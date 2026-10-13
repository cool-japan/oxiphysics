//! Physics regression tests for the particle-laden LBM helpers: Schiller-Naumann
//! drag, net weight/buoyancy, quiescent-fluid settling and soft-sphere contact
//! forces.

use super::*;

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1e-300)
}

// ---------------------------------------------------------------------------
// Schiller-Naumann drag
// ---------------------------------------------------------------------------

#[test]
fn schiller_naumann_force_matches_stokes_force_at_low_reynolds() {
    // Water, d = 10 um, slip speed 1 mm/s -> Re_p = 1e-2.
    let (rho_f, nu, d) = (1000.0, 1e-6, 1e-5);
    let fluid_vel = [1e-3, 0.0, 0.0];
    let particle_vel = [0.0; 3];
    let sn = drag_force_schiller_naumann(&fluid_vel, &particle_vel, d, rho_f, nu);
    let stokes = stokes_drag_force(&fluid_vel, &particle_vel, d, rho_f * nu);
    assert!(
        rel_err(sn[0], stokes[0]) < 0.01,
        "S-N {} vs Stokes {}",
        sn[0],
        stokes[0]
    );
    assert_eq!(sn[1], 0.0);
    assert_eq!(sn[2], 0.0);
}

#[test]
fn schiller_naumann_force_equals_closed_form_expression() {
    // Slip speed |u| = 2 along a general direction; Re_p = 20.
    let (rho_f, nu, d) = (1000.0, 1e-6, 1e-5);
    let fluid_vel = [1.2, 1.6, 0.0];
    let particle_vel = [0.0; 3];
    let speed = 2.0_f64;
    let re_p = speed * d / nu;
    let expected_mag = 3.0 * PI * rho_f * nu * d * speed * (1.0 + 0.15 * re_p.powf(0.687));
    let f = drag_force_schiller_naumann(&fluid_vel, &particle_vel, d, rho_f, nu);
    let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
    assert!(
        rel_err(mag, expected_mag) < 1e-10,
        "|F| = {mag}, expected {expected_mag}"
    );
    // The force points along the slip velocity.
    assert!(rel_err(f[0] / mag, 0.6) < 1e-12);
    assert!(rel_err(f[1] / mag, 0.8) < 1e-12);
    assert_eq!(f[2], 0.0);
}

#[test]
fn schiller_naumann_force_doubles_with_slip_speed_at_low_reynolds() {
    let (rho_f, nu, d) = (1000.0, 1e-6, 1e-5);
    let slow = drag_force_schiller_naumann(&[1e-3, 0.0, 0.0], &[0.0; 3], d, rho_f, nu);
    let fast = drag_force_schiller_naumann(&[2e-3, 0.0, 0.0], &[0.0; 3], d, rho_f, nu);
    let ratio = fast[0] / slow[0];
    // Re_p = 1e-2 and 2e-2: the force is linear in the slip speed up to the
    // small (0.15 Re^0.687) inertial correction.
    assert!(
        (ratio - 2.0).abs() < 1e-2,
        "force ratio for doubled slip is {ratio}, expected ~2"
    );
    let correction = |re: f64| 1.0 + 0.15 * re.powf(0.687);
    let expected_ratio = 2.0 * correction(2e-2) / correction(1e-2);
    assert!(
        rel_err(ratio, expected_ratio) < 1e-10,
        "force ratio {ratio} vs closed form {expected_ratio}"
    );
}

#[test]
fn schiller_naumann_force_is_zero_without_slip() {
    let f = drag_force_schiller_naumann(&[0.3, -0.2, 0.1], &[0.3, -0.2, 0.1], 1e-3, 1000.0, 1e-6);
    assert_eq!(f, [0.0; 3]);
}

#[test]
fn schiller_naumann_force_scales_quadratically_in_newton_regime() {
    // d = 1 mm, water: slip 2 m/s -> Re_p = 2000 and 4 m/s -> Re_p = 4000, both
    // above the correlation limit, so C_D = 0.44 and F ~ |u|^2.
    let (rho_f, nu, d) = (1000.0, 1e-6, 1e-3);
    let f1 = drag_force_schiller_naumann(&[2.0, 0.0, 0.0], &[0.0; 3], d, rho_f, nu);
    let f2 = drag_force_schiller_naumann(&[4.0, 0.0, 0.0], &[0.0; 3], d, rho_f, nu);
    assert!(rel_err(f2[0] / f1[0], 4.0) < 1e-12);
    let area = PI * d * d / 4.0;
    let expected = 0.5 * NEWTON_REGIME_DRAG_COEFFICIENT * rho_f * area * 4.0;
    assert!(rel_err(f1[0], expected) < 1e-12, "F = {}", f1[0]);
}

#[test]
fn schiller_naumann_coefficient_is_capped_in_newton_regime() {
    assert_eq!(
        schiller_naumann_drag_coefficient(1e3),
        NEWTON_REGIME_DRAG_COEFFICIENT
    );
    assert_eq!(
        schiller_naumann_drag_coefficient(5e5),
        NEWTON_REGIME_DRAG_COEFFICIENT
    );
    // Continuity across the limit: the correlation itself is within 0.5 % of
    // the plateau at Re = 1000.
    let below = schiller_naumann_drag_coefficient(SCHILLER_NAUMANN_RE_LIMIT * (1.0 - 1e-9));
    assert!(rel_err(below, NEWTON_REGIME_DRAG_COEFFICIENT) < 5e-3);
}

// ---------------------------------------------------------------------------
// Net weight / buoyancy and settling
// ---------------------------------------------------------------------------

#[test]
fn heavy_particle_settles_at_stokes_velocity_in_quiescent_fluid() {
    let field = VelocityField2D::uniform(0.0, 0.0, 20, 20, 0.1);
    let fluid = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
    let (d, rho_p) = (1e-4, 2500.0);
    let mut particles = vec![LagrangianParticle::new([1.0, 1.0, 0.0], d, rho_p, 0)];
    let tau_p = particle_relaxation_time(d, rho_p, fluid.mu);
    let dt = tau_p / 20.0;
    // 400 steps = 20 relaxation times: the transient has decayed to ~1e-9.
    for _ in 0..400 {
        track_particles_2d(&mut particles, &field, &fluid, dt);
    }
    let v_stokes = (rho_p - fluid.rho_f) * 9.81 * d * d / (18.0 * fluid.mu);
    let p = &particles[0];
    assert!(p.vel[1] < 0.0, "a heavy grain must sink, vy = {}", p.vel[1]);
    assert!(
        rel_err(-p.vel[1], v_stokes) < 1e-6,
        "vy = {}, Stokes settling velocity = -{v_stokes}",
        p.vel[1]
    );
    assert!(p.vel[0].abs() < 1e-14);
    assert!(p.pos[1] < 1.0, "the grain must have moved towards -y");
    assert!((p.pos[0] - 1.0).abs() < 1e-14);
}

#[test]
fn light_particle_rises_at_stokes_velocity_in_quiescent_fluid() {
    let field = VelocityField2D::uniform(0.0, 0.0, 20, 20, 0.1);
    let fluid = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
    let (d, rho_p) = (1e-4, 500.0);
    let mut particles = vec![LagrangianParticle::new([1.0, 1.0, 0.0], d, rho_p, 0)];
    let tau_p = particle_relaxation_time(d, rho_p, fluid.mu);
    let dt = tau_p / 20.0;
    for _ in 0..400 {
        track_particles_2d(&mut particles, &field, &fluid, dt);
    }
    let v_rise = (fluid.rho_f - rho_p) * 9.81 * d * d / (18.0 * fluid.mu);
    let p = &particles[0];
    assert!(p.vel[1] > 0.0);
    assert!(rel_err(p.vel[1], v_rise) < 1e-6, "vy = {}", p.vel[1]);
}

#[test]
fn neutrally_buoyant_particle_stays_at_rest() {
    let field = VelocityField2D::uniform(0.0, 0.0, 20, 20, 0.1);
    let fluid = FluidParams::new(1000.0, 1e-6, [0.0, -9.81, 0.0]);
    let mut particles = vec![LagrangianParticle::new([1.0, 1.0, 0.0], 1e-4, 1000.0, 0)];
    for _ in 0..100 {
        track_particles_2d(&mut particles, &field, &fluid, 1e-4);
    }
    assert_eq!(particles[0].vel, [0.0; 3]);
    assert_eq!(particles[0].pos, [1.0, 1.0, 0.0]);
}

// ---------------------------------------------------------------------------
// Soft-sphere DEM: stateless contact
// ---------------------------------------------------------------------------

fn overlapping_pair() -> (LagrangianParticle, LagrangianParticle) {
    // Two d = 0.1 spheres whose centres are 0.08 apart -> overlap 0.02.
    let p1 = LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 0);
    let p2 = LagrangianParticle::new([0.08, 0.0, 0.0], 0.1, 1000.0, 1);
    (p1, p2)
}

#[test]
fn dem_normal_damping_adds_repulsion_while_approaching() {
    let params = DemParams::new(1000.0, 500.0, 10.0, 0.0, 0.3);
    let (mut p1, p2) = overlapping_pair();

    // Static: spring only, 1000 * 0.02 = 20 pushing p1 towards -x.
    let (f_static, _) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(-f_static[0], 20.0) < 1e-12);

    // p1 moves towards p2 (+x): approach rate 1 -> 20 + 10 * 1 = 30.
    p1.vel = [1.0, 0.0, 0.0];
    let (f_approach, _) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(-f_approach[0], 30.0) < 1e-12);

    // p1 moves away from p2 (-x): 20 - 10 * 1 = 10.
    p1.vel = [-1.0, 0.0, 0.0];
    let (f_recede, _) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(-f_recede[0], 10.0) < 1e-12);

    // Fast separation would make the damping term exceed the spring: the
    // contact is never attractive.
    p1.vel = [-10.0, 0.0, 0.0];
    let (f_fast, f_fast2) = soft_sphere_collision(&p1, &p2, &params);
    assert_eq!(f_fast, [0.0; 3]);
    assert_eq!(f_fast2, [0.0; 3]);
}

#[test]
fn dem_forces_obey_newtons_third_law() {
    let params = DemParams::new(2000.0, 800.0, 12.0, 6.0, 0.4);
    let mut p1 = LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 0);
    let mut p2 = LagrangianParticle::new([0.05, 0.04, 0.03], 0.1, 1000.0, 1);
    p1.vel = [0.3, -0.2, 0.5];
    p2.vel = [-0.1, 0.4, 0.2];
    let (f1, f2) = soft_sphere_collision(&p1, &p2, &params);
    assert!(f1.iter().any(|v| v.abs() > 0.0));
    for k in 0..3 {
        assert!((f1[k] + f2[k]).abs() < 1e-12, "component {k}");
    }
}

#[test]
fn dem_friction_opposes_the_sliding_of_the_first_particle() {
    let params = DemParams::new(1000.0, 500.0, 0.0, 5.0, 0.3);
    let (mut p1, p2) = overlapping_pair();

    // p1 slides along +y relative to p2: normal 20 (towards -x), viscous
    // friction 5 * 0.1 = 0.5 < mu * Fn = 6, acting along -y on p1.
    p1.vel = [0.0, 0.1, 0.0];
    let (f1, f2) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(-f1[0], 20.0) < 1e-12);
    assert!(
        f1[1] < 0.0,
        "friction must oppose +y sliding, got {}",
        f1[1]
    );
    assert!(rel_err(-f1[1], 0.5) < 1e-12);
    assert!(rel_err(f2[1], 0.5) < 1e-12, "p2 is dragged along +y");

    // Opposite sliding direction flips the friction force.
    p1.vel = [0.0, -0.1, 0.0];
    let (g1, g2) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(g1[1], 0.5) < 1e-12);
    assert!(rel_err(-g2[1], 0.5) < 1e-12);
}

#[test]
fn dem_friction_is_capped_by_the_coulomb_limit() {
    let params = DemParams::new(1000.0, 500.0, 0.0, 5.0, 0.3);
    let (mut p1, p2) = overlapping_pair();
    // Viscous friction 5 * 10 = 50 would exceed mu * Fn = 0.3 * 20 = 6.
    p1.vel = [0.0, 10.0, 0.0];
    let (f1, _) = soft_sphere_collision(&p1, &p2, &params);
    assert!(rel_err(-f1[1], 6.0) < 1e-12, "capped friction {}", f1[1]);
}

#[test]
fn dem_friction_vanishes_without_slip_even_with_overlap() {
    // The old model produced a tangential force proportional to the overlap
    // even for pure normal contact.
    let params = DemParams::new(1000.0, 500.0, 5.0, 5.0, 0.3);
    let (p1, p2) = overlapping_pair();
    let (f1, _) = soft_sphere_collision(&p1, &p2, &params);
    assert_eq!(f1[1], 0.0);
    assert_eq!(f1[2], 0.0);
}

/// Integrate a head-on collision of two equal spheres and return the
/// restitution coefficient (separation speed / approach speed).
fn head_on_restitution(params: &DemParams) -> f64 {
    let d = 0.1;
    let mut p = vec![
        LagrangianParticle::new([-0.06, 0.0, 0.0], d, 1000.0, 0),
        LagrangianParticle::new([0.06, 0.0, 0.0], d, 1000.0, 1),
    ];
    p[0].vel = [0.5, 0.0, 0.0];
    p[1].vel = [-0.5, 0.0, 0.0];
    let approach_speed = 1.0;
    let dt = 1e-5;
    let mut touched = false;
    for _ in 0..200_000 {
        for q in p.iter_mut() {
            q.reset_forces();
        }
        compute_dem_collisions(&mut p, params);
        let in_contact = p[0].force != [0.0; 3];
        touched |= in_contact;
        for q in p.iter_mut() {
            euler_step(q, dt);
        }
        if touched && !in_contact {
            break;
        }
    }
    assert!(touched, "the spheres never touched");
    (p[1].vel[0] - p[0].vel[0]) / approach_speed
}

#[test]
fn dem_head_on_collision_dissipates_energy() {
    // m = rho * pi d^3 / 6, m_eff = m / 2; kn = 1e4, cn = 20 -> zeta = 0.195.
    let params = DemParams::new(1.0e4, 0.0, 20.0, 0.0, 0.0);
    let e = head_on_restitution(&params);
    let m_eff = 0.5 * 1000.0 * PI * 0.1_f64.powi(3) / 6.0;
    let zeta = params.cn / (2.0 * (params.kn * m_eff).sqrt());
    let e_linear = (-PI * zeta / (1.0 - zeta * zeta).sqrt()).exp();
    assert!(e > 0.0 && e < 1.0, "restitution must be in (0, 1): {e}");
    assert!(
        (e - e_linear).abs() < 0.05,
        "restitution {e} vs linear spring-dashpot value {e_linear}"
    );
}

#[test]
fn dem_undamped_head_on_collision_is_elastic() {
    let params = DemParams::new(1.0e4, 0.0, 0.0, 0.0, 0.0);
    let e = head_on_restitution(&params);
    assert!((e - 1.0).abs() < 5e-3, "elastic restitution {e}");
}

// ---------------------------------------------------------------------------
// Soft-sphere DEM: Cundall-Strack tangential history
// ---------------------------------------------------------------------------

#[test]
fn dem_history_tangential_force_grows_then_saturates_at_the_coulomb_limit() {
    // Fn = 20, mu = 0.3 -> Coulomb limit 6.  Slip 0.2 over dt = 0.01 adds
    // xi = 0.002 per step, i.e. 0.4 of spring force per step (kt = 200).
    let params = DemParams::new(1000.0, 200.0, 0.0, 0.0, 0.3);
    let (mut p1, p2) = overlapping_pair();
    p1.vel = [0.0, 0.2, 0.0];
    let dt = 0.01;
    let mut spring = [0.0_f64; 3];
    for step in 1..=40 {
        let (f1, f2) = soft_sphere_collision_with_history(&p1, &p2, &params, dt, &mut spring);
        let expected = (0.4 * step as f64).min(6.0);
        assert!(rel_err(-f1[0], 20.0) < 1e-12);
        assert!(
            rel_err(-f1[1], expected) < 1e-9,
            "step {step}: tangential force {} expected -{expected}",
            f1[1]
        );
        for k in 0..3 {
            assert!((f1[k] + f2[k]).abs() < 1e-12);
        }
    }
    // Sliding: the spring is truncated to the displacement that produces the
    // capped force.
    assert!(rel_err(spring[1], 6.0 / 200.0) < 1e-9, "xi = {}", spring[1]);
    assert!(spring[0].abs() < 1e-15);
}

#[test]
fn dem_history_spring_is_reprojected_onto_the_contact_plane() {
    let params = DemParams::new(1000.0, 200.0, 0.0, 0.0, 0.3);
    let (p1, p2) = overlapping_pair();
    // A stored displacement with a component along the contact normal (x).
    let mut spring = [0.01, 0.02, 0.0];
    let _ = soft_sphere_collision_with_history(&p1, &p2, &params, 0.01, &mut spring);
    assert!(
        spring[0].abs() < 1e-15,
        "normal part removed: {}",
        spring[0]
    );
    assert!((spring[1] - 0.02).abs() < 1e-15);
}

#[test]
fn dem_history_resets_when_the_contact_breaks() {
    let params = DemParams::new(1000.0, 200.0, 0.0, 0.0, 0.3);
    let (mut p1, mut p2) = overlapping_pair();
    p1.vel = [0.0, 0.2, 0.0];
    let mut spring = [0.0_f64; 3];
    let _ = soft_sphere_collision_with_history(&p1, &p2, &params, 0.01, &mut spring);
    assert!(spring[1] > 0.0);
    p2.pos[0] = 1.0;
    let (f1, f2) = soft_sphere_collision_with_history(&p1, &p2, &params, 0.01, &mut spring);
    assert_eq!(f1, [0.0; 3]);
    assert_eq!(f2, [0.0; 3]);
    assert_eq!(spring, [0.0; 3]);
}

#[test]
fn dem_history_without_tangential_stiffness_reduces_to_the_stateless_model() {
    let params = DemParams::new(1000.0, 0.0, 10.0, 0.0, 0.3);
    let (mut p1, p2) = overlapping_pair();
    p1.vel = [0.4, 0.3, -0.2];
    let mut spring = [0.0_f64; 3];
    let (h1, h2) = soft_sphere_collision_with_history(&p1, &p2, &params, 0.01, &mut spring);
    let (s1, s2) = soft_sphere_collision(&p1, &p2, &params);
    for k in 0..3 {
        assert!((h1[k] - s1[k]).abs() < 1e-12, "force 1, component {k}");
        assert!((h2[k] - s2[k]).abs() < 1e-12, "force 2, component {k}");
    }
}

#[test]
fn dem_contact_history_tracks_only_touching_pairs() {
    let params = DemParams::new(1000.0, 200.0, 0.0, 0.0, 0.3);
    let mut particles = vec![
        LagrangianParticle::new([0.0, 0.0, 0.0], 0.1, 1000.0, 10),
        LagrangianParticle::new([0.08, 0.0, 0.0], 0.1, 1000.0, 11),
        LagrangianParticle::new([5.0, 0.0, 0.0], 0.1, 1000.0, 12),
    ];
    particles[0].vel = [0.0, 0.2, 0.0];
    let mut history = DemContactHistory::new();
    assert!(history.is_empty());

    compute_dem_collisions_with_history(&mut particles, &params, &mut history, 0.01);
    assert_eq!(history.len(), 1);
    let xi = history
        .tangential_displacement(11, 10)
        .expect("pair (10, 11) is in contact");
    assert!(rel_err(xi[1], 0.002) < 1e-12);
    assert_eq!(history.tangential_displacement(10, 12), None);
    // Forces were accumulated on the touching pair only.
    assert!(particles[0].force[0] < 0.0 && particles[1].force[0] > 0.0);
    assert_eq!(particles[2].force, [0.0; 3]);

    // Pull the pair apart: the spring is forgotten.
    particles[1].pos[0] = 1.0;
    for p in particles.iter_mut() {
        p.reset_forces();
    }
    compute_dem_collisions_with_history(&mut particles, &params, &mut history, 0.01);
    assert!(history.is_empty());
    history.clear();
    assert_eq!(history.len(), 0);
}
