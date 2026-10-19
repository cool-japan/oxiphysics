//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::types::{DemContact, DemParticle, DemSimulation, GranularContactModel};

#[cfg(test)]
mod tests_dem_contact {
    use super::*;
    fn make_contact(overlap: f64) -> DemContact {
        DemContact::new(0, 1, overlap, [1.0, 0.0, 0.0], 1.0e9, 0.05, 0.5, 0.01)
    }
    #[test]
    fn hertz_contact_force_zero_overlap_returns_zero() {
        let c = make_contact(0.0);
        let f = c.compute_hertz_contact_force();
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn hertz_contact_force_negative_overlap_returns_zero() {
        let c = make_contact(-0.001);
        let f = c.compute_hertz_contact_force();
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn hertz_contact_force_positive_overlap_repulsive_on_i() {
        let c = make_contact(0.001);
        let f = c.compute_hertz_contact_force();
        assert!(
            f[0] < 0.0,
            "Force on i should be in −x (repulsive), got {}",
            f[0]
        );
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn hertz_contact_force_scales_as_delta_to_3_2() {
        let c1 = make_contact(0.001);
        let c2 = make_contact(0.002);
        let f1 = c1.compute_hertz_contact_force()[0].abs();
        let f2 = c2.compute_hertz_contact_force()[0].abs();
        let ratio = f2 / f1;
        let expected = 2.0_f64.powf(1.5);
        assert!(
            (ratio - expected).abs() < 1e-9,
            "ratio={ratio} expected={expected}"
        );
    }
    #[test]
    fn rolling_resistance_zero_relative_angular_velocity() {
        let c = make_contact(0.001);
        let tau = c.compute_rolling_resistance([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(tau, [0.0; 3], "Zero relative ω → zero torque");
    }
    #[test]
    fn rolling_resistance_zero_overlap_returns_zero() {
        let c = make_contact(0.0);
        let tau = c.compute_rolling_resistance([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        assert_eq!(tau, [0.0; 3]);
    }
    #[test]
    fn rolling_resistance_direction_opposes_relative_spin() {
        let c = make_contact(0.001);
        let tau = c.compute_rolling_resistance([0.0, 0.0, 1.0], [0.0; 3]);
        assert!(
            tau[2] < 0.0,
            "Torque should oppose relative spin, tau_z={}",
            tau[2]
        );
    }
    #[test]
    fn coordination_number_no_contact() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let p0 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p1 = DemParticle::new_sphere([10.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle(p0);
        sim.add_particle(p1);
        let z = sim.compute_coordination_number(0.0);
        assert_eq!(z, 0.0, "Far-apart particles have Z=0, got {z}");
    }
    #[test]
    fn coordination_number_one_contact_two_particles() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let p0 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p1 = DemParticle::new_sphere([0.1, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle(p0);
        sim.add_particle(p1);
        let z = sim.compute_coordination_number(0.0);
        assert!((z - 1.0).abs() < 1e-12, "Expected Z̄=1.0, got {z}");
    }
    #[test]
    fn coordination_number_empty_simulation_is_zero() {
        let sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let z = sim.compute_coordination_number(0.0);
        assert_eq!(z, 0.0);
    }
}

#[cfg(test)]
mod tests_tangential_friction {
    //! Regression tests of the tangential (friction) contact forces of the three
    //! DEM drivers: friction opposes the sliding, is capped by the Coulomb cone,
    //! keeps a per-contact spring history, spins the particles up and conserves
    //! the linear and the angular momentum.
    use crate::granular::functions::*;
    use crate::granular::types::*;

    const DT: f64 = 1e-5;

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }

    fn spring_dashpot(kt: f64, gamma_t: f64) -> ContactModel {
        ContactModel::SpringDashpot {
            kn: 1e5,
            kt,
            gamma_n: 0.0,
            gamma_t,
        }
    }

    /// Unit-mass spheres of radius 0.1 overlapping by 0.02 along y: the first
    /// sits at the origin and slides with −x velocity, the second sits above it
    /// and slides with +x velocity (relative sliding speed 2, normal force 2000 N
    /// for `kn = 1e5`).
    const PAIR_POSITIONS: [[f64; 3]; 2] = [[0.0, 0.0, 0.0], [0.0, 0.18, 0.0]];
    const PAIR_VELOCITIES: [[f64; 3]; 2] = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

    fn sliding_granular_sim(model: ContactModel, mu: f64) -> GranularSim {
        let mut sim = GranularSim::new(model, mu, [0.0; 3]);
        for (pos, vel) in PAIR_POSITIONS.iter().zip(PAIR_VELOCITIES.iter()) {
            sim.add_particle(*pos, *vel, 0.1, 1.0);
        }
        sim
    }

    fn sliding_extended_sim(
        model: ContactModel,
        mu: f64,
        rolling: RollingResistanceModel,
        mu_r: f64,
        masses: [f64; 2],
    ) -> GranularSimExtended {
        let mut sim = GranularSimExtended::new(model, mu, [0.0; 3], rolling, mu_r);
        for ((pos, vel), mass) in PAIR_POSITIONS
            .iter()
            .zip(PAIR_VELOCITIES.iter())
            .zip(masses.iter())
        {
            sim.add_particle(*pos, *vel, 0.1, *mass);
        }
        sim
    }

    /// [`DemSimulation`] of the sliding pair with a linear normal spring
    /// (`kn = 1e5`), the tangential spring `kt` and the friction coefficient `mu`.
    fn sliding_dem_sim(kt: f64, mu: f64) -> DemSimulation {
        let mut sim = DemSimulation::new(
            DT,
            [0.0; 3],
            1e5,
            kt,
            0.0,
            mu,
            GranularContactModel::LinearSpring,
            1e9,
        );
        for (pos, vel) in PAIR_POSITIONS.iter().zip(PAIR_VELOCITIES.iter()) {
            sim.add_particle(DemParticle::new_sphere(*pos, *vel, 0.1, 1.0));
        }
        sim
    }

    /// Total angular momentum about the origin, `Σ m x × v + I ω`.
    fn extended_angular_momentum(sim: &GranularSimExtended) -> [f64; 3] {
        let mut total = [0.0; 3];
        for i in 0..sim.len() {
            let orbital = cross(sim.sim.positions[i], sim.sim.velocities[i]);
            let inertia = 0.4 * sim.sim.masses[i] * sim.sim.radii[i].powi(2);
            for (k, t) in total.iter_mut().enumerate() {
                *t += sim.sim.masses[i] * orbital[k] + inertia * sim.angular_velocities[i][k];
            }
        }
        total
    }

    fn extended_slip_speed(sim: &GranularSimExtended) -> f64 {
        let (x, v, w, r) = (
            &sim.sim.positions,
            &sim.sim.velocities,
            &sim.angular_velocities,
            &sim.sim.radii,
        );
        let d = [x[1][0] - x[0][0], x[1][1] - x[0][1], x[1][2] - x[0][2]];
        let dist = norm(d);
        let normal = [d[0] / dist, d[1] / dist, d[2] / dist];
        let lever = r[0] - 0.5 * (r[0] + r[1] - dist);
        let slip = contact_point_relative_velocity(v[0], v[1], w[0], w[1], lever, lever, normal);
        norm(tangential_component(slip, normal))
    }

    // ── Pure helpers ─────────────────────────────────────────────────────────

    #[test]
    fn contact_point_velocity_vanishes_for_a_rolling_sphere() {
        // Sphere i of radius r rolls to +x without slipping on the fixed sphere j
        // below it: ω = −(v/r) ẑ, normal i → j = −ŷ, lever = r.
        let (v, r) = (0.7, 0.1);
        let slip = contact_point_relative_velocity(
            [v, 0.0, 0.0],
            [0.0; 3],
            [0.0, 0.0, -v / r],
            [0.0; 3],
            r,
            0.0,
            [0.0, -1.0, 0.0],
        );
        assert!(norm(slip) < 1e-12, "rolling contact slips: {slip:?}");
        // The same translation without spin slips with the full velocity.
        let skid = contact_point_relative_velocity(
            [v, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            r,
            0.0,
            [0.0, -1.0, 0.0],
        );
        assert!((skid[0] - v).abs() < 1e-14);
    }

    #[test]
    fn tangential_component_removes_the_normal_part() {
        let t = tangential_component([1.0, 2.0, 3.0], [0.0, 1.0, 0.0]);
        assert_eq!(t, [1.0, 0.0, 3.0]);
    }

    // ── TangentialContactLaw ─────────────────────────────────────────────────

    #[test]
    fn tangential_law_in_the_stick_regime_is_spring_plus_dashpot() {
        let law = TangentialContactLaw {
            kt: 5e4,
            gamma_t: 50.0,
            mu: 0.5,
        };
        let mut spring = [0.0; 3];
        let normal = [0.0, 1.0, 0.0];
        let v_t = [-2.0, 0.0, 0.0];
        let f1 = law.force(&mut spring, normal, v_t, 2000.0, DT);
        // δ = v Δt, F = −(k δ + γ v) = +(1 + 100) x̂ on i.
        assert!((spring[0] + 2.0 * DT).abs() < 1e-18);
        assert!((f1[0] - 101.0).abs() < 1e-10, "F = {f1:?}");
        let f2 = law.force(&mut spring, normal, v_t, 2000.0, DT);
        assert!((spring[0] + 4.0 * DT).abs() < 1e-18);
        assert!((f2[0] - 102.0).abs() < 1e-10, "F = {f2:?}");
    }

    #[test]
    fn tangential_law_slip_caps_the_force_and_keeps_the_spring_on_the_cone() {
        let law = TangentialContactLaw {
            kt: 5e4,
            gamma_t: 1000.0,
            mu: 0.3,
        };
        let mut spring = [0.0; 3];
        let normal = [0.0, 1.0, 0.0];
        let v_t = [-2.0, 0.0, 0.0];
        for _ in 0..4 {
            let f = law.force(&mut spring, normal, v_t, 2000.0, DT);
            assert!((norm(f) - 600.0).abs() < 1e-9, "|F| = {}", norm(f));
            assert!(f[0] > 0.0, "friction must oppose the −x sliding of i");
            // The stored spring + damper force sits exactly on the cone.
            let stored = [
                -(law.kt * spring[0] + law.gamma_t * v_t[0]),
                -(law.kt * spring[1] + law.gamma_t * v_t[1]),
                -(law.kt * spring[2] + law.gamma_t * v_t[2]),
            ];
            assert!((norm(stored) - 600.0).abs() < 1e-9);
        }
    }

    #[test]
    fn tangential_law_without_normal_force_gives_no_friction() {
        let law = TangentialContactLaw {
            kt: 5e4,
            gamma_t: 0.0,
            mu: 0.5,
        };
        let mut spring = [1e-3, 0.0, 0.0];
        let f = law.force(&mut spring, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 0.0, DT);
        assert_eq!(f, [0.0; 3]);
        assert_eq!(spring, [0.0; 3], "a contact without load keeps no history");
        let f = law.force(&mut spring, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 100.0, DT);
        assert!(
            (f[0] + law.kt * DT).abs() < 1e-9,
            "history restarts at zero"
        );
    }

    #[test]
    fn tangential_law_zero_load_clears_the_history_with_a_damper_too() {
        // The slip branch alone would leave the spring at −(γ/k) v_t for γ > 0;
        // a contact that carries no load must not keep any history at all, and a
        // frictionless contact (μ = 0) neither.
        for (mu, fn_mag) in [(0.5, 0.0), (0.0, 100.0)] {
            let law = TangentialContactLaw {
                kt: 5e4,
                gamma_t: 50.0,
                mu,
            };
            let mut spring = [1e-3, 0.0, 2e-4];
            let f = law.force(&mut spring, [0.0, 1.0, 0.0], [-2.0, 0.0, 0.0], fn_mag, DT);
            assert_eq!(f, [0.0; 3], "μ = {mu}, F_n = {fn_mag}");
            assert_eq!(spring, [0.0; 3], "μ = {mu}, F_n = {fn_mag}");
        }
        // The next loaded step restarts from a fresh spring: k δ + γ v = (2 k Δt + 100).
        let law = TangentialContactLaw {
            kt: 5e4,
            gamma_t: 50.0,
            mu: 0.5,
        };
        let mut spring = [0.0; 3];
        law.force(&mut spring, [0.0, 1.0, 0.0], [-2.0, 0.0, 0.0], 0.0, DT);
        let f = law.force(&mut spring, [0.0, 1.0, 0.0], [-2.0, 0.0, 0.0], 2000.0, DT);
        assert!((f[0] - (5e4 * 2.0 * DT + 100.0)).abs() < 1e-9, "F = {f:?}");
    }

    #[test]
    fn tangential_law_without_spring_is_a_capped_dashpot() {
        let law = TangentialContactLaw {
            kt: 0.0,
            gamma_t: 10.0,
            mu: 0.5,
        };
        let mut spring = [0.0; 3];
        let f = law.force(&mut spring, [0.0, 1.0, 0.0], [-3.0, 0.0, 0.0], 100.0, DT);
        assert!(
            (f[0] - 30.0).abs() < 1e-12,
            "viscous friction 10·3, got {f:?}"
        );
        let capped = law.force(&mut spring, [0.0, 1.0, 0.0], [-30.0, 0.0, 0.0], 100.0, DT);
        assert!(
            (capped[0] - 50.0).abs() < 1e-12,
            "μ F_n cap, got {capped:?}"
        );
        assert_eq!(spring, [0.0; 3]);
    }

    #[test]
    fn tangential_law_carries_the_spring_into_the_new_tangent_plane() {
        // A spring displacement along z is rotated with a normal that tilts from y
        // towards z: it stays perpendicular to the new normal and keeps its length.
        let law = TangentialContactLaw {
            kt: 1e3,
            gamma_t: 0.0,
            mu: 10.0,
        };
        let mut spring = [0.0, 0.0, 2e-3];
        let tilt = 0.3_f64;
        let normal = [0.0, tilt.cos(), tilt.sin()];
        law.force(&mut spring, normal, [0.0; 3], 100.0, DT);
        assert!(dot(spring, normal).abs() < 1e-15, "spring not tangential");
        assert!((norm(spring) - 2e-3).abs() < 1e-15, "spring length changed");
        // A spring parallel to the new normal has no tangential part left.
        let mut parallel = [0.0, 1e-3, 0.0];
        law.force(&mut parallel, [0.0, 1.0, 0.0], [0.0; 3], 100.0, DT);
        assert_eq!(parallel, [0.0; 3]);
    }

    // ── GranularSim ──────────────────────────────────────────────────────────

    #[test]
    fn granular_sim_friction_opposes_the_relative_sliding() {
        let mut sim = sliding_granular_sim(spring_dashpot(5e4, 50.0), 0.5);
        let contacts = sim.detect_contacts();
        let loads = sim.resolve_contacts(&contacts, DT);
        assert_eq!(loads.len(), 1);
        // Particle 0 slides with −2 x̂ relative to particle 1: F_t = +(kt·2Δt + γt·2).
        let expected = 5e4 * 2.0 * DT + 50.0 * 2.0;
        assert!((loads[0].tangential[0] - expected).abs() < 1e-9);
        assert!(loads[0].tangential[1].abs() < 1e-12 && loads[0].tangential[2].abs() < 1e-12);
        assert!((loads[0].normal - 2000.0).abs() < 1e-9);
        assert_eq!((loads[0].i, loads[0].j), (0, 1));
        // The forces are the contact force and its reaction (no gravity here).
        assert!((sim.forces[0][0] - expected).abs() < 1e-9);
        assert!((sim.forces[1][0] + expected).abs() < 1e-9);
        assert!((sim.forces[0][1] + 2000.0).abs() < 1e-9);
        for k in 0..3 {
            assert!((sim.forces[0][k] + sim.forces[1][k]).abs() < 1e-9);
        }
    }

    #[test]
    fn granular_sim_friction_is_limited_by_the_coulomb_cone() {
        let mut sim = sliding_granular_sim(spring_dashpot(5e4, 1000.0), 0.3);
        let contacts = sim.detect_contacts();
        let loads = sim.resolve_contacts(&contacts, DT);
        assert!((norm(loads[0].tangential) - 0.3 * loads[0].normal).abs() < 1e-9);
        assert!((loads[0].tangential[0] - 600.0).abs() < 1e-9);
    }

    #[test]
    fn granular_sim_sliding_decreases_and_momentum_is_conserved() {
        let mut sim = sliding_granular_sim(spring_dashpot(5e4, 1000.0), 0.3);
        let slip_before = sim.velocities[1][0] - sim.velocities[0][0];
        assert!((slip_before - 2.0).abs() < 1e-14);
        for _ in 0..5 {
            sim.step(DT);
        }
        let slip_after = sim.velocities[1][0] - sim.velocities[0][0];
        // Each particle loses μ F_n Δt / m ≈ 600 · 5·10⁻⁵ = 0.03 of its speed.
        assert!(
            (slip_after - 1.94).abs() < 0.01,
            "relative sliding speed {slip_before} -> {slip_after}"
        );
        assert!((sim.velocities[0][0] + sim.velocities[1][0]).abs() < 1e-12);
        assert!((sim.velocities[0][1] + sim.velocities[1][1]).abs() < 1e-12);
        assert!(sim.velocities[1][0] > 0.9 && sim.velocities[0][0] < -0.9);
    }

    #[test]
    fn granular_sim_hertzian_friction_uses_the_mindlin_stiffness() {
        let model = ContactModel::Hertzian {
            e_star: 1e6,
            g_star: 4e5,
            r_star: 0.05,
        };
        let mut sim = sliding_granular_sim(model, 0.3);
        let contacts = sim.detect_contacts();
        let loads = sim.resolve_contacts(&contacts, DT);
        // k_t = 8 G* √(R* δ_n), F_t = k_t · 2Δt, F_n = (4/3) E* √R* δ^{3/2}.
        let kt = 8.0 * 4e5 * (0.05_f64 * 0.02).sqrt();
        assert!((loads[0].tangential[0] - kt * 2.0 * DT).abs() < 1e-9 * kt * DT);
        let fn_expected = (4.0 / 3.0) * 1e6 * 0.05_f64.sqrt() * 0.02_f64.powf(1.5);
        assert!((loads[0].normal - fn_expected).abs() < 1e-9 * fn_expected);
        assert!(loads[0].tangential[0] > 0.0);
    }

    #[test]
    fn granular_sim_contact_history_accumulates_and_is_dropped_on_separation() {
        let mut sim = sliding_granular_sim(spring_dashpot(5e4, 0.0), 0.5);
        let contacts = sim.detect_contacts();
        sim.resolve_contacts(&contacts, DT);
        assert_eq!(sim.tangential_springs.len(), 1);
        assert!((sim.tangential_springs[&(0, 1)][0] + 2.0 * DT).abs() < 1e-18);
        sim.resolve_contacts(&contacts, DT);
        assert!((sim.tangential_springs[&(0, 1)][0] + 4.0 * DT).abs() < 1e-18);

        // The pair separates: its history is dropped ...
        sim.positions[1] = [0.0, 1.0, 0.0];
        let apart = sim.detect_contacts();
        assert!(apart.is_empty());
        assert!(sim.resolve_contacts(&apart, DT).is_empty());
        assert!(sim.tangential_springs.is_empty());

        // ... and a new contact starts with a fresh spring.
        sim.positions[1] = [0.0, 0.18, 0.0];
        let again = sim.detect_contacts();
        sim.resolve_contacts(&again, DT);
        assert!((sim.tangential_springs[&(0, 1)][0] + 2.0 * DT).abs() < 1e-18);
    }

    #[test]
    fn granular_sim_without_sliding_has_no_friction() {
        let mut sim = GranularSim::new(spring_dashpot(5e4, 50.0), 0.5, [0.0; 3]);
        sim.add_particle([0.0, 0.0, 0.0], [0.3, 0.0, 0.0], 0.1, 1.0);
        sim.add_particle([0.0, 0.18, 0.0], [0.3, 0.0, 0.0], 0.1, 1.0);
        let contacts = sim.detect_contacts();
        let loads = sim.resolve_contacts(&contacts, DT);
        assert_eq!(loads[0].tangential, [0.0; 3]);
        assert_eq!(sim.tangential_springs[&(0, 1)], [0.0; 3]);
    }

    // ── GranularSimExtended ──────────────────────────────────────────────────

    #[test]
    fn extended_sim_friction_spins_the_particles_and_conserves_momentum() {
        let mut sim = sliding_extended_sim(
            spring_dashpot(5e4, 50.0),
            0.5,
            RollingResistanceModel::None,
            0.0,
            [1.0, 1.0],
        );
        let slip_before = extended_slip_speed(&sim);
        assert!((slip_before - 2.0).abs() < 1e-12);
        let angular_before = extended_angular_momentum(&sim);
        for _ in 0..20 {
            sim.step(DT);
        }
        // The lower sphere (slides −x, friction +x at its top) and the upper
        // sphere (slides +x, friction −x at its bottom) both spin about −z.
        for (i, w) in sim.angular_velocities.iter().enumerate() {
            assert!(w[2] < -0.1, "particle {i} spin {w:?}");
            assert!(w[0].abs() < 1e-12 && w[1].abs() < 1e-12);
        }
        let slip_after = extended_slip_speed(&sim);
        assert!(
            slip_after < 0.95 * slip_before,
            "sliding speed {slip_before} -> {slip_after}"
        );
        // The contact forces are internal: linear and angular momentum stay put.
        for k in 0..3 {
            let momentum: f64 = (0..2)
                .map(|i| sim.sim.masses[i] * sim.sim.velocities[i][k])
                .sum();
            assert!(momentum.abs() < 1e-12, "momentum[{k}] = {momentum}");
        }
        let angular_after = extended_angular_momentum(&sim);
        for k in 0..3 {
            assert!(
                (angular_after[k] - angular_before[k]).abs() < 1e-12,
                "angular momentum[{k}]: {} -> {}",
                angular_before[k],
                angular_after[k]
            );
        }
        assert!((angular_before[2] + 0.18).abs() < 1e-12);
    }

    #[test]
    fn extended_sim_hertzian_friction_spins_the_particles_with_the_mindlin_stiffness() {
        // Hertz–Mindlin pair (no damper): k_t = 8 G* √(R* δ) ≈ 1.0e5 N/m, so the
        // friction force grows as k_t v_t Δt per step in the stick regime and the
        // torques ℓ n × F_t spin both spheres about −z while the sliding decays.
        let model = ContactModel::Hertzian {
            e_star: 1e6,
            g_star: 4e5,
            r_star: 0.05,
        };
        let mut sim =
            sliding_extended_sim(model, 0.5, RollingResistanceModel::None, 0.0, [1.0, 1.0]);
        let slip_before = extended_slip_speed(&sim);
        let angular_before = extended_angular_momentum(&sim);
        for _ in 0..20 {
            sim.step(DT);
        }
        for (i, w) in sim.angular_velocities.iter().enumerate() {
            assert!(w[2] < -0.05, "particle {i} spin {w:?}");
            assert!(w[0].abs() < 1e-12 && w[1].abs() < 1e-12);
        }
        assert!(
            extended_slip_speed(&sim) < 0.995 * slip_before,
            "sliding speed {slip_before} -> {}",
            extended_slip_speed(&sim)
        );
        let angular_after = extended_angular_momentum(&sim);
        for k in 0..3 {
            assert!(
                (angular_after[k] - angular_before[k]).abs() < 1e-12,
                "angular momentum[{k}]: {} -> {}",
                angular_before[k],
                angular_after[k]
            );
        }
        // The 21st step of the stick regime (the initial slip speed 2 only decays
        // by about a percent): the spring has accumulated ≈ 21 v_t Δt, so the
        // friction force on particle 0 is ≈ k_t · 2 · 21 Δt ≈ 42 N along +x (against
        // its −x sliding), far below the Coulomb limit μ F_n ≈ 420 N.
        let contacts = sim.sim.detect_contacts();
        let loads = sim
            .sim
            .resolve_contacts_with_spin(&contacts, DT, &sim.angular_velocities);
        assert_eq!(loads.len(), 1);
        let kt = 8.0 * 4e5 * (0.05_f64 * 0.02).sqrt();
        let stick_force = kt * 2.0 * 21.0 * DT;
        let friction = norm(loads[0].tangential);
        assert!(
            friction > 0.95 * stick_force && friction < 1.001 * stick_force,
            "|F_t| = {friction}, stick-regime estimate {stick_force}"
        );
        assert!(loads[0].tangential[0] > 0.0, "{:?}", loads[0].tangential);
        assert!(friction < 0.2 * 0.5 * loads[0].normal);
    }

    #[test]
    fn extended_sim_rolling_resistance_uses_the_contact_normal_force() {
        // Twisting spin about the contact axis causes no sliding, so only the
        // rolling-resistance torque acts.  Gravity along the contact axis changes
        // the total force on a particle but not the contact's normal force.
        let model = spring_dashpot(5e4, 0.0);
        let mut sim = GranularSimExtended::new(
            model,
            0.3,
            [-100.0, 0.0, 0.0],
            RollingResistanceModel::Constant,
            0.01,
        );
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.18, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.angular_velocities[0] = [10.0, 0.0, 0.0];
        let dt = 1e-4;
        sim.step(dt);
        // F_n = kn δ = 2000 N, R_eff = 0.05, |τ| = μ_r R_eff F_n = 1 N·m,
        // I = 0.4 m r² = 0.004 kg m²: Δω = τ Δt / I = 0.025.
        let delta_omega = 1.0 / 0.004 * dt;
        assert!(
            (sim.angular_velocities[0][0] - (10.0 - delta_omega)).abs() < 1e-12,
            "ω₀ = {:?}",
            sim.angular_velocities[0]
        );
        assert!((sim.angular_velocities[1][0] - delta_omega).abs() < 1e-12);
        assert!(sim.angular_velocities[0][1].abs() < 1e-14);
        assert!(sim.angular_velocities[0][2].abs() < 1e-14);
        assert!(
            sim.sim.tangential_springs[&(0, 1)]
                .iter()
                .all(|&d| d.abs() < 1e-18)
        );
    }

    #[test]
    fn extended_sim_rolling_resistance_damps_the_relative_spin_of_a_sliding_pair() {
        // Different masses give the two spheres different spin-up rates, so the
        // relative spin ω₀ − ω₁ is nonzero and the rolling resistance acts on it.
        let run = |rolling: RollingResistanceModel, mu_r: f64| {
            let mut sim =
                sliding_extended_sim(spring_dashpot(5e4, 50.0), 0.5, rolling, mu_r, [1.0, 4.0]);
            let angular_before = extended_angular_momentum(&sim);
            for _ in 0..60 {
                sim.step(DT);
            }
            let angular_after = extended_angular_momentum(&sim);
            for k in 0..3 {
                assert!(
                    (angular_after[k] - angular_before[k]).abs() < 1e-12,
                    "angular momentum[{k}] not conserved with {rolling:?}"
                );
            }
            (sim.angular_velocities[0][2] - sim.angular_velocities[1][2]).abs()
        };
        let free = run(RollingResistanceModel::None, 0.0);
        let damped = run(RollingResistanceModel::Viscous, 0.1);
        assert!(free > 1.0, "relative spin without resistance {free}");
        assert!(
            damped < 0.9 * free,
            "rolling resistance must reduce the relative spin: {free} -> {damped}"
        );
    }

    // ── GranularSimulation (Hertz–Mindlin) and DemSimulation ─────────────────

    #[test]
    fn granular_simulation_friction_opposes_sliding_and_conserves_angular_momentum() {
        let contact = HertzContact::new(1e6, 4e5);
        let mut sim = GranularSimulation::new(contact, [0.0; 3], 0.3, 0.0, DT);
        for (pos, vel) in PAIR_POSITIONS.iter().zip(PAIR_VELOCITIES.iter()) {
            sim.add_particle(*pos, *vel, 0.1, 1.0);
        }
        let angular_momentum = |sim: &GranularSimulation| {
            let mut total = [0.0; 3];
            for p in &sim.particles {
                let orbital = cross(p.position, p.velocity);
                for (k, t) in total.iter_mut().enumerate() {
                    *t += p.mass * orbital[k] + p.moment_of_inertia * p.angular_velocity[k];
                }
            }
            total
        };
        let slip_speed = |sim: &GranularSimulation| {
            let (a, b) = (&sim.particles[0], &sim.particles[1]);
            let d = [
                b.position[0] - a.position[0],
                b.position[1] - a.position[1],
                b.position[2] - a.position[2],
            ];
            let dist = norm(d);
            let normal = [d[0] / dist, d[1] / dist, d[2] / dist];
            let lever = a.radius - 0.5 * (a.radius + b.radius - dist);
            let slip = contact_point_relative_velocity(
                a.velocity,
                b.velocity,
                a.angular_velocity,
                b.angular_velocity,
                lever,
                lever,
                normal,
            );
            norm(tangential_component(slip, normal))
        };
        let slip_before = slip_speed(&sim);
        let l_before = angular_momentum(&sim);
        assert!((slip_before - 2.0).abs() < 1e-12);
        sim.step(DT);
        // First step: F_t = k_t · v_t Δt = 8 G* √(R* δ) · 2Δt on particle 0, along +x.
        let kt = 8.0 * 4e5 * (0.05_f64 * 0.02).sqrt();
        let expected_dv = kt * 2.0 * DT * DT;
        assert!(
            (sim.particles[0].velocity[0] - (-1.0 + expected_dv)).abs() < 1e-12,
            "v₀ₓ = {}",
            sim.particles[0].velocity[0]
        );
        for _ in 0..99 {
            sim.step(DT);
        }
        for (i, p) in sim.particles.iter().enumerate() {
            assert!(
                p.angular_velocity[2] < -0.5,
                "particle {i} spin {:?}",
                p.angular_velocity
            );
        }
        assert!(
            slip_speed(&sim) < 0.9 * slip_before,
            "sliding speed {slip_before} -> {}",
            slip_speed(&sim)
        );
        let l_after = angular_momentum(&sim);
        for k in 0..3 {
            assert!(
                (l_after[k] - l_before[k]).abs() < 1e-12,
                "angular momentum[{k}]: {} -> {}",
                l_before[k],
                l_after[k]
            );
        }
        let momentum = sim.total_momentum();
        assert!(
            momentum.iter().all(|m| m.abs() < 1e-12),
            "momentum {momentum:?}"
        );
        assert!(!sim.tangential_springs.is_empty());
    }

    #[test]
    fn granular_simulation_drops_the_history_of_separated_pairs() {
        let contact = HertzContact::new(1e6, 4e5);
        let mut sim = GranularSimulation::new(contact, [0.0; 3], 0.3, 0.0, DT);
        for (pos, vel) in PAIR_POSITIONS.iter().zip(PAIR_VELOCITIES.iter()) {
            sim.add_particle(*pos, *vel, 0.1, 1.0);
        }
        sim.step(DT);
        assert_eq!(sim.tangential_springs.len(), 1);
        sim.particles[1].position = [0.0, 2.0, 0.0];
        sim.step(DT);
        assert!(sim.tangential_springs.is_empty());
    }

    #[test]
    fn dem_simulation_friction_opposes_sliding_and_conserves_angular_momentum() {
        let mut sim = sliding_dem_sim(5e4, 0.3);
        let angular_momentum = |sim: &DemSimulation| {
            let mut total = [0.0; 3];
            for p in &sim.particles {
                let orbital = cross(p.position, p.velocity);
                for (k, t) in total.iter_mut().enumerate() {
                    *t += p.mass * orbital[k] + p.moment_of_inertia * p.angular_velocity[k];
                }
            }
            total
        };
        let l_before = angular_momentum(&sim);
        sim.step();
        // First step: F_t = k_t · v_t Δt = 5e4 · 2Δt on particle 0, along +x.
        let expected_dv = 5e4 * 2.0 * DT * DT;
        assert!(
            (sim.particles[0].velocity[0] - (-1.0 + expected_dv)).abs() < 1e-12,
            "v₀ₓ = {} (friction must decelerate the sliding)",
            sim.particles[0].velocity[0]
        );
        assert!((sim.particles[1].velocity[0] - (1.0 - expected_dv)).abs() < 1e-12);
        for _ in 0..99 {
            sim.step();
        }
        for (i, p) in sim.particles.iter().enumerate() {
            assert!(
                p.angular_velocity[2] < -0.5,
                "particle {i} spin {:?}",
                p.angular_velocity
            );
            assert!(p.angular_velocity[0].abs() < 1e-12 && p.angular_velocity[1].abs() < 1e-12);
        }
        // The spheres also push each other apart along x (the contact normal tilts
        // while they slide past each other), so the friction-induced change of the
        // relative translational velocity is measured against a frictionless run.
        let relative_x = sim.particles[1].velocity[0] - sim.particles[0].velocity[0];
        let mut frictionless = sliding_dem_sim(5e4, 0.0);
        for _ in 0..100 {
            frictionless.step();
        }
        let reference_x =
            frictionless.particles[1].velocity[0] - frictionless.particles[0].velocity[0];
        let reduction = reference_x - relative_x;
        // Stick regime (|F_t| ≤ 600 N): the slip velocity oscillates with
        // ω₀² = k_t (2/m + 2 ℓ²/I), ℓ = r − δ/2 = 0.09, I = 0.4 m r² = 0.004, so the
        // friction impulse after t is ∫F dt = k_t v₀ (1 − cos ω₀ t) / ω₀², v₀ = 2, and
        // the relative velocity of the two unit masses drops by twice that.
        let omega0_sq = 5e4 * (2.0 + 2.0 * 0.09_f64.powi(2) / 0.004);
        let t = 100.0 * DT;
        let expected_reduction = 2.0 * 5e4 * 2.0 * (1.0 - (omega0_sq.sqrt() * t).cos()) / omega0_sq;
        assert!(
            reduction > 0.05,
            "friction must decelerate the translational sliding: {reference_x} -> {relative_x}"
        );
        assert!(
            (reduction - expected_reduction).abs() < 0.1 * expected_reduction,
            "friction impulse {reduction}, stick-regime closed form {expected_reduction}"
        );
        let l_after = angular_momentum(&sim);
        for k in 0..3 {
            assert!(
                (l_after[k] - l_before[k]).abs() < 1e-12,
                "angular momentum[{k}]: {} -> {}",
                l_before[k],
                l_after[k]
            );
        }
        let momentum: f64 = sim.particles.iter().map(|p| p.mass * p.velocity[0]).sum();
        assert!(momentum.abs() < 1e-12);
        assert_eq!(sim.tangential_springs.len(), 1);
    }

    #[test]
    fn dem_simulation_friction_is_capped_by_the_coulomb_cone() {
        // A very stiff tangential spring saturates at μ F_n = 0.3 · 2000 N: the
        // first step changes the sliding velocity by exactly μ F_n Δt / m.
        let mut sim = sliding_dem_sim(1e12, 0.3);
        sim.step();
        let expected = -1.0 + 0.3 * 2000.0 * DT;
        assert!(
            (sim.particles[0].velocity[0] - expected).abs() < 1e-12,
            "v₀ₓ = {}",
            sim.particles[0].velocity[0]
        );
    }
}
