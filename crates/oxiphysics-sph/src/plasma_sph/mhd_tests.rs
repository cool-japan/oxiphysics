// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the SPH operators and the time integration of
//! [`MhdSph`].
//!
//! The operators are evaluated at the centre of a cubic particle lattice of
//! spacing `s` with the smoothing length `h = 1.3 s` (kernel radius `2.6 s`), so
//! that the centre particle has a complete kernel support, and are compared with
//! the analytic derivatives of smooth fields.  The density is the SPH summation,
//! so `Σⱼ Vⱼ W(rᵢⱼ) = 1` and the estimates are accurate to a few 10⁻³.

use super::*;

/// Lattice spacing \[m\].
const SPACING: f64 = 0.01;
/// Smoothing length in lattice spacings.
const H_OVER_S: f64 = 1.3;
/// Rest density of the lattice \[kg m⁻³\].
const REST_DENSITY: f64 = 1.0;

/// A cold (`u = 0`) cubic lattice of `n³` particles at rest, centred on the
/// origin, with consistent summation densities.
fn lattice(n: usize) -> MhdSph {
    let mut solver = MhdSph::new(0.0, 5.0 / 3.0);
    let half = (n as f64 - 1.0) / 2.0;
    let mass = REST_DENSITY * SPACING.powi(3);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let position = [
                    (ix as f64 - half) * SPACING,
                    (iy as f64 - half) * SPACING,
                    (iz as f64 - half) * SPACING,
                ];
                let mut p =
                    PlasmaParticle::new(position, [0.0; 3], mass, E_CHARGE, 0.0, REST_DENSITY);
                p.smoothing_length = H_OVER_S * SPACING;
                solver.add_particle(p);
            }
        }
    }
    solver.update_thermodynamics();
    solver
}

/// Index of the particle closest to `target`.
fn nearest(solver: &MhdSph, target: [f64; 3]) -> usize {
    let mut best = 0;
    let mut best_d2 = f64::INFINITY;
    for (i, p) in solver.particles.iter().enumerate() {
        let d = sub3(p.position, target);
        let d2 = dot3(d, d);
        if d2 < best_d2 {
            best = i;
            best_d2 = d2;
        }
    }
    best
}

/// Index of the central particle of a lattice with an odd number of particles
/// per edge.
fn centre(solver: &MhdSph) -> usize {
    nearest(solver, [0.0; 3])
}

/// Set the magnetic field of every particle from its position.
fn set_field(solver: &mut MhdSph, field: impl Fn([f64; 3]) -> [f64; 3]) {
    for p in solver.particles.iter_mut() {
        p.magnetic_field = field(p.position);
    }
}

/// Set the velocity of every particle from its position.
fn set_velocity(solver: &mut MhdSph, velocity: impl Fn([f64; 3]) -> [f64; 3]) {
    for p in solver.particles.iter_mut() {
        p.velocity = velocity(p.position);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ∇×B and the Lorentz force
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn curl_of_b_equals_dby_dx_along_plus_z() {
    // B = (0, x, 0) has ∇×B = (0, 0, +1); the old operator returned −1.
    let mut solver = lattice(7);
    set_field(&mut solver, |x| [0.0, x[0], 0.0]);
    let curl = solver.curl_b(centre(&solver));
    assert!((curl[2] - 1.0).abs() < 0.01, "curl_z = {}", curl[2]);
    assert!(curl[0].abs() < 1e-9, "curl_x = {}", curl[0]);
    assert!(curl[1].abs() < 1e-9, "curl_y = {}", curl[1]);
}

#[test]
fn curl_of_b_is_the_true_curl_for_a_rotating_field() {
    // B = (−y, x, 0) has ∇×B = (0, 0, 2).
    let mut solver = lattice(7);
    set_field(&mut solver, |x| [-x[1], x[0], 0.0]);
    let curl = solver.curl_b(centre(&solver));
    assert!((curl[2] - 2.0).abs() < 0.02, "curl_z = {}", curl[2]);
    assert!(curl[0].abs() < 1e-9 && curl[1].abs() < 1e-9);
}

#[test]
fn curl_of_a_uniform_field_vanishes() {
    let mut solver = lattice(5);
    set_field(&mut solver, |_| [0.3, -0.2, 0.7]);
    for i in 0..solver.particles.len() {
        assert_eq!(solver.curl_b(i), [0.0; 3]);
    }
}

#[test]
fn lorentz_force_points_down_the_magnetic_pressure_gradient() {
    // B = (0, 1 + x, 0): the magnetic pressure B²/(2μ₀) grows with x, so the
    // plasma is pushed towards −x with a = −(dB_y/dx) B_y/(μ₀ρ).
    let mut solver = lattice(7);
    set_field(&mut solver, |x| [0.0, 1.0 + x[0], 0.0]);
    let c = centre(&solver);
    let a = solver.lorentz_force(c);
    let expected = -1.0 / (MU_0 * solver.particles[c].density);
    assert!(a[0] < 0.0, "a_x = {} must point towards −x", a[0]);
    assert!(
        (a[0] / expected - 1.0).abs() < 0.01,
        "a_x = {} vs {expected}",
        a[0]
    );
    assert!(a[1].abs() < 1e-6 * expected.abs(), "a_y = {}", a[1]);
    assert!(a[2].abs() < 1e-6 * expected.abs(), "a_z = {}", a[2]);
}

#[test]
fn lorentz_force_of_a_force_free_uniform_field_vanishes() {
    let mut solver = lattice(5);
    set_field(&mut solver, |_| [0.0, 0.0, 1.0]);
    for i in 0..solver.particles.len() {
        assert_eq!(solver.lorentz_force(i), [0.0; 3]);
    }
}

#[test]
fn operators_with_out_of_range_index_return_zero() {
    let mut solver = lattice(3);
    set_field(&mut solver, |x| [0.0, x[0], 0.0]);
    let missing = solver.particles.len() + 7;
    assert_eq!(solver.curl_b(missing), [0.0; 3]);
    assert_eq!(solver.divergence_b(missing), 0.0);
    assert_eq!(solver.induction_rhs(missing), [0.0; 3]);
    assert_eq!(solver.lorentz_force(missing), [0.0; 3]);
    assert_eq!(solver.pressure_acceleration(missing), [0.0; 3]);
    assert_eq!(solver.internal_energy_rate(missing), 0.0);
    assert_eq!(solver.resistive_diffusion(missing), [0.0; 3]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Induction equation
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn induction_of_uniform_flow_and_field_vanishes() {
    let mut solver = lattice(5);
    set_field(&mut solver, |_| [0.0, 1.0, 0.0]);
    set_velocity(&mut solver, |_| [3.0, -1.0, 2.0]);
    for i in 0..solver.particles.len() {
        assert_eq!(solver.induction_rhs(i), [0.0; 3]);
    }
}

#[test]
fn induction_winds_a_sheared_field_up_with_the_full_rate() {
    // v = (y, 0, 0), B = (0, 1, 0): dB_x/dt = (B·∇)v_x = +1.  The old estimate
    // was about −½ of the curl.
    let mut solver = lattice(7);
    set_velocity(&mut solver, |x| [x[1], 0.0, 0.0]);
    set_field(&mut solver, |_| [0.0, 1.0, 0.0]);
    let rate = solver.induction_rhs(centre(&solver));
    assert!((rate[0] - 1.0).abs() < 0.01, "dBx/dt = {}", rate[0]);
    assert!(rate[1].abs() < 1e-9, "dBy/dt = {}", rate[1]);
    assert!(rate[2].abs() < 1e-9, "dBz/dt = {}", rate[2]);
}

#[test]
fn induction_freezes_the_flux_of_a_perpendicular_field_in_a_compression() {
    // v = (α x, 0, 0) compresses the plasma along x (∇·v = α); a field B ⊥ x
    // is frozen into it: dB_y/dt = −α B_y.
    let alpha = 2.0;
    let mut solver = lattice(7);
    set_velocity(&mut solver, |x| [alpha * x[0], 0.0, 0.0]);
    set_field(&mut solver, |_| [0.0, 1.5, 0.0]);
    let rate = solver.induction_rhs(centre(&solver));
    let expected = -alpha * 1.5;
    assert!(
        (rate[1] / expected - 1.0).abs() < 0.01,
        "dBy/dt = {} vs {expected}",
        rate[1]
    );
    assert!(rate[0].abs() < 1e-9 && rate[2].abs() < 1e-9);
}

#[test]
fn induction_leaves_a_field_along_a_one_dimensional_compression_unchanged() {
    // For B ∥ v = (α x, 0, 0): (B·∇)v − B(∇·v) = αB − αB = 0.
    let mut solver = lattice(7);
    set_velocity(&mut solver, |x| [2.0 * x[0], 0.0, 0.0]);
    set_field(&mut solver, |_| [1.0, 0.0, 0.0]);
    let rate = solver.induction_rhs(centre(&solver));
    assert!(rate[0].abs() < 1e-9, "dBx/dt = {}", rate[0]);
}

#[test]
fn induction_is_the_lagrangian_rate_not_the_eulerian_one() {
    // v = (y, 0, 0), B = (0, 1 + c x, 0).  At (0, s, 0) the Eulerian
    // ∇×(v×B) = (1, −c s, 0), but the field moves with the particle: the
    // advection −(v·∇)B = (0, −c s, 0) must not be applied on top of the motion,
    // so dB_y/dt ≈ 0 and dB_x/dt ≈ B_y = 1.
    let c = 40.0;
    let mut solver = lattice(7);
    set_velocity(&mut solver, |x| [x[1], 0.0, 0.0]);
    set_field(&mut solver, |x| [0.0, 1.0 + c * x[0], 0.0]);
    let i = nearest(&solver, [0.0, SPACING, 0.0]);
    let rate = solver.induction_rhs(i);
    let eulerian_by = -c * SPACING;
    assert!(rate[1].abs() < 0.05, "dBy/dt = {}", rate[1]);
    assert!(
        (rate[1] - eulerian_by).abs() > 0.3,
        "dBy/dt = {} must differ from the Eulerian {eulerian_by}",
        rate[1]
    );
    assert!((rate[0] - 1.0).abs() < 0.05, "dBx/dt = {}", rate[0]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Resistive diffusion
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn resistive_diffusion_of_a_quadratic_field_is_the_magnetic_diffusivity_times_two() {
    // B_x = x² has ∇²B_x = 2; the diffusivity is η/μ₀ = 2.5 m² s⁻¹.
    let mut solver = lattice(7);
    solver.resistivity = 2.5 * MU_0;
    set_field(&mut solver, |x| [x[0] * x[0], 0.0, 0.0]);
    let rate = solver.resistive_diffusion(centre(&solver));
    assert!((rate[0] / 5.0 - 1.0).abs() < 0.03, "rate_x = {}", rate[0]);
    assert!(rate[1].abs() < 1e-9 && rate[2].abs() < 1e-9);
}

#[test]
fn resistive_diffusion_of_a_linear_field_vanishes() {
    let mut solver = lattice(7);
    solver.resistivity = 3.0 * MU_0;
    set_field(&mut solver, |x| [x[0], 2.0 * x[1], -x[2]]);
    let rate = solver.resistive_diffusion(centre(&solver));
    for (k, r) in rate.iter().enumerate() {
        assert!(r.abs() < 1e-9, "component {k} = {r}");
    }
}

#[test]
fn resistive_diffusion_flattens_a_field_peak() {
    // A Gaussian peak has a negative Laplacian at its maximum.
    let mut solver = lattice(7);
    solver.resistivity = MU_0;
    let width = 2.0 * SPACING;
    set_field(&mut solver, |x| {
        let r2 = dot3(x, x);
        [(-r2 / (2.0 * width * width)).exp(), 0.0, 0.0]
    });
    let rate = solver.resistive_diffusion(centre(&solver));
    // ∇²exp(−r²/2w²) = −3/w² at the centre.
    let expected = -3.0 / (width * width);
    assert!(rate[0] < 0.0, "peak must diffuse, got {}", rate[0]);
    assert!(
        (rate[0] / expected - 1.0).abs() < 0.1,
        "rate_x = {} vs {expected}",
        rate[0]
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Gas pressure and the energy equation
// ─────────────────────────────────────────────────────────────────────────────

/// A lattice whose gas is hot in a Gaussian core: `u = 10³ (1 + 4 e^{−r²/2w²})`.
fn hot_core_lattice(n: usize) -> MhdSph {
    let mut solver = lattice(n);
    let width = 2.0 * SPACING;
    for p in solver.particles.iter_mut() {
        let r2 = dot3(p.position, p.position);
        p.internal_energy = 1.0e3 * (1.0 + 4.0 * (-r2 / (2.0 * width * width)).exp());
    }
    solver.update_thermodynamics();
    solver
}

/// Sound-crossing based time step of the hot-core lattice.
fn hot_core_dt(solver: &MhdSph) -> f64 {
    let u_max = solver
        .particles
        .iter()
        .map(|p| p.internal_energy)
        .fold(0.0_f64, f64::max);
    let c_s = (solver.gamma * (solver.gamma - 1.0) * u_max).sqrt();
    0.1 * SPACING / c_s
}

#[test]
fn add_particle_derives_the_internal_energy_from_the_temperature() {
    let gamma = 1.4;
    let mut solver = MhdSph::new(0.0, gamma);
    let p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1.0e6, 3.0e-9);
    let pressure = p.pressure;
    solver.add_particle(p);
    let q = &solver.particles[0];
    assert!(q.internal_energy > 0.0);
    let p_eos = (gamma - 1.0) * q.density * q.internal_energy;
    assert!(
        ((p_eos - pressure) / pressure).abs() < 1e-12,
        "(γ−1)ρu = {p_eos} vs reported pressure {pressure}"
    );
}

#[test]
fn add_particle_keeps_an_explicit_internal_energy() {
    let mut solver = MhdSph::new(0.0, 5.0 / 3.0);
    let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1.0e6, 1.0e-9);
    p.internal_energy = 123.0;
    solver.add_particle(p);
    assert_eq!(solver.particles[0].internal_energy, 123.0);
}

#[test]
fn pressure_follows_the_ideal_gas_law_with_the_solver_gamma() {
    let mut a = hot_core_lattice(5);
    let mut b = hot_core_lattice(5);
    b.gamma = 1.4;
    a.update_thermodynamics();
    b.update_thermodynamics();
    let c = centre(&a);
    let (pa, pb) = (a.particles[c].pressure, b.particles[c].pressure);
    let ratio = pb / pa;
    assert!(
        (ratio - 0.4 / (2.0 / 3.0)).abs() < 1e-12,
        "P(γ=1.4)/P(γ=5/3) = {ratio}"
    );
    for p in &a.particles {
        let expected = (a.gamma - 1.0) * p.density * p.internal_energy;
        assert!((p.pressure - expected).abs() <= 1e-12 * expected);
    }
}

#[test]
fn summation_density_of_the_lattice_interior_is_the_rest_density() {
    let solver = lattice(7);
    let rho = solver.particles[centre(&solver)].density;
    assert!((rho / REST_DENSITY - 1.0).abs() < 0.01, "ρ = {rho}");
    let corner = solver.particles[0].density;
    assert!(
        corner < 0.6 * rho,
        "corner ρ = {corner} lacks neighbours (centre ρ = {rho})"
    );
}

#[test]
fn gas_pressure_pushes_a_hot_core_outwards() {
    let mut solver = hot_core_lattice(7);
    let dt = hot_core_dt(&solver);
    solver.euler_step(dt);
    let mut checked = 0;
    for p in &solver.particles {
        let r = len3(p.position);
        if r < 0.5 * SPACING {
            continue;
        }
        let radial = dot3(p.velocity, scale3(p.position, 1.0 / r));
        assert!(
            radial > 0.05,
            "particle at r = {r} moves inwards or not at all: v_r = {radial}"
        );
        checked += 1;
    }
    assert_eq!(checked, 7 * 7 * 7 - 1);
    assert!(solver.total_kinetic_energy() > 0.0);
    // The step's acceleration is stored; by symmetry it vanishes at the centre.
    let c = centre(&solver);
    let a_max = solver
        .particles
        .iter()
        .map(|p| len3(p.accel))
        .fold(0.0_f64, f64::max);
    assert!(a_max > 1.0e3, "stored accelerations are missing: {a_max}");
    assert!(
        len3(solver.particles[c].accel) < 1e-9 * a_max,
        "centre accel = {:?}",
        solver.particles[c].accel
    );
}

#[test]
fn pure_gas_conserves_the_total_energy() {
    let mut solver = hot_core_lattice(7);
    let dt = hot_core_dt(&solver);
    let energy = |s: &MhdSph| s.total_kinetic_energy() + s.total_thermal_energy();
    let e0 = energy(&solver);
    for _ in 0..20 {
        solver.euler_step(dt);
    }
    let e1 = energy(&solver);
    assert!(
        ((e1 - e0) / e0).abs() < 1e-10,
        "energy drifted from {e0} to {e1}"
    );
    // The gas has really exchanged energy: a good part of the thermal energy
    // has become kinetic.
    assert!(
        solver.total_kinetic_energy() > 0.05 * e0,
        "kinetic energy {} stays tiny",
        solver.total_kinetic_energy()
    );
}

#[test]
fn expansion_cools_the_hot_core() {
    let mut solver = hot_core_lattice(7);
    let dt = hot_core_dt(&solver);
    let c = centre(&solver);
    let u_core = solver.particles[c].internal_energy;
    for _ in 0..10 {
        solver.euler_step(dt);
    }
    assert!(
        solver.particles[c].internal_energy < u_core,
        "the expanding core must cool: {} vs {u_core}",
        solver.particles[c].internal_energy
    );
}

#[test]
fn euler_step_refreshes_density_pressure_and_temperature() {
    let mut solver = hot_core_lattice(5);
    for p in solver.particles.iter_mut() {
        p.density = 7.0;
        p.pressure = -1.0;
    }
    solver.euler_step(1e-9);
    for p in &solver.particles {
        assert!(
            p.density > 0.0 && p.density < 1.5 * REST_DENSITY,
            "ρ = {}",
            p.density
        );
        let expected = (solver.gamma - 1.0) * p.density * p.internal_energy;
        assert!((p.pressure - expected).abs() <= 1e-12 * expected);
        let temperature = (solver.gamma - 1.0) * p.internal_energy * M_PROTON / K_B;
        assert!((p.temperature - temperature).abs() <= 1e-12 * temperature);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Time integration with a magnetic field
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn euler_step_accelerates_the_plasma_with_the_lorentz_force() {
    let mut solver = lattice(7);
    set_field(&mut solver, |x| [0.0, 1.0 + x[0], 0.0]);
    let c = centre(&solver);
    let a_expected = -1.0 / (MU_0 * solver.particles[c].density);
    let dt = 1.0e-7;
    solver.euler_step(dt);
    let v = solver.particles[c].velocity;
    assert!(
        (v[0] / (a_expected * dt) - 1.0).abs() < 0.01,
        "v_x = {} vs {}",
        v[0],
        a_expected * dt
    );
    assert!(v[1].abs() < 1e-6 * v[0].abs() && v[2].abs() < 1e-6 * v[0].abs());
    assert!((solver.particles[c].accel[0] / a_expected - 1.0).abs() < 0.01);
}

#[test]
fn euler_step_winds_up_a_sheared_field() {
    let mut solver = lattice(7);
    set_velocity(&mut solver, |x| [x[1], 0.0, 0.0]);
    set_field(&mut solver, |_| [0.0, 1.0, 0.0]);
    let dt = 1.0e-3;
    solver.euler_step(dt);
    let c = centre(&solver);
    // The centre particle moves with the shear; its field gained dt·dBx/dt.
    let b = solver.particles[c].magnetic_field;
    assert!((b[0] / dt - 1.0).abs() < 0.01, "Bx/dt = {}", b[0] / dt);
    assert!((b[1] - 1.0).abs() < 1e-9, "By = {}", b[1]);
}

#[test]
fn euler_step_keeps_a_uniform_field_in_a_uniform_flow() {
    let mut solver = lattice(5);
    set_velocity(&mut solver, |_| [3.0, 0.0, 0.0]);
    set_field(&mut solver, |_| [0.0, 1.0, 0.0]);
    let before: Vec<_> = solver.particles.iter().map(|p| p.position).collect();
    let dt = 1.0e-3;
    solver.euler_step(dt);
    for (p, x0) in solver.particles.iter().zip(before) {
        assert_eq!(p.magnetic_field, [0.0, 1.0, 0.0]);
        assert_eq!(p.psi, 0.0);
        assert!((p.position[0] - x0[0] - 3.0 * dt).abs() < 1e-15);
        assert!((p.position[1] - x0[1]).abs() < 1e-15);
    }
}

#[test]
fn euler_step_of_a_single_particle_stays_finite() {
    let mut solver = MhdSph::new(1e-7, 5.0 / 3.0);
    let mut p = PlasmaParticle::new([0.0; 3], [1.0, 0.0, 0.0], M_PROTON, E_CHARGE, 1e4, 1e-6);
    p.magnetic_field = [0.0, 0.0, 1e-3];
    solver.add_particle(p);
    solver.euler_step(0.01);
    let q = &solver.particles[0];
    assert!((q.position[0] - 0.01).abs() < 1e-15);
    assert_eq!(q.velocity, [1.0, 0.0, 0.0]);
    assert_eq!(q.magnetic_field, [0.0, 0.0, 1e-3]);
    assert!(q.density.is_finite() && q.density > 0.0);
    assert!(q.internal_energy.is_finite() && q.internal_energy > 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Dedner divergence cleaning
// ─────────────────────────────────────────────────────────────────────────────

/// A 9³ lattice carrying the radial field `B = 10⁻³ (x/s) e^{−r²/(2 (2s)²)}`,
/// whose divergence is concentrated near the centre, and the indices of the
/// particles in the inner half of the block (complete kernel support).
fn radial_field_lattice() -> (MhdSph, Vec<usize>) {
    let mut solver = lattice(9);
    set_field(&mut solver, |x| {
        let w = (-dot3(x, x) / (2.0 * (2.0 * SPACING).powi(2))).exp();
        scale3(x, 1.0e-3 / SPACING * w)
    });
    let inner: Vec<usize> = solver
        .particles
        .iter()
        .enumerate()
        .filter(|(_, p)| len3(p.position) < 2.0 * SPACING)
        .map(|(i, _)| i)
        .collect();
    (solver, inner)
}

fn max_abs_div_b(solver: &MhdSph, indices: &[usize]) -> f64 {
    indices
        .iter()
        .map(|&i| solver.divergence_b(i).abs())
        .fold(0.0_f64, f64::max)
}

#[test]
fn divergence_of_a_linear_field_is_its_trace() {
    // B = (x, 2y, −3z) has ∇·B = 0.
    let mut solver = lattice(7);
    set_field(&mut solver, |x| [x[0], 2.0 * x[1], -3.0 * x[2]]);
    let div = solver.divergence_b(centre(&solver));
    assert!(div.abs() < 1e-3, "∇·B = {div}");
    set_field(&mut solver, |x| [x[0], 2.0 * x[1], x[2]]);
    let div = solver.divergence_b(centre(&solver));
    assert!((div - 4.0).abs() < 0.04, "∇·B = {div}");
}

#[test]
fn dedner_cleaning_sources_psi_from_the_divergence() {
    // ψ starts at zero; the old step never wrote to it.
    let (mut solver, _) = radial_field_lattice();
    solver.dedner_ch = 1.0e3;
    let dt = 0.1 * H_OVER_S * SPACING / solver.dedner_ch;
    assert!(solver.particles.iter().all(|p| p.psi == 0.0));
    solver.dedner_cleaning_step(dt);
    let psi_max = solver
        .particles
        .iter()
        .map(|p| p.psi.abs())
        .fold(0.0_f64, f64::max);
    assert!(psi_max > 1e-3, "ψ was not sourced: max |ψ| = {psi_max}");
    // ψ ≈ −c_h² ∇·B dt: the sign follows the divergence (positive near the
    // centre for the outward field).
    let c = centre(&solver);
    assert!(
        solver.particles[c].psi < 0.0,
        "ψ = {}",
        solver.particles[c].psi
    );
}

#[test]
fn dedner_cleaning_reduces_the_divergence_of_the_field() {
    let (mut solver, inner) = radial_field_lattice();
    solver.dedner_ch = 1.0e3;
    let dt = 0.1 * H_OVER_S * SPACING / solver.dedner_ch;
    let initial = max_abs_div_b(&solver, &inner);
    assert!(initial > 0.1, "test field has no divergence: {initial}");
    for _ in 0..50 {
        solver.dedner_cleaning_step(dt);
    }
    let cleaned = max_abs_div_b(&solver, &inner);
    assert!(
        cleaned < 0.15 * initial,
        "max |∇·B| only went from {initial} to {cleaned}"
    );
    assert!(solver.particles.iter().all(|p| p.psi.is_finite()));
    assert!(
        solver
            .particles
            .iter()
            .all(|p| len3(p.magnetic_field).is_finite())
    );
}

#[test]
fn dedner_cleaning_does_not_rescale_the_field_of_an_isolated_particle() {
    // The old step shifted |B| by ψ dt c_h along B̂ (here it even reversed B).
    let mut solver = MhdSph::new(0.0, 5.0 / 3.0);
    let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-6);
    p.psi = 1.0;
    p.magnetic_field = [0.0, 0.0, 1e-3];
    solver.add_particle(p);
    solver.dedner_cleaning_step(1e-6);
    assert_eq!(solver.particles[0].magnetic_field, [0.0, 0.0, 1e-3]);
}

#[test]
fn dedner_psi_damping_is_stable_for_a_step_longer_than_the_damping_time() {
    // τ = h/(σ c_h) = 1.25e-7 s; a step of 1e-3 s is 8000 τ.
    let mut solver = MhdSph::new(0.0, 5.0 / 3.0);
    let mut p = PlasmaParticle::new([0.0; 3], [0.0; 3], M_PROTON, E_CHARGE, 1e4, 1e-6);
    p.psi = 1.0;
    solver.add_particle(p);
    solver.dedner_cleaning_step(1.0e-3);
    let psi = solver.particles[0].psi;
    assert!(psi.abs() < 1e-12, "ψ = {psi}");
}

#[test]
fn dedner_cleaning_ignores_a_non_positive_time_step() {
    let (mut solver, _) = radial_field_lattice();
    let before: Vec<_> = solver.particles.iter().map(|p| p.magnetic_field).collect();
    solver.dedner_cleaning_step(0.0);
    solver.dedner_cleaning_step(-1.0);
    solver.dedner_cleaning_step(f64::NAN);
    for (p, b) in solver.particles.iter().zip(before) {
        assert_eq!(p.magnetic_field, b);
        assert_eq!(p.psi, 0.0);
    }
}

#[test]
fn euler_step_runs_the_cleaning_with_the_updated_field() {
    // The induction update of a compressed, sheared field creates ∇·B ≠ 0 at the
    // block edge; after a few steps ψ is non-zero.
    let mut solver = lattice(7);
    solver.dedner_ch = 1.0e3;
    set_velocity(&mut solver, |x| [x[1], 0.0, 0.0]);
    set_field(&mut solver, |_| [0.0, 1.0, 0.0]);
    let dt = 0.05 * H_OVER_S * SPACING / solver.dedner_ch;
    for _ in 0..3 {
        solver.euler_step(dt);
    }
    let psi_max = solver
        .particles
        .iter()
        .map(|p| p.psi.abs())
        .fold(0.0_f64, f64::max);
    assert!(psi_max > 0.0, "cleaning never saw a divergence");
}
