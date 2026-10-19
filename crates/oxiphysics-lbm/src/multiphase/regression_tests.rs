//! Regression tests for the Rayleigh-Plesset bubble model and the
//! multi-component Shan-Chen force.

use super::types::{BubbleState, MultiComponentSC};

/// Conserved energy of the inviscid Rayleigh-Plesset equation with an
/// isothermal gas core:
/// `2*pi*rho*R^3*R_t^2 - 4*pi*P_gas0*R_eq^3*ln(R) + 4*pi/3*P_inf*R^3 + 4*pi*sigma*R^2`.
fn rayleigh_plesset_energy(b: &BubbleState) -> f64 {
    let p_gas0 = b.p_inf + 2.0 * b.sigma / b.r_eq;
    let pi = std::f64::consts::PI;
    2.0 * pi * b.rho_liquid * b.radius.powi(3) * b.radius_rate * b.radius_rate
        - 4.0 * pi * p_gas0 * b.r_eq.powi(3) * b.radius.ln()
        + 4.0 * pi / 3.0 * b.p_inf * b.radius.powi(3)
        + 4.0 * pi * b.sigma * b.radius * b.radius
}

/// Deterministic pseudo-random field in `[0.2, 1.0)` (64-bit LCG).
fn pseudo_random_field(seed: u64, n: usize) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            0.2 + 0.8 * ((state >> 11) as f64 / (1u64 << 53) as f64)
        })
        .collect()
}

#[test]
fn rayleigh_plesset_acceleration_matches_hand_evaluation() {
    // R_eq = 1 mm, R = 2 mm, R_t = 0.5 m/s, water-like liquid.
    //   p_gas0 = 1e5 + 2*0.072/1e-3 = 100144 Pa, p_gas = p_gas0 / 8 = 12518 Pa
    //   p_wall = 12518 - 1e5 - 2*0.072/2e-3 - 4*1e-3*0.5/2e-3 = -87555 Pa
    //   R_tt   = (p_wall / rho - 1.5 * R_t^2) / R = (-87.555 - 0.375) / 2e-3
    let mut b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 1e-3);
    b.radius = 2e-3;
    b.radius_rate = 0.5;
    let expected = -43965.0;
    let a = b.acceleration();
    assert!(
        ((a - expected) / expected).abs() < 1e-12,
        "R_tt = {a}, expected {expected}"
    );
}

#[test]
fn rayleigh_plesset_inertia_term_is_not_divided_by_density() {
    // At R = R_eq with sigma = mu = 0 the pressure balance vanishes, so the
    // acceleration is purely the inertia term -1.5 * R_t^2 / R.
    let r_eq = 1e-3;
    let rate = 2.0;
    let expected = -1.5 * rate * rate / r_eq;
    for &rho in &[1.0, 1000.0, 13_600.0] {
        let mut b = BubbleState::new(r_eq, 1e5, rho, 0.0, 0.0);
        b.radius_rate = rate;
        let a = b.acceleration();
        assert!(
            ((a - expected) / expected).abs() < 1e-12,
            "rho = {rho}: R_tt = {a}, expected {expected}"
        );
    }
}

#[test]
fn rayleigh_plesset_viscous_term_is_four_mu_rate_over_rho_r2() {
    let mut inviscid = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 0.0);
    let mut viscous = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 1e-3);
    for b in [&mut inviscid, &mut viscous] {
        b.radius = 2e-3;
        b.radius_rate = 0.5;
    }
    // -4 mu R_t / (rho R^2) = -4e-3 * 0.5 / (1000 * 4e-6) = -0.5
    let delta = viscous.acceleration() - inviscid.acceleration();
    assert!((delta + 0.5).abs() < 1e-9, "viscous contribution = {delta}");
}

#[test]
fn rayleigh_plesset_rk4_conserves_energy_of_inviscid_bubble() {
    let mut b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 0.0);
    b.radius = 1.5e-3;
    let e0 = rayleigh_plesset_energy(&b);
    let mut max_drift = 0.0_f64;
    let mut min_radius = b.radius;
    // ~1.1 oscillation periods of a large-amplitude oscillation.
    for _ in 0..40_000 {
        b.step(1e-8);
        max_drift = max_drift.max(((rayleigh_plesset_energy(&b) - e0) / e0).abs());
        min_radius = min_radius.min(b.radius);
    }
    assert!(
        min_radius < 0.8e-3,
        "bubble must actually collapse below equilibrium, min R = {min_radius}"
    );
    assert!(max_drift < 1e-11, "relative energy drift = {max_drift}");
}

#[test]
fn rayleigh_plesset_step_is_fourth_order() {
    let r0 = 1.3e-3;
    let horizon = 1.5e-4;
    let integrate = |n: usize| {
        let mut b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 1e-3);
        b.radius = r0;
        let dt = horizon / n as f64;
        for _ in 0..n {
            b.step(dt);
        }
        b.radius
    };
    let reference = integrate(40_000);
    let errors: Vec<f64> = [50usize, 100, 200]
        .iter()
        .map(|&n| (integrate(n) - reference).abs())
        .collect();
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (12.0..20.0).contains(&ratio),
            "halving dt must reduce the error ~16x, got {ratio} (errors {errors:?})"
        );
    }
}

#[test]
fn rayleigh_plesset_small_oscillation_half_period_matches_linear_theory() {
    // Linearised isothermal bubble with surface tension:
    //   omega0^2 = (3 * p_gas0 - 2 * sigma / R_eq) / (rho * R_eq^2)
    let (r_eq, p_inf, rho, sigma): (f64, f64, f64, f64) = (1e-3, 1e5, 1000.0, 0.072);
    let p_gas0 = p_inf + 2.0 * sigma / r_eq;
    let omega0 = ((3.0 * p_gas0 - 2.0 * sigma / r_eq) / (rho * r_eq * r_eq)).sqrt();
    let half_period = std::f64::consts::PI / omega0;
    let mut b = BubbleState::new(r_eq, p_inf, rho, sigma, 0.0);
    b.radius = r_eq * (1.0 + 1e-4);
    let dt = 1e-9;
    let mut t = 0.0;
    let mut crossing = None;
    for _ in 0..400_000 {
        let rate_before = b.radius_rate;
        b.step(dt);
        t += dt;
        if rate_before < 0.0 && b.radius_rate >= 0.0 {
            let frac = rate_before / (rate_before - b.radius_rate);
            crossing = Some(t - dt + frac * dt);
            break;
        }
    }
    let measured = crossing.expect("the bubble must turn around within one period");
    assert!(
        ((measured - half_period) / half_period).abs() < 5e-4,
        "measured half period {measured}, linear theory {half_period}"
    );
}

#[test]
fn rayleigh_plesset_viscosity_dissipates_energy_and_damps_amplitude() {
    let mut b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 0.5);
    b.radius = 1.2e-3;
    let dt = 1e-8;
    let period_steps = 36_260;
    let mut previous_energy = rayleigh_plesset_energy(&b);
    let energy_scale = previous_energy.abs();
    let mut peaks = [0.0_f64; 2];
    for (k, peak) in peaks.iter_mut().enumerate() {
        let mut max_radius = 0.0_f64;
        for _ in 0..period_steps {
            b.step(dt);
            max_radius = max_radius.max(b.radius);
            let energy = rayleigh_plesset_energy(&b);
            assert!(
                energy <= previous_energy + 1e-13 * energy_scale,
                "viscous energy must not increase (period {k}): {previous_energy} -> {energy}"
            );
            previous_energy = energy;
        }
        *peak = max_radius;
    }
    let first_excursion = peaks[0] - 1e-3;
    let second_excursion = peaks[1] - 1e-3;
    assert!(
        second_excursion < 0.85 * first_excursion && second_excursion > 0.4 * first_excursion,
        "excursions {first_excursion} -> {second_excursion}"
    );
}

#[test]
fn multicomponent_sc_like_force_on_a_uses_psi_a() {
    let (nx, ny) = (8usize, 8usize);
    let mc = MultiComponentSC::new(-0.7, -0.9, -1.5, 1.0);
    // Smooth field with a non-zero gradient in both directions at (x, y).
    let psi_a: Vec<f64> = (0..nx * ny)
        .map(|k| {
            let phase_x = std::f64::consts::TAU * (k % nx) as f64 / nx as f64;
            let phase_y = std::f64::consts::TAU * (k / nx) as f64 / ny as f64;
            1.0 + 0.3 * phase_x.sin() + 0.2 * (phase_y + 0.5).cos()
        })
        .collect();
    let psi_b = vec![0.0_f64; nx * ny];
    let (x, y) = (1usize, 3usize);
    // Independent 3x3 stencil evaluation of  -g_aa * psi_a(x) * sum w psi_a(x+c) c.
    let weight = |dx: i64, dy: i64| match (dx != 0, dy != 0) {
        (false, false) => 0.0,
        (true, true) => 1.0 / 36.0,
        _ => 1.0 / 9.0,
    };
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    for dy in -1i64..=1 {
        for dx in -1i64..=1 {
            let xi = (x as i64 + dx).rem_euclid(nx as i64) as usize;
            let yi = (y as i64 + dy).rem_euclid(ny as i64) as usize;
            sum_x += weight(dx, dy) * psi_a[yi * nx + xi] * dx as f64;
            sum_y += weight(dx, dy) * psi_a[yi * nx + xi] * dy as f64;
        }
    }
    let psi_here = psi_a[y * nx + x];
    let expected = (-mc.g_aa * psi_here * sum_x, -mc.g_aa * psi_here * sum_y);
    let (fx, fy) = mc.interaction_force(&psi_a, &psi_b, x, y, nx, ny);
    assert!(
        expected.0.abs() > 1e-3 && expected.1.abs() > 1e-3,
        "test needs a non-trivial gradient in both directions: {expected:?}"
    );
    assert!(
        (fx - expected.0).abs() < 1e-14,
        "fx = {fx}, expected {}",
        expected.0
    );
    assert!(
        (fy - expected.1).abs() < 1e-14,
        "fy = {fy}, expected {}",
        expected.1
    );
    // With psi_B == 0 component B feels nothing.
    let (bx, by) = mc.interaction_force_b(&psi_a, &psi_b, x, y, nx, ny);
    assert_eq!((bx, by), (0.0, 0.0));
}

#[test]
fn multicomponent_sc_like_force_on_b_uses_g_bb_and_psi_b() {
    let (nx, ny) = (8usize, 6usize);
    let psi_a = vec![0.0_f64; nx * ny];
    let psi_b: Vec<f64> = (0..nx * ny)
        .map(|k| 0.5 + 0.2 * ((k % nx) as f64 * 0.7).cos() + 0.1 * (k / nx) as f64)
        .collect();
    let with_gbb = MultiComponentSC::new(0.0, -0.9, -1.5, 1.0);
    let without_gbb = MultiComponentSC::new(0.0, 0.0, -1.5, 1.0);
    let (x, y) = (3usize, 2usize);
    let (fx, fy) = with_gbb.interaction_force_b(&psi_a, &psi_b, x, y, nx, ny);
    let (f0x, f0y) = without_gbb.interaction_force_b(&psi_a, &psi_b, x, y, nx, ny);
    assert!(
        f0x == 0.0 && f0y == 0.0,
        "g_bb = 0 and psi_A = 0 give no force"
    );
    assert!(fx.abs() + fy.abs() > 1e-3, "g_bb must act on component B");
    // Mirror symmetry: F_B(psi_A, psi_B; g_aa, g_bb) = F_A(psi_B, psi_A; g_bb, g_aa).
    let mirror = MultiComponentSC::new(-0.9, 0.0, -1.5, 1.0);
    let (mx, my) = mirror.interaction_force(&psi_b, &psi_a, x, y, nx, ny);
    assert!((fx - mx).abs() < 1e-15 && (fy - my).abs() < 1e-15);
}

#[test]
fn multicomponent_sc_forces_conserve_total_momentum() {
    let (nx, ny) = (7usize, 6usize);
    let psi_a = pseudo_random_field(11, nx * ny);
    let psi_b = pseudo_random_field(23, nx * ny);
    let models = [
        MultiComponentSC::new(-0.8, -0.4, -1.3, 1.0),
        MultiComponentSC::new(-0.8, -0.4, 0.0, 1.0),
        MultiComponentSC::new(0.0, 0.0, -1.3, 1.0),
    ];
    for mc in &models {
        let mut total_a = (0.0, 0.0);
        let mut total_b = (0.0, 0.0);
        for y in 0..ny {
            for x in 0..nx {
                let ((ax, ay), (bx, by)) = mc.interaction_forces(&psi_a, &psi_b, x, y, nx, ny);
                total_a = (total_a.0 + ax, total_a.1 + ay);
                total_b = (total_b.0 + bx, total_b.1 + by);
            }
        }
        // Cross forces are equal and opposite, self forces sum to zero within a component.
        assert!(
            (total_a.0 + total_b.0).abs() < 1e-12 && (total_a.1 + total_b.1).abs() < 1e-12,
            "net force {:?} for {mc:?}",
            (total_a.0 + total_b.0, total_a.1 + total_b.1)
        );
        if mc.g_cc == 0.0 {
            assert!(
                total_a.0.abs() < 1e-12 && total_a.1.abs() < 1e-12,
                "self-interaction of A must conserve momentum: {total_a:?}"
            );
            assert!(
                total_b.0.abs() < 1e-12 && total_b.1.abs() < 1e-12,
                "self-interaction of B must conserve momentum: {total_b:?}"
            );
        }
    }
}
