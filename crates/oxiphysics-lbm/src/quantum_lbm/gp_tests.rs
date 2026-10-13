//! Physics regression tests of the Gross-Pitaevskii split-step solver.
//!
//! The reference solutions are exact: the free-particle dispersion relation,
//! the phase of a uniform condensate (`theta = -g rho t / hbar`), the
//! stationary bright soliton of the attractive 1-D GP equation and the ground
//! state of a harmonic trap.

use super::*;

/// Positions `x_j = (j - nx/2) dx` centred on the origin.
fn centred_positions(nx: usize, dx: f64) -> Vec<f64> {
    (0..nx).map(|j| (j as f64 - (nx / 2) as f64) * dx).collect()
}

fn solver(nx: usize, hbar: f64, mass: f64, coupling: f64, dx: f64, dt: f64) -> GrossPitaevskiiLbm {
    GrossPitaevskiiLbm::new(nx, QuantumLbmParams::new(hbar, mass, 1.0, coupling), dx, dt)
}

#[test]
fn particle_number_is_conserved_to_rounding_for_nonlinear_trapped_state() {
    let nx = 64;
    let dx = 0.5;
    let length = nx as f64 * dx;
    let mut gp = solver(nx, 1.0, 1.0, 0.8, dx, 0.01);
    let x = centred_positions(nx, dx);
    gp.potential = x.iter().map(|&xi| 0.5 * 0.09 * xi * xi).collect();
    let amplitude: Vec<f64> = (0..nx)
        .map(|j| {
            let s = j as f64 * dx / length;
            1.0 + 0.3 * (2.0 * PI * s).cos() + 0.2 * (4.0 * PI * s).sin()
        })
        .collect();
    let phase: Vec<f64> = (0..nx)
        .map(|j| 0.7 * (2.0 * PI * j as f64 * dx / length).sin())
        .collect();
    gp.set_wavefunction(amplitude, phase);
    let n0 = gp.total_number();
    for _ in 0..1000 {
        gp.step();
    }
    let n1 = gp.total_number();
    assert!(
        ((n1 - n0) / n0).abs() < 1e-10,
        "particle number drifted: {n0} -> {n1}"
    );
    assert!(gp.amplitude.iter().all(|a| a.is_finite() && *a >= 0.0));
    // The state really evolved (it is not a trivially frozen field).
    assert!(
        gp.amplitude[3]
            != 1.0
                + 0.3 * (2.0 * PI * 3.0 * dx / length).cos()
                + 0.2 * (4.0 * PI * 3.0 * dx / length).sin()
    );
}

#[test]
fn uniform_condensate_phase_advances_as_minus_g_rho_t_over_hbar() {
    let (hbar, mass, g) = (1.3, 0.7, 0.4);
    let rho: f64 = 2.5;
    let (dt, steps) = (0.02, 200);
    let mut gp = solver(16, hbar, mass, g, 0.5, dt);
    gp.set_wavefunction(vec![rho.sqrt(); 16], vec![0.0; 16]);
    for _ in 0..steps {
        gp.step();
    }
    let t = dt * steps as f64;
    let expected = -g * rho * t / hbar;
    for (j, &theta) in gp.phase.iter().enumerate() {
        let diff = wrap_phase(theta - expected);
        assert!(diff.abs() < 1e-10, "site {j}: phase {theta} vs {expected}");
        assert!((gp.amplitude[j] - rho.sqrt()).abs() < 1e-12);
    }
}

#[test]
fn uniform_potential_shifts_the_phase_by_minus_v_t_over_hbar() {
    let hbar = 0.9;
    let mut gp = solver(8, hbar, 1.0, 0.0, 1.0, 0.05);
    gp.potential = vec![0.6; 8];
    for _ in 0..40 {
        gp.step();
    }
    let expected = -0.6 * 2.0 / hbar;
    for &theta in &gp.phase {
        assert!(wrap_phase(theta - expected).abs() < 1e-12);
    }
}

#[test]
fn plane_waves_follow_the_free_particle_dispersion() {
    // psi = exp(i (k x - hbar k^2 t / 2m)).  Modes with negative wavenumbers
    // sit in the upper half of the FFT output; using the wrong |k| there would
    // change the phase, so both signs and an odd (non power of two) length
    // are checked.
    let (hbar, mass) = (1.1, 0.8);
    for &(nx, mode) in &[(32usize, 5i32), (32, -3), (30, -4), (30, 7)] {
        let dx = 0.25;
        let length = nx as f64 * dx;
        let k = 2.0 * PI * f64::from(mode) / length;
        let (dt, steps) = (0.05, 100);
        let mut gp = solver(nx, hbar, mass, 0.0, dx, dt);
        let phase0: Vec<f64> = (0..nx).map(|j| wrap_phase(k * j as f64 * dx)).collect();
        gp.set_wavefunction(vec![1.0; nx], phase0);

        let velocity = gp.superfluid_velocity();
        for v in &velocity {
            assert!(
                (v - hbar * k / mass).abs() < 1e-9,
                "nx={nx} mode={mode}: wrapped phase must give v = hbar k / m, got {v}"
            );
        }

        for _ in 0..steps {
            gp.step();
        }
        let t = dt * steps as f64;
        for j in 0..nx {
            let expected = k * j as f64 * dx - hbar * k * k * t / (2.0 * mass);
            let diff = wrap_phase(gp.phase[j] - expected);
            assert!(diff.abs() < 1e-9, "nx={nx} mode={mode} site {j}: {diff}");
            assert!((gp.amplitude[j] - 1.0).abs() < 1e-11);
        }
    }
}

#[test]
fn free_gaussian_packet_moves_with_the_group_velocity() {
    let (hbar, mass) = (1.0, 2.0);
    let nx = 512;
    let dx = 0.125;
    let k0 = 1.5;
    let (x0, sigma) = (-8.0, 1.5);
    let mut gp = solver(nx, hbar, mass, 0.0, dx, 0.02);
    let x = centred_positions(nx, dx);
    let amplitude: Vec<f64> = x
        .iter()
        .map(|&xi| (-(xi - x0).powi(2) / (4.0 * sigma * sigma)).exp())
        .collect();
    let phase: Vec<f64> = x.iter().map(|&xi| wrap_phase(k0 * xi)).collect();
    gp.set_wavefunction(amplitude, phase);
    let centre = |gp: &GrossPitaevskiiLbm| -> f64 {
        let n = gp.total_number();
        x.iter()
            .zip(gp.amplitude.iter())
            .map(|(&xi, a)| xi * a * a * dx)
            .sum::<f64>()
            / n
    };
    let c0 = centre(&gp);
    assert!((c0 - x0).abs() < 1e-9);
    for _ in 0..500 {
        gp.step();
    }
    let t = 10.0;
    let expected = x0 + hbar * k0 / mass * t;
    let c1 = centre(&gp);
    assert!(
        (c1 - expected).abs() < 1e-6,
        "packet centre {c1} vs {expected}"
    );
}

#[test]
fn attractive_condensate_supports_a_stationary_bright_soliton() {
    // For hbar = m = 1 and g = -1 the profile psi = sech(x) exp(i t / 2) solves
    // the GP equation exactly.
    let nx = 256;
    let dx = 40.0 / nx as f64;
    let (dt, steps) = (0.002, 1500);
    let mut gp = solver(nx, 1.0, 1.0, -1.0, dx, dt);
    let x = centred_positions(nx, dx);
    let sech: Vec<f64> = x.iter().map(|&xi| 1.0 / xi.cosh()).collect();
    gp.set_wavefunction(sech.clone(), vec![0.0; nx]);
    let n0 = gp.total_number();
    let e0 = gp.total_energy();
    // Analytic soliton: N = 2 L_x with L_x = 1, E = -1/3 (kinetic 1/3, interaction -2/3).
    assert!((n0 - 2.0).abs() < 1e-9, "N = {n0}");
    assert!((e0 + 1.0 / 3.0).abs() < 1e-8, "E = {e0}");
    for _ in 0..steps {
        gp.step();
    }
    let t = dt * steps as f64;
    let mut worst_amplitude: f64 = 0.0;
    let mut worst_phase: f64 = 0.0;
    for j in 0..nx {
        worst_amplitude = worst_amplitude.max((gp.amplitude[j] - sech[j]).abs());
        if sech[j] > 0.1 {
            worst_phase = worst_phase.max(wrap_phase(gp.phase[j] - 0.5 * t).abs());
        }
    }
    assert!(worst_amplitude < 1e-4, "profile error {worst_amplitude}");
    assert!(worst_phase < 1e-4, "phase error {worst_phase}");
    assert!(((gp.total_number() - n0) / n0).abs() < 1e-12);
    assert!(((gp.total_energy() - e0) / e0).abs() < 1e-5);
}

#[test]
fn harmonic_oscillator_ground_state_only_acquires_a_phase() {
    // hbar = m = omega = 1: psi_0 = pi^(-1/4) exp(-x^2/2), E_0 = 1/2.
    let nx = 128;
    let dx = 24.0 / nx as f64;
    let (dt, steps) = (0.002, 1000);
    let mut gp = solver(nx, 1.0, 1.0, 0.0, dx, dt);
    let x = centred_positions(nx, dx);
    gp.potential = x.iter().map(|&xi| 0.5 * xi * xi).collect();
    let ground: Vec<f64> = x
        .iter()
        .map(|&xi| PI.powf(-0.25) * (-0.5 * xi * xi).exp())
        .collect();
    gp.set_wavefunction(ground.clone(), vec![0.0; nx]);
    assert!(
        (gp.total_energy() - 0.5).abs() < 1e-9,
        "E0 = {}",
        gp.total_energy()
    );
    for _ in 0..steps {
        gp.step();
    }
    let t = dt * steps as f64;
    for j in 0..nx {
        assert!((gp.amplitude[j] - ground[j]).abs() < 1e-5);
        if ground[j] > 0.05 {
            assert!(wrap_phase(gp.phase[j] + 0.5 * t).abs() < 1e-4);
        }
    }
}

#[test]
fn total_energy_is_conserved_by_the_splitting() {
    let nx = 64;
    let dx = 0.5;
    let length = nx as f64 * dx;
    let mut gp = solver(nx, 1.0, 1.0, 1.2, dx, 0.005);
    let x = centred_positions(nx, dx);
    gp.potential = x.iter().map(|&xi| 0.5 * 0.04 * xi * xi).collect();
    let amplitude: Vec<f64> = (0..nx)
        .map(|j| 1.0 + 0.25 * (2.0 * PI * j as f64 * dx / length).cos())
        .collect();
    let phase: Vec<f64> = (0..nx)
        .map(|j| 0.5 * (2.0 * PI * j as f64 * dx / length).sin())
        .collect();
    gp.set_wavefunction(amplitude, phase);
    let e0 = gp.total_energy();
    let mut worst: f64 = 0.0;
    for step in 0..600 {
        gp.step();
        if step % 50 == 49 {
            worst = worst.max(((gp.total_energy() - e0) / e0).abs());
        }
    }
    assert!(worst < 2e-3, "energy drift {worst}");
}

#[test]
fn plane_wave_energy_matches_the_analytic_value() {
    let (hbar, mass, g) = (1.0, 2.0, 0.6);
    let nx = 32;
    let dx = 0.5;
    let length = nx as f64 * dx;
    let k = 2.0 * PI * 3.0 / length;
    let mut gp = solver(nx, hbar, mass, g, dx, 0.01);
    let phase: Vec<f64> = (0..nx).map(|j| wrap_phase(k * j as f64 * dx)).collect();
    gp.set_wavefunction(vec![1.0; nx], phase);
    let expected = (0.5 * hbar * hbar * k * k / mass + 0.5 * g) * length;
    assert!((gp.total_energy() - expected).abs() < 1e-10);
}

#[test]
fn zero_time_step_is_an_exact_identity() {
    // Guards the normalisation of the forward/inverse transform pair.
    let nx = 20;
    let mut gp = solver(nx, 1.0, 1.0, 0.5, 0.3, 0.0);
    let amplitude: Vec<f64> = (0..nx).map(|j| 0.5 + 0.1 * j as f64).collect();
    let phase: Vec<f64> = (0..nx).map(|j| wrap_phase(0.37 * j as f64 - 1.0)).collect();
    gp.set_wavefunction(amplitude.clone(), phase.clone());
    gp.step();
    for j in 0..nx {
        assert!((gp.amplitude[j] - amplitude[j]).abs() < 1e-13);
        assert!(wrap_phase(gp.phase[j] - phase[j]).abs() < 1e-13);
    }
}

#[test]
fn wrapped_phase_differences_give_the_right_velocity() {
    let nx = 16;
    let dx = 1.0;
    let mut gp = solver(nx, 2.0, 4.0, 0.0, dx, 0.1);
    // Phase ramp of 1.2 rad per site, wrapped many times around the ring
    // (16 * 1.2 is not a multiple of 2 pi: the seam is the only defect).
    let phase: Vec<f64> = (0..nx).map(|j| wrap_phase(1.2 * j as f64)).collect();
    gp.set_wavefunction(vec![1.0; nx], phase);
    let v = gp.superfluid_velocity();
    for (j, vj) in v.iter().enumerate() {
        if j != 0 && j != nx - 1 {
            assert!((vj - 2.0 / 4.0 * 1.2).abs() < 1e-12, "site {j}: {vj}");
        }
    }
}

#[test]
fn invalid_configurations_are_reported_and_leave_the_state_untouched() {
    let base = solver(8, 1.0, 1.0, 0.1, 1.0, 0.1);

    let mut gp = base.clone();
    gp.potential.pop();
    let err = gp.try_step().expect_err("short potential");
    assert!(
        matches!(err, Error::ComponentCountMismatch { .. }),
        "{err:?}"
    );

    let mut gp = base.clone();
    gp.amplitude.push(1.0);
    assert!(matches!(
        gp.try_step().expect_err("long amplitude"),
        Error::ComponentCountMismatch { .. }
    ));

    for (hbar, mass, dx, dt) in [
        (0.0, 1.0, 1.0, 0.1),
        (1.0, -1.0, 1.0, 0.1),
        (1.0, 1.0, 0.0, 0.1),
        (1.0, 1.0, 1.0, f64::NAN),
        (f64::INFINITY, 1.0, 1.0, 0.1),
    ] {
        let mut gp = solver(8, hbar, mass, 0.1, dx, dt);
        gp.phase[3] = 0.25;
        let before = (gp.amplitude.clone(), gp.phase.clone());
        let err = gp.try_step().expect_err("invalid parameter");
        assert!(matches!(err, Error::ParameterOutOfBounds { .. }), "{err:?}");
        gp.step();
        assert_eq!((gp.amplitude.clone(), gp.phase.clone()), before);
    }

    let mut empty = solver(0, 1.0, 1.0, 0.1, 1.0, 0.1);
    assert!(matches!(
        empty.try_step().expect_err("no sites"),
        Error::InvalidGridDimensions { .. }
    ));
    empty.step();
    assert_eq!(empty.total_energy(), 0.0);
    assert!(empty.superfluid_velocity().is_empty());
}
