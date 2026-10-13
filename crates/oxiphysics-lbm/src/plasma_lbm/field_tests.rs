//! Tests of the electrostatic coupling, the Boris push and the continuity update
//! of the two-fluid plasma solver, including a Langmuir-oscillation check.

use super::*;

fn rel(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1e-300)
}

fn neutral_sim(nx: usize, ny: usize, n0: f64, bz: f64) -> PlasmaLbm {
    let ion = PlasmaSpecies::new(1.0, 1836.0, 1.0e3, [0.0; 3]);
    let elec = PlasmaSpecies::new(-1.0, 1.0, 1.0e4, [0.0; 3]);
    PlasmaLbm::new(nx, ny, n0, n0, ion, elec, bz)
}

/// `rho0 cos(kx x + ky y)` on the grid with the integer mode numbers `(mx, my)`;
/// returns the field and the wavenumbers.
fn cosine_charge(
    nx: usize,
    ny: usize,
    dx: f64,
    modes: (usize, usize),
    rho0: f64,
) -> (Vec<f64>, f64, f64) {
    let kx = 2.0 * PI * modes.0 as f64 / (nx as f64 * dx);
    let ky = 2.0 * PI * modes.1 as f64 / (ny as f64 * dx);
    let mut rho = vec![0.0; nx * ny];
    for iy in 0..ny {
        for ix in 0..nx {
            rho[iy * nx + ix] = rho0 * (kx * ix as f64 * dx + ky * iy as f64 * dx).cos();
        }
    }
    (rho, kx, ky)
}

// ---------------------------------------------------------------------------
// Poisson solve
// ---------------------------------------------------------------------------

#[test]
fn poisson_field_of_a_sinusoidal_charge_is_analytic() {
    let (nx, ny, dx) = (16, 8, 2.0e-6);
    let rho0 = 1.0e-3; // C m^-3
    for &modes in &[(1, 0), (2, 0), (0, 1), (2, 1), (3, 2)] {
        let (rho, kx, ky) = cosine_charge(nx, ny, dx, modes, rho0);
        let (ex, ey) = poisson_electric_field_periodic(&rho, nx, ny, dx).expect("valid input");
        let k2 = kx * kx + ky * ky;
        let e0 = rho0 / (EPSILON_0 * k2.sqrt());
        for iy in 0..ny {
            for ix in 0..nx {
                let theta = kx * ix as f64 * dx + ky * iy as f64 * dx;
                // phi = rho0 cos(theta) / (eps0 k^2),  E = -grad phi
                let amplitude = rho0 * theta.sin() / (EPSILON_0 * k2);
                let idx = iy * nx + ix;
                assert!(
                    (ex[idx] - amplitude * kx).abs() < 1e-10 * e0,
                    "modes {modes:?} Ex[{ix},{iy}] = {} expected {}",
                    ex[idx],
                    amplitude * kx
                );
                assert!(
                    (ey[idx] - amplitude * ky).abs() < 1e-10 * e0,
                    "modes {modes:?} Ey[{ix},{iy}] = {} expected {}",
                    ey[idx],
                    amplitude * ky
                );
            }
        }
    }
}

#[test]
fn poisson_field_satisfies_gauss_law_to_second_order() {
    // div E = rho / eps0: the central-difference divergence of the field of a
    // well-resolved mode (64 cells per wavelength) reproduces the charge density
    // up to the O((k dx)^2) truncation error of the difference.
    let (nx, ny, dx) = (64, 2, 1.0e-6);
    let rho0 = 5.0e-4;
    let (rho, _, _) = cosine_charge(nx, ny, dx, (1, 0), rho0);
    let (ex, _) = poisson_electric_field_periodic(&rho, nx, ny, dx).expect("valid input");
    for ix in 0..nx {
        let divergence = (ex[(ix + 1) % nx] - ex[(ix + nx - 1) % nx]) / (2.0 * dx);
        let expected = rho[ix] / EPSILON_0;
        assert!(
            (divergence - expected).abs() < 5e-3 * rho0 / EPSILON_0,
            "div E at {ix} = {divergence}, rho/eps0 = {expected}"
        );
    }
}

#[test]
fn poisson_field_of_a_uniform_charge_vanishes() {
    let (nx, ny, dx) = (8, 6, 1.0e-6);
    let rho = vec![0.16; nx * ny];
    let (ex, ey) = poisson_electric_field_periodic(&rho, nx, ny, dx).expect("valid input");
    assert!(ex.iter().all(|e| e.abs() < 1e-6), "Ex = {ex:?}");
    assert!(ey.iter().all(|e| e.abs() < 1e-6), "Ey = {ey:?}");
}

#[test]
fn poisson_nyquist_mode_produces_no_field() {
    // (-1)^ix is the Nyquist mode of an even-length axis; its derivative is not
    // representable by a real field and must not leak into E.
    let (nx, ny, dx) = (8, 4, 1.0e-6);
    let mut rho = vec![0.0; nx * ny];
    for iy in 0..ny {
        for ix in 0..nx {
            rho[iy * nx + ix] = if ix % 2 == 0 { 1.0e-3 } else { -1.0e-3 };
        }
    }
    let (ex, ey) = poisson_electric_field_periodic(&rho, nx, ny, dx).expect("valid input");
    assert!(ex.iter().all(|e| e.abs() < 1e-6), "Ex = {ex:?}");
    assert!(ey.iter().all(|e| e.abs() < 1e-6), "Ey = {ey:?}");
}

#[test]
fn poisson_rejects_invalid_input() {
    let rho = vec![0.0; 12];
    assert!(poisson_electric_field_periodic(&rho, 0, 4, 1e-6).is_err());
    assert!(poisson_electric_field_periodic(&rho, 4, 4, 1e-6).is_err());
    assert!(poisson_electric_field_periodic(&rho, 4, 3, 0.0).is_err());
    assert!(poisson_electric_field_periodic(&rho, 4, 3, f64::NAN).is_err());
    let mut bad = rho.clone();
    bad[5] = f64::INFINITY;
    assert!(poisson_electric_field_periodic(&bad, 4, 3, 1e-6).is_err());
    assert!(poisson_electric_field_periodic(&rho, 4, 3, 1e-6).is_ok());
}

#[test]
fn single_cell_grid_has_no_field() {
    let (ex, ey) = poisson_electric_field_periodic(&[0.5], 1, 1, 1e-6).expect("valid input");
    assert_eq!((ex[0], ey[0]), (0.0, 0.0));
}

#[test]
fn step_stores_the_electrostatic_field_of_the_charge_separation() {
    let (nx, ny) = (16, 4);
    let n0 = 1.0e18;
    let alpha = 1.0e-2;
    let mut sim = neutral_sim(nx, ny, n0, 0.0).with_grid_spacing(2.0e-6);
    let kx = 2.0 * PI / (nx as f64 * sim.dx);
    for iy in 0..ny {
        for ix in 0..nx {
            sim.electrons.density[iy * nx + ix] =
                n0 * (1.0 + alpha * (kx * ix as f64 * sim.dx).cos());
        }
    }
    sim.try_step(1.0e-16).expect("valid step");
    // rho_c = e (n_i - n_e) = -e alpha n0 cos(kx x);  E_x = rho0 sin(kx x) / (eps0 kx)
    let rho0 = -ELEM_CHARGE * alpha * n0;
    let e0 = rho0.abs() / (EPSILON_0 * kx);
    for iy in 0..ny {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            let expected = rho0 * (kx * ix as f64 * sim.dx).sin() / (EPSILON_0 * kx);
            assert!(
                (sim.electric_field.ex[idx] - expected).abs() < 1e-9 * e0,
                "Ex[{ix},{iy}] = {} expected {expected}",
                sim.electric_field.ex[idx]
            );
            assert!(
                sim.electric_field.ey[idx].abs() < 1e-9 * e0,
                "Ey must vanish"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Boris push
// ---------------------------------------------------------------------------

#[test]
fn boris_rotates_without_changing_the_speed() {
    let (q, m, bz) = (-ELEM_CHARGE, ELECTRON_MASS, 1.0);
    let omega_c = ELEM_CHARGE * bz / ELECTRON_MASS;
    let dt = 1.0e-12; // omega_c dt = 0.176
    let v0 = [1.0e6, 2.0e5];

    // Rotation angle of the first steps: 2 atan(omega_c dt / 2) per step,
    // counter-clockwise for an electron in B_z > 0.
    let mut v = v0;
    let steps = 10;
    for _ in 0..steps {
        v = boris_push_2d(v, [0.0, 0.0], bz, q, m, dt);
    }
    let angle = v[1].atan2(v[0]) - v0[1].atan2(v0[0]);
    let expected = steps as f64 * 2.0 * (0.5 * omega_c * dt).atan();
    assert!(
        (angle - expected).abs() < 1e-9,
        "rotation {angle} vs {expected}"
    );

    // Long runs conserve |v| to rounding error.
    let speed0 = (v0[0] * v0[0] + v0[1] * v0[1]).sqrt();
    let mut v = v0;
    for _ in 0..100_000 {
        v = boris_push_2d(v, [0.0, 0.0], bz, q, m, dt);
    }
    let speed = (v[0] * v[0] + v[1] * v[1]).sqrt();
    assert!(
        (speed - speed0).abs() < 1e-10 * speed0,
        "|v| drifted from {speed0} to {speed}"
    );
}

#[test]
fn boris_keeps_the_exb_drift_for_both_charge_signs() {
    let (bz, e0) = (1.0, 1.0e3);
    let drift = [e0 / bz, 0.0]; // E x B / B^2 for E = (0, E0), B = (0, 0, Bz)
    for &(q, m) in &[
        (ELEM_CHARGE, 1836.0 * ELECTRON_MASS),
        (-ELEM_CHARGE, ELECTRON_MASS),
    ] {
        let mut v = drift;
        for _ in 0..1000 {
            v = boris_push_2d(v, [0.0, e0], bz, q, m, 1.0e-12);
        }
        assert!(
            (v[0] - drift[0]).abs() < 1e-9 * drift[0] && v[1].abs() < 1e-9 * drift[0],
            "drift lost for q = {q}: v = {v:?}"
        );
    }
}

#[test]
fn boris_without_magnetic_field_is_a_uniform_acceleration() {
    let (q, m, dt) = (ELEM_CHARGE, 1836.0 * ELECTRON_MASS, 1.0e-9);
    let e = [3.0e2, -1.0e2];
    let v = boris_push_2d([10.0, 20.0], e, 0.0, q, m, dt);
    assert!(rel(v[0], 10.0 + q / m * e[0] * dt) < 1e-14);
    assert!(rel(v[1], 20.0 + q / m * e[1] * dt) < 1e-14);
}

#[test]
fn boris_rotation_sense_matches_the_magnetic_force() {
    // For a positive charge and B_z > 0 the force on v = (1, 0) points to -y.
    let bf = MagneticField::new(1.0);
    let force = bf.magnetic_force(ELEM_CHARGE, [1.0, 0.0, 0.0]);
    let v = boris_push_2d(
        [1.0, 0.0],
        [0.0, 0.0],
        1.0,
        ELEM_CHARGE,
        1836.0 * ELECTRON_MASS,
        1.0e-9,
    );
    assert!(force[1] < 0.0 && v[1] < 0.0, "force {force:?}, v = {v:?}");
}

#[test]
fn plasma_step_conserves_the_speed_of_gyrating_species() {
    let mut sim = neutral_sim(1, 1, 1.0e18, 1.0);
    sim.electrons.vx[0] = 1.0e6;
    sim.electrons.vy[0] = -3.0e5;
    sim.ions.vx[0] = 2.0e4;
    sim.ions.vy[0] = 1.0e4;
    let speed = |vx: f64, vy: f64| (vx * vx + vy * vy).sqrt();
    let s_e = speed(sim.electrons.vx[0], sim.electrons.vy[0]);
    let s_i = speed(sim.ions.vx[0], sim.ions.vy[0]);
    for _ in 0..5_000 {
        sim.try_step(1.0e-12).expect("valid step");
    }
    let s_e1 = speed(sim.electrons.vx[0], sim.electrons.vy[0]);
    let s_i1 = speed(sim.ions.vx[0], sim.ions.vy[0]);
    assert!(
        (s_e1 - s_e).abs() < 1e-10 * s_e,
        "electron speed {s_e} -> {s_e1}"
    );
    assert!(
        (s_i1 - s_i).abs() < 1e-10 * s_i,
        "ion speed {s_i} -> {s_i1}"
    );
}

// ---------------------------------------------------------------------------
// Continuity
// ---------------------------------------------------------------------------

#[test]
fn upwind_continuity_conserves_particles_and_transports_a_bump() {
    let (nx, ny) = (64, 4);
    let dx = 1.0e-6;
    let speed = 2.0e3; // m/s
    let dt = 1.0e-9; // 2 cells per step -> CFL sub-cycling
    let mut density = vec![0.0; nx * ny];
    for iy in 0..ny {
        for ix in 0..nx {
            let d = ix as f64 - 20.0;
            density[iy * nx + ix] = 1.0e18 * (1.0 + (-d * d / 8.0).exp());
        }
    }
    let vx = vec![speed; nx * ny];
    let vy = vec![0.0; nx * ny];
    let moment = |n: &[f64]| -> (f64, f64) {
        let mut total = 0.0;
        let mut first = 0.0;
        for iy in 0..ny {
            for ix in 0..nx {
                let excess = n[iy * nx + ix] - 1.0e18;
                total += excess;
                first += ix as f64 * excess;
            }
        }
        (total, first / total)
    };
    let (m0, x0) = moment(&density);
    let sum0: f64 = density.iter().sum();
    advance_density_upwind(&mut density, &vx, &vy, nx, ny, dx, dt).expect("valid step");
    let (m1, x1) = moment(&density);
    let sum1: f64 = density.iter().sum();
    assert!(
        (sum1 - sum0).abs() < 1e-12 * sum0,
        "particle number {sum0} -> {sum1}"
    );
    assert!((m1 - m0).abs() < 1e-9 * m0.abs());
    // The excess density moves by v dt / dx = 2 cells.
    assert!((x1 - x0 - 2.0).abs() < 1e-8, "shift {}", x1 - x0);
    assert!(
        density.iter().all(|&n| n >= 0.0),
        "density must stay non-negative"
    );
}

#[test]
fn upwind_continuity_transports_in_both_directions() {
    let (nx, ny) = (12, 12);
    let dx = 1.0e-6;
    let mut density = vec![1.0e18; nx * ny];
    density[5 * nx + 5] = 3.0e18;
    let vx = vec![-1.0e3; nx * ny];
    let vy = vec![5.0e2; nx * ny];
    let sum0: f64 = density.iter().sum();
    for _ in 0..4 {
        advance_density_upwind(&mut density, &vx, &vy, nx, ny, dx, 1.0e-9).expect("valid step");
    }
    let sum1: f64 = density.iter().sum();
    assert!((sum1 - sum0).abs() < 1e-12 * sum0);
    // Mass flows towards -x and +y: the maximum has left its start cell.
    let (mut best, mut best_idx) = (0.0, 0);
    for (k, &n) in density.iter().enumerate() {
        if n > best {
            best = n;
            best_idx = k;
        }
    }
    let (bx, by) = (best_idx % nx, best_idx / nx);
    assert!(bx < 5 && by > 5, "peak at ({bx}, {by})");
}

#[test]
fn step_reports_a_time_step_beyond_the_substep_limit_and_changes_nothing() {
    let mut sim = neutral_sim(8, 8, 1.0e18, 0.0).with_grid_spacing(1.0e-6);
    for v in sim.electrons.vx.iter_mut() {
        *v = 1.0e9;
    }
    let before = sim.clone();
    let result = sim.try_step(1.0e-3);
    assert!(
        matches!(result, Err(Error::TimeStepTooLarge { .. })),
        "{result:?}"
    );
    assert_eq!(sim.electrons.vx, before.electrons.vx);
    assert_eq!(sim.electrons.density, before.electrons.density);
    assert_eq!(sim.ions.density, before.ions.density);
}

// ---------------------------------------------------------------------------
// Whole-solver behaviour
// ---------------------------------------------------------------------------

#[test]
fn quasi_neutral_plasma_stays_at_rest() {
    let mut sim = neutral_sim(8, 8, 1.0e18, 0.0);
    for _ in 0..20 {
        sim.try_step(1.0e-12).expect("valid step");
    }
    assert!(
        sim.electrons.vx.iter().all(|&v| v.abs() < 1e-9),
        "electron vx"
    );
    assert!(sim.ions.vy.iter().all(|&v| v.abs() < 1e-9), "ion vy");
    assert!(
        sim.electrons
            .density
            .iter()
            .all(|&n| rel(n, 1.0e18) < 1e-12)
    );
}

#[test]
fn langmuir_oscillation_has_the_plasma_frequency() {
    let (nx, ny) = (32, 4);
    let n0 = 1.0e18;
    let dx = 1.0e-6;
    let alpha = 1.0e-3;
    let mut sim = neutral_sim(nx, ny, n0, 0.0).with_grid_spacing(dx);
    let k = 2.0 * PI / (nx as f64 * dx);
    for iy in 0..ny {
        for ix in 0..nx {
            sim.electrons.density[iy * nx + ix] = n0 * (1.0 + alpha * (k * ix as f64 * dx).cos());
        }
    }
    let wp = plasma_frequency(n0, ELEM_CHARGE, ELECTRON_MASS);
    let dt = 0.02 / wp;
    let amplitude = |s: &PlasmaLbm| -> f64 {
        let mut sum = 0.0;
        for iy in 0..ny {
            for ix in 0..nx {
                sum += s.electrons.density[iy * nx + ix] * (k * ix as f64 * dx).cos();
            }
        }
        2.0 * sum / (nx * ny) as f64
    };
    let mut series = vec![amplitude(&sim)];
    for _ in 0..300 {
        sim.try_step(dt).expect("valid step");
        series.push(amplitude(&sim));
    }
    assert!(rel(series[0], alpha * n0) < 1e-12);

    // Zero crossings of the cosine-like oscillation are half a period apart.
    let mut crossings = Vec::new();
    for j in 1..series.len() {
        if series[j - 1] * series[j] < 0.0 {
            let frac = series[j - 1] / (series[j - 1] - series[j]);
            crossings.push((j as f64 - 1.0 + frac) * dt);
        }
    }
    assert!(
        crossings.len() >= 2,
        "found {} zero crossings",
        crossings.len()
    );
    let omega = PI / (crossings[1] - crossings[0]);
    assert!(
        rel(omega, wp) < 0.02,
        "oscillation frequency {omega} rad/s vs plasma frequency {wp} rad/s"
    );
    // The amplitude neither decays nor grows.
    let lowest = series.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(
        (-lowest - alpha * n0).abs() < 0.03 * alpha * n0,
        "swing to {lowest}, initial amplitude {}",
        alpha * n0
    );
    // The field stays one-dimensional.
    let ex_max = sim
        .electric_field
        .ex
        .iter()
        .fold(0.0_f64, |m, &e| m.max(e.abs()));
    let ey_max = sim
        .electric_field
        .ey
        .iter()
        .fold(0.0_f64, |m, &e| m.max(e.abs()));
    assert!(ey_max < 1e-8 * ex_max, "Ey = {ey_max}, Ex = {ex_max}");
}

#[test]
fn invalid_steps_are_rejected_and_leave_the_state_untouched() {
    let mut sim = neutral_sim(4, 4, 1.0e18, 0.1);
    sim.electrons.vx[3] = 123.0;
    let before = sim.clone();
    for &dt in &[0.0, -1.0e-12, f64::NAN, f64::INFINITY] {
        assert!(sim.try_step(dt).is_err(), "dt = {dt}");
        sim.step(dt);
        assert_eq!(sim.electrons.vx, before.electrons.vx);
        assert_eq!(sim.ions.density, before.ions.density);
    }
    let mut bad_dx = before.clone();
    bad_dx.dx = 0.0;
    assert!(bad_dx.try_step(1.0e-12).is_err());
    let mut bad_mass = before.clone();
    bad_mass.electrons.species.mass_ratio = 0.0;
    assert!(bad_mass.try_step(1.0e-12).is_err());
    let mut bad_len = before.clone();
    bad_len.ions.density.pop();
    assert!(bad_len.try_step(1.0e-12).is_err());
    let mut nan_field = before.clone();
    nan_field.ions.vy[1] = f64::NAN;
    assert!(nan_field.try_step(1.0e-12).is_err());
}

#[test]
fn default_grid_spacing_is_the_electron_debye_length() {
    let sim = neutral_sim(4, 4, 1.0e18, 0.0);
    let debye = debye_length(1.0e18, 1.0e4, ELEM_CHARGE);
    assert!(rel(sim.dx, debye) < 1e-12);
    let sim = sim.with_grid_spacing(3.5e-6);
    assert_eq!(sim.dx, 3.5e-6);
    // No electrons -> no Debye length -> fall back to 1 micrometre.
    let ion = PlasmaSpecies::new(1.0, 1836.0, 1.0e3, [0.0; 3]);
    let elec = PlasmaSpecies::new(-1.0, 1.0, 1.0e4, [0.0; 3]);
    let empty = PlasmaLbm::new(2, 2, 1.0e18, 0.0, ion, elec, 0.0);
    assert_eq!(empty.dx, 1.0e-6);
}

#[test]
fn charge_density_uses_the_species_charges() {
    let ion = PlasmaSpecies::new(2.0, 7344.0, 1.0e3, [0.0; 3]);
    let elec = PlasmaSpecies::new(-1.0, 1.0, 1.0e4, [0.0; 3]);
    let sim = PlasmaLbm::new(2, 2, 1.0e18, 2.0e18, ion, elec, 0.0);
    // Doubly charged ions at n_i and singly charged electrons at 2 n_i: neutral.
    for &rho in &sim.charge_density() {
        assert!(rho.abs() < 1e-12 * ELEM_CHARGE * 1.0e18, "rho = {rho}");
    }
}
