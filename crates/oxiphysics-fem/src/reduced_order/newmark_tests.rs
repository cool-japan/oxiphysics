// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the reduced-order Newmark integrator
//! ([`RomSystem::step_newmark`]).
//!
//! The average-acceleration scheme (`β = ¼, γ = ½`) is the trapezoidal rule of
//! the first-order system `(q, q̇)`, so for an undamped single mode it advances
//! the state by the exact rotation `θ = 2·atan(ω·dt/2)` per step and conserves
//! the energy; these closed forms are checked to rounding accuracy.

use super::*;

const BETA: f64 = 0.25;
const GAMMA: f64 = 0.5;

/// One-mode ROM (`Φ = [1]`) with `M_r = m`, `K_r = k`, `C_r = c`.
fn single_mode(m: f64, k: f64, c: f64) -> RomSystem {
    let basis = PodBasis {
        modes: vec![vec![1.0]],
        singular_values: vec![1.0],
        n_modes: 1,
    };
    let mut rom = RomSystem::project_system(&basis, &[vec![m]], &[vec![k]]);
    rom.reduced_damping = vec![vec![c]];
    rom
}

/// Exact discrete rotation angle of the average-acceleration scheme.
fn discrete_angle(omega: f64, dt: f64) -> f64 {
    2.0 * (0.5 * omega * dt).atan()
}

#[test]
fn undamped_mode_follows_the_exact_discrete_cosine_and_conserves_energy() {
    let omega = 2.0 * PI;
    let (m, k) = (1.0, omega * omega);
    let dt = 1.0 / 200.0;
    let mut rom = single_mode(m, k, 0.0);
    rom.set_initial_state(&[1.0], &[0.0], &[0.0])
        .expect("initial state");
    assert!((rom.reduced_acceleration()[0] + omega * omega).abs() < 1e-9);
    let theta = discrete_angle(omega, dt);
    let initial_energy = 0.5 * k;
    for step in 1..=100 {
        let u = rom
            .step_newmark(&[0.0], dt, BETA, GAMMA)
            .expect("Newmark step");
        let q = u[0];
        let t = step as f64 * dt;
        // Closed form of the scheme ...
        assert!(
            (q - (step as f64 * theta).cos()).abs() < 1e-12,
            "step {step}: {q}"
        );
        // ... which follows q0 cos(omega t) to the expected period-elongation error.
        assert!(
            (q - (omega * t).cos()).abs() < 1e-3,
            "step {step}: {q} vs {}",
            (omega * t).cos()
        );
        // The acceleration is stored and satisfies the equation of motion.
        let a = rom.reduced_acceleration()[0];
        assert!((a + omega * omega * q).abs() < 1e-9, "step {step}: a = {a}");
        let v = rom.reduced_velocity()[0];
        let energy = 0.5 * m * v * v + 0.5 * k * q * q;
        assert!(
            (energy - initial_energy).abs() < 1e-10 * initial_energy,
            "step {step}: energy {energy}"
        );
    }
}

#[test]
fn damped_mode_matches_the_analytic_decay() {
    let omega = 2.0 * PI;
    let zeta = 0.05;
    let (m, k) = (2.0, 2.0 * omega * omega);
    let c = 2.0 * zeta * omega * m;
    let omega_d = omega * (1.0 - zeta * zeta).sqrt();
    let exact = |t: f64| {
        (-zeta * omega * t).exp()
            * ((omega_d * t).cos() + zeta * omega / omega_d * (omega_d * t).sin())
    };
    let run = |n_steps: usize| -> f64 {
        let dt = 1.0 / n_steps as f64;
        let mut rom = single_mode(m, k, c);
        rom.set_initial_state(&[1.0], &[0.0], &[0.0])
            .expect("initial state");
        let mut max_error = 0.0_f64;
        for step in 1..=n_steps {
            let u = rom
                .step_newmark(&[0.0], dt, BETA, GAMMA)
                .expect("Newmark step");
            max_error = max_error.max((u[0] - exact(step as f64 * dt)).abs());
        }
        max_error
    };
    let coarse = run(200);
    let fine = run(400);
    assert!(coarse < 2e-3, "coarse error {coarse}");
    assert!(fine < 6e-4, "fine error {fine}");
    // Second-order accuracy: halving dt divides the error by about four.
    let ratio = coarse / fine;
    assert!((3.5..4.5).contains(&ratio), "convergence ratio {ratio}");
}

#[test]
fn step_load_oscillates_about_the_static_deflection() {
    let omega = 3.0;
    let (m, k, force) = (2.0, 2.0 * 9.0, 5.0);
    let dt = 0.01;
    let theta = discrete_angle(omega, dt);
    let mut rom = single_mode(m, k, 0.0);
    // f(0) = F is applied from the start: the consistent acceleration is F / m.
    rom.set_initial_state(&[0.0], &[0.0], &[force])
        .expect("initial state");
    assert!((rom.reduced_acceleration()[0] - force / m).abs() < 1e-12);
    let static_deflection = force / k;
    for step in 1..=300 {
        let u = rom
            .step_newmark(&[force], dt, BETA, GAMMA)
            .expect("Newmark step");
        let expected = static_deflection * (1.0 - (step as f64 * theta).cos());
        assert!(
            (u[0] - expected).abs() < 1e-12,
            "step {step}: {} vs {expected}",
            u[0]
        );
    }
}

#[test]
fn constant_force_at_rest_is_in_equilibrium_at_the_static_deflection() {
    // Starting at q = F/k with a = 0 the scheme must stay put (the old effective
    // force dropped the M q~ and C terms and drifted away).
    let (m, k, c, force) = (1.5, 40.0, 0.7, 10.0);
    let mut rom = single_mode(m, k, c);
    rom.set_initial_state(&[force / k], &[0.0], &[force])
        .expect("initial state");
    assert!(rom.reduced_acceleration()[0].abs() < 1e-12);
    for _ in 0..50 {
        let u = rom
            .step_newmark(&[force], 0.02, BETA, GAMMA)
            .expect("Newmark step");
        assert!((u[0] - force / k).abs() < 1e-12);
    }
}

/// Orthonormal basis of three mixed sine modes of an 8-DOF chain.
fn chain_basis() -> PodBasis {
    let n = 8;
    let sine = |j: usize| -> Vec<f64> {
        (1..=n)
            .map(|i| (j as f64 * PI * i as f64 / (n + 1) as f64).sin())
            .collect()
    };
    let (s1, s2, s3) = (sine(1), sine(2), sine(3));
    let mixed = [
        (0..n).map(|i| s1[i] + 0.3 * s2[i]).collect::<Vec<f64>>(),
        (0..n).map(|i| s2[i] - 0.2 * s3[i]).collect::<Vec<f64>>(),
        s3.clone(),
    ];
    let mut modes: Vec<Vec<f64>> = Vec::new();
    for v in mixed {
        let mut w = v;
        for phi in &modes {
            let proj = dot(&w, phi);
            for (wi, pi) in w.iter_mut().zip(phi.iter()) {
                *wi -= proj * pi;
            }
        }
        let nrm = norm(&w);
        modes.push(vec_scale(&w, 1.0 / nrm));
    }
    PodBasis {
        modes,
        singular_values: vec![1.0; 3],
        n_modes: 3,
    }
}

fn chain_matrices(n: usize, stiffness: f64) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut mass = vec![vec![0.0; n]; n];
    let mut k = vec![vec![0.0; n]; n];
    for i in 0..n {
        mass[i][i] = 1.0 + 0.1 * i as f64;
        k[i][i] = 2.0 * stiffness;
        if i + 1 < n {
            k[i][i + 1] = -stiffness;
            k[i + 1][i] = -stiffness;
        }
    }
    (mass, k)
}

#[test]
fn coupled_modes_conserve_energy_and_agree_with_runge_kutta() {
    let basis = chain_basis();
    let (mass, stiffness) = chain_matrices(8, 100.0);
    let mut rom = RomSystem::project_system(&basis, &mass, &stiffness);
    // Initial displacement inside the span of the basis, zero velocity.
    let u0 = basis.reconstruct(&[0.01, -0.02, 0.005]);
    let zero = vec![0.0; 8];
    rom.set_initial_state(&u0, &zero, &zero)
        .expect("initial state");
    let q0 = rom.reduced_displacement().to_vec();

    // Reference: classical RK4 on the reduced first-order system.
    let m_r = rom.reduced_mass.clone();
    let k_r = rom.reduced_stiffness.clone();
    let m_factor = cholesky_spd(&m_r).expect("M_r is positive definite");
    let accel = |q: &[f64]| -> Vec<f64> {
        let kq = mat_vec(&k_r, q);
        cholesky_solve(&m_factor, &kq.iter().map(|x| -x).collect::<Vec<f64>>())
    };
    let rk4_step = |q: &[f64], v: &[f64], h: f64| -> (Vec<f64>, Vec<f64>) {
        let add = |a: &[f64], b: &[f64], s: f64| -> Vec<f64> {
            a.iter().zip(b.iter()).map(|(x, y)| x + s * y).collect()
        };
        let (k1q, k1v) = (v.to_vec(), accel(q));
        let (k2q, k2v) = (add(v, &k1v, 0.5 * h), accel(&add(q, &k1q, 0.5 * h)));
        let (k3q, k3v) = (add(v, &k2v, 0.5 * h), accel(&add(q, &k2q, 0.5 * h)));
        let (k4q, k4v) = (add(v, &k3v, h), accel(&add(q, &k3q, h)));
        let q_new: Vec<f64> = (0..q.len())
            .map(|i| q[i] + h / 6.0 * (k1q[i] + 2.0 * k2q[i] + 2.0 * k3q[i] + k4q[i]))
            .collect();
        let v_new: Vec<f64> = (0..q.len())
            .map(|i| v[i] + h / 6.0 * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]))
            .collect();
        (q_new, v_new)
    };

    let energy = |q: &[f64], v: &[f64]| -> f64 {
        0.5 * dot(v, &mat_vec(&m_r, v)) + 0.5 * dot(q, &mat_vec(&k_r, q))
    };
    let initial_energy = energy(&q0, &[0.0; 3]);
    assert!(initial_energy > 0.0);

    let dt = 2.0e-3;
    let sub_steps = 20;
    let (mut q_ref, mut v_ref) = (q0, vec![0.0; 3]);
    for step in 1..=250 {
        let u = rom
            .step_newmark(&zero, dt, BETA, GAMMA)
            .expect("Newmark step");
        assert_eq!(u.len(), 8);
        for _ in 0..sub_steps {
            let (q, v) = rk4_step(&q_ref, &v_ref, dt / sub_steps as f64);
            q_ref = q;
            v_ref = v;
        }
        let e = energy(rom.reduced_displacement(), rom.reduced_velocity());
        assert!(
            (e - initial_energy).abs() < 1e-10 * initial_energy,
            "step {step}: energy {e} vs {initial_energy}"
        );
        let scale = q_ref.iter().map(|x| x.abs()).fold(0.0, f64::max).max(1e-3);
        for (a, b) in rom.reduced_displacement().iter().zip(q_ref.iter()) {
            assert!(
                (a - b).abs() < 2e-3 * scale,
                "step {step}: Newmark {a} vs RK4 {b}"
            );
        }
    }
}

#[test]
fn damped_coupled_modes_dissipate_energy_monotonically() {
    let basis = chain_basis();
    let (mass, stiffness) = chain_matrices(8, 100.0);
    let mut rom = RomSystem::project_system(&basis, &mass, &stiffness);
    rom.set_rayleigh_damping(2.0, 0.002);
    let u0 = basis.reconstruct(&[0.01, 0.0, 0.0]);
    let zero = vec![0.0; 8];
    rom.set_initial_state(&u0, &zero, &zero)
        .expect("initial state");
    let energy = |rom: &RomSystem| -> f64 {
        let (q, v) = (rom.reduced_displacement(), rom.reduced_velocity());
        0.5 * dot(v, &mat_vec(&rom.reduced_mass, v))
            + 0.5 * dot(q, &mat_vec(&rom.reduced_stiffness, q))
    };
    let mut previous = energy(&rom);
    let first = previous;
    for step in 1..=400 {
        rom.step_newmark(&zero, 5.0e-3, BETA, GAMMA)
            .expect("Newmark step");
        let e = energy(&rom);
        assert!(
            e <= previous * (1.0 + 1e-12),
            "step {step}: {e} > {previous}"
        );
        previous = e;
    }
    assert!(previous < 0.2 * first, "energy fell only to {previous}");
}

#[test]
fn cached_factor_is_reused_and_invalidated_by_matrix_edits() {
    let mut rom = single_mode(1.0, 100.0, 0.0);
    rom.set_initial_state(&[0.3], &[0.0], &[0.0])
        .expect("initial state");
    assert!(rom.newmark_factor.is_none());
    rom.step_newmark(&[0.0], 0.01, BETA, GAMMA).expect("step");
    let first = rom.newmark_factor.clone().expect("factor cached");
    rom.step_newmark(&[0.0], 0.01, BETA, GAMMA).expect("step");
    let second = rom.newmark_factor.clone().expect("factor cached");
    assert_eq!(first.k_eff, second.k_eff);
    assert_eq!(first.chol, second.chol);

    // A different time step needs a different effective stiffness ...
    let mut other_dt = rom.clone();
    other_dt
        .step_newmark(&[0.0], 0.02, BETA, GAMMA)
        .expect("step");
    let refactored = other_dt.newmark_factor.expect("factor cached");
    assert!(refactored.k_eff != second.k_eff);

    // ... and so does an edit of the (public) reduced stiffness: the step must use
    // the new value, checked against the scalar closed form.
    let mut edited = rom.clone();
    edited.reduced_stiffness[0][0] = 250.0;
    let (q, v, a) = (edited.disp[0], edited.vel[0], edited.acc[0]);
    let dt = 0.01;
    let q_pred = q + dt * v + dt * dt * (0.5 - BETA) * a;
    let c1 = 1.0 / (BETA * dt * dt);
    let expected = c1 * q_pred / (c1 + 250.0);
    let u = edited.step_newmark(&[0.0], dt, BETA, GAMMA).expect("step");
    assert!((u[0] - expected).abs() < 1e-14, "{} vs {expected}", u[0]);
}

#[test]
fn invalid_arguments_are_rejected_and_leave_the_state_untouched() {
    let mut rom = single_mode(1.0, 100.0, 0.0);
    rom.set_initial_state(&[0.3], &[0.1], &[0.0])
        .expect("initial state");
    let (q, v, a) = (rom.disp.clone(), rom.vel.clone(), rom.acc.clone());
    let invalid = |result: Result<Vec<f64>, FemError>| {
        assert!(
            matches!(result, Err(FemError::InvalidInput(_))),
            "unexpected result {result:?}"
        );
    };
    invalid(rom.step_newmark(&[0.0], 0.0, BETA, GAMMA));
    invalid(rom.step_newmark(&[0.0], -0.01, BETA, GAMMA));
    invalid(rom.step_newmark(&[0.0], f64::NAN, BETA, GAMMA));
    invalid(rom.step_newmark(&[0.0], 0.01, 0.0, GAMMA));
    invalid(rom.step_newmark(&[0.0], 0.01, BETA, -0.5));
    invalid(rom.step_newmark(&[0.0], 0.01, BETA, f64::INFINITY));
    invalid(rom.step_newmark(&[0.0, 1.0], 0.01, BETA, GAMMA));
    assert_eq!(
        (rom.disp.clone(), rom.vel.clone(), rom.acc.clone()),
        (q, v, a)
    );

    // Negative stiffness beyond M/(beta dt^2): the effective stiffness is indefinite.
    let mut unstable = single_mode(1.0, -1.0e6, 0.0);
    let result = unstable.step_newmark(&[0.0], 0.1, BETA, GAMMA);
    assert_eq!(result, Err(FemError::SingularSystem));
    assert_eq!(unstable.disp, vec![0.0]);

    // Public matrices of the wrong size are detected instead of indexing out of range.
    let mut broken = single_mode(1.0, 100.0, 0.0);
    broken.reduced_mass = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    invalid(broken.step_newmark(&[0.0], 0.01, BETA, GAMMA));
}

#[test]
fn initial_state_uses_the_equation_of_motion() {
    // M = 2, K = 8, C = 1:  a0 = (f0 - C v0 - K q0) / M = (3 - 0.25 - 4) / 2.
    let mut rom = single_mode(2.0, 8.0, 1.0);
    rom.set_initial_state(&[0.5], &[0.25], &[3.0])
        .expect("initial state");
    assert_eq!(rom.reduced_displacement(), &[0.5]);
    assert_eq!(rom.reduced_velocity(), &[0.25]);
    assert!((rom.reduced_acceleration()[0] + 0.625).abs() < 1e-14);

    assert!(matches!(
        rom.set_initial_state(&[0.5, 0.0], &[0.25], &[3.0]),
        Err(FemError::InvalidInput(_))
    ));
    let mut massless = single_mode(0.0, 8.0, 0.0);
    assert_eq!(
        massless.set_initial_state(&[0.5], &[0.0], &[0.0]),
        Err(FemError::SingularSystem)
    );
}

#[test]
fn empty_basis_gives_an_empty_displacement() {
    let basis = PodBasis {
        modes: Vec::new(),
        singular_values: Vec::new(),
        n_modes: 0,
    };
    let mut rom = RomSystem::project_system(&basis, &[], &[]);
    let u = rom
        .step_newmark(&[], 0.01, BETA, GAMMA)
        .expect("empty system");
    assert!(u.is_empty());
}
