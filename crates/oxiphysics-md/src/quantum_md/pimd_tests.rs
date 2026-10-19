// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the path-integral kinetic-energy estimators and the
//! ring-polymer integrator of [`super`].
//!
//! The estimators are checked on *normal-mode snapshots*: bead positions whose
//! normal-mode amplitudes all equal their canonical rms value in a harmonic
//! external potential.  Both estimators are quadratic forms that are diagonal in
//! the normal modes, so such a snapshot reproduces their canonical average
//! exactly and the test is deterministic (no sampling noise).  With the
//! convention of the module (full potential on every bead, `β_P = β / P`) the
//! variance of mode `n` of a particle in `V = m ω² r² / 2` is
//! `1 / (β_P m (ω_n² + ω²))`.

use super::*;

/// Positions of a ring polymer whose normal-mode amplitudes are the canonical rms
/// values in `V = m ω² r² / 2` (`omega = 0`: free ring polymer, zero centroid).
fn canonical_snapshot(
    n_beads: usize,
    mass: f64,
    beta: f64,
    hbar: f64,
    omega: f64,
) -> Vec<[f64; 3]> {
    let transform = NormalModeTransform::new(n_beads);
    let modes: Vec<[f64; 3]> = (0..n_beads)
        .map(|n| {
            let w_n = mode_frequency(n, n_beads, beta, hbar);
            let stiffness = w_n * w_n + omega * omega;
            if stiffness == 0.0 {
                return [0.0; 3];
            }
            let amplitude = (n_beads as f64 / (beta * mass * stiffness)).sqrt();
            [amplitude; 3]
        })
        .collect();
    transform.inverse(&modes)
}

/// External forces `−m ω² r` of per-atom harmonic potentials.
fn harmonic_forces(sys: &QuantumSystem, omegas: &[f64]) -> Vec<Vec<[f64; 3]>> {
    sys.beads
        .iter()
        .zip(omegas.iter())
        .map(|(bead, &omega)| {
            bead.positions
                .iter()
                .map(|&r| v3_scale(-bead.mass * omega * omega, r))
                .collect()
        })
        .collect()
}

/// Exact canonical kinetic energy `3 (ℏω/4) coth(βℏω/2)` of a 3-D oscillator.
fn exact_oscillator_kinetic_energy(beta: f64, hbar: f64, omega: f64) -> f64 {
    0.75 * hbar * omega / (0.5 * beta * hbar * omega).tanh()
}

/// Finite-`P` canonical average of both estimators for one 3-D oscillator:
/// `3/(2β) [1 + Σ_{n=1}^{P−1} ω² / (ω_n² + ω²)]`.
fn finite_p_oscillator_kinetic_energy(n_beads: usize, beta: f64, hbar: f64, omega: f64) -> f64 {
    let correction: f64 = (1..n_beads)
        .map(|n| {
            let w_n = mode_frequency(n, n_beads, beta, hbar);
            omega * omega / (w_n * w_n + omega * omega)
        })
        .sum();
    3.0 / (2.0 * beta) * (1.0 + correction)
}

#[test]
fn primitive_estimator_of_a_free_ring_polymer_is_the_classical_free_value() {
    // Free particle: <T> = 3/(2β) for every P.  The pre-fix spring term was P times
    // too large, which made the estimator negative already for a handful of beads.
    let (beta, mass, hbar) = (4.0, 1.7, 1.3);
    for &p in &[2_usize, 3, 8, 33] {
        let mut sys = QuantumSystem::new(1, p, beta, mass, hbar);
        sys.beads[0].positions = canonical_snapshot(p, mass, beta, hbar, 0.0);
        let t_prim = primitive_kinetic_energy(&sys);
        assert!(
            (t_prim - 1.5 / beta).abs() < 1e-12,
            "P = {p}: T_prim = {t_prim}, expected {}",
            1.5 / beta
        );
    }
}

#[test]
fn both_estimators_reproduce_the_harmonic_oscillator_kinetic_energy() {
    // Two atoms with different mass, ℏ_eff and frequency: beta*hbar*omega = 4 and 3.36.
    let p = 32;
    let beta = 4.0;
    let params = [(1.0, 1.0, 1.0), (2.5, 1.4, 0.6)]; // (mass, hbar, omega)
    let mut sys = QuantumSystem::new(2, p, beta, 1.0, 1.0);
    for (bead, &(mass, hbar, omega)) in sys.beads.iter_mut().zip(params.iter()) {
        bead.mass = mass;
        bead.hbar_eff = hbar;
        bead.positions = canonical_snapshot(p, mass, beta, hbar, omega);
    }
    let omegas = [params[0].2, params[1].2];
    let forces = harmonic_forces(&sys, &omegas);
    let t_virial = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
    let t_primitive = primitive_kinetic_energy(&sys);

    let finite_p: f64 = params
        .iter()
        .map(|&(_, hbar, omega)| finite_p_oscillator_kinetic_energy(p, beta, hbar, omega))
        .sum();
    let exact: f64 = params
        .iter()
        .map(|&(_, hbar, omega)| exact_oscillator_kinetic_energy(beta, hbar, omega))
        .sum();
    assert!(
        (t_virial - finite_p).abs() < 1e-10 * finite_p,
        "virial {t_virial} vs finite-P average {finite_p}"
    );
    assert!(
        (t_primitive - finite_p).abs() < 1e-10 * finite_p,
        "primitive {t_primitive} vs finite-P average {finite_p}"
    );
    // P = 32 is converged to a fraction of a percent against the exact quantum value.
    assert!(
        (t_virial - exact).abs() / exact < 0.01,
        "{t_virial} vs {exact}"
    );
    assert!(
        (t_primitive - exact).abs() / exact < 0.01,
        "{t_primitive} vs {exact}"
    );
    // Zero-point motion: far above the classical equipartition value 3N/(2β).
    let classical = 3.0 * 2.0 / (2.0 * beta);
    assert!(
        t_virial > 1.6 * classical,
        "{t_virial} vs classical {classical}"
    );
}

#[test]
fn estimators_converge_to_the_quantum_value_as_the_bead_number_grows() {
    let (beta, mass, hbar, omega) = (4.0, 1.0, 1.0, 1.0);
    let exact = exact_oscillator_kinetic_energy(beta, hbar, omega);
    let mut previous_error = f64::INFINITY;
    for &p in &[4_usize, 8, 16, 32, 64] {
        let mut sys = QuantumSystem::new(1, p, beta, mass, hbar);
        sys.beads[0].positions = canonical_snapshot(p, mass, beta, hbar, omega);
        let forces = harmonic_forces(&sys, &[omega]);
        let t_virial = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
        let error = (t_virial - exact).abs();
        assert!(
            error < previous_error,
            "P = {p}: error {error} did not shrink"
        );
        // The Trotter error falls like 1/P²: doubling P cuts it by about four.
        if previous_error.is_finite() {
            assert!(
                previous_error / error > 3.0,
                "P = {p}: ratio {}",
                previous_error / error
            );
        }
        previous_error = error;
    }
    assert!(previous_error / exact < 2e-3);
}

#[test]
fn estimators_agree_for_odd_and_tiny_bead_numbers() {
    // The normal-mode frequencies and the snapshot construction hold for every P,
    // including P = 1, 2 and odd P (cosine and sine modes pair up as n <-> P - n).
    let (beta, mass, hbar, omega) = (3.0, 1.3, 1.1, 0.8);
    for &p in &[1_usize, 2, 3, 5, 7, 15, 33] {
        let mut sys = QuantumSystem::new(1, p, beta, mass, hbar);
        sys.beads[0].positions = canonical_snapshot(p, mass, beta, hbar, omega);
        let forces = harmonic_forces(&sys, &[omega]);
        let t_virial = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
        let t_primitive = primitive_kinetic_energy(&sys);
        let finite_p = finite_p_oscillator_kinetic_energy(p, beta, hbar, omega);
        assert!(
            (t_virial - finite_p).abs() < 1e-10 * finite_p,
            "P = {p}: virial {t_virial} vs {finite_p}"
        );
        assert!(
            (t_primitive - finite_p).abs() < 1e-10 * finite_p,
            "P = {p}: primitive {t_primitive} vs {finite_p}"
        );
    }
}

#[test]
fn virial_estimator_is_translation_invariant_and_odd_in_the_force() {
    let (p, beta, mass, hbar, omega) = (16, 3.0, 1.2, 1.0, 0.9);
    let mut sys = QuantumSystem::new(1, p, beta, mass, hbar);
    sys.beads[0].positions = canonical_snapshot(p, mass, beta, hbar, omega);
    let forces = harmonic_forces(&sys, &[omega]);
    let reference = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");

    // A rigid translation changes the centroid but not the bead spread.
    for pos in sys.beads[0].positions.iter_mut() {
        *pos = v3_add(*pos, [0.7, -1.1, 0.4]);
    }
    let shifted_forces = harmonic_forces(&sys, &[omega]);
    let shifted = quantum_kinetic_energy(&sys, &shifted_forces).expect("forces match the system");
    assert!((shifted - reference).abs() < 1e-10 * reference);

    // The correction term changes sign with the force.
    let classical = 1.5 / beta;
    let negated: Vec<Vec<[f64; 3]>> = shifted_forces
        .iter()
        .map(|bead| bead.iter().map(|&f| v3_scale(-1.0, f)).collect())
        .collect();
    let flipped = quantum_kinetic_energy(&sys, &negated).expect("forces match the system");
    assert!(((flipped - classical) + (shifted - classical)).abs() < 1e-10 * reference);
    assert!(shifted > classical && flipped < classical);
}

#[test]
fn virial_estimator_rejects_forces_of_the_wrong_shape() {
    let sys = QuantumSystem::new(2, 4, 1.0, 1.0, 1.0);
    let too_few_atoms = vec![vec![[0.0; 3]; 4]];
    assert!(quantum_kinetic_energy(&sys, &too_few_atoms).is_err());
    let too_few_beads = vec![vec![[0.0; 3]; 3], vec![[0.0; 3]; 4]];
    assert!(quantum_kinetic_energy(&sys, &too_few_beads).is_err());
    assert!(QuantumDiagnostics::from_system(&sys, &too_few_beads).is_err());
    let matching = vec![vec![[0.0; 3]; 4]; 2];
    assert!(quantum_kinetic_energy(&sys, &matching).is_ok());
}

#[test]
fn spring_quantities_use_each_atoms_own_effective_hbar() {
    // ring_polymer_spring_constant ignored hbar_eff (hard-coded ℏ = 1).
    let positions = [
        [0.3, 0.0, 0.0],
        [0.0, 0.4, 0.0],
        [-0.2, 0.1, 0.0],
        [0.0, -0.3, 0.5],
    ];
    let single_atom = |hbar: f64| {
        let mut sys = QuantumSystem::new(1, 4, 1.5, 1.0, hbar);
        sys.beads[0].positions = positions.to_vec();
        sys
    };
    let mut pair = QuantumSystem::new(2, 4, 1.5, 1.0, 1.0);
    pair.beads[1].hbar_eff = 2.0;
    for bead in pair.beads.iter_mut() {
        bead.positions = positions.to_vec();
    }

    // ℏ = 2 gives springs four times softer: E_pair = E(ℏ=1) (1 + 1/4).
    let unit_energy = ring_polymer_spring_energy(&single_atom(1.0));
    let pair_energy = ring_polymer_spring_energy(&pair);
    assert!((pair_energy - 1.25 * unit_energy).abs() < 1e-12 * unit_energy);

    // The primitive estimator's spring term (classical term minus T_prim) scales alike.
    let spring_term = |sys: &QuantumSystem| {
        3.0 * sys.n_atoms as f64 * 4.0 / (2.0 * sys.beta) - primitive_kinetic_energy(sys)
    };
    let unit_term = spring_term(&single_atom(1.0));
    assert!((spring_term(&pair) - 1.25 * unit_term).abs() < 1e-12 * unit_term.abs());

    // Direct value for one atom: Σ_links k / (2P) |Δr|² with k = m (P / (β ℏ))².
    let k = ring_polymer_spring_constant(1.0, 4, 1.5, 1.0);
    let link_sum: f64 = (0..4)
        .map(|j| v3_norm_sq(v3_sub(positions[(j + 1) % 4], positions[j])))
        .sum();
    assert!((unit_term - k / (2.0 * 4.0) * link_sum).abs() < 1e-12 * unit_term.abs());
}

// ============================================================================
// Integrator
// ============================================================================

/// Deterministic, non-trivial initial condition for `n_atoms` atoms of different masses.
fn displaced_system(n_atoms: usize, n_beads: usize, beta: f64) -> QuantumSystem {
    let mut sys = QuantumSystem::new(n_atoms, n_beads, beta, 1.3, 1.0);
    for (i, bead) in sys.beads.iter_mut().enumerate() {
        bead.mass = 1.3 + 0.4 * i as f64;
        for j in 0..n_beads {
            let x = (j + 3 * i) as f64;
            bead.positions[j] = [
                0.3 * (0.7 * x + 0.2).sin() + 0.5,
                -0.2 * (1.3 * x).cos(),
                0.1 * (0.4 * x + 1.0).sin(),
            ];
            bead.velocities[j] = [
                0.4 * (0.9 * x).cos(),
                0.3 * (0.5 * x + 0.3).sin(),
                -0.2 * (1.1 * x).cos(),
            ];
        }
    }
    sys
}

/// `H_P = Σ m v²/2 + spring energy + Σ_p V(m, r_p)` for the ring-polymer Hamiltonian.
fn ring_polymer_hamiltonian(sys: &QuantumSystem, potential: impl Fn(f64, [f64; 3]) -> f64) -> f64 {
    let kinetic: f64 = sys
        .beads
        .iter()
        .map(|b| {
            b.velocities
                .iter()
                .map(|&v| 0.5 * b.mass * v3_norm_sq(v))
                .sum::<f64>()
        })
        .sum();
    let external: f64 = sys
        .beads
        .iter()
        .map(|b| {
            b.positions
                .iter()
                .map(|&r| potential(b.mass, r))
                .sum::<f64>()
        })
        .sum();
    kinetic + ring_polymer_spring_energy(sys) + external
}

/// Fill `forces` with `f(m, r)` for every bead.
fn fill_forces(
    sys: &QuantumSystem,
    forces: &mut [Vec<[f64; 3]>],
    force: impl Fn(f64, [f64; 3]) -> [f64; 3],
) {
    for (bead, bead_forces) in sys.beads.iter().zip(forces.iter_mut()) {
        for (&r, f) in bead.positions.iter().zip(bead_forces.iter_mut()) {
            *f = force(bead.mass, r);
        }
    }
}

fn harmonic_potential(omega: f64) -> impl Fn(f64, [f64; 3]) -> f64 {
    move |mass, r| 0.5 * mass * omega * omega * v3_norm_sq(r)
}

fn harmonic_force(omega: f64) -> impl Fn(f64, [f64; 3]) -> [f64; 3] {
    move |mass, r| v3_scale(-mass * omega * omega, r)
}

fn quartic_potential(kappa: f64) -> impl Fn(f64, [f64; 3]) -> f64 {
    move |_, r| 0.25 * kappa * v3_norm_sq(r).powi(2)
}

fn quartic_force(kappa: f64) -> impl Fn(f64, [f64; 3]) -> [f64; 3] {
    move |_, r| v3_scale(-kappa * v3_norm_sq(r), r)
}

/// Run `steps` steps of size `dt`; returns the largest relative energy deviation.
fn run_and_track_energy(
    sys: &mut QuantumSystem,
    steps: usize,
    dt: f64,
    potential: impl Fn(f64, [f64; 3]) -> f64,
    force: impl Fn(f64, [f64; 3]) -> [f64; 3],
) -> f64 {
    let mut forces = vec![vec![[0.0; 3]; sys.n_beads]; sys.n_atoms];
    fill_forces(sys, &mut forces, &force);
    let initial = ring_polymer_hamiltonian(sys, &potential);
    let mut worst = 0.0_f64;
    for _ in 0..steps {
        step_pimd(sys, &mut forces, dt, |s, f| fill_forces(s, f, &force))
            .expect("force shape matches the system");
        let energy = ring_polymer_hamiltonian(sys, &potential);
        worst = worst.max((energy - initial).abs() / initial.abs());
    }
    worst
}

#[test]
fn free_ring_polymer_conserves_energy_for_arbitrarily_large_steps() {
    // The springs are propagated exactly in normal-mode space, so even a step with
    // ω_max dt = 5 (a Verlet step on the springs would be unstable beyond 2) keeps the
    // free ring-polymer energy to rounding error.
    let beta = 2.0;
    let mut sys = displaced_system(2, 8, beta);
    let omega_max = mode_frequency(4, 8, beta, 1.0);
    let dt = 5.0 / omega_max;
    let worst = run_and_track_energy(&mut sys, 2000, dt, |_, _| 0.0, |_, _| [0.0; 3]);
    assert!(worst < 1e-11, "free energy drift {worst}");
}

#[test]
fn free_step_is_the_exact_normal_mode_rotation() {
    let (p, beta, dt) = (6, 1.7, 0.37);
    let mut sys = displaced_system(1, p, beta);
    let transform = NormalModeTransform::new(p);
    let q0 = transform.forward(&sys.beads[0].positions);
    let v0 = transform.forward(&sys.beads[0].velocities);
    let frequencies = normal_mode_frequencies(p, beta, sys.beads[0].hbar_eff);
    let mut q1 = vec![[0.0; 3]; p];
    let mut v1 = vec![[0.0; 3]; p];
    for n in 0..p {
        for d in 0..3 {
            let w = frequencies[n];
            if w == 0.0 {
                q1[n][d] = q0[n][d] + dt * v0[n][d];
                v1[n][d] = v0[n][d];
            } else {
                q1[n][d] = q0[n][d] * (w * dt).cos() + v0[n][d] / w * (w * dt).sin();
                v1[n][d] = -q0[n][d] * w * (w * dt).sin() + v0[n][d] * (w * dt).cos();
            }
        }
    }
    let expected_positions = transform.inverse(&q1);
    let expected_velocities = transform.inverse(&v1);
    let mut forces = vec![vec![[0.0; 3]; p]];
    step_pimd(&mut sys, &mut forces, dt, |_, _| {}).expect("force shape matches the system");
    for j in 0..p {
        for d in 0..3 {
            assert!((sys.beads[0].positions[j][d] - expected_positions[j][d]).abs() < 1e-12);
            assert!((sys.beads[0].velocities[j][d] - expected_velocities[j][d]).abs() < 1e-12);
        }
    }
}

#[test]
fn constant_force_accelerates_the_whole_ring_uniformly() {
    // A uniform force acts on the centroid only: relative to the free ring polymer every
    // bead is displaced by f t²/(2m) and every velocity changes by f t/m (exact for
    // constant forces).
    let (p, beta, dt, steps) = (8, 2.0, 0.05, 400);
    let mut forced = displaced_system(2, p, beta);
    let mut free = forced.clone();
    let uniform = [0.3, -0.2, 0.5];
    let mut forces: Vec<Vec<[f64; 3]>> = vec![vec![uniform; p]; 2];
    let mut no_forces = vec![vec![[0.0; 3]; p]; 2];
    for _ in 0..steps {
        step_pimd(&mut forced, &mut forces, dt, |_, f| {
            for bead_forces in f.iter_mut() {
                bead_forces.fill(uniform);
            }
        })
        .expect("force shape matches the system");
        step_pimd(&mut free, &mut no_forces, dt, |_, _| {})
            .expect("force shape matches the system");
    }
    let t = steps as f64 * dt;
    for (bead_forced, bead_free) in forced.beads.iter().zip(free.beads.iter()) {
        for j in 0..p {
            for (d, &force_component) in uniform.iter().enumerate() {
                let shift = 0.5 * force_component / bead_forced.mass * t * t;
                let kick = force_component / bead_forced.mass * t;
                let dx = bead_forced.positions[j][d] - bead_free.positions[j][d];
                let dv = bead_forced.velocities[j][d] - bead_free.velocities[j][d];
                assert!((dx - shift).abs() < 1e-9, "bead {j} d {d}: {dx} vs {shift}");
                assert!((dv - kick).abs() < 1e-9, "bead {j} d {d}: {dv} vs {kick}");
            }
        }
    }
}

#[test]
fn second_half_kick_uses_the_refreshed_forces() {
    // P = 1 degenerates to textbook velocity Verlet for V = k x²/2:
    //   x1 = x0 + dt v0 + dt² a0 / 2,  v1 = v0 + dt (a0 + a1) / 2,  a1 = a(x1).
    // The old code reused a0 in the second kick.
    let (mass, stiffness, dt) = (2.0, 3.0, 0.1);
    let (x0, v0) = (1.0, 0.2);
    let mut sys = QuantumSystem::new(1, 1, 1.0, mass, 1.0);
    sys.beads[0].positions[0] = [x0, 0.0, 0.0];
    sys.beads[0].velocities[0] = [v0, 0.0, 0.0];
    let acceleration = |x: f64| -stiffness * x / mass;
    let mut forces = vec![vec![[-stiffness * x0, 0.0, 0.0]]];
    let mut callbacks = 0;
    let mut position_seen_by_callback = f64::NAN;
    step_pimd(&mut sys, &mut forces, dt, |s, f| {
        callbacks += 1;
        position_seen_by_callback = s.beads[0].positions[0][0];
        f[0][0] = [-stiffness * s.beads[0].positions[0][0], 0.0, 0.0];
    })
    .expect("force shape matches the system");
    let a0 = acceleration(x0);
    let x1 = x0 + dt * v0 + 0.5 * dt * dt * a0;
    let a1 = acceleration(x1);
    let v1 = v0 + 0.5 * dt * (a0 + a1);
    assert_eq!(callbacks, 1, "the force routine must run once per step");
    assert!(
        (position_seen_by_callback - x1).abs() < 1e-15,
        "forces refreshed at old positions"
    );
    assert!((sys.beads[0].positions[0][0] - x1).abs() < 1e-15);
    assert!((sys.beads[0].velocities[0][0] - v1).abs() < 1e-15);
    // The cached forces now belong to the new positions.
    assert!((forces[0][0][0] + stiffness * x1).abs() < 1e-15);
    // The stale second kick would have produced this velocity instead.
    let stale = v0 + dt * a0;
    assert!((sys.beads[0].velocities[0][0] - stale).abs() > 1e-4);
}

#[test]
fn harmonic_external_potential_conserves_the_ring_polymer_hamiltonian() {
    let beta = 2.0;
    let omega = 1.7;
    let mut sys = displaced_system(2, 8, beta);
    let worst = run_and_track_energy(
        &mut sys,
        10_000,
        0.001,
        harmonic_potential(omega),
        harmonic_force(omega),
    );
    assert!(worst < 1e-5, "relative energy deviation {worst}");
}

#[test]
fn quartic_external_potential_conserves_the_ring_polymer_hamiltonian() {
    let beta = 2.0;
    let kappa = 0.8;
    let mut sys = displaced_system(2, 8, beta);
    let worst = run_and_track_energy(
        &mut sys,
        10_000,
        0.002,
        quartic_potential(kappa),
        quartic_force(kappa),
    );
    assert!(worst < 1e-5, "relative energy deviation {worst}");
}

/// Largest bead-position error of `steps` steps of size `dt` against the exact
/// normal-mode solution `q_n(t) = q_n cos Ω_n t + (v_n/Ω_n) sin Ω_n t`,
/// `Ω_n² = ω_n² + ω²`, of the ring polymer in a harmonic potential.
fn harmonic_trajectory_error(steps: usize, dt: f64) -> f64 {
    let (p, beta, omega) = (8, 2.0, 1.7);
    let mut sys = displaced_system(1, p, beta);
    let transform = NormalModeTransform::new(p);
    let q0 = transform.forward(&sys.beads[0].positions);
    let v0 = transform.forward(&sys.beads[0].velocities);
    let frequencies = normal_mode_frequencies(p, beta, sys.beads[0].hbar_eff);
    let mut forces = vec![vec![[0.0; 3]; p]];
    fill_forces(&sys, &mut forces, harmonic_force(omega));
    for _ in 0..steps {
        step_pimd(&mut sys, &mut forces, dt, |s, f| {
            fill_forces(s, f, harmonic_force(omega))
        })
        .expect("force shape matches the system");
    }
    let t = steps as f64 * dt;
    let exact_modes: Vec<[f64; 3]> = frequencies
        .iter()
        .zip(q0.iter().zip(v0.iter()))
        .map(|(&w_n, (q_n, v_n))| {
            let big_omega = (w_n.powi(2) + omega * omega).sqrt();
            let (sin_wt, cos_wt) = (big_omega * t).sin_cos();
            [0, 1, 2].map(|d| q_n[d] * cos_wt + v_n[d] / big_omega * sin_wt)
        })
        .collect();
    let exact = transform.inverse(&exact_modes);
    sys.beads[0]
        .positions
        .iter()
        .zip(exact.iter())
        .flat_map(|(simulated, reference)| {
            simulated
                .iter()
                .zip(reference.iter())
                .map(|(a, b)| (a - b).abs())
        })
        .fold(0.0_f64, f64::max)
}

#[test]
fn harmonic_trajectory_matches_the_normal_mode_solution_with_second_order_accuracy() {
    let coarse = harmonic_trajectory_error(200, 0.02);
    let fine = harmonic_trajectory_error(400, 0.01);
    let finer = harmonic_trajectory_error(800, 0.005);
    assert!(finer < 1e-4, "error at dt = 0.005 is {finer}");
    let ratio_1 = coarse / fine;
    let ratio_2 = fine / finer;
    assert!(
        (3.5..4.5).contains(&ratio_1),
        "error ratio {ratio_1} (expected ~4)"
    );
    assert!(
        (3.5..4.5).contains(&ratio_2),
        "error ratio {ratio_2} (expected ~4)"
    );
}

#[test]
fn integrator_is_time_reversible() {
    // The symmetric splitting (kick - exact free rotation - kick) is time reversible:
    // flipping all velocities and integrating for as long again retraces the trajectory.
    let (beta, dt, steps) = (2.0, 0.01, 300);
    let mut sys = displaced_system(2, 8, beta);
    let initial = sys.clone();
    let mut forces = vec![vec![[0.0; 3]; 8]; 2];
    fill_forces(&sys, &mut forces, quartic_force(0.8));
    for _ in 0..steps {
        step_pimd(&mut sys, &mut forces, dt, |s, f| {
            fill_forces(s, f, quartic_force(0.8))
        })
        .expect("force shape matches the system");
    }
    let turning_point = sys.clone();
    for bead in sys.beads.iter_mut() {
        for v in bead.velocities.iter_mut() {
            *v = v3_scale(-1.0, *v);
        }
    }
    for _ in 0..steps {
        step_pimd(&mut sys, &mut forces, dt, |s, f| {
            fill_forces(s, f, quartic_force(0.8))
        })
        .expect("force shape matches the system");
    }
    let mut drift = 0.0_f64;
    for (end, start) in sys.beads.iter().zip(initial.beads.iter()) {
        for (r_end, r_start) in end.positions.iter().zip(start.positions.iter()) {
            for (a, b) in r_end.iter().zip(r_start.iter()) {
                drift = drift.max((a - b).abs());
            }
        }
        for (v_end, v_start) in end.velocities.iter().zip(start.velocities.iter()) {
            for (a, b) in v_end.iter().zip(v_start.iter()) {
                drift = drift.max((a + b).abs());
            }
        }
    }
    assert!(
        drift < 1e-10,
        "trajectory did not retrace, deviation {drift}"
    );
    // The forward trajectory itself moved far from the start (the test is not vacuous).
    let travelled = turning_point.beads[0].positions[0]
        .iter()
        .zip(initial.beads[0].positions[0].iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        travelled > 0.1,
        "forward trajectory only moved by {travelled}"
    );
}

#[test]
fn mismatched_force_shapes_are_rejected_without_touching_the_system() {
    let mut sys = displaced_system(2, 4, 1.0);
    let before = sys.clone();
    let mut wrong_atoms = vec![vec![[0.0; 3]; 4]];
    assert!(step_pimd(&mut sys, &mut wrong_atoms, 0.01, |_, _| {}).is_err());
    let mut wrong_beads = vec![vec![[0.0; 3]; 4], vec![[0.0; 3]; 3]];
    assert!(step_pimd(&mut sys, &mut wrong_beads, 0.01, |_, _| {}).is_err());
    for (a, b) in sys.beads.iter().zip(before.beads.iter()) {
        assert_eq!(a.positions, b.positions);
        assert_eq!(a.velocities, b.velocities);
    }
    // A force routine that breaks the shape is reported as well.
    let mut forces = vec![vec![[0.0; 3]; 4]; 2];
    let result = step_pimd(&mut sys, &mut forces, 0.01, |_, f| f[1].truncate(2));
    assert!(result.is_err());
}

#[test]
fn normal_mode_cache_follows_a_change_of_the_bead_number() {
    let beta = 1.5;
    let mut sys = displaced_system(1, 4, beta);
    let mut forces = vec![vec![[0.0; 3]; 4]];
    step_pimd(&mut sys, &mut forces, 0.05, |_, _| {}).expect("force shape matches the system");
    // Re-discretise the same particle with six beads.
    let reference = displaced_system(1, 6, beta);
    sys.n_beads = 6;
    sys.beads = reference.beads.clone();
    let e0 = ring_polymer_hamiltonian(&sys, |_, _| 0.0);
    let mut forces = vec![vec![[0.0; 3]; 6]];
    for _ in 0..50 {
        step_pimd(&mut sys, &mut forces, 0.05, |_, _| {}).expect("force shape matches the system");
    }
    let e1 = ring_polymer_hamiltonian(&sys, |_, _| 0.0);
    assert!((e1 - e0).abs() < 1e-11 * e0.abs());
    assert_eq!(sys.beads[0].positions.len(), 6);
}

#[test]
fn rpmd_state_tracks_time_and_conserves_energy_with_refreshed_forces() {
    let (beta, omega, dt) = (2.0, 1.2, 0.002);
    let mut state = RingPolymerMdState::new(1, 8, beta, 1.0, 1.0);
    let reference = displaced_system(1, 8, beta);
    state.system.beads[0].positions = reference.beads[0].positions.clone();
    state.system.beads[0].velocities = reference.beads[0].velocities.clone();
    // Fresh state starts with zero forces; they must be refreshed explicitly.
    assert!(state.forces.iter().flatten().all(|f| *f == [0.0; 3]));
    state.refresh_forces(|s, f| fill_forces(s, f, harmonic_force(omega)));
    assert!(state.forces.iter().flatten().any(|f| f[0] != 0.0));

    let potential = harmonic_potential(omega);
    let initial = ring_polymer_hamiltonian(&state.system, &potential);
    for _ in 0..1500 {
        state
            .step(dt, |s, f| fill_forces(s, f, harmonic_force(omega)))
            .expect("force shape matches the system");
    }
    assert!((state.time - 1500.0 * dt).abs() < 1e-12);
    let final_energy = ring_polymer_hamiltonian(&state.system, &potential);
    assert!((final_energy - initial).abs() < 1e-5 * initial.abs());
    assert!((state.kinetic_energy() + state.spring_energy()) > 0.0);
}

#[test]
fn diagnostics_report_the_virial_estimator() {
    let (p, beta, mass, hbar, omega) = (16, 3.0, 1.0, 1.0, 1.0);
    let mut sys = QuantumSystem::new(1, p, beta, mass, hbar);
    sys.beads[0].positions = canonical_snapshot(p, mass, beta, hbar, omega);
    let forces = harmonic_forces(&sys, &[omega]);
    let diagnostics = QuantumDiagnostics::from_system(&sys, &forces).expect("forces match");
    let virial = quantum_kinetic_energy(&sys, &forces).expect("forces match");
    assert!((diagnostics.quantum_ke - virial).abs() < 1e-15);
    assert!((diagnostics.spring_energy - ring_polymer_spring_energy(&sys)).abs() < 1e-15);
    assert!(diagnostics.quantum_ke > 1.5 / beta);
}
