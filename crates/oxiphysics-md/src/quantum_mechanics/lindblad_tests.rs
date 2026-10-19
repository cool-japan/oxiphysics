// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Lindblad dissipator of [`DensityMatrix`].
//!
//! The pre-fix `lindblad_dissipation(gamma, k)` wrote only row `k` of
//! `−½{L†L, ρ}` (the column was never touched, so `D_kj ≠ conj(D_jk)`) and added
//! the *rate* to `ρ` without multiplying by a time step.  All quantities below
//! are in SI units, like the rest of the density-matrix code (`ħ` in J·s).

use super::*;

/// Largest violation `|ρ_ij − conj(ρ_ji)|` of Hermiticity.
fn hermiticity_defect(rho: &[Vec<(f64, f64)>]) -> f64 {
    let mut worst = 0.0_f64;
    for (i, row) in rho.iter().enumerate() {
        for (j, &(re_ij, im_ij)) in row.iter().enumerate() {
            let (re_ji, im_ji) = rho[j][i];
            worst = worst.max((re_ij - re_ji).abs()).max((im_ij + im_ji).abs());
        }
    }
    worst
}

/// Pure state `|ψ⟩ = Σ c_i |i⟩` (normalised here) as the density matrix `ρ_ij = c_i c_j*`.
fn pure_state(amplitudes: &[(f64, f64)]) -> DensityMatrix {
    let norm = amplitudes
        .iter()
        .map(|&(re, im)| re * re + im * im)
        .sum::<f64>()
        .sqrt();
    let c: Vec<(f64, f64)> = amplitudes
        .iter()
        .map(|&(re, im)| (re / norm, im / norm))
        .collect();
    let n = c.len();
    let mut dm = DensityMatrix::ground_state(n);
    for i in 0..n {
        for j in 0..n {
            let (a_re, a_im) = c[i];
            let (b_re, b_im) = c[j];
            dm.rho[i][j] = (a_re * b_re + a_im * b_im, a_im * b_re - a_re * b_im);
        }
    }
    dm
}

/// Four-level pure state with unequal weights and complex phases on every coherence.
fn coherent_four_level_state() -> DensityMatrix {
    pure_state(&[(0.35, 0.0), (0.4, 0.3), (-0.3, 0.45), (0.2, -0.5)])
}

/// Two-level system prepared in the excited state `|1⟩⟨1|`.
fn excited_two_level_state() -> DensityMatrix {
    let mut dm = DensityMatrix::ground_state(2);
    dm.rho[0][0] = (0.0, 0.0);
    dm.rho[1][1] = (1.0, 0.0);
    dm
}

fn assert_close(actual: (f64, f64), expected: (f64, f64), tol: f64, what: &str) {
    assert!(
        (actual.0 - expected.0).abs() <= tol && (actual.1 - expected.1).abs() <= tol,
        "{what}: got {actual:?}, expected {expected:?}"
    );
}

#[test]
fn excited_state_population_decays_as_exp_minus_gamma_t() {
    // γ = 1e13 s⁻¹, 300 steps of 1 fs: ρ_11(t) = e^{−γ t}, ρ_00 = 1 − ρ_11.
    let gamma = 1.0e13;
    let dt = 1.0e-15;
    let mut dm = excited_two_level_state();
    for step in 1..=300 {
        dm.lindblad_dissipation(gamma, 1, dt);
        let t = step as f64 * dt;
        let survival = (-gamma * t).exp();
        assert!(
            (dm.population(1) - survival).abs() < 1e-13,
            "step {step}: ρ_11 = {}, e^(−γt) = {survival}",
            dm.population(1)
        );
        assert!((dm.population(0) - (1.0 - survival)).abs() < 1e-13);
        assert!((dm.trace() - 1.0).abs() < 1e-14);
    }
    assert!(hermiticity_defect(&dm.rho) < 1e-15);
    // e^{−3} after three lifetimes.
    assert!((dm.population(1) - (-3.0_f64).exp()).abs() < 1e-12);
}

#[test]
fn exact_dissipator_steps_compose() {
    // The exact map is a semigroup: one step of 2 dt equals two steps of dt.
    let gamma = 4.0e12;
    let mut once = coherent_four_level_state();
    let mut twice = once.clone();
    once.lindblad_dissipation(gamma, 2, 8.0e-14);
    twice.lindblad_dissipation(gamma, 2, 4.0e-14);
    twice.lindblad_dissipation(gamma, 2, 4.0e-14);
    for i in 0..4 {
        for j in 0..4 {
            assert_close(once.rho[i][j], twice.rho[i][j], 1e-14, "composition");
        }
    }
}

#[test]
fn dissipator_keeps_complex_coherences_hermitian() {
    // Complex coherences are what exposes the old row-only update: D_kj and D_jk
    // were then not complex conjugates of each other.
    let original = coherent_four_level_state();
    assert!(hermiticity_defect(&original.rho) < 1e-16);
    let gamma_dt = 0.7;
    let mut dm = original.clone();
    dm.lindblad_dissipation(gamma_dt, 2, 1.0);
    assert!(
        hermiticity_defect(&dm.rho) < 1e-15,
        "ρ is no longer Hermitian"
    );
    assert!((dm.trace() - 1.0).abs() < 1e-15);

    let population_decay = (-gamma_dt).exp();
    let coherence_decay = (-0.5 * gamma_dt).exp();
    // Decaying state: ρ_22 shrinks, the loss lands on ρ_00.
    assert!((dm.rho[2][2].0 - original.rho[2][2].0 * population_decay).abs() < 1e-15);
    let gained = original.rho[2][2].0 * (1.0 - population_decay);
    assert!((dm.rho[0][0].0 - (original.rho[0][0].0 + gained)).abs() < 1e-15);
    // Every coherence of state 2 (row and column) decays with half the rate.
    for j in [0_usize, 1, 3] {
        let (re, im) = original.rho[2][j];
        assert_close(
            dm.rho[2][j],
            (re * coherence_decay, im * coherence_decay),
            1e-15,
            "row 2",
        );
        let (re, im) = original.rho[j][2];
        assert_close(
            dm.rho[j][2],
            (re * coherence_decay, im * coherence_decay),
            1e-15,
            "column 2",
        );
    }
    // Elements not involving state 2 or the arrival state 0 are untouched.
    for &(i, j) in &[
        (1_usize, 1_usize),
        (1, 3),
        (3, 1),
        (3, 3),
        (0, 1),
        (1, 0),
        (0, 3),
    ] {
        assert_close(dm.rho[i][j], original.rho[i][j], 1e-16, "spectator element");
    }
}

#[test]
fn lowest_state_channel_is_pure_dephasing() {
    // L = |0⟩⟨0|: populations are conserved, only coherences with state 0 decay.
    let original = coherent_four_level_state();
    let mut dm = original.clone();
    let gamma_dt = 1.2;
    dm.lindblad_dissipation(gamma_dt, 0, 1.0);
    let decay = (-0.5 * gamma_dt).exp();
    for i in 0..4 {
        assert_eq!(dm.rho[i][i], original.rho[i][i], "population {i} changed");
    }
    for j in 1..4 {
        let (re, im) = original.rho[0][j];
        assert_close(dm.rho[0][j], (re * decay, im * decay), 1e-15, "ρ_0j");
        let (re, im) = original.rho[j][0];
        assert_close(dm.rho[j][0], (re * decay, im * decay), 1e-15, "ρ_j0");
    }
    for &(i, j) in &[(1_usize, 2_usize), (2, 1), (1, 3), (3, 2)] {
        assert_eq!(dm.rho[i][j], original.rho[i][j], "ρ_{i}{j} changed");
    }
    assert!(hermiticity_defect(&dm.rho) < 1e-15);
}

#[test]
fn very_large_gamma_dt_stays_physical() {
    // A first-order Euler update with γ dt = 50 would give ρ_11 = −49.
    let mut dm = excited_two_level_state();
    dm.lindblad_dissipation(50.0, 1, 1.0);
    assert!(dm.population(1) >= 0.0 && dm.population(1) < 1e-20);
    assert!((dm.population(0) - 1.0).abs() < 1e-15);
    assert!((dm.trace() - 1.0).abs() < 1e-15);
}

#[test]
fn zero_time_step_and_invalid_index_leave_the_state_unchanged() {
    let original = coherent_four_level_state();
    let mut dm = original.clone();
    dm.lindblad_dissipation(1.0e13, 2, 0.0);
    assert_eq!(dm.rho, original.rho);
    dm.lindblad_dissipation(1.0e13, 4, 1.0e-12);
    assert_eq!(dm.rho, original.rho);
    assert!(
        dm.lindblad_rate(1.0, 9)
            .iter()
            .flatten()
            .all(|&element| element == (0.0, 0.0))
    );
}

#[test]
fn rate_matrix_is_hermitian_traceless_and_has_the_textbook_elements() {
    let dm = coherent_four_level_state();
    let gamma = 2.5;
    let k = 1;
    let d = dm.lindblad_rate(gamma, k);
    assert!(hermiticity_defect(&d) < 1e-15, "D must be Hermitian");
    let trace: f64 = (0..4).map(|i| d[i][i].0).sum();
    assert!(trace.abs() < 1e-15, "Tr D = {trace}");
    // D_00 = γ ρ_11, D_11 = −γ ρ_11 (gain and loss of the decaying state).
    assert!((d[0][0].0 - gamma * dm.rho[1][1].0).abs() < 1e-15);
    assert!((d[1][1].0 + gamma * dm.rho[1][1].0).abs() < 1e-15);
    // The column element that the old code never wrote: D_j1 = −½ γ ρ_j1.
    for j in [0_usize, 2, 3] {
        let (re, im) = dm.rho[j][1];
        assert_close(
            d[j][1],
            (-0.5 * gamma * re, -0.5 * gamma * im),
            1e-15,
            "D_j1",
        );
        let (re, im) = dm.rho[1][j];
        assert_close(
            d[1][j],
            (-0.5 * gamma * re, -0.5 * gamma * im),
            1e-15,
            "D_1j",
        );
    }
    // Everything away from row/column 1 and from the (0,0) gain term vanishes.
    for &(i, j) in &[
        (2_usize, 2_usize),
        (2, 3),
        (3, 2),
        (3, 3),
        (0, 2),
        (2, 0),
        (0, 3),
    ] {
        assert_eq!(d[i][j], (0.0, 0.0), "D_{i}{j}");
    }
}

#[test]
fn rate_matrix_is_the_time_derivative_of_the_exact_step() {
    let original = coherent_four_level_state();
    let gamma = 3.0;
    for k in 0..4 {
        let rate = original.lindblad_rate(gamma, k);
        let dt = 1.0e-7;
        let mut stepped = original.clone();
        stepped.lindblad_dissipation(gamma, k, dt);
        for (i, (rate_row, (stepped_row, original_row))) in rate
            .iter()
            .zip(stepped.rho.iter().zip(original.rho.iter()))
            .enumerate()
        {
            for (j, (&expected, (&after, &before))) in rate_row
                .iter()
                .zip(stepped_row.iter().zip(original_row.iter()))
                .enumerate()
            {
                let numerical = ((after.0 - before.0) / dt, (after.1 - before.1) / dt);
                assert_close(numerical, expected, 1e-5, &format!("k = {k}, ({i},{j})"));
            }
        }
    }
}

/// Two-level Hamiltonian `ħΩ/2 (|0⟩⟨1| + |1⟩⟨0|)` of a resonantly driven transition.
fn driven_two_level_hamiltonian(rabi_frequency: f64) -> HamiltonianMatrix {
    let mut h = HamiltonianMatrix::zeros(2);
    h.add_coupling(0, 1, 0.5 * HBAR * rabi_frequency);
    h
}

#[test]
fn evolve_lindblad_pure_decay_matches_the_exponential() {
    let h = HamiltonianMatrix::zeros(2);
    let gamma = 1.0e13;
    let dt = 1.0e-15;
    let mut dm = excited_two_level_state();
    let steps = 300;
    for _ in 0..steps {
        dm.evolve_lindblad(&h, &[(gamma, 1)], dt);
    }
    // Three lifetimes: e^{-3} = 0.0498.
    let survival = (-gamma * steps as f64 * dt).exp();
    assert!((survival - (-3.0_f64).exp()).abs() < 1e-12);
    assert!(
        (dm.population(1) - survival).abs() < 1e-9,
        "ρ_11 = {}, e^(−γt) = {survival}",
        dm.population(1)
    );
    assert!((dm.trace() - 1.0).abs() < 1e-12);
    assert!(hermiticity_defect(&dm.rho) < 1e-12);
}

#[test]
fn evolve_lindblad_matches_exact_dissipator_for_a_single_small_step() {
    // With H = 0 one RK4 step agrees with the exact map up to its O((γ dt)^5) error.
    let h = HamiltonianMatrix::zeros(4);
    let gamma = 1.0e13;
    let dt = 5.0e-16;
    let mut rk4 = coherent_four_level_state();
    let mut exact = rk4.clone();
    rk4.evolve_lindblad(&h, &[(gamma, 3)], dt);
    exact.lindblad_dissipation(gamma, 3, dt);
    for i in 0..4 {
        for j in 0..4 {
            assert_close(rk4.rho[i][j], exact.rho[i][j], 1e-12, "RK4 vs exact");
        }
    }
}

#[test]
fn driven_decaying_two_level_atom_reaches_the_optical_bloch_steady_state() {
    // On resonance: ρ_11(∞) = (Ω²/4) / (γ²/4 + Ω²/2) = Ω² / (γ² + 2Ω²).
    let gamma = 1.0e13;
    let rabi = 5.0e12;
    let h = driven_two_level_hamiltonian(rabi);
    let mut dm = DensityMatrix::ground_state(2);
    let dt = 2.0e-16;
    for _ in 0..20_000 {
        dm.evolve_lindblad(&h, &[(gamma, 1)], dt);
    }
    let expected = rabi * rabi / (gamma * gamma + 2.0 * rabi * rabi);
    assert!(
        (dm.population(1) - expected).abs() < 1e-6,
        "steady-state ρ_11 = {}, expected {expected}",
        dm.population(1)
    );
    assert!((dm.population(0) - (1.0 - expected)).abs() < 1e-6);
    assert!((dm.trace() - 1.0).abs() < 1e-10);
    assert!(hermiticity_defect(&dm.rho) < 1e-10);
    // Dissipation makes the steady state mixed: purity strictly below one.
    assert!(dm.purity() < 0.99);
}

#[test]
fn undamped_driven_two_level_atom_performs_rabi_oscillations() {
    // Without decay the coherent drive alone transfers the population:
    // ρ_11(t) = sin²(Ω t / 2).  `evolve_unitary` is `evolve_lindblad` without channels.
    let rabi = 4.0e12;
    let h = driven_two_level_hamiltonian(rabi);
    let mut with_lindblad = DensityMatrix::ground_state(2);
    let mut unitary = DensityMatrix::ground_state(2);
    let dt = 5.0e-16;
    let steps = 1000;
    for _ in 0..steps {
        with_lindblad.evolve_lindblad(&h, &[], dt);
        unitary.evolve_unitary(&h, dt);
    }
    let t = steps as f64 * dt;
    let expected = (0.5 * rabi * t).sin().powi(2);
    assert!(expected > 0.5, "the test should cover a sizeable transfer");
    for dm in [&with_lindblad, &unitary] {
        assert!(
            (dm.population(1) - expected).abs() < 1e-8,
            "ρ_11 = {}, sin² = {expected}",
            dm.population(1)
        );
        assert!((dm.trace() - 1.0).abs() < 1e-12);
        // A pure state stays pure under unitary dynamics.
        assert!((dm.purity() - 1.0).abs() < 1e-8);
    }
    assert_eq!(with_lindblad.rho, unitary.rho);
}

#[test]
fn coupled_three_level_cascade_keeps_trace_hermiticity_and_positivity() {
    // |2> --γ2--> |0> and |1> --γ1--> |0> with a coherent coupling 1-2: two channels
    // acting on a state with complex coherences.
    let mut h = HamiltonianMatrix::diagonal(&[0.0, 6.0e-21, 7.5e-21]);
    h.add_coupling(1, 2, 2.0e-21);
    let channels = [(4.0e12, 1_usize), (7.0e12, 2_usize)];
    let mut dm = pure_state(&[(0.1, 0.0), (0.5, 0.3), (-0.4, 0.6)]);
    let dt = 1.0e-16;
    for step in 0..40_000 {
        dm.evolve_lindblad(&h, &channels, dt);
        if step % 5000 == 0 {
            assert!(
                (dm.trace() - 1.0).abs() < 1e-12,
                "step {step}: {}",
                dm.trace()
            );
            assert!(hermiticity_defect(&dm.rho) < 1e-12, "step {step}");
            assert!(
                (0..3).all(|i| dm.population(i) > -1e-9 && dm.population(i) < 1.0 + 1e-9),
                "step {step}: populations {:?}",
                (0..3).map(|i| dm.population(i)).collect::<Vec<_>>()
            );
        }
    }
    // t = 4 ps = 16 / γ_min: everything has decayed into |0> (e^{-16} ~ 1e-7).
    assert!(dm.population(0) > 1.0 - 1e-5, "ρ_00 = {}", dm.population(0));
    assert!((dm.trace() - 1.0).abs() < 1e-12);
    assert!(hermiticity_defect(&dm.rho) < 1e-12);
}
