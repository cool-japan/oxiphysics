// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Asymptotic uncertainty of MBAR free-energy differences.
//!
//! Shirts & Chodera, *J. Chem. Phys.* **129**, 124105 (2008), Appendix D: with the
//! `N × K` weight matrix `W_nk = exp(f_k − u_k(x_n)) / Σ_l N_l exp(f_l − u_l(x_n))`
//! and `N = diag(N_k)`, the asymptotic covariance of the estimated reduced free
//! energies `f_k` is
//!
//! ```text
//! Θ = Wᵀ (I_N − W N Wᵀ)⁺ W
//! ```
//!
//! where `⁺` is the Moore–Penrose pseudo-inverse.  Forming the `N × N` matrix is
//! infeasible for realistic sample counts, so the formula is evaluated in the
//! equivalent `K × K` form: with the eigendecomposition `Wᵀ W = U Σ² Uᵀ`,
//!
//! ```text
//! Θ = U Σ (I_K − Σ Uᵀ N U Σ)⁺ Σ Uᵀ
//! ```
//!
//! The variance of a difference follows from `Var(f_k − f_l) = Θ_kk + Θ_ll − 2 Θ_kl`.
//! The pseudo-inverse is needed because only differences of the `f_k` are
//! defined: `I_K − Σ Uᵀ N U Σ` has one zero eigenvalue (a common shift of all
//! `f_k`), which cancels in every difference.

use crate::error::MdError;
use oxiphysics_core::numerical_methods::symmetric_eigen_n;

/// Eigenvalues of `I_K − Σ Uᵀ N U Σ` below this magnitude are treated as zero.
const MBAR_PSEUDO_INVERSE_CUTOFF: f64 = 1.0e-10;

/// `ln Σ exp(v_i)` over the finite and `−∞` entries of `values` (`−∞` for none).
fn log_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values.iter().map(|&v| (v - max).exp()).sum::<f64>().ln()
}

/// Gram matrix `Wᵀ W` of the MBAR weight matrix, accumulated sample by sample so
/// the `N × K` matrix itself is never stored.
fn mbar_weight_gram(u_kn: &[Vec<f64>], n_k: &[usize], f: &[f64]) -> Result<Vec<Vec<f64>>, MdError> {
    let n_states = u_kn.len();
    let n_total: usize = n_k.iter().sum();
    let ln_n: Vec<f64> = n_k
        .iter()
        .map(|&n| {
            if n > 0 {
                (n as f64).ln()
            } else {
                f64::NEG_INFINITY
            }
        })
        .collect();
    let mut gram = vec![vec![0.0_f64; n_states]; n_states];
    let mut log_numerator = vec![0.0_f64; n_states];
    let mut log_terms = vec![0.0_f64; n_states];
    for sample in 0..n_total {
        for (state, u) in u_kn.iter().map(|row| row[sample]).enumerate() {
            log_numerator[state] = f[state] - u;
            log_terms[state] = ln_n[state] + log_numerator[state];
        }
        if log_numerator.iter().any(|v| v.is_nan()) {
            return Err(MdError::General(format!(
                "mbar_df_uncertainty: NaN reduced potential in sample {sample}"
            )));
        }
        let log_denominator = log_sum_exp(&log_terms);
        if log_denominator == f64::NEG_INFINITY {
            // The sample is impossible in every sampled state: zero weight.
            continue;
        }
        let weights: Vec<f64> = log_numerator
            .iter()
            .map(|&v| (v - log_denominator).exp())
            .collect();
        for a in 0..n_states {
            for b in a..n_states {
                let product = weights[a] * weights[b];
                gram[a][b] += product;
                if a != b {
                    gram[b][a] += product;
                }
            }
        }
    }
    Ok(gram)
}

/// MBAR uncertainty estimate via error propagation (asymptotic standard error).
///
/// Returns the **variance** of the reduced free-energy difference
/// `f_l − f_k` between states `k` and `l` (the square of the standard error, in
/// units of `(k_B T)²` when `u_kn` is the reduced potential `u = β U`),
/// evaluated with the asymptotic covariance of Shirts & Chodera (2008); see the
/// module documentation for the formula.  The result is symmetric in `k` and
/// `l`, vanishes for `k == l`, and does not depend on the additive gauge of `f`.
///
/// # Arguments
/// * `u_kn`  – K×N reduced potential matrix; `u_kn[k][n]` is the reduced
///   potential of sample `n` evaluated at state `k`, with the samples of all
///   states concatenated (`N = Σ n_k`), as for [`mbar_free_energies`].
/// * `n_k`   – samples per state.
/// * `f`     – converged MBAR free energies, e.g. from [`mbar_free_energies`].
/// * `k`, `l` – state indices.
///
/// # Errors
/// [`MdError::General`] if the dimensions are inconsistent (`u_kn`, `n_k`, `f`
/// disagree on the number of states, a row of `u_kn` is shorter than
/// `Σ n_k`), `k` or `l` is not a state, there are no samples, `f` is not
/// finite, or a reduced potential is NaN.
///
/// [`mbar_free_energies`]: super::functions::mbar_free_energies
pub fn mbar_df_uncertainty(
    u_kn: &[Vec<f64>],
    n_k: &[usize],
    f: &[f64],
    k: usize,
    l: usize,
) -> Result<f64, MdError> {
    let n_states = u_kn.len();
    if n_states == 0 {
        return Err(MdError::General(
            "mbar_df_uncertainty: no states in u_kn".into(),
        ));
    }
    if n_k.len() != n_states || f.len() != n_states {
        return Err(MdError::General(format!(
            "mbar_df_uncertainty: u_kn has {n_states} states but n_k has {} and f has {} entries",
            n_k.len(),
            f.len()
        )));
    }
    if k >= n_states || l >= n_states {
        return Err(MdError::General(format!(
            "mbar_df_uncertainty: state indices ({k}, {l}) out of range for {n_states} states"
        )));
    }
    let n_total: usize = n_k.iter().sum();
    if n_total == 0 {
        return Err(MdError::General(
            "mbar_df_uncertainty: no samples (all n_k are zero)".into(),
        ));
    }
    if let Some(short) = u_kn.iter().position(|row| row.len() < n_total) {
        return Err(MdError::General(format!(
            "mbar_df_uncertainty: row {short} of u_kn has {} samples, expected {n_total}",
            u_kn[short].len()
        )));
    }
    if f.iter().any(|v| !v.is_finite()) {
        return Err(MdError::General(
            "mbar_df_uncertainty: free energies f must be finite".into(),
        ));
    }
    if k == l {
        return Ok(0.0);
    }

    // Wᵀ W = U Σ² Uᵀ  (eigenvectors `u[a][state]`, ascending eigenvalues).
    let gram = mbar_weight_gram(u_kn, n_k, f)?;
    let (gram_eigenvalues, u) = symmetric_eigen_n(&gram);
    let sigma: Vec<f64> = gram_eigenvalues
        .iter()
        .map(|&v| v.max(0.0).sqrt())
        .collect();

    // Q = I − Σ Uᵀ N U Σ
    let mut q = vec![vec![0.0_f64; n_states]; n_states];
    for a in 0..n_states {
        for b in a..n_states {
            let weighted_overlap: f64 = (0..n_states)
                .map(|state| n_k[state] as f64 * u[a][state] * u[b][state])
                .sum();
            let value = if a == b { 1.0 } else { 0.0 } - sigma[a] * sigma[b] * weighted_overlap;
            q[a][b] = value;
            q[b][a] = value;
        }
    }
    let (q_eigenvalues, q_vectors) = symmetric_eigen_n(&q);

    // Var(f_k − f_l) = yᵀ Q⁺ y   with   y_a = σ_a (U_ka − U_la).
    let y: Vec<f64> = (0..n_states)
        .map(|a| sigma[a] * (u[a][k] - u[a][l]))
        .collect();
    let mut variance = 0.0_f64;
    for (mu, vector) in q_eigenvalues.iter().zip(q_vectors.iter()) {
        if mu.abs() < MBAR_PSEUDO_INVERSE_CUTOFF {
            continue;
        }
        let projection: f64 = vector.iter().zip(y.iter()).map(|(v, yv)| v * yv).sum();
        variance += projection * projection / mu;
    }
    Ok(variance.max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::free_energy::bar::{bar_uncertainty, fep_bar};
    use crate::free_energy::functions::mbar_free_energies;
    use crate::free_energy::test_support::{mean_and_std, standard_normal};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    /// One-dimensional harmonic states with reduced potentials
    /// `u_k(x) = spring_k (x - center_k)² / 2` (β = 1).
    struct HarmonicStates {
        spring: Vec<f64>,
        center: Vec<f64>,
    }

    impl HarmonicStates {
        /// Exact reduced free energies relative to state 0:
        /// `f_k - f_0 = ln(sqrt(spring_k / spring_0))`.
        fn exact_f(&self) -> Vec<f64> {
            self.spring
                .iter()
                .map(|&k| 0.5 * (k / self.spring[0]).ln())
                .collect()
        }

        /// Draw `per_state` exact Boltzmann samples from every state and return
        /// the concatenated `u_kn` matrix together with `n_k`.
        fn sample(&self, rng: &mut StdRng, per_state: usize) -> (Vec<Vec<f64>>, Vec<usize>) {
            let n_states = self.spring.len();
            let mut x = Vec::with_capacity(n_states * per_state);
            for s in 0..n_states {
                let sigma = 1.0 / self.spring[s].sqrt();
                for _ in 0..per_state {
                    x.push(self.center[s] + sigma * standard_normal(rng));
                }
            }
            let u_kn = (0..n_states)
                .map(|l| {
                    x.iter()
                        .map(|&xi| 0.5 * self.spring[l] * (xi - self.center[l]).powi(2))
                        .collect()
                })
                .collect();
            (u_kn, vec![per_state; n_states])
        }
    }

    /// Self-consistent MBAR solver with the denominators evaluated once per
    /// sweep (the library routine recomputes them for every state).
    fn solve_mbar(u_kn: &[Vec<f64>], n_k: &[usize]) -> Vec<f64> {
        let n_states = u_kn.len();
        let n_total: usize = n_k.iter().sum();
        let ln_n: Vec<f64> = n_k.iter().map(|&c| (c as f64).ln()).collect();
        let mut f = vec![0.0_f64; n_states];
        for _ in 0..50_000 {
            let log_denominator: Vec<f64> = (0..n_total)
                .map(|s| {
                    let terms: Vec<f64> =
                        (0..n_states).map(|l| ln_n[l] + f[l] - u_kn[l][s]).collect();
                    log_sum_exp(&terms)
                })
                .collect();
            let mut f_new: Vec<f64> = (0..n_states)
                .map(|i| {
                    let terms: Vec<f64> = (0..n_total)
                        .map(|s| -u_kn[i][s] - log_denominator[s])
                        .collect();
                    -log_sum_exp(&terms)
                })
                .collect();
            let shift = f_new[0];
            f_new.iter_mut().for_each(|v| *v -= shift);
            let delta = f
                .iter()
                .zip(f_new.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            f = f_new;
            if delta < 1.0e-12 {
                break;
            }
        }
        f
    }

    fn three_harmonic_states() -> HarmonicStates {
        HarmonicStates {
            spring: vec![1.0, 2.0, 4.0],
            center: vec![0.0, 0.5, 1.0],
        }
    }

    #[test]
    fn variance_vanishes_for_identical_states_of_the_pair() {
        let states = three_harmonic_states();
        let mut rng = StdRng::seed_from_u64(3);
        let (u_kn, n_k) = states.sample(&mut rng, 60);
        let f = solve_mbar(&u_kn, &n_k);
        for state in 0..3 {
            let v = mbar_df_uncertainty(&u_kn, &n_k, &f, state, state).unwrap();
            assert_eq!(v, 0.0, "state {state}");
        }
    }

    #[test]
    fn variance_is_positive_symmetric_gauge_invariant_and_obeys_the_triangle_inequality() {
        let states = three_harmonic_states();
        let mut rng = StdRng::seed_from_u64(4);
        let (u_kn, n_k) = states.sample(&mut rng, 150);
        let f = solve_mbar(&u_kn, &n_k);
        let var = |k, l| mbar_df_uncertainty(&u_kn, &n_k, &f, k, l).unwrap();
        for &(k, l) in &[(0, 1), (1, 2), (0, 2)] {
            assert!(var(k, l) > 1.0e-6, "({k}, {l}) = {}", var(k, l));
            assert!((var(k, l) - var(l, k)).abs() < 1.0e-12 * var(k, l));
        }
        // Θ is a covariance matrix, so standard deviations obey the triangle inequality.
        let (s01, s12, s02) = (var(0, 1).sqrt(), var(1, 2).sqrt(), var(0, 2).sqrt());
        assert!(s02 <= s01 + s12 + 1.0e-9);
        assert!(s01 <= s02 + s12 + 1.0e-9);
        assert!(s12 <= s01 + s02 + 1.0e-9);
        // Shifting every free energy by a constant changes nothing.
        let shifted: Vec<f64> = f.iter().map(|v| v + 7.25).collect();
        let v_shifted = mbar_df_uncertainty(&u_kn, &n_k, &shifted, 0, 2).unwrap();
        assert!((v_shifted - var(0, 2)).abs() < 1.0e-9 * var(0, 2));
    }

    #[test]
    fn k_by_k_reduction_matches_the_brute_force_n_by_n_pseudo_inverse() {
        // Θ = Wᵀ (I - W N Wᵀ)⁺ W with an explicit N×N pseudo-inverse.
        let states = three_harmonic_states();
        let mut rng = StdRng::seed_from_u64(21);
        let (u_kn, n_k) = states.sample(&mut rng, 5);
        let n_total: usize = n_k.iter().sum();
        let f = mbar_free_energies(&u_kn, &n_k, 1.0e-13, 200_000);
        assert!(f.iter().all(|v| v.is_finite()));

        let ln_n: Vec<f64> = n_k.iter().map(|&c| (c as f64).ln()).collect();
        let w: Vec<Vec<f64>> = (0..n_total)
            .map(|s| {
                let terms: Vec<f64> = (0..3).map(|l| ln_n[l] + f[l] - u_kn[l][s]).collect();
                let log_denominator = log_sum_exp(&terms);
                (0..3)
                    .map(|k| (f[k] - u_kn[k][s] - log_denominator).exp())
                    .collect()
            })
            .collect();
        // M = I - W N Wᵀ
        let mut m = vec![vec![0.0_f64; n_total]; n_total];
        for a in 0..n_total {
            for b in 0..n_total {
                let wnw: f64 = (0..3).map(|k| n_k[k] as f64 * w[a][k] * w[b][k]).sum();
                m[a][b] = if a == b { 1.0 } else { 0.0 } - wnw;
            }
        }
        let (mu, vectors) = symmetric_eigen_n(&m);
        // Θ = Wᵀ M⁺ W
        let mut theta = vec![vec![0.0_f64; 3]; 3];
        for (value, vector) in mu.iter().zip(vectors.iter()) {
            if value.abs() < 1.0e-9 {
                continue;
            }
            let projected: Vec<f64> = (0..3)
                .map(|k| (0..n_total).map(|s| w[s][k] * vector[s]).sum())
                .collect();
            for a in 0..3 {
                for b in 0..3 {
                    theta[a][b] += projected[a] * projected[b] / value;
                }
            }
        }
        for &(k, l) in &[(0usize, 1usize), (1, 2), (0, 2)] {
            let brute = theta[k][k] + theta[l][l] - 2.0 * theta[k][l];
            let reduced = mbar_df_uncertainty(&u_kn, &n_k, &f, k, l).unwrap();
            assert!(
                (brute - reduced).abs() < 1.0e-6 * brute.abs().max(1.0e-3),
                "({k}, {l}): brute force {brute} vs K x K form {reduced}"
            );
        }
    }

    #[test]
    fn predicted_variance_matches_the_scatter_over_replicate_datasets() {
        let states = three_harmonic_states();
        let exact = states.exact_f();
        let mut rng = StdRng::seed_from_u64(20_260_613);
        let n_replicates = 200;
        let mut estimates: Vec<Vec<f64>> =
            (0..3).map(|_| Vec::with_capacity(n_replicates)).collect();
        let mut predicted: Vec<Vec<f64>> =
            (0..3).map(|_| Vec::with_capacity(n_replicates)).collect();
        for _ in 0..n_replicates {
            let (u_kn, n_k) = states.sample(&mut rng, 100);
            let f = solve_mbar(&u_kn, &n_k);
            for (slot, &(k, l)) in [(1usize, 0usize), (2, 1), (2, 0)].iter().enumerate() {
                estimates[slot].push(f[k] - f[l]);
                predicted[slot].push(mbar_df_uncertainty(&u_kn, &n_k, &f, k, l).unwrap());
            }
        }
        for (slot, &(k, l)) in [(1usize, 0usize), (2, 1), (2, 0)].iter().enumerate() {
            let (mean, empirical_std) = mean_and_std(&estimates[slot]);
            let predicted_std = (predicted[slot].iter().sum::<f64>() / n_replicates as f64).sqrt();
            let truth = exact[k] - exact[l];
            assert!(
                (mean - truth).abs() < 4.0 * empirical_std / (n_replicates as f64).sqrt() + 1.0e-3,
                "({k}, {l}): replicate mean {mean} vs exact {truth}"
            );
            let relative = (predicted_std - empirical_std).abs() / empirical_std;
            assert!(
                relative < 0.25,
                "({k}, {l}): predicted std {predicted_std} vs replicate std {empirical_std}"
            );
        }
    }

    #[test]
    fn two_states_reproduce_bar_free_energy_and_bennett_variance() {
        let states = HarmonicStates {
            spring: vec![1.0, 2.0],
            center: vec![0.0, 0.5],
        };
        let mut rng = StdRng::seed_from_u64(77);
        let per_state = 4000;
        let (u_kn, n_k) = states.sample(&mut rng, per_state);
        let f = solve_mbar(&u_kn, &n_k);
        // Two-state MBAR *is* BAR.
        let forward: Vec<f64> = (0..per_state).map(|s| u_kn[1][s] - u_kn[0][s]).collect();
        let backward: Vec<f64> = (per_state..2 * per_state)
            .map(|s| u_kn[0][s] - u_kn[1][s])
            .collect();
        let df_bar = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        assert!((f[1] - df_bar).abs() < 1.0e-6, "{} vs {df_bar}", f[1]);
        // The two asymptotic variance estimators agree up to sampling noise.
        let mbar_variance = mbar_df_uncertainty(&u_kn, &n_k, &f, 1, 0).unwrap();
        let bennett_sigma = bar_uncertainty(&forward, &backward, df_bar, 1.0);
        let relative = (mbar_variance.sqrt() - bennett_sigma).abs() / bennett_sigma;
        assert!(
            relative < 0.1,
            "MBAR sigma {} vs Bennett sigma {bennett_sigma}",
            mbar_variance.sqrt()
        );
    }

    #[test]
    fn uncertainty_falls_as_one_over_the_sample_count() {
        let states = three_harmonic_states();
        let mut rng = StdRng::seed_from_u64(8);
        let mut mean_variance = |per_state: usize, n_runs: usize| -> f64 {
            let mut total = 0.0;
            for _ in 0..n_runs {
                let (u_kn, n_k) = states.sample(&mut rng, per_state);
                let f = solve_mbar(&u_kn, &n_k);
                total += mbar_df_uncertainty(&u_kn, &n_k, &f, 0, 2).unwrap();
            }
            total / n_runs as f64
        };
        let small = mean_variance(60, 40);
        let large = mean_variance(240, 40);
        let ratio = small / large;
        assert!(
            (3.0..5.5).contains(&ratio),
            "4x more samples should cut the variance ~4x, ratio {ratio}"
        );
    }

    #[test]
    fn states_with_identical_or_constant_offset_potentials_have_no_uncertainty() {
        let mut rng = StdRng::seed_from_u64(5);
        let x: Vec<f64> = (0..80).map(|_| standard_normal(&mut rng)).collect();
        let u0: Vec<f64> = x.iter().map(|v| 0.5 * v * v).collect();
        let n_k = vec![40, 40];
        // Same potential twice.
        let u_kn = vec![u0.clone(), u0.clone()];
        let v = mbar_df_uncertainty(&u_kn, &n_k, &[0.0, 0.0], 1, 0).unwrap();
        assert!(v.abs() < 1.0e-12, "identical states: {v}");
        // Same potential up to a constant: f_1 - f_0 = 2.5 exactly.
        let u1: Vec<f64> = u0.iter().map(|v| v + 2.5).collect();
        let u_kn = vec![u0, u1];
        let f = solve_mbar(&u_kn, &n_k);
        assert!((f[1] - 2.5).abs() < 1.0e-9, "{f:?}");
        let v = mbar_df_uncertainty(&u_kn, &n_k, &f, 1, 0).unwrap();
        assert!(v.abs() < 1.0e-12, "offset states: {v}");
    }

    #[test]
    fn unsampled_states_are_supported() {
        // State 2 has no samples of its own but is still a valid target state.
        let states = three_harmonic_states();
        let mut rng = StdRng::seed_from_u64(31);
        let (mut u_kn, _) = states.sample(&mut rng, 120);
        // Keep only the samples of states 0 and 1.
        for row in u_kn.iter_mut() {
            row.truncate(240);
        }
        let n_k = vec![120, 120, 0];
        let mut ln_n = [f64::NEG_INFINITY; 3];
        ln_n[0] = 120.0_f64.ln();
        ln_n[1] = 120.0_f64.ln();
        // Self-consistent iteration restricted to the sampled states, state 2 follows.
        let mut f = vec![0.0_f64; 3];
        for _ in 0..20_000 {
            let log_denominator: Vec<f64> = (0..240)
                .map(|s| {
                    let terms: Vec<f64> = (0..3).map(|l| ln_n[l] + f[l] - u_kn[l][s]).collect();
                    log_sum_exp(&terms)
                })
                .collect();
            let mut f_new: Vec<f64> = (0..3)
                .map(|i| {
                    let terms: Vec<f64> =
                        (0..240).map(|s| -u_kn[i][s] - log_denominator[s]).collect();
                    -log_sum_exp(&terms)
                })
                .collect();
            let shift = f_new[0];
            f_new.iter_mut().for_each(|v| *v -= shift);
            let delta = (0..3).map(|i| (f[i] - f_new[i]).abs()).fold(0.0, f64::max);
            f = f_new;
            if delta < 1.0e-12 {
                break;
            }
        }
        let v = mbar_df_uncertainty(&u_kn, &n_k, &f, 2, 0).unwrap();
        assert!(v.is_finite() && v > 0.0, "{v}");
    }

    #[test]
    fn invalid_input_is_an_error() {
        let u_kn = vec![vec![0.1, 0.2, 0.3, 0.4], vec![0.4, 0.3, 0.2, 0.1]];
        let n_k = vec![2, 2];
        let f = vec![0.0, 0.1];
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &f, 1, 0).is_ok());
        assert!(mbar_df_uncertainty(&[], &[], &[], 0, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &[2], &f, 1, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &[0.0], 1, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &f, 2, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &f, 0, 9).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &[0, 0], &f, 1, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &[3, 3], &f, 1, 0).is_err());
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &[0.0, f64::INFINITY], 1, 0).is_err());
        let nan = vec![vec![0.1, f64::NAN, 0.3, 0.4], u_kn[1].clone()];
        assert!(mbar_df_uncertainty(&nan, &n_k, &f, 1, 0).is_err());
        // Index checks come first: a bad index is an error even for k == l.
        assert!(mbar_df_uncertainty(&u_kn, &n_k, &f, 5, 5).is_err());
    }
}
