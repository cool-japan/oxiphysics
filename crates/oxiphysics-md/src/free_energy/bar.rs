// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bennett acceptance ratio (BAR) free-energy estimator and its asymptotic
//! uncertainty.
//!
//! # Conventions
//!
//! Two thermodynamic states 0 and 1 have potentials `U₀` and `U₁`, and
//! `ΔU = U₁ − U₀`.  The estimator consumes
//!
//! * `du_forward`: `ΔU = U₁ − U₀` evaluated on `n₀` configurations sampled from
//!   state 0 (the forward work `W_F`),
//! * `du_backward`: `U₀ − U₁` evaluated on `n₁` configurations sampled from
//!   state 1 (the reverse work `W_R`),
//!
//! and returns `ΔF = F₁ − F₀`.
//!
//! # Equation
//!
//! With the Fermi function `f(x) = 1 / (1 + eˣ)` and the shift
//! `C = ΔF − kT ln(n₀ / n₁)`, Bennett's (1976) self-consistent equation is
//!
//! ```text
//! Σ_{i ∈ 0} f((W_F,i − C) / kT)  =  Σ_{j ∈ 1} f((W_R,j + C) / kT)
//! ```
//!
//! Writing `num` and `den` for the right- and left-hand side, the residual
//! `g(ΔF) = ln num − ln den` is strictly decreasing in `ΔF`, which makes the
//! root unique and lets the solver bracket it.  All sums are carried out in log
//! space, so very large energy differences (hundreds of `kT`) neither overflow
//! nor underflow.

use super::functions::fep_zwanzig;

/// Upper bound on the number of root-finder iterations.
const BAR_MAX_ITERATIONS: usize = 400;

/// Upper bound on the number of bracket-doubling steps.
const BAR_MAX_BRACKET_EXPANSIONS: usize = 200;

/// Natural logarithm of the Fermi function `f(x) = 1 / (1 + eˣ)`.
///
/// Evaluated as `-softplus(x)`, so it neither overflows for large positive `x`
/// nor loses precision for large negative `x`.  `ln f(+∞) = −∞` and
/// `ln f(−∞) = 0`; a NaN argument yields NaN.
#[inline]
pub(crate) fn log_fermi(x: f64) -> f64 {
    -(x.max(0.0) + (-x.abs()).exp().ln_1p())
}

/// Log-space summary of `Σ_i f(a_i)` with `a_i = beta · du_i + offset`.
#[derive(Debug, Clone, Copy)]
struct FermiSum {
    /// `ln Σ_i f(a_i)`; `−∞` when every term vanishes, NaN for NaN input.
    log_sum: f64,
    /// `Σ_i f_i (1 − f_i) / Σ_i f_i`: the mean of `1 − f` weighted by `f`.
    mean_complement: f64,
}

/// Evaluate [`FermiSum`] for the samples `du`.
fn fermi_sum(du: &[f64], beta: f64, offset: f64) -> FermiSum {
    if du.iter().any(|d| d.is_nan()) {
        return FermiSum {
            log_sum: f64::NAN,
            mean_complement: f64::NAN,
        };
    }
    let max_log = du
        .iter()
        .map(|&d| log_fermi(beta * d + offset))
        .fold(f64::NEG_INFINITY, f64::max);
    if max_log == f64::NEG_INFINITY {
        return FermiSum {
            log_sum: f64::NEG_INFINITY,
            mean_complement: 0.0,
        };
    }
    let mut total = 0.0_f64;
    let mut complement = 0.0_f64;
    for &d in du {
        let a = beta * d + offset;
        let weight = (log_fermi(a) - max_log).exp();
        total += weight;
        // 1 - f(a) = f(-a)
        complement += weight * log_fermi(-a).exp();
    }
    FermiSum {
        log_sum: max_log + total.ln(),
        mean_complement: complement / total,
    }
}

/// `ln ⟨f(beta · du + offset)⟩` over the samples `du` (their arithmetic mean
/// inside the logarithm).
///
/// Returns `−∞` for an empty slice or when every term vanishes.
pub(crate) fn log_mean_fermi(du: &[f64], beta: f64, offset: f64) -> f64 {
    if du.is_empty() {
        return f64::NEG_INFINITY;
    }
    fermi_sum(du, beta, offset).log_sum - (du.len() as f64).ln()
}

/// Residual of the BAR equation at a trial free energy.
#[derive(Debug, Clone, Copy)]
struct BarResidual {
    /// `ln num − ln den`; positive when the trial `ΔF` is too small.
    value: f64,
    /// `−∂value / ∂(β ΔF)`, always in `(0, 2)` for finite data.
    curvature: f64,
}

fn bar_residual(
    du_forward: &[f64],
    du_backward: &[f64],
    beta: f64,
    ln_ratio: f64,
    df: f64,
) -> BarResidual {
    let c = beta * df - ln_ratio;
    let num = fermi_sum(du_backward, beta, c);
    let den = fermi_sum(du_forward, beta, -c);
    BarResidual {
        value: num.log_sum - den.log_sum,
        curvature: num.mean_complement + den.mean_complement,
    }
}

fn slice_mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Solve the two-sided BAR equation for `ΔF` (both sample sets non-empty and
/// `kt` positive and finite).
fn solve_bar(du_forward: &[f64], du_backward: &[f64], kt: f64, tol: f64) -> f64 {
    let beta = 1.0 / kt;
    let ln_ratio = (du_forward.len() as f64 / du_backward.len() as f64).ln();
    let residual = |df: f64| bar_residual(du_forward, du_backward, beta, ln_ratio, df);

    // The Jensen bounds -<W_R> <= ΔF <= <W_F> bracket the root in expectation;
    // the loops below widen the interval until the sample equation really does
    // change sign inside it.
    let mean_forward = slice_mean(du_forward);
    let mean_backward = slice_mean(du_backward);
    let mut lo = (-mean_backward).min(mean_forward);
    let mut hi = (-mean_backward).max(mean_forward);
    if !(lo.is_finite() && hi.is_finite()) {
        lo = -kt;
        hi = kt;
    }
    let mut width = (hi - lo).max(kt);
    let mut bracketed_lo = false;
    for _ in 0..BAR_MAX_BRACKET_EXPANSIONS {
        let g = residual(lo).value;
        if g.is_nan() {
            return f64::NAN;
        }
        if g >= 0.0 {
            bracketed_lo = true;
            break;
        }
        lo -= width;
        width *= 2.0;
    }
    let mut width = (hi - lo).max(kt);
    let mut bracketed_hi = false;
    for _ in 0..BAR_MAX_BRACKET_EXPANSIONS {
        let g = residual(hi).value;
        if g.is_nan() {
            return f64::NAN;
        }
        if g <= 0.0 {
            bracketed_hi = true;
            break;
        }
        hi += width;
        width *= 2.0;
    }
    if !(bracketed_lo && bracketed_hi) {
        return f64::NAN;
    }

    // Safeguarded Newton iteration: g is strictly decreasing, so every
    // evaluation tightens one end of the bracket, and a Newton step that leaves
    // the bracket is replaced by bisection.
    let mut x = 0.5 * (lo + hi);
    for _ in 0..BAR_MAX_ITERATIONS {
        let r = residual(x);
        if r.value.is_nan() {
            return f64::NAN;
        }
        if r.value == 0.0 {
            return x;
        }
        if r.value > 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let newton = if r.curvature > 0.0 {
            x + kt * r.value / r.curvature
        } else {
            f64::NAN
        };
        let next = if newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };
        if next == x || (next - x).abs() <= tol || hi - lo <= tol {
            return next;
        }
        x = next;
    }
    x
}

/// Bennett Acceptance Ratio (BAR) free energy estimate `ΔF = F₁ − F₀`.
///
/// Solves Bennett's self-consistent equation (see the module documentation)
/// with a bracketed Newton iteration until the step falls below `tol`.  The
/// Zwanzig estimate is used as a fallback when samples from only one state are
/// available.
///
/// # Arguments
/// * `du_forward`  – ΔU = U₁ − U₀ samples from state 0.
/// * `du_backward` – ΔU = U₀ − U₁ samples from state 1.
/// * `kt`          – thermal energy k_B T.
/// * `tol`         – convergence tolerance (absolute, in units of energy).
///
/// Returns `0` when both sample sets are empty and NaN when `kt` is not
/// positive and finite, when the input contains NaN, or when the two sample
/// sets share no phase-space overlap at all (every Fermi weight vanishes).
pub fn fep_bar(du_forward: &[f64], du_backward: &[f64], kt: f64, tol: f64) -> f64 {
    if du_forward.is_empty() && du_backward.is_empty() {
        return 0.0;
    }
    if du_backward.is_empty() {
        return fep_zwanzig(du_forward, kt);
    }
    if du_forward.is_empty() {
        return -fep_zwanzig(du_backward, kt);
    }
    if !(kt > 0.0 && kt.is_finite()) {
        return f64::NAN;
    }
    solve_bar(du_forward, du_backward, kt, tol)
}

/// `Var(f) / ⟨f⟩²` of the Fermi weights `f(beta · du + offset)`.
///
/// The ratio is invariant under a common rescaling of `f`, so the weights are
/// normalised to their maximum in log space and cannot underflow.
fn fermi_relative_variance(du: &[f64], beta: f64, offset: f64) -> f64 {
    let log_f: Vec<f64> = du.iter().map(|&d| log_fermi(beta * d + offset)).collect();
    let max_log = log_f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max_log.is_finite() || log_f.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let n = log_f.len() as f64;
    let scaled: Vec<f64> = log_f.iter().map(|&v| (v - max_log).exp()).collect();
    let mean = scaled.iter().sum::<f64>() / n;
    let var = scaled.iter().map(|&u| (u - mean).powi(2)).sum::<f64>() / n;
    var / (mean * mean)
}

/// Asymptotic standard error of the BAR estimate (Bennett 1976, Shirts *et
/// al.* 2003, Eq. 10).
///
/// ```text
/// σ²(ΔF) = kT² [ (⟨f²⟩₀ / ⟨f⟩₀² − 1) / n₀  +  (⟨f²⟩₁ / ⟨f⟩₁² − 1) / n₁ ]
/// ```
///
/// where `f₀ = f((ΔU − C) / kT)` is evaluated on the forward samples,
/// `f₁ = f((W_R + C) / kT)` on the backward samples and
/// `C = ΔF − kT ln(n₀ / n₁)`.  The expression is the variance of the logarithm
/// of the ratio of the two Fermi averages, and it is exact to leading order at
/// the self-consistent `ΔF`, because the estimator is stationary with respect
/// to `C` there.  It vanishes when the Fermi weights are constant (for example
/// a pure constant shift of the potential).
///
/// # Arguments
/// * `du_forward`  – ΔU = U₁ − U₀ samples from state 0.
/// * `du_backward` – ΔU = U₀ − U₁ samples from state 1.
/// * `df`          – the BAR free energy, normally the result of [`fep_bar`].
/// * `kt`          – thermal energy k_B T.
///
/// Returns `0` when either sample set has fewer than two entries or `kt` is
/// zero, and NaN when the weights cannot be evaluated.
pub fn bar_uncertainty(du_forward: &[f64], du_backward: &[f64], df: f64, kt: f64) -> f64 {
    let n0 = du_forward.len();
    let n1 = du_backward.len();
    if n0 < 2 || n1 < 2 || kt == 0.0 {
        return 0.0;
    }
    let beta = 1.0 / kt;
    let c = beta * df - (n0 as f64 / n1 as f64).ln();
    let var0 = fermi_relative_variance(du_forward, beta, -c) / n0 as f64;
    let var1 = fermi_relative_variance(du_backward, beta, c) / n1 as f64;
    kt.abs() * (var0 + var1).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::free_energy::test_support::{mean_and_std, standard_normal};
    use crate::free_energy::{LambdaWindow, alchemical_bar_profile, total_fep_bar};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    /// Forward and reverse work samples of a Gaussian work distribution.
    ///
    /// `W_F ~ N(mu, sigma²)` and `-W_R ~ N(mu - sigma²/kt, sigma²)`, which by
    /// Crooks' theorem corresponds to `ΔF = mu - sigma²/(2 kt)`.
    fn gaussian_work(
        rng: &mut StdRng,
        mu: f64,
        sigma: f64,
        kt: f64,
        n_forward: usize,
        n_backward: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let forward: Vec<f64> = (0..n_forward)
            .map(|_| mu + sigma * standard_normal(rng))
            .collect();
        let mu_reverse = mu - sigma * sigma / kt;
        let backward: Vec<f64> = (0..n_backward)
            .map(|_| -(mu_reverse + sigma * standard_normal(rng)))
            .collect();
        (forward, backward)
    }

    #[test]
    fn log_fermi_matches_direct_evaluation_and_saturates() {
        for &x in &[-30.0, -2.5, -0.1, 0.0, 0.1, 2.5, 30.0] {
            let direct = (1.0 / (1.0 + f64::exp(x))).ln();
            assert!(
                (log_fermi(x) - direct).abs() < 1.0e-12,
                "x = {x}: {} vs {direct}",
                log_fermi(x)
            );
        }
        assert_eq!(log_fermi(f64::INFINITY), f64::NEG_INFINITY);
        assert_eq!(log_fermi(f64::NEG_INFINITY), 0.0);
        // The naive 1/(1+exp(x)) would be exactly zero here.
        assert!((log_fermi(1000.0) + 1000.0).abs() < 1.0e-12);
        assert!(log_fermi(f64::NAN).is_nan());
    }

    #[test]
    fn bar_recovers_constant_shift_with_equal_sample_counts() {
        // Pure shift U1 = U0 + 1: the exact answer is 1.  The previous
        // implementation converged to 0.
        let forward = vec![1.0_f64; 100];
        let backward = vec![-1.0_f64; 100];
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        assert!((df - 1.0).abs() < 1.0e-8, "expected 1, got {df}");
    }

    #[test]
    fn bar_recovers_constant_shift_with_unbalanced_sample_counts() {
        let forward = vec![1.0_f64; 100];
        let backward = vec![-1.0_f64; 50];
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        assert!((df - 1.0).abs() < 1.0e-8, "expected 1, got {df}");
        let df_swapped = fep_bar(&backward, &forward, 1.0, 1.0e-12);
        assert!(
            (df_swapped + 1.0).abs() < 1.0e-8,
            "expected -1, got {df_swapped}"
        );
    }

    #[test]
    fn bar_constant_shift_is_independent_of_kt_and_counts() {
        for &(kt, shift, n0, n1) in &[(2.5, 3.0, 30_usize, 70_usize), (0.4, -1.3, 90, 10)] {
            let forward = vec![shift; n0];
            let backward = vec![-shift; n1];
            let df = fep_bar(&forward, &backward, kt, 1.0e-12);
            assert!(
                (df - shift).abs() < 1.0e-8,
                "kt = {kt}, shift = {shift}: got {df}"
            );
        }
    }

    #[test]
    fn bar_handles_hundreds_of_kt_without_overflow() {
        let forward = vec![1000.0_f64; 10];
        let backward = vec![-1000.0_f64; 10];
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-10);
        assert!(
            (df - 1000.0).abs() < 1.0e-6,
            "expected 1000 for a 1000 kT shift, got {df}"
        );
    }

    #[test]
    fn bar_is_antisymmetric_under_state_exchange() {
        let mut rng = StdRng::seed_from_u64(11);
        let (forward, backward) = gaussian_work(&mut rng, 2.0, 1.5, 1.0, 300, 200);
        let a = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        let b = fep_bar(&backward, &forward, 1.0, 1.0e-12);
        assert!((a + b).abs() < 1.0e-9, "{a} vs {b}");
    }

    #[test]
    fn bar_lies_between_the_one_sided_jensen_bounds() {
        let mut rng = StdRng::seed_from_u64(5);
        let (forward, backward) = gaussian_work(&mut rng, 3.0, 2.0, 1.0, 2000, 2000);
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-10);
        let upper = forward.iter().sum::<f64>() / forward.len() as f64;
        let lower = -backward.iter().sum::<f64>() / backward.len() as f64;
        assert!(lower < df && df < upper, "{lower} < {df} < {upper}");
    }

    #[test]
    fn bar_recovers_gaussian_crooks_free_energy() {
        let mut rng = StdRng::seed_from_u64(2024);
        let (mu, sigma, kt) = (3.0, 2.0, 1.0);
        let (forward, backward) = gaussian_work(&mut rng, mu, sigma, kt, 20_000, 20_000);
        let df = fep_bar(&forward, &backward, kt, 1.0e-10);
        let exact = mu - sigma * sigma / (2.0 * kt);
        assert!(
            (df - exact).abs() < 0.05,
            "BAR = {df}, Crooks value = {exact}"
        );
    }

    #[test]
    fn bar_recovers_gaussian_crooks_free_energy_for_unequal_counts_and_kt() {
        let mut rng = StdRng::seed_from_u64(77);
        let (mu, sigma, kt) = (6.0, 2.5, 2.0);
        let (forward, backward) = gaussian_work(&mut rng, mu, sigma, kt, 30_000, 10_000);
        let df = fep_bar(&forward, &backward, kt, 1.0e-10);
        let exact = mu - sigma * sigma / (2.0 * kt);
        assert!(
            (df - exact).abs() < 0.08,
            "BAR = {df}, Crooks value = {exact}"
        );
    }

    #[test]
    fn bar_one_sided_and_empty_fallbacks() {
        assert_eq!(fep_bar(&[], &[], 1.0, 1.0e-10), 0.0);
        let shifted = vec![0.75_f64; 40];
        assert!((fep_bar(&shifted, &[], 1.0, 1.0e-10) - 0.75).abs() < 1.0e-10);
        assert!((fep_bar(&[], &shifted, 1.0, 1.0e-10) + 0.75).abs() < 1.0e-10);
    }

    #[test]
    fn bar_rejects_non_positive_temperature_and_nan_input() {
        let a = [0.5_f64; 4];
        assert!(fep_bar(&a, &a, 0.0, 1.0e-10).is_nan());
        assert!(fep_bar(&a, &a, -1.0, 1.0e-10).is_nan());
        assert!(fep_bar(&[f64::NAN, 1.0], &[0.0, 1.0], 1.0, 1.0e-10).is_nan());
    }

    #[test]
    fn bar_with_non_overlapping_infinite_energies_is_nan() {
        // Every forward sample has an infinite energy difference and so does
        // every backward sample: no information about ΔF at all.
        let inf = [f64::INFINITY; 3];
        assert!(fep_bar(&inf, &inf, 1.0, 1.0e-10).is_nan());
    }

    #[test]
    fn bar_ignores_forward_samples_with_infinite_energy() {
        // A forward sample that would put the system into a forbidden region of
        // state 1 (ΔU = +∞) simply carries zero Fermi weight.
        let forward = vec![1.0, 1.0, 1.0, f64::INFINITY];
        let backward = vec![-1.0, -1.0, -1.0];
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        assert!(df.is_finite());
        assert!((df - 1.0).abs() < 1.0, "unexpected estimate {df}");
    }

    #[test]
    fn bar_residual_vanishes_at_the_returned_free_energy() {
        let mut rng = StdRng::seed_from_u64(9);
        let (forward, backward) = gaussian_work(&mut rng, 1.0, 1.0, 1.0, 400, 250);
        let kt = 1.0;
        let df = fep_bar(&forward, &backward, kt, 1.0e-13);
        let ln_ratio = (400.0_f64 / 250.0).ln();
        let r = bar_residual(&forward, &backward, 1.0 / kt, ln_ratio, df);
        assert!(r.value.abs() < 1.0e-9, "residual {}", r.value);
        assert!(r.curvature > 0.0 && r.curvature < 2.0);
    }

    #[test]
    fn bar_uncertainty_vanishes_for_constant_weights() {
        let forward = vec![0.5_f64; 50];
        let backward = vec![-0.5_f64; 50];
        let df = fep_bar(&forward, &backward, 1.0, 1.0e-12);
        let sigma = bar_uncertainty(&forward, &backward, df, 1.0);
        assert!(sigma.abs() < 1.0e-6, "expected 0, got {sigma}");
    }

    #[test]
    fn bar_uncertainty_needs_two_samples_per_state() {
        assert_eq!(bar_uncertainty(&[1.0], &[0.0, 1.0], 0.0, 1.0), 0.0);
        assert_eq!(bar_uncertainty(&[1.0, 2.0], &[0.0], 0.0, 1.0), 0.0);
        assert_eq!(bar_uncertainty(&[1.0, 2.0], &[0.0, 1.0], 0.0, 0.0), 0.0);
    }

    #[test]
    fn bar_uncertainty_scales_with_kt_and_inverse_sqrt_of_sample_count() {
        let mut rng = StdRng::seed_from_u64(31);
        let (fwd_small, bwd_small) = gaussian_work(&mut rng, 3.0, 2.0, 1.0, 400, 400);
        let (fwd_large, bwd_large) = gaussian_work(&mut rng, 3.0, 2.0, 1.0, 6400, 6400);
        let s_small = bar_uncertainty(
            &fwd_small,
            &bwd_small,
            fep_bar(&fwd_small, &bwd_small, 1.0, 1.0e-12),
            1.0,
        );
        let s_large = bar_uncertainty(
            &fwd_large,
            &bwd_large,
            fep_bar(&fwd_large, &bwd_large, 1.0, 1.0e-12),
            1.0,
        );
        let ratio = s_small / s_large;
        assert!(
            (3.0..5.5).contains(&ratio),
            "16x more samples should shrink sigma by ~4, got ratio {ratio}"
        );
    }

    #[test]
    fn bar_uncertainty_matches_replicate_scatter() {
        let (mu, sigma, kt) = (3.0, 2.0, 1.0);
        let mut rng = StdRng::seed_from_u64(1234);
        let n_replicates = 200;
        let mut estimates = Vec::with_capacity(n_replicates);
        let mut variances = Vec::with_capacity(n_replicates);
        for _ in 0..n_replicates {
            let (forward, backward) = gaussian_work(&mut rng, mu, sigma, kt, 400, 400);
            let df = fep_bar(&forward, &backward, kt, 1.0e-10);
            let predicted = bar_uncertainty(&forward, &backward, df, kt);
            estimates.push(df);
            variances.push(predicted * predicted);
        }
        let (mean_df, empirical_std) = mean_and_std(&estimates);
        let predicted_std = (variances.iter().sum::<f64>() / n_replicates as f64).sqrt();
        let exact = mu - sigma * sigma / (2.0 * kt);
        assert!(
            (mean_df - exact).abs() < 0.05,
            "replicate mean {mean_df} vs Crooks value {exact}"
        );
        let relative = (predicted_std - empirical_std).abs() / empirical_std;
        assert!(
            relative < 0.25,
            "predicted sigma {predicted_std} vs replicate scatter {empirical_std}"
        );
    }

    #[test]
    fn window_chain_sums_constant_shifts() {
        let n_windows = 4;
        let shift = 0.5;
        let mut windows: Vec<LambdaWindow> = (0..n_windows)
            .map(|i| LambdaWindow::new(i as f64 / (n_windows - 1) as f64))
            .collect();
        for i in 0..n_windows - 1 {
            for _ in 0..40 {
                windows[i].add_forward(shift);
                windows[i + 1].add_backward(-shift);
            }
        }
        let first = windows[0].bar_estimate(&windows[1], 1.0);
        assert!((first - shift).abs() < 1.0e-8, "window estimate {first}");
        let total = total_fep_bar(&windows, 1.0);
        let expected = shift * (n_windows - 1) as f64;
        assert!(
            (total - expected).abs() < 1.0e-8,
            "expected {expected}, got {total}"
        );
    }

    #[test]
    fn alchemical_profile_accumulates_constant_shifts() {
        let forwards = vec![vec![0.5_f64; 20]; 3];
        let backwards = vec![vec![-0.5_f64; 20]; 3];
        let profile = alchemical_bar_profile(&forwards, &backwards, 1.0, 1.0e-12);
        assert_eq!(profile.len(), 4);
        for (i, v) in profile.iter().enumerate() {
            assert!((v - 0.5 * i as f64).abs() < 1.0e-8, "profile[{i}] = {v}");
        }
    }
}
