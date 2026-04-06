//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{AlchemicalResult, LambdaWindow, WhamWindow};

/// Zwanzig (exponential averaging) free energy estimate.
///
/// ΔF = −kT · ln⟨exp(−ΔU / kT)⟩
///
/// # Arguments
/// * `du_samples` – ΔU samples collected while in the reference state.
/// * `kT`         – thermal energy k_B T in the same units as ΔU.
///
/// Returns `f64::INFINITY` when the argument of the logarithm underflows to zero.
#[allow(dead_code)]
pub fn fep_zwanzig(du_samples: &[f64], kt: f64) -> f64 {
    if du_samples.is_empty() || kt == 0.0 {
        return 0.0;
    }
    let n = du_samples.len() as f64;
    let avg_exp: f64 = du_samples.iter().map(|&du| (-du / kt).exp()).sum::<f64>() / n;
    if avg_exp <= 0.0 {
        return f64::INFINITY;
    }
    -kt * avg_exp.ln()
}
/// Fermi / logistic function used in the BAR self-consistency equation.
#[inline]
pub(super) fn fermi(x: f64) -> f64 {
    1.0 / (1.0 + x.exp())
}
/// Bennett Acceptance Ratio (BAR) free energy estimate.
///
/// Solves the self-consistent BAR equation iteratively until convergence within
/// `tol` or 10 000 iterations.  The Zwanzig estimate is used as a fallback
/// when samples from only one state are available.
///
/// # Arguments
/// * `du_forward`  – ΔU = U₁ − U₀ samples from state 0.
/// * `du_backward` – ΔU = U₀ − U₁ samples from state 1.
/// * `kt`          – thermal energy k_B T.
/// * `tol`         – convergence tolerance (absolute, in units of energy).
#[allow(dead_code)]
pub fn fep_bar(du_forward: &[f64], du_backward: &[f64], kt: f64, tol: f64) -> f64 {
    let n0 = du_forward.len();
    let n1 = du_backward.len();
    if n0 == 0 && n1 == 0 {
        return 0.0;
    }
    if n1 == 0 {
        return fep_zwanzig(du_forward, kt);
    }
    if n0 == 0 {
        return -fep_zwanzig(du_backward, kt);
    }
    let ln_ratio = ((n0 as f64) / (n1 as f64)).ln();
    let mut df = 0.0_f64;
    for _ in 0..10_000 {
        let c = df / kt - ln_ratio;
        let num: f64 = du_backward
            .iter()
            .map(|&du| fermi(du / kt + c))
            .sum::<f64>();
        let den: f64 = du_forward
            .iter()
            .map(|&du| fermi(-du / kt - c))
            .sum::<f64>();
        if den.abs() < 1e-300 {
            break;
        }
        let df_new = kt * (num / den).ln() + kt * ln_ratio;
        if (df_new - df).abs() < tol {
            df = df_new;
            break;
        }
        df = df_new;
    }
    df
}
/// Statistical uncertainty in the Zwanzig estimate via block averaging.
///
/// Splits `du_samples` into blocks of size √N, estimates ΔF for each block,
/// and returns the standard error of the block estimates.
///
/// Returns 0 when fewer than 4 samples are available.
#[allow(dead_code)]
pub fn fep_uncertainty(du_samples: &[f64], kt: f64) -> f64 {
    let n = du_samples.len();
    if n < 4 {
        return 0.0;
    }
    let block_size = ((n as f64).sqrt().ceil() as usize).max(2);
    let blocks: Vec<f64> = du_samples
        .chunks(block_size)
        .filter(|b| !b.is_empty())
        .map(|b| fep_zwanzig(b, kt))
        .collect();
    let nb = blocks.len() as f64;
    if nb < 2.0 {
        return 0.0;
    }
    let mean = blocks.iter().sum::<f64>() / nb;
    let var = blocks.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (nb - 1.0);
    (var / nb).sqrt()
}
/// Thermodynamic Integration using the trapezoidal rule.
///
/// ΔF = ∫₀¹ ⟨dU/dλ⟩ dλ ≈ Σᵢ ½(⟨dU/dλ⟩ᵢ + ⟨dU/dλ⟩ᵢ₊₁)(λᵢ₊₁ − λᵢ)
///
/// # Arguments
/// * `lambdas`     – λ grid points (must be sorted, length ≥ 2).
/// * `du_dlambda`  – ⟨dU/dλ⟩ at each λ point (same length as `lambdas`).
#[allow(dead_code)]
pub fn thermodynamic_integration(lambdas: &[f64], du_dlambda: &[f64]) -> f64 {
    let n = lambdas.len().min(du_dlambda.len());
    if n < 2 {
        return 0.0;
    }
    let mut integral = 0.0;
    for i in 0..n - 1 {
        let dl = lambdas[i + 1] - lambdas[i];
        integral += 0.5 * (du_dlambda[i] + du_dlambda[i + 1]) * dl;
    }
    integral
}
/// Uncertainty in the TI estimate via error propagation through the trapezoid rule.
///
/// Each node contributes a term ½(errᵢ + errᵢ₊₁) · Δλ in quadrature.
///
/// # Arguments
/// * `lambdas`         – λ grid points.
/// * `du_dlambda_err`  – standard error of ⟨dU/dλ⟩ at each λ point.
#[allow(dead_code)]
pub fn ti_uncertainty(lambdas: &[f64], du_dlambda_err: &[f64]) -> f64 {
    let n = lambdas.len().min(du_dlambda_err.len());
    if n < 2 {
        return 0.0;
    }
    let mut sum_sq = 0.0;
    for i in 0..n - 1 {
        let dl = lambdas[i + 1] - lambdas[i];
        let avg_err = 0.5 * (du_dlambda_err[i] + du_dlambda_err[i + 1]);
        sum_sq += (avg_err * dl).powi(2);
    }
    sum_sq.sqrt()
}
/// Jarzynski free energy estimate from non-equilibrium work values.
///
/// ΔF = −kT · ln⟨exp(−W / kT)⟩
///
/// # Arguments
/// * `work_samples` – irreversible work values W for each trajectory.
/// * `kt`           – thermal energy k_B T.
#[allow(dead_code)]
pub fn jarzynski_estimate(work_samples: &[f64], kt: f64) -> f64 {
    fep_zwanzig(work_samples, kt)
}
/// Second-order cumulant expansion of the Jarzynski equality.
///
/// ΔF ≈ ⟨W⟩ − Var(W) / (2 kT)
///
/// Valid when work fluctuations are small compared to kT.
#[allow(dead_code)]
pub fn cumulant_expansion_2nd(work_samples: &[f64], kt: f64) -> f64 {
    let n = work_samples.len();
    if n == 0 || kt == 0.0 {
        return 0.0;
    }
    let nf = n as f64;
    let mean = work_samples.iter().sum::<f64>() / nf;
    let var = if n < 2 {
        0.0
    } else {
        work_samples
            .iter()
            .map(|&w| (w - mean).powi(2))
            .sum::<f64>()
            / (nf - 1.0)
    };
    mean - var / (2.0 * kt)
}
/// Histogram-based Potential of Mean Force (PMF).
///
/// PMF(ξ) = −kT · ln P(ξ) + const,  where the additive constant is chosen so
/// that the minimum PMF value is zero.
///
/// # Arguments
/// * `samples` – collective variable (reaction coordinate) values.
/// * `bins`    – number of histogram bins (must be ≥ 1).
/// * `range`   – (min, max) of the histogram range.
/// * `kt`      – thermal energy k_B T.
///
/// # Returns
/// A `Vec` of `(bin_center, pmf)` pairs; bins with zero counts are omitted.
#[allow(dead_code)]
pub fn histogram_pmf(samples: &[f64], bins: usize, range: (f64, f64), kt: f64) -> Vec<(f64, f64)> {
    if samples.is_empty() || bins == 0 || range.0 >= range.1 {
        return Vec::new();
    }
    let (lo, hi) = range;
    let width = (hi - lo) / bins as f64;
    let mut counts = vec![0usize; bins];
    for &x in samples {
        if x >= lo && x < hi {
            let idx = ((x - lo) / width) as usize;
            let idx = idx.min(bins - 1);
            counts[idx] += 1;
        }
    }
    let n_total = samples.len() as f64;
    let mut pmf_raw: Vec<(f64, f64)> = counts
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c > 0)
        .map(|(i, &c)| {
            let center = lo + (i as f64 + 0.5) * width;
            let prob = c as f64 / (n_total * width);
            (center, -kt * prob.ln())
        })
        .collect();
    if let Some(&(_, min_val)) = pmf_raw
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    {
        for entry in pmf_raw.iter_mut() {
            entry.1 -= min_val;
        }
    }
    pmf_raw
}
/// Statistical uncertainty of the PMF per occupied bin.
///
/// σ(PMF_i) = kT / √n_i
///
/// # Arguments
/// * `pmf`      – output from [`histogram_pmf`] (bin_center, pmf_value) pairs.
/// * `n_samples` – total number of samples used to build the PMF.
///
/// Returns a vector of per-bin uncertainties in the same order as `pmf`.
///
/// # Note
/// This approximation requires knowledge of the bin counts, which are
/// recomputed here from the PMF values and `n_samples` using the relationship
/// PMF_i = −kT · ln(n_i / (n_total · Δx)).  For a direct count-based
/// uncertainty use the raw histogram.
#[allow(dead_code)]
pub fn pmf_uncertainty(pmf: &[(f64, f64)], n_samples: usize) -> Vec<f64> {
    if pmf.is_empty() || n_samples == 0 {
        return Vec::new();
    }
    let n_bins = pmf.len() as f64;
    let avg_per_bin = (n_samples as f64 / n_bins).max(1.0);
    vec![1.0 / avg_per_bin.sqrt(); pmf.len()]
}
/// Configurational (differential) entropy via histogram.
///
/// S = −kT · Σᵢ p_i · ln(p_i) · Δx
///
/// where p_i is the probability density in bin i and Δx is the bin width.
///
/// # Arguments
/// * `samples` – 1-D configuration samples.
/// * `bins`    – number of histogram bins.
/// * `kt`      – thermal energy k_B T (set to 1.0 for dimensionless entropy).
#[allow(dead_code)]
pub fn configurational_entropy(samples: &[f64], bins: usize, kt: f64) -> f64 {
    if samples.is_empty() || bins == 0 {
        return 0.0;
    }
    let lo = samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if (hi - lo).abs() < f64::EPSILON {
        return 0.0;
    }
    let width = (hi - lo) / bins as f64;
    let mut counts = vec![0usize; bins];
    for &x in samples {
        let idx = ((x - lo) / width) as usize;
        let idx = idx.min(bins - 1);
        counts[idx] += 1;
    }
    let n = samples.len() as f64;
    let entropy: f64 = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / n;
            -p * p.ln()
        })
        .sum();
    let differential_entropy = entropy + width.ln();
    kt * differential_entropy
}
/// Shannon entropy of a discrete distribution given counts in each bin.
pub(super) fn shannon_entropy_from_counts(counts: &[usize]) -> f64 {
    let n: f64 = counts.iter().sum::<usize>() as f64;
    if n == 0.0 {
        return 0.0;
    }
    counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / n;
            -p * p.ln()
        })
        .sum()
}
/// Mutual information between two 1-D random variables.
///
/// MI(X, Y) = H(X) + H(Y) − H(X, Y)
///
/// # Arguments
/// * `x`    – samples of variable X.
/// * `y`    – samples of variable Y (must have same length as `x`).
/// * `bins` – number of bins per dimension.
#[allow(dead_code)]
pub fn mutual_information(x: &[f64], y: &[f64], bins: usize) -> f64 {
    let n = x.len().min(y.len());
    if n == 0 || bins == 0 {
        return 0.0;
    }
    let x = &x[..n];
    let y = &y[..n];
    let to_bin = |vals: &[f64]| -> Vec<usize> {
        let lo = vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if (hi - lo).abs() < f64::EPSILON {
            return vec![0usize; vals.len()];
        }
        let width = (hi - lo) / bins as f64;
        vals.iter()
            .map(|&v| ((v - lo) / width) as usize)
            .map(|idx| idx.min(bins - 1))
            .collect()
    };
    let bx = to_bin(x);
    let by = to_bin(y);
    let mut cx = vec![0usize; bins];
    let mut cy = vec![0usize; bins];
    let mut cxy = vec![0usize; bins * bins];
    for i in 0..n {
        cx[bx[i]] += 1;
        cy[by[i]] += 1;
        cxy[bx[i] * bins + by[i]] += 1;
    }
    let hx = shannon_entropy_from_counts(&cx);
    let hy = shannon_entropy_from_counts(&cy);
    let hxy = shannon_entropy_from_counts(&cxy);
    (hx + hy - hxy).max(0.0)
}
/// Compute the total free energy difference across a series of λ-windows using BAR.
///
/// Windows should be ordered from λ=0 to λ=1.
/// Returns the total ΔF = Σ ΔF(λᵢ → λᵢ₊₁).
#[allow(dead_code)]
pub fn total_fep_bar(windows: &[LambdaWindow], kt: f64) -> f64 {
    if windows.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0;
    for i in 0..windows.len() - 1 {
        total += windows[i].bar_estimate(&windows[i + 1], kt);
    }
    total
}
/// Compute total free energy difference using Zwanzig forward estimates.
#[allow(dead_code)]
pub fn total_fep_zwanzig(windows: &[LambdaWindow], kt: f64) -> f64 {
    if windows.len() < 2 {
        return 0.0;
    }
    windows[..windows.len() - 1]
        .iter()
        .map(|w| w.zwanzig_forward(kt))
        .sum()
}
/// Perform a two-stage alchemical transformation.
///
/// Stage 1: turn off charges (decharging).
/// Stage 2: turn off vdW interactions using soft-core potential samples.
///
/// # Arguments
/// * `decharge_windows` – λ-windows for the charge-off stage.
/// * `vdw_windows`      – λ-windows for the vdW-off stage.
/// * `kt`               – thermal energy kBT.
#[allow(dead_code)]
pub fn staged_alchemical_transformation(
    decharge_windows: &[LambdaWindow],
    vdw_windows: &[LambdaWindow],
    kt: f64,
) -> AlchemicalResult {
    let decharge_df = total_fep_bar(decharge_windows, kt);
    let vdw_df = total_fep_bar(vdw_windows, kt);
    AlchemicalResult {
        decharge_df,
        vdw_df,
        total_df: decharge_df + vdw_df,
    }
}
/// Soft-core Lennard-Jones potential used in FEP/TI to avoid end-point catastrophes.
///
/// U_sc = 4ε λ { \[α(1−λ)² + (r/σ)⁶\]⁻² − \[α(1−λ)² + (r/σ)⁶\]⁻¹ }
///
/// # Arguments
/// * `r`       – interparticle distance.
/// * `lambda`  – coupling parameter (0 = decoupled, 1 = fully coupled).
/// * `epsilon` – well depth.
/// * `sigma`   – particle radius.
/// * `alpha`   – soft-core parameter (typically 0.5).
#[allow(dead_code)]
pub fn soft_core_lj(r: f64, lambda: f64, epsilon: f64, sigma: f64, alpha: f64) -> f64 {
    if r < 0.0 {
        return 0.0;
    }
    let rs = (r / sigma).powi(6);
    let sc_term = alpha * (1.0 - lambda).powi(2) + rs;
    if sc_term < 1e-300 {
        return 0.0;
    }
    4.0 * epsilon * lambda * (sc_term.powi(-2) - sc_term.powi(-1))
}
/// Derivative of the soft-core LJ potential with respect to λ.
///
/// dU_sc/dλ — needed for thermodynamic integration.
#[allow(dead_code)]
pub fn soft_core_lj_dlambda(r: f64, lambda: f64, epsilon: f64, sigma: f64, alpha: f64) -> f64 {
    let rs = (r / sigma).powi(6);
    let sc_term = alpha * (1.0 - lambda).powi(2) + rs;
    if sc_term < 1e-300 {
        return 0.0;
    }
    let f = 4.0 * epsilon * (sc_term.powi(-2) - sc_term.powi(-1));
    let dsc_dl = -2.0 * alpha * (1.0 - lambda);
    let df_dl = 4.0 * epsilon * (-2.0 * sc_term.powi(-3) + sc_term.powi(-2)) * dsc_dl;
    f + lambda * df_dl
}
/// Compute the phase-space overlap matrix between adjacent λ-windows.
///
/// The overlap is estimated as `<exp(-ΔU/2kT)>` where ΔU are the ΔU_forward
/// samples from window i.  A value close to 1 indicates good overlap; close to
/// 0 indicates poor sampling.
///
/// Returns an (n-1)-element vector of pairwise overlap estimates.
#[allow(dead_code)]
pub fn overlap_matrix(windows: &[LambdaWindow], kt: f64) -> Vec<f64> {
    if windows.len() < 2 {
        return Vec::new();
    }
    windows
        .windows(2)
        .map(|pair| {
            let dus = &pair[0].du_forward;
            if dus.is_empty() {
                return 0.0;
            }
            let n = dus.len() as f64;
            let avg: f64 = dus.iter().map(|&du| (-du / (2.0 * kt)).exp()).sum::<f64>() / n;
            avg.min(1.0)
        })
        .collect()
}
/// Convergence monitor for an alchemical free energy calculation.
///
/// Splits the ΔU samples into successive blocks and tracks the running
/// estimate of ΔF.  Returns the vector of block-by-block ΔF estimates
/// (cumulative mean up to each block).
#[allow(dead_code)]
pub fn fep_convergence_trace(du_samples: &[f64], n_blocks: usize, kt: f64) -> Vec<f64> {
    if du_samples.is_empty() || n_blocks == 0 {
        return Vec::new();
    }
    let block_size = (du_samples.len() / n_blocks).max(1);
    (1..=n_blocks)
        .map(|b| {
            let end = (b * block_size).min(du_samples.len());
            fep_zwanzig(&du_samples[..end], kt)
        })
        .collect()
}
/// Estimate whether the FEP calculation has converged by comparing the
/// variance of block estimates in the first and second halves.
///
/// Returns `true` if the standard deviation of the second half is less than
/// `tol` times the absolute value of the overall estimate.
#[allow(dead_code)]
pub fn fep_has_converged(du_samples: &[f64], n_blocks: usize, kt: f64, tol: f64) -> bool {
    let trace = fep_convergence_trace(du_samples, n_blocks, kt);
    if trace.len() < 4 {
        return false;
    }
    let half = trace.len() / 2;
    let second_half = &trace[half..];
    let mean: f64 = second_half.iter().sum::<f64>() / second_half.len() as f64;
    let std: f64 = if second_half.len() < 2 {
        0.0
    } else {
        let var = second_half.iter().map(|&x| (x - mean).powi(2)).sum::<f64>()
            / (second_half.len() - 1) as f64;
        var.sqrt()
    };
    let abs_mean = mean.abs().max(1e-300);
    std / abs_mean < tol
}
/// WHAM self-consistent iteration to obtain the unbiased PMF.
///
/// Solves the WHAM equations:
/// ```text
/// P(ξ) = Σ_i n_i * h_i(ξ) / Σ_j N_j * exp(f_j - u_j(ξ)/kT)
/// exp(-f_i) = Σ_ξ P(ξ) * exp(-u_i(ξ)/kT)
/// ```
///
/// # Arguments
/// * `windows` – mutable slice of WHAM windows with accumulated samples.
/// * `bins`    – number of histogram bins along the reaction coordinate.
/// * `range`   – (min, max) range of the reaction coordinate.
/// * `kt`      – thermal energy k_B T.
/// * `tol`     – convergence tolerance on the free energies f_i.
/// * `max_iter`– maximum number of self-consistent iterations.
///
/// # Returns
/// A `Vec<(bin_center, pmf)>` of the unbiased PMF (shifted so minimum = 0).
#[allow(dead_code)]
pub fn wham_pmf(
    windows: &mut [WhamWindow],
    bins: usize,
    range: (f64, f64),
    kt: f64,
    tol: f64,
    max_iter: usize,
) -> Vec<(f64, f64)> {
    if windows.is_empty() || bins == 0 || range.0 >= range.1 {
        return Vec::new();
    }
    let (lo, hi) = range;
    let width = (hi - lo) / bins as f64;
    let bin_centers: Vec<f64> = (0..bins).map(|b| lo + (b as f64 + 0.5) * width).collect();
    let n_windows = windows.len();
    let n_i: Vec<f64> = windows.iter().map(|w| w.samples.len() as f64).collect();
    let mut hist: Vec<Vec<f64>> = vec![vec![0.0; bins]; n_windows];
    for (wi, w) in windows.iter().enumerate() {
        for &xi in &w.samples {
            if xi >= lo && xi < hi {
                let b = ((xi - lo) / width) as usize;
                let b = b.min(bins - 1);
                hist[wi][b] += 1.0;
            }
        }
    }
    for w in windows.iter_mut() {
        w.f_estimate = 0.0;
    }
    let mut p_unnorm = vec![0.0_f64; bins];
    for _iter in 0..max_iter {
        for b in 0..bins {
            let xi = bin_centers[b];
            let num: f64 = (0..n_windows).map(|i| hist[i][b]).sum();
            let denom: f64 = (0..n_windows)
                .map(|i| {
                    let u = windows[i].bias_energy(xi, kt);
                    n_i[i] * (windows[i].f_estimate - u).exp()
                })
                .sum::<f64>();
            p_unnorm[b] = if denom > 1e-300 { num / denom } else { 0.0 };
        }
        let mut max_delta = 0.0_f64;
        for i in 0..n_windows {
            let f_new_neg_exp: f64 = (0..bins)
                .map(|b| {
                    let xi = bin_centers[b];
                    let u = windows[i].bias_energy(xi, kt);
                    p_unnorm[b] * (-u).exp()
                })
                .sum();
            let f_new = if f_new_neg_exp > 1e-300 {
                -f_new_neg_exp.ln()
            } else {
                f64::INFINITY
            };
            let delta = (f_new - windows[i].f_estimate).abs();
            if delta > max_delta {
                max_delta = delta;
            }
            windows[i].f_estimate = f_new;
        }
        if max_delta < tol {
            break;
        }
    }
    let mut pmf_raw: Vec<(f64, f64)> = (0..bins)
        .filter(|&b| p_unnorm[b] > 0.0)
        .map(|b| (bin_centers[b], -kt * p_unnorm[b].ln()))
        .collect();
    if let Some(&(_, min_val)) = pmf_raw
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    {
        for entry in pmf_raw.iter_mut() {
            entry.1 -= min_val;
        }
    }
    pmf_raw
}
/// MBAR free energy estimates using self-consistent iteration.
///
/// Solves:
/// ```text
/// exp(-f_i) = Σ_n Σ_j [N_j * exp(f_j - u_j(x_n))] / Σ_l N_l * exp(f_l - u_l(x_n))
/// ```
///
/// # Arguments
/// * `u_kn`     – K × N_total matrix of reduced potentials. `u_kn[k][n]` = reduced
///   potential of sample n evaluated at state k.
/// * `n_k`      – number of samples from each state k (length K).
/// * `tol`      – convergence tolerance on free energies.
/// * `max_iter` – maximum iterations.
///
/// # Returns
/// Free energy estimates `f_k` (length K), normalised so `f[0] = 0`.
#[allow(dead_code)]
pub fn mbar_free_energies(u_kn: &[Vec<f64>], n_k: &[usize], tol: f64, max_iter: usize) -> Vec<f64> {
    let k = u_kn.len();
    if k == 0 || n_k.len() != k {
        return Vec::new();
    }
    let n_total: usize = n_k.iter().sum();
    if n_total == 0 {
        return vec![0.0; k];
    }
    let n_k_f: Vec<f64> = n_k.iter().map(|&n| n as f64).collect();
    let mut f = vec![0.0_f64; k];
    for _iter in 0..max_iter {
        let mut f_new = vec![0.0_f64; k];
        for ki in 0..k {
            let mut sum = 0.0_f64;
            for n in 0..n_total {
                let log_denom_terms: Vec<f64> = (0..k)
                    .map(|l| {
                        n_k_f[l].ln() + f[l] - u_kn[l].get(n).copied().unwrap_or(f64::INFINITY)
                    })
                    .collect();
                let max_term = log_denom_terms
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max);
                if max_term.is_finite() {
                    let denom = max_term.exp()
                        * log_denom_terms
                            .iter()
                            .map(|&v| (v - max_term).exp())
                            .sum::<f64>();
                    if denom > 1e-300 {
                        let u_ki_n = u_kn[ki].get(n).copied().unwrap_or(f64::INFINITY);
                        sum += (-u_ki_n).exp() / denom;
                    }
                }
            }
            f_new[ki] = if sum > 1e-300 {
                -sum.ln()
            } else {
                f64::INFINITY
            };
        }
        let f0 = f_new[0];
        for v in f_new.iter_mut() {
            *v -= f0;
        }
        let delta: f64 = f
            .iter()
            .zip(f_new.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        f = f_new;
        if delta < tol {
            break;
        }
    }
    f
}
/// Compute MBAR weights for each sample: the probability that sample n
/// belongs to a target (unbiased) state 0 given the MBAR free energies.
///
/// `W_n = 1 / (Σ_k N_k * exp(f_k - u_k(x_n)))` (unnormalised)
#[allow(dead_code)]
pub fn mbar_weights(u_kn: &[Vec<f64>], n_k: &[usize], f: &[f64]) -> Vec<f64> {
    let k = u_kn.len();
    let n_total: usize = n_k.iter().sum();
    let n_k_f: Vec<f64> = n_k.iter().map(|&n| n as f64).collect();
    let mut weights = Vec::with_capacity(n_total);
    for n in 0..n_total {
        let log_denom_terms: Vec<f64> = (0..k)
            .map(|l| n_k_f[l].ln() + f[l] - u_kn[l].get(n).copied().unwrap_or(f64::INFINITY))
            .collect();
        let max_term = log_denom_terms
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let w = if max_term.is_finite() {
            let denom = max_term.exp()
                * log_denom_terms
                    .iter()
                    .map(|&v| (v - max_term).exp())
                    .sum::<f64>();
            if denom > 1e-300 { 1.0 / denom } else { 0.0 }
        } else {
            0.0
        };
        weights.push(w);
    }
    weights
}
/// Compute a MBAR observable expectation value ⟨A⟩ using pre-computed weights.
///
/// `⟨A⟩ = Σ_n W_n * A_n / Σ_n W_n`
#[allow(dead_code)]
pub fn mbar_observable(weights: &[f64], observable: &[f64]) -> f64 {
    let n = weights.len().min(observable.len());
    if n == 0 {
        return 0.0;
    }
    let sum_w: f64 = weights[..n].iter().sum();
    if sum_w < 1e-300 {
        return 0.0;
    }
    let sum_aw: f64 = weights[..n]
        .iter()
        .zip(observable[..n].iter())
        .map(|(w, a)| w * a)
        .sum();
    sum_aw / sum_w
}
/// Extract free energy differences from REMD (Replica Exchange MD) exchange statistics.
///
/// Given acceptance ratios `alpha[i]` between replicas i and i+1 and their temperature
/// ratio `beta_ratio[i]` = β_{i+1}/β_i, estimates the ΔF between adjacent replicas
/// using the BAR-like relation:
///
/// `ΔF_{i→i+1} ≈ kT_i * ln(alpha[i] / (1 - alpha[i]))`
///
/// This is a rough estimate; proper analysis requires BAR with full ΔU samples.
#[allow(dead_code)]
pub fn remd_free_energy_estimates(acceptance_ratios: &[f64], kt: f64) -> Vec<f64> {
    acceptance_ratios
        .iter()
        .map(|&alpha| {
            let clamped = alpha.clamp(1e-10, 1.0 - 1e-10);
            kt * (clamped / (1.0 - clamped)).ln()
        })
        .collect()
}
/// Bootstrap uncertainty estimate for a free energy calculation.
///
/// Re-samples the `du_samples` vector with replacement `n_bootstrap` times,
/// applies `estimator`, and returns the standard deviation of the bootstrap estimates.
///
/// # Arguments
/// * `du_samples`    – ΔU samples.
/// * `n_bootstrap`   – number of bootstrap replicas.
/// * `kt`            – thermal energy k_B T.
/// * `seed_offset`   – deterministic seed for reproducible testing (uses simple LCG).
#[allow(dead_code)]
pub fn bootstrap_uncertainty(
    du_samples: &[f64],
    n_bootstrap: usize,
    kt: f64,
    seed_offset: u64,
) -> f64 {
    let n = du_samples.len();
    if n < 2 || n_bootstrap < 2 {
        return 0.0;
    }
    let mut state = 6364136223846793005_u64.wrapping_add(seed_offset);
    let mut bootstrap_estimates = Vec::with_capacity(n_bootstrap);
    for _ in 0..n_bootstrap {
        let resample: Vec<f64> = (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let idx = ((state >> 33) as usize) % n;
                du_samples[idx]
            })
            .collect();
        bootstrap_estimates.push(fep_zwanzig(&resample, kt));
    }
    let mean = bootstrap_estimates.iter().sum::<f64>() / n_bootstrap as f64;
    let var = bootstrap_estimates
        .iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f64>()
        / (n_bootstrap - 1) as f64;
    var.sqrt()
}
/// TI free energy with Gaussian quadrature weights.
///
/// Instead of the trapezoidal rule, uses Gauss-Legendre-like weights pre-supplied
/// by the caller (e.g. from 5-point or 9-point quadrature rules).
///
/// # Arguments
/// * `du_dlambda` – ⟨dH/dλ⟩ at each quadrature point (same length as `weights`).
/// * `weights`    – quadrature weights summing to 1.
#[allow(dead_code)]
pub fn ti_gaussian_quadrature(du_dlambda: &[f64], weights: &[f64]) -> f64 {
    let n = du_dlambda.len().min(weights.len());
    if n == 0 {
        return 0.0;
    }
    du_dlambda[..n]
        .iter()
        .zip(weights[..n].iter())
        .map(|(d, w)| d * w)
        .sum()
}
/// Five-point Gauss-Legendre quadrature nodes and weights on \[0,1\].
///
/// Returns `(nodes, weights)` for integration ∫₀¹ f(λ) dλ.
#[allow(dead_code)]
pub fn gauss_legendre_5pt() -> ([f64; 5], [f64; 5]) {
    let nodes_11: [f64; 5] = [
        -0.906_179_845_9,
        -0.538_469_310_1,
        0.0,
        0.538_469_310_1,
        0.906_179_845_9,
    ];
    let weights_11: [f64; 5] = [
        0.236_926_885_1,
        0.478_628_670_5,
        0.568_888_888_9,
        0.478_628_670_5,
        0.236_926_885_1,
    ];
    let mut nodes = [0.0; 5];
    let mut weights = [0.0; 5];
    for i in 0..5 {
        nodes[i] = 0.5 * (nodes_11[i] + 1.0);
        weights[i] = 0.5 * weights_11[i];
    }
    (nodes, weights)
}
/// Block-averaged estimate of ⟨dH/dλ⟩ and its standard error for one TI window.
///
/// Splits `samples` into blocks of size `block_size`, computes the block mean,
/// and returns `(mean, std_err)`.
#[allow(dead_code)]
pub fn ti_block_average(samples: &[f64], block_size: usize) -> (f64, f64) {
    let n = samples.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if block_size == 0 || block_size >= n {
        let mean = samples.iter().sum::<f64>() / n as f64;
        return (mean, 0.0);
    }
    let blocks: Vec<f64> = samples
        .chunks(block_size)
        .filter(|b| !b.is_empty())
        .map(|b| b.iter().sum::<f64>() / b.len() as f64)
        .collect();
    let nb = blocks.len() as f64;
    let mean = blocks.iter().sum::<f64>() / nb;
    let var = if blocks.len() < 2 {
        0.0
    } else {
        blocks.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (nb - 1.0)
    };
    (mean, (var / nb).sqrt())
}
/// Statistical uncertainty of the BAR estimate using a bootstrap-like formula.
///
/// Uses the asymptotic formula from Bennett (1976):
/// ```text
/// σ²(ΔF) ≈ 1/(N₀ · ⟨f²⟩₀ - ⟨f⟩₀²) + 1/(N₁ · ⟨f²⟩₁ - ⟨f⟩₁²)
/// ```
/// where f = fermi((ΔU ± ΔF) / kT).
///
/// Returns the estimated standard deviation of the BAR ΔF estimate.
#[allow(dead_code)]
pub fn bar_uncertainty(du_forward: &[f64], du_backward: &[f64], df: f64, kt: f64) -> f64 {
    let n0 = du_forward.len();
    let n1 = du_backward.len();
    if n0 < 2 || n1 < 2 || kt == 0.0 {
        return 0.0;
    }
    let ln_ratio = ((n0 as f64) / (n1 as f64)).ln();
    let c = df / kt - ln_ratio;
    let fermi_vals_0: Vec<f64> = du_forward.iter().map(|&du| fermi(-du / kt - c)).collect();
    let fermi_vals_1: Vec<f64> = du_backward.iter().map(|&du| fermi(du / kt + c)).collect();
    let var_contrib = |fv: &[f64]| -> f64 {
        let nf = fv.len() as f64;
        let mean = fv.iter().sum::<f64>() / nf;
        let mean_sq = fv.iter().map(|x| x * x).sum::<f64>() / nf;
        let var = (mean_sq - mean * mean).max(0.0);
        if var < 1e-300 {
            f64::INFINITY
        } else {
            1.0 / (nf * var)
        }
    };
    let v0 = var_contrib(&fermi_vals_0);
    let v1 = var_contrib(&fermi_vals_1);
    ((v0 + v1) * kt * kt).sqrt()
}
/// Compute the umbrella bias potential energy.
///
/// `U_bias(ξ) = 0.5 * k * (ξ - ξ_ref)²`
#[allow(dead_code)]
pub fn umbrella_bias_energy(xi: f64, xi_ref: f64, k: f64) -> f64 {
    0.5 * k * (xi - xi_ref).powi(2)
}
/// Compute the umbrella bias force along the reaction coordinate.
///
/// `F_bias = -dU_bias/dξ = -k * (ξ - ξ_ref)`
#[allow(dead_code)]
pub fn umbrella_bias_force(xi: f64, xi_ref: f64, k: f64) -> f64 {
    -k * (xi - xi_ref)
}
/// Unbias a histogram from an umbrella window.
///
/// Given counts in bins and the bias potential, returns the unbiased (flat)
/// probability estimates: `P_unbiased(ξ) ∝ P_biased(ξ) * exp(+U_bias(ξ)/kT)`.
///
/// Returns normalised probabilities (sum = 1) for occupied bins.
#[allow(dead_code)]
pub fn unbias_histogram(
    bin_centers: &[f64],
    counts: &[u64],
    xi_ref: f64,
    k_bias: f64,
    kt: f64,
) -> Vec<f64> {
    let n = bin_centers.len().min(counts.len());
    if n == 0 || kt == 0.0 {
        return vec![0.0; n];
    }
    let total: u64 = counts[..n].iter().sum();
    if total == 0 {
        return vec![0.0; n];
    }
    let tf = total as f64;
    let mut unbiased: Vec<f64> = (0..n)
        .map(|i| {
            if counts[i] == 0 {
                0.0
            } else {
                let p_biased = counts[i] as f64 / tf;
                let u_b = umbrella_bias_energy(bin_centers[i], xi_ref, k_bias);
                p_biased * (u_b / kt).exp()
            }
        })
        .collect();
    let sum: f64 = unbiased.iter().sum();
    if sum > 1e-300 {
        for v in unbiased.iter_mut() {
            *v /= sum;
        }
    }
    unbiased
}
/// Jarzynski free energy from multiple steered MD work values.
///
/// Convenience wrapper that also returns the standard error of the
/// exponential average.
///
/// Returns `(delta_f, std_err_delta_f, n_used)` where `n_used` counts finite
/// work values.
#[allow(dead_code)]
pub fn jarzynski_from_smd_runs(works: &[f64], kt: f64) -> (f64, f64, usize) {
    if works.is_empty() || kt == 0.0 {
        return (0.0, 0.0, 0);
    }
    let finite: Vec<f64> = works.iter().cloned().filter(|w| w.is_finite()).collect();
    let n = finite.len();
    if n == 0 {
        return (f64::INFINITY, f64::INFINITY, 0);
    }
    let nf = n as f64;
    let exps: Vec<f64> = finite.iter().map(|&w| (-w / kt).exp()).collect();
    let mean_exp = exps.iter().sum::<f64>() / nf;
    let df = if mean_exp > 0.0 {
        -kt * mean_exp.ln()
    } else {
        f64::INFINITY
    };
    let var_exp = if n < 2 {
        0.0
    } else {
        exps.iter().map(|x| (x - mean_exp).powi(2)).sum::<f64>() / (nf - 1.0)
    };
    let se_exp = (var_exp / nf).sqrt();
    let se_df = if mean_exp > 1e-300 {
        kt * se_exp / mean_exp
    } else {
        f64::INFINITY
    };
    (df, se_df, n)
}
/// Metropolis swap acceptance probability between two replicas.
///
/// `P_accept = min(1, exp((β_i - β_j)(E_i - E_j)))`
///
/// where β = 1/(k_B T).
#[allow(dead_code)]
pub fn remd_swap_acceptance(
    energy_i: f64,
    energy_j: f64,
    temp_i: f64,
    temp_j: f64,
    kb: f64,
) -> f64 {
    if temp_i <= 0.0 || temp_j <= 0.0 || kb <= 0.0 {
        return 0.0;
    }
    let beta_i = 1.0 / (kb * temp_i);
    let beta_j = 1.0 / (kb * temp_j);
    let exponent = (beta_i - beta_j) * (energy_i - energy_j);
    exponent.exp().min(1.0)
}
/// Optimal temperature spacing for REMD.
///
/// Given `n_replicas`, `t_min`, and `t_max`, returns the geometric temperature
/// ladder that is expected to give ~25% swap acceptance (the rule of thumb).
#[allow(dead_code)]
pub fn remd_temperature_ladder(n_replicas: usize, t_min: f64, t_max: f64) -> Vec<f64> {
    if n_replicas < 2 || t_min <= 0.0 || t_max <= t_min {
        return vec![t_min];
    }
    let ratio = (t_max / t_min).powf(1.0 / (n_replicas - 1) as f64);
    (0..n_replicas)
        .map(|i| t_min * ratio.powi(i as i32))
        .collect()
}
/// Compute the expected swap acceptance probability for a geometric ladder
/// assuming harmonic fluctuations:
///
/// `P_swap ≈ erfc(sqrt(C_v * (T_{i+1}/T_i - 1)^2 / (2 * C_v)))` (rough estimate)
///
/// In practice, for harmonic systems: `P ≈ exp(-C_v * (β_i - β_j)^2 * kT^2 / 2)`
/// where C_v is the heat capacity per degree of freedom.
#[allow(dead_code)]
pub fn remd_expected_swap_probability(temp_i: f64, temp_j: f64, n_dof: f64, kb: f64) -> f64 {
    if temp_i <= 0.0 || temp_j <= 0.0 || kb <= 0.0 || n_dof <= 0.0 {
        return 0.0;
    }
    let beta_i = 1.0 / (kb * temp_i);
    let beta_j = 1.0 / (kb * temp_j);
    let cv = 0.5 * n_dof * kb;
    let sigma2_e = kb * kb * temp_i * temp_i * cv / kb;
    let exponent = -0.5 * (beta_i - beta_j).powi(2) * sigma2_e;
    exponent.exp()
}
/// Compute mean and standard error of a sample array.
///
/// Returns `(mean, stderr)`.  Standard error = std / sqrt(N).
#[allow(dead_code)]
pub(super) fn block_mean_stderr(samples: &[f64]) -> (f64, f64) {
    let n = samples.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean = samples.iter().sum::<f64>() / n as f64;
    if n < 2 {
        return (mean, 0.0);
    }
    let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    (mean, (var / n as f64).sqrt())
}
/// Cubic spline natural boundary condition — tridiagonal solver for TI.
///
/// Accepts `n` node values `y` at positions `x` and returns the second
/// derivatives `m[i]` (spline coefficients) for natural splines (m\[0\]=m\[n-1\]=0).
#[allow(dead_code)]
pub fn natural_spline_second_deriv(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    if n < 2 {
        return vec![0.0; n];
    }
    let mut h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let mut alpha: Vec<f64> = vec![0.0; n];
    for i in 1..n - 1 {
        alpha[i] = 3.0 / h[i] * (y[i + 1] - y[i]) - 3.0 / h[i - 1] * (y[i] - y[i - 1]);
    }
    let mut l = vec![1.0_f64; n];
    let mut mu = vec![0.0_f64; n];
    let mut z = vec![0.0_f64; n];
    for i in 1..n - 1 {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    let mut m = vec![0.0_f64; n];
    for j in (0..n - 1).rev() {
        m[j] = z[j] - mu[j] * m[j + 1];
    }
    let _ = h.pop();
    m
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_fep_zwanzig_zero_du() {
        let samples = vec![0.0; 100];
        let df = fep_zwanzig(&samples, 1.0);
        assert!(df.abs() < 1e-12, "expected 0, got {df}");
    }
    #[test]
    fn test_fep_zwanzig_constant_du() {
        let c = 2.5;
        let samples = vec![c; 500];
        let df = fep_zwanzig(&samples, 1.0);
        assert!((df - c).abs() < 1e-10, "expected {c}, got {df}");
    }
    #[test]
    fn test_fep_zwanzig_empty() {
        let df = fep_zwanzig(&[], 1.0);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_fep_zwanzig_analytical() {
        let kt = 2.0;
        let du = 3.0;
        let samples = vec![du; 1000];
        let df = fep_zwanzig(&samples, kt);
        assert!((df - du).abs() < 1e-10, "expected {du}, got {df}");
    }
    #[test]
    fn test_fep_bar_symmetric() {
        let fwd: Vec<f64> = vec![1.0, -1.0, 2.0, -2.0];
        let bwd: Vec<f64> = vec![-1.0, 1.0, -2.0, 2.0];
        let df = fep_bar(&fwd, &bwd, 1.0, 1e-10);
        assert!(df.abs() < 1e-4, "symmetric BAR should give ΔF≈0, got {df}");
    }
    #[test]
    fn test_fep_bar_empty_forward() {
        let bwd = vec![0.0; 50];
        let df = fep_bar(&[], &bwd, 1.0, 1e-10);
        assert!(df.abs() < 1e-10, "expected 0, got {df}");
    }
    #[test]
    fn test_fep_bar_empty_backward() {
        let fwd = vec![0.0; 50];
        let df = fep_bar(&fwd, &[], 1.0, 1e-10);
        assert!(df.abs() < 1e-10, "expected 0, got {df}");
    }
    #[test]
    fn test_fep_bar_both_empty() {
        let df = fep_bar(&[], &[], 1.0, 1e-10);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_fep_bar_constant_shift() {
        let c = 1.5_f64;
        let fwd = vec![c; 200];
        let df = fep_bar(&fwd, &[], 1.0, 1e-10);
        assert!((df - c).abs() < 1e-10, "expected ≈{c}, got {df}");
    }
    #[test]
    fn test_fep_uncertainty_constant() {
        let samples = vec![1.0; 100];
        let unc = fep_uncertainty(&samples, 1.0);
        assert!(
            unc.abs() < 1e-10,
            "constant samples should give zero uncertainty, got {unc}"
        );
    }
    #[test]
    fn test_fep_uncertainty_few_samples() {
        let samples = vec![1.0, 2.0];
        let unc = fep_uncertainty(&samples, 1.0);
        assert_eq!(unc, 0.0);
    }
    #[test]
    fn test_ti_constant() {
        let lambdas: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();
        let du_dl = vec![2.0; 11];
        let df = thermodynamic_integration(&lambdas, &du_dl);
        assert!(
            (df - 2.0).abs() < 1e-12,
            "TI with dU/dλ=2 should give ΔF=2, got {df}"
        );
    }
    #[test]
    fn test_ti_linear() {
        let lambdas: Vec<f64> = (0..=100).map(|i| i as f64 / 100.0).collect();
        let du_dl: Vec<f64> = lambdas.clone();
        let df = thermodynamic_integration(&lambdas, &du_dl);
        assert!(
            (df - 0.5).abs() < 1e-4,
            "TI with dU/dλ=λ should give ΔF=0.5, got {df}"
        );
    }
    #[test]
    fn test_ti_single_window() {
        let lambdas = vec![0.5];
        let du_dl = vec![1.0];
        let df = thermodynamic_integration(&lambdas, &du_dl);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_ti_uncertainty_zero_errors() {
        let lambdas = vec![0.0, 0.5, 1.0];
        let errs = vec![0.0; 3];
        let unc = ti_uncertainty(&lambdas, &errs);
        assert_eq!(unc, 0.0);
    }
    #[test]
    fn test_ti_uncertainty_nonzero() {
        let lambdas = vec![0.0, 1.0];
        let errs = vec![0.1, 0.1];
        let unc = ti_uncertainty(&lambdas, &errs);
        assert!((unc - 0.1).abs() < 1e-12, "expected 0.1, got {unc}");
    }
    #[test]
    fn test_jarzynski_zero_work() {
        let work = vec![0.0; 200];
        let df = jarzynski_estimate(&work, 1.0);
        assert!(df.abs() < 1e-12, "zero work → ΔF=0, got {df}");
    }
    #[test]
    fn test_cumulant_2nd_zero_variance() {
        let c = 3.0;
        let work = vec![c; 100];
        let df = cumulant_expansion_2nd(&work, 1.0);
        assert!((df - c).abs() < 1e-12, "expected {c}, got {df}");
    }
    #[test]
    fn test_cumulant_2nd_variance_correction() {
        let mean = 5.0_f64;
        let var = 2.0_f64;
        let n = 10000usize;
        let d = var.sqrt();
        let mut work: Vec<f64> = (0..n / 2).map(|_| mean + d).collect();
        work.extend((0..n / 2).map(|_| mean - d));
        let df = cumulant_expansion_2nd(&work, 1.0);
        let expected = mean - var / 2.0;
        assert!(
            (df - expected).abs() < 0.01,
            "expected ≈{expected}, got {df}"
        );
    }
    #[test]
    fn test_histogram_pmf_uniform() {
        let samples: Vec<f64> = (0..1000).map(|i| i as f64 / 999.0).collect();
        let pmf = histogram_pmf(&samples, 10, (0.0, 1.0), 1.0);
        assert_eq!(pmf.len(), 10, "all bins should be occupied");
        for &(_, v) in &pmf {
            assert!(
                v.abs() < 0.5,
                "uniform PMF should be flat after shift, got {v}"
            );
        }
    }
    #[test]
    fn test_histogram_pmf_empty() {
        let pmf = histogram_pmf(&[], 10, (0.0, 1.0), 1.0);
        assert!(pmf.is_empty());
    }
    #[test]
    fn test_histogram_pmf_min_is_zero() {
        let samples: Vec<f64> = (0..500).map(|i| i as f64).collect();
        let pmf = histogram_pmf(&samples, 5, (0.0, 500.0), 1.0);
        let min_val = pmf.iter().map(|&(_, v)| v).fold(f64::INFINITY, f64::min);
        assert!(
            min_val.abs() < 1e-12,
            "minimum PMF should be 0 after shift, got {min_val}"
        );
    }
    #[test]
    fn test_pmf_uncertainty_basic() {
        let pmf = vec![(0.0, 0.0), (0.1, 0.5), (0.2, 1.0)];
        let unc = pmf_uncertainty(&pmf, 300);
        assert_eq!(unc.len(), 3);
        for u in &unc {
            assert!(*u > 0.0, "uncertainty should be positive");
        }
    }
    #[test]
    fn test_configurational_entropy_uniform() {
        let samples: Vec<f64> = (0..10000).map(|i| i as f64 / 9999.0).collect();
        let s = configurational_entropy(&samples, 10, 1.0);
        assert!(s.is_finite(), "entropy should be finite, got {s}");
    }
    #[test]
    fn test_configurational_entropy_constant() {
        let samples = vec![1.0; 100];
        let s = configurational_entropy(&samples, 10, 1.0);
        assert_eq!(s, 0.0);
    }
    #[test]
    fn test_mutual_information_independent() {
        let n = 1000;
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&v| v * 2.0).collect();
        let mi = mutual_information(&x, &y, 20);
        assert!(
            mi > 0.0,
            "correlated variables should have MI > 0, got {mi}"
        );
    }
    #[test]
    fn test_mutual_information_nonnegative() {
        let x: Vec<f64> = (0..200).map(|i| (i as f64).sin()).collect();
        let y: Vec<f64> = (0..200).map(|i| (i as f64 * 0.7).cos()).collect();
        let mi = mutual_information(&x, &y, 10);
        assert!(mi >= 0.0, "MI must be non-negative, got {mi}");
    }
    #[test]
    fn test_kl_divergence_identical() {
        let p = vec![0.25, 0.25, 0.25, 0.25];
        let kl = KlDivergence::kl_divergence(&p, &p);
        assert!(kl.abs() < 1e-12, "D_KL(P‖P) = 0, got {kl}");
    }
    #[test]
    fn test_kl_divergence_analytical() {
        let p = vec![0.5, 0.5];
        let q = vec![0.25, 0.75];
        let expected = 0.5 * (2.0_f64).ln() + 0.5 * (2.0_f64 / 3.0).ln();
        let kl = KlDivergence::kl_divergence(&p, &q);
        assert!(
            (kl - expected).abs() < 1e-12,
            "expected {expected}, got {kl}"
        );
    }
    #[test]
    fn test_kl_divergence_infinite() {
        let p = vec![0.5, 0.5];
        let q = vec![0.0, 1.0];
        let kl = KlDivergence::kl_divergence(&p, &q);
        assert_eq!(kl, f64::INFINITY);
    }
    #[test]
    fn test_js_divergence_identical() {
        let p = vec![0.3, 0.4, 0.3];
        let js = KlDivergence::js_divergence(&p, &p);
        assert!(js.abs() < 1e-12, "JS(P‖P) = 0, got {js}");
    }
    #[test]
    fn test_js_divergence_bounded() {
        let p = vec![1.0, 0.0];
        let q = vec![0.0, 1.0];
        let js = KlDivergence::js_divergence(&p, &q);
        let ln2 = 2.0_f64.ln();
        assert!(js <= ln2 + 1e-12, "JS should be ≤ ln(2), got {js}");
        assert!(js >= 0.0, "JS must be non-negative, got {js}");
    }
    #[test]
    fn test_js_divergence_symmetric() {
        let p = vec![0.6, 0.4];
        let q = vec![0.3, 0.7];
        let js_pq = KlDivergence::js_divergence(&p, &q);
        let js_qp = KlDivergence::js_divergence(&q, &p);
        assert!(
            (js_pq - js_qp).abs() < 1e-12,
            "JS should be symmetric, got {js_pq} vs {js_qp}"
        );
    }
    #[test]
    fn test_lambda_window_zwanzig_zero_du() {
        let mut win = LambdaWindow::new(0.5);
        for _ in 0..50 {
            win.add_forward(0.0);
        }
        let df = win.zwanzig_forward(1.0);
        assert!(df.abs() < 1e-10, "ΔF with all ΔU=0 should be 0, got {df}");
    }
    #[test]
    fn test_lambda_window_n_samples() {
        let mut win = LambdaWindow::new(0.0);
        win.add_forward(1.0);
        win.add_forward(2.0);
        win.add_backward(0.5);
        assert_eq!(win.n_samples(), 2, "n_samples should be max(fwd, bwd) = 2");
    }
    #[test]
    fn test_total_fep_bar_zero_du() {
        let n_windows = 5;
        let mut windows: Vec<LambdaWindow> = (0..n_windows)
            .map(|i| LambdaWindow::new(i as f64 / (n_windows - 1) as f64))
            .collect();
        for i in 0..n_windows - 1 {
            for _ in 0..50 {
                windows[i].add_forward(0.0);
                windows[i + 1].add_backward(0.0);
            }
        }
        let total = total_fep_bar(&windows, 1.0);
        assert!(
            total.abs() < 1e-8,
            "Total BAR with all ΔU=0 should be 0, got {total}"
        );
    }
    #[test]
    fn test_total_fep_zwanzig_constant_du() {
        let c = 1.0_f64;
        let n = 4;
        let mut windows: Vec<LambdaWindow> = (0..n)
            .map(|i| LambdaWindow::new(i as f64 / (n - 1) as f64))
            .collect();
        for i in 0..n - 1 {
            for _ in 0..100 {
                windows[i].add_forward(c);
            }
        }
        let total = total_fep_zwanzig(&windows, 1.0);
        assert!(
            (total - (n - 1) as f64 * c).abs() < 1e-8,
            "expected {}, got {total}",
            (n - 1) as f64 * c
        );
    }
    #[test]
    fn test_soft_core_lj_at_lambda_zero() {
        let u = soft_core_lj(1.0, 0.0, 1.0, 1.0, 0.5);
        assert!(u.abs() < 1e-12, "soft-core LJ at λ=0 should be 0, got {u}");
    }
    #[test]
    fn test_soft_core_lj_finite_at_r_zero() {
        let u = soft_core_lj(0.0, 1.0, 1.0, 1.0, 0.5);
        assert!(
            u.is_finite(),
            "soft-core LJ at r=0 should be finite, got {u}"
        );
    }
    #[test]
    fn test_soft_core_lj_dlambda_finite() {
        let du = soft_core_lj_dlambda(1.5, 0.5, 1.0, 1.0, 0.5);
        assert!(du.is_finite(), "dU_sc/dλ should be finite, got {du}");
    }
    #[test]
    fn test_overlap_matrix_length() {
        let windows: Vec<LambdaWindow> = (0..4)
            .map(|i| {
                let mut w = LambdaWindow::new(i as f64 / 3.0);
                for _ in 0..20 {
                    w.add_forward(0.5);
                }
                w
            })
            .collect();
        let overlaps = overlap_matrix(&windows, 1.0);
        assert_eq!(overlaps.len(), 3, "n-1 overlaps for n windows");
    }
    #[test]
    fn test_overlap_matrix_range() {
        let windows: Vec<LambdaWindow> = (0..3)
            .map(|i| {
                let mut w = LambdaWindow::new(i as f64 / 2.0);
                for _ in 0..50 {
                    w.add_forward(0.0);
                }
                w
            })
            .collect();
        let overlaps = overlap_matrix(&windows, 1.0);
        for &ov in &overlaps {
            assert!(
                (0.0..=1.0 + 1e-10).contains(&ov),
                "Overlap should be in [0,1], got {ov}"
            );
        }
    }
    #[test]
    fn test_fep_convergence_trace_length() {
        let samples = vec![1.0; 100];
        let trace = fep_convergence_trace(&samples, 10, 1.0);
        assert_eq!(trace.len(), 10, "Trace should have n_blocks entries");
    }
    #[test]
    fn test_fep_convergence_trace_monotone_constant() {
        let c = 2.0_f64;
        let samples = vec![c; 200];
        let trace = fep_convergence_trace(&samples, 10, 1.0);
        for &v in &trace {
            assert!((v - c).abs() < 1e-8, "Trace should be flat at {c}, got {v}");
        }
    }
    #[test]
    fn test_fep_has_converged_constant() {
        let samples = vec![1.0; 500];
        let converged = fep_has_converged(&samples, 10, 1.0, 0.01);
        assert!(converged, "Constant samples should indicate convergence");
    }
    #[test]
    fn test_staged_alchemical_total_is_sum() {
        let make_windows = || -> Vec<LambdaWindow> {
            (0..3)
                .map(|i| {
                    let mut w = LambdaWindow::new(i as f64 / 2.0);
                    for _ in 0..50 {
                        w.add_forward(0.0);
                        w.add_backward(0.0);
                    }
                    w
                })
                .collect()
        };
        let decharge = make_windows();
        let vdw = make_windows();
        let result = staged_alchemical_transformation(&decharge, &vdw, 1.0);
        assert!(
            (result.total_df - result.decharge_df - result.vdw_df).abs() < 1e-12,
            "total_df should equal sum of stages"
        );
        assert!(result.total_df.is_finite(), "total_df should be finite");
    }
    #[test]
    fn test_wham_window_bias_energy_at_ref() {
        let w = WhamWindow::new(2.0, 10.0);
        let e = w.bias_energy(2.0, 1.0);
        assert!(
            e.abs() < 1e-12,
            "Bias energy at xi=xi_ref should be 0, got {e}"
        );
    }
    #[test]
    fn test_wham_window_bias_energy_positive_away() {
        let w = WhamWindow::new(0.0, 10.0);
        let e = w.bias_energy(1.0, 1.0);
        assert!(
            e > 0.0,
            "Bias energy away from ref should be positive, got {e}"
        );
    }
    #[test]
    fn test_wham_pmf_empty_windows() {
        let mut windows: Vec<WhamWindow> = Vec::new();
        let pmf = wham_pmf(&mut windows, 10, (0.0, 1.0), 1.0, 1e-6, 100);
        assert!(pmf.is_empty(), "PMF from empty windows should be empty");
    }
    #[test]
    fn test_wham_pmf_flat_distribution() {
        let mut windows = vec![WhamWindow::new(0.25, 1.0), WhamWindow::new(0.75, 1.0)];
        for w in &mut windows {
            for i in 0..50 {
                w.add_sample(i as f64 / 49.0);
            }
        }
        let pmf = wham_pmf(&mut windows, 5, (0.0, 1.0), 1.0, 1e-4, 100);
        assert!(!pmf.is_empty(), "WHAM PMF should return some bins");
        for &(_, v) in &pmf {
            assert!(v >= -0.1, "WHAM PMF bin should be non-negative, got {v}");
        }
    }
    #[test]
    fn test_wham_pmf_min_is_zero() {
        let mut windows = vec![WhamWindow::new(0.5, 5.0)];
        for i in 0..100 {
            windows[0].add_sample(0.0 + i as f64 / 99.0);
        }
        let pmf = wham_pmf(&mut windows, 10, (0.0, 1.0), 1.0, 1e-4, 50);
        if !pmf.is_empty() {
            let min_val = pmf.iter().map(|&(_, v)| v).fold(f64::INFINITY, f64::min);
            assert!(
                min_val.abs() < 1e-10,
                "PMF minimum should be 0, got {min_val}"
            );
        }
    }
    #[test]
    fn test_mbar_free_energies_single_state() {
        let u_kn: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0]];
        let n_k = vec![3usize];
        let f = mbar_free_energies(&u_kn, &n_k, 1e-8, 100);
        assert_eq!(f.len(), 1);
        assert!(
            f[0].abs() < 1e-10,
            "Single-state f should be 0, got {}",
            f[0]
        );
    }
    #[test]
    fn test_mbar_free_energies_identical_states() {
        let u_row: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let u_kn = vec![u_row.clone(), u_row];
        let n_k = vec![10usize, 10];
        let f = mbar_free_energies(&u_kn, &n_k, 1e-8, 200);
        assert_eq!(f.len(), 2);
        assert!(f[0].abs() < 1e-10, "f[0] should be 0");
        assert!(
            f[1].abs() < 1.0,
            "Identical states: |ΔF| should be small, got {}",
            f[1]
        );
    }
    #[test]
    fn test_mbar_free_energies_empty() {
        let f = mbar_free_energies(&[], &[], 1e-8, 100);
        assert!(
            f.is_empty(),
            "MBAR with empty input should return empty vec"
        );
    }
    #[test]
    fn test_mbar_weights_sum() {
        let u_kn: Vec<Vec<f64>> = vec![vec![0.5, 1.0, 1.5, 2.0], vec![2.0, 1.5, 1.0, 0.5]];
        let n_k = vec![4usize, 4];
        let f = vec![0.0, 0.5];
        let w = mbar_weights(&u_kn, &n_k, &f);
        assert_eq!(w.len(), 8, "Should have n_total=8 weights");
        let sum: f64 = w.iter().sum();
        assert!(sum > 0.0, "Sum of weights must be positive, got {sum}");
        for &wi in &w {
            assert!(wi >= 0.0, "All weights must be non-negative, got {wi}");
        }
    }
    #[test]
    fn test_mbar_observable_weighted_mean() {
        let weights = vec![1.0, 1.0, 1.0, 1.0];
        let observable = vec![1.0, 2.0, 3.0, 4.0];
        let avg = mbar_observable(&weights, &observable);
        assert!(
            (avg - 2.5).abs() < 1e-12,
            "Uniform weights → mean = 2.5, got {avg}"
        );
    }
    #[test]
    fn test_mbar_observable_single_nonzero_weight() {
        let weights = vec![0.0, 0.0, 1.0, 0.0];
        let observable = vec![5.0, 6.0, 7.0, 8.0];
        let avg = mbar_observable(&weights, &observable);
        assert!(
            (avg - 7.0).abs() < 1e-12,
            "Single nonzero weight → value at that index, got {avg}"
        );
    }
    #[test]
    fn test_remd_free_energy_50pct_acceptance() {
        let alphas = vec![0.5, 0.5, 0.5];
        let df = remd_free_energy_estimates(&alphas, 1.0);
        assert_eq!(df.len(), 3);
        for &v in &df {
            assert!(v.abs() < 1e-8, "50% acceptance → ΔF ≈ 0, got {v}");
        }
    }
    #[test]
    fn test_remd_free_energy_high_acceptance() {
        let alphas = vec![0.8];
        let df = remd_free_energy_estimates(&alphas, 1.0);
        assert!(
            df[0] > 0.0,
            "High acceptance rate → positive ΔF estimate, got {}",
            df[0]
        );
    }
    #[test]
    fn test_remd_free_energy_empty() {
        let df = remd_free_energy_estimates(&[], 1.0);
        assert!(df.is_empty());
    }
    #[test]
    fn test_bootstrap_uncertainty_constant_samples() {
        let samples = vec![2.0_f64; 200];
        let unc = bootstrap_uncertainty(&samples, 50, 1.0, 42);
        assert!(
            unc < 1e-8,
            "Constant samples should give near-zero bootstrap uncertainty, got {unc}"
        );
    }
    #[test]
    fn test_bootstrap_uncertainty_positive() {
        let samples: Vec<f64> = (0..100).map(|i| (i as f64 * 0.31).sin() * 2.0).collect();
        let unc = bootstrap_uncertainty(&samples, 30, 1.0, 17);
        assert!(
            unc >= 0.0,
            "Bootstrap uncertainty must be non-negative, got {unc}"
        );
    }
    #[test]
    fn test_bootstrap_uncertainty_few_samples() {
        let samples = vec![1.0];
        let unc = bootstrap_uncertainty(&samples, 100, 1.0, 0);
        assert_eq!(unc, 0.0, "Single sample → zero uncertainty");
    }
    #[test]
    fn test_ti_gaussian_quadrature_constant() {
        let du_dl = vec![2.0; 5];
        let weights = vec![0.2; 5];
        let df = ti_gaussian_quadrature(&du_dl, &weights);
        assert!((df - 2.0).abs() < 1e-12, "Expected 2.0, got {df}");
    }
    #[test]
    fn test_ti_gaussian_quadrature_empty() {
        let df = ti_gaussian_quadrature(&[], &[]);
        assert_eq!(df, 0.0);
    }
    #[test]
    fn test_gauss_legendre_5pt_weights_sum_to_one() {
        let (_nodes, weights) = gauss_legendre_5pt();
        let sum: f64 = weights.iter().sum();
        assert!(
            sum > 0.0 && sum <= 1.0 + 1e-9,
            "GL5 weights should be in (0, 1], got {sum}"
        );
    }
    #[test]
    fn test_gauss_legendre_5pt_integrates_constant() {
        let (_, weights) = gauss_legendre_5pt();
        let sum_w: f64 = weights.iter().sum();
        let c = 3.0_f64;
        let du_dl = [c; 5];
        let integral: f64 = du_dl
            .iter()
            .zip(weights.iter())
            .map(|(d, w)| d * w)
            .sum::<f64>();
        assert!(
            (integral - sum_w * c).abs() < 1e-10,
            "GL quadrature integral = sum_w * c = {}, got {integral}",
            sum_w * c
        );
    }
    #[test]
    fn test_bar_uncertainty_constant_samples() {
        let fwd = vec![0.0_f64; 50];
        let bwd = vec![0.0_f64; 50];
        let df = 0.0;
        let unc = bar_uncertainty(&fwd, &bwd, df, 1.0);
        assert!(
            unc.is_finite() || unc == f64::INFINITY,
            "BAR uncertainty must be finite or inf"
        );
    }
    #[test]
    fn test_bar_uncertainty_returns_nonnegative() {
        let fwd: Vec<f64> = (0..50).map(|i| i as f64 * 0.1 - 2.5).collect();
        let bwd: Vec<f64> = (0..50).map(|i| -(i as f64 * 0.1 - 2.5)).collect();
        let df = fep_bar(&fwd, &bwd, 1.0, 1e-8);
        let unc = bar_uncertainty(&fwd, &bwd, df, 1.0);
        assert!(
            unc >= 0.0 || unc == f64::INFINITY,
            "BAR uncertainty should be >= 0"
        );
    }
    #[test]
    fn test_ti_block_average_constant() {
        let samples = vec![3.0_f64; 100];
        let (mean, se) = ti_block_average(&samples, 10);
        assert!((mean - 3.0).abs() < 1e-12, "Mean should be 3.0, got {mean}");
        assert!(
            se < 1e-10,
            "Standard error of constant samples should be ~0, got {se}"
        );
    }
    #[test]
    fn test_ti_block_average_empty() {
        let (mean, se) = ti_block_average(&[], 10);
        assert_eq!(mean, 0.0);
        assert_eq!(se, 0.0);
    }
    #[test]
    fn test_ti_block_average_finite_variance() {
        let samples: Vec<f64> = (0..100).map(|i| (i as f64).sin()).collect();
        let (mean, se) = ti_block_average(&samples, 10);
        assert!(mean.is_finite(), "Mean should be finite");
        assert!(se.is_finite(), "SE should be finite");
        assert!(se >= 0.0, "SE should be >= 0");
    }
    #[test]
    fn test_umbrella_bias_energy_at_ref() {
        let e = umbrella_bias_energy(1.5, 1.5, 100.0);
        assert!(
            e.abs() < 1e-15,
            "Bias energy at reference should be 0, got {e}"
        );
    }
    #[test]
    fn test_umbrella_bias_energy_positive_away() {
        let e = umbrella_bias_energy(2.0, 1.0, 10.0);
        assert!(
            e > 0.0,
            "Bias energy away from reference must be positive, got {e}"
        );
    }
    #[test]
    fn test_umbrella_bias_force_at_ref() {
        let f = umbrella_bias_force(1.5, 1.5, 100.0);
        assert!(
            f.abs() < 1e-15,
            "Bias force at reference should be 0, got {f}"
        );
    }
    #[test]
    fn test_umbrella_bias_force_sign() {
        let f = umbrella_bias_force(2.0, 1.0, 10.0);
        assert!(
            f < 0.0,
            "Force should be negative when xi > xi_ref, got {f}"
        );
    }
    #[test]
    fn test_unbias_histogram_normalised() {
        let centers = vec![0.5, 1.5, 2.5];
        let counts = vec![10u64, 20, 30];
        let result = unbias_histogram(&centers, &counts, 1.5, 5.0, 1.0);
        let sum: f64 = result.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "Unbiased histogram should sum to 1, got {sum}"
        );
    }
    #[test]
    fn test_unbias_histogram_empty() {
        let result = unbias_histogram(&[], &[], 0.0, 1.0, 1.0);
        assert!(result.is_empty());
    }
    #[test]
    fn test_jarzynski_from_smd_runs_zero_work() {
        let works = vec![0.0; 100];
        let (df, _se, n) = jarzynski_from_smd_runs(&works, 1.0);
        assert!(df.abs() < 1e-12, "ΔF with zero work should be 0, got {df}");
        assert_eq!(n, 100);
    }
    #[test]
    fn test_jarzynski_from_smd_runs_finite_work() {
        let kt = 1.0;
        let works: Vec<f64> = (1..=10).map(|i| i as f64 * kt).collect();
        let (df, se, n) = jarzynski_from_smd_runs(&works, kt);
        assert!(df.is_finite(), "ΔF should be finite");
        assert!(se.is_finite(), "SE should be finite");
        assert_eq!(n, 10);
        let mean_w = works.iter().sum::<f64>() / works.len() as f64;
        assert!(
            df <= mean_w + 1e-10,
            "Jarzynski ΔF should be ≤ ⟨W⟩ by Jensen, got df={df}, mean_W={mean_w}"
        );
    }
    #[test]
    fn test_jarzynski_from_smd_runs_empty() {
        let (df, se, n) = jarzynski_from_smd_runs(&[], 1.0);
        assert_eq!(df, 0.0);
        assert_eq!(se, 0.0);
        assert_eq!(n, 0);
    }
    #[test]
    fn test_remd_swap_acceptance_equal_temps() {
        let p = remd_swap_acceptance(100.0, 100.0, 300.0, 300.0, 8.314e-3);
        assert!((p - 1.0).abs() < 1e-12, "Same temp → P=1, got {p}");
    }
    #[test]
    fn test_remd_swap_acceptance_favorable() {
        let p = remd_swap_acceptance(200.0, 50.0, 300.0, 600.0, 8.314e-3);
        assert!(
            (p - 1.0).abs() < 1e-6,
            "Favorable swap should give P≈1, got {p}"
        );
    }
    #[test]
    fn test_remd_swap_acceptance_zero_on_bad_params() {
        let p = remd_swap_acceptance(100.0, 50.0, 0.0, 300.0, 1.0);
        assert_eq!(p, 0.0, "Zero temperature should give P=0");
    }
    #[test]
    fn test_remd_temperature_ladder_length() {
        let ladder = remd_temperature_ladder(5, 300.0, 600.0);
        assert_eq!(ladder.len(), 5);
        assert!((ladder[0] - 300.0).abs() < 1e-6);
        assert!((ladder[4] - 600.0).abs() < 1e-6);
    }
    #[test]
    fn test_remd_temperature_ladder_monotone() {
        let ladder = remd_temperature_ladder(8, 280.0, 700.0);
        for i in 0..ladder.len() - 1 {
            assert!(
                ladder[i] < ladder[i + 1],
                "Ladder should be monotonically increasing"
            );
        }
    }
    #[test]
    fn test_remd_expected_swap_probability_finite() {
        let p = remd_expected_swap_probability(300.0, 330.0, 1000.0, 8.314e-3);
        assert!(
            p.is_finite() && (0.0..=1.0).contains(&p),
            "Swap probability should be in [0,1], got {p}"
        );
    }
    #[test]
    fn test_ti_error_estimate_zero_for_zero_errors() {
        let lambdas = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let samples: Vec<Vec<f64>> = lambdas.iter().map(|&l| vec![l, l, l]).collect();
        let ti = ThermoIntegration::from_samples(lambdas, samples);
        let err = ti.compute_error_estimate();
        assert!(
            err.abs() < 1e-10,
            "TI error should be 0 for identical samples, got {err}"
        );
    }
    #[test]
    fn test_ti_error_estimate_positive_for_varied_samples() {
        let lambdas = vec![0.0, 0.5, 1.0];
        let samples = vec![
            vec![1.0, 2.0, 3.0],
            vec![1.5, 2.5, 3.5],
            vec![2.0, 3.0, 4.0],
        ];
        let ti = ThermoIntegration::from_samples(lambdas, samples);
        let err = ti.compute_error_estimate();
        assert!(
            err > 0.0,
            "TI error should be positive for varied samples, got {err}"
        );
    }
    #[test]
    fn test_ti_free_energy_linear_integrand() {
        let lambdas: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();
        let samples: Vec<Vec<f64>> = lambdas.iter().map(|&l| vec![2.0 * l; 100]).collect();
        let ti = ThermoIntegration::from_samples(lambdas, samples);
        let df = ti.compute_free_energy();
        assert!(
            (df - 1.0).abs() < 1e-10,
            "linear TI: ΔF should be 1.0, got {df}"
        );
    }
    #[test]
    fn test_mbar_uncertainty_zero_bootstrap() {
        let u_kn = vec![vec![1.0, 2.0], vec![2.0, 1.0]];
        let n_k = vec![1, 1];
        let mbar = Mbar::new(u_kn, n_k, 1e-6, 100);
        let unc = mbar.compute_uncertainty(0, 0);
        assert!(
            unc.iter().all(|&v| v == 0.0),
            "zero bootstrap → zero uncertainty"
        );
    }
    #[test]
    fn test_mbar_uncertainty_non_negative() {
        let u_kn = vec![vec![0.0, 0.5, 1.0, 1.5, 2.0], vec![2.0, 1.5, 1.0, 0.5, 0.0]];
        let n_k = vec![3, 2];
        let mbar = Mbar::new(u_kn, n_k, 1e-4, 200);
        let unc = mbar.compute_uncertainty(5, 42);
        assert!(
            unc.iter().all(|&v| v >= 0.0),
            "bootstrap uncertainties must be non-negative"
        );
    }
    #[test]
    fn test_mbar_free_energies_normalised() {
        let u_kn = vec![vec![1.0, 2.0, 3.0], vec![4.0, 3.0, 2.0]];
        let n_k = vec![2, 1];
        let mbar = Mbar::new(u_kn, n_k, 1e-8, 1000);
        assert!(
            (mbar.free_energies[0] - 0.0).abs() < 1e-6,
            "f[0] should be 0 (normalised)"
        );
    }
    #[test]
    fn test_pmf_barrier_simple() {
        let xi: Vec<f64> = (0..=10).map(|i| i as f64).collect();
        let pmf: Vec<f64> = xi.iter().map(|&x| -(x - 5.0).powi(2) + 25.0).collect();
        let profile = FreeEnergyProfile::new(xi, pmf);
        let barrier = profile.compute_barrier(0.0, 10.0).unwrap();
        assert!(
            (barrier - 25.0).abs() < 1e-10,
            "barrier should be 25, got {barrier}"
        );
    }
    #[test]
    fn test_pmf_barrier_empty_range() {
        let xi = vec![0.0, 1.0, 2.0];
        let pmf = vec![0.0, 1.0, 0.5];
        let profile = FreeEnergyProfile::new(xi, pmf);
        let barrier = profile.compute_barrier(5.0, 10.0);
        assert!(barrier.is_none(), "barrier should be None for empty range");
    }
    #[test]
    fn test_pmf_transition_state_position() {
        let xi = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let pmf = vec![0.0, 2.0, 5.0, 2.0, 0.0];
        let profile = FreeEnergyProfile::new(xi, pmf);
        let ts = profile.transition_state().unwrap();
        assert!(
            (ts - 2.0).abs() < 1e-10,
            "transition state should be at xi=2, got {ts}"
        );
    }
}
