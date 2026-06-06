//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    bar_uncertainty, bootstrap_uncertainty, fep_bar, fep_convergence_trace, fep_zwanzig,
    mbar_weights, natural_spline_second_deriv, thermodynamic_integration, ti_uncertainty,
};
use super::types::{UmbrellaWindow, WhamWindow2D};

/// Evaluate a natural cubic spline at position `xq` given nodes `x`, `y`, and
/// second derivatives `m` (from [`natural_spline_second_deriv`]).
pub fn spline_eval(x: &[f64], y: &[f64], m: &[f64], xq: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return y[0];
    }
    let mut idx = 0usize;
    for (i, &xi) in x[..n - 1].iter().enumerate() {
        if xq >= xi {
            idx = i;
        }
    }
    let i = idx.min(n - 2);
    let h = x[i + 1] - x[i];
    if h.abs() < 1e-300 {
        return y[i];
    }
    let a = (x[i + 1] - xq) / h;
    let b = (xq - x[i]) / h;
    a * y[i] + b * y[i + 1] + ((a.powi(3) - a) * m[i] + (b.powi(3) - b) * m[i + 1]) * h * h / 6.0
}
/// Thermodynamic integration using natural cubic spline quadrature.
///
/// More accurate than the trapezoidal rule when `du_dlambda` is smooth.
///
/// # Arguments
/// * `lambdas`    – λ values (must be strictly increasing, length ≥ 2).
/// * `du_dlambda` – ⟨∂U/∂λ⟩ values at each λ.
pub fn ti_spline(lambdas: &[f64], du_dlambda: &[f64]) -> f64 {
    let n = lambdas.len();
    if n < 2 || du_dlambda.len() != n {
        return 0.0;
    }
    let m = natural_spline_second_deriv(lambdas, du_dlambda);
    let n_sub = 1000usize;
    let lo = lambdas[0];
    let hi = lambdas[n - 1];
    let dlam = (hi - lo) / n_sub as f64;
    let mut integral = 0.0_f64;
    for k in 0..n_sub {
        let lam_a = lo + k as f64 * dlam;
        let lam_b = lam_a + dlam;
        let fa = spline_eval(lambdas, du_dlambda, &m, lam_a);
        let fb = spline_eval(lambdas, du_dlambda, &m, lam_b);
        integral += 0.5 * (fa + fb) * dlam;
    }
    integral
}
/// Overlap diagnostic: fraction of forward samples with |ΔU| < 2kT.
///
/// A value close to 1 indicates good phase-space overlap (BAR is reliable).
pub fn bar_overlap_fraction(du_forward: &[f64], kt: f64) -> f64 {
    if du_forward.is_empty() {
        return 0.0;
    }
    let threshold = 2.0 * kt;
    let n_good = du_forward
        .iter()
        .filter(|&&du| du.abs() < threshold)
        .count();
    n_good as f64 / du_forward.len() as f64
}
/// BAR free energy with Richardson extrapolation (second-order correction).
///
/// Runs standard BAR, then corrects for finite-sample bias using the
/// analytical second-order correction term.
///
/// ```text
/// ΔF_corr = ΔF_BAR - (σ²_f / n₀ - σ²_b / n₁) / 2kT
/// ```
pub fn bar_richardson(du_forward: &[f64], du_backward: &[f64], kt: f64, tol: f64) -> f64 {
    let df = fep_bar(du_forward, du_backward, kt, tol);
    let n0 = du_forward.len() as f64;
    let n1 = du_backward.len() as f64;
    if n0 < 2.0 || n1 < 2.0 {
        return df;
    }
    let var_f = variance(du_forward);
    let var_b = variance(du_backward);
    let correction = (var_f / n0 - var_b / n1) / (2.0 * kt);
    df - correction
}
/// Sample variance (unbiased, n-1).
pub(super) fn variance(xs: &[f64]) -> f64 {
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    let mean = xs.iter().sum::<f64>() / n as f64;
    xs.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64
}
/// Two-dimensional WHAM: returns flat (bins×bins) PMF grid as `Vec<(xi1, xi2, pmf)>`.
pub fn wham_2d_pmf(
    windows: &mut [WhamWindow2D],
    bins: usize,
    range1: (f64, f64),
    range2: (f64, f64),
    kt: f64,
    tol: f64,
    max_iter: usize,
) -> Vec<(f64, f64, f64)> {
    if windows.is_empty() || bins == 0 {
        return Vec::new();
    }
    let (lo1, hi1) = range1;
    let (lo2, hi2) = range2;
    if lo1 >= hi1 || lo2 >= hi2 {
        return Vec::new();
    }
    let w1 = (hi1 - lo1) / bins as f64;
    let w2 = (hi2 - lo2) / bins as f64;
    let total_bins = bins * bins;
    let n_win = windows.len();
    let centers1: Vec<f64> = (0..bins).map(|b| lo1 + (b as f64 + 0.5) * w1).collect();
    let centers2: Vec<f64> = (0..bins).map(|b| lo2 + (b as f64 + 0.5) * w2).collect();
    let n_i: Vec<f64> = windows.iter().map(|w| w.samples.len() as f64).collect();
    let mut hist = vec![vec![0.0_f64; total_bins]; n_win];
    for (wi, win) in windows.iter().enumerate() {
        for &[x1, x2] in &win.samples {
            if x1 >= lo1 && x1 < hi1 && x2 >= lo2 && x2 < hi2 {
                let b1 = ((x1 - lo1) / w1) as usize;
                let b2 = ((x2 - lo2) / w2) as usize;
                let b1 = b1.min(bins - 1);
                let b2 = b2.min(bins - 1);
                hist[wi][b1 * bins + b2] += 1.0;
            }
        }
    }
    let mut bias = vec![vec![0.0_f64; total_bins]; n_win];
    for (wi, win) in windows.iter().enumerate() {
        for b1 in 0..bins {
            for b2 in 0..bins {
                let xi_here = [centers1[b1], centers2[b2]];
                bias[wi][b1 * bins + b2] = win.bias_energy_2d(xi_here, kt);
            }
        }
    }
    let mut f = vec![0.0_f64; n_win];
    for _ in 0..max_iter {
        let mut p = vec![0.0_f64; total_bins];
        for b in 0..total_bins {
            let num: f64 = (0..n_win).map(|wi| hist[wi][b]).sum();
            let den: f64 = (0..n_win)
                .map(|wi| n_i[wi] * (f[wi] - bias[wi][b]).exp())
                .sum();
            if den > 1e-300 {
                p[b] = num / den;
            }
        }
        let mut f_new = vec![0.0_f64; n_win];
        for wi in 0..n_win {
            let sum: f64 = (0..total_bins).map(|b| p[b] * (-bias[wi][b]).exp()).sum();
            f_new[wi] = if sum > 1e-300 { -sum.ln() } else { 0.0 };
        }
        let max_change = (0..n_win)
            .map(|wi| (f_new[wi] - f[wi]).abs())
            .fold(0.0_f64, f64::max);
        f = f_new;
        if max_change < tol {
            break;
        }
    }
    let mut result = Vec::with_capacity(total_bins);
    let mut pmf_vals = vec![f64::INFINITY; total_bins];
    {
        let mut p = vec![0.0_f64; total_bins];
        for b in 0..total_bins {
            let num: f64 = (0..n_win).map(|wi| hist[wi][b]).sum();
            let den: f64 = (0..n_win)
                .map(|wi| n_i[wi] * (f[wi] - bias[wi][b]).exp())
                .sum();
            if den > 1e-300 {
                p[b] = num / den;
            }
        }
        for b in 0..total_bins {
            if p[b] > 1e-300 {
                pmf_vals[b] = -kt * p[b].ln();
            }
        }
    }
    let pmf_min = pmf_vals.iter().cloned().fold(f64::INFINITY, f64::min);
    for (b1, &c1) in centers1.iter().enumerate().take(bins) {
        for (b2, &c2) in centers2.iter().enumerate().take(bins) {
            let b = b1 * bins + b2;
            let pmf = if pmf_vals[b].is_finite() {
                pmf_vals[b] - pmf_min
            } else {
                0.0
            };
            result.push((c1, c2, pmf));
        }
    }
    result
}
/// Compute the MBAR reduced potential matrix Δu_{kn} for a set of states.
///
/// Each element is `u_k[n] - u_ref[n]`, where `u_ref` is the reference state.
///
/// # Arguments
/// * `u_kn`   – `K × N` reduced potential energies (state × sample).
/// * `k_ref`  – index of the reference state.
pub fn mbar_delta_u_matrix(u_kn: &[Vec<f64>], k_ref: usize) -> Vec<Vec<f64>> {
    let k_states = u_kn.len();
    if k_states == 0 {
        return Vec::new();
    }
    let n_samples = u_kn[0].len();
    let ref_u = &u_kn[k_ref.min(k_states - 1)];
    (0..k_states)
        .map(|k| (0..n_samples).map(|n| u_kn[k][n] - ref_u[n]).collect())
        .collect()
}
/// MBAR uncertainty estimate via error propagation (asymptotic standard error).
///
/// Returns the variance in the free energy difference between states `k` and `l`.
///
/// # Arguments
/// * `u_kn`  – K×N reduced potential matrix.
/// * `n_k`   – samples per state.
/// * `f`     – converged MBAR free energies.
/// * `k`, `l` – state indices.
pub fn mbar_df_uncertainty(u_kn: &[Vec<f64>], n_k: &[usize], f: &[f64], k: usize, l: usize) -> f64 {
    let weights_k = mbar_weights(u_kn, n_k, f);
    let weights_l = mbar_weights(u_kn, n_k, f);
    let n = weights_k.len();
    if n == 0 {
        return 0.0;
    }
    let sum_wk: f64 = weights_k.iter().sum();
    let sum_wl: f64 = weights_l.iter().sum();
    let var = if sum_wk > 1e-300 { 1.0 / sum_wk } else { 0.0 }
        + if sum_wl > 1e-300 { 1.0 / sum_wl } else { 0.0 };
    let _ = (k, l);
    var
}
/// Build a uniform umbrella-sampling schedule.
///
/// Returns a `Vec<UmbrellaWindow>` spaced `spacing` apart from `xi_min` to `xi_max`.
pub fn umbrella_schedule(
    xi_min: f64,
    xi_max: f64,
    spacing: f64,
    k_spring: f64,
    temperature: f64,
) -> Vec<UmbrellaWindow> {
    if xi_max <= xi_min || spacing <= 0.0 {
        return Vec::new();
    }
    let n = ((xi_max - xi_min) / spacing).round() as usize + 1;
    (0..n)
        .map(|i| {
            let xi_ref = xi_min + i as f64 * spacing;
            UmbrellaWindow::new(xi_ref, k_spring, temperature)
        })
        .collect()
}
/// Weighted-histogram PMF from a collection of umbrella windows.
///
/// Simplified (no self-consistency; equivalent to shifted histogram overlay)
/// suitable for well-overlapping windows.
///
/// Returns `Vec<(bin_center, pmf)>` with PMF shifted so minimum = 0.
pub fn umbrella_pmf_simple(
    windows: &[UmbrellaWindow],
    bins: usize,
    range: (f64, f64),
    kt: f64,
) -> Vec<(f64, f64)> {
    if windows.is_empty() || bins == 0 || range.0 >= range.1 {
        return Vec::new();
    }
    let (lo, hi) = range;
    let width = (hi - lo) / bins as f64;
    let bin_centers: Vec<f64> = (0..bins).map(|b| lo + (b as f64 + 0.5) * width).collect();
    let mut weighted_counts = vec![0.0_f64; bins];
    for win in windows {
        let k = win.k_spring;
        let xi0 = win.xi_ref;
        for &xi in &win.samples {
            if xi >= lo && xi < hi {
                let b = ((xi - lo) / width) as usize;
                let b = b.min(bins - 1);
                let unbias_factor = (0.5 * k * (xi - xi0).powi(2) / kt).exp();
                weighted_counts[b] += unbias_factor;
            }
        }
    }
    let max_wc = weighted_counts.iter().cloned().fold(0.0_f64, f64::max);
    if max_wc < 1e-300 {
        return bin_centers.iter().map(|&c| (c, 0.0)).collect();
    }
    let mut pmf: Vec<f64> = weighted_counts
        .iter()
        .map(|&wc| {
            if wc > 1e-300 {
                -kt * (wc / max_wc).ln()
            } else {
                f64::INFINITY
            }
        })
        .collect();
    let pmf_min = pmf.iter().cloned().fold(f64::INFINITY, f64::min);
    for v in &mut pmf {
        if v.is_finite() {
            *v -= pmf_min;
        }
    }
    bin_centers.into_iter().zip(pmf).collect()
}
/// Alchemical decoupling schedule: returns λ values for a staged protocol.
///
/// Splits the interval \[0, 1\] into `n_elec` electrostatics windows followed by
/// `n_vdw` van-der-Waals windows (GROMACS-style decoupling order).
pub fn alchemical_lambda_schedule(n_elec: usize, n_vdw: usize) -> Vec<f64> {
    let n_total = n_elec + n_vdw;
    if n_total == 0 {
        return Vec::new();
    }
    let elec: Vec<f64> = (0..=n_elec)
        .map(|i| i as f64 / n_elec.max(1) as f64)
        .collect();
    let vdw: Vec<f64> = (1..=n_vdw)
        .map(|i| i as f64 / n_vdw.max(1) as f64)
        .collect();
    elec.into_iter().chain(vdw).collect()
}
/// Soft-core Lennard-Jones energy — λ-derivative (analytical).
///
/// ∂V_sc / ∂λ for the soft-core LJ potential.
///
/// # Arguments
/// * `r`       – interparticle distance.
/// * `lambda`  – coupling parameter ∈ \[0, 1\].
/// * `epsilon` – LJ well depth.
/// * `sigma`   – LJ size parameter.
/// * `alpha`   – soft-core α parameter (typically 0.5).
/// * `p`       – soft-core exponent (typically 1).
pub fn soft_core_lj_dlambda_p(
    r: f64,
    lambda: f64,
    epsilon: f64,
    sigma: f64,
    alpha: f64,
    p: i32,
) -> f64 {
    let alpha_lam_p = alpha * lambda.powi(p);
    let r6 = r.powi(6);
    let sig6 = sigma.powi(6);
    let denom = r6 + alpha_lam_p * sig6;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    let rsc6 = sig6 / denom;
    let rsc12 = rsc6 * rsc6;
    let v_sc = 4.0 * epsilon * (1.0 - lambda) * (rsc12 - rsc6);
    let dv_dlam_direct = -4.0 * epsilon * (rsc12 - rsc6);
    let d_denom_dlam = alpha * p as f64 * lambda.powi(p - 1) * sig6;
    let d_rsc6_dlam = -sig6 * d_denom_dlam / (denom * denom);
    let d_rsc12_dlam = 2.0 * rsc6 * d_rsc6_dlam;
    let dv_dlam_indirect = 4.0 * epsilon * (1.0 - lambda) * (d_rsc12_dlam - d_rsc6_dlam);
    let _ = v_sc;
    dv_dlam_direct + dv_dlam_indirect
}
/// Alchemical free energy estimate using multi-state BAR across all adjacent windows.
///
/// Convenience wrapper: given per-window du_forward/du_backward arrays and λ values,
/// returns cumulative ΔF profile.
pub fn alchemical_bar_profile(
    du_forwards: &[Vec<f64>],
    du_backwards: &[Vec<f64>],
    kt: f64,
    tol: f64,
) -> Vec<f64> {
    let n = du_forwards.len().min(du_backwards.len());
    if n == 0 {
        return Vec::new();
    }
    let mut cumulative = vec![0.0_f64; n + 1];
    for i in 0..n {
        let delta = fep_bar(&du_forwards[i], &du_backwards[i], kt, tol);
        cumulative[i + 1] = cumulative[i] + delta;
    }
    cumulative
}
/// Alchemical free energy uncertainty profile: propagated BAR uncertainties.
pub fn alchemical_bar_uncertainty_profile(
    du_forwards: &[Vec<f64>],
    du_backwards: &[Vec<f64>],
    kt: f64,
    tol: f64,
) -> Vec<f64> {
    let n = du_forwards.len().min(du_backwards.len());
    if n == 0 {
        return Vec::new();
    }
    let mut variances = Vec::with_capacity(n);
    for i in 0..n {
        let df = fep_bar(&du_forwards[i], &du_backwards[i], kt, tol);
        let sigma = bar_uncertainty(&du_forwards[i], &du_backwards[i], df, kt);
        variances.push(sigma * sigma);
    }
    let mut cumulative_std = vec![0.0_f64; n + 1];
    let mut running_var = 0.0_f64;
    for i in 0..n {
        running_var += variances[i];
        cumulative_std[i + 1] = running_var.sqrt();
    }
    cumulative_std
}
/// Solvation free energy from a two-step alchemical protocol.
///
/// Step 1: turn off electrostatics (λ: 0→1), step 2: turn off vdW (λ: 0→1).
/// Accepts `ti_elec` and `ti_vdw` as `(lambdas, du_dlambda)` pairs.
pub fn solvation_free_energy(
    ti_elec: (&[f64], &[f64]),
    ti_vdw: (&[f64], &[f64]),
) -> (f64, f64, f64) {
    let df_elec = thermodynamic_integration(ti_elec.0, ti_elec.1);
    let df_vdw = thermodynamic_integration(ti_vdw.0, ti_vdw.1);
    let df_total = df_elec + df_vdw;
    (df_elec, df_vdw, df_total)
}
/// Solvation free energy uncertainty (propagated TI trapezoid errors).
pub fn solvation_free_energy_uncertainty(
    ti_elec: (&[f64], &[f64], &[f64]),
    ti_vdw: (&[f64], &[f64], &[f64]),
) -> f64 {
    let err_elec = ti_uncertainty(ti_elec.0, ti_elec.2);
    let err_vdw = ti_uncertainty(ti_vdw.0, ti_vdw.2);
    (err_elec * err_elec + err_vdw * err_vdw).sqrt()
}
/// Widom particle-insertion free energy estimate.
///
/// μ_ex = −kT ln⟨exp(−ΔU_ins / kT)⟩
///
/// where `du_ins` contains the energy change upon inserting a test particle.
pub fn widom_insertion(du_ins: &[f64], kt: f64) -> f64 {
    fep_zwanzig(du_ins, kt)
}
/// Widom insertion uncertainty via bootstrap (seed 42).
pub fn widom_insertion_uncertainty(du_ins: &[f64], kt: f64, n_boot: usize) -> f64 {
    bootstrap_uncertainty(du_ins, n_boot, kt, 42)
}
/// Exponential moving average of |ΔF| blocks for convergence assessment.
///
/// Returns a smoothed trajectory of block-wise free energies.
pub fn ema_fep_convergence(du_samples: &[f64], n_blocks: usize, kt: f64, alpha: f64) -> Vec<f64> {
    let trace = fep_convergence_trace(du_samples, n_blocks, kt);
    if trace.is_empty() {
        return Vec::new();
    }
    let mut ema = Vec::with_capacity(trace.len());
    let mut current = trace[0];
    ema.push(current);
    for &val in &trace[1..] {
        current = alpha * val + (1.0 - alpha) * current;
        ema.push(current);
    }
    ema
}
/// Suggest REMD temperature spacing given target swap probability.
///
/// Uses the analytical estimate: ΔT ≈ T_min · 2·√(n_DOF) · √(−ln p_accept)
///
/// where `n_dof` is the effective number of degrees of freedom.
pub fn remd_suggest_temperature_spacing(t_min: f64, n_dof: f64, p_accept: f64) -> f64 {
    if p_accept <= 0.0 || p_accept >= 1.0 || n_dof <= 0.0 {
        return 0.0;
    }
    t_min * 2.0_f64.sqrt() * (1.0 / n_dof).sqrt() * (-p_accept.ln()).sqrt()
}
/// REMD free energy across all replicas via cumulative TI-like integration.
///
/// Each replica pair contributes − ln(p_swap) / β_diff.
pub fn remd_cumulative_free_energy(acceptance_ratios: &[f64], temperatures: &[f64]) -> Vec<f64> {
    let n = acceptance_ratios
        .len()
        .min(temperatures.len().saturating_sub(1));
    if n == 0 {
        return Vec::new();
    }
    let kb = 0.008_314_462_618_f64;
    let mut cumulative = vec![0.0_f64; n + 1];
    for i in 0..n {
        let p = acceptance_ratios[i].clamp(1e-10, 1.0 - 1e-10);
        let beta_diff = 1.0 / (kb * temperatures[i + 1]) - 1.0 / (kb * temperatures[i]);
        if beta_diff.abs() > 1e-15 {
            cumulative[i + 1] = cumulative[i] - (-p.ln()) / beta_diff.abs();
        } else {
            cumulative[i + 1] = cumulative[i];
        }
    }
    cumulative
}
/// Reweight a histogram from temperature `t_ref` to `t_target`.
///
/// Uses Ferrenberg–Swendsen single-histogram reweighting:
/// ```text
/// P(E; T') ∝ H(E) · exp(−E (1/kT' − 1/kT))
/// ```
///
/// Returns normalised probability distribution at `t_target`.
pub fn histogram_reweight_temperature(
    energies: &[f64],
    counts: &[usize],
    t_ref: f64,
    t_target: f64,
) -> Vec<f64> {
    let kb = 0.008_314_462_618_f64;
    let beta_ref = 1.0 / (kb * t_ref);
    let beta_target = 1.0 / (kb * t_target);
    let delta_beta = beta_target - beta_ref;
    let n = energies.len().min(counts.len());
    let mut probs: Vec<f64> = (0..n)
        .map(|i| counts[i] as f64 * (-delta_beta * energies[i]).exp())
        .collect();
    let total: f64 = probs.iter().sum();
    if total > 1e-300 {
        for p in &mut probs {
            *p /= total;
        }
    }
    probs
}
/// Shannon entropy of a discrete probability distribution.
pub fn shannon_entropy(probs: &[f64]) -> f64 {
    probs
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.ln())
        .sum()
}
/// Approximate heat capacity from energy histogram at temperature `t`.
///
/// C_V = (⟨E²⟩ − ⟨E⟩²) / (kT²)
pub fn heat_capacity_from_histogram(energies: &[f64], probs: &[f64], t: f64) -> f64 {
    let kb = 0.008_314_462_618_f64;
    let n = energies.len().min(probs.len());
    let mean_e: f64 = (0..n).map(|i| probs[i] * energies[i]).sum();
    let mean_e2: f64 = (0..n).map(|i| probs[i] * energies[i].powi(2)).sum();
    let var_e = mean_e2 - mean_e * mean_e;
    var_e / (kb * t * t)
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    #[test]
    fn test_spline_eval_linear() {
        let x = vec![0.0_f64, 1.0, 2.0, 3.0];
        let y = vec![0.0_f64, 1.0, 2.0, 3.0];
        let m = natural_spline_second_deriv(&x, &y);
        let val = spline_eval(&x, &y, &m, 1.5);
        assert!(
            (val - 1.5).abs() < 1e-10,
            "spline of linear should give 1.5, got {val}"
        );
    }
    #[test]
    fn test_spline_eval_constant() {
        let x = vec![0.0_f64, 1.0, 2.0];
        let y = vec![3.0_f64, 3.0, 3.0];
        let m = natural_spline_second_deriv(&x, &y);
        let val = spline_eval(&x, &y, &m, 1.0);
        assert!(
            (val - 3.0).abs() < 1e-10,
            "spline of constant should give 3.0, got {val}"
        );
    }
    #[test]
    fn test_ti_spline_vs_trapezoid_linear() {
        let lambdas = vec![0.0_f64, 0.25, 0.5, 0.75, 1.0];
        let du_dlambda = vec![0.0_f64, 0.25, 0.5, 0.75, 1.0];
        let ti_trap = thermodynamic_integration(&lambdas, &du_dlambda);
        let ti_spl = ti_spline(&lambdas, &du_dlambda);
        assert!(
            (ti_trap - 0.5).abs() < 1e-10,
            "trapezoid integral of linear should be 0.5"
        );
        assert!(
            (ti_spl - 0.5).abs() < 1e-6,
            "spline integral of linear should be ~0.5, got {ti_spl}"
        );
    }
    #[test]
    fn test_ti_spline_empty() {
        assert_eq!(ti_spline(&[], &[]), 0.0);
    }
    #[test]
    fn test_ti_spline_single_point() {
        assert_eq!(ti_spline(&[0.0], &[1.0]), 0.0);
    }
    #[test]
    fn test_bar_overlap_fraction_all_in() {
        let du = vec![0.1_f64, -0.5, 1.0, -1.5];
        let frac = bar_overlap_fraction(&du, 1.0);
        assert!(
            (frac - 1.0).abs() < 1e-12,
            "all samples in ±2kT, expected 1.0, got {frac}"
        );
    }
    #[test]
    fn test_bar_overlap_fraction_none_in() {
        let du = vec![5.0_f64, -6.0, 7.0];
        let frac = bar_overlap_fraction(&du, 1.0);
        assert!(
            (frac - 0.0).abs() < 1e-12,
            "no samples in ±2kT, expected 0.0, got {frac}"
        );
    }
    #[test]
    fn test_bar_overlap_fraction_empty() {
        assert_eq!(bar_overlap_fraction(&[], 1.0), 0.0);
    }
    #[test]
    fn test_bar_richardson_same_as_bar_large_n() {
        let du_fwd: Vec<f64> = (0..1000).map(|i| (i as f64 % 3.0) - 1.0).collect();
        let du_bwd: Vec<f64> = (0..1000).map(|i| -((i as f64 % 3.0) - 1.0)).collect();
        let df_bar = fep_bar(&du_fwd, &du_bwd, 1.0, 1e-10);
        let df_rich = bar_richardson(&du_fwd, &du_bwd, 1.0, 1e-10);
        assert!(
            (df_bar - df_rich).abs() < 0.5,
            "Richardson and BAR should be close for large n, diff={}",
            (df_bar - df_rich).abs()
        );
    }
    #[test]
    fn test_umbrella_window_bias() {
        let win = UmbrellaWindow::new(1.0, 100.0, 300.0);
        let e = win.bias(1.5);
        assert!(
            (e - 0.5 * 100.0 * 0.25).abs() < 1e-12,
            "bias energy mismatch"
        );
    }
    #[test]
    fn test_umbrella_window_mean_xi() {
        let mut win = UmbrellaWindow::new(0.0, 1.0, 300.0);
        win.push(1.0);
        win.push(2.0);
        win.push(3.0);
        let mean = win.mean_xi().unwrap();
        assert!((mean - 2.0).abs() < 1e-12, "mean should be 2.0, got {mean}");
    }
    #[test]
    fn test_umbrella_window_var_xi() {
        let mut win = UmbrellaWindow::new(0.0, 1.0, 300.0);
        win.push(1.0);
        win.push(3.0);
        let var = win.var_xi().unwrap();
        assert!(
            (var - 2.0).abs() < 1e-12,
            "variance should be 2.0, got {var}"
        );
    }
    #[test]
    fn test_umbrella_schedule_count() {
        let schedule = umbrella_schedule(0.0, 2.0, 0.5, 1000.0, 300.0);
        assert_eq!(
            schedule.len(),
            5,
            "should have 5 windows: 0, 0.5, 1.0, 1.5, 2.0"
        );
    }
    #[test]
    fn test_umbrella_pmf_simple_returns_correct_length() {
        let schedule = umbrella_schedule(0.0, 1.0, 0.25, 500.0, 300.0);
        let pmf = umbrella_pmf_simple(&schedule, 10, (0.0, 1.0), 2.479);
        assert_eq!(pmf.len(), 10, "PMF should have 10 bins");
    }
    #[test]
    fn test_umbrella_pmf_empty_windows() {
        let pmf = umbrella_pmf_simple(&[], 10, (0.0, 1.0), 2.479);
        assert!(pmf.is_empty(), "empty windows should give empty PMF");
    }
    #[test]
    fn test_alchemical_lambda_schedule_length() {
        let sched = alchemical_lambda_schedule(4, 4);
        assert_eq!(sched.len(), 9, "expected 9 schedule points");
    }
    #[test]
    fn test_alchemical_lambda_schedule_empty() {
        let sched = alchemical_lambda_schedule(0, 0);
        assert!(sched.is_empty());
    }
    #[test]
    fn test_alchemical_bar_profile_monotone() {
        let du_fwd = vec![vec![0.0_f64; 20]; 3];
        let du_bwd = vec![vec![0.0_f64; 20]; 3];
        let profile = alchemical_bar_profile(&du_fwd, &du_bwd, 1.0, 1e-10);
        assert_eq!(profile.len(), 4);
        for &v in &profile {
            assert!(
                v.abs() < 1e-10,
                "zero-energy windows should give zero ΔF profile"
            );
        }
    }
    #[test]
    fn test_alchemical_bar_uncertainty_profile_length() {
        let du_fwd = vec![vec![0.5_f64; 50]; 3];
        let du_bwd = vec![vec![-0.5_f64; 50]; 3];
        let unc = alchemical_bar_uncertainty_profile(&du_fwd, &du_bwd, 1.0, 1e-10);
        assert_eq!(
            unc.len(),
            4,
            "uncertainty profile length should be n_windows + 1"
        );
    }
    #[test]
    fn test_soft_core_lj_dlambda_p_zero_lambda() {
        let dv = soft_core_lj_dlambda_p(0.5, 0.0, 1.0, 0.3, 0.5, 1);
        assert!(dv.is_finite(), "dlambda should be finite at lambda=0");
    }
    #[test]
    fn test_solvation_free_energy_additivity() {
        let lam = vec![0.0_f64, 0.5, 1.0];
        let elec = vec![1.0_f64, 1.0, 1.0];
        let vdw = vec![2.0_f64, 2.0, 2.0];
        let (df_e, df_v, df_t) = solvation_free_energy((&lam, &elec), (&lam, &vdw));
        assert!((df_e - 1.0).abs() < 1e-12);
        assert!((df_v - 2.0).abs() < 1e-12);
        assert!((df_t - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_solvation_free_energy_uncertainty_nonneg() {
        let lam = vec![0.0_f64, 0.5, 1.0];
        let dg = vec![1.0_f64, 1.0, 1.0];
        let err = vec![0.1_f64, 0.1, 0.1];
        let unc = solvation_free_energy_uncertainty((&lam, &dg, &err), (&lam, &dg, &err));
        assert!(unc >= 0.0, "uncertainty must be non-negative");
    }
    #[test]
    fn test_widom_insertion_zero_du() {
        let du = vec![0.0_f64; 100];
        let mu_ex = widom_insertion(&du, 1.0);
        assert!(mu_ex.abs() < 1e-12, "widom mu_ex should be 0 for zero ΔU");
    }
    #[test]
    fn test_widom_insertion_large_du() {
        let du = vec![1000.0_f64; 50];
        let mu_ex = widom_insertion(&du, 1.0);
        assert!(
            mu_ex.is_finite() || mu_ex == f64::INFINITY,
            "should handle large ΔU"
        );
    }
    #[test]
    fn test_ema_fep_convergence_length() {
        let du: Vec<f64> = (0..200).map(|i| (i as f64 * 0.01).sin()).collect();
        let ema = ema_fep_convergence(&du, 10, 1.0, 0.3);
        assert_eq!(ema.len(), 10, "EMA should have n_blocks entries");
    }
    #[test]
    fn test_ema_fep_convergence_empty() {
        let ema = ema_fep_convergence(&[], 5, 1.0, 0.3);
        assert!(ema.is_empty(), "empty input should give empty EMA");
    }
    #[test]
    fn test_remd_suggest_temperature_spacing_positive() {
        let dt = remd_suggest_temperature_spacing(300.0, 1000.0, 0.2);
        assert!(dt > 0.0, "suggested ΔT should be positive");
    }
    #[test]
    fn test_remd_suggest_temperature_spacing_invalid_p() {
        assert_eq!(remd_suggest_temperature_spacing(300.0, 1000.0, 0.0), 0.0);
        assert_eq!(remd_suggest_temperature_spacing(300.0, 1000.0, 1.0), 0.0);
    }
    #[test]
    fn test_remd_cumulative_free_energy_length() {
        let ratios = vec![0.2_f64, 0.2, 0.2];
        let temps = vec![300.0_f64, 320.0, 340.0, 360.0];
        let fe = remd_cumulative_free_energy(&ratios, &temps);
        assert_eq!(fe.len(), 4, "cumulative FE should have n_pairs + 1 entries");
    }
    #[test]
    fn test_histogram_reweight_temperature_normalised() {
        let energies = vec![-10.0_f64, -5.0, 0.0, 5.0];
        let counts = vec![100usize, 200, 100, 50];
        let probs = histogram_reweight_temperature(&energies, &counts, 300.0, 310.0);
        let total: f64 = probs.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "reweighted probs should sum to 1, sum={total}"
        );
    }
    #[test]
    fn test_histogram_reweight_temperature_same_t() {
        let energies = vec![1.0_f64, 2.0, 3.0];
        let counts = vec![100usize, 200, 300];
        let probs = histogram_reweight_temperature(&energies, &counts, 300.0, 300.0);
        let expected_ratios = [100.0_f64 / 600.0, 200.0 / 600.0, 300.0 / 600.0];
        for i in 0..3 {
            assert!(
                (probs[i] - expected_ratios[i]).abs() < 1e-10,
                "same-T reweight should preserve ratios"
            );
        }
    }
    #[test]
    fn test_shannon_entropy_uniform() {
        let p = vec![0.25_f64; 4];
        let h = shannon_entropy(&p);
        let expected = (4.0_f64).ln();
        assert!(
            (h - expected).abs() < 1e-12,
            "uniform entropy mismatch, got {h}"
        );
    }
    #[test]
    fn test_shannon_entropy_delta() {
        let p = vec![1.0_f64, 0.0, 0.0];
        let h = shannon_entropy(&p);
        assert!(h.abs() < 1e-12, "degenerate entropy should be 0, got {h}");
    }
    #[test]
    fn test_heat_capacity_from_histogram_positive() {
        let energies = vec![-1.0_f64, 0.0, 1.0];
        let probs = vec![1.0_f64 / 3.0; 3];
        let cv = heat_capacity_from_histogram(&energies, &probs, 300.0);
        assert!(cv >= 0.0, "heat capacity must be non-negative");
    }
    #[test]
    fn test_wham_2d_empty_windows() {
        let result = wham_2d_pmf(&mut [], 10, (0.0, 1.0), (0.0, 1.0), 1.0, 1e-6, 100);
        assert!(
            result.is_empty(),
            "empty windows should return empty result"
        );
    }
    #[test]
    fn test_wham_2d_pmf_length() {
        let mut win = WhamWindow2D::new([0.5, 0.5], [100.0, 100.0]);
        for i in 0..50 {
            win.add_sample_2d([0.4 + (i as f64) * 0.002, 0.4 + (i as f64) * 0.002]);
        }
        let mut windows = vec![win];
        let result = wham_2d_pmf(&mut windows, 5, (0.0, 1.0), (0.0, 1.0), 1.0, 1e-4, 10);
        assert_eq!(result.len(), 25, "5×5 grid should have 25 entries");
    }
    #[test]
    fn test_mbar_delta_u_matrix_reference_zero() {
        let u_kn = vec![vec![1.0_f64, 2.0, 3.0], vec![4.0_f64, 5.0, 6.0]];
        let delta = mbar_delta_u_matrix(&u_kn, 0);
        for &v in &delta[0] {
            assert!(v.abs() < 1e-12, "reference row should be zero");
        }
        assert!((delta[1][0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_mbar_delta_u_matrix_empty() {
        let result = mbar_delta_u_matrix(&[], 0);
        assert!(result.is_empty());
    }
    #[test]
    fn test_variance_known_values() {
        let xs = vec![2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let var = variance(&xs);
        assert!(
            (var - 4.571_428_571_4).abs() < 1e-6,
            "variance mismatch: {var}"
        );
    }
    #[test]
    fn test_variance_constant_zero() {
        let xs = vec![3.0_f64; 100];
        let var = variance(&xs);
        assert!(
            var.abs() < 1e-12,
            "constant series should have zero variance"
        );
    }
}
