//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BasquinCurve, DamageToleranceResult, ParisLaw, RainflowCycle};

/// Rainflow cycle counting using the ASTM E1049 simplified 3-point algorithm.
///
/// Returns a list of `(range, mean)` pairs for each counted cycle.
/// Half-cycles (residual cycles) are combined as full cycles where possible.
///
/// # Arguments
/// * `signal` - Time-domain load or strain signal (peaks/valleys, or raw signal)
///
/// # Returns
/// A `Vec<(f64, f64)>` of `(range, mean)` tuples for each full cycle.
pub fn rain_flow_count(signal: &[f64]) -> Vec<(f64, f64)> {
    if signal.len() < 3 {
        return Vec::new();
    }
    let tp = extract_turning_points(signal);
    let mut cycles: Vec<(f64, f64)> = Vec::new();
    let mut stack: Vec<f64> = Vec::new();
    let mut i = 0;
    while i < tp.len() {
        stack.push(tp[i]);
        i += 1;
        loop {
            let n = stack.len();
            if n < 3 {
                break;
            }
            let x0 = stack[n - 3];
            let x1 = stack[n - 2];
            let x2 = stack[n - 1];
            let range_01 = (x1 - x0).abs();
            let range_12 = (x2 - x1).abs();
            if range_12 >= range_01 {
                let range = range_01;
                let mean = (x0 + x1) / 2.0;
                cycles.push((range, mean));
                stack.remove(n - 3);
                stack.remove(n - 3);
            } else {
                break;
            }
        }
    }
    let residuals = stack;
    let mut j = 0;
    while j + 1 < residuals.len() {
        let x0 = residuals[j];
        let x1 = residuals[j + 1];
        let range = (x1 - x0).abs();
        let mean = (x0 + x1) / 2.0;
        cycles.push((range, mean));
        j += 2;
    }
    cycles
}
/// Rainflow counting that also returns half-cycle information.
///
/// Returns `Vec<(range, mean, count)>` where `count` is 1.0 for full cycles
/// and 0.5 for half cycles.
pub fn rain_flow_count_with_half_cycles(signal: &[f64]) -> Vec<(f64, f64, f64)> {
    if signal.len() < 3 {
        return Vec::new();
    }
    let tp = extract_turning_points(signal);
    let mut cycles: Vec<(f64, f64, f64)> = Vec::new();
    let mut stack: Vec<f64> = Vec::new();
    let mut i = 0;
    while i < tp.len() {
        stack.push(tp[i]);
        i += 1;
        loop {
            let n = stack.len();
            if n < 3 {
                break;
            }
            let x0 = stack[n - 3];
            let x1 = stack[n - 2];
            let x2 = stack[n - 1];
            let range_01 = (x1 - x0).abs();
            let range_12 = (x2 - x1).abs();
            if range_12 >= range_01 {
                let range = range_01;
                let mean = (x0 + x1) / 2.0;
                cycles.push((range, mean, 1.0));
                stack.remove(n - 3);
                stack.remove(n - 3);
            } else {
                break;
            }
        }
    }
    let residuals = stack;
    let mut j = 0;
    while j + 1 < residuals.len() {
        let x0 = residuals[j];
        let x1 = residuals[j + 1];
        let range = (x1 - x0).abs();
        let mean = (x0 + x1) / 2.0;
        cycles.push((range, mean, 0.5));
        j += 2;
    }
    cycles
}
/// Compute Palmgren-Miner damage from rainflow-counted cycles using an S-N curve.
///
/// For each cycle (range, mean), the stress amplitude = range/2 is used to
/// compute cycles to failure N_f from the S-N curve, and damage D += 1/N_f.
pub fn miner_damage_from_rainflow(cycles: &[(f64, f64)], sn: &BasquinCurve) -> f64 {
    let mut damage = 0.0;
    for &(range, _mean) in cycles {
        let stress_amp = range / 2.0;
        let nf = sn.cycles_to_failure(stress_amp);
        if nf.is_finite() && nf > 0.0 {
            damage += 1.0 / nf;
        }
    }
    damage
}
/// Extract turning points (peaks and valleys) from a signal.
///
/// The first and last points are always included.
pub(super) fn extract_turning_points(signal: &[f64]) -> Vec<f64> {
    if signal.len() <= 2 {
        return signal.to_vec();
    }
    let mut tp = vec![signal[0]];
    for k in 1..signal.len() - 1 {
        let prev = signal[k - 1];
        let curr = signal[k];
        let next = signal[k + 1];
        if ((curr >= prev && curr >= next) || (curr <= prev && curr <= next))
            && *tp.last().expect("collection should not be empty") != curr
        {
            tp.push(curr);
        }
    }
    let last = *signal.last().expect("collection should not be empty");
    if *tp.last().expect("collection should not be empty") != last {
        tp.push(last);
    }
    tp
}
/// Compute stress concentration factor for a notched specimen.
///
/// Kt is the theoretical stress concentration factor.
/// Kf is the fatigue notch factor: Kf = 1 + q*(Kt - 1)
/// where q is the notch sensitivity (0 ≤ q ≤ 1).
pub fn fatigue_notch_factor(kt: f64, notch_sensitivity: f64) -> f64 {
    1.0 + notch_sensitivity * (kt - 1.0)
}
/// Peterson's notch sensitivity estimate.
///
/// q = 1 / (1 + a/r)
///
/// where `a` is the material constant (mm) and `r` is the notch radius (mm).
/// For steel, a ≈ 0.0254 * (2070 / σ_u)^1.8 \[MPa, mm\].
pub fn peterson_notch_sensitivity(notch_radius_mm: f64, material_constant_mm: f64) -> f64 {
    1.0 / (1.0 + material_constant_mm / notch_radius_mm)
}
/// Estimate Peterson's material constant `a` for steel.
///
/// a (mm) ≈ 0.0254 * (2070 / σ_u_mpa)^1.8
pub fn peterson_material_constant_steel(ultimate_strength_mpa: f64) -> f64 {
    0.0254 * (2070.0 / ultimate_strength_mpa).powf(1.8)
}
/// Fit a Basquin S-N curve (log-log linear) to a set of (N, sigma) data points
/// using ordinary least squares in log-space.
///
/// The model is: log(σ) = log(A) + B * log(N)
/// i.e., σ = A * N^B
///
/// Returns (A, B) on success, None if fewer than 2 data points.
pub fn fit_sn_curve_least_squares(data: &[(f64, f64)]) -> Option<(f64, f64)> {
    let n = data.len();
    if n < 2 {
        return None;
    }
    let log_data: Vec<(f64, f64)> = data
        .iter()
        .filter(|&&(cycles, sigma)| cycles > 0.0 && sigma > 0.0)
        .map(|&(cycles, sigma)| (cycles.ln(), sigma.ln()))
        .collect();
    let m = log_data.len();
    if m < 2 {
        return None;
    }
    let mean_x: f64 = log_data.iter().map(|(x, _)| x).sum::<f64>() / m as f64;
    let mean_y: f64 = log_data.iter().map(|(_, y)| y).sum::<f64>() / m as f64;
    let ssxy: f64 = log_data
        .iter()
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum();
    let ssxx: f64 = log_data.iter().map(|(x, _)| (x - mean_x).powi(2)).sum();
    if ssxx.abs() < 1e-30 {
        return None;
    }
    let b_exp = ssxy / ssxx;
    let log_a = mean_y - b_exp * mean_x;
    let a = log_a.exp();
    Some((a, b_exp))
}
/// Fit S-N curve and return a `BasquinCurve` (with endurance limit estimated at N=1e7).
pub fn fit_basquin_from_data(data: &[(f64, f64)], endurance_fraction: f64) -> Option<BasquinCurve> {
    let (a, b_exp) = fit_sn_curve_least_squares(data)?;
    let endurance = a * (1e7_f64).powf(b_exp) * endurance_fraction;
    Some(BasquinCurve::new(a, b_exp, endurance.max(0.0)))
}
/// R² coefficient of determination for a fitted S-N curve.
pub fn sn_r_squared(data: &[(f64, f64)], a: f64, b_exp: f64) -> f64 {
    let valid: Vec<(f64, f64)> = data
        .iter()
        .filter(|&&(n, s)| n > 0.0 && s > 0.0)
        .copied()
        .collect();
    if valid.len() < 2 {
        return 0.0;
    }
    let mean_s: f64 = valid.iter().map(|(_, s)| s.ln()).sum::<f64>() / valid.len() as f64;
    let ss_tot: f64 = valid.iter().map(|(_, s)| (s.ln() - mean_s).powi(2)).sum();
    let ss_res: f64 = valid
        .iter()
        .map(|(n, s)| (s.ln() - (a * n.powf(b_exp)).ln()).powi(2))
        .sum();
    if ss_tot < 1e-30 {
        return 1.0;
    }
    1.0 - ss_res / ss_tot
}
/// Gerber mean stress correction.
///
/// σ_ar = σ_a / (1 - (σ_m / σ_u)²)
///
/// Gerber is less conservative than Goodman for compressive mean stress.
pub fn gerber_equivalent_amplitude(
    stress_amplitude: f64,
    mean_stress: f64,
    ultimate_strength: f64,
) -> f64 {
    let ratio = mean_stress / ultimate_strength;
    let denom = 1.0 - ratio * ratio;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    stress_amplitude / denom
}
/// Soderberg mean stress correction.
///
/// σ_ar = σ_a / (1 - σ_m / σ_ys)
pub fn soderberg_equivalent_amplitude(
    stress_amplitude: f64,
    mean_stress: f64,
    yield_strength: f64,
) -> f64 {
    let denom = 1.0 - mean_stress / yield_strength;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    stress_amplitude / denom
}
/// Goodman mean stress correction (see also `GoodmanDiagram::equivalent_amplitude`).
///
/// σ_ar = σ_a / (1 - σ_m / σ_u)
pub fn goodman_equivalent_amplitude(
    stress_amplitude: f64,
    mean_stress: f64,
    ultimate_strength: f64,
) -> f64 {
    let denom = 1.0 - mean_stress / ultimate_strength;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    stress_amplitude / denom
}
/// Walker mean stress correction.
///
/// σ_ar = σ_max^(1-γ) * σ_a^γ
/// where γ is the Walker exponent (typically 0.3..0.7 for metals).
pub fn walker_equivalent_amplitude(
    stress_max: f64,
    stress_amplitude: f64,
    walker_exponent: f64,
) -> f64 {
    if stress_max <= 0.0 {
        return stress_amplitude;
    }
    stress_max.powf(1.0 - walker_exponent) * stress_amplitude.powf(walker_exponent)
}
/// Compute the von Mises equivalent stress amplitude for a multiaxial
/// stress state.
///
/// The stress amplitude components are given in Voigt notation:
/// \[σ_x_amp, σ_y_amp, σ_z_amp, τ_xy_amp, τ_yz_amp, τ_xz_amp\]
///
/// σ_vm_amp = sqrt(((σx-σy)² + (σy-σz)² + (σz-σx)²)/2 + 3*(τxy²+τyz²+τxz²))
pub fn von_mises_amplitude(sigma_amp: [f64; 6]) -> f64 {
    let sx = sigma_amp[0];
    let sy = sigma_amp[1];
    let sz = sigma_amp[2];
    let txy = sigma_amp[3];
    let tyz = sigma_amp[4];
    let txz = sigma_amp[5];
    let term1 = (sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2);
    let term2 = 6.0 * (txy.powi(2) + tyz.powi(2) + txz.powi(2));
    ((term1 + term2) / 2.0).sqrt()
}
/// Compute the multiaxial Sines fatigue criterion value.
///
/// Sines criterion: f = J2_a + α * I1_mean ≤ β
/// where J2_a = sqrt(J2_amplitude) and I1_mean = hydrostatic mean stress.
///
/// Returns the Sines criterion value (failure if >= beta).
pub fn sines_criterion(stress_amp: [f64; 6], mean_stress: [f64; 6], alpha: f64) -> f64 {
    let vm_amp = von_mises_amplitude(stress_amp) / 3.0_f64.sqrt();
    let i1_mean = mean_stress[0] + mean_stress[1] + mean_stress[2];
    vm_amp + alpha * i1_mean
}
/// Critical plane method: find the plane with maximum damage.
///
/// Scans through `n_planes` candidate planes and returns the maximum
/// shear amplitude on any plane.
pub fn critical_plane_max_shear(sigma_amp: [f64; 6], n_planes: usize) -> f64 {
    let mut max_shear = 0.0_f64;
    let n = n_planes.max(1);
    for k in 0..n {
        let theta = std::f64::consts::PI * k as f64 / n as f64;
        let nx = theta.cos();
        let ny = theta.sin();
        let tx = sigma_amp[0] * nx + sigma_amp[3] * ny;
        let ty = sigma_amp[3] * nx + sigma_amp[1] * ny;
        let tn = tx * nx + ty * ny;
        let shear = ((tx - tn * nx).powi(2) + (ty - tn * ny).powi(2)).sqrt();
        max_shear = max_shear.max(shear);
    }
    max_shear
}
/// Perform a damage tolerance analysis using Paris law.
pub fn damage_tolerance_analysis(
    paris: &ParisLaw,
    a_initial: f64,
    a_critical: f64,
    delta_k_fn: &impl Fn(f64) -> f64,
    da_step: f64,
    inspection_interval: f64,
    safety_factor: f64,
) -> DamageToleranceResult {
    let cycles_to_critical = paris.cycles_to_grow(a_initial, a_critical, delta_k_fn, da_step);
    DamageToleranceResult {
        a_initial,
        a_critical,
        inspection_interval,
        cycles_to_critical,
        safety_factor,
    }
}
/// Non-linear damage accumulation (Marco-Starkey modification).
///
/// D_total = Σ (n_i/N_i)^α
/// where α is typically > 1 for load-order effects.
pub fn marco_starkey_damage(blocks: &[(f64, f64)], exponent: f64) -> f64 {
    blocks.iter().map(|&(n, nf)| (n / nf).powf(exponent)).sum()
}
/// Corten-Dolan cumulative damage model.
///
/// D = Σ n_i * σ_i^d / (N_1 * σ_1^d)
/// where σ_1 is the reference stress, N_1 is cycles at reference, d is the exponent.
pub fn corten_dolan_damage(
    blocks: &[(f64, f64, f64)],
    sigma_ref: f64,
    n_ref: f64,
    d_exp: f64,
) -> f64 {
    if n_ref <= 0.0 || sigma_ref <= 0.0 {
        return 0.0;
    }
    let numerator: f64 = blocks
        .iter()
        .map(|&(n, sigma, _)| n * sigma.powf(d_exp))
        .sum();
    numerator / (n_ref * sigma_ref.powf(d_exp))
}
/// Theoretical stress concentration factor Kt for an edge notch.
///
/// Uses the Inglis approximation: Kt ≈ 1 + 2 * sqrt(half_width / notch_radius).
///
/// # Arguments
/// * `notch_radius` - Notch tip radius ρ (m).
/// * `half_width`   - Half the notch depth (or semi-axis) a (m).
pub fn stress_concentration_kt(notch_radius: f64, half_width: f64) -> f64 {
    if notch_radius <= 0.0 {
        return f64::INFINITY;
    }
    1.0 + 2.0 * (half_width / notch_radius).sqrt()
}
/// Effective stress range accounting for the R-ratio.
///
/// Δσ_eff = σ_max * (1 - R) where R = σ_min / σ_max.
/// Alternatively computed directly from σ_max and σ_min.
///
/// # Arguments
/// * `sigma_max` - Maximum stress in cycle (Pa).
/// * `sigma_min` - Minimum stress in cycle (Pa).
/// * `r_ratio`   - Stress ratio R = σ_min / σ_max (used for cross-check only).
pub fn effective_stress_range(sigma_max: f64, sigma_min: f64, _r_ratio: f64) -> f64 {
    (sigma_max - sigma_min).abs()
}
/// Paris law crack growth rate: da/dN = C * ΔK^m.
///
/// # Arguments
/// * `c`        - Paris coefficient C.
/// * `m`        - Paris exponent m.
/// * `delta_k`  - Stress-intensity factor range ΔK (Pa·√m).
pub fn paris_law_crack_growth(c: f64, m: f64, delta_k: f64) -> f64 {
    c * delta_k.powf(m)
}
/// Threshold stress-intensity factor range at a given R-ratio.
///
/// Uses the empirical Walker-type relationship:
/// ΔK_th(R) = ΔK_th0 * (1 - R)^(1 - m_w)   with m_w = 0.5 (default).
///
/// For simplicity a linear reduction is used here:
/// ΔK_th(R) = ΔK_th0 * (1 - R).clamp(0, 1)
///
/// # Arguments
/// * `delta_k_th0` - Threshold at R=0 (Pa·√m).
/// * `r_ratio`     - Stress ratio R.
pub fn threshold_sif_range(delta_k_th0: f64, r_ratio: f64) -> f64 {
    delta_k_th0 * (1.0 - r_ratio).clamp(0.0, 1.0)
}
/// Sequenced Miner's rule: applies damage blocks in order, tracking residual life.
///
/// Returns the damage history at each step.
pub fn sequenced_miner_damage(blocks: &[(f64, f64)]) -> Vec<f64> {
    let mut d = 0.0;
    let mut history = Vec::with_capacity(blocks.len());
    for &(n_applied, n_failure) in blocks {
        if n_failure > 0.0 {
            d += n_applied / n_failure;
        }
        history.push(d);
    }
    history
}
/// Compute the rainflow cycle histogram.
///
/// Groups cycles from rainflow counting into `n_bins` linearly-spaced
/// stress-range bins and returns `(bin_midpoints, cycle_counts)`.
pub fn rainflow_histogram(cycles: &[(f64, f64)], n_bins: usize) -> (Vec<f64>, Vec<f64>) {
    if cycles.is_empty() || n_bins == 0 {
        return (vec![], vec![]);
    }
    let max_range = cycles.iter().map(|&(r, _)| r).fold(0.0_f64, f64::max);
    if max_range <= 0.0 {
        return (vec![0.0; n_bins], vec![0.0; n_bins]);
    }
    let bin_width = max_range / n_bins as f64;
    let mut counts = vec![0.0_f64; n_bins];
    let mut midpoints = vec![0.0_f64; n_bins];
    for (k, mid) in midpoints.iter_mut().enumerate() {
        *mid = (k as f64 + 0.5) * bin_width;
    }
    for &(range, _mean) in cycles {
        let idx = ((range / bin_width) as usize).min(n_bins - 1);
        counts[idx] += 1.0;
    }
    (midpoints, counts)
}
/// Compute damage-equivalent stress range (DESR) for variable amplitude loading.
///
/// DESR = (Σ n_i * Δσ_i^m / N_total)^(1/m)
///
/// where n_i is the number of cycles at range Δσ_i and m is the Paris/S-N
/// exponent. Represents the constant-amplitude range giving the same damage.
pub fn damage_equivalent_stress_range(cycles: &[(f64, f64)], m_exp: f64) -> f64 {
    if cycles.is_empty() {
        return 0.0;
    }
    let n_total = cycles.len() as f64;
    let sum: f64 = cycles.iter().map(|&(range, _)| range.powf(m_exp)).sum();
    (sum / n_total).powf(1.0 / m_exp)
}
/// Compute root-mean-square (RMS) stress range from rainflow cycles.
pub fn rms_stress_range(cycles: &[(f64, f64)]) -> f64 {
    if cycles.is_empty() {
        return 0.0;
    }
    let mean_sq: f64 = cycles.iter().map(|&(r, _)| r * r).sum::<f64>() / cycles.len() as f64;
    mean_sq.sqrt()
}
/// Irregularity factor I = number of upward zero crossings / number of peaks.
///
/// I = 1 for narrow-band loading, 0 < I < 1 for broad-band.
/// Computed from a time-domain signal.
pub fn irregularity_factor(signal: &[f64]) -> f64 {
    if signal.len() < 3 {
        return 0.0;
    }
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    let mut upward_crossings = 0_usize;
    let mut peaks = 0_usize;
    for i in 1..signal.len() - 1 {
        let prev = signal[i - 1] - mean;
        let curr = signal[i] - mean;
        let next = signal[i + 1] - mean;
        if prev < 0.0 && curr >= 0.0 {
            upward_crossings += 1;
        }
        if curr >= prev && curr >= next {
            peaks += 1;
        }
    }
    if peaks == 0 {
        return 0.0;
    }
    upward_crossings as f64 / peaks as f64
}
/// Dirlik's probability density function (closed-form rainflow damage).
///
/// Dirlik's formula gives the cycle range distribution for broad-band random
/// loading as a combination of exponential and Rayleigh distributions.
/// The damage per unit time is estimated without actual rainflow counting.
///
/// Inputs:
/// - m0: variance of the process (0th spectral moment).
/// - m1, m2, m4: 1st, 2nd, 4th spectral moments.
/// - sigma_f, b_exp: Basquin curve parameters (σ_a = σ_f * N^b_exp).
///
/// Returns the expected fatigue damage rate (D per unit time).
pub fn dirlik_damage_rate(m0: f64, m1: f64, m2: f64, m4: f64, c_coeff: f64, m_exp: f64) -> f64 {
    if m0 <= 0.0 || m2 <= 0.0 || m4 <= 0.0 {
        return 0.0;
    }
    let sigma_rms = m0.sqrt();
    let e_p = (m4 / m2).sqrt();
    let _e_0 = (m2 / m0).sqrt();
    let xm = m1 / m0 * (m2 / m4).sqrt();
    let gamma = m2 / (m0 * m4).sqrt();
    let d1 = 2.0 * (xm - gamma * gamma) / (1.0 + gamma * gamma);
    let r = (gamma - xm - d1 * d1) / (1.0 - gamma - d1 + d1 * d1);
    let d2 = (1.0 - gamma - d1 + d1 * d1) / (1.0 - r);
    let d3 = 1.0 - d1 - d2;
    let q = 1.25 * (gamma - d3 - d2 * r) / d1;
    let n_bins = 200_usize;
    let s_max = 6.0 * sigma_rms;
    let ds = s_max / n_bins as f64;
    let mut integral = 0.0_f64;
    for i in 0..n_bins {
        let s = (i as f64 + 0.5) * ds;
        let z = s / sigma_rms;
        let pdf = (d1 / q * (-z / q).exp()
            + d2 * z / (r * r) * (-z * z / (2.0 * r * r)).exp()
            + d3 * z * (-z * z / 2.0).exp())
            / sigma_rms;
        integral += s.powf(m_exp) * pdf * ds;
    }
    c_coeff * e_p * integral
}
/// Standard normal CDF Φ(z) using the rational approximation (Abramowitz & Stegun 26.2.17).
pub(super) fn standard_normal_cdf(z: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.2316419 * z.abs());
    let poly = t
        * (0.319_381_53
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let phi = 1.0 - ((-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()) * poly;
    if z >= 0.0 { phi } else { 1.0 - phi }
}
/// Normal quantile function (inverse CDF) z = Φ⁻¹(p).
///
/// Uses Beasley-Springer-Moro approximation for p ∈ (0, 1).
pub(super) fn normal_quantile(p: f64) -> f64 {
    let p_clamped = p.clamp(1e-12, 1.0 - 1e-12);
    let t = if p_clamped < 0.5 {
        (-2.0 * p_clamped.ln()).sqrt()
    } else {
        (-2.0 * (1.0 - p_clamped).ln()).sqrt()
    };
    let c0 = 2.515_517;
    let c1 = 0.802_853;
    let c2 = 0.010_328;
    let d1 = 1.432_788;
    let d2 = 0.189_269;
    let d3 = 0.001_308;
    let num = c0 + c1 * t + c2 * t * t;
    let den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;
    let z = t - num / den;
    if p_clamped < 0.5 { -z } else { z }
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::fatigue::CoffinManson;

    use crate::fatigue::FatigueLifePredictor;
    use crate::fatigue::GerberDiagram;
    use crate::fatigue::GoodmanDiagram;
    use crate::fatigue::GoodmanDiagramNew;
    use crate::fatigue::MinerRule;
    use crate::fatigue::MinerWithCrit;
    use crate::fatigue::MorrowCorrection;
    use crate::fatigue::NeuberRule;
    use crate::fatigue::PalmgrenMinor;

    use crate::fatigue::RambergOsgood;
    use crate::fatigue::SNScatterBand;
    use crate::fatigue::SNcurve;
    use crate::fatigue::SoderbergDiagram;
    use crate::fatigue::SwtParameter;

    /// Steel-like Coffin-Manson-Basquin parameters (approximate SAE 1015 steel)
    fn steel_sn() -> SNcurve {
        SNcurve::new(1000.0e6, -0.1, 0.2, -0.5, 200.0e9)
    }
    #[test]
    fn test_basquin_elastic_strain_roundtrip() {
        let sn = steel_sn();
        let two_n = 1000.0_f64;
        let eps_e = sn.elastic_strain_amplitude(two_n);
        let eps_p = sn.plastic_strain_amplitude(two_n);
        let eps_total = sn.total_strain_amplitude(two_n);
        let expected_e = (1000.0e6 / 200.0e9) * 1000.0_f64.powf(-0.1);
        assert!(
            (eps_e - expected_e).abs() < 1.0e-15,
            "Elastic strain mismatch: got {eps_e}, expected {expected_e}"
        );
        let expected_p = 0.2 * 1000.0_f64.powf(-0.5);
        assert!(
            (eps_p - expected_p).abs() < 1.0e-15,
            "Plastic strain mismatch: got {eps_p}, expected {expected_p}"
        );
        assert!(
            (eps_total - (expected_e + expected_p)).abs() < 1.0e-15,
            "Total strain mismatch"
        );
    }
    #[test]
    fn test_basquin_stress_from_n() {
        let sn = steel_sn();
        let n_cycles = 1.0e5_f64;
        let sigma_a = sn.stress_amplitude_from_n(n_cycles);
        let expected = 1000.0e6 * (2.0e5_f64).powf(-0.1);
        assert!(
            (sigma_a - expected).abs() / expected < 1.0e-10,
            "Basquin stress mismatch: got {sigma_a}, expected {expected}"
        );
    }
    #[test]
    fn test_cycles_to_failure_roundtrip() {
        let sn = steel_sn();
        let two_n_ref = 2000.0_f64;
        let eps = sn.total_strain_amplitude(two_n_ref);
        let n_back = sn.cycles_to_failure_strain(eps);
        let two_n_back = 2.0 * n_back;
        assert!(
            (two_n_back - two_n_ref).abs() / two_n_ref < 1.0e-6,
            "Roundtrip 2N mismatch: computed {two_n_back}, expected {two_n_ref}"
        );
    }
    #[test]
    fn test_goodman_zero_mean_stress() {
        let gd = GoodmanDiagram::new(600.0e6, 400.0e6);
        let endurance = 250.0e6_f64;
        let fs = 2.0_f64;
        let allowable = gd.allowable_amplitude(0.0, endurance, fs);
        let expected = endurance / fs;
        assert!(
            (allowable - expected).abs() < 1.0e-6,
            "Goodman zero-mean mismatch: got {allowable}, expected {expected}"
        );
    }
    #[test]
    fn test_goodman_is_safe() {
        let gd = GoodmanDiagram::new(600.0e6, 400.0e6);
        let endurance = 250.0e6_f64;
        let sigma_m = 300.0e6_f64;
        let sigma_a = endurance * (1.0 - sigma_m / 600.0e6);
        assert!(
            gd.is_safe(sigma_m, sigma_a, endurance),
            "Point on Goodman line should be safe (<=1)"
        );
        assert!(
            !gd.is_safe(sigma_m, sigma_a * 1.01, endurance),
            "Point above Goodman line should not be safe"
        );
    }
    #[test]
    fn test_soderberg_amplitude() {
        let gd = GoodmanDiagram::new(600.0e6, 400.0e6);
        let endurance = 250.0e6_f64;
        let sigma_m = 200.0e6_f64;
        let amp = gd.soderberg_amplitude(sigma_m, endurance);
        let expected = endurance * (1.0 - sigma_m / 400.0e6);
        assert!(
            (amp - expected).abs() < 1.0e-6,
            "Soderberg amplitude mismatch: got {amp}, expected {expected}"
        );
    }
    #[test]
    fn test_miner_damage_accumulation() {
        let mut miner = PalmgrenMinor::new();
        assert_eq!(miner.total_damage(), 0.0);
        assert!(!miner.is_failed());
        miner.apply_cycle_block(500.0, 1000.0);
        assert!((miner.total_damage() - 0.5).abs() < 1.0e-15);
        assert!(!miner.is_failed());
        miner.apply_cycle_block(500.0, 1000.0);
        assert!((miner.total_damage() - 1.0).abs() < 1.0e-15);
        assert!(miner.is_failed(), "D=1.0 should indicate failure");
    }
    #[test]
    fn test_miner_reset() {
        let mut miner = PalmgrenMinor::new();
        miner.apply_cycle_block(1000.0, 1000.0);
        assert!(miner.is_failed());
        miner.reset();
        assert_eq!(miner.total_damage(), 0.0);
        assert!(!miner.is_failed());
    }
    #[test]
    fn test_rain_flow_count_simple() {
        let signal = [0.0_f64, 1.0, -1.0, 0.0];
        let cycles = rain_flow_count(&signal);
        assert!(
            !cycles.is_empty(),
            "rain_flow_count should return at least one cycle for [0,1,-1,0]"
        );
        assert_eq!(cycles.len(), 2, "Expected 2 cycles, got: {cycles:?}");
        assert!(
            cycles.iter().all(|&(r, _)| (r - 1.0).abs() < 1.0e-10),
            "Expected all cycles with range 1.0, got: {cycles:?}"
        );
    }
    #[test]
    fn test_rain_flow_count_empty() {
        let signal = [0.0_f64, 1.0];
        let cycles = rain_flow_count(&signal);
        assert!(
            cycles.is_empty(),
            "Signal with fewer than 3 points should return no cycles"
        );
    }
    #[test]
    fn test_rain_flow_count_constant() {
        let signal = [1.0_f64; 10];
        let cycles = rain_flow_count(&signal);
        assert!(
            cycles.is_empty(),
            "Constant signal should produce no cycles"
        );
    }
    #[test]
    fn test_basquin_curve_stress_at_n() {
        let bc = BasquinCurve::new(1500.0e6, -0.1, 200.0e6);
        let sigma = bc.stress_at_n(1e6);
        let expected = 1500.0e6 * (1e6_f64).powf(-0.1);
        assert!((sigma - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_basquin_curve_roundtrip() {
        let bc = BasquinCurve::new(1500.0e6, -0.1, 200.0e6);
        let n = 1e5_f64;
        let sigma = bc.stress_at_n(n);
        let n_back = bc.cycles_to_failure(sigma);
        assert!(
            (n_back - n).abs() / n < 1e-6,
            "n_back={n_back}, expected={n}"
        );
    }
    #[test]
    fn test_basquin_below_endurance() {
        let bc = BasquinCurve::new(1500.0e6, -0.1, 200.0e6);
        let n = bc.cycles_to_failure(100.0e6);
        assert!(
            n.is_infinite(),
            "Below endurance limit should give infinite life"
        );
    }
    #[test]
    fn test_basquin_from_two_points() {
        let n1 = 1e3_f64;
        let s1 = 500.0e6;
        let n2 = 1e6;
        let s2 = 300.0e6;
        let bc = BasquinCurve::from_two_points(n1, s1, n2, s2, 100.0e6);
        assert!((bc.stress_at_n(n1) - s1).abs() / s1 < 1e-10);
        assert!((bc.stress_at_n(n2) - s2).abs() / s2 < 1e-10);
    }
    #[test]
    fn test_coffin_manson_strain() {
        let cm = CoffinManson::new(0.2, -0.5);
        let two_n = 1000.0_f64;
        let eps = cm.plastic_strain_amplitude(two_n);
        let expected = 0.2 * 1000.0_f64.powf(-0.5);
        assert!((eps - expected).abs() < 1e-15);
    }
    #[test]
    fn test_coffin_manson_roundtrip() {
        let cm = CoffinManson::new(0.2, -0.5);
        let eps = 0.005;
        let two_n = cm.reversals_to_failure(eps);
        let eps_back = cm.plastic_strain_amplitude(two_n);
        assert!((eps_back - eps).abs() / eps < 1e-10);
    }
    #[test]
    fn test_coffin_manson_cycles() {
        let cm = CoffinManson::new(0.2, -0.5);
        let eps = 0.01;
        let n_cycles = cm.cycles_to_failure(eps);
        let two_n = cm.reversals_to_failure(eps);
        assert!((n_cycles - two_n / 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_transition_life() {
        let cm = CoffinManson::new(0.2, -0.5);
        let two_nt = cm.transition_life(1000.0e6, -0.1, 200.0e9);
        let sn = steel_sn();
        let eps_e = sn.elastic_strain_amplitude(two_nt);
        let eps_p = cm.plastic_strain_amplitude(two_nt);
        assert!(
            (eps_e - eps_p).abs() / eps_e < 0.01,
            "At transition: eps_e={eps_e}, eps_p={eps_p}"
        );
    }
    #[test]
    fn test_morrow_zero_mean() {
        let mc = MorrowCorrection::new(1000.0e6, -0.1, 200.0e9);
        let two_n = 1000.0_f64;
        let eps_no_mean = mc.corrected_strain_amplitude(0.0, two_n);
        let eps_basquin = (1000.0e6 / 200.0e9) * two_n.powf(-0.1);
        assert!((eps_no_mean - eps_basquin).abs() < 1e-15);
    }
    #[test]
    fn test_morrow_tensile_mean_reduces_life() {
        let mc = MorrowCorrection::new(1000.0e6, -0.1, 200.0e9);
        let two_n = 1000.0_f64;
        let eps_zero = mc.corrected_strain_amplitude(0.0, two_n);
        let eps_tensile = mc.corrected_strain_amplitude(200.0e6, two_n);
        assert!(
            eps_tensile < eps_zero,
            "Tensile mean stress should reduce allowable strain"
        );
    }
    #[test]
    fn test_morrow_effective_sigma_f() {
        let mc = MorrowCorrection::new(1000.0e6, -0.1, 200.0e9);
        let eff = mc.effective_sigma_f(300.0e6);
        assert!((eff - 700.0e6).abs() < 1e-6);
    }
    #[test]
    fn test_morrow_cycles_to_failure() {
        let mc = MorrowCorrection::new(1000.0e6, -0.1, 200.0e9);
        let eps = 0.001;
        let n_zero = mc.cycles_to_failure(eps, 0.0);
        let n_tensile = mc.cycles_to_failure(eps, 200.0e6);
        assert!(
            n_tensile < n_zero,
            "Tensile mean stress should reduce fatigue life"
        );
    }
    #[test]
    fn test_swt_compute() {
        let swt = SwtParameter::compute(500.0e6, 0.001);
        assert!((swt - 500.0e3).abs() < 1e-6);
    }
    #[test]
    fn test_swt_compressive_zero() {
        let swt = SwtParameter::compute(-100.0e6, 0.001);
        assert_eq!(swt, 0.0, "Compressive max stress should give zero SWT");
    }
    #[test]
    fn test_swt_material_capacity() {
        let swt = SwtParameter::new(1000.0e6, -0.1, 0.2, -0.5, 200.0e9);
        let cap = swt.material_swt(1000.0);
        assert!(cap > 0.0, "Material SWT capacity should be positive");
    }
    #[test]
    fn test_swt_cycles_to_failure() {
        let swt = SwtParameter::new(1000.0e6, -0.1, 0.2, -0.5, 200.0e9);
        let cap = swt.material_swt(1000.0);
        let n = swt.cycles_to_failure(cap);
        assert!(
            (n - 500.0).abs() / 500.0 < 0.01,
            "SWT roundtrip: expected ~500, got {n}"
        );
    }
    #[test]
    fn test_swt_infinite_life() {
        let swt = SwtParameter::new(1000.0e6, -0.1, 0.2, -0.5, 200.0e9);
        let n = swt.cycles_to_failure(0.0);
        assert!(n.is_infinite());
    }
    #[test]
    fn test_goodman_equivalent_amplitude() {
        let gd = GoodmanDiagram::new(600.0e6, 400.0e6);
        let eq = gd.equivalent_amplitude(0.0, 200.0e6);
        assert!((eq - 200.0e6).abs() < 1e-6);
        let eq2 = gd.equivalent_amplitude(300.0e6, 100.0e6);
        assert!(
            eq2 > 100.0e6,
            "Equivalent amplitude should increase with mean stress"
        );
    }
    #[test]
    fn test_goodman_factor_of_safety() {
        let gd = GoodmanDiagram::new(600.0e6, 400.0e6);
        let endurance = 250.0e6_f64;
        let fs = gd.factor_of_safety(0.0, 125.0e6, endurance);
        assert!((fs - 2.0).abs() < 1e-6, "FS={fs}");
    }
    #[test]
    fn test_miner_apply_spectrum() {
        let mut miner = PalmgrenMinor::new();
        let blocks = vec![(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)];
        miner.apply_spectrum(&blocks);
        assert!((miner.total_damage() - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_miner_remaining_life() {
        let mut miner = PalmgrenMinor::new();
        miner.apply_cycle_block(300.0, 1000.0);
        assert!((miner.remaining_life() - 0.7).abs() < 1e-12);
    }
    #[test]
    fn test_miner_remaining_cycles() {
        let mut miner = PalmgrenMinor::new();
        miner.apply_cycle_block(500.0, 1000.0);
        let remaining = miner.remaining_cycles(2000.0);
        assert!((remaining - 1000.0).abs() < 1e-12);
    }
    #[test]
    fn test_rainflow_with_half_cycles() {
        let signal = [0.0_f64, 2.0, -1.0, 3.0, -2.0, 1.0];
        let cycles = rain_flow_count_with_half_cycles(&signal);
        assert!(!cycles.is_empty());
        for &(_, _, count) in &cycles {
            assert!(count == 0.5 || count == 1.0, "Invalid cycle count: {count}");
        }
    }
    #[test]
    fn test_miner_damage_from_rainflow() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let cycles = vec![(400.0e6, 200.0e6), (200.0e6, 100.0e6)];
        let damage = miner_damage_from_rainflow(&cycles, &sn);
        assert!(damage > 0.0, "Should accumulate some damage");
    }
    #[test]
    fn test_fatigue_notch_factor_no_sensitivity() {
        let kf = fatigue_notch_factor(3.0, 0.0);
        assert!((kf - 1.0).abs() < 1e-12, "q=0 → Kf=1");
    }
    #[test]
    fn test_fatigue_notch_factor_full_sensitivity() {
        let kf = fatigue_notch_factor(3.0, 1.0);
        assert!((kf - 3.0).abs() < 1e-12, "q=1 → Kf=Kt");
    }
    #[test]
    fn test_fatigue_notch_factor_partial() {
        let kf = fatigue_notch_factor(3.0, 0.5);
        assert!((kf - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_peterson_notch_sensitivity() {
        let q = peterson_notch_sensitivity(1.0, 0.1);
        assert!((q - 1.0 / 1.1).abs() < 1e-12);
    }
    #[test]
    fn test_peterson_material_constant() {
        let a = peterson_material_constant_steel(500.0);
        assert!(a > 0.0, "Material constant should be positive");
        let a2 = peterson_material_constant_steel(1000.0);
        assert!(a2 < a, "Higher strength steel should have smaller a");
    }
    #[test]
    fn test_sn_curve_least_squares_two_points() {
        let n1 = 1e3_f64;
        let s1 = 500.0e6;
        let n2 = 1e6_f64;
        let s2 = 300.0e6;
        let (a, b_exp) = fit_sn_curve_least_squares(&[(n1, s1), (n2, s2)]).unwrap();
        assert!(
            (a * n1.powf(b_exp) - s1).abs() / s1 < 1e-6,
            "Point 1 not reproduced"
        );
        assert!(
            (a * n2.powf(b_exp) - s2).abs() / s2 < 1e-6,
            "Point 2 not reproduced"
        );
    }
    #[test]
    fn test_sn_curve_least_squares_insufficient_data() {
        let result = fit_sn_curve_least_squares(&[(1e6, 200.0e6)]);
        assert!(result.is_none(), "Single point should return None");
    }
    #[test]
    fn test_sn_curve_r_squared_perfect_fit() {
        let n1 = 1e3_f64;
        let s1 = 500.0e6;
        let n2 = 1e6_f64;
        let s2 = 300.0e6;
        let (a, b_exp) = fit_sn_curve_least_squares(&[(n1, s1), (n2, s2)]).unwrap();
        let r2 = sn_r_squared(&[(n1, s1), (n2, s2)], a, b_exp);
        assert!(
            (r2 - 1.0).abs() < 1e-8,
            "Perfect 2-point fit should have R²=1: {r2}"
        );
    }
    #[test]
    fn test_fit_basquin_from_data() {
        let data = vec![
            (1e3, 600.0e6),
            (1e4, 450.0e6),
            (1e5, 350.0e6),
            (1e6, 280.0e6),
        ];
        let bc = fit_basquin_from_data(&data, 0.5).unwrap();
        assert!(bc.b_exp < 0.0, "Basquin exponent should be negative");
        assert!(bc.a > 0.0, "Basquin coefficient should be positive");
    }
    #[test]
    fn test_gerber_zero_mean() {
        let sigma_ar = gerber_equivalent_amplitude(200.0e6, 0.0, 600.0e6);
        assert!(
            (sigma_ar - 200.0e6).abs() < 1e-6,
            "Zero mean: Gerber = amplitude"
        );
    }
    #[test]
    fn test_gerber_vs_goodman_conservatism() {
        let sigma_a = 150.0e6;
        let sigma_m = 200.0e6;
        let sigma_u = 600.0e6;
        let goodman = goodman_equivalent_amplitude(sigma_a, sigma_m, sigma_u);
        let gerber = gerber_equivalent_amplitude(sigma_a, sigma_m, sigma_u);
        assert!(
            goodman > gerber,
            "Goodman should be more conservative: G={goodman} < Gerber={gerber}"
        );
    }
    #[test]
    fn test_soderberg_equivalent_amplitude() {
        let sigma_a = 100.0e6;
        let sigma_m = 200.0e6;
        let sigma_ys = 400.0e6;
        let sigma_ar = soderberg_equivalent_amplitude(sigma_a, sigma_m, sigma_ys);
        assert!(
            (sigma_ar - 200.0e6).abs() < 1e-6,
            "Soderberg mismatch: {sigma_ar}"
        );
    }
    #[test]
    fn test_walker_correction() {
        let sigma_max = 400.0e6;
        let sigma_a = 200.0e6;
        let gamma = 0.5;
        let sigma_ar = walker_equivalent_amplitude(sigma_max, sigma_a, gamma);
        assert!(
            sigma_ar > 0.0,
            "Walker equivalent amplitude should be positive"
        );
        let sigma_ar_r_neg1 = walker_equivalent_amplitude(sigma_a, sigma_a, gamma);
        assert!(
            (sigma_ar_r_neg1 - sigma_a).abs() / sigma_a < 1e-6,
            "At R=-1, Walker should equal amplitude"
        );
    }
    #[test]
    fn test_von_mises_amplitude_uniaxial() {
        let sigma_amp = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = von_mises_amplitude(sigma_amp);
        assert!(
            (vm - 300.0e6).abs() < 1e-6,
            "Uniaxial vm amplitude = sigma_x: {vm}"
        );
    }
    #[test]
    fn test_von_mises_amplitude_hydrostatic() {
        let p = 100.0e6;
        let sigma_amp = [p, p, p, 0.0, 0.0, 0.0];
        let vm = von_mises_amplitude(sigma_amp);
        assert!(vm.abs() < 1e-6, "Hydrostatic has zero VM amplitude: {vm}");
    }
    #[test]
    fn test_von_mises_amplitude_pure_shear() {
        let tau = 100.0e6;
        let sigma_amp = [0.0, 0.0, 0.0, tau, 0.0, 0.0];
        let vm = von_mises_amplitude(sigma_amp);
        let expected = 3.0_f64.sqrt() * tau;
        assert!(
            (vm - expected).abs() / expected < 1e-6,
            "Pure shear vm: {vm} vs {expected}"
        );
    }
    #[test]
    fn test_sines_criterion_uniaxial() {
        let sigma_amp = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mean = [0.0; 6];
        let val = sines_criterion(sigma_amp, mean, 0.0);
        assert!(
            val > 0.0,
            "Sines criterion should be positive for non-zero amplitude"
        );
    }
    #[test]
    fn test_critical_plane_max_shear() {
        let sigma_amp = [200.0e6, -100.0e6, 0.0, 50.0e6, 0.0, 0.0];
        let max_shear = critical_plane_max_shear(sigma_amp, 36);
        assert!(max_shear > 0.0, "Max shear should be positive");
    }
    #[test]
    fn test_paris_law_below_threshold() {
        let paris = ParisLaw::new(1e-12, 3.0, 5.0, 100.0);
        let rate = paris.crack_growth_rate(2.0);
        assert_eq!(rate, 0.0, "Below threshold: rate should be zero");
    }
    #[test]
    fn test_paris_law_above_threshold() {
        let paris = ParisLaw::new(1e-12, 3.0, 5.0, 100.0);
        let rate = paris.crack_growth_rate(10.0);
        let expected = 1e-12 * 10.0_f64.powi(3);
        assert!(
            (rate - expected).abs() / expected < 1e-10,
            "Paris rate mismatch: {rate} vs {expected}"
        );
    }
    #[test]
    fn test_paris_law_fast_fracture() {
        let paris = ParisLaw::new(1e-12, 3.0, 5.0, 100.0);
        let rate = paris.crack_growth_rate(100.0);
        assert!(rate.is_infinite(), "At K_c: rate should be infinite");
    }
    #[test]
    fn test_paris_law_cycles_to_grow() {
        let paris = ParisLaw::new(3.8e-11, 3.0, 1.0e6, 50.0e6);
        let delta_sigma = 100.0e6;
        let delta_k_fn = |a: f64| ParisLaw::delta_k_center_crack(delta_sigma, a, 1.0);
        let cycles = paris.cycles_to_grow(0.001, 0.01, &delta_k_fn, 0.0001);
        assert!(cycles.is_some(), "Should compute finite cycles");
        assert!(cycles.unwrap() > 0.0, "Cycles should be positive");
    }
    #[test]
    fn test_paris_law_delta_k_formula() {
        let delta_sigma = 100.0e6;
        let a = 0.01;
        let y = 1.0;
        let dk = ParisLaw::delta_k_center_crack(delta_sigma, a, y);
        let expected = 100.0e6 * (std::f64::consts::PI * 0.01).sqrt();
        assert!(
            (dk - expected).abs() / expected < 1e-10,
            "ΔK mismatch: {dk} vs {expected}"
        );
    }
    #[test]
    fn test_damage_tolerance_analysis_safe() {
        let paris = ParisLaw::new(3.8e-11, 3.0, 1.0, 50.0);
        let delta_sigma = 50.0e6;
        let delta_k_fn = |a: f64| ParisLaw::delta_k_center_crack(delta_sigma, a, 1.0);
        let result =
            damage_tolerance_analysis(&paris, 0.001, 0.01, &delta_k_fn, 0.0001, 1000.0, 2.0);
        assert!(result.a_critical > result.a_initial);
    }
    #[test]
    fn test_damage_tolerance_result_safe_check() {
        let result = DamageToleranceResult {
            a_initial: 0.001,
            a_critical: 0.01,
            inspection_interval: 1000.0,
            cycles_to_critical: Some(5000.0),
            safety_factor: 2.0,
        };
        assert!(result.is_safe(), "5000 cycles > 1000 * 2 = 2000");
    }
    #[test]
    fn test_damage_tolerance_result_unsafe() {
        let result = DamageToleranceResult {
            a_initial: 0.001,
            a_critical: 0.01,
            inspection_interval: 5000.0,
            cycles_to_critical: Some(3000.0),
            safety_factor: 2.0,
        };
        assert!(!result.is_safe(), "3000 < 5000 * 2 = 10000");
    }
    #[test]
    fn test_marco_starkey_damage() {
        let blocks = vec![(100.0, 1000.0), (200.0, 2000.0)];
        let d = marco_starkey_damage(&blocks, 1.0);
        assert!(
            (d - 0.2).abs() < 1e-12,
            "Marco-Starkey at alpha=1 = Miner: {d}"
        );
    }
    #[test]
    fn test_marco_starkey_vs_miner_nonlinear() {
        let blocks = vec![(100.0, 1000.0), (200.0, 2000.0)];
        let d_miner = marco_starkey_damage(&blocks, 1.0);
        let d_nonlinear = marco_starkey_damage(&blocks, 2.0);
        assert!(
            d_nonlinear < d_miner,
            "Nonlinear Marco-Starkey < Miner for alpha>1"
        );
    }
    #[test]
    fn test_sequenced_miner_damage_history() {
        let blocks = vec![(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)];
        let history = sequenced_miner_damage(&blocks);
        assert_eq!(history.len(), 3);
        assert!((history[0] - 0.1).abs() < 1e-12);
        assert!((history[1] - 0.2).abs() < 1e-12);
        assert!((history[2] - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_corten_dolan_damage() {
        let blocks = vec![(100.0, 400.0e6, 1000.0), (200.0, 200.0e6, 5000.0)];
        let d = corten_dolan_damage(&blocks, 400.0e6, 1000.0, 2.0);
        let expected = (100.0 * (400.0e6_f64).powi(2) + 200.0 * (200.0e6_f64).powi(2))
            / (1000.0 * (400.0e6_f64).powi(2));
        assert!(
            (d - expected).abs() < 1e-10,
            "Corten-Dolan: {d} vs {expected}"
        );
    }
    #[test]
    fn test_goodman_diagram_new_amplitude_decreases_with_mean() {
        let gd = GoodmanDiagramNew::new(600.0e6, 300.0e6);
        let a0 = gd.goodman_allowed_amplitude(0.0);
        let a1 = gd.goodman_allowed_amplitude(100.0e6);
        let a2 = gd.goodman_allowed_amplitude(300.0e6);
        assert!(
            (a0 - 300.0e6).abs() < 1.0,
            "Zero mean should give endurance limit"
        );
        assert!(a1 < a0, "Amplitude should decrease with mean stress");
        assert!(a2 < a1, "Amplitude should keep decreasing");
    }
    #[test]
    fn test_gerber_diagram_amplitude_decreases_with_mean() {
        let gd = GerberDiagram::new(600.0e6, 300.0e6);
        let a0 = gd.gerber_allowed_amplitude(0.0);
        let a1 = gd.gerber_allowed_amplitude(300.0e6);
        assert!(
            (a0 - 300.0e6).abs() < 1.0,
            "Zero mean should give endurance limit"
        );
        assert!(a1 < a0, "Gerber amplitude should decrease with mean stress");
    }
    #[test]
    fn test_soderberg_diagram_amplitude_decreases_with_mean() {
        let sd = SoderbergDiagram::new(400.0e6, 200.0e6);
        let a0 = sd.soderberg_allowed_amplitude(0.0);
        let a1 = sd.soderberg_allowed_amplitude(200.0e6);
        assert!(
            (a0 - 200.0e6).abs() < 1.0,
            "Zero mean should give endurance limit"
        );
        assert!(
            a1 < a0,
            "Soderberg amplitude should decrease with mean stress"
        );
        let a_yield = sd.soderberg_allowed_amplitude(400.0e6);
        assert!(
            a_yield.abs() < 1.0,
            "At yield stress, allowable amplitude = 0"
        );
    }
    #[test]
    fn test_miner_rule_failure_at_total_one() {
        let mr = MinerRule::new(vec![0.3, 0.4, 0.3]);
        assert!((mr.total_damage() - 1.0).abs() < 1e-12);
        assert!(mr.is_failed(), "D=1 should indicate failure");
    }
    #[test]
    fn test_miner_rule_no_failure_below_one() {
        let mr = MinerRule::new(vec![0.2, 0.3, 0.4]);
        assert!(!mr.is_failed(), "D=0.9 should not indicate failure");
    }
    #[test]
    fn test_paris_law_crack_growth_positive() {
        let rate = paris_law_crack_growth(1e-12, 3.0, 10.0);
        let expected = 1e-12 * 10.0_f64.powi(3);
        assert!(
            (rate - expected).abs() / expected < 1e-10,
            "Paris rate: {rate}"
        );
        assert!(rate > 0.0);
    }
    #[test]
    fn test_stress_concentration_kt() {
        let kt = stress_concentration_kt(1.0e-3, 4.0e-3);
        let expected = 1.0 + 2.0 * (4.0e-3_f64 / 1.0e-3).sqrt();
        assert!(
            (kt - expected).abs() < 1e-10,
            "Kt mismatch: {kt} vs {expected}"
        );
        assert!(kt > 1.0);
    }
    #[test]
    fn test_threshold_sif_range() {
        let dkth = threshold_sif_range(5.0, 0.0);
        assert!((dkth - 5.0).abs() < 1e-12, "At R=0, threshold = dkth0");
        let dkth_r05 = threshold_sif_range(5.0, 0.5);
        assert!((dkth_r05 - 2.5).abs() < 1e-12, "At R=0.5, threshold = 2.5");
    }
    #[test]
    fn test_fatigue_life_predictor() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let predictor = FatigueLifePredictor::new(sn, 600.0e6);
        let n_zero = predictor.predict_cycles(200.0e6, 0.0);
        let n_mean = predictor.predict_cycles(200.0e6, 200.0e6);
        assert!(n_zero > 0.0, "Life must be positive: {n_zero}");
        assert!(n_mean < n_zero, "Mean stress should reduce fatigue life");
    }
    #[test]
    fn test_effective_stress_range() {
        let range = effective_stress_range(300.0e6, -100.0e6, -1.0 / 3.0);
        assert!(
            (range - 400.0e6).abs() < 1.0,
            "Stress range = sigma_max - sigma_min"
        );
    }
    #[test]
    fn test_neuber_rule_product() {
        let neuber = NeuberRule::new(2.0, 200.0e9, 1200.0e6, 0.15);
        let product = neuber.neuber_product(100.0e6);
        let expected = (2.0 * 100.0e6_f64).powi(2) / 200.0e9;
        assert!(
            (product - expected).abs() / expected < 1e-10,
            "Neuber product mismatch"
        );
    }
    #[test]
    fn test_neuber_ramberg_osgood_strain_positive() {
        let neuber = NeuberRule::new(2.0, 200.0e9, 1200.0e6, 0.15);
        let eps = neuber.ramberg_osgood_strain(200.0e6);
        assert!(eps > 200.0e6 / 200.0e9, "Total strain > elastic component");
    }
    #[test]
    fn test_neuber_local_stress_below_elastic_for_small_nom() {
        let neuber = NeuberRule::new(2.0, 200.0e9, 1200.0e6, 0.15);
        let sigma_nom = 10.0e6;
        let sigma_local = neuber.local_stress(sigma_nom, 100);
        let sigma_elastic = neuber.kt * sigma_nom;
        assert!(
            (sigma_local - sigma_elastic).abs() / sigma_elastic < 0.05,
            "Small nom stress: local ≈ Kt*nom, got {sigma_local} vs {sigma_elastic}"
        );
    }
    #[test]
    fn test_neuber_local_stress_less_than_elastic_at_high_nom() {
        let neuber = NeuberRule::new(3.0, 200.0e9, 800.0e6, 0.15);
        let sigma_nom = 200.0e6;
        let sigma_local = neuber.local_stress(sigma_nom, 100);
        let sigma_elastic = neuber.kt * sigma_nom;
        assert!(
            sigma_local < sigma_elastic,
            "With plasticity: σ_local < Kt*σ_nom ({sigma_local} < {sigma_elastic})"
        );
    }
    #[test]
    fn test_neuber_cyclic_stress_strain() {
        let neuber = NeuberRule::new(2.5, 200.0e9, 1000.0e6, 0.15);
        let (sigma_a, eps_a) = neuber.cyclic_local_stress_strain(80.0e6, 100);
        assert!(sigma_a > 0.0, "Cyclic local stress should be positive");
        assert!(eps_a > 0.0, "Cyclic local strain should be positive");
        let product = sigma_a * eps_a;
        let target = neuber.neuber_product(80.0e6);
        assert!(
            (product - target).abs() / target < 1e-4,
            "Neuber product not satisfied: {product} vs {target}"
        );
    }
    #[test]
    fn test_rainflow_histogram_basic() {
        let cycles = vec![(100.0e6, 0.0), (200.0e6, 10.0e6), (50.0e6, -25.0e6)];
        let (midpoints, counts) = rainflow_histogram(&cycles, 4);
        assert_eq!(midpoints.len(), 4, "Should have 4 bins");
        assert_eq!(counts.len(), 4, "Should have 4 bins");
        let total: f64 = counts.iter().sum();
        assert!((total - 3.0).abs() < 1e-10, "All cycles should be counted");
    }
    #[test]
    fn test_rainflow_histogram_empty() {
        let (mid, cnt) = rainflow_histogram(&[], 5);
        assert!(mid.is_empty() && cnt.is_empty());
    }
    #[test]
    fn test_damage_equivalent_stress_range_uniform() {
        let cycles = vec![(100.0e6, 0.0); 10];
        let desr = damage_equivalent_stress_range(&cycles, 3.0);
        assert!(
            (desr - 100.0e6).abs() / 100.0e6 < 1e-10,
            "Uniform cycles: DESR = range"
        );
    }
    #[test]
    fn test_damage_equivalent_stress_range_positive() {
        let cycles = vec![(50.0e6, 0.0), (100.0e6, 10.0e6), (150.0e6, -50.0e6)];
        let desr = damage_equivalent_stress_range(&cycles, 3.0);
        assert!(desr > 0.0, "DESR should be positive");
        assert!(
            (50.0e6..=150.0e6).contains(&desr),
            "DESR should be between min and max range"
        );
    }
    #[test]
    fn test_rms_stress_range_uniform() {
        let cycles = vec![(100.0e6, 0.0); 5];
        let rms = rms_stress_range(&cycles);
        assert!(
            (rms - 100.0e6).abs() / 100.0e6 < 1e-10,
            "Uniform: RMS = range"
        );
    }
    #[test]
    fn test_rms_stress_range_empty() {
        let rms = rms_stress_range(&[]);
        assert_eq!(rms, 0.0);
    }
    #[test]
    fn test_irregularity_factor_narrow_band() {
        let n = 200_usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 20.0).sin())
            .collect();
        let i_f = irregularity_factor(&signal);
        assert!(i_f > 0.5, "Sinusoid should have I close to 1, got {i_f}");
    }
    #[test]
    fn test_irregularity_factor_short_signal() {
        let result = irregularity_factor(&[1.0, 2.0]);
        assert_eq!(result, 0.0, "Short signal should give 0");
    }
    #[test]
    fn test_sn_scatter_band_median_life() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let band = SNScatterBand::new(sn.clone(), 10.0);
        let sigma = 300.0e6;
        let n_median_direct = sn.cycles_to_failure(sigma);
        let n_at_50 = band.life_at_survival(sigma, 0.5);
        assert!(
            (n_at_50 - n_median_direct).abs() / n_median_direct < 0.01,
            "50% survival should give median life: {n_at_50} vs {n_median_direct}"
        );
    }
    #[test]
    fn test_sn_scatter_band_high_survival_longer_life() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let band = SNScatterBand::new(sn, 10.0);
        let sigma = 300.0e6;
        let n_90 = band.life_at_survival(sigma, 0.9);
        let n_10 = band.life_at_survival(sigma, 0.1);
        assert!(
            n_10 > n_90,
            "10% survival life > 90% survival life (less conservative)"
        );
    }
    #[test]
    fn test_sn_scatter_band_survival_probability() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let band = SNScatterBand::new(sn.clone(), 10.0);
        let sigma = 300.0e6;
        let n_median = sn.cycles_to_failure(sigma);
        let ps = band.survival_probability(sigma, n_median);
        assert!((ps - 0.5).abs() < 0.05, "At median life, P_s ≈ 0.5: {ps}");
    }
    #[test]
    fn test_sn_scatter_band_below_endurance() {
        let sn = BasquinCurve::new(1500.0e6, -0.1, 100.0e6);
        let band = SNScatterBand::new(sn, 10.0);
        let n = band.life_at_survival(50.0e6, 0.5);
        assert!(n.is_infinite(), "Below endurance limit → infinite life");
    }
    #[test]
    fn test_ramberg_osgood_strain_elastic_dominated() {
        let ro = RambergOsgood::new(200.0e9, 1500.0e6, 0.15);
        let sigma = 1.0e6;
        let eps = ro.strain(sigma);
        let eps_e = ro.elastic_strain(sigma);
        assert!(
            (eps - eps_e).abs() / eps_e < 0.01,
            "Very small stress → elastic dominated"
        );
    }
    #[test]
    fn test_ramberg_osgood_plastic_strain_increases_with_stress() {
        let ro = RambergOsgood::new(200.0e9, 800.0e6, 0.15);
        let ep1 = ro.plastic_strain(100.0e6).abs();
        let ep2 = ro.plastic_strain(200.0e6).abs();
        assert!(ep2 > ep1, "Plastic strain should increase with stress");
    }
    #[test]
    fn test_ramberg_osgood_stress_from_strain_roundtrip() {
        let ro = RambergOsgood::new(200.0e9, 1000.0e6, 0.15);
        let sigma_orig = 300.0e6;
        let eps = ro.strain(sigma_orig);
        let sigma_back = ro.stress_from_strain(eps, 100);
        assert!(
            (sigma_back - sigma_orig).abs() / sigma_orig < 1e-5,
            "R-O roundtrip: {sigma_back} vs {sigma_orig}"
        );
    }
    #[test]
    fn test_ramberg_osgood_energy_per_cycle_positive() {
        let ro = RambergOsgood::new(200.0e9, 1000.0e6, 0.15);
        let w = ro.energy_per_cycle(400.0e6);
        assert!(w > 0.0, "Hysteresis energy per cycle should be positive");
    }
    #[test]
    fn test_ramberg_osgood_energy_increases_with_delta_sigma() {
        let ro = RambergOsgood::new(200.0e9, 1000.0e6, 0.15);
        let w1 = ro.energy_per_cycle(200.0e6);
        let w2 = ro.energy_per_cycle(400.0e6);
        assert!(w2 > w1, "More stress range → more hysteresis energy");
    }
    #[test]
    fn test_miner_with_crit_default_no_failure() {
        let mut miner = MinerWithCrit::new(1.0);
        miner.add_block(0.5, 1.0);
        assert!(!miner.is_failed(), "D=0.5 < D_crit=1.0 → no failure");
    }
    #[test]
    fn test_miner_with_crit_failure_at_crit() {
        let mut miner = MinerWithCrit::new(0.7);
        miner.add_block(0.7, 1.0);
        assert!(miner.is_failed(), "D = D_crit → failure");
    }
    #[test]
    fn test_miner_with_crit_remaining_life() {
        let miner = MinerWithCrit {
            damage: 0.3,
            d_crit: 1.0,
        };
        let rem = miner.remaining_life_fraction();
        assert!(
            (rem - 0.7).abs() < 1e-12,
            "Remaining life fraction = 0.7, got {rem}"
        );
    }
    #[test]
    fn test_miner_with_crit_cycles_to_failure() {
        let miner = MinerWithCrit {
            damage: 0.5,
            d_crit: 1.0,
        };
        let n_f = 10_000.0;
        let n_remaining = miner.cycles_to_failure(n_f);
        assert!(
            (n_remaining - 5000.0).abs() < 1e-10,
            "Remaining cycles: {n_remaining}"
        );
    }
    #[test]
    fn test_standard_normal_cdf_at_zero() {
        let phi = standard_normal_cdf(0.0);
        assert!((phi - 0.5).abs() < 1e-4, "Φ(0) ≈ 0.5: {phi}");
    }
    #[test]
    fn test_standard_normal_cdf_symmetry() {
        let phi_pos = standard_normal_cdf(1.0);
        let phi_neg = standard_normal_cdf(-1.0);
        assert!(
            (phi_pos + phi_neg - 1.0).abs() < 1e-5,
            "Symmetry: Φ(z) + Φ(-z) = 1"
        );
    }
    #[test]
    fn test_normal_quantile_at_half() {
        let z = normal_quantile(0.5);
        assert!(z.abs() < 0.01, "Φ⁻¹(0.5) ≈ 0, got {z}");
    }
    #[test]
    fn test_normal_quantile_roundtrip() {
        let p = 0.95;
        let z = normal_quantile(p);
        let p_back = standard_normal_cdf(z);
        assert!(
            (p_back - p).abs() < 0.002,
            "Normal quantile roundtrip: {p_back} vs {p}"
        );
    }
}
/// Perform rainflow cycle counting on a stress/strain time series.
///
/// Implements the ASTM E1049 four-point algorithm.  Returns a list of
/// [`RainflowCycle`] structs with `count = 1.0` for full cycles and
/// `count = 0.5` for residual half-cycles at the end.
pub fn rainflow_count(signal: &[f64]) -> Vec<RainflowCycle> {
    if signal.len() < 2 {
        return Vec::new();
    }
    let mut pts: Vec<f64> = Vec::with_capacity(signal.len());
    pts.push(signal[0]);
    for i in 1..signal.len() - 1 {
        let prev = signal[i - 1];
        let curr = signal[i];
        let next = signal[i + 1];
        if (curr >= prev && curr >= next) || (curr <= prev && curr <= next) {
            pts.push(curr);
        }
    }
    pts.push(*signal.last().expect("collection should not be empty"));
    let mut stack: Vec<f64> = Vec::new();
    let mut cycles: Vec<RainflowCycle> = Vec::new();
    for &pt in &pts {
        stack.push(pt);
        loop {
            let n = stack.len();
            if n < 4 {
                break;
            }
            let x = (stack[n - 2] - stack[n - 3]).abs();
            let y = (stack[n - 1] - stack[n - 2]).abs();
            let z = (stack[n - 3] - stack[n - 4]).abs();
            if y >= x && x <= z {
                let mean = (stack[n - 3] + stack[n - 2]) / 2.0;
                cycles.push(RainflowCycle {
                    range: x,
                    mean,
                    count: 1.0,
                });
                stack.remove(n - 3);
                stack.remove(n - 3);
            } else {
                break;
            }
        }
    }
    let n = stack.len();
    for i in 0..n.saturating_sub(1) {
        let range = (stack[i + 1] - stack[i]).abs();
        let mean = (stack[i + 1] + stack[i]) / 2.0;
        cycles.push(RainflowCycle {
            range,
            mean,
            count: 0.5,
        });
    }
    cycles
}
/// Compute Miner cumulative damage from a slice of (n_applied, N_failure) pairs.
///
/// D = Σ n_i / N_fi
pub fn miners_rule_damage(blocks: &[(f64, f64)]) -> f64 {
    blocks
        .iter()
        .map(|&(n, nf)| if nf > 0.0 { n / nf } else { 0.0 })
        .sum()
}
/// True when Miner damage from given blocks reaches `d_crit`.
pub fn miners_rule_failed(blocks: &[(f64, f64)], d_crit: f64) -> bool {
    miners_rule_damage(blocks) >= d_crit
}
/// Remaining allowable cycles for a constant-amplitude block given existing damage.
///
/// Returns `f64::INFINITY` when no additional damage can accumulate.
pub fn miners_rule_remaining(existing_damage: f64, n_f: f64, d_crit: f64) -> f64 {
    let remaining = (d_crit - existing_damage).max(0.0);
    remaining * n_f
}
#[cfg(test)]
mod tests_basquin_expansion {

    use crate::fatigue::BasquinStressLife;

    use crate::fatigue::CoffinMansonLcf;
    use crate::fatigue::CycleDamageAccumulator;

    use crate::fatigue::RainflowCycle;

    use crate::fatigue::miners_rule_damage;
    use crate::fatigue::miners_rule_failed;
    use crate::fatigue::miners_rule_remaining;
    use crate::fatigue::rainflow_count;
    #[test]
    fn test_basquin_sl_stress_amplitude_at_two_reversals() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 200.0e6);
        let sa = sn.stress_amplitude(2.0);
        let expected = 1000.0e6 * 2.0_f64.powf(-0.1);
        assert!(
            (sa - expected).abs() < 1.0,
            "stress amplitude mismatch: {sa}"
        );
    }
    #[test]
    fn test_basquin_sl_cycles_to_failure_below_endurance_is_inf() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 200.0e6);
        let n_f = sn.cycles_to_failure(100.0e6);
        assert!(n_f.is_infinite(), "below endurance → infinite life");
    }
    #[test]
    fn test_basquin_sl_cycles_to_failure_monotone() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let n1 = sn.cycles_to_failure(400.0e6);
        let n2 = sn.cycles_to_failure(600.0e6);
        assert!(n2 < n1, "higher stress → fewer cycles to failure");
    }
    #[test]
    fn test_basquin_sl_damage_per_cycle_above_endurance() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let d = sn.damage_per_cycle(400.0e6);
        assert!(d > 0.0, "damage > 0 above endurance: {d}");
    }
    #[test]
    fn test_basquin_sl_damage_per_cycle_below_endurance_zero() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 200.0e6);
        let d = sn.damage_per_cycle(100.0e6);
        assert_eq!(d, 0.0, "damage = 0 below endurance");
    }
    #[test]
    fn test_basquin_sl_endurance_ratio() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 500.0e6);
        let ratio = sn.endurance_ratio();
        assert!(
            (ratio - 0.5).abs() < 1e-12,
            "endurance ratio = 0.5: {ratio}"
        );
    }
    #[test]
    fn test_basquin_sl_safety_factor_reduces_life() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let n_no_sf = sn.cycles_to_failure(400.0e6);
        let n_with_sf = sn.cycles_with_safety_factor(400.0e6, 1.5);
        assert!(n_with_sf < n_no_sf, "safety factor should reduce life");
    }
    #[test]
    fn test_basquin_sl_transition_cycles_positive() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 200.0e6);
        let n_t = sn.transition_cycles();
        assert!(n_t > 0.0 && n_t.is_finite(), "transition cycles: {n_t}");
    }
    #[test]
    fn test_coffin_manson_lcf_plastic_strain_amplitude_positive() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let eps = cm.plastic_strain_amplitude(2.0);
        assert!(eps > 0.0, "plastic strain amplitude > 0: {eps}");
    }
    #[test]
    fn test_coffin_manson_lcf_cycles_to_failure_positive() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let n_f = cm.cycles_to_failure(0.01);
        assert!(n_f > 0.0 && n_f.is_finite(), "cycles to failure: {n_f}");
    }
    #[test]
    fn test_coffin_manson_lcf_higher_strain_fewer_cycles() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let n1 = cm.cycles_to_failure(0.005);
        let n2 = cm.cycles_to_failure(0.02);
        assert!(n2 < n1, "higher strain amplitude → fewer cycles");
    }
    #[test]
    fn test_coffin_manson_lcf_transition_reversals_positive() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let n_t = cm.transition_reversals(1000.0e6, 200.0e9, -0.1);
        assert!(n_t > 0.0, "transition reversals > 0: {n_t}");
    }
    #[test]
    fn test_coffin_manson_lcf_damage_per_cycle() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let d = cm.damage_per_cycle(0.01);
        assert!(d > 0.0, "damage per cycle > 0: {d}");
    }
    #[test]
    fn test_coffin_manson_lcf_roundtrip() {
        let cm = CoffinMansonLcf::new(0.5, -0.6);
        let eps = 0.01;
        let n_f = cm.cycles_to_failure(eps);
        let eps_back = cm.plastic_strain_amplitude(2.0 * n_f);
        assert!(
            (eps_back - eps).abs() / eps < 1e-8,
            "roundtrip mismatch: {eps_back} vs {eps}"
        );
    }
    #[test]
    fn test_rainflow_cycle_amplitude() {
        let cyc = RainflowCycle {
            range: 200.0e6,
            mean: 50.0e6,
            count: 1.0,
        };
        assert!(
            (cyc.amplitude() - 100.0e6).abs() < 1.0,
            "amplitude = range/2"
        );
    }
    #[test]
    fn test_rainflow_cycle_stress_ratio_positive_r() {
        let cyc = RainflowCycle {
            range: 200.0e6,
            mean: 50.0e6,
            count: 1.0,
        };
        let r = cyc.stress_ratio();
        assert!((r - (-1.0 / 3.0)).abs() < 1e-10, "stress ratio: {r}");
    }
    #[test]
    fn test_rainflow_count_empty_signal() {
        let cycles = rainflow_count(&[]);
        assert!(cycles.is_empty(), "empty signal → no cycles");
    }
    #[test]
    fn test_rainflow_count_single_point() {
        let cycles = rainflow_count(&[1.0]);
        assert!(cycles.is_empty(), "single point → no cycles");
    }
    #[test]
    fn test_rainflow_count_simple_triangle() {
        let signal = vec![0.0, 10.0, 0.0];
        let cycles = rainflow_count(&signal);
        assert!(!cycles.is_empty(), "should have at least one half-cycle");
    }
    #[test]
    fn test_rainflow_count_positive_ranges() {
        let signal = vec![0.0, 5.0, -5.0, 8.0, -3.0, 6.0, -2.0];
        let cycles = rainflow_count(&signal);
        for c in &cycles {
            assert!(
                c.range >= 0.0,
                "cycle range must be non-negative: {}",
                c.range
            );
        }
    }
    #[test]
    fn test_rainflow_count_total_weighted_cycles() {
        let signal = vec![0.0, 4.0, -4.0, 4.0, -4.0, 0.0];
        let cycles = rainflow_count(&signal);
        let total: f64 = cycles.iter().map(|c| c.count).sum();
        assert!(
            total >= 1.0,
            "should count at least one weighted cycle: {total}"
        );
    }
    #[test]
    fn test_cycle_damage_accumulator_zero_initial() {
        let acc = CycleDamageAccumulator::new(1.0);
        assert_eq!(acc.damage, 0.0);
        assert!(!acc.is_failed());
    }
    #[test]
    fn test_cycle_damage_accumulator_single_block() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        let n_f = sn.cycles_to_failure(400.0e6);
        acc.apply_cycles(n_f / 2.0, 400.0e6, &sn);
        assert!(
            (acc.damage - 0.5).abs() < 1e-10,
            "damage = 0.5: {}",
            acc.damage
        );
    }
    #[test]
    fn test_cycle_damage_accumulator_failure_at_d_crit() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        let n_f = sn.cycles_to_failure(400.0e6);
        acc.apply_cycles(n_f, 400.0e6, &sn);
        assert!(acc.is_failed(), "D = 1 ≥ D_crit → failed");
    }
    #[test]
    fn test_cycle_damage_accumulator_below_endurance_no_damage() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 200.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        acc.apply_cycles(1_000_000.0, 100.0e6, &sn);
        assert_eq!(acc.damage, 0.0, "below endurance → no damage");
    }
    #[test]
    fn test_cycle_damage_accumulator_remaining_life() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        let n_f = sn.cycles_to_failure(400.0e6);
        acc.apply_cycles(n_f * 0.3, 400.0e6, &sn);
        let frac = acc.remaining_life_fraction();
        assert!(
            (frac - 0.7).abs() < 1e-10,
            "remaining fraction = 0.7: {frac}"
        );
    }
    #[test]
    fn test_cycle_damage_accumulator_reset() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        acc.apply_cycles(100.0, 400.0e6, &sn);
        acc.reset();
        assert_eq!(acc.damage, 0.0, "reset should zero damage");
        assert_eq!(acc.block_count(), 0, "reset should clear history");
    }
    #[test]
    fn test_cycle_damage_accumulator_block_count() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        acc.apply_cycles(10.0, 300.0e6, &sn);
        acc.apply_cycles(10.0, 400.0e6, &sn);
        acc.apply_cycles(10.0, 500.0e6, &sn);
        assert_eq!(acc.block_count(), 3, "three blocks applied");
    }
    #[test]
    fn test_cycle_damage_accumulator_remaining_cycles_estimate() {
        let sn = BasquinStressLife::new(1000.0e6, -0.1, 50.0e6);
        let mut acc = CycleDamageAccumulator::new(1.0);
        let n_f = sn.cycles_to_failure(400.0e6);
        acc.apply_cycles(n_f * 0.5, 400.0e6, &sn);
        let rem = acc.remaining_cycles(400.0e6, &sn);
        assert!(
            (rem - n_f * 0.5).abs() / (n_f * 0.5) < 1e-9,
            "remaining cycles ≈ N_f/2: {rem}"
        );
    }
    #[test]
    fn test_miners_rule_damage_zero_for_empty() {
        let d = miners_rule_damage(&[]);
        assert_eq!(d, 0.0, "empty blocks → D = 0");
    }
    #[test]
    fn test_miners_rule_damage_single_block() {
        let blocks = vec![(500.0, 1000.0)];
        let d = miners_rule_damage(&blocks);
        assert!((d - 0.5).abs() < 1e-12, "D = 0.5: {d}");
    }
    #[test]
    fn test_miners_rule_damage_accumulates() {
        let blocks = vec![(200.0, 1000.0), (300.0, 1000.0)];
        let d = miners_rule_damage(&blocks);
        assert!((d - 0.5).abs() < 1e-12, "D = 0.5: {d}");
    }
    #[test]
    fn test_miners_rule_failed_true() {
        let blocks = vec![(1000.0, 1000.0)];
        assert!(miners_rule_failed(&blocks, 1.0), "D ≥ D_crit → failed");
    }
    #[test]
    fn test_miners_rule_failed_false() {
        let blocks = vec![(400.0, 1000.0)];
        assert!(!miners_rule_failed(&blocks, 1.0), "D < D_crit → not failed");
    }
    #[test]
    fn test_miners_rule_remaining() {
        let rem = miners_rule_remaining(0.4, 1000.0, 1.0);
        assert!((rem - 600.0).abs() < 1e-10, "remaining cycles = 600: {rem}");
    }
    #[test]
    fn test_miners_rule_remaining_no_overshoot() {
        let rem = miners_rule_remaining(1.2, 1000.0, 1.0);
        assert_eq!(rem, 0.0, "negative remaining life should clamp to 0");
    }
}
