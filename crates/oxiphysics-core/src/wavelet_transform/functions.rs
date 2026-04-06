//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::{
    CwtResult, DwtDecomposition, DwtLevel, ModwtDecomposition, MotherWavelet,
    MultiresolutionAnalysis, Scalogram, SwtDecomposition, ThresholdMode, WaveletFamily,
    WaveletPacketNode, WaveletPacketTree,
};

/// Haar low-pass filter coefficients (db1).
pub(super) const HAAR_LO: [f64; 2] = [
    std::f64::consts::FRAC_1_SQRT_2,
    std::f64::consts::FRAC_1_SQRT_2,
];
/// Haar high-pass filter coefficients (db1).
pub(super) const HAAR_HI: [f64; 2] = [
    std::f64::consts::FRAC_1_SQRT_2,
    -std::f64::consts::FRAC_1_SQRT_2,
];
/// Daubechies-2 low-pass filter coefficients.
pub(super) const DB2_LO: [f64; 4] = [
    0.4829629131445341,
    0.8365163037378079,
    0.2241438680420134,
    -0.1294095225512604,
];
/// Daubechies-2 high-pass filter coefficients.
pub(super) const DB2_HI: [f64; 4] = [
    -0.1294095225512604,
    -0.2241438680420134,
    0.8365163037378079,
    -0.4829629131445341,
];
/// Daubechies-3 low-pass filter coefficients.
pub(super) const DB3_LO: [f64; 6] = [
    0.3326705529500826,
    0.8068915093110925,
    0.4598775021184915,
    -0.1350110200102546,
    -0.0854412738820267,
    0.0352262918857095,
];
/// Daubechies-3 high-pass filter coefficients.
pub(super) const DB3_HI: [f64; 6] = [
    0.0352262918857095,
    0.0854412738820267,
    -0.1350110200102546,
    -0.4598775021184915,
    0.8068915093110925,
    -0.3326705529500826,
];
/// Daubechies-4 low-pass filter coefficients.
pub(super) const DB4_LO: [f64; 8] = [
    0.2303778133088965,
    0.7148465705529156,
    0.6308807679298589,
    -0.0279837694169839,
    -0.1870348117190931,
    0.0308413818355607,
    0.0328830116668852,
    -0.0105974017850690,
];
/// Daubechies-4 high-pass filter coefficients.
pub(super) const DB4_HI: [f64; 8] = [
    -0.0105974017850690,
    -0.0328830116668852,
    0.0308413818355607,
    0.1870348117190931,
    -0.0279837694169839,
    -0.6308807679298589,
    0.7148465705529156,
    -0.2303778133088965,
];
/// Daubechies-5 low-pass filter coefficients.
pub(super) const DB5_LO: [f64; 10] = [
    0.160_102_397_974_193,
    0.6038292697971898,
    0.7243085284377729,
    0.1384281459013204,
    -0.2422948870663824,
    -0.0322448695846381,
    0.0775714938400459,
    -0.0062414902127983,
    -0.0125807519990820,
    0.0033357252854738,
];
/// Daubechies-5 high-pass filter coefficients.
pub(super) const DB5_HI: [f64; 10] = [
    0.0033357252854738,
    0.0125807519990820,
    -0.0062414902127983,
    -0.0775714938400459,
    -0.0322448695846381,
    0.2422948870663824,
    0.1384281459013204,
    -0.7243085284377729,
    0.6038292697971898,
    -0.160_102_397_974_193,
];
/// Daubechies-6 low-pass filter coefficients.
pub(super) const DB6_LO: [f64; 12] = [
    0.1115407433501095,
    0.4946238903984533,
    0.7511339080210959,
    0.3152503517091982,
    -0.226_264_693_965_44,
    -0.1297668675672625,
    0.0975016055873225,
    0.0275228655303053,
    -0.0315820393174862,
    0.0005538422011614,
    0.0047772575109455,
    -0.0010773010853085,
];
/// Daubechies-6 high-pass filter coefficients.
pub(super) const DB6_HI: [f64; 12] = [
    -0.0010773010853085,
    -0.0047772575109455,
    0.0005538422011614,
    0.0315820393174862,
    0.0275228655303053,
    -0.0975016055873225,
    -0.1297668675672625,
    0.226_264_693_965_44,
    0.3152503517091982,
    -0.7511339080210959,
    0.4946238903984533,
    -0.1115407433501095,
];
/// Convolve `signal` with `filter` and downsample by 2 (periodic extension).
pub fn convolve_downsample(signal: &[f64], filter: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let flen = filter.len();
    let out_len = (n + flen - 1) / 2;
    let mut result = Vec::with_capacity(out_len);
    for i in 0..out_len {
        let mut sum = 0.0;
        for (j, &h) in filter.iter().enumerate() {
            let idx = (2 * i + flen - 1).wrapping_sub(j) % n;
            sum += h * signal[idx];
        }
        result.push(sum);
    }
    result
}
/// Upsample by 2 (zero-insertion) and convolve with `filter` (periodic extension).
pub fn upsample_convolve(coeffs: &[f64], filter: &[f64], target_len: usize) -> Vec<f64> {
    let _flen = filter.len();
    let mut result = vec![0.0; target_len];
    for (k, &c) in coeffs.iter().enumerate() {
        for (j, &h) in filter.iter().enumerate() {
            let idx = (2 * k + j) % target_len;
            result[idx] += c * h;
        }
    }
    result
}
/// Perform a single-level forward DWT.
pub fn dwt_single(signal: &[f64], wavelet: WaveletFamily) -> DwtLevel {
    let lo = wavelet.lo_dec();
    let hi = wavelet.hi_dec();
    DwtLevel {
        approx: convolve_downsample(signal, lo),
        detail: convolve_downsample(signal, hi),
    }
}
/// Perform a single-level inverse DWT.
pub fn idwt_single(level: &DwtLevel, wavelet: WaveletFamily, target_len: usize) -> Vec<f64> {
    let lo_r = wavelet.lo_rec();
    let hi_r = wavelet.hi_rec();
    let a = upsample_convolve(&level.approx, &lo_r, target_len);
    let d = upsample_convolve(&level.detail, &hi_r, target_len);
    a.iter().zip(d.iter()).map(|(&x, &y)| x + y).collect()
}
/// Perform multi-level forward DWT decomposition.
///
/// `levels` is the number of decomposition levels. If the signal is too
/// short for the requested number of levels, fewer levels are computed.
pub fn dwt(signal: &[f64], wavelet: WaveletFamily, levels: usize) -> DwtDecomposition {
    let mut current = signal.to_vec();
    let mut details = Vec::with_capacity(levels);
    let mut lengths = Vec::with_capacity(levels);
    let filter_len = wavelet.filter_length();
    for _level in 0..levels {
        if current.len() < filter_len {
            break;
        }
        lengths.push(current.len());
        let decomp = dwt_single(&current, wavelet);
        details.push(decomp.detail);
        current = decomp.approx;
    }
    DwtDecomposition {
        details,
        approx: current,
        wavelet,
        lengths,
    }
}
/// Perform multi-level inverse DWT reconstruction.
pub fn idwt(decomp: &DwtDecomposition) -> Vec<f64> {
    let mut current = decomp.approx.clone();
    let n_levels = decomp.details.len();
    for i in (0..n_levels).rev() {
        let target_len = decomp.lengths[i];
        let level = DwtLevel {
            approx: current,
            detail: decomp.details[i].clone(),
        };
        current = idwt_single(&level, decomp.wavelet, target_len);
    }
    current
}
/// Perform multiresolution analysis using the given wavelet.
///
/// Returns approximations and detail contributions at each level.
pub fn multiresolution_analysis(
    signal: &[f64],
    wavelet: WaveletFamily,
    levels: usize,
) -> MultiresolutionAnalysis {
    let decomp = dwt(signal, wavelet, levels);
    let n_levels = decomp.details.len();
    let mut approximations = Vec::with_capacity(n_levels + 1);
    let mut detail_contributions = Vec::with_capacity(n_levels);
    let mut current_approx = decomp.approx.clone();
    approximations.push(current_approx.clone());
    for i in (0..n_levels).rev() {
        let target_len = decomp.lengths[i];
        let zero_approx = vec![0.0; current_approx.len()];
        let detail_level = DwtLevel {
            approx: zero_approx,
            detail: decomp.details[i].clone(),
        };
        let detail_contrib = idwt_single(&detail_level, wavelet, target_len);
        detail_contributions.push(detail_contrib);
        let full_level = DwtLevel {
            approx: current_approx,
            detail: decomp.details[i].clone(),
        };
        current_approx = idwt_single(&full_level, wavelet, target_len);
        approximations.push(current_approx.clone());
    }
    detail_contributions.reverse();
    MultiresolutionAnalysis {
        approximations,
        detail_contributions,
    }
}
/// Compute Shannon entropy of a coefficient vector (normalized).
pub fn shannon_entropy(coeffs: &[f64]) -> f64 {
    let total_energy: f64 = coeffs.iter().map(|&c| c * c).sum();
    if total_energy < 1e-30 {
        return 0.0;
    }
    let mut entropy = 0.0;
    for &c in coeffs {
        let p = (c * c) / total_energy;
        if p > 1e-30 {
            entropy -= p * p.ln();
        }
    }
    entropy
}
/// Perform wavelet packet decomposition.
///
/// Builds a full binary tree of wavelet packet nodes up to `max_level` levels.
pub fn wavelet_packet_decompose(
    signal: &[f64],
    wavelet: WaveletFamily,
    max_level: usize,
) -> WaveletPacketTree {
    let root = WaveletPacketNode {
        coefficients: signal.to_vec(),
        level: 0,
        position: 0,
        entropy: shannon_entropy(signal),
    };
    let mut nodes: Vec<Vec<WaveletPacketNode>> = vec![vec![root]];
    let filter_len = wavelet.filter_length();
    for level in 0..max_level {
        let mut next_level = Vec::new();
        for node in &nodes[level] {
            if node.coefficients.len() < filter_len {
                continue;
            }
            let lo = wavelet.lo_dec();
            let hi = wavelet.hi_dec();
            let approx = convolve_downsample(&node.coefficients, lo);
            let detail = convolve_downsample(&node.coefficients, hi);
            next_level.push(WaveletPacketNode {
                entropy: shannon_entropy(&approx),
                coefficients: approx,
                level: level + 1,
                position: 2 * node.position,
            });
            next_level.push(WaveletPacketNode {
                entropy: shannon_entropy(&detail),
                coefficients: detail,
                level: level + 1,
                position: 2 * node.position + 1,
            });
        }
        if next_level.is_empty() {
            break;
        }
        nodes.push(next_level);
    }
    WaveletPacketTree {
        nodes,
        wavelet,
        max_level,
    }
}
/// Select best basis from wavelet packet tree using minimum entropy criterion.
///
/// Returns indices (level, position) of the selected nodes.
pub fn best_basis_selection(tree: &WaveletPacketTree) -> Vec<(usize, usize)> {
    let max_lvl = tree.nodes.len() - 1;
    if max_lvl == 0 {
        return vec![(0, 0)];
    }
    let mut selected: Vec<Vec<bool>> = tree
        .nodes
        .iter()
        .map(|level| vec![true; level.len()])
        .collect();
    for level in (0..max_lvl).rev() {
        for (i, node) in tree.nodes[level].iter().enumerate() {
            let child_base = 2 * i;
            if level + 1 < tree.nodes.len() && child_base + 1 < tree.nodes[level + 1].len() {
                let child_entropy = tree.nodes[level + 1][child_base].entropy
                    + tree.nodes[level + 1][child_base + 1].entropy;
                if node.entropy <= child_entropy {
                    selected[level][i] = true;
                    deselect_subtree(&mut selected, level + 1, child_base);
                    deselect_subtree(&mut selected, level + 1, child_base + 1);
                } else {
                    selected[level][i] = false;
                }
            }
        }
    }
    let mut result = Vec::new();
    for (level, level_selected) in selected.iter().enumerate() {
        for (pos, &sel) in level_selected.iter().enumerate() {
            if sel {
                result.push((level, pos));
            }
        }
    }
    result
}
/// Recursively deselect a subtree rooted at (level, position).
pub fn deselect_subtree(selected: &mut [Vec<bool>], level: usize, pos: usize) {
    if level >= selected.len() || pos >= selected[level].len() {
        return;
    }
    selected[level][pos] = false;
    if level + 1 < selected.len() {
        deselect_subtree(selected, level + 1, 2 * pos);
        deselect_subtree(selected, level + 1, 2 * pos + 1);
    }
}
/// Compute the continuous wavelet transform.
///
/// # Arguments
/// * `signal` - Input signal.
/// * `scales` - Array of scales at which to compute the CWT.
/// * `wavelet` - Mother wavelet to use.
/// * `dt` - Sampling period.
pub fn cwt(signal: &[f64], scales: &[f64], wavelet: MotherWavelet, dt: f64) -> CwtResult {
    let n = signal.len();
    let mut coefficients = Vec::with_capacity(scales.len());
    for &scale in scales {
        let mut row = Vec::with_capacity(n);
        let norm = (dt / scale).sqrt();
        for b in 0..n {
            let mut sum = 0.0;
            for (k, &s) in signal.iter().enumerate() {
                let t = ((k as f64) - (b as f64)) * dt / scale;
                sum += s * wavelet.evaluate(t);
            }
            row.push(sum * norm * dt);
        }
        coefficients.push(row);
    }
    CwtResult {
        coefficients,
        scales: scales.to_vec(),
        wavelet,
    }
}
/// Compute the inverse CWT using the Morlet reconstruction formula.
///
/// This is an approximate reconstruction using the admissibility constant.
pub fn icwt(cwt_result: &CwtResult, dt: f64) -> Vec<f64> {
    let n = if cwt_result.coefficients.is_empty() {
        0
    } else {
        cwt_result.coefficients[0].len()
    };
    let mut signal = vec![0.0; n];
    if cwt_result.scales.is_empty() || n == 0 {
        return signal;
    }
    let c_psi = 1.0;
    for (si, &scale) in cwt_result.scales.iter().enumerate() {
        let norm = 1.0 / (scale * scale);
        for b in 0..n {
            signal[b] += cwt_result.coefficients[si][b] * norm * dt;
        }
    }
    let dj = if cwt_result.scales.len() > 1 {
        (cwt_result.scales[1] / cwt_result.scales[0]).ln()
    } else {
        1.0
    };
    for val in &mut signal {
        *val *= dj / c_psi;
    }
    signal
}
/// Compute the scalogram from a CWT result.
pub fn scalogram(cwt_result: &CwtResult) -> Scalogram {
    let mut energy = Vec::with_capacity(cwt_result.coefficients.len());
    let mut scale_energy = Vec::with_capacity(cwt_result.scales.len());
    for row in &cwt_result.coefficients {
        let e_row: Vec<f64> = row.iter().map(|&c| c * c).collect();
        let total: f64 = e_row.iter().sum();
        scale_energy.push(total);
        energy.push(e_row);
    }
    Scalogram {
        energy,
        scales: cwt_result.scales.clone(),
        scale_energy,
    }
}
/// Compute the global wavelet spectrum (time-averaged energy at each scale).
pub fn global_wavelet_spectrum(scalo: &Scalogram) -> Vec<f64> {
    scalo
        .energy
        .iter()
        .map(|row| {
            if row.is_empty() {
                0.0
            } else {
                row.iter().sum::<f64>() / row.len() as f64
            }
        })
        .collect()
}
/// Apply thresholding to a coefficient vector.
pub fn apply_threshold(coeffs: &[f64], threshold: f64, mode: ThresholdMode) -> Vec<f64> {
    match mode {
        ThresholdMode::Hard => coeffs
            .iter()
            .map(|&c| if c.abs() < threshold { 0.0 } else { c })
            .collect(),
        ThresholdMode::Soft => coeffs
            .iter()
            .map(|&c| {
                if c.abs() < threshold {
                    0.0
                } else {
                    c.signum() * (c.abs() - threshold)
                }
            })
            .collect(),
    }
}
/// Compute the universal (VisuShrink) threshold.
///
/// `sigma` is the noise standard deviation.
/// `n` is the signal length.
pub fn universal_threshold(sigma: f64, n: usize) -> f64 {
    sigma * (2.0 * (n as f64).ln()).sqrt()
}
/// Estimate noise standard deviation from the finest detail coefficients.
///
/// Uses the median absolute deviation (MAD) estimator.
pub fn estimate_noise_sigma(detail_coeffs: &[f64]) -> f64 {
    if detail_coeffs.is_empty() {
        return 0.0;
    }
    let mut abs_coeffs: Vec<f64> = detail_coeffs.iter().map(|c| c.abs()).collect();
    abs_coeffs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if abs_coeffs.len().is_multiple_of(2) {
        (abs_coeffs[abs_coeffs.len() / 2 - 1] + abs_coeffs[abs_coeffs.len() / 2]) / 2.0
    } else {
        abs_coeffs[abs_coeffs.len() / 2]
    };
    median / 0.6745
}
/// Denoise a signal using wavelet thresholding.
///
/// # Arguments
/// * `signal` - Noisy input signal.
/// * `wavelet` - Wavelet family to use.
/// * `levels` - Number of decomposition levels.
/// * `mode` - Thresholding mode (hard or soft).
/// * `threshold` - Optional threshold value. If `None`, uses universal threshold.
#[allow(clippy::too_many_arguments)]
pub fn wavelet_denoise(
    signal: &[f64],
    wavelet: WaveletFamily,
    levels: usize,
    mode: ThresholdMode,
    threshold: Option<f64>,
) -> Vec<f64> {
    let mut decomp = dwt(signal, wavelet, levels);
    let thresh = threshold.unwrap_or_else(|| {
        if decomp.details.is_empty() {
            0.0
        } else {
            let sigma = estimate_noise_sigma(&decomp.details[0]);
            universal_threshold(sigma, signal.len())
        }
    });
    for detail in &mut decomp.details {
        *detail = apply_threshold(detail, thresh, mode);
    }
    idwt(&decomp)
}
/// Energy of wavelet coefficients at a given level.
pub fn level_energy(coeffs: &[f64]) -> f64 {
    coeffs.iter().map(|&c| c * c).sum()
}
/// Compute energy distribution across all DWT levels.
///
/// Returns a vector of (level_index, energy) pairs, plus the approximation energy.
pub fn energy_distribution(decomp: &DwtDecomposition) -> (Vec<f64>, f64) {
    let detail_energies: Vec<f64> = decomp.details.iter().map(|d| level_energy(d)).collect();
    let approx_energy = level_energy(&decomp.approx);
    (detail_energies, approx_energy)
}
/// Compute the relative energy (percentage) at each level.
pub fn relative_energy(decomp: &DwtDecomposition) -> Vec<f64> {
    let (detail_energies, approx_energy) = energy_distribution(decomp);
    let total: f64 = detail_energies.iter().sum::<f64>() + approx_energy;
    if total < 1e-30 {
        return vec![0.0; detail_energies.len() + 1];
    }
    let mut result: Vec<f64> = detail_energies.iter().map(|&e| e / total).collect();
    result.push(approx_energy / total);
    result
}
/// Compute the wavelet entropy (measure of signal complexity).
pub fn wavelet_entropy(decomp: &DwtDecomposition) -> f64 {
    let (detail_energies, approx_energy) = energy_distribution(decomp);
    let total: f64 = detail_energies.iter().sum::<f64>() + approx_energy;
    if total < 1e-30 {
        return 0.0;
    }
    let mut entropy = 0.0;
    for &e in detail_energies
        .iter()
        .chain(std::iter::once(&approx_energy))
    {
        let p = e / total;
        if p > 1e-30 {
            entropy -= p * p.ln();
        }
    }
    entropy
}
/// Perform a stationary wavelet transform (SWT).
///
/// Also known as the "algorithme a trous" (with holes).
pub fn swt(signal: &[f64], wavelet: WaveletFamily, levels: usize) -> SwtDecomposition {
    let n = signal.len();
    let lo = wavelet.lo_dec();
    let hi = wavelet.hi_dec();
    let flen = lo.len();
    let mut current = signal.to_vec();
    let mut details = Vec::with_capacity(levels);
    for level in 0..levels {
        let step = 1_usize << level;
        let mut approx = vec![0.0; n];
        let mut detail = vec![0.0; n];
        for i in 0..n {
            let mut sum_lo = 0.0;
            let mut sum_hi = 0.0;
            for j in 0..flen {
                let idx = (i + j * step) % n;
                sum_lo += lo[j] * current[idx];
                sum_hi += hi[j] * current[idx];
            }
            approx[i] = sum_lo;
            detail[i] = sum_hi;
        }
        details.push(detail);
        current = approx;
    }
    SwtDecomposition {
        details,
        approx: current,
        wavelet,
    }
}
/// Compute wavelet cross-spectrum between two CWT results.
///
/// Returns the cross-spectrum matrix (scale x time).
pub fn wavelet_cross_spectrum(cwt_x: &CwtResult, cwt_y: &CwtResult) -> Vec<Vec<f64>> {
    let n_scales = cwt_x.scales.len().min(cwt_y.scales.len());
    let mut cross = Vec::with_capacity(n_scales);
    for s in 0..n_scales {
        let n_time = cwt_x.coefficients[s].len().min(cwt_y.coefficients[s].len());
        let row: Vec<f64> = (0..n_time)
            .map(|t| cwt_x.coefficients[s][t] * cwt_y.coefficients[s][t])
            .collect();
        cross.push(row);
    }
    cross
}
/// Compute wavelet coherence between two signals.
///
/// Returns values between 0 and 1 at each (scale, time) point.
pub fn wavelet_coherence(cwt_x: &CwtResult, cwt_y: &CwtResult, smoothing: usize) -> Vec<Vec<f64>> {
    let cross = wavelet_cross_spectrum(cwt_x, cwt_y);
    let scalo_x = scalogram(cwt_x);
    let scalo_y = scalogram(cwt_y);
    let n_scales = cross.len();
    let mut coherence = Vec::with_capacity(n_scales);
    for s in 0..n_scales {
        let n_time = cross[s].len();
        let mut coh_row = Vec::with_capacity(n_time);
        for t in 0..n_time {
            let lo = t.saturating_sub(smoothing);
            let hi_t = (t + smoothing + 1).min(n_time);
            let mut sum_cross = 0.0;
            let mut sum_xx = 0.0;
            let mut sum_yy = 0.0;
            for k in lo..hi_t {
                sum_cross += cross[s][k];
                sum_xx += scalo_x.energy[s][k];
                sum_yy += scalo_y.energy[s][k];
            }
            let denom = (sum_xx * sum_yy).sqrt();
            let c = if denom > 1e-30 {
                (sum_cross.abs() / denom).min(1.0)
            } else {
                0.0
            };
            coh_row.push(c);
        }
        coherence.push(coh_row);
    }
    coherence
}
/// Compute wavelet variance (energy per scale, normalized by signal length).
pub fn wavelet_variance(cwt_result: &CwtResult) -> Vec<f64> {
    cwt_result
        .coefficients
        .iter()
        .map(|row| {
            let n = row.len().max(1) as f64;
            row.iter().map(|&c| c * c).sum::<f64>() / n
        })
        .collect()
}
/// Compute instantaneous wavelet power at each (scale, time) point.
pub fn wavelet_power(cwt_result: &CwtResult) -> Vec<Vec<f64>> {
    cwt_result
        .coefficients
        .iter()
        .map(|row| row.iter().map(|&c| c * c).collect())
        .collect()
}
/// Generate logarithmically spaced scales for CWT.
///
/// # Arguments
/// * `s0` - Smallest scale.
/// * `num_scales` - Number of scales.
/// * `dj` - Scale spacing in octaves (e.g. 0.25 for 4 scales per octave).
pub fn log_scales(s0: f64, num_scales: usize, dj: f64) -> Vec<f64> {
    (0..num_scales)
        .map(|j| s0 * (2.0_f64).powf(j as f64 * dj))
        .collect()
}
/// Generate linearly spaced scales for CWT.
pub fn linear_scales(s_min: f64, s_max: f64, num_scales: usize) -> Vec<f64> {
    if num_scales <= 1 {
        return vec![s_min];
    }
    let step = (s_max - s_min) / (num_scales - 1) as f64;
    (0..num_scales).map(|i| s_min + i as f64 * step).collect()
}
/// Detect ridges in a CWT scalogram (local maxima along scale axis).
///
/// Returns positions of ridge points as (scale_index, time_index).
pub fn detect_ridges(scalo: &Scalogram) -> Vec<(usize, usize)> {
    let n_scales = scalo.energy.len();
    if n_scales < 3 {
        return Vec::new();
    }
    let mut ridges = Vec::new();
    let n_time = if scalo.energy.is_empty() {
        0
    } else {
        scalo.energy[0].len()
    };
    for t in 0..n_time {
        for s in 1..n_scales - 1 {
            if s < scalo.energy.len()
                && t < scalo.energy[s].len()
                && t < scalo.energy[s - 1].len()
                && t < scalo.energy[s + 1].len()
            {
                let e = scalo.energy[s][t];
                if e > scalo.energy[s - 1][t] && e > scalo.energy[s + 1][t] && e > 1e-30 {
                    ridges.push((s, t));
                }
            }
        }
    }
    ridges
}
/// Compress a signal by keeping only the top `fraction` of wavelet coefficients.
///
/// Returns the denoised/compressed signal.
pub fn wavelet_compress(
    signal: &[f64],
    wavelet: WaveletFamily,
    levels: usize,
    fraction: f64,
) -> Vec<f64> {
    let mut decomp = dwt(signal, wavelet, levels);
    let mut all_mags: Vec<f64> = decomp
        .details
        .iter()
        .flat_map(|d| d.iter().map(|&c| c.abs()))
        .chain(decomp.approx.iter().map(|&c| c.abs()))
        .collect();
    all_mags.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let keep_count = ((all_mags.len() as f64 * fraction).ceil() as usize).min(all_mags.len());
    let threshold = if keep_count > 0 && keep_count <= all_mags.len() {
        all_mags[keep_count.saturating_sub(1)]
    } else {
        0.0
    };
    for detail in &mut decomp.details {
        for c in detail.iter_mut() {
            if c.abs() < threshold {
                *c = 0.0;
            }
        }
    }
    for c in decomp.approx.iter_mut() {
        if c.abs() < threshold {
            *c = 0.0;
        }
    }
    idwt(&decomp)
}
/// Compute the compression ratio (fraction of non-zero coefficients).
pub fn compression_ratio(decomp: &DwtDecomposition) -> f64 {
    let total: usize = decomp.details.iter().map(|d| d.len()).sum::<usize>() + decomp.approx.len();
    let nonzero: usize = decomp
        .details
        .iter()
        .flat_map(|d| d.iter())
        .chain(decomp.approx.iter())
        .filter(|&&c| c.abs() > 1e-30)
        .count();
    if total == 0 {
        0.0
    } else {
        nonzero as f64 / total as f64
    }
}
/// Compute the MODWT (non-decimated DWT with rescaled filters).
pub fn modwt(signal: &[f64], wavelet: WaveletFamily, levels: usize) -> ModwtDecomposition {
    let n = signal.len();
    let lo_orig = wavelet.lo_dec();
    let hi_orig = wavelet.hi_dec();
    let flen = lo_orig.len();
    let scale = std::f64::consts::FRAC_1_SQRT_2;
    let lo: Vec<f64> = lo_orig.iter().map(|&h| h * scale).collect();
    let hi: Vec<f64> = hi_orig.iter().map(|&h| h * scale).collect();
    let mut current = signal.to_vec();
    let mut details = Vec::with_capacity(levels);
    for level in 0..levels {
        let step = 1_usize << level;
        let mut approx = vec![0.0; n];
        let mut detail = vec![0.0; n];
        for i in 0..n {
            let mut s_lo = 0.0;
            let mut s_hi = 0.0;
            for j in 0..flen {
                let idx = (i + n - j * step) % n;
                s_lo += lo[j] * current[idx];
                s_hi += hi[j] * current[idx];
            }
            approx[i] = s_lo;
            detail[i] = s_hi;
        }
        details.push(detail);
        current = approx;
    }
    ModwtDecomposition {
        details,
        approx: current,
        wavelet,
    }
}
/// Compute the cone of influence for CWT.
///
/// The cone of influence defines the region where edge effects become important.
/// Returns the e-folding time for each scale.
pub fn cone_of_influence(scales: &[f64], wavelet: MotherWavelet) -> Vec<f64> {
    scales
        .iter()
        .map(|&s| match wavelet {
            MotherWavelet::Morlet { .. } => s * 2.0_f64.sqrt(),
            MotherWavelet::MexicanHat => s * (2.0_f64.sqrt()),
        })
        .collect()
}
/// Convert CWT scale to approximate pseudo-frequency.
///
/// f = center_frequency / (scale * dt)
pub fn scale_to_frequency(scale: f64, dt: f64, wavelet: MotherWavelet) -> f64 {
    let center_freq = match wavelet {
        MotherWavelet::Morlet { omega0 } => omega0 / (2.0 * PI),
        MotherWavelet::MexicanHat => 2.0 / (PI * (2.0_f64 / 3.0).sqrt()),
    };
    center_freq / (scale * dt)
}
/// Convert frequency to CWT scale.
pub fn frequency_to_scale(freq: f64, dt: f64, wavelet: MotherWavelet) -> f64 {
    let center_freq = match wavelet {
        MotherWavelet::Morlet { omega0 } => omega0 / (2.0 * PI),
        MotherWavelet::MexicanHat => 2.0 / (PI * (2.0_f64 / 3.0).sqrt()),
    };
    center_freq / (freq * dt)
}
/// Compute wavelet-based signal features for classification/analysis.
///
/// Returns: mean energy, std energy, max coefficient, min coefficient per level.
pub fn wavelet_features(decomp: &DwtDecomposition) -> Vec<[f64; 4]> {
    let mut features = Vec::new();
    for detail in &decomp.details {
        if detail.is_empty() {
            features.push([0.0; 4]);
            continue;
        }
        let n = detail.len() as f64;
        let energy: f64 = detail.iter().map(|&c| c * c).sum();
        let mean_energy = energy / n;
        let mean = detail.iter().sum::<f64>() / n;
        let var = detail.iter().map(|&c| (c - mean).powi(2)).sum::<f64>() / n;
        let std_energy = var.sqrt();
        let max_c = detail.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_c = detail.iter().cloned().fold(f64::INFINITY, f64::min);
        features.push([mean_energy, std_energy, max_c, min_c]);
    }
    if !decomp.approx.is_empty() {
        let n = decomp.approx.len() as f64;
        let energy: f64 = decomp.approx.iter().map(|&c| c * c).sum();
        let mean = decomp.approx.iter().sum::<f64>() / n;
        let var = decomp
            .approx
            .iter()
            .map(|&c| (c - mean).powi(2))
            .sum::<f64>()
            / n;
        let max_c = decomp
            .approx
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_c = decomp.approx.iter().cloned().fold(f64::INFINITY, f64::min);
        features.push([energy / n, var.sqrt(), max_c, min_c]);
    }
    features
}
/// Compute the reconstruction error (L2 norm) between original and reconstructed signals.
pub fn reconstruction_error(original: &[f64], reconstructed: &[f64]) -> f64 {
    let n = original.len().min(reconstructed.len());
    let err: f64 = (0..n)
        .map(|i| (original[i] - reconstructed[i]).powi(2))
        .sum();
    (err / n as f64).sqrt()
}
/// Compute the signal-to-noise ratio (SNR) of the reconstruction.
pub fn reconstruction_snr(original: &[f64], reconstructed: &[f64]) -> f64 {
    let n = original.len().min(reconstructed.len());
    let signal_power: f64 = (0..n).map(|i| original[i].powi(2)).sum();
    let noise_power: f64 = (0..n)
        .map(|i| (original[i] - reconstructed[i]).powi(2))
        .sum();
    if noise_power < 1e-30 {
        return f64::INFINITY;
    }
    10.0 * (signal_power / noise_power).log10()
}
/// Compute the BayesShrink adaptive threshold for each level.
///
/// The BayesShrink threshold is: sigma_noise^2 / sigma_signal
/// where sigma_signal = max(0, sigma_x^2 - sigma_noise^2).
pub fn bayes_shrink_threshold(detail_coeffs: &[f64], noise_sigma: f64) -> f64 {
    let n = detail_coeffs.len() as f64;
    if n < 1.0 {
        return 0.0;
    }
    let sigma_x_sq = detail_coeffs.iter().map(|&c| c * c).sum::<f64>() / n;
    let sigma_n_sq = noise_sigma * noise_sigma;
    let sigma_s_sq = (sigma_x_sq - sigma_n_sq).max(0.0);
    if sigma_s_sq < 1e-30 {
        detail_coeffs
            .iter()
            .map(|c| c.abs())
            .fold(0.0_f64, f64::max)
    } else {
        sigma_n_sq / sigma_s_sq.sqrt()
    }
}
/// Denoise using BayesShrink (adaptive, level-dependent threshold).
pub fn bayes_shrink_denoise(
    signal: &[f64],
    wavelet: WaveletFamily,
    levels: usize,
    mode: ThresholdMode,
) -> Vec<f64> {
    let mut decomp = dwt(signal, wavelet, levels);
    let noise_sigma = if decomp.details.is_empty() {
        0.0
    } else {
        estimate_noise_sigma(&decomp.details[0])
    };
    for detail in &mut decomp.details {
        let thresh = bayes_shrink_threshold(detail, noise_sigma);
        *detail = apply_threshold(detail, thresh, mode);
    }
    idwt(&decomp)
}
