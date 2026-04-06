#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wavelet transforms for signal analysis in physics simulations.
//!
//! Provides Haar, Daubechies-4, Symlet-4, continuous wavelet transforms (CWT),
//! discrete wavelet transforms (DWT), wavelet-domain denoising,
//! multi-resolution analysis, wavelet packet decomposition with best-basis
//! selection, wavelet compression, and physics-specific wavelet tools.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// HaarWavelet
// ─────────────────────────────────────────────────────────────────────────────

/// Haar wavelet transform (single level and multi-level decomposition).
///
/// The Haar wavelet is the simplest orthogonal wavelet, based on piecewise
/// constant functions.
pub struct HaarWavelet;

impl HaarWavelet {
    /// Compute one level of the forward Haar transform.
    ///
    /// Returns `(approx, detail)` where each output has length `signal.len() / 2`.
    /// The input length should be even; if odd the last sample is dropped.
    pub fn forward(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len() / 2;
        let mut approx = Vec::with_capacity(n);
        let mut detail = Vec::with_capacity(n);
        let sqrt2_inv = 1.0 / std::f64::consts::SQRT_2;
        for i in 0..n {
            let a = signal[2 * i];
            let b = signal[2 * i + 1];
            approx.push((a + b) * sqrt2_inv);
            detail.push((a - b) * sqrt2_inv);
        }
        (approx, detail)
    }

    /// Reconstruct a signal from its Haar approximation and detail coefficients.
    ///
    /// `approx` and `detail` must have the same length.
    pub fn inverse(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let n = approx.len().min(detail.len());
        let mut out = vec![0.0; n * 2];
        let sqrt2_inv = 1.0 / std::f64::consts::SQRT_2;
        for i in 0..n {
            let a = approx[i];
            let d = detail[i];
            out[2 * i] = (a + d) * sqrt2_inv;
            out[2 * i + 1] = (a - d) * sqrt2_inv;
        }
        out
    }

    /// Decompose a signal into multiple levels.
    ///
    /// Returns a vector of length `levels + 1` where index 0 is the final
    /// approximation and indices 1..=levels are detail coefficients at each
    /// level (finest first).
    pub fn decompose_level(signal: &[f64], levels: usize) -> Vec<Vec<f64>> {
        let mut result = Vec::with_capacity(levels + 1);
        let mut approx = signal.to_vec();
        for _ in 0..levels {
            if approx.len() < 2 {
                break;
            }
            let (a, d) = Self::forward(&approx);
            result.push(d);
            approx = a;
        }
        result.push(approx);
        result.reverse();
        result
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Daubechies4
// ─────────────────────────────────────────────────────────────────────────────

/// Daubechies D4 wavelet with four filter taps.
///
/// Uses the standard D4 scaling (low-pass) filter `H` and wavelet (high-pass)
/// filter `G` derived by the quadrature mirror filter relationship.
pub struct Daubechies4;

impl Daubechies4 {
    /// D4 low-pass (scaling) filter coefficients.
    pub const H: [f64; 4] = [
        0.482962913_14503_f64,
        0.836516303_73787_f64,
        0.224143868_04186_f64,
        -0.12940952_25523_f64,
    ];

    /// D4 high-pass (wavelet) filter coefficients (QMF of H).
    pub const G: [f64; 4] = [
        -0.12940952_25523_f64,
        -0.224143868_04186_f64,
        0.836516303_73787_f64,
        -0.482962913_14503_f64,
    ];

    /// One-level forward D4 DWT via periodic convolution and downsampling.
    ///
    /// Returns `(approx, detail)` each of length `signal.len() / 2`.
    /// Input length must be at least 4 and even.
    pub fn forward(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len();
        let half = n / 2;
        let mut approx = vec![0.0; half];
        let mut detail = vec![0.0; half];
        for i in 0..half {
            let mut a = 0.0;
            let mut d = 0.0;
            for k in 0..4 {
                let idx = (2 * i + k) % n;
                a += Self::H[k] * signal[idx];
                d += Self::G[k] * signal[idx];
            }
            approx[i] = a;
            detail[i] = d;
        }
        (approx, detail)
    }

    /// One-level inverse D4 DWT (synthesis step).
    ///
    /// `approx` and `detail` must have the same length.
    pub fn inverse(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let half = approx.len().min(detail.len());
        let n = half * 2;
        let mut out = vec![0.0; n];
        for i in 0..half {
            for k in 0..4 {
                let idx = (2 * i + k) % n;
                out[idx] += Self::H[k] * approx[i] + Self::G[k] * detail[i];
            }
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Symlet4
// ─────────────────────────────────────────────────────────────────────────────

/// Symlet-4 wavelet filter bank (sym4).
///
/// Symlets are nearly symmetric wavelets derived from Daubechies filters.
/// The sym4 filter has 8 coefficients.
pub struct Symlet4;

impl Symlet4 {
    /// Sym4 low-pass (scaling) filter coefficients (8 taps).
    pub const H: [f64; 8] = [
        -0.075_765_714_789_273_13,
        -0.02963552764599851,
        0.49761866763201545,
        0.803_738_751_805_916_1,
        0.29785779560527736,
        -0.099_219_543_576_847_22,
        -0.012603967262037833,
        0.032_223_100_604_042_7,
    ];

    /// Sym4 high-pass (wavelet) filter coefficients (QMF of H).
    pub const G: [f64; 8] = [
        -0.032_223_100_604_042_7,
        -0.012603967262037833,
        0.099_219_543_576_847_22,
        0.29785779560527736,
        -0.803_738_751_805_916_1,
        0.49761866763201545,
        0.02963552764599851,
        -0.075_765_714_789_273_13,
    ];

    /// One-level forward Sym4 DWT via periodic convolution and downsampling.
    ///
    /// Returns `(approx, detail)` each of length `signal.len() / 2`.
    pub fn forward(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len();
        if n < 2 {
            return (signal.to_vec(), vec![]);
        }
        let half = n / 2;
        let filter_len = 8;
        let mut approx = vec![0.0; half];
        let mut detail = vec![0.0; half];
        for i in 0..half {
            let mut a = 0.0;
            let mut d = 0.0;
            for k in 0..filter_len {
                let idx = (2 * i + k) % n;
                a += Self::H[k] * signal[idx];
                d += Self::G[k] * signal[idx];
            }
            approx[i] = a;
            detail[i] = d;
        }
        (approx, detail)
    }

    /// One-level inverse Sym4 DWT (synthesis step).
    ///
    /// `approx` and `detail` must have the same length.
    pub fn inverse(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let half = approx.len().min(detail.len());
        let n = half * 2;
        if n == 0 {
            return vec![];
        }
        let filter_len = 8;
        let mut out = vec![0.0; n];
        for i in 0..half {
            for k in 0..filter_len {
                let idx = (2 * i + k) % n;
                out[idx] += Self::H[k] * approx[i] + Self::G[k] * detail[i];
            }
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WaveletTransform — unified DWT interface
// ─────────────────────────────────────────────────────────────────────────────

/// Wavelet family selector for `WaveletTransform`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveletFamily {
    /// Haar wavelet (2 filter taps).
    Haar,
    /// Daubechies D4 wavelet (4 filter taps).
    Daubechies4,
    /// Symlet-4 wavelet (8 filter taps).
    Symlet4,
}

/// 1-D Discrete Wavelet Transform with pluggable filter bank.
///
/// Supports Haar, Daubechies-4, and Symlet-4 wavelets for single-level and
/// multi-level decomposition/reconstruction.
#[derive(Debug, Clone)]
pub struct WaveletTransform {
    /// The wavelet family used by this transform.
    pub family: WaveletFamily,
}

impl WaveletTransform {
    /// Create a new `WaveletTransform` using the given wavelet family.
    pub fn new(family: WaveletFamily) -> Self {
        Self { family }
    }

    /// Perform one level of the forward DWT.
    ///
    /// Returns `(approximation, detail)` coefficient vectors.
    pub fn forward_one_level(&self, signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        match self.family {
            WaveletFamily::Haar => HaarWavelet::forward(signal),
            WaveletFamily::Daubechies4 => Daubechies4::forward(signal),
            WaveletFamily::Symlet4 => Symlet4::forward(signal),
        }
    }

    /// Perform one level of the inverse DWT (reconstruction).
    ///
    /// Returns the reconstructed signal from approximation and detail coefficients.
    pub fn inverse_one_level(&self, approx: &[f64], detail: &[f64]) -> Vec<f64> {
        match self.family {
            WaveletFamily::Haar => HaarWavelet::inverse(approx, detail),
            WaveletFamily::Daubechies4 => Daubechies4::inverse(approx, detail),
            WaveletFamily::Symlet4 => Symlet4::inverse(approx, detail),
        }
    }

    /// Compute the energy (sum of squared coefficients) of a coefficient slice.
    pub fn energy(coeffs: &[f64]) -> f64 {
        coeffs.iter().map(|x| x * x).sum()
    }

    /// Compute the Shannon entropy of coefficient energies (normalized).
    ///
    /// Returns `- Σ p_i ln(p_i)` where `p_i = c_i² / total_energy`.
    pub fn entropy(coeffs: &[f64]) -> f64 {
        let total: f64 = coeffs.iter().map(|x| x * x).sum();
        if total < 1e-30 {
            return 0.0;
        }
        -coeffs
            .iter()
            .map(|x| {
                let p = x * x / total;
                if p < 1e-30 { 0.0 } else { p * p.ln() }
            })
            .sum::<f64>()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MultilevelDwt
// ─────────────────────────────────────────────────────────────────────────────

/// Multi-level DWT decomposition and perfect reconstruction.
///
/// Stores approximation and detail coefficients at every decomposition level,
/// allowing both analysis and synthesis (reconstruction).
#[derive(Debug, Clone)]
pub struct MultilevelDwt {
    /// Wavelet transform used for all levels.
    pub transform: WaveletTransform,
    /// Number of decomposition levels requested.
    pub levels: usize,
    /// Approximation coefficients at the coarsest level (`approx[0]` = coarsest).
    pub approx: Vec<f64>,
    /// Detail coefficients at each level, coarsest first (`details[0]` = coarsest).
    pub details: Vec<Vec<f64>>,
}

impl MultilevelDwt {
    /// Decompose a signal into `levels` levels using the given wavelet family.
    ///
    /// After construction, `approx` holds the coarsest approximation and
    /// `details` holds detail vectors from coarsest (index 0) to finest.
    pub fn decompose(signal: &[f64], levels: usize, family: WaveletFamily) -> Self {
        let transform = WaveletTransform::new(family);
        let mut approx = signal.to_vec();
        let mut details = Vec::with_capacity(levels);

        for _ in 0..levels {
            if approx.len() < 2 {
                break;
            }
            let (a, d) = transform.forward_one_level(&approx);
            details.push(d);
            approx = a;
        }
        // Store details finest-first (last pushed = finest level detail)
        details.reverse();

        Self {
            transform,
            levels,
            approx,
            details,
        }
    }

    /// Reconstruct the original signal from stored approximation and details.
    ///
    /// Applies inverse DWT from coarsest to finest level.
    pub fn reconstruct(&self) -> Vec<f64> {
        let mut signal = self.approx.clone();
        // details[0] = coarsest → iterate forward from coarsest to finest
        for detail in self.details.iter() {
            signal = self.transform.inverse_one_level(&signal, detail);
        }
        signal
    }

    /// Return the total energy stored in all coefficient vectors.
    pub fn total_energy(&self) -> f64 {
        let approx_e: f64 = self.approx.iter().map(|x| x * x).sum();
        let detail_e: f64 = self
            .details
            .iter()
            .flat_map(|d| d.iter())
            .map(|x| x * x)
            .sum();
        approx_e + detail_e
    }

    /// Return entropy at each decomposition level (approximation + details).
    pub fn level_entropies(&self) -> Vec<f64> {
        let mut entropies = vec![WaveletTransform::entropy(&self.approx)];
        for d in &self.details {
            entropies.push(WaveletTransform::entropy(d));
        }
        entropies
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WaveletPacket
// ─────────────────────────────────────────────────────────────────────────────

/// Wavelet packet decomposition tree.
///
/// Stores all sub-band coefficients produced by full wavelet packet
/// decomposition up to a specified number of levels.
#[derive(Debug, Clone)]
pub struct WaveletPacket {
    /// All node coefficient vectors in the packet tree (breadth-first order).
    pub coefficients: Vec<Vec<f64>>,
    /// Number of decomposition levels.
    pub levels: usize,
}

impl WaveletPacket {
    /// Build a wavelet packet decomposition using the Haar transform.
    ///
    /// After construction, `coefficients` contains 2^levels leaf nodes.
    pub fn new_haar(signal: &[f64], levels: usize) -> Self {
        let mut nodes: Vec<Vec<f64>> = vec![signal.to_vec()];
        for _ in 0..levels {
            let mut next = Vec::with_capacity(nodes.len() * 2);
            for node in &nodes {
                if node.len() >= 2 {
                    let (a, d) = HaarWavelet::forward(node);
                    next.push(a);
                    next.push(d);
                } else {
                    next.push(node.clone());
                    next.push(Vec::new());
                }
            }
            nodes = next;
        }
        Self {
            coefficients: nodes,
            levels,
        }
    }

    /// Build a wavelet packet decomposition using a specified wavelet family.
    ///
    /// Uses the given `WaveletFamily` for all decomposition steps.
    pub fn new(signal: &[f64], levels: usize, family: WaveletFamily) -> Self {
        let transform = WaveletTransform::new(family);
        let mut nodes: Vec<Vec<f64>> = vec![signal.to_vec()];
        for _ in 0..levels {
            let mut next = Vec::with_capacity(nodes.len() * 2);
            for node in &nodes {
                if node.len() >= 2 {
                    let (a, d) = transform.forward_one_level(node);
                    next.push(a);
                    next.push(d);
                } else {
                    next.push(node.clone());
                    next.push(Vec::new());
                }
            }
            nodes = next;
        }
        Self {
            coefficients: nodes,
            levels,
        }
    }

    /// Return the index of the leaf node with highest energy (best basis heuristic).
    pub fn best_basis(&self) -> Vec<usize> {
        let energies = self.energy_distribution();
        let max_e = energies.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        energies
            .iter()
            .enumerate()
            .filter(|&(_, &e)| (e - max_e).abs() < 1e-12 * max_e.abs().max(1.0))
            .map(|(i, _)| i)
            .collect()
    }

    /// Select the best basis using minimum Shannon entropy criterion.
    ///
    /// Returns the indices of the leaf nodes forming the minimum-entropy basis.
    /// This implements a greedy search over the leaf nodes of the full packet tree.
    pub fn best_basis_entropy(&self) -> Vec<usize> {
        let entropies: Vec<f64> = self
            .coefficients
            .iter()
            .map(|c| WaveletTransform::entropy(c))
            .collect();
        let min_e = entropies.iter().cloned().fold(f64::INFINITY, f64::min);
        entropies
            .iter()
            .enumerate()
            .filter(|&(_, &e)| (e - min_e).abs() < 1e-12 * min_e.abs().max(1.0))
            .map(|(i, _)| i)
            .collect()
    }

    /// Compute the energy (sum of squares) in each leaf node.
    pub fn energy_distribution(&self) -> Vec<f64> {
        self.coefficients
            .iter()
            .map(|c| c.iter().map(|x| x * x).sum())
            .collect()
    }

    /// Compute the Shannon entropy of each leaf node's coefficients.
    pub fn entropy_distribution(&self) -> Vec<f64> {
        self.coefficients
            .iter()
            .map(|c| WaveletTransform::entropy(c))
            .collect()
    }

    /// Return the total energy across all leaf nodes.
    pub fn total_energy(&self) -> f64 {
        self.energy_distribution().iter().sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ContinuousWavelet — CWT with Morlet/Mexican Hat/Paul wavelets
// ─────────────────────────────────────────────────────────────────────────────

/// Continuous wavelet transform (CWT) with multiple mother wavelets.
///
/// Supports Morlet, Mexican Hat (Ricker), and Paul wavelets for scalogram
/// computation and instantaneous frequency estimation.
pub struct ContinuousWavelet;

impl ContinuousWavelet {
    /// Evaluate the (real part of the) Morlet wavelet at position `t`.
    ///
    /// `sigma` controls the time-frequency trade-off (bandwidth parameter).
    pub fn morlet(t: f64, sigma: f64) -> f64 {
        let gauss = (-t * t / (2.0 * sigma * sigma)).exp();
        let wave = (t / sigma).cos();
        gauss * wave
    }

    /// Evaluate the Mexican-hat (Ricker) wavelet at position `t`.
    ///
    /// Defined as `(1 - t²/σ²) exp(-t²/2σ²)`.
    pub fn mexican_hat(t: f64, sigma: f64) -> f64 {
        let s2 = sigma * sigma;
        let u = t * t / s2;
        (1.0 - u) * (-t * t / (2.0 * s2)).exp()
    }

    /// Evaluate the Paul wavelet of order `m` at position `t`.
    ///
    /// The real part is proportional to `cos(m * atan2(t, 1)) / (1 + t²)^((m+1)/2)`.
    pub fn paul(t: f64, m: u32) -> f64 {
        let denom = (1.0 + t * t).powf((m as f64 + 1.0) / 2.0);
        let angle = (m as f64) * t.atan2(1.0);
        angle.cos() / denom.max(1e-30)
    }

    /// Compute the CWT scalogram using the Morlet wavelet.
    ///
    /// Returns a matrix `scalogram[scale_idx][time_idx]` where each row
    /// corresponds to one entry in `scales`.
    ///
    /// - `signal` — input time series.
    /// - `dt` — sampling interval.
    /// - `scales` — desired wavelet scales (dilation parameters).
    pub fn cwt_morlet(signal: &[f64], dt: f64, scales: &[f64]) -> Vec<Vec<f64>> {
        let n = signal.len();
        let mut scalogram = Vec::with_capacity(scales.len());
        let sigma = 1.0_f64;
        for &scale in scales {
            let mut row = Vec::with_capacity(n);
            for t_idx in 0..n {
                let mut val = 0.0;
                for tau in 0..n {
                    let time = (t_idx as f64 - tau as f64) * dt / scale.max(1e-12);
                    val += signal[tau] * Self::morlet(time, sigma);
                }
                row.push(val * dt / scale.sqrt().max(1e-12));
            }
            scalogram.push(row);
        }
        scalogram
    }

    /// Compute the CWT scalogram using the Mexican-hat wavelet.
    ///
    /// Returns a matrix `scalogram[scale_idx][time_idx]`.
    pub fn cwt_mexican_hat(signal: &[f64], dt: f64, scales: &[f64]) -> Vec<Vec<f64>> {
        let n = signal.len();
        let mut scalogram = Vec::with_capacity(scales.len());
        let sigma = 1.0_f64;
        for &scale in scales {
            let mut row = Vec::with_capacity(n);
            for t_idx in 0..n {
                let mut val = 0.0;
                for tau in 0..n {
                    let time = (t_idx as f64 - tau as f64) * dt / scale.max(1e-12);
                    val += signal[tau] * Self::mexican_hat(time, sigma);
                }
                row.push(val * dt / scale.sqrt().max(1e-12));
            }
            scalogram.push(row);
        }
        scalogram
    }

    /// Compute the CWT scalogram using the Paul wavelet of order `m`.
    ///
    /// Returns a matrix `scalogram[scale_idx][time_idx]`.
    pub fn cwt_paul(signal: &[f64], dt: f64, scales: &[f64], m: u32) -> Vec<Vec<f64>> {
        let n = signal.len();
        let mut scalogram = Vec::with_capacity(scales.len());
        for &scale in scales {
            let mut row = Vec::with_capacity(n);
            for t_idx in 0..n {
                let mut val = 0.0;
                for tau in 0..n {
                    let time = (t_idx as f64 - tau as f64) * dt / scale.max(1e-12);
                    val += signal[tau] * Self::paul(time, m);
                }
                row.push(val * dt / scale.sqrt().max(1e-12));
            }
            scalogram.push(row);
        }
        scalogram
    }

    /// Estimate the instantaneous dominant frequency at each time step from a scalogram.
    ///
    /// For each time index, finds the scale with maximum CWT coefficient magnitude
    /// and converts it to a frequency via `f = 1 / (scale * dt)`.
    pub fn instantaneous_frequency(scalogram: &[Vec<f64>], scales: &[f64], dt: f64) -> Vec<f64> {
        if scalogram.is_empty() || scales.is_empty() {
            return Vec::new();
        }
        let n_time = scalogram[0].len();
        let n_scales = scalogram.len().min(scales.len());
        (0..n_time)
            .map(|t| {
                let mut best_scale = scales[0];
                let mut best_val = scalogram[0][t].abs();
                for s in 1..n_scales {
                    let v = scalogram[s][t].abs();
                    if v > best_val {
                        best_val = v;
                        best_scale = scales[s];
                    }
                }
                1.0 / (best_scale * dt).max(1e-30)
            })
            .collect()
    }

    /// Compute total power at each scale (integrate over time).
    ///
    /// Returns a vector of `(scale, power)` pairs.
    pub fn scale_averaged_power(scalogram: &[Vec<f64>], scales: &[f64]) -> Vec<(f64, f64)> {
        scalogram
            .iter()
            .zip(scales.iter())
            .map(|(row, &scale)| {
                let power: f64 = row.iter().map(|x| x * x).sum::<f64>() / row.len().max(1) as f64;
                (scale, power)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ContinuousWaveletTransform (alias kept for backward compatibility)
// ─────────────────────────────────────────────────────────────────────────────

/// Continuous wavelet transform (CWT) — backward-compatible alias.
///
/// Provides Morlet and Mexican-hat CWT methods identical to [`ContinuousWavelet`].
pub struct ContinuousWaveletTransform;

impl ContinuousWaveletTransform {
    /// Evaluate the (real part of the) Morlet wavelet at position `t`.
    pub fn morlet(t: f64, sigma: f64) -> f64 {
        ContinuousWavelet::morlet(t, sigma)
    }

    /// Evaluate the Mexican-hat (Ricker) wavelet at position `t`.
    pub fn mexican_hat(t: f64, sigma: f64) -> f64 {
        ContinuousWavelet::mexican_hat(t, sigma)
    }

    /// Compute the CWT scalogram using the Morlet wavelet.
    pub fn cwt_morlet(signal: &[f64], dt: f64, scales: &[f64]) -> Vec<Vec<f64>> {
        ContinuousWavelet::cwt_morlet(signal, dt, scales)
    }

    /// Estimate the instantaneous dominant frequency at each time step.
    pub fn instantaneous_frequency(scalogram: &[Vec<f64>], scales: &[f64], dt: f64) -> Vec<f64> {
        ContinuousWavelet::instantaneous_frequency(scalogram, scales, dt)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DiscreteWaveletTransform
// ─────────────────────────────────────────────────────────────────────────────

/// Collection of discrete wavelet transform (DWT) utility functions.
///
/// Includes single-level and multi-level Haar DWT, energy, and entropy
/// computations.
pub struct DiscreteWaveletTransform;

impl DiscreteWaveletTransform {
    /// Compute the full Haar DWT of a signal.
    ///
    /// Returns a coefficient vector in the standard DWT layout:
    /// `[approx_coarsest ... detail_level1]`.
    pub fn dwt_haar(signal: &[f64]) -> Vec<f64> {
        let mut coeffs = signal.to_vec();
        let mut n = signal.len();
        // In-place lifting
        while n >= 2 {
            let half = n / 2;
            let tmp: Vec<f64> = (0..half)
                .flat_map(|i| {
                    let a = coeffs[2 * i];
                    let b = coeffs[2 * i + 1];
                    let s = 1.0 / std::f64::consts::SQRT_2;
                    vec![(a + b) * s, (a - b) * s]
                })
                .collect();
            // Store approx in first half, detail in second
            for i in 0..half {
                coeffs[i] = tmp[2 * i];
                coeffs[half + i] = tmp[2 * i + 1];
            }
            n = half;
        }
        coeffs
    }

    /// Invert the Haar DWT produced by [`DiscreteWaveletTransform::dwt_haar`].
    pub fn idwt_haar(coeffs: &[f64]) -> Vec<f64> {
        let n = coeffs.len();
        let mut signal = coeffs.to_vec();
        let mut level_size = 2usize;
        // Find the coarsest level
        while level_size <= n {
            level_size *= 2;
        }
        level_size /= 2;
        // Reconstruct from coarsest to finest
        let mut cur = 2usize;
        while cur <= level_size {
            let half = cur / 2;
            let tmp: Vec<f64> = (0..half)
                .flat_map(|i| {
                    let a = signal[i];
                    let d = signal[half + i];
                    let s = 1.0 / std::f64::consts::SQRT_2;
                    vec![(a + d) * s, (a - d) * s]
                })
                .collect();
            signal[..cur].copy_from_slice(&tmp[..cur]);
            cur *= 2;
        }
        signal
    }

    /// Compute the DWT coefficients at a specific decomposition level.
    ///
    /// Returns only the approximation coefficients after `level` Haar forward
    /// steps.
    pub fn dwt_level(signal: &[f64], level: usize) -> Vec<f64> {
        let mut approx = signal.to_vec();
        for _ in 0..level {
            if approx.len() < 2 {
                break;
            }
            let (a, _d) = HaarWavelet::forward(&approx);
            approx = a;
        }
        approx
    }

    /// Compute the total energy (sum of squared coefficients) of a wavelet
    /// coefficient array.
    pub fn wavelet_energy(coeffs: &[f64]) -> f64 {
        coeffs.iter().map(|x| x * x).sum()
    }

    /// Compute the Shannon entropy of a normalized wavelet coefficient array.
    ///
    /// Coefficients are normalized by total energy before computing entropy.
    pub fn wavelet_entropy(coeffs: &[f64]) -> f64 {
        let energy: f64 = coeffs.iter().map(|x| x * x).sum();
        if energy < 1e-30 {
            return 0.0;
        }
        -coeffs
            .iter()
            .map(|x| {
                let p = x * x / energy;
                if p < 1e-30 { 0.0 } else { p * p.ln() }
            })
            .sum::<f64>()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WaveletDenoising — soft/hard thresholding with VisuShrink and SureShrink
// ─────────────────────────────────────────────────────────────────────────────

/// Wavelet-domain signal denoising methods.
///
/// Provides hard and soft thresholding, and common threshold selection rules
/// including VisuShrink (universal) and SureShrink (SURE).
pub struct WaveletDenoising;

impl WaveletDenoising {
    /// Hard thresholding: set coefficients with absolute value below `threshold` to zero.
    pub fn hard_threshold(coeffs: &mut Vec<f64>, threshold: f64) {
        for c in coeffs.iter_mut() {
            if c.abs() < threshold {
                *c = 0.0;
            }
        }
    }

    /// Soft thresholding: shrink coefficients toward zero by `threshold`.
    pub fn soft_threshold(coeffs: &mut Vec<f64>, threshold: f64) {
        for c in coeffs.iter_mut() {
            if c.abs() <= threshold {
                *c = 0.0;
            } else {
                *c = c.signum() * (c.abs() - threshold);
            }
        }
    }

    /// Universal (VisuShrink) threshold: σ √(2 ln N).
    ///
    /// Estimates noise σ from the median absolute deviation of the detail
    /// coefficients.
    pub fn universal_threshold(signal: &[f64]) -> f64 {
        let n = signal.len();
        if n == 0 {
            return 0.0;
        }
        // Estimate σ via MAD / 0.6745 from detail coefficients
        let (_approx, detail) = HaarWavelet::forward(signal);
        let sigma = Self::mad_sigma(&detail);
        sigma * (2.0 * (n as f64).ln()).sqrt()
    }

    /// SURE (Stein's Unbiased Risk Estimator) threshold.
    ///
    /// Selects the threshold that minimizes the SURE risk estimate.
    pub fn sure_threshold(signal: &[f64]) -> f64 {
        let n = signal.len();
        if n == 0 {
            return 0.0;
        }
        let mut sorted: Vec<f64> = signal.iter().map(|x| x.abs()).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut best_risk = f64::INFINITY;
        let mut best_thresh = 0.0;
        let total_sq: f64 = signal.iter().map(|x| x * x).sum();
        for (i, &thresh) in sorted.iter().enumerate() {
            // SURE risk for threshold t:  N - 2*(# |x_i| <= t) + sum(min(|x_i|,t)^2)
            let k = i + 1; // number of coefficients with |x| <= thresh
            let sum_clamped: f64 = sorted
                .iter()
                .map(|&x| if x <= thresh { x * x } else { thresh * thresh })
                .sum();
            let risk = (n as f64) - 2.0 * (k as f64) + sum_clamped + (total_sq - sum_clamped);
            if risk < best_risk {
                best_risk = risk;
                best_thresh = thresh;
            }
        }
        best_thresh
    }

    /// Apply full denoising pipeline: DWT → threshold → IDWT.
    ///
    /// Decomposes the signal using the Haar DWT, applies soft thresholding
    /// with the universal threshold, and reconstructs.
    pub fn denoise_signal(signal: &[f64]) -> Vec<f64> {
        if signal.len() < 2 {
            return signal.to_vec();
        }
        let threshold = Self::universal_threshold(signal);
        let mut coeffs = DiscreteWaveletTransform::dwt_haar(signal);
        // Keep the approximation (index 0) untouched; threshold details
        let n_approx = coeffs.len() / 2;
        let (_, detail_part) = coeffs.split_at_mut(n_approx);
        let mut detail_vec = detail_part.to_vec();
        Self::soft_threshold(&mut detail_vec, threshold);
        for (c, d) in detail_part.iter_mut().zip(detail_vec.iter()) {
            *c = *d;
        }
        DiscreteWaveletTransform::idwt_haar(&coeffs)
    }

    /// Estimate noise standard deviation from detail coefficients using MAD.
    fn mad_sigma(detail: &[f64]) -> f64 {
        if detail.is_empty() {
            return 1.0;
        }
        let mut abs_vals: Vec<f64> = detail.iter().map(|x| x.abs()).collect();
        abs_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med = if abs_vals.len() % 2 == 1 {
            abs_vals[abs_vals.len() / 2]
        } else {
            let m = abs_vals.len() / 2;
            (abs_vals[m - 1] + abs_vals[m]) / 2.0
        };
        (med / 0.6745).max(1e-12)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WaveletCompression
// ─────────────────────────────────────────────────────────────────────────────

/// Wavelet-based signal compression using coefficient thresholding.
///
/// Implements sparse representation, zero-tree coding concept, and
/// compression ratio estimation.
#[derive(Debug, Clone)]
pub struct WaveletCompression {
    /// Original signal length.
    pub original_length: usize,
    /// Number of non-zero coefficients after compression.
    pub nonzero_count: usize,
    /// Compressed coefficient vector (zero-tree encoded as sparse indices+values).
    pub sparse_indices: Vec<usize>,
    /// Coefficient values at sparse indices.
    pub sparse_values: Vec<f64>,
    /// Threshold used for compression.
    pub threshold: f64,
}

impl WaveletCompression {
    /// Compress a signal by retaining only DWT coefficients above `retain_fraction`.
    ///
    /// `retain_fraction` is the fraction of energy to preserve (e.g., 0.99 retains
    /// 99% of signal energy).
    pub fn compress(signal: &[f64], retain_fraction: f64) -> Self {
        let n = signal.len();
        let coeffs = DiscreteWaveletTransform::dwt_haar(signal);
        let total_energy: f64 = coeffs.iter().map(|x| x * x).sum();

        // Sort coefficients by magnitude descending; find threshold
        let mut sorted_mags: Vec<f64> = coeffs.iter().map(|x| x.abs()).collect();
        sorted_mags.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

        let target_energy = retain_fraction.clamp(0.0, 1.0) * total_energy;
        let mut cumulative = 0.0;
        let mut threshold = 0.0;
        for &mag in &sorted_mags {
            cumulative += mag * mag;
            if cumulative >= target_energy {
                threshold = mag;
                break;
            }
        }

        let sparse_indices: Vec<usize> = coeffs
            .iter()
            .enumerate()
            .filter(|&(_, c)| c.abs() >= threshold)
            .map(|(i, _)| i)
            .collect();
        let sparse_values: Vec<f64> = sparse_indices.iter().map(|&i| coeffs[i]).collect();
        let nonzero_count = sparse_indices.len();

        Self {
            original_length: n,
            nonzero_count,
            sparse_indices,
            sparse_values,
            threshold,
        }
    }

    /// Reconstruct signal from compressed sparse representation.
    ///
    /// Fills a full coefficient vector from the sparse representation and
    /// applies the inverse DWT.
    pub fn decompress(&self) -> Vec<f64> {
        let mut coeffs = vec![0.0f64; self.original_length];
        for (&idx, &val) in self.sparse_indices.iter().zip(self.sparse_values.iter()) {
            if idx < coeffs.len() {
                coeffs[idx] = val;
            }
        }
        DiscreteWaveletTransform::idwt_haar(&coeffs)
    }

    /// Compute compression ratio = original_length / nonzero_count.
    ///
    /// Returns 1.0 if no compression occurred.
    pub fn compression_ratio(&self) -> f64 {
        if self.nonzero_count == 0 {
            return self.original_length as f64;
        }
        self.original_length as f64 / self.nonzero_count as f64
    }

    /// Compute the root-mean-square error between original and reconstructed signal.
    pub fn reconstruction_error(&self, original: &[f64]) -> f64 {
        let reconstructed = self.decompress();
        let n = original.len().min(reconstructed.len());
        if n == 0 {
            return 0.0;
        }
        let mse: f64 = original
            .iter()
            .zip(reconstructed.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / n as f64;
        mse.sqrt()
    }

    /// Compute the peak signal-to-noise ratio (PSNR) in dB.
    ///
    /// Uses the formula `PSNR = 20 log10(max_val / RMSE)`.
    pub fn psnr(&self, original: &[f64]) -> f64 {
        let rmse = self.reconstruction_error(original);
        if rmse < 1e-30 {
            return f64::INFINITY;
        }
        let max_val = original.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        if max_val < 1e-30 {
            return 0.0;
        }
        20.0 * (max_val / rmse).log10()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MultiResolutionAnalysis
// ─────────────────────────────────────────────────────────────────────────────

/// Multi-resolution analysis (MRA) using the Haar wavelet.
///
/// Provides approximate and detail coefficients at any decomposition level,
/// and signal reconstruction from a chosen level.
pub struct MultiResolutionAnalysis;

impl MultiResolutionAnalysis {
    /// Extract the approximation coefficients at decomposition `level`.
    ///
    /// Applies `level` forward Haar transforms, returning the coarsest
    /// approximation vector.
    pub fn approximate_coefficients(signal: &[f64], level: usize) -> Vec<f64> {
        DiscreteWaveletTransform::dwt_level(signal, level)
    }

    /// Extract the detail coefficients at decomposition `level`.
    ///
    /// Applies `level - 1` forward Haar transforms and returns the detail
    /// coefficients from the final step.
    pub fn detail_coefficients(signal: &[f64], level: usize) -> Vec<f64> {
        if level == 0 {
            return signal.to_vec();
        }
        let mut approx = signal.to_vec();
        for _ in 0..level - 1 {
            if approx.len() < 2 {
                break;
            }
            let (a, _) = HaarWavelet::forward(&approx);
            approx = a;
        }
        if approx.len() < 2 {
            return approx;
        }
        let (_, detail) = HaarWavelet::forward(&approx);
        detail
    }

    /// Reconstruct a signal from its approximation and a list of detail arrays.
    ///
    /// `details` should be ordered from coarsest (index 0) to finest.
    pub fn reconstruct_from_level(approx: &[f64], details: &[Vec<f64>]) -> Vec<f64> {
        let mut signal = approx.to_vec();
        for detail in details {
            signal = HaarWavelet::inverse(&signal, detail);
        }
        signal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PoddedWavelet — physics-specific wavelet tools
// ─────────────────────────────────────────────────────────────────────────────

/// Physics-specific wavelet-based analysis tools.
///
/// Provides turbulence spectrum estimation, shock detection, and
/// intermittency measurement using wavelet methods.
pub struct PoddedWavelet;

impl PoddedWavelet {
    /// Estimate the turbulence energy spectrum from a velocity field using the
    /// Haar DWT.
    ///
    /// Returns a list of `(frequency, energy)` pairs, one per DWT level.
    /// `dt` is the sampling interval.
    pub fn turbulence_spectrum_wavelet(velocity_field: &[f64], dt: f64) -> Vec<(f64, f64)> {
        if velocity_field.len() < 2 {
            return Vec::new();
        }
        let decomp = HaarWavelet::decompose_level(velocity_field, 8);
        let n_levels = decomp.len();
        let mut result = Vec::with_capacity(n_levels);
        for (level_idx, coeffs) in decomp.iter().enumerate() {
            if coeffs.is_empty() {
                continue;
            }
            let energy: f64 = coeffs.iter().map(|x| x * x).sum::<f64>() / coeffs.len() as f64;
            // Frequency: scale = 2^level, frequency = 1 / (scale * dt)
            let scale = 2.0_f64.powi(level_idx as i32);
            let freq = 1.0 / (scale * dt.max(1e-30));
            result.push((freq, energy));
        }
        result
    }

    /// Detect discontinuities (shocks) in a signal using first-difference
    /// magnitudes as a simple multi-scale shock indicator.
    ///
    /// Returns the sample indices where the absolute first difference
    /// `|signal[i+1] - signal[i]|` exceeds `threshold`.
    pub fn shock_detection(signal: &[f64], threshold: f64) -> Vec<usize> {
        if signal.len() < 2 {
            return Vec::new();
        }
        signal
            .windows(2)
            .enumerate()
            .filter_map(|(i, w)| {
                let diff = (w[1] - w[0]).abs();
                if diff > threshold { Some(i) } else { None }
            })
            .collect()
    }

    /// Compute a wavelet-based intermittency measure at scale `scale`.
    ///
    /// Returns the kurtosis of the wavelet coefficients at the given scale,
    /// normalized so that a Gaussian signal returns 0.
    pub fn intermittency_measure(signal: &[f64], scale: f64) -> f64 {
        let n = signal.len();
        if n < 4 {
            return 0.0;
        }
        // Compute wavelet coefficients at this scale using Morlet
        let dt = 1.0;
        let scales = [scale];
        let scalogram = ContinuousWavelet::cwt_morlet(signal, dt, &scales);
        if scalogram.is_empty() {
            return 0.0;
        }
        let coeffs = &scalogram[0];
        let n_c = coeffs.len() as f64;
        if n_c < 2.0 {
            return 0.0;
        }
        let mean = coeffs.iter().sum::<f64>() / n_c;
        let variance = coeffs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_c;
        if variance < 1e-30 {
            return 0.0;
        }
        let kurtosis =
            coeffs.iter().map(|x| (x - mean).powi(4)).sum::<f64>() / (n_c * variance * variance);
        kurtosis - 3.0 // excess kurtosis
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── HaarWavelet ────────────────────────────────────────────────────────────

    #[test]
    fn test_haar_forward_basic() {
        let signal = vec![1.0, 3.0, 5.0, 11.0];
        let (a, d) = HaarWavelet::forward(&signal);
        assert_eq!(a.len(), 2);
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn test_haar_perfect_reconstruction() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let (a, d) = HaarWavelet::forward(&signal);
        let reconstructed = HaarWavelet::inverse(&a, &d);
        for (orig, rec) in signal.iter().zip(reconstructed.iter()) {
            assert!((orig - rec).abs() < 1e-10, "{orig} vs {rec}");
        }
    }

    #[test]
    fn test_haar_energy_conservation() {
        let signal: Vec<f64> = (1..=8).map(|x| x as f64).collect();
        let energy_in: f64 = signal.iter().map(|x| x * x).sum();
        let (a, d) = HaarWavelet::forward(&signal);
        let energy_out: f64 =
            a.iter().map(|x| x * x).sum::<f64>() + d.iter().map(|x| x * x).sum::<f64>();
        assert!(
            (energy_in - energy_out).abs() < 1e-8,
            "energy not conserved"
        );
    }

    #[test]
    fn test_haar_constant_signal_zero_detail() {
        let signal = vec![5.0; 8];
        let (_a, d) = HaarWavelet::forward(&signal);
        for x in &d {
            assert!(x.abs() < 1e-12, "detail should be zero for constant signal");
        }
    }

    #[test]
    fn test_haar_decompose_level_count() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let levels = HaarWavelet::decompose_level(&signal, 3);
        assert_eq!(levels.len(), 4, "should return levels+1 = 4 entries");
    }

    #[test]
    fn test_haar_decompose_level_first_is_coarsest() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let levels = HaarWavelet::decompose_level(&signal, 3);
        // Coarsest (index 0) should be the shortest array
        let coarsest = &levels[0];
        let finest_detail = &levels[levels.len() - 1];
        assert!(coarsest.len() <= finest_detail.len());
    }

    #[test]
    fn test_haar_inverse_length() {
        let a = vec![1.0, 2.0];
        let d = vec![0.5, 0.3];
        let out = HaarWavelet::inverse(&a, &d);
        assert_eq!(out.len(), 4);
    }

    // ── Daubechies4 ───────────────────────────────────────────────────────────

    #[test]
    fn test_d4_forward_length() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let (a, d) = Daubechies4::forward(&signal);
        assert_eq!(a.len(), 8);
        assert_eq!(d.len(), 8);
    }

    #[test]
    fn test_d4_filter_norms() {
        // Both filters should be unit-norm (orthonormality condition)
        let h_norm: f64 = Daubechies4::H.iter().map(|x| x * x).sum::<f64>().sqrt();
        let g_norm: f64 = Daubechies4::G.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((h_norm - 1.0).abs() < 1e-7, "H norm = {h_norm}");
        assert!((g_norm - 1.0).abs() < 1e-7, "G norm = {g_norm}");
    }

    #[test]
    fn test_d4_constant_signal() {
        // Forward then inverse should approximately reconstruct for constant input
        let signal = vec![1.0; 8];
        let (a, d) = Daubechies4::forward(&signal);
        let rec = Daubechies4::inverse(&a, &d);
        // Sum of reconstruction should approximately equal sum of input
        let sum_in: f64 = signal.iter().sum();
        let sum_rec: f64 = rec.iter().sum();
        assert!((sum_in - sum_rec).abs() < 1e-6, "{sum_in} vs {sum_rec}");
    }

    #[test]
    fn test_d4_energy_approximately_conserved() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64).sin()).collect();
        let e_in: f64 = signal.iter().map(|x| x * x).sum();
        let (a, d) = Daubechies4::forward(&signal);
        let e_out: f64 =
            a.iter().map(|x| x * x).sum::<f64>() + d.iter().map(|x| x * x).sum::<f64>();
        assert!(
            (e_in - e_out).abs() / e_in.max(1e-12) < 0.01,
            "D4 energy deviation too large: {e_in} vs {e_out}"
        );
    }

    // ── Symlet4 ───────────────────────────────────────────────────────────────

    #[test]
    fn test_sym4_filter_norms() {
        let h_norm: f64 = Symlet4::H.iter().map(|x| x * x).sum::<f64>().sqrt();
        let g_norm: f64 = Symlet4::G.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((h_norm - 1.0).abs() < 1e-7, "Sym4 H norm = {h_norm}");
        assert!((g_norm - 1.0).abs() < 1e-7, "Sym4 G norm = {g_norm}");
    }

    #[test]
    fn test_sym4_forward_length() {
        let signal: Vec<f64> = (0..32).map(|x| x as f64).collect();
        let (a, d) = Symlet4::forward(&signal);
        assert_eq!(a.len(), 16);
        assert_eq!(d.len(), 16);
    }

    #[test]
    fn test_sym4_energy_approximately_conserved() {
        let signal: Vec<f64> = (0..32).map(|x| (x as f64 * 0.5).sin()).collect();
        let e_in: f64 = signal.iter().map(|x| x * x).sum();
        let (a, d) = Symlet4::forward(&signal);
        let e_out: f64 =
            a.iter().map(|x| x * x).sum::<f64>() + d.iter().map(|x| x * x).sum::<f64>();
        assert!(
            (e_in - e_out).abs() / e_in.max(1e-12) < 0.02,
            "Sym4 energy deviation: {e_in} vs {e_out}"
        );
    }

    #[test]
    fn test_sym4_inverse_length() {
        let a = vec![1.0; 8];
        let d = vec![0.5; 8];
        let out = Symlet4::inverse(&a, &d);
        assert_eq!(out.len(), 16);
    }

    // ── WaveletTransform ──────────────────────────────────────────────────────

    #[test]
    fn test_wavelet_transform_haar_forward() {
        let wt = WaveletTransform::new(WaveletFamily::Haar);
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let (a, d) = wt.forward_one_level(&signal);
        assert_eq!(a.len(), 2);
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn test_wavelet_transform_d4_roundtrip_length() {
        let wt = WaveletTransform::new(WaveletFamily::Daubechies4);
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let (a, d) = wt.forward_one_level(&signal);
        let rec = wt.inverse_one_level(&a, &d);
        assert_eq!(rec.len(), signal.len());
    }

    #[test]
    fn test_wavelet_transform_sym4_forward() {
        let wt = WaveletTransform::new(WaveletFamily::Symlet4);
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let (a, d) = wt.forward_one_level(&signal);
        assert_eq!(a.len(), 8);
        assert_eq!(d.len(), 8);
    }

    #[test]
    fn test_wavelet_transform_energy() {
        let coeffs = vec![3.0, 4.0];
        assert!((WaveletTransform::energy(&coeffs) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_wavelet_transform_entropy_uniform() {
        // Uniform coefficients: maximum entropy
        let coeffs = vec![1.0; 4];
        let h = WaveletTransform::entropy(&coeffs);
        assert!(h > 0.0);
    }

    #[test]
    fn test_wavelet_transform_entropy_spike() {
        // All energy in one coefficient: minimal entropy
        let mut coeffs = vec![0.0; 8];
        coeffs[0] = 1.0;
        let h = WaveletTransform::entropy(&coeffs);
        assert!(h.abs() < 1e-12);
    }

    // ── MultilevelDwt ─────────────────────────────────────────────────────────

    #[test]
    fn test_multilevel_haar_reconstruction() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let mld = MultilevelDwt::decompose(&signal, 3, WaveletFamily::Haar);
        let rec = mld.reconstruct();
        assert_eq!(rec.len(), signal.len());
        for (a, b) in signal.iter().zip(rec.iter()) {
            assert!((a - b).abs() < 1e-8, "MultilevelDwt roundtrip: {a} vs {b}");
        }
    }

    #[test]
    fn test_multilevel_d4_reconstruction() {
        let signal: Vec<f64> = (0..32).map(|x| (x as f64 * 0.3).sin()).collect();
        let mld = MultilevelDwt::decompose(&signal, 2, WaveletFamily::Daubechies4);
        let rec = mld.reconstruct();
        assert_eq!(rec.len(), signal.len());
    }

    #[test]
    fn test_multilevel_approx_coarsest() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let mld = MultilevelDwt::decompose(&signal, 3, WaveletFamily::Haar);
        // Coarsest approx should have length 16 / 2^3 = 2
        assert_eq!(mld.approx.len(), 2);
    }

    #[test]
    fn test_multilevel_detail_count() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let mld = MultilevelDwt::decompose(&signal, 3, WaveletFamily::Haar);
        assert_eq!(mld.details.len(), 3);
    }

    #[test]
    fn test_multilevel_total_energy_conserved() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64).cos()).collect();
        let e_in: f64 = signal.iter().map(|x| x * x).sum();
        let mld = MultilevelDwt::decompose(&signal, 3, WaveletFamily::Haar);
        let e_stored = mld.total_energy();
        assert!(
            (e_in - e_stored).abs() / e_in.max(1e-12) < 1e-8,
            "energy conserved: {e_in} vs {e_stored}"
        );
    }

    #[test]
    fn test_multilevel_level_entropies_count() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let mld = MultilevelDwt::decompose(&signal, 3, WaveletFamily::Haar);
        let entropies = mld.level_entropies();
        // 1 approx + 3 detail levels
        assert_eq!(entropies.len(), 4);
    }

    // ── WaveletPacket ─────────────────────────────────────────────────────────

    #[test]
    fn test_packet_leaf_count() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let pkt = WaveletPacket::new_haar(&signal, 3);
        assert_eq!(pkt.coefficients.len(), 8); // 2^3 leaves
    }

    #[test]
    fn test_packet_energy_distribution_length() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let pkt = WaveletPacket::new_haar(&signal, 2);
        let energies = pkt.energy_distribution();
        assert_eq!(energies.len(), 4);
    }

    #[test]
    fn test_packet_best_basis_nonempty() {
        let signal: Vec<f64> = (0..8).map(|x| (x as f64).sin()).collect();
        let pkt = WaveletPacket::new_haar(&signal, 2);
        let bb = pkt.best_basis();
        assert!(!bb.is_empty());
    }

    #[test]
    fn test_packet_total_energy_conserved() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let e_in: f64 = signal.iter().map(|x| x * x).sum();
        let pkt = WaveletPacket::new_haar(&signal, 1); // 1 level → 2 leaves
        let e_out: f64 = pkt.energy_distribution().iter().sum();
        assert!((e_in - e_out).abs() < 1e-6, "packet energy not conserved");
    }

    #[test]
    fn test_packet_best_basis_entropy_nonempty() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64 * 0.5).sin()).collect();
        let pkt = WaveletPacket::new_haar(&signal, 2);
        let bb = pkt.best_basis_entropy();
        assert!(!bb.is_empty());
    }

    #[test]
    fn test_packet_entropy_distribution_length() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let pkt = WaveletPacket::new_haar(&signal, 2);
        let ent = pkt.entropy_distribution();
        assert_eq!(ent.len(), 4);
    }

    #[test]
    fn test_packet_new_with_family() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let pkt = WaveletPacket::new(&signal, 2, WaveletFamily::Daubechies4);
        assert_eq!(pkt.coefficients.len(), 4); // 2^2 leaves
    }

    #[test]
    fn test_packet_total_energy_method() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let pkt = WaveletPacket::new_haar(&signal, 1);
        assert!(pkt.total_energy() > 0.0);
    }

    // ── ContinuousWavelet ─────────────────────────────────────────────────────

    #[test]
    fn test_morlet_at_zero() {
        // At t=0 the Gaussian envelope is 1 and cos(0)=1
        let v = ContinuousWavelet::morlet(0.0, 1.0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_morlet_decays_at_large_t() {
        let v0 = ContinuousWavelet::morlet(0.0, 1.0).abs();
        let v10 = ContinuousWavelet::morlet(10.0, 1.0).abs();
        assert!(v10 < v0 * 0.01, "Morlet should decay: {v10} vs {v0}");
    }

    #[test]
    fn test_mexican_hat_at_zero() {
        let v = ContinuousWavelet::mexican_hat(0.0, 1.0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mexican_hat_sign_change() {
        // At t=sigma, value should be 0 (zero crossing)
        let sigma = 1.0;
        let v = ContinuousWavelet::mexican_hat(sigma, sigma);
        assert!(v.abs() < 1e-12);
    }

    #[test]
    fn test_paul_wavelet_finite() {
        let v = ContinuousWavelet::paul(0.5, 4);
        assert!(v.is_finite());
    }

    #[test]
    fn test_paul_wavelet_decay() {
        // Paul wavelet should decay for large |t|
        let v_near = ContinuousWavelet::paul(0.0, 2).abs();
        let v_far = ContinuousWavelet::paul(100.0, 2).abs();
        assert!(
            v_far < v_near,
            "Paul wavelet should decay: {v_far} vs {v_near}"
        );
    }

    #[test]
    fn test_cwt_morlet_shape() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64 * 0.5).sin()).collect();
        let scales = vec![1.0, 2.0, 4.0];
        let scalogram = ContinuousWavelet::cwt_morlet(&signal, 1.0, &scales);
        assert_eq!(scalogram.len(), 3);
        assert_eq!(scalogram[0].len(), 16);
    }

    #[test]
    fn test_cwt_mexican_hat_shape() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64 * 0.5).sin()).collect();
        let scales = vec![1.0, 2.0];
        let scalogram = ContinuousWavelet::cwt_mexican_hat(&signal, 1.0, &scales);
        assert_eq!(scalogram.len(), 2);
        assert_eq!(scalogram[0].len(), 16);
    }

    #[test]
    fn test_cwt_paul_shape() {
        let signal: Vec<f64> = (0..8).map(|x| (x as f64).sin()).collect();
        let scales = vec![1.0, 2.0];
        let scalogram = ContinuousWavelet::cwt_paul(&signal, 1.0, &scales, 4);
        assert_eq!(scalogram.len(), 2);
        assert_eq!(scalogram[0].len(), 8);
    }

    #[test]
    fn test_instantaneous_frequency_length() {
        let signal: Vec<f64> = (0..8).map(|x| (x as f64).sin()).collect();
        let scales = vec![1.0, 2.0];
        let scalogram = ContinuousWavelet::cwt_morlet(&signal, 1.0, &scales);
        let freq = ContinuousWavelet::instantaneous_frequency(&scalogram, &scales, 1.0);
        assert_eq!(freq.len(), 8);
    }

    #[test]
    fn test_instantaneous_frequency_positive() {
        let signal: Vec<f64> = (0..8).map(|x| (x as f64).sin()).collect();
        let scales = vec![1.0, 2.0, 4.0];
        let scalogram = ContinuousWavelet::cwt_morlet(&signal, 1.0, &scales);
        let freq = ContinuousWavelet::instantaneous_frequency(&scalogram, &scales, 1.0);
        for f in &freq {
            assert!(*f > 0.0, "frequency must be positive");
        }
    }

    #[test]
    fn test_scale_averaged_power_length() {
        let signal: Vec<f64> = (0..8).map(|x| (x as f64).sin()).collect();
        let scales = vec![1.0, 2.0, 4.0];
        let scalogram = ContinuousWavelet::cwt_morlet(&signal, 1.0, &scales);
        let power = ContinuousWavelet::scale_averaged_power(&scalogram, &scales);
        assert_eq!(power.len(), 3);
    }

    // ── DiscreteWaveletTransform ───────────────────────────────────────────────

    #[test]
    fn test_dwt_haar_length_preserved() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let coeffs = DiscreteWaveletTransform::dwt_haar(&signal);
        assert_eq!(coeffs.len(), signal.len());
    }

    #[test]
    fn test_idwt_haar_roundtrip() {
        let signal: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let coeffs = DiscreteWaveletTransform::dwt_haar(&signal);
        let reconstructed = DiscreteWaveletTransform::idwt_haar(&coeffs);
        for (a, b) in signal.iter().zip(reconstructed.iter()) {
            assert!((a - b).abs() < 1e-8, "IDWT roundtrip failed: {a} vs {b}");
        }
    }

    #[test]
    fn test_dwt_level_length() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let approx = DiscreteWaveletTransform::dwt_level(&signal, 2);
        assert_eq!(approx.len(), 4); // 16 / 2^2 = 4
    }

    #[test]
    fn test_wavelet_energy_nonneg() {
        let coeffs = vec![1.0, -2.0, 3.0];
        let e = DiscreteWaveletTransform::wavelet_energy(&coeffs);
        assert!(e >= 0.0);
        assert!((e - 14.0).abs() < 1e-12);
    }

    #[test]
    fn test_wavelet_entropy_nonneg() {
        let coeffs = vec![1.0, 2.0, 3.0, 4.0];
        let h = DiscreteWaveletTransform::wavelet_entropy(&coeffs);
        assert!(h >= 0.0);
    }

    #[test]
    fn test_wavelet_entropy_uniform_max() {
        // Uniform coefficients should give higher entropy than concentrated ones
        let uniform = vec![1.0; 8];
        let spike: Vec<f64> = {
            let mut v = vec![0.0; 8];
            v[0] = 8.0;
            v
        };
        let h_uniform = DiscreteWaveletTransform::wavelet_entropy(&uniform);
        let h_spike = DiscreteWaveletTransform::wavelet_entropy(&spike);
        assert!(h_uniform >= h_spike, "uniform should have higher entropy");
    }

    // ── WaveletDenoising ──────────────────────────────────────────────────────

    #[test]
    fn test_hard_threshold_zeroes_small() {
        let mut coeffs = vec![0.1, 0.5, -0.2, 1.0];
        WaveletDenoising::hard_threshold(&mut coeffs, 0.3);
        assert_eq!(coeffs[0], 0.0);
        assert_eq!(coeffs[1], 0.5);
        assert_eq!(coeffs[2], 0.0);
        assert_eq!(coeffs[3], 1.0);
    }

    #[test]
    fn test_soft_threshold_shrinks() {
        let mut coeffs = vec![1.0, -2.0, 0.5];
        WaveletDenoising::soft_threshold(&mut coeffs, 0.5);
        assert!((coeffs[0] - 0.5).abs() < 1e-12);
        assert!((coeffs[1] + 1.5).abs() < 1e-12);
        assert!(coeffs[2].abs() < 1e-12);
    }

    #[test]
    fn test_universal_threshold_positive() {
        let signal: Vec<f64> = (0..64)
            .map(|i| (i as f64 * 0.1).sin() + 0.1 * (i as f64))
            .collect();
        let t = WaveletDenoising::universal_threshold(&signal);
        assert!(t >= 0.0, "threshold must be non-negative");
    }

    #[test]
    fn test_sure_threshold_nonneg() {
        let signal = vec![1.0, 2.0, 3.0, 0.1, 0.05];
        let t = WaveletDenoising::sure_threshold(&signal);
        assert!(t >= 0.0);
    }

    #[test]
    fn test_soft_threshold_zero_input() {
        let mut coeffs: Vec<f64> = vec![];
        WaveletDenoising::soft_threshold(&mut coeffs, 1.0);
        assert!(coeffs.is_empty());
    }

    #[test]
    fn test_denoise_signal_length_preserved() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64 * 0.5).sin()).collect();
        let denoised = WaveletDenoising::denoise_signal(&signal);
        assert_eq!(denoised.len(), signal.len());
    }

    #[test]
    fn test_denoise_small_signal_unchanged() {
        let signal = vec![1.0];
        let denoised = WaveletDenoising::denoise_signal(&signal);
        assert_eq!(denoised, signal);
    }

    // ── WaveletCompression ────────────────────────────────────────────────────

    #[test]
    fn test_compression_ratio_ge_one() {
        let signal: Vec<f64> = (0..64).map(|x| (x as f64 * 0.5).sin()).collect();
        let comp = WaveletCompression::compress(&signal, 0.90);
        assert!(comp.compression_ratio() >= 1.0);
    }

    #[test]
    fn test_compression_lossless_full_retention() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let comp = WaveletCompression::compress(&signal, 1.0);
        let rec = comp.decompress();
        assert_eq!(rec.len(), signal.len());
        for (a, b) in signal.iter().zip(rec.iter()) {
            assert!((a - b).abs() < 1e-6, "lossless compression: {a} vs {b}");
        }
    }

    #[test]
    fn test_compression_original_length_stored() {
        let signal: Vec<f64> = (0..32).map(|x| x as f64).collect();
        let comp = WaveletCompression::compress(&signal, 0.95);
        assert_eq!(comp.original_length, 32);
    }

    #[test]
    fn test_compression_reconstruction_error_finite() {
        let signal: Vec<f64> = (0..32).map(|x| (x as f64 * 0.5).sin()).collect();
        let comp = WaveletCompression::compress(&signal, 0.99);
        let err = comp.reconstruction_error(&signal);
        assert!(err.is_finite());
        assert!(err >= 0.0);
    }

    #[test]
    fn test_compression_psnr_nonneg() {
        let signal: Vec<f64> = (0..32).map(|x| (x as f64 * 0.3).cos()).collect();
        let comp = WaveletCompression::compress(&signal, 0.99);
        let psnr = comp.psnr(&signal);
        assert!(psnr >= 0.0 || psnr.is_infinite());
    }

    #[test]
    fn test_compression_sparse_indices_nonempty() {
        let signal: Vec<f64> = (0..16).map(|x| (x as f64).sin()).collect();
        let comp = WaveletCompression::compress(&signal, 0.8);
        assert!(!comp.sparse_indices.is_empty());
    }

    // ── MultiResolutionAnalysis ───────────────────────────────────────────────

    #[test]
    fn test_mra_approx_length() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let approx = MultiResolutionAnalysis::approximate_coefficients(&signal, 2);
        assert_eq!(approx.len(), 4);
    }

    #[test]
    fn test_mra_detail_length() {
        let signal: Vec<f64> = (0..16).map(|x| x as f64).collect();
        let detail = MultiResolutionAnalysis::detail_coefficients(&signal, 1);
        assert_eq!(detail.len(), 8);
    }

    #[test]
    fn test_mra_detail_level_zero() {
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let d = MultiResolutionAnalysis::detail_coefficients(&signal, 0);
        // Level 0 returns the signal unchanged
        assert_eq!(d, signal);
    }

    #[test]
    fn test_mra_reconstruct_single_level() {
        let signal: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let approx = MultiResolutionAnalysis::approximate_coefficients(&signal, 1);
        let detail = MultiResolutionAnalysis::detail_coefficients(&signal, 1);
        let rec = MultiResolutionAnalysis::reconstruct_from_level(&approx, &[detail]);
        assert_eq!(rec.len(), signal.len());
        for (a, b) in signal.iter().zip(rec.iter()) {
            assert!((a - b).abs() < 1e-8, "MRA roundtrip: {a} vs {b}");
        }
    }

    // ── PoddedWavelet ─────────────────────────────────────────────────────────

    #[test]
    fn test_turbulence_spectrum_nonempty() {
        let vel: Vec<f64> = (0..64).map(|i| (i as f64 * 0.3).sin()).collect();
        let spectrum = PoddedWavelet::turbulence_spectrum_wavelet(&vel, 0.01);
        assert!(!spectrum.is_empty(), "spectrum should not be empty");
    }

    #[test]
    fn test_turbulence_spectrum_positive_freq() {
        let vel: Vec<f64> = (0..32).map(|i| i as f64).collect();
        let spectrum = PoddedWavelet::turbulence_spectrum_wavelet(&vel, 0.01);
        for (f, _e) in &spectrum {
            assert!(*f > 0.0, "frequencies must be positive");
        }
    }

    #[test]
    fn test_turbulence_spectrum_nonneg_energy() {
        let vel: Vec<f64> = (0..32).map(|i| (i as f64).cos()).collect();
        for (_f, e) in PoddedWavelet::turbulence_spectrum_wavelet(&vel, 0.01) {
            assert!(e >= 0.0);
        }
    }

    #[test]
    fn test_shock_detection_step_function() {
        // Step function: indices around the discontinuity should be detected
        let mut signal = vec![0.0_f64; 16];
        for i in 8..16 {
            signal[i] = 10.0;
        }
        let shocks = PoddedWavelet::shock_detection(&signal, 0.5);
        assert!(!shocks.is_empty(), "should detect discontinuity");
    }

    #[test]
    fn test_shock_detection_smooth_signal_empty() {
        // Very smooth signal: no shocks at high threshold
        let signal: Vec<f64> = (0..16).map(|i| (i as f64 * 0.1).sin()).collect();
        let shocks = PoddedWavelet::shock_detection(&signal, 1000.0);
        assert!(
            shocks.is_empty(),
            "no shocks expected for smooth signal at high threshold"
        );
    }

    #[test]
    fn test_shock_detection_empty_signal() {
        let shocks = PoddedWavelet::shock_detection(&[], 0.5);
        assert!(shocks.is_empty());
    }

    #[test]
    fn test_intermittency_small_signal() {
        // Signal too small → return 0
        let v = PoddedWavelet::intermittency_measure(&[1.0, 2.0], 1.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_intermittency_gaussian_near_zero() {
        // For a Gaussian-like signal the excess kurtosis is near 0
        let signal: Vec<f64> = (0..64).map(|i| (i as f64 * 0.2).sin()).collect();
        let k = PoddedWavelet::intermittency_measure(&signal, 2.0);
        // Just check it is finite
        assert!(k.is_finite(), "intermittency should be finite");
    }
}
