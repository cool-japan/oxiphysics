#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Signal analysis: autocorrelation, cross-correlation, PSD, Hilbert transform,
//! empirical mode decomposition, peak detection, and signal filtering.

#![allow(dead_code)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Helper: FFT (Cooley-Tukey radix-2)
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct C64 {
    re: f64,
    im: f64,
}

impl C64 {
    #[inline]
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    #[inline]
    fn exp_i(theta: f64) -> Self {
        Self {
            re: theta.cos(),
            im: theta.sin(),
        }
    }
}

impl std::ops::Add for C64 {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        Self::new(self.re + r.re, self.im + r.im)
    }
}
impl std::ops::Sub for C64 {
    type Output = Self;
    fn sub(self, r: Self) -> Self {
        Self::new(self.re - r.re, self.im - r.im)
    }
}
impl std::ops::Mul for C64 {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self::new(
            self.re * r.re - self.im * r.im,
            self.re * r.im + self.im * r.re,
        )
    }
}

fn fft(a: &mut Vec<C64>, invert: bool) {
    let n = a.len();
    // Bit-reverse permutation
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2usize;
    while len <= n {
        let ang = if invert {
            2.0 * PI / len as f64
        } else {
            -2.0 * PI / len as f64
        };
        let wlen = C64::exp_i(ang);
        let mut i = 0;
        while i < n {
            let mut w = C64::new(1.0, 0.0);
            for jj in 0..(len / 2) {
                let u = a[i + jj];
                let v = a[i + jj + len / 2] * w;
                a[i + jj] = u + v;
                a[i + jj + len / 2] = u - v;
                w = w * wlen;
            }
            i += len;
        }
        len <<= 1;
    }
    if invert {
        let nf = n as f64;
        for x in a.iter_mut() {
            x.re /= nf;
            x.im /= nf;
        }
    }
}

/// Next power of two >= n.
fn next_pow2(n: usize) -> usize {
    let mut p = 1;
    while p < n {
        p <<= 1;
    }
    p
}

// ─────────────────────────────────────────────────────────────────────────────
// AutoCorrelation
// ─────────────────────────────────────────────────────────────────────────────

/// Autocorrelation estimator for a real-valued signal.
///
/// Supports biased and unbiased estimates and optional normalization.
#[derive(Debug, Clone)]
pub struct AutoCorrelation {
    /// The autocorrelation values at lags 0, 1, ..., max_lag.
    pub values: Vec<f64>,
    /// Maximum lag computed.
    pub max_lag: usize,
    /// Whether the unbiased (1/(N-k)) estimator was used.
    pub unbiased: bool,
}

impl AutoCorrelation {
    /// Compute autocorrelation of `signal` up to `max_lag`.
    ///
    /// If `unbiased` is true, use the `1/(N-k)` normalization; otherwise use `1/N`.
    /// If `normalize` is true, divide all values by the lag-0 value.
    pub fn compute(signal: &[f64], max_lag: usize, unbiased: bool, normalize: bool) -> Self {
        let n = signal.len();
        let mean = signal.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = signal.iter().map(|&x| x - mean).collect();
        let lag_limit = max_lag.min(n.saturating_sub(1));
        let mut values = vec![0.0_f64; lag_limit + 1];
        for k in 0..=lag_limit {
            let sum: f64 = (0..(n - k)).map(|i| centered[i] * centered[i + k]).sum();
            let denom = if unbiased { (n - k) as f64 } else { n as f64 };
            values[k] = sum / denom;
        }
        if normalize && values[0].abs() > 1e-15 {
            let r0 = values[0];
            for v in values.iter_mut() {
                *v /= r0;
            }
        }
        Self {
            values,
            max_lag: lag_limit,
            unbiased,
        }
    }

    /// Return the autocorrelation at lag `k`, or 0 if out of range.
    pub fn at_lag(&self, k: usize) -> f64 {
        self.values.get(k).copied().unwrap_or(0.0)
    }

    /// Return the dominant lag (lag > 0 with maximum autocorrelation).
    pub fn dominant_lag(&self) -> usize {
        self.values[1..]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i + 1)
            .unwrap_or(0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CrossCorrelation
// ─────────────────────────────────────────────────────────────────────────────

/// Cross-correlation between two real-valued signals.
///
/// Computes R_xy(k) = Σ x\[n\] * y\[n+k\] for k in \[-max_lag, max_lag\].
#[derive(Debug, Clone)]
pub struct CrossCorrelation {
    /// Cross-correlation values; index 0 corresponds to lag -max_lag.
    pub values: Vec<f64>,
    /// Maximum lag in each direction.
    pub max_lag: usize,
}

impl CrossCorrelation {
    /// Compute cross-correlation between `x` and `y` up to `max_lag`.
    ///
    /// Both signals are mean-subtracted before computation.
    pub fn compute(x: &[f64], y: &[f64], max_lag: usize) -> Self {
        let n = x.len().min(y.len());
        let mx = x[..n].iter().sum::<f64>() / n as f64;
        let my = y[..n].iter().sum::<f64>() / n as f64;
        let lag_limit = max_lag.min(n.saturating_sub(1));
        let total_lags = 2 * lag_limit + 1;
        let mut values = vec![0.0_f64; total_lags];
        for (idx, lag) in (-(lag_limit as i64)..=(lag_limit as i64)).enumerate() {
            let sum: f64 = (0..n)
                .filter_map(|i| {
                    let j = i as i64 + lag;
                    if j >= 0 && (j as usize) < n {
                        Some((x[i] - mx) * (y[j as usize] - my))
                    } else {
                        None
                    }
                })
                .sum();
            values[idx] = sum / n as f64;
        }
        Self {
            values,
            max_lag: lag_limit,
        }
    }

    /// Return the lag at which cross-correlation is maximum.
    ///
    /// Negative lags mean x leads y; positive lags mean y leads x.
    pub fn peak_lag(&self) -> i64 {
        let (idx, _) = self
            .values
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.abs()
                    .partial_cmp(&b.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or((self.max_lag, &0.0));
        idx as i64 - self.max_lag as i64
    }

    /// Return the cross-correlation value at a given lag.
    pub fn at_lag(&self, lag: i64) -> f64 {
        let idx = lag + self.max_lag as i64;
        if idx < 0 || idx as usize >= self.values.len() {
            return 0.0;
        }
        self.values[idx as usize]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PowerSpectralDensity
// ─────────────────────────────────────────────────────────────────────────────

/// Power Spectral Density estimate via Welch's method.
///
/// The signal is divided into overlapping segments, each windowed with a Hann
/// window and FFT'd, and the resulting periodograms are averaged.
#[derive(Debug, Clone)]
pub struct PowerSpectralDensity {
    /// Frequency bins (Hz).
    pub frequencies: Vec<f64>,
    /// PSD values at each frequency bin.
    pub power: Vec<f64>,
    /// Sampling frequency (Hz).
    pub fs: f64,
    /// FFT window length used.
    pub window_len: usize,
}

impl PowerSpectralDensity {
    /// Estimate PSD using Welch's method.
    ///
    /// `fs` is the sampling frequency.
    /// `window_len` is the segment length (will be rounded to power-of-two).
    /// `overlap` in \[0, 1) fraction of the window to overlap.
    pub fn welch(signal: &[f64], fs: f64, window_len: usize, overlap: f64) -> Self {
        let n = signal.len();
        let seg_len = next_pow2(window_len.min(n));
        let step = ((seg_len as f64 * (1.0 - overlap.clamp(0.0, 0.99))) as usize).max(1);
        let hann: Vec<f64> = (0..seg_len)
            .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (seg_len - 1) as f64).cos()))
            .collect();
        let win_power: f64 = hann.iter().map(|&w| w * w).sum::<f64>();

        let half = seg_len / 2 + 1;
        let mut psd_sum = vec![0.0_f64; half];
        let mut n_segs = 0usize;

        let mut start = 0;
        while start + seg_len <= n {
            let seg = &signal[start..start + seg_len];
            let mut buf: Vec<C64> = seg
                .iter()
                .zip(hann.iter())
                .map(|(&x, &w)| C64::new(x * w, 0.0))
                .collect();
            fft(&mut buf, false);
            for k in 0..half {
                let power = buf[k].re * buf[k].re + buf[k].im * buf[k].im;
                // Double for all bins except DC and Nyquist
                psd_sum[k] += if k == 0 || k == half - 1 {
                    power
                } else {
                    2.0 * power
                };
            }
            n_segs += 1;
            start += step;
        }
        if n_segs == 0 {
            n_segs = 1;
        }
        let scale = 1.0 / (fs * win_power * n_segs as f64);
        let power: Vec<f64> = psd_sum.iter().map(|&p| p * scale).collect();
        let frequencies: Vec<f64> = (0..half).map(|k| k as f64 * fs / seg_len as f64).collect();
        Self {
            frequencies,
            power,
            fs,
            window_len: seg_len,
        }
    }

    /// Return the frequency resolution (Hz per bin).
    pub fn freq_resolution(&self) -> f64 {
        if self.window_len > 0 {
            self.fs / self.window_len as f64
        } else {
            0.0
        }
    }

    /// Return the dominant frequency (frequency bin with maximum power).
    pub fn dominant_frequency(&self) -> f64 {
        self.frequencies
            .iter()
            .zip(self.power.iter())
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(&f, _)| f)
            .unwrap_or(0.0)
    }

    /// Compute total power (integral of PSD over all frequencies).
    pub fn total_power(&self) -> f64 {
        let df = self.freq_resolution();
        self.power.iter().sum::<f64>() * df
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HilbertTransform
// ─────────────────────────────────────────────────────────────────────────────

/// Hilbert transform of a real signal, yielding the analytic signal.
///
/// The analytic signal is z\[n\] = x\[n\] + i*H{x}\[n\], where H{x} is the
/// Hilbert transform.  From this, instantaneous amplitude, phase, and
/// frequency can be extracted.
#[derive(Debug, Clone)]
pub struct HilbertTransform {
    /// Real part of the analytic signal (original signal).
    pub real_part: Vec<f64>,
    /// Imaginary part of the analytic signal (Hilbert transform).
    pub imaginary_part: Vec<f64>,
    /// Instantaneous amplitude (envelope).
    pub amplitude: Vec<f64>,
    /// Instantaneous phase in radians.
    pub phase: Vec<f64>,
}

impl HilbertTransform {
    /// Compute the Hilbert transform of `signal` using FFT.
    ///
    /// Returns the analytic signal components.
    pub fn compute(signal: &[f64]) -> Self {
        let n = signal.len();
        let n_fft = next_pow2(n);
        let mut buf: Vec<C64> = signal
            .iter()
            .map(|&x| C64::new(x, 0.0))
            .chain(std::iter::repeat(C64::new(0.0, 0.0)))
            .take(n_fft)
            .collect();
        fft(&mut buf, false);

        // Apply the one-sided filter: H[k] = 2 for 1..N/2, 1 at DC and Nyquist, 0 otherwise
        for k in 0..n_fft {
            if k == 0 || (n_fft.is_multiple_of(2) && k == n_fft / 2) {
                // keep as-is
            } else if k < n_fft / 2 {
                buf[k].re *= 2.0;
                buf[k].im *= 2.0;
            } else {
                buf[k] = C64::new(0.0, 0.0);
            }
        }
        fft(&mut buf, true); // IFFT

        let real_part: Vec<f64> = buf[..n].iter().map(|c| c.re).collect();
        let imaginary_part: Vec<f64> = buf[..n].iter().map(|c| c.im).collect();
        let amplitude: Vec<f64> = real_part
            .iter()
            .zip(imaginary_part.iter())
            .map(|(&r, &i)| (r * r + i * i).sqrt())
            .collect();
        let phase: Vec<f64> = real_part
            .iter()
            .zip(imaginary_part.iter())
            .map(|(&r, &i)| i.atan2(r))
            .collect();
        Self {
            real_part,
            imaginary_part,
            amplitude,
            phase,
        }
    }

    /// Compute instantaneous frequency from the phase using finite differences.
    ///
    /// `fs` is the sampling frequency in Hz.
    /// Phase unwrapping is applied before differencing.
    pub fn instantaneous_frequency(&self, fs: f64) -> Vec<f64> {
        let n = self.phase.len();
        if n < 2 {
            return vec![0.0; n];
        }
        // Unwrap phase
        let mut unwrapped = self.phase.clone();
        for i in 1..n {
            let diff = unwrapped[i] - unwrapped[i - 1];
            if diff > PI {
                for v in unwrapped[i..].iter_mut() {
                    *v -= 2.0 * PI;
                }
            } else if diff < -PI {
                for v in unwrapped[i..].iter_mut() {
                    *v += 2.0 * PI;
                }
            }
        }
        let mut freq = vec![0.0_f64; n];
        for i in 1..(n - 1) {
            freq[i] = (unwrapped[i + 1] - unwrapped[i - 1]) * fs / (4.0 * PI);
        }
        freq[0] = freq[1];
        freq[n - 1] = freq[n - 2];
        freq
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EMD (Empirical Mode Decomposition)
// ─────────────────────────────────────────────────────────────────────────────

/// Empirical Mode Decomposition (EMD) into Intrinsic Mode Functions (IMFs).
///
/// The sifting algorithm extracts IMFs by subtracting the mean envelope
/// (computed from local extrema interpolation) until a stopping criterion
/// is met.  The Hilbert-Huang spectrum can then be computed from the IMFs.
#[derive(Debug, Clone)]
pub struct EMD {
    /// Extracted Intrinsic Mode Functions (IMFs).
    pub imfs: Vec<Vec<f64>>,
    /// Residual signal after all IMFs are extracted.
    pub residual: Vec<f64>,
}

impl EMD {
    /// Decompose `signal` into IMFs via the sifting algorithm.
    ///
    /// `max_imfs` limits the number of IMFs extracted.
    /// `max_sift_iter` limits the sifting iterations per IMF.
    /// `sift_threshold` is the stopping threshold for sifting.
    pub fn decompose(
        signal: &[f64],
        max_imfs: usize,
        max_sift_iter: usize,
        sift_threshold: f64,
    ) -> Self {
        let n = signal.len();
        let mut residual = signal.to_vec();
        let mut imfs = Vec::new();

        for _ in 0..max_imfs {
            if residual.iter().all(|&x| x.abs() < sift_threshold) {
                break;
            }
            let (maxima_idx, _) = find_local_maxima(&residual);
            let (minima_idx, _) = find_local_minima(&residual);
            if maxima_idx.len() < 2 || minima_idx.len() < 2 {
                break;
            }
            let mut h = residual.clone();
            for _ in 0..max_sift_iter {
                let upper = cubic_spline_envelope(&h, &find_local_maxima(&h).0, true);
                let lower = cubic_spline_envelope(&h, &find_local_minima(&h).0, false);
                let mean_env: Vec<f64> = upper
                    .iter()
                    .zip(lower.iter())
                    .map(|(u, l)| (u + l) / 2.0)
                    .collect();
                let prev = h.clone();
                for i in 0..n {
                    h[i] -= mean_env[i];
                }
                // SD stopping criterion
                let sd: f64 = prev
                    .iter()
                    .zip(h.iter())
                    .map(|(p, c)| (p - c) * (p - c) / (p * p + 1e-15))
                    .sum::<f64>()
                    / n as f64;
                if sd < sift_threshold {
                    break;
                }
            }
            for i in 0..n {
                residual[i] -= h[i];
            }
            imfs.push(h);
        }
        Self { imfs, residual }
    }

    /// Compute the Hilbert-Huang spectrum (instantaneous frequency vs time for each IMF).
    ///
    /// `fs` is the sampling frequency.
    /// Returns one frequency vector per IMF.
    pub fn hilbert_huang_spectrum(&self, fs: f64) -> Vec<Vec<f64>> {
        self.imfs
            .iter()
            .map(|imf| {
                let ht = HilbertTransform::compute(imf);
                ht.instantaneous_frequency(fs)
            })
            .collect()
    }
}

/// Find local maxima of a signal.  Returns (indices, values).
fn find_local_maxima(x: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let n = x.len();
    let mut idx = Vec::new();
    let mut vals = Vec::new();
    for i in 1..(n.saturating_sub(1)) {
        if x[i] > x[i - 1] && x[i] > x[i + 1] {
            idx.push(i);
            vals.push(x[i]);
        }
    }
    (idx, vals)
}

/// Find local minima of a signal.  Returns (indices, values).
fn find_local_minima(x: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let n = x.len();
    let mut idx = Vec::new();
    let mut vals = Vec::new();
    for i in 1..(n.saturating_sub(1)) {
        if x[i] < x[i - 1] && x[i] < x[i + 1] {
            idx.push(i);
            vals.push(x[i]);
        }
    }
    (idx, vals)
}

/// Build a piecewise-linear (degenerate cubic spline) envelope from extrema.
fn cubic_spline_envelope(signal: &[f64], extrema_idx: &[usize], upper: bool) -> Vec<f64> {
    let n = signal.len();
    if extrema_idx.is_empty() {
        let fill = if upper {
            signal.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        } else {
            signal.iter().cloned().fold(f64::INFINITY, f64::min)
        };
        return vec![fill; n];
    }
    // Add boundary extrema
    let mut xs: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    xs.push(0.0);
    ys.push(signal[extrema_idx[0]]);
    for &i in extrema_idx {
        xs.push(i as f64);
        ys.push(signal[i]);
    }
    xs.push((n - 1) as f64);
    ys.push(signal[extrema_idx[extrema_idx.len() - 1]]);

    // Piecewise linear interpolation
    let mut env = vec![0.0_f64; n];
    for i in 0..n {
        let xi = i as f64;
        // Find surrounding xs
        let mut lo = 0usize;
        let mut hi = xs.len() - 1;
        for k in 0..(xs.len() - 1) {
            if xs[k] <= xi && xi <= xs[k + 1] {
                lo = k;
                hi = k + 1;
                break;
            }
        }
        let t = if (xs[hi] - xs[lo]).abs() < 1e-15 {
            0.0
        } else {
            (xi - xs[lo]) / (xs[hi] - xs[lo])
        };
        env[i] = ys[lo] + t * (ys[hi] - ys[lo]);
    }
    env
}

// ─────────────────────────────────────────────────────────────────────────────
// PeakDetection
// ─────────────────────────────────────────────────────────────────────────────

/// Detected peak in a signal.
#[derive(Debug, Clone, PartialEq)]
pub struct Peak {
    /// Sample index of the peak.
    pub index: usize,
    /// Signal value at the peak.
    pub value: f64,
    /// Prominence of the peak (height above surrounding baseline).
    pub prominence: f64,
    /// Estimated half-peak width (in samples).
    pub width: f64,
}

/// Peak detection algorithm.
#[derive(Debug, Clone)]
pub struct PeakDetection {
    /// Detected peaks.
    pub peaks: Vec<Peak>,
}

impl PeakDetection {
    /// Detect local maxima in `signal` above `min_height`.
    ///
    /// `min_distance` enforces a minimum sample separation between peaks.
    /// `min_prominence` filters out peaks below this prominence.
    pub fn find_peaks(
        signal: &[f64],
        min_height: f64,
        min_distance: usize,
        min_prominence: f64,
    ) -> Self {
        let n = signal.len();
        // Collect raw local maxima
        let mut candidates: Vec<usize> = (1..(n.saturating_sub(1)))
            .filter(|&i| {
                signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] >= min_height
            })
            .collect();

        // Apply min_distance: keep only the tallest in each window
        if min_distance > 0 {
            candidates = apply_min_distance(signal, &candidates, min_distance);
        }

        // Compute prominence and width
        let peaks: Vec<Peak> = candidates
            .iter()
            .filter_map(|&i| {
                let prom = prominence(signal, i);
                if prom < min_prominence {
                    return None;
                }
                let w = half_peak_width(signal, i, signal[i] - prom / 2.0);
                Some(Peak {
                    index: i,
                    value: signal[i],
                    prominence: prom,
                    width: w,
                })
            })
            .collect();

        Self { peaks }
    }

    /// Return peak indices sorted by decreasing prominence.
    pub fn sorted_by_prominence(&self) -> Vec<usize> {
        let mut sorted = self.peaks.clone();
        sorted.sort_by(|a, b| {
            b.prominence
                .partial_cmp(&a.prominence)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        sorted.iter().map(|p| p.index).collect()
    }
}

fn apply_min_distance(signal: &[f64], candidates: &[usize], min_dist: usize) -> Vec<usize> {
    let mut result: Vec<usize> = Vec::new();
    for &c in candidates {
        if result
            .iter()
            .all(|&r| (c as i64 - r as i64).unsigned_abs() as usize >= min_dist)
        {
            result.push(c);
        } else {
            // Replace if higher
            if let Some(close) = result
                .iter_mut()
                .find(|r| ((c as i64 - **r as i64).unsigned_abs() as usize) < min_dist)
                && signal[c] > signal[*close]
            {
                *close = c;
            }
        }
    }
    result
}

fn prominence(signal: &[f64], peak_idx: usize) -> f64 {
    let n = signal.len();
    let val = signal[peak_idx];
    // Find the minimum between peak and nearest higher peak on each side
    let left_min = (0..peak_idx)
        .rev()
        .take(100)
        .map(|i| signal[i])
        .fold(f64::INFINITY, f64::min);
    let right_min = ((peak_idx + 1)..n.min(peak_idx + 101))
        .map(|i| signal[i])
        .fold(f64::INFINITY, f64::min);
    let base = left_min.min(right_min);
    if base.is_infinite() { val } else { val - base }
}

fn half_peak_width(signal: &[f64], peak_idx: usize, level: f64) -> f64 {
    let n = signal.len();
    let mut left = peak_idx;
    let mut right = peak_idx;
    while left > 0 && signal[left] > level {
        left -= 1;
    }
    while right < n - 1 && signal[right] > level {
        right += 1;
    }
    (right - left) as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// SignalFilter
// ─────────────────────────────────────────────────────────────────────────────

/// Collection of signal filtering methods.
#[derive(Debug, Clone)]
pub struct SignalFilter;

impl SignalFilter {
    /// Apply a simple moving average filter.
    ///
    /// The output at index `i` is the mean of `signal[i-window/2 .. i+window/2]`
    /// (with boundary clamping).
    pub fn moving_average(signal: &[f64], window: usize) -> Vec<f64> {
        let n = signal.len();
        let hw = window / 2;
        (0..n)
            .map(|i| {
                let lo = i.saturating_sub(hw);
                let hi = (i + hw + 1).min(n);
                let count = hi - lo;
                signal[lo..hi].iter().sum::<f64>() / count as f64
            })
            .collect()
    }

    /// Apply exponential smoothing (single exponential).
    ///
    /// `alpha` in (0, 1]: smoothing factor.  Larger alpha = less smoothing.
    pub fn exponential_smoothing(signal: &[f64], alpha: f64) -> Vec<f64> {
        let n = signal.len();
        if n == 0 {
            return vec![];
        }
        let alpha = alpha.clamp(1e-6, 1.0);
        let mut out = vec![0.0_f64; n];
        out[0] = signal[0];
        for i in 1..n {
            out[i] = alpha * signal[i] + (1.0 - alpha) * out[i - 1];
        }
        out
    }

    /// Apply a Savitzky-Golay filter (polynomial smoothing).
    ///
    /// `window` must be odd and >= `poly_order + 1`.
    /// Uses a polynomial of degree `poly_order` fitted to each local window.
    ///
    /// This implementation uses a direct least-squares fit per window.
    pub fn savitzky_golay(signal: &[f64], window: usize, poly_order: usize) -> Vec<f64> {
        let n = signal.len();
        if n == 0 || window > n {
            return signal.to_vec();
        }
        let window = if window.is_multiple_of(2) {
            window + 1
        } else {
            window
        };
        let hw = window / 2;
        let poly_order = poly_order.min(window - 1);

        let mut out = vec![0.0_f64; n];
        for i in 0..n {
            let lo = i.saturating_sub(hw);
            let hi = (i + hw + 1).min(n);
            let seg = &signal[lo..hi];
            let m = seg.len();
            let center = (i - lo) as f64;
            // Build Vandermonde matrix and solve for central coefficient
            // Using a simple analytical approximation: weighted polynomial fit
            let xs: Vec<f64> = (0..m).map(|k| k as f64 - center).collect();
            out[i] = poly_fit_center(seg, &xs, poly_order);
        }
        out
    }

    /// Apply a simple FIR low-pass filter using a Hann-windowed sinc kernel.
    ///
    /// `cutoff` is the normalized cutoff frequency (0..0.5).
    /// `num_taps` is the number of filter taps (odd).
    pub fn lowpass_fir(signal: &[f64], cutoff: f64, num_taps: usize) -> Vec<f64> {
        let num_taps = if num_taps.is_multiple_of(2) {
            num_taps + 1
        } else {
            num_taps
        };
        let hw = num_taps / 2;
        // Build sinc kernel with Hann window
        let kernel: Vec<f64> = (0..num_taps)
            .map(|i| {
                let k = i as f64 - hw as f64;
                let sinc = if k.abs() < 1e-10 {
                    2.0 * cutoff
                } else {
                    (2.0 * PI * cutoff * k).sin() / (PI * k)
                };
                let hann = 0.5 * (1.0 - (2.0 * PI * i as f64 / (num_taps - 1) as f64).cos());
                sinc * hann
            })
            .collect();
        // Normalize
        let sum: f64 = kernel.iter().sum();
        let kernel: Vec<f64> = kernel.iter().map(|&k| k / sum.max(1e-15)).collect();
        // Convolve (valid mode padded with edge values)
        let n = signal.len();
        (0..n)
            .map(|i| {
                kernel
                    .iter()
                    .enumerate()
                    .map(|(k, &w)| {
                        let j = i as i64 + k as i64 - hw as i64;
                        let j = j.clamp(0, n as i64 - 1) as usize;
                        w * signal[j]
                    })
                    .sum()
            })
            .collect()
    }
}

/// Fit a polynomial of degree `p` to `(xs, ys)` and return the value at x=0.
fn poly_fit_center(ys: &[f64], xs: &[f64], p: usize) -> f64 {
    let m = ys.len();
    let deg = p + 1;
    // Build Vandermonde matrix V (m x deg)
    let mut vt: Vec<Vec<f64>> = vec![vec![0.0; m]; deg];
    for (j, &x) in xs.iter().enumerate() {
        let mut xp = 1.0_f64;
        for i in 0..deg {
            vt[i][j] = xp;
            xp *= x;
        }
    }
    // Compute V^T V and V^T y
    let mut vtv = vec![vec![0.0_f64; deg]; deg];
    let mut vty = vec![0.0_f64; deg];
    for i in 0..deg {
        for k in 0..m {
            vty[i] += vt[i][k] * ys[k];
            for j in 0..deg {
                vtv[i][j] += vt[i][k] * vt[j][k];
            }
        }
    }
    // Solve via Gaussian elimination
    let coeffs = gaussian_elimination(&mut vtv, &mut vty);
    // Evaluate at x=0: only constant term
    coeffs.first().copied().unwrap_or(ys[m / 2])
}

/// Gaussian elimination for Ax = b (modifies A and b in place).
fn gaussian_elimination(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot
        let pivot = (col..n)
            .max_by(|&i, &j| {
                a[i][col]
                    .abs()
                    .partial_cmp(&a[j][col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(col);
        a.swap(col, pivot);
        b.swap(col, pivot);
        if a[col][col].abs() < 1e-15 {
            continue;
        }
        let diag = a[col][col];
        for row in (col + 1)..n {
            let factor = a[row][col] / diag;
            for k in col..n {
                let v = a[col][k];
                a[row][k] -= factor * v;
            }
            let bv = b[col];
            b[row] -= factor * bv;
        }
    }
    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = if a[i][i].abs() > 1e-15 {
            (b[i] - s) / a[i][i]
        } else {
            0.0
        };
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn sine_wave(n: usize, freq: f64, fs: f64) -> Vec<f64> {
        (0..n)
            .map(|i| (2.0 * PI * freq * i as f64 / fs).sin())
            .collect()
    }

    fn ramp(n: usize) -> Vec<f64> {
        (0..n).map(|i| i as f64).collect()
    }

    // ── AutoCorrelation ───────────────────────────────────────────────────────

    #[test]
    fn test_autocorr_lag0_equals_variance() {
        let sig = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let ac = AutoCorrelation::compute(&sig, 2, false, false);
        let mean = sig.iter().sum::<f64>() / sig.len() as f64;
        let var: f64 = sig.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / sig.len() as f64;
        assert!((ac.at_lag(0) - var).abs() < 1e-10);
    }

    #[test]
    fn test_autocorr_normalized_lag0_is_one() {
        let sig: Vec<f64> = (0..20).map(|i| (i as f64).sin()).collect();
        let ac = AutoCorrelation::compute(&sig, 5, false, true);
        assert!((ac.at_lag(0) - 1.0_f64).abs() < 1e-10);
    }

    #[test]
    fn test_autocorr_periodic_signal() {
        let n = 128;
        let sig = sine_wave(n, 8.0, n as f64);
        let ac = AutoCorrelation::compute(&sig, n / 2, false, true);
        // For a pure sine, autocorr at lag n/freq should be close to 1
        let period = (n as f64 / 8.0) as usize;
        assert!(ac.at_lag(period).abs() > 0.5);
    }

    #[test]
    fn test_autocorr_unbiased_vs_biased() {
        let sig = vec![1.0_f64, 2.0, 3.0, 4.0];
        let biased = AutoCorrelation::compute(&sig, 1, false, false);
        let unbiased = AutoCorrelation::compute(&sig, 1, true, false);
        // Biased uses N, unbiased uses N-k; at lag 1 unbiased should be larger
        assert!(unbiased.at_lag(1).abs() >= biased.at_lag(1).abs());
    }

    #[test]
    fn test_autocorr_constant_signal_zero_lag_positive() {
        let sig = vec![3.0_f64; 10];
        let ac = AutoCorrelation::compute(&sig, 3, false, false);
        assert!(ac.at_lag(0) >= 0.0);
    }

    #[test]
    fn test_autocorr_dominant_lag_returns_index() {
        let sig: Vec<f64> = (0..20).map(|i| (i as f64).sin()).collect();
        let ac = AutoCorrelation::compute(&sig, 10, false, false);
        let dl = ac.dominant_lag();
        assert!((1..=10).contains(&dl));
    }

    // ── CrossCorrelation ──────────────────────────────────────────────────────

    #[test]
    fn test_crosscorr_identical_signals() {
        let sig = vec![1.0_f64, 2.0, 3.0, 2.0, 1.0];
        let cc = CrossCorrelation::compute(&sig, &sig, 2);
        // Peak should be at lag 0
        assert_eq!(cc.peak_lag(), 0);
    }

    #[test]
    fn test_crosscorr_shifted_signal() {
        let n = 32;
        // Use delay=2 with period=8 and max_lag=5 so the aliased peak at
        // -(period - delay) = -6 falls outside the search range, giving a
        // unique positive peak at lag = -delay = -2.
        let delay = 2usize;
        let x: Vec<f64> = (0..n).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let y: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * (i + delay) as f64 / 8.0).sin())
            .collect();
        let cc = CrossCorrelation::compute(&x, &y, 5);
        let lag = cc.peak_lag();
        // R_xy(k) = -½ sin(2πk/8); maximum positive value at k = -delay = -2
        assert!(
            (lag + delay as i64).abs() <= 1,
            "lag={lag}, expected near -{delay}"
        );
    }

    #[test]
    fn test_crosscorr_at_lag() {
        let x = vec![1.0_f64, 0.0, -1.0, 0.0, 1.0];
        let y = vec![0.0_f64, 1.0, 0.0, -1.0, 0.0];
        let cc = CrossCorrelation::compute(&x, &y, 2);
        // Value at lag 0 should be finite
        assert!(cc.at_lag(0).is_finite());
    }

    #[test]
    fn test_crosscorr_out_of_range_lag() {
        let x = vec![1.0_f64, 2.0];
        let y = vec![1.0_f64, 2.0];
        let cc = CrossCorrelation::compute(&x, &y, 1);
        assert_eq!(cc.at_lag(100), 0.0);
    }

    #[test]
    fn test_crosscorr_values_length() {
        let x = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let y = vec![5.0_f64, 4.0, 3.0, 2.0, 1.0];
        let cc = CrossCorrelation::compute(&x, &y, 3);
        assert_eq!(cc.values.len(), 7); // 2*3+1
    }

    // ── PowerSpectralDensity ──────────────────────────────────────────────────

    #[test]
    fn test_psd_frequencies_start_at_zero() {
        let sig = vec![0.0_f64; 64];
        let psd = PowerSpectralDensity::welch(&sig, 100.0, 32, 0.5);
        assert!((psd.frequencies[0]).abs() < 1e-12);
    }

    #[test]
    fn test_psd_power_nonneg() {
        let sig: Vec<f64> = (0..128)
            .map(|i| (2.0 * PI * 10.0 * i as f64 / 128.0).sin())
            .collect();
        let psd = PowerSpectralDensity::welch(&sig, 128.0, 32, 0.5);
        assert!(psd.power.iter().all(|&p| p >= 0.0));
    }

    #[test]
    fn test_psd_dominant_freq() {
        let fs = 256.0_f64;
        let f0 = 16.0_f64;
        let sig: Vec<f64> = (0..256)
            .map(|i| (2.0 * PI * f0 * i as f64 / fs).sin())
            .collect();
        let psd = PowerSpectralDensity::welch(&sig, fs, 64, 0.5);
        let dominant = psd.dominant_frequency();
        assert!(
            (dominant - f0).abs() < 5.0,
            "expected ~{f0} Hz, got {dominant}"
        );
    }

    #[test]
    fn test_psd_freq_resolution() {
        let sig = vec![0.0_f64; 64];
        let psd = PowerSpectralDensity::welch(&sig, 100.0, 32, 0.0);
        let df = psd.freq_resolution();
        assert!(df > 0.0);
    }

    #[test]
    fn test_psd_total_power_nonneg() {
        let sig: Vec<f64> = (0..64).map(|i| i as f64).collect();
        let psd = PowerSpectralDensity::welch(&sig, 64.0, 16, 0.0);
        assert!(psd.total_power() >= 0.0);
    }

    // ── HilbertTransform ──────────────────────────────────────────────────────

    #[test]
    fn test_hilbert_real_part_matches_input() {
        let sig: Vec<f64> = (0..32).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let ht = HilbertTransform::compute(&sig);
        // Real part may deviate slightly due to zero-padding
        for (r, s) in ht.real_part.iter().zip(sig.iter()) {
            assert!((r - s).abs() < 0.2, "real part deviation: {} vs {}", r, s);
        }
    }

    #[test]
    fn test_hilbert_amplitude_nonneg() {
        let sig: Vec<f64> = (0..64).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let ht = HilbertTransform::compute(&sig);
        assert!(ht.amplitude.iter().all(|&a| a >= 0.0));
    }

    #[test]
    fn test_hilbert_phase_range() {
        let sig: Vec<f64> = (0..64).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let ht = HilbertTransform::compute(&sig);
        assert!(
            ht.phase
                .iter()
                .all(|&p| (-PI - 1e-9..=PI + 1e-9).contains(&p))
        );
    }

    #[test]
    fn test_hilbert_instantaneous_freq() {
        let fs = 64.0_f64;
        let f0 = 4.0_f64;
        let sig: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * f0 * i as f64 / fs).sin())
            .collect();
        let ht = HilbertTransform::compute(&sig);
        let ifreq = ht.instantaneous_frequency(fs);
        // Central samples should be close to f0
        let mid_freq: f64 = ifreq[16..48].iter().sum::<f64>() / 32.0;
        assert!(
            (mid_freq - f0).abs() < 2.0,
            "instantaneous freq {mid_freq} vs expected {f0}"
        );
    }

    #[test]
    fn test_hilbert_analytic_length() {
        let n = 50;
        let sig = vec![1.0_f64; n];
        let ht = HilbertTransform::compute(&sig);
        assert_eq!(ht.real_part.len(), n);
        assert_eq!(ht.imaginary_part.len(), n);
        assert_eq!(ht.amplitude.len(), n);
        assert_eq!(ht.phase.len(), n);
    }

    // ── EMD ───────────────────────────────────────────────────────────────────

    #[test]
    fn test_emd_sum_of_imfs_plus_residual() {
        let sig: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * i as f64 / 8.0).sin() + 0.5 * (2.0 * PI * i as f64 / 16.0).sin())
            .collect();
        let emd = EMD::decompose(&sig, 3, 10, 0.05);
        let mut recon = emd.residual.clone();
        for imf in &emd.imfs {
            for (i, &v) in imf.iter().enumerate() {
                recon[i] += v;
            }
        }
        let err: f64 = sig
            .iter()
            .zip(recon.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>();
        assert!(err.sqrt() < 1.0, "reconstruction error: {err}");
    }

    #[test]
    fn test_emd_imfs_not_empty_for_rich_signal() {
        let sig: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * i as f64 / 8.0).sin() + 0.3 * (i as f64).cos())
            .collect();
        let emd = EMD::decompose(&sig, 4, 10, 0.01);
        assert!(!emd.imfs.is_empty());
    }

    #[test]
    fn test_emd_hilbert_huang_spectrum_shape() {
        let sig: Vec<f64> = (0..32).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let emd = EMD::decompose(&sig, 2, 5, 0.1);
        let hhs = emd.hilbert_huang_spectrum(32.0);
        for row in &hhs {
            assert_eq!(row.len(), 32);
        }
    }

    // ── PeakDetection ─────────────────────────────────────────────────────────

    #[test]
    fn test_peaks_simple_sine() {
        let n = 64;
        let sig: Vec<f64> = (0..n).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let pd = PeakDetection::find_peaks(&sig, 0.5, 4, 0.0);
        assert!(!pd.peaks.is_empty());
    }

    #[test]
    fn test_peaks_heights_above_threshold() {
        let sig = vec![0.0_f64, 1.0, 0.5, 2.0, 0.3, 1.5, 0.0];
        let pd = PeakDetection::find_peaks(&sig, 0.8, 1, 0.0);
        for p in &pd.peaks {
            assert!(p.value >= 0.8_f64);
        }
    }

    #[test]
    fn test_peaks_no_peaks_in_flat_signal() {
        let sig = vec![1.0_f64; 20];
        let pd = PeakDetection::find_peaks(&sig, 0.5, 1, 0.0);
        assert!(pd.peaks.is_empty());
    }

    #[test]
    fn test_peaks_prominence_positive() {
        let sig = vec![0.0_f64, 1.0, 0.0, 2.0, 0.0];
        let pd = PeakDetection::find_peaks(&sig, 0.0, 1, 0.0);
        for p in &pd.peaks {
            assert!(p.prominence >= 0.0);
        }
    }

    #[test]
    fn test_peaks_sorted_by_prominence() {
        let sig = vec![0.0_f64, 1.0, 0.0, 3.0, 0.0, 2.0, 0.0];
        let pd = PeakDetection::find_peaks(&sig, 0.0, 1, 0.0);
        let sorted = pd.sorted_by_prominence();
        // Index of highest peak (3.0 at index 3) should come first
        if !sorted.is_empty() {
            assert_eq!(sorted[0], 3);
        }
    }

    #[test]
    fn test_peaks_width_nonneg() {
        let sig: Vec<f64> = (0..32).map(|i| (2.0 * PI * i as f64 / 8.0).sin()).collect();
        let pd = PeakDetection::find_peaks(&sig, 0.0, 2, 0.0);
        for p in &pd.peaks {
            assert!(p.width >= 0.0);
        }
    }

    // ── SignalFilter ──────────────────────────────────────────────────────────

    #[test]
    fn test_moving_average_constant() {
        let sig = vec![5.0_f64; 20];
        let filtered = SignalFilter::moving_average(&sig, 5);
        for &v in &filtered {
            assert!((v - 5.0_f64).abs() < 1e-12);
        }
    }

    #[test]
    fn test_moving_average_length_preserved() {
        let sig = ramp(30);
        let out = SignalFilter::moving_average(&sig, 5);
        assert_eq!(out.len(), 30);
    }

    #[test]
    fn test_exp_smoothing_constant_signal() {
        let sig = vec![3.0_f64; 10];
        let out = SignalFilter::exponential_smoothing(&sig, 0.3);
        // Should converge to 3
        assert!((out.last().unwrap() - 3.0_f64).abs() < 0.1);
    }

    #[test]
    fn test_exp_smoothing_length_preserved() {
        let sig = ramp(15);
        let out = SignalFilter::exponential_smoothing(&sig, 0.5);
        assert_eq!(out.len(), 15);
    }

    #[test]
    fn test_exp_smoothing_smooths_noise() {
        let sig: Vec<f64> = (0..20)
            .map(|i| if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 })
            .collect();
        let out = SignalFilter::exponential_smoothing(&sig, 0.1);
        // Output should have smaller range than input
        let range_out = out.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - out.iter().cloned().fold(f64::INFINITY, f64::min);
        let range_in = 2.0_f64;
        assert!(range_out < range_in);
    }

    #[test]
    fn test_savitzky_golay_constant() {
        let sig = vec![7.0_f64; 20];
        let out = SignalFilter::savitzky_golay(&sig, 5, 2);
        for &v in &out {
            assert!((v - 7.0_f64).abs() < 1e-8);
        }
    }

    #[test]
    fn test_savitzky_golay_length_preserved() {
        let sig = ramp(25);
        let out = SignalFilter::savitzky_golay(&sig, 5, 2);
        assert_eq!(out.len(), 25);
    }

    #[test]
    fn test_lowpass_fir_length_preserved() {
        let sig: Vec<f64> = (0..64).map(|i| i as f64).collect();
        let out = SignalFilter::lowpass_fir(&sig, 0.1, 15);
        assert_eq!(out.len(), 64);
    }

    #[test]
    fn test_lowpass_fir_attenuates_high_freq() {
        let n = 256;
        let fs = 256.0_f64;
        // Mix of low (4 Hz) and high (80 Hz) frequency
        let sig: Vec<f64> = (0..n)
            .map(|i| {
                (2.0 * PI * 4.0 * i as f64 / fs).sin() + (2.0 * PI * 80.0 * i as f64 / fs).sin()
            })
            .collect();
        let cutoff = 0.1_f64; // normalized ~25.6 Hz
        let out = SignalFilter::lowpass_fir(&sig, cutoff, 31);
        // Power of output should be less than input
        let p_in: f64 = sig.iter().map(|x| x * x).sum::<f64>();
        let p_out: f64 = out.iter().map(|x| x * x).sum::<f64>();
        assert!(p_out < p_in);
    }
}
