//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use crate::complex::Complex;
use std::f64::consts::PI;

/// Supported window functions for spectral analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WindowType {
    /// Hann (von Hann) window — good general-purpose window.
    Hann,
    /// Hamming window — slightly higher sidelobes than Hann.
    Hamming,
    /// Blackman window — very low sidelobes.
    Blackman,
    /// Flat-top window — accurate amplitude measurements.
    FlatTop,
    /// Rectangular (boxcar) window — no windowing.
    Rectangular,
}
impl WindowType {
    /// Generate the window coefficients of length `n`.
    pub fn generate(&self, n: usize) -> Vec<f64> {
        match self {
            WindowType::Hann => hann_window(n),
            WindowType::Hamming => hamming_window(n),
            WindowType::Blackman => blackman_window(n),
            WindowType::FlatTop => flat_top_window(n),
            WindowType::Rectangular => vec![1.0_f64; n],
        }
    }
}
/// Arbitrary FIR (Finite Impulse Response) filter.
///
/// Supports both batch processing via [`FirFilter::apply`] and online
/// sample-by-sample processing via [`FirFilter::step`].
///
/// # Example
/// ```no_run
/// # use oxiphysics_core::signal::FirFilter;
/// let coeffs = vec![0.25, 0.5, 0.25]; // simple low-pass
/// let mut fir = FirFilter::new(coeffs);
/// let out = fir.apply(&[1.0, 0.0, 0.0, 0.0, 0.0]);
/// ```
pub struct FirFilter {
    /// Filter coefficients (tap weights), index 0 = most recent.
    pub coeffs: Vec<f64>,
    /// Internal delay line (ring buffer).
    pub(super) delay: Vec<f64>,
    /// Current write position in the ring buffer.
    pub(super) pos: usize,
}
impl FirFilter {
    /// Create a new `FirFilter` from the given tap coefficients.
    pub fn new(coeffs: Vec<f64>) -> Self {
        let n = coeffs.len();
        Self {
            coeffs,
            delay: vec![0.0_f64; n],
            pos: 0,
        }
    }
    /// Reset the internal delay line to zero.
    pub fn reset(&mut self) {
        for v in self.delay.iter_mut() {
            *v = 0.0;
        }
        self.pos = 0;
    }
    /// Process a single sample and return the filtered output.
    pub fn step(&mut self, sample: f64) -> f64 {
        let n = self.coeffs.len();
        if n == 0 {
            return sample;
        }
        self.delay[self.pos] = sample;
        let mut acc = 0.0_f64;
        for (k, &c) in self.coeffs.iter().enumerate() {
            let idx = (self.pos + n - k) % n;
            acc += c * self.delay[idx];
        }
        self.pos = (self.pos + 1) % n;
        acc
    }
    /// Filter an entire signal block and return the output vector.
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&x| self.step(x)).collect()
    }
}
/// Short-Time Fourier Transform (STFT).
///
/// Splits `signal` into overlapping frames, applies a window function, computes
/// the FFT of each frame, and returns the 2-D spectrogram.
pub struct Stft {
    /// Frame (window) length in samples.
    pub frame_len: usize,
    /// Hop size in samples (frame_len - overlap).
    pub hop_size: usize,
}
impl Stft {
    /// Create a new `Stft` with given frame length and hop size.
    pub fn new(frame_len: usize, hop_size: usize) -> Self {
        Self {
            frame_len,
            hop_size,
        }
    }
    /// Compute the STFT of `signal`.
    ///
    /// Returns a `Vec` of complex spectra, one per frame. Each spectrum has
    /// length `next_pow2(frame_len)`.
    pub fn compute(&self, signal: &[f64]) -> Vec<Vec<Complex>> {
        if signal.is_empty() || self.frame_len == 0 || self.hop_size == 0 {
            return Vec::new();
        }
        let window = hann_window(self.frame_len);
        let n = signal.len();
        let mut frames = Vec::new();
        let mut start = 0;
        while start + self.frame_len <= n {
            let frame: Vec<Complex> = (0..self.frame_len)
                .map(|i| Complex::new(signal[start + i] * window[i], 0.0))
                .collect();
            frames.push(fft(&frame));
            start += self.hop_size;
        }
        frames
    }
    /// Compute the STFT magnitude spectrogram (|X\[k\]| per frame).
    pub fn magnitude_spectrogram(&self, signal: &[f64]) -> Vec<Vec<f64>> {
        self.compute(signal)
            .into_iter()
            .map(|frame| {
                frame
                    .iter()
                    .map(|c| (c.re * c.re + c.im * c.im).sqrt())
                    .collect()
            })
            .collect()
    }
}
/// Welch's method for estimating the power spectral density (PSD).
///
/// Divides the signal into overlapping segments, applies a window function,
/// computes the FFT of each segment, and averages the squared magnitudes.
pub struct WelchPsd {
    /// Length of each segment (should be a power of 2).
    pub segment_len: usize,
    /// Hop (overlap step) between consecutive segments.
    pub hop: usize,
    /// Window function to apply to each segment.
    pub window_type: WindowType,
}
impl WelchPsd {
    /// Create a new Welch PSD estimator.
    pub fn new(segment_len: usize, hop: usize, window_type: WindowType) -> Self {
        Self {
            segment_len,
            hop,
            window_type,
        }
    }
    /// Compute the one-sided PSD estimate of `signal`.
    ///
    /// Returns a vector of length `segment_len / 2 + 1` containing the
    /// power at each frequency bin from DC to Nyquist.
    pub fn compute(&self, signal: &[f64]) -> Vec<f64> {
        let n = self.segment_len;
        let m = n / 2 + 1;
        let window = self.window_type.generate(n);
        let win_power: f64 = window.iter().map(|w| w * w).sum::<f64>();
        let normalisation = if win_power > 0.0 { win_power } else { 1.0 };
        let mut accum = vec![0.0_f64; m];
        let mut n_frames = 0usize;
        let mut start = 0usize;
        while start + n <= signal.len() {
            let frame: Vec<f64> = (0..n).map(|i| signal[start + i] * window[i]).collect();
            let spectrum = fft_real(&frame);
            for k in 0..m {
                let mag2 = spectrum[k].norm_sq();
                let scale = if k == 0 || k == m - 1 { 1.0 } else { 2.0 };
                accum[k] += scale * mag2 / normalisation;
            }
            n_frames += 1;
            start += self.hop;
        }
        if n_frames > 0 {
            for p in &mut accum {
                *p /= n_frames as f64;
            }
        }
        accum
    }
}
/// Spectrum peak detector: finds local maxima above a threshold.
pub struct SpectrumPeakDetector;
impl SpectrumPeakDetector {
    /// Find local maxima in `spectrum` that exceed `threshold`.
    ///
    /// # Arguments
    /// * `spectrum`      - Power spectrum or PSD vector.
    /// * `min_distance`  - Minimum separation (in bins) between peaks.
    /// * `threshold`     - Minimum power for a peak to be reported.
    ///
    /// Returns a sorted list of bin indices of detected peaks.
    pub fn find_peaks(spectrum: &[f64], min_distance: usize, threshold: f64) -> Vec<usize> {
        let n = spectrum.len();
        if n < 3 {
            return Vec::new();
        }
        let mut candidates: Vec<(usize, f64)> = Vec::new();
        for i in 1..n - 1 {
            if spectrum[i] > threshold
                && spectrum[i] >= spectrum[i - 1]
                && spectrum[i] >= spectrum[i + 1]
            {
                candidates.push((i, spectrum[i]));
            }
        }
        candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut peaks: Vec<usize> = Vec::new();
        for (idx, _power) in candidates {
            let too_close = peaks
                .iter()
                .any(|&p| (p as isize - idx as isize).unsigned_abs() < min_distance);
            if !too_close {
                peaks.push(idx);
            }
        }
        peaks.sort_unstable();
        peaks
    }
}
/// Biquad IIR filter coefficients (Direct Form II).
///
/// Transfer function: `H(z) = (b0 + b1*z⁻¹ + b2*z⁻²) / (1 + a1*z⁻¹ + a2*z⁻²)`
#[derive(Debug, Clone, Copy)]
pub struct BiquadCoeff {
    /// Feed-forward coefficient b0.
    pub b0: f64,
    /// Feed-forward coefficient b1.
    pub b1: f64,
    /// Feed-forward coefficient b2.
    pub b2: f64,
    /// Feedback coefficient a1.
    pub a1: f64,
    /// Feedback coefficient a2.
    pub a2: f64,
}
impl BiquadCoeff {
    /// Create a new set of biquad coefficients.
    pub fn new(b0: f64, b1: f64, b2: f64, a1: f64, a2: f64) -> Self {
        Self { b0, b1, b2, a1, a2 }
    }
}
/// A higher-level signal processor that bundles STFT-based spectral analysis,
/// Mel filterbank computation, and full autocorrelation functions.
pub struct SignalProcessor {
    /// Frame length (window size) for STFT-based methods.
    pub frame_len: usize,
    /// Hop size (step between frames) for STFT-based methods.
    pub hop_size: usize,
    /// Sample rate in Hz.
    pub sample_rate: f64,
}
impl SignalProcessor {
    /// Create a new `SignalProcessor`.
    pub fn new(frame_len: usize, hop_size: usize, sample_rate: f64) -> Self {
        Self {
            frame_len,
            hop_size,
            sample_rate,
        }
    }
    /// Compute a power spectrogram via STFT.
    ///
    /// Each row of the returned matrix corresponds to one time frame;
    /// each column is a frequency bin (0..frame_len/2+1).
    /// Values are power (`|X[k]|²` / `frame_len`).
    pub fn compute_spectrogram(&self, signal: &[f64]) -> Vec<Vec<f64>> {
        if signal.is_empty() || self.frame_len == 0 || self.hop_size == 0 {
            return Vec::new();
        }
        let stft = Stft::new(self.frame_len, self.hop_size);
        let complex_frames = stft.compute(signal);
        let n_bins = self.frame_len / 2 + 1;
        let norm = (self.frame_len as f64).powi(2);
        complex_frames
            .into_iter()
            .map(|frame| {
                (0..n_bins.min(frame.len()))
                    .map(|k| {
                        let re = frame[k].re;
                        let im = frame[k].im;
                        (re * re + im * im) / norm
                    })
                    .collect()
            })
            .collect()
    }
    /// Compute a Mel-scale filterbank matrix.
    ///
    /// Returns a `n_mels × n_fft_bins` matrix where each row is one triangular
    /// Mel filter.  Frequency bins span `[f_min, f_max]` (in Hz) converted
    /// through the standard Mel formula `m = 2595 * log10(1 + f/700)`.
    ///
    /// # Arguments
    /// * `n_mels`    - Number of Mel filters.
    /// * `n_fft`     - Number of FFT frequency bins (typically `frame_len / 2 + 1`).
    /// * `f_min`     - Lowest frequency (Hz).
    /// * `f_max`     - Highest frequency (Hz, ≤ Nyquist).
    pub fn compute_mel_filterbank(
        &self,
        n_mels: usize,
        n_fft: usize,
        f_min: f64,
        f_max: f64,
    ) -> Vec<Vec<f64>> {
        if n_mels == 0 || n_fft == 0 {
            return Vec::new();
        }
        let hz_to_mel = |f: f64| 2595.0 * (1.0 + f / 700.0).log10();
        let mel_to_hz = |m: f64| 700.0 * (10.0_f64.powf(m / 2595.0) - 1.0);
        let mel_min = hz_to_mel(f_min);
        let mel_max = hz_to_mel(f_max);
        let mel_points: Vec<f64> = (0..=n_mels + 1)
            .map(|i| mel_min + i as f64 * (mel_max - mel_min) / (n_mels + 1) as f64)
            .collect();
        let fft_bins: Vec<f64> = mel_points
            .iter()
            .map(|&m| {
                let hz = mel_to_hz(m);
                hz * (n_fft as f64 - 1.0) * 2.0 / self.sample_rate
            })
            .collect();
        (0..n_mels)
            .map(|m| {
                (0..n_fft)
                    .map(|k| {
                        let k_f = k as f64;
                        let lo = fft_bins[m];
                        let centre = fft_bins[m + 1];
                        let hi = fft_bins[m + 2];
                        if k_f < lo || k_f > hi {
                            0.0
                        } else if k_f <= centre {
                            if (centre - lo).abs() < f64::EPSILON {
                                1.0
                            } else {
                                (k_f - lo) / (centre - lo)
                            }
                        } else if (hi - centre).abs() < f64::EPSILON {
                            1.0
                        } else {
                            (hi - k_f) / (hi - centre)
                        }
                    })
                    .collect()
            })
            .collect()
    }
    /// Compute the full normalized autocorrelation function (ACF) up to `max_lag`.
    ///
    /// Returns a vector of length `max_lag + 1` where element `k` is the
    /// normalized autocorrelation at lag `k`.  At lag 0 the value is always 1.
    /// The normalisation uses `r(0)` (zero-lag covariance) so all values lie in
    /// `[-1, 1]`.
    pub fn compute_autocorrelation(&self, signal: &[f64], max_lag: usize) -> Vec<f64> {
        let n = signal.len();
        if n == 0 {
            return Vec::new();
        }
        let mean = signal.iter().sum::<f64>() / n as f64;
        let variance: f64 = signal.iter().map(|&x| (x - mean).powi(2)).sum();
        if variance < f64::EPSILON {
            return vec![1.0; (max_lag + 1).min(n)];
        }
        let lags = (max_lag + 1).min(n);
        (0..lags)
            .map(|lag| {
                let corr: f64 = signal[..n - lag]
                    .iter()
                    .zip(signal[lag..].iter())
                    .map(|(&a, &b)| (a - mean) * (b - mean))
                    .sum();
                corr / variance
            })
            .collect()
    }
}
/// 1-D (scalar) Kalman filter for tracking a slowly varying scalar state.
///
/// Model: `x[k] = x[k-1] + process_noise`, `z[k] = x[k] + measurement_noise`
pub struct KalmanFilter1D {
    /// Current state estimate.
    pub x: f64,
    /// Estimate error covariance.
    pub p: f64,
    /// Process noise variance.
    pub q: f64,
    /// Measurement noise variance.
    pub r: f64,
}
impl KalmanFilter1D {
    /// Create a new `KalmanFilter1D`.
    ///
    /// # Arguments
    /// * `x0` — Initial state estimate.
    /// * `p0` — Initial estimate error covariance.
    /// * `q` — Process noise variance.
    /// * `r` — Measurement noise variance.
    pub fn new(x0: f64, p0: f64, q: f64, r: f64) -> Self {
        Self { x: x0, p: p0, q, r }
    }
    /// Prediction step: project state and covariance forward.
    pub fn predict(&mut self) {
        self.p += self.q;
    }
    /// Update step: incorporate a new measurement `z`.
    pub fn update(&mut self, z: f64) {
        let k = self.p / (self.p + self.r);
        self.x += k * (z - self.x);
        self.p *= 1.0 - k;
    }
    /// Combined predict + update for convenience.
    pub fn step(&mut self, z: f64) -> f64 {
        self.predict();
        self.update(z);
        self.x
    }
}
/// Second-order Butterworth band-pass filter using cascaded biquad sections.
///
/// Designed via bilinear transform. Cascades `order` biquad sections, each
/// targeting the band `[f_low, f_high]`.
pub struct ButterworthBandpass {
    /// Biquad sections (cascaded).
    pub(super) sections: Vec<IirFilter>,
}
impl ButterworthBandpass {
    /// Design a Butterworth band-pass filter.
    ///
    /// # Arguments
    /// * `order`    - Number of biquad sections (total filter order = 2*order).
    /// * `f_low`    - Lower 3-dB cutoff frequency in Hz.
    /// * `f_high`   - Upper 3-dB cutoff frequency in Hz.
    /// * `sample_rate` - Sampling frequency in Hz.
    pub fn new(order: usize, f_low: f64, f_high: f64, sample_rate: f64) -> Self {
        let mut sections = Vec::new();
        let fs = sample_rate;
        let order = order.max(1);
        let w_low = 2.0 * fs * (PI * f_low / fs).tan();
        let w_high = 2.0 * fs * (PI * f_high / fs).tan();
        let bw = w_high - w_low;
        let w0 = (w_low * w_high).sqrt();
        for k in 0..order {
            let theta = PI * (2 * k + 1) as f64 / (2 * order) as f64;
            let alpha = bw / (2.0 * w0) * theta.sin();
            let w_c = 2.0 * PI * ((f_low + f_high) / 2.0) / fs;
            let q = (w0 / bw) * theta.sin();
            let r = (-PI * (f_high - f_low) / (q * fs)).exp();
            let cos_wc = w_c.cos();
            let sin_wc = w_c.sin();
            let bw_rad = 2.0 * PI * (f_high - f_low) / fs;
            let alpha_bw = sin_wc * (bw_rad / 2.0).sinh().min(0.5 * alpha).max(0.01);
            let b0 = alpha_bw;
            let b1 = 0.0;
            let b2 = -alpha_bw;
            let a0 = 1.0 + alpha_bw;
            let a1 = -2.0 * cos_wc;
            let a2 = 1.0 - alpha_bw;
            let _ = r;
            let coeff = BiquadCoeff::new(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0);
            sections.push(IirFilter::new(coeff));
        }
        Self { sections }
    }
    /// Apply the band-pass filter to a signal in-place (cascade of biquads).
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        let mut out: Vec<f64> = signal.to_vec();
        for sec in &mut self.sections {
            out = sec.apply(&out);
        }
        out
    }
    /// Reset all biquad states.
    pub fn reset(&mut self) {
        for sec in &mut self.sections {
            sec.reset();
        }
    }
}
/// General N-dimensional linear Kalman filter.
///
/// State model:    `x[k] = F * x[k-1] + process_noise`
/// Measurement:    `z[k] = H * x[k]   + measurement_noise`
///
/// All matrices are stored in row-major flat `Vec`f64` with dimension metadata.
pub struct KalmanFilterND {
    /// State dimension n.
    pub n: usize,
    /// Measurement dimension m.
    pub m: usize,
    /// State estimate vector (length n).
    pub x: Vec<f64>,
    /// Error covariance matrix P (n×n, row-major).
    pub p: Vec<f64>,
    /// State transition matrix F (n×n, row-major).
    pub f: Vec<f64>,
    /// Measurement matrix H (m×n, row-major).
    pub h: Vec<f64>,
    /// Process noise covariance Q (n×n, row-major).
    pub q: Vec<f64>,
    /// Measurement noise covariance R (m×m, row-major).
    pub r: Vec<f64>,
}
impl KalmanFilterND {
    /// Create a new `KalmanFilterND`.
    ///
    /// # Arguments
    /// * `n` — State dimension.
    /// * `m` — Measurement dimension.
    /// * `x0` — Initial state (length n).
    /// * `p0` — Initial error covariance (n×n, row-major).
    /// * `f` — State transition matrix (n×n, row-major).
    /// * `h` — Measurement matrix (m×n, row-major).
    /// * `q` — Process noise covariance (n×n, row-major).
    /// * `r` — Measurement noise covariance (m×m, row-major).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        n: usize,
        m: usize,
        x0: Vec<f64>,
        p0: Vec<f64>,
        f: Vec<f64>,
        h: Vec<f64>,
        q: Vec<f64>,
        r: Vec<f64>,
    ) -> Self {
        Self {
            n,
            m,
            x: x0,
            p: p0,
            f,
            h,
            q,
            r,
        }
    }
    fn mat_mul(a: &[f64], b: &[f64], rows_a: usize, cols_a: usize, cols_b: usize) -> Vec<f64> {
        let mut c = vec![0.0_f64; rows_a * cols_b];
        for i in 0..rows_a {
            for j in 0..cols_b {
                for k in 0..cols_a {
                    c[i * cols_b + j] += a[i * cols_a + k] * b[k * cols_b + j];
                }
            }
        }
        c
    }
    fn mat_add(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
    }
    fn mat_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
    }
    fn mat_transpose(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
        let mut t = vec![0.0_f64; rows * cols];
        for i in 0..rows {
            for j in 0..cols {
                t[j * rows + i] = a[i * cols + j];
            }
        }
        t
    }
    fn mat_vec_mul(a: &[f64], x: &[f64], rows: usize, cols: usize) -> Vec<f64> {
        (0..rows)
            .map(|i| (0..cols).map(|j| a[i * cols + j] * x[j]).sum())
            .collect()
    }
    fn mat_inv_2x2(a: &[f64]) -> Option<Vec<f64>> {
        let det = a[0] * a[3] - a[1] * a[2];
        if det.abs() < 1e-15 {
            return None;
        }
        Some(vec![a[3] / det, -a[1] / det, -a[2] / det, a[0] / det])
    }
    /// Invert a square matrix using Gaussian elimination.
    fn mat_inv(a: &[f64], n: usize) -> Option<Vec<f64>> {
        if n == 1 {
            if a[0].abs() < 1e-15 {
                return None;
            }
            return Some(vec![1.0 / a[0]]);
        }
        if n == 2 {
            return Self::mat_inv_2x2(a);
        }
        let mut aug = vec![0.0_f64; n * 2 * n];
        for i in 0..n {
            for j in 0..n {
                aug[i * 2 * n + j] = a[i * n + j];
            }
            aug[i * 2 * n + n + i] = 1.0;
        }
        for col in 0..n {
            let pivot_row = (col..n).max_by(|&r1, &r2| {
                aug[r1 * 2 * n + col]
                    .abs()
                    .partial_cmp(&aug[r2 * 2 * n + col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })?;
            if pivot_row != col {
                for j in 0..2 * n {
                    aug.swap(col * 2 * n + j, pivot_row * 2 * n + j);
                }
            }
            let pivot = aug[col * 2 * n + col];
            if pivot.abs() < 1e-15 {
                return None;
            }
            for j in 0..2 * n {
                aug[col * 2 * n + j] /= pivot;
            }
            for row in 0..n {
                if row != col {
                    let factor = aug[row * 2 * n + col];
                    for j in 0..2 * n {
                        let v = aug[col * 2 * n + j];
                        aug[row * 2 * n + j] -= factor * v;
                    }
                }
            }
        }
        let mut inv = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                inv[i * n + j] = aug[i * 2 * n + n + j];
            }
        }
        Some(inv)
    }
    fn identity(n: usize) -> Vec<f64> {
        let mut m = vec![0.0_f64; n * n];
        for i in 0..n {
            m[i * n + i] = 1.0;
        }
        m
    }
    /// Prediction step: `x = F*x`, `P = F*P*F' + Q`.
    pub fn predict(&mut self) {
        let n = self.n;
        self.x = Self::mat_vec_mul(&self.f, &self.x, n, n);
        let fp = Self::mat_mul(&self.f, &self.p, n, n, n);
        let ft = Self::mat_transpose(&self.f, n, n);
        let fpft = Self::mat_mul(&fp, &ft, n, n, n);
        self.p = Self::mat_add(&fpft, &self.q);
    }
    /// Update step: incorporate measurement `z` (length m).
    pub fn update(&mut self, z: &[f64]) {
        let n = self.n;
        let m = self.m;
        let ht = Self::mat_transpose(&self.h, m, n);
        let hp = Self::mat_mul(&self.h, &self.p, m, n, n);
        let hpht = Self::mat_mul(&hp, &ht, m, n, m);
        let s = Self::mat_add(&hpht, &self.r);
        let s_inv = match Self::mat_inv(&s, m) {
            Some(inv) => inv,
            None => return,
        };
        let pht = Self::mat_mul(&self.p, &ht, n, n, m);
        let k = Self::mat_mul(&pht, &s_inv, n, m, m);
        let hx = Self::mat_vec_mul(&self.h, &self.x, m, n);
        let innov: Vec<f64> = z.iter().zip(hx.iter()).map(|(zi, hi)| zi - hi).collect();
        let k_innov = Self::mat_vec_mul(&k, &innov, n, m);
        for i in 0..n {
            self.x[i] += k_innov[i];
        }
        let kh = Self::mat_mul(&k, &self.h, n, m, n);
        let i_kh = Self::mat_sub(&Self::identity(n), &kh);
        self.p = Self::mat_mul(&i_kh, &self.p, n, n, n);
    }
}
/// Biquad IIR filter (Direct Form II transposed).
///
/// Processes signals using a single second-order section. For higher-order
/// filters, cascade multiple `IirFilter` instances.
pub struct IirFilter {
    /// Biquad coefficients.
    pub coeff: BiquadCoeff,
    /// Delay state w1.
    pub(super) w1: f64,
    /// Delay state w2.
    pub(super) w2: f64,
}
impl IirFilter {
    /// Create a new `IirFilter` from the given biquad coefficients.
    pub fn new(coeff: BiquadCoeff) -> Self {
        Self {
            coeff,
            w1: 0.0,
            w2: 0.0,
        }
    }
    /// Reset the internal state.
    pub fn reset(&mut self) {
        self.w1 = 0.0;
        self.w2 = 0.0;
    }
    /// Process a single sample (Direct Form II transposed).
    pub fn step(&mut self, x: f64) -> f64 {
        let y = self.coeff.b0 * x + self.w1;
        self.w1 = self.coeff.b1 * x - self.coeff.a1 * y + self.w2;
        self.w2 = self.coeff.b2 * x - self.coeff.a2 * y;
        y
    }
    /// Filter an entire signal and return the output vector.
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&x| self.step(x)).collect()
    }
}
/// nth-order Butterworth low-pass filter.
///
/// Designed via the bilinear transform. The filter is implemented as a
/// cascade of biquad sections (and optionally one first-order section for
/// odd orders).
pub struct ButterworthLowPass {
    /// Cascade of biquad sections.
    pub(super) sections: Vec<IirFilter>,
    /// Optional first-order section for odd filter orders.
    pub(super) first_order: Option<(f64, f64, f64)>,
    /// First-order delay state.
    pub(super) fo_state: (f64, f64),
}
impl ButterworthLowPass {
    /// Design an nth-order Butterworth low-pass filter.
    ///
    /// # Arguments
    /// * `order` — Filter order (1..=8 supported).
    /// * `cutoff_hz` — Cutoff frequency in Hz.
    /// * `sample_rate_hz` — Sample rate in Hz.
    pub fn new(order: usize, cutoff_hz: f64, sample_rate_hz: f64) -> Self {
        let wc = 2.0 * sample_rate_hz * (PI * cutoff_hz / sample_rate_hz).tan();
        let n_pairs = order / 2;
        let has_odd = order % 2 == 1;
        let mut sections = Vec::new();
        for k in 0..n_pairs {
            let theta = PI * (2 * (k + 1) + order - 1) as f64 / (2 * order) as f64;
            let re = theta.cos();
            let im = theta.sin();
            let a_s2 = 1.0 / (wc * wc);
            let a_s1 = -2.0 * re / wc;
            let k_bl = 2.0 * sample_rate_hz;
            let k2 = k_bl * k_bl;
            let d = a_s2 * k2 + a_s1 * k_bl + 1.0;
            let b0 = 1.0 / d;
            let b1 = 2.0 * b0;
            let b2 = b0;
            let ia1 = (2.0 * (1.0 - a_s2 * k2)) / d;
            let ia2 = (a_s2 * k2 - a_s1 * k_bl + 1.0) / d;
            let coeff = BiquadCoeff::new(b0, b1, b2, ia1, ia2);
            sections.push(IirFilter::new(coeff));
            let _ = im;
        }
        let first_order = if has_odd {
            let alpha = wc / (wc + 2.0 * sample_rate_hz);
            let a1 = 1.0 - 2.0 * alpha;
            Some((alpha, alpha, -a1))
        } else {
            None
        };
        Self {
            sections,
            first_order,
            fo_state: (0.0, 0.0),
        }
    }
    /// Reset all internal filter states.
    pub fn reset(&mut self) {
        for s in self.sections.iter_mut() {
            s.reset();
        }
        self.fo_state = (0.0, 0.0);
    }
    /// Process a single sample through the filter cascade.
    pub fn step(&mut self, x: f64) -> f64 {
        let mut y = x;
        for s in self.sections.iter_mut() {
            y = s.step(y);
        }
        if let Some((b0, b1, a1)) = self.first_order {
            let x_prev = self.fo_state.0;
            let y_prev = self.fo_state.1;
            let y_new = b0 * y + b1 * x_prev + a1 * y_prev;
            self.fo_state = (y, y_new);
            y = y_new;
        }
        y
    }
    /// Filter an entire signal and return the output.
    pub fn apply(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&x| self.step(x)).collect()
    }
}
/// Compute the analytic signal and its envelope via the Hilbert transform.
///
/// Uses the FFT-based approach:
/// 1. Compute FFT(x)
/// 2. Zero negative frequencies, double positive frequencies
/// 3. IFFT to obtain the analytic signal
/// 4. Return the instantaneous envelope `|analytic signal|`
pub struct HilbertTransform;
impl HilbertTransform {
    /// Compute the analytic signal of `x`.
    ///
    /// Returns a complex vector of the same length as `x` (zero-padded to next
    /// power of 2 internally, then truncated to `x.len()`).
    pub fn analytic_signal(x: &[f64]) -> Vec<Complex> {
        let n_orig = x.len();
        if n_orig == 0 {
            return Vec::new();
        }
        let n = next_pow2(n_orig);
        let mut spec: Vec<Complex> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
        spec.resize(n, Complex::zero());
        fft_inplace_impl(&mut spec, false);
        let half = n / 2;
        for k in 1..half {
            spec[k] = spec[k] * Complex::new(2.0, 0.0);
        }
        for k in (half + 1)..n {
            spec[k] = Complex::zero();
        }
        fft_inplace_impl(&mut spec, true);
        let scale = 1.0 / n as f64;
        spec.iter_mut().for_each(|c| {
            c.re *= scale;
            c.im *= scale;
        });
        spec.truncate(n_orig);
        spec
    }
    /// Compute the instantaneous envelope (magnitude of the analytic signal).
    pub fn envelope(x: &[f64]) -> Vec<f64> {
        Self::analytic_signal(x)
            .iter()
            .map(|c| (c.re * c.re + c.im * c.im).sqrt())
            .collect()
    }
}
