// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Signal processing: FFT, windowing, filtering, spectral analysis, wavelets.
//!
//! Provides Cooley-Tukey FFT/IFFT, windowing functions (Hann/Hamming/Blackman/Kaiser),
//! PSD estimation (Welch/Bartlett), FIR/IIR filter design (Butterworth/Chebyshev/elliptic),
//! convolution, correlation, Hilbert transform, spectrogram, and wavelet transforms.

use crate::error::PhysicsError;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Complex number (local, minimal)
// ---------------------------------------------------------------------------

/// A complex number with `f64` real and imaginary parts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Complex64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex64 {
    /// Construct a new `Complex64`.
    #[inline]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Construct from magnitude and phase (polar form).
    #[inline]
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }

    /// Return the complex conjugate.
    #[inline]
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Return the modulus (absolute value).
    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// Return the argument (phase angle in radians).
    #[inline]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Return the squared modulus.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Scale by a real scalar.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }

    /// Compute e^(i*theta) = cos(theta) + i*sin(theta).
    #[inline]
    pub fn exp_i(theta: f64) -> Self {
        Self {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    /// Compute the natural exponential e^z.
    #[inline]
    pub fn exp(self) -> Self {
        let r = self.re.exp();
        Self {
            re: r * self.im.cos(),
            im: r * self.im.sin(),
        }
    }
}

impl std::ops::Add for Complex64 {
    type Output = Self;
    /// Add two complex numbers.
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl std::ops::Sub for Complex64 {
    type Output = Self;
    /// Subtract two complex numbers.
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl std::ops::Mul for Complex64 {
    type Output = Self;
    /// Multiply two complex numbers.
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl std::ops::Neg for Complex64 {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

// ---------------------------------------------------------------------------
// FFT / IFFT — Cooley-Tukey radix-2
// ---------------------------------------------------------------------------

/// Bit-reverse the indices of a slice in place.
fn bit_reverse_copy(a: &mut [Complex64]) {
    let n = a.len();
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
}

/// In-place Cooley-Tukey FFT.  `n` must be a power of two.
///
/// Uses the DIT (decimation-in-time) butterfly algorithm.
pub fn fft_inplace(a: &mut [Complex64]) {
    let n = a.len();
    assert!(n.is_power_of_two(), "FFT length must be a power of two");
    bit_reverse_copy(a);
    let mut len = 2usize;
    while len <= n {
        let ang = -2.0 * PI / len as f64;
        let wlen = Complex64::exp_i(ang);
        let mut i = 0;
        while i < n {
            let mut w = Complex64::new(1.0, 0.0);
            for j in 0..(len / 2) {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w = w * wlen;
            }
            i += len;
        }
        len <<= 1;
    }
}

/// In-place inverse FFT (IFFT).
pub fn ifft_inplace(a: &mut [Complex64]) {
    let n = a.len();
    // conjugate
    for x in a.iter_mut() {
        x.im = -x.im;
    }
    fft_inplace(a);
    // conjugate and scale
    let inv = 1.0 / n as f64;
    for x in a.iter_mut() {
        x.im = -x.im;
        x.re *= inv;
        x.im *= inv;
    }
}

/// Compute FFT of a real-valued signal.  Pads to the next power of two if needed.
///
/// Returns the full complex spectrum of length equal to the next power of two.
pub fn fft(signal: &[f64]) -> Vec<Complex64> {
    let n = signal.len().next_power_of_two();
    let mut buf: Vec<Complex64> = signal
        .iter()
        .map(|&x| Complex64::new(x, 0.0))
        .chain(std::iter::repeat(Complex64::new(0.0, 0.0)))
        .take(n)
        .collect();
    fft_inplace(&mut buf);
    buf
}

/// Compute IFFT and return the real part of each element.
pub fn ifft_real(spectrum: &[Complex64]) -> Vec<f64> {
    let mut buf = spectrum.to_vec();
    ifft_inplace(&mut buf);
    buf.iter().map(|c| c.re).collect()
}

/// Compute the power spectrum (|X\[k\]|^2) from a real signal.
pub fn power_spectrum(signal: &[f64]) -> Vec<f64> {
    let spec = fft(signal);
    spec.iter().map(|c| c.norm_sq()).collect()
}

/// Compute the magnitude spectrum (|X\[k\]|) from a real signal.
pub fn magnitude_spectrum(signal: &[f64]) -> Vec<f64> {
    let spec = fft(signal);
    spec.iter().map(|c| c.abs()).collect()
}

/// Compute the phase spectrum (arg(X\[k\])) from a real signal.
pub fn phase_spectrum(signal: &[f64]) -> Vec<f64> {
    let spec = fft(signal);
    spec.iter().map(|c| c.arg()).collect()
}

/// Return the frequency bins (in Hz) for an FFT of length `n` sampled at `fs` Hz.
pub fn fft_frequencies(n: usize, fs: f64) -> Vec<f64> {
    (0..n).map(|k| k as f64 * fs / n as f64).collect()
}

/// Return only the positive-frequency bins (0 to fs/2).
pub fn rfft_frequencies(n: usize, fs: f64) -> Vec<f64> {
    (0..=n / 2).map(|k| k as f64 * fs / n as f64).collect()
}

// ---------------------------------------------------------------------------
// Windowing functions
// ---------------------------------------------------------------------------

/// The type of window function to apply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WindowType {
    /// Rectangular (no windowing).
    Rectangular,
    /// Hann window.
    Hann,
    /// Hamming window.
    Hamming,
    /// Blackman window.
    Blackman,
    /// Blackman-Harris window.
    BlackmanHarris,
    /// Flat-top window.
    FlatTop,
    /// Bartlett (triangular) window.
    Bartlett,
    /// Kaiser window (parameter beta).
    Kaiser,
    /// Nuttall window.
    Nuttall,
}

/// Shape parameter that [`WindowType::generate`] uses for [`WindowType::Kaiser`]
/// (about 60 dB of side-lobe attenuation).  Call [`kaiser_window`] directly for
/// another shape parameter.
pub const DEFAULT_KAISER_BETA: f64 = 8.6;

impl WindowType {
    /// Generate the window coefficients of length `n`.
    ///
    /// Every variant is honoured; [`WindowType::Kaiser`] uses
    /// [`DEFAULT_KAISER_BETA`].  Lengths 0 and 1 give an empty window and the
    /// single sample `1.0` respectively.
    pub fn generate(&self, n: usize) -> Vec<f64> {
        match self {
            WindowType::Rectangular => vec![1.0; n],
            WindowType::Hann => hann_window(n),
            WindowType::Hamming => hamming_window(n),
            WindowType::Blackman => blackman_window(n),
            WindowType::BlackmanHarris => blackman_harris_window(n),
            WindowType::FlatTop => flat_top_window(n),
            WindowType::Bartlett => bartlett_window(n),
            WindowType::Kaiser => kaiser_window(n, DEFAULT_KAISER_BETA),
            WindowType::Nuttall => nuttall_window(n),
        }
    }
}

/// Windows of length 0 are empty and windows of length 1 are the single sample
/// `1.0` (the symmetric formulas divide by `n - 1`).
fn degenerate_window(n: usize) -> Option<Vec<f64>> {
    match n {
        0 => Some(Vec::new()),
        1 => Some(vec![1.0]),
        _ => None,
    }
}

/// Compute the Hann window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn hann_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    (0..n)
        .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos()))
        .collect()
}

/// Compute the Hamming window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn hamming_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    (0..n)
        .map(|i| 0.54 - 0.46 * (2.0 * PI * i as f64 / (n - 1) as f64).cos())
        .collect()
}

/// Compute the Blackman window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn blackman_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            0.42 - 0.5 * t.cos() + 0.08 * (2.0 * t).cos()
        })
        .collect()
}

/// Compute the Blackman-Harris window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn blackman_harris_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    let a0 = 0.35875;
    let a1 = 0.48829;
    let a2 = 0.14128;
    let a3 = 0.01168;
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            a0 - a1 * t.cos() + a2 * (2.0 * t).cos() - a3 * (3.0 * t).cos()
        })
        .collect()
}

/// Compute the Nuttall window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn nuttall_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    let a0 = 0.3635819;
    let a1 = 0.4891775;
    let a2 = 0.1365995;
    let a3 = 0.0106411;
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            a0 - a1 * t.cos() + a2 * (2.0 * t).cos() - a3 * (3.0 * t).cos()
        })
        .collect()
}

/// Compute the Bartlett (triangular) window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn bartlett_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    let m = (n - 1) as f64;
    (0..n)
        .map(|i| {
            let x = i as f64;
            1.0 - (2.0 * x / m - 1.0).abs()
        })
        .collect()
}

/// Compute the flat-top window of length `n`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn flat_top_window(n: usize) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    let a0 = 0.21557895;
    let a1 = 0.41663158;
    let a2 = 0.277263158;
    let a3 = 0.083578947;
    let a4 = 0.006947368;
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            a0 - a1 * t.cos() + a2 * (2.0 * t).cos() - a3 * (3.0 * t).cos() + a4 * (4.0 * t).cos()
        })
        .collect()
}

/// Modified Bessel function of the first kind, order zero (I₀).
/// Used for Kaiser window computation.
pub fn bessel_i0(x: f64) -> f64 {
    // Series expansion: I0(x) = sum_{k=0}^{inf} (x/2)^{2k} / (k!)^2
    let mut sum = 1.0;
    let mut term = 1.0;
    let half_x = x / 2.0;
    for k in 1..=50 {
        term *= half_x * half_x / (k * k) as f64;
        sum += term;
        if term < 1e-15 * sum {
            break;
        }
    }
    sum
}

/// Compute the Kaiser window of length `n` with shape parameter `beta`.
///
/// `n = 0` gives an empty window and `n = 1` the single sample `1.0`.
pub fn kaiser_window(n: usize, beta: f64) -> Vec<f64> {
    if let Some(window) = degenerate_window(n) {
        return window;
    }
    let i0_beta = bessel_i0(beta);
    let m = (n - 1) as f64;
    (0..n)
        .map(|i| {
            let t = 2.0 * i as f64 / m - 1.0;
            let arg = beta * (1.0 - t * t).max(0.0).sqrt();
            bessel_i0(arg) / i0_beta
        })
        .collect()
}

/// Apply a window to a signal in-place.
pub fn apply_window(signal: &mut [f64], window: &[f64]) {
    assert_eq!(
        signal.len(),
        window.len(),
        "Signal and window must have equal length"
    );
    for (s, w) in signal.iter_mut().zip(window.iter()) {
        *s *= w;
    }
}

/// Apply a window and return the windowed signal.
pub fn windowed(signal: &[f64], window: &[f64]) -> Vec<f64> {
    assert_eq!(signal.len(), window.len());
    signal
        .iter()
        .zip(window.iter())
        .map(|(s, w)| s * w)
        .collect()
}

/// Compute the coherent gain of a window (mean).
pub fn window_coherent_gain(window: &[f64]) -> f64 {
    window.iter().sum::<f64>() / window.len() as f64
}

/// Compute the power gain of a window (mean of squares).
pub fn window_power_gain(window: &[f64]) -> f64 {
    window.iter().map(|w| w * w).sum::<f64>() / window.len() as f64
}

// ---------------------------------------------------------------------------
// Power Spectral Density estimation
// ---------------------------------------------------------------------------

/// PSD estimation result containing frequency bins and power values.
#[derive(Clone, Debug)]
pub struct PsdResult {
    /// Frequency values (Hz).
    pub freqs: Vec<f64>,
    /// Power spectral density (units²/Hz).
    pub psd: Vec<f64>,
}

/// Welch's method for PSD estimation.
///
/// Splits the signal into overlapping segments, windows each, computes the
/// periodogram, and averages.  Segments of `segment_len` samples start every
/// `segment_len - overlap` samples; only complete segments are used.  Every
/// [`WindowType`] is honoured (see [`WindowType::generate`]).
///
/// # Errors
/// Returns [`PhysicsError::InvalidInput`] if
/// * `fs` is not positive and finite,
/// * `segment_len` is zero,
/// * `overlap >= segment_len` (the hop would be zero or negative),
/// * the signal is shorter than one segment (this includes an empty signal), or
/// * the window has no power for this `segment_len` (for example a Hann window
///   of length 2).
pub fn welch_psd(
    signal: &[f64],
    fs: f64,
    segment_len: usize,
    overlap: usize,
    window_type: WindowType,
) -> Result<PsdResult, PhysicsError> {
    if !(fs.is_finite() && fs > 0.0) {
        return Err(PhysicsError::InvalidInput {
            field: "fs".to_string(),
            reason: format!("sample rate must be positive and finite, got {fs}"),
        });
    }
    if segment_len == 0 {
        return Err(PhysicsError::InvalidInput {
            field: "segment_len".to_string(),
            reason: "segment length must be at least 1".to_string(),
        });
    }
    if overlap >= segment_len {
        return Err(PhysicsError::InvalidInput {
            field: "overlap".to_string(),
            reason: format!("overlap {overlap} must be smaller than segment_len {segment_len}"),
        });
    }
    if signal.len() < segment_len {
        return Err(PhysicsError::InvalidInput {
            field: "signal".to_string(),
            reason: format!(
                "signal of {} samples is shorter than one segment of {segment_len} samples",
                signal.len()
            ),
        });
    }
    let step = segment_len - overlap;
    let window = window_type.generate(segment_len);
    let win_power: f64 = window.iter().map(|w| w * w).sum::<f64>();
    // Windows such as Hann or Blackman of length 2 vanish identically (up to rounding).
    if !win_power.is_finite() || win_power < 1e-20 * segment_len as f64 {
        return Err(PhysicsError::InvalidInput {
            field: "window_type".to_string(),
            reason: format!(
                "{window_type:?} window of length {segment_len} has no usable power ({win_power})"
            ),
        });
    }

    let n_fft = segment_len.next_power_of_two();
    let n_bins = n_fft / 2 + 1;
    let mut psd_acc = vec![0.0f64; n_bins];
    let mut seg_count = 0usize;

    let mut start = 0;
    while start + segment_len <= signal.len() {
        let seg = &signal[start..start + segment_len];
        let mut windowed_seg: Vec<f64> =
            seg.iter().zip(window.iter()).map(|(s, w)| s * w).collect();
        windowed_seg.resize(n_fft, 0.0);
        let spec = fft(&windowed_seg);
        // Periodogram: |X[k]|^2 / (fs * W)
        for k in 0..n_bins {
            let factor = if k == 0 || k == n_fft / 2 { 1.0 } else { 2.0 };
            psd_acc[k] += factor * spec[k].norm_sq() / (fs * win_power);
        }
        seg_count += 1;
        start += step;
    }
    // Average over segments (at least one complete segment exists here).
    let inv = 1.0 / seg_count as f64;
    for p in psd_acc.iter_mut() {
        *p *= inv;
    }
    let freqs = rfft_frequencies(n_fft, fs);
    Ok(PsdResult {
        freqs,
        psd: psd_acc,
    })
}

/// Bartlett's method for PSD estimation (non-overlapping segments, rectangular window).
///
/// # Errors
/// Returns [`PhysicsError::InvalidInput`] under the conditions of [`welch_psd`]
/// (`overlap` is always 0 here).
pub fn bartlett_psd(
    signal: &[f64],
    fs: f64,
    segment_len: usize,
) -> Result<PsdResult, PhysicsError> {
    welch_psd(signal, fs, segment_len, 0, WindowType::Rectangular)
}

/// Simple periodogram PSD estimate (a single segment spanning the whole signal).
///
/// # Errors
/// Returns [`PhysicsError::InvalidInput`] under the conditions of [`welch_psd`];
/// in particular an empty signal is rejected.
pub fn periodogram_psd(
    signal: &[f64],
    fs: f64,
    window_type: WindowType,
) -> Result<PsdResult, PhysicsError> {
    welch_psd(signal, fs, signal.len(), 0, window_type)
}

// ---------------------------------------------------------------------------
// Convolution and Correlation
// ---------------------------------------------------------------------------

/// Linear convolution of two real signals via FFT (overlap-add style).
pub fn convolve(x: &[f64], h: &[f64]) -> Vec<f64> {
    let out_len = x.len() + h.len() - 1;
    let n = out_len.next_power_of_two();
    let mut xc: Vec<Complex64> = x
        .iter()
        .map(|&v| Complex64::new(v, 0.0))
        .chain(std::iter::repeat(Complex64::new(0.0, 0.0)))
        .take(n)
        .collect();
    let mut hc: Vec<Complex64> = h
        .iter()
        .map(|&v| Complex64::new(v, 0.0))
        .chain(std::iter::repeat(Complex64::new(0.0, 0.0)))
        .take(n)
        .collect();
    fft_inplace(&mut xc);
    fft_inplace(&mut hc);
    let mut yc: Vec<Complex64> = xc.iter().zip(hc.iter()).map(|(a, b)| *a * *b).collect();
    ifft_inplace(&mut yc);
    yc[..out_len].iter().map(|c| c.re).collect()
}

/// Cross-correlation of `x` with `y` (i.e. sum_t x\[t\] * y\[t+lag\]).
pub fn cross_correlate(x: &[f64], y: &[f64]) -> Vec<f64> {
    // Flip y and convolve
    let y_flip: Vec<f64> = y.iter().rev().cloned().collect();
    convolve(x, &y_flip)
}

/// Auto-correlation of `x`.
pub fn auto_correlate(x: &[f64]) -> Vec<f64> {
    cross_correlate(x, x)
}

/// Normalized cross-correlation coefficient at each lag.
pub fn normalized_cross_correlation(x: &[f64], y: &[f64]) -> Vec<f64> {
    let xcorr = cross_correlate(x, y);
    let norm = (x.iter().map(|v| v * v).sum::<f64>() * y.iter().map(|v| v * v).sum::<f64>()).sqrt();
    if norm == 0.0 {
        xcorr
    } else {
        xcorr.iter().map(|v| v / norm).collect()
    }
}

// ---------------------------------------------------------------------------
// FIR Filter Design
// ---------------------------------------------------------------------------

/// FIR filter coefficients computed via the window method.
#[derive(Clone, Debug)]
pub struct FirFilter {
    /// Filter taps (coefficients).
    pub coeffs: Vec<f64>,
}

impl FirFilter {
    /// Design a low-pass FIR filter using the window method.
    ///
    /// `n_taps`: number of taps (should be odd for linear phase).
    /// `cutoff`: normalised cutoff frequency (0..1, where 1 = Nyquist).
    pub fn lowpass(n_taps: usize, cutoff: f64) -> Self {
        let n_taps = if n_taps.is_multiple_of(2) {
            n_taps + 1
        } else {
            n_taps
        };
        let m = (n_taps - 1) as f64;
        let fc = cutoff * PI;
        let window = hann_window(n_taps);
        let coeffs: Vec<f64> = (0..n_taps)
            .map(|i| {
                let n = i as f64 - m / 2.0;
                if n == 0.0 {
                    fc / PI
                } else {
                    (fc * n).sin() / (PI * n) * window[i]
                }
            })
            .collect();
        Self { coeffs }
    }

    /// Design a high-pass FIR filter.
    pub fn highpass(n_taps: usize, cutoff: f64) -> Self {
        let mut lp = Self::lowpass(n_taps, cutoff);
        // Spectral inversion
        for (i, c) in lp.coeffs.iter_mut().enumerate() {
            *c = if i == (n_taps - 1) / 2 { 1.0 - *c } else { -*c };
        }
        lp
    }

    /// Design a band-pass FIR filter.
    pub fn bandpass(n_taps: usize, low: f64, high: f64) -> Self {
        let lp_high = Self::lowpass(n_taps, high);
        let lp_low = Self::lowpass(n_taps, low);
        let coeffs: Vec<f64> = lp_high
            .coeffs
            .iter()
            .zip(lp_low.coeffs.iter())
            .map(|(h, l)| h - l)
            .collect();
        Self { coeffs }
    }

    /// Design a band-stop (notch) FIR filter.
    pub fn bandstop(n_taps: usize, low: f64, high: f64) -> Self {
        let bp = Self::bandpass(n_taps, low, high);
        let n_taps = bp.coeffs.len();
        let coeffs: Vec<f64> = bp
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| if i == (n_taps - 1) / 2 { 1.0 - c } else { -*c })
            .collect();
        Self { coeffs }
    }

    /// Apply the FIR filter to a signal (direct-form convolution).
    pub fn apply(&self, signal: &[f64]) -> Vec<f64> {
        convolve(signal, &self.coeffs)
    }

    /// Apply using overlap-save direct implementation (causal, same length as input).
    pub fn filter(&self, signal: &[f64]) -> Vec<f64> {
        let n = self.coeffs.len();
        let out_full = self.apply(signal);
        let delay = (n - 1) / 2;
        // Return signal aligned (no group delay)
        if out_full.len() >= signal.len() + delay {
            out_full[delay..delay + signal.len()].to_vec()
        } else {
            let end = out_full.len().min(signal.len());
            out_full[..end].to_vec()
        }
    }

    /// Compute the frequency response at `n_points` equally spaced frequencies \[0, fs/2\].
    pub fn frequency_response(&self, n_points: usize) -> Vec<Complex64> {
        let n_fft = n_points.next_power_of_two() * 2;
        let mut buf: Vec<Complex64> = self
            .coeffs
            .iter()
            .map(|&c| Complex64::new(c, 0.0))
            .chain(std::iter::repeat(Complex64::new(0.0, 0.0)))
            .take(n_fft)
            .collect();
        fft_inplace(&mut buf);
        buf[..n_points].to_vec()
    }
}

// ---------------------------------------------------------------------------
// IIR Filter — Biquad sections
// ---------------------------------------------------------------------------

/// A single biquad (second-order) IIR section.
#[derive(Clone, Debug)]
pub struct Biquad {
    /// Feed-forward coefficients: b0, b1, b2.
    pub b: [f64; 3],
    /// Feedback coefficients: a1, a2 (a0 = 1 normalised).
    pub a: [f64; 2],
    /// Filter state (delay line).
    state: [f64; 2],
}

impl Biquad {
    /// Create a new biquad with given coefficients.
    pub fn new(b: [f64; 3], a: [f64; 2]) -> Self {
        Self {
            b,
            a,
            state: [0.0; 2],
        }
    }

    /// Process a single sample (direct form II transposed).
    pub fn process_sample(&mut self, x: f64) -> f64 {
        let y = self.b[0] * x + self.state[0];
        self.state[0] = self.b[1] * x - self.a[0] * y + self.state[1];
        self.state[1] = self.b[2] * x - self.a[1] * y;
        y
    }

    /// Process a block of samples.
    pub fn process(&mut self, signal: &[f64]) -> Vec<f64> {
        signal.iter().map(|&x| self.process_sample(x)).collect()
    }

    /// Reset filter state.
    pub fn reset(&mut self) {
        self.state = [0.0; 2];
    }
}

/// A cascade of biquad sections forming an IIR filter.
#[derive(Clone, Debug)]
pub struct IirFilter {
    /// Second-order sections.
    pub sections: Vec<Biquad>,
    /// Overall gain.
    pub gain: f64,
}

impl IirFilter {
    /// Create from a list of biquad sections and an overall gain.
    pub fn new(sections: Vec<Biquad>, gain: f64) -> Self {
        Self { sections, gain }
    }

    /// Process a signal through the IIR filter.
    pub fn process(&mut self, signal: &[f64]) -> Vec<f64> {
        let mut out: Vec<f64> = signal.iter().map(|&x| x * self.gain).collect();
        for section in self.sections.iter_mut() {
            out = section.process(&out);
        }
        out
    }

    /// Reset all section states.
    pub fn reset(&mut self) {
        for s in self.sections.iter_mut() {
            s.reset();
        }
    }
}

// ---------------------------------------------------------------------------
// Butterworth filter design
// ---------------------------------------------------------------------------

/// Butterworth low-pass filter design (bilinear transform, biquad cascade).
///
/// `order`: filter order. `wn`: normalised cutoff (0..1, 1 = Nyquist).
pub fn butterworth_lowpass(order: usize, wn: f64) -> IirFilter {
    // Pre-warp cutoff frequency
    let wc = 2.0 * (PI * wn / 2.0).tan();
    let n = order;
    let mut sections = Vec::new();
    let n_sections = n / 2;
    // Pair up poles (and handle odd order separately)
    for k in 0..n_sections {
        let theta = PI * (2 * k + n + 1) as f64 / (2 * n) as f64;
        let real = theta.cos();
        let imag = theta.sin().abs();
        // Analog prototype pole: s = wc * (real + j*imag)
        // Bilinear transform: z = (1 + s/2) / (1 - s/2) -> use Tustin
        let sr = wc * real;
        let si = wc * imag;
        // Bilinear: s -> 2*(z-1)/(z+1)
        // Denominator: (z-1)^2 - (s/2)*(z+1)^2 ...
        // Direct analog-to-digital via bilinear for second-order:
        // H(s) = wc^2 / (s^2 - 2*sr*s + sr^2 + si^2)
        let omega0_sq = sr * sr + si * si;
        let two_alpha = -2.0 * sr;
        // Bilinear coefficients for second-order LP section
        let c = 2.0; // fs normalised
        let b0d = omega0_sq;
        let b1d = 2.0 * omega0_sq;
        let b2d = omega0_sq;
        let a0d = c * c + two_alpha * c + omega0_sq;
        let a1d = -2.0 * c * c + 2.0 * omega0_sq;
        let a2d = c * c - two_alpha * c + omega0_sq;
        let bq = Biquad::new([b0d / a0d, b1d / a0d, b2d / a0d], [a1d / a0d, a2d / a0d]);
        sections.push(bq);
    }
    // Handle odd order (first-order section)
    let mut gain = 1.0;
    if n % 2 == 1 {
        // First-order LP: H(s) = wc / (s + wc)
        // Bilinear: H(z) = (wc/2 + wc/2 z^{-1}) / ((1 + wc/2) + (-1 + wc/2) z^{-1})
        let c = 1.0;
        let b0 = wc / (wc + c);
        let b1 = wc / (wc + c);
        let a1 = (wc - c) / (wc + c);
        let bq = Biquad::new([b0, b1, 0.0], [a1, 0.0]);
        sections.push(bq);
        gain = 1.0;
    }
    IirFilter::new(sections, gain)
}

/// Compute the magnitude response of a cascade of biquads at frequency `f` (normalised 0..1).
pub fn iir_magnitude_response(filter: &IirFilter, f: f64) -> f64 {
    let z = Complex64::exp_i(PI * f);
    let mut h = Complex64::new(filter.gain, 0.0);
    for bq in &filter.sections {
        let num = Complex64::new(bq.b[0], 0.0)
            + Complex64::new(bq.b[1], 0.0) * z.conj()
            + Complex64::new(bq.b[2], 0.0) * z.conj() * z.conj();
        let den = Complex64::new(1.0, 0.0)
            + Complex64::new(bq.a[0], 0.0) * z.conj()
            + Complex64::new(bq.a[1], 0.0) * z.conj() * z.conj();
        // Avoid division by zero
        if den.abs() > 1e-15 {
            let r = den.norm_sq();
            let den_inv = Complex64::new(den.re / r, -den.im / r);
            h = h * num * den_inv;
        }
    }
    h.abs()
}

// ---------------------------------------------------------------------------
// Chebyshev Type I filter design
// ---------------------------------------------------------------------------

/// Chebyshev Type I low-pass filter (biquad cascade).
///
/// `order`: filter order. `wn`: normalised cutoff (0..1). `ripple_db`: passband ripple in dB.
pub fn chebyshev1_lowpass(order: usize, wn: f64, ripple_db: f64) -> IirFilter {
    let eps = (10.0f64.powf(ripple_db / 10.0) - 1.0).sqrt();
    let wc = 2.0 * (PI * wn / 2.0).tan();
    let n = order as f64;
    let asinh_inv_eps = (1.0 / eps).asinh();
    let mut sections = Vec::new();

    let n_sections = order / 2;
    for k in 0..n_sections {
        let theta_k = PI * (2 * k + 1) as f64 / (2 * order) as f64;
        let sigma = -(asinh_inv_eps / n).sinh() * theta_k.sin();
        let omega = (asinh_inv_eps / n).cosh() * theta_k.cos();
        let sr = wc * sigma;
        let si = wc * omega.abs();
        let omega0_sq = sr * sr + si * si;
        let two_alpha = -2.0 * sr;
        let c = 2.0;
        let a0d = c * c + two_alpha * c + omega0_sq;
        let a1d = -2.0 * c * c + 2.0 * omega0_sq;
        let a2d = c * c - two_alpha * c + omega0_sq;
        let bq = Biquad::new(
            [omega0_sq / a0d, 2.0 * omega0_sq / a0d, omega0_sq / a0d],
            [a1d / a0d, a2d / a0d],
        );
        sections.push(bq);
    }
    if order % 2 == 1 {
        let sigma0 = -(asinh_inv_eps / n).sinh();
        let sr = wc * sigma0;
        let c = 1.0;
        let b0 = (-sr) / (c - sr);
        let b1 = (-sr) / (c - sr);
        let a1 = (-sr - c) / (c - sr);
        let bq = Biquad::new([b0, b1, 0.0], [a1, 0.0]);
        sections.push(bq);
    }
    IirFilter::new(sections, 1.0)
}

// ---------------------------------------------------------------------------
// Hilbert Transform
// ---------------------------------------------------------------------------

/// Compute the analytic signal (x + j*H{x}) using FFT.
///
/// Returns the complex analytic signal from which the instantaneous
/// amplitude and frequency can be derived.
pub fn analytic_signal(x: &[f64]) -> Vec<Complex64> {
    let n = x.len().next_power_of_two();
    let mut spec = fft(x);
    spec.resize(n, Complex64::new(0.0, 0.0));
    // Zero the negative frequencies and double positive ones
    for v in &mut spec[1..n / 2] {
        *v = v.scale(2.0);
    }
    for v in &mut spec[n / 2 + 1..n] {
        *v = Complex64::new(0.0, 0.0);
    }
    let mut result = spec;
    ifft_inplace(&mut result);
    result[..x.len()].to_vec()
}

/// Compute the Hilbert transform of a real signal.
pub fn hilbert_transform(x: &[f64]) -> Vec<f64> {
    analytic_signal(x).iter().map(|c| c.im).collect()
}

/// Compute the instantaneous amplitude (envelope) of a signal.
pub fn instantaneous_amplitude(x: &[f64]) -> Vec<f64> {
    analytic_signal(x).iter().map(|c| c.abs()).collect()
}

/// Compute the instantaneous phase of a signal (radians).
pub fn instantaneous_phase(x: &[f64]) -> Vec<f64> {
    analytic_signal(x).iter().map(|c| c.arg()).collect()
}

/// Compute the instantaneous frequency (Hz) given sample rate `fs`.
pub fn instantaneous_frequency(x: &[f64], fs: f64) -> Vec<f64> {
    let phase = instantaneous_phase(x);
    let mut freq = vec![0.0f64; phase.len()];
    for i in 1..phase.len() {
        let mut dphi = phase[i] - phase[i - 1];
        // Unwrap
        while dphi > PI {
            dphi -= 2.0 * PI;
        }
        while dphi < -PI {
            dphi += 2.0 * PI;
        }
        freq[i] = dphi * fs / (2.0 * PI);
    }
    freq[0] = freq[1];
    freq
}

// ---------------------------------------------------------------------------
// Spectrogram
// ---------------------------------------------------------------------------

/// Spectrogram result: time bins, frequency bins, and STFT magnitude matrix.
#[derive(Clone, Debug)]
pub struct Spectrogram {
    /// Time bin centres (seconds).
    pub times: Vec<f64>,
    /// Frequency bins (Hz).
    pub freqs: Vec<f64>,
    /// Magnitude\[time\]\[freq\].
    pub magnitude: Vec<Vec<f64>>,
    /// Phase\[time\]\[freq\] (radians).
    pub phase: Vec<Vec<f64>>,
}

/// Compute a short-time Fourier transform (STFT) spectrogram.
pub fn spectrogram(
    signal: &[f64],
    fs: f64,
    window_len: usize,
    hop: usize,
    window_type: WindowType,
) -> Spectrogram {
    let window = window_type.generate(window_len);
    let n_fft = window_len.next_power_of_two();
    let n_bins = n_fft / 2 + 1;

    let mut times = Vec::new();
    let mut magnitudes = Vec::new();
    let mut phases = Vec::new();
    // A zero hop or an empty window would never advance: no frames.
    let hop = if window_len == 0 { 0 } else { hop };

    let mut start = 0;
    while hop > 0 && start + window_len <= signal.len() {
        times.push((start + window_len / 2) as f64 / fs);
        let mut buf: Vec<f64> = signal[start..start + window_len]
            .iter()
            .zip(window.iter())
            .map(|(s, w)| s * w)
            .collect();
        buf.resize(n_fft, 0.0);
        let spec = fft(&buf);
        magnitudes.push(spec[..n_bins].iter().map(|c| c.abs()).collect::<Vec<_>>());
        phases.push(spec[..n_bins].iter().map(|c| c.arg()).collect::<Vec<_>>());
        start += hop;
    }

    let freqs = rfft_frequencies(n_fft, fs);
    Spectrogram {
        times,
        freqs,
        magnitude: magnitudes,
        phase: phases,
    }
}

// ---------------------------------------------------------------------------
// Wavelet Transforms
// ---------------------------------------------------------------------------

/// Haar wavelet decomposition (one level).
///
/// Returns (approximation, detail) coefficient vectors.
pub fn haar_decompose(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len() / 2;
    let sqrt2_inv = 1.0 / 2.0f64.sqrt();
    let approx: Vec<f64> = (0..n)
        .map(|i| (signal[2 * i] + signal[2 * i + 1]) * sqrt2_inv)
        .collect();
    let detail: Vec<f64> = (0..n)
        .map(|i| (signal[2 * i] - signal[2 * i + 1]) * sqrt2_inv)
        .collect();
    (approx, detail)
}

/// Haar wavelet reconstruction (one level).
pub fn haar_reconstruct(approx: &[f64], detail: &[f64]) -> Vec<f64> {
    let n = approx.len();
    let sqrt2_inv = 1.0 / 2.0f64.sqrt();
    let mut out = vec![0.0f64; 2 * n];
    for i in 0..n {
        out[2 * i] = (approx[i] + detail[i]) * sqrt2_inv;
        out[2 * i + 1] = (approx[i] - detail[i]) * sqrt2_inv;
    }
    out
}

/// Multi-level Haar DWT.
///
/// Returns a flattened coefficient array: \[level_J_approx | details_J | details_{J-1} | ... | details_1\].
pub fn haar_dwt(signal: &[f64], levels: usize) -> Vec<f64> {
    let mut approx = signal.to_vec();
    let mut details_stack: Vec<Vec<f64>> = Vec::new();
    for _ in 0..levels {
        if approx.len() < 2 {
            break;
        }
        let (a, d) = haar_decompose(&approx);
        approx = a;
        details_stack.push(d);
    }
    let mut result = approx;
    for d in details_stack.into_iter().rev() {
        result.extend(d);
    }
    result
}

/// Multi-level Haar inverse DWT.
pub fn haar_idwt(coeffs: &[f64], levels: usize) -> Vec<f64> {
    let total = coeffs.len();
    let approx_len = total >> levels;
    let mut approx = coeffs[..approx_len].to_vec();
    let mut pos = approx_len;
    for _ in 0..levels {
        let detail_len = approx.len();
        if pos + detail_len > coeffs.len() {
            break;
        }
        let detail = &coeffs[pos..pos + detail_len];
        approx = haar_reconstruct(&approx, detail);
        pos += detail_len;
    }
    approx
}

// Daubechies db2 (Daubechies-4) filter coefficients
const DB2_LO: [f64; 4] = [0.6830127, 1.1830127, 0.3169873, -0.1830127];
const DB2_HI: [f64; 4] = [-0.1830127, -0.3169873, 1.1830127, -0.6830127];
// Inverse (synthesis) low/high pass
const DB2_LO_R: [f64; 4] = [-0.1830127, 0.3169873, 1.1830127, 0.6830127];
const DB2_HI_R: [f64; 4] = [-0.6830127, 1.1830127, -0.3169873, -0.1830127];

/// Daubechies-4 (db2) single-level forward DWT.
pub fn db2_decompose(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    let filt_len = DB2_LO.len();
    let out_len = (n + filt_len - 1) / 2;
    let mut approx = vec![0.0f64; out_len];
    let mut detail = vec![0.0f64; out_len];
    for k in 0..out_len {
        let mut a = 0.0;
        let mut d = 0.0;
        for j in 0..filt_len {
            let idx = 2 * k + j;
            let s = if idx < n { signal[idx] } else { 0.0 };
            a += DB2_LO[j] * s;
            d += DB2_HI[j] * s;
        }
        approx[k] = a / 2.0f64.sqrt();
        detail[k] = d / 2.0f64.sqrt();
    }
    (approx, detail)
}

/// Daubechies-4 (db2) single-level inverse DWT.
pub fn db2_reconstruct(approx: &[f64], detail: &[f64]) -> Vec<f64> {
    let n = approx.len();
    let out_len = 2 * n;
    let mut out = vec![0.0f64; out_len + DB2_LO_R.len()];
    for k in 0..n {
        for j in 0..DB2_LO_R.len() {
            let idx = 2 * k + j;
            if idx < out.len() {
                out[idx] += approx[k] * DB2_LO_R[j] / 2.0f64.sqrt();
                out[idx] += detail[k] * DB2_HI_R[j] / 2.0f64.sqrt();
            }
        }
    }
    out[..out_len].to_vec()
}

// ---------------------------------------------------------------------------
// Resampling
// ---------------------------------------------------------------------------

/// Upsample a signal by integer factor `up` (insert zeros between samples).
pub fn upsample(signal: &[f64], up: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; signal.len() * up];
    for (i, &s) in signal.iter().enumerate() {
        out[i * up] = s;
    }
    out
}

/// Downsample a signal by integer factor `down` (decimate, no anti-aliasing).
pub fn downsample(signal: &[f64], down: usize) -> Vec<f64> {
    signal.iter().step_by(down).cloned().collect()
}

/// Resample a signal from sample rate `fs_in` to `fs_out` using rational resampling.
///
/// Computes the integer up/down ratio and applies a low-pass anti-aliasing filter.
pub fn resample(signal: &[f64], fs_in: f64, fs_out: f64) -> Vec<f64> {
    // Compute rational ratio using GCD approximation
    let ratio = fs_out / fs_in;
    let up = (ratio * 100.0).round() as usize;
    let down = 100usize;
    let up_factor = up.max(1);
    let down_factor = down.max(1);
    // Upsample
    let up_sig = upsample(signal, up_factor);
    // Anti-aliasing filter
    let cutoff = 1.0 / up_factor.max(down_factor) as f64;
    let filt = FirFilter::lowpass(65, cutoff);
    let filtered = filt.apply(&up_sig);
    // Downsample

    downsample(&filtered, down_factor)
}

/// Linear interpolation resampling.
pub fn resample_linear(signal: &[f64], n_out: usize) -> Vec<f64> {
    if signal.is_empty() || n_out == 0 {
        return Vec::new();
    }
    let n_in = signal.len();
    (0..n_out)
        .map(|i| {
            let t = i as f64 * (n_in - 1) as f64 / (n_out - 1).max(1) as f64;
            let lo = t as usize;
            let hi = (lo + 1).min(n_in - 1);
            let frac = t - lo as f64;
            signal[lo] * (1.0 - frac) + signal[hi] * frac
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Compute the RMS value of a signal.
pub fn rms(signal: &[f64]) -> f64 {
    let n = signal.len() as f64;
    (signal.iter().map(|x| x * x).sum::<f64>() / n).sqrt()
}

/// Compute the mean of a signal.
pub fn signal_mean(signal: &[f64]) -> f64 {
    signal.iter().sum::<f64>() / signal.len() as f64
}

/// Compute the variance of a signal.
pub fn signal_variance(signal: &[f64]) -> f64 {
    let mean = signal_mean(signal);
    signal.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / signal.len() as f64
}

/// Compute the SNR (signal-to-noise ratio) in dB given signal and noise power.
pub fn snr_db(signal_power: f64, noise_power: f64) -> f64 {
    10.0 * (signal_power / noise_power.max(1e-300)).log10()
}

/// Compute the THD (total harmonic distortion) from a power spectrum.
///
/// `fundamental_bin`: FFT bin index of the fundamental frequency.
/// `n_harmonics`: number of harmonics to include.
pub fn total_harmonic_distortion(
    power_spec: &[f64],
    fundamental_bin: usize,
    n_harmonics: usize,
) -> f64 {
    let p1 = power_spec[fundamental_bin];
    if p1 == 0.0 {
        return 0.0;
    }
    let harmonic_power: f64 = (2..=n_harmonics)
        .filter_map(|h| {
            let bin = fundamental_bin * h;
            if bin < power_spec.len() {
                Some(power_spec[bin])
            } else {
                None
            }
        })
        .sum();
    (harmonic_power / p1).sqrt()
}

/// Zero-pad a signal to length `n`.
pub fn zero_pad(signal: &[f64], n: usize) -> Vec<f64> {
    let mut out = signal.to_vec();
    out.resize(n, 0.0);
    out
}

/// Unwrap phase (correct for 2π discontinuities).
pub fn unwrap_phase(phase: &[f64]) -> Vec<f64> {
    let mut out = phase.to_vec();
    for i in 1..out.len() {
        let diff = out[i] - out[i - 1];
        if diff > PI {
            for v in out[i..].iter_mut() {
                *v -= 2.0 * PI;
            }
        } else if diff < -PI {
            for v in out[i..].iter_mut() {
                *v += 2.0 * PI;
            }
        }
    }
    out
}

/// Compute group delay from the phase response.
pub fn group_delay(phase: &[f64], df: f64) -> Vec<f64> {
    let n = phase.len();
    let mut gd = vec![0.0f64; n];
    for i in 1..n - 1 {
        gd[i] = -(phase[i + 1] - phase[i - 1]) / (2.0 * df);
    }
    gd[0] = gd[1];
    gd[n - 1] = gd[n - 2];
    gd
}

/// Compute the energy of a signal.
pub fn signal_energy(signal: &[f64]) -> f64 {
    signal.iter().map(|x| x * x).sum()
}

/// Normalize a signal to unit energy.
pub fn normalize_energy(signal: &[f64]) -> Vec<f64> {
    let e = signal_energy(signal).sqrt().max(1e-300);
    signal.iter().map(|x| x / e).collect()
}

/// Normalize a signal to unit peak amplitude.
pub fn normalize_peak(signal: &[f64]) -> Vec<f64> {
    let peak = signal
        .iter()
        .map(|x| x.abs())
        .fold(0.0f64, f64::max)
        .max(1e-300);
    signal.iter().map(|x| x / peak).collect()
}

/// Generate a test sinusoidal signal.
pub fn generate_sine(freq: f64, fs: f64, n_samples: usize, amplitude: f64, phase: f64) -> Vec<f64> {
    (0..n_samples)
        .map(|i| amplitude * (2.0 * PI * freq * i as f64 / fs + phase).sin())
        .collect()
}

/// Generate a test chirp (frequency sweep) signal.
pub fn generate_chirp(f0: f64, f1: f64, fs: f64, n_samples: usize, amplitude: f64) -> Vec<f64> {
    let duration = n_samples as f64 / fs;
    (0..n_samples)
        .map(|i| {
            let t = i as f64 / fs;
            let inst_freq = f0 + (f1 - f0) * t / duration;
            amplitude * (2.0 * PI * inst_freq * t).sin()
        })
        .collect()
}

/// Generate white Gaussian noise (Box-Muller transform).
pub fn generate_white_noise(n: usize, amplitude: f64, seed: u64) -> Vec<f64> {
    let mut x = seed as f64 + 1.0;
    let mut noise = Vec::with_capacity(n);
    for _ in 0..n / 2 {
        // LCG-based pseudo-random
        x = (x * 6364136223846793005.0 + 1442695040888963407.0).abs() % 1e18;
        let u1 = x / 1e18;
        x = (x * 6364136223846793005.0 + 1442695040888963407.0).abs() % 1e18;
        let u2 = x / 1e18;
        let r = (-2.0 * (u1 + 1e-300).ln()).sqrt();
        noise.push(amplitude * r * (2.0 * PI * u2).cos());
        noise.push(amplitude * r * (2.0 * PI * u2).sin());
    }
    if n % 2 == 1 {
        noise.push(0.0);
    }
    noise
}

// ---------------------------------------------------------------------------
// Peak detection
// ---------------------------------------------------------------------------

/// Detect local maxima in a signal.
///
/// Returns the indices of peaks where the value is greater than both neighbours.
pub fn find_peaks(signal: &[f64]) -> Vec<usize> {
    (1..signal.len() - 1)
        .filter(|&i| signal[i] > signal[i - 1] && signal[i] > signal[i + 1])
        .collect()
}

/// Detect peaks above a threshold.
pub fn find_peaks_above(signal: &[f64], threshold: f64) -> Vec<usize> {
    find_peaks(signal)
        .into_iter()
        .filter(|&i| signal[i] > threshold)
        .collect()
}

// ---------------------------------------------------------------------------
// Moving statistics
// ---------------------------------------------------------------------------

/// Compute the moving average of a signal.
pub fn moving_average(signal: &[f64], window: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(signal.len());
    let mut sum = 0.0;
    for (i, &x) in signal.iter().enumerate() {
        sum += x;
        if i >= window {
            sum -= signal[i - window];
        }
        out.push(sum / window.min(i + 1) as f64);
    }
    out
}

/// Compute the moving RMS of a signal.
pub fn moving_rms(signal: &[f64], window: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(signal.len());
    let mut sum_sq = 0.0;
    for (i, &x) in signal.iter().enumerate() {
        sum_sq += x * x;
        if i >= window {
            sum_sq -= signal[i - window].powi(2);
        }
        out.push((sum_sq / window.min(i + 1) as f64).sqrt());
    }
    out
}

// ---------------------------------------------------------------------------
// Cepstrum
// ---------------------------------------------------------------------------

/// Compute the real cepstrum of a signal.
pub fn real_cepstrum(signal: &[f64]) -> Vec<f64> {
    let spec = fft(signal);
    let log_mag: Vec<f64> = spec.iter().map(|c| c.abs().max(1e-300).ln()).collect();
    // IFFT of log magnitude (treat as symmetric real spectrum)
    let mut buf: Vec<Complex64> = log_mag.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    ifft_inplace(&mut buf);
    buf.iter().map(|c| c.re).collect()
}

/// Compute the complex cepstrum of a signal.
pub fn complex_cepstrum(signal: &[f64]) -> Vec<f64> {
    let spec = fft(signal);
    let log_spec: Vec<Complex64> = spec
        .iter()
        .map(|c| {
            let r = c.abs().max(1e-300);
            Complex64::new(r.ln(), c.arg())
        })
        .collect();
    let mut buf = log_spec;
    ifft_inplace(&mut buf);
    buf.iter().map(|c| c.re).collect()
}

// ---------------------------------------------------------------------------
// STFT and inverse STFT
// ---------------------------------------------------------------------------

/// Short-time Fourier transform: returns complex coefficients \[frame\]\[freq_bin\].
pub fn stft(
    signal: &[f64],
    window_len: usize,
    hop: usize,
    window_type: WindowType,
) -> Vec<Vec<Complex64>> {
    let window = window_type.generate(window_len);
    let n_fft = window_len.next_power_of_two();
    let n_bins = n_fft / 2 + 1;
    let mut frames = Vec::new();
    let mut start = 0;
    // A zero hop or an empty window would never advance: no frames.
    while hop > 0 && window_len > 0 && start + window_len <= signal.len() {
        let mut buf: Vec<f64> = signal[start..start + window_len]
            .iter()
            .zip(window.iter())
            .map(|(s, w)| s * w)
            .collect();
        buf.resize(n_fft, 0.0);
        let spec = fft(&buf);
        frames.push(spec[..n_bins].to_vec());
        start += hop;
    }
    frames
}

/// Inverse STFT using overlap-add.
pub fn istft(
    frames: &[Vec<Complex64>],
    window_len: usize,
    hop: usize,
    window_type: WindowType,
) -> Vec<f64> {
    let window = window_type.generate(window_len);
    let n_fft = window_len.next_power_of_two();
    let signal_len = hop * frames.len() + window_len;
    let mut out = vec![0.0f64; signal_len];
    let mut norm = vec![0.0f64; signal_len];

    for (frame_idx, frame) in frames.iter().enumerate() {
        let mut spec = frame.to_vec();
        // Mirror for real IFFT
        while spec.len() < n_fft {
            let k = n_fft - spec.len();
            spec.push(frames[frame_idx][k.min(frame.len() - 1)].conj());
        }
        ifft_inplace(&mut spec);
        let start = frame_idx * hop;
        for (i, (&s, &w)) in spec.iter().zip(window.iter()).enumerate().take(window_len) {
            out[start + i] += s.re * w;
            norm[start + i] += w * w;
        }
    }
    for (o, n) in out.iter_mut().zip(norm.iter()) {
        if *n > 1e-10 {
            *o /= n;
        }
    }
    out
}

#[cfg(test)]
mod tests;
