//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use crate::complex::Complex;
use std::f64::consts::PI;

use super::types::WindowType;

/// Next power of two >= n.
pub(super) fn next_pow2(n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}
/// Bit-reversal permutation in-place.
pub(super) fn bit_reverse_permute(buf: &mut [Complex]) {
    let n = buf.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            buf.swap(i, j);
        }
    }
}
/// Cooley-Tukey radix-2 DIT butterfly in-place.
/// `inverse = true` uses +2π twiddle factor (for IFFT).
pub(super) fn fft_inplace_impl(buf: &mut Vec<Complex>, inverse: bool) {
    let n = buf.len();
    if n <= 1 {
        return;
    }
    bit_reverse_permute(buf);
    let sign = if inverse { 1.0_f64 } else { -1.0_f64 };
    let mut len = 2usize;
    while len <= n {
        let half = len / 2;
        let angle = sign * 2.0 * PI / (len as f64);
        let w_base = Complex::new(angle.cos(), angle.sin());
        let mut k = 0usize;
        while k < n {
            let mut w = Complex::one();
            for l in 0..half {
                let u = buf[k + l];
                let v = buf[k + l + half] * w;
                buf[k + l] = u + v;
                buf[k + l + half] = u - v;
                w = w * w_base;
            }
            k += len;
        }
        len <<= 1;
    }
}
/// Cooley-Tukey radix-2 DIT FFT.
///
/// Pads `input` to the next power of 2 and returns the full two-sided
/// complex spectrum of that length.
pub fn fft(input: &[Complex]) -> Vec<Complex> {
    let n = next_pow2(input.len().max(1));
    let mut buf: Vec<Complex> = input.to_vec();
    buf.resize(n, Complex::zero());
    fft_inplace_impl(&mut buf, false);
    buf
}
/// Inverse FFT.
///
/// Pads `input` to the next power of 2, computes the IFFT, and scales by
/// 1/N. Returns the full complex output.
pub fn ifft(input: &[Complex]) -> Vec<Complex> {
    let n = next_pow2(input.len().max(1));
    let mut buf: Vec<Complex> = input.to_vec();
    buf.resize(n, Complex::zero());
    fft_inplace_impl(&mut buf, true);
    let scale = 1.0 / (n as f64);
    buf.iter()
        .map(|c| Complex::new(c.re * scale, c.im * scale))
        .collect()
}
/// Real-to-complex FFT.
///
/// Converts `input` to complex (zero imaginary part), pads to next power of 2,
/// and returns the full complex spectrum.
pub fn fft_real(input: &[f64]) -> Vec<Complex> {
    let complex_input: Vec<Complex> = input.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fft(&complex_input)
}
/// Single-sided power spectrum: |FFT(x)|² for bins 0..=N/2.
///
/// Returns a vector of length `N/2 + 1` where `N` is the next power of 2
/// >= `input.len()`.
pub fn power_spectrum(input: &[f64]) -> Vec<f64> {
    let spectrum = fft_real(input);
    let n = spectrum.len();
    let half = n / 2 + 1;
    spectrum.iter().take(half).map(|c| c.norm_sq()).collect()
}
/// Frequency (Hz) of the dominant peak in the power spectrum (excluding DC).
pub fn dominant_frequency(input: &[f64], sample_rate: f64) -> f64 {
    let ps = power_spectrum(input);
    let n_padded = next_pow2(input.len().max(1));
    let mut peak_bin = 1usize;
    let mut peak_val = 0.0_f64;
    for (k, &p) in ps.iter().enumerate().skip(1) {
        if p > peak_val {
            peak_val = p;
            peak_bin = k;
        }
    }
    (peak_bin as f64) * sample_rate / (n_padded as f64)
}
/// Linear convolution of `signal` with `kernel`.
///
/// Output length is `signal.len() + kernel.len() - 1`.
pub fn convolve(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() {
        return Vec::new();
    }
    let out_len = signal.len() + kernel.len() - 1;
    let mut out = vec![0.0_f64; out_len];
    for (i, &s) in signal.iter().enumerate() {
        for (j, &k) in kernel.iter().enumerate() {
            out[i + j] += s * k;
        }
    }
    out
}
/// Cross-correlation of `a` with `b`.
///
/// Computes `corr[k] = Σ_n a[n] * b[n - k]` for all valid lags.
/// Output length is `a.len() + b.len() - 1`.
pub fn correlate(a: &[f64], b: &[f64]) -> Vec<f64> {
    let b_rev: Vec<f64> = b.iter().rev().copied().collect();
    convolve(a, &b_rev)
}
/// Causal moving average filter.
///
/// Each output sample is the mean of the current and previous `window-1`
/// input samples. Output has the same length as `signal`.
pub fn moving_average(signal: &[f64], window: usize) -> Vec<f64> {
    if signal.is_empty() || window == 0 {
        return signal.to_vec();
    }
    let w = window.min(signal.len());
    let mut out = Vec::with_capacity(signal.len());
    let mut sum = 0.0_f64;
    for (i, &x) in signal.iter().enumerate() {
        sum += x;
        if i >= w {
            sum -= signal[i - w];
        }
        let count = (i + 1).min(w);
        out.push(sum / count as f64);
    }
    out
}
/// 1-D Gaussian filter.
///
/// Builds a Gaussian kernel of radius `3*sigma` (truncated), normalises it,
/// and convolves it with `signal`. Output length equals `signal.len()` using
/// "same" mode (centre-aligned, zero-padded boundaries).
pub fn gaussian_filter_1d(signal: &[f64], sigma: f64) -> Vec<f64> {
    if sigma <= 0.0 || signal.is_empty() {
        return signal.to_vec();
    }
    let radius = (3.0 * sigma).ceil() as usize;
    let kernel_len = 2 * radius + 1;
    let mut kernel = vec![0.0_f64; kernel_len];
    let denom = 2.0 * sigma * sigma;
    for i in 0..kernel_len {
        let x = i as f64 - radius as f64;
        kernel[i] = (-x * x / denom).exp();
    }
    let k_sum: f64 = kernel.iter().sum();
    for k in kernel.iter_mut() {
        *k /= k_sum;
    }
    let n = signal.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let mut acc = 0.0_f64;
        let mut weight = 0.0_f64;
        for (j, &kv) in kernel.iter().enumerate() {
            let src = i as isize + j as isize - radius as isize;
            if src >= 0 && (src as usize) < n {
                acc += signal[src as usize] * kv;
                weight += kv;
            }
        }
        out[i] = if weight > 0.0 {
            acc / weight
        } else {
            signal[i]
        };
    }
    out
}
/// 1-D median filter with a running sort over `window` samples.
///
/// Output has the same length as `signal` (boundary-replicated padding).
pub fn median_filter_1d(signal: &[f64], window: usize) -> Vec<f64> {
    if signal.is_empty() || window == 0 {
        return signal.to_vec();
    }
    let w = window.min(signal.len());
    let half = w / 2;
    let n = signal.len();
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let start = i.saturating_sub(half);
        let end = (i + half + 1).min(n);
        let mut window_vals: Vec<f64> = signal[start..end].to_vec();
        window_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = window_vals.len() / 2;
        out[i] = if window_vals.len() % 2 == 1 {
            window_vals[mid]
        } else {
            (window_vals[mid - 1] + window_vals[mid]) * 0.5
        };
    }
    out
}
/// Hann window of length `n`.
pub fn hann_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos()))
        .collect()
}
/// Hamming window of length `n`.
pub fn hamming_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|i| 0.54 - 0.46 * (2.0 * PI * i as f64 / (n - 1) as f64).cos())
        .collect()
}
/// Blackman window of length `n`.
pub fn blackman_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            0.42 - 0.5 * t.cos() + 0.08 * (2.0 * t).cos()
        })
        .collect()
}
/// Flat-top window of length `n`.
pub fn flat_top_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    let a0 = 0.21557895_f64;
    let a1 = 0.41663158_f64;
    let a2 = 0.277263158_f64;
    let a3 = 0.083578947_f64;
    let a4 = 0.006947368_f64;
    (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / (n - 1) as f64;
            a0 - a1 * t.cos() + a2 * (2.0 * t).cos() - a3 * (3.0 * t).cos() + a4 * (4.0 * t).cos()
        })
        .collect()
}
/// Root-mean-square of `signal`.
pub fn rms(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let mean_sq = signal.iter().map(|&x| x * x).sum::<f64>() / signal.len() as f64;
    mean_sq.sqrt()
}
/// Signal-to-noise ratio in dB: `10 * log10(|signal|² / |noise|²)`.
pub fn snr_db(signal: &[f64], noise: &[f64]) -> f64 {
    let sig_power: f64 = signal.iter().map(|&x| x * x).sum();
    let noise_power: f64 = noise.iter().map(|&x| x * x).sum();
    if noise_power == 0.0 {
        return f64::INFINITY;
    }
    10.0 * (sig_power / noise_power).log10()
}
/// Number of zero crossings in `signal`.
pub fn zero_crossings(signal: &[f64]) -> usize {
    if signal.len() < 2 {
        return 0;
    }
    signal.windows(2).filter(|w| w[0] * w[1] < 0.0).count()
}
/// Peak-to-peak amplitude of `signal` (max - min).
pub fn peak_to_peak(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }
    let max = signal.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = signal.iter().cloned().fold(f64::INFINITY, f64::min);
    max - min
}
/// Autocorrelation at a single `lag`.
///
/// Returns the normalised autocorrelation coefficient r(lag).
pub fn autocorrelation(signal: &[f64], lag: usize) -> f64 {
    let n = signal.len();
    if n == 0 || lag >= n {
        return 0.0;
    }
    let mean = signal.iter().sum::<f64>() / n as f64;
    let variance: f64 = signal.iter().map(|&x| (x - mean).powi(2)).sum();
    if variance == 0.0 {
        return 1.0;
    }
    let corr: f64 = signal[..n - lag]
        .iter()
        .zip(signal[lag..].iter())
        .map(|(&a, &b)| (a - mean) * (b - mean))
        .sum();
    corr / variance
}
/// Design a windowed-sinc low-pass FIR filter.
///
/// # Arguments
/// * `num_taps`   - Number of filter coefficients (should be odd for type I).
/// * `cutoff_norm`- Normalised cutoff frequency in `(0, 0.5]` (fraction of
///   sampling rate; 0.5 = Nyquist).
/// * `window`     - Window function to apply to the sinc.
///
/// Returns a vector of `num_taps` coefficients.
pub fn fir_windowed_sinc(num_taps: usize, cutoff_norm: f64, window: WindowType) -> Vec<f64> {
    if num_taps == 0 {
        return Vec::new();
    }
    let m = num_taps as isize - 1;
    let win = window.generate(num_taps);
    let mut h: Vec<f64> = (0..num_taps as isize)
        .map(|n| {
            let x = n - m / 2;
            let sinc = if x == 0 {
                2.0 * cutoff_norm
            } else {
                (2.0 * PI * cutoff_norm * x as f64).sin() / (PI * x as f64)
            };
            sinc * win[n as usize]
        })
        .collect();
    let sum: f64 = h.iter().sum();
    if sum.abs() > 1e-15 {
        for c in &mut h {
            *c /= sum;
        }
    }
    h
}
/// Compute the full linear cross-correlation of two real signals.
///
/// Equivalent to `correlate(a, b)` but exposed separately for clarity.
/// Returns a vector of length `len(a) + len(b) - 1`.
///
/// `cc[k]` = Σ_n a\[n\] * b\[n - (k - (len(b)-1))\]  (with zero padding).
pub fn cross_correlate_full(a: &[f64], b: &[f64]) -> Vec<f64> {
    let b_rev: Vec<f64> = b.iter().rev().cloned().collect();
    let na = a.len();
    let nb = b.len();
    let out_len = na + nb - 1;
    let mut out = vec![0.0_f64; out_len];
    for i in 0..na {
        for j in 0..nb {
            out[i + j] += a[i] * b_rev[j];
        }
    }
    out
}
/// Compute Savitzky-Golay smoothing coefficients for a symmetric window.
///
/// Uses a least-squares polynomial fit of degree `poly_order` over a window
/// of `2*half_window + 1` points.  Returns the `2*half_window + 1` coefficients.
///
/// # Arguments
/// * `half_window` - Number of points on each side of the centre.
/// * `poly_order`  - Polynomial order (must be < window size).
pub(super) fn savgol_coeffs(half_window: usize, poly_order: usize) -> Vec<f64> {
    let m = half_window as isize;
    let window = 2 * half_window + 1;
    let deg = poly_order.min(window - 1);
    let rows = window;
    let cols = deg + 1;
    let mut v = vec![0.0_f64; rows * cols];
    for (r, i) in (-m..=m).enumerate() {
        let x = i as f64;
        let mut xpow = 1.0;
        for c in 0..cols {
            v[r * cols + c] = xpow;
            xpow *= x;
        }
    }
    let mut vtv = vec![0.0_f64; cols * cols];
    for i in 0..cols {
        for j in 0..cols {
            let mut s = 0.0;
            for r in 0..rows {
                s += v[r * cols + i] * v[r * cols + j];
            }
            vtv[i * cols + j] = s;
        }
    }
    let mut aug = vec![0.0_f64; cols * 2 * cols];
    for i in 0..cols {
        for j in 0..cols {
            aug[i * 2 * cols + j] = vtv[i * cols + j];
        }
        aug[i * 2 * cols + cols + i] = 1.0;
    }
    for col in 0..cols {
        let mut max_row = col;
        for row in col + 1..cols {
            if aug[row * 2 * cols + col].abs() > aug[max_row * 2 * cols + col].abs() {
                max_row = row;
            }
        }
        for j in 0..2 * cols {
            aug.swap(col * 2 * cols + j, max_row * 2 * cols + j);
        }
        let pivot = aug[col * 2 * cols + col];
        if pivot.abs() < 1e-15 {
            continue;
        }
        for j in 0..2 * cols {
            aug[col * 2 * cols + j] /= pivot;
        }
        for row in 0..cols {
            if row == col {
                continue;
            }
            let factor = aug[row * 2 * cols + col];
            for j in 0..2 * cols {
                let val = aug[col * 2 * cols + j] * factor;
                aug[row * 2 * cols + j] -= val;
            }
        }
    }
    let mut vtv_inv = vec![0.0_f64; cols * cols];
    for i in 0..cols {
        for j in 0..cols {
            vtv_inv[i * cols + j] = aug[i * 2 * cols + cols + j];
        }
    }
    let mut coeffs = vec![0.0_f64; rows];
    for r in 0..rows {
        let mut val = 0.0;
        for j in 0..cols {
            val += v[r * cols + j] * vtv_inv[j * cols];
        }
        coeffs[r] = val;
    }
    coeffs
}
/// Savitzky-Golay smoothing filter.
///
/// Smooths `signal` using a polynomial of degree `poly_order` fitted over a
/// sliding window of `2*half_window + 1` samples.  Endpoints are handled by
/// replicating the signal edge values.
///
/// # Arguments
/// * `signal`      - Input signal.
/// * `half_window` - Half-width of the smoothing window (window = 2*half_window+1).
/// * `poly_order`  - Polynomial order for local fitting (2 or 4 recommended).
pub fn savitzky_golay(signal: &[f64], half_window: usize, poly_order: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return vec![];
    }
    let coeffs = savgol_coeffs(half_window, poly_order);
    let m = half_window as isize;
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let mut val = 0.0;
        for (k, &c) in coeffs.iter().enumerate() {
            let idx = i as isize + k as isize - m;
            let idx = idx.clamp(0, n as isize - 1) as usize;
            val += c * signal[idx];
        }
        out[i] = val;
    }
    out
}
/// Compute a one-sided power spectral density estimate using the periodogram method.
///
/// Applies an optional window (none = rectangular) and returns the one-sided
/// PSD vector of length `N/2 + 1` and the corresponding frequency bins in Hz.
///
/// # Arguments
/// * `signal`      - Input signal.
/// * `sample_rate` - Sampling frequency in Hz.
/// * `window`      - Optional window function (length must match signal length).
///
/// Returns `(psd, freqs)`.
pub fn periodogram_psd(
    signal: &[f64],
    sample_rate: f64,
    window: Option<&[f64]>,
) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    if n == 0 {
        return (vec![], vec![]);
    }
    let windowed: Vec<f64> = if let Some(w) = window {
        signal.iter().zip(w.iter()).map(|(s, wi)| s * wi).collect()
    } else {
        signal.to_vec()
    };
    let w_norm: f64 = if let Some(w) = window {
        w.iter().map(|wi| wi * wi).sum::<f64>()
    } else {
        n as f64
    };
    let complex_in: Vec<Complex> = windowed.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let spectrum = fft(&complex_in);
    let n_out = n / 2 + 1;
    let df = sample_rate / n as f64;
    let mut psd = vec![0.0_f64; n_out];
    let mut freqs = vec![0.0_f64; n_out];
    for k in 0..n_out {
        let power = spectrum[k].norm_sq() / (w_norm * sample_rate);
        psd[k] = if k == 0 || k == n / 2 {
            power
        } else {
            2.0 * power
        };
        freqs[k] = k as f64 * df;
    }
    (psd, freqs)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::signal::BiquadCoeff;
    use crate::signal::ButterworthBandpass;
    use crate::signal::ButterworthLowPass;
    use crate::signal::FirFilter;
    use crate::signal::HilbertTransform;
    use crate::signal::IirFilter;
    use crate::signal::KalmanFilter1D;
    use crate::signal::KalmanFilterND;
    use crate::signal::SignalProcessor;
    use crate::signal::SpectrumPeakDetector;
    use crate::signal::Stft;
    use crate::signal::WelchPsd;
    pub(super) const EPS: f64 = 1e-9;
    pub(super) const LOOSE: f64 = 1e-6;
    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }
    fn close_complex(a: Complex, b: Complex, eps: f64) -> bool {
        close(a.re, b.re, eps) && close(a.im, b.im, eps)
    }
    #[test]
    fn fft_ifft_roundtrip_impulse() {
        let n = 8;
        let impulse: Vec<Complex> = (0..n)
            .map(|i| {
                if i == 0 {
                    Complex::one()
                } else {
                    Complex::zero()
                }
            })
            .collect();
        let spectrum = fft(&impulse);
        let recovered = ifft(&spectrum);
        assert!(
            close(recovered[0].re, 1.0, LOOSE),
            "recovered[0].re = {}",
            recovered[0].re
        );
        for i in 1..n {
            assert!(
                close(recovered[i].re, 0.0, LOOSE),
                "recovered[{i}].re = {}",
                recovered[i].re
            );
        }
    }
    #[test]
    fn fft_ifft_roundtrip_random_like() {
        let signal: Vec<Complex> = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(-1.0, 0.0),
            Complex::new(0.5, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(-2.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.5, 0.0),
        ];
        let recovered = ifft(&fft(&signal));
        for (i, (&orig, rec)) in signal.iter().zip(recovered.iter()).enumerate() {
            assert!(close(orig.re, rec.re, LOOSE), "re mismatch at {i}");
            assert!(close(orig.im, rec.im, LOOSE), "im mismatch at {i}");
        }
    }
    #[test]
    fn fft_real_dc_only() {
        let n = 8;
        let signal = vec![1.0_f64; n];
        let spectrum = fft_real(&signal);
        assert!(
            close(spectrum[0].re, n as f64, LOOSE),
            "DC bin = {}",
            spectrum[0].re
        );
        for k in 1..spectrum.len() {
            assert!(
                spectrum[k].norm() < LOOSE,
                "bin {k} magnitude = {}",
                spectrum[k].norm()
            );
        }
    }
    #[test]
    fn fft_linearity() {
        let a: Vec<Complex> = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];
        let b: Vec<Complex> = vec![
            Complex::new(0.5, 0.0),
            Complex::new(-1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(0.0, 0.0),
        ];
        let alpha = 3.0_f64;
        let beta = -2.0_f64;
        let ab_sum: Vec<Complex> = a
            .iter()
            .zip(b.iter())
            .map(|(&ai, &bi)| Complex::new(alpha * ai.re + beta * bi.re, 0.0))
            .collect();
        let fa = fft(&a);
        let fb = fft(&b);
        let fab = fft(&ab_sum);
        for k in 0..fa.len() {
            let expected = Complex::new(
                alpha * fa[k].re + beta * fb[k].re,
                alpha * fa[k].im + beta * fb[k].im,
            );
            assert!(
                close_complex(fab[k], expected, LOOSE),
                "linearity failed at bin {k}"
            );
        }
    }
    #[test]
    fn fft_parseval_theorem() {
        let signal: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let time_power: f64 = signal.iter().map(|&x| x * x).sum();
        let spectrum = fft_real(&signal);
        let n = spectrum.len() as f64;
        let freq_power: f64 = spectrum.iter().map(|c| c.norm_sq()).sum::<f64>() / n;
        assert!(
            close(time_power, freq_power, 1e-6),
            "Parseval: {time_power} vs {freq_power}"
        );
    }
    #[test]
    fn fft_sine_peak_bin() {
        let n = 16usize;
        let freq_bin = 3usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * freq_bin as f64 * i as f64 / n as f64).sin())
            .collect();
        let ps = power_spectrum(&signal);
        let peak = ps
            .iter()
            .enumerate()
            .skip(1)
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        assert_eq!(
            peak, freq_bin,
            "Expected peak at bin {freq_bin}, got {peak}"
        );
    }
    #[test]
    fn power_spectrum_length() {
        let signal = vec![1.0_f64; 7];
        let ps = power_spectrum(&signal);
        assert_eq!(ps.len(), 5, "power_spectrum length = {}", ps.len());
    }
    #[test]
    fn power_spectrum_nonnegative() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        for &p in power_spectrum(&signal).iter() {
            assert!(p >= 0.0, "negative power {p}");
        }
    }
    #[test]
    fn dominant_frequency_10hz() {
        let sample_rate = 100.0_f64;
        let freq = 10.0_f64;
        let n = 128usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * freq * i as f64 / sample_rate).sin())
            .collect();
        let dom = dominant_frequency(&signal, sample_rate);
        assert!(
            (dom - freq).abs() < 1.0,
            "Expected ~{freq} Hz, got {dom} Hz"
        );
    }
    #[test]
    fn convolve_identity_kernel() {
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let result = convolve(&signal, &[1.0]);
        assert_eq!(result.len(), signal.len());
        for (a, b) in signal.iter().zip(result.iter()) {
            assert!(close(*a, *b, EPS));
        }
    }
    #[test]
    fn convolve_length() {
        let sig = vec![1.0; 5];
        let ker = vec![1.0; 3];
        let out = convolve(&sig, &ker);
        assert_eq!(out.len(), 5 + 3 - 1, "convolve output length");
    }
    #[test]
    fn convolve_box_filter() {
        let s = vec![1.0, 1.0, 1.0];
        let k = vec![1.0, 1.0, 1.0];
        let out = convolve(&s, &k);
        let expected = [1.0, 2.0, 3.0, 2.0, 1.0];
        assert_eq!(out.len(), expected.len());
        for (a, b) in out.iter().zip(expected.iter()) {
            assert!(close(*a, *b, EPS), "{a} vs {b}");
        }
    }
    #[test]
    fn correlate_autocorr_peak_at_zero() {
        let signal = vec![1.0, 2.0, 1.0, 0.0, -1.0];
        let cc = correlate(&signal, &signal);
        let centre = signal.len() - 1;
        let peak_idx = cc
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        assert_eq!(
            peak_idx, centre,
            "Auto-corr peak at {peak_idx}, expected {centre}"
        );
    }
    #[test]
    fn moving_average_constant_signal() {
        let signal = vec![5.0_f64; 20];
        let out = moving_average(&signal, 5);
        for &y in &out {
            assert!(close(y, 5.0, EPS), "Expected 5.0, got {y}");
        }
    }
    #[test]
    fn moving_average_same_length() {
        let signal: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let out = moving_average(&signal, 3);
        assert_eq!(out.len(), signal.len());
    }
    #[test]
    fn moving_average_window_1_identity() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let out = moving_average(&signal, 1);
        for (a, b) in signal.iter().zip(out.iter()) {
            assert!(close(*a, *b, EPS));
        }
    }
    #[test]
    fn gaussian_filter_constant_unchanged() {
        let signal = vec![3.0_f64; 20];
        let out = gaussian_filter_1d(&signal, 1.0);
        for &y in &out {
            assert!(close(y, 3.0, 1e-10), "Expected 3.0, got {y}");
        }
    }
    #[test]
    fn gaussian_filter_same_length() {
        let signal: Vec<f64> = (0..15).map(|i| i as f64).collect();
        let out = gaussian_filter_1d(&signal, 2.0);
        assert_eq!(out.len(), signal.len());
    }
    #[test]
    fn median_filter_constant() {
        let signal = vec![7.0_f64; 15];
        let out = median_filter_1d(&signal, 5);
        for &y in &out {
            assert!(close(y, 7.0, EPS), "Expected 7.0, got {y}");
        }
    }
    #[test]
    fn median_filter_removes_spike() {
        let mut signal = vec![1.0_f64; 11];
        signal[5] = 100.0;
        let out = median_filter_1d(&signal, 3);
        assert!(out[5] < 50.0, "Spike not removed: out[5] = {}", out[5]);
    }
    #[test]
    fn hann_window_endpoints_zero() {
        let w = hann_window(16);
        assert!(close(w[0], 0.0, EPS), "hann[0] = {}", w[0]);
        assert!(
            close(*w.last().unwrap(), 0.0, EPS),
            "hann[last] = {}",
            w.last().unwrap()
        );
    }
    #[test]
    fn hamming_window_length() {
        assert_eq!(hamming_window(32).len(), 32);
    }
    #[test]
    fn blackman_window_endpoints_zero() {
        let w = blackman_window(8);
        assert!(close(w[0], 0.0, 1e-10), "blackman[0] = {}", w[0]);
        assert!(close(*w.last().unwrap(), 0.0, 1e-10));
    }
    #[test]
    fn flat_top_window_length() {
        assert_eq!(flat_top_window(64).len(), 64);
    }
    #[test]
    fn hann_window_symmetry() {
        let n = 17;
        let w = hann_window(n);
        for i in 0..n {
            assert!(close(w[i], w[n - 1 - i], EPS), "asymmetric at {i}");
        }
    }
    #[test]
    fn rms_sine_half_sqrt2() {
        let n = 10_000usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let r = rms(&signal);
        let expected = 1.0 / 2.0_f64.sqrt();
        assert!((r - expected).abs() < 1e-4, "Expected {expected}, got {r}");
    }
    #[test]
    fn rms_constant() {
        let signal = vec![3.0_f64; 100];
        assert!(close(rms(&signal), 3.0, EPS));
    }
    #[test]
    fn snr_db_equal_power() {
        let s = vec![1.0_f64; 10];
        let n = vec![1.0_f64; 10];
        assert!(close(snr_db(&s, &n), 0.0, EPS));
    }
    #[test]
    fn snr_db_zero_noise() {
        let s = vec![1.0_f64; 5];
        let n = vec![0.0_f64; 5];
        assert_eq!(snr_db(&s, &n), f64::INFINITY);
    }
    #[test]
    fn zero_crossings_sine() {
        let n = 1000usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 4.0 * i as f64 / n as f64).sin())
            .collect();
        let zc = zero_crossings(&signal);
        assert!(
            (zc as isize - 8).abs() <= 2,
            "Expected ~8 zero crossings, got {zc}"
        );
    }
    #[test]
    fn peak_to_peak_sine() {
        let n = 10_000usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let p2p = peak_to_peak(&signal);
        assert!((p2p - 2.0).abs() < 1e-3, "Expected ~2.0, got {p2p}");
    }
    #[test]
    fn peak_to_peak_constant() {
        let signal = vec![5.0_f64; 100];
        assert!(close(peak_to_peak(&signal), 0.0, EPS));
    }
    #[test]
    fn autocorrelation_lag0_is_one() {
        let signal: Vec<f64> = (0..50).map(|i| (i as f64 * 0.1).sin()).collect();
        assert!(close(autocorrelation(&signal, 0), 1.0, LOOSE));
    }
    #[test]
    fn autocorrelation_decreases_with_lag() {
        let signal: Vec<f64> = (0..100)
            .map(|i| (i as f64 * 1.234).sin() + (i as f64 * 5.678).sin())
            .collect();
        let r0 = autocorrelation(&signal, 0);
        let r10 = autocorrelation(&signal, 10).abs();
        assert!(r0 >= r10, "autocorr(0)={r0} should >= autocorr(10)={r10}");
    }
    #[test]
    fn fir_filter_impulse_response() {
        let coeffs = vec![1.0_f64, 0.5, 0.25];
        let mut fir = FirFilter::new(coeffs.clone());
        let impulse: Vec<f64> = std::iter::once(1.0)
            .chain(std::iter::repeat_n(0.0, 4))
            .collect();
        let out = fir.apply(&impulse);
        assert!((out[0] - 1.0).abs() < 1e-12, "out[0]={}", out[0]);
        assert!((out[1] - 0.5).abs() < 1e-12, "out[1]={}", out[1]);
        assert!((out[2] - 0.25).abs() < 1e-12, "out[2]={}", out[2]);
    }
    #[test]
    fn fir_filter_constant_signal() {
        let coeffs = vec![0.25_f64; 4];
        let mut fir = FirFilter::new(coeffs);
        let signal = vec![1.0_f64; 20];
        let out = fir.apply(&signal);
        for &y in out.iter().skip(3) {
            assert!((y - 1.0).abs() < 1e-12, "y={}", y);
        }
    }
    #[test]
    fn fir_filter_step_response_length() {
        let mut fir = FirFilter::new(vec![1.0, -1.0]);
        let step: Vec<f64> = std::iter::repeat_n(1.0_f64, 10).collect();
        let out = fir.apply(&step);
        assert_eq!(out.len(), 10);
    }
    #[test]
    fn iir_filter_passthrough_all_pass() {
        let coeff = BiquadCoeff::new(1.0, 0.0, 0.0, 0.0, 0.0);
        let mut iir = IirFilter::new(coeff);
        let signal = vec![1.0_f64, 2.0, 3.0, 4.0];
        let out = iir.apply(&signal);
        for (a, b) in signal.iter().zip(out.iter()) {
            assert!(close(*a, *b, EPS), "{a} vs {b}");
        }
    }
    #[test]
    fn iir_filter_dc_gain() {
        let b0 = 0.25_f64;
        let b1 = 0.5;
        let b2 = 0.25;
        let a1 = -0.3;
        let a2 = 0.1;
        let dc_gain = (b0 + b1 + b2) / (1.0 + a1 + a2);
        let coeff = BiquadCoeff::new(b0, b1, b2, a1, a2);
        let mut iir = IirFilter::new(coeff);
        let n = 2000;
        let signal = vec![1.0_f64; n];
        let out = iir.apply(&signal);
        let last = out[n - 1];
        assert!(
            (last - dc_gain).abs() < 1e-6,
            "DC gain: expected {dc_gain}, got {last}"
        );
    }
    #[test]
    fn butterworth_lp_dc_passes() {
        let mut bp = ButterworthLowPass::new(2, 100.0, 1000.0);
        let signal = vec![1.0_f64; 500];
        let out = bp.apply(&signal);
        let last = out[499];
        assert!((last - 1.0).abs() < 0.05, "DC gain = {last}");
    }
    #[test]
    fn butterworth_lp_attenuates_high_freq() {
        let sample_rate = 1000.0_f64;
        let cutoff = 50.0_f64;
        let test_freq = 400.0_f64;
        let n = 2000;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * test_freq * i as f64 / sample_rate).sin())
            .collect();
        let mut bp = ButterworthLowPass::new(4, cutoff, sample_rate);
        let out = bp.apply(&signal);
        let out_rms = rms(&out[200..]);
        assert!(
            out_rms < 0.3,
            "High-frequency signal not attenuated: rms={out_rms}"
        );
    }
    #[test]
    fn kalman_1d_tracks_constant() {
        let mut kf = KalmanFilter1D::new(0.0, 1.0, 0.01, 0.1);
        let true_val = 5.0_f64;
        let mut estimate = 0.0_f64;
        for _ in 0..200 {
            estimate = kf.step(true_val);
        }
        assert!(
            (estimate - true_val).abs() < 0.1,
            "Kalman 1D: estimate={estimate}"
        );
    }
    #[test]
    fn kalman_1d_uncertainty_decreases() {
        let mut kf = KalmanFilter1D::new(0.0, 10.0, 0.001, 1.0);
        let p0 = kf.p;
        kf.predict();
        kf.update(1.0);
        assert!(kf.p < p0, "Uncertainty should decrease after update");
    }
    #[test]
    fn kalman_nd_constant_velocity_tracking() {
        let n = 2;
        let m = 1;
        let dt = 0.1_f64;
        let f_mat = vec![1.0, dt, 0.0, 1.0];
        let h_mat = vec![1.0, 0.0];
        let q = vec![0.0001, 0.0, 0.0, 0.0001];
        let r = vec![0.5];
        let p0 = vec![1.0, 0.0, 0.0, 1.0];
        let x0 = vec![0.0, 0.0];
        let mut kf = KalmanFilterND::new(n, m, x0, p0, f_mat, h_mat, q, r);
        let v_true = 1.0_f64;
        for step in 0..50 {
            kf.predict();
            let pos_meas = v_true * (step + 1) as f64 * dt;
            kf.update(&[pos_meas]);
        }
        assert!(
            (kf.x[1] - v_true).abs() < 0.2,
            "Velocity estimate = {}",
            kf.x[1]
        );
    }
    #[test]
    fn hilbert_envelope_sine() {
        let n = 128;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 4.0 * i as f64 / n as f64).sin())
            .collect();
        let env = HilbertTransform::envelope(&signal);
        let mid_rms: f64 = {
            let mid = &env[16..112];
            (mid.iter().map(|v| (v - 1.0) * (v - 1.0)).sum::<f64>() / mid.len() as f64).sqrt()
        };
        assert!(mid_rms < 0.1, "Hilbert envelope error: mid_rms={mid_rms}");
    }
    #[test]
    fn hilbert_envelope_nonnegative() {
        let signal: Vec<f64> = (0..64).map(|i| (i as f64 * 0.5).sin()).collect();
        let env = HilbertTransform::envelope(&signal);
        for &v in &env {
            assert!(v >= 0.0, "Envelope should be non-negative, got {v}");
        }
    }
    #[test]
    fn stft_frame_count() {
        let signal = vec![0.0_f64; 100];
        let stft = Stft::new(32, 16);
        let frames = stft.compute(&signal);
        assert_eq!(frames.len(), 5, "Frame count = {}", frames.len());
    }
    #[test]
    fn stft_frame_length() {
        let signal = vec![1.0_f64; 64];
        let stft = Stft::new(16, 8);
        let frames = stft.compute(&signal);
        for f in &frames {
            assert_eq!(f.len(), 16, "Frame length = {}", f.len());
        }
    }
    #[test]
    fn stft_magnitude_spectrogram_nonneg() {
        let signal: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
        let stft = Stft::new(64, 32);
        let mag = stft.magnitude_spectrogram(&signal);
        for frame in &mag {
            for &v in frame {
                assert!(v >= 0.0, "Magnitude should be non-negative, got {v}");
            }
        }
    }
    #[test]
    fn stft_dc_signal() {
        let signal = vec![1.0_f64; 128];
        let stft = Stft::new(32, 16);
        let frames = stft.compute(&signal);
        if !frames.is_empty() {
            let dc_mag = frames[0][0].norm();
            let other_max = frames[0][1..]
                .iter()
                .map(|c| c.norm())
                .fold(0.0_f64, f64::max);
            assert!(
                dc_mag > other_max,
                "DC energy should dominate: dc={dc_mag}, max_other={other_max}"
            );
        }
    }
    #[test]
    fn fir_lowpass_frequency_selectivity() {
        let n_taps = 20;
        let coeffs = vec![1.0 / n_taps as f64; n_taps];
        let mut fir = FirFilter::new(coeffs);
        let n = 500;
        let low_freq_sig: Vec<f64> = (0..n).map(|i| (2.0 * PI * 0.01 * i as f64).sin()).collect();
        let out_low = fir.apply(&low_freq_sig);
        fir.reset();
        let high_freq_sig: Vec<f64> = (0..n).map(|i| (2.0 * PI * 0.4 * i as f64).sin()).collect();
        let out_high = fir.apply(&high_freq_sig);
        let rms_low = rms(&out_low[n_taps..]);
        let rms_high = rms(&out_high[n_taps..]);
        assert!(
            rms_low > rms_high * 2.0,
            "Low-pass selectivity: rms_low={rms_low}, rms_high={rms_high}"
        );
    }
    #[test]
    fn welch_psd_length() {
        let signal: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
        let psd = WelchPsd::new(64, 32, WindowType::Hann).compute(&signal);
        assert!(!psd.is_empty(), "Welch PSD should be non-empty");
    }
    #[test]
    fn welch_psd_nonnegative() {
        let signal: Vec<f64> = (0..512).map(|i| (i as f64 * 0.2).sin()).collect();
        let psd = WelchPsd::new(128, 64, WindowType::Hamming).compute(&signal);
        for &p in &psd {
            assert!(p >= 0.0, "Welch PSD should be non-negative, got {p}");
        }
    }
    #[test]
    fn welch_psd_dc_signal() {
        let signal = vec![1.0_f64; 256];
        let psd = WelchPsd::new(64, 32, WindowType::Hann).compute(&signal);
        if psd.len() > 1 {
            assert!(
                psd[0]
                    > *psd[1..]
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap_or(&0.0),
                "DC power should dominate"
            );
        }
    }
    #[test]
    fn fir_windowed_sinc_length() {
        let h = fir_windowed_sinc(31, 0.2, WindowType::Hann);
        assert_eq!(h.len(), 31, "FIR filter should have exactly 31 taps");
    }
    #[test]
    fn fir_windowed_sinc_dc_gain_near_one() {
        let h = fir_windowed_sinc(63, 0.1, WindowType::Hamming);
        let dc_gain: f64 = h.iter().sum();
        assert!((dc_gain - 1.0).abs() < 0.05, "DC gain = {dc_gain}");
    }
    #[test]
    fn fir_windowed_sinc_symmetry() {
        let n = 33;
        let h = fir_windowed_sinc(n, 0.2, WindowType::Blackman);
        for i in 0..n / 2 {
            assert!((h[i] - h[n - 1 - i]).abs() < 1e-10, "asymmetric at {i}");
        }
    }
    #[test]
    fn peak_detector_finds_dominant() {
        let n = 128usize;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 8.0 * i as f64 / n as f64).sin())
            .collect();
        let ps = power_spectrum(&signal);
        let peaks = SpectrumPeakDetector::find_peaks(&ps, 3, ps[8] * 0.5);
        assert!(peaks.contains(&8), "peak at bin 8 not found: {:?}", peaks);
    }
    #[test]
    fn peak_detector_no_peaks_in_flat() {
        let ps = vec![1.0_f64; 32];
        let peaks = SpectrumPeakDetector::find_peaks(&ps, 3, 2.0);
        assert!(
            peaks.is_empty(),
            "should find no peaks in flat spectrum above threshold"
        );
    }
    #[test]
    fn cross_correlation_lag_detection() {
        let n = 64;
        let signal_a: Vec<f64> = (0..n).map(|i| (i as f64 * 0.5).sin()).collect();
        let delay = 5usize;
        let signal_b: Vec<f64> = (0..n)
            .map(|i| if i >= delay { signal_a[i - delay] } else { 0.0 })
            .collect();
        let cc = cross_correlate_full(&signal_a, &signal_b);
        let centre = n - 1;
        let peak_idx = cc
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        let expected_peak = centre as isize - delay as isize;
        let detected_lag = peak_idx as isize;
        assert!(
            (detected_lag - expected_peak).abs() <= 2,
            "Expected peak near {expected_peak}, got {detected_lag}"
        );
    }
    #[test]
    fn cross_correlation_self_is_autocorrelation() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let cc = cross_correlate_full(&signal, &signal);
        let ac = correlate(&signal, &signal);
        assert_eq!(cc.len(), ac.len());
        for (a, b) in cc.iter().zip(ac.iter()) {
            assert!((a - b).abs() < 1e-10, "cc={a} vs ac={b}");
        }
    }
    #[test]
    fn butterworth_bandpass_passes_center() {
        let fs = 1000.0_f64;
        let f_low = 100.0_f64;
        let f_high = 200.0_f64;
        let f_center = 150.0_f64;
        let n = 2000;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * f_center * i as f64 / fs).sin())
            .collect();
        let mut bp = ButterworthBandpass::new(2, f_low, f_high, fs);
        let out = bp.apply(&signal);
        let out_rms = rms(&out[200..]);
        assert!(
            out_rms > 0.2,
            "bandpass should pass center freq, rms={out_rms}"
        );
    }
    #[test]
    fn butterworth_bandpass_rejects_dc() {
        let fs = 1000.0_f64;
        let signal = vec![1.0_f64; 1000];
        let mut bp = ButterworthBandpass::new(2, 100.0, 300.0, fs);
        let out = bp.apply(&signal);
        let out_rms = rms(&out[200..]);
        assert!(out_rms < 0.3, "bandpass should reject DC, rms={out_rms}");
    }
    #[test]
    fn savgol_preserves_constant_signal() {
        let signal = vec![5.0_f64; 50];
        let out = savitzky_golay(&signal, 3, 2);
        for (i, &v) in out.iter().enumerate() {
            assert!((v - 5.0).abs() < 1e-8, "SG constant at {i}: {v}");
        }
    }
    #[test]
    fn savgol_preserves_linear_signal() {
        let signal: Vec<f64> = (0..40).map(|i| i as f64 * 2.0 + 1.0).collect();
        let out = savitzky_golay(&signal, 3, 2);
        for i in 3..37 {
            assert!(
                (out[i] - signal[i]).abs() < 1e-6,
                "SG linear at {i}: {}",
                out[i] - signal[i]
            );
        }
    }
    #[test]
    fn savgol_smooths_noisy_signal() {
        let n = 100;
        let clean: Vec<f64> = (0..n).map(|i| (i as f64 / n as f64).powi(2)).collect();
        let noisy: Vec<f64> = clean
            .iter()
            .enumerate()
            .map(|(i, &c)| c + 0.05 * (i as f64 * 1.23456).sin())
            .collect();
        let smoothed = savitzky_golay(&noisy, 4, 2);
        let err_noisy: f64 = noisy
            .iter()
            .zip(clean.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        let err_smooth: f64 = smoothed
            .iter()
            .zip(clean.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        assert!(
            err_smooth < err_noisy,
            "SG should reduce noise: {err_smooth} vs {err_noisy}"
        );
    }
    #[test]
    fn savgol_empty_input() {
        let out = savitzky_golay(&[], 3, 2);
        assert!(out.is_empty());
    }
    #[test]
    fn savgol_single_element() {
        let out = savitzky_golay(&[7.0], 3, 2);
        assert!(!out.is_empty());
        assert!((out[0] - 7.0).abs() < 1e-6);
    }
    #[test]
    fn savgol_output_length_matches_input() {
        let signal: Vec<f64> = (0..64).map(|i| i as f64).collect();
        let out = savitzky_golay(&signal, 5, 3);
        assert_eq!(out.len(), signal.len());
    }
    #[test]
    fn periodogram_psd_length() {
        let n = 64;
        let signal = vec![1.0_f64; n];
        let (psd, freqs) = periodogram_psd(&signal, 1000.0, None);
        assert_eq!(psd.len(), n / 2 + 1);
        assert_eq!(freqs.len(), n / 2 + 1);
    }
    #[test]
    fn periodogram_psd_nonnegative() {
        let signal: Vec<f64> = (0..64)
            .map(|i| (2.0 * PI * 10.0 * i as f64 / 64.0).sin())
            .collect();
        let (psd, _) = periodogram_psd(&signal, 64.0, None);
        for &p in &psd {
            assert!(p >= 0.0, "PSD must be non-negative, got {p}");
        }
    }
    #[test]
    fn periodogram_psd_dc_peak() {
        let n = 64;
        let signal = vec![2.0_f64; n];
        let (psd, _) = periodogram_psd(&signal, 1000.0, None);
        let max_ac = psd[1..].iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            psd[0] > max_ac * 10.0,
            "DC should dominate: psd[0]={} max_ac={max_ac}",
            psd[0]
        );
    }
    #[test]
    fn periodogram_psd_with_hann_window() {
        let n = 128;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 5.0 * i as f64 / n as f64).sin())
            .collect();
        let w = hann_window(n);
        let (psd, _) = periodogram_psd(&signal, n as f64, Some(&w));
        assert_eq!(psd.len(), n / 2 + 1);
        for &p in &psd {
            assert!(p.is_finite(), "PSD must be finite");
        }
    }
    #[test]
    fn periodogram_psd_empty() {
        let (psd, freqs) = periodogram_psd(&[], 1000.0, None);
        assert!(psd.is_empty());
        assert!(freqs.is_empty());
    }
    #[test]
    fn periodogram_freqs_start_at_zero() {
        let signal: Vec<f64> = (0..32).map(|_| 1.0).collect();
        let (_, freqs) = periodogram_psd(&signal, 100.0, None);
        assert!(freqs[0].abs() < 1e-12, "first freq bin should be 0 Hz");
    }
    #[test]
    fn hann_window_endpoints_near_zero() {
        let w = hann_window(16);
        assert!(w[0].abs() < 1e-10, "Hann window should start near 0");
    }
    #[test]
    fn hamming_window_length_32() {
        let w = hamming_window(32);
        assert_eq!(w.len(), 32);
    }
    #[test]
    fn blackman_window_sum_positive() {
        let w = blackman_window(64);
        let sum: f64 = w.iter().sum();
        assert!(sum > 0.0, "Blackman window sum should be positive");
    }
    #[test]
    fn blackman_window_endpoints() {
        let w = blackman_window(16);
        assert!(w[0].abs() < 1e-10, "Blackman should start near 0");
    }
    #[test]
    fn power_spectrum_nonnegative_extended() {
        let signal: Vec<f64> = (0..32).map(|i| (i as f64 * 0.5).sin()).collect();
        let ps = power_spectrum(&signal);
        for &p in &ps {
            assert!(p >= 0.0, "power spectrum must be non-negative: {p}");
        }
    }
    #[test]
    fn fft_real_length() {
        let signal = vec![1.0_f64; 16];
        let spectrum = fft_real(&signal);
        assert_eq!(spectrum.len(), 16);
    }
    #[test]
    fn butterworth_lowpass_output_length() {
        let signal: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let mut lp = ButterworthLowPass::new(2, 100.0, 1000.0);
        let out = lp.apply(&signal);
        assert_eq!(out.len(), signal.len());
    }
    #[test]
    fn butterworth_lowpass_dc_passes() {
        let fs = 1000.0_f64;
        let n = 500;
        let signal = vec![1.0_f64; n];
        let mut lp = ButterworthLowPass::new(2, 200.0, fs);
        let out = lp.apply(&signal);
        let out_mean: f64 = out[200..].iter().sum::<f64>() / (n - 200) as f64;
        assert!(
            (out_mean - 1.0).abs() < 0.05,
            "DC should pass: mean={out_mean}"
        );
    }
    #[test]
    fn butterworth_lowpass_high_freq_attenuated() {
        let fs = 1000.0_f64;
        let n = 1000;
        let f_high = 400.0_f64;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * f_high * i as f64 / fs).sin())
            .collect();
        let mut lp = ButterworthLowPass::new(4, 100.0, fs);
        let out = lp.apply(&signal);
        let out_rms = rms(&out[200..]);
        assert!(
            out_rms < 0.1,
            "high freq should be attenuated: rms={out_rms}"
        );
    }
    #[test]
    fn moving_average_constant_signal_v2() {
        let signal = vec![3.0_f64; 20];
        let out = moving_average(&signal, 5);
        for (i, &v) in out.iter().enumerate() {
            assert!((v - 3.0).abs() < 1e-10, "MA constant at {i}: {v}");
        }
    }
    #[test]
    fn moving_average_output_length() {
        let signal: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let out = moving_average(&signal, 5);
        assert_eq!(out.len(), signal.len());
    }
    #[test]
    fn rms_of_ones() {
        let s = vec![1.0_f64; 100];
        assert!((rms(&s) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn rms_of_sine() {
        let n = 1024;
        let s: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let r = rms(&s);
        assert!((r - 1.0 / 2.0_f64.sqrt()).abs() < 0.01, "rms of sine = {r}");
    }
    #[test]
    fn zero_crossings_sine_4_periods() {
        let n = 256;
        let s: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * 4.0 * i as f64 / n as f64).sin())
            .collect();
        let zc = zero_crossings(&s);
        assert!(
            (6..=10).contains(&zc),
            "expected ~8 zero crossings, got {zc}"
        );
    }
    #[test]
    fn peak_to_peak_constant_signal() {
        let s = vec![5.0_f64; 20];
        assert!(peak_to_peak(&s).abs() < 1e-12);
    }
    #[test]
    fn autocorrelation_at_zero_lag() {
        let s: Vec<f64> = (0..32).map(|i| i as f64 + 1.0).collect();
        let ac0 = autocorrelation(&s, 0);
        assert!(
            (ac0 - 1.0).abs() < 1e-12,
            "autocorrelation at lag=0 should be 1.0, got {ac0}"
        );
    }
    #[test]
    fn convolve_identity_kernel_extended() {
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let kernel = vec![1.0];
        let out = convolve(&signal, &kernel);
        assert_eq!(out.len(), signal.len());
        for (a, b) in out.iter().zip(signal.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn correlate_self_peak_at_centre() {
        let s: Vec<f64> = (0..16).map(|i| (i as f64 * 0.5).sin()).collect();
        let cc = correlate(&s, &s);
        let n = s.len();
        let centre = n - 1;
        let peak_idx = cc
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        assert_eq!(
            peak_idx, centre,
            "autocorr peak at {peak_idx}, expected {centre}"
        );
    }
    #[test]
    fn spectrogram_empty_signal_returns_empty() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let sg = sp.compute_spectrogram(&[]);
        assert!(sg.is_empty(), "empty signal should give empty spectrogram");
    }
    #[test]
    fn spectrogram_bin_count_matches_frame_len() {
        let frame_len = 64;
        let sp = SignalProcessor::new(frame_len, 32, 16000.0);
        let signal: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
        let sg = sp.compute_spectrogram(&signal);
        assert!(!sg.is_empty(), "spectrogram should not be empty");
        let expected_bins = frame_len / 2 + 1;
        for (fi, frame) in sg.iter().enumerate() {
            assert_eq!(
                frame.len(),
                expected_bins,
                "frame {fi} has {} bins, expected {expected_bins}",
                frame.len()
            );
        }
    }
    #[test]
    fn spectrogram_values_are_non_negative() {
        let sp = SignalProcessor::new(32, 16, 8000.0);
        let signal: Vec<f64> = (0..128).map(|i| (i as f64).sin()).collect();
        let sg = sp.compute_spectrogram(&signal);
        for (fi, frame) in sg.iter().enumerate() {
            for (ki, &v) in frame.iter().enumerate() {
                assert!(v >= 0.0, "spectrogram[{fi}][{ki}] = {v} < 0");
            }
        }
    }
    #[test]
    fn mel_filterbank_shape_is_correct() {
        let sp = SignalProcessor::new(512, 256, 22050.0);
        let fb = sp.compute_mel_filterbank(40, 257, 0.0, 8000.0);
        assert_eq!(fb.len(), 40, "should have 40 Mel filters");
        for (i, row) in fb.iter().enumerate() {
            assert_eq!(row.len(), 257, "filter {i} should have 257 bins");
        }
    }
    #[test]
    fn mel_filterbank_values_in_unit_range() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let fb = sp.compute_mel_filterbank(20, 129, 80.0, 7600.0);
        for (i, row) in fb.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!(
                    (0.0..=1.0 + 1e-12).contains(&v),
                    "mel_filterbank[{i}][{j}] = {v} out of [0,1]"
                );
            }
        }
    }
    #[test]
    fn mel_filterbank_empty_for_zero_mels() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let fb = sp.compute_mel_filterbank(0, 128, 0.0, 8000.0);
        assert!(fb.is_empty(), "0 mel filters should give empty filterbank");
    }
    #[test]
    fn acf_at_zero_lag_is_one() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let s: Vec<f64> = (0..64).map(|i| (i as f64 * 0.3).sin()).collect();
        let acf = sp.compute_autocorrelation(&s, 0);
        assert!(
            (acf[0] - 1.0).abs() < 1e-12,
            "ACF[0] should be 1, got {}",
            acf[0]
        );
    }
    #[test]
    fn acf_length_matches_max_lag_plus_one() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let s = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let max_lag = 3;
        let acf = sp.compute_autocorrelation(&s, max_lag);
        assert_eq!(
            acf.len(),
            max_lag + 1,
            "ACF length should be {}, got {}",
            max_lag + 1,
            acf.len()
        );
    }
    #[test]
    fn acf_values_bounded_by_one() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let s: Vec<f64> = (0..32).map(|i| (i as f64 * PI / 8.0).sin()).collect();
        let acf = sp.compute_autocorrelation(&s, 10);
        for (k, &v) in acf.iter().enumerate() {
            assert!(
                (-1.0 - 1e-10..=1.0 + 1e-10).contains(&v),
                "ACF[{k}] = {v} out of [-1, 1]"
            );
        }
    }
    #[test]
    fn acf_constant_signal_returns_ones() {
        let sp = SignalProcessor::new(256, 128, 16000.0);
        let s = vec![3.7_f64; 20];
        let acf = sp.compute_autocorrelation(&s, 5);
        for (k, &v) in acf.iter().enumerate() {
            assert!(
                (v - 1.0).abs() < 1e-12,
                "constant signal ACF[{k}] should be 1, got {v}"
            );
        }
    }
}
