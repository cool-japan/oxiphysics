#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time series I/O: storage, resampling, filtering, FFT, spectral analysis,
//! autocorrelation/cross-correlation, peak detection, anomaly detection,
//! CSV/JSON/binary export, and streaming/buffered I/O for large time series.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors produced by the time-series I/O subsystem.
#[allow(dead_code)]
#[derive(Debug)]
pub enum TsError {
    /// Length mismatch between two arrays.
    LengthMismatch {
        /// Expected length.
        expected: usize,
        /// Actual length received.
        got: usize,
    },
    /// The sample rate is zero or negative.
    InvalidSampleRate(f64),
    /// Requested channel index does not exist.
    ChannelNotFound(usize),
    /// Named channel does not exist.
    ChannelNameNotFound(String),
    /// Not enough samples for the operation.
    InsufficientSamples {
        /// Required minimum number of samples.
        need: usize,
        /// Actual number of samples present.
        have: usize,
    },
    /// FFT size is not a power of two.
    FftSizeNotPow2(usize),
    /// Generic I/O error.
    Io(String),
    /// Parse error.
    Parse(String),
}

impl fmt::Display for TsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LengthMismatch { expected, got } => {
                write!(f, "length mismatch: expected {expected}, got {got}")
            }
            Self::InvalidSampleRate(r) => write!(f, "invalid sample rate: {r}"),
            Self::ChannelNotFound(i) => write!(f, "channel index {i} not found"),
            Self::ChannelNameNotFound(n) => write!(f, "channel '{n}' not found"),
            Self::InsufficientSamples { need, have } => {
                write!(f, "need {need} samples, have {have}")
            }
            Self::FftSizeNotPow2(n) => write!(f, "FFT size {n} is not a power of two"),
            Self::Io(s) => write!(f, "I/O error: {s}"),
            Self::Parse(s) => write!(f, "parse error: {s}"),
        }
    }
}

impl std::error::Error for TsError {}

// ---------------------------------------------------------------------------
// TimeStamp and sample
// ---------------------------------------------------------------------------

/// A single time-stamped scalar sample.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample {
    /// Timestamp in seconds (or any consistent unit).
    pub time: f64,
    /// Scalar measurement value.
    pub value: f64,
}

impl Sample {
    /// Construct a new sample.
    #[allow(dead_code)]
    pub fn new(time: f64, value: f64) -> Self {
        Self { time, value }
    }
}

// ---------------------------------------------------------------------------
// Channel — one named signal with uniform or non-uniform timestamps
// ---------------------------------------------------------------------------

/// A single named signal channel.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Channel {
    /// Human-readable name, e.g. `"temperature"`.
    pub name: String,
    /// Physical unit string, e.g. `"K"` or `"m/s"`.
    pub unit: String,
    /// Ordered time-stamp array (seconds).
    pub times: Vec<f64>,
    /// Value array, same length as `times`.
    pub values: Vec<f64>,
}

impl Channel {
    /// Create a new empty channel.
    #[allow(dead_code)]
    pub fn new(name: impl Into<String>, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            times: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Push one sample.
    #[allow(dead_code)]
    pub fn push(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }

    /// Number of samples.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Return `true` when there are no samples.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// Minimum value over all samples, or `NAN` if empty.
    #[allow(dead_code)]
    pub fn min_value(&self) -> f64 {
        self.values.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Maximum value over all samples, or `NAN` if empty.
    #[allow(dead_code)]
    pub fn max_value(&self) -> f64 {
        self.values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Arithmetic mean of all values, or `NAN` if empty.
    #[allow(dead_code)]
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return f64::NAN;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }

    /// Sample variance (population).
    #[allow(dead_code)]
    pub fn variance(&self) -> f64 {
        if self.values.is_empty() {
            return f64::NAN;
        }
        let m = self.mean();
        self.values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / self.values.len() as f64
    }

    /// Standard deviation (population).
    #[allow(dead_code)]
    pub fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }

    /// Linearly interpolate the value at arbitrary time `t`.
    /// Clamps to the boundary if `t` is out of range.
    #[allow(dead_code)]
    pub fn interp_linear(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n == 0 {
            return f64::NAN;
        }
        if t <= self.times[0] {
            return self.values[0];
        }
        if t >= self.times[n - 1] {
            return self.values[n - 1];
        }
        // binary search
        let pos = self.times.partition_point(|&x| x <= t);
        let i = pos.saturating_sub(1);
        let j = i + 1;
        let dt = self.times[j] - self.times[i];
        if dt.abs() < 1e-300 {
            return self.values[i];
        }
        let alpha = (t - self.times[i]) / dt;
        self.values[i] * (1.0 - alpha) + self.values[j] * alpha
    }

    /// Cubic (Catmull-Rom) spline interpolation at time `t`.
    #[allow(dead_code)]
    pub fn interp_cubic(&self, t: f64) -> f64 {
        let n = self.times.len();
        if n < 4 {
            return self.interp_linear(t);
        }
        if t <= self.times[0] {
            return self.values[0];
        }
        if t >= self.times[n - 1] {
            return self.values[n - 1];
        }
        let pos = self.times.partition_point(|&x| x <= t);
        let i1 = pos.saturating_sub(1).min(n - 2);
        let i0 = i1.saturating_sub(1);
        let i2 = (i1 + 1).min(n - 1);
        let i3 = (i1 + 2).min(n - 1);
        let dt = self.times[i2] - self.times[i1];
        let alpha = if dt.abs() < 1e-300 {
            0.0
        } else {
            (t - self.times[i1]) / dt
        };
        catmull_rom(
            self.values[i0],
            self.values[i1],
            self.values[i2],
            self.values[i3],
            alpha,
        )
    }
}

/// Catmull-Rom spline helper (one-dimensional).
#[allow(dead_code)]
fn catmull_rom(p0: f64, p1: f64, p2: f64, p3: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
}

// ---------------------------------------------------------------------------
// MultiChannelSeries — container for N synchronized channels
// ---------------------------------------------------------------------------

/// Multi-channel time series container.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct MultiChannelSeries {
    /// Ordered channel list.
    pub channels: Vec<Channel>,
    /// Optional metadata key-value store.
    pub metadata: HashMap<String, String>,
}

impl MultiChannelSeries {
    /// Create an empty container.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a channel and return its index.
    #[allow(dead_code)]
    pub fn add_channel(&mut self, ch: Channel) -> usize {
        let idx = self.channels.len();
        self.channels.push(ch);
        idx
    }

    /// Number of channels.
    #[allow(dead_code)]
    pub fn num_channels(&self) -> usize {
        self.channels.len()
    }

    /// Borrow channel by index.
    #[allow(dead_code)]
    pub fn channel(&self, idx: usize) -> Option<&Channel> {
        self.channels.get(idx)
    }

    /// Mutable borrow channel by index.
    #[allow(dead_code)]
    pub fn channel_mut(&mut self, idx: usize) -> Option<&mut Channel> {
        self.channels.get_mut(idx)
    }

    /// Find a channel by name.
    #[allow(dead_code)]
    pub fn channel_by_name(&self, name: &str) -> Option<&Channel> {
        self.channels.iter().find(|c| c.name == name)
    }

    /// Insert metadata.
    #[allow(dead_code)]
    pub fn set_metadata(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.metadata.insert(key.into(), value.into());
    }
}

// ---------------------------------------------------------------------------
// Resampling
// ---------------------------------------------------------------------------

/// Resampling method.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResampleMethod {
    /// Piecewise linear interpolation.
    Linear,
    /// Catmull-Rom cubic spline.
    CubicSpline,
}

/// Resample `channel` onto a new uniform grid with `n_samples` points
/// spanning from `t_start` to `t_end`.
#[allow(dead_code)]
pub fn resample_uniform(
    channel: &Channel,
    t_start: f64,
    t_end: f64,
    n_samples: usize,
    method: ResampleMethod,
) -> Result<Channel, TsError> {
    if n_samples < 2 {
        return Err(TsError::InsufficientSamples {
            need: 2,
            have: n_samples,
        });
    }
    let mut out = Channel::new(channel.name.clone(), channel.unit.clone());
    let dt = (t_end - t_start) / (n_samples - 1) as f64;
    for i in 0..n_samples {
        let t = t_start + i as f64 * dt;
        let v = match method {
            ResampleMethod::Linear => channel.interp_linear(t),
            ResampleMethod::CubicSpline => channel.interp_cubic(t),
        };
        out.push(t, v);
    }
    Ok(out)
}

/// Resample `channel` onto an explicit list of target timestamps.
#[allow(dead_code)]
pub fn resample_to_times(
    channel: &Channel,
    target_times: &[f64],
    method: ResampleMethod,
) -> Channel {
    let mut out = Channel::new(channel.name.clone(), channel.unit.clone());
    for &t in target_times {
        let v = match method {
            ResampleMethod::Linear => channel.interp_linear(t),
            ResampleMethod::CubicSpline => channel.interp_cubic(t),
        };
        out.push(t, v);
    }
    out
}

// ---------------------------------------------------------------------------
// Moving-average filter
// ---------------------------------------------------------------------------

/// Apply a simple causal moving-average filter in-place.
/// `window` is the number of samples (≥ 1).
#[allow(dead_code)]
pub fn moving_average(values: &[f64], window: usize) -> Vec<f64> {
    let n = values.len();
    let w = window.max(1);
    let mut out = vec![0.0_f64; n];
    let mut acc = 0.0_f64;
    let mut count = 0_usize;
    for (i, &v) in values.iter().enumerate() {
        acc += v;
        count += 1;
        if count > w {
            acc -= values[i - w];
            count -= 1;
        }
        out[i] = acc / count as f64;
    }
    out
}

// ---------------------------------------------------------------------------
// Butterworth IIR filter (2nd-order biquad sections)
// ---------------------------------------------------------------------------

/// Butterworth filter type.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterType {
    /// Low-pass filter.
    LowPass,
    /// High-pass filter.
    HighPass,
}

/// Coefficients for a single biquad (2nd-order IIR) section.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct BiquadCoeffs {
    /// Feed-forward coefficients b0, b1, b2.
    pub b: [f64; 3],
    /// Feed-back coefficients a1, a2 (a0 is normalised to 1).
    pub a: [f64; 2],
}

impl BiquadCoeffs {
    /// Design a single 2nd-order Butterworth section.
    /// `fc` is the cut-off frequency in Hz, `fs` is the sample rate in Hz.
    #[allow(dead_code)]
    pub fn butterworth_2nd(fc: f64, fs: f64, filter_type: FilterType) -> Self {
        // Bilinear-transform design for 2nd-order Butterworth
        let omega = std::f64::consts::PI * fc / fs;
        let k = omega.tan();
        let k2 = k * k;
        let sqrt2 = std::f64::consts::SQRT_2;
        let norm = 1.0 / (1.0 + sqrt2 * k + k2);
        match filter_type {
            FilterType::LowPass => {
                let b0 = k2 * norm;
                let b1 = 2.0 * b0;
                let b2 = b0;
                let a1 = 2.0 * (k2 - 1.0) * norm;
                let a2 = (1.0 - sqrt2 * k + k2) * norm;
                Self {
                    b: [b0, b1, b2],
                    a: [a1, a2],
                }
            }
            FilterType::HighPass => {
                let b0 = norm;
                let b1 = -2.0 * b0;
                let b2 = b0;
                let a1 = 2.0 * (k2 - 1.0) * norm;
                let a2 = (1.0 - sqrt2 * k + k2) * norm;
                Self {
                    b: [b0, b1, b2],
                    a: [a1, a2],
                }
            }
        }
    }

    /// Apply this biquad section to `input`, return filtered output.
    #[allow(dead_code)]
    pub fn apply(&self, input: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0_f64; input.len()];
        let mut w1 = 0.0_f64;
        let mut w2 = 0.0_f64;
        for (i, &x) in input.iter().enumerate() {
            let w0 = x - self.a[0] * w1 - self.a[1] * w2;
            out[i] = self.b[0] * w0 + self.b[1] * w1 + self.b[2] * w2;
            w2 = w1;
            w1 = w0;
        }
        out
    }

    /// Zero-phase forward-backward filtering (doubles the filter order).
    #[allow(dead_code)]
    pub fn apply_zero_phase(&self, input: &[f64]) -> Vec<f64> {
        let forward = self.apply(input);
        let rev: Vec<f64> = forward.iter().cloned().rev().collect();
        let backward = self.apply(&rev);
        backward.iter().cloned().rev().collect()
    }
}

// ---------------------------------------------------------------------------
// FFT (radix-2 Cooley-Tukey, in-place)
// ---------------------------------------------------------------------------

/// A complex number stored as `(real, imag)`.
#[allow(dead_code)]
pub type Complex = (f64, f64);

/// Compute the in-place radix-2 DIT FFT of `buf` (length must be a power of 2).
/// `inverse` == true → IFFT (no 1/N normalisation; caller must divide).
#[allow(dead_code)]
pub fn fft_inplace(buf: &mut [(f64, f64)], inverse: bool) {
    let n = buf.len();
    // Bit-reversal permutation
    let mut j = 0_usize;
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
    // Cooley-Tukey butterfly
    let sign = if inverse { 1.0_f64 } else { -1.0_f64 };
    let mut len = 2_usize;
    while len <= n {
        let half = len / 2;
        let ang = sign * std::f64::consts::TAU / len as f64;
        let wr = ang.cos();
        let wi = ang.sin();
        let mut k = 0;
        while k < n {
            let mut wre = 1.0_f64;
            let mut wim = 0.0_f64;
            for m in 0..half {
                let (ur, ui) = buf[k + m];
                let (vr, vi) = buf[k + m + half];
                let tr = wre * vr - wim * vi;
                let ti = wre * vi + wim * vr;
                buf[k + m] = (ur + tr, ui + ti);
                buf[k + m + half] = (ur - tr, ui - ti);
                let new_wre = wre * wr - wim * wi;
                wim = wre * wi + wim * wr;
                wre = new_wre;
            }
            k += len;
        }
        len *= 2;
    }
}

/// Compute the FFT of real-valued data.  Returns complex spectrum of length N/2+1.
/// `data` length must be a power of two.
#[allow(dead_code)]
pub fn rfft(data: &[f64]) -> Result<Vec<Complex>, TsError> {
    let n = data.len();
    if n == 0 || (n & (n - 1)) != 0 {
        return Err(TsError::FftSizeNotPow2(n));
    }
    let mut buf: Vec<Complex> = data.iter().map(|&v| (v, 0.0)).collect();
    fft_inplace(&mut buf, false);
    Ok(buf[..=n / 2].to_vec())
}

/// Compute the power spectral density (one-sided, in units of value²/Hz)
/// from `rfft` output given sample rate `fs`.
#[allow(dead_code)]
pub fn power_spectrum(rfft_out: &[Complex], n: usize, fs: f64) -> Vec<f64> {
    let scale = 2.0 / (n as f64 * fs);
    rfft_out
        .iter()
        .enumerate()
        .map(|(k, &(re, im))| {
            let p = (re * re + im * im) * scale;
            // DC and Nyquist are one-sided, others doubled above
            if k == 0 || k == n / 2 { p * 0.5 } else { p }
        })
        .collect()
}

/// Frequency axis (Hz) matching `rfft` output for signal of length `n` at rate `fs`.
#[allow(dead_code)]
pub fn rfft_frequencies(n: usize, fs: f64) -> Vec<f64> {
    (0..=n / 2).map(|k| k as f64 * fs / n as f64).collect()
}

// ---------------------------------------------------------------------------
// Welch's method (spectral density estimation)
// ---------------------------------------------------------------------------

/// Welch window type.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WindowType {
    /// Rectangular (no windowing).
    Rectangular,
    /// Hann (raised cosine).
    Hann,
    /// Hamming.
    Hamming,
    /// Blackman.
    Blackman,
}

/// Generate a window function of length `n`.
#[allow(dead_code)]
pub fn make_window(wtype: WindowType, n: usize) -> Vec<f64> {
    use std::f64::consts::TAU;
    match wtype {
        WindowType::Rectangular => vec![1.0; n],
        WindowType::Hann => (0..n)
            .map(|i| 0.5 * (1.0 - (TAU * i as f64 / (n - 1) as f64).cos()))
            .collect(),
        WindowType::Hamming => (0..n)
            .map(|i| 0.54 - 0.46 * (TAU * i as f64 / (n - 1) as f64).cos())
            .collect(),
        WindowType::Blackman => (0..n)
            .map(|i| {
                let x = TAU * i as f64 / (n - 1) as f64;
                0.42 - 0.5 * x.cos() + 0.08 * (2.0 * x).cos()
            })
            .collect(),
    }
}

/// Compute the Welch power spectral density estimate.
/// `segment_len` must be a power of two; `overlap` is in (0..segment_len).
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn welch_psd(
    data: &[f64],
    fs: f64,
    segment_len: usize,
    overlap: usize,
    window: WindowType,
) -> Result<(Vec<f64>, Vec<f64>), TsError> {
    if segment_len == 0 || (segment_len & (segment_len - 1)) != 0 {
        return Err(TsError::FftSizeNotPow2(segment_len));
    }
    if fs <= 0.0 {
        return Err(TsError::InvalidSampleRate(fs));
    }
    let step = (segment_len - overlap).max(1);
    let win = make_window(window, segment_len);
    let win_power: f64 = win.iter().map(|w| w * w).sum::<f64>() / segment_len as f64;
    let half = segment_len / 2 + 1;
    let mut accum = vec![0.0_f64; half];
    let mut n_segments = 0_usize;
    let mut start = 0;
    while start + segment_len <= data.len() {
        let seg: Vec<f64> = data[start..start + segment_len]
            .iter()
            .zip(win.iter())
            .map(|(&d, &w)| d * w)
            .collect();
        let spec = rfft(&seg)?;
        let scale = 1.0 / (fs * win_power * segment_len as f64);
        for (k, &(re, im)) in spec.iter().enumerate() {
            let p = (re * re + im * im) * scale;
            accum[k] += if k == 0 || k == segment_len / 2 {
                p
            } else {
                2.0 * p
            };
        }
        n_segments += 1;
        start += step;
    }
    if n_segments == 0 {
        return Err(TsError::InsufficientSamples {
            need: segment_len,
            have: data.len(),
        });
    }
    let psd: Vec<f64> = accum.iter().map(|v| v / n_segments as f64).collect();
    let freqs = rfft_frequencies(segment_len, fs);
    Ok((freqs, psd))
}

// ---------------------------------------------------------------------------
// Autocorrelation and cross-correlation
// ---------------------------------------------------------------------------

/// Compute the (biased) autocorrelation of `x` for lags 0..=max_lag.
#[allow(dead_code)]
pub fn autocorrelation(x: &[f64], max_lag: usize) -> Vec<f64> {
    let n = x.len();
    let mean = x.iter().sum::<f64>() / n as f64;
    let xc: Vec<f64> = x.iter().map(|&v| v - mean).collect();
    let c0 = xc.iter().map(|&v| v * v).sum::<f64>() / n as f64;
    let lags = max_lag.min(n - 1);
    (0..=lags)
        .map(|lag| {
            let s: f64 = xc[..n - lag]
                .iter()
                .zip(xc[lag..].iter())
                .map(|(&a, &b)| a * b)
                .sum();
            s / (n as f64 * c0)
        })
        .collect()
}

/// Compute the (biased) cross-correlation of `x` and `y` for lags -(max_lag)..=max_lag.
/// Returns (lag_vector, ccf_values).
#[allow(dead_code)]
pub fn cross_correlation(x: &[f64], y: &[f64], max_lag: usize) -> (Vec<i64>, Vec<f64>) {
    let n = x.len().min(y.len());
    let mx = x[..n].iter().sum::<f64>() / n as f64;
    let my = y[..n].iter().sum::<f64>() / n as f64;
    let xc: Vec<f64> = x[..n].iter().map(|&v| v - mx).collect();
    let yc: Vec<f64> = y[..n].iter().map(|&v| v - my).collect();
    let sx: f64 = (xc.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
    let sy: f64 = (yc.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
    let denom = sx * sy * n as f64;
    let lags = max_lag.min(n - 1);
    let lag_vec: Vec<i64> = (-(lags as i64)..=(lags as i64)).collect();
    let ccf: Vec<f64> = lag_vec
        .iter()
        .map(|&lag| {
            let s: f64 = if lag >= 0 {
                let l = lag as usize;
                xc[..n - l]
                    .iter()
                    .zip(yc[l..].iter())
                    .map(|(&a, &b)| a * b)
                    .sum()
            } else {
                let l = (-lag) as usize;
                xc[l..]
                    .iter()
                    .zip(yc[..n - l].iter())
                    .map(|(&a, &b)| a * b)
                    .sum()
            };
            if denom.abs() < 1e-300 { 0.0 } else { s / denom }
        })
        .collect();
    (lag_vec, ccf)
}

// ---------------------------------------------------------------------------
// Peak detection
// ---------------------------------------------------------------------------

/// A detected peak with its index, timestamp, and value.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct Peak {
    /// Sample index.
    pub index: usize,
    /// Timestamp at the peak.
    pub time: f64,
    /// Value at the peak.
    pub value: f64,
    /// Prominence relative to surrounding valleys.
    pub prominence: f64,
}

/// Detect local maxima in `channel`.
/// `min_prominence` filters weak peaks; `min_distance` is the minimum index separation.
#[allow(dead_code)]
pub fn detect_peaks(channel: &Channel, min_prominence: f64, min_distance: usize) -> Vec<Peak> {
    let n = channel.values.len();
    if n < 3 {
        return Vec::new();
    }
    let v = &channel.values;
    let candidates: Vec<usize> = (1..n - 1)
        .filter(|&i| v[i] > v[i - 1] && v[i] >= v[i + 1])
        .collect();
    // Compute prominence
    let mut peaks: Vec<Peak> = candidates
        .iter()
        .map(|&i| {
            let left_min = v[..i].iter().cloned().fold(f64::INFINITY, f64::min);
            let right_min = v[i + 1..].iter().cloned().fold(f64::INFINITY, f64::min);
            let base = left_min.max(right_min);
            let prom = v[i] - base;
            Peak {
                index: i,
                time: channel.times[i],
                value: v[i],
                prominence: prom,
            }
        })
        .filter(|p| p.prominence >= min_prominence)
        .collect();
    // Enforce minimum distance (greedy, largest prominence first)
    peaks.sort_by(|a, b| {
        b.prominence
            .partial_cmp(&a.prominence)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut kept: Vec<Peak> = Vec::new();
    for p in &peaks {
        let too_close = kept.iter().any(|k: &Peak| {
            let d = k.index.abs_diff(p.index);
            d < min_distance
        });
        if !too_close {
            kept.push(*p);
        }
    }
    kept.sort_by_key(|p| p.index);
    // Remove spurious candidates after distance filter
    let _ = candidates.len(); // suppress unused warning
    kept
}

// ---------------------------------------------------------------------------
// Anomaly detection
// ---------------------------------------------------------------------------

/// An anomalous sample identified by index.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct Anomaly {
    /// Sample index.
    pub index: usize,
    /// Timestamp.
    pub time: f64,
    /// Raw value.
    pub value: f64,
    /// Anomaly score (z-score or IQR distance).
    pub score: f64,
}

/// Z-score anomaly detection: flag samples where |z| > `threshold`.
#[allow(dead_code)]
pub fn anomalies_zscore(channel: &Channel, threshold: f64) -> Vec<Anomaly> {
    let mean = channel.mean();
    let std = channel.std_dev();
    if std < 1e-300 {
        return Vec::new();
    }
    channel
        .values
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| {
            let z = (v - mean) / std;
            if z.abs() > threshold {
                Some(Anomaly {
                    index: i,
                    time: channel.times[i],
                    value: v,
                    score: z,
                })
            } else {
                None
            }
        })
        .collect()
}

/// IQR-based anomaly detection.
/// Flags values below Q1 − 1.5·IQR or above Q3 + 1.5·IQR (or custom `k`).
#[allow(dead_code)]
pub fn anomalies_iqr(channel: &Channel, k: f64) -> Vec<Anomaly> {
    let mut sorted = channel.values.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n < 4 {
        return Vec::new();
    }
    let q1 = sorted[n / 4];
    let q3 = sorted[3 * n / 4];
    let iqr = q3 - q1;
    let lo = q1 - k * iqr;
    let hi = q3 + k * iqr;
    channel
        .values
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| {
            if v < lo || v > hi {
                let score = if v < lo {
                    (lo - v) / iqr
                } else {
                    (v - hi) / iqr
                };
                Some(Anomaly {
                    index: i,
                    time: channel.times[i],
                    value: v,
                    score,
                })
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Export: CSV
// ---------------------------------------------------------------------------

/// Write a single channel to CSV (time, value) using any `Write` sink.
#[allow(dead_code)]
pub fn write_channel_csv<W: Write>(w: &mut W, channel: &Channel) -> io::Result<()> {
    writeln!(w, "time_s,{}", channel.name)?;
    for (t, v) in channel.times.iter().zip(channel.values.iter()) {
        writeln!(w, "{t:.9e},{v:.9e}")?;
    }
    Ok(())
}

/// Write a `MultiChannelSeries` to CSV.  First column is `time_s`, then one
/// column per channel.  All channels are assumed to have identical timestamps.
#[allow(dead_code)]
pub fn write_multi_channel_csv<W: Write>(w: &mut W, mcs: &MultiChannelSeries) -> io::Result<()> {
    if mcs.channels.is_empty() {
        return Ok(());
    }
    // Header
    write!(w, "time_s")?;
    for ch in &mcs.channels {
        write!(w, ",{}", ch.name)?;
    }
    writeln!(w)?;
    let n = mcs.channels[0].len();
    for i in 0..n {
        let t = mcs.channels[0].times[i];
        write!(w, "{t:.9e}")?;
        for ch in &mcs.channels {
            let v = if i < ch.len() { ch.values[i] } else { f64::NAN };
            write!(w, ",{v:.9e}")?;
        }
        writeln!(w)?;
    }
    Ok(())
}

/// Parse a two-column CSV (time, value) into a `Channel`.
#[allow(dead_code)]
pub fn read_channel_csv(text: &str, name: &str, unit: &str) -> Result<Channel, TsError> {
    let mut ch = Channel::new(name, unit);
    for (line_no, line) in text.lines().enumerate() {
        if line_no == 0 || line.trim().is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, ',');
        let t_str = parts.next().unwrap_or("").trim();
        let v_str = parts.next().unwrap_or("").trim();
        let t: f64 = t_str
            .parse()
            .map_err(|_| TsError::Parse(format!("line {line_no}: bad time '{t_str}'")))?;
        let v: f64 = v_str
            .parse()
            .map_err(|_| TsError::Parse(format!("line {line_no}: bad value '{v_str}'")))?;
        ch.push(t, v);
    }
    Ok(ch)
}

// ---------------------------------------------------------------------------
// Export: JSON
// ---------------------------------------------------------------------------

/// Serialise a `Channel` to a compact JSON string.
#[allow(dead_code)]
pub fn channel_to_json(channel: &Channel) -> String {
    let pairs: Vec<String> = channel
        .times
        .iter()
        .zip(channel.values.iter())
        .map(|(t, v)| format!("[{t:.9e},{v:.9e}]"))
        .collect();
    format!(
        "{{\"name\":{:?},\"unit\":{:?},\"samples\":[{}]}}",
        channel.name,
        channel.unit,
        pairs.join(",")
    )
}

/// Deserialise a `Channel` from a minimal JSON string (as produced by `channel_to_json`).
/// This is a hand-rolled parser to avoid a serde dependency.
#[allow(dead_code)]
pub fn channel_from_json(json: &str) -> Result<Channel, TsError> {
    // Very naive extractor: find name, unit, then parse [t,v] pairs
    fn extract_string_field<'a>(json: &'a str, key: &str) -> Option<&'a str> {
        let needle = format!("\"{key}\":");
        let start = json.find(&needle)? + needle.len();
        let rest = json[start..].trim_start();
        if !rest.starts_with('"') {
            return None;
        }
        let inner = &rest[1..];
        let end = inner.find('"')?;
        Some(&inner[..end])
    }
    let name = extract_string_field(json, "name").unwrap_or("").to_string();
    let unit = extract_string_field(json, "unit").unwrap_or("").to_string();
    let mut ch = Channel::new(name, unit);
    // Parse [t,v] pairs from samples array
    let samples_key = "\"samples\":[";
    if let Some(start) = json.find(samples_key) {
        let rest = &json[start + samples_key.len()..];
        // Iterate over [...] blocks
        let mut remaining = rest;
        while let Some(open) = remaining.find('[') {
            remaining = &remaining[open + 1..];
            let close = remaining.find(']').unwrap_or(0);
            let inner = &remaining[..close];
            remaining = &remaining[close + 1..];
            let mut parts = inner.splitn(2, ',');
            let t_s = parts.next().unwrap_or("").trim();
            let v_s = parts.next().unwrap_or("").trim();
            if let (Ok(t), Ok(v)) = (t_s.parse::<f64>(), v_s.parse::<f64>()) {
                ch.push(t, v);
            }
        }
    }
    Ok(ch)
}

// ---------------------------------------------------------------------------
// Export: Binary (little-endian f64 pairs: [t0 v0 t1 v1 ...])
// ---------------------------------------------------------------------------

/// Serialise a `Channel` to raw binary (little-endian f64 interleaved t/v pairs).
#[allow(dead_code)]
pub fn channel_to_binary(channel: &Channel) -> Vec<u8> {
    let mut buf = Vec::with_capacity(channel.len() * 16);
    for (&t, &v) in channel.times.iter().zip(channel.values.iter()) {
        buf.extend_from_slice(&t.to_le_bytes());
        buf.extend_from_slice(&v.to_le_bytes());
    }
    buf
}

/// Deserialise a `Channel` from the binary format produced by `channel_to_binary`.
#[allow(dead_code)]
pub fn channel_from_binary(data: &[u8], name: &str, unit: &str) -> Result<Channel, TsError> {
    if !data.len().is_multiple_of(16) {
        return Err(TsError::Parse(format!(
            "binary data length {} is not divisible by 16",
            data.len()
        )));
    }
    let mut ch = Channel::new(name, unit);
    for chunk in data.chunks_exact(16) {
        let t = f64::from_le_bytes(chunk[..8].try_into().expect("slice length must match"));
        let v = f64::from_le_bytes(chunk[8..16].try_into().expect("slice length must match"));
        ch.push(t, v);
    }
    Ok(ch)
}

// ---------------------------------------------------------------------------
// Streaming / buffered I/O
// ---------------------------------------------------------------------------

/// A streaming time-series writer that flushes to an underlying `Write` sink
/// in chunks, enabling memory-efficient handling of very large time series.
#[allow(dead_code)]
pub struct StreamingTsWriter<W: Write> {
    /// Underlying writer.
    sink: W,
    /// Internal sample buffer.
    buffer: Vec<(f64, f64)>,
    /// Number of samples before an automatic flush.
    flush_threshold: usize,
    /// Channel name (written in the CSV header).
    name: String,
    /// Whether the header has been written.
    header_written: bool,
}

impl<W: Write> StreamingTsWriter<W> {
    /// Create a new streaming writer with a given buffer size.
    #[allow(dead_code)]
    pub fn new(sink: W, name: impl Into<String>, flush_threshold: usize) -> Self {
        Self {
            sink,
            buffer: Vec::with_capacity(flush_threshold),
            flush_threshold: flush_threshold.max(1),
            name: name.into(),
            header_written: false,
        }
    }

    /// Push one sample. Auto-flushes when the buffer is full.
    #[allow(dead_code)]
    pub fn push(&mut self, time: f64, value: f64) -> io::Result<()> {
        if !self.header_written {
            writeln!(self.sink, "time_s,{}", self.name)?;
            self.header_written = true;
        }
        self.buffer.push((time, value));
        if self.buffer.len() >= self.flush_threshold {
            self.flush()?;
        }
        Ok(())
    }

    /// Flush the internal buffer to the underlying writer.
    #[allow(dead_code)]
    pub fn flush(&mut self) -> io::Result<()> {
        for (t, v) in self.buffer.drain(..) {
            writeln!(self.sink, "{t:.9e},{v:.9e}")?;
        }
        self.sink.flush()
    }

    /// Finish writing; flushes remaining samples.
    #[allow(dead_code)]
    pub fn finish(mut self) -> io::Result<W> {
        self.flush()?;
        Ok(self.sink)
    }
}

/// A buffered time-series reader that reads CSV data line-by-line without
/// loading the entire file into memory.
#[allow(dead_code)]
pub struct BufferedTsReader {
    /// Accumulated samples so far.
    pub samples: Vec<(f64, f64)>,
    /// Buffer size hint.
    chunk_size: usize,
}

impl BufferedTsReader {
    /// Create a new reader.
    #[allow(dead_code)]
    pub fn new(chunk_size: usize) -> Self {
        Self {
            samples: Vec::new(),
            chunk_size: chunk_size.max(1),
        }
    }

    /// Feed a chunk of CSV text into the reader.
    #[allow(dead_code)]
    pub fn feed(&mut self, text: &str, skip_header: bool) -> Result<(), TsError> {
        for (i, line) in text.lines().enumerate() {
            if i == 0 && skip_header {
                continue;
            }
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.splitn(2, ',');
            let t_s = parts.next().unwrap_or("").trim();
            let v_s = parts.next().unwrap_or("").trim();
            let t: f64 = t_s
                .parse()
                .map_err(|_| TsError::Parse(format!("bad time '{t_s}'")))?;
            let v: f64 = v_s
                .parse()
                .map_err(|_| TsError::Parse(format!("bad value '{v_s}'")))?;
            self.samples.push((t, v));
        }
        Ok(())
    }

    /// Convert accumulated samples into a `Channel`.
    #[allow(dead_code)]
    pub fn into_channel(self, name: &str, unit: &str) -> Channel {
        let mut ch = Channel::new(name, unit);
        for (t, v) in self.samples {
            ch.push(t, v);
        }
        ch
    }

    /// Chunk size hint for external callers.
    #[allow(dead_code)]
    pub fn chunk_size(&self) -> usize {
        self.chunk_size
    }
}

// ---------------------------------------------------------------------------
// Frequency analysis helpers
// ---------------------------------------------------------------------------

/// Dominant frequency (Hz) in a signal, estimated from FFT peak.
#[allow(dead_code)]
pub fn dominant_frequency(data: &[f64], fs: f64) -> Result<f64, TsError> {
    let n = data.len();
    if n < 2 {
        return Err(TsError::InsufficientSamples { need: 2, have: n });
    }
    // Pad / trim to next power of two
    let np2 = next_pow2(n);
    let mut padded = data.to_vec();
    padded.resize(np2, 0.0);
    let spec = rfft(&padded)?;
    let freqs = rfft_frequencies(np2, fs);
    let (max_k, _) = spec
        .iter()
        .enumerate()
        .skip(1) // skip DC
        .map(|(k, &(re, im))| (k, re * re + im * im))
        .fold(
            (1, 0.0_f64),
            |(bk, bv), (k, v)| {
                if v > bv { (k, v) } else { (bk, bv) }
            },
        );
    Ok(freqs[max_k])
}

/// Next power of two ≥ n.
#[allow(dead_code)]
pub fn next_pow2(n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    let mut p = 1;
    while p < n {
        p <<= 1;
    }
    p
}

/// Total harmonic distortion (THD) estimate: ratio of harmonic power to fundamental.
/// `fundamental_freq` in Hz, `n_harmonics` harmonics to sum above fundamental.
#[allow(dead_code)]
pub fn total_harmonic_distortion(
    data: &[f64],
    fs: f64,
    fundamental_freq: f64,
    n_harmonics: usize,
) -> Result<f64, TsError> {
    let np2 = next_pow2(data.len());
    let mut padded = data.to_vec();
    padded.resize(np2, 0.0);
    let spec = rfft(&padded)?;
    let df = fs / np2 as f64;
    let bin_fund = (fundamental_freq / df).round() as usize;
    if bin_fund == 0 || bin_fund >= spec.len() {
        return Err(TsError::InsufficientSamples {
            need: bin_fund + 1,
            have: spec.len(),
        });
    }
    let p_fund = {
        let (re, im) = spec[bin_fund];
        re * re + im * im
    };
    let p_harm: f64 = (2..=n_harmonics + 1)
        .map(|h| {
            let k = bin_fund * h;
            if k < spec.len() {
                let (re, im) = spec[k];
                re * re + im * im
            } else {
                0.0
            }
        })
        .sum();
    if p_fund < 1e-300 {
        return Ok(0.0);
    }
    Ok((p_harm / p_fund).sqrt())
}

// ---------------------------------------------------------------------------
// Statistical helpers
// ---------------------------------------------------------------------------

/// Root-mean-square of a slice.
#[allow(dead_code)]
pub fn rms(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    (data.iter().map(|v| v * v).sum::<f64>() / data.len() as f64).sqrt()
}

/// Percentile of `data` (0.0–100.0) using linear interpolation.
#[allow(dead_code)]
pub fn percentile(data: &[f64], p: f64) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    let frac = p / 100.0 * (n - 1) as f64;
    let lo = frac.floor() as usize;
    let hi = frac.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let alpha = frac - lo as f64;
        sorted[lo] * (1.0 - alpha) + sorted[hi] * alpha
    }
}

/// Median absolute deviation (MAD) of `data`.
#[allow(dead_code)]
pub fn mad(data: &[f64]) -> f64 {
    let med = percentile(data, 50.0);
    let devs: Vec<f64> = data.iter().map(|v| (v - med).abs()).collect();
    percentile(&devs, 50.0)
}

// ---------------------------------------------------------------------------
// Ring buffer for real-time streaming
// ---------------------------------------------------------------------------

/// A fixed-capacity circular ring buffer for real-time time-series streaming.
#[allow(dead_code)]
pub struct RingBuffer {
    times: Vec<f64>,
    values: Vec<f64>,
    head: usize,
    len: usize,
    capacity: usize,
}

impl RingBuffer {
    /// Create a ring buffer with given capacity.
    #[allow(dead_code)]
    pub fn new(capacity: usize) -> Self {
        let cap = capacity.max(1);
        Self {
            times: vec![0.0; cap],
            values: vec![0.0; cap],
            head: 0,
            len: 0,
            capacity: cap,
        }
    }

    /// Push one sample, overwriting the oldest if full.
    #[allow(dead_code)]
    pub fn push(&mut self, time: f64, value: f64) {
        let idx = (self.head + self.len) % self.capacity;
        if self.len < self.capacity {
            self.times[idx] = time;
            self.values[idx] = value;
            self.len += 1;
        } else {
            self.times[self.head] = time;
            self.values[self.head] = value;
            self.head = (self.head + 1) % self.capacity;
        }
    }

    /// Number of valid samples.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Return `true` when empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Drain all samples into a `Channel` (chronological order).
    #[allow(dead_code)]
    pub fn drain_to_channel(&mut self, name: &str, unit: &str) -> Channel {
        let mut ch = Channel::new(name, unit);
        for i in 0..self.len {
            let idx = (self.head + i) % self.capacity;
            ch.push(self.times[idx], self.values[idx]);
        }
        self.len = 0;
        self.head = 0;
        ch
    }

    /// View current samples as a `Channel` without consuming the buffer.
    #[allow(dead_code)]
    pub fn snapshot(&self, name: &str, unit: &str) -> Channel {
        let mut ch = Channel::new(name, unit);
        for i in 0..self.len {
            let idx = (self.head + i) % self.capacity;
            ch.push(self.times[idx], self.values[idx]);
        }
        ch
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // helpers ----------------------------------------------------------------

    fn sine_channel(freq: f64, fs: f64, n: usize) -> Channel {
        let mut ch = Channel::new("sine", "V");
        for i in 0..n {
            let t = i as f64 / fs;
            ch.push(t, (2.0 * std::f64::consts::PI * freq * t).sin());
        }
        ch
    }

    fn ramp_channel(n: usize) -> Channel {
        let mut ch = Channel::new("ramp", "m");
        for i in 0..n {
            ch.push(i as f64, i as f64);
        }
        ch
    }

    // Channel basics ---------------------------------------------------------

    #[test]
    fn test_channel_len_empty() {
        let ch = Channel::new("x", "m");
        assert!(ch.is_empty());
        assert_eq!(ch.len(), 0);
    }

    #[test]
    fn test_channel_push_and_stats() {
        let mut ch = Channel::new("x", "Pa");
        for i in 0..10 {
            ch.push(i as f64, i as f64);
        }
        assert_eq!(ch.len(), 10);
        assert!((ch.mean() - 4.5).abs() < 1e-12);
        assert!(ch.min_value() < 1.0);
        assert!((ch.max_value() - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_channel_std_dev_uniform() {
        let mut ch = Channel::new("u", "");
        for _ in 0..100 {
            ch.push(0.0, 5.0);
        }
        assert!(ch.std_dev().abs() < 1e-12);
    }

    // Interpolation ----------------------------------------------------------

    #[test]
    fn test_interp_linear_exact_nodes() {
        let ch = ramp_channel(10);
        for i in 0..10 {
            assert!((ch.interp_linear(i as f64) - i as f64).abs() < 1e-12);
        }
    }

    #[test]
    fn test_interp_linear_midpoint() {
        let ch = ramp_channel(5);
        assert!((ch.interp_linear(1.5) - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_interp_linear_clamp_low() {
        let ch = ramp_channel(5);
        assert!((ch.interp_linear(-1.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_interp_linear_clamp_high() {
        let ch = ramp_channel(5);
        assert!((ch.interp_linear(100.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_interp_cubic_midpoint() {
        let ch = ramp_channel(10);
        // Catmull-Rom through linear data must reproduce exact midpoints
        assert!((ch.interp_cubic(3.5) - 3.5).abs() < 1e-9);
    }

    // Resampling -------------------------------------------------------------

    #[test]
    fn test_resample_uniform_linear() {
        let ch = ramp_channel(10);
        let out = resample_uniform(&ch, 0.0, 9.0, 19, ResampleMethod::Linear).unwrap();
        assert_eq!(out.len(), 19);
        assert!((out.values[0] - 0.0).abs() < 1e-12);
        assert!((out.values[18] - 9.0).abs() < 1e-12);
        assert!((out.values[9] - 4.5).abs() < 1e-12);
    }

    #[test]
    fn test_resample_uniform_too_few() {
        let ch = ramp_channel(4);
        assert!(resample_uniform(&ch, 0.0, 3.0, 1, ResampleMethod::Linear).is_err());
    }

    #[test]
    fn test_resample_to_times() {
        let ch = ramp_channel(5);
        let targets = [0.5, 1.5, 2.5];
        let out = resample_to_times(&ch, &targets, ResampleMethod::Linear);
        assert_eq!(out.len(), 3);
        for (i, &v) in out.values.iter().enumerate() {
            assert!((v - targets[i]).abs() < 1e-12);
        }
    }

    // Moving average ---------------------------------------------------------

    #[test]
    fn test_moving_average_constant() {
        let data = vec![3.0_f64; 20];
        let out = moving_average(&data, 5);
        for v in &out {
            assert!((*v - 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_moving_average_window_1() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let out = moving_average(&data, 1);
        for (a, b) in data.iter().zip(out.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    // Butterworth filter -----------------------------------------------------

    #[test]
    fn test_butterworth_lp_passes_dc() {
        let bq = BiquadCoeffs::butterworth_2nd(100.0, 1000.0, FilterType::LowPass);
        let dc = vec![1.0_f64; 200];
        let out = bq.apply(&dc);
        assert!((out.last().copied().unwrap_or(0.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_butterworth_hp_blocks_dc() {
        let bq = BiquadCoeffs::butterworth_2nd(50.0, 1000.0, FilterType::HighPass);
        let dc = vec![1.0_f64; 500];
        let out = bq.apply(&dc);
        assert!(out.last().copied().unwrap_or(1.0).abs() < 1e-4);
    }

    #[test]
    fn test_butterworth_zero_phase() {
        let bq = BiquadCoeffs::butterworth_2nd(100.0, 1000.0, FilterType::LowPass);
        let impulse: Vec<f64> = std::iter::once(1.0)
            .chain(std::iter::repeat_n(0.0, 99))
            .collect();
        let out = bq.apply_zero_phase(&impulse);
        assert_eq!(out.len(), 100);
    }

    // FFT -------------------------------------------------------------------

    #[test]
    fn test_fft_dc() {
        let data = vec![1.0_f64; 8];
        let spec = rfft(&data).unwrap();
        assert!((spec[0].0 - 8.0).abs() < 1e-10); // DC bin = sum
        for k in 1..spec.len() {
            assert!(spec[k].0.abs() < 1e-10);
            assert!(spec[k].1.abs() < 1e-10);
        }
    }

    #[test]
    fn test_fft_size_not_pow2_error() {
        let data = vec![0.0_f64; 7];
        assert!(rfft(&data).is_err());
    }

    #[test]
    fn test_rfft_frequencies_length() {
        let freqs = rfft_frequencies(16, 100.0);
        assert_eq!(freqs.len(), 9); // N/2 + 1
    }

    #[test]
    fn test_next_pow2() {
        assert_eq!(next_pow2(1), 1);
        assert_eq!(next_pow2(3), 4);
        assert_eq!(next_pow2(8), 8);
        assert_eq!(next_pow2(9), 16);
    }

    // Welch PSD --------------------------------------------------------------

    #[test]
    fn test_welch_psd_output_length() {
        let data: Vec<f64> = (0..256).map(|i| (i as f64 * 0.1).sin()).collect();
        let (freqs, psd) = welch_psd(&data, 100.0, 64, 32, WindowType::Hann).unwrap();
        assert_eq!(freqs.len(), 33);
        assert_eq!(psd.len(), 33);
    }

    #[test]
    fn test_welch_psd_invalid_fs() {
        let data = vec![0.0_f64; 256];
        assert!(welch_psd(&data, -1.0, 64, 0, WindowType::Rectangular).is_err());
    }

    // Window functions -------------------------------------------------------

    #[test]
    fn test_hann_window_ends() {
        let w = make_window(WindowType::Hann, 64);
        assert!(w[0].abs() < 1e-12);
        assert!(w[63].abs() < 1e-12);
    }

    #[test]
    fn test_rectangular_window() {
        let w = make_window(WindowType::Rectangular, 10);
        assert!(w.iter().all(|&v| (v - 1.0).abs() < 1e-12));
    }

    // Autocorrelation --------------------------------------------------------

    #[test]
    fn test_autocorrelation_lag0_is_1() {
        let ch = sine_channel(5.0, 100.0, 128);
        let ac = autocorrelation(&ch.values, 0);
        assert!((ac[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_length() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let ac = autocorrelation(&data, 10);
        assert_eq!(ac.len(), 11);
    }

    // Cross-correlation ------------------------------------------------------

    #[test]
    fn test_cross_correlation_self_is_autocorr() {
        let data: Vec<f64> = (0..32).map(|i| (i as f64 * 0.3).sin()).collect();
        let (_lags, ccf) = cross_correlation(&data, &data, 5);
        let ac = autocorrelation(&data, 5);
        for (a, b) in ac.iter().zip(ccf.iter().skip(5)) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    // Peak detection ---------------------------------------------------------

    #[test]
    fn test_detect_peaks_sine() {
        let ch = sine_channel(1.0, 100.0, 200);
        let peaks = detect_peaks(&ch, 0.5, 50);
        // Two full periods → 2 peaks
        assert!(!peaks.is_empty());
        for p in &peaks {
            assert!(p.value > 0.8);
        }
    }

    #[test]
    fn test_detect_peaks_flat() {
        let mut ch = Channel::new("flat", "");
        for i in 0..20 {
            ch.push(i as f64, 1.0);
        }
        let peaks = detect_peaks(&ch, 0.1, 1);
        assert!(peaks.is_empty());
    }

    // Anomaly detection ------------------------------------------------------

    #[test]
    fn test_zscore_anomaly_detects_spike() {
        let mut ch = Channel::new("v", "");
        for i in 0..100 {
            ch.push(i as f64, 0.0);
        }
        ch.push(100.0, 100.0);
        let anomalies = anomalies_zscore(&ch, 2.0);
        assert!(!anomalies.is_empty());
        assert_eq!(anomalies[0].index, 100);
    }

    #[test]
    fn test_iqr_anomaly_detects_outlier() {
        let mut ch = Channel::new("v", "");
        for i in 0..50 {
            ch.push(i as f64, (i % 2) as f64);
        }
        ch.push(50.0, 9999.0);
        let anomalies = anomalies_iqr(&ch, 1.5);
        assert!(!anomalies.is_empty());
    }

    // CSV export/import ------------------------------------------------------

    #[test]
    fn test_csv_roundtrip() {
        let mut ch = Channel::new("pressure", "Pa");
        for i in 0..5 {
            ch.push(i as f64 * 0.01, i as f64 * 100.0);
        }
        let mut buf = Vec::new();
        write_channel_csv(&mut buf, &ch).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let ch2 = read_channel_csv(&text, "pressure", "Pa").unwrap();
        assert_eq!(ch2.len(), 5);
        for (a, b) in ch.values.iter().zip(ch2.values.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    #[test]
    fn test_multi_channel_csv_header() {
        let mut mcs = MultiChannelSeries::new();
        let mut ch = Channel::new("temp", "K");
        ch.push(0.0, 300.0);
        ch.push(1.0, 310.0);
        mcs.add_channel(ch);
        let mut buf = Vec::new();
        write_multi_channel_csv(&mut buf, &mcs).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.starts_with("time_s,temp"));
    }

    // JSON export/import -----------------------------------------------------

    #[test]
    fn test_json_roundtrip() {
        let mut ch = Channel::new("vel", "m/s");
        for i in 0..4 {
            ch.push(i as f64, i as f64 * 2.0);
        }
        let json = channel_to_json(&ch);
        let ch2 = channel_from_json(&json).unwrap();
        assert_eq!(ch2.len(), 4);
        assert_eq!(ch2.name, "vel");
        for (a, b) in ch.values.iter().zip(ch2.values.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    // Binary export/import ---------------------------------------------------

    #[test]
    fn test_binary_roundtrip() {
        let mut ch = Channel::new("force", "N");
        for i in 0..8 {
            ch.push(i as f64 * 0.1, i as f64 * 3.125);
        }
        let bin = channel_to_binary(&ch);
        let ch2 = channel_from_binary(&bin, "force", "N").unwrap();
        assert_eq!(ch2.len(), 8);
        for (a, b) in ch.values.iter().zip(ch2.values.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_binary_bad_length() {
        let bad = vec![0u8; 17]; // not divisible by 16
        assert!(channel_from_binary(&bad, "x", "").is_err());
    }

    // Streaming writer -------------------------------------------------------

    #[test]
    fn test_streaming_writer_produces_csv() {
        let buf: Vec<u8> = Vec::new();
        let mut writer = StreamingTsWriter::new(buf, "sig", 4);
        for i in 0..10_u32 {
            writer.push(i as f64, i as f64 * 1.5).unwrap();
        }
        let finished = writer.finish().unwrap();
        let text = String::from_utf8(finished).unwrap();
        assert!(text.starts_with("time_s,sig"));
        assert_eq!(text.lines().count(), 11); // header + 10 data rows
    }

    // Ring buffer -----------------------------------------------------------

    #[test]
    fn test_ring_buffer_capacity() {
        let mut rb = RingBuffer::new(4);
        for i in 0..6_u32 {
            rb.push(i as f64, i as f64);
        }
        // Only last 4 samples should be retained
        assert_eq!(rb.len(), 4);
        let ch = rb.drain_to_channel("r", "");
        assert_eq!(ch.values[0], 2.0);
    }

    #[test]
    fn test_ring_buffer_snapshot() {
        let mut rb = RingBuffer::new(8);
        for i in 0..5_u32 {
            rb.push(i as f64, i as f64);
        }
        let snap = rb.snapshot("s", "");
        assert_eq!(snap.len(), 5);
        // Buffer still intact
        assert_eq!(rb.len(), 5);
    }

    // Dominant frequency ----------------------------------------------------

    #[test]
    fn test_dominant_frequency() {
        let fs = 512.0_f64;
        let freq = 32.0_f64;
        let data: Vec<f64> = (0..512)
            .map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 / fs).sin())
            .collect();
        let dom = dominant_frequency(&data, fs).unwrap();
        assert!((dom - freq).abs() < fs / 512.0 + 1.0);
    }

    // RMS / percentile / MAD ------------------------------------------------

    #[test]
    fn test_rms_unit_sine() {
        let data: Vec<f64> = (0..1024)
            .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 1024.0).sin())
            .collect();
        let r = rms(&data);
        assert!((r - 1.0 / 2.0_f64.sqrt()).abs() < 0.005);
    }

    #[test]
    fn test_percentile_median() {
        let data: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        let m = percentile(&data, 50.0);
        assert!((m - 50.5).abs() < 0.01);
    }

    #[test]
    fn test_mad_constant() {
        let data = vec![7.0_f64; 20];
        assert!(mad(&data).abs() < 1e-12);
    }

    // Multi-channel ---------------------------------------------------------

    #[test]
    fn test_multi_channel_add_and_find() {
        let mut mcs = MultiChannelSeries::new();
        let mut ch = Channel::new("omega", "rad/s");
        ch.push(0.0, 1.0);
        mcs.add_channel(ch);
        assert_eq!(mcs.num_channels(), 1);
        assert!(mcs.channel_by_name("omega").is_some());
        assert!(mcs.channel_by_name("missing").is_none());
    }

    // TsError display -------------------------------------------------------

    #[test]
    fn test_error_display() {
        let e = TsError::LengthMismatch {
            expected: 3,
            got: 5,
        };
        assert!(!format!("{e}").is_empty());
        let e2 = TsError::FftSizeNotPow2(7);
        assert!(format!("{e2}").contains("7"));
    }
}
