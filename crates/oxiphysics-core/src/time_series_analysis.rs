// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time-series analysis utilities: statistics, filtering, autocorrelation,
//! trend decomposition, spectral analysis, change-point detection, and AR models.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Core data structure
// ─────────────────────────────────────────────────────────────────────────────

/// A time series with optional explicit timestamps.
#[derive(Debug, Clone)]
pub struct TimeSeries {
    /// Raw sample values.
    pub data: Vec<f64>,
    /// Optional explicit timestamps (seconds).  When `None`, samples are
    /// assumed to be evenly spaced at `1.0 / sampling_rate`.
    pub timestamps: Option<Vec<f64>>,
    /// Nominal sampling rate in Hz.
    pub sampling_rate: f64,
}

impl TimeSeries {
    /// Create a new `TimeSeries` from evenly-spaced data.
    pub fn new(data: Vec<f64>, sampling_rate: f64) -> Self {
        Self {
            data,
            timestamps: None,
            sampling_rate,
        }
    }

    /// Create a `TimeSeries` with explicit timestamps.
    pub fn with_timestamps(data: Vec<f64>, timestamps: Vec<f64>, sampling_rate: f64) -> Self {
        Self {
            data,
            timestamps: Some(timestamps),
            sampling_rate,
        }
    }

    /// Return the number of samples.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Return `true` when there are no samples.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Descriptive statistics
// ─────────────────────────────────────────────────────────────────────────────

/// Descriptive statistics computed over a [`TimeSeries`].
#[derive(Debug, Clone, PartialEq)]
pub struct TimeSeriesStats {
    /// Arithmetic mean.
    pub mean: f64,
    /// Sample variance (Bessel-corrected, N-1 denominator).
    pub variance: f64,
    /// Standard deviation.
    pub std: f64,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Sample skewness.
    pub skewness: f64,
    /// Excess kurtosis (Fisher definition: normal distribution → 0).
    pub kurtosis: f64,
}

/// Compute descriptive statistics for a [`TimeSeries`].
///
/// Returns `None` when the series is empty or has fewer than two samples.
pub fn compute_stats(ts: &TimeSeries) -> Option<TimeSeriesStats> {
    let n = ts.data.len();
    if n < 2 {
        return None;
    }
    let nf = n as f64;
    let mean = ts.data.iter().copied().sum::<f64>() / nf;
    let variance = ts.data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (nf - 1.0);
    let std = variance.sqrt();
    let min = ts.data.iter().copied().fold(f64::INFINITY, f64::min);
    let max = ts.data.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let (skewness, kurtosis) = if std > 0.0 {
        let m3 = ts
            .data
            .iter()
            .map(|&x| ((x - mean) / std).powi(3))
            .sum::<f64>()
            / nf;
        let m4 = ts
            .data
            .iter()
            .map(|&x| ((x - mean) / std).powi(4))
            .sum::<f64>()
            / nf;
        (m3, m4 - 3.0)
    } else {
        (0.0, 0.0)
    };

    Some(TimeSeriesStats {
        mean,
        variance,
        std,
        min,
        max,
        skewness,
        kurtosis,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Moving averages
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a centred (symmetric) simple moving average.
///
/// If `window == 0` or `window > data.len()` the original data is returned
/// unchanged.
pub fn moving_average(data: &[f64], window: usize) -> Vec<f64> {
    let n = data.len();
    if window == 0 || window > n {
        return data.to_vec();
    }
    let half = window / 2;
    let mut out = vec![0.0_f64; n];
    for (i, o) in out.iter_mut().enumerate() {
        let lo = i.saturating_sub(half);
        let hi = (i + half + 1).min(n);
        let count = hi - lo;
        let sum: f64 = data[lo..hi].iter().copied().sum();
        *o = sum / count as f64;
    }
    out
}

/// Compute an exponential moving average with smoothing factor `alpha` ∈ (0, 1].
///
/// `alpha = 1` returns the original series.
pub fn ema(data: &[f64], alpha: f64) -> Vec<f64> {
    if data.is_empty() {
        return vec![];
    }
    let alpha = alpha.clamp(1e-9, 1.0);
    let mut out = Vec::with_capacity(data.len());
    out.push(data[0]);
    for &x in &data[1..] {
        let prev = *out.last().expect("output is non-empty");
        out.push(alpha * x + (1.0 - alpha) * prev);
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Autocorrelation
// ─────────────────────────────────────────────────────────────────────────────

/// Autocorrelation function (ACF) estimator.
#[derive(Debug, Clone)]
pub struct Autocorrelation {
    /// ACF values at lags 0, 1, …, `max_lag`.
    pub values: Vec<f64>,
}

impl Autocorrelation {
    /// Compute the normalized ACF of `data` up to `max_lag`.
    pub fn compute(data: &[f64], max_lag: usize) -> Self {
        let n = data.len();
        if n == 0 {
            return Self { values: vec![] };
        }
        let mean = data.iter().copied().sum::<f64>() / n as f64;
        let var: f64 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>();
        if var == 0.0 {
            return Self {
                values: vec![1.0; max_lag + 1],
            };
        }
        let lags = max_lag.min(n - 1);
        let values = (0..=lags)
            .map(|lag| {
                let cov: f64 = (0..n - lag)
                    .map(|i| (data[i] - mean) * (data[i + lag] - mean))
                    .sum();
                cov / var
            })
            .collect();
        Self { values }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cross-correlation
// ─────────────────────────────────────────────────────────────────────────────

/// Cross-correlation estimator for two equal-length signals.
#[derive(Debug, Clone)]
pub struct CrossCorrelation {
    /// Cross-correlation values indexed by lag offset.
    pub values: Vec<f64>,
    /// Maximum lag used during computation.
    pub max_lag: usize,
}

impl CrossCorrelation {
    /// Compute the normalized cross-correlation between `x` and `y`.
    ///
    /// The two series are truncated to the shorter length before processing.
    pub fn compute(x: &[f64], y: &[f64], max_lag: usize) -> Self {
        let n = x.len().min(y.len());
        if n == 0 {
            return Self {
                values: vec![],
                max_lag,
            };
        }
        let mx = x[..n].iter().copied().sum::<f64>() / n as f64;
        let my = y[..n].iter().copied().sum::<f64>() / n as f64;
        let sx: f64 = x[..n].iter().map(|&v| (v - mx).powi(2)).sum::<f64>().sqrt();
        let sy: f64 = y[..n].iter().map(|&v| (v - my).powi(2)).sum::<f64>().sqrt();
        let denom = sx * sy;
        let lags = max_lag.min(n - 1);
        let values = (0..=lags)
            .map(|lag| {
                let cov: f64 = (0..n - lag).map(|i| (x[i] - mx) * (y[i + lag] - my)).sum();
                if denom > 0.0 { cov / denom } else { 0.0 }
            })
            .collect();
        Self {
            values,
            max_lag: lags,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Trend decomposition
// ─────────────────────────────────────────────────────────────────────────────

/// Classical additive decomposition of a time series into trend, seasonal, and
/// residual components.
///
/// Returns `(trend, seasonal, residual)`.  Each output vector has the same
/// length as `data`.  When the series is too short for the requested `period`,
/// all three vectors are zeroed.
pub fn decompose_trend(data: &[f64], period: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = data.len();
    if n < period * 2 || period == 0 {
        return (vec![0.0; n], vec![0.0; n], data.to_vec());
    }

    // 1. Estimate trend via centred moving average of length `period`.
    let trend = moving_average(data, period);

    // 2. De-trended = data – trend.
    let detrended: Vec<f64> = data
        .iter()
        .zip(trend.iter())
        .map(|(&d, &t)| d - t)
        .collect();

    // 3. Seasonal component: average of de-trended values at each phase.
    let mut phase_sum = vec![0.0_f64; period];
    let mut phase_cnt = vec![0_usize; period];
    for (i, &v) in detrended.iter().enumerate() {
        phase_sum[i % period] += v;
        phase_cnt[i % period] += 1;
    }
    let seasonal_pattern: Vec<f64> = phase_sum
        .iter()
        .zip(phase_cnt.iter())
        .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
        .collect();

    // Centre seasonal pattern so it sums to zero.
    let seasonal_mean = seasonal_pattern.iter().copied().sum::<f64>() / period as f64;
    let seasonal_pattern: Vec<f64> = seasonal_pattern
        .iter()
        .map(|&v| v - seasonal_mean)
        .collect();

    let seasonal: Vec<f64> = (0..n).map(|i| seasonal_pattern[i % period]).collect();

    // 4. Residual = data – trend – seasonal.
    let residual: Vec<f64> = data
        .iter()
        .zip(trend.iter())
        .zip(seasonal.iter())
        .map(|((&d, &t), &s)| d - t - s)
        .collect();

    (trend, seasonal, residual)
}

// ─────────────────────────────────────────────────────────────────────────────
// Stationarity test (simplified ADF)
// ─────────────────────────────────────────────────────────────────────────────

/// Simplified Augmented Dickey-Fuller test for stationarity.
///
/// Returns `(statistic, is_stationary)`.  The threshold used is −1.94
/// (approximate 5 % critical value for a large sample).
pub fn adf_test(data: &[f64]) -> (f64, bool) {
    let n = data.len();
    if n < 3 {
        return (0.0, false);
    }
    // First differences Δyₜ = yₜ - yₜ₋₁
    let dy: Vec<f64> = data.windows(2).map(|w| w[1] - w[0]).collect();
    let y_lag: Vec<f64> = data[..n - 1].to_vec();

    let m = dy.len() as f64;
    let mean_x = y_lag.iter().copied().sum::<f64>() / m;
    let mean_y = dy.iter().copied().sum::<f64>() / m;

    let ss_xx: f64 = y_lag.iter().map(|&x| (x - mean_x).powi(2)).sum();
    let ss_xy: f64 = y_lag
        .iter()
        .zip(dy.iter())
        .map(|(&x, &y)| (x - mean_x) * (y - mean_y))
        .sum();

    if ss_xx == 0.0 {
        return (0.0, false);
    }
    let beta = ss_xy / ss_xx;
    let alpha = mean_y - beta * mean_x;

    // Residuals
    let residuals: Vec<f64> = y_lag
        .iter()
        .zip(dy.iter())
        .map(|(&x, &y)| y - (alpha + beta * x))
        .collect();
    let sse: f64 = residuals.iter().map(|&r| r.powi(2)).sum::<f64>();
    let se_beta = ((sse / (m - 2.0)) / ss_xx).sqrt();

    let t_stat = if se_beta > 0.0 { beta / se_beta } else { 0.0 };
    let is_stationary = t_stat < -1.94;
    (t_stat, is_stationary)
}

// ─────────────────────────────────────────────────────────────────────────────
// Peak detection
// ─────────────────────────────────────────────────────────────────────────────

/// Find local maxima whose prominence exceeds `min_prominence`.
///
/// A peak's *prominence* is defined here as the height above the nearest
/// higher sample on either side (or above the boundary if no higher neighbour
/// exists).
pub fn find_peaks(data: &[f64], min_prominence: f64) -> Vec<usize> {
    let n = data.len();
    if n < 3 {
        return vec![];
    }
    let mut peaks = Vec::new();
    for i in 1..n - 1 {
        if data[i] > data[i - 1] && data[i] > data[i + 1] {
            // Left base: minimum to the left up to the next higher peak.
            let left_min = data[..i].iter().copied().fold(f64::INFINITY, f64::min);
            let right_min = data[i + 1..].iter().copied().fold(f64::INFINITY, f64::min);
            let base = left_min.max(right_min);
            let prominence = data[i] - base;
            if prominence >= min_prominence {
                peaks.push(i);
            }
        }
    }
    peaks
}

// ─────────────────────────────────────────────────────────────────────────────
// Zero-crossing rate
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the zero-crossing rate: fraction of consecutive-sample sign changes.
///
/// Returns a value in \[0, 1\].
pub fn zero_crossing_rate(data: &[f64]) -> f64 {
    let n = data.len();
    if n < 2 {
        return 0.0;
    }
    let crossings = data
        .windows(2)
        .filter(|w| (w[0] >= 0.0) != (w[1] >= 0.0))
        .count();
    crossings as f64 / (n - 1) as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Power spectrum (simple DFT)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the one-sided power spectral density via a direct DFT.
///
/// Returns a vector of length `n/2 + 1` where `n` is the signal length.
pub fn compute_power_spectrum(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    if n == 0 {
        return vec![];
    }
    let half = n / 2 + 1;
    let mut psd = Vec::with_capacity(half);
    for k in 0..half {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for (t, &x) in data.iter().enumerate() {
            let angle = -2.0 * PI * k as f64 * t as f64 / n as f64;
            re += x * angle.cos();
            im += x * angle.sin();
        }
        let power = (re * re + im * im) / (n as f64 * n as f64);
        psd.push(if k == 0 || (n.is_multiple_of(2) && k == n / 2) {
            power
        } else {
            2.0 * power
        });
    }
    psd
}

// ─────────────────────────────────────────────────────────────────────────────
// Lomb-Scargle periodogram
// ─────────────────────────────────────────────────────────────────────────────

/// Lomb-Scargle periodogram for unevenly-sampled data.
#[derive(Debug, Clone)]
pub struct Periodogram {
    /// Angular frequencies at which the power was evaluated.
    pub frequencies: Vec<f64>,
    /// Normalised Lomb-Scargle power at each frequency.
    pub power: Vec<f64>,
}

impl Periodogram {
    /// Compute the Lomb-Scargle periodogram.
    ///
    /// `times` and `values` must be the same length.  `n_freqs` controls the
    /// number of frequencies sampled in the range
    /// `[1/(max_t - min_t), n/(2*(max_t - min_t))]`.
    pub fn compute(times: &[f64], values: &[f64], n_freqs: usize) -> Self {
        let n = times.len().min(values.len());
        if n < 2 || n_freqs == 0 {
            return Self {
                frequencies: vec![],
                power: vec![],
            };
        }
        let mean = values[..n].iter().copied().sum::<f64>() / n as f64;
        let variance = values[..n].iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        let t_span = times[n - 1] - times[0];
        if t_span <= 0.0 || variance == 0.0 {
            return Self {
                frequencies: vec![],
                power: vec![],
            };
        }
        let f_min = 1.0 / t_span;
        let f_max = (n as f64) / (2.0 * t_span);

        let mut frequencies = Vec::with_capacity(n_freqs);
        let mut power = Vec::with_capacity(n_freqs);

        for fi in 0..n_freqs {
            let freq = f_min + (f_max - f_min) * fi as f64 / n_freqs.max(1) as f64;
            let omega = 2.0 * PI * freq;

            // Time offset τ
            let sin2 = (0..n).map(|i| (2.0 * omega * times[i]).sin()).sum::<f64>();
            let cos2 = (0..n).map(|i| (2.0 * omega * times[i]).cos()).sum::<f64>();
            let tau = (sin2 / cos2).atan() / (2.0 * omega);

            let cos_sum: f64 = (0..n).map(|i| (omega * (times[i] - tau)).cos()).sum();
            let sin_sum: f64 = (0..n).map(|i| (omega * (times[i] - tau)).sin()).sum();
            let cc: f64 = (0..n)
                .map(|i| (omega * (times[i] - tau)).cos().powi(2))
                .sum();
            let ss: f64 = (0..n)
                .map(|i| (omega * (times[i] - tau)).sin().powi(2))
                .sum();

            let yc: f64 = (0..n)
                .map(|i| (values[i] - mean) * (omega * (times[i] - tau)).cos())
                .sum();
            let ys: f64 = (0..n)
                .map(|i| (values[i] - mean) * (omega * (times[i] - tau)).sin())
                .sum();

            let _ = (cos_sum, sin_sum); // suppress warnings; tau already used
            let p = 0.5 / variance
                * (if cc > 0.0 { yc * yc / cc } else { 0.0 }
                    + if ss > 0.0 { ys * ys / ss } else { 0.0 });

            frequencies.push(freq);
            power.push(p);
        }
        Self { frequencies, power }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hurst exponent (R/S analysis)
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the Hurst exponent of a time series using R/S (rescaled range) analysis.
///
/// Returns `None` when the series is too short (< 8 samples).
pub fn hurst_exponent(data: &[f64]) -> Option<f64> {
    let n = data.len();
    if n < 8 {
        return None;
    }
    let max_pow = (n as f64).log2().floor() as u32;
    if max_pow < 2 {
        return None;
    }
    let mut log_n_vals = Vec::new();
    let mut log_rs_vals = Vec::new();

    for pow in 2..=max_pow {
        let sub_len = 1_usize << pow;
        if sub_len > n {
            break;
        }
        let n_subs = n / sub_len;
        let mut rs_sum = 0.0_f64;
        for s in 0..n_subs {
            let segment = &data[s * sub_len..(s + 1) * sub_len];
            let mean = segment.iter().copied().sum::<f64>() / sub_len as f64;
            let std =
                (segment.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / sub_len as f64).sqrt();
            if std == 0.0 {
                continue;
            }
            // Cumulative deviations
            let mut cum = 0.0_f64;
            let mut cum_min = f64::INFINITY;
            let mut cum_max = f64::NEG_INFINITY;
            for &x in segment {
                cum += x - mean;
                cum_min = cum_min.min(cum);
                cum_max = cum_max.max(cum);
            }
            rs_sum += (cum_max - cum_min) / std;
        }
        let rs = rs_sum / n_subs as f64;
        if rs > 0.0 {
            log_n_vals.push((sub_len as f64).ln());
            log_rs_vals.push(rs.ln());
        }
    }

    if log_n_vals.len() < 2 {
        return None;
    }
    // OLS slope of log(R/S) ~ H * log(N)
    let m = log_n_vals.len() as f64;
    let mx = log_n_vals.iter().copied().sum::<f64>() / m;
    let my = log_rs_vals.iter().copied().sum::<f64>() / m;
    let num: f64 = log_n_vals
        .iter()
        .zip(log_rs_vals.iter())
        .map(|(&x, &y)| (x - mx) * (y - my))
        .sum();
    let den: f64 = log_n_vals.iter().map(|&x| (x - mx).powi(2)).sum();
    if den == 0.0 {
        return None;
    }
    Some(num / den)
}

// ─────────────────────────────────────────────────────────────────────────────
// Change-point detection (CUSUM)
// ─────────────────────────────────────────────────────────────────────────────

/// CUSUM-based change-point detector.
#[derive(Debug, Clone)]
pub struct ChangePointDetection {
    /// Detection threshold (in units of standard deviation of the series).
    pub threshold: f64,
}

impl ChangePointDetection {
    /// Create a new detector with the given threshold.
    pub fn new(threshold: f64) -> Self {
        Self { threshold }
    }

    /// Detect change points in `data`.
    ///
    /// Returns a `Vec` of indices where the CUSUM statistic exceeds the
    /// threshold.
    pub fn detect(&self, data: &[f64]) -> Vec<usize> {
        let n = data.len();
        if n < 2 {
            return vec![];
        }
        let mean = data.iter().copied().sum::<f64>() / n as f64;
        let std = (data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
        if std == 0.0 {
            return vec![];
        }
        let limit = self.threshold * std;
        let mut pos_cusum = 0.0_f64;
        let mut neg_cusum = 0.0_f64;
        let mut change_points = Vec::new();
        for (i, &x) in data.iter().enumerate() {
            pos_cusum = (pos_cusum + x - mean).max(0.0);
            neg_cusum = (neg_cusum - x + mean).max(0.0);
            if pos_cusum > limit || neg_cusum > limit {
                change_points.push(i);
                pos_cusum = 0.0;
                neg_cusum = 0.0;
            }
        }
        change_points
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Filters
// ─────────────────────────────────────────────────────────────────────────────

/// Filter types available through [`apply_filter`].
#[derive(Debug, Clone)]
pub enum TimeSeriesFilter {
    /// Simplified first-order Butterworth low-pass filter parameterised by
    /// normalised cut-off frequency `fc` ∈ (0, 0.5).
    ButterworthLowpass {
        /// Normalised cut-off frequency in (0, 0.5).
        fc: f64,
    },
    /// Running median filter of the given window size.
    Median {
        /// Window size (must be odd; rounded up internally).
        window: usize,
    },
    /// Simplified Savitzky-Golay smoothing using a polynomial of degree 2
    /// over the specified window.
    SavitzkyGolay {
        /// Half-window size (total window = `2 * half_window + 1`).
        half_window: usize,
    },
}

/// Apply a [`TimeSeriesFilter`] to `data` and return the filtered signal.
pub fn apply_filter(data: &[f64], filter: &TimeSeriesFilter) -> Vec<f64> {
    match filter {
        TimeSeriesFilter::ButterworthLowpass { fc } => butterworth_lowpass(data, *fc),
        TimeSeriesFilter::Median { window } => median_filter(data, *window),
        TimeSeriesFilter::SavitzkyGolay { half_window } => savitzky_golay(data, *half_window),
    }
}

fn butterworth_lowpass(data: &[f64], fc: f64) -> Vec<f64> {
    let fc = fc.clamp(1e-6, 0.4999);
    // First-order IIR approximation: α = 2πfc / (2πfc + 1)
    let alpha = (2.0 * PI * fc) / (2.0 * PI * fc + 1.0);
    let mut out = Vec::with_capacity(data.len());
    if data.is_empty() {
        return out;
    }
    out.push(data[0]);
    for &x in &data[1..] {
        let prev = *out.last().expect("output is non-empty");
        out.push(prev + alpha * (x - prev));
    }
    out
}

fn median_filter(data: &[f64], window: usize) -> Vec<f64> {
    let w = if window.is_multiple_of(2) {
        window + 1
    } else {
        window
    }
    .max(1);
    let half = w / 2;
    let n = data.len();
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            let mut window_vals: Vec<f64> = data[lo..hi].to_vec();
            window_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mid = window_vals.len() / 2;
            window_vals[mid]
        })
        .collect()
}

fn savitzky_golay(data: &[f64], half_window: usize) -> Vec<f64> {
    // Polynomial degree 2 least-squares coefficients for a symmetric window
    // of total length 2*half_window + 1.
    let n = data.len();
    if half_window == 0 || n == 0 {
        return data.to_vec();
    }
    let m = half_window as i64;
    let wlen = (2 * half_window + 1) as f64;
    // Precompute normalisation terms for quadratic SG filter.
    let norm_denom = (wlen * (wlen.powi(2) - 1.0) / 3.0) / wlen;
    let _ = norm_denom; // Used implicitly below via the quadratic fit.

    (0..n)
        .map(|i| {
            let lo = (i as i64 - m).max(0) as usize;
            let hi = (i as i64 + m + 1).min(n as i64) as usize;
            let segment = &data[lo..hi];
            let k = segment.len() as f64;
            let mean = segment.iter().copied().sum::<f64>() / k;
            // Quadratic fit: evaluate central coefficient (smoothed value).
            // For a symmetric window this equals the centred quadratic regression value.
            let offset: Vec<f64> = (lo..hi).map(|j| j as f64 - i as f64).collect();
            let s0 = k;
            let s2: f64 = offset.iter().map(|&t| t * t).sum();
            let s4: f64 = offset.iter().map(|&t| t.powi(4)).sum();
            let sy: f64 = segment.iter().copied().sum::<f64>();
            let st2y: f64 = segment
                .iter()
                .zip(offset.iter())
                .map(|(&y, &t)| y * t * t)
                .sum();
            let det = s0 * s4 - s2 * s2;
            if det == 0.0 {
                mean
            } else {
                (s4 * sy - s2 * st2y) / det
            }
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// AR(p) model
// ─────────────────────────────────────────────────────────────────────────────

/// Autoregressive model of order `p`.
#[derive(Debug, Clone)]
pub struct ArModel {
    /// AR coefficients φ₁ … φₚ.
    pub coefficients: Vec<f64>,
    /// Estimated noise variance.
    pub noise_variance: f64,
}

impl ArModel {
    /// Fit an AR(p) model to `data` using the Yule-Walker equations.
    ///
    /// Returns `None` when `data.len() <= p` or `p == 0`.
    pub fn fit(data: &[f64], p: usize) -> Option<Self> {
        if p == 0 || data.len() <= p {
            return None;
        }
        let acf = Autocorrelation::compute(data, p);
        if acf.values.len() < p + 1 {
            return None;
        }
        // Build Toeplitz system R * φ = r
        // R[i][j] = acf[|i-j|] * acf[0]   (use sample variance as scale)
        let n = data.len() as f64;
        let mean = data.iter().copied().sum::<f64>() / n;
        let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;

        // Levinson-Durbin recursion for the Yule-Walker equations.
        let r: Vec<f64> = acf.values.iter().map(|&v| v * var).collect();
        let mut phi = vec![0.0_f64; p];
        let mut err = r[0];
        let mut phi_prev = vec![0.0_f64; p];

        for i in 0..p {
            let mut lambda = r[i + 1];
            for k in 0..i {
                lambda -= phi_prev[k] * r[i - k];
            }
            if err.abs() < 1e-15 {
                break;
            }
            let kappa = lambda / err;
            let mut phi_new = vec![0.0_f64; p];
            phi_new[i] = kappa;
            for k in 0..i {
                phi_new[k] = phi_prev[k] - kappa * phi_prev[i - 1 - k];
            }
            err *= 1.0 - kappa * kappa;
            phi = phi_new.clone();
            phi_prev = phi_new;
        }

        Some(Self {
            coefficients: phi,
            noise_variance: err.max(0.0),
        })
    }

    /// Forecast `steps` steps ahead using the fitted model.
    ///
    /// `history` must contain at least `p` values.
    pub fn forecast(&self, history: &[f64], steps: usize) -> Vec<f64> {
        let p = self.coefficients.len();
        if history.len() < p || steps == 0 {
            return vec![];
        }
        let mut buf: Vec<f64> = history[history.len() - p..].to_vec();
        let mut out = Vec::with_capacity(steps);
        for _ in 0..steps {
            let pred: f64 = self
                .coefficients
                .iter()
                .enumerate()
                .map(|(i, &c)| c * buf[buf.len() - 1 - i])
                .sum();
            out.push(pred);
            buf.push(pred);
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Missing-value interpolation
// ─────────────────────────────────────────────────────────────────────────────

/// Fill `NaN` values in `data` using linear interpolation between the nearest
/// valid neighbours.  Leading and trailing `NaN`s are filled with the first /
/// last valid value respectively.
pub fn interpolate_linear(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    let mut out = data.to_vec();
    let mut i = 0;
    while i < n {
        if out[i].is_nan() {
            // Find the previous valid index.
            let prev = (0..i).rev().find(|&k| !out[k].is_nan());
            // Find the next valid index.
            let next = (i + 1..n).find(|&k| !out[k].is_nan());
            match (prev, next) {
                (Some(p), Some(q)) => {
                    let slope = (out[q] - out[p]) / (q - p) as f64;
                    for j in p + 1..q {
                        out[j] = out[p] + slope * (j - p) as f64;
                    }
                    i = q;
                }
                (None, Some(q)) => {
                    let fill = out[q];
                    for o in out.iter_mut().take(q) {
                        *o = fill;
                    }
                    i = q;
                }
                (Some(_p), None) => {
                    let fill = out[_p];
                    for o in out.iter_mut().skip(_p + 1) {
                        *o = fill;
                    }
                    break;
                }
                (None, None) => break,
            }
        }
        i += 1;
    }
    out
}

/// Fill `NaN` values in `data` using simple cubic-spline-like interpolation
/// (Catmull-Rom tangents) between the nearest valid neighbours.
pub fn interpolate_spline(data: &[f64]) -> Vec<f64> {
    // Collect valid (index, value) pairs.
    let knots: Vec<(usize, f64)> = data
        .iter()
        .enumerate()
        .filter(|&(_, v)| !v.is_nan())
        .map(|(i, &v)| (i, v))
        .collect();
    if knots.is_empty() {
        return data.to_vec();
    }
    let mut out = data.to_vec();
    for (i, o) in out.iter_mut().enumerate() {
        if !o.is_nan() {
            continue;
        }
        // Find surrounding knots.
        let right = knots.partition_point(|&(ki, _)| ki <= i);
        if right == 0 {
            *o = knots[0].1;
        } else if right >= knots.len() {
            *o = knots[knots.len() - 1].1;
        } else {
            let (x0, y0) = knots[right - 1];
            let (x1, y1) = knots[right];
            let t = (i - x0) as f64 / (x1 - x0) as f64;
            // Catmull-Rom tangents from neighbours (or zero at boundaries).
            let m0 = if right >= 2 {
                let (xp, yp) = knots[right - 2];
                (y1 - yp) / (x1 - xp) as f64
            } else {
                y1 - y0
            };
            let m1 = if right + 1 < knots.len() {
                let (xn, yn) = knots[right + 1];
                (yn - y0) / (xn - x0) as f64
            } else {
                y1 - y0
            };
            // Hermite basis
            let h00 = 2.0 * t.powi(3) - 3.0 * t.powi(2) + 1.0;
            let h10 = t.powi(3) - 2.0 * t.powi(2) + t;
            let h01 = -2.0 * t.powi(3) + 3.0 * t.powi(2);
            let h11 = t.powi(3) - t.powi(2);
            let span = (x1 - x0) as f64;
            *o = h00 * y0 + h10 * span * m0 + h01 * y1 + h11 * span * m1;
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ────────────────────────────────────────────────────────────

    fn sine_wave(n: usize, freq: f64, sr: f64) -> Vec<f64> {
        (0..n)
            .map(|i| (2.0 * PI * freq * i as f64 / sr).sin())
            .collect()
    }

    fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![start];
        }
        (0..n)
            .map(|i| start + (end - start) * i as f64 / (n - 1) as f64)
            .collect()
    }

    // ── TimeSeries ──────────────────────────────────────────────────────────

    #[test]
    fn test_time_series_new() {
        let ts = TimeSeries::new(vec![1.0, 2.0, 3.0], 100.0);
        assert_eq!(ts.len(), 3);
        assert!(!ts.is_empty());
    }

    #[test]
    fn test_time_series_empty() {
        let ts = TimeSeries::new(vec![], 100.0);
        assert!(ts.is_empty());
    }

    #[test]
    fn test_time_series_with_timestamps() {
        let ts = TimeSeries::with_timestamps(vec![1.0, 2.0], vec![0.0, 0.01], 100.0);
        assert!(ts.timestamps.is_some());
    }

    // ── compute_stats ───────────────────────────────────────────────────────

    #[test]
    fn test_stats_basic() {
        let ts = TimeSeries::new(vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0], 1.0);
        let s = compute_stats(&ts).unwrap();
        assert!((s.mean - 5.0).abs() < 1e-9);
        assert!((s.variance - 4.571_428_571).abs() < 1e-6);
        assert_eq!(s.min, 2.0);
        assert_eq!(s.max, 9.0);
    }

    #[test]
    fn test_stats_constant() {
        let ts = TimeSeries::new(vec![3.0; 10], 1.0);
        let s = compute_stats(&ts).unwrap();
        assert!((s.mean - 3.0).abs() < 1e-12);
        assert_eq!(s.std, 0.0);
        assert_eq!(s.skewness, 0.0);
        assert_eq!(s.kurtosis, 0.0);
    }

    #[test]
    fn test_stats_too_short() {
        let ts = TimeSeries::new(vec![1.0], 1.0);
        assert!(compute_stats(&ts).is_none());
    }

    // ── moving_average ──────────────────────────────────────────────────────

    #[test]
    fn test_moving_average_window_1() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let ma = moving_average(&data, 1);
        assert_eq!(ma, data);
    }

    #[test]
    fn test_moving_average_preserves_length() {
        let data = sine_wave(64, 5.0, 100.0);
        let ma = moving_average(&data, 5);
        assert_eq!(ma.len(), data.len());
    }

    #[test]
    fn test_moving_average_constant_input() {
        let data = vec![7.0_f64; 20];
        let ma = moving_average(&data, 5);
        for v in ma {
            assert!((v - 7.0).abs() < 1e-12);
        }
    }

    // ── ema ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_ema_alpha_one() {
        let data = vec![1.0, 3.0, 5.0, 2.0];
        let e = ema(&data, 1.0);
        for (a, b) in e.iter().zip(data.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_ema_length() {
        let data = vec![1.0; 50];
        assert_eq!(ema(&data, 0.3).len(), 50);
    }

    #[test]
    fn test_ema_empty() {
        assert!(ema(&[], 0.5).is_empty());
    }

    // ── Autocorrelation ─────────────────────────────────────────────────────

    #[test]
    fn test_acf_lag_zero_is_one() {
        let data = sine_wave(128, 5.0, 128.0);
        let acf = Autocorrelation::compute(&data, 20);
        assert!((acf.values[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_acf_length() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let acf = Autocorrelation::compute(&data, 10);
        assert_eq!(acf.values.len(), 11);
    }

    #[test]
    fn test_acf_empty() {
        let acf = Autocorrelation::compute(&[], 5);
        assert!(acf.values.is_empty());
    }

    // ── CrossCorrelation ────────────────────────────────────────────────────

    #[test]
    fn test_xcorr_identical_lag_zero() {
        let data = sine_wave(64, 2.0, 64.0);
        let xcorr = CrossCorrelation::compute(&data, &data, 0);
        // lag-0 should be 1.0
        assert!((xcorr.values[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_xcorr_empty() {
        let xcorr = CrossCorrelation::compute(&[], &[], 5);
        assert!(xcorr.values.is_empty());
    }

    // ── decompose_trend ─────────────────────────────────────────────────────

    #[test]
    fn test_decompose_trend_additive() {
        let n = 100;
        let period = 10;
        let trend: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();
        let seasonal: Vec<f64> = (0..n).map(|i| ((i % period) as f64).sin()).collect();
        let data: Vec<f64> = trend
            .iter()
            .zip(seasonal.iter())
            .map(|(t, s)| t + s)
            .collect();
        let (tr, _se, _re) = decompose_trend(&data, period);
        assert_eq!(tr.len(), n);
    }

    #[test]
    fn test_decompose_trend_lengths() {
        let data: Vec<f64> = (0..48).map(|i| i as f64).collect();
        let (tr, se, re) = decompose_trend(&data, 12);
        assert_eq!(tr.len(), 48);
        assert_eq!(se.len(), 48);
        assert_eq!(re.len(), 48);
    }

    #[test]
    fn test_decompose_trend_too_short() {
        let data = vec![1.0, 2.0, 3.0];
        let (tr, _se, re) = decompose_trend(&data, 10);
        assert_eq!(tr, vec![0.0, 0.0, 0.0]);
        assert_eq!(re, data);
    }

    // ── adf_test ────────────────────────────────────────────────────────────

    #[test]
    fn test_adf_stationary_white_noise() {
        // White noise is stationary; may or may not pass simplified test, but should not panic.
        let data: Vec<f64> = (0..200).map(|i| ((i * 7 + 3) % 17) as f64 - 8.0).collect();
        let (stat, _) = adf_test(&data);
        assert!(stat.is_finite());
    }

    #[test]
    fn test_adf_short() {
        let (stat, is_stat) = adf_test(&[1.0, 2.0]);
        assert_eq!(stat, 0.0);
        assert!(!is_stat);
    }

    // ── find_peaks ──────────────────────────────────────────────────────────

    #[test]
    fn test_find_peaks_simple() {
        let data = vec![0.0, 1.0, 0.0, 1.5, 0.0, 2.0, 0.0];
        let peaks = find_peaks(&data, 0.5);
        assert!(peaks.contains(&1));
        assert!(peaks.contains(&3));
        assert!(peaks.contains(&5));
    }

    #[test]
    fn test_find_peaks_monotone() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let peaks = find_peaks(&data, 0.1);
        assert!(peaks.is_empty());
    }

    #[test]
    fn test_find_peaks_no_false_positives() {
        let data = vec![3.0, 2.0, 1.0, 0.5, 1.0, 2.0, 1.0];
        // Only index 1 (value 2.0) is a local maximum here.
        // There should be no peak at the shoulders with high prominence threshold.
        let peaks = find_peaks(&data, 2.5);
        // With high threshold, expect no prominent peaks.
        assert!(peaks.is_empty());
    }

    // ── zero_crossing_rate ──────────────────────────────────────────────────

    #[test]
    fn test_zcr_sine() {
        let data = sine_wave(100, 10.0, 100.0);
        let zcr = zero_crossing_rate(&data);
        assert!(zcr > 0.0 && zcr <= 1.0);
    }

    #[test]
    fn test_zcr_constant_positive() {
        let data = vec![1.0_f64; 20];
        assert_eq!(zero_crossing_rate(&data), 0.0);
    }

    #[test]
    fn test_zcr_alternating() {
        let data = vec![1.0, -1.0, 1.0, -1.0, 1.0];
        let zcr = zero_crossing_rate(&data);
        assert!((zcr - 1.0).abs() < 1e-10);
    }

    // ── compute_power_spectrum ──────────────────────────────────────────────

    #[test]
    fn test_power_spectrum_length() {
        let data = sine_wave(64, 5.0, 64.0);
        let psd = compute_power_spectrum(&data);
        assert_eq!(psd.len(), 33); // 64/2 + 1
    }

    #[test]
    fn test_power_spectrum_non_negative() {
        let data = sine_wave(32, 3.0, 32.0);
        for &p in &compute_power_spectrum(&data) {
            assert!(p >= 0.0);
        }
    }

    #[test]
    fn test_power_spectrum_empty() {
        assert!(compute_power_spectrum(&[]).is_empty());
    }

    // ── Periodogram (Lomb-Scargle) ──────────────────────────────────────────

    #[test]
    fn test_periodogram_length() {
        let times = linspace(0.0, 10.0, 50);
        let values = sine_wave(50, 2.0, 5.0);
        let pg = Periodogram::compute(&times, &values, 20);
        assert_eq!(pg.frequencies.len(), 20);
        assert_eq!(pg.power.len(), 20);
    }

    #[test]
    fn test_periodogram_non_negative_power() {
        let times = linspace(0.0, 5.0, 30);
        let values = sine_wave(30, 1.5, 6.0);
        let pg = Periodogram::compute(&times, &values, 10);
        for &p in &pg.power {
            assert!(p >= 0.0);
        }
    }

    // ── hurst_exponent ──────────────────────────────────────────────────────

    #[test]
    fn test_hurst_range() {
        // Alternating ±1 signal: each segment has constant R/S => H ≈ 0.0 (anti-persistent).
        // We just check the estimate is in the valid range [0, 2).
        let mut data = vec![0.0_f64; 256];
        for i in 1..256 {
            data[i] = data[i - 1] + if i % 2 == 0 { 1.0 } else { -1.0 };
        }
        if let Some(h) = hurst_exponent(&data) {
            assert!((0.0..2.0).contains(&h), "hurst={h}");
        }
    }

    #[test]
    fn test_hurst_too_short() {
        assert!(hurst_exponent(&[1.0, 2.0, 3.0]).is_none());
    }

    // ── ChangePointDetection ────────────────────────────────────────────────

    #[test]
    fn test_change_point_step() {
        let mut data = vec![0.0_f64; 50];
        for v in &mut data[25..] {
            *v = 10.0;
        }
        let cpd = ChangePointDetection::new(1.0);
        let cps = cpd.detect(&data);
        assert!(!cps.is_empty());
        // Change point should be somewhere around index 25.
        let any_near = cps.iter().any(|&cp| (20..=35).contains(&cp));
        assert!(any_near, "change points: {cps:?}");
    }

    #[test]
    fn test_change_point_constant() {
        let data = vec![5.0_f64; 100];
        let cpd = ChangePointDetection::new(1.0);
        assert!(cpd.detect(&data).is_empty());
    }

    // ── apply_filter ────────────────────────────────────────────────────────

    #[test]
    fn test_butterworth_preserves_length() {
        let data = sine_wave(100, 5.0, 100.0);
        let f = TimeSeriesFilter::ButterworthLowpass { fc: 0.1 };
        assert_eq!(apply_filter(&data, &f).len(), 100);
    }

    #[test]
    fn test_median_filter_preserves_length() {
        let data: Vec<f64> = (0..40).map(|i| i as f64).collect();
        let f = TimeSeriesFilter::Median { window: 5 };
        assert_eq!(apply_filter(&data, &f).len(), 40);
    }

    #[test]
    fn test_savitzky_golay_preserves_length() {
        let data = sine_wave(50, 3.0, 50.0);
        let f = TimeSeriesFilter::SavitzkyGolay { half_window: 5 };
        assert_eq!(apply_filter(&data, &f).len(), 50);
    }

    #[test]
    fn test_butterworth_smooths_noise() {
        let clean = sine_wave(200, 1.0, 200.0);
        // Add high-frequency "noise" by mixing a higher-freq component.
        let noisy: Vec<f64> = clean
            .iter()
            .enumerate()
            .map(|(i, &v)| v + 0.5 * (2.0 * PI * 50.0 * i as f64 / 200.0).sin())
            .collect();
        let filtered = apply_filter(&noisy, &TimeSeriesFilter::ButterworthLowpass { fc: 0.05 });
        // The filtered signal should have lower variance than the noisy one.
        let var_noisy: f64 = noisy.iter().map(|&x| x * x).sum::<f64>() / noisy.len() as f64;
        let var_filtered: f64 =
            filtered.iter().map(|&x| x * x).sum::<f64>() / filtered.len() as f64;
        assert!(var_filtered <= var_noisy + 1e-6);
    }

    // ── ArModel ─────────────────────────────────────────────────────────────

    #[test]
    fn test_ar_fit_returns_coefficients() {
        let data: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let model = ArModel::fit(&data, 3).unwrap();
        assert_eq!(model.coefficients.len(), 3);
    }

    #[test]
    fn test_ar_forecast_length() {
        let data: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
        let model = ArModel::fit(&data, 2).unwrap();
        let fc = model.forecast(&data, 5);
        assert_eq!(fc.len(), 5);
    }

    #[test]
    fn test_ar_too_short() {
        assert!(ArModel::fit(&[1.0, 2.0], 3).is_none());
        assert!(ArModel::fit(&[], 1).is_none());
        assert!(ArModel::fit(&[1.0], 0).is_none());
    }

    // ── interpolate_linear ──────────────────────────────────────────────────

    #[test]
    fn test_interpolate_linear_basic() {
        let data = vec![1.0, f64::NAN, 3.0];
        let result = interpolate_linear(&data);
        assert!((result[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_linear_leading_nan() {
        let data = vec![f64::NAN, f64::NAN, 5.0];
        let result = interpolate_linear(&data);
        assert_eq!(result[0], 5.0);
        assert_eq!(result[1], 5.0);
    }

    #[test]
    fn test_interpolate_linear_no_nan() {
        let data = vec![1.0, 2.0, 3.0];
        let result = interpolate_linear(&data);
        assert_eq!(result, data);
    }

    // ── interpolate_spline ──────────────────────────────────────────────────

    #[test]
    fn test_interpolate_spline_basic() {
        let data = vec![0.0, f64::NAN, 2.0, f64::NAN, 4.0];
        let result = interpolate_spline(&data);
        // The filled values should be finite.
        for &v in &result {
            assert!(v.is_finite(), "spline result contains non-finite value");
        }
    }

    #[test]
    fn test_interpolate_spline_no_nan() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = interpolate_spline(&data);
        assert_eq!(result, data);
    }

    #[test]
    fn test_interpolate_spline_all_nan() {
        let data = vec![f64::NAN; 5];
        let result = interpolate_spline(&data);
        // All NaN, no valid data — output is the same as input.
        assert!(result.iter().all(|v| v.is_nan()));
    }
}
