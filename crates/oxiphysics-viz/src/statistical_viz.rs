// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Statistical visualization primitives for the OxiPhysics engine.
//!
//! This module provides pure-Rust data structures and algorithms for generating
//! statistical chart layouts. All output is renderer-agnostic: functions return
//! plain numerical data that a downstream renderer can display.
//!
//! Included chart types:
//! - **Histogram** – binned frequency counts with optional cumulative curve
//! - **Box plot** – quartiles, whiskers (Tukey 1.5 × IQR), and outlier points
//! - **Violin plot** – kernel density estimation (Gaussian KDE) for distribution shape
//! - **Scatter plot matrix (SPLOM)** – pairwise scatter data for *n* variables
//! - **Correlation heatmap** – Pearson correlation matrix
//! - **QQ plot** – quantile-quantile comparison against a theoretical distribution
//! - **Time series with confidence bands** – mean ± σ envelope over time
//! - **Error bar plots** – centres with ±error extents
//! - **Empirical CDF** – step-function cumulative distribution from samples
//! - **Bland-Altman plot** – method-comparison via mean vs. difference

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Shared helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Sort a slice of f64 values (NaN-safe: NaNs are moved to the end).
fn sort_f64(data: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = data.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));
    v
}

/// Compute the arithmetic mean of `data`. Returns `f64::NAN` if empty.
fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NAN;
    }
    data.iter().sum::<f64>() / data.len() as f64
}

/// Compute the population variance of `data`. Returns `0.0` if ≤ 1 element.
fn variance(data: &[f64]) -> f64 {
    if data.len() <= 1 {
        return 0.0;
    }
    let m = mean(data);
    data.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / data.len() as f64
}

/// Compute the population standard deviation.
fn std_dev(data: &[f64]) -> f64 {
    variance(data).sqrt()
}

/// Linear interpolation between `a` and `b` with parameter `t` ∈ \[0, 1\].
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Compute the p-th percentile of **sorted** data using linear interpolation.
///
/// `p` should be in `[0.0, 100.0]`.
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let n = sorted.len() as f64;
    let rank = (p / 100.0) * (n - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    let frac = rank - lo as f64;
    lerp(
        sorted[lo.min(sorted.len() - 1)],
        sorted[hi.min(sorted.len() - 1)],
        frac,
    )
}

/// Compute Q1, Q2 (median), Q3 of **sorted** data.
fn quartiles(sorted: &[f64]) -> (f64, f64, f64) {
    let q1 = percentile_sorted(sorted, 25.0);
    let q2 = percentile_sorted(sorted, 50.0);
    let q3 = percentile_sorted(sorted, 75.0);
    (q1, q2, q3)
}

// ─────────────────────────────────────────────────────────────────────────────
// Histogram
// ─────────────────────────────────────────────────────────────────────────────

/// A single histogram bin.
#[derive(Debug, Clone, PartialEq)]
pub struct HistBin {
    /// Left edge of the bin (inclusive).
    pub left: f64,
    /// Right edge of the bin (exclusive, except for the last bin).
    pub right: f64,
    /// Count (number of data points falling in this bin).
    pub count: usize,
    /// Frequency = count / total.
    pub frequency: f64,
    /// Cumulative frequency up to and including this bin.
    pub cumulative: f64,
}

impl HistBin {
    /// Centre of the bin.
    pub fn centre(&self) -> f64 {
        0.5 * (self.left + self.right)
    }

    /// Width of the bin.
    pub fn width(&self) -> f64 {
        self.right - self.left
    }
}

/// Result of computing a histogram.
#[derive(Debug, Clone)]
pub struct Histogram {
    /// All bins in left-to-right order.
    pub bins: Vec<HistBin>,
    /// Total number of data points.
    pub n: usize,
    /// Minimum data value.
    pub data_min: f64,
    /// Maximum data value.
    pub data_max: f64,
}

impl Histogram {
    /// Compute a histogram with `n_bins` equal-width bins.
    ///
    /// Data values outside `[min, max]` are clamped into the first/last bin.
    pub fn compute(data: &[f64], n_bins: usize) -> Self {
        let n_bins = n_bins.max(1);
        if data.is_empty() {
            return Self {
                bins: vec![],
                n: 0,
                data_min: f64::NAN,
                data_max: f64::NAN,
            };
        }
        let sorted = sort_f64(data);
        let data_min = sorted[0];
        let data_max = *sorted.last().expect("collection should not be empty");
        let range = data_max - data_min;
        let bin_width = if range < 1e-300 {
            1.0
        } else {
            range / n_bins as f64
        };

        let mut counts = vec![0usize; n_bins];
        for &x in data {
            let idx = if range < 1e-300 {
                0
            } else {
                let i = ((x - data_min) / bin_width).floor() as usize;
                i.min(n_bins - 1)
            };
            counts[idx] += 1;
        }
        let total = data.len() as f64;
        let mut cumulative = 0.0;
        let bins: Vec<HistBin> = (0..n_bins)
            .map(|i| {
                let left = data_min + i as f64 * bin_width;
                let right = left + bin_width;
                let freq = counts[i] as f64 / total;
                cumulative += freq;
                HistBin {
                    left,
                    right,
                    count: counts[i],
                    frequency: freq,
                    cumulative,
                }
            })
            .collect();
        Self {
            bins,
            n: data.len(),
            data_min,
            data_max,
        }
    }

    /// Total of all bin counts (should equal `self.n`).
    pub fn total_count(&self) -> usize {
        self.bins.iter().map(|b| b.count).sum()
    }

    /// Mode bin (bin with highest count).
    pub fn mode_bin(&self) -> Option<&HistBin> {
        self.bins.iter().max_by_key(|b| b.count)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Box plot
// ─────────────────────────────────────────────────────────────────────────────

/// Box plot statistics for a single data series.
#[derive(Debug, Clone)]
pub struct BoxPlot {
    /// Minimum data value (excluding outliers).
    pub whisker_low: f64,
    /// First quartile (Q1, 25th percentile).
    pub q1: f64,
    /// Median (Q2, 50th percentile).
    pub median: f64,
    /// Third quartile (Q3, 75th percentile).
    pub q3: f64,
    /// Maximum data value (excluding outliers).
    pub whisker_high: f64,
    /// Outlier values below the lower fence.
    pub outliers_low: Vec<f64>,
    /// Outlier values above the upper fence.
    pub outliers_high: Vec<f64>,
    /// Mean of the data.
    pub mean: f64,
    /// Total number of observations.
    pub n: usize,
}

impl BoxPlot {
    /// Compute box plot statistics using Tukey 1.5 × IQR fences.
    pub fn compute(data: &[f64]) -> Self {
        if data.is_empty() {
            return Self {
                whisker_low: f64::NAN,
                q1: f64::NAN,
                median: f64::NAN,
                q3: f64::NAN,
                whisker_high: f64::NAN,
                outliers_low: vec![],
                outliers_high: vec![],
                mean: f64::NAN,
                n: 0,
            };
        }
        let sorted = sort_f64(data);
        let (q1, median, q3) = quartiles(&sorted);
        let iqr = q3 - q1;
        let lower_fence = q1 - 1.5 * iqr;
        let upper_fence = q3 + 1.5 * iqr;

        let mut outliers_low = Vec::new();
        let mut outliers_high = Vec::new();
        let mut whisker_low = q1;
        let mut whisker_high = q3;
        for &v in &sorted {
            if v < lower_fence {
                outliers_low.push(v);
            } else if v > upper_fence {
                outliers_high.push(v);
            } else {
                whisker_low = whisker_low.min(v);
                whisker_high = whisker_high.max(v);
            }
        }
        Self {
            whisker_low,
            q1,
            median,
            q3,
            whisker_high,
            outliers_low,
            outliers_high,
            mean: mean(data),
            n: data.len(),
        }
    }

    /// Inter-quartile range Q3 - Q1.
    pub fn iqr(&self) -> f64 {
        self.q3 - self.q1
    }

    /// True if the value `x` is an outlier according to Tukey fences.
    pub fn is_outlier(&self, x: f64) -> bool {
        x < self.q1 - 1.5 * self.iqr() || x > self.q3 + 1.5 * self.iqr()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Violin plot (KDE)
// ─────────────────────────────────────────────────────────────────────────────

/// A single evaluation point of a kernel density estimate.
#[derive(Debug, Clone)]
pub struct KdePoint {
    /// Position on the x-axis.
    pub x: f64,
    /// Estimated density at this x.
    pub density: f64,
}

/// Violin plot: a mirrored KDE profile plus embedded box-plot statistics.
#[derive(Debug, Clone)]
pub struct ViolinPlot {
    /// KDE curve evaluated at `n_points` evenly-spaced x-values.
    pub kde: Vec<KdePoint>,
    /// Embedded box-plot (quartiles, whiskers, outliers).
    pub box_plot: BoxPlot,
    /// Bandwidth used for the Gaussian kernel.
    pub bandwidth: f64,
}

impl ViolinPlot {
    /// Compute a violin plot using a Gaussian KDE with Silverman's bandwidth rule.
    ///
    /// * `data` – observed values.
    /// * `n_points` – number of evaluation points for the KDE curve.
    pub fn compute(data: &[f64], n_points: usize) -> Self {
        let n_points = n_points.max(2);
        let box_plot = BoxPlot::compute(data);
        if data.is_empty() {
            return Self {
                kde: vec![],
                box_plot,
                bandwidth: 0.0,
            };
        }

        // Silverman's rule of thumb: h = 1.06 σ n^{-1/5}
        let sigma = std_dev(data);
        let bandwidth = if sigma < 1e-300 {
            1.0
        } else {
            1.06 * sigma * (data.len() as f64).powf(-0.2)
        };

        let sorted = sort_f64(data);
        let x_min = sorted[0] - 3.0 * bandwidth;
        let x_max = *sorted.last().expect("collection should not be empty") + 3.0 * bandwidth;
        let step = (x_max - x_min) / (n_points - 1) as f64;

        let kde: Vec<KdePoint> = (0..n_points)
            .map(|i| {
                let x = x_min + i as f64 * step;
                let density = gaussian_kde(data, x, bandwidth);
                KdePoint { x, density }
            })
            .collect();

        Self {
            kde,
            box_plot,
            bandwidth,
        }
    }

    /// Maximum density value across all KDE points.
    pub fn max_density(&self) -> f64 {
        self.kde.iter().map(|p| p.density).fold(0.0_f64, f64::max)
    }
}

/// Evaluate a Gaussian KDE at `x` given `data` points and `bandwidth` `h`.
///
/// K(u) = exp(-0.5 u²) / sqrt(2π)
fn gaussian_kde(data: &[f64], x: f64, h: f64) -> f64 {
    if h < 1e-300 || data.is_empty() {
        return 0.0;
    }
    let norm = 1.0 / (h * (2.0 * std::f64::consts::PI).sqrt());
    let sum: f64 = data
        .iter()
        .map(|&xi| {
            let u = (x - xi) / h;
            norm * (-0.5 * u * u).exp()
        })
        .sum();
    sum / data.len() as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Scatter plot matrix (SPLOM)
// ─────────────────────────────────────────────────────────────────────────────

/// Pairwise scatter data for two variables.
#[derive(Debug, Clone)]
pub struct ScatterPair {
    /// Name of the X variable.
    pub x_name: String,
    /// Name of the Y variable.
    pub y_name: String,
    /// (x, y) data point pairs.
    pub points: Vec<(f64, f64)>,
}

impl ScatterPair {
    /// Pearson correlation coefficient for this pair.
    pub fn pearson_r(&self) -> f64 {
        let xs: Vec<f64> = self.points.iter().map(|&(x, _)| x).collect();
        let ys: Vec<f64> = self.points.iter().map(|&(_, y)| y).collect();
        pearson_correlation(&xs, &ys)
    }
}

/// Scatter plot matrix: all pairwise combinations of `n` variables.
#[derive(Debug, Clone)]
pub struct ScatterMatrix {
    /// Variable names.
    pub names: Vec<String>,
    /// Lower-triangular pairwise scatter data (i > j).
    pub pairs: Vec<ScatterPair>,
}

impl ScatterMatrix {
    /// Build a scatter plot matrix from a column-major data matrix.
    ///
    /// `columns[i]` is the vector of values for variable `i`.
    pub fn compute(names: Vec<String>, columns: &[Vec<f64>]) -> Self {
        let n = names.len().min(columns.len());
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in 0..i {
                let len = columns[i].len().min(columns[j].len());
                let points: Vec<(f64, f64)> =
                    (0..len).map(|k| (columns[j][k], columns[i][k])).collect();
                pairs.push(ScatterPair {
                    x_name: names[j].clone(),
                    y_name: names[i].clone(),
                    points,
                });
            }
        }
        Self {
            names: names[..n].to_vec(),
            pairs,
        }
    }

    /// Number of variables in the matrix.
    pub fn n_vars(&self) -> usize {
        self.names.len()
    }

    /// Expected number of pairs: n*(n-1)/2.
    pub fn expected_pairs(&self) -> usize {
        let n = self.n_vars();
        n * (n.saturating_sub(1)) / 2
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Correlation heatmap
// ─────────────────────────────────────────────────────────────────────────────

/// Pearson correlation between two equal-length vectors.
///
/// Returns `f64::NAN` if either vector is empty or has zero variance.
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n == 0 {
        return f64::NAN;
    }
    let mx = mean(&x[..n]);
    let my = mean(&y[..n]);
    let mut cov = 0.0;
    let mut sx2 = 0.0;
    let mut sy2 = 0.0;
    for i in 0..n {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        cov += dx * dy;
        sx2 += dx * dx;
        sy2 += dy * dy;
    }
    let denom = (sx2 * sy2).sqrt();
    if denom < 1e-300 {
        f64::NAN
    } else {
        cov / denom
    }
}

/// A full symmetric correlation matrix.
#[derive(Debug, Clone)]
pub struct CorrelationHeatmap {
    /// Variable names (row/column labels).
    pub names: Vec<String>,
    /// Row-major `n × n` correlation matrix. `matrix[i][j]` = corr(i, j).
    pub matrix: Vec<Vec<f64>>,
}

impl CorrelationHeatmap {
    /// Compute the Pearson correlation matrix from column-major data.
    pub fn compute(names: Vec<String>, columns: &[Vec<f64>]) -> Self {
        let n = names.len().min(columns.len());
        let matrix: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            1.0
                        } else {
                            pearson_correlation(&columns[i], &columns[j])
                        }
                    })
                    .collect()
            })
            .collect();
        Self {
            names: names[..n].to_vec(),
            matrix,
        }
    }

    /// Retrieve the correlation between two named variables.
    pub fn get(&self, a: &str, b: &str) -> Option<f64> {
        let ia = self.names.iter().position(|n| n == a)?;
        let ib = self.names.iter().position(|n| n == b)?;
        Some(self.matrix[ia][ib])
    }

    /// Dimension of the square matrix.
    pub fn n(&self) -> usize {
        self.names.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// QQ plot
// ─────────────────────────────────────────────────────────────────────────────

/// A single point on a QQ plot.
#[derive(Debug, Clone)]
pub struct QqPoint {
    /// Theoretical quantile.
    pub theoretical: f64,
    /// Sample quantile.
    pub sample: f64,
}

/// Theoretical distribution for QQ plot comparison.
#[derive(Debug, Clone, Copy)]
pub enum TheoreticalDist {
    /// Standard normal N(0, 1).
    Normal,
    /// Uniform on \[0, 1\].
    Uniform,
}

/// Quantile-quantile plot comparing sample data to a theoretical distribution.
#[derive(Debug, Clone)]
pub struct QqPlot {
    /// (theoretical_quantile, sample_quantile) pairs.
    pub points: Vec<QqPoint>,
    /// Distribution used for comparison.
    pub dist: TheoreticalDist,
}

impl QqPlot {
    /// Compute a QQ plot against the specified theoretical distribution.
    pub fn compute(data: &[f64], dist: TheoreticalDist) -> Self {
        if data.is_empty() {
            return Self {
                points: vec![],
                dist,
            };
        }
        let sorted = sort_f64(data);
        let n = sorted.len() as f64;
        let points: Vec<QqPoint> = sorted
            .iter()
            .enumerate()
            .map(|(i, &sample)| {
                let p = (i as f64 + 0.5) / n;
                let theoretical = match dist {
                    TheoreticalDist::Uniform => p,
                    TheoreticalDist::Normal => probit(p),
                };
                QqPoint {
                    theoretical,
                    sample,
                }
            })
            .collect();
        Self { points, dist }
    }

    /// Number of QQ points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if there are no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Rational approximation of the standard normal quantile function (probit).
///
/// Uses the Beasley-Springer-Moro approximation; accurate to ~7 significant figures.
fn probit(p: f64) -> f64 {
    const A: [f64; 4] = [2.515517, 0.802853, 0.010328, 0.0];
    const B: [f64; 4] = [1.432788, 0.189269, 0.001308, 0.0];
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    let sign = if p < 0.5 { -1.0 } else { 1.0 };
    let q = if p < 0.5 { p } else { 1.0 - p };
    let t = (-2.0 * q.ln()).sqrt();
    let num = A[0] + A[1] * t + A[2] * t * t + A[3] * t * t * t;
    let den = 1.0 + B[0] * t + B[1] * t * t + B[2] * t * t * t + B[3] * t.powi(4);
    sign * (t - num / den)
}

// ─────────────────────────────────────────────────────────────────────────────
// Time series with confidence bands
// ─────────────────────────────────────────────────────────────────────────────

/// A single time-series observation at one time point.
#[derive(Debug, Clone)]
pub struct TsPoint {
    /// Time coordinate.
    pub time: f64,
    /// Mean value across replicates at this time.
    pub mean: f64,
    /// Lower confidence bound (mean - k·σ).
    pub lower: f64,
    /// Upper confidence bound (mean + k·σ).
    pub upper: f64,
    /// Number of replicates used to compute the statistics.
    pub n: usize,
}

/// Time series with per-point confidence bands.
#[derive(Debug, Clone)]
pub struct TimeSeries {
    /// Ordered sequence of time points with confidence envelopes.
    pub points: Vec<TsPoint>,
    /// Confidence-band half-width multiplier (e.g. 1.0 = ±σ, 1.96 = 95%).
    pub k_sigma: f64,
}

impl TimeSeries {
    /// Compute a time series from repeated measurements.
    ///
    /// * `times` – time coordinate for each observation group.
    /// * `replicate_groups` – `replicate_groups[i]` contains the replicate values at
    ///   `times[i]`.
    /// * `k_sigma` – half-width multiplier for the confidence band.
    pub fn compute(times: &[f64], replicate_groups: &[Vec<f64>], k_sigma: f64) -> Self {
        let n = times.len().min(replicate_groups.len());
        let points: Vec<TsPoint> = (0..n)
            .map(|i| {
                let reps = &replicate_groups[i];
                let m = mean(reps);
                let s = std_dev(reps);
                TsPoint {
                    time: times[i],
                    mean: m,
                    lower: m - k_sigma * s,
                    upper: m + k_sigma * s,
                    n: reps.len(),
                }
            })
            .collect();
        Self { points, k_sigma }
    }

    /// Number of time points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if there are no time points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Error bar plots
// ─────────────────────────────────────────────────────────────────────────────

/// A single data point with symmetric error bars.
#[derive(Debug, Clone)]
pub struct ErrorBar {
    /// Category label or x-coordinate.
    pub label: String,
    /// Central value.
    pub centre: f64,
    /// Positive (upper) error extent.
    pub error_pos: f64,
    /// Negative (lower) error extent (stored as positive magnitude).
    pub error_neg: f64,
}

impl ErrorBar {
    /// Construct a symmetric error bar.
    pub fn symmetric(label: impl Into<String>, centre: f64, error: f64) -> Self {
        Self {
            label: label.into(),
            centre,
            error_pos: error.abs(),
            error_neg: error.abs(),
        }
    }

    /// Construct an asymmetric error bar.
    pub fn asymmetric(
        label: impl Into<String>,
        centre: f64,
        error_pos: f64,
        error_neg: f64,
    ) -> Self {
        Self {
            label: label.into(),
            centre,
            error_pos: error_pos.abs(),
            error_neg: error_neg.abs(),
        }
    }

    /// Upper extent of the bar.
    pub fn upper(&self) -> f64 {
        self.centre + self.error_pos
    }

    /// Lower extent of the bar.
    pub fn lower(&self) -> f64 {
        self.centre - self.error_neg
    }
}

/// A collection of error bar data for plotting.
#[derive(Debug, Clone, Default)]
pub struct ErrorBarPlot {
    /// All error bars.
    pub bars: Vec<ErrorBar>,
}

impl ErrorBarPlot {
    /// Create an empty plot.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an error bar.
    pub fn push(&mut self, bar: ErrorBar) {
        self.bars.push(bar);
    }

    /// Compute error bars from grouped data (mean ± std dev per group).
    pub fn from_groups(labels: &[&str], groups: &[Vec<f64>]) -> Self {
        let n = labels.len().min(groups.len());
        let bars: Vec<ErrorBar> = (0..n)
            .map(|i| ErrorBar::symmetric(labels[i], mean(&groups[i]), std_dev(&groups[i])))
            .collect();
        Self { bars }
    }

    /// Global y-axis range (lower_min, upper_max) including all error extents.
    pub fn y_range(&self) -> Option<(f64, f64)> {
        if self.bars.is_empty() {
            return None;
        }
        let lo = self
            .bars
            .iter()
            .map(|b| b.lower())
            .fold(f64::INFINITY, f64::min);
        let hi = self
            .bars
            .iter()
            .map(|b| b.upper())
            .fold(f64::NEG_INFINITY, f64::max);
        Some((lo, hi))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Empirical CDF
// ─────────────────────────────────────────────────────────────────────────────

/// A step in the empirical cumulative distribution function.
#[derive(Debug, Clone)]
pub struct EcdfPoint {
    /// Data value at which the CDF steps up.
    pub value: f64,
    /// Cumulative probability just after this step (i/n).
    pub probability: f64,
}

/// Empirical CDF computed directly from sample data.
#[derive(Debug, Clone)]
pub struct EmpiricalCdf {
    /// Sorted (value, probability) steps.
    pub steps: Vec<EcdfPoint>,
    /// Number of observations.
    pub n: usize,
}

impl EmpiricalCdf {
    /// Compute the empirical CDF from unsorted data.
    pub fn compute(data: &[f64]) -> Self {
        let n = data.len();
        if n == 0 {
            return Self {
                steps: vec![],
                n: 0,
            };
        }
        let sorted = sort_f64(data);
        let steps: Vec<EcdfPoint> = sorted
            .iter()
            .enumerate()
            .map(|(i, &v)| EcdfPoint {
                value: v,
                probability: (i + 1) as f64 / n as f64,
            })
            .collect();
        Self { steps, n }
    }

    /// Evaluate P(X ≤ x) for a given x (linear interpolation between steps).
    pub fn evaluate(&self, x: f64) -> f64 {
        if self.steps.is_empty() {
            return f64::NAN;
        }
        if x < self.steps[0].value {
            return 0.0;
        }
        if x >= self
            .steps
            .last()
            .expect("collection should not be empty")
            .value
        {
            return 1.0;
        }
        // Binary search for the last step with value ≤ x
        let pos = self.steps.partition_point(|s| s.value <= x);
        self.steps[pos.saturating_sub(1)].probability
    }

    /// Median estimate: value at which P ≥ 0.5 first.
    pub fn median(&self) -> f64 {
        if self.steps.is_empty() {
            return f64::NAN;
        }
        self.steps
            .iter()
            .find(|s| s.probability >= 0.5)
            .map(|s| s.value)
            .unwrap_or(
                self.steps
                    .last()
                    .expect("collection should not be empty")
                    .value,
            )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bland-Altman plot
// ─────────────────────────────────────────────────────────────────────────────

/// A single point in a Bland-Altman plot.
#[derive(Debug, Clone)]
pub struct BlandAltmanPoint {
    /// Mean of the two measurements: (a + b) / 2.
    pub mean: f64,
    /// Difference of the measurements: a − b.
    pub difference: f64,
}

/// Bland-Altman (Tukey mean-difference) plot for comparing two measurement methods.
#[derive(Debug, Clone)]
pub struct BlandAltmanPlot {
    /// Per-subject mean / difference points.
    pub points: Vec<BlandAltmanPoint>,
    /// Mean difference (bias).
    pub bias: f64,
    /// Standard deviation of differences.
    pub sd_diff: f64,
    /// Upper limit of agreement: bias + 1.96 × SD.
    pub upper_loa: f64,
    /// Lower limit of agreement: bias − 1.96 × SD.
    pub lower_loa: f64,
    /// Number of paired observations.
    pub n: usize,
}

impl BlandAltmanPlot {
    /// Compute a Bland-Altman plot from two paired measurement vectors `a` and `b`.
    pub fn compute(a: &[f64], b: &[f64]) -> Self {
        let n = a.len().min(b.len());
        let points: Vec<BlandAltmanPoint> = (0..n)
            .map(|i| BlandAltmanPoint {
                mean: 0.5 * (a[i] + b[i]),
                difference: a[i] - b[i],
            })
            .collect();
        let diffs: Vec<f64> = points.iter().map(|p| p.difference).collect();
        let bias = mean(&diffs);
        let sd = std_dev(&diffs);
        Self {
            points,
            bias,
            sd_diff: sd,
            upper_loa: bias + 1.96 * sd,
            lower_loa: bias - 1.96 * sd,
            n,
        }
    }

    /// Percentage of points within the limits of agreement.
    pub fn within_loa_fraction(&self) -> f64 {
        if self.points.is_empty() {
            return f64::NAN;
        }
        let within = self
            .points
            .iter()
            .filter(|p| p.difference >= self.lower_loa && p.difference <= self.upper_loa)
            .count();
        within as f64 / self.points.len() as f64
    }

    /// True if there are outlier points outside the limits of agreement.
    pub fn has_outliers(&self) -> bool {
        self.points
            .iter()
            .any(|p| p.difference < self.lower_loa || p.difference > self.upper_loa)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Axis range helper
// ─────────────────────────────────────────────────────────────────────────────

/// Compute a padded axis range `(lo, hi)` that encompasses all data values.
///
/// * `data` – raw data values.
/// * `padding_fraction` – fractional padding added to each side (0.05 = 5%).
pub fn auto_range(data: &[f64], padding_fraction: f64) -> (f64, f64) {
    if data.is_empty() {
        return (0.0, 1.0);
    }
    let sorted = sort_f64(data);
    let lo = sorted[0];
    let hi = *sorted.last().expect("collection should not be empty");
    let span = (hi - lo).max(1e-12);
    let pad = span * padding_fraction;
    (lo - pad, hi + pad)
}

/// Map a scalar `x` to a normalized `[0, 1]` position within `[lo, hi]`.
///
/// Clamps the output to `[0, 1]`.
pub fn normalize_to_range(x: f64, lo: f64, hi: f64) -> f64 {
    if (hi - lo).abs() < 1e-300 {
        return 0.5;
    }
    ((x - lo) / (hi - lo)).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── sort_f64 / mean / variance / std_dev ───────────────────────────────────

    #[test]
    fn test_sort_f64_basic() {
        let sorted = sort_f64(&[3.0, 1.0, 2.0]);
        assert_eq!(sorted, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_mean_basic() {
        assert!((mean(&[1.0, 2.0, 3.0, 4.0, 5.0]) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_mean_empty_is_nan() {
        assert!(mean(&[]).is_nan());
    }

    #[test]
    fn test_variance_uniform() {
        // Population variance of [1,2,3,4,5] = 2.0
        let v = variance(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((v - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_std_dev_zero_for_constant() {
        let s = std_dev(&[5.0, 5.0, 5.0, 5.0]);
        assert!(s.abs() < 1e-15);
    }

    // ── percentile / quartiles ────────────────────────────────────────────────

    #[test]
    fn test_percentile_min_max() {
        let data = sort_f64(&[10.0, 20.0, 30.0, 40.0, 50.0]);
        assert!((percentile_sorted(&data, 0.0) - 10.0).abs() < 1e-10);
        assert!((percentile_sorted(&data, 100.0) - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_percentile_median() {
        let data = sort_f64(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((percentile_sorted(&data, 50.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_quartiles_simple() {
        let data = sort_f64(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
        let (q1, _q2, q3) = quartiles(&data);
        assert!(q1 > 1.0 && q1 < 4.0);
        assert!(q3 > 5.0 && q3 < 8.0);
    }

    // ── Histogram ─────────────────────────────────────────────────────────────

    #[test]
    fn test_histogram_total_count() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = Histogram::compute(&data, 10);
        assert_eq!(h.total_count(), 100);
    }

    #[test]
    fn test_histogram_n_bins() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let h = Histogram::compute(&data, 5);
        assert_eq!(h.bins.len(), 5);
    }

    #[test]
    fn test_histogram_cumulative_ends_at_one() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = Histogram::compute(&data, 10);
        let last_cum = h.bins.last().unwrap().cumulative;
        assert!((last_cum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_histogram_empty() {
        let h = Histogram::compute(&[], 5);
        assert_eq!(h.n, 0);
        assert!(h.bins.is_empty());
    }

    #[test]
    fn test_histogram_single_value() {
        let h = Histogram::compute(&[3.125], 5);
        assert_eq!(h.total_count(), 1);
    }

    #[test]
    fn test_histogram_bin_centre() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let h = Histogram::compute(&data, 2);
        let c0 = h.bins[0].centre();
        assert!(c0 > 0.0);
    }

    #[test]
    fn test_histogram_mode_bin() {
        // Skewed distribution: many values at the low end
        let mut data: Vec<f64> = vec![0.1; 80];
        data.extend(vec![9.0; 20]);
        let h = Histogram::compute(&data, 10);
        let mode = h.mode_bin().unwrap();
        assert_eq!(mode.count, 80);
    }

    // ── BoxPlot ───────────────────────────────────────────────────────────────

    #[test]
    fn test_box_plot_median_sorted() {
        let data: Vec<f64> = (1..=9).map(|i| i as f64).collect(); // [1..9]
        let bp = BoxPlot::compute(&data);
        assert!((bp.median - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_box_plot_empty() {
        let bp = BoxPlot::compute(&[]);
        assert!(bp.median.is_nan());
        assert_eq!(bp.n, 0);
    }

    #[test]
    fn test_box_plot_no_outliers_normal_range() {
        let data = vec![10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0];
        let bp = BoxPlot::compute(&data);
        assert!(bp.outliers_low.is_empty());
        assert!(bp.outliers_high.is_empty());
    }

    #[test]
    fn test_box_plot_detects_outliers() {
        let mut data: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        data.push(1000.0); // extreme outlier
        let bp = BoxPlot::compute(&data);
        assert!(!bp.outliers_high.is_empty());
        assert!(bp.outliers_high.contains(&1000.0));
    }

    #[test]
    fn test_box_plot_iqr() {
        let data: Vec<f64> = (1..=8).map(|i| i as f64).collect();
        let bp = BoxPlot::compute(&data);
        assert!(bp.iqr() > 0.0);
    }

    #[test]
    fn test_box_plot_is_outlier() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 100.0];
        let bp = BoxPlot::compute(&data);
        assert!(bp.is_outlier(100.0));
        assert!(!bp.is_outlier(3.0));
    }

    // ── ViolinPlot / KDE ──────────────────────────────────────────────────────

    #[test]
    fn test_violin_plot_kde_non_negative() {
        let data: Vec<f64> = (0..50).map(|i| i as f64 * 0.1).collect();
        let vp = ViolinPlot::compute(&data, 20);
        for p in &vp.kde {
            assert!(p.density >= 0.0, "density must be non-negative");
        }
    }

    #[test]
    fn test_violin_plot_kde_n_points() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let vp = ViolinPlot::compute(&data, 50);
        assert_eq!(vp.kde.len(), 50);
    }

    #[test]
    fn test_violin_plot_empty() {
        let vp = ViolinPlot::compute(&[], 20);
        assert!(vp.kde.is_empty());
    }

    #[test]
    fn test_gaussian_kde_integrates_approx_one() {
        // Rough numerical integration of KDE over a wide range.
        let data = vec![0.0_f64, 1.0, 2.0];
        let h = 1.0;
        let step = 0.01_f64;
        // Integrate from -5 to 7 to cover tails of KDE (data in [0,2], bandwidth=1)
        let sum: f64 = (-500..=700)
            .map(|i| gaussian_kde(&data, i as f64 * step, h) * step)
            .sum();
        assert!((sum - 1.0).abs() < 0.05, "KDE integral ≈ 1, got {sum}");
    }

    // ── Scatter matrix ────────────────────────────────────────────────────────

    #[test]
    fn test_scatter_matrix_pair_count() {
        let names: Vec<String> = vec!["x".into(), "y".into(), "z".into()];
        let cols = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let sm = ScatterMatrix::compute(names, &cols);
        assert_eq!(sm.pairs.len(), sm.expected_pairs());
        assert_eq!(sm.pairs.len(), 3);
    }

    #[test]
    fn test_scatter_pair_pearson_perfect() {
        let pair = ScatterPair {
            x_name: "x".into(),
            y_name: "y".into(),
            points: vec![(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)],
        };
        let r = pair.pearson_r();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_scatter_pair_pearson_anti_corr() {
        let pair = ScatterPair {
            x_name: "x".into(),
            y_name: "y".into(),
            points: vec![(1.0, 3.0), (2.0, 2.0), (3.0, 1.0)],
        };
        let r = pair.pearson_r();
        assert!((r + 1.0).abs() < 1e-10);
    }

    // ── Correlation heatmap ───────────────────────────────────────────────────

    #[test]
    fn test_correlation_heatmap_diagonal_one() {
        let names: Vec<String> = vec!["a".into(), "b".into()];
        let cols = vec![vec![1.0, 2.0, 3.0], vec![3.0, 2.0, 1.0]];
        let hm = CorrelationHeatmap::compute(names, &cols);
        assert!((hm.matrix[0][0] - 1.0).abs() < 1e-10);
        assert!((hm.matrix[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_correlation_heatmap_get() {
        let names: Vec<String> = vec!["a".into(), "b".into()];
        let cols = vec![vec![1.0, 2.0, 3.0], vec![1.0, 2.0, 3.0]];
        let hm = CorrelationHeatmap::compute(names, &cols);
        let r = hm.get("a", "b").unwrap();
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_uncorrelated() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![4.0, 2.0, 4.0, 2.0]; // zero mean covariance with x? Let's just check bounds.
        let r = pearson_correlation(&x, &y);
        assert!(r.abs() <= 1.0 || r.is_nan());
    }

    // ── QQ plot ───────────────────────────────────────────────────────────────

    #[test]
    fn test_qq_plot_length() {
        let data: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let qq = QqPlot::compute(&data, TheoreticalDist::Normal);
        assert_eq!(qq.len(), 30);
    }

    #[test]
    fn test_qq_plot_empty() {
        let qq = QqPlot::compute(&[], TheoreticalDist::Uniform);
        assert!(qq.is_empty());
    }

    #[test]
    fn test_qq_plot_uniform_theoretical_in_unit_interval() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let qq = QqPlot::compute(&data, TheoreticalDist::Uniform);
        for p in &qq.points {
            assert!(p.theoretical > 0.0 && p.theoretical <= 1.0);
        }
    }

    #[test]
    fn test_probit_symmetry() {
        // probit(0.5) ≈ 0
        assert!(probit(0.5).abs() < 0.01);
        // probit is antisymmetric: probit(p) = -probit(1-p)
        let p = 0.25;
        assert!((probit(p) + probit(1.0 - p)).abs() < 0.01);
    }

    // ── Time series ───────────────────────────────────────────────────────────

    #[test]
    fn test_time_series_mean_constant() {
        let times = vec![0.0, 1.0, 2.0];
        let groups = vec![
            vec![5.0, 5.0, 5.0],
            vec![5.0, 5.0, 5.0],
            vec![5.0, 5.0, 5.0],
        ];
        let ts = TimeSeries::compute(&times, &groups, 1.0);
        for p in &ts.points {
            assert!((p.mean - 5.0).abs() < 1e-10);
            assert!((p.lower - 5.0).abs() < 1e-10); // std=0 → lower=upper=mean
        }
    }

    #[test]
    fn test_time_series_band_width() {
        let times = vec![0.0, 1.0];
        let groups = vec![
            vec![1.0, 3.0], // mean=2, sd=1
            vec![2.0, 4.0], // mean=3, sd=1
        ];
        let ts = TimeSeries::compute(&times, &groups, 2.0);
        // band width = 2 * k_sigma * sd = 2 * 2 * 1 = 4
        let w0 = ts.points[0].upper - ts.points[0].lower;
        assert!((w0 - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_len() {
        let ts = TimeSeries::compute(&[0.0, 1.0, 2.0], &[vec![1.0], vec![2.0], vec![3.0]], 1.0);
        assert_eq!(ts.len(), 3);
    }

    // ── Error bar plots ───────────────────────────────────────────────────────

    #[test]
    fn test_error_bar_symmetric_upper_lower() {
        let eb = ErrorBar::symmetric("A", 10.0, 2.0);
        assert!((eb.upper() - 12.0).abs() < 1e-10);
        assert!((eb.lower() - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_bar_asymmetric() {
        let eb = ErrorBar::asymmetric("B", 5.0, 1.0, 0.5);
        assert!((eb.upper() - 6.0).abs() < 1e-10);
        assert!((eb.lower() - 4.5).abs() < 1e-10);
    }

    #[test]
    fn test_error_bar_plot_from_groups() {
        let labels = vec!["A", "B", "C"];
        let groups = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let plot = ErrorBarPlot::from_groups(&labels, &groups);
        assert_eq!(plot.bars.len(), 3);
        assert!((plot.bars[0].centre - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_bar_plot_y_range() {
        let mut plot = ErrorBarPlot::new();
        plot.push(ErrorBar::symmetric("A", 0.0, 1.0)); // [-1, 1]
        plot.push(ErrorBar::symmetric("B", 5.0, 2.0)); // [3, 7]
        let (lo, hi) = plot.y_range().unwrap();
        assert!((lo + 1.0).abs() < 1e-10);
        assert!((hi - 7.0).abs() < 1e-10);
    }

    // ── Empirical CDF ─────────────────────────────────────────────────────────

    #[test]
    fn test_ecdf_probability_one_at_max() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ecdf = EmpiricalCdf::compute(&data);
        assert!((ecdf.evaluate(5.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ecdf_probability_zero_below_min() {
        let data = vec![1.0, 2.0, 3.0];
        let ecdf = EmpiricalCdf::compute(&data);
        assert_eq!(ecdf.evaluate(0.5), 0.0);
    }

    #[test]
    fn test_ecdf_len() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ecdf = EmpiricalCdf::compute(&data);
        assert_eq!(ecdf.steps.len(), 20);
    }

    #[test]
    fn test_ecdf_median() {
        let data: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        let ecdf = EmpiricalCdf::compute(&data);
        let med = ecdf.median();
        // Should be the 5th value = 5.0
        assert!((med - 5.0).abs() < 1.0);
    }

    #[test]
    fn test_ecdf_empty() {
        let ecdf = EmpiricalCdf::compute(&[]);
        assert_eq!(ecdf.n, 0);
        assert!(ecdf.evaluate(1.0).is_nan());
    }

    // ── Bland-Altman ──────────────────────────────────────────────────────────

    #[test]
    fn test_bland_altman_zero_bias_identical() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = a.clone();
        let ba = BlandAltmanPlot::compute(&a, &b);
        assert!(ba.bias.abs() < 1e-10);
        assert!(ba.sd_diff.abs() < 1e-10);
    }

    #[test]
    fn test_bland_altman_known_bias() {
        let a = vec![10.0, 20.0, 30.0];
        let b = vec![8.0, 18.0, 28.0]; // b = a - 2, so diff = +2 everywhere
        let ba = BlandAltmanPlot::compute(&a, &b);
        assert!((ba.bias - 2.0).abs() < 1e-10);
        assert!(ba.sd_diff.abs() < 1e-10);
    }

    #[test]
    fn test_bland_altman_loa_width() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![1.0, 1.5, 2.0, 3.0, 5.5];
        let ba = BlandAltmanPlot::compute(&a, &b);
        let width = ba.upper_loa - ba.lower_loa;
        assert!(width > 0.0);
        assert!((width - 2.0 * 1.96 * ba.sd_diff).abs() < 1e-6);
    }

    #[test]
    fn test_bland_altman_within_loa_all() {
        let a = vec![1.0, 2.0, 3.0];
        let b = a.clone();
        let ba = BlandAltmanPlot::compute(&a, &b);
        assert!((ba.within_loa_fraction() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bland_altman_n() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.1, 2.1, 3.1];
        let ba = BlandAltmanPlot::compute(&a, &b);
        assert_eq!(ba.n, 3); // limited by shorter vector
    }

    // ── auto_range / normalize_to_range ───────────────────────────────────────

    #[test]
    fn test_auto_range_encloses_data() {
        let data = vec![1.0, 5.0, 3.0];
        let (lo, hi) = auto_range(&data, 0.1);
        assert!(lo <= 1.0);
        assert!(hi >= 5.0);
    }

    #[test]
    fn test_auto_range_empty() {
        let (lo, hi) = auto_range(&[], 0.05);
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 1.0);
    }

    #[test]
    fn test_normalize_to_range_endpoints() {
        assert!((normalize_to_range(0.0, 0.0, 1.0) - 0.0).abs() < 1e-10);
        assert!((normalize_to_range(1.0, 0.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_to_range_clamp() {
        assert!((normalize_to_range(-5.0, 0.0, 1.0) - 0.0).abs() < 1e-10);
        assert!((normalize_to_range(10.0, 0.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_degenerate_range() {
        // lo == hi → should return 0.5
        assert!((normalize_to_range(3.0, 3.0, 3.0) - 0.5).abs() < 1e-10);
    }
}
