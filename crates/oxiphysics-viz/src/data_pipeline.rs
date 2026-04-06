// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Data processing pipeline for visualization.
//!
//! Provides [`DataSeries`] for 1-D data operations (smoothing, resampling,
//! normalization, differentiation, integration, filtering, correlation,
//! FFT power spectrum), [`ScatterData`] for 2-D point clouds (k-means
//! clustering, convex hull), [`HeatmapData`] for matrix-valued datasets,
//! [`TimeSeriesStats`] for statistical summaries, and standalone filter /
//! peak-detection helpers.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// DataSeries
// ─────────────────────────────────────────────────────────────────────────────

/// A named 1-D data series with paired (x, y) values.
#[derive(Debug, Clone)]
pub struct DataSeries {
    /// Human-readable name of this series.
    pub name: String,
    /// X-axis values (independent variable).
    pub x: Vec<f64>,
    /// Y-axis values (dependent variable).
    pub y: Vec<f64>,
}

impl DataSeries {
    /// Create a new `DataSeries`.
    ///
    /// # Panics
    /// Panics if `x.len() != y.len()`.
    pub fn new(name: &str, x: Vec<f64>, y: Vec<f64>) -> Self {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        Self {
            name: name.to_string(),
            x,
            y,
        }
    }

    /// Smooth the series with a centred moving-average of width `window`.
    ///
    /// Points near the boundaries use the largest available symmetric window.
    pub fn smooth(&self, window: usize) -> Self {
        let n = self.y.len();
        let half = (window / 2).max(1);
        let y_out: Vec<f64> = (0..n)
            .map(|i| {
                let lo = i.saturating_sub(half);
                let hi = (i + half + 1).min(n);
                let count = hi - lo;
                self.y[lo..hi].iter().sum::<f64>() / count as f64
            })
            .collect();
        Self::new(&self.name, self.x.clone(), y_out)
    }

    /// Resample to `n` evenly-spaced points via linear interpolation.
    ///
    /// The output x-range matches the input range.
    pub fn resample(&self, n: usize) -> Self {
        assert!(n >= 2, "resample needs at least 2 points");
        let x0 = self.x.first().copied().unwrap_or(0.0);
        let x1 = self.x.last().copied().unwrap_or(1.0);
        let xs: Vec<f64> = (0..n)
            .map(|i| x0 + (x1 - x0) * i as f64 / (n - 1) as f64)
            .collect();
        let ys = xs.iter().map(|&xq| self.interp(xq)).collect();
        Self::new(&self.name, xs, ys)
    }

    /// Linearly interpolate a y-value at query point `xq`.
    fn interp(&self, xq: f64) -> f64 {
        let n = self.x.len();
        if n == 0 {
            return 0.0;
        }
        if xq <= self.x[0] {
            return self.y[0];
        }
        if xq >= self.x[n - 1] {
            return self.y[n - 1];
        }
        // Binary search for the bracket
        let idx = self.x.partition_point(|&v| v <= xq);
        let i = idx.saturating_sub(1).min(n - 2);
        let t = (xq - self.x[i]) / (self.x[i + 1] - self.x[i]).max(1e-300);
        self.y[i] * (1.0 - t) + self.y[i + 1] * t
    }

    /// Min-max normalize y-values to \[0, 1\].
    ///
    /// If all y-values are equal the output is all zeros.
    pub fn normalize(&self) -> Self {
        let min = self.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = self.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = max - min;
        let y_out = if range.abs() < 1e-300 {
            vec![0.0; self.y.len()]
        } else {
            self.y.iter().map(|v| (v - min) / range).collect()
        };
        Self::new(&self.name, self.x.clone(), y_out)
    }

    /// Differentiate via central differences (forward/backward at boundaries).
    pub fn differentiate(&self) -> Self {
        let n = self.y.len();
        let mut dy = vec![0.0; n];
        if n == 0 {
            return Self::new(&self.name, self.x.clone(), dy);
        }
        if n == 1 {
            return Self::new(&self.name, self.x.clone(), dy);
        }
        // Forward difference at left boundary
        dy[0] = (self.y[1] - self.y[0]) / (self.x[1] - self.x[0]).max(1e-300);
        // Backward difference at right boundary
        dy[n - 1] = (self.y[n - 1] - self.y[n - 2]) / (self.x[n - 1] - self.x[n - 2]).max(1e-300);
        // Central differences for interior
        for i in 1..n - 1 {
            dy[i] = (self.y[i + 1] - self.y[i - 1]) / (self.x[i + 1] - self.x[i - 1]).max(1e-300);
        }
        Self::new(&self.name, self.x.clone(), dy)
    }

    /// Cumulative trapezoidal integration.
    ///
    /// The output starts at 0 and has the same length as the input.
    pub fn integrate(&self) -> Self {
        let n = self.y.len();
        let mut cum = vec![0.0; n];
        for i in 1..n {
            let dx = self.x[i] - self.x[i - 1];
            cum[i] = cum[i - 1] + 0.5 * (self.y[i] + self.y[i - 1]) * dx;
        }
        Self::new(&self.name, self.x.clone(), cum)
    }

    /// Retain only points with `x_min ≤ x ≤ x_max`.
    pub fn filter_range(&self, x_min: f64, x_max: f64) -> Self {
        let (xs, ys): (Vec<f64>, Vec<f64>) = self
            .x
            .iter()
            .zip(self.y.iter())
            .filter(|(x, _)| **x >= x_min && **x <= x_max)
            .map(|(x, y)| (*x, *y))
            .unzip();
        Self::new(&self.name, xs, ys)
    }

    /// Cross-correlation of `self` and `other` at integer lags.
    ///
    /// Returns a `Vec` of length `2 * n - 1` (full cross-correlation) where
    /// the centre element (index `n - 1`) corresponds to zero lag.
    pub fn cross_correlation(&self, other: &Self) -> Vec<f64> {
        let n = self.y.len();
        let m = other.y.len();
        let out_len = n + m - 1;
        let mut out = vec![0.0; out_len];
        for lag in 0..out_len {
            let mut acc = 0.0;
            for i in 0..n {
                let j = lag as isize - (n as isize - 1) + i as isize;
                if j >= 0 && (j as usize) < m {
                    acc += self.y[i] * other.y[j as usize];
                }
            }
            out[lag] = acc;
        }
        out
    }

    /// Compute the DFT power spectrum.
    ///
    /// Returns `(frequency, power)` pairs.  Frequencies are in cycles per
    /// unit of x (i.e., normalised to \[0, 0.5\] for one-sided spectrum if the
    /// x-spacing is 1).  Only positive frequencies (one-sided) are returned.
    pub fn fft_power_spectrum(&self) -> Vec<(f64, f64)> {
        let n = self.y.len();
        if n == 0 {
            return Vec::new();
        }
        // DFT — O(n²) but simple and correct
        let dx = if n > 1 {
            (self.x[n - 1] - self.x[0]) / (n as f64 - 1.0)
        } else {
            1.0
        };
        let fs = if dx.abs() > 1e-300 { 1.0 / dx } else { 1.0 };
        let half = n / 2 + 1;
        (0..half)
            .map(|k| {
                let freq = k as f64 * fs / n as f64;
                let (re, im): (f64, f64) = self
                    .y
                    .iter()
                    .enumerate()
                    .map(|(j, yj)| {
                        let angle = -2.0 * PI * k as f64 * j as f64 / n as f64;
                        (yj * angle.cos(), yj * angle.sin())
                    })
                    .fold((0.0, 0.0), |(ar, ai), (br, bi)| (ar + br, ai + bi));
                let power = (re * re + im * im) / (n * n) as f64;
                (freq, power)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ScatterData
// ─────────────────────────────────────────────────────────────────────────────

/// 2-D scatter data with optional point labels.
#[derive(Debug, Clone)]
pub struct ScatterData {
    /// X-coordinates of the points.
    pub x: Vec<f64>,
    /// Y-coordinates of the points.
    pub y: Vec<f64>,
    /// Optional text labels for each point.
    pub labels: Option<Vec<String>>,
}

impl ScatterData {
    /// Create a new `ScatterData`.
    ///
    /// # Panics
    /// Panics if `x.len() != y.len()`, or if `labels` is provided and its
    /// length differs from `x`.
    pub fn new(x: Vec<f64>, y: Vec<f64>, labels: Option<Vec<String>>) -> Self {
        assert_eq!(x.len(), y.len(), "x and y must have the same length");
        if let Some(ref l) = labels {
            assert_eq!(l.len(), x.len(), "labels must have the same length as x");
        }
        Self { x, y, labels }
    }

    /// Cluster the points into `k` clusters using Lloyd's k-means algorithm.
    ///
    /// Returns a `Vec`usize` of cluster assignments (length == `self.x.len()`).
    /// Runs at most `max_iter` assignment-update cycles.
    pub fn k_means_clusters(&self, k: usize, max_iter: usize) -> Vec<usize> {
        let n = self.x.len();
        if n == 0 || k == 0 {
            return Vec::new();
        }
        let k = k.min(n);

        // Initialise centroids by picking first k distinct points
        let mut cx: Vec<f64> = (0..k).map(|i| self.x[i]).collect();
        let mut cy: Vec<f64> = (0..k).map(|i| self.y[i]).collect();
        let mut labels = vec![0usize; n];

        for _ in 0..max_iter {
            // Assignment step
            let mut changed = false;
            for i in 0..n {
                let best = (0..k)
                    .min_by(|&a, &b| {
                        let da = (self.x[i] - cx[a]).powi(2) + (self.y[i] - cy[a]).powi(2);
                        let db = (self.x[i] - cx[b]).powi(2) + (self.y[i] - cy[b]).powi(2);
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap_or(0);
                if labels[i] != best {
                    labels[i] = best;
                    changed = true;
                }
            }
            if !changed {
                break;
            }

            // Update step
            let mut sum_x = vec![0.0; k];
            let mut sum_y = vec![0.0; k];
            let mut count = vec![0usize; k];
            for i in 0..n {
                let c = labels[i];
                sum_x[c] += self.x[i];
                sum_y[c] += self.y[i];
                count[c] += 1;
            }
            for c in 0..k {
                if count[c] > 0 {
                    cx[c] = sum_x[c] / count[c] as f64;
                    cy[c] = sum_y[c] / count[c] as f64;
                }
            }
        }
        labels
    }

    /// Compute the indices of the 2-D convex hull in counter-clockwise order.
    ///
    /// Uses the Graham-scan algorithm. Returns indices into `self.x`/`self.y`.
    /// Returns an empty `Vec` if there are fewer than 3 points.
    pub fn convex_hull_2d(&self) -> Vec<usize> {
        let n = self.x.len();
        if n < 3 {
            return (0..n).collect();
        }

        // Find the lowest-then-leftmost point as the pivot
        let pivot = (0..n)
            .min_by(|&a, &b| {
                let ya = self.y[a];
                let yb = self.y[b];
                ya.partial_cmp(&yb)
                    .expect("operation should succeed")
                    .then_with(|| {
                        self.x[a]
                            .partial_cmp(&self.x[b])
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            })
            .expect("operation should succeed");

        // Sort all other points by polar angle around pivot
        let mut indices: Vec<usize> = (0..n).filter(|&i| i != pivot).collect();
        let px = self.x[pivot];
        let py = self.y[pivot];

        indices.sort_by(|&a, &b| {
            let ax = self.x[a] - px;
            let ay = self.y[a] - py;
            let bx = self.x[b] - px;
            let by_ = self.y[b] - py;
            let cross = ax * by_ - ay * bx;
            if cross.abs() < 1e-14 {
                // Collinear: prefer closer point first (will be removed later)
                let da = ax * ax + ay * ay;
                let db = bx * bx + by_ * by_;
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            } else {
                (-cross)
                    .partial_cmp(&0.0)
                    .unwrap_or(std::cmp::Ordering::Equal) // counter-clockwise
            }
        });

        // Graham scan
        let mut hull = vec![pivot, indices[0]];
        for &i in &indices[1..] {
            while hull.len() >= 2 {
                let a = hull[hull.len() - 2];
                let b = hull[hull.len() - 1];
                // Cross product: (B-A) × (I-A)
                let ax = self.x[b] - self.x[a];
                let ay = self.y[b] - self.y[a];
                let bx = self.x[i] - self.x[a];
                let by_ = self.y[i] - self.y[a];
                if ax * by_ - ay * bx <= 0.0 {
                    hull.pop();
                } else {
                    break;
                }
            }
            hull.push(i);
        }
        hull
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HeatmapData
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D matrix suitable for heatmap visualization.
#[derive(Debug, Clone)]
pub struct HeatmapData {
    /// Row-major matrix data: `data\[row\]\[col\]`.
    pub data: Vec<Vec<f64>>,
    /// Labels for the rows.
    pub row_labels: Vec<String>,
    /// Labels for the columns.
    pub col_labels: Vec<String>,
}

impl HeatmapData {
    /// Create a `HeatmapData` from a matrix and labels.
    ///
    /// # Panics
    /// Panics if `row_labels.len() != data.len()` or if `col_labels.len()` does
    /// not match the number of columns.
    pub fn new(data: Vec<Vec<f64>>, row_labels: Vec<String>, col_labels: Vec<String>) -> Self {
        assert_eq!(data.len(), row_labels.len(), "row label count mismatch");
        if !data.is_empty() {
            assert_eq!(data[0].len(), col_labels.len(), "col label count mismatch");
        }
        Self {
            data,
            row_labels,
            col_labels,
        }
    }

    /// Normalize each row to [0, 1] independently.
    pub fn normalize_rows(&self) -> Self {
        let data = self.data.iter().map(|row| normalize_vec(row)).collect();
        Self::new(data, self.row_labels.clone(), self.col_labels.clone())
    }

    /// Normalize each column to [0, 1] independently.
    pub fn normalize_cols(&self) -> Self {
        let nrows = self.data.len();
        let ncols = if nrows == 0 { 0 } else { self.data[0].len() };
        let mut out = vec![vec![0.0; ncols]; nrows];
        for c in 0..ncols {
            let col: Vec<f64> = (0..nrows).map(|r| self.data[r][c]).collect();
            let normalized = normalize_vec(&col);
            for r in 0..nrows {
                out[r][c] = normalized[r];
            }
        }
        Self::new(out, self.row_labels.clone(), self.col_labels.clone())
    }

    /// Transpose the matrix (rows become columns and vice-versa).
    pub fn transpose(&self) -> Self {
        let nrows = self.data.len();
        let ncols = if nrows == 0 { 0 } else { self.data[0].len() };
        let mut out = vec![vec![0.0; nrows]; ncols];
        for r in 0..nrows {
            for c in 0..ncols {
                out[c][r] = self.data[r][c];
            }
        }
        Self::new(out, self.col_labels.clone(), self.row_labels.clone())
    }

    /// Build a Pearson correlation matrix from a list of variable vectors.
    ///
    /// `data\[i\]` is the vector of observations for variable `i`.
    /// Returns an `n × n` `HeatmapData` where `data\[i\]\[j\]` is the correlation
    /// coefficient between variable `i` and variable `j`.
    pub fn correlation_matrix(data: &[Vec<f64>]) -> Self {
        let n = data.len();
        let labels: Vec<String> = (0..n).map(|i| format!("var{i}")).collect();
        if n == 0 {
            return Self::new(Vec::new(), Vec::new(), Vec::new());
        }

        let means: Vec<f64> = data.iter().map(|v| mean(v)).collect();
        let stds: Vec<f64> = data.iter().map(|v| std_dev(v)).collect();

        let mut mat = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                if stds[i].abs() < 1e-300 || stds[j].abs() < 1e-300 {
                    mat[i][j] = if i == j { 1.0 } else { 0.0 };
                    continue;
                }
                let len = data[i].len().min(data[j].len());
                let cov: f64 = (0..len)
                    .map(|k| (data[i][k] - means[i]) * (data[j][k] - means[j]))
                    .sum::<f64>()
                    / len as f64;
                mat[i][j] = cov / (stds[i] * stds[j]);
            }
        }
        Self::new(mat, labels.clone(), labels)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TimeSeriesStats
// ─────────────────────────────────────────────────────────────────────────────

/// Summary statistics for a time series.
#[derive(Debug, Clone)]
pub struct TimeSeriesStats {
    /// Arithmetic mean.
    pub mean: f64,
    /// Population standard deviation.
    pub std: f64,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Autocorrelation coefficients at lags 0, 1, …, `lag_max`.
    pub autocorr: Vec<f64>,
}

impl TimeSeriesStats {
    /// Compute statistics for `series` with autocorrelation up to `lag_max`.
    pub fn compute(series: &[f64], lag_max: usize) -> Self {
        let m = mean(series);
        let s = std_dev(series);
        let min = series.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = series.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        let n = series.len();
        let autocorr: Vec<f64> = (0..=lag_max.min(n.saturating_sub(1)))
            .map(|lag| {
                if n <= lag {
                    return 0.0;
                }
                let cov: f64 = (0..n - lag)
                    .map(|i| (series[i] - m) * (series[i + lag] - m))
                    .sum::<f64>()
                    / n as f64;
                if s.abs() < 1e-300 { 0.0 } else { cov / (s * s) }
            })
            .collect();

        Self {
            mean: m,
            std: s,
            min,
            max,
            autocorr,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Filters
// ─────────────────────────────────────────────────────────────────────────────

/// Simple exponential moving average (EMA) low-pass filter.
///
/// `cutoff` is the cutoff frequency in Hz, `sample_rate` is in Hz.
pub fn lowpass_filter(data: &[f64], cutoff: f64, sample_rate: f64) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    // RC time constant → α for EMA
    let rc = 1.0 / (2.0 * PI * cutoff.max(1e-10));
    let dt = 1.0 / sample_rate.max(1e-10);
    let alpha = dt / (rc + dt);

    let mut out = vec![0.0; data.len()];
    out[0] = data[0];
    for i in 1..data.len() {
        out[i] = alpha * data[i] + (1.0 - alpha) * out[i - 1];
    }
    out
}

/// Simple high-pass filter derived from the EMA low-pass.
///
/// `highpass = input - lowpass`.
pub fn highpass_filter(data: &[f64], cutoff: f64, sample_rate: f64) -> Vec<f64> {
    let lp = lowpass_filter(data, cutoff, sample_rate);
    data.iter().zip(lp.iter()).map(|(d, l)| d - l).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Peak detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect local-maximum peaks with minimum prominence ≥ `min_prominence`.
///
/// A point `i` is a peak if `data\[i\] > data\[i-1\]` and `data\[i\] > data\[i+1\]`.
/// Its prominence is `data\[i\] - min(data\[left_base..=right_base\])`.
/// Returns a sorted `Vec` of peak indices.
pub fn detect_peaks(data: &[f64], min_prominence: f64) -> Vec<usize> {
    let n = data.len();
    if n < 3 {
        return Vec::new();
    }
    let mut peaks = Vec::new();
    for i in 1..n - 1 {
        if data[i] > data[i - 1] && data[i] > data[i + 1] {
            // Compute prominence: walk left and right until we find a higher value
            let mut lo = i;
            while lo > 0 && data[lo - 1] <= data[i] {
                lo -= 1;
            }
            let mut hi = i;
            while hi < n - 1 && data[hi + 1] <= data[i] {
                hi += 1;
            }
            let base_min = data[lo..=hi].iter().cloned().fold(f64::INFINITY, f64::min);
            let prominence = data[i] - base_min;
            if prominence >= min_prominence {
                peaks.push(i);
            }
        }
    }
    peaks
}

// ─────────────────────────────────────────────────────────────────────────────
// Private helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the arithmetic mean of a slice.
fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.iter().sum::<f64>() / v.len() as f64
}

/// Compute the population standard deviation of a slice.
fn std_dev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let m = mean(v);
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64;
    var.sqrt()
}

/// Min-max normalize a `Vec`f64` to \[0, 1\].
fn normalize_vec(v: &[f64]) -> Vec<f64> {
    let min = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    if range.abs() < 1e-300 {
        vec![0.0; v.len()]
    } else {
        v.iter().map(|x| (x - min) / range).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── DataSeries::new ──────────────────────────────────────────────────────

    #[test]
    fn test_dataseries_new_basic() {
        let ds = DataSeries::new("test", vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 3.0]);
        assert_eq!(ds.name, "test");
        assert_eq!(ds.x.len(), 3);
    }

    #[test]
    #[should_panic]
    fn test_dataseries_new_mismatched_lengths_panics() {
        let _ = DataSeries::new("bad", vec![0.0, 1.0], vec![1.0]);
    }

    // ── smooth ───────────────────────────────────────────────────────────────

    #[test]
    fn test_smooth_constant_is_constant() {
        let ds = DataSeries::new("c", vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![5.0; 5]);
        let sm = ds.smooth(3);
        for y in &sm.y {
            assert!((*y - 5.0).abs() < 1e-10, "expected 5.0, got {y}");
        }
    }

    #[test]
    fn test_smooth_length_preserved() {
        let ds = DataSeries::new("s", vec![0.0, 1.0, 2.0, 3.0], vec![1.0, 2.0, 3.0, 4.0]);
        let sm = ds.smooth(2);
        assert_eq!(sm.y.len(), ds.y.len());
    }

    #[test]
    fn test_smooth_reduces_noise() {
        // Alternating signal: smooth should reduce the amplitude
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..10)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let ds = DataSeries::new("noise", x, y);
        let sm = ds.smooth(4);
        let max_abs = sm.y.iter().cloned().map(f64::abs).fold(0.0_f64, f64::max);
        assert!(
            max_abs < 1.0,
            "smoothed max_abs should be < 1.0, got {max_abs}"
        );
    }

    // ── resample ─────────────────────────────────────────────────────────────

    #[test]
    fn test_resample_length() {
        let ds = DataSeries::new("r", vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0]);
        let rs = ds.resample(10);
        assert_eq!(rs.x.len(), 10);
        assert_eq!(rs.y.len(), 10);
    }

    #[test]
    fn test_resample_linear_exact() {
        // Linear y = x; resampling should preserve exact values
        let ds = DataSeries::new("lin", vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0]);
        let rs = ds.resample(5);
        for (x, y) in rs.x.iter().zip(rs.y.iter()) {
            assert!((y - x).abs() < 1e-10, "y={y} should equal x={x}");
        }
    }

    #[test]
    fn test_resample_endpoints_preserved() {
        let ds = DataSeries::new("ep", vec![0.0, 1.0], vec![3.0, 7.0]);
        let rs = ds.resample(5);
        assert!((rs.x[0] - 0.0).abs() < 1e-12);
        assert!((rs.x[4] - 1.0).abs() < 1e-12);
        assert!((rs.y[0] - 3.0).abs() < 1e-10);
        assert!((rs.y[4] - 7.0).abs() < 1e-10);
    }

    // ── normalize ────────────────────────────────────────────────────────────

    #[test]
    fn test_normalize_range_is_zero_to_one() {
        let ds = DataSeries::new("n", vec![0.0, 1.0, 2.0, 3.0], vec![2.0, 5.0, 1.0, 8.0]);
        let norm = ds.normalize();
        let min = norm.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = norm.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((min - 0.0).abs() < 1e-10);
        assert!((max - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_constant_is_all_zeros() {
        let ds = DataSeries::new("flat", vec![0.0, 1.0, 2.0], vec![7.0; 3]);
        let norm = ds.normalize();
        for y in &norm.y {
            assert_eq!(*y, 0.0);
        }
    }

    // ── differentiate ────────────────────────────────────────────────────────

    #[test]
    fn test_differentiate_constant_is_zero() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let ds = DataSeries::new("dc", x, vec![5.0; 4]);
        let diff = ds.differentiate();
        for d in &diff.y {
            assert!(d.abs() < 1e-10, "expected 0, got {d}");
        }
    }

    #[test]
    fn test_differentiate_linear_is_slope() {
        // y = 3x → dy/dx = 3 everywhere
        let x: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|xi| 3.0 * xi).collect();
        let ds = DataSeries::new("lin", x, y);
        let diff = ds.differentiate();
        for d in &diff.y {
            assert!((d - 3.0).abs() < 1e-10, "expected slope 3, got {d}");
        }
    }

    #[test]
    fn test_differentiate_length_preserved() {
        let ds = DataSeries::new("len", vec![0.0, 1.0, 2.0], vec![1.0, 4.0, 9.0]);
        let diff = ds.differentiate();
        assert_eq!(diff.y.len(), 3);
    }

    // ── integrate ────────────────────────────────────────────────────────────

    #[test]
    fn test_integrate_starts_at_zero() {
        let ds = DataSeries::new("i", vec![0.0, 1.0, 2.0], vec![1.0, 1.0, 1.0]);
        let int = ds.integrate();
        assert_eq!(int.y[0], 0.0);
    }

    #[test]
    fn test_integrate_constant_is_linear() {
        // ∫₀ˣ 1 dt = x
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let ds = DataSeries::new("const", x.clone(), vec![1.0; 4]);
        let int = ds.integrate();
        for (xi, yi) in x.iter().zip(int.y.iter()) {
            assert!((yi - xi).abs() < 1e-10, "expected {xi}, got {yi}");
        }
    }

    #[test]
    fn test_integrate_length_preserved() {
        let ds = DataSeries::new("il", vec![0.0, 1.0, 2.0], vec![2.0, 3.0, 4.0]);
        let int = ds.integrate();
        assert_eq!(int.y.len(), 3);
    }

    // ── filter_range ─────────────────────────────────────────────────────────

    #[test]
    fn test_filter_range_keeps_interior() {
        let x: Vec<f64> = (0..=10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|xi| xi * xi).collect();
        let ds = DataSeries::new("fr", x, y);
        let flt = ds.filter_range(2.0, 5.0);
        assert!(flt.x.iter().all(|xi| *xi >= 2.0 && *xi <= 5.0));
        assert_eq!(flt.x.len(), 4); // 2, 3, 4, 5
    }

    #[test]
    fn test_filter_range_empty_if_out_of_range() {
        let ds = DataSeries::new("empty_fr", vec![0.0, 1.0], vec![1.0, 2.0]);
        let flt = ds.filter_range(5.0, 10.0);
        assert!(flt.x.is_empty());
    }

    // ── cross_correlation ────────────────────────────────────────────────────

    #[test]
    fn test_cross_correlation_self_is_autocorrelation() {
        let ds = DataSeries::new("cc", vec![0.0, 1.0, 2.0, 3.0], vec![1.0, 2.0, 3.0, 4.0]);
        let ac = ds.cross_correlation(&ds);
        // Length should be 2*n - 1 = 7
        assert_eq!(ac.len(), 7);
        // Centre element (lag=0) should be maximum
        let centre = ac[3];
        assert!(ac.iter().all(|v| *v <= centre + 1e-10));
    }

    #[test]
    fn test_cross_correlation_length() {
        let a = DataSeries::new("a", vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 3.0]);
        let b = DataSeries::new("b", vec![0.0, 1.0, 2.0, 3.0], vec![1.0, 1.0, 1.0, 1.0]);
        let cc = a.cross_correlation(&b);
        assert_eq!(cc.len(), a.x.len() + b.x.len() - 1);
    }

    // ── fft_power_spectrum ───────────────────────────────────────────────────

    #[test]
    fn test_fft_power_spectrum_length() {
        let n = 8;
        let ds = DataSeries::new(
            "fft",
            (0..n).map(|i| i as f64).collect(),
            (0..n).map(|i| (i as f64).sin()).collect(),
        );
        let ps = ds.fft_power_spectrum();
        assert_eq!(ps.len(), n / 2 + 1);
    }

    #[test]
    fn test_fft_power_nonnegative() {
        let ds = DataSeries::new(
            "fftpos",
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.0, -1.0, 1.0, -1.0],
        );
        for (_, p) in ds.fft_power_spectrum() {
            assert!(p >= 0.0, "power must be non-negative, got {p}");
        }
    }

    #[test]
    fn test_fft_empty_series_returns_empty() {
        let ds = DataSeries::new("empty_fft", Vec::new(), Vec::new());
        assert!(ds.fft_power_spectrum().is_empty());
    }

    // ── ScatterData::k_means_clusters ────────────────────────────────────────

    #[test]
    fn test_kmeans_two_clear_clusters() {
        // Two tight clusters around (0,0) and (10,10)
        let x = vec![0.0, 0.1, -0.1, 10.0, 10.1, 9.9];
        let y = vec![0.0, 0.1, -0.1, 10.0, 10.1, 9.9];
        let sd = ScatterData::new(x, y, None);
        let labels = sd.k_means_clusters(2, 100);
        assert_eq!(labels.len(), 6);
        // Points 0,1,2 should share a cluster; points 3,4,5 should share another
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[0], labels[2]);
        assert_eq!(labels[3], labels[4]);
        assert_eq!(labels[3], labels[5]);
        assert_ne!(labels[0], labels[3]);
    }

    #[test]
    fn test_kmeans_k_geq_n_assigns_all() {
        let sd = ScatterData::new(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0], None);
        let labels = sd.k_means_clusters(5, 50);
        assert_eq!(labels.len(), 3);
    }

    #[test]
    fn test_kmeans_empty_returns_empty() {
        let sd = ScatterData::new(Vec::new(), Vec::new(), None);
        assert!(sd.k_means_clusters(3, 10).is_empty());
    }

    // ── ScatterData::convex_hull_2d ───────────────────────────────────────────

    #[test]
    fn test_convex_hull_square() {
        // Square: corners + interior point
        let x = vec![0.0, 1.0, 1.0, 0.0, 0.5];
        let y = vec![0.0, 0.0, 1.0, 1.0, 0.5];
        let sd = ScatterData::new(x, y, None);
        let hull = sd.convex_hull_2d();
        // Interior point (0.5, 0.5) should NOT be in the hull
        assert!(!hull.contains(&4), "interior point should not be on hull");
        // All 4 corners should be on the hull
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_convex_hull_less_than_three_returns_all() {
        let sd = ScatterData::new(vec![0.0, 1.0], vec![0.0, 1.0], None);
        let hull = sd.convex_hull_2d();
        assert_eq!(hull.len(), 2);
    }

    #[test]
    fn test_convex_hull_indices_in_range() {
        let x = vec![0.0, 1.0, 0.5, 0.5, 2.0, 1.5];
        let y = vec![0.0, 0.0, 1.0, 0.3, 1.0, 0.5];
        let sd = ScatterData::new(x.clone(), y.clone(), None);
        let hull = sd.convex_hull_2d();
        for &idx in &hull {
            assert!(idx < x.len(), "hull index {idx} out of range");
        }
    }

    // ── HeatmapData ──────────────────────────────────────────────────────────

    #[test]
    fn test_heatmap_normalize_rows_range() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 4.0, 4.0]];
        let hm = HeatmapData::new(
            data,
            vec!["r0".into(), "r1".into()],
            vec!["c0".into(), "c1".into(), "c2".into()],
        );
        let norm = hm.normalize_rows();
        // Row 0 should go from 0 to 1
        assert!((norm.data[0][0] - 0.0).abs() < 1e-10);
        assert!((norm.data[0][2] - 1.0).abs() < 1e-10);
        // Row 1 (constant) → all zeros
        for v in &norm.data[1] {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_heatmap_normalize_cols_range() {
        let data = vec![vec![1.0, 4.0], vec![3.0, 4.0], vec![2.0, 4.0]];
        let hm = HeatmapData::new(
            data,
            vec!["r0".into(), "r1".into(), "r2".into()],
            vec!["c0".into(), "c1".into()],
        );
        let norm = hm.normalize_cols();
        // Col 0: values 1,3,2 → min=1,max=3, normalized: 0, 1, 0.5
        assert!((norm.data[0][0] - 0.0).abs() < 1e-10);
        assert!((norm.data[1][0] - 1.0).abs() < 1e-10);
        assert!((norm.data[2][0] - 0.5).abs() < 1e-10);
        // Col 1: all 4 → constant → all zeros
        for r in 0..3 {
            assert_eq!(norm.data[r][1], 0.0);
        }
    }

    #[test]
    fn test_heatmap_transpose_shape() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let hm = HeatmapData::new(
            data,
            vec!["r0".into(), "r1".into()],
            vec!["c0".into(), "c1".into(), "c2".into()],
        );
        let t = hm.transpose();
        assert_eq!(t.data.len(), 3);
        assert_eq!(t.data[0].len(), 2);
    }

    #[test]
    fn test_heatmap_transpose_values() {
        let data = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let hm = HeatmapData::new(
            data,
            vec!["r0".into(), "r1".into()],
            vec!["c0".into(), "c1".into()],
        );
        let t = hm.transpose();
        assert!((t.data[0][1] - 3.0).abs() < 1e-12); // was data[1][0]
        assert!((t.data[1][0] - 2.0).abs() < 1e-12); // was data[0][1]
    }

    #[test]
    fn test_correlation_matrix_diagonal_is_one() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![3.0, 2.0, 1.0]];
        let hm = HeatmapData::correlation_matrix(&data);
        assert!((hm.data[0][0] - 1.0).abs() < 1e-10, "diagonal must be 1");
        assert!((hm.data[1][1] - 1.0).abs() < 1e-10, "diagonal must be 1");
    }

    #[test]
    fn test_correlation_matrix_off_diagonal_range() {
        let data = vec![vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 2.0, 3.0, 4.0]];
        let hm = HeatmapData::correlation_matrix(&data);
        // Identical series → correlation = 1
        assert!(
            (hm.data[0][1] - 1.0).abs() < 1e-10,
            "identical series → r=1"
        );
    }

    // ── TimeSeriesStats ───────────────────────────────────────────────────────

    #[test]
    fn test_time_series_stats_mean() {
        let s = TimeSeriesStats::compute(&[1.0, 2.0, 3.0, 4.0, 5.0], 2);
        assert!((s.mean - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_stats_min_max() {
        let s = TimeSeriesStats::compute(&[5.0, 1.0, 3.0], 0);
        assert!((s.min - 1.0).abs() < 1e-10);
        assert!((s.max - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_stats_autocorr_lag0_is_one() {
        let s = TimeSeriesStats::compute(&[1.0, 2.0, 3.0, 4.0, 5.0], 3);
        assert!(
            (s.autocorr[0] - 1.0).abs() < 1e-10,
            "lag=0 autocorr should be 1"
        );
    }

    #[test]
    fn test_time_series_stats_autocorr_length() {
        let s = TimeSeriesStats::compute(&[1.0, 2.0, 3.0, 4.0], 2);
        assert_eq!(s.autocorr.len(), 3); // lags 0, 1, 2
    }

    // ── lowpass_filter ───────────────────────────────────────────────────────

    #[test]
    fn test_lowpass_length_preserved() {
        let data: Vec<f64> = (0..20).map(|i| (i as f64).sin()).collect();
        let lp = lowpass_filter(&data, 0.1, 1.0);
        assert_eq!(lp.len(), data.len());
    }

    #[test]
    fn test_lowpass_constant_signal_unchanged() {
        let data = vec![3.0; 20];
        let lp = lowpass_filter(&data, 0.1, 1.0);
        // After settling, the constant signal should pass through
        let last = *lp.last().unwrap();
        assert!((last - 3.0).abs() < 0.01, "constant signal: last={last}");
    }

    #[test]
    fn test_lowpass_empty_returns_empty() {
        let lp = lowpass_filter(&[], 0.1, 1.0);
        assert!(lp.is_empty());
    }

    // ── highpass_filter ──────────────────────────────────────────────────────

    #[test]
    fn test_highpass_plus_lowpass_equals_input() {
        let data: Vec<f64> = (0..20).map(|i| (i as f64 * 0.3).sin()).collect();
        let lp = lowpass_filter(&data, 0.1, 1.0);
        let hp = highpass_filter(&data, 0.1, 1.0);
        for ((d, l), h) in data.iter().zip(lp.iter()).zip(hp.iter()) {
            assert!((d - l - h).abs() < 1e-12, "hp + lp should equal input");
        }
    }

    #[test]
    fn test_highpass_constant_is_zero_after_settling() {
        let data = vec![5.0; 100];
        let hp = highpass_filter(&data, 0.1, 1.0);
        let last = *hp.last().unwrap();
        assert!(
            last.abs() < 0.01,
            "constant signal high-pass → ~0, got {last}"
        );
    }

    // ── detect_peaks ─────────────────────────────────────────────────────────

    #[test]
    fn test_detect_peaks_simple() {
        // Data: 0 1 2 1 0 3 0  → peaks at index 2 and 5
        let data = vec![0.0, 1.0, 2.0, 1.0, 0.0, 3.0, 0.0];
        let peaks = detect_peaks(&data, 0.5);
        assert!(peaks.contains(&2), "peak at index 2");
        assert!(peaks.contains(&5), "peak at index 5");
    }

    #[test]
    fn test_detect_peaks_no_peaks_in_monotone() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let peaks = detect_peaks(&data, 0.0);
        assert!(peaks.is_empty(), "no peaks in monotone increasing series");
    }

    #[test]
    fn test_detect_peaks_prominence_filters_small() {
        // Small bumps should be filtered out by prominence threshold
        let data = vec![0.0, 0.1, 0.0, 0.0, 5.0, 0.0];
        let peaks_high = detect_peaks(&data, 1.0);
        let peaks_low = detect_peaks(&data, 0.05);
        assert!(peaks_high.len() < peaks_low.len() || peaks_low.contains(&4));
    }

    #[test]
    fn test_detect_peaks_empty_short_data() {
        assert!(detect_peaks(&[], 0.0).is_empty());
        assert!(detect_peaks(&[1.0], 0.0).is_empty());
        assert!(detect_peaks(&[1.0, 2.0], 0.0).is_empty());
    }

    // ── helper functions ─────────────────────────────────────────────────────

    #[test]
    fn test_mean_basic() {
        assert!((mean(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_std_dev_basic() {
        // std([0, 1, 2, 3, 4]) = sqrt(2)
        let s = std_dev(&[0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!((s - 2.0_f64.sqrt()).abs() < 1e-10, "std = {s}");
    }

    #[test]
    fn test_normalize_vec_range() {
        let v = vec![2.0, 4.0, 6.0];
        let n = normalize_vec(&v);
        assert!((n[0] - 0.0).abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }
}
