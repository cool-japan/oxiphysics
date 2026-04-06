#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uncertainty visualization primitives for physics simulation post-processing.
//!
//! Provides:
//! - [`ConfidenceInterval`] — mean ± CI, percentile bands, bootstrap CI
//! - [`ErrorBand`] — upper/lower bound arrays with shading parameters
//! - [`BoxPlot`] — quartiles, whiskers, outliers, notched variant
//! - [`ViolinPlot`] — Gaussian KDE, Scott/Silverman bandwidth, violin shape
//! - [`SensitivityChart`] — tornado plot, spider plot, Sobol index bar chart
//! - [`EnsemblePlot`] — spaghetti plot, ensemble mean/spread, percentile envelopes
//! - Legacy types: [`UncertaintyBand`], [`ErrorBar`], [`ConfidenceEllipse`],
//!   [`SensitivityPlot`] plus free functions [`fan_chart`], [`reliability_diagram`],
//!   [`calibration_error`], [`ensemble_spread`].

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Internal math helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Standard normal quantile (probit) via Beasley-Springer-Moro rational approximation.
/// `p` must lie in (0, 1).
fn probit(p: f64) -> f64 {
    let a = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    let b = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    let c = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    let d = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];

    let p_low = 0.024_25;
    let p_high = 1.0 - p_low;

    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }

    if p < p_low {
        let q = (-2.0 * p.ln()).sqrt();
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else if p <= p_high {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    }
}

/// Gaussian kernel: K(x; h) = exp(-x²/(2h²)) / (h·√(2π)).
#[inline]
fn gaussian_kernel(x: f64, h: f64) -> f64 {
    let z = x / h;
    (-0.5 * z * z).exp() / (h * (2.0 * PI).sqrt())
}

/// Compute the sample mean and (population) standard deviation of a slice.
fn mean_std(data: &[f64]) -> (f64, f64) {
    if data.is_empty() {
        return (0.0, 0.0);
    }
    let n = data.len() as f64;
    let mu = data.iter().sum::<f64>() / n;
    let var = data.iter().map(|&x| (x - mu).powi(2)).sum::<f64>() / n;
    (mu, var.sqrt())
}

/// Sort a slice and return the value at quantile `p` ∈ \[0, 1\].
fn quantile(data: &[f64], p: f64) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = (p * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

/// Sort a vec in place.
fn sort_vec(v: &mut Vec<f64>) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}

// ─────────────────────────────────────────────────────────────────────────────
// ConfidenceInterval
// ─────────────────────────────────────────────────────────────────────────────

/// A confidence interval computed from a data sample.
///
/// Supports normal-theory CI (`mean ± z·SE`), percentile bands, and
/// bootstrap (percentile) confidence intervals.
#[derive(Debug, Clone)]
pub struct ConfidenceInterval {
    /// Sample mean.
    pub mean: f64,
    /// Lower bound of the confidence interval.
    pub lower: f64,
    /// Upper bound of the confidence interval.
    pub upper: f64,
    /// Coverage level (e.g. `0.95` for a 95% CI).
    pub level: f64,
    /// Original data used to construct this CI.
    pub data: Vec<f64>,
}

impl ConfidenceInterval {
    /// Construct a normal-theory CI for the mean.
    ///
    /// Uses z-score corresponding to `level` (e.g. 0.95 → z ≈ 1.96).
    pub fn normal(data: &[f64], level: f64) -> Self {
        let n = data.len() as f64;
        let (mu, sigma) = mean_std(data);
        let se = if n > 1.0 {
            sigma / (n - 1.0).sqrt()
        } else {
            sigma
        };
        let alpha = 1.0 - level;
        let z = probit(1.0 - alpha / 2.0);
        Self {
            mean: mu,
            lower: mu - z * se,
            upper: mu + z * se,
            level,
            data: data.to_vec(),
        }
    }

    /// Construct a percentile CI from sorted quantiles of the data.
    ///
    /// `lower = quantile(alpha/2)`, `upper = quantile(1 - alpha/2)`.
    pub fn percentile(data: &[f64], level: f64) -> Self {
        let (mu, _) = mean_std(data);
        let alpha = 1.0 - level;
        let lo = quantile(data, alpha / 2.0);
        let hi = quantile(data, 1.0 - alpha / 2.0);
        Self {
            mean: mu,
            lower: lo,
            upper: hi,
            level,
            data: data.to_vec(),
        }
    }

    /// Bootstrap percentile CI using `n_bootstrap` resamples.
    ///
    /// Each resample draws `data.len()` values with replacement and records
    /// its mean. The CI bounds are the `alpha/2` and `1-alpha/2` quantiles
    /// of the bootstrap distribution.
    pub fn bootstrap(data: &[f64], level: f64, n_bootstrap: usize, seed: u64) -> Self {
        let (mu, _) = mean_std(data);
        if data.is_empty() || n_bootstrap == 0 {
            return Self {
                mean: mu,
                lower: mu,
                upper: mu,
                level,
                data: data.to_vec(),
            };
        }
        let n = data.len();
        let mut rng = Lcg64::new(seed);
        let mut boot_means: Vec<f64> = (0..n_bootstrap)
            .map(|_| {
                let sum: f64 = (0..n).map(|_| data[rng.next_usize() % n]).sum();
                sum / n as f64
            })
            .collect();
        sort_vec(&mut boot_means);
        let alpha = 1.0 - level;
        let lo = quantile(&boot_means, alpha / 2.0);
        let hi = quantile(&boot_means, 1.0 - alpha / 2.0);
        Self {
            mean: mu,
            lower: lo,
            upper: hi,
            level,
            data: data.to_vec(),
        }
    }

    /// Width of the confidence interval.
    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }

    /// Returns `true` if the value `v` lies within the CI.
    pub fn contains(&self, v: f64) -> bool {
        v >= self.lower && v <= self.upper
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ErrorBand
// ─────────────────────────────────────────────────────────────────────────────

/// A shaded error band defined by upper and lower bound arrays.
///
/// Suitable for plotting as a filled region around a mean curve.
#[derive(Debug, Clone)]
pub struct ErrorBand {
    /// X positions.
    pub x: Vec<f64>,
    /// Upper bound values.
    pub upper: Vec<f64>,
    /// Lower bound values.
    pub lower: Vec<f64>,
    /// RGBA fill color (each component in \[0, 1\]).
    pub fill_color: [f64; 4],
    /// Opacity of the shaded region (0 = fully transparent, 1 = opaque).
    pub opacity: f64,
    /// Optional label for legend display.
    pub label: Option<String>,
}

impl ErrorBand {
    /// Create a symmetric error band: lower = center − half_width, upper = center + half_width.
    pub fn symmetric(x: Vec<f64>, center: Vec<f64>, half_width: Vec<f64>) -> Self {
        let lower: Vec<f64> = center
            .iter()
            .zip(half_width.iter())
            .map(|(&c, &h)| c - h)
            .collect();
        let upper: Vec<f64> = center
            .iter()
            .zip(half_width.iter())
            .map(|(&c, &h)| c + h)
            .collect();
        Self {
            x,
            upper,
            lower,
            fill_color: [0.2, 0.5, 0.9, 0.3],
            opacity: 0.3,
            label: None,
        }
    }

    /// Create an explicit error band from separate lower and upper arrays.
    pub fn explicit(x: Vec<f64>, lower: Vec<f64>, upper: Vec<f64>) -> Self {
        Self {
            x,
            upper,
            lower,
            fill_color: [0.2, 0.5, 0.9, 0.3],
            opacity: 0.3,
            label: None,
        }
    }

    /// Builder: set the fill color \[r, g, b, a\].
    pub fn with_color(mut self, rgba: [f64; 4]) -> Self {
        self.fill_color = rgba;
        self.opacity = rgba[3];
        self
    }

    /// Builder: set an optional legend label.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Number of points in the band.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Returns `true` if the band has no points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Compute the area between upper and lower bounds (trapezoidal rule).
    pub fn area(&self) -> f64 {
        let n = self.x.len().min(self.upper.len()).min(self.lower.len());
        if n < 2 {
            return 0.0;
        }
        (0..n - 1)
            .map(|i| {
                let dx = self.x[i + 1] - self.x[i];
                let h0 = self.upper[i] - self.lower[i];
                let h1 = self.upper[i + 1] - self.lower[i + 1];
                0.5 * (h0 + h1) * dx.abs()
            })
            .sum()
    }

    /// Return the point-wise midline: (upper + lower) / 2.
    pub fn midline(&self) -> Vec<f64> {
        self.upper
            .iter()
            .zip(self.lower.iter())
            .map(|(&u, &l)| 0.5 * (u + l))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BoxPlot
// ─────────────────────────────────────────────────────────────────────────────

/// Five-number summary and outlier data for a box plot.
///
/// Supports standard and notched box plots.
#[derive(Debug, Clone)]
pub struct BoxPlot {
    /// Group label.
    pub label: String,
    /// Minimum (lower whisker tip).
    pub whisker_low: f64,
    /// First quartile (Q1).
    pub q1: f64,
    /// Median (Q2).
    pub median: f64,
    /// Third quartile (Q3).
    pub q3: f64,
    /// Maximum (upper whisker tip).
    pub whisker_high: f64,
    /// Outlier values (below lower fence or above upper fence).
    pub outliers: Vec<f64>,
    /// Notch half-width: `1.58 * IQR / sqrt(n)`. `None` if not notched.
    pub notch_half_width: Option<f64>,
}

impl BoxPlot {
    /// Build a [`BoxPlot`] from raw data.
    ///
    /// Uses Tukey's fences: `Q1 − 1.5·IQR` and `Q3 + 1.5·IQR`.
    pub fn from_data(label: impl Into<String>, data: &[f64]) -> Self {
        Self::from_data_with_options(label, data, false)
    }

    /// Build a notched [`BoxPlot`] from raw data.
    pub fn notched(label: impl Into<String>, data: &[f64]) -> Self {
        Self::from_data_with_options(label, data, true)
    }

    fn from_data_with_options(label: impl Into<String>, data: &[f64], notched: bool) -> Self {
        let n = data.len();
        let label = label.into();
        if n == 0 {
            return Self {
                label,
                whisker_low: 0.0,
                q1: 0.0,
                median: 0.0,
                q3: 0.0,
                whisker_high: 0.0,
                outliers: vec![],
                notch_half_width: None,
            };
        }
        let q1 = quantile(data, 0.25);
        let med = quantile(data, 0.5);
        let q3 = quantile(data, 0.75);
        let iqr = q3 - q1;
        let fence_low = q1 - 1.5 * iqr;
        let fence_high = q3 + 1.5 * iqr;

        let mut sorted = data.to_vec();
        sort_vec(&mut sorted);

        let whisker_low = sorted
            .iter()
            .cloned()
            .find(|&x| x >= fence_low)
            .unwrap_or(q1);
        let whisker_high = sorted
            .iter()
            .cloned()
            .rev()
            .find(|&x| x <= fence_high)
            .unwrap_or(q3);
        let outliers: Vec<f64> = data
            .iter()
            .cloned()
            .filter(|&x| x < fence_low || x > fence_high)
            .collect();

        let notch_half_width = if notched && n > 0 {
            Some(1.58 * iqr / (n as f64).sqrt())
        } else {
            None
        };

        Self {
            label,
            whisker_low,
            q1,
            median: med,
            q3,
            whisker_high,
            outliers,
            notch_half_width,
        }
    }

    /// Interquartile range: Q3 − Q1.
    pub fn iqr(&self) -> f64 {
        self.q3 - self.q1
    }

    /// Returns `true` if this is a notched box plot.
    pub fn is_notched(&self) -> bool {
        self.notch_half_width.is_some()
    }

    /// Number of outliers.
    pub fn outlier_count(&self) -> usize {
        self.outliers.len()
    }

    /// Notch lower bound: `median − notch_half_width`.
    pub fn notch_lower(&self) -> Option<f64> {
        self.notch_half_width.map(|h| self.median - h)
    }

    /// Notch upper bound: `median + notch_half_width`.
    pub fn notch_upper(&self) -> Option<f64> {
        self.notch_half_width.map(|h| self.median + h)
    }

    /// Five-number summary as an array `[min, Q1, median, Q3, max]`.
    pub fn five_number_summary(&self) -> [f64; 5] {
        [
            self.whisker_low,
            self.q1,
            self.median,
            self.q3,
            self.whisker_high,
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ViolinPlot (expanded)
// ─────────────────────────────────────────────────────────────────────────────

/// Data container for a violin plot with Gaussian KDE and bandwidth selection.
///
/// The violin shape is defined by mirroring the KDE density estimate about
/// the centre axis.
#[derive(Debug, Clone)]
pub struct ViolinPlot {
    /// Group label.
    pub label: String,
    /// Raw data observations.
    pub data: Vec<f64>,
    /// Bandwidth override. If `None`, Scott's rule is used.
    pub bandwidth: Option<f64>,
}

impl ViolinPlot {
    /// Create a new `ViolinPlot` from observations.
    pub fn new(data: Vec<f64>) -> Self {
        Self {
            label: String::new(),
            data,
            bandwidth: None,
        }
    }

    /// Create a labelled violin plot.
    pub fn labelled(label: impl Into<String>, data: Vec<f64>) -> Self {
        Self {
            label: label.into(),
            data,
            bandwidth: None,
        }
    }

    /// Builder: set a fixed bandwidth.
    pub fn with_bandwidth(mut self, h: f64) -> Self {
        self.bandwidth = Some(h);
        self
    }

    /// Scott's rule bandwidth: `h = n^{-1/5} · σ̂`.
    pub fn scott_bandwidth(&self) -> f64 {
        let (_, std) = mean_std(&self.data);
        if std < 1e-14 {
            return 1.0;
        }
        std * (self.data.len() as f64).powf(-0.2)
    }

    /// Silverman's rule bandwidth: `h = 0.9 · min(σ, IQR/1.34) · n^{-1/5}`.
    pub fn silverman_bandwidth(&self) -> f64 {
        let n = self.data.len();
        if n == 0 {
            return 1.0;
        }
        let (_, sigma) = mean_std(&self.data);
        let iqr = quantile(&self.data, 0.75) - quantile(&self.data, 0.25);
        let scale = sigma.min(iqr / 1.34);
        if scale < 1e-14 {
            return 1.0;
        }
        0.9 * scale * (n as f64).powf(-0.2)
    }

    /// Effective bandwidth (override if set, else Scott's rule).
    pub fn effective_bandwidth(&self) -> f64 {
        self.bandwidth.unwrap_or_else(|| self.scott_bandwidth())
    }

    /// Evaluate the Gaussian KDE at `x`.
    pub fn kde(&self, x: f64, h: f64) -> f64 {
        if self.data.is_empty() || h <= 0.0 {
            return 0.0;
        }
        let n = self.data.len() as f64;
        self.data
            .iter()
            .map(|&xi| gaussian_kernel(x - xi, h))
            .sum::<f64>()
            / n
    }

    /// Evaluate the KDE over a grid of `n_grid` points spanning `[lo, hi]`.
    ///
    /// Returns `(grid_x, densities)`.
    pub fn kde_grid(&self, lo: f64, hi: f64, n_grid: usize) -> (Vec<f64>, Vec<f64>) {
        if n_grid == 0 {
            return (vec![], vec![]);
        }
        let h = self.effective_bandwidth();
        let step = (hi - lo) / (n_grid - 1).max(1) as f64;
        let xs: Vec<f64> = (0..n_grid).map(|i| lo + i as f64 * step).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| self.kde(x, h)).collect();
        (xs, ys)
    }

    /// Generate violin polygon points: returns `(left_x, y_positions)` and
    /// `(right_x, y_positions)` for the two sides.
    ///
    /// `center_x` is the x-position of the violin center.
    /// `scale` controls maximum half-width of the violin.
    pub fn violin_shape(
        &self,
        lo: f64,
        hi: f64,
        n_grid: usize,
        center_x: f64,
        scale: f64,
    ) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
        let (ys, densities) = self.kde_grid(lo, hi, n_grid);
        let max_d = densities.iter().cloned().fold(0.0_f64, f64::max).max(1e-14);
        let left: Vec<[f64; 2]> = ys
            .iter()
            .zip(densities.iter())
            .map(|(&y, &d)| [center_x - scale * d / max_d, y])
            .collect();
        let right: Vec<[f64; 2]> = ys
            .iter()
            .zip(densities.iter())
            .map(|(&y, &d)| [center_x + scale * d / max_d, y])
            .collect();
        (left, right)
    }

    /// Median of the data.
    pub fn median(&self) -> f64 {
        quantile(&self.data, 0.5)
    }

    /// Interquartile range: Q3 − Q1.
    pub fn iqr(&self) -> f64 {
        quantile(&self.data, 0.75) - quantile(&self.data, 0.25)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SensitivityChart
// ─────────────────────────────────────────────────────────────────────────────

/// Sobol sensitivity visualization: tornado plot, spider plot, and bar chart.
#[derive(Debug, Clone, Default)]
pub struct SensitivityChart {
    /// Names of the input factors.
    pub factor_names: Vec<String>,
    /// First-order Sobol indices S₁ ∈ \[0, 1\].
    pub s1: Vec<f64>,
    /// Total-order Sobol indices Sₜ ∈ \[0, 1\].
    pub st: Vec<f64>,
    /// Second-order interaction indices S₂ (optional, stored as flat upper triangle).
    pub s2: Vec<f64>,
}

impl SensitivityChart {
    /// Create an empty `SensitivityChart`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a factor with S₁ and Sₜ indices.
    pub fn add_factor(&mut self, name: impl Into<String>, s1: f64, st: f64) {
        self.factor_names.push(name.into());
        self.s1.push(s1.clamp(0.0, 1.0));
        self.st.push(st.clamp(0.0, 1.0));
    }

    /// Number of factors.
    pub fn len(&self) -> usize {
        self.factor_names.len()
    }

    /// Returns `true` if no factors have been added.
    pub fn is_empty(&self) -> bool {
        self.factor_names.is_empty()
    }

    /// Return a tornado plot representation: factors sorted by |S₁| descending,
    /// as `(name, s1, st)` triples.
    pub fn tornado_order(&self) -> Vec<(&str, f64, f64)> {
        let mut triples: Vec<(&str, f64, f64)> = self
            .factor_names
            .iter()
            .map(|s| s.as_str())
            .zip(self.s1.iter().cloned())
            .zip(self.st.iter().cloned())
            .map(|((name, s1), st)| (name, s1, st))
            .collect();
        triples.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        triples
    }

    /// Spider/radar plot vertices at angles evenly distributed over 0..2π.
    ///
    /// Returns `(angles_rad, s1_values)` for drawing a polygon.
    pub fn spider_vertices(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.factor_names.len();
        if n == 0 {
            return (vec![], vec![]);
        }
        let angles: Vec<f64> = (0..n).map(|i| 2.0 * PI * i as f64 / n as f64).collect();
        (angles, self.s1.clone())
    }

    /// Total first-order sensitivity (sum of all S₁).
    pub fn total_s1(&self) -> f64 {
        self.s1.iter().sum()
    }

    /// Total total-order sensitivity (sum of all Sₜ).
    pub fn total_st(&self) -> f64 {
        self.st.iter().sum()
    }

    /// Compute the interaction index for each factor: `Sₜᵢ − S₁ᵢ`.
    pub fn interaction_indices(&self) -> Vec<f64> {
        self.st
            .iter()
            .zip(self.s1.iter())
            .map(|(&t, &f)| (t - f).max(0.0))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EnsemblePlot
// ─────────────────────────────────────────────────────────────────────────────

/// Ensemble data container for spaghetti, envelope, and spread visualization.
///
/// `members[m][t]` is member `m` at time step `t`.
#[derive(Debug, Clone)]
pub struct EnsemblePlot {
    /// Time axis (shared across all members).
    pub times: Vec<f64>,
    /// Ensemble members: each inner `Vec` must have the same length as `times`.
    pub members: Vec<Vec<f64>>,
    /// Optional member labels.
    pub member_labels: Vec<String>,
}

impl EnsemblePlot {
    /// Create an empty ensemble plot.
    pub fn new(times: Vec<f64>) -> Self {
        Self {
            times,
            members: Vec::new(),
            member_labels: Vec::new(),
        }
    }

    /// Add an ensemble member.
    pub fn add_member(&mut self, values: Vec<f64>) {
        self.members.push(values);
        self.member_labels
            .push(format!("member_{}", self.members.len()));
    }

    /// Add a labelled member.
    pub fn add_labelled_member(&mut self, label: impl Into<String>, values: Vec<f64>) {
        self.members.push(values);
        self.member_labels.push(label.into());
    }

    /// Number of ensemble members.
    pub fn n_members(&self) -> usize {
        self.members.len()
    }

    /// Number of time steps.
    pub fn n_times(&self) -> usize {
        self.times.len()
    }

    /// Compute the ensemble mean at each time step.
    pub fn ensemble_mean(&self) -> Vec<f64> {
        let m = self.members.len();
        if m == 0 {
            return vec![0.0; self.times.len()];
        }
        let t = self.times.len();
        (0..t)
            .map(|i| self.members.iter().filter_map(|v| v.get(i)).sum::<f64>() / m as f64)
            .collect()
    }

    /// Compute the ensemble spread (population std) at each time step.
    pub fn ensemble_spread(&self) -> Vec<f64> {
        let m = self.members.len();
        if m == 0 {
            return vec![0.0; self.times.len()];
        }
        let t = self.times.len();
        (0..t)
            .map(|i| {
                let vals: Vec<f64> = self
                    .members
                    .iter()
                    .filter_map(|v| v.get(i).cloned())
                    .collect();
                let (_, std) = mean_std(&vals);
                std
            })
            .collect()
    }

    /// Compute percentile envelopes at each time step.
    ///
    /// Returns `(lower_envelope, upper_envelope)` for the given coverage level.
    /// E.g. `level = 0.9` → 5th and 95th percentile.
    pub fn percentile_envelope(&self, level: f64) -> (Vec<f64>, Vec<f64>) {
        let t = self.times.len();
        let alpha = 1.0 - level;
        let p_lo = alpha / 2.0;
        let p_hi = 1.0 - alpha / 2.0;
        let lo: Vec<f64> = (0..t)
            .map(|i| {
                let vals: Vec<f64> = self
                    .members
                    .iter()
                    .filter_map(|v| v.get(i).cloned())
                    .collect();
                quantile(&vals, p_lo)
            })
            .collect();
        let hi: Vec<f64> = (0..t)
            .map(|i| {
                let vals: Vec<f64> = self
                    .members
                    .iter()
                    .filter_map(|v| v.get(i).cloned())
                    .collect();
                quantile(&vals, p_hi)
            })
            .collect();
        (lo, hi)
    }

    /// Compute multiple percentile envelopes.
    ///
    /// Returns a `Vec` of `(level, lower, upper)` triples.
    pub fn multi_envelope(&self, levels: &[f64]) -> Vec<(f64, Vec<f64>, Vec<f64>)> {
        levels
            .iter()
            .map(|&l| {
                let (lo, hi) = self.percentile_envelope(l);
                (l, lo, hi)
            })
            .collect()
    }

    /// Maximum spread (maximum of `ensemble_spread` across all time steps).
    pub fn max_spread(&self) -> f64 {
        self.ensemble_spread().into_iter().fold(0.0_f64, f64::max)
    }

    /// Spaghetti plot data: returns a reference to all member series.
    pub fn spaghetti(&self) -> &[Vec<f64>] {
        &self.members
    }

    /// Build an [`ErrorBand`] from the ensemble's percentile envelope at `level`.
    pub fn to_error_band(&self, level: f64) -> ErrorBand {
        let (lo, hi) = self.percentile_envelope(level);
        ErrorBand::explicit(self.times.clone(), lo, hi)
            .with_label(format!("{:.0}% envelope", level * 100.0))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Legacy: UncertaintyBand
// ─────────────────────────────────────────────────────────────────────────────

/// A mean curve with associated per-point standard deviations.
///
/// Supports generating confidence bands (σ-based) and percentile bands
/// (quantile-based, assuming Gaussian distribution).
#[derive(Debug, Clone)]
pub struct UncertaintyBand {
    /// X positions (independent variable).
    pub x: Vec<f64>,
    /// Mean values at each x position.
    pub y_mean: Vec<f64>,
    /// Standard deviation at each x position.
    pub y_std: Vec<f64>,
}

impl UncertaintyBand {
    /// Construct a new `UncertaintyBand`.
    ///
    /// All three slices must have the same length.
    pub fn new(x: Vec<f64>, y_mean: Vec<f64>, y_std: Vec<f64>) -> Self {
        Self { x, y_mean, y_std }
    }

    /// Number of points in the band.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Returns `true` if the band has no points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Return `(lower, upper)` confidence band at `±sigma` standard deviations.
    pub fn confidence_band(&self, sigma: f64) -> (Vec<f64>, Vec<f64>) {
        let lower = self
            .y_mean
            .iter()
            .zip(self.y_std.iter())
            .map(|(&m, &s)| m - sigma * s)
            .collect();
        let upper = self
            .y_mean
            .iter()
            .zip(self.y_std.iter())
            .map(|(&m, &s)| m + sigma * s)
            .collect();
        (lower, upper)
    }

    /// Return `(lower, upper)` percentile band assuming Gaussian per-point distributions.
    ///
    /// `p_low` and `p_high` are probabilities in (0, 1).
    pub fn percentile_band(&self, p_low: f64, p_high: f64) -> (Vec<f64>, Vec<f64>) {
        let z_low = probit(p_low.clamp(1e-6, 1.0 - 1e-6));
        let z_high = probit(p_high.clamp(1e-6, 1.0 - 1e-6));
        let lower = self
            .y_mean
            .iter()
            .zip(self.y_std.iter())
            .map(|(&m, &s)| m + z_low * s)
            .collect();
        let upper = self
            .y_mean
            .iter()
            .zip(self.y_std.iter())
            .map(|(&m, &s)| m + z_high * s)
            .collect();
        (lower, upper)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Legacy: ErrorBar (point)
// ─────────────────────────────────────────────────────────────────────────────

/// A single asymmetric error bar data point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrorBar {
    /// X position.
    pub x: f64,
    /// Central Y value.
    pub y: f64,
    /// Lower Y bound.
    pub y_lower: f64,
    /// Upper Y bound.
    pub y_upper: f64,
}

impl ErrorBar {
    /// Create a symmetric error bar: `y_lower = y − err`, `y_upper = y + err`.
    pub fn symmetric(x: f64, y: f64, err: f64) -> Self {
        Self {
            x,
            y,
            y_lower: y - err,
            y_upper: y + err,
        }
    }

    /// Total span of the error bar: `y_upper − y_lower`.
    pub fn span(&self) -> f64 {
        self.y_upper - self.y_lower
    }
}

/// Convert raw `(x, y, y_lower, y_upper)` tuples to a `Vec`ErrorBar`.
pub fn error_bars(data: &[(f64, f64, f64, f64)]) -> Vec<ErrorBar> {
    data.iter()
        .map(|&(x, y, y_lower, y_upper)| ErrorBar {
            x,
            y,
            y_lower,
            y_upper,
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Legacy: ConfidenceEllipse
// ─────────────────────────────────────────────────────────────────────────────

/// A 2-D confidence ellipse defined by a center and covariance matrix.
#[derive(Debug, Clone)]
pub struct ConfidenceEllipse {
    /// Center of the ellipse [x, y].
    pub center: [f64; 2],
    /// 2×2 covariance matrix [[c00, c01\], [c10, c11]].
    pub cov: [[f64; 2]; 2],
}

impl ConfidenceEllipse {
    /// Create a new confidence ellipse.
    pub fn new(center: [f64; 2], cov: [[f64; 2]; 2]) -> Self {
        Self { center, cov }
    }

    /// Compute `n` points on the confidence ellipse for the given `sigma` level.
    pub fn ellipse_points(&self, n: usize, sigma: f64) -> Vec<[f64; 2]> {
        if n == 0 {
            return vec![];
        }
        let a = self.cov[0][0];
        let b = self.cov[0][1];
        let dd = self.cov[1][1];
        let trace = a + dd;
        let det = a * dd - b * b;
        let disc = ((trace / 2.0).powi(2) - det).max(0.0).sqrt();
        let lam1 = trace / 2.0 + disc;
        let lam2 = trace / 2.0 - disc;
        let (e1x, e1y) = if b.abs() > 1e-14 {
            let vx = lam1 - dd;
            let vy = b;
            let norm = (vx * vx + vy * vy).sqrt();
            (vx / norm, vy / norm)
        } else if a >= dd {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };
        let (e2x, e2y) = (-e1y, e1x);
        let r1 = sigma * lam1.max(0.0).sqrt();
        let r2 = sigma * lam2.max(0.0).sqrt();
        (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                let u = r1 * theta.cos();
                let v = r2 * theta.sin();
                [
                    self.center[0] + u * e1x + v * e2x,
                    self.center[1] + u * e1y + v * e2y,
                ]
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Legacy: SensitivityPlot
// ─────────────────────────────────────────────────────────────────────────────

/// Sobol-style first-order sensitivity index visualization data (legacy).
#[derive(Debug, Clone, Default)]
pub struct SensitivityPlot {
    /// Names of the input factors.
    pub factor_names: Vec<String>,
    /// First-order Sobol sensitivity indices S₁ ∈ [0, 1].
    pub s1_indices: Vec<f64>,
    /// Total-order sensitivity indices Sₜ ∈ [0, 1].
    pub st_indices: Vec<f64>,
}

impl SensitivityPlot {
    /// Create an empty `SensitivityPlot`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a factor with its first-order and total-order sensitivity indices.
    pub fn add_factor(&mut self, name: &str, s1: f64, st: f64) {
        self.factor_names.push(name.to_string());
        self.s1_indices.push(s1.clamp(0.0, 1.0));
        self.st_indices.push(st.clamp(0.0, 1.0));
    }

    /// Return pairs `(factor_name, s1_index)` sorted by S₁ descending.
    pub fn ranked_s1(&self) -> Vec<(&str, f64)> {
        let mut pairs: Vec<(&str, f64)> = self
            .factor_names
            .iter()
            .map(|s| s.as_str())
            .zip(self.s1_indices.iter().cloned())
            .collect();
        pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        pairs
    }

    /// Sum of all S₁ indices.
    pub fn total_first_order(&self) -> f64 {
        self.s1_indices.iter().sum()
    }

    /// Number of factors registered.
    pub fn len(&self) -> usize {
        self.factor_names.len()
    }

    /// Returns `true` if no factors have been added.
    pub fn is_empty(&self) -> bool {
        self.factor_names.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions (legacy + new)
// ─────────────────────────────────────────────────────────────────────────────

/// Generate fan chart bands from a mean forecast and standard deviation.
///
/// For each percentile `q` in `percentiles`, returns a vector of `(time, upper_value)` pairs
/// for the upper boundary of that band.
pub fn fan_chart(mean: &[f64], std: &[f64], percentiles: &[f64]) -> Vec<Vec<(f64, f64)>> {
    let n = mean.len().min(std.len());
    percentiles
        .iter()
        .map(|&p| {
            let q = ((1.0 + p) / 2.0).clamp(0.5, 1.0 - 1e-9);
            let z = probit(q);
            (0..n).map(|i| (i as f64, mean[i] + z * std[i])).collect()
        })
        .collect()
}

/// Compute a reliability (calibration) diagram.
///
/// Groups `predicted_probs` into `bins` equal-width bins in \[0, 1\].
pub fn reliability_diagram(
    predicted_probs: &[f64],
    outcomes: &[bool],
    bins: usize,
) -> Vec<(f64, f64)> {
    if predicted_probs.is_empty() || outcomes.is_empty() || bins == 0 {
        return vec![];
    }
    let n = predicted_probs.len().min(outcomes.len());
    let bin_width = 1.0 / bins as f64;
    let mut sum_pred = vec![0.0_f64; bins];
    let mut sum_pos = vec![0.0_f64; bins];
    let mut count = vec![0usize; bins];
    for i in 0..n {
        let p = predicted_probs[i].clamp(0.0, 1.0);
        let bin = ((p / bin_width) as usize).min(bins - 1);
        sum_pred[bin] += p;
        if outcomes[i] {
            sum_pos[bin] += 1.0;
        }
        count[bin] += 1;
    }
    (0..bins)
        .filter_map(|b| {
            if count[b] == 0 {
                None
            } else {
                Some((sum_pred[b] / count[b] as f64, sum_pos[b] / count[b] as f64))
            }
        })
        .collect()
}

/// Expected Calibration Error (ECE).
pub fn calibration_error(predicted_probs: &[f64], outcomes: &[bool], bins: usize) -> f64 {
    if predicted_probs.is_empty() || outcomes.is_empty() || bins == 0 {
        return 0.0;
    }
    let n_total = predicted_probs.len().min(outcomes.len());
    let bin_width = 1.0 / bins as f64;
    let mut sum_pred = vec![0.0_f64; bins];
    let mut sum_pos = vec![0.0_f64; bins];
    let mut count = vec![0usize; bins];
    for i in 0..n_total {
        let p = predicted_probs[i].clamp(0.0, 1.0);
        let bin = ((p / bin_width) as usize).min(bins - 1);
        sum_pred[bin] += p;
        if outcomes[i] {
            sum_pos[bin] += 1.0;
        }
        count[bin] += 1;
    }
    (0..bins)
        .filter(|&b| count[b] > 0)
        .map(|b| {
            let n_b = count[b] as f64;
            let conf = sum_pred[b] / n_b;
            let acc = sum_pos[b] / n_b;
            (n_b / n_total as f64) * (conf - acc).abs()
        })
        .sum()
}

/// Compute ensemble spread (std) across members at each time step.
pub fn ensemble_spread(ensemble: &[Vec<f64>]) -> Vec<f64> {
    if ensemble.is_empty() {
        return vec![];
    }
    let n_members = ensemble.len() as f64;
    let n_steps = ensemble.iter().map(|m| m.len()).min().unwrap_or(0);
    if n_steps == 0 {
        return vec![];
    }
    (0..n_steps)
        .map(|t| {
            let values: Vec<f64> = ensemble.iter().map(|m| m[t]).collect();
            let mu = values.iter().sum::<f64>() / n_members;
            let var = values.iter().map(|&v| (v - mu).powi(2)).sum::<f64>() / n_members;
            var.sqrt()
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple LCG RNG (no external dependency needed for bootstrap)
// ─────────────────────────────────────────────────────────────────────────────

/// A minimal 64-bit linear congruential generator for bootstrap resampling.
struct Lcg64 {
    state: u64,
}

impl Lcg64 {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        // Numerical Recipes LCG constants
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_usize(&mut self) -> usize {
        self.next_u64() as usize
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ConfidenceInterval ────────────────────────────────────────────────────

    #[test]
    fn test_ci_normal_contains_mean() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ci = ConfidenceInterval::normal(&data, 0.95);
        assert!(ci.lower < ci.mean && ci.mean < ci.upper);
    }

    #[test]
    fn test_ci_normal_95_percent_width() {
        let data = vec![0.0_f64; 100];
        // Zero variance data: CI degenerates to a point.
        let ci = ConfidenceInterval::normal(&data, 0.95);
        assert!(ci.width().abs() < 1e-10);
    }

    #[test]
    fn test_ci_percentile_level() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let ci = ConfidenceInterval::percentile(&data, 0.90);
        assert!((ci.level - 0.90).abs() < 1e-10);
        assert!(ci.lower < ci.upper);
    }

    #[test]
    fn test_ci_percentile_bounds_in_data_range() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let ci = ConfidenceInterval::percentile(&data, 0.80);
        let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(ci.lower >= min);
        assert!(ci.upper <= max);
    }

    #[test]
    fn test_ci_bootstrap_width_positive() {
        let data: Vec<f64> = (0..30).map(|i| i as f64 * 0.5).collect();
        let ci = ConfidenceInterval::bootstrap(&data, 0.95, 500, 42);
        assert!(ci.width() > 0.0, "bootstrap CI should have positive width");
    }

    #[test]
    fn test_ci_bootstrap_empty_data() {
        let ci = ConfidenceInterval::bootstrap(&[], 0.95, 100, 1);
        assert_eq!(ci.width(), 0.0);
    }

    #[test]
    fn test_ci_contains() {
        let ci = ConfidenceInterval {
            mean: 5.0,
            lower: 3.0,
            upper: 7.0,
            level: 0.95,
            data: vec![],
        };
        assert!(ci.contains(5.0));
        assert!(!ci.contains(8.0));
        assert!(!ci.contains(2.9));
    }

    #[test]
    fn test_ci_width_symmetric() {
        let data = vec![0.0_f64; 1];
        let ci = ConfidenceInterval::normal(&data, 0.95);
        assert!((ci.upper - ci.mean - (ci.mean - ci.lower)).abs() < 1e-10);
    }

    // ── ErrorBand ─────────────────────────────────────────────────────────────

    #[test]
    fn test_error_band_symmetric() {
        let x = vec![0.0, 1.0, 2.0];
        let center = vec![1.0, 2.0, 3.0];
        let half = vec![0.5, 0.5, 0.5];
        let band = ErrorBand::symmetric(x, center, half);
        assert!((band.upper[0] - 1.5).abs() < 1e-10);
        assert!((band.lower[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_error_band_explicit() {
        let band = ErrorBand::explicit(vec![0.0, 1.0], vec![0.0, 1.0], vec![1.0, 2.0]);
        assert_eq!(band.len(), 2);
        assert!(!band.is_empty());
    }

    #[test]
    fn test_error_band_area_positive() {
        let x = vec![0.0, 1.0, 2.0];
        let band = ErrorBand::explicit(x, vec![0.0; 3], vec![1.0; 3]);
        assert!((band.area() - 2.0).abs() < 1e-10, "area={}", band.area());
    }

    #[test]
    fn test_error_band_midline() {
        let band = ErrorBand::explicit(vec![0.0, 1.0], vec![1.0, 3.0], vec![3.0, 5.0]);
        let mid = band.midline();
        assert!((mid[0] - 2.0).abs() < 1e-10);
        assert!((mid[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_band_with_color() {
        let band =
            ErrorBand::symmetric(vec![0.0], vec![0.0], vec![1.0]).with_color([1.0, 0.0, 0.0, 0.5]);
        assert!((band.fill_color[0] - 1.0).abs() < 1e-10);
        assert!((band.opacity - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_error_band_label() {
        let band = ErrorBand::symmetric(vec![0.0], vec![0.0], vec![1.0]).with_label("95% CI");
        assert_eq!(band.label.as_deref(), Some("95% CI"));
    }

    #[test]
    fn test_error_band_empty() {
        let band = ErrorBand::explicit(vec![], vec![], vec![]);
        assert!(band.is_empty());
        assert_eq!(band.area(), 0.0);
    }

    // ── BoxPlot ───────────────────────────────────────────────────────────────

    #[test]
    fn test_boxplot_from_data_five_numbers() {
        let data: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        let bp = BoxPlot::from_data("test", &data);
        assert!(bp.q1 < bp.median, "Q1 < median");
        assert!(bp.median < bp.q3, "median < Q3");
        assert!(bp.whisker_low <= bp.q1);
        assert!(bp.whisker_high >= bp.q3);
    }

    #[test]
    fn test_boxplot_iqr_positive() {
        let data: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let bp = BoxPlot::from_data("T", &data);
        assert!(bp.iqr() > 0.0);
    }

    #[test]
    fn test_boxplot_outliers_detected() {
        let mut data: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        data.push(1000.0); // extreme outlier
        let bp = BoxPlot::from_data("T", &data);
        assert!(!bp.outliers.is_empty(), "should detect outlier 1000.0");
        assert!(bp.outliers.contains(&1000.0));
    }

    #[test]
    fn test_boxplot_no_outliers_uniform() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let bp = BoxPlot::from_data("T", &data);
        assert!(bp.outliers.is_empty(), "no outliers in uniform data");
    }

    #[test]
    fn test_boxplot_notched() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let bp = BoxPlot::notched("T", &data);
        assert!(bp.is_notched());
        assert!(bp.notch_lower().is_some());
        assert!(bp.notch_upper().is_some());
        assert!(bp.notch_lower().unwrap() < bp.median);
        assert!(bp.notch_upper().unwrap() > bp.median);
    }

    #[test]
    fn test_boxplot_empty_data() {
        let bp = BoxPlot::from_data("T", &[]);
        assert_eq!(bp.median, 0.0);
        assert_eq!(bp.iqr(), 0.0);
    }

    #[test]
    fn test_boxplot_five_number_summary() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let bp = BoxPlot::from_data("T", &data);
        let fns = bp.five_number_summary();
        assert_eq!(fns[2], bp.median);
        assert!(fns[0] <= fns[1]);
        assert!(fns[3] <= fns[4]);
    }

    // ── ViolinPlot ────────────────────────────────────────────────────────────

    #[test]
    fn test_violin_silverman_bandwidth_positive() {
        let vp = ViolinPlot::new(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(vp.silverman_bandwidth() > 0.0);
    }

    #[test]
    fn test_violin_scott_bandwidth_positive() {
        let vp = ViolinPlot::new(vec![0.0, 1.0, 2.0]);
        assert!(vp.scott_bandwidth() > 0.0);
    }

    #[test]
    fn test_violin_bandwidth_override() {
        let vp = ViolinPlot::new(vec![0.0, 1.0]).with_bandwidth(2.5);
        assert!((vp.effective_bandwidth() - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_violin_kde_positive_near_data() {
        let vp = ViolinPlot::new(vec![0.0, 1.0, 2.0]);
        assert!(vp.kde(1.0, 0.5) > 0.0);
    }

    #[test]
    fn test_violin_kde_zero_bandwidth() {
        let vp = ViolinPlot::new(vec![1.0, 2.0, 3.0]);
        assert_eq!(vp.kde(1.5, 0.0), 0.0);
    }

    #[test]
    fn test_violin_kde_grid_length() {
        let vp = ViolinPlot::new(vec![0.0, 1.0, 2.0]);
        let (xs, ys) = vp.kde_grid(0.0, 2.0, 10);
        assert_eq!(xs.len(), 10);
        assert_eq!(ys.len(), 10);
    }

    #[test]
    fn test_violin_shape_left_right() {
        let vp = ViolinPlot::new(vec![0.0, 1.0, 2.0]);
        let (left, right) = vp.violin_shape(0.0, 2.0, 20, 5.0, 1.0);
        assert_eq!(left.len(), 20);
        assert_eq!(right.len(), 20);
        // Right side x > center_x for positive density
        for (l, r) in left.iter().zip(right.iter()) {
            assert!(r[0] >= l[0], "right should be >= left");
        }
    }

    #[test]
    fn test_violin_labelled() {
        let vp = ViolinPlot::labelled("Group A", vec![1.0, 2.0]);
        assert_eq!(vp.label, "Group A");
    }

    #[test]
    fn test_violin_median() {
        let vp = ViolinPlot::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((vp.median() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_violin_iqr_positive() {
        let vp = ViolinPlot::new(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert!(vp.iqr() > 0.0);
    }

    // ── SensitivityChart ──────────────────────────────────────────────────────

    #[test]
    fn test_sensitivity_chart_add_factor() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("x1", 0.4, 0.5);
        sc.add_factor("x2", 0.3, 0.35);
        assert_eq!(sc.len(), 2);
    }

    #[test]
    fn test_sensitivity_chart_tornado_order() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("x1", 0.1, 0.15);
        sc.add_factor("x2", 0.6, 0.7);
        sc.add_factor("x3", 0.3, 0.35);
        let order = sc.tornado_order();
        assert_eq!(order[0].0, "x2", "highest S1 should be first");
        assert_eq!(order[2].0, "x1", "lowest S1 should be last");
    }

    #[test]
    fn test_sensitivity_chart_spider_vertices() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("a", 0.5, 0.6);
        sc.add_factor("b", 0.3, 0.4);
        sc.add_factor("c", 0.2, 0.25);
        let (angles, vals) = sc.spider_vertices();
        assert_eq!(angles.len(), 3);
        assert_eq!(vals.len(), 3);
        assert!((angles[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_sensitivity_chart_total_s1() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("x1", 0.3, 0.4);
        sc.add_factor("x2", 0.4, 0.5);
        assert!((sc.total_s1() - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_sensitivity_chart_interaction_indices() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("x1", 0.3, 0.5);
        sc.add_factor("x2", 0.4, 0.4);
        let inter = sc.interaction_indices();
        assert!(
            (inter[0] - 0.2).abs() < 1e-10,
            "x1 interaction = 0.2, got {}",
            inter[0]
        );
        assert!(
            (inter[1] - 0.0).abs() < 1e-10,
            "x2 interaction = 0, got {}",
            inter[1]
        );
    }

    #[test]
    fn test_sensitivity_chart_empty() {
        let sc = SensitivityChart::new();
        assert!(sc.is_empty());
        assert_eq!(sc.total_s1(), 0.0);
        let (a, v) = sc.spider_vertices();
        assert!(a.is_empty());
        assert!(v.is_empty());
    }

    #[test]
    fn test_sensitivity_chart_clamp() {
        let mut sc = SensitivityChart::new();
        sc.add_factor("x", 1.5, -0.1);
        assert!((sc.s1[0] - 1.0).abs() < 1e-10);
        assert!(sc.st[0].abs() < 1e-10);
    }

    // ── EnsemblePlot ──────────────────────────────────────────────────────────

    #[test]
    fn test_ensemble_mean_single_member() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0]);
        ep.add_member(vec![1.0, 2.0, 3.0]);
        let mean = ep.ensemble_mean();
        assert_eq!(mean, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_ensemble_mean_two_members() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0]);
        ep.add_member(vec![0.0, 0.0]);
        ep.add_member(vec![2.0, 2.0]);
        let mean = ep.ensemble_mean();
        assert!((mean[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ensemble_spread_identical() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0]);
        ep.add_member(vec![5.0, 5.0, 5.0]);
        ep.add_member(vec![5.0, 5.0, 5.0]);
        let spread = ep.ensemble_spread();
        for &s in &spread {
            assert!(s.abs() < 1e-10, "identical members: spread={s}");
        }
    }

    #[test]
    fn test_ensemble_percentile_envelope_ordering() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0]);
        for v in 0..5 {
            ep.add_member(vec![v as f64, v as f64, v as f64]);
        }
        let (lo, hi) = ep.percentile_envelope(0.8);
        for (l, h) in lo.iter().zip(hi.iter()) {
            assert!(l <= h, "lower should be <= upper: {l} vs {h}");
        }
    }

    #[test]
    fn test_ensemble_n_members_and_n_times() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0, 3.0]);
        ep.add_member(vec![1.0; 4]);
        ep.add_member(vec![2.0; 4]);
        assert_eq!(ep.n_members(), 2);
        assert_eq!(ep.n_times(), 4);
    }

    #[test]
    fn test_ensemble_max_spread_positive() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0]);
        ep.add_member(vec![0.0, 0.0, 0.0]);
        ep.add_member(vec![10.0, 10.0, 10.0]);
        assert!(ep.max_spread() > 0.0);
    }

    #[test]
    fn test_ensemble_multi_envelope_levels() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0]);
        for v in 0..10 {
            ep.add_member(vec![v as f64, v as f64]);
        }
        let envelopes = ep.multi_envelope(&[0.5, 0.9]);
        assert_eq!(envelopes.len(), 2);
        // 90% envelope should be wider than 50%.
        let (_, lo50, hi50) = &envelopes[0];
        let (_, lo90, hi90) = &envelopes[1];
        assert!(hi90[0] - lo90[0] >= hi50[0] - lo50[0] - 1e-9);
    }

    #[test]
    fn test_ensemble_to_error_band() {
        let mut ep = EnsemblePlot::new(vec![0.0, 1.0, 2.0]);
        for v in 0..5 {
            ep.add_member(vec![v as f64; 3]);
        }
        let band = ep.to_error_band(0.9);
        assert_eq!(band.len(), 3);
    }

    #[test]
    fn test_ensemble_spaghetti() {
        let mut ep = EnsemblePlot::new(vec![0.0]);
        ep.add_member(vec![1.0]);
        ep.add_member(vec![2.0]);
        let spaghetti = ep.spaghetti();
        assert_eq!(spaghetti.len(), 2);
    }

    #[test]
    fn test_ensemble_labelled_member() {
        let mut ep = EnsemblePlot::new(vec![0.0]);
        ep.add_labelled_member("control", vec![1.0]);
        assert_eq!(ep.member_labels[0], "control");
    }

    // ── Legacy: UncertaintyBand ───────────────────────────────────────────────

    #[test]
    fn test_uncertainty_band_confidence_band_zero_sigma() {
        let band = UncertaintyBand::new(vec![0.0, 1.0], vec![1.0, 2.0], vec![0.5, 0.5]);
        let (lo, hi) = band.confidence_band(0.0);
        assert_eq!(lo, vec![1.0, 2.0]);
        assert_eq!(hi, vec![1.0, 2.0]);
    }

    #[test]
    fn test_uncertainty_band_percentile_95() {
        let band = UncertaintyBand::new(vec![0.0], vec![0.0], vec![1.0]);
        let (lo, hi) = band.percentile_band(0.025, 0.975);
        assert!(lo[0] < -1.9 && lo[0] > -2.1, "lo={}", lo[0]);
        assert!(hi[0] > 1.9 && hi[0] < 2.1, "hi={}", hi[0]);
    }

    #[test]
    fn test_uncertainty_band_len_and_is_empty() {
        let band = UncertaintyBand::new(vec![1.0, 2.0, 3.0], vec![0.0; 3], vec![1.0; 3]);
        assert_eq!(band.len(), 3);
        assert!(!band.is_empty());
    }

    // ── Legacy: ErrorBar ──────────────────────────────────────────────────────

    #[test]
    fn test_error_bar_symmetric() {
        let eb = ErrorBar::symmetric(1.0, 2.0, 0.5);
        assert!((eb.y_lower - 1.5).abs() < 1e-10);
        assert!((eb.y_upper - 2.5).abs() < 1e-10);
        assert!((eb.span() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_bars_conversion() {
        let data = vec![(1.0, 2.0, 1.5, 2.5), (3.0, 4.0, 3.5, 4.5)];
        let bars = error_bars(&data);
        assert_eq!(bars.len(), 2);
        assert_eq!(bars[0].x, 1.0);
    }

    // ── Legacy: ConfidenceEllipse ─────────────────────────────────────────────

    #[test]
    fn test_confidence_ellipse_unit_circle() {
        let ellipse = ConfidenceEllipse::new([0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]]);
        let pts = ellipse.ellipse_points(100, 1.0);
        assert_eq!(pts.len(), 100);
        for p in &pts {
            let r = (p[0].powi(2) + p[1].powi(2)).sqrt();
            assert!((r - 1.0).abs() < 0.01, "radius={r}");
        }
    }

    #[test]
    fn test_confidence_ellipse_zero_points() {
        let ellipse = ConfidenceEllipse::new([0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]]);
        assert!(ellipse.ellipse_points(0, 1.0).is_empty());
    }

    // ── Legacy: SensitivityPlot ───────────────────────────────────────────────

    #[test]
    fn test_sensitivity_plot_add_factor() {
        let mut sp = SensitivityPlot::new();
        sp.add_factor("x1", 0.5, 0.6);
        sp.add_factor("x2", 0.3, 0.4);
        assert_eq!(sp.len(), 2);
        assert!(!sp.is_empty());
    }

    #[test]
    fn test_sensitivity_plot_ranked_s1() {
        let mut sp = SensitivityPlot::new();
        sp.add_factor("x1", 0.2, 0.3);
        sp.add_factor("x2", 0.6, 0.7);
        sp.add_factor("x3", 0.1, 0.2);
        let ranked = sp.ranked_s1();
        assert_eq!(ranked[0].0, "x2");
        assert_eq!(ranked[2].0, "x3");
    }

    // ── Legacy free functions ─────────────────────────────────────────────────

    #[test]
    fn test_fan_chart_basic() {
        let mean = vec![1.0, 2.0, 3.0];
        let std_arr = vec![0.5, 0.5, 0.5];
        let bands = fan_chart(&mean, &std_arr, &[0.5, 0.9]);
        assert_eq!(bands.len(), 2);
        assert_eq!(bands[0].len(), 3);
    }

    #[test]
    fn test_fan_chart_wider_at_higher_percentile() {
        let mean = vec![0.0, 0.0];
        let std_arr = vec![1.0, 1.0];
        let bands = fan_chart(&mean, &std_arr, &[0.5, 0.95]);
        assert!(bands[1][0].1 > bands[0][0].1);
    }

    #[test]
    fn test_reliability_diagram_empty() {
        assert!(reliability_diagram(&[], &[], 10).is_empty());
    }

    #[test]
    fn test_calibration_error_empty() {
        assert_eq!(calibration_error(&[], &[], 10), 0.0);
    }

    #[test]
    fn test_calibration_error_perfect() {
        let preds = vec![0.5; 100];
        let outcomes: Vec<bool> = (0..100).map(|i| i % 2 == 0).collect();
        let ece = calibration_error(&preds, &outcomes, 10);
        assert!(ece < 0.01);
    }

    #[test]
    fn test_ensemble_spread_identical_members() {
        let ensemble = vec![vec![1.0, 2.0, 3.0]; 4];
        let spread = ensemble_spread(&ensemble);
        for &s in &spread {
            assert!(s.abs() < 1e-10);
        }
    }

    #[test]
    fn test_ensemble_spread_empty() {
        assert!(ensemble_spread(&[]).is_empty());
    }

    // ── probit helper ─────────────────────────────────────────────────────────

    #[test]
    fn test_probit_midpoint_zero() {
        assert!(probit(0.5).abs() < 1e-3);
    }

    #[test]
    fn test_probit_symmetry() {
        let z1 = probit(0.025);
        let z2 = probit(0.975);
        assert!((z1 + z2).abs() < 1e-3);
    }

    #[test]
    fn test_probit_monotone() {
        assert!(probit(0.1) < probit(0.9));
    }

    // ── LCG RNG sanity ────────────────────────────────────────────────────────

    #[test]
    fn test_lcg_produces_different_values() {
        let mut rng = Lcg64::new(123);
        let v1 = rng.next_u64();
        let v2 = rng.next_u64();
        assert_ne!(v1, v2);
    }

    #[test]
    fn test_lcg_different_seeds() {
        let mut r1 = Lcg64::new(1);
        let mut r2 = Lcg64::new(2);
        assert_ne!(r1.next_u64(), r2.next_u64());
    }
}
