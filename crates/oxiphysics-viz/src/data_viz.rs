// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Data visualization utilities: charts, plots, histograms, scatter plots.
//!
//! Pure CPU data structures suitable for downstream rendering. No GPU or
//! windowing dependency is introduced here.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// RGBA color
// ---------------------------------------------------------------------------

/// RGBA color with components in \[0, 1\].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red channel.
    pub r: f64,
    /// Green channel.
    pub g: f64,
    /// Blue channel.
    pub b: f64,
    /// Alpha channel.
    pub a: f64,
}

impl Rgba {
    /// Construct an RGBA colour from four `f64` components.
    pub fn new(r: f64, g: f64, b: f64, a: f64) -> Self {
        Self { r, g, b, a }
    }

    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }

    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }

    /// Opaque red.
    pub fn red() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Opaque green.
    pub fn green() -> Self {
        Self::new(0.0, 1.0, 0.0, 1.0)
    }

    /// Opaque blue.
    pub fn blue() -> Self {
        Self::new(0.0, 0.0, 1.0, 1.0)
    }

    /// Return `true` if all components are within \[0, 1\].
    pub fn is_valid(&self) -> bool {
        let in_range = |v: f64| (0.0..=1.0).contains(&v);
        in_range(self.r) && in_range(self.g) && in_range(self.b) && in_range(self.a)
    }
}

// ---------------------------------------------------------------------------
// ColorMapping — scalar-to-RGBA colormap
// ---------------------------------------------------------------------------

/// Built-in colormap variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColormapKind {
    /// Sequential blue-green-yellow colormap (perceptually uniform).
    Viridis,
    /// Diverging blue-white-red colormap.
    Bwr,
    /// Classic rainbow jet colormap.
    Jet,
    /// Greyscale colormap.
    Greys,
    /// Hot (black-red-yellow-white) colormap.
    Hot,
}

/// Maps a scalar value to an [`Rgba`] colour via a built-in colormap.
#[derive(Debug, Clone)]
pub struct ColorMapping {
    /// The colormap to use.
    pub kind: ColormapKind,
    /// Minimum of the input range.
    pub vmin: f64,
    /// Maximum of the input range.
    pub vmax: f64,
}

impl ColorMapping {
    /// Construct a new `ColorMapping`.
    pub fn new(kind: ColormapKind, vmin: f64, vmax: f64) -> Self {
        Self { kind, vmin, vmax }
    }

    /// Map `value` to an [`Rgba`] colour.  Values outside \[vmin, vmax\] are
    /// clamped.
    pub fn map(&self, value: f64) -> Rgba {
        let t = if (self.vmax - self.vmin).abs() < 1e-15 {
            0.5
        } else {
            ((value - self.vmin) / (self.vmax - self.vmin)).clamp(0.0, 1.0)
        };
        match self.kind {
            ColormapKind::Viridis => viridis(t),
            ColormapKind::Bwr => bwr(t),
            ColormapKind::Jet => jet(t),
            ColormapKind::Greys => Rgba::new(t, t, t, 1.0),
            ColormapKind::Hot => hot(t),
        }
    }
}

fn viridis(t: f64) -> Rgba {
    // Simplified 4-point viridis approximation
    let r = (0.267_f64 + 2.1_f64 * t - 2.5_f64 * t * t + 1.1_f64 * t * t * t).clamp(0.0, 1.0);
    let g = (0.004 + 1.39 * t - 0.56 * t * t).clamp(0.0, 1.0);
    let b = (0.329 + 1.72 * t * (1.0 - t) * 1.5 - 0.8 * t).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}

fn bwr(t: f64) -> Rgba {
    // Blue-White-Red
    if t < 0.5 {
        let s = 2.0 * t;
        Rgba::new(s, s, 1.0, 1.0)
    } else {
        let s = 2.0 * (1.0 - t);
        Rgba::new(1.0, s, s, 1.0)
    }
}

fn jet(t: f64) -> Rgba {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}

fn hot(t: f64) -> Rgba {
    let r = (t * 3.0).clamp(0.0, 1.0);
    let g = ((t * 3.0) - 1.0).clamp(0.0, 1.0);
    let b = ((t * 3.0) - 2.0).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}

// ---------------------------------------------------------------------------
// PlotRange
// ---------------------------------------------------------------------------

/// Axis range, either automatically computed or manually specified.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlotRange {
    /// Auto-compute from the data.
    Auto,
    /// Manually specified \[min, max\].
    Manual(f64, f64),
}

impl PlotRange {
    /// Resolve to concrete `(min, max)` given the data slice.  Returns
    /// `(0.0, 1.0)` if the slice is empty or min == max.
    pub fn resolve(&self, data: &[f64]) -> (f64, f64) {
        match self {
            PlotRange::Manual(lo, hi) => (*lo, *hi),
            PlotRange::Auto => {
                if data.is_empty() {
                    return (0.0, 1.0);
                }
                let lo = data.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                if (hi - lo).abs() < 1e-15 {
                    (lo - 1.0, lo + 1.0)
                } else {
                    (lo, hi)
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// PlotData — time-series line chart
// ---------------------------------------------------------------------------

/// Data for a single line in a time-series plot.
#[derive(Debug, Clone)]
pub struct PlotData {
    /// X axis values.
    pub x: Vec<f64>,
    /// Y axis values (same length as `x`).
    pub y: Vec<f64>,
    /// Series label.
    pub label: String,
    /// Line colour.
    pub color: Rgba,
    /// Axis range for X.
    pub x_range: PlotRange,
    /// Axis range for Y.
    pub y_range: PlotRange,
}

impl PlotData {
    /// Create a new `PlotData` series.
    pub fn new(x: Vec<f64>, y: Vec<f64>, label: impl Into<String>, color: Rgba) -> Self {
        Self {
            x,
            y,
            label: label.into(),
            color,
            x_range: PlotRange::Auto,
            y_range: PlotRange::Auto,
        }
    }

    /// Number of data points.
    pub fn len(&self) -> usize {
        self.x.len().min(self.y.len())
    }

    /// Return `true` if there are no data points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Resolved X range.
    pub fn resolved_x_range(&self) -> (f64, f64) {
        self.x_range.resolve(&self.x)
    }

    /// Resolved Y range.
    pub fn resolved_y_range(&self) -> (f64, f64) {
        self.y_range.resolve(&self.y)
    }
}

// ---------------------------------------------------------------------------
// ScatterPlot
// ---------------------------------------------------------------------------

/// A single point in a scatter plot.
#[derive(Debug, Clone, Copy)]
pub struct ScatterPoint {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Optional Z coordinate for 3-D scatter.
    pub z: Option<f64>,
    /// Point size (pixels).
    pub size: f64,
    /// Point colour.
    pub color: Rgba,
}

impl ScatterPoint {
    /// Create a 2-D scatter point.
    pub fn new_2d(x: f64, y: f64, size: f64, color: Rgba) -> Self {
        Self {
            x,
            y,
            z: None,
            size,
            color,
        }
    }

    /// Create a 3-D scatter point.
    pub fn new_3d(x: f64, y: f64, z: f64, size: f64, color: Rgba) -> Self {
        Self {
            x,
            y,
            z: Some(z),
            size,
            color,
        }
    }
}

/// 2-D / 3-D scatter plot dataset.
#[derive(Debug, Clone)]
pub struct ScatterPlot {
    /// Collection of scatter points.
    pub points: Vec<ScatterPoint>,
    /// Plot label / title.
    pub label: String,
    /// Colour mapping used if colouring by a continuous scalar.
    pub color_mapping: Option<ColorMapping>,
}

impl ScatterPlot {
    /// Create an empty `ScatterPlot`.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            points: Vec::new(),
            label: label.into(),
            color_mapping: None,
        }
    }

    /// Add a 2-D point.
    pub fn push_2d(&mut self, x: f64, y: f64, size: f64, color: Rgba) {
        self.points.push(ScatterPoint::new_2d(x, y, size, color));
    }

    /// Add a 3-D point.
    pub fn push_3d(&mut self, x: f64, y: f64, z: f64, size: f64, color: Rgba) {
        self.points.push(ScatterPoint::new_3d(x, y, z, size, color));
    }

    /// Colour all points from a continuous scalar array via the stored
    /// `color_mapping`.  `values` must have the same length as `points`.
    pub fn apply_color_mapping(&mut self, values: &[f64]) {
        if let Some(ref cm) = self.color_mapping {
            for (pt, &v) in self.points.iter_mut().zip(values.iter()) {
                pt.color = cm.map(v);
            }
        }
    }

    /// Return `true` if all point colours are within the valid \[0,1\] range.
    pub fn all_colors_valid(&self) -> bool {
        self.points.iter().all(|p| p.color.is_valid())
    }
}

// ---------------------------------------------------------------------------
// Histogram bin selection rules
// ---------------------------------------------------------------------------

/// Method for choosing the number of bins in a histogram.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinMethod {
    /// Sturges: k = ceil(log2(n) + 1).
    Sturges,
    /// Scott: h = 3.5 * σ / n^(1/3).
    Scott,
    /// Freedman-Diaconis: h = 2 * IQR / n^(1/3).
    FreedmanDiaconis,
    /// Fixed number of bins supplied by the caller.
    Fixed(usize),
}

/// Normalization mode for a histogram.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistNorm {
    /// Raw bin counts.
    Count,
    /// Probability density (area integrates to 1).
    Pdf,
    /// Cumulative distribution function.
    Cdf,
}

/// A 1-D histogram.
#[derive(Debug, Clone)]
pub struct Histogram {
    /// Left edges of each bin (length = `counts.len()`).
    pub bin_edges: Vec<f64>,
    /// Right edge of the last bin.
    pub bin_right: f64,
    /// Bin heights (after normalization).
    pub counts: Vec<f64>,
    /// Normalization mode used.
    pub norm: HistNorm,
    /// Total number of data points.
    pub total: usize,
}

impl Histogram {
    /// Build a histogram from `data` using the supplied `method` and
    /// `normalization`.
    pub fn from_data(data: &[f64], method: BinMethod, norm: HistNorm) -> Self {
        if data.is_empty() {
            return Self {
                bin_edges: vec![0.0],
                bin_right: 1.0,
                counts: vec![0.0],
                norm,
                total: 0,
            };
        }

        let n = data.len();
        let lo = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = if (hi - lo) < 1e-15 { 1.0 } else { hi - lo };

        // Compute number of bins
        let num_bins: usize = match method {
            BinMethod::Fixed(k) => k.max(1),
            BinMethod::Sturges => ((n as f64).log2().ceil() as usize + 1).max(1),
            BinMethod::Scott => {
                let mean = data.iter().sum::<f64>() / n as f64;
                let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
                let std = var.sqrt().max(1e-15);
                let h = 3.5 * std / (n as f64).cbrt();
                ((range / h).ceil() as usize).max(1)
            }
            BinMethod::FreedmanDiaconis => {
                let mut sorted = data.to_vec();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let q1 = sorted[n / 4];
                let q3 = sorted[(3 * n) / 4];
                let iqr = (q3 - q1).max(1e-15);
                let h = 2.0 * iqr / (n as f64).cbrt();
                ((range / h).ceil() as usize).max(1)
            }
        };

        let width = range / num_bins as f64;
        let mut raw = vec![0usize; num_bins];
        for &v in data {
            let idx = (((v - lo) / width) as usize).min(num_bins - 1);
            raw[idx] += 1;
        }

        let bin_edges: Vec<f64> = (0..num_bins).map(|i| lo + i as f64 * width).collect();
        let bin_right = lo + num_bins as f64 * width;

        let heights = match norm {
            HistNorm::Count => raw.iter().map(|&c| c as f64).collect::<Vec<_>>(),
            HistNorm::Pdf => {
                let denom = n as f64 * width;
                raw.iter().map(|&c| c as f64 / denom).collect()
            }
            HistNorm::Cdf => {
                let mut cum = 0.0;
                raw.iter()
                    .map(|&c| {
                        cum += c as f64 / n as f64;
                        cum
                    })
                    .collect()
            }
        };

        Self {
            bin_edges,
            bin_right,
            counts: heights,
            norm,
            total: n,
        }
    }

    /// Number of bins.
    pub fn num_bins(&self) -> usize {
        self.counts.len()
    }

    /// Bin width (assumes uniform bins).
    pub fn bin_width(&self) -> f64 {
        if self.bin_edges.len() < 2 {
            self.bin_right - self.bin_edges[0]
        } else {
            self.bin_edges[1] - self.bin_edges[0]
        }
    }

    /// Sum of raw counts (for Count normalisation only).
    pub fn count_sum(&self) -> f64 {
        self.counts.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// BarChart
// ---------------------------------------------------------------------------

/// A single bar in a bar chart.
#[derive(Debug, Clone)]
pub struct Bar {
    /// Category label.
    pub label: String,
    /// Bar height / value.
    pub value: f64,
    /// Standard error / half-width of error bar.
    pub error: Option<f64>,
    /// Bar colour.
    pub color: Rgba,
}

impl Bar {
    /// Create a bar without an error bar.
    pub fn new(label: impl Into<String>, value: f64, color: Rgba) -> Self {
        Self {
            label: label.into(),
            value,
            error: None,
            color,
        }
    }

    /// Create a bar with an error bar.
    pub fn with_error(label: impl Into<String>, value: f64, error: f64, color: Rgba) -> Self {
        Self {
            label: label.into(),
            value,
            error: Some(error),
            color,
        }
    }
}

/// A dataset for grouped or stacked bar charts.
#[derive(Debug, Clone)]
pub struct BarChart {
    /// Groups of bars. Outer vec = groups; inner vec = series per group.
    pub groups: Vec<Vec<Bar>>,
    /// Whether bars within a group should be stacked instead of side-by-side.
    pub stacked: bool,
    /// Chart title.
    pub title: String,
}

impl BarChart {
    /// Create an empty bar chart.
    pub fn new(title: impl Into<String>, stacked: bool) -> Self {
        Self {
            groups: Vec::new(),
            title: title.into(),
            stacked,
        }
    }

    /// Add a group of bars.
    pub fn push_group(&mut self, group: Vec<Bar>) {
        self.groups.push(group);
    }

    /// Compute stacked heights for each group (sum of values in each series).
    pub fn stacked_totals(&self) -> Vec<f64> {
        self.groups
            .iter()
            .map(|g| g.iter().map(|b| b.value).sum())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// HeatmapData
// ---------------------------------------------------------------------------

/// Range normalization strategy for heatmap colour mapping.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NormalizeRange {
    /// Normalize using the global minimum and maximum.
    MinMax,
    /// Normalize using given percentiles (e.g., 2nd and 98th).
    Percentile(f64, f64),
    /// Symmetric around zero: range = \[-abs_max, +abs_max\].
    Symmetric,
}

/// A 2-D grid of scalar values suitable for colormapped rendering.
#[derive(Debug, Clone)]
pub struct HeatmapData {
    /// Row-major scalar grid.  `data[row * cols + col]`.
    pub data: Vec<f64>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Colour mapping.
    pub mapping: ColorMapping,
    /// Normalization strategy.
    pub normalize: NormalizeRange,
}

impl HeatmapData {
    /// Construct a heatmap from row-major data.
    pub fn new(
        data: Vec<f64>,
        rows: usize,
        cols: usize,
        mapping: ColorMapping,
        normalize: NormalizeRange,
    ) -> Self {
        Self {
            data,
            rows,
            cols,
            mapping,
            normalize,
        }
    }

    /// Return the value at `(row, col)`.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.cols + col]
    }

    /// Render the heatmap into a flat RGBA byte buffer (4 bytes per pixel,
    /// row-major, top-left first).
    pub fn to_rgba_bytes(&self) -> Vec<u8> {
        let (lo, hi) = self.resolve_range();
        let mut out = Vec::with_capacity(self.rows * self.cols * 4);
        let mut cm = self.mapping.clone();
        cm.vmin = lo;
        cm.vmax = hi;
        for &v in &self.data {
            let c = cm.map(v);
            out.push((c.r * 255.0).round() as u8);
            out.push((c.g * 255.0).round() as u8);
            out.push((c.b * 255.0).round() as u8);
            out.push((c.a * 255.0).round() as u8);
        }
        out
    }

    fn resolve_range(&self) -> (f64, f64) {
        match self.normalize {
            NormalizeRange::MinMax => {
                let lo = self.data.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                (lo, hi)
            }
            NormalizeRange::Percentile(p_lo, p_hi) => {
                let mut sorted = self.data.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let n = sorted.len();
                let idx_lo = ((p_lo / 100.0) * (n - 1) as f64).round() as usize;
                let idx_hi = ((p_hi / 100.0) * (n - 1) as f64).round() as usize;
                (sorted[idx_lo], sorted[idx_hi])
            }
            NormalizeRange::Symmetric => {
                let abs_max = self.data.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
                (-abs_max, abs_max)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ContourPlot — Marching Squares
// ---------------------------------------------------------------------------

/// A single iso-contour line segment (a pair of 2-D points).
#[derive(Debug, Clone, Copy)]
pub struct ContourSegment {
    /// Start point `[x, y]`.
    pub p0: [f64; 2],
    /// End point `[x, y]`.
    pub p1: [f64; 2],
}

/// 2-D contour plot: extracts iso-contour lines from a scalar grid using the
/// Marching Squares algorithm.
#[derive(Debug, Clone)]
pub struct ContourPlot {
    /// Scalar grid (row-major, `data[row * cols + col]`).
    pub data: Vec<f64>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Physical extent: `[x_min, x_max, y_min, y_max]`.
    pub extent: [f64; 4],
}

impl ContourPlot {
    /// Construct a `ContourPlot`.
    pub fn new(data: Vec<f64>, rows: usize, cols: usize, extent: [f64; 4]) -> Self {
        Self {
            data,
            rows,
            cols,
            extent,
        }
    }

    fn val(&self, r: usize, c: usize) -> f64 {
        self.data[r * self.cols + c]
    }

    fn to_xy(&self, r: f64, c: f64) -> [f64; 2] {
        let [x0, x1, y0, y1] = self.extent;
        let x = x0 + c / (self.cols - 1).max(1) as f64 * (x1 - x0);
        let y = y0 + r / (self.rows - 1).max(1) as f64 * (y1 - y0);
        [x, y]
    }

    fn interp(v0: f64, v1: f64, iso: f64) -> f64 {
        if (v1 - v0).abs() < 1e-15 {
            0.5
        } else {
            (iso - v0) / (v1 - v0)
        }
    }

    /// Extract all contour segments at the given iso-value using Marching
    /// Squares.
    pub fn extract(&self, iso: f64) -> Vec<ContourSegment> {
        let mut segs = Vec::new();
        for r in 0..self.rows.saturating_sub(1) {
            for c in 0..self.cols.saturating_sub(1) {
                let bl = self.val(r, c);
                let br = self.val(r, c + 1);
                let tl = self.val(r + 1, c);
                let tr = self.val(r + 1, c + 1);

                let idx = ((bl >= iso) as u8)
                    | (((br >= iso) as u8) << 1)
                    | (((tr >= iso) as u8) << 2)
                    | (((tl >= iso) as u8) << 3);

                // Edge midpoints (fractional row/col)
                let left = || {
                    let t = Self::interp(bl, tl, iso);
                    self.to_xy(r as f64 + t, c as f64)
                };
                let right = || {
                    let t = Self::interp(br, tr, iso);
                    self.to_xy(r as f64 + t, (c + 1) as f64)
                };
                let bottom = || {
                    let t = Self::interp(bl, br, iso);
                    self.to_xy(r as f64, c as f64 + t)
                };
                let top = || {
                    let t = Self::interp(tl, tr, iso);
                    self.to_xy((r + 1) as f64, c as f64 + t)
                };

                match idx {
                    0 | 15 => {}
                    1 | 14 => segs.push(ContourSegment {
                        p0: left(),
                        p1: bottom(),
                    }),
                    2 | 13 => segs.push(ContourSegment {
                        p0: bottom(),
                        p1: right(),
                    }),
                    3 | 12 => segs.push(ContourSegment {
                        p0: left(),
                        p1: right(),
                    }),
                    4 | 11 => segs.push(ContourSegment {
                        p0: right(),
                        p1: top(),
                    }),
                    5 => {
                        segs.push(ContourSegment {
                            p0: left(),
                            p1: top(),
                        });
                        segs.push(ContourSegment {
                            p0: bottom(),
                            p1: right(),
                        });
                    }
                    6 | 9 => segs.push(ContourSegment {
                        p0: bottom(),
                        p1: top(),
                    }),
                    7 | 8 => segs.push(ContourSegment {
                        p0: left(),
                        p1: top(),
                    }),
                    10 => {
                        segs.push(ContourSegment {
                            p0: left(),
                            p1: bottom(),
                        });
                        segs.push(ContourSegment {
                            p0: right(),
                            p1: top(),
                        });
                    }
                    _ => {}
                }
            }
        }
        segs
    }
}

// ---------------------------------------------------------------------------
// VectorField2d
// ---------------------------------------------------------------------------

/// A single arrow in a 2-D vector field visualization.
#[derive(Debug, Clone, Copy)]
pub struct Arrow2d {
    /// Arrow origin `[x, y]`.
    pub origin: [f64; 2],
    /// Arrow direction (not normalized) `[dx, dy]`.
    pub direction: [f64; 2],
    /// Display colour.
    pub color: Rgba,
}

/// 2-D vector field visualization (arrow plot).
#[derive(Debug, Clone)]
pub struct VectorField2d {
    /// Collection of arrows.
    pub arrows: Vec<Arrow2d>,
    /// Scale factor applied to arrow lengths for display.
    pub scale: f64,
}

impl VectorField2d {
    /// Create an empty `VectorField2d` with the given display scale.
    pub fn new(scale: f64) -> Self {
        Self {
            arrows: Vec::new(),
            scale,
        }
    }

    /// Sample a vector field `f(x,y) -> [dx,dy]` on a regular grid.
    pub fn from_grid<F>(
        x_range: (f64, f64),
        y_range: (f64, f64),
        nx: usize,
        ny: usize,
        scale: f64,
        cmap: Option<&ColorMapping>,
        f: F,
    ) -> Self
    where
        F: Fn(f64, f64) -> [f64; 2],
    {
        let mut arrows = Vec::with_capacity(nx * ny);
        let mut mags: Vec<f64> = Vec::new();
        let mut raw_arrows: Vec<([f64; 2], [f64; 2])> = Vec::new();
        for iy in 0..ny {
            for ix in 0..nx {
                let x = x_range.0 + ix as f64 / (nx - 1).max(1) as f64 * (x_range.1 - x_range.0);
                let y = y_range.0 + iy as f64 / (ny - 1).max(1) as f64 * (y_range.1 - y_range.0);
                let d = f(x, y);
                let mag = (d[0] * d[0] + d[1] * d[1]).sqrt();
                mags.push(mag);
                raw_arrows.push(([x, y], d));
            }
        }
        let max_mag = mags.iter().cloned().fold(0.0_f64, f64::max).max(1e-15);
        for (i, (origin, direction)) in raw_arrows.into_iter().enumerate() {
            let color = cmap.map(|cm| cm.map(mags[i])).unwrap_or(Rgba::new(
                mags[i] / max_mag,
                0.2,
                1.0 - mags[i] / max_mag,
                1.0,
            ));
            arrows.push(Arrow2d {
                origin,
                direction,
                color,
            });
        }
        Self { arrows, scale }
    }
}

// ---------------------------------------------------------------------------
// ParallelCoordinates
// ---------------------------------------------------------------------------

/// A single observation (row) in a parallel-coordinates plot.
#[derive(Debug, Clone)]
pub struct PcRow {
    /// Values along each axis (one per axis).
    pub values: Vec<f64>,
    /// Line colour for this observation.
    pub color: Rgba,
    /// Optional label / class name.
    pub label: Option<String>,
}

/// Parallel-coordinates plot for multivariate data.
#[derive(Debug, Clone)]
pub struct ParallelCoordinates {
    /// Axis names.
    pub axes: Vec<String>,
    /// Observations.
    pub rows: Vec<PcRow>,
    /// Whether to normalize each axis independently to \[0, 1\].
    pub normalize: bool,
}

impl ParallelCoordinates {
    /// Create an empty `ParallelCoordinates` with given axis names.
    pub fn new(axes: Vec<String>, normalize: bool) -> Self {
        Self {
            axes,
            rows: Vec::new(),
            normalize,
        }
    }

    /// Add an observation.
    pub fn push(&mut self, values: Vec<f64>, color: Rgba, label: Option<String>) {
        self.rows.push(PcRow {
            values,
            color,
            label,
        });
    }

    /// Return normalized values for a row in \[0, 1\] per axis.
    pub fn normalized_row(&self, row_idx: usize) -> Vec<f64> {
        let row = &self.rows[row_idx];
        (0..self.axes.len())
            .map(|ax| {
                let col: Vec<f64> = self.rows.iter().map(|r| r.values[ax]).collect();
                let lo = col.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = col.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                if (hi - lo) < 1e-15 {
                    0.5
                } else {
                    (row.values[ax] - lo) / (hi - lo)
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// RadarChart
// ---------------------------------------------------------------------------

/// A single series in a radar / spider chart.
#[derive(Debug, Clone)]
pub struct RadarSeries {
    /// Series name.
    pub name: String,
    /// Values for each axis (same length as `RadarChart::axes`).
    pub values: Vec<f64>,
    /// Display colour.
    pub color: Rgba,
}

/// Radar (spider) chart for multivariate comparison.
#[derive(Debug, Clone)]
pub struct RadarChart {
    /// Axis names (one per spoke).
    pub axes: Vec<String>,
    /// Data series.
    pub series: Vec<RadarSeries>,
    /// Common scale maximum (all axes share the same max for fair comparison).
    pub scale_max: f64,
}

impl RadarChart {
    /// Create a new `RadarChart`.
    pub fn new(axes: Vec<String>, scale_max: f64) -> Self {
        Self {
            axes,
            series: Vec::new(),
            scale_max,
        }
    }

    /// Add a data series.
    pub fn push(&mut self, name: impl Into<String>, values: Vec<f64>, color: Rgba) {
        self.series.push(RadarSeries {
            name: name.into(),
            values,
            color,
        });
    }

    /// Compute the `[x, y]` polygon points for a series (unit circle layout).
    pub fn polygon_points(&self, series_idx: usize) -> Vec<[f64; 2]> {
        let n = self.axes.len();
        let s = &self.series[series_idx];
        (0..n)
            .map(|i| {
                let angle = 2.0 * PI * i as f64 / n as f64 - PI / 2.0;
                let r = (s.values[i] / self.scale_max.max(1e-15)).clamp(0.0, 1.0);
                [r * angle.cos(), r * angle.sin()]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// TimeSeriesPlot
// ---------------------------------------------------------------------------

/// A marker annotation on a time-series plot.
#[derive(Debug, Clone)]
pub struct TimeMarker {
    /// Time (X axis) position.
    pub t: f64,
    /// Annotation text.
    pub label: String,
    /// Marker colour.
    pub color: Rgba,
}

/// Time-indexed plot with support for multiple channels, markers, and zoom.
#[derive(Debug, Clone)]
pub struct TimeSeriesPlot {
    /// Named channels: `(channel_name, timestamps, values)`.
    pub channels: Vec<(String, Vec<f64>, Vec<f64>)>,
    /// Annotated time markers.
    pub markers: Vec<TimeMarker>,
    /// Visible time window `[t_min, t_max]` (None = auto).
    pub time_window: Option<(f64, f64)>,
    /// Visible value window `[v_min, v_max]` (None = auto).
    pub value_window: Option<(f64, f64)>,
}

impl TimeSeriesPlot {
    /// Create an empty `TimeSeriesPlot`.
    pub fn new() -> Self {
        Self {
            channels: Vec::new(),
            markers: Vec::new(),
            time_window: None,
            value_window: None,
        }
    }

    /// Add a channel.
    pub fn add_channel(&mut self, name: impl Into<String>, t: Vec<f64>, v: Vec<f64>) {
        self.channels.push((name.into(), t, v));
    }

    /// Add a time marker.
    pub fn add_marker(&mut self, t: f64, label: impl Into<String>, color: Rgba) {
        self.markers.push(TimeMarker {
            t,
            label: label.into(),
            color,
        });
    }

    /// Set the visible time window.
    pub fn zoom_time(&mut self, t_min: f64, t_max: f64) {
        self.time_window = Some((t_min, t_max));
    }

    /// Return only the (t, v) pairs of `channel_idx` that fall within the
    /// current time window.
    pub fn visible_data(&self, channel_idx: usize) -> Vec<(f64, f64)> {
        let (ref _name, ref ts, ref vs) = self.channels[channel_idx];
        let (t0, t1) = self
            .time_window
            .unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
        ts.iter()
            .zip(vs.iter())
            .filter(|&(&t, _)| t >= t0 && t <= t1)
            .map(|(&t, &v)| (t, v))
            .collect()
    }
}

impl Default for TimeSeriesPlot {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Rgba ──────────────────────────────────────────────────────────────

    #[test]
    fn rgba_constructors_valid() {
        assert!(Rgba::black().is_valid());
        assert!(Rgba::white().is_valid());
        assert!(Rgba::red().is_valid());
        assert!(Rgba::green().is_valid());
        assert!(Rgba::blue().is_valid());
    }

    #[test]
    fn rgba_out_of_range_invalid() {
        let c = Rgba::new(1.5, 0.0, 0.0, 1.0);
        assert!(!c.is_valid());
    }

    // ── ColorMapping ──────────────────────────────────────────────────────

    #[test]
    fn color_mapping_clamps_out_of_range() {
        let cm = ColorMapping::new(ColormapKind::Jet, 0.0, 1.0);
        let c = cm.map(5.0); // should clamp to 1.0
        assert!(c.is_valid());
    }

    #[test]
    fn color_mapping_viridis_midpoint_is_valid() {
        let cm = ColorMapping::new(ColormapKind::Viridis, 0.0, 1.0);
        let c = cm.map(0.5);
        assert!(c.is_valid());
    }

    #[test]
    fn color_mapping_bwr_endpoints() {
        let cm = ColorMapping::new(ColormapKind::Bwr, 0.0, 1.0);
        let blue = cm.map(0.0);
        assert!(blue.b > blue.r, "bwr low end should be blue");
        let red = cm.map(1.0);
        assert!(red.r > red.b, "bwr high end should be red");
    }

    #[test]
    fn color_mapping_equal_range_uses_midpoint() {
        let cm = ColorMapping::new(ColormapKind::Greys, 5.0, 5.0);
        let c = cm.map(5.0);
        // midpoint grey = (0.5, 0.5, 0.5, 1.0)
        assert!((c.r - 0.5).abs() < 1e-9);
    }

    #[test]
    fn scatter_all_colors_valid_after_mapping() {
        let mut sp = ScatterPlot::new("test");
        let cm = ColorMapping::new(ColormapKind::Jet, 0.0, 10.0);
        sp.color_mapping = Some(cm);
        for i in 0..5 {
            sp.push_2d(i as f64, i as f64, 3.0, Rgba::black());
        }
        let values: Vec<f64> = (0..5).map(|i| i as f64 * 2.0).collect();
        sp.apply_color_mapping(&values);
        assert!(
            sp.all_colors_valid(),
            "all mapped colours should be in [0,1]"
        );
    }

    // ── PlotRange ─────────────────────────────────────────────────────────

    #[test]
    fn plot_range_auto_empty() {
        let r = PlotRange::Auto.resolve(&[]);
        assert_eq!(r, (0.0, 1.0));
    }

    #[test]
    fn plot_range_auto_single_value() {
        let (lo, hi) = PlotRange::Auto.resolve(&[3.0]);
        assert!(lo < 3.0);
        assert!(hi > 3.0);
    }

    #[test]
    fn plot_range_manual_ignores_data() {
        let r = PlotRange::Manual(-5.0, 5.0).resolve(&[100.0, -100.0]);
        assert_eq!(r, (-5.0, 5.0));
    }

    // ── PlotData ──────────────────────────────────────────────────────────

    #[test]
    fn plot_data_len_min_of_x_y() {
        let pd = PlotData::new(vec![1.0, 2.0, 3.0], vec![4.0, 5.0], "s", Rgba::black());
        assert_eq!(pd.len(), 2);
    }

    #[test]
    fn plot_data_is_empty() {
        let pd = PlotData::new(vec![], vec![], "empty", Rgba::black());
        assert!(pd.is_empty());
    }

    #[test]
    fn plot_data_resolved_ranges_valid() {
        let pd = PlotData::new(vec![1.0, 2.0, 3.0], vec![0.0, 5.0, 2.0], "s", Rgba::black());
        let (xlo, xhi) = pd.resolved_x_range();
        assert!(xlo <= 1.0);
        assert!(xhi >= 3.0);
        let (ylo, yhi) = pd.resolved_y_range();
        assert!(ylo <= 0.0);
        assert!(yhi >= 5.0);
    }

    // ── Histogram ─────────────────────────────────────────────────────────

    #[test]
    fn histogram_count_sums_to_n() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = Histogram::from_data(&data, BinMethod::Sturges, HistNorm::Count);
        assert_eq!(h.count_sum() as usize, 100, "total count must equal N");
    }

    #[test]
    fn histogram_values_in_correct_bins() {
        let data = vec![0.5, 1.5, 2.5, 3.5];
        let h = Histogram::from_data(&data, BinMethod::Fixed(4), HistNorm::Count);
        // Each of 4 bins should have exactly 1 element
        for &c in &h.counts {
            assert_eq!(c as usize, 1, "each bin should contain exactly 1 value");
        }
    }

    #[test]
    fn histogram_pdf_integrates_to_one() {
        let data: Vec<f64> = (0..200).map(|i| (i as f64) * 0.5).collect();
        let h = Histogram::from_data(&data, BinMethod::Sturges, HistNorm::Pdf);
        let area: f64 = h.counts.iter().map(|c| c * h.bin_width()).sum();
        assert!(
            (area - 1.0).abs() < 1e-9,
            "PDF must integrate to 1, got {:.6}",
            area
        );
    }

    #[test]
    fn histogram_cdf_ends_at_one() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let h = Histogram::from_data(&data, BinMethod::Sturges, HistNorm::Cdf);
        let last = *h.counts.last().unwrap();
        assert!(
            (last - 1.0).abs() < 1e-9,
            "CDF last bin should be 1.0, got {:.6}",
            last
        );
    }

    #[test]
    fn histogram_fixed_bins_count() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = Histogram::from_data(&data, BinMethod::Fixed(10), HistNorm::Count);
        assert_eq!(h.num_bins(), 10);
    }

    #[test]
    fn histogram_scott_nonempty() {
        let data: Vec<f64> = (0..50).map(|i| (i as f64 - 25.0) * 0.1).collect();
        let h = Histogram::from_data(&data, BinMethod::Scott, HistNorm::Count);
        assert!(h.num_bins() > 0);
        assert_eq!(h.count_sum() as usize, 50);
    }

    #[test]
    fn histogram_fd_nonempty() {
        let data: Vec<f64> = (0..60).map(|i| i as f64).collect();
        let h = Histogram::from_data(&data, BinMethod::FreedmanDiaconis, HistNorm::Count);
        assert!(h.num_bins() > 0);
        assert_eq!(h.count_sum() as usize, 60);
    }

    // ── BarChart ──────────────────────────────────────────────────────────

    #[test]
    fn bar_chart_stacked_totals_sum_correctly() {
        let mut bc = BarChart::new("test", true);
        bc.push_group(vec![
            Bar::new("a", 3.0, Rgba::red()),
            Bar::new("b", 7.0, Rgba::blue()),
        ]);
        bc.push_group(vec![
            Bar::new("c", 5.0, Rgba::red()),
            Bar::new("d", 5.0, Rgba::blue()),
        ]);
        let totals = bc.stacked_totals();
        assert_eq!(totals.len(), 2);
        assert!(
            (totals[0] - 10.0).abs() < 1e-9,
            "group 0 total should be 10"
        );
        assert!(
            (totals[1] - 10.0).abs() < 1e-9,
            "group 1 total should be 10"
        );
    }

    #[test]
    fn bar_with_error_stores_error() {
        let b = Bar::with_error("x", 5.0, 0.5, Rgba::black());
        assert!(b.error.is_some());
        assert!((b.error.unwrap() - 0.5).abs() < 1e-9);
    }

    // ── HeatmapData ───────────────────────────────────────────────────────

    #[test]
    fn heatmap_size_matches_grid() {
        let rows = 4;
        let cols = 5;
        let data: Vec<f64> = (0..rows * cols).map(|i| i as f64).collect();
        let cm = ColorMapping::new(ColormapKind::Viridis, 0.0, 20.0);
        let hm = HeatmapData::new(data, rows, cols, cm, NormalizeRange::MinMax);
        let bytes = hm.to_rgba_bytes();
        assert_eq!(
            bytes.len(),
            rows * cols * 4,
            "RGBA byte buffer length must match"
        );
    }

    #[test]
    fn heatmap_symmetric_range_is_symmetric() {
        let data = vec![-5.0, 0.0, 3.0];
        let cm = ColorMapping::new(ColormapKind::Bwr, -5.0, 5.0);
        let hm = HeatmapData::new(data, 1, 3, cm, NormalizeRange::Symmetric);
        let (lo, hi) = hm.resolve_range();
        assert!(
            (lo + hi).abs() < 1e-9,
            "symmetric range must be centred at 0"
        );
    }

    #[test]
    fn heatmap_percentile_range_ordered() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let cm = ColorMapping::new(ColormapKind::Greys, 0.0, 1.0);
        let hm = HeatmapData::new(data, 10, 10, cm, NormalizeRange::Percentile(2.0, 98.0));
        let (lo, hi) = hm.resolve_range();
        assert!(lo < hi, "percentile lo must be less than hi");
    }

    // ── ContourPlot ───────────────────────────────────────────────────────

    #[test]
    fn contour_empty_for_iso_above_all_values() {
        // 2×2 grid all zeros: iso = 1.0 → no crossings
        let data = vec![0.0, 0.0, 0.0, 0.0];
        let cp = ContourPlot::new(data, 2, 2, [0.0, 1.0, 0.0, 1.0]);
        let segs = cp.extract(1.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn contour_circle_has_segments() {
        // Build a 21×21 grid: value = -(x^2 + y^2 - 0.25); iso at 0 = circle r=0.5
        let n = 21usize;
        let mut data = Vec::with_capacity(n * n);
        for r in 0..n {
            for c in 0..n {
                let x = -1.0 + 2.0 * c as f64 / (n - 1) as f64;
                let y = -1.0 + 2.0 * r as f64 / (n - 1) as f64;
                data.push(0.25 - (x * x + y * y));
            }
        }
        let cp = ContourPlot::new(data, n, n, [-1.0, 1.0, -1.0, 1.0]);
        let segs = cp.extract(0.0);
        assert!(!segs.is_empty(), "circle should produce contour segments");
    }

    #[test]
    fn contour_single_crossing_produces_segment() {
        // 2×2 grid: left column < iso, right column > iso
        let data = vec![0.0, 2.0, 0.0, 2.0];
        let cp = ContourPlot::new(data, 2, 2, [0.0, 1.0, 0.0, 1.0]);
        let segs = cp.extract(1.0);
        assert_eq!(segs.len(), 1);
    }

    // ── VectorField2d ─────────────────────────────────────────────────────

    #[test]
    fn vector_field_arrow_count_matches_grid() {
        let vf = VectorField2d::from_grid((0.0, 1.0), (0.0, 1.0), 4, 3, 1.0, None, |x, y| [y, -x]);
        assert_eq!(vf.arrows.len(), 12, "4×3 grid should give 12 arrows");
    }

    #[test]
    fn vector_field_colors_valid() {
        let cm = ColorMapping::new(ColormapKind::Hot, 0.0, 2.0);
        let vf =
            VectorField2d::from_grid((-1.0, 1.0), (-1.0, 1.0), 3, 3, 1.0, Some(&cm), |x, y| {
                [x, y]
            });
        for a in &vf.arrows {
            assert!(a.color.is_valid(), "arrow colour out of [0,1] range");
        }
    }

    // ── ParallelCoordinates ───────────────────────────────────────────────

    #[test]
    fn parallel_coordinates_normalized_in_unit_range() {
        let axes: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let mut pc = ParallelCoordinates::new(axes, true);
        pc.push(vec![1.0, 10.0, 100.0], Rgba::red(), None);
        pc.push(vec![2.0, 20.0, 200.0], Rgba::blue(), None);
        let nv = pc.normalized_row(0);
        for &v in &nv {
            assert!(
                (0.0..=1.0).contains(&v),
                "normalized value {:.6} out of [0,1]",
                v
            );
        }
    }

    #[test]
    fn parallel_coordinates_first_row_normalizes_to_zero() {
        let axes: Vec<String> = ["x"].iter().map(|s| s.to_string()).collect();
        let mut pc = ParallelCoordinates::new(axes, true);
        pc.push(vec![0.0], Rgba::black(), None);
        pc.push(vec![10.0], Rgba::white(), None);
        let nv = pc.normalized_row(0);
        assert!((nv[0] - 0.0).abs() < 1e-9, "min row should normalize to 0");
    }

    // ── RadarChart ────────────────────────────────────────────────────────

    #[test]
    fn radar_polygon_points_same_scale() {
        let axes: Vec<String> = ["s", "e", "c", "r", "e"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut rc = RadarChart::new(axes, 10.0);
        rc.push("hero", vec![10.0, 10.0, 10.0, 10.0, 10.0], Rgba::red());
        let pts = rc.polygon_points(0);
        // All values at max → all points at unit circle
        for p in &pts {
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-9,
                "point radius should be 1.0, got {:.6}",
                r
            );
        }
    }

    #[test]
    fn radar_chart_all_axes_same_scale() {
        let axes: Vec<String> = ["a", "b"].iter().map(|s| s.to_string()).collect();
        let mut rc = RadarChart::new(axes, 5.0);
        rc.push("s1", vec![5.0, 5.0], Rgba::blue());
        rc.push("s2", vec![2.5, 2.5], Rgba::green());
        let pts1 = rc.polygon_points(0);
        let pts2 = rc.polygon_points(1);
        let r1: Vec<f64> = pts1
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .collect();
        let r2: Vec<f64> = pts2
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt())
            .collect();
        for (&ra, &rb) in r1.iter().zip(r2.iter()) {
            assert!(
                ra > rb - 1e-9,
                "series 1 at max should exceed series 2 at half"
            );
        }
    }

    // ── TimeSeriesPlot ────────────────────────────────────────────────────

    #[test]
    fn time_series_visible_data_filters_by_window() {
        let mut ts = TimeSeriesPlot::new();
        let t: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let v: Vec<f64> = (0..10).map(|i| i as f64).collect();
        ts.add_channel("ch", t, v);
        ts.zoom_time(2.0, 5.0);
        let vis = ts.visible_data(0);
        assert_eq!(vis.len(), 4, "should have 4 points in [2,5]");
    }

    #[test]
    fn time_series_no_window_returns_all() {
        let mut ts = TimeSeriesPlot::new();
        ts.add_channel("ch", vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0]);
        let vis = ts.visible_data(0);
        assert_eq!(vis.len(), 3);
    }

    #[test]
    fn time_series_marker_added() {
        let mut ts = TimeSeriesPlot::new();
        ts.add_marker(1.5, "event", Rgba::red());
        assert_eq!(ts.markers.len(), 1);
        assert!((ts.markers[0].t - 1.5).abs() < 1e-9);
    }

    // ── ScatterPlot ───────────────────────────────────────────────────────

    #[test]
    fn scatter_3d_point_stored() {
        let mut sp = ScatterPlot::new("s");
        sp.push_3d(1.0, 2.0, 3.0, 4.0, Rgba::blue());
        assert_eq!(sp.points.len(), 1);
        assert!(sp.points[0].z.is_some());
        assert!((sp.points[0].z.unwrap() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn scatter_2d_z_is_none() {
        let mut sp = ScatterPlot::new("s");
        sp.push_2d(0.0, 0.0, 1.0, Rgba::black());
        assert!(sp.points[0].z.is_none());
    }

    // ── integration / round-trip ──────────────────────────────────────────

    #[test]
    fn histogram_then_pdf_then_cdf_consistent() {
        let data: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let h_count = Histogram::from_data(&data, BinMethod::Fixed(5), HistNorm::Count);
        let h_pdf = Histogram::from_data(&data, BinMethod::Fixed(5), HistNorm::Pdf);
        let h_cdf = Histogram::from_data(&data, BinMethod::Fixed(5), HistNorm::Cdf);
        assert_eq!(h_count.num_bins(), h_pdf.num_bins());
        assert_eq!(h_pdf.num_bins(), h_cdf.num_bins());
        // CDF last bin = 1.0
        assert!((h_cdf.counts.last().unwrap() - 1.0).abs() < 1e-9);
    }
}
