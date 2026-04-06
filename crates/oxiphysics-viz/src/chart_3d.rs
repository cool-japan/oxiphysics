// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! 3D chart rendering data structures.
//!
//! Provides 3D bar charts, surface plots, 3D scatter plots, vector field
//! plots, colormap utilities, and axis line-segment generation.

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// BarChart3D
// ─────────────────────────────────────────────────────────────────────────────

/// 3D bar chart data: bars stored as `(x_index, y_index, value)` triples.
#[derive(Debug, Clone, Default)]
pub struct BarChart3D {
    /// Bars as `(x_col, y_row, value)`.
    pub bars: Vec<(f64, f64, f64)>,
    /// Labels for the X axis categories.
    pub x_labels: Vec<String>,
    /// Labels for the Y axis categories.
    pub y_labels: Vec<String>,
}

impl BarChart3D {
    /// Create a new empty bar chart.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a bar at `(x, y)` with the given `value`.
    pub fn add_bar(&mut self, x: f64, y: f64, value: f64) {
        self.bars.push((x, y, value));
    }

    /// Number of bars in this chart.
    pub fn bar_count(&self) -> usize {
        self.bars.len()
    }

    /// Maximum value across all bars, or `f64::NEG_INFINITY` if empty.
    pub fn max_value(&self) -> f64 {
        self.bars
            .iter()
            .map(|&(_, _, v)| v)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum value across all bars, or `f64::INFINITY` if empty.
    pub fn min_value(&self) -> f64 {
        self.bars
            .iter()
            .map(|&(_, _, v)| v)
            .fold(f64::INFINITY, f64::min)
    }

    /// Normalize bar values to `[0, 1]` relative to the chart's min/max.
    ///
    /// Returns `(x_index_as_usize, y_index_as_usize, normalized_value)`.
    /// If all values are equal, all normalized values are `0.5`.
    pub fn normalize(&self) -> Vec<(usize, usize, f64)> {
        let lo = self.min_value();
        let hi = self.max_value();
        let range = hi - lo;
        self.bars
            .iter()
            .map(|&(x, y, v)| {
                let norm = if range.abs() < 1e-15 {
                    0.5
                } else {
                    (v - lo) / range
                };
                (x as usize, y as usize, norm)
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SurfacePlot
// ─────────────────────────────────────────────────────────────────────────────

/// A surface plot `z = f(x, y)` sampled on a regular grid.
#[derive(Debug, Clone)]
pub struct SurfacePlot {
    /// X-axis sample points.
    pub x: Vec<f64>,
    /// Y-axis sample points.
    pub y: Vec<f64>,
    /// Z values: `z[iy][ix]` is the height at `(x[ix], y[iy])`.
    pub z: Vec<Vec<f64>>,
}

impl SurfacePlot {
    /// Build a surface plot by evaluating `f` on an `nx × ny` grid.
    ///
    /// `x_range` and `y_range` are `(min, max)` inclusive bounds.
    pub fn from_fn(
        nx: usize,
        ny: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
        f: impl Fn(f64, f64) -> f64,
    ) -> Self {
        let nx = nx.max(2);
        let ny = ny.max(2);
        let x: Vec<f64> = (0..nx)
            .map(|i| x_range.0 + (x_range.1 - x_range.0) * i as f64 / (nx - 1) as f64)
            .collect();
        let y: Vec<f64> = (0..ny)
            .map(|j| y_range.0 + (y_range.1 - y_range.0) * j as f64 / (ny - 1) as f64)
            .collect();
        let z: Vec<Vec<f64>> = y
            .iter()
            .map(|&yv| x.iter().map(|&xv| f(xv, yv)).collect())
            .collect();
        Self { x, y, z }
    }

    /// Number of samples along the X axis.
    pub fn nx(&self) -> usize {
        self.x.len()
    }

    /// Number of samples along the Y axis.
    pub fn ny(&self) -> usize {
        self.y.len()
    }

    /// Minimum Z value over the entire grid.
    pub fn z_min(&self) -> f64 {
        self.z
            .iter()
            .flat_map(|row| row.iter())
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }

    /// Maximum Z value over the entire grid.
    pub fn z_max(&self) -> f64 {
        self.z
            .iter()
            .flat_map(|row| row.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Triangulate the surface into triangle index triples `[i0, i1, i2]`.
    ///
    /// Each grid cell `(ix, iy)` → `(ix, iy+1)` → `(ix+1, iy+1)` → `(ix+1, iy)` is
    /// split into two triangles. Total triangle count = `2 * (nx-1) * (ny-1)`.
    ///
    /// Vertex index: `iy * nx + ix`.
    pub fn triangulate(&self) -> Vec<[usize; 3]> {
        let nx = self.nx();
        let ny = self.ny();
        if nx < 2 || ny < 2 {
            return vec![];
        }
        let mut tris = Vec::with_capacity(2 * (nx - 1) * (ny - 1));
        for iy in 0..(ny - 1) {
            for ix in 0..(nx - 1) {
                let i00 = iy * nx + ix;
                let i10 = iy * nx + ix + 1;
                let i01 = (iy + 1) * nx + ix;
                let i11 = (iy + 1) * nx + ix + 1;
                tris.push([i00, i10, i11]);
                tris.push([i00, i11, i01]);
            }
        }
        tris
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ScatterPlot3D
// ─────────────────────────────────────────────────────────────────────────────

/// A 3D scatter plot with per-point colors and sizes.
#[derive(Debug, Clone, Default)]
pub struct ScatterPlot3D {
    /// 3D point positions.
    pub points: Vec<[f64; 3]>,
    /// Per-point RGBA colors.
    pub colors: Vec<[f32; 4]>,
    /// Per-point billboard sizes.
    pub sizes: Vec<f32>,
}

impl ScatterPlot3D {
    /// Create an empty scatter plot.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a point with the given position, color, and size.
    pub fn add_point(&mut self, xyz: [f64; 3], color: [f32; 4], size: f32) {
        self.points.push(xyz);
        self.colors.push(color);
        self.sizes.push(size);
    }

    /// Number of points in this scatter plot.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Compute the axis-aligned bounding box `(min_xyz, max_xyz)`.
    ///
    /// Returns `([0,0,0], [0,0,0])` if the plot is empty.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.points.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = self.points[0];
        let mut mx = self.points[0];
        for p in &self.points {
            for k in 0..3 {
                if p[k] < mn[k] {
                    mn[k] = p[k];
                }
                if p[k] > mx[k] {
                    mx[k] = p[k];
                }
            }
        }
        (mn, mx)
    }

    /// Compute the centroid (mean position) of all points.
    ///
    /// Returns `[0, 0, 0]` if empty.
    pub fn center(&self) -> [f64; 3] {
        if self.points.is_empty() {
            return [0.0; 3];
        }
        let n = self.points.len() as f64;
        let mut sum = [0.0_f64; 3];
        for p in &self.points {
            sum[0] += p[0];
            sum[1] += p[1];
            sum[2] += p[2];
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VectorFieldPlot
// ─────────────────────────────────────────────────────────────────────────────

/// A vector field visualization: arrows at given origins.
#[derive(Debug, Clone)]
pub struct VectorFieldPlot {
    /// Arrow origins.
    pub origins: Vec<[f64; 3]>,
    /// Arrow vectors (unscaled).
    pub vectors: Vec<[f64; 3]>,
    /// Scale factor applied when rendering.
    pub scale: f64,
}

impl VectorFieldPlot {
    /// Create an empty vector field plot with the given scale.
    pub fn new(scale: f64) -> Self {
        Self {
            origins: Vec::new(),
            vectors: Vec::new(),
            scale,
        }
    }

    /// Add an arrow at `origin` with direction/magnitude `vec`.
    pub fn add_arrow(&mut self, origin: [f64; 3], vec: [f64; 3]) {
        self.origins.push(origin);
        self.vectors.push(vec);
    }

    /// Number of arrows in this vector field.
    pub fn arrow_count(&self) -> usize {
        self.origins.len()
    }

    /// Maximum vector magnitude across all arrows.
    ///
    /// Returns `0.0` if there are no arrows.
    pub fn max_magnitude(&self) -> f64 {
        self.vectors
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }

    /// Return arrows as `(origin, tip)` pairs with vectors normalized to unit length and scaled.
    ///
    /// Zero-length vectors are left at the origin (tip == origin).
    pub fn normalized_arrows(&self) -> Vec<([f64; 3], [f64; 3])> {
        self.origins
            .iter()
            .zip(self.vectors.iter())
            .map(|(o, v)| {
                let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if mag < 1e-15 {
                    (*o, *o)
                } else {
                    let tip = [
                        o[0] + v[0] / mag * self.scale,
                        o[1] + v[1] / mag * self.scale,
                        o[2] + v[2] / mag * self.scale,
                    ];
                    (*o, tip)
                }
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// color_by_value
// ─────────────────────────────────────────────────────────────────────────────

/// Map a scalar `value` in `[min, max]` to an RGBA color using a named colormap.
///
/// Supported colormap names (case-insensitive): `"viridis"`, `"jet"`, `"hot"`,
/// `"cool"`, `"gray"`. Unknown names fall back to grayscale.
pub fn color_by_value(value: f64, min: f64, max: f64, colormap: &str) -> [f32; 4] {
    let t = if (max - min).abs() < 1e-15 {
        0.5
    } else {
        ((value - min) / (max - min)).clamp(0.0, 1.0)
    } as f32;

    match colormap.to_lowercase().as_str() {
        "viridis" => viridis_color(t),
        "jet" => jet_color(t),
        "hot" => hot_color(t),
        "cool" => cool_color(t),
        "gray" | "grey" => [t, t, t, 1.0],
        _ => [t, t, t, 1.0],
    }
}

/// Viridis colormap approximation.
fn viridis_color(t: f32) -> [f32; 4] {
    // Three-stop approximation: dark purple → teal → yellow-green
    let r = (0.267 + 0.003 * t + 1.394 * t * t - 0.664 * t * t * t).clamp(0.0, 1.0);
    let g = (0.005 + 1.237 * t - 0.660 * t * t + 0.418 * t * t * t).clamp(0.0, 1.0);
    let b = (0.329 + 1.294 * t - 2.638 * t * t + 1.415 * t * t * t).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// Jet colormap: blue → cyan → green → yellow → red.
///
/// t=0 → pure blue (0,0,1), t=0.5 → green, t=1 → pure red (1,0,0).
fn jet_color(t: f32) -> [f32; 4] {
    let r = if t < 0.5 {
        0.0
    } else {
        ((t - 0.5) * 4.0).clamp(0.0, 1.0)
    };
    let g = if t < 0.25 {
        (t * 4.0).clamp(0.0, 1.0)
    } else if t < 0.75 {
        1.0
    } else {
        ((1.0 - t) * 4.0).clamp(0.0, 1.0)
    };
    let b = (1.0 - (t - 0.25).max(0.0) * 4.0).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// Hot colormap: black → red → yellow → white.
fn hot_color(t: f32) -> [f32; 4] {
    let r = (t * 3.0).clamp(0.0, 1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// Cool colormap: cyan → magenta.
fn cool_color(t: f32) -> [f32; 4] {
    [t, 1.0 - t, 1.0, 1.0]
}

// ─────────────────────────────────────────────────────────────────────────────
// chart_axes_3d
// ─────────────────────────────────────────────────────────────────────────────

/// Generate axis line segments for a 3D chart.
///
/// Returns tick mark line segments for X, Y, and Z axes.  Each tick is a short
/// segment perpendicular to the axis.  `n_ticks` ticks are placed on each axis.
pub fn chart_axes_3d(
    x_range: (f64, f64),
    y_range: (f64, f64),
    z_range: (f64, f64),
    n_ticks: usize,
) -> Vec<([f64; 3], [f64; 3])> {
    let n = n_ticks.max(2);
    let mut segs: Vec<([f64; 3], [f64; 3])> = Vec::new();

    // X axis spine
    segs.push(([x_range.0, 0.0, 0.0], [x_range.1, 0.0, 0.0]));
    // Y axis spine
    segs.push(([0.0, y_range.0, 0.0], [0.0, y_range.1, 0.0]));
    // Z axis spine
    segs.push(([0.0, 0.0, z_range.0], [0.0, 0.0, z_range.1]));

    let tick_len = 0.05;

    // X ticks
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let x = x_range.0 + t * (x_range.1 - x_range.0);
        segs.push(([x, -tick_len, 0.0], [x, tick_len, 0.0]));
    }

    // Y ticks
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let y = y_range.0 + t * (y_range.1 - y_range.0);
        segs.push(([-tick_len, y, 0.0], [tick_len, y, 0.0]));
    }

    // Z ticks
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let z = z_range.0 + t * (z_range.1 - z_range.0);
        segs.push(([0.0, -tick_len, z], [0.0, tick_len, z]));
    }

    segs
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BarChart3D ────────────────────────────────────────────────────────────

    #[test]
    fn test_bar_chart_empty() {
        let chart = BarChart3D::new();
        assert_eq!(chart.bar_count(), 0);
    }

    #[test]
    fn test_bar_chart_add_bar() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 5.0);
        assert_eq!(chart.bar_count(), 1);
    }

    #[test]
    fn test_bar_chart_max_value() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 3.0);
        chart.add_bar(1.0, 0.0, 7.0);
        chart.add_bar(0.0, 1.0, 1.0);
        assert!((chart.max_value() - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_bar_chart_min_value() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 3.0);
        chart.add_bar(1.0, 0.0, 7.0);
        chart.add_bar(0.0, 1.0, 1.0);
        assert!((chart.min_value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bar_chart_normalize_count() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 1.0);
        chart.add_bar(1.0, 0.0, 2.0);
        chart.add_bar(2.0, 0.0, 3.0);
        let norm = chart.normalize();
        assert_eq!(norm.len(), 3);
    }

    #[test]
    fn test_bar_chart_normalize_range() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 0.0);
        chart.add_bar(1.0, 0.0, 10.0);
        let norm = chart.normalize();
        // min normalizes to 0.0, max to 1.0
        let values: Vec<f64> = norm.iter().map(|&(_, _, v)| v).collect();
        let min_v = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_v = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(min_v.abs() < 1e-12);
        assert!((max_v - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bar_chart_normalize_all_equal() {
        let mut chart = BarChart3D::new();
        chart.add_bar(0.0, 0.0, 5.0);
        chart.add_bar(1.0, 0.0, 5.0);
        let norm = chart.normalize();
        for (_, _, v) in &norm {
            assert!((v - 0.5).abs() < 1e-12);
        }
    }

    // ── SurfacePlot ───────────────────────────────────────────────────────────

    #[test]
    fn test_surface_plot_dimensions() {
        let sp = SurfacePlot::from_fn(5, 4, (0.0, 1.0), (0.0, 1.0), |x, y| x + y);
        assert_eq!(sp.nx(), 5);
        assert_eq!(sp.ny(), 4);
    }

    #[test]
    fn test_surface_plot_z_shape() {
        let sp = SurfacePlot::from_fn(3, 3, (0.0, 1.0), (0.0, 1.0), |x, y| x * y);
        assert_eq!(sp.z.len(), 3); // ny rows
        assert_eq!(sp.z[0].len(), 3); // nx cols
    }

    #[test]
    fn test_surface_plot_z_min_max() {
        let sp = SurfacePlot::from_fn(10, 10, (0.0, 1.0), (0.0, 1.0), |x, _y| x);
        assert!(sp.z_min() >= 0.0 - 1e-12);
        assert!(sp.z_max() <= 1.0 + 1e-12);
    }

    #[test]
    fn test_surface_plot_triangulate_count() {
        let nx = 5usize;
        let ny = 4usize;
        let sp = SurfacePlot::from_fn(nx, ny, (0.0, 1.0), (0.0, 1.0), |_, _| 0.0);
        let tris = sp.triangulate();
        let expected = 2 * (nx - 1) * (ny - 1);
        assert_eq!(tris.len(), expected);
    }

    #[test]
    fn test_surface_plot_triangulate_2x2() {
        // 2x2 grid → 2*(1)*(1) = 2 triangles
        let sp = SurfacePlot::from_fn(2, 2, (0.0, 1.0), (0.0, 1.0), |_, _| 0.0);
        assert_eq!(sp.triangulate().len(), 2);
    }

    #[test]
    fn test_surface_plot_triangulate_large() {
        let sp = SurfacePlot::from_fn(20, 15, (0.0, 1.0), (0.0, 1.0), |x, y| x.sin() * y.cos());
        let expected = 2 * 19 * 14;
        assert_eq!(sp.triangulate().len(), expected);
    }

    #[test]
    fn test_surface_plot_x_range() {
        let sp = SurfacePlot::from_fn(5, 3, (-1.0, 1.0), (0.0, 1.0), |_, _| 0.0);
        assert!((sp.x[0] - (-1.0)).abs() < 1e-12);
        assert!((sp.x[4] - 1.0).abs() < 1e-12);
    }

    // ── ScatterPlot3D ─────────────────────────────────────────────────────────

    #[test]
    fn test_scatter_empty() {
        let s = ScatterPlot3D::new();
        assert_eq!(s.point_count(), 0);
    }

    #[test]
    fn test_scatter_add_point() {
        let mut s = ScatterPlot3D::new();
        s.add_point([1.0, 2.0, 3.0], [1.0, 0.0, 0.0, 1.0], 0.1);
        assert_eq!(s.point_count(), 1);
    }

    #[test]
    fn test_scatter_bounding_box() {
        let mut s = ScatterPlot3D::new();
        s.add_point([0.0, 0.0, 0.0], [1.0; 4], 0.1);
        s.add_point([3.0, -1.0, 5.0], [1.0; 4], 0.1);
        let (mn, mx) = s.bounding_box();
        assert!((mn[0] - 0.0).abs() < 1e-12);
        assert!((mx[0] - 3.0).abs() < 1e-12);
        assert!((mn[1] - (-1.0)).abs() < 1e-12);
        assert!((mx[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_scatter_bounding_box_empty() {
        let s = ScatterPlot3D::new();
        let (mn, mx) = s.bounding_box();
        for k in 0..3 {
            assert!(mn[k].abs() < 1e-12);
            assert!(mx[k].abs() < 1e-12);
        }
    }

    #[test]
    fn test_scatter_center() {
        let mut s = ScatterPlot3D::new();
        s.add_point([0.0, 0.0, 0.0], [1.0; 4], 0.1);
        s.add_point([2.0, 4.0, 6.0], [1.0; 4], 0.1);
        let c = s.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_scatter_center_empty() {
        let s = ScatterPlot3D::new();
        let c = s.center();
        assert!(c[0].abs() < 1e-12);
    }

    // ── VectorFieldPlot ───────────────────────────────────────────────────────

    #[test]
    fn test_vector_field_empty() {
        let vf = VectorFieldPlot::new(1.0);
        assert_eq!(vf.arrow_count(), 0);
        assert!(vf.max_magnitude().abs() < 1e-12);
    }

    #[test]
    fn test_vector_field_add_arrow() {
        let mut vf = VectorFieldPlot::new(1.0);
        vf.add_arrow([0.0; 3], [1.0, 0.0, 0.0]);
        assert_eq!(vf.arrow_count(), 1);
    }

    #[test]
    fn test_vector_field_max_magnitude() {
        let mut vf = VectorFieldPlot::new(1.0);
        vf.add_arrow([0.0; 3], [3.0, 4.0, 0.0]); // magnitude = 5
        vf.add_arrow([0.0; 3], [1.0, 0.0, 0.0]); // magnitude = 1
        assert!((vf.max_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_normalized_arrows_length() {
        let mut vf = VectorFieldPlot::new(2.0);
        vf.add_arrow([0.0; 3], [3.0, 0.0, 0.0]);
        let arrows = vf.normalized_arrows();
        let (o, tip) = arrows[0];
        // Normalized direction is [1,0,0], scaled by 2.0
        let dx = tip[0] - o[0];
        assert!((dx - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_field_zero_vector() {
        let mut vf = VectorFieldPlot::new(1.0);
        vf.add_arrow([1.0, 2.0, 3.0], [0.0; 3]);
        let arrows = vf.normalized_arrows();
        let (o, tip) = arrows[0];
        // tip == origin for zero vector
        for k in 0..3 {
            assert!((tip[k] - o[k]).abs() < 1e-12);
        }
    }

    // ── color_by_value ────────────────────────────────────────────────────────

    #[test]
    fn test_color_by_value_jet_min() {
        let c = color_by_value(0.0, 0.0, 1.0, "jet");
        // jet(0) is blue: r≈0, b≈1
        assert!(c[2] > 0.9, "jet(0) should be blue, b={}", c[2]);
        assert!(c[3].abs() - 1.0 < 1e-6);
    }

    #[test]
    fn test_color_by_value_jet_max() {
        let c = color_by_value(1.0, 0.0, 1.0, "jet");
        // jet(1) is red: r≈1, b≈0
        assert!(c[0] > 0.9, "jet(1) should be red, r={}", c[0]);
    }

    #[test]
    fn test_color_by_value_gray() {
        let c = color_by_value(0.5, 0.0, 1.0, "gray");
        assert!((c[0] - 0.5).abs() < 1e-6);
        assert!((c[1] - 0.5).abs() < 1e-6);
        assert!((c[2] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_color_by_value_equal_range() {
        // All values equal → t = 0.5 fallback
        let c = color_by_value(5.0, 5.0, 5.0, "gray");
        assert!((c[0] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_color_by_value_unknown_colormap() {
        // Unknown colormap: gray fallback
        let c = color_by_value(1.0, 0.0, 1.0, "unknown_map");
        assert!((c[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_color_by_value_hot_zero() {
        let c = color_by_value(0.0, 0.0, 1.0, "hot");
        assert!(c[0].abs() < 1e-6); // black at start
    }

    #[test]
    fn test_color_by_value_cool_endpoints() {
        let c0 = color_by_value(0.0, 0.0, 1.0, "cool");
        let c1 = color_by_value(1.0, 0.0, 1.0, "cool");
        // At t=0: cyan [0,1,1]; at t=1: magenta [1,0,1]
        assert!(c0[1] > 0.9);
        assert!(c1[0] > 0.9);
    }

    // ── chart_axes_3d ─────────────────────────────────────────────────────────

    #[test]
    fn test_chart_axes_3d_segment_count() {
        // 3 axis spines + 3 * n_ticks tick marks
        let n = 5usize;
        let segs = chart_axes_3d((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), n);
        let expected = 3 + 3 * n;
        assert_eq!(segs.len(), expected);
    }

    #[test]
    fn test_chart_axes_3d_min_ticks() {
        // n_ticks=0 → clamped to 2
        let segs = chart_axes_3d((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), 0);
        assert_eq!(segs.len(), 3 + 3 * 2);
    }

    #[test]
    fn test_chart_axes_3d_x_spine_endpoints() {
        let segs = chart_axes_3d((-1.0, 2.0), (0.0, 1.0), (0.0, 1.0), 2);
        // First segment should be the X spine: from (-1,0,0) to (2,0,0)
        let (start, end) = segs[0];
        assert!((start[0] - (-1.0)).abs() < 1e-12);
        assert!((end[0] - 2.0).abs() < 1e-12);
    }
}
