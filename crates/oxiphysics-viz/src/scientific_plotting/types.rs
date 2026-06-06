//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

/// Panel index `(row, col)` → axis range `((x_min, x_max), (y_min, y_max))`.
pub type PanelAxisRanges = HashMap<(usize, usize), ((f64, f64), (f64, f64))>;

/// A vertex on a 3D surface.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceVertex {
    /// 3D position.
    pub position: [f64; 3],
    /// Surface normal.
    pub normal: [f64; 3],
    /// Scalar value for colormap.
    pub value: f64,
}
/// Colormap name for scalar-to-color mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlotColormap {
    /// Perceptually uniform, dark-to-bright.
    #[default]
    Viridis,
    /// Rainbow (blue→cyan→green→yellow→red).
    Jet,
    /// Diverging: blue–white–red.
    BlueWhiteRed,
    /// Grey scale.
    Greys,
    /// Heat / plasma.
    Plasma,
    /// Inferno.
    Inferno,
    /// Coolwarm.
    Coolwarm,
}
impl PlotColormap {
    /// Map a value `t ∈ [0, 1]` to a color.
    pub fn map(&self, t: f32) -> PlotColor {
        let t = t.clamp(0.0, 1.0);
        match self {
            PlotColormap::Jet => {
                let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
                let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
                let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
                PlotColor::rgb(r, g, b)
            }
            PlotColormap::BlueWhiteRed => {
                if t < 0.5 {
                    PlotColor::lerp(PlotColor::blue(), PlotColor::white(), t * 2.0)
                } else {
                    PlotColor::lerp(PlotColor::white(), PlotColor::red(), (t - 0.5) * 2.0)
                }
            }
            PlotColormap::Greys => PlotColor::rgb(t, t, t),
            PlotColormap::Plasma => {
                let r = (0.05 + 0.95 * t).clamp(0.0, 1.0);
                let g = (0.02 + 0.3 * (t * std::f32::consts::PI).sin()).clamp(0.0, 1.0);
                let b = (0.55 - 0.55 * t).clamp(0.0, 1.0);
                PlotColor::rgb(r, g, b)
            }
            PlotColormap::Inferno => {
                let r = (t * 1.2).clamp(0.0, 1.0);
                let g = (t * t).clamp(0.0, 1.0);
                let b = (0.3 * (1.0 - t)).clamp(0.0, 1.0);
                PlotColor::rgb(r, g, b)
            }
            PlotColormap::Coolwarm => {
                if t < 0.5 {
                    let s = t * 2.0;
                    PlotColor::lerp(PlotColor::rgb(0.23, 0.30, 0.75), PlotColor::white(), s)
                } else {
                    let s = (t - 0.5) * 2.0;
                    PlotColor::lerp(PlotColor::white(), PlotColor::rgb(0.70, 0.09, 0.09), s)
                }
            }
            PlotColormap::Viridis => {
                let r = (-2.0 * (t - 0.9) * (t - 0.9)).exp() * 0.9;
                let g = (-(3.5 * t - 2.5) * (3.5 * t - 2.5)).exp() * 0.8
                    + (-(3.0 * t - 0.5) * (3.0 * t - 0.5)).exp() * 0.5;
                let b = (-(3.0 * t) * (3.0 * t)).exp() * 0.5 + 0.1;
                PlotColor::rgb(r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0))
            }
        }
    }
}
/// Axis configuration for a 2D/3D plot.
#[derive(Debug, Clone)]
pub struct AxisConfig {
    /// Axis label text.
    pub label: String,
    /// Optional fixed minimum for the axis range. `None` → auto.
    pub min: Option<f64>,
    /// Optional fixed maximum for the axis range. `None` → auto.
    pub max: Option<f64>,
    /// Scale type for this axis.
    pub scale: AxisScale,
    /// Number of major tick marks (0 = auto).
    pub num_ticks: usize,
    /// Whether to show a grid line for this axis.
    pub show_grid: bool,
}
impl AxisConfig {
    /// Create an axis with a given label.
    pub fn with_label(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            ..Default::default()
        }
    }
    /// Set fixed data range.
    pub fn with_range(mut self, min: f64, max: f64) -> Self {
        self.min = Some(min);
        self.max = Some(max);
        self
    }
    /// Set axis scale.
    pub fn with_scale(mut self, scale: AxisScale) -> Self {
        self.scale = scale;
        self
    }
    /// Compute auto range from a data slice, respecting overrides.
    pub fn effective_range(&self, data: &[f64]) -> (f64, f64) {
        let lo = self
            .min
            .unwrap_or_else(|| data.iter().cloned().fold(f64::INFINITY, f64::min));
        let hi = self
            .max
            .unwrap_or_else(|| data.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        if (hi - lo).abs() < 1e-14 {
            (lo - 1.0, hi + 1.0)
        } else {
            (lo, hi)
        }
    }
}
/// Lighting model for a 3D surface.
#[derive(Debug, Clone)]
pub struct SurfaceLighting {
    /// Ambient coefficient (0–1).
    pub ambient: f32,
    /// Diffuse coefficient (0–1).
    pub diffuse: f32,
    /// Specular coefficient (0–1).
    pub specular: f32,
    /// Shininess exponent.
    pub shininess: f32,
    /// Light direction (unit vector).
    pub light_dir: [f32; 3],
}
impl SurfaceLighting {
    /// Apply Phong lighting to a color given a normal.
    pub fn apply(&self, color: PlotColor, normal: [f32; 3]) -> PlotColor {
        let n = normal;
        let l = self.light_dir;
        let ndotl = (n[0] * l[0] + n[1] * l[1] + n[2] * l[2]).clamp(0.0, 1.0);
        let intensity = self.ambient + self.diffuse * ndotl;
        PlotColor::rgb(
            (color.r * intensity).clamp(0.0, 1.0),
            (color.g * intensity).clamp(0.0, 1.0),
            (color.b * intensity).clamp(0.0, 1.0),
        )
    }
}
/// Options for publication-ready figure export.
#[derive(Debug, Clone)]
pub struct ExportOptions {
    /// Output format.
    pub format: ExportFormat,
    /// Dots per inch.
    pub dpi: u32,
    /// Whether to embed fonts.
    pub embed_fonts: bool,
    /// Whether to use a transparent background.
    pub transparent_background: bool,
    /// Crop whitespace from the output.
    pub tight_layout: bool,
}
/// Marker style for scatter points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MarkerStyle {
    /// Filled circle.
    #[default]
    Circle,
    /// Square.
    Square,
    /// Up-pointing triangle.
    TriangleUp,
    /// Down-pointing triangle.
    TriangleDown,
    /// Diamond (rotated square).
    Diamond,
    /// Plus sign.
    Plus,
    /// Cross (×).
    Cross,
    /// No marker.
    None,
}
/// A 2D line plot series (x values and corresponding y values).
#[derive(Debug, Clone)]
pub struct LineSeries {
    /// X data values.
    pub x: Vec<f64>,
    /// Y data values (must be the same length as `x`).
    pub y: Vec<f64>,
    /// Visual style for this series.
    pub style: LineStyle2D,
}
impl LineSeries {
    /// Create a series from x/y arrays.
    pub fn new(x: Vec<f64>, y: Vec<f64>, style: LineStyle2D) -> Self {
        debug_assert_eq!(x.len(), y.len(), "x and y must have equal length");
        Self { x, y, style }
    }
    /// Number of data points.
    pub fn len(&self) -> usize {
        self.x.len()
    }
    /// Return `true` if the series has no data.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
    /// Compute the y range `(min, max)` over the series.
    pub fn y_range(&self) -> (f64, f64) {
        let mn = self.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Compute the x range.
    pub fn x_range(&self) -> (f64, f64) {
        let mn = self.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Resample this series to `n` evenly-spaced x values using linear interpolation.
    pub fn resample(&self, n: usize) -> Self {
        if n == 0 || self.x.is_empty() {
            return Self::new(vec![], vec![], self.style.clone());
        }
        let (x0, x1) = self.x_range();
        let step = (x1 - x0) / (n - 1).max(1) as f64;
        let mut new_x = Vec::with_capacity(n);
        let mut new_y = Vec::with_capacity(n);
        for i in 0..n {
            let xi = x0 + i as f64 * step;
            new_x.push(xi);
            new_y.push(self.interpolate_at(xi));
        }
        Self::new(new_x, new_y, self.style.clone())
    }
    /// Linearly interpolate y at x value `xi` (clamped to data range).
    pub fn interpolate_at(&self, xi: f64) -> f64 {
        if self.x.len() < 2 {
            return self.y.first().cloned().unwrap_or(0.0);
        }
        let idx = self.x.partition_point(|&v| v <= xi);
        if idx == 0 {
            return self.y[0];
        }
        if idx >= self.x.len() {
            return *self.y.last().expect("collection should not be empty");
        }
        let x0 = self.x[idx - 1];
        let x1 = self.x[idx];
        let y0 = self.y[idx - 1];
        let y1 = self.y[idx];
        if (x1 - x0).abs() < 1e-15 {
            return y0;
        }
        y0 + (xi - x0) / (x1 - x0) * (y1 - y0)
    }
}
/// Represents one computed contour level.
#[derive(Debug, Clone)]
pub struct ContourLevel {
    /// The scalar value for this contour.
    pub value: f64,
    /// Color for this level.
    pub color: PlotColor,
    /// Polyline segments making up this contour (pairs of connected points).
    pub segments: Vec<([f64; 2], [f64; 2])>,
}
/// A 2D contour plot from a regular grid of scalar values.
#[derive(Debug, Clone)]
pub struct ContourPlot {
    /// Grid values as a flat row-major array of size `nx × ny`.
    pub values: Vec<f64>,
    /// Number of columns (x direction).
    pub nx: usize,
    /// Number of rows (y direction).
    pub ny: usize,
    /// X extent `[x_min, x_max]`.
    pub x_range: [f64; 2],
    /// Y extent `[y_min, y_max]`.
    pub y_range: [f64; 2],
    /// Contour levels to draw.
    pub levels: Vec<f64>,
    /// Whether to draw filled contour bands.
    pub filled: bool,
    /// Colormap for filled contours.
    pub colormap: PlotColormap,
    /// X axis configuration.
    pub x_axis: AxisConfig,
    /// Y axis configuration.
    pub y_axis: AxisConfig,
    /// Plot title.
    pub title: String,
    /// Whether to show a colorbar.
    pub show_colorbar: bool,
}
impl ContourPlot {
    /// Create a contour plot from a regular grid.
    ///
    /// `values` must have length `nx * ny`.
    pub fn from_grid(
        values: Vec<f64>,
        nx: usize,
        ny: usize,
        x_range: [f64; 2],
        y_range: [f64; 2],
        title: impl Into<String>,
    ) -> Self {
        assert_eq!(values.len(), nx * ny, "values must have length nx*ny");
        Self {
            values,
            nx,
            ny,
            x_range,
            y_range,
            levels: vec![],
            filled: true,
            colormap: PlotColormap::Viridis,
            x_axis: AxisConfig::default(),
            y_axis: AxisConfig::default(),
            title: title.into(),
            show_colorbar: true,
        }
    }
    /// Build a contour plot by evaluating a function on a regular grid.
    pub fn from_function<F>(
        f: F,
        nx: usize,
        ny: usize,
        x_range: [f64; 2],
        y_range: [f64; 2],
        title: impl Into<String>,
    ) -> Self
    where
        F: Fn(f64, f64) -> f64,
    {
        let dx = (x_range[1] - x_range[0]) / (nx - 1).max(1) as f64;
        let dy = (y_range[1] - y_range[0]) / (ny - 1).max(1) as f64;
        let mut values = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            let y = y_range[0] + iy as f64 * dy;
            for ix in 0..nx {
                let x = x_range[0] + ix as f64 * dx;
                values.push(f(x, y));
            }
        }
        Self::from_grid(values, nx, ny, x_range, y_range, title)
    }
    /// Set evenly-spaced contour levels from the data range.
    pub fn auto_levels(&mut self, n: usize) {
        let mn = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let step = (mx - mn) / (n + 1) as f64;
        self.levels = (1..=n).map(|i| mn + i as f64 * step).collect();
    }
    /// Get the grid value at `(ix, iy)`.
    pub fn value_at(&self, ix: usize, iy: usize) -> f64 {
        self.values[iy * self.nx + ix]
    }
    /// Compute the scalar data range.
    pub fn data_range(&self) -> (f64, f64) {
        let mn = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Perform simple marching-squares to extract polyline segments for one level.
    pub fn extract_level(&self, level: f64) -> ContourLevel {
        let mn = self.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = self
            .values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let t = if (mx - mn).abs() < 1e-14 {
            0.5
        } else {
            ((level - mn) / (mx - mn)).clamp(0.0, 1.0) as f32
        };
        let color = self.colormap.map(t);
        let dx = (self.x_range[1] - self.x_range[0]) / (self.nx - 1).max(1) as f64;
        let dy = (self.y_range[1] - self.y_range[0]) / (self.ny - 1).max(1) as f64;
        let mut segments = vec![];
        for iy in 0..self.ny.saturating_sub(1) {
            for ix in 0..self.nx.saturating_sub(1) {
                let v00 = self.value_at(ix, iy);
                let v10 = self.value_at(ix + 1, iy);
                let v01 = self.value_at(ix, iy + 1);
                let v11 = self.value_at(ix + 1, iy + 1);
                let x0 = self.x_range[0] + ix as f64 * dx;
                let y0 = self.y_range[0] + iy as f64 * dy;
                let x1 = x0 + dx;
                let y1 = y0 + dy;
                let above = [v00 >= level, v10 >= level, v01 >= level, v11 >= level];
                let case = (above[0] as u8)
                    | ((above[1] as u8) << 1)
                    | ((above[2] as u8) << 2)
                    | ((above[3] as u8) << 3);
                let lerp = |a: f64, b: f64| -> f64 {
                    if (b - a).abs() < 1e-15 {
                        0.5
                    } else {
                        ((level - a) / (b - a)).clamp(0.0, 1.0)
                    }
                };
                let e_bottom = [x0 + lerp(v00, v10) * dx, y0];
                let e_top = [x0 + lerp(v01, v11) * dx, y1];
                let e_left = [x0, y0 + lerp(v00, v01) * dy];
                let e_right = [x1, y0 + lerp(v10, v11) * dy];
                match case {
                    0 | 15 => {}
                    1 | 14 => segments.push((e_bottom, e_left)),
                    2 | 13 => segments.push((e_bottom, e_right)),
                    3 | 12 => segments.push((e_left, e_right)),
                    4 | 11 => segments.push((e_top, e_left)),
                    5 | 10 => {
                        segments.push((e_bottom, e_right));
                        segments.push((e_top, e_left));
                    }
                    6 | 9 => segments.push((e_bottom, e_top)),
                    7 | 8 => segments.push((e_top, e_right)),
                    _ => {}
                }
            }
        }
        ContourLevel {
            value: level,
            color,
            segments,
        }
    }
    /// Extract all configured levels.
    pub fn extract_all_levels(&self) -> Vec<ContourLevel> {
        self.levels.iter().map(|&l| self.extract_level(l)).collect()
    }
}
/// A 3D surface plot built from a function or data grid.
#[derive(Debug, Clone)]
pub struct SurfacePlot3D {
    /// Surface vertices (row-major: `ny × nx`).
    pub vertices: Vec<SurfaceVertex>,
    /// Triangle index triples.
    pub triangles: Vec<[usize; 3]>,
    /// Number of columns (x direction).
    pub nx: usize,
    /// Number of rows (y direction).
    pub ny: usize,
    /// Colormap for scalar values.
    pub colormap: PlotColormap,
    /// Lighting model.
    pub lighting: SurfaceLighting,
    /// Whether to draw a wireframe overlay.
    pub show_wireframe: bool,
    /// Plot title.
    pub title: String,
    /// Whether to perform hidden-line removal (conceptual flag for renderers).
    pub hidden_line_removal: bool,
}
impl SurfacePlot3D {
    /// Build a surface from a function `z = f(x, y)`.
    pub fn from_function<F>(
        f: F,
        nx: usize,
        ny: usize,
        x_range: [f64; 2],
        y_range: [f64; 2],
        title: impl Into<String>,
    ) -> Self
    where
        F: Fn(f64, f64) -> f64,
    {
        let dx = (x_range[1] - x_range[0]) / (nx - 1).max(1) as f64;
        let dy = (y_range[1] - y_range[0]) / (ny - 1).max(1) as f64;
        let mut vertices = Vec::with_capacity(nx * ny);
        let mut values: Vec<f64> = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            let y = y_range[0] + iy as f64 * dy;
            for ix in 0..nx {
                let x = x_range[0] + ix as f64 * dx;
                let z = f(x, y);
                values.push(z);
                vertices.push(SurfaceVertex {
                    position: [x, y, z],
                    normal: [0.0, 0.0, 1.0],
                    value: z,
                });
            }
        }
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                let zl = if ix > 0 {
                    values[iy * nx + ix - 1]
                } else {
                    values[idx]
                };
                let zr = if ix + 1 < nx {
                    values[iy * nx + ix + 1]
                } else {
                    values[idx]
                };
                let zd = if iy > 0 {
                    values[(iy - 1) * nx + ix]
                } else {
                    values[idx]
                };
                let zu = if iy + 1 < ny {
                    values[(iy + 1) * nx + ix]
                } else {
                    values[idx]
                };
                let dzdx = (zr - zl) / (2.0 * dx.max(1e-15));
                let dzdy = (zu - zd) / (2.0 * dy.max(1e-15));
                let len = (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
                vertices[idx].normal = [(-dzdx / len), (-dzdy / len), 1.0 / len];
            }
        }
        let mut triangles = Vec::with_capacity((nx - 1) * (ny - 1) * 2);
        for iy in 0..ny.saturating_sub(1) {
            for ix in 0..nx.saturating_sub(1) {
                let i00 = iy * nx + ix;
                let i10 = i00 + 1;
                let i01 = i00 + nx;
                let i11 = i01 + 1;
                triangles.push([i00, i10, i11]);
                triangles.push([i00, i11, i01]);
            }
        }
        Self {
            vertices,
            triangles,
            nx,
            ny,
            colormap: PlotColormap::Viridis,
            lighting: SurfaceLighting::default(),
            show_wireframe: false,
            title: title.into(),
            hidden_line_removal: true,
        }
    }
    /// Compute the z range of the surface.
    pub fn z_range(&self) -> (f64, f64) {
        let mn = self
            .vertices
            .iter()
            .map(|v| v.value)
            .fold(f64::INFINITY, f64::min);
        let mx = self
            .vertices
            .iter()
            .map(|v| v.value)
            .fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Apply the colormap to a vertex, returning its display color.
    pub fn color_for_vertex(&self, i: usize, z_min: f64, z_max: f64) -> PlotColor {
        let z = self.vertices[i].value;
        let range = z_max - z_min;
        let t = if range.abs() < 1e-14 {
            0.5
        } else {
            ((z - z_min) / range).clamp(0.0, 1.0) as f32
        };
        self.colormap.map(t)
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }
}
/// Tag identifying the type of plot held in a panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PlotType {
    /// 2D line plot.
    Line,
    /// Scatter plot.
    Scatter,
    /// Contour plot.
    Contour,
    /// 3D surface.
    Surface,
    /// Phase portrait.
    Phase,
    /// Empty panel (placeholder).
    Empty,
}
/// An RGBA color (each component in `[0.0, 1.0]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlotColor {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}
impl PlotColor {
    /// Create a fully-opaque color.
    pub fn rgb(r: f32, g: f32, b: f32) -> Self {
        Self { r, g, b, a: 1.0 }
    }
    /// Create a color with alpha.
    pub fn rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }
    /// Solid black.
    pub fn black() -> Self {
        Self::rgb(0.0, 0.0, 0.0)
    }
    /// Solid white.
    pub fn white() -> Self {
        Self::rgb(1.0, 1.0, 1.0)
    }
    /// Solid red.
    pub fn red() -> Self {
        Self::rgb(1.0, 0.0, 0.0)
    }
    /// Solid green.
    pub fn green() -> Self {
        Self::rgb(0.0, 0.8, 0.0)
    }
    /// Solid blue.
    pub fn blue() -> Self {
        Self::rgb(0.0, 0.0, 1.0)
    }
    /// Orange.
    pub fn orange() -> Self {
        Self::rgb(1.0, 0.5, 0.0)
    }
    /// Cyan.
    pub fn cyan() -> Self {
        Self::rgb(0.0, 1.0, 1.0)
    }
    /// Magenta.
    pub fn magenta() -> Self {
        Self::rgb(1.0, 0.0, 1.0)
    }
    /// Blend two colors at parameter `t ∈ [0, 1]`.
    pub fn lerp(a: PlotColor, b: PlotColor, t: f32) -> Self {
        Self {
            r: a.r + t * (b.r - a.r),
            g: a.g + t * (b.g - a.g),
            b: a.b + t * (b.b - a.b),
            a: a.a + t * (b.a - a.a),
        }
    }
}
/// Style for a single line series.
#[derive(Debug, Clone)]
pub struct LineStyle2D {
    /// Line color.
    pub color: PlotColor,
    /// Line pattern.
    pub style: LineStyle,
    /// Line width in points.
    pub width: f32,
    /// Marker drawn at each data point.
    pub marker: MarkerStyle,
    /// Marker size in points.
    pub marker_size: f32,
    /// Label for the legend.
    pub label: String,
}
/// Axis scale type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AxisScale {
    /// Linear scale (default).
    #[default]
    Linear,
    /// Base-10 logarithmic scale.
    Log10,
    /// Natural-logarithm scale.
    Ln,
    /// Square-root scale.
    Sqrt,
}
/// Error bar data for a scatter point.
#[derive(Debug, Clone, Copy)]
pub struct ErrorBar {
    /// Symmetric horizontal error (±).
    pub x_err: f64,
    /// Symmetric vertical error (±).
    pub y_err: f64,
}
/// A nullcline for a 2D phase portrait: a curve where dX/dt = 0 or dY/dt = 0.
#[derive(Debug, Clone)]
pub struct Nullcline {
    /// Points on the nullcline.
    pub points: Vec<[f64; 2]>,
    /// Color.
    pub color: PlotColor,
    /// Whether this is an x-nullcline (true) or y-nullcline (false).
    pub is_x_null: bool,
}
/// A single trajectory in phase space.
#[derive(Debug, Clone)]
pub struct PhaseTrajectory {
    /// State-space points (each entry is a state vector; 2D or 3D).
    pub points: Vec<Vec<f64>>,
    /// Display color.
    pub color: PlotColor,
    /// Line width.
    pub line_width: f32,
    /// Label.
    pub label: String,
    /// Whether to draw an arrowhead at the end.
    pub show_arrow: bool,
}
impl PhaseTrajectory {
    /// Create a 2D phase trajectory from x/y arrays.
    pub fn from_xy(x: Vec<f64>, y: Vec<f64>, color: PlotColor) -> Self {
        let points = x.into_iter().zip(y).map(|(xi, yi)| vec![xi, yi]).collect();
        Self {
            points,
            color,
            line_width: 1.5,
            label: String::new(),
            show_arrow: true,
        }
    }
    /// Create a 3D phase trajectory from x/y/z arrays.
    pub fn from_xyz(x: Vec<f64>, y: Vec<f64>, z: Vec<f64>, color: PlotColor) -> Self {
        let points = x
            .into_iter()
            .zip(y)
            .zip(z)
            .map(|((xi, yi), zi)| vec![xi, yi, zi])
            .collect();
        Self {
            points,
            color,
            line_width: 1.5,
            label: String::new(),
            show_arrow: true,
        }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Dimensionality (2 or 3, or 0 if empty).
    pub fn dim(&self) -> usize {
        self.points.first().map_or(0, |p| p.len())
    }
    /// Arc length of the trajectory.
    pub fn arc_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| {
                let a = &w[0];
                let b = &w[1];
                let d: f64 = a
                    .iter()
                    .zip(b.iter())
                    .map(|(ai, bi)| (bi - ai) * (bi - ai))
                    .sum();
                d.sqrt()
            })
            .sum()
    }
}
/// File format for figure export.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    /// Scalable Vector Graphics.
    Svg,
    /// Portable Network Graphics.
    Png,
    /// Portable Document Format.
    Pdf,
    /// Encapsulated PostScript.
    Eps,
}
/// A 2D scatter plot with optional colormap, error bars, and log scale.
#[derive(Debug, Clone)]
pub struct ScatterPlot {
    /// Data points.
    pub points: Vec<ScatterPoint>,
    /// Default marker style.
    pub marker: MarkerStyle,
    /// Default marker size.
    pub marker_size: f32,
    /// Default color (used when `color_value` is `None`).
    pub default_color: PlotColor,
    /// Colormap for scalar coloring.
    pub colormap: PlotColormap,
    /// Whether colormap coloring is enabled.
    pub use_colormap: bool,
    /// Scalar range `(min, max)` for colormap. `None` → auto.
    pub color_range: Option<(f64, f64)>,
    /// X axis configuration.
    pub x_axis: AxisConfig,
    /// Y axis configuration.
    pub y_axis: AxisConfig,
    /// Plot title.
    pub title: String,
    /// Whether to show a colorbar.
    pub show_colorbar: bool,
    /// Whether to draw error bars.
    pub show_error_bars: bool,
}
impl ScatterPlot {
    /// Create an empty scatter plot.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            points: vec![],
            marker: MarkerStyle::Circle,
            marker_size: 6.0,
            default_color: PlotColor::blue(),
            colormap: PlotColormap::Viridis,
            use_colormap: false,
            color_range: None,
            x_axis: AxisConfig::default(),
            y_axis: AxisConfig::default(),
            title: title.into(),
            show_colorbar: true,
            show_error_bars: true,
        }
    }
    /// Add a point.
    pub fn add_point(&mut self, p: ScatterPoint) {
        self.points.push(p);
    }
    /// Add raw (x, y) pairs.
    pub fn add_xy(&mut self, x: &[f64], y: &[f64]) {
        assert_eq!(x.len(), y.len());
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            self.points.push(ScatterPoint::new(xi, yi));
        }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if there are no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Resolve the color for point `i` given the effective color range.
    pub fn color_for(&self, i: usize, color_min: f64, color_max: f64) -> PlotColor {
        if !self.use_colormap {
            return self.default_color;
        }
        if let Some(cv) = self.points[i].color_value {
            let range = color_max - color_min;
            let t = if range.abs() < 1e-14 {
                0.5
            } else {
                ((cv - color_min) / range).clamp(0.0, 1.0) as f32
            };
            self.colormap.map(t)
        } else {
            self.default_color
        }
    }
    /// Compute the auto color range from all `color_value` fields.
    pub fn auto_color_range(&self) -> (f64, f64) {
        let vals: Vec<f64> = self.points.iter().filter_map(|p| p.color_value).collect();
        if vals.is_empty() {
            return (0.0, 1.0);
        }
        let mn = vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Compute the x data range.
    pub fn x_range(&self) -> (f64, f64) {
        let xs: Vec<f64> = self.points.iter().map(|p| p.x).collect();
        let mn = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
    /// Compute the y data range.
    pub fn y_range(&self) -> (f64, f64) {
        let ys: Vec<f64> = self.points.iter().map(|p| p.y).collect();
        let mn = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (mn, mx)
    }
}
/// A complete line plot with multiple series, axis labels, and title.
#[derive(Debug, Clone, Default)]
pub struct LinePlot {
    /// Series to draw.
    pub series: Vec<LineSeries>,
    /// Horizontal axis configuration.
    pub x_axis: AxisConfig,
    /// Vertical axis configuration.
    pub y_axis: AxisConfig,
    /// Plot title.
    pub title: String,
    /// Whether to show a legend.
    pub show_legend: bool,
    /// Width in logical pixels.
    pub width: u32,
    /// Height in logical pixels.
    pub height: u32,
    /// Background color.
    pub background: PlotColor,
}
impl LinePlot {
    /// Create a new empty line plot.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            show_legend: true,
            width: 800,
            height: 600,
            background: PlotColor::white(),
            ..Default::default()
        }
    }
    /// Add a series to the plot.
    pub fn add_series(&mut self, series: LineSeries) {
        self.series.push(series);
    }
    /// Number of series.
    pub fn num_series(&self) -> usize {
        self.series.len()
    }
    /// Compute the overall x range from all series data.
    pub fn auto_x_range(&self) -> (f64, f64) {
        let all_x: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.x.iter().cloned())
            .collect();
        self.x_axis.effective_range(&all_x)
    }
    /// Compute the overall y range from all series data.
    pub fn auto_y_range(&self) -> (f64, f64) {
        let all_y: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.y.iter().cloned())
            .collect();
        self.y_axis.effective_range(&all_y)
    }
    /// Add a series from raw (x, y) pairs.
    pub fn add_xy(&mut self, x: Vec<f64>, y: Vec<f64>, color: PlotColor, label: impl Into<String>) {
        let style = LineStyle2D {
            color,
            label: label.into(),
            ..Default::default()
        };
        self.add_series(LineSeries::new(x, y, style));
    }
    /// Set axis labels.
    pub fn set_labels(&mut self, x_label: impl Into<String>, y_label: impl Into<String>) {
        self.x_axis.label = x_label.into();
        self.y_axis.label = y_label.into();
    }
}
/// A single scatter data point.
#[derive(Debug, Clone)]
pub struct ScatterPoint {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Scalar value used for colormap coloring (optional; `None` → default color).
    pub color_value: Option<f64>,
    /// Marker size override (optional; `None` → default size).
    pub size: Option<f32>,
    /// Optional error bars.
    pub error_bar: Option<ErrorBar>,
}
impl ScatterPoint {
    /// Create a basic point without error bars.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            x,
            y,
            color_value: None,
            size: None,
            error_bar: None,
        }
    }
    /// Create a point with a colormap scalar.
    pub fn with_color_value(mut self, v: f64) -> Self {
        self.color_value = Some(v);
        self
    }
    /// Create a point with error bars.
    pub fn with_error(mut self, x_err: f64, y_err: f64) -> Self {
        self.error_bar = Some(ErrorBar { x_err, y_err });
        self
    }
    /// Create a point with custom size.
    pub fn with_size(mut self, size: f32) -> Self {
        self.size = Some(size);
        self
    }
}
/// Minimal SVG exporter for a `LinePlot` (renders line geometry only).
pub struct SvgExporter;
impl SvgExporter {
    /// Export a `LinePlot` as a minimal SVG string.
    pub fn export_line_plot(plot: &LinePlot) -> String {
        let w = plot.width;
        let h = plot.height;
        let (xmin, xmax) = plot.auto_x_range();
        let (ymin, ymax) = plot.auto_y_range();
        let xrange = (xmax - xmin).max(1e-14);
        let yrange = (ymax - ymin).max(1e-14);
        let to_sx = |x: f64| -> f32 { ((x - xmin) / xrange * (w as f64 - 80.0) + 40.0) as f32 };
        let to_sy = |y: f64| -> f32 {
            (h as f32) - (((y - ymin) / yrange * (h as f64 - 80.0) + 40.0) as f32)
        };
        let mut svg =
            format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\">\n");
        let bg = plot.background;
        svg.push_str(&format!(
            "  <rect width=\"{w}\" height=\"{h}\" fill=\"rgb({},{},{})\" />\n",
            (bg.r * 255.0) as u8,
            (bg.g * 255.0) as u8,
            (bg.b * 255.0) as u8
        ));
        if !plot.title.is_empty() {
            svg.push_str(&format!(
                "  <text x=\"{}\" y=\"20\" font-size=\"14\" text-anchor=\"middle\">{}</text>\n",
                w / 2,
                plot.title
            ));
        }
        for series in &plot.series {
            if series.len() < 2 {
                continue;
            }
            let c = series.style.color;
            let stroke = format!(
                "rgb({},{},{})",
                (c.r * 255.0) as u8,
                (c.g * 255.0) as u8,
                (c.b * 255.0) as u8
            );
            svg.push_str(&format!(
                "  <polyline stroke=\"{}\" stroke-width=\"{}\" fill=\"none\" points=\"",
                stroke, series.style.width
            ));
            for i in 0..series.len() {
                svg.push_str(&format!(
                    "{:.1},{:.1} ",
                    to_sx(series.x[i]),
                    to_sy(series.y[i])
                ));
            }
            svg.push_str("\" />\n");
        }
        svg.push_str("</svg>\n");
        svg
    }
}
/// 2D vector field function signature: `(x, y) → (dx, dy)`.
pub type VectorField2dFn = fn(f64, f64) -> (f64, f64);

/// An equilibrium point for a dynamical system.
#[derive(Debug, Clone)]
pub struct EquilibriumPoint {
    /// Location in phase space.
    pub position: Vec<f64>,
    /// Stability type description (e.g. "stable node").
    pub stability: String,
    /// Display color.
    pub color: PlotColor,
    /// Marker size.
    pub size: f32,
}
/// 2D or 3D phase portrait: trajectories, nullclines, and equilibrium points.
#[derive(Debug, Clone)]
pub struct PhasePortrait {
    /// Phase trajectories.
    pub trajectories: Vec<PhaseTrajectory>,
    /// Nullclines (2D only).
    pub nullclines: Vec<Nullcline>,
    /// Equilibrium points.
    pub equilibria: Vec<EquilibriumPoint>,
    /// Axis labels.
    pub axis_labels: Vec<String>,
    /// Plot title.
    pub title: String,
    /// Whether to draw a vector field background.
    pub show_vector_field: bool,
    /// Resolution of the background vector field grid.
    pub vector_field_resolution: usize,
    /// 2D vector field function: (x, y) → (dx, dy).
    pub vector_field_2d: Option<VectorField2dFn>,
}
impl PhasePortrait {
    /// Create an empty phase portrait.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            trajectories: vec![],
            nullclines: vec![],
            equilibria: vec![],
            axis_labels: vec!["x".to_string(), "y".to_string()],
            title: title.into(),
            show_vector_field: false,
            vector_field_resolution: 20,
            vector_field_2d: None,
        }
    }
    /// Add a phase trajectory.
    pub fn add_trajectory(&mut self, t: PhaseTrajectory) {
        self.trajectories.push(t);
    }
    /// Add a nullcline.
    pub fn add_nullcline(&mut self, nc: Nullcline) {
        self.nullclines.push(nc);
    }
    /// Add an equilibrium point.
    pub fn add_equilibrium(&mut self, eq: EquilibriumPoint) {
        self.equilibria.push(eq);
    }
    /// Number of trajectories.
    pub fn num_trajectories(&self) -> usize {
        self.trajectories.len()
    }
    /// Integrate a 2D system `(dx/dt, dy/dt) = f(x, y)` from an initial condition.
    pub fn integrate_2d<F>(
        f: F,
        x0: f64,
        y0: f64,
        dt: f64,
        steps: usize,
        color: PlotColor,
    ) -> PhaseTrajectory
    where
        F: Fn(f64, f64) -> (f64, f64),
    {
        let mut xs = Vec::with_capacity(steps + 1);
        let mut ys = Vec::with_capacity(steps + 1);
        xs.push(x0);
        ys.push(y0);
        let mut x = x0;
        let mut y = y0;
        for _ in 0..steps {
            let (dx, dy) = f(x, y);
            x += dt * dx;
            y += dt * dy;
            xs.push(x);
            ys.push(y);
        }
        PhaseTrajectory::from_xy(xs, ys, color)
    }
    /// Integrate a 3D system from an initial condition.
    pub fn integrate_3d<F>(
        f: F,
        x0: f64,
        y0: f64,
        z0: f64,
        dt: f64,
        steps: usize,
        color: PlotColor,
    ) -> PhaseTrajectory
    where
        F: Fn(f64, f64, f64) -> (f64, f64, f64),
    {
        let mut xs = Vec::with_capacity(steps + 1);
        let mut ys = Vec::with_capacity(steps + 1);
        let mut zs = Vec::with_capacity(steps + 1);
        xs.push(x0);
        ys.push(y0);
        zs.push(z0);
        let mut x = x0;
        let mut y = y0;
        let mut z = z0;
        for _ in 0..steps {
            let (dx, dy, dz) = f(x, y, z);
            x += dt * dx;
            y += dt * dy;
            z += dt * dz;
            xs.push(x);
            ys.push(y);
            zs.push(z);
        }
        PhaseTrajectory::from_xyz(xs, ys, zs, color)
    }
    /// Compute x-nullcline by scanning for sign changes in `dx/dt`.
    pub fn compute_nullcline_x<F>(
        f: F,
        x_range: [f64; 2],
        y_range: [f64; 2],
        resolution: usize,
    ) -> Nullcline
    where
        F: Fn(f64, f64) -> (f64, f64),
    {
        let mut points = vec![];
        let dx = (x_range[1] - x_range[0]) / (resolution - 1).max(1) as f64;
        let dy_step = (y_range[1] - y_range[0]) / (resolution - 1).max(1) as f64;
        for iy in 0..resolution {
            let y = y_range[0] + iy as f64 * dy_step;
            for ix in 0..resolution.saturating_sub(1) {
                let x0 = x_range[0] + ix as f64 * dx;
                let x1 = x0 + dx;
                let (fx0, _) = f(x0, y);
                let (fx1, _) = f(x1, y);
                if fx0 * fx1 <= 0.0 {
                    let t = if (fx1 - fx0).abs() < 1e-14 {
                        0.5
                    } else {
                        -fx0 / (fx1 - fx0)
                    };
                    points.push([x0 + t * dx, y]);
                }
            }
        }
        Nullcline {
            points,
            color: PlotColor::red(),
            is_x_null: true,
        }
    }
}
/// A single panel within a multi-panel figure.
#[derive(Debug, Clone)]
pub struct PlotPanel {
    /// Row index (0-based).
    pub row: usize,
    /// Column index (0-based).
    pub col: usize,
    /// Optional plot title (overrides child title if set).
    pub title: Option<String>,
    /// Whether this panel shares its x axis with the panel above.
    pub shared_x: bool,
    /// Whether this panel shares its y axis with the panel to the left.
    pub shared_y: bool,
    /// Type tag for what kind of plot this panel holds.
    pub plot_type: PlotType,
}
/// Layout manager for a multi-panel scientific figure.
///
/// Panels are arranged in a grid of `rows × cols`.
#[derive(Debug, Clone, Default)]
pub struct PlotLayout {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Figure-level title.
    pub title: String,
    /// Panels in this layout.
    pub panels: Vec<PlotPanel>,
    /// Per-panel axis ranges: key = `(row, col)`, value = `((x0,x1),(y0,y1))`.
    pub axis_ranges: PanelAxisRanges,
    /// Figure width in points/pixels.
    pub width: u32,
    /// Figure height in points/pixels.
    pub height: u32,
    /// Horizontal spacing between subplots (fraction of panel width).
    pub hspace: f32,
    /// Vertical spacing between subplots (fraction of panel height).
    pub vspace: f32,
}
impl PlotLayout {
    /// Create a new layout with the given grid dimensions.
    pub fn new(rows: usize, cols: usize, title: impl Into<String>) -> Self {
        Self {
            rows,
            cols,
            title: title.into(),
            panels: vec![],
            axis_ranges: HashMap::new(),
            width: 1200,
            height: 900,
            hspace: 0.05,
            vspace: 0.05,
        }
    }
    /// Add a panel at `(row, col)`.
    pub fn add_panel(&mut self, row: usize, col: usize, plot_type: PlotType) {
        self.panels.push(PlotPanel {
            row,
            col,
            title: None,
            shared_x: false,
            shared_y: false,
            plot_type,
        });
    }
    /// Set the axis range for panel `(row, col)`.
    pub fn set_axis_range(
        &mut self,
        row: usize,
        col: usize,
        x_range: (f64, f64),
        y_range: (f64, f64),
    ) {
        self.axis_ranges.insert((row, col), (x_range, y_range));
    }
    /// Mark a panel as sharing the x axis with the panel above it.
    pub fn share_x_axis(&mut self, row: usize, col: usize) {
        if let Some(p) = self
            .panels
            .iter_mut()
            .find(|p| p.row == row && p.col == col)
        {
            p.shared_x = true;
        }
    }
    /// Mark a panel as sharing the y axis with the panel to its left.
    pub fn share_y_axis(&mut self, row: usize, col: usize) {
        if let Some(p) = self
            .panels
            .iter_mut()
            .find(|p| p.row == row && p.col == col)
        {
            p.shared_y = true;
        }
    }
    /// Number of panels.
    pub fn num_panels(&self) -> usize {
        self.panels.len()
    }
    /// Compute the normalized pixel bounds `[left, top, right, bottom]` for panel `(row, col)`.
    ///
    /// All values are in `[0.0, 1.0]` (fraction of figure size).
    pub fn panel_bounds(&self, row: usize, col: usize) -> [f32; 4] {
        let pw = 1.0 / self.cols as f32 - self.hspace;
        let ph = 1.0 / self.rows as f32 - self.vspace;
        let x0 = col as f32 / self.cols as f32 + self.hspace * 0.5;
        let y0 = row as f32 / self.rows as f32 + self.vspace * 0.5;
        [x0, y0, x0 + pw, y0 + ph]
    }
    /// Generate an SVG-like layout description string (for debugging/export).
    pub fn describe(&self) -> String {
        let mut s = format!(
            "PlotLayout '{}' {}×{} panels ({}×{} figure)\n",
            self.title, self.rows, self.cols, self.width, self.height
        );
        for p in &self.panels {
            let bounds = self.panel_bounds(p.row, p.col);
            s.push_str(&format!(
                "  [{},{}] {:?} bounds=[{:.2},{:.2},{:.2},{:.2}]\n",
                p.row, p.col, p.plot_type, bounds[0], bounds[1], bounds[2], bounds[3]
            ));
        }
        s
    }
    /// Fill all grid positions with a default `Empty` panel.
    pub fn fill_empty(&mut self) {
        let existing: std::collections::HashSet<(usize, usize)> =
            self.panels.iter().map(|p| (p.row, p.col)).collect();
        for r in 0..self.rows {
            for c in 0..self.cols {
                if !existing.contains(&(r, c)) {
                    self.add_panel(r, c, PlotType::Empty);
                }
            }
        }
    }
}
/// Line style for plots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LineStyle {
    /// Solid line.
    #[default]
    Solid,
    /// Dashed line (equal dash/gap).
    Dashed,
    /// Dotted line (small dots).
    Dotted,
    /// Dash-dot alternating pattern.
    DashDot,
    /// No line (invisible).
    None,
}
