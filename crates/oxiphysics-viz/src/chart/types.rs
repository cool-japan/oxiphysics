//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

/// Residual convergence chart (semi-log).
#[allow(dead_code)]
pub struct ConvergencePlot {
    /// Iteration indices.
    pub iterations: Vec<usize>,
    /// Residual values.
    pub residuals: Vec<f64>,
    /// Convergence rate (estimated geometric ratio).
    pub convergence_rate: f64,
    /// Target residual.
    pub target: f64,
    /// Chart title.
    pub title: String,
    /// Multiple residual channels (name, data).
    pub channels: Vec<(String, Vec<f64>)>,
}
impl ConvergencePlot {
    /// Create a convergence plot.
    pub fn new(title: &str, target: f64) -> Self {
        Self {
            iterations: Vec::new(),
            residuals: Vec::new(),
            convergence_rate: 0.0,
            target,
            title: title.to_string(),
            channels: Vec::new(),
        }
    }
    /// Append a residual value.
    pub fn push(&mut self, iteration: usize, residual: f64) {
        self.iterations.push(iteration);
        self.residuals.push(residual);
        self.update_rate();
    }
    /// Update the estimated convergence rate.
    fn update_rate(&mut self) {
        let n = self.residuals.len();
        if n < 2 {
            self.convergence_rate = 0.0;
            return;
        }
        let r_prev = self.residuals[n - 2];
        let r_curr = self.residuals[n - 1];
        if r_prev > 1e-15 {
            self.convergence_rate = r_curr / r_prev;
        }
    }
    /// Returns `true` when the residual is below the target.
    pub fn converged(&self) -> bool {
        self.residuals
            .last()
            .map(|&r| r < self.target)
            .unwrap_or(false)
    }
    /// Estimate how many more iterations until convergence.
    pub fn estimated_iters_to_convergence(&self) -> Option<usize> {
        let last = *self.residuals.last()?;
        if self.convergence_rate >= 1.0 || self.convergence_rate <= 0.0 {
            return None;
        }
        if last <= self.target {
            return Some(0);
        }
        let ratio = (self.target / last).ln() / self.convergence_rate.ln();
        Some(ratio.ceil() as usize)
    }
    /// Add a named residual channel.
    pub fn add_channel(&mut self, name: &str, residuals: Vec<f64>) {
        self.channels.push((name.to_string(), residuals));
    }
}
/// Multi-series line chart with configurable legend and smooth-curve overlay.
#[allow(dead_code)]
pub struct MultiLinePlot {
    /// Series data.
    pub series: Vec<Series>,
    /// X axis configuration.
    pub x_axis: AxisConfig,
    /// Y axis configuration.
    pub y_axis: AxisConfig,
    /// Chart title.
    pub title: String,
    /// Show legend flag.
    pub show_legend: bool,
}
impl MultiLinePlot {
    /// Creates an empty multi-series line chart.
    pub fn new(title: &str) -> Self {
        Self {
            series: Vec::new(),
            x_axis: AxisConfig::new("x"),
            y_axis: AxisConfig::new("y"),
            title: title.to_string(),
            show_legend: true,
        }
    }
    /// Adds a data series to the chart.
    pub fn add_series(&mut self, s: Series) {
        self.series.push(s);
    }
    /// Computes the overall x range across all series.
    pub fn x_range(&self) -> (f64, f64) {
        let all_x: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.x.iter().cloned())
            .collect();
        if all_x.is_empty() {
            return (0.0, 1.0);
        }
        self.x_axis.effective_range(&all_x)
    }
    /// Computes the overall y range across all series.
    pub fn y_range(&self) -> (f64, f64) {
        let all_y: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.y.iter().cloned())
            .collect();
        if all_y.is_empty() {
            return (0.0, 1.0);
        }
        self.y_axis.effective_range(&all_y)
    }
    /// Returns legend entries as (label, color) pairs.
    pub fn legend_entries(&self) -> Vec<(&str, [f32; 4])> {
        self.series
            .iter()
            .map(|s| (s.label.as_str(), s.color))
            .collect()
    }
    /// Applies moving-average smoothing to series `idx` and returns the result.
    pub fn smoothed_series(&self, idx: usize, window: usize) -> Vec<f64> {
        if idx >= self.series.len() {
            return Vec::new();
        }
        moving_average(&self.series[idx].y, window)
    }
    /// Exports the chart as a minimal SVG string.
    pub fn to_svg(&self, width: usize, height: usize) -> String {
        let exp = SvgExporter::new(width, height);
        let mut out = exp.export_header();
        let (x0, x1) = self.x_range();
        let (y0, y1) = self.y_range();
        for s in &self.series {
            let color_str = format!(
                "rgba({},{},{},{})",
                (s.color[0] * 255.0) as u8,
                (s.color[1] * 255.0) as u8,
                (s.color[2] * 255.0) as u8,
                s.color[3]
            );
            let pts: Vec<String> =
                s.x.iter()
                    .zip(s.y.iter())
                    .map(|(&xi, &yi)| {
                        let px = normalize(xi, x0, x1) * width as f64;
                        let py = (1.0 - normalize(yi, y0, y1)) * height as f64;
                        format!("{:.1},{:.1}", px, py)
                    })
                    .collect();
            if !pts.is_empty() {
                out.push_str(&format!(
                    "<polyline points=\"{}\" stroke=\"{}\" fill=\"none\" stroke-width=\"1.5\"/>",
                    pts.join(" "),
                    color_str
                ));
            }
        }
        out.push_str("</svg>");
        out
    }
}
/// A single scatter point.
#[allow(dead_code)]
pub struct ScatterPoint {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Scalar value for colour mapping.
    pub value: f64,
    /// Marker size (pixels).
    pub size: f32,
}
/// Error bar data for scatter or line charts.
#[allow(dead_code)]
pub struct ErrorBars {
    /// X positions.
    pub x: Vec<f64>,
    /// Y centre values.
    pub y: Vec<f64>,
    /// Error half-widths (symmetric ±).
    pub errors: Vec<f64>,
    /// Cap width in display units.
    pub cap_width: f64,
    /// Bar color.
    pub color: [f32; 4],
}
impl ErrorBars {
    /// Creates new error bars.
    pub fn new(x: Vec<f64>, y: Vec<f64>, errors: Vec<f64>) -> Self {
        Self {
            x,
            y,
            errors,
            cap_width: 4.0,
            color: [0.1, 0.1, 0.1, 1.0],
        }
    }
    /// Returns the number of error bars.
    pub fn len(&self) -> usize {
        self.x.len().min(self.y.len()).min(self.errors.len())
    }
    /// Returns `true` if there are no error bars.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Returns (lower, upper) bounds for error bar at index `i`.
    pub fn error_range(&self, i: usize) -> (f64, f64) {
        let y = self.y[i];
        let e = self.errors[i];
        (y - e, y + e)
    }
    /// Generates SVG elements for all error bars.
    pub fn to_svg_elements(
        &self,
        x_lo: f64,
        x_hi: f64,
        y_lo: f64,
        y_hi: f64,
        w: f64,
        h: f64,
    ) -> Vec<String> {
        let mut elems = Vec::new();
        for i in 0..self.len() {
            let px = normalize(self.x[i], x_lo, x_hi) * w;
            let py_c = (1.0 - normalize(self.y[i], y_lo, y_hi)) * h;
            let py_lo = (1.0 - normalize(self.y[i] - self.errors[i], y_lo, y_hi)) * h;
            let py_hi = (1.0 - normalize(self.y[i] + self.errors[i], y_lo, y_hi)) * h;
            let cw = self.cap_width * 0.5;
            let r = (self.color[0] * 255.0) as u8;
            let g_c = (self.color[1] * 255.0) as u8;
            let b = (self.color[2] * 255.0) as u8;
            let col = format!("rgb({},{},{})", r, g_c, b);
            elems
                .push(
                    format!(
                        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\"/>",
                        px, py_lo, px, py_hi, col
                    ),
                );
            elems
                .push(
                    format!(
                        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\"/>",
                        px - cw, py_lo, px + cw, py_lo, col
                    ),
                );
            elems
                .push(
                    format!(
                        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\"/>",
                        px - cw, py_hi, px + cw, py_hi, col
                    ),
                );
            let _ = py_c;
        }
        elems
    }
}
/// A single bubble in a 3D bubble chart.
#[allow(dead_code)]
pub struct BubblePoint {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z value → bubble size.
    pub z: f64,
    /// Optional label.
    pub label: Option<String>,
    /// Colour.
    pub color: [f32; 4],
}
/// A text annotation placed at a specific data-coordinate position.
#[allow(dead_code)]
pub struct Annotation {
    /// Annotation text.
    pub text: String,
    /// X position in data coordinates.
    pub x: f64,
    /// Y position in data coordinates.
    pub y: f64,
    /// Font size in points.
    pub font_size: usize,
    /// Text color as CSS string.
    pub color: String,
    /// Horizontal alignment: "left", "center", or "right".
    pub align: String,
}
impl Annotation {
    /// Creates an annotation at (x, y) with default style.
    pub fn new(text: &str, x: f64, y: f64) -> Self {
        Self {
            text: text.to_string(),
            x,
            y,
            font_size: 12,
            color: "black".to_string(),
            align: "center".to_string(),
        }
    }
    /// Sets the font size and returns `self` for chaining.
    pub fn with_font_size(mut self, size: usize) -> Self {
        self.font_size = size;
        self
    }
    /// Sets the color and returns `self` for chaining.
    pub fn with_color(mut self, color: &str) -> Self {
        self.color = color.to_string();
        self
    }
    /// Renders the annotation as an SVG ``text` element.
    ///
    /// `x_px` and `y_px` are the pixel coordinates of the annotation.
    pub fn to_svg_at(&self, x_px: f64, y_px: f64) -> String {
        format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{}\" fill=\"{}\" text-anchor=\"{}\">{}</text>",
            x_px, y_px, self.font_size, self.color, self.align, self.text
        )
    }
    /// Renders the annotation using its own (x, y) as pixel coords.
    pub fn to_svg(&self) -> String {
        self.to_svg_at(self.x, self.y)
    }
}
/// Polar plot (rose diagram) for directional or cyclic data.
#[allow(dead_code)]
pub struct PolarPlot {
    /// (angle_rad, radius) pairs.
    pub points: Vec<(f64, f64)>,
    /// Chart title.
    pub title: String,
    /// Maximum radius for scaling.
    pub r_max: f64,
}
impl PolarPlot {
    /// Creates an empty polar plot.
    pub fn new(title: &str) -> Self {
        Self {
            points: Vec::new(),
            title: title.to_string(),
            r_max: 1.0,
        }
    }
    /// Adds a (angle, radius) data point.
    pub fn add_point(&mut self, angle_rad: f64, r: f64) {
        self.points.push((angle_rad, r));
        if r > self.r_max {
            self.r_max = r;
        }
    }
    /// Converts a polar point to Cartesian (x, y) in [-1, 1]².
    pub fn to_cartesian(&self, angle_rad: f64, r: f64) -> (f64, f64) {
        let rn = r / self.r_max.max(1e-15);
        (rn * angle_rad.cos(), rn * angle_rad.sin())
    }
    /// Determines the rose-diagram bin index for a given angle and `n_bins`.
    pub fn rose_bin(&self, angle_rad: f64, n_bins: usize) -> usize {
        use std::f64::consts::TAU;
        let a = angle_rad.rem_euclid(TAU);
        let bin = (a / TAU * n_bins as f64) as usize;
        bin.min(n_bins - 1)
    }
    /// Builds a rose histogram with `n_bins` angular bins.
    ///
    /// Returns bin counts.
    pub fn rose_histogram(&self, n_bins: usize) -> Vec<usize> {
        let mut bins = vec![0usize; n_bins];
        for &(angle, _r) in &self.points {
            bins[self.rose_bin(angle, n_bins)] += 1;
        }
        bins
    }
}
/// Contour plot using the marching squares algorithm.
///
/// Extracts iso-contour line segments from a 2D scalar field.
#[allow(dead_code)]
pub struct ContourPlot {
    /// Grid width (columns).
    pub width: usize,
    /// Grid height (rows).
    pub height: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Scalar field data (row-major, size = width * height).
    pub data: Vec<f64>,
    /// Chart title.
    pub title: String,
}
impl ContourPlot {
    /// Creates a new contour plot.
    pub fn new(title: &str, data: Vec<f64>, width: usize, height: usize, dx: f64) -> Self {
        Self {
            width,
            height,
            dx,
            data,
            title: title.to_string(),
        }
    }
    /// Generates `n` evenly-spaced threshold values spanning the data range.
    pub fn auto_thresholds(&self, n: usize) -> Vec<f64> {
        let lo = slice_min(&self.data);
        let hi = slice_max(&self.data);
        if n == 0 {
            return Vec::new();
        }
        (0..n)
            .map(|i| lo + (hi - lo) * i as f64 / (n - 1).max(1) as f64)
            .collect()
    }
    /// Linearly interpolates along an edge to find the crossing point.
    fn interp(v0: f64, v1: f64, iso: f64, p0: f64, p1: f64) -> f64 {
        if (v1 - v0).abs() < 1e-15 {
            return (p0 + p1) * 0.5;
        }
        p0 + (iso - v0) / (v1 - v0) * (p1 - p0)
    }
    /// Runs marching squares for a single iso-value, returning line segments.
    ///
    /// Each segment is `\[(x0,y0), (x1,y1)\]` in data coordinates.
    pub fn marching_squares(&self, iso: f64) -> Vec<[(f64, f64); 2]> {
        let w = self.width;
        let h = self.height;
        let dx = self.dx;
        let mut segments = Vec::new();
        for row in 0..h - 1 {
            for col in 0..w - 1 {
                let v00 = self.data[row * w + col];
                let v10 = self.data[row * w + col + 1];
                let v01 = self.data[(row + 1) * w + col];
                let v11 = self.data[(row + 1) * w + col + 1];
                let x0 = col as f64 * dx;
                let x1 = (col + 1) as f64 * dx;
                let y0 = row as f64 * dx;
                let y1 = (row + 1) as f64 * dx;
                let code = ((v00 >= iso) as u8)
                    | (((v10 >= iso) as u8) << 1)
                    | (((v01 >= iso) as u8) << 2)
                    | (((v11 >= iso) as u8) << 3);
                let bottom = (Self::interp(v00, v10, iso, x0, x1), y0);
                let top = (Self::interp(v01, v11, iso, x0, x1), y1);
                let left = (x0, Self::interp(v00, v01, iso, y0, y1));
                let right = (x1, Self::interp(v10, v11, iso, y0, y1));
                match code {
                    1 | 14 => segments.push([left, bottom]),
                    2 | 13 => segments.push([bottom, right]),
                    3 | 12 => segments.push([left, right]),
                    4 | 11 => segments.push([left, top]),
                    6 | 9 => segments.push([bottom, top]),
                    7 | 8 => segments.push([top, right]),
                    5 => {
                        segments.push([left, bottom]);
                        segments.push([top, right]);
                    }
                    10 => {
                        segments.push([bottom, right]);
                        segments.push([left, top]);
                    }
                    _ => {}
                }
            }
        }
        segments
    }
}
/// Vehicle telemetry chart.
#[allow(dead_code)]
pub struct SensorPlot {
    /// Telemetry samples.
    pub samples: Vec<TelemetrySample>,
    /// Chart title.
    pub title: String,
    /// Visible channels.
    pub channels: Vec<String>,
}
impl SensorPlot {
    /// Create a sensor plot.
    pub fn new(title: &str) -> Self {
        Self {
            samples: Vec::new(),
            title: title.to_string(),
            channels: vec![
                "speed".to_string(),
                "g_lat".to_string(),
                "g_long".to_string(),
            ],
        }
    }
    /// Append a telemetry sample.
    pub fn push(&mut self, sample: TelemetrySample) {
        self.samples.push(sample);
    }
    /// Compute maximum lateral g-force over the session.
    pub fn max_lateral_g(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.g_lateral.abs())
            .fold(0.0, f64::max)
    }
    /// Compute maximum speed over the session (m/s).
    pub fn max_speed(&self) -> f64 {
        self.samples.iter().map(|s| s.speed).fold(0.0, f64::max)
    }
    /// Extract speed time series.
    pub fn speed_series(&self) -> (Vec<f64>, Vec<f64>) {
        let t: Vec<f64> = self.samples.iter().map(|s| s.time).collect();
        let v: Vec<f64> = self.samples.iter().map(|s| s.speed).collect();
        (t, v)
    }
}
/// Logarithmic-scale axis configuration.
#[allow(dead_code)]
pub struct LogAxis {
    /// Lower bound (must be > 0).
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// Number of ticks.
    pub n_ticks: usize,
}
impl LogAxis {
    /// Creates a logarithmic axis.
    pub fn new(lo: f64, hi: f64, n_ticks: usize) -> Self {
        Self {
            lo: lo.max(1e-15),
            hi,
            n_ticks,
        }
    }
    /// Generates `n_ticks` evenly-spaced positions in log space.
    pub fn tick_positions(&self) -> Vec<f64> {
        let log_lo = self.lo.log10();
        let log_hi = self.hi.log10();
        (0..self.n_ticks)
            .map(|i| {
                let t = if self.n_ticks > 1 {
                    i as f64 / (self.n_ticks - 1) as f64
                } else {
                    0.0
                };
                10.0_f64.powf(log_lo + t * (log_hi - log_lo))
            })
            .collect()
    }
    /// Maps a value to normalised [0, 1] in log space.
    pub fn normalize_log(&self, v: f64) -> f64 {
        let log_lo = self.lo.log10();
        let log_hi = self.hi.log10();
        let log_v = v.max(self.lo).log10();
        ((log_v - log_lo) / (log_hi - log_lo)).clamp(0.0, 1.0)
    }
}
/// 2D phase portrait: x vs dx/dt.
#[allow(dead_code)]
pub struct PhasePortrait {
    /// Trajectory x values.
    pub x: Vec<f64>,
    /// Trajectory dx/dt values.
    pub dxdt: Vec<f64>,
    /// Fixed points (equilibria).
    pub fixed_points: Vec<[f64; 2]>,
    /// Separatrix points (saddle connections).
    pub separatrix: Vec<Vec<[f64; 2]>>,
    /// Chart title.
    pub title: String,
}
impl PhasePortrait {
    /// Create a phase portrait.
    pub fn new(title: &str) -> Self {
        Self {
            x: Vec::new(),
            dxdt: Vec::new(),
            fixed_points: Vec::new(),
            separatrix: Vec::new(),
            title: title.to_string(),
        }
    }
    /// Add a trajectory point.
    pub fn add_point(&mut self, x: f64, dxdt: f64) {
        self.x.push(x);
        self.dxdt.push(dxdt);
    }
    /// Add a fixed point.
    pub fn add_fixed_point(&mut self, x: f64, dxdt: f64) {
        self.fixed_points.push([x, dxdt]);
    }
    /// Classify a fixed point as stable / unstable / saddle.
    ///
    /// Uses the trace and determinant of the Jacobian.
    /// `j` = [df/dx, df/dy, dg/dx, dg/dy] (2×2 row-major).
    pub fn classify_fixed_point(j: [f64; 4]) -> &'static str {
        let trace = j[0] + j[3];
        let det = j[0] * j[3] - j[1] * j[2];
        if det < 0.0 {
            "saddle"
        } else if trace < 0.0 {
            "stable"
        } else if trace > 0.0 {
            "unstable"
        } else {
            "centre"
        }
    }
    /// Compute the bounding box of the trajectory.
    pub fn bounding_box(&self) -> ([f64; 2], [f64; 2]) {
        if self.x.is_empty() {
            return ([0.0, 0.0], [1.0, 1.0]);
        }
        (
            [slice_min(&self.x), slice_min(&self.dxdt)],
            [slice_max(&self.x), slice_max(&self.dxdt)],
        )
    }
}
/// Axis configuration for a chart.
#[allow(dead_code)]
pub struct AxisConfig {
    /// Axis label.
    pub label: String,
    /// Use logarithmic scale.
    pub log_scale: bool,
    /// Number of major grid lines.
    pub grid_lines: usize,
    /// Number of tick marks.
    pub ticks: usize,
    /// Manual range (min, max), or `None` to auto-range.
    pub range: Option<(f64, f64)>,
}
impl AxisConfig {
    /// Create an axis configuration with defaults.
    pub fn new(label: &str) -> Self {
        Self {
            label: label.to_string(),
            log_scale: false,
            grid_lines: 5,
            ticks: 5,
            range: None,
        }
    }
    /// Compute the data range for this axis given a slice of values.
    pub fn effective_range(&self, data: &[f64]) -> (f64, f64) {
        if let Some(r) = self.range {
            return r;
        }
        let lo = slice_min(data);
        let hi = slice_max(data);
        if (hi - lo).abs() < 1e-12 {
            (lo - 1.0, lo + 1.0)
        } else {
            (lo, hi)
        }
    }
    /// Generate tick positions within the data range.
    pub fn tick_positions(&self, lo: f64, hi: f64) -> Vec<f64> {
        (0..=self.ticks)
            .map(|i| lo + (hi - lo) * i as f64 / self.ticks as f64)
            .collect()
    }
}
/// A named data series for plotting.
#[allow(dead_code)]
pub struct Series {
    /// Series label.
    pub label: String,
    /// X values.
    pub x: Vec<f64>,
    /// Y values.
    pub y: Vec<f64>,
    /// Line colour (RGBA 0..1).
    pub color: [f32; 4],
}
impl Series {
    /// Create a series.
    pub fn new(label: &str, x: Vec<f64>, y: Vec<f64>, color: [f32; 4]) -> Self {
        Self {
            label: label.to_string(),
            x,
            y,
            color,
        }
    }
    /// Number of data points.
    pub fn len(&self) -> usize {
        self.x.len().min(self.y.len())
    }
    /// Returns `true` if the series has no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
/// Multi-series line chart.
#[allow(dead_code)]
pub struct LinePlot {
    /// Data series.
    pub series: Vec<Series>,
    /// X axis configuration.
    pub x_axis: AxisConfig,
    /// Y axis configuration.
    pub y_axis: AxisConfig,
    /// Chart title.
    pub title: String,
    /// Show legend.
    pub show_legend: bool,
}
impl LinePlot {
    /// Create an empty line plot.
    pub fn new(title: &str) -> Self {
        Self {
            series: Vec::new(),
            x_axis: AxisConfig::new("x"),
            y_axis: AxisConfig::new("y"),
            title: title.to_string(),
            show_legend: true,
        }
    }
    /// Add a series.
    pub fn add_series(&mut self, s: Series) {
        self.series.push(s);
    }
    /// Compute the overall x range across all series.
    pub fn x_range(&self) -> (f64, f64) {
        let all_x: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.x.iter().cloned())
            .collect();
        if all_x.is_empty() {
            return (0.0, 1.0);
        }
        self.x_axis.effective_range(&all_x)
    }
    /// Compute the overall y range across all series.
    pub fn y_range(&self) -> (f64, f64) {
        let all_y: Vec<f64> = self
            .series
            .iter()
            .flat_map(|s| s.y.iter().cloned())
            .collect();
        if all_y.is_empty() {
            return (0.0, 1.0);
        }
        self.y_axis.effective_range(&all_y)
    }
    /// Normalize a data point to `\[0, 1\]` chart space.
    pub fn normalize_point(&self, x: f64, y: f64) -> [f64; 2] {
        let (x0, x1) = self.x_range();
        let (y0, y1) = self.y_range();
        [normalize(x, x0, x1), normalize(y, y0, y1)]
    }
}
/// 1D histogram plot.
#[allow(dead_code)]
pub struct HistogramPlot {
    /// Number of bins.
    pub bins: usize,
    /// Bin edges (length = bins + 1).
    pub bin_edges: Vec<f64>,
    /// Bin counts.
    pub counts: Vec<usize>,
    /// Whether to show density (normalised) vs count.
    pub density: bool,
    /// Whether to overlay a normal distribution.
    pub normal_overlay: bool,
    /// Chart title.
    pub title: String,
}
impl HistogramPlot {
    /// Create a histogram from data.
    pub fn new(title: &str, data: &[f64], bins: usize) -> Self {
        let lo = slice_min(data);
        let hi = slice_max(data);
        let range = (hi - lo).max(1e-12);
        let bin_width = range / bins as f64;
        let bin_edges: Vec<f64> = (0..=bins).map(|i| lo + i as f64 * bin_width).collect();
        let mut counts = vec![0usize; bins];
        for &v in data {
            let b = ((v - lo) / bin_width) as usize;
            let b = b.min(bins - 1);
            counts[b] += 1;
        }
        Self {
            bins,
            bin_edges,
            counts,
            density: false,
            normal_overlay: false,
            title: title.to_string(),
        }
    }
    /// Compute density values (counts / (n * bin_width)).
    pub fn densities(&self) -> Vec<f64> {
        let n: usize = self.counts.iter().sum();
        let bw = if self.bin_edges.len() > 1 {
            self.bin_edges[1] - self.bin_edges[0]
        } else {
            1.0
        };
        self.counts
            .iter()
            .map(|&c| c as f64 / (n.max(1) as f64 * bw))
            .collect()
    }
    /// Compute normal distribution overlay values at bin centres.
    ///
    /// Uses sample mean and variance from the bin counts.
    pub fn normal_pdf_overlay(&self) -> Vec<f64> {
        let densities = self.densities();
        let n = self.bins;
        let centres: Vec<f64> = (0..n)
            .map(|i| (self.bin_edges[i] + self.bin_edges[i + 1]) / 2.0)
            .collect();
        let mean = centres
            .iter()
            .zip(densities.iter())
            .map(|(c, d)| c * d)
            .sum::<f64>();
        let var = centres
            .iter()
            .zip(densities.iter())
            .map(|(c, d)| (c - mean).powi(2) * d)
            .sum::<f64>()
            .max(1e-12);
        let std = var.sqrt();
        use std::f64::consts::PI;
        centres
            .iter()
            .map(|&c| {
                let z = (c - mean) / std;
                (-0.5 * z * z).exp() / (std * (2.0 * PI).sqrt())
            })
            .collect()
    }
}
/// Scatter plot renderer.
#[allow(dead_code)]
pub struct ScatterPlot {
    /// Data points.
    pub points: Vec<ScatterPoint>,
    /// Marker style.
    pub marker: MarkerStyle,
    /// Scalar value range for colour mapping.
    pub value_range: (f64, f64),
    /// X axis configuration.
    pub x_axis: AxisConfig,
    /// Y axis configuration.
    pub y_axis: AxisConfig,
    /// Chart title.
    pub title: String,
}
impl ScatterPlot {
    /// Create a scatter plot.
    pub fn new(title: &str) -> Self {
        Self {
            points: Vec::new(),
            marker: MarkerStyle::Circle,
            value_range: (0.0, 1.0),
            x_axis: AxisConfig::new("x"),
            y_axis: AxisConfig::new("y"),
            title: title.to_string(),
        }
    }
    /// Add a point.
    pub fn add_point(&mut self, x: f64, y: f64, value: f64, size: f32) {
        self.points.push(ScatterPoint { x, y, value, size });
    }
    /// Map a scalar value to a colour (cool-to-warm).
    pub fn value_to_color(&self, value: f64) -> [f32; 4] {
        let t = normalize(value, self.value_range.0, self.value_range.1) as f32;
        [t, 0.2, 1.0 - t, 1.0]
    }
    /// Return the number of points within a rectangular region.
    pub fn count_in_region(&self, x0: f64, x1: f64, y0: f64, y1: f64) -> usize {
        self.points
            .iter()
            .filter(|p| p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
            .count()
    }
}
/// 3D bubble chart (x, y, z = size).
#[allow(dead_code)]
pub struct BubbleChart {
    /// Data points.
    pub points: Vec<BubblePoint>,
    /// Maximum bubble display size (pixels).
    pub max_size: f32,
    /// Minimum bubble display size (pixels).
    pub min_size: f32,
    /// Chart title.
    pub title: String,
}
impl BubbleChart {
    /// Create a bubble chart.
    pub fn new(title: &str) -> Self {
        Self {
            points: Vec::new(),
            max_size: 50.0,
            min_size: 5.0,
            title: title.to_string(),
        }
    }
    /// Add a point.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64, color: [f32; 4]) {
        self.points.push(BubblePoint {
            x,
            y,
            z,
            label: None,
            color,
        });
    }
    /// Compute display radius for a z value.
    pub fn display_size(&self, z: f64) -> f32 {
        let z_values: Vec<f64> = self.points.iter().map(|p| p.z).collect();
        if z_values.is_empty() {
            return self.min_size;
        }
        let lo = slice_min(&z_values);
        let hi = slice_max(&z_values);
        let t = normalize(z, lo, hi) as f32;
        self.min_size + (self.max_size - self.min_size) * t
    }
    /// Return all points sorted by z descending (paint order).
    pub fn sorted_by_z_desc(&self) -> Vec<&BubblePoint> {
        let mut refs: Vec<&BubblePoint> = self.points.iter().collect();
        refs.sort_by(|a, b| b.z.partial_cmp(&a.z).unwrap_or(std::cmp::Ordering::Equal));
        refs
    }
}
/// Parallel coordinates chart for multi-dimensional physics data.
#[allow(dead_code)]
pub struct ParallelCoord {
    /// Axis names.
    pub axes: Vec<String>,
    /// Data rows (one entry per sample, one value per axis).
    pub data: Vec<Vec<f64>>,
    /// Axis ranges (auto-computed if None).
    pub axis_ranges: Vec<Option<(f64, f64)>>,
    /// Line colour per row.
    pub colors: Vec<[f32; 4]>,
    /// Chart title.
    pub title: String,
}
impl ParallelCoord {
    /// Create a parallel coordinates chart.
    pub fn new(title: &str, axes: Vec<String>) -> Self {
        let n = axes.len();
        Self {
            axes,
            data: Vec::new(),
            axis_ranges: vec![None; n],
            colors: Vec::new(),
            title: title.to_string(),
        }
    }
    /// Add a data row.
    pub fn add_row(&mut self, row: Vec<f64>, color: [f32; 4]) {
        self.data.push(row);
        self.colors.push(color);
    }
    /// Compute effective range for axis `i`.
    pub fn axis_range(&self, i: usize) -> (f64, f64) {
        if let Some(Some(r)) = self.axis_ranges.get(i) {
            return *r;
        }
        let values: Vec<f64> = self
            .data
            .iter()
            .filter_map(|row| row.get(i).cloned())
            .collect();
        if values.is_empty() {
            return (0.0, 1.0);
        }
        let lo = slice_min(&values);
        let hi = slice_max(&values);
        if (hi - lo).abs() < 1e-12 {
            (lo - 1.0, lo + 1.0)
        } else {
            (lo, hi)
        }
    }
    /// Normalise row `row_idx` to [0, 1] per axis.
    pub fn normalized_row(&self, row_idx: usize) -> Vec<f64> {
        let row = &self.data[row_idx];
        row.iter()
            .enumerate()
            .map(|(i, &v)| {
                let (lo, hi) = self.axis_range(i);
                normalize(v, lo, hi)
            })
            .collect()
    }
}
/// Multi-line energy time series chart.
#[allow(dead_code)]
pub struct EnergyTimeSeries {
    /// Time values.
    pub time: Vec<f64>,
    /// Kinetic energy.
    pub ke: Vec<f64>,
    /// Potential energy.
    pub pe: Vec<f64>,
    /// Total energy.
    pub total: Vec<f64>,
    /// Chart title.
    pub title: String,
    /// Extra named energy channels (e.g., "elastic").
    pub extra: Vec<(String, Vec<f64>)>,
}
impl EnergyTimeSeries {
    /// Create an energy time series chart.
    pub fn new(title: &str) -> Self {
        Self {
            time: Vec::new(),
            ke: Vec::new(),
            pe: Vec::new(),
            total: Vec::new(),
            title: title.to_string(),
            extra: Vec::new(),
        }
    }
    /// Append an energy snapshot.
    pub fn push(&mut self, t: f64, ke: f64, pe: f64) {
        self.time.push(t);
        self.ke.push(ke);
        self.pe.push(pe);
        self.total.push(ke + pe);
    }
    /// Compute relative energy drift (max(total) - min(total)) / mean(total).
    pub fn energy_drift(&self) -> f64 {
        if self.total.is_empty() {
            return 0.0;
        }
        let lo = slice_min(&self.total);
        let hi = slice_max(&self.total);
        let mean = self.total.iter().sum::<f64>() / self.total.len() as f64;
        if mean.abs() < 1e-12 {
            0.0
        } else {
            (hi - lo) / mean.abs()
        }
    }
    /// Add an extra energy channel.
    pub fn add_channel(&mut self, name: &str, data: Vec<f64>) {
        self.extra.push((name.to_string(), data));
    }
}
/// A single telemetry sample.
#[allow(dead_code)]
pub struct TelemetrySample {
    /// Timestamp (seconds).
    pub time: f64,
    /// Front-left tire longitudinal force (N).
    pub tire_fl: f64,
    /// Front-right tire longitudinal force (N).
    pub tire_fr: f64,
    /// Rear-left tire longitudinal force (N).
    pub tire_rl: f64,
    /// Rear-right tire longitudinal force (N).
    pub tire_rr: f64,
    /// Lateral g-force (m/s²).
    pub g_lateral: f64,
    /// Longitudinal g-force (m/s²).
    pub g_longitudinal: f64,
    /// Speed (m/s).
    pub speed: f64,
}
/// Grouped / stacked bar chart.
#[allow(dead_code)]
pub struct BarChart {
    /// Group labels.
    pub group_labels: Vec<String>,
    /// Data: each entry is a group with one value per category.
    pub groups: Vec<Vec<f64>>,
    /// Category labels.
    pub category_labels: Vec<String>,
    /// Whether to render as stacked (true) or grouped (false).
    pub stacked: bool,
    /// Chart title.
    pub title: String,
    /// Bar colors per group.
    pub colors: Vec<[f32; 4]>,
}
impl BarChart {
    /// Creates an empty bar chart.
    pub fn new(title: &str) -> Self {
        Self {
            group_labels: Vec::new(),
            groups: Vec::new(),
            category_labels: Vec::new(),
            stacked: false,
            title: title.to_string(),
            colors: Vec::new(),
        }
    }
    /// Adds a group with values for each category.
    pub fn add_group(&mut self, label: &str, values: Vec<f64>) {
        self.group_labels.push(label.to_string());
        self.groups.push(values);
        let n = self.colors.len() as f32 / 8.0;
        self.colors
            .push([0.2 + n * 0.1, 0.4 + n * 0.05, 0.8 - n * 0.1, 1.0]);
    }
    /// Returns the number of categories (length of the first group, or 0).
    pub fn n_categories(&self) -> usize {
        self.groups.first().map(|g| g.len()).unwrap_or(0)
    }
    /// Computes the sum of values in group `idx`.
    pub fn group_total(&self, idx: usize) -> f64 {
        self.groups.get(idx).map(|g| g.iter().sum()).unwrap_or(0.0)
    }
    /// Computes the maximum stacked height across all categories.
    pub fn stacked_max(&self) -> f64 {
        let n_cat = self.n_categories();
        (0..n_cat)
            .map(|cat| {
                self.groups
                    .iter()
                    .map(|g| g.get(cat).cloned().unwrap_or(0.0))
                    .sum::<f64>()
            })
            .fold(0.0_f64, f64::max)
    }
    /// Computes the maximum grouped bar height.
    pub fn grouped_max(&self) -> f64 {
        self.groups
            .iter()
            .flat_map(|g| g.iter().cloned())
            .fold(0.0_f64, f64::max)
    }
}
/// Axis tick marks with formatted string labels.
#[allow(dead_code)]
pub struct AxisTicksLabel {
    /// Tick positions in data coordinates.
    pub positions: Vec<f64>,
    /// Formatted labels.
    pub labels: Vec<String>,
}
impl AxisTicksLabel {
    /// Generates `n` linearly-spaced ticks between `lo` and `hi`.
    pub fn linear(lo: f64, hi: f64, n: usize) -> Self {
        let positions: Vec<f64> = (0..n)
            .map(|i| {
                lo + if n > 1 {
                    (hi - lo) * i as f64 / (n - 1) as f64
                } else {
                    0.0
                }
            })
            .collect();
        let labels = positions.iter().map(|&p| format!("{:.3}", p)).collect();
        Self { positions, labels }
    }
    /// Generates ticks at powers of 10 between `lo` and `hi`.
    pub fn logarithmic(lo: f64, hi: f64) -> Self {
        let lo_log = lo.max(1e-15).log10().floor() as i32;
        let hi_log = hi.max(1e-15).log10().ceil() as i32;
        let positions: Vec<f64> = (lo_log..=hi_log)
            .map(|e| 10.0_f64.powi(e))
            .filter(|&v| v >= lo && v <= hi)
            .collect();
        let labels = positions.iter().map(|&p| format!("{:.2e}", p)).collect();
        Self { positions, labels }
    }
}
/// SVG string exporter for chart rendering.
///
/// Produces minimal SVG markup suitable for embedding in HTML or saving to disk.
#[allow(dead_code)]
pub struct SvgExporter {
    /// Canvas width in pixels.
    pub width: usize,
    /// Canvas height in pixels.
    pub height: usize,
}
impl SvgExporter {
    /// Creates a new SVG exporter with the given dimensions.
    pub fn new(width: usize, height: usize) -> Self {
        Self { width, height }
    }
    /// Generates the SVG opening tag with viewBox.
    pub fn export_header(&self) -> String {
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">",
            self.width, self.height, self.width, self.height
        )
    }
    /// Generates an SVG ``line` element.
    pub fn line_element(
        &self,
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        color: &str,
        width: f64,
    ) -> String {
        format!(
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{:.2}\"/>",
            x1, y1, x2, y2, color, width
        )
    }
    /// Generates an SVG ``circle` element.
    pub fn circle_element(&self, cx: f64, cy: f64, r: f64, color: &str) -> String {
        format!(
            "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"{}\"/>",
            cx, cy, r, color
        )
    }
    /// Generates an SVG ``rect` element.
    pub fn rect_element(&self, x: f64, y: f64, w: f64, h: f64, color: &str) -> String {
        format!(
            "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>",
            x, y, w, h, color
        )
    }
    /// Generates an SVG ``text` element.
    pub fn text_element(&self, x: f64, y: f64, text: &str, font_size: usize) -> String {
        format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"{}\">{}</text>",
            x, y, font_size, text
        )
    }
    /// Generates a ``polyline` element from a list of (x, y) pairs.
    pub fn polyline_element(&self, pts: &[(f64, f64)], color: &str, stroke_width: f64) -> String {
        let pts_str: String = pts
            .iter()
            .map(|(x, y)| format!("{:.2},{:.2}", x, y))
            .collect::<Vec<_>>()
            .join(" ");
        format!(
            "<polyline points=\"{}\" stroke=\"{}\" fill=\"none\" stroke-width=\"{:.2}\"/>",
            pts_str, color, stroke_width
        )
    }
}
/// Circular ring buffer for real-time chart updates.
///
/// Maintains the last `capacity` samples in insertion order.
#[allow(dead_code)]
pub struct RealTimeBuffer {
    /// Maximum number of samples to retain.
    pub capacity: usize,
    /// Stored samples (most-recent at the end).
    pub samples: std::collections::VecDeque<f64>,
}
impl RealTimeBuffer {
    /// Creates a new real-time buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            samples: std::collections::VecDeque::with_capacity(capacity),
        }
    }
    /// Pushes a new sample, evicting the oldest if at capacity.
    pub fn push_sample(&mut self, value: f64) {
        if self.samples.len() >= self.capacity {
            self.samples.pop_front();
        }
        self.samples.push_back(value);
    }
    /// Returns the most-recently added sample, or `None` if empty.
    pub fn latest(&self) -> Option<f64> {
        self.samples.back().cloned()
    }
    /// Returns the arithmetic mean of all stored samples.
    pub fn mean(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().sum::<f64>() / self.samples.len() as f64
    }
    /// Returns the minimum stored sample value.
    pub fn min_val(&self) -> f64 {
        self.samples.iter().cloned().fold(f64::MAX, f64::min)
    }
    /// Returns the maximum stored sample value.
    pub fn max_val(&self) -> f64 {
        self.samples.iter().cloned().fold(f64::MIN, f64::max)
    }
    /// Copies current samples into a `Vec`f64` in chronological order.
    pub fn to_vec(&self) -> Vec<f64> {
        self.samples.iter().cloned().collect()
    }
}
/// 2D heatmap with optional contour overlays.
#[allow(dead_code)]
pub struct HeatmapPlot {
    /// Grid width (columns).
    pub width: usize,
    /// Grid height (rows).
    pub height: usize,
    /// Scalar data (row-major).
    pub data: Vec<f64>,
    /// X axis label.
    pub x_label: String,
    /// Y axis label.
    pub y_label: String,
    /// Colormap name (stored for reference).
    pub colormap: String,
    /// Number of contour levels.
    pub contour_levels: usize,
    /// Chart title.
    pub title: String,
}
impl HeatmapPlot {
    /// Create a heatmap.
    pub fn new(title: &str, data: Vec<f64>, width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            data,
            x_label: "x".to_string(),
            y_label: "y".to_string(),
            colormap: "plasma".to_string(),
            contour_levels: 8,
            title: title.to_string(),
        }
    }
    /// Map a scalar value to a colour using a plasma-like palette.
    pub fn map_color(&self, value: f64) -> [f32; 4] {
        let lo = slice_min(&self.data);
        let hi = slice_max(&self.data);
        let t = normalize(value, lo, hi) as f32;
        let r = t;
        let g = 0.2 * (1.0 - t);
        let b = 1.0 - t;
        [r, g, b, 1.0]
    }
    /// Compute contour thresholds.
    pub fn contour_thresholds(&self) -> Vec<f64> {
        let lo = slice_min(&self.data);
        let hi = slice_max(&self.data);
        (0..=self.contour_levels)
            .map(|i| lo + (hi - lo) * i as f64 / self.contour_levels as f64)
            .collect()
    }
    /// Sample value at (row, col) with bilinear interpolation.
    pub fn sample(&self, row: f64, col: f64) -> f64 {
        let r = row.clamp(0.0, (self.height - 1) as f64);
        let c = col.clamp(0.0, (self.width - 1) as f64);
        let ri = r as usize;
        let ci = c as usize;
        let ri1 = (ri + 1).min(self.height - 1);
        let ci1 = (ci + 1).min(self.width - 1);
        let tr = r - ri as f64;
        let tc = c - ci as f64;
        let v00 = self.data[ri * self.width + ci];
        let v10 = self.data[ri1 * self.width + ci];
        let v01 = self.data[ri * self.width + ci1];
        let v11 = self.data[ri1 * self.width + ci1];
        (1.0 - tr) * (1.0 - tc) * v00
            + tr * (1.0 - tc) * v10
            + (1.0 - tr) * tc * v01
            + tr * tc * v11
    }
}
/// Scatter plot marker shapes.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MarkerStyle {
    /// Circular marker.
    Circle,
    /// Cross marker.
    Cross,
    /// Square marker.
    Square,
    /// Triangle marker.
    Triangle,
}
