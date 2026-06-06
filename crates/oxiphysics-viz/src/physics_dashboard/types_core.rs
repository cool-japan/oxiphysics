//! Core types for physics dashboard
//!
//! Foundation types with no dependencies on other dashboard types.

/// Colormap legend data for dashboard rendering.
#[derive(Debug, Clone)]
pub struct ColormapLegend {
    /// Minimum scalar value.
    pub vmin: f64,
    /// Maximum scalar value.
    pub vmax: f64,
    /// Number of legend ticks.
    pub n_ticks: usize,
    /// Label format string.
    pub label_format: String,
    /// Title of the legend.
    pub title: String,
    /// Whether to use log scale.
    pub log_scale: bool,
}
impl ColormapLegend {
    /// Create with linear scale.
    pub fn linear(vmin: f64, vmax: f64, title: &str) -> Self {
        Self {
            vmin,
            vmax,
            n_ticks: 5,
            label_format: "{:.3g}".to_string(),
            title: title.to_string(),
            log_scale: false,
        }
    }
    /// Tick values for the legend.
    pub fn tick_values(&self) -> Vec<f64> {
        if self.log_scale && self.vmin > 0.0 {
            let log_min = self.vmin.log10();
            let log_max = self.vmax.log10();
            (0..self.n_ticks)
                .map(|i| {
                    let t = i as f64 / (self.n_ticks - 1).max(1) as f64;
                    10_f64.powf(log_min + t * (log_max - log_min))
                })
                .collect()
        } else {
            (0..self.n_ticks)
                .map(|i| {
                    let t = i as f64 / (self.n_ticks - 1).max(1) as f64;
                    self.vmin + t * (self.vmax - self.vmin)
                })
                .collect()
        }
    }
    /// Normalize a value to \[0,1\] for colormap lookup.
    pub fn normalize(&self, v: f64) -> f64 {
        if self.log_scale && self.vmin > 0.0 && self.vmax > 0.0 && v > 0.0 {
            let log_v = v.log10();
            let log_min = self.vmin.log10();
            let log_max = self.vmax.log10();
            ((log_v - log_min) / (log_max - log_min).max(1e-30)).clamp(0.0, 1.0)
        } else {
            ((v - self.vmin) / (self.vmax - self.vmin).max(1e-30)).clamp(0.0, 1.0)
        }
    }
    /// Update range to fit data.
    pub fn fit_range(&mut self, data: &[f64]) {
        if data.is_empty() {
            return;
        }
        self.vmin = data.iter().cloned().fold(f64::INFINITY, f64::min);
        self.vmax = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if (self.vmax - self.vmin).abs() < 1e-30 {
            self.vmax = self.vmin + 1.0;
        }
    }
}
/// A 2-D heatmap of scalar field values (pressure, stress, temperature, …).
///
/// Supports colour mapping, normalization per cell, and bilinear downsampling
/// for efficient rendering at reduced resolution.
#[derive(Debug, Clone)]
pub struct HeatmapRenderer {
    /// Grid data: `grid[row][col]`.
    pub grid: Vec<Vec<f64>>,
    /// Physical label for this field.
    pub field_label: String,
}
impl HeatmapRenderer {
    /// Build a heatmap from a 2-D grid of values.
    pub fn from_grid(grid: Vec<Vec<f64>>, field_label: &str) -> Self {
        Self {
            grid,
            field_label: field_label.to_string(),
        }
    }
    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.grid.len()
    }
    /// Number of columns (of the first row, or 0 if empty).
    pub fn cols(&self) -> usize {
        self.grid.first().map(|r| r.len()).unwrap_or(0)
    }
    /// Global minimum value across all cells.
    pub fn global_min(&self) -> f64 {
        self.grid
            .iter()
            .flat_map(|r| r.iter())
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
    /// Global maximum value across all cells.
    pub fn global_max(&self) -> f64 {
        self.grid
            .iter()
            .flat_map(|r| r.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Normalised value at cell (row, col) in \[0, 1\] relative to global min/max.
    pub fn cell_normalized(&self, row: usize, col: usize) -> f64 {
        let mn = self.global_min();
        let mx = self.global_max();
        let span = (mx - mn).max(1e-15);
        (self.grid[row][col] - mn) / span
    }
    /// Nearest-neighbour resample to `out_rows × out_cols`.
    pub fn resample_to(&self, out_rows: usize, out_cols: usize) -> HeatmapRenderer {
        let src_r = self.rows().max(1);
        let src_c = self.cols().max(1);
        let new_grid: Vec<Vec<f64>> = (0..out_rows)
            .map(|r| {
                (0..out_cols)
                    .map(|c| {
                        let sr = (r * src_r / out_rows.max(1)).min(src_r - 1);
                        let sc = (c * src_c / out_cols.max(1)).min(src_c - 1);
                        self.grid[sr][sc]
                    })
                    .collect()
            })
            .collect();
        HeatmapRenderer::from_grid(new_grid, &self.field_label)
    }
    /// Mean value across the entire grid.
    pub fn mean(&self) -> f64 {
        let n: usize = self.grid.iter().map(|r| r.len()).sum();
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = self.grid.iter().flat_map(|r| r.iter()).sum();
        sum / n as f64
    }
}
/// A single run result in a parameter sweep.
#[derive(Debug, Clone)]
pub struct SweepResult {
    /// Parameter values for this run.
    pub parameters: Vec<f64>,
    /// Result metrics for this run.
    pub metrics: Vec<f64>,
}
/// A sample on a material's stress-strain curve.
#[derive(Debug, Clone, Copy)]
pub struct StressStrainPoint {
    /// Strain (dimensionless).
    pub strain: f64,
    /// Stress (Pa).
    pub stress: f64,
}
/// A time series of scalar values from a simulation run.
#[derive(Debug, Clone)]
pub struct SimTrace {
    /// Name of this trace.
    pub name: String,
    /// Time stamps.
    pub times: Vec<f64>,
    /// Values at each time stamp.
    pub values: Vec<f64>,
}
impl SimTrace {
    /// Create a new trace.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            times: Vec::new(),
            values: Vec::new(),
        }
    }
    /// Add a sample.
    pub fn push(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }
    /// Interpolate value at `t` (linear).
    pub fn sample_at(&self, t: f64) -> Option<f64> {
        let n = self.times.len();
        if n == 0 {
            return None;
        }
        let idx = self.times.partition_point(|&ti| ti <= t);
        if idx == 0 {
            return Some(self.values[0]);
        }
        if idx >= n {
            return Some(*self.values.last().expect("collection should not be empty"));
        }
        let t0 = self.times[idx - 1];
        let t1 = self.times[idx];
        let v0 = self.values[idx - 1];
        let v1 = self.values[idx];
        let frac = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
        Some(v0 + (v1 - v0) * frac)
    }
    /// Mean value.
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }
    /// Standard deviation.
    pub fn std_dev(&self) -> f64 {
        if self.values.len() < 2 {
            return 0.0;
        }
        let m = self.mean();
        let var = self.values.iter().map(|&v| (v - m).powi(2)).sum::<f64>()
            / (self.values.len() - 1) as f64;
        var.sqrt()
    }
}
/// Minimal SVG plot builder.
#[derive(Debug, Clone, Default)]
pub struct SvgPlotBuilder {
    /// Plot elements accumulated.
    elements: Vec<String>,
    /// SVG width (px).
    pub width: u32,
    /// SVG height (px).
    pub height: u32,
    /// Plot area margins \[left, right, top, bottom\].
    pub margins: [u32; 4],
}
impl SvgPlotBuilder {
    /// Create a builder.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            elements: Vec::new(),
            width,
            height,
            margins: [60, 20, 20, 40],
        }
    }
    /// Plot width (inside margins).
    pub fn plot_width(&self) -> u32 {
        self.width.saturating_sub(self.margins[0] + self.margins[1])
    }
    /// Plot height (inside margins).
    pub fn plot_height(&self) -> u32 {
        self.height
            .saturating_sub(self.margins[2] + self.margins[3])
    }
    /// Map data coordinate to SVG pixel.
    pub fn map_x(&self, x: f64, x_min: f64, x_max: f64) -> f64 {
        let range = (x_max - x_min).max(1e-30);
        self.margins[0] as f64 + (x - x_min) / range * self.plot_width() as f64
    }
    /// Map a data y-value to pixel y-coordinate within the plot area.
    pub fn map_y(&self, y: f64, y_min: f64, y_max: f64) -> f64 {
        let range = (y_max - y_min).max(1e-30);
        (self.margins[2] + self.plot_height()) as f64
            - (y - y_min) / range * self.plot_height() as f64
    }
    /// Add a polyline from data series.
    pub fn add_line(&mut self, xs: &[f64], ys: &[f64], color: &str, stroke_width: f64) {
        let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let points: String = xs
            .iter()
            .zip(ys.iter())
            .map(|(&x, &y)| {
                format!(
                    "{:.1},{:.1}",
                    self.map_x(x, x_min, x_max),
                    self.map_y(y, y_min, y_max)
                )
            })
            .collect::<Vec<_>>()
            .join(" ");
        self.elements
            .push(
                format!(
                    r#"<polyline points="{points}" fill="none" stroke="{color}" stroke-width="{stroke_width:.1}"/>"#
                ),
            );
    }
    /// Add a text label.
    pub fn add_text(&mut self, x: f64, y: f64, text: &str, color: &str, font_size: f64) {
        self.elements.push(format!(
            r#"<text x="{x:.1}" y="{y:.1}" fill="{color}" font-size="{font_size:.1}">{text}</text>"#
        ));
    }
    /// Add a background rectangle.
    pub fn add_background(&mut self, color: &str) {
        self.elements.push(format!(
            r#"<rect width="{}" height="{}" fill="{color}"/>"#,
            self.width, self.height
        ));
    }
    /// Render to SVG string.
    pub fn render(&self) -> String {
        let body = self.elements.join("\n  ");
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">{}</svg>"#,
            self.width, self.height, body
        )
    }
    /// Clear all elements.
    pub fn clear(&mut self) {
        self.elements.clear();
    }
}
/// A single snapshot of simulation state for recording.
#[derive(Debug, Clone)]
pub struct SimSnapshot {
    /// Simulation time (s).
    pub time: f64,
    /// Kinetic energy (J).
    pub ke: f64,
    /// Potential energy (J).
    pub pe: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Number of active particles/bodies.
    pub n_active: usize,
    /// Maximum force magnitude (N).
    pub max_force: f64,
    /// RMS velocity (m/s).
    pub rms_velocity: f64,
}
impl SimSnapshot {
    /// Total mechanical energy.
    pub fn total_energy(&self) -> f64 {
        self.ke + self.pe
    }
}
/// Snapshot of a rigid body's state.
#[derive(Debug, Clone)]
pub struct BodyState {
    /// Body index.
    pub index: usize,
    /// Mass (kg).
    pub mass: f64,
    /// Inertia tensor diagonal (kg·m²).
    pub inertia: [f64; 3],
    /// Linear velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
    /// Net force (N).
    pub net_force: [f64; 3],
    /// Net torque (N·m).
    pub net_torque: [f64; 3],
    /// Whether the body is sleeping.
    pub is_sleeping: bool,
}
impl BodyState {
    /// Kinetic energy of this body (translational + rotational).
    pub fn kinetic_energy(&self) -> f64 {
        let v = self.velocity;
        let w = self.angular_velocity;
        let ke_trans = 0.5 * self.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        let ke_rot = 0.5
            * (self.inertia[0] * w[0] * w[0]
                + self.inertia[1] * w[1] * w[1]
                + self.inertia[2] * w[2] * w[2]);
        ke_trans + ke_rot
    }
}
/// A single frame's timing breakdown (milliseconds).
#[derive(Debug, Clone, Copy, Default)]
pub struct FrameTiming {
    /// Broad-phase collision detection time (ms).
    pub broadphase_ms: f64,
    /// Narrow-phase collision detection time (ms).
    pub narrowphase_ms: f64,
    /// Constraint solver time (ms).
    pub solver_ms: f64,
    /// Integration time (ms).
    pub integration_ms: f64,
    /// Memory used (bytes).
    pub memory_bytes: usize,
}
impl FrameTiming {
    /// Total physics step time.
    pub fn total_ms(&self) -> f64 {
        self.broadphase_ms + self.narrowphase_ms + self.solver_ms + self.integration_ms
    }
}
/// Interactive chart viewport: zoom in/out and pan left/right.
///
/// Maintains an axis-aligned view rectangle and applies scale/translate
/// operations while clamping to a user-defined data domain.
#[derive(Debug, Clone)]
pub struct ChartZoomPan {
    /// Current view x-minimum.
    pub xmin: f64,
    /// Current view x-maximum.
    pub xmax: f64,
    /// Current view y-minimum.
    pub ymin: f64,
    /// Current view y-maximum.
    pub ymax: f64,
    /// Initial data domain x-minimum (zoom-out limit).
    domain_xmin: f64,
    /// Initial data domain x-maximum.
    domain_xmax: f64,
    /// Initial data domain y-minimum.
    domain_ymin: f64,
    /// Initial data domain y-maximum.
    domain_ymax: f64,
}
impl ChartZoomPan {
    /// Create a new chart view covering the full data domain.
    pub fn new(xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> Self {
        Self {
            xmin,
            xmax,
            ymin,
            ymax,
            domain_xmin: xmin,
            domain_xmax: xmax,
            domain_ymin: ymin,
            domain_ymax: ymax,
        }
    }
    /// Current view as (xmin, xmax, ymin, ymax).
    pub fn view(&self) -> (f64, f64, f64, f64) {
        (self.xmin, self.xmax, self.ymin, self.ymax)
    }
    /// Zoom by `factor` centred on the current view centre.
    ///
    /// `factor` < 1.0 zooms in; `factor` > 1.0 zooms out (clamped to domain).
    pub fn zoom(&mut self, factor: f64) {
        let cx = (self.xmin + self.xmax) * 0.5;
        let cy = (self.ymin + self.ymax) * 0.5;
        let hw = (self.xmax - self.xmin) * 0.5 * factor;
        let hh = (self.ymax - self.ymin) * 0.5 * factor;
        self.xmin = (cx - hw).max(self.domain_xmin);
        self.xmax = (cx + hw).min(self.domain_xmax);
        self.ymin = (cy - hh).max(self.domain_ymin);
        self.ymax = (cy + hh).min(self.domain_ymax);
    }
    /// Pan the view by (dx, dy) in data coordinates.
    pub fn pan(&mut self, dx: f64, dy: f64) {
        let span_x = self.xmax - self.xmin;
        let span_y = self.ymax - self.ymin;
        self.xmin = (self.xmin + dx).max(self.domain_xmin);
        self.xmax = self.xmin + span_x;
        if self.xmax > self.domain_xmax {
            self.xmax = self.domain_xmax;
            self.xmin = self.domain_xmax - span_x;
        }
        self.ymin = (self.ymin + dy).max(self.domain_ymin);
        self.ymax = self.ymin + span_y;
        if self.ymax > self.domain_ymax {
            self.ymax = self.domain_ymax;
            self.ymin = self.domain_ymax - span_y;
        }
    }
    /// Reset view to the full domain.
    pub fn reset(&mut self) {
        self.xmin = self.domain_xmin;
        self.xmax = self.domain_xmax;
        self.ymin = self.domain_ymin;
        self.ymax = self.domain_ymax;
    }
    /// Map a data x-value to a \[0, 1\] normalized screen coordinate.
    pub fn normalize_x(&self, x: f64) -> f64 {
        let span = (self.xmax - self.xmin).max(1e-15);
        ((x - self.xmin) / span).clamp(0.0, 1.0)
    }
    /// Map a data y-value to a \[0, 1\] normalized screen coordinate.
    pub fn normalize_y(&self, y: f64) -> f64 {
        let span = (self.ymax - self.ymin).max(1e-15);
        ((y - self.ymin) / span).clamp(0.0, 1.0)
    }
}
/// Overall simulation status panel.
#[derive(Debug, Clone, Default)]
pub struct SimDashboard {
    /// Current simulation time (seconds).
    pub sim_time: f64,
    /// Current frames per second.
    pub fps: f64,
    /// Number of active rigid bodies.
    pub body_count: usize,
    /// Total kinetic energy (J).
    pub total_ke: f64,
    /// Total potential energy (J).
    pub total_pe: f64,
    /// Whether the simulation is paused.
    pub paused: bool,
    /// Step counter.
    pub step_count: u64,
}
impl SimDashboard {
    /// Create a new dashboard.
    pub fn new() -> Self {
        Self::default()
    }
    /// Total mechanical energy.
    pub fn total_energy(&self) -> f64 {
        self.total_ke + self.total_pe
    }
    /// Advance the simulation time by `dt` and increment step counter.
    pub fn step(&mut self, dt: f64) {
        if !self.paused {
            self.sim_time += dt;
            self.step_count += 1;
        }
    }
    /// Return a formatted status string.
    pub fn status_string(&self) -> String {
        format!(
            "t={:.3}s fps={:.1} bodies={} KE={:.3}J PE={:.3}J E={:.3}J",
            self.sim_time,
            self.fps,
            self.body_count,
            self.total_ke,
            self.total_pe,
            self.total_energy()
        )
    }
}
/// Severity of a physics alert.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum AlertSeverity {
    /// Informational message.
    Info,
    /// Warning: simulation may be inaccurate.
    Warning,
    /// Error: simulation is likely unstable.
    Error,
}
/// Type of a joint constraint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Fixed joint (0 DOF).
    Fixed,
    /// Revolute/hinge joint (1 rotational DOF).
    Revolute,
    /// Prismatic/slider joint (1 linear DOF).
    Prismatic,
    /// Ball-and-socket joint (3 rotational DOF).
    Ball,
    /// Generic spring.
    Spring,
}
/// A single parameter sweep result entry.
#[derive(Debug, Clone)]
pub struct SweepEntry {
    /// Parameter value.
    pub param_value: f64,
    /// Scalar result.
    pub result: f64,
    /// Optional error bar (std dev or range).
    pub error: Option<f64>,
    /// Whether this run converged.
    pub converged: bool,
}
/// Configuration entry for a single named dashboard panel.
#[derive(Debug, Clone)]
pub struct PanelConfig {
    /// Unique name identifier for the panel.
    pub name: String,
    /// Whether the panel is currently visible.
    pub visible: bool,
    /// Preferred width in pixels (0 = auto).
    pub width_px: u32,
    /// Preferred height in pixels (0 = auto).
    pub height_px: u32,
}
impl PanelConfig {
    /// Create a new visible panel configuration.
    pub fn new(name: &str, width_px: u32, height_px: u32) -> Self {
        Self {
            name: name.to_string(),
            visible: true,
            width_px,
            height_px,
        }
    }
}
/// A rolling time-series buffer for plotting physics quantities vs time.
///
/// Maintains at most `capacity` (time, value) samples, evicting the oldest
/// when full. Supports range queries and linear interpolation.
#[derive(Debug, Clone)]
pub struct TimeSeriesBuffer {
    /// Sorted (time, value) pairs.
    pub(crate) data: std::collections::VecDeque<(f64, f64)>,
    /// Maximum number of samples retained.
    capacity: usize,
}
impl TimeSeriesBuffer {
    /// Create a buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            data: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
        }
    }
    /// Append a new (time, value) sample, evicting the oldest if full.
    pub fn push(&mut self, time: f64, value: f64) {
        if self.data.len() >= self.capacity {
            self.data.pop_front();
        }
        self.data.push_back((time, value));
    }
    /// Number of samples currently stored.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    /// All stored time values.
    pub fn times(&self) -> Vec<f64> {
        self.data.iter().map(|&(t, _)| t).collect()
    }
    /// All stored values.
    pub fn values(&self) -> Vec<f64> {
        self.data.iter().map(|&(_, v)| v).collect()
    }
    /// Return (time, value) pairs in the closed time range \[t_lo, t_hi\].
    pub fn range(&self, t_lo: f64, t_hi: f64) -> Vec<(f64, f64)> {
        self.data
            .iter()
            .filter(|&&(t, _)| t >= t_lo && t <= t_hi)
            .copied()
            .collect()
    }
    /// Linearly interpolate value at the given time.
    ///
    /// Returns the nearest boundary value if `t` is outside the data range.
    pub fn interpolate_at(&self, t: f64) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let mut prev: Option<(f64, f64)> = None;
        for &(ti, vi) in &self.data {
            if ti >= t {
                match prev {
                    None => return vi,
                    Some((tp, vp)) => {
                        let dt = ti - tp;
                        if dt < 1e-15 {
                            return vi;
                        }
                        return vp + (vi - vp) * (t - tp) / dt;
                    }
                }
            }
            prev = Some((ti, vi));
        }
        self.data.back().map(|&(_, v)| v).unwrap_or(0.0)
    }
    /// Minimum value in the buffer.
    pub fn min_value(&self) -> f64 {
        self.data
            .iter()
            .map(|&(_, v)| v)
            .fold(f64::INFINITY, f64::min)
    }
    /// Maximum value in the buffer.
    pub fn max_value(&self) -> f64 {
        self.data
            .iter()
            .map(|&(_, v)| v)
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Mean value across all samples.
    pub fn mean_value(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().map(|&(_, v)| v).sum::<f64>() / self.data.len() as f64
    }
    /// Reference to the internal deque of (time, value) pairs.
    pub fn data(&self) -> &std::collections::VecDeque<(f64, f64)> {
        &self.data
    }
    /// First (oldest) sample, if any.
    pub fn front(&self) -> Option<&(f64, f64)> {
        self.data.front()
    }
    /// Last (newest) sample, if any.
    pub fn back(&self) -> Option<&(f64, f64)> {
        self.data.back()
    }
}
/// A time window for chart display.
#[derive(Debug, Clone, Copy)]
pub struct TimeRange {
    /// Start time (seconds).
    pub start: f64,
    /// End time (seconds).
    pub end: f64,
}
impl TimeRange {
    /// Create a new time range.
    pub fn new(start: f64, end: f64) -> Self {
        Self { start, end }
    }
    /// Duration of the range.
    pub fn duration(&self) -> f64 {
        (self.end - self.start).max(0.0)
    }
    /// Map a time value to a normalized `[0, 1]` position.
    pub fn normalize(&self, t: f64) -> f64 {
        let d = self.duration();
        if d < 1e-15 {
            return 0.5;
        }
        ((t - self.start) / d).clamp(0.0, 1.0)
    }
}
/// A multi-channel data streaming buffer for real-time physics metrics.
///
/// Each sample carries a scalar value and a step counter. The buffer is
/// bounded: old samples are evicted when capacity is reached.
#[derive(Debug, Clone)]
pub struct DataStreamBuffer {
    /// (value, step) ring.
    data: std::collections::VecDeque<(f64, u64)>,
    /// Maximum retained samples.
    capacity: usize,
    /// Number of channels (metadata only; all values stored flat).
    pub channels: usize,
}
impl DataStreamBuffer {
    /// Create a streaming buffer with the given capacity and channel count.
    pub fn new(capacity: usize, channels: usize) -> Self {
        Self {
            data: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
            channels: channels.max(1),
        }
    }
    /// Push a (value, step) sample.
    pub fn push_sample(&mut self, value: f64, step: u64) {
        if self.data.len() >= self.capacity {
            self.data.pop_front();
        }
        self.data.push_back((value, step));
    }
    /// Current number of stored samples.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Returns `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    /// The N most recent samples (value, step).
    pub fn latest_n(&self, n: usize) -> Vec<(f64, u64)> {
        let skip = self.data.len().saturating_sub(n);
        self.data.iter().skip(skip).copied().collect()
    }
    /// Mean value across all stored samples.
    pub fn mean_value(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().map(|&(v, _)| v).sum::<f64>() / self.data.len() as f64
    }
    /// Standard deviation of values.
    pub fn std_value(&self) -> f64 {
        if self.data.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_value();
        let var = self
            .data
            .iter()
            .map(|&(v, _)| (v - mean).powi(2))
            .sum::<f64>()
            / (self.data.len() - 1) as f64;
        var.sqrt()
    }
    /// Minimum value stored.
    pub fn min_value(&self) -> f64 {
        self.data
            .iter()
            .map(|&(v, _)| v)
            .fold(f64::INFINITY, f64::min)
    }
    /// Maximum value stored.
    pub fn max_value(&self) -> f64 {
        self.data
            .iter()
            .map(|&(v, _)| v)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
/// A single panel cell in a multi-panel grid layout.
#[derive(Debug, Clone)]
pub struct PanelCell {
    /// Grid row index (0-based).
    pub row: usize,
    /// Grid column index (0-based).
    pub col: usize,
    /// Panel pixel x-origin.
    pub x_px: f64,
    /// Panel pixel y-origin.
    pub y_px: f64,
    /// Panel pixel width.
    pub w_px: f64,
    /// Panel pixel height.
    pub h_px: f64,
    /// Optional name/label for this panel.
    pub name: String,
}
impl PanelCell {
    /// Return the (width, height) in pixels.
    pub fn cell_size(&self) -> (f64, f64) {
        (self.w_px, self.h_px)
    }
    /// Returns `true` if the given screen coordinate falls inside this panel.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.x_px && x < self.x_px + self.w_px && y >= self.y_px && y < self.y_px + self.h_px
    }
}
/// Convergence criterion type.
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceCriterion {
    /// Absolute tolerance on residual.
    Absolute(f64),
    /// Relative tolerance (residual / initial residual).
    Relative(f64),
    /// Both absolute and relative.
    Both(f64, f64),
}
impl ConvergenceCriterion {
    /// Check if residual satisfies criterion.
    pub fn satisfied(&self, residual: f64, initial_residual: f64) -> bool {
        match self {
            ConvergenceCriterion::Absolute(tol) => residual <= *tol,
            ConvergenceCriterion::Relative(tol) => {
                if initial_residual.abs() < 1e-30 {
                    residual <= 1e-12
                } else {
                    residual / initial_residual <= *tol
                }
            }
            ConvergenceCriterion::Both(abs, rel) => {
                residual <= *abs
                    || (initial_residual.abs() > 1e-30 && residual / initial_residual <= *rel)
            }
        }
    }
}
/// Real-time physics metrics: FPS, step time, body count, memory.
///
/// Maintains a rolling window of per-frame measurements and provides
/// instant statistics for HUD display.
#[derive(Debug, Clone)]
pub struct RealTimeMetrics {
    /// Ring buffer of frame step times (ms).
    step_times: Vec<f64>,
    /// Write cursor into the ring buffer.
    cursor: usize,
    /// Number of valid entries (≤ capacity).
    count: usize,
    /// Ring buffer capacity.
    capacity: usize,
}
impl RealTimeMetrics {
    /// Create a real-time metrics tracker with the given rolling window size.
    pub fn new(capacity: usize) -> Self {
        Self {
            step_times: vec![0.0; capacity.max(1)],
            cursor: 0,
            count: 0,
            capacity: capacity.max(1),
        }
    }
    /// Record a new frame step time in milliseconds.
    pub fn record_frame(&mut self, step_ms: f64) {
        self.step_times[self.cursor] = step_ms;
        self.cursor = (self.cursor + 1) % self.capacity;
        if self.count < self.capacity {
            self.count += 1;
        }
    }
    /// Number of recorded frames in the window.
    pub fn len(&self) -> usize {
        self.count
    }
    /// Returns `true` if no frames have been recorded.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
    fn valid_samples(&self) -> &[f64] {
        &self.step_times[..self.count]
    }
    /// Mean step time (ms).
    pub fn mean_step_ms(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        self.valid_samples().iter().sum::<f64>() / self.count as f64
    }
    /// Maximum step time (ms) in the current window.
    pub fn max_step_ms(&self) -> f64 {
        self.valid_samples()
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum step time (ms) in the current window.
    pub fn min_step_ms(&self) -> f64 {
        self.valid_samples()
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
    /// Estimated frames per second (1000 / mean_step_ms).
    pub fn fps(&self) -> f64 {
        let ms = self.mean_step_ms();
        if ms < 1e-10 {
            f64::INFINITY
        } else {
            1000.0 / ms
        }
    }
    /// Standard deviation of step times (ms).
    pub fn std_step_ms(&self) -> f64 {
        if self.count < 2 {
            return 0.0;
        }
        let mean = self.mean_step_ms();
        let var = self
            .valid_samples()
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / (self.count - 1) as f64;
        var.sqrt()
    }
    /// Format a one-line HUD status string.
    pub fn hud_string(&self) -> String {
        format!(
            "FPS:{:.1}  step:{:.2}ms±{:.2}ms  max:{:.2}ms",
            self.fps(),
            self.mean_step_ms(),
            self.std_step_ms(),
            self.max_step_ms()
        )
    }
}
/// Detailed timing breakdown for one simulation step (microseconds).
#[derive(Debug, Clone, Copy, Default)]
pub struct StepBreakdown {
    /// Broad-phase collision detection (µs).
    pub broadphase_us: u64,
    /// Narrow-phase collision detection (µs).
    pub narrowphase_us: u64,
    /// Constraint solver (µs).
    pub solver_us: u64,
    /// Integrator (µs).
    pub integration_us: u64,
    /// Rendering / visualization (µs).
    pub render_us: u64,
}
impl StepBreakdown {
    /// Total step time (µs).
    pub fn total_us(&self) -> u64 {
        self.broadphase_us
            + self.narrowphase_us
            + self.solver_us
            + self.integration_us
            + self.render_us
    }
    /// Total step time in milliseconds.
    pub fn total_ms(&self) -> f64 {
        self.total_us() as f64 / 1000.0
    }
    /// Fraction of time in each phase: \[bp, np, sv, it, rn\].
    pub fn fractions(&self) -> [f64; 5] {
        let total = self.total_us() as f64;
        if total < 1.0 {
            return [0.2; 5];
        }
        [
            self.broadphase_us as f64 / total,
            self.narrowphase_us as f64 / total,
            self.solver_us as f64 / total,
            self.integration_us as f64 / total,
            self.render_us as f64 / total,
        ]
    }
}
/// Histogram renderer for physics scalar distributions (forces, velocities, …).
///
/// Bins data into equal-width intervals and supports normalized PDF rendering.
#[derive(Debug, Clone)]
pub struct HistogramRenderer {
    /// Bin counts.
    pub bins: Vec<u32>,
    /// Left edge of first bin.
    pub range_min: f64,
    /// Right edge of last bin.
    pub range_max: f64,
    /// Label for the x-axis.
    pub label: String,
}
impl HistogramRenderer {
    /// Build a histogram from raw data with `num_bins` equal-width bins.
    pub fn from_data(data: &[f64], num_bins: usize, range_min: f64, range_max: f64) -> Self {
        let nb = num_bins.max(1);
        let mut bins = vec![0u32; nb];
        let span = (range_max - range_min).max(1e-15);
        for &x in data {
            let idx = ((x - range_min) / span * nb as f64) as usize;
            let idx = idx.min(nb - 1);
            bins[idx] += 1;
        }
        Self {
            bins,
            range_min,
            range_max,
            label: String::new(),
        }
    }
    /// Create an empty histogram.
    pub fn new(num_bins: usize, range_min: f64, range_max: f64, label: &str) -> Self {
        Self {
            bins: vec![0u32; num_bins.max(1)],
            range_min,
            range_max,
            label: label.to_string(),
        }
    }
    /// Bin width.
    pub fn bin_width(&self) -> f64 {
        (self.range_max - self.range_min).max(1e-15) / self.bins.len() as f64
    }
    /// Centre of bin `i`.
    pub fn bin_center(&self, i: usize) -> f64 {
        self.range_min + (i as f64 + 0.5) * self.bin_width()
    }
    /// Total sample count.
    pub fn total(&self) -> u32 {
        self.bins.iter().sum()
    }
    /// Normalized histogram (probability density).
    ///
    /// Each bar height = count / (total * bin_width), so the sum integrates to 1.
    pub fn normalized(&self) -> Vec<f64> {
        let total = self.total() as f64;
        let bw = self.bin_width();
        if total < 1e-15 || bw < 1e-15 {
            return vec![0.0; self.bins.len()];
        }
        self.bins.iter().map(|&c| c as f64 / (total * bw)).collect()
    }
    /// Peak bin index.
    pub fn peak_bin(&self) -> usize {
        self.bins
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Add a new value to the histogram.
    pub fn add(&mut self, value: f64) {
        let nb = self.bins.len();
        let span = (self.range_max - self.range_min).max(1e-15);
        let idx = ((value - self.range_min) / span * nb as f64) as usize;
        let idx = idx.min(nb - 1);
        self.bins[idx] += 1;
    }
    /// Reset all bin counts to zero.
    pub fn clear(&mut self) {
        for b in &mut self.bins {
            *b = 0;
        }
    }
}
/// A recorded energy sample.
#[derive(Debug, Clone, Copy)]
pub struct EnergySample {
    /// Simulation time (s).
    pub time: f64,
    /// Kinetic energy (J).
    pub ke: f64,
    /// Potential energy (J).
    pub pe: f64,
}
impl EnergySample {
    /// Total energy.
    pub fn total(&self) -> f64 {
        self.ke + self.pe
    }
}
/// 2-D phase-space scatter plot (position vs velocity, or q vs p).
///
/// Stores up to `capacity` (x, y) points, evicting the oldest to maintain
/// a fixed-size display window for real-time visualization.
#[derive(Debug, Clone)]
pub struct PhaseSpaceDiagram {
    /// Rolling buffer of (x, y) phase-space points.
    points: std::collections::VecDeque<(f64, f64)>,
    /// Maximum points retained.
    capacity: usize,
}
impl PhaseSpaceDiagram {
    /// Create a new phase-space diagram.
    pub fn new(capacity: usize) -> Self {
        Self {
            points: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
        }
    }
    /// Add a (position, velocity) or (q, p) point to the diagram.
    pub fn add_point(&mut self, x: f64, y: f64) {
        if self.points.len() >= self.capacity {
            self.points.pop_front();
        }
        self.points.push_back((x, y));
    }
    /// Number of points in the diagram.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Returns `true` if the diagram has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// All stored points as a slice.
    pub fn points(&self) -> Vec<(f64, f64)> {
        self.points.iter().copied().collect()
    }
    /// Bounding box (xmin, xmax, ymin, ymax) of all points.
    ///
    /// Returns (0,0,0,0) if the diagram is empty.
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        if self.points.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let xmin = self
            .points
            .iter()
            .map(|&(x, _)| x)
            .fold(f64::INFINITY, f64::min);
        let xmax = self
            .points
            .iter()
            .map(|&(x, _)| x)
            .fold(f64::NEG_INFINITY, f64::max);
        let ymin = self
            .points
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::INFINITY, f64::min);
        let ymax = self
            .points
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::NEG_INFINITY, f64::max);
        (xmin, xmax, ymin, ymax)
    }
    /// Centroid (mean x, mean y) of all points.
    pub fn centroid(&self) -> (f64, f64) {
        if self.points.is_empty() {
            return (0.0, 0.0);
        }
        let n = self.points.len() as f64;
        let mx = self.points.iter().map(|&(x, _)| x).sum::<f64>() / n;
        let my = self.points.iter().map(|&(_, y)| y).sum::<f64>() / n;
        (mx, my)
    }
}
/// Panel position in a grid layout.
#[derive(Debug, Clone, Copy)]
pub struct PanelPosition {
    /// Row index (0-based).
    pub row: usize,
    /// Column index (0-based).
    pub col: usize,
    /// Row span.
    pub row_span: usize,
    /// Column span.
    pub col_span: usize,
}
impl PanelPosition {
    /// Single-cell position.
    pub fn single(row: usize, col: usize) -> Self {
        Self {
            row,
            col,
            row_span: 1,
            col_span: 1,
        }
    }
    /// Wide panel spanning multiple columns.
    pub fn wide(row: usize, col: usize, col_span: usize) -> Self {
        Self {
            row,
            col,
            row_span: 1,
            col_span,
        }
    }
}
/// A contact count sample.
#[derive(Debug, Clone, Copy)]
pub struct ContactSample {
    /// Simulation time.
    pub time: f64,
    /// Number of active contacts.
    pub count: usize,
    /// Mean contact force magnitude (N).
    pub mean_force: f64,
    /// Maximum contact force magnitude (N).
    pub max_force: f64,
}
/// Histogram of particle velocity magnitudes.
#[derive(Debug, Clone)]
pub struct VelocityHistogram {
    /// Bin counts.
    pub counts: Vec<u64>,
    /// Bin edges (n_bins + 1 values).
    pub edges: Vec<f64>,
    /// Total sample count.
    pub total: u64,
}
impl VelocityHistogram {
    /// Create with given number of bins and range.
    pub fn new(n_bins: usize, v_min: f64, v_max: f64) -> Self {
        let edges: Vec<f64> = (0..=n_bins)
            .map(|i| v_min + (v_max - v_min) * i as f64 / n_bins as f64)
            .collect();
        Self {
            counts: vec![0; n_bins],
            edges,
            total: 0,
        }
    }
    /// Add a velocity value to the histogram.
    pub fn add(&mut self, v: f64) {
        let n = self.counts.len();
        let v_min = *self.edges.first().unwrap_or(&0.0);
        let v_max = *self.edges.last().unwrap_or(&1.0);
        if v < v_min || v > v_max {
            return;
        }
        let range = (v_max - v_min).max(1e-30);
        let idx = ((v - v_min) / range * n as f64) as usize;
        let idx = idx.min(n - 1);
        self.counts[idx] += 1;
        self.total += 1;
    }
    /// Add a batch of velocity values.
    pub fn add_batch(&mut self, velocities: &[f64]) {
        for &v in velocities {
            self.add(v);
        }
    }
    /// Normalized probability density for each bin.
    pub fn pdf(&self) -> Vec<f64> {
        let n = self.counts.len();
        if self.total == 0 || n == 0 {
            return vec![0.0; n];
        }
        self.counts
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                let dv = self.edges[i + 1] - self.edges[i];
                if dv > 1e-30 {
                    c as f64 / (self.total as f64 * dv)
                } else {
                    0.0
                }
            })
            .collect()
    }
    /// Maxwell-Boltzmann PDF at speed v for temperature T (3D).
    pub fn maxwell_boltzmann(v: f64, kbt: f64, mass: f64) -> f64 {
        if kbt <= 0.0 || mass <= 0.0 || v < 0.0 {
            return 0.0;
        }
        let a = (mass / (2.0 * kbt)).sqrt();
        4.0 * std::f64::consts::PI * a * a * a / std::f64::consts::PI.sqrt()
            * v
            * v
            * (-0.5 * mass * v * v / kbt).exp()
    }
    /// KL divergence from Maxwell-Boltzmann (thermalization check).
    pub fn kl_from_maxwell_boltzmann(&self, kbt: f64, mass: f64) -> f64 {
        let pdf = self.pdf();
        let mut kl = 0.0;
        for (i, &p_data) in pdf.iter().enumerate() {
            let v_mid = 0.5 * (self.edges[i] + self.edges[i + 1]);
            let dv = self.edges[i + 1] - self.edges[i];
            let p_mb = Self::maxwell_boltzmann(v_mid, kbt, mass) * dv;
            if p_data > 1e-30 && p_mb > 1e-30 {
                kl += p_data * dv * (p_data / p_mb).ln();
            }
        }
        kl
    }
    /// Mean velocity.
    pub fn mean(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let n = self.counts.len();
        (0..n)
            .map(|i| {
                let v_mid = 0.5 * (self.edges[i] + self.edges[i + 1]);
                v_mid * self.counts[i] as f64
            })
            .sum::<f64>()
            / self.total as f64
    }
    /// Standard deviation of velocity distribution.
    pub fn std_dev(&self) -> f64 {
        let mean = self.mean();
        if self.total == 0 {
            return 0.0;
        }
        let n = self.counts.len();
        let var = (0..n)
            .map(|i| {
                let v_mid = 0.5 * (self.edges[i] + self.edges[i + 1]);
                let d = v_mid - mean;
                d * d * self.counts[i] as f64
            })
            .sum::<f64>()
            / self.total as f64;
        var.sqrt()
    }
}
