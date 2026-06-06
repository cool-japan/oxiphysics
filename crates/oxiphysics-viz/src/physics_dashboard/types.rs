//! Composite types for physics dashboard
//!
//! Types that depend on core dashboard types.

use super::types_core::*;

/// Body inspector panel.
#[derive(Debug, Clone, Default)]
pub struct BodyInspector {
    /// Currently selected body state.
    pub selected: Option<BodyState>,
}
impl BodyInspector {
    /// Create a new body inspector.
    pub fn new() -> Self {
        Self::default()
    }
    /// Select a body for inspection.
    pub fn select(&mut self, state: BodyState) {
        self.selected = Some(state);
    }
    /// Deselect the body.
    pub fn deselect(&mut self) {
        self.selected = None;
    }
    /// Kinetic energy of the selected body.
    pub fn selected_ke(&self) -> f64 {
        self.selected
            .as_ref()
            .map(|s| s.kinetic_energy())
            .unwrap_or(0.0)
    }
}
/// Linear momentum \[x,y,z\] vs time recorder for conservation verification.
///
/// Measures drift in total system momentum to detect numerical errors or
/// improperly applied external forces.
#[derive(Debug, Clone)]
pub struct MomentumTimeSeries {
    /// Time series buffer for |p| (magnitude of momentum).
    pub mag_series: TimeSeriesBuffer,
    /// Stored \[px,py,pz\] vectors.
    pub vectors: std::collections::VecDeque<[f64; 3]>,
    /// Maximum samples retained.
    pub(super) capacity: usize,
}
impl MomentumTimeSeries {
    /// Create a momentum time series with the given sample capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            mag_series: TimeSeriesBuffer::new(capacity),
            vectors: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
        }
    }
    /// Record total system momentum at simulation time `t`.
    pub fn record(&mut self, t: f64, p: [f64; 3]) {
        let mag = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        self.mag_series.push(t, mag);
        if self.vectors.len() >= self.capacity {
            self.vectors.pop_front();
        }
        self.vectors.push_back(p);
    }
    /// Maximum change in momentum magnitude relative to the first recorded value.
    pub fn drift(&self) -> f64 {
        if self.mag_series.is_empty() {
            return 0.0;
        }
        let first = self.mag_series.front().map(|&(_, v)| v).unwrap_or(0.0);
        let last = self.mag_series.back().map(|&(_, v)| v).unwrap_or(0.0);
        if first.abs() < 1e-15 {
            return (last - first).abs();
        }
        ((last - first) / first).abs()
    }
    /// Returns `true` if momentum is conserved within tolerance.
    pub fn is_conserved(&self, tol: f64) -> bool {
        self.drift() <= tol
    }
}
/// Energy vs time recorder for energy conservation plots.
///
/// Stores kinetic, potential, and total mechanical energy at each simulation step.
#[derive(Debug, Clone)]
pub struct EnergyTimeSeries {
    /// Time-indexed kinetic energy buffer.
    pub ke_series: TimeSeriesBuffer,
    /// Time-indexed potential energy buffer.
    pub pe_series: TimeSeriesBuffer,
    /// Maximum retention window in seconds.
    pub window_s: f64,
}
impl EnergyTimeSeries {
    /// Create an energy time series with the given time window (seconds).
    ///
    /// Capacity is set to `window_s * 1000` samples (1 ms resolution).
    pub fn new(window_s: f64) -> Self {
        let cap = (window_s * 1000.0).max(100.0) as usize;
        Self {
            ke_series: TimeSeriesBuffer::new(cap),
            pe_series: TimeSeriesBuffer::new(cap),
            window_s,
        }
    }
    /// Record a step with kinetic and potential energy at simulation time `t`.
    pub fn record(&mut self, t: f64, ke: f64, pe: f64) {
        self.ke_series.push(t, ke);
        self.pe_series.push(t, pe);
    }
    /// The N most recent kinetic energy samples.
    pub fn recent_ke(&self, n: usize) -> Vec<(f64, f64)> {
        let all: Vec<(f64, f64)> = self.ke_series.data().iter().copied().collect();
        let start = all.len().saturating_sub(n);
        all[start..].to_vec()
    }
    /// Total mechanical energy (KE + PE) at each recorded time.
    pub fn total_series(&self) -> Vec<(f64, f64)> {
        let kes: Vec<(f64, f64)> = self.ke_series.data().iter().copied().collect();
        let pes: Vec<(f64, f64)> = self.pe_series.data().iter().copied().collect();
        kes.iter()
            .zip(pes.iter())
            .map(|(&(t, ke), &(_tp, pe))| (t, ke + pe))
            .collect()
    }
    /// Energy drift: (E_last − E_first) / E_first (relative change).
    pub fn relative_drift(&self) -> f64 {
        let totals = self.total_series();
        if totals.len() < 2 {
            return 0.0;
        }
        let e0 = totals[0].1;
        let e1 = totals.last().expect("collection should not be empty").1;
        if e0.abs() < 1e-15 {
            return 0.0;
        }
        (e1 - e0) / e0.abs()
    }
}
/// A panel in the dashboard grid.
#[derive(Debug, Clone)]
pub struct DashboardPanel {
    /// Panel title.
    pub title: String,
    /// Grid position.
    pub position: PanelPosition,
    /// Content type.
    pub content: PanelContent,
    /// Background color (RGBA, 0–255).
    pub bg_color: [u8; 4],
    /// Border color.
    pub border_color: [u8; 4],
    /// Whether the panel is visible.
    pub visible: bool,
}
impl DashboardPanel {
    /// Create a titled panel at position.
    pub fn new(title: &str, position: PanelPosition, content: PanelContent) -> Self {
        Self {
            title: title.to_string(),
            position,
            content,
            bg_color: [20, 20, 30, 255],
            border_color: [80, 80, 100, 255],
            visible: true,
        }
    }
    /// Create an energy plot panel.
    pub fn energy_panel(row: usize, col: usize) -> Self {
        Self::new(
            "Energy (KE/PE/Total)",
            PanelPosition::single(row, col),
            PanelContent::EnergyPlot,
        )
    }
    /// Create a velocity histogram panel.
    pub fn velocity_panel(row: usize, col: usize) -> Self {
        Self::new(
            "Velocity Distribution",
            PanelPosition::single(row, col),
            PanelContent::VelocityHistogram,
        )
    }
}
/// Full dashboard configuration: panel list, theme, update rate.
#[derive(Debug, Clone)]
pub struct DashboardConfig {
    /// Ordered list of panel configurations.
    pub panels: Vec<PanelConfig>,
    /// Target update rate for real-time display (Hz).
    pub update_hz: f64,
    /// Dark/light theme selector.
    pub dark_theme: bool,
}
impl DashboardConfig {
    /// Create a standard 4-panel physics dashboard configuration.
    ///
    /// Panels: `"energy"`, `"momentum"`, `"performance"`, `"contacts"`.
    pub fn default_4panel() -> Self {
        let panels = vec![
            PanelConfig::new("energy", 400, 300),
            PanelConfig::new("momentum", 400, 300),
            PanelConfig::new("performance", 400, 300),
            PanelConfig::new("contacts", 400, 300),
        ];
        Self {
            panels,
            update_hz: 60.0,
            dark_theme: true,
        }
    }
    /// Set visibility of the named panel (no-op if not found).
    pub fn set_panel_visible(&mut self, name: &str, visible: bool) {
        if let Some(p) = self.panels.iter_mut().find(|p| p.name == name) {
            p.visible = visible;
        }
    }
    /// Returns `true` if the named panel exists and is visible.
    pub fn is_panel_visible(&self, name: &str) -> bool {
        self.panels.iter().any(|p| p.name == name && p.visible)
    }
    /// Number of currently visible panels.
    pub fn visible_count(&self) -> usize {
        self.panels.iter().filter(|p| p.visible).count()
    }
    /// Add a new panel configuration.
    pub fn add_panel(&mut self, cfg: PanelConfig) {
        self.panels.push(cfg);
    }
}
/// A multi-panel display layout, dividing screen space into a regular grid.
///
/// Each panel occupies one cell; the grid can be filled with different
/// visualization types (energy, phase space, heatmap, etc.).
#[derive(Debug, Clone)]
pub struct MultiPanelLayout {
    /// All panel cells in row-major order.
    pub panels: Vec<PanelCell>,
    /// Number of rows in the grid.
    pub num_rows: usize,
    /// Number of columns in the grid.
    pub num_cols: usize,
    /// Total canvas width (px).
    pub canvas_w: f64,
    /// Total canvas height (px).
    pub canvas_h: f64,
}
impl MultiPanelLayout {
    /// Build a regular `rows × cols` grid over the given canvas dimensions.
    pub fn grid(rows: usize, cols: usize, canvas_w: f64, canvas_h: f64) -> Self {
        let rows = rows.max(1);
        let cols = cols.max(1);
        let cw = canvas_w / cols as f64;
        let ch = canvas_h / rows as f64;
        let panels: Vec<PanelCell> = (0..rows)
            .flat_map(|r| {
                (0..cols).map(move |c| PanelCell {
                    row: r,
                    col: c,
                    x_px: c as f64 * cw,
                    y_px: r as f64 * ch,
                    w_px: cw,
                    h_px: ch,
                    name: String::new(),
                })
            })
            .collect();
        Self {
            panels,
            num_rows: rows,
            num_cols: cols,
            canvas_w,
            canvas_h,
        }
    }
    /// Find the panel at screen position (x, y).
    pub fn panel_at(&self, x: f64, y: f64) -> Option<&PanelCell> {
        self.panels.iter().find(|p| p.contains(x, y))
    }
    /// Assign a name to the panel at grid position (row, col).
    pub fn set_name(&mut self, row: usize, col: usize, name: &str) {
        if let Some(p) = self
            .panels
            .iter_mut()
            .find(|p| p.row == row && p.col == col)
        {
            p.name = name.to_string();
        }
    }
    /// Find a panel by name.
    pub fn panel_named(&self, name: &str) -> Option<&PanelCell> {
        self.panels.iter().find(|p| p.name == name)
    }
}
/// Energy panel: stacked area chart with anomaly detection.
#[derive(Debug, Clone)]
pub struct EnergyPanel {
    /// History of energy samples.
    pub samples: Vec<EnergySample>,
    /// Maximum samples to retain.
    pub max_samples: usize,
    /// Energy anomaly threshold (fraction increase per step).
    pub anomaly_threshold: f64,
    /// Number of energy anomalies detected.
    pub anomaly_count: u32,
}
impl EnergyPanel {
    /// Create a new energy panel.
    pub fn new(max_samples: usize, anomaly_threshold: f64) -> Self {
        Self {
            samples: Vec::new(),
            max_samples,
            anomaly_threshold,
            anomaly_count: 0,
        }
    }
    /// Add an energy sample.
    pub fn add_sample(&mut self, sample: EnergySample) {
        if let Some(prev) = self.samples.last()
            && prev.total() > 1e-10
        {
            let delta = (sample.total() - prev.total()) / prev.total();
            if delta > self.anomaly_threshold {
                self.anomaly_count += 1;
            }
        }
        self.samples.push(sample);
        if self.samples.len() > self.max_samples {
            self.samples.remove(0);
        }
    }
    /// Mean kinetic energy over the sample window.
    pub fn mean_ke(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().map(|s| s.ke).sum::<f64>() / self.samples.len() as f64
    }
    /// Maximum total energy in the window.
    pub fn max_total_energy(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.total())
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum total energy in the window.
    pub fn min_total_energy(&self) -> f64 {
        self.samples
            .iter()
            .map(|s| s.total())
            .fold(f64::INFINITY, f64::min)
    }
}
/// State of a joint constraint for inspection.
#[derive(Debug, Clone)]
pub struct ConstraintState {
    /// Constraint index.
    pub index: usize,
    /// Joint type.
    pub joint_type: JointType,
    /// Current angle or displacement.
    pub current_value: f64,
    /// Current constraint force/torque.
    pub current_force: f64,
    /// Lower limit.
    pub limit_lo: f64,
    /// Upper limit.
    pub limit_hi: f64,
    /// Motor target speed.
    pub motor_speed: f64,
    /// Whether motor is active.
    pub motor_active: bool,
}
impl ConstraintState {
    /// Whether the joint is at its lower limit.
    pub fn at_lower_limit(&self) -> bool {
        self.current_value <= self.limit_lo + 1e-6
    }
    /// Whether the joint is at its upper limit.
    pub fn at_upper_limit(&self) -> bool {
        self.current_value >= self.limit_hi - 1e-6
    }
    /// Normalized position within limits `[0, 1]`.
    pub fn normalized_position(&self) -> f64 {
        let range = self.limit_hi - self.limit_lo;
        if range < 1e-15 {
            return 0.5;
        }
        ((self.current_value - self.limit_lo) / range).clamp(0.0, 1.0)
    }
}
/// Constraint inspector panel.
#[derive(Debug, Clone, Default)]
pub struct ConstraintInspector {
    /// Currently selected constraint state.
    pub selected: Option<ConstraintState>,
}
impl ConstraintInspector {
    /// Create a new constraint inspector.
    pub fn new() -> Self {
        Self::default()
    }
    /// Select a constraint.
    pub fn select(&mut self, state: ConstraintState) {
        self.selected = Some(state);
    }
}
/// Physics anomaly alert system.
#[derive(Debug, Clone)]
pub struct AlertSystem {
    /// All alerts raised.
    pub alerts: Vec<PhysicsAlert>,
    /// Energy spike threshold (fractional increase per step).
    pub energy_spike_threshold: f64,
    /// Maximum allowed penetration depth (m).
    pub max_penetration_depth: f64,
    /// Previous total energy (for spike detection).
    pub prev_energy: f64,
}
impl AlertSystem {
    /// Create a new alert system.
    pub fn new(energy_spike_threshold: f64, max_penetration_depth: f64) -> Self {
        Self {
            alerts: Vec::new(),
            energy_spike_threshold,
            max_penetration_depth,
            prev_energy: 0.0,
        }
    }
    /// Check energy for spikes and raise an alert if necessary.
    pub fn check_energy(&mut self, total_energy: f64, sim_time: f64) {
        if self.prev_energy > 1e-10 {
            let delta = (total_energy - self.prev_energy) / self.prev_energy;
            if delta > self.energy_spike_threshold {
                self.alerts.push(PhysicsAlert::new(
                    format!(
                        "Energy spike: {:.1}% increase at t={:.3}s",
                        delta * 100.0,
                        sim_time
                    ),
                    AlertSeverity::Warning,
                    sim_time,
                ));
            }
        }
        self.prev_energy = total_energy;
    }
    /// Check penetration depth and raise an error if it exceeds threshold.
    pub fn check_penetration(&mut self, depth: f64, sim_time: f64) {
        if depth > self.max_penetration_depth {
            self.alerts.push(PhysicsAlert::new(
                format!(
                    "Penetration depth {:.4}m exceeds threshold {:.4}m",
                    depth, self.max_penetration_depth
                ),
                AlertSeverity::Error,
                sim_time,
            ));
        }
    }
    /// Acknowledge all alerts.
    pub fn acknowledge_all(&mut self) {
        for a in &mut self.alerts {
            a.acknowledged = true;
        }
    }
    /// Count of unacknowledged alerts at or above the given severity.
    pub fn unacknowledged_count(&self, min_severity: AlertSeverity) -> usize {
        self.alerts
            .iter()
            .filter(|a| !a.acknowledged && a.severity >= min_severity)
            .count()
    }
    /// Clear all acknowledged alerts.
    pub fn clear_acknowledged(&mut self) {
        self.alerts.retain(|a| !a.acknowledged);
    }
}
/// Panel content type.
#[derive(Debug, Clone, PartialEq)]
pub enum PanelContent {
    /// Energy time series plot.
    EnergyPlot,
    /// Velocity histogram.
    VelocityHistogram,
    /// Pressure/temperature history.
    ThermoHistory,
    /// Stress heatmap.
    StressHeatmap,
    /// Convergence plot.
    ConvergencePlot,
    /// Parameter sweep scatter.
    ParamSweepViz,
    /// Custom text label.
    Label(String),
    /// Empty panel.
    Empty,
}
/// Ring buffer for simulation snapshots.
#[derive(Debug, Clone)]
pub struct SnapshotBuffer {
    /// Circular buffer of snapshots.
    pub(super) buffer: Vec<SimSnapshot>,
    /// Maximum capacity.
    pub(super) capacity: usize,
    /// Current write head.
    pub(super) head: usize,
    /// Current count.
    pub(super) count: usize,
}
impl SnapshotBuffer {
    /// Create buffer with given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
            head: 0,
            count: 0,
        }
    }
    /// Push a new snapshot.
    pub fn push(&mut self, snap: SimSnapshot) {
        if self.buffer.len() < self.capacity {
            self.buffer.push(snap);
        } else {
            self.buffer[self.head] = snap;
        }
        self.head = (self.head + 1) % self.capacity;
        self.count = (self.count + 1).min(self.capacity);
    }
    /// Get snapshot by age (0 = most recent).
    pub fn get_by_age(&self, age: usize) -> Option<&SimSnapshot> {
        if age >= self.count {
            return None;
        }
        let actual_cap = self.buffer.len();
        if actual_cap == 0 {
            return None;
        }
        let idx = (self.head + self.capacity - 1 - age) % actual_cap;
        self.buffer.get(idx)
    }
    /// Ordered slice of all snapshots (oldest first).
    pub fn ordered(&self) -> Vec<&SimSnapshot> {
        let n = self.count;
        (0..n)
            .rev()
            .filter_map(|age| self.get_by_age(age))
            .collect()
    }
    /// Current number of snapshots.
    pub fn len(&self) -> usize {
        self.count
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
    /// Most recent snapshot.
    pub fn latest(&self) -> Option<&SimSnapshot> {
        self.get_by_age(0)
    }
    /// Energy drift: (E_latest - E_initial) / E_initial.
    pub fn energy_drift(&self) -> Option<f64> {
        let ordered = self.ordered();
        if ordered.len() < 2 {
            return None;
        }
        let e0 = ordered[0].total_energy();
        let e1 = ordered[ordered.len() - 1].total_energy();
        if e0.abs() < 1e-30 {
            return None;
        }
        Some((e1 - e0) / e0.abs())
    }
}
/// Material properties for inspection.
#[derive(Debug, Clone)]
pub struct MaterialState {
    /// Material name.
    pub name: String,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson ratio.
    pub poisson_ratio: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Yield stress (Pa).
    pub yield_stress: f64,
    /// Stress-strain curve data points.
    pub stress_strain_curve: Vec<StressStrainPoint>,
}
impl MaterialState {
    /// Create a new material state.
    pub fn new(
        name: impl Into<String>,
        youngs_modulus: f64,
        poisson_ratio: f64,
        density: f64,
    ) -> Self {
        Self {
            name: name.into(),
            youngs_modulus,
            poisson_ratio,
            density,
            yield_stress: youngs_modulus * 0.001,
            stress_strain_curve: Vec::new(),
        }
    }
    /// Generate a linear elastic stress-strain curve with `n` points.
    pub fn generate_linear_curve(&mut self, max_strain: f64, n: usize) {
        self.stress_strain_curve = (0..n)
            .map(|i| {
                let e = max_strain * i as f64 / (n - 1) as f64;
                let s = self.youngs_modulus * e;
                StressStrainPoint {
                    strain: e,
                    stress: s,
                }
            })
            .collect();
    }
    /// Interpolate stress at a given strain.
    pub fn stress_at_strain(&self, strain: f64) -> f64 {
        if self.stress_strain_curve.is_empty() {
            return self.youngs_modulus * strain;
        }
        let last = self
            .stress_strain_curve
            .last()
            .expect("collection should not be empty");
        if strain >= last.strain {
            return last.stress;
        }
        let idx = self
            .stress_strain_curve
            .partition_point(|p| p.strain <= strain);
        if idx == 0 {
            return 0.0;
        }
        let lo = &self.stress_strain_curve[idx - 1];
        let hi = &self.stress_strain_curve[idx];
        let t = (strain - lo.strain) / (hi.strain - lo.strain);
        lo.stress + (hi.stress - lo.stress) * t
    }
}
/// Performance monitoring panel.
#[derive(Debug, Clone)]
pub struct PerformancePanel {
    /// Frame timing history.
    pub frames: Vec<FrameTiming>,
    /// Maximum frame history size.
    pub max_frames: usize,
    /// Peak memory usage (bytes).
    pub peak_memory: usize,
}
impl PerformancePanel {
    /// Create a new performance panel.
    pub fn new(max_frames: usize) -> Self {
        Self {
            frames: Vec::new(),
            max_frames,
            peak_memory: 0,
        }
    }
    /// Record a frame timing.
    pub fn record_frame(&mut self, timing: FrameTiming) {
        if timing.memory_bytes > self.peak_memory {
            self.peak_memory = timing.memory_bytes;
        }
        self.frames.push(timing);
        if self.frames.len() > self.max_frames {
            self.frames.remove(0);
        }
    }
    /// Average total frame time over the window (ms).
    pub fn avg_total_ms(&self) -> f64 {
        if self.frames.is_empty() {
            return 0.0;
        }
        self.frames.iter().map(|f| f.total_ms()).sum::<f64>() / self.frames.len() as f64
    }
    /// Average frames per second.
    pub fn avg_fps(&self) -> f64 {
        let ms = self.avg_total_ms();
        if ms < 1e-10 {
            f64::INFINITY
        } else {
            1000.0 / ms
        }
    }
    /// Fraction of time spent in each phase (broadphase, narrowphase, solver, integration).
    pub fn phase_fractions(&self) -> [f64; 4] {
        if self.frames.is_empty() {
            return [0.25; 4];
        }
        let n = self.frames.len() as f64;
        let bp = self.frames.iter().map(|f| f.broadphase_ms).sum::<f64>() / n;
        let np = self.frames.iter().map(|f| f.narrowphase_ms).sum::<f64>() / n;
        let sv = self.frames.iter().map(|f| f.solver_ms).sum::<f64>() / n;
        let it = self.frames.iter().map(|f| f.integration_ms).sum::<f64>() / n;
        let total = bp + np + sv + it;
        if total < 1e-10 {
            return [0.25; 4];
        }
        [bp / total, np / total, sv / total, it / total]
    }
}
/// A single physics alert.
#[derive(Debug, Clone)]
pub struct PhysicsAlert {
    /// Alert message.
    pub message: String,
    /// Severity level.
    pub severity: AlertSeverity,
    /// Simulation time when alert was raised.
    pub time: f64,
    /// Whether the alert has been acknowledged.
    pub acknowledged: bool,
}
impl PhysicsAlert {
    /// Create a new alert.
    pub fn new(message: impl Into<String>, severity: AlertSeverity, time: f64) -> Self {
        Self {
            message: message.into(),
            severity,
            time,
            acknowledged: false,
        }
    }
}
/// Contact panel: contact count time series and force distribution.
#[derive(Debug, Clone)]
pub struct ContactPanel {
    /// Contact history.
    pub samples: Vec<ContactSample>,
    /// Maximum samples.
    pub max_samples: usize,
    /// Force histogram bins (bin counts).
    pub histogram: Vec<u32>,
    /// Maximum force for histogram range.
    pub hist_max_force: f64,
}
impl ContactPanel {
    /// Create a new contact panel.
    pub fn new(max_samples: usize, num_bins: usize, hist_max_force: f64) -> Self {
        Self {
            samples: Vec::new(),
            max_samples,
            histogram: vec![0u32; num_bins],
            hist_max_force,
        }
    }
    /// Add a contact sample and update histogram.
    pub fn add_sample(&mut self, sample: ContactSample) {
        let nb = self.histogram.len();
        if nb > 0 && sample.max_force > 0.0 {
            let bin = ((sample.mean_force / self.hist_max_force) * nb as f64) as usize;
            let bin = bin.min(nb - 1);
            self.histogram[bin] += 1;
        }
        self.samples.push(sample);
        if self.samples.len() > self.max_samples {
            self.samples.remove(0);
        }
    }
    /// Mean contact count over the window.
    pub fn mean_contact_count(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        self.samples.iter().map(|s| s.count as f64).sum::<f64>() / self.samples.len() as f64
    }
    /// Clear histogram.
    pub fn clear_histogram(&mut self) {
        for b in &mut self.histogram {
            *b = 0;
        }
    }
}
/// Overlay comparison of two simulation runs.
#[derive(Debug, Clone)]
pub struct SimCompare {
    /// Reference trace (run A).
    pub trace_a: SimTrace,
    /// Comparison trace (run B).
    pub trace_b: SimTrace,
}
impl SimCompare {
    /// Create a new comparison.
    pub fn new(trace_a: SimTrace, trace_b: SimTrace) -> Self {
        Self { trace_a, trace_b }
    }
    /// Compute delta trace (B - A) at times from trace_a.
    pub fn delta_trace(&self) -> Vec<(f64, f64)> {
        self.trace_a
            .times
            .iter()
            .zip(self.trace_a.values.iter())
            .filter_map(|(&t, &va)| {
                let vb = self.trace_b.sample_at(t)?;
                Some((t, vb - va))
            })
            .collect()
    }
    /// Mean absolute difference.
    pub fn mean_abs_diff(&self) -> f64 {
        let deltas = self.delta_trace();
        if deltas.is_empty() {
            return 0.0;
        }
        deltas.iter().map(|(_, d)| d.abs()).sum::<f64>() / deltas.len() as f64
    }
    /// Root mean square difference.
    pub fn rms_diff(&self) -> f64 {
        let deltas = self.delta_trace();
        if deltas.is_empty() {
            return 0.0;
        }
        (deltas.iter().map(|(_, d)| d * d).sum::<f64>() / deltas.len() as f64).sqrt()
    }
}
/// Multi-run parameter sweep visualisation.
#[derive(Debug, Clone)]
pub struct ParameterSweep {
    /// Sweep results.
    pub results: Vec<SweepResult>,
    /// Parameter names.
    pub param_names: Vec<String>,
    /// Metric names.
    pub metric_names: Vec<String>,
}
impl ParameterSweep {
    /// Create a new parameter sweep study.
    pub fn new(param_names: Vec<String>, metric_names: Vec<String>) -> Self {
        Self {
            results: Vec::new(),
            param_names,
            metric_names,
        }
    }
    /// Add a run result.
    pub fn add_result(&mut self, result: SweepResult) {
        self.results.push(result);
    }
    /// Return a scatter matrix entry: (param_i, param_j) → list of (xi, xj) pairs.
    pub fn scatter_pair(&self, param_i: usize, param_j: usize) -> Vec<(f64, f64)> {
        self.results
            .iter()
            .filter_map(|r| {
                let xi = r.parameters.get(param_i).copied()?;
                let xj = r.parameters.get(param_j).copied()?;
                Some((xi, xj))
            })
            .collect()
    }
    /// Compute a metric heatmap cell value (mean of metric `m` for runs with param_i ≈ val).
    pub fn heatmap_value(
        &self,
        param_i: usize,
        val: f64,
        metric_m: usize,
        tol: f64,
    ) -> Option<f64> {
        let matching: Vec<f64> = self
            .results
            .iter()
            .filter_map(|r| {
                let p = r.parameters.get(param_i).copied()?;
                if (p - val).abs() <= tol {
                    r.metrics.get(metric_m).copied()
                } else {
                    None
                }
            })
            .collect();
        if matching.is_empty() {
            return None;
        }
        Some(matching.iter().sum::<f64>() / matching.len() as f64)
    }
}
/// Performance profiler: rolling window of `StepBreakdown` records.
///
/// Computes average breakdowns and estimated FPS from the rolling window.
#[derive(Debug, Clone)]
pub struct PerformanceProfiler {
    /// Ring buffer of step breakdowns.
    pub(super) records: std::collections::VecDeque<StepBreakdown>,
    /// Maximum history size.
    pub(super) capacity: usize,
}
impl PerformanceProfiler {
    /// Create a performance profiler with the given history capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            records: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
        }
    }
    /// Record one step's timing breakdown.
    pub fn record(&mut self, bd: StepBreakdown) {
        if self.records.len() >= self.capacity {
            self.records.pop_front();
        }
        self.records.push_back(bd);
    }
    /// Number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }
    /// Returns `true` if no records have been stored.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    /// Average step breakdown over the rolling window.
    pub fn average_breakdown(&self) -> StepBreakdown {
        if self.records.is_empty() {
            return StepBreakdown::default();
        }
        let n = self.records.len() as u64;
        StepBreakdown {
            broadphase_us: self.records.iter().map(|r| r.broadphase_us).sum::<u64>() / n,
            narrowphase_us: self.records.iter().map(|r| r.narrowphase_us).sum::<u64>() / n,
            solver_us: self.records.iter().map(|r| r.solver_us).sum::<u64>() / n,
            integration_us: self.records.iter().map(|r| r.integration_us).sum::<u64>() / n,
            render_us: self.records.iter().map(|r| r.render_us).sum::<u64>() / n,
        }
    }
    /// Estimated FPS from the mean total step time.
    pub fn estimated_fps(&self) -> f64 {
        let avg = self.average_breakdown();
        let total_us = avg.total_us();
        if total_us == 0 {
            return f64::INFINITY;
        }
        1_000_000.0 / total_us as f64
    }
    /// Peak total step time (µs) in the window.
    pub fn peak_total_us(&self) -> u64 {
        self.records.iter().map(|r| r.total_us()).max().unwrap_or(0)
    }
}
/// Solver convergence history tracker.
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor {
    /// Residual history.
    pub residuals: Vec<f64>,
    /// Iteration labels.
    pub iterations: Vec<usize>,
    /// Convergence criterion.
    pub criterion: ConvergenceCriterion,
    /// Whether solver has converged.
    pub converged: bool,
    /// Convergence iteration.
    pub converge_iter: Option<usize>,
    /// Initial residual.
    pub initial_residual: f64,
    /// Maximum allowed iterations.
    pub max_iter: usize,
}
impl ConvergenceMonitor {
    /// Create with given criterion and maximum iterations.
    pub fn new(criterion: ConvergenceCriterion, max_iter: usize) -> Self {
        Self {
            residuals: Vec::new(),
            iterations: Vec::new(),
            criterion,
            converged: false,
            converge_iter: None,
            initial_residual: 0.0,
            max_iter,
        }
    }
    /// Record a residual value and check convergence.
    pub fn record(&mut self, residual: f64) -> bool {
        let iter = self.residuals.len();
        if iter == 0 {
            self.initial_residual = residual;
        }
        self.residuals.push(residual);
        self.iterations.push(iter);
        if !self.converged && self.criterion.satisfied(residual, self.initial_residual) {
            self.converged = true;
            self.converge_iter = Some(iter);
        }
        self.converged
    }
    /// Convergence rate (slope of log residual per iteration).
    pub fn convergence_rate(&self) -> Option<f64> {
        let n = self.residuals.len();
        if n < 2 {
            return None;
        }
        let r0 = self.residuals[0];
        let r1 = self.residuals[n - 1];
        if r0 <= 0.0 || r1 <= 0.0 {
            return None;
        }
        Some((r1.ln() - r0.ln()) / (n - 1) as f64)
    }
    /// Estimated iterations to reach tolerance `tol`.
    pub fn estimated_iters_to(&self, tol: f64) -> Option<usize> {
        let rate = self.convergence_rate()?;
        if rate >= 0.0 {
            return None;
        }
        let last_r = self.residuals.last()?;
        if *last_r <= tol {
            return Some(0);
        }
        let n_more = ((tol.ln() - last_r.ln()) / rate).ceil() as usize;
        Some(n_more)
    }
    /// Generate ASCII sparkline of residual history.
    pub fn sparkline(&self) -> String {
        let chars = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
        if self.residuals.is_empty() {
            return String::new();
        }
        let max_r = self
            .residuals
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_r = self.residuals.iter().cloned().fold(f64::INFINITY, f64::min);
        let range = (max_r - min_r).max(1e-30);
        self.residuals
            .iter()
            .map(|&r| {
                let t = (r - min_r) / range;
                let idx = (t * 7.0).round() as usize;
                chars[idx.min(7)]
            })
            .collect()
    }
}
/// Material inspector panel.
#[derive(Debug, Clone, Default)]
pub struct MaterialInspector {
    /// Currently selected element's material.
    pub selected: Option<MaterialState>,
}
impl MaterialInspector {
    /// Create a new material inspector.
    pub fn new() -> Self {
        Self::default()
    }
    /// Select a material for inspection.
    pub fn select(&mut self, state: MaterialState) {
        self.selected = Some(state);
    }
}
/// Collection of parameter sweep results.
#[derive(Debug, Clone, Default)]
pub struct ParamSweepViz {
    /// Parameter name.
    pub param_name: String,
    /// Result name.
    pub result_name: String,
    /// Entries sorted by param_value.
    pub entries: Vec<SweepEntry>,
}
impl ParamSweepViz {
    /// Create empty sweep.
    pub fn new(param_name: &str, result_name: &str) -> Self {
        Self {
            param_name: param_name.to_string(),
            result_name: result_name.to_string(),
            entries: Vec::new(),
        }
    }
    /// Add an entry.
    pub fn add(&mut self, param_value: f64, result: f64, error: Option<f64>, converged: bool) {
        self.entries.push(SweepEntry {
            param_value,
            result,
            error,
            converged,
        });
        self.entries.sort_by(|a, b| {
            a.param_value
                .partial_cmp(&b.param_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Parameter values as Vec.
    pub fn param_values(&self) -> Vec<f64> {
        self.entries.iter().map(|e| e.param_value).collect()
    }
    /// Result values as Vec.
    pub fn result_values(&self) -> Vec<f64> {
        self.entries.iter().map(|e| e.result).collect()
    }
    /// Find optimal parameter value (minimizes result).
    pub fn optimal_min(&self) -> Option<f64> {
        self.entries
            .iter()
            .filter(|e| e.converged)
            .min_by(|a, b| {
                a.result
                    .partial_cmp(&b.result)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|e| e.param_value)
    }
    /// Find optimal parameter value (maximizes result).
    pub fn optimal_max(&self) -> Option<f64> {
        self.entries
            .iter()
            .filter(|e| e.converged)
            .max_by(|a, b| {
                a.result
                    .partial_cmp(&b.result)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|e| e.param_value)
    }
    /// Convergence fraction.
    pub fn convergence_fraction(&self) -> f64 {
        if self.entries.is_empty() {
            return 0.0;
        }
        self.entries.iter().filter(|e| e.converged).count() as f64 / self.entries.len() as f64
    }
}
/// Multi-panel simulation dashboard layout.
#[derive(Debug, Clone, Default)]
pub struct DashboardLayout {
    /// Grid rows.
    pub rows: usize,
    /// Grid columns.
    pub cols: usize,
    /// Panels in the layout.
    pub panels: Vec<DashboardPanel>,
    /// Dashboard title.
    pub title: String,
    /// Background color.
    pub bg_color: [u8; 4],
}
impl DashboardLayout {
    /// Create a layout grid.
    pub fn new(rows: usize, cols: usize, title: &str) -> Self {
        Self {
            rows,
            cols,
            panels: Vec::new(),
            title: title.to_string(),
            bg_color: [10, 10, 15, 255],
        }
    }
    /// Add a panel (if within bounds).
    pub fn add_panel(&mut self, panel: DashboardPanel) -> bool {
        if panel.position.row + panel.position.row_span <= self.rows
            && panel.position.col + panel.position.col_span <= self.cols
        {
            self.panels.push(panel);
            true
        } else {
            false
        }
    }
    /// Standard physics dashboard with 4 panels (2×2).
    pub fn standard_4panel(sim_name: &str) -> Self {
        let mut layout = Self::new(2, 2, sim_name);
        layout.add_panel(DashboardPanel::energy_panel(0, 0));
        layout.add_panel(DashboardPanel::velocity_panel(0, 1));
        layout.add_panel(DashboardPanel::new(
            "Convergence",
            PanelPosition::single(1, 0),
            PanelContent::ConvergencePlot,
        ));
        layout.add_panel(DashboardPanel::new(
            "Thermodynamics",
            PanelPosition::single(1, 1),
            PanelContent::ThermoHistory,
        ));
        layout
    }
    /// Number of visible panels.
    pub fn visible_count(&self) -> usize {
        self.panels.iter().filter(|p| p.visible).count()
    }
}
