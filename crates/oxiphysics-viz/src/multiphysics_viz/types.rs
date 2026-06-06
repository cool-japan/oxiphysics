//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{vec3_dot, vec3_norm, vec3_sub};

/// A keyframe in a temperature-color ramp.
#[derive(Debug, Clone)]
pub struct TempColorKey {
    /// Temperature (K).
    pub temperature: f64,
    /// Color at this temperature.
    pub color: Rgba,
}
impl TempColorKey {
    /// Construct a new `TempColorKey`.
    pub fn new(temperature: f64, color: Rgba) -> Self {
        Self { temperature, color }
    }
}
/// A multiphysics convergence dashboard aggregating multiple domains.
#[derive(Debug, Clone, Default)]
pub struct MultiphysicsConvergenceDashboard {
    /// Convergence records for each physics domain.
    pub domains: Vec<DomainConvergence>,
    /// Current coupling iteration number.
    pub iteration: usize,
    /// Maximum coupling iterations allowed.
    pub max_iterations: usize,
}
impl MultiphysicsConvergenceDashboard {
    /// Construct an empty dashboard.
    pub fn new(max_iterations: usize) -> Self {
        Self {
            domains: Vec::new(),
            max_iterations,
            iteration: 0,
        }
    }
    /// Register a new physics domain.
    pub fn add_domain(&mut self, domain: DomainConvergence) {
        self.domains.push(domain);
    }
    /// Advance the coupling iteration counter.
    pub fn next_iteration(&mut self) {
        self.iteration += 1;
    }
    /// Return `true` if all registered domains have converged.
    pub fn all_converged(&self) -> bool {
        !self.domains.is_empty() && self.domains.iter().all(|d| d.converged())
    }
    /// Return `true` if the iteration limit has been reached.
    pub fn iteration_limit_reached(&self) -> bool {
        self.iteration >= self.max_iterations
    }
    /// Return `true` if either all domains converged or the limit is reached.
    pub fn should_stop(&self) -> bool {
        self.all_converged() || self.iteration_limit_reached()
    }
    /// Build a human-readable summary string.
    pub fn summary(&self) -> String {
        let mut s = format!("Iteration {}/{}\n", self.iteration, self.max_iterations);
        for d in &self.domains {
            let last = d.residuals.last().copied().unwrap_or(f64::NAN);
            let status = if d.converged() {
                "CONVERGED"
            } else {
                "running"
            };
            s.push_str(&format!(
                "  {} residual={:.3e} tol={:.3e} [{}]\n",
                d.name, last, d.tolerance, status
            ));
        }
        s
    }
}
/// Configuration for a coupled thermal-mechanical overlay.
///
/// Blends a thermal field (temperature colormap) and a mechanical field
/// (displacement magnitude colormap) into a single RGBA color per node.
#[derive(Debug, Clone)]
pub struct ThermalMechanicalOverlayConfig {
    /// Minimum temperature for color mapping (K).
    pub temp_min: f64,
    /// Maximum temperature for color mapping (K).
    pub temp_max: f64,
    /// Minimum displacement magnitude for color mapping (m).
    pub disp_min: f64,
    /// Maximum displacement magnitude for color mapping (m).
    pub disp_max: f64,
    /// Blend weight for thermal channel in `[0, 1]`; mechanical weight = 1 - thermal_weight.
    pub thermal_weight: f64,
}
// Default impl is in trait_impls.rs to avoid duplication
/// A stack of clipping planes; a point is visible only if it passes all planes.
#[derive(Debug, Clone, Default)]
pub struct ClippingPlaneStack {
    /// Active clipping planes.
    pub planes: Vec<ClippingPlane>,
}
impl ClippingPlaneStack {
    /// Construct an empty stack.
    pub fn new() -> Self {
        Self { planes: Vec::new() }
    }
    /// Add a clipping plane.
    pub fn add(&mut self, plane: ClippingPlane) {
        self.planes.push(plane);
    }
    /// Return `true` if `point` is visible (passes all planes).
    pub fn visible(&self, point: [f64; 3]) -> bool {
        self.planes.iter().all(|p| p.visible(point))
    }
    /// Filter a list of nodes, returning those visible through all clipping planes.
    pub fn filter_nodes<'a>(&self, nodes: &'a [FieldNode]) -> Vec<&'a FieldNode> {
        nodes.iter().filter(|n| self.visible(n.position)).collect()
    }
    /// Filter line segments, retaining those where both endpoints are visible.
    pub fn filter_segments<'a>(&self, segs: &'a [LineSegment]) -> Vec<&'a LineSegment> {
        segs.iter()
            .filter(|s| self.visible(s.start) && self.visible(s.end))
            .collect()
    }
}
/// Synchronization policy when physics domains have different time steps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SyncPolicy {
    /// Use the coarsest time step (minimum frequency).
    Coarsest,
    /// Use the finest time step (maximum frequency).
    Finest,
    /// Use a user-specified fixed time step.
    Fixed(f64),
}
/// Interpolation method for cross-mesh field transfer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeshInterpolationMethod {
    /// Nearest-neighbor: use the value of the closest source node.
    NearestNeighbor,
    /// Inverse-distance weighted interpolation (IDW).
    InverseDistance {
        /// Power parameter for IDW (typically 2).
        power: f64,
        /// Number of nearest neighbors to use.
        k_neighbors: usize,
    },
    /// Bilinear interpolation within the source triangle (fallback to IDW).
    Barycentric,
}
/// A colored 3-D line segment for visualization output.
#[derive(Debug, Clone)]
pub struct LineSegment {
    /// Start point.
    pub start: [f64; 3],
    /// End point.
    pub end: [f64; 3],
    /// Color of the segment.
    pub color: Rgba,
}
impl LineSegment {
    /// Construct a new `LineSegment`.
    pub fn new(start: [f64; 3], end: [f64; 3], color: Rgba) -> Self {
        Self { start, end, color }
    }
}
/// Integration method for field-line tracing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FieldLineMethod {
    /// First-order Euler integration.
    Euler,
    /// Fourth-order Runge-Kutta integration.
    Rk4,
}
/// A piecewise-linear temperature → color ramp.
///
/// Keys must be sorted by ascending temperature; behavior is undefined otherwise.
/// Use [`TempColorRamp::sorted`] to ensure ordering.
#[derive(Debug, Clone)]
pub struct TempColorRamp {
    /// Sorted keyframes defining the ramp.
    pub keys: Vec<TempColorKey>,
}
impl TempColorRamp {
    /// Construct a ramp from unsorted keys (sorts them in place).
    pub fn sorted(mut keys: Vec<TempColorKey>) -> Self {
        keys.sort_by(|a, b| {
            a.temperature
                .partial_cmp(&b.temperature)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Self { keys }
    }
    /// Map a temperature to a color via piecewise-linear interpolation.
    pub fn evaluate(&self, temp: f64) -> Rgba {
        if self.keys.is_empty() {
            return Rgba::white();
        }
        if temp
            <= self
                .keys
                .first()
                .expect("collection should not be empty")
                .temperature
        {
            return self
                .keys
                .first()
                .expect("collection should not be empty")
                .color;
        }
        if temp
            >= self
                .keys
                .last()
                .expect("collection should not be empty")
                .temperature
        {
            return self
                .keys
                .last()
                .expect("collection should not be empty")
                .color;
        }
        for i in 0..self.keys.len() - 1 {
            let t0 = self.keys[i].temperature;
            let t1 = self.keys[i + 1].temperature;
            if temp >= t0 && temp <= t1 {
                let frac = ((temp - t0) / (t1 - t0)).clamp(0.0, 1.0) as f32;
                return Rgba::lerp(self.keys[i].color, self.keys[i + 1].color, frac);
            }
        }
        Rgba::white()
    }
    /// Generate per-node colors from a temperature field.
    pub fn colorize_temperatures(&self, temperatures: &[f64]) -> Vec<Rgba> {
        temperatures.iter().map(|&t| self.evaluate(t)).collect()
    }
}
/// An axis-aligned or arbitrarily-oriented clipping plane.
#[derive(Debug, Clone)]
pub struct ClippingPlane {
    /// A point on the plane.
    pub origin: [f64; 3],
    /// Outward normal of the plane (points toward the clipped-away side).
    pub normal: [f64; 3],
}
impl ClippingPlane {
    /// Construct a new `ClippingPlane`.
    pub fn new(origin: [f64; 3], normal: [f64; 3]) -> Self {
        let normal = vec3_norm(normal);
        Self { origin, normal }
    }
    /// Return the signed distance from `point` to the plane.
    /// Positive = in front (visible side); negative = behind (clipped side).
    pub fn signed_distance(&self, point: [f64; 3]) -> f64 {
        vec3_dot(vec3_sub(point, self.origin), self.normal)
    }
    /// Return `true` if `point` is on the visible side (signed distance >= 0).
    pub fn visible(&self, point: [f64; 3]) -> bool {
        self.signed_distance(point) >= 0.0
    }
}
/// A synchronized multiphysics timeline that aligns multiple domain time series.
#[derive(Debug, Clone)]
pub struct MultiphysicsTimeline {
    /// Registered physics domain time series.
    pub domains: Vec<PhysicsTimeSeries>,
    /// Synchronization policy.
    pub policy: SyncPolicy,
}
impl MultiphysicsTimeline {
    /// Construct an empty timeline with the given sync policy.
    pub fn new(policy: SyncPolicy) -> Self {
        Self {
            domains: Vec::new(),
            policy,
        }
    }
    /// Register a domain time series.
    pub fn add_domain(&mut self, series: PhysicsTimeSeries) {
        self.domains.push(series);
    }
    /// Compute the global start time (minimum across all domains).
    pub fn global_start(&self) -> f64 {
        self.domains
            .iter()
            .filter_map(|d| d.times.first().copied())
            .fold(f64::INFINITY, f64::min)
    }
    /// Compute the global end time (maximum across all domains).
    pub fn global_end(&self) -> f64 {
        self.domains
            .iter()
            .filter_map(|d| d.times.last().copied())
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Compute the synchronized time step under the current policy.
    pub fn sync_step(&self) -> f64 {
        match self.policy {
            SyncPolicy::Fixed(dt) => dt,
            SyncPolicy::Coarsest => self
                .domains
                .iter()
                .filter_map(|d| {
                    if d.times.len() < 2 {
                        return None;
                    }
                    let dts: Vec<f64> = d.times.windows(2).map(|w| w[1] - w[0]).collect();
                    dts.iter().copied().reduce(f64::max)
                })
                .fold(0.0_f64, f64::max)
                .max(1e-15),
            SyncPolicy::Finest => self
                .domains
                .iter()
                .filter_map(|d| {
                    if d.times.len() < 2 {
                        return None;
                    }
                    let dts: Vec<f64> = d.times.windows(2).map(|w| w[1] - w[0]).collect();
                    dts.iter().copied().reduce(f64::min)
                })
                .fold(f64::INFINITY, f64::min)
                .max(1e-15),
        }
    }
    /// Return interpolated values from all domains at `query_time`.
    pub fn sample_all(&self, query_time: f64) -> Vec<Option<f64>> {
        self.domains
            .iter()
            .map(|d| d.interpolate_at(query_time))
            .collect()
    }
    /// Generate the synchronized time grid as a `Vec`f64`.
    pub fn sync_grid(&self) -> Vec<f64> {
        let start = self.global_start();
        let end = self.global_end();
        let dt = self.sync_step();
        if start >= end || dt < 1e-15 {
            return vec![start];
        }
        let n = ((end - start) / dt).ceil() as usize + 1;
        (0..n).map(|i| (start + i as f64 * dt).min(end)).collect()
    }
}
/// Parameters for electromagnetic field-line tracing.
#[derive(Debug, Clone)]
pub struct EmFieldLineConfig {
    /// Integration step size (spatial units).
    pub step_size: f64,
    /// Maximum number of integration steps.
    pub max_steps: usize,
    /// Stop integration when field magnitude falls below this threshold.
    pub min_magnitude: f64,
    /// Integration method.
    pub method: FieldLineMethod,
}
// Default impl is in trait_impls.rs to avoid duplication
/// A single traced electromagnetic field line.
#[derive(Debug, Clone)]
pub struct EmFieldLine {
    /// Sequence of positions along the field line.
    pub points: Vec<[f64; 3]>,
    /// Field magnitude sampled at each point.
    pub magnitudes: Vec<f64>,
}
impl EmFieldLine {
    pub(super) fn new() -> Self {
        Self {
            points: Vec::new(),
            magnitudes: Vec::new(),
        }
    }
}
/// Convergence history for a single physics domain.
#[derive(Debug, Clone)]
pub struct DomainConvergence {
    /// Name of the physics domain (e.g. `"thermal"`, `"mechanical"`).
    pub name: String,
    /// Residual norm at each coupling iteration.
    pub residuals: Vec<f64>,
    /// Convergence tolerance.
    pub tolerance: f64,
}
impl DomainConvergence {
    /// Construct a new `DomainConvergence` record.
    pub fn new(name: impl Into<String>, tolerance: f64) -> Self {
        Self {
            name: name.into(),
            residuals: Vec::new(),
            tolerance,
        }
    }
    /// Append a new residual sample.
    pub fn push(&mut self, residual: f64) {
        self.residuals.push(residual);
    }
    /// Return `true` if the last residual is below the tolerance.
    pub fn converged(&self) -> bool {
        self.residuals.last().is_some_and(|&r| r <= self.tolerance)
    }
    /// Return the convergence rate (ratio of last two residuals), or `None` if fewer than 2 samples.
    pub fn rate(&self) -> Option<f64> {
        let n = self.residuals.len();
        if n < 2 {
            return None;
        }
        let prev = self.residuals[n - 2];
        if prev.abs() < 1e-300 {
            return None;
        }
        Some(self.residuals[n - 1] / prev)
    }
}
/// Configuration for coupled field arrow glyphs.
#[derive(Debug, Clone)]
pub struct CoupledArrowConfig {
    /// Scale factor for heat flux arrows.
    pub heat_flux_scale: f64,
    /// Scale factor for displacement arrows.
    pub displacement_scale: f64,
    /// Color for heat flux arrows.
    pub heat_flux_color: Rgba,
    /// Color for displacement arrows.
    pub displacement_color: Rgba,
    /// Minimum vector magnitude below which no arrow is drawn.
    pub min_magnitude: f64,
}
// Default impl is in trait_impls.rs to avoid duplication
/// A snapshot of FSI coupling state at one time step.
#[derive(Debug, Clone)]
pub struct FsiSnapshot {
    /// Simulation time (s).
    pub time: f64,
    /// Fluid node positions (deformed).
    pub fluid_positions: Vec<[f64; 3]>,
    /// Fluid velocity vectors.
    pub fluid_velocities: Vec<[f64; 3]>,
    /// Structural node positions (deformed).
    pub structural_positions: Vec<[f64; 3]>,
    /// Structural displacement vectors.
    pub structural_displacements: Vec<[f64; 3]>,
    /// Fluid-structure interface node indices in `fluid_positions`.
    pub fsi_interface_indices: Vec<usize>,
}
impl FsiSnapshot {
    /// Construct an empty `FsiSnapshot` at time `t`.
    pub fn new(time: f64) -> Self {
        Self {
            time,
            fluid_positions: Vec::new(),
            fluid_velocities: Vec::new(),
            structural_positions: Vec::new(),
            structural_displacements: Vec::new(),
            fsi_interface_indices: Vec::new(),
        }
    }
}
/// RGBA color with f32 components in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}
impl Rgba {
    /// Construct a new `Rgba` color.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
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
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }
    /// Opaque yellow.
    pub fn yellow() -> Self {
        Self::new(1.0, 1.0, 0.0, 1.0)
    }
    /// Opaque cyan.
    pub fn cyan() -> Self {
        Self::new(0.0, 1.0, 1.0, 1.0)
    }
    /// Opaque magenta.
    pub fn magenta() -> Self {
        Self::new(1.0, 0.0, 1.0, 1.0)
    }
    /// Linear interpolation between two colors.
    pub fn lerp(a: Self, b: Self, t: f32) -> Self {
        Self {
            r: a.r + t * (b.r - a.r),
            g: a.g + t * (b.g - a.g),
            b: a.b + t * (b.b - a.b),
            a: a.a + t * (b.a - a.a),
        }
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// A node in a physics mesh that carries scalar and vector field data.
#[derive(Debug, Clone)]
pub struct FieldNode {
    /// Spatial position.
    pub position: [f64; 3],
    /// Scalar value (temperature, pressure, density, etc.).
    pub scalar: f64,
    /// Vector value (velocity, displacement, heat flux, etc.).
    pub vector: [f64; 3],
}
impl FieldNode {
    /// Construct a new `FieldNode`.
    pub fn new(position: [f64; 3], scalar: f64, vector: [f64; 3]) -> Self {
        Self {
            position,
            scalar,
            vector,
        }
    }
}
/// A single phase-boundary sample point in 3-D.
#[derive(Debug, Clone)]
pub struct InterfacePoint {
    /// Spatial position of the interface sample.
    pub position: [f64; 3],
    /// Outward normal direction at this point.
    pub normal: [f64; 3],
    /// Phase indicator (0 = phase A, 1 = phase B).
    pub phase: u8,
}
impl InterfacePoint {
    /// Construct a new `InterfacePoint`.
    pub fn new(position: [f64; 3], normal: [f64; 3], phase: u8) -> Self {
        Self {
            position,
            normal,
            phase,
        }
    }
}
/// A time series record for one physics domain.
#[derive(Debug, Clone)]
pub struct PhysicsTimeSeries {
    /// Name of the domain.
    pub name: String,
    /// Time stamps of recorded frames (seconds).
    pub times: Vec<f64>,
    /// Scalar quantity sampled at each frame (e.g. max temperature, max stress).
    pub values: Vec<f64>,
}
impl PhysicsTimeSeries {
    /// Construct an empty time series.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            times: Vec::new(),
            values: Vec::new(),
        }
    }
    /// Append a `(time, value)` sample.
    pub fn push(&mut self, time: f64, value: f64) {
        self.times.push(time);
        self.values.push(value);
    }
    /// Interpolate the value at `query_time` using linear interpolation.
    pub fn interpolate_at(&self, query_time: f64) -> Option<f64> {
        let n = self.times.len();
        if n == 0 {
            return None;
        }
        if query_time <= self.times[0] {
            return Some(self.values[0]);
        }
        if query_time >= self.times[n - 1] {
            return Some(self.values[n - 1]);
        }
        for i in 0..n - 1 {
            if query_time >= self.times[i] && query_time <= self.times[i + 1] {
                let t = (query_time - self.times[i]) / (self.times[i + 1] - self.times[i]);
                return Some(self.values[i] + t * (self.values[i + 1] - self.values[i]));
            }
        }
        None
    }
}
