// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiphysics field overlay visualization for the OxiPhysics engine.
//!
//! This module provides utilities for visualizing coupled multiphysics simulations,
//! including:
//! - Coupled thermal-mechanical field overlays
//! - Coupled field arrow glyphs (heat flux + deformation vectors)
//! - Interface tracking and phase boundary visualization
//! - Fluid-structure interaction (FSI) visualization
//! - Electromagnetic field lines
//! - Temperature-dependent color transitions
//! - Multiphysics convergence dashboards
//! - Field interpolation between different meshes
//! - Clipping planes for multiphysics fields
//! - Timeline synchronization across physics domains

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Basic 3-D vector arithmetic helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-D vectors.
#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors.
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by a scalar.
#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-D vectors.
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean length of a 3-D vector.
#[inline]
fn vec3_len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-D vector; returns `[0,0,0]` for near-zero vectors.
#[inline]
fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

/// Cross product of two 3-D vectors.
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Linear interpolation between two 3-D vectors.
#[inline]
fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

// ---------------------------------------------------------------------------
// RGBA color helper
// ---------------------------------------------------------------------------

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

impl Default for Rgba {
    fn default() -> Self {
        Self::white()
    }
}

// ---------------------------------------------------------------------------
// LineSegment — a colored line primitive
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// FieldNode — a scalar or vector sample on a mesh node
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// 1. Coupled thermal-mechanical field overlay
// ---------------------------------------------------------------------------

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

impl Default for ThermalMechanicalOverlayConfig {
    fn default() -> Self {
        Self {
            temp_min: 0.0,
            temp_max: 1000.0,
            disp_min: 0.0,
            disp_max: 0.01,
            thermal_weight: 0.5,
        }
    }
}

/// Map a normalized `t` in `[0,1]` to a hot colormap (black → red → yellow → white).
fn hot_colormap(t: f64) -> Rgba {
    let t = t.clamp(0.0, 1.0) as f32;
    let r = (t * 3.0).min(1.0);
    let g = (t * 3.0 - 1.0).clamp(0.0, 1.0);
    let b = (t * 3.0 - 2.0).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}

/// Map a normalized `t` in `[0,1]` to a cool colormap (cyan → magenta).
fn cool_colormap(t: f64) -> Rgba {
    let t = t.clamp(0.0, 1.0) as f32;
    Rgba::new(t, 1.0 - t, 1.0, 1.0)
}

/// Normalize a scalar to `[0,1]` given min/max.
fn normalize_scalar(val: f64, min: f64, max: f64) -> f64 {
    if (max - min).abs() < 1e-15 {
        0.5
    } else {
        ((val - min) / (max - min)).clamp(0.0, 1.0)
    }
}

/// Compute per-node colors for a coupled thermal-mechanical overlay.
///
/// For each node the thermal color (hot colormap from temperature) and
/// mechanical color (cool colormap from displacement magnitude) are blended
/// according to `config.thermal_weight`.
pub fn thermal_mechanical_overlay(
    temperature_nodes: &[FieldNode],
    displacement_nodes: &[FieldNode],
    config: &ThermalMechanicalOverlayConfig,
) -> Vec<Rgba> {
    let n = temperature_nodes.len().min(displacement_nodes.len());
    let mut colors = Vec::with_capacity(n);
    for i in 0..n {
        let t_norm = normalize_scalar(
            temperature_nodes[i].scalar,
            config.temp_min,
            config.temp_max,
        );
        let d_norm = normalize_scalar(
            vec3_len(displacement_nodes[i].vector),
            config.disp_min,
            config.disp_max,
        );
        let thermal_color = hot_colormap(t_norm);
        let mech_color = cool_colormap(d_norm);
        let w = config.thermal_weight.clamp(0.0, 1.0) as f32;
        colors.push(Rgba::lerp(mech_color, thermal_color, w));
    }
    colors
}

// ---------------------------------------------------------------------------
// 2. Coupled field arrows (heat flux + deformation)
// ---------------------------------------------------------------------------

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

impl Default for CoupledArrowConfig {
    fn default() -> Self {
        Self {
            heat_flux_scale: 1.0,
            displacement_scale: 1.0,
            heat_flux_color: Rgba::new(1.0, 0.4, 0.0, 1.0),
            displacement_color: Rgba::new(0.2, 0.6, 1.0, 1.0),
            min_magnitude: 1e-12,
        }
    }
}

/// Generate arrow [`LineSegment`]s for two coupled vector fields at the same nodes.
///
/// Returns a flat list: for each node two segments are generated — one for
/// `heat_flux` (orange) and one for `displacement` (blue), filtered by
/// `config.min_magnitude`.
pub fn coupled_field_arrows(
    nodes: &[FieldNode],
    heat_flux: &[[f64; 3]],
    displacement: &[[f64; 3]],
    config: &CoupledArrowConfig,
) -> Vec<LineSegment> {
    let n = nodes.len().min(heat_flux.len()).min(displacement.len());
    let mut segments = Vec::with_capacity(n * 2);
    for i in 0..n {
        let origin = nodes[i].position;

        let hf = heat_flux[i];
        if vec3_len(hf) >= config.min_magnitude {
            let tip = vec3_add(origin, vec3_scale(hf, config.heat_flux_scale));
            segments.push(LineSegment::new(origin, tip, config.heat_flux_color));
        }

        let disp = displacement[i];
        if vec3_len(disp) >= config.min_magnitude {
            let tip = vec3_add(origin, vec3_scale(disp, config.displacement_scale));
            segments.push(LineSegment::new(origin, tip, config.displacement_color));
        }
    }
    segments
}

// ---------------------------------------------------------------------------
// 3. Interface tracking (phase boundaries)
// ---------------------------------------------------------------------------

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

/// Extract interface points from a scalar volume field using marching-squares
/// in 1-D stencil (simplified): a point is on the interface if its scalar
/// value straddles `iso_value` between adjacent samples.
///
/// `scalars` and `positions` must have the same length.  Returns the
/// interface samples sorted by position along the first coordinate axis.
pub fn extract_interface_points(
    positions: &[[f64; 3]],
    scalars: &[f64],
    iso_value: f64,
) -> Vec<InterfacePoint> {
    let n = positions.len().min(scalars.len());
    let mut points = Vec::new();
    if n < 2 {
        return points;
    }

    for i in 0..n - 1 {
        let s0 = scalars[i];
        let s1 = scalars[i + 1];
        let crosses = (s0 < iso_value) != (s1 < iso_value);
        if crosses {
            let t = (iso_value - s0) / (s1 - s0 + 1e-300);
            let pos = vec3_lerp(positions[i], positions[i + 1], t);
            let dir = vec3_sub(positions[i + 1], positions[i]);
            // normal perpendicular in XY plane
            let normal = vec3_norm([dir[1], -dir[0], dir[2]]);
            let phase = if s0 < iso_value { 0 } else { 1 };
            points.push(InterfacePoint::new(pos, normal, phase));
        }
    }
    points
}

/// Generate [`LineSegment`]s that display phase interface normals as short ticks.
pub fn interface_normal_lines(
    interface: &[InterfacePoint],
    tick_length: f64,
    color_a: Rgba,
    color_b: Rgba,
) -> Vec<LineSegment> {
    interface
        .iter()
        .map(|ip| {
            let tip = vec3_add(ip.position, vec3_scale(ip.normal, tick_length));
            let color = if ip.phase == 0 { color_a } else { color_b };
            LineSegment::new(ip.position, tip, color)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 4. FSI visualization (fluid + structural deformation)
// ---------------------------------------------------------------------------

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

/// Generate velocity arrow lines for the fluid domain.
pub fn fsi_fluid_arrows(
    snapshot: &FsiSnapshot,
    scale: f64,
    min_mag: f64,
    color: Rgba,
) -> Vec<LineSegment> {
    snapshot
        .fluid_positions
        .iter()
        .zip(snapshot.fluid_velocities.iter())
        .filter(|(_, v)| vec3_len(**v) >= min_mag)
        .map(|(p, v)| LineSegment::new(*p, vec3_add(*p, vec3_scale(*v, scale)), color))
        .collect()
}

/// Generate displacement arrow lines for the structural domain.
pub fn fsi_structural_arrows(
    snapshot: &FsiSnapshot,
    scale: f64,
    min_mag: f64,
    color: Rgba,
) -> Vec<LineSegment> {
    snapshot
        .structural_positions
        .iter()
        .zip(snapshot.structural_displacements.iter())
        .filter(|(_, d)| vec3_len(**d) >= min_mag)
        .map(|(p, d)| LineSegment::new(*p, vec3_add(*p, vec3_scale(*d, scale)), color))
        .collect()
}

/// Generate coupling lines between FSI interface nodes and their structural counterparts.
///
/// Uses the `fsi_interface_indices` in the snapshot to pair fluid boundary
/// nodes with structural nodes (by matching index into `structural_positions`).
pub fn fsi_coupling_lines(snapshot: &FsiSnapshot, color: Rgba) -> Vec<LineSegment> {
    snapshot
        .fsi_interface_indices
        .iter()
        .filter_map(|&fi| {
            let fp = snapshot.fluid_positions.get(fi)?;
            let sp = snapshot.structural_positions.get(fi)?;
            Some(LineSegment::new(*fp, *sp, color))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 5. Electromagnetic field lines
// ---------------------------------------------------------------------------

/// Integration method for field-line tracing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FieldLineMethod {
    /// First-order Euler integration.
    Euler,
    /// Fourth-order Runge-Kutta integration.
    Rk4,
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

impl Default for EmFieldLineConfig {
    fn default() -> Self {
        Self {
            step_size: 0.05,
            max_steps: 500,
            min_magnitude: 1e-10,
            method: FieldLineMethod::Rk4,
        }
    }
}

/// A single traced electromagnetic field line.
#[derive(Debug, Clone)]
pub struct EmFieldLine {
    /// Sequence of positions along the field line.
    pub points: Vec<[f64; 3]>,
    /// Field magnitude sampled at each point.
    pub magnitudes: Vec<f64>,
}

impl EmFieldLine {
    fn new() -> Self {
        Self {
            points: Vec::new(),
            magnitudes: Vec::new(),
        }
    }
}

/// Trace one electromagnetic field line starting from `seed`.
///
/// `field_fn` takes a position and returns the field vector at that point.
pub fn trace_em_field_line<F>(
    seed: [f64; 3],
    field_fn: &F,
    config: &EmFieldLineConfig,
) -> EmFieldLine
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    let mut line = EmFieldLine::new();
    let mut pos = seed;

    // Push seed point
    line.points.push(pos);
    line.magnitudes.push(vec3_len(field_fn(pos)));

    for _ in 0..config.max_steps {
        let b = field_fn(pos);
        let mag = vec3_len(b);
        if mag < config.min_magnitude {
            break;
        }

        let dir = vec3_norm(b);
        pos = match config.method {
            FieldLineMethod::Euler => vec3_add(pos, vec3_scale(dir, config.step_size)),
            FieldLineMethod::Rk4 => {
                let h = config.step_size;
                let k1 = vec3_norm(field_fn(pos));
                let k2 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k1, h * 0.5))));
                let k3 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k2, h * 0.5))));
                let k4 = vec3_norm(field_fn(vec3_add(pos, vec3_scale(k3, h))));
                let avg = [
                    (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
                    (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
                    (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
                ];
                vec3_add(pos, vec3_scale(avg, h))
            }
        };
        let b_new = field_fn(pos);
        line.points.push(pos);
        line.magnitudes.push(vec3_len(b_new));
    }
    line
}

/// Trace multiple field lines from a list of seed points.
pub fn trace_em_field_lines<F>(
    seeds: &[[f64; 3]],
    field_fn: &F,
    config: &EmFieldLineConfig,
) -> Vec<EmFieldLine>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    seeds
        .iter()
        .map(|&s| trace_em_field_line(s, field_fn, config))
        .collect()
}

/// Convert traced field lines to [`LineSegment`]s, colored by local magnitude.
///
/// `mag_min` and `mag_max` define the normalization range for coloring;
/// the cool colormap (cyan–magenta) is used.
pub fn em_field_lines_to_segments(
    lines: &[EmFieldLine],
    mag_min: f64,
    mag_max: f64,
) -> Vec<LineSegment> {
    let mut segs = Vec::new();
    for line in lines {
        let pts = &line.points;
        let mags = &line.magnitudes;
        for i in 0..pts.len().saturating_sub(1) {
            let t = normalize_scalar(mags[i], mag_min, mag_max);
            let color = cool_colormap(t);
            segs.push(LineSegment::new(pts[i], pts[i + 1], color));
        }
    }
    segs
}

// ---------------------------------------------------------------------------
// 6. Temperature-dependent color transitions
// ---------------------------------------------------------------------------

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

/// Build a standard iron-body black-body temperature ramp
/// (black → red → orange → yellow → white).
pub fn blackbody_ramp() -> TempColorRamp {
    TempColorRamp::sorted(vec![
        TempColorKey::new(0.0, Rgba::new(0.0, 0.0, 0.0, 1.0)),
        TempColorKey::new(300.0, Rgba::new(0.1, 0.0, 0.0, 1.0)),
        TempColorKey::new(800.0, Rgba::new(0.8, 0.1, 0.0, 1.0)),
        TempColorKey::new(1200.0, Rgba::new(1.0, 0.5, 0.0, 1.0)),
        TempColorKey::new(2000.0, Rgba::new(1.0, 1.0, 0.5, 1.0)),
        TempColorKey::new(5000.0, Rgba::new(1.0, 1.0, 1.0, 1.0)),
    ])
}

// ---------------------------------------------------------------------------
// 7. Multiphysics convergence dashboard
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// 8. Field interpolation between different meshes
// ---------------------------------------------------------------------------

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

/// Transfer a scalar field from a source mesh to a target mesh.
///
/// `source_positions` and `source_values` define the source field.
/// `target_positions` defines the query points.
/// Returns one interpolated scalar per target node.
pub fn interpolate_scalar_field(
    source_positions: &[[f64; 3]],
    source_values: &[f64],
    target_positions: &[[f64; 3]],
    method: MeshInterpolationMethod,
) -> Vec<f64> {
    let n_src = source_positions.len().min(source_values.len());
    target_positions
        .iter()
        .map(|&qp| {
            match method {
                MeshInterpolationMethod::NearestNeighbor => {
                    let mut best_dist = f64::INFINITY;
                    let mut best_val = 0.0;
                    for i in 0..n_src {
                        let d = vec3_len(vec3_sub(source_positions[i], qp));
                        if d < best_dist {
                            best_dist = d;
                            best_val = source_values[i];
                        }
                    }
                    best_val
                }
                MeshInterpolationMethod::InverseDistance { power, k_neighbors } => {
                    // Collect distances
                    let mut dist_val: Vec<(f64, f64)> = (0..n_src)
                        .map(|i| {
                            (
                                vec3_len(vec3_sub(source_positions[i], qp)),
                                source_values[i],
                            )
                        })
                        .collect();
                    dist_val
                        .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                    let k = k_neighbors.min(n_src).max(1);
                    let mut w_sum = 0.0;
                    let mut v_sum = 0.0;
                    for (d, v) in dist_val.iter().take(k) {
                        if *d < 1e-14 {
                            return *v;
                        }
                        let w = 1.0 / d.powf(power);
                        w_sum += w;
                        v_sum += w * v;
                    }
                    if w_sum < 1e-300 { 0.0 } else { v_sum / w_sum }
                }
                MeshInterpolationMethod::Barycentric => {
                    // Fallback to IDW-2 with k=3 for now
                    let mut dist_val: Vec<(f64, f64)> = (0..n_src)
                        .map(|i| {
                            (
                                vec3_len(vec3_sub(source_positions[i], qp)),
                                source_values[i],
                            )
                        })
                        .collect();
                    dist_val
                        .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                    let k = 3_usize.min(n_src);
                    let mut w_sum = 0.0;
                    let mut v_sum = 0.0;
                    for (d, v) in dist_val.iter().take(k) {
                        if *d < 1e-14 {
                            return *v;
                        }
                        let w = 1.0 / (d * d);
                        w_sum += w;
                        v_sum += w * v;
                    }
                    if w_sum < 1e-300 { 0.0 } else { v_sum / w_sum }
                }
            }
        })
        .collect()
}

/// Transfer a vector field from a source mesh to a target mesh using IDW-2.
pub fn interpolate_vector_field(
    source_positions: &[[f64; 3]],
    source_vectors: &[[f64; 3]],
    target_positions: &[[f64; 3]],
    k_neighbors: usize,
) -> Vec<[f64; 3]> {
    let n_src = source_positions.len().min(source_vectors.len());
    target_positions
        .iter()
        .map(|&qp| {
            let mut dist_vec: Vec<(f64, [f64; 3])> = (0..n_src)
                .map(|i| {
                    (
                        vec3_len(vec3_sub(source_positions[i], qp)),
                        source_vectors[i],
                    )
                })
                .collect();
            dist_vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let k = k_neighbors.min(n_src).max(1);
            let mut w_sum = 0.0;
            let mut v_sum = [0.0_f64; 3];
            for (d, v) in dist_vec.iter().take(k) {
                if *d < 1e-14 {
                    return *v;
                }
                let w = 1.0 / (d * d);
                w_sum += w;
                v_sum = vec3_add(v_sum, vec3_scale(*v, w));
            }
            if w_sum < 1e-300 {
                [0.0; 3]
            } else {
                vec3_scale(v_sum, 1.0 / w_sum)
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 9. Clipping planes for multiphysics
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// 10. Timeline synchronization across physics domains
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Multiphysics viz utility functions
// ---------------------------------------------------------------------------

/// Generate a color-coded scalar field visualization as colored [`LineSegment`]s
/// (vertical bars from node to colored tip), useful for 2-D scalar profiles.
pub fn scalar_profile_bars(
    positions: &[[f64; 3]],
    scalars: &[f64],
    min_val: f64,
    max_val: f64,
    bar_height_scale: f64,
) -> Vec<LineSegment> {
    positions
        .iter()
        .zip(scalars.iter())
        .map(|(&pos, &val)| {
            let t = normalize_scalar(val, min_val, max_val);
            let color = hot_colormap(t);
            let tip = [pos[0], pos[1] + t * bar_height_scale, pos[2]];
            LineSegment::new(pos, tip, color)
        })
        .collect()
}

/// Generate iso-contour lines in 2-D (XZ plane, y=0) using simple marching squares
/// on a regular grid.
///
/// `grid_values` is row-major with `nx` columns and `nz` rows.
/// `cell_size` is the spacing between grid points.
/// `origin` is the bottom-left corner of the grid.
pub fn iso_contour_lines(
    grid_values: &[f64],
    nx: usize,
    nz: usize,
    iso_value: f64,
    cell_size: f64,
    origin: [f64; 3],
    color: Rgba,
) -> Vec<LineSegment> {
    let mut segs = Vec::new();
    if nx < 2 || nz < 2 {
        return segs;
    }

    let idx = |ix: usize, iz: usize| iz * nx + ix;

    for iz in 0..nz - 1 {
        for ix in 0..nx - 1 {
            let s00 = *grid_values.get(idx(ix, iz)).unwrap_or(&0.0);
            let s10 = *grid_values.get(idx(ix + 1, iz)).unwrap_or(&0.0);
            let s01 = *grid_values.get(idx(ix, iz + 1)).unwrap_or(&0.0);
            let s11 = *grid_values.get(idx(ix + 1, iz + 1)).unwrap_or(&0.0);

            let above = [
                s00 >= iso_value,
                s10 >= iso_value,
                s11 >= iso_value,
                s01 >= iso_value,
            ];
            let code: u8 = (above[0] as u8)
                | ((above[1] as u8) << 1)
                | ((above[2] as u8) << 2)
                | ((above[3] as u8) << 3);

            if code == 0 || code == 15 {
                continue;
            }

            let x0 = origin[0] + ix as f64 * cell_size;
            let z0 = origin[2] + iz as f64 * cell_size;
            let x1 = x0 + cell_size;
            let z1 = z0 + cell_size;
            let y = origin[1];

            // Helper: interpolate along an edge
            let interp = |a: f64, b: f64, pa: [f64; 3], pb: [f64; 3]| -> [f64; 3] {
                let denom = b - a;
                if denom.abs() < 1e-14 {
                    return pa;
                }
                let t = (iso_value - a) / denom;
                vec3_lerp(pa, pb, t)
            };

            // Edge midpoints
            let e_bottom = interp(s00, s10, [x0, y, z0], [x1, y, z0]);
            let e_right = interp(s10, s11, [x1, y, z0], [x1, y, z1]);
            let e_top = interp(s11, s01, [x1, y, z1], [x0, y, z1]);
            let e_left = interp(s01, s00, [x0, y, z1], [x0, y, z0]);

            // Marching squares table (simplified — connect crossing edges)
            let push_seg = |segs: &mut Vec<LineSegment>, a: [f64; 3], b: [f64; 3]| {
                segs.push(LineSegment::new(a, b, color));
            };

            match code {
                1 | 14 => push_seg(&mut segs, e_bottom, e_left),
                2 | 13 => push_seg(&mut segs, e_bottom, e_right),
                3 | 12 => push_seg(&mut segs, e_left, e_right),
                4 | 11 => push_seg(&mut segs, e_right, e_top),
                5 => {
                    push_seg(&mut segs, e_bottom, e_right);
                    push_seg(&mut segs, e_left, e_top);
                }
                6 | 9 => push_seg(&mut segs, e_bottom, e_top),
                7 | 8 => push_seg(&mut segs, e_left, e_top),
                10 => {
                    push_seg(&mut segs, e_bottom, e_left);
                    push_seg(&mut segs, e_right, e_top);
                }
                _ => {}
            }
        }
    }
    segs
}

/// Compute the multiphysics residual norm — the maximum absolute residual across
/// multiple field arrays.
///
/// Used to assess overall coupling convergence.
pub fn multiphysics_residual_norm(fields: &[&[f64]]) -> f64 {
    fields
        .iter()
        .flat_map(|f| f.iter())
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max)
}

/// Compute per-node coupling error between two sets of field values.
///
/// Returns the absolute pointwise difference.
pub fn coupling_error(old_values: &[f64], new_values: &[f64]) -> Vec<f64> {
    old_values
        .iter()
        .zip(new_values.iter())
        .map(|(a, b)| (b - a).abs())
        .collect()
}

/// Smooth a 1-D scalar field using a simple box filter of half-width `radius`.
pub fn smooth_scalar_field(values: &[f64], radius: usize) -> Vec<f64> {
    let n = values.len();
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(radius);
            let hi = (i + radius + 1).min(n);
            let sum: f64 = values[lo..hi].iter().sum();
            sum / (hi - lo) as f64
        })
        .collect()
}

/// Apply a clipping plane to a scalar field, returning only the visible nodes
/// and their values.
pub fn clip_scalar_field<'a>(
    nodes: &'a [[f64; 3]],
    values: &'a [f64],
    plane: &ClippingPlane,
) -> (Vec<&'a [f64; 3]>, Vec<f64>) {
    let mut out_nodes = Vec::new();
    let mut out_vals = Vec::new();
    for (p, v) in nodes.iter().zip(values.iter()) {
        if plane.visible(*p) {
            out_nodes.push(p);
            out_vals.push(*v);
        }
    }
    (out_nodes, out_vals)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- vec3 helpers ----

    #[test]
    fn test_vec3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = vec3_add(a, b);
        assert!((c[0] - 5.0).abs() < 1e-12);
        assert!((c[1] - 7.0).abs() < 1e-12);
        assert!((c[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_dot() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(vec3_dot(a, b).abs() < 1e-12);
        assert!((vec3_dot(a, a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_norm_unit() {
        let v = [3.0, 4.0, 0.0];
        let u = vec3_norm(v);
        assert!((vec3_len(u) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_norm_zero() {
        let v = [0.0, 0.0, 0.0];
        let u = vec3_norm(v);
        assert_eq!(u, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_vec3_cross() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }

    // ---- Rgba ----

    #[test]
    fn test_rgba_lerp_midpoint() {
        let r = Rgba::red();
        let b = Rgba::blue();
        let m = Rgba::lerp(r, b, 0.5);
        assert!((m.r - 0.5).abs() < 1e-6);
        assert!((m.b - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_rgba_lerp_endpoints() {
        let r = Rgba::red();
        let b = Rgba::blue();
        let a = Rgba::lerp(r, b, 0.0);
        let c = Rgba::lerp(r, b, 1.0);
        assert!((a.r - 1.0).abs() < 1e-6);
        assert!((c.b - 1.0).abs() < 1e-6);
    }

    // ---- thermal-mechanical overlay ----

    #[test]
    fn test_thermal_mechanical_overlay_length() {
        let config = ThermalMechanicalOverlayConfig::default();
        let t_nodes: Vec<FieldNode> = (0..5)
            .map(|i| FieldNode::new([i as f64, 0.0, 0.0], 300.0 + i as f64 * 100.0, [0.0; 3]))
            .collect();
        let d_nodes: Vec<FieldNode> = (0..5)
            .map(|i| FieldNode::new([i as f64, 0.0, 0.0], 0.0, [0.001 * i as f64, 0.0, 0.0]))
            .collect();
        let colors = thermal_mechanical_overlay(&t_nodes, &d_nodes, &config);
        assert_eq!(colors.len(), 5);
    }

    #[test]
    fn test_thermal_mechanical_overlay_pure_thermal() {
        let config = ThermalMechanicalOverlayConfig {
            thermal_weight: 1.0,
            ..Default::default()
        };
        let t_nodes = vec![FieldNode::new([0.0; 3], 1000.0, [0.0; 3])]; // max temp
        let d_nodes = vec![FieldNode::new([0.0; 3], 0.0, [0.0; 3])];
        let colors = thermal_mechanical_overlay(&t_nodes, &d_nodes, &config);
        // hot colormap at t=1 → white
        assert!(colors[0].r > 0.9);
    }

    // ---- coupled field arrows ----

    #[test]
    fn test_coupled_arrows_count() {
        let config = CoupledArrowConfig::default();
        let nodes = vec![
            FieldNode::new([0.0, 0.0, 0.0], 0.0, [0.0; 3]),
            FieldNode::new([1.0, 0.0, 0.0], 0.0, [0.0; 3]),
        ];
        let hf = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let disp = vec![[0.0, 1.0, 0.0], [0.0, 0.5, 0.0]];
        let segs = coupled_field_arrows(&nodes, &hf, &disp, &config);
        // 2 nodes × 2 arrows each = 4
        assert_eq!(segs.len(), 4);
    }

    #[test]
    fn test_coupled_arrows_below_threshold_filtered() {
        let config = CoupledArrowConfig {
            min_magnitude: 1e3,
            ..Default::default()
        };
        let nodes = vec![FieldNode::new([0.0; 3], 0.0, [0.0; 3])];
        let hf = vec![[1.0, 0.0, 0.0]];
        let disp = vec![[0.0, 1.0, 0.0]];
        let segs = coupled_field_arrows(&nodes, &hf, &disp, &config);
        assert_eq!(segs.len(), 0);
    }

    // ---- interface extraction ----

    #[test]
    fn test_extract_interface_single_crossing() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let scalars: Vec<f64> = (0..10).map(|i| i as f64 - 4.5).collect();
        let pts = extract_interface_points(&positions, &scalars, 0.0);
        assert_eq!(pts.len(), 1);
    }

    #[test]
    fn test_extract_interface_no_crossing() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let scalars = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let pts = extract_interface_points(&positions, &scalars, 0.0);
        assert_eq!(pts.len(), 0);
    }

    #[test]
    fn test_interface_normal_lines() {
        let pts = vec![
            InterfacePoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0),
            InterfacePoint::new([1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 1),
        ];
        let segs = interface_normal_lines(&pts, 0.1, Rgba::red(), Rgba::blue());
        assert_eq!(segs.len(), 2);
        // First segment starts at origin
        assert!((segs[0].start[0]).abs() < 1e-10);
    }

    // ---- FSI ----

    #[test]
    fn test_fsi_fluid_arrows_count() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.fluid_positions = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        snap.fluid_velocities = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let segs = fsi_fluid_arrows(&snap, 1.0, 1e-10, Rgba::cyan());
        assert_eq!(segs.len(), 2);
    }

    #[test]
    fn test_fsi_coupling_lines_count() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.fluid_positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        snap.structural_positions = vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        snap.fsi_interface_indices = vec![0, 1];
        let segs = fsi_coupling_lines(&snap, Rgba::yellow());
        assert_eq!(segs.len(), 2);
    }

    #[test]
    fn test_fsi_structural_arrows_filtered() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.structural_positions = vec![[0.0; 3]];
        snap.structural_displacements = vec![[0.0; 3]]; // zero displacement
        let segs = fsi_structural_arrows(&snap, 1.0, 1e-10, Rgba::green());
        assert_eq!(segs.len(), 0);
    }

    // ---- EM field lines ----

    #[test]
    fn test_em_field_line_uniform_euler() {
        let config = EmFieldLineConfig {
            step_size: 0.1,
            max_steps: 10,
            min_magnitude: 1e-15,
            method: FieldLineMethod::Euler,
        };
        let field = |_: [f64; 3]| [1.0, 0.0, 0.0];
        let line = trace_em_field_line([0.0, 0.0, 0.0], &field, &config);
        assert_eq!(line.points.len(), 11);
        // Should travel ~1 unit along x
        assert!((line.points.last().unwrap()[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_em_field_line_rk4() {
        let config = EmFieldLineConfig {
            step_size: 0.1,
            max_steps: 5,
            min_magnitude: 1e-15,
            method: FieldLineMethod::Rk4,
        };
        let field = |_: [f64; 3]| [0.0, 1.0, 0.0];
        let line = trace_em_field_line([0.0, 0.0, 0.0], &field, &config);
        assert_eq!(line.points.len(), 6);
        assert!((line.points.last().unwrap()[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_em_field_lines_to_segments() {
        let config = EmFieldLineConfig::default();
        let field = |_: [f64; 3]| [1.0, 0.0, 0.0];
        let seeds = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let lines = trace_em_field_lines(&seeds, &field, &config);
        let segs = em_field_lines_to_segments(&lines, 0.0, 1.0);
        // Each line produces (max_steps) segments
        assert!(!segs.is_empty());
    }

    // ---- temperature ramp ----

    #[test]
    fn test_blackbody_ramp_endpoints() {
        let ramp = blackbody_ramp();
        let cold = ramp.evaluate(0.0);
        let hot = ramp.evaluate(5000.0);
        // cold is black, hot is white
        assert!(cold.r < 0.3);
        assert!(hot.r > 0.9);
        assert!(hot.g > 0.9);
        assert!(hot.b > 0.9);
    }

    #[test]
    fn test_blackbody_ramp_midpoint() {
        let ramp = blackbody_ramp();
        let mid = ramp.evaluate(800.0);
        // at 800 K should be reddish
        assert!(mid.r > mid.b);
    }

    #[test]
    fn test_colorize_temperatures_length() {
        let ramp = blackbody_ramp();
        let temps = vec![300.0, 600.0, 900.0, 1200.0];
        let colors = ramp.colorize_temperatures(&temps);
        assert_eq!(colors.len(), 4);
    }

    // ---- convergence dashboard ----

    #[test]
    fn test_dashboard_all_converged() {
        let mut dash = MultiphysicsConvergenceDashboard::new(100);
        let mut d1 = DomainConvergence::new("thermal", 1e-6);
        d1.push(1e-7);
        let mut d2 = DomainConvergence::new("mechanical", 1e-6);
        d2.push(5e-8);
        dash.add_domain(d1);
        dash.add_domain(d2);
        assert!(dash.all_converged());
    }

    #[test]
    fn test_dashboard_not_converged() {
        let mut dash = MultiphysicsConvergenceDashboard::new(100);
        let mut d = DomainConvergence::new("fluid", 1e-6);
        d.push(1e-3);
        dash.add_domain(d);
        assert!(!dash.all_converged());
    }

    #[test]
    fn test_dashboard_iteration_limit() {
        let mut dash = MultiphysicsConvergenceDashboard::new(10);
        for _ in 0..10 {
            dash.next_iteration();
        }
        assert!(dash.iteration_limit_reached());
    }

    #[test]
    fn test_convergence_rate() {
        let mut d = DomainConvergence::new("x", 1e-8);
        d.push(1.0);
        d.push(0.1);
        let rate = d.rate().unwrap();
        assert!((rate - 0.1).abs() < 1e-12);
    }

    // ---- mesh interpolation ----

    #[test]
    fn test_interpolate_scalar_nearest() {
        let src_pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_val = vec![10.0, 20.0, 30.0];
        let tgt_pos = vec![[0.1, 0.0, 0.0]];
        let result = interpolate_scalar_field(
            &src_pos,
            &src_val,
            &tgt_pos,
            MeshInterpolationMethod::NearestNeighbor,
        );
        assert!((result[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_scalar_idw() {
        let src_pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_val = vec![0.0, 2.0];
        let tgt_pos = vec![[1.0, 0.0, 0.0]]; // equidistant → should give average
        let result = interpolate_scalar_field(
            &src_pos,
            &src_val,
            &tgt_pos,
            MeshInterpolationMethod::InverseDistance {
                power: 2.0,
                k_neighbors: 2,
            },
        );
        assert!((result[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_vector_field() {
        let src_pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_vec = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tgt_pos = vec![[1.0, 0.0, 0.0]];
        let result = interpolate_vector_field(&src_pos, &src_vec, &tgt_pos, 2);
        assert!((result[0][0] - 1.0).abs() < 1e-10);
    }

    // ---- clipping planes ----

    #[test]
    fn test_clipping_plane_visible() {
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(plane.visible([1.0, 0.0, 0.0]));
        assert!(!plane.visible([-1.0, 0.0, 0.0]));
    }

    #[test]
    fn test_clipping_plane_signed_distance() {
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let d = plane.signed_distance([0.0, 3.0, 0.0]);
        assert!((d - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_clipping_stack_filter_nodes() {
        let mut stack = ClippingPlaneStack::new();
        stack.add(ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])); // x >= 0
        let nodes = vec![
            FieldNode::new([1.0, 0.0, 0.0], 0.0, [0.0; 3]),
            FieldNode::new([-1.0, 0.0, 0.0], 0.0, [0.0; 3]),
        ];
        let visible = stack.filter_nodes(&nodes);
        assert_eq!(visible.len(), 1);
        assert!((visible[0].position[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clipping_stack_filter_segments() {
        let mut stack = ClippingPlaneStack::new();
        stack.add(ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
        let segs = vec![
            LineSegment::new([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], Rgba::white()), // visible
            LineSegment::new([-1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], Rgba::white()), // clipped
            LineSegment::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], Rgba::white()), // straddles
        ];
        let visible = stack.filter_segments(&segs);
        assert_eq!(visible.len(), 1);
    }

    // ---- timeline synchronization ----

    #[test]
    fn test_timeline_global_bounds() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.1));
        let mut s1 = PhysicsTimeSeries::new("thermal");
        s1.push(0.0, 0.0);
        s1.push(1.0, 1.0);
        let mut s2 = PhysicsTimeSeries::new("fluid");
        s2.push(0.5, 0.0);
        s2.push(2.0, 2.0);
        tl.add_domain(s1);
        tl.add_domain(s2);
        assert!((tl.global_start() - 0.0).abs() < 1e-12);
        assert!((tl.global_end() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_timeline_fixed_sync_step() {
        let tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.25));
        assert!((tl.sync_step() - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_timeline_sample_all() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.1));
        let mut s = PhysicsTimeSeries::new("domain");
        s.push(0.0, 10.0);
        s.push(1.0, 20.0);
        tl.add_domain(s);
        let vals = tl.sample_all(0.5);
        assert_eq!(vals.len(), 1);
        assert!((vals[0].unwrap() - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_interpolate_out_of_range() {
        let mut s = PhysicsTimeSeries::new("x");
        s.push(1.0, 100.0);
        s.push(2.0, 200.0);
        // Before first point
        assert!((s.interpolate_at(0.0).unwrap() - 100.0).abs() < 1e-10);
        // After last point
        assert!((s.interpolate_at(3.0).unwrap() - 200.0).abs() < 1e-10);
    }

    // ---- scalar profile bars ----

    #[test]
    fn test_scalar_profile_bars_count() {
        let pos: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, 0.0, 0.0]).collect();
        let vals: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let bars = scalar_profile_bars(&pos, &vals, 0.0, 7.0, 1.0);
        assert_eq!(bars.len(), 8);
    }

    #[test]
    fn test_scalar_profile_bar_height() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let vals = vec![1.0]; // max value
        let bars = scalar_profile_bars(&pos, &vals, 0.0, 1.0, 2.0);
        assert_eq!(bars.len(), 1);
        // tip y = 0 + 1.0 * 2.0 = 2.0
        assert!((bars[0].end[1] - 2.0).abs() < 1e-10);
    }

    // ---- iso-contour lines ----

    #[test]
    fn test_iso_contour_2x2_single_crossing() {
        // 2×2 grid with iso at 0.5: all cells have crossing
        let grid = vec![0.0, 0.0, 1.0, 1.0];
        let segs = iso_contour_lines(&grid, 2, 2, 0.5, 1.0, [0.0; 3], Rgba::white());
        assert!(!segs.is_empty());
    }

    #[test]
    fn test_iso_contour_uniform_no_crossing() {
        // Uniform field above iso_value → no contour
        let grid = vec![1.0, 1.0, 1.0, 1.0];
        let segs = iso_contour_lines(&grid, 2, 2, 0.5, 1.0, [0.0; 3], Rgba::white());
        assert_eq!(segs.len(), 0);
    }

    // ---- multiphysics residual norm ----

    #[test]
    fn test_multiphysics_residual_norm() {
        let f1 = vec![0.1, -0.5, 0.3];
        let f2 = vec![1.2, 0.0];
        let norm = multiphysics_residual_norm(&[&f1, &f2]);
        assert!((norm - 1.2).abs() < 1e-12);
    }

    #[test]
    fn test_coupling_error() {
        let old = vec![1.0, 2.0, 3.0];
        let new = vec![1.1, 1.8, 3.5];
        let err = coupling_error(&old, &new);
        assert!((err[0] - 0.1).abs() < 1e-12);
        assert!((err[1] - 0.2).abs() < 1e-12);
        assert!((err[2] - 0.5).abs() < 1e-12);
    }

    // ---- smooth scalar field ----

    #[test]
    fn test_smooth_scalar_field_uniform() {
        let vals = vec![2.0; 10];
        let smoothed = smooth_scalar_field(&vals, 2);
        for v in &smoothed {
            assert!((v - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_smooth_scalar_field_reduces_peak() {
        let mut vals = vec![0.0; 11];
        vals[5] = 10.0; // spike
        let smoothed = smooth_scalar_field(&vals, 2);
        assert!(smoothed[5] < 10.0, "smoothing should reduce peak");
    }

    // ---- clip scalar field ----

    #[test]
    fn test_clip_scalar_field() {
        let nodes = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let vals = vec![10.0, 20.0, 30.0];
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let (vis_nodes, vis_vals) = clip_scalar_field(&nodes, &vals, &plane);
        assert_eq!(vis_nodes.len(), 2);
        assert_eq!(vis_vals.len(), 2);
        // values 10 and 30 should be retained
        assert!(vis_vals.contains(&10.0));
        assert!(vis_vals.contains(&30.0));
    }

    // ---- FieldNode constructor ----

    #[test]
    fn test_field_node_construction() {
        let n = FieldNode::new([1.0, 2.0, 3.0], 42.0, [0.1, 0.2, 0.3]);
        assert!((n.scalar - 42.0).abs() < 1e-12);
        assert!((n.position[0] - 1.0).abs() < 1e-12);
        assert!((n.vector[1] - 0.2).abs() < 1e-12);
    }

    // ---- normalize_scalar ----

    #[test]
    fn test_normalize_scalar_clamped() {
        assert!((normalize_scalar(-1.0, 0.0, 1.0) - 0.0).abs() < 1e-12);
        assert!((normalize_scalar(2.0, 0.0, 1.0) - 1.0).abs() < 1e-12);
        assert!((normalize_scalar(0.5, 0.0, 1.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_scalar_degenerate() {
        // min == max → 0.5
        assert!((normalize_scalar(5.0, 5.0, 5.0) - 0.5).abs() < 1e-12);
    }

    // ---- hot / cool colormap ----

    #[test]
    fn test_hot_colormap_zero() {
        let c = hot_colormap(0.0);
        assert!(c.r < 1e-6 && c.g < 1e-6 && c.b < 1e-6);
    }

    #[test]
    fn test_hot_colormap_one() {
        let c = hot_colormap(1.0);
        assert!((c.r - 1.0).abs() < 1e-6);
        assert!((c.g - 1.0).abs() < 1e-6);
        assert!((c.b - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_cool_colormap_range() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let c = cool_colormap(t);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }

    // ---- FsiSnapshot new ----

    #[test]
    fn test_fsi_snapshot_initial_state() {
        let snap = FsiSnapshot::new(1.5);
        assert!((snap.time - 1.5).abs() < 1e-12);
        assert!(snap.fluid_positions.is_empty());
        assert!(snap.structural_positions.is_empty());
    }

    // ---- MultiphysicsTimeline sync grid ----

    #[test]
    fn test_timeline_sync_grid_length() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.5));
        let mut s = PhysicsTimeSeries::new("d");
        s.push(0.0, 0.0);
        s.push(2.0, 1.0);
        tl.add_domain(s);
        let grid = tl.sync_grid();
        // [0.0, 0.5, 1.0, 1.5, 2.0] → 5 points
        assert_eq!(grid.len(), 5);
    }
}
