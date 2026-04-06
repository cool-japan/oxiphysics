#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cable and rope constraint solvers.
//!
//! This module provides a comprehensive set of tools for simulating cables,
//! ropes, and cable-pulley systems in physics engines:
//!
//! - **Catenary solver**: Computes the shape of a hanging cable under gravity.
//! - **Elastic cable (Hookean)**: Spring-like cable with configurable stiffness.
//! - **Inextensible cable (position-based)**: Rigid cable constraint via PBD.
//! - **Cable-pulley systems**: Multi-segment cable routed through pulleys.
//! - **Cable winding**: Drum/winch winding simulation.
//! - **Slack detection**: Determines if a cable is taut or slack.
//! - **Cable tension propagation**: Propagates tension through multi-segment cables.
//! - **Multi-segment cables**: Discretized cable as chain of segments.
//! - **Cable-rigid body attachment**: Attach cables to rigid body anchor points.

#![allow(dead_code)]

use std::f64::consts::PI;

// ── Vector helpers (no nalgebra) ────────────────────────────────────────────

/// 3D dot product.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 3D cross product.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Length of a 3-vector.
fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3-vector. Returns zero if length < epsilon.
fn vec3_normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vec3_len(v);
    if len < 1e-15 {
        return [0.0, 0.0, 0.0];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

/// Subtract two 3-vectors (a - b).
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-vectors.
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-vector.
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Linear interpolation between two 3-vectors.
fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Distance between two points.
fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_len(vec3_sub(b, a))
}

// ── Catenary Solver ─────────────────────────────────────────────────────────

/// Parameters for a catenary cable.
#[derive(Debug, Clone, Copy)]
pub struct CatenaryParams {
    /// Left anchor point.
    pub anchor_a: [f64; 3],
    /// Right anchor point.
    pub anchor_b: [f64; 3],
    /// Total cable length (must be >= distance between anchors).
    pub cable_length: f64,
    /// Linear weight density (weight per unit length, N/m).
    pub weight_density: f64,
}

/// Result of a catenary computation.
#[derive(Debug, Clone)]
pub struct CatenaryResult {
    /// The catenary parameter `a` (horizontal tension / weight density).
    pub catenary_a: f64,
    /// Horizontal tension component.
    pub horizontal_tension: f64,
    /// Minimum sag (vertical distance below straight line at midpoint).
    pub sag: f64,
    /// Sampled points along the catenary curve.
    pub points: Vec<[f64; 3]>,
}

/// Solve the catenary equation for a cable hanging between two anchor points
/// under uniform gravity.
///
/// The catenary shape is `y = a * cosh((x - x0) / a) + c`, where `a = H / w`,
/// `H` is horizontal tension, `w` is weight per unit length.
///
/// Uses Newton's method to find the catenary parameter `a`.
///
/// # Arguments
/// * `params` - Catenary cable parameters.
/// * `num_samples` - Number of sample points along the curve.
///
/// # Returns
/// `Some(CatenaryResult)` on success, `None` if cable is too short.
pub fn solve_catenary(params: &CatenaryParams, num_samples: usize) -> Option<CatenaryResult> {
    let dx = params.anchor_b[0] - params.anchor_a[0];
    let dy = params.anchor_b[1] - params.anchor_a[1];
    let dz = params.anchor_b[2] - params.anchor_a[2];

    // Project to 2D: horizontal distance and vertical difference
    let horizontal_dist = (dx * dx + dz * dz).sqrt();
    let vertical_diff = dy;
    let span = (horizontal_dist * horizontal_dist + vertical_diff * vertical_diff).sqrt();

    if params.cable_length < span - 1e-10 {
        return None; // Cable too short
    }

    // If cable length ~ span, it's a straight line
    if (params.cable_length - span).abs() < 1e-8 {
        let mut points = Vec::with_capacity(num_samples);
        for i in 0..num_samples {
            let t = i as f64 / (num_samples - 1).max(1) as f64;
            points.push(vec3_lerp(params.anchor_a, params.anchor_b, t));
        }
        return Some(CatenaryResult {
            catenary_a: f64::INFINITY,
            horizontal_tension: params.weight_density * params.cable_length * 0.5,
            sag: 0.0,
            points,
        });
    }

    // Bisection method to find catenary parameter `a`
    // Equation: 2*a*sinh(d/(2*a)) = sqrt(L^2 - v^2)
    // where d = horizontal_dist, L = cable_length, v = vertical_diff
    let l = params.cable_length;
    let d = horizontal_dist.max(1e-10);
    let target = (l * l - vertical_diff * vertical_diff).max(0.0).sqrt();

    // f(a) = 2*a*sinh(d/(2*a)) - target
    // As a -> 0+, 2*a*sinh(d/(2a)) -> infinity (cable has large sag)
    // As a -> infinity, 2*a*sinh(d/(2a)) -> d (straight line)
    // We need f(a) = 0, and target > d (cable longer than span)

    let catenary_f = |a_val: f64| -> f64 { 2.0 * a_val * (d / (2.0 * a_val)).sinh() - target };

    // Bisection: find bracket [a_lo, a_hi] such that f changes sign
    let mut a_lo = 1e-6;
    let mut a_hi = d * 100.0;

    // Ensure bracket is valid
    for _ in 0..50 {
        if catenary_f(a_lo) * catenary_f(a_hi) < 0.0 {
            break;
        }
        a_lo *= 0.5;
        a_hi *= 2.0;
    }

    let mut a = (a_lo + a_hi) / 2.0;
    for _iter in 0..100 {
        a = (a_lo + a_hi) / 2.0;
        let fval = catenary_f(a);
        if fval.abs() < 1e-12 {
            break;
        }
        if fval * catenary_f(a_lo) < 0.0 {
            a_hi = a;
        } else {
            a_lo = a;
        }
        if (a_hi - a_lo) < 1e-14 {
            break;
        }
    }

    let horizontal_tension = a * params.weight_density;

    // Compute the x-offset for the catenary vertex
    let _x0 = d / 2.0 - a * ((l + vertical_diff) / (l - vertical_diff)).ln() / 2.0;

    // Sag at midpoint (approximate)
    let mid_y = a * (d / (2.0 * a)).cosh();
    let end_avg_y = (a * 0.0_f64.cosh() + a * (d / a).cosh()) / 2.0;
    let sag = (end_avg_y - mid_y).abs();

    // Generate sample points (3D interpolation)
    let dir_h = if horizontal_dist > 1e-10 {
        [dx / horizontal_dist, 0.0, dz / horizontal_dist]
    } else {
        [1.0, 0.0, 0.0]
    };

    let mut points = Vec::with_capacity(num_samples);
    let n = num_samples.max(2);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let x_local = t * d;
        // Catenary y = a * cosh((x - d/2) / a) - a * cosh(d/(2a)) + linear_offset
        let y_cat = a * ((x_local - d / 2.0) / a).cosh() - a * (d / (2.0 * a)).cosh();
        let y_offset = params.anchor_a[1] + t * vertical_diff + y_cat;

        let px = params.anchor_a[0] + dir_h[0] * x_local;
        let py = y_offset;
        let pz = params.anchor_a[2] + dir_h[2] * x_local;
        points.push([px, py, pz]);
    }

    Some(CatenaryResult {
        catenary_a: a,
        horizontal_tension,
        sag,
        points,
    })
}

// ── Elastic (Hookean) Cable ─────────────────────────────────────────────────

/// An elastic cable constraint following Hooke's law.
///
/// Force = stiffness * max(0, current_length - rest_length) along cable direction.
/// Supports damping and maximum tension clamping.
#[derive(Debug, Clone, Copy)]
pub struct ElasticCable {
    /// Anchor point A position.
    pub anchor_a: [f64; 3],
    /// Anchor point B position.
    pub anchor_b: [f64; 3],
    /// Rest length of the cable.
    pub rest_length: f64,
    /// Stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N*s/m).
    pub damping: f64,
    /// Maximum allowable tension (N). 0 = unlimited.
    pub max_tension: f64,
}

/// Result of elastic cable force computation.
#[derive(Debug, Clone, Copy)]
pub struct ElasticCableForce {
    /// Force vector on anchor A (pulls A toward B when cable stretched).
    pub force_on_a: [f64; 3],
    /// Force vector on anchor B (pulls B toward A when cable stretched).
    pub force_on_b: [f64; 3],
    /// Scalar tension magnitude.
    pub tension: f64,
    /// Current length of the cable.
    pub current_length: f64,
    /// Cable stretch (current_length - rest_length), clamped to >= 0.
    pub stretch: f64,
}

impl ElasticCable {
    /// Create a new elastic cable.
    pub fn new(anchor_a: [f64; 3], anchor_b: [f64; 3], rest_length: f64, stiffness: f64) -> Self {
        Self {
            anchor_a,
            anchor_b,
            rest_length,
            stiffness,
            damping: 0.0,
            max_tension: 0.0,
        }
    }

    /// Set the damping coefficient.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping;
        self
    }

    /// Set maximum tension.
    pub fn with_max_tension(mut self, max_tension: f64) -> Self {
        self.max_tension = max_tension;
        self
    }

    /// Compute the cable force given current velocities of the anchor points.
    ///
    /// The cable only exerts force when stretched beyond its rest length (tension only,
    /// no compression).
    pub fn compute_force(&self, vel_a: [f64; 3], vel_b: [f64; 3]) -> ElasticCableForce {
        let delta = vec3_sub(self.anchor_b, self.anchor_a);
        let current_length = vec3_len(delta);

        if current_length < 1e-15 {
            return ElasticCableForce {
                force_on_a: [0.0; 3],
                force_on_b: [0.0; 3],
                tension: 0.0,
                current_length: 0.0,
                stretch: 0.0,
            };
        }

        let direction = vec3_normalize(delta);
        let stretch = (current_length - self.rest_length).max(0.0);

        // Spring force
        let mut tension = self.stiffness * stretch;

        // Damping force (only along cable direction)
        let relative_vel = vec3_sub(vel_b, vel_a);
        let vel_along = dot3(relative_vel, direction);
        tension += self.damping * vel_along;

        // Clamp tension (cable can only pull, not push)
        tension = tension.max(0.0);
        if self.max_tension > 0.0 {
            tension = tension.min(self.max_tension);
        }

        let force_on_a = vec3_scale(direction, tension);
        let force_on_b = vec3_scale(direction, -tension);

        ElasticCableForce {
            force_on_a,
            force_on_b,
            tension,
            current_length,
            stretch,
        }
    }
}

// ── Inextensible Cable (Position-Based) ─────────────────────────────────────

/// Position-based dynamics constraint for an inextensible cable.
///
/// Projects particle positions to satisfy the distance constraint.
/// The cable has a fixed maximum length and cannot stretch.
#[derive(Debug, Clone, Copy)]
pub struct InextensibleCable {
    /// Maximum cable length.
    pub max_length: f64,
    /// Compliance (inverse stiffness, 0 = perfectly rigid).
    pub compliance: f64,
    /// Number of solver iterations.
    pub iterations: usize,
}

/// A position correction from the PBD solver.
#[derive(Debug, Clone, Copy)]
pub struct PbdCorrection {
    /// Position correction for particle A.
    pub delta_a: [f64; 3],
    /// Position correction for particle B.
    pub delta_b: [f64; 3],
    /// The constraint error (positive = violation).
    pub constraint_error: f64,
}

impl InextensibleCable {
    /// Create a new inextensible cable constraint.
    pub fn new(max_length: f64) -> Self {
        Self {
            max_length,
            compliance: 0.0,
            iterations: 10,
        }
    }

    /// Set compliance (XPBD).
    pub fn with_compliance(mut self, compliance: f64) -> Self {
        self.compliance = compliance;
        self
    }

    /// Set iteration count.
    pub fn with_iterations(mut self, iterations: usize) -> Self {
        self.iterations = iterations;
        self
    }

    /// Project positions of two particles to satisfy the cable constraint.
    ///
    /// `inv_mass_a` and `inv_mass_b` are the inverse masses of the particles.
    /// Returns the position correction for both particles.
    pub fn project(
        &self,
        pos_a: [f64; 3],
        pos_b: [f64; 3],
        inv_mass_a: f64,
        inv_mass_b: f64,
        dt: f64,
    ) -> PbdCorrection {
        let delta = vec3_sub(pos_b, pos_a);
        let dist = vec3_len(delta);

        // Only constrain if stretched beyond max_length (cable, not rod)
        if dist <= self.max_length || dist < 1e-15 {
            return PbdCorrection {
                delta_a: [0.0; 3],
                delta_b: [0.0; 3],
                constraint_error: 0.0,
            };
        }

        let n = vec3_normalize(delta);
        let c = dist - self.max_length; // constraint error

        let w_sum = inv_mass_a + inv_mass_b;
        let alpha = self.compliance / (dt * dt);

        if (w_sum + alpha).abs() < 1e-15 {
            return PbdCorrection {
                delta_a: [0.0; 3],
                delta_b: [0.0; 3],
                constraint_error: c,
            };
        }

        let lambda = c / (w_sum + alpha);
        let delta_a = vec3_scale(n, inv_mass_a * lambda);
        let delta_b = vec3_scale(n, -inv_mass_b * lambda);

        PbdCorrection {
            delta_a,
            delta_b,
            constraint_error: c,
        }
    }

    /// Iteratively project a chain of particles connected by cable segments.
    ///
    /// `positions` - Mutable slice of particle positions.
    /// `inv_masses` - Inverse masses for each particle.
    /// `segment_lengths` - Maximum length for each cable segment (N-1 segments for N particles).
    /// `dt` - Time step.
    pub fn project_chain(
        &self,
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        segment_lengths: &[f64],
        dt: f64,
    ) {
        let n = positions.len();
        if n < 2 || segment_lengths.len() != n - 1 {
            return;
        }

        for _iter in 0..self.iterations {
            for i in 0..n - 1 {
                let delta = vec3_sub(positions[i + 1], positions[i]);
                let dist = vec3_len(delta);
                let max_len = segment_lengths[i];

                if dist <= max_len || dist < 1e-15 {
                    continue;
                }

                let dir = vec3_normalize(delta);
                let c = dist - max_len;
                let w_sum = inv_masses[i] + inv_masses[i + 1];
                let alpha = self.compliance / (dt * dt);

                if (w_sum + alpha).abs() < 1e-15 {
                    continue;
                }

                let lam = c / (w_sum + alpha);
                positions[i] = vec3_add(positions[i], vec3_scale(dir, inv_masses[i] * lam));
                positions[i + 1] =
                    vec3_add(positions[i + 1], vec3_scale(dir, -inv_masses[i + 1] * lam));
            }
        }
    }
}

// ── Cable-Pulley System ─────────────────────────────────────────────────────

/// A pulley in a cable-pulley system.
#[derive(Debug, Clone, Copy)]
pub struct Pulley {
    /// Center position of the pulley.
    pub center: [f64; 3],
    /// Radius of the pulley.
    pub radius: f64,
    /// Friction coefficient at the pulley (0 = frictionless).
    pub friction: f64,
    /// Pulley axis direction (unit vector).
    pub axis: [f64; 3],
}

/// A cable-pulley system consisting of a cable routed through multiple pulleys.
#[derive(Debug, Clone)]
pub struct CablePulleySystem {
    /// Start anchor point.
    pub start_anchor: [f64; 3],
    /// End anchor point.
    pub end_anchor: [f64; 3],
    /// Pulleys in order from start to end.
    pub pulleys: Vec<Pulley>,
    /// Total cable length.
    pub total_length: f64,
}

/// Result of cable-pulley analysis.
#[derive(Debug, Clone)]
pub struct CablePulleyResult {
    /// Tension at each segment (between pulleys).
    pub segment_tensions: Vec<f64>,
    /// Wrap angle at each pulley (radians).
    pub wrap_angles: Vec<f64>,
    /// Total cable path length through the system.
    pub path_length: f64,
    /// Whether the cable is taut.
    pub is_taut: bool,
}

impl CablePulleySystem {
    /// Create a new cable-pulley system.
    pub fn new(start_anchor: [f64; 3], end_anchor: [f64; 3], total_length: f64) -> Self {
        Self {
            start_anchor,
            end_anchor,
            pulleys: Vec::new(),
            total_length,
        }
    }

    /// Add a pulley to the system.
    pub fn add_pulley(&mut self, pulley: Pulley) {
        self.pulleys.push(pulley);
    }

    /// Compute the straight-line path length through all pulleys.
    ///
    /// This is the sum of distances: start -> pulley1 -> pulley2 -> ... -> end.
    pub fn compute_path_length(&self) -> f64 {
        let mut length = 0.0;
        let mut prev = self.start_anchor;

        for pulley in &self.pulleys {
            length += vec3_dist(prev, pulley.center);
            prev = pulley.center;
        }
        length += vec3_dist(prev, self.end_anchor);

        // Add wrap lengths at each pulley (approximate)
        for (i, pulley) in self.pulleys.iter().enumerate() {
            let wrap = self.compute_wrap_angle(i);
            length += pulley.radius * wrap;
        }

        length
    }

    /// Compute the wrap angle at a given pulley index.
    fn compute_wrap_angle(&self, pulley_idx: usize) -> f64 {
        if pulley_idx >= self.pulleys.len() {
            return 0.0;
        }

        let pulley = &self.pulleys[pulley_idx];
        let prev = if pulley_idx == 0 {
            self.start_anchor
        } else {
            self.pulleys[pulley_idx - 1].center
        };
        let next = if pulley_idx == self.pulleys.len() - 1 {
            self.end_anchor
        } else {
            self.pulleys[pulley_idx + 1].center
        };

        let to_prev = vec3_normalize(vec3_sub(prev, pulley.center));
        let to_next = vec3_normalize(vec3_sub(next, pulley.center));

        let cos_angle = dot3(to_prev, to_next).clamp(-1.0, 1.0);
        PI - cos_angle.acos()
    }

    /// Analyze the cable-pulley system.
    ///
    /// Computes tensions in each segment, wrap angles, and tautness.
    /// `input_tension` is the tension applied at the start anchor.
    pub fn analyze(&self, input_tension: f64) -> CablePulleyResult {
        let n_pulleys = self.pulleys.len();
        let n_segments = n_pulleys + 1;

        let mut wrap_angles = Vec::with_capacity(n_pulleys);
        let mut segment_tensions = Vec::with_capacity(n_segments);

        // Compute wrap angles
        for i in 0..n_pulleys {
            wrap_angles.push(self.compute_wrap_angle(i));
        }

        // Tension propagation using Euler-Eytelwein (capstan) equation:
        // T_out = T_in * exp(mu * theta)
        let mut current_tension = input_tension;
        segment_tensions.push(current_tension);

        for i in 0..n_pulleys {
            let mu = self.pulleys[i].friction;
            let theta = wrap_angles[i];
            // Friction reduces tension as cable passes over pulley
            current_tension *= (-mu * theta).exp();
            segment_tensions.push(current_tension);
        }

        let path_length = self.compute_path_length();
        let is_taut = path_length >= self.total_length - 1e-6;

        CablePulleyResult {
            segment_tensions,
            wrap_angles,
            path_length,
            is_taut,
        }
    }
}

// ── Cable Winding ───────────────────────────────────────────────────────────

/// Represents a cable winding drum/winch.
#[derive(Debug, Clone, Copy)]
pub struct CableWindingDrum {
    /// Drum radius (m).
    pub drum_radius: f64,
    /// Drum width (m).
    pub drum_width: f64,
    /// Cable diameter (m).
    pub cable_diameter: f64,
    /// Current wound length (m).
    pub wound_length: f64,
    /// Maximum cable capacity (m).
    pub max_capacity: f64,
    /// Current winding speed (m/s).
    pub winding_speed: f64,
}

/// State of the winding drum.
#[derive(Debug, Clone, Copy)]
pub struct WindingState {
    /// Number of complete layers wound.
    pub layers: u32,
    /// Effective drum radius with wound cable.
    pub effective_radius: f64,
    /// Remaining capacity (m).
    pub remaining_capacity: f64,
    /// Current tension in cable at drum (N).
    pub drum_tension: f64,
    /// Required torque at drum (N*m).
    pub required_torque: f64,
}

impl CableWindingDrum {
    /// Create a new cable winding drum.
    pub fn new(drum_radius: f64, drum_width: f64, cable_diameter: f64) -> Self {
        // Compute max capacity based on geometry
        let cables_per_layer = (drum_width / cable_diameter).floor() as u32;
        // Approximate: assume ~20 layers max
        let max_layers = 20u32;
        let mut capacity = 0.0;
        for layer in 0..max_layers {
            let r = drum_radius + (layer as f64 + 0.5) * cable_diameter;
            let circumference = 2.0 * PI * r;
            capacity += circumference * cables_per_layer as f64;
        }

        Self {
            drum_radius,
            drum_width,
            cable_diameter,
            wound_length: 0.0,
            max_capacity: capacity,
            winding_speed: 0.0,
        }
    }

    /// Set the current wound length.
    pub fn with_wound_length(mut self, length: f64) -> Self {
        self.wound_length = length;
        self
    }

    /// Compute the current winding state.
    ///
    /// `cable_tension` is the tension in the incoming cable.
    pub fn compute_state(&self, cable_tension: f64) -> WindingState {
        let cables_per_layer = (self.drum_width / self.cable_diameter).floor().max(1.0);
        let circumference_first = 2.0 * PI * (self.drum_radius + self.cable_diameter * 0.5);

        // Determine current layer
        let mut remaining = self.wound_length;
        let mut layer = 0u32;
        let current_radius;

        loop {
            let r = self.drum_radius + (layer as f64 + 0.5) * self.cable_diameter;
            let layer_capacity = 2.0 * PI * r * cables_per_layer;
            if remaining <= layer_capacity || layer >= 100 {
                current_radius = r;
                break;
            }
            remaining -= layer_capacity;
            layer += 1;
        }

        let _circumference = circumference_first; // suppress
        let effective_radius = current_radius;
        let remaining_capacity = (self.max_capacity - self.wound_length).max(0.0);
        let required_torque = cable_tension * effective_radius;

        WindingState {
            layers: layer,
            effective_radius,
            remaining_capacity,
            drum_tension: cable_tension,
            required_torque,
        }
    }

    /// Update the wound length given a time step.
    pub fn update(&mut self, dt: f64) {
        self.wound_length += self.winding_speed * dt;
        self.wound_length = self.wound_length.clamp(0.0, self.max_capacity);
    }
}

// ── Slack Detection ─────────────────────────────────────────────────────────

/// Cable slack state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CableSlackState {
    /// Cable is taut (stretched to or beyond rest length).
    Taut,
    /// Cable is slack (shorter than rest length).
    Slack,
    /// Cable is at rest length (within tolerance).
    Neutral,
}

/// Detect cable slack state.
///
/// # Arguments
/// * `endpoint_a` - Position of anchor A.
/// * `endpoint_b` - Position of anchor B.
/// * `cable_length` - Rest/max length of the cable.
/// * `tolerance` - Tolerance for neutral state detection.
///
/// # Returns
/// The slack state of the cable.
pub fn detect_slack(
    endpoint_a: [f64; 3],
    endpoint_b: [f64; 3],
    cable_length: f64,
    tolerance: f64,
) -> CableSlackState {
    let dist = vec3_dist(endpoint_a, endpoint_b);
    if dist > cable_length + tolerance {
        CableSlackState::Taut
    } else if dist < cable_length - tolerance {
        CableSlackState::Slack
    } else {
        CableSlackState::Neutral
    }
}

/// Compute the slack amount (positive = slack, negative = overstretched).
pub fn compute_slack_amount(endpoint_a: [f64; 3], endpoint_b: [f64; 3], cable_length: f64) -> f64 {
    cable_length - vec3_dist(endpoint_a, endpoint_b)
}

// ── Cable Tension Propagation ───────────────────────────────────────────────

/// A cable segment in a multi-segment cable system.
#[derive(Debug, Clone, Copy)]
pub struct CableSegment {
    /// Start position.
    pub start: [f64; 3],
    /// End position.
    pub end: [f64; 3],
    /// Segment rest length.
    pub rest_length: f64,
    /// Segment stiffness.
    pub stiffness: f64,
    /// Segment damping.
    pub damping: f64,
    /// Linear density (kg/m).
    pub linear_density: f64,
}

/// Result of tension propagation through a multi-segment cable.
#[derive(Debug, Clone)]
pub struct TensionPropagationResult {
    /// Tension at each segment.
    pub tensions: Vec<f64>,
    /// Total cable length.
    pub total_length: f64,
    /// Total cable stretch.
    pub total_stretch: f64,
    /// Whether the entire cable is taut.
    pub all_taut: bool,
}

/// Propagate tension through a chain of cable segments.
///
/// Computes the tension in each segment based on positions and velocities
/// of the nodes, considering stiffness and damping.
///
/// `node_positions` - Positions of N nodes.
/// `node_velocities` - Velocities of N nodes.
/// `segments` - N-1 cable segments connecting consecutive nodes.
pub fn propagate_tension(
    node_positions: &[[f64; 3]],
    node_velocities: &[[f64; 3]],
    segments: &[CableSegment],
) -> TensionPropagationResult {
    let n = node_positions.len();
    if n < 2 || segments.len() != n - 1 {
        return TensionPropagationResult {
            tensions: Vec::new(),
            total_length: 0.0,
            total_stretch: 0.0,
            all_taut: false,
        };
    }

    let mut tensions = Vec::with_capacity(segments.len());
    let mut total_length = 0.0;
    let mut total_stretch = 0.0;
    let mut all_taut = true;

    for (i, seg) in segments.iter().enumerate() {
        let delta = vec3_sub(node_positions[i + 1], node_positions[i]);
        let dist = vec3_len(delta);
        total_length += dist;

        let stretch = (dist - seg.rest_length).max(0.0);
        total_stretch += stretch;

        if stretch < 1e-10 {
            tensions.push(0.0);
            all_taut = false;
            continue;
        }

        let dir = vec3_normalize(delta);
        let rel_vel = vec3_sub(node_velocities[i + 1], node_velocities[i]);
        let vel_along = dot3(rel_vel, dir);

        let tension = (seg.stiffness * stretch + seg.damping * vel_along).max(0.0);
        tensions.push(tension);
    }

    TensionPropagationResult {
        tensions,
        total_length,
        total_stretch,
        all_taut,
    }
}

// ── Multi-Segment Cable ─────────────────────────────────────────────────────

/// A discretized cable represented as a chain of mass points connected by
/// spring-damper segments.
#[derive(Debug, Clone)]
pub struct MultiSegmentCable {
    /// Node positions.
    pub positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Inverse masses (0 = fixed).
    pub inv_masses: Vec<f64>,
    /// Segment rest lengths.
    pub rest_lengths: Vec<f64>,
    /// Global stiffness (N/m).
    pub stiffness: f64,
    /// Global damping coefficient.
    pub damping: f64,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

impl MultiSegmentCable {
    /// Create a multi-segment cable between two endpoints.
    ///
    /// `n_segments` - Number of segments (n_segments+1 nodes).
    /// `mass_per_meter` - Linear mass density.
    /// `stiffness` - Spring stiffness.
    /// `damping` - Damping coefficient.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: [f64; 3],
        end: [f64; 3],
        n_segments: usize,
        mass_per_meter: f64,
        stiffness: f64,
        damping: f64,
        gravity: [f64; 3],
    ) -> Self {
        let n_nodes = n_segments + 1;
        let total_len = vec3_dist(start, end);
        let seg_len = total_len / n_segments as f64;
        let seg_mass = mass_per_meter * seg_len;

        let mut positions = Vec::with_capacity(n_nodes);
        let mut velocities = Vec::with_capacity(n_nodes);
        let mut inv_masses = Vec::with_capacity(n_nodes);
        let mut rest_lengths = Vec::with_capacity(n_segments);

        for i in 0..n_nodes {
            let t = i as f64 / n_segments as f64;
            positions.push(vec3_lerp(start, end, t));
            velocities.push([0.0; 3]);

            // First and last nodes are fixed by default
            if i == 0 || i == n_nodes - 1 {
                inv_masses.push(0.0);
            } else {
                inv_masses.push(1.0 / seg_mass);
            }
        }

        for _i in 0..n_segments {
            rest_lengths.push(seg_len);
        }

        Self {
            positions,
            velocities,
            inv_masses,
            rest_lengths,
            stiffness,
            damping,
            gravity,
        }
    }

    /// Step the cable simulation forward by `dt`.
    ///
    /// Uses symplectic Euler integration with spring-damper forces.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];

        // Gravity
        for i in 0..n {
            if self.inv_masses[i] > 0.0 {
                let mass = 1.0 / self.inv_masses[i];
                forces[i] = vec3_scale(self.gravity, mass);
            }
        }

        // Spring-damper forces
        for i in 0..n - 1 {
            let delta = vec3_sub(self.positions[i + 1], self.positions[i]);
            let dist = vec3_len(delta);
            if dist < 1e-15 {
                continue;
            }

            let dir = vec3_normalize(delta);
            let stretch = dist - self.rest_lengths[i];
            // Cable only pulls (tension), not pushes
            let spring_force = if stretch > 0.0 {
                self.stiffness * stretch
            } else {
                0.0
            };

            let rel_vel = vec3_sub(self.velocities[i + 1], self.velocities[i]);
            let vel_along = dot3(rel_vel, dir);
            let damp_force = self.damping * vel_along;

            let total_force = spring_force + damp_force;
            let f = vec3_scale(dir, total_force);

            forces[i] = vec3_add(forces[i], f);
            forces[i + 1] = vec3_sub(forces[i + 1], f);
        }

        // Symplectic Euler integration
        for i in 0..n {
            if self.inv_masses[i] <= 0.0 {
                continue;
            }
            let acc = vec3_scale(forces[i], self.inv_masses[i]);
            self.velocities[i] = vec3_add(self.velocities[i], vec3_scale(acc, dt));
            self.positions[i] = vec3_add(self.positions[i], vec3_scale(self.velocities[i], dt));
        }
    }

    /// Get the total current length of the cable.
    pub fn current_length(&self) -> f64 {
        let mut len = 0.0;
        for i in 0..self.positions.len() - 1 {
            len += vec3_dist(self.positions[i], self.positions[i + 1]);
        }
        len
    }

    /// Get the total rest length.
    pub fn rest_length(&self) -> f64 {
        self.rest_lengths.iter().sum()
    }

    /// Compute tensions in each segment.
    pub fn segment_tensions(&self) -> Vec<f64> {
        let n = self.positions.len();
        let mut tensions = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let delta = vec3_sub(self.positions[i + 1], self.positions[i]);
            let dist = vec3_len(delta);
            let stretch = (dist - self.rest_lengths[i]).max(0.0);
            tensions.push(self.stiffness * stretch);
        }
        tensions
    }

    /// Get the number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.positions.len()
    }

    /// Get the number of segments.
    pub fn num_segments(&self) -> usize {
        self.rest_lengths.len()
    }
}

// ── Cable-Rigid Body Attachment ─────────────────────────────────────────────

/// Describes how a cable attaches to a rigid body.
#[derive(Debug, Clone, Copy)]
pub struct CableAttachment {
    /// Body-local attachment point (offset from center of mass).
    pub local_offset: [f64; 3],
    /// Body position (center of mass, world space).
    pub body_position: [f64; 3],
    /// Body orientation as a 3x3 rotation matrix (row-major).
    pub body_rotation: [[f64; 3]; 3],
    /// Body inverse mass (0 = static).
    pub inv_mass: f64,
}

impl CableAttachment {
    /// Create a new cable attachment.
    pub fn new(local_offset: [f64; 3], body_position: [f64; 3], inv_mass: f64) -> Self {
        Self {
            local_offset,
            body_position,
            body_rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            inv_mass,
        }
    }

    /// Set the body rotation matrix.
    pub fn with_rotation(mut self, rotation: [[f64; 3]; 3]) -> Self {
        self.body_rotation = rotation;
        self
    }

    /// Compute the world-space attachment point.
    pub fn world_position(&self) -> [f64; 3] {
        let rotated = mat3_vec(self.body_rotation, self.local_offset);
        vec3_add(self.body_position, rotated)
    }
}

/// Multiply a 3x3 matrix by a 3-vector.
fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Compute the cable force on a rigid body at an attachment point.
///
/// Returns `(force, torque)` in world space.
pub fn cable_force_on_body(
    attachment: &CableAttachment,
    cable_endpoint: [f64; 3],
    cable_tension: f64,
) -> ([f64; 3], [f64; 3]) {
    let attach_world = attachment.world_position();
    let delta = vec3_sub(cable_endpoint, attach_world);
    let dist = vec3_len(delta);

    if dist < 1e-15 || cable_tension <= 0.0 {
        return ([0.0; 3], [0.0; 3]);
    }

    let direction = vec3_normalize(delta);
    let force = vec3_scale(direction, cable_tension);

    // Torque = r x F, where r is the lever arm from body center to attach point
    let lever = vec3_sub(attach_world, attachment.body_position);
    let torque = cross3(lever, force);

    (force, torque)
}

// ── Cable Energy Computations ───────────────────────────────────────────────

/// Compute the elastic potential energy in a cable segment.
///
/// PE = 0.5 * k * stretch^2 (only when stretched beyond rest length).
pub fn elastic_potential_energy(current_length: f64, rest_length: f64, stiffness: f64) -> f64 {
    let stretch = (current_length - rest_length).max(0.0);
    0.5 * stiffness * stretch * stretch
}

/// Compute the gravitational potential energy of a cable.
///
/// Assumes uniform linear density.
///
/// `positions` - Node positions of the cable.
/// `linear_density` - Mass per unit length (kg/m).
/// `gravity_magnitude` - Magnitude of gravitational acceleration (m/s^2).
/// `reference_height` - y-coordinate of the reference plane.
pub fn gravitational_potential_energy(
    positions: &[[f64; 3]],
    linear_density: f64,
    gravity_magnitude: f64,
    reference_height: f64,
) -> f64 {
    if positions.len() < 2 {
        return 0.0;
    }

    let mut total_pe = 0.0;
    for i in 0..positions.len() - 1 {
        let seg_len = vec3_dist(positions[i], positions[i + 1]);
        let seg_mass = linear_density * seg_len;
        let mid_height = (positions[i][1] + positions[i + 1][1]) / 2.0;
        total_pe += seg_mass * gravity_magnitude * (mid_height - reference_height);
    }
    total_pe
}

/// Compute the kinetic energy of a cable.
///
/// `velocities` - Node velocities.
/// `masses` - Node masses.
pub fn kinetic_energy(velocities: &[[f64; 3]], masses: &[f64]) -> f64 {
    let n = velocities.len().min(masses.len());
    let mut ke = 0.0;
    for i in 0..n {
        let v2 = dot3(velocities[i], velocities[i]);
        ke += 0.5 * masses[i] * v2;
    }
    ke
}

// ── Cable Natural Frequency ─────────────────────────────────────────────────

/// Compute the fundamental natural frequency of a taut cable (Hz).
///
/// f = (1 / (2 * L)) * sqrt(T / mu)
///
/// `length` - Cable length (m).
/// `tension` - Cable tension (N).
/// `linear_density` - Mass per unit length (kg/m).
pub fn cable_natural_frequency(length: f64, tension: f64, linear_density: f64) -> f64 {
    if length < 1e-15 || linear_density < 1e-15 || tension < 0.0 {
        return 0.0;
    }
    (1.0 / (2.0 * length)) * (tension / linear_density).sqrt()
}

/// Compute the nth harmonic frequency of a cable (Hz).
pub fn cable_harmonic_frequency(
    length: f64,
    tension: f64,
    linear_density: f64,
    harmonic: u32,
) -> f64 {
    cable_natural_frequency(length, tension, linear_density) * harmonic as f64
}

// ── Cable Wave Speed ────────────────────────────────────────────────────────

/// Compute the transverse wave speed in a cable.
///
/// c = sqrt(T / mu)
pub fn cable_wave_speed(tension: f64, linear_density: f64) -> f64 {
    if linear_density < 1e-15 || tension < 0.0 {
        return 0.0;
    }
    (tension / linear_density).sqrt()
}

// ── Cable Sag Computation ───────────────────────────────────────────────────

/// Compute the sag of a cable under uniform load (parabolic approximation).
///
/// sag = w * L^2 / (8 * H)
///
/// `span` - Horizontal span (m).
/// `weight_per_length` - Weight per unit length (N/m).
/// `horizontal_tension` - Horizontal component of cable tension (N).
pub fn parabolic_sag(span: f64, weight_per_length: f64, horizontal_tension: f64) -> f64 {
    if horizontal_tension < 1e-15 {
        return f64::INFINITY;
    }
    weight_per_length * span * span / (8.0 * horizontal_tension)
}

/// Compute the cable length for a parabolic cable profile.
///
/// L ≈ span * (1 + 8/3 * (sag/span)^2) (approximate for small sag/span ratios).
pub fn parabolic_cable_length(span: f64, sag: f64) -> f64 {
    let ratio = sag / span.max(1e-15);
    span * (1.0 + 8.0 / 3.0 * ratio * ratio)
}

// ── Cable Breaking Strength ─────────────────────────────────────────────────

/// Check if a cable would break given the current tension.
///
/// `tension` - Current tension (N).
/// `breaking_strength` - Cable breaking strength (N).
/// `safety_factor` - Required safety factor (typically 2-5).
///
/// Returns true if the cable would break (tension exceeds allowable).
pub fn would_break(tension: f64, breaking_strength: f64, safety_factor: f64) -> bool {
    let allowable = breaking_strength / safety_factor.max(1.0);
    tension > allowable
}

/// Compute the safety margin of a cable.
///
/// Returns the ratio of allowable tension to current tension.
/// Values > 1.0 are safe, < 1.0 means overloaded.
pub fn safety_margin(tension: f64, breaking_strength: f64, safety_factor: f64) -> f64 {
    if tension < 1e-15 {
        return f64::INFINITY;
    }
    let allowable = breaking_strength / safety_factor.max(1.0);
    allowable / tension
}

// ── Cable Drag ──────────────────────────────────────────────────────────────

/// Compute aerodynamic drag force on a cable segment.
///
/// F_drag = 0.5 * rho * Cd * d * L * V^2
///
/// `wind_velocity` - Wind velocity vector (m/s).
/// `cable_direction` - Unit direction of cable segment.
/// `cable_diameter` - Cable diameter (m).
/// `segment_length` - Cable segment length (m).
/// `air_density` - Air density (kg/m^3, typically 1.225).
/// `drag_coefficient` - Drag coefficient (typically 1.0-1.2 for cables).
#[allow(clippy::too_many_arguments)]
pub fn cable_drag_force(
    wind_velocity: [f64; 3],
    cable_direction: [f64; 3],
    cable_diameter: f64,
    segment_length: f64,
    air_density: f64,
    drag_coefficient: f64,
) -> [f64; 3] {
    // Component of wind perpendicular to cable
    let wind_along = dot3(wind_velocity, cable_direction);
    let wind_perp = [
        wind_velocity[0] - wind_along * cable_direction[0],
        wind_velocity[1] - wind_along * cable_direction[1],
        wind_velocity[2] - wind_along * cable_direction[2],
    ];

    let v_perp = vec3_len(wind_perp);
    if v_perp < 1e-15 {
        return [0.0; 3];
    }

    let f_mag =
        0.5 * air_density * drag_coefficient * cable_diameter * segment_length * v_perp * v_perp;
    let wind_perp_dir = vec3_normalize(wind_perp);
    vec3_scale(wind_perp_dir, f_mag)
}

// ── Cable Temperature Effects ───────────────────────────────────────────────

/// Compute cable length change due to thermal expansion.
///
/// delta_L = alpha * L * delta_T
///
/// `original_length` - Original cable length (m).
/// `thermal_coefficient` - Coefficient of thermal expansion (1/K).
/// `temperature_change` - Temperature change (K).
pub fn thermal_length_change(
    original_length: f64,
    thermal_coefficient: f64,
    temperature_change: f64,
) -> f64 {
    thermal_coefficient * original_length * temperature_change
}

/// Compute the adjusted rest length of a cable due to temperature.
pub fn thermal_adjusted_length(
    original_length: f64,
    thermal_coefficient: f64,
    temperature_change: f64,
) -> f64 {
    original_length
        + thermal_length_change(original_length, thermal_coefficient, temperature_change)
}

// ── Cable Creep Model ───────────────────────────────────────────────────────

/// Simple cable creep model.
///
/// Creep strain rate = A * sigma^n
///
/// `stress` - Current cable stress (Pa).
/// `creep_coefficient` - Material creep coefficient A.
/// `creep_exponent` - Stress exponent n (typically 1-5).
pub fn creep_strain_rate(stress: f64, creep_coefficient: f64, creep_exponent: f64) -> f64 {
    if stress < 0.0 {
        return 0.0;
    }
    creep_coefficient * stress.powf(creep_exponent)
}

/// Compute the elongation due to creep over a time interval.
///
/// `original_length` - Original length (m).
/// `stress` - Current stress (Pa).
/// `creep_coefficient` - Material creep coefficient.
/// `creep_exponent` - Stress exponent.
/// `dt` - Time interval (s).
pub fn creep_elongation(
    original_length: f64,
    stress: f64,
    creep_coefficient: f64,
    creep_exponent: f64,
    dt: f64,
) -> f64 {
    let strain_rate = creep_strain_rate(stress, creep_coefficient, creep_exponent);
    strain_rate * original_length * dt
}

// ── Catenary with Concentrated Load ─────────────────────────────────────────

/// Compute the deflection of a cable with a concentrated point load.
///
/// For a cable spanning distance `span` with a point load `P` at distance `a`
/// from the left support, the vertical deflection at the load point is:
///
/// delta = P * a * (span - a) / (span * H)
///
/// where H is the horizontal tension.
///
/// `span` - Horizontal span (m).
/// `load_position` - Distance from left support to load (m).
/// `load` - Point load magnitude (N, downward positive).
/// `horizontal_tension` - Horizontal cable tension (N).
pub fn concentrated_load_deflection(
    span: f64,
    load_position: f64,
    load: f64,
    horizontal_tension: f64,
) -> f64 {
    if horizontal_tension < 1e-15 || span < 1e-15 {
        return 0.0;
    }
    let a = load_position.clamp(0.0, span);
    let b = span - a;
    load * a * b / (span * horizontal_tension)
}

// ── Cable Vibration Damping ─────────────────────────────────────────────────

/// Compute the logarithmic decrement for cable vibration damping.
///
/// delta = 2 * pi * zeta / sqrt(1 - zeta^2)
///
/// `damping_ratio` - Damping ratio (0-1).
pub fn logarithmic_decrement(damping_ratio: f64) -> f64 {
    let zeta = damping_ratio.clamp(0.0, 0.999);
    2.0 * PI * zeta / (1.0 - zeta * zeta).sqrt()
}

/// Compute the number of cycles to decay to a fraction of initial amplitude.
///
/// n = ln(1/fraction) / delta
///
/// `damping_ratio` - Damping ratio.
/// `fraction` - Target fraction of initial amplitude (0-1).
pub fn cycles_to_decay(damping_ratio: f64, fraction: f64) -> f64 {
    let delta = logarithmic_decrement(damping_ratio);
    if delta < 1e-15 || fraction <= 0.0 {
        return f64::INFINITY;
    }
    (1.0 / fraction).ln() / delta
}

// ── Cable Ice Loading ───────────────────────────────────────────────────────

/// Compute the additional weight per unit length due to ice accretion on a cable.
///
/// Assumes cylindrical ice coating around the cable.
///
/// `cable_diameter` - Cable outer diameter (m).
/// `ice_thickness` - Radial ice thickness (m).
/// `ice_density` - Ice density (kg/m^3, typically 900).
/// `gravity` - Gravitational acceleration (m/s^2).
pub fn ice_loading_weight(
    cable_diameter: f64,
    ice_thickness: f64,
    ice_density: f64,
    gravity: f64,
) -> f64 {
    let r_cable = cable_diameter / 2.0;
    let r_ice = r_cable + ice_thickness;
    let ice_area = PI * (r_ice * r_ice - r_cable * r_cable);
    ice_area * ice_density * gravity
}

// ── Cable Galloping ─────────────────────────────────────────────────────────

/// Estimate the onset wind speed for cable galloping (Den Hartog criterion).
///
/// Galloping can occur when dCL/dalpha + CD < 0
///
/// For a simplified model, critical wind speed:
/// Vc = 4 * m * omega * zeta / (rho * D * (dCL_dalpha + CD))
///
/// `mass_per_length` - Cable mass per unit length (kg/m).
/// `natural_freq` - Natural frequency (rad/s).
/// `damping_ratio` - Structural damping ratio.
/// `air_density` - Air density (kg/m^3).
/// `cable_diameter` - Cable diameter (m).
/// `aerodynamic_coeff` - Combined aero coefficient (dCL/dalpha + CD), negative for galloping.
#[allow(clippy::too_many_arguments)]
pub fn galloping_onset_speed(
    mass_per_length: f64,
    natural_freq: f64,
    damping_ratio: f64,
    air_density: f64,
    cable_diameter: f64,
    aerodynamic_coeff: f64,
) -> f64 {
    if aerodynamic_coeff.abs() < 1e-15 || air_density < 1e-15 || cable_diameter < 1e-15 {
        return f64::INFINITY;
    }
    let numerator = 4.0 * mass_per_length * natural_freq * damping_ratio;
    let denominator = air_density * cable_diameter * aerodynamic_coeff.abs();
    (numerator / denominator).abs()
}

// ── Cable Fatigue ───────────────────────────────────────────────────────────

/// Estimate cable fatigue life using a simple S-N curve model.
///
/// N = C / S^m
///
/// `stress_amplitude` - Stress amplitude (Pa).
/// `fatigue_coefficient` - Material constant C.
/// `fatigue_exponent` - Material exponent m (typically 3-5 for steel cables).
pub fn fatigue_life_cycles(
    stress_amplitude: f64,
    fatigue_coefficient: f64,
    fatigue_exponent: f64,
) -> f64 {
    if stress_amplitude < 1e-15 {
        return f64::INFINITY;
    }
    fatigue_coefficient / stress_amplitude.powf(fatigue_exponent)
}

/// Compute cumulative fatigue damage using Miner's rule.
///
/// D = sum(ni / Ni) for each stress level.
///
/// `cycle_counts` - Number of cycles at each stress level.
/// `allowable_cycles` - Allowable cycles at each stress level (from S-N curve).
///
/// Returns damage fraction (>= 1.0 means failure).
pub fn miners_rule_damage(cycle_counts: &[f64], allowable_cycles: &[f64]) -> f64 {
    let n = cycle_counts.len().min(allowable_cycles.len());
    let mut damage = 0.0;
    for i in 0..n {
        if allowable_cycles[i] > 0.0 {
            damage += cycle_counts[i] / allowable_cycles[i];
        }
    }
    damage
}

// ════════════════════════════════════════════════════════════════════════════
// Tests
// ════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catenary_straight_line() {
        let params = CatenaryParams {
            anchor_a: [0.0, 0.0, 0.0],
            anchor_b: [10.0, 0.0, 0.0],
            cable_length: 10.0, // exactly the span
            weight_density: 1.0,
        };
        let result = solve_catenary(&params, 5).unwrap();
        assert_eq!(result.points.len(), 5);
        assert!(result.sag.abs() < 1e-6, "sag={}", result.sag);
    }

    #[test]
    fn test_catenary_too_short() {
        let params = CatenaryParams {
            anchor_a: [0.0, 0.0, 0.0],
            anchor_b: [10.0, 0.0, 0.0],
            cable_length: 5.0, // shorter than span
            weight_density: 1.0,
        };
        assert!(solve_catenary(&params, 10).is_none());
    }

    #[test]
    fn test_catenary_with_sag() {
        let params = CatenaryParams {
            anchor_a: [0.0, 0.0, 0.0],
            anchor_b: [10.0, 0.0, 0.0],
            cable_length: 12.0, // longer than span
            weight_density: 1.0,
        };
        let result = solve_catenary(&params, 20).unwrap();
        assert_eq!(result.points.len(), 20);
        assert!(result.catenary_a > 0.0);
        assert!(result.horizontal_tension > 0.0);
    }

    #[test]
    fn test_elastic_cable_no_stretch() {
        let cable = ElasticCable::new(
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            10.0, // rest length > current length
            100.0,
        );
        let result = cable.compute_force([0.0; 3], [0.0; 3]);
        assert!(result.tension.abs() < 1e-10, "tension={}", result.tension);
        assert!((result.stretch).abs() < 1e-10);
    }

    #[test]
    fn test_elastic_cable_stretched() {
        let cable = ElasticCable::new(
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            5.0, // rest length < current length
            100.0,
        );
        let result = cable.compute_force([0.0; 3], [0.0; 3]);
        // stretch = 10 - 5 = 5, tension = 100 * 5 = 500
        assert!(
            (result.tension - 500.0).abs() < 1e-6,
            "tension={}",
            result.tension
        );
        assert!((result.stretch - 5.0).abs() < 1e-10);
        // Force on A should point toward B (positive x)
        assert!(result.force_on_a[0] > 0.0);
        // Force on B should point toward A (negative x)
        assert!(result.force_on_b[0] < 0.0);
    }

    #[test]
    fn test_elastic_cable_with_damping() {
        let cable =
            ElasticCable::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 5.0, 100.0).with_damping(10.0);
        // B moving away from A (positive damping contribution)
        let result = cable.compute_force([0.0; 3], [1.0, 0.0, 0.0]);
        assert!(
            result.tension > 500.0,
            "tension should be > 500 with damping"
        );
    }

    #[test]
    fn test_elastic_cable_max_tension() {
        let cable = ElasticCable::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 5.0, 100.0)
            .with_max_tension(200.0);
        let result = cable.compute_force([0.0; 3], [0.0; 3]);
        assert!(
            (result.tension - 200.0).abs() < 1e-6,
            "tension={}",
            result.tension
        );
    }

    #[test]
    fn test_inextensible_cable_no_violation() {
        let cable = InextensibleCable::new(10.0);
        let correction = cable.project([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0, 1.0, 0.01);
        assert!(correction.constraint_error.abs() < 1e-10);
        assert!(vec3_len(correction.delta_a) < 1e-10);
    }

    #[test]
    fn test_inextensible_cable_violation() {
        let cable = InextensibleCable::new(5.0);
        let correction = cable.project([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 1.0, 1.0, 0.01);
        assert!(correction.constraint_error > 0.0);
        // A should move toward B
        assert!(correction.delta_a[0] > 0.0);
        // B should move toward A
        assert!(correction.delta_b[0] < 0.0);
    }

    #[test]
    fn test_inextensible_cable_chain() {
        let cable = InextensibleCable::new(5.0).with_iterations(50);
        let mut positions = vec![
            [0.0, 0.0, 0.0],
            [5.5, 0.0, 0.0], // slightly too far from [0]
            [10.0, 0.0, 0.0],
        ];
        let inv_masses = vec![0.0, 1.0, 0.0]; // endpoints fixed
        let segment_lengths = vec![5.0, 5.0];
        cable.project_chain(&mut positions, &inv_masses, &segment_lengths, 0.01);
        // Middle particle should be corrected to satisfy both constraints
        let d0 = vec3_dist(positions[0], positions[1]);
        let d1 = vec3_dist(positions[1], positions[2]);
        assert!(d0 <= 5.01, "d0={d0}");
        assert!(d1 <= 5.01, "d1={d1}");
    }

    #[test]
    fn test_slack_detection_taut() {
        let state = detect_slack([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 8.0, 0.01);
        assert_eq!(state, CableSlackState::Taut);
    }

    #[test]
    fn test_slack_detection_slack() {
        let state = detect_slack([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 10.0, 0.01);
        assert_eq!(state, CableSlackState::Slack);
    }

    #[test]
    fn test_slack_detection_neutral() {
        let state = detect_slack([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 10.0, 0.1);
        assert_eq!(state, CableSlackState::Neutral);
    }

    #[test]
    fn test_slack_amount() {
        let slack = compute_slack_amount([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 10.0);
        assert!((slack - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_tension_propagation() {
        let positions = vec![[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [12.0, 0.0, 0.0]];
        let velocities = vec![[0.0; 3]; 3];
        let segments = vec![
            CableSegment {
                start: positions[0],
                end: positions[1],
                rest_length: 5.0,
                stiffness: 100.0,
                damping: 0.0,
                linear_density: 1.0,
            },
            CableSegment {
                start: positions[1],
                end: positions[2],
                rest_length: 5.0,
                stiffness: 100.0,
                damping: 0.0,
                linear_density: 1.0,
            },
        ];
        let result = propagate_tension(&positions, &velocities, &segments);
        assert_eq!(result.tensions.len(), 2);
        // Each segment stretched by 1.0, tension = 100 * 1 = 100
        assert!((result.tensions[0] - 100.0).abs() < 1e-6);
        assert!((result.tensions[1] - 100.0).abs() < 1e-6);
        assert!(result.all_taut);
    }

    #[test]
    fn test_multi_segment_cable_creation() {
        let cable = MultiSegmentCable::new(
            [0.0, 10.0, 0.0],
            [10.0, 10.0, 0.0],
            5,
            1.0,
            1000.0,
            10.0,
            [0.0, -9.81, 0.0],
        );
        assert_eq!(cable.num_nodes(), 6);
        assert_eq!(cable.num_segments(), 5);
        assert!((cable.rest_length() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_segment_cable_step() {
        let mut cable = MultiSegmentCable::new(
            [0.0, 10.0, 0.0],
            [10.0, 10.0, 0.0],
            4,
            1.0,
            10000.0,
            100.0,
            [0.0, -9.81, 0.0],
        );
        // Step a few times — cable should sag under gravity
        for _i in 0..100 {
            cable.step(0.001);
        }
        // Middle nodes should have dropped below the initial height
        let mid = cable.positions[2][1];
        assert!(mid < 10.0, "mid_y={mid} should be < 10.0");
    }

    #[test]
    fn test_cable_pulley_system() {
        let mut system = CablePulleySystem::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 15.0);
        system.add_pulley(Pulley {
            center: [5.0, 5.0, 0.0],
            radius: 0.1,
            friction: 0.0,
            axis: [0.0, 0.0, 1.0],
        });
        let result = system.analyze(100.0);
        assert_eq!(result.segment_tensions.len(), 2);
        // No friction: tension should be same on both sides
        assert!((result.segment_tensions[0] - 100.0).abs() < 1e-6);
        assert!((result.segment_tensions[1] - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_cable_pulley_with_friction() {
        let mut system = CablePulleySystem::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 15.0);
        system.add_pulley(Pulley {
            center: [5.0, 5.0, 0.0],
            radius: 0.1,
            friction: 0.3,
            axis: [0.0, 0.0, 1.0],
        });
        let result = system.analyze(100.0);
        // With friction, output tension < input tension
        assert!(result.segment_tensions[1] < result.segment_tensions[0]);
    }

    #[test]
    fn test_cable_winding_drum() {
        let drum = CableWindingDrum::new(0.5, 0.3, 0.01);
        assert!(drum.max_capacity > 0.0);
        let state = drum.compute_state(1000.0);
        assert_eq!(state.layers, 0);
        assert!(state.required_torque > 0.0);
    }

    #[test]
    fn test_cable_winding_update() {
        let mut drum = CableWindingDrum::new(0.5, 0.3, 0.01);
        drum.winding_speed = 1.0;
        drum.update(5.0);
        assert!((drum.wound_length - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cable_attachment() {
        let attachment = CableAttachment::new([1.0, 0.0, 0.0], [5.0, 5.0, 0.0], 1.0);
        let world_pos = attachment.world_position();
        assert!((world_pos[0] - 6.0).abs() < 1e-10);
        assert!((world_pos[1] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cable_force_on_body() {
        let attachment = CableAttachment::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        let (force, _torque) = cable_force_on_body(&attachment, [10.0, 0.0, 0.0], 100.0);
        assert!((force[0] - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_elastic_potential_energy() {
        let pe = elastic_potential_energy(10.0, 8.0, 100.0);
        // 0.5 * 100 * 2^2 = 200
        assert!((pe - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_elastic_pe_no_stretch() {
        let pe = elastic_potential_energy(5.0, 8.0, 100.0);
        assert!(pe.abs() < 1e-10);
    }

    #[test]
    fn test_cable_natural_frequency() {
        let f = cable_natural_frequency(1.0, 100.0, 0.01);
        // f = 1/(2*1) * sqrt(100/0.01) = 0.5 * 100 = 50 Hz
        assert!((f - 50.0).abs() < 1e-6, "f={f}");
    }

    #[test]
    fn test_cable_wave_speed() {
        let c = cable_wave_speed(100.0, 0.01);
        assert!((c - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_parabolic_sag() {
        let sag = parabolic_sag(100.0, 10.0, 5000.0);
        // 10 * 100^2 / (8 * 5000) = 100000/40000 = 2.5
        assert!((sag - 2.5).abs() < 1e-6, "sag={sag}");
    }

    #[test]
    fn test_would_break() {
        assert!(would_break(1000.0, 2000.0, 3.0));
        assert!(!would_break(500.0, 2000.0, 3.0));
    }

    #[test]
    fn test_safety_margin() {
        let margin = safety_margin(500.0, 2000.0, 2.0);
        // allowable = 1000, margin = 1000/500 = 2.0
        assert!((margin - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_length_change() {
        let dl = thermal_length_change(100.0, 12e-6, 50.0);
        // 12e-6 * 100 * 50 = 0.06
        assert!((dl - 0.06).abs() < 1e-10);
    }

    #[test]
    fn test_logarithmic_decrement() {
        let delta = logarithmic_decrement(0.05);
        // 2*pi*0.05/sqrt(1-0.0025) ≈ 0.31438
        assert!((delta - 0.31438).abs() < 0.001, "delta={delta}");
    }

    #[test]
    fn test_ice_loading() {
        let w = ice_loading_weight(0.03, 0.01, 900.0, 9.81);
        assert!(w > 0.0);
    }

    #[test]
    fn test_fatigue_life() {
        let n = fatigue_life_cycles(100e6, 1e30, 3.0);
        assert!(n > 0.0);
        assert!(n.is_finite());
    }

    #[test]
    fn test_miners_rule() {
        let counts = [1000.0, 2000.0];
        let allowable = [10000.0, 5000.0];
        let damage = miners_rule_damage(&counts, &allowable);
        // 1000/10000 + 2000/5000 = 0.1 + 0.4 = 0.5
        assert!((damage - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_concentrated_load_deflection() {
        let defl = concentrated_load_deflection(10.0, 5.0, 100.0, 500.0);
        // 100 * 5 * 5 / (10 * 500) = 2500/5000 = 0.5
        assert!((defl - 0.5).abs() < 1e-10, "defl={defl}");
    }

    #[test]
    fn test_cable_drag_force_perpendicular() {
        let f = cable_drag_force(
            [10.0, 0.0, 0.0], // wind along x
            [0.0, 0.0, 1.0],  // cable along z
            0.05,
            10.0,
            1.225,
            1.0,
        );
        // Force should be along x (perpendicular to cable)
        assert!(f[0] > 0.0);
        assert!(f[2].abs() < 1e-10);
    }

    #[test]
    fn test_gravitational_pe() {
        let positions = vec![[0.0, 10.0, 0.0], [5.0, 10.0, 0.0], [10.0, 10.0, 0.0]];
        let pe = gravitational_potential_energy(&positions, 1.0, 9.81, 0.0);
        // Total length = 10m, total mass = 10kg, avg height = 10m
        // PE = 10 * 9.81 * 10 = 981
        assert!((pe - 981.0).abs() < 1.0, "pe={pe}");
    }

    #[test]
    fn test_kinetic_energy() {
        let velocities = vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let masses = vec![2.0, 3.0];
        let ke = kinetic_energy(&velocities, &masses);
        // 0.5*2*1 + 0.5*3*4 = 1 + 6 = 7
        assert!((ke - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_creep_strain_rate() {
        let rate = creep_strain_rate(100.0, 1e-10, 3.0);
        // 1e-10 * 100^3 = 1e-10 * 1e6 = 1e-4
        assert!((rate - 1e-4).abs() < 1e-10);
    }

    #[test]
    fn test_galloping_onset() {
        let v = galloping_onset_speed(10.0, std::f64::consts::TAU, 0.01, 1.225, 0.05, -0.5);
        assert!(v > 0.0);
        assert!(v.is_finite());
    }
}
