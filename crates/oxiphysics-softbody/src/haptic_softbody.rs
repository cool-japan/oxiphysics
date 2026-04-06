#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Haptic rendering for deformable soft bodies.
//!
//! Implements the god-object / proxy-based haptic rendering pipeline,
//! virtual coupling for haptic thread simulation, force feedback computation,
//! phantom haptic device model, needle insertion force modelling, deformable
//! tissue haptic simulation, virtual fixture constraints, latency compensation,
//! and multi-point haptic contact.
//!
//! # Overview
//!
//! Haptic rendering requires a high-frequency inner loop (~1 kHz) that
//! computes contact forces between a proxy (constrained version of the haptic
//! device cursor) and the deformable surface.  The outer simulation loop runs
//! at the normal physics rate (~60 Hz) and updates the soft-body deformation.
//! The two loops are coupled through a [`VirtualCoupling`] spring-damper.
//!
//! # References
//!
//! * Zilles & Salisbury (1995) – god-object method.
//! * Ruspini et al. (1997) – proxy-based haptic rendering.
//! * McNeely et al. (1999) – voxel-based haptic rendering.
//! * Misra et al. (2010) – mechanics-based models for needle insertion.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / n)
    }
}

#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Minimum haptic loop frequency (Hz) required for stable force feedback.
pub const HAPTIC_MIN_FREQUENCY_HZ: f64 = 500.0;

/// Typical phantom haptic device maximum force (N).
pub const PHANTOM_MAX_FORCE_N: f64 = 8.5;

/// Typical haptic rendering loop period at 1 kHz (s).
pub const HAPTIC_LOOP_DT: f64 = 1.0e-3;

// ---------------------------------------------------------------------------
// ProxyState
// ---------------------------------------------------------------------------

/// The proxy (god-object) represents the closest reachable position on or
/// outside the deformable surface.
///
/// The proxy position is constrained to lie outside (or on) the deformable
/// mesh, while the haptic cursor may be inside the object.  The difference
/// between cursor and proxy drives the force feedback signal.
#[derive(Debug, Clone)]
pub struct ProxyState {
    /// Current proxy position in world coordinates (m).
    pub position: [f64; 3],
    /// Proxy velocity (m/s), estimated by finite difference.
    pub velocity: [f64; 3],
    /// Surface normal at the proxy contact point (unit vector).
    pub surface_normal: [f64; 3],
    /// Penetration depth of the haptic cursor into the surface (m, ≥ 0).
    pub penetration_depth: f64,
    /// Whether the proxy is currently in contact with the surface.
    pub in_contact: bool,
}

impl ProxyState {
    /// Create a new proxy at `position` with no contact.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            surface_normal: [0.0, 1.0, 0.0],
            penetration_depth: 0.0,
            in_contact: false,
        }
    }

    /// Update the proxy position given the haptic cursor position, a surface
    /// contact point, and its outward normal.
    ///
    /// The god-object method projects the cursor onto the constraint manifold
    /// (the closest point on the surface).
    pub fn update(
        &mut self,
        cursor_pos: [f64; 3],
        contact_point: [f64; 3],
        contact_normal: [f64; 3],
        dt: f64,
    ) {
        let prev = self.position;
        // Signed distance from surface: negative means inside object.
        let d = dot3(sub3(cursor_pos, contact_point), contact_normal);
        if d < 0.0 {
            // Cursor is inside: clamp proxy to surface.
            self.position = add3(cursor_pos, scale3(contact_normal, -d));
            self.surface_normal = contact_normal;
            self.penetration_depth = (-d).max(0.0);
            self.in_contact = true;
        } else {
            self.position = cursor_pos;
            self.penetration_depth = 0.0;
            self.in_contact = false;
        }
        if dt > 1e-12 {
            self.velocity = scale3(sub3(self.position, prev), 1.0 / dt);
        }
    }
}

// ---------------------------------------------------------------------------
// HapticCursor
// ---------------------------------------------------------------------------

/// Haptic cursor representing the physical haptic device tip position.
///
/// Position is typically updated at the haptic device servo rate (≥ 1 kHz).
#[derive(Debug, Clone)]
pub struct HapticCursor {
    /// Current cursor position in world coordinates (m).
    pub position: [f64; 3],
    /// Current cursor velocity (m/s).
    pub velocity: [f64; 3],
    /// Rendered force fed back to the user (N), updated each haptic frame.
    pub feedback_force: [f64; 3],
}

impl HapticCursor {
    /// Construct a cursor at `position` with zero velocity.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            feedback_force: [0.0; 3],
        }
    }

    /// Update cursor state from device readings.
    pub fn set_state(&mut self, position: [f64; 3], velocity: [f64; 3]) {
        self.position = position;
        self.velocity = velocity;
    }
}

// ---------------------------------------------------------------------------
// ContactForceModel
// ---------------------------------------------------------------------------

/// Contact force model for haptic rendering using penalty-based stiffness.
///
/// The rendered force is:
///
/// ```text
/// F = k_s * penetration * n̂  -  k_d * (v_proxy · n̂) * n̂
/// ```
///
/// where `n̂` is the surface normal, `k_s` is surface stiffness, and
/// `k_d` is viscous damping.
#[derive(Debug, Clone)]
pub struct ContactForceModel {
    /// Surface stiffness (N/m).
    pub stiffness: f64,
    /// Viscous damping (N·s/m).
    pub damping: f64,
    /// Maximum renderable force magnitude (N) — clamps to device limit.
    pub max_force: f64,
}

impl ContactForceModel {
    /// Create a new contact force model.
    pub fn new(stiffness: f64, damping: f64, max_force: f64) -> Self {
        Self {
            stiffness,
            damping,
            max_force,
        }
    }

    /// Compute the haptic feedback force given proxy state.
    ///
    /// Returns the force vector (N) to be rendered to the user.
    pub fn compute_force(&self, proxy: &ProxyState) -> [f64; 3] {
        if !proxy.in_contact {
            return [0.0; 3];
        }
        let n = proxy.surface_normal;
        // Spring term.
        let f_spring = scale3(n, self.stiffness * proxy.penetration_depth);
        // Damping term (opposes normal velocity).
        let v_normal = dot3(proxy.velocity, n);
        let f_damp = scale3(n, -self.damping * v_normal);
        let f_total = add3(f_spring, f_damp);
        // Clamp to device limits.
        let mag = norm3(f_total);
        if mag > self.max_force {
            scale3(f_total, self.max_force / mag)
        } else {
            f_total
        }
    }

    /// Compute stiffness as the derivative of force w.r.t. penetration depth
    /// (i.e., just returns `stiffness`).
    pub fn local_stiffness(&self) -> f64 {
        self.stiffness
    }
}

// ---------------------------------------------------------------------------
// VirtualCoupling
// ---------------------------------------------------------------------------

/// Virtual coupling between the haptic device and the simulated tool.
///
/// Because the haptic loop runs at ~1 kHz and the physics loop at ~60 Hz,
/// a spring-damper "virtual coupling" element connects them.  The haptic
/// loop computes a coupling force based on the difference between device
/// position and simulated tool position.
///
/// The coupling force is:
///
/// ```text
/// F_coupling = K_c * (x_device - x_tool)  +  B_c * (v_device - v_tool)
/// ```
#[derive(Debug, Clone)]
pub struct VirtualCoupling {
    /// Coupling stiffness K_c (N/m).
    pub stiffness: f64,
    /// Coupling damping B_c (N·s/m).
    pub damping: f64,
    /// Simulated tool position (m).
    pub tool_position: [f64; 3],
    /// Simulated tool velocity (m/s).
    pub tool_velocity: [f64; 3],
}

impl VirtualCoupling {
    /// Create a new virtual coupling element.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self {
            stiffness,
            damping,
            tool_position: [0.0; 3],
            tool_velocity: [0.0; 3],
        }
    }

    /// Compute the coupling force from device state to tool.
    ///
    /// * `device_pos` – haptic device tip position (m).
    /// * `device_vel` – haptic device tip velocity (m/s).
    ///
    /// Returns the force (N) applied to the simulated tool.
    pub fn coupling_force(&self, device_pos: [f64; 3], device_vel: [f64; 3]) -> [f64; 3] {
        let dp = sub3(device_pos, self.tool_position);
        let dv = sub3(device_vel, self.tool_velocity);
        add3(scale3(dp, self.stiffness), scale3(dv, self.damping))
    }

    /// Update tool state by integrating the coupling force.
    ///
    /// Uses explicit Euler with mass `m_tool` (kg) and time step `dt` (s).
    pub fn integrate_tool(
        &mut self,
        device_pos: [f64; 3],
        device_vel: [f64; 3],
        m_tool: f64,
        dt: f64,
    ) {
        let f = self.coupling_force(device_pos, device_vel);
        let accel = scale3(f, 1.0 / m_tool);
        self.tool_velocity = add3(self.tool_velocity, scale3(accel, dt));
        self.tool_position = add3(self.tool_position, scale3(self.tool_velocity, dt));
    }

    /// Compute the stability bound for the coupling stiffness.
    ///
    /// For a mass-spring system the maximum stable stiffness is:
    /// `K_max = 2 * m / dt^2`.
    pub fn max_stable_stiffness(m_tool: f64, dt: f64) -> f64 {
        2.0 * m_tool / (dt * dt)
    }
}

// ---------------------------------------------------------------------------
// PhantomHapticDevice
// ---------------------------------------------------------------------------

/// Simplified model of a Phantom Omni/Premium haptic device.
///
/// Models the kinematics as a 6-DOF device with position workspace limits
/// and force saturation.  Provides a simplified gravity compensation model.
#[derive(Debug, Clone)]
pub struct PhantomHapticDevice {
    /// Workspace half-extents \[x, y, z\] (m).
    pub workspace: [f64; 3],
    /// Maximum output force (N).
    pub max_force: f64,
    /// Maximum output torque (N·m).
    pub max_torque: f64,
    /// Device stylus position (m).
    pub position: [f64; 3],
    /// Device stylus velocity (m/s).
    pub velocity: [f64; 3],
    /// Device stylus orientation as a rotation vector (rad).
    pub orientation: [f64; 3],
    /// Gravity compensation enabled flag.
    pub gravity_compensation: bool,
    /// Effective stylus mass for gravity compensation (kg).
    pub stylus_mass: f64,
}

impl PhantomHapticDevice {
    /// Create a new Phantom device model with default workspace ±150 mm.
    pub fn new() -> Self {
        Self {
            workspace: [0.15, 0.15, 0.15],
            max_force: PHANTOM_MAX_FORCE_N,
            max_torque: 0.188,
            position: [0.0; 3],
            velocity: [0.0; 3],
            orientation: [0.0; 3],
            gravity_compensation: true,
            stylus_mass: 0.050,
        }
    }

    /// Update device state from sensor readings.
    pub fn update_state(&mut self, position: [f64; 3], velocity: [f64; 3], orientation: [f64; 3]) {
        // Clamp to workspace.
        for i in 0..3 {
            self.position[i] = clamp(position[i], -self.workspace[i], self.workspace[i]);
        }
        self.velocity = velocity;
        self.orientation = orientation;
    }

    /// Compute the gravity compensation torque vector.
    ///
    /// Returns approximate compensation force (N) to counteract stylus weight.
    pub fn gravity_compensation_force(&self) -> [f64; 3] {
        if self.gravity_compensation {
            [0.0, self.stylus_mass * 9.81, 0.0]
        } else {
            [0.0; 3]
        }
    }

    /// Saturate force command to device limits.
    pub fn saturate_force(&self, force: [f64; 3]) -> [f64; 3] {
        let mag = norm3(force);
        if mag > self.max_force {
            scale3(force, self.max_force / mag)
        } else {
            force
        }
    }

    /// Check whether the device position is inside the workspace.
    pub fn in_workspace(&self, pos: [f64; 3]) -> bool {
        pos[0].abs() <= self.workspace[0]
            && pos[1].abs() <= self.workspace[1]
            && pos[2].abs() <= self.workspace[2]
    }
}

impl Default for PhantomHapticDevice {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// NeedleInsertionModel
// ---------------------------------------------------------------------------

/// Force model for needle insertion into soft tissue.
///
/// Models the three main phases of needle-tissue interaction:
/// 1. Pre-puncture: tissue deformation before puncture (elastic loading).
/// 2. Puncture: tissue rupture event (force drop).
/// 3. Post-puncture: needle advancement (cutting + friction forces).
///
/// Reference: Okamura et al. (2004), Misra et al. (2010).
#[derive(Debug, Clone)]
pub struct NeedleInsertionModel {
    /// Tissue stiffness before puncture (N/m).
    pub pre_puncture_stiffness: f64,
    /// Maximum force at puncture (N).
    pub puncture_force: f64,
    /// Post-puncture cutting force (N) per unit insertion depth.
    pub cutting_stiffness: f64,
    /// Needle-tissue friction coefficient (dimensionless).
    pub friction_coeff: f64,
    /// Needle shaft diameter (m).
    pub needle_diameter: f64,
    /// Current insertion depth (m).
    pub insertion_depth: f64,
    /// Whether the needle has punctured the surface.
    pub punctured: bool,
    /// Depth at which puncture occurred (m).
    pub puncture_depth: f64,
}

impl NeedleInsertionModel {
    /// Create a new needle insertion model for soft tissue.
    pub fn new(
        pre_puncture_stiffness: f64,
        puncture_force: f64,
        cutting_stiffness: f64,
        friction_coeff: f64,
        needle_diameter: f64,
    ) -> Self {
        Self {
            pre_puncture_stiffness,
            puncture_force,
            cutting_stiffness,
            friction_coeff,
            needle_diameter,
            insertion_depth: 0.0,
            punctured: false,
            puncture_depth: 0.0,
        }
    }

    /// Compute the axial insertion force at the current depth.
    ///
    /// Returns force (N) in the needle axis direction (positive = resistance).
    pub fn axial_force(&self) -> f64 {
        if !self.punctured {
            // Pre-puncture: linear spring loading.
            self.pre_puncture_stiffness * self.insertion_depth
        } else {
            // Post-puncture: cutting + friction.
            let depth_beyond = (self.insertion_depth - self.puncture_depth).max(0.0);
            let f_cut = self.cutting_stiffness * depth_beyond;
            // Friction force proportional to needle circumference × contact length.
            let circumference = std::f64::consts::PI * self.needle_diameter;
            let f_friction = self.friction_coeff * circumference * depth_beyond * 1000.0; // Pa/m × m
            f_cut + f_friction
        }
    }

    /// Advance the needle by `delta` metres and update puncture state.
    ///
    /// If pre-puncture force exceeds `puncture_force`, the tissue punctures.
    pub fn advance(&mut self, delta: f64) {
        self.insertion_depth += delta.max(0.0);
        if !self.punctured {
            let f = self.pre_puncture_stiffness * self.insertion_depth;
            if f >= self.puncture_force {
                self.punctured = true;
                self.puncture_depth = self.insertion_depth;
            }
        }
    }

    /// Retract the needle by `delta` metres (hysteresis: friction reverses).
    pub fn retract(&mut self, delta: f64) {
        self.insertion_depth = (self.insertion_depth - delta.max(0.0)).max(0.0);
        if self.insertion_depth < 1e-6 {
            self.punctured = false;
            self.puncture_depth = 0.0;
        }
    }

    /// Compute the lateral tissue reaction force (N) for a needle bent by
    /// `lateral_deflection` (m) at the tip.
    pub fn lateral_force(&self, lateral_deflection: f64) -> f64 {
        // Simplified beam model: EI/L^2 × deflection.
        // Use pre-puncture stiffness as a proxy for tissue lateral stiffness.
        0.1 * self.pre_puncture_stiffness * lateral_deflection
    }
}

// ---------------------------------------------------------------------------
// DeformableTissueHaptics
// ---------------------------------------------------------------------------

/// Haptic simulation layer for deformable tissue.
///
/// Combines the proxy-based haptic rendering with a linearised tissue
/// deformation model.  The tissue is modelled as a grid of Kelvin–Voigt
/// elements (spring + damper in parallel).
#[derive(Debug, Clone)]
pub struct DeformableTissueHaptics {
    /// Number of tissue nodes in X.
    pub nx: usize,
    /// Number of tissue nodes in Y.
    pub ny: usize,
    /// Tissue Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Tissue viscosity (Pa·s).
    pub viscosity: f64,
    /// Node spacing (m).
    pub node_spacing: f64,
    /// Node displacements (m) — linearised vertical deformation.
    pub displacements: Vec<f64>,
    /// Node velocities (m/s).
    pub velocities: Vec<f64>,
    /// Node rest positions (x, y) in the surface plane.
    pub rest_positions: Vec<[f64; 2]>,
}

impl DeformableTissueHaptics {
    /// Create a new tissue haptic model on an `nx × ny` grid.
    pub fn new(
        nx: usize,
        ny: usize,
        node_spacing: f64,
        youngs_modulus: f64,
        viscosity: f64,
    ) -> Self {
        let n = nx * ny;
        let mut rest_positions = Vec::with_capacity(n);
        for j in 0..ny {
            for i in 0..nx {
                rest_positions.push([i as f64 * node_spacing, j as f64 * node_spacing]);
            }
        }
        Self {
            nx,
            ny,
            youngs_modulus,
            viscosity,
            node_spacing,
            displacements: vec![0.0; n],
            velocities: vec![0.0; n],
            rest_positions,
        }
    }

    /// Apply a point load `force_z` (N) at world position `(px, py)`.
    ///
    /// Distributes force to the four surrounding nodes using bilinear weights.
    pub fn apply_point_load(&mut self, px: f64, py: f64, force_z: f64) {
        let inv_h = 1.0 / self.node_spacing;
        let ix = (px * inv_h).floor() as isize;
        let iy = (py * inv_h).floor() as isize;
        let tx = (px * inv_h) - ix as f64;
        let ty = (py * inv_h) - iy as f64;
        let nx = self.nx as isize;
        let ny = self.ny as isize;
        // Four corner weights (bilinear interpolation).
        let corners = [
            (ix, iy, (1.0 - tx) * (1.0 - ty)),
            (ix + 1, iy, tx * (1.0 - ty)),
            (ix, iy + 1, (1.0 - tx) * ty),
            (ix + 1, iy + 1, tx * ty),
        ];
        let area = self.node_spacing * self.node_spacing;
        for (ci, cj, w) in corners {
            if ci >= 0 && ci < nx && cj >= 0 && cj < ny {
                let idx = cj as usize * self.nx + ci as usize;
                // Distributed force per unit area.
                let pressure = force_z * w / area;
                // Kelvin–Voigt equilibrium: stress = E * strain ≈ E * disp / h.
                let stiffness = self.youngs_modulus / self.node_spacing;
                self.displacements[idx] += pressure / stiffness;
            }
        }
    }

    /// Step the tissue dynamics forward by `dt` (s).
    ///
    /// Each node is modelled as an independent Kelvin–Voigt element that
    /// relaxes back towards zero displacement with time constant `η/E`.
    pub fn step(&mut self, dt: f64) {
        let tau = self.viscosity / self.youngs_modulus;
        let decay = (-dt / tau).exp();
        for i in 0..self.displacements.len() {
            self.displacements[i] *= decay;
            self.velocities[i] = self.displacements[i] * (decay - 1.0) / dt;
        }
    }

    /// Query the surface height at position `(px, py)` by bilinear
    /// interpolation of node displacements.
    pub fn surface_height(&self, px: f64, py: f64) -> f64 {
        let inv_h = 1.0 / self.node_spacing;
        let ix = (px * inv_h).floor() as isize;
        let iy = (py * inv_h).floor() as isize;
        let tx = (px * inv_h) - ix as f64;
        let ty = (py * inv_h) - iy as f64;
        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let corners = [
            (ix, iy, (1.0 - tx) * (1.0 - ty)),
            (ix + 1, iy, tx * (1.0 - ty)),
            (ix, iy + 1, (1.0 - tx) * ty),
            (ix + 1, iy + 1, tx * ty),
        ];
        let mut h = 0.0;
        for (ci, cj, w) in corners {
            if ci >= 0 && ci < nx && cj >= 0 && cj < ny {
                let idx = cj as usize * self.nx + ci as usize;
                h += w * self.displacements[idx];
            }
        }
        h
    }
}

// ---------------------------------------------------------------------------
// VirtualFixture
// ---------------------------------------------------------------------------

/// Virtual fixture constraint for haptic guidance.
///
/// Virtual fixtures are software-generated constraints that either:
/// * Guide the user along a desired path (assistive fixture).
/// * Prevent motion into forbidden regions (forbidden-region fixture).
///
/// Reference: Rosenberg (1993).
#[derive(Debug, Clone)]
pub struct VirtualFixture {
    /// Fixture type.
    pub fixture_type: FixtureType,
    /// Fixture stiffness for guidance / barrier (N/m).
    pub stiffness: f64,
    /// Fixture damping (N·s/m).
    pub damping: f64,
}

/// Type of virtual fixture constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum FixtureType {
    /// Forbidden region — half-space defined by a point and outward normal.
    ForbiddenHalfSpace {
        /// Point on the boundary plane (m).
        point: [f64; 3],
        /// Outward normal (unit vector, pointing into forbidden region).
        normal: [f64; 3],
    },
    /// Linear guide — constrains motion to a line.
    LinearGuide {
        /// Point on the guide line (m).
        origin: [f64; 3],
        /// Unit direction of the guide line.
        direction: [f64; 3],
    },
    /// Planar guide — constrains motion to a plane.
    PlanarGuide {
        /// Point on the guide plane (m).
        origin: [f64; 3],
        /// Plane normal (unit vector).
        normal: [f64; 3],
    },
}

impl VirtualFixture {
    /// Create a new virtual fixture.
    pub fn new(fixture_type: FixtureType, stiffness: f64, damping: f64) -> Self {
        Self {
            fixture_type,
            stiffness,
            damping,
        }
    }

    /// Compute the constraint force to apply to the haptic cursor.
    ///
    /// * `pos` – current cursor position (m).
    /// * `vel` – current cursor velocity (m/s).
    ///
    /// Returns the force (N) in world coordinates.
    pub fn constraint_force(&self, pos: [f64; 3], vel: [f64; 3]) -> [f64; 3] {
        match &self.fixture_type {
            FixtureType::ForbiddenHalfSpace { point, normal } => {
                let d = dot3(sub3(pos, *point), *normal);
                if d > 0.0 {
                    // Inside forbidden region — push back.
                    let f_spring = scale3(*normal, -self.stiffness * d);
                    let v_n = dot3(vel, *normal);
                    let f_damp = scale3(*normal, -self.damping * v_n);
                    add3(f_spring, f_damp)
                } else {
                    [0.0; 3]
                }
            }
            FixtureType::LinearGuide { origin, direction } => {
                // Project position onto the line; error is the off-line component.
                let dp = sub3(pos, *origin);
                let along = dot3(dp, *direction);
                let on_line = add3(*origin, scale3(*direction, along));
                let error = sub3(pos, on_line);
                let v_perp = sub3(vel, scale3(*direction, dot3(vel, *direction)));
                let f_spring = scale3(error, -self.stiffness);
                let f_damp = scale3(v_perp, -self.damping);
                add3(f_spring, f_damp)
            }
            FixtureType::PlanarGuide { origin, normal } => {
                let d = dot3(sub3(pos, *origin), *normal);
                let f_spring = scale3(*normal, -self.stiffness * d);
                let v_n = dot3(vel, *normal);
                let f_damp = scale3(*normal, -self.damping * v_n);
                add3(f_spring, f_damp)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// LatencyCompensator
// ---------------------------------------------------------------------------

/// Latency compensator using a first-order predictor.
///
/// In networked haptic systems the round-trip delay can be several
/// milliseconds.  This predictor extrapolates position forward by the
/// estimated latency using constant-velocity assumption.
#[derive(Debug, Clone)]
pub struct LatencyCompensator {
    /// Estimated one-way network latency (s).
    pub latency: f64,
    /// History buffer of (time, position) pairs.
    history: Vec<(f64, [f64; 3])>,
    /// Maximum history length.
    max_history: usize,
}

impl LatencyCompensator {
    /// Create a new latency compensator.
    pub fn new(latency: f64) -> Self {
        Self {
            latency,
            history: Vec::with_capacity(64),
            max_history: 64,
        }
    }

    /// Record a position sample at time `t`.
    pub fn record(&mut self, t: f64, pos: [f64; 3]) {
        if self.history.len() >= self.max_history {
            self.history.remove(0);
        }
        self.history.push((t, pos));
    }

    /// Predict the position at `t + latency` using linear extrapolation.
    pub fn predict(&self, t: f64) -> [f64; 3] {
        if self.history.len() < 2 {
            return self.history.last().map(|e| e.1).unwrap_or([0.0; 3]);
        }
        let n = self.history.len();
        let (t1, p1) = self.history[n - 2];
        let (t2, p2) = self.history[n - 1];
        let dt_hist = t2 - t1;
        if dt_hist < 1e-12 {
            return p2;
        }
        let vel = scale3(sub3(p2, p1), 1.0 / dt_hist);
        let dt_pred = (t + self.latency) - t2;
        add3(p2, scale3(vel, dt_pred))
    }
}

// ---------------------------------------------------------------------------
// MultiPointHapticContact
// ---------------------------------------------------------------------------

/// Multi-point haptic contact for tools with multiple contact zones.
///
/// Aggregates force contributions from several contact points (e.g., a
/// forceps with two jaws, or a needle with shaft contacts).
#[derive(Debug, Clone)]
pub struct MultiPointHapticContact {
    /// Individual contact points and their local force models.
    pub contacts: Vec<ContactPoint>,
    /// Maximum combined force (N) — clamps total feedback.
    pub max_combined_force: f64,
}

/// A single contact point with its own stiffness/damping.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// Local identifier for this contact zone.
    pub id: usize,
    /// Contact stiffness (N/m).
    pub stiffness: f64,
    /// Contact damping (N·s/m).
    pub damping: f64,
    /// Current penetration depth (m).
    pub penetration: f64,
    /// Surface normal at this contact.
    pub normal: [f64; 3],
    /// Contact velocity (m/s).
    pub contact_velocity: [f64; 3],
    /// Whether this contact is active.
    pub active: bool,
}

impl ContactPoint {
    /// Create a new contact point.
    pub fn new(id: usize, stiffness: f64, damping: f64) -> Self {
        Self {
            id,
            stiffness,
            damping,
            penetration: 0.0,
            normal: [0.0, 1.0, 0.0],
            contact_velocity: [0.0; 3],
            active: false,
        }
    }

    /// Compute the force contribution from this contact point.
    pub fn force(&self) -> [f64; 3] {
        if !self.active || self.penetration <= 0.0 {
            return [0.0; 3];
        }
        let f_spring = scale3(self.normal, self.stiffness * self.penetration);
        let v_n = dot3(self.contact_velocity, self.normal);
        let f_damp = scale3(self.normal, -self.damping * v_n);
        add3(f_spring, f_damp)
    }
}

impl MultiPointHapticContact {
    /// Create a multi-point haptic contact manager.
    pub fn new(max_combined_force: f64) -> Self {
        Self {
            contacts: Vec::new(),
            max_combined_force,
        }
    }

    /// Add a contact point.
    pub fn add_contact(&mut self, cp: ContactPoint) {
        self.contacts.push(cp);
    }

    /// Compute the total haptic feedback force from all active contacts.
    pub fn total_force(&self) -> [f64; 3] {
        let mut total = [0.0; 3];
        for cp in &self.contacts {
            total = add3(total, cp.force());
        }
        let mag = norm3(total);
        if mag > self.max_combined_force {
            scale3(total, self.max_combined_force / mag)
        } else {
            total
        }
    }

    /// Return the number of currently active contact points.
    pub fn active_count(&self) -> usize {
        self.contacts.iter().filter(|c| c.active).count()
    }
}

// ---------------------------------------------------------------------------
// HapticRenderingLoop
// ---------------------------------------------------------------------------

/// High-frequency haptic rendering loop controller.
///
/// Manages the 1 kHz haptic servo loop, integrating proxy updates,
/// virtual coupling, contact force computation, and latency compensation.
#[derive(Debug, Clone)]
pub struct HapticRenderingLoop {
    /// Haptic cursor.
    pub cursor: HapticCursor,
    /// Proxy state.
    pub proxy: ProxyState,
    /// Contact force model.
    pub contact_model: ContactForceModel,
    /// Virtual coupling to physics simulation.
    pub coupling: VirtualCoupling,
    /// Current simulation time (s).
    pub time: f64,
    /// Haptic loop time step (s).
    pub dt: f64,
    /// Iteration counter (for diagnostics).
    pub iteration: u64,
}

impl HapticRenderingLoop {
    /// Create a new haptic rendering loop.
    pub fn new(
        stiffness: f64,
        damping: f64,
        coupling_stiffness: f64,
        coupling_damping: f64,
        dt: f64,
    ) -> Self {
        Self {
            cursor: HapticCursor::new([0.0; 3]),
            proxy: ProxyState::new([0.0; 3]),
            contact_model: ContactForceModel::new(stiffness, damping, PHANTOM_MAX_FORCE_N),
            coupling: VirtualCoupling::new(coupling_stiffness, coupling_damping),
            time: 0.0,
            dt,
            iteration: 0,
        }
    }

    /// Execute one haptic servo tick.
    ///
    /// * `device_pos` – raw haptic device position (m).
    /// * `device_vel` – raw haptic device velocity (m/s).
    /// * `contact_point` – closest point on deformable surface (m).
    /// * `contact_normal` – outward normal at contact (unit vector).
    ///
    /// Returns the force vector (N) to command to the haptic device.
    pub fn tick(
        &mut self,
        device_pos: [f64; 3],
        device_vel: [f64; 3],
        contact_point: [f64; 3],
        contact_normal: [f64; 3],
    ) -> [f64; 3] {
        self.cursor.set_state(device_pos, device_vel);
        self.proxy
            .update(device_pos, contact_point, contact_normal, self.dt);
        let f_contact = self.contact_model.compute_force(&self.proxy);
        // Virtual coupling force from device to simulated tool.
        let f_coupling = self.coupling.coupling_force(device_pos, device_vel);
        // Total force = contact + coupling (coupling feeds back to user).
        let f_total = add3(f_contact, f_coupling);
        self.cursor.feedback_force = self.contact_model.compute_force(&self.proxy); // store for diagnostics
        self.time += self.dt;
        self.iteration += 1;
        f_total
    }
}

// ---------------------------------------------------------------------------
// HapticFrequencyAnalysis
// ---------------------------------------------------------------------------

/// Utilities for analysing haptic rendering frequency requirements.
///
/// A haptic system is stable if the loop frequency satisfies:
/// `f_haptic > f_cutoff = sqrt(k_s / m_device) / (2 * π)`
/// to avoid the Z-width problem.
pub struct HapticFrequencyAnalysis;

impl HapticFrequencyAnalysis {
    /// Compute the minimum stable haptic loop frequency for a given
    /// surface stiffness and device effective mass.
    ///
    /// Returns minimum frequency (Hz).
    pub fn min_stable_frequency(surface_stiffness: f64, device_mass: f64) -> f64 {
        let omega = (surface_stiffness / device_mass).sqrt();
        omega / (2.0 * std::f64::consts::PI)
    }

    /// Compute the maximum renderable stiffness for a given loop frequency
    /// and device mass (`Z-width` bound).
    ///
    /// `k_max = m * (2 * π * f)^2`
    pub fn max_renderable_stiffness(frequency_hz: f64, device_mass: f64) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * frequency_hz;
        device_mass * omega * omega
    }

    /// Estimate the required loop period `dt` to stably render stiffness `k`
    /// with device mass `m`.
    pub fn required_dt(surface_stiffness: f64, device_mass: f64) -> f64 {
        let f_min = Self::min_stable_frequency(surface_stiffness, device_mass);
        if f_min < 1e-6 {
            1.0
        } else {
            1.0 / (2.0 * f_min)
        }
    }
}

// ---------------------------------------------------------------------------
// TissueLayerModel
// ---------------------------------------------------------------------------

/// Multi-layer tissue model for haptic simulation of layered anatomy.
///
/// Models skin, fat, muscle, and organ layers each with distinct stiffness
/// and thickness parameters.
#[derive(Debug, Clone)]
pub struct TissueLayerModel {
    /// Tissue layers ordered from surface inward.
    pub layers: Vec<TissueLayer>,
}

/// A single tissue layer.
#[derive(Debug, Clone)]
pub struct TissueLayer {
    /// Layer name (e.g., "skin", "fat", "muscle").
    pub name: &'static str,
    /// Layer thickness (m).
    pub thickness: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Viscosity (Pa·s).
    pub viscosity: f64,
    /// Current compression (m).
    pub compression: f64,
}

impl TissueLayer {
    /// Create a new tissue layer.
    pub fn new(name: &'static str, thickness: f64, youngs_modulus: f64, viscosity: f64) -> Self {
        Self {
            name,
            thickness,
            youngs_modulus,
            viscosity,
            compression: 0.0,
        }
    }

    /// Compute the reaction force for a given compression `delta` (m).
    pub fn reaction_force(&self, delta: f64) -> f64 {
        let strain = (delta / self.thickness).min(1.0);
        // Nonlinear stiffness: stiffer near full compression.
        let k_eff = self.youngs_modulus * (1.0 + 10.0 * strain * strain);
        k_eff * delta
    }
}

impl TissueLayerModel {
    /// Create a model with the given layers.
    pub fn new(layers: Vec<TissueLayer>) -> Self {
        Self { layers }
    }

    /// Create a default 4-layer anatomical model (skin / fat / muscle / organ).
    pub fn default_anatomical() -> Self {
        Self {
            layers: vec![
                TissueLayer::new("skin", 0.002, 100_000.0, 50.0),
                TissueLayer::new("fat", 0.020, 5_000.0, 20.0),
                TissueLayer::new("muscle", 0.030, 50_000.0, 200.0),
                TissueLayer::new("organ", 0.050, 20_000.0, 100.0),
            ],
        }
    }

    /// Compute the total reaction force for a tool tip penetrating to depth
    /// `total_depth` (m) from the surface.
    ///
    /// Compresses layers sequentially until total depth is consumed.
    pub fn total_reaction_force(&self, total_depth: f64) -> f64 {
        let mut remaining = total_depth;
        let mut force = 0.0;
        for layer in &self.layers {
            if remaining <= 0.0 {
                break;
            }
            let delta = remaining.min(layer.thickness);
            force += layer.reaction_force(delta);
            remaining -= delta;
        }
        force
    }

    /// Compute the depth of the tool tip within a specific layer index.
    pub fn depth_in_layer(&self, total_depth: f64) -> (usize, f64) {
        let mut remaining = total_depth;
        for (i, layer) in self.layers.iter().enumerate() {
            if remaining <= layer.thickness {
                return (i, remaining);
            }
            remaining -= layer.thickness;
        }
        (self.layers.len().saturating_sub(1), remaining)
    }
}

// ---------------------------------------------------------------------------
// HapticSimulationStats
// ---------------------------------------------------------------------------

/// Runtime statistics for the haptic simulation loop.
#[derive(Debug, Clone, Default)]
pub struct HapticSimulationStats {
    /// Total number of haptic ticks executed.
    pub tick_count: u64,
    /// Maximum observed force magnitude (N).
    pub max_force_n: f64,
    /// Average force magnitude (N) over the session.
    pub avg_force_n: f64,
    /// Number of force clamp events (force reached device limit).
    pub clamp_events: u64,
    /// Simulation time elapsed (s).
    pub elapsed_time: f64,
}

impl HapticSimulationStats {
    /// Record a new force sample.
    pub fn record_force(&mut self, force_mag: f64, clamped: bool) {
        self.tick_count += 1;
        if force_mag > self.max_force_n {
            self.max_force_n = force_mag;
        }
        self.avg_force_n += (force_mag - self.avg_force_n) / self.tick_count as f64;
        if clamped {
            self.clamp_events += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Unit Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // -----------------------------------------------------------------------
    // Math helper tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dot3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot3(a, b) - 32.0).abs() < EPS);
    }

    #[test]
    fn test_cross3() {
        let i = [1.0_f64, 0.0, 0.0];
        let j = [0.0_f64, 1.0, 0.0];
        let k = cross3(i, j);
        assert!(k[0].abs() < EPS && k[1].abs() < EPS && (k[2] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_normalize3_unit_length() {
        let v = [3.0_f64, 4.0, 0.0];
        let u = normalize3(v);
        let n = norm3(u);
        assert!((n - 1.0).abs() < EPS, "norm={n}");
    }

    #[test]
    fn test_normalize3_zero_returns_zero() {
        let v = [0.0_f64; 3];
        let u = normalize3(v);
        assert_eq!(u, [0.0; 3]);
    }

    // -----------------------------------------------------------------------
    // ProxyState tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_proxy_no_penetration() {
        let mut proxy = ProxyState::new([0.0; 3]);
        // Cursor above surface: no contact.
        proxy.update([0.0, 0.1, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.001);
        assert!(!proxy.in_contact, "should not be in contact");
        assert!(proxy.penetration_depth < EPS);
    }

    #[test]
    fn test_proxy_penetration_detected() {
        let mut proxy = ProxyState::new([0.0; 3]);
        // Cursor 5 mm inside surface.
        proxy.update([0.0, -0.005, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.001);
        assert!(proxy.in_contact, "should be in contact");
        assert!(
            (proxy.penetration_depth - 0.005).abs() < 1e-6,
            "depth={}",
            proxy.penetration_depth
        );
    }

    #[test]
    fn test_proxy_position_clamped_to_surface() {
        let mut proxy = ProxyState::new([0.0; 3]);
        proxy.update([0.0, -0.01, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.001);
        // Proxy should be at y=0 (surface), not inside.
        assert!(
            proxy.position[1].abs() < EPS,
            "proxy y={}",
            proxy.position[1]
        );
    }

    // -----------------------------------------------------------------------
    // ContactForceModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_contact_force_no_contact() {
        let model = ContactForceModel::new(1000.0, 10.0, 8.5);
        let proxy = ProxyState::new([0.0; 3]);
        let f = model.compute_force(&proxy);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_contact_force_proportional_to_penetration() {
        let model = ContactForceModel::new(1000.0, 0.0, 100.0);
        let mut proxy = ProxyState::new([0.0; 3]);
        proxy.in_contact = true;
        proxy.surface_normal = [0.0, 1.0, 0.0];
        proxy.penetration_depth = 0.005;
        proxy.velocity = [0.0; 3];
        let f = model.compute_force(&proxy);
        // F = 1000 * 0.005 = 5 N in +Y direction.
        assert!((f[1] - 5.0).abs() < EPS, "f_y={}", f[1]);
    }

    #[test]
    fn test_contact_force_clamped_to_max() {
        let model = ContactForceModel::new(100_000.0, 0.0, 8.5);
        let mut proxy = ProxyState::new([0.0; 3]);
        proxy.in_contact = true;
        proxy.surface_normal = [0.0, 1.0, 0.0];
        proxy.penetration_depth = 1.0; // huge penetration
        proxy.velocity = [0.0; 3];
        let f = model.compute_force(&proxy);
        let mag = norm3(f);
        assert!(
            (mag - 8.5).abs() < 1e-6,
            "should be clamped to 8.5 N, got {mag}"
        );
    }

    // -----------------------------------------------------------------------
    // VirtualCoupling tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_virtual_coupling_zero_error_zero_force() {
        let vc = VirtualCoupling::new(500.0, 10.0);
        let f = vc.coupling_force([0.0; 3], [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_virtual_coupling_spring_force() {
        let vc = VirtualCoupling::new(500.0, 0.0);
        // Device 10 mm ahead of tool.
        let f = vc.coupling_force([0.01, 0.0, 0.0], [0.0; 3]);
        // F = 500 * 0.01 = 5 N in +X.
        assert!((f[0] - 5.0).abs() < EPS, "f_x={}", f[0]);
    }

    #[test]
    fn test_virtual_coupling_max_stable_stiffness() {
        let k_max = VirtualCoupling::max_stable_stiffness(0.05, 0.001);
        let expected = 2.0 * 0.05 / (0.001 * 0.001);
        assert!((k_max - expected).abs() < EPS);
    }

    #[test]
    fn test_virtual_coupling_integrate_tool_moves() {
        let mut vc = VirtualCoupling::new(1000.0, 5.0);
        let init_pos = vc.tool_position;
        vc.integrate_tool([0.01, 0.0, 0.0], [0.0; 3], 0.05, 0.001);
        assert!(
            norm3(sub3(vc.tool_position, init_pos)) > 0.0,
            "tool should have moved"
        );
    }

    // -----------------------------------------------------------------------
    // PhantomHapticDevice tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_phantom_in_workspace() {
        let device = PhantomHapticDevice::new();
        assert!(device.in_workspace([0.0; 3]));
        assert!(device.in_workspace([0.14, 0.14, 0.14]));
        assert!(!device.in_workspace([0.2, 0.0, 0.0]));
    }

    #[test]
    fn test_phantom_saturate_force() {
        let device = PhantomHapticDevice::new();
        let f_big = [100.0, 0.0, 0.0];
        let f_sat = device.saturate_force(f_big);
        let mag = norm3(f_sat);
        assert!((mag - PHANTOM_MAX_FORCE_N).abs() < 1e-6, "mag={mag}");
    }

    #[test]
    fn test_phantom_gravity_compensation_enabled() {
        let device = PhantomHapticDevice::new();
        let f_gc = device.gravity_compensation_force();
        assert!(f_gc[1] > 0.0, "gravity comp should push up");
    }

    // -----------------------------------------------------------------------
    // NeedleInsertionModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_needle_pre_puncture_linear() {
        let mut needle = NeedleInsertionModel::new(5000.0, 2.0, 100.0, 0.3, 0.001);
        needle.advance(0.0002); // 0.2 mm, force = 5000 * 0.0002 = 1 N < 2 N
        assert!(!needle.punctured);
        let f = needle.axial_force();
        assert!((f - 5000.0 * 0.0002).abs() < 1e-6);
    }

    #[test]
    fn test_needle_puncture_event() {
        let mut needle = NeedleInsertionModel::new(5000.0, 0.5, 100.0, 0.3, 0.001);
        // Advance enough that force exceeds puncture_force.
        needle.advance(0.001); // force = 5000 * 0.001 = 5 N > 0.5 N
        assert!(needle.punctured, "needle should have punctured");
    }

    #[test]
    fn test_needle_retract_resets_puncture() {
        let mut needle = NeedleInsertionModel::new(5000.0, 0.5, 100.0, 0.3, 0.001);
        needle.advance(0.01);
        assert!(needle.punctured);
        needle.retract(0.02); // retract past start
        assert!(!needle.punctured, "retraction past start resets puncture");
    }

    #[test]
    fn test_needle_lateral_force() {
        let needle = NeedleInsertionModel::new(5000.0, 2.0, 100.0, 0.3, 0.001);
        let f = needle.lateral_force(0.001);
        assert!(f > 0.0, "lateral force should be positive");
    }

    // -----------------------------------------------------------------------
    // DeformableTissueHaptics tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_tissue_haptics_no_load_zero_height() {
        let tissue = DeformableTissueHaptics::new(5, 5, 0.01, 10_000.0, 50.0);
        let h = tissue.surface_height(0.02, 0.02);
        assert!(h.abs() < EPS, "no load → zero height, got {h}");
    }

    #[test]
    fn test_tissue_haptics_apply_load_increases_displacement() {
        let mut tissue = DeformableTissueHaptics::new(10, 10, 0.01, 10_000.0, 50.0);
        tissue.apply_point_load(0.04, 0.04, -10.0);
        let h = tissue.surface_height(0.04, 0.04);
        assert!(h < 0.0, "negative load should depress tissue, h={h}");
    }

    #[test]
    fn test_tissue_haptics_step_decays_displacement() {
        let mut tissue = DeformableTissueHaptics::new(10, 10, 0.01, 10_000.0, 10.0);
        tissue.apply_point_load(0.04, 0.04, -10.0);
        let h0 = tissue.surface_height(0.04, 0.04);
        tissue.step(0.1);
        let h1 = tissue.surface_height(0.04, 0.04);
        assert!(
            h1.abs() < h0.abs(),
            "displacement should decay after step: h0={h0}, h1={h1}"
        );
    }

    // -----------------------------------------------------------------------
    // VirtualFixture tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_virtual_fixture_forbidden_outside_region_no_force() {
        let vf = VirtualFixture::new(
            FixtureType::ForbiddenHalfSpace {
                point: [0.0; 3],
                normal: [0.0, 1.0, 0.0],
            },
            1000.0,
            10.0,
        );
        // Position below the plane → outside forbidden region → zero force.
        let f = vf.constraint_force([0.0, -0.1, 0.0], [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_virtual_fixture_forbidden_inside_region_pushes_out() {
        let vf = VirtualFixture::new(
            FixtureType::ForbiddenHalfSpace {
                point: [0.0; 3],
                normal: [0.0, 1.0, 0.0],
            },
            1000.0,
            0.0,
        );
        // Position 5 mm inside forbidden region.
        let f = vf.constraint_force([0.0, 0.005, 0.0], [0.0; 3]);
        // Force should push in -Y direction (away from forbidden region).
        assert!(f[1] < 0.0, "should push out, f_y={}", f[1]);
        assert!((f[1] + 5.0).abs() < EPS, "f_y = -1000 * 0.005 = -5 N");
    }

    #[test]
    fn test_virtual_fixture_linear_guide() {
        let vf = VirtualFixture::new(
            FixtureType::LinearGuide {
                origin: [0.0; 3],
                direction: [1.0, 0.0, 0.0],
            },
            500.0,
            0.0,
        );
        // Position off the X-axis.
        let f = vf.constraint_force([1.0, 0.01, 0.0], [0.0; 3]);
        // Force should pull back towards the X-axis in -Y direction.
        assert!(f[1] < 0.0, "should pull back to guide line, f_y={}", f[1]);
    }

    #[test]
    fn test_virtual_fixture_planar_guide() {
        let vf = VirtualFixture::new(
            FixtureType::PlanarGuide {
                origin: [0.0; 3],
                normal: [0.0, 1.0, 0.0],
            },
            800.0,
            0.0,
        );
        // Position 2 mm above the plane.
        let f = vf.constraint_force([0.0, 0.002, 0.0], [0.0; 3]);
        assert!(f[1] < 0.0, "should push back to plane, f_y={}", f[1]);
        assert!((f[1] + 1.6).abs() < EPS, "f_y = -800 * 0.002 = -1.6 N");
    }

    // -----------------------------------------------------------------------
    // LatencyCompensator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_latency_compensator_constant_velocity() {
        let mut lc = LatencyCompensator::new(0.010); // 10 ms latency
        lc.record(0.000, [0.0, 0.0, 0.0]);
        lc.record(0.001, [0.001, 0.0, 0.0]); // 1 m/s in X
        let pred = lc.predict(0.001);
        // Expected: position at t=0.011 → x = 0.001 + 0.001 * 0.010/0.001 = 0.011
        assert!((pred[0] - 0.011).abs() < 1e-6, "pred_x={}", pred[0]);
    }

    #[test]
    fn test_latency_compensator_single_sample_returns_last() {
        let mut lc = LatencyCompensator::new(0.005);
        lc.record(0.0, [1.0, 2.0, 3.0]);
        let pred = lc.predict(0.0);
        assert_eq!(pred, [1.0, 2.0, 3.0]);
    }

    // -----------------------------------------------------------------------
    // MultiPointHapticContact tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_multipoint_no_active_contacts_zero_force() {
        let mp = MultiPointHapticContact::new(20.0);
        let f = mp.total_force();
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_multipoint_two_active_contacts_sum() {
        let mut mp = MultiPointHapticContact::new(100.0);
        let mut cp1 = ContactPoint::new(0, 1000.0, 0.0);
        cp1.active = true;
        cp1.normal = [0.0, 1.0, 0.0];
        cp1.penetration = 0.002;
        let mut cp2 = ContactPoint::new(1, 1000.0, 0.0);
        cp2.active = true;
        cp2.normal = [0.0, 1.0, 0.0];
        cp2.penetration = 0.003;
        mp.add_contact(cp1);
        mp.add_contact(cp2);
        let f = mp.total_force();
        // 1000 * (0.002 + 0.003) = 5 N in +Y.
        assert!((f[1] - 5.0).abs() < EPS, "f_y={}", f[1]);
    }

    #[test]
    fn test_multipoint_force_clamped() {
        let mut mp = MultiPointHapticContact::new(3.0);
        let mut cp = ContactPoint::new(0, 100_000.0, 0.0);
        cp.active = true;
        cp.normal = [0.0, 1.0, 0.0];
        cp.penetration = 1.0;
        mp.add_contact(cp);
        let f = mp.total_force();
        let mag = norm3(f);
        assert!(
            (mag - 3.0).abs() < 1e-6,
            "should be clamped to 3 N, got {mag}"
        );
    }

    #[test]
    fn test_multipoint_active_count() {
        let mut mp = MultiPointHapticContact::new(10.0);
        let mut cp1 = ContactPoint::new(0, 100.0, 0.0);
        cp1.active = true;
        let cp2 = ContactPoint::new(1, 100.0, 0.0);
        mp.add_contact(cp1);
        mp.add_contact(cp2);
        assert_eq!(mp.active_count(), 1);
    }

    // -----------------------------------------------------------------------
    // HapticFrequencyAnalysis tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_min_stable_frequency() {
        // k=1000 N/m, m=0.05 kg → f_min = sqrt(1000/0.05)/(2π) ≈ 22.5 Hz
        let f = HapticFrequencyAnalysis::min_stable_frequency(1000.0, 0.05);
        assert!((f - 22.508).abs() < 0.01, "f_min={f}");
    }

    #[test]
    fn test_max_renderable_stiffness_at_1khz() {
        // At 1 kHz, m=0.05 kg: k_max = 0.05 * (2π * 1000)^2 ≈ 1.97e6 N/m
        let k = HapticFrequencyAnalysis::max_renderable_stiffness(1000.0, 0.05);
        assert!(k > 1e6, "k_max={k}");
    }

    #[test]
    fn test_required_dt_higher_stiffness_needs_smaller_dt() {
        let dt_soft = HapticFrequencyAnalysis::required_dt(1_000.0, 0.05);
        let dt_hard = HapticFrequencyAnalysis::required_dt(10_000.0, 0.05);
        assert!(
            dt_hard < dt_soft,
            "harder tissue requires smaller dt: dt_soft={dt_soft}, dt_hard={dt_hard}"
        );
    }

    // -----------------------------------------------------------------------
    // TissueLayerModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_tissue_layer_zero_depth_zero_force() {
        let model = TissueLayerModel::default_anatomical();
        let f = model.total_reaction_force(0.0);
        assert!(f.abs() < EPS, "f={f}");
    }

    #[test]
    fn test_tissue_layer_shallow_penetration_skin_only() {
        let model = TissueLayerModel::default_anatomical();
        let f_shallow = model.total_reaction_force(0.001); // 1 mm (< skin 2 mm)
        let (layer_idx, _depth) = model.depth_in_layer(0.001);
        assert_eq!(layer_idx, 0, "should be in skin layer");
        assert!(f_shallow > 0.0);
    }

    #[test]
    fn test_tissue_layer_force_increases_with_depth() {
        let model = TissueLayerModel::default_anatomical();
        let f1 = model.total_reaction_force(0.005);
        let f2 = model.total_reaction_force(0.010);
        assert!(f2 > f1, "deeper penetration should give more force");
    }

    #[test]
    fn test_tissue_layer_depth_in_layer_correct() {
        let model = TissueLayerModel::default_anatomical();
        // Skin is 2 mm thick, fat starts at 2 mm.
        let (layer_idx, _depth) = model.depth_in_layer(0.025); // 25 mm → into muscle
        // skin=2mm, fat=20mm, muscle starts at 22mm → 25mm is in muscle (layer 2)
        assert_eq!(layer_idx, 2, "should be in muscle layer");
    }

    // -----------------------------------------------------------------------
    // HapticSimulationStats tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_haptic_stats_max_force_tracked() {
        let mut stats = HapticSimulationStats::default();
        stats.record_force(1.0, false);
        stats.record_force(5.0, false);
        stats.record_force(3.0, false);
        assert!((stats.max_force_n - 5.0).abs() < EPS);
    }

    #[test]
    fn test_haptic_stats_clamp_events_counted() {
        let mut stats = HapticSimulationStats::default();
        stats.record_force(1.0, false);
        stats.record_force(8.5, true);
        stats.record_force(8.5, true);
        assert_eq!(stats.clamp_events, 2);
    }

    #[test]
    fn test_haptic_stats_avg_force() {
        let mut stats = HapticSimulationStats::default();
        stats.record_force(2.0, false);
        stats.record_force(4.0, false);
        assert!(
            (stats.avg_force_n - 3.0).abs() < 1e-6,
            "avg={}",
            stats.avg_force_n
        );
    }

    // -----------------------------------------------------------------------
    // HapticRenderingLoop integration test
    // -----------------------------------------------------------------------

    #[test]
    fn test_haptic_loop_tick_no_contact_zero_contact_force() {
        let mut loop_ctrl = HapticRenderingLoop::new(1000.0, 10.0, 500.0, 5.0, 0.001);
        // Cursor above surface — no contact.
        let f = loop_ctrl.tick([0.0, 0.01, 0.0], [0.0; 3], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        // Contact force = 0, coupling force = 500 * 0.01 = 5 N (device ahead of tool).
        assert!(norm3(f) > 0.0);
    }

    #[test]
    fn test_haptic_loop_tick_advances_time() {
        let mut loop_ctrl = HapticRenderingLoop::new(1000.0, 10.0, 500.0, 5.0, 0.001);
        loop_ctrl.tick([0.0; 3], [0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]);
        loop_ctrl.tick([0.0; 3], [0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]);
        assert!((loop_ctrl.time - 0.002).abs() < EPS);
        assert_eq!(loop_ctrl.iteration, 2);
    }
}
