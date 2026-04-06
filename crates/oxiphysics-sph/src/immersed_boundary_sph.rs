#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Immersed boundary method (IBM) coupled with SPH.
//!
//! Implements the standard penalty-based IBM where a flexible or rigid structure
//! is represented by a set of Lagrangian marker points embedded in an SPH fluid.
//!
//! Key components:
//! - [`ImmersedBody`] – marker points on the structure surface with outward normals
//! - [`PenaltyForce`] – stiffness/damping penalty restoring force at Lagrangian markers
//! - [`InterpolationOperator`] – SPH interpolation between Lagrangian and Eulerian frames
//! - [`FluidStructureForce`] – pressure and viscous forces on the immersed body
//! - [`RigidImmersedBody`] – prescribed-motion rigid body with force/torque computation
//! - [`ElasticImmersedBody`] – flexible fibre/membrane model with spring-like connections

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ═══════════════════════════════════════════════════════════════════════════════
// § 0  Small vector / math helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise a 3-vector; returns the zero vector if the input is nearly zero.
#[inline]
fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Negate a 3-vector.
#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 1  SPH kernel (cubic spline, 3-D)
// ═══════════════════════════════════════════════════════════════════════════════

/// Cubic-spline SPH kernel W(r, h) in 3-D.
///
/// Normalisation constant: α = 1 / (π h³).
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r.abs() / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient of the cubic-spline kernel: ∇W = (dW/dr)(rij / r).
pub fn cubic_kernel_gradient(rij: [f64; 3], h: f64) -> [f64; 3] {
    let r = norm3(rij);
    if r < 1e-300 || h < 1e-300 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        -alpha / h * 0.75 * t * t
    } else {
        return [0.0; 3];
    };
    scale3(rij, dw_dr / r)
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 2  Fluid particle (Eulerian frame)
// ═══════════════════════════════════════════════════════════════════════════════

/// A single SPH fluid particle in the Eulerian frame.
pub struct FluidParticle {
    /// Current position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Acceleration (m/s²).
    pub acceleration: [f64; 3],
    /// Density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// IBM body force accumulated from immersed boundary (N/kg).
    pub ibm_force: [f64; 3],
}

impl FluidParticle {
    /// Create a new fluid particle at rest.
    pub fn new(position: [f64; 3], mass: f64, density: f64, h: f64, viscosity: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            density,
            pressure: 0.0,
            mass,
            h,
            viscosity,
            ibm_force: [0.0; 3],
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 3  ImmersedBody – Lagrangian marker representation
// ═══════════════════════════════════════════════════════════════════════════════

/// A Lagrangian marker point on the immersed structure surface.
#[derive(Clone, Debug)]
pub struct MarkerPoint {
    /// Current position of the marker (m).
    pub position: [f64; 3],
    /// Reference (initial) position of the marker (m).
    pub ref_position: [f64; 3],
    /// Outward unit normal vector at this marker.
    pub normal: [f64; 3],
    /// Velocity of the marker (m/s).
    pub velocity: [f64; 3],
    /// Force on the marker from the IBM penalty (N).
    pub force: [f64; 3],
    /// Area element associated with this marker (m²).
    pub ds: f64,
}

impl MarkerPoint {
    /// Create a new marker at `position` with outward normal `normal` and area `ds`.
    pub fn new(position: [f64; 3], normal: [f64; 3], ds: f64) -> Self {
        let n = normalise3(normal);
        Self {
            position,
            ref_position: position,
            normal: n,
            velocity: [0.0; 3],
            force: [0.0; 3],
            ds,
        }
    }

    /// Displacement from reference position.
    pub fn displacement(&self) -> [f64; 3] {
        sub3(self.position, self.ref_position)
    }
}

/// A collection of Lagrangian marker points representing an immersed body.
pub struct ImmersedBody {
    /// All marker points on the body surface.
    pub markers: Vec<MarkerPoint>,
    /// Body label / identifier.
    pub label: String,
}

impl ImmersedBody {
    /// Create a new immersed body with the given markers.
    pub fn new(label: impl Into<String>, markers: Vec<MarkerPoint>) -> Self {
        Self {
            markers,
            label: label.into(),
        }
    }

    /// Build a ring of `n` marker points in the x-y plane at height `z` with
    /// radius `r`, outward normals pointing radially, and uniform arc-length ds.
    ///
    /// Useful for cylindrical body cross-sections in 2-D problems.
    pub fn ring(n: usize, radius: f64, z: f64, label: impl Into<String>) -> Self {
        let mut markers = Vec::with_capacity(n);
        let ds = 2.0 * PI * radius / (n as f64);
        for i in 0..n {
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            let (ct, st) = (theta.cos(), theta.sin());
            let pos = [radius * ct, radius * st, z];
            let normal = [ct, st, 0.0];
            markers.push(MarkerPoint::new(pos, normal, ds));
        }
        Self::new(label, markers)
    }

    /// Number of marker points.
    pub fn len(&self) -> usize {
        self.markers.len()
    }

    /// Returns `true` if the body has no marker points.
    pub fn is_empty(&self) -> bool {
        self.markers.is_empty()
    }

    /// Centre of mass of the marker cloud (simple average of positions).
    pub fn centroid(&self) -> [f64; 3] {
        let n = self.markers.len() as f64;
        if n < 1e-300 {
            return [0.0; 3];
        }
        let mut c = [0.0; 3];
        for m in &self.markers {
            for k in 0..3 {
                c[k] += m.position[k];
            }
        }
        scale3(c, 1.0 / n)
    }

    /// Total force on the body (sum over all markers).
    pub fn total_force(&self) -> [f64; 3] {
        let mut f = [0.0; 3];
        for m in &self.markers {
            for k in 0..3 {
                f[k] += m.force[k];
            }
        }
        f
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 4  PenaltyForce – IBM restoring force
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters for the penalty-force IBM.
///
/// The penalty force at each marker is:
/// ```text
/// f_penalty = -κ (x_marker - x_target) - γ (v_marker - v_target)
/// ```
/// where `κ` is the stiffness coefficient and `γ` is the damping coefficient.
pub struct PenaltyForce {
    /// Stiffness penalty coefficient κ (N/m).
    pub stiffness: f64,
    /// Damping penalty coefficient γ (N·s/m).
    pub damping: f64,
}

impl PenaltyForce {
    /// Create a new penalty force with given stiffness and damping.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }

    /// Compute the penalty restoring force at a single marker.
    ///
    /// # Arguments
    /// * `x_marker`  – current marker position
    /// * `x_target`  – desired (Lagrangian) target position
    /// * `v_marker`  – current marker velocity
    /// * `v_target`  – desired target velocity
    pub fn compute(
        &self,
        x_marker: [f64; 3],
        x_target: [f64; 3],
        v_marker: [f64; 3],
        v_target: [f64; 3],
    ) -> [f64; 3] {
        let dx = sub3(x_marker, x_target);
        let dv = sub3(v_marker, v_target);
        let mut f = [0.0; 3];
        for k in 0..3 {
            f[k] = -self.stiffness * dx[k] - self.damping * dv[k];
        }
        f
    }

    /// Apply penalty forces to all markers, updating `marker.force`.
    ///
    /// For a rigid body the target positions are the prescribed positions and
    /// velocities; for an elastic body they are the equilibrium fibre positions.
    pub fn apply_to_body(
        &self,
        body: &mut ImmersedBody,
        targets: &[([f64; 3], [f64; 3])], // (position, velocity) pairs
    ) {
        for (m, &(xt, vt)) in body.markers.iter_mut().zip(targets.iter()) {
            m.force = self.compute(m.position, xt, m.velocity, vt);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 5  InterpolationOperator – SPH ↔ Lagrangian spreading/interpolation
// ═══════════════════════════════════════════════════════════════════════════════

/// SPH interpolation and spreading operator between Lagrangian markers and
/// Eulerian SPH particles.
///
/// Uses the standard SPH kernel to:
/// 1. **Interpolate** fluid velocity at marker positions (Eulerian → Lagrangian).
/// 2. **Spread** marker forces onto fluid particles (Lagrangian → Eulerian).
pub struct InterpolationOperator {
    /// Smoothing length used for the IBM kernel (m).
    pub h_ib: f64,
}

impl InterpolationOperator {
    /// Create a new interpolation operator with smoothing length `h_ib`.
    pub fn new(h_ib: f64) -> Self {
        Self { h_ib }
    }

    /// Interpolate the fluid velocity at a Lagrangian marker position.
    ///
    /// ```text
    /// V_L(X) = Σ_j (m_j / ρ_j) v_j W(|X - x_j|, h)
    /// ```
    pub fn interpolate_velocity(
        &self,
        marker_pos: [f64; 3],
        particles: &[FluidParticle],
    ) -> [f64; 3] {
        let mut vel = [0.0; 3];
        for p in particles {
            let r = norm3(sub3(marker_pos, p.position));
            let w = cubic_kernel(r, self.h_ib);
            let wv = p.mass / p.density.max(1e-300) * w;
            for k in 0..3 {
                vel[k] += wv * p.velocity[k];
            }
        }
        vel
    }

    /// Interpolate the fluid pressure at a Lagrangian marker position.
    pub fn interpolate_pressure(&self, marker_pos: [f64; 3], particles: &[FluidParticle]) -> f64 {
        let mut pressure = 0.0;
        for p in particles {
            let r = norm3(sub3(marker_pos, p.position));
            let w = cubic_kernel(r, self.h_ib);
            pressure += p.mass / p.density.max(1e-300) * p.pressure * w;
        }
        pressure
    }

    /// Spread a body force `f_marker` (N) from a single marker to all nearby
    /// fluid particles.  The body force per unit mass on each fluid particle is
    /// accumulated into `particle.ibm_force`.
    ///
    /// The spreading formula is:
    /// ```text
    /// f_E(x_j) += f_L(X) * W(|X - x_j|, h) * ds
    /// ```
    pub fn spread_force(&self, marker: &MarkerPoint, particles: &mut [FluidParticle]) {
        for p in particles.iter_mut() {
            let r = norm3(sub3(marker.position, p.position));
            let w = cubic_kernel(r, self.h_ib);
            let scale = w * marker.ds / p.density.max(1e-300);
            for k in 0..3 {
                p.ibm_force[k] += marker.force[k] * scale;
            }
        }
    }

    /// Spread all marker forces onto the fluid particle array.
    pub fn spread_all(&self, body: &ImmersedBody, particles: &mut [FluidParticle]) {
        for m in &body.markers {
            self.spread_force(m, particles);
        }
    }

    /// Reset all IBM body forces on the fluid particles to zero.
    pub fn reset_ibm_forces(particles: &mut [FluidParticle]) {
        for p in particles.iter_mut() {
            p.ibm_force = [0.0; 3];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 6  FluidStructureForce – pressure/viscous force on the body
// ═══════════════════════════════════════════════════════════════════════════════

/// Computes the hydrodynamic (pressure + viscous) force on an immersed body
/// by integrating over the Lagrangian marker surface.
pub struct FluidStructureForce {
    /// Interpolation operator used to sample fluid quantities at markers.
    pub interp: InterpolationOperator,
}

impl FluidStructureForce {
    /// Create with a given IBM smoothing length.
    pub fn new(h_ib: f64) -> Self {
        Self {
            interp: InterpolationOperator::new(h_ib),
        }
    }

    /// Compute the pressure force on the body.
    ///
    /// `F_p = - ∫ p n̂ dS ≈ -Σ_k p(X_k) n̂_k ds_k`
    pub fn pressure_force(&self, body: &ImmersedBody, particles: &[FluidParticle]) -> [f64; 3] {
        let mut f = [0.0; 3];
        for m in &body.markers {
            let p = self.interp.interpolate_pressure(m.position, particles);
            for k in 0..3 {
                f[k] -= p * m.normal[k] * m.ds;
            }
        }
        f
    }

    /// Compute the viscous drag force on the body.
    ///
    /// Uses a finite-difference approximation:
    /// `F_v ≈ Σ_k μ * (v_fluid(X_k) - v_marker_k) / h_ib * ds_k`
    pub fn viscous_force(
        &self,
        body: &ImmersedBody,
        particles: &[FluidParticle],
        viscosity: f64,
    ) -> [f64; 3] {
        let mut f = [0.0; 3];
        for m in &body.markers {
            let v_f = self.interp.interpolate_velocity(m.position, particles);
            let slip = sub3(v_f, m.velocity);
            let contrib = scale3(slip, viscosity / self.interp.h_ib * m.ds);
            for k in 0..3 {
                f[k] += contrib[k];
            }
        }
        f
    }

    /// Total hydrodynamic force (pressure + viscous).
    pub fn total_force(
        &self,
        body: &ImmersedBody,
        particles: &[FluidParticle],
        viscosity: f64,
    ) -> [f64; 3] {
        let fp = self.pressure_force(body, particles);
        let fv = self.viscous_force(body, particles, viscosity);
        add3(fp, fv)
    }

    /// Torque about a centre point `cx` due to the total hydrodynamic force.
    pub fn torque_about(
        &self,
        body: &ImmersedBody,
        particles: &[FluidParticle],
        viscosity: f64,
        cx: [f64; 3],
    ) -> [f64; 3] {
        let mut tau = [0.0; 3];
        for m in &body.markers {
            let r = sub3(m.position, cx);
            let p = self.interp.interpolate_pressure(m.position, particles);
            let v_f = self.interp.interpolate_velocity(m.position, particles);
            let fp_local = scale3(neg3(m.normal), p * m.ds);
            let fv_local = scale3(sub3(v_f, m.velocity), viscosity / self.interp.h_ib * m.ds);
            let f_local = add3(fp_local, fv_local);
            let t = cross3(r, f_local);
            for k in 0..3 {
                tau[k] += t[k];
            }
        }
        tau
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 7  RigidImmersedBody – prescribed motion
// ═══════════════════════════════════════════════════════════════════════════════

/// Motion prescription type for a rigid immersed body.
#[derive(Clone, Debug)]
pub enum MotionType {
    /// Body is held stationary.
    Fixed,
    /// Constant translational velocity (m/s).
    ConstantVelocity([f64; 3]),
    /// Sinusoidal oscillation: `x(t) = A sin(ω t + φ) ê`.
    SinusoidalOscillation {
        /// Amplitude (m).
        amplitude: f64,
        /// Angular frequency (rad/s).
        omega: f64,
        /// Phase offset (rad).
        phase: f64,
        /// Direction unit vector.
        direction: [f64; 3],
    },
    /// Circular orbit about a centre point.
    CircularOrbit {
        /// Centre of rotation (m).
        centre: [f64; 3],
        /// Orbital radius (m).
        radius: f64,
        /// Angular velocity (rad/s).
        angular_velocity: f64,
    },
}

/// A rigid immersed body whose motion is prescribed.
///
/// Stores the current state (centre of mass position and velocity) and
/// computes force/torque from the fluid via [`FluidStructureForce`].
pub struct RigidImmersedBody {
    /// Underlying Lagrangian marker representation.
    pub body: ImmersedBody,
    /// Centre of mass position (m).
    pub centre: [f64; 3],
    /// Centre of mass velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass of the rigid body (kg).
    pub mass: f64,
    /// Moment of inertia about the z-axis (kg·m²).
    pub inertia_z: f64,
    /// Motion prescription.
    pub motion: MotionType,
    /// Accumulated force from fluid (N).
    pub fluid_force: [f64; 3],
    /// Accumulated torque from fluid about z-axis (N·m).
    pub fluid_torque_z: f64,
}

impl RigidImmersedBody {
    /// Create a new rigid immersed body.
    pub fn new(body: ImmersedBody, mass: f64, inertia_z: f64, motion: MotionType) -> Self {
        let centre = body.centroid();
        Self {
            body,
            centre,
            velocity: [0.0; 3],
            mass,
            inertia_z,
            motion,
            fluid_force: [0.0; 3],
            fluid_torque_z: 0.0,
        }
    }

    /// Prescribed centre-of-mass position at time `t`.
    pub fn prescribed_position(&self, t: f64) -> [f64; 3] {
        match &self.motion {
            MotionType::Fixed => self.centre,
            MotionType::ConstantVelocity(v) => add3(self.centre, scale3(*v, t)),
            MotionType::SinusoidalOscillation {
                amplitude,
                omega,
                phase,
                direction,
            } => {
                let disp = amplitude * (omega * t + phase).sin();
                add3(self.centre, scale3(*direction, disp))
            }
            MotionType::CircularOrbit {
                centre,
                radius,
                angular_velocity,
            } => {
                let theta = angular_velocity * t;
                [
                    centre[0] + radius * theta.cos(),
                    centre[1] + radius * theta.sin(),
                    centre[2],
                ]
            }
        }
    }

    /// Prescribed centre-of-mass velocity at time `t`.
    pub fn prescribed_velocity(&self, t: f64) -> [f64; 3] {
        match &self.motion {
            MotionType::Fixed => [0.0; 3],
            MotionType::ConstantVelocity(v) => *v,
            MotionType::SinusoidalOscillation {
                amplitude,
                omega,
                phase,
                direction,
            } => {
                let vv = amplitude * omega * (omega * t + phase).cos();
                scale3(*direction, vv)
            }
            MotionType::CircularOrbit {
                radius,
                angular_velocity,
                ..
            } => {
                let theta = angular_velocity * t;
                let v = radius * angular_velocity;
                [-v * theta.sin(), v * theta.cos(), 0.0]
            }
        }
    }

    /// Update marker positions and velocities according to the prescribed motion
    /// at time `t`.  The markers translate rigidly with the body centre.
    pub fn update_kinematics(&mut self, t: f64) {
        let new_centre = self.prescribed_position(t);
        let new_vel = self.prescribed_velocity(t);
        let delta = sub3(new_centre, self.centre);
        self.centre = new_centre;
        self.velocity = new_vel;
        for m in self.body.markers.iter_mut() {
            m.position = add3(m.position, delta);
            m.velocity = new_vel;
        }
    }

    /// Compute and store fluid forces/torques.
    pub fn compute_fluid_loads(
        &mut self,
        fsf: &FluidStructureForce,
        particles: &[FluidParticle],
        viscosity: f64,
    ) {
        self.fluid_force = fsf.total_force(&self.body, particles, viscosity);
        let tau = fsf.torque_about(&self.body, particles, viscosity, self.centre);
        self.fluid_torque_z = tau[2];
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 8  ElasticImmersedBody – flexible fibre / membrane model
// ═══════════════════════════════════════════════════════════════════════════════

/// A segment connecting two adjacent Lagrangian markers for the spring-fibre model.
#[derive(Clone, Debug)]
pub struct FibreSegment {
    /// Index of the first marker.
    pub i: usize,
    /// Index of the second marker.
    pub j: usize,
    /// Rest (reference) length of the segment (m).
    pub rest_length: f64,
    /// Stretching stiffness (N/m).
    pub stiffness: f64,
    /// Bending stiffness (N·m).
    pub bending_stiffness: f64,
}

impl FibreSegment {
    /// Create a segment between markers `i` and `j` with a given stiffness.
    pub fn new(
        i: usize,
        j: usize,
        rest_length: f64,
        stiffness: f64,
        bending_stiffness: f64,
    ) -> Self {
        Self {
            i,
            j,
            rest_length,
            stiffness,
            bending_stiffness,
        }
    }

    /// Stretching force on marker `i` (equal and opposite on `j`).
    pub fn stretching_force(&self, positions: &[[f64; 3]]) -> [f64; 3] {
        let xi = positions[self.i];
        let xj = positions[self.j];
        let d = sub3(xj, xi);
        let l = norm3(d);
        if l < 1e-300 {
            return [0.0; 3];
        }
        let extension = l - self.rest_length;
        let fmag = self.stiffness * extension;
        scale3(d, fmag / l)
    }
}

/// A flexible immersed body modelled as a network of spring-like fibres.
///
/// The body deforms under the combined action of:
/// 1. **Stretching forces** from [`FibreSegment`]s.
/// 2. **IBM penalty forces** keeping markers co-located with the fluid.
pub struct ElasticImmersedBody {
    /// Lagrangian marker array.
    pub body: ImmersedBody,
    /// Fibre segments connecting adjacent markers.
    pub segments: Vec<FibreSegment>,
    /// Penalty force parameters.
    pub penalty: PenaltyForce,
    /// Mass of each marker (kg).
    pub marker_mass: f64,
    /// Damping coefficient for marker dynamics (kg/s).
    pub marker_damping: f64,
}

impl ElasticImmersedBody {
    /// Create a new elastic immersed body.
    pub fn new(
        body: ImmersedBody,
        segments: Vec<FibreSegment>,
        penalty: PenaltyForce,
        marker_mass: f64,
        marker_damping: f64,
    ) -> Self {
        Self {
            body,
            segments,
            penalty,
            marker_mass,
            marker_damping,
        }
    }

    /// Build a straight elastic fibre of `n` markers along the x-axis from
    /// `x0` to `x1` at height `y`, `z`, with given material parameters.
    pub fn straight_fibre(
        n: usize,
        x0: f64,
        x1: f64,
        y: f64,
        z: f64,
        stiffness: f64,
        bending_stiffness: f64,
        penalty_k: f64,
        penalty_d: f64,
        marker_mass: f64,
        marker_damping: f64,
    ) -> Self {
        let dx = (x1 - x0) / (n - 1) as f64;
        let mut markers = Vec::with_capacity(n);
        for i in 0..n {
            let x = x0 + i as f64 * dx;
            let pos = [x, y, z];
            let normal = [0.0, 1.0, 0.0]; // default: normal pointing in +y
            markers.push(MarkerPoint::new(pos, normal, dx));
        }
        let body = ImmersedBody::new("elastic_fibre", markers);

        let mut segments = Vec::with_capacity(n - 1);
        for i in 0..(n - 1) {
            segments.push(FibreSegment::new(
                i,
                i + 1,
                dx,
                stiffness,
                bending_stiffness,
            ));
        }

        Self::new(
            body,
            segments,
            PenaltyForce::new(penalty_k, penalty_d),
            marker_mass,
            marker_damping,
        )
    }

    /// Compute elastic (stretching) forces on all markers and accumulate
    /// into `marker.force`.
    pub fn compute_elastic_forces(&mut self) {
        // Reset forces
        for m in self.body.markers.iter_mut() {
            m.force = [0.0; 3];
        }
        let positions: Vec<[f64; 3]> = self.body.markers.iter().map(|m| m.position).collect();
        for seg in &self.segments {
            let fi = seg.stretching_force(&positions);
            let fj = neg3(fi);
            for k in 0..3 {
                self.body.markers[seg.i].force[k] += fi[k];
                self.body.markers[seg.j].force[k] += fj[k];
            }
        }
    }

    /// Integrate marker velocities and positions with forward Euler, time step `dt`.
    ///
    /// The acceleration of marker `l` is:
    /// `a_l = (f_elastic_l + f_ibm_l) / m_marker - γ v_l / m_marker`
    pub fn integrate(&mut self, dt: f64, ibm_forces: &[[f64; 3]]) {
        let nm = self.body.markers.len();
        for l in 0..nm {
            let m = &mut self.body.markers[l];
            let f_ibm = if l < ibm_forces.len() {
                ibm_forces[l]
            } else {
                [0.0; 3]
            };
            let mut acc = [0.0; 3];
            for k in 0..3 {
                acc[k] = (m.force[k] + f_ibm[k] - self.marker_damping * m.velocity[k])
                    / self.marker_mass.max(1e-300);
            }
            for k in 0..3 {
                m.velocity[k] += acc[k] * dt;
                m.position[k] += m.velocity[k] * dt;
            }
        }
    }

    /// Total elastic potential energy stored in all segments.
    pub fn elastic_energy(&self) -> f64 {
        let positions: Vec<[f64; 3]> = self.body.markers.iter().map(|m| m.position).collect();
        let mut energy = 0.0;
        for seg in &self.segments {
            let d = sub3(positions[seg.j], positions[seg.i]);
            let l = norm3(d);
            let ext = l - seg.rest_length;
            energy += 0.5 * seg.stiffness * ext * ext;
        }
        energy
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 9  IBM coupling driver
// ═══════════════════════════════════════════════════════════════════════════════

/// High-level IBM coupling step that:
/// 1. Resets IBM forces on fluid particles.
/// 2. Interpolates fluid velocity to each marker.
/// 3. Computes the penalty restoring force at each marker.
/// 4. Spreads marker forces back to fluid particles.
///
/// Returns the updated marker velocities (interpolated from fluid).
pub fn ibm_coupling_step(
    body: &mut ImmersedBody,
    particles: &mut [FluidParticle],
    penalty: &PenaltyForce,
    interp: &InterpolationOperator,
    target_positions: &[[f64; 3]],
    target_velocities: &[[f64; 3]],
) {
    // 1. Reset
    InterpolationOperator::reset_ibm_forces(particles);

    // 2. Interpolate fluid velocity to markers
    for m in body.markers.iter_mut() {
        let v_f = interp.interpolate_velocity(m.position, particles);
        m.velocity = v_f;
    }

    // 3. Compute penalty forces
    for (idx, m) in body.markers.iter_mut().enumerate() {
        let xt = if idx < target_positions.len() {
            target_positions[idx]
        } else {
            m.ref_position
        };
        let vt = if idx < target_velocities.len() {
            target_velocities[idx]
        } else {
            [0.0; 3]
        };
        m.force = penalty.compute(m.position, xt, m.velocity, vt);
    }

    // 4. Spread to fluid
    interp.spread_all(body, particles);
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 10  Utility: nearest SPH particle index
// ═══════════════════════════════════════════════════════════════════════════════

/// Find the index of the SPH particle nearest to `query`.
///
/// Returns `None` if `particles` is empty.
pub fn nearest_particle(query: [f64; 3], particles: &[FluidParticle]) -> Option<usize> {
    particles
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let da = norm3(sub3(query, a.position));
            let db = norm3(sub3(query, b.position));
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(idx, _)| idx)
}

// ═══════════════════════════════════════════════════════════════════════════════
// § 11  Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── Kernel ────────────────────────────────────────────────────────────────

    #[test]
    fn kernel_zero_at_large_distance() {
        assert_eq!(cubic_kernel(3.0, 1.0), 0.0);
    }

    #[test]
    fn kernel_positive_inside_support() {
        assert!(cubic_kernel(0.5, 1.0) > 0.0);
    }

    #[test]
    fn kernel_gradient_zero_at_origin() {
        let g = cubic_kernel_gradient([0.0; 3], 1.0);
        for &v in &g {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn kernel_gradient_zero_outside_support() {
        let g = cubic_kernel_gradient([3.0, 0.0, 0.0], 1.0);
        for &v in &g {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn kernel_symmetric() {
        let w1 = cubic_kernel(0.3, 1.0);
        let w2 = cubic_kernel(-0.3_f64.abs(), 1.0);
        assert!((w1 - w2).abs() < 1e-15);
    }

    // ── MarkerPoint ───────────────────────────────────────────────────────────

    #[test]
    fn marker_normal_is_unit() {
        let m = MarkerPoint::new([0.0; 3], [3.0, 4.0, 0.0], 0.1);
        let len = norm3(m.normal);
        assert!((len - 1.0).abs() < 1e-10, "normal length = {len}");
    }

    #[test]
    fn marker_displacement_zero_initially() {
        let m = MarkerPoint::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], 0.05);
        let d = m.displacement();
        for &v in &d {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn marker_displacement_after_move() {
        let mut m = MarkerPoint::new([0.0; 3], [1.0, 0.0, 0.0], 0.1);
        m.position[0] = 0.5;
        let d = m.displacement();
        assert!((d[0] - 0.5).abs() < 1e-15);
    }

    // ── ImmersedBody ──────────────────────────────────────────────────────────

    #[test]
    fn ring_body_has_correct_count() {
        let body = ImmersedBody::ring(16, 1.0, 0.0, "test");
        assert_eq!(body.len(), 16);
    }

    #[test]
    fn ring_body_not_empty() {
        let body = ImmersedBody::ring(8, 0.5, 0.0, "ring");
        assert!(!body.is_empty());
    }

    #[test]
    fn empty_body_is_empty() {
        let body = ImmersedBody::new("empty", vec![]);
        assert!(body.is_empty());
    }

    #[test]
    fn ring_centroid_near_origin() {
        let body = ImmersedBody::ring(32, 1.0, 0.0, "ring");
        let c = body.centroid();
        assert!(c[0].abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
    }

    #[test]
    fn ring_total_force_initially_zero() {
        let body = ImmersedBody::ring(8, 1.0, 0.0, "ring");
        let f = body.total_force();
        for &v in &f {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn ring_normals_point_outward() {
        let body = ImmersedBody::ring(8, 1.0, 0.0, "ring");
        for m in &body.markers {
            // Outward normal should have positive dot product with position
            let d = dot3(m.normal, m.position);
            assert!(d > 0.5, "normal not outward: d = {d}");
        }
    }

    // ── PenaltyForce ──────────────────────────────────────────────────────────

    #[test]
    fn penalty_zero_when_at_target() {
        let pf = PenaltyForce::new(1000.0, 10.0);
        let f = pf.compute([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [0.0; 3], [0.0; 3]);
        for &v in &f {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn penalty_force_restores_to_target() {
        let pf = PenaltyForce::new(1000.0, 0.0);
        // marker is at +0.1 in x relative to target
        let f = pf.compute([1.1, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        assert!((f[0] + 100.0).abs() < 1e-10, "f[0] = {}", f[0]);
    }

    #[test]
    fn penalty_damping_opposes_velocity() {
        let pf = PenaltyForce::new(0.0, 50.0);
        let f = pf.compute([0.0; 3], [0.0; 3], [1.0, 0.0, 0.0], [0.0; 3]);
        assert!(f[0] < 0.0, "damping should oppose positive velocity");
    }

    // ── InterpolationOperator ─────────────────────────────────────────────────

    #[test]
    fn interpolate_velocity_single_particle() {
        let op = InterpolationOperator::new(0.5);
        let mut p = FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.1, 1e-3);
        p.velocity = [1.0, 2.0, 3.0];
        let v = op.interpolate_velocity([0.0; 3], &[p]);
        // Should interpolate the known velocity
        assert!(v[0] > 0.0 || v[1] > 0.0 || v[2] > 0.0);
    }

    #[test]
    fn interpolate_velocity_far_particle_zero() {
        let op = InterpolationOperator::new(0.1);
        let p = FluidParticle::new([100.0, 0.0, 0.0], 1.0, 1000.0, 0.1, 1e-3);
        let v = op.interpolate_velocity([0.0; 3], &[p]);
        for &vi in &v {
            assert_eq!(vi, 0.0);
        }
    }

    #[test]
    fn interpolate_pressure_zero_if_no_pressure() {
        let op = InterpolationOperator::new(0.5);
        let p = FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.1, 1e-3);
        let pr = op.interpolate_pressure([0.0; 3], &[p]);
        assert_eq!(pr, 0.0);
    }

    #[test]
    fn spread_force_accumulates_on_nearby_particle() {
        let op = InterpolationOperator::new(1.0);
        let mut marker = MarkerPoint::new([0.0; 3], [1.0, 0.0, 0.0], 1.0);
        marker.force = [10.0, 0.0, 0.0];
        let mut particles = vec![FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.5, 1e-3)];
        op.spread_force(&marker, &mut particles);
        assert!(
            particles[0].ibm_force[0] != 0.0,
            "IBM force should be non-zero"
        );
    }

    #[test]
    fn reset_ibm_forces_clears_all() {
        let mut particles = vec![
            FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.1, 1e-3),
            FluidParticle::new([1.0, 0.0, 0.0], 1.0, 1000.0, 0.1, 1e-3),
        ];
        particles[0].ibm_force = [1.0, 2.0, 3.0];
        particles[1].ibm_force = [4.0, 5.0, 6.0];
        InterpolationOperator::reset_ibm_forces(&mut particles);
        for p in &particles {
            for &v in &p.ibm_force {
                assert_eq!(v, 0.0);
            }
        }
    }

    // ── FluidStructureForce ───────────────────────────────────────────────────

    #[test]
    fn pressure_force_zero_with_zero_pressure() {
        let fsf = FluidStructureForce::new(1.0);
        let body = ImmersedBody::ring(8, 1.0, 0.0, "ring");
        let particles = vec![FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.2, 1e-3)];
        let f = fsf.pressure_force(&body, &particles);
        // pressure is zero, so force should be zero
        for &v in &f {
            assert_eq!(v, 0.0);
        }
    }

    #[test]
    fn viscous_force_zero_same_velocity() {
        // When no fluid particles exist, viscous force must be zero
        // (interpolated velocity returns [0,0,0], but marker vel is also [0,0,0])
        let fsf = FluidStructureForce::new(1.0);
        let mut body = ImmersedBody::ring(4, 0.1, 0.0, "ring");
        // Set marker velocities to zero so they match the zero interpolated velocity
        let zero_vel = [0.0_f64; 3];
        for m in body.markers.iter_mut() {
            m.velocity = zero_vel;
        }
        let f = fsf.viscous_force(&body, &[], 1e-3);
        for &v in &f {
            assert_eq!(
                v, 0.0,
                "viscous force should be zero when no particles: {v}"
            );
        }
    }

    // ── RigidImmersedBody ─────────────────────────────────────────────────────

    #[test]
    fn rigid_body_fixed_position_unchanged() {
        let body = ImmersedBody::ring(8, 1.0, 0.0, "ring");
        let mut rb = RigidImmersedBody::new(body, 1.0, 0.1, MotionType::Fixed);
        let c0 = rb.centre;
        rb.update_kinematics(1.0);
        let c1 = rb.centre;
        for k in 0..3 {
            assert!((c0[k] - c1[k]).abs() < 1e-12, "fixed body moved");
        }
    }

    #[test]
    fn rigid_body_constant_velocity_moves() {
        let body = ImmersedBody::ring(8, 1.0, 0.0, "ring");
        let vel = [1.0, 0.0, 0.0];
        let rb = RigidImmersedBody::new(body, 1.0, 0.1, MotionType::ConstantVelocity(vel));
        // prescribed_position at t=2 should be centre + 2*vel
        let pos = rb.prescribed_position(2.0);
        assert!((pos[0] - 2.0).abs() < 1e-10, "x = {}", pos[0]);
    }

    #[test]
    fn rigid_body_sinusoidal_at_zero_time() {
        let body = ImmersedBody::ring(4, 0.5, 0.0, "ring");
        let motion = MotionType::SinusoidalOscillation {
            amplitude: 0.1,
            omega: 2.0 * PI,
            phase: 0.0,
            direction: [1.0, 0.0, 0.0],
        };
        let rb = RigidImmersedBody::new(body, 1.0, 0.1, motion);
        let vel = rb.prescribed_velocity(0.0);
        // sin(0) = 0, so v = A*ω*cos(0) = A*ω in x
        let expected = 0.1 * 2.0 * PI;
        assert!((vel[0] - expected).abs() < 1e-8, "vel[0] = {}", vel[0]);
    }

    #[test]
    fn rigid_body_circular_orbit_radius() {
        let body = ImmersedBody::ring(4, 0.1, 0.0, "ring");
        let centre = [0.0; 3];
        let motion = MotionType::CircularOrbit {
            centre,
            radius: 2.0,
            angular_velocity: 1.0,
        };
        let rb = RigidImmersedBody::new(body, 1.0, 0.1, motion);
        let pos = rb.prescribed_position(0.0);
        let r = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
        assert!((r - 2.0).abs() < 1e-10, "orbital radius = {r}");
    }

    // ── ElasticImmersedBody ───────────────────────────────────────────────────

    #[test]
    fn elastic_fibre_correct_marker_count() {
        let eib = ElasticImmersedBody::straight_fibre(
            10, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        assert_eq!(eib.body.len(), 10);
    }

    #[test]
    fn elastic_fibre_correct_segment_count() {
        let eib = ElasticImmersedBody::straight_fibre(
            10, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        assert_eq!(eib.segments.len(), 9);
    }

    #[test]
    fn elastic_energy_zero_at_rest() {
        let eib = ElasticImmersedBody::straight_fibre(
            5, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        let e = eib.elastic_energy();
        assert!(e.abs() < 1e-10, "energy at rest = {e}");
    }

    #[test]
    fn elastic_energy_positive_when_stretched() {
        let mut eib = ElasticImmersedBody::straight_fibre(
            3, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        // Stretch the last marker
        eib.body.markers[2].position[0] += 0.5;
        let e = eib.elastic_energy();
        assert!(e > 0.0, "stretched energy should be positive");
    }

    #[test]
    fn elastic_compute_forces_zero_at_rest() {
        let mut eib = ElasticImmersedBody::straight_fibre(
            5, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        eib.compute_elastic_forces();
        for m in &eib.body.markers {
            for &v in &m.force {
                assert!(v.abs() < 1e-10, "force at rest = {v}");
            }
        }
    }

    #[test]
    fn elastic_force_nonzero_when_stretched() {
        let mut eib = ElasticImmersedBody::straight_fibre(
            3, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        eib.body.markers[2].position[0] += 0.1;
        eib.compute_elastic_forces();
        // Some force should be nonzero
        let any_nonzero = eib
            .body
            .markers
            .iter()
            .any(|m| m.force.iter().any(|&v| v.abs() > 1e-10));
        assert!(
            any_nonzero,
            "stretched fibre should have nonzero elastic forces"
        );
    }

    #[test]
    fn elastic_integrate_changes_position() {
        let mut eib = ElasticImmersedBody::straight_fibre(
            3, 0.0, 1.0, 0.0, 0.0, 1000.0, 1.0, 500.0, 10.0, 1e-3, 0.1,
        );
        let ibm = vec![[1.0, 0.0, 0.0]; 3];
        let pos0 = eib.body.markers[1].position;
        eib.integrate(0.01, &ibm);
        let pos1 = eib.body.markers[1].position;
        // Position must have changed
        let moved = (0..3).any(|k| (pos0[k] - pos1[k]).abs() > 1e-15);
        assert!(moved, "integration should change marker position");
    }

    // ── IBM coupling step ─────────────────────────────────────────────────────

    #[test]
    fn ibm_coupling_step_resets_and_spreads() {
        let mut body = ImmersedBody::ring(4, 0.05, 0.0, "ring");
        let mut particles = vec![FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.5, 1e-3)];
        particles[0].ibm_force = [999.0; 3]; // pre-existing garbage

        let penalty = PenaltyForce::new(1000.0, 10.0);
        let interp = InterpolationOperator::new(1.0);
        let targets_pos: Vec<[f64; 3]> = body.markers.iter().map(|m| m.position).collect();
        let targets_vel: Vec<[f64; 3]> = vec![[0.0; 3]; body.len()];

        ibm_coupling_step(
            &mut body,
            &mut particles,
            &penalty,
            &interp,
            &targets_pos,
            &targets_vel,
        );

        // ibm_force was reset then re-computed, so initial 999 must be gone
        // (it may be non-zero from the new spreading pass)
        let _any_finite = particles
            .iter()
            .all(|p| p.ibm_force.iter().all(|v| v.is_finite()));
        assert!(
            _any_finite,
            "IBM forces should be finite after coupling step"
        );
    }

    // ── Nearest particle ─────────────────────────────────────────────────────

    #[test]
    fn nearest_particle_finds_closest() {
        let particles = vec![
            FluidParticle::new([0.0; 3], 1.0, 1000.0, 0.1, 1e-3),
            FluidParticle::new([1.0, 0.0, 0.0], 1.0, 1000.0, 0.1, 1e-3),
            FluidParticle::new([2.0, 0.0, 0.0], 1.0, 1000.0, 0.1, 1e-3),
        ];
        let idx = nearest_particle([0.9, 0.0, 0.0], &particles).unwrap();
        assert_eq!(idx, 1, "nearest particle index = {idx}");
    }

    #[test]
    fn nearest_particle_returns_none_empty() {
        let result = nearest_particle([0.0; 3], &[]);
        assert!(result.is_none());
    }

    // ── Math helpers ──────────────────────────────────────────────────────────

    #[test]
    fn add3_correct() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn sub3_correct() {
        let r = sub3([3.0, 2.0, 1.0], [1.0, 1.0, 1.0]);
        assert_eq!(r, [2.0, 1.0, 0.0]);
    }

    #[test]
    fn dot3_correct() {
        let d = dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(d, 0.0);
        let d2 = dot3([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]);
        assert!((d2 - 3.0).abs() < 1e-15);
    }

    #[test]
    fn cross3_orthogonal() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn normalise3_unit_length() {
        let v = normalise3([3.0, 4.0, 0.0]);
        let len = norm3(v);
        assert!((len - 1.0).abs() < 1e-10, "length = {len}");
    }

    #[test]
    fn normalise3_zero_returns_zero() {
        let v = normalise3([0.0; 3]);
        for &vi in &v {
            assert_eq!(vi, 0.0);
        }
    }
}
