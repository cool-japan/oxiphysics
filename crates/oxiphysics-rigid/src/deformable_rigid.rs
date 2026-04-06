#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Deformable-rigid body coupling: elastic shells, contact compliance,
//! rigid-flexi force/torque coupling, impact deformation, soft contact
//! patches, and 6-D wrench space.
//!
//! # Overview
//!
//! - [`DeformableBody`] — elastic shell with rigid core, compliance matrix,
//!   and modal deformation representation.
//! - [`ContactCompliance`] — Hertz contact model, compliant contact, and
//!   penetration depth utilities.
//! - [`RigidFlexiCoupling`] — maps deformation modal amplitudes to resultant
//!   force/torque on the rigid body via modal superposition.
//! - [`ImpactDeformation`] — energy absorption, permanent plastic deformation,
//!   and restitution coefficient computation.
//! - [`SoftContactPatch`] — contact area, pressure distribution, and traction
//!   integration over an elliptical contact patch.
//! - [`WrenchSpace`] — 6-D wrench representation, grasp map, and contact
//!   wrench cone.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper math (no nalgebra in this crate)
// ---------------------------------------------------------------------------

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn v3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n < 1e-12 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

// ---------------------------------------------------------------------------
// DeformableBody
// ---------------------------------------------------------------------------

/// Number of retained vibration modes in the reduced-order model.
pub const NUM_MODES: usize = 6;

/// Elastic shell body coupled to a rigid core.
///
/// The deformation field is represented via a truncated modal basis.
/// `q[i]` is the generalised coordinate (amplitude) of mode `i`.
/// `q_dot[i]` is the corresponding modal velocity.
#[derive(Debug, Clone)]
pub struct DeformableBody {
    /// Mass of the rigid core (kg).
    pub core_mass: f64,
    /// Shell thickness (m).
    pub shell_thickness: f64,
    /// Young's modulus of the shell material (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio of the shell material (dimensionless).
    pub poissons_ratio: f64,
    /// Shell mid-surface area (m²).
    pub shell_area: f64,
    /// Modal stiffness coefficients (N/m), one per mode.
    pub modal_stiffness: [f64; NUM_MODES],
    /// Modal damping coefficients (N·s/m), one per mode.
    pub modal_damping: [f64; NUM_MODES],
    /// Modal mass coefficients (kg), one per mode.
    pub modal_mass: [f64; NUM_MODES],
    /// Current modal amplitudes (m).
    pub q: [f64; NUM_MODES],
    /// Current modal velocities (m/s).
    pub q_dot: [f64; NUM_MODES],
    /// 6×6 compliance matrix (flexibility, m/N for translational entries).
    ///
    /// Stored in row-major order as a flat array of 36 entries.
    pub compliance: [f64; 36],
}

impl DeformableBody {
    /// Construct a deformable body with given core mass and shell parameters.
    ///
    /// Modal properties are initialised using a simple beam analogy.
    pub fn new(
        core_mass: f64,
        shell_thickness: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
        shell_area: f64,
    ) -> Self {
        // Simple modal property estimate: k_i ∝ E * h³ / (A * i²)
        let h = shell_thickness;
        let bending_stiffness =
            youngs_modulus * h * h * h / (12.0 * (1.0 - poissons_ratio * poissons_ratio));
        let mut modal_stiffness = [0.0_f64; NUM_MODES];
        let mut modal_damping = [0.0_f64; NUM_MODES];
        let mut modal_mass = [0.0_f64; NUM_MODES];
        for i in 0..NUM_MODES {
            let n = (i + 1) as f64;
            modal_stiffness[i] = bending_stiffness * shell_area * n * n;
            modal_mass[i] = core_mass / (NUM_MODES as f64);
            // 5 % critical damping per mode
            let omega_n = (modal_stiffness[i] / modal_mass[i]).sqrt();
            modal_damping[i] = 2.0 * 0.05 * modal_mass[i] * omega_n;
        }

        // Diagonal compliance: 1/k for diagonal, off-diagonal zero
        let mut compliance = [0.0_f64; 36];
        for i in 0..6 {
            let k = if i < NUM_MODES {
                modal_stiffness[i]
            } else {
                1.0e6
            };
            compliance[i * 6 + i] = 1.0 / k;
        }

        Self {
            core_mass,
            shell_thickness,
            youngs_modulus,
            poissons_ratio,
            shell_area,
            modal_stiffness,
            modal_damping,
            modal_mass,
            q: [0.0; NUM_MODES],
            q_dot: [0.0; NUM_MODES],
            compliance,
        }
    }

    /// Integrate modal equations of motion for one time step `dt` under
    /// generalised modal forces `f_modal`.
    ///
    /// Uses semi-implicit Euler.
    pub fn integrate(&mut self, f_modal: [f64; NUM_MODES], dt: f64) {
        for i in 0..NUM_MODES {
            let f_damp = -self.modal_damping[i] * self.q_dot[i];
            let f_spring = -self.modal_stiffness[i] * self.q[i];
            let q_ddot = (f_modal[i] + f_spring + f_damp) / self.modal_mass[i];
            self.q_dot[i] += q_ddot * dt;
            self.q[i] += self.q_dot[i] * dt;
        }
    }

    /// Compute bending stiffness D = E h³ / (12 (1 - ν²)).
    pub fn bending_stiffness(&self) -> f64 {
        self.youngs_modulus * self.shell_thickness.powi(3)
            / (12.0 * (1.0 - self.poissons_ratio * self.poissons_ratio))
    }

    /// Compute membrane stiffness C = E h / (1 - ν²).
    pub fn membrane_stiffness(&self) -> f64 {
        self.youngs_modulus * self.shell_thickness
            / (1.0 - self.poissons_ratio * self.poissons_ratio)
    }

    /// Return the current elastic strain energy (J).
    pub fn strain_energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..NUM_MODES {
            e += 0.5 * self.modal_stiffness[i] * self.q[i] * self.q[i];
        }
        e
    }

    /// Apply the compliance matrix to a 6-D wrench, returning the 6-D
    /// deformation vector.
    pub fn apply_compliance(&self, wrench: [f64; 6]) -> [f64; 6] {
        let mut deform = [0.0_f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                deform[i] += self.compliance[i * 6 + j] * wrench[j];
            }
        }
        deform
    }

    /// Reset all modal amplitudes and velocities to zero.
    pub fn reset_deformation(&mut self) {
        self.q = [0.0; NUM_MODES];
        self.q_dot = [0.0; NUM_MODES];
    }
}

// ---------------------------------------------------------------------------
// ContactCompliance
// ---------------------------------------------------------------------------

/// Contact compliance model combining Hertz elastic contact and a linear
/// compliant contact regularisation.
#[derive(Debug, Clone)]
pub struct ContactCompliance {
    /// Combined effective elastic modulus E* (Pa).
    pub effective_modulus: f64,
    /// Effective contact radius R* (m).
    pub effective_radius: f64,
    /// Linear stiffness for the compliant contact regularisation (N/m).
    pub stiffness: f64,
    /// Linear damping coefficient (N·s/m).
    pub damping: f64,
    /// Maximum allowed penetration before clamping (m).
    pub max_penetration: f64,
}

impl ContactCompliance {
    /// Construct a contact compliance from material and geometry data.
    ///
    /// * `e1`, `nu1` — Young's modulus and Poisson's ratio of body 1.
    /// * `e2`, `nu2` — Young's modulus and Poisson's ratio of body 2.
    /// * `r1`, `r2` — principal radii of curvature (m).
    pub fn new(e1: f64, nu1: f64, e2: f64, nu2: f64, r1: f64, r2: f64) -> Self {
        let inv_e_star = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
        let effective_modulus = 1.0 / inv_e_star;
        let effective_radius = if r1 > 0.0 && r2 > 0.0 {
            1.0 / (1.0 / r1 + 1.0 / r2)
        } else if r1 > 0.0 {
            r1
        } else {
            r2
        };
        let stiffness = (4.0 / 3.0) * effective_modulus * effective_radius.sqrt();
        let damping = 2.0 * 0.05 * (stiffness * 1.0).sqrt(); // approximate
        Self {
            effective_modulus,
            effective_radius,
            stiffness,
            damping,
            max_penetration: 1e-2,
        }
    }

    /// Compute Hertz contact force for penetration depth `delta` (m).
    ///
    /// F = (4/3) E* √R* δ^(3/2)
    pub fn hertz_force(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let delta = delta.min(self.max_penetration);
        (4.0 / 3.0) * self.effective_modulus * self.effective_radius.sqrt() * delta.powf(1.5)
    }

    /// Compute Hertz contact stiffness ∂F/∂δ at current penetration `delta`.
    pub fn hertz_stiffness(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let delta = delta.min(self.max_penetration);
        2.0 * self.effective_modulus * self.effective_radius.sqrt() * delta.sqrt()
    }

    /// Compliant contact force: linear spring-damper for small δ, Hertz for large δ.
    pub fn compliant_force(&self, delta: f64, delta_dot: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let elastic = self.hertz_force(delta);
        let viscous = (self.damping * delta_dot).max(0.0);
        elastic + viscous
    }

    /// Compute contact radius from Hertz theory: a = (3 F R* / (4 E*))^(1/3).
    pub fn hertz_contact_radius(&self, normal_force: f64) -> f64 {
        if normal_force <= 0.0 {
            return 0.0;
        }
        (3.0 * normal_force * self.effective_radius / (4.0 * self.effective_modulus))
            .powf(1.0 / 3.0)
    }

    /// Penetration depth from Hertz contact force (inverse of hertz_force).
    pub fn penetration_from_force(&self, force: f64) -> f64 {
        if force <= 0.0 {
            return 0.0;
        }
        let c = (4.0 / 3.0) * self.effective_modulus * self.effective_radius.sqrt();
        (force / c).powf(2.0 / 3.0)
    }
}

// ---------------------------------------------------------------------------
// RigidFlexiCoupling
// ---------------------------------------------------------------------------

/// Maps modal deformation amplitudes to resultant force and torque on the
/// rigid body via modal superposition.
#[derive(Debug, Clone)]
pub struct RigidFlexiCoupling {
    /// Number of deformation modes retained.
    pub n_modes: usize,
    /// Modal force participation vectors (n_modes × 3) stored in row-major order.
    /// `force_participation[i*3 .. i*3+3]` is the force direction for mode i.
    pub force_participation: Vec<f64>,
    /// Modal torque participation vectors (n_modes × 3), row-major.
    pub torque_participation: Vec<f64>,
    /// Natural frequencies ω_i (rad/s) for each mode.
    pub natural_frequencies: Vec<f64>,
    /// Damping ratios ζ_i for each mode.
    pub damping_ratios: Vec<f64>,
}

impl RigidFlexiCoupling {
    /// Construct a coupling with `n_modes` retained modes.
    ///
    /// All participation vectors default to unit vectors along X; callers
    /// should overwrite `force_participation` and `torque_participation`.
    pub fn new(n_modes: usize) -> Self {
        let mut force_participation = vec![0.0_f64; n_modes * 3];
        let torque_participation = vec![0.0_f64; n_modes * 3];
        let mut natural_frequencies = vec![0.0_f64; n_modes];
        let damping_ratios = vec![0.05_f64; n_modes];
        for i in 0..n_modes {
            force_participation[i * 3] = 1.0;
            natural_frequencies[i] = 100.0 * (i + 1) as f64; // placeholder Hz
        }
        Self {
            n_modes,
            force_participation,
            torque_participation,
            natural_frequencies,
            damping_ratios,
        }
    }

    /// Compute the resultant force on the rigid body from modal amplitudes `q`.
    ///
    /// F = Σ_i φ_i^F · k_i · q_i
    pub fn resultant_force(&self, q: &[f64], modal_stiffness: &[f64]) -> [f64; 3] {
        let mut force = [0.0_f64; 3];
        let n = q.len().min(self.n_modes);
        for i in 0..n {
            let k_q = modal_stiffness.get(i).copied().unwrap_or(1.0) * q[i];
            let phi = &self.force_participation[i * 3..i * 3 + 3];
            force[0] += phi[0] * k_q;
            force[1] += phi[1] * k_q;
            force[2] += phi[2] * k_q;
        }
        force
    }

    /// Compute the resultant torque on the rigid body from modal amplitudes `q`.
    pub fn resultant_torque(&self, q: &[f64], modal_stiffness: &[f64]) -> [f64; 3] {
        let mut torque = [0.0_f64; 3];
        let n = q.len().min(self.n_modes);
        for i in 0..n {
            let k_q = modal_stiffness.get(i).copied().unwrap_or(1.0) * q[i];
            let phi = &self.torque_participation[i * 3..i * 3 + 3];
            torque[0] += phi[0] * k_q;
            torque[1] += phi[1] * k_q;
            torque[2] += phi[2] * k_q;
        }
        torque
    }

    /// Compute modal generalised forces from an external force `f` applied
    /// at the surface, projected through participation vectors.
    pub fn project_force_to_modes(&self, f: [f64; 3]) -> Vec<f64> {
        let mut q_force = vec![0.0_f64; self.n_modes];
        for i in 0..self.n_modes {
            let phi = &self.force_participation[i * 3..i * 3 + 3];
            q_force[i] = phi[0] * f[0] + phi[1] * f[1] + phi[2] * f[2];
        }
        q_force
    }

    /// Reconstruct the physical deformation vector at a point from modal amplitudes.
    ///
    /// The displacement at the point is approximated as u = Σ_i φ_i^F q_i.
    pub fn reconstruct_displacement(&self, q: &[f64]) -> [f64; 3] {
        let mut disp = [0.0_f64; 3];
        let n = q.len().min(self.n_modes);
        for i in 0..n {
            let phi = &self.force_participation[i * 3..i * 3 + 3];
            disp[0] += phi[0] * q[i];
            disp[1] += phi[1] * q[i];
            disp[2] += phi[2] * q[i];
        }
        disp
    }
}

// ---------------------------------------------------------------------------
// ImpactDeformation
// ---------------------------------------------------------------------------

/// Models energy absorption and permanent plastic deformation during impact.
#[derive(Debug, Clone)]
pub struct ImpactDeformation {
    /// Yield force (N): below this force the response is elastic.
    pub yield_force: f64,
    /// Plastic stiffness after yield (N/m).
    pub plastic_stiffness: f64,
    /// Elastic stiffness (N/m).
    pub elastic_stiffness: f64,
    /// Accumulated permanent deformation (m).
    pub permanent_deformation: f64,
    /// Maximum deformation reached (m).
    pub max_deformation: f64,
    /// Energy absorbed so far (J).
    pub energy_absorbed: f64,
}

impl ImpactDeformation {
    /// Construct an impact deformation model.
    ///
    /// # Arguments
    /// * `elastic_stiffness` — linear stiffness in the elastic regime (N/m).
    /// * `yield_force` — force at which plastic yielding begins (N).
    /// * `plastic_stiffness` — post-yield hardening stiffness (N/m).
    pub fn new(elastic_stiffness: f64, yield_force: f64, plastic_stiffness: f64) -> Self {
        Self {
            yield_force,
            plastic_stiffness,
            elastic_stiffness,
            permanent_deformation: 0.0,
            max_deformation: 0.0,
            energy_absorbed: 0.0,
        }
    }

    /// Compute contact force and update state for current deformation `delta`.
    ///
    /// Uses an elastic-plastic model: elastic below yield, linear hardening above.
    /// Returns the contact force (N).
    pub fn update(&mut self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let yield_disp = self.yield_force / self.elastic_stiffness;
        let force = if delta <= yield_disp {
            self.elastic_stiffness * delta
        } else {
            self.yield_force + self.plastic_stiffness * (delta - yield_disp)
        };
        // Track energy absorbed (trapezoidal integration approximation)
        if delta > self.max_deformation {
            let d_delta = delta - self.max_deformation;
            self.energy_absorbed += force * d_delta * 0.5;
            self.max_deformation = delta;
        }
        // Update permanent deformation
        if delta > yield_disp {
            self.permanent_deformation = (delta - yield_disp).max(self.permanent_deformation);
        }
        force
    }

    /// Compute coefficient of restitution from absorbed vs. input kinetic energy.
    ///
    /// e = √(1 - E_absorbed / E_kinetic)
    pub fn restitution_coefficient(&self, kinetic_energy_in: f64) -> f64 {
        if kinetic_energy_in <= 0.0 {
            return 0.0;
        }
        let ratio = (self.energy_absorbed / kinetic_energy_in).min(1.0);
        (1.0 - ratio).sqrt()
    }

    /// Check whether the body has undergone any permanent plastic deformation.
    pub fn has_permanent_deformation(&self) -> bool {
        self.permanent_deformation > 1e-12
    }

    /// Reset accumulated state (useful for unit tests or re-use).
    pub fn reset(&mut self) {
        self.permanent_deformation = 0.0;
        self.max_deformation = 0.0;
        self.energy_absorbed = 0.0;
    }
}

// ---------------------------------------------------------------------------
// SoftContactPatch
// ---------------------------------------------------------------------------

/// Elliptical soft contact patch.
///
/// Models the contact between two deformable surfaces as an ellipse with
/// semi-axes `a` (major) and `b` (minor). Pressure distribution follows
/// the Hertzian parabolic profile p(r) = p0 √(1 - (r/a)²).
#[derive(Debug, Clone)]
pub struct SoftContactPatch {
    /// Semi-major axis of the contact ellipse (m).
    pub semi_major: f64,
    /// Semi-minor axis of the contact ellipse (m).
    pub semi_minor: f64,
    /// Peak contact pressure at the ellipse centre (Pa).
    pub peak_pressure: f64,
    /// Contact centre in body-local coordinates (m).
    pub centre: [f64; 3],
    /// Contact normal (unit vector pointing from surface into body).
    pub normal: [f64; 3],
    /// Friction coefficient for traction computation.
    pub friction: f64,
    /// Tangential slip velocity at the patch (m/s).
    pub slip_velocity: [f64; 3],
}

impl SoftContactPatch {
    /// Construct from Hertz theory given normal force and material properties.
    ///
    /// * `normal_force` — applied normal force (N).
    /// * `e_star` — combined effective modulus E* (Pa).
    /// * `r_star` — effective contact radius R* (m).
    /// * `friction` — Coulomb friction coefficient.
    pub fn from_hertz(normal_force: f64, e_star: f64, r_star: f64, friction: f64) -> Self {
        if normal_force <= 0.0 || e_star <= 0.0 || r_star <= 0.0 {
            return Self {
                semi_major: 0.0,
                semi_minor: 0.0,
                peak_pressure: 0.0,
                centre: [0.0; 3],
                normal: [0.0, 1.0, 0.0],
                friction,
                slip_velocity: [0.0; 3],
            };
        }
        let a = (3.0 * normal_force * r_star / (4.0 * e_star)).powf(1.0 / 3.0);
        let area = PI * a * a;
        let p0 = (3.0 * normal_force) / (2.0 * area);
        Self {
            semi_major: a,
            semi_minor: a,
            peak_pressure: p0,
            centre: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            friction,
            slip_velocity: [0.0; 3],
        }
    }

    /// Contact area (m²).
    pub fn area(&self) -> f64 {
        PI * self.semi_major * self.semi_minor
    }

    /// Total normal force integrated over the patch (N).
    ///
    /// For a Hertzian parabolic profile: F = (2/3) p0 A.
    pub fn total_normal_force(&self) -> f64 {
        (2.0 / 3.0) * self.peak_pressure * self.area()
    }

    /// Integrate traction vector over the contact patch.
    ///
    /// Returns the resultant traction force vector (N) in the plane
    /// perpendicular to the contact normal.
    pub fn integrate_traction(&self) -> [f64; 3] {
        let slip_mag = v3_norm(self.slip_velocity);
        if slip_mag < 1e-9 {
            return [0.0; 3];
        }
        let slip_dir = v3_normalise(self.slip_velocity);
        let total_normal = self.total_normal_force();
        let traction_mag = self.friction * total_normal;
        v3_scale(slip_dir, traction_mag)
    }

    /// Pressure at radial distance `r` from the patch centre (Pa).
    ///
    /// Uses the Hertzian parabolic profile.
    pub fn pressure_at(&self, r: f64) -> f64 {
        let a = self.semi_major;
        if a < 1e-15 || r >= a {
            return 0.0;
        }
        self.peak_pressure * (1.0 - (r / a).powi(2)).sqrt()
    }

    /// Average pressure over the contact patch (Pa).
    pub fn average_pressure(&self) -> f64 {
        let area = self.area();
        if area < 1e-20 {
            return 0.0;
        }
        self.total_normal_force() / area
    }
}

// ---------------------------------------------------------------------------
// WrenchSpace
// ---------------------------------------------------------------------------

/// A 6-D wrench: force (3 components) and moment (3 components) about a point.
#[derive(Debug, Clone, Copy)]
pub struct Wrench {
    /// Force component (N).
    pub force: [f64; 3],
    /// Moment component (N·m).
    pub moment: [f64; 3],
}

impl Wrench {
    /// Construct a zero wrench.
    pub fn zero() -> Self {
        Self {
            force: [0.0; 3],
            moment: [0.0; 3],
        }
    }

    /// Construct from separate force and moment vectors.
    pub fn new(force: [f64; 3], moment: [f64; 3]) -> Self {
        Self { force, moment }
    }

    /// Add two wrenches component-wise.
    pub fn add(&self, other: &Wrench) -> Wrench {
        Wrench {
            force: v3_add(self.force, other.force),
            moment: v3_add(self.moment, other.moment),
        }
    }

    /// Scale all components by scalar `s`.
    pub fn scale(&self, s: f64) -> Wrench {
        Wrench {
            force: v3_scale(self.force, s),
            moment: v3_scale(self.moment, s),
        }
    }

    /// Return the wrench as a flat 6-element array \[fx,fy,fz,mx,my,mz\].
    pub fn as_array(&self) -> [f64; 6] {
        [
            self.force[0],
            self.force[1],
            self.force[2],
            self.moment[0],
            self.moment[1],
            self.moment[2],
        ]
    }

    /// Magnitude of the force component.
    pub fn force_magnitude(&self) -> f64 {
        v3_norm(self.force)
    }

    /// Magnitude of the moment component.
    pub fn moment_magnitude(&self) -> f64 {
        v3_norm(self.moment)
    }

    /// Transport the wrench from its current reference point to a new point `p`.
    ///
    /// New moment = old moment + r × force where r = contact_point - p.
    pub fn transport_to(&self, contact_point: [f64; 3], new_point: [f64; 3]) -> Wrench {
        let r = v3_sub(contact_point, new_point);
        let cross = v3_cross(r, self.force);
        Wrench {
            force: self.force,
            moment: v3_add(self.moment, cross),
        }
    }
}

/// Grasp map: maps contact wrenches to a net object wrench.
///
/// Stores `n_contacts` contact frames (normal + two tangent directions)
/// and assembles the 6 × (n_contacts · 3) grasp matrix G.
#[derive(Debug, Clone)]
pub struct GraspMap {
    /// Number of contact points.
    pub n_contacts: usize,
    /// Contact points in object frame (m).
    pub contact_points: Vec<[f64; 3]>,
    /// Contact normals (unit vectors, into object).
    pub contact_normals: Vec<[f64; 3]>,
    /// First tangent directions at each contact.
    pub tangent1: Vec<[f64; 3]>,
    /// Second tangent directions at each contact.
    pub tangent2: Vec<[f64; 3]>,
}

impl GraspMap {
    /// Construct from contact data.
    ///
    /// Tangent vectors are computed automatically from each normal via the
    /// Gram-Schmidt process.
    pub fn new(contact_points: Vec<[f64; 3]>, contact_normals: Vec<[f64; 3]>) -> Self {
        let n = contact_points.len();
        let mut tangent1 = Vec::with_capacity(n);
        let mut tangent2 = Vec::with_capacity(n);
        for i in 0..n {
            let n_hat = v3_normalise(contact_normals[i]);
            // Choose an auxiliary vector not parallel to n_hat
            let aux = if n_hat[0].abs() < 0.9 {
                [1.0_f64, 0.0, 0.0]
            } else {
                [0.0_f64, 1.0, 0.0]
            };
            let t1 = v3_normalise(v3_sub(aux, v3_scale(n_hat, v3_dot(aux, n_hat))));
            let t2 = v3_cross(n_hat, t1);
            tangent1.push(t1);
            tangent2.push(t2);
        }
        Self {
            n_contacts: n,
            contact_points,
            contact_normals,
            tangent1,
            tangent2,
        }
    }

    /// Map a set of contact force vectors (one per contact) to a net wrench
    /// about the origin.
    ///
    /// Each `contact_forces[i]` is expressed in world coordinates.
    pub fn map_to_wrench(&self, contact_forces: &[[f64; 3]]) -> Wrench {
        let mut net = Wrench::zero();
        let n = contact_forces.len().min(self.n_contacts);
        for i in 0..n {
            let f = contact_forces[i];
            let r = self.contact_points[i];
            let torque = v3_cross(r, f);
            net.force = v3_add(net.force, f);
            net.moment = v3_add(net.moment, torque);
        }
        net
    }

    /// Check whether a given wrench is inside the contact wrench cone (CWC).
    ///
    /// A wrench is feasible if it can be generated by non-negative normal
    /// forces at all contacts with friction cone constraints.  This is a
    /// simplified linearised check: verifies that the normal component at
    /// each contact is non-negative and friction limits are satisfied.
    pub fn in_contact_wrench_cone(&self, wrench: &Wrench, friction_coeff: f64) -> bool {
        // For each contact, compute the component of the net force along the normal
        // This is a simplified check — a full CWC requires convex optimisation
        let n_hat = &self.contact_normals;
        if n_hat.is_empty() {
            return false;
        }
        for i in 0..self.n_contacts {
            let fn_i = v3_dot(wrench.force, n_hat[i]);
            if fn_i < 0.0 {
                return false;
            }
            // Friction cone: |tangential force| ≤ μ * fn
            let ft1 = v3_dot(wrench.force, self.tangent1[i]).abs();
            let ft2 = v3_dot(wrench.force, self.tangent2[i]).abs();
            if ft1 + ft2 > friction_coeff * fn_i + 1e-9 {
                return false;
            }
        }
        true
    }
}

/// 6-D wrench space utilities.
#[derive(Debug, Clone)]
pub struct WrenchSpace {
    /// Grasp map associated with this wrench space.
    pub grasp_map: GraspMap,
}

impl WrenchSpace {
    /// Construct from a set of contact points and normals.
    pub fn new(contact_points: Vec<[f64; 3]>, contact_normals: Vec<[f64; 3]>) -> Self {
        Self {
            grasp_map: GraspMap::new(contact_points, contact_normals),
        }
    }

    /// Return the number of contact points.
    pub fn n_contacts(&self) -> usize {
        self.grasp_map.n_contacts
    }

    /// Compute the wrench from contact forces.
    pub fn compute_wrench(&self, contact_forces: &[[f64; 3]]) -> Wrench {
        self.grasp_map.map_to_wrench(contact_forces)
    }

    /// Check if a wrench is feasible given the friction coefficient.
    pub fn is_feasible(&self, wrench: &Wrench, friction_coeff: f64) -> bool {
        self.grasp_map
            .in_contact_wrench_cone(wrench, friction_coeff)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── DeformableBody ───────────────────────────────────────────────────

    #[test]
    fn deformable_body_bending_stiffness_positive() {
        let body = DeformableBody::new(10.0, 0.005, 70e9, 0.33, 0.5);
        assert!(body.bending_stiffness() > 0.0);
    }

    #[test]
    fn deformable_body_membrane_stiffness_positive() {
        let body = DeformableBody::new(10.0, 0.005, 70e9, 0.33, 0.5);
        assert!(body.membrane_stiffness() > 0.0);
    }

    #[test]
    fn deformable_body_initial_strain_energy_zero() {
        let body = DeformableBody::new(5.0, 0.003, 200e9, 0.3, 1.0);
        assert!(body.strain_energy().abs() < 1e-10);
    }

    #[test]
    fn deformable_body_integrate_moves_modes() {
        let mut body = DeformableBody::new(1.0, 0.01, 10e9, 0.3, 0.2);
        let mut f = [0.0_f64; NUM_MODES];
        f[0] = 100.0;
        body.integrate(f, 0.001);
        assert!(
            body.q[0].abs() > 0.0 || body.q_dot[0].abs() > 0.0,
            "mode should have moved"
        );
    }

    #[test]
    fn deformable_body_reset_zeros_modes() {
        let mut body = DeformableBody::new(1.0, 0.01, 10e9, 0.3, 0.2);
        body.q[0] = 1.0;
        body.q_dot[0] = 2.0;
        body.reset_deformation();
        assert_eq!(body.q, [0.0; NUM_MODES]);
        assert_eq!(body.q_dot, [0.0; NUM_MODES]);
    }

    #[test]
    fn deformable_body_compliance_matrix_diagonal() {
        let body = DeformableBody::new(5.0, 0.005, 70e9, 0.3, 0.5);
        // diagonal entries should be positive
        for i in 0..6 {
            assert!(
                body.compliance[i * 6 + i] > 0.0,
                "compliance[{i},{i}] should be positive"
            );
        }
    }

    #[test]
    fn deformable_body_apply_compliance_scales_wrench() {
        let body = DeformableBody::new(5.0, 0.005, 70e9, 0.3, 0.5);
        let wrench = [1.0_f64; 6];
        let deform = body.apply_compliance(wrench);
        // All deformations should be non-zero
        for i in 0..6 {
            assert!(deform[i].abs() > 0.0, "deform[{i}] should be non-zero");
        }
    }

    #[test]
    fn deformable_body_strain_energy_increases_with_amplitude() {
        let mut body = DeformableBody::new(1.0, 0.01, 10e9, 0.3, 0.2);
        body.q[0] = 1e-4;
        let e1 = body.strain_energy();
        body.q[0] = 2e-4;
        let e2 = body.strain_energy();
        assert!(e2 > e1, "strain energy should increase with amplitude");
    }

    #[test]
    fn deformable_body_num_modes_constant() {
        assert_eq!(NUM_MODES, 6);
    }

    // ── ContactCompliance ────────────────────────────────────────────────

    #[test]
    fn contact_compliance_hertz_force_zero_for_negative_delta() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        assert_eq!(cc.hertz_force(-0.001), 0.0);
    }

    #[test]
    fn contact_compliance_hertz_force_positive_for_positive_delta() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        let f = cc.hertz_force(0.001);
        assert!(f > 0.0, "Hertz force should be positive: {f}");
    }

    #[test]
    fn contact_compliance_hertz_stiffness_zero_for_negative_delta() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        assert_eq!(cc.hertz_stiffness(-0.001), 0.0);
    }

    #[test]
    fn contact_compliance_hertz_stiffness_increases_with_delta() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        let k1 = cc.hertz_stiffness(0.0005);
        let k2 = cc.hertz_stiffness(0.001);
        assert!(k2 > k1, "stiffness should increase with penetration");
    }

    #[test]
    fn contact_compliance_compliant_force_ge_hertz() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        let f_h = cc.hertz_force(0.001);
        let f_c = cc.compliant_force(0.001, 0.01);
        assert!(
            f_c >= f_h,
            "compliant force >= hertz force when delta_dot>0"
        );
    }

    #[test]
    fn contact_compliance_contact_radius_positive() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        let a = cc.hertz_contact_radius(1000.0);
        assert!(a > 0.0, "contact radius should be positive: {a}");
    }

    #[test]
    fn contact_compliance_penetration_from_force_roundtrip() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        let delta_in = 0.0005;
        let force = cc.hertz_force(delta_in);
        let delta_out = cc.penetration_from_force(force);
        assert!(
            (delta_out - delta_in).abs() < 1e-10,
            "roundtrip: delta_in={delta_in}, delta_out={delta_out}"
        );
    }

    #[test]
    fn contact_compliance_effective_modulus_finite() {
        let cc = ContactCompliance::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01);
        assert!(cc.effective_modulus.is_finite() && cc.effective_modulus > 0.0);
    }

    // ── RigidFlexiCoupling ───────────────────────────────────────────────

    #[test]
    fn rigid_flexi_coupling_resultant_force_zero_at_rest() {
        let coupling = RigidFlexiCoupling::new(NUM_MODES);
        let q = [0.0_f64; NUM_MODES];
        let k = [1.0_f64; NUM_MODES];
        let f = coupling.resultant_force(&q, &k);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn rigid_flexi_coupling_resultant_force_nonzero_when_displaced() {
        let coupling = RigidFlexiCoupling::new(NUM_MODES);
        let mut q = [0.0_f64; NUM_MODES];
        q[0] = 0.001;
        let k = [1e6_f64; NUM_MODES];
        let f = coupling.resultant_force(&q, &k);
        assert!(v3_norm(f) > 0.0, "force should be nonzero: {f:?}");
    }

    #[test]
    fn rigid_flexi_coupling_project_force_to_modes_length() {
        let coupling = RigidFlexiCoupling::new(NUM_MODES);
        let projected = coupling.project_force_to_modes([10.0, 0.0, 0.0]);
        assert_eq!(projected.len(), NUM_MODES);
    }

    #[test]
    fn rigid_flexi_coupling_reconstruct_displacement_zero_at_rest() {
        let coupling = RigidFlexiCoupling::new(NUM_MODES);
        let q = vec![0.0_f64; NUM_MODES];
        let d = coupling.reconstruct_displacement(&q);
        assert_eq!(d, [0.0; 3]);
    }

    #[test]
    fn rigid_flexi_coupling_resultant_torque_zero_at_rest() {
        let coupling = RigidFlexiCoupling::new(NUM_MODES);
        let q = [0.0_f64; NUM_MODES];
        let k = [1.0_f64; NUM_MODES];
        let tau = coupling.resultant_torque(&q, &k);
        assert_eq!(tau, [0.0; 3]);
    }

    // ── ImpactDeformation ────────────────────────────────────────────────

    #[test]
    fn impact_deformation_zero_force_at_zero_delta() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        assert_eq!(imp.update(0.0), 0.0);
    }

    #[test]
    fn impact_deformation_elastic_regime() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        // Within elastic regime: δ < F_yield / k_elastic = 5000/1e6 = 0.005 m
        let delta = 0.003;
        let force = imp.update(delta);
        let expected = 1e6 * delta;
        assert!(
            (force - expected).abs() < 1e-3,
            "force={force}, expected={expected}"
        );
    }

    #[test]
    fn impact_deformation_plastic_regime() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        // Beyond yield: δ > 0.005
        let force = imp.update(0.01);
        assert!(
            force > 5000.0,
            "force in plastic regime should exceed yield force"
        );
    }

    #[test]
    fn impact_deformation_energy_absorbed_increases() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        imp.update(0.005);
        let e1 = imp.energy_absorbed;
        imp.update(0.01);
        let e2 = imp.energy_absorbed;
        assert!(e2 >= e1, "energy absorbed should not decrease");
    }

    #[test]
    fn impact_deformation_restitution_zero_when_all_absorbed() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        imp.energy_absorbed = 100.0;
        let e = imp.restitution_coefficient(100.0);
        assert!(e.abs() < 1e-10, "e={e}");
    }

    #[test]
    fn impact_deformation_restitution_one_when_nothing_absorbed() {
        let imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        let e = imp.restitution_coefficient(100.0);
        assert!((e - 1.0).abs() < 1e-10, "e={e}");
    }

    #[test]
    fn impact_deformation_has_permanent_deformation_after_plastic() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        imp.update(0.02);
        assert!(imp.has_permanent_deformation());
    }

    #[test]
    fn impact_deformation_no_permanent_deformation_elastic_only() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        imp.update(0.001);
        assert!(!imp.has_permanent_deformation());
    }

    #[test]
    fn impact_deformation_reset_clears_state() {
        let mut imp = ImpactDeformation::new(1e6, 5000.0, 2e5);
        imp.update(0.02);
        imp.reset();
        assert_eq!(imp.permanent_deformation, 0.0);
        assert_eq!(imp.energy_absorbed, 0.0);
        assert_eq!(imp.max_deformation, 0.0);
    }

    // ── SoftContactPatch ─────────────────────────────────────────────────

    #[test]
    fn soft_contact_patch_area_positive_for_valid_contact() {
        let patch = SoftContactPatch::from_hertz(1000.0, 200e9, 0.01, 0.3);
        assert!(patch.area() > 0.0, "area={}", patch.area());
    }

    #[test]
    fn soft_contact_patch_area_zero_for_zero_force() {
        let patch = SoftContactPatch::from_hertz(0.0, 200e9, 0.01, 0.3);
        assert_eq!(patch.area(), 0.0);
    }

    #[test]
    fn soft_contact_patch_total_normal_force_matches_input() {
        // For a circular Hertz patch, total_normal_force ≈ normal_force
        let normal_force = 1000.0;
        let patch = SoftContactPatch::from_hertz(normal_force, 200e9, 0.01, 0.3);
        let computed = patch.total_normal_force();
        // Allow 1% relative error due to simplified formula
        assert!(
            (computed - normal_force).abs() / normal_force < 0.02,
            "computed={computed}, expected≈{normal_force}"
        );
    }

    #[test]
    fn soft_contact_patch_pressure_at_centre_is_peak() {
        let patch = SoftContactPatch::from_hertz(1000.0, 200e9, 0.01, 0.3);
        let p0 = patch.pressure_at(0.0);
        assert!(
            (p0 - patch.peak_pressure).abs() < 1e-6,
            "p0={p0}, peak={}",
            patch.peak_pressure
        );
    }

    #[test]
    fn soft_contact_patch_pressure_zero_outside_contact() {
        let patch = SoftContactPatch::from_hertz(1000.0, 200e9, 0.01, 0.3);
        let p = patch.pressure_at(patch.semi_major + 1e-5);
        assert_eq!(p, 0.0, "pressure outside contact should be zero");
    }

    #[test]
    fn soft_contact_patch_integrate_traction_zero_no_slip() {
        let patch = SoftContactPatch {
            semi_major: 0.001,
            semi_minor: 0.001,
            peak_pressure: 1e6,
            centre: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            friction: 0.5,
            slip_velocity: [0.0; 3],
        };
        let traction = patch.integrate_traction();
        assert_eq!(traction, [0.0; 3]);
    }

    #[test]
    fn soft_contact_patch_average_pressure_less_than_peak() {
        let patch = SoftContactPatch::from_hertz(1000.0, 200e9, 0.01, 0.3);
        assert!(patch.average_pressure() < patch.peak_pressure);
    }

    // ── Wrench ───────────────────────────────────────────────────────────

    #[test]
    fn wrench_zero_is_zero() {
        let w = Wrench::zero();
        assert_eq!(w.force, [0.0; 3]);
        assert_eq!(w.moment, [0.0; 3]);
    }

    #[test]
    fn wrench_add() {
        let a = Wrench::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let b = Wrench::new([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]);
        let c = a.add(&b);
        assert_eq!(c.force, [2.0, 3.0, 4.0]);
        assert_eq!(c.moment, [5.0, 6.0, 7.0]);
    }

    #[test]
    fn wrench_scale() {
        let w = Wrench::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let ws = w.scale(2.0);
        assert_eq!(ws.force, [2.0, 4.0, 6.0]);
        assert_eq!(ws.moment, [8.0, 10.0, 12.0]);
    }

    #[test]
    fn wrench_as_array() {
        let w = Wrench::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let arr = w.as_array();
        assert_eq!(arr, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn wrench_force_magnitude() {
        let w = Wrench::new([3.0, 4.0, 0.0], [0.0; 3]);
        assert!((w.force_magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn wrench_transport_pure_force_generates_moment() {
        // Force [0,0,F] at x=1,y=0,z=0. Transport to origin.
        // moment = r × F = [1,0,0] × [0,0,F] = [0*F - 0*0, 0*0 - 1*F, 1*0 - 0*0]
        //        = [0, -F, 0]
        let w = Wrench::new([0.0, 0.0, 10.0], [0.0; 3]);
        let tw = w.transport_to([1.0, 0.0, 0.0], [0.0; 3]);
        assert!(
            (tw.moment[1] - (-10.0)).abs() < 1e-10,
            "moment.y={}",
            tw.moment[1]
        );
    }

    // ── GraspMap & WrenchSpace ───────────────────────────────────────────

    #[test]
    fn grasp_map_two_contacts_map_to_wrench() {
        let pts = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let normals = vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let gm = GraspMap::new(pts, normals);
        let forces = vec![[0.0, 5.0, 0.0], [0.0, 5.0, 0.0]];
        let w = gm.map_to_wrench(&forces);
        assert!(
            (w.force[1] - 10.0).abs() < 1e-10,
            "net fy should be 10, got {}",
            w.force[1]
        );
    }

    #[test]
    fn wrench_space_n_contacts() {
        let ws = WrenchSpace::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        );
        assert_eq!(ws.n_contacts(), 2);
    }

    #[test]
    fn wrench_space_feasible_normal_force() {
        let ws = WrenchSpace::new(vec![[0.0; 3]], vec![[0.0, 1.0, 0.0]]);
        // A force purely along the normal should be feasible
        let w = Wrench::new([0.0, 10.0, 0.0], [0.0; 3]);
        assert!(ws.is_feasible(&w, 0.5));
    }

    #[test]
    fn wrench_space_infeasible_opposing_normal() {
        let ws = WrenchSpace::new(vec![[0.0; 3]], vec![[0.0, 1.0, 0.0]]);
        // Force opposite the contact normal: infeasible
        let w = Wrench::new([0.0, -10.0, 0.0], [0.0; 3]);
        assert!(!ws.is_feasible(&w, 0.5));
    }

    #[test]
    fn grasp_map_tangent_orthogonal_to_normal() {
        let pts = vec![[0.0; 3]];
        let normals = vec![[0.0, 0.0, 1.0]];
        let gm = GraspMap::new(pts, normals);
        let dot = v3_dot(gm.tangent1[0], gm.contact_normals[0]);
        assert!(
            dot.abs() < 1e-10,
            "tangent1 should be perpendicular to normal, dot={dot}"
        );
    }
}
