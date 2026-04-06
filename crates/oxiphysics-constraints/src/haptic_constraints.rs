// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Haptic feedback constraints module.
//!
//! Implements proxy-based haptic rendering (Zilles-Salisbury), friction models,
//! texture rendering, impedance/admittance control, passivity-based stability
//! analysis, and multi-point haptic interaction for grasping.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vec3H — lightweight 3-D vector (no nalgebra)
// ---------------------------------------------------------------------------

/// A simple 3-D vector for haptic computations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3H {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3H {
    /// Create a new vector.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Squared magnitude.
    pub fn len_sq(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Euclidean magnitude.
    pub fn len(self) -> f64 {
        self.len_sq().sqrt()
    }

    /// Normalise. Returns zero vector if magnitude is near-zero.
    pub fn normalised(self) -> Self {
        let l = self.len();
        if l < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / l, self.y / l, self.z / l)
        }
    }

    /// Dot product.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Clamp magnitude.
    pub fn clamped(self, max_len: f64) -> Self {
        let l = self.len();
        if l > max_len && l > 1e-15 {
            let s = max_len / l;
            Self::new(self.x * s, self.y * s, self.z * s)
        } else {
            self
        }
    }

    /// Component-wise absolute value.
    pub fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Distance to another point.
    pub fn dist(self, other: Self) -> f64 {
        (self - other).len()
    }
}

impl std::ops::Add for Vec3H {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl std::ops::AddAssign for Vec3H {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl std::ops::Sub for Vec3H {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl std::ops::Neg for Vec3H {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl std::ops::Mul<f64> for Vec3H {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl std::ops::Div<f64> for Vec3H {
    type Output = Self;
    fn div(self, s: f64) -> Self {
        Self::new(self.x / s, self.y / s, self.z / s)
    }
}

// ---------------------------------------------------------------------------
// Surface representation for proxy method
// ---------------------------------------------------------------------------

/// An implicit surface represented as a plane (point + normal).
#[derive(Debug, Clone, Copy)]
pub struct SurfacePlane {
    /// A point on the surface.
    pub point: Vec3H,
    /// Outward normal (unit vector).
    pub normal: Vec3H,
}

impl SurfacePlane {
    /// Create a new surface plane.
    pub fn new(point: Vec3H, normal: Vec3H) -> Self {
        Self {
            point,
            normal: normal.normalised(),
        }
    }

    /// Signed distance from a point to the plane (positive = outside).
    pub fn signed_distance(&self, p: Vec3H) -> f64 {
        (p - self.point).dot(self.normal)
    }

    /// Project a point onto the plane surface.
    pub fn project(&self, p: Vec3H) -> Vec3H {
        let d = self.signed_distance(p);
        p - self.normal * d
    }
}

/// A sphere surface.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceSphere {
    /// Centre of the sphere.
    pub centre: Vec3H,
    /// Radius.
    pub radius: f64,
}

impl SurfaceSphere {
    /// Create a new sphere surface.
    pub fn new(centre: Vec3H, radius: f64) -> Self {
        Self { centre, radius }
    }

    /// Signed distance from a point (positive = outside).
    pub fn signed_distance(&self, p: Vec3H) -> f64 {
        (p - self.centre).len() - self.radius
    }

    /// Project a point onto the sphere surface.
    pub fn project(&self, p: Vec3H) -> Vec3H {
        let dir = (p - self.centre).normalised();
        self.centre + dir * self.radius
    }

    /// Surface normal at a point (pointing outward).
    pub fn normal_at(&self, p: Vec3H) -> Vec3H {
        (p - self.centre).normalised()
    }
}

// ---------------------------------------------------------------------------
// HapticProxy
// ---------------------------------------------------------------------------

/// Haptic Interaction Point (HIP) / proxy that tracks the constrained position.
#[derive(Debug, Clone)]
pub struct HapticProxy {
    /// Current proxy (constrained) position on or outside the surface.
    pub proxy_pos: Vec3H,
    /// Current device (unconstrained) position reported by the haptic device.
    pub device_pos: Vec3H,
    /// Virtual coupling stiffness (N/m).
    pub stiffness: f64,
    /// Virtual coupling damping (Ns/m).
    pub damping: f64,
    /// Device velocity (estimated or measured).
    pub device_vel: Vec3H,
    /// Previous device position (for velocity estimation).
    prev_device_pos: Vec3H,
}

impl HapticProxy {
    /// Create a new haptic proxy at the given starting position.
    pub fn new(start_pos: Vec3H, stiffness: f64, damping: f64) -> Self {
        Self {
            proxy_pos: start_pos,
            device_pos: start_pos,
            stiffness,
            damping,
            device_vel: Vec3H::zero(),
            prev_device_pos: start_pos,
        }
    }

    /// Update the device position and estimate velocity.
    pub fn update_device(&mut self, new_device_pos: Vec3H, dt: f64) {
        self.prev_device_pos = self.device_pos;
        self.device_pos = new_device_pos;
        if dt > 1e-15 {
            self.device_vel = (self.device_pos - self.prev_device_pos) / dt;
        }
    }

    /// Compute the spring-damper force from proxy to device: F = K*(proxy - device) + B*v.
    pub fn compute_force(&self) -> Vec3H {
        let spring = (self.proxy_pos - self.device_pos) * self.stiffness;
        let damp = self.device_vel * (-self.damping);
        spring + damp
    }

    /// Penetration depth (negative means no contact).
    pub fn penetration_depth(&self) -> f64 {
        (self.device_pos - self.proxy_pos).len()
    }

    /// Whether the device has penetrated the surface (proxy != device).
    pub fn in_contact(&self) -> bool {
        self.penetration_depth() > 1e-10
    }
}

// ---------------------------------------------------------------------------
// ProxyMethod (Zilles-Salisbury)
// ---------------------------------------------------------------------------

/// Zilles-Salisbury proxy method for god-object haptic rendering.
///
/// The proxy is constrained to remain on the surface of virtual objects.
/// When the device penetrates, the proxy stays on the surface and a spring
/// force is computed between proxy and device.
#[derive(Debug, Clone)]
pub struct ProxyMethod {
    /// The proxy state.
    pub proxy: HapticProxy,
    /// Maximum force magnitude to output (safety clamp).
    pub max_force: f64,
}

impl ProxyMethod {
    /// Create a new proxy method renderer.
    pub fn new(start_pos: Vec3H, stiffness: f64, damping: f64, max_force: f64) -> Self {
        Self {
            proxy: HapticProxy::new(start_pos, stiffness, damping),
            max_force,
        }
    }

    /// Update against a planar surface: constrain proxy if device has penetrated.
    pub fn update_plane(&mut self, device_pos: Vec3H, dt: f64, surface: &SurfacePlane) -> Vec3H {
        self.proxy.update_device(device_pos, dt);
        let sd = surface.signed_distance(device_pos);
        if sd < 0.0 {
            // Device is inside → project proxy onto surface
            self.proxy.proxy_pos = surface.project(device_pos);
        } else {
            // Device is outside → proxy follows device (free space)
            self.proxy.proxy_pos = device_pos;
        }
        self.proxy.compute_force().clamped(self.max_force)
    }

    /// Update against a sphere surface.
    pub fn update_sphere(&mut self, device_pos: Vec3H, dt: f64, surface: &SurfaceSphere) -> Vec3H {
        self.proxy.update_device(device_pos, dt);
        let sd = surface.signed_distance(device_pos);
        if sd < 0.0 {
            self.proxy.proxy_pos = surface.project(device_pos);
        } else {
            self.proxy.proxy_pos = device_pos;
        }
        self.proxy.compute_force().clamped(self.max_force)
    }

    /// Raw force magnitude currently being rendered.
    pub fn force_magnitude(&self) -> f64 {
        self.proxy.compute_force().len()
    }
}

// ---------------------------------------------------------------------------
// FrictionHaptic
// ---------------------------------------------------------------------------

/// Coulomb friction model with stick-slip behaviour for haptic rendering.
#[derive(Debug, Clone)]
pub struct FrictionHaptic {
    /// Static friction coefficient.
    pub mu_s: f64,
    /// Dynamic/kinetic friction coefficient.
    pub mu_k: f64,
    /// Whether currently in stick state.
    pub is_stuck: bool,
    /// Anchor point (position where stick began).
    pub anchor: Vec3H,
    /// Lateral stiffness during stick phase (N/m).
    pub lateral_stiffness: f64,
}

impl FrictionHaptic {
    /// Create a new friction model.
    pub fn new(mu_s: f64, mu_k: f64, lateral_stiffness: f64) -> Self {
        Self {
            mu_s,
            mu_k,
            is_stuck: true,
            anchor: Vec3H::zero(),
            lateral_stiffness,
        }
    }

    /// Set the anchor point (called on first contact).
    pub fn set_anchor(&mut self, pos: Vec3H) {
        self.anchor = pos;
        self.is_stuck = true;
    }

    /// Compute friction force given the proxy position on the surface,
    /// the normal force magnitude, and the tangential displacement.
    ///
    /// Returns the friction force vector (tangential to surface).
    pub fn compute_friction(
        &mut self,
        proxy_pos: Vec3H,
        normal_force_mag: f64,
        surface_normal: Vec3H,
    ) -> Vec3H {
        // Tangential displacement from anchor
        let disp = proxy_pos - self.anchor;
        let normal_component = surface_normal * disp.dot(surface_normal);
        let tangential_disp = disp - normal_component;
        let tang_dist = tangential_disp.len();

        let max_static_force = self.mu_s * normal_force_mag;
        let spring_force = tang_dist * self.lateral_stiffness;

        if self.is_stuck {
            if spring_force > max_static_force && tang_dist > 1e-15 {
                // Transition to slip
                self.is_stuck = false;
                let kinetic_force = self.mu_k * normal_force_mag;
                -tangential_disp.normalised() * kinetic_force
            } else {
                // Stick: spring-like restoring force
                -tangential_disp * self.lateral_stiffness
            }
        } else {
            // Kinetic friction
            let kinetic_force = self.mu_k * normal_force_mag;
            if tang_dist < 1e-15 {
                // Re-enter stick when stationary
                self.is_stuck = true;
                self.anchor = proxy_pos;
                Vec3H::zero()
            } else {
                // Check if we should re-stick (velocity very low)
                if spring_force < max_static_force * 0.5 {
                    self.is_stuck = true;
                    self.anchor = proxy_pos;
                    Vec3H::zero()
                } else {
                    -tangential_disp.normalised() * kinetic_force
                }
            }
        }
    }

    /// Reset friction state.
    pub fn reset(&mut self) {
        self.is_stuck = true;
        self.anchor = Vec3H::zero();
    }
}

// ---------------------------------------------------------------------------
// HapticTexture
// ---------------------------------------------------------------------------

/// Haptic texture rendering using bump-map / sinusoidal perturbation.
#[derive(Debug, Clone)]
pub struct HapticTexture {
    /// Spatial frequency in X (cycles per metre).
    pub freq_x: f64,
    /// Spatial frequency in Y (cycles per metre).
    pub freq_y: f64,
    /// Bump amplitude (metres).
    pub amplitude: f64,
    /// Texture stiffness gain.
    pub stiffness: f64,
}

impl HapticTexture {
    /// Create a new texture.
    pub fn new(freq_x: f64, freq_y: f64, amplitude: f64, stiffness: f64) -> Self {
        Self {
            freq_x,
            freq_y,
            amplitude,
            stiffness,
        }
    }

    /// Evaluate the height field at a surface point (x, y).
    pub fn height(&self, x: f64, y: f64) -> f64 {
        self.amplitude
            * ((2.0 * PI * self.freq_x * x).sin() + (2.0 * PI * self.freq_y * y).sin())
            * 0.5
    }

    /// Compute the texture gradient at a surface point.
    pub fn gradient(&self, x: f64, y: f64) -> (f64, f64) {
        let dhdx =
            self.amplitude * 0.5 * 2.0 * PI * self.freq_x * (2.0 * PI * self.freq_x * x).cos();
        let dhdy =
            self.amplitude * 0.5 * 2.0 * PI * self.freq_y * (2.0 * PI * self.freq_y * y).cos();
        (dhdx, dhdy)
    }

    /// Compute the texture force given the proxy position on a surface
    /// with local tangent vectors `tx`, `ty` and normal `n`.
    ///
    /// The force perturbs the normal direction based on the bump gradient.
    pub fn compute_force(
        &self,
        proxy_x: f64,
        proxy_y: f64,
        surface_normal: Vec3H,
        tangent_x: Vec3H,
        tangent_y: Vec3H,
        normal_force_mag: f64,
    ) -> Vec3H {
        let (gx, gy) = self.gradient(proxy_x, proxy_y);
        // Perturbed normal
        let perturbed = surface_normal + tangent_x * (-gx) + tangent_y * (-gy);
        let perturbed_n = perturbed.normalised();
        // Force proportional to normal force magnitude in the perturbed direction
        let delta_n = perturbed_n - surface_normal;
        delta_n * normal_force_mag * self.stiffness
    }

    /// Wavelength in X direction.
    pub fn wavelength_x(&self) -> f64 {
        if self.freq_x.abs() < 1e-15 {
            f64::INFINITY
        } else {
            1.0 / self.freq_x
        }
    }

    /// Wavelength in Y direction.
    pub fn wavelength_y(&self) -> f64 {
        if self.freq_y.abs() < 1e-15 {
            f64::INFINITY
        } else {
            1.0 / self.freq_y
        }
    }
}

// ---------------------------------------------------------------------------
// ImpedanceControl
// ---------------------------------------------------------------------------

/// Impedance-type haptic rendering: the device controls position, the
/// controller outputs force.
///
/// Renders a virtual mass-spring-damper: `F = -K_d * x - B_d * v - M_d * a`.
#[derive(Debug, Clone)]
pub struct ImpedanceControl {
    /// Desired inertia (kg).
    pub m_d: f64,
    /// Desired damping (Ns/m).
    pub b_d: f64,
    /// Desired stiffness (N/m).
    pub k_d: f64,
    /// Equilibrium (target) position.
    pub target_pos: Vec3H,
    /// Previous velocity for acceleration estimation.
    prev_vel: Vec3H,
}

impl ImpedanceControl {
    /// Create a new impedance controller.
    pub fn new(m_d: f64, b_d: f64, k_d: f64, target_pos: Vec3H) -> Self {
        Self {
            m_d,
            b_d,
            k_d,
            target_pos,
            prev_vel: Vec3H::zero(),
        }
    }

    /// Compute the impedance force given device position and velocity.
    ///
    /// `F = -K_d * (pos - target) - B_d * vel - M_d * accel`
    pub fn compute_force(&mut self, pos: Vec3H, vel: Vec3H, dt: f64) -> Vec3H {
        let x_err = pos - self.target_pos;
        let accel = if dt > 1e-15 {
            (vel - self.prev_vel) / dt
        } else {
            Vec3H::zero()
        };
        self.prev_vel = vel;
        -(x_err * self.k_d) - (vel * self.b_d) - (accel * self.m_d)
    }

    /// Set the equilibrium position.
    pub fn set_target(&mut self, target: Vec3H) {
        self.target_pos = target;
    }

    /// Reset internal state.
    pub fn reset(&mut self) {
        self.prev_vel = Vec3H::zero();
    }

    /// Compute natural frequency of the impedance model.
    pub fn natural_frequency(&self) -> f64 {
        if self.m_d > 1e-15 {
            (self.k_d / self.m_d).sqrt()
        } else {
            0.0
        }
    }

    /// Compute damping ratio.
    pub fn damping_ratio(&self) -> f64 {
        let wn = self.natural_frequency();
        if wn > 1e-15 && self.m_d > 1e-15 {
            self.b_d / (2.0 * self.m_d * wn)
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// AdmittanceControl
// ---------------------------------------------------------------------------

/// Admittance-type haptic rendering: the device controls force, the
/// controller outputs position.
///
/// Integrates the equation: `M_d * a + B_d * v + K_d * x = F_ext`.
#[derive(Debug, Clone)]
pub struct AdmittanceControl {
    /// Virtual mass.
    pub m_d: f64,
    /// Virtual damping.
    pub b_d: f64,
    /// Virtual stiffness.
    pub k_d: f64,
    /// Current virtual position.
    pub pos: Vec3H,
    /// Current virtual velocity.
    pub vel: Vec3H,
    /// Equilibrium position.
    pub equilibrium: Vec3H,
}

impl AdmittanceControl {
    /// Create a new admittance controller.
    pub fn new(m_d: f64, b_d: f64, k_d: f64, start_pos: Vec3H) -> Self {
        Self {
            m_d,
            b_d,
            k_d,
            pos: start_pos,
            vel: Vec3H::zero(),
            equilibrium: start_pos,
        }
    }

    /// Step the admittance model given an external force and time step.
    ///
    /// Returns the new desired position for the haptic device.
    pub fn step(&mut self, f_ext: Vec3H, dt: f64) -> Vec3H {
        if self.m_d < 1e-15 {
            return self.pos;
        }
        let x_err = self.pos - self.equilibrium;
        let accel = (f_ext - self.vel * self.b_d - x_err * self.k_d) / self.m_d;
        self.vel += accel * dt;
        self.pos += self.vel * dt;
        self.pos
    }

    /// Reset state.
    pub fn reset(&mut self, pos: Vec3H) {
        self.pos = pos;
        self.vel = Vec3H::zero();
        self.equilibrium = pos;
    }

    /// Current kinetic energy of the virtual mass.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.m_d * self.vel.len_sq()
    }

    /// Current potential energy (spring).
    pub fn potential_energy(&self) -> f64 {
        let x = (self.pos - self.equilibrium).len();
        0.5 * self.k_d * x * x
    }

    /// Total mechanical energy.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
}

// ---------------------------------------------------------------------------
// StabilityAnalysis
// ---------------------------------------------------------------------------

/// Passivity-based stability analysis for haptic rendering.
///
/// Uses an energy observer / energy tank to monitor whether the haptic
/// system is passive (energy-dissipating).
#[derive(Debug, Clone)]
pub struct StabilityAnalysis {
    /// Observed energy in the system.
    pub energy_observed: f64,
    /// Energy tank capacity (upper bound for stored energy).
    pub tank_capacity: f64,
    /// Current energy in the tank.
    pub tank_energy: f64,
    /// Whether the system is currently deemed passive.
    pub is_passive: bool,
    /// Accumulated energy input.
    pub energy_input: f64,
    /// Accumulated energy output.
    pub energy_output: f64,
}

impl StabilityAnalysis {
    /// Create a new stability analyser with a given tank capacity.
    pub fn new(tank_capacity: f64) -> Self {
        Self {
            energy_observed: 0.0,
            tank_capacity,
            tank_energy: tank_capacity,
            is_passive: true,
            energy_input: 0.0,
            energy_output: 0.0,
        }
    }

    /// Update the energy observer with the power exchange for this time step.
    ///
    /// `power = F . v` where F is the force rendered and v is device velocity.
    /// Positive power = energy flowing out of the haptic device (user doing work).
    /// Negative power = energy flowing into the device (haptic generating energy).
    pub fn update(&mut self, force: Vec3H, velocity: Vec3H, dt: f64) {
        let power = force.dot(velocity);
        let energy_delta = power * dt;

        self.energy_observed += energy_delta;

        if energy_delta > 0.0 {
            // Energy entering the system (user input)
            self.energy_input += energy_delta;
            self.tank_energy = (self.tank_energy + energy_delta).min(self.tank_capacity);
        } else {
            // Energy leaving (haptic output)
            self.energy_output += -energy_delta;
            self.tank_energy += energy_delta; // reduces tank
        }

        self.is_passive = self.tank_energy >= 0.0;
    }

    /// Check the passivity condition: total energy output <= total energy input.
    pub fn check_passivity(&self) -> bool {
        self.energy_output <= self.energy_input + 1e-10
    }

    /// Compute the Z-width: ratio of maximum stiffness to minimum damping
    /// that maintains passivity at a given sample rate.
    ///
    /// For a virtual wall: `K_max = B / T` where T is the sample period.
    pub fn z_width(damping: f64, sample_period: f64) -> f64 {
        if sample_period > 1e-15 {
            damping / sample_period
        } else {
            f64::INFINITY
        }
    }

    /// Maximum stable stiffness for a sampled haptic system.
    ///
    /// `K_max = 2 * B / T - (B^2 * T) / (2 * M)`
    /// simplified for zero-order hold at sample period `T`.
    #[allow(clippy::too_many_arguments)]
    pub fn max_stable_stiffness(damping: f64, mass: f64, sample_period: f64) -> f64 {
        if sample_period < 1e-15 {
            return f64::INFINITY;
        }
        let term1 = 2.0 * damping / sample_period;
        let term2 = if mass > 1e-15 {
            damping * damping * sample_period / (2.0 * mass)
        } else {
            0.0
        };
        (term1 - term2).max(0.0)
    }

    /// Scale a force to ensure passivity (clamp if tank is empty).
    pub fn passivity_clamp(&mut self, force: Vec3H, velocity: Vec3H, dt: f64) -> Vec3H {
        let power = force.dot(velocity);
        let energy_needed = -power * dt; // energy that would leave the tank
        if energy_needed > 0.0 && energy_needed > self.tank_energy {
            // Scale down the force to not exceed available energy
            if power.abs() < 1e-15 || dt < 1e-15 {
                return force;
            }
            let scale = self.tank_energy / energy_needed;
            force * scale.max(0.0)
        } else {
            force
        }
    }

    /// Reset the analyser.
    pub fn reset(&mut self) {
        self.energy_observed = 0.0;
        self.tank_energy = self.tank_capacity;
        self.is_passive = true;
        self.energy_input = 0.0;
        self.energy_output = 0.0;
    }
}

// ---------------------------------------------------------------------------
// MultiPointHaptic
// ---------------------------------------------------------------------------

/// A single haptic contact finger.
#[derive(Debug, Clone)]
pub struct HapticFinger {
    /// Finger identifier.
    pub id: usize,
    /// Current finger position.
    pub pos: Vec3H,
    /// Proxy position (constrained).
    pub proxy_pos: Vec3H,
    /// Finger stiffness.
    pub stiffness: f64,
    /// Contact normal (if in contact).
    pub contact_normal: Vec3H,
    /// Whether currently in contact.
    pub in_contact: bool,
}

impl HapticFinger {
    /// Create a new finger.
    pub fn new(id: usize, pos: Vec3H, stiffness: f64) -> Self {
        Self {
            id,
            pos,
            proxy_pos: pos,
            stiffness,
            contact_normal: Vec3H::zero(),
            in_contact: false,
        }
    }

    /// Compute the contact force for this finger.
    pub fn compute_force(&self) -> Vec3H {
        if self.in_contact {
            (self.proxy_pos - self.pos) * self.stiffness
        } else {
            Vec3H::zero()
        }
    }
}

/// Multi-finger haptic interaction for grasp rendering.
#[derive(Debug, Clone)]
pub struct MultiPointHaptic {
    /// The fingers.
    pub fingers: Vec<HapticFinger>,
    /// Maximum total grasp force.
    pub max_grasp_force: f64,
}

impl MultiPointHaptic {
    /// Create a new multi-point haptic system.
    pub fn new(max_grasp_force: f64) -> Self {
        Self {
            fingers: Vec::new(),
            max_grasp_force,
        }
    }

    /// Add a finger.
    pub fn add_finger(&mut self, finger: HapticFinger) {
        self.fingers.push(finger);
    }

    /// Number of fingers.
    pub fn num_fingers(&self) -> usize {
        self.fingers.len()
    }

    /// Update finger positions against a sphere surface.
    pub fn update_sphere(&mut self, surface: &SurfaceSphere) {
        for finger in &mut self.fingers {
            let sd = surface.signed_distance(finger.pos);
            if sd < 0.0 {
                finger.proxy_pos = surface.project(finger.pos);
                finger.contact_normal = surface.normal_at(finger.proxy_pos);
                finger.in_contact = true;
            } else {
                finger.proxy_pos = finger.pos;
                finger.in_contact = false;
            }
        }
    }

    /// Update finger positions against a plane surface.
    pub fn update_plane(&mut self, surface: &SurfacePlane) {
        for finger in &mut self.fingers {
            let sd = surface.signed_distance(finger.pos);
            if sd < 0.0 {
                finger.proxy_pos = surface.project(finger.pos);
                finger.contact_normal = surface.normal;
                finger.in_contact = true;
            } else {
                finger.proxy_pos = finger.pos;
                finger.in_contact = false;
            }
        }
    }

    /// Compute total grasp force (sum of all finger forces).
    pub fn total_force(&self) -> Vec3H {
        let mut total = Vec3H::zero();
        for f in &self.fingers {
            total += f.compute_force();
        }
        total.clamped(self.max_grasp_force)
    }

    /// Compute grasp wrench (force + torque about a reference point).
    pub fn grasp_wrench(&self, ref_point: Vec3H) -> (Vec3H, Vec3H) {
        let mut force = Vec3H::zero();
        let mut torque = Vec3H::zero();
        for f in &self.fingers {
            let fi = f.compute_force();
            force += fi;
            let r = f.pos - ref_point;
            torque += r.cross(fi);
        }
        (force.clamped(self.max_grasp_force), torque)
    }

    /// Number of fingers currently in contact.
    pub fn contact_count(&self) -> usize {
        self.fingers.iter().filter(|f| f.in_contact).count()
    }

    /// Whether we have a stable grasp (at least 2 opposing contacts).
    pub fn is_grasp_stable(&self) -> bool {
        let contacting: Vec<&HapticFinger> = self.fingers.iter().filter(|f| f.in_contact).collect();
        if contacting.len() < 2 {
            return false;
        }
        // Check for opposing normals
        for i in 0..contacting.len() {
            for j in (i + 1)..contacting.len() {
                let dot = contacting[i]
                    .contact_normal
                    .dot(contacting[j].contact_normal);
                if dot < -0.5 {
                    return true; // opposing normals → stable
                }
            }
        }
        false
    }

    /// Compute the grasp centroid (mean of contacting finger positions).
    pub fn grasp_centroid(&self) -> Vec3H {
        let contacting: Vec<&HapticFinger> = self.fingers.iter().filter(|f| f.in_contact).collect();
        if contacting.is_empty() {
            return Vec3H::zero();
        }
        let mut sum = Vec3H::zero();
        for f in &contacting {
            sum += f.pos;
        }
        sum / contacting.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Force interpolation utility
// ---------------------------------------------------------------------------

/// Linearly interpolate between two forces.
pub fn lerp_force(a: Vec3H, b: Vec3H, t: f64) -> Vec3H {
    a * (1.0 - t) + b * t
}

/// Compute the work done by a force over a displacement.
pub fn work(force: Vec3H, displacement: Vec3H) -> f64 {
    force.dot(displacement)
}

/// Compute instantaneous power: P = F . v.
pub fn power(force: Vec3H, velocity: Vec3H) -> f64 {
    force.dot(velocity)
}

/// Compute the mechanical impedance magnitude |Z| = |F| / |v|.
pub fn impedance_magnitude(force_mag: f64, velocity_mag: f64) -> f64 {
    if velocity_mag > 1e-15 {
        force_mag / velocity_mag
    } else {
        f64::INFINITY
    }
}

/// Virtual wall: compute contact force for a 1-D virtual wall at position `wall_x`.
///
/// Returns force in the x direction.
pub fn virtual_wall_force(
    device_x: f64,
    wall_x: f64,
    stiffness: f64,
    damping: f64,
    vel_x: f64,
) -> f64 {
    let penetration = wall_x - device_x;
    if penetration > 0.0 {
        stiffness * penetration - damping * vel_x
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Vec3H ----

    #[test]
    fn vec3h_basic_ops() {
        let a = Vec3H::new(1.0, 2.0, 3.0);
        let b = Vec3H::new(4.0, 5.0, 6.0);
        let c = a + b;
        assert!((c.x - 5.0).abs() < 1e-10);
        assert!((c.y - 7.0).abs() < 1e-10);
        assert!((c.z - 9.0).abs() < 1e-10);
    }

    #[test]
    fn vec3h_dot_and_cross() {
        let a = Vec3H::new(1.0, 0.0, 0.0);
        let b = Vec3H::new(0.0, 1.0, 0.0);
        assert!(a.dot(b).abs() < 1e-10);
        let c = a.cross(b);
        assert!((c.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn vec3h_normalised() {
        let v = Vec3H::new(3.0, 4.0, 0.0);
        let n = v.normalised();
        assert!((n.len() - 1.0).abs() < 1e-10);
    }

    // ---- SurfacePlane ----

    #[test]
    fn plane_signed_distance() {
        let plane = SurfacePlane::new(Vec3H::zero(), Vec3H::new(0.0, 1.0, 0.0));
        assert!((plane.signed_distance(Vec3H::new(0.0, 5.0, 0.0)) - 5.0).abs() < 1e-10);
        assert!((plane.signed_distance(Vec3H::new(0.0, -3.0, 0.0)) + 3.0).abs() < 1e-10);
    }

    #[test]
    fn plane_project() {
        let plane = SurfacePlane::new(Vec3H::zero(), Vec3H::new(0.0, 1.0, 0.0));
        let proj = plane.project(Vec3H::new(3.0, 5.0, 7.0));
        assert!((proj.y).abs() < 1e-10);
        assert!((proj.x - 3.0).abs() < 1e-10);
    }

    // ---- SurfaceSphere ----

    #[test]
    fn sphere_signed_distance() {
        let sphere = SurfaceSphere::new(Vec3H::zero(), 2.0);
        assert!((sphere.signed_distance(Vec3H::new(5.0, 0.0, 0.0)) - 3.0).abs() < 1e-10);
        assert!((sphere.signed_distance(Vec3H::new(1.0, 0.0, 0.0)) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_project() {
        let sphere = SurfaceSphere::new(Vec3H::zero(), 3.0);
        let proj = sphere.project(Vec3H::new(6.0, 0.0, 0.0));
        assert!((proj.x - 3.0).abs() < 1e-10);
        assert!(proj.y.abs() < 1e-10);
    }

    // ---- HapticProxy ----

    #[test]
    fn proxy_spring_force() {
        let mut proxy = HapticProxy::new(Vec3H::zero(), 1000.0, 0.0);
        proxy.device_pos = Vec3H::new(0.0, -0.01, 0.0);
        proxy.proxy_pos = Vec3H::zero(); // on surface
        let f = proxy.compute_force();
        // F = K * (proxy - device) = 1000 * (0 - (-0.01)) = 10 in y
        assert!((f.y - 10.0).abs() < 1e-6, "f.y={}", f.y);
    }

    #[test]
    fn proxy_no_force_in_free_space() {
        let proxy = HapticProxy::new(Vec3H::new(1.0, 2.0, 3.0), 1000.0, 0.0);
        // proxy_pos == device_pos
        let f = proxy.compute_force();
        assert!(f.len() < 1e-10);
    }

    #[test]
    fn proxy_in_contact_detection() {
        let mut proxy = HapticProxy::new(Vec3H::zero(), 1000.0, 0.0);
        assert!(!proxy.in_contact());
        proxy.device_pos = Vec3H::new(0.0, -0.01, 0.0);
        assert!(proxy.in_contact());
    }

    // ---- ProxyMethod ----

    #[test]
    fn proxy_method_plane_contact() {
        let mut pm = ProxyMethod::new(Vec3H::new(0.0, 1.0, 0.0), 500.0, 0.0, 100.0);
        let plane = SurfacePlane::new(Vec3H::zero(), Vec3H::new(0.0, 1.0, 0.0));
        // Device penetrates the plane
        let f = pm.update_plane(Vec3H::new(0.0, -0.05, 0.0), 0.001, &plane);
        // Force should push upward (positive y)
        assert!(f.y > 0.0, "f.y={}", f.y);
    }

    #[test]
    fn proxy_method_no_force_outside() {
        let mut pm = ProxyMethod::new(Vec3H::new(0.0, 1.0, 0.0), 500.0, 0.0, 100.0);
        let plane = SurfacePlane::new(Vec3H::zero(), Vec3H::new(0.0, 1.0, 0.0));
        let f = pm.update_plane(Vec3H::new(0.0, 2.0, 0.0), 0.001, &plane);
        assert!(f.len() < 1e-10);
    }

    #[test]
    fn proxy_method_sphere_contact() {
        let mut pm = ProxyMethod::new(Vec3H::new(3.0, 0.0, 0.0), 500.0, 0.0, 1000.0);
        let sphere = SurfaceSphere::new(Vec3H::zero(), 2.0);
        let f = pm.update_sphere(Vec3H::new(1.0, 0.0, 0.0), 0.001, &sphere);
        // Device is inside sphere → force pushes outward (positive x)
        assert!(f.x > 0.0, "f.x={}", f.x);
    }

    // ---- FrictionHaptic ----

    #[test]
    fn friction_stick_phase() {
        let mut fh = FrictionHaptic::new(0.5, 0.3, 1000.0);
        fh.set_anchor(Vec3H::zero());
        let normal = Vec3H::new(0.0, 1.0, 0.0);
        // Small displacement (within static friction cone)
        let f = fh.compute_friction(Vec3H::new(0.001, 0.0, 0.0), 10.0, normal);
        assert!(fh.is_stuck);
        // Should be restoring force in -x direction
        assert!(f.x < 0.0, "f.x={}", f.x);
    }

    #[test]
    fn friction_slip_transition() {
        let mut fh = FrictionHaptic::new(0.3, 0.2, 1000.0);
        fh.set_anchor(Vec3H::zero());
        let normal = Vec3H::new(0.0, 1.0, 0.0);
        // Large displacement to exceed static friction
        // max_static = 0.3 * 10 = 3.0, spring_force = 10.0 * 1000 = 10000
        let f = fh.compute_friction(Vec3H::new(10.0, 0.0, 0.0), 10.0, normal);
        assert!(!fh.is_stuck, "should have transitioned to slip");
        // Kinetic friction force
        assert!(f.x < 0.0);
    }

    #[test]
    fn friction_zero_normal_force() {
        let mut fh = FrictionHaptic::new(0.5, 0.3, 1000.0);
        fh.set_anchor(Vec3H::zero());
        let normal = Vec3H::new(0.0, 1.0, 0.0);
        // Even with large displacement, zero normal force means tiny friction threshold
        let _f = fh.compute_friction(Vec3H::new(5.0, 0.0, 0.0), 0.0, normal);
        // With zero normal force, static limit=0, so it slips immediately
        assert!(!fh.is_stuck);
    }

    // ---- HapticTexture ----

    #[test]
    fn texture_height_sinusoidal() {
        let tex = HapticTexture::new(1.0, 1.0, 0.001, 1.0);
        let h0 = tex.height(0.0, 0.0);
        assert!(h0.abs() < 1e-10); // sin(0) = 0
        let h_quarter = tex.height(0.25, 0.0);
        // sin(PI/2) = 1, so h = 0.001 * 1 * 0.5 = 0.0005
        assert!((h_quarter - 0.0005).abs() < 1e-6, "h={h_quarter}");
    }

    #[test]
    fn texture_wavelength() {
        let tex = HapticTexture::new(10.0, 5.0, 0.001, 1.0);
        assert!((tex.wavelength_x() - 0.1).abs() < 1e-10);
        assert!((tex.wavelength_y() - 0.2).abs() < 1e-10);
    }

    // ---- ImpedanceControl ----

    #[test]
    fn impedance_force_at_equilibrium() {
        let mut ic = ImpedanceControl::new(0.0, 0.0, 100.0, Vec3H::zero());
        let f = ic.compute_force(Vec3H::zero(), Vec3H::zero(), 0.001);
        assert!(f.len() < 1e-10);
    }

    #[test]
    fn impedance_spring_force() {
        let mut ic = ImpedanceControl::new(0.0, 0.0, 100.0, Vec3H::zero());
        let f = ic.compute_force(Vec3H::new(0.1, 0.0, 0.0), Vec3H::zero(), 0.001);
        // F = -K * x = -100 * 0.1 = -10
        assert!((f.x + 10.0).abs() < 1e-6, "f.x={}", f.x);
    }

    #[test]
    fn impedance_damping_force() {
        let mut ic = ImpedanceControl::new(0.0, 10.0, 0.0, Vec3H::zero());
        let f = ic.compute_force(Vec3H::zero(), Vec3H::new(1.0, 0.0, 0.0), 0.001);
        // F = -B * v = -10 * 1 = -10
        assert!((f.x + 10.0).abs() < 1e-6, "f.x={}", f.x);
    }

    #[test]
    fn impedance_natural_frequency() {
        let ic = ImpedanceControl::new(1.0, 0.0, 100.0, Vec3H::zero());
        assert!((ic.natural_frequency() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn impedance_damping_ratio() {
        // Critical damping: zeta = B / (2*sqrt(K*M)) = 20 / (2*10) = 1.0
        let ic = ImpedanceControl::new(1.0, 20.0, 100.0, Vec3H::zero());
        assert!((ic.damping_ratio() - 1.0).abs() < 1e-10);
    }

    // ---- AdmittanceControl ----

    #[test]
    fn admittance_stationary_no_force() {
        let mut ac = AdmittanceControl::new(1.0, 10.0, 100.0, Vec3H::zero());
        let pos = ac.step(Vec3H::zero(), 0.001);
        assert!(pos.len() < 1e-6);
    }

    #[test]
    fn admittance_force_causes_motion() {
        let mut ac = AdmittanceControl::new(1.0, 0.0, 0.0, Vec3H::zero());
        let _pos = ac.step(Vec3H::new(10.0, 0.0, 0.0), 0.01);
        // After one step: accel = 10/1 = 10, vel = 10*0.01 = 0.1, pos = 0.1*0.01 = 0.001
        assert!(ac.vel.x > 0.0);
        assert!(ac.pos.x > 0.0);
    }

    #[test]
    fn admittance_energy_conservation_like() {
        let mut ac = AdmittanceControl::new(1.0, 0.0, 100.0, Vec3H::zero());
        ac.pos = Vec3H::new(0.1, 0.0, 0.0); // displaced
        let e0 = ac.total_energy();
        // Step with no external force → energy should be conserved (approximately)
        for _ in 0..10 {
            ac.step(Vec3H::zero(), 0.0001);
        }
        let e1 = ac.total_energy();
        // With no damping, energy should be approximately conserved for small dt
        assert!((e1 - e0).abs() / (e0 + 1e-15) < 0.01, "e0={e0}, e1={e1}");
    }

    // ---- StabilityAnalysis ----

    #[test]
    fn passivity_condition_holds_for_passive_system() {
        let mut sa = StabilityAnalysis::new(10.0);
        // User pushes device (positive work → energy input)
        let f = Vec3H::new(5.0, 0.0, 0.0);
        let v = Vec3H::new(1.0, 0.0, 0.0);
        sa.update(f, v, 0.001);
        assert!(sa.check_passivity());
        assert!(sa.is_passive);
    }

    #[test]
    fn passivity_violated_when_generating_energy() {
        let mut sa = StabilityAnalysis::new(0.001); // tiny tank
        // Force opposes velocity → system generates energy
        let f = Vec3H::new(-100.0, 0.0, 0.0);
        let v = Vec3H::new(1.0, 0.0, 0.0);
        // Multiple steps to drain the tank
        for _ in 0..100 {
            sa.update(f, v, 0.001);
        }
        assert!(!sa.is_passive, "tank should be drained");
    }

    #[test]
    fn z_width_computation() {
        let zw = StabilityAnalysis::z_width(10.0, 0.001);
        assert!((zw - 10000.0).abs() < 1e-6);
    }

    #[test]
    fn max_stable_stiffness_basic() {
        let k_max = StabilityAnalysis::max_stable_stiffness(10.0, 1.0, 0.001);
        // term1 = 2*10/0.001 = 20000
        // term2 = 100 * 0.001 / 2 = 0.05
        assert!((k_max - 19999.95).abs() < 0.1, "k_max={k_max}");
    }

    #[test]
    fn passivity_clamp_reduces_force() {
        let mut sa = StabilityAnalysis::new(0.001);
        // Force generating large negative power
        let f = Vec3H::new(-1000.0, 0.0, 0.0);
        let v = Vec3H::new(1.0, 0.0, 0.0);
        let clamped = sa.passivity_clamp(f, v, 1.0);
        // The clamped force should be much smaller than the original
        assert!(clamped.len() < f.len());
    }

    // ---- MultiPointHaptic ----

    #[test]
    fn multipoint_no_contact() {
        let mut mph = MultiPointHaptic::new(100.0);
        mph.add_finger(HapticFinger::new(0, Vec3H::new(5.0, 0.0, 0.0), 500.0));
        let sphere = SurfaceSphere::new(Vec3H::zero(), 2.0);
        mph.update_sphere(&sphere);
        assert_eq!(mph.contact_count(), 0);
        assert!(mph.total_force().len() < 1e-10);
    }

    #[test]
    fn multipoint_single_contact() {
        let mut mph = MultiPointHaptic::new(1000.0);
        mph.add_finger(HapticFinger::new(0, Vec3H::new(1.0, 0.0, 0.0), 500.0));
        let sphere = SurfaceSphere::new(Vec3H::zero(), 2.0);
        mph.update_sphere(&sphere);
        assert_eq!(mph.contact_count(), 1);
        let f = mph.total_force();
        assert!(f.x > 0.0, "should push outward");
    }

    #[test]
    fn multipoint_grasp_stability() {
        let mut mph = MultiPointHaptic::new(1000.0);
        // Two fingers on opposite sides of a sphere
        mph.add_finger(HapticFinger::new(0, Vec3H::new(1.0, 0.0, 0.0), 500.0));
        mph.add_finger(HapticFinger::new(1, Vec3H::new(-1.0, 0.0, 0.0), 500.0));
        let sphere = SurfaceSphere::new(Vec3H::zero(), 2.0);
        mph.update_sphere(&sphere);
        assert_eq!(mph.contact_count(), 2);
        assert!(mph.is_grasp_stable());
    }

    #[test]
    fn multipoint_grasp_centroid() {
        let mut mph = MultiPointHaptic::new(1000.0);
        mph.add_finger(HapticFinger::new(0, Vec3H::new(1.0, 0.0, 0.0), 500.0));
        mph.add_finger(HapticFinger::new(1, Vec3H::new(-1.0, 0.0, 0.0), 500.0));
        let sphere = SurfaceSphere::new(Vec3H::zero(), 3.0);
        mph.update_sphere(&sphere);
        let c = mph.grasp_centroid();
        assert!(c.len() < 1e-10, "centroid should be at origin");
    }

    // ---- Utility functions ----

    #[test]
    fn work_computation() {
        let f = Vec3H::new(10.0, 0.0, 0.0);
        let d = Vec3H::new(2.0, 0.0, 0.0);
        assert!((work(f, d) - 20.0).abs() < 1e-10);
    }

    #[test]
    fn power_computation() {
        let f = Vec3H::new(5.0, 0.0, 0.0);
        let v = Vec3H::new(3.0, 0.0, 0.0);
        assert!((power(f, v) - 15.0).abs() < 1e-10);
    }

    #[test]
    fn virtual_wall_no_contact() {
        let f = virtual_wall_force(1.0, 0.0, 1000.0, 10.0, 0.0);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn virtual_wall_contact() {
        let f = virtual_wall_force(-0.01, 0.0, 1000.0, 0.0, 0.0);
        // penetration = 0.0 - (-0.01) = 0.01 > 0, F = 1000 * 0.01 = 10
        assert!((f - 10.0).abs() < 1e-6, "f={f}");
    }

    #[test]
    fn impedance_magnitude_basic() {
        let z = impedance_magnitude(10.0, 2.0);
        assert!((z - 5.0).abs() < 1e-10);
    }

    #[test]
    fn lerp_force_midpoint() {
        let a = Vec3H::new(0.0, 0.0, 0.0);
        let b = Vec3H::new(10.0, 0.0, 0.0);
        let mid = lerp_force(a, b, 0.5);
        assert!((mid.x - 5.0).abs() < 1e-10);
    }
}
