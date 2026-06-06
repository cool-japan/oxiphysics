// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular matter simulation using Discrete Element Method (DEM).
//!
//! Implements sphere-sphere and sphere-wall contacts with Hertz-Mindlin
//! contact mechanics, angular integration, and granular flow analysis.

use std::f64::consts::PI;

// ── vector helpers (no nalgebra) ──────────────────────────────────────────────

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n > 1e-15 {
        vec3_scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_zeros() -> [f64; 3] {
    [0.0, 0.0, 0.0]
}

// ── GranularParticle ──────────────────────────────────────────────────────────

/// A single spherical granular particle for DEM simulation.
///
/// Stores position, linear and angular velocity, radius, and mass.
/// Forces and torques are accumulated each step and then cleared.
#[derive(Debug, Clone)]
pub struct GranularParticle {
    /// World-space position \[m\].
    pub position: [f64; 3],
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Accumulated force this step \[N\].
    pub force: [f64; 3],
    /// Accumulated torque this step \[N·m\].
    pub torque: [f64; 3],
}

impl GranularParticle {
    /// Create a new particle with given position, radius and density.
    ///
    /// Mass is computed from `4/3 π r³ ρ`. Initial velocities are zero.
    pub fn new(position: [f64; 3], radius: f64, density: f64) -> Self {
        let mass = (4.0 / 3.0) * PI * radius.powi(3) * density;
        Self {
            position,
            velocity: vec3_zeros(),
            angular_velocity: vec3_zeros(),
            radius,
            mass,
            force: vec3_zeros(),
            torque: vec3_zeros(),
        }
    }

    /// Moment of inertia for a solid sphere: `2/5 m r²`.
    pub fn inertia(&self) -> f64 {
        (2.0 / 5.0) * self.mass * self.radius * self.radius
    }

    /// Translational kinetic energy `½ m v²`.
    pub fn kinetic_energy_translational(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }

    /// Rotational kinetic energy `½ I ω²`.
    pub fn kinetic_energy_rotational(&self) -> f64 {
        0.5 * self.inertia() * vec3_dot(self.angular_velocity, self.angular_velocity)
    }

    /// Clear accumulated force and torque.
    pub fn clear_accumulators(&mut self) {
        self.force = vec3_zeros();
        self.torque = vec3_zeros();
    }

    /// Apply gravity `[0, -g, 0]`.
    pub fn apply_gravity(&mut self, g: f64) {
        self.force[1] -= self.mass * g;
    }
}

// ── DemContact ────────────────────────────────────────────────────────────────

/// Hertz-Mindlin contact parameters for sphere-sphere interaction.
///
/// Normal force uses Hertz contact (`kn * δ^1.5 - cn * δ̇`),
/// tangential uses linear spring-dashpot with Coulomb friction limit.
#[derive(Debug, Clone)]
pub struct DemContact {
    /// Normal spring stiffness \[N/m^1.5\].
    pub kn: f64,
    /// Normal damping coefficient \[N·s/m\].
    pub cn: f64,
    /// Tangential spring stiffness \[N/m\].
    pub kt: f64,
    /// Tangential damping coefficient \[N·s/m\].
    pub ct: f64,
    /// Coulomb friction coefficient (static ≈ dynamic for simplicity).
    pub mu: f64,
}

impl DemContact {
    /// Create a contact model for glass beads (typical DEM values).
    pub fn glass_beads() -> Self {
        Self {
            kn: 1.0e6,
            cn: 50.0,
            kt: 7.0e5,
            ct: 30.0,
            mu: 0.3,
        }
    }

    /// Create a contact model for sand-like particles.
    pub fn sand() -> Self {
        Self {
            kn: 5.0e5,
            cn: 30.0,
            kt: 3.5e5,
            ct: 20.0,
            mu: 0.5,
        }
    }

    /// Compute and apply contact forces between two particles.
    ///
    /// Updates `force` and `torque` accumulators on both particles.
    pub fn apply_contact(&self, a: &mut GranularParticle, b: &mut GranularParticle) {
        let delta_pos = vec3_sub(b.position, a.position);
        let dist = vec3_norm(delta_pos);
        let overlap = a.radius + b.radius - dist;
        if overlap <= 0.0 {
            return;
        }

        let n = vec3_normalize(delta_pos); // unit normal a → b

        // Relative velocity at contact point
        let r_a = vec3_scale(n, a.radius);
        let r_b = vec3_scale(n, -b.radius);
        let v_contact_a = vec3_add(a.velocity, vec3_cross(a.angular_velocity, r_a));
        let v_contact_b = vec3_add(b.velocity, vec3_cross(b.angular_velocity, r_b));
        let v_rel = vec3_sub(v_contact_a, v_contact_b);

        // Normal component of relative velocity
        let v_n_scalar = vec3_dot(v_rel, n);
        let v_n = vec3_scale(n, v_n_scalar);

        // Hertz normal force (overlap^1.5)
        // v_n_scalar > 0 means particles approaching (closing velocity along n: a→b)
        let fn_mag = self.kn * overlap.powf(1.5) - self.cn * v_n_scalar;
        let fn_mag = fn_mag.max(0.0); // no tensile force
        // Repulsive: a is pushed in -n direction (away from b)
        let f_normal_on_a = vec3_scale(n, -fn_mag);

        // Tangential relative velocity
        let v_t = vec3_sub(v_rel, v_n);
        let f_t_raw = vec3_sub(
            vec3_scale(v_t, -self.ct),
            vec3_scale(vec3_normalize(v_t), self.kt * overlap),
        );
        // Coulomb limit
        let ft_max = self.mu * fn_mag;
        let ft_norm = vec3_norm(f_t_raw);
        let f_tangential = if ft_norm > ft_max && ft_norm > 1e-15 {
            vec3_scale(f_t_raw, ft_max / ft_norm)
        } else {
            f_t_raw
        };

        let f_total_on_a = vec3_add(f_normal_on_a, f_tangential);

        // Apply equal and opposite forces
        a.force = vec3_add(a.force, f_total_on_a);
        b.force = vec3_sub(b.force, f_total_on_a);

        // Torques: τ = r_contact × F_tangential
        // r_a = radius * n (vector from a centre to contact point)
        let torque_a = vec3_cross(r_a, f_tangential);
        // force on b is -f_tangential; r_b points from b centre to contact = -n * r_b
        let torque_b = vec3_cross(r_b, vec3_scale(f_tangential, -1.0));
        a.torque = vec3_add(a.torque, torque_a);
        b.torque = vec3_add(b.torque, torque_b);
    }
}

// ── AngularIntegration ────────────────────────────────────────────────────────

/// Utilities for integrating angular velocity and orientation from torque.
pub struct AngularIntegration;

impl AngularIntegration {
    /// Update angular velocity using Euler integration.
    ///
    /// `ω_new = ω + (τ / I) * dt`
    pub fn integrate(particle: &mut GranularParticle, dt: f64) {
        let alpha = vec3_scale(particle.torque, 1.0 / particle.inertia());
        particle.angular_velocity = vec3_add(particle.angular_velocity, vec3_scale(alpha, dt));
    }

    /// Apply angular damping: `ω *= exp(-γ dt)`.
    pub fn damp(particle: &mut GranularParticle, damping: f64, dt: f64) {
        let factor = (-damping * dt).exp();
        particle.angular_velocity = vec3_scale(particle.angular_velocity, factor);
    }
}

// ── WallBoundary ──────────────────────────────────────────────────────────────

/// An infinite planar wall for particle-wall contact.
///
/// The wall is defined by a point on the plane and an inward normal.
#[derive(Debug, Clone)]
pub struct WallBoundary {
    /// A point on the wall plane.
    pub point: [f64; 3],
    /// Inward normal (pointing into the simulation domain).
    pub normal: [f64; 3],
    /// Contact stiffness \[N/m^1.5\].
    pub kn: f64,
    /// Contact damping \[N·s/m\].
    pub cn: f64,
    /// Wall friction coefficient.
    pub mu: f64,
}

impl WallBoundary {
    /// Create a floor wall at `y = y0`.
    pub fn floor(y0: f64) -> Self {
        Self {
            point: [0.0, y0, 0.0],
            normal: [0.0, 1.0, 0.0],
            kn: 1.0e6,
            cn: 50.0,
            mu: 0.3,
        }
    }

    /// Signed distance from particle centre to wall (positive = inside domain).
    pub fn signed_distance(&self, pos: [f64; 3]) -> f64 {
        vec3_dot(vec3_sub(pos, self.point), self.normal)
    }

    /// Apply wall contact force to `particle` if it overlaps.
    pub fn apply_contact(&self, particle: &mut GranularParticle) {
        let dist = self.signed_distance(particle.position);
        let overlap = particle.radius - dist;
        if overlap <= 0.0 {
            return;
        }
        let n = self.normal;
        let v_n_scalar = vec3_dot(particle.velocity, n);
        let fn_mag = (self.kn * overlap.powf(1.5) - self.cn * v_n_scalar).max(0.0);
        let f_normal = vec3_scale(n, fn_mag);

        // Tangential velocity component at contact
        let v_t = vec3_sub(particle.velocity, vec3_scale(n, v_n_scalar));
        let ft_limit = self.mu * fn_mag;
        let vt_norm = vec3_norm(v_t);
        let f_tangential = if vt_norm > 1e-15 {
            let ft_mag = (1000.0 * vt_norm).min(ft_limit); // viscous wall friction
            vec3_scale(vec3_normalize(v_t), -ft_mag)
        } else {
            vec3_zeros()
        };

        particle.force = vec3_add(particle.force, vec3_add(f_normal, f_tangential));

        // Torque from tangential wall force
        let r = vec3_scale(n, -particle.radius);
        particle.torque = vec3_add(particle.torque, vec3_cross(r, f_tangential));
    }
}

// ── DemSimulation ─────────────────────────────────────────────────────────────

/// A full DEM simulation holding particles, walls, and contact parameters.
///
/// Each call to [`DemSimulation::step`] performs:
/// 1. Clear force accumulators
/// 2. Apply gravity
/// 3. Particle-particle contacts (brute force O(n²) for simplicity)
/// 4. Particle-wall contacts
/// 5. Semi-implicit Euler integration
#[derive(Debug)]
pub struct DemSimulation {
    /// All particles in the simulation.
    pub particles: Vec<GranularParticle>,
    /// Wall boundaries.
    pub walls: Vec<WallBoundary>,
    /// Contact law for particle-particle interactions.
    pub contact: DemContact,
    /// Gravitational acceleration \[m/s²\].
    pub gravity: f64,
    /// Simulation time \[s\].
    pub time: f64,
}

impl DemSimulation {
    /// Create a new empty simulation.
    pub fn new(contact: DemContact, gravity: f64) -> Self {
        Self {
            particles: Vec::new(),
            walls: Vec::new(),
            contact,
            gravity,
            time: 0.0,
        }
    }

    /// Add a particle and return its index.
    pub fn add_particle(&mut self, p: GranularParticle) -> usize {
        let idx = self.particles.len();
        self.particles.push(p);
        idx
    }

    /// Add a wall boundary.
    pub fn add_wall(&mut self, w: WallBoundary) {
        self.walls.push(w);
    }

    /// Advance the simulation by one time step `dt`.
    pub fn step(&mut self, dt: f64) {
        // 1. Clear accumulators
        for p in &mut self.particles {
            p.clear_accumulators();
            p.apply_gravity(self.gravity);
        }

        // 2. Particle-particle contacts (O(n²) neighbor check)
        let n = self.particles.len();
        for i in 0..n {
            for j in (i + 1)..n {
                // Safety: split borrow
                let (left, right) = self.particles.split_at_mut(j);
                self.contact.apply_contact(&mut left[i], &mut right[0]);
            }
        }

        // 3. Particle-wall contacts
        let walls: Vec<WallBoundary> = self.walls.clone();
        for p in &mut self.particles {
            for w in &walls {
                w.apply_contact(p);
            }
        }

        // 4. Integrate
        for p in &mut self.particles {
            // Linear
            let a = vec3_scale(p.force, 1.0 / p.mass);
            p.velocity = vec3_add(p.velocity, vec3_scale(a, dt));
            p.position = vec3_add(p.position, vec3_scale(p.velocity, dt));
            // Angular
            AngularIntegration::integrate(p, dt);
        }

        self.time += dt;
    }

    /// Return the total kinetic energy (translational + rotational) of all particles.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.kinetic_energy_translational() + p.kinetic_energy_rotational())
            .sum()
    }

    /// Return the number of particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }
}

// ── GranularFlow ──────────────────────────────────────────────────────────────

/// Measures granular flow rate through a planar cross-section.
///
/// The cross-section is defined by a signed-distance plane; particles
/// that cross from negative to positive side are counted.
#[derive(Debug, Clone)]
pub struct GranularFlow {
    /// Plane point defining the cross-section.
    pub plane_point: [f64; 3],
    /// Plane normal (flow direction).
    pub plane_normal: [f64; 3],
    /// Number of particles that have crossed.
    pub crossing_count: usize,
    /// Previous signed distances per particle (keyed by particle index).
    prev_distances: Vec<Option<f64>>,
}

impl GranularFlow {
    /// Create a new flow counter for the given plane.
    pub fn new(plane_point: [f64; 3], plane_normal: [f64; 3]) -> Self {
        Self {
            plane_point,
            plane_normal,
            crossing_count: 0,
            prev_distances: Vec::new(),
        }
    }

    /// Update crossing count given current particle positions.
    ///
    /// Returns the number of new crossings detected this step.
    pub fn update(&mut self, particles: &[GranularParticle]) -> usize {
        // Grow buffer if needed
        while self.prev_distances.len() < particles.len() {
            self.prev_distances.push(None);
        }
        let mut new_crossings = 0usize;
        for (i, p) in particles.iter().enumerate() {
            let d = vec3_dot(vec3_sub(p.position, self.plane_point), self.plane_normal);
            if let Some(prev_d) = self.prev_distances[i]
                && prev_d < 0.0
                && d >= 0.0
            {
                new_crossings += 1;
                self.crossing_count += 1;
            }
            self.prev_distances[i] = Some(d);
        }
        new_crossings
    }

    /// Compute average velocity along the normal direction for all particles.
    pub fn mean_normal_velocity(&self, particles: &[GranularParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let sum: f64 = particles
            .iter()
            .map(|p| vec3_dot(p.velocity, self.plane_normal))
            .sum();
        sum / particles.len() as f64
    }
}

// ── AngleOfRepose ─────────────────────────────────────────────────────────────

/// Compute the angle of repose from a granular pile.
///
/// The angle is estimated from pile geometry: `arctan(height / half_base)`.
#[derive(Debug, Clone)]
pub struct AngleOfRepose {
    /// Measured pile height \[m\].
    pub height: f64,
    /// Measured pile half-base width \[m\].
    pub half_base: f64,
}

impl AngleOfRepose {
    /// Create from measured pile dimensions.
    pub fn new(height: f64, half_base: f64) -> Self {
        Self { height, half_base }
    }

    /// Compute angle of repose in radians.
    pub fn angle_rad(&self) -> f64 {
        (self.height / self.half_base.max(1e-15)).atan()
    }

    /// Compute angle of repose in degrees.
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad().to_degrees()
    }

    /// Estimate from a list of particle positions above a floor at `y_floor`.
    ///
    /// Finds the highest point and the horizontal extent; returns the angle.
    pub fn from_particles(particles: &[GranularParticle], y_floor: f64) -> Self {
        if particles.is_empty() {
            return Self::new(0.0, 1.0);
        }
        let height = particles
            .iter()
            .map(|p| p.position[1] - y_floor)
            .fold(f64::NEG_INFINITY, f64::max);
        let x_min = particles
            .iter()
            .map(|p| p.position[0])
            .fold(f64::INFINITY, f64::min);
        let x_max = particles
            .iter()
            .map(|p| p.position[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let half_base = ((x_max - x_min) / 2.0).max(1e-15);
        Self::new(height.max(0.0), half_base)
    }
}

// ── PackingDensity ────────────────────────────────────────────────────────────

/// Random loose packing analysis: volume fraction and coordination number.
#[derive(Debug, Clone)]
pub struct PackingDensity {
    /// Domain volume \[m³\].
    pub domain_volume: f64,
}

impl PackingDensity {
    /// Create for a given domain volume.
    pub fn new(domain_volume: f64) -> Self {
        Self { domain_volume }
    }

    /// Compute volume fraction (solid volume / domain volume).
    pub fn volume_fraction(&self, particles: &[GranularParticle]) -> f64 {
        let solid: f64 = particles
            .iter()
            .map(|p| (4.0 / 3.0) * PI * p.radius.powi(3))
            .sum();
        solid / self.domain_volume.max(1e-30)
    }

    /// Compute mean coordination number (average contacts per particle).
    ///
    /// Two particles are in contact if their centre distance ≤ sum of radii + `tol`.
    pub fn coordination_number(&self, particles: &[GranularParticle], tol: f64) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let n = particles.len();
        let mut total_contacts = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = vec3_norm(vec3_sub(particles[i].position, particles[j].position));
                if d < particles[i].radius + particles[j].radius + tol {
                    total_contacts += 2; // count for both i and j
                }
            }
        }
        total_contacts as f64 / n as f64
    }
}

// ── AvalancheDynamics ─────────────────────────────────────────────────────────

/// Detects and tracks avalanche events in a granular assembly.
///
/// An avalanche is triggered when the mean kinetic energy exceeds
/// `energy_threshold` and ends when it drops below `calm_threshold`.
#[derive(Debug, Clone)]
pub struct AvalancheDynamics {
    /// Energy threshold for avalanche onset \[J\].
    pub energy_threshold: f64,
    /// Energy threshold for avalanche cessation \[J\].
    pub calm_threshold: f64,
    /// Whether an avalanche is currently active.
    pub in_avalanche: bool,
    /// Number of detected avalanche events.
    pub event_count: usize,
    /// Sizes (particle count) of completed avalanches.
    pub event_sizes: Vec<usize>,
    /// Maximum energy recorded during current event.
    current_peak_energy: f64,
}

impl AvalancheDynamics {
    /// Create with given onset / calm thresholds.
    pub fn new(energy_threshold: f64, calm_threshold: f64) -> Self {
        Self {
            energy_threshold,
            calm_threshold,
            in_avalanche: false,
            event_count: 0,
            event_sizes: Vec::new(),
            current_peak_energy: 0.0,
        }
    }

    /// Update avalanche state given current particle system.
    ///
    /// Returns `true` if a new avalanche event just started.
    pub fn update(&mut self, particles: &[GranularParticle]) -> bool {
        let ke: f64 = particles
            .iter()
            .map(|p| p.kinetic_energy_translational())
            .sum();
        let mut new_event = false;
        if !self.in_avalanche && ke > self.energy_threshold {
            self.in_avalanche = true;
            self.current_peak_energy = ke;
            new_event = true;
        } else if self.in_avalanche {
            if ke > self.current_peak_energy {
                self.current_peak_energy = ke;
            }
            if ke < self.calm_threshold {
                self.in_avalanche = false;
                self.event_count += 1;
                // Rough size proxy: particles with KE > calm_threshold / n
                let n_active = particles
                    .iter()
                    .filter(|p| p.kinetic_energy_translational() > self.calm_threshold * 0.1)
                    .count();
                self.event_sizes.push(n_active);
                self.current_peak_energy = 0.0;
            }
        }
        new_event
    }

    /// Return mean avalanche size across all completed events.
    pub fn mean_event_size(&self) -> f64 {
        if self.event_sizes.is_empty() {
            return 0.0;
        }
        self.event_sizes.iter().sum::<usize>() as f64 / self.event_sizes.len() as f64
    }
}

// ── GranularTemperature ───────────────────────────────────────────────────────

/// Granular temperature: measure of velocity fluctuations around mean flow.
///
/// Translational: `T_t = (1/3) * <(v - `v`)²>`
/// Rotational: `T_r = (1/3) * <(ω - <ω>)²> * I/m`
#[derive(Debug, Clone)]
pub struct GranularTemperature;

impl GranularTemperature {
    /// Compute translational granular temperature \[m²/s²\].
    pub fn translational(particles: &[GranularParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let n = particles.len() as f64;
        let mean_v = particles
            .iter()
            .fold(vec3_zeros(), |acc, p| vec3_add(acc, p.velocity));
        let mean_v = vec3_scale(mean_v, 1.0 / n);
        let var: f64 = particles
            .iter()
            .map(|p| {
                let dv = vec3_sub(p.velocity, mean_v);
                vec3_dot(dv, dv)
            })
            .sum::<f64>()
            / n;
        var / 3.0
    }

    /// Compute rotational granular temperature \[m²/s²\] (scaled by I/m = 2r²/5).
    pub fn rotational(particles: &[GranularParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let n = particles.len() as f64;
        let mean_w = particles
            .iter()
            .fold(vec3_zeros(), |acc, p| vec3_add(acc, p.angular_velocity));
        let mean_w = vec3_scale(mean_w, 1.0 / n);
        let var: f64 = particles
            .iter()
            .map(|p| {
                let dw = vec3_sub(p.angular_velocity, mean_w);
                let im_ratio = (2.0 / 5.0) * p.radius * p.radius; // I/m
                im_ratio * vec3_dot(dw, dw)
            })
            .sum::<f64>()
            / n;
        var / 3.0
    }

    /// Total granular temperature (translational + rotational).
    pub fn total(particles: &[GranularParticle]) -> f64 {
        Self::translational(particles) + Self::rotational(particles)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 3], r: f64) -> GranularParticle {
        GranularParticle::new(pos, r, 1000.0)
    }

    // ── GranularParticle ──────────────────────────────────────────────────

    #[test]
    fn particle_mass_from_density() {
        let p = make_particle([0.0; 3], 0.01);
        let expected = (4.0 / 3.0) * PI * 0.01_f64.powi(3) * 1000.0;
        assert!((p.mass - expected).abs() < 1e-15);
    }

    #[test]
    fn particle_inertia_solid_sphere() {
        let p = make_particle([0.0; 3], 0.01);
        let expected = (2.0 / 5.0) * p.mass * 0.01_f64.powi(2);
        // Tolerance relative to the magnitude (~1.7e-7); use 1e-20 to allow FP rounding
        assert!((p.inertia() - expected).abs() < 1e-20);
    }

    #[test]
    fn particle_kinetic_energy_translational() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.velocity = [1.0, 0.0, 0.0];
        let ke = p.kinetic_energy_translational();
        assert!((ke - 0.5 * p.mass).abs() < 1e-12);
    }

    #[test]
    fn particle_kinetic_energy_rotational() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.angular_velocity = [0.0, 1.0, 0.0];
        let ke = p.kinetic_energy_rotational();
        assert!(ke > 0.0);
    }

    #[test]
    fn particle_clear_accumulators() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.force = [1.0, 2.0, 3.0];
        p.torque = [4.0, 5.0, 6.0];
        p.clear_accumulators();
        assert_eq!(p.force, [0.0; 3]);
        assert_eq!(p.torque, [0.0; 3]);
    }

    #[test]
    fn particle_apply_gravity_adds_downward_force() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.apply_gravity(9.81);
        assert!(p.force[1] < 0.0);
        assert!((p.force[1] + p.mass * 9.81).abs() < 1e-15);
    }

    // ── DemContact ────────────────────────────────────────────────────────

    #[test]
    fn contact_no_overlap_no_force() {
        let contact = DemContact::glass_beads();
        let mut a = make_particle([0.0, 0.0, 0.0], 0.01);
        let mut b = make_particle([0.1, 0.0, 0.0], 0.01); // far apart
        contact.apply_contact(&mut a, &mut b);
        assert_eq!(a.force, [0.0; 3]);
        assert_eq!(b.force, [0.0; 3]);
    }

    #[test]
    fn contact_overlap_generates_repulsive_force() {
        let contact = DemContact::glass_beads();
        let mut a = make_particle([0.0, 0.0, 0.0], 0.01);
        let mut b = make_particle([0.015, 0.0, 0.0], 0.01); // overlap = 0.005
        contact.apply_contact(&mut a, &mut b);
        // a should be pushed in -x, b in +x
        assert!(a.force[0] < 0.0);
        assert!(b.force[0] > 0.0);
    }

    #[test]
    fn contact_forces_equal_and_opposite() {
        let contact = DemContact::sand();
        let mut a = make_particle([0.0, 0.0, 0.0], 0.01);
        let mut b = make_particle([0.015, 0.0, 0.0], 0.01);
        contact.apply_contact(&mut a, &mut b);
        for i in 0..3 {
            assert!((a.force[i] + b.force[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn contact_sand_parameters() {
        let c = DemContact::sand();
        assert!(c.mu > 0.0);
        assert!(c.kn > 0.0);
    }

    // ── WallBoundary ──────────────────────────────────────────────────────

    #[test]
    fn wall_signed_distance_above_floor() {
        let w = WallBoundary::floor(0.0);
        let p = make_particle([0.0, 0.5, 0.0], 0.01);
        assert!((w.signed_distance(p.position) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn wall_signed_distance_below_floor() {
        let w = WallBoundary::floor(0.0);
        let p = make_particle([0.0, -0.1, 0.0], 0.01);
        assert!(w.signed_distance(p.position) < 0.0);
    }

    #[test]
    fn wall_contact_particle_above_no_force() {
        let w = WallBoundary::floor(0.0);
        let mut p = make_particle([0.0, 0.5, 0.0], 0.01);
        w.apply_contact(&mut p);
        assert_eq!(p.force, [0.0; 3]);
    }

    #[test]
    fn wall_contact_particle_penetrating_gets_pushed_up() {
        let w = WallBoundary::floor(0.0);
        let mut p = make_particle([0.0, 0.005, 0.0], 0.01); // overlap = 0.005
        w.apply_contact(&mut p);
        assert!(p.force[1] > 0.0, "should push upward: {:?}", p.force);
    }

    // ── AngularIntegration ────────────────────────────────────────────────

    #[test]
    fn angular_integration_updates_omega() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.torque = [0.0, 1e-6, 0.0];
        let omega_before = p.angular_velocity[1];
        AngularIntegration::integrate(&mut p, 0.01);
        assert!(p.angular_velocity[1] > omega_before);
    }

    #[test]
    fn angular_damping_reduces_omega() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.angular_velocity = [0.0, 10.0, 0.0];
        AngularIntegration::damp(&mut p, 5.0, 0.1);
        assert!(p.angular_velocity[1] < 10.0);
    }

    // ── DemSimulation ─────────────────────────────────────────────────────

    #[test]
    fn simulation_step_applies_gravity() {
        let mut sim = DemSimulation::new(DemContact::glass_beads(), 9.81);
        sim.add_particle(make_particle([0.0, 1.0, 0.0], 0.01));
        let y0 = sim.particles[0].position[1];
        sim.step(0.01);
        let y1 = sim.particles[0].position[1];
        assert!(y1 < y0, "particle should fall under gravity");
    }

    #[test]
    fn simulation_particle_count() {
        let mut sim = DemSimulation::new(DemContact::glass_beads(), 9.81);
        sim.add_particle(make_particle([0.0, 0.0, 0.0], 0.01));
        sim.add_particle(make_particle([0.1, 0.0, 0.0], 0.01));
        assert_eq!(sim.particle_count(), 2);
    }

    #[test]
    fn simulation_wall_stops_falling_particle() {
        let mut sim = DemSimulation::new(DemContact::glass_beads(), 9.81);
        sim.add_wall(WallBoundary::floor(0.0));
        sim.add_particle(make_particle([0.0, 0.05, 0.0], 0.01));
        // Step enough to hit the floor
        for _ in 0..500 {
            sim.step(0.0001);
        }
        let y = sim.particles[0].position[1];
        // Should be resting near radius height
        assert!(y > 0.0, "particle above floor: y={y}");
        assert!(y < 0.1, "particle not too high: y={y}");
    }

    #[test]
    fn simulation_kinetic_energy_initially_zero() {
        let sim = DemSimulation::new(DemContact::glass_beads(), 9.81);
        assert_eq!(sim.total_kinetic_energy(), 0.0);
    }

    #[test]
    fn simulation_kinetic_energy_grows_under_gravity() {
        let mut sim = DemSimulation::new(DemContact::glass_beads(), 9.81);
        sim.add_particle(make_particle([0.0, 1.0, 0.0], 0.01));
        let ke0 = sim.total_kinetic_energy();
        for _ in 0..10 {
            sim.step(0.001);
        }
        assert!(sim.total_kinetic_energy() > ke0);
    }

    // ── GranularFlow ──────────────────────────────────────────────────────

    #[test]
    fn granular_flow_counts_crossings() {
        let mut flow = GranularFlow::new([0.0; 3], [1.0, 0.0, 0.0]);
        let mut particles = vec![make_particle([-0.01, 0.0, 0.0], 0.005)];
        flow.update(&particles);
        // Move past the plane
        particles[0].position[0] = 0.01;
        let crossings = flow.update(&particles);
        assert_eq!(crossings, 1);
        assert_eq!(flow.crossing_count, 1);
    }

    #[test]
    fn granular_flow_no_crossing_same_side() {
        let mut flow = GranularFlow::new([0.0; 3], [1.0, 0.0, 0.0]);
        let particles = vec![make_particle([0.1, 0.0, 0.0], 0.005)];
        flow.update(&particles);
        let crossings = flow.update(&particles);
        assert_eq!(crossings, 0);
    }

    #[test]
    fn granular_flow_mean_normal_velocity() {
        let flow = GranularFlow::new([0.0; 3], [1.0, 0.0, 0.0]);
        let mut p = make_particle([0.0; 3], 0.005);
        p.velocity = [2.0, 0.0, 0.0];
        let mean = flow.mean_normal_velocity(&[p]);
        assert!((mean - 2.0).abs() < 1e-12);
    }

    // ── AngleOfRepose ─────────────────────────────────────────────────────

    #[test]
    fn angle_of_repose_45_deg() {
        let a = AngleOfRepose::new(1.0, 1.0);
        assert!((a.angle_deg() - 45.0).abs() < 1e-10);
    }

    #[test]
    fn angle_of_repose_shallow_pile() {
        let a = AngleOfRepose::new(0.5, 2.0);
        assert!(a.angle_deg() < 30.0);
    }

    #[test]
    fn angle_of_repose_from_particles() {
        let mut particles = Vec::new();
        // Pile: base from x=-1 to x=1, peak at y=0.5
        for xi in &[-1.0_f64, 0.0, 1.0] {
            let mut p = make_particle([*xi, 0.5 * (1.0 - xi.abs()), 0.0], 0.05);
            p.position[1] += 0.05; // radius offset
            particles.push(p);
        }
        let a = AngleOfRepose::from_particles(&particles, 0.0);
        assert!(a.angle_deg() > 0.0);
        assert!(a.angle_deg() < 90.0);
    }

    // ── PackingDensity ────────────────────────────────────────────────────

    #[test]
    fn packing_density_single_sphere() {
        let domain_vol = 1.0;
        let pd = PackingDensity::new(domain_vol);
        let particles = vec![make_particle([0.0; 3], 0.5)]; // radius 0.5 in unit cube
        let vf = pd.volume_fraction(&particles);
        let expected = (4.0 / 3.0) * PI * 0.5_f64.powi(3);
        assert!((vf - expected).abs() < 1e-12);
    }

    #[test]
    fn packing_coordination_no_contact() {
        let pd = PackingDensity::new(1.0);
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], 0.01),
            make_particle([1.0, 0.0, 0.0], 0.01),
        ];
        let cn = pd.coordination_number(&particles, 0.001);
        assert_eq!(cn, 0.0);
    }

    #[test]
    fn packing_coordination_in_contact() {
        let pd = PackingDensity::new(1.0);
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], 0.01),
            make_particle([0.02, 0.0, 0.0], 0.01), // touching
        ];
        let cn = pd.coordination_number(&particles, 0.001);
        assert!(cn > 0.0);
    }

    // ── AvalancheDynamics ─────────────────────────────────────────────────

    #[test]
    fn avalanche_detects_onset() {
        let mut av = AvalancheDynamics::new(0.01, 0.001);
        let mut particles = vec![make_particle([0.0; 3], 0.01)];
        particles[0].velocity = [10.0, 0.0, 0.0]; // high KE
        let started = av.update(&particles);
        assert!(started);
        assert!(av.in_avalanche);
    }

    #[test]
    fn avalanche_no_onset_below_threshold() {
        let mut av = AvalancheDynamics::new(100.0, 1.0);
        let particles = vec![make_particle([0.0; 3], 0.01)];
        let started = av.update(&particles);
        assert!(!started);
    }

    #[test]
    fn avalanche_detects_cessation() {
        let mut av = AvalancheDynamics::new(0.01, 0.001);
        let mut particles = vec![make_particle([0.0; 3], 0.01)];
        particles[0].velocity = [10.0, 0.0, 0.0];
        av.update(&particles);
        // Stop
        particles[0].velocity = [0.0; 3];
        av.update(&particles);
        assert!(!av.in_avalanche);
        assert_eq!(av.event_count, 1);
    }

    #[test]
    fn avalanche_mean_event_size_empty() {
        let av = AvalancheDynamics::new(0.01, 0.001);
        assert_eq!(av.mean_event_size(), 0.0);
    }

    // ── GranularTemperature ───────────────────────────────────────────────

    #[test]
    fn granular_temperature_zero_for_identical_velocities() {
        let mut p1 = make_particle([0.0; 3], 0.01);
        let mut p2 = make_particle([1.0, 0.0, 0.0], 0.01);
        p1.velocity = [1.0, 0.0, 0.0];
        p2.velocity = [1.0, 0.0, 0.0];
        let t = GranularTemperature::translational(&[p1, p2]);
        assert!(t.abs() < 1e-14);
    }

    #[test]
    fn granular_temperature_nonzero_for_varying_velocities() {
        let mut p1 = make_particle([0.0; 3], 0.01);
        let mut p2 = make_particle([1.0, 0.0, 0.0], 0.01);
        p1.velocity = [2.0, 0.0, 0.0];
        p2.velocity = [-2.0, 0.0, 0.0];
        let t = GranularTemperature::translational(&[p1, p2]);
        assert!(t > 0.0);
    }

    #[test]
    fn granular_temperature_empty_particles() {
        assert_eq!(GranularTemperature::translational(&[]), 0.0);
        assert_eq!(GranularTemperature::rotational(&[]), 0.0);
        assert_eq!(GranularTemperature::total(&[]), 0.0);
    }

    #[test]
    fn granular_temperature_rotational_nonzero() {
        let mut p1 = make_particle([0.0; 3], 0.01);
        let mut p2 = make_particle([1.0, 0.0, 0.0], 0.01);
        p1.angular_velocity = [0.0, 5.0, 0.0];
        p2.angular_velocity = [0.0, -5.0, 0.0];
        let t = GranularTemperature::rotational(&[p1, p2]);
        assert!(t > 0.0);
    }

    #[test]
    fn granular_temperature_total_is_sum() {
        let mut p1 = make_particle([0.0; 3], 0.01);
        let mut p2 = make_particle([1.0, 0.0, 0.0], 0.01);
        p1.velocity = [2.0, 0.0, 0.0];
        p2.velocity = [-2.0, 0.0, 0.0];
        p1.angular_velocity = [0.0, 3.0, 0.0];
        p2.angular_velocity = [0.0, -3.0, 0.0];
        let particles = [p1, p2];
        let total = GranularTemperature::total(&particles);
        let trans = GranularTemperature::translational(&particles);
        let rot = GranularTemperature::rotational(&particles);
        assert!((total - trans - rot).abs() < 1e-14);
    }

    // ── vec3 helpers ──────────────────────────────────────────────────────

    #[test]
    fn vec3_cross_product() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn vec3_normalize_unit_vector() {
        let v = [3.0, 4.0, 0.0];
        let n = vec3_normalize(v);
        assert!((vec3_norm(n) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn vec3_normalize_zero_returns_zero() {
        let z = vec3_normalize([0.0; 3]);
        assert_eq!(z, [0.0; 3]);
    }
}
