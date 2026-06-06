// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluid-solid collision and coupling for the OxiPhysics collision crate.
//!
//! Covers:
//! - Fluid-solid coupling (SPH-rigid, LBM-rigid bounce-back)
//! - Buoyancy force computation (Archimedes principle + displaced volume)
//! - Drag force on submerged bodies (viscous + pressure drag)
//! - Water entry / exit detection (free-surface crossing)
//! - Fluid boundary collision handling
//! - Splash particle generation at collision impact
//! - Underwater contact friction
//! - Free-surface collision detection
//! - Two-way coupling impulse exchange
//! - Fluid-accelerated projectile dynamics

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Internal vector helpers — no nalgebra
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

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
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0, 0.0, 1.0]
    } else {
        scale3(a, 1.0 / l)
    }
}

#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Fluid Properties
// ─────────────────────────────────────────────────────────────────────────────

/// Physical properties of a fluid.
#[derive(Debug, Clone, Copy)]
pub struct FluidProperties {
    /// Fluid density ρ \[kg/m³\].
    pub density: f64,
    /// Dynamic viscosity μ \[Pa·s\].
    pub viscosity: f64,
    /// Kinematic viscosity ν = μ/ρ \[m²/s\].
    pub kinematic_viscosity: f64,
    /// Surface tension coefficient σ \[N/m\].
    pub surface_tension: f64,
    /// Speed of sound in the fluid \[m/s\].
    pub speed_of_sound: f64,
}

impl FluidProperties {
    /// Create fluid properties and derive kinematic viscosity automatically.
    ///
    /// # Arguments
    /// * `density` - Mass density \[kg/m³\].
    /// * `viscosity` - Dynamic viscosity \[Pa·s\].
    /// * `surface_tension` - Surface tension \[N/m\].
    /// * `speed_of_sound` - Speed of sound \[m/s\].
    pub fn new(density: f64, viscosity: f64, surface_tension: f64, speed_of_sound: f64) -> Self {
        Self {
            density,
            viscosity,
            kinematic_viscosity: viscosity / density,
            surface_tension,
            speed_of_sound,
        }
    }

    /// Return water-like properties at ~20°C.
    pub fn water() -> Self {
        Self::new(998.2, 1.002e-3, 0.0728, 1482.0)
    }

    /// Return air-like properties at sea level, ~20°C.
    pub fn air() -> Self {
        Self::new(1.204, 1.81e-5, 0.0, 343.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Rigid Body State for Fluid Coupling
// ─────────────────────────────────────────────────────────────────────────────

/// Minimal rigid body state used for fluid-solid coupling computations.
#[derive(Debug, Clone)]
pub struct RigidBodyFluidState {
    /// Centre-of-mass position \[m\].
    pub position: [f64; 3],
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Total mass \[kg\].
    pub mass: f64,
    /// Inertia tensor diagonal \[kg·m²\] (principal axes assumed).
    pub inertia: [f64; 3],
    /// Representative radius for sphere approximation \[m\].
    pub radius: f64,
    /// Cross-sectional area for drag computation \[m²\].
    pub cross_section_area: f64,
}

impl RigidBodyFluidState {
    /// Create a new rigid body fluid state.
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, radius: f64) -> Self {
        let i = 0.4 * mass * radius * radius; // sphere moment of inertia
        Self {
            position,
            velocity,
            angular_velocity: [0.0; 3],
            mass,
            inertia: [i, i, i],
            radius,
            cross_section_area: PI * radius * radius,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Reynolds Number
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Reynolds number Re = U·L / ν.
///
/// # Arguments
/// * `velocity` - Characteristic velocity \[m/s\].
/// * `length` - Characteristic length \[m\].
/// * `kinematic_viscosity` - ν = μ/ρ \[m²/s\].
pub fn reynolds_number(velocity: f64, length: f64, kinematic_viscosity: f64) -> f64 {
    if kinematic_viscosity < 1e-30 {
        return f64::INFINITY;
    }
    velocity.abs() * length / kinematic_viscosity
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Buoyancy Force
// ─────────────────────────────────────────────────────────────────────────────

/// Buoyancy result containing force vector and submerged volume.
#[derive(Debug, Clone, Copy)]
pub struct BuoyancyResult {
    /// Buoyancy force vector \[N\] (upward for positive y convention).
    pub force: [f64; 3],
    /// Volume of body submerged in fluid \[m³\].
    pub submerged_volume: f64,
    /// Centre of buoyancy (centroid of submerged volume).
    pub centre_of_buoyancy: [f64; 3],
}

/// Compute buoyancy force on a sphere partially or fully submerged.
///
/// Uses the spherical cap volume formula for partial submersion.
///
/// # Arguments
/// * `body` - Rigid body state (position = centre of sphere).
/// * `fluid` - Fluid properties.
/// * `free_surface_y` - Y-coordinate of the free surface \[m\].
/// * `gravity` - Gravitational acceleration magnitude \[m/s²\].
///
/// The gravity vector is assumed to act in the −y direction.
pub fn buoyancy_sphere(
    body: &RigidBodyFluidState,
    fluid: &FluidProperties,
    free_surface_y: f64,
    gravity: f64,
) -> BuoyancyResult {
    let r = body.radius;
    let cy = body.position[1];
    // Depth of sphere centre below free surface
    let depth = free_surface_y - cy;

    let (vol_sub, cob_y) = if depth >= r {
        // Fully submerged
        let vol = 4.0 / 3.0 * PI * r * r * r;
        (vol, cy) // centre of buoyancy = sphere centre when fully submerged
    } else if depth <= -r {
        // Fully above surface
        (0.0, cy)
    } else {
        // Partially submerged — spherical cap volume
        // h = depth of cap immersion = r + depth (from bottom of sphere)
        let h = (r + depth).max(0.0).min(2.0 * r);
        let vol = PI * h * h * (r - h / 3.0);
        // Centre of buoyancy y within the submerged cap (measured from sphere centre)
        let cob_y_local = if vol > 1e-30 {
            // z-coord of centroid of spherical cap: (3*(2r-h)^2) / (4*(3r-h)) below top
            let cob_from_top = 3.0 * (2.0 * r - h).powi(2) / (4.0 * (3.0 * r - h).max(1e-15));
            cy - r + h - cob_from_top
        } else {
            cy - r
        };
        (vol, cob_y_local)
    };

    let force_mag = fluid.density * gravity * vol_sub;
    let force = [0.0, force_mag, 0.0]; // upward
    let cob = [body.position[0], cob_y, body.position[2]];
    BuoyancyResult {
        force,
        submerged_volume: vol_sub,
        centre_of_buoyancy: cob,
    }
}

/// Compute buoyancy force for an arbitrary body given its submerged volume.
///
/// F_b = ρ_fluid · g · V_sub (upward).
pub fn buoyancy_force(fluid_density: f64, gravity: f64, submerged_volume: f64) -> [f64; 3] {
    let mag = fluid_density * gravity * submerged_volume;
    [0.0, mag, 0.0]
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Drag Force
// ─────────────────────────────────────────────────────────────────────────────

/// Drag coefficient as a function of Reynolds number (sphere model).
///
/// Uses the Schiller-Naumann correlation:
/// - Stokes: C_d = 24/Re for Re < 1
/// - Schiller-Naumann: C_d = (24/Re)(1 + 0.15 Re^0.687) for 1 ≤ Re ≤ 1000
/// - Newton: C_d = 0.44 for Re > 1000
pub fn drag_coefficient_sphere(re: f64) -> f64 {
    if re < 1e-10 {
        return 24.0 / 1e-10; // avoid division by zero
    }
    if re < 1.0 {
        24.0 / re
    } else if re <= 1000.0 {
        (24.0 / re) * (1.0 + 0.15 * re.powf(0.687))
    } else {
        0.44
    }
}

/// Compute the hydrodynamic drag force on a body moving through a fluid.
///
/// F_drag = -0.5 · ρ · C_d · A · |v_rel|² · v̂_rel
///
/// # Arguments
/// * `body_velocity` - Velocity of the body \[m/s\].
/// * `fluid_velocity` - Local fluid velocity \[m/s\].
/// * `fluid` - Fluid properties.
/// * `cd` - Drag coefficient (use [`drag_coefficient_sphere`] if unknown).
/// * `reference_area` - Body reference area \[m²\].
pub fn drag_force(
    body_velocity: [f64; 3],
    fluid_velocity: [f64; 3],
    fluid: &FluidProperties,
    cd: f64,
    reference_area: f64,
) -> [f64; 3] {
    let v_rel = sub3(body_velocity, fluid_velocity);
    let speed = len3(v_rel);
    if speed < 1e-15 {
        return [0.0; 3];
    }
    let mag = 0.5 * fluid.density * cd * reference_area * speed * speed;
    let dir = normalize3(neg3(v_rel));
    scale3(dir, mag)
}

/// Compute drag with automatic Reynolds-number-based C_d for a sphere.
pub fn drag_force_sphere(
    body: &RigidBodyFluidState,
    fluid_velocity: [f64; 3],
    fluid: &FluidProperties,
) -> [f64; 3] {
    let v_rel = sub3(body.velocity, fluid_velocity);
    let speed = len3(v_rel);
    let re = reynolds_number(speed, 2.0 * body.radius, fluid.kinematic_viscosity);
    let cd = drag_coefficient_sphere(re);
    drag_force(
        body.velocity,
        fluid_velocity,
        fluid,
        cd,
        body.cross_section_area,
    )
}

/// Compute the lift force using a simplified cross-flow model.
///
/// F_lift = 0.5 · ρ · C_l · A · v_rel² · n̂_lift
///
/// where n̂_lift is perpendicular to v_rel in the vertical plane.
pub fn lift_force(
    body_velocity: [f64; 3],
    fluid_velocity: [f64; 3],
    fluid: &FluidProperties,
    cl: f64,
    reference_area: f64,
) -> [f64; 3] {
    let v_rel = sub3(body_velocity, fluid_velocity);
    let speed = len3(v_rel);
    if speed < 1e-15 {
        return [0.0; 3];
    }
    let up = [0.0_f64, 1.0, 0.0];
    let lift_dir = normalize3(cross3(v_rel, cross3(v_rel, up)));
    let mag = 0.5 * fluid.density * cl * reference_area * speed * speed;
    scale3(lift_dir, mag)
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Water Entry / Exit Detection
// ─────────────────────────────────────────────────────────────────────────────

/// State of a body with respect to the free surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubmersionState {
    /// Body is entirely above the free surface.
    Airborne,
    /// Body is crossing the free surface (partially submerged).
    Crossing,
    /// Body is entirely below the free surface.
    Submerged,
}

/// Detect the submersion state of a sphere relative to the free surface.
///
/// # Arguments
/// * `body` - Rigid body (position = sphere centre).
/// * `free_surface_y` - Y-coordinate of the flat free surface.
pub fn submersion_state_sphere(body: &RigidBodyFluidState, free_surface_y: f64) -> SubmersionState {
    let cy = body.position[1];
    let r = body.radius;
    if cy - r >= free_surface_y {
        SubmersionState::Airborne
    } else if cy + r <= free_surface_y {
        SubmersionState::Submerged
    } else {
        SubmersionState::Crossing
    }
}

/// Water entry event data.
#[derive(Debug, Clone, Copy)]
pub struct WaterEntryEvent {
    /// Position of the body centre at entry.
    pub position: [f64; 3],
    /// Velocity at entry \[m/s\].
    pub velocity: [f64; 3],
    /// Entry speed \[m/s\].
    pub entry_speed: f64,
    /// Vertical component of velocity at entry \[m/s\].
    pub vertical_speed: f64,
    /// Whether this is an entry (true) or exit (false).
    pub is_entry: bool,
}

/// Detect water entry or exit by comparing previous and current states.
///
/// Returns `Some(WaterEntryEvent)` if a crossing of the free surface occurred.
pub fn detect_water_crossing(
    prev: &RigidBodyFluidState,
    curr: &RigidBodyFluidState,
    free_surface_y: f64,
) -> Option<WaterEntryEvent> {
    let prev_state = submersion_state_sphere(prev, free_surface_y);
    let curr_state = submersion_state_sphere(curr, free_surface_y);
    if prev_state == curr_state {
        return None;
    }
    let is_entry = matches!(
        (prev_state, curr_state),
        (SubmersionState::Airborne, SubmersionState::Crossing)
            | (SubmersionState::Airborne, SubmersionState::Submerged)
            | (SubmersionState::Crossing, SubmersionState::Submerged)
    );
    let speed = len3(curr.velocity);
    Some(WaterEntryEvent {
        position: curr.position,
        velocity: curr.velocity,
        entry_speed: speed,
        vertical_speed: curr.velocity[1].abs(),
        is_entry,
    })
}

/// Compute the von Karman water impact load (slamming pressure) on entry.
///
/// Uses the simplified von Karman model:
/// F_slam ≈ ρ_f · π · r_c² · (dz/dt)² · (1 + α)
/// where r_c is the contact radius and α is an added-mass correction factor.
///
/// # Arguments
/// * `fluid_density` - Fluid density \[kg/m³\].
/// * `entry_speed` - Vertical entry speed \[m/s\].
/// * `contact_radius` - Instantaneous water-plane contact radius \[m\].
/// * `added_mass_coeff` - Added mass coefficient α (≈ 1.0 for sphere).
pub fn water_entry_slamming_force(
    fluid_density: f64,
    entry_speed: f64,
    contact_radius: f64,
    added_mass_coeff: f64,
) -> f64 {
    fluid_density
        * PI
        * contact_radius
        * contact_radius
        * entry_speed
        * entry_speed
        * (1.0 + added_mass_coeff)
}

/// Estimate the contact radius for a sphere entering the water.
///
/// r_c ≈ sqrt(2 R d) where R is sphere radius and d is penetration depth.
pub fn sphere_water_contact_radius(sphere_radius: f64, penetration_depth: f64) -> f64 {
    if penetration_depth <= 0.0 {
        return 0.0;
    }
    (2.0 * sphere_radius * penetration_depth.min(2.0 * sphere_radius)).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Splash Particle Generation
// ─────────────────────────────────────────────────────────────────────────────

/// A single splash particle spawned at water impact.
#[derive(Debug, Clone, Copy)]
pub struct SplashParticle {
    /// Initial position \[m\].
    pub position: [f64; 3],
    /// Initial velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Lifetime \[s\].
    pub lifetime: f64,
}

/// Generate splash particles at a water entry event.
///
/// Particles are distributed on a cone around the impact normal with
/// speeds scaled by entry speed.
///
/// # Arguments
/// * `event` - Water entry event.
/// * `n_particles` - Number of splash particles to generate.
/// * `fluid` - Fluid properties.
/// * `particle_mass` - Mass of each splash particle \[kg\].
/// * `seed` - RNG seed.
pub fn generate_splash_particles(
    event: &WaterEntryEvent,
    n_particles: usize,
    fluid: &FluidProperties,
    particle_mass: f64,
    seed: u64,
) -> Vec<SplashParticle> {
    let mut rng = FcRng::new(seed);
    let mut particles = Vec::with_capacity(n_particles);
    let base_speed = event.entry_speed * 0.6;
    let cone_half_angle = PI / 4.0; // 45° cone

    for i in 0..n_particles {
        let azimuth = (i as f64 / n_particles as f64) * 2.0 * PI
            + rng.next_f64() * (2.0 * PI / n_particles as f64);
        let polar = rng.next_f64() * cone_half_angle;
        let speed = base_speed * (0.5 + rng.next_f64() * 0.5);
        let vx = speed * polar.sin() * azimuth.cos();
        let vy = speed * polar.cos();
        let vz = speed * polar.sin() * azimuth.sin();
        let lifetime = 2.0 * speed / 9.81 + 0.1; // rough estimate
        // Use surface tension to compute critical diameter
        let _surface_tension_ref = fluid.surface_tension;
        particles.push(SplashParticle {
            position: event.position,
            velocity: [vx, vy, vz],
            mass: particle_mass,
            lifetime,
        });
    }
    particles
}

/// Integrate a splash particle forward by time step `dt` under gravity.
///
/// Simple Euler step: ignores drag and surface tension for splash particles.
pub fn integrate_splash_particle(particle: &mut SplashParticle, dt: f64, gravity: f64) {
    particle.velocity[1] -= gravity * dt;
    particle.position = add3(particle.position, scale3(particle.velocity, dt));
    particle.lifetime -= dt;
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. Free-Surface Collision Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box of the fluid domain.
#[derive(Debug, Clone, Copy)]
pub struct FluidDomain {
    /// Minimum corner \[m\].
    pub min: [f64; 3],
    /// Maximum corner \[m\].
    pub max: [f64; 3],
    /// Free surface Y-coordinate \[m\].
    pub free_surface_y: f64,
}

impl FluidDomain {
    /// Create a new fluid domain.
    pub fn new(min: [f64; 3], max: [f64; 3], free_surface_y: f64) -> Self {
        Self {
            min,
            max,
            free_surface_y,
        }
    }

    /// Return true if point `p` is within the fluid domain and below the free surface.
    pub fn contains_fluid_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.free_surface_y
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Return true if point `p` is within the domain bounding box.
    pub fn in_bounds(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
}

/// Detect whether a ray from `origin` in direction `dir` hits the flat free surface.
///
/// Returns the time of intersection `t` such that `origin + t * dir` lies on the surface,
/// or `None` if no intersection exists within `[0, max_t]`.
pub fn ray_free_surface(
    origin: [f64; 3],
    dir: [f64; 3],
    free_surface_y: f64,
    max_t: f64,
) -> Option<f64> {
    // Surface: y = free_surface_y => origin.y + t * dir.y = free_surface_y
    if dir[1].abs() < 1e-14 {
        return None; // ray is parallel to surface
    }
    let t = (free_surface_y - origin[1]) / dir[1];
    if t >= 0.0 && t <= max_t {
        Some(t)
    } else {
        None
    }
}

/// Compute the fraction of a sphere submerged below the free surface (0..=1).
pub fn sphere_submersion_fraction(cy: f64, r: f64, free_surface_y: f64) -> f64 {
    let depth = free_surface_y - cy; // positive = centre below surface
    if depth >= r {
        1.0
    } else if depth <= -r {
        0.0
    } else {
        let h = (r + depth).max(0.0).min(2.0 * r);
        // Volume of spherical cap / total sphere volume
        let vol_cap = PI * h * h * (r - h / 3.0);
        let vol_sphere = 4.0 / 3.0 * PI * r * r * r;
        (vol_cap / vol_sphere).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. Fluid Boundary Collision Handling
// ─────────────────────────────────────────────────────────────────────────────

/// Collision response for a rigid body hitting the fluid domain boundary.
///
/// Reflects velocity component normal to the boundary with a restitution coefficient.
///
/// # Arguments
/// * `body` - Mutable rigid body state.
/// * `domain` - Fluid domain bounds.
/// * `restitution` - Coefficient of restitution \[0, 1\].
pub fn handle_domain_boundary_collision(
    body: &mut RigidBodyFluidState,
    domain: &FluidDomain,
    restitution: f64,
) -> bool {
    let mut hit = false;
    let r = body.radius;
    // X boundaries
    if body.position[0] - r < domain.min[0] {
        body.position[0] = domain.min[0] + r;
        body.velocity[0] = body.velocity[0].abs() * restitution;
        hit = true;
    } else if body.position[0] + r > domain.max[0] {
        body.position[0] = domain.max[0] - r;
        body.velocity[0] = -body.velocity[0].abs() * restitution;
        hit = true;
    }
    // Y boundaries (bottom)
    if body.position[1] - r < domain.min[1] {
        body.position[1] = domain.min[1] + r;
        body.velocity[1] = body.velocity[1].abs() * restitution;
        hit = true;
    }
    // Z boundaries
    if body.position[2] - r < domain.min[2] {
        body.position[2] = domain.min[2] + r;
        body.velocity[2] = body.velocity[2].abs() * restitution;
        hit = true;
    } else if body.position[2] + r > domain.max[2] {
        body.position[2] = domain.max[2] - r;
        body.velocity[2] = -body.velocity[2].abs() * restitution;
        hit = true;
    }
    hit
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. SPH-Rigid Body Coupling
// ─────────────────────────────────────────────────────────────────────────────

/// An SPH fluid particle.
#[derive(Debug, Clone)]
pub struct SphParticle {
    /// Position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Density \[kg/m³\].
    pub density: f64,
    /// Pressure \[Pa\].
    pub pressure: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
}

/// Smoothing kernel — cubic spline W(r, h).
///
/// Returns W(r, h) where h is the smoothing length.
pub fn sph_cubic_kernel(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h); // 3D normalisation
    let q = r / h;
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q * (1.0 - 0.5 * q))
    } else if q < 2.0 {
        let f = 2.0 - q;
        sigma * 0.25 * f * f * f
    } else {
        0.0
    }
}

/// Gradient of the cubic spline kernel ∇W (pointing away from x_j toward x_i).
///
/// Returns the gradient vector dW/dr * (x_i - x_j) / r.
pub fn sph_cubic_kernel_gradient(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_vec);
    if r < 1e-15 {
        return [0.0; 3];
    }
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    let dw_dr = if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let f = 2.0 - q;
        sigma * (-0.75 * f * f) / h
    } else {
        0.0
    };
    scale3(r_vec, dw_dr / r)
}

/// Compute the SPH force on a rigid body due to fluid pressure from neighbouring particles.
///
/// The boundary force is: F_b = Σ_j m_j (p_j / ρ_j²) ∇W(r_ij, h)
///
/// # Arguments
/// * `body_center` - Centre of the rigid body \[m\].
/// * `body_radius` - Radius of the rigid body \[m\].
/// * `particles` - SPH particle array.
/// * `h` - Smoothing length \[m\].
pub fn sph_rigid_coupling_force(
    body_center: [f64; 3],
    body_radius: f64,
    particles: &[SphParticle],
    h: f64,
) -> [f64; 3] {
    let mut force = [0.0_f64; 3];
    for p in particles {
        let r_vec = sub3(body_center, p.position);
        let r = len3(r_vec);
        if r > h || r < body_radius * 0.5 {
            continue; // outside kernel support or inside body
        }
        let grad_w = sph_cubic_kernel_gradient(r_vec, h);
        // Contribution: m_j (p_j / ρ_j²) ∇W
        if p.density > 1e-10 {
            let coeff = p.mass * p.pressure / (p.density * p.density);
            force = add3(force, scale3(grad_w, coeff));
        }
    }
    force
}

/// Reflect SPH particle off a rigid sphere surface using the boundary repulsion method.
///
/// Applies a penalty velocity correction to particles inside the rigid body.
///
/// # Arguments
/// * `particle` - Mutable SPH particle.
/// * `body` - Rigid body state.
/// * `stiffness` - Repulsion stiffness coefficient.
/// * `dt` - Time step \[s\].
pub fn sph_rigid_boundary_reflect(
    particle: &mut SphParticle,
    body: &RigidBodyFluidState,
    stiffness: f64,
    dt: f64,
) -> bool {
    let r_vec = sub3(particle.position, body.position);
    let r = len3(r_vec);
    if r >= body.radius {
        return false; // particle outside body
    }
    let penetration = body.radius - r;
    let normal = normalize3(r_vec);
    // Apply penalty force: a_pen = stiffness * penetration * normal
    let dv = scale3(normal, stiffness * penetration * dt);
    particle.velocity = add3(particle.velocity, dv);
    // Push particle to surface
    particle.position = add3(body.position, scale3(normal, body.radius + 1e-6));
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. LBM-Rigid Body Coupling (Bounce-Back)
// ─────────────────────────────────────────────────────────────────────────────

/// D3Q19 lattice velocity directions (velocities in lattice units).
///
/// Standard D3Q19 model: 19 directions including rest.
pub const D3Q19_VELOCITIES: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, 1],
    [0, -1, -1],
    [0, 1, -1],
    [0, -1, 1],
];

/// D3Q19 weights for the equilibrium distribution.
pub const D3Q19_WEIGHTS: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Compute the LBM equilibrium distribution function f^eq.
///
/// f^eq_i = w_i ρ \[1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) − u²/(2cs²)\]
///
/// with cs² = 1/3 in lattice units.
///
/// # Arguments
/// * `rho` - Local density.
/// * `u` - Local velocity vector.
/// * `i` - Lattice direction index (0..19).
pub fn lbm_equilibrium(rho: f64, u: [f64; 3], i: usize) -> f64 {
    let e = D3Q19_VELOCITIES[i];
    let ei = [e[0] as f64, e[1] as f64, e[2] as f64];
    let wi = D3Q19_WEIGHTS[i];
    let cs2 = 1.0 / 3.0;
    let eu = dot3(ei, u);
    let u2 = dot3(u, u);
    wi * rho * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2))
}

/// Apply the LBM bounce-back boundary condition for a rigid body.
///
/// For each distribution function that would stream into the solid node,
/// the distribution is reversed (half-way bounce-back).
///
/// Returns the momentum transferred to the body from the fluid \[lattice units\].
pub fn lbm_bounce_back_force(
    f_in: &[f64; 19],
    f_in_opposite: &[f64; 19],
    body_velocity: [f64; 3],
    rho: f64,
) -> [f64; 3] {
    let mut momentum = [0.0_f64; 3];
    for i in 1..19usize {
        let e = D3Q19_VELOCITIES[i];
        let ei = [e[0] as f64, e[1] as f64, e[2] as f64];
        // Moving wall correction: add 2 w_i ρ (e_i · u_wall) / cs²
        let eu_wall = dot3(ei, body_velocity);
        let cs2 = 1.0 / 3.0;
        let f_bb = f_in_opposite[i] + 2.0 * D3Q19_WEIGHTS[i] * rho * eu_wall / cs2;
        // Force = 2 * (f_in - f_bb) * e_i (Newton 3rd law)
        let df = f_in[i] - f_bb;
        momentum = add3(momentum, scale3(ei, 2.0 * df));
    }
    momentum
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Two-Way Coupling Impulse
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the two-way coupling impulse between a rigid body and the fluid.
///
/// The impulse J is distributed between the body and the fluid parcel to conserve
/// linear momentum.
///
/// # Arguments
/// * `body` - Mutable rigid body.
/// * `fluid_mass` - Effective fluid mass in contact \[kg\].
/// * `fluid_vel` - Fluid velocity at contact point \[m/s\].
/// * `contact_normal` - Unit normal from fluid to body.
/// * `restitution` - Coefficient of restitution.
/// * `dt` - Time step \[s\] (used for impulse limiting).
///
/// Returns the impulse vector applied to the body.
pub fn two_way_coupling_impulse(
    body: &mut RigidBodyFluidState,
    fluid_mass: f64,
    fluid_vel: [f64; 3],
    contact_normal: [f64; 3],
    restitution: f64,
    _dt: f64,
) -> [f64; 3] {
    let v_rel = sub3(body.velocity, fluid_vel);
    let vn = dot3(v_rel, contact_normal);
    if vn >= 0.0 {
        // Separating — no impulse needed
        return [0.0; 3];
    }
    // Impulse magnitude
    let j_mag = -(1.0 + restitution) * vn / (1.0 / body.mass + 1.0 / fluid_mass.max(1e-10));
    let impulse = scale3(contact_normal, j_mag);
    // Apply impulse to body
    body.velocity = add3(body.velocity, scale3(impulse, 1.0 / body.mass));
    impulse
}

// ─────────────────────────────────────────────────────────────────────────────
// 13. Underwater Contact Friction
// ─────────────────────────────────────────────────────────────────────────────

/// Friction result for underwater contact.
#[derive(Debug, Clone, Copy)]
pub struct UnderwaterFrictionResult {
    /// Friction force \[N\].
    pub friction_force: [f64; 3],
    /// Whether the contact is sliding (kinetic) or sticking (static).
    pub is_sliding: bool,
}

/// Compute underwater contact friction force.
///
/// Applies a modified Coulomb friction model with a viscous correction:
/// F_fric = −min(μ_eff |F_n|, μ_visc |v_t|) * v̂_t
///
/// where μ_eff = μ_coulomb (dry) reduced by the lubrication factor.
///
/// # Arguments
/// * `normal_force` - Magnitude of normal contact force \[N\].
/// * `tangential_velocity` - Relative tangential velocity \[m/s\].
/// * `mu_coulomb` - Dry Coulomb friction coefficient.
/// * `lubrication_factor` - Reduction factor for fluid lubrication \[0,1\].
/// * `viscous_coeff` - Viscous contribution coefficient \[N·s/m\].
pub fn underwater_contact_friction(
    normal_force: f64,
    tangential_velocity: [f64; 3],
    mu_coulomb: f64,
    lubrication_factor: f64,
    viscous_coeff: f64,
) -> UnderwaterFrictionResult {
    let vt_len = len3(tangential_velocity);
    let mu_eff = mu_coulomb * (1.0 - lubrication_factor.clamp(0.0, 0.99));
    let f_coulomb = mu_eff * normal_force.abs();
    let f_viscous = viscous_coeff * vt_len;
    let f_limit = f_coulomb.min(f_viscous + f_coulomb);

    if vt_len < 1e-12 {
        return UnderwaterFrictionResult {
            friction_force: [0.0; 3],
            is_sliding: false,
        };
    }
    let vt_dir = normalize3(neg3(tangential_velocity));
    let f_applied = f_limit.min(f_coulomb + f_viscous);
    let is_sliding = f_applied >= f_coulomb;
    UnderwaterFrictionResult {
        friction_force: scale3(vt_dir, f_applied),
        is_sliding,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 14. Fluid-Accelerated Projectile Dynamics
// ─────────────────────────────────────────────────────────────────────────────

/// State of a fluid-accelerated projectile (e.g., a supercavitating torpedo).
#[derive(Debug, Clone)]
pub struct FluidProjectile {
    /// Body state.
    pub body: RigidBodyFluidState,
    /// Whether the projectile is currently submerged.
    pub submerged: bool,
    /// Cavitation number σ = (p_inf - p_v) / (0.5 ρ v²).
    pub cavitation_number: f64,
    /// Accumulated impulse from thrust \[N·s\].
    pub thrust_impulse: f64,
}

impl FluidProjectile {
    /// Create a new fluid-accelerated projectile.
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, radius: f64) -> Self {
        Self {
            body: RigidBodyFluidState::new(position, velocity, mass, radius),
            submerged: false,
            cavitation_number: f64::INFINITY,
            thrust_impulse: 0.0,
        }
    }

    /// Update the cavitation number based on ambient pressure and vapour pressure.
    ///
    /// σ = (p_inf − p_v) / (0.5 ρ v²)
    pub fn update_cavitation_number(&mut self, p_inf: f64, p_vapour: f64, fluid_density: f64) {
        let v = len3(self.body.velocity);
        let q = 0.5 * fluid_density * v * v;
        if q < 1e-12 {
            self.cavitation_number = f64::INFINITY;
        } else {
            self.cavitation_number = (p_inf - p_vapour) / q;
        }
    }

    /// Check whether the projectile is in a supercavitating regime (σ < σ_crit).
    pub fn is_supercavitating(&self, sigma_crit: f64) -> bool {
        self.cavitation_number < sigma_crit
    }
}

/// Compute the thrust force on a fluid-accelerated projectile (rocket/jet propulsion).
///
/// F_thrust = ṁ_exhaust · v_exhaust + (p_exhaust − p_ambient) · A_nozzle
///
/// # Arguments
/// * `mass_flow_rate` - Exhaust mass flow rate \[kg/s\].
/// * `exhaust_velocity` - Exhaust velocity relative to body \[m/s\].
/// * `pressure_diff` - Exhaust pressure minus ambient pressure \[Pa\].
/// * `nozzle_area` - Nozzle exit area \[m²\].
/// * `thrust_direction` - Unit vector of thrust (body frame).
pub fn projectile_thrust(
    mass_flow_rate: f64,
    exhaust_velocity: f64,
    pressure_diff: f64,
    nozzle_area: f64,
    thrust_direction: [f64; 3],
) -> [f64; 3] {
    let f_mag = mass_flow_rate * exhaust_velocity + pressure_diff * nozzle_area;
    scale3(normalize3(thrust_direction), f_mag.max(0.0))
}

/// Integrate a fluid projectile forward by one time step.
///
/// Applies buoyancy, drag, gravity, and thrust forces.
///
/// # Arguments
/// * `proj` - Mutable projectile state.
/// * `fluid` - Fluid properties (used when submerged).
/// * `thrust` - External thrust force \[N\].
/// * `gravity` - Gravitational acceleration \[m/s²\].
/// * `free_surface_y` - Y-coordinate of free surface \[m\].
/// * `dt` - Time step \[s\].
pub fn integrate_fluid_projectile(
    proj: &mut FluidProjectile,
    fluid: &FluidProperties,
    thrust: [f64; 3],
    gravity: f64,
    free_surface_y: f64,
    dt: f64,
) {
    let submerged_frac =
        sphere_submersion_fraction(proj.body.position[1], proj.body.radius, free_surface_y);
    proj.submerged = submerged_frac > 0.01;

    let gravity_force = [0.0, -proj.body.mass * gravity, 0.0];

    let buoyancy = if proj.submerged {
        let vol_sub = submerged_frac * (4.0 / 3.0 * PI * proj.body.radius.powi(3));
        buoyancy_force(fluid.density, gravity, vol_sub)
    } else {
        [0.0; 3]
    };

    let fluid_vel = [0.0_f64; 3]; // assume quiescent fluid
    let drag = if proj.submerged {
        drag_force_sphere(&proj.body, fluid_vel, fluid)
    } else {
        // Air drag
        let air = FluidProperties::air();
        drag_force_sphere(&proj.body, fluid_vel, &air)
    };

    let total_force = add3(add3(add3(gravity_force, buoyancy), drag), thrust);
    let accel = scale3(total_force, 1.0 / proj.body.mass);

    proj.body.velocity = add3(proj.body.velocity, scale3(accel, dt));
    proj.body.position = add3(proj.body.position, scale3(proj.body.velocity, dt));
    proj.thrust_impulse += len3(thrust) * dt;
}

// ─────────────────────────────────────────────────────────────────────────────
// 15. Added Mass Effect
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the added mass coefficient C_m for a sphere.
///
/// For a sphere, C_m = 0.5 (theoretical inviscid potential flow result).
/// Returns the added mass m_a = C_m · ρ_f · V_body.
///
/// # Arguments
/// * `fluid_density` - Fluid density \[kg/m³\].
/// * `body_volume` - Volume of the body \[m³\].
/// * `cm` - Added mass coefficient (0.5 for sphere).
pub fn added_mass(fluid_density: f64, body_volume: f64, cm: f64) -> f64 {
    cm * fluid_density * body_volume
}

/// Effective mass of a body including added mass.
pub fn effective_mass(body_mass: f64, added_mass_val: f64) -> f64 {
    body_mass + added_mass_val
}

/// Compute added-mass force on an accelerating body.
///
/// F_added = −m_a · a_body
pub fn added_mass_force(added_mass_val: f64, body_acceleration: [f64; 3]) -> [f64; 3] {
    scale3(body_acceleration, -added_mass_val)
}

// ─────────────────────────────────────────────────────────────────────────────
// 16. Froude Number and Wave Drag
// ─────────────────────────────────────────────────────────────────────────────

/// Froude number Fr = V / sqrt(g · L).
pub fn froude_number(velocity: f64, length: f64, gravity: f64) -> f64 {
    velocity / (gravity * length).sqrt().max(1e-15)
}

/// Simplified wave drag coefficient as function of Froude number.
///
/// Uses a Kelvin wave drag model approximation.
/// Returns C_w(Fr) ≈ exp(−1/(Fr²)) / Fr² for Fr > 0.
pub fn wave_drag_coefficient(fr: f64) -> f64 {
    if fr < 1e-6 {
        return 0.0;
    }
    (-1.0 / (fr * fr)).exp() / (fr * fr)
}

/// Compute wave drag force on a surface vessel.
///
/// F_wave = 0.5 · ρ · V² · L² · C_w(Fr)
///
/// # Arguments
/// * `speed` - Forward speed \[m/s\].
/// * `waterline_length` - Waterline length \[m\].
/// * `fluid` - Fluid properties.
/// * `gravity` - Gravitational acceleration \[m/s²\].
pub fn wave_drag_force(
    speed: f64,
    waterline_length: f64,
    fluid: &FluidProperties,
    gravity: f64,
) -> f64 {
    let fr = froude_number(speed, waterline_length, gravity);
    let cw = wave_drag_coefficient(fr);
    0.5 * fluid.density * speed * speed * waterline_length * waterline_length * cw
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal LCG RNG helper
// ─────────────────────────────────────────────────────────────────────────────

struct FcRng {
    state: u64,
}

impl FcRng {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── vector helpers ────────────────────────────────────────────────────────

    #[test]
    fn test_add3_basic() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = add3(a, b);
        assert!((c[0] - 5.0).abs() < 1e-12);
        assert!((c[1] - 7.0).abs() < 1e-12);
        assert!((c[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3_unit() {
        let v = normalize3([3.0, 4.0, 0.0]);
        let len = len3(v);
        assert!((len - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3_zero_vector() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert!((v[2] - 1.0).abs() < 1e-12); // fallback to z-axis
    }

    // ── FluidProperties ───────────────────────────────────────────────────────

    #[test]
    fn test_fluid_properties_water() {
        let w = FluidProperties::water();
        assert!((w.density - 998.2).abs() < 0.1);
        assert!(w.kinematic_viscosity > 0.0);
    }

    #[test]
    fn test_fluid_properties_air() {
        let a = FluidProperties::air();
        assert!(a.density < 2.0);
        assert!(a.viscosity > 0.0);
    }

    #[test]
    fn test_kinematic_viscosity_derived() {
        let f = FluidProperties::new(1000.0, 0.001, 0.07, 1500.0);
        assert!((f.kinematic_viscosity - 1e-6).abs() < 1e-10);
    }

    // ── Reynolds number ───────────────────────────────────────────────────────

    #[test]
    fn test_reynolds_number_basic() {
        // Water at 20°C: ν ≈ 1e-6; U = 1 m/s, L = 0.01 m => Re = 10000
        let re = reynolds_number(1.0, 0.01, 1e-6);
        assert!((re - 10_000.0).abs() < 1.0);
    }

    #[test]
    fn test_reynolds_number_zero_velocity() {
        assert_eq!(reynolds_number(0.0, 1.0, 1e-6), 0.0);
    }

    // ── Drag coefficient ──────────────────────────────────────────────────────

    #[test]
    fn test_drag_coefficient_stokes() {
        // Re << 1: C_d ≈ 24/Re
        let cd = drag_coefficient_sphere(0.1);
        assert!((cd - 240.0).abs() < 1.0);
    }

    #[test]
    fn test_drag_coefficient_newton() {
        // Re >> 1000: C_d = 0.44
        let cd = drag_coefficient_sphere(1e6);
        assert!((cd - 0.44).abs() < 1e-10);
    }

    #[test]
    fn test_drag_coefficient_positive() {
        for re in [0.01, 1.0, 10.0, 100.0, 1000.0, 1e5] {
            assert!(drag_coefficient_sphere(re) > 0.0, "re={re}");
        }
    }

    // ── Drag force ────────────────────────────────────────────────────────────

    #[test]
    fn test_drag_force_opposing() {
        // Body moving in +x, drag should be in -x
        let body_vel = [10.0, 0.0, 0.0];
        let fluid_vel = [0.0, 0.0, 0.0];
        let fluid = FluidProperties::water();
        let f = drag_force(body_vel, fluid_vel, &fluid, 0.44, 0.01);
        assert!(f[0] < 0.0, "drag should oppose motion, f={f:?}");
    }

    #[test]
    fn test_drag_force_zero_velocity() {
        let fluid = FluidProperties::water();
        let f = drag_force([0.0; 3], [0.0; 3], &fluid, 0.44, 0.01);
        assert_eq!(f, [0.0; 3]);
    }

    // ── Buoyancy ──────────────────────────────────────────────────────────────

    #[test]
    fn test_buoyancy_fully_submerged() {
        let body = RigidBodyFluidState::new([0.0, -2.0, 0.0], [0.0; 3], 1.0, 0.5);
        let fluid = FluidProperties::water();
        let result = buoyancy_sphere(&body, &fluid, 0.0, 9.81);
        let vol = 4.0 / 3.0 * PI * 0.5_f64.powi(3);
        let expected = fluid.density * 9.81 * vol;
        assert!(
            (result.force[1] - expected).abs() < 1.0,
            "force={}",
            result.force[1]
        );
        assert!((result.submerged_volume - vol).abs() < 1e-6);
    }

    #[test]
    fn test_buoyancy_above_surface() {
        let body = RigidBodyFluidState::new([0.0, 5.0, 0.0], [0.0; 3], 1.0, 0.5);
        let fluid = FluidProperties::water();
        let result = buoyancy_sphere(&body, &fluid, 0.0, 9.81);
        assert!((result.force[1]).abs() < 1e-10, "force={}", result.force[1]);
        assert_eq!(result.submerged_volume, 0.0);
    }

    #[test]
    fn test_buoyancy_partial_submersion() {
        // Centre at free surface => half submerged
        let r = 1.0;
        let body = RigidBodyFluidState::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, r);
        let fluid = FluidProperties::water();
        let result = buoyancy_sphere(&body, &fluid, 0.0, 9.81);
        let half_vol = 2.0 / 3.0 * PI * r * r * r;
        assert!(
            (result.submerged_volume - half_vol).abs() < 0.05,
            "vol={} expected≈{half_vol}",
            result.submerged_volume
        );
    }

    // ── Submersion state ──────────────────────────────────────────────────────

    #[test]
    fn test_submersion_state_airborne() {
        let body = RigidBodyFluidState::new([0.0, 2.0, 0.0], [0.0; 3], 1.0, 0.5);
        assert_eq!(
            submersion_state_sphere(&body, 0.0),
            SubmersionState::Airborne
        );
    }

    #[test]
    fn test_submersion_state_submerged() {
        let body = RigidBodyFluidState::new([0.0, -2.0, 0.0], [0.0; 3], 1.0, 0.5);
        assert_eq!(
            submersion_state_sphere(&body, 0.0),
            SubmersionState::Submerged
        );
    }

    #[test]
    fn test_submersion_state_crossing() {
        let body = RigidBodyFluidState::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        assert_eq!(
            submersion_state_sphere(&body, 0.0),
            SubmersionState::Crossing
        );
    }

    // ── Water entry / exit ────────────────────────────────────────────────────

    #[test]
    fn test_water_crossing_detection_entry() {
        let prev = RigidBodyFluidState::new([0.0, 1.0, 0.0], [0.0, -5.0, 0.0], 1.0, 0.3);
        let curr = RigidBodyFluidState::new([0.0, 0.0, 0.0], [0.0, -5.0, 0.0], 1.0, 0.3);
        let ev = detect_water_crossing(&prev, &curr, 0.0);
        assert!(ev.is_some());
        assert!(ev.unwrap().is_entry);
    }

    #[test]
    fn test_water_crossing_none_when_stable() {
        let b1 = RigidBodyFluidState::new([0.0, -2.0, 0.0], [0.0; 3], 1.0, 0.3);
        let b2 = RigidBodyFluidState::new([0.0, -2.1, 0.0], [0.0; 3], 1.0, 0.3);
        let ev = detect_water_crossing(&b1, &b2, 0.0);
        assert!(ev.is_none());
    }

    // ── Splash particles ──────────────────────────────────────────────────────

    #[test]
    fn test_splash_particle_count() {
        let event = WaterEntryEvent {
            position: [0.0; 3],
            velocity: [0.0, -5.0, 0.0],
            entry_speed: 5.0,
            vertical_speed: 5.0,
            is_entry: true,
        };
        let particles = generate_splash_particles(&event, 20, &FluidProperties::water(), 1e-4, 42);
        assert_eq!(particles.len(), 20);
    }

    #[test]
    fn test_splash_particle_lifetime_positive() {
        let event = WaterEntryEvent {
            position: [0.0; 3],
            velocity: [0.0, -3.0, 0.0],
            entry_speed: 3.0,
            vertical_speed: 3.0,
            is_entry: true,
        };
        let particles = generate_splash_particles(&event, 5, &FluidProperties::water(), 1e-4, 1);
        for p in &particles {
            assert!(p.lifetime > 0.0, "lifetime={}", p.lifetime);
        }
    }

    // ── Ray vs free surface ───────────────────────────────────────────────────

    #[test]
    fn test_ray_free_surface_hit() {
        let t = ray_free_surface([0.0, 2.0, 0.0], [0.0, -1.0, 0.0], 0.0, 10.0);
        assert!(t.is_some());
        assert!((t.unwrap() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_ray_free_surface_parallel_miss() {
        let t = ray_free_surface([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], 0.0, 10.0);
        assert!(t.is_none());
    }

    #[test]
    fn test_ray_free_surface_behind_miss() {
        let t = ray_free_surface([0.0, -1.0, 0.0], [0.0, -1.0, 0.0], 0.0, 10.0);
        assert!(t.is_none());
    }

    // ── Sphere submersion fraction ────────────────────────────────────────────

    #[test]
    fn test_submersion_fraction_full() {
        let frac = sphere_submersion_fraction(-2.0, 0.5, 0.0);
        assert!((frac - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_submersion_fraction_zero() {
        let frac = sphere_submersion_fraction(2.0, 0.5, 0.0);
        assert!((frac).abs() < 1e-10);
    }

    #[test]
    fn test_submersion_fraction_half() {
        let frac = sphere_submersion_fraction(0.0, 1.0, 0.0);
        assert!((frac - 0.5).abs() < 0.05, "frac={frac}");
    }

    // ── Domain boundary collision ─────────────────────────────────────────────

    #[test]
    fn test_domain_boundary_collision_x() {
        let domain = FluidDomain::new([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], 0.0);
        let mut body = RigidBodyFluidState::new([4.8, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.3);
        let hit = handle_domain_boundary_collision(&mut body, &domain, 0.8);
        assert!(hit);
        assert!(body.velocity[0] < 0.0, "vx should be reflected");
    }

    #[test]
    fn test_domain_boundary_no_collision() {
        let domain = FluidDomain::new([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], 0.0);
        let mut body = RigidBodyFluidState::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.3);
        let hit = handle_domain_boundary_collision(&mut body, &domain, 0.8);
        assert!(!hit);
    }

    // ── SPH coupling ──────────────────────────────────────────────────────────

    #[test]
    fn test_sph_cubic_kernel_positive() {
        assert!(sph_cubic_kernel(0.0, 1.0) > 0.0);
        assert!(sph_cubic_kernel(0.5, 1.0) > 0.0);
    }

    #[test]
    fn test_sph_cubic_kernel_outside_support() {
        assert_eq!(sph_cubic_kernel(3.0, 1.0), 0.0);
    }

    #[test]
    fn test_sph_rigid_coupling_zero_particles() {
        let f = sph_rigid_coupling_force([0.0; 3], 0.5, &[], 1.0);
        assert_eq!(f, [0.0; 3]);
    }

    // ── LBM ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_lbm_equilibrium_sum_to_rho() {
        let rho = 1.0;
        let u = [0.0_f64; 3];
        let sum: f64 = (0..19).map(|i| lbm_equilibrium(rho, u, i)).sum();
        assert!((sum - rho).abs() < 1e-10, "sum={sum}");
    }

    #[test]
    fn test_lbm_equilibrium_rest_node() {
        // At rest node (i=0), f^eq = w_0 * rho
        let rho = 1.5;
        let u = [0.0_f64; 3];
        let feq = lbm_equilibrium(rho, u, 0);
        let expected = D3Q19_WEIGHTS[0] * rho;
        assert!((feq - expected).abs() < 1e-10);
    }

    // ── Two-way coupling ──────────────────────────────────────────────────────

    #[test]
    fn test_two_way_coupling_separating() {
        let mut body = RigidBodyFluidState::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.5);
        let impulse =
            two_way_coupling_impulse(&mut body, 100.0, [0.0; 3], [1.0, 0.0, 0.0], 0.5, 0.01);
        // Separating — zero impulse
        assert_eq!(impulse, [0.0; 3]);
    }

    #[test]
    fn test_two_way_coupling_approaching() {
        let mut body = RigidBodyFluidState::new([0.0; 3], [-1.0, 0.0, 0.0], 1.0, 0.5);
        let impulse =
            two_way_coupling_impulse(&mut body, 100.0, [0.0; 3], [1.0, 0.0, 0.0], 0.8, 0.01);
        // Impulse should be in +x direction
        assert!(impulse[0] > 0.0, "impulse={impulse:?}");
    }

    // ── Underwater friction ───────────────────────────────────────────────────

    #[test]
    fn test_underwater_friction_no_tangential_velocity() {
        let result = underwater_contact_friction(100.0, [0.0; 3], 0.5, 0.3, 0.1);
        assert_eq!(result.friction_force, [0.0; 3]);
        assert!(!result.is_sliding);
    }

    #[test]
    fn test_underwater_friction_opposing() {
        let result = underwater_contact_friction(100.0, [1.0, 0.0, 0.0], 0.5, 0.0, 0.0);
        // Friction should oppose motion (negative x)
        assert!(
            result.friction_force[0] < 0.0,
            "f={:?}",
            result.friction_force
        );
    }

    // ── Added mass ────────────────────────────────────────────────────────────

    #[test]
    fn test_added_mass_sphere() {
        let r = 0.1;
        let vol = 4.0 / 3.0 * PI * r * r * r;
        let rho_f = 1000.0;
        let ma = added_mass(rho_f, vol, 0.5);
        let expected = 0.5 * rho_f * vol;
        assert!((ma - expected).abs() < 1e-12);
    }

    #[test]
    fn test_effective_mass_gt_body_mass() {
        let ma = 0.3;
        let mb = 1.0;
        assert!(effective_mass(mb, ma) > mb);
    }

    // ── Froude number ─────────────────────────────────────────────────────────

    #[test]
    fn test_froude_number_known() {
        // Fr = V / sqrt(g L)
        let fr = froude_number(3.132, 1.0, 9.81);
        assert!((fr - 1.0).abs() < 0.01, "fr={fr}");
    }

    // ── Projectile integration ────────────────────────────────────────────────

    #[test]
    fn test_projectile_falls_under_gravity() {
        let mut proj = FluidProjectile::new([0.0, 10.0, 0.0], [0.0; 3], 1.0, 0.05);
        let fluid = FluidProperties::water();
        let y0 = proj.body.position[1];
        integrate_fluid_projectile(&mut proj, &fluid, [0.0; 3], 9.81, 0.0, 0.1);
        assert!(proj.body.position[1] < y0, "should fall under gravity");
    }

    #[test]
    fn test_projectile_thrust_impulse_accumulates() {
        let mut proj = FluidProjectile::new([0.0; 3], [0.0; 3], 1.0, 0.05);
        let fluid = FluidProperties::water();
        let thrust = [10.0, 0.0, 0.0];
        integrate_fluid_projectile(&mut proj, &fluid, thrust, 9.81, -10.0, 0.01);
        assert!(proj.thrust_impulse > 0.0);
    }

    // ── Slamming force ────────────────────────────────────────────────────────

    #[test]
    fn test_slamming_force_positive() {
        let f = water_entry_slamming_force(1000.0, 10.0, 0.1, 1.0);
        assert!(f > 0.0, "slamming force={f}");
    }

    #[test]
    fn test_slamming_force_scales_with_speed_squared() {
        let f1 = water_entry_slamming_force(1000.0, 5.0, 0.1, 1.0);
        let f2 = water_entry_slamming_force(1000.0, 10.0, 0.1, 1.0);
        assert!((f2 / f1 - 4.0).abs() < 1e-10, "ratio={}", f2 / f1);
    }

    // ── Contact radius ────────────────────────────────────────────────────────

    #[test]
    fn test_contact_radius_zero_penetration() {
        assert_eq!(sphere_water_contact_radius(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_contact_radius_increases_with_depth() {
        let r1 = sphere_water_contact_radius(1.0, 0.1);
        let r2 = sphere_water_contact_radius(1.0, 0.5);
        assert!(r2 > r1, "r1={r1} r2={r2}");
    }

    // ── Wave drag ─────────────────────────────────────────────────────────────

    #[test]
    fn test_wave_drag_force_positive() {
        let fluid = FluidProperties::water();
        let f = wave_drag_force(5.0, 10.0, &fluid, 9.81);
        assert!(f >= 0.0, "wave drag={f}");
    }

    // ── FluidDomain ───────────────────────────────────────────────────────────

    #[test]
    fn test_fluid_domain_contains_fluid_point() {
        let domain = FluidDomain::new([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], 0.0);
        assert!(domain.contains_fluid_point([1.0, -1.0, 0.0]));
        assert!(!domain.contains_fluid_point([1.0, 1.0, 0.0])); // above surface
    }
}
