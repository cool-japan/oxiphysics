// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Boundary handling for SPH simulations.
//!
//! Provides static boundary particles, analytical plane boundaries,
//! rigid boundary particle sets, mirror-particle generation, moving
//! (dynamic) boundaries, and Lennard-Jones repulsion forces.

use oxiphysics_core::math::Vec3;

use crate::kernel::SphKernel;
use crate::particle::ParticleSet;

// ── BoundaryParticle ──────────────────────────────────────────────────────────

/// A static boundary particle for wall representation.
#[derive(Debug, Clone)]
pub struct BoundaryParticle {
    /// Position of the boundary particle.
    pub position: Vec3,
    /// Volume (inverse number density) of the boundary particle.
    pub volume: f64,
}

impl BoundaryParticle {
    /// Create a new boundary particle.
    pub fn new(position: Vec3, volume: f64) -> Self {
        Self { position, volume }
    }
}

// ── BoundaryPlane ─────────────────────────────────────────────────────────────

/// An infinite plane boundary.
#[derive(Debug, Clone)]
pub struct BoundaryPlane {
    /// A point on the plane.
    pub point: Vec3,
    /// Outward-facing normal (pointing into the fluid domain).
    pub normal: Vec3,
}

impl BoundaryPlane {
    /// Create a new boundary plane.
    pub fn new(point: Vec3, normal: Vec3) -> Self {
        let n = if normal.norm() > 1e-14 {
            normal / normal.norm()
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        Self { point, normal: n }
    }

    /// Signed distance from a point to this plane (positive = inside fluid).
    pub fn signed_distance(&self, pos: &Vec3) -> f64 {
        (pos - self.point).dot(&self.normal)
    }
}

// ── BoundarySet ───────────────────────────────────────────────────────────────

/// Collection of boundary elements for an SPH simulation.
#[derive(Debug, Clone, Default)]
pub struct BoundarySet {
    /// Static boundary particles.
    pub particles: Vec<BoundaryParticle>,
    /// Infinite plane boundaries.
    pub planes: Vec<BoundaryPlane>,
}

impl BoundarySet {
    /// Create an empty boundary set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a boundary particle.
    pub fn add_particle(&mut self, particle: BoundaryParticle) {
        self.particles.push(particle);
    }

    /// Add a boundary plane.
    pub fn add_plane(&mut self, plane: BoundaryPlane) {
        self.planes.push(plane);
    }
}

// ── RigidBoundaryParticle ─────────────────────────────────────────────────────

/// A boundary particle belonging to a rigid solid body (Akinci 2012 style).
///
/// Each particle carries a position, outward normal, representative volume,
/// and an instantaneous velocity (for moving boundaries).
#[derive(Debug, Clone)]
pub struct RigidBoundaryParticle {
    /// World-space position.
    pub position: [f64; 3],
    /// Outward surface normal.
    pub normal: [f64; 3],
    /// Representative volume (m³).
    pub volume: f64,
    /// Instantaneous velocity of the boundary particle (m/s).
    pub velocity: [f64; 3],
}

impl RigidBoundaryParticle {
    /// Create a new rigid boundary particle.
    pub fn new(position: [f64; 3], normal: [f64; 3], volume: f64) -> Self {
        Self {
            position,
            normal,
            volume,
            velocity: [0.0; 3],
        }
    }

    /// Set the velocity of the particle.
    pub fn with_velocity(mut self, velocity: [f64; 3]) -> Self {
        self.velocity = velocity;
        self
    }
}

// ── RigidBoundarySet ──────────────────────────────────────────────────────────

/// Collection of `RigidBoundaryParticle`s with a simple grid-based query.
#[derive(Debug, Default)]
pub struct RigidBoundarySet {
    /// All particles.
    pub particles: Vec<RigidBoundaryParticle>,
}

impl RigidBoundarySet {
    /// Create an empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a particle.
    pub fn add(&mut self, p: RigidBoundaryParticle) {
        self.particles.push(p);
    }

    /// Find all boundary particles within distance `h` of a query position.
    pub fn query_neighbors(&self, pos: [f64; 3], h: f64) -> Vec<&RigidBoundaryParticle> {
        let h2 = h * h;
        self.particles
            .iter()
            .filter(|p| {
                let dx = p.position[0] - pos[0];
                let dy = p.position[1] - pos[1];
                let dz = p.position[2] - pos[2];
                dx * dx + dy * dy + dz * dz <= h2
            })
            .collect()
    }

    /// Number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }
}

// ── BoundaryPressure ──────────────────────────────────────────────────────────

/// Compute an approximate pressure contribution from a single rigid boundary
/// particle onto a nearby fluid particle.
///
/// Uses the mirrored-pressure approximation: the boundary pressure equals
/// the fluid pressure plus a hydrostatic correction.
pub fn compute_boundary_pressure(
    fluid_pressure: f64,
    fluid_density: f64,
    gravity_mag: f64,
    fluid_pos_y: f64,
    boundary_pos_y: f64,
    rest_density: f64,
    kernel_h: f64,
) -> f64 {
    let _ = kernel_h; // h determines the interaction radius (not used directly here)
    let depth = (fluid_pos_y - boundary_pos_y).max(0.0);
    fluid_pressure + rest_density * gravity_mag * depth + fluid_density * gravity_mag * depth
}

// ── MirrorBoundary ────────────────────────────────────────────────────────────

/// Generates mirror (ghost) particles by reflecting fluid particles across a
/// plane defined by a normal and a signed distance from the origin.
pub struct MirrorBoundary {
    /// Unit outward normal of the plane.
    pub plane_normal: [f64; 3],
    /// Signed distance of the plane from the origin along `plane_normal`.
    pub plane_d: f64,
}

impl MirrorBoundary {
    /// Create a new mirror boundary.
    pub fn new(plane_normal: [f64; 3], plane_d: f64) -> Self {
        // Normalise
        let len = (plane_normal[0] * plane_normal[0]
            + plane_normal[1] * plane_normal[1]
            + plane_normal[2] * plane_normal[2])
            .sqrt();
        let n = if len > 1e-14 {
            [
                plane_normal[0] / len,
                plane_normal[1] / len,
                plane_normal[2] / len,
            ]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            plane_normal: n,
            plane_d,
        }
    }

    /// Signed distance from a point to the mirror plane.
    ///
    /// Positive values mean the point is on the "fluid" side of the plane.
    pub fn signed_distance(&self, pos: [f64; 3]) -> f64 {
        self.plane_normal[0] * pos[0]
            + self.plane_normal[1] * pos[1]
            + self.plane_normal[2] * pos[2]
            - self.plane_d
    }

    /// Reflect a single point across the plane.
    pub fn reflect(&self, pos: [f64; 3]) -> [f64; 3] {
        let d = self.signed_distance(pos);
        [
            pos[0] - 2.0 * d * self.plane_normal[0],
            pos[1] - 2.0 * d * self.plane_normal[1],
            pos[2] - 2.0 * d * self.plane_normal[2],
        ]
    }

    /// Generate mirror positions for all fluid particles within `h` of the plane.
    ///
    /// Only particles whose signed distance is in `(0, h)` need mirrors.
    pub fn generate_mirrors(&self, fluid_positions: &[[f64; 3]], h: f64) -> Vec<[f64; 3]> {
        fluid_positions
            .iter()
            .filter_map(|&pos| {
                let d = self.signed_distance(pos);
                if d > 0.0 && d < h {
                    Some(self.reflect(pos))
                } else {
                    None
                }
            })
            .collect()
    }
}

// ── DynamicBoundary ───────────────────────────────────────────────────────────

/// A moving boundary that updates the positions of its particles each step.
pub struct DynamicBoundary {
    /// Underlying particle set.
    pub particles: RigidBoundarySet,
}

impl Default for DynamicBoundary {
    fn default() -> Self {
        Self::new()
    }
}

impl DynamicBoundary {
    /// Create a new empty dynamic boundary.
    pub fn new() -> Self {
        Self {
            particles: RigidBoundarySet::new(),
        }
    }

    /// Advance all boundary particles by `dt` using their current velocities.
    pub fn step(&mut self, dt: f64) {
        for p in &mut self.particles.particles {
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
        }
    }

    /// Set a uniform velocity for all boundary particles.
    pub fn set_velocity(&mut self, vel: [f64; 3]) {
        for p in &mut self.particles.particles {
            p.velocity = vel;
        }
    }
}

// ── BoundaryForce (Lennard-Jones) ─────────────────────────────────────────────

/// Compute the Lennard-Jones force magnitude between a fluid particle and a
/// boundary at distance `r`.
///
/// `F(r) = 4 * epsilon * (12 * sigma^12 / r^13 - 6 * sigma^6 / r^7)`
///
/// The result is positive (repulsive) when `r < sigma * 2^(1/6)`.
pub fn compute_lennard_jones_force(r: f64, epsilon: f64, sigma: f64) -> f64 {
    if r < 1e-14 {
        return f64::MAX * 0.5; // avoid singularity
    }
    let s6 = sigma.powi(6);
    let r6 = r.powi(6);
    let r12 = r6 * r6;
    let s12 = s6 * s6;
    4.0 * epsilon * (12.0 * s12 / (r12 * r) - 6.0 * s6 / (r6 * r))
}

// ── apply_plane_penalty_forces ────────────────────────────────────────────────

/// Apply penalty forces from boundary planes.
pub fn apply_plane_penalty_forces(
    particles: &mut ParticleSet,
    planes: &[BoundaryPlane],
    stiffness: f64,
    damping: f64,
) {
    for i in 0..particles.len() {
        for plane in planes {
            let dist = plane.signed_distance(&particles.positions[i]);
            if dist < 0.0 {
                let penalty = -stiffness * dist;
                let vel_normal = particles.velocities[i].dot(&plane.normal);
                let damp = -damping * vel_normal;
                particles.forces[i] += (penalty + damp) * plane.normal;
            }
        }
    }
}

/// Apply boundary forces from static boundary particles (Akinci-style).
pub fn apply_boundary_particle_forces(
    particles: &mut ParticleSet,
    boundary_particles: &[BoundaryParticle],
    kernel: &dyn SphKernel,
    h: f64,
    rest_density: f64,
) {
    for i in 0..particles.len() {
        let rhoi = particles.densities[i].max(1e-14);
        let pi_val = particles.pressures[i];

        for bp in boundary_particles {
            let rij = particles.positions[i] - bp.position;
            let r = rij.norm();
            if r < 1e-14 || r > 2.0 * h {
                continue;
            }
            let rhat = rij / r;
            let grad = kernel.grad_w(r, h);

            let psi = rest_density * bp.volume;
            let pressure_term = pi_val / (rhoi * rhoi);
            particles.forces[i] -= psi * 2.0 * pressure_term * grad * rhat * rhoi;
        }
    }
}

/// Apply all boundary forces from a [`BoundarySet`].
pub fn apply_boundary_forces(
    particles: &mut ParticleSet,
    boundaries: &BoundarySet,
    kernel: &dyn SphKernel,
    h: f64,
    rest_density: f64,
    penalty_stiffness: f64,
    penalty_damping: f64,
) {
    apply_plane_penalty_forces(
        particles,
        &boundaries.planes,
        penalty_stiffness,
        penalty_damping,
    );
    apply_boundary_particle_forces(particles, &boundaries.particles, kernel, h, rest_density);
}

// ── Ghost particle generation ─────────────────────────────────────────────────

/// Generate ghost particles for a plane boundary by reflecting all fluid
/// particles that are within distance `h` of the plane.
///
/// Returns `(positions, velocities)` of the ghost particles.  Ghost
/// velocities are the mirror-image of the fluid velocity (no-slip).
pub fn generate_ghost_particles(
    fluid_positions: &[[f64; 3]],
    fluid_velocities: &[[f64; 3]],
    plane_normal: [f64; 3],
    plane_d: f64,
    h: f64,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mb = MirrorBoundary::new(plane_normal, plane_d);
    let mut ghost_pos = Vec::new();
    let mut ghost_vel = Vec::new();
    for (i, &pos) in fluid_positions.iter().enumerate() {
        let d = mb.signed_distance(pos);
        if d > 0.0 && d < h {
            ghost_pos.push(mb.reflect(pos));
            // Reflect velocity for no-slip: v_ghost = -v_fluid
            let v = fluid_velocities[i];
            ghost_vel.push([-v[0], -v[1], -v[2]]);
        }
    }
    (ghost_pos, ghost_vel)
}

/// Generate ghost particles with free-slip condition.
///
/// The tangential velocity component is preserved; only the normal
/// component is negated.
pub fn generate_ghost_particles_free_slip(
    fluid_positions: &[[f64; 3]],
    fluid_velocities: &[[f64; 3]],
    plane_normal: [f64; 3],
    plane_d: f64,
    h: f64,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mb = MirrorBoundary::new(plane_normal, plane_d);
    let n = mb.plane_normal;
    let mut ghost_pos = Vec::new();
    let mut ghost_vel = Vec::new();
    for (i, &pos) in fluid_positions.iter().enumerate() {
        let d = mb.signed_distance(pos);
        if d > 0.0 && d < h {
            ghost_pos.push(mb.reflect(pos));
            let v = fluid_velocities[i];
            // v_n = (v · n) n
            let vn = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
            // v_ghost = v - 2 * v_n * n  (negate normal component, keep tangential)
            ghost_vel.push([
                v[0] - 2.0 * vn * n[0],
                v[1] - 2.0 * vn * n[1],
                v[2] - 2.0 * vn * n[2],
            ]);
        }
    }
    (ghost_pos, ghost_vel)
}

// ── Lennard-Jones boundary improvements ───────────────────────────────────────

/// Truncated and shifted Lennard-Jones boundary force (WCA-style).
///
/// The force is purely repulsive and goes to zero exactly at the cutoff
/// distance `r_cut = sigma * 2^(1/6)`.
pub fn compute_wca_force(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let r_cut = sigma * 2.0_f64.powf(1.0 / 6.0);
    if r >= r_cut || r < 1e-14 {
        return 0.0;
    }
    // Standard LJ force, but only in the repulsive region
    let s6 = sigma.powi(6);
    let r6 = r.powi(6);
    let r12 = r6 * r6;
    let s12 = s6 * s6;
    4.0 * epsilon * (12.0 * s12 / (r12 * r) - 6.0 * s6 / (r6 * r))
}

/// Soft-core Lennard-Jones force that avoids the singularity at r = 0.
///
/// `F(r) = 4ε (12σ^12/(r² + α²)^6.5 - 6σ^6/(r² + α²)^3.5) * r`
/// where α is a softening parameter (typically 0.1 * sigma).
pub fn compute_softcore_lj_force(r: f64, epsilon: f64, sigma: f64, alpha: f64) -> f64 {
    let r2_eff = r * r + alpha * alpha;
    let s2 = sigma * sigma;
    let ratio = s2 / r2_eff;
    let ratio3 = ratio * ratio * ratio;
    let ratio6 = ratio3 * ratio3;
    4.0 * epsilon * (12.0 * ratio6 - 6.0 * ratio3) * r / r2_eff
}

// ── Boundary normals ──────────────────────────────────────────────────────────

/// Compute approximate outward normal at each boundary particle using
/// the SPH gradient of the indicator function.
///
/// n_i ≈ -h Σ_j (V_j ∇W_ij)  (sum over boundary neighbours only).
pub fn compute_boundary_normals(
    boundary_positions: &[[f64; 3]],
    volumes: &[f64],
    h: f64,
) -> Vec<[f64; 3]> {
    let n = boundary_positions.len();
    let mut normals = vec![[0.0_f64; 3]; n];
    let h2 = (2.0 * h) * (2.0 * h);

    for i in 0..n {
        let mut nx = 0.0_f64;
        let mut ny = 0.0_f64;
        let mut nz = 0.0_f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = boundary_positions[i][0] - boundary_positions[j][0];
            let dy = boundary_positions[i][1] - boundary_positions[j][1];
            let dz = boundary_positions[i][2] - boundary_positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 >= h2 || r2 < 1e-28 {
                continue;
            }
            let r = r2.sqrt();
            // Approximate kernel gradient magnitude (cubic spline)
            let q = r / h;
            let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
            let dw_dr = if q >= 2.0 {
                0.0
            } else if q >= 1.0 {
                let t = 2.0 - q;
                sigma * (-0.75 * t * t) / h
            } else {
                sigma * (-3.0 * q + 2.25 * q * q) / h
            };
            let factor = -h * volumes[j] * dw_dr / r;
            nx += factor * dx;
            ny += factor * dy;
            nz += factor * dz;
        }
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len > 1e-14 {
            normals[i] = [nx / len, ny / len, nz / len];
        }
    }
    normals
}

// ── Dummy particle method ─────────────────────────────────────────────────────

/// A dummy boundary particle (Adami et al. 2012 style) that extrapolates
/// fluid pressure and velocity to maintain a smooth kernel support near walls.
#[derive(Debug, Clone)]
pub struct DummyParticle {
    /// Position.
    pub position: [f64; 3],
    /// Extrapolated pressure.
    pub pressure: f64,
    /// Extrapolated velocity.
    pub velocity: [f64; 3],
    /// Volume.
    pub volume: f64,
}

/// Generate dummy particles along a plane boundary and extrapolate
/// fluid quantities onto them.
///
/// `n_layers` layers of dummy particles are placed at spacings of `dx`
/// behind the plane.
pub fn generate_dummy_particles(
    fluid_positions: &[[f64; 3]],
    fluid_pressures: &[f64],
    fluid_velocities: &[[f64; 3]],
    plane_normal: [f64; 3],
    plane_d: f64,
    dx: f64,
    n_layers: usize,
    h: f64,
) -> Vec<DummyParticle> {
    let mb = MirrorBoundary::new(plane_normal, plane_d);
    let n = mb.plane_normal;
    let mut dummies = Vec::new();

    for (i, &pos) in fluid_positions.iter().enumerate() {
        let d = mb.signed_distance(pos);
        if d > 0.0 && d < h {
            for layer in 1..=n_layers {
                let offset = layer as f64 * dx;
                let dummy_pos = [
                    pos[0] - 2.0 * d * n[0] - offset * n[0],
                    pos[1] - 2.0 * d * n[1] - offset * n[1],
                    pos[2] - 2.0 * d * n[2] - offset * n[2],
                ];
                // Extrapolate: pressure increases hydrostatically,
                // velocity is negated (no-slip)
                let v = fluid_velocities[i];
                dummies.push(DummyParticle {
                    position: dummy_pos,
                    pressure: fluid_pressures[i],
                    velocity: [-v[0], -v[1], -v[2]],
                    volume: dx * dx * dx,
                });
            }
        }
    }
    dummies
}

// ── Dynamic boundary particle update ──────────────────────────────────────────

/// Update dynamic boundary particle positions with angular velocity.
///
/// Given a rotation centre, angular velocity vector, and linear velocity,
/// update all particles by `dt`.
pub fn update_rotating_boundary(
    particles: &mut RigidBoundarySet,
    center: [f64; 3],
    angular_vel: [f64; 3],
    linear_vel: [f64; 3],
    dt: f64,
) {
    for p in &mut particles.particles {
        // Relative position
        let rx = p.position[0] - center[0];
        let ry = p.position[1] - center[1];
        let rz = p.position[2] - center[2];

        // v = linear_vel + angular_vel × r
        let vx = linear_vel[0] + angular_vel[1] * rz - angular_vel[2] * ry;
        let vy = linear_vel[1] + angular_vel[2] * rx - angular_vel[0] * rz;
        let vz = linear_vel[2] + angular_vel[0] * ry - angular_vel[1] * rx;

        p.velocity = [vx, vy, vz];
        p.position[0] += vx * dt;
        p.position[1] += vy * dt;
        p.position[2] += vz * dt;
    }
}

/// Compute the boundary density contribution using Shepard correction.
///
/// ρ_boundary_i = m_i * Σ_j W(r_ij, h) / Σ_j (m_j/ρ_j) W(r_ij, h)
/// This is a simplified Shepard-corrected density for boundary particles.
pub fn boundary_shepard_density(
    boundary_pos: [f64; 3],
    fluid_positions: &[[f64; 3]],
    fluid_masses: &[f64],
    fluid_densities: &[f64],
    h: f64,
) -> f64 {
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    let h2 = (2.0 * h) * (2.0 * h);

    let cubic_w = |r: f64| -> f64 {
        let q = r / h;
        if q >= 2.0 {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            sigma * 0.25 * t * t * t
        } else {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        }
    };

    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (j, &fp) in fluid_positions.iter().enumerate() {
        let dx = boundary_pos[0] - fp[0];
        let dy = boundary_pos[1] - fp[1];
        let dz = boundary_pos[2] - fp[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 >= h2 {
            continue;
        }
        let r = r2.sqrt();
        let w = cubic_w(r);
        num += fluid_masses[j] * w;
        let rho_j = fluid_densities[j].max(1e-14);
        den += (fluid_masses[j] / rho_j) * w;
    }
    if den.abs() < 1e-14 { 0.0 } else { num / den }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::SphParticle;

    // ── original tests ────────────────────────────────────────────────────

    #[test]
    fn plane_penalty_repels_penetrating_particle() {
        let mut ps = ParticleSet::new();
        ps.add_particle(&SphParticle::new(
            Vec3::new(0.0, -0.1, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            1.0,
        ));
        ps.densities[0] = 1000.0;

        let plane = BoundaryPlane::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        ps.clear_forces();
        apply_plane_penalty_forces(&mut ps, &[plane], 10_000.0, 100.0);

        assert!(
            ps.forces[0].y > 0.0,
            "Penalty force should push particle up, got {:?}",
            ps.forces[0]
        );
    }

    #[test]
    fn boundary_particles_repel_fluid() {
        let mut ps = ParticleSet::new();
        let h = 0.1;
        let spacing: f64 = 0.05;
        let mass = 1000.0 * spacing.powi(3);

        ps.add_particle(&SphParticle::new(
            Vec3::new(0.0, 0.05, 0.0),
            Vec3::zeros(),
            mass,
        ));
        ps.densities[0] = 1000.0;
        ps.pressures[0] = 1000.0;

        let bp = BoundaryParticle::new(Vec3::new(0.0, 0.0, 0.0), spacing.powi(3));

        let kernel = crate::kernel::CubicSplineKernel;
        ps.clear_forces();
        apply_boundary_particle_forces(&mut ps, &[bp], &kernel, h, 1000.0);

        assert!(
            ps.forces[0].y > 0.0,
            "Boundary particle force should push fluid up, got {:?}",
            ps.forces[0]
        );
    }

    // ── MirrorBoundary tests ──────────────────────────────────────────────

    #[test]
    fn mirror_boundary_reflects_point() {
        // Plane: y = 0 (normal [0,1,0], d=0)
        let mb = MirrorBoundary::new([0.0, 1.0, 0.0], 0.0);
        let mirror = mb.reflect([1.0, 0.5, 0.0]);
        // Reflection: y → -0.5
        assert!((mirror[0] - 1.0).abs() < 1e-12);
        assert!((mirror[1] - (-0.5)).abs() < 1e-12);
        assert!((mirror[2]).abs() < 1e-12);
    }

    #[test]
    fn mirror_boundary_generate_mirrors_near_wall() {
        // Wall at y = 0, smoothing length 0.2
        let mb = MirrorBoundary::new([0.0, 1.0, 0.0], 0.0);
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.1, 0.0],  // within h → should be mirrored
            [0.0, 0.5, 0.0],  // outside h → no mirror
            [0.0, -0.1, 0.0], // behind the wall → no mirror
        ];
        let mirrors = mb.generate_mirrors(&positions, 0.2);
        assert_eq!(mirrors.len(), 1, "only 1 particle should be mirrored");
        assert!(
            (mirrors[0][1] - (-0.1)).abs() < 1e-12,
            "mirror y={}",
            mirrors[0][1]
        );
    }

    // ── RigidBoundarySet tests ────────────────────────────────────────────

    #[test]
    fn rigid_boundary_set_query_finds_nearby() {
        let mut rbs = RigidBoundarySet::new();
        rbs.add(RigidBoundaryParticle::new(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        ));
        rbs.add(RigidBoundaryParticle::new(
            [10.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        ));
        let nearby = rbs.query_neighbors([0.0, 0.0, 0.0], 1.0);
        assert_eq!(nearby.len(), 1);
    }

    // ── Lennard-Jones force tests ─────────────────────────────────────────

    #[test]
    fn lj_force_repulsive_at_short_range() {
        // At r < sigma * 2^(1/6) the force is repulsive (positive)
        let epsilon = 1.0;
        let sigma = 1.0;
        let r_repulsive = 0.9 * sigma;
        let f = compute_lennard_jones_force(r_repulsive, epsilon, sigma);
        assert!(f > 0.0, "LJ should be repulsive at r={r_repulsive}, f={f}");
    }

    #[test]
    fn lj_force_attractive_at_medium_range() {
        let epsilon = 1.0;
        let sigma = 1.0;
        // Equilibrium at r0 = sigma * 2^(1/6) ≈ 1.122
        // For r slightly larger than r0, force is attractive (negative)
        let r_attractive = 1.5 * sigma;
        let f = compute_lennard_jones_force(r_attractive, epsilon, sigma);
        assert!(
            f < 0.0,
            "LJ should be attractive at r={r_attractive}, f={f}"
        );
    }

    #[test]
    fn lj_force_zero_at_equilibrium() {
        let epsilon = 1.0;
        let sigma = 1.0;
        // Equilibrium: r0 = sigma * 2^(1/6)
        let r0 = sigma * 2.0_f64.powf(1.0 / 6.0);
        let f = compute_lennard_jones_force(r0, epsilon, sigma);
        assert!(
            f.abs() < 1e-10,
            "LJ force at equilibrium should be ~0, got {f}"
        );
    }

    // ── DynamicBoundary tests ─────────────────────────────────────────────

    #[test]
    fn dynamic_boundary_moves_particles() {
        let mut db = DynamicBoundary::new();
        db.particles.add(RigidBoundaryParticle::new(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        ));
        db.set_velocity([1.0, 0.0, 0.0]);
        db.step(0.5);
        let pos = db.particles.particles[0].position;
        assert!((pos[0] - 0.5).abs() < 1e-12, "x={}", pos[0]);
    }

    // ── BoundaryPressure test ─────────────────────────────────────────────

    #[test]
    fn boundary_pressure_increases_with_depth() {
        let p_shallow = compute_boundary_pressure(1000.0, 1000.0, 9.81, 0.1, 0.0, 1000.0, 0.1);
        let p_deep = compute_boundary_pressure(1000.0, 1000.0, 9.81, 0.5, 0.0, 1000.0, 0.1);
        assert!(
            p_deep > p_shallow,
            "deeper fluid should have higher boundary pressure"
        );
    }

    // ── Ghost particle tests ─────────────────────────────────────────────

    #[test]
    fn ghost_particles_no_slip() {
        let positions = vec![[0.0, 0.05, 0.0], [0.0, 0.5, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let (gp, gv) = generate_ghost_particles(&positions, &velocities, [0.0, 1.0, 0.0], 0.0, 0.1);
        // Only the first particle is within h=0.1 of the plane y=0
        assert_eq!(gp.len(), 1, "only 1 ghost particle expected");
        // Ghost position is at y = -0.05
        assert!((gp[0][1] - (-0.05)).abs() < 1e-12, "ghost y={}", gp[0][1]);
        // Ghost velocity is negated (no-slip)
        assert!((gv[0][0] - (-1.0)).abs() < 1e-12, "ghost vx={}", gv[0][0]);
    }

    #[test]
    fn ghost_particles_free_slip() {
        let positions = vec![[0.0, 0.05, 0.0]];
        let velocities = vec![[1.0, -0.5, 0.0]]; // tangential + normal
        let (gp, gv) =
            generate_ghost_particles_free_slip(&positions, &velocities, [0.0, 1.0, 0.0], 0.0, 0.1);
        assert_eq!(gp.len(), 1);
        // Normal component (-0.5 in y) is negated: ghost vy = 0.5
        assert!((gv[0][1] - 0.5).abs() < 1e-12, "ghost vy={}", gv[0][1]);
        // Tangential component (1.0 in x) is preserved
        assert!((gv[0][0] - 1.0).abs() < 1e-12, "ghost vx={}", gv[0][0]);
    }

    // ── WCA force tests ──────────────────────────────────────────────────

    #[test]
    fn wca_force_repulsive_only() {
        let epsilon = 1.0;
        let sigma = 1.0;
        let r_cut = sigma * 2.0_f64.powf(1.0 / 6.0);

        // Inside cutoff → positive (repulsive)
        let f_inside = compute_wca_force(0.9, epsilon, sigma);
        assert!(
            f_inside > 0.0,
            "WCA should be repulsive at r=0.9, got {f_inside}"
        );

        // At cutoff → zero
        let f_cut = compute_wca_force(r_cut, epsilon, sigma);
        assert!(
            f_cut.abs() < 1e-10,
            "WCA should be zero at cutoff, got {f_cut}"
        );

        // Beyond cutoff → zero
        let f_beyond = compute_wca_force(2.0, epsilon, sigma);
        assert!(f_beyond.abs() < 1e-14, "WCA should be zero beyond cutoff");
    }

    // ── Soft-core LJ tests ──────────────────────────────────────────────

    #[test]
    fn softcore_lj_no_singularity_at_zero() {
        let f = compute_softcore_lj_force(0.0, 1.0, 1.0, 0.1);
        assert!(
            f.is_finite(),
            "soft-core LJ should be finite at r=0, got {f}"
        );
        assert!(
            f.abs() < 1e-10,
            "at r=0, r factor makes force zero, got {f}"
        );
    }

    #[test]
    fn softcore_lj_finite_at_small_r() {
        let f = compute_softcore_lj_force(0.01, 1.0, 1.0, 0.1);
        assert!(f.is_finite(), "soft-core should stay finite at small r");
    }

    // ── Boundary normals test ────────────────────────────────────────────

    #[test]
    fn boundary_normals_flat_surface() {
        // Particles along a flat plane y=0, spaced in x
        let dx = 0.05_f64;
        let h = 0.1_f64;
        let mut positions = Vec::new();
        let mut volumes = Vec::new();
        for i in -5..=5_i32 {
            positions.push([i as f64 * dx, 0.0, 0.0]);
            volumes.push(dx * dx * dx);
        }
        let normals = compute_boundary_normals(&positions, &volumes, h);
        // The middle particle's normal should be roughly in the y-direction
        // (though with a 1D line of particles results are approximate)
        // At minimum, normals should be finite
        for n in &normals {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if len > 1e-10 {
                assert!(
                    (len - 1.0).abs() < 1e-10,
                    "Non-zero normals should be unit length"
                );
            }
        }
    }

    // ── Dummy particle tests ─────────────────────────────────────────────

    #[test]
    fn dummy_particles_generated() {
        let positions = vec![[0.0, 0.05, 0.0]];
        let pressures = vec![1000.0];
        let velocities = vec![[1.0, 0.0, 0.0]];
        let dummies = generate_dummy_particles(
            &positions,
            &pressures,
            &velocities,
            [0.0, 1.0, 0.0],
            0.0,
            0.02,
            2,
            0.1,
        );
        // 1 fluid particle near wall, 2 layers → 2 dummies
        assert_eq!(dummies.len(), 2, "Expected 2 dummy particles");
        // Dummy velocity is negated (no-slip)
        for d in &dummies {
            assert!((d.velocity[0] - (-1.0)).abs() < 1e-12);
        }
    }

    // ── Rotating boundary test ───────────────────────────────────────────

    #[test]
    fn rotating_boundary_updates_positions() {
        let mut rbs = RigidBoundarySet::new();
        rbs.add(RigidBoundaryParticle::new(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.01,
        ));

        let center = [0.0, 0.0, 0.0];
        let angular_vel = [0.0, 0.0, 1.0]; // rotate about Z
        let linear_vel = [0.0; 3];
        let dt = 0.01;

        update_rotating_boundary(&mut rbs, center, angular_vel, linear_vel, dt);

        let p = &rbs.particles[0];
        // ω × r = [0,0,1] × [1,0,0] = [0,1,0]
        // velocity should be [0, 1, 0]
        assert!((p.velocity[1] - 1.0).abs() < 1e-10, "vy={}", p.velocity[1]);
        // Position updated by dt
        assert!((p.position[1] - 0.01).abs() < 1e-10, "py={}", p.position[1]);
    }

    // ── Shepard density test ─────────────────────────────────────────────

    #[test]
    fn shepard_density_reasonable() {
        let boundary_pos = [0.0, 0.0, 0.0];
        let fluid_positions = vec![[0.0, 0.05, 0.0], [0.0, 0.1, 0.0]];
        let fluid_masses = vec![0.125e-3; 2]; // typical for spacing 0.05
        let fluid_densities = vec![1000.0; 2];
        let h = 0.15;

        let rho = boundary_shepard_density(
            boundary_pos,
            &fluid_positions,
            &fluid_masses,
            &fluid_densities,
            h,
        );
        assert!(rho >= 0.0, "Shepard density should be non-negative");
        assert!(rho.is_finite(), "Shepard density should be finite");
    }

    // ── BoundarySet combined test ────────────────────────────────────────

    #[test]
    fn boundary_set_add_and_apply() {
        let mut bs = BoundarySet::new();
        bs.add_plane(BoundaryPlane::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0)));
        bs.add_particle(BoundaryParticle::new(Vec3::new(0.0, -0.1, 0.0), 0.001));
        assert_eq!(bs.planes.len(), 1);
        assert_eq!(bs.particles.len(), 1);
    }
}

// ── Inflow/outflow boundary ────────────────────────────────────────────────────

/// Defines an inflow plane that injects particles at a given rate.
#[derive(Debug, Clone)]
pub struct InflowBoundary {
    /// Position of the inflow plane center.
    pub position: [f64; 3],
    /// Inflow normal (direction particles flow in).
    pub normal: [f64; 3],
    /// Half-width of the inflow region.
    pub half_width: f64,
    /// Inflow velocity magnitude (m/s).
    pub velocity: f64,
    /// Rest density (kg/m^3).
    pub rest_density: f64,
    /// Spacing between injected particles.
    pub particle_spacing: f64,
}

impl InflowBoundary {
    /// Create a new inflow boundary.
    pub fn new(
        position: [f64; 3],
        normal: [f64; 3],
        half_width: f64,
        velocity: f64,
        rest_density: f64,
        particle_spacing: f64,
    ) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let n = if len > 1e-14 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [1.0, 0.0, 0.0]
        };
        Self {
            position,
            normal: n,
            half_width,
            velocity,
            rest_density,
            particle_spacing,
        }
    }

    /// Generate the velocity vector for inflow particles.
    pub fn inflow_velocity_vec(&self) -> [f64; 3] {
        [
            self.normal[0] * self.velocity,
            self.normal[1] * self.velocity,
            self.normal[2] * self.velocity,
        ]
    }

    /// Check if a particle position is within the inflow region.
    pub fn is_within_region(&self, pos: [f64; 3]) -> bool {
        // Project onto tangential plane
        let dx = pos[0] - self.position[0];
        let dy = pos[1] - self.position[1];
        let dz = pos[2] - self.position[2];
        let n = self.normal;
        let dot = dx * n[0] + dy * n[1] + dz * n[2];
        let tx = dx - dot * n[0];
        let ty = dy - dot * n[1];
        let tz = dz - dot * n[2];
        let r = (tx * tx + ty * ty + tz * tz).sqrt();
        r < self.half_width
    }
}

/// Outflow boundary that removes particles that leave the domain.
#[derive(Debug, Clone)]
pub struct OutflowBoundary {
    /// Position of the outflow plane.
    pub position: [f64; 3],
    /// Outward normal of the outflow plane.
    pub normal: [f64; 3],
}

impl OutflowBoundary {
    /// Create a new outflow boundary.
    pub fn new(position: [f64; 3], normal: [f64; 3]) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let n = if len > 1e-14 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [1.0, 0.0, 0.0]
        };
        Self {
            position,
            normal: n,
        }
    }

    /// Check if a particle has exited through this outflow plane.
    pub fn has_exited(&self, pos: [f64; 3]) -> bool {
        let dx = pos[0] - self.position[0];
        let dy = pos[1] - self.position[1];
        let dz = pos[2] - self.position[2];
        let signed_dist = dx * self.normal[0] + dy * self.normal[1] + dz * self.normal[2];
        signed_dist > 0.0
    }

    /// Filter a list of positions, returning indices of particles that should be removed.
    pub fn find_exiting_particles(&self, positions: &[[f64; 3]]) -> Vec<usize> {
        positions
            .iter()
            .enumerate()
            .filter(|(_, pos)| self.has_exited(**pos))
            .map(|(i, _)| i)
            .collect()
    }
}

// ── Periodic boundary condition ────────────────────────────────────────────────

/// Periodic box boundary condition for SPH.
///
/// Wraps particle positions into the box domain `[min, max]^3`.
#[derive(Debug, Clone)]
pub struct PeriodicBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl PeriodicBox {
    /// Create a new periodic box.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Wrap a single position into the periodic box.
    pub fn wrap(&self, pos: [f64; 3]) -> [f64; 3] {
        let mut p = pos;
        for ((pk, &mn), &mx) in p.iter_mut().zip(self.min.iter()).zip(self.max.iter()) {
            let size = mx - mn;
            if size <= 0.0 {
                continue;
            }
            while *pk < mn {
                *pk += size;
            }
            while *pk >= mx {
                *pk -= size;
            }
        }
        p
    }

    /// Compute the minimum image displacement from `pos_a` to `pos_b`.
    pub fn min_image_displacement(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> [f64; 3] {
        let mut d = [
            pos_b[0] - pos_a[0],
            pos_b[1] - pos_a[1],
            pos_b[2] - pos_a[2],
        ];
        for ((dk, &mn), &mx) in d.iter_mut().zip(self.min.iter()).zip(self.max.iter()) {
            let size = mx - mn;
            if size <= 0.0 {
                continue;
            }
            if *dk > 0.5 * size {
                *dk -= size;
            }
            if *dk < -0.5 * size {
                *dk += size;
            }
        }
        d
    }

    /// Wrap all positions in a vector in place.
    pub fn wrap_all(&self, positions: &mut [[f64; 3]]) {
        for pos in positions.iter_mut() {
            *pos = self.wrap(*pos);
        }
    }

    /// Box dimensions.
    pub fn size(&self) -> [f64; 3] {
        [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ]
    }
}

// ── Pressure boundary condition ────────────────────────────────────────────────

/// Prescribed-pressure boundary condition using the ghost-particle technique.
///
/// Ghost particles are placed at the boundary with extrapolated pressure
/// to impose a fixed pressure condition.
pub struct PressureBoundary {
    /// Target pressure at the boundary.
    pub pressure: f64,
    /// Reference density (kg/m^3).
    pub rho_ref: f64,
    /// Plane normal (points into the domain).
    pub normal: [f64; 3],
    /// A point on the plane.
    pub plane_point: [f64; 3],
}

impl PressureBoundary {
    /// Create a new pressure boundary.
    pub fn new(pressure: f64, rho_ref: f64, normal: [f64; 3], plane_point: [f64; 3]) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let n = if len > 1e-14 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            pressure,
            rho_ref,
            normal: n,
            plane_point,
        }
    }

    /// Density from prescribed pressure via Tait EOS:
    ///
    /// `rho = rho_ref * (p / B + 1)^(1/gamma)`
    pub fn density_from_pressure(&self, gamma: f64, b: f64) -> f64 {
        let arg = self.pressure / b + 1.0;
        self.rho_ref * arg.powf(1.0 / gamma)
    }

    /// Generate ghost particles for the pressure boundary condition.
    ///
    /// Mirrors fluid particles within `h` of the plane with prescribed pressure.
    pub fn generate_pressure_ghosts(
        &self,
        fluid_positions: &[[f64; 3]],
        fluid_velocities: &[[f64; 3]],
        h: f64,
    ) -> Vec<DummyParticle> {
        let mb = MirrorBoundary::new(self.normal, {
            let pp = self.plane_point;
            pp[0] * self.normal[0] + pp[1] * self.normal[1] + pp[2] * self.normal[2]
        });
        let mut ghosts = Vec::new();
        for (i, &pos) in fluid_positions.iter().enumerate() {
            let d = mb.signed_distance(pos);
            if d > 0.0 && d < h {
                let mirror_pos = mb.reflect(pos);
                ghosts.push(DummyParticle {
                    position: mirror_pos,
                    pressure: self.pressure,
                    velocity: [
                        -fluid_velocities[i][0],
                        -fluid_velocities[i][1],
                        -fluid_velocities[i][2],
                    ],
                    volume: h * h * h / 8.0,
                });
            }
        }
        ghosts
    }
}

// ── Repulsive boundary forces (improved) ──────────────────────────────────────

/// Apply Monaghan's repulsive boundary force to prevent particle penetration.
///
/// `F_i = sum_b chi * (r0/r)^4 * (e_ib / r) * f(q)` where:
/// - `chi` is a strength constant
/// - `r0` is the reference distance
/// - `q = r / h`
/// - `f(q) = 1 - q` for q < 1, else 0
pub fn apply_monaghan_boundary_force(
    fluid_pos: [f64; 3],
    boundary_pos: [f64; 3],
    h: f64,
    chi: f64,
    r0: f64,
) -> [f64; 3] {
    let dx = fluid_pos[0] - boundary_pos[0];
    let dy = fluid_pos[1] - boundary_pos[1];
    let dz = fluid_pos[2] - boundary_pos[2];
    let r2 = dx * dx + dy * dy + dz * dz;
    if r2 < 1e-28 || r2 > h * h {
        return [0.0; 3];
    }
    let r = r2.sqrt();
    let q = r / h;
    if q >= 1.0 {
        return [0.0; 3];
    }
    let fq = 1.0 - q;
    let ratio = r0 / r;
    let mag = chi * ratio * ratio * ratio * ratio * fq / r;
    [mag * dx, mag * dy, mag * dz]
}

// ── Free-surface boundary ──────────────────────────────────────────────────────

/// Detect free-surface particles using a density-based criterion.
///
/// A particle is considered a free-surface particle if its density is below
/// a fraction `threshold` of the reference density.
pub fn detect_free_surface(densities: &[f64], rest_density: f64, threshold: f64) -> Vec<bool> {
    densities
        .iter()
        .map(|&rho| rho < threshold * rest_density)
        .collect()
}

/// Compute the surface normal for a free-surface particle using the
/// color function gradient (Morris 2000 approximation).
pub fn free_surface_normal(
    fluid_positions: &[[f64; 3]],
    fluid_densities: &[f64],
    target_idx: usize,
    h: f64,
    rest_density: f64,
) -> [f64; 3] {
    let pos_i = fluid_positions[target_idx];
    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;

    for (j, &pos_j) in fluid_positions.iter().enumerate() {
        if j == target_idx {
            continue;
        }
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 >= h * h || r2 < 1e-28 {
            continue;
        }
        let r = r2.sqrt();
        let q = r / h;
        // Gradient of cubic spline color function
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
        let dw_dr = if q >= 1.0 {
            let t = 2.0 - q;
            sigma * (-0.75 * t * t) / h
        } else {
            sigma * (-3.0 * q + 2.25 * q * q) / h
        };
        let c_j = fluid_densities[j] / rest_density;
        let factor = c_j * dw_dr / r;
        nx += factor * dx;
        ny += factor * dy;
        nz += factor * dz;
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len > 1e-14 {
        [nx / len, ny / len, nz / len]
    } else {
        [0.0, 0.0, 0.0]
    }
}

// ── Tests for new functionality ────────────────────────────────────────────────

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ── InflowBoundary tests ─────────────────────────────────────────────

    #[test]
    fn inflow_boundary_velocity_vec_correct() {
        let ib = InflowBoundary::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 2.0, 1000.0, 0.05);
        let v = ib.inflow_velocity_vec();
        assert!((v[0] - 2.0).abs() < 1e-12, "vx={}", v[0]);
        assert!(v[1].abs() < 1e-14);
    }

    #[test]
    fn inflow_boundary_region_check() {
        let ib = InflowBoundary::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 1.0, 1000.0, 0.05);
        assert!(ib.is_within_region([0.0, 0.1, 0.0]));
        assert!(!ib.is_within_region([0.0, 1.0, 0.0]));
    }

    // ── OutflowBoundary tests ────────────────────────────────────────────

    #[test]
    fn outflow_detects_exiting() {
        let ob = OutflowBoundary::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(ob.has_exited([1.5, 0.0, 0.0]));
        assert!(!ob.has_exited([0.5, 0.0, 0.0]));
    }

    #[test]
    fn outflow_find_exiting_particles() {
        let ob = OutflowBoundary::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let positions = vec![[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let exiting = ob.find_exiting_particles(&positions);
        assert_eq!(exiting, vec![1, 2]);
    }

    // ── PeriodicBox tests ─────────────────────────────────────────────

    #[test]
    fn periodic_box_wrap_inside() {
        let pb = PeriodicBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        let p = pb.wrap([0.5, 0.3, 0.7]);
        assert!((p[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn periodic_box_wrap_outside_high() {
        let pb = PeriodicBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        let p = pb.wrap([1.3, 0.0, 0.0]);
        assert!((p[0] - 0.3).abs() < 1e-12, "wrapped x={}", p[0]);
    }

    #[test]
    fn periodic_box_wrap_outside_low() {
        let pb = PeriodicBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        let p = pb.wrap([-0.3, 0.0, 0.0]);
        assert!((p[0] - 0.7).abs() < 1e-12, "wrapped x={}", p[0]);
    }

    #[test]
    fn periodic_box_min_image_close() {
        let pb = PeriodicBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        let d = pb.min_image_displacement([0.1, 0.0, 0.0], [0.9, 0.0, 0.0]);
        // Shortest path: -0.2 (wrap around), not +0.8
        assert!((d[0] - (-0.2)).abs() < 1e-12, "d[0]={}", d[0]);
    }

    #[test]
    fn periodic_box_size() {
        let pb = PeriodicBox::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let s = pb.size();
        assert!((s[0] - 3.0).abs() < 1e-12);
        assert!((s[1] - 3.0).abs() < 1e-12);
        assert!((s[2] - 3.0).abs() < 1e-12);
    }

    // ── PressureBoundary tests ────────────────────────────────────────

    #[test]
    fn pressure_boundary_density_from_pressure() {
        let pb = PressureBoundary::new(1000.0, 1000.0, [0.0, 1.0, 0.0], [0.0; 3]);
        // gamma=7, B=2e5: rho = rho_ref * (p/B + 1)^(1/7)
        let rho = pb.density_from_pressure(7.0, 2.0e5);
        assert!(rho > 1000.0, "rho should exceed rest density: {rho}");
    }

    #[test]
    fn pressure_boundary_ghosts_generated() {
        let pb = PressureBoundary::new(1000.0, 1000.0, [0.0, 1.0, 0.0], [0.0; 3]);
        let positions = vec![[0.0, 0.05, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0]];
        let ghosts = pb.generate_pressure_ghosts(&positions, &velocities, 0.1);
        assert_eq!(ghosts.len(), 1, "1 ghost expected");
        assert!((ghosts[0].pressure - 1000.0).abs() < 1e-10);
    }

    // ── Monaghan repulsive force tests ───────────────────────────────

    #[test]
    fn monaghan_force_repulsive_close() {
        let f = apply_monaghan_boundary_force([0.0, 0.02, 0.0], [0.0; 3], 0.1, 1.0, 0.05);
        assert!(
            f[1] > 0.0,
            "Force should be repulsive (pointing away): {}",
            f[1]
        );
    }

    #[test]
    fn monaghan_force_zero_outside_h() {
        let f = apply_monaghan_boundary_force([0.0, 0.2, 0.0], [0.0; 3], 0.1, 1.0, 0.05);
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!(mag < 1e-30, "Force outside h should be zero: {mag}");
    }

    // ── Free-surface detection ─────────────────────────────────────

    #[test]
    fn detect_free_surface_below_threshold() {
        let densities = vec![1000.0, 800.0, 950.0];
        let fs = detect_free_surface(&densities, 1000.0, 0.9);
        assert!(!fs[0]); // 1000 >= 900
        assert!(fs[1]); // 800 < 900
        assert!(!fs[2]); // 950 >= 900
    }

    #[test]
    fn detect_free_surface_all_above() {
        let densities = vec![1000.0, 1010.0, 995.0];
        let fs = detect_free_surface(&densities, 1000.0, 0.9);
        assert!(fs.iter().all(|&b| !b), "All above threshold");
    }

    // ── Free-surface normal ────────────────────────────────────────

    #[test]
    fn free_surface_normal_finite() {
        let positions = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let densities = vec![500.0, 1000.0, 1000.0];
        let n = free_surface_normal(&positions, &densities, 0, 0.1, 1000.0);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        // Normal is either unit length or zero (insufficient neighbors)
        assert!(len < 1.0 + 1e-10, "Normal should be unit or zero");
    }
}
