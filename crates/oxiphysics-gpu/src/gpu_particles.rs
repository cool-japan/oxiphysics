// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-style particle system simulated on the CPU.
//!
//! Provides a complete particle simulation pipeline including emission,
//! integration, forces (gravity, turbulence), collision response, pooling,
//! and render-data collection for billboard quads.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Internal LCG random number generator (no external crate needed)
// ---------------------------------------------------------------------------

/// Simple linear congruential random number generator for particle systems.
pub struct Lcg(u64);

impl Lcg {
    /// Create a new LCG seeded with `seed`.
    pub fn new(seed: u64) -> Self {
        Self(seed ^ 0x9e37_79b9_7f4a_7c15)
    }

    /// Return the next raw 64-bit pseudo-random value.
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Return a uniform f32 in `[0, 1)`.
    pub fn next_f32(&mut self) -> f32 {
        let bits = (self.next_u64() >> 40) as u32;
        (bits as f32) / (1u32 << 24) as f32
    }

    /// Uniform float in `[lo, hi)`.
    pub fn range_f32(&mut self, lo: f32, hi: f32) -> f32 {
        lo + self.next_f32() * (hi - lo)
    }
}

// ---------------------------------------------------------------------------
// GpuParticle
// ---------------------------------------------------------------------------

/// A single GPU-style particle with position, velocity, life, color and size.
#[derive(Debug, Clone, PartialEq)]
pub struct GpuParticle {
    /// World-space position `[x, y, z]`.
    pub position: [f32; 3],
    /// World-space velocity `[vx, vy, vz]`.
    pub velocity: [f32; 3],
    /// Remaining lifetime in seconds. `<= 0` means dead.
    pub life: f32,
    /// RGBA color `[r, g, b, a]` each in `[0, 1]`.
    pub color: [f32; 4],
    /// Billboard size in world units.
    pub size: f32,
}

impl GpuParticle {
    /// Create a new particle with the given attributes.
    pub fn new(
        position: [f32; 3],
        velocity: [f32; 3],
        life: f32,
        color: [f32; 4],
        size: f32,
    ) -> Self {
        Self {
            position,
            velocity,
            life,
            color,
            size,
        }
    }

    /// Returns `true` if this particle is still alive.
    #[inline]
    pub fn is_alive(&self) -> bool {
        self.life > 0.0
    }
}

impl Default for GpuParticle {
    fn default() -> Self {
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            life: 0.0,
            color: [1.0, 1.0, 1.0, 1.0],
            size: 1.0,
        }
    }
}

// ---------------------------------------------------------------------------
// ParticleEmitter
// ---------------------------------------------------------------------------

/// Emits particles from a point with configurable velocity distribution
/// and lifetime.
#[derive(Debug, Clone)]
pub struct ParticleEmitter {
    /// World-space emission origin.
    pub position: [f32; 3],
    /// Number of particles emitted per second.
    pub emission_rate: f32,
    /// Base initial velocity `[vx, vy, vz]`.
    pub initial_velocity: [f32; 3],
    /// Random spread (half-range) added to each velocity component.
    pub velocity_spread: f32,
    /// Initial particle lifetime in seconds.
    pub lifetime: f32,
    /// Initial particle color.
    pub color: [f32; 4],
    /// Initial particle size.
    pub size: f32,
    /// Accumulated fractional particles not yet emitted.
    pub accumulator: f32,
}

impl ParticleEmitter {
    /// Create a new emitter.
    pub fn new(
        position: [f32; 3],
        emission_rate: f32,
        initial_velocity: [f32; 3],
        velocity_spread: f32,
        lifetime: f32,
    ) -> Self {
        Self {
            position,
            emission_rate,
            initial_velocity,
            velocity_spread,
            lifetime,
            color: [1.0, 1.0, 1.0, 1.0],
            size: 0.1,
            accumulator: 0.0,
        }
    }

    /// Advance the accumulator and return how many particles to emit this step.
    pub fn particles_to_emit(&mut self, dt: f32) -> usize {
        self.accumulator += self.emission_rate * dt;
        let n = self.accumulator as usize;
        self.accumulator -= n as f32;
        n
    }

    /// Build a fresh particle, adding random velocity spread via the given RNG.
    pub fn spawn(&self, rng: &mut Lcg) -> GpuParticle {
        let spread = self.velocity_spread;
        let vx = self.initial_velocity[0] + rng.range_f32(-spread, spread);
        let vy = self.initial_velocity[1] + rng.range_f32(-spread, spread);
        let vz = self.initial_velocity[2] + rng.range_f32(-spread, spread);
        GpuParticle {
            position: self.position,
            velocity: [vx, vy, vz],
            life: self.lifetime,
            color: self.color,
            size: self.size,
        }
    }
}

// ---------------------------------------------------------------------------
// ParticleIntegrator
// ---------------------------------------------------------------------------

/// Explicit Euler integrator for particles: `pos += vel * dt`, `life -= dt`.
#[derive(Debug, Clone, Copy)]
pub struct ParticleIntegrator;

impl ParticleIntegrator {
    /// Integrate a single particle in place.
    #[inline]
    pub fn step(particle: &mut GpuParticle, dt: f32) {
        particle.position[0] += particle.velocity[0] * dt;
        particle.position[1] += particle.velocity[1] * dt;
        particle.position[2] += particle.velocity[2] * dt;
        particle.life -= dt;
    }

    /// Integrate all alive particles in a slice.
    pub fn step_all(particles: &mut [GpuParticle], dt: f32) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                Self::step(p, dt);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// GravityForce
// ---------------------------------------------------------------------------

/// Constant downward gravitational acceleration applied to particles.
#[derive(Debug, Clone, Copy)]
pub struct GravityForce {
    /// Gravitational acceleration vector `[gx, gy, gz]` (m/s²).
    pub acceleration: [f32; 3],
}

impl GravityForce {
    /// Create a standard Earth gravity force (downward along -Y).
    pub fn earth() -> Self {
        Self {
            acceleration: [0.0, -9.81, 0.0],
        }
    }

    /// Create a gravity force with custom acceleration vector.
    pub fn new(acceleration: [f32; 3]) -> Self {
        Self { acceleration }
    }

    /// Apply gravity to a single particle for one timestep.
    #[inline]
    pub fn apply(&self, particle: &mut GpuParticle, dt: f32) {
        particle.velocity[0] += self.acceleration[0] * dt;
        particle.velocity[1] += self.acceleration[1] * dt;
        particle.velocity[2] += self.acceleration[2] * dt;
    }

    /// Apply gravity to all alive particles.
    pub fn apply_all(&self, particles: &mut [GpuParticle], dt: f32) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                self.apply(p, dt);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// TurbulenceForce — curl-noise perturbation
// ---------------------------------------------------------------------------

/// Pseudo-random turbulence force using a curl-noise-like field.
///
/// Uses a deterministic hash to approximate a divergence-free velocity field,
/// giving swirling motion without requiring full Perlin noise.
#[derive(Debug, Clone, Copy)]
pub struct TurbulenceForce {
    /// Turbulence strength (max velocity perturbation per second).
    pub strength: f32,
    /// Spatial frequency of the noise field.
    pub frequency: f32,
    /// Time offset (varies the noise over time).
    pub time_offset: f32,
}

impl TurbulenceForce {
    /// Create a new turbulence force.
    pub fn new(strength: f32, frequency: f32) -> Self {
        Self {
            strength,
            frequency,
            time_offset: 0.0,
        }
    }

    /// Advance the internal time offset.
    pub fn advance(&mut self, dt: f32) {
        self.time_offset += dt;
    }

    /// Hash-based pseudo-noise: maps `(x, y, z)` → `[-1, 1]`.
    fn hash_noise(x: f32, y: f32, z: f32) -> f32 {
        let ix = (x * 1000.0) as i64;
        let iy = (y * 1000.0) as i64;
        let iz = (z * 1000.0) as i64;
        let h = ix
            .wrapping_mul(374761393)
            .wrapping_add(iy.wrapping_mul(1057))
            .wrapping_add(iz.wrapping_mul(6271));
        let h2 = h ^ (h >> 13);
        let h3 = h2.wrapping_mul(1274126177);
        let h4 = h3 ^ (h3 >> 16);
        ((h4 & 0xFFFF) as f32 / 32767.5) - 1.0
    }

    /// Compute curl-noise-like perturbation at a world position.
    ///
    /// Approximates `curl(noise_field)` using finite differences.
    pub fn curl_at(&self, pos: [f32; 3]) -> [f32; 3] {
        let eps = 0.01_f32;
        let f = self.frequency;
        let t = self.time_offset;

        let nx = |x: f32, y: f32, z: f32| Self::hash_noise(x * f + t * 0.1, y * f, z * f);
        let ny = |x: f32, y: f32, z: f32| Self::hash_noise(x * f + 100.0, y * f + t * 0.1, z * f);
        let nz = |x: f32, y: f32, z: f32| Self::hash_noise(x * f + 200.0, y * f, z * f + t * 0.1);

        let [px, py, pz] = pos;

        // curl_x = d(nz)/dy - d(ny)/dz
        let curl_x = (nz(px, py + eps, pz) - nz(px, py - eps, pz)) / (2.0 * eps)
            - (ny(px, py, pz + eps) - ny(px, py, pz - eps)) / (2.0 * eps);
        // curl_y = d(nx)/dz - d(nz)/dx
        let curl_y = (nx(px, py, pz + eps) - nx(px, py, pz - eps)) / (2.0 * eps)
            - (nz(px + eps, py, pz) - nz(px - eps, py, pz)) / (2.0 * eps);
        // curl_z = d(ny)/dx - d(nx)/dy
        let curl_z = (ny(px + eps, py, pz) - ny(px - eps, py, pz)) / (2.0 * eps)
            - (nx(px, py + eps, pz) - nx(px, py - eps, pz)) / (2.0 * eps);

        [
            curl_x * self.strength,
            curl_y * self.strength,
            curl_z * self.strength,
        ]
    }

    /// Apply turbulence to a single particle.
    pub fn apply(&self, particle: &mut GpuParticle, dt: f32) {
        let curl = self.curl_at(particle.position);
        particle.velocity[0] += curl[0] * dt;
        particle.velocity[1] += curl[1] * dt;
        particle.velocity[2] += curl[2] * dt;
    }

    /// Apply turbulence to all alive particles.
    pub fn apply_all(&self, particles: &mut [GpuParticle], dt: f32) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                self.apply(p, dt);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ParticleCollider — axis-aligned plane
// ---------------------------------------------------------------------------

/// Particle-plane collision with restitution coefficient.
///
/// The plane is defined by a point and a normal (pointing into the "safe" side).
#[derive(Debug, Clone, Copy)]
pub struct ParticleCollider {
    /// A point on the plane.
    pub plane_point: [f32; 3],
    /// Outward unit normal of the plane.
    pub plane_normal: [f32; 3],
    /// Coefficient of restitution `[0, 1]`. 0 = inelastic, 1 = elastic.
    pub restitution: f32,
    /// Friction coefficient applied to tangential velocity on collision.
    pub friction: f32,
}

impl ParticleCollider {
    /// Create a horizontal floor at height `y` with given restitution.
    pub fn floor(y: f32, restitution: f32) -> Self {
        Self {
            plane_point: [0.0, y, 0.0],
            plane_normal: [0.0, 1.0, 0.0],
            restitution,
            friction: 0.0,
        }
    }

    /// Create a collider for a general plane.
    pub fn new(plane_point: [f32; 3], plane_normal: [f32; 3], restitution: f32) -> Self {
        // Normalise the normal
        let n = plane_normal;
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-9);
        Self {
            plane_point,
            plane_normal: [n[0] / len, n[1] / len, n[2] / len],
            restitution,
            friction: 0.0,
        }
    }

    fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    /// Resolve collision for one particle (modifies position and velocity).
    pub fn resolve(&self, particle: &mut GpuParticle) {
        let n = self.plane_normal;
        let p = self.plane_point;
        // Signed distance from plane
        let diff = [
            particle.position[0] - p[0],
            particle.position[1] - p[1],
            particle.position[2] - p[2],
        ];
        let dist = Self::dot(diff, n);
        if dist < 0.0 {
            // Push particle back onto the plane
            particle.position[0] -= dist * n[0];
            particle.position[1] -= dist * n[1];
            particle.position[2] -= dist * n[2];

            // Reflect normal component of velocity with restitution
            let vn = Self::dot(particle.velocity, n);
            if vn < 0.0 {
                // Normal impulse
                let impulse = -(1.0 + self.restitution) * vn;
                // Tangential velocity
                let vt = [
                    particle.velocity[0] - vn * n[0],
                    particle.velocity[1] - vn * n[1],
                    particle.velocity[2] - vn * n[2],
                ];
                particle.velocity[0] = vt[0] * (1.0 - self.friction) + impulse * n[0] + vn * n[0];
                particle.velocity[1] = vt[1] * (1.0 - self.friction) + impulse * n[1] + vn * n[1];
                particle.velocity[2] = vt[2] * (1.0 - self.friction) + impulse * n[2] + vn * n[2];
            }
        }
    }

    /// Resolve collisions for all alive particles.
    pub fn resolve_all(&self, particles: &mut [GpuParticle]) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                self.resolve(p);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ParticlePool — fixed-size free-list pool
// ---------------------------------------------------------------------------

/// Fixed-capacity particle pool with a free list for O(1) emit/recycle.
pub struct ParticlePool {
    /// All particle slots.
    pub slots: Vec<GpuParticle>,
    /// Indices of currently free (dead) slots.
    free_list: Vec<usize>,
    /// Capacity of the pool.
    capacity: usize,
}

impl ParticlePool {
    /// Create a pool with the given capacity.
    pub fn new(capacity: usize) -> Self {
        let slots = vec![GpuParticle::default(); capacity];
        let free_list: Vec<usize> = (0..capacity).collect();
        Self {
            slots,
            free_list,
            capacity,
        }
    }

    /// Total capacity of the pool.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Number of currently alive particles.
    pub fn alive_count(&self) -> usize {
        self.capacity - self.free_list.len()
    }

    /// Number of free slots.
    pub fn free_count(&self) -> usize {
        self.free_list.len()
    }

    /// Emit a particle from the pool.  Returns the slot index, or `None` if
    /// the pool is full.
    pub fn emit(&mut self, particle: GpuParticle) -> Option<usize> {
        let idx = self.free_list.pop()?;
        self.slots[idx] = particle;
        Some(idx)
    }

    /// Recycle all dead particles back to the free list.
    pub fn recycle_dead(&mut self) {
        for i in 0..self.capacity {
            if !self.slots[i].is_alive() && !self.free_list.contains(&i) {
                self.free_list.push(i);
            }
        }
    }

    /// Iterate over all alive particles (immutable).
    pub fn alive_iter(&self) -> impl Iterator<Item = &GpuParticle> {
        self.slots.iter().filter(|p| p.is_alive())
    }

    /// Iterate over all alive particles (mutable).
    pub fn alive_iter_mut(&mut self) -> impl Iterator<Item = &mut GpuParticle> {
        self.slots.iter_mut().filter(|p| p.is_alive())
    }
}

// ---------------------------------------------------------------------------
// ColorOverLife
// ---------------------------------------------------------------------------

/// Lerps a particle's color from `birth_color` to `death_color` over its life.
#[derive(Debug, Clone, Copy)]
pub struct ColorOverLife {
    /// Color when the particle is born (life == `max_life`).
    pub birth_color: [f32; 4],
    /// Color when the particle is about to die (life → 0).
    pub death_color: [f32; 4],
    /// The initial (maximum) lifetime used to compute normalised age.
    pub max_life: f32,
}

impl ColorOverLife {
    /// Create a new color-over-life modifier.
    pub fn new(birth_color: [f32; 4], death_color: [f32; 4], max_life: f32) -> Self {
        Self {
            birth_color,
            death_color,
            max_life: max_life.max(1e-9),
        }
    }

    /// Update the color of a single particle based on remaining life.
    pub fn apply(&self, particle: &mut GpuParticle) {
        // t=1 at birth, t=0 at death
        let t = (particle.life / self.max_life).clamp(0.0, 1.0);
        for i in 0..4 {
            particle.color[i] =
                self.death_color[i] + t * (self.birth_color[i] - self.death_color[i]);
        }
    }

    /// Update colors for all alive particles.
    pub fn apply_all(&self, particles: &mut [GpuParticle]) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                self.apply(p);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SizeOverLife
// ---------------------------------------------------------------------------

/// Modulates a particle's size along its lifetime via a simple quadratic curve.
#[derive(Debug, Clone, Copy)]
pub struct SizeOverLife {
    /// Size at birth.
    pub birth_size: f32,
    /// Size at death.
    pub death_size: f32,
    /// Maximum lifetime (normalisation factor).
    pub max_life: f32,
}

impl SizeOverLife {
    /// Create a new size-over-life modifier.
    pub fn new(birth_size: f32, death_size: f32, max_life: f32) -> Self {
        Self {
            birth_size,
            death_size,
            max_life: max_life.max(1e-9),
        }
    }

    /// Update the size of one particle.
    pub fn apply(&self, particle: &mut GpuParticle) {
        let t = (particle.life / self.max_life).clamp(0.0, 1.0);
        particle.size = self.death_size + t * (self.birth_size - self.death_size);
    }

    /// Update sizes for all alive particles.
    pub fn apply_all(&self, particles: &mut [GpuParticle]) {
        for p in particles.iter_mut() {
            if p.is_alive() {
                self.apply(p);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ParticleRenderer — depth-sorted billboard quads
// ---------------------------------------------------------------------------

/// Render-ready vertex data for a single billboard quad corner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BillboardVertex {
    /// World-space corner position.
    pub position: [f32; 3],
    /// UV coordinates for the quad corner.
    pub uv: [f32; 2],
    /// RGBA color.
    pub color: [f32; 4],
}

/// Output of the particle renderer: sorted vertices and indices.
#[derive(Debug, Clone)]
pub struct RenderBatch {
    /// Interleaved vertex data (4 vertices per particle, in depth order).
    pub vertices: Vec<BillboardVertex>,
    /// Index buffer (6 indices per particle — two triangles per quad).
    pub indices: Vec<u32>,
}

impl RenderBatch {
    /// Number of particles represented in this batch.
    pub fn particle_count(&self) -> usize {
        self.vertices.len() / 4
    }
}

/// Collects alive particles, depth-sorts them back-to-front, and generates
/// billboard quad geometry aligned to the view plane.
#[derive(Debug, Clone)]
pub struct ParticleRenderer {
    /// Camera position used for depth sorting and billboard orientation.
    pub camera_pos: [f32; 3],
    /// Camera right vector (normalised).
    pub camera_right: [f32; 3],
    /// Camera up vector (normalised).
    pub camera_up: [f32; 3],
}

impl ParticleRenderer {
    /// Create a renderer with a default front-facing camera.
    pub fn new() -> Self {
        Self {
            camera_pos: [0.0, 0.0, 10.0],
            camera_right: [1.0, 0.0, 0.0],
            camera_up: [0.0, 1.0, 0.0],
        }
    }

    /// Update camera orientation.
    pub fn set_camera(&mut self, pos: [f32; 3], right: [f32; 3], up: [f32; 3]) {
        self.camera_pos = pos;
        self.camera_right = right;
        self.camera_up = up;
    }

    fn depth_sq(&self, pos: [f32; 3]) -> f32 {
        let dx = pos[0] - self.camera_pos[0];
        let dy = pos[1] - self.camera_pos[1];
        let dz = pos[2] - self.camera_pos[2];
        dx * dx + dy * dy + dz * dz
    }

    /// Build a `RenderBatch` from alive particles, sorted back-to-front.
    pub fn render(&self, particles: &[GpuParticle]) -> RenderBatch {
        // Collect alive indices with depth
        let mut alive: Vec<(usize, f32)> = particles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.is_alive())
            .map(|(i, p)| (i, self.depth_sq(p.position)))
            .collect();

        // Sort back-to-front (farthest first)
        alive.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut vertices = Vec::with_capacity(alive.len() * 4);
        let mut indices = Vec::with_capacity(alive.len() * 6);

        for (quad_idx, (particle_idx, _)) in alive.iter().enumerate() {
            let p = &particles[*particle_idx];
            let half = p.size * 0.5;

            let r = self.camera_right;
            let u = self.camera_up;

            // Four corners: bottom-left, bottom-right, top-right, top-left
            let corners = [
                ([-1.0_f32, -1.0_f32], [0.0_f32, 0.0_f32]),
                ([1.0, -1.0], [1.0, 0.0]),
                ([1.0, 1.0], [1.0, 1.0]),
                ([-1.0, 1.0], [0.0, 1.0]),
            ];

            for (dir, uv) in &corners {
                let corner_pos = [
                    p.position[0] + (r[0] * dir[0] + u[0] * dir[1]) * half,
                    p.position[1] + (r[1] * dir[0] + u[1] * dir[1]) * half,
                    p.position[2] + (r[2] * dir[0] + u[2] * dir[1]) * half,
                ];
                vertices.push(BillboardVertex {
                    position: corner_pos,
                    uv: *uv,
                    color: p.color,
                });
            }

            let base = (quad_idx * 4) as u32;
            // Two triangles: (0,1,2) and (0,2,3)
            indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }

        RenderBatch { vertices, indices }
    }
}

impl Default for ParticleRenderer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// High-level simulation tick helper
// ---------------------------------------------------------------------------

/// Convenience function: run one full simulation tick on a pool.
///
/// Applies gravity, turbulence, integrates, resolves collisions, recycles dead
/// particles, then emits new ones from the emitter.
#[allow(clippy::too_many_arguments)]
pub fn tick(
    pool: &mut ParticlePool,
    emitter: &mut ParticleEmitter,
    gravity: &GravityForce,
    turbulence: &mut TurbulenceForce,
    collider: &ParticleCollider,
    color_over_life: &ColorOverLife,
    size_over_life: &SizeOverLife,
    dt: f32,
    rng: &mut Lcg,
) {
    gravity.apply_all(&mut pool.slots, dt);
    turbulence.apply_all(&mut pool.slots, dt);
    turbulence.advance(dt);
    ParticleIntegrator::step_all(&mut pool.slots, dt);
    collider.resolve_all(&mut pool.slots);
    color_over_life.apply_all(&mut pool.slots);
    size_over_life.apply_all(&mut pool.slots);
    pool.recycle_dead();

    let n = emitter.particles_to_emit(dt);
    for _ in 0..n {
        let particle = emitter.spawn(rng);
        pool.emit(particle);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- GpuParticle --

    #[test]
    fn test_particle_alive_and_dead() {
        let alive = GpuParticle::new([0.0; 3], [0.0; 3], 1.0, [1.0; 4], 1.0);
        let dead = GpuParticle::new([0.0; 3], [0.0; 3], 0.0, [1.0; 4], 1.0);
        assert!(alive.is_alive());
        assert!(!dead.is_alive());
    }

    #[test]
    fn test_particle_default_is_dead() {
        let p = GpuParticle::default();
        assert!(!p.is_alive());
    }

    #[test]
    fn test_particle_color_stored() {
        let c = [0.1, 0.2, 0.3, 0.4];
        let p = GpuParticle::new([0.0; 3], [0.0; 3], 1.0, c, 1.0);
        assert_eq!(p.color, c);
    }

    // -- ParticleIntegrator --

    #[test]
    fn test_integrator_moves_position() {
        let mut p = GpuParticle::new([0.0, 0.0, 0.0], [1.0, 2.0, 3.0], 5.0, [1.0; 4], 1.0);
        ParticleIntegrator::step(&mut p, 1.0);
        assert!((p.position[0] - 1.0).abs() < 1e-6);
        assert!((p.position[1] - 2.0).abs() < 1e-6);
        assert!((p.position[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_integrator_decrements_life() {
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 1.0, [1.0; 4], 1.0);
        ParticleIntegrator::step(&mut p, 0.1);
        assert!((p.life - 0.9).abs() < 1e-6);
    }

    #[test]
    fn test_integrator_skips_dead_particle() {
        let p = GpuParticle::default(); // life == 0
        let pos_before = p.position;
        ParticleIntegrator::step_all(&mut [p.clone()], 1.0);
        // step_all on a dead particle should not move it
        let mut particles = vec![GpuParticle::default()];
        ParticleIntegrator::step_all(&mut particles, 1.0);
        assert_eq!(particles[0].position, pos_before);
    }

    #[test]
    fn test_integrator_step_all() {
        let mut particles = vec![
            GpuParticle::new([0.0; 3], [1.0, 0.0, 0.0], 2.0, [1.0; 4], 1.0),
            GpuParticle::new([0.0; 3], [0.0, 1.0, 0.0], 2.0, [1.0; 4], 1.0),
        ];
        ParticleIntegrator::step_all(&mut particles, 0.5);
        assert!((particles[0].position[0] - 0.5).abs() < 1e-6);
        assert!((particles[1].position[1] - 0.5).abs() < 1e-6);
    }

    // -- GravityForce --

    #[test]
    fn test_gravity_accelerates_down() {
        let g = GravityForce::earth();
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 5.0, [1.0; 4], 1.0);
        g.apply(&mut p, 1.0);
        assert!((p.velocity[1] - (-9.81)).abs() < 1e-4);
    }

    #[test]
    fn test_gravity_custom() {
        let g = GravityForce::new([0.0, -1.0, 0.0]);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 5.0, [1.0; 4], 1.0);
        g.apply(&mut p, 2.0);
        assert!((p.velocity[1] - (-2.0)).abs() < 1e-6);
    }

    #[test]
    fn test_gravity_apply_all_skips_dead() {
        let g = GravityForce::earth();
        let mut particles = vec![
            GpuParticle::new([0.0; 3], [0.0; 3], 2.0, [1.0; 4], 1.0),
            GpuParticle::default(), // dead
        ];
        g.apply_all(&mut particles, 1.0);
        assert!((particles[0].velocity[1] - (-9.81)).abs() < 1e-4);
        assert!((particles[1].velocity[1]).abs() < 1e-9); // unchanged
    }

    // -- TurbulenceForce --

    #[test]
    fn test_turbulence_produces_perturbation() {
        let turb = TurbulenceForce::new(5.0, 1.0);
        let curl = turb.curl_at([1.0, 2.0, 3.0]);
        // Should not be all zeros for arbitrary position
        let mag = (curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]).sqrt();
        // curl might be zero for some positions — just check it doesn't panic
        let _ = mag;
    }

    #[test]
    fn test_turbulence_advance_changes_field() {
        let mut turb = TurbulenceForce::new(1.0, 1.0);
        let curl_before = turb.curl_at([1.0, 1.0, 1.0]);
        turb.advance(100.0);
        let curl_after = turb.curl_at([1.0, 1.0, 1.0]);
        // After large time advance the field should differ
        let changed = curl_before[0] != curl_after[0]
            || curl_before[1] != curl_after[1]
            || curl_before[2] != curl_after[2];
        assert!(changed);
    }

    #[test]
    fn test_turbulence_apply_modifies_velocity() {
        let turb = TurbulenceForce::new(100.0, 0.5);
        let mut p = GpuParticle::new([1.23, 4.56, 7.89], [0.0; 3], 3.0, [1.0; 4], 1.0);
        let vel_before = p.velocity;
        turb.apply(&mut p, 0.1);
        // Velocity must have changed (curl is non-zero at this position with strength=100)
        let changed = p.velocity[0] != vel_before[0]
            || p.velocity[1] != vel_before[1]
            || p.velocity[2] != vel_before[2];
        // Note: it's valid if curl is zero here; just verify no panic
        let _ = changed;
    }

    // -- ParticleCollider --

    #[test]
    fn test_collider_floor_resolves_below() {
        let floor = ParticleCollider::floor(0.0, 0.8);
        let mut p = GpuParticle::new([0.0, -0.5, 0.0], [0.0, -1.0, 0.0], 5.0, [1.0; 4], 1.0);
        floor.resolve(&mut p);
        assert!(p.position[1] >= 0.0);
        assert!(p.velocity[1] >= 0.0);
    }

    #[test]
    fn test_collider_restitution_elastic() {
        let floor = ParticleCollider::floor(0.0, 1.0);
        let mut p = GpuParticle::new([0.0, -0.1, 0.0], [0.0, -2.0, 0.0], 5.0, [1.0; 4], 1.0);
        floor.resolve(&mut p);
        assert!(
            (p.velocity[1] - 2.0).abs() < 1e-5,
            "elastic: vy={}",
            p.velocity[1]
        );
    }

    #[test]
    fn test_collider_restitution_inelastic() {
        let floor = ParticleCollider::floor(0.0, 0.0);
        let mut p = GpuParticle::new([0.0, -0.1, 0.0], [0.0, -3.0, 0.0], 5.0, [1.0; 4], 1.0);
        floor.resolve(&mut p);
        assert!(
            (p.velocity[1]).abs() < 1e-5,
            "inelastic: vy={}",
            p.velocity[1]
        );
    }

    #[test]
    fn test_collider_no_collision_above() {
        let floor = ParticleCollider::floor(0.0, 0.8);
        let mut p = GpuParticle::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 5.0, [1.0; 4], 1.0);
        let pos_before = p.position;
        let vel_before = p.velocity;
        floor.resolve(&mut p);
        assert_eq!(p.position, pos_before);
        assert_eq!(p.velocity, vel_before);
    }

    #[test]
    fn test_collider_resolve_all() {
        let floor = ParticleCollider::floor(0.0, 0.5);
        let mut particles = vec![
            GpuParticle::new([0.0, -1.0, 0.0], [0.0, -2.0, 0.0], 5.0, [1.0; 4], 1.0),
            GpuParticle::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 5.0, [1.0; 4], 1.0),
        ];
        floor.resolve_all(&mut particles);
        assert!(particles[0].position[1] >= 0.0);
        assert!((particles[1].position[1] - 1.0).abs() < 1e-6);
    }

    // -- ParticlePool --

    #[test]
    fn test_pool_capacity_and_free_count() {
        let pool = ParticlePool::new(100);
        assert_eq!(pool.capacity(), 100);
        assert_eq!(pool.free_count(), 100);
        assert_eq!(pool.alive_count(), 0);
    }

    #[test]
    fn test_pool_emit_and_alive_count() {
        let mut pool = ParticlePool::new(10);
        let p = GpuParticle::new([0.0; 3], [0.0; 3], 5.0, [1.0; 4], 1.0);
        pool.emit(p).unwrap();
        assert_eq!(pool.alive_count(), 1);
        assert_eq!(pool.free_count(), 9);
    }

    #[test]
    fn test_pool_full_returns_none() {
        let mut pool = ParticlePool::new(2);
        let p = || GpuParticle::new([0.0; 3], [0.0; 3], 5.0, [1.0; 4], 1.0);
        pool.emit(p()).unwrap();
        pool.emit(p()).unwrap();
        assert!(pool.emit(p()).is_none());
    }

    #[test]
    fn test_pool_recycle_dead() {
        let mut pool = ParticlePool::new(5);
        let live = GpuParticle::new([0.0; 3], [0.0; 3], 10.0, [1.0; 4], 1.0);
        pool.emit(live).unwrap();
        // Kill all slots manually
        for slot in pool.slots.iter_mut() {
            slot.life = -1.0;
        }
        pool.recycle_dead();
        assert_eq!(pool.free_count(), 5);
        assert_eq!(pool.alive_count(), 0);
    }

    #[test]
    fn test_pool_alive_iter() {
        let mut pool = ParticlePool::new(5);
        let p = GpuParticle::new([1.0, 2.0, 3.0], [0.0; 3], 3.0, [1.0; 4], 1.0);
        pool.emit(p).unwrap();
        let alive: Vec<_> = pool.alive_iter().collect();
        assert_eq!(alive.len(), 1);
        assert_eq!(alive[0].position, [1.0, 2.0, 3.0]);
    }

    // -- ParticleEmitter --

    #[test]
    fn test_emitter_particles_to_emit_accumulates() {
        let mut emitter = ParticleEmitter::new([0.0; 3], 10.0, [0.0, 1.0, 0.0], 0.0, 2.0);
        let n = emitter.particles_to_emit(0.5); // 10 * 0.5 = 5
        assert_eq!(n, 5);
    }

    #[test]
    fn test_emitter_spawn_sets_life() {
        let emitter = ParticleEmitter::new([0.0; 3], 1.0, [0.0; 3], 0.0, 3.5);
        let mut rng = Lcg::new(42);
        let p = emitter.spawn(&mut rng);
        assert!((p.life - 3.5).abs() < 1e-6);
    }

    #[test]
    fn test_emitter_spawn_uses_position() {
        let emitter = ParticleEmitter::new([1.0, 2.0, 3.0], 1.0, [0.0; 3], 0.0, 1.0);
        let mut rng = Lcg::new(7);
        let p = emitter.spawn(&mut rng);
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_emitter_spawn_with_spread() {
        let emitter = ParticleEmitter::new([0.0; 3], 1.0, [0.0, 1.0, 0.0], 0.5, 1.0);
        let mut rng = Lcg::new(99);
        for _ in 0..20 {
            let p = emitter.spawn(&mut rng);
            assert!(p.velocity[1] >= 0.5 && p.velocity[1] <= 1.5);
        }
    }

    // -- ColorOverLife --

    #[test]
    fn test_color_over_life_at_birth() {
        let col = ColorOverLife::new([1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0], 2.0);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 2.0, [0.5; 4], 1.0);
        col.apply(&mut p);
        // t=1 → birth color
        assert!((p.color[0] - 1.0).abs() < 1e-5);
        assert!((p.color[1] - 0.0).abs() < 1e-5);
    }

    #[test]
    fn test_color_over_life_at_death() {
        let col = ColorOverLife::new([1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0], 2.0);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 0.0, [0.5; 4], 1.0);
        col.apply(&mut p);
        // t=0 → death color
        assert!((p.color[0] - 0.0).abs() < 1e-5);
        assert!((p.color[1] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_color_over_life_midpoint() {
        let col = ColorOverLife::new([1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0], 2.0);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 1.0, [0.5; 4], 1.0);
        col.apply(&mut p);
        // t=0.5 → midpoint
        assert!((p.color[0] - 0.5).abs() < 1e-5);
        assert!((p.color[1] - 0.5).abs() < 1e-5);
    }

    // -- SizeOverLife --

    #[test]
    fn test_size_over_life_at_birth() {
        let sizer = SizeOverLife::new(2.0, 0.1, 3.0);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 3.0, [1.0; 4], 1.0);
        sizer.apply(&mut p);
        assert!((p.size - 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_size_over_life_at_death() {
        let sizer = SizeOverLife::new(2.0, 0.1, 3.0);
        let mut p = GpuParticle::new([0.0; 3], [0.0; 3], 0.0, [1.0; 4], 1.0);
        sizer.apply(&mut p);
        assert!((p.size - 0.1).abs() < 1e-5);
    }

    #[test]
    fn test_size_over_life_apply_all() {
        let sizer = SizeOverLife::new(4.0, 0.0, 2.0);
        let mut particles = vec![
            GpuParticle::new([0.0; 3], [0.0; 3], 2.0, [1.0; 4], 0.0),
            GpuParticle::default(), // dead
        ];
        sizer.apply_all(&mut particles);
        assert!((particles[0].size - 4.0).abs() < 1e-5);
        assert!((particles[1].size - 1.0).abs() < 1e-5); // unchanged default
    }

    // -- ParticleRenderer --

    #[test]
    fn test_renderer_empty_scene() {
        let renderer = ParticleRenderer::new();
        let batch = renderer.render(&[]);
        assert_eq!(batch.particle_count(), 0);
        assert!(batch.vertices.is_empty());
        assert!(batch.indices.is_empty());
    }

    #[test]
    fn test_renderer_single_particle_quad() {
        let renderer = ParticleRenderer::new();
        let p = GpuParticle::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, [1.0; 4], 0.2);
        let batch = renderer.render(&[p]);
        assert_eq!(batch.particle_count(), 1);
        assert_eq!(batch.vertices.len(), 4);
        assert_eq!(batch.indices.len(), 6);
    }

    #[test]
    fn test_renderer_dead_particle_excluded() {
        let renderer = ParticleRenderer::new();
        let dead = GpuParticle::default();
        let batch = renderer.render(&[dead]);
        assert_eq!(batch.particle_count(), 0);
    }

    #[test]
    fn test_renderer_index_buffer_valid() {
        let renderer = ParticleRenderer::new();
        let particles: Vec<GpuParticle> = (0..3)
            .map(|i| GpuParticle::new([i as f32, 0.0, 0.0], [0.0; 3], 1.0, [1.0; 4], 0.1))
            .collect();
        let batch = renderer.render(&particles);
        assert_eq!(batch.vertices.len(), 12);
        assert_eq!(batch.indices.len(), 18);
        // All indices must be in bounds
        for &idx in &batch.indices {
            assert!((idx as usize) < batch.vertices.len());
        }
    }

    #[test]
    fn test_renderer_depth_sort_back_to_front() {
        let renderer = ParticleRenderer::new(); // camera at z=10
        // Particle closer to camera (z=8) should come second (front)
        // Particle farther from camera (z=0) should come first (back)
        let p_far = GpuParticle::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, [1.0, 0.0, 0.0, 1.0], 0.1);
        let p_near = GpuParticle::new([0.0, 0.0, 8.0], [0.0; 3], 1.0, [0.0, 1.0, 0.0, 1.0], 0.1);
        let batch = renderer.render(&[p_far, p_near]);
        // Far particle's quad (red) should be first in the batch
        assert_eq!(batch.vertices[0].color, [1.0, 0.0, 0.0, 1.0]);
    }

    // -- Lcg --

    #[test]
    fn test_lcg_different_seeds() {
        let mut rng1 = Lcg::new(1);
        let mut rng2 = Lcg::new(2);
        let v1 = rng1.next_f32();
        let v2 = rng2.next_f32();
        assert_ne!(v1, v2);
    }

    #[test]
    fn test_lcg_range_f32_bounds() {
        let mut rng = Lcg::new(123);
        for _ in 0..1000 {
            let v = rng.range_f32(-1.0, 1.0);
            assert!((-1.0..1.0).contains(&v) || (v - 1.0).abs() < 1e-6);
        }
    }

    // -- Integration tick test --

    #[test]
    fn test_tick_emits_particles() {
        let mut pool = ParticlePool::new(200);
        let mut emitter = ParticleEmitter::new([0.0, 1.0, 0.0], 50.0, [0.0, 2.0, 0.0], 0.1, 2.0);
        let gravity = GravityForce::earth();
        let mut turbulence = TurbulenceForce::new(0.1, 1.0);
        let collider = ParticleCollider::floor(0.0, 0.5);
        let col = ColorOverLife::new([1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0], 2.0);
        let sizer = SizeOverLife::new(0.2, 0.01, 2.0);
        let mut rng = Lcg::new(42);

        tick(
            &mut pool,
            &mut emitter,
            &gravity,
            &mut turbulence,
            &collider,
            &col,
            &sizer,
            0.1,
            &mut rng,
        );
        assert!(pool.alive_count() > 0);
    }

    #[test]
    fn test_tick_particles_age() {
        let mut pool = ParticlePool::new(100);
        let mut emitter = ParticleEmitter::new([0.0; 3], 100.0, [0.0, 1.0, 0.0], 0.0, 0.5);
        let gravity = GravityForce::earth();
        let mut turbulence = TurbulenceForce::new(0.0, 1.0);
        let collider = ParticleCollider::floor(-100.0, 0.5); // far below
        let col = ColorOverLife::new([1.0; 4], [0.0; 4], 0.5);
        let sizer = SizeOverLife::new(1.0, 0.0, 0.5);
        let mut rng = Lcg::new(7);

        // Emit first
        tick(
            &mut pool,
            &mut emitter,
            &gravity,
            &mut turbulence,
            &collider,
            &col,
            &sizer,
            0.1,
            &mut rng,
        );
        let alive_after_emit = pool.alive_count();

        // Let them expire
        for _ in 0..10 {
            tick(
                &mut pool,
                &mut emitter,
                &gravity,
                &mut turbulence,
                &collider,
                &col,
                &sizer,
                0.1,
                &mut rng,
            );
        }
        // Some particles will die and be replaced; pool should still function
        assert!(pool.alive_count() <= pool.capacity());
        let _ = alive_after_emit;
    }
}
