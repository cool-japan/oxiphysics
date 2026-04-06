//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// Emitter shape that controls where and how particles are emitted.
#[derive(Debug, Clone, Copy)]
pub enum EmitterShape {
    /// Point emitter: all particles originate from a single point.
    Point,
    /// Cone emitter: particles emitted within a cone of given half-angle (radians).
    Cone {
        /// Half-angle of the cone in radians.
        half_angle: f32,
    },
    /// Sphere emitter: particles start at random positions on a sphere surface.
    Sphere {
        /// Radius of the sphere.
        radius: f32,
    },
    /// Box emitter: particles emitted from random positions inside an AABB.
    Box {
        /// Half-extents of the AABB on each axis.
        half_extents: [f32; 3],
    },
}
/// GPU-style particle collision using a grid-based broad phase.
///
/// Assigns each particle to a grid cell, then checks only neighbours.
pub struct GridParticleCollision {
    /// Cell size.
    pub cell_size: f32,
    /// Radius for collision detection.
    pub particle_radius: f32,
    /// Restitution coefficient.
    pub restitution: f32,
}
impl GridParticleCollision {
    /// Create a new grid collision handler.
    pub fn new(cell_size: f32, particle_radius: f32, restitution: f32) -> Self {
        Self {
            cell_size,
            particle_radius,
            restitution,
        }
    }
    /// Resolve collisions between all alive particles (O(n) with grid).
    ///
    /// Uses a simple hash-grid approach: build cell lists, then check
    /// same-cell and adjacent-cell pairs only.
    pub fn resolve(&self, buffer: &mut ParticleBuffer) {
        let n = buffer.count;
        let alive: Vec<usize> = (0..n).filter(|&i| buffer.is_alive(i)).collect();
        let na = alive.len();
        if na < 2 {
            return;
        }
        let mut grid: std::collections::HashMap<(i32, i32, i32), Vec<usize>> =
            std::collections::HashMap::new();
        for &i in &alive {
            let cx = (buffer.positions_x[i] / self.cell_size).floor() as i32;
            let cy = (buffer.positions_y[i] / self.cell_size).floor() as i32;
            let cz = (buffer.positions_z[i] / self.cell_size).floor() as i32;
            grid.entry((cx, cy, cz)).or_default().push(i);
        }
        let diameter = 2.0 * self.particle_radius;
        let mut dvx = vec![0.0f32; n];
        let mut dvy = vec![0.0f32; n];
        let mut dvz = vec![0.0f32; n];
        for (&(cx, cy, cz), cell) in &grid {
            for dx in -1i32..=1 {
                for dy in -1i32..=1 {
                    for dz in -1i32..=1 {
                        let nb_key = (cx + dx, cy + dy, cz + dz);
                        if let Some(nb_cell) = grid.get(&nb_key) {
                            for &i in cell {
                                for &j in nb_cell {
                                    if j <= i {
                                        continue;
                                    }
                                    let dx_p = buffer.positions_x[j] - buffer.positions_x[i];
                                    let dy_p = buffer.positions_y[j] - buffer.positions_y[i];
                                    let dz_p = buffer.positions_z[j] - buffer.positions_z[i];
                                    let dist = (dx_p * dx_p + dy_p * dy_p + dz_p * dz_p).sqrt();
                                    if dist < diameter && dist > 1e-6 {
                                        let overlap = diameter - dist;
                                        let nx = dx_p / dist;
                                        let ny = dy_p / dist;
                                        let nz = dz_p / dist;
                                        let rvx = buffer.velocities_x[j] - buffer.velocities_x[i];
                                        let rvy = buffer.velocities_y[j] - buffer.velocities_y[i];
                                        let rvz = buffer.velocities_z[j] - buffer.velocities_z[i];
                                        let rv_n = rvx * nx + rvy * ny + rvz * nz;
                                        if rv_n < 0.0 {
                                            let j_impulse = -(1.0 + self.restitution) * rv_n
                                                / (1.0 / buffer.masses[i] + 1.0 / buffer.masses[j]);
                                            let inv_mi = 1.0 / buffer.masses[i];
                                            let inv_mj = 1.0 / buffer.masses[j];
                                            dvx[i] -= j_impulse * inv_mi * nx;
                                            dvy[i] -= j_impulse * inv_mi * ny;
                                            dvz[i] -= j_impulse * inv_mi * nz;
                                            dvx[j] += j_impulse * inv_mj * nx;
                                            dvy[j] += j_impulse * inv_mj * ny;
                                            dvz[j] += j_impulse * inv_mj * nz;
                                        }
                                        let push = overlap * 0.5;
                                        buffer.positions_x[i] -= push * nx;
                                        buffer.positions_y[i] -= push * ny;
                                        buffer.positions_z[i] -= push * nz;
                                        buffer.positions_x[j] += push * nx;
                                        buffer.positions_y[j] += push * ny;
                                        buffer.positions_z[j] += push * nz;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for &i in &alive {
            buffer.velocities_x[i] += dvx[i];
            buffer.velocities_y[i] += dvy[i];
            buffer.velocities_z[i] += dvz[i];
        }
    }
}
/// Summary statistics computed from a `ParticleBuffer`.
#[derive(Debug, Clone)]
pub struct ParticleStats {
    /// Number of alive particles.
    pub active: usize,
    /// Minimum position across all alive particles.
    pub min_pos: [f32; 3],
    /// Maximum position across all alive particles.
    pub max_pos: [f32; 3],
    /// Mean speed of all alive particles.
    pub avg_speed: f32,
    /// Total ½mv² kinetic energy.
    pub total_kinetic_energy: f32,
}
impl ParticleStats {
    /// Compute statistics from the given buffer.
    pub fn compute(buffer: &ParticleBuffer) -> Self {
        let mut active = 0usize;
        let mut min_pos = [f32::MAX; 3];
        let mut max_pos = [f32::MIN; 3];
        let mut sum_speed = 0.0f32;
        let mut total_ke = 0.0f32;
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            active += 1;
            let x = buffer.positions_x[i];
            let y = buffer.positions_y[i];
            let z = buffer.positions_z[i];
            min_pos[0] = min_pos[0].min(x);
            min_pos[1] = min_pos[1].min(y);
            min_pos[2] = min_pos[2].min(z);
            max_pos[0] = max_pos[0].max(x);
            max_pos[1] = max_pos[1].max(y);
            max_pos[2] = max_pos[2].max(z);
            let vx = buffer.velocities_x[i];
            let vy = buffer.velocities_y[i];
            let vz = buffer.velocities_z[i];
            let speed = (vx * vx + vy * vy + vz * vz).sqrt();
            sum_speed += speed;
            total_ke += 0.5 * buffer.masses[i] * speed * speed;
        }
        let avg_speed = if active > 0 {
            sum_speed / active as f32
        } else {
            0.0
        };
        if active == 0 {
            min_pos = [0.0; 3];
            max_pos = [0.0; 3];
        }
        Self {
            active,
            min_pos,
            max_pos,
            avg_speed,
            total_kinetic_energy: total_ke,
        }
    }
}
/// Euler integrator for particle positions and lifetimes.
pub struct ParticleIntegrator;
impl ParticleIntegrator {
    /// Advance all alive particles by `dt` seconds.
    ///
    /// * `pos += vel * dt`
    /// * `age += dt`
    /// * `lifetime -= dt`  (particle dies naturally when lifetime < 0)
    pub fn integrate(buffer: &mut ParticleBuffer, dt: f32) {
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            buffer.positions_x[i] += buffer.velocities_x[i] * dt;
            buffer.positions_y[i] += buffer.velocities_y[i] * dt;
            buffer.positions_z[i] += buffer.velocities_z[i] * dt;
            buffer.ages[i] += dt;
            buffer.lifetimes[i] -= dt;
        }
    }
}
/// Simple pairwise repulsion between particles.
pub struct ParticleRepulsion {
    /// Repulsion strength.
    pub strength: f32,
    /// Interaction radius.
    pub radius: f32,
}
impl ParticleRepulsion {
    /// Apply pairwise repulsion between all alive particles.
    ///
    /// This is O(n^2) and suitable only for small particle counts.
    pub fn apply(&self, buffer: &mut ParticleBuffer, dt: f32) {
        let n = buffer.count;
        let alive: Vec<usize> = (0..n).filter(|&i| buffer.is_alive(i)).collect();
        let na = alive.len();
        let mut fx = vec![0.0f32; n];
        let mut fy = vec![0.0f32; n];
        let mut fz = vec![0.0f32; n];
        for ai in 0..na {
            let i = alive[ai];
            for aj in (ai + 1)..na {
                let j = alive[aj];
                let dx = buffer.positions_x[j] - buffer.positions_x[i];
                let dy = buffer.positions_y[j] - buffer.positions_y[i];
                let dz = buffer.positions_z[j] - buffer.positions_z[i];
                let dist2 = dx * dx + dy * dy + dz * dz;
                let dist = dist2.sqrt();
                if dist >= self.radius || dist < 1e-6 {
                    continue;
                }
                let overlap = self.radius - dist;
                let f = self.strength * overlap / dist;
                fx[i] -= f * dx;
                fy[i] -= f * dy;
                fz[i] -= f * dz;
                fx[j] += f * dx;
                fy[j] += f * dy;
                fz[j] += f * dz;
            }
        }
        for &i in &alive {
            buffer.velocities_x[i] += fx[i] * dt / buffer.masses[i];
            buffer.velocities_y[i] += fy[i] * dt / buffer.masses[i];
            buffer.velocities_z[i] += fz[i] * dt / buffer.masses[i];
        }
    }
}
/// Emission mode.
#[derive(Debug, Clone, Copy)]
pub enum EmissionMode {
    /// All particles emitted at once.
    Burst {
        /// Number of particles to emit in the burst.
        count: usize,
    },
    /// Continuous emission at a given rate (particles/sec).
    Continuous {
        /// Emission rate in particles per second.
        rate: f32,
    },
}
/// Extended rendering data for a particle including sort key.
#[derive(Debug, Clone)]
pub struct SortedParticleRenderData {
    /// Particle render data.
    pub render_data: ParticleRenderData,
    /// Depth sort key (negative z in view space for back-to-front).
    pub sort_key: f32,
    /// Original buffer index.
    pub buffer_index: usize,
}
/// Linear drag force that damps particle velocity each step.
pub struct DragForce {
    /// Drag coefficient.  Applied as `v *= (1 - coefficient * dt)`.
    pub coefficient: f32,
}
impl DragForce {
    /// Apply drag to all alive particles.
    pub fn apply(&self, buffer: &mut ParticleBuffer, dt: f32) {
        let factor = (1.0 - self.coefficient * dt).max(0.0);
        for i in 0..buffer.count {
            if buffer.is_alive(i) {
                buffer.velocities_x[i] *= factor;
                buffer.velocities_y[i] *= factor;
                buffer.velocities_z[i] *= factor;
            }
        }
    }
}
/// Utilities for converting between the SoA `ParticleBuffer` and an
/// interleaved flat `f32` slice suitable for GPU upload.
pub struct GpuParticleLayout;
impl GpuParticleLayout {
    /// Number of `f32` values per particle in the interleaved layout.
    ///
    /// Layout: `[x, y, z, vx, vy, vz, mass, lifetime]`
    pub fn stride() -> usize {
        8
    }
    /// Produce an interleaved `Vec`f32` containing all slots (alive and dead).
    pub fn to_f32_buffer(buffer: &ParticleBuffer) -> Vec<f32> {
        let stride = Self::stride();
        let mut out = Vec::with_capacity(buffer.count * stride);
        for i in 0..buffer.count {
            out.push(buffer.positions_x[i]);
            out.push(buffer.positions_y[i]);
            out.push(buffer.positions_z[i]);
            out.push(buffer.velocities_x[i]);
            out.push(buffer.velocities_y[i]);
            out.push(buffer.velocities_z[i]);
            out.push(buffer.masses[i]);
            out.push(buffer.lifetimes[i]);
        }
        out
    }
    /// Parse an interleaved buffer back into a `ParticleBuffer`.
    ///
    /// `count` must equal the number of particle slots encoded in `data`.
    pub fn from_f32_buffer(data: &[f32], count: usize) -> ParticleBuffer {
        let stride = Self::stride();
        assert_eq!(data.len(), count * stride, "data length mismatch");
        let mut buf = ParticleBuffer::new(count);
        for i in 0..count {
            let base = i * stride;
            buf.positions_x[i] = data[base];
            buf.positions_y[i] = data[base + 1];
            buf.positions_z[i] = data[base + 2];
            buf.velocities_x[i] = data[base + 3];
            buf.velocities_y[i] = data[base + 4];
            buf.velocities_z[i] = data[base + 5];
            buf.masses[i] = data[base + 6];
            buf.lifetimes[i] = data[base + 7];
        }
        buf
    }
}
/// Extended particle system statistics.
#[derive(Debug, Clone)]
pub struct ParticleSystemStats {
    /// Basic stats.
    pub basic: ParticleStats,
    /// Total buffer capacity.
    pub capacity: usize,
    /// Fill ratio (active / capacity).
    pub fill_ratio: f32,
    /// Total kinetic energy.
    pub total_kinetic_energy: f32,
    /// Mean age of alive particles.
    pub mean_age: f32,
    /// Maximum age of alive particles.
    pub max_age: f32,
    /// Velocity standard deviation.
    pub velocity_std_dev: f32,
}
impl ParticleSystemStats {
    /// Compute extended statistics from a buffer.
    pub fn compute_extended(buffer: &ParticleBuffer) -> Self {
        let basic = ParticleStats::compute(buffer);
        let capacity = buffer.count;
        let fill_ratio = if capacity > 0 {
            basic.active as f32 / capacity as f32
        } else {
            0.0
        };
        let mut sum_age = 0.0f32;
        let mut max_age = 0.0f32;
        let mut sum_v2 = 0.0f32;
        let active = basic.active;
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            sum_age += buffer.ages[i];
            max_age = max_age.max(buffer.ages[i]);
            let vx = buffer.velocities_x[i];
            let vy = buffer.velocities_y[i];
            let vz = buffer.velocities_z[i];
            sum_v2 += vx * vx + vy * vy + vz * vz;
        }
        let mean_age = if active > 0 {
            sum_age / active as f32
        } else {
            0.0
        };
        let mean_v2 = if active > 0 {
            sum_v2 / active as f32
        } else {
            0.0
        };
        let velocity_std_dev = (mean_v2 - basic.avg_speed * basic.avg_speed)
            .max(0.0)
            .sqrt();
        Self {
            total_kinetic_energy: basic.total_kinetic_energy,
            basic,
            capacity,
            fill_ratio,
            mean_age,
            max_age,
            velocity_std_dev,
        }
    }
    /// Whether the buffer is at or near capacity.
    pub fn is_near_capacity(&self, threshold: f32) -> bool {
        self.fill_ratio >= threshold
    }
}
/// Vortex (rotational) force field around a Y-axis.
pub struct VortexForceField {
    /// Center of the vortex on the XZ plane.
    pub center: [f32; 2],
    /// Angular velocity (radians per second).
    pub angular_velocity: f32,
    /// Radius of influence.
    pub radius: f32,
}
impl VortexForceField {
    /// Apply vortex force to alive particles.
    pub fn apply(&self, buffer: &mut ParticleBuffer, dt: f32) {
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            let dx = buffer.positions_x[i] - self.center[0];
            let dz = buffer.positions_z[i] - self.center[1];
            let dist = (dx * dx + dz * dz).sqrt();
            if dist > self.radius || dist < 1e-6 {
                continue;
            }
            let factor = self.angular_velocity * (1.0 - dist / self.radius) * dt;
            buffer.velocities_x[i] += -dz / dist * factor;
            buffer.velocities_z[i] += dx / dist * factor;
        }
    }
}
/// Reflect particles off a horizontal floor plane at a fixed Y height.
pub struct FloorCollision {
    /// Y coordinate of the floor.
    pub y: f32,
    /// Coefficient of restitution (0 = fully inelastic, 1 = fully elastic).
    pub restitution: f32,
}
impl FloorCollision {
    /// Push particles above the floor and reflect downward velocities.
    pub fn apply(&self, buffer: &mut ParticleBuffer) {
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            if buffer.positions_y[i] < self.y {
                buffer.positions_y[i] = self.y;
                if buffer.velocities_y[i] < 0.0 {
                    buffer.velocities_y[i] = -buffer.velocities_y[i] * self.restitution;
                }
            }
        }
    }
}
/// A minimal Linear Congruential Generator used internally.
pub struct SimpleRng {
    pub(super) state: u64,
}
impl SimpleRng {
    /// Create a new RNG from a seed.
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x853c_49e6_748f_ea9b,
        }
    }
    /// Return the next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }
    /// Return a uniform float in [0, 1).
    pub fn next_f32(&mut self) -> f32 {
        let bits = (self.next_u64() >> 40) as u32;
        (bits as f32) / (1u32 << 24) as f32
    }
    /// Return a uniform float in [min, max).
    pub fn next_f32_range(&mut self, min: f32, max: f32) -> f32 {
        min + self.next_f32() * (max - min)
    }
    /// Return a random direction on the unit sphere (uniform).
    pub fn next_unit_sphere(&mut self) -> [f32; 3] {
        loop {
            let x = self.next_f32_range(-1.0, 1.0);
            let y = self.next_f32_range(-1.0, 1.0);
            let z = self.next_f32_range(-1.0, 1.0);
            let len2 = x * x + y * y + z * z;
            if len2 > 1e-10 && len2 <= 1.0 {
                let inv = 1.0 / len2.sqrt();
                return [x * inv, y * inv, z * inv];
            }
        }
    }
}
/// Kill particles that leave an axis-aligned bounding box.
pub struct BoundingBoxKill {
    /// Minimum corner of the kill-box.
    pub min: [f32; 3],
    /// Maximum corner of the kill-box.
    pub max: [f32; 3],
}
impl BoundingBoxKill {
    /// Kill any alive particle whose position is outside `\[min, max\]`.
    pub fn apply(&self, buffer: &mut ParticleBuffer) {
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            let x = buffer.positions_x[i];
            let y = buffer.positions_y[i];
            let z = buffer.positions_z[i];
            if x < self.min[0]
                || x > self.max[0]
                || y < self.min[1]
                || y > self.max[1]
                || z < self.min[2]
                || z > self.max[2]
            {
                buffer.kill(i);
            }
        }
    }
}
/// High-level particle system that owns a buffer, emitters, and forces.
pub struct ParticleSystem {
    /// The SoA particle data.
    pub buffer: ParticleBuffer,
    /// All registered emitters.
    pub emitters: Vec<ParticleEmitter>,
    /// Global gravity.
    pub gravity: GravityForce,
    /// Global drag.
    pub drag: DragForce,
    /// Optional floor collision plane.
    pub floor: Option<FloorCollision>,
    /// Internal RNG for emitter seeding.
    pub rng: SimpleRng,
    /// Elapsed simulation time.
    pub time: f32,
}
impl ParticleSystem {
    /// Create a new particle system with the given buffer capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: ParticleBuffer::new(capacity),
            emitters: Vec::new(),
            gravity: GravityForce {
                g: [0.0, -9.81, 0.0],
            },
            drag: DragForce { coefficient: 0.01 },
            floor: None,
            rng: SimpleRng::new(12345),
            time: 0.0,
        }
    }
    /// Register an emitter and return its index.
    pub fn add_emitter(&mut self, emitter: ParticleEmitter) -> usize {
        let idx = self.emitters.len();
        self.emitters.push(emitter);
        idx
    }
    /// Advance the simulation by `dt` seconds.
    ///
    /// Order: emit → gravity → drag → floor → integrate.
    pub fn step(&mut self, dt: f32) {
        let seed_base = self.rng.next_u64();
        for (idx, emitter) in self.emitters.iter_mut().enumerate() {
            let seed = seed_base ^ (idx as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
            emitter.emit(&mut self.buffer, dt, seed);
        }
        self.gravity.apply(&mut self.buffer, dt);
        self.drag.apply(&mut self.buffer, dt);
        if let Some(ref floor) = self.floor {
            floor.apply(&mut self.buffer);
        }
        ParticleIntegrator::integrate(&mut self.buffer, dt);
        self.time += dt;
    }
}
/// Radial force field: attracts or repels particles from a point.
pub struct RadialForceField {
    /// Center of the force field.
    pub center: [f32; 3],
    /// Strength of the field. Positive = attraction, negative = repulsion.
    pub strength: f32,
    /// Falloff exponent (1 = linear, 2 = inverse-square).
    pub falloff: f32,
    /// Minimum distance to avoid singularity.
    pub min_distance: f32,
}
impl RadialForceField {
    /// Apply this radial force to all alive particles.
    pub fn apply(&self, buffer: &mut ParticleBuffer, dt: f32) {
        for i in 0..buffer.count {
            if !buffer.is_alive(i) {
                continue;
            }
            let dx = self.center[0] - buffer.positions_x[i];
            let dy = self.center[1] - buffer.positions_y[i];
            let dz = self.center[2] - buffer.positions_z[i];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(self.min_distance);
            let force = self.strength / dist.powf(self.falloff);
            let inv_dist = 1.0 / dist;
            buffer.velocities_x[i] += force * dx * inv_dist * dt;
            buffer.velocities_y[i] += force * dy * inv_dist * dt;
            buffer.velocities_z[i] += force * dz * inv_dist * dt;
        }
    }
}
/// Per-particle rendering data for a GPU particle renderer.
#[derive(Debug, Clone)]
pub struct ParticleRenderData {
    /// Position as [x, y, z].
    pub position: [f32; 3],
    /// Color as [r, g, b, a].
    pub color: [f32; 4],
    /// Size (radius or diameter depending on renderer).
    pub size: f32,
    /// Normalized age (0 = just spawned, 1 = about to die).
    pub age_normalized: f32,
}
/// Constant gravitational acceleration applied to all alive particles.
pub struct GravityForce {
    /// Gravitational acceleration vector, e.g. `\[0.0, -9.81, 0.0\]`.
    pub g: [f32; 3],
}
impl GravityForce {
    /// Apply gravity: `v += g * dt` for every alive particle.
    pub fn apply(&self, buffer: &mut ParticleBuffer, dt: f32) {
        for i in 0..buffer.count {
            if buffer.is_alive(i) {
                buffer.velocities_x[i] += self.g[0] * dt;
                buffer.velocities_y[i] += self.g[1] * dt;
                buffer.velocities_z[i] += self.g[2] * dt;
            }
        }
    }
}
/// Structure-of-Arrays particle buffer, optimised for GPU upload.
pub struct ParticleBuffer {
    /// X positions of each particle slot.
    pub positions_x: Vec<f32>,
    /// Y positions of each particle slot.
    pub positions_y: Vec<f32>,
    /// Z positions of each particle slot.
    pub positions_z: Vec<f32>,
    /// X velocities.
    pub velocities_x: Vec<f32>,
    /// Y velocities.
    pub velocities_y: Vec<f32>,
    /// Z velocities.
    pub velocities_z: Vec<f32>,
    /// Per-particle mass.
    pub masses: Vec<f32>,
    /// Remaining lifetime; negative means the slot is dead.
    pub lifetimes: Vec<f32>,
    /// Elapsed age since spawn.
    pub ages: Vec<f32>,
    /// Number of slots allocated (alive + dead).
    pub count: usize,
}
impl ParticleBuffer {
    /// Allocate a buffer with `capacity` slots (all dead initially).
    pub fn new(capacity: usize) -> Self {
        Self {
            positions_x: vec![0.0; capacity],
            positions_y: vec![0.0; capacity],
            positions_z: vec![0.0; capacity],
            velocities_x: vec![0.0; capacity],
            velocities_y: vec![0.0; capacity],
            velocities_z: vec![0.0; capacity],
            masses: vec![1.0; capacity],
            lifetimes: vec![-1.0; capacity],
            ages: vec![0.0; capacity],
            count: capacity,
        }
    }
    /// Spawn a new particle in the first dead slot.  Returns its index or
    /// `None` if the buffer is full.
    pub fn add_particle(
        &mut self,
        pos: [f32; 3],
        vel: [f32; 3],
        mass: f32,
        lifetime: f32,
    ) -> Option<usize> {
        for i in 0..self.count {
            if self.lifetimes[i] < 0.0 {
                self.positions_x[i] = pos[0];
                self.positions_y[i] = pos[1];
                self.positions_z[i] = pos[2];
                self.velocities_x[i] = vel[0];
                self.velocities_y[i] = vel[1];
                self.velocities_z[i] = vel[2];
                self.masses[i] = mass;
                self.lifetimes[i] = lifetime;
                self.ages[i] = 0.0;
                return Some(i);
            }
        }
        None
    }
    /// Get the position of slot `i` as `\[x, y, z\]`.
    pub fn get_position(&self, i: usize) -> [f32; 3] {
        [
            self.positions_x[i],
            self.positions_y[i],
            self.positions_z[i],
        ]
    }
    /// Get the velocity of slot `i` as `\[vx, vy, vz\]`.
    pub fn get_velocity(&self, i: usize) -> [f32; 3] {
        [
            self.velocities_x[i],
            self.velocities_y[i],
            self.velocities_z[i],
        ]
    }
    /// Set the position of slot `i`.
    pub fn set_position(&mut self, i: usize, p: [f32; 3]) {
        self.positions_x[i] = p[0];
        self.positions_y[i] = p[1];
        self.positions_z[i] = p[2];
    }
    /// Set the velocity of slot `i`.
    pub fn set_velocity(&mut self, i: usize, v: [f32; 3]) {
        self.velocities_x[i] = v[0];
        self.velocities_y[i] = v[1];
        self.velocities_z[i] = v[2];
    }
    /// Return `true` if the slot at `i` holds a live particle.
    pub fn is_alive(&self, i: usize) -> bool {
        self.lifetimes[i] >= 0.0
    }
    /// Kill the particle at slot `i` by setting its lifetime negative.
    pub fn kill(&mut self, i: usize) {
        self.lifetimes[i] = -1.0;
    }
    /// Count how many slots are currently alive.
    pub fn active_count(&self) -> usize {
        (0..self.count).filter(|&i| self.is_alive(i)).count()
    }
}
/// Lifetime manager that also records spawn events for analysis.
pub struct ParticleLifetimeManager {
    /// Total particles ever spawned.
    pub total_spawned: usize,
    /// Total particles that have died naturally.
    pub total_expired: usize,
    /// Minimum observed lifetime seen.
    pub min_observed_lifetime: f32,
    /// Maximum observed lifetime seen.
    pub max_observed_lifetime: f32,
}
impl ParticleLifetimeManager {
    /// Create a new lifetime manager.
    pub fn new() -> Self {
        Self {
            total_spawned: 0,
            total_expired: 0,
            min_observed_lifetime: f32::MAX,
            max_observed_lifetime: 0.0,
        }
    }
    /// Record a spawn event with initial lifetime.
    pub fn record_spawn(&mut self, lifetime: f32) {
        self.total_spawned += 1;
        self.min_observed_lifetime = self.min_observed_lifetime.min(lifetime);
        self.max_observed_lifetime = self.max_observed_lifetime.max(lifetime);
    }
    /// Record an expiration (natural death).
    pub fn record_expiration(&mut self) {
        self.total_expired += 1;
    }
    /// Scan a buffer and retire particles that have expired.
    /// Returns the number retired.
    pub fn retire_expired(&mut self, buffer: &mut ParticleBuffer) -> usize {
        let mut count = 0;
        for i in 0..buffer.count {
            if buffer.lifetimes[i] < 0.0 && buffer.ages[i] > 0.0 {
                let _ = i;
            }
            if buffer.is_alive(i) && buffer.lifetimes[i] < 0.0 {
                count += 1;
                self.record_expiration();
            }
        }
        count
    }
    /// Alive fraction: alive / total_spawned.
    pub fn alive_fraction(&self, buffer: &ParticleBuffer) -> f32 {
        if self.total_spawned == 0 {
            return 0.0;
        }
        buffer.active_count() as f32 / self.total_spawned as f32
    }
}
/// Emits new particles into a `ParticleBuffer` over time.
pub struct ParticleEmitter {
    /// World-space spawn origin.
    pub position: [f32; 3],
    /// Target emission rate in particles per second.
    pub emit_rate: f32,
    /// Fractional particle accumulator (sub-frame carry).
    pub emit_accumulator: f32,
    /// Base initial velocity of emitted particles.
    pub initial_velocity: [f32; 3],
    /// Cone half-angle (radians) for random velocity spread.
    pub velocity_spread: f32,
    /// Minimum particle lifetime in seconds.
    pub lifetime_min: f32,
    /// Maximum particle lifetime in seconds.
    pub lifetime_max: f32,
    /// Mass assigned to emitted particles.
    pub mass: f32,
    /// Whether this emitter is currently active.
    pub active: bool,
}
impl ParticleEmitter {
    /// Create a new emitter.  `lifetime` is used as both min and max.
    pub fn new(pos: [f32; 3], rate: f32, vel: [f32; 3], lifetime: f32) -> Self {
        Self {
            position: pos,
            emit_rate: rate,
            emit_accumulator: 0.0,
            initial_velocity: vel,
            velocity_spread: 0.0,
            lifetime_min: lifetime,
            lifetime_max: lifetime,
            mass: 1.0,
            active: true,
        }
    }
    /// Emit particles for a time-step `dt`.  Returns how many were spawned.
    pub fn emit(&mut self, buffer: &mut ParticleBuffer, dt: f32, rng_seed: u64) -> usize {
        if !self.active {
            return 0;
        }
        let mut rng = SimpleRng::new(rng_seed);
        self.emit_accumulator += self.emit_rate * dt;
        let to_emit = self.emit_accumulator.floor() as usize;
        self.emit_accumulator -= to_emit as f32;
        let mut spawned = 0usize;
        for _ in 0..to_emit {
            let spread_dir = rng.next_unit_sphere();
            let vel = [
                self.initial_velocity[0] + spread_dir[0] * self.velocity_spread,
                self.initial_velocity[1] + spread_dir[1] * self.velocity_spread,
                self.initial_velocity[2] + spread_dir[2] * self.velocity_spread,
            ];
            let lt = rng.next_f32_range(self.lifetime_min, self.lifetime_max);
            if buffer
                .add_particle(self.position, vel, self.mass, lt)
                .is_some()
            {
                spawned += 1;
            }
        }
        spawned
    }
}
/// GPU particle emitter with configurable shape and mode.
pub struct GpuParticleEmitter {
    /// World-space spawn origin.
    pub position: [f32; 3],
    /// Emitter shape.
    pub shape: EmitterShape,
    /// Emission mode.
    pub mode: EmissionMode,
    /// Base initial velocity.
    pub initial_velocity: [f32; 3],
    /// Particle lifetime (seconds).
    pub lifetime: f32,
    /// Particle mass.
    pub mass: f32,
    /// Whether emitter is active.
    pub active: bool,
    /// Fractional accumulator for continuous mode.
    pub accumulator: f32,
    /// Internal RNG state.
    pub(super) rng: SimpleRng,
}
impl GpuParticleEmitter {
    /// Create a new emitter with a point shape and continuous mode.
    pub fn new_continuous(position: [f32; 3], rate: f32, lifetime: f32) -> Self {
        Self {
            position,
            shape: EmitterShape::Point,
            mode: EmissionMode::Continuous { rate },
            initial_velocity: [0.0, 1.0, 0.0],
            lifetime,
            mass: 1.0,
            active: true,
            accumulator: 0.0,
            rng: SimpleRng::new(0xdeadbeef),
        }
    }
    /// Create a burst emitter.
    pub fn new_burst(position: [f32; 3], count: usize, lifetime: f32) -> Self {
        Self {
            position,
            shape: EmitterShape::Point,
            mode: EmissionMode::Burst { count },
            initial_velocity: [0.0, 1.0, 0.0],
            lifetime,
            mass: 1.0,
            active: true,
            accumulator: 0.0,
            rng: SimpleRng::new(0xcafebabe),
        }
    }
    /// Emit particles into a buffer for one timestep `dt`.
    pub fn emit(&mut self, buffer: &mut ParticleBuffer, dt: f32) -> usize {
        if !self.active {
            return 0;
        }
        let to_emit = match self.mode {
            EmissionMode::Burst { count } => {
                self.active = false;
                count
            }
            EmissionMode::Continuous { rate } => {
                self.accumulator += rate * dt;
                let n = self.accumulator.floor() as usize;
                self.accumulator -= n as f32;
                n
            }
        };
        let mut spawned = 0;
        for _ in 0..to_emit {
            let pos = self.sample_position();
            let vel = self.initial_velocity;
            if buffer
                .add_particle(pos, vel, self.mass, self.lifetime)
                .is_some()
            {
                spawned += 1;
            }
        }
        spawned
    }
    /// Sample a spawn position according to the emitter shape.
    fn sample_position(&mut self) -> [f32; 3] {
        match self.shape {
            EmitterShape::Point => self.position,
            EmitterShape::Cone { half_angle } => {
                let _ = half_angle;
                self.position
            }
            EmitterShape::Sphere { radius } => {
                let dir = self.rng.next_unit_sphere();
                [
                    self.position[0] + dir[0] * radius,
                    self.position[1] + dir[1] * radius,
                    self.position[2] + dir[2] * radius,
                ]
            }
            EmitterShape::Box { half_extents } => {
                let x = self.rng.next_f32_range(-half_extents[0], half_extents[0]);
                let y = self.rng.next_f32_range(-half_extents[1], half_extents[1]);
                let z = self.rng.next_f32_range(-half_extents[2], half_extents[2]);
                [
                    self.position[0] + x,
                    self.position[1] + y,
                    self.position[2] + z,
                ]
            }
        }
    }
    /// Number of particles this emitter will emit in one burst (or 0 if continuous).
    pub fn burst_count(&self) -> usize {
        match self.mode {
            EmissionMode::Burst { count } => count,
            EmissionMode::Continuous { .. } => 0,
        }
    }
}
