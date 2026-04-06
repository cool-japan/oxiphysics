//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use oxiphysics_core::math::{Real, Vec3};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{Rgba, halton, value_noise_4d};

/// A turbulence field that applies pseudo-random velocity perturbations to
/// particles, simulating atmospheric or fluid turbulence.
///
/// The turbulence is computed using a layered (octave) approach inspired by
/// value noise, using deterministic integer-based hash functions.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TurbulenceField {
    /// Turbulence strength (m/s).
    pub strength: f64,
    /// Spatial frequency of the turbulence (1/m). Higher = finer detail.
    pub frequency: f64,
    /// Number of octaves to layer. More = richer structure.
    pub octaves: usize,
    /// Lacunarity: frequency multiplier per octave.
    pub lacunarity: f64,
    /// Persistence: amplitude multiplier per octave (0..1).
    pub persistence: f64,
    /// Time offset for animated turbulence.
    pub time_offset: f64,
}
#[allow(dead_code)]
impl TurbulenceField {
    /// Create a new turbulence field with default parameters.
    pub fn new(strength: f64, frequency: f64) -> Self {
        Self {
            strength,
            frequency,
            octaves: 3,
            lacunarity: 2.0,
            persistence: 0.5,
            time_offset: 0.0,
        }
    }
    /// Sample the turbulence velocity at position `pos` and time `t`.
    ///
    /// Returns a 3-vector perturbation velocity (m/s).
    pub fn sample(&self, pos: [f64; 3], t: f64) -> [f64; 3] {
        let mut vx = 0.0_f64;
        let mut vy = 0.0_f64;
        let mut vz = 0.0_f64;
        let mut amplitude = 1.0_f64;
        let mut freq = self.frequency;
        let mut amp_sum = 0.0_f64;
        for oct in 0..self.octaves {
            let seed = oct as u64;
            let px = pos[0] * freq;
            let py = pos[1] * freq;
            let pz = pos[2] * freq;
            let pt = (t + self.time_offset) * freq;
            vx += amplitude * value_noise_4d(px, py, pz, pt, seed);
            vy += amplitude * value_noise_4d(px + 31.4, py + 17.3, pz + 5.7, pt, seed + 1);
            vz += amplitude * value_noise_4d(px + 2.7, py + 41.1, pz + 93.2, pt, seed + 2);
            amp_sum += amplitude;
            amplitude *= self.persistence;
            freq *= self.lacunarity;
        }
        let inv_amp = if amp_sum > 1e-15 {
            self.strength / amp_sum
        } else {
            0.0
        };
        [vx * inv_amp, vy * inv_amp, vz * inv_amp]
    }
    /// Apply turbulence velocity impulse to all dynamic particles in a set.
    ///
    /// `dt` is the time step; the turbulence velocity is treated as an
    /// acceleration impulse: `Δv = F/m * dt` where `F = turbulence_vel * m`.
    pub fn apply(&self, ps: &mut ParticleSet, t: f64, dt: f64) {
        for p in &mut ps.particles {
            if p.is_static() {
                continue;
            }
            let turb = self.sample(p.position, t);
            p.velocity[0] += turb[0] * dt;
            p.velocity[1] += turb[1] * dt;
            p.velocity[2] += turb[2] * dt;
        }
    }
    /// Advance the time offset by `dt` (for animated turbulence).
    pub fn advance_time(&mut self, dt: f64) {
        self.time_offset += dt;
    }
}
/// A very simple spatial hash mapping 3-D grid cells → particle indices.
///
/// Cell size is fixed at construction time.
#[derive(Debug, Clone)]
pub(super) struct SpatialHash {
    /// Side length of one grid cell.
    pub(super) cell_size: f64,
    /// The map: cell key → list of particle indices.
    pub(super) table: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
}
impl SpatialHash {
    fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            table: std::collections::HashMap::new(),
        }
    }
    fn cell_of(&self, pos: [f64; 3]) -> (i64, i64, i64) {
        (
            pos[0].div_euclid(self.cell_size).floor() as i64,
            pos[1].div_euclid(self.cell_size).floor() as i64,
            pos[2].div_euclid(self.cell_size).floor() as i64,
        )
    }
    fn insert(&mut self, pos: [f64; 3], idx: usize) {
        let cell = self.cell_of(pos);
        self.table.entry(cell).or_default().push(idx);
    }
    #[allow(dead_code)]
    fn clear(&mut self) {
        self.table.clear();
    }
    /// Return all particle indices in cells overlapping the given sphere.
    fn query(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let r_cells = (radius / self.cell_size).ceil() as i64;
        let base = self.cell_of(center);
        let mut result = Vec::new();
        for dx in -r_cells..=r_cells {
            for dy in -r_cells..=r_cells {
                for dz in -r_cells..=r_cells {
                    let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                    if let Some(ids) = self.table.get(&key) {
                        result.extend_from_slice(ids);
                    }
                }
            }
        }
        result
    }
}
/// A single particle in a soft body simulation.
#[derive(Debug, Clone)]
pub struct SoftParticle {
    /// Current position.
    pub position: Vec3,
    /// Position from the previous time step (used by PBD).
    pub prev_position: Vec3,
    /// Current velocity.
    pub velocity: Vec3,
    /// Inverse mass (0.0 means the particle is fixed/static).
    pub inverse_mass: Real,
    /// Accumulated external force for this frame.
    pub external_force: Vec3,
}
impl SoftParticle {
    /// Create a new particle at the given position with the specified mass.
    ///
    /// If `mass` is zero or negative the particle is treated as static
    /// (infinite mass).
    pub fn new(position: Vec3, mass: Real) -> Self {
        let inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };
        Self {
            position,
            prev_position: position,
            velocity: Vec3::zeros(),
            inverse_mass: inv_mass,
            external_force: Vec3::zeros(),
        }
    }
    /// Create a static (immovable) particle at the given position.
    pub fn new_static(position: Vec3) -> Self {
        Self::new(position, 0.0)
    }
    /// Returns `true` when the particle has zero inverse mass (immovable).
    pub fn is_static(&self) -> bool {
        self.inverse_mass == 0.0
    }
}
/// Bit-flags for a [`Particle`] in a [`ParticleSet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ParticleFlags(pub u32);
impl ParticleFlags {
    /// The particle is static (pinned).
    pub const STATIC: u32 = 1 << 0;
    /// The particle is asleep (not being integrated).
    pub const ASLEEP: u32 = 1 << 1;
    /// The particle is a ghost / sensor (no collision response).
    pub const GHOST: u32 = 1 << 2;
    /// Create empty flags.
    #[allow(dead_code)]
    pub fn none() -> Self {
        Self(0)
    }
    /// Test whether a flag is set.
    #[allow(dead_code)]
    pub fn has(&self, flag: u32) -> bool {
        self.0 & flag != 0
    }
    /// Set a flag.
    #[allow(dead_code)]
    pub fn set(&mut self, flag: u32) {
        self.0 |= flag;
    }
    /// Clear a flag.
    #[allow(dead_code)]
    pub fn clear(&mut self, flag: u32) {
        self.0 &= !flag;
    }
}
/// Emits particles from the surface or volume of an axis-aligned box.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BoxEmitter {
    /// Minimum corner of the box.
    pub min: [f64; 3],
    /// Maximum corner of the box.
    pub max: [f64; 3],
    /// Emission velocity for all particles.
    pub velocity: [f64; 3],
    /// Particle mass.
    pub particle_mass: f64,
    /// Particle radius.
    pub particle_radius: f64,
    /// Emission rate (particles/s).
    pub rate: f64,
    /// Accumulated partial-particle counter.
    pub(super) accumulator: f64,
}
#[allow(dead_code)]
impl BoxEmitter {
    /// Create a new box emitter with given bounds.
    pub fn new(min: [f64; 3], max: [f64; 3], velocity: [f64; 3]) -> Self {
        Self {
            min,
            max,
            velocity,
            particle_mass: 1.0,
            particle_radius: 0.05,
            rate: 10.0,
            accumulator: 0.0,
        }
    }
    /// Emit `n` particles uniformly distributed inside the box.
    ///
    /// Uses deterministic low-discrepancy positions (Halton sequence base 2,3,5).
    pub fn emit_n(&self, ps: &mut ParticleSet, n: usize) {
        if n == 0 {
            return;
        }
        for k in 0..n {
            let hx = halton(k + 1, 2);
            let hy = halton(k + 1, 3);
            let hz = halton(k + 1, 5);
            let pos = [
                self.min[0] + hx * (self.max[0] - self.min[0]),
                self.min[1] + hy * (self.max[1] - self.min[1]),
                self.min[2] + hz * (self.max[2] - self.min[2]),
            ];
            let mut p = Particle::new(pos, self.particle_mass, self.particle_radius);
            p.velocity = self.velocity;
            ps.add_particle(p);
        }
    }
    /// Rate-based emit.
    pub fn emit(&mut self, ps: &mut ParticleSet, dt: f64) -> usize {
        self.accumulator += self.rate * dt;
        let n = self.accumulator.floor() as usize;
        self.accumulator -= n as f64;
        self.emit_n(ps, n);
        n
    }
    /// Volume of the emission box.
    pub fn volume(&self) -> f64 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        dx * dy * dz
    }
}
/// A particle in a [`ParticleSet`].
///
/// Uses plain `[f64; 3]` arrays (no nalgebra dependency in this struct).
#[derive(Debug, Clone)]
pub struct Particle {
    /// Current position (m).
    pub position: [f64; 3],
    /// Current velocity (m/s).
    pub velocity: [f64; 3],
    /// Particle mass (kg). Must be > 0 unless the particle is static.
    pub mass: f64,
    /// Particle radius (m) – used for spatial queries.
    pub radius: f64,
    /// Bit-flags.
    pub flags: ParticleFlags,
}
impl Particle {
    /// Create a new dynamic particle.
    #[allow(dead_code)]
    pub fn new(position: [f64; 3], mass: f64, radius: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            radius,
            flags: ParticleFlags::none(),
        }
    }
    /// Create a static (pinned) particle.
    #[allow(dead_code)]
    pub fn new_static(position: [f64; 3], radius: f64) -> Self {
        let mut p = Self::new(position, 1.0, radius);
        p.flags.set(ParticleFlags::STATIC);
        p
    }
    /// Returns `true` when the static flag is set.
    #[allow(dead_code)]
    pub fn is_static(&self) -> bool {
        self.flags.has(ParticleFlags::STATIC)
    }
    /// Kinetic energy of this particle (½mv²).
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        if self.is_static() {
            return 0.0;
        }
        let v2 = self.velocity[0] * self.velocity[0]
            + self.velocity[1] * self.velocity[1]
            + self.velocity[2] * self.velocity[2];
        0.5 * self.mass * v2
    }
    /// Linear momentum of this particle (mv).
    #[allow(dead_code)]
    pub fn momentum(&self) -> [f64; 3] {
        if self.is_static() {
            return [0.0; 3];
        }
        [
            self.mass * self.velocity[0],
            self.mass * self.velocity[1],
            self.mass * self.velocity[2],
        ]
    }
}
/// Emits particles along a line segment.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LineEmitter {
    /// Start point of the line.
    pub start: [f64; 3],
    /// End point of the line.
    pub end: [f64; 3],
    /// Emission velocity for all particles.
    pub velocity: [f64; 3],
    /// Particle mass.
    pub particle_mass: f64,
    /// Particle radius.
    pub particle_radius: f64,
}
#[allow(dead_code)]
impl LineEmitter {
    /// Create a new line emitter.
    pub fn new(start: [f64; 3], end: [f64; 3], velocity: [f64; 3]) -> Self {
        Self {
            start,
            end,
            velocity,
            particle_mass: 1.0,
            particle_radius: 0.1,
        }
    }
    /// Emit `n` particles evenly spaced along the line.
    pub fn emit(&self, ps: &mut ParticleSet, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        for i in 0..n {
            let t = if n == 1 {
                0.5
            } else {
                i as f64 / (n - 1) as f64
            };
            let pos = [
                self.start[0] + t * (self.end[0] - self.start[0]),
                self.start[1] + t * (self.end[1] - self.start[1]),
                self.start[2] + t * (self.end[2] - self.start[2]),
            ];
            let mut p = Particle::new(pos, self.particle_mass, self.particle_radius);
            p.velocity = self.velocity;
            ps.add_particle(p);
        }
        n
    }
}
/// Handles simple elastic collision response between particles.
#[allow(dead_code)]
pub struct ParticleCollisionHandler {
    /// Coefficient of restitution (0 = perfectly inelastic, 1 = perfectly elastic).
    pub restitution: f64,
}
#[allow(dead_code)]
impl ParticleCollisionHandler {
    /// Create a new collision handler.
    pub fn new(restitution: f64) -> Self {
        Self { restitution }
    }
    /// Resolve collisions between all overlapping particles.
    ///
    /// Returns the number of collisions resolved.
    pub fn resolve_collisions(&self, ps: &mut ParticleSet) -> usize {
        let n = ps.particles.len();
        let mut collision_count = 0;
        let data: Vec<([f64; 3], f64, f64, bool)> = ps
            .particles
            .iter()
            .map(|p| (p.position, p.radius, p.mass, p.is_static()))
            .collect();
        for i in 0..n {
            for j in (i + 1)..n {
                if data[i].3 && data[j].3 {
                    continue;
                }
                let dx = [
                    data[j].0[0] - data[i].0[0],
                    data[j].0[1] - data[i].0[1],
                    data[j].0[2] - data[i].0[2],
                ];
                let dist2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let min_dist = data[i].1 + data[j].1;
                if dist2 < min_dist * min_dist && dist2 > 1e-20 {
                    let dist = dist2.sqrt();
                    let normal = [dx[0] / dist, dx[1] / dist, dx[2] / dist];
                    let vi = ps.particles[i].velocity;
                    let vj = ps.particles[j].velocity;
                    let rel_vel = [vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]];
                    let vn =
                        rel_vel[0] * normal[0] + rel_vel[1] * normal[1] + rel_vel[2] * normal[2];
                    if vn > 0.0 {
                        continue;
                    }
                    let mi = data[i].2;
                    let mj = data[j].2;
                    let inv_mi = if data[i].3 { 0.0 } else { 1.0 / mi };
                    let inv_mj = if data[j].3 { 0.0 } else { 1.0 / mj };
                    let inv_total = inv_mi + inv_mj;
                    if inv_total < 1e-20 {
                        continue;
                    }
                    let j_impulse = -(1.0 + self.restitution) * vn / inv_total;
                    if !data[i].3 {
                        ps.particles[i].velocity[0] -= j_impulse * inv_mi * normal[0];
                        ps.particles[i].velocity[1] -= j_impulse * inv_mi * normal[1];
                        ps.particles[i].velocity[2] -= j_impulse * inv_mi * normal[2];
                    }
                    if !data[j].3 {
                        ps.particles[j].velocity[0] += j_impulse * inv_mj * normal[0];
                        ps.particles[j].velocity[1] += j_impulse * inv_mj * normal[1];
                        ps.particles[j].velocity[2] += j_impulse * inv_mj * normal[2];
                    }
                    collision_count += 1;
                }
            }
        }
        collision_count
    }
}
/// Evaluates a colour gradient over a particle's normalised age.
///
/// Keys must be sorted by `t` in ascending order.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ColorOverLifetime {
    /// Sorted keyframes.
    pub keys: Vec<ColorKey>,
}
#[allow(dead_code)]
impl ColorOverLifetime {
    /// Create a gradient from a list of `(t, rgba)` keyframes.
    ///
    /// Keyframes are automatically sorted by `t`.
    pub fn new(mut keys: Vec<(f32, Rgba)>) -> Self {
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self {
            keys: keys
                .into_iter()
                .map(|(t, color)| ColorKey { t, color })
                .collect(),
        }
    }
    /// White-to-transparent fade (useful default).
    pub fn white_to_transparent() -> Self {
        Self::new(vec![
            (0.0, [1.0, 1.0, 1.0, 1.0]),
            (1.0, [1.0, 1.0, 1.0, 0.0]),
        ])
    }
    /// Fire gradient: yellow → orange → red → black.
    pub fn fire() -> Self {
        Self::new(vec![
            (0.0, [1.0, 1.0, 0.0, 1.0]),
            (0.33, [1.0, 0.5, 0.0, 1.0]),
            (0.66, [1.0, 0.0, 0.0, 0.8]),
            (1.0, [0.1, 0.0, 0.0, 0.0]),
        ])
    }
    /// Evaluate the colour at normalised age `t` by linear interpolation.
    ///
    /// `t` is clamped to \[0, 1\].
    pub fn evaluate(&self, t: f32) -> Rgba {
        let t = t.clamp(0.0, 1.0);
        if self.keys.is_empty() {
            return [1.0, 1.0, 1.0, 1.0];
        }
        if t <= self.keys[0].t {
            return self.keys[0].color;
        }
        if t >= self.keys[self.keys.len() - 1].t {
            return self.keys[self.keys.len() - 1].color;
        }
        for i in 0..self.keys.len() - 1 {
            let k0 = &self.keys[i];
            let k1 = &self.keys[i + 1];
            if t >= k0.t && t <= k1.t {
                let span = k1.t - k0.t;
                let alpha = if span < 1e-7 { 0.0 } else { (t - k0.t) / span };
                return [
                    k0.color[0] + alpha * (k1.color[0] - k0.color[0]),
                    k0.color[1] + alpha * (k1.color[1] - k0.color[1]),
                    k0.color[2] + alpha * (k1.color[2] - k0.color[2]),
                    k0.color[3] + alpha * (k1.color[3] - k0.color[3]),
                ];
            }
        }
        self.keys[self.keys.len() - 1].color
    }
    /// Return the alpha value at normalised age `t`.
    pub fn alpha_at(&self, t: f32) -> f32 {
        self.evaluate(t)[3]
    }
}
/// Emits particles uniformly from within a cone volume.
///
/// Particles originate at `apex` and travel within `half_angle` of `axis`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConeEmitter {
    /// Cone apex (emission origin).
    pub apex: [f64; 3],
    /// Central axis of the cone (normalized).
    pub axis: [f64; 3],
    /// Half-angle of the cone spread (radians). 0 = pencil beam.
    pub half_angle: f64,
    /// Emission speed (m/s).
    pub speed: f64,
    /// Particle mass for emitted particles.
    pub particle_mass: f64,
    /// Particle radius for emitted particles.
    pub particle_radius: f64,
    /// Emission rate (particles/s).
    pub rate: f64,
    /// Accumulated partial-particle counter.
    pub(super) accumulator: f64,
}
#[allow(dead_code)]
impl ConeEmitter {
    /// Create a new cone emitter.
    pub fn new(apex: [f64; 3], axis: [f64; 3], half_angle: f64, speed: f64) -> Self {
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let axis = if len > 1e-12 {
            [axis[0] / len, axis[1] / len, axis[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            apex,
            axis,
            half_angle: half_angle.clamp(0.0, std::f64::consts::FRAC_PI_2),
            speed,
            particle_mass: 1.0,
            particle_radius: 0.1,
            rate: 10.0,
            accumulator: 0.0,
        }
    }
    /// Emit `n` particles with uniform cone-distributed velocities.
    ///
    /// Uses deterministic Halton-like sampling so tests are reproducible.
    pub fn emit_n(&self, ps: &mut ParticleSet, n: usize) {
        if n == 0 {
            return;
        }
        let ax = self.axis;
        let perp0 = if ax[0].abs() < 0.9 {
            let t = [1.0, 0.0, 0.0];
            let d = t[0] * ax[0] + t[1] * ax[1] + t[2] * ax[2];
            let raw = [t[0] - d * ax[0], t[1] - d * ax[1], t[2] - d * ax[2]];
            let l = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2])
                .sqrt()
                .max(1e-15);
            [raw[0] / l, raw[1] / l, raw[2] / l]
        } else {
            let t = [0.0, 1.0, 0.0];
            let d = t[0] * ax[0] + t[1] * ax[1] + t[2] * ax[2];
            let raw = [t[0] - d * ax[0], t[1] - d * ax[1], t[2] - d * ax[2]];
            let l = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2])
                .sqrt()
                .max(1e-15);
            [raw[0] / l, raw[1] / l, raw[2] / l]
        };
        let perp1 = [
            ax[1] * perp0[2] - ax[2] * perp0[1],
            ax[2] * perp0[0] - ax[0] * perp0[2],
            ax[0] * perp0[1] - ax[1] * perp0[0],
        ];
        for k in 0..n {
            let phi = 2.0 * std::f64::consts::PI * (k as f64 / n as f64);
            let cos_theta_max = self.half_angle.cos();
            let frac = (k as f64 + 0.5) / n as f64;
            let cos_theta = 1.0 - frac * (1.0 - cos_theta_max);
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(0.0);
            let dir = [
                ax[0] * cos_theta
                    + perp0[0] * sin_theta * phi.cos()
                    + perp1[0] * sin_theta * phi.sin(),
                ax[1] * cos_theta
                    + perp0[1] * sin_theta * phi.cos()
                    + perp1[1] * sin_theta * phi.sin(),
                ax[2] * cos_theta
                    + perp0[2] * sin_theta * phi.cos()
                    + perp1[2] * sin_theta * phi.sin(),
            ];
            let mut p = Particle::new(self.apex, self.particle_mass, self.particle_radius);
            p.velocity = [
                dir[0] * self.speed,
                dir[1] * self.speed,
                dir[2] * self.speed,
            ];
            ps.add_particle(p);
        }
    }
    /// Rate-based emit: call each frame with `dt`.
    ///
    /// Returns number of particles emitted this frame.
    pub fn emit(&mut self, ps: &mut ParticleSet, dt: f64) -> usize {
        self.accumulator += self.rate * dt;
        let n = self.accumulator.floor() as usize;
        self.accumulator -= n as f64;
        self.emit_n(ps, n);
        n
    }
    /// Solid angle of the cone: Ω = 2π(1 − cos θ).
    pub fn solid_angle(&self) -> f64 {
        2.0 * std::f64::consts::PI * (1.0 - self.half_angle.cos())
    }
}
/// A keyframe in a size-over-lifetime curve.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct SizeKey {
    /// Normalised age position (0 = birth, 1 = death).
    pub t: f32,
    /// Particle size multiplier at this key.
    pub size: f32,
}
/// A point attractor that pulls particles toward a point.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PointAttractor {
    /// Position of the attractor.
    pub position: [f64; 3],
    /// Strength of attraction (positive = attract, negative = repel).
    pub strength: f64,
    /// Maximum radius of influence.
    pub radius: f64,
}
#[allow(dead_code)]
impl PointAttractor {
    /// Create a new attractor.
    pub fn new(position: [f64; 3], strength: f64, radius: f64) -> Self {
        Self {
            position,
            strength,
            radius,
        }
    }
    /// Apply attraction force to all particles in the set.
    pub fn apply(&self, ps: &mut ParticleSet, dt: f64) {
        let r2_max = self.radius * self.radius;
        for p in &mut ps.particles {
            if p.is_static() {
                continue;
            }
            let dx = [
                self.position[0] - p.position[0],
                self.position[1] - p.position[1],
                self.position[2] - p.position[2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 > r2_max || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            let force_mag = self.strength / (r2 + 1e-6);
            let inv_m = 1.0 / p.mass;
            p.velocity[0] += force_mag * dx[0] / r * inv_m * dt;
            p.velocity[1] += force_mag * dx[1] / r * inv_m * dt;
            p.velocity[2] += force_mag * dx[2] / r * inv_m * dt;
        }
    }
}
/// A container for an SPH-like soft body simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SphSoftBody {
    /// Particles in the simulation.
    pub particles: Vec<SphSoftBodyParticle>,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Pressure stiffness k (Pa·m³/kg).
    pub stiffness: f64,
    /// Viscosity coefficient μ.
    pub viscosity: f64,
}
#[allow(dead_code)]
impl SphSoftBody {
    /// Create a new SPH soft body system.
    pub fn new(rest_density: f64, stiffness: f64, viscosity: f64) -> Self {
        Self {
            particles: Vec::new(),
            rest_density,
            stiffness,
            viscosity,
        }
    }
    /// Add a particle.
    pub fn add_particle(&mut self, p: SphSoftBodyParticle) {
        self.particles.push(p);
    }
    /// SPH cubic spline kernel W(q) with q = r/h.
    fn kernel_w(&self, r: f64, h: f64) -> f64 {
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
        let q = r / h;
        if q < 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q < 2.0 {
            sigma * 0.25 * (2.0 - q).powi(3)
        } else {
            0.0
        }
    }
    /// Gradient of SPH kernel: dW/dr.
    fn kernel_dw(&self, r: f64, h: f64) -> f64 {
        let sigma = 1.0 / (std::f64::consts::PI * h.powi(4));
        let q = r / h;
        if q < 1.0 {
            sigma * (-3.0 * q + 2.25 * q * q)
        } else if q < 2.0 {
            -sigma * 0.75 * (2.0 - q).powi(2)
        } else {
            0.0
        }
    }
    /// Compute density for all particles.
    pub fn compute_densities(&mut self) {
        let n = self.particles.len();
        for i in 0..n {
            let h = self.particles[i].h;
            let mut rho = 0.0;
            for j in 0..n {
                let dx = [
                    self.particles[i].position[0] - self.particles[j].position[0],
                    self.particles[i].position[1] - self.particles[j].position[1],
                    self.particles[i].position[2] - self.particles[j].position[2],
                ];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                rho += self.particles[j].mass * self.kernel_w(r, h);
            }
            self.particles[i].density = rho;
        }
    }
    /// Compute pressure from density via equation of state: P = k*(ρ/ρ₀ - 1).
    pub fn compute_pressures(&mut self) {
        for p in &mut self.particles {
            p.pressure = self.stiffness * (p.density / self.rest_density - 1.0);
        }
    }
    /// Compute pressure gradient forces.
    pub fn compute_forces(&mut self) {
        let n = self.particles.len();
        for p in &mut self.particles {
            p.force = [0.0; 3];
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = [
                    self.particles[i].position[0] - self.particles[j].position[0],
                    self.particles[i].position[1] - self.particles[j].position[1],
                    self.particles[i].position[2] - self.particles[j].position[2],
                ];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                if r < 1e-15 {
                    continue;
                }
                let h = (self.particles[i].h + self.particles[j].h) * 0.5;
                let dw = self.kernel_dw(r, h);
                let pi = self.particles[i].pressure;
                let pj = self.particles[j].pressure;
                let rhoi = self.particles[i].density.max(1e-15);
                let rhoj = self.particles[j].density.max(1e-15);
                let mj = self.particles[j].mass;
                let mi = self.particles[i].mass;
                let fmag = -mj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * dw;
                for d in 0..3 {
                    let rhat = dx[d] / r;
                    self.particles[i].force[d] += fmag * rhat * mi;
                    self.particles[j].force[d] -= fmag * rhat * mj;
                }
            }
        }
    }
    /// Integrate using explicit Euler.
    pub fn step(&mut self, dt: f64) {
        self.compute_densities();
        self.compute_pressures();
        self.compute_forces();
        for p in &mut self.particles {
            if p.pinned {
                continue;
            }
            let inv_m = 1.0 / p.mass;
            for d in 0..3 {
                p.velocity[d] += p.force[d] * inv_m * dt;
                p.position[d] += p.velocity[d] * dt;
            }
        }
    }
    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| !p.pinned)
            .map(|p| {
                let v2 = p.velocity[0] * p.velocity[0]
                    + p.velocity[1] * p.velocity[1]
                    + p.velocity[2] * p.velocity[2];
                0.5 * p.mass * v2
            })
            .sum()
    }
}
/// Emits particles from a point source with configurable velocity and spread.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PointEmitter {
    /// Emission origin.
    pub origin: [f64; 3],
    /// Base emission direction (normalized).
    pub direction: [f64; 3],
    /// Emission speed (m/s).
    pub speed: f64,
    /// Half-angle cone spread (radians, 0 = no spread).
    pub spread: f64,
    /// Particle mass for emitted particles.
    pub particle_mass: f64,
    /// Particle radius for emitted particles.
    pub particle_radius: f64,
    /// Emission rate (particles per second).
    pub rate: f64,
    /// Accumulated time since last emission.
    pub(super) accumulator: f64,
}
#[allow(dead_code)]
impl PointEmitter {
    /// Create a new point emitter.
    pub fn new(origin: [f64; 3], direction: [f64; 3], speed: f64) -> Self {
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        let dir = if len > 1e-12 {
            [direction[0] / len, direction[1] / len, direction[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            origin,
            direction: dir,
            speed,
            spread: 0.0,
            particle_mass: 1.0,
            particle_radius: 0.1,
            rate: 10.0,
            accumulator: 0.0,
        }
    }
    /// Emit particles into the given `ParticleSet`, advancing by `dt` seconds.
    ///
    /// Returns the number of particles emitted.
    pub fn emit(&mut self, ps: &mut ParticleSet, dt: f64) -> usize {
        self.accumulator += dt * self.rate;
        let n = self.accumulator.floor() as usize;
        self.accumulator -= n as f64;
        for _ in 0..n {
            let vel = [
                self.direction[0] * self.speed,
                self.direction[1] * self.speed,
                self.direction[2] * self.speed,
            ];
            let mut p = Particle::new(self.origin, self.particle_mass, self.particle_radius);
            p.velocity = vel;
            ps.add_particle(p);
        }
        n
    }
}
/// Converts particle data to a regular grid (scalar field).
#[allow(dead_code)]
pub struct ParticleToGrid {
    /// Grid origin (minimum corner).
    pub origin: [f64; 3],
    /// Cell size.
    pub cell_size: f64,
    /// Grid dimensions \[nx, ny, nz\].
    pub dims: [usize; 3],
}
#[allow(dead_code)]
impl ParticleToGrid {
    /// Create a new particle-to-grid converter.
    pub fn new(origin: [f64; 3], cell_size: f64, dims: [usize; 3]) -> Self {
        Self {
            origin,
            cell_size,
            dims,
        }
    }
    /// Scatter particle masses onto the grid using Cloud-In-Cell (CIC) weighting.
    ///
    /// Returns a flat array of size `nx * ny * nz`.
    pub fn scatter_mass(&self, particles: &[Particle]) -> Vec<f64> {
        let [nx, ny, nz] = self.dims;
        let mut grid = vec![0.0_f64; nx * ny * nz];
        for p in particles {
            let fx = (p.position[0] - self.origin[0]) / self.cell_size;
            let fy = (p.position[1] - self.origin[1]) / self.cell_size;
            let fz = (p.position[2] - self.origin[2]) / self.cell_size;
            let ix = fx.floor() as isize;
            let iy = fy.floor() as isize;
            let iz = fz.floor() as isize;
            let tx = fx - ix as f64;
            let ty = fy - iy as f64;
            let tz = fz - iz as f64;
            for di in 0..2_isize {
                for dj in 0..2_isize {
                    for dk in 0..2_isize {
                        let gi = ix + di;
                        let gj = iy + dj;
                        let gk = iz + dk;
                        if gi < 0 || gj < 0 || gk < 0 {
                            continue;
                        }
                        let gi = gi as usize;
                        let gj = gj as usize;
                        let gk = gk as usize;
                        if gi >= nx || gj >= ny || gk >= nz {
                            continue;
                        }
                        let wx = if di == 0 { 1.0 - tx } else { tx };
                        let wy = if dj == 0 { 1.0 - ty } else { ty };
                        let wz = if dk == 0 { 1.0 - tz } else { tz };
                        let idx = gi * ny * nz + gj * nz + gk;
                        grid[idx] += p.mass * wx * wy * wz;
                    }
                }
            }
        }
        grid
    }
    /// Total grid size (number of cells).
    pub fn total_cells(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
}
/// A keyframe in a colour-over-lifetime gradient.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct ColorKey {
    /// Normalised age position (0 = birth, 1 = death).
    pub t: f32,
    /// RGBA colour at this key.
    pub color: Rgba,
}
/// Evaluates particle size over its normalised lifetime.
///
/// Useful for making particles grow, shrink, or pulsate.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SizeOverLifetime {
    /// Sorted keyframes.
    pub keys: Vec<SizeKey>,
}
#[allow(dead_code)]
impl SizeOverLifetime {
    /// Create a size curve from a list of `(t, size)` keyframes.
    pub fn new(mut keys: Vec<(f32, f32)>) -> Self {
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self {
            keys: keys
                .into_iter()
                .map(|(t, size)| SizeKey { t, size })
                .collect(),
        }
    }
    /// Constant size = 1 for the entire lifetime.
    pub fn constant() -> Self {
        Self::new(vec![(0.0, 1.0), (1.0, 1.0)])
    }
    /// Grow from 0 to 1, hold, then shrink back to 0.
    pub fn grow_shrink() -> Self {
        Self::new(vec![(0.0, 0.0), (0.2, 1.0), (0.8, 1.0), (1.0, 0.0)])
    }
    /// Evaluate the size multiplier at normalised age `t`.
    pub fn evaluate(&self, t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        if self.keys.is_empty() {
            return 1.0;
        }
        if t <= self.keys[0].t {
            return self.keys[0].size;
        }
        if t >= self.keys[self.keys.len() - 1].t {
            return self.keys[self.keys.len() - 1].size;
        }
        for i in 0..self.keys.len() - 1 {
            let k0 = &self.keys[i];
            let k1 = &self.keys[i + 1];
            if t >= k0.t && t <= k1.t {
                let span = k1.t - k0.t;
                let alpha = if span < 1e-7 { 0.0 } else { (t - k0.t) / span };
                return k0.size + alpha * (k1.size - k0.size);
            }
        }
        self.keys[self.keys.len() - 1].size
    }
    /// Minimum size over the whole lifetime.
    pub fn min_size(&self) -> f32 {
        self.keys
            .iter()
            .map(|k| k.size)
            .fold(f32::INFINITY, f32::min)
    }
    /// Maximum size over the whole lifetime.
    pub fn max_size(&self) -> f32 {
        self.keys.iter().map(|k| k.size).fold(0.0_f32, f32::max)
    }
}
/// A collection of particles that form a deformable body.
#[derive(Debug, Clone)]
pub struct SoftBody {
    /// The particles that make up this body.
    pub particles: Vec<SoftParticle>,
    /// Triangle connectivity (indices into `particles`).
    pub triangles: Vec<[usize; 3]>,
    /// Tetrahedron connectivity (indices into `particles`).
    pub tetrahedra: Vec<[usize; 4]>,
    /// Velocity damping factor applied each sub-step (0 = no damping, 1 = full).
    pub damping: Real,
}
impl SoftBody {
    /// Create an empty soft body.
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
            triangles: Vec::new(),
            tetrahedra: Vec::new(),
            damping: 0.01,
        }
    }
    /// Create a soft body from existing particles.
    pub fn from_particles(particles: Vec<SoftParticle>) -> Self {
        Self {
            particles,
            triangles: Vec::new(),
            tetrahedra: Vec::new(),
            damping: 0.01,
        }
    }
    /// Apply an external force (e.g. gravity) to all dynamic particles.
    pub fn apply_force(&mut self, force: &Vec3) {
        for p in &mut self.particles {
            if !p.is_static() {
                p.external_force += force;
            }
        }
    }
    /// Clear all accumulated external forces.
    pub fn clear_forces(&mut self) {
        for p in &mut self.particles {
            p.external_force = Vec3::zeros();
        }
    }
}
/// Sample a scalar field at arbitrary positions using particle data.
#[allow(dead_code)]
pub struct ParticleFieldSampler;
#[allow(dead_code)]
impl ParticleFieldSampler {
    /// Sample the density field at position `pos` using SPH kernel.
    ///
    /// Uses a simple Gaussian-like falloff: w(r) = exp(-r²/(2h²)).
    pub fn sample_density(particles: &[Particle], pos: [f64; 3], h: f64) -> f64 {
        let h2 = h * h;
        let mut density = 0.0;
        for p in particles {
            let dx = [
                pos[0] - p.position[0],
                pos[1] - p.position[1],
                pos[2] - p.position[2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 < 9.0 * h2 {
                let w = (-r2 / (2.0 * h2)).exp();
                density += p.mass * w;
            }
        }
        density
    }
    /// Sample the velocity field at position `pos`.
    pub fn sample_velocity(particles: &[Particle], pos: [f64; 3], h: f64) -> [f64; 3] {
        let h2 = h * h;
        let mut vel = [0.0_f64; 3];
        let mut w_sum = 0.0_f64;
        for p in particles {
            let dx = [
                pos[0] - p.position[0],
                pos[1] - p.position[1],
                pos[2] - p.position[2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 < 9.0 * h2 {
                let w = (-r2 / (2.0 * h2)).exp();
                vel[0] += w * p.velocity[0];
                vel[1] += w * p.velocity[1];
                vel[2] += w * p.velocity[2];
                w_sum += w;
            }
        }
        if w_sum > 1e-30 {
            vel[0] /= w_sum;
            vel[1] /= w_sum;
            vel[2] /= w_sum;
        }
        vel
    }
}
/// Simple k-means-style particle cluster.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParticleCluster {
    /// Cluster centroid.
    pub centroid: [f64; 3],
    /// Indices of particles in this cluster.
    pub members: Vec<usize>,
    /// Total mass of the cluster.
    pub total_mass: f64,
}
#[allow(dead_code)]
impl ParticleCluster {
    /// Create a cluster with initial centroid and no members.
    pub fn new(centroid: [f64; 3]) -> Self {
        Self {
            centroid,
            members: Vec::new(),
            total_mass: 0.0,
        }
    }
    /// Update the centroid based on current member positions.
    pub fn update_centroid(&mut self, particles: &[Particle]) {
        if self.members.is_empty() {
            return;
        }
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        let mut tm = 0.0_f64;
        for &i in &self.members {
            cx += particles[i].position[0] * particles[i].mass;
            cy += particles[i].position[1] * particles[i].mass;
            cz += particles[i].position[2] * particles[i].mass;
            tm += particles[i].mass;
        }
        if tm > 1e-15 {
            self.centroid = [cx / tm, cy / tm, cz / tm];
            self.total_mass = tm;
        }
    }
}
/// A pool of particles with lifetime management.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct LifeParticlePool {
    /// Active particles.
    pub particles: Vec<LifeParticle>,
}
#[allow(dead_code)]
impl LifeParticlePool {
    /// Create an empty pool.
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
        }
    }
    /// Spawn a new particle.
    pub fn spawn(&mut self, particle: Particle, lifetime: f64) {
        self.particles.push(LifeParticle::new(particle, lifetime));
    }
    /// Step the pool: age all particles, remove dead ones.
    pub fn step(&mut self, dt: f64) {
        for lp in &mut self.particles {
            if lp.is_alive() {
                lp.age(dt);
                if !lp.particle.is_static() {
                    for d in 0..3 {
                        lp.particle.position[d] += lp.particle.velocity[d] * dt;
                    }
                }
            }
        }
        self.particles.retain(|lp| lp.is_alive());
    }
    /// Number of alive particles.
    pub fn alive_count(&self) -> usize {
        self.particles.len()
    }
}
/// Particle with a finite lifetime.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LifeParticle {
    /// Underlying particle data.
    pub particle: Particle,
    /// Remaining lifetime (s). 0 or below means dead.
    pub lifetime: f64,
    /// Maximum lifetime at birth (for normalizing age).
    pub max_lifetime: f64,
}
#[allow(dead_code)]
impl LifeParticle {
    /// Create a new life particle.
    pub fn new(particle: Particle, lifetime: f64) -> Self {
        Self {
            particle,
            lifetime,
            max_lifetime: lifetime,
        }
    }
    /// Age fraction: 0 = just born, 1 = just died.
    pub fn age_fraction(&self) -> f64 {
        if self.max_lifetime < 1e-15 {
            return 1.0;
        }
        1.0 - (self.lifetime / self.max_lifetime).clamp(0.0, 1.0)
    }
    /// Returns `true` if the particle is still alive.
    pub fn is_alive(&self) -> bool {
        self.lifetime > 0.0
    }
    /// Advance lifetime by `dt` seconds.
    pub fn age(&mut self, dt: f64) {
        self.lifetime -= dt;
    }
}
/// Emits particles from the surface or volume of a sphere.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SphereEmitter {
    /// Centre of the emission sphere.
    pub center: [f64; 3],
    /// Radius of the emission sphere (m).
    pub radius: f64,
    /// Emission speed (m/s). Velocity direction = outward normal.
    pub speed: f64,
    /// If `true`, emit from volume; if `false`, emit from surface only.
    pub volume_emit: bool,
    /// Particle mass.
    pub particle_mass: f64,
    /// Particle radius.
    pub particle_radius: f64,
    /// Emission rate (particles/s).
    pub rate: f64,
    /// Accumulated partial-particle counter.
    pub(super) accumulator: f64,
}
#[allow(dead_code)]
impl SphereEmitter {
    /// Create a new sphere emitter (surface mode by default).
    pub fn new(center: [f64; 3], radius: f64, speed: f64) -> Self {
        Self {
            center,
            radius,
            speed,
            volume_emit: false,
            particle_mass: 1.0,
            particle_radius: 0.05,
            rate: 10.0,
            accumulator: 0.0,
        }
    }
    /// Emit `n` particles uniformly distributed on/in the sphere.
    ///
    /// Uses the Fibonacci sphere sampling for uniform surface distribution.
    pub fn emit_n(&self, ps: &mut ParticleSet, n: usize) {
        if n == 0 {
            return;
        }
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        for k in 0..n {
            let y = 1.0 - 2.0 * (k as f64 + 0.5) / n as f64;
            let r_xy = (1.0 - y * y).sqrt().max(0.0);
            let theta = golden * k as f64;
            let nx = r_xy * theta.cos();
            let nz = r_xy * theta.sin();
            let normal = [nx, y, nz];
            let r_emit = if self.volume_emit {
                let u = (k as f64 + 0.5) / n as f64;
                self.radius * u.cbrt()
            } else {
                self.radius
            };
            let pos = [
                self.center[0] + normal[0] * r_emit,
                self.center[1] + normal[1] * r_emit,
                self.center[2] + normal[2] * r_emit,
            ];
            let mut p = Particle::new(pos, self.particle_mass, self.particle_radius);
            p.velocity = [
                normal[0] * self.speed,
                normal[1] * self.speed,
                normal[2] * self.speed,
            ];
            ps.add_particle(p);
        }
    }
    /// Rate-based emit.
    pub fn emit(&mut self, ps: &mut ParticleSet, dt: f64) -> usize {
        self.accumulator += self.rate * dt;
        let n = self.accumulator.floor() as usize;
        self.accumulator -= n as f64;
        self.emit_n(ps, n);
        n
    }
    /// Surface area of the emission sphere.
    pub fn surface_area(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }
}
/// SPH-like soft body with density-based pressure forces.
///
/// Each `SphSoftBodyParticle` tracks position, velocity, density, pressure, and
/// accumulated force. The system integrates using Leapfrog and applies
/// neighbour-based pressure gradients.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SphSoftBodyParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing radius (m).
    pub h: f64,
    /// SPH density estimate (kg/m^3).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Accumulated force for this frame.
    pub force: [f64; 3],
    /// Is this particle pinned (static)?
    pub pinned: bool,
}
#[allow(dead_code)]
impl SphSoftBodyParticle {
    /// Create a new SPH soft body particle.
    pub fn new(position: [f64; 3], mass: f64, h: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            h,
            density: 0.0,
            pressure: 0.0,
            force: [0.0; 3],
            pinned: false,
        }
    }
    /// Squared distance to another particle.
    pub fn dist_sq(&self, other: &Self) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        dx * dx + dy * dy + dz * dz
    }
}
/// A managed collection of [`Particle`]s with spatial hashing for fast
/// neighbour queries.
///
/// # Example
/// ```rust,ignore
/// let mut ps = ParticleSet::new(0.5);
/// ps.add_particle(Particle::new([0.0, 0.0, 0.0], 1.0, 0.1));
/// ps.add_particle(Particle::new([0.3, 0.0, 0.0], 1.0, 0.1));
/// let neighbours = ps.neighbors_within([0.0, 0.0, 0.0], 0.5);
/// ```
#[derive(Debug, Clone)]
pub struct ParticleSet {
    /// The particles.
    pub particles: Vec<Particle>,
    /// Spatial hash (rebuilt on demand).
    pub(super) hash: SpatialHash,
    /// Cell size used for the spatial hash (should match the typical query radius).
    pub(super) cell_size: f64,
    /// Whether the spatial hash is up to date.
    pub(super) hash_dirty: bool,
}
impl ParticleSet {
    /// Create an empty particle set.
    ///
    /// `cell_size` controls the spatial hash resolution.  A good value is the
    /// typical neighbour-query radius.
    #[allow(dead_code)]
    pub fn new(cell_size: f64) -> Self {
        Self {
            particles: Vec::new(),
            hash: SpatialHash::new(cell_size),
            cell_size,
            hash_dirty: false,
        }
    }
    /// Add a particle and return its index.
    #[allow(dead_code)]
    pub fn add_particle(&mut self, p: Particle) -> usize {
        let idx = self.particles.len();
        self.hash.insert(p.position, idx);
        self.particles.push(p);
        idx
    }
    /// Remove the particle at `idx` by swapping with the last element.
    ///
    /// Note: this changes the index of the last particle.  Returns the swapped
    /// particle's old index (which is now `idx`), or `None` if `idx` was the
    /// last particle.
    #[allow(dead_code)]
    pub fn remove_particle(&mut self, idx: usize) -> Option<usize> {
        let n = self.particles.len();
        if idx >= n {
            return None;
        }
        self.particles.swap_remove(idx);
        self.hash_dirty = true;
        if idx < self.particles.len() {
            Some(n - 1)
        } else {
            None
        }
    }
    /// Rebuild the spatial hash from scratch.
    #[allow(dead_code)]
    pub fn rebuild_hash(&mut self) {
        self.hash = SpatialHash::new(self.cell_size);
        for (i, p) in self.particles.iter().enumerate() {
            self.hash.insert(p.position, i);
        }
        self.hash_dirty = false;
    }
    /// Ensure the spatial hash is up to date.
    fn ensure_hash(&mut self) {
        if self.hash_dirty {
            self.rebuild_hash();
        }
    }
    /// Return the indices of all particles within `radius` of `pos`.
    ///
    /// The query is exact: only particles whose centre is within `radius` are
    /// returned (the hash is used as a broad-phase).
    #[allow(dead_code)]
    pub fn neighbors_within(&mut self, pos: [f64; 3], radius: f64) -> Vec<usize> {
        self.ensure_hash();
        let candidates = self.hash.query(pos, radius);
        let r2 = radius * radius;
        candidates
            .into_iter()
            .filter(|&i| {
                let p = &self.particles[i].position;
                let dx = p[0] - pos[0];
                let dy = p[1] - pos[1];
                let dz = p[2] - pos[2];
                dx * dx + dy * dy + dz * dz <= r2
            })
            .collect()
    }
    /// Apply an external force to particle `idx` for one frame.
    ///
    /// The force is integrated by the caller; here we directly modify velocity
    /// using an implicit Euler half-step `Δv = F/m * dt`.
    #[allow(dead_code)]
    pub fn apply_external_force(&mut self, idx: usize, force: [f64; 3], dt: f64) {
        let p = &mut self.particles[idx];
        if p.is_static() {
            return;
        }
        let inv_m = 1.0 / p.mass;
        p.velocity[0] += force[0] * inv_m * dt;
        p.velocity[1] += force[1] * inv_m * dt;
        p.velocity[2] += force[2] * inv_m * dt;
    }
    /// Apply an instantaneous impulse `J` (kg·m/s) to particle `idx`.
    ///
    /// `Δv = J / m`
    #[allow(dead_code)]
    pub fn apply_impulse(&mut self, idx: usize, impulse: [f64; 3]) {
        let p = &mut self.particles[idx];
        if p.is_static() {
            return;
        }
        let inv_m = 1.0 / p.mass;
        p.velocity[0] += impulse[0] * inv_m;
        p.velocity[1] += impulse[1] * inv_m;
        p.velocity[2] += impulse[2] * inv_m;
    }
    /// Total kinetic energy of the system (½ Σ mᵢ |vᵢ|²).
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
    /// Total linear momentum of the system (Σ mᵢ vᵢ).
    #[allow(dead_code)]
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut m = [0.0_f64; 3];
        for p in &self.particles {
            let pm = p.momentum();
            m[0] += pm[0];
            m[1] += pm[1];
            m[2] += pm[2];
        }
        m
    }
    /// Step all dynamic particles forward by `dt` (simple explicit Euler).
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        for p in &mut self.particles {
            if p.is_static() || p.flags.has(ParticleFlags::ASLEEP) {
                continue;
            }
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
        }
        self.hash_dirty = true;
    }
}
