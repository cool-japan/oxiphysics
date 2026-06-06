// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle effect rendering data structures.
//!
//! Provides single-particle effects, fire and smoke emitters, explosion bursts,
//! a pool-based particle manager, billboard quad generation, and back-to-front
//! sorting for correct alpha blending.

use rand::RngExt;

// ─────────────────────────────────────────────────────────────────────────────
// ParticleEffect
// ─────────────────────────────────────────────────────────────────────────────

/// A single particle effect with position, velocity, color, and lifetime state.
#[derive(Debug, Clone)]
pub struct ParticleEffect {
    /// Total lifetime of this particle in seconds.
    pub lifetime: f64,
    /// Current age of this particle in seconds.
    pub age: f64,
    /// World-space position `[x, y, z]`.
    pub pos: [f32; 3],
    /// Velocity vector `[vx, vy, vz]`.
    pub vel: [f32; 3],
    /// RGBA color.
    pub color: [f32; 4],
    /// Billboard size in world units.
    pub size: f32,
    /// Current opacity (0 = transparent, 1 = opaque).
    pub opacity: f32,
}

impl ParticleEffect {
    /// Create a new particle.
    ///
    /// The particle starts at `pos`, moves with `vel`, uses `color`, lives for
    /// `lifetime` seconds, and has the given billboard `size`.
    pub fn new(pos: [f32; 3], vel: [f32; 3], color: [f32; 4], lifetime: f64, size: f32) -> Self {
        Self {
            lifetime,
            age: 0.0,
            pos,
            vel,
            color,
            size,
            opacity: 1.0,
        }
    }

    /// Advance the particle by `dt` seconds: move position, age, and fade opacity.
    pub fn update(&mut self, dt: f32) {
        self.age += dt as f64;
        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
        self.pos[2] += self.vel[2] * dt;
        let t = self.normalized_age();
        self.opacity = (1.0 - t).max(0.0);
    }

    /// Return `true` if the particle's age is less than its lifetime.
    pub fn is_alive(&self) -> bool {
        self.age < self.lifetime
    }

    /// Current alpha value (same as `opacity`).
    pub fn alpha(&self) -> f32 {
        self.opacity
    }

    /// Normalized age in `[0, 1]` (0 = just spawned, 1 = at end of lifetime).
    pub fn normalized_age(&self) -> f32 {
        if self.lifetime <= 0.0 {
            1.0
        } else {
            (self.age / self.lifetime).min(1.0) as f32
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FireEmitter
// ─────────────────────────────────────────────────────────────────────────────

/// A fire particle emitter.
///
/// Emits particles that rise with randomised velocity; color interpolates from
/// hot (near base) to cool (at tip) over each particle's lifetime.
#[derive(Debug, Clone)]
pub struct FireEmitter {
    /// World-space origin of the fire `[x, y, z]`.
    pub pos: [f32; 3],
    /// Particles emitted per second.
    pub rate: f32,
    /// Intensity multiplier affecting velocity and scale.
    pub intensity: f32,
    /// RGBA color for the hottest part of the flame.
    pub color_hot: [f32; 4],
    /// RGBA color for the coolest part of the flame (tips/edges).
    pub color_cool: [f32; 4],
}

impl FireEmitter {
    /// Create a new fire emitter at `pos` with given `intensity`.
    ///
    /// Hot color defaults to orange-yellow, cool color to dark red.
    pub fn new(pos: [f32; 3], intensity: f32) -> Self {
        Self {
            pos,
            rate: 20.0 * intensity,
            intensity,
            color_hot: [1.0, 0.8, 0.1, 1.0],
            color_cool: [0.6, 0.05, 0.0, 0.0],
        }
    }

    /// Emit particles for a time step `dt`.
    ///
    /// The number of particles emitted is `floor(rate * dt)` (minimum 0).
    pub fn emit(&self, dt: f32) -> Vec<ParticleEffect> {
        use rand::RngExt;
        let mut rng = rand::rng();
        let count = (self.rate * dt).floor() as usize;
        let mut particles = Vec::with_capacity(count);
        for _ in 0..count {
            let vx: f32 = rng.random_range(-0.5_f32..0.5_f32) * self.intensity;
            let vy: f32 = rng.random_range(1.0_f32..3.0_f32) * self.intensity;
            let vz: f32 = rng.random_range(-0.5_f32..0.5_f32) * self.intensity;
            let lifetime = rng.random_range(0.5_f64..1.5_f64);
            let size = rng.random_range(0.05_f32..0.2_f32) * self.intensity;
            particles.push(ParticleEffect::new(
                self.pos,
                [vx, vy, vz],
                self.color_hot,
                lifetime,
                size,
            ));
        }
        particles
    }

    /// Interpolate color between hot and cool based on normalized age `t ∈ [0, 1]`.
    pub fn color_at_age(&self, t: f32) -> [f32; 4] {
        let t = t.clamp(0.0, 1.0);
        [
            self.color_hot[0] * (1.0 - t) + self.color_cool[0] * t,
            self.color_hot[1] * (1.0 - t) + self.color_cool[1] * t,
            self.color_hot[2] * (1.0 - t) + self.color_cool[2] * t,
            self.color_hot[3] * (1.0 - t) + self.color_cool[3] * t,
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SmokeEmitter
// ─────────────────────────────────────────────────────────────────────────────

/// A smoke particle emitter.
///
/// Emits slow-rising, wide-spreading grey particles.
#[derive(Debug, Clone)]
pub struct SmokeEmitter {
    /// World-space origin `[x, y, z]`.
    pub pos: [f32; 3],
    /// Particles emitted per second.
    pub rate: f32,
    /// Horizontal spread radius.
    pub spread: f32,
    /// Upward drift speed.
    pub rise_speed: f32,
}

impl SmokeEmitter {
    /// Create a smoke emitter at `pos` with given `spread` radius.
    pub fn new(pos: [f32; 3], spread: f32) -> Self {
        Self {
            pos,
            rate: 5.0,
            spread,
            rise_speed: 0.3,
        }
    }

    /// Emit smoke particles for a time step `dt`.
    pub fn emit(&self, dt: f32) -> Vec<ParticleEffect> {
        let mut rng = rand::rng();
        let count = (self.rate * dt).floor().max(0.0) as usize;
        let mut particles = Vec::with_capacity(count);
        for _ in 0..count {
            let vx: f32 = rng.random_range(-self.spread..self.spread);
            let vy: f32 = rng.random_range(self.rise_speed * 0.5..self.rise_speed * 1.5);
            let vz: f32 = rng.random_range(-self.spread..self.spread);
            let lifetime = rng.random_range(2.0_f64..4.0_f64);
            let size = rng.random_range(0.1_f32..0.4_f32);
            let grey: f32 = rng.random_range(0.5_f32..0.8_f32);
            particles.push(ParticleEffect::new(
                self.pos,
                [vx, vy, vz],
                [grey, grey, grey, 0.7],
                lifetime,
                size,
            ));
        }
        particles
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ExplosionEffect
// ─────────────────────────────────────────────────────────────────────────────

/// A burst explosion that spawns `n_particles` all at once from a center point.
#[derive(Debug, Clone)]
pub struct ExplosionEffect {
    /// Center of the explosion in world space.
    pub center: [f32; 3],
    /// Radius of the explosion sphere.
    pub radius: f32,
    /// Number of particles to spawn.
    pub n_particles: usize,
}

impl ExplosionEffect {
    /// Create a new explosion effect.
    pub fn new(center: [f32; 3], radius: f32, n_particles: usize) -> Self {
        Self {
            center,
            radius,
            n_particles,
        }
    }

    /// Spawn all particles for this explosion burst.
    ///
    /// Particles are launched radially outward with random velocity magnitudes.
    pub fn spawn_particles(&self) -> Vec<ParticleEffect> {
        use std::f32::consts::PI as PIf;
        let mut rng = rand::rng();
        let mut particles = Vec::with_capacity(self.n_particles);
        for _ in 0..self.n_particles {
            // Random direction on unit sphere
            let theta: f32 = rng.random_range(0.0_f32..PIf);
            let phi: f32 = rng.random_range(0.0_f32..(2.0 * PIf));
            let speed: f32 = rng.random_range(0.5_f32..self.radius);
            let vx = speed * theta.sin() * phi.cos();
            let vy = speed * theta.sin() * phi.sin();
            let vz = speed * theta.cos();
            let lifetime = rng.random_range(0.3_f64..1.0_f64);
            let size = rng.random_range(0.02_f32..0.1_f32);
            let r: f32 = rng.random_range(0.8_f32..1.0_f32);
            let g: f32 = rng.random_range(0.2_f32..0.6_f32);
            particles.push(ParticleEffect::new(
                self.center,
                [vx, vy, vz],
                [r, g, 0.0, 1.0],
                lifetime,
                size,
            ));
        }
        particles
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticlePool
// ─────────────────────────────────────────────────────────────────────────────

/// A bounded pool of particle effects.
///
/// When the pool is full, new particles are silently dropped.
#[derive(Debug, Clone)]
pub struct ParticlePool {
    /// Active and dead particles stored together.
    pub particles: Vec<ParticleEffect>,
    /// Maximum number of particles this pool may hold.
    pub max_particles: usize,
}

impl ParticlePool {
    /// Create a new empty pool with the given capacity.
    pub fn new(max_particles: usize) -> Self {
        Self {
            particles: Vec::with_capacity(max_particles),
            max_particles,
        }
    }

    /// Add a particle to the pool. If the pool is at capacity, the particle is dropped.
    pub fn add(&mut self, p: ParticleEffect) {
        if self.particles.len() < self.max_particles {
            self.particles.push(p);
        }
    }

    /// Advance all particles by `dt` and remove dead ones.
    pub fn update(&mut self, dt: f32) {
        for p in &mut self.particles {
            p.update(dt);
        }
        self.particles.retain(|p| p.is_alive());
    }

    /// Number of alive particles currently in the pool.
    pub fn alive_count(&self) -> usize {
        self.particles.iter().filter(|p| p.is_alive()).count()
    }

    /// Number of dead particles still in the pool (before the next update).
    pub fn dead_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.is_alive()).count()
    }

    /// Maximum capacity of this pool.
    pub fn capacity(&self) -> usize {
        self.max_particles
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Billboard helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the four corner positions of a billboard quad centered at `pos`.
///
/// The quad lies in the XY plane (no camera orientation applied) and has side
/// length `size`. Corners are returned in order: bottom-left, bottom-right,
/// top-right, top-left.
pub fn billboard_vertices(pos: [f32; 3], size: f32) -> [[f32; 3]; 4] {
    let half = size / 2.0;
    [
        [pos[0] - half, pos[1] - half, pos[2]], // bottom-left
        [pos[0] + half, pos[1] - half, pos[2]], // bottom-right
        [pos[0] + half, pos[1] + half, pos[2]], // top-right
        [pos[0] - half, pos[1] + half, pos[2]], // top-left
    ]
}

/// Sort a slice of particles back-to-front relative to `camera_pos` for correct alpha blending.
///
/// Uses the squared distance to the camera: farthest particles come first.
pub fn particle_sort_back_to_front(particles: &mut [ParticleEffect], camera_pos: [f32; 3]) {
    particles.sort_by(|a, b| {
        let dist_sq = |p: &ParticleEffect| {
            let dx = p.pos[0] - camera_pos[0];
            let dy = p.pos[1] - camera_pos[1];
            let dz = p.pos[2] - camera_pos[2];
            dx * dx + dy * dy + dz * dz
        };
        dist_sq(b)
            .partial_cmp(&dist_sq(a))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ParticleEffect ────────────────────────────────────────────────────────

    #[test]
    fn test_particle_starts_alive() {
        let p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1);
        assert!(p.is_alive());
    }

    #[test]
    fn test_particle_dies_after_lifetime() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 0.5, 0.1);
        p.update(0.6); // age > lifetime
        assert!(!p.is_alive());
    }

    #[test]
    fn test_particle_age_advances() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1);
        p.update(0.3);
        assert!((p.age - 0.3).abs() < 1e-7);
    }

    #[test]
    fn test_particle_position_advances() {
        let mut p = ParticleEffect::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0; 4], 5.0, 0.1);
        p.update(2.0);
        assert!((p.pos[0] - 2.0).abs() < 1e-6);
        assert!(p.pos[1].abs() < 1e-6);
    }

    #[test]
    fn test_particle_opacity_fades() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1);
        p.update(0.5); // halfway
        let alpha = p.alpha();
        assert!(alpha > 0.0 && alpha < 1.0, "opacity={alpha}");
    }

    #[test]
    fn test_particle_opacity_zero_at_end() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1);
        p.update(1.0);
        assert!(p.alpha() < 1e-6, "opacity at end={}", p.alpha());
    }

    #[test]
    fn test_particle_normalized_age_zero_at_start() {
        let p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1);
        assert!(p.normalized_age().abs() < 1e-6);
    }

    #[test]
    fn test_particle_normalized_age_half() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1);
        p.update(1.0);
        assert!((p.normalized_age() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_particle_normalized_age_clamps_to_one() {
        let mut p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1);
        p.update(5.0);
        assert!((p.normalized_age() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_particle_zero_lifetime() {
        let p = ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 0.0, 0.1);
        // normalized_age returns 1.0 for zero lifetime
        assert!((p.normalized_age() - 1.0).abs() < 1e-6);
    }

    // ── FireEmitter ───────────────────────────────────────────────────────────

    #[test]
    fn test_fire_emitter_emit_large_dt() {
        let emitter = FireEmitter::new([0.0; 3], 1.0);
        let particles = emitter.emit(1.0);
        assert!(
            !particles.is_empty(),
            "Fire should emit particles for dt=1.0"
        );
    }

    #[test]
    fn test_fire_emitter_emit_zero_dt() {
        let emitter = FireEmitter::new([0.0; 3], 1.0);
        let particles = emitter.emit(0.0);
        assert!(particles.is_empty(), "Zero dt should emit no particles");
    }

    #[test]
    fn test_fire_emitter_color_at_age_zero() {
        let emitter = FireEmitter::new([0.0; 3], 1.0);
        let c = emitter.color_at_age(0.0);
        // At t=0 the color should equal color_hot
        assert!((c[0] - emitter.color_hot[0]).abs() < 1e-6);
        assert!((c[1] - emitter.color_hot[1]).abs() < 1e-6);
    }

    #[test]
    fn test_fire_emitter_color_at_age_one() {
        let emitter = FireEmitter::new([0.0; 3], 1.0);
        let c = emitter.color_at_age(1.0);
        assert!((c[0] - emitter.color_cool[0]).abs() < 1e-6);
    }

    #[test]
    fn test_fire_emitter_color_midpoint() {
        let emitter = FireEmitter::new([0.0; 3], 1.0);
        let c = emitter.color_at_age(0.5);
        let expected = (emitter.color_hot[0] + emitter.color_cool[0]) / 2.0;
        assert!((c[0] - expected).abs() < 1e-6);
    }

    #[test]
    fn test_fire_emitter_particles_alive() {
        let emitter = FireEmitter::new([0.0; 3], 2.0);
        let particles = emitter.emit(2.0);
        assert!(!particles.is_empty());
        for p in &particles {
            assert!(p.is_alive(), "freshly emitted particles should be alive");
        }
    }

    // ── SmokeEmitter ──────────────────────────────────────────────────────────

    #[test]
    fn test_smoke_emitter_large_dt() {
        let emitter = SmokeEmitter::new([0.0; 3], 0.5);
        let particles = emitter.emit(2.0);
        assert!(!particles.is_empty(), "Smoke should emit for large dt");
    }

    #[test]
    fn test_smoke_emitter_zero_dt() {
        let emitter = SmokeEmitter::new([0.0; 3], 0.5);
        let particles = emitter.emit(0.0);
        assert!(particles.is_empty());
    }

    // ── ExplosionEffect ───────────────────────────────────────────────────────

    #[test]
    fn test_explosion_spawns_n_particles() {
        let expl = ExplosionEffect::new([0.0; 3], 5.0, 50);
        let particles = expl.spawn_particles();
        assert_eq!(particles.len(), 50);
    }

    #[test]
    fn test_explosion_zero_particles() {
        let expl = ExplosionEffect::new([0.0; 3], 1.0, 0);
        let particles = expl.spawn_particles();
        assert!(particles.is_empty());
    }

    #[test]
    fn test_explosion_particles_alive() {
        let expl = ExplosionEffect::new([0.0; 3], 3.0, 10);
        let particles = expl.spawn_particles();
        for p in &particles {
            assert!(p.is_alive());
        }
    }

    #[test]
    fn test_explosion_particles_start_at_center() {
        let center = [1.0_f32, 2.0, 3.0];
        let expl = ExplosionEffect::new(center, 2.0, 5);
        let particles = expl.spawn_particles();
        for p in &particles {
            assert!((p.pos[0] - center[0]).abs() < 1e-6);
            assert!((p.pos[1] - center[1]).abs() < 1e-6);
            assert!((p.pos[2] - center[2]).abs() < 1e-6);
        }
    }

    // ── ParticlePool ──────────────────────────────────────────────────────────

    #[test]
    fn test_pool_capacity() {
        let pool = ParticlePool::new(100);
        assert_eq!(pool.capacity(), 100);
    }

    #[test]
    fn test_pool_add_within_limit() {
        let mut pool = ParticlePool::new(5);
        for _ in 0..5 {
            pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1));
        }
        assert_eq!(pool.particles.len(), 5);
    }

    #[test]
    fn test_pool_add_drops_over_limit() {
        let mut pool = ParticlePool::new(3);
        for _ in 0..10 {
            pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 1.0, 0.1));
        }
        assert_eq!(pool.particles.len(), 3);
    }

    #[test]
    fn test_pool_update_ages_particles() {
        let mut pool = ParticlePool::new(10);
        pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1));
        pool.update(0.5);
        // After update, particle still alive but age advanced (still in pool)
        assert_eq!(pool.alive_count(), 1);
    }

    #[test]
    fn test_pool_update_removes_dead() {
        let mut pool = ParticlePool::new(10);
        pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 0.5, 0.1));
        pool.update(1.0); // kills the particle
        assert_eq!(pool.alive_count(), 0);
        assert_eq!(pool.particles.len(), 0);
    }

    #[test]
    fn test_pool_alive_dead_count() {
        let mut pool = ParticlePool::new(10);
        pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1));
        pool.add(ParticleEffect::new([0.0; 3], [0.0; 3], [1.0; 4], 2.0, 0.1));
        assert_eq!(pool.alive_count(), 2);
        assert_eq!(pool.dead_count(), 0);
    }

    // ── billboard_vertices ────────────────────────────────────────────────────

    #[test]
    fn test_billboard_four_corners() {
        let verts = billboard_vertices([0.0; 3], 2.0);
        assert_eq!(verts.len(), 4);
    }

    #[test]
    fn test_billboard_centered_at_pos() {
        let pos = [3.0_f32, 1.0, 0.0];
        let verts = billboard_vertices(pos, 2.0);
        let cx = verts.iter().map(|v| v[0]).sum::<f32>() / 4.0;
        let cy = verts.iter().map(|v| v[1]).sum::<f32>() / 4.0;
        assert!((cx - pos[0]).abs() < 1e-6);
        assert!((cy - pos[1]).abs() < 1e-6);
    }

    #[test]
    fn test_billboard_forms_square() {
        let verts = billboard_vertices([0.0; 3], 1.0);
        // width = top-right x - bottom-left x = 1.0
        let width = verts[1][0] - verts[0][0];
        let height = verts[2][1] - verts[1][1];
        assert!((width - 1.0).abs() < 1e-6);
        assert!((height - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_billboard_all_same_z() {
        let pos = [0.0, 0.0, 5.0_f32];
        let verts = billboard_vertices(pos, 1.0);
        for v in &verts {
            assert!((v[2] - 5.0).abs() < 1e-6);
        }
    }

    // ── particle_sort_back_to_front ───────────────────────────────────────────

    #[test]
    fn test_sort_back_to_front_ordering() {
        let cam = [0.0_f32; 3];
        let mut particles = vec![
            ParticleEffect::new([1.0, 0.0, 0.0], [0.0; 3], [1.0; 4], 1.0, 0.1), // dist=1
            ParticleEffect::new([3.0, 0.0, 0.0], [0.0; 3], [1.0; 4], 1.0, 0.1), // dist=3
            ParticleEffect::new([2.0, 0.0, 0.0], [0.0; 3], [1.0; 4], 1.0, 0.1), // dist=2
        ];
        particle_sort_back_to_front(&mut particles, cam);
        assert!((particles[0].pos[0] - 3.0).abs() < 1e-6); // farthest first
        assert!((particles[2].pos[0] - 1.0).abs() < 1e-6); // nearest last
    }

    #[test]
    fn test_sort_empty_slice() {
        let cam = [0.0_f32; 3];
        let mut particles: Vec<ParticleEffect> = vec![];
        particle_sort_back_to_front(&mut particles, cam);
        assert!(particles.is_empty());
    }

    #[test]
    fn test_sort_single_particle() {
        let cam = [0.0_f32; 3];
        let mut particles = vec![ParticleEffect::new(
            [1.0, 0.0, 0.0],
            [0.0; 3],
            [1.0; 4],
            1.0,
            0.1,
        )];
        particle_sort_back_to_front(&mut particles, cam);
        assert_eq!(particles.len(), 1);
    }
}
