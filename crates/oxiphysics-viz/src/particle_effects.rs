// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle effect rendering data structures.
//!
//! Provides single-particle effects, fire and smoke emitters, explosion bursts,
//! a pool-based particle manager, billboard quad generation, and back-to-front
//! sorting for correct alpha blending.
//!
//! The emitters draw their randomness from the thread-local generator by
//! default.  Every random entry point also has a `*_with_rng` variant taking
//! any [`rand::Rng`], so seeding e.g. a `StdRng` makes the output reproducible.

use rand::RngExt;

// ─────────────────────────────────────────────────────────────────────────────
// Sampling helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Draw a uniform sample from the interval spanned by `a` and `b` (in either
/// order) without ever panicking.
///
/// `rand` rejects empty or non-finite ranges, which is what a zero spread, a
/// zero rise speed or a tiny explosion radius would otherwise produce.  Here a
/// degenerate interval (`a == b`) yields that single value, an interval whose
/// span overflows `f32` yields its midpoint, and non-finite bounds yield `0.0`.
fn sample_between<R: rand::Rng + ?Sized>(rng: &mut R, a: f32, b: f32) -> f32 {
    if !a.is_finite() || !b.is_finite() {
        return 0.0;
    }
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    if hi > lo && (hi - lo).is_finite() {
        rng.random_range(lo..hi)
    } else {
        0.5 * lo + 0.5 * hi
    }
}

/// Largest number of particles a single emitter call produces.
///
/// It only guards against runaway allocations from absurd `rate * dt`
/// products (a stalled clock, a huge rate); real frames stay far below it.
pub const MAX_PARTICLES_PER_EMIT: usize = 1 << 20;

/// Advance a fractional emission `accumulator` by `rate * dt` and return how
/// many whole particles are due.
///
/// `rate` (particles per second) and `dt` (seconds) are treated as zero when
/// negative or NaN, so a bad frame time can neither emit a negative count nor
/// rewind the accumulator.  The fractional remainder stays in `accumulator`:
/// that carry is what lets a low rate (5 particles/s at 60 Hz is only 0.083
/// particles per frame) emit at all.  A non-finite total clears the accumulator
/// and emits nothing; a total of [`MAX_PARTICLES_PER_EMIT`] or more emits that
/// many and discards the excess.
fn take_whole_particles(accumulator: &mut f32, rate: f32, dt: f32) -> usize {
    let total = *accumulator + rate.max(0.0) * dt.max(0.0);
    if !total.is_finite() {
        *accumulator = 0.0;
        return 0;
    }
    let whole = total.floor();
    if whole >= MAX_PARTICLES_PER_EMIT as f32 {
        *accumulator = 0.0;
        return MAX_PARTICLES_PER_EMIT;
    }
    *accumulator = total - whole;
    whole as usize
}

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
    /// Fraction of a particle carried over between [`emit`](Self::emit) calls
    /// (always in `[0, 1)` after a call).
    pub emit_accumulator: f32,
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
            emit_accumulator: 0.0,
        }
    }

    /// Emit particles for a time step `dt` using the thread-local generator.
    ///
    /// See [`emit_with_rng`](Self::emit_with_rng) for the emission count rule.
    pub fn emit(&mut self, dt: f32) -> Vec<ParticleEffect> {
        let mut rng = rand::rng();
        self.emit_with_rng(dt, &mut rng)
    }

    /// Emit particles for a time step `dt`, drawing randomness from `rng`.
    ///
    /// `rate * dt` is added to [`emit_accumulator`](Self::emit_accumulator) and
    /// the whole part is emitted, so the long-run emission rate equals `rate`
    /// at any frame rate (at 60 Hz a rate of 20/s emits one particle every
    /// third frame, rather than none).  A negative or NaN `rate` or `dt`
    /// emits nothing, and one call never emits more than
    /// [`MAX_PARTICLES_PER_EMIT`] particles.
    pub fn emit_with_rng<R: rand::Rng + ?Sized>(
        &mut self,
        dt: f32,
        rng: &mut R,
    ) -> Vec<ParticleEffect> {
        let count = take_whole_particles(&mut self.emit_accumulator, self.rate, dt);
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
    /// Fraction of a particle carried over between [`emit`](Self::emit) calls
    /// (always in `[0, 1)` after a call).
    pub emit_accumulator: f32,
}

impl SmokeEmitter {
    /// Create a smoke emitter at `pos` with given `spread` radius.
    pub fn new(pos: [f32; 3], spread: f32) -> Self {
        Self {
            pos,
            rate: 5.0,
            spread,
            rise_speed: 0.3,
            emit_accumulator: 0.0,
        }
    }

    /// Emit smoke particles for a time step `dt` using the thread-local
    /// generator.
    ///
    /// See [`emit_with_rng`](Self::emit_with_rng) for the emission count rule.
    pub fn emit(&mut self, dt: f32) -> Vec<ParticleEffect> {
        let mut rng = rand::rng();
        self.emit_with_rng(dt, &mut rng)
    }

    /// Emit smoke particles for a time step `dt`, drawing randomness from `rng`.
    ///
    /// `rate * dt` is added to [`emit_accumulator`](Self::emit_accumulator) and
    /// the whole part is emitted, so the long-run emission rate equals `rate`
    /// at any frame rate.  A zero `spread` or `rise_speed` is valid and yields
    /// particles without horizontal or vertical drift respectively.
    pub fn emit_with_rng<R: rand::Rng + ?Sized>(
        &mut self,
        dt: f32,
        rng: &mut R,
    ) -> Vec<ParticleEffect> {
        let count = take_whole_particles(&mut self.emit_accumulator, self.rate, dt);
        let mut particles = Vec::with_capacity(count);
        for _ in 0..count {
            let vx: f32 = sample_between(rng, -self.spread, self.spread);
            let vy: f32 = sample_between(rng, self.rise_speed * 0.5, self.rise_speed * 1.5);
            let vz: f32 = sample_between(rng, -self.spread, self.spread);
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

    /// Spawn all particles for this explosion burst using the thread-local
    /// generator.
    ///
    /// See [`spawn_particles_with_rng`](Self::spawn_particles_with_rng).
    pub fn spawn_particles(&self) -> Vec<ParticleEffect> {
        let mut rng = rand::rng();
        self.spawn_particles_with_rng(&mut rng)
    }

    /// Spawn all particles for this explosion burst, drawing randomness from
    /// `rng`.
    ///
    /// Particles are launched radially outward in directions distributed
    /// uniformly over the sphere (`cos(theta)` and the azimuth are sampled
    /// uniformly) with speeds in `[0.5, radius)`.  A `radius` of at most `0.5`
    /// leaves no such interval, so every particle then moves at `radius`
    /// (clamped to be non-negative).
    pub fn spawn_particles_with_rng<R: rand::Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> Vec<ParticleEffect> {
        use std::f32::consts::PI as PIf;
        let max_speed = if self.radius.is_finite() {
            self.radius.max(0.0)
        } else {
            0.0
        };
        let min_speed = 0.5_f32.min(max_speed);
        let mut particles = Vec::with_capacity(self.n_particles);
        for _ in 0..self.n_particles {
            // Uniform direction on the unit sphere: cos(theta) is uniform in
            // [-1, 1), the azimuth uniform in [0, 2 pi).
            let cos_theta = sample_between(rng, -1.0, 1.0);
            let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
            let phi = sample_between(rng, 0.0, 2.0 * PIf);
            let speed = sample_between(rng, min_speed, max_speed);
            let vx = speed * sin_theta * phi.cos();
            let vy = speed * sin_theta * phi.sin();
            let vz = speed * cos_theta;
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
    use rand::SeedableRng;
    use rand::rngs::StdRng;

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
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
        let particles = emitter.emit(1.0);
        assert!(
            !particles.is_empty(),
            "Fire should emit particles for dt=1.0"
        );
    }

    #[test]
    fn test_fire_emitter_emit_zero_dt() {
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
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
        let mut emitter = FireEmitter::new([0.0; 3], 2.0);
        let particles = emitter.emit(2.0);
        assert!(!particles.is_empty());
        for p in &particles {
            assert!(p.is_alive(), "freshly emitted particles should be alive");
        }
    }

    #[test]
    fn test_fire_emitter_emits_at_real_frame_rate() {
        // 20 particles/s at 60 Hz is only 0.33 particles per frame; the
        // remainder must be carried so that one second yields ~20 particles.
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
        let mut total = 0;
        for _ in 0..60 {
            total += emitter.emit(1.0 / 60.0).len();
        }
        assert!(
            (19..=21).contains(&total),
            "60 frames of 1/60 s at 20/s must emit 19..=21 particles, got {total}"
        );
    }

    #[test]
    fn test_fire_emitter_carries_fractional_remainder() {
        // 2 particles/s with dt = 0.25 s (exact in binary) adds 0.5 per call.
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
        emitter.rate = 2.0;
        let counts: Vec<usize> = (0..4).map(|_| emitter.emit(0.25).len()).collect();
        assert_eq!(counts, vec![0, 1, 0, 1]);
        assert!(emitter.emit_accumulator.abs() < 1e-6);

        // 20/s with a 0.125 s step: 2.5 particles -> 2 now, 0.5 carried over;
        // a further 0.0625 s adds 1.25 -> 1 more, 0.75 carried over.
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
        assert_eq!(emitter.emit(0.125).len(), 2);
        assert!((emitter.emit_accumulator - 0.5).abs() < 1e-6);
        assert_eq!(emitter.emit(0.0625).len(), 1);
        assert!((emitter.emit_accumulator - 0.75).abs() < 1e-6);
    }

    #[test]
    fn test_fire_emitter_ignores_negative_and_nan_dt() {
        let mut emitter = FireEmitter::new([0.0; 3], 1.0);
        assert!(emitter.emit(-1.0).is_empty());
        assert!(emitter.emit(f32::NAN).is_empty());
        assert_eq!(emitter.emit_accumulator, 0.0);
        // A negative rate also emits nothing.
        emitter.rate = -5.0;
        assert!(emitter.emit(10.0).is_empty());
        // An infinite step clears the accumulator instead of flooding.
        emitter.rate = 20.0;
        assert!(emitter.emit(f32::INFINITY).is_empty());
        assert_eq!(emitter.emit_accumulator, 0.0);
    }

    #[test]
    fn test_take_whole_particles_caps_runaway_counts() {
        // A finite but absurd product must not request an impossible allocation.
        let mut acc = 0.25_f32;
        assert_eq!(
            take_whole_particles(&mut acc, 1.0e30, 1.0),
            MAX_PARTICLES_PER_EMIT
        );
        assert_eq!(acc, 0.0, "the excess is discarded, not carried over");
        // Realistic counts are untouched.
        let mut acc = 0.0_f32;
        assert_eq!(take_whole_particles(&mut acc, 1000.0, 2.5), 2500);
        assert_eq!(acc, 0.0);
    }

    #[test]
    fn test_fire_emitter_with_rng_is_reproducible() {
        let mut a = FireEmitter::new([1.0, 2.0, 3.0], 1.5);
        let mut b = a.clone();
        let mut rng_a = StdRng::seed_from_u64(7);
        let mut rng_b = StdRng::seed_from_u64(7);
        let pa = a.emit_with_rng(1.0, &mut rng_a);
        let pb = b.emit_with_rng(1.0, &mut rng_b);
        assert_eq!(pa.len(), 30);
        assert_eq!(pa.len(), pb.len());
        for (x, y) in pa.iter().zip(pb.iter()) {
            assert_eq!(x.vel, y.vel);
            assert_eq!(x.size, y.size);
            assert_eq!(x.lifetime, y.lifetime);
        }
    }

    // ── SmokeEmitter ──────────────────────────────────────────────────────────

    #[test]
    fn test_smoke_emitter_emits_at_real_frame_rate() {
        // 5 particles/s at 60 Hz is 0.083 per frame and used to emit nothing.
        let mut emitter = SmokeEmitter::new([0.0; 3], 0.5);
        let mut total = 0;
        for _ in 0..60 {
            total += emitter.emit(1.0 / 60.0).len();
        }
        assert!(
            (4..=6).contains(&total),
            "60 frames of 1/60 s at 5/s must emit 4..=6 particles, got {total}"
        );
    }

    #[test]
    fn test_smoke_emitter_zero_spread_and_rise_speed_do_not_panic() {
        let mut emitter = SmokeEmitter::new([0.0; 3], 0.0);
        let particles = emitter.emit(2.0);
        assert_eq!(particles.len(), 10);
        for p in &particles {
            assert_eq!(p.vel[0], 0.0, "zero spread gives no horizontal drift");
            assert_eq!(p.vel[2], 0.0, "zero spread gives no horizontal drift");
            assert!(p.vel[1] > 0.0, "default rise speed still lifts the smoke");
        }
        emitter.rise_speed = 0.0;
        let particles = emitter.emit(2.0);
        assert_eq!(particles.len(), 10);
        for p in &particles {
            assert_eq!(p.vel[1], 0.0, "zero rise speed gives no vertical drift");
        }
    }

    #[test]
    fn test_smoke_emitter_negative_and_nan_parameters_do_not_panic() {
        let mut emitter = SmokeEmitter::new([0.0; 3], -0.5);
        emitter.rise_speed = -1.0;
        let particles = emitter.emit(2.0);
        assert_eq!(particles.len(), 10);
        for p in &particles {
            assert!(p.vel[0].abs() <= 0.5 && p.vel[2].abs() <= 0.5);
            assert!(
                (-1.5..=-0.5).contains(&p.vel[1]),
                "negative rise speed sinks the smoke: {}",
                p.vel[1]
            );
        }
        emitter.spread = f32::NAN;
        emitter.rise_speed = f32::NAN;
        for p in emitter.emit(2.0) {
            assert!(p.vel.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn test_smoke_emitter_with_rng_is_reproducible() {
        let mut a = SmokeEmitter::new([0.0; 3], 0.5);
        let mut b = a.clone();
        let pa = a.emit_with_rng(2.0, &mut StdRng::seed_from_u64(99));
        let pb = b.emit_with_rng(2.0, &mut StdRng::seed_from_u64(99));
        assert_eq!(pa.len(), 10);
        assert_eq!(
            pa.iter().map(|p| p.vel).collect::<Vec<_>>(),
            pb.iter().map(|p| p.vel).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_smoke_emitter_large_dt() {
        let mut emitter = SmokeEmitter::new([0.0; 3], 0.5);
        let particles = emitter.emit(2.0);
        assert!(!particles.is_empty(), "Smoke should emit for large dt");
    }

    #[test]
    fn test_smoke_emitter_zero_dt() {
        let mut emitter = SmokeEmitter::new([0.0; 3], 0.5);
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

    #[test]
    fn test_explosion_small_radius_does_not_panic() {
        // radius 0.3 < 0.5 used to produce the empty range 0.5..0.3.
        let expl = ExplosionEffect::new([0.0; 3], 0.3, 25);
        let particles = expl.spawn_particles();
        assert_eq!(particles.len(), 25);
        for p in &particles {
            let speed = (p.vel[0].powi(2) + p.vel[1].powi(2) + p.vel[2].powi(2)).sqrt();
            assert!(
                (speed - 0.3).abs() < 1e-5,
                "speed must fall back to the radius, got {speed}"
            );
        }
    }

    #[test]
    fn test_explosion_degenerate_radii_do_not_panic() {
        for radius in [0.0_f32, 0.5, -3.0, f32::NAN, f32::INFINITY] {
            let expl = ExplosionEffect::new([0.0; 3], radius, 8);
            let particles = expl.spawn_particles();
            assert_eq!(particles.len(), 8);
            for p in &particles {
                assert!(p.vel.iter().all(|v| v.is_finite()), "radius {radius}");
            }
        }
    }

    #[test]
    fn test_explosion_speed_within_radius() {
        let expl = ExplosionEffect::new([0.0; 3], 5.0, 500);
        let particles = expl.spawn_particles_with_rng(&mut StdRng::seed_from_u64(3));
        for p in &particles {
            let speed = (p.vel[0].powi(2) + p.vel[1].powi(2) + p.vel[2].powi(2)).sqrt();
            assert!(
                (0.5 - 1e-4..=5.0 + 1e-4).contains(&speed),
                "speed {speed} outside [0.5, 5)"
            );
        }
    }

    #[test]
    fn test_explosion_directions_are_uniform_on_sphere() {
        // For directions uniform on the sphere cos(theta) is uniform on
        // [-1, 1): E[cos] = 0 and E[cos^2] = 1/3.  Sampling theta uniformly
        // instead would give E[cos^2] = 1/2 (clustering at the poles).
        let expl = ExplosionEffect::new([0.0; 3], 5.0, 20_000);
        let particles = expl.spawn_particles_with_rng(&mut StdRng::seed_from_u64(11));
        let mut mean_cos = 0.0_f64;
        let mut mean_cos2 = 0.0_f64;
        for p in &particles {
            let speed = (p.vel[0].powi(2) + p.vel[1].powi(2) + p.vel[2].powi(2)).sqrt();
            let cos_theta = f64::from(p.vel[2] / speed);
            mean_cos += cos_theta;
            mean_cos2 += cos_theta * cos_theta;
        }
        let n = particles.len() as f64;
        mean_cos /= n;
        mean_cos2 /= n;
        assert!(mean_cos.abs() < 0.03, "mean cos(theta) = {mean_cos}");
        assert!(
            (mean_cos2 - 1.0 / 3.0).abs() < 0.02,
            "mean cos^2(theta) = {mean_cos2}, expected 1/3"
        );
    }

    #[test]
    fn test_explosion_with_rng_is_reproducible() {
        let expl = ExplosionEffect::new([1.0, 0.0, -1.0], 4.0, 40);
        let a = expl.spawn_particles_with_rng(&mut StdRng::seed_from_u64(5));
        let b = expl.spawn_particles_with_rng(&mut StdRng::seed_from_u64(5));
        let c = expl.spawn_particles_with_rng(&mut StdRng::seed_from_u64(6));
        let vels = |ps: &[ParticleEffect]| ps.iter().map(|p| p.vel).collect::<Vec<_>>();
        assert_eq!(vels(&a), vels(&b));
        assert_ne!(vels(&a), vels(&c));
    }

    // ── sample_between ────────────────────────────────────────────────────────

    #[test]
    fn test_sample_between_handles_degenerate_intervals() {
        let mut rng = StdRng::seed_from_u64(1);
        assert_eq!(sample_between(&mut rng, 2.0, 2.0), 2.0);
        assert_eq!(sample_between(&mut rng, 0.0, 0.0), 0.0);
        // Non-finite bounds and an overflowing span never panic.
        assert_eq!(sample_between(&mut rng, f32::NAN, 1.0), 0.0);
        assert_eq!(sample_between(&mut rng, 0.0, f32::INFINITY), 0.0);
        assert_eq!(sample_between(&mut rng, -3.0e38, 3.0e38), 0.0);
        // Reversed bounds sample the same interval.
        for _ in 0..200 {
            let v = sample_between(&mut rng, 5.0, 2.0);
            assert!((2.0..=5.0).contains(&v), "{v}");
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
