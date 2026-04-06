// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly particle system.
//!
//! Provides a complete CPU-side particle simulation designed for WASM export.
//! Uses plain `f64`/`f32` arrays — no nalgebra or heavy dependencies.

#![allow(dead_code)]
#![allow(missing_docs)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Minimal LCG RNG (no external crates)
// ---------------------------------------------------------------------------

pub(crate) struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    pub(crate) fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x853c_49e6_748f_ea9b,
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        let bits = self.next_u64() >> 11;
        bits as f64 / (1u64 << 53) as f64
    }

    fn range_f64(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
    }

    fn unit_sphere(&mut self) -> [f64; 3] {
        loop {
            let x = self.range_f64(-1.0, 1.0);
            let y = self.range_f64(-1.0, 1.0);
            let z = self.range_f64(-1.0, 1.0);
            let r2 = x * x + y * y + z * z;
            if r2 <= 1.0 && r2 > 1e-15 {
                let r = r2.sqrt();
                return [x / r, y / r, z / r];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// WasmParticleConfig
// ---------------------------------------------------------------------------

/// Configuration for the particle system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleConfig {
    /// Maximum number of particles alive simultaneously.
    pub max_particles: usize,
    /// Gravity vector `[gx, gy, gz]`.
    pub gravity: [f64; 3],
    /// Velocity damping factor per second (0 = no damping, 1 = full stop).
    pub damping: f64,
    /// Radius used for inter-particle collision checks.
    pub collision_radius: f64,
    /// Minimum and maximum particle lifetime in seconds.
    pub lifetime_range: [f64; 2],
}

impl WasmParticleConfig {
    /// Default config with Earth gravity.
    pub fn default_earth() -> Self {
        Self {
            max_particles: 1000,
            gravity: [0.0, -9.81, 0.0],
            damping: 0.01,
            collision_radius: 0.05,
            lifetime_range: [1.0, 5.0],
        }
    }
}

impl Default for WasmParticleConfig {
    fn default() -> Self {
        Self::default_earth()
    }
}

// ---------------------------------------------------------------------------
// WasmParticle
// ---------------------------------------------------------------------------

/// A single particle in the system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticle {
    /// World-space position.
    pub position: [f64; 3],
    /// Velocity in m/s.
    pub velocity: [f64; 3],
    /// Current age in seconds.
    pub age: f64,
    /// Lifetime in seconds (particle dies when age >= lifetime).
    pub lifetime: f64,
    /// RGBA color `[r, g, b, a]` in \[0, 1\].
    pub color: [f32; 4],
    /// Billboard size in world units.
    pub size: f64,
    /// Whether the particle is alive.
    pub alive: bool,
}

impl WasmParticle {
    /// Create a new particle at the given position.
    pub fn new(position: [f64; 3], velocity: [f64; 3], lifetime: f64) -> Self {
        Self {
            position,
            velocity,
            age: 0.0,
            lifetime,
            color: [1.0, 1.0, 1.0, 1.0],
            size: 0.1,
            alive: true,
        }
    }

    /// Normalised age in \[0, 1\].
    pub fn normalized_age(&self) -> f64 {
        if self.lifetime <= 0.0 {
            1.0
        } else {
            (self.age / self.lifetime).clamp(0.0, 1.0)
        }
    }

    /// True when particle has exceeded its lifetime.
    pub fn is_dead(&self) -> bool {
        self.age >= self.lifetime
    }
}

// ---------------------------------------------------------------------------
// WasmEmitterShape
// ---------------------------------------------------------------------------

/// Shape from which particles are emitted.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmEmitterShape {
    /// Emit from a single point.
    Point,
    /// Emit from within a sphere of the given radius.
    Sphere(f64),
    /// Emit from within an axis-aligned box with the given half-extents.
    Box([f64; 3]),
    /// Emit from a cone: `angle` is half-angle in radians, `radius` is base disc radius.
    Cone { angle: f64, radius: f64 },
    /// Emit from a disc of the given radius in the XZ plane.
    Disc(f64),
}

impl WasmEmitterShape {
    /// Sample a spawn position from this shape.
    fn sample_position(&self, rng: &mut SimpleRng, origin: [f64; 3]) -> [f64; 3] {
        match self {
            WasmEmitterShape::Point => origin,
            WasmEmitterShape::Sphere(r) => {
                let dir = rng.unit_sphere();
                let t = rng.range_f64(0.0, 1.0).cbrt();
                add3(origin, scale3(dir, t * r))
            }
            WasmEmitterShape::Box(h) => {
                let dx = rng.range_f64(-h[0], h[0]);
                let dy = rng.range_f64(-h[1], h[1]);
                let dz = rng.range_f64(-h[2], h[2]);
                add3(origin, [dx, dy, dz])
            }
            WasmEmitterShape::Cone { angle, radius } => {
                let t = rng.range_f64(0.0, 1.0).sqrt() * radius;
                let theta = rng.range_f64(0.0, std::f64::consts::TAU);
                let h = t / angle.tan().max(1e-6);
                add3(origin, [t * theta.cos(), h, t * theta.sin()])
            }
            WasmEmitterShape::Disc(r) => {
                let t = rng.range_f64(0.0, 1.0).sqrt() * r;
                let theta = rng.range_f64(0.0, std::f64::consts::TAU);
                add3(origin, [t * theta.cos(), 0.0, t * theta.sin()])
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WasmEmitter
// ---------------------------------------------------------------------------

/// Controls where and how particles are spawned.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmEmitter {
    /// Shape from which particles emerge.
    pub shape: WasmEmitterShape,
    /// Particles per second emission rate.
    pub rate: f64,
    /// Number of particles to emit in a single burst (0 = continuous only).
    pub burst_count: usize,
    /// Min/max initial speed range.
    pub velocity_range: [f64; 2],
    /// Origin of the emitter in world space.
    pub origin: [f64; 3],
    /// Accumulated fractional particles not yet emitted.
    accumulator: f64,
}

impl WasmEmitter {
    /// Create a new emitter.
    pub fn new(shape: WasmEmitterShape, rate: f64, origin: [f64; 3]) -> Self {
        Self {
            shape,
            rate,
            burst_count: 0,
            velocity_range: [1.0, 5.0],
            origin,
            accumulator: 0.0,
        }
    }

    /// Advance the emitter and return the positions + velocities of new particles.
    pub(crate) fn emit(&mut self, dt: f64, rng: &mut SimpleRng) -> Vec<([f64; 3], [f64; 3])> {
        self.accumulator += self.rate * dt;
        let count = self.accumulator as usize + self.burst_count;
        self.accumulator -= (self.accumulator as usize) as f64;
        self.burst_count = 0; // reset burst

        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            let pos = self.shape.sample_position(rng, self.origin);
            let dir = rng.unit_sphere();
            let speed = rng.range_f64(self.velocity_range[0], self.velocity_range[1]);
            out.push((pos, scale3(dir, speed)));
        }
        out
    }
}

// ---------------------------------------------------------------------------
// WasmParticleForce
// ---------------------------------------------------------------------------

/// External forces that can be applied to particles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmParticleForce {
    /// Gravitational acceleration.
    Gravity([f64; 3]),
    /// Constant wind velocity added to particle velocity.
    Wind([f64; 3]),
    /// Noise-based turbulence.
    Turbulence { freq: f64, amplitude: f64 },
    /// Point attractor pulling particles toward a position.
    Attractor { pos: [f64; 3], strength: f64 },
}

impl WasmParticleForce {
    /// Compute the acceleration contribution for a particle at `pos`.
    pub fn acceleration(&self, pos: [f64; 3], time: f64) -> [f64; 3] {
        match self {
            WasmParticleForce::Gravity(g) => *g,
            WasmParticleForce::Wind(w) => *w,
            WasmParticleForce::Turbulence { freq, amplitude } => {
                // Simple pseudo-turbulence using sin/cos
                let t = time * freq;
                let ax = amplitude * (t * 1.3 + pos[0]).sin();
                let ay = amplitude * (t * 0.7 + pos[1]).cos();
                let az = amplitude * (t * 1.1 + pos[2]).sin();
                [ax, ay, az]
            }
            WasmParticleForce::Attractor {
                pos: apos,
                strength,
            } => {
                let d = sub3(*apos, pos);
                let dist = len3(d) + 1e-6;
                scale3(normalize3(d), strength / (dist * dist))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WasmParticleCollider
// ---------------------------------------------------------------------------

/// Collision shape for particles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmColliderShape {
    /// Sphere collider at a position with a radius.
    Sphere { center: [f64; 3], radius: f64 },
    /// Infinite plane with normal and offset.
    Plane { normal: [f64; 3], offset: f64 },
}

/// Particle collider with restitution and friction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleCollider {
    /// The collision shape.
    pub shape: WasmColliderShape,
    /// Coefficient of restitution (0 = inelastic, 1 = elastic).
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
}

impl WasmParticleCollider {
    /// Create a plane collider.
    pub fn plane(normal: [f64; 3], offset: f64, restitution: f64) -> Self {
        Self {
            shape: WasmColliderShape::Plane {
                normal: normalize3(normal),
                offset,
            },
            restitution,
            friction: 0.1,
        }
    }

    /// Create a sphere collider.
    pub fn sphere(center: [f64; 3], radius: f64, restitution: f64) -> Self {
        Self {
            shape: WasmColliderShape::Sphere { center, radius },
            restitution,
            friction: 0.1,
        }
    }

    /// Resolve collision for a particle: returns the corrected (position, velocity).
    pub fn resolve(&self, pos: [f64; 3], vel: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        match &self.shape {
            WasmColliderShape::Plane { normal, offset } => {
                let n = *normal;
                let dist = dot3(pos, n) - offset;
                if dist < 0.0 {
                    let new_pos = add3(pos, scale3(n, -dist + 1e-5));
                    let vn = dot3(vel, n);
                    if vn < 0.0 {
                        let normal_v = scale3(n, vn);
                        let tangent_v = sub3(vel, normal_v);
                        let new_vel = add3(
                            scale3(n, -vn * self.restitution),
                            scale3(tangent_v, 1.0 - self.friction),
                        );
                        (new_pos, new_vel)
                    } else {
                        (new_pos, vel)
                    }
                } else {
                    (pos, vel)
                }
            }
            WasmColliderShape::Sphere { center, radius } => {
                let d = sub3(pos, *center);
                let dist = len3(d);
                if dist < *radius {
                    let n = if dist > 1e-15 {
                        scale3(d, 1.0 / dist)
                    } else {
                        [0.0, 1.0, 0.0]
                    };
                    let new_pos = add3(*center, scale3(n, radius + 1e-5));
                    let vn = dot3(vel, n);
                    if vn < 0.0 {
                        let normal_v = scale3(n, vn);
                        let tangent_v = sub3(vel, normal_v);
                        let new_vel = add3(
                            scale3(n, -vn * self.restitution),
                            scale3(tangent_v, 1.0 - self.friction),
                        );
                        (new_pos, new_vel)
                    } else {
                        (new_pos, vel)
                    }
                } else {
                    (pos, vel)
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WasmTrailSystem
// ---------------------------------------------------------------------------

/// Trail point: position + time stamp.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmTrailPoint {
    /// World position.
    pub position: [f64; 3],
    /// Time the trail point was recorded.
    pub time: f64,
}

/// Per-particle trail data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmTrail {
    /// Trail positions in order from oldest to newest.
    pub points: Vec<WasmTrailPoint>,
    /// Maximum trail length (points).
    pub max_length: usize,
    /// Whether trail fades toward the end.
    pub fading: bool,
}

impl WasmTrail {
    /// Create a new trail.
    pub fn new(max_length: usize, fading: bool) -> Self {
        Self {
            points: Vec::with_capacity(max_length),
            max_length,
            fading,
        }
    }

    /// Add a new trail point at the current time.
    pub fn push(&mut self, pos: [f64; 3], time: f64) {
        self.points.push(WasmTrailPoint {
            position: pos,
            time,
        });
        while self.points.len() > self.max_length {
            self.points.remove(0);
        }
    }

    /// Compute the width at normalised position `t` in \[0,1\] (0 = oldest).
    pub fn width_at(&self, t: f64, base_width: f64) -> f64 {
        if self.fading {
            base_width * t
        } else {
            base_width
        }
    }
}

/// Manages trails for multiple particles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmTrailSystem {
    /// Per-particle trails (indexed by particle ID).
    pub trails: Vec<WasmTrail>,
    /// Maximum trail length.
    pub max_length: usize,
    /// Whether trails fade.
    pub fading: bool,
}

impl WasmTrailSystem {
    /// Create a trail system.
    pub fn new(max_particles: usize, max_length: usize, fading: bool) -> Self {
        let trails = (0..max_particles)
            .map(|_| WasmTrail::new(max_length, fading))
            .collect();
        Self {
            trails,
            max_length,
            fading,
        }
    }

    /// Update trail for particle `id`.
    pub fn update(&mut self, id: usize, pos: [f64; 3], time: f64) {
        if let Some(trail) = self.trails.get_mut(id) {
            trail.push(pos, time);
        }
    }

    /// Clear trail for particle `id`.
    pub fn clear(&mut self, id: usize) {
        if let Some(trail) = self.trails.get_mut(id) {
            trail.points.clear();
        }
    }
}

// ---------------------------------------------------------------------------
// WasmParticleRenderer
// ---------------------------------------------------------------------------

/// Color/alpha curve key-frame.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmColorKey {
    /// Time in normalised \[0, 1\] particle age.
    pub time: f32,
    /// RGBA color.
    pub color: [f32; 4],
}

/// Rendering parameters for particles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleRenderer {
    /// Billboard size over lifetime (time, size) key-frames.
    pub size_curve: Vec<[f32; 2]>,
    /// Color over lifetime key-frames.
    pub color_curve: Vec<WasmColorKey>,
    /// Alpha over lifetime: (time, alpha) key-frames.
    pub alpha_curve: Vec<[f32; 2]>,
}

impl WasmParticleRenderer {
    /// Create a default renderer (constant white, size 0.1).
    pub fn default_white() -> Self {
        Self {
            size_curve: vec![[0.0, 0.1], [1.0, 0.0]],
            color_curve: vec![
                WasmColorKey {
                    time: 0.0,
                    color: [1.0, 1.0, 1.0, 1.0],
                },
                WasmColorKey {
                    time: 1.0,
                    color: [1.0, 1.0, 1.0, 0.0],
                },
            ],
            alpha_curve: vec![[0.0, 1.0], [1.0, 0.0]],
        }
    }

    /// Sample size at normalised age `t`.
    pub fn size_at(&self, t: f32) -> f32 {
        sample_curve_f32(&self.size_curve, t)
    }

    /// Sample alpha at normalised age `t`.
    pub fn alpha_at(&self, t: f32) -> f32 {
        sample_curve_f32(&self.alpha_curve, t)
    }

    /// Sample color at normalised age `t` (linear interpolation).
    pub fn color_at(&self, t: f32) -> [f32; 4] {
        if self.color_curve.is_empty() {
            return [1.0, 1.0, 1.0, 1.0];
        }
        let n = self.color_curve.len();
        if n == 1 || t <= self.color_curve[0].time {
            return self.color_curve[0].color;
        }
        if t >= self.color_curve[n - 1].time {
            return self.color_curve[n - 1].color;
        }
        for i in 0..n - 1 {
            let a = &self.color_curve[i];
            let b = &self.color_curve[i + 1];
            if t >= a.time && t <= b.time {
                let range = b.time - a.time;
                let f = if range > 1e-6 {
                    (t - a.time) / range
                } else {
                    0.0
                };
                return [
                    a.color[0] + (b.color[0] - a.color[0]) * f,
                    a.color[1] + (b.color[1] - a.color[1]) * f,
                    a.color[2] + (b.color[2] - a.color[2]) * f,
                    a.color[3] + (b.color[3] - a.color[3]) * f,
                ];
            }
        }
        self.color_curve[n - 1].color
    }
}

fn sample_curve_f32(curve: &[[f32; 2]], t: f32) -> f32 {
    if curve.is_empty() {
        return 0.0;
    }
    let n = curve.len();
    if n == 1 || t <= curve[0][0] {
        return curve[0][1];
    }
    if t >= curve[n - 1][0] {
        return curve[n - 1][1];
    }
    for i in 0..n - 1 {
        let ta = curve[i][0];
        let tb = curve[i + 1][0];
        if t >= ta && t <= tb {
            let range = tb - ta;
            let f = if range > 1e-6 { (t - ta) / range } else { 0.0 };
            return curve[i][1] + (curve[i + 1][1] - curve[i][1]) * f;
        }
    }
    curve[n - 1][1]
}

// ---------------------------------------------------------------------------
// WasmFireParticle
// ---------------------------------------------------------------------------

/// Fire-specific particle data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmFireParticle {
    /// Base particle data.
    pub particle: WasmParticle,
    /// Temperature in arbitrary units (high = yellow/white flame, low = red/smoke).
    pub temperature: f64,
    /// Heat-rise buoyancy factor.
    pub heat_rise: f64,
    /// Transition to smoke: 0 = pure fire, 1 = pure smoke.
    pub smoke_factor: f64,
}

impl WasmFireParticle {
    /// Create a new fire particle.
    pub fn new(position: [f64; 3], temperature: f64, lifetime: f64) -> Self {
        let vel = [0.0, heat_rise_vel(temperature), 0.0];
        let mut p = WasmParticle::new(position, vel, lifetime);
        p.color = fire_color(temperature);
        Self {
            particle: p,
            temperature,
            heat_rise: heat_rise_vel(temperature),
            smoke_factor: 0.0,
        }
    }

    /// Step this fire particle by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.particle.age += dt;
        let t = self.particle.normalized_age();
        // Temperature cools over time
        self.temperature *= 1.0 - t * 0.8;
        self.smoke_factor = t * t;
        // Update rise based on temperature
        self.heat_rise = heat_rise_vel(self.temperature);
        self.particle.velocity[1] = self.heat_rise;
        // Integrate position
        for i in 0..3 {
            self.particle.position[i] += self.particle.velocity[i] * dt;
        }
        self.particle.color = lerp_color(
            fire_color(self.temperature),
            [0.3, 0.3, 0.3, 0.5],
            self.smoke_factor as f32,
        );
        if self.particle.is_dead() {
            self.particle.alive = false;
        }
    }
}

fn heat_rise_vel(temperature: f64) -> f64 {
    (temperature * 0.05).max(0.0)
}

fn fire_color(temp: f64) -> [f32; 4] {
    // Simplified black-body: low = red, mid = orange, high = yellow/white
    let t = (temp / 1000.0).clamp(0.0, 1.0) as f32;
    [1.0, t * 0.7, t * t * 0.3, 1.0 - t * 0.2]
}

fn lerp_color(a: [f32; 4], b: [f32; 4], t: f32) -> [f32; 4] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
        a[3] + (b[3] - a[3]) * t,
    ]
}

// ---------------------------------------------------------------------------
// WasmParticleSystem
// ---------------------------------------------------------------------------

/// Complete particle system managing all particles and forces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmParticleSystem {
    /// Configuration.
    pub config: WasmParticleConfig,
    /// Active particles (some may be dead but not yet recycled).
    pub particles: Vec<WasmParticle>,
    /// Registered forces.
    pub forces: Vec<WasmParticleForce>,
    /// Registered colliders.
    pub colliders: Vec<WasmParticleCollider>,
    /// Current simulation time.
    pub time: f64,
    /// RNG seed (stored for determinism; actual RNG is local per step).
    rng_seed: u64,
}

impl WasmParticleSystem {
    /// Create a new particle system.
    pub fn new(config: WasmParticleConfig) -> Self {
        let max = config.max_particles;
        Self {
            config,
            particles: Vec::with_capacity(max),
            forces: vec![],
            colliders: vec![],
            time: 0.0,
            rng_seed: 42,
        }
    }

    /// Create with Earth gravity pre-configured.
    pub fn earth() -> Self {
        let mut sys = Self::new(WasmParticleConfig::default_earth());
        sys.forces
            .push(WasmParticleForce::Gravity([0.0, -9.81, 0.0]));
        sys
    }

    /// Spawn `n` particles at the origin with random velocities.
    pub fn spawn(&mut self, n: usize) {
        let mut rng = SimpleRng::new(self.rng_seed);
        self.rng_seed = rng.next_u64();
        let alive = self.particles.iter().filter(|p| p.alive).count();
        let to_spawn = n.min(self.config.max_particles.saturating_sub(alive));
        for _ in 0..to_spawn {
            let dir = rng.unit_sphere();
            let speed = rng.range_f64(1.0, 5.0);
            let lt = rng.range_f64(self.config.lifetime_range[0], self.config.lifetime_range[1]);
            let mut p = WasmParticle::new([0.0; 3], scale3(dir, speed), lt);
            p.size = rng.range_f64(0.05, 0.2);
            // Recycle a dead slot if available
            let slot = self.particles.iter().position(|p| !p.alive);
            if let Some(i) = slot {
                self.particles[i] = p;
            } else {
                self.particles.push(p);
            }
        }
    }

    /// Spawn particles from a given emitter.
    pub fn spawn_from_emitter(&mut self, emitter: &mut WasmEmitter, dt: f64) {
        let mut rng = SimpleRng::new(self.rng_seed);
        self.rng_seed = rng.next_u64();
        let spawns = emitter.emit(dt, &mut rng);
        for (pos, vel) in spawns {
            let alive = self.particles.iter().filter(|p| p.alive).count();
            if alive >= self.config.max_particles {
                break;
            }
            let lt = rng.range_f64(self.config.lifetime_range[0], self.config.lifetime_range[1]);
            let p = WasmParticle::new(pos, vel, lt);
            let slot = self.particles.iter().position(|p| !p.alive);
            if let Some(i) = slot {
                self.particles[i] = p;
            } else {
                self.particles.push(p);
            }
        }
    }

    /// Add a force.
    pub fn add_force(&mut self, force: WasmParticleForce) {
        self.forces.push(force);
    }

    /// Add a collider.
    pub fn add_collider(&mut self, collider: WasmParticleCollider) {
        self.colliders.push(collider);
    }

    /// Step the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.time += dt;
        let damping = 1.0 - self.config.damping * dt;
        for p in &mut self.particles {
            if !p.alive {
                continue;
            }
            p.age += dt;
            if p.is_dead() {
                p.alive = false;
                continue;
            }
            // Accumulate forces
            let mut acc = [0.0f64; 3];
            for force in &self.forces {
                let fa = force.acceleration(p.position, self.time);
                acc = add3(acc, fa);
            }
            // Integrate velocity and position (semi-implicit Euler)
            p.velocity = add3(p.velocity, scale3(acc, dt));
            p.velocity = scale3(p.velocity, damping.max(0.0));
            p.position = add3(p.position, scale3(p.velocity, dt));
            // Resolve collisions
            for collider in &self.colliders {
                let (np, nv) = collider.resolve(p.position, p.velocity);
                p.position = np;
                p.velocity = nv;
            }
        }
    }

    /// Return flat position array `[x0,y0,z0, x1,y1,z1, ...]` for alive particles.
    pub fn get_positions(&self) -> Vec<f64> {
        let mut out = Vec::new();
        for p in &self.particles {
            if p.alive {
                out.push(p.position[0]);
                out.push(p.position[1]);
                out.push(p.position[2]);
            }
        }
        out
    }

    /// Return flat color array `[r0,g0,b0,a0, r1,g1,b1,a1, ...]` for alive particles.
    pub fn get_colors(&self) -> Vec<f32> {
        let mut out = Vec::new();
        for p in &self.particles {
            if p.alive {
                out.extend_from_slice(&p.color);
            }
        }
        out
    }

    /// Return the number of alive particles.
    pub fn alive_count(&self) -> usize {
        self.particles.iter().filter(|p| p.alive).count()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Config tests ---

    #[test]
    fn test_default_config() {
        let cfg = WasmParticleConfig::default();
        assert_eq!(cfg.max_particles, 1000);
        assert!((cfg.gravity[1] - (-9.81)).abs() < 1e-10);
    }

    // --- Particle tests ---

    #[test]
    fn test_particle_new() {
        let p = WasmParticle::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], 5.0);
        assert!(p.alive);
        assert!((p.age).abs() < 1e-10);
        assert!((p.lifetime - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_normalized_age() {
        let mut p = WasmParticle::new([0.0; 3], [0.0; 3], 4.0);
        p.age = 2.0;
        assert!((p.normalized_age() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_particle_is_dead() {
        let mut p = WasmParticle::new([0.0; 3], [0.0; 3], 1.0);
        p.age = 1.5;
        assert!(p.is_dead());
    }

    // --- Emitter shape tests ---

    #[test]
    fn test_emitter_point_shape() {
        let mut rng = SimpleRng::new(1);
        let pos = WasmEmitterShape::Point.sample_position(&mut rng, [1.0, 2.0, 3.0]);
        assert_eq!(pos, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_emitter_sphere_shape_in_bounds() {
        let mut rng = SimpleRng::new(2);
        let shape = WasmEmitterShape::Sphere(5.0);
        for _ in 0..20 {
            let p = shape.sample_position(&mut rng, [0.0; 3]);
            let dist = len3(p);
            assert!(dist <= 5.0 + 1e-9, "dist={dist}");
        }
    }

    #[test]
    fn test_emitter_disc_shape_y_zero() {
        let mut rng = SimpleRng::new(3);
        let shape = WasmEmitterShape::Disc(2.0);
        for _ in 0..10 {
            let p = shape.sample_position(&mut rng, [0.0; 3]);
            assert!((p[1]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_emitter_emit_count() {
        let mut emitter = WasmEmitter::new(WasmEmitterShape::Point, 10.0, [0.0; 3]);
        let mut rng = SimpleRng::new(42);
        let spawns = emitter.emit(1.0, &mut rng);
        assert_eq!(spawns.len(), 10);
    }

    #[test]
    fn test_emitter_burst() {
        let mut emitter = WasmEmitter::new(WasmEmitterShape::Point, 0.0, [0.0; 3]);
        emitter.burst_count = 5;
        let mut rng = SimpleRng::new(1);
        let spawns = emitter.emit(0.016, &mut rng);
        assert_eq!(spawns.len(), 5);
    }

    // --- Force tests ---

    #[test]
    fn test_gravity_force() {
        let f = WasmParticleForce::Gravity([0.0, -9.81, 0.0]);
        let acc = f.acceleration([0.0; 3], 0.0);
        assert!((acc[1] - (-9.81)).abs() < 1e-10);
    }

    #[test]
    fn test_wind_force() {
        let f = WasmParticleForce::Wind([3.0, 0.0, 0.0]);
        let acc = f.acceleration([0.0; 3], 0.0);
        assert!((acc[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_attractor_force_direction() {
        let f = WasmParticleForce::Attractor {
            pos: [1.0, 0.0, 0.0],
            strength: 1.0,
        };
        let acc = f.acceleration([0.0; 3], 0.0);
        // Should pull toward +x
        assert!(acc[0] > 0.0);
    }

    #[test]
    fn test_turbulence_force_nonzero() {
        let f = WasmParticleForce::Turbulence {
            freq: 1.0,
            amplitude: 0.5,
        };
        let acc = f.acceleration([0.1, 0.2, 0.3], 1.0);
        // At least one component should be nonzero
        let mag: f64 = acc.iter().map(|&x| x * x).sum::<f64>().sqrt();
        assert!(mag > 0.0);
    }

    // --- Collider tests ---

    #[test]
    fn test_plane_collider_bounce() {
        let col = WasmParticleCollider::plane([0.0, 1.0, 0.0], 0.0, 0.8);
        // Particle below plane moving downward
        let (new_pos, new_vel) = col.resolve([0.0, -0.1, 0.0], [0.0, -1.0, 0.0]);
        assert!(new_pos[1] >= 0.0);
        assert!(new_vel[1] > 0.0); // bounced upward
    }

    #[test]
    fn test_sphere_collider_push_out() {
        let col = WasmParticleCollider::sphere([0.0; 3], 1.0, 0.5);
        let (new_pos, _) = col.resolve([0.3, 0.0, 0.0], [0.0; 3]);
        assert!(len3(new_pos) >= 1.0 - 1e-5);
    }

    // --- Trail tests ---

    #[test]
    fn test_trail_push_capped() {
        let mut trail = WasmTrail::new(3, false);
        for i in 0..10 {
            trail.push([i as f64, 0.0, 0.0], i as f64);
        }
        assert!(trail.points.len() <= 3);
    }

    #[test]
    fn test_trail_fading_width() {
        let trail = WasmTrail::new(5, true);
        let w0 = trail.width_at(0.0, 1.0);
        let w1 = trail.width_at(1.0, 1.0);
        assert!((w0).abs() < 1e-10);
        assert!((w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_trail_system_update() {
        let mut ts = WasmTrailSystem::new(5, 10, false);
        ts.update(0, [1.0, 0.0, 0.0], 0.5);
        assert_eq!(ts.trails[0].points.len(), 1);
    }

    #[test]
    fn test_trail_system_clear() {
        let mut ts = WasmTrailSystem::new(3, 5, false);
        ts.update(1, [0.0; 3], 0.0);
        ts.clear(1);
        assert!(ts.trails[1].points.is_empty());
    }

    // --- Renderer tests ---

    #[test]
    fn test_renderer_size_at_edges() {
        let r = WasmParticleRenderer::default_white();
        let s0 = r.size_at(0.0);
        assert!((s0 - 0.1).abs() < 1e-6);
        let s1 = r.size_at(1.0);
        assert!((s1).abs() < 1e-6);
    }

    #[test]
    fn test_renderer_alpha_at() {
        let r = WasmParticleRenderer::default_white();
        assert!((r.alpha_at(0.0) - 1.0).abs() < 1e-6);
        assert!((r.alpha_at(1.0)).abs() < 1e-6);
    }

    #[test]
    fn test_renderer_color_at_midpoint() {
        let r = WasmParticleRenderer::default_white();
        let c = r.color_at(0.5);
        assert!(c[3] >= 0.0 && c[3] <= 1.0);
    }

    // --- Fire particle tests ---

    #[test]
    fn test_fire_particle_step() {
        let mut fp = WasmFireParticle::new([0.0; 3], 1000.0, 5.0);
        let start_y = fp.particle.position[1];
        fp.step(0.1);
        assert!(fp.particle.position[1] > start_y); // rose due to heat
        assert!(fp.particle.age > 0.0);
    }

    #[test]
    fn test_fire_particle_cools() {
        let mut fp = WasmFireParticle::new([0.0; 3], 1000.0, 2.0);
        fp.particle.age = 1.9; // near end of life
        fp.step(0.05);
        assert!(fp.temperature < 1000.0);
    }

    // --- WasmParticleSystem tests ---

    #[test]
    fn test_system_spawn_and_count() {
        let mut sys = WasmParticleSystem::earth();
        sys.spawn(10);
        assert_eq!(sys.alive_count(), 10);
    }

    #[test]
    fn test_system_step_moves_particles() {
        let mut sys = WasmParticleSystem::earth();
        sys.spawn(5);
        let before: Vec<[f64; 3]> = sys
            .particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| p.position)
            .collect();
        sys.step(0.1);
        let after: Vec<[f64; 3]> = sys
            .particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| p.position)
            .collect();
        // At least some should have moved
        let moved = before
            .iter()
            .zip(after.iter())
            .any(|(b, a)| len3(sub3(*a, *b)) > 1e-10);
        assert!(moved);
    }

    #[test]
    fn test_system_get_positions_length() {
        let mut sys = WasmParticleSystem::earth();
        sys.spawn(7);
        let pos = sys.get_positions();
        assert_eq!(pos.len(), 7 * 3);
    }

    #[test]
    fn test_system_get_colors_length() {
        let mut sys = WasmParticleSystem::earth();
        sys.spawn(4);
        let colors = sys.get_colors();
        assert_eq!(colors.len(), 4 * 4);
    }

    #[test]
    fn test_system_max_particles_respected() {
        let cfg = WasmParticleConfig {
            max_particles: 5,
            ..Default::default()
        };
        let mut sys = WasmParticleSystem::new(cfg);
        sys.spawn(100);
        assert!(sys.alive_count() <= 5);
    }

    #[test]
    fn test_system_particles_die_over_time() {
        let cfg = WasmParticleConfig {
            lifetime_range: [0.1, 0.1],
            gravity: [0.0; 3],
            ..Default::default()
        };
        let mut sys = WasmParticleSystem::new(cfg);
        sys.spawn(5);
        // Step past lifetime
        for _ in 0..20 {
            sys.step(0.02);
        }
        assert_eq!(sys.alive_count(), 0);
    }

    #[test]
    fn test_system_with_plane_collider() {
        let mut sys = WasmParticleSystem::earth();
        sys.add_collider(WasmParticleCollider::plane([0.0, 1.0, 0.0], 0.0, 0.5));
        sys.spawn(5);
        for _ in 0..100 {
            sys.step(0.016);
        }
        for p in sys.particles.iter().filter(|p| p.alive) {
            assert!(p.position[1] >= -0.01, "y={}", p.position[1]);
        }
    }
}
