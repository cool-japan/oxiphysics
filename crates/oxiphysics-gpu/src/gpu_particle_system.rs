// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated particle system (CPU mock implementation).
//!
//! Provides a particle emitter, Euler integrator, and depth-sort for
//! alpha-blended rendering.  All operations are implemented on the CPU as a
//! reference / fallback backend.

/// A single particle in the GPU particle system.
#[derive(Debug, Clone, PartialEq)]
pub struct GpuParticle {
    /// World-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]` in world units per second.
    pub velocity: [f64; 3],
    /// Remaining lifetime in seconds (`<= 0` means dead).
    pub lifetime: f64,
    /// RGBA colour, each component in `[0, 1]`.
    pub color: [f32; 4],
}

impl GpuParticle {
    /// Create a new particle at the given position.
    pub fn new(position: [f64; 3], velocity: [f64; 3], lifetime: f64, color: [f32; 4]) -> Self {
        Self {
            position,
            velocity,
            lifetime,
            color,
        }
    }

    /// Returns `true` if the particle is still alive.
    pub fn is_alive(&self) -> bool {
        self.lifetime > 0.0
    }
}

/// Configuration for a point-emitter.
#[derive(Debug, Clone)]
pub struct EmitterConfig {
    /// World-space origin of the emitter.
    pub origin: [f64; 3],
    /// Initial speed applied along the emission direction.
    pub initial_speed: f64,
    /// Spread half-angle in radians (0 = directional).
    pub spread_radians: f64,
    /// Lifetime in seconds for freshly spawned particles.
    pub particle_lifetime: f64,
    /// RGBA colour assigned to every new particle.
    pub color: [f32; 4],
}

impl Default for EmitterConfig {
    fn default() -> Self {
        Self {
            origin: [0.0; 3],
            initial_speed: 1.0,
            spread_radians: 0.3,
            particle_lifetime: 2.0,
            color: [1.0, 1.0, 1.0, 1.0],
        }
    }
}

/// GPU particle system managing emission and integration of many particles.
#[derive(Debug, Clone)]
pub struct GpuParticleSystem {
    /// Emitter configuration.
    pub config: EmitterConfig,
    /// Currently active particles.
    pub particles: Vec<GpuParticle>,
    /// Maximum number of particles allowed simultaneously.
    pub max_particles: usize,
}

impl GpuParticleSystem {
    /// Create a new particle system with the given emitter config and capacity.
    pub fn new(config: EmitterConfig, max_particles: usize) -> Self {
        Self {
            config,
            particles: Vec::with_capacity(max_particles),
            max_particles,
        }
    }

    /// Return the number of currently alive particles.
    pub fn active_count(&self) -> usize {
        self.particles.len()
    }
}

// ── Core GPU-mock operations ──────────────────────────────────────────────────

/// Spawn `n` particles from the emitter into the system.
///
/// Uses a deterministic direction fan spread around `+Z` so results are
/// reproducible without a random number generator.  Respects `max_particles`.
pub fn gpu_emit_particles(system: &mut GpuParticleSystem, n: usize) {
    let cfg = &system.config;
    let slots = system.max_particles.saturating_sub(system.particles.len());
    let to_spawn = n.min(slots);

    for i in 0..to_spawn {
        // Deterministic spread: fan angle around +Z axis.
        let angle = if to_spawn > 1 {
            let t = i as f64 / (to_spawn - 1) as f64;
            (t - 0.5) * 2.0 * cfg.spread_radians
        } else {
            0.0
        };
        let vx = angle.sin() * cfg.initial_speed;
        let vz = angle.cos() * cfg.initial_speed;
        let velocity = [vx, 0.0, vz];
        system.particles.push(GpuParticle::new(
            cfg.origin,
            velocity,
            cfg.particle_lifetime,
            cfg.color,
        ));
    }
}

/// Advance all particles by `dt` seconds using symplectic Euler integration.
///
/// Also decrements each particle's `lifetime` by `dt`.
pub fn gpu_integrate_particles(system: &mut GpuParticleSystem, dt: f64) {
    for p in &mut system.particles {
        p.position[0] += p.velocity[0] * dt;
        p.position[1] += p.velocity[1] * dt;
        p.position[2] += p.velocity[2] * dt;
        p.lifetime -= dt;
    }
}

/// Remove particles whose `lifetime <= 0`, compacting the particle list.
///
/// This mirrors the GPU stream-compaction pattern.
pub fn gpu_kill_dead_particles(system: &mut GpuParticleSystem) {
    system.particles.retain(|p| p.is_alive());
}

/// Sort active particles by depth along `camera_dir` (back-to-front) for
/// correct alpha blending.
///
/// `camera_dir` should be the normalised view direction (world space).  The
/// sort is stable so that ties preserve the original emission order.
pub fn gpu_sort_by_depth(system: &mut GpuParticleSystem, camera_dir: [f64; 3]) {
    system.particles.sort_by(|a, b| {
        let da = dot3(a.position, camera_dir);
        let db = dot3(b.position, camera_dir);
        // Back-to-front: largest depth first.
        db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
    });
}

/// Emit all `n` particles at once (burst emission).
///
/// Equivalent to calling `gpu_emit_particles` once with `n`.
pub fn spawn_burst(system: &mut GpuParticleSystem, n: usize) {
    gpu_emit_particles(system, n);
}

// ── Internal helpers ──────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_system(max: usize) -> GpuParticleSystem {
        GpuParticleSystem::new(EmitterConfig::default(), max)
    }

    // --- GpuParticle ---

    #[test]
    fn test_particle_is_alive_positive_lifetime() {
        let p = GpuParticle::new([0.0; 3], [0.0; 3], 1.0, [1.0; 4]);
        assert!(p.is_alive());
    }

    #[test]
    fn test_particle_is_dead_zero_lifetime() {
        let p = GpuParticle::new([0.0; 3], [0.0; 3], 0.0, [1.0; 4]);
        assert!(!p.is_alive());
    }

    #[test]
    fn test_particle_is_dead_negative_lifetime() {
        let p = GpuParticle::new([0.0; 3], [0.0; 3], -1.0, [1.0; 4]);
        assert!(!p.is_alive());
    }

    #[test]
    fn test_particle_new_fields() {
        let pos = [1.0, 2.0, 3.0];
        let vel = [0.1, 0.2, 0.3];
        let lt = 5.0;
        let col = [0.5, 0.5, 0.5, 1.0];
        let p = GpuParticle::new(pos, vel, lt, col);
        assert_eq!(p.position, pos);
        assert_eq!(p.velocity, vel);
        assert!((p.lifetime - lt).abs() < 1e-12);
        assert_eq!(p.color, col);
    }

    // --- EmitterConfig ---

    #[test]
    fn test_emitter_default_speed_positive() {
        let cfg = EmitterConfig::default();
        assert!(cfg.initial_speed > 0.0);
    }

    #[test]
    fn test_emitter_default_lifetime_positive() {
        let cfg = EmitterConfig::default();
        assert!(cfg.particle_lifetime > 0.0);
    }

    // --- GpuParticleSystem ---

    #[test]
    fn test_system_starts_empty() {
        let sys = default_system(100);
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn test_system_max_particles_stored() {
        let sys = default_system(42);
        assert_eq!(sys.max_particles, 42);
    }

    // --- gpu_emit_particles ---

    #[test]
    fn test_emit_spawns_n_particles() {
        let mut sys = default_system(100);
        gpu_emit_particles(&mut sys, 10);
        assert_eq!(sys.active_count(), 10);
    }

    #[test]
    fn test_emit_respects_max_particles() {
        let mut sys = default_system(5);
        gpu_emit_particles(&mut sys, 100);
        assert_eq!(sys.active_count(), 5);
    }

    #[test]
    fn test_emit_zero_particles() {
        let mut sys = default_system(100);
        gpu_emit_particles(&mut sys, 0);
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn test_emit_single_particle_at_origin() {
        let mut sys = default_system(10);
        gpu_emit_particles(&mut sys, 1);
        assert_eq!(sys.particles[0].position, [0.0; 3]);
    }

    #[test]
    fn test_emit_particles_have_positive_lifetime() {
        let mut sys = default_system(10);
        gpu_emit_particles(&mut sys, 5);
        for p in &sys.particles {
            assert!(p.lifetime > 0.0);
        }
    }

    // --- gpu_integrate_particles ---

    #[test]
    fn test_integrate_moves_position() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([0.0; 3], [1.0, 0.0, 0.0], 5.0, [1.0; 4]));
        gpu_integrate_particles(&mut sys, 1.0);
        assert!((sys.particles[0].position[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_integrate_decrements_lifetime() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([0.0; 3], [0.0; 3], 3.0, [1.0; 4]));
        gpu_integrate_particles(&mut sys, 1.0);
        assert!((sys.particles[0].lifetime - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_integrate_zero_dt_no_movement() {
        let mut sys = default_system(10);
        sys.particles.push(GpuParticle::new(
            [1.0, 2.0, 3.0],
            [5.0, 5.0, 5.0],
            1.0,
            [1.0; 4],
        ));
        gpu_integrate_particles(&mut sys, 0.0);
        assert_eq!(sys.particles[0].position, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_integrate_multiple_steps() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([0.0; 3], [2.0, 0.0, 0.0], 10.0, [1.0; 4]));
        gpu_integrate_particles(&mut sys, 0.5);
        gpu_integrate_particles(&mut sys, 0.5);
        assert!((sys.particles[0].position[0] - 2.0).abs() < 1e-10);
    }

    // --- gpu_kill_dead_particles ---

    #[test]
    fn test_kill_removes_dead_particles() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([0.0; 3], [0.0; 3], -1.0, [1.0; 4]));
        sys.particles
            .push(GpuParticle::new([0.0; 3], [0.0; 3], 1.0, [1.0; 4]));
        gpu_kill_dead_particles(&mut sys);
        assert_eq!(sys.active_count(), 1);
        assert!(sys.particles[0].is_alive());
    }

    #[test]
    fn test_kill_all_dead() {
        let mut sys = default_system(10);
        for _ in 0..5 {
            sys.particles
                .push(GpuParticle::new([0.0; 3], [0.0; 3], -1.0, [1.0; 4]));
        }
        gpu_kill_dead_particles(&mut sys);
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn test_kill_none_dead() {
        let mut sys = default_system(10);
        gpu_emit_particles(&mut sys, 5);
        gpu_kill_dead_particles(&mut sys);
        assert_eq!(sys.active_count(), 5);
    }

    #[test]
    fn test_integrate_then_kill() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, [1.0; 4]));
        gpu_integrate_particles(&mut sys, 1.0); // lifetime becomes -0.5
        gpu_kill_dead_particles(&mut sys);
        assert_eq!(sys.active_count(), 0);
    }

    // --- gpu_sort_by_depth ---

    #[test]
    fn test_sort_by_depth_back_to_front() {
        let mut sys = default_system(10);
        // Two particles along Z axis
        sys.particles
            .push(GpuParticle::new([0.0, 0.0, 1.0], [0.0; 3], 1.0, [1.0; 4]));
        sys.particles
            .push(GpuParticle::new([0.0, 0.0, 5.0], [0.0; 3], 1.0, [1.0; 4]));
        let cam = [0.0, 0.0, 1.0]; // camera looks along +Z
        gpu_sort_by_depth(&mut sys, cam);
        // Particle at z=5 should come first (further from camera = drawn first)
        assert!((sys.particles[0].position[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_sort_by_depth_single_particle() {
        let mut sys = default_system(10);
        sys.particles
            .push(GpuParticle::new([1.0, 2.0, 3.0], [0.0; 3], 1.0, [1.0; 4]));
        gpu_sort_by_depth(&mut sys, [0.0, 0.0, 1.0]);
        assert_eq!(sys.active_count(), 1);
    }

    #[test]
    fn test_sort_by_depth_preserves_count() {
        let mut sys = default_system(20);
        gpu_emit_particles(&mut sys, 10);
        gpu_sort_by_depth(&mut sys, [1.0, 0.0, 0.0]);
        assert_eq!(sys.active_count(), 10);
    }

    // --- spawn_burst ---

    #[test]
    fn test_spawn_burst_emits_all_at_once() {
        let mut sys = default_system(50);
        spawn_burst(&mut sys, 20);
        assert_eq!(sys.active_count(), 20);
    }

    #[test]
    fn test_spawn_burst_respects_max() {
        let mut sys = default_system(5);
        spawn_burst(&mut sys, 100);
        assert_eq!(sys.active_count(), 5);
    }

    // --- Integration scenario ---

    #[test]
    fn test_full_lifecycle() {
        let mut sys = default_system(100);
        spawn_burst(&mut sys, 30);
        assert_eq!(sys.active_count(), 30);
        // Integrate past the lifetime
        let dt = EmitterConfig::default().particle_lifetime + 0.1;
        gpu_integrate_particles(&mut sys, dt);
        gpu_kill_dead_particles(&mut sys);
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn test_emission_incremental() {
        let mut sys = default_system(100);
        gpu_emit_particles(&mut sys, 10);
        gpu_emit_particles(&mut sys, 10);
        assert_eq!(sys.active_count(), 20);
    }

    #[test]
    fn test_particle_color_propagated() {
        let cfg = EmitterConfig {
            color: [1.0, 0.0, 0.0, 1.0],
            ..Default::default()
        };
        let mut sys = GpuParticleSystem::new(cfg, 10);
        gpu_emit_particles(&mut sys, 1);
        assert_eq!(sys.particles[0].color, [1.0, 0.0, 0.0, 1.0]);
    }
}
