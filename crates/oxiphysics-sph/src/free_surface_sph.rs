// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Free-surface SPH.
//!
//! Implements free-surface detection, normal estimation, ghost-particle
//! reflection, wave-height profiling, and volume-conservation monitoring.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// § 1  Particle and parameter structs
// ─────────────────────────────────────────────────────────────────────────────

/// A single particle for free-surface SPH simulations.
#[derive(Debug, Clone)]
pub struct FreeSurfaceParticle {
    /// Position \[x, y, z\] \[m\].
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] \[m/s\].
    pub velocity: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// SPH density estimate \[kg/m³\].
    pub density: f64,
    /// Pressure \[Pa\].
    pub pressure: f64,
    /// Whether this particle is on the free surface.
    pub is_free_surface: bool,
    /// Outward unit normal at the free surface.
    pub normal: [f64; 3],
}

impl FreeSurfaceParticle {
    /// Create a new interior particle.
    pub fn new(position: [f64; 3], mass: f64, density: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density,
            pressure: 0.0,
            is_free_surface: false,
            normal: [0.0; 3],
        }
    }

    /// Squared distance to another particle.
    pub fn dist2(&self, other: &FreeSurfaceParticle) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        dx * dx + dy * dy + dz * dz
    }
}

/// Parameters controlling free-surface detection and wave analysis.
#[derive(Debug, Clone)]
pub struct FreeSurfaceParams {
    /// Density threshold factor β: particle is free surface if ρ < β ρ₀.
    pub threshold_factor: f64,
    /// SPH kernel support radius h \[m\].
    pub kernel_radius: f64,
    /// Gravitational acceleration vector \[m/s²\].
    pub gravity: [f64; 3],
}

impl FreeSurfaceParams {
    /// Create a new `FreeSurfaceParams` with standard downward gravity.
    pub fn new(threshold_factor: f64, kernel_radius: f64) -> Self {
        Self {
            threshold_factor,
            kernel_radius,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Create with a custom gravity vector.
    pub fn with_gravity(threshold_factor: f64, kernel_radius: f64, gravity: [f64; 3]) -> Self {
        Self {
            threshold_factor,
            kernel_radius,
            gravity,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Kernel helper (cubic spline, 3-D)
// ─────────────────────────────────────────────────────────────────────────────

fn cubic_kernel_grad_3d(rij: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if r < 1e-30 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dq = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        alpha * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    };
    let dw_dr = dw_dq / h;
    [dw_dr * rij[0] / r, dw_dr * rij[1] / r, dw_dr * rij[2] / r]
}

// ─────────────────────────────────────────────────────────────────────────────
// § 2  Free-surface detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect free-surface particles using the density deficiency criterion.
///
/// A particle is marked as free-surface if its density is less than
/// `params.threshold_factor * ρ_max`, where ρ_max is the maximum density
/// in the particle set.
///
/// # Arguments
/// * `particles` – mutable particle slice
/// * `params`    – detection parameters
pub fn detect_free_surface(particles: &mut [FreeSurfaceParticle], params: &FreeSurfaceParams) {
    if particles.is_empty() {
        return;
    }
    let rho_max = particles
        .iter()
        .map(|p| p.density)
        .fold(f64::NEG_INFINITY, f64::max);
    let threshold = params.threshold_factor * rho_max;
    for p in particles.iter_mut() {
        p.is_free_surface = p.density < threshold;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 3  Free-surface normal estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the outward normal at a free-surface particle via SPH colour-function gradient.
///
/// Uses a smoothed colour function (1 = fluid, 0 = air) and computes ∇C.
///
/// # Arguments
/// * `particles` – all particles
/// * `i`         – query particle index
/// * `h`         – smoothing length \[m\]
pub fn free_surface_normal_estimation(
    particles: &[FreeSurfaceParticle],
    i: usize,
    h: f64,
) -> [f64; 3] {
    let xi = particles[i].position;
    let mut grad = [0.0_f64; 3];
    for p in particles {
        let rij = [
            xi[0] - p.position[0],
            xi[1] - p.position[1],
            xi[2] - p.position[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        let dw = cubic_kernel_grad_3d(rij, r, h);
        let vol_j = p.mass / p.density.max(1e-30);
        // Colour function: 1 for fluid particles, 0 for free-surface
        let c_j = if p.is_free_surface { 0.0 } else { 1.0 };
        grad[0] += vol_j * c_j * dw[0];
        grad[1] += vol_j * c_j * dw[1];
        grad[2] += vol_j * c_j * dw[2];
    }
    grad
}

// ─────────────────────────────────────────────────────────────────────────────
// § 4  Zero pressure at free surface
// ─────────────────────────────────────────────────────────────────────────────

/// Enforce p = 0 boundary condition at free-surface particles.
///
/// Sets the pressure of every particle flagged as `is_free_surface` to zero.
pub fn pressure_zero_free_surface(particles: &mut [FreeSurfaceParticle]) {
    for p in particles.iter_mut() {
        if p.is_free_surface {
            p.pressure = 0.0;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 5  Ghost particle reflection
// ─────────────────────────────────────────────────────────────────────────────

/// Create a ghost particle by mirroring `particle` across the surface defined
/// by its outward `normal`.
///
/// The ghost is placed at `x_ghost = x + 2 h_mirror normal` with mirrored
/// velocity and the same mass/density.  `h_mirror = kernel_radius` is used as
/// the mirror distance.
///
/// # Arguments
/// * `particle`      – original free-surface particle
/// * `normal`        – outward surface normal (need not be normalised)
/// * `kernel_radius` – mirror distance \[m\]
pub fn ghost_particle_reflection(
    particle: &FreeSurfaceParticle,
    normal: [f64; 3],
    kernel_radius: f64,
) -> FreeSurfaceParticle {
    let n_mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let n_hat = if n_mag > 1e-30 {
        [normal[0] / n_mag, normal[1] / n_mag, normal[2] / n_mag]
    } else {
        [0.0, 1.0, 0.0] // default upward
    };
    let d = 2.0 * kernel_radius;
    let ghost_pos = [
        particle.position[0] + d * n_hat[0],
        particle.position[1] + d * n_hat[1],
        particle.position[2] + d * n_hat[2],
    ];
    // Mirror the normal component of velocity
    let v = particle.velocity;
    let v_dot_n = v[0] * n_hat[0] + v[1] * n_hat[1] + v[2] * n_hat[2];
    let ghost_vel = [
        v[0] - 2.0 * v_dot_n * n_hat[0],
        v[1] - 2.0 * v_dot_n * n_hat[1],
        v[2] - 2.0 * v_dot_n * n_hat[2],
    ];
    FreeSurfaceParticle {
        position: ghost_pos,
        velocity: ghost_vel,
        mass: particle.mass,
        density: particle.density,
        pressure: 0.0,
        is_free_surface: false,
        normal: n_hat,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 6  Wave height profile
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the wave height profile along the x-axis.
///
/// Divides \[x_min, x_max\] into `n_bins` columns; for each bin, finds the
/// maximum y-coordinate among free-surface particles.  Returns 0.0 for empty
/// bins.
///
/// # Arguments
/// * `particles` – all particles
/// * `x_range`   – (x_min, x_max) extent of the domain
/// * `n_bins`    – number of horizontal bins
pub fn wave_height_profile(
    particles: &[FreeSurfaceParticle],
    x_range: (f64, f64),
    n_bins: usize,
) -> Vec<f64> {
    if n_bins == 0 {
        return vec![];
    }
    let (x_min, x_max) = x_range;
    let dx = (x_max - x_min) / n_bins as f64;
    let mut heights = vec![f64::NEG_INFINITY; n_bins];
    for p in particles {
        if !p.is_free_surface {
            continue;
        }
        let x = p.position[0];
        if x < x_min || x >= x_max {
            continue;
        }
        let bin = ((x - x_min) / dx) as usize;
        let bin = bin.min(n_bins - 1);
        if p.position[1] > heights[bin] {
            heights[bin] = p.position[1];
        }
    }
    heights
        .iter()
        .map(|&h| if h == f64::NEG_INFINITY { 0.0 } else { h })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 7  Wave speed estimate
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the shallow-water wave speed from a height profile.
///
/// `c = sqrt(g * h_mean)` where h_mean is the mean height in the profile.
///
/// # Arguments
/// * `height_profile` – wave height per bin \[m\]
/// * `dx`             – bin width \[m\] (informational; not used in formula)
/// * `g`              – gravitational acceleration magnitude \[m/s²\]
pub fn wave_speed_estimate(height_profile: &[f64], _dx: f64, g: f64) -> f64 {
    if height_profile.is_empty() {
        return 0.0;
    }
    let h_mean = height_profile.iter().sum::<f64>() / height_profile.len() as f64;
    if h_mean <= 0.0 {
        return 0.0;
    }
    (g * h_mean).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 8  Splash detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect airborne (splash) particles: free-surface particles above `min_height`.
///
/// Returns the indices of particles satisfying both conditions.
///
/// # Arguments
/// * `particles`  – all particles
/// * `min_height` – minimum y-coordinate \[m\] above which a particle is "airborne"
pub fn splash_detection(particles: &[FreeSurfaceParticle], min_height: f64) -> Vec<usize> {
    particles
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_free_surface && p.position[1] > min_height)
        .map(|(i, _)| i)
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 9  Volume conservation error
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the relative volume conservation error.
///
/// Volume is approximated as `V = Σ m_i / ρ_i`.
///
/// Returns `(V_current - V_initial) / V_initial`.
///
/// # Arguments
/// * `particles`      – all particles
/// * `initial_volume` – initial fluid volume V₀ \[m³\]
pub fn volume_conservation_error(particles: &[FreeSurfaceParticle], initial_volume: f64) -> f64 {
    let v_current: f64 = particles
        .iter()
        .map(|p| p.mass / p.density.max(1e-30))
        .sum();
    if initial_volume.abs() < 1e-30 {
        return 0.0;
    }
    (v_current - initial_volume) / initial_volume
}

// ─────────────────────────────────────────────────────────────────────────────
// § 10  Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FreeSurfaceParticle ──────────────────────────────────────────────────

    #[test]
    fn test_fsp_new() {
        let p = FreeSurfaceParticle::new([1.0, 2.0, 3.0], 0.5, 1000.0);
        assert!((p.position[0] - 1.0).abs() < 1e-12);
        assert!((p.mass - 0.5).abs() < 1e-12);
        assert!(!p.is_free_surface);
        assert_eq!(p.pressure, 0.0);
    }

    #[test]
    fn test_fsp_dist2() {
        let p1 = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        let p2 = FreeSurfaceParticle::new([3.0, 4.0, 0.0], 1.0, 1000.0);
        assert!((p1.dist2(&p2) - 25.0).abs() < 1e-12);
    }

    // ── FreeSurfaceParams ────────────────────────────────────────────────────

    #[test]
    fn test_params_new() {
        let params = FreeSurfaceParams::new(0.85, 0.1);
        assert!((params.threshold_factor - 0.85).abs() < 1e-12);
        assert!((params.kernel_radius - 0.1).abs() < 1e-12);
        assert!((params.gravity[1] + 9.81).abs() < 1e-12);
    }

    #[test]
    fn test_params_with_gravity() {
        let params = FreeSurfaceParams::with_gravity(0.8, 0.05, [0.0, -3.7, 0.0]);
        assert!((params.gravity[1] + 3.7).abs() < 1e-12);
    }

    // ── detect_free_surface ──────────────────────────────────────────────────

    fn make_particles() -> Vec<FreeSurfaceParticle> {
        let mut v: Vec<FreeSurfaceParticle> = (0..8)
            .map(|i| FreeSurfaceParticle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 1000.0))
            .collect();
        // Make top 2 particles have lower density (free surface)
        v[6].density = 600.0;
        v[7].density = 500.0;
        v
    }

    #[test]
    fn test_detect_free_surface_marks_low_density() {
        let mut particles = make_particles();
        let params = FreeSurfaceParams::new(0.75, 0.1);
        detect_free_surface(&mut particles, &params);
        assert!(particles[6].is_free_surface);
        assert!(particles[7].is_free_surface);
        assert!(!particles[0].is_free_surface);
    }

    #[test]
    fn test_detect_free_surface_empty() {
        let mut particles: Vec<FreeSurfaceParticle> = vec![];
        let params = FreeSurfaceParams::new(0.85, 0.1);
        detect_free_surface(&mut particles, &params); // should not panic
    }

    #[test]
    fn test_detect_free_surface_all_same_density() {
        let mut particles: Vec<FreeSurfaceParticle> = (0..5)
            .map(|_| FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0))
            .collect();
        let params = FreeSurfaceParams::new(0.85, 0.1);
        detect_free_surface(&mut particles, &params);
        // ρ = 0.85 * ρ_max = 850 < 1000 → none are free surface
        assert!(particles.iter().all(|p| !p.is_free_surface));
    }

    // ── pressure_zero_free_surface ───────────────────────────────────────────

    #[test]
    fn test_pressure_zero_free_surface() {
        let mut particles = make_particles();
        particles[0].is_free_surface = true;
        particles[0].pressure = 500.0;
        particles[1].pressure = 200.0;
        pressure_zero_free_surface(&mut particles);
        assert_eq!(particles[0].pressure, 0.0);
        assert!((particles[1].pressure - 200.0).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_zero_not_free_surface_unchanged() {
        let mut particles: Vec<FreeSurfaceParticle> = (0..3)
            .map(|_| FreeSurfaceParticle::new([0.0; 3], 1.0, 1000.0))
            .collect();
        for p in &mut particles {
            p.pressure = 100.0;
        }
        pressure_zero_free_surface(&mut particles);
        // None are free surface, so pressures are unchanged
        for p in &particles {
            assert!((p.pressure - 100.0).abs() < 1e-12);
        }
    }

    // ── ghost_particle_reflection ────────────────────────────────────────────

    #[test]
    fn test_ghost_particle_position() {
        let mut p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        p.velocity = [1.0, 2.0, 0.0];
        let normal = [0.0, 1.0, 0.0]; // upward normal
        let ghost = ghost_particle_reflection(&p, normal, 0.1);
        // Ghost should be 2*h=0.2 above particle
        assert!((ghost.position[1] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_ghost_velocity_mirror() {
        let mut p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        p.velocity = [0.0, -1.0, 0.0]; // downward velocity
        let normal = [0.0, 1.0, 0.0]; // upward normal
        let ghost = ghost_particle_reflection(&p, normal, 0.1);
        // Normal component of velocity reflected: -1 → +1
        assert!((ghost.velocity[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ghost_zero_normal_default() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        // Zero normal → default upward ghost
        let ghost = ghost_particle_reflection(&p, [0.0, 0.0, 0.0], 0.1);
        assert!(ghost.position[1] > 0.0);
    }

    #[test]
    fn test_ghost_mass_density_preserved() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 2.0, 500.0);
        let ghost = ghost_particle_reflection(&p, [1.0, 0.0, 0.0], 0.1);
        assert!((ghost.mass - 2.0).abs() < 1e-12);
        assert!((ghost.density - 500.0).abs() < 1e-12);
    }

    // ── wave_height_profile ──────────────────────────────────────────────────

    #[test]
    fn test_wave_height_profile_length() {
        let particles: Vec<FreeSurfaceParticle> = vec![];
        let h = wave_height_profile(&particles, (0.0, 1.0), 10);
        assert_eq!(h.len(), 10);
    }

    #[test]
    fn test_wave_height_profile_zero_bins() {
        let particles: Vec<FreeSurfaceParticle> = vec![];
        let h = wave_height_profile(&particles, (0.0, 1.0), 0);
        assert!(h.is_empty());
    }

    #[test]
    fn test_wave_height_profile_values() {
        let mut p = FreeSurfaceParticle::new([0.5, 2.0, 0.0], 1.0, 1000.0);
        p.is_free_surface = true;
        let particles = vec![p];
        let h = wave_height_profile(&particles, (0.0, 1.0), 5);
        // bin 2 (x=0.5 in [0,1] with 5 bins → bin index 2) should have height 2.0
        assert!((h[2] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_wave_height_profile_non_free_ignored() {
        let p = FreeSurfaceParticle::new([0.5, 5.0, 0.0], 1.0, 1000.0);
        // not marked as free surface
        let particles = vec![p];
        let h = wave_height_profile(&particles, (0.0, 1.0), 5);
        assert_eq!(h[2], 0.0);
    }

    // ── wave_speed_estimate ──────────────────────────────────────────────────

    #[test]
    fn test_wave_speed_empty() {
        assert_eq!(wave_speed_estimate(&[], 0.1, 9.81), 0.0);
    }

    #[test]
    fn test_wave_speed_positive() {
        let h = vec![1.0, 1.0, 1.0];
        let c = wave_speed_estimate(&h, 0.1, 9.81);
        assert!((c - 9.81_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_wave_speed_scales_with_g() {
        let h = vec![1.0; 4];
        let c1 = wave_speed_estimate(&h, 0.1, 9.81);
        let c2 = wave_speed_estimate(&h, 0.1, 4.0 * 9.81);
        assert!((c2 / c1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_wave_speed_zero_height() {
        let h = vec![0.0; 3];
        assert_eq!(wave_speed_estimate(&h, 0.1, 9.81), 0.0);
    }

    // ── splash_detection ─────────────────────────────────────────────────────

    #[test]
    fn test_splash_detection_empty() {
        let particles: Vec<FreeSurfaceParticle> = vec![];
        let s = splash_detection(&particles, 0.5);
        assert!(s.is_empty());
    }

    #[test]
    fn test_splash_detection_none_airborne() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        let particles = vec![p];
        let s = splash_detection(&particles, 0.5);
        assert!(s.is_empty());
    }

    #[test]
    fn test_splash_detection_airborne() {
        let mut p = FreeSurfaceParticle::new([0.0, 1.0, 0.0], 1.0, 1000.0);
        p.is_free_surface = true;
        let particles = vec![p];
        let s = splash_detection(&particles, 0.5);
        assert_eq!(s, vec![0]);
    }

    #[test]
    fn test_splash_detection_mixed() {
        let mut p1 = FreeSurfaceParticle::new([0.0, 2.0, 0.0], 1.0, 1000.0);
        p1.is_free_surface = true;
        let mut p2 = FreeSurfaceParticle::new([0.0, 0.1, 0.0], 1.0, 1000.0);
        p2.is_free_surface = true;
        let p3 = FreeSurfaceParticle::new([0.0, 3.0, 0.0], 1.0, 1000.0);
        // p3 is not free surface
        let particles = vec![p1, p2, p3];
        let s = splash_detection(&particles, 0.5);
        // Only p1 is free-surface AND above 0.5
        assert_eq!(s, vec![0]);
    }

    // ── volume_conservation_error ────────────────────────────────────────────

    #[test]
    fn test_volume_error_zero() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        // V = 1/1000 = 0.001
        let err = volume_conservation_error(&[p], 0.001);
        assert!(err.abs() < 1e-10);
    }

    #[test]
    fn test_volume_error_positive() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 2.0, 1000.0);
        // V = 0.002, initial = 0.001 → err = +1.0
        let err = volume_conservation_error(&[p], 0.001);
        assert!((err - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_error_zero_initial() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        let err = volume_conservation_error(&[p], 0.0);
        assert_eq!(err, 0.0);
    }

    #[test]
    fn test_volume_conservation_empty() {
        let err = volume_conservation_error(&[], 0.001);
        assert!(err < 0.0 || err == -1.0);
    }

    // ── free_surface_normal_estimation ──────────────────────────────────────

    #[test]
    fn test_normal_estimation_returns_3d() {
        let particles: Vec<FreeSurfaceParticle> = (0..4)
            .map(|i| FreeSurfaceParticle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 1000.0))
            .collect();
        let n = free_surface_normal_estimation(&particles, 0, 0.5);
        assert!(n[0].is_finite() && n[1].is_finite() && n[2].is_finite());
    }

    #[test]
    fn test_normal_estimation_single_particle() {
        let p = FreeSurfaceParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0);
        let n = free_surface_normal_estimation(&[p], 0, 0.1);
        // Self-contribution: gradient of self = 0 (dW at r=0 = 0)
        assert!(n[0].is_finite());
    }
}
