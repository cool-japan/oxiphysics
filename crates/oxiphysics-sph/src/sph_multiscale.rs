#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale SPH methods.
//!
//! Implements:
//! - [`MultiScaleSph`] — particle sets at multiple resolution levels
//! - \[`coarsen_region()`\] — merge fine particles into one coarse particle
//! - \[`refine_region()`\] — split a coarse particle into fine particles
//! - \[`interpolation_zone()`\] — identify the overlap region for scale coupling
//! - \[`transition_kernel()`\] — smooth blending weight between scales
//! - \[`energy_conservation_check()`\] — verify kinetic energy is conserved across a scale change

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── Math helpers ─────────────────────────────────────────────────────────────

/// Euclidean distance between two 3-D points.
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Dot product of two 3-D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Squared magnitude of a 3-D vector.
fn sq_mag3(v: [f64; 3]) -> f64 {
    dot3(v, v)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 1  Particle
// ═══════════════════════════════════════════════════════════════════════════

/// A single SPH particle carrying position, velocity, mass, and smoothing length.
#[derive(Debug, Clone, PartialEq)]
pub struct Particle {
    /// Position \[x, y, z\].
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\].
    pub vel: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Smoothing length.
    pub h: f64,
    /// Particle density.
    pub density: f64,
}

impl Particle {
    /// Create a new particle.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64, h: f64, density: f64) -> Self {
        Self {
            pos,
            vel,
            mass,
            h,
            density,
        }
    }

    /// Kinetic energy `½ m |v|²`.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * sq_mag3(self.vel)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 2  MultiScaleSph
// ═══════════════════════════════════════════════════════════════════════════

/// Multiscale SPH container holding both fine and coarse particle sets.
///
/// Fine particles (small `h`) resolve high-gradient regions; coarse
/// particles (large `h`) cover the bulk of the domain cheaply.
#[derive(Debug, Clone)]
pub struct MultiScaleSph {
    /// Fine-resolution particle set.
    pub fine: Vec<Particle>,
    /// Coarse-resolution particle set.
    pub coarse: Vec<Particle>,
    /// Ratio of coarse smoothing length to fine smoothing length.
    pub scale_ratio: f64,
    /// Width of the transition / overlap zone.
    pub overlap_width: f64,
}

impl MultiScaleSph {
    /// Create a new [`MultiScaleSph`] container.
    ///
    /// # Arguments
    /// * `fine` — fine-resolution particles
    /// * `coarse` — coarse-resolution particles
    /// * `scale_ratio` — `h_coarse / h_fine` (typically 2–4)
    /// * `overlap_width` — width of the blending zone
    pub fn new(
        fine: Vec<Particle>,
        coarse: Vec<Particle>,
        scale_ratio: f64,
        overlap_width: f64,
    ) -> Self {
        Self {
            fine,
            coarse,
            scale_ratio,
            overlap_width,
        }
    }

    /// Total kinetic energy of all particles in both levels.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.fine
            .iter()
            .chain(self.coarse.iter())
            .map(|p| p.kinetic_energy())
            .sum()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 3  coarsen_region
// ═══════════════════════════════════════════════════════════════════════════

/// Merge all fine particles within `radius` of `center` into a single coarse
/// particle by conserving total mass, momentum, and centre-of-mass position.
///
/// Returns the merged coarse particle, or `None` if no fine particles lie
/// within `radius`.
///
/// The original fine particles are **not** removed from any list; the caller
/// is responsible for bookkeeping.
pub fn coarsen_region(fine: &[Particle], center: [f64; 3], radius: f64) -> Option<Particle> {
    let candidates: Vec<&Particle> = fine
        .iter()
        .filter(|p| dist3(p.pos, center) <= radius)
        .collect();
    if candidates.is_empty() {
        return None;
    }

    let total_mass: f64 = candidates.iter().map(|p| p.mass).sum();
    let mut com = [0.0_f64; 3];
    let mut mom = [0.0_f64; 3];
    for p in &candidates {
        for d in 0..3 {
            com[d] += p.mass * p.pos[d];
            mom[d] += p.mass * p.vel[d];
        }
    }
    for d in 0..3 {
        com[d] /= total_mass;
        mom[d] /= total_mass; // momentum → velocity
    }
    // Coarse smoothing length: scale with cube root of number of particles merged
    let h_coarse = candidates[0].h * (candidates.len() as f64).cbrt();
    Some(Particle::new(
        com,
        mom,
        total_mass,
        h_coarse,
        total_mass / (4.0 / 3.0 * std::f64::consts::PI * h_coarse.powi(3)),
    ))
}

// ═══════════════════════════════════════════════════════════════════════════
// § 4  refine_region
// ═══════════════════════════════════════════════════════════════════════════

/// Split a coarse particle into `n_fine` fine particles arranged on a regular
/// lattice around the coarse particle's position.
///
/// Mass and momentum are conserved: each child gets `mass / n_fine` and the
/// same velocity as the parent.  The fine smoothing length is
/// `h_coarse / scale_ratio`.
///
/// `n_fine` is clamped to 1 if zero or negative.
pub fn refine_region(coarse: &Particle, n_fine: usize, scale_ratio: f64) -> Vec<Particle> {
    let n = n_fine.max(1);
    let m_child = coarse.mass / n as f64;
    let h_child = coarse.h / scale_ratio;
    let spacing = h_child * 0.5;

    let side = (n as f64).cbrt().ceil() as usize;
    let mut children = Vec::with_capacity(n);
    let mut count = 0;
    'outer: for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                if count >= n {
                    break 'outer;
                }
                let offset = [
                    (ix as f64 - side as f64 * 0.5) * spacing,
                    (iy as f64 - side as f64 * 0.5) * spacing,
                    (iz as f64 - side as f64 * 0.5) * spacing,
                ];
                let pos = [
                    coarse.pos[0] + offset[0],
                    coarse.pos[1] + offset[1],
                    coarse.pos[2] + offset[2],
                ];
                children.push(Particle::new(
                    pos,
                    coarse.vel,
                    m_child,
                    h_child,
                    coarse.density,
                ));
                count += 1;
            }
        }
    }
    children
}

// ═══════════════════════════════════════════════════════════════════════════
// § 5  interpolation_zone
// ═══════════════════════════════════════════════════════════════════════════

/// Return the indices of particles (from `particles`) that lie within the
/// transition / overlap zone defined by `[inner_radius, outer_radius]` from
/// `center`.
///
/// Particles inside `inner_radius` are pure fine-scale; particles outside
/// `outer_radius` are pure coarse-scale; those in between form the blending zone.
pub fn interpolation_zone(
    particles: &[Particle],
    center: [f64; 3],
    inner_radius: f64,
    outer_radius: f64,
) -> Vec<usize> {
    particles
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            let d = dist3(p.pos, center);
            d >= inner_radius && d <= outer_radius
        })
        .map(|(i, _)| i)
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// § 6  transition_kernel
// ═══════════════════════════════════════════════════════════════════════════

/// Smooth blending weight ∈ \[0, 1\] for a particle at distance `r` from the
/// fine-scale region boundary.
///
/// Uses a quintic polynomial that has zero first and second derivatives at
/// both ends of the transition zone `[r_inner, r_outer]`:
///
/// ```text
/// t = (r − r_inner) / (r_outer − r_inner)
/// w = 1 − 10 t³ + 15 t⁴ − 6 t⁵
/// ```
///
/// Returns 1 in the fine region, 0 in the coarse region.
pub fn transition_kernel(r: f64, r_inner: f64, r_outer: f64) -> f64 {
    if r <= r_inner {
        return 1.0;
    }
    if r >= r_outer {
        return 0.0;
    }
    let t = (r - r_inner) / (r_outer - r_inner);
    1.0 - 10.0 * t.powi(3) + 15.0 * t.powi(4) - 6.0 * t.powi(5)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 7  energy_conservation_check
// ═══════════════════════════════════════════════════════════════════════════

/// Verify that kinetic energy is approximately conserved after a coarsening or
/// refinement operation.
///
/// Returns `(energy_before, energy_after, relative_error)`.  A relative error
/// below `tol` indicates good conservation.
pub fn energy_conservation_check(
    particles_before: &[Particle],
    particles_after: &[Particle],
) -> (f64, f64, f64) {
    let ke_before: f64 = particles_before.iter().map(|p| p.kinetic_energy()).sum();
    let ke_after: f64 = particles_after.iter().map(|p| p.kinetic_energy()).sum();
    let rel = if ke_before.abs() > 1e-300 {
        (ke_after - ke_before).abs() / ke_before
    } else {
        (ke_after - ke_before).abs()
    };
    (ke_before, ke_after, rel)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 8  Additional helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Interpolate the velocity at position `r` using SPH kernel-weighted average
/// of the supplied `particles`.
///
/// Uses the standard cubic-spline kernel in 3D.
pub fn interpolate_velocity(r: [f64; 3], particles: &[Particle]) -> [f64; 3] {
    let mut weighted_vel = [0.0_f64; 3];
    let mut weight_sum = 0.0_f64;
    for p in particles {
        let d = dist3(r, p.pos);
        let q = d / (p.h + 1e-300);
        let w = cubic_spline_kernel(q);
        let vol = p.mass / (p.density + 1e-300);
        let wv = w * vol;
        for d in 0..3 {
            weighted_vel[d] += wv * p.vel[d];
        }
        weight_sum += wv;
    }
    if weight_sum > 1e-300 {
        for d in 0..3 {
            weighted_vel[d] /= weight_sum;
        }
    }
    weighted_vel
}

/// Cubic-spline SPH kernel value for normalized distance `q = r/h`.
fn cubic_spline_kernel(q: f64) -> f64 {
    let sigma = 1.0 / std::f64::consts::PI;
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── helpers ──────────────────────────────────────────────────────────

    fn make_particle(pos: [f64; 3], vel: [f64; 3], mass: f64, h: f64) -> Particle {
        Particle::new(
            pos,
            vel,
            mass,
            h,
            mass / (4.0 / 3.0 * std::f64::consts::PI * h.powi(3)),
        )
    }

    // ── Particle tests ────────────────────────────────────────────────────

    #[test]
    fn test_particle_kinetic_energy_zero_vel() {
        let p = make_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1);
        assert!(p.kinetic_energy() < 1e-15);
    }

    #[test]
    fn test_particle_kinetic_energy_nonzero() {
        // KE = 0.5 * 2.0 * (3^2+4^2+0^2) = 25
        let p = make_particle([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 2.0, 0.1);
        assert!((p.kinetic_energy() - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_particle_clone() {
        let p = make_particle([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 0.5, 0.05);
        let q = p.clone();
        assert_eq!(p, q);
    }

    // ── MultiScaleSph tests ───────────────────────────────────────────────

    #[test]
    fn test_multiscale_total_ke() {
        let fine = vec![make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.1)];
        let coarse = vec![make_particle([10.0, 0.0, 0.0], [2.0, 0.0, 0.0], 4.0, 0.4)];
        let ms = MultiScaleSph::new(fine, coarse, 4.0, 0.2);
        // KE_fine = 0.5, KE_coarse = 8.0
        assert!((ms.total_kinetic_energy() - 8.5).abs() < 1e-10);
    }

    #[test]
    fn test_multiscale_new_stores_fields() {
        let ms = MultiScaleSph::new(vec![], vec![], 2.0, 0.5);
        assert_eq!(ms.scale_ratio, 2.0);
        assert_eq!(ms.overlap_width, 0.5);
    }

    // ── coarsen_region tests ──────────────────────────────────────────────

    #[test]
    fn test_coarsen_empty_region() {
        let fine = vec![make_particle([100.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1)];
        let result = coarsen_region(&fine, [0.0, 0.0, 0.0], 0.5);
        assert!(result.is_none());
    }

    #[test]
    fn test_coarsen_single_particle() {
        let fine = vec![make_particle([0.0, 0.0, 0.0], [1.0, 2.0, 3.0], 1.0, 0.1)];
        let result = coarsen_region(&fine, [0.0, 0.0, 0.0], 1.0).unwrap();
        assert!((result.mass - 1.0).abs() < 1e-10);
        assert!((result.vel[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_coarsen_conserves_mass() {
        let fine = vec![
            make_particle([0.1, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.1),
            make_particle([-0.1, 0.0, 0.0], [3.0, 0.0, 0.0], 3.0, 0.1),
        ];
        let coarse = coarsen_region(&fine, [0.0, 0.0, 0.0], 0.5).unwrap();
        assert!((coarse.mass - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_coarsen_conserves_momentum() {
        let fine = vec![
            make_particle([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], 1.0, 0.1),
            make_particle([0.2, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.1),
        ];
        let coarse = coarsen_region(&fine, [0.1, 0.0, 0.0], 0.5).unwrap();
        // total momentum: 4+2=6, total mass: 2 → vel = 3
        assert!((coarse.vel[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_coarsen_centre_of_mass() {
        let fine = vec![
            make_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1),
            make_particle([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1),
        ];
        let coarse = coarsen_region(&fine, [1.0, 0.0, 0.0], 2.0).unwrap();
        assert!((coarse.pos[0] - 1.0).abs() < 1e-10);
    }

    // ── refine_region tests ───────────────────────────────────────────────

    #[test]
    fn test_refine_conserves_mass() {
        let coarse = make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 8.0, 0.4);
        let children = refine_region(&coarse, 8, 2.0);
        let total_mass: f64 = children.iter().map(|p| p.mass).sum();
        assert!((total_mass - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_refine_child_count() {
        let coarse = make_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 4.0, 0.4);
        let children = refine_region(&coarse, 4, 2.0);
        assert_eq!(children.len(), 4);
    }

    #[test]
    fn test_refine_child_smoothing_length() {
        let coarse = make_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.4);
        let children = refine_region(&coarse, 2, 2.0);
        for c in &children {
            assert!((c.h - 0.2).abs() < 1e-10);
        }
    }

    #[test]
    fn test_refine_velocity_inherited() {
        let coarse = make_particle([0.0, 0.0, 0.0], [5.0, -3.0, 2.0], 1.0, 0.1);
        let children = refine_region(&coarse, 3, 2.0);
        for c in &children {
            assert!((c.vel[0] - 5.0).abs() < 1e-10);
            assert!((c.vel[1] + 3.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_refine_n_zero_gives_one_child() {
        let coarse = make_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1);
        let children = refine_region(&coarse, 0, 2.0);
        assert_eq!(children.len(), 1);
    }

    // ── interpolation_zone tests ──────────────────────────────────────────

    #[test]
    fn test_interpolation_zone_empty() {
        let particles = vec![make_particle([10.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1)];
        let zone = interpolation_zone(&particles, [0.0, 0.0, 0.0], 0.5, 1.0);
        assert!(zone.is_empty());
    }

    #[test]
    fn test_interpolation_zone_all_inside() {
        // All particles at distance 0.75, inner=0.5, outer=1.0 → all in zone
        let particles = vec![
            make_particle([0.75, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1),
            make_particle([-0.75, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1),
        ];
        let zone = interpolation_zone(&particles, [0.0, 0.0, 0.0], 0.5, 1.0);
        assert_eq!(zone.len(), 2);
    }

    #[test]
    fn test_interpolation_zone_excludes_interior() {
        let particles = vec![
            make_particle([0.3, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1), // inside inner
            make_particle([0.75, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1), // in zone
        ];
        let zone = interpolation_zone(&particles, [0.0, 0.0, 0.0], 0.5, 1.0);
        assert_eq!(zone, vec![1]);
    }

    // ── transition_kernel tests ───────────────────────────────────────────

    #[test]
    fn test_transition_kernel_inner() {
        assert!((transition_kernel(0.0, 0.5, 1.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_transition_kernel_outer() {
        assert!(transition_kernel(1.5, 0.5, 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_transition_kernel_midpoint() {
        let w = transition_kernel(0.75, 0.5, 1.0);
        assert!(w > 0.0 && w < 1.0, "expected w in (0,1), got {w}");
    }

    #[test]
    fn test_transition_kernel_monotone() {
        let r_inner = 0.0;
        let r_outer = 1.0;
        let mut prev = 1.0;
        for i in 1..=10 {
            let r = i as f64 * 0.1;
            let w = transition_kernel(r, r_inner, r_outer);
            assert!(w <= prev + 1e-12, "not monotone at r={r}");
            prev = w;
        }
    }

    // ── energy_conservation_check tests ──────────────────────────────────

    #[test]
    fn test_energy_conservation_identical() {
        let ps = vec![make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.1)];
        let (_kb, _ka, rel) = energy_conservation_check(&ps, &ps);
        assert!(rel < 1e-12);
    }

    #[test]
    fn test_energy_conservation_nonzero_error() {
        let before = vec![make_particle([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.1)];
        let after = vec![make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.1)];
        let (kb, ka, rel) = energy_conservation_check(&before, &after);
        assert!(kb > ka);
        assert!(rel > 0.0);
    }

    #[test]
    fn test_energy_conservation_after_refinement() {
        let coarse = make_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 0.4);
        let ke_before = coarse.kinetic_energy(); // 0.5 * 4 * 1 = 2.0
        let children = refine_region(&coarse, 4, 2.0);
        let ke_after: f64 = children.iter().map(|p| p.kinetic_energy()).sum();
        let rel_err = (ke_after - ke_before).abs() / ke_before;
        assert!(rel_err < 1e-10, "energy rel error: {rel_err}");
    }

    // ── cubic_spline_kernel tests ─────────────────────────────────────────

    #[test]
    fn test_cubic_spline_kernel_zero() {
        let w = cubic_spline_kernel(0.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_spline_kernel_beyond_support() {
        assert!(cubic_spline_kernel(2.5).abs() < 1e-15);
    }

    #[test]
    fn test_cubic_spline_kernel_monotone_decrease() {
        let mut prev = cubic_spline_kernel(0.0);
        for i in 1..=20 {
            let q = i as f64 * 0.1;
            let w = cubic_spline_kernel(q);
            assert!(w <= prev + 1e-12, "not monotone at q={q}");
            prev = w;
        }
    }

    // ── interpolate_velocity tests ────────────────────────────────────────

    #[test]
    fn test_interpolate_velocity_single_particle() {
        let p = make_particle([0.0, 0.0, 0.0], [1.0, 2.0, 3.0], 1.0, 1.0);
        let vel = interpolate_velocity([0.0, 0.0, 0.0], &[p]);
        // Should recover the particle velocity (approximately)
        assert!((vel[0] - 1.0).abs() < 0.01, "vx={}", vel[0]);
    }

    #[test]
    fn test_interpolate_velocity_empty() {
        let vel = interpolate_velocity([0.0, 0.0, 0.0], &[]);
        assert_eq!(vel, [0.0, 0.0, 0.0]);
    }

    // ── dist3 helper tests ────────────────────────────────────────────────

    #[test]
    fn test_dist3_same_point() {
        assert!(dist3([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) < 1e-15);
    }

    #[test]
    fn test_dist3_known_distance() {
        assert!((dist3([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]) - 5.0).abs() < 1e-10);
    }
}
