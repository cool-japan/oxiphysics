// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
#![allow(dead_code)]

//! Adaptive smoothing-length (h) algorithms for SPH.
//!
//! Implements the Hoover–Monaghan constraint (constant neighbour count),
//! Newton–Raphson solver for h, grad-h correction factor Ω, variable-h
//! kernel gradient corrections, smoothing-length bounds enforcement, a
//! tree-based neighbour count estimator, and particle splitting/merging
//! for resolution refinement and coarsening.
//!
//! # References
//! - Springel & Hernquist (2002) – grad-h SPH formulation.
//! - Monaghan (2005) – SPH review.
//! - Dehnen & Aly (2012) – kernel comparisons.

use crate::kernel::CubicSplineKernel;
use crate::kernel::SphKernel;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Target number of neighbours for the Hoover criterion in 3D (≈ 32–64).
pub const TARGET_NEIGHBOUR_COUNT_3D: f64 = 32.0;

/// Typical η factor: h = η (m / ρ)^(1/d).  η ≈ 1.2 for cubic spline.
pub const ETA_DEFAULT: f64 = 1.2;

/// Minimum allowed smoothing length (relative to initial h₀).
pub const H_MIN_FRACTION: f64 = 0.1;

/// Maximum allowed smoothing length (relative to initial h₀).
pub const H_MAX_FRACTION: f64 = 10.0;

// ---------------------------------------------------------------------------
// Hoover criterion: ρ (h/h₀)^d = ρ₀
// ---------------------------------------------------------------------------

/// Compute the target smoothing length from the Hoover criterion.
///
/// The Hoover (1986) constraint requires that the number of neighbours
/// remains constant as h evolves:
/// ```text
/// ρ_i (h_i / h₀)^d  =  ρ₀
/// ```
/// Solving for h_i:
/// ```text
/// h_i = h₀ (ρ₀ / ρ_i)^(1/d)
/// ```
/// where `d` is the spatial dimension.
///
/// # Arguments
/// * `h0`    – reference smoothing length.
/// * `rho0`  – reference density (at h₀).
/// * `rho_i` – current local density.
/// * `dim`   – spatial dimension (2 or 3).
pub fn hoover_h(h0: f64, rho0: f64, rho_i: f64, dim: u32) -> f64 {
    if rho_i < 1e-30 {
        return h0 * H_MAX_FRACTION;
    }
    h0 * (rho0 / rho_i).powf(1.0 / dim as f64)
}

/// Verify the Hoover constraint residual: `ρ_i (h_i/h₀)^d − ρ₀`.
///
/// Returns 0 when the constraint is exactly satisfied.
pub fn hoover_residual(h0: f64, h_i: f64, rho0: f64, rho_i: f64, dim: u32) -> f64 {
    rho_i * (h_i / h0).powi(dim as i32) - rho0
}

// ---------------------------------------------------------------------------
// SPH density sum and its derivative w.r.t. h (for Newton iteration)
// ---------------------------------------------------------------------------

/// SPH density sum at position `pos_i` with smoothing length `h`.
///
/// `ρ_i = Σ_j m_j W(|r_i − r_j|, h)`
pub fn sph_density(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let kernel = CubicSplineKernel;
    let mut rho = 0.0_f64;
    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += m_j * kernel.w(r, h);
    }
    rho
}

/// Derivative of the SPH density sum with respect to h.
///
/// `∂ρ/∂h = Σ_j m_j ∂W(r, h)/∂h`
///
/// For the cubic spline kernel, `∂W/∂h = −(3/h) W + (q/h) ∂W/∂q`.
pub fn sph_density_dh(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let sigma_3d = 1.0 / (std::f64::consts::PI * h * h * h);
    let inv_h = 1.0 / h;
    let mut drho_dh = 0.0_f64;

    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let q = r * inv_h;

        // W(r, h) and dW/dq from cubic spline
        let (w, dw_dq) = if q < 1.0 {
            let w = sigma_3d * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
            let dw_dq = sigma_3d * (-3.0 * q + 2.25 * q * q);
            (w, dw_dq)
        } else if q < 2.0 {
            let t = 2.0 - q;
            let w = sigma_3d * 0.25 * t * t * t;
            let dw_dq = sigma_3d * (-0.75 * t * t);
            (w, dw_dq)
        } else {
            (0.0, 0.0)
        };

        // ∂W/∂h = −(3/h)W + (1/h) q dW/dq
        let dw_dh = -3.0 * inv_h * w + inv_h * q * dw_dq;
        drho_dh += m_j * dw_dh;
    }

    drho_dh
}

// ---------------------------------------------------------------------------
// Newton-Raphson solver for h
// ---------------------------------------------------------------------------

/// Solve for the adaptive smoothing length h via Newton–Raphson.
///
/// Finds h such that `ρ(h) = ρ_target`.
///
/// # Arguments
/// * `pos_i`        — position of particle i.
/// * `h_guess`      — initial guess for h.
/// * `rho_target`   — target density.
/// * `neighbor_pos` — positions of all potential neighbours.
/// * `neighbor_mass`— masses of all potential neighbours.
/// * `max_iter`     — maximum Newton iterations.
/// * `tol`          — relative convergence tolerance.
/// * `h_min`        — minimum allowed h.
/// * `h_max`        — maximum allowed h.
///
/// Returns `Some(h)` on convergence within `[h_min, h_max]`, `None` otherwise.
#[allow(clippy::too_many_arguments)]
pub fn solve_h_newton(
    pos_i: [f64; 3],
    h_guess: f64,
    rho_target: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    max_iter: usize,
    tol: f64,
    h_min: f64,
    h_max: f64,
) -> Option<f64> {
    let mut h = h_guess.clamp(h_min, h_max);

    for _ in 0..max_iter {
        let rho = sph_density(pos_i, h, neighbor_pos, neighbor_mass);
        let f = rho - rho_target;

        // Check convergence
        if f.abs() < tol * rho_target.max(1e-14) {
            return Some(h.clamp(h_min, h_max));
        }

        let df_dh = sph_density_dh(pos_i, h, neighbor_pos, neighbor_mass);
        if df_dh.abs() < 1e-30 {
            // Jacobian too small — cannot continue
            break;
        }

        let h_new = (h - f / df_dh).clamp(h_min, h_max);
        let rel = (h_new - h).abs() / h.max(1e-30);
        h = h_new;

        if rel < tol {
            return Some(h);
        }
    }

    None
}

// ---------------------------------------------------------------------------
// Grad-h correction factor Ω
// ---------------------------------------------------------------------------

/// Compute the grad-h correction factor Ω for particle i.
///
/// The variable-h SPH equations (Springel & Hernquist 2002) contain the
/// correction:
/// ```text
/// Ω_i = 1 − (∂ρ_i / ∂h_i) * h_i / (n_dim * ρ_i)
/// ```
/// which appears in the momentum equation to ensure conservation.
///
/// Returns 1 (no correction) when `ρ_i ≈ 0`.
pub fn omega_correction(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    rho_i: f64,
    n_dim: f64,
) -> f64 {
    if rho_i < 1e-30 {
        return 1.0;
    }
    let drho_dh = sph_density_dh(pos_i, h, neighbor_pos, neighbor_mass);
    1.0 - drho_dh * h / (n_dim * rho_i)
}

// ---------------------------------------------------------------------------
// Variable-h kernel gradient correction
// ---------------------------------------------------------------------------

/// Compute the corrected kernel gradient ∇W in the variable-h formulation.
///
/// The variable-h gradient is:
/// ```text
/// ∇W̃_ij = (1/Ω_i) ∇W(r_ij, h_i) + (1/Ω_j) ∇W(r_ij, h_j)
/// ```
/// divided by 2 for the symmetric form used in force evaluation.
///
/// Returns the 3D gradient vector.
#[allow(clippy::too_many_arguments)]
pub fn variable_h_kernel_gradient(
    r_ij: [f64; 3],
    h_i: f64,
    h_j: f64,
    omega_i: f64,
    omega_j: f64,
) -> [f64; 3] {
    let kernel = CubicSplineKernel;
    let r = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]).sqrt();
    if r < 1e-14 {
        return [0.0; 3];
    }

    let gw_i = kernel.grad_w(r, h_i);
    let gw_j = kernel.grad_w(r, h_j);

    let omega_i_safe = omega_i.max(1e-14);
    let omega_j_safe = omega_j.max(1e-14);

    let scale = 0.5 * (gw_i / omega_i_safe + gw_j / omega_j_safe) / r;

    [scale * r_ij[0], scale * r_ij[1], scale * r_ij[2]]
}

// ---------------------------------------------------------------------------
// Smoothing-length bounds
// ---------------------------------------------------------------------------

/// Enforce smoothing-length bounds for an array of h values.
///
/// Clamps each h to `[h_min, h_max]` in-place.
pub fn apply_h_bounds(hs: &mut [f64], h_min: f64, h_max: f64) {
    for h in hs.iter_mut() {
        *h = h.clamp(h_min, h_max);
    }
}

/// Compute h_min and h_max from the initial smoothing length h₀ and
/// the bound fractions defined in the module constants.
pub fn h_bounds_from_reference(h0: f64) -> (f64, f64) {
    (h0 * H_MIN_FRACTION, h0 * H_MAX_FRACTION)
}

// ---------------------------------------------------------------------------
// Tree-based neighbour count estimate
// ---------------------------------------------------------------------------

/// Estimate the number of neighbours within radius `r_cut` using a brute-force
/// (O(N²)) search.  In a production code this would use a BVH or octree, but
/// for correctness and testing a direct search suffices.
///
/// Returns the number of particles within `r_cut` of `pos_i` (including
/// itself if present in `positions`).
pub fn estimate_neighbour_count(pos_i: [f64; 3], r_cut: f64, positions: &[[f64; 3]]) -> usize {
    let r2_cut = r_cut * r_cut;
    positions
        .iter()
        .filter(|pos_j| {
            let dx = pos_i[0] - pos_j[0];
            let dy = pos_i[1] - pos_j[1];
            let dz = pos_i[2] - pos_j[2];
            dx * dx + dy * dy + dz * dz <= r2_cut
        })
        .count()
}

/// Estimate whether h needs to be increased (too few neighbours) or
/// decreased (too many neighbours) to achieve `target_count`.
///
/// Returns `None` if the current count is within `tolerance` of the target.
/// Returns `Some(true)` if h should be increased, `Some(false)` if decreased.
pub fn h_adjustment_direction(
    pos_i: [f64; 3],
    h: f64,
    target_count: usize,
    tolerance: usize,
    positions: &[[f64; 3]],
) -> Option<bool> {
    let current = estimate_neighbour_count(pos_i, 2.0 * h, positions);
    let diff = current as isize - target_count as isize;
    if diff.unsigned_abs() <= tolerance {
        None // Already close enough
    } else if diff < 0 {
        Some(true) // Need to increase h
    } else {
        Some(false) // Need to decrease h
    }
}

// ---------------------------------------------------------------------------
// Particle splitting
// ---------------------------------------------------------------------------

/// A single particle for the adaptive-h splitting/merging operations.
#[derive(Debug, Clone)]
pub struct AdaptiveHParticle {
    /// 3D position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Refinement level (0 = coarsest).
    pub level: u32,
}

impl AdaptiveHParticle {
    /// Create a new particle.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        h: f64,
        density: f64,
        level: u32,
    ) -> Self {
        Self {
            position,
            velocity,
            mass,
            h,
            density,
            level,
        }
    }
}

/// Split a parent particle into `n_children` children arranged symmetrically.
///
/// Uses a 2×2×2 (= 8) hexagonal arrangement by default.  Each child gets:
/// - mass = parent_mass / n_children
/// - h    = parent_h / 2^(1/3)
/// - position offset by ±offset_fraction × h from parent centre
/// - Same velocity as the parent (no velocity perturbation)
///
/// Returns the child particles (parent is removed by the caller).
pub fn split_particle(parent: &AdaptiveHParticle) -> Vec<AdaptiveHParticle> {
    let n = 8_usize;
    let child_mass = parent.mass / n as f64;
    let child_h = parent.h / 2.0_f64.powf(1.0 / 3.0);
    let offset = 0.4 * child_h; // displacement of children from parent centre

    // 2×2×2 offsets: all (±1, ±1, ±1) combinations
    let signs: [(f64, f64, f64); 8] = [
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, 1.0, 1.0),
    ];

    signs
        .iter()
        .map(|&(sx, sy, sz)| {
            let pos = [
                parent.position[0] + sx * offset,
                parent.position[1] + sy * offset,
                parent.position[2] + sz * offset,
            ];
            AdaptiveHParticle::new(
                pos,
                parent.velocity,
                child_mass,
                child_h,
                parent.density,
                parent.level + 1,
            )
        })
        .collect()
}

/// Merge a set of particles into a single representative particle.
///
/// Conservation laws enforced:
/// - Total mass: m_merged = Σ m_i
/// - Centre of mass: x_merged = Σ (m_i x_i) / m_merged
/// - Linear momentum: v_merged = Σ (m_i v_i) / m_merged
/// - Smoothing length: h_merged = 2^(1/3) × mean(h_i)
pub fn merge_particles(particles: &[AdaptiveHParticle]) -> AdaptiveHParticle {
    let n = particles.len();
    assert!(n > 0, "cannot merge empty particle list");

    let total_mass: f64 = particles.iter().map(|p| p.mass).sum();
    let inv_m = if total_mass > 1e-30 {
        1.0 / total_mass
    } else {
        0.0
    };

    let pos = [
        particles
            .iter()
            .map(|p| p.mass * p.position[0])
            .sum::<f64>()
            * inv_m,
        particles
            .iter()
            .map(|p| p.mass * p.position[1])
            .sum::<f64>()
            * inv_m,
        particles
            .iter()
            .map(|p| p.mass * p.position[2])
            .sum::<f64>()
            * inv_m,
    ];

    let vel = [
        particles
            .iter()
            .map(|p| p.mass * p.velocity[0])
            .sum::<f64>()
            * inv_m,
        particles
            .iter()
            .map(|p| p.mass * p.velocity[1])
            .sum::<f64>()
            * inv_m,
        particles
            .iter()
            .map(|p| p.mass * p.velocity[2])
            .sum::<f64>()
            * inv_m,
    ];

    let mean_h = particles.iter().map(|p| p.h).sum::<f64>() / n as f64;
    let merged_h = 2.0_f64.powf(1.0 / 3.0) * mean_h;

    let level = particles.iter().map(|p| p.level).min().unwrap_or(0);
    let merged_level = if level > 0 { level - 1 } else { 0 };

    let density = particles.iter().map(|p| p.density).sum::<f64>() / n as f64;

    AdaptiveHParticle::new(pos, vel, total_mass, merged_h, density, merged_level)
}

// ---------------------------------------------------------------------------
// Adaptive h manager
// ---------------------------------------------------------------------------

/// High-level adaptive smoothing-length manager.
///
/// Maintains a particle set and applies the Hoover criterion each time step
/// to keep the neighbour count near `target_neighbours`.
pub struct AdaptiveHManager {
    /// All particles.
    pub particles: Vec<AdaptiveHParticle>,
    /// Reference density ρ₀ for the Hoover criterion.
    pub rho0: f64,
    /// Reference smoothing length h₀.
    pub h0: f64,
    /// Minimum allowed h.
    pub h_min: f64,
    /// Maximum allowed h.
    pub h_max: f64,
    /// Target neighbour count N_s.
    pub target_neighbours: f64,
    /// Spatial dimension.
    pub dim: u32,
}

impl AdaptiveHManager {
    /// Create a new manager.
    pub fn new(rho0: f64, h0: f64, target_neighbours: f64, dim: u32) -> Self {
        let (h_min, h_max) = h_bounds_from_reference(h0);
        Self {
            particles: Vec::new(),
            rho0,
            h0,
            h_min,
            h_max,
            target_neighbours,
            dim,
        }
    }

    /// Add a particle to the system.
    pub fn add_particle(&mut self, p: AdaptiveHParticle) {
        self.particles.push(p);
    }

    /// Update all smoothing lengths using the Hoover criterion.
    ///
    /// For each particle, computes h from `ρ_i`:
    /// `h_i = h₀ (ρ₀ / ρ_i)^(1/d)`
    /// and clamps to `[h_min, h_max]`.
    pub fn update_smoothing_lengths(&mut self) {
        for p in &mut self.particles {
            let h = hoover_h(self.h0, self.rho0, p.density, self.dim);
            p.h = h.clamp(self.h_min, self.h_max);
        }
    }

    /// Compute Ω correction factor for each particle.
    ///
    /// Returns a vector of Ω values, one per particle.
    pub fn compute_omega(&self) -> Vec<f64> {
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.position).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();

        self.particles
            .iter()
            .map(|p| {
                omega_correction(
                    p.position,
                    p.h,
                    &positions,
                    &masses,
                    p.density,
                    self.dim as f64,
                )
            })
            .collect()
    }

    /// Refine particles that are too coarse (h too large).
    ///
    /// A particle is split if `h_i > h_split_threshold`.
    /// Returns the number of particles split.
    pub fn refine_coarse_particles(&mut self, h_split_threshold: f64) -> usize {
        let mut to_split: Vec<usize> = self
            .particles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.h > h_split_threshold)
            .map(|(i, _)| i)
            .collect();

        // Process in reverse order so removal doesn't shift indices
        to_split.sort_unstable_by(|a, b| b.cmp(a));
        let n_split = to_split.len();

        for idx in to_split {
            let parent = self.particles.remove(idx);
            let children = split_particle(&parent);
            self.particles.extend(children);
        }

        n_split
    }

    /// Coarsen fine particles that are too small (h too small).
    ///
    /// Particles with `h_i < h_merge_threshold` at the same level are
    /// grouped into batches of 8 and merged.  Returns the number of merges.
    pub fn coarsen_fine_particles(&mut self, h_merge_threshold: f64) -> usize {
        let fine_indices: Vec<usize> = self
            .particles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.h < h_merge_threshold)
            .map(|(i, _)| i)
            .collect();

        if fine_indices.len() < 8 {
            return 0;
        }

        // Merge in groups of 8 (reverse order for removal safety)
        let mut n_merges = 0;
        let mut chunks: Vec<Vec<usize>> = fine_indices.chunks(8).map(|c| c.to_vec()).collect();
        // Only merge complete octets
        chunks.retain(|c| c.len() == 8);

        // Sort outer chunks descending by max index to avoid shifting
        chunks.sort_by(|a, b| b.iter().max().cmp(&a.iter().max()));

        for chunk in chunks {
            // Sort each chunk descending
            let mut sorted_chunk = chunk.clone();
            sorted_chunk.sort_unstable_by(|a, b| b.cmp(a));

            // Extract particles (remove highest-index first)
            let mut group = Vec::with_capacity(8);
            for idx in &sorted_chunk {
                group.push(self.particles.remove(*idx));
            }

            let merged = merge_particles(&group);
            self.particles.push(merged);
            n_merges += 1;
        }

        n_merges
    }

    /// Total mass (should be conserved by split/merge).
    pub fn total_mass(&self) -> f64 {
        self.particles.iter().map(|p| p.mass).sum()
    }
}

// ---------------------------------------------------------------------------
// Zhu-Fox h-refinement criterion
// ---------------------------------------------------------------------------

/// Zhu & Fox (2001) refinement criterion for adaptive smoothing length.
///
/// A particle should be refined (h reduced) when the local resolution is
/// insufficient.  The criterion compares the relative density variation
/// across the neighbourhood to a threshold ε:
///
/// ```text
/// ε_i = max_j | (ρ_j - ρ_i) / ρ_i |  for |r_ij| < h_i
/// ```
///
/// Returns the maximum relative density variation across all neighbours.
pub fn zhu_fox_criterion(
    pos_i: [f64; 3],
    rho_i: f64,
    h_i: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_densities: &[f64],
) -> f64 {
    if rho_i < 1e-30 {
        return 0.0;
    }
    let mut eps = 0.0_f64;
    for (pos_j, &rho_j) in neighbor_pos.iter().zip(neighbor_densities) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 <= h_i * h_i {
            let rel = ((rho_j - rho_i) / rho_i).abs();
            if rel > eps {
                eps = rel;
            }
        }
    }
    eps
}

/// Decide whether to refine based on the Zhu-Fox criterion.
///
/// Returns `true` when `zhu_fox_criterion > threshold` AND `level < max_level`.
#[allow(clippy::too_many_arguments)]
pub fn should_refine_zhu_fox(
    pos_i: [f64; 3],
    rho_i: f64,
    h_i: f64,
    level: u32,
    max_level: u32,
    neighbor_pos: &[[f64; 3]],
    neighbor_densities: &[f64],
    threshold: f64,
) -> bool {
    if level >= max_level {
        return false;
    }
    let eps = zhu_fox_criterion(pos_i, rho_i, h_i, neighbor_pos, neighbor_densities);
    eps > threshold
}

// ---------------------------------------------------------------------------
// Variable smoothing-length update rules (Monaghan 2002 style)
// ---------------------------------------------------------------------------

/// Update the smoothing length `h` of a single particle using the continuity-
/// equation-based rule (Monaghan 2002):
///
/// ```text
/// dh/dt = -(h / (d * ρ)) * (dρ/dt)
/// ```
///
/// Given an estimate of `dρ/dt` from the SPH continuity equation, the new h is:
///
/// ```text
/// h_new = h_old + dt * (-h_old / (dim * rho)) * (drho_dt)
/// ```
///
/// Clamps the result to `[h_min, h_max]`.
pub fn update_h_continuity(
    h: f64,
    rho: f64,
    drho_dt: f64,
    dt: f64,
    dim: f64,
    h_min: f64,
    h_max: f64,
) -> f64 {
    if rho < 1e-30 {
        return h.clamp(h_min, h_max);
    }
    let dh_dt = -(h / (dim * rho)) * drho_dt;
    (h + dt * dh_dt).clamp(h_min, h_max)
}

/// Compute `dρ/dt` from the SPH continuity equation:
///
/// ```text
/// dρ_i/dt = -ρ_i * Σ_j (m_j / ρ_j) * (v_i - v_j) · ∇W_ij
/// ```
///
/// This uses the anti-symmetric velocity-gradient form.
#[allow(clippy::too_many_arguments)]
pub fn sph_drho_dt(
    pos_i: [f64; 3],
    vel_i: [f64; 3],
    rho_i: f64,
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_vel: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_rho: &[f64],
) -> f64 {
    let kernel = crate::kernel::CubicSplineKernel;
    let mut sum = 0.0_f64;
    for (((pos_j, vel_j), &m_j), &rho_j) in neighbor_pos
        .iter()
        .zip(neighbor_vel.iter())
        .zip(neighbor_mass.iter())
        .zip(neighbor_rho.iter())
    {
        let rij = [
            pos_i[0] - pos_j[0],
            pos_i[1] - pos_j[1],
            pos_i[2] - pos_j[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        if r < 1e-14 {
            continue;
        }
        let gw = kernel.grad_w(r, h);
        // ∇W_ij = (dW/dr) * r_hat
        let vij = [
            vel_i[0] - vel_j[0],
            vel_i[1] - vel_j[1],
            vel_i[2] - vel_j[2],
        ];
        // v_ij · ∇W_ij = v_ij · r_hat * (dW/dr)
        let v_dot_rhat = (vij[0] * rij[0] + vij[1] * rij[1] + vij[2] * rij[2]) / r;
        let rho_j_safe = rho_j.max(1e-14);
        sum += (m_j / rho_j_safe) * v_dot_rhat * gw;
    }
    -rho_i * sum
}

// ---------------------------------------------------------------------------
// Gather / scatter symmetric formulation
// ---------------------------------------------------------------------------

/// Gather-scatter symmetric kernel gradient for SPH (Monaghan 2005).
///
/// The symmetric form averages contributions from the two particles:
/// ```text
/// ∇W̃_ij = (∇W(r, h_i) + ∇W(r, h_j)) / 2
/// ```
///
/// This is used in the momentum equation to ensure linear-momentum conservation
/// even when smoothing lengths differ between particles.
///
/// Returns the 3-vector: `∇W̃_ij = (dW/dr @ h_i + dW/dr @ h_j) / 2 * r_hat`
pub fn gather_scatter_gradient(r_ij: [f64; 3], h_i: f64, h_j: f64) -> [f64; 3] {
    use crate::kernel::{CubicSplineKernel, SphKernel};
    let kernel = CubicSplineKernel;
    let r = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]).sqrt();
    if r < 1e-14 {
        return [0.0; 3];
    }
    let gw_i = kernel.grad_w(r, h_i);
    let gw_j = kernel.grad_w(r, h_j);
    let half_sum = 0.5 * (gw_i + gw_j);
    let inv_r = 1.0 / r;
    [
        r_ij[0] * half_sum * inv_r,
        r_ij[1] * half_sum * inv_r,
        r_ij[2] * half_sum * inv_r,
    ]
}

/// Compute the SPH density using the *scatter* form:
///
/// ```text
/// ρ_i = Σ_j m_j W(r_ij, h_j)
/// ```
///
/// In the scatter formulation each *source* particle j spreads its mass
/// with its own kernel width h_j (rather than using h_i of the receiver).
pub fn sph_density_scatter(
    pos_i: [f64; 3],
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_h: &[f64],
) -> f64 {
    let kernel = CubicSplineKernel;
    let mut rho = 0.0;
    for ((pos_j, &m_j), &h_j) in neighbor_pos.iter().zip(neighbor_mass).zip(neighbor_h) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += m_j * kernel.w(r, h_j);
    }
    rho
}

/// Compute the SPH density using the *gather* form:
///
/// ```text
/// ρ_i = Σ_j m_j W(r_ij, h_i)
/// ```
///
/// In the gather formulation the *receiver* particle i uses its own kernel
/// width h_i for all contributions.
pub fn sph_density_gather(
    pos_i: [f64; 3],
    h_i: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let kernel = CubicSplineKernel;
    let mut rho = 0.0;
    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += m_j * kernel.w(r, h_i);
    }
    rho
}

/// Symmetric SPH density: average of gather and scatter.
///
/// ```text
/// ρ_i^sym = (ρ_i^gather + ρ_i^scatter) / 2
/// ```
pub fn sph_density_symmetric(
    pos_i: [f64; 3],
    h_i: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_h: &[f64],
) -> f64 {
    let gather = sph_density_gather(pos_i, h_i, neighbor_pos, neighbor_mass);
    let scatter = sph_density_scatter(pos_i, neighbor_pos, neighbor_mass, neighbor_h);
    0.5 * (gather + scatter)
}

// ---------------------------------------------------------------------------
// Multi-resolution particle splitting (tetrahedral pattern)
// ---------------------------------------------------------------------------

/// Split a parent particle into 4 children arranged in a regular tetrahedron.
///
/// This 4-particle split is useful when only moderate refinement is needed
/// (vs the full 8-particle hexahedral split).  The tetrahedron vertices are
/// displaced from the parent centre by `offset = 0.5 * child_h`.
///
/// Conservation:
/// - mass: each child gets parent_mass / 4
/// - h   : child_h = parent_h / 4^(1/3) ≈ parent_h / 1.587
pub fn split_particle_4(parent: &AdaptiveHParticle) -> Vec<AdaptiveHParticle> {
    let n = 4_usize;
    let child_mass = parent.mass / n as f64;
    let child_h = parent.h / (n as f64).powf(1.0 / 3.0);
    let offset = 0.5 * child_h;

    // Regular tetrahedron vertices (unit edge length, scaled to `offset`)
    let sqrt2 = std::f64::consts::SQRT_2;
    let s = offset / sqrt2;
    let offsets: [[f64; 3]; 4] = [[s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]];

    offsets
        .iter()
        .map(|off| {
            let pos = [
                parent.position[0] + off[0],
                parent.position[1] + off[1],
                parent.position[2] + off[2],
            ];
            AdaptiveHParticle::new(
                pos,
                parent.velocity,
                child_mass,
                child_h,
                parent.density,
                parent.level + 1,
            )
        })
        .collect()
}

/// Merge two particles (2-to-1 merge) conserving mass, momentum and h.
///
/// Uses the same conservation laws as [`merge_particles`].
pub fn merge_pair(a: &AdaptiveHParticle, b: &AdaptiveHParticle) -> AdaptiveHParticle {
    merge_particles(&[a.clone(), b.clone()])
}

// ---------------------------------------------------------------------------
// h-ratio diagnostic
// ---------------------------------------------------------------------------

/// Compute the maximum h-ratio across all particles: `max_i h_i / h_j`
/// over all pairs (i, j) where j is a neighbour of i.
///
/// A large ratio (>> 1) indicates a rapid spatial resolution change and may
/// cause accuracy problems.
pub fn max_h_ratio(positions: &[[f64; 3]], hs: &[f64], radius: f64) -> f64 {
    let n = positions.len();
    let mut max_ratio = 1.0_f64;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < radius * radius {
                let ratio = hs[i] / hs[j].max(1e-30);
                if ratio > max_ratio {
                    max_ratio = ratio;
                }
            }
        }
    }
    max_ratio
}

// ---------------------------------------------------------------------------
// Variable-h pressure force (Springel & Hernquist symmetric)
// ---------------------------------------------------------------------------

/// Compute the symmetric variable-h pressure acceleration on particle i.
///
/// Uses the Springel & Hernquist (2002) conservative form:
/// ```text
/// a_i = -Σ_j m_j [ f_i P_i / (Ω_i ρ_i²) ∇W(r_ij, h_i)
///                 + f_j P_j / (Ω_j ρ_j²) ∇W(r_ij, h_j) ]
/// ```
/// where `f_i = (1 + (h_i / (d * ρ_i)) * dρ_i/dh)^{-1}` ≈ 1/Ω_i.
///
/// Returns a 3-vector acceleration.
#[allow(clippy::too_many_arguments)]
pub fn variable_h_pressure_accel(
    pos_i: [f64; 3],
    pressure_i: f64,
    rho_i: f64,
    h_i: f64,
    omega_i: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_pressure: &[f64],
    neighbor_rho: &[f64],
    neighbor_h: &[f64],
    neighbor_omega: &[f64],
) -> [f64; 3] {
    let kernel = CubicSplineKernel;
    let mut accel = [0.0_f64; 3];
    let rho_i_sq = (rho_i * rho_i).max(1e-28);
    let omega_i_safe = omega_i.max(1e-14);

    for i in 0..neighbor_pos.len() {
        let pos_j = neighbor_pos[i];
        let m_j = neighbor_mass[i];
        let p_j = neighbor_pressure[i];
        let rho_j = neighbor_rho[i].max(1e-14);
        let h_j = neighbor_h[i];
        let omega_j = neighbor_omega[i].max(1e-14);

        let rij = [
            pos_i[0] - pos_j[0],
            pos_i[1] - pos_j[1],
            pos_i[2] - pos_j[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        if r < 1e-14 {
            continue;
        }
        let gw_i = kernel.grad_w(r, h_i);
        let gw_j = kernel.grad_w(r, h_j);
        let factor_i = pressure_i / (omega_i_safe * rho_i_sq);
        let factor_j = p_j / (omega_j * rho_j * rho_j);
        let scale = -m_j * (factor_i * gw_i + factor_j * gw_j) / r;
        accel[0] += scale * rij[0];
        accel[1] += scale * rij[1];
        accel[2] += scale * rij[2];
    }
    accel
}

// ---------------------------------------------------------------------------
// KernelMassConservation check
// ---------------------------------------------------------------------------

/// Verify that the kernel sum Σ_j m_j W(r_ij, h) ≈ ρ₀ for an ideal uniform
/// distribution.  Returns the relative error |Σ - ρ₀| / ρ₀.
pub fn kernel_mass_conservation_error(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    rho0: f64,
) -> f64 {
    let rho = sph_density(pos_i, h, neighbor_pos, neighbor_mass);
    (rho - rho0).abs() / rho0.max(1e-30)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Hoover criterion tests
    // -----------------------------------------------------------------------

    #[test]
    fn hoover_h_at_reference_density() {
        // At ρ_i = ρ₀, hoover_h should return h₀.
        let h = hoover_h(0.1, 1000.0, 1000.0, 3);
        assert!(
            (h - 0.1).abs() < 1e-14,
            "At reference density, h should equal h₀: {h}"
        );
    }

    #[test]
    fn hoover_h_denser_fluid_gives_smaller_h() {
        // Higher density → smaller h (need fewer radius to find same number of neighbours).
        let h_dense = hoover_h(0.1, 1000.0, 2000.0, 3);
        let h_ref = 0.1_f64;
        assert!(
            h_dense < h_ref,
            "Denser fluid should have smaller h: {h_dense} vs {h_ref}"
        );
    }

    #[test]
    fn hoover_h_rarefied_fluid_gives_larger_h() {
        // Lower density → larger h.
        let h_rare = hoover_h(0.1, 1000.0, 500.0, 3);
        assert!(
            h_rare > 0.1,
            "Rarefied fluid should have larger h: {h_rare}"
        );
    }

    #[test]
    fn hoover_residual_zero_at_solution() {
        let h0 = 0.1;
        let rho0 = 1000.0;
        let rho_i = 800.0;
        let h_sol = hoover_h(h0, rho0, rho_i, 3);
        let res = hoover_residual(h0, h_sol, rho0, rho_i, 3);
        assert!(
            res.abs() < 1e-10,
            "Residual at Hoover solution should be ~0: {res}"
        );
    }

    // -----------------------------------------------------------------------
    // SPH density sum tests
    // -----------------------------------------------------------------------

    #[test]
    fn sph_density_single_neighbour_at_origin() {
        let pos_i = [0.0_f64; 3];
        let pos_j = [[0.0_f64; 3]];
        let mass_j = [1.0_f64];
        let h = 0.1;
        let kernel = CubicSplineKernel;
        let rho = sph_density(pos_i, h, &pos_j, &mass_j);
        let expected = kernel.w(0.0, h);
        assert!(
            (rho - expected).abs() < 1e-14,
            "Self contribution: {rho} vs {expected}"
        );
    }

    #[test]
    fn sph_density_zero_outside_support() {
        let pos_i = [0.0_f64; 3];
        let pos_j = [[5.0_f64, 0.0, 0.0]]; // r = 5 >> 2h
        let mass_j = [1.0_f64];
        let h = 0.1;
        let rho = sph_density(pos_i, h, &pos_j, &mass_j);
        assert!(
            rho.abs() < 1e-14,
            "Particle outside support should contribute 0: {rho}"
        );
    }

    #[test]
    fn sph_density_increases_with_more_neighbours() {
        let pos_i = [0.0_f64; 3];
        let h = 0.1;
        let pos_1 = vec![[0.05_f64, 0.0, 0.0]];
        let pos_2 = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let mass = vec![1.0_f64];
        let mass2 = vec![1.0_f64, 1.0_f64];

        let rho1 = sph_density(pos_i, h, &pos_1, &mass);
        let rho2 = sph_density(pos_i, h, &pos_2, &mass2);
        assert!(rho2 > rho1, "More neighbours → higher density");
    }

    // -----------------------------------------------------------------------
    // Newton-Raphson solver tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve_h_newton_converges_to_known_density() {
        let pos_i = [0.0_f64; 3];
        let h0 = 0.1;
        let positions = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0], [0.0, 0.05, 0.0]];
        let masses = vec![0.001_f64; 3];

        // Compute the density at h0 to use as target
        let rho_target = sph_density(pos_i, h0, &positions, &masses);

        let h_solved = solve_h_newton(
            pos_i, h0, rho_target, &positions, &masses, 100, 1e-8, 0.01, 1.0,
        );

        assert!(h_solved.is_some(), "Newton solver should converge");
        let h_s = h_solved.unwrap();
        let rho_check = sph_density(pos_i, h_s, &positions, &masses);
        let rel_err = (rho_check - rho_target).abs() / rho_target.max(1e-14);
        assert!(
            rel_err < 1e-6,
            "Solved h gives wrong density: target={rho_target}, got={rho_check}"
        );
    }

    #[test]
    fn solve_h_newton_clamps_to_bounds() {
        let pos_i = [0.0_f64; 3];
        let h_min = 0.05;
        let h_max = 0.5;
        // Use a very low target density (only possible with very large h)
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let masses = vec![0.0001_f64];
        let h_solved = solve_h_newton(
            pos_i, 0.1, 1e10, // unreachably large target
            &positions, &masses, 50, 1e-6, h_min, h_max,
        );
        // Should either converge within bounds or return None (but if Some, must be in bounds)
        if let Some(h) = h_solved {
            assert!(
                h >= h_min && h <= h_max,
                "Solved h must be in [h_min, h_max]: {h}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Omega correction tests
    // -----------------------------------------------------------------------

    #[test]
    fn omega_correction_near_one_for_few_neighbours() {
        // With no neighbours, drho/dh = 0 → Ω = 1
        let pos_i = [0.0_f64; 3];
        let omega = omega_correction(pos_i, 0.1, &[], &[], 1000.0, 3.0);
        assert!((omega - 1.0).abs() < 1e-14, "No neighbours → Ω = 1");
    }

    #[test]
    fn omega_correction_finite_for_uniform_cluster() {
        let pos_i = [0.0_f64; 3];
        let positions: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64 * 0.02 - 0.09, 0.0, 0.0])
            .collect();
        let masses = vec![0.001_f64; 10];
        let rho = sph_density(pos_i, 0.1, &positions, &masses);
        let omega = omega_correction(pos_i, 0.1, &positions, &masses, rho, 3.0);
        assert!(omega.is_finite(), "Ω must be finite: {omega}");
    }

    // -----------------------------------------------------------------------
    // Variable-h kernel gradient tests
    // -----------------------------------------------------------------------

    #[test]
    fn variable_h_gradient_zero_at_coincident_positions() {
        let r_ij = [0.0_f64; 3];
        let grad = variable_h_kernel_gradient(r_ij, 0.1, 0.1, 1.0, 1.0);
        let mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
        assert!(mag < 1e-14, "Zero separation → zero gradient");
    }

    #[test]
    fn variable_h_gradient_symmetric_equal_h() {
        // With h_i = h_j and Ω_i = Ω_j, the symmetric gradient should be the same
        // as the standard kernel gradient.
        let r_ij = [0.05_f64, 0.0, 0.0];
        let h = 0.1;
        let kernel = CubicSplineKernel;
        let r = 0.05;
        let gw = kernel.grad_w(r, h);

        let grad = variable_h_kernel_gradient(r_ij, h, h, 1.0, 1.0);
        let expected = gw / r * r_ij[0];
        assert!(
            (grad[0] - expected).abs() < 1e-12,
            "Equal-h gradient: {}, expected {}",
            grad[0],
            expected
        );
    }

    // -----------------------------------------------------------------------
    // H bounds tests
    // -----------------------------------------------------------------------

    #[test]
    fn apply_h_bounds_clamps_values() {
        let mut hs = vec![0.001, 0.1, 5.0, 0.5];
        apply_h_bounds(&mut hs, 0.01, 1.0);
        assert!((hs[0] - 0.01).abs() < 1e-14, "Too small clamped to h_min");
        assert!((hs[1] - 0.1).abs() < 1e-14, "In range unchanged");
        assert!((hs[2] - 1.0).abs() < 1e-14, "Too large clamped to h_max");
        assert!((hs[3] - 0.5).abs() < 1e-14, "In range unchanged");
    }

    #[test]
    fn h_bounds_from_reference_uses_fractions() {
        let h0 = 0.1;
        let (h_min, h_max) = h_bounds_from_reference(h0);
        assert!((h_min - h0 * H_MIN_FRACTION).abs() < 1e-14);
        assert!((h_max - h0 * H_MAX_FRACTION).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // Neighbour count estimate tests
    // -----------------------------------------------------------------------

    #[test]
    fn neighbour_count_single_at_origin() {
        let pos_i = [0.0_f64; 3];
        let positions = vec![[0.0_f64; 3]]; // the particle itself
        let count = estimate_neighbour_count(pos_i, 0.1, &positions);
        assert_eq!(count, 1, "Should count itself");
    }

    #[test]
    fn neighbour_count_excludes_outside_radius() {
        let pos_i = [0.0_f64; 3];
        let positions = vec![[1.0_f64, 0.0, 0.0]]; // far away
        let count = estimate_neighbour_count(pos_i, 0.1, &positions);
        assert_eq!(count, 0, "Far particle should not be counted");
    }

    #[test]
    fn neighbour_count_uniform_grid() {
        let pos_i = [0.0_f64; 3];
        // 3×3×3 grid of particles at spacing 0.05 → all within 2h = 0.2 (h=0.1)
        let positions: Vec<[f64; 3]> = (0..3)
            .flat_map(|ix| {
                (0..3).flat_map(move |iy| {
                    (0..3).map(move |iz| {
                        [
                            (ix as f64 - 1.0) * 0.05,
                            (iy as f64 - 1.0) * 0.05,
                            (iz as f64 - 1.0) * 0.05,
                        ]
                    })
                })
            })
            .collect();
        let count = estimate_neighbour_count(pos_i, 0.15, &positions);
        assert!(count > 0, "Should find at least some neighbours: {count}");
    }

    // -----------------------------------------------------------------------
    // Particle split/merge tests
    // -----------------------------------------------------------------------

    #[test]
    fn split_particle_gives_eight_children() {
        let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 0);
        let children = split_particle(&parent);
        assert_eq!(children.len(), 8, "Should produce 8 children");
    }

    #[test]
    fn split_conserves_mass() {
        let parent_mass = 2.5;
        let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], parent_mass, 0.1, 1000.0, 0);
        let children = split_particle(&parent);
        let total_child_mass: f64 = children.iter().map(|c| c.mass).sum();
        assert!(
            (total_child_mass - parent_mass).abs() < 1e-14,
            "Split must conserve mass: {total_child_mass} vs {parent_mass}"
        );
    }

    #[test]
    fn split_centre_of_mass_conserved() {
        let parent_pos = [1.0_f64, 2.0, 3.0];
        let parent = AdaptiveHParticle::new(parent_pos, [0.0; 3], 1.0, 0.1, 1000.0, 0);
        let children = split_particle(&parent);
        let n = children.len() as f64;
        let cm = [
            children.iter().map(|c| c.position[0]).sum::<f64>() / n,
            children.iter().map(|c| c.position[1]).sum::<f64>() / n,
            children.iter().map(|c| c.position[2]).sum::<f64>() / n,
        ];
        for k in 0..3 {
            assert!(
                (cm[k] - parent_pos[k]).abs() < 1e-12,
                "Centre of mass dim {k}: {}, expected {}",
                cm[k],
                parent_pos[k]
            );
        }
    }

    #[test]
    fn merge_particles_conserves_mass() {
        let particles: Vec<AdaptiveHParticle> = (0..8)
            .map(|i| AdaptiveHParticle::new([i as f64 * 0.01; 3], [0.0; 3], 0.125, 0.05, 1000.0, 1))
            .collect();
        let total_before: f64 = particles.iter().map(|p| p.mass).sum();
        let merged = merge_particles(&particles);
        assert!(
            (merged.mass - total_before).abs() < 1e-14,
            "Merge must conserve mass: {} vs {}",
            merged.mass,
            total_before
        );
    }

    #[test]
    fn merge_particles_centre_of_mass() {
        // Two equal-mass particles at x=0 and x=2 → merged at x=1
        let p1 = AdaptiveHParticle::new([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.05, 1000.0, 1);
        let p2 = AdaptiveHParticle::new([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.05, 1000.0, 1);
        let merged = merge_particles(&[p1, p2]);
        assert!(
            (merged.position[0] - 1.0).abs() < 1e-14,
            "CoM should be at x=1"
        );
    }

    #[test]
    fn split_then_merge_preserves_mass() {
        let parent = AdaptiveHParticle::new([0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.1, 1000.0, 0);
        let original_mass = parent.mass;
        let children = split_particle(&parent);
        let merged = merge_particles(&children);
        assert!(
            (merged.mass - original_mass).abs() < 1e-14,
            "Round-trip mass: {} vs {}",
            merged.mass,
            original_mass
        );
    }

    // -----------------------------------------------------------------------
    // AdaptiveHManager tests
    // -----------------------------------------------------------------------

    #[test]
    fn manager_update_smoothing_lengths_hoover() {
        let h0 = 0.1;
        let rho0 = 1000.0;
        let mut mgr = AdaptiveHManager::new(rho0, h0, 32.0, 3);

        // Add particle at double the reference density → h should halve (3D)
        mgr.add_particle(AdaptiveHParticle::new(
            [0.0; 3],
            [0.0; 3],
            1.0,
            h0,
            2.0 * rho0,
            0,
        ));
        mgr.update_smoothing_lengths();

        let expected_h = hoover_h(h0, rho0, 2.0 * rho0, 3);
        let actual_h = mgr.particles[0].h;
        assert!(
            (actual_h - expected_h).abs() < 1e-12,
            "Hoover h mismatch: {actual_h} vs {expected_h}"
        );
    }

    #[test]
    fn manager_total_mass_conserved_after_refine() {
        let mut mgr = AdaptiveHManager::new(1000.0, 0.1, 32.0, 3);
        for i in 0..4 {
            mgr.add_particle(AdaptiveHParticle::new(
                [i as f64 * 0.1, 0.0, 0.0],
                [0.0; 3],
                1.0,
                0.5, // Large h → will be split
                1000.0,
                0,
            ));
        }
        let mass_before = mgr.total_mass();
        mgr.refine_coarse_particles(0.3); // threshold < 0.5 → all 4 split
        let mass_after = mgr.total_mass();
        assert!(
            (mass_after - mass_before).abs() < 1e-12,
            "Mass after refine: {mass_after} vs {mass_before}"
        );
    }

    #[test]
    fn manager_omega_all_finite() {
        let mut mgr = AdaptiveHManager::new(1000.0, 0.1, 32.0, 3);
        for i in 0..5 {
            mgr.add_particle(AdaptiveHParticle::new(
                [i as f64 * 0.05, 0.0, 0.0],
                [0.0; 3],
                0.001,
                0.1,
                1000.0,
                0,
            ));
        }
        // Update densities from SPH sum
        let positions: Vec<[f64; 3]> = mgr.particles.iter().map(|p| p.position).collect();
        let masses: Vec<f64> = mgr.particles.iter().map(|p| p.mass).collect();
        for p in &mut mgr.particles {
            p.density = sph_density(p.position, p.h, &positions, &masses).max(1e-14);
        }

        let omegas = mgr.compute_omega();
        for (i, &w) in omegas.iter().enumerate() {
            assert!(w.is_finite(), "Omega[{i}] should be finite: {w}");
        }
    }

    // -----------------------------------------------------------------------
    // Zhu-Fox criterion tests
    // -----------------------------------------------------------------------

    #[test]
    fn zhu_fox_uniform_density_zero() {
        let pos_i = [0.0_f64; 3];
        let rho = 1000.0;
        let h = 0.1;
        // All neighbours at same density → epsilon = 0
        let neighbour_pos = vec![[0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let neighbour_rho = vec![rho, rho];
        let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
        assert!(eps.abs() < 1e-14, "Uniform density → eps=0, got {eps}");
    }

    #[test]
    fn zhu_fox_varying_density_positive() {
        let pos_i = [0.0_f64; 3];
        let rho = 1000.0;
        let h = 0.1;
        let neighbour_pos = vec![[0.05, 0.0, 0.0]];
        let neighbour_rho = vec![1500.0_f64]; // 50% higher
        let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
        assert!(
            (eps - 0.5).abs() < 1e-10,
            "Should report 50% variation, got {eps}"
        );
    }

    #[test]
    fn zhu_fox_far_neighbour_ignored() {
        let pos_i = [0.0_f64; 3];
        let rho = 1000.0;
        let h = 0.05; // small h
        // Neighbour is at r = 0.1 > h, so outside the criterion radius
        let neighbour_pos = vec![[0.1, 0.0, 0.0]];
        let neighbour_rho = vec![5000.0_f64];
        let eps = zhu_fox_criterion(pos_i, rho, h, &neighbour_pos, &neighbour_rho);
        assert!(eps.abs() < 1e-14, "Far neighbour should be ignored: {eps}");
    }

    #[test]
    fn should_refine_zhu_fox_at_max_level_returns_false() {
        let pos_i = [0.0_f64; 3];
        let rho = 1000.0;
        let result = should_refine_zhu_fox(pos_i, rho, 0.1, 5, 5, &[], &[], 0.1);
        assert!(!result, "At max_level should never refine");
    }

    #[test]
    fn should_refine_zhu_fox_large_variation_returns_true() {
        let pos_i = [0.0_f64; 3];
        let rho = 1000.0;
        let h = 0.1;
        let neighbour_pos = vec![[0.05, 0.0, 0.0]];
        let neighbour_rho = vec![2000.0_f64]; // 100% variation
        let result =
            should_refine_zhu_fox(pos_i, rho, h, 0, 5, &neighbour_pos, &neighbour_rho, 0.5);
        assert!(result, "100% variation > threshold 50% → should refine");
    }

    // -----------------------------------------------------------------------
    // sph_drho_dt tests
    // -----------------------------------------------------------------------

    #[test]
    fn drho_dt_zero_for_zero_velocity() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [0.0; 3];
        let rho_i = 1000.0;
        let h = 0.1;
        let neighbour_pos = vec![[0.05, 0.0, 0.0]];
        let neighbour_vel = vec![[0.0_f64; 3]]; // same velocity → zero relative
        let neighbour_mass = vec![0.001_f64];
        let neighbour_rho = vec![1000.0_f64];
        let drho = sph_drho_dt(
            pos_i,
            vel_i,
            rho_i,
            h,
            &neighbour_pos,
            &neighbour_vel,
            &neighbour_mass,
            &neighbour_rho,
        );
        assert!(
            drho.abs() < 1e-14,
            "Zero relative velocity → dρ/dt=0, got {drho}"
        );
    }

    #[test]
    fn drho_dt_finite_for_nonzero_velocity() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let rho_i = 1000.0;
        let h = 0.1;
        let neighbour_pos = vec![[0.05, 0.0, 0.0]];
        let neighbour_vel = vec![[0.0_f64; 3]];
        let neighbour_mass = vec![0.001_f64];
        let neighbour_rho = vec![1000.0_f64];
        let drho = sph_drho_dt(
            pos_i,
            vel_i,
            rho_i,
            h,
            &neighbour_pos,
            &neighbour_vel,
            &neighbour_mass,
            &neighbour_rho,
        );
        assert!(drho.is_finite(), "dρ/dt should be finite");
    }

    // -----------------------------------------------------------------------
    // update_h_continuity tests
    // -----------------------------------------------------------------------

    #[test]
    fn update_h_continuity_no_change_for_zero_drho() {
        let h = 0.1;
        let h_new = update_h_continuity(h, 1000.0, 0.0, 0.001, 3.0, 0.01, 1.0);
        assert!((h_new - h).abs() < 1e-14, "Zero dρ/dt → h unchanged");
    }

    #[test]
    fn update_h_continuity_increases_for_expanding_flow() {
        // Expanding flow: dρ/dt < 0 → dh/dt > 0 → h increases
        let h = 0.1;
        let h_new = update_h_continuity(h, 1000.0, -100.0, 0.001, 3.0, 0.01, 1.0);
        assert!(
            h_new > h,
            "Expanding flow should increase h: {h_new} vs {h}"
        );
    }

    #[test]
    fn update_h_continuity_clamps_to_bounds() {
        let h = 0.1;
        // drho/dt very large → h would go very small
        let h_new = update_h_continuity(h, 1000.0, 1e10, 0.001, 3.0, 0.05, 1.0);
        assert!(h_new >= 0.05, "Should be clamped to h_min: {h_new}");
    }

    // -----------------------------------------------------------------------
    // Gather / scatter density tests
    // -----------------------------------------------------------------------

    #[test]
    fn gather_density_equals_standard() {
        let pos_i = [0.0_f64; 3];
        let h = 0.1;
        let positions = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let masses = vec![0.001_f64, 0.001];
        let rho_gather = sph_density_gather(pos_i, h, &positions, &masses);
        let rho_std = sph_density(pos_i, h, &positions, &masses);
        assert!(
            (rho_gather - rho_std).abs() < 1e-12,
            "Gather == standard for uniform h"
        );
    }

    #[test]
    fn scatter_density_positive() {
        let pos_i = [0.0_f64; 3];
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let masses = vec![0.001_f64];
        let hs = vec![0.1_f64];
        let rho = sph_density_scatter(pos_i, &positions, &masses, &hs);
        assert!(rho > 0.0, "Scatter density should be positive: {rho}");
    }

    #[test]
    fn symmetric_density_between_gather_and_scatter() {
        let pos_i = [0.0_f64; 3];
        let h_i = 0.1;
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let masses = vec![0.001_f64];
        let hs_j = vec![0.08_f64]; // different h for neighbour
        let rho_g = sph_density_gather(pos_i, h_i, &positions, &masses);
        let rho_s = sph_density_scatter(pos_i, &positions, &masses, &hs_j);
        let rho_sym = sph_density_symmetric(pos_i, h_i, &positions, &masses, &hs_j);
        let expected = 0.5 * (rho_g + rho_s);
        assert!(
            (rho_sym - expected).abs() < 1e-14,
            "Symmetric = avg(gather, scatter)"
        );
    }

    // -----------------------------------------------------------------------
    // Gather-scatter gradient test
    // -----------------------------------------------------------------------

    #[test]
    fn gather_scatter_gradient_zero_at_origin() {
        let r_ij = [0.0_f64; 3];
        let g = gather_scatter_gradient(r_ij, 0.1, 0.1);
        let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!(mag < 1e-14, "Zero separation → zero gradient");
    }

    #[test]
    fn gather_scatter_gradient_symmetric_equal_h() {
        // With h_i = h_j, gather-scatter == standard gradient
        let r_ij = [0.05_f64, 0.0, 0.0];
        let h = 0.1;
        let kernel = CubicSplineKernel;
        let r = 0.05;
        let gw = kernel.grad_w(r, h);
        let g = gather_scatter_gradient(r_ij, h, h);
        let expected_x = gw / r * r_ij[0];
        assert!(
            (g[0] - expected_x).abs() < 1e-12,
            "Equal-h gather-scatter == standard"
        );
    }

    // -----------------------------------------------------------------------
    // 4-particle split tests
    // -----------------------------------------------------------------------

    #[test]
    fn split_4_gives_four_children() {
        let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 0);
        let children = split_particle_4(&parent);
        assert_eq!(children.len(), 4);
    }

    #[test]
    fn split_4_conserves_mass() {
        let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 3.0, 0.1, 1000.0, 0);
        let children = split_particle_4(&parent);
        let total: f64 = children.iter().map(|c| c.mass).sum();
        assert!((total - 3.0).abs() < 1e-14, "Split-4 should conserve mass");
    }

    #[test]
    fn split_4_increments_level() {
        let parent = AdaptiveHParticle::new([0.0; 3], [0.0; 3], 1.0, 0.1, 1000.0, 2);
        let children = split_particle_4(&parent);
        for c in &children {
            assert_eq!(c.level, 3, "Children should be at level parent.level + 1");
        }
    }

    // -----------------------------------------------------------------------
    // merge_pair test
    // -----------------------------------------------------------------------

    #[test]
    fn merge_pair_conserves_mass() {
        let a = AdaptiveHParticle::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 0.1, 1000.0, 1);
        let b = AdaptiveHParticle::new([0.1; 3], [0.0, 1.0, 0.0], 0.5, 0.1, 1000.0, 1);
        let merged = merge_pair(&a, &b);
        assert!((merged.mass - 1.0).abs() < 1e-14);
    }

    #[test]
    fn merge_pair_momentum_conserved() {
        let a = AdaptiveHParticle::new([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.1, 1000.0, 1);
        let b = AdaptiveHParticle::new([0.0; 3], [0.0, 0.0, 0.0], 1.0, 0.1, 1000.0, 1);
        let merged = merge_pair(&a, &b);
        // Momentum: m_a * v_a + m_b * v_b = 1*2 + 1*0 = 2 kg m/s
        // merged mass = 2, merged velocity should be 1.0 m/s
        assert!(
            (merged.velocity[0] - 1.0).abs() < 1e-14,
            "Velocity after merge: {}",
            merged.velocity[0]
        );
    }

    // -----------------------------------------------------------------------
    // max_h_ratio test
    // -----------------------------------------------------------------------

    #[test]
    fn max_h_ratio_uniform_is_one() {
        let positions = vec![[0.0_f64; 3], [0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let hs = vec![0.1_f64; 3];
        let ratio = max_h_ratio(&positions, &hs, 0.2);
        assert!(
            (ratio - 1.0).abs() < 1e-10,
            "Uniform h → ratio 1.0, got {ratio}"
        );
    }

    #[test]
    fn max_h_ratio_detects_difference() {
        let positions = vec![[0.0_f64; 3], [0.05, 0.0, 0.0]];
        let hs = vec![0.2_f64, 0.1_f64]; // ratio 2
        let ratio = max_h_ratio(&positions, &hs, 0.2);
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Expected ratio 2.0, got {ratio}"
        );
    }

    // -----------------------------------------------------------------------
    // kernel_mass_conservation_error test
    // -----------------------------------------------------------------------

    #[test]
    fn kernel_mass_conservation_finite() {
        let pos_i = [0.0_f64; 3];
        let h = 0.1;
        let positions: Vec<[f64; 3]> = (0..8)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 8.0;
                [0.07 * angle.cos(), 0.07 * angle.sin(), 0.0]
            })
            .collect();
        let masses = vec![0.001_f64; 8];
        let err = kernel_mass_conservation_error(pos_i, h, &positions, &masses, 1000.0);
        assert!(
            err.is_finite(),
            "Mass conservation error should be finite: {err}"
        );
    }

    // -----------------------------------------------------------------------
    // variable_h_pressure_accel test
    // -----------------------------------------------------------------------

    #[test]
    fn variable_h_pressure_accel_zero_for_zero_pressure() {
        let pos_i = [0.0_f64; 3];
        let accel = variable_h_pressure_accel(
            pos_i,
            0.0,
            1000.0,
            0.1,
            1.0,
            &[[0.05, 0.0, 0.0_f64]],
            &[0.001],
            &[0.0], // neighbour pressure = 0
            &[1000.0],
            &[0.1],
            &[1.0],
        );
        let mag = (accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]).sqrt();
        assert!(mag < 1e-14, "Zero pressure → zero force: {mag}");
    }

    #[test]
    fn variable_h_pressure_accel_finite_for_nonzero_pressure() {
        let pos_i = [0.0_f64; 3];
        let accel = variable_h_pressure_accel(
            pos_i,
            1000.0,
            1000.0,
            0.1,
            1.0,
            &[[0.05, 0.0, 0.0_f64]],
            &[0.001],
            &[500.0],
            &[1000.0],
            &[0.1],
            &[1.0],
        );
        for a in accel {
            assert!(
                a.is_finite(),
                "Acceleration component should be finite: {a}"
            );
        }
    }
}
