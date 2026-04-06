#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle shifting for SPH regularity.
//!
//! Particle shifting techniques prevent disorder (tensile instability,
//! void formation, particle clustering) in Lagrangian SPH simulations
//! by periodically redistributing particles toward a more uniform arrangement.

/// Parameters controlling the particle shifting algorithm.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ParticleShift {
    /// Shifting coefficient (typically 0.01–0.5).
    pub coefficient: f64,
    /// Smoothing length reference for normalising shift magnitude.
    pub smoothing_length: f64,
    /// Maximum allowed shift per time step as a fraction of `smoothing_length`.
    pub max_shift_fraction: f64,
    /// Artificial tensile correction parameter `epsilon` (typically ~0.2).
    pub tensile_epsilon: f64,
    /// Exponent for tensile correction (typically 4).
    pub tensile_exponent: i32,
    /// Spatial dimension (2 or 3).
    pub dim: usize,
}

impl ParticleShift {
    /// Create a new `ParticleShift` with default parameters.
    ///
    /// # Arguments
    /// * `coefficient` - shifting coefficient A in the concentration gradient formula
    /// * `smoothing_length` - particle smoothing length h
    /// * `dim` - spatial dimension
    pub fn new(coefficient: f64, smoothing_length: f64, dim: usize) -> Self {
        Self {
            coefficient,
            smoothing_length,
            max_shift_fraction: 0.5,
            tensile_epsilon: 0.2,
            tensile_exponent: 4,
            dim,
        }
    }

    /// Maximum shift magnitude for a single step.
    pub fn max_shift(&self) -> f64 {
        self.max_shift_fraction * self.smoothing_length
    }
}

/// Compute the concentration-gradient shifting vector for a particle.
///
/// Uses the formulation:
/// `delta_r_i = -A * h^2 * sum_j (1 + R*(W_ij/W_0)^n) * grad_W_ij`
/// where the tensile correction is omitted in this base version (R=0).
///
/// # Arguments
/// * `grad_w` - kernel gradient contributions from each neighbour, shape `[n_neigh][dim]`
/// * `mass_ratio` - mass_j / mass_i for each neighbour (use `1.0` for equal masses)
/// * `coefficient` - shifting coefficient A
/// * `h` - smoothing length
/// * `dim` - spatial dimension (must match `grad_w[i].len()`)
///
/// Returns the shifting vector of length `dim`.
pub fn shifting_vector(
    grad_w: &[Vec<f64>],
    mass_ratio: &[f64],
    coefficient: f64,
    h: f64,
    dim: usize,
) -> Vec<f64> {
    assert_eq!(grad_w.len(), mass_ratio.len());
    let mut delta = vec![0.0f64; dim];
    let scale = coefficient * h * h;
    for (j, gw) in grad_w.iter().enumerate() {
        for d in 0..dim {
            delta[d] -= scale * mass_ratio[j] * gw[d];
        }
    }
    delta
}

/// Fickian diffusion-based particle shifting.
///
/// Computes a diffusion-type shift based on the inter-particle concentration
/// gradient, designed to equidistribute particles:
/// `delta_r_i = D * dt * sum_j (r_j - r_i) * W_ij * V_j`
///
/// where `D` is the diffusion coefficient and `V_j` is the volume of neighbour `j`.
///
/// # Arguments
/// * `pos_i` - position of particle i, length `dim`
/// * `pos_j` - positions of neighbours, shape `[n_neigh][dim]`
/// * `w_ij` - kernel values at each neighbour
/// * `vol_j` - volumes of neighbours
/// * `diffusion_coeff` - Fickian diffusion coefficient D
/// * `dt` - time step
/// * `dim` - spatial dimension
#[allow(clippy::too_many_arguments)]
pub fn fickian_shift(
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    vol_j: &[f64],
    diffusion_coeff: f64,
    dt: f64,
    dim: usize,
) -> Vec<f64> {
    assert_eq!(pos_j.len(), w_ij.len());
    assert_eq!(pos_j.len(), vol_j.len());
    let mut delta = vec![0.0f64; dim];
    for j in 0..pos_j.len() {
        let wv = w_ij[j] * vol_j[j];
        for d in 0..dim {
            delta[d] += (pos_j[j][d] - pos_i[d]) * wv;
        }
    }
    for d in 0..dim {
        delta[d] *= diffusion_coeff * dt;
    }
    delta
}

/// Artificial tensile instability correction term.
///
/// Adds a positive small correction to the shifting vector to prevent
/// particle clustering under tension. Based on:
/// `R_ij = epsilon * (W_ij / W_0)^n`
/// where `W_0 = W(0)` is the kernel value at zero distance.
///
/// # Arguments
/// * `w_ij` - kernel value between particle i and neighbour j
/// * `w_0` - kernel value at zero distance
/// * `epsilon` - tensile correction coefficient (typically 0.2)
/// * `n_exp` - tensile exponent (typically 4)
pub fn tensile_correction(w_ij: f64, w_0: f64, epsilon: f64, n_exp: i32) -> f64 {
    if w_0.abs() < 1e-20 {
        return 0.0;
    }
    let ratio = w_ij / w_0;
    epsilon * ratio.powi(n_exp)
}

/// Apply computed shift vectors to update particle positions.
///
/// Clamps each shift to `max_shift` to maintain stability.
///
/// # Arguments
/// * `positions` - particle positions, shape `[n_particles][dim]` (modified in place)
/// * `shifts` - shift vectors, shape `[n_particles][dim]`
/// * `max_shift` - maximum allowed magnitude of shift per particle
/// * `dim` - spatial dimension
pub fn apply_shift(positions: &mut [Vec<f64>], shifts: &[Vec<f64>], max_shift: f64, dim: usize) {
    assert_eq!(positions.len(), shifts.len());
    for (pos, shift) in positions.iter_mut().zip(shifts.iter()) {
        // Compute magnitude of shift
        let mag_sq: f64 = shift.iter().map(|&s| s * s).sum();
        let mag = mag_sq.sqrt();
        let scale = if mag > max_shift && mag > 1e-20 {
            max_shift / mag
        } else {
            1.0
        };
        for d in 0..dim {
            pos[d] += scale * shift[d];
        }
    }
}

/// Check and correct particle positions against a box boundary.
///
/// If a particle has been shifted outside `[lo, hi]` in any dimension,
/// its position is reflected back to the boundary. Returns the number
/// of particles whose shift was clamped.
///
/// # Arguments
/// * `positions` - particle positions, shape `[n_particles][dim]` (modified in place)
/// * `lo` - lower bound per dimension
/// * `hi` - upper bound per dimension
/// * `dim` - spatial dimension
pub fn shift_boundary_check(
    positions: &mut [Vec<f64>],
    lo: &[f64],
    hi: &[f64],
    dim: usize,
) -> usize {
    assert_eq!(lo.len(), dim);
    assert_eq!(hi.len(), dim);
    let mut n_clamped = 0;
    for pos in positions.iter_mut() {
        let mut clamped = false;
        for d in 0..dim {
            if pos[d] < lo[d] {
                pos[d] = lo[d];
                clamped = true;
            } else if pos[d] > hi[d] {
                pos[d] = hi[d];
                clamped = true;
            }
        }
        if clamped {
            n_clamped += 1;
        }
    }
    n_clamped
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── ParticleShift ──────────────────────────────────────────────────────

    #[test]
    fn test_particle_shift_new() {
        let ps = ParticleShift::new(0.2, 0.1, 2);
        assert!((ps.coefficient - 0.2).abs() < 1e-14);
        assert!((ps.smoothing_length - 0.1).abs() < 1e-14);
        assert_eq!(ps.dim, 2);
    }

    #[test]
    fn test_particle_shift_max_shift() {
        let ps = ParticleShift::new(0.2, 0.1, 2);
        assert!((ps.max_shift() - 0.05).abs() < 1e-14);
    }

    #[test]
    fn test_particle_shift_defaults() {
        let ps = ParticleShift::new(0.1, 0.05, 3);
        assert!((ps.max_shift_fraction - 0.5).abs() < 1e-14);
        assert_eq!(ps.tensile_exponent, 4);
    }

    #[test]
    fn test_particle_shift_clone() {
        let ps = ParticleShift::new(0.3, 0.2, 2);
        let ps2 = ps.clone();
        assert!((ps2.coefficient - 0.3).abs() < 1e-14);
    }

    // ── shifting_vector ────────────────────────────────────────────────────

    #[test]
    fn test_shifting_vector_no_neighbours() {
        let delta = shifting_vector(&[], &[], 0.2, 0.1, 2);
        assert_eq!(delta.len(), 2);
        assert!(delta[0].abs() < 1e-15);
        assert!(delta[1].abs() < 1e-15);
    }

    #[test]
    fn test_shifting_vector_single_neighbour() {
        let grad_w = vec![vec![1.0, 0.0]];
        let mass_ratio = vec![1.0];
        let delta = shifting_vector(&grad_w, &mass_ratio, 1.0, 1.0, 2);
        // delta[0] = -1.0 * 1.0 * 1.0^2 * 1.0 = -1
        assert!((delta[0] - (-1.0)).abs() < 1e-14);
        assert!(delta[1].abs() < 1e-14);
    }

    #[test]
    fn test_shifting_vector_symmetry_cancels() {
        // Two neighbours with opposite gradients → shift should cancel
        let grad_w = vec![vec![1.0, 0.0], vec![-1.0, 0.0]];
        let mass_ratio = vec![1.0, 1.0];
        let delta = shifting_vector(&grad_w, &mass_ratio, 0.2, 0.1, 2);
        assert!(delta[0].abs() < 1e-13);
    }

    #[test]
    fn test_shifting_vector_scales_with_h_squared() {
        let grad_w = vec![vec![1.0, 0.0]];
        let mass_ratio = vec![1.0];
        let d1 = shifting_vector(&grad_w, &mass_ratio, 1.0, 1.0, 2);
        let d2 = shifting_vector(&grad_w, &mass_ratio, 1.0, 2.0, 2);
        assert!((d2[0] / d1[0] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_shifting_vector_3d() {
        let grad_w = vec![vec![1.0, 2.0, 3.0]];
        let mass_ratio = vec![1.0];
        let delta = shifting_vector(&grad_w, &mass_ratio, 1.0, 1.0, 3);
        assert_eq!(delta.len(), 3);
        assert!((delta[0] - (-1.0)).abs() < 1e-14);
        assert!((delta[1] - (-2.0)).abs() < 1e-14);
        assert!((delta[2] - (-3.0)).abs() < 1e-14);
    }

    // ── fickian_shift ──────────────────────────────────────────────────────

    #[test]
    fn test_fickian_no_neighbours() {
        let pos_i = [0.0, 0.0];
        let delta = fickian_shift(&pos_i, &[], &[], &[], 0.1, 0.01, 2);
        assert!(delta[0].abs() < 1e-15);
    }

    #[test]
    fn test_fickian_single_neighbour_positive() {
        let pos_i = [0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0]];
        let w_ij = vec![1.0];
        let vol_j = vec![1.0];
        let delta = fickian_shift(&pos_i, &pos_j, &w_ij, &vol_j, 1.0, 1.0, 2);
        // delta[0] = 1.0 * (1.0 - 0.0) * 1.0 * 1.0 = 1.0
        assert!((delta[0] - 1.0).abs() < 1e-14);
        assert!(delta[1].abs() < 1e-14);
    }

    #[test]
    fn test_fickian_zero_dt() {
        let pos_i = [0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0]];
        let w_ij = vec![1.0];
        let vol_j = vec![1.0];
        let delta = fickian_shift(&pos_i, &pos_j, &w_ij, &vol_j, 1.0, 0.0, 2);
        assert!(delta[0].abs() < 1e-15);
    }

    #[test]
    fn test_fickian_scales_with_diffusion_coeff() {
        let pos_i = [0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0]];
        let w_ij = vec![1.0];
        let vol_j = vec![1.0];
        let d1 = fickian_shift(&pos_i, &pos_j, &w_ij, &vol_j, 1.0, 1.0, 2);
        let d2 = fickian_shift(&pos_i, &pos_j, &w_ij, &vol_j, 2.0, 1.0, 2);
        assert!((d2[0] / d1[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_fickian_symmetric_neighbours_cancel() {
        let pos_i = [0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0], vec![-1.0, 0.0]];
        let w_ij = vec![1.0, 1.0];
        let vol_j = vec![1.0, 1.0];
        let delta = fickian_shift(&pos_i, &pos_j, &w_ij, &vol_j, 1.0, 1.0, 2);
        assert!(delta[0].abs() < 1e-14);
    }

    // ── tensile_correction ─────────────────────────────────────────────────

    #[test]
    fn test_tensile_correction_unit() {
        // w_ij = w_0 => R = epsilon * 1^n = epsilon
        let r = tensile_correction(1.0, 1.0, 0.2, 4);
        assert!((r - 0.2).abs() < 1e-14);
    }

    #[test]
    fn test_tensile_correction_zero_w() {
        let r = tensile_correction(0.0, 1.0, 0.2, 4);
        assert!(r.abs() < 1e-14);
    }

    #[test]
    fn test_tensile_correction_zero_w0() {
        let r = tensile_correction(1.0, 0.0, 0.2, 4);
        assert!(r.abs() < 1e-14);
    }

    #[test]
    fn test_tensile_correction_half_ratio() {
        // ratio = 0.5, n=4 => 0.2 * 0.5^4 = 0.2 * 0.0625 = 0.0125
        let r = tensile_correction(0.5, 1.0, 0.2, 4);
        assert!((r - 0.0125).abs() < 1e-14);
    }

    #[test]
    fn test_tensile_correction_non_negative() {
        for i in 0..10 {
            let ratio = i as f64 * 0.1;
            let r = tensile_correction(ratio, 1.0, 0.2, 4);
            assert!(r >= 0.0);
        }
    }

    // ── apply_shift ────────────────────────────────────────────────────────

    #[test]
    fn test_apply_shift_moves_particle() {
        let mut pos = vec![vec![0.0, 0.0]];
        let shifts = vec![vec![0.1, 0.2]];
        apply_shift(&mut pos, &shifts, 1.0, 2);
        assert!((pos[0][0] - 0.1).abs() < 1e-14);
        assert!((pos[0][1] - 0.2).abs() < 1e-14);
    }

    #[test]
    fn test_apply_shift_clamped() {
        let mut pos = vec![vec![0.0, 0.0]];
        // shift of magnitude sqrt(2) = 1.414, max_shift = 1.0 => clamped
        let shifts = vec![vec![1.0, 1.0]];
        apply_shift(&mut pos, &shifts, 1.0, 2);
        let mag: f64 = pos[0].iter().map(|&p| p * p).sum::<f64>().sqrt();
        assert!(mag <= 1.0 + 1e-12);
    }

    #[test]
    fn test_apply_shift_small_not_clamped() {
        let mut pos = vec![vec![0.0, 0.0]];
        let shifts = vec![vec![0.01, 0.0]];
        apply_shift(&mut pos, &shifts, 1.0, 2);
        assert!((pos[0][0] - 0.01).abs() < 1e-14);
    }

    #[test]
    fn test_apply_shift_multiple_particles() {
        let mut pos = vec![vec![0.0], vec![1.0]];
        let shifts = vec![vec![0.1], vec![-0.1]];
        apply_shift(&mut pos, &shifts, 1.0, 1);
        assert!((pos[0][0] - 0.1).abs() < 1e-14);
        assert!((pos[1][0] - 0.9).abs() < 1e-14);
    }

    // ── shift_boundary_check ───────────────────────────────────────────────

    #[test]
    fn test_boundary_check_no_violation() {
        let mut pos = vec![vec![0.5, 0.5]];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let n = shift_boundary_check(&mut pos, &lo, &hi, 2);
        assert_eq!(n, 0);
        assert!((pos[0][0] - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_boundary_check_clamps_lo() {
        let mut pos = vec![vec![-0.1, 0.5]];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let n = shift_boundary_check(&mut pos, &lo, &hi, 2);
        assert_eq!(n, 1);
        assert!((pos[0][0] - 0.0).abs() < 1e-14);
    }

    #[test]
    fn test_boundary_check_clamps_hi() {
        let mut pos = vec![vec![1.5, 0.5]];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let n = shift_boundary_check(&mut pos, &lo, &hi, 2);
        assert_eq!(n, 1);
        assert!((pos[0][0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_boundary_check_multiple_particles() {
        let mut pos = vec![vec![-0.1, 0.5], vec![0.5, 1.5], vec![0.3, 0.3]];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let n = shift_boundary_check(&mut pos, &lo, &hi, 2);
        assert_eq!(n, 2);
    }

    #[test]
    fn test_boundary_check_on_boundary_not_clamped() {
        let mut pos = vec![vec![0.0, 1.0]];
        let lo = vec![0.0, 0.0];
        let hi = vec![1.0, 1.0];
        let n = shift_boundary_check(&mut pos, &lo, &hi, 2);
        assert_eq!(n, 0);
    }
}
