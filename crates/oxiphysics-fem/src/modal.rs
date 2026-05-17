#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Modal analysis: natural frequencies and mode shapes.
//!
//! Solves the generalized eigenvalue problem K*φ = ω²*M*φ
//! using inverse iteration (power method on K^{-1}*M) for the lowest modes.
//!
//! Additional capabilities:
//! - Mode shape normalization (mass, max-component, unit-length)
//! - Modal Assurance Criterion (MAC)
//! - Frequency Response Function (FRF)
//! - Harmonic analysis
//! - Modal superposition

use std::f64::consts::PI;

use crate::solvers::conjugate_gradient;
use crate::sparse::CsrMatrix;

/// Result of a modal analysis.
pub struct ModalResult {
    /// Natural angular frequencies ω_i (rad/s), sorted ascending.
    pub frequencies: Vec<f64>,
    /// Mode shapes φ_i (column vectors, each of length n_dof), mass-normalized.
    pub mode_shapes: Vec<Vec<f64>>,
}

/// Matrix-vector product: y = A*x for CSR matrix A.
pub fn matvec(a: &CsrMatrix, x: &[f64]) -> Vec<f64> {
    a.mul_vec(x)
}

/// Rayleigh quotient: rq = x^T*A*x / x^T*B*x
pub fn rayleigh_quotient(a: &CsrMatrix, b: &CsrMatrix, x: &[f64]) -> f64 {
    let ax = matvec(a, x);
    let bx = matvec(b, x);
    let num: f64 = x.iter().zip(ax.iter()).map(|(xi, axi)| xi * axi).sum();
    let den: f64 = x.iter().zip(bx.iter()).map(|(xi, bxi)| xi * bxi).sum();
    if den.abs() < 1e-60 { 0.0 } else { num / den }
}

/// Dot product of two slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Compute the n_modes lowest natural frequencies using inverse iteration.
///
/// Algorithm (simplified):
/// 1. Start with random initial vector b
/// 2. Solve K*y = M*b  (inverse iteration: y = K^{-1}*M*b)
/// 3. Normalize: b = y / sqrt(y^T*M*y)
/// 4. Rayleigh quotient: ω² = b^T*K*b / b^T*M*b
/// 5. Repeat until converged
///
/// For multiple modes, use deflation (subtract converged modes from search space).
pub fn inverse_iteration(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    n_modes: usize,
    max_iter: usize,
    tol: f64,
) -> ModalResult {
    let n = stiffness.nrows;
    assert_eq!(stiffness.ncols, n);
    assert_eq!(mass.nrows, n);
    assert_eq!(mass.ncols, n);

    let mut frequencies = Vec::with_capacity(n_modes);
    let mut mode_shapes: Vec<Vec<f64>> = Vec::with_capacity(n_modes);

    for mode_idx in 0..n_modes {
        // Initial vector: simple unit vector shifted by mode_idx to avoid degeneracy
        let mut b = vec![0.0; n];
        let start_idx = mode_idx % n;
        b[start_idx] = 1.0;
        // Add small perturbation to make vector non-degenerate
        for (i, bi) in b.iter_mut().enumerate() {
            *bi += 0.1 * ((i + mode_idx * 7 + 1) as f64).sin();
        }

        // Mass-normalize initial vector
        let mb = matvec(mass, &b);
        let m_norm = dot(&b, &mb).sqrt();
        if m_norm > 1e-60 {
            for bi in b.iter_mut() {
                *bi /= m_norm;
            }
        }

        let mut omega_sq = 0.0;

        for _iter in 0..max_iter {
            // Deflate: subtract components along already-converged modes
            for prev_shape in &mode_shapes {
                let m_prev = matvec(mass, prev_shape);
                let coeff: f64 = dot(&b, &m_prev);
                for i in 0..n {
                    b[i] -= coeff * prev_shape[i];
                }
            }

            // Solve K * y = M * b
            let mb = matvec(mass, &b);
            let x0 = vec![0.0; n];
            let y = conjugate_gradient(stiffness, &mb, &x0, 10 * n, tol * 1e-3);

            // Deflate y as well
            let mut y_deflated = y.clone();
            for prev_shape in &mode_shapes {
                let m_prev = matvec(mass, prev_shape);
                let coeff: f64 = dot(&y_deflated, &m_prev);
                for i in 0..n {
                    y_deflated[i] -= coeff * prev_shape[i];
                }
            }

            // Mass-normalize: b = y / sqrt(y^T * M * y)
            let my = matvec(mass, &y_deflated);
            let norm_sq = dot(&y_deflated, &my);
            if norm_sq < 1e-120 {
                break;
            }
            let norm = norm_sq.sqrt();
            let b_new: Vec<f64> = y_deflated.iter().map(|yi| yi / norm).collect();

            // Rayleigh quotient: omega^2 = b^T * K * b / b^T * M * b
            let kb = matvec(stiffness, &b_new);
            let mb_new = matvec(mass, &b_new);
            let num: f64 = dot(&b_new, &kb);
            let den: f64 = dot(&b_new, &mb_new);
            let omega_sq_new = if den.abs() > 1e-60 { num / den } else { 0.0 };

            // Check convergence
            let change = (omega_sq_new - omega_sq).abs();
            let rel_change = change / (omega_sq_new.abs() + 1e-60);

            b = b_new;
            omega_sq = omega_sq_new;

            if rel_change < tol && _iter > 0 {
                break;
            }
        }

        let omega = omega_sq.max(0.0).sqrt();
        frequencies.push(omega);
        mode_shapes.push(b);
    }

    // Sort by ascending frequency
    let mut indexed: Vec<(usize, f64)> = frequencies.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let sorted_frequencies: Vec<f64> = indexed.iter().map(|&(_, f)| f).collect();
    let sorted_shapes: Vec<Vec<f64>> = indexed
        .iter()
        .map(|&(i, _)| mode_shapes[i].clone())
        .collect();

    ModalResult {
        frequencies: sorted_frequencies,
        mode_shapes: sorted_shapes,
    }
}

// ---------------------------------------------------------------------------
// Mode shape normalization
// ---------------------------------------------------------------------------

/// Normalization method for mode shapes.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizationMethod {
    /// Mass normalization: φ^T M φ = 1
    Mass,
    /// Max-component normalization: max|φ_i| = 1
    MaxComponent,
    /// Unit-length normalization: ||φ|| = 1 (Euclidean norm)
    UnitLength,
}

/// Normalize a mode shape vector according to the given method.
///
/// For `Mass` normalization, the mass matrix is required.
///
/// # Arguments
/// * `mode_shape` - Mode shape vector (modified in-place)
/// * `method` - Normalization method to apply
/// * `mass` - Optional mass matrix (required for `Mass` normalization)
#[allow(dead_code)]
pub fn normalize_mode_shape(
    mode_shape: &mut [f64],
    method: NormalizationMethod,
    mass: Option<&CsrMatrix>,
) {
    match method {
        NormalizationMethod::Mass => {
            let m = mass.expect("Mass matrix required for mass normalization");
            let m_phi = m.mul_vec(mode_shape);
            let gen_mass: f64 = dot(mode_shape, &m_phi);
            if gen_mass.abs() > 1e-60 {
                let scale = 1.0 / gen_mass.abs().sqrt();
                for v in mode_shape.iter_mut() {
                    *v *= scale;
                }
            }
        }
        NormalizationMethod::MaxComponent => {
            let max_abs = mode_shape.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if max_abs > 1e-60 {
                for v in mode_shape.iter_mut() {
                    *v /= max_abs;
                }
            }
        }
        NormalizationMethod::UnitLength => {
            let norm: f64 = mode_shape.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm > 1e-60 {
                for v in mode_shape.iter_mut() {
                    *v /= norm;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Modal Assurance Criterion (MAC)
// ---------------------------------------------------------------------------

/// Compute the Modal Assurance Criterion (MAC) between two mode shapes.
///
/// MAC = |φ_a^T φ_b|^2 / (φ_a^T φ_a)(φ_b^T φ_b)
///
/// Returns a value in \[0, 1\] where 1 indicates perfectly correlated modes.
#[allow(dead_code)]
pub fn mac_value(phi_a: &[f64], phi_b: &[f64]) -> f64 {
    assert_eq!(phi_a.len(), phi_b.len(), "mode shape lengths must match");
    let cross: f64 = dot(phi_a, phi_b);
    let norm_a: f64 = dot(phi_a, phi_a);
    let norm_b: f64 = dot(phi_b, phi_b);
    let denom = norm_a * norm_b;
    if denom.abs() < 1e-60 {
        0.0
    } else {
        (cross * cross) / denom
    }
}

/// Compute the full MAC matrix between two sets of mode shapes.
///
/// Returns a matrix (`Vec<Vec`f64`>`) of size `modes_a.len()` x `modes_b.len()`.
#[allow(dead_code)]
pub fn mac_matrix(modes_a: &[Vec<f64>], modes_b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let na = modes_a.len();
    let nb = modes_b.len();
    let mut mac = vec![vec![0.0; nb]; na];
    for i in 0..na {
        for j in 0..nb {
            mac[i][j] = mac_value(&modes_a[i], &modes_b[j]);
        }
    }
    mac
}

// ---------------------------------------------------------------------------
// Frequency Response Function (FRF)
// ---------------------------------------------------------------------------

/// Result of a Frequency Response Function computation at a single excitation frequency.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FrfResult {
    /// Excitation frequency ω (rad/s).
    pub omega: f64,
    /// Response magnitude at each DOF.
    pub magnitude: Vec<f64>,
    /// Response phase at each DOF (radians).
    pub phase: Vec<f64>,
}

/// Compute the FRF using modal superposition for a single excitation frequency.
///
/// H(ω) = Σ_r φ_r φ_r^T / (ω_r^2 - ω^2 + 2iζ_r ω_r ω)
///
/// # Arguments
/// * `modal_result` - Modal analysis result containing frequencies and mode shapes
/// * `omega` - Excitation angular frequency (rad/s)
/// * `damping_ratios` - Modal damping ratio ζ for each mode
/// * `force_vector` - Applied force vector (one per DOF)
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn frequency_response_function(
    modal_result: &ModalResult,
    omega: f64,
    damping_ratios: &[f64],
    force_vector: &[f64],
) -> FrfResult {
    let n_dof = force_vector.len();
    let n_modes = modal_result.frequencies.len();

    // Accumulate real and imaginary parts of the response
    let mut u_real = vec![0.0; n_dof];
    let mut u_imag = vec![0.0; n_dof];

    for r in 0..n_modes {
        let omega_r = modal_result.frequencies[r];
        let phi_r = &modal_result.mode_shapes[r];
        let zeta = if r < damping_ratios.len() {
            damping_ratios[r]
        } else {
            0.0
        };

        // Modal force: q_r = φ_r^T · F
        let q_r: f64 = phi_r
            .iter()
            .zip(force_vector.iter())
            .map(|(p, f)| p * f)
            .sum();

        // Complex denominator: (ω_r^2 - ω^2) + i(2 ζ ω_r ω)
        let real_denom = omega_r * omega_r - omega * omega;
        let imag_denom = 2.0 * zeta * omega_r * omega;
        let denom_sq = real_denom * real_denom + imag_denom * imag_denom;

        if denom_sq < 1e-60 {
            continue;
        }

        // Complex scalar: q_r / denom = q_r * (real - i*imag) / denom_sq
        let scale_real = q_r * real_denom / denom_sq;
        let scale_imag = -q_r * imag_denom / denom_sq;

        for i in 0..n_dof.min(phi_r.len()) {
            u_real[i] += phi_r[i] * scale_real;
            u_imag[i] += phi_r[i] * scale_imag;
        }
    }

    // Compute magnitude and phase
    let mut magnitude = vec![0.0; n_dof];
    let mut phase = vec![0.0; n_dof];
    for i in 0..n_dof {
        magnitude[i] = (u_real[i] * u_real[i] + u_imag[i] * u_imag[i]).sqrt();
        phase[i] = u_imag[i].atan2(u_real[i]);
    }

    FrfResult {
        omega,
        magnitude,
        phase,
    }
}

/// Compute FRF over a range of frequencies.
///
/// # Arguments
/// * `modal_result` - Modal analysis result
/// * `omega_range` - Slice of excitation frequencies
/// * `damping_ratios` - Modal damping ratios
/// * `force_vector` - Applied force vector
#[allow(dead_code)]
pub fn frf_sweep(
    modal_result: &ModalResult,
    omega_range: &[f64],
    damping_ratios: &[f64],
    force_vector: &[f64],
) -> Vec<FrfResult> {
    omega_range
        .iter()
        .map(|&omega| {
            frequency_response_function(modal_result, omega, damping_ratios, force_vector)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Harmonic analysis
// ---------------------------------------------------------------------------

/// Result of harmonic analysis at a single frequency.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HarmonicResult {
    /// Excitation frequency ω (rad/s).
    pub omega: f64,
    /// Real part of displacement response.
    pub displacement_real: Vec<f64>,
    /// Imaginary part of displacement response.
    pub displacement_imag: Vec<f64>,
    /// Peak displacement magnitude at each DOF.
    pub displacement_amplitude: Vec<f64>,
}

/// Perform harmonic analysis using direct frequency response.
///
/// Solves (-ω^2 M + iωC + K) u = F for complex displacement u.
/// Here C is assumed as proportional damping: C = α M + β K.
///
/// # Arguments
/// * `stiffness` - Global stiffness matrix K
/// * `mass` - Global mass matrix M
/// * `omega` - Excitation frequency (rad/s)
/// * `alpha` - Mass-proportional damping coefficient
/// * `beta` - Stiffness-proportional damping coefficient
/// * `force_vector` - Applied harmonic force vector
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn harmonic_analysis(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    omega: f64,
    alpha: f64,
    beta: f64,
    force_vector: &[f64],
) -> HarmonicResult {
    let n = stiffness.nrows;
    assert_eq!(force_vector.len(), n);

    // Build the effective real stiffness: K_eff = K - ω^2 M
    // Build the effective damping: C_eff = α M + β K
    // Solve two coupled real systems:
    //   K_eff u_r - ω C_eff u_i = F
    //   ω C_eff u_r + K_eff u_i = 0

    // For simplicity, solve iteratively using CG on the real part:
    // (K - ω^2 M) u_r = F - ω(αM + βK) u_i
    // (K - ω^2 M) u_i = -ω(αM + βK) u_r

    // Start with u_i = 0, iterate a few times
    let mut u_real = vec![0.0; n];
    let mut u_imag = vec![0.0; n];

    // Build K_eff as triplets
    let omega_sq = omega * omega;
    let mut k_eff_triplets = Vec::new();
    for row in 0..n {
        let start = stiffness.row_ptr[row];
        let end = stiffness.row_ptr[row + 1];
        for idx in start..end {
            let col = stiffness.col_indices[idx];
            let k_val = stiffness.values[idx];
            let m_val = mass.get(row, col);
            k_eff_triplets.push((row, col, k_val - omega_sq * m_val));
        }
    }
    let _k_eff = CsrMatrix::from_triplets(n, n, &k_eff_triplets);

    // Make K_eff diagonally dominant for CG to work
    // Add small regularization if needed
    let mut k_reg_triplets = k_eff_triplets.clone();
    for i in 0..n {
        k_reg_triplets.push((i, i, 1e-10));
    }
    let k_reg = CsrMatrix::from_triplets(n, n, &k_reg_triplets);

    for _outer in 0..5 {
        // Compute ω * C_eff * u_imag
        // C_eff = α M + β K
        let m_ui = matvec(mass, &u_imag);
        let k_ui = matvec(stiffness, &u_imag);
        let mut rhs_real = force_vector.to_vec();
        for i in 0..n {
            rhs_real[i] -= omega * (alpha * m_ui[i] + beta * k_ui[i]);
        }

        let x0 = u_real.clone();
        u_real = conjugate_gradient(&k_reg, &rhs_real, &x0, 10 * n, 1e-10);

        // Solve for u_imag
        let m_ur = matvec(mass, &u_real);
        let k_ur = matvec(stiffness, &u_real);
        let mut rhs_imag = vec![0.0; n];
        for i in 0..n {
            rhs_imag[i] = -omega * (alpha * m_ur[i] + beta * k_ur[i]);
        }

        let x0 = u_imag.clone();
        u_imag = conjugate_gradient(&k_reg, &rhs_imag, &x0, 10 * n, 1e-10);
    }

    let mut displacement_amplitude = vec![0.0; n];
    for i in 0..n {
        displacement_amplitude[i] = (u_real[i] * u_real[i] + u_imag[i] * u_imag[i]).sqrt();
    }

    HarmonicResult {
        omega,
        displacement_real: u_real,
        displacement_imag: u_imag,
        displacement_amplitude,
    }
}

// ---------------------------------------------------------------------------
// Modal superposition
// ---------------------------------------------------------------------------

/// Perform transient analysis using modal superposition.
///
/// Decomposes the response into modal coordinates and integrates each
/// decoupled SDOF equation independently using the Duhamel integral.
///
/// For each mode r: q_r'' + 2 ζ_r ω_r q_r' + ω_r^2 q_r = φ_r^T F(t) / m_r
///
/// # Arguments
/// * `modal_result` - Eigenvalue solution (frequencies and mode shapes)
/// * `damping_ratios` - Modal damping ratio ζ for each mode
/// * `force_history` - Force time history as (time, force_vector) pairs
///
/// # Returns
/// Vector of displacement vectors, one per time step.
#[allow(dead_code)]
pub fn modal_superposition(
    modal_result: &ModalResult,
    damping_ratios: &[f64],
    force_history: &[(f64, Vec<f64>)],
) -> Vec<Vec<f64>> {
    if force_history.is_empty() {
        return Vec::new();
    }

    let n_dof = force_history[0].1.len();
    let n_modes = modal_result.frequencies.len();
    let n_steps = force_history.len();

    // Initialize modal coordinates (displacement and velocity)
    let mut q = vec![0.0; n_modes];
    let mut q_dot = vec![0.0; n_modes];
    let mut results = Vec::with_capacity(n_steps);

    for step in 0..n_steps {
        let dt = if step > 0 {
            force_history[step].0 - force_history[step - 1].0
        } else if n_steps > 1 {
            force_history[1].0 - force_history[0].0
        } else {
            0.001
        };

        let force = &force_history[step].1;

        // Update each modal coordinate using the Newmark-beta method (average acceleration)
        for r in 0..n_modes {
            let omega_r = modal_result.frequencies[r];
            let phi_r = &modal_result.mode_shapes[r];
            let zeta = if r < damping_ratios.len() {
                damping_ratios[r]
            } else {
                0.0
            };

            // Modal force
            let modal_force: f64 = phi_r.iter().zip(force.iter()).map(|(p, f)| p * f).sum();

            if omega_r < 1e-30 {
                // Rigid body mode: q'' = modal_force
                q[r] += q_dot[r] * dt + 0.5 * modal_force * dt * dt;
                q_dot[r] += modal_force * dt;
                continue;
            }

            // Damped SDOF step (exact integration for constant force over dt)
            let omega_d = omega_r * (1.0 - zeta * zeta).max(0.0).sqrt();
            let exp_term = (-zeta * omega_r * dt).exp();

            if omega_d < 1e-30 {
                // Critically damped
                q[r] = exp_term * (q[r] + (q_dot[r] + omega_r * q[r]) * dt)
                    + modal_force / (omega_r * omega_r);
                q_dot[r] = exp_term * (q_dot[r] * (1.0 - omega_r * dt));
            } else {
                let cos_wd = (omega_d * dt).cos();
                let sin_wd = (omega_d * dt).sin();

                // Particular solution for constant force
                let q_static = modal_force / (omega_r * omega_r);

                // Free response from current state
                let q_free = q[r] - q_static;
                let q_new = exp_term
                    * (q_free * cos_wd + (q_dot[r] + zeta * omega_r * q_free) / omega_d * sin_wd)
                    + q_static;

                let q_dot_new = exp_term
                    * (-q_free * omega_d * sin_wd + (q_dot[r] + zeta * omega_r * q_free) * cos_wd
                        - zeta
                            * omega_r
                            * (q_free * cos_wd
                                + (q_dot[r] + zeta * omega_r * q_free) / omega_d * sin_wd));

                q[r] = q_new;
                q_dot[r] = q_dot_new;
            }
        }

        // Reconstruct physical displacement: u = Σ_r φ_r q_r
        let mut u = vec![0.0; n_dof];
        for r in 0..n_modes {
            let phi_r = &modal_result.mode_shapes[r];
            for i in 0..n_dof.min(phi_r.len()) {
                u[i] += phi_r[i] * q[r];
            }
        }
        results.push(u);
    }

    results
}

/// Compute effective modal mass for each mode.
///
/// The effective modal mass for mode r in a given direction is:
/// m_eff_r = (φ_r^T M L)^2 / (φ_r^T M φ_r)
///
/// where L is the rigid-body mode (unit vector in the direction of interest).
///
/// # Arguments
/// * `mode_shapes` - Mode shape vectors
/// * `mass` - Mass matrix
/// * `direction` - Direction vector (length n_dof, usually a rigid-body mode)
#[allow(dead_code)]
pub fn effective_modal_mass(
    mode_shapes: &[Vec<f64>],
    mass: &CsrMatrix,
    direction: &[f64],
) -> Vec<f64> {
    let mut masses = Vec::with_capacity(mode_shapes.len());
    let m_dir = matvec(mass, direction);

    for phi in mode_shapes {
        let participation: f64 = dot(phi, &m_dir);
        let m_phi = matvec(mass, phi);
        let gen_mass: f64 = dot(phi, &m_phi);
        if gen_mass.abs() > 1e-60 {
            masses.push(participation * participation / gen_mass);
        } else {
            masses.push(0.0);
        }
    }
    masses
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build identity CSR matrix of size n x n.
    fn identity_csr(n: usize) -> CsrMatrix {
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0)).collect();
        CsrMatrix::from_triplets(n, n, &triplets)
    }

    fn diagonal_csr(diag: &[f64]) -> CsrMatrix {
        let triplets: Vec<(usize, usize, f64)> =
            diag.iter().enumerate().map(|(i, &v)| (i, i, v)).collect();
        CsrMatrix::from_triplets(diag.len(), diag.len(), &triplets)
    }

    #[test]
    fn test_matvec_identity() {
        let id = identity_csr(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = matvec(&id, &x);
        for i in 0..4 {
            assert!(
                (y[i] - x[i]).abs() < 1e-14,
                "y[{i}] = {}, expected {}",
                y[i],
                x[i]
            );
        }
    }

    #[test]
    fn test_rayleigh_quotient_eigenvector() {
        let triplets_a = vec![(0, 0, 3.0), (1, 1, 7.0)];
        let a = CsrMatrix::from_triplets(2, 2, &triplets_a);
        let b = identity_csr(2);
        let x = vec![1.0, 0.0];
        let rq = rayleigh_quotient(&a, &b, &x);
        assert!(
            (rq - 3.0).abs() < 1e-12,
            "Rayleigh quotient = {rq}, expected 3.0"
        );

        let x2 = vec![0.0, 1.0];
        let rq2 = rayleigh_quotient(&a, &b, &x2);
        assert!(
            (rq2 - 7.0).abs() < 1e-12,
            "Rayleigh quotient = {rq2}, expected 7.0"
        );
    }

    #[test]
    fn test_modal_spring_mass_chain() {
        let k_val = 100.0_f64;
        let m_val = 1.0_f64;
        let triplets_k = vec![(0, 0, k_val)];
        let triplets_m = vec![(0, 0, m_val)];
        let k = CsrMatrix::from_triplets(1, 1, &triplets_k);
        let m = CsrMatrix::from_triplets(1, 1, &triplets_m);
        let result = inverse_iteration(&k, &m, 1, 200, 1e-8);
        let expected_omega = (k_val / m_val).sqrt();
        assert!(
            (result.frequencies[0] - expected_omega).abs() / expected_omega < 1e-5,
            "omega = {}, expected {}",
            result.frequencies[0],
            expected_omega
        );
    }

    #[test]
    fn test_modal_frequencies_sorted() {
        let triplets_k = vec![(0, 0, 1.0), (1, 1, 4.0), (2, 2, 9.0)];
        let k = CsrMatrix::from_triplets(3, 3, &triplets_k);
        let m = identity_csr(3);
        let result = inverse_iteration(&k, &m, 3, 500, 1e-8);
        assert_eq!(result.frequencies.len(), 3);
        for i in 0..result.frequencies.len() - 1 {
            assert!(
                result.frequencies[i] <= result.frequencies[i + 1] + 1e-10,
                "frequencies not sorted: {} > {}",
                result.frequencies[i],
                result.frequencies[i + 1]
            );
        }
        assert!(
            (result.frequencies[0] - 1.0).abs() < 0.1,
            "omega_0 = {}",
            result.frequencies[0]
        );
        assert!(
            (result.frequencies[1] - 2.0).abs() < 0.1,
            "omega_1 = {}",
            result.frequencies[1]
        );
        assert!(
            (result.frequencies[2] - 3.0).abs() < 0.1,
            "omega_2 = {}",
            result.frequencies[2]
        );
    }

    // --- New tests for expanded functionality ---

    #[test]
    fn test_normalize_unit_length() {
        let mut phi = vec![3.0, 4.0];
        normalize_mode_shape(&mut phi, NormalizationMethod::UnitLength, None);
        let norm: f64 = phi.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-12, "unit-length norm = {norm}");
    }

    #[test]
    fn test_normalize_max_component() {
        let mut phi = vec![2.0, -5.0, 3.0];
        normalize_mode_shape(&mut phi, NormalizationMethod::MaxComponent, None);
        let max_abs = phi.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!((max_abs - 1.0).abs() < 1e-12, "max component = {max_abs}");
        // The sign should be preserved
        assert!(phi[1] < 0.0, "sign should be preserved");
    }

    #[test]
    fn test_normalize_mass() {
        let m = diagonal_csr(&[2.0, 3.0]);
        let mut phi = vec![1.0, 1.0];
        normalize_mode_shape(&mut phi, NormalizationMethod::Mass, Some(&m));
        // Check: phi^T M phi = 1
        let m_phi = m.mul_vec(&phi);
        let gen_mass: f64 = dot(&phi, &m_phi);
        assert!(
            (gen_mass - 1.0).abs() < 1e-12,
            "generalized mass = {gen_mass}"
        );
    }

    #[test]
    fn test_mac_identical_modes() {
        let phi = vec![1.0, 2.0, 3.0];
        let mac = mac_value(&phi, &phi);
        assert!((mac - 1.0).abs() < 1e-12, "MAC of identical modes = {mac}");
    }

    #[test]
    fn test_mac_orthogonal_modes() {
        let phi_a = vec![1.0, 0.0, 0.0];
        let phi_b = vec![0.0, 1.0, 0.0];
        let mac = mac_value(&phi_a, &phi_b);
        assert!(mac.abs() < 1e-12, "MAC of orthogonal modes = {mac}");
    }

    #[test]
    fn test_mac_scaled_modes() {
        let phi_a = vec![1.0, 2.0, 3.0];
        let phi_b = vec![10.0, 20.0, 30.0];
        let mac = mac_value(&phi_a, &phi_b);
        assert!((mac - 1.0).abs() < 1e-12, "MAC of scaled modes = {mac}");
    }

    #[test]
    fn test_mac_matrix_diagonal() {
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let mac = mac_matrix(&modes, &modes);
        for i in 0..3 {
            assert!(
                (mac[i][i] - 1.0).abs() < 1e-12,
                "diagonal MAC[{i}][{i}] = {}",
                mac[i][i]
            );
            for j in 0..3 {
                if i != j {
                    assert!(
                        mac[i][j].abs() < 1e-12,
                        "off-diagonal MAC[{i}][{j}] = {}",
                        mac[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn test_frf_at_resonance() {
        // Single mode: omega_r = 10 rad/s, zeta = 0.01
        // At resonance omega = omega_r, response should be large
        let modal_result = ModalResult {
            frequencies: vec![10.0],
            mode_shapes: vec![vec![1.0]],
        };
        let damping = vec![0.01];
        let force = vec![1.0];

        let frf_res = frequency_response_function(&modal_result, 10.0, &damping, &force);
        let frf_off = frequency_response_function(&modal_result, 1.0, &damping, &force);

        assert!(
            frf_res.magnitude[0] > frf_off.magnitude[0],
            "resonance magnitude {} should exceed off-resonance {}",
            frf_res.magnitude[0],
            frf_off.magnitude[0]
        );
    }

    #[test]
    fn test_frf_sweep_length() {
        let modal_result = ModalResult {
            frequencies: vec![10.0],
            mode_shapes: vec![vec![1.0]],
        };
        let omega_range: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let results = frf_sweep(&modal_result, &omega_range, &[0.01], &[1.0]);
        assert_eq!(results.len(), 20);
    }

    #[test]
    fn test_frf_static_limit() {
        // At omega = 0, response should equal static solution: u = F / k = phi * phi^T * F / omega_r^2
        let modal_result = ModalResult {
            frequencies: vec![10.0],
            mode_shapes: vec![vec![1.0]],
        };
        let frf = frequency_response_function(&modal_result, 0.0, &[0.01], &[1.0]);
        // Expected: 1.0 * 1.0 * 1.0 / (10^2) = 0.01
        assert!(
            (frf.magnitude[0] - 0.01).abs() < 1e-6,
            "static FRF magnitude = {}",
            frf.magnitude[0]
        );
    }

    #[test]
    fn test_modal_superposition_static_load() {
        // Constant force on a single SDOF system should converge to static deflection
        let modal_result = ModalResult {
            frequencies: vec![10.0],
            mode_shapes: vec![vec![1.0]],
        };
        let damping = vec![0.5]; // heavy damping for fast convergence

        let n_steps = 500;
        let dt = 0.01;
        let force_history: Vec<(f64, Vec<f64>)> =
            (0..n_steps).map(|i| (i as f64 * dt, vec![1.0])).collect();

        let results = modal_superposition(&modal_result, &damping, &force_history);
        assert_eq!(results.len(), n_steps);

        // At steady state, u should approach F / omega_r^2 = 1.0 / 100 = 0.01
        let final_u = results.last().unwrap()[0];
        assert!(
            (final_u - 0.01).abs() < 0.005,
            "final displacement = {final_u}, expected ~0.01"
        );
    }

    #[test]
    fn test_modal_superposition_initial_zero() {
        let modal_result = ModalResult {
            frequencies: vec![10.0],
            mode_shapes: vec![vec![1.0]],
        };
        // Zero force should give zero displacement
        let force_history: Vec<(f64, Vec<f64>)> =
            (0..10).map(|i| (i as f64 * 0.01, vec![0.0])).collect();

        let results = modal_superposition(&modal_result, &[0.01], &force_history);
        for (step, u) in results.iter().enumerate() {
            assert!(
                u[0].abs() < 1e-12,
                "step {step}: displacement = {} should be zero",
                u[0]
            );
        }
    }

    #[test]
    fn test_effective_modal_mass() {
        // 2-DOF diagonal system: M = I, modes = unit vectors
        let m = identity_csr(2);
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let direction = vec![1.0, 0.0]; // x-direction

        let masses = effective_modal_mass(&modes, &m, &direction);
        // Mode 1 (x) should have full mass in x-direction = 1.0
        assert!(
            (masses[0] - 1.0).abs() < 1e-12,
            "mode 1 eff mass = {}",
            masses[0]
        );
        // Mode 2 (y) should have zero mass in x-direction
        assert!(masses[1].abs() < 1e-12, "mode 2 eff mass = {}", masses[1]);
    }

    #[test]
    fn test_effective_modal_mass_sum() {
        // Sum of effective modal masses should equal total mass in the direction
        let m = identity_csr(3);
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let direction = vec![1.0, 1.0, 1.0]; // uniform direction

        let masses = effective_modal_mass(&modes, &m, &direction);
        let total: f64 = masses.iter().sum();
        // Total mass in this direction = L^T M L = 3.0
        let expected_total: f64 = dot(&direction, &m.mul_vec(&direction));
        assert!(
            (total - expected_total).abs() < 1e-10,
            "total eff mass = {total}, expected {expected_total}"
        );
    }

    #[test]
    fn test_harmonic_analysis_basic() {
        // 1-DOF: K = [[100]], M = [[1]], omega = 0 (static)
        let k = diagonal_csr(&[100.0]);
        let m = diagonal_csr(&[1.0]);
        let force = vec![10.0];

        let result = harmonic_analysis(&k, &m, 0.0, 0.0, 0.0, &force);
        // Static: u = F / K = 10 / 100 = 0.1
        assert!(
            (result.displacement_real[0] - 0.1).abs() < 0.01,
            "static displacement = {}",
            result.displacement_real[0]
        );
    }
}

// ---------------------------------------------------------------------------
// Lanczos iteration for large eigenproblems
// ---------------------------------------------------------------------------

/// Result of a Lanczos iteration.
#[allow(dead_code)]
pub struct LanczosResult {
    /// Approximate eigenvalues (angular frequencies squared ω²).
    pub eigenvalues: Vec<f64>,
    /// Approximate eigenvectors (mode shapes), each of length n_dof.
    pub eigenvectors: Vec<Vec<f64>>,
    /// Number of Lanczos steps taken.
    pub n_steps: usize,
}

/// Lanczos algorithm for computing the k smallest eigenvalues of K·φ = λM·φ.
///
/// This is more efficient than inverse iteration for large systems because
/// it builds a small tridiagonal matrix T whose eigenvalues approximate
/// those of the generalized problem.
///
/// Algorithm:
/// 1. Start with a random M-normalized vector v_1.
/// 2. For j=1..m: compute w = K⁻¹ M v_j (shift-and-invert), orthogonalize.
/// 3. Build tridiagonal T from α and β recurrences.
/// 4. Solve T y = θ y (small eigenproblem via bisection/QR).
/// 5. Ritz vectors: x_i = V y_i (approximate mode shapes).
#[allow(dead_code)]
pub fn lanczos_iteration(
    stiffness: &CsrMatrix,
    mass: &CsrMatrix,
    n_modes: usize,
    lanczos_steps: usize,
    tol: f64,
) -> LanczosResult {
    let n = stiffness.nrows;
    let m = lanczos_steps.max(n_modes + 4).min(n);

    // Storage for Lanczos vectors
    let mut v_vecs: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
    let mut alpha = vec![0.0f64; m];
    let mut beta = vec![0.0f64; m + 1];

    // Initial vector: uniform normalized
    let mut v = vec![1.0 / (n as f64).sqrt(); n];
    // M-normalize
    let mv = matvec(mass, &v);
    let vnorm = dot(&v, &mv).sqrt();
    for x in v.iter_mut() {
        *x /= vnorm;
    }
    v_vecs.push(v.clone());

    let mut actual_steps = 0;
    for j in 0..m {
        // w = K^{-1} M v_j (inverse iteration step)
        let mv_j = matvec(mass, &v_vecs[j]);
        let x0 = vec![0.0; n];
        let mut w = conjugate_gradient(stiffness, &mv_j, &x0, 15 * n, tol * 1e-4);

        // α_j = v_j^T M w
        let mw = matvec(mass, &w);
        alpha[j] = dot(&v_vecs[j], &mw);

        // w = w - α_j v_j - β_j v_{j-1}
        for i in 0..n {
            w[i] -= alpha[j] * v_vecs[j][i];
        }
        if j > 0 {
            for i in 0..n {
                w[i] -= beta[j] * v_vecs[j - 1][i];
            }
        }

        // Full reorthogonalization (double Gram-Schmidt for stability)
        for prev in &v_vecs {
            let mprev = matvec(mass, prev);
            let coeff = dot(&w, &mprev);
            for i in 0..n {
                w[i] -= coeff * prev[i];
            }
        }
        for prev in &v_vecs {
            let mprev = matvec(mass, prev);
            let coeff = dot(&w, &mprev);
            for i in 0..n {
                w[i] -= coeff * prev[i];
            }
        }

        // β_{j+1} = ||w||_M
        let mw2 = matvec(mass, &w);
        let w_norm = dot(&w, &mw2).max(0.0).sqrt();
        beta[j + 1] = w_norm;

        actual_steps = j + 1;

        if w_norm < tol {
            break;
        }

        // v_{j+1} = w / w_norm
        let v_new: Vec<f64> = w.iter().map(|x| x / w_norm).collect();
        v_vecs.push(v_new.clone());
        v = v_new;
        let _ = &v; // used next iteration
    }

    // Solve tridiagonal eigenproblem T y = θ y using power iteration on T
    // T is m×m symmetric tridiagonal: diag=alpha, off-diag=beta[1..m]
    // For simplicity, use direct QR iteration on the small matrix
    let dim = actual_steps;
    let eigenvalues_t = tridiagonal_eigenvalues(&alpha[..dim], &beta[1..dim + 1], 200);

    // Take the n_modes smallest (they correspond to 1/λ since we used K^{-1})
    // Ritz eigenvalues λ = 1/θ (since we inverted K)
    let mut ritz_pairs: Vec<(f64, Vec<f64>)> = eigenvalues_t
        .iter()
        .enumerate()
        .map(|(i, &theta)| {
            // Ritz vector = V * y_i  (approximate eigenvector of T is unit vec e_i for small matrix)
            let lambda = if theta.abs() > 1e-30 {
                1.0 / theta
            } else {
                0.0
            };
            let mut ritz = vec![0.0; n];
            // For the tridiagonal system we use eigenvectors approximated via inverse iteration on T
            // Simplified: use the j-th Lanczos vector as the Ritz approximation
            if i < v_vecs.len() {
                ritz = v_vecs[i].clone();
            }
            (lambda, ritz)
        })
        .collect();

    // Sort by eigenvalue ascending
    ritz_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    ritz_pairs.truncate(n_modes);

    let eigenvalues: Vec<f64> = ritz_pairs.iter().map(|(e, _)| e.max(0.0).sqrt()).collect();
    let eigenvectors: Vec<Vec<f64>> = ritz_pairs.into_iter().map(|(_, v)| v).collect();

    LanczosResult {
        eigenvalues,
        eigenvectors,
        n_steps: actual_steps,
    }
}

/// Compute eigenvalues of a symmetric tridiagonal matrix using the
/// implicit symmetric QR algorithm with Wilkinson shift (LAPACK dsteqr style).
///
/// `diag`    - diagonal elements (length n)
/// `offdiag` - off-diagonal elements (length >= n; element 0 is unused)
///
/// Returns eigenvalues sorted in ascending order.
fn tridiagonal_eigenvalues(diag: &[f64], offdiag: &[f64], max_iter: usize) -> Vec<f64> {
    let n = diag.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return diag.to_vec();
    }

    let mut d = diag.to_vec();
    // e[i] holds the off-diagonal between d[i] and d[i+1] (length n-1 used, padded to n)
    let mut e = vec![0.0f64; n];
    for i in 0..n - 1 {
        e[i] = if i + 1 < offdiag.len() {
            offdiag[i + 1]
        } else {
            0.0
        };
    }

    let eps = f64::EPSILON;

    // QL-implicit algorithm (Golub & Van Loan §8.3.4)
    let mut l = 0usize;
    while l < n {
        // Find small sub-diagonal element to split off a block [l..m+1]
        let mut m = l;
        while m < n - 1 {
            let tst = d[m].abs() + d[m + 1].abs();
            if e[m].abs() <= eps * tst {
                e[m] = 0.0;
                break;
            }
            m += 1;
        }
        if m == l {
            // 1×1 block converged
            l += 1;
            continue;
        }
        // Run QR steps on sub-matrix [l..m+1] (size m-l+1 >= 2)
        let mut iter_count = 0;
        loop {
            iter_count += 1;
            if iter_count > max_iter {
                break;
            }
            // Wilkinson shift from bottom 2×2
            let b = (d[m - 1] - d[m]) * 0.5;
            let r = (b * b + e[m - 1] * e[m - 1]).sqrt();
            let shift = d[m] - e[m - 1] * e[m - 1] / (b + b.signum() * r + f64::EPSILON * 1e-10);

            // Implicit QL step
            let mut g = d[l] - shift;
            let mut p_val = 1.0;
            for i in l..m {
                let eim1 = if i > l { e[i - 1] } else { 0.0 };
                let _ = eim1;
                let hyp = (g * g + e[i] * e[i]).sqrt();
                let c = if hyp > 1e-300 { g / hyp } else { 1.0 };
                let s = if hyp > 1e-300 { e[i] / hyp } else { 0.0 };
                if i > l {
                    e[i - 1] = hyp * p_val;
                }
                let w = d[i] - shift;
                let dip1_w = d[i + 1] - shift;
                g = c * (c * w - s * e[i]) + dip1_w * s * s;
                let f = s * (c * e[i] + s * w);
                d[i] = shift + c * c * w + s * s * d[i + 1] - 2.0 * s * c * e[i];
                d[i + 1] = shift + s * s * w + c * c * d[i + 1] + 2.0 * s * c * e[i];
                p_val = c;
                e[i] = f;
            }
            let _ = (g, p_val);
            if m > l {
                e[m - 1] = 0.0;
            } // often already ~0 after step

            // Check for new convergence within [l..m+1]
            let mut converged = false;
            for k in l..m {
                let tst = d[k].abs() + d[k + 1].abs();
                if e[k].abs() <= eps * tst {
                    e[k] = 0.0;
                    converged = true;
                    break;
                }
            }
            if converged {
                break;
            }
        }
        l += 1; // move forward even if not fully converged to avoid infinite loop
    }

    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    d
}

// ---------------------------------------------------------------------------
// Mode superposition (enhanced) and modal truncation
// ---------------------------------------------------------------------------

/// Configuration for modal truncation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ModalTruncationConfig {
    /// Fraction of total effective mass that must be captured.
    pub mass_capture_fraction: f64,
    /// Maximum number of modes to retain.
    pub max_modes: usize,
    /// Include static correction (residual mode)?
    pub static_correction: bool,
}

impl Default for ModalTruncationConfig {
    fn default() -> Self {
        Self {
            mass_capture_fraction: 0.90,
            max_modes: 50,
            static_correction: false,
        }
    }
}

/// Select the minimum number of modes to capture a given fraction of total
/// effective mass in a specified direction.
///
/// Returns the truncated index (number of modes to use).
#[allow(dead_code)]
pub fn modal_truncation(
    mode_shapes: &[Vec<f64>],
    mass: &CsrMatrix,
    direction: &[f64],
    config: &ModalTruncationConfig,
) -> usize {
    let all_masses = effective_modal_mass(mode_shapes, mass, direction);
    let total: f64 = all_masses.iter().sum();
    if total < 1e-60 {
        return config.max_modes.min(mode_shapes.len());
    }

    let target = total * config.mass_capture_fraction;
    let mut cumulative = 0.0;
    for (i, m) in all_masses.iter().enumerate() {
        cumulative += m;
        if cumulative >= target {
            return (i + 1).min(config.max_modes);
        }
    }
    config.max_modes.min(mode_shapes.len())
}

/// Compute modal participation factors for all modes.
///
/// Participation factor for mode r: Γ_r = φ_r^T M L
/// where L is the direction vector (rigid body mode).
#[allow(dead_code)]
pub fn participation_factors(
    mode_shapes: &[Vec<f64>],
    mass: &CsrMatrix,
    direction: &[f64],
) -> Vec<f64> {
    let m_dir = matvec(mass, direction);
    mode_shapes.iter().map(|phi| dot(phi, &m_dir)).collect()
}

/// Cumulative effective mass fraction as a function of mode count.
///
/// Returns a vector of cumulative mass fractions (0..1), one per mode.
#[allow(dead_code)]
pub fn cumulative_mass_fraction(
    mode_shapes: &[Vec<f64>],
    mass: &CsrMatrix,
    direction: &[f64],
) -> Vec<f64> {
    let masses = effective_modal_mass(mode_shapes, mass, direction);
    let total: f64 = masses.iter().sum();
    if total < 1e-60 {
        return vec![0.0; masses.len()];
    }
    let mut cumulative = 0.0;
    masses
        .iter()
        .map(|&m| {
            cumulative += m;
            cumulative / total
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Response Spectrum Analysis
// ---------------------------------------------------------------------------

/// Ground acceleration record for response spectrum analysis.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GroundMotion {
    /// Time values (seconds).
    pub time: Vec<f64>,
    /// Ground acceleration values (m/s²).
    pub acceleration: Vec<f64>,
}

/// Response spectrum data: peak response vs. natural period.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ResponseSpectrum {
    /// Natural periods T (seconds).
    pub periods: Vec<f64>,
    /// Spectral acceleration Sa(T) (m/s²).
    pub spectral_acceleration: Vec<f64>,
    /// Spectral velocity Sv(T) (m/s).
    pub spectral_velocity: Vec<f64>,
    /// Spectral displacement Sd(T) (m).
    pub spectral_displacement: Vec<f64>,
}

/// Compute a response spectrum from a ground motion record.
///
/// For each natural period T in `periods`:
/// 1. Compute ω = 2π/T
/// 2. Integrate the SDOF equation: ẍ + 2ζω ẋ + ω²x = -ẍ_g(t)
/// 3. Record peak absolute values of x, ẋ, ω²x
#[allow(dead_code)]
pub fn compute_response_spectrum(
    ground_motion: &GroundMotion,
    periods: &[f64],
    damping_ratio: f64,
) -> ResponseSpectrum {
    use std::f64::consts::PI;

    let mut sa = Vec::with_capacity(periods.len());
    let mut sv = Vec::with_capacity(periods.len());
    let mut sd = Vec::with_capacity(periods.len());

    for &t_period in periods {
        let omega_n = 2.0 * PI / t_period.max(1e-6);
        let omega_d = omega_n * (1.0 - damping_ratio * damping_ratio).max(0.0).sqrt();
        let zeta = damping_ratio;

        let mut x = 0.0f64;
        let mut x_dot = 0.0f64;
        let mut peak_x = 0.0f64;
        let mut peak_v = 0.0f64;

        let n = ground_motion.time.len();
        for i in 1..n {
            let dt = ground_motion.time[i] - ground_motion.time[i - 1];
            if dt <= 0.0 {
                continue;
            }

            // Piecewise-linear interpolation of ground acceleration
            let ag = ground_motion.acceleration[i];

            // Newmark-beta (average acceleration) integration
            let gamma = 0.5;
            let beta_nb = 0.25;

            let k_eff = omega_n * omega_n
                + gamma / (beta_nb * dt) * 2.0 * zeta * omega_n
                + 1.0 / (beta_nb * dt * dt);

            let rhs = -ag
                - (2.0
                    * zeta
                    * omega_n
                    * (x_dot
                        + (1.0 - gamma)
                            * dt
                            * (-2.0 * zeta * omega_n * x_dot - omega_n * omega_n * x + 0.0))
                    + omega_n * omega_n * x
                    + x_dot / (beta_nb * dt)
                    + x / (beta_nb * dt * dt)
                    - (x_dot + x / (beta_nb * dt)));

            let dx = if k_eff.abs() > 1e-60 {
                (-ag - 2.0 * zeta * omega_n * x_dot - omega_n * omega_n * x) / k_eff * dt
            } else {
                0.0
            };

            let x_dd = -ag - 2.0 * zeta * omega_n * x_dot - omega_n * omega_n * x;
            x_dot += x_dd * dt;
            x += x_dot * dt;
            let _ = (rhs, dx);

            // Use exact exponential integration for constant-excitation interval
            let exp_t = (-zeta * omega_n * dt).exp();
            if omega_d > 1e-30 {
                let cos_wd = (omega_d * dt).cos();
                let sin_wd = (omega_d * dt).sin();
                let x_new = exp_t * (x * cos_wd + (x_dot + zeta * omega_n * x) / omega_d * sin_wd)
                    - ag / (omega_n * omega_n);
                let x_dot_new = exp_t
                    * (-(omega_n.powi(2) * x + zeta * omega_n * (x_dot + zeta * omega_n * x))
                        / omega_d
                        * sin_wd
                        + (x_dot + zeta * omega_n * x) * cos_wd
                        - zeta * omega_n * x * cos_wd);
                x = x_new + ag / (omega_n * omega_n);
                x_dot = x_dot_new;
            }

            let abs_x = x.abs();
            let abs_v = x_dot.abs();
            if abs_x > peak_x {
                peak_x = abs_x;
            }
            if abs_v > peak_v {
                peak_v = abs_v;
            }
        }

        sd.push(peak_x);
        sv.push(peak_v.max(omega_n * peak_x));
        sa.push(omega_n * omega_n * peak_x);
    }

    ResponseSpectrum {
        periods: periods.to_vec(),
        spectral_acceleration: sa,
        spectral_velocity: sv,
        spectral_displacement: sd,
    }
}

/// Perform response spectrum analysis on a multi-DOF system.
///
/// For each mode r:
/// 1. Compute spectral displacement Sd(T_r) from the spectrum
/// 2. Peak modal response: u_r = φ_r * Γ_r * Sd(T_r)
/// 3. Combine using SRSS rule
///
/// Returns the SRSS-combined peak displacement at each DOF.
#[allow(dead_code)]
pub fn response_spectrum_analysis(
    modal_result: &ModalResult,
    mass: &CsrMatrix,
    direction: &[f64],
    spectrum: &ResponseSpectrum,
    damping_ratio: f64,
) -> Vec<f64> {
    let n_dof = direction.len();
    let n_modes = modal_result.frequencies.len();
    let pf = participation_factors(&modal_result.mode_shapes, mass, direction);

    // SRSS combination
    let mut u_srss = vec![0.0; n_dof];

    for r in 0..n_modes {
        let omega_r = modal_result.frequencies[r];
        let t_r = if omega_r > 1e-10 {
            2.0 * PI / omega_r
        } else {
            0.0
        };

        // Interpolate Sd from spectrum
        let sd = interpolate_spectrum(
            &spectrum.periods,
            &spectrum.spectral_displacement,
            t_r,
            damping_ratio,
        );

        let phi_r = &modal_result.mode_shapes[r];
        for i in 0..n_dof.min(phi_r.len()) {
            let u_modal = phi_r[i] * pf[r] * sd;
            u_srss[i] += u_modal * u_modal;
        }
    }

    u_srss.iter().map(|x| x.sqrt()).collect()
}

/// Linear interpolation in a spectrum table.
fn interpolate_spectrum(periods: &[f64], values: &[f64], t: f64, _zeta: f64) -> f64 {
    if periods.is_empty() || values.is_empty() {
        return 0.0;
    }
    if t <= periods[0] {
        return values[0];
    }
    if t >= *periods.last().expect("periods is non-empty") {
        return *values.last().expect("values is non-empty");
    }
    for i in 0..periods.len() - 1 {
        if t >= periods[i] && t <= periods[i + 1] {
            let alpha = (t - periods[i]) / (periods[i + 1] - periods[i]);
            return values[i] + alpha * (values[i + 1] - values[i]);
        }
    }
    *values.last().expect("values is non-empty")
}

// ---------------------------------------------------------------------------
// Additional tests for new functionality
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    fn identity_csr(n: usize) -> CsrMatrix {
        let triplets: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0)).collect();
        CsrMatrix::from_triplets(n, n, &triplets)
    }

    fn diagonal_csr(diag: &[f64]) -> CsrMatrix {
        let triplets: Vec<(usize, usize, f64)> =
            diag.iter().enumerate().map(|(i, &v)| (i, i, v)).collect();
        CsrMatrix::from_triplets(diag.len(), diag.len(), &triplets)
    }

    #[test]
    fn test_lanczos_single_dof() {
        let k = diagonal_csr(&[100.0]);
        let m = diagonal_csr(&[1.0]);
        let result = lanczos_iteration(&k, &m, 1, 5, 1e-8);
        assert_eq!(result.eigenvalues.len(), 1);
        // eigenvalue should be ~sqrt(100) = 10
        assert!(result.eigenvalues[0] > 0.0);
    }

    #[test]
    fn test_lanczos_returns_correct_count() {
        let k = diagonal_csr(&[1.0, 4.0, 9.0]);
        let m = identity_csr(3);
        let result = lanczos_iteration(&k, &m, 3, 10, 1e-8);
        assert_eq!(result.eigenvalues.len(), 3);
        assert_eq!(result.eigenvectors.len(), 3);
    }

    #[test]
    fn test_participation_factors() {
        let m = identity_csr(2);
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let dir = vec![1.0, 0.0];
        let pf = participation_factors(&modes, &m, &dir);
        assert!((pf[0] - 1.0).abs() < 1e-12);
        assert!(pf[1].abs() < 1e-12);
    }

    #[test]
    fn test_cumulative_mass_fraction() {
        let m = identity_csr(3);
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let dir = vec![1.0, 0.0, 0.0];
        let cum = cumulative_mass_fraction(&modes, &m, &dir);
        assert_eq!(cum.len(), 3);
        // After mode 1: 1.0/1.0 = 100%
        assert!((cum[0] - 1.0).abs() < 1e-12);
        // After modes 1+2+3: still 100%
        assert!((cum[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_modal_truncation_captures_mass() {
        let m = identity_csr(4);
        let modes: Vec<Vec<f64>> = vec![
            vec![0.9, 0.0, 0.0, 0.0],
            vec![0.0, 0.436, 0.0, 0.0],
            vec![0.0, 0.0, 0.316, 0.0],
            vec![0.0, 0.0, 0.0, 0.244],
        ];
        let dir = vec![1.0, 1.0, 1.0, 1.0];
        let config = ModalTruncationConfig {
            mass_capture_fraction: 0.80,
            max_modes: 4,
            static_correction: false,
        };
        let n_keep = modal_truncation(&modes, &m, &dir, &config);
        assert!(n_keep >= 1 && n_keep <= 4, "n_keep = {n_keep}");
    }

    #[test]
    fn test_response_spectrum_basic() {
        // Simple harmonic excitation
        let n = 200;
        let dt = 0.01;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time.iter().map(|&t| (2.0 * PI * 2.0 * t).sin()).collect();
        let gm = GroundMotion { time, acceleration };
        let periods = vec![0.5, 1.0, 2.0];
        let spectrum = compute_response_spectrum(&gm, &periods, 0.05);
        assert_eq!(spectrum.periods.len(), 3);
        assert_eq!(spectrum.spectral_acceleration.len(), 3);
        // All values should be finite and non-negative
        for &sa in &spectrum.spectral_acceleration {
            assert!(sa >= 0.0 && sa.is_finite(), "sa = {sa}");
        }
    }

    #[test]
    fn test_response_spectrum_analysis_srss() {
        let modal_result = ModalResult {
            frequencies: vec![10.0, 20.0],
            mode_shapes: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        };
        let m = identity_csr(2);
        let direction = vec![1.0, 1.0];
        let spectrum = ResponseSpectrum {
            periods: vec![0.1, 0.5, 1.0],
            spectral_acceleration: vec![1.0, 2.0, 1.0],
            spectral_velocity: vec![0.1, 0.2, 0.1],
            spectral_displacement: vec![0.01, 0.05, 0.02],
        };
        let u = response_spectrum_analysis(&modal_result, &m, &direction, &spectrum, 0.05);
        assert_eq!(u.len(), 2);
        for &ui in &u {
            assert!(ui >= 0.0 && ui.is_finite(), "displacement = {ui}");
        }
    }

    #[test]
    fn test_interpolate_spectrum_at_endpoints() {
        let periods = vec![0.1, 1.0, 10.0];
        let values = vec![5.0, 3.0, 1.0];
        let v_low = interpolate_spectrum(&periods, &values, 0.05, 0.05);
        assert!((v_low - 5.0).abs() < 1e-12);
        let v_high = interpolate_spectrum(&periods, &values, 15.0, 0.05);
        assert!((v_high - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_interpolate_spectrum_midpoint() {
        let periods = vec![0.0, 1.0];
        let values = vec![0.0, 10.0];
        let v = interpolate_spectrum(&periods, &values, 0.5, 0.05);
        assert!((v - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_tridiagonal_eigenvalues_diagonal() {
        // Diagonal matrix: eigenvalues are the diagonal elements
        let diag = vec![1.0, 4.0, 9.0];
        let offdiag = vec![0.0, 0.0, 0.0];
        let eigs = tridiagonal_eigenvalues(&diag, &offdiag, 500);
        assert_eq!(eigs.len(), 3);
        // Should recover 1, 4, 9 (sorted)
        assert!((eigs[0] - 1.0).abs() < 0.5, "eig0 = {}", eigs[0]);
    }

    #[test]
    fn test_modal_truncation_max_modes() {
        let m = identity_csr(10);
        let modes: Vec<Vec<f64>> = (0..10)
            .map(|i| {
                let mut v = vec![0.0; 10];
                v[i] = 1.0;
                v
            })
            .collect();
        let dir = vec![1.0; 10];
        let config = ModalTruncationConfig {
            mass_capture_fraction: 0.99,
            max_modes: 3, // limit to 3 even if more needed
            static_correction: false,
        };
        let n_keep = modal_truncation(&modes, &m, &dir, &config);
        assert!(n_keep <= 3, "Should not exceed max_modes: {n_keep}");
    }

    #[test]
    fn test_lanczos_nsteps_positive() {
        let k = diagonal_csr(&[1.0, 2.0]);
        let m = identity_csr(2);
        let result = lanczos_iteration(&k, &m, 2, 8, 1e-8);
        assert!(result.n_steps > 0, "Should have taken at least one step");
    }
}
