// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FFT-based spectral homogenization (Moulinec–Suquet 1998).
//!
//! Computes the effective fourth-order stiffness tensor `C_eff` of a periodic
//! two-phase (or arbitrary voxel-wise) representative volume element (RVE) by
//! solving the periodic Lippmann–Schwinger equation in Fourier space.
//!
//! Two solvers are provided:
//!
//! * [`fft::lippmann_schwinger`] — the classical Moulinec–Suquet basic
//!   fixed-point scheme (robust, `O(contrast)` iterations).
//! * [`fft::eyre_milton`] — the Eyre–Milton (1999) accelerated polarization
//!   scheme (`O(√contrast)` iterations at high phase contrast).
//!
//! The effective stiffness is assembled by [`compute_c_eff`], which applies six
//! independent unit macroscopic-strain loadings and reads off the volume-average
//! stress for each.
//!
//! # Conventions
//!
//! All 6-component vectors use the Voigt ordering of
//! [`crate::elastic::LinearElastic::stress_strain_matrix_3d`]:
//! `[xx, yy, zz, yz, xz, xy]`, with **engineering shear strain** (γ = 2ε) on the
//! shear entries. Stiffness fields are flat `Vec<[[f64;6];6]>` of length
//! `nx·ny·nz` indexed row-major as `ix·ny·nz + iy·nz + iz`.
//!
//! # Example
//!
//! ```
//! use oxiphysics_materials::elastic::LinearElastic;
//! use oxiphysics_materials::homogenization::compute_c_eff;
//!
//! // Homogeneous 2×2×2 RVE: C_eff must equal the single-phase stiffness.
//! let c = LinearElastic::new(70.0e9, 0.3).stress_strain_matrix_3d();
//! let field = vec![c; 8];
//! let (lambda0, mu0) = {
//!     let m = LinearElastic::new(70.0e9, 0.3);
//!     (m.lame_lambda(), m.lame_mu())
//! };
//! let result = compute_c_eff(&field, 2, 2, 2, (lambda0, mu0), 1e-10, 100).unwrap();
//! assert!(result.converged);
//! assert!((result.c_eff[0][0] - c[0][0]).abs() / c[0][0] < 1e-9);
//! ```

pub mod fft;

/// Errors raised by the FFT spectral homogenization solvers.
#[derive(Debug, thiserror::Error)]
pub enum HomogenizationError {
    /// The fixed-point iteration did not reach the requested tolerance.
    #[error(
        "FFT homogenization did not converge in {max_iter} iterations; final residual {residual:.3e}"
    )]
    NotConverged {
        /// Iteration budget that was exhausted.
        max_iter: usize,
        /// Final normalized equilibrium residual.
        residual: f64,
    },
    /// The grid description is inconsistent (zero size, length mismatch, …).
    #[error("Invalid grid: {0}")]
    InvalidGrid(String),
}

/// Result of an effective-stiffness homogenization run.
#[derive(Debug, Clone)]
pub struct HomogenizationResult {
    /// Effective fourth-order stiffness in 6×6 Voigt form.
    pub c_eff: [[f64; 6]; 6],
    /// Number of solver iterations used for each of the six unit loadings.
    pub iterations: [usize; 6],
    /// Whether every loading converged within the iteration budget.
    pub converged: bool,
}

impl HomogenizationResult {
    /// Effective bulk modulus from the isotropic part of `C_eff`.
    ///
    /// `K = (C₁₁ + C₂₂ + C₃₃ + 2(C₁₂ + C₁₃ + C₂₃)) / 9`.
    pub fn effective_bulk_modulus(&self) -> f64 {
        let c = &self.c_eff;
        let trace_normal = c[0][0] + c[1][1] + c[2][2];
        let off = c[0][1] + c[0][2] + c[1][2];
        (trace_normal + 2.0 * off) / 9.0
    }

    /// Effective shear modulus estimate from the deviatoric part of `C_eff`.
    ///
    /// Uses the isotropic projection
    /// `G = (C₁₁+C₂₂+C₃₃ − (C₁₂+C₁₃+C₂₃) + 3(C₄₄+C₅₅+C₆₆)) / 15`.
    pub fn effective_shear_modulus(&self) -> f64 {
        let c = &self.c_eff;
        let normal = c[0][0] + c[1][1] + c[2][2];
        let off = c[0][1] + c[0][2] + c[1][2];
        let shear = c[3][3] + c[4][4] + c[5][5];
        (normal - off + 3.0 * shear) / 15.0
    }
}

/// Selectable iteration scheme for [`compute_c_eff_with`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scheme {
    /// Classical Moulinec–Suquet basic fixed-point scheme.
    Basic,
    /// Eyre–Milton accelerated polarization scheme.
    Accelerated,
}

/// Assemble the effective stiffness tensor using the basic Moulinec–Suquet
/// scheme.
///
/// Convenience wrapper over [`compute_c_eff_with`] with [`Scheme::Basic`].
pub fn compute_c_eff(
    c_field: &[[[f64; 6]; 6]],
    nx: usize,
    ny: usize,
    nz: usize,
    ref_moduli: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> Result<HomogenizationResult, HomogenizationError> {
    compute_c_eff_with(
        c_field,
        [nx, ny, nz],
        ref_moduli,
        tol,
        max_iter,
        Scheme::Basic,
    )
}

/// Assemble the effective stiffness tensor `C_eff` by running six unit
/// macroscopic-strain loadings through the chosen spectral solver.
///
/// Column `m` of `C_eff` is the volume-average stress `⟨σ⟩` produced by the
/// loading `ε̄ = e_m` (the `m`-th engineering-strain unit vector).
///
/// # Arguments
/// * `c_field`    — flat voxel stiffness field, length `nx·ny·nz`.
/// * `dims`       — grid dimensions `[nx, ny, nz]`.
/// * `ref_moduli` — reference medium `(λ₀, μ₀)`; for best convergence pick the
///   midpoint of the phase stiffness spectrum (see [`reference_from_field`]).
/// * `tol`        — relative equilibrium-residual tolerance.
/// * `max_iter`   — per-loading iteration budget.
/// * `scheme`     — basic or accelerated iteration.
pub fn compute_c_eff_with(
    c_field: &[[[f64; 6]; 6]],
    dims: [usize; 3],
    ref_moduli: (f64, f64),
    tol: f64,
    max_iter: usize,
    scheme: Scheme,
) -> Result<HomogenizationResult, HomogenizationError> {
    let [nx, ny, nz] = dims;
    let n_total = nx.checked_mul(ny).and_then(|v| v.checked_mul(nz));
    let n_total = match n_total {
        Some(v) if v > 0 => v,
        _ => {
            return Err(HomogenizationError::InvalidGrid(format!(
                "grid {nx}×{ny}×{nz} has zero or overflowing voxel count"
            )));
        }
    };
    if c_field.len() != n_total {
        return Err(HomogenizationError::InvalidGrid(format!(
            "c_field length {} does not match grid {nx}×{ny}×{nz} = {n_total}",
            c_field.len()
        )));
    }

    let inv_v = 1.0 / n_total as f64;

    let mut c_eff = [[0.0_f64; 6]; 6];
    let mut iterations = [0_usize; 6];
    let mut converged = true;

    for m in 0..6 {
        let mut mean_strain = [0.0_f64; 6];
        mean_strain[m] = 1.0;

        // The inner solvers always return the final field together with the
        // convergence flag, so a slowly-converging loading still contributes a
        // best-effort averaged stress while flipping `converged` to false.
        let outcome = match scheme {
            Scheme::Basic => fft::lippmann_schwinger_inner(
                c_field,
                dims,
                mean_strain,
                ref_moduli,
                tol,
                max_iter,
            )?,
            Scheme::Accelerated => {
                fft::eyre_milton_inner(c_field, dims, mean_strain, ref_moduli, tol, max_iter)?
            }
        };

        if !outcome.converged {
            converged = false;
        }
        iterations[m] = outcome.iterations;

        // Column m = ⟨σ⟩ = (1/V) Σ_x C(x):ε(x).
        for (c, eps) in c_field.iter().zip(outcome.field.iter()) {
            let sigma = matvec6(c, eps);
            for (row, s) in sigma.iter().enumerate() {
                c_eff[row][m] += s;
            }
        }
        for row_arr in c_eff.iter_mut() {
            row_arr[m] *= inv_v;
        }
        // Guard: if the averaged column contains non-finite values, the solver
        // diverged.  Return an error rather than silently propagating NaN.
        for row_arr in c_eff.iter() {
            if !row_arr[m].is_finite() {
                return Err(HomogenizationError::NotConverged {
                    max_iter,
                    residual: outcome.residual,
                });
            }
        }
    }

    Ok(HomogenizationResult {
        c_eff,
        iterations,
        converged,
    })
}

/// Apply a 6×6 Voigt stiffness to an engineering-strain Voigt vector.
#[inline]
fn matvec6(c: &[[f64; 6]; 6], eps: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, oi) in out.iter_mut().enumerate() {
        let row = &c[i];
        *oi = row[0] * eps[0]
            + row[1] * eps[1]
            + row[2] * eps[2]
            + row[3] * eps[3]
            + row[4] * eps[4]
            + row[5] * eps[5];
    }
    out
}

/// Suggest reference Lamé moduli `(λ₀, μ₀)` for a two-phase problem from the
/// phase moduli, using the spectral midpoint that minimizes the iteration count.
///
/// `λ₀ = (λ_min + λ_max)/2`, `μ₀ = (μ_min + μ_max)/2` over the two phases.
pub fn reference_from_two_phase(lambda_a: f64, mu_a: f64, lambda_b: f64, mu_b: f64) -> (f64, f64) {
    let lambda0 = 0.5 * (lambda_a.min(lambda_b) + lambda_a.max(lambda_b));
    let mu0 = 0.5 * (mu_a.min(mu_b) + mu_a.max(mu_b));
    (lambda0, mu0)
}

/// Suggest reference Lamé moduli `(λ₀, μ₀)` for use with the **basic**
/// Moulinec–Suquet scheme on a two-phase problem.
///
/// The basic fixed-point map has spectral radius < 1 only when the reference
/// medium is stiffer than **both** phases in both λ and μ.  Using the maximum
/// of the two-phase values (with a small safety margin) guarantees convergence
/// for any phase contrast.
///
/// Contrast with [`reference_from_two_phase`], which computes the midpoint —
/// the correct choice for the Eyre–Milton accelerated scheme but not for the
/// basic scheme at high contrast.
pub fn reference_stiffer_than_both(
    lambda_a: f64,
    mu_a: f64,
    lambda_b: f64,
    mu_b: f64,
) -> (f64, f64) {
    let safety = 1.01_f64;
    let lambda0 = lambda_a.max(lambda_b) * safety;
    let mu0 = mu_a.max(mu_b) * safety;
    (lambda0, mu0)
}

/// Suggest reference Lamé moduli `(λ₀, μ₀)` for the **Eyre–Milton** accelerated
/// scheme on a two-phase problem with high phase contrast.
///
/// The optimal reference for Eyre–Milton is the **geometric mean** of the phase
/// moduli, which gives `O(√contrast)` convergence.  Contrast with
/// [`reference_from_two_phase`], which uses the arithmetic mean and performs
/// poorly at contrast > 100.
///
/// `λ₀ = √(λ_min · λ_max)`, `μ₀ = √(μ_min · μ_max)` over the two phases.
/// If either minimum is zero (or negative), falls back to the arithmetic mean.
pub fn reference_geometric_mean(lambda_a: f64, mu_a: f64, lambda_b: f64, mu_b: f64) -> (f64, f64) {
    let lambda_min = lambda_a.min(lambda_b);
    let lambda_max = lambda_a.max(lambda_b);
    let mu_min = mu_a.min(mu_b);
    let mu_max = mu_a.max(mu_b);
    let lambda0 = if lambda_min > 0.0 {
        (lambda_min * lambda_max).sqrt()
    } else {
        0.5 * (lambda_min + lambda_max)
    };
    let mu0 = if mu_min > 0.0 {
        (mu_min * mu_max).sqrt()
    } else {
        0.5 * (mu_min + mu_max)
    };
    (lambda0, mu0)
}

/// Suggest reference Lamé moduli from an explicit voxel stiffness field by
/// scanning the diagonal entries for the min/max P-wave and shear moduli.
///
/// Returns `(λ₀, μ₀)` at the midpoint of the observed `C₁₁` and `C₄₄` ranges,
/// a robust default reference for arbitrary microstructures.
pub fn reference_from_field(c_field: &[[[f64; 6]; 6]]) -> (f64, f64) {
    let mut m_min = f64::INFINITY;
    let mut m_max = f64::NEG_INFINITY;
    let mut mu_min = f64::INFINITY;
    let mut mu_max = f64::NEG_INFINITY;
    for c in c_field {
        let m_mod = c[0][0]; // λ + 2μ
        let mu = c[3][3]; // μ
        m_min = m_min.min(m_mod);
        m_max = m_max.max(m_mod);
        mu_min = mu_min.min(mu);
        mu_max = mu_max.max(mu);
    }
    if !m_min.is_finite() {
        return (1.0, 1.0);
    }
    let mu0 = 0.5 * (mu_min + mu_max);
    let m0 = 0.5 * (m_min + m_max);
    let lambda0 = m0 - 2.0 * mu0;
    (lambda0, mu0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::combination::{
        hashin_shtrikman_bulk_lower, hashin_shtrikman_bulk_upper, mori_tanaka_homogenization,
    };
    use crate::elastic::LinearElastic;

    /// Build a flat homogeneous stiffness field on an `n³` grid.
    fn homogeneous_field(c: [[f64; 6]; 6], n: usize) -> Vec<[[f64; 6]; 6]> {
        vec![c; n * n * n]
    }

    /// Build a two-phase checkerboard field at (roughly) the target volume
    /// fraction `vf` of the inclusion phase, using a sphere-packing-free
    /// deterministic assignment: a voxel belongs to the inclusion if its
    /// linear index falls in the first `round(vf·N)` of a bit-reversal-like
    /// scramble, which spreads the inclusion uniformly through the cell.
    fn two_phase_field(
        c_matrix: [[f64; 6]; 6],
        c_incl: [[f64; 6]; 6],
        n: usize,
        vf: f64,
    ) -> (Vec<[[f64; 6]; 6]>, f64) {
        let total = n * n * n;
        let n_incl = ((vf * total as f64).round() as usize).min(total);
        let mut field = vec![c_matrix; total];
        // Golden-ratio low-discrepancy stride over [0,total) spreads the
        // inclusion voxels uniformly. For a power-of-two grid the odd stride is
        // coprime with `total`, so each index is visited exactly once before any
        // repeat — no clustering and no need to track visited voxels.
        let stride = (((total as f64) * 0.618_033_988_75) as usize | 1) % total;
        let stride = if stride == 0 { 1 } else { stride };
        let mut idx = 0usize;
        for _ in 0..n_incl {
            field[idx] = c_incl;
            idx = (idx + stride) % total;
        }
        let actual_vf = n_incl as f64 / total as f64;
        (field, actual_vf)
    }

    #[test]
    fn test_homogeneous_medium() {
        let mat = LinearElastic::new(70.0e9, 0.3);
        let c = mat.stress_strain_matrix_3d();
        let field = homogeneous_field(c, 4);
        let ref_moduli = (mat.lame_lambda(), mat.lame_mu());

        let result = compute_c_eff(&field, 4, 4, 4, ref_moduli, 1e-12, 200)
            .expect("homogeneous solve must succeed");
        assert!(result.converged, "homogeneous case must converge");

        let mut max_err = 0.0_f64;
        for (ci_row, eff_row) in c.iter().zip(result.c_eff.iter()) {
            for (c_val, eff_val) in ci_row.iter().zip(eff_row.iter()) {
                let scale = c_val.abs().max(1.0);
                let err = (eff_val - c_val).abs() / scale;
                max_err = max_err.max(err);
            }
        }
        assert!(
            max_err < 1e-9,
            "homogeneous C_eff error {max_err:.3e} exceeds 1e-9"
        );
    }

    #[test]
    fn test_homogeneous_medium_accelerated() {
        let mat = LinearElastic::new(120.0e9, 0.25);
        let c = mat.stress_strain_matrix_3d();
        let field = homogeneous_field(c, 4);
        let ref_moduli = (mat.lame_lambda(), mat.lame_mu());

        let result = compute_c_eff_with(
            &field,
            [4, 4, 4],
            ref_moduli,
            1e-12,
            200,
            Scheme::Accelerated,
        )
        .expect("accelerated homogeneous solve must succeed");
        assert!(result.converged);
        let mut max_err = 0.0_f64;
        for (ci_row, eff_row) in c.iter().zip(result.c_eff.iter()) {
            for (c_val, eff_val) in ci_row.iter().zip(eff_row.iter()) {
                let scale = c_val.abs().max(1.0);
                max_err = max_err.max((eff_val - c_val).abs() / scale);
            }
        }
        assert!(max_err < 1e-9, "accelerated homogeneous err {max_err:.3e}");
    }

    #[test]
    fn test_two_phase_hs_bounds() {
        // Matrix: K_m via E_m, nu_m; inclusion 10× stiffer in bulk.
        let matrix = LinearElastic::new(3.0e9, 0.35); // soft polymer-like
        let incl = LinearElastic::new(30.0e9, 0.25); // stiffer phase, K ratio ~10
        let c_m = matrix.stress_strain_matrix_3d();
        let c_i = incl.stress_strain_matrix_3d();

        let n = 8;
        let vf = 0.3;
        let (field, actual_vf) = two_phase_field(c_m, c_i, n, vf);

        let ref_moduli = reference_stiffer_than_both(
            matrix.lame_lambda(),
            matrix.lame_mu(),
            incl.lame_lambda(),
            incl.lame_mu(),
        );

        let result = compute_c_eff(&field, n, n, n, ref_moduli, 1e-6, 500)
            .expect("two-phase solve must succeed");

        let k_eff = result.effective_bulk_modulus();
        let g_eff = result.effective_shear_modulus();

        // Phase moduli.
        let k_m = matrix.bulk_modulus();
        let g_m = matrix.shear_modulus();
        let k_i = incl.bulk_modulus();
        let g_i = incl.shear_modulus();

        // HS bounds for bulk (v1 = matrix fraction). Soft phase = matrix.
        let v_matrix = 1.0 - actual_vf;
        let k_hs_lower = hashin_shtrikman_bulk_lower(v_matrix, k_m, g_m, k_i);
        let k_hs_upper = hashin_shtrikman_bulk_upper(v_matrix, k_m, k_i, g_i);
        let (lo, hi) = (k_hs_lower.min(k_hs_upper), k_hs_lower.max(k_hs_upper));

        // Allow a small discretization slack on a coarse 8³ grid.
        let slack = 0.05 * (hi - lo).max(hi.abs() * 0.02);
        assert!(
            k_eff >= lo - slack && k_eff <= hi + slack,
            "K_eff {k_eff:.4e} outside HS bounds [{lo:.4e}, {hi:.4e}] (slack {slack:.2e}), vf={actual_vf}"
        );

        // Shear must be bracketed by the soft/stiff phase shear moduli at least.
        assert!(
            g_eff >= g_m * 0.9 && g_eff <= g_i * 1.1,
            "G_eff {g_eff:.4e} not between phase shear moduli [{g_m:.4e}, {g_i:.4e}]"
        );
    }

    #[test]
    fn test_eshelby_dilute_limit() {
        // Low volume fraction: C_eff bulk modulus should approach the
        // Mori-Tanaka / Eshelby dilute estimate.
        let matrix = LinearElastic::new(3.0e9, 0.35);
        let incl = LinearElastic::new(30.0e9, 0.25);
        let c_m = matrix.stress_strain_matrix_3d();
        let c_i = incl.stress_strain_matrix_3d();

        let n = 8;
        let vf = 0.05;
        let (field, actual_vf) = two_phase_field(c_m, c_i, n, vf);

        let ref_moduli = reference_stiffer_than_both(
            matrix.lame_lambda(),
            matrix.lame_mu(),
            incl.lame_lambda(),
            incl.lame_mu(),
        );
        let result = compute_c_eff(&field, n, n, n, ref_moduli, 1e-6, 500)
            .expect("dilute solve must succeed");
        let k_eff = result.effective_bulk_modulus();

        let (k_mt, _g_mt) = mori_tanaka_homogenization(
            actual_vf,
            incl.bulk_modulus(),
            incl.shear_modulus(),
            matrix.bulk_modulus(),
            matrix.shear_modulus(),
        );

        let rel = (k_eff - k_mt).abs() / k_mt;
        assert!(
            rel < 0.10,
            "dilute K_eff {k_eff:.4e} vs Mori-Tanaka {k_mt:.4e}, rel err {rel:.3} > 10%"
        );
    }

    #[test]
    fn test_invalid_grid_rejected() {
        let mat = LinearElastic::new(70.0e9, 0.3);
        let c = mat.stress_strain_matrix_3d();
        let field = vec![c; 8];
        // Wrong dimensions (2×2×3 = 12 ≠ 8).
        let err = compute_c_eff(&field, 2, 2, 3, (1.0, 1.0), 1e-6, 10);
        assert!(matches!(err, Err(HomogenizationError::InvalidGrid(_))));
    }

    #[test]
    fn test_high_contrast_accelerated_converges() {
        // Phase contrast 1000 in stiffness; accelerated scheme must converge.
        let matrix = LinearElastic::new(1.0e6, 0.3);
        let incl = LinearElastic::new(1.0e9, 0.3); // 1000× stiffer
        let c_m = matrix.stress_strain_matrix_3d();
        let c_i = incl.stress_strain_matrix_3d();

        let n = 8;
        let (field, _vf) = two_phase_field(c_m, c_i, n, 0.2);
        let ref_moduli = reference_geometric_mean(
            matrix.lame_lambda(),
            matrix.lame_mu(),
            incl.lame_lambda(),
            incl.lame_mu(),
        );

        let result = compute_c_eff_with(
            &field,
            [n, n, n],
            ref_moduli,
            1e-4,
            500,
            Scheme::Accelerated,
        )
        .expect("accelerated high-contrast solve must succeed");
        assert!(
            result.converged,
            "accelerated scheme failed to converge at contrast 1000: iters {:?}",
            result.iterations
        );
        // Effective bulk modulus must lie between the two phase bulk moduli.
        let k_eff = result.effective_bulk_modulus();
        assert!(
            k_eff > matrix.bulk_modulus() * 0.9 && k_eff < incl.bulk_modulus() * 1.1,
            "high-contrast K_eff {k_eff:.4e} out of phase range"
        );
    }
}



