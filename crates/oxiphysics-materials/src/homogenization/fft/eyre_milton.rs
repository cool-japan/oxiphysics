// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eyre–Milton (1999) accelerated polarization scheme for the periodic
//! Lippmann–Schwinger equation.
//!
//! The basic Moulinec–Suquet scheme converges at a rate governed by the phase
//! contrast: roughly `O(contrast)` iterations. The Eyre–Milton accelerated
//! scheme reformulates the problem on the **polarization field**
//! `p(x) = (C(x) − C⁰):ε(x)` and reaches the same fixed point in
//! `O(√contrast)` iterations by choosing the reference stiffness C⁰ at the
//! geometric/arithmetic midpoint of the local stiffness spectrum and using an
//! over-relaxation factor.
//!
//! # Scheme
//!
//! With an isotropic reference C⁰ = (λ₀, μ₀) and relaxation factor `ω`:
//!
//! 1. ε(x) = ε̄ − Γ⁰ * p(x)            (strain from polarization; DC enforces ε̄)
//! 2. q(x) = (C(x) − C⁰):ε(x)          (target polarization)
//! 3. p(x) ← p(x) + ω·(q(x) − p(x))     (relaxed update)
//!
//! repeated until the equilibrium residual `‖ξ·σ̂(ξ)‖` falls below `tol`.
//!
//! The reference modulus is taken at the midpoint of the bounds `λ⁻, λ⁺` on the
//! eigenvalues of `C(x)` (Eyre–Milton's optimal choice), and the over-relaxation
//! factor `ω = 2` realizes the accelerated convergence. Step 1 uses the same
//! `oxifft` transforms and Green operator as the basic scheme.

use super::green_operator::{apply_gamma0, voigt_to_pair};
use super::lippmann_schwinger::LsOutcome;
use crate::homogenization::HomogenizationError;

/// Apply a 6×6 Voigt matrix to a Voigt vector.
#[inline]
fn matvec6(c: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, oi) in out.iter_mut().enumerate() {
        let row = &c[i];
        *oi = row[0] * v[0]
            + row[1] * v[1]
            + row[2] * v[2]
            + row[3] * v[3]
            + row[4] * v[4]
            + row[5] * v[5];
    }
    out
}

/// Wrapped integer frequency for FFT index `i` on a grid of size `n`.
#[inline]
fn wrapped_freq(i: usize, n: usize) -> f64 {
    if i <= n / 2 {
        i as f64
    } else {
        i as f64 - n as f64
    }
}

/// Build the isotropic reference 6×6 stiffness from Lamé parameters.
fn isotropic_c(lambda0: f64, mu0: f64) -> [[f64; 6]; 6] {
    let mut c = [[0.0_f64; 6]; 6];
    let diag = lambda0 + 2.0 * mu0;
    c[0][0] = diag;
    c[1][1] = diag;
    c[2][2] = diag;
    c[0][1] = lambda0;
    c[0][2] = lambda0;
    c[1][0] = lambda0;
    c[1][2] = lambda0;
    c[2][0] = lambda0;
    c[2][1] = lambda0;
    c[3][3] = mu0;
    c[4][4] = mu0;
    c[5][5] = mu0;
    c
}

/// Solve the periodic Lippmann–Schwinger problem by the Eyre–Milton accelerated
/// polarization scheme.
///
/// # Arguments
/// * `c_field`     — flat voxel stiffness field, length `nx*ny*nz`.
/// * `dims`        — grid dimensions `[nx, ny, nz]`.
/// * `mean_strain` — prescribed macroscopic strain ε̄ (engineering-strain Voigt).
/// * `ref_moduli`  — reference medium `(λ₀, μ₀)`. For best acceleration this
///   should sit near the midpoint of the phase stiffness spectrum; callers using
///   [`crate::homogenization::compute_c_eff`] receive an auto-selected value.
/// * `tol`         — relative equilibrium-residual tolerance.
/// * `max_iter`    — maximum number of iterations.
///
/// Returns the converged strain field and the iteration count, or
/// [`HomogenizationError::NotConverged`] if the residual never reached `tol`.
pub fn eyre_milton(
    c_field: &[[[f64; 6]; 6]],
    dims: [usize; 3],
    mean_strain: [f64; 6],
    ref_moduli: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<[f64; 6]>, usize), HomogenizationError> {
    let outcome = eyre_milton_inner(c_field, dims, mean_strain, ref_moduli, tol, max_iter)?;
    if outcome.converged {
        Ok((outcome.field, outcome.iterations))
    } else {
        Err(HomogenizationError::NotConverged {
            max_iter,
            residual: outcome.residual,
        })
    }
}

/// Core sweep returning the final field, residual and convergence flag.
pub(crate) fn eyre_milton_inner(
    c_field: &[[[f64; 6]; 6]],
    dims: [usize; 3],
    mean_strain: [f64; 6],
    ref_moduli: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> Result<LsOutcome, HomogenizationError> {
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

    let (lambda0, mu0) = ref_moduli;
    let c0 = isotropic_c(lambda0, mu0);
    let n_f = n_total as f64;

    // Over-relaxation factor (Eyre–Milton accelerated convergence).
    let omega = 2.0_f64;

    // Polarization field p(x) = (C(x) − C⁰):ε(x), initialized from ε = ε̄.
    let mut polar = vec![[0.0_f64; 6]; n_total];
    for (idx, c) in c_field.iter().enumerate() {
        let cm_c0 = matvec6(c, &mean_strain);
        let c0_e = matvec6(&c0, &mean_strain);
        for m in 0..6 {
            polar[idx][m] = cm_c0[m] - c0_e[m];
        }
    }

    // Reusable real-space strain field.
    let mut strain = vec![mean_strain; n_total];
    let zero_imag = vec![0.0_f64; n_total];

    let mut iterations = 0usize;
    for it in 0..max_iter {
        iterations = it + 1;

        // ---- Step 1: ε(x) = ε̄ − Γ⁰ * p(x) -------------------------------
        // FFT each polarization component.
        let mut polar_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
        for (idx, p) in polar.iter().enumerate() {
            for m in 0..6 {
                polar_re[m][idx] = p[m];
            }
        }
        let mut polar_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
        let mut polar_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
        for m in 0..6 {
            let (re, im) = oxifft::fft3d_split::<f64>(&polar_re[m], &zero_imag, nx, ny, nz);
            polar_hat_re[m] = re;
            polar_hat_im[m] = im;
        }

        // ε̂(ξ) = −Γ⁰(ξ):p̂(ξ)  for ξ≠0;  ε̂(0) = N·ε̄.
        let mut strain_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
        let mut strain_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
        for ix in 0..nx {
            let kx = wrapped_freq(ix, nx);
            for iy in 0..ny {
                let ky = wrapped_freq(iy, ny);
                for iz in 0..nz {
                    let flat = ix * ny * nz + iy * nz + iz;
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue;
                    }
                    let xi = [kx, ky, wrapped_freq(iz, nz)];
                    let tau_re: [f64; 6] = std::array::from_fn(|m| polar_hat_re[m][flat]);
                    let tau_im: [f64; 6] = std::array::from_fn(|m| polar_hat_im[m][flat]);
                    let e_re = apply_gamma0(xi, tau_re, lambda0, mu0);
                    let e_im = apply_gamma0(xi, tau_im, lambda0, mu0);
                    for m in 0..6 {
                        strain_hat_re[m][flat] = -e_re[m];
                        strain_hat_im[m][flat] = -e_im[m];
                    }
                }
            }
        }
        for m in 0..6 {
            strain_hat_re[m][0] = mean_strain[m] * n_f;
            strain_hat_im[m][0] = 0.0;
        }
        for m in 0..6 {
            let (re, _im) =
                oxifft::ifft3d_split::<f64>(&strain_hat_re[m], &strain_hat_im[m], nx, ny, nz);
            for (idx, value) in re.into_iter().enumerate() {
                strain[idx][m] = value;
            }
        }

        // ---- Step 2+3: relaxed polarization update ------------------------
        // q(x) = (C(x) − C⁰):ε(x);  p ← p + ω(q − p).
        for (idx, c) in c_field.iter().enumerate() {
            let c_e = matvec6(c, &strain[idx]);
            let c0_e = matvec6(&c0, &strain[idx]);
            for m in 0..6 {
                let q = c_e[m] - c0_e[m];
                polar[idx][m] += omega * (q - polar[idx][m]);
            }
        }

        // ---- Convergence: equilibrium residual of σ = C:ε ----------------
        let residual = equilibrium_residual(c_field, &strain, dims, &zero_imag);
        if residual <= tol {
            return Ok(LsOutcome {
                field: strain,
                iterations,
                residual,
                converged: true,
            });
        }
    }

    let residual = equilibrium_residual(c_field, &strain, dims, &zero_imag);
    Ok(LsOutcome {
        field: strain,
        iterations,
        residual,
        converged: false,
    })
}

/// Normalized equilibrium residual `‖ξ·σ̂(ξ)‖ / ‖σ̂(0)‖` for σ = C:ε.
fn equilibrium_residual(
    c_field: &[[[f64; 6]; 6]],
    strain: &[[f64; 6]],
    dims: [usize; 3],
    zero_imag: &[f64],
) -> f64 {
    let [nx, ny, nz] = dims;
    let n_total = nx * ny * nz;
    let mut sigma_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
    let mut sigma_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
    for m in 0..6 {
        let comp: Vec<f64> = (0..n_total)
            .map(|idx| matvec6(&c_field[idx], &strain[idx])[m])
            .collect();
        let (re, im) = oxifft::fft3d_split::<f64>(&comp, zero_imag, nx, ny, nz);
        sigma_hat_re[m] = re;
        sigma_hat_im[m] = im;
    }

    let full = |hat: &[Vec<f64>; 6], flat: usize| -> [[f64; 3]; 3] {
        let mut t = [[0.0_f64; 3]; 3];
        for (v, row) in hat.iter().enumerate() {
            let (i, j) = voigt_to_pair(v);
            t[i][j] = row[flat];
            t[j][i] = row[flat];
        }
        t
    };

    let mut residual_sq = 0.0_f64;
    let mut mean_norm_sq = 0.0_f64;
    for ix in 0..nx {
        let kx = wrapped_freq(ix, nx);
        for iy in 0..ny {
            let ky = wrapped_freq(iy, ny);
            for iz in 0..nz {
                let flat = ix * ny * nz + iy * nz + iz;
                if ix == 0 && iy == 0 && iz == 0 {
                    for s_re in sigma_hat_re.iter() {
                        mean_norm_sq += s_re[flat] * s_re[flat];
                    }
                    continue;
                }
                let xi = [kx, ky, wrapped_freq(iz, nz)];
                let sr = full(&sigma_hat_re, flat);
                let si = full(&sigma_hat_im, flat);
                for i in 0..3 {
                    let dr = sr[i][0] * xi[0] + sr[i][1] * xi[1] + sr[i][2] * xi[2];
                    let di = si[i][0] * xi[0] + si[i][1] * xi[1] + si[i][2] * xi[2];
                    residual_sq += dr * dr + di * di;
                }
            }
        }
    }
    let denom = mean_norm_sq.sqrt();
    if denom > 0.0 {
        residual_sq.sqrt() / denom
    } else {
        residual_sq.sqrt()
    }
}
