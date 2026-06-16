// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Basic Moulinec–Suquet fixed-point scheme for the periodic
//! Lippmann–Schwinger equation.
//!
//! Given a voxel-wise stiffness field C(x) (6×6 Voigt per voxel) and a
//! prescribed macroscopic strain ε̄, the local strain field ε(x) satisfies the
//! periodic Lippmann–Schwinger equation
//!
//! ```text
//! ε(x) = ε̄ − Γ⁰ * (C(x):ε(x) − C⁰:ε(x))      [polarization form]
//! ```
//!
//! The classical basic scheme iterates directly on the strain field:
//!
//! 1. σ(x) = C(x):ε(x)
//! 2. σ̂(ξ) = FFT[σ](ξ)
//! 3. ε̂(ξ) = ε̂(ξ) − Γ⁰(ξ):σ̂(ξ)        for ξ ≠ 0
//! 4. ε̂(0) = ε̄
//! 5. ε(x) = IFFT[ε̂](x)
//!
//! repeated until the equilibrium residual `‖div σ‖` (evaluated in Fourier
//! space as `‖ξ·σ̂(ξ)‖`) drops below the tolerance.
//!
//! The transforms come from `oxifft::fft3d_split` / `oxifft::ifft3d_split`,
//! which use the convention: forward FFT is the unnormalized DFT sum and the
//! inverse FFT divides by `N = nx·ny·nz`, so `ifft(fft(x)) = x`.

use super::green_operator::{apply_gamma0, voigt_to_pair};
use crate::homogenization::HomogenizationError;

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

/// Wrapped integer frequency for FFT index `i` on a grid of size `n`.
///
/// Maps `i ∈ [0, n)` to `i` for the lower half and `i − n` for the upper half,
/// matching the standard FFT layout (and `pme.rs`).
#[inline]
fn wrapped_freq(i: usize, n: usize) -> f64 {
    if i <= n / 2 {
        i as f64
    } else {
        i as f64 - n as f64
    }
}

/// Outcome of an internal solver sweep, exposing the field regardless of
/// whether the tolerance was met.
pub(crate) struct LsOutcome {
    /// Converged-or-final strain field (one `[f64;6]` per voxel).
    pub field: Vec<[f64; 6]>,
    /// Iterations performed.
    pub iterations: usize,
    /// Final normalized equilibrium residual.
    pub residual: f64,
    /// Whether `residual <= tol`.
    pub converged: bool,
}

/// Solve the periodic Lippmann–Schwinger problem by the basic Moulinec–Suquet
/// fixed-point scheme.
///
/// # Arguments
/// * `c_field`     — flat voxel stiffness field, length `nx*ny*nz`, row-major
///   index `ix*ny*nz + iy*nz + iz`. Each entry is a 6×6 Voigt stiffness.
/// * `dims`        — grid dimensions `[nx, ny, nz]`.
/// * `mean_strain` — prescribed macroscopic strain ε̄ (engineering-strain Voigt).
/// * `ref_moduli`  — reference medium `(λ₀, μ₀)`.
/// * `tol`         — relative equilibrium-residual tolerance.
/// * `max_iter`    — maximum number of fixed-point iterations.
///
/// Returns the converged strain field (one `[f64;6]` per voxel) and the number
/// of iterations performed, or [`HomogenizationError::NotConverged`] if the
/// residual never reached `tol`.
pub fn lippmann_schwinger(
    c_field: &[[[f64; 6]; 6]],
    dims: [usize; 3],
    mean_strain: [f64; 6],
    ref_moduli: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<[f64; 6]>, usize), HomogenizationError> {
    let outcome = lippmann_schwinger_inner(c_field, dims, mean_strain, ref_moduli, tol, max_iter)?;
    if outcome.converged {
        Ok((outcome.field, outcome.iterations))
    } else {
        Err(HomogenizationError::NotConverged {
            max_iter,
            residual: outcome.residual,
        })
    }
}

/// Core sweep that always returns the final field together with the residual
/// and convergence flag. Used directly by the effective-stiffness assembler so
/// that a (possibly unconverged) field can still be averaged.
pub(crate) fn lippmann_schwinger_inner(
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

    let n_f = n_total as f64;

    // Initialize strain_hat in Fourier space: DC = mean_strain * N, all others 0.
    // This avoids re-FFTing the strain field every iteration (saves 6 FFTs/iter).
    let mut strain_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
    let mut strain_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
    for m in 0..6 {
        strain_hat_re[m][0] = mean_strain[m] * n_f;
    }

    // Initial real-space strain field from the initial strain_hat (uniform mean strain).
    let mut strain = vec![[0.0_f64; 6]; n_total];
    for m in 0..6 {
        let (re, _im) =
            oxifft::ifft3d_split::<f64>(&strain_hat_re[m], &strain_hat_im[m], nx, ny, nz);
        for (idx, value) in re.into_iter().enumerate() {
            strain[idx][m] = value;
        }
    }

    // Reusable split-complex buffers: 6 Voigt components of the stress field.
    let mut sigma_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
    let zero_imag = vec![0.0_f64; n_total];

    let mut iterations = 0usize;
    let mut best_residual = f64::INFINITY;
    for it in 0..max_iter {
        iterations = it + 1;

        // 1. Local stress σ(x) = C(x):ε(x), stored component-wise.
        for (idx, (c, eps)) in c_field.iter().zip(strain.iter()).enumerate() {
            let s = matvec6(c, eps);
            for m in 0..6 {
                sigma_re[m][idx] = s[m];
            }
        }

        // 2. Forward FFT of each stress component (6 FFTs).
        let mut sigma_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
        let mut sigma_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
        for m in 0..6 {
            let (re, im) = oxifft::fft3d_split::<f64>(&sigma_re[m], &zero_imag, nx, ny, nz);
            sigma_hat_re[m] = re;
            sigma_hat_im[m] = im;
        }

        // 3. Apply Γ⁰ update in Fourier space and accumulate equilibrium residual.
        //    ε̂(ξ) ← ε̂(ξ) − Γ⁰(ξ):σ̂(ξ)   for ξ ≠ 0  (increment form).
        //    strain_hat is maintained across iterations — no re-FFT of strain needed.
        let mut residual_sq = 0.0_f64;
        let mut mean_stress_norm_sq = 0.0_f64;

        for ix in 0..nx {
            let kx = wrapped_freq(ix, nx);
            for iy in 0..ny {
                let ky = wrapped_freq(iy, ny);
                for iz in 0..nz {
                    let flat = ix * ny * nz + iy * nz + iz;
                    let kz = wrapped_freq(iz, nz);

                    if ix == 0 && iy == 0 && iz == 0 {
                        // DC: mean stress, used to normalize the residual.
                        for s_re in sigma_hat_re.iter() {
                            mean_stress_norm_sq += s_re[flat] * s_re[flat];
                        }
                        continue;
                    }

                    // Equilibrium residual: r_i = Σ_j ξ_j σ̂_ij(ξ).
                    // Build the full symmetric stress tensor from the 6 Voigt
                    // components (real and imaginary parts handled separately).
                    let xi = [kx, ky, kz];
                    accumulate_residual(&xi, &sigma_hat_re, &sigma_hat_im, flat, &mut residual_sq);

                    // Apply Γ⁰ to the real and imaginary stress separately
                    // (Γ⁰ is real, so it acts independently on each part).
                    let tau_re: [f64; 6] = std::array::from_fn(|m| sigma_hat_re[m][flat]);
                    let tau_im: [f64; 6] = std::array::from_fn(|m| sigma_hat_im[m][flat]);
                    let eps_re = apply_gamma0(xi, tau_re, lambda0, mu0);
                    let eps_im = apply_gamma0(xi, tau_im, lambda0, mu0);
                    // Increment form: ε̂ ← ε̂ − Γ⁰:σ̂  (NOT replacement).
                    for m in 0..6 {
                        strain_hat_re[m][flat] -= eps_re[m];
                        strain_hat_im[m][flat] -= eps_im[m];
                    }
                }
            }
        }

        // 4. Enforce the prescribed mean strain at the DC component.
        //    The inverse FFT divides by N, so the DC coefficient must be
        //    N·ε̄ to recover a real-space mean of ε̄.
        for m in 0..6 {
            strain_hat_re[m][0] = mean_strain[m] * n_f;
            strain_hat_im[m][0] = 0.0;
        }

        // 5. Inverse FFT back to the real-space strain field (6 IFFTs).
        for m in 0..6 {
            let (re, _im) =
                oxifft::ifft3d_split::<f64>(&strain_hat_re[m], &strain_hat_im[m], nx, ny, nz);
            for (idx, value) in re.into_iter().enumerate() {
                strain[idx][m] = value;
            }
        }

        // 7. Convergence test on the normalized equilibrium residual.
        let denom = mean_stress_norm_sq.sqrt();
        let residual = if denom > 0.0 {
            residual_sq.sqrt() / denom
        } else {
            residual_sq.sqrt()
        };
        // Divergence guard: break early if residual is non-finite or has grown
        // catastrophically (spectral radius >= 1 → geometric blow-up → NaN).
        if !residual.is_finite() || residual > best_residual * 1.0e6 {
            break;
        }
        if residual < best_residual {
            best_residual = residual;
        }
        if residual <= tol {
            return Ok(LsOutcome {
                field: strain,
                iterations,
                residual,
                converged: true,
            });
        }
    }

    // Did not converge: recompute the final residual for the report.
    let residual = final_residual(c_field, &strain, nx, ny, nz);
    Ok(LsOutcome {
        field: strain,
        iterations,
        residual,
        converged: false,
    })
}

/// Accumulate the squared equilibrium residual `Σ_i |Σ_j ξ_j σ̂_ij|²` at one
/// frequency into `residual_sq`. The stress tensor is reconstructed from its 6
/// Voigt components; real and imaginary contributions add in quadrature.
fn accumulate_residual(
    xi: &[f64; 3],
    sigma_hat_re: &[Vec<f64>; 6],
    sigma_hat_im: &[Vec<f64>; 6],
    flat: usize,
    residual_sq: &mut f64,
) {
    // Full 3×3 stress from Voigt (real and imag).
    let full = |hat: &[Vec<f64>; 6]| -> [[f64; 3]; 3] {
        let mut t = [[0.0_f64; 3]; 3];
        for (v, row) in hat.iter().enumerate() {
            let (i, j) = voigt_to_pair(v);
            let val = row[flat];
            t[i][j] = val;
            t[j][i] = val;
        }
        t
    };
    let sr = full(sigma_hat_re);
    let si = full(sigma_hat_im);
    for i in 0..3 {
        let div_re = sr[i][0] * xi[0] + sr[i][1] * xi[1] + sr[i][2] * xi[2];
        let div_im = si[i][0] * xi[0] + si[i][1] * xi[1] + si[i][2] * xi[2];
        *residual_sq += div_re * div_re + div_im * div_im;
    }
}

/// Recompute the normalized equilibrium residual for a converged-or-not strain
/// field (used to report the final residual on non-convergence).
fn final_residual(
    c_field: &[[[f64; 6]; 6]],
    strain: &[[f64; 6]],
    nx: usize,
    ny: usize,
    nz: usize,
) -> f64 {
    let n_total = nx * ny * nz;
    let zero_imag = vec![0.0_f64; n_total];
    let mut sigma_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
    let mut sigma_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::new());
    for m in 0..6 {
        let comp: Vec<f64> = (0..n_total)
            .map(|idx| matvec6(&c_field[idx], &strain[idx])[m])
            .collect();
        let (re, im) = oxifft::fft3d_split::<f64>(&comp, &zero_imag, nx, ny, nz);
        sigma_hat_re[m] = re;
        sigma_hat_im[m] = im;
    }
    let mut residual_sq = 0.0_f64;
    let mut mean_stress_norm_sq = 0.0_f64;
    for ix in 0..nx {
        let kx = wrapped_freq(ix, nx);
        for iy in 0..ny {
            let ky = wrapped_freq(iy, ny);
            for iz in 0..nz {
                let flat = ix * ny * nz + iy * nz + iz;
                let kz = wrapped_freq(iz, nz);
                if ix == 0 && iy == 0 && iz == 0 {
                    for s_re in sigma_hat_re.iter() {
                        mean_stress_norm_sq += s_re[flat] * s_re[flat];
                    }
                    continue;
                }
                accumulate_residual(
                    &[kx, ky, kz],
                    &sigma_hat_re,
                    &sigma_hat_im,
                    flat,
                    &mut residual_sq,
                );
            }
        }
    }
    let denom = mean_stress_norm_sq.sqrt();
    if denom > 0.0 {
        residual_sq.sqrt() / denom
    } else {
        residual_sq.sqrt()
    }
}
