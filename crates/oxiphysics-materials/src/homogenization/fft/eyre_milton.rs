// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polarization-based ADMM solver for the periodic Lippmann-Schwinger equation.
//!
//! Implements the Brisard-Dormieux (2010) polarization scheme, which alternates
//! between a Fourier strain-update step (applying the reference Green operator
//! Γ⁰) and a local per-voxel polarization-update step.
//!
//! Given reference medium (λ₀, μ₀) the iteration is:
//!
//! ```text
//! τ₀ = 0
//! ε^{n+1}(ξ) = ε̄ − Γ⁰(ξ) τ̂^n(ξ)          (Fourier update; DC = N·ε̄)
//! σ^{n+1}    = C : ε^{n+1}                   (local stress)
//! rhs        = σ^{n+1} − C₀ : ε^{n+1} − τⁿ (residual polarization)
//! ψ          = (C + C₀)⁻¹ : rhs              (local solve per voxel)
//! τ^{n+1}    = τⁿ + 2·C₀ : ψ               (polarization update, factor 2)
//! ```
//!
//! With a geometric-mean reference `C₀ = (√(λ_min·λ_max), √(μ_min·μ_max))` and
//! phase contrast κ, the spectral radius of the error propagator is
//! approximately `|c₀/(c+c₀)|` for each mode/phase, giving O(√κ) convergence.

use super::green_operator::apply_gamma0;
use super::lippmann_schwinger::LsOutcome;
use crate::homogenization::HomogenizationError;

// ─────────────────────────────────── helpers ──────────────────────────────────

#[inline]
fn wrapped_freq(i: usize, n: usize) -> f64 {
    if i <= n / 2 {
        i as f64
    } else {
        i as f64 - n as f64
    }
}

#[inline]
fn matvec6(c: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i] += c[i][j] * v[j];
        }
    }
    out
}

/// Apply the isotropic reference stiffness C₀ : ε.
///
/// For Voigt engineering strain `eps` (shear entries γ = 2ε):
///   σ_nn = λ₀·tr(ε) + 2μ₀·ε_nn
///   σ_ij = μ₀·γ_ij   (i≠j)
fn apply_c0(lambda0: f64, mu0: f64, eps: &[f64; 6]) -> [f64; 6] {
    let trace = eps[0] + eps[1] + eps[2];
    [
        lambda0 * trace + 2.0 * mu0 * eps[0],
        lambda0 * trace + 2.0 * mu0 * eps[1],
        lambda0 * trace + 2.0 * mu0 * eps[2],
        mu0 * eps[3],
        mu0 * eps[4],
        mu0 * eps[5],
    ]
}

/// Solve `(C_iso) : ψ = v` for an isotropic stiffness with Lamé moduli
/// `(lam_eff, mu_eff)`.
///
/// Input `v` is a Voigt stress (shear entries are plain tensor components).
/// Output is a Voigt engineering strain (shear entries γ = 2ε).
#[inline]
fn isotropic_solve(lam_eff: f64, mu_eff: f64, v: [f64; 6]) -> [f64; 6] {
    let tr_v = v[0] + v[1] + v[2];
    let denom = 2.0 * mu_eff * (3.0 * lam_eff + 2.0 * mu_eff);
    let lam_factor = if denom.abs() > 1e-300 { lam_eff / denom } else { 0.0 };
    [
        v[0] / (2.0 * mu_eff) - lam_factor * tr_v,
        v[1] / (2.0 * mu_eff) - lam_factor * tr_v,
        v[2] / (2.0 * mu_eff) - lam_factor * tr_v,
        v[3] / mu_eff,
        v[4] / mu_eff,
        v[5] / mu_eff,
    ]
}

/// General 6×6 linear solve via Gauss elimination with partial pivoting.
///
/// Solves `M·x = b`, returning `Some(x)` or `None` if the system is singular.
fn solve_6x6(m: &[[f64; 6]; 6], b: &[f64; 6]) -> Option<[f64; 6]> {
    const N: usize = 6;
    let mut aug = [[0.0_f64; 7]; N];
    for i in 0..N {
        for j in 0..N {
            aug[i][j] = m[i][j];
        }
        aug[i][N] = b[i];
    }
    for col in 0..N {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (offset, aug_row) in aug[(col + 1)..N].iter().enumerate() {
            let row = col + 1 + offset;
            let v = aug_row[col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-300 {
            return None;
        }
        aug.swap(col, max_row);
        let inv_pivot = 1.0 / aug[col][col];
        for row in (col + 1)..N {
            let factor = aug[row][col] * inv_pivot;
            let col_row: [f64; 7] = aug[col];
            for (k_off, &col_val) in col_row[col..=N].iter().enumerate() {
                let k = col + k_off;
                aug[row][k] -= factor * col_val;
            }
        }
    }
    let mut x = [0.0_f64; N];
    for i in (0..N).rev() {
        x[i] = aug[i][N];
        for j in (i + 1)..N {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    Some(x)
}

/// Solve `(C(x) + C₀) : ψ = rhs` at a single voxel.
///
/// Uses the fast isotropic path when `C(x)` is isotropic (detected from its
/// structure), falling back to Gauss elimination for anisotropic phases.
fn local_solve(c: &[[f64; 6]; 6], lambda0: f64, mu0: f64, rhs: [f64; 6]) -> [f64; 6] {
    // Fast path: extract Lamé params assuming isotropic structure.
    // For isotropic C: C[5][5] = mu_x, C[0][0] = lambda_x + 2*mu_x, C[0][1] = lambda_x.
    let mu_x = c[5][5];
    let lam_x = c[0][0] - 2.0 * mu_x;
    let iso_ok = (c[0][1] - lam_x).abs() < 1e-6 * (lam_x.abs() + 1.0)
        && (c[1][2] - lam_x).abs() < 1e-6 * (lam_x.abs() + 1.0)
        && c[0][3].abs() < 1e-6 * (mu_x.abs() + 1.0)
        && (c[3][3] - mu_x).abs() < 1e-6 * (mu_x.abs() + 1.0)
        && (c[4][4] - mu_x).abs() < 1e-6 * (mu_x.abs() + 1.0);
    if iso_ok {
        return isotropic_solve(lam_x + lambda0, mu_x + mu0, rhs);
    }
    // General path: build M = C + C0 and solve via Gauss elimination.
    let mut m = *c;
    let lam2mu0 = lambda0 + 2.0 * mu0;
    m[0][0] += lam2mu0;
    m[1][1] += lam2mu0;
    m[2][2] += lam2mu0;
    m[0][1] += lambda0;
    m[1][0] += lambda0;
    m[0][2] += lambda0;
    m[2][0] += lambda0;
    m[1][2] += lambda0;
    m[2][1] += lambda0;
    m[3][3] += mu0;
    m[4][4] += mu0;
    m[5][5] += mu0;
    solve_6x6(&m, &rhs).unwrap_or([0.0; 6])
}

// ─────────────────────────────── public API ───────────────────────────────────

/// Solve the periodic Lippmann-Schwinger problem using the Brisard-Dormieux
/// polarization ADMM scheme with the supplied reference medium.
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

/// Core ADMM polarization iteration.
///
/// Alternates between:
/// 1. Fourier step: apply Γ⁰ to update strain from current polarization.
/// 2. Local step: update polarization τ via the per-voxel proximal operator.
///
/// The iteration terminates when the relative polarization update
/// ‖Δτ‖/‖σ₀‖ falls below `tol`, where ‖σ₀‖ = N·‖C₀ : ε̄‖₂ is a fixed
/// reference scale computed once before the loop.
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
                "grid {nx}x{ny}x{nz} has zero or overflowing voxel count"
            )));
        }
    };
    if c_field.len() != n_total {
        return Err(HomogenizationError::InvalidGrid(format!(
            "c_field length {} does not match grid {nx}x{ny}x{nz} = {n_total}",
            c_field.len()
        )));
    }

    let (lambda0, mu0) = ref_moduli;
    let zero_imag = vec![0.0_f64; n_total];
    let n_f64 = n_total as f64;

    // Fixed reference scale = N * ||C₀ : ε̄||₂
    let c0_eps_mean = apply_c0(lambda0, mu0, &mean_strain);
    let sigma0_scale = {
        let sq: f64 = c0_eps_mean.iter().map(|v| v * v).sum();
        n_f64 * sq.sqrt().max(1e-300)
    };

    // ── Scratch buffers ───────────────────────────────────────────────────────
    // Real-space polarization τ(x), initialised to zero.
    let mut tau: Vec<[f64; 6]> = vec![[0.0_f64; 6]; n_total];
    // Fourier τ̂ components (initialised to zero = transform of zero field).
    let mut tau_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
    let mut tau_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);

    // Real-space strain ε(x), initialised to uniform mean strain.
    let mut eps: Vec<[f64; 6]> = vec![mean_strain; n_total];

    // Fourier strain ε̂ — recomputed each iteration.
    let mut eps_hat_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);
    let mut eps_hat_im: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);

    // Scratch for FFT of each τ component.
    let mut tau_re_comp: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);

    // Real-space sigma σ(x) = C : ε(x)
    let mut sigma_re: [Vec<f64>; 6] = std::array::from_fn(|_| vec![0.0; n_total]);

    let mut iterations = 0usize;
    let mut tau_update_sq = 0.0_f64;
    let mut best_residual = f64::INFINITY;

    for it in 0..max_iter {
        iterations = it + 1;

        // ── Step 1: Fourier strain update ε̂ = ε̄ − Γ⁰ · τ̂ ──────────────────
        // FFT τ → τ̂
        for m in 0..6 {
            for (idx, t) in tau.iter().enumerate() {
                tau_re_comp[m][idx] = t[m];
            }
            let (re, im) = oxifft::fft3d_split::<f64>(&tau_re_comp[m], &zero_imag, nx, ny, nz);
            tau_hat_re[m] = re;
            tau_hat_im[m] = im;
        }

        // DC: ε̂(0) = N · ε̄
        for m in 0..6 {
            eps_hat_re[m][0] = n_f64 * mean_strain[m];
            eps_hat_im[m][0] = 0.0;
        }

        // Non-DC: ε̂(ξ) = −Γ⁰(ξ) · τ̂(ξ)
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue;
                    }
                    let flat = ix * ny * nz + iy * nz + iz;
                    let xi = [
                        wrapped_freq(ix, nx),
                        wrapped_freq(iy, ny),
                        wrapped_freq(iz, nz),
                    ];
                    let tr: [f64; 6] = std::array::from_fn(|m| tau_hat_re[m][flat]);
                    let ti: [f64; 6] = std::array::from_fn(|m| tau_hat_im[m][flat]);
                    let gr = apply_gamma0(xi, tr, lambda0, mu0);
                    let gi = apply_gamma0(xi, ti, lambda0, mu0);
                    for m in 0..6 {
                        eps_hat_re[m][flat] = -gr[m];
                        eps_hat_im[m][flat] = -gi[m];
                    }
                }
            }
        }

        // IFFT ε̂ → ε
        for m in 0..6 {
            let (re, _) =
                oxifft::ifft3d_split::<f64>(&eps_hat_re[m], &eps_hat_im[m], nx, ny, nz);
            for (idx, val) in re.into_iter().enumerate() {
                eps[idx][m] = val;
            }
        }

        // ── Step 2: Compute σ = C : ε in real space ──────────────────────────
        for (idx, (c, e)) in c_field.iter().zip(eps.iter()).enumerate() {
            let s = matvec6(c, e);
            for m in 0..6 {
                sigma_re[m][idx] = s[m];
            }
        }

        // ── Step 3: Local polarization update τ^{n+1} = τⁿ + 2·C₀ : ψ ────────
        tau_update_sq = 0.0_f64;
        for (idx, (c, e)) in c_field.iter().zip(eps.iter()).enumerate() {
            let sigma_x: [f64; 6] = std::array::from_fn(|m| sigma_re[m][idx]);
            let c0_eps = apply_c0(lambda0, mu0, e);
            let mut rhs = [0.0_f64; 6];
            for m in 0..6 {
                rhs[m] = sigma_x[m] - c0_eps[m] - tau[idx][m];
            }
            let psi = local_solve(c, lambda0, mu0, rhs);
            let c0_psi = apply_c0(lambda0, mu0, &psi);
            for m in 0..6 {
                let delta = 2.0 * c0_psi[m];
                tau_update_sq += delta * delta;
                tau[idx][m] += delta;
            }
        }

        // ── Step 4: Check convergence via polarization-update criterion ─────────
        let residual = tau_update_sq.sqrt() / sigma0_scale;
        // Divergence guard: break early on NaN or catastrophic growth.
        if !residual.is_finite() || residual > best_residual * 1.0e6 {
            break;
        }
        if residual < best_residual {
            best_residual = residual;
        }
        if residual <= tol {
            return Ok(LsOutcome {
                field: eps,
                iterations,
                residual,
                converged: true,
            });
        }
    }

    // Not converged — return the final field with the last polarization-update residual.
    let final_residual = tau_update_sq.sqrt() / sigma0_scale;
    Ok(LsOutcome {
        field: eps,
        iterations,
        residual: final_residual,
        converged: false,
    })
}
