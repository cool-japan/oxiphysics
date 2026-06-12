// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LJ-PME (Lennard-Jones Particle Mesh Ewald) long-range dispersion.
//!
//! Implements the reciprocal-space r⁻⁶ dispersion sum following:
//! in 't Veld, Ismail & Bhatt, J. Chem. Phys. 127, 144711 (2007).
//!
//! # Algorithm
//!
//! The r⁻⁶ potential is split using the Ewald decomposition:
//! ```text
//! 1/r⁶ = g₆(αr)/r⁶  [real space]  +  (1 - g₆(αr))/r⁶  [reciprocal space]
//! ```
//! where `g₆(x) = exp(-x²)(1 + x² + x⁴/2)`.
//!
//! The "charges" spread onto the grid are `√C₆ᵢ` (geometric combining rule),
//! so that `C₆ᵢⱼ = (√C₆ᵢ)(√C₆ⱼ)`.
//!
//! ## Reciprocal-space Green's function
//! ```text
//! Ĝ_disp(k) = (π^(3/2) / (2α³V)) * (k²/(4α²)) * exp(-k²/(4α²)) / |B(k)|²
//! ```
//! Note: the B-spline correction divides here (not multiplies), because we
//! are deconvolving the gridding kernel in k-space.
//!
//! ## Self-energy correction
//! ```text
//! E_self = -(π^(3/2) * α³ / 6) * Σᵢ C₆ᵢ
//! ```
//!
//! # Units
//! - Distance: Å
//! - C₆:       kJ mol⁻¹ Å⁶  (GROMACS-compatible)
//! - Energy:   kJ mol⁻¹
//! - Forces:   kJ mol⁻¹ Å⁻¹

use std::f64::consts::PI;
use thiserror::Error;

// ---------------------------------------------------------------------------
// B-spline order (must match pme.rs)
// ---------------------------------------------------------------------------

const BSPLINE_ORDER: usize = 4;

// ---------------------------------------------------------------------------
// B-spline helpers (private; mirrors pme.rs implementations exactly)
// ---------------------------------------------------------------------------

#[inline]
fn bspline(p: usize, x: f64) -> f64 {
    if p == 1 {
        return if (0.0..1.0).contains(&x) { 1.0 } else { 0.0 };
    }
    let pf = p as f64;
    let term1 = if x > 0.0 {
        x / (pf - 1.0) * bspline(p - 1, x)
    } else {
        0.0
    };
    let term2 = if pf - x > 0.0 {
        (pf - x) / (pf - 1.0) * bspline(p - 1, x - 1.0)
    } else {
        0.0
    };
    term1 + term2
}

#[inline]
fn bspline_deriv(p: usize, x: f64) -> f64 {
    if p <= 1 {
        return 0.0;
    }
    bspline(p - 1, x) - bspline(p - 1, x - 1.0)
}

fn bspline_weights(
    u: f64,
    n_grid: usize,
    order: usize,
) -> ([f64; BSPLINE_ORDER], [usize; BSPLINE_ORDER]) {
    let scaled = u * (n_grid as f64);
    let i0 = scaled.floor() as i64;
    let mut weights = [0.0f64; BSPLINE_ORDER];
    let mut grid_pts = [0usize; BSPLINE_ORDER];
    for j in 0..order {
        let raw = i0 - (order as i64 - 1) + j as i64;
        let x = scaled - raw as f64;
        weights[j] = bspline(order, x);
        grid_pts[j] = raw.rem_euclid(n_grid as i64) as usize;
    }
    (weights, grid_pts)
}

fn bspline_deriv_weights(u: f64, n_grid: usize, order: usize) -> [f64; BSPLINE_ORDER] {
    let scaled = u * (n_grid as f64);
    let i0 = scaled.floor() as i64;
    let mut dweights = [0.0f64; BSPLINE_ORDER];
    let scale = n_grid as f64;
    for (j, dw) in dweights.iter_mut().enumerate().take(order) {
        let raw = i0 - (order as i64 - 1) + j as i64;
        let x = scaled - raw as f64;
        *dw = bspline_deriv(order, x) * scale;
    }
    dweights
}

fn bspline_correction_sq(m: i64, n: usize, p: usize) -> f64 {
    let n_f = n as f64;
    let two_pi_m_over_n = 2.0 * PI * (m as f64) / n_f;
    let mut den_re = 0.0f64;
    let mut den_im = 0.0f64;
    for l in 0..(p - 1) {
        let bval = bspline(p, (l + 1) as f64);
        let angle = two_pi_m_over_n * (l as f64);
        den_re += bval * angle.cos();
        den_im += bval * angle.sin();
    }
    let num_angle = two_pi_m_over_n * ((p - 1) as f64);
    let num_re = num_angle.cos();
    let num_im = num_angle.sin();
    let num_sq = num_re * num_re + num_im * num_im;
    let den_sq = den_re * den_re + den_im * den_im;
    if den_sq < 1e-30 {
        return 0.0;
    }
    num_sq / den_sq
}

#[inline]
fn min_image(dx: f64, l: f64) -> f64 {
    dx - l * (dx / l).round()
}

// ---------------------------------------------------------------------------
// Ewald r⁻⁶ splitting function  g₆(x) = exp(-x²)(1 + x² + x⁴/2)
// ---------------------------------------------------------------------------

/// Short-range complement function for r⁻⁶ Ewald splitting.
///
/// The real-space contribution for pair (i,j) is:
/// `-C₆ᵢⱼ * g₆(α·r) / r⁶`
///
/// As r→0, g₆→1 (real-space recovers full potential).
/// As r→∞, g₆→0 (real-space sum decays rapidly).
#[inline]
pub fn g6(x: f64) -> f64 {
    let x2 = x * x;
    (-x2).exp() * (1.0 + x2 + x2 * x2 / 2.0)
}

/// Derivative dg₆/dx = -x⁵ · exp(-x²).
#[inline]
pub fn dg6_dx(x: f64) -> f64 {
    let x2 = x * x;
    -(x2 * x2 * x) * (-x2).exp()
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Parameters for the LJ-PME dispersion calculation.
#[derive(Debug, Clone)]
pub struct LjPmeParams {
    /// Ewald splitting parameter α (Å⁻¹). Larger α → more work in reciprocal space.
    pub alpha: f64,
    /// PME grid dimensions [mx, my, mz].
    pub grid: [usize; 3],
    /// Real-space cutoff radius (Å).
    pub r_cut: f64,
}

impl LjPmeParams {
    /// Create new LJ-PME parameters.
    ///
    /// # Arguments
    /// * `alpha`  — Ewald splitting parameter (Å⁻¹); typical range 0.3–2.0
    /// * `grid`   — PME grid dimensions; should be even and have small prime factors
    /// * `r_cut`  — real-space cutoff radius (Å)
    pub fn new(alpha: f64, grid: [usize; 3], r_cut: f64) -> Self {
        Self { alpha, grid, r_cut }
    }
}

/// Errors from the LJ-PME calculation.
#[derive(Debug, Error)]
pub enum LjPmeError {
    /// Input length mismatch between positions and C6 arrays.
    #[error("positions length {pos} != c6 length {c6}")]
    LengthMismatch {
        /// Number of positions provided.
        pos: usize,
        /// Number of C₆ coefficients provided.
        c6: usize,
    },
    /// Invalid parameter value.
    #[error("invalid LJ-PME parameter: {0}")]
    InvalidParam(String),
}

// ---------------------------------------------------------------------------
// Spreading √C₆ onto the PME grid
// ---------------------------------------------------------------------------

/// Spread `√C₆ᵢ` values onto a 3D grid using B-spline order 4.
///
/// Grid layout: row-major with index `ix * ny * nz + iy * nz + iz`.
fn spread_c6_sqrt_bspline(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    nx: usize,
    ny: usize,
    nz: usize,
) -> Vec<f64> {
    let mut rho = vec![0.0f64; nx * ny * nz];

    for (pos, &c6i) in positions.iter().zip(c6.iter()) {
        let sqrt_c6i = c6i.abs().sqrt();
        let ux = pos[0].rem_euclid(box_len[0]) / box_len[0];
        let uy = pos[1].rem_euclid(box_len[1]) / box_len[1];
        let uz = pos[2].rem_euclid(box_len[2]) / box_len[2];

        let (wx, gx) = bspline_weights(ux, nx, BSPLINE_ORDER);
        let (wy, gy) = bspline_weights(uy, ny, BSPLINE_ORDER);
        let (wz, gz) = bspline_weights(uz, nz, BSPLINE_ORDER);

        for jx in 0..BSPLINE_ORDER {
            for jy in 0..BSPLINE_ORDER {
                let wxy = wx[jx] * wy[jy];
                for jz in 0..BSPLINE_ORDER {
                    let idx = gx[jx] * ny * nz + gy[jy] * nz + gz[jz];
                    rho[idx] += sqrt_c6i * wxy * wz[jz];
                }
            }
        }
    }
    rho
}

// ---------------------------------------------------------------------------
// Reciprocal-space energy and potential grid
// ---------------------------------------------------------------------------

/// Compute the LJ-PME reciprocal-space energy (kJ mol⁻¹) and the potential
/// grid for force gathering.
///
/// Returns `(energy, pot_re)` where `pot_re` is the real-space potential grid
/// after IFFT of `Ĝ_disp(k) · ρ̂_disp(k)`.
fn lj_pme_reciprocal_energy_and_potential(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    alpha: f64,
    grid: [usize; 3],
) -> (f64, Vec<f64>) {
    let [mx, my, mz] = grid;
    let volume = box_len[0] * box_len[1] * box_len[2];
    let n_total = mx * my * mz;

    // 1. Spread √C₆ᵢ onto grid
    let rho = spread_c6_sqrt_bspline(positions, c6, box_len, mx, my, mz);

    // 2. 3D FFT via OxiFFT (imaginary part of real input = 0)
    let rho_imag = vec![0.0f64; n_total];
    let (rho_re, rho_im) = oxifft::fft3d_split::<f64>(&rho, &rho_imag, mx, my, mz);

    // 3. Build Ĝ_disp(k) and compute energy; multiply to get φ̂_disp(k)
    //
    // Ĝ_disp(k) = (π^(3/2) / (2α³V)) * (k² / (4α²)) * exp(-k² / (4α²)) / |B(k)|²
    //
    // Note: division by |B(k)|² deconvolves the B-spline spreading kernel.
    let mut phi_re = vec![0.0f64; n_total];
    let mut phi_im = vec![0.0f64; n_total];

    let two_pi = 2.0 * PI;
    let alpha3 = alpha * alpha * alpha;
    // Prefactor = π^(3/2) / (2 α³ V), multiplied by 0.5 for the 1/2 in front of Σ|S(k)|²
    let g_prefactor = PI.powf(1.5) / (2.0 * alpha3 * volume);
    let four_alpha_sq = 4.0 * alpha * alpha;

    let mut energy = 0.0f64;

    for ix in 0..mx {
        let nx = if ix <= mx / 2 { ix as i64 } else { ix as i64 - mx as i64 };
        let kx = two_pi * nx as f64 / box_len[0];
        let bc_x = bspline_correction_sq(nx, mx, BSPLINE_ORDER);

        for iy in 0..my {
            let ny = if iy <= my / 2 { iy as i64 } else { iy as i64 - my as i64 };
            let ky = two_pi * ny as f64 / box_len[1];
            let bc_y = bspline_correction_sq(ny, my, BSPLINE_ORDER);

            for iz in 0..mz {
                if ix == 0 && iy == 0 && iz == 0 {
                    continue; // k=0 term vanishes (no monopole dispersion)
                }
                let nz = if iz <= mz / 2 { iz as i64 } else { iz as i64 - mz as i64 };
                let kz = two_pi * nz as f64 / box_len[2];
                let bc_z = bspline_correction_sq(nz, mz, BSPLINE_ORDER);

                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 < 1e-20 {
                    continue;
                }

                let bspline_corr = bc_x * bc_y * bc_z;
                if bspline_corr < 1e-30 {
                    continue;
                }

                // Dispersion Green's function
                let u = k2 / four_alpha_sq;
                let gk = g_prefactor * u * (-u).exp() / bspline_corr;

                let flat = ix * my * mz + iy * mz + iz;
                // ρ̂(k) = rho_re + i*rho_im  (output of FFT of real input)
                let r_re = rho_re[flat];
                let r_im = rho_im[flat];

                // |ρ̂(k)|² contributes to energy
                let rho_sq = r_re * r_re + r_im * r_im;
                energy += gk * rho_sq;

                // φ̂(k) = G(k) * ρ̂(k)  for potential → forces
                phi_re[flat] = gk * r_re;
                phi_im[flat] = gk * r_im;
            }
        }
    }

    // Factor of 1/2 for double-counting (each pair counted twice)
    energy *= 0.5;

    // 4. Inverse 3D FFT to get real-space potential grid φ(r)
    let (pot_re, _pot_im) = oxifft::ifft3d_split::<f64>(&phi_re, &phi_im, mx, my, mz);

    (energy, pot_re)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the LJ-PME reciprocal-space dispersion energy (kJ mol⁻¹).
///
/// This is the long-range correction to the Lennard-Jones r⁻⁶ interaction,
/// computed via Particle Mesh Ewald in reciprocal space.
///
/// # Arguments
/// * `positions`  — atom positions (Å), length N
/// * `c6`         — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N; must be ≥ 0
/// * `box_len`    — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `params`     — LJ-PME parameters (α, grid, cutoff)
///
/// # Errors
/// Returns `LjPmeError::LengthMismatch` if `positions.len() != c6.len()`,
/// or `LjPmeError::InvalidParam` if α or box dimensions are non-positive.
pub fn lj_pme_reciprocal_energy(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<f64, LjPmeError> {
    validate_inputs(positions, c6, box_len, params)?;
    let (energy, _pot) =
        lj_pme_reciprocal_energy_and_potential(positions, c6, box_len, params.alpha, params.grid);
    Ok(energy)
}

/// Compute the LJ-PME self-energy correction (kJ mol⁻¹).
///
/// Removes the spurious self-interaction of each atom with its own periodic
/// images introduced by the Ewald splitting:
/// ```text
/// E_self = -(π^(3/2) * α³ / 6) * Σᵢ C₆ᵢ
/// ```
///
/// # Arguments
/// * `c6`    — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N
/// * `alpha` — Ewald splitting parameter (Å⁻¹)
pub fn lj_pme_self_energy(c6: &[f64], alpha: f64) -> f64 {
    let alpha3 = alpha * alpha * alpha;
    let sum_c6: f64 = c6.iter().sum();
    -(PI.powf(1.5) * alpha3 / 6.0) * sum_c6
}

/// Compute the real-space LJ dispersion energy (kJ mol⁻¹) with Ewald damping.
///
/// ```text
/// E_real = -Σ_{i<j}  C₆ᵢⱼ · g₆(α·rᵢⱼ) / rᵢⱼ⁶
/// ```
/// where `C₆ᵢⱼ = √(C₆ᵢ · C₆ⱼ)` (geometric combining rule) and
/// `g₆(x) = exp(-x²)(1 + x² + x⁴/2)`.
///
/// # Arguments
/// * `positions`  — atom positions (Å), length N
/// * `c6`         — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N
/// * `box_len`    — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `params`     — LJ-PME parameters (uses `alpha` and `r_cut`)
///
/// # Errors
/// Returns `LjPmeError::LengthMismatch` if `positions.len() != c6.len()`.
pub fn lj_pme_real_space_energy(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<f64, LjPmeError> {
    validate_inputs(positions, c6, box_len, params)?;
    let alpha = params.alpha;
    let r_cut = params.r_cut;
    let r_cut_sq = r_cut * r_cut;
    let n = positions.len();
    let mut energy = 0.0f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = min_image(positions[j][0] - positions[i][0], box_len[0]);
            let dy = min_image(positions[j][1] - positions[i][1], box_len[1]);
            let dz = min_image(positions[j][2] - positions[i][2], box_len[2]);
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > r_cut_sq || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            let c6ij = (c6[i] * c6[j]).sqrt();
            let r6 = r2 * r2 * r2;
            energy -= c6ij * g6(alpha * r) / r6;
        }
    }
    Ok(energy)
}

/// Compute the LJ-PME reciprocal-space forces (kJ mol⁻¹ Å⁻¹).
///
/// Gathers forces from the potential grid via B-spline interpolation of
/// the potential gradient.
///
/// # Arguments
/// * `positions`  — atom positions (Å), length N
/// * `c6`         — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N
/// * `box_len`    — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `params`     — LJ-PME parameters
///
/// # Errors
/// Returns `LjPmeError` on invalid inputs.
pub fn lj_pme_reciprocal_forces(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<Vec<[f64; 3]>, LjPmeError> {
    validate_inputs(positions, c6, box_len, params)?;
    let [mx, my, mz] = params.grid;
    let n_atoms = positions.len();
    let n_total = mx * my * mz;

    let (_energy, pot_re) =
        lj_pme_reciprocal_energy_and_potential(positions, c6, box_len, params.alpha, params.grid);

    // Gather forces via B-spline gradient of the potential
    // F_i^α = -√C₆ᵢ * Σ_m (dM/du_α) * M_y * M_z * φ(m) * (n_total / L_α)
    // (The n_total factor cancels the IFFT normalization 1/N.)
    let ntot_f64 = n_total as f64;
    let mut forces = vec![[0.0f64; 3]; n_atoms];

    for (atom_idx, (pos, &c6i)) in positions.iter().zip(c6.iter()).enumerate() {
        let sqrt_c6i = c6i.abs().sqrt();
        let ux = pos[0].rem_euclid(box_len[0]) / box_len[0];
        let uy = pos[1].rem_euclid(box_len[1]) / box_len[1];
        let uz = pos[2].rem_euclid(box_len[2]) / box_len[2];

        let (wx, gx) = bspline_weights(ux, mx, BSPLINE_ORDER);
        let (wy, gy) = bspline_weights(uy, my, BSPLINE_ORDER);
        let (wz, gz) = bspline_weights(uz, mz, BSPLINE_ORDER);

        let dwx = bspline_deriv_weights(ux, mx, BSPLINE_ORDER);
        let dwy = bspline_deriv_weights(uy, my, BSPLINE_ORDER);
        let dwz = bspline_deriv_weights(uz, mz, BSPLINE_ORDER);

        let mut fx = 0.0f64;
        let mut fy = 0.0f64;
        let mut fz = 0.0f64;

        for jx in 0..BSPLINE_ORDER {
            for jy in 0..BSPLINE_ORDER {
                for jz in 0..BSPLINE_ORDER {
                    let flat = gx[jx] * my * mz + gy[jy] * mz + gz[jz];
                    let phi = pot_re[flat];
                    fx += dwx[jx] * wy[jy] * wz[jz] * phi;
                    fy += wx[jx] * dwy[jy] * wz[jz] * phi;
                    fz += wx[jx] * wy[jy] * dwz[jz] * phi;
                }
            }
        }

        forces[atom_idx][0] = -sqrt_c6i * ntot_f64 * fx / box_len[0];
        forces[atom_idx][1] = -sqrt_c6i * ntot_f64 * fy / box_len[1];
        forces[atom_idx][2] = -sqrt_c6i * ntot_f64 * fz / box_len[2];
    }

    Ok(forces)
}

// ---------------------------------------------------------------------------
// Unified energy + forces entry point
// ---------------------------------------------------------------------------

/// Compute the total LJ-PME dispersion energy and forces.
///
/// Combines real-space, reciprocal-space, and self-energy contributions:
/// ```text
/// E_total = E_real + E_recip + E_self
/// ```
///
/// # Arguments
/// * `positions` — atom positions (Å), length N
/// * `c6`        — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N; must be ≥ 0
/// * `box_len`   — orthorhombic box side lengths [Lx, Ly, Lz] (Å)
/// * `params`    — LJ-PME parameters (α, grid, r_cut)
///
/// # Returns
/// `(energy, forces)` where energy is in kJ mol⁻¹ and forces in kJ mol⁻¹ Å⁻¹.
///
/// # Errors
/// Returns [`LjPmeError`] if input lengths are inconsistent or parameters are invalid.
pub fn lj_pme_energy_and_forces(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<(f64, Vec<[f64; 3]>), LjPmeError> {
    validate_inputs(positions, c6, box_len, params)?;
    let [mx, my, mz] = params.grid;
    let n_atoms = positions.len();
    let alpha = params.alpha;
    let r_cut = params.r_cut;

    // Real-space energy and forces
    let r_cut_sq = r_cut * r_cut;
    let n = positions.len();
    let mut e_real = 0.0f64;
    let mut forces = vec![[0.0f64; 3]; n_atoms];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = min_image(positions[j][0] - positions[i][0], box_len[0]);
            let dy = min_image(positions[j][1] - positions[i][1], box_len[1]);
            let dz = min_image(positions[j][2] - positions[i][2], box_len[2]);
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > r_cut_sq || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            let r6 = r2 * r2 * r2;
            let c6ij = (c6[i] * c6[j]).abs().sqrt();
            let s = alpha * r;
            let g = g6(s);
            e_real -= c6ij * g / r6;

            // Force: F_i = -(dU/dr) * (r_j - r_i)/r
            // dU/dr = -c6ij * (dg6_dx(s)*alpha/r^6 - 6*g/r^7)
            //       = -c6ij * (dg6_dx(s)*alpha - 6*g/r) / r^6
            let dg = dg6_dx(s);
            let du_dr = -c6ij * (dg * alpha - 6.0 * g / r) / r6;
            let f_scalar = -du_dr / r;
            forces[i][0] += f_scalar * dx;
            forces[i][1] += f_scalar * dy;
            forces[i][2] += f_scalar * dz;
            forces[j][0] -= f_scalar * dx;
            forces[j][1] -= f_scalar * dy;
            forces[j][2] -= f_scalar * dz;
        }
    }

    // Reciprocal-space energy + potential grid for forces
    let (e_recip, pot_re) =
        lj_pme_reciprocal_energy_and_potential(positions, c6, box_len, alpha, params.grid);

    // Self-energy
    let e_self = lj_pme_self_energy(c6, alpha);

    // Gather reciprocal forces
    let n_total = mx * my * mz;
    let ntot_f64 = n_total as f64;

    for (atom_idx, (pos, &c6i)) in positions.iter().zip(c6.iter()).enumerate() {
        let sqrt_c6i = c6i.abs().sqrt();
        let ux = pos[0].rem_euclid(box_len[0]) / box_len[0];
        let uy = pos[1].rem_euclid(box_len[1]) / box_len[1];
        let uz = pos[2].rem_euclid(box_len[2]) / box_len[2];

        let (wx, gx) = bspline_weights(ux, mx, BSPLINE_ORDER);
        let (wy, gy) = bspline_weights(uy, my, BSPLINE_ORDER);
        let (wz, gz) = bspline_weights(uz, mz, BSPLINE_ORDER);

        let dwx = bspline_deriv_weights(ux, mx, BSPLINE_ORDER);
        let dwy = bspline_deriv_weights(uy, my, BSPLINE_ORDER);
        let dwz = bspline_deriv_weights(uz, mz, BSPLINE_ORDER);

        let mut fx = 0.0f64;
        let mut fy = 0.0f64;
        let mut fz = 0.0f64;

        for jx in 0..BSPLINE_ORDER {
            for jy in 0..BSPLINE_ORDER {
                for jz in 0..BSPLINE_ORDER {
                    let flat = gx[jx] * my * mz + gy[jy] * mz + gz[jz];
                    let phi = pot_re[flat];
                    fx += dwx[jx] * wy[jy] * wz[jz] * phi;
                    fy += wx[jx] * dwy[jy] * wz[jz] * phi;
                    fz += wx[jx] * wy[jy] * dwz[jz] * phi;
                }
            }
        }

        forces[atom_idx][0] += -sqrt_c6i * ntot_f64 * fx / box_len[0];
        forces[atom_idx][1] += -sqrt_c6i * ntot_f64 * fy / box_len[1];
        forces[atom_idx][2] += -sqrt_c6i * ntot_f64 * fz / box_len[2];
    }

    let energy = e_real + e_recip + e_self;
    Ok((energy, forces))
}

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------

fn validate_inputs(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<(), LjPmeError> {
    if positions.len() != c6.len() {
        return Err(LjPmeError::LengthMismatch {
            pos: positions.len(),
            c6: c6.len(),
        });
    }
    if params.alpha <= 0.0 {
        return Err(LjPmeError::InvalidParam(format!(
            "alpha must be positive, got {}",
            params.alpha
        )));
    }
    if box_len[0] <= 0.0 || box_len[1] <= 0.0 || box_len[2] <= 0.0 {
        return Err(LjPmeError::InvalidParam(format!(
            "box dimensions must be positive, got {:?}",
            box_len
        )));
    }
    for &d in &params.grid {
        if d == 0 {
            return Err(LjPmeError::InvalidParam(
                "grid dimensions must be positive".to_string(),
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> LjPmeParams {
        LjPmeParams::new(0.35, [16, 16, 16], 10.0)
    }

    // -----------------------------------------------------------------------
    // Test 1: Self-energy matches analytical formula
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_self_energy() {
        let c6 = [10.0f64, 20.0, 15.0];
        let alpha = 0.35_f64;
        let sum_c6: f64 = c6.iter().sum();
        let expected = -(PI.powf(1.5) * alpha.powi(3) / 6.0) * sum_c6;
        let got = lj_pme_self_energy(&c6, alpha);
        assert!(
            (got - expected).abs() < 1e-12,
            "self-energy mismatch: got {got:.12} expected {expected:.12}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: Splitting function g₆ boundary values
    // -----------------------------------------------------------------------
    #[test]
    fn test_g6_limits() {
        // g₆(0) = exp(0) * (1 + 0 + 0) = 1
        assert!((g6(0.0) - 1.0).abs() < 1e-15, "g6(0) should be 1");
        // g₆(∞) → 0  (g₆(10) ≈ exp(-100)·5101 ≈ 1.9e-40, so use 1e-38 as guard)
        assert!(g6(10.0) < 1e-38, "g6(large) should be ~0");
        // g₆ ∈ (0, 1] for x > 0
        for &x in &[0.1, 0.5, 1.0, 2.0, 3.0] {
            let v = g6(x);
            assert!(v > 0.0 && v <= 1.0, "g6({x}) = {v} out of (0,1]");
        }
    }

    // -----------------------------------------------------------------------
    // Test 3: Reciprocal energy returns finite, non-NaN values
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_reciprocal_energy_finite() {
        let box_len = [15.0f64; 3];
        let positions = [[2.0, 3.0, 4.0], [8.0, 9.0, 10.0], [1.0, 12.0, 5.0]];
        let c6 = [1000.0f64, 800.0, 1200.0]; // kJ/mol Å⁶ (argon-like)
        let params = default_params();

        let e = lj_pme_reciprocal_energy(&positions, &c6, box_len, &params)
            .expect("reciprocal energy failed");
        assert!(!e.is_nan(), "energy is NaN");
        assert!(!e.is_infinite(), "energy is infinite");
    }

    // -----------------------------------------------------------------------
    // Test 4: Length mismatch returns error
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_length_mismatch() {
        let positions = vec![[1.0, 2.0, 3.0]; 3];
        let c6 = vec![1.0f64; 2]; // wrong length
        let params = default_params();
        let result = lj_pme_reciprocal_energy(&positions, &c6, [10.0; 3], &params);
        assert!(
            matches!(result, Err(LjPmeError::LengthMismatch { pos: 3, c6: 2 })),
            "Expected LengthMismatch, got: {result:?}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 5: Invalid alpha returns error
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_invalid_alpha() {
        let positions = vec![[1.0, 2.0, 3.0]];
        let c6 = vec![100.0f64];
        let params = LjPmeParams::new(-0.1, [8, 8, 8], 5.0);
        let result = lj_pme_reciprocal_energy(&positions, &c6, [10.0; 3], &params);
        assert!(matches!(result, Err(LjPmeError::InvalidParam(_))));
    }

    // -----------------------------------------------------------------------
    // Test 6: Real-space energy is negative for attractive dispersion
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_real_space_energy_sign() {
        let box_len = [20.0f64; 3];
        // Two atoms 4 Å apart (well within cutoff)
        let positions = [[5.0, 5.0, 5.0], [9.0, 5.0, 5.0]];
        let c6 = [1000.0f64, 1000.0]; // large C₆ → strong attraction
        let params = LjPmeParams::new(0.35, [16, 16, 16], 10.0);
        let e = lj_pme_real_space_energy(&positions, &c6, box_len, &params)
            .expect("real-space energy failed");
        assert!(e < 0.0, "dispersion real-space energy should be negative, got {e}");
    }

    // -----------------------------------------------------------------------
    // Test 7: Force-energy consistency via finite differences
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_forces_finite_diff() {
        let box_len = [12.0f64; 3];
        let positions = [
            [1.5, 2.0, 3.0],
            [6.0, 5.5, 4.0],
            [2.5, 9.0, 8.5],
        ];
        let c6 = [800.0f64, 1200.0, 1000.0];
        let params = default_params();
        let h = 1e-3_f64;

        let forces = lj_pme_reciprocal_forces(&positions, &c6, box_len, &params)
            .expect("forces failed");

        for atom in 0..positions.len() {
            for dim in 0..3usize {
                let mut pos_plus = positions.to_vec();
                let mut pos_minus = positions.to_vec();
                pos_plus[atom][dim] += h;
                pos_minus[atom][dim] -= h;

                let e_plus = lj_pme_reciprocal_energy(&pos_plus, &c6, box_len, &params)
                    .expect("energy+ failed");
                let e_minus = lj_pme_reciprocal_energy(&pos_minus, &c6, box_len, &params)
                    .expect("energy- failed");

                let f_fd = -(e_plus - e_minus) / (2.0 * h);
                let f_anal = forces[atom][dim];
                let err = (f_anal - f_fd).abs();
                let scale = f_fd.abs().max(1.0);
                assert!(
                    err / scale < 2e-2,
                    "Force mismatch atom={atom} dim={dim}: analytical={f_anal:.6} fd={f_fd:.6} rel_err={:.2e}",
                    err / scale
                );
            }
        }
    }
}
