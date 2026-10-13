// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LJ-PME (Lennard-Jones Particle Mesh Ewald) long-range dispersion.
//!
//! Implements the r⁻⁶ dispersion Ewald sum (real space + reciprocal space +
//! self term) following Essmann et al., J. Chem. Phys. 103, 8577 (1995), §V, and
//! in 't Veld, Ismail & Grest, J. Chem. Phys. 127, 144711 (2007).
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
//! ## Reciprocal-space energy
//! The Fourier transform of the smooth part `(1 - g₆(αr))/r⁶` is
//! `π^(3/2) α³ f₆(b)` with `b = |k|/(2α)` and
//! `f₆(b) = ⅓[(1 - 2b²) exp(-b²) + 2b³ √π erfc(b)]`, so that
//! ```text
//! E_recip = -(π^(3/2) α³ / 2V) Σ_k f₆(|k|/2α) |S(k)|²,   S(k) = Σᵢ √C₆ᵢ exp(i k·rᵢ)
//! ```
//! The sum includes `k = 0` (`f₆(0) = 1/3`): unlike the Coulomb case the
//! dispersion integral of the smooth part does not vanish.  On the mesh,
//! `|S(k)|²` is replaced by `B(k) |ρ̂(k)|²`, where `ρ̂` is the DFT of the spread
//! `√C₆ᵢ` and `B(k)` is the squared inverse modulus of the B-spline
//! interpolation.  As in `pme.rs`, the B-spline correction *multiplies*.
//!
//! ## Self-energy correction
//! The smooth part stays finite at `r → 0`, `(1 - g₆(αr))/r⁶ → α⁶/6`, so the
//! reciprocal sum contains a spurious `i = j` term that must be removed:
//! ```text
//! E_self = +(α⁶ / 12) Σᵢ C₆ᵢ
//! ```
//!
//! # Units
//! - Distance: Å
//! - C₆:       kJ mol⁻¹ Å⁶
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
// Reciprocal-space kernel  f₆(b) = ⅓[(1 - 2b²) e^(-b²) + 2b³ √π erfc(b)]
// ---------------------------------------------------------------------------

/// `√π`.
const SQRT_PI: f64 = 1.772_453_850_905_516;

/// Below this argument [`erfc_precise`] sums the power series of `erf`; from
/// here on it evaluates a continued fraction.
const ERFC_SERIES_LIMIT: f64 = 1.5;

/// Number of continued-fraction levels used by [`erfc_precise`].
const ERFC_CF_LEVELS: usize = 80;

/// Beyond this argument `erfc` underflows and `f6` is below `1e-38`.
const ERFC_ZERO_BEYOND: f64 = 27.3;

/// `f6(b)` is smaller than `1e-38` for `b` above this value and is returned as zero.
const F6_ZERO_BEYOND: f64 = 9.5;

/// Complementary error function with a relative accuracy of about `1e-13`.
///
/// The Abramowitz-Stegun 7.1.26 approximation used elsewhere in the crate is
/// only good to `1.5e-7` *absolute*, which the cancellation inside [`f6`]
/// amplifies beyond use.  For `0 <= x < 1.5` this sums the all-positive series
/// `erf(x) = 2/√π e^(-x²) Σₙ 2ⁿ x^(2n+1) / (2n+1)!!`; above that it evaluates
/// the continued fraction
/// `erfc(x) = e^(-x²)/√π / (x + (1/2)/(x + (2/2)/(x + (3/2)/(x + …))))`
/// from the tail.  Negative arguments use `erfc(-x) = 2 - erfc(x)`.
fn erfc_precise(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_precise(-x);
    }
    if x.is_nan() {
        return f64::NAN;
    }
    if x < ERFC_SERIES_LIMIT {
        let x2 = x * x;
        let mut term = x;
        let mut sum = x;
        for n in 1..=200_u32 {
            term *= 2.0 * x2 / f64::from(2 * n + 1);
            sum += term;
            if term <= 1.0e-17 * sum {
                break;
            }
        }
        return 1.0 - std::f64::consts::FRAC_2_SQRT_PI * (-x2).exp() * sum;
    }
    if x > ERFC_ZERO_BEYOND {
        return 0.0;
    }
    let mut tail = x;
    for level in (1..=ERFC_CF_LEVELS).rev() {
        tail = x + 0.5 * level as f64 / tail;
    }
    (-x * x).exp() / (SQRT_PI * tail)
}

/// Reciprocal-space kernel of the r⁻⁶ Ewald split,
/// `f₆(b) = ⅓[(1 - 2b²) e^(-b²) + 2b³ √π erfc(b)]` with `b = |k|/(2α)`.
///
/// The Fourier transform of the smooth part `(1 - g₆(αr))/r⁶` of the split is
/// `π^(3/2) α³ f₆(b)`.  `f₆(0) = 1/3` and `f₆` decays like `e^(-b²)/(2b²)`; it is
/// returned as zero for `b > 9.5`, where it is below `1e-38`.
pub fn f6(b: f64) -> f64 {
    let b = b.abs();
    if b > F6_ZERO_BEYOND {
        return 0.0;
    }
    let b2 = b * b;
    ((1.0 - 2.0 * b2) * (-b2).exp() + 2.0 * b2 * b * SQRT_PI * erfc_precise(b)) / 3.0
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

/// Per-axis reciprocal wavevector components `k = 2π n / L` and squared B-spline
/// moduli `|b(n)|²` for the FFT frequencies of an `m`-point axis of edge `edge`.
fn axis_modes(m: usize, edge: f64) -> Vec<(f64, f64)> {
    (0..m)
        .map(|index| {
            let n = if index <= m / 2 {
                index as i64
            } else {
                index as i64 - m as i64
            };
            (
                2.0 * PI * n as f64 / edge,
                bspline_correction_sq(n, m, BSPLINE_ORDER),
            )
        })
        .collect()
}

/// Compute the LJ-PME reciprocal-space energy (kJ mol⁻¹) and the potential
/// grid for force gathering.
///
/// Returns `(energy, pot_re)` where `pot_re` is the real-space potential grid
/// after IFFT of `G(k) · ρ̂(k)`, with the Green's function
/// `G(k) = -(π^(3/2) α³ / V) f₆(|k|/2α) B(k)` and `E = ½ Σ_k G(k) |ρ̂(k)|²`.
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

    // 3. Green's function and energy:
    //      G(k) = -(π^(3/2) α³ / V) · f₆(|k| / 2α) · B(k),   E = ½ Σ_k G(k) |ρ̂(k)|²
    //    B(k) = |b_x|² |b_y|² |b_z|² multiplies (it is the *inverse* squared modulus of
    //    the spline interpolation, as in pme.rs).  The k = 0 term is kept:
    //    f₆(0) = 1/3, the dispersion integral of the smooth part does not vanish.
    let modes_x = axis_modes(mx, box_len[0]);
    let modes_y = axis_modes(my, box_len[1]);
    let modes_z = axis_modes(mz, box_len[2]);
    let g_prefactor = -PI.powf(1.5) * alpha * alpha * alpha / volume;
    let two_alpha = 2.0 * alpha;

    let mut phi_re = vec![0.0f64; n_total];
    let mut phi_im = vec![0.0f64; n_total];
    let mut energy = 0.0f64;

    for (ix, &(kx, bx)) in modes_x.iter().enumerate() {
        for (iy, &(ky, by)) in modes_y.iter().enumerate() {
            let kxy_sq = kx * kx + ky * ky;
            let bxy = bx * by;
            for (iz, &(kz, bz)) in modes_z.iter().enumerate() {
                let b = (kxy_sq + kz * kz).sqrt() / two_alpha;
                let gk = g_prefactor * f6(b) * bxy * bz;

                let flat = ix * my * mz + iy * mz + iz;
                // ρ̂(k) = rho_re + i*rho_im  (output of FFT of real input)
                let r_re = rho_re[flat];
                let r_im = rho_im[flat];

                // |ρ̂(k)|² contributes to the energy
                energy += gk * (r_re * r_re + r_im * r_im);

                // φ̂(k) = G(k) * ρ̂(k)  for potential → forces
                phi_re[flat] = gk * r_re;
                phi_im[flat] = gk * r_im;
            }
        }
    }

    // Factor of 1/2: every pair is counted twice in |ρ̂|².
    energy *= 0.5;

    // 4. Inverse 3D FFT to get real-space potential grid φ(r)
    let (pot_re, _pot_im) = oxifft::ifft3d_split::<f64>(&phi_re, &phi_im, mx, my, mz);

    (energy, pot_re)
}

/// Gather the reciprocal-space forces (kJ mol⁻¹ Å⁻¹) from the potential grid
/// `pot_re` by B-spline interpolation of its gradient.
///
/// `F_i^α = -√C₆ᵢ · (n_total / L_α) · Σ_m (dM/du_α) M_β M_γ φ(m)`; the `n_total`
/// factor undoes the `1/N` normalisation of the inverse FFT.  This is the exact
/// gradient of the energy returned by `lj_pme_reciprocal_energy_and_potential`.
fn gather_reciprocal_forces(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    grid: [usize; 3],
    pot_re: &[f64],
) -> Vec<[f64; 3]> {
    let [mx, my, mz] = grid;
    let ntot_f64 = (mx * my * mz) as f64;
    let mut forces = vec![[0.0f64; 3]; positions.len()];

    for (force, (pos, &c6i)) in forces.iter_mut().zip(positions.iter().zip(c6.iter())) {
        let sqrt_c6i = c6i.sqrt();
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

        force[0] = -sqrt_c6i * ntot_f64 * fx / box_len[0];
        force[1] = -sqrt_c6i * ntot_f64 * fy / box_len[1];
        force[2] = -sqrt_c6i * ntot_f64 * fz / box_len[2];
    }

    forces
}

/// Real-space dispersion energy and forces with Ewald damping.
///
/// `E_real = -Σ_{i<j} C₆ᵢⱼ g₆(α rᵢⱼ) / rᵢⱼ⁶` over minimum-image pairs closer than
/// `r_cut`; the forces are its exact negative gradient.
fn real_space_energy_and_forces(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    alpha: f64,
    r_cut: f64,
) -> (f64, Vec<[f64; 3]>) {
    let r_cut_sq = r_cut * r_cut;
    let n = positions.len();
    let mut energy = 0.0f64;
    let mut forces = vec![[0.0f64; 3]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            // d points from atom i to atom j.
            let dx = min_image(positions[j][0] - positions[i][0], box_len[0]);
            let dy = min_image(positions[j][1] - positions[i][1], box_len[1]);
            let dz = min_image(positions[j][2] - positions[i][2], box_len[2]);
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > r_cut_sq || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            let r6 = r2 * r2 * r2;
            let c6ij = (c6[i] * c6[j]).sqrt();
            let s = alpha * r;
            let g = g6(s);
            energy -= c6ij * g / r6;

            // U(r) = -C g₆(αr) / r⁶
            // dU/dr = -C (α g₆'(αr) - 6 g₆(αr) / r) / r⁶   (positive: U rises towards 0)
            let dg = dg6_dx(s);
            let du_dr = -c6ij * (dg * alpha - 6.0 * g / r) / r6;
            // F_i = -dU/dr_i = +(dU/dr) (r_j - r_i) / r: an attractive pair
            // (dU/dr > 0) pulls atom i towards atom j.
            let f_scalar = du_dr / r;
            forces[i][0] += f_scalar * dx;
            forces[i][1] += f_scalar * dy;
            forces[i][2] += f_scalar * dz;
            forces[j][0] -= f_scalar * dx;
            forces[j][1] -= f_scalar * dy;
            forces[j][2] -= f_scalar * dz;
        }
    }

    (energy, forces)
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
/// or `LjPmeError::InvalidParam` if α or box dimensions are non-positive, the grid
/// has an empty axis, a C₆ is negative or non-finite, or a coordinate is not finite.
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
/// The smooth (reciprocal) part of the r⁻⁶ split stays finite at zero
/// separation, `(1 - g₆(αr))/r⁶ → α⁶/6`, so the reciprocal sum contains a
/// spurious interaction of every atom with itself, which is removed by:
/// ```text
/// E_self = +(α⁶ / 12) * Σᵢ C₆ᵢ
/// ```
///
/// # Arguments
/// * `c6`    — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N
/// * `alpha` — Ewald splitting parameter (Å⁻¹)
pub fn lj_pme_self_energy(c6: &[f64], alpha: f64) -> f64 {
    let alpha2 = alpha * alpha;
    let alpha6 = alpha2 * alpha2 * alpha2;
    let sum_c6: f64 = c6.iter().sum();
    alpha6 / 12.0 * sum_c6
}

/// Compute the real-space LJ dispersion energy (kJ mol⁻¹) with Ewald damping.
///
/// ```text
/// E_real = -Σ_{i<j}  C₆ᵢⱼ · g₆(α·rᵢⱼ) / rᵢⱼ⁶
/// ```
/// where `C₆ᵢⱼ = √(C₆ᵢ · C₆ⱼ)` (geometric combining rule) and
/// `g₆(x) = exp(-x²)(1 + x² + x⁴/2)`.  Pairs use the minimum image and must be
/// closer than `params.r_cut`.
///
/// # Arguments
/// * `positions`  — atom positions (Å), length N
/// * `c6`         — per-atom C₆ coefficients (kJ mol⁻¹ Å⁶), length N
/// * `box_len`    — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `params`     — LJ-PME parameters (uses `alpha` and `r_cut`)
///
/// # Errors
/// Returns `LjPmeError::LengthMismatch` if `positions.len() != c6.len()`, or
/// `LjPmeError::InvalidParam` for an invalid `alpha`, box, C₆ or an `r_cut` that
/// exceeds half of the shortest box edge (the minimum image would then miss pairs).
pub fn lj_pme_real_space_energy(
    positions: &[[f64; 3]],
    c6: &[f64],
    box_len: [f64; 3],
    params: &LjPmeParams,
) -> Result<f64, LjPmeError> {
    validate_inputs(positions, c6, box_len, params)?;
    validate_real_space_cutoff(box_len, params)?;
    let (energy, _forces) =
        real_space_energy_and_forces(positions, c6, box_len, params.alpha, params.r_cut);
    Ok(energy)
}

/// Compute the LJ-PME reciprocal-space forces (kJ mol⁻¹ Å⁻¹).
///
/// Gathers forces from the potential grid via B-spline interpolation of
/// the potential gradient.  They are the exact negative gradient of
/// [`lj_pme_reciprocal_energy`].
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
    let (_energy, pot_re) =
        lj_pme_reciprocal_energy_and_potential(positions, c6, box_len, params.alpha, params.grid);
    Ok(gather_reciprocal_forces(
        positions,
        c6,
        box_len,
        params.grid,
        &pot_re,
    ))
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
/// which reproduces the periodic lattice sum `-½ Σ_{i,j,n}' C₆ᵢⱼ / |rᵢⱼ + nL|⁶`
/// independently of the splitting parameter `α` (up to the cutoff and mesh error).
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
    validate_real_space_cutoff(box_len, params)?;
    let alpha = params.alpha;

    // Real-space energy and forces
    let (e_real, mut forces) =
        real_space_energy_and_forces(positions, c6, box_len, alpha, params.r_cut);

    // Reciprocal-space energy + potential grid for forces
    let (e_recip, pot_re) =
        lj_pme_reciprocal_energy_and_potential(positions, c6, box_len, alpha, params.grid);

    // Self-energy
    let e_self = lj_pme_self_energy(c6, alpha);

    // Gather reciprocal forces
    let recip_forces = gather_reciprocal_forces(positions, c6, box_len, params.grid, &pot_re);
    for (force, recip) in forces.iter_mut().zip(recip_forces.iter()) {
        for axis in 0..3 {
            force[axis] += recip[axis];
        }
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
    if !params.alpha.is_finite() || params.alpha <= 0.0 {
        return Err(LjPmeError::InvalidParam(format!(
            "alpha must be positive, got {}",
            params.alpha
        )));
    }
    if box_len.iter().any(|l| !l.is_finite() || *l <= 0.0) {
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
    if let Some(bad) = c6.iter().find(|c| !c.is_finite() || **c < 0.0) {
        return Err(LjPmeError::InvalidParam(format!(
            "C6 coefficients must be finite and non-negative, got {bad}"
        )));
    }
    if positions.iter().flatten().any(|x| !x.is_finite()) {
        return Err(LjPmeError::InvalidParam(
            "atom positions must be finite".to_string(),
        ));
    }
    Ok(())
}

/// The real-space sum uses the minimum image, which is complete only up to half
/// of the shortest box edge.
fn validate_real_space_cutoff(box_len: [f64; 3], params: &LjPmeParams) -> Result<(), LjPmeError> {
    let half_edge = 0.5 * box_len.iter().copied().fold(f64::INFINITY, f64::min);
    if !params.r_cut.is_finite() || params.r_cut <= 0.0 || params.r_cut > half_edge {
        return Err(LjPmeError::InvalidParam(format!(
            "r_cut must lie in (0, L_min/2 = {half_edge}], got {}",
            params.r_cut
        )));
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

    /// Cubic box edge of the eight-atom reference system (Å).
    const REF_BOX: f64 = 15.0;

    /// Eight atoms with argon-like C₆ in a 15 Å cubic box; no pair is closer than 3.3 Å.
    const REF_POSITIONS: [[f64; 3]; 8] = [
        [3.762, 14.201, 2.84],
        [2.689, 5.248, 3.458],
        [10.057, 1.726, 13.445],
        [12.872, 0.042, 8.122],
        [1.603, 3.869, 6.253],
        [6.804, 7.022, 13.913],
        [3.882, 2.818, 10.058],
        [14.199, 13.842, 13.204],
    ];

    /// Per-atom C₆ of the eight-atom reference system (kJ mol⁻¹ Å⁶).
    const REF_C6: [f64; 8] = [651.5, 1349.4, 1119.4, 1297.2, 926.5, 775.5, 1234.4, 1129.3];

    /// Direct periodic lattice sum `-½ Σ_{i,j} Σ_n' √(C₆ᵢ C₆ⱼ) / |rᵢⱼ + n L|⁶` over the
    /// images `|n_α| <= n_max`; the `i = j, n = 0` term is excluded.
    fn direct_dispersion_energy(positions: &[[f64; 3]], c6: &[f64], edge: f64, n_max: i32) -> f64 {
        let mut energy = 0.0;
        for (i, pi_pos) in positions.iter().enumerate() {
            for (j, pj_pos) in positions.iter().enumerate() {
                let c6ij = (c6[i] * c6[j]).sqrt();
                for nx in -n_max..=n_max {
                    for ny in -n_max..=n_max {
                        for nz in -n_max..=n_max {
                            if i == j && nx == 0 && ny == 0 && nz == 0 {
                                continue;
                            }
                            let dx = pj_pos[0] - pi_pos[0] + f64::from(nx) * edge;
                            let dy = pj_pos[1] - pi_pos[1] + f64::from(ny) * edge;
                            let dz = pj_pos[2] - pi_pos[2] + f64::from(nz) * edge;
                            let r2 = dx * dx + dy * dy + dz * dz;
                            energy -= 0.5 * c6ij / (r2 * r2 * r2);
                        }
                    }
                }
            }
        }
        energy
    }

    /// Exact (non-mesh) reciprocal sum `-(π^(3/2) α³ / 2V) Σ_k f₆(|k|/2α) |S(k)|²`
    /// over the integer vectors `|n_α| <= k_max`, with `S(k) = Σ √C₆ᵢ e^(i k·rᵢ)`.
    fn exact_reciprocal_energy(
        positions: &[[f64; 3]],
        c6: &[f64],
        edge: f64,
        alpha: f64,
        k_max: i32,
    ) -> f64 {
        let volume = edge * edge * edge;
        let unit = 2.0 * PI / edge;
        let mut sum = 0.0;
        for nx in -k_max..=k_max {
            for ny in -k_max..=k_max {
                for nz in -k_max..=k_max {
                    let k = [
                        unit * f64::from(nx),
                        unit * f64::from(ny),
                        unit * f64::from(nz),
                    ];
                    let k_norm = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
                    let mut s_re = 0.0;
                    let mut s_im = 0.0;
                    for (pos, &c6i) in positions.iter().zip(c6.iter()) {
                        let phase = k[0] * pos[0] + k[1] * pos[1] + k[2] * pos[2];
                        s_re += c6i.sqrt() * phase.cos();
                        s_im += c6i.sqrt() * phase.sin();
                    }
                    sum += f6(k_norm / (2.0 * alpha)) * (s_re * s_re + s_im * s_im);
                }
            }
        }
        -(PI.powf(1.5) * alpha.powi(3) / (2.0 * volume)) * sum
    }

    fn reference_params(alpha: f64, grid: usize) -> LjPmeParams {
        LjPmeParams::new(alpha, [grid; 3], 7.4)
    }

    fn reference_total_energy(alpha: f64, grid: usize) -> f64 {
        let (energy, _forces) = lj_pme_energy_and_forces(
            &REF_POSITIONS,
            &REF_C6,
            [REF_BOX; 3],
            &reference_params(alpha, grid),
        )
        .expect("reference system is valid");
        energy
    }

    // -----------------------------------------------------------------------
    // Test 1: Self-energy matches the analytical formula +α⁶/12 ΣC₆
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_self_energy() {
        let c6 = [10.0f64, 20.0, 15.0];
        let alpha = 0.35_f64;
        let sum_c6: f64 = c6.iter().sum();
        let expected = alpha.powi(6) / 12.0 * sum_c6;
        let got = lj_pme_self_energy(&c6, alpha);
        assert!(
            (got - expected).abs() < 1e-15,
            "self-energy mismatch: got {got:.12e} expected {expected:.12e}"
        );
        assert!(got > 0.0, "the removed i = j term is repulsive: {got}");
        // E_self carries the units of C₆ (kJ/mol Å⁶) times α⁶ (Å⁻⁶): it scales as α⁶.
        let doubled = lj_pme_self_energy(&c6, 2.0 * alpha);
        assert!(
            (doubled / got - 64.0).abs() < 1e-12,
            "ratio {}",
            doubled / got
        );
        assert_eq!(lj_pme_self_energy(&[], alpha), 0.0);
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
    // Test 3: erfc and the reciprocal kernel f₆ against high-precision references
    // -----------------------------------------------------------------------
    #[test]
    fn test_erfc_precise_matches_reference_values() {
        // 40-digit mpmath references, rounded to 15 significant digits.
        let reference = [
            (0.0, 1.0),
            (0.1, 0.887_537_083_981_715),
            (0.5, 0.479_500_122_186_953),
            (1.0, 0.157_299_207_050_285),
            (1.4999, 0.033_906_748_337_704_7),
            (1.5, 0.033_894_853_524_689_3),
            (2.0, 0.004_677_734_981_047_27),
            (3.0, 2.209_049_699_858_54e-5),
            (5.0, 1.537_459_794_428_03e-12),
            (10.0, 2.088_487_583_762_54e-45),
        ];
        for (x, expected) in reference {
            let got = erfc_precise(x);
            assert!(
                ((got - expected) / expected).abs() < 1e-12,
                "erfc({x}) = {got:e}, expected {expected:e}"
            );
        }
        // Reflection and the far tail.
        assert!((erfc_precise(-0.7) + erfc_precise(0.7) - 2.0).abs() < 1e-14);
        assert_eq!(erfc_precise(40.0), 0.0);
        assert!(erfc_precise(f64::NAN).is_nan());
    }

    #[test]
    fn test_f6_matches_reference_values() {
        // f₆(b) = ⅓[(1 - 2b²) e^(-b²) + 2b³ √π erfc(b)] (40-digit mpmath references).
        let reference = [
            (0.0, 1.0 / 3.0),
            (0.25, 0.287_356_706_617_941),
            (0.5, 0.200_624_450_351_895),
            (1.0, 0.063_243_909_796_627_2),
            (2.0, 0.001_482_545_956_541_14),
            (3.0, 5.456_733_366_220_39e-6),
            (4.0, 3.061_995_026_527_69e-9),
        ];
        for (b, expected) in reference {
            let got = f6(b);
            assert!(
                ((got - expected) / expected).abs() < 1e-9,
                "f6({b}) = {got:e}, expected {expected:e}"
            );
        }
        // Even in b, monotonically decreasing and vanishing for large arguments.
        assert_eq!(f6(-1.3), f6(1.3));
        let mut previous = f6(0.0);
        for step in 1..=60 {
            let value = f6(0.1 * f64::from(step));
            assert!(
                value < previous && value > 0.0,
                "f6 not decreasing at step {step}"
            );
            previous = value;
        }
        assert_eq!(f6(12.0), 0.0);
    }

    // -----------------------------------------------------------------------
    // Test 4: Reciprocal energy returns finite, non-NaN values
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
        // The smooth part of an attractive r⁻⁶ potential is attractive.
        assert!(
            e < 0.0,
            "reciprocal dispersion energy must be negative: {e}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 5: Length mismatch returns error
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
    // Test 6: Invalid alpha returns error
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_invalid_alpha() {
        let positions = vec![[1.0, 2.0, 3.0]];
        let c6 = vec![100.0f64];
        let params = LjPmeParams::new(-0.1, [8, 8, 8], 5.0);
        let result = lj_pme_reciprocal_energy(&positions, &c6, [10.0; 3], &params);
        assert!(matches!(result, Err(LjPmeError::InvalidParam(_))));
    }

    #[test]
    fn test_lj_pme_rejects_invalid_inputs() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let params = LjPmeParams::new(0.35, [8, 8, 8], 4.0);
        // Negative C₆ has no square root.
        let negative = lj_pme_reciprocal_energy(&positions, &[1.0, -1.0], [10.0; 3], &params);
        assert!(
            matches!(negative, Err(LjPmeError::InvalidParam(_))),
            "{negative:?}"
        );
        // Non-finite coordinates.
        let bad_positions = vec![[f64::NAN, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let nan = lj_pme_reciprocal_energy(&bad_positions, &[1.0, 1.0], [10.0; 3], &params);
        assert!(matches!(nan, Err(LjPmeError::InvalidParam(_))), "{nan:?}");
        // The minimum image is incomplete beyond half the shortest edge.
        let long_cut = LjPmeParams::new(0.35, [8, 8, 8], 5.5);
        let too_long = lj_pme_real_space_energy(&positions, &[1.0, 1.0], [10.0; 3], &long_cut);
        assert!(
            matches!(too_long, Err(LjPmeError::InvalidParam(_))),
            "{too_long:?}"
        );
        let total = lj_pme_energy_and_forces(&positions, &[1.0, 1.0], [10.0; 3], &long_cut);
        assert!(
            matches!(total, Err(LjPmeError::InvalidParam(_))),
            "{total:?}"
        );
        // r_cut == L/2 is allowed.
        let half = LjPmeParams::new(0.35, [8, 8, 8], 5.0);
        assert!(lj_pme_real_space_energy(&positions, &[1.0, 1.0], [10.0; 3], &half).is_ok());
    }

    // -----------------------------------------------------------------------
    // Test 7: Real-space energy is negative for attractive dispersion
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
        assert!(
            e < 0.0,
            "dispersion real-space energy should be negative, got {e}"
        );
        // Closed form: -C g₆(αr) / r⁶ for the single pair.
        let expected = -1000.0 * g6(0.35 * 4.0) / 4.0_f64.powi(6);
        assert!(
            (e - expected).abs() < 1e-12 * expected.abs(),
            "{e} vs {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 8: The mesh reciprocal energy reproduces the exact Ewald k-sum
    // -----------------------------------------------------------------------
    #[test]
    fn test_reciprocal_energy_matches_exact_ewald_k_sum() {
        let alpha = 0.45;
        let exact = exact_reciprocal_energy(&REF_POSITIONS, &REF_C6, REF_BOX, alpha, 14);
        let mesh = lj_pme_reciprocal_energy(
            &REF_POSITIONS,
            &REF_C6,
            [REF_BOX; 3],
            &reference_params(alpha, 48),
        )
        .expect("reciprocal energy");
        assert!(exact < -5.0, "exact reciprocal energy {exact}");
        assert!(
            ((mesh - exact) / exact).abs() < 1e-4,
            "PME reciprocal energy {mesh} vs exact k-sum {exact}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 9: real + reciprocal + self equals the direct lattice sum
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_energy_matches_direct_lattice_sum() {
        // 8 atoms in a 15 Å box; direct sum over 17³ periodic images (the r⁻⁶ tail
        // beyond them is below 1e-5 of the total).
        let direct = direct_dispersion_energy(&REF_POSITIONS, &REF_C6, REF_BOX, 8);
        assert!(direct < -1.0 && direct > -2.0, "direct sum {direct}");
        let total = reference_total_energy(0.45, 48);
        assert!(
            ((total - direct) / direct).abs() < 2e-4,
            "LJ-PME total {total} vs direct lattice sum {direct}"
        );
    }

    #[test]
    fn test_total_energy_is_independent_of_the_splitting_parameter() {
        let direct = direct_dispersion_energy(&REF_POSITIONS, &REF_C6, REF_BOX, 8);
        let mut totals = Vec::new();
        for alpha in [0.45, 0.5, 0.55] {
            let total = reference_total_energy(alpha, 48);
            assert!(
                ((total - direct) / direct).abs() < 1e-3,
                "alpha = {alpha}: total {total} vs direct {direct}"
            );
            totals.push(total);
        }
        // The pieces themselves vary by orders of magnitude over this alpha range.
        let self_small = lj_pme_self_energy(&REF_C6, 0.45);
        let self_large = lj_pme_self_energy(&REF_C6, 0.55);
        assert!(self_large > 3.0 * self_small);
        let spread = totals.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - totals.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(spread < 1e-3 * direct.abs(), "alpha-dependence {spread}");
    }

    #[test]
    fn test_total_energy_converges_with_the_mesh() {
        let direct = direct_dispersion_energy(&REF_POSITIONS, &REF_C6, REF_BOX, 8);
        let errors: Vec<f64> = [24usize, 32, 48]
            .iter()
            .map(|&grid| ((reference_total_energy(0.45, grid) - direct) / direct).abs())
            .collect();
        assert!(
            errors[0] > errors[1] && errors[1] > errors[2],
            "errors must shrink with the mesh: {errors:?}"
        );
        assert!(errors[2] < 2e-4, "fine-mesh error {}", errors[2]);
    }

    #[test]
    fn test_energy_scales_linearly_with_c6() {
        let scaled: Vec<f64> = REF_C6.iter().map(|c| 3.0 * c).collect();
        let (e1, _) = lj_pme_energy_and_forces(
            &REF_POSITIONS,
            &REF_C6,
            [REF_BOX; 3],
            &reference_params(0.5, 24),
        )
        .expect("energy");
        let (e3, _) = lj_pme_energy_and_forces(
            &REF_POSITIONS,
            &scaled,
            [REF_BOX; 3],
            &reference_params(0.5, 24),
        )
        .expect("energy");
        assert!(((e3 / e1) - 3.0).abs() < 1e-12, "{}", e3 / e1);
    }

    #[test]
    fn test_energy_is_invariant_under_a_uniform_translation() {
        let shifted: Vec<[f64; 3]> = REF_POSITIONS
            .iter()
            .map(|p| [p[0] + 4.321, p[1] - 2.5, p[2] + 17.0])
            .collect();
        let params = reference_params(0.45, 48);
        let (e0, _) = lj_pme_energy_and_forces(&REF_POSITIONS, &REF_C6, [REF_BOX; 3], &params)
            .expect("energy");
        let (e1, _) =
            lj_pme_energy_and_forces(&shifted, &REF_C6, [REF_BOX; 3], &params).expect("energy");
        // PME energy depends on the position relative to the mesh only at the
        // interpolation-error level (below 1e-4 of the total for this mesh).
        assert!(((e1 - e0) / e0).abs() < 4e-4, "{e0} vs {e1}");
    }

    // -----------------------------------------------------------------------
    // Test 10: Force-energy consistency via finite differences
    // -----------------------------------------------------------------------
    #[test]
    fn test_lj_pme_forces_finite_diff() {
        let box_len = [12.0f64; 3];
        let positions = [[1.5, 2.0, 3.0], [6.0, 5.5, 4.0], [2.5, 9.0, 8.5]];
        let c6 = [800.0f64, 1200.0, 1000.0];
        let params = default_params();
        let h = 1e-4_f64;

        let forces =
            lj_pme_reciprocal_forces(&positions, &c6, box_len, &params).expect("forces failed");

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
                let scale = f_fd.abs().max(1e-3);
                assert!(
                    err / scale < 1e-5,
                    "Force mismatch atom={atom} dim={dim}: analytical={f_anal:.9} fd={f_fd:.9} rel_err={:.2e}",
                    err / scale
                );
            }
        }
    }

    #[test]
    fn test_real_space_forces_are_the_negative_energy_gradient() {
        // Four atoms inside the cutoff; the pair 0-1 is only 3.4 Å apart so its force
        // is large, and the pair 2-3 straddles the periodic boundary.
        let box_len = [14.0f64; 3];
        let positions = [
            [2.0, 2.0, 2.0],
            [5.4, 2.0, 2.0],
            [0.5, 7.0, 7.0],
            [13.0, 8.0, 6.5],
        ];
        let c6 = [900.0f64, 1100.0, 1000.0, 700.0];
        let (alpha, r_cut) = (0.4, 6.9);
        let (_energy, forces) =
            real_space_energy_and_forces(&positions, &c6, box_len, alpha, r_cut);
        let h = 1e-5;
        for atom in 0..positions.len() {
            for dim in 0..3usize {
                let mut plus = positions;
                let mut minus = positions;
                plus[atom][dim] += h;
                minus[atom][dim] -= h;
                let e_plus = real_space_energy_and_forces(&plus, &c6, box_len, alpha, r_cut).0;
                let e_minus = real_space_energy_and_forces(&minus, &c6, box_len, alpha, r_cut).0;
                let numeric = -(e_plus - e_minus) / (2.0 * h);
                let analytic = forces[atom][dim];
                assert!(
                    (analytic - numeric).abs() < 1e-7 * numeric.abs().max(1e-3),
                    "atom {atom} axis {dim}: analytic {analytic} vs finite difference {numeric}"
                );
            }
        }
        // Atom 0 is attracted towards atom 1 (+x).
        assert!(forces[0][0] > 0.0 && forces[1][0] < 0.0, "{:?}", forces);
        // Newton's third law.
        for axis in 0..3 {
            let net: f64 = forces.iter().map(|f| f[axis]).sum();
            assert!(
                net.abs() < 1e-12,
                "net real-space force {net} on axis {axis}"
            );
        }
    }

    #[test]
    fn test_total_forces_are_the_negative_gradient_of_the_total_energy() {
        // Real-space, mesh and self terms together: the unified entry point must
        // return forces that are the exact derivative of the energy it returns.
        let params = LjPmeParams::new(0.55, [24, 24, 24], 7.0);
        let box_len = [REF_BOX; 3];
        let (_energy, forces) =
            lj_pme_energy_and_forces(&REF_POSITIONS, &REF_C6, box_len, &params).expect("forces");
        let h = 1e-4;
        for atom in 0..REF_POSITIONS.len() {
            for dim in 0..3usize {
                let mut plus = REF_POSITIONS;
                let mut minus = REF_POSITIONS;
                plus[atom][dim] += h;
                minus[atom][dim] -= h;
                let e_plus = lj_pme_energy_and_forces(&plus, &REF_C6, box_len, &params)
                    .expect("energy+")
                    .0;
                let e_minus = lj_pme_energy_and_forces(&minus, &REF_C6, box_len, &params)
                    .expect("energy-")
                    .0;
                let numeric = -(e_plus - e_minus) / (2.0 * h);
                let analytic = forces[atom][dim];
                assert!(
                    (analytic - numeric).abs() < 2e-6 * numeric.abs().max(1e-2),
                    "atom {atom} axis {dim}: analytic {analytic} vs finite difference {numeric}"
                );
            }
        }
    }
}
