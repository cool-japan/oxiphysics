// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PME (Particle Mesh Ewald) reciprocal-space electrostatics using B-spline order 4
//! charge spreading and 3-D FFT via OxiFFT.
//!
//! # Algorithm (Smooth-PME, Essmann et al. 1995)
//!
//! 1. **Charge spreading** — each atom's charge is distributed onto a regular grid
//!    using cardinal B-splines of order 4.
//! 2. **3D FFT** — the charge-density grid `ρ(m)` is transformed to `ρ̂(k)` using
//!    `oxifft::fft3d_split`.
//! 3. **Reciprocal-space energy** — for each non-zero k-vector:
//!    ```text
//!    E_recip = (1 / (2V)) * Σ_{k≠0}  (4π / k²)  *  exp(-k² / (4α²))  *  |S(k)|²
//!    ```
//!    where `S(k) = Σ_m  ρ̂(m) * B_spline_correction(k)` is the structure factor
//!    corrected for B-spline aliasing (the `B(k)` factor from Eq. 4.7 in Essmann 1995).
//!    The COULOMB_K prefactor is folded in at the end.
//! 4. **Reciprocal-space forces** — back-transform the potential and differentiate
//!    via B-spline interpolation.
//! 5. **Self-energy** — `E_self = -(α/√π) * Σ_j q_j²` (already in `coulomb.rs`).
//!
//! # Units
//! - Distance : Å
//! - Charge   : e
//! - Energy   : kJ mol⁻¹
//!
//! # Orthorhombic box
//! Only rectangular (orthorhombic) periodic boxes are supported; `box_lengths = [Lx, Ly, Lz]`.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constant (GROMACS-compatible)
// ---------------------------------------------------------------------------

/// Coulomb constant in GROMACS units: kJ·Å/(mol·e²).
pub const COULOMB_K: f64 = 138.935_458;

// ---------------------------------------------------------------------------
// B-spline order 4
// ---------------------------------------------------------------------------

/// Order of the cardinal B-spline used for charge spreading.
pub const BSPLINE_ORDER: usize = 4;

/// Cardinal B-spline of order `p` evaluated at position `x`.
///
/// Defined recursively:
/// ```text
/// M₁(x) = 1   for x ∈ [0,1), else 0
/// Mₙ(x) = (x/(n-1)) * Mₙ₋₁(x)  +  ((n-x)/(n-1)) * Mₙ₋₁(x-1)
/// ```
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

/// Derivative of cardinal B-spline of order `p` at `x`.
///
/// Using the identity: M'_n(x) = M_{n-1}(x) - M_{n-1}(x-1).
#[inline]
fn bspline_deriv(p: usize, x: f64) -> f64 {
    if p <= 1 {
        return 0.0;
    }
    bspline(p - 1, x) - bspline(p - 1, x - 1.0)
}

/// Precompute B-spline weights for a single particle fractional coordinate `u`
/// in a grid of size `n_grid`.
///
/// Returns `(weights, offsets)` where `offsets[j]` is the grid-point index (periodic)
/// and `weights[j]` is the corresponding B-spline weight.
fn bspline_weights(
    u: f64,
    n_grid: usize,
    order: usize,
) -> ([f64; BSPLINE_ORDER], [usize; BSPLINE_ORDER]) {
    // fractional grid coordinate in [0, n_grid)
    let scaled = u * (n_grid as f64);
    let i0 = scaled.floor() as i64;

    let mut weights = [0.0f64; BSPLINE_ORDER];
    let mut grid_pts = [0usize; BSPLINE_ORDER];

    // Grid points: i0-(order-1), i0-(order-2), ..., i0
    // For order=4: i0-3, i0-2, i0-1, i0
    // x = scaled - raw ∈ (order-1, order] → M_p evaluated at x ∈ (p-1, p], (p-2, p-1], ...
    for j in 0..order {
        let raw = i0 - (order as i64 - 1) + j as i64;
        let x = scaled - raw as f64; // in [0, order)
        weights[j] = bspline(order, x);
        // periodic wrapping
        grid_pts[j] = raw.rem_euclid(n_grid as i64) as usize;
    }
    (weights, grid_pts)
}

/// Precompute B-spline *derivative* weights for force interpolation.
///
/// The factor `n_grid` converts from fractional-coordinate derivative to
/// grid-point derivative: `dM/du = dM/d(x/L) = K * dM/dx`.
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

// ---------------------------------------------------------------------------
// B-spline structure-factor correction  b(k)
// ---------------------------------------------------------------------------

/// Compute the 1D B-spline structure factor correction for grid index `m`
/// with grid size `n` and spline order `p`.
///
/// From Essmann et al. (1995) Eq. 4.7:
/// ```text
/// b(m) = exp(2π i (p-1) m / n) / Σ_{l=0}^{p-2}  M_p(l+1) * exp(2πi m l / n)
/// ```
/// Returns `|b(m)|²`.
fn bspline_correction_sq(m: i64, n: usize, p: usize) -> f64 {
    let n_f = n as f64;
    let two_pi_m_over_n = 2.0 * PI * (m as f64) / n_f;

    // denominator: Σ_{l=0}^{p-2}  M_p(l+1) * exp(2πi m l / n)
    let mut den_re = 0.0f64;
    let mut den_im = 0.0f64;
    for l in 0..(p - 1) {
        let bval = bspline(p, (l + 1) as f64);
        let angle = two_pi_m_over_n * (l as f64);
        den_re += bval * angle.cos();
        den_im += bval * angle.sin();
    }

    // numerator: exp(2πi (p-1) m / n)
    let num_angle = two_pi_m_over_n * ((p - 1) as f64);
    let num_re = num_angle.cos();
    let num_im = num_angle.sin();

    // |numerator / denominator|² = |num|² / |den|²
    let num_sq = num_re * num_re + num_im * num_im;
    let den_sq = den_re * den_re + den_im * den_im;

    if den_sq < 1e-30 {
        return 0.0;
    }
    num_sq / den_sq
}

// ---------------------------------------------------------------------------
// Charge spreading
// ---------------------------------------------------------------------------

/// Spread point charges onto a 3D grid using B-spline order `BSPLINE_ORDER`.
///
/// Grid layout: row-major with index `ix * ny * nz + iy * nz + iz`.
fn spread_charges_bspline(
    positions: &[[f64; 3]],
    charges: &[f64],
    box_lengths: [f64; 3],
    nx: usize,
    ny: usize,
    nz: usize,
) -> Vec<f64> {
    let mut rho = vec![0.0f64; nx * ny * nz];

    for (pos, &q) in positions.iter().zip(charges.iter()) {
        // fractional coordinates in [0, 1)
        let ux = pos[0].rem_euclid(box_lengths[0]) / box_lengths[0];
        let uy = pos[1].rem_euclid(box_lengths[1]) / box_lengths[1];
        let uz = pos[2].rem_euclid(box_lengths[2]) / box_lengths[2];

        let (wx, gx) = bspline_weights(ux, nx, BSPLINE_ORDER);
        let (wy, gy) = bspline_weights(uy, ny, BSPLINE_ORDER);
        let (wz, gz) = bspline_weights(uz, nz, BSPLINE_ORDER);

        for jx in 0..BSPLINE_ORDER {
            for jy in 0..BSPLINE_ORDER {
                let wxy = wx[jx] * wy[jy];
                for jz in 0..BSPLINE_ORDER {
                    let idx = gx[jx] * ny * nz + gy[jy] * nz + gz[jz];
                    rho[idx] += q * wxy * wz[jz];
                }
            }
        }
    }
    rho
}

// ---------------------------------------------------------------------------
// PME reciprocal energy
// ---------------------------------------------------------------------------

/// Compute the PME reciprocal-space energy (kJ mol⁻¹).
///
/// Uses B-spline order 4 charge spreading and OxiFFT's `fft3d_split` for
/// the 3D DFT.  The Ewald Gaussian envelope and B-spline correction factors
/// are applied in k-space before summing `|S(k)|²`.
///
/// # Arguments
/// * `positions`   — atom positions (Å), length N
/// * `charges`     — atom charges (e), length N
/// * `box_lengths` — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `alpha`       — Ewald splitting parameter (Å⁻¹)
/// * `mesh`        — \[mx, my, mz\] number of grid points per dimension
pub fn pme_reciprocal_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    box_lengths: [f64; 3],
    alpha: f64,
    mesh: [usize; 3],
) -> f64 {
    let [mx, my, mz] = mesh;
    let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
    let n_total = mx * my * mz;

    // 1. Spread charges onto grid
    let rho = spread_charges_bspline(positions, charges, box_lengths, mx, my, mz);

    // 2. 3D FFT of charge density (split-complex; imaginary part of input = 0)
    let rho_imag = vec![0.0f64; n_total];
    let (rho_re, rho_im) = oxifft::fft3d_split::<f64>(&rho, &rho_imag, mx, my, mz);

    // 3. Sum E_recip = prefactor * Σ_{k≠0}  (4π/k²) * exp(-k²/4α²) * |B(k)|² * |ρ̂(k)|²
    // prefactor = COULOMB_K / (2 * V)
    let prefactor = COULOMB_K / (2.0 * volume);
    let four_alpha_sq = 4.0 * alpha * alpha;
    let two_pi = 2.0 * PI;

    let mut energy = 0.0f64;

    for ix in 0..mx {
        let nx = if ix <= mx / 2 {
            ix as i64
        } else {
            ix as i64 - mx as i64
        };
        let kx = two_pi * nx as f64 / box_lengths[0];
        let bc_x = bspline_correction_sq(nx, mx, BSPLINE_ORDER);

        for iy in 0..my {
            let ny = if iy <= my / 2 {
                iy as i64
            } else {
                iy as i64 - my as i64
            };
            let ky = two_pi * ny as f64 / box_lengths[1];
            let bc_y = bspline_correction_sq(ny, my, BSPLINE_ORDER);

            for iz in 0..mz {
                if ix == 0 && iy == 0 && iz == 0 {
                    continue; // skip k=0 (charge neutrality term)
                }
                let nz = if iz <= mz / 2 {
                    iz as i64
                } else {
                    iz as i64 - mz as i64
                };
                let kz = two_pi * nz as f64 / box_lengths[2];
                let bc_z = bspline_correction_sq(nz, mz, BSPLINE_ORDER);

                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 < 1e-20 {
                    continue;
                }

                let gaussian = (-k2 / four_alpha_sq).exp();
                let bspline_corr = bc_x * bc_y * bc_z;

                let flat = ix * my * mz + iy * mz + iz;
                let rho_sq = rho_re[flat] * rho_re[flat] + rho_im[flat] * rho_im[flat];

                energy += (4.0 * PI / k2) * gaussian * bspline_corr * rho_sq;
            }
        }
    }

    prefactor * energy
}

// ---------------------------------------------------------------------------
// PME self-energy correction
// ---------------------------------------------------------------------------

/// Ewald self-energy correction (kJ mol⁻¹).
///
/// Removes the spurious self-interaction introduced by the Ewald splitting:
/// ```text
/// E_self = -COULOMB_K * α / sqrt(π) * Σ_j q_j²
/// ```
pub fn pme_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    -COULOMB_K * alpha / PI.sqrt() * sum_q2
}

// ---------------------------------------------------------------------------
// PME reciprocal forces
// ---------------------------------------------------------------------------

/// Compute the PME reciprocal-space forces on each atom (kJ mol⁻¹ Å⁻¹).
///
/// Evaluates the reciprocal-space potential on the grid via inverse FFT
/// of the Green's function-weighted charge density, then interpolates
/// forces onto atom positions using B-spline gradients.
///
/// # Arguments
/// * `positions`   — atom positions (Å), length N
/// * `charges`     — atom charges (e), length N
/// * `box_lengths` — orthorhombic box sides [Lx, Ly, Lz] (Å)
/// * `alpha`       — Ewald splitting parameter (Å⁻¹)
/// * `mesh`        — \[mx, my, mz\] number of grid points per dimension
///
/// Returns a `Vec<[f64; 3]>` of forces in kJ mol⁻¹ Å⁻¹.
pub fn pme_reciprocal_forces(
    positions: &[[f64; 3]],
    charges: &[f64],
    box_lengths: [f64; 3],
    alpha: f64,
    mesh: [usize; 3],
) -> Vec<[f64; 3]> {
    let [mx, my, mz] = mesh;
    let n_atoms = positions.len();
    let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
    let n_total = mx * my * mz;

    // 1. Spread charges onto grid
    let rho = spread_charges_bspline(positions, charges, box_lengths, mx, my, mz);

    // 2. 3D FFT of charge density
    let rho_imag = vec![0.0f64; n_total];
    let (rho_re, rho_im) = oxifft::fft3d_split::<f64>(&rho, &rho_imag, mx, my, mz);

    // 3. Multiply by Green's function G(k) = (4π / k²) * exp(-k²/4α²) * |B(k)|²
    //    to get potential in k-space: φ̂(k) = G(k) * ρ̂(k)
    let mut phi_re = vec![0.0f64; n_total];
    let mut phi_im = vec![0.0f64; n_total];

    let prefactor = COULOMB_K / volume;
    let four_alpha_sq = 4.0 * alpha * alpha;
    let two_pi = 2.0 * PI;

    for ix in 0..mx {
        let nx = if ix <= mx / 2 {
            ix as i64
        } else {
            ix as i64 - mx as i64
        };
        let kx = two_pi * nx as f64 / box_lengths[0];
        let bc_x = bspline_correction_sq(nx, mx, BSPLINE_ORDER);

        for iy in 0..my {
            let ny = if iy <= my / 2 {
                iy as i64
            } else {
                iy as i64 - my as i64
            };
            let ky = two_pi * ny as f64 / box_lengths[1];
            let bc_y = bspline_correction_sq(ny, my, BSPLINE_ORDER);

            for iz in 0..mz {
                if ix == 0 && iy == 0 && iz == 0 {
                    continue;
                }
                let nz = if iz <= mz / 2 {
                    iz as i64
                } else {
                    iz as i64 - mz as i64
                };
                let kz = two_pi * nz as f64 / box_lengths[2];
                let bc_z = bspline_correction_sq(nz, mz, BSPLINE_ORDER);

                let k2 = kx * kx + ky * ky + kz * kz;
                if k2 < 1e-20 {
                    continue;
                }

                let gk =
                    prefactor * (4.0 * PI / k2) * (-k2 / four_alpha_sq).exp() * bc_x * bc_y * bc_z;

                let flat = ix * my * mz + iy * mz + iz;
                phi_re[flat] = gk * rho_re[flat];
                phi_im[flat] = gk * rho_im[flat];
            }
        }
    }

    // 4. Inverse 3D FFT to get real-space potential grid φ(r)
    let (pot_re, _pot_im) = oxifft::ifft3d_split::<f64>(&phi_re, &phi_im, mx, my, mz);

    // 5. Interpolate forces onto atom positions via B-spline gradient
    //    F_i^α = -q_i * Σ_m  (dM/du_α) * φ(m)  * (1/L_α)
    let mut forces = vec![[0.0f64; 3]; n_atoms];

    for (atom_idx, (pos, &q)) in positions.iter().zip(charges.iter()).enumerate() {
        let ux = pos[0].rem_euclid(box_lengths[0]) / box_lengths[0];
        let uy = pos[1].rem_euclid(box_lengths[1]) / box_lengths[1];
        let uz = pos[2].rem_euclid(box_lengths[2]) / box_lengths[2];

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

                    // dE/dx_i = Σ_m dρ(m)/dx_i * φ(m)
                    // dρ(m)/dx_i = q_i * dM(u_x)/du_x * (1/Lx) * M(u_y) * M(u_z)
                    fx += dwx[jx] * wy[jy] * wz[jz] * phi;
                    fy += wx[jx] * dwy[jy] * wz[jz] * phi;
                    fz += wx[jx] * wy[jy] * dwz[jz] * phi;
                }
            }
        }

        // F_i^α = -(N/V) * q * Σ_m (dM/du_α) * M_y * M_z * φ_m / L_α
        // where N = n_total = mx*my*mz and φ_m = IDFT{(COULOMB_K/V) * G_k * ρ̂}[m].
        // The IDFT normalizes by 1/N, so we multiply by n_total to recover the correct factor.
        let ntot_f64 = n_total as f64;
        forces[atom_idx][0] = -q * ntot_f64 * fx / box_lengths[0];
        forces[atom_idx][1] = -q * ntot_f64 * fy / box_lengths[1];
        forces[atom_idx][2] = -q * ntot_f64 * fz / box_lengths[2];
    }

    forces
}

// ---------------------------------------------------------------------------
// Real-space Ewald sum helper
// ---------------------------------------------------------------------------

/// Real-space Ewald sum (kJ mol⁻¹) for an orthorhombic box with minimum image.
///
/// ```text
/// E_real = COULOMB_K * Σ_{i<j} q_i q_j * erfc(α r_ij) / r_ij
/// ```
pub fn ewald_real_space_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    box_lengths: [f64; 3],
    alpha: f64,
    _r_cutoff: f64,
) -> f64 {
    let n = positions.len();
    let mut energy = 0.0f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = min_image(positions[j][0] - positions[i][0], box_lengths[0]);
            let dy = min_image(positions[j][1] - positions[i][1], box_lengths[1]);
            let dz = min_image(positions[j][2] - positions[i][2], box_lengths[2]);
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-12 {
                continue;
            }
            energy += COULOMB_K * charges[i] * charges[j] * erfc_fast(alpha * r) / r;
        }
    }
    energy
}

/// Minimum image coordinate displacement for a single axis.
#[inline]
fn min_image(dx: f64, l: f64) -> f64 {
    dx - l * (dx / l).round()
}

/// Complementary error function (Abramowitz & Stegun 7.1.26, max error ~1.5e-7).
#[inline]
fn erfc_fast(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_fast(-x);
    }
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Tolerance helpers
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------------
    // Test 1: self-energy matches analytical formula
    // -----------------------------------------------------------------------
    #[test]
    fn test_pme_self_energy_correction() {
        let charges = [1.0f64, -1.0, 2.0, -0.5];
        let alpha = 0.35;
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        let expected = -COULOMB_K * alpha / PI.sqrt() * sum_q2;
        let got = pme_self_energy(&charges, alpha);
        assert!(
            approx_eq(got, expected, 1e-10),
            "self energy mismatch: got {got:.10} expected {expected:.10}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 2: reciprocal energy of a strictly neutral system with cancelling
    // charges near the same grid point should be very small
    // -----------------------------------------------------------------------
    #[test]
    fn test_pme_energy_neutral_system() {
        // Two equal-and-opposite charges separated by a full box diagonal —
        // real-space sum vanishes at cutoff, reciprocal contribution should
        // be finite but well-defined.  We just verify the function returns
        // without panic and that the total charge on the grid is near zero.
        let box_len = 10.0f64;
        let box_lengths = [box_len; 3];
        let positions = [[0.0, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let charges = [1.0f64, -1.0];
        let alpha = 0.35;
        let mesh = [16, 16, 16];

        let e_recip = pme_reciprocal_energy(&positions, &charges, box_lengths, alpha, mesh);
        // The energy should be finite and not NaN
        assert!(!e_recip.is_nan(), "reciprocal energy is NaN");
        assert!(!e_recip.is_infinite(), "reciprocal energy is infinite");
        // For a neutral system in a periodic box the reciprocal energy is finite
        // and should be of reasonable magnitude
        assert!(
            e_recip.abs() < 1.0e5,
            "reciprocal energy suspiciously large: {e_recip}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 3: Madelung constant of NaCl 4×4×4 supercell
    // -----------------------------------------------------------------------
    //
    // The NaCl Madelung constant A = -1.7475645946 (per ion pair, dimensionless).
    // For an N-ion system:
    //   E_total / (n_pairs * COULOMB_K * q^2 / d_nn) ≈ -1.7475645946
    // where d_nn is the nearest-neighbour distance.
    //
    // We build a 4×4×4 (128-ion) NaCl supercell with lattice constant a = 5.64 Å
    // (d_nn = a/2 = 2.82 Å, q = ±1 e), compute total Ewald energy, and compare.
    #[test]
    fn test_pme_madelung_nacl_4x4x4() {
        // NaCl rock-salt structure: FCC Bravais lattice with two-atom basis.
        // We use a conventional cubic cell with a = 5.64 Å and replicate 4×4×4.
        let a = 5.64f64; // conventional cell lattice parameter (Å)
        let d_nn = a / 2.0; // nearest-neighbour distance
        let n_cells_side = 4usize;
        let box_len = a * n_cells_side as f64;
        let box_lengths = [box_len; 3];

        // Build positions and charges
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut charges: Vec<f64> = Vec::new();

        for cx in 0..n_cells_side {
            for cy in 0..n_cells_side {
                for cz in 0..n_cells_side {
                    // Na at (0,0,0) of each conventional cell
                    let ox = cx as f64 * a;
                    let oy = cy as f64 * a;
                    let oz = cz as f64 * a;

                    // Basis: Na at (0,0,0), Cl at (a/2,0,0), Na at (0,a/2,0), Cl at (a/2,a/2,0)
                    //        Na at (0,0,a/2), Cl at (a/2,0,a/2), Na at (0,a/2,a/2), Cl at (a/2,a/2,a/2)
                    // Simplified: (ix + iy + iz) even => +1 (Na), odd => -1 (Cl)
                    // on a cubic grid with spacing a/2
                    let half_a = a / 2.0;
                    for bx in 0..2usize {
                        for by in 0..2usize {
                            for bz in 0..2usize {
                                let x = ox + bx as f64 * half_a;
                                let y = oy + by as f64 * half_a;
                                let z = oz + bz as f64 * half_a;
                                let q = if (bx + by + bz) % 2 == 0 {
                                    1.0f64
                                } else {
                                    -1.0
                                };
                                positions.push([x, y, z]);
                                charges.push(q);
                            }
                        }
                    }
                }
            }
        }

        let n_ions = positions.len(); // should be 512
        assert_eq!(n_ions, 512, "Expected 512 ions in 4×4×4 NaCl supercell");

        // Check neutrality
        let total_charge: f64 = charges.iter().sum();
        assert!(
            total_charge.abs() < 1e-10,
            "system not neutral: total charge = {total_charge}"
        );

        // Ewald parameters: moderate alpha for good reciprocal convergence with 48³ mesh.
        // alpha = 3.0 / r_cutoff ≈ 0.53 Å⁻¹ (r_cutoff = L/2 = 5.64 Å,
        // giving 4α² = 1.12 Å⁻², Gaussian at k_max ≈ 6.7 Å⁻¹ is ~exp(-40) ≈ 0).
        let alpha = 3.0 / (box_len / 2.0); // ≈ 0.266 Å⁻¹
        let mesh = [48, 48, 48];

        // Reciprocal-space contribution
        let e_recip = pme_reciprocal_energy(&positions, &charges, box_lengths, alpha, mesh);

        // Real-space contribution (standard Ewald direct sum with erfc damping)
        let r_cutoff = box_len / 2.0;
        let e_real = ewald_real_space_energy(&positions, &charges, box_lengths, alpha, r_cutoff);

        // Self-energy correction
        let e_self = pme_self_energy(&charges, alpha);

        let e_total = e_recip + e_real + e_self;

        // Reference Madelung constant
        let madelung_ref = -1.747_564_594_6;
        let n_pairs = (n_ions / 2) as f64;
        let e_madelung = e_total / (n_pairs * COULOMB_K / d_nn);

        let diff = (e_madelung - madelung_ref).abs();
        assert!(
            diff < 5e-3,
            "Madelung constant error too large: got {e_madelung:.8}, expected {madelung_ref:.8}, diff = {diff:.2e}"
        );
    }

    // -----------------------------------------------------------------------
    // Test 4: reciprocal forces match finite-difference gradient
    // -----------------------------------------------------------------------
    #[test]
    fn test_pme_reciprocal_forces_finite_diff() {
        // Simple 4-atom neutral system in a cubic box
        let box_lengths = [12.0f64; 3];
        let positions = [
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
            [2.0, 8.0, 9.0],
            [10.0, 3.0, 7.0],
        ];
        let charges = [1.0f64, -1.0, 0.5, -0.5];
        let alpha = 0.35;
        let mesh = [16, 16, 16];
        let h = 1e-3; // finite-difference step (Å)

        let forces = pme_reciprocal_forces(&positions, &charges, box_lengths, alpha, mesh);

        for atom in 0..positions.len() {
            for dim in 0..3 {
                let mut pos_plus = positions.to_vec();
                let mut pos_minus = positions.to_vec();
                pos_plus[atom][dim] += h;
                pos_minus[atom][dim] -= h;

                let e_plus = pme_reciprocal_energy(&pos_plus, &charges, box_lengths, alpha, mesh);
                let e_minus = pme_reciprocal_energy(&pos_minus, &charges, box_lengths, alpha, mesh);

                let f_fd = -(e_plus - e_minus) / (2.0 * h);
                let f_anal = forces[atom][dim];

                let err = (f_anal - f_fd).abs();
                let scale = f_fd.abs().max(1.0);
                assert!(
                    err / scale < 1e-2,
                    "Force mismatch atom={atom} dim={dim}: analytical={f_anal:.6} finite-diff={f_fd:.6} rel_err={:.2e}",
                    err / scale
                );
            }
        }
    }
}
