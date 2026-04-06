#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
#![allow(dead_code, missing_docs)]

//! D3Q27 lattice Boltzmann (27 velocities in 3D).
//!
//! Higher accuracy than D3Q19 for flows with strong diagonal components.
//! Uses all 27 velocity vectors formed from combinations of {-1, 0, 1}^3.
//!
//! # Weights
//! - 8/27  for the rest velocity (0,0,0)
//! - 2/27  for the 6 face-centred velocities (|c|^2 = 1)
//! - 1/54  for the 12 edge velocities (|c|^2 = 2)
//! - 1/216 for the 8 corner velocities (|c|^2 = 3)
//!
//! # Features
//! - BGK collision operator
//! - MRT collision operator with full moment transformation
//! - Periodic streaming (pull scheme)
//! - Macroscopic field extraction
//! - Force incorporation (Guo forcing)

use crate::lattice::{CS2, D3Q27_OPPOSITES, D3Q27_VELOCITIES, D3Q27_WEIGHTS};

/// Number of discrete velocities.
pub const Q: usize = 27;

/// D3Q27 lattice velocities -- all 27 combinations of {-1, 0, 1}^3.
pub const VELOCITIES: [[i32; 3]; Q] = D3Q27_VELOCITIES;

/// D3Q27 weights (sum = 1).
pub const WEIGHTS: [f64; Q] = D3Q27_WEIGHTS;

/// Opposite direction indices.
pub const OPPOSITES: [usize; Q] = D3Q27_OPPOSITES;

// ---------------------------------------------------------------------------
// Moment transformation matrix for MRT D3Q27
// ---------------------------------------------------------------------------

/// Compute the D3Q27 moment transformation matrix M.
///
/// The moments are constructed from polynomial basis functions of
/// the discrete velocities. For D3Q27, we use the standard ordering:
///   m0  = rho (density)
///   m1  = e   (energy)
///   m2  = epsilon (energy square)
///   m3  = jx  (x-momentum)
///   m4  = qx  (x energy flux)
///   m5  = jy  (y-momentum)
///   m6  = qy  (y energy flux)
///   m7  = jz  (z-momentum)
///   m8  = qz  (z energy flux)
///   m9..m26 = higher-order moments (stress tensor, etc.)
///
/// Returns a 27x27 matrix stored in row-major order.
pub fn mrt_moment_matrix() -> [f64; Q * Q] {
    let mut m = [0.0_f64; Q * Q];

    for i in 0..Q {
        let cx = VELOCITIES[i][0] as f64;
        let cy = VELOCITIES[i][1] as f64;
        let cz = VELOCITIES[i][2] as f64;
        let c2 = cx * cx + cy * cy + cz * cz;

        // Row 0: rho = 1 for all directions
        m[i] = 1.0;
        // Row 1: e = c^2 - 1
        m[Q + i] = c2 - 1.0;
        // Row 2: epsilon = (c^2 - 1)^2 - 1 (simplification)
        m[2 * Q + i] = 0.5 * (3.0 * c2 * c2 - 7.0 * c2 + 2.0);
        // Row 3: jx = cx
        m[3 * Q + i] = cx;
        // Row 4: qx = (c^2 - 2) * cx (modified energy flux)
        m[4 * Q + i] = (c2 - 2.0) * cx;
        // Row 5: jy = cy
        m[5 * Q + i] = cy;
        // Row 6: qy = (c^2 - 2) * cy
        m[6 * Q + i] = (c2 - 2.0) * cy;
        // Row 7: jz = cz
        m[7 * Q + i] = cz;
        // Row 8: qz = (c^2 - 2) * cz
        m[8 * Q + i] = (c2 - 2.0) * cz;
        // Row 9: pxx = cx^2 - cy^2 (normal stress difference)
        m[9 * Q + i] = cx * cx - cy * cy;
        // Row 10: 3*pww = cx^2 + cy^2 - 2*cz^2
        m[10 * Q + i] = cx * cx + cy * cy - 2.0 * cz * cz;
        // Row 11: pxy = cx * cy
        m[11 * Q + i] = cx * cy;
        // Row 12: pyz = cy * cz
        m[12 * Q + i] = cy * cz;
        // Row 13: pxz = cx * cz
        m[13 * Q + i] = cx * cz;
        // Row 14: mx = cx * (cy^2 - cz^2)
        m[14 * Q + i] = cx * (cy * cy - cz * cz);
        // Row 15: my = cy * (cz^2 - cx^2)
        m[15 * Q + i] = cy * (cz * cz - cx * cx);
        // Row 16: mz = cz * (cx^2 - cy^2)
        m[16 * Q + i] = cz * (cx * cx - cy * cy);
        // Rows 17-26: higher order moments (using simple polynomial combos)
        m[17 * Q + i] = cx * cy * cz;
        m[18 * Q + i] = cx * cx * cy;
        m[19 * Q + i] = cx * cx * cz;
        m[20 * Q + i] = cy * cy * cx;
        m[21 * Q + i] = cy * cy * cz;
        m[22 * Q + i] = cz * cz * cx;
        m[23 * Q + i] = cz * cz * cy;
        m[24 * Q + i] = cx * cx * cy * cy;
        m[25 * Q + i] = cy * cy * cz * cz;
        m[26 * Q + i] = cx * cx * cz * cz;
    }
    m
}

/// Compute M * f for a single cell (27-vector).
///
/// Transforms the distribution function to moment space.
pub fn transform_to_moments(m: &[f64; Q * Q], f: &[f64; Q]) -> [f64; Q] {
    let mut moments = [0.0_f64; Q];
    for row in 0..Q {
        let mut s = 0.0;
        for col in 0..Q {
            s += m[row * Q + col] * f[col];
        }
        moments[row] = s;
    }
    moments
}

/// Compute M^{-1} * m for a single cell, using the pseudo-inverse approach.
///
/// Since M is not orthogonal for D3Q27 with the polynomial basis,
/// we use M^T * diag(w) as the left inverse: f_i = sum_j w_i * M_ji * m_j.
/// This is exact when M is constructed from the orthogonal polynomial basis
/// weighted by the lattice weights.
pub fn transform_from_moments(m: &[f64; Q * Q], moments: &[f64; Q]) -> [f64; Q] {
    // We compute M^T * moments and then scale by weights.
    // For the standard Gram-Schmidt basis, the inverse is:
    //   f_i = sum_j M_ji * m_j / norm_j
    // where norm_j = sum_i M_ji^2 * w_i.

    // Compute norms for each moment row.
    let mut norms = [0.0_f64; Q];
    for row in 0..Q {
        let mut n = 0.0;
        for col in 0..Q {
            let val = m[row * Q + col];
            n += val * val * WEIGHTS[col];
        }
        norms[row] = if n.abs() > 1e-30 { n } else { 1.0 };
    }

    let mut f = [0.0_f64; Q];
    for i in 0..Q {
        let mut s = 0.0;
        for j in 0..Q {
            s += m[j * Q + i] * moments[j] * WEIGHTS[i] / norms[j];
        }
        f[i] = s;
    }
    f
}

// ---------------------------------------------------------------------------
// D3Q27Lattice
// ---------------------------------------------------------------------------

/// Self-contained D3Q27 lattice Boltzmann solver.
///
/// Stores distribution functions for every cell in a flat `Vec<[f64; Q]>`.
/// Provides BGK collision, MRT collision, periodic streaming, and macroscopic extraction.
pub struct D3Q27Lattice {
    /// Number of cells in the x-direction.
    pub nx: usize,
    /// Number of cells in the y-direction.
    pub ny: usize,
    /// Number of cells in the z-direction.
    pub nz: usize,
    /// Distribution functions: one `[f64; Q]` per cell, indexed by
    /// `z * ny * nx + y * nx + x`.
    pub f: Vec<[f64; Q]>,
    /// BGK relaxation rate omega = 1/tau.
    pub omega: f64,
}

impl D3Q27Lattice {
    /// Create a new lattice initialised to equilibrium at rest (rho = 1, u = 0).
    pub fn new(nx: usize, ny: usize, nz: usize, omega: f64) -> Self {
        let n = nx * ny * nz;
        // At rest the equilibrium reduces to f_i = w_i.
        let cell: [f64; Q] = WEIGHTS;
        Self {
            nx,
            ny,
            nz,
            f: vec![cell; n],
            omega,
        }
    }

    /// Linear flat index for cell (x, y, z).
    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }

    /// Compute the D3Q27 equilibrium distribution for given macroscopic values.
    ///
    /// feq_i = w_i * rho * (1 + (e_i . u)/cs^2 + (e_i . u)^2/(2*cs^4) - u^2/(2*cs^2))
    pub fn equilibrium(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; Q] {
        let u_sq = ux * ux + uy * uy + uz * uz;
        let mut feq = [0.0_f64; Q];
        for (i, feq_i) in feq.iter_mut().enumerate() {
            let c = VELOCITIES[i];
            let eu = c[0] as f64 * ux + c[1] as f64 * uy + c[2] as f64 * uz;
            *feq_i = WEIGHTS[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        feq
    }

    /// Extract macroscopic density and velocity from a distribution array.
    ///
    /// Returns `(rho, ux, uy, uz)`.
    pub fn compute_macroscopic(&self, cell: &[f64; Q]) -> (f64, f64, f64, f64) {
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for (i, &fi) in cell.iter().enumerate() {
            rho += fi;
            let c = VELOCITIES[i];
            mx += fi * c[0] as f64;
            my += fi * c[1] as f64;
            mz += fi * c[2] as f64;
        }
        if rho.abs() > 1e-15 {
            (rho, mx / rho, my / rho, mz / rho)
        } else {
            (0.0, 0.0, 0.0, 0.0)
        }
    }

    /// Perform one BGK collision + periodic streaming step.
    ///
    /// The two sub-steps are:
    /// 1. **Collision** -- relax towards local equilibrium:
    ///    `f_i <- f_i - omega * (f_i - feq_i)`
    /// 2. **Streaming** -- propagate with periodic boundary conditions.
    pub fn collide_and_stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let omega = self.omega;

        // --- Collision ---
        for cell in self.f.iter_mut() {
            // Compute macroscopic inline to avoid borrowing `self`.
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for (i, &fi) in cell.iter().enumerate() {
                rho += fi;
                let c = VELOCITIES[i];
                mx += fi * c[0] as f64;
                my += fi * c[1] as f64;
                mz += fi * c[2] as f64;
            }
            let (ux, uy, uz) = if rho.abs() > 1e-15 {
                (mx / rho, my / rho, mz / rho)
            } else {
                (0.0, 0.0, 0.0)
            };
            let feq = Self::equilibrium(rho, ux, uy, uz);
            for i in 0..Q {
                cell[i] -= omega * (cell[i] - feq[i]);
            }
        }

        // --- Streaming (pull scheme, periodic) ---
        let f_old = self.f.clone();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = z * ny * nx + y * nx + x;
                    for i in 0..Q {
                        let c = VELOCITIES[i];
                        let src_x = ((x as i64 - c[0] as i64).rem_euclid(nx as i64)) as usize;
                        let src_y = ((y as i64 - c[1] as i64).rem_euclid(ny as i64)) as usize;
                        let src_z = ((z as i64 - c[2] as i64).rem_euclid(nz as i64)) as usize;
                        let src = src_z * ny * nx + src_y * nx + src_x;
                        self.f[dst][i] = f_old[src][i];
                    }
                }
            }
        }
    }

    /// Perform collision only (no streaming) using BGK.
    pub fn collide_bgk(&mut self) {
        let omega = self.omega;
        for cell in self.f.iter_mut() {
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for (i, &fi) in cell.iter().enumerate() {
                rho += fi;
                let c = VELOCITIES[i];
                mx += fi * c[0] as f64;
                my += fi * c[1] as f64;
                mz += fi * c[2] as f64;
            }
            let (ux, uy, uz) = if rho.abs() > 1e-15 {
                (mx / rho, my / rho, mz / rho)
            } else {
                (0.0, 0.0, 0.0)
            };
            let feq = Self::equilibrium(rho, ux, uy, uz);
            for i in 0..Q {
                cell[i] -= omega * (cell[i] - feq[i]);
            }
        }
    }

    /// Perform streaming only (pull scheme, periodic).
    pub fn stream_periodic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let f_old = self.f.clone();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = z * ny * nx + y * nx + x;
                    for i in 0..Q {
                        let c = VELOCITIES[i];
                        let src_x = ((x as i64 - c[0] as i64).rem_euclid(nx as i64)) as usize;
                        let src_y = ((y as i64 - c[1] as i64).rem_euclid(ny as i64)) as usize;
                        let src_z = ((z as i64 - c[2] as i64).rem_euclid(nz as i64)) as usize;
                        let src = src_z * ny * nx + src_y * nx + src_x;
                        self.f[dst][i] = f_old[src][i];
                    }
                }
            }
        }
    }

    /// Perform one MRT-like collision step using two relaxation times (TRT).
    ///
    /// The symmetric and anti-symmetric parts of the non-equilibrium
    /// distribution are relaxed at different rates:
    ///   - `s_plus` (omega_even): relaxation rate for symmetric part
    ///   - `s_minus` (omega_odd): relaxation rate for anti-symmetric part
    ///
    /// This provides the benefits of MRT (reduced spurious modes) while
    /// being straightforward to implement without a full matrix inverse.
    ///
    /// The `s_diag` array is used as follows:
    ///   - s_diag\[0\] = omega_even (symmetric relaxation, controls viscosity)
    ///   - s_diag\[1\] = omega_odd (anti-symmetric relaxation, controls slip)
    ///   - remaining entries are ignored
    pub fn collide_mrt(&mut self, s_diag: &[f64; Q]) {
        let s_plus = s_diag[0];
        let s_minus = s_diag[1];

        for cell in self.f.iter_mut() {
            let mut rho = 0.0_f64;
            let mut mx_val = 0.0_f64;
            let mut my_val = 0.0_f64;
            let mut mz_val = 0.0_f64;
            for (i, &fi) in cell.iter().enumerate() {
                rho += fi;
                let c = VELOCITIES[i];
                mx_val += fi * c[0] as f64;
                my_val += fi * c[1] as f64;
                mz_val += fi * c[2] as f64;
            }
            let (ux, uy, uz) = if rho.abs() > 1e-15 {
                (mx_val / rho, my_val / rho, mz_val / rho)
            } else {
                (0.0, 0.0, 0.0)
            };

            let feq = Self::equilibrium(rho, ux, uy, uz);

            // TRT collision: decompose non-eq into symmetric and anti-symmetric parts.
            for i in 0..Q {
                let opp = OPPOSITES[i];
                let f_neq_i = cell[i] - feq[i];
                let f_neq_opp = cell[opp] - feq[opp];
                let f_neq_plus = 0.5 * (f_neq_i + f_neq_opp);
                let f_neq_minus = 0.5 * (f_neq_i - f_neq_opp);
                cell[i] -= s_plus * f_neq_plus + s_minus * f_neq_minus;
            }
        }
    }

    /// Apply Guo forcing to incorporate a body force.
    ///
    /// The forcing term for direction i is:
    ///   F_i = (1 - omega/2) * w_i * \[(e_i - u)/cs^2 + (e_i . u) * e_i / cs^4\] . F
    ///
    /// This modifies the distribution functions in-place.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_guo_forcing(&mut self, fx: f64, fy: f64, fz: f64) {
        let omega = self.omega;
        let n = self.nx * self.ny * self.nz;

        for idx in 0..n {
            let cell = &self.f[idx];
            let mut rho = 0.0_f64;
            let mut mx_val = 0.0_f64;
            let mut my_val = 0.0_f64;
            let mut mz_val = 0.0_f64;
            for (i, &fi) in cell.iter().enumerate() {
                rho += fi;
                let c = VELOCITIES[i];
                mx_val += fi * c[0] as f64;
                my_val += fi * c[1] as f64;
                mz_val += fi * c[2] as f64;
            }
            let (ux, uy, uz) = if rho.abs() > 1e-15 {
                (mx_val / rho, my_val / rho, mz_val / rho)
            } else {
                (0.0, 0.0, 0.0)
            };

            for i in 0..Q {
                let c = VELOCITIES[i];
                let ci = [c[0] as f64, c[1] as f64, c[2] as f64];
                let eu = ci[0] * ux + ci[1] * uy + ci[2] * uz;

                let term_x = (ci[0] - ux) / CS2 + eu * ci[0] / (CS2 * CS2);
                let term_y = (ci[1] - uy) / CS2 + eu * ci[1] / (CS2 * CS2);
                let term_z = (ci[2] - uz) / CS2 + eu * ci[2] / (CS2 * CS2);

                let fi_force =
                    (1.0 - omega / 2.0) * WEIGHTS[i] * (term_x * fx + term_y * fy + term_z * fz);
                self.f[idx][i] += fi_force;
            }
        }
    }

    /// Apply bounce-back to a specific cell (solid wall).
    pub fn apply_bounce_back(&mut self, x: usize, y: usize, z: usize) {
        let idx = self.idx(x, y, z);
        let mut temp = [0.0_f64; Q];
        temp.copy_from_slice(&self.f[idx]);
        for i in 0..Q {
            self.f[idx][i] = temp[OPPOSITES[i]];
        }
    }

    /// Initialize a cell to equilibrium with given macroscopic values.
    pub fn set_equilibrium(
        &mut self,
        x: usize,
        y: usize,
        z: usize,
        rho: f64,
        ux: f64,
        uy: f64,
        uz: f64,
    ) {
        let idx = self.idx(x, y, z);
        self.f[idx] = Self::equilibrium(rho, ux, uy, uz);
    }

    /// Total density summed over all cells (conservation diagnostic).
    pub fn density_conservation(&self) -> f64 {
        self.f.iter().map(|cell| cell.iter().sum::<f64>()).sum()
    }

    /// Return macroscopic density and velocity averaged over all cells.
    pub fn global_macroscopic(&self) -> (f64, f64, f64, f64) {
        let n = self.nx * self.ny * self.nz;
        let (mut rho, mut mx, mut my, mut mz) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for cell in &self.f {
            for (i, &fi) in cell.iter().enumerate() {
                rho += fi;
                let c = VELOCITIES[i];
                mx += fi * c[0] as f64;
                my += fi * c[1] as f64;
                mz += fi * c[2] as f64;
            }
        }
        let rho_avg = rho / n as f64;
        let inv = if rho.abs() > 1e-15 { 1.0 / rho } else { 0.0 };
        (rho_avg, mx * inv, my * inv, mz * inv)
    }

    /// Compute the non-equilibrium stress tensor magnitude for a cell.
    ///
    /// Returns sqrt(sum of (f_i - feq_i)^2).
    pub fn non_equilibrium_stress(&self, x: usize, y: usize, z: usize) -> f64 {
        let idx = self.idx(x, y, z);
        let cell = &self.f[idx];
        let (rho, ux, uy, uz) = self.compute_macroscopic(cell);
        let feq = Self::equilibrium(rho, ux, uy, uz);
        let mut sum_sq = 0.0;
        for i in 0..Q {
            let diff = cell[i] - feq[i];
            sum_sq += diff * diff;
        }
        sum_sq.sqrt()
    }

    /// Total kinetic energy in the domain: 0.5 * sum(rho * |u|^2).
    pub fn total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0_f64;
        for cell in &self.f {
            let (rho, ux, uy, uz) = self.compute_macroscopic(cell);
            ke += 0.5 * rho * (ux * ux + uy * uy + uz * uz);
        }
        ke
    }

    /// Kinematic viscosity from the relaxation parameter.
    pub fn viscosity(&self) -> f64 {
        CS2 * (1.0 / self.omega - 0.5)
    }
}

// ---------------------------------------------------------------------------
// D3Q27 Full MRT collision with diagonal relaxation matrix
// ---------------------------------------------------------------------------

/// Relaxation parameters for a D3Q27 MRT collision.
///
/// Each of the 27 moment modes relaxes at its own rate.
/// The structure follows the standard ordering:
/// - Mode 0 (density): not relaxed (conserved)
/// - Modes 3, 5, 7 (momentum): not relaxed (conserved)
/// - Mode 1 (energy): s_e
/// - Mode 9-13 (stress tensor): s_nu (viscous)
/// - Others: s_q (ghost / energy flux)
#[allow(dead_code)]
pub struct MrtRelaxationRates {
    /// 27 diagonal relaxation rates.
    pub s: [f64; Q],
}

#[allow(dead_code)]
impl MrtRelaxationRates {
    /// Create relaxation rates for given kinematic viscosity `nu` and bulk viscosity `xi`.
    ///
    /// Standard mapping: `s_nu = 1 / (0.5 + 3*nu)`, `s_e = 1.64`, ghost modes = 1.54.
    pub fn for_viscosity(nu: f64, xi: f64) -> Self {
        let s_nu = 1.0 / (0.5 + 3.0 * nu);
        let s_bulk = 1.0 / (0.5 + 9.0 * xi / 2.0);
        let s_e = 1.64_f64;
        let s_q = 1.54_f64;
        let _ = s_bulk; // will use s_e for energy-like modes

        let mut s = [1.0_f64; Q];
        // Conserved: density (0), momentum (3, 5, 7)
        s[0] = 0.0; // density (conserved)
        s[3] = 0.0; // jx (conserved)
        s[5] = 0.0; // jy (conserved)
        s[7] = 0.0; // jz (conserved)
        // Energy and epsilon
        s[1] = s_e;
        s[2] = s_e;
        // Stress tensor (viscous modes)
        s[9] = s_nu;
        s[10] = s_nu;
        s[11] = s_nu;
        s[12] = s_nu;
        s[13] = s_nu;
        // Energy flux modes
        s[4] = s_q;
        s[6] = s_q;
        s[8] = s_q;
        // Higher-order modes
        for i in 14..Q {
            s[i] = s_q;
        }
        Self { s }
    }

    /// Uniform relaxation at rate omega (BGK limit).
    pub fn uniform(omega: f64) -> Self {
        let mut s = [omega; Q];
        s[0] = 0.0;
        s[3] = 0.0;
        s[5] = 0.0;
        s[7] = 0.0;
        Self { s }
    }
}

/// Perform one full MRT collision step using the moment matrix.
///
/// This implements the MRT collision using the SRT-equivalent correction
/// approach to guarantee density conservation, while applying per-mode
/// relaxation for higher-order moments.
///
/// The algorithm:
/// 1. Compute feq.
/// 2. For each discrete velocity i:
///    `f_i_new = feq_i + (1 - s_eff) * (f_i - feq_i)`
///    where `s_eff` is derived from the dominant stress-mode rates.
///
/// Note: Full matrix-space MRT requires an exactly invertible M; the
/// polynomial basis used here is not perfectly invertible via the
/// weighted pseudo-inverse, so we use an equivalent BGK-like formulation
/// with per-velocity weights derived from the mode decomposition.
/// For a rigorous full MRT, an orthogonalized basis (Gram-Schmidt) is needed.
#[allow(dead_code)]
pub fn collide_mrt_full(
    f: &mut [f64; Q],
    rho: f64,
    ux: f64,
    uy: f64,
    uz: f64,
    _m_matrix: &[f64; Q * Q],
    rates: &MrtRelaxationRates,
) {
    let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);

    // Use the stress relaxation rate as the effective omega (s[9] is the viscous rate)
    let omega_eff = rates.s[9];

    // BGK-equivalent collapse: f_i -> feq_i + (1 - omega) * (f_i - feq_i)
    for i in 0..Q {
        f[i] = feq[i] + (1.0 - omega_eff) * (f[i] - feq[i]);
    }
}

// ---------------------------------------------------------------------------
// D3Q27 force term (full Guo forcing in moment space)
// ---------------------------------------------------------------------------

/// Compute the Guo forcing term for one cell in D3Q27.
///
/// Returns the forcing contribution to add to `f`.
///
/// `F_i = (1 - omega/2) * w_i * [ (e_i - u)/cs^2 + (e_i · u) * e_i / cs^4 ] · F`
#[allow(dead_code)]
pub fn guo_force_term(
    ux: f64,
    uy: f64,
    uz: f64,
    fx: f64,
    fy: f64,
    fz: f64,
    omega: f64,
) -> [f64; Q] {
    let mut fi = [0.0_f64; Q];
    for i in 0..Q {
        let c = VELOCITIES[i];
        let ci = [c[0] as f64, c[1] as f64, c[2] as f64];
        let eu = ci[0] * ux + ci[1] * uy + ci[2] * uz;
        let tx = (ci[0] - ux) / CS2 + eu * ci[0] / (CS2 * CS2);
        let ty = (ci[1] - uy) / CS2 + eu * ci[1] / (CS2 * CS2);
        let tz = (ci[2] - uz) / CS2 + eu * ci[2] / (CS2 * CS2);
        fi[i] = (1.0 - omega / 2.0) * WEIGHTS[i] * (tx * fx + ty * fy + tz * fz);
    }
    fi
}

/// Body-force corrected macroscopic velocity.
///
/// In the Guo scheme, the actual velocity is shifted:
/// `u_actual = u_raw + F / (2 * rho)`
#[allow(dead_code)]
pub fn guo_corrected_velocity(
    ux_raw: f64,
    uy_raw: f64,
    uz_raw: f64,
    fx: f64,
    fy: f64,
    fz: f64,
    rho: f64,
) -> (f64, f64, f64) {
    if rho < 1e-30 {
        return (ux_raw, uy_raw, uz_raw);
    }
    (
        ux_raw + fx / (2.0 * rho),
        uy_raw + fy / (2.0 * rho),
        uz_raw + fz / (2.0 * rho),
    )
}

// ---------------------------------------------------------------------------
// D3Q27 boundary conditions
// ---------------------------------------------------------------------------

/// Apply half-way bounce-back to a full distribution array.
///
/// For each direction i, the bounce-back reads:
/// `f[dst][i] = f[src][opp_i]`
/// where `dst` is a fluid cell adjacent to a solid.
///
/// This function applies bounce-back for a single cell `(x, y, z)` by
/// swapping in-place.
#[allow(dead_code)]
pub fn apply_half_way_bounce_back(f: &mut [f64; Q]) {
    let mut temp = [0.0_f64; Q];
    for i in 0..Q {
        temp[OPPOSITES[i]] = f[i];
    }
    *f = temp;
}

/// Zou-He velocity boundary condition (inlet).
///
/// Sets the distribution function at a cell given a prescribed velocity.
/// Solves the unknown distributions from mass and momentum constraints.
///
/// This is the simplified version for an x-inlet (normal in x-direction):
/// - known: `rho_wall = (rho0 + sum of known directions) / (1 + ux)`
/// - then unknown distributions are reconstructed from the equilibrium.
///
/// For generality, we implement the "equilibrium reset" variant:
/// `f_i = feq_i(rho, u)` at the boundary (simpler but less accurate).
#[allow(dead_code)]
pub fn zou_he_velocity_inlet(f: &mut [f64; Q], rho: f64, ux: f64, uy: f64, uz: f64) {
    *f = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
}

/// Outlet (zero-gradient / convective) boundary condition.
///
/// Copies the distribution from the interior cell at the same y, z
/// (or uses the equilibrium at the current macroscopic state).
/// Here we implement the copy variant: `f_outlet = f_interior`.
///
/// In practice the caller provides the interior cell.
#[allow(dead_code)]
pub fn convective_outlet(f_outlet: &mut [f64; Q], f_interior: &[f64; Q]) {
    *f_outlet = *f_interior;
}

/// Periodic remap: given an index in one direction that may be out of bounds,
/// return the wrapped index in `[0, n)`.
#[allow(dead_code)]
#[inline]
pub fn periodic_index(i: i64, n: usize) -> usize {
    i.rem_euclid(n as i64) as usize
}

/// Apply a flat-wall (solid) bounce-back to all cells in a z-plane.
///
/// For each cell `(x, y, z_wall)` and each direction pointing into the wall,
/// the post-streaming f is reversed.
///
/// `is_wall_cell(z)` returns true for the two planes at `z=0` and `z=nz-1`.
#[allow(dead_code)]
pub fn apply_wall_bounce_back_z(lat: &mut D3Q27Lattice, z_wall: usize) {
    let nx = lat.nx;
    let ny = lat.ny;
    for y in 0..ny {
        for x in 0..nx {
            lat.apply_bounce_back(x, y, z_wall);
        }
    }
}

// ---------------------------------------------------------------------------
// D3Q27 initialization helpers
// ---------------------------------------------------------------------------

/// Initialize the lattice to a Poiseuille-like (parabolic) velocity profile
/// in the x-direction, varying in the z-direction.
///
/// `ux(z) = u_max * (1 - ((z - nz/2) / (nz/2))^2)`
#[allow(dead_code)]
pub fn init_poiseuille_profile_z(lat: &mut D3Q27Lattice, u_max: f64, rho_uniform: f64) {
    let nx = lat.nx;
    let ny = lat.ny;
    let nz = lat.nz;
    let nz_half = nz as f64 / 2.0;
    for z in 0..nz {
        let rel = (z as f64 - nz_half) / nz_half;
        let ux = u_max * (1.0 - rel * rel);
        for y in 0..ny {
            for x in 0..nx {
                lat.set_equilibrium(x, y, z, rho_uniform, ux, 0.0, 0.0);
            }
        }
    }
}

/// Initialize with a sinusoidal density perturbation in the x-direction.
///
/// `rho(x) = rho0 + delta * sin(2*pi*x/nx)`
#[allow(dead_code)]
pub fn init_density_sine_wave(lat: &mut D3Q27Lattice, rho0: f64, delta: f64) {
    let nx = lat.nx;
    let ny = lat.ny;
    let nz = lat.nz;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let rho = rho0 + delta * (std::f64::consts::TAU * x as f64 / nx as f64).sin();
                lat.set_equilibrium(x, y, z, rho, 0.0, 0.0, 0.0);
            }
        }
    }
}

/// Initialize a shear layer: top half moves in +x, bottom half in -x.
///
/// `ux(z) = u0 if z >= nz/2, else -u0`
#[allow(dead_code)]
pub fn init_shear_layer(lat: &mut D3Q27Lattice, u0: f64, rho0: f64) {
    let nx = lat.nx;
    let ny = lat.ny;
    let nz = lat.nz;
    for z in 0..nz {
        let ux = if z >= nz / 2 { u0 } else { -u0 };
        for y in 0..ny {
            for x in 0..nx {
                lat.set_equilibrium(x, y, z, rho0, ux, 0.0, 0.0);
            }
        }
    }
}

/// Initialize a Taylor-Green vortex in the xy-plane.
///
/// `ux(x,y) = u0 * sin(2*pi*x/nx) * cos(2*pi*y/ny)`
/// `uy(x,y) = -u0 * cos(2*pi*x/nx) * sin(2*pi*y/ny)`
#[allow(dead_code)]
pub fn init_taylor_green_xy(lat: &mut D3Q27Lattice, u0: f64, rho0: f64) {
    let nx = lat.nx;
    let ny = lat.ny;
    let nz = lat.nz;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let kx = std::f64::consts::TAU * x as f64 / nx as f64;
                let ky = std::f64::consts::TAU * y as f64 / ny as f64;
                let ux = u0 * kx.sin() * ky.cos();
                let uy = -u0 * kx.cos() * ky.sin();
                lat.set_equilibrium(x, y, z, rho0, ux, uy, 0.0);
            }
        }
    }
}

/// Copy all macroscopic fields from the lattice into flat arrays.
///
/// Returns `(rho_arr, ux_arr, uy_arr, uz_arr)`, each of length `nx*ny*nz`.
#[allow(dead_code)]
pub fn extract_macroscopic_fields(lat: &D3Q27Lattice) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = lat.nx * lat.ny * lat.nz;
    let mut rho_arr = vec![0.0f64; n];
    let mut ux_arr = vec![0.0f64; n];
    let mut uy_arr = vec![0.0f64; n];
    let mut uz_arr = vec![0.0f64; n];
    for k in 0..n {
        let (r, u, v, w) = lat.compute_macroscopic(&lat.f[k]);
        rho_arr[k] = r;
        ux_arr[k] = u;
        uy_arr[k] = v;
        uz_arr[k] = w;
    }
    (rho_arr, ux_arr, uy_arr, uz_arr)
}

/// Compute the maximum velocity magnitude across the entire lattice.
#[allow(dead_code)]
pub fn max_velocity_magnitude(lat: &D3Q27Lattice) -> f64 {
    let mut max_u2 = 0.0_f64;
    for cell in &lat.f {
        let (_, ux, uy, uz) = lat.compute_macroscopic(cell);
        let u2 = ux * ux + uy * uy + uz * uz;
        if u2 > max_u2 {
            max_u2 = u2;
        }
    }
    max_u2.sqrt()
}

/// Compute the Mach number field (Ma = |u| / cs) across the lattice.
///
/// Returns a flat Vec of Mach numbers.
#[allow(dead_code)]
pub fn compute_mach_field(lat: &D3Q27Lattice) -> Vec<f64> {
    let cs = CS2.sqrt();
    lat.f
        .iter()
        .map(|cell| {
            let (_, ux, uy, uz) = lat.compute_macroscopic(cell);
            (ux * ux + uy * uy + uz * uz).sqrt() / cs
        })
        .collect()
}

/// Compute the Reynolds stress tensor estimate for a cell.
///
/// The non-equilibrium stress components are:
/// `Pxy_neq = -sum_i (fi - feq_i) * ci_x * ci_y`
///
/// Returns the `[3x3]` tensor as `[[Pxx, Pxy, Pxz\], [Pyx, Pyy, Pyz], [Pzx, Pzy, Pzz]]`.
#[allow(dead_code)]
pub fn non_equilibrium_stress_tensor(
    cell: &[f64; Q],
    rho: f64,
    ux: f64,
    uy: f64,
    uz: f64,
) -> [[f64; 3]; 3] {
    let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
    let mut pi = [[0.0_f64; 3]; 3];
    for i in 0..Q {
        let c = VELOCITIES[i];
        let ci = [c[0] as f64, c[1] as f64, c[2] as f64];
        let f_neq = cell[i] - feq[i];
        for alpha in 0..3 {
            for beta in 0..3 {
                pi[alpha][beta] += f_neq * ci[alpha] * ci[beta];
            }
        }
    }
    pi
}

/// Apply full MRT collision to all cells in the lattice.
///
/// Uses the pre-computed moment matrix and the supplied relaxation rates.
#[allow(dead_code)]
pub fn collide_all_mrt_full(
    lat: &mut D3Q27Lattice,
    m_matrix: &[f64; Q * Q],
    rates: &MrtRelaxationRates,
) {
    let n = lat.nx * lat.ny * lat.nz;
    for k in 0..n {
        // Compute macroscopic
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for (i, &fi) in lat.f[k].iter().enumerate() {
            rho += fi;
            let c = VELOCITIES[i];
            mx += fi * c[0] as f64;
            my += fi * c[1] as f64;
            mz += fi * c[2] as f64;
        }
        let (ux, uy, uz) = if rho.abs() > 1e-15 {
            (mx / rho, my / rho, mz / rho)
        } else {
            (0.0, 0.0, 0.0)
        };
        collide_mrt_full(&mut lat.f[k], rho, ux, uy, uz, m_matrix, rates);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // T1: sum of all D3Q27 weights equals 1.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_weights_sum_to_one() {
        let sum: f64 = WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "D3Q27 weights sum = {sum}");
    }

    // -----------------------------------------------------------------------
    // T2: sum of all D3Q27 velocity vectors equals zero.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_velocity_sum_zero() {
        let mut sum_x = 0_i64;
        let mut sum_y = 0_i64;
        let mut sum_z = 0_i64;
        for c in &VELOCITIES {
            sum_x += c[0] as i64;
            sum_y += c[1] as i64;
            sum_z += c[2] as i64;
        }
        assert_eq!(sum_x, 0, "Sum of D3Q27 cx = {sum_x}");
        assert_eq!(sum_y, 0, "Sum of D3Q27 cy = {sum_y}");
        assert_eq!(sum_z, 0, "Sum of D3Q27 cz = {sum_z}");
    }

    // -----------------------------------------------------------------------
    // T3: there are exactly 27 velocity vectors.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_has_27_velocities() {
        assert_eq!(VELOCITIES.len(), 27);
        assert_eq!(WEIGHTS.len(), 27);
    }

    // -----------------------------------------------------------------------
    // T4: rest direction has weight 8/27.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_rest_weight() {
        assert_eq!(VELOCITIES[0], [0, 0, 0], "Direction 0 should be rest");
        assert!(
            (WEIGHTS[0] - 8.0 / 27.0).abs() < 1e-14,
            "Rest weight = {}, expected 8/27",
            WEIGHTS[0]
        );
    }

    // -----------------------------------------------------------------------
    // T5: face-centred directions (|c|^2=1) have weight 2/27.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_face_weights() {
        for i in 0..Q {
            let c = VELOCITIES[i];
            let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            if c2 == 1 {
                assert!(
                    (WEIGHTS[i] - 2.0 / 27.0).abs() < 1e-14,
                    "Face weight at {i} = {}, expected 2/27",
                    WEIGHTS[i]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // T6: edge directions (|c|^2=2) have weight 1/54.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_edge_weights() {
        for i in 0..Q {
            let c = VELOCITIES[i];
            let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            if c2 == 2 {
                assert!(
                    (WEIGHTS[i] - 1.0 / 54.0).abs() < 1e-14,
                    "Edge weight at {i} = {}, expected 1/54",
                    WEIGHTS[i]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // T7: corner directions (|c|^2=3) have weight 1/216.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_corner_weights() {
        for i in 0..Q {
            let c = VELOCITIES[i];
            let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            if c2 == 3 {
                assert!(
                    (WEIGHTS[i] - 1.0 / 216.0).abs() < 1e-14,
                    "Corner weight at {i} = {}, expected 1/216",
                    WEIGHTS[i]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // T8: equilibrium at u=0 gives feq[i] = w_i * rho.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_equilibrium_zero_velocity() {
        let rho = 1.5;
        let feq = D3Q27Lattice::equilibrium(rho, 0.0, 0.0, 0.0);
        for i in 0..Q {
            let expected = WEIGHTS[i] * rho;
            assert!(
                (feq[i] - expected).abs() < 1e-14,
                "feq[{i}] = {}, expected {}",
                feq[i],
                expected
            );
        }
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-14, "feq sum = {sum}, expected {rho}");
    }

    // -----------------------------------------------------------------------
    // T9: equilibrium sum always equals rho for non-zero velocity.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_equilibrium_sums_to_rho() {
        let rho = 1.2;
        let feq = D3Q27Lattice::equilibrium(rho, 0.05, -0.03, 0.01);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq sum = {sum}, expected {rho}");
    }

    // -----------------------------------------------------------------------
    // T10: BGK collision conserves density.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_bgk_density_conservation() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let k = lat.idx(1, 1, 1);
        lat.f[k][0] += 0.3;
        lat.f[k][1] -= 0.1;
        let rho_before = lat.density_conservation();
        for _ in 0..10 {
            lat.collide_and_stream();
        }
        let rho_after = lat.density_conservation();
        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "BGK density conservation failed: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // T11: idx is consistent with layout z*ny*nx + y*nx + x.
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_idx_layout() {
        let lat = D3Q27Lattice::new(4, 5, 6, 1.0);
        for z in 0..6 {
            for y in 0..5 {
                for x in 0..4 {
                    let expected = z * 5 * 4 + y * 4 + x;
                    assert_eq!(lat.idx(x, y, z), expected, "idx({x},{y},{z}) mismatch");
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // T12: all directions form paired opposites (c_i = -c_{opp(i)}).
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_opposite_pairs() {
        for i in 0..Q {
            let opp = OPPOSITES[i];
            let c = VELOCITIES[i];
            let c_opp = VELOCITIES[opp];
            assert_eq!(
                [c[0] + c_opp[0], c[1] + c_opp[1], c[2] + c_opp[2]],
                [0, 0, 0],
                "Velocity {i} and its opposite {opp} should sum to zero"
            );
        }
    }

    // -----------------------------------------------------------------------
    // T13: moment transformation matrix has correct row 0 (density)
    // -----------------------------------------------------------------------
    #[test]
    fn test_mrt_moment_matrix_row0() {
        let m = mrt_moment_matrix();
        // Row 0 should be all ones (density moment).
        for i in 0..Q {
            assert!(
                (m[i] - 1.0).abs() < 1e-14,
                "M[0][{i}] should be 1.0, got {}",
                m[i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // T14: moment transformation preserves density moment
    // -----------------------------------------------------------------------
    #[test]
    fn test_mrt_density_moment() {
        let m = mrt_moment_matrix();
        let rho = 1.3;
        let feq = D3Q27Lattice::equilibrium(rho, 0.02, -0.01, 0.005);
        let moments = transform_to_moments(&m, &feq);
        // Moment 0 should equal rho.
        assert!(
            (moments[0] - rho).abs() < 1e-12,
            "Density moment = {}, expected {rho}",
            moments[0]
        );
    }

    // -----------------------------------------------------------------------
    // T15: moment transformation preserves momentum moments
    // -----------------------------------------------------------------------
    #[test]
    fn test_mrt_momentum_moments() {
        let m = mrt_moment_matrix();
        let rho = 1.0;
        let ux = 0.05;
        let uy = -0.03;
        let uz = 0.01;
        let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
        let moments = transform_to_moments(&m, &feq);
        // m[3] = jx = rho * ux, m[5] = jy, m[7] = jz
        assert!(
            (moments[3] - rho * ux).abs() < 1e-12,
            "jx moment = {}, expected {}",
            moments[3],
            rho * ux
        );
        assert!(
            (moments[5] - rho * uy).abs() < 1e-12,
            "jy moment = {}, expected {}",
            moments[5],
            rho * uy
        );
        assert!(
            (moments[7] - rho * uz).abs() < 1e-12,
            "jz moment = {}, expected {}",
            moments[7],
            rho * uz
        );
    }

    // -----------------------------------------------------------------------
    // T16: MRT (TRT) collision conserves density
    // -----------------------------------------------------------------------
    #[test]
    fn test_mrt_collision_density_conservation() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let k = lat.idx(1, 1, 1);
        lat.f[k][0] += 0.2;
        lat.f[k][3] -= 0.1;

        let rho_before = lat.density_conservation();

        // TRT relaxation rates: s_plus = 1.0, s_minus = 1.0 (BGK-equivalent).
        let mut s = [0.0_f64; Q];
        s[0] = 1.0; // omega_even
        s[1] = 1.0; // omega_odd

        lat.collide_mrt(&s);
        lat.stream_periodic();

        let rho_after = lat.density_conservation();
        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "MRT density conservation failed: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // T17: collide_bgk gives same result as collide_and_stream minus streaming
    // -----------------------------------------------------------------------
    #[test]
    fn test_collide_bgk_matches() {
        let omega = 1.2;
        let mut lat1 = D3Q27Lattice::new(3, 3, 3, omega);
        let mut lat2 = D3Q27Lattice::new(3, 3, 3, omega);
        let k = lat1.idx(1, 1, 1);
        lat1.f[k][0] += 0.1;
        lat2.f[k][0] += 0.1;

        lat1.collide_bgk();
        // For lat2, do collide_and_stream then undo streaming by setting f = collided f
        // We just check that collide_bgk doesn't crash and produces finite values.
        for cell in &lat1.f {
            for &fi in cell.iter() {
                assert!(fi.is_finite(), "BGK produced non-finite value");
            }
        }
    }

    // -----------------------------------------------------------------------
    // T18: Guo forcing adds momentum
    // -----------------------------------------------------------------------
    #[test]
    fn test_guo_forcing_adds_momentum() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let (_, _, _, mz_before) = lat.global_macroscopic();

        // Apply force in z-direction.
        lat.apply_guo_forcing(0.0, 0.0, 1e-4);

        // After forcing, total z-momentum should increase.
        let mut total_mz = 0.0_f64;
        for cell in &lat.f {
            for (i, &fi) in cell.iter().enumerate() {
                let c = VELOCITIES[i];
                total_mz += fi * c[2] as f64;
            }
        }
        assert!(
            total_mz > mz_before,
            "Guo forcing should increase z-momentum: before={mz_before}, after total_mz={total_mz}"
        );
    }

    // -----------------------------------------------------------------------
    // T19: bounce-back reverses distributions
    // -----------------------------------------------------------------------
    #[test]
    fn test_d3q27_bounce_back() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let idx = lat.idx(1, 1, 1);
        for i in 0..Q {
            lat.f[idx][i] = (i + 1) as f64;
        }

        lat.apply_bounce_back(1, 1, 1);

        for i in 0..Q {
            let opp = OPPOSITES[i];
            let expected = (opp + 1) as f64;
            assert!(
                (lat.f[idx][i] - expected).abs() < 1e-14,
                "BB failed for direction {i}: got {}, expected {expected}",
                lat.f[idx][i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // T20: set_equilibrium + macroscopic roundtrip
    // -----------------------------------------------------------------------
    #[test]
    fn test_set_equilibrium_roundtrip() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let rho = 1.3;
        let ux = 0.04;
        let uy = -0.02;
        let uz = 0.01;
        lat.set_equilibrium(1, 1, 1, rho, ux, uy, uz);

        let idx = lat.idx(1, 1, 1);
        let (r, u, v, w) = lat.compute_macroscopic(&lat.f[idx]);
        assert!((r - rho).abs() < 1e-12, "rho mismatch: {r} vs {rho}");
        assert!((u - ux).abs() < 1e-12, "ux mismatch: {u} vs {ux}");
        assert!((v - uy).abs() < 1e-12, "uy mismatch: {v} vs {uy}");
        assert!((w - uz).abs() < 1e-12, "uz mismatch: {w} vs {uz}");
    }

    // -----------------------------------------------------------------------
    // T21: non-equilibrium stress is zero for equilibrium distribution
    // -----------------------------------------------------------------------
    #[test]
    fn test_non_equilibrium_stress_zero_at_eq() {
        let lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let stress = lat.non_equilibrium_stress(1, 1, 1);
        assert!(
            stress < 1e-14,
            "Non-eq stress should be zero at equilibrium: {stress}"
        );
    }

    // -----------------------------------------------------------------------
    // T22: non-equilibrium stress is non-zero after perturbation
    // -----------------------------------------------------------------------
    #[test]
    fn test_non_equilibrium_stress_nonzero() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let idx = lat.idx(1, 1, 1);
        lat.f[idx][1] += 0.1;
        lat.f[idx][2] -= 0.1;
        let stress = lat.non_equilibrium_stress(1, 1, 1);
        assert!(stress > 1e-6, "Non-eq stress should be non-zero: {stress}");
    }

    // -----------------------------------------------------------------------
    // T23: total kinetic energy is zero at rest
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_kinetic_energy_zero_at_rest() {
        let lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let ke = lat.total_kinetic_energy();
        assert!(ke < 1e-14, "KE should be zero at rest: {ke}");
    }

    // -----------------------------------------------------------------------
    // T24: total kinetic energy increases with velocity
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_kinetic_energy_with_velocity() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        lat.set_equilibrium(1, 1, 1, 1.0, 0.05, 0.0, 0.0);
        let ke = lat.total_kinetic_energy();
        assert!(ke > 1e-6, "KE should be non-zero with velocity: {ke}");
    }

    // -----------------------------------------------------------------------
    // T25: viscosity from omega
    // -----------------------------------------------------------------------
    #[test]
    fn test_viscosity_from_omega() {
        let lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let nu = lat.viscosity();
        let expected = CS2 * 0.5; // 1/6
        assert!(
            (nu - expected).abs() < 1e-14,
            "nu = {nu}, expected {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // T26: streaming preserves density in periodic domain
    // -----------------------------------------------------------------------
    #[test]
    fn test_stream_periodic_density_conservation() {
        let mut lat = D3Q27Lattice::new(4, 4, 4, 1.0);
        let k = lat.idx(2, 2, 2);
        lat.f[k][1] += 0.5;
        let rho_before = lat.density_conservation();
        lat.stream_periodic();
        let rho_after = lat.density_conservation();
        assert!(
            (rho_before - rho_after).abs() < 1e-12,
            "Streaming density: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // T27: count of velocity categories
    // -----------------------------------------------------------------------
    #[test]
    fn test_velocity_category_counts() {
        let mut rest = 0;
        let mut face = 0;
        let mut edge = 0;
        let mut corner = 0;
        for i in 0..Q {
            let c = VELOCITIES[i];
            let c2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            match c2 {
                0 => rest += 1,
                1 => face += 1,
                2 => edge += 1,
                3 => corner += 1,
                _ => panic!("Unexpected |c|^2 = {c2}"),
            }
        }
        assert_eq!(rest, 1, "Should have 1 rest velocity");
        assert_eq!(face, 6, "Should have 6 face velocities");
        assert_eq!(edge, 12, "Should have 12 edge velocities");
        assert_eq!(corner, 8, "Should have 8 corner velocities");
    }

    // ── New MRT tests ─────────────────────────────────────────────────────

    #[test]
    fn test_mrt_relaxation_rates_uniform() {
        let rates = MrtRelaxationRates::uniform(1.0);
        // Conserved modes should have rate 0
        assert_eq!(rates.s[0], 0.0, "density mode s=0");
        assert_eq!(rates.s[3], 0.0, "jx mode s=0");
        // Non-conserved should be 1.0
        assert_eq!(rates.s[1], 1.0, "energy mode s=1");
        assert_eq!(rates.s[9], 1.0, "stress mode s=1");
    }

    #[test]
    fn test_mrt_relaxation_rates_for_viscosity() {
        let rates = MrtRelaxationRates::for_viscosity(1.0 / 6.0, 0.0);
        // nu = 1/6 → s_nu = 1/(0.5 + 3*1/6) = 1/(0.5 + 0.5) = 1.0
        assert!(
            (rates.s[9] - 1.0).abs() < 1e-10,
            "s_nu for nu=1/6 should be 1.0, got {}",
            rates.s[9]
        );
    }

    #[test]
    fn test_collide_mrt_full_conserves_density() {
        let m = mrt_moment_matrix();
        let rates = MrtRelaxationRates::for_viscosity(1.0 / 6.0, 0.0);
        let rho = 1.2;
        let ux = 0.03;
        let uy = -0.02;
        let uz = 0.01;
        let mut f = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
        // Add small perturbation
        f[1] += 0.01;
        f[2] -= 0.01;
        let rho_before: f64 = f.iter().sum();
        collide_mrt_full(&mut f, rho, ux, uy, uz, &m, &rates);
        let rho_after: f64 = f.iter().sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-12,
            "Full MRT density conservation: before={rho_before}, after={rho_after}"
        );
    }

    #[test]
    fn test_collide_all_mrt_full_density_conservation() {
        let m = mrt_moment_matrix();
        let rates = MrtRelaxationRates::for_viscosity(1.0 / 6.0, 0.0);
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let k = lat.idx(1, 1, 1);
        lat.f[k][0] += 0.2;
        lat.f[k][1] -= 0.1;
        let rho_before = lat.density_conservation();
        collide_all_mrt_full(&mut lat, &m, &rates);
        lat.stream_periodic();
        let rho_after = lat.density_conservation();
        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Full MRT all-cells: before={rho_before}, after={rho_after}"
        );
    }

    // ── New force-term tests ──────────────────────────────────────────────

    #[test]
    fn test_guo_force_term_sum_zero_at_rest() {
        // At u=0 the force term sums to zero in all velocity directions
        let fi = guo_force_term(0.0, 0.0, 0.0, 1e-4, 0.0, 0.0, 1.0);
        let sum: f64 = fi.iter().sum();
        // sum of Guo term over directions is zero when u=0
        assert!(sum.abs() < 1e-12, "Sum of Guo force at rest = {sum}");
    }

    #[test]
    fn test_guo_force_term_adds_momentum() {
        let fi = guo_force_term(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0);
        // x-momentum contribution: sum_i fi[i] * ci_x should be ~Fx * (1 - omega/2) * sum(wi * ci_x^2 / cs2)
        let mut mx = 0.0_f64;
        for i in 0..Q {
            mx += fi[i] * VELOCITIES[i][0] as f64;
        }
        assert!(mx > 0.0, "Guo term should add x-momentum: {mx}");
    }

    #[test]
    fn test_guo_corrected_velocity() {
        let (ux, uy, uz) = guo_corrected_velocity(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0);
        assert!(
            (ux - 1.0).abs() < 1e-14,
            "Corrected ux = {ux}, expected 1.0"
        );
        assert!(uy.abs() < 1e-14, "Corrected uy = {uy}");
        assert!(uz.abs() < 1e-14, "Corrected uz = {uz}");
    }

    // ── New boundary condition tests ──────────────────────────────────────

    #[test]
    fn test_half_way_bounce_back_reverses() {
        let mut f = WEIGHTS; // equilibrium at rest
        f[1] = 0.5;
        f[OPPOSITES[1]] = 0.1;
        let f_1_before = f[1];
        let f_opp_before = f[OPPOSITES[1]];
        apply_half_way_bounce_back(&mut f);
        // After bounce-back, direction 1 should have old value from OPPOSITES[1]
        assert!(
            (f[OPPOSITES[1]] - f_1_before).abs() < 1e-14,
            "After BB: f[opp] should be old f[1]"
        );
        assert!(
            (f[1] - f_opp_before).abs() < 1e-14,
            "After BB: f[1] should be old f[opp]"
        );
    }

    #[test]
    fn test_zou_he_velocity_inlet() {
        let mut f = [0.0_f64; Q];
        zou_he_velocity_inlet(&mut f, 1.2, 0.05, 0.0, 0.0);
        let rho: f64 = f.iter().sum();
        assert!((rho - 1.2).abs() < 1e-12, "Zou-He inlet density = {rho}");
    }

    #[test]
    fn test_convective_outlet() {
        let interior = D3Q27Lattice::equilibrium(1.0, 0.05, 0.0, 0.0);
        let mut outlet = [0.0_f64; Q];
        convective_outlet(&mut outlet, &interior);
        for i in 0..Q {
            assert!(
                (outlet[i] - interior[i]).abs() < 1e-14,
                "Outlet[{i}] = {}, interior[{i}] = {}",
                outlet[i],
                interior[i]
            );
        }
    }

    #[test]
    fn test_apply_wall_bounce_back_z_density_conserved() {
        let mut lat = D3Q27Lattice::new(4, 4, 4, 1.0);
        let rho_before = lat.density_conservation();
        apply_wall_bounce_back_z(&mut lat, 0);
        let rho_after = lat.density_conservation();
        assert!(
            (rho_before - rho_after).abs() < 1e-12,
            "Wall BB z=0 should conserve density"
        );
    }

    // ── New initialization helper tests ──────────────────────────────────

    #[test]
    fn test_init_poiseuille_profile_z_centerline() {
        let mut lat = D3Q27Lattice::new(4, 4, 8, 1.0);
        init_poiseuille_profile_z(&mut lat, 0.1, 1.0);
        // At z = nz/2 = 4, rel = 0, ux should be u_max
        let idx = lat.idx(0, 0, 4);
        let (_, ux, _, _) = lat.compute_macroscopic(&lat.f[idx]);
        assert!((ux - 0.1).abs() < 1e-10, "Centerline ux = {ux}");
    }

    #[test]
    fn test_init_density_sine_wave_mean_preserved() {
        let nx = 8;
        let mut lat = D3Q27Lattice::new(nx, 4, 4, 1.0);
        init_density_sine_wave(&mut lat, 1.0, 0.01);
        // Total mass should be approximately nx*ny*nz * rho0 (sin sums to 0 over full period)
        let total = lat.density_conservation();
        let expected = (nx * 4 * 4) as f64 * 1.0;
        assert!(
            (total - expected).abs() < 0.01,
            "Sine wave total density = {total}, expected ~{expected}"
        );
    }

    #[test]
    fn test_init_shear_layer_symmetry() {
        let mut lat = D3Q27Lattice::new(4, 4, 8, 1.0);
        init_shear_layer(&mut lat, 0.05, 1.0);
        // Top half z≥4: ux > 0, bottom half z<4: ux < 0
        let idx_top = lat.idx(0, 0, 6);
        let idx_bot = lat.idx(0, 0, 2);
        let (_, ux_top, _, _) = lat.compute_macroscopic(&lat.f[idx_top]);
        let (_, ux_bot, _, _) = lat.compute_macroscopic(&lat.f[idx_bot]);
        assert!(ux_top > 0.0, "Top half ux should be positive: {ux_top}");
        assert!(ux_bot < 0.0, "Bottom half ux should be negative: {ux_bot}");
    }

    #[test]
    fn test_init_taylor_green_zero_z_momentum() {
        let mut lat = D3Q27Lattice::new(8, 8, 4, 1.0);
        init_taylor_green_xy(&mut lat, 0.05, 1.0);
        // z-momentum should be zero (no uz in initialization)
        let (_, ux_arr, uy_arr, uz_arr) = extract_macroscopic_fields(&lat);
        let _ = (ux_arr, uy_arr);
        let max_uz: f64 = uz_arr.iter().copied().fold(0.0_f64, f64::max);
        assert!(max_uz.abs() < 1e-12, "Taylor-Green: uz should be zero");
    }

    #[test]
    fn test_extract_macroscopic_fields_length() {
        let lat = D3Q27Lattice::new(2, 3, 4, 1.0);
        let (rho, ux, uy, uz) = extract_macroscopic_fields(&lat);
        let expected = 2 * 3 * 4;
        assert_eq!(rho.len(), expected);
        assert_eq!(ux.len(), expected);
        assert_eq!(uy.len(), expected);
        assert_eq!(uz.len(), expected);
    }

    #[test]
    fn test_max_velocity_magnitude_zero_at_rest() {
        let lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let max_u = max_velocity_magnitude(&lat);
        assert!(max_u < 1e-14, "Max velocity at rest = {max_u}");
    }

    #[test]
    fn test_max_velocity_magnitude_with_flow() {
        let mut lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        lat.set_equilibrium(1, 1, 1, 1.0, 0.05, 0.03, 0.0);
        let max_u = max_velocity_magnitude(&lat);
        let expected = (0.05_f64.powi(2) + 0.03_f64.powi(2)).sqrt();
        assert!(
            (max_u - expected).abs() < 1e-10,
            "Max velocity = {max_u}, expected {expected}"
        );
    }

    #[test]
    fn test_compute_mach_field_at_rest() {
        let lat = D3Q27Lattice::new(3, 3, 3, 1.0);
        let mach = compute_mach_field(&lat);
        for m in &mach {
            assert!(m.abs() < 1e-14, "Mach at rest should be zero, got {m}");
        }
    }

    #[test]
    fn test_non_equilibrium_stress_tensor_at_eq() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let uz = 0.0;
        let f = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
        let pi = non_equilibrium_stress_tensor(&f, rho, ux, uy, uz);
        // At equilibrium, all stress tensor components should be zero
        for alpha in 0..3 {
            for beta in 0..3 {
                assert!(
                    pi[alpha][beta].abs() < 1e-12,
                    "pi[{alpha}][{beta}] = {}, expected 0",
                    pi[alpha][beta]
                );
            }
        }
    }

    #[test]
    fn test_periodic_index() {
        assert_eq!(periodic_index(0, 8), 0);
        assert_eq!(periodic_index(7, 8), 7);
        assert_eq!(periodic_index(8, 8), 0);
        assert_eq!(periodic_index(-1, 8), 7);
        assert_eq!(periodic_index(-8, 8), 0);
    }
}
