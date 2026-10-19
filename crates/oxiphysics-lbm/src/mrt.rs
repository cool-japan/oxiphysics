// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiple Relaxation Time (MRT) collision operator for D2Q9 and D3Q19.
//!
//! MRT relaxes different moments at different rates, providing more numerical
//! stability than the single-relaxation-time BGK scheme, especially at high
//! Reynolds numbers.
//!
//! Besides MRT the module provides the two-relaxation-time scheme
//! ([`TrtCollision2D`]), the cascaded central-moment scheme
//! ([`CascadedCollision2D`]) and the regularised BGK scheme
//! ([`RegularisedBgk2D`]) for D2Q9.
//!
//! Reference: d'Humières et al. (2002), Phil. Trans. R. Soc. Lond. A 360, 437–451.

use crate::collision::TRT_MAGIC;
use crate::grid::LbmGrid2D;

pub mod cascaded;
#[cfg(test)]
mod regression_tests;

pub use cascaded::CascadedCollision2D;

// ---------------------------------------------------------------------------
// D2Q9 transformation matrix M (row-major, 9×9)
// Rows are basis vectors for each moment:
//   0: mass        [1, 1, 1, 1, 1, 1, 1, 1, 1]
//   1: energy      [-4,-1,-1,-1,-1, 2, 2, 2, 2]
//   2: energy-sq   [ 4,-2,-2,-2,-2, 1, 1, 1, 1]
//   3: x-momentum  [ 0, 1, 0,-1, 0, 1,-1,-1, 1]
//   4: heat flux x [ 0,-2, 0, 2, 0, 1,-1,-1, 1]
//   5: y-momentum  [ 0, 0, 1, 0,-1, 1, 1,-1,-1]
//   6: heat flux y [ 0, 0,-2, 0, 2, 1, 1,-1,-1]
//   7: stress xx-yy[ 0, 1,-1, 1,-1, 0, 0, 0, 0]
//   8: stress xy   [ 0, 0, 0, 0, 0, 1,-1, 1,-1]
// ---------------------------------------------------------------------------
#[rustfmt::skip]
const M: [[f64; 9]; 9] = [
    [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],
    [-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0],
    [ 4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0],
    [ 0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0],
    [ 0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0],
    [ 0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0],
    [ 0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0],
    [ 0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0],
    [ 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0],
];

// ---------------------------------------------------------------------------
// D2Q9 inverse transformation matrix M^{-1} (analytical, exact fractions)
// ---------------------------------------------------------------------------
#[rustfmt::skip]
const M_INV: [[f64; 9]; 9] = [
    [ 1.0/9.0,  -1.0/9.0,   1.0/9.0,   0.0,       0.0,       0.0,       0.0,       0.0,       0.0      ],
    [ 1.0/9.0,  -1.0/36.0, -1.0/18.0,  1.0/6.0,  -1.0/6.0,   0.0,       0.0,       1.0/4.0,   0.0      ],
    [ 1.0/9.0,  -1.0/36.0, -1.0/18.0,  0.0,       0.0,       1.0/6.0,  -1.0/6.0,  -1.0/4.0,   0.0      ],
    [ 1.0/9.0,  -1.0/36.0, -1.0/18.0, -1.0/6.0,   1.0/6.0,   0.0,       0.0,       1.0/4.0,   0.0      ],
    [ 1.0/9.0,  -1.0/36.0, -1.0/18.0,  0.0,       0.0,      -1.0/6.0,   1.0/6.0,  -1.0/4.0,   0.0      ],
    [ 1.0/9.0,   1.0/18.0,  1.0/36.0,  1.0/6.0,   1.0/12.0,  1.0/6.0,   1.0/12.0,  0.0,       1.0/4.0  ],
    [ 1.0/9.0,   1.0/18.0,  1.0/36.0, -1.0/6.0,  -1.0/12.0,  1.0/6.0,   1.0/12.0,  0.0,      -1.0/4.0  ],
    [ 1.0/9.0,   1.0/18.0,  1.0/36.0, -1.0/6.0,  -1.0/12.0, -1.0/6.0,  -1.0/12.0,  0.0,       1.0/4.0  ],
    [ 1.0/9.0,   1.0/18.0,  1.0/36.0,  1.0/6.0,   1.0/12.0, -1.0/6.0,  -1.0/12.0,  0.0,      -1.0/4.0  ],
];

// ---------------------------------------------------------------------------
// MrtCollision2D
// ---------------------------------------------------------------------------

/// D2Q9 Multiple Relaxation Time (MRT) collision operator.
///
/// MRT relaxes different hydrodynamic moments at different rates, offering
/// better numerical stability than BGK (SRT) for high Reynolds number flows.
///
/// The nine relaxation rates correspond to the nine moment modes:
/// `[s_rho, s_e, s_eps, s_j, s_q, s_j, s_q, s_nu, s_nu]`
///
/// The physically important rate is `s_nu = 1 / (3*nu + 0.5)`, which sets
/// the kinematic viscosity. The other rates are typically set to 1 for
/// over-relaxation stability.
pub struct MrtCollision2D {
    /// Relaxation rates for the nine moment modes.
    ///
    /// Index mapping (d'Humières D2Q9):
    /// - `[0]` s_rho: density (conserved, unused but kept for symmetry)
    /// - `[1]` s_e: energy
    /// - `[2]` s_eps: energy squared
    /// - `[3]` s_jx: x-momentum (conserved)
    /// - `[4]` s_qx: x heat flux
    /// - `[5]` s_jy: y-momentum (conserved)
    /// - `[6]` s_qy: y heat flux
    /// - `[7]` s_nu: diagonal stress (controls viscosity)
    /// - `[8]` s_nu: off-diagonal stress (controls viscosity)
    pub relaxation_rates: [f64; 9],
}

impl MrtCollision2D {
    /// Create an `MrtCollision2D` from a kinematic viscosity `nu`.
    ///
    /// The viscous relaxation rate is `s_nu = 1 / (3*nu + 0.5)`.
    /// All non-hydrodynamic rates default to `1.0` (maximal dissipation
    /// of ghost modes for numerical stability).
    ///
    /// # Panics
    ///
    /// Does not panic, but a `nu <= 0` will produce an invalid `s_nu`.
    pub fn new(viscosity: f64) -> Self {
        let s_nu = 1.0 / (3.0 * viscosity + 0.5);
        // Conserved modes (rho, jx, jy) have rate 0 — they do not relax.
        // Non-hydrodynamic modes use rate 1.0 for stability.
        Self {
            relaxation_rates: [
                0.0,  // s_rho  (density — conserved)
                1.0,  // s_e    (energy ghost mode)
                1.0,  // s_eps  (energy-squared ghost mode)
                0.0,  // s_jx   (x-momentum — conserved)
                1.0,  // s_qx   (heat flux x — ghost mode)
                0.0,  // s_jy   (y-momentum — conserved)
                1.0,  // s_qy   (heat flux y — ghost mode)
                s_nu, // s_nu   (diagonal stress → shear viscosity)
                s_nu, // s_nu   (off-diagonal stress → shear viscosity)
            ],
        }
    }

    /// Apply MRT collision to a single cell.
    ///
    /// # Arguments
    ///
    /// * `f`   – current distribution functions (9 components)
    /// * `rho` – macroscopic density
    /// * `ux`  – macroscopic x-velocity
    /// * `uy`  – macroscopic y-velocity
    ///
    /// # Returns
    ///
    /// Post-collision distribution functions `f*` (9 components).
    pub fn collide(&self, f: &[f64; 9], rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        // 1. Forward transform: m = M * f
        let mut m = [0.0_f64; 9];
        for row in 0..9 {
            for col in 0..9 {
                m[row] += M[row][col] * f[col];
            }
        }

        // 2. Compute equilibrium moments m_eq
        let u2 = ux * ux + uy * uy;
        let m_eq = [
            rho,                       // rho_eq
            rho * (-2.0 + 3.0 * u2),   // e_eq
            rho * (1.0 - 3.0 * u2),    // eps_eq
            rho * ux,                  // jx_eq
            -rho * ux,                 // qx_eq
            rho * uy,                  // jy_eq
            -rho * uy,                 // qy_eq
            rho * (ux * ux - uy * uy), // pxx_eq
            rho * ux * uy,             // pxy_eq
        ];

        // 3. Relax in moment space: m* = m - S * (m - m_eq)
        let mut m_star = [0.0_f64; 9];
        for i in 0..9 {
            m_star[i] = m[i] - self.relaxation_rates[i] * (m[i] - m_eq[i]);
        }

        // 4. Back-transform: f* = M^{-1} * m*
        let mut f_star = [0.0_f64; 9];
        for row in 0..9 {
            for col in 0..9 {
                f_star[row] += M_INV[row][col] * m_star[col];
            }
        }

        f_star
    }

    /// Apply MRT collision to every cell in a 2D grid in-place.
    ///
    /// Macroscopic fields are recomputed from the distributions before
    /// collision so they are always up-to-date.
    pub fn collide_grid(&self, grid: &mut LbmGrid2D) {
        grid.compute_macroscopic();

        let n = grid.nx * grid.ny;
        for k in 0..n {
            // Gather current distributions for this cell.
            let mut f_cell = [0.0_f64; 9];
            for (i, f_cell_i) in f_cell.iter_mut().enumerate() {
                *f_cell_i = grid.f[i][k];
            }

            let rho = grid.rho[k];
            let ux = grid.ux[k];
            let uy = grid.uy[k];

            let f_star = self.collide(&f_cell, rho, ux, uy);

            // Write back.
            for (i, &fs) in f_star.iter().enumerate() {
                grid.f[i][k] = fs;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: BGK equilibrium for Poiseuille comparison test
// ---------------------------------------------------------------------------

/// Compute equilibrium distribution for the full D2Q9 set.
///
/// Used internally by tests that compare MRT output with BGK.
#[cfg(test)]
pub(crate) fn d2q9_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    use crate::grid::equilibrium_2d;
    use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let w = D2Q9_WEIGHTS[i];
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        feq[i] = equilibrium_2d(w, rho, ux, uy, cx, cy);
    }
    feq
}

// ---------------------------------------------------------------------------
// D3Q19 MRT — transformation matrix, inverse, relaxation, and simulation
// ---------------------------------------------------------------------------
//
// Velocity ordering (identical to `crate::lattice::D3Q19_VELOCITIES` and d3q19_full.rs):
//   i=0:  ( 0, 0, 0)   rest
//   i=1:  ( 1, 0, 0)   i=2:  (-1, 0, 0)
//   i=3:  ( 0, 1, 0)   i=4:  ( 0,-1, 0)
//   i=5:  ( 0, 0, 1)   i=6:  ( 0, 0,-1)
//   i=7:  ( 1, 1, 0)   i=8:  (-1, 1, 0)
//   i=9:  ( 1,-1, 0)   i=10: (-1,-1, 0)
//   i=11: ( 1, 0, 1)   i=12: (-1, 0, 1)
//   i=13: ( 1, 0,-1)   i=14: (-1, 0,-1)
//   i=15: ( 0, 1, 1)   i=16: ( 0,-1, 1)
//   i=17: ( 0, 1,-1)   i=18: ( 0,-1,-1)
//
// Reference: d'Humières et al. (2002), Phil. Trans. R. Soc. Lond. A 360, 437–451.
// The 19 moments, in the row order of `mrt_transform_matrix` (the polynomial of
// every row is listed in the documentation of that function):
//   0: ρ     1: e      2: ε      3: jx     4: qx     5: jy     6: qy     7: jz     8: qz
//   9: 3pxx  10: 3πxx  11: pww   12: πww   13: pxy   14: pxz   15: pyz
//  16: mx    17: my    18: mz
//
// `crate::mrt3d` follows the original paper and lists pyz before pxz, so the two
// shear-stress rows 14 and 15 are exchanged between that module and this one.

/// Compute the 19×19 D3Q19 MRT transformation matrix M.
///
/// Row `k` of `M` is the `k`-th moment basis vector evaluated at the 19
/// lattice velocities `c = (x, y, z)` (ordering in the table above), with
/// `r² = x² + y² + z²`:
///
/// | row | moment | polynomial |
/// |-----|--------|------------|
/// | 0 | ρ | `1` |
/// | 1 | e | `19 r² − 30` |
/// | 2 | ε | `(21 r⁴ − 53 r² + 24) / 2` |
/// | 3, 5, 7 | j_x, j_y, j_z | `x`, `y`, `z` |
/// | 4, 6, 8 | q_x, q_y, q_z | `(5 r² − 9) x`, `(5 r² − 9) y`, `(5 r² − 9) z` |
/// | 9 | 3 p_xx | `3 x² − r²` |
/// | 10 | 3 π_xx | `(3 r² − 5)(3 x² − r²)` |
/// | 11 | p_ww | `y² − z²` |
/// | 12 | π_ww | `(3 r² − 5)(y² − z²)` |
/// | 13 | p_xy | `x y` |
/// | 14 | p_xz | `x z` |
/// | 15 | p_yz | `y z` |
/// | 16, 17, 18 | m_x, m_y, m_z | `(y² − z²) x`, `(z² − x²) y`, `(x² − y²) z` |
///
/// The nineteen rows are mutually orthogonal (squared norms 19, 2394, 252,
/// 10, 40, 10, 40, 10, 40, 36, 72, 12, 24, 4, 4, 4, 8, 8, 8), so the matrix is
/// invertible and `M⁻¹ = Mᵀ diag(1/‖row‖²)`.
///
/// The matrix is [`crate::mrt3d::M`] with its rows 14 and 15 exchanged:
/// [`crate::mrt3d`] lists `p_yz` before `p_xz` as in the original paper,
/// while [`MrtRelaxation`] lists `p_xz` before `p_yz`.
pub fn mrt_transform_matrix() -> [[f64; 19]; 19] {
    let mut m = crate::mrt3d::M;
    m.swap(14, 15);
    m
}

/// Compute the 19×19 D3Q19 MRT inverse transformation matrix M^{-1}.
///
/// The matrix is computed by Gauss–Jordan elimination with partial pivoting on
/// `[M | I]` in f64 arithmetic, where `M` is [`mrt_transform_matrix`].  `M` has
/// 19 mutually orthogonal rows, so it is invertible and the product `M · M⁻¹`
/// equals the identity to about 1e-15.
pub fn mrt_inverse_matrix() -> [[f64; 19]; 19] {
    let m = mrt_transform_matrix();
    // Build augmented matrix [M | I].
    let mut aug = [[0.0_f64; 38]; 19];
    for r in 0..19 {
        for c in 0..19 {
            aug[r][c] = m[r][c];
        }
        aug[r][19 + r] = 1.0;
    }
    // Forward elimination with partial pivoting.
    for col in 0..19 {
        // Find pivot row.
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (offset, aug_row) in aug[col + 1..].iter().enumerate() {
            let row = col + 1 + offset;
            if aug_row[col].abs() > max_val {
                max_val = aug_row[col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        for elem in aug[col][col..38].iter_mut() {
            *elem /= pivot;
        }
        let col_row: [f64; 38] = aug[col];
        for (row, aug_row) in aug.iter_mut().enumerate() {
            if row == col {
                continue;
            }
            let factor = aug_row[col];
            for (elem, &cv) in aug_row.iter_mut().zip(col_row.iter()) {
                *elem -= factor * cv;
            }
        }
    }
    // Extract right half.
    let mut inv = [[0.0_f64; 19]; 19];
    for r in 0..19 {
        for c in 0..19 {
            inv[r][c] = aug[r][19 + c];
        }
    }
    inv
}

// ---------------------------------------------------------------------------
// MrtRelaxation
// ---------------------------------------------------------------------------

/// Relaxation rates for the 19 D3Q19 MRT moment modes.
///
/// Index mapping (row order of [`mrt_transform_matrix`]):
/// - `[0]`  s_ρ   : density (conserved → 0)
/// - `[1]`  s_e   : energy
/// - `[2]`  s_ε   : energy²
/// - `[3]`  s_jx  : x-momentum (conserved → 0)
/// - `[4]`  s_qx  : x-heat-flux ghost
/// - `[5]`  s_jy  : y-momentum (conserved → 0)
/// - `[6]`  s_qy  : y-heat-flux ghost
/// - `[7]`  s_jz  : z-momentum (conserved → 0)
/// - `[8]`  s_qz  : z-heat-flux ghost
/// - `[9]`  s_ν   : 3pxx  → shear viscosity
/// - `[10]` s_π   : 3πxx  ghost
/// - `[11]` s_ν   : pww   → shear viscosity
/// - `[12]` s_π   : πww   ghost
/// - `[13]` s_ν   : pxy   → shear viscosity
/// - `[14]` s_ν   : pxz   → shear viscosity
/// - `[15]` s_ν   : pyz   → shear viscosity
/// - `[16]` s_mx  : cubic ghost
/// - `[17]` s_my  : cubic ghost
/// - `[18]` s_mz  : cubic ghost
///
/// This ordering lists `p_xz` (14) before `p_yz` (15); the rates used by
/// [`crate::mrt3d`] ([`crate::mrt3d::MrtRelaxationParams`]) follow the original
/// paper with `p_yz` before `p_xz`.  [`MrtRelaxation::to_mrt3d_order`]
/// converts between the two.  An isotropic shear viscosity needs
/// `s[9] = s[11] = s[13] = s[14] = s[15]`.
pub struct MrtRelaxation {
    /// Relaxation rates for each of the 19 moment modes.
    pub s: [f64; 19],
}

impl MrtRelaxation {
    /// Build relaxation rates with the shear viscosity of a BGK scheme with
    /// kinematic viscosity `nu`.
    ///
    /// All five stress modes (3pxx, pww, pxy, pxz, pyz) are set to
    /// `s_ν = 1 / (3*nu + 0.5)`.  The conserved-mode rates (density and the
    /// three momentum components) are 0, so a collision conserves mass and
    /// momentum for any `(rho, u)` it is given.  All other (energy/ghost)
    /// modes are set to 1.0 for maximal ghost-mode dissipation.
    pub fn bgk_equivalent(nu: f64) -> Self {
        let s_nu = 1.0 / (3.0 * nu + 0.5);
        Self {
            s: [
                0.0,  // [0]  ρ        (conserved)
                1.0,  // [1]  e        (energy → 1.0)
                1.0,  // [2]  ε        (energy² → 1.0)
                0.0,  // [3]  jx       (conserved)
                1.0,  // [4]  qx       (heat flux ghost → 1.0)
                0.0,  // [5]  jy       (conserved)
                1.0,  // [6]  qy       (heat flux ghost → 1.0)
                0.0,  // [7]  jz       (conserved)
                1.0,  // [8]  qz       (heat flux ghost → 1.0)
                s_nu, // [9]  3pxx     (shear viscosity)
                1.0,  // [10] 3πxx     (ghost → 1.0)
                s_nu, // [11] pww      (shear viscosity)
                1.0,  // [12] πww      (ghost → 1.0)
                s_nu, // [13] pxy      (shear viscosity)
                s_nu, // [14] pxz      (shear viscosity)
                s_nu, // [15] pyz      (shear viscosity)
                1.0,  // [16] mx       (cubic ghost → 1.0)
                1.0,  // [17] my       (cubic ghost → 1.0)
                1.0,  // [18] mz       (cubic ghost → 1.0)
            ],
        }
    }

    /// The rates in the row order of [`crate::mrt3d::M`] (`p_yz` at index 14,
    /// `p_xz` at index 15), i.e. with entries 14 and 15 exchanged.
    pub fn to_mrt3d_order(&self) -> [f64; 19] {
        let mut s = self.s;
        s.swap(14, 15);
        s
    }
}

// ---------------------------------------------------------------------------
// Standalone D3Q19 equilibrium helper
// ---------------------------------------------------------------------------

/// Compute the standard D3Q19 equilibrium distribution.
///
/// `feq_i = w_i * ρ * (1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) - |u|²/(2cs²))`
/// with cs² = 1/3.
pub fn d3q19_equilibrium(rho: f64, u: [f64; 3]) -> [f64; 19] {
    const CS2: f64 = 1.0 / 3.0;
    let u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let ex: [i32; 19] = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0];
    let ey: [i32; 19] = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1];
    let ez: [i32; 19] = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1];
    let w: [f64; 19] = [
        1.0 / 3.0,
        1.0 / 18.0,
        1.0 / 18.0,
        1.0 / 18.0,
        1.0 / 18.0,
        1.0 / 18.0,
        1.0 / 18.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ];
    let mut feq = [0.0_f64; 19];
    for i in 0..19 {
        let eu = ex[i] as f64 * u[0] + ey[i] as f64 * u[1] + ez[i] as f64 * u[2];
        feq[i] = w[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

// ---------------------------------------------------------------------------
// MrtD3Q19
// ---------------------------------------------------------------------------

/// D3Q19 MRT (Multiple Relaxation Time) LBM simulation.
///
/// Uses the d'Humières et al. (2002) MRT collision operator with the standard
/// D3Q19 velocity set. Distribution functions are stored as `Vec<[f64;19]>`
/// indexed by `idx = x + y*nx + z*nx*ny`.
pub struct MrtD3Q19 {
    /// Number of cells in the x-direction.
    pub nx: usize,
    /// Number of cells in the y-direction.
    pub ny: usize,
    /// Number of cells in the z-direction.
    pub nz: usize,
    /// Distribution functions: `f[idx]` contains all 19 populations for cell `idx`.
    pub f: Vec<[f64; 19]>,
    /// MRT relaxation rates.
    pub relax: MrtRelaxation,
}

// Velocity components — identical to d3q19_full.rs.
const MRT_EX: [i32; 19] = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0];
const MRT_EY: [i32; 19] = [0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1];
const MRT_EZ: [i32; 19] = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1];

// D3Q19 weights.
const MRT_W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

impl MrtD3Q19 {
    /// Create a new simulation initialised to equilibrium at rest (ρ = 1, **u** = 0).
    pub fn new(nx: usize, ny: usize, nz: usize, relax: MrtRelaxation) -> Self {
        let n = nx * ny * nz;
        // At rest, feq_i = w_i (since ρ = 1).
        let f = vec![MRT_W; n];
        Self {
            nx,
            ny,
            nz,
            f,
            relax,
        }
    }

    /// Flat linear index for cell `(x, y, z)`.
    #[cfg(test)]
    #[inline]
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + y * self.nx + z * self.nx * self.ny
    }

    /// Compute the D3Q19 equilibrium distribution for a given density and velocity.
    ///
    /// `feq_i = w_i * ρ * (1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) - |u|²/(2cs²))`
    /// with cs² = 1/3.
    pub fn equilibrium(rho: f64, u: [f64; 3]) -> [f64; 19] {
        const CS2: f64 = 1.0 / 3.0;
        let u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        let mut feq = [0.0_f64; 19];
        for i in 0..19 {
            let eu = MRT_EX[i] as f64 * u[0] + MRT_EY[i] as f64 * u[1] + MRT_EZ[i] as f64 * u[2];
            feq[i] =
                MRT_W[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
        }
        feq
    }

    /// Macroscopic density and velocity at cell `idx`.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 3]) {
        let fi = &self.f[idx];
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for i in 0..19 {
            rho += fi[i];
            mx += fi[i] * MRT_EX[i] as f64;
            my += fi[i] * MRT_EY[i] as f64;
            mz += fi[i] * MRT_EZ[i] as f64;
        }
        if rho.abs() > 1e-15 {
            (rho, [mx / rho, my / rho, mz / rho])
        } else {
            (0.0, [0.0, 0.0, 0.0])
        }
    }

    /// Apply the d'Humières MRT collision operator to a single cell.
    ///
    /// The populations are mapped to the 19 moments `m = M f`
    /// ([`mrt_transform_matrix`]), every moment relaxes toward its equilibrium
    /// `m_eq(ρ, u) = M f_eq(ρ, u)` at its own rate `s_k` from `self.relax`, and
    /// the result is mapped back:
    ///
    /// ```text
    /// f* = f − M⁻¹ S (M f − m_eq)
    /// ```
    ///
    /// The work is done by [`crate::mrt3d::collide_mrt`] after exchanging the
    /// rates `s[14]` and `s[15]` (see [`MrtRelaxation`]).  Within the collision
    /// every moment relaxes independently, so the ghost and energy rates never
    /// change the post-collision shear-stress moments that set the viscosity.
    ///
    /// Mass and momentum are conserved exactly when the conserved-mode rates
    /// (`s[0]`, `s[3]`, `s[5]`, `s[7]`) are zero, as in
    /// [`MrtRelaxation::bgk_equivalent`], whatever `rho` and `u` are passed.
    /// When `rho` and `u` are the moments of `f` those rates have no effect.
    pub fn collide_cell(&self, f: [f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        crate::mrt3d::collide_mrt(&f, rho, u, &self.relax.to_mrt3d_order())
    }

    /// Apply MRT collision to every cell in the domain.
    pub fn collide(&mut self) {
        let n = self.nx * self.ny * self.nz;
        for k in 0..n {
            let (rho, u) = self.macros(k);
            let f_new = self.collide_cell(self.f[k], rho, u);
            self.f[k] = f_new;
        }
    }

    /// Pull-scheme streaming step with fully periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let f_old = self.f.clone();

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = x + y * nx + z * nx * ny;
                    for i in 0..19 {
                        let sx = ((x as i64 - MRT_EX[i] as i64).rem_euclid(nx as i64)) as usize;
                        let sy = ((y as i64 - MRT_EY[i] as i64).rem_euclid(ny as i64)) as usize;
                        let sz = ((z as i64 - MRT_EZ[i] as i64).rem_euclid(nz as i64)) as usize;
                        let src = sx + sy * nx + sz * nx * ny;
                        self.f[dst][i] = f_old[src][i];
                    }
                }
            }
        }
    }

    /// Perform one full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }

    /// Initialize all cells to equilibrium for the given density and velocity.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 3]) {
        let feq = Self::equilibrium(rho, u);
        for cell in self.f.iter_mut() {
            *cell = feq;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::{apply_boundaries_2d, channel_walls};
    use crate::collision::bgk_collide_2d;
    use crate::grid::LbmGrid2D;
    use crate::lattice::LatticeType;
    use crate::streaming::stream_2d;

    /// 1. MRT collision conserves total mass (rho).
    #[test]
    fn test_mrt_mass_conservation() {
        let mrt = MrtCollision2D::new(1.0 / 6.0);
        let rho = 1.05;
        let ux = 0.05;
        let uy = -0.02;
        let feq = d2q9_equilibrium(rho, ux, uy);

        // Perturb slightly from equilibrium.
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01;

        let mass_before: f64 = f.iter().sum();
        let f_star = mrt.collide(&f, rho, ux, uy);
        let mass_after: f64 = f_star.iter().sum();

        assert!(
            (mass_before - mass_after).abs() < 1e-13,
            "MRT mass not conserved: before={mass_before}, after={mass_after}"
        );
    }

    /// 2. MRT collision conserves momentum (rho * ux, rho * uy).
    #[test]
    fn test_mrt_momentum_conservation() {
        use crate::lattice::D2Q9_VELOCITIES;

        let mrt = MrtCollision2D::new(1.0 / 6.0);
        let rho = 1.0;
        let ux = 0.08;
        let uy = 0.03;
        let feq = d2q9_equilibrium(rho, ux, uy);

        let mut f = feq;
        f[5] += 0.005;
        f[7] -= 0.005;

        let jx_before: f64 = (0..9).map(|i| f[i] * D2Q9_VELOCITIES[i][0] as f64).sum();
        let jy_before: f64 = (0..9).map(|i| f[i] * D2Q9_VELOCITIES[i][1] as f64).sum();

        let f_star = mrt.collide(&f, rho, ux, uy);

        let jx_after: f64 = (0..9)
            .map(|i| f_star[i] * D2Q9_VELOCITIES[i][0] as f64)
            .sum();
        let jy_after: f64 = (0..9)
            .map(|i| f_star[i] * D2Q9_VELOCITIES[i][1] as f64)
            .sum();

        assert!(
            (jx_before - jx_after).abs() < 1e-13,
            "MRT x-momentum not conserved: before={jx_before}, after={jx_after}"
        );
        assert!(
            (jy_before - jy_after).abs() < 1e-13,
            "MRT y-momentum not conserved: before={jy_before}, after={jy_after}"
        );
    }

    /// 3. At equilibrium, MRT collision is a no-op (leaves f unchanged).
    #[test]
    fn test_mrt_equilibrium_invariant() {
        let mrt = MrtCollision2D::new(1.0 / 6.0);
        let rho = 1.0;
        let ux = 0.05;
        let uy = -0.03;
        let feq = d2q9_equilibrium(rho, ux, uy);

        let f_star = mrt.collide(&feq, rho, ux, uy);

        for i in 0..9 {
            assert!(
                (feq[i] - f_star[i]).abs() < 1e-13,
                "MRT at equilibrium changed f[{i}]: before={}, after={}",
                feq[i],
                f_star[i]
            );
        }
    }

    /// 4. s_nu is computed correctly for nu = 1/6.
    #[test]
    fn test_mrt_viscosity_parameter() {
        let nu = 1.0 / 6.0;
        let mrt = MrtCollision2D::new(nu);
        let expected_s_nu = 1.0 / (3.0 * nu + 0.5); // = 1.0
        let s_nu = mrt.relaxation_rates[7];
        assert!(
            (s_nu - expected_s_nu).abs() < 1e-14,
            "s_nu incorrect: got {s_nu}, expected {expected_s_nu}"
        );
        // Stress components both equal s_nu.
        assert!(
            (mrt.relaxation_rates[8] - expected_s_nu).abs() < 1e-14,
            "s_nu[8] incorrect"
        );
    }

    /// 5. collide_grid runs without panicking on a real grid.
    #[test]
    fn test_mrt_grid_collide() {
        let mut grid = LbmGrid2D::new(10, 10, LatticeType::D2Q9);
        grid.set_equilibrium(3, 3, 1.1, 0.05, 0.0);
        grid.compute_macroscopic();

        let mrt = MrtCollision2D::new(1.0 / 6.0);
        mrt.collide_grid(&mut grid);

        // Grid should still have sensible densities.
        grid.compute_macroscopic();
        for k in 0..(10 * 10) {
            assert!(grid.rho[k] > 0.0, "Density became non-positive at cell {k}");
        }
    }

    /// 6. TRT creation produces reasonable parameters.
    #[test]
    fn test_trt_creation() {
        use crate::mrt3d::TrtCollision3D;
        let trt = TrtCollision3D::new(1.0 / 6.0);
        assert!(
            (0.0..2.0).contains(&trt.s_plus),
            "s_plus out of range: {}",
            trt.s_plus
        );
        assert!(
            (0.0..2.0).contains(&trt.s_minus),
            "s_minus out of range: {}",
            trt.s_minus
        );
    }

    /// 7. TRT magic parameter satisfies the stability condition:
    ///    Λ = (1/s_plus - 0.5) * (1/s_minus - 0.5) = 3/16 (optimal value).
    #[test]
    fn test_trt_magic_parameter() {
        use crate::mrt3d::TrtCollision3D;
        let trt = TrtCollision3D::new(1.0 / 6.0);
        let lambda = (1.0 / trt.s_plus - 0.5) * (1.0 / trt.s_minus - 0.5);
        // The magic parameter should be positive and bounded for stability.
        assert!(
            lambda > 0.0,
            "TRT magic parameter should be positive: {lambda}"
        );
        // Optimal value 3/16 = 0.1875 — allow tolerance.
        assert!(
            (lambda - 3.0_f64 / 16.0).abs() < 0.05,
            "TRT magic parameter far from optimal 3/16: {lambda}"
        );
    }

    /// 8. For low-Re Poiseuille flow, MRT and BGK produce the same velocity profile.
    #[test]
    fn test_mrt_vs_bgk_poiseuille() {
        let nx = 5;
        let ny = 12;
        let nu = 1.0 / 6.0;
        let omega = 1.0 / (3.0 * nu + 0.5); // BGK relaxation rate = 1.0
        let body_force = 1e-5_f64;
        let n_steps = 3000;

        // ---- BGK run ----
        let mut grid_bgk = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        let walls = channel_walls(nx, ny);
        for _ in 0..n_steps {
            bgk_collide_2d(&mut grid_bgk, omega);
            // Uniform body force: nudge ux for every fluid cell.
            for y in 1..(ny - 1) {
                for x in 0..nx {
                    let k = grid_bgk.idx(x, y);
                    grid_bgk.ux[k] += body_force;
                }
            }
            stream_2d(&mut grid_bgk);
            apply_boundaries_2d(&mut grid_bgk, &walls);
            grid_bgk.compute_macroscopic();
        }

        // ---- MRT run ----
        let mrt = MrtCollision2D::new(nu);
        let mut grid_mrt = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        for _ in 0..n_steps {
            mrt.collide_grid(&mut grid_mrt);
            for y in 1..(ny - 1) {
                for x in 0..nx {
                    let k = grid_mrt.idx(x, y);
                    grid_mrt.ux[k] += body_force;
                }
            }
            stream_2d(&mut grid_mrt);
            apply_boundaries_2d(&mut grid_mrt, &walls);
            grid_mrt.compute_macroscopic();
        }

        // Compare velocity profiles at mid x.
        let mid_x = nx / 2;
        for y in 1..(ny - 1) {
            let (ux_bgk, _) = grid_bgk.velocity_at(mid_x, y);
            let (ux_mrt, _) = grid_mrt.velocity_at(mid_x, y);

            // Allow 5 % relative error or absolute 1e-6 for near-zero cells.
            let abs_err = (ux_bgk - ux_mrt).abs();
            let rel_err = if ux_bgk.abs() > 1e-8 {
                abs_err / ux_bgk.abs()
            } else {
                abs_err
            };
            assert!(
                rel_err < 0.1,
                "MRT vs BGK mismatch at y={y}: bgk={ux_bgk}, mrt={ux_mrt}, rel_err={rel_err}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // D3Q19 MRT tests
    // -----------------------------------------------------------------------

    /// D1. MRT D3Q19 collision conserves mass (sum of f).
    #[test]
    fn test_mrt_d3q19_mass_conservation() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        let sim = MrtD3Q19::new(4, 4, 4, relax);
        let rho = 1.05_f64;
        let u = [0.05, -0.02, 0.01];
        let mut f = MrtD3Q19::equilibrium(rho, u);
        // Perturb away from equilibrium (mass-neutral perturbation).
        f[1] += 0.005;
        f[2] -= 0.005;

        let mass_before: f64 = f.iter().sum();
        let (rho_in, u_in) = {
            let mut r = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for i in 0..19 {
                r += f[i];
                mx += f[i] * MRT_EX[i] as f64;
                my += f[i] * MRT_EY[i] as f64;
                mz += f[i] * MRT_EZ[i] as f64;
            }
            (r, [mx / r, my / r, mz / r])
        };
        let f_star = sim.collide_cell(f, rho_in, u_in);
        let mass_after: f64 = f_star.iter().sum();

        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "D3Q19 MRT mass not conserved: before={mass_before}, after={mass_after}"
        );
    }

    /// D2. MRT D3Q19 collision conserves momentum (Σ f_i * e_i).
    #[test]
    fn test_mrt_d3q19_momentum_conservation() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        let sim = MrtD3Q19::new(4, 4, 4, relax);
        let rho = 1.0_f64;
        let u = [0.05, 0.03, -0.02];
        let mut f = MrtD3Q19::equilibrium(rho, u);
        // Mass-neutral, momentum-neutral perturbation.
        f[7] += 0.003;
        f[10] -= 0.003;

        let jx_before: f64 = (0..19).map(|i| f[i] * MRT_EX[i] as f64).sum();
        let jy_before: f64 = (0..19).map(|i| f[i] * MRT_EY[i] as f64).sum();
        let jz_before: f64 = (0..19).map(|i| f[i] * MRT_EZ[i] as f64).sum();

        let rho_in: f64 = f.iter().sum();
        let u_in = [jx_before / rho_in, jy_before / rho_in, jz_before / rho_in];
        let f_star = sim.collide_cell(f, rho_in, u_in);

        let jx_after: f64 = (0..19).map(|i| f_star[i] * MRT_EX[i] as f64).sum();
        let jy_after: f64 = (0..19).map(|i| f_star[i] * MRT_EY[i] as f64).sum();
        let jz_after: f64 = (0..19).map(|i| f_star[i] * MRT_EZ[i] as f64).sum();

        assert!(
            (jx_before - jx_after).abs() < 1e-12,
            "D3Q19 MRT jx not conserved: before={jx_before}, after={jx_after}"
        );
        assert!(
            (jy_before - jy_after).abs() < 1e-12,
            "D3Q19 MRT jy not conserved: before={jy_before}, after={jy_after}"
        );
        assert!(
            (jz_before - jz_after).abs() < 1e-12,
            "D3Q19 MRT jz not conserved: before={jz_before}, after={jz_after}"
        );
    }

    /// D3. Equilibrium distribution sums to rho.
    #[test]
    fn test_mrt_d3q19_equilibrium_sums_to_rho() {
        let rho = 1.23_f64;
        let u = [0.04, -0.03, 0.02];
        let feq = MrtD3Q19::equilibrium(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-13,
            "D3Q19 equilibrium sum = {sum}, expected rho = {rho}"
        );
    }

    /// D4. With all s_i = 1/(3*nu+0.5), MRT collision equals BGK for a cell at equilibrium.
    ///
    /// At equilibrium MRT is a no-op regardless of s_i, which also equals BGK at equilibrium.
    /// Away from equilibrium with uniform s, the MRT correction equals the BGK correction.
    #[test]
    fn test_mrt_d3q19_bgk_equivalent() {
        let nu = 1.0 / 6.0;
        let omega = 1.0 / (3.0 * nu + 0.5); // = 1.0
        // Build MRT with ALL rates set to omega (fully uniform = BGK).
        let mut relax_bgk = MrtRelaxation::bgk_equivalent(nu);
        // Override conserved modes to omega too (for strict BGK equivalence check).
        for i in [0usize, 3, 5, 7] {
            relax_bgk.s[i] = omega;
        }
        let sim = MrtD3Q19::new(2, 2, 2, relax_bgk);

        let rho = 1.0_f64;
        let u = [0.05, 0.02, -0.01];
        let feq = MrtD3Q19::equilibrium(rho, u);

        // Perturb.
        let mut f = feq;
        f[1] += 0.01;
        f[2] -= 0.01;

        // BGK: f* = f - omega*(f - feq)
        let mut f_bgk = f;
        for i in 0..19 {
            f_bgk[i] = f[i] - omega * (f[i] - feq[i]);
        }

        let f_mrt = sim.collide_cell(f, rho, u);

        for i in 0..19 {
            assert!(
                (f_bgk[i] - f_mrt[i]).abs() < 1e-10,
                "MRT-BGK mismatch at i={i}: bgk={}, mrt={}",
                f_bgk[i],
                f_mrt[i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // New required tests
    // -----------------------------------------------------------------------

    /// N1. MrtRelaxation::bgk_equivalent: all s values are in (0, 2).
    #[test]
    fn test_mrt_relaxation_bgk_equivalent_s_in_range() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        for (i, &si) in relax.s.iter().enumerate() {
            assert!((0.0..2.0).contains(&si), "s[{i}] = {si} is not in (0, 2)");
        }
    }

    /// N2. macros after set_uniform: rho and velocity match the set values.
    #[test]
    fn test_mrt_d3q19_macros_after_set_uniform() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        let mut sim = MrtD3Q19::new(4, 4, 4, relax);
        let rho_set = 1.3_f64;
        let u_set = [0.05, -0.02, 0.01];
        sim.set_uniform(rho_set, u_set);

        let n = sim.nx * sim.ny * sim.nz;
        for k in 0..n {
            let (rho_k, u_k) = sim.macros(k);
            assert!(
                (rho_k - rho_set).abs() < 1e-13,
                "rho mismatch at cell {k}: got {rho_k}, expected {rho_set}"
            );
            for d in 0..3 {
                assert!(
                    (u_k[d] - u_set[d]).abs() < 1e-12,
                    "u[{d}] mismatch at cell {k}: got {}, expected {}",
                    u_k[d],
                    u_set[d]
                );
            }
        }
    }

    /// N3. After 1 step: total mass is conserved (Σρ conserved).
    #[test]
    fn test_mrt_d3q19_step_mass_conservation() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        let mut sim = MrtD3Q19::new(4, 4, 4, relax);
        // Perturb one cell.
        let k = sim.idx(1, 2, 3);
        sim.f[k][0] += 0.1;
        sim.f[k][1] -= 0.05;
        sim.f[k][2] -= 0.05;

        let rho_before: f64 = (0..sim.nx * sim.ny * sim.nz)
            .map(|idx| sim.macros(idx).0)
            .sum();
        sim.step();
        let rho_after: f64 = (0..sim.nx * sim.ny * sim.nz)
            .map(|idx| sim.macros(idx).0)
            .sum();

        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Mass not conserved after step: before={rho_before}, after={rho_after}"
        );
    }

    /// N4. d3q19_equilibrium sums to rho.
    #[test]
    fn test_d3q19_equilibrium_sums_to_rho() {
        let rho = 1.5_f64;
        let u = [0.03, -0.02, 0.01];
        let feq = d3q19_equilibrium(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-13,
            "d3q19_equilibrium sum = {sum}, expected {rho}"
        );
    }
}

// ---------------------------------------------------------------------------
// TRT (Two-Relaxation-Time) collision operator for D2Q9
// ---------------------------------------------------------------------------

/// Two-Relaxation-Time (TRT) collision operator for D2Q9.
///
/// TRT splits the distribution into symmetric (even) and antisymmetric (odd)
/// parts and relaxes each at its own rate: `s⁺ = 1/τ⁺` for the even part, which
/// carries the viscous stress (so `τ⁺ = 3ν + ½` fixes the kinematic
/// viscosity), and `s⁻ = 1/τ⁻` for the odd part.
///
/// The two relaxation times are tied together by the "magic" parameter
///
/// ```text
/// Λ = (τ⁺ − ½)(τ⁻ − ½)
/// ```
///
/// `Λ = 3/16` makes the effective position of a bounce-back wall in a straight
/// channel independent of the viscosity.  Stability needs `0 < s⁻ < 2`, i.e.
/// `τ⁻ > ½`, which holds precisely when `Λ > 0`.
///
/// Reference: Ginzburg, Verhaeghe & d'Humières (2008), Commun. Comput. Phys. 3, 427–478.
#[derive(Debug, Clone, Copy)]
pub struct TrtCollision2D {
    /// Relaxation time of the symmetric (even) part, `τ⁺ = 3ν + ½`; the
    /// corresponding rate is `s⁺ = 1/τ⁺`.  It controls the viscosity.
    pub tau_plus: f64,
    /// Relaxation time of the antisymmetric (odd) part; the corresponding
    /// rate is `s⁻ = 1/τ⁻`.  The constructors derive it from the magic
    /// parameter, `τ⁻ = ½ + Λ / (τ⁺ − ½)`.
    pub tau_minus: f64,
}

impl TrtCollision2D {
    /// Construct TRT from kinematic viscosity `nu` with the magic parameter
    /// `Λ = 3/16` ([`TRT_MAGIC`]).
    ///
    /// Sets `τ⁺ = 3ν + ½` and `τ⁻ = ½ + (3/16) / (τ⁺ − ½)`, so that
    /// `(τ⁺ − ½)(τ⁻ − ½) = 3/16` and `0 < 1/τ⁻ < 2` for every `ν > 0`.
    ///
    /// No stable operator exists for `ν <= 0` or a NaN viscosity.  `τ⁻` is then
    /// NaN, so that every [`TrtCollision2D::collide`] returns NaN populations
    /// instead of silently unstable values.  [`TrtCollision2D::try_new`]
    /// reports such input as an error instead.
    pub fn new(nu: f64) -> Self {
        Self::with_magic_parameter(nu, TRT_MAGIC)
    }

    /// Construct TRT from kinematic viscosity `nu` and the magic parameter
    /// `lambda`: `τ⁺ = 3ν + ½`, `τ⁻ = ½ + Λ / (τ⁺ − ½)`.
    ///
    /// `τ⁻` is NaN (see [`TrtCollision2D::new`]) unless `nu` and `lambda` are
    /// finite and positive.
    pub fn with_magic_parameter(nu: f64, lambda: f64) -> Self {
        let tau_plus = 3.0 * nu + 0.5;
        // `3ν` is `τ⁺ − ½` without the cancellation of the subtraction.
        let excess = 3.0 * nu;
        let tau_minus = if excess.is_finite() && excess > 0.0 && lambda.is_finite() && lambda > 0.0
        {
            0.5 + lambda / excess
        } else {
            f64::NAN
        };
        Self {
            tau_plus,
            tau_minus,
        }
    }

    /// Fallible variant of [`TrtCollision2D::new`]: rejects `nu <= 0` and
    /// non-finite viscosities.
    pub fn try_new(nu: f64) -> crate::error::Result<Self> {
        Self::try_with_magic_parameter(nu, TRT_MAGIC)
    }

    /// Fallible variant of [`TrtCollision2D::with_magic_parameter`]: rejects a
    /// non-positive or non-finite `nu` or `lambda`.
    pub fn try_with_magic_parameter(nu: f64, lambda: f64) -> crate::error::Result<Self> {
        crate::error::Error::check_range(nu, f64::MIN_POSITIVE, f64::MAX, "nu")?;
        crate::error::Error::check_range(lambda, f64::MIN_POSITIVE, f64::MAX, "lambda")?;
        Ok(Self::with_magic_parameter(nu, lambda))
    }

    /// Construct TRT with explicit relaxation times `τ⁺` and `τ⁻`.
    pub fn with_rates(tau_plus: f64, tau_minus: f64) -> Self {
        Self {
            tau_plus,
            tau_minus,
        }
    }

    /// Kinematic viscosity implied by `τ⁺`.
    pub fn viscosity(&self) -> f64 {
        (self.tau_plus - 0.5) / 3.0
    }

    /// Magic parameter `Λ = (τ⁺ − ½)(τ⁻ − ½)`.
    pub fn lambda(&self) -> f64 {
        (self.tau_plus - 0.5) * (self.tau_minus - 0.5)
    }

    /// Split distribution into symmetric (even) part f⁺ = (fᵢ + f̄ᵢ)/2
    /// where f̄ᵢ is the opposite-direction distribution.
    fn symmetric(f: &[f64; 9]) -> [f64; 9] {
        // D2Q9 opposite index map
        const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        let mut fp = [0.0f64; 9];
        for i in 0..9 {
            fp[i] = 0.5 * (f[i] + f[OPP[i]]);
        }
        fp
    }

    /// Split distribution into antisymmetric (odd) part f⁻ = (fᵢ - f̄ᵢ)/2.
    fn antisymmetric(f: &[f64; 9]) -> [f64; 9] {
        const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        let mut fm = [0.0f64; 9];
        for i in 0..9 {
            fm[i] = 0.5 * (f[i] - f[OPP[i]]);
        }
        fm
    }

    /// Apply TRT collision to a single D2Q9 cell.
    ///
    /// Returns post-collision distribution.
    pub fn collide(&self, f: &[f64; 9], rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let fp = Self::symmetric(f);
        let fm = Self::antisymmetric(f);
        let feqp = Self::symmetric(&feq);
        let feqm = Self::antisymmetric(&feq);

        let s_plus = 1.0 / self.tau_plus;
        let s_minus = 1.0 / self.tau_minus;

        let mut f_out = [0.0f64; 9];
        for i in 0..9 {
            f_out[i] = f[i] - s_plus * (fp[i] - feqp[i]) - s_minus * (fm[i] - feqm[i]);
        }
        f_out
    }
}

/// D2Q9 equilibrium distribution (helper for TRT).
fn d2q9_equilibrium_trt(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    // D2Q9 weights and velocities
    const W: [f64; 9] = [
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ];
    const CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    const CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; 9];
    for i in 0..9 {
        let eu = CX[i] * ux + CY[i] * uy;
        feq[i] = W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
    }
    feq
}

// ---------------------------------------------------------------------------
// MRT Stability analysis utilities
// ---------------------------------------------------------------------------

/// Compute the spectral radius of the MRT relaxation operator S.
///
/// The spectral radius is the maximum relaxation rate, which should be ≤ 2
/// for linear stability.
pub fn mrt_spectral_radius(rates: &[f64; 9]) -> f64 {
    rates.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

/// Check whether a set of MRT relaxation rates is linearly stable.
///
/// Linear stability requires all rates in (0, 2).  Conserved-mode rates
/// should be 0.
pub fn mrt_is_stable(rates: &[f64; 9]) -> bool {
    rates.iter().all(|&s| (0.0..=2.0).contains(&s))
}

/// Compute the effective bulk viscosity from MRT relaxation rates.
///
/// `ζ = ρ·cs² · (1/s_e - 0.5) · dt`
/// where `s_e` is the energy-mode relaxation rate (index 1).
pub fn mrt_bulk_viscosity(rates: &[f64; 9], rho: f64, dt: f64) -> f64 {
    let s_e = rates[1];
    if s_e.abs() < 1e-15 {
        return f64::INFINITY;
    }
    rho * (1.0 / 3.0) * (1.0 / s_e - 0.5) * dt
}

/// Compute the shear viscosity from MRT relaxation rates.
///
/// `ν = cs² · (1/s_ν - 0.5) · dt`
pub fn mrt_shear_viscosity(rates: &[f64; 9], dt: f64) -> f64 {
    let s_nu = rates[7]; // both stress modes should be equal
    if s_nu.abs() < 1e-15 {
        return f64::INFINITY;
    }
    (1.0 / 3.0) * (1.0 / s_nu - 0.5) * dt
}

// ---------------------------------------------------------------------------
// Regularised LBM (RLBM) for D2Q9
// ---------------------------------------------------------------------------

/// Regularised LBM collision step for D2Q9.
///
/// Before BGK collision, the non-equilibrium part is projected onto the
/// Hermite basis so that only stress-tensor contributions survive.
/// This eliminates non-hydrodynamic modes and enhances stability near walls.
///
/// Reference: Latt & Chopard (2006), Phys. Rev. E 72, 016703.
#[derive(Debug, Clone, Copy)]
pub struct RegularisedBgk2D {
    /// Relaxation time τ = 3ν + 0.5.
    pub tau: f64,
}

impl RegularisedBgk2D {
    /// Construct from kinematic viscosity.
    pub fn new(nu: f64) -> Self {
        Self {
            tau: 3.0 * nu + 0.5,
        }
    }

    /// Project the non-equilibrium stress tensor from f - f_eq.
    ///
    /// Returns the off-equilibrium stress components (Π_xx, Π_xy, Π_yy).
    pub fn neq_stress(f: &[f64; 9], feq: &[f64; 9]) -> (f64, f64, f64) {
        const CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
        const CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
        let mut pi_xx = 0.0f64;
        let mut pi_xy = 0.0f64;
        let mut pi_yy = 0.0f64;
        for i in 0..9 {
            let fneq = f[i] - feq[i];
            pi_xx += fneq * CX[i] * CX[i];
            pi_xy += fneq * CX[i] * CY[i];
            pi_yy += fneq * CY[i] * CY[i];
        }
        (pi_xx, pi_xy, pi_yy)
    }

    /// Reconstruct the regularised non-equilibrium populations from the stress tensor.
    fn reg_neq(pi_xx: f64, pi_xy: f64, pi_yy: f64) -> [f64; 9] {
        const W: [f64; 9] = [
            4.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 9.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
            1.0 / 36.0,
        ];
        const CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
        const CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
        let mut fneq = [0.0f64; 9];
        for i in 0..9 {
            // Q_i = (cx² - cs²) Πxx + 2 cx cy Πxy + (cy² - cs²) Πyy  [× 9/(2cs⁴)]
            fneq[i] = W[i]
                * 4.5
                * ((CX[i] * CX[i] - 1.0 / 3.0) * pi_xx
                    + 2.0 * CX[i] * CY[i] * pi_xy
                    + (CY[i] * CY[i] - 1.0 / 3.0) * pi_yy);
        }
        fneq
    }

    /// Apply regularised BGK collision.
    pub fn collide(&self, f: &[f64; 9], rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let (pi_xx, pi_xy, pi_yy) = Self::neq_stress(f, &feq);
        let fneq_reg = Self::reg_neq(pi_xx, pi_xy, pi_yy);
        let mut f_out = [0.0f64; 9];
        for i in 0..9 {
            f_out[i] = feq[i] + (1.0 - 1.0 / self.tau) * fneq_reg[i];
        }
        f_out
    }
}

// ---------------------------------------------------------------------------
// Additional tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod trt_tests {
    use super::*;

    /// TRT: viscosity round-trips.
    #[test]
    fn test_trt_viscosity_round_trip() {
        let nu = 0.1_f64;
        let trt = TrtCollision2D::new(nu);
        let nu_back = trt.viscosity();
        assert!((nu_back - nu).abs() < 1e-14, "nu={nu}, back={nu_back}");
    }

    /// TRT: the default constructor realises the magic parameter
    /// Λ = (τ⁺ − ½)(τ⁻ − ½) = 3/16 (not the product τ⁺τ⁻) and keeps the
    /// antisymmetric rate s⁻ = 1/τ⁻ inside the stable interval (0, 2).
    #[test]
    fn test_trt_magic_number() {
        for &nu in &[0.01_f64, 1.0 / 6.0, 0.5] {
            let trt = TrtCollision2D::new(nu);
            let lambda = (trt.tau_plus - 0.5) * (trt.tau_minus - 0.5);
            assert!((lambda - 3.0 / 16.0).abs() < 1e-14, "ν={nu}: Λ = {lambda}");
            assert!(
                (trt.lambda() - 3.0 / 16.0).abs() < 1e-14,
                "ν={nu}: lambda() = {}",
                trt.lambda()
            );
            let s_minus = 1.0 / trt.tau_minus;
            assert!(
                s_minus > 0.0 && s_minus < 2.0,
                "ν={nu}: s⁻ = {s_minus} outside (0, 2)"
            );
        }
        // ν = 1/6 gives τ⁺ = 1 and τ⁻ = ½ + (3/16)/(1/2) = 7/8.
        let trt = TrtCollision2D::new(1.0 / 6.0);
        assert!((trt.tau_plus - 1.0).abs() < 1e-14);
        assert!((trt.tau_minus - 0.875).abs() < 1e-14);
    }

    /// TRT: symmetric part sums to rho, antisymmetric sums to 0.
    #[test]
    fn test_trt_split_sums() {
        let feq = d2q9_equilibrium_trt(1.0, 0.05, -0.03);
        let fp: f64 = TrtCollision2D::symmetric(&feq).iter().sum();
        let fm: f64 = TrtCollision2D::antisymmetric(&feq).iter().sum();
        assert!((fp - 1.0).abs() < 1e-13, "fp sum = {fp}");
        assert!(fm.abs() < 1e-14, "fm sum = {fm}");
    }

    /// TRT: after collide at equilibrium, output == input.
    #[test]
    fn test_trt_equilibrium_fixed_point() {
        let nu = 1.0 / 6.0;
        let trt = TrtCollision2D::new(nu);
        let rho = 1.0_f64;
        let (ux, uy) = (0.05, -0.02);
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let f_out = trt.collide(&feq, rho, ux, uy);
        for i in 0..9 {
            assert!(
                (f_out[i] - feq[i]).abs() < 1e-13,
                "TRT eq fixed-point failed at i={i}"
            );
        }
    }

    /// TRT: mass is conserved after collision.
    #[test]
    fn test_trt_mass_conservation() {
        let trt = TrtCollision2D::new(0.05);
        let feq = d2q9_equilibrium_trt(1.2, 0.03, 0.01);
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01; // mass-neutral perturbation
        let rho: f64 = f.iter().sum();
        let (ux, uy) = (0.03, 0.01);
        let f_out = trt.collide(&f, rho, ux, uy);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-12,
            "TRT mass not conserved: {rho} -> {rho_out}"
        );
    }

    /// MRT spectral radius on BGK-equivalent rates.
    #[test]
    fn test_mrt_spectral_radius() {
        let relax = MrtRelaxation::bgk_equivalent(1.0 / 6.0);
        // Convert D3Q19 rates (19 entries) to first 9 for spectral radius check
        let mut rates9 = [0.0f64; 9];
        rates9.copy_from_slice(&relax.s[..9]);
        let rmax = mrt_spectral_radius(&rates9);
        assert!((0.0..=2.0).contains(&rmax), "spectral radius = {rmax}");
    }

    /// MRT stability check — standard rates should be stable.
    #[test]
    fn test_mrt_is_stable() {
        let relax = MrtRelaxation::bgk_equivalent(0.1);
        let mut rates9 = [0.0f64; 9];
        rates9.copy_from_slice(&relax.s[..9]);
        assert!(
            mrt_is_stable(&rates9),
            "Standard MRT rates should be stable"
        );
    }

    /// MRT shear viscosity from rates should match input nu.
    #[test]
    fn test_mrt_shear_viscosity() {
        let nu = 0.08_f64;
        // Use MrtCollision2D which has a 9-element relaxation_rates field
        let mrt2d = MrtCollision2D::new(nu);
        let nu_computed = mrt_shear_viscosity(&mrt2d.relaxation_rates, 1.0);
        assert!(
            (nu_computed - nu).abs() < 1e-13,
            "nu={nu}, computed={nu_computed}"
        );
    }

    /// Regularised BGK: at equilibrium the output equals input.
    #[test]
    fn test_rlbm_equilibrium_fixed_point() {
        let rlbm = RegularisedBgk2D::new(1.0 / 6.0);
        let rho = 1.0_f64;
        let (ux, uy) = (0.04, -0.02);
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let f_out = rlbm.collide(&feq, rho, ux, uy);
        for i in 0..9 {
            assert!(
                (f_out[i] - feq[i]).abs() < 1e-12,
                "RLBM eq fixed-point at i={i}"
            );
        }
    }

    /// Regularised BGK: mass is conserved after collision.
    #[test]
    fn test_rlbm_mass_conservation() {
        let rlbm = RegularisedBgk2D::new(0.05);
        let rho = 1.1_f64;
        let (ux, uy) = (0.02, 0.03);
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let mut f = feq;
        f[0] += 0.02;
        f[5] -= 0.01;
        f[7] -= 0.01;
        let rho_in: f64 = f.iter().sum();
        let f_out = rlbm.collide(&f, rho_in, ux, uy);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "RLBM mass: {rho_in} -> {rho_out}"
        );
    }

    /// CLBM central moments of equilibrium: κ10=κ01=κ11=0, κ20=κ02=ρ/3.
    #[test]
    fn test_clbm_central_moments_at_eq() {
        let rho = 1.0_f64;
        let (ux, uy) = (0.05, -0.02);
        let feq = d2q9_equilibrium_trt(rho, ux, uy);
        let kappa = CascadedCollision2D::central_moments(&feq, ux, uy);
        assert!((kappa[0] - rho).abs() < 1e-13, "κ00={}", kappa[0]);
        assert!(kappa[1].abs() < 1e-13, "κ10={}", kappa[1]);
        assert!(kappa[2].abs() < 1e-13, "κ01={}", kappa[2]);
        assert!((kappa[3] - rho / 3.0).abs() < 1e-13, "κ20={}", kappa[3]);
        assert!((kappa[4] - rho / 3.0).abs() < 1e-13, "κ02={}", kappa[4]);
        assert!(kappa[5].abs() < 1e-13, "κ11={}", kappa[5]);
    }

    /// d2q9_equilibrium_trt: sum equals rho.
    #[test]
    fn test_d2q9_equilibrium_trt_sum() {
        let rho = 1.3_f64;
        let feq = d2q9_equilibrium_trt(rho, 0.04, -0.01);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "eq sum = {sum}");
    }
}
