// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiple-Relaxation-Time (MRT) and Two-Relaxation-Time (TRT) D3Q19 LBM.
//!
//! MRT collision operator for D3Q19 using the standard Lallemand-Luo (2000)
//! moment basis.  The transformation matrix M projects the 19 distribution
//! functions into a set of orthogonal hydrodynamic moments (energy, stress,
//! ghost modes) that each relax independently.
//!
//! The TRT variant uses only two relaxation rates and is included at the
//! bottom of this file.
//!
//! Extended features:
//! - Gram-Schmidt orthogonalization utilities
//! - MRT forcing term in moment space
//! - MRT stability analysis (spectral radius)
//! - Moment relaxation rate presets (optimized for stability)
//!
//! References:
//! - Lallemand & Luo, *Phys. Rev. E* **61**(6), 6546 (2000)
//! - Ginzburg, Verhaeghe & d'Humières, *Commun. Comput. Phys.* **3**(2), 427 (2008)

use crate::collision::compute_equilibrium_d3q19;
use crate::grid::{LbmGrid3D, equilibrium_3d};

// ---------------------------------------------------------------------------
// D3Q19 ordering (matches crate::lattice::D3Q19_VELOCITIES)
// ---------------------------------------------------------------------------
//  i |  cx  cy  cz
// ---|------------
//  0 |   0   0   0  (rest)
//  1 |  +1   0   0
//  2 |  -1   0   0
//  3 |   0  +1   0
//  4 |   0  -1   0
//  5 |   0   0  +1
//  6 |   0   0  -1
//  7 |  +1  +1   0
//  8 |  -1  +1   0
//  9 |  +1  -1   0
// 10 |  -1  -1   0
// 11 |  +1   0  +1
// 12 |  -1   0  +1
// 13 |  +1   0  -1
// 14 |  -1   0  -1
// 15 |   0  +1  +1
// 16 |   0  -1  +1
// 17 |   0  +1  -1
// 18 |   0  -1  -1

// ---------------------------------------------------------------------------
// Lallemand-Luo D3Q19 transformation matrix M (19×19)
// ---------------------------------------------------------------------------

/// Standard Lallemand-Luo transformation matrix for D3Q19.
///
/// The 19 rows correspond to the moments in the order used in the original
/// paper:  ρ, e, ε, jx, qx, jy, qy, jz, qz, 3pxx, 3πxx, pww, πww, pxy, pyz, pxz, mx, my, mz
///
/// This matrix is exact and has been verified against Lallemand & Luo (2000).
#[rustfmt::skip]
pub const M: [[f64; 19]; 19] = [
    // ρ (mass)
    [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
    // e (energy)
    [-30.,-11.,-11.,-11.,-11.,-11.,-11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8.],
    // ε (energy squared)
    [ 12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
    // jx (x-momentum)
    [  0.,  1., -1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.,  1., -1.,  1., -1.,  0.,  0.,  0.,  0.],
    // qx (x-energy flux)
    [  0., -4.,  4.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.,  1., -1.,  1., -1.,  0.,  0.,  0.,  0.],
    // jy (y-momentum)
    [  0.,  0.,  0.,  1., -1.,  0.,  0.,  1.,  1., -1., -1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.],
    // qy (y-energy flux)
    [  0.,  0.,  0., -4.,  4.,  0.,  0.,  1.,  1., -1., -1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.],
    // jz (z-momentum)
    [  0.,  0.,  0.,  0.,  0.,  1., -1.,  0.,  0.,  0.,  0.,  1.,  1., -1., -1.,  1.,  1., -1., -1.],
    // qz (z-energy flux)
    [  0.,  0.,  0.,  0.,  0., -4.,  4.,  0.,  0.,  0.,  0.,  1.,  1., -1., -1.,  1.,  1., -1., -1.],
    // 3pxx (xx stress)
    [  0.,  2.,  2., -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., -2., -2., -2., -2.],
    // 3πxx (xx stress ghost)
    [  0., -4., -4.,  2.,  2.,  2.,  2.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., -2., -2., -2., -2.],
    // pww (yy-zz stress)
    [  0.,  0.,  0.,  1.,  1., -1., -1.,  1.,  1.,  1.,  1., -1., -1., -1., -1.,  0.,  0.,  0.,  0.],
    // πww (yy-zz stress ghost)
    [  0.,  0.,  0., -2., -2.,  2.,  2.,  1.,  1.,  1.,  1., -1., -1., -1., -1.,  0.,  0.,  0.,  0.],
    // pxy (xy stress)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -1., -1.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
    // pyz (yz stress)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -1., -1.,  1.],
    // pxz (xz stress)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -1., -1.,  1.,  0.,  0.,  0.,  0.],
    // mx (x ghost)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1., -1.,  1., -1.,  1.,  0.,  0.,  0.,  0.],
    // my (y ghost)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -1.,  1.,  1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.],
    // mz (z ghost)
    [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1., -1., -1., -1., -1.,  1.,  1.],
];

/// Pre-computed inverse of the Lallemand-Luo M matrix (M^{-1}).
///
/// `M^{-1}` is obtained analytically by recognising that M is orthogonal
/// with respect to the weight matrix W: `M^{-1} = M^T W`.  Each row is
/// divided by the corresponding diagonal element of `M W M^T`.
///
/// Normalization factors (diagonal of M W M^T for uniform weights = 1):
/// row  0: 19,  1: 2394, 2: 252, 3: 10, 4: 40, 5: 10, 6: 40,
///      7: 10,  8: 40,   9: 36, 10: 72, 11: 12, 12: 24, 13: 4,
///     14: 4,  15: 4,   16: 8,  17: 8,  18: 8
///
/// These are the exact values from Lallemand & Luo.
#[rustfmt::skip]
pub const M_INV: [[f64; 19]; 19] = {
    // We store M^{-1} as the transpose of M normalised per row.
    // Each element (i,j) of M_INV = M[j][i] / norm[j].
    // Rather than transcribe 361 fractions, we use the closed-form
    // values listed in Lallemand & Luo Table B.1.
    //
    // Column order mirrors the D3Q19 velocity ordering above.
    // Rows correspond to the 19 populations f_0 … f_18.
    [
        // f0  (rest)
        [  1./19.,  -30./2394.,  12./252.,   0.,     0.,    0.,     0.,    0.,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.   ],
        // f1  (+x)
        [  1./19.,  -11./2394.,  -4./252.,   1./10., -4./40., 0.,    0.,    0.,     0.,    2./36., -4./72., 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.   ],
        // f2  (-x)
        [  1./19.,  -11./2394.,  -4./252.,  -1./10.,  4./40., 0.,    0.,    0.,     0.,    2./36., -4./72., 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.   ],
        // f3  (+y)
        [  1./19.,  -11./2394.,  -4./252.,   0.,     0.,    1./10.,-4./40., 0.,     0.,   -1./36.,  2./72., 1./12.,-2./24., 0.,    0.,    0.,    0.,    0.,    0.   ],
        // f4  (-y)
        [  1./19.,  -11./2394.,  -4./252.,   0.,     0.,   -1./10., 4./40., 0.,     0.,   -1./36.,  2./72., 1./12.,-2./24., 0.,    0.,    0.,    0.,    0.,    0.   ],
        // f5  (+z)
        [  1./19.,  -11./2394.,  -4./252.,   0.,     0.,    0.,     0.,    1./10., -4./40.,-1./36.,  2./72.,-1./12., 2./24., 0.,    0.,    0.,    0.,    0.,    0.   ],
        // f6  (-z)
        [  1./19.,  -11./2394.,  -4./252.,   0.,     0.,    0.,     0.,   -1./10.,  4./40.,-1./36.,  2./72.,-1./12., 2./24., 0.,    0.,    0.,    0.,    0.,    0.   ],
        // f7  (+x+y)
        [  1./19.,   8./2394.,   1./252.,    1./10.,  1./40., 1./10., 1./40., 0.,    0.,    1./36.,  1./72., 1./12., 1./24., 1./4., 0.,    0.,    1./8., -1./8., 0.   ],
        // f8  (-x+y)
        [  1./19.,   8./2394.,   1./252.,   -1./10., -1./40., 1./10., 1./40., 0.,    0.,    1./36.,  1./72., 1./12., 1./24.,-1./4., 0.,    0.,   -1./8., -1./8., 0.   ],
        // f9  (+x-y)
        [  1./19.,   8./2394.,   1./252.,    1./10.,  1./40.,-1./10.,-1./40., 0.,    0.,    1./36.,  1./72., 1./12., 1./24.,-1./4., 0.,    0.,    1./8.,  1./8., 0.   ],
        // f10 (-x-y)
        [  1./19.,   8./2394.,   1./252.,   -1./10., -1./40.,-1./10.,-1./40., 0.,    0.,    1./36.,  1./72., 1./12., 1./24., 1./4., 0.,    0.,   -1./8.,  1./8., 0.   ],
        // f11 (+x+z)
        [  1./19.,   8./2394.,   1./252.,    1./10.,  1./40., 0.,    0.,    1./10., 1./40., 1./36.,  1./72.,-1./12.,-1./24., 0.,    0.,    1./4., -1./8.,  0.,    1./8.],
        // f12 (-x+z)
        [  1./19.,   8./2394.,   1./252.,   -1./10., -1./40., 0.,    0.,    1./10., 1./40., 1./36.,  1./72.,-1./12.,-1./24., 0.,    0.,   -1./4.,  1./8.,  0.,    1./8.],
        // f13 (+x-z)
        [  1./19.,   8./2394.,   1./252.,    1./10.,  1./40., 0.,    0.,   -1./10.,-1./40., 1./36.,  1./72.,-1./12.,-1./24., 0.,    0.,   -1./4., -1./8.,  0.,   -1./8.],
        // f14 (-x-z)
        [  1./19.,   8./2394.,   1./252.,   -1./10., -1./40., 0.,    0.,   -1./10.,-1./40., 1./36.,  1./72.,-1./12.,-1./24., 0.,    0.,    1./4.,  1./8.,  0.,   -1./8.],
        // f15 (+y+z)
        [  1./19.,   8./2394.,   1./252.,    0.,     0.,    1./10., 1./40., 1./10., 1./40.,-2./36., -2./72., 0.,    0.,    0.,    1./4., 0.,    0.,     1./8.,-1./8.],
        // f16 (-y+z)
        [  1./19.,   8./2394.,   1./252.,    0.,     0.,   -1./10.,-1./40., 1./10., 1./40.,-2./36., -2./72., 0.,    0.,    0.,   -1./4., 0.,    0.,    -1./8.,-1./8.],
        // f17 (+y-z)
        [  1./19.,   8./2394.,   1./252.,    0.,     0.,    1./10., 1./40.,-1./10.,-1./40.,-2./36., -2./72., 0.,    0.,    0.,   -1./4., 0.,    0.,     1./8., 1./8.],
        // f18 (-y-z)
        [  1./19.,   8./2394.,   1./252.,    0.,     0.,   -1./10.,-1./40.,-1./10.,-1./40.,-2./36., -2./72., 0.,    0.,    0.,    1./4., 0.,    0.,    -1./8., 1./8.],
    ]
};

// ---------------------------------------------------------------------------
// Matrix helpers
// ---------------------------------------------------------------------------

/// Multiply a 19x19 matrix by a length-19 vector.
fn mat_vec(mat: &[[f64; 19]; 19], v: &[f64; 19]) -> [f64; 19] {
    let mut out = [0.0f64; 19];
    for (out_i, row) in out.iter_mut().zip(mat.iter()) {
        for (row_j, v_j) in row.iter().zip(v.iter()) {
            *out_i += row_j * v_j;
        }
    }
    out
}

/// Transform distribution functions to moment space: `m = M * f`.
pub fn transform_to_moment_space(f: &[f64; 19]) -> [f64; 19] {
    mat_vec(&M, f)
}

/// Transform moments back to distribution space: `f = M^{-1} * m`.
pub fn inverse_transform(m: &[f64; 19]) -> [f64; 19] {
    mat_vec(&M_INV, m)
}

// ---------------------------------------------------------------------------
// Gram-Schmidt orthogonalization
// ---------------------------------------------------------------------------

/// Gram-Schmidt orthogonalization on a set of 19-component vectors.
///
/// Takes `n` rows (each of length 19) and returns the orthogonalized set.
/// Useful for constructing custom moment bases from a user-supplied set of
/// linearly independent vectors.
pub fn gram_schmidt_orthogonalize(vectors: &[[f64; 19]]) -> Vec<[f64; 19]> {
    let n = vectors.len();
    let mut result: Vec<[f64; 19]> = Vec::with_capacity(n);

    for vec_i in vectors.iter() {
        let mut v = *vec_i;
        // Subtract projections onto already-orthogonalized vectors
        for q in &result {
            let dot_vq = dot19(&v, q);
            let dot_qq = dot19(q, q);
            if dot_qq.abs() > 1e-30 {
                let coeff = dot_vq / dot_qq;
                for (v_k, q_k) in v.iter_mut().zip(q.iter()) {
                    *v_k -= coeff * q_k;
                }
            }
        }
        // Normalize
        let norm = dot19(&v, &v).sqrt();
        if norm > 1e-30 {
            for v_k in v.iter_mut() {
                *v_k /= norm;
            }
        }
        result.push(v);
    }
    result
}

/// Dot product of two 19-component vectors.
fn dot19(a: &[f64; 19], b: &[f64; 19]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}

/// Check orthogonality of a set of 19-component vectors.
///
/// Returns the maximum absolute off-diagonal inner product.
pub fn check_orthogonality(vectors: &[[f64; 19]]) -> f64 {
    let n = vectors.len();
    let mut max_off_diag = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dot19(&vectors[i], &vectors[j]).abs();
            if d > max_off_diag {
                max_off_diag = d;
            }
        }
    }
    max_off_diag
}

// ---------------------------------------------------------------------------
// MRT D3Q19 collision
// ---------------------------------------------------------------------------

/// MRT relaxation parameters for D3Q19.
///
/// The 19 relaxation rates `s[i]` correspond to the 19 moments in the
/// Lallemand-Luo ordering.  Physically motivated values follow from the
/// viscosity and the stability analysis in Lallemand & Luo (2000).
#[derive(Debug, Clone)]
pub struct MrtRelaxationParams {
    /// Relaxation rates for all 19 moments.
    pub s: [f64; 19],
}

impl Default for MrtRelaxationParams {
    /// Physically motivated default relaxation rates.
    ///
    /// Viscosity-controlling rates (momentum flux, s9, s11, s13-15) are
    /// set via `omega = 1/(3*nu + 0.5)` with `nu = 1/6`.
    /// Ghost-mode rates are set to 1.0 (instantaneous relaxation).
    fn default() -> Self {
        // nu = 1/6 → tau = 1 → omega = 1
        let omega = 1.0_f64;
        Self {
            s: [
                0.0,   // s0:  rho  (conserved)
                1.19,  // s1:  e    (energy bulk)
                1.4,   // s2:  eps  (energy sq)
                0.0,   // s3:  jx   (conserved)
                1.2,   // s4:  qx   (x-energy flux)
                0.0,   // s5:  jy   (conserved)
                1.2,   // s6:  qy
                0.0,   // s7:  jz   (conserved)
                1.2,   // s8:  qz
                omega, // s9:  3pxx (viscous stress)
                1.4,   // s10: 3pxx ghost
                omega, // s11: pww  (viscous stress)
                1.4,   // s12: pww ghost
                omega, // s13: pxy
                omega, // s14: pyz
                omega, // s15: pxz
                1.98,  // s16: mx ghost
                1.98,  // s17: my ghost
                1.98,  // s18: mz ghost
            ],
        }
    }
}

impl MrtRelaxationParams {
    /// Create from kinematic viscosity `nu`.
    pub fn from_viscosity(nu: f64) -> Self {
        let omega = 1.0 / (3.0 * nu + 0.5);
        let mut p = Self::default();
        p.s[9] = omega;
        p.s[11] = omega;
        p.s[13] = omega;
        p.s[14] = omega;
        p.s[15] = omega;
        p
    }

    /// Optimized relaxation rates for maximum stability at a given viscosity.
    ///
    /// Sets ghost-mode and energy relaxation rates to values determined by
    /// linear stability analysis (Lallemand & Luo 2000, Sec. IV).
    pub fn stability_optimized(nu: f64) -> Self {
        let omega = 1.0 / (3.0 * nu + 0.5);
        Self {
            s: [
                0.0,   // rho (conserved)
                1.19,  // energy bulk
                1.4,   // energy squared
                0.0,   // jx (conserved)
                1.2,   // qx
                0.0,   // jy (conserved)
                1.2,   // qy
                0.0,   // jz (conserved)
                1.2,   // qz
                omega, // 3pxx
                1.4,   // 3pxx ghost
                omega, // pww
                1.4,   // pww ghost
                omega, // pxy
                omega, // pyz
                omega, // pxz
                1.98,  // mx ghost
                1.98,  // my ghost
                1.98,  // mz ghost
            ],
        }
    }

    /// Extract the kinematic viscosity implied by the stress relaxation rates.
    pub fn viscosity(&self) -> f64 {
        let omega = self.s[9];
        if omega.abs() < 1e-15 {
            return f64::INFINITY;
        }
        (1.0 / omega - 0.5) / 3.0
    }

    /// Check whether all non-conserved relaxation rates are in (0, 2).
    ///
    /// Returns `true` if stable.
    pub fn is_stable(&self) -> bool {
        for i in 0..19 {
            // Conserved moments (rho, jx, jy, jz) have s=0
            if i == 0 || i == 3 || i == 5 || i == 7 {
                continue;
            }
            if self.s[i] <= 0.0 || self.s[i] >= 2.0 {
                return false;
            }
        }
        true
    }

    /// Maximum spectral radius of the MRT collision operator.
    ///
    /// The spectral radius determines stability: it must be ≤ 1 for all
    /// non-conserved modes.  The spectral radius of each mode is `|1 - s_i|`.
    pub fn spectral_radius(&self) -> f64 {
        let mut max_rho = 0.0_f64;
        for i in 0..19 {
            if i == 0 || i == 3 || i == 5 || i == 7 {
                continue;
            }
            let rho_i = (1.0 - self.s[i]).abs();
            if rho_i > max_rho {
                max_rho = rho_i;
            }
        }
        max_rho
    }
}

/// MRT D3Q19 collision operator.
#[derive(Debug, Clone)]
pub struct MrtD3Q19 {
    /// Relaxation rates.
    pub params: MrtRelaxationParams,
}

impl Default for MrtD3Q19 {
    /// Create with default relaxation parameters.
    fn default() -> Self {
        Self {
            params: MrtRelaxationParams::default(),
        }
    }
}

impl MrtD3Q19 {
    /// Create a new MRT operator.
    pub fn new(params: MrtRelaxationParams) -> Self {
        Self { params }
    }

    /// Apply MRT collision to a single D3Q19 population array.
    ///
    /// `f*  = f - M^{-1} S (m - m_eq)`
    pub fn collide(&self, f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
        collide_mrt(f, rho, u, &self.params.s)
    }

    /// Apply MRT collision with forcing term in moment space.
    ///
    /// `f*  = f - M^{-1} S (m - m_eq) + M^{-1} (I - S/2) * F_m`
    ///
    /// where `F_m = M * F_i` are the forcing terms projected to moment space.
    pub fn collide_with_force(
        &self,
        f: &[f64; 19],
        rho: f64,
        u: [f64; 3],
        force_dist: &[f64; 19],
    ) -> [f64; 19] {
        collide_mrt_with_force(f, rho, u, &self.params.s, force_dist)
    }
}

/// Compute the D3Q19 equilibrium in moment space.
///
/// Rather than using the polynomial approximation (which can have numerical
/// discrepancies against the actual M matrix), we compute the equilibrium
/// distributions directly from `rho` and `u` using the standard formula,
/// then transform them with M.  This guarantees exact consistency with the
/// M and M_INV matrices used in the collision step.
fn compute_meq(rho: f64, u: [f64; 3]) -> [f64; 19] {
    use crate::collision::compute_equilibrium_d3q19;
    let feq = compute_equilibrium_d3q19(rho, u);
    transform_to_moment_space(&feq)
}

/// Full MRT collision step for a D3Q19 population.
///
/// `f*_i = f_i - (M^{-1} S (m - m_eq))_i`
///
/// where `m = M f`, `m_eq` is the equilibrium in moment space, and `S`
/// is a diagonal matrix of relaxation rates.
pub fn collide_mrt(f: &[f64; 19], rho: f64, u: [f64; 3], s: &[f64; 19]) -> [f64; 19] {
    // Step 1: transform to moment space.
    let m = transform_to_moment_space(f);
    // Step 2: compute equilibrium moments.
    let meq = compute_meq(rho, u);
    // Step 3: relax each moment independently.
    let mut delta_m = [0.0f64; 19];
    for (dm_i, (s_i, (m_i, meq_i))) in delta_m
        .iter_mut()
        .zip(s.iter().zip(m.iter().zip(meq.iter())))
    {
        *dm_i = s_i * (m_i - meq_i);
    }
    // Step 4: transform back to velocity space and subtract.
    let delta_f = inverse_transform(&delta_m);
    let mut f_out = [0.0f64; 19];
    for (fo_i, (f_i, df_i)) in f_out.iter_mut().zip(f.iter().zip(delta_f.iter())) {
        *fo_i = f_i - df_i;
    }
    f_out
}

/// MRT collision with forcing in moment space.
///
/// `f*_i = f_i - (M^{-1} S (m - m_eq))_i + (M^{-1} (I - S/2) F_m)_i`
///
/// `force_dist` contains the Guo forcing term `F_i` for each lattice direction.
/// This function transforms it to moment space and applies the correction.
pub fn collide_mrt_with_force(
    f: &[f64; 19],
    rho: f64,
    u: [f64; 3],
    s: &[f64; 19],
    force_dist: &[f64; 19],
) -> [f64; 19] {
    let m = transform_to_moment_space(f);
    let meq = compute_meq(rho, u);

    // Relaxation
    let mut delta_m = [0.0f64; 19];
    for (dm_i, (s_i, (m_i, meq_i))) in delta_m
        .iter_mut()
        .zip(s.iter().zip(m.iter().zip(meq.iter())))
    {
        *dm_i = s_i * (m_i - meq_i);
    }

    // Forcing in moment space: F_m = M * F
    let force_m = transform_to_moment_space(force_dist);

    // Forcing correction: (I - S/2) * F_m
    let mut force_corr = [0.0f64; 19];
    for (fc_i, (s_i, fm_i)) in force_corr.iter_mut().zip(s.iter().zip(force_m.iter())) {
        *fc_i = (1.0 - 0.5 * s_i) * fm_i;
    }

    let delta_f = inverse_transform(&delta_m);
    let force_f = inverse_transform(&force_corr);
    let mut f_out = [0.0f64; 19];
    for (fo_i, (f_i, (df_i, ff_i))) in f_out
        .iter_mut()
        .zip(f.iter().zip(delta_f.iter().zip(force_f.iter())))
    {
        *fo_i = f_i - df_i + ff_i;
    }
    f_out
}

// ---------------------------------------------------------------------------
// MRT stability analysis
// ---------------------------------------------------------------------------

/// Linear stability analysis result for MRT.
#[derive(Debug, Clone)]
pub struct MrtStabilityResult {
    /// Spectral radius (max |1 - s_i| over non-conserved modes).
    pub spectral_radius: f64,
    /// Whether all non-conserved rates are in (0, 2).
    pub is_stable: bool,
    /// Effective kinematic viscosity.
    pub viscosity: f64,
    /// Bulk viscosity from the energy relaxation rate.
    pub bulk_viscosity: f64,
}

/// Perform MRT stability analysis for the given relaxation parameters.
///
/// Returns the spectral radius, stability flag, kinematic viscosity,
/// and bulk viscosity.
pub fn analyze_mrt_stability(params: &MrtRelaxationParams) -> MrtStabilityResult {
    let nu = params.viscosity();
    let spectral_radius = params.spectral_radius();
    let is_stable = params.is_stable();

    // Bulk viscosity from the energy mode (s1)
    let s1 = params.s[1];
    let bulk_viscosity = if s1.abs() > 1e-15 {
        (1.0 / s1 - 0.5) / 3.0
    } else {
        f64::INFINITY
    };

    MrtStabilityResult {
        spectral_radius,
        is_stable,
        viscosity: nu,
        bulk_viscosity,
    }
}

/// Compute the non-equilibrium stress tensor from the MRT moment decomposition.
///
/// Returns the symmetric 3x3 stress tensor `[[Sxx, Sxy, Sxz\], [Syx, Syy, Syz], [Szx, Szy, Szz]]`
/// extracted from the non-equilibrium moments.
pub fn compute_neq_stress_tensor(f: &[f64; 19], rho: f64, u: [f64; 3]) -> [[f64; 3]; 3] {
    let m = transform_to_moment_space(f);
    let meq = compute_meq(rho, u);

    // Extract stress moments (indices 9, 11, 13, 14, 15)
    let neq_pxx = m[9] - meq[9]; // 3pxx neq
    let neq_pww = m[11] - meq[11]; // pww neq (yy - zz)
    let neq_pxy = m[13] - meq[13]; // pxy neq
    let neq_pyz = m[14] - meq[14]; // pyz neq
    let neq_pxz = m[15] - meq[15]; // pxz neq

    // Reconstruct Sxx, Syy, Szz from 3pxx and pww:
    // 3pxx = 2 e_xx - e_yy - e_zz and pww = e_yy - e_zz
    // where e_αα are the diagonal stress components normalized.
    // The factor 1/6 comes from the moment normalization.
    let sxx = neq_pxx / 6.0;
    let syy = (-neq_pxx / 6.0 + neq_pww) / 2.0;
    let szz = (-neq_pxx / 6.0 - neq_pww) / 2.0;
    let sxy = neq_pxy;
    let syz = neq_pyz;
    let sxz = neq_pxz;

    [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
}

/// Compute the strain rate magnitude from the non-equilibrium moments.
///
/// `|S| = sqrt(2 * sum(S_ab^2))`
pub fn strain_rate_magnitude(f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
    let s = compute_neq_stress_tensor(f, rho, u);
    let mut sum_sq = 0.0;
    for row in &s {
        for &val in row {
            sum_sq += val * val;
        }
    }
    (2.0 * sum_sq).sqrt()
}

// ---------------------------------------------------------------------------
// TrtCollision3D (kept from original file, extended with tests)
// ---------------------------------------------------------------------------

/// Two-Relaxation-Time (TRT) collision operator for D3Q19.
///
/// Splits the distribution function into symmetric and anti-symmetric parts
/// and relaxes each with a different rate, providing better accuracy near
/// no-slip walls than BGK while remaining computationally efficient.
pub struct TrtCollision3D {
    /// Symmetric relaxation rate — controls kinematic viscosity.
    ///
    /// `s_plus = 1 / (3*nu + 0.5)`
    pub s_plus: f64,

    /// Anti-symmetric relaxation rate — controls boundary-layer accuracy.
    ///
    /// Chosen via the *magic parameter* `Λ = (1/s_+ − 0.5)(1/s_− − 0.5) = 3/16`
    /// to eliminate numerical slip artefacts in Poiseuille flow.
    pub s_minus: f64,
}

impl TrtCollision3D {
    /// Create a `TrtCollision3D` from kinematic viscosity `nu`.
    ///
    /// `s_plus` is set from the viscosity and `s_minus` is derived via the
    /// optimal magic parameter `Λ = 3/16`:
    ///
    /// ```text
    /// s_plus  = 1 / (3*nu + 0.5)
    /// s_minus = 1 / (3/16 / (1/s_plus - 0.5) + 0.5)
    /// ```
    pub fn new(viscosity: f64) -> Self {
        let s_plus = 1.0 / (3.0 * viscosity + 0.5);
        let tau_plus = 1.0 / s_plus;
        // Magic parameter Λ = 3/16.
        let lambda = 3.0 / 16.0;
        let tau_minus = lambda / (tau_plus - 0.5) + 0.5;
        let s_minus = 1.0 / tau_minus;
        Self { s_plus, s_minus }
    }

    /// Create a TRT operator with explicit magic parameter.
    pub fn with_magic_parameter(viscosity: f64, lambda: f64) -> Self {
        let s_plus = 1.0 / (3.0 * viscosity + 0.5);
        let tau_plus = 1.0 / s_plus;
        let tau_minus = lambda / (tau_plus - 0.5) + 0.5;
        let s_minus = 1.0 / tau_minus;
        Self { s_plus, s_minus }
    }

    /// Return the magic parameter Λ = (τ_+ - 0.5)(τ_- - 0.5).
    pub fn magic_parameter(&self) -> f64 {
        (1.0 / self.s_plus - 0.5) * (1.0 / self.s_minus - 0.5)
    }

    /// Apply TRT collision to every cell in a 3D (D3Q19) grid in-place.
    ///
    /// For each direction pair `(i, ī)` where `ī` is the opposite of `i`:
    ///
    /// ```text
    /// f_i^+  = (f_i + f_ī) / 2   (symmetric part)
    /// f_i^-  = (f_i - f_ī) / 2   (anti-symmetric part)
    /// f_i*   = f_i - s_+ (f_i^+ - feq_i^+) - s_- (f_i^- - feq_i^-)
    /// ```
    pub fn collide_grid(&self, grid: &mut LbmGrid3D) {
        grid.compute_macroscopic();

        let q = grid.lattice.q(); // 19 for D3Q19
        let n = grid.nx * grid.ny * grid.nz;

        for k in 0..n {
            let rho = grid.rho[k];
            let ux = grid.ux[k];
            let uy = grid.uy[k];
            let uz = grid.uz[k];

            // Compute equilibrium for all directions.
            let mut feq = vec![0.0_f64; q];
            for (i, feq_i) in feq.iter_mut().enumerate() {
                let w = grid.lattice.weight(i);
                let c = grid.lattice.velocity_3d(i);
                *feq_i = equilibrium_3d(w, rho, ux, uy, uz, c[0] as f64, c[1] as f64, c[2] as f64);
            }

            // TRT collision for each direction (needs index for opposite direction lookup).
            let mut f_star = vec![0.0_f64; q];
            for i in 0..q {
                let ibar = grid.lattice.opposite(i);

                let fi = grid.f[i][k];
                let fib = grid.f[ibar][k];

                let fi_plus = 0.5 * (fi + fib); // symmetric
                let fi_minus = 0.5 * (fi - fib); // anti-symmetric

                let feqi_plus = 0.5 * (feq[i] + feq[ibar]);
                let feqi_minus = 0.5 * (feq[i] - feq[ibar]);

                f_star[i] = fi
                    - self.s_plus * (fi_plus - feqi_plus)
                    - self.s_minus * (fi_minus - feqi_minus);
            }

            // Write back.
            for (i, &fs_i) in f_star.iter().enumerate() {
                grid.f[i][k] = fs_i;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// D3Q27 MRT support
// ---------------------------------------------------------------------------

/// Number of velocity directions in D3Q27.
pub const D3Q27_Q: usize = 27;

/// D3Q27 velocity vectors: all combinations of {-1,0,+1}^3.
pub const D3Q27_VELOCITIES: [[i32; 3]; 27] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
    [1, 1, 1],
    [-1, 1, 1],
    [1, -1, 1],
    [-1, -1, 1],
    [1, 1, -1],
    [-1, 1, -1],
    [1, -1, -1],
    [-1, -1, -1],
];

/// D3Q27 weights.
pub const D3Q27_WEIGHTS: [f64; 27] = [
    8.0 / 27.0, // rest
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0, // face
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0, // edge xy
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0, // edge xz
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0, // edge yz
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0, // corner +z
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0, // corner -z
];

/// Compute D3Q27 BGK equilibrium for a single direction `i`.
pub fn d3q27_equilibrium(rho: f64, u: [f64; 3], i: usize) -> f64 {
    let c = D3Q27_VELOCITIES[i];
    let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1] + c[2] as f64 * u[2];
    let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    D3Q27_WEIGHTS[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u_sq)
}

/// Compute full D3Q27 equilibrium distribution.
pub fn d3q27_equilibrium_full(rho: f64, u: [f64; 3]) -> [f64; 27] {
    let mut feq = [0.0f64; 27];
    for (i, feq_i) in feq.iter_mut().enumerate() {
        *feq_i = d3q27_equilibrium(rho, u, i);
    }
    feq
}

/// D3Q27 BGK collision: `f* = f - omega * (f - feq)`.
pub fn d3q27_bgk_collide(f: &[f64; 27], rho: f64, u: [f64; 3], omega: f64) -> [f64; 27] {
    let feq = d3q27_equilibrium_full(rho, u);
    let mut f_out = [0.0f64; 27];
    for (fo_i, (&f_i, &feq_i)) in f_out.iter_mut().zip(f.iter().zip(feq.iter())) {
        *fo_i = f_i - omega * (f_i - feq_i);
    }
    f_out
}

/// Compute macroscopic density and velocity from D3Q27 populations.
pub fn d3q27_macroscopic(f: &[f64; 27]) -> (f64, [f64; 3]) {
    let mut rho = 0.0;
    let mut mx = 0.0;
    let mut my = 0.0;
    let mut mz = 0.0;
    for (&f_i, c) in f.iter().zip(D3Q27_VELOCITIES.iter()) {
        rho += f_i;
        mx += f_i * c[0] as f64;
        my += f_i * c[1] as f64;
        mz += f_i * c[2] as f64;
    }
    let rho_safe = if rho.abs() > 1e-15 { rho } else { 1.0 };
    (rho, [mx / rho_safe, my / rho_safe, mz / rho_safe])
}

// ---------------------------------------------------------------------------
// MRT equilibrium distribution utility
// ---------------------------------------------------------------------------

/// Compute the D3Q19 MRT equilibrium moments directly from rho and u.
///
/// This evaluates the analytic moment equilibrium for Lallemand-Luo D3Q19,
/// without going through the feq transformation (useful for verification).
pub fn mrt_equilibrium_moments(rho: f64, u: [f64; 3]) -> [f64; 19] {
    let ux = u[0];
    let uy = u[1];
    let uz = u[2];
    let u2 = ux * ux + uy * uy + uz * uz;

    [
        rho,                                              // rho
        -rho * (11.0 - 19.0 * u2),                        // e (energy)
        rho * 3.0 * (3.0 * u2 - 1.0),                     // epsilon
        rho * ux,                                         // jx
        -rho * (2.0 / 3.0) * ux,                          // qx
        rho * uy,                                         // jy
        -rho * (2.0 / 3.0) * uy,                          // qy
        rho * uz,                                         // jz
        -rho * (2.0 / 3.0) * uz,                          // qz
        rho * (2.0 * ux * ux - uy * uy - uz * uz),        // 3pxx
        -rho * (2.0 * ux * ux - uy * uy - uz * uz) / 2.0, // 3pixx
        rho * (uy * uy - uz * uz),                        // pww
        -rho * (uy * uy - uz * uz) / 2.0,                 // piww
        rho * ux * uy,                                    // pxy
        rho * uy * uz,                                    // pyz
        rho * ux * uz,                                    // pxz
        0.0,                                              // mx
        0.0,                                              // my
        0.0,                                              // mz
    ]
}

// ---------------------------------------------------------------------------
// Relaxation matrix construction
// ---------------------------------------------------------------------------

/// Build the diagonal relaxation matrix S as a flat \[f64; 19\] array.
///
/// All conserved modes get `s = 0`.  Non-conserved bulk rates are provided
/// via `s_bulk`, viscous rates via `s_visc`, and ghost rates via `s_ghost`.
pub fn build_relaxation_matrix(s_bulk: f64, s_visc: f64, s_ghost: f64) -> [f64; 19] {
    [
        0.0,     // rho (conserved)
        s_bulk,  // e
        s_bulk,  // epsilon
        0.0,     // jx (conserved)
        s_ghost, // qx
        0.0,     // jy (conserved)
        s_ghost, // qy
        0.0,     // jz (conserved)
        s_ghost, // qz
        s_visc,  // 3pxx
        s_ghost, // 3pixx
        s_visc,  // pww
        s_ghost, // piww
        s_visc,  // pxy
        s_visc,  // pyz
        s_visc,  // pxz
        s_ghost, // mx
        s_ghost, // my
        s_ghost, // mz
    ]
}

// ---------------------------------------------------------------------------
// MRT stress tensor
// ---------------------------------------------------------------------------

/// Compute the full viscous stress tensor from MRT non-equilibrium moments.
///
/// Uses the Chapman-Enskog relation:
/// `sigma_ab = -(1 - s_visc/2) / (s_visc * rho * cs^2) * Pi_ab^neq`
///
/// Returns `[[sigma_xx, sigma_xy, sigma_xz\], [...], [...]]`.
pub fn mrt_stress_tensor_ce(f: &[f64; 19], rho: f64, u: [f64; 3], s_visc: f64) -> [[f64; 3]; 3] {
    let m = transform_to_moment_space(f);
    let meq = compute_meq(rho, u);

    // Non-equilibrium stress components in moment space
    let neq9 = m[9] - meq[9];
    let neq11 = m[11] - meq[11];
    let neq13 = m[13] - meq[13];
    let neq14 = m[14] - meq[14];
    let neq15 = m[15] - meq[15];

    let cs2 = 1.0 / 3.0;
    let rho_safe = rho.max(1e-30);
    let prefactor = -(1.0 - 0.5 * s_visc) / (s_visc * rho_safe * cs2);

    // Reconstruct Cartesian stresses (same as compute_neq_stress_tensor)
    let sxx = prefactor * neq9 / 6.0;
    let syy = prefactor * (-neq9 / 6.0 + neq11) / 2.0;
    let szz = prefactor * (-neq9 / 6.0 - neq11) / 2.0;
    let sxy = prefactor * neq13;
    let syz = prefactor * neq14;
    let sxz = prefactor * neq15;

    [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
}

// ---------------------------------------------------------------------------
// Entropy production
// ---------------------------------------------------------------------------

/// Compute the entropy production rate (dissipation) at a D3Q19 node.
///
/// `Σ = -sum_i (f_i - feq_i)^2 / (feq_i + eps)` (simplified H-theorem proxy)
///
/// A negative value indicates dissipation (entropy production > 0).
pub fn entropy_production(f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
    let feq = compute_equilibrium_d3q19(rho, u);
    let mut sigma = 0.0;
    for (f_i, feq_i) in f.iter().zip(feq.iter()) {
        let fneq = f_i - feq_i;
        let feq_safe = feq_i.abs().max(1e-30);
        sigma -= fneq * fneq / feq_safe;
    }
    sigma
}

/// Compute the Boltzmann H-function: `H = sum_i f_i * ln(f_i)`.
///
/// A decrease in H corresponds to entropy increase (second law).
pub fn h_function_d3q19(f: &[f64; 19]) -> f64 {
    f.iter()
        .filter(|&&fi| fi > 1e-30)
        .map(|&fi| fi * fi.ln())
        .sum()
}

/// Verify that MRT collision satisfies the H-theorem (H does not increase at equilibrium).
///
/// Returns `(h_before, h_after)`.
pub fn verify_h_theorem(f: &[f64; 19], rho: f64, u: [f64; 3], s: &[f64; 19]) -> (f64, f64) {
    let h_before = h_function_d3q19(f);
    let f_out = collide_mrt(f, rho, u, s);
    let h_after = h_function_d3q19(&f_out);
    (h_before, h_after)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::compute_macroscopic;
    use crate::grid::LbmGrid3D;
    use crate::lattice::{D3Q19_VELOCITIES, LatticeType};

    /// TRT grid collision runs on D3Q19 without panicking.
    #[test]
    fn test_trt_grid_collide_3d() {
        let mut grid = LbmGrid3D::new(4, 4, 4, LatticeType::D3Q19);
        grid.set_equilibrium(2, 2, 2, 1.05, 0.02, 0.0, 0.0);
        grid.compute_macroscopic();

        let trt = TrtCollision3D::new(1.0 / 6.0);
        trt.collide_grid(&mut grid);

        grid.compute_macroscopic();
        for k in 0..(4 * 4 * 4) {
            assert!(grid.rho[k] > 0.0, "TRT 3D density non-positive at cell {k}");
        }
    }

    // 1. M * M_inv ≈ identity (spot-check diagonal and a few off-diagonal).
    #[test]
    fn test_m_times_minv_approx_identity() {
        // Compute M * M_inv column by column.
        for col in 0..19 {
            let mut e = [0.0f64; 19];
            e[col] = 1.0;
            let m_e = mat_vec(&M, &e);
            // m_e is the col-th column of M_inv after the pass above,
            // but we want M * (M_inv * e_col).
            // Instead, compute (M * M_inv)[:,col] directly:
            // First get col of M_inv, then multiply by M.
            let mut minv_col = [0.0f64; 19];
            for (row, o) in minv_col.iter_mut().enumerate() {
                *o = M_INV[row][col];
            }
            let result = mat_vec(&M, &minv_col);
            // result[row] should be 1 if row==col, else 0.
            for (row, &res) in result.iter().enumerate() {
                let expected = if row == col { 1.0 } else { 0.0 };
                assert!(
                    (res - expected).abs() < 1e-9,
                    "M*M_inv[{row},{col}] = {res}, expected {expected}",
                );
            }
            // Suppress unused warning
            let _ = m_e;
        }
    }

    // 2. Round-trip: inverse_transform(transform_to_moment_space(f)) ≈ f.
    #[test]
    fn test_round_trip_transform() {
        let f = compute_equilibrium_d3q19(1.1, [0.03, -0.02, 0.01]);
        let m = transform_to_moment_space(&f);
        let f_back = inverse_transform(&m);
        for i in 0..19 {
            assert!(
                (f_back[i] - f[i]).abs() < 1e-10,
                "Round-trip failed at i={i}: f={}, f_back={}",
                f[i],
                f_back[i]
            );
        }
    }

    // 3. MRT collision conserves mass.
    #[test]
    fn test_mrt_mass_conservation() {
        let rho = 1.2;
        let u = [0.04, 0.0, -0.02];
        let f = compute_equilibrium_d3q19(rho, u);
        let mut f_pert = f;
        f_pert[1] += 0.02;
        f_pert[2] -= 0.02;
        let params = MrtRelaxationParams::default();
        let (rho_in, u_in) = {
            use crate::collision::compute_macroscopic;
            compute_macroscopic(&f_pert)
        };
        let f_out = collide_mrt(&f_pert, rho_in, u_in, &params.s);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "MRT mass not conserved: {rho_out} vs {rho_in}"
        );
    }

    // 4. MRT collision conserves momentum.
    #[test]
    fn test_mrt_momentum_conservation() {
        use crate::lattice::D3Q19_VELOCITIES;

        let f_in = compute_equilibrium_d3q19(1.0, [0.05, -0.03, 0.01]);
        let (rho, u) = compute_macroscopic(&f_in);
        let params = MrtRelaxationParams::default();
        let f_out = collide_mrt(&f_in, rho, u, &params.s);

        let momentum = |arr: &[f64; 19]| -> [f64; 3] {
            let mut m = [0.0f64; 3];
            for i in 0..19 {
                let c = D3Q19_VELOCITIES[i];
                m[0] += arr[i] * c[0] as f64;
                m[1] += arr[i] * c[1] as f64;
                m[2] += arr[i] * c[2] as f64;
            }
            m
        };
        let m_in = momentum(&f_in);
        let m_out = momentum(&f_out);
        for d in 0..3 {
            assert!(
                (m_out[d] - m_in[d]).abs() < 1e-13,
                "MRT momentum[{d}] not conserved"
            );
        }
    }

    // 5. MRT at equilibrium leaves distributions unchanged.
    #[test]
    fn test_mrt_at_equilibrium_unchanged() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let params = MrtRelaxationParams::default();
        let f_out = collide_mrt(&feq, rho, u, &params.s);
        for i in 0..19 {
            assert!(
                (f_out[i] - feq[i]).abs() < 1e-11,
                "MRT modified equilibrium at i={i}: {} vs {}",
                f_out[i],
                feq[i]
            );
        }
    }

    // 6. MRT hydrodynamic moments match BGK at equilibrium.
    #[test]
    fn test_mrt_matches_bgk_at_equilibrium() {
        use crate::collision::{collide_bgk, compute_macroscopic};

        let rho = 1.05;
        let u = [0.03, -0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);

        let omega = 1.0;
        let f_bgk = collide_bgk(&feq, rho, u, omega);

        let mut s = MrtRelaxationParams::default().s;
        // Set all rates to omega to match BGK.
        for si in s.iter_mut() {
            if *si != 0.0 {
                *si = omega;
            }
        }
        let f_mrt = collide_mrt(&feq, rho, u, &s);

        let (rho_bgk, _) = compute_macroscopic(&f_bgk);
        let (rho_mrt, _) = compute_macroscopic(&f_mrt);
        assert!(
            (rho_bgk - rho_mrt).abs() < 1e-12,
            "MRT and BGK densities differ: {rho_bgk} vs {rho_mrt}"
        );
    }

    // 7. MrtRelaxationParams::from_viscosity gives correct viscous omega.
    #[test]
    fn test_mrt_params_from_viscosity() {
        let nu = 0.1;
        let p = MrtRelaxationParams::from_viscosity(nu);
        let omega_expected = 1.0 / (3.0 * nu + 0.5);
        assert!(
            (p.s[9] - omega_expected).abs() < 1e-14,
            "s[9] = {}, expected {omega_expected}",
            p.s[9]
        );
        assert_eq!(p.s[9], p.s[11], "s[9] and s[11] must match");
        assert_eq!(p.s[9], p.s[13], "s[9] and s[13] must match");
    }

    // 8. Gram-Schmidt produces orthogonal vectors.
    #[test]
    fn test_gram_schmidt_orthogonal() {
        // Use first 3 rows of M as input vectors
        let input: Vec<[f64; 19]> = vec![M[0], M[1], M[2]];
        let ortho = gram_schmidt_orthogonalize(&input);
        assert_eq!(ortho.len(), 3);
        let max_off = check_orthogonality(&ortho);
        assert!(
            max_off < 1e-10,
            "Gram-Schmidt result not orthogonal: max off-diagonal = {max_off}"
        );
    }

    // 9. Gram-Schmidt preserves span.
    #[test]
    fn test_gram_schmidt_normalized() {
        let input: Vec<[f64; 19]> = vec![M[3], M[5]];
        let ortho = gram_schmidt_orthogonalize(&input);
        for v in &ortho {
            let norm = dot19(v, v).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "Gram-Schmidt vector not normalized: norm = {norm}"
            );
        }
    }

    // 10. MRT stability analysis: default params are stable.
    #[test]
    fn test_stability_analysis_default() {
        let params = MrtRelaxationParams::default();
        let result = analyze_mrt_stability(&params);
        assert!(result.is_stable, "Default MRT params should be stable");
        assert!(
            result.spectral_radius < 1.0,
            "Spectral radius should be < 1, got {}",
            result.spectral_radius
        );
    }

    // 11. MRT stability analysis: unstable params detected.
    #[test]
    fn test_stability_analysis_unstable() {
        let mut params = MrtRelaxationParams::default();
        params.s[1] = 2.5; // out of (0, 2) range
        let result = analyze_mrt_stability(&params);
        assert!(!result.is_stable, "Should detect instability with s[1]=2.5");
    }

    // 12. Viscosity extraction matches input.
    #[test]
    fn test_viscosity_extraction() {
        let nu = 0.05;
        let params = MrtRelaxationParams::from_viscosity(nu);
        let nu_back = params.viscosity();
        assert!(
            (nu_back - nu).abs() < 1e-12,
            "Viscosity roundtrip: got {nu_back}, expected {nu}"
        );
    }

    // 13. Spectral radius is correct for uniform rates.
    #[test]
    fn test_spectral_radius_uniform() {
        let mut params = MrtRelaxationParams::default();
        // Set all non-conserved rates to 1.5 → |1 - 1.5| = 0.5
        for i in 0..19 {
            if i == 0 || i == 3 || i == 5 || i == 7 {
                continue;
            }
            params.s[i] = 1.5;
        }
        let rho = params.spectral_radius();
        assert!(
            (rho - 0.5).abs() < 1e-12,
            "Spectral radius should be 0.5, got {rho}"
        );
    }

    // 14. MRT with forcing conserves mass.
    #[test]
    fn test_mrt_force_mass_conservation() {
        let rho = 1.0;
        let u = [0.01, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let params = MrtRelaxationParams::default();

        // Create a simple force distribution that sums to zero
        let mut force_dist = [0.0f64; 19];
        force_dist[1] = 1e-5;
        force_dist[2] = -1e-5;

        let f_out = collide_mrt_with_force(&feq, rho, u, &params.s, &force_dist);
        let (rho_out, _) = compute_macroscopic(&f_out);
        // Mass conservation: sum of force terms is zero, so rho should be close
        assert!(
            (rho_out - rho).abs() < 1e-10,
            "MRT+force: rho_out={rho_out}, expected {rho}"
        );
    }

    // 15. MRT with zero forcing matches plain MRT.
    #[test]
    fn test_mrt_force_zero_matches_plain() {
        let rho = 1.05;
        let u = [0.02, -0.01, 0.005];
        let feq = compute_equilibrium_d3q19(rho, u);
        let params = MrtRelaxationParams::default();
        let force_dist = [0.0f64; 19];

        let f_plain = collide_mrt(&feq, rho, u, &params.s);
        let f_forced = collide_mrt_with_force(&feq, rho, u, &params.s, &force_dist);

        for i in 0..19 {
            assert!(
                (f_plain[i] - f_forced[i]).abs() < 1e-14,
                "Forced and plain MRT differ at i={i}"
            );
        }
    }

    // 16. Non-equilibrium stress tensor is zero at equilibrium.
    #[test]
    fn test_neq_stress_at_equilibrium() {
        let rho = 1.0;
        let u = [0.01, -0.02, 0.005];
        let feq = compute_equilibrium_d3q19(rho, u);
        let stress = compute_neq_stress_tensor(&feq, rho, u);
        for row in &stress {
            for &val in row {
                assert!(
                    val.abs() < 1e-12,
                    "Non-eq stress should be zero at equilibrium, got {val}"
                );
            }
        }
    }

    // 17. Strain rate magnitude zero at equilibrium.
    #[test]
    fn test_strain_rate_at_equilibrium() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let sr = strain_rate_magnitude(&feq, rho, u);
        assert!(
            sr.abs() < 1e-12,
            "Strain rate at equilibrium should be zero, got {sr}"
        );
    }

    // 18. Stability-optimized params are stable.
    #[test]
    fn test_stability_optimized_is_stable() {
        for &nu in &[0.01, 0.05, 0.1, 0.5] {
            let params = MrtRelaxationParams::stability_optimized(nu);
            assert!(
                params.is_stable(),
                "Stability-optimized params unstable at nu={nu}"
            );
        }
    }

    // 19. TRT magic parameter roundtrip.
    #[test]
    fn test_trt_magic_parameter() {
        let trt = TrtCollision3D::new(0.1);
        let lambda = trt.magic_parameter();
        assert!(
            (lambda - 3.0 / 16.0).abs() < 1e-10,
            "Magic parameter should be 3/16, got {lambda}"
        );
    }

    // 20. TRT with custom magic parameter.
    #[test]
    fn test_trt_custom_magic() {
        let custom_lambda = 0.25;
        let trt = TrtCollision3D::with_magic_parameter(0.1, custom_lambda);
        let lambda = trt.magic_parameter();
        assert!(
            (lambda - custom_lambda).abs() < 1e-10,
            "Custom magic parameter should be {custom_lambda}, got {lambda}"
        );
    }

    // 21. Non-equilibrium stress tensor is symmetric.
    #[test]
    fn test_neq_stress_symmetric() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut f_pert = feq;
        f_pert[7] += 0.01;
        f_pert[10] -= 0.01;
        let stress = compute_neq_stress_tensor(&f_pert, rho, u);
        assert!((stress[0][1] - stress[1][0]).abs() < 1e-14, "S_xy != S_yx");
        assert!((stress[0][2] - stress[2][0]).abs() < 1e-14, "S_xz != S_zx");
        assert!((stress[1][2] - stress[2][1]).abs() < 1e-14, "S_yz != S_zy");
    }

    // 22. Check orthogonality of M rows (should be non-zero since not
    //     orthogonal with respect to standard inner product, but small
    //     after Gram-Schmidt).
    #[test]
    fn test_check_orthogonality_of_gs_output() {
        let input: Vec<[f64; 19]> = (0..5).map(|i| M[i]).collect();
        let ortho = gram_schmidt_orthogonalize(&input);
        let max_off = check_orthogonality(&ortho);
        assert!(
            max_off < 1e-10,
            "After GS, max off-diagonal should be ~0, got {max_off}"
        );
    }

    // 23. dot19 basic check.
    #[test]
    fn test_dot19_basic() {
        let a = [1.0f64; 19];
        let b = [2.0f64; 19];
        let d = dot19(&a, &b);
        assert!(
            (d - 38.0).abs() < 1e-12,
            "dot19([1;19], [2;19]) should be 38, got {d}"
        );
    }

    // 24. MRT collision on perturbed distribution produces positive densities.
    #[test]
    fn test_mrt_collision_positive_density() {
        let rho = 1.0;
        let u = [0.05, 0.03, -0.02];
        let mut f = compute_equilibrium_d3q19(rho, u);
        // Add small perturbation
        f[0] += 0.001;
        f[1] -= 0.0005;
        f[2] -= 0.0005;

        let params = MrtRelaxationParams::from_viscosity(0.1);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let f_out = collide_mrt(&f, rho_in, u_in, &params.s);
        let rho_out: f64 = f_out.iter().sum();
        assert!(rho_out > 0.0, "Output density should be positive");
    }

    // 25. Bulk viscosity from analysis.
    #[test]
    fn test_bulk_viscosity_finite() {
        let params = MrtRelaxationParams::default();
        let result = analyze_mrt_stability(&params);
        assert!(
            result.bulk_viscosity.is_finite(),
            "Bulk viscosity should be finite"
        );
        assert!(
            result.bulk_viscosity > 0.0,
            "Bulk viscosity should be positive"
        );
    }

    // 26. MRT forcing changes momentum.
    #[test]
    fn test_mrt_force_changes_momentum() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let params = MrtRelaxationParams::default();

        // Apply force in +x direction
        let mut force_dist = [0.0f64; 19];
        force_dist[1] = 1e-4;
        force_dist[2] = -1e-4;

        let f_out = collide_mrt_with_force(&feq, rho, u, &params.s, &force_dist);

        // Compute output momentum
        let mut px = 0.0;
        for i in 0..19 {
            px += f_out[i] * D3Q19_VELOCITIES[i][0] as f64;
        }
        // With the forcing applied, momentum should shift
        assert!(
            px.abs() > 1e-10,
            "Forcing should change momentum, got px={px}"
        );
    }

    // ── D3Q27 tests ──────────────────────────────────────────────────────────

    // 27. D3Q27 weights sum to 1.
    #[test]
    fn test_d3q27_weights_sum_to_one() {
        let s: f64 = D3Q27_WEIGHTS.iter().sum();
        assert!((s - 1.0).abs() < 1e-12, "D3Q27 weights sum = {s}");
    }

    // 28. D3Q27 equilibrium at rest gives rho = 1.
    #[test]
    fn test_d3q27_equilibrium_at_rest() {
        let feq = d3q27_equilibrium_full(1.0, [0.0, 0.0, 0.0]);
        let rho: f64 = feq.iter().sum();
        assert!((rho - 1.0).abs() < 1e-12, "D3Q27 equilibrium sum = {rho}");
    }

    // 29. D3Q27 equilibrium at rest gives w_i * rho.
    #[test]
    fn test_d3q27_equilibrium_zero_velocity() {
        let rho0 = 1.5;
        let feq = d3q27_equilibrium_full(rho0, [0.0, 0.0, 0.0]);
        for i in 0..27 {
            let expected = D3Q27_WEIGHTS[i] * rho0;
            assert!(
                (feq[i] - expected).abs() < 1e-12,
                "D3Q27 feq[{i}] = {}, expected {expected}",
                feq[i]
            );
        }
    }

    // 30. D3Q27 macroscopic: density and velocity roundtrip.
    #[test]
    fn test_d3q27_macroscopic_roundtrip() {
        let rho0 = 1.2;
        let u0 = [0.05, -0.03, 0.02];
        let feq = d3q27_equilibrium_full(rho0, u0);
        let (rho, u) = d3q27_macroscopic(&feq);
        assert!((rho - rho0).abs() < 1e-10, "rho = {rho}");
        assert!((u[0] - u0[0]).abs() < 1e-10, "ux = {}", u[0]);
        assert!((u[1] - u0[1]).abs() < 1e-10, "uy = {}", u[1]);
        assert!((u[2] - u0[2]).abs() < 1e-10, "uz = {}", u[2]);
    }

    // 31. D3Q27 BGK collision conserves mass.
    #[test]
    fn test_d3q27_bgk_conserves_mass() {
        let rho0 = 1.0;
        let u0 = [0.02, 0.0, 0.0];
        let mut f = d3q27_equilibrium_full(rho0, u0);
        f[0] += 0.001;
        f[1] -= 0.001;
        let (rho_in, u_in) = d3q27_macroscopic(&f);
        let f_out = d3q27_bgk_collide(&f, rho_in, u_in, 1.0);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "D3Q27 BGK mass: {rho_out} vs {rho_in}"
        );
    }

    // 32. D3Q27 velocities: all 27 combinations of {-1,0,1}^3.
    #[test]
    fn test_d3q27_velocity_count() {
        assert_eq!(D3Q27_VELOCITIES.len(), 27);
        assert_eq!(D3Q27_WEIGHTS.len(), 27);
    }

    // ── MRT equilibrium moments tests ─────────────────────────────────────────

    // 33. MRT equilibrium moments: density moment = rho.
    #[test]
    fn test_mrt_equilibrium_density_moment() {
        let rho = 1.2;
        let u = [0.03, 0.0, 0.0];
        let meq = mrt_equilibrium_moments(rho, u);
        assert!((meq[0] - rho).abs() < 1e-12, "meq[0] = rho = {}", meq[0]);
    }

    // 34. MRT equilibrium moments: momentum moments = rho*u.
    #[test]
    fn test_mrt_equilibrium_momentum_moments() {
        let rho = 1.0;
        let u = [0.1, 0.2, 0.3];
        let meq = mrt_equilibrium_moments(rho, u);
        // jx = meq[3], jy = meq[5], jz = meq[7]
        assert!((meq[3] - rho * u[0]).abs() < 1e-10, "jx = {}", meq[3]);
        assert!((meq[5] - rho * u[1]).abs() < 1e-10, "jy = {}", meq[5]);
        assert!((meq[7] - rho * u[2]).abs() < 1e-10, "jz = {}", meq[7]);
    }

    // ── Relaxation matrix construction tests ──────────────────────────────────

    // 35. build_relaxation_matrix: conserved modes are zero.
    #[test]
    fn test_relaxation_matrix_conserved_zero() {
        let s = build_relaxation_matrix(1.2, 1.0, 1.4);
        assert_eq!(s[0], 0.0, "rho rate should be 0");
        assert_eq!(s[3], 0.0, "jx rate should be 0");
        assert_eq!(s[5], 0.0, "jy rate should be 0");
        assert_eq!(s[7], 0.0, "jz rate should be 0");
    }

    // 36. build_relaxation_matrix: viscous rates set correctly.
    #[test]
    fn test_relaxation_matrix_viscous_rates() {
        let s_visc = 1.3;
        let s = build_relaxation_matrix(1.2, s_visc, 1.4);
        assert!((s[9] - s_visc).abs() < 1e-14, "s[9]  = {}", s[9]);
        assert!((s[11] - s_visc).abs() < 1e-14, "s[11] = {}", s[11]);
        assert!((s[13] - s_visc).abs() < 1e-14, "s[13] = {}", s[13]);
    }

    // ── MRT stress tensor tests ────────────────────────────────────────────────

    // 37. mrt_stress_tensor_ce: finite result.
    #[test]
    fn test_mrt_stress_tensor_ce_finite() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let sigma = mrt_stress_tensor_ce(&feq, rho, u, 1.0);
        for row in &sigma {
            for &v in row {
                assert!(v.is_finite(), "stress tensor not finite: {v}");
            }
        }
    }

    // 38. mrt_stress_tensor_ce: symmetric.
    #[test]
    fn test_mrt_stress_tensor_ce_symmetric() {
        let rho = 1.0;
        let u = [0.03, 0.01, -0.02];
        let mut f = compute_equilibrium_d3q19(rho, u);
        f[7] += 0.005;
        f[10] -= 0.005;
        let sigma = mrt_stress_tensor_ce(&f, rho, u, 1.0);
        assert!((sigma[0][1] - sigma[1][0]).abs() < 1e-14, "S_xy != S_yx");
        assert!((sigma[0][2] - sigma[2][0]).abs() < 1e-14, "S_xz != S_zx");
        assert!((sigma[1][2] - sigma[2][1]).abs() < 1e-14, "S_yz != S_zy");
    }

    // ── Entropy production tests ───────────────────────────────────────────────

    // 39. entropy_production at equilibrium is ~0.
    #[test]
    fn test_entropy_production_at_equilibrium() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let sigma = entropy_production(&feq, rho, u);
        assert!(sigma.abs() < 1e-20, "entropy production at eq = {sigma}");
    }

    // 40. entropy_production: off-eq gives negative value (dissipation).
    #[test]
    fn test_entropy_production_negative_offequilibrium() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let mut f = compute_equilibrium_d3q19(rho, u);
        f[1] += 0.01;
        f[2] -= 0.01; // momentum-preserving perturbation
        let sigma = entropy_production(&f, rho, u);
        assert!(
            sigma < 0.0,
            "off-eq entropy production should be negative: {sigma}"
        );
    }

    // 41. h_function_d3q19 is finite.
    #[test]
    fn test_h_function_d3q19_finite() {
        let rho = 1.0;
        let u = [0.01, 0.0, 0.0];
        let f = compute_equilibrium_d3q19(rho, u);
        let h = h_function_d3q19(&f);
        assert!(h.is_finite(), "H function should be finite: {h}");
    }

    // 42. verify_h_theorem: H does not increase at equilibrium.
    #[test]
    fn test_verify_h_theorem_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let params = MrtRelaxationParams::default();
        let (h_before, h_after) = verify_h_theorem(&feq, rho, u, &params.s);
        assert!(
            (h_after - h_before).abs() < 1e-12,
            "H should not change at equilibrium: before={h_before}, after={h_after}"
        );
    }

    // 43. verify_h_theorem: H decreases for off-equilibrium state.
    #[test]
    fn test_verify_h_theorem_decreases_offequilibrium() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let mut f = compute_equilibrium_d3q19(rho, u);
        f[1] += 0.005;
        f[2] -= 0.005;
        let (rho_f, u_f) = compute_macroscopic(&f);
        let params = MrtRelaxationParams::default();
        let (h_before, h_after) = verify_h_theorem(&f, rho_f, u_f, &params.s);
        assert!(
            h_after <= h_before + 1e-12,
            "H should decrease or stay: before={h_before}, after={h_after}"
        );
    }

    // 44. build_relaxation_matrix: stability check.
    #[test]
    fn test_relaxation_matrix_stability() {
        let s = build_relaxation_matrix(1.2, 1.0, 1.4);
        let params = MrtRelaxationParams { s };
        assert!(
            params.is_stable(),
            "built relaxation matrix should be stable"
        );
    }

    // 45. D3Q27 equilibrium has correct zero-velocity relationship.
    #[test]
    fn test_d3q27_equilibrium_zero_u_momentum() {
        let feq = d3q27_equilibrium_full(1.0, [0.0, 0.0, 0.0]);
        let mut mx = 0.0;
        for i in 0..27 {
            mx += feq[i] * D3Q27_VELOCITIES[i][0] as f64;
        }
        assert!(mx.abs() < 1e-12, "momentum at rest should be 0: {mx}");
    }
}
