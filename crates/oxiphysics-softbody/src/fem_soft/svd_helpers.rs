// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rotation-variant (signed) 3x3 SVD for invertible FEM.
//!
//! Wraps the reflection-free `oxiphysics_core::svd3` and applies the
//! Irving et al. (2004) / Stomakhin et al. (2012) reflection convention so
//! that `det(U) = det(V) = +1`. Under this convention a negative smallest
//! singular value `sigma[2] < 0` encodes element inversion (`det(F) < 0`),
//! which lets corotational and invertible material models recover from a
//! fully inverted tetrahedron without producing NaN/Inf.

use oxiphysics_core::math::{Mat3, Vec3};
use oxiphysics_core::svd3;

/// Convert a row-major `[[f64;3];3]` into a column-major nalgebra `Mat3`.
fn array_to_mat3(a: [[f64; 3]; 3]) -> Mat3 {
    Mat3::new(
        a[0][0], a[0][1], a[0][2], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2],
    )
}

/// Convert a `Mat3` into a row-major `[[f64;3];3]`.
fn mat3_to_array(m: Mat3) -> [[f64; 3]; 3] {
    [
        [m[(0, 0)], m[(0, 1)], m[(0, 2)]],
        [m[(1, 0)], m[(1, 1)], m[(1, 2)]],
        [m[(2, 0)], m[(2, 1)], m[(2, 2)]],
    ]
}

/// Signed SVD with reflection convention.
///
/// Returns `(u, sigma, vt)` such that `F = U * diag(sigma) * Vt` with
/// `det(U) = det(Vt) = +1`. `sigma[2]` may be negative, encoding inversion.
/// All matrices are returned as row-major `[[f64;3];3]`; `vt` is V-transpose.
pub fn signed_svd3(f: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3]) {
    let m = array_to_mat3(f);
    let (mut u, sigma_v, mut v) = svd3(m);
    let mut sigma = [sigma_v.x, sigma_v.y, sigma_v.z];

    // Enforce det(U) = +1 by flipping the column of the smallest singular value.
    if u.determinant() < 0.0 {
        for r in 0..3 {
            u[(r, 2)] = -u[(r, 2)];
        }
        sigma[2] = -sigma[2];
    }
    // Enforce det(V) = +1 likewise.
    if v.determinant() < 0.0 {
        for r in 0..3 {
            v[(r, 2)] = -v[(r, 2)];
        }
        sigma[2] = -sigma[2];
    }

    let vt = v.transpose();
    (mat3_to_array(u), sigma, mat3_to_array(vt))
}

/// Reconstruct `F = U * diag(sigma) * Vt` from a signed SVD (row-major arrays).
pub fn reconstruct_from_svd(u: [[f64; 3]; 3], sigma: [f64; 3], vt: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let um = array_to_mat3(u);
    let vtm = array_to_mat3(vt);
    let s = Mat3::from_diagonal(&Vec3::new(sigma[0], sigma[1], sigma[2]));
    let recon = um * s * vtm;
    mat3_to_array(recon)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem_soft::math_helpers::det3x3;

    fn frob_diff(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
        let mut s = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let d = a[i][j] - b[i][j];
                s += d * d;
            }
        }
        s.sqrt()
    }

    // Deterministic pseudo-random matrices (LCG) so the test is reproducible.
    fn lcg_matrix(seed: &mut u64) -> [[f64; 3]; 3] {
        let mut next = || {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        [
            [next(), next(), next()],
            [next(), next(), next()],
            [next(), next(), next()],
        ]
    }

    #[test]
    fn test_signed_svd3_reconstruction() {
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut max_err = 0.0f64;
        for _ in 0..10 {
            let f = lcg_matrix(&mut seed);
            let (u, sigma, vt) = signed_svd3(f);
            let recon = reconstruct_from_svd(u, sigma, vt);
            let err = frob_diff(f, recon);
            max_err = max_err.max(err);
            assert!(err < 1e-10, "reconstruction error {err}");
            assert!((det3x3(u) - 1.0).abs() < 1e-12, "det(U) = {}", det3x3(u));
            assert!((det3x3(vt) - 1.0).abs() < 1e-12, "det(Vt) = {}", det3x3(vt));
        }
        // Surface the worst-case error for the report.
        assert!(max_err < 1e-10);
    }

    #[test]
    fn test_signed_svd3_inverted() {
        // Start from a well-conditioned matrix, flip one column to make det < 0.
        let mut f = [[2.0, 0.1, 0.0], [0.0, 1.5, 0.2], [0.1, 0.0, 1.2]];
        for r in 0..3 {
            f[r][0] = -f[r][0];
        }
        assert!(det3x3(f) < 0.0, "test matrix must be inverted");
        let (u, sigma, vt) = signed_svd3(f);
        assert!(
            sigma[2] < 0.0,
            "smallest singular value must be negative: {sigma:?}"
        );
        assert!((det3x3(u) - 1.0).abs() < 1e-12);
        assert!((det3x3(vt) - 1.0).abs() < 1e-12);
        let recon = reconstruct_from_svd(u, sigma, vt);
        assert!(frob_diff(f, recon) < 1e-10);
    }
}
