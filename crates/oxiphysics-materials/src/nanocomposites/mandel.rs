// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense 6×6 linear algebra in Mandel (Kelvin) notation.
//!
//! A symmetric second-order tensor `a` is stored as
//! `[a11, a22, a33, √2·a23, √2·a13, √2·a12]`.  A fourth-order tensor with the
//! minor symmetries becomes a 6×6 matrix whose normal-normal block equals the
//! tensor components, whose normal-shear blocks carry a factor `√2` and whose
//! shear-shear block carries a factor `2` (an isotropic stiffness has shear
//! diagonal `2G`).  In this orthonormal basis the double contraction of two
//! fourth-order tensors is the ordinary matrix product and the tensor inverse
//! is the matrix inverse, so Eshelby-type homogenization schemes reduce to
//! plain 6×6 linear algebra with no factors of two to track.

/// A 6×6 matrix in Mandel notation.
pub type Mat6 = [[f64; 6]; 6];

/// The 6×6 identity matrix (the fourth-order symmetric identity tensor).
pub fn identity6() -> Mat6 {
    let mut m = [[0.0; 6]; 6];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    m
}

/// Matrix product `a · b` (double contraction of the fourth-order tensors).
pub fn mul6(a: &Mat6, b: &Mat6) -> Mat6 {
    let mut c = [[0.0; 6]; 6];
    for (c_row, a_row) in c.iter_mut().zip(a.iter()) {
        for (a_ik, b_row) in a_row.iter().zip(b.iter()) {
            for (c_ij, b_kj) in c_row.iter_mut().zip(b_row.iter()) {
                *c_ij += a_ik * b_kj;
            }
        }
    }
    c
}

/// Entry-wise sum `a + b`.
pub fn add6(a: &Mat6, b: &Mat6) -> Mat6 {
    let mut c = *a;
    for (c_row, b_row) in c.iter_mut().zip(b.iter()) {
        for (c_ij, b_ij) in c_row.iter_mut().zip(b_row.iter()) {
            *c_ij += b_ij;
        }
    }
    c
}

/// Entry-wise difference `a - b`.
pub fn sub6(a: &Mat6, b: &Mat6) -> Mat6 {
    let mut c = *a;
    for (c_row, b_row) in c.iter_mut().zip(b.iter()) {
        for (c_ij, b_ij) in c_row.iter_mut().zip(b_row.iter()) {
            *c_ij -= b_ij;
        }
    }
    c
}

/// Scalar multiple `s · a`.
pub fn scale6(a: &Mat6, s: f64) -> Mat6 {
    let mut c = *a;
    for entry in c.iter_mut().flatten() {
        *entry *= s;
    }
    c
}

/// Symmetric part `(a + aᵀ)/2`.
pub fn symmetrize6(a: &Mat6) -> Mat6 {
    let mut c = [[0.0; 6]; 6];
    for (i, row) in c.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = 0.5 * (a[i][j] + a[j][i]);
        }
    }
    c
}

/// Inverse of a 6×6 matrix by Gauss-Jordan elimination with partial pivoting.
///
/// Returns `None` if the matrix contains non-finite entries or is singular to
/// working precision (a pivot smaller than `1e-14` of the largest entry).
pub fn inverse6(m: &Mat6) -> Option<Mat6> {
    let scale = m
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, value| acc.max(value.abs()));
    if !(scale.is_finite() && scale > 0.0) {
        return None;
    }
    let mut a = *m;
    let mut inv = identity6();
    for col in 0..6 {
        let mut pivot_row = col;
        let mut pivot_abs = a[col][col].abs();
        for (row, a_row) in a.iter().enumerate().skip(col + 1) {
            if a_row[col].abs() > pivot_abs {
                pivot_abs = a_row[col].abs();
                pivot_row = row;
            }
        }
        if pivot_abs <= 1.0e-14 * scale {
            return None;
        }
        a.swap(col, pivot_row);
        inv.swap(col, pivot_row);
        let pivot = a[col][col];
        for value in a[col].iter_mut().chain(inv[col].iter_mut()) {
            *value /= pivot;
        }
        let pivot_a = a[col];
        let pivot_inv = inv[col];
        for (row, (a_row, inv_row)) in a.iter_mut().zip(inv.iter_mut()).enumerate() {
            if row == col {
                continue;
            }
            let factor = a_row[col];
            if factor == 0.0 {
                continue;
            }
            for (a_v, p_v) in a_row.iter_mut().zip(pivot_a.iter()) {
                *a_v -= factor * p_v;
            }
            for (i_v, p_v) in inv_row.iter_mut().zip(pivot_inv.iter()) {
                *i_v -= factor * p_v;
            }
        }
    }
    if inv.iter().flatten().all(|v| v.is_finite()) {
        Some(inv)
    } else {
        None
    }
}

/// Whether the symmetric part of `m` is positive definite (Cholesky test with
/// a pivot tolerance of `1e-13` relative to the largest diagonal entry).
pub fn is_positive_definite(m: &Mat6) -> bool {
    let sym = symmetrize6(m);
    let scale = (0..6).map(|i| sym[i][i].abs()).fold(0.0_f64, f64::max);
    if !(scale.is_finite() && scale > 0.0) {
        return false;
    }
    let mut l = [[0.0_f64; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let partial: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(a, b)| a * b)
                .sum();
            let value = sym[i][j] - partial;
            if i == j {
                if value <= 1.0e-13 * scale {
                    return false;
                }
                l[i][i] = value.sqrt();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    true
}

/// Mandel stiffness matrix of an isotropic material with Young's modulus
/// `young` \[Pa\] and Poisson's ratio `nu`.
///
/// The normal block is `λ + 2G` on the diagonal and `λ` off it; the shear
/// diagonal is `2G`.
pub fn isotropic_stiffness(young: f64, nu: f64) -> Mat6 {
    let lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let shear = young / (2.0 * (1.0 + nu));
    let mut c = [[0.0; 6]; 6];
    for (i, row) in c.iter_mut().enumerate() {
        if i < 3 {
            for entry in row.iter_mut().take(3) {
                *entry = lambda;
            }
            row[i] += 2.0 * shear;
        } else {
            row[i] = 2.0 * shear;
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_abs_diff(a: &Mat6, b: &Mat6) -> f64 {
        a.iter()
            .flatten()
            .zip(b.iter().flatten())
            .fold(0.0_f64, |acc, (x, y)| acc.max((x - y).abs()))
    }

    fn sample_spd() -> Mat6 {
        // B Bᵀ + I for a fixed, well-conditioned B.
        let mut b = [[0.0; 6]; 6];
        for (i, row) in b.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = ((i * 7 + j * 3 + 1) % 11) as f64 / 11.0 - 0.4;
            }
        }
        let mut bt = [[0.0; 6]; 6];
        for (i, row) in bt.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = b[j][i];
            }
        }
        add6(&mul6(&b, &bt), &identity6())
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let m = sample_spd();
        let inv = inverse6(&m).expect("sample matrix is non-singular");
        assert!(max_abs_diff(&mul6(&m, &inv), &identity6()) < 1e-12);
        assert!(max_abs_diff(&mul6(&inv, &m), &identity6()) < 1e-12);
    }

    #[test]
    fn inverse_handles_row_swaps_and_wide_dynamic_range() {
        // Zero on the diagonal forces pivoting; the scaled block mimics Pa vs 1/Pa entries.
        let mut m = [[0.0; 6]; 6];
        m[0][1] = 2.0e9;
        m[1][0] = 3.0e9;
        m[2][2] = 5.0e11;
        m[3][4] = 1.0e9;
        m[4][3] = 1.0e9;
        m[5][5] = 7.0e10;
        m[2][0] = 1.0e9;
        let inv = inverse6(&m).expect("permuted block matrix is invertible");
        let product = mul6(&m, &inv);
        assert!(max_abs_diff(&product, &identity6()) < 1e-12);
    }

    #[test]
    fn singular_and_non_finite_matrices_have_no_inverse() {
        let mut singular = identity6();
        singular[3] = singular[2];
        assert!(inverse6(&singular).is_none());
        assert!(inverse6(&[[0.0; 6]; 6]).is_none());
        let mut nan = identity6();
        nan[1][1] = f64::NAN;
        assert!(inverse6(&nan).is_none());
    }

    #[test]
    fn arithmetic_helpers_are_consistent() {
        let a = sample_spd();
        let b = identity6();
        assert!(max_abs_diff(&sub6(&add6(&a, &b), &b), &a) < 1e-15);
        assert!(max_abs_diff(&scale6(&b, 3.0), &add6(&b, &add6(&b, &b))) < 1e-15);
        let mut skew = a;
        skew[0][1] += 0.5;
        let sym = symmetrize6(&skew);
        assert!((sym[0][1] - sym[1][0]).abs() < 1e-15);
        assert!(max_abs_diff(&symmetrize6(&a), &a) < 1e-15);
    }

    #[test]
    fn positive_definiteness_test() {
        assert!(is_positive_definite(&sample_spd()));
        assert!(is_positive_definite(&isotropic_stiffness(3.5e9, 0.35)));
        let mut negative = sample_spd();
        negative[2][2] = -1.0;
        assert!(!is_positive_definite(&negative));
        assert!(!is_positive_definite(&[[0.0; 6]; 6]));
        // Indefinite: ±1 on the leading 2×2 block.
        let mut indefinite = identity6();
        indefinite[0][1] = 2.0;
        indefinite[1][0] = 2.0;
        assert!(!is_positive_definite(&indefinite));
    }

    #[test]
    fn isotropic_stiffness_matches_lame_constants_and_inverts_to_compliance() {
        let (young, nu) = (3.5e9, 0.35);
        let c = isotropic_stiffness(young, nu);
        let lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let g = young / (2.0 * (1.0 + nu));
        assert!((c[0][0] - (lambda + 2.0 * g)).abs() < 1e-3);
        assert!((c[0][1] - lambda).abs() < 1e-3);
        assert!((c[3][3] - 2.0 * g).abs() < 1e-3);
        assert_eq!(c[0][3], 0.0);
        let s = inverse6(&c).expect("isotropic stiffness is invertible");
        assert!((s[0][0] - 1.0 / young).abs() < 1e-15 * 1e3);
        assert!((s[0][1] + nu / young).abs() < 1e-15 * 1e3);
        assert!((s[3][3] - 1.0 / (2.0 * g)).abs() < 1e-15 * 1e3);
    }
}
