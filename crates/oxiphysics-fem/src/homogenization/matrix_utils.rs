#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! 6x6 matrix utilities for Voigt notation tensors.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// 6×6 matrix utilities
// ---------------------------------------------------------------------------

/// Multiply two 6×6 matrices: result = a * b.
pub fn mat6_mul(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Add two 6×6 matrices.
pub fn mat6_add(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Subtract two 6×6 matrices: result = a - b.
pub fn mat6_sub(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    c
}

/// Scale a 6×6 matrix by scalar s.
pub fn mat6_scale(a: &[[f64; 6]; 6], s: f64) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}

/// Return the 6×6 identity matrix.
pub fn mat6_identity() -> [[f64; 6]; 6] {
    let mut m = [[0.0f64; 6]; 6];
    for i in 0..6 {
        m[i][i] = 1.0;
    }
    m
}

/// Transpose a 6×6 matrix.
pub fn mat6_transpose(a: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut t = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            t[i][j] = a[j][i];
        }
    }
    t
}

/// Max absolute value of (a\[i\]\[j\] - a\[j\]\[i\]) over all i, j.
pub fn mat6_symmetry_error(a: &[[f64; 6]; 6]) -> f64 {
    let mut max_err = 0.0f64;
    for i in 0..6 {
        for j in 0..6 {
            let err = (a[i][j] - a[j][i]).abs();
            if err > max_err {
                max_err = err;
            }
        }
    }
    max_err
}

/// Frobenius norm of a 6×6 matrix.
pub fn mat6_frobenius_norm(a: &[[f64; 6]; 6]) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..6 {
        for j in 0..6 {
            sum += a[i][j] * a[i][j];
        }
    }
    sum.sqrt()
}

/// Matrix-vector product for 6×6 matrix and 6-vector.
pub fn mat6_vec_mul(a: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut result = [0.0f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            result[i] += a[i][j] * v[j];
        }
    }
    result
}

/// Dot product of two 6-vectors.
pub fn vec6_dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..6 {
        sum += a[i] * b[i];
    }
    sum
}

/// Invert a 6×6 matrix via Gaussian elimination with partial pivoting.
/// Returns `None` if the matrix is singular.
pub fn mat6_inv(a: &[[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
    const N: usize = 6;
    let mut aug = [[0.0f64; 12]; N];
    for i in 0..N {
        for j in 0..N {
            aug[i][j] = a[i][j];
        }
        aug[i][N + i] = 1.0;
    }

    for col in 0..N {
        // Partial pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..N {
            let v = aug[row][col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-14 {
            return None;
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for j in 0..12 {
            aug[col][j] /= pivot;
        }
        for row in 0..N {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in 0..12 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    let mut inv = [[0.0f64; 6]; 6];
    for i in 0..N {
        for j in 0..N {
            inv[i][j] = aug[i][N + j];
        }
    }
    Some(inv)
}
