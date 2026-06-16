// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SVD-based per-triangle strain limiting for cloth simulation.
//!
//! Uses the 2×2 deformation gradient (computed in the local tangent frame of
//! each rest triangle) to clamp singular values into `[min_stretch, max_stretch]`
//! and redistribute position corrections in a momentum-conserving way.

/// Per-triangle strain-limiting parameters.
#[derive(Debug, Clone)]
pub struct StrainLimitConfig {
    /// Maximum allowed singular value (stretch). Default `1.05`.
    pub max_stretch: f64,
    /// Minimum allowed singular value (compression). Default `0.95`.
    pub min_stretch: f64,
    /// Number of Gauss–Seidel sweeps per call. Default `1`.
    pub iterations: u32,
}

impl Default for StrainLimitConfig {
    fn default() -> Self {
        Self {
            max_stretch: 1.05,
            min_stretch: 0.95,
            iterations: 1,
        }
    }
}

// ─────────────────────────── helpers ─────────────────────────────────────────

/// Normalise a 3-vector; returns the zero vector if the norm is too small.
#[inline]
fn safe_normalize(v: [f64; 3]) -> [f64; 3] {
    let n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if n2 < 1e-30 {
        return [0.0; 3];
    }
    let inv = 1.0 / n2.sqrt();
    [v[0] * inv, v[1] * inv, v[2] * inv]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// 2×2 matrix-vector multiply: M is row-major [[row0], [row1]].
#[inline]
fn mat2_mul_vec2(m: [[f64; 2]; 2], v: [f64; 2]) -> [f64; 2] {
    [
        m[0][0] * v[0] + m[0][1] * v[1],
        m[1][0] * v[0] + m[1][1] * v[1],
    ]
}

/// 2×2 matrix multiply: A * B (both row-major).
#[inline]
fn mat2_mul(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// Transpose of a 2×2 row-major matrix.
#[inline]
fn mat2_transpose(m: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [[m[0][0], m[1][0]], [m[0][1], m[1][1]]]
}

/// Inverse of a 2×2 row-major matrix; returns `None` if singular.
#[inline]
fn mat2_inverse(m: [[f64; 2]; 2]) -> Option<[[f64; 2]; 2]> {
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    if det.abs() < 1e-30 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [m[1][1] * inv_det, -m[0][1] * inv_det],
        [-m[1][0] * inv_det, m[0][0] * inv_det],
    ])
}

// ─────────────────────────── 2×2 SVD ─────────────────────────────────────────

/// Analytic 2×2 SVD.
///
/// Returns `(U, sigma, Vt)` (both row-major 2×2) such that
/// `F = U * diag(sigma) * Vt`, with `sigma[0] >= sigma[1] >= 0`.
fn svd2x2(f: [[f64; 2]; 2]) -> ([[f64; 2]; 2], [f64; 2], [[f64; 2]; 2]) {
    // Compute A = F^T * F (symmetric 2x2).
    let ft = mat2_transpose(f);
    let a = mat2_mul(ft, f);

    // Eigendecomposition of symmetric A = [[a,b],[b,d]].
    let aa = a[0][0];
    let b = a[0][1]; // == a[1][0]
    let d = a[1][1];

    let trace_half = (aa + d) * 0.5;
    let disc = ((aa - d) * 0.5) * ((aa - d) * 0.5) + b * b;
    let disc_sqrt = disc.sqrt();

    let lam1 = (trace_half + disc_sqrt).max(0.0);
    let lam2 = (trace_half - disc_sqrt).max(0.0);

    let sigma1 = lam1.sqrt();
    let sigma2 = lam2.sqrt();

    // Eigenvectors of A (= right-singular vectors V).
    // For eigenvalue lam1:
    let vt = if b.abs() > 1e-14 {
        // First eigenvector: [b, lam1 - aa] normalized.
        let v1 = safe_normalize_2d([b, lam1 - aa]);
        let v2 = [-v1[1], v1[0]]; // perpendicular
        // V columns are v1, v2; Vt rows are v1^T, v2^T.
        [[v1[0], v1[1]], [v2[0], v2[1]]]
    } else {
        // Already diagonal.
        if aa >= d {
            [[1.0_f64, 0.0], [0.0, 1.0_f64]]
        } else {
            [[0.0_f64, 1.0], [1.0_f64, 0.0]]
        }
    };

    // U = F * V * Sigma^{-1}.  Column k of U = F * v_k / sigma_k.
    // v_k is row k of Vt.
    let u = compute_u(f, vt, sigma1, sigma2);

    ([u[0], u[1]], [sigma1, sigma2], vt)
}

#[inline]
fn safe_normalize_2d(v: [f64; 2]) -> [f64; 2] {
    let n2 = v[0] * v[0] + v[1] * v[1];
    if n2 < 1e-30 {
        return [1.0, 0.0];
    }
    let inv = 1.0 / n2.sqrt();
    [v[0] * inv, v[1] * inv]
}

/// Compute U = F V Sigma^{-1}, returned as row-major [[row0],[row1]].
#[inline]
fn compute_u(f: [[f64; 2]; 2], vt: [[f64; 2]; 2], s1: f64, s2: f64) -> [[f64; 2]; 2] {
    // Column 0 of V = row 0 of Vt = vt[0].
    let v0 = [vt[0][0], vt[0][1]];
    // Column 1 of V = row 1 of Vt = vt[1].
    let v1 = [vt[1][0], vt[1][1]];

    // u_col0 = F * v0 / s1
    let fv0 = mat2_mul_vec2(f, v0);
    let u_col0 = if s1 > 1e-30 {
        [fv0[0] / s1, fv0[1] / s1]
    } else {
        safe_normalize_2d(fv0)
    };

    // u_col1 = F * v1 / s2
    let fv1 = mat2_mul_vec2(f, v1);
    let u_col1 = if s2 > 1e-30 {
        [fv1[0] / s2, fv1[1] / s2]
    } else {
        // Perpendicular to u_col0.
        [-u_col0[1], u_col0[0]]
    };

    // Return as row-major: U[row][col].
    [[u_col0[0], u_col1[0]], [u_col0[1], u_col1[1]]]
}

// ─────────────────────────── main API ────────────────────────────────────────

/// Apply SVD-based per-triangle strain limiting.
///
/// For each triangle the 2-D deformation gradient is computed in the local
/// tangent frame of the rest triangle, decomposed with an analytic 2×2 SVD,
/// and any singular values outside `[min_stretch, max_stretch]` are clamped.
/// The resulting position corrections are distributed among the three vertices
/// in proportion to their inverse masses (momentum-conserving).
///
/// # Arguments
///
/// * `positions`      – current 3-D positions, modified in-place.
/// * `rest_positions` – reference (rest) 3-D positions, not modified.
/// * `triangles`      – index triples `[i0, i1, i2]` into the position arrays.
/// * `inv_masses`     – per-particle inverse mass (`0` = pinned / static).
/// * `config`         – algorithm parameters.
pub fn apply_strain_limiting(
    positions: &mut [[f64; 3]],
    rest_positions: &[[f64; 3]],
    triangles: &[[usize; 3]],
    inv_masses: &[f64],
    config: &StrainLimitConfig,
) {
    for _ in 0..config.iterations {
        for tri in triangles {
            let [i0, i1, i2] = *tri;

            // Rest edge vectors.
            let q0 = rest_positions[i0];
            let q1 = rest_positions[i1];
            let q2 = rest_positions[i2];
            let e0_rest = sub3(q1, q0);
            let e1_rest = sub3(q2, q0);

            // Degenerate rest triangle guard.
            let area_vec = cross3(e0_rest, e1_rest);
            let area2 = dot3(area_vec, area_vec);
            if area2 < 1e-24 {
                // rest triangle is degenerate — skip
                continue;
            }

            // Build local orthonormal 2-D frame in the rest triangle's plane.
            let t1 = safe_normalize(e0_rest);
            let e1_proj_t1 = dot3(e1_rest, t1);
            let e1_orth = sub3(e1_rest, scale3(t1, e1_proj_t1));
            let t2 = safe_normalize(e1_orth);

            // Dm (rest shape matrix, 2×2, column-major stored as row-major here):
            //   col0 = [e0_rest · t1,  e0_rest · t2]  = [|e0_rest|, 0]
            //   col1 = [e1_rest · t1,  e1_rest · t2]
            let dm: [[f64; 2]; 2] = [
                [dot3(e0_rest, t1), dot3(e1_rest, t1)],
                [dot3(e0_rest, t2), dot3(e1_rest, t2)],
            ];

            let dm_inv = match mat2_inverse(dm) {
                Some(inv) => inv,
                None => continue,
            };

            // Current edge vectors.
            let p0 = positions[i0];
            let p1 = positions[i1];
            let p2 = positions[i2];
            let f0 = sub3(p1, p0);
            let f1 = sub3(p2, p0);

            // Ds (current shape matrix in local frame).
            let ds: [[f64; 2]; 2] = [[dot3(f0, t1), dot3(f1, t1)], [dot3(f0, t2), dot3(f1, t2)]];

            // Deformation gradient F_2d = Ds * Dm_inv.
            let f_2d = mat2_mul(ds, dm_inv);

            // 2×2 SVD.
            let (u, sigma, vt) = svd2x2(f_2d);

            // Clamp singular values.
            let sigma0_clamped = sigma[0].clamp(config.min_stretch, config.max_stretch);
            let sigma1_clamped = sigma[1].clamp(config.min_stretch, config.max_stretch);

            if (sigma0_clamped - sigma[0]).abs() < 1e-15
                && (sigma1_clamped - sigma[1]).abs() < 1e-15
            {
                // No clamping needed.
                continue;
            }

            // Reconstruct F' = U * diag(sigma') * Vt.
            let sigma_clamped_mat: [[f64; 2]; 2] = [[sigma0_clamped, 0.0], [0.0, sigma1_clamped]];
            let f_prime = mat2_mul(mat2_mul(u, sigma_clamped_mat), vt);

            // Delta_F = F' - F_2d (in local frame).
            let delta_f: [[f64; 2]; 2] = [
                [f_prime[0][0] - f_2d[0][0], f_prime[0][1] - f_2d[0][1]],
                [f_prime[1][0] - f_2d[1][0], f_prime[1][1] - f_2d[1][1]],
            ];

            // Delta_Ds = Delta_F * Dm   (convert back from F-space to edge space).
            let delta_ds = mat2_mul(delta_f, dm);

            // Convert corrections from 2-D local frame back to 3-D.
            // delta_ds[row][col]: row=local coord (t1/t2), col=edge index (0=e0, 1=e1).
            let delta_e0 = add3(scale3(t1, delta_ds[0][0]), scale3(t2, delta_ds[1][0]));
            let delta_e1 = add3(scale3(t1, delta_ds[0][1]), scale3(t2, delta_ds[1][1]));

            let w0 = inv_masses[i0];
            let w1 = inv_masses[i1];
            let w2 = inv_masses[i2];

            // Distribute delta_e0 between i0 and i1 (momentum-conserving).
            let (c0a, c1a) = {
                let denom = w0 + w1;
                if denom < 1e-30 {
                    ([0.0_f64; 3], [0.0_f64; 3])
                } else {
                    (scale3(delta_e0, -w0 / denom), scale3(delta_e0, w1 / denom))
                }
            };

            // Distribute delta_e1 between i0 and i2 (momentum-conserving).
            let (c0b, c2b) = {
                let denom = w0 + w2;
                if denom < 1e-30 {
                    ([0.0_f64; 3], [0.0_f64; 3])
                } else {
                    (scale3(delta_e1, -w0 / denom), scale3(delta_e1, w2 / denom))
                }
            };

            // Apply.
            positions[i0] = add3(positions[i0], add3(c0a, c0b));
            positions[i1] = add3(positions[i1], c1a);
            positions[i2] = add3(positions[i2], c2b);
        }
    }
}

/// Compute cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svd2x2_identity() {
        let f = [[1.0_f64, 0.0], [0.0, 1.0_f64]];
        let (u, sigma, vt) = svd2x2(f);
        assert!((sigma[0] - 1.0).abs() < 1e-12);
        assert!((sigma[1] - 1.0).abs() < 1e-12);
        // Reconstruct.
        let s_mat = [[sigma[0], 0.0], [0.0, sigma[1]]];
        let recon = mat2_mul(mat2_mul(u, s_mat), vt);
        for i in 0..2 {
            for j in 0..2 {
                assert!((recon[i][j] - f[i][j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn svd2x2_diagonal() {
        let f = [[3.0_f64, 0.0], [0.0, 2.0_f64]];
        let (u, sigma, vt) = svd2x2(f);
        assert!((sigma[0] - 3.0).abs() < 1e-10, "sigma[0]={}", sigma[0]);
        assert!((sigma[1] - 2.0).abs() < 1e-10, "sigma[1]={}", sigma[1]);
        let s_mat = [[sigma[0], 0.0], [0.0, sigma[1]]];
        let recon = mat2_mul(mat2_mul(u, s_mat), vt);
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    (recon[i][j] - f[i][j]).abs() < 1e-10,
                    "recon[{i}][{j}]={}",
                    recon[i][j]
                );
            }
        }
    }
}
