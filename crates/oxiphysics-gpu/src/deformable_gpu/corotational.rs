//! 3×3 matrix helpers and the polar decomposition behind the corotational
//! tetrahedral element ([`CorotationalElement`](super::types::CorotationalElement)).
//!
//! All matrices are row-major `[[f64; 3]; 3]`.

use super::functions::{add3, cross3, dot3, length3, scale3, sub3};

/// Row-major 3×3 matrix.
pub(super) type Mat3 = [[f64; 3]; 3];

/// The 3×3 identity matrix.
pub(super) const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Relative off-diagonal tolerance of the one-sided Jacobi sweeps.
const JACOBI_TOL: f64 = 1.0e-15;
/// Upper bound on Jacobi sweeps (a 3×3 matrix converges in well under ten).
const JACOBI_MAX_SWEEPS: usize = 32;
/// Singular values below this fraction of the largest one count as zero.
const RANK_TOL: f64 = 1.0e-13;

/// Matrix product `a · b`.
pub(super) fn mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    out
}

/// Transpose.
pub(super) fn transpose(a: &Mat3) -> Mat3 {
    [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]],
    ]
}

/// Trace.
pub(super) fn trace(a: &Mat3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

/// Squared Frobenius norm.
pub(super) fn frobenius_sq(a: &Mat3) -> f64 {
    a.iter().flatten().map(|v| v * v).sum()
}

/// Column `index` of `a`.
fn column(a: &Mat3, index: usize) -> [f64; 3] {
    [a[0][index], a[1][index], a[2][index]]
}

/// Outer product `u vᵀ`.
fn outer(u: [f64; 3], v: [f64; 3]) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = u[i] * v[j];
        }
    }
    out
}

/// Rest-pose data of a tetrahedron with vertices `rest = [a, b, c, d]`: its
/// volume `|det Dm| / 6` and `Dm⁻ᵀ` (row-major), where
/// `Dm = [b − a | c − a | d − a]` is the rest edge matrix.
///
/// Returns `None` if the edge vectors are linearly dependent (relative
/// determinant `|det Dm| <= 1e-12 |b − a| |c − a| |d − a|`, or a `NaN`), i.e. if
/// the rest tetrahedron is degenerate and `Dm` cannot be inverted.
pub(super) fn rest_shape(rest: &[[f64; 3]; 4]) -> Option<(f64, Mat3)> {
    let ra = rest[0];
    let c0 = sub3(rest[1], ra);
    let c1 = sub3(rest[2], ra);
    let c2 = sub3(rest[3], ra);
    let det = dot3(c0, cross3(c1, c2));
    let edge_scale = length3(c0) * length3(c1) * length3(c2);
    if det.is_nan() || det.abs() <= 1.0e-12 * edge_scale {
        return None;
    }
    // Rows of Dm⁻¹ are (c1×c2, c2×c0, c0×c1)/det, so its columns -- the rows of
    // Dm⁻ᵀ returned below -- are assembled component by component.
    let inv_det = 1.0 / det;
    let r0 = scale3(cross3(c1, c2), inv_det);
    let r1 = scale3(cross3(c2, c0), inv_det);
    let r2 = scale3(cross3(c0, c1), inv_det);
    Some((
        det.abs() / 6.0,
        [
            [r0[0], r1[0], r2[0]],
            [r0[1], r1[1], r2[1]],
            [r0[2], r1[2], r2[2]],
        ],
    ))
}

/// Deformation gradient `F = Ds · Dm⁻¹` of the tetrahedron whose current vertex
/// positions are `positions = [a, b, c, d]`, with `Ds = [b − a | c − a | d − a]`
/// and `dm_inv_t = Dm⁻ᵀ` as returned by [`rest_shape`].
pub(super) fn deformation_gradient(dm_inv_t: &Mat3, positions: &[[f64; 3]; 4]) -> Mat3 {
    let pa = positions[0];
    let edges = [
        sub3(positions[1], pa),
        sub3(positions[2], pa),
        sub3(positions[3], pa),
    ];
    let mut f = [[0.0; 3]; 3];
    for (a, row) in f.iter_mut().enumerate() {
        for (b, entry) in row.iter_mut().enumerate() {
            let dm_col = dm_inv_t[b];
            *entry = edges[0][a] * dm_col[0] + edges[1][a] * dm_col[1] + edges[2][a] * dm_col[2];
        }
    }
    f
}

/// Nodal forces `[f_a, f_b, f_c, f_d]` of a tetrahedron with first
/// Piola-Kirchhoff stress `piola` (evaluated at its deformation gradient).
///
/// The forces on `b`, `c`, `d` are the columns of `−V₀ · P · Dm⁻ᵀ` and node `a`
/// takes the negative sum, so the element exerts no net force.  For a
/// hyperelastic `P = ∂ψ/∂F` they are the exact negative gradient of the strain
/// energy `V₀ ψ` with respect to the vertex positions.
pub(super) fn nodal_forces(piola: &Mat3, dm_inv_t: &Mat3, rest_volume: f64) -> [[f64; 3]; 4] {
    let mut forces = [[0.0f64; 3]; 4];
    for column in 0..3 {
        for i in 0..3 {
            let h: f64 = (0..3).map(|l| piola[i][l] * dm_inv_t[l][column]).sum();
            forces[column + 1][i] = -rest_volume * h;
        }
    }
    forces[0] = scale3(add3(add3(forces[1], forces[2]), forces[3]), -1.0);
    forces
}

/// A unit vector perpendicular to the unit vector `u`.
fn unit_perpendicular(u: [f64; 3]) -> [f64; 3] {
    let axis = (0..3)
        .min_by(|&i, &j| u[i].abs().total_cmp(&u[j].abs()))
        .unwrap_or(0);
    let mut e = [0.0; 3];
    e[axis] = 1.0;
    let p = cross3(u, e);
    scale3(p, 1.0 / length3(p))
}

/// Rotation factor `R` of the polar decomposition `F = R S` — the proper
/// rotation (`det R = +1`) closest to `F` in the Frobenius norm.
///
/// `F` is decomposed as `F = U Σ Vᵀ` with the one-sided Jacobi method, which
/// keeps small singular values accurate; `R = U diag(1, 1, det(U Vᵀ)) Vᵀ`.  An
/// inverted `F` (`det F < 0`) therefore yields a rotation too — the smallest
/// singular value simply becomes negative — and `R` stays continuous as `F`
/// passes through singularity.  `F = 0` and rank-deficient `F` return a valid
/// rotation (the identity for `F = 0`).
pub(super) fn polar_rotation(f: &Mat3) -> Mat3 {
    let mut a = *f;
    let mut v = IDENTITY;
    for _ in 0..JACOBI_MAX_SWEEPS {
        let mut rotated = false;
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
            for row in &a {
                alpha += row[p] * row[p];
                beta += row[q] * row[q];
                gamma += row[p] * row[q];
            }
            if gamma.abs() <= JACOBI_TOL * (alpha * beta).sqrt() {
                continue;
            }
            rotated = true;
            let zeta = (beta - alpha) / (2.0 * gamma);
            let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = c * t;
            for r in 0..3 {
                let (ap, aq) = (a[r][p], a[r][q]);
                a[r][p] = c * ap - s * aq;
                a[r][q] = s * ap + c * aq;
                let (vp, vq) = (v[r][p], v[r][q]);
                v[r][p] = c * vp - s * vq;
                v[r][q] = s * vp + c * vq;
            }
        }
        if !rotated {
            break;
        }
    }

    // A = F V now has orthogonal columns σᵢ uᵢ; order them by decreasing σ.
    let norms = [
        length3(column(&a, 0)),
        length3(column(&a, 1)),
        length3(column(&a, 2)),
    ];
    let mut order = [0usize, 1, 2];
    order.sort_by(|&i, &j| norms[j].total_cmp(&norms[i]));
    let [i0, i1, i2] = order;
    let sigma_max = norms[i0];
    if sigma_max.is_nan() || sigma_max <= 0.0 {
        return IDENTITY;
    }
    let u0 = scale3(column(&a, i0), 1.0 / sigma_max);
    let u1 = if norms[i1] > RANK_TOL * sigma_max {
        let raw = scale3(column(&a, i1), 1.0 / norms[i1]);
        // One Gram-Schmidt step removes the residual non-orthogonality.
        let ortho = sub3(raw, scale3(u0, dot3(raw, u0)));
        scale3(ortho, 1.0 / length3(ortho))
    } else {
        unit_perpendicular(u0)
    };
    let (v0, v1, v2) = (column(&v, i0), column(&v, i1), column(&v, i2));
    let det_v = dot3(v0, cross3(v1, v2));
    let w = scale3(cross3(u0, u1), det_v);
    let (r0, r1, r2) = (outer(u0, v0), outer(u1, v1), outer(w, v2));
    let mut rotation = [[0.0; 3]; 3];
    for (i, row) in rotation.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = r0[i][j] + r1[i][j] + r2[i][j];
        }
    }
    rotation
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rotation by `angle` about the (not necessarily unit) `axis` (Rodrigues).
    fn rotation(axis: [f64; 3], angle: f64) -> Mat3 {
        let k = scale3(axis, 1.0 / length3(axis));
        let (s, c) = angle.sin_cos();
        let kx = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]];
        let kk = mul(&kx, &kx);
        let mut r = IDENTITY;
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] += s * kx[i][j] + (1.0 - c) * kk[i][j];
            }
        }
        r
    }
    fn max_abs_diff(a: &Mat3, b: &Mat3) -> f64 {
        let mut worst = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                worst = worst.max((a[i][j] - b[i][j]).abs());
            }
        }
        worst
    }
    fn det(a: &Mat3) -> f64 {
        dot3(a[0], cross3(a[1], a[2]))
    }
    fn assert_proper_rotation(r: &Mat3) {
        let rtr = mul(&transpose(r), r);
        assert!(
            max_abs_diff(&rtr, &IDENTITY) < 1.0e-12,
            "RᵀR = {rtr:?} is not the identity"
        );
        assert!((det(r) - 1.0).abs() < 1.0e-12, "det R = {}", det(r));
    }
    #[test]
    fn polar_rotation_of_a_rotation_is_itself() {
        let q = rotation([1.0, 2.0, 3.0], 0.7);
        let r = polar_rotation(&q);
        assert!(max_abs_diff(&r, &q) < 1.0e-12, "{r:?} vs {q:?}");
    }
    #[test]
    fn polar_rotation_recovers_the_rotation_of_a_stretched_matrix() {
        let q = rotation([-0.3, 0.5, 1.0], 2.1);
        let s = [[2.0, 0.3, 0.0], [0.3, 1.5, 0.1], [0.0, 0.1, 0.8]];
        let f = mul(&q, &s);
        let r = polar_rotation(&f);
        assert!(max_abs_diff(&r, &q) < 1.0e-12, "{r:?} vs {q:?}");
        // F = R S with S = RᵀF symmetric.
        let recovered = mul(&transpose(&r), &f);
        assert!(max_abs_diff(&recovered, &transpose(&recovered)) < 1.0e-12);
    }
    #[test]
    fn polar_rotation_is_the_closest_rotation_even_for_an_inverted_matrix() {
        let q = rotation([0.2, -1.0, 0.4], 1.3);
        let f = mul(&q, &[[1.2, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, -0.5]]);
        assert!(det(&f) < 0.0);
        let r = polar_rotation(&f);
        assert_proper_rotation(&r);
        // The optimum of tr(R'ᵀ F) over rotations is σ1 + σ2 - σ3 = 1.6.
        let best = trace(&mul(&transpose(&r), &f));
        assert!((best - 1.6).abs() < 1.0e-12, "tr(RᵀF) = {best}");
        for (axis, angle) in [
            ([1.0, 0.0, 0.0], 0.3),
            ([0.0, 1.0, 1.0], 1.0),
            ([1.0, 2.0, 3.0], 2.5),
        ] {
            let other = rotation(axis, angle);
            assert!(trace(&mul(&transpose(&other), &f)) <= best + 1.0e-12);
        }
    }
    #[test]
    fn polar_rotation_of_zero_is_identity() {
        assert_eq!(polar_rotation(&[[0.0; 3]; 3]), IDENTITY);
    }
    #[test]
    fn polar_rotation_of_rank_deficient_matrices_is_a_rotation() {
        // Rank one: a bᵀ.
        let rank_one = [[2.0, 4.0, 6.0], [1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]];
        assert_proper_rotation(&polar_rotation(&rank_one));
        // Rank two: a flat matrix.
        let rank_two = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]];
        let r = polar_rotation(&rank_two);
        assert_proper_rotation(&r);
        assert!(max_abs_diff(&r, &IDENTITY) < 1.0e-12);
    }
    #[test]
    fn polar_rotation_is_continuous_through_singularity() {
        let q = rotation([1.0, 1.0, 0.0], 0.9);
        let above = polar_rotation(&mul(
            &q,
            &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0e-6]],
        ));
        let below = polar_rotation(&mul(
            &q,
            &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0e-6]],
        ));
        assert_proper_rotation(&above);
        assert_proper_rotation(&below);
        assert!(max_abs_diff(&above, &below) < 1.0e-5);
    }
    #[test]
    fn polar_rotation_handles_repeated_singular_values() {
        let q = rotation([0.0, 0.0, 1.0], -0.4);
        let f = mul(&q, &[[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]);
        let r = polar_rotation(&f);
        assert!(max_abs_diff(&r, &q) < 1.0e-12);
    }
}
