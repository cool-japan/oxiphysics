// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Flat-array vector, quaternion and 3x3 matrix helpers for the rigid-body engine.
//!
//! The engine keeps every quantity as a plain `[f64; N]` so that nothing
//! heavier than the standard library is needed at the WASM boundary.  This
//! module collects the handful of linear-algebra primitives shared by the
//! narrow phase, the ray caster and the mass-property code.

/// A 3-vector stored as a flat array `[x, y, z]`.
pub(super) type Vec3 = [f64; 3];

/// A row-major 3x3 matrix, indexed `m[row][col]`.
pub(super) type Mat3 = [[f64; 3]; 3];

/// The 3x3 identity matrix.
pub(super) const MAT3_IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// The 3x3 zero matrix.
pub(super) const MAT3_ZERO: Mat3 = [[0.0; 3]; 3];

/// Component-wise sum `a + b`.
#[inline]
pub(super) fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Component-wise difference `a - b`.
#[inline]
pub(super) fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Component-wise negation.
#[inline]
pub(super) fn neg(a: Vec3) -> Vec3 {
    [-a[0], -a[1], -a[2]]
}

/// Scalar multiple `a * s`.
#[inline]
pub(super) fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Fused `a + b * s`.
#[inline]
pub(super) fn madd(a: Vec3, b: Vec3, s: f64) -> Vec3 {
    [a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s]
}

/// Dot product.
#[inline]
pub(super) fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product `a x b`.
#[inline]
pub(super) fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared Euclidean length.
#[inline]
pub(super) fn length_sq(a: Vec3) -> f64 {
    dot(a, a)
}

/// Euclidean length.
#[inline]
pub(super) fn length(a: Vec3) -> f64 {
    length_sq(a).sqrt()
}

/// Squared distance between two points.
#[inline]
pub(super) fn dist_sq(a: Vec3, b: Vec3) -> f64 {
    length_sq(sub(a, b))
}

/// Unit vector along `a`, or `None` when `a` is too short (or non-finite) to
/// define a direction.
#[inline]
pub(super) fn normalized(a: Vec3) -> Option<Vec3> {
    let len_sq = length_sq(a);
    if len_sq.is_finite() && len_sq > 1e-300 {
        Some(scale(a, 1.0 / len_sq.sqrt()))
    } else {
        None
    }
}

/// Unit vector along the part of `v` that is perpendicular to the unit vector
/// `axis`, or `None` when `v` is parallel to `axis` (to rounding accuracy).
///
/// For a `v` that is exactly or almost exactly along `axis` the perpendicular
/// remainder `v - axis (axis . v)` is pure rounding noise whose direction has
/// nothing to do with the geometry (it is not even perpendicular to `axis`),
/// so remainders below `1e-12 |v|` are rejected, and the accepted ones are
/// re-orthogonalised against `axis` before they are normalised.
pub(super) fn unit_perpendicular_part(v: Vec3, axis: Vec3) -> Option<Vec3> {
    let v_len_sq = length_sq(v);
    let radial = sub(v, scale(axis, dot(v, axis)));
    if !(v_len_sq > 0.0 && v_len_sq.is_finite()) || length_sq(radial) <= 1e-24 * v_len_sq {
        return None;
    }
    normalized(sub(radial, scale(axis, dot(radial, axis))))
}

/// A unit vector perpendicular to the (non-zero) vector `a`.
///
/// The result is deterministic: the world axis least aligned with `a` is
/// crossed with it.
pub(super) fn any_perpendicular(a: Vec3) -> Vec3 {
    let abs = [a[0].abs(), a[1].abs(), a[2].abs()];
    let helper = if abs[0] <= abs[1] && abs[0] <= abs[2] {
        [1.0, 0.0, 0.0]
    } else if abs[1] <= abs[2] {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalized(cross(a, helper)).unwrap_or([1.0, 0.0, 0.0])
}

/// Hamilton product `a * b` of two `[x, y, z, w]` quaternions.
#[inline]
pub(super) fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
        a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
        a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2],
    ]
}

/// Normalise a quaternion; `None` when it is non-finite or (nearly) zero.
pub(super) fn quat_normalized(q: [f64; 4]) -> Option<[f64; 4]> {
    let n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if n2.is_finite() && n2 > 1e-300 {
        let inv = 1.0 / n2.sqrt();
        Some([q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv])
    } else {
        None
    }
}

/// Unit quaternion for a rotation of `angle` radians about the unit `axis`.
pub(super) fn quat_from_axis_angle(axis: Vec3, angle: f64) -> [f64; 4] {
    let half = 0.5 * angle;
    let s = half.sin();
    [axis[0] * s, axis[1] * s, axis[2] * s, half.cos()]
}

/// Rotation matrix of a unit quaternion `[x, y, z, w]`.
///
/// Column `i` of the result is the image of the `i`-th world axis, so the
/// matrix maps body-frame vectors to world-frame vectors.
pub(super) fn quat_to_mat3(q: [f64; 4]) -> Mat3 {
    let [x, y, z, w] = q;
    let (xx, yy, zz) = (x * x, y * y, z * z);
    let (xy, xz, yz) = (x * y, x * z, y * z);
    let (wx, wy, wz) = (w * x, w * y, w * z);
    [
        [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy)],
        [2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
        [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)],
    ]
}

/// Matrix-vector product `m * v`.
#[inline]
pub(super) fn mat_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// Transposed matrix-vector product `m^T * v`.
#[inline]
pub(super) fn mat_t_vec(m: &Mat3, v: Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

/// Column `i` of `m` as a vector (`i` must be `0..3`).
#[inline]
pub(super) fn mat_col(m: &Mat3, i: usize) -> Vec3 {
    [m[0][i], m[1][i], m[2][i]]
}

/// Matrix transpose.
pub(super) fn mat_transpose(m: &Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Matrix product `a * b`.
pub(super) fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let bt = mat_transpose(b);
    [
        [dot(a[0], bt[0]), dot(a[0], bt[1]), dot(a[0], bt[2])],
        [dot(a[1], bt[0]), dot(a[1], bt[1]), dot(a[1], bt[2])],
        [dot(a[2], bt[0]), dot(a[2], bt[1]), dot(a[2], bt[2])],
    ]
}

/// Element-wise matrix sum.
pub(super) fn mat_add(a: &Mat3, b: &Mat3) -> Mat3 {
    [add(a[0], b[0]), add(a[1], b[1]), add(a[2], b[2])]
}

/// Matrix scaled by `s`.
pub(super) fn mat_scale(a: &Mat3, s: f64) -> Mat3 {
    [scale(a[0], s), scale(a[1], s), scale(a[2], s)]
}

/// Diagonal matrix with the given diagonal entries.
pub(super) fn mat_diag(d: Vec3) -> Mat3 {
    [[d[0], 0.0, 0.0], [0.0, d[1], 0.0], [0.0, 0.0, d[2]]]
}

/// Skew-symmetric (cross-product) matrix `[v]x` with `[v]x * u == v x u`.
pub(super) fn mat_skew(v: Vec3) -> Mat3 {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}

/// Outer product `a * b^T`.
pub(super) fn mat_outer(a: Vec3, b: Vec3) -> Mat3 {
    [scale(b, a[0]), scale(b, a[1]), scale(b, a[2])]
}

/// Inverse of a general 3x3 matrix; `None` when it is (numerically) singular.
pub(super) fn mat_inverse(m: &Mat3) -> Option<Mat3> {
    let c0 = cross(m[1], m[2]);
    let c1 = cross(m[2], m[0]);
    let c2 = cross(m[0], m[1]);
    let det = dot(m[0], c0);
    let scale_ref = length(m[0]) * length(m[1]) * length(m[2]);
    if !det.is_finite() || !scale_ref.is_finite() || det.abs() <= 1e-14 * scale_ref {
        return None;
    }
    let inv_det = 1.0 / det;
    // Rows of the inverse are the cofactor vectors c0, c1, c2 stacked as
    // columns: inv = [c0 c1 c2] / det.
    Some([
        [c0[0] * inv_det, c1[0] * inv_det, c2[0] * inv_det],
        [c0[1] * inv_det, c1[1] * inv_det, c2[1] * inv_det],
        [c0[2] * inv_det, c1[2] * inv_det, c2[2] * inv_det],
    ])
}

/// Clamp `v` component-wise into the box `[-half, half]`.
#[inline]
pub(super) fn clamp_to_box(v: Vec3, half: Vec3) -> Vec3 {
    [
        v[0].clamp(-half[0], half[0]),
        v[1].clamp(-half[1], half[1]),
        v[2].clamp(-half[2], half[2]),
    ]
}

/// Closest points between the segments `p1..q1` and `p2..q2`.
///
/// Returns `(c1, c2)` with `c1` on the first and `c2` on the second segment
/// (Ericson, *Real-Time Collision Detection*, 5.1.9).  Degenerate
/// (zero-length) segments are handled.
pub(super) fn closest_points_segments(p1: Vec3, q1: Vec3, p2: Vec3, q2: Vec3) -> (Vec3, Vec3) {
    const EPS: f64 = 1e-14;
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);

    let (s, t);
    if a <= EPS && e <= EPS {
        s = 0.0;
        t = 0.0;
    } else if a <= EPS {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= EPS {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let mut s_candidate = if denom > EPS * a * e {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t_candidate = (b * s_candidate + f) / e;
            if t_candidate < 0.0 {
                t_candidate = 0.0;
                s_candidate = (-c / a).clamp(0.0, 1.0);
            } else if t_candidate > 1.0 {
                t_candidate = 1.0;
                s_candidate = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = s_candidate;
            t = t_candidate;
        }
    }
    (madd(p1, d1, s), madd(p2, d2, t))
}

/// Closest point to `p` on the segment `a..b`.
pub(super) fn closest_point_on_segment(p: Vec3, a: Vec3, b: Vec3) -> Vec3 {
    let ab = sub(b, a);
    let len_sq = dot(ab, ab);
    if len_sq <= 1e-300 {
        return a;
    }
    let t = (dot(sub(p, a), ab) / len_sq).clamp(0.0, 1.0);
    madd(a, ab, t)
}
