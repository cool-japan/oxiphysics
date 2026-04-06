// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vec3 / Quaternion conversion helpers for the Python binding layer.

#![allow(missing_docs)]

// ===========================================================================
// Vec3 / Quaternion conversion helpers (for Python binding layer)
// ===========================================================================

/// Convert a `[f64; 3]` array to `PyVec3`.
#[inline]
#[allow(dead_code)]
pub fn array_to_vec3(arr: [f64; 3]) -> crate::types::PyVec3 {
    crate::types::PyVec3::from_array(arr)
}

/// Convert `PyVec3` to a `[f64; 3]` array.
#[inline]
#[allow(dead_code)]
pub fn vec3_to_array(v: crate::types::PyVec3) -> [f64; 3] {
    v.to_array()
}

/// Normalize a quaternion `[x, y, z, w]` and return it.
///
/// If the quaternion is near-zero, returns the identity `[0,0,0,1]`.
#[inline]
#[allow(dead_code)]
pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if norm < 1e-12 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm]
    }
}

/// Multiply two quaternions q1 * q2 (Hamilton product).
///
/// Both inputs are `[x, y, z, w]`.
#[allow(dead_code)]
pub fn quat_mul(q1: [f64; 4], q2: [f64; 4]) -> [f64; 4] {
    let (ax, ay, az, aw) = (q1[0], q1[1], q1[2], q1[3]);
    let (bx, by, bz, bw) = (q2[0], q2[1], q2[2], q2[3]);
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}

/// Conjugate (inverse for unit quaternion) of `[x, y, z, w]`.
#[allow(dead_code)]
pub fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}

/// Rotate a vector `v` by unit quaternion `q`.
///
/// Uses sandwich product: v' = q * \[v,0\] * q*.
#[allow(dead_code)]
pub fn quat_rotate_vec(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let (qx, qy, qz, qw) = (q[0], q[1], q[2], q[3]);
    let (vx, vy, vz) = (v[0], v[1], v[2]);
    // t = 2 * cross(q.xyz, v)
    let tx = 2.0 * (qy * vz - qz * vy);
    let ty = 2.0 * (qz * vx - qx * vz);
    let tz = 2.0 * (qx * vy - qy * vx);
    [
        vx + qw * tx + qy * tz - qz * ty,
        vy + qw * ty + qz * tx - qx * tz,
        vz + qw * tz + qx * ty - qy * tx,
    ]
}

/// Create a quaternion from an axis-angle representation.
///
/// `axis` need not be normalised; `angle` is in radians.
#[allow(dead_code)]
pub fn quat_from_axis_angle(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if len < 1e-12 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    let half = angle * 0.5;
    let s = half.sin() / len;
    [axis[0] * s, axis[1] * s, axis[2] * s, half.cos()]
}
