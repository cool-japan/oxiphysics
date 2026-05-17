// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Primitive signed distance functions for common shapes.

use super::helpers::{clamp, dot2, dot3, length2, norm3, scale2, scale3, sub2, sub3};

/// Signed distance from point `p` to a sphere centred at the origin with radius `r`.
///
/// Negative inside, positive outside.
pub fn sdf_sphere(p: [f64; 3], r: f64) -> f64 {
    norm3(p) - r
}

/// Signed distance from point `p` to an axis-aligned box centred at the origin
/// with half-extents `b`.
pub fn sdf_box(p: [f64; 3], b: [f64; 3]) -> f64 {
    let q = [p[0].abs() - b[0], p[1].abs() - b[1], p[2].abs() - b[2]];
    let pos_part = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
    norm3(pos_part) + q[0].max(q[1]).max(q[2]).min(0.0)
}

/// Signed distance from point `p` to a rounded box with half-extents `b` and
/// corner radius `r`.
pub fn sdf_rounded_box(p: [f64; 3], b: [f64; 3], r: f64) -> f64 {
    sdf_box(p, b) - r
}

/// Signed distance from `p` to a capsule defined by segment endpoints `a`, `b`
/// with radius `r`.
pub fn sdf_capsule(p: [f64; 3], a: [f64; 3], b: [f64; 3], r: f64) -> f64 {
    let pa = sub3(p, a);
    let ba = sub3(b, a);
    let h = clamp(dot3(pa, ba) / dot3(ba, ba).max(1e-30), 0.0, 1.0);
    norm3(sub3(pa, scale3(ba, h))) - r
}

/// Signed distance from `p` to an infinite cylinder along the y-axis with radius `r`.
pub fn sdf_cylinder_infinite(p: [f64; 3], r: f64) -> f64 {
    let xz = [p[0], p[2]];
    length2(xz) - r
}

/// Signed distance from `p` to a finite cylinder centred at the origin,
/// aligned along y-axis, with radius `r` and half-height `h`.
pub fn sdf_cylinder(p: [f64; 3], r: f64, h: f64) -> f64 {
    let d_xy = length2([p[0], p[2]]) - r;
    let d_z = p[1].abs() - h;
    let outer = [d_xy.max(0.0), d_z.max(0.0)];
    length2(outer) + d_xy.min(0.0).max(d_z.min(0.0))
}

/// Signed distance from `p` to a torus in the xz-plane with major radius `r1`
/// and minor radius `r2`.
pub fn sdf_torus(p: [f64; 3], r1: f64, r2: f64) -> f64 {
    let xz = [p[0], p[2]];
    let q = [length2(xz) - r1, p[1]];
    length2(q) - r2
}

/// Signed distance from `p` to an infinite plane defined by normal `n` (must be
/// unit length) and scalar offset `d` (plane equation: n·x = d).
pub fn sdf_plane(p: [f64; 3], n: [f64; 3], d: f64) -> f64 {
    dot3(p, n) - d
}

/// Signed distance from `p` to a cone with apex at the origin pointing along +y,
/// half-angle `angle` (radians), height `h`.
pub fn sdf_cone(p: [f64; 3], angle: f64, h: f64) -> f64 {
    let q = length2([p[0], p[2]]);
    let k = [angle.sin(), angle.cos()];
    let w = [q, -p[1]]; // q, -y
    let a = sub2(
        w,
        scale2(k, clamp(dot2(w, k) / dot2(k, k).max(1e-30), 0.0, h)),
    );
    let b = sub2(w, [k[0] * h, -k[1] * h]);
    let s = if w[1] * k[0] - w[0] * k[1] > 0.0 {
        -1.0
    } else {
        1.0
    };
    let la = length2(a);
    let lb = length2(b);
    s * la.min(lb)
}

/// Signed distance from `p` to a line segment defined by endpoints `a` and `b`.
pub fn sdf_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let pa = sub3(p, a);
    let ba = sub3(b, a);
    let h = clamp(dot3(pa, ba) / dot3(ba, ba).max(1e-30), 0.0, 1.0);
    norm3(sub3(pa, scale3(ba, h)))
}

/// Signed distance from `p` to an ellipsoid with semi-axes `r`.
///
/// Uses an approximate formula (Quilez 2018).
pub fn sdf_ellipsoid(p: [f64; 3], r: [f64; 3]) -> f64 {
    let k0 = norm3([p[0] / r[0], p[1] / r[1], p[2] / r[2]]);
    let k1 = norm3([
        p[0] / (r[0] * r[0]),
        p[1] / (r[1] * r[1]),
        p[2] / (r[2] * r[2]),
    ]);
    if k1 < 1e-30 {
        return -r[0].min(r[1]).min(r[2]);
    }
    k0 * (k0 - 1.0) / k1
}
