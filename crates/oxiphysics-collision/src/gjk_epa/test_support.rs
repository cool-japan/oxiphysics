// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact oracles shared by the regression tests of the GJK and EPA modules.
//!
//! Everything here is brute force on purpose: a separating-axis gap over every
//! triangle normal and edge pair of two point clouds, and vertex/triangle plus
//! edge/edge enumeration for the distance between two polyhedra. The oracles
//! are slow but share no failure mode with GJK or EPA.

use super::{add3, cross3, dot3, len3, scale3, sub3};

/// Deterministic xorshift generator.
pub(crate) struct Rng(pub(crate) u64);

impl Rng {
    /// Uniform sample in `[0, 1)`.
    pub(crate) fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform sample in `[lo, hi)`.
    pub(crate) fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Vector with every component uniform in `[lo, hi)`.
    pub(crate) fn vec(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

/// Row-major rotation by `angle` radians about the (not necessarily unit) `axis`.
pub(crate) fn rotation_about(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let length = len3(axis);
    let [x, y, z] = scale3(axis, 1.0 / length);
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

/// Apply a row-major 3x3 matrix to a vector.
pub(crate) fn rotate(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [dot3(m[0], v), dot3(m[1], v), dot3(m[2], v)]
}

/// The eight corners of the box with half extents `half`, rotated by `rotation`
/// about its centre `center`.
pub(crate) fn box_corners(
    center: [f64; 3],
    half: [f64; 3],
    rotation: [[f64; 3]; 3],
) -> Vec<[f64; 3]> {
    let mut corners = Vec::with_capacity(8);
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let local = [sx * half[0], sy * half[1], sz * half[2]];
                corners.push(add3(center, rotate(rotation, local)));
            }
        }
    }
    corners
}

/// The point cloud moved by `offset`.
pub(crate) fn translated(cloud: &[[f64; 3]], offset: [f64; 3]) -> Vec<[f64; 3]> {
    cloud.iter().map(|p| add3(*p, offset)).collect()
}

/// Vertex of the point cloud furthest along `direction`.
pub(crate) fn cloud_support(cloud: &[[f64; 3]], direction: [f64; 3]) -> [f64; 3] {
    let mut best = cloud[0];
    let mut best_dot = dot3(best, direction);
    for &p in &cloud[1..] {
        let d = dot3(p, direction);
        if d > best_dot {
            best_dot = d;
            best = p;
        }
    }
    best
}

/// Smallest and largest projection of the cloud on the unit axis `u`.
pub(crate) fn interval(cloud: &[[f64; 3]], u: [f64; 3]) -> (f64, f64) {
    cloud
        .iter()
        .map(|p| dot3(*p, u))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
            (lo.min(x), hi.max(x))
        })
}

/// Distance from point `p` to segment `a`-`b`.
fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let t = (dot3(sub3(p, a), ab) / dot3(ab, ab)).clamp(0.0, 1.0);
    len3(sub3(p, add3(a, scale3(ab, t))))
}

/// Closest-point distance between segments `p1`-`q1` and `p2`-`q2`.
fn segment_segment_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    let d1 = sub3(q1, p1);
    let d2 = sub3(q2, p2);
    let r = sub3(p1, p2);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    let c = dot3(d1, r);
    let b = dot3(d1, d2);
    let denom = a * e - b * b;
    let mut s = if denom > 1e-12 {
        ((b * f - c * e) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    len3(sub3(add3(p1, scale3(d1, s)), add3(p2, scale3(d2, t))))
}

/// Distance from `p` to the triangle `abc`: the plane distance when the foot of
/// the perpendicular lies inside the triangle, otherwise the distance to the
/// nearest edge.
fn point_triangle_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let edges = point_segment_distance(p, a, b)
        .min(point_segment_distance(p, b, c))
        .min(point_segment_distance(p, c, a));
    let normal = cross3(sub3(b, a), sub3(c, a));
    let normal_sq = dot3(normal, normal);
    if normal_sq < 1e-24 {
        return edges;
    }
    let offset = dot3(sub3(p, a), normal) / normal_sq;
    let foot = sub3(p, scale3(normal, offset));
    let inside = [(a, b), (b, c), (c, a)]
        .iter()
        .all(|&(from, to)| dot3(cross3(sub3(to, from), sub3(foot, from)), normal) >= 0.0);
    if inside {
        offset.abs() * normal_sq.sqrt()
    } else {
        edges
    }
}

/// Exact distance between the convex hulls of two disjoint point sets.
///
/// The closest pair of two disjoint polyhedra is a vertex and a face or two
/// edges. Every face is covered by the triangles of its vertices and every edge
/// is a vertex pair, so the minimum over all vertex/triangle and edge/edge
/// combinations (chords through the interior can only be farther away) is exact.
pub(crate) fn polyhedron_distance(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut best = f64::INFINITY;
    for (own, other) in [(a, b), (b, a)] {
        for &p in own {
            for i in 0..other.len() {
                for j in (i + 1)..other.len() {
                    for k in (j + 1)..other.len() {
                        best = best.min(point_triangle_distance(p, other[i], other[j], other[k]));
                    }
                }
            }
        }
    }
    for i in 0..a.len() {
        for j in (i + 1)..a.len() {
            for k in 0..b.len() {
                for l in (k + 1)..b.len() {
                    best = best.min(segment_segment_distance(a[i], a[j], b[k], b[l]));
                }
            }
        }
    }
    best
}

/// Largest separating-axis gap of two convex hulls.
pub(crate) struct Separation {
    /// `max(lo_b - hi_a, lo_a - hi_b)` along the best axis: positive exactly
    /// when the hulls are disjoint, otherwise minus the penetration depth.
    pub(crate) gap: f64,
    /// Unit direction, along the best axis, that points from hull B toward hull
    /// A. For penetrating hulls it is the direction in which A must move by
    /// `-gap` to separate (the contact normal "from B to A").
    pub(crate) from_b_to_a: [f64; 3],
}

/// Candidate separating axes of two convex hulls given as point clouds: the
/// normals of all triangles spanned by three points of either cloud and the
/// cross products of all point-pair directions, a superset of the face normals
/// and edge pairs of the hulls.
pub(crate) fn polyhedron_axes(a: &[[f64; 3]], b: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let mut axes = Vec::new();
    for cloud in [a, b] {
        for i in 0..cloud.len() {
            for j in (i + 1)..cloud.len() {
                for k in (j + 1)..cloud.len() {
                    axes.push(cross3(sub3(cloud[j], cloud[i]), sub3(cloud[k], cloud[i])));
                }
            }
        }
    }
    for i in 0..a.len() {
        for j in (i + 1)..a.len() {
            for k in 0..b.len() {
                for l in (k + 1)..b.len() {
                    axes.push(cross3(sub3(a[j], a[i]), sub3(b[l], b[k])));
                }
            }
        }
    }
    axes
}

/// The 15 separating axes of two rotated boxes: the three face normals of each
/// and the nine cross products of their edge directions.
pub(crate) fn box_axes(rotation_a: [[f64; 3]; 3], rotation_b: [[f64; 3]; 3]) -> Vec<[f64; 3]> {
    let column = |m: [[f64; 3]; 3], k: usize| [m[0][k], m[1][k], m[2][k]];
    let mut axes = Vec::with_capacity(15);
    for k in 0..3 {
        axes.push(column(rotation_a, k));
        axes.push(column(rotation_b, k));
    }
    for i in 0..3 {
        for j in 0..3 {
            axes.push(cross3(column(rotation_a, i), column(rotation_b, j)));
        }
    }
    axes
}

/// Largest separating-axis gap of two point clouds over the candidate `axes`
/// (axes of near-zero length are skipped).
pub(crate) fn axis_separation(a: &[[f64; 3]], b: &[[f64; 3]], axes: &[[f64; 3]]) -> Separation {
    let mut best = Separation {
        gap: f64::NEG_INFINITY,
        from_b_to_a: [0.0; 3],
    };
    for &axis in axes {
        let length = len3(axis);
        if length < 1e-9 {
            continue;
        }
        let u = scale3(axis, 1.0 / length);
        let (a_lo, a_hi) = interval(a, u);
        let (b_lo, b_hi) = interval(b, u);
        let a_below_b = b_lo - a_hi;
        let a_above_b = a_lo - b_hi;
        let gap = a_below_b.max(a_above_b);
        if gap > best.gap {
            best.gap = gap;
            best.from_b_to_a = if a_below_b >= a_above_b {
                scale3(u, -1.0)
            } else {
                u
            };
        }
    }
    best
}

/// Separating-axis gap of two convex hulls given as point clouds, over every
/// axis of [`polyhedron_axes`].
pub(crate) fn polyhedron_separation(a: &[[f64; 3]], b: &[[f64; 3]]) -> Separation {
    axis_separation(a, b, &polyhedron_axes(a, b))
}
