// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact distance queries against axis-aligned boxes.
//!
//! These are the closed-form building blocks of the capsule and box casts:
//! the distance between a segment and a box, the distance between two boxes,
//! and the contact normals to fall back on when the shapes overlap and the
//! separation vector vanishes.

use super::advancement::separation_direction;
use super::{add3, len3, scale3, sub3};

/// The point of the box `[lo, hi]` closest to `p`.
pub(super) fn clamp_to_box(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> [f64; 3] {
    [
        p[0].clamp(lo[0], hi[0]),
        p[1].clamp(lo[1], hi[1]),
        p[2].clamp(lo[2], hi[2]),
    ]
}

/// Outward unit normal of the face of the box `[lo, hi]` that `p` is closest
/// to (or has crossed furthest, if `p` is outside the box).
///
/// This is the shortest way out of the box for a point inside it, which makes it
/// the natural contact normal when the closest point on the box coincides with
/// the query point.
pub(super) fn box_exit_normal(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> [f64; 3] {
    let mut best_axis = 0usize;
    let mut best_sign = 1.0;
    let mut best_distance = f64::INFINITY;
    for axis in 0..3 {
        let to_lower = p[axis] - lo[axis];
        let to_upper = hi[axis] - p[axis];
        if to_lower < best_distance {
            best_distance = to_lower;
            best_axis = axis;
            best_sign = -1.0;
        }
        if to_upper < best_distance {
            best_distance = to_upper;
            best_axis = axis;
            best_sign = 1.0;
        }
    }
    let mut normal = [0.0; 3];
    normal[best_axis] = best_sign;
    normal
}

/// The eight corners of the box with half extents `half`, relative to its centre.
pub(super) fn box_corners(half: [f64; 3]) -> [[f64; 3]; 8] {
    std::array::from_fn(|i| {
        [
            if i & 1 == 0 { -half[0] } else { half[0] },
            if i & 2 == 0 { -half[1] } else { half[1] },
            if i & 4 == 0 { -half[2] } else { half[2] },
        ]
    })
}

/// Squared distance from `p` to the box `[lo, hi]`.
fn point_box_distance_sq(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let diff = sub3(p, clamp_to_box(p, lo, hi));
    diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]
}

/// Closest pair between a segment and a box.
pub(super) struct SegmentBoxDistance {
    /// Distance between the segment and the box (zero if they intersect).
    pub distance: f64,
    /// Closest point on the segment (a point inside the box if they intersect).
    pub segment_point: [f64; 3],
    /// Closest point of the box.
    pub box_point: [f64; 3],
}

/// Exact distance between the segment `p0..p1` and the box `[lo, hi]`.
///
/// With `p(s) = p0 + s (p1 - p0)` the squared distance
/// `f(s) = |p(s) - clamp(p(s), lo, hi)|^2` is a sum of squared per-axis
/// excesses over the slab boundaries, hence convex and piecewise quadratic with
/// breakpoints at the parameters where the segment crosses a slab boundary.
/// Between two consecutive breakpoints every axis is either inside its slab
/// (no contribution) or on a fixed side of it, so `f` is a single quadratic
/// whose constrained minimiser has a closed form. The global minimum is the best
/// of those piece minima.
pub(super) fn segment_box_distance(
    p0: [f64; 3],
    p1: [f64; 3],
    lo: [f64; 3],
    hi: [f64; 3],
) -> SegmentBoxDistance {
    let d = sub3(p1, p0);
    let mut breaks = [0.0_f64; 8];
    breaks[1] = 1.0;
    let mut count = 2;
    for axis in 0..3 {
        if d[axis] != 0.0 {
            for bound in [lo[axis], hi[axis]] {
                let s = (bound - p0[axis]) / d[axis];
                if s > 0.0 && s < 1.0 {
                    breaks[count] = s;
                    count += 1;
                }
            }
        }
    }
    let breaks = &mut breaks[..count];
    breaks.sort_unstable_by(f64::total_cmp);

    let point_at = |s: f64| add3(p0, scale3(d, s));
    let mut best_s = 0.0;
    let mut best_sq = f64::INFINITY;
    for piece in breaks.windows(2) {
        let (s_lo, s_hi) = (piece[0], piece[1]);
        if s_hi <= s_lo {
            continue;
        }
        let s_mid = 0.5 * (s_lo + s_hi);
        let mid = point_at(s_mid);
        // On this piece axis `k` contributes `(p0[k] + s d[k] - bound[k])^2` if
        // the segment is outside its slab, where `bound` is the slab side it is
        // on; setting the derivative of the sum to zero gives the minimiser.
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for axis in 0..3 {
            let bound = if mid[axis] < lo[axis] {
                lo[axis]
            } else if mid[axis] > hi[axis] {
                hi[axis]
            } else {
                continue;
            };
            numerator += (p0[axis] - bound) * d[axis];
            denominator += d[axis] * d[axis];
        }
        let s = if denominator > 0.0 {
            (-numerator / denominator).clamp(s_lo, s_hi)
        } else {
            s_mid
        };
        let sq = point_box_distance_sq(point_at(s), lo, hi);
        if sq < best_sq {
            best_sq = sq;
            best_s = s;
        }
    }
    let segment_point = point_at(best_s);
    let box_point = clamp_to_box(segment_point, lo, hi);
    SegmentBoxDistance {
        distance: len3(sub3(segment_point, box_point)),
        segment_point,
        box_point,
    }
}

/// Closest pair between two axis-aligned boxes.
pub(super) struct AabbSeparation {
    /// Distance between the boxes (zero if they overlap or touch).
    pub distance: f64,
    /// Unit vector from box B towards box A along the closest points, or `None`
    /// if the boxes overlap or touch.
    pub direction: Option<[f64; 3]>,
    /// Closest point of box A.
    pub point_a: [f64; 3],
    /// Closest point of box B.
    pub point_b: [f64; 3],
}

/// Exact distance and closest points of the boxes `[lo_a, hi_a]` and `[lo_b, hi_b]`.
///
/// On an axis where the boxes are separated the closest points are the facing
/// faces; on an axis where their extents overlap both closest points sit in the
/// middle of the overlap.
pub(super) fn aabb_separation(
    lo_a: [f64; 3],
    hi_a: [f64; 3],
    lo_b: [f64; 3],
    hi_b: [f64; 3],
) -> AabbSeparation {
    let axes: [(f64, f64); 3] = std::array::from_fn(|k| {
        if lo_a[k] > hi_b[k] {
            (lo_a[k], hi_b[k])
        } else if lo_b[k] > hi_a[k] {
            (hi_a[k], lo_b[k])
        } else {
            let middle = 0.5 * (lo_a[k].max(lo_b[k]) + hi_a[k].min(hi_b[k]));
            (middle, middle)
        }
    });
    let point_a = [axes[0].0, axes[1].0, axes[2].0];
    let point_b = [axes[0].1, axes[1].1, axes[2].1];
    AabbSeparation {
        distance: len3(sub3(point_a, point_b)),
        direction: separation_direction(point_a, point_b),
        point_a,
        point_b,
    }
}

/// Contact normal (pointing from box B to box A) of two overlapping boxes: the
/// axis of least penetration, signed by the offset of the centres.
pub(super) fn aabb_overlap_normal(
    center_a: [f64; 3],
    half_a: [f64; 3],
    center_b: [f64; 3],
    half_b: [f64; 3],
) -> [f64; 3] {
    let mut best_axis = 0usize;
    let mut best_sign = 1.0;
    let mut best_penetration = f64::INFINITY;
    for axis in 0..3 {
        let offset = center_a[axis] - center_b[axis];
        let penetration = half_a[axis] + half_b[axis] - offset.abs();
        if penetration < best_penetration {
            best_penetration = penetration;
            best_axis = axis;
            best_sign = if offset < 0.0 { -1.0 } else { 1.0 };
        }
    }
    let mut normal = [0.0; 3];
    normal[best_axis] = best_sign;
    normal
}
