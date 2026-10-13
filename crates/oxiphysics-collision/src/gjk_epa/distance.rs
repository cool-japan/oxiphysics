// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GJK distance iteration between two convex shapes.
//!
//! The distance between convex shapes `A` and `B` is the distance from the
//! origin to their Minkowski difference `A - B`. GJK approximates that set by
//! the convex hull of at most four support points (the *simplex*):
//!
//! 1. take the point `v` of the simplex nearest to the origin together with its
//!    barycentric weights (see [`super::simplex_closest`]),
//! 2. drop the vertices that do not support `v`,
//! 3. fetch the support point `w` of `A - B` in direction `-v`,
//! 4. stop when `w` cannot improve on `v`, i.e. `|v|^2 - v.w <= tol |v|^2`,
//!    otherwise add `w` to the simplex and repeat.
//!
//! Every simplex vertex remembers the support points on `A` and `B` it was
//! built from, so the barycentric weights of `v` also give the pair of closest
//! points `(sum l_i sa_i, sum l_i sb_i)`, whose difference is exactly `v`.

use super::simplex_closest::closest_point_to_origin;
use super::{ConvexShape, add3, dot3, neg3, scale3, sub3};

/// Hard cap on the number of GJK iterations.
const MAX_ITERATIONS: usize = 128;

/// Relative tolerance of the termination test `|v|^2 - v.w <= GAP_TOLERANCE |v|^2`.
///
/// The returned distance is within `GAP_TOLERANCE * |v|` of the true distance.
/// For smooth shapes the distance error is quadratic in the position error of
/// the witness points, so the witnesses are accurate to roughly
/// `sqrt(GAP_TOLERANCE) * |v|`.
const GAP_TOLERANCE: f64 = 1e-12;

/// Shapes whose closest point of the Minkowski difference is closer to the
/// origin than this fraction of the simplex extent count as touching.
const TOUCH_TOLERANCE: f64 = 1e-12;

/// A vertex of the Minkowski-difference simplex and the support points it was
/// built from.
#[derive(Debug, Clone, Copy)]
struct SimplexVertex {
    /// Point of `A - B`: `sa - sb`.
    w: [f64; 3],
    /// Support point on shape `A`.
    sa: [f64; 3],
    /// Support point on shape `B`.
    sb: [f64; 3],
}

/// Outcome of a GJK distance query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum DistanceOutcome {
    /// The shapes overlap or touch (the origin lies in the Minkowski difference
    /// or closer to it than numerical resolution allows).
    Overlapping,
    /// The shapes are disjoint.
    Separated {
        /// Euclidean distance between the shapes.
        distance: f64,
        /// Closest point on shape `A`.
        point_a: [f64; 3],
        /// Closest point on shape `B`.
        point_b: [f64; 3],
    },
}

/// Best feature found so far.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    dist_sq: f64,
    point_a: [f64; 3],
    point_b: [f64; 3],
}

impl Candidate {
    fn into_outcome(self) -> DistanceOutcome {
        DistanceOutcome::Separated {
            distance: self.dist_sq.sqrt(),
            point_a: self.point_a,
            point_b: self.point_b,
        }
    }
}

/// Closest points of the simplex `weights` blend on shapes `A` and `B`.
fn witnesses(simplex: &[SimplexVertex], weights: &[f64; 4]) -> ([f64; 3], [f64; 3]) {
    let mut point_a = [0.0; 3];
    let mut point_b = [0.0; 3];
    for (vertex, &weight) in simplex.iter().zip(weights.iter()) {
        point_a = add3(point_a, scale3(vertex.sa, weight));
        point_b = add3(point_b, scale3(vertex.sb, weight));
    }
    (point_a, point_b)
}

/// Run the GJK distance iteration on two convex shapes.
///
/// Overlap is decided by the iteration itself (origin enclosed by the simplex,
/// or the simplex touching the origin); no separate intersection test is used.
/// If the iteration limit or numerical stagnation is reached before the
/// termination test passes, the best (smallest) feature found so far is
/// returned, which is always an upper bound on the true distance.
pub(super) fn solve(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> DistanceOutcome {
    let support = |direction: [f64; 3]| -> SimplexVertex {
        let sa = shape_a.support(direction);
        let sb = shape_b.support(neg3(direction));
        SimplexVertex {
            w: sub3(sa, sb),
            sa,
            sb,
        }
    };

    let mut simplex: Vec<SimplexVertex> = Vec::with_capacity(4);
    simplex.push(support([1.0, 0.0, 0.0]));
    let mut best: Option<Candidate> = None;

    for _ in 0..MAX_ITERATIONS {
        let count = simplex.len().min(4);
        let mut points = [[0.0_f64; 3]; 4];
        for (slot, vertex) in points.iter_mut().zip(simplex.iter()) {
            *slot = vertex.w;
        }
        let closest = closest_point_to_origin(&points[..count]);
        if closest.encloses_origin {
            return DistanceOutcome::Overlapping;
        }

        let v = closest.point;
        let v_sq = dot3(v, v);
        let extent_sq = simplex
            .iter()
            .map(|vertex| dot3(vertex.w, vertex.w))
            .fold(0.0_f64, f64::max);
        if v_sq <= TOUCH_TOLERANCE * TOUCH_TOLERANCE * extent_sq {
            return DistanceOutcome::Overlapping;
        }

        // The closest point of a growing hull can only get closer; a larger value
        // is rounding noise, so the previous (better) feature stands.
        if let Some(previous) = best
            && v_sq >= previous.dist_sq
        {
            return previous.into_outcome();
        }
        let (point_a, point_b) = witnesses(&simplex, &closest.weights);
        let candidate = Candidate {
            dist_sq: v_sq,
            point_a,
            point_b,
        };
        best = Some(candidate);

        // Keep only the vertices that support the closest point.
        let mut index = 0;
        simplex.retain(|_| {
            let keep = closest.weights.get(index).is_some_and(|w| *w > 0.0);
            index += 1;
            keep
        });

        let next = support(neg3(v));
        let gap = v_sq - dot3(v, next.w);
        let repeated = simplex.iter().any(|vertex| vertex.w == next.w);
        if gap <= GAP_TOLERANCE * v_sq || repeated {
            return candidate.into_outcome();
        }
        simplex.push(next);
    }

    match best {
        Some(candidate) => candidate.into_outcome(),
        None => DistanceOutcome::Overlapping,
    }
}
