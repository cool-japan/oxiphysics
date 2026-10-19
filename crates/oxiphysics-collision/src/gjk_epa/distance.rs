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
//!
//! [`solve`] runs the iteration on two [`ConvexShape`]s; [`solve_with`] takes
//! any support function of a Minkowski difference (for example one that moves
//! the shapes with time) and also returns the final simplex, so that the
//! extended GJK solver, the time-of-impact search and the EPA seeding all share
//! this single, tested implementation.

use super::simplex_closest::closest_point_to_origin;
use super::{ConvexShape, add3, dot3, neg3, scale3, sub3};

/// Default hard cap on the number of GJK iterations.
pub(crate) const MAX_ITERATIONS: usize = 128;

/// Default relative tolerance of the termination test
/// `|v|^2 - v.w <= GAP_TOLERANCE |v|^2`.
///
/// The returned distance is within `GAP_TOLERANCE * |v|` of the true distance.
/// For smooth shapes the distance error is quadratic in the position error of
/// the witness points, so the witnesses are accurate to roughly
/// `sqrt(GAP_TOLERANCE) * |v|`.
pub(crate) const GAP_TOLERANCE: f64 = 1e-12;

/// Shapes whose closest point of the Minkowski difference is closer to the
/// origin than this fraction of the simplex extent count as touching.
const TOUCH_TOLERANCE: f64 = 1e-12;

/// A vertex of the Minkowski-difference simplex and the support points it was
/// built from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SimplexVertex {
    /// Point of `A - B`: `sa - sb`.
    pub(crate) w: [f64; 3],
    /// Support point on shape `A`.
    pub(crate) sa: [f64; 3],
    /// Support point on shape `B`.
    pub(crate) sb: [f64; 3],
}

impl SimplexVertex {
    /// Placeholder used to fill the unused slots of a fixed-size simplex.
    const ORIGIN: Self = Self {
        w: [0.0; 3],
        sa: [0.0; 3],
        sb: [0.0; 3],
    };
}

/// Tunable parameters of the GJK distance iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DistanceSettings {
    /// Hard cap on the number of iterations (at least one is always performed).
    pub(crate) max_iterations: usize,
    /// Relative tolerance of the termination test
    /// `|v|^2 - v.w <= gap_tolerance |v|^2`.
    pub(crate) gap_tolerance: f64,
}

impl Default for DistanceSettings {
    fn default() -> Self {
        Self {
            max_iterations: MAX_ITERATIONS,
            gap_tolerance: GAP_TOLERANCE,
        }
    }
}

impl DistanceSettings {
    /// Settings with a caller-chosen iteration cap and relative tolerance.
    ///
    /// A tolerance that is negative or not finite falls back to
    /// [`GAP_TOLERANCE`].
    pub(crate) fn new(max_iterations: usize, gap_tolerance: f64) -> Self {
        let gap_tolerance = if gap_tolerance.is_finite() && gap_tolerance >= 0.0 {
            gap_tolerance
        } else {
            GAP_TOLERANCE
        };
        Self {
            max_iterations,
            gap_tolerance,
        }
    }
}

/// Outcome of a GJK distance query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum DistanceOutcome {
    /// The shapes overlap or touch (the origin lies in the Minkowski difference
    /// or closer to it than numerical resolution allows).
    ///
    /// The two points are the witnesses of the final simplex. When the origin is
    /// enclosed they coincide and are a point shared by both shapes.
    Overlapping {
        /// Witness point on shape `A`.
        point_a: [f64; 3],
        /// Witness point on shape `B`.
        point_b: [f64; 3],
    },
    /// The shapes are disjoint.
    Separated {
        /// Euclidean distance between the shapes (an upper bound that is exact
        /// to the relative tolerance of the termination test once converged).
        distance: f64,
        /// Certified lower bound of the distance.
        ///
        /// The plane through the last support point, perpendicular to the
        /// closest point `v`, separates the shapes by `v.w / |v|`. It never
        /// exceeds the true distance and equals `distance` up to the tolerance
        /// when the iteration converged; a conservative advancement must step
        /// on this value, not on `distance`.
        lower_bound: f64,
        /// Closest point on shape `A`.
        point_a: [f64; 3],
        /// Closest point on shape `B`.
        point_b: [f64; 3],
    },
}

/// Result of [`solve_with`]: the outcome together with the final simplex.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct DistanceSolution {
    /// Distance and witnesses, or the overlap verdict.
    pub(crate) outcome: DistanceOutcome,
    /// Final simplex. For a separated pair these are the vertices that support
    /// the closest point; for an overlapping pair they are the simplex that
    /// enclosed (or touched) the origin, a tetrahedron when it was enclosed.
    pub(crate) simplex: Vec<SimplexVertex>,
    /// Number of iterations performed.
    pub(crate) iterations: usize,
}

/// Best feature found so far.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    dist_sq: f64,
    lower_bound: f64,
    point_a: [f64; 3],
    point_b: [f64; 3],
    vertices: [SimplexVertex; 4],
    count: usize,
}

impl Candidate {
    fn into_solution(self, iterations: usize) -> DistanceSolution {
        DistanceSolution {
            outcome: DistanceOutcome::Separated {
                distance: self.dist_sq.sqrt(),
                lower_bound: self.lower_bound,
                point_a: self.point_a,
                point_b: self.point_b,
            },
            simplex: self.vertices[..self.count].to_vec(),
            iterations,
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

/// Verdict for a simplex whose closest point is the origin (or encloses it).
fn overlapping(
    simplex: Vec<SimplexVertex>,
    weights: &[f64; 4],
    iterations: usize,
) -> DistanceSolution {
    let (point_a, point_b) = witnesses(&simplex, weights);
    DistanceSolution {
        outcome: DistanceOutcome::Overlapping { point_a, point_b },
        simplex,
        iterations,
    }
}

/// Run the GJK distance iteration on a Minkowski-difference support function.
///
/// `support(d)` returns the vertex of `A - B` furthest along `d`, together with
/// the support points on `A` (along `d`) and on `B` (along `-d`) it was built
/// from. The first query uses `first_direction`; a zero or non-finite value
/// falls back to `+x`. A good `first_direction` (the last separating axis) is
/// the warm start of coherent queries.
///
/// Overlap is decided by the iteration itself (origin enclosed by the simplex,
/// or the simplex touching the origin); no separate intersection test is used.
/// If the iteration limit or numerical stagnation is reached before the
/// termination test passes, the best (smallest) feature found so far is
/// returned, which is always an upper bound on the true distance.
pub(crate) fn solve_with(
    support: &dyn Fn([f64; 3]) -> SimplexVertex,
    first_direction: [f64; 3],
    settings: &DistanceSettings,
) -> DistanceSolution {
    let direction = if first_direction.iter().all(|c| c.is_finite())
        && dot3(first_direction, first_direction) > 0.0
    {
        first_direction
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut simplex: Vec<SimplexVertex> = Vec::with_capacity(4);
    simplex.push(support(direction));
    let mut best: Option<Candidate> = None;
    let mut iterations = 0;

    for _ in 0..settings.max_iterations.max(1) {
        iterations += 1;
        let count = simplex.len().min(4);
        let mut points = [[0.0_f64; 3]; 4];
        for (slot, vertex) in points.iter_mut().zip(simplex.iter()) {
            *slot = vertex.w;
        }
        let closest = closest_point_to_origin(&points[..count]);
        if closest.encloses_origin {
            return overlapping(simplex, &closest.weights, iterations);
        }

        let v = closest.point;
        let v_sq = dot3(v, v);
        let extent_sq = simplex
            .iter()
            .map(|vertex| dot3(vertex.w, vertex.w))
            .fold(0.0_f64, f64::max);
        if v_sq <= TOUCH_TOLERANCE * TOUCH_TOLERANCE * extent_sq {
            return overlapping(simplex, &closest.weights, iterations);
        }

        // The closest point of a growing hull can only get closer; a larger value
        // is rounding noise, so the previous (better) feature stands.
        if let Some(previous) = best
            && v_sq >= previous.dist_sq
        {
            return previous.into_solution(iterations);
        }
        let (point_a, point_b) = witnesses(&simplex, &closest.weights);

        // Keep only the vertices that support the closest point.
        let mut index = 0;
        simplex.retain(|_| {
            let keep = closest.weights.get(index).is_some_and(|w| *w > 0.0);
            index += 1;
            keep
        });
        let mut vertices = [SimplexVertex::ORIGIN; 4];
        for (slot, vertex) in vertices.iter_mut().zip(simplex.iter()) {
            *slot = *vertex;
        }

        let next = support(neg3(v));
        let along = dot3(v, next.w);
        let candidate = Candidate {
            dist_sq: v_sq,
            lower_bound: (along / v_sq.sqrt()).max(0.0),
            point_a,
            point_b,
            vertices,
            count: simplex.len().min(4),
        };
        best = Some(candidate);

        let gap = v_sq - along;
        let repeated = simplex.iter().any(|vertex| vertex.w == next.w);
        if gap <= settings.gap_tolerance * v_sq || repeated {
            return candidate.into_solution(iterations);
        }
        simplex.push(next);
    }

    match best {
        Some(candidate) => candidate.into_solution(iterations),
        None => overlapping(simplex, &[0.0; 4], iterations),
    }
}

/// Run the GJK distance iteration on two convex shapes.
///
/// See [`solve_with`]; this is the Minkowski difference `A - B` of the two
/// shapes, started along `+x` with the default tolerances.
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
    solve_with(&support, [1.0, 0.0, 0.0], &DistanceSettings::default()).outcome
}
