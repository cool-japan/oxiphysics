// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lower-bound GJK distance between two convex shapes given by support
//! functions.
//!
//! The distance of two convex shapes is the distance of the origin from their
//! Minkowski difference `C = A - B`. GJK (Gilbert, Johnson and Keerthi, 1988)
//! keeps a simplex of points of `C` and repeats two steps:
//!
//! 1. take the point `v` of the simplex hull that is nearest to the origin
//!    (and drop the simplex vertices that do not carry it);
//! 2. query the support point `w` of `C` in the direction `-v`.
//!
//! `|v|` is an *upper* bound of the distance, because the simplex lies inside
//! `C`. Every point `x` of `C` satisfies `x . v >= w . v`, so the plane through
//! `w` with normal `v` separates `C` from the origin by `v . w / |v|`: that is
//! a *lower* bound of the distance. The iteration ends when the two bounds
//! meet, and the lower bound is returned: conservative advancement advances
//! time by `distance / speed`, which is only safe if the distance is never
//! overestimated.

use super::{add3, cross3, dot3, neg3, scale3, sub3};

/// Maximum number of support queries of one distance computation.
const MAX_ITERATIONS: usize = 64;

/// Relative part of the gap `|v| - lower_bound` below which the distance counts
/// as converged.
const RELATIVE_TOLERANCE: f64 = 1e-9;

/// Absolute part of the converged gap. Round-off in the nearest simplex point
/// limits the lower bound of a very small distance between large shapes to
/// about `epsilon * size^2 / distance`, so a purely relative test would never
/// be met there.
const ABSOLUTE_TOLERANCE: f64 = 1e-10;

/// Squared distance of the nearest simplex point below which the origin counts
/// as lying on the simplex, i.e. the shapes touch or overlap.
const TOUCHING_DISTANCE_SQ: f64 = 1e-24;

/// Squared edge length below which a segment counts as a single point.
const DEGENERATE_LENGTH_SQ: f64 = 1e-300;

/// Distance between the convex shapes `A` and `B` given by their support
/// functions (`support_a(d)` is the point of `A` farthest along `d`).
///
/// The result never exceeds the true distance and converges to it: iteration
/// stops when the upper and lower bounds of the distance are within `1e-10`
/// plus a relative `1e-9` of each other, when round-off stalls the simplex, or
/// after 64 support queries, in which case the best lower bound found so far is
/// returned. It is `0.0` when the shapes touch or overlap.
pub fn gjk_distance_approx(
    support_a: impl Fn([f64; 3]) -> [f64; 3],
    support_b: impl Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    // Support point of the Minkowski difference A - B along `d`.
    let support = |d: [f64; 3]| -> [f64; 3] { sub3(support_a(d), support_b(neg3(d))) };

    let first = support([1.0, 0.0, 0.0]);
    let mut simplex = Simplex::new(first);
    let mut nearest = first;
    let mut lower_bound = 0.0_f64;

    for _ in 0..MAX_ITERATIONS {
        let distance_sq = dot3(nearest, nearest);
        if distance_sq <= TOUCHING_DISTANCE_SQ {
            // The origin lies on the simplex hull: touching or overlapping.
            return 0.0;
        }
        let distance = distance_sq.sqrt();
        let toward_origin = scale3(nearest, -1.0 / distance);
        let candidate = support(toward_origin);

        // Plane through `candidate` with normal `nearest`: separates C from the
        // origin by `candidate . nearest / |nearest|`.
        lower_bound = lower_bound.max(-dot3(candidate, toward_origin));
        let converged = distance - lower_bound <= ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * distance;
        if converged || simplex.contains(&candidate) {
            return lower_bound;
        }

        simplex.push(candidate);
        match simplex.reduce_to_nearest() {
            // Adding a point below the separating plane must bring the nearest
            // point closer; if it does not, round-off has stalled the iteration.
            Some(point) if dot3(point, point) >= distance_sq => return lower_bound,
            Some(point) => nearest = point,
            // The tetrahedron encloses the origin.
            None => return 0.0,
        }
    }
    lower_bound
}

/// Up to four points of the Minkowski difference; the newest point is last.
struct Simplex {
    points: [[f64; 3]; 4],
    len: usize,
}

/// Vertices (indices into the points that were examined) whose hull holds the
/// point nearest to the origin: a vertex, an edge or a face.
#[derive(Clone, Copy)]
struct Feature {
    indices: [usize; 3],
    len: usize,
}

impl Feature {
    fn vertex(i: usize) -> Self {
        Self {
            indices: [i, 0, 0],
            len: 1,
        }
    }

    fn edge(i: usize, j: usize) -> Self {
        Self {
            indices: [i, j, 0],
            len: 2,
        }
    }

    fn face(i: usize, j: usize, k: usize) -> Self {
        Self {
            indices: [i, j, k],
            len: 3,
        }
    }

    /// Re-index the feature through `map` (local index to outer index).
    fn mapped(self, map: &[usize]) -> Self {
        let mut indices = [0; 3];
        for (slot, &local) in indices.iter_mut().zip(&self.indices[..self.len]) {
            *slot = map[local];
        }
        Self {
            indices,
            len: self.len,
        }
    }
}

impl Simplex {
    fn new(first: [f64; 3]) -> Self {
        let mut points = [[0.0; 3]; 4];
        points[0] = first;
        Self { points, len: 1 }
    }

    /// Whether `point` is already a vertex of the simplex.
    fn contains(&self, point: &[f64; 3]) -> bool {
        self.points[..self.len].iter().any(|p| p == point)
    }

    /// Append `point` as the newest vertex.
    fn push(&mut self, point: [f64; 3]) {
        if self.len < self.points.len() {
            self.points[self.len] = point;
            self.len += 1;
        }
    }

    /// Replace the simplex by the vertices of `feature`, preserving their order.
    fn keep(&mut self, feature: Feature) {
        let mut kept = [[0.0; 3]; 4];
        for (slot, &index) in kept.iter_mut().zip(&feature.indices[..feature.len]) {
            *slot = self.points[index];
        }
        self.points = kept;
        self.len = feature.len;
    }

    /// Reduce the simplex to the feature that holds the point nearest to the
    /// origin and return that point, or `None` when a tetrahedron encloses the
    /// origin.
    fn reduce_to_nearest(&mut self) -> Option<[f64; 3]> {
        let p = self.points;
        let (point, feature) = match self.len {
            1 => return Some(p[0]),
            2 => closest_on_segment(p[0], p[1]),
            3 => closest_on_triangle(p[0], p[1], p[2]),
            _ => closest_on_tetrahedron(&p)?,
        };
        self.keep(feature);
        Some(point)
    }
}

/// Point of the segment `a`-`b` nearest to the origin.
fn closest_on_segment(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], Feature) {
    let ab = sub3(b, a);
    let length_sq = dot3(ab, ab);
    if length_sq <= DEGENERATE_LENGTH_SQ {
        return (b, Feature::vertex(1));
    }
    let t = -dot3(a, ab) / length_sq;
    if t <= 0.0 {
        (a, Feature::vertex(0))
    } else if t >= 1.0 {
        (b, Feature::vertex(1))
    } else {
        (add3(a, scale3(ab, t)), Feature::edge(0, 1))
    }
}

/// Point of the triangle `a`, `b`, `c` nearest to the origin (Ericson, Real-Time
/// Collision Detection, 5.1.5, with the query point at the origin).
fn closest_on_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], Feature) {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ao = neg3(a);
    let d1 = dot3(ab, ao);
    let d2 = dot3(ac, ao);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, Feature::vertex(0));
    }

    let bo = neg3(b);
    let d3 = dot3(ab, bo);
    let d4 = dot3(ac, bo);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, Feature::vertex(1));
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (add3(a, scale3(ab, v)), Feature::edge(0, 1));
    }

    let co = neg3(c);
    let d5 = dot3(ab, co);
    let d6 = dot3(ac, co);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, Feature::vertex(2));
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (add3(a, scale3(ac, w)), Feature::edge(0, 2));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (add3(b, scale3(sub3(c, b), w)), Feature::edge(1, 2));
    }

    let sum = va + vb + vc;
    if sum.abs() <= DEGENERATE_LENGTH_SQ {
        // Collinear vertices: the triangle is a segment, take the best edge.
        return nearest_of_edges(&[a, b, c], &[(0, 1), (0, 2), (1, 2)]);
    }
    let denom = 1.0 / sum;
    let v = vb * denom;
    let w = vc * denom;
    (
        add3(a, add3(scale3(ab, v), scale3(ac, w))),
        Feature::face(0, 1, 2),
    )
}

/// The closest point to the origin among several segments of `points`, given as
/// index pairs.
fn nearest_of_edges(points: &[[f64; 3]], edges: &[(usize, usize)]) -> ([f64; 3], Feature) {
    let mut best: Option<(f64, [f64; 3], Feature)> = None;
    for &(i, j) in edges {
        let (point, feature) = closest_on_segment(points[i], points[j]);
        let feature = feature.mapped(&[i, j]);
        let distance_sq = dot3(point, point);
        if best.is_none_or(|(best_sq, _, _)| distance_sq < best_sq) {
            best = Some((distance_sq, point, feature));
        }
    }
    match best {
        Some((_, point, feature)) => (point, feature),
        None => (points[0], Feature::vertex(0)),
    }
}

/// Point of the tetrahedron `p` nearest to the origin, or `None` when the origin
/// lies inside it.
fn closest_on_tetrahedron(p: &[[f64; 3]; 4]) -> Option<([f64; 3], Feature)> {
    const FACES: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    // Vertex of the tetrahedron opposite each face of `FACES`.
    const OPPOSITE: [usize; 4] = [3, 2, 1, 0];

    let extent_sq = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
        .map(|(i, j)| {
            let edge = sub3(p[j], p[i]);
            dot3(edge, edge)
        })
        .fold(0.0_f64, f64::max);
    let volume6 = dot3(sub3(p[1], p[0]), cross3(sub3(p[2], p[0]), sub3(p[3], p[0])));
    // A (nearly) flat tetrahedron cannot enclose the origin.
    let flat = volume6.abs() <= 1e-12 * extent_sq * extent_sq.sqrt();

    if !flat {
        let enclosed = FACES.iter().zip(OPPOSITE).all(|(face, opposite)| {
            let [a, b, c] = [p[face[0]], p[face[1]], p[face[2]]];
            let normal = cross3(sub3(b, a), sub3(c, a));
            // The origin is on the same side of the face as the opposite vertex.
            dot3(normal, neg3(a)) * dot3(normal, sub3(p[opposite], a)) >= 0.0
        });
        if enclosed {
            return None;
        }
    }

    let mut best: Option<(f64, [f64; 3], Feature)> = None;
    for face in &FACES {
        let (point, feature) = closest_on_triangle(p[face[0]], p[face[1]], p[face[2]]);
        let feature = feature.mapped(face);
        let distance_sq = dot3(point, point);
        if best.is_none_or(|(best_sq, _, _)| distance_sq < best_sq) {
            best = Some((distance_sq, point, feature));
        }
    }
    best.map(|(_, point, feature)| (point, feature))
}
