// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact distance between the convex hulls of two finite point sets.
//!
//! This is the Gilbert-Johnson-Keerthi algorithm specialised to polytopes given
//! by their vertices. The support function of such a shape is an exact `argmax`
//! over the vertex list, the Minkowski difference of the two hulls is the hull
//! of all pairwise differences, and the algorithm terminates after finitely many
//! iterations with the exact closest pair of points (up to floating-point
//! rounding), including for degenerate hulls (a single point, a segment, a flat
//! polygon).
//!
//! The closest point of the current simplex (at most four vertices of the
//! Minkowski difference) to the origin is found by examining every sub-simplex
//! and keeping the one whose closest point passes the optimality certificate
//! `c . q >= c . c` for all simplex vertices `q`. Termination follows van den
//! Bergen: the iteration stops when the support point in direction `-c` no
//! longer improves on the closest point `c` of the current simplex.

use super::{add3, dot3, len3, scale3, sub3};

/// Upper bound on GJK iterations. A polytope with `n` vertices needs far fewer.
const MAX_ITERATIONS: usize = 128;

/// One vertex of the Minkowski difference `A - B`, remembering its sources.
#[derive(Clone, Copy)]
struct CsoVertex {
    /// Index of the support point in shape A.
    index_a: usize,
    /// Index of the support point in shape B.
    index_b: usize,
    /// Support point of A (already translated).
    point_a: [f64; 3],
    /// Support point of B.
    point_b: [f64; 3],
    /// `point_a - point_b`.
    point: [f64; 3],
}

/// Distance query result for two convex hulls.
#[derive(Debug, Clone)]
pub(super) struct HullDistance {
    /// Euclidean distance between the hulls (zero if they intersect).
    pub distance: f64,
    /// Closest point on hull A (a point common to both hulls if they intersect).
    pub point_a: [f64; 3],
    /// Closest point on hull B (a point common to both hulls if they intersect).
    pub point_b: [f64; 3],
    /// True if the hulls intersect or touch.
    pub intersecting: bool,
}

/// Convex combination behind the witnesses of a [`HullDistance`], as
/// `(index in A, index in B, weight)`. The weights are positive and sum to one;
/// `point_a = sum(w * (A[index_a] + offset))` and `point_b = sum(w * B[index_b])`.
pub(super) type SupportCombination = Vec<(usize, usize, f64)>;

/// Index of the first point of `points` that maximises `p . dir`.
fn support_index(points: &[[f64; 3]], dir: [f64; 3]) -> usize {
    points
        .iter()
        .enumerate()
        .fold((0usize, f64::NEG_INFINITY), |(best, best_value), (i, p)| {
            let value = dot3(*p, dir);
            if value > best_value {
                (i, value)
            } else {
                (best, best_value)
            }
        })
        .0
}

/// The point of `points` that is extreme in direction `dir` (`None` if empty).
pub(super) fn support_point(points: &[[f64; 3]], dir: [f64; 3]) -> Option<[f64; 3]> {
    points.get(support_index(points, dir)).copied()
}

/// Mean of the points (`None` if empty).
pub(super) fn centroid(points: &[[f64; 3]]) -> Option<[f64; 3]> {
    if points.is_empty() {
        return None;
    }
    let sum = points.iter().fold([0.0; 3], |acc, p| add3(acc, *p));
    Some(scale3(sum, 1.0 / points.len() as f64))
}

/// Solve the Gram system `G mu = rhs` of size `k` (1 to 3) for the barycentric
/// offsets of a sub-simplex. Returns `None` for a (numerically) singular system,
/// i.e. an affinely dependent sub-simplex, which a smaller sub-simplex covers.
fn solve_gram(g: &[[f64; 3]; 3], rhs: [f64; 3], k: usize) -> Option<[f64; 3]> {
    const SINGULAR: f64 = 1e-20;
    match k {
        1 => {
            if g[0][0].is_nan() || g[0][0] <= 0.0 {
                return None;
            }
            Some([rhs[0] / g[0][0], 0.0, 0.0])
        }
        2 => {
            let det = g[0][0] * g[1][1] - g[0][1] * g[0][1];
            if det.is_nan() || det <= SINGULAR * g[0][0] * g[1][1] {
                return None;
            }
            Some([
                (rhs[0] * g[1][1] - g[0][1] * rhs[1]) / det,
                (g[0][0] * rhs[1] - g[0][1] * rhs[0]) / det,
                0.0,
            ])
        }
        3 => {
            let a00 = g[1][1] * g[2][2] - g[1][2] * g[1][2];
            let a01 = g[0][2] * g[1][2] - g[0][1] * g[2][2];
            let a02 = g[0][1] * g[1][2] - g[0][2] * g[1][1];
            let a11 = g[0][0] * g[2][2] - g[0][2] * g[0][2];
            let a12 = g[0][1] * g[0][2] - g[0][0] * g[1][2];
            let a22 = g[0][0] * g[1][1] - g[0][1] * g[0][1];
            let det = g[0][0] * a00 + g[0][1] * a01 + g[0][2] * a02;
            if det.is_nan() || det <= SINGULAR * g[0][0] * g[1][1] * g[2][2] {
                return None;
            }
            Some([
                (a00 * rhs[0] + a01 * rhs[1] + a02 * rhs[2]) / det,
                (a01 * rhs[0] + a11 * rhs[1] + a12 * rhs[2]) / det,
                (a02 * rhs[0] + a12 * rhs[1] + a22 * rhs[2]) / det,
            ])
        }
        _ => None,
    }
}

/// Barycentric weights (aligned with `verts`) of the point of the affine hull of
/// the sub-simplex `members` that is closest to the origin, or `None` when the
/// sub-simplex is affinely dependent or the point lies outside the sub-simplex.
fn subset_weights(verts: &[CsoVertex], members: &[usize]) -> Option<[f64; 4]> {
    let first = verts.get(*members.first()?)?.point;
    let k = members.len() - 1;
    let mut edges = [[0.0; 3]; 3];
    for (edge, &m) in edges.iter_mut().zip(members.iter().skip(1)) {
        *edge = sub3(verts.get(m)?.point, first);
    }
    let mut gram = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for (i, edge_i) in edges.iter().enumerate().take(k) {
        rhs[i] = -dot3(*edge_i, first);
        for (j, edge_j) in edges.iter().enumerate().take(k) {
            gram[i][j] = dot3(*edge_i, *edge_j);
        }
    }
    let mut weights = [0.0; 4];
    if k == 0 {
        *weights.get_mut(members[0])? = 1.0;
        return Some(weights);
    }
    let mu = solve_gram(&gram, rhs, k)?;
    let mut sum = 0.0;
    for (slot, &m) in members.iter().enumerate().skip(1) {
        *weights.get_mut(m)? = mu[slot - 1];
        sum += mu[slot - 1];
    }
    *weights.get_mut(members[0])? = 1.0 - sum;
    // Weights that are negative beyond rounding put the point outside the face.
    if weights.iter().any(|&w| w < -1e-12) {
        return None;
    }
    let mut total = 0.0;
    for w in weights.iter_mut() {
        *w = w.max(0.0);
        total += *w;
    }
    if total.is_nan() || total <= 0.0 {
        return None;
    }
    for w in weights.iter_mut() {
        *w /= total;
    }
    Some(weights)
}

/// Closest point to the origin of the convex hull of `verts` (one to four
/// vertices), with barycentric weights aligned with `verts`.
fn closest_to_origin(verts: &[CsoVertex]) -> ([f64; 3], [f64; 4]) {
    let n = verts.len().min(4);
    let verts = &verts[..n];
    let max_norm = verts.iter().map(|v| len3(v.point)).fold(0.0, f64::max);

    let mut best: Option<([f64; 3], [f64; 4], f64)> = None;
    // Sub-simplices in order of increasing size, so that of two equally good
    // candidates the one with fewer vertices wins.
    for size in 1..=n {
        for mask in 1u32..(1u32 << n) {
            if mask.count_ones() as usize != size {
                continue;
            }
            let mut members = [0usize; 4];
            let mut count = 0;
            for i in 0..n {
                if mask & (1 << i) != 0 {
                    members[count] = i;
                    count += 1;
                }
            }
            let Some(weights) = subset_weights(verts, &members[..count]) else {
                continue;
            };
            let mut point = [0.0; 3];
            for (v, &w) in verts.iter().zip(weights.iter()) {
                point = add3(point, scale3(v.point, w));
            }
            let point_sq = dot3(point, point);
            // Optimality certificate: no simplex vertex lies behind the plane
            // through `point` perpendicular to it.
            let slack = 1e-12 * (len3(point) * max_norm) + 1e-300;
            if verts.iter().any(|v| dot3(v.point, point) < point_sq - slack) {
                continue;
            }
            let better = match best {
                None => true,
                Some((_, _, best_sq)) => point_sq < best_sq * (1.0 - 1e-12),
            };
            if better {
                best = Some((point, weights, point_sq));
            }
        }
    }
    if let Some((point, weights, _)) = best {
        return (point, weights);
    }
    // Numerical fallback: the vertex closest to the origin.
    let nearest = verts
        .iter()
        .enumerate()
        .min_by(|a, b| dot3(a.1.point, a.1.point).total_cmp(&dot3(b.1.point, b.1.point)));
    let mut weights = [0.0; 4];
    match nearest {
        Some((index, vertex)) => {
            weights[index] = 1.0;
            (vertex.point, weights)
        }
        None => ([0.0; 3], weights),
    }
}

/// Assemble the witnesses of the active vertices `(vertex, weight)`.
fn assemble(
    active: &[(CsoVertex, f64)],
    closest_sq: f64,
    intersecting: bool,
) -> (HullDistance, SupportCombination) {
    let mut point_a = [0.0; 3];
    let mut point_b = [0.0; 3];
    let mut combination = Vec::with_capacity(active.len());
    for &(v, w) in active {
        point_a = add3(point_a, scale3(v.point_a, w));
        point_b = add3(point_b, scale3(v.point_b, w));
        combination.push((v.index_a, v.index_b, w));
    }
    let result = HullDistance {
        distance: if intersecting { 0.0 } else { closest_sq.sqrt() },
        point_a,
        point_b,
        intersecting,
    };
    (result, combination)
}

/// Like [`hull_distance`], additionally returning the convex combination of
/// support points behind the witnesses.
pub(super) fn hull_distance_with_support(
    points_a: &[[f64; 3]],
    offset_a: [f64; 3],
    points_b: &[[f64; 3]],
) -> Option<(HullDistance, SupportCombination)> {
    let centroid_a = add3(centroid(points_a)?, offset_a);
    let centroid_b = centroid(points_b)?;

    let support = |dir: [f64; 3]| -> Option<CsoVertex> {
        let index_a = support_index(points_a, dir);
        let index_b = support_index(points_b, scale3(dir, -1.0));
        let point_a = add3(*points_a.get(index_a)?, offset_a);
        let point_b = *points_b.get(index_b)?;
        Some(CsoVertex {
            index_a,
            index_b,
            point_a,
            point_b,
            point: sub3(point_a, point_b),
        })
    };

    let mut direction = sub3(centroid_b, centroid_a);
    if !direction.iter().all(|c| c.is_finite()) || dot3(direction, direction) <= 0.0 {
        direction = [1.0, 0.0, 0.0];
    }

    let mut simplex: Vec<CsoVertex> = Vec::with_capacity(4);
    simplex.push(support(direction)?);
    let mut scale = len3(simplex[0].point);
    let mut last: Option<(HullDistance, SupportCombination)> = None;

    for _ in 0..MAX_ITERATIONS {
        let (closest, weights) = closest_to_origin(&simplex);
        // Keep only the vertices that carry the closest point.
        let mut active: Vec<(CsoVertex, f64)> = simplex
            .iter()
            .zip(weights.iter())
            .filter(|&(_, &w)| w > 0.0)
            .map(|(v, &w)| (*v, w))
            .collect();
        if active.is_empty() {
            active.push((simplex[0], 1.0));
        }
        simplex = active.iter().map(|&(v, _)| v).collect();

        let closest_sq = dot3(closest, closest);
        let zero = 1e-12 * scale.max(1e-300);
        if closest_sq <= zero * zero {
            return Some(assemble(&active, closest_sq, true));
        }
        last = Some(assemble(&active, closest_sq, false));

        let candidate = support(scale3(closest, -1.0))?;
        scale = scale.max(len3(candidate.point));
        let gap = closest_sq - dot3(closest, candidate.point);
        let tolerance =
            1e-13 * closest_sq + 16.0 * f64::EPSILON * len3(closest) * len3(candidate.point);
        let repeated = simplex
            .iter()
            .any(|v| v.index_a == candidate.index_a && v.index_b == candidate.index_b);
        if gap <= tolerance || repeated || simplex.len() >= 4 {
            return last;
        }
        simplex.push(candidate);
    }
    last
}

/// Distance between the convex hull of `points_a`, translated by `offset_a`, and
/// the convex hull of `points_b`.
///
/// Returns `None` if either point set is empty.
pub(super) fn hull_distance(
    points_a: &[[f64; 3]],
    offset_a: [f64; 3],
    points_b: &[[f64; 3]],
) -> Option<HullDistance> {
    hull_distance_with_support(points_a, offset_a, points_b).map(|(result, _)| result)
}
