// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Welzl's minimum enclosing sphere.
//!
//! The textbook formulation recurses once per input point and therefore needs a
//! call stack as deep as the point list. Here the recursion over the point list is
//! unrolled into a loop; only the support set (at most four points on the
//! sphere surface) is handled recursively, so the call depth is bounded by four.

/// Relative slack applied when testing whether a point lies inside a sphere, so
/// that points on the surface (cospherical input, the support points themselves)
/// are not rejected because of round-off.
const SPHERE_CONTAINMENT_SLACK: f64 = 1e-12;

/// A triangle counts as degenerate (collinear) when the sine squared of its
/// largest angle drops below this value.
const TRIANGLE_DEGENERACY: f64 = 1e-16;

/// A tetrahedron counts as degenerate (coplanar) when `|det| / (|u| |v| |w|)` of
/// its edge vectors at the best-conditioned corner drops below this value.
const TETRAHEDRON_DEGENERACY: f64 = 1e-10;

/// Seed of the pseudo-random permutation applied by [`minimum_enclosing_sphere`].
const WELZL_SHUFFLE_SEED: u64 = 0x9E37_79B9_7F4A_7C15;

/// Minimum enclosing sphere of `points` as `(center, radius)`.
///
/// An empty input gives the radius-zero sphere at the origin. The points are
/// processed in a pseudo-random order (a copy of the input permuted with a
/// fixed-seed LCG), which gives the expected linear running time of the
/// randomised algorithm and keeps the result deterministic.
pub(super) fn minimum_enclosing_sphere(points: &[[f64; 3]]) -> ([f64; 3], f64) {
    let mut pts: Vec<[f64; 3]> = points.to_vec();
    shuffle_points(&mut pts, WELZL_SHUFFLE_SEED);
    let mut support = [[0.0; 3]; 4];
    let n = pts.len();
    let Some((center, _)) = welzl_sphere(&pts, &mut support, 0, n) else {
        return ([0.0; 3], 0.0);
    };
    // Re-measure the radius so that round-off in the support sphere can never
    // leave an input point outside.
    let radius = points
        .iter()
        .map(|&p| sphere_dist3(p, center))
        .fold(0.0, f64::max);
    (center, radius)
}

/// Smallest sphere that has every point of `boundary` on its surface.
///
/// `boundary` holds at most four points (extra points are ignored): one point
/// gives the radius-zero sphere, two points the sphere on their diameter, three
/// points the sphere whose centre lies in their plane and four points the
/// circumsphere of the tetrahedron. An empty slice gives the radius-zero sphere
/// at the origin.
fn min_sphere_with_boundary(boundary: &[[f64; 3]]) -> ([f64; 3], f64) {
    match boundary {
        [] => ([0.0; 3], 0.0),
        [a] => (*a, 0.0),
        [a, b] => diametral_sphere(*a, *b),
        [a, b, c] => circumsphere_3(*a, *b, *c),
        [a, b, c, d, ..] => circumsphere_4(*a, *b, *c, *d),
    }
}

fn sphere_dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Sphere that has the segment `ab` as a diameter.
fn diametral_sphere(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], f64) {
    (midpoint3sph(a, b), 0.5 * sphere_dist3(a, b))
}

/// Whether `p` lies inside the sphere `(centre, radius)`, up to
/// [`SPHERE_CONTAINMENT_SLACK`].
fn sphere_contains(sphere: ([f64; 3], f64), p: [f64; 3]) -> bool {
    let (center, radius) = sphere;
    sphere_dist3(p, center) <= radius * (1.0 + SPHERE_CONTAINMENT_SLACK)
}

/// Welzl's algorithm with a bounded support set.
///
/// Returns the smallest sphere that encloses the first `n` points of `pts` and
/// has the first `support_len` entries of `support` on its surface, or `None`
/// when there is neither a point nor a support entry. The recursion over the
/// point list of the textbook formulation is unrolled into a loop, so the call
/// depth is at most four regardless of the number of points.
fn welzl_sphere(
    pts: &[[f64; 3]],
    support: &mut [[f64; 3]; 4],
    support_len: usize,
    n: usize,
) -> Option<([f64; 3], f64)> {
    let mut sphere = if support_len == 0 {
        None
    } else {
        Some(min_sphere_with_boundary(&support[..support_len]))
    };
    if support_len == support.len() {
        return sphere;
    }
    for (i, &p) in pts.iter().enumerate().take(n) {
        if sphere.is_some_and(|s| sphere_contains(s, p)) {
            continue;
        }
        support[support_len] = p;
        sphere = welzl_sphere(pts, support, support_len + 1, i);
    }
    sphere
}

/// Deterministic Fisher-Yates shuffle driven by a 64-bit LCG.
fn shuffle_points(pts: &mut [[f64; 3]], seed: u64) {
    let mut state = seed;
    for i in (1..pts.len()).rev() {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        // The high bits of an LCG are the well distributed ones.
        let j = ((state >> 11) % (i as u64 + 1)) as usize;
        pts.swap(i, j);
    }
}

/// Circumscribed sphere of the triangle `abc`, whose centre lies in the plane
/// of the triangle. Returns `None` when the triangle is degenerate (collinear
/// or coincident vertices).
fn try_circumsphere_3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<([f64; 3], f64)> {
    // Use the vertex opposite the longest edge as the apex: its angle is the
    // largest one, so the cross product below is well conditioned unless the
    // three vertices are (nearly) collinear.
    let ab_sq = dot3sph(sub3sph(b, a), sub3sph(b, a));
    let ac_sq = dot3sph(sub3sph(c, a), sub3sph(c, a));
    let bc_sq = dot3sph(sub3sph(c, b), sub3sph(c, b));
    let (a, b, c, ab_sq, ac_sq) = if bc_sq >= ab_sq && bc_sq >= ac_sq {
        (a, b, c, ab_sq, ac_sq)
    } else if ac_sq >= ab_sq {
        (b, c, a, bc_sq, ab_sq)
    } else {
        (c, a, b, ac_sq, bc_sq)
    };
    let ab = sub3sph(b, a);
    let ac = sub3sph(c, a);
    let normal = cross3sph(ab, ac);
    let normal_sq = dot3sph(normal, normal);
    // |ab x ac|^2 = |ab|^2 |ac|^2 sin^2(angle at a); non-finite input is rejected
    // as well.
    if !normal_sq.is_finite() || normal_sq <= TRIANGLE_DEGENERACY * ab_sq * ac_sq {
        return None;
    }
    // centre = a + (|ac|^2 (n x ab) + |ab|^2 (ac x n)) / (2 |n|^2) with n = ab x ac
    let first = scale3sph(cross3sph(normal, ab), ac_sq);
    let second = scale3sph(cross3sph(ac, normal), ab_sq);
    let k = 0.5 / normal_sq;
    let center = [
        a[0] + (first[0] + second[0]) * k,
        a[1] + (first[1] + second[1]) * k,
        a[2] + (first[2] + second[2]) * k,
    ];
    let radius = sphere_dist3(a, center)
        .max(sphere_dist3(b, center))
        .max(sphere_dist3(c, center));
    Some((center, radius))
}

/// Edge vectors from corner `apex` of a tetrahedron to the three other corners.
fn tetrahedron_edges(corners: &[[f64; 3]; 4], apex: usize) -> [[f64; 3]; 3] {
    let origin = corners[apex % 4];
    [1, 2, 3].map(|step| sub3sph(corners[(apex + step) % 4], origin))
}

/// `|u . (v x w)| / (|u| |v| |w|)`: the sine-like measure of how far the three
/// edge vectors are from being coplanar.
fn relative_triple_product(edges: &[[f64; 3]; 3]) -> f64 {
    let [u, v, w] = *edges;
    let scale = dot3sph(u, u).sqrt() * dot3sph(v, v).sqrt() * dot3sph(w, w).sqrt();
    dot3sph(u, cross3sph(v, w)).abs() / scale
}

/// Circumscribed sphere of the tetrahedron `abcd`. Returns `None` when the four
/// points are (nearly) coplanar.
fn try_circumsphere_4(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<([f64; 3], f64)> {
    let corners = [a, b, c, d];
    // Cramer's rule below is best conditioned at the corner whose edge vectors
    // span the most volume relative to their lengths.
    let mut best: Option<(usize, f64)> = None;
    for apex in 0..corners.len() {
        let measure = relative_triple_product(&tetrahedron_edges(&corners, apex));
        if measure.is_finite() && best.is_none_or(|(_, current)| measure > current) {
            best = Some((apex, measure));
        }
    }
    let (apex, measure) = best?;
    if measure <= TETRAHEDRON_DEGENERACY {
        return None;
    }
    let origin = corners[apex];
    let [u, v, w] = tetrahedron_edges(&corners, apex);
    // 2 [u; v; w] x = [|u|^2; |v|^2; |w|^2] solved with Cramer's rule gives
    // x = (|u|^2 (v x w) + |v|^2 (w x u) + |w|^2 (u x v)) / (2 u . (v x w)).
    let det = dot3sph(u, cross3sph(v, w));
    let tu = scale3sph(cross3sph(v, w), dot3sph(u, u));
    let tv = scale3sph(cross3sph(w, u), dot3sph(v, v));
    let tw = scale3sph(cross3sph(u, v), dot3sph(w, w));
    let k = 0.5 / det;
    let center = [
        origin[0] + (tu[0] + tv[0] + tw[0]) * k,
        origin[1] + (tu[1] + tv[1] + tw[1]) * k,
        origin[2] + (tu[2] + tv[2] + tw[2]) * k,
    ];
    let radius = corners
        .iter()
        .map(|&p| sphere_dist3(p, center))
        .fold(0.0, f64::max);
    Some((center, radius))
}

/// Smallest sphere enclosing a handful of degenerate (collinear or coplanar)
/// points.
///
/// The centre of a minimum enclosing sphere is the midpoint of two of the points
/// or the circumcentre of three of them, so every such candidate is tried and
/// the one with the smallest enclosing radius wins.
fn enclosing_sphere_of_candidates(points: &[[f64; 3]]) -> ([f64; 3], f64) {
    let mut best: Option<([f64; 3], f64)> = None;
    let mut consider = |center: [f64; 3]| {
        let radius = points
            .iter()
            .map(|&p| sphere_dist3(p, center))
            .fold(0.0, f64::max);
        if best.is_none_or(|(_, current)| radius < current) {
            best = Some((center, radius));
        }
    };
    for (i, &p) in points.iter().enumerate() {
        for (j, &q) in points.iter().enumerate().skip(i + 1) {
            consider(midpoint3sph(p, q));
            for &r in points.iter().skip(j + 1) {
                if let Some((center, _)) = try_circumsphere_3(p, q, r) {
                    consider(center);
                }
            }
        }
    }
    best.unwrap_or_else(|| (points.first().copied().unwrap_or([0.0; 3]), 0.0))
}

/// Sphere through the three points `a`, `b`, `c` with its centre in their plane
/// (the circumsphere of the triangle).
///
/// Collinear points cannot lie on a common sphere; in that case the smallest
/// sphere enclosing all three is returned instead.
fn circumsphere_3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> ([f64; 3], f64) {
    try_circumsphere_3(a, b, c).unwrap_or_else(|| enclosing_sphere_of_candidates(&[a, b, c]))
}

/// Circumsphere of the tetrahedron `abcd`: the unique sphere through all four
/// points.
///
/// Coplanar points have no such sphere unless they are concyclic; in that case
/// the smallest sphere enclosing all four is returned instead.
fn circumsphere_4(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> ([f64; 3], f64) {
    try_circumsphere_4(a, b, c, d).unwrap_or_else(|| enclosing_sphere_of_candidates(&[a, b, c, d]))
}

fn midpoint3sph(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        (a[0] + b[0]) * 0.5,
        (a[1] + b[1]) * 0.5,
        (a[2] + b[2]) * 0.5,
    ]
}
fn sub3sph(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross3sph(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot3sph(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scale3sph(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[cfg(test)]
mod tests;
