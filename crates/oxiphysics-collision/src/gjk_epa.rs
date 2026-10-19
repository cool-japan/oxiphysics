// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GJK (Gilbert-Johnson-Keerthi) and EPA (Expanding Polytope Algorithm)
//! for convex shape collision detection using plain `[f64; 3]` arrays.

use std::collections::BTreeMap;

pub(crate) mod distance;
#[cfg(test)]
mod distance_tests;
#[cfg(test)]
mod epa_contact_tests;
pub(crate) mod simplex_closest;
#[cfg(test)]
pub(crate) mod test_support;

// ─── Helper math ────────────────────────────────────────────────────────────

/// Dot product of two 3-vectors stored as `[f64; 3]`.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        scale3(a, 1.0 / l)
    }
}

#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ─── ConvexShape trait ───────────────────────────────────────────────────────

/// A convex shape that exposes a support function used by GJK/EPA.
pub trait ConvexShape {
    /// Returns the furthest point on the shape in the given direction.
    fn support(&self, direction: [f64; 3]) -> [f64; 3];
}

// ─── GjkSphere ───────────────────────────────────────────────────────────────

/// A sphere represented by center and radius, used with GJK/EPA.
pub struct GjkSphere {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}

impl ConvexShape for GjkSphere {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let n = normalize3(direction);
        add3(self.center, scale3(n, self.radius))
    }
}

// ─── GjkBox ──────────────────────────────────────────────────────────────────

/// An axis-aligned box represented by center and half-extents, used with GJK/EPA.
pub struct GjkBox {
    /// Center of the box.
    pub center: [f64; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}

impl ConvexShape for GjkBox {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        [
            self.center[0] + self.half_extents[0] * direction[0].signum(),
            self.center[1] + self.half_extents[1] * direction[1].signum(),
            self.center[2] + self.half_extents[2] * direction[2].signum(),
        ]
    }
}

// ─── GjkCapsule ──────────────────────────────────────────────────────────────

/// A capsule (cylinder with hemispherical caps) used with GJK/EPA.
pub struct GjkCapsule {
    /// One end of the capsule axis.
    pub p1: [f64; 3],
    /// Other end of the capsule axis.
    pub p2: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}

impl ConvexShape for GjkCapsule {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let n = normalize3(direction);
        let d1 = dot3(self.p1, n);
        let d2 = dot3(self.p2, n);
        let base = if d1 > d2 { self.p1 } else { self.p2 };
        add3(base, scale3(n, self.radius))
    }
}

// ─── Minkowski difference support ────────────────────────────────────────────

/// Support function for the Minkowski difference of two convex shapes.
pub fn minkowski_support(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    dir: [f64; 3],
) -> [f64; 3] {
    sub3(shape_a.support(dir), shape_b.support(neg3(dir)))
}

/// A point of the Minkowski difference `A - B` together with the support
/// points on `A` and `B` it is the difference of.
///
/// Keeping the pair lets EPA turn the barycentric coordinates of a face into
/// contact points on the two shapes: `sum(l_i * a_i)` and `sum(l_i * b_i)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportPoint {
    /// Point of the Minkowski difference: `a - b`.
    pub w: [f64; 3],
    /// Support point on shape A.
    pub a: [f64; 3],
    /// Support point on shape B.
    pub b: [f64; 3],
}

/// Support point of the Minkowski difference of two convex shapes along `dir`,
/// together with the support points on both shapes.
pub fn minkowski_support_point(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    dir: [f64; 3],
) -> SupportPoint {
    let a = shape_a.support(dir);
    let b = shape_b.support(neg3(dir));
    SupportPoint {
        w: sub3(a, b),
        a,
        b,
    }
}

// ─── Simplex ─────────────────────────────────────────────────────────────────

/// A simplex of up to 4 points in the Minkowski difference space.
#[derive(Clone)]
pub struct Simplex {
    /// Points of the simplex (at most 4), each with its support points on both
    /// shapes.
    pub points: Vec<SupportPoint>,
}

impl Simplex {
    fn new() -> Self {
        Self {
            points: Vec::with_capacity(4),
        }
    }

    fn push(&mut self, p: SupportPoint) {
        self.points.push(p);
    }

    fn len(&self) -> usize {
        self.points.len()
    }
}

// ─── do_simplex ──────────────────────────────────────────────────────────────

/// Process the simplex; returns true if origin is contained.
/// Also updates `direction` with the next search direction.
pub fn do_simplex(simplex: &mut Simplex, direction: &mut [f64; 3]) -> bool {
    match simplex.len() {
        2 => handle_line(simplex, direction),
        3 => handle_triangle(simplex, direction),
        4 => handle_tetrahedron(simplex, direction),
        _ => false,
    }
}

fn handle_line(simplex: &mut Simplex, direction: &mut [f64; 3]) -> bool {
    let b = simplex.points[0];
    let a = simplex.points[1];
    let ab = sub3(b.w, a.w);
    let ao = neg3(a.w);
    if dot3(ab, ao) > 0.0 {
        // direction = AB × (AO × AB) (triple product)
        let t = cross3(cross3(ab, ao), ab);
        if len3(t) < 1e-12 {
            // AO parallel to AB — pick perpendicular
            *direction = [1.0, 0.0, 0.0];
        } else {
            *direction = t;
        }
    } else {
        simplex.points = vec![a];
        *direction = ao;
    }
    false
}

fn handle_triangle(simplex: &mut Simplex, direction: &mut [f64; 3]) -> bool {
    let c = simplex.points[0];
    let b = simplex.points[1];
    let a = simplex.points[2];

    let ab = sub3(b.w, a.w);
    let ac = sub3(c.w, a.w);
    let ao = neg3(a.w);
    let abc = cross3(ab, ac);

    // Is origin outside edge AB (away from C)?
    if dot3(cross3(abc, ac), ao) > 0.0 {
        if dot3(ac, ao) > 0.0 {
            simplex.points = vec![c, a];
            *direction = cross3(cross3(ac, ao), ac);
        } else {
            simplex.points = vec![b, a];
            return handle_line(simplex, direction);
        }
    } else if dot3(cross3(ab, abc), ao) > 0.0 {
        simplex.points = vec![b, a];
        return handle_line(simplex, direction);
    } else {
        // Inside triangle; origin above or below?
        if dot3(abc, ao) > 0.0 {
            *direction = abc;
        } else {
            simplex.points = vec![b, c, a];
            *direction = neg3(abc);
        }
    }
    false
}

fn handle_tetrahedron(simplex: &mut Simplex, direction: &mut [f64; 3]) -> bool {
    let d = simplex.points[0];
    let c = simplex.points[1];
    let b = simplex.points[2];
    let a = simplex.points[3];

    let ab = sub3(b.w, a.w);
    let ac = sub3(c.w, a.w);
    let ad = sub3(d.w, a.w);
    let ao = neg3(a.w);

    let abc = cross3(ab, ac);
    let acd = cross3(ac, ad);
    let adb = cross3(ad, ab);

    if dot3(abc, ao) > 0.0 {
        simplex.points = vec![c, b, a];
        return handle_triangle(simplex, direction);
    }
    if dot3(acd, ao) > 0.0 {
        simplex.points = vec![d, c, a];
        return handle_triangle(simplex, direction);
    }
    if dot3(adb, ao) > 0.0 {
        simplex.points = vec![b, d, a];
        return handle_triangle(simplex, direction);
    }
    true
}

// ─── GJK intersection test ───────────────────────────────────────────────────

/// Run the GJK boolean iteration; returns the final simplex when the shapes
/// intersect and `None` when a separating direction was found.
fn gjk_simplex(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> Option<Simplex> {
    let mut direction = [1.0_f64, 0.0, 0.0];
    let mut simplex = Simplex::new();

    let first = minkowski_support_point(shape_a, shape_b, direction);
    simplex.push(first);
    direction = neg3(first.w);

    for _ in 0..64 {
        if len3(direction) < 1e-12 {
            return Some(simplex);
        }
        let next = minkowski_support_point(shape_a, shape_b, direction);
        if dot3(next.w, direction) < 0.0 {
            return None;
        }
        simplex.push(next);
        if do_simplex(&mut simplex, &mut direction) {
            return Some(simplex);
        }
    }
    None
}

/// Returns true if the two convex shapes intersect using GJK.
pub fn gjk_intersect(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> bool {
    gjk_simplex(shape_a, shape_b).is_some()
}

/// Euclidean distance between two convex shapes.
///
/// Returns `Some(distance)` when the shapes are disjoint and `None` when they
/// overlap or touch (distance below `1e-12` of the Minkowski-difference extent).
///
/// This runs the Gilbert-Johnson-Keerthi distance iteration: the closest point
/// `v` of the Minkowski difference `A - B` to the origin is refined with
/// support points until `|v|^2 - v.w <= 1e-12 |v|^2`, so the result is exact
/// (to that relative tolerance) rather than a support-vertex estimate. Use
/// [`gjk_distance_with_witnesses`] to also obtain the closest points.
pub fn gjk_distance(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> Option<f64> {
    match distance::solve(shape_a, shape_b) {
        distance::DistanceOutcome::Overlapping { .. } => None,
        distance::DistanceOutcome::Separated { distance, .. } => Some(distance),
    }
}

// ─── EPA ─────────────────────────────────────────────────────────────────────

/// Hard cap on the number of EPA expansion steps.
const EPA_MAX_ITERATIONS: usize = 128;

/// EPA stops once the support point along the normal of the nearest face lies
/// within this fraction of the polytope's extent beyond the face plane.
const EPA_TOLERANCE: f64 = 1e-9;

/// A support point must lie beyond a face plane by more than this fraction of
/// the polytope's extent for the face to count as visible from it.
const EPA_VISIBILITY: f64 = 1e-10;

/// A face of the polytope used in EPA.
#[derive(Debug, Clone)]
pub struct EpaFace {
    /// First vertex.
    pub a: SupportPoint,
    /// Second vertex.
    pub b: SupportPoint,
    /// Third vertex.
    pub c: SupportPoint,
    /// Outward-pointing unit normal (away from the interior of the polytope).
    pub normal: [f64; 3],
    /// Signed distance from the origin to the face plane along `normal`. It is
    /// negative while the origin is still outside the polytope.
    pub dist: f64,
}

impl EpaFace {
    /// Face through three vertices whose outward side is the one that faces away
    /// from `interior`, a point inside the polytope. The vertices are reordered
    /// so that they wind counter-clockwise seen from outside.
    ///
    /// Returns `None` for a degenerate (zero-area) triangle.
    fn oriented(
        a: SupportPoint,
        b: SupportPoint,
        c: SupportPoint,
        interior: [f64; 3],
    ) -> Option<Self> {
        let ab = sub3(b.w, a.w);
        let ac = sub3(c.w, a.w);
        let n = cross3(ab, ac);
        let length = len3(n);
        let reach = dot3(ab, ab).max(dot3(ac, ac));
        if !length.is_finite() || length <= 1e-13 * reach {
            return None;
        }
        let mut normal = scale3(n, 1.0 / length);
        let (mut b, mut c) = (b, c);
        if dot3(normal, sub3(interior, a.w)) > 0.0 {
            normal = neg3(normal);
            std::mem::swap(&mut b, &mut c);
        }
        Some(Self {
            a,
            b,
            c,
            normal,
            dist: dot3(normal, a.w),
        })
    }

    /// The three vertices of the face.
    fn corners(&self) -> [SupportPoint; 3] {
        [self.a, self.b, self.c]
    }
}

/// Result of the EPA penetration depth query.
pub struct EpaResult {
    /// Penetration depth (positive means overlap).
    pub depth: f64,
    /// Collision normal pointing from shape B to shape A: the direction in which
    /// shape A has to move to separate from shape B.
    pub normal: [f64; 3],
    /// Contact point on shape A: its deepest point inside shape B.
    pub contact_a: [f64; 3],
    /// Contact point on shape B: its deepest point inside shape A.
    ///
    /// `contact_a - contact_b == -normal * depth`.
    pub contact_b: [f64; 3],
}

/// Outcome of [`expand_polytope`]: the nearest facet of the Minkowski difference
/// `A - B`.
pub(crate) struct Penetration {
    /// Penetration depth: the distance of the nearest facet from the origin.
    pub(crate) depth: f64,
    /// Outward unit normal of that facet, from shape A toward shape B.
    pub(crate) normal: [f64; 3],
    /// Contact point on shape A: its deepest point inside shape B.
    pub(crate) contact_a: [f64; 3],
    /// Contact point on shape B: its deepest point inside shape A.
    pub(crate) contact_b: [f64; 3],
    /// Number of expansion steps performed.
    pub(crate) iterations: usize,
}

/// Turn the polytope's nearest boundary region into a [`Penetration`].
///
/// Several faces can lie in the nearest plane (the triangles of one facet of a
/// polyhedral Minkowski difference) and only one of them contains the foot of
/// the perpendicular from the origin. Among all faces whose plane is within a
/// relative `1e-9` of the nearest one, the face whose triangle is closest to
/// the origin is used. The contact points are the blends of the support points
/// of its vertices with the barycentric weights of the point of the triangle
/// nearest to the origin; that is the foot `n * dist` whenever the foot lies in
/// the triangle, which makes `contact_a - contact_b` equal `normal * depth`
/// exactly.
///
/// Returns `None` when `faces` is empty or the origin lies outside the
/// polytope by more than `tolerance`, which proves that the shapes are
/// disjoint.
fn penetration(faces: &[EpaFace], tolerance: f64, iterations: usize) -> Option<Penetration> {
    let nearest_plane = faces.iter().map(|f| f.dist).fold(f64::INFINITY, f64::min);
    if !nearest_plane.is_finite() || nearest_plane < -tolerance {
        return None;
    }
    let slack = 1e-9 * nearest_plane.abs().max(1.0);
    let mut best: Option<(f64, &EpaFace, simplex_closest::SimplexClosest)> = None;
    for face in faces.iter().filter(|f| f.dist <= nearest_plane + slack) {
        let corners = face.corners();
        let closest =
            simplex_closest::closest_point_to_origin(&[corners[0].w, corners[1].w, corners[2].w]);
        let reach = dot3(closest.point, closest.point);
        if best
            .as_ref()
            .is_none_or(|(smallest, _, _)| reach < *smallest)
        {
            best = Some((reach, face, closest));
        }
    }
    let (_, face, closest) = best?;
    let mut contact_a = [0.0; 3];
    let mut contact_b = [0.0; 3];
    for (corner, &weight) in face.corners().iter().zip(closest.weights.iter()) {
        contact_a = add3(contact_a, scale3(corner.a, weight));
        contact_b = add3(contact_b, scale3(corner.b, weight));
    }
    Some(Penetration {
        depth: face.dist.max(0.0),
        normal: face.normal,
        contact_a,
        contact_b,
        iterations,
    })
}

/// Key of an undirected edge: the bit patterns of its end points, ordered.
type EdgeKey = ([u64; 3], [u64; 3]);

fn edge_key(p: &SupportPoint, q: &SupportPoint) -> EdgeKey {
    let bits = |v: [f64; 3]| [v[0].to_bits(), v[1].to_bits(), v[2].to_bits()];
    let (kp, kq) = (bits(p.w), bits(q.w));
    if kp <= kq { (kp, kq) } else { (kq, kp) }
}

/// Extent of a set of Minkowski-difference points: the largest distance from
/// the origin, but at least one. Tolerances of the expansion are fractions of it.
pub(crate) fn seed_scale(points: &[SupportPoint]) -> f64 {
    points.iter().map(|p| len3(p.w)).fold(1.0_f64, f64::max)
}

/// Grow the polytope `seed` of the Minkowski difference `A - B` until its
/// nearest facet is a facet of `A - B`, and report that facet as a
/// [`Penetration`].
///
/// `support(d)` returns the point of `A - B` furthest along `d` together with
/// the support points on `A` and `B`. `seed` must hold four non-coplanar points
/// of `A - B`; it does not have to enclose the origin: faces are oriented
/// against the centroid of the seed, so the polytope first grows toward the
/// origin and then converges on the nearest facet.
///
/// Every step takes the face nearest to the origin (the most negative signed
/// distance first), asks for the support point along its normal, removes all
/// faces that can see that point, and joins the point to the horizon (the
/// edges of exactly one removed face). The iteration ends when the support point
/// lies within `tolerance` beyond the nearest face (never less than `1e-13`
/// of the seed extent, which is the rounding noise of that test) or after
/// `max_iterations` steps, in which case the nearest facet reached so far is
/// reported and its depth is a lower bound of the true depth.
///
/// Returns `None` when the seed is flat, or when the polytope cannot be grown
/// around the origin, which proves that the shapes do not overlap.
pub(crate) fn expand_polytope(
    seed: [SupportPoint; 4],
    support: &dyn Fn([f64; 3]) -> SupportPoint,
    tolerance: f64,
    max_iterations: usize,
) -> Option<Penetration> {
    let scale = seed_scale(&seed);
    let interior = scale3(
        add3(add3(seed[0].w, seed[1].w), add3(seed[2].w, seed[3].w)),
        0.25,
    );
    let volume = dot3(
        cross3(sub3(seed[1].w, seed[0].w), sub3(seed[2].w, seed[0].w)),
        sub3(seed[3].w, seed[0].w),
    )
    .abs();
    if volume.is_nan() || volume <= 1e-12 * scale * scale * scale {
        return None;
    }

    let mut faces: Vec<EpaFace> = Vec::with_capacity(32);
    for [i, j, k] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
        faces.push(EpaFace::oriented(seed[i], seed[j], seed[k], interior)?);
    }

    let tolerance = tolerance.max(1e-13 * scale);
    let margin = EPA_VISIBILITY * scale;
    let mut iterations = 0;
    while iterations < max_iterations {
        let nearest = faces
            .iter()
            .min_by(|x, y| x.dist.total_cmp(&y.dist))?
            .clone();
        let point = support(nearest.normal);
        if dot3(point.w, nearest.normal) - nearest.dist < tolerance {
            break;
        }
        iterations += 1;

        // Remove the faces that can see the new point; the horizon is made of the
        // edges that belong to exactly one of them.
        let (visible, kept): (Vec<EpaFace>, Vec<EpaFace>) = faces
            .into_iter()
            .partition(|f| dot3(f.normal, sub3(point.w, f.a.w)) > margin);
        faces = kept;
        if visible.is_empty() {
            // The point is beyond the nearest plane by less than the visibility
            // margin: the polytope cannot grow there any more.
            break;
        }
        let mut edges: BTreeMap<EdgeKey, (usize, SupportPoint, SupportPoint)> = BTreeMap::new();
        for face in &visible {
            for (p, q) in [(face.a, face.b), (face.b, face.c), (face.c, face.a)] {
                edges.entry(edge_key(&p, &q)).or_insert((0, p, q)).0 += 1;
            }
        }
        for (uses, p, q) in edges.into_values() {
            if uses == 1
                && let Some(face) = EpaFace::oriented(p, q, point, interior)
            {
                faces.push(face);
            }
        }
    }
    penetration(&faces, tolerance, iterations)
}

/// Expand polytope algorithm: finds penetration depth given initial simplex from GJK.
///
/// The simplex must hold four non-coplanar points of the Minkowski difference
/// `A - B` (a degenerate GJK simplex can be completed with the support points
/// along fixed directions, as [`gjk_epa_contact`] does). It does not have to
/// enclose the origin: faces are oriented against the centroid of the seed, so
/// the polytope first grows toward the origin and then converges on the
/// nearest facet.
///
/// Every step takes the face nearest to the origin, asks for the support point
/// along its normal, removes all faces that can see that point, and joins the
/// point to the horizon (the edges of exactly one removed face). The iteration
/// ends when the support point lies within a tolerance of the nearest face.
///
/// The returned normal points from shape B to shape A and the contact points
/// are the deepest points of each shape inside the other (see [`EpaResult`]).
///
/// Returns `None` when the simplex has fewer than four points, is flat, or
/// when the polytope cannot be grown around the origin, i.e. the shapes do not
/// overlap.
pub fn epa(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    simplex: Simplex,
) -> Option<EpaResult> {
    let seed = <[SupportPoint; 4]>::try_from(simplex.points.get(..4)?).ok()?;
    let oracle = |direction: [f64; 3]| minkowski_support_point(shape_a, shape_b, direction);
    let found = expand_polytope(
        seed,
        &oracle,
        EPA_TOLERANCE * seed_scale(&seed),
        EPA_MAX_ITERATIONS,
    )?;
    Some(EpaResult {
        depth: found.depth,
        normal: neg3(found.normal),
        contact_a: found.contact_a,
        contact_b: found.contact_b,
    })
}

// ─── Full GJK + EPA pipeline ─────────────────────────────────────────────────

/// Run GJK to detect intersection, then EPA to find penetration depth and contact.
/// Returns `None` if the shapes are not intersecting, or if their Minkowski
/// difference is flat (for example two zero-thickness shapes in one plane), which
/// leaves no volume to penetrate.
///
/// The contact points `contact_a` and `contact_b` are the deepest points of each
/// shape inside the other: blends of the support points that built the nearest
/// facet of the Minkowski difference, so `contact_a - contact_b` is the
/// penetration vector `-normal * depth`.
pub fn gjk_epa_contact(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> Option<EpaResult> {
    let simplex = gjk_simplex(shape_a, shape_b)?;
    let simplex = ensure_tetrahedron(shape_a, shape_b, simplex);
    epa(shape_a, shape_b, simplex)
}

/// Directions along which a degenerate GJK simplex is completed: the six axes
/// and the eight cube diagonals.
const COMPLETION_DIRECTIONS: [[f64; 3]; 14] = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
    [1.0, -1.0, -1.0],
    [-1.0, -1.0, -1.0],
];

/// Complete the GJK simplex to four non-coplanar points of the Minkowski
/// difference, for EPA.
///
/// GJK can stop with fewer than four points (the origin lies on a point,
/// segment or triangle of the simplex) or with repeated points. The distinct,
/// affinely independent points of `known` are kept and the remaining ones are
/// chosen among the support points along [`COMPLETION_DIRECTIONS`]: the
/// farthest from the first point, then the farthest from the line through the
/// two, then the one that rises the most above their plane (also trying the
/// support points along the plane normal). Every added point is a support point
/// of the Minkowski difference, so the result lies inside it. If the
/// Minkowski difference is itself flat there is no fourth point and the
/// result stays short, which [`expand_polytope`] callers report as no
/// penetration.
pub(crate) fn complete_tetrahedron(
    known: &[SupportPoint],
    support: &dyn Fn([f64; 3]) -> SupportPoint,
) -> Vec<SupportPoint> {
    let scale = seed_scale(known);
    let epsilon = 1e-9 * scale;

    // Affinely independent subset of the points GJK found.
    let mut basis: Vec<SupportPoint> = Vec::with_capacity(4);
    for &p in known {
        let independent = match basis.len() {
            0 => true,
            1 => len3(sub3(p.w, basis[0].w)) > epsilon,
            2 => {
                len3(cross3(sub3(basis[1].w, basis[0].w), sub3(p.w, basis[0].w))) > epsilon * scale
            }
            3 => {
                let normal = cross3(sub3(basis[1].w, basis[0].w), sub3(basis[2].w, basis[0].w));
                dot3(normal, sub3(p.w, basis[0].w)).abs() > epsilon * scale * scale
            }
            _ => false,
        };
        if independent {
            basis.push(p);
        }
    }
    if basis.len() == 4 {
        return basis;
    }

    let probes: Vec<SupportPoint> = COMPLETION_DIRECTIONS
        .iter()
        .map(|&d| support(normalize3(d)))
        .collect();
    let farthest = |measure: &dyn Fn(&SupportPoint) -> f64| {
        probes
            .iter()
            .copied()
            .map(|p| (measure(&p), p))
            .filter(|(m, _)| m.is_finite())
            .max_by(|x, y| x.0.total_cmp(&y.0))
    };

    if basis.is_empty()
        && let Some(&first) = probes.first()
    {
        basis.push(first);
    }
    if basis.len() == 1 {
        let origin = basis[0].w;
        match farthest(&|p| len3(sub3(p.w, origin))) {
            Some((reach, p)) if reach > epsilon => basis.push(p),
            _ => return basis,
        }
    }
    if basis.len() == 2 {
        let (origin, edge) = (basis[0].w, sub3(basis[1].w, basis[0].w));
        match farthest(&|p| len3(cross3(edge, sub3(p.w, origin)))) {
            Some((reach, p)) if reach > epsilon * scale => basis.push(p),
            _ => return basis,
        }
    }
    if basis.len() == 3 {
        let origin = basis[0].w;
        let normal = normalize3(cross3(
            sub3(basis[1].w, basis[0].w),
            sub3(basis[2].w, basis[0].w),
        ));
        let mut candidates = probes.clone();
        candidates.push(support(normal));
        candidates.push(support(neg3(normal)));
        let best = candidates
            .into_iter()
            .map(|p| (dot3(normal, sub3(p.w, origin)).abs(), p))
            .filter(|(m, _)| m.is_finite())
            .max_by(|x, y| x.0.total_cmp(&y.0));
        if let Some((rise, p)) = best
            && rise > epsilon
        {
            basis.push(p);
        }
    }
    basis
}

/// Complete a GJK simplex to a tetrahedron for [`epa`]; see
/// [`complete_tetrahedron`].
fn ensure_tetrahedron(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    simplex: Simplex,
) -> Simplex {
    let oracle = |direction: [f64; 3]| minkowski_support_point(shape_a, shape_b, direction);
    Simplex {
        points: complete_tetrahedron(&simplex.points, &oracle),
    }
}

// ─── GJK Distance query (improved) ───────────────────────────────────────────

/// Closest point on a line segment AB to the origin.
pub fn closest_point_segment_to_origin(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = sub3(b, a);
    let t = dot3(neg3(a), ab) / dot3(ab, ab).max(1e-30);
    let t = t.clamp(0.0, 1.0);
    add3(a, scale3(ab, t))
}

/// Closest point on triangle ABC to the origin (point projection).
///
/// Degenerate triangles (coincident or collinear vertices) are handled by
/// falling back to the nearest point of their edges instead of dividing by a
/// vanishing area.
pub fn closest_point_triangle_to_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    simplex_closest::closest_point_to_origin(&[a, b, c]).point
}

/// GJK distance query that also returns the closest points on both shapes.
///
/// Returns `Some((dist, witness_a, witness_b))` — the minimum distance and the
/// closest point on each shape — when the shapes do **not** intersect, and
/// `None` when they overlap or touch. The witnesses are the barycentric blend
/// of the support points of the final simplex, so
/// `witness_a - witness_b` is the closest point of the Minkowski difference to
/// the origin and `|witness_a - witness_b| == dist`. When the closest features
/// are parallel faces or edges the witness pair is one of infinitely many
/// equally close pairs. On smooth shapes (spheres, capsules) the distance is
/// accurate to about `1e-12` relative while the witnesses are accurate to
/// about its square root (`~1e-6` relative), because the distance is
/// quadratic in the witness position error.
pub fn gjk_distance_with_witnesses(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
) -> Option<(f64, [f64; 3], [f64; 3])> {
    match distance::solve(shape_a, shape_b) {
        distance::DistanceOutcome::Overlapping { .. } => None,
        distance::DistanceOutcome::Separated {
            distance,
            point_a,
            point_b,
            ..
        } => Some((distance, point_a, point_b)),
    }
}

// ─── GJK with rotation ───────────────────────────────────────────────────────

/// Apply a 3×3 rotation matrix (row-major) to a 3D vector.
pub fn rotate_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Rotate a convex shape by `rot` (row-major 3×3 rotation matrix) and
/// translate it by `translation`.  Returns a new support function wrapper.
pub struct TransformedShape<'a> {
    /// Underlying shape.
    pub shape: &'a dyn ConvexShape,
    /// Row-major 3×3 rotation matrix.
    pub rotation: [[f64; 3]; 3],
    /// Translation vector.
    pub translation: [f64; 3],
}

impl<'a> ConvexShape for TransformedShape<'a> {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        // Transform direction into local frame: R^T * dir
        let local_dir = [
            self.rotation[0][0] * direction[0]
                + self.rotation[1][0] * direction[1]
                + self.rotation[2][0] * direction[2],
            self.rotation[0][1] * direction[0]
                + self.rotation[1][1] * direction[1]
                + self.rotation[2][1] * direction[2],
            self.rotation[0][2] * direction[0]
                + self.rotation[1][2] * direction[1]
                + self.rotation[2][2] * direction[2],
        ];
        let local_pt = self.shape.support(local_dir);
        // Transform back to world frame and add translation
        add3(rotate_vec(self.rotation, local_pt), self.translation)
    }
}

/// Test intersection between two shapes given explicit poses.
pub fn gjk_intersect_transformed(
    shape_a: &dyn ConvexShape,
    rot_a: [[f64; 3]; 3],
    pos_a: [f64; 3],
    shape_b: &dyn ConvexShape,
    rot_b: [[f64; 3]; 3],
    pos_b: [f64; 3],
) -> bool {
    let ta = TransformedShape {
        shape: shape_a,
        rotation: rot_a,
        translation: pos_a,
    };
    let tb = TransformedShape {
        shape: shape_b,
        rotation: rot_b,
        translation: pos_b,
    };
    gjk_intersect(&ta, &tb)
}

// ─── EPA polytope expansion helpers ──────────────────────────────────────────

/// Count the number of faces in the EPA polytope after `n` expansion steps.
///
/// The initial tetrahedron has 4 faces.  Each expansion step removes at least
/// 1 face visible from the new support point and adds at least 3 new faces
/// along the silhouette edges.  This is a conservative *lower-bound* estimate.
pub fn epa_face_count_lower_bound(expansion_steps: usize) -> usize {
    4 + 2 * expansion_steps
}

/// A pair of closest points on two shapes — the result of a closest-point query.
#[derive(Debug, Clone)]
pub struct ClosestPoints {
    /// Closest point on shape A.
    pub point_a: [f64; 3],
    /// Closest point on shape B.
    pub point_b: [f64; 3],
    /// Distance between the two points.
    pub distance: f64,
}

/// Compute the closest points between two non-intersecting convex shapes.
///
/// Returns `None` if the shapes are intersecting.
pub fn closest_points_between_shapes(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
) -> Option<ClosestPoints> {
    gjk_distance_with_witnesses(shape_a, shape_b).map(|(dist, pa, pb)| ClosestPoints {
        point_a: pa,
        point_b: pb,
        distance: dist,
    })
}

// ─── GjkConvexHull ────────────────────────────────────────────────────────────

/// A convex hull defined by an explicit list of vertices.
pub struct GjkConvexHull {
    /// Vertices of the convex hull.
    pub vertices: Vec<[f64; 3]>,
}

impl GjkConvexHull {
    /// Construct from a list of points.
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        Self { vertices }
    }
}

impl ConvexShape for GjkConvexHull {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        let mut best = self.vertices[0];
        let mut best_dot = dot3(best, direction);
        for &v in &self.vertices[1..] {
            let d = dot3(v, direction);
            if d > best_dot {
                best_dot = d;
                best = v;
            }
        }
        best
    }
}

/// A finite right circular cylinder with its exact support function.
pub struct GjkCylinder {
    /// Center of the bottom cap.
    pub base: [f64; 3],
    /// Axis direction (unit vector).
    pub axis: [f64; 3],
    /// Height of the cylinder.
    pub height: f64,
    /// Radius of the cylinder.
    pub radius: f64,
}

/// A unit vector perpendicular to the unit vector `axis`.
///
/// The cross product with the coordinate axis least aligned with `axis` is well
/// conditioned: its length is at least `sqrt(2/3)`.
fn perpendicular_to(axis: [f64; 3]) -> [f64; 3] {
    let helper = if axis[0].abs() <= axis[1].abs() && axis[0].abs() <= axis[2].abs() {
        [1.0, 0.0, 0.0]
    } else if axis[1].abs() <= axis[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize3(cross3(axis, helper))
}

impl ConvexShape for GjkCylinder {
    fn support(&self, direction: [f64; 3]) -> [f64; 3] {
        // Project direction onto the axis and radial plane.
        let axis = normalize3(self.axis);
        let d_axial = dot3(direction, axis);
        // Radial component of direction
        let radial = sub3(direction, scale3(axis, d_axial));
        let radial_len = len3(radial);
        let radial_pt = if radial_len > 1e-12 {
            scale3(radial, self.radius / radial_len)
        } else {
            // The direction runs along the axis, so every rim point of the cap
            // supports it. A fixed offset such as `(radius, 0, 0)` is a rim
            // point only when the axis is perpendicular to x; take a point
            // perpendicular to the axis instead.
            scale3(perpendicular_to(axis), self.radius)
        };
        // Choose top or bottom cap
        let cap_offset = if d_axial >= 0.0 {
            scale3(axis, self.height)
        } else {
            [0.0; 3]
        };
        add3(add3(self.base, radial_pt), cap_offset)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn sphere(cx: f64, cy: f64, cz: f64, r: f64) -> GjkSphere {
        GjkSphere {
            center: [cx, cy, cz],
            radius: r,
        }
    }

    fn aabb(cx: f64, cy: f64, cz: f64, hx: f64, hy: f64, hz: f64) -> GjkBox {
        GjkBox {
            center: [cx, cy, cz],
            half_extents: [hx, hy, hz],
        }
    }

    fn capsule(p1: [f64; 3], p2: [f64; 3], r: f64) -> GjkCapsule {
        GjkCapsule { p1, p2, radius: r }
    }

    // ── math helpers ────────────────────────────────────────────────────────

    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_cross3() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3() {
        let n = normalize3([3.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3_zero() {
        let n = normalize3([0.0, 0.0, 0.0]);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_len3() {
        assert!((len3([3.0, 4.0, 0.0]) - 5.0).abs() < 1e-12);
    }

    // ── support functions ────────────────────────────────────────────────────

    #[test]
    fn test_sphere_support() {
        let s = sphere(0.0, 0.0, 0.0, 2.0);
        let p = s.support([1.0, 0.0, 0.0]);
        assert!((p[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_box_support() {
        let b = aabb(0.0, 0.0, 0.0, 1.0, 2.0, 3.0);
        let p = b.support([1.0, 1.0, 1.0]);
        assert_eq!(p, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_capsule_support() {
        let c = capsule([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let p = c.support([0.0, 1.0, 0.0]);
        assert!((p[1] - 1.5).abs() < 1e-10);
    }

    // ── sphere-sphere ────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_sphere_overlap() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(1.5, 0.0, 0.0, 1.0);
        assert!(gjk_intersect(&a, &b));
    }

    #[test]
    fn test_sphere_sphere_touching() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(2.0, 0.0, 0.0, 1.0);
        // Touching spheres: the boolean test is borderline and may report either,
        // but the distance query must agree with contact: no gap at all.
        let _ = gjk_intersect(&a, &b);
        match gjk_distance(&a, &b) {
            None => {}
            Some(d) => assert!(d < 1e-7, "tangent spheres reported a gap of {}", d),
        }
    }

    #[test]
    fn test_sphere_sphere_separated() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(3.0, 0.0, 0.0, 1.0);
        assert!(!gjk_intersect(&a, &b));
    }

    #[test]
    fn test_sphere_sphere_distance_separated() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(4.0, 0.0, 0.0, 1.0);
        // Surface gap: 4 - 1 - 1 = 2.
        let d = gjk_distance(&a, &b).expect("disjoint spheres have a distance");
        assert!((d - 2.0).abs() < 1e-6, "distance={}", d);
    }

    #[test]
    fn test_sphere_sphere_distance_intersecting() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(1.0, 0.0, 0.0, 1.0);
        assert!(gjk_distance(&a, &b).is_none());
    }

    // ── box-box ──────────────────────────────────────────────────────────────

    #[test]
    fn test_box_box_overlap() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(1.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(gjk_intersect(&a, &b));
    }

    #[test]
    fn test_box_box_separated() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(3.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(!gjk_intersect(&a, &b));
    }

    #[test]
    fn test_box_box_same_position() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(gjk_intersect(&a, &b));
    }

    #[test]
    fn test_box_box_distance() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        // Face gap: 4 - 1 - 1 = 2.
        let d = gjk_distance(&a, &b).expect("disjoint boxes have a distance");
        assert!((d - 2.0).abs() < 1e-6, "distance={}", d);
    }

    // ── sphere-box ───────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_box_overlap() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = aabb(1.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(gjk_intersect(&a, &b));
    }

    #[test]
    fn test_sphere_box_separated() {
        let a = sphere(0.0, 0.0, 0.0, 0.5);
        let b = aabb(3.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(!gjk_intersect(&a, &b));
    }

    #[test]
    fn test_sphere_box_distance() {
        let a = sphere(0.0, 0.0, 0.0, 0.5);
        let b = aabb(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        // Gap: box face at x = 3, sphere surface at x = 0.5.
        let d = gjk_distance(&a, &b).expect("disjoint shapes have a distance");
        assert!((d - 2.5).abs() < 1e-6, "distance={}", d);
    }

    // ── capsule tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_capsule_sphere_overlap() {
        let a = capsule([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        // capsule extends to x=1.0 (radius), sphere center at x=1.2 with r=0.5: overlap
        let b = sphere(1.2, 0.0, 0.0, 0.5);
        assert!(gjk_intersect(&a, &b));
    }

    #[test]
    fn test_capsule_sphere_separated() {
        let a = capsule([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let b = sphere(5.0, 0.0, 0.0, 0.5);
        assert!(!gjk_intersect(&a, &b));
    }

    #[test]
    fn test_capsule_box_overlap() {
        let a = capsule([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        let b = aabb(1.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(gjk_intersect(&a, &b));
    }

    // ── EPA contact ───────────────────────────────────────────────────────────

    #[test]
    fn test_epa_sphere_sphere_depth() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(1.0, 0.0, 0.0, 1.0);
        assert!(gjk_intersect(&a, &b), "spheres must intersect");
        // Radii 1 + 1 minus the centre distance 1.
        let r = gjk_epa_contact(&a, &b).expect("overlapping spheres have a contact");
        assert!((r.depth - 1.0).abs() < 1e-6, "depth {}", r.depth);
    }

    #[test]
    fn test_epa_sphere_sphere_no_contact_when_separated() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(3.0, 0.0, 0.0, 1.0);
        let result = gjk_epa_contact(&a, &b);
        assert!(result.is_none());
    }

    #[test]
    fn test_epa_box_box_depth() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(1.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let r = gjk_epa_contact(&a, &b).expect("overlapping boxes have a contact");
        // Half extents 1 + 1 minus the centre distance 1.5.
        assert!((r.depth - 0.5).abs() < 1e-9, "depth {}", r.depth);
    }

    #[test]
    fn test_epa_normal_unit_length() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(1.5, 0.0, 0.0, 1.0);
        let r = gjk_epa_contact(&a, &b).expect("overlapping spheres have a contact");
        let nl = len3(r.normal);
        assert!((nl - 1.0).abs() < 1e-6, "normal length={}", nl);
    }

    // ── minkowski_support ────────────────────────────────────────────────────

    #[test]
    fn test_minkowski_support_spheres() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(3.0, 0.0, 0.0, 1.0);
        // In direction [1,0,0]: sa = [1,0,0], sb_neg = sb.support([-1,0,0]) = [2,0,0]
        // minkowski = [1,0,0] - [2,0,0] = [-1,0,0]
        let m = minkowski_support(&a, &b, [1.0, 0.0, 0.0]);
        assert!((m[0] - (-1.0)).abs() < 1e-10, "m[0]={}", m[0]);
    }

    // ── closest_point_segment_to_origin ──────────────────────────────────────

    #[test]
    fn test_closest_point_segment_origin_inside_segment() {
        let p = closest_point_segment_to_origin([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(
            p[0].abs() < 1e-12,
            "closest point on segment to origin: {:?}",
            p
        );
    }

    #[test]
    fn test_closest_point_segment_origin_past_end() {
        // Origin is past B end → closest is B
        let p = closest_point_segment_to_origin([2.0, 0.0, 0.0], [3.0, 0.0, 0.0]);
        assert!((p[0] - 2.0).abs() < 1e-12, "expected [2,0,0], got {:?}", p);
    }

    // ── closest_point_triangle_to_origin ─────────────────────────────────────

    #[test]
    fn test_closest_point_triangle_origin_inside() {
        // Triangle in xz-plane containing origin: origin projects to itself
        let a = [-1.0, 0.0, -1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [0.0, 0.0, 1.0];
        let p = closest_point_triangle_to_origin(a, b, c);
        // origin is inside the triangle → closest point is origin itself
        assert!(
            len3(p) < 0.5,
            "closest point should be near origin, got {:?}",
            p
        );
    }

    #[test]
    fn test_closest_point_triangle_degenerate_inputs_stay_finite() {
        // Coincident vertices: the triangle is the segment a-c.
        let a = [1.0, 1.0, 0.0];
        let c = [-1.0, 1.0, 0.0];
        let p = closest_point_triangle_to_origin(a, a, c);
        assert!(len3(sub3(p, [0.0, 1.0, 0.0])) < 1e-12, "got {:?}", p);
        // Collinear vertices.
        let p =
            closest_point_triangle_to_origin([-3.0, 2.0, 0.0], [1.0, 2.0, 0.0], [5.0, 2.0, 0.0]);
        assert!(len3(sub3(p, [0.0, 2.0, 0.0])) < 1e-12, "got {:?}", p);
        // All three vertices equal.
        let p = closest_point_triangle_to_origin([0.0, 3.0, 4.0], [0.0, 3.0, 4.0], [0.0, 3.0, 4.0]);
        assert!(len3(sub3(p, [0.0, 3.0, 4.0])) < 1e-12, "got {:?}", p);
    }

    #[test]
    fn test_closest_point_triangle_interior_projection_is_exact() {
        // The origin projects into the interior of a triangle in the plane z = 3.
        let p =
            closest_point_triangle_to_origin([-1.0, -1.0, 3.0], [2.0, -1.0, 3.0], [-1.0, 2.0, 3.0]);
        assert!(len3(sub3(p, [0.0, 0.0, 3.0])) < 1e-12, "got {:?}", p);
    }

    #[test]
    fn test_closest_point_triangle_origin_outside_near_vertex() {
        // Triangle far from origin: closest should be vertex A
        let a = [10.0, 0.0, 0.0];
        let b = [11.0, 0.0, 0.0];
        let c = [10.0, 1.0, 0.0];
        let p = closest_point_triangle_to_origin(a, b, c);
        assert!((p[0] - 10.0).abs() < 1e-10, "expected x~10, got {:?}", p);
    }

    // ── gjk_distance_with_witnesses ───────────────────────────────────────────

    #[test]
    fn test_gjk_distance_witnesses_separated_spheres() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(5.0, 0.0, 0.0, 1.0);
        let result = gjk_distance_with_witnesses(&a, &b);
        assert!(result.is_some(), "separated spheres should have a distance");
        let (dist, pa, pb) = result.unwrap();
        assert!((dist - 3.0).abs() < 1e-6, "distance={}", dist);
        assert!((pa[0] - 1.0).abs() < 1e-5, "witness_a={:?}", pa);
        assert!((pb[0] - 4.0).abs() < 1e-5, "witness_b={:?}", pb);
        assert!((len3(sub3(pa, pb)) - dist).abs() < 1e-9);
    }

    #[test]
    fn test_gjk_distance_witnesses_intersecting_returns_none() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(0.5, 0.0, 0.0, 1.0);
        assert!(gjk_distance_with_witnesses(&a, &b).is_none());
    }

    // ── TransformedShape ─────────────────────────────────────────────────────

    fn identity_rot() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    #[test]
    fn test_transformed_shape_identity_matches_direct() {
        let s = sphere(0.0, 0.0, 0.0, 1.0);
        let t = TransformedShape {
            shape: &s,
            rotation: identity_rot(),
            translation: [0.0; 3],
        };
        let p_direct = s.support([1.0, 0.0, 0.0]);
        let p_trans = t.support([1.0, 0.0, 0.0]);
        assert!((p_direct[0] - p_trans[0]).abs() < 1e-10);
    }

    #[test]
    fn test_transformed_shape_translation() {
        let s = sphere(0.0, 0.0, 0.0, 1.0);
        let t = TransformedShape {
            shape: &s,
            rotation: identity_rot(),
            translation: [3.0, 0.0, 0.0],
        };
        let p = t.support([1.0, 0.0, 0.0]);
        // sphere center at (3,0,0) with r=1 → support in +x is (4,0,0)
        assert!((p[0] - 4.0).abs() < 1e-10, "p={:?}", p);
    }

    #[test]
    fn test_gjk_intersect_transformed_separated() {
        let s = sphere(0.0, 0.0, 0.0, 1.0);
        let intersects = gjk_intersect_transformed(
            &s,
            identity_rot(),
            [0.0; 3],
            &s,
            identity_rot(),
            [5.0, 0.0, 0.0],
        );
        assert!(!intersects);
    }

    #[test]
    fn test_gjk_intersect_transformed_overlapping() {
        let s = sphere(0.0, 0.0, 0.0, 1.0);
        let intersects = gjk_intersect_transformed(
            &s,
            identity_rot(),
            [0.0; 3],
            &s,
            identity_rot(),
            [1.5, 0.0, 0.0],
        );
        assert!(intersects);
    }

    // ── closest_points_between_shapes ────────────────────────────────────────

    #[test]
    fn test_closest_points_separated_boxes() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let cp = closest_points_between_shapes(&a, &b);
        assert!(cp.is_some(), "separated boxes should have closest points");
        let cp = cp.unwrap();
        assert!((cp.distance - 3.0).abs() < 1e-6, "distance={}", cp.distance);
        assert!(
            (cp.point_a[0] - 1.0).abs() < 1e-9,
            "point_a={:?}",
            cp.point_a
        );
        assert!(
            (cp.point_b[0] - 4.0).abs() < 1e-9,
            "point_b={:?}",
            cp.point_b
        );
    }

    #[test]
    fn test_closest_points_intersecting_returns_none() {
        let a = aabb(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = aabb(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        assert!(closest_points_between_shapes(&a, &b).is_none());
    }

    // ── GjkConvexHull ─────────────────────────────────────────────────────────

    #[test]
    fn test_convex_hull_support_in_x() {
        let hull = GjkConvexHull::new(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
        let p = hull.support([1.0, 0.0, 0.0]);
        assert!((p[0] - 2.0).abs() < 1e-12, "p={:?}", p);
    }

    #[test]
    fn test_convex_hull_intersects_sphere() {
        // Tetrahedron hull centered near origin
        let hull = GjkConvexHull::new(vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let s = sphere(0.0, 0.0, 0.0, 0.5);
        assert!(gjk_intersect(&hull, &s));
    }

    // ── GjkCylinder ──────────────────────────────────────────────────────────

    #[test]
    fn test_cylinder_support_top() {
        let cyl = GjkCylinder {
            base: [0.0, 0.0, 0.0],
            axis: [0.0, 1.0, 0.0],
            height: 2.0,
            radius: 0.5,
        };
        // Direction straight up: support should be near (0, 2, 0) + radial
        let p = cyl.support([0.0, 1.0, 0.0]);
        assert!(p[1] >= 1.9, "top support y should be near 2, got {}", p[1]);
    }

    #[test]
    fn test_cylinder_intersects_sphere_at_side() {
        let cyl = GjkCylinder {
            base: [0.0, 0.0, 0.0],
            axis: [0.0, 1.0, 0.0],
            height: 2.0,
            radius: 1.0,
        };
        let s = sphere(1.5, 1.0, 0.0, 1.0);
        assert!(gjk_intersect(&cyl, &s));
    }

    #[test]
    fn test_cylinder_separated_from_sphere() {
        let cyl = GjkCylinder {
            base: [0.0, 0.0, 0.0],
            axis: [0.0, 1.0, 0.0],
            height: 1.0,
            radius: 0.5,
        };
        let s = sphere(10.0, 0.5, 0.0, 0.5);
        assert!(!gjk_intersect(&cyl, &s));
    }

    #[test]
    fn test_cylinder_support_along_the_axis_is_a_point_of_the_cylinder() {
        // A direction parallel to the axis has no radial part: every rim point of
        // the cap supports it, but the point returned must belong to the cylinder.
        // The old fallback `(radius, 0, 0)` lies beyond the cap unless the axis is
        // perpendicular to x.
        let (height, radius) = (3.0, 0.75);
        for axis in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 2.0, -2.0],
            [-0.3, 0.1, 0.9],
        ] {
            let unit = normalize3(axis);
            let cyl = GjkCylinder {
                base: [0.5, -1.0, 2.0],
                axis: unit,
                height,
                radius,
            };
            for sign in [1.0, -1.0] {
                let p = cyl.support(scale3(unit, sign));
                let relative = sub3(p, cyl.base);
                let along = dot3(relative, unit);
                let rim = len3(sub3(relative, scale3(unit, along)));
                let cap = if sign > 0.0 { height } else { 0.0 };
                assert!(
                    (along - cap).abs() < 1e-12,
                    "axis {axis:?} sign {sign}: axial position {along}, cap at {cap}"
                );
                assert!(
                    (rim - radius).abs() < 1e-12,
                    "axis {axis:?} sign {sign}: distance from the axis {rim}"
                );
            }
        }
    }

    #[test]
    fn test_cylinder_along_x_is_not_reported_in_contact_beyond_its_end_cap() {
        // The cylinder fills x in [0, 2] around the x axis. GJK starts with the
        // direction +x, which is parallel to the axis: the support point beyond the
        // cap made a sphere 1.0 away look like it touched the cylinder.
        let cyl = GjkCylinder {
            base: [0.0; 3],
            axis: [1.0, 0.0, 0.0],
            height: 2.0,
            radius: 1.0,
        };
        let apart = sphere(3.5, 0.0, 0.0, 0.5);
        assert!(!gjk_intersect(&cyl, &apart));
        assert!(gjk_epa_contact(&cyl, &apart).is_none());
        let d = gjk_distance(&cyl, &apart).expect("the shapes are apart");
        assert!((d - 1.0).abs() < 1e-6, "distance {d}");
        let (dist, on_cylinder, on_sphere) =
            gjk_distance_with_witnesses(&cyl, &apart).expect("the shapes are apart");
        assert!((dist - 1.0).abs() < 1e-6);
        assert!((on_cylinder[0] - 2.0).abs() < 1e-6, "{on_cylinder:?}");
        assert!((on_sphere[0] - 3.0).abs() < 1e-6, "{on_sphere:?}");

        // Overlapping the end cap by 0.2 gives a depth of exactly 0.2 along -x.
        let overlapping = sphere(2.3, 0.0, 0.0, 0.5);
        assert!(gjk_intersect(&cyl, &overlapping));
        let contact = gjk_epa_contact(&cyl, &overlapping).expect("the shapes overlap");
        assert!(
            (contact.depth - 0.2).abs() < 1e-6,
            "depth {}",
            contact.depth
        );
        assert!(
            (contact.normal[0] + 1.0).abs() < 1e-6,
            "{:?}",
            contact.normal
        );
    }

    // ── epa_face_count_lower_bound ────────────────────────────────────────────

    #[test]
    fn test_epa_face_count_lower_bound() {
        assert_eq!(epa_face_count_lower_bound(0), 4);
        assert_eq!(epa_face_count_lower_bound(3), 10);
    }

    // ── rotate_vec ────────────────────────────────────────────────────────────

    #[test]
    fn test_rotate_vec_identity() {
        let v = [1.0, 2.0, 3.0];
        let r = rotate_vec([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], v);
        assert!((r[0] - v[0]).abs() < 1e-12);
        assert!((r[1] - v[1]).abs() < 1e-12);
        assert!((r[2] - v[2]).abs() < 1e-12);
    }

    #[test]
    fn test_rotate_vec_90_around_z() {
        // Rotation 90° around Z: (1,0,0) → (0,1,0)
        let rot = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let r = rotate_vec(rot, [1.0, 0.0, 0.0]);
        assert!((r[0]).abs() < 1e-12, "r[0]={}", r[0]);
        assert!((r[1] - 1.0).abs() < 1e-12, "r[1]={}", r[1]);
    }
}
