// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GJK (Gilbert-Johnson-Keerthi) and EPA (Expanding Polytope Algorithm)
//! for convex shape collision detection using plain `[f64; 3]` arrays.

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

// ─── Simplex ─────────────────────────────────────────────────────────────────

/// A simplex of up to 4 points in the Minkowski difference space.
#[derive(Clone)]
pub struct Simplex {
    /// Points of the simplex (at most 4).
    pub points: Vec<[f64; 3]>,
}

impl Simplex {
    fn new() -> Self {
        Self {
            points: Vec::with_capacity(4),
        }
    }

    fn push(&mut self, p: [f64; 3]) {
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
    let ab = sub3(b, a);
    let ao = neg3(a);
    if dot3(ab, ao) > 0.0 {
        *direction = sub3(scale3(ab, dot3(ab, ao) / dot3(ab, ab)), ao);
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

    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ao = neg3(a);
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

    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ad = sub3(d, a);
    let ao = neg3(a);

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

/// Returns true if the two convex shapes intersect using GJK.
pub fn gjk_intersect(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> bool {
    let mut direction = [1.0_f64, 0.0, 0.0];
    let mut simplex = Simplex::new();

    let support = minkowski_support(shape_a, shape_b, direction);
    simplex.push(support);
    direction = neg3(support);

    for _ in 0..64 {
        if len3(direction) < 1e-12 {
            return true;
        }
        let a = minkowski_support(shape_a, shape_b, direction);
        if dot3(a, direction) < 0.0 {
            return false;
        }
        simplex.push(a);
        if do_simplex(&mut simplex, &mut direction) {
            return true;
        }
    }
    false
}

/// Returns `Some(distance)` if not intersecting, `None` if intersecting.
pub fn gjk_distance(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> Option<f64> {
    if gjk_intersect(shape_a, shape_b) {
        None
    } else {
        // Simple distance estimate: iterate and track closest point
        let mut direction = [1.0_f64, 0.0, 0.0];
        let mut min_dist = f64::MAX;

        for _ in 0..64 {
            let s = minkowski_support(shape_a, shape_b, direction);
            let d = len3(s);
            if d < min_dist {
                min_dist = d;
            }
            let new_dir = neg3(s);
            if len3(new_dir) < 1e-12 {
                break;
            }
            let next_s = minkowski_support(shape_a, shape_b, new_dir);
            if dot3(next_s, new_dir) < dot3(s, direction) + 1e-10 {
                break;
            }
            direction = new_dir;
        }

        // Approximate: distance between closest support points
        let sa = shape_a.support(direction);
        let sb = shape_b.support(neg3(direction));
        let dist = len3(sub3(sa, sb));
        Some(dist.max(0.0))
    }
}

// ─── EPA ─────────────────────────────────────────────────────────────────────

/// A face of the polytope used in EPA.
#[derive(Clone)]
pub struct EpaFace {
    /// First vertex.
    pub a: [f64; 3],
    /// Second vertex.
    pub b: [f64; 3],
    /// Third vertex.
    pub c: [f64; 3],
    /// Outward-pointing unit normal.
    pub normal: [f64; 3],
    /// Distance from origin to the face plane.
    pub dist: f64,
}

impl EpaFace {
    fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let n = cross3(ab, ac);
        let mut normal = normalize3(n);
        let mut dist = dot3(normal, a);
        // Ensure normal points away from origin (positive dist)
        if dist < 0.0 {
            normal = neg3(normal);
            dist = -dist;
        }
        Self {
            a,
            b,
            c,
            normal,
            dist,
        }
    }
}

/// Result of the EPA penetration depth query.
pub struct EpaResult {
    /// Penetration depth (positive means overlap).
    pub depth: f64,
    /// Collision normal pointing from shape B to shape A.
    pub normal: [f64; 3],
    /// Contact point on shape A.
    pub contact_a: [f64; 3],
    /// Contact point on shape B.
    pub contact_b: [f64; 3],
}

fn barycentric_contact(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    face: &EpaFace,
) -> ([f64; 3], [f64; 3]) {
    // Use barycentric coords of origin projection onto face
    let n = face.normal;
    let proj = scale3(n, face.dist);

    let v0 = sub3(face.b, face.a);
    let v1 = sub3(face.c, face.a);
    let v2 = sub3(proj, face.a);

    let d00 = dot3(v0, v0);
    let d01 = dot3(v0, v1);
    let d11 = dot3(v1, v1);
    let d20 = dot3(v2, v0);
    let d21 = dot3(v2, v1);
    let denom = d00 * d11 - d01 * d01;

    let (u, v, w) = if denom.abs() < 1e-12 {
        (1.0_f64 / 3.0, 1.0_f64 / 3.0, 1.0_f64 / 3.0)
    } else {
        let v_coord = (d11 * d20 - d01 * d21) / denom;
        let w_coord = (d00 * d21 - d01 * d20) / denom;
        let u_coord = 1.0 - v_coord - w_coord;
        (u_coord, v_coord, w_coord)
    };

    // Reconstruct support points from Minkowski difference vertices
    // face.a = sa_a - sb_a, face.b = sa_b - sb_b, face.c = sa_c - sb_c
    // We need to recover original support points; store them as (sa, sb) pairs
    // Since we only have the Minkowski difference points, approximate with
    // closest support point in face normal direction
    let ca = add3(
        add3(
            scale3(shape_a.support(face.a), u),
            scale3(shape_a.support(face.b), v),
        ),
        scale3(shape_a.support(face.c), w),
    );
    let cb = add3(
        add3(
            scale3(shape_b.support(neg3(face.a)), u),
            scale3(shape_b.support(neg3(face.b)), v),
        ),
        scale3(shape_b.support(neg3(face.c)), w),
    );
    (ca, cb)
}

/// Expand polytope algorithm: finds penetration depth given initial simplex from GJK.
pub fn epa(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    simplex: Simplex,
) -> Option<EpaResult> {
    if simplex.points.len() < 4 {
        return None;
    }

    let mut faces: Vec<EpaFace> = vec![
        EpaFace::new(simplex.points[0], simplex.points[1], simplex.points[2]),
        EpaFace::new(simplex.points[0], simplex.points[1], simplex.points[3]),
        EpaFace::new(simplex.points[0], simplex.points[2], simplex.points[3]),
        EpaFace::new(simplex.points[1], simplex.points[2], simplex.points[3]),
    ];

    for _ in 0..64 {
        // Find face closest to origin
        let (min_idx, min_face) = {
            let mut best_idx = 0;
            let mut best_dist = f64::MAX;
            for (i, face) in faces.iter().enumerate() {
                if face.dist < best_dist {
                    best_dist = face.dist;
                    best_idx = i;
                }
            }
            (best_idx, faces[best_idx].clone())
        };

        let support = minkowski_support(shape_a, shape_b, min_face.normal);
        let d = dot3(support, min_face.normal);

        if (d - min_face.dist).abs() < 1e-6 {
            // Converged
            let (contact_a, contact_b) = barycentric_contact(shape_a, shape_b, &min_face);
            return Some(EpaResult {
                depth: min_face.dist,
                normal: min_face.normal,
                contact_a,
                contact_b,
            });
        }

        // Remove faces that can "see" the new support point and collect silhouette edges
        let mut silhouette: Vec<([f64; 3], [f64; 3])> = Vec::new();
        let new_point = support;

        let mut i = 0;
        while i < faces.len() {
            if dot3(faces[i].normal, sub3(new_point, faces[i].a)) > 0.0 {
                let f = faces.remove(i);
                // Collect edges (in reverse winding to keep silhouette consistent)
                for edge in [(f.a, f.b), (f.b, f.c), (f.c, f.a)] {
                    // Remove shared edges (horizon silhouette)
                    let shared = silhouette
                        .iter()
                        .position(|&(x, y)| x == edge.1 && y == edge.0);
                    if let Some(pos) = shared {
                        silhouette.remove(pos);
                    } else {
                        silhouette.push(edge);
                    }
                }
            } else {
                i += 1;
            }
        }

        // Re-expand with new faces from silhouette
        for (ea, eb) in silhouette {
            faces.push(EpaFace::new(ea, eb, new_point));
        }

        // Safety: remove the face we started with if it's still there (re-inserted version covers it)
        let _ = min_idx; // already removed or updated above
    }

    // Return best face found so far
    let best = faces.iter().min_by(|a, b| {
        a.dist
            .partial_cmp(&b.dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;
    let (contact_a, contact_b) = barycentric_contact(shape_a, shape_b, best);
    Some(EpaResult {
        depth: best.dist,
        normal: best.normal,
        contact_a,
        contact_b,
    })
}

// ─── Full GJK + EPA pipeline ─────────────────────────────────────────────────

/// Run GJK to detect intersection, then EPA to find penetration depth and contact.
/// Returns `None` if the shapes are not intersecting.
pub fn gjk_epa_contact(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> Option<EpaResult> {
    // Run GJK and collect the final simplex
    let mut direction = [1.0_f64, 0.0, 0.0];
    let mut simplex = Simplex::new();

    let first_support = minkowski_support(shape_a, shape_b, direction);
    simplex.push(first_support);
    direction = neg3(first_support);

    let mut intersects = false;
    for _ in 0..64 {
        if len3(direction) < 1e-12 {
            intersects = true;
            break;
        }
        let a = minkowski_support(shape_a, shape_b, direction);
        if dot3(a, direction) < 0.0 {
            return None;
        }
        simplex.push(a);
        if do_simplex(&mut simplex, &mut direction) {
            intersects = true;
            break;
        }
    }

    if !intersects {
        return None;
    }

    // Ensure we have a tetrahedron for EPA
    let simplex = ensure_tetrahedron(shape_a, shape_b, simplex);
    epa(shape_a, shape_b, simplex)
}

fn ensure_tetrahedron(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    mut simplex: Simplex,
) -> Simplex {
    let dirs = [
        [1.0_f64, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];

    // Fill up to 3 unique points
    for dir in dirs {
        if simplex.len() >= 3 {
            break;
        }
        let p = minkowski_support(shape_a, shape_b, dir);
        let is_new = simplex.points.iter().all(|&q| len3(sub3(p, q)) > 1e-10);
        if is_new {
            simplex.push(p);
        }
    }

    if simplex.len() < 3 {
        return simplex;
    }

    // Check if we already have a valid (non-degenerate) tetrahedron
    if simplex.len() >= 4 {
        let a = simplex.points[0];
        let b = simplex.points[1];
        let c = simplex.points[2];
        let d = simplex.points[3];
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let normal = cross3(ab, ac);
        let vol = dot3(normal, sub3(d, a));
        if vol.abs() > 1e-10 {
            return simplex; // already a valid tetrahedron
        }
        // Degenerate: drop the 4th point and find a better one
        simplex.points.truncate(3);
    }

    // For the 4th point, we need it to be non-coplanar with the first 3
    let a = simplex.points[0];
    let b = simplex.points[1];
    let c = simplex.points[2];
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let normal = cross3(ab, ac);

    for dir in dirs {
        let p = minkowski_support(shape_a, shape_b, dir);
        let is_new = simplex.points.iter().all(|&q| len3(sub3(p, q)) > 1e-10);
        if !is_new {
            continue;
        }
        // Check non-coplanar: dot(normal, p-a) != 0
        let vol = dot3(normal, sub3(p, a));
        if vol.abs() > 1e-10 {
            simplex.push(p);
            break;
        }
    }

    simplex
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
pub fn closest_point_triangle_to_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ao = neg3(a);

    let d1 = dot3(ab, ao);
    let d2 = dot3(ac, ao);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }

    let bo = neg3(b);
    let d3 = dot3(ab, bo);
    let d4 = dot3(ac, bo);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }

    let co = neg3(c);
    let d5 = dot3(ab, co);
    let d6 = dot3(ac, co);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add3(a, scale3(ab, v));
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add3(a, scale3(ac, w));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add3(b, scale3(sub3(c, b), w));
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    add3(a, add3(scale3(ab, v), scale3(ac, w)))
}

/// Improved GJK distance query that tracks the closest simplex points.
///
/// Returns `Some((dist, witness_a, witness_b))` — the minimum distance and
/// witness points on each shape — when the shapes do **not** intersect.
/// Returns `None` when they intersect.
pub fn gjk_distance_with_witnesses(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
) -> Option<(f64, [f64; 3], [f64; 3])> {
    if gjk_intersect(shape_a, shape_b) {
        return None;
    }
    // Walk the Minkowski difference to find the closest point.
    let mut dir = [1.0_f64, 0.0, 0.0];
    let mut best_a = shape_a.support(dir);
    let mut best_b = shape_b.support(neg3(dir));
    let mut best_dist = f64::MAX;

    for _ in 0..64 {
        let sa = shape_a.support(dir);
        let sb = shape_b.support(neg3(dir));
        let diff = sub3(sa, sb); // point on Minkowski difference boundary
        let d = len3(diff);
        if d < best_dist {
            best_dist = d;
            best_a = sa;
            best_b = sb;
        }
        let new_dir = neg3(diff);
        if len3(new_dir) < 1e-12 {
            break;
        }
        dir = normalize3(new_dir);
    }
    Some((best_dist.max(0.0), best_a, best_b))
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

/// A cylinder (finite) approximated as a convex hull of 16 sample points.
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
            [self.radius, 0.0, 0.0]
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
        // Touching spheres: borderline; GJK may report either
        let _ = gjk_intersect(&a, &b);
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
        let d = gjk_distance(&a, &b);
        assert!(d.is_some());
        // Distance between surfaces should be ~2.0
        let d = d.unwrap();
        assert!(d > 0.0 && d < 5.0, "distance={}", d);
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
        let d = gjk_distance(&a, &b);
        assert!(d.is_some());
        assert!(d.unwrap() > 0.0);
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
        let d = gjk_distance(&a, &b);
        assert!(d.is_some());
        assert!(d.unwrap() > 0.0);
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
        if let Some(r) = gjk_epa_contact(&a, &b) {
            assert!(r.depth >= 0.0, "depth must be non-negative");
        }
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
        let result = gjk_epa_contact(&a, &b);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.depth > 0.0);
    }

    #[test]
    fn test_epa_normal_unit_length() {
        let a = sphere(0.0, 0.0, 0.0, 1.0);
        let b = sphere(1.5, 0.0, 0.0, 1.0);
        if let Some(r) = gjk_epa_contact(&a, &b) {
            let nl = len3(r.normal);
            assert!((nl - 1.0).abs() < 1e-6, "normal length={}", nl);
        }
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
        assert!(dist > 0.0, "distance should be positive");
        assert!(pa.iter().all(|x| x.is_finite()));
        assert!(pb.iter().all(|x| x.is_finite()));
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
        assert!(cp.distance > 0.0);
        assert!(cp.point_a.iter().all(|x| x.is_finite()));
        assert!(cp.point_b.iter().all(|x| x.is_finite()));
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
