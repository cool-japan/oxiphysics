// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced contact manifold generation for collision detection.
//!
//! Provides:
//! - Polyhedron face-face, face-edge, edge-edge contact
//! - Box-box specialized manifold
//! - Capsule-capsule contact
//! - Mesh vs convex contact
//! - 4-point manifold reduction
//! - Persistent manifold with frame-to-frame tracking
//! - Feature contact (vertex/edge/face type tracking)
//! - Contact normal smoothing on mesh geometry
//! - Contact caching with LRU
//! - Contact event system (begin/end events)

mod box_box;
pub use box_box::BoxBoxManifold;

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared norm of a 3-vector.
#[inline]
pub fn norm_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

/// Norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    norm_sq3(a).sqrt()
}

/// Normalize a 3-vector.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-15 {
        return [0.0, 1.0, 0.0];
    }
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Negate a 3-vector.
#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ---------------------------------------------------------------------------
// Contact point
// ---------------------------------------------------------------------------

/// A single contact point in a manifold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactPoint {
    /// World-space contact position on body A.
    pub pos_a: [f64; 3],
    /// World-space contact position on body B.
    pub pos_b: [f64; 3],
    /// Contact normal pointing from B to A.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Feature ID on body A.
    pub feature_a: u32,
    /// Feature ID on body B.
    pub feature_b: u32,
    /// Warm-start accumulated normal impulse.
    pub impulse: f64,
}

impl ContactPoint {
    /// Create a new contact point.
    pub fn new(pos_a: [f64; 3], pos_b: [f64; 3], normal: [f64; 3], depth: f64) -> Self {
        Self {
            pos_a,
            pos_b,
            normal,
            depth,
            feature_a: 0,
            feature_b: 0,
            impulse: 0.0,
        }
    }

    /// Contact point midpoint.
    pub fn midpoint(&self) -> [f64; 3] {
        scale3(add3(self.pos_a, self.pos_b), 0.5)
    }
}

// ---------------------------------------------------------------------------
// Contact manifold
// ---------------------------------------------------------------------------

/// Maximum contact points in a manifold.
pub const MAX_MANIFOLD_POINTS: usize = 4;

/// A contact manifold between two bodies.
#[derive(Debug, Clone)]
pub struct ContactManifoldEx {
    /// Contact points (up to MAX_MANIFOLD_POINTS).
    pub points: Vec<ContactPoint>,
    /// Shared contact normal.
    pub normal: [f64; 3],
    /// Body A ID.
    pub body_a: u32,
    /// Body B ID.
    pub body_b: u32,
    /// Frame counter when last updated.
    pub last_frame: u64,
}

impl ContactManifoldEx {
    /// Create an empty contact manifold.
    pub fn new(body_a: u32, body_b: u32) -> Self {
        Self {
            points: Vec::new(),
            normal: [0.0, 1.0, 0.0],
            body_a,
            body_b,
            last_frame: 0,
        }
    }

    /// Number of contact points.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Add a contact point (up to MAX_MANIFOLD_POINTS).
    pub fn add_point(&mut self, cp: ContactPoint) {
        if self.points.len() < MAX_MANIFOLD_POINTS {
            self.points.push(cp);
        }
    }

    /// Maximum penetration depth.
    pub fn max_depth(&self) -> f64 {
        self.points.iter().map(|p| p.depth).fold(0.0_f64, f64::max)
    }

    /// Is this manifold active (has contact points with positive depth)?
    pub fn is_active(&self) -> bool {
        self.points.iter().any(|p| p.depth > 0.0)
    }

    /// Remove all points with non-positive depth.
    pub fn prune_inactive(&mut self) {
        self.points.retain(|p| p.depth > -0.01);
    }
}

// ---------------------------------------------------------------------------
// Clip polygon by plane
// ---------------------------------------------------------------------------

/// Clip a polygon against a half-plane defined by a support point and normal.
///
/// Points on the positive side of the plane (dot(p - s, n) >= 0) are kept.
/// Returns clipped polygon vertices.
pub fn clip_polygon_by_plane(
    polygon: &[[f64; 3]],
    plane_point: [f64; 3],
    plane_normal: [f64; 3],
) -> Vec<[f64; 3]> {
    if polygon.is_empty() {
        return Vec::new();
    }
    let mut result = Vec::new();
    let n = polygon.len();

    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let da = dot3(sub3(a, plane_point), plane_normal);
        let db = dot3(sub3(b, plane_point), plane_normal);

        if da >= 0.0 {
            result.push(a);
        }
        if (da >= 0.0) != (db >= 0.0) {
            // Edge crosses plane, compute intersection
            let t = da / (da - db);
            let p = add3(a, scale3(sub3(b, a), t));
            result.push(p);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// SAT helpers
// ---------------------------------------------------------------------------

/// SAT face axis test between two convex bodies.
///
/// Tests separation along a face normal axis.
/// Returns the separation distance (negative = overlap).
pub fn sat_face_axis(axis: [f64; 3], ref_pts: &[[f64; 3]], inc_pts: &[[f64; 3]]) -> f64 {
    let ref_max = ref_pts
        .iter()
        .map(|&p| dot3(p, axis))
        .fold(f64::NEG_INFINITY, f64::max);
    let inc_min = inc_pts
        .iter()
        .map(|&p| dot3(p, axis))
        .fold(f64::INFINITY, f64::min);
    inc_min - ref_max
}

/// SAT edge axis test between two edges.
///
/// Tests separation along the cross product of two edge directions.
pub fn sat_edge_axis(
    edge_a: [f64; 3],
    edge_b: [f64; 3],
    ref_center: [f64; 3],
    ref_pts: &[[f64; 3]],
    inc_pts: &[[f64; 3]],
) -> f64 {
    let axis = cross3(edge_a, edge_b);
    let len = norm3(axis);
    if len < 1e-10 {
        return f64::NEG_INFINITY;
    } // Parallel edges
    let n = scale3(axis, 1.0 / len);

    // Orient normal toward ref_center
    let n = if dot3(sub3(ref_center, inc_pts[0]), n) < 0.0 {
        neg3(n)
    } else {
        n
    };

    let ref_max = ref_pts
        .iter()
        .map(|&p| dot3(p, n))
        .fold(f64::NEG_INFINITY, f64::max);
    let inc_min = inc_pts
        .iter()
        .map(|&p| dot3(p, n))
        .fold(f64::INFINITY, f64::min);
    inc_min - ref_max
}

/// Find the incident face of a convex body most anti-parallel to a reference normal.
///
/// Returns the face index and its normal.
pub fn find_incident_face(normals: &[[f64; 3]], ref_normal: [f64; 3]) -> (usize, [f64; 3]) {
    let mut best_idx = 0;
    let mut best_dot = f64::INFINITY;
    for (i, &n) in normals.iter().enumerate() {
        let d = dot3(n, ref_normal);
        if d < best_dot {
            best_dot = d;
            best_idx = i;
        }
    }
    (best_idx, normals[best_idx])
}

// ---------------------------------------------------------------------------
// PolyhedronContact
// ---------------------------------------------------------------------------

/// Face-face, face-edge, and edge-edge contact generation for convex polyhedra.
///
/// Uses the Sutherland-Hodgman clipping algorithm for face-face contacts.
#[derive(Debug, Clone)]
pub struct PolyhedronContact {
    /// Maximum contact points to generate.
    pub max_contacts: usize,
}

impl PolyhedronContact {
    /// Create a new polyhedron contact generator.
    pub fn new(max_contacts: usize) -> Self {
        Self { max_contacts }
    }

    /// Compute face-face contact between two convex polyhedra.
    ///
    /// # Arguments
    /// * `ref_face` – polygon of reference face vertices
    /// * `inc_face` – polygon of incident face vertices
    /// * `ref_normal` – reference face outward normal
    /// * `depth` – maximum clipping depth threshold
    pub fn face_face_contact(
        &self,
        ref_face: &[[f64; 3]],
        inc_face: &[[f64; 3]],
        ref_normal: [f64; 3],
        depth: f64,
    ) -> Vec<ContactPoint> {
        if ref_face.is_empty() || inc_face.is_empty() {
            return Vec::new();
        }

        // Clip incident face against reference face side planes
        let mut clipped = inc_face.to_vec();
        let n = ref_face.len();
        for i in 0..n {
            if clipped.is_empty() {
                break;
            }
            let edge_start = ref_face[i];
            let edge_end = ref_face[(i + 1) % n];
            let edge_dir = sub3(edge_end, edge_start);
            // Side plane normal (pointing inward)
            let side_normal = normalize3(cross3(ref_normal, edge_dir));
            clipped = clip_polygon_by_plane(&clipped, edge_start, side_normal);
        }

        // Generate contact points from clipped polygon
        let mut contacts = Vec::new();
        let ref_plane_d = dot3(ref_face[0], ref_normal);
        for p in &clipped {
            let pen = ref_plane_d - dot3(*p, ref_normal);
            if pen > -depth {
                let pos_b = *p;
                let pos_a = add3(*p, scale3(ref_normal, pen));
                contacts.push(ContactPoint::new(pos_a, pos_b, ref_normal, pen));
                if contacts.len() >= self.max_contacts {
                    break;
                }
            }
        }
        contacts
    }

    /// Compute edge-edge contact distance and closest points.
    ///
    /// Returns (closest_point_a, closest_point_b, distance).
    pub fn edge_edge_closest(
        &self,
        pa: [f64; 3],
        da: [f64; 3],
        pb: [f64; 3],
        db: [f64; 3],
    ) -> ([f64; 3], [f64; 3], f64) {
        // Squared edge length (and denominator) below which an edge counts as
        // a point (or the two edges as parallel).
        const DEGENERATE: f64 = 1e-14;

        let r = sub3(pa, pb);
        let a = dot3(da, da);
        let e = dot3(db, db);
        let f = dot3(db, r);

        if a < DEGENERATE && e < DEGENERATE {
            return (pa, pb, norm3(r));
        }

        // Parameters of the closest points, `s` along edge A and `t` along
        // edge B, both in [0, 1] (Ericson, Real-Time Collision Detection 5.1.9).
        let (s, t) = if a < DEGENERATE {
            (0.0, (f / e).clamp(0.0, 1.0))
        } else {
            let c = dot3(da, r);
            if e < DEGENERATE {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else {
                let b = dot3(da, db);
                let denom = a * e - b * b;
                // Closest points of the two infinite lines, clamped to edge A.
                // Parallel edges have no unique pair: take the middle of A.
                let s_line = if denom.abs() > DEGENERATE {
                    ((b * f - c * e) / denom).clamp(0.0, 1.0)
                } else {
                    0.5
                };
                let t_line = (b * s_line + f) / e;
                if t_line < 0.0 {
                    // Beyond the start of B: pin t and re-solve s.
                    ((-c / a).clamp(0.0, 1.0), 0.0)
                } else if t_line > 1.0 {
                    // Beyond the end of B: pin t and re-solve s.
                    (((b - c) / a).clamp(0.0, 1.0), 1.0)
                } else {
                    (s_line, t_line)
                }
            }
        };

        let ca = add3(pa, scale3(da, s));
        let cb = add3(pb, scale3(db, t));
        (ca, cb, norm3(sub3(ca, cb)))
    }
}

// ---------------------------------------------------------------------------
// CapsuleCapsuleManifold
// ---------------------------------------------------------------------------

/// Capsule-capsule contact manifold generator.
///
/// Handles the segment distance computation and generates 1 or 2 contact points.
#[derive(Debug, Clone)]
pub struct CapsuleCapsuleManifold;

impl CapsuleCapsuleManifold {
    /// Compute capsule-capsule contact.
    ///
    /// Returns up to 2 contact points.
    ///
    /// # Arguments
    /// * `pa_start`, `pa_end` – endpoints of capsule A segment
    /// * `ra` – radius of capsule A
    /// * `pb_start`, `pb_end` – endpoints of capsule B segment
    /// * `rb` – radius of capsule B
    pub fn compute(
        pa_start: [f64; 3],
        pa_end: [f64; 3],
        ra: f64,
        pb_start: [f64; 3],
        pb_end: [f64; 3],
        rb: f64,
    ) -> Vec<ContactPoint> {
        let da = sub3(pa_end, pa_start);
        let db = sub3(pb_end, pb_start);
        let r = sub3(pa_start, pb_start);

        let a = dot3(da, da);
        let e = dot3(db, db);
        let f = dot3(db, r);

        let s;
        let t;
        if a <= 1e-14 && e <= 1e-14 {
            s = 0.0;
            t = 0.0;
        } else if a <= 1e-14 {
            s = 0.0;
            t = (f / e).clamp(0.0, 1.0);
        } else {
            let c = dot3(da, r);
            if e <= 1e-14 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else {
                let b = dot3(da, db);
                let denom = a * e - b * b;
                if denom.abs() > 1e-14 {
                    s = ((b * f - c * e) / denom).clamp(0.0, 1.0);
                } else {
                    s = 0.5;
                }
                t = ((b * s + f) / e).clamp(0.0, 1.0);
            }
        }

        let ca = add3(pa_start, scale3(da, s));
        let cb = add3(pb_start, scale3(db, t));
        let diff = sub3(ca, cb);
        let dist = norm3(diff);
        let total_radius = ra + rb;

        if dist >= total_radius {
            return Vec::new();
        }

        let normal = if dist > 1e-10 {
            normalize3(diff)
        } else {
            [0.0, 1.0, 0.0]
        };
        let depth = total_radius - dist;
        let pos_a = sub3(ca, scale3(normal, ra));
        let pos_b = add3(cb, scale3(normal, rb));

        let mut contacts = vec![ContactPoint::new(pos_a, pos_b, normal, depth)];

        // Check if segments are parallel and long — generate second point
        let da_len = norm3(da);
        let db_len = norm3(db);
        if da_len > 1e-5 && db_len > 1e-5 {
            let da_n = scale3(da, 1.0 / da_len);
            let db_n = scale3(db, 1.0 / db_len);
            let parallel = dot3(da_n, db_n).abs();
            if parallel > 0.9 && (s > 0.1 && s < 0.9 || t > 0.1 && t < 0.9) {
                // Add second contact at the other end
                let s2 = if s < 0.5 { 1.0 } else { 0.0 };
                let t2 = ((dot3(da_n, db_n) * (s2 - s) + t) * db_len / db_len).clamp(0.0, 1.0);
                let ca2 = add3(pa_start, scale3(da, s2));
                let cb2 = add3(pb_start, scale3(db, t2));
                let diff2 = sub3(ca2, cb2);
                let dist2 = norm3(diff2);
                if dist2 < total_radius {
                    let depth2 = total_radius - dist2;
                    let n2 = if dist2 > 1e-10 {
                        normalize3(diff2)
                    } else {
                        normal
                    };
                    let pa2 = sub3(ca2, scale3(n2, ra));
                    let pb2 = add3(cb2, scale3(n2, rb));
                    contacts.push(ContactPoint::new(pa2, pb2, n2, depth2));
                }
            }
        }

        contacts
    }
}

// ---------------------------------------------------------------------------
// MeshConvexManifold
// ---------------------------------------------------------------------------

/// Triangle mesh vs convex shape contact manifold.
///
/// Tests each triangle of the mesh against the convex shape.
#[derive(Debug, Clone)]
pub struct MeshConvexManifold {
    /// Maximum contact points.
    pub max_contacts: usize,
    /// Contact distance tolerance.
    pub tolerance: f64,
}

impl MeshConvexManifold {
    /// Create a new mesh-convex manifold generator.
    pub fn new(max_contacts: usize, tolerance: f64) -> Self {
        Self {
            max_contacts,
            tolerance,
        }
    }

    /// Compute contact between a single triangle and a sphere.
    ///
    /// Returns contact point if penetrating.
    pub fn triangle_sphere_contact(
        &self,
        tri: ([f64; 3], [f64; 3], [f64; 3]),
        tri_normal: [f64; 3],
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> Option<ContactPoint> {
        // Find closest point on triangle to sphere center
        let cp = Self::closest_on_triangle(sphere_center, tri.0, tri.1, tri.2);
        let diff = sub3(sphere_center, cp);
        let dist = norm3(diff);
        if dist >= sphere_radius {
            return None;
        }
        let normal = if dist > 1e-10 {
            normalize3(diff)
        } else {
            tri_normal
        };
        let depth = sphere_radius - dist;
        let pos_a = sub3(sphere_center, scale3(normal, sphere_radius));
        Some(ContactPoint::new(pos_a, cp, normal, depth))
    }

    fn closest_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let ap = sub3(p, a);
        let d1 = dot3(ab, ap);
        let d2 = dot3(ac, ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return a;
        }

        let bp = sub3(p, b);
        let d3 = dot3(ab, bp);
        let d4 = dot3(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            return b;
        }

        let cp_pt = sub3(p, c);
        let d5 = dot3(ab, cp_pt);
        let d6 = dot3(ac, cp_pt);
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
        if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return add3(b, scale3(sub3(c, b), w));
        }
        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        add3(a, add3(scale3(ab, v), scale3(ac, w)))
    }

    /// Generate contacts between a triangle list and a sphere.
    pub fn mesh_sphere_contacts(
        &self,
        triangles: &[([f64; 3], [f64; 3], [f64; 3])],
        normals: &[[f64; 3]],
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> Vec<ContactPoint> {
        let mut contacts = Vec::new();
        for (i, &tri) in triangles.iter().enumerate() {
            let n = normals.get(i).copied().unwrap_or([0.0, 1.0, 0.0]);
            if let Some(cp) = self.triangle_sphere_contact(tri, n, sphere_center, sphere_radius) {
                contacts.push(cp);
                if contacts.len() >= self.max_contacts {
                    break;
                }
            }
        }
        contacts
    }
}

// ---------------------------------------------------------------------------
// ManifoldReduction
// ---------------------------------------------------------------------------

/// Reduce a large contact point set to a 4-point manifold.
///
/// Uses the maximum area criterion to keep the 4 most widely spread points.
#[derive(Debug, Clone)]
pub struct ManifoldReduction;

impl ManifoldReduction {
    /// Reduce contact points to at most 4 points using maximum area criterion.
    ///
    /// The kept points are those selected by [`ManifoldReduction::reduce_indices`].
    pub fn reduce(points: &[ContactPoint]) -> Vec<ContactPoint> {
        Self::reduce_indices(points)
            .into_iter()
            .filter_map(|i| points.get(i).copied())
            .collect()
    }

    /// Indices (into `points`) of the points that [`ManifoldReduction::reduce`]
    /// keeps, in selection order.
    ///
    /// With at most [`MAX_MANIFOLD_POINTS`] input points every index is
    /// returned in order. Otherwise four distinct points are chosen:
    /// 1. the point with maximum penetration depth;
    /// 2. the point furthest from point 1;
    /// 3. the point furthest from the segment 1-2 (maximum triangle area);
    /// 4. the point that adds the largest area to triangle 1-2-3, which also
    ///    works when all points are coplanar (the usual face-face case).
    ///
    /// Ties are broken towards the lowest index, so the selection is
    /// deterministic.
    pub fn reduce_indices(points: &[ContactPoint]) -> Vec<usize> {
        if points.len() <= MAX_MANIFOLD_POINTS {
            return (0..points.len()).collect();
        }

        let mut chosen: Vec<usize> = Vec::with_capacity(MAX_MANIFOLD_POINTS);

        // Step 1: deepest point.
        let mut deepest = 0;
        let mut deepest_key = Self::score_key(points[0].depth);
        for (i, point) in points.iter().enumerate().skip(1) {
            let key = Self::score_key(point.depth);
            if key > deepest_key {
                deepest = i;
                deepest_key = key;
            }
        }
        chosen.push(deepest);
        let p0 = points[deepest].pos_a;

        // Step 2: furthest from p0.
        let second = Self::best_unchosen(points, &chosen, |p| norm_sq3(sub3(p.pos_a, p0)));
        chosen.push(second);
        let p1 = points[second].pos_a;

        // Step 3: furthest from the segment p0-p1.
        let third = Self::best_unchosen(points, &chosen, |p| {
            Self::point_to_line_dist_sq(p.pos_a, p0, p1)
        });
        chosen.push(third);
        let p2 = points[third].pos_a;

        // Step 4: largest total area of the triangles it spans with the edges
        // of the triangle p0-p1-p2.
        let fourth = Self::best_unchosen(points, &chosen, |p| {
            let q = p.pos_a;
            norm3(cross3(sub3(p1, p0), sub3(q, p0)))
                + norm3(cross3(sub3(p2, p1), sub3(q, p1)))
                + norm3(cross3(sub3(p0, p2), sub3(q, p2)))
        });
        chosen.push(fourth);

        chosen
    }

    /// Comparison key of a selection score: NaN (a point with non-finite
    /// coordinates) ranks below every number, so it is only picked when
    /// nothing else is left.
    fn score_key(value: f64) -> f64 {
        if value.is_nan() {
            f64::NEG_INFINITY
        } else {
            value
        }
    }

    /// Index of the point with the largest `score` among those not yet in
    /// `chosen` (lowest index wins ties; a NaN score never beats a number).
    fn best_unchosen(
        points: &[ContactPoint],
        chosen: &[usize],
        score: impl Fn(&ContactPoint) -> f64,
    ) -> usize {
        let mut best: Option<(usize, f64)> = None;
        for (i, point) in points.iter().enumerate() {
            if chosen.contains(&i) {
                continue;
            }
            let value = Self::score_key(score(point));
            let better = match best {
                None => true,
                Some((_, best_value)) => value > best_value,
            };
            if better {
                best = Some((i, value));
            }
        }
        best.map(|(i, _)| i).unwrap_or(0)
    }

    fn point_to_line_dist_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
        let ab = sub3(b, a);
        let ap = sub3(p, a);
        let t = dot3(ap, ab) / dot3(ab, ab).max(1e-15);
        let closest = add3(a, scale3(ab, t.clamp(0.0, 1.0)));
        norm_sq3(sub3(p, closest))
    }
}

// ---------------------------------------------------------------------------
// PersistentManifold
// ---------------------------------------------------------------------------

/// Persistent contact manifold with frame-to-frame point matching.
///
/// Maintains contact points across frames for warm-starting and stability.
#[derive(Debug, Clone)]
pub struct PersistentManifold {
    /// Current manifold.
    pub manifold: ContactManifoldEx,
    /// Distance threshold for point matching (m).
    pub match_threshold: f64,
    /// Maximum lifetime of an unmatched point (frames).
    pub max_lifetime: usize,
    /// Frame counter for each point.
    pub point_ages: Vec<usize>,
}

impl PersistentManifold {
    /// Create a new persistent manifold.
    pub fn new(body_a: u32, body_b: u32, match_threshold: f64) -> Self {
        Self {
            manifold: ContactManifoldEx::new(body_a, body_b),
            match_threshold,
            max_lifetime: 5,
            point_ages: Vec::new(),
        }
    }

    /// Update manifold with new contact points.
    ///
    /// Matches new points to existing points by proximity and transfers warm-start impulses.
    ///
    /// A point that is matched again starts a new life (age 0). An existing
    /// point without a match is kept for at most [`PersistentManifold::max_lifetime`]
    /// further updates and then dropped, so contacts of bodies that have
    /// separated expire instead of lingering forever. After reduction to
    /// [`MAX_MANIFOLD_POINTS`] points, `point_ages[k]` is the age of
    /// `manifold.points[k]`.
    pub fn update(&mut self, new_points: &[ContactPoint], frame: u64) {
        let threshold_sq = self.match_threshold * self.match_threshold;
        let mut matched = vec![false; new_points.len()];
        let mut new_manifold_pts = Vec::new();
        let mut new_ages = Vec::new();

        // Try to match new points to existing points
        for (idx, old_pt) in self.manifold.points.iter().enumerate() {
            let old_age = self.point_ages.get(idx).copied().unwrap_or(0);
            let mut best_idx = None;
            let mut best_dist = f64::MAX;
            for (j, new_pt) in new_points.iter().enumerate() {
                if matched[j] {
                    continue;
                }
                let d = norm_sq3(sub3(old_pt.pos_a, new_pt.pos_a));
                if d < threshold_sq && d < best_dist {
                    best_dist = d;
                    best_idx = Some(j);
                }
            }

            if let Some(j) = best_idx {
                matched[j] = true;
                let mut updated = new_points[j];
                // Transfer warm-start impulse
                updated.impulse = old_pt.impulse;
                new_manifold_pts.push(updated);
                new_ages.push(0);
            } else if old_age < self.max_lifetime {
                // Keep old point if not too old
                new_manifold_pts.push(*old_pt);
                new_ages.push(old_age + 1);
            }
        }

        // Add unmatched new points
        for (j, &pt) in new_points.iter().enumerate() {
            if !matched[j] {
                new_manifold_pts.push(pt);
                new_ages.push(0);
            }
        }

        // Reduce to max 4 points, carrying every kept point's age along.
        let kept = ManifoldReduction::reduce_indices(&new_manifold_pts);
        self.manifold.points = kept
            .iter()
            .filter_map(|&i| new_manifold_pts.get(i).copied())
            .collect();
        self.point_ages = kept
            .iter()
            .filter_map(|&i| new_ages.get(i).copied())
            .collect();
        self.manifold.last_frame = frame;
    }

    /// Number of currently tracked contact points.
    pub fn num_points(&self) -> usize {
        self.manifold.num_points()
    }

    /// Get accumulated impulse for warm-starting.
    pub fn warm_start_impulse(&self, point_idx: usize) -> f64 {
        self.manifold
            .points
            .get(point_idx)
            .map(|p| p.impulse)
            .unwrap_or(0.0)
    }
}

// ---------------------------------------------------------------------------
// FeatureContact
// ---------------------------------------------------------------------------

/// Contact feature type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureType {
    /// Vertex contact.
    Vertex,
    /// Edge contact.
    Edge,
    /// Face contact.
    Face,
}

/// Contact feature with type tracking for warm-starting.
#[derive(Debug, Clone)]
pub struct FeatureContact {
    /// Contact point.
    pub contact: ContactPoint,
    /// Feature type on body A.
    pub feature_type_a: FeatureType,
    /// Feature type on body B.
    pub feature_type_b: FeatureType,
    /// Previous feature type on A (for transition detection).
    pub prev_feature_a: FeatureType,
    /// Previous feature type on B.
    pub prev_feature_b: FeatureType,
}

impl FeatureContact {
    /// Create a new feature contact.
    pub fn new(
        contact: ContactPoint,
        feature_type_a: FeatureType,
        feature_type_b: FeatureType,
    ) -> Self {
        Self {
            contact,
            feature_type_a,
            feature_type_b,
            prev_feature_a: feature_type_a,
            prev_feature_b: feature_type_b,
        }
    }

    /// Check if a feature type transition occurred.
    pub fn has_transition(&self) -> bool {
        self.feature_type_a != self.prev_feature_a || self.feature_type_b != self.prev_feature_b
    }

    /// Update previous feature types.
    pub fn update_prev(&mut self) {
        self.prev_feature_a = self.feature_type_a;
        self.prev_feature_b = self.feature_type_b;
    }

    /// Classify feature type based on contact normal and geometry.
    ///
    /// Returns FeatureType based on alignment between normal and face/edge axes.
    pub fn classify_feature(normal: [f64; 3], face_normals: &[[f64; 3]]) -> FeatureType {
        let mut best_dot = 0.0f64;
        for &fn_ in face_normals {
            let d = dot3(normal, fn_).abs();
            if d > best_dot {
                best_dot = d;
            }
        }
        if best_dot > 0.95 {
            FeatureType::Face
        } else if best_dot > 0.7 {
            FeatureType::Edge
        } else {
            FeatureType::Vertex
        }
    }
}

// ---------------------------------------------------------------------------
// ContactNormalSmoothing
// ---------------------------------------------------------------------------

/// Contact normal smoother for mesh collision.
///
/// Smooths contact normals using Phong normal interpolation
/// or Garland-Heckbert normal averaging.
#[derive(Debug, Clone)]
pub struct ContactNormalSmoothing {
    /// Method for normal smoothing.
    pub method: NormalSmoothingMethod,
    /// Angle threshold for crease detection (radians).
    pub crease_angle: f64,
}

/// Normal smoothing method.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalSmoothingMethod {
    /// Phong normal interpolation (barycentric interpolation of vertex normals).
    Phong,
    /// Garland-Heckbert area-weighted normal averaging.
    Garland,
    /// No smoothing (flat normals).
    Flat,
}

impl ContactNormalSmoothing {
    /// Create a new contact normal smoother.
    pub fn new(method: NormalSmoothingMethod, crease_angle: f64) -> Self {
        Self {
            method,
            crease_angle,
        }
    }

    /// Smooth contact normal at barycentric coordinates (u, v) on a triangle.
    ///
    /// # Arguments
    /// * `face_normal` – flat face normal
    /// * `vertex_normals` – per-vertex normals (a, b, c)
    /// * `u`, `v` – barycentric coordinates (w = 1-u-v)
    pub fn smooth_normal(
        &self,
        face_normal: [f64; 3],
        vertex_normals: ([f64; 3], [f64; 3], [f64; 3]),
        u: f64,
        v: f64,
    ) -> [f64; 3] {
        match self.method {
            NormalSmoothingMethod::Flat => face_normal,
            NormalSmoothingMethod::Phong => {
                let w = 1.0 - u - v;
                let (na, nb, nc) = vertex_normals;
                let interpolated = [
                    w * na[0] + u * nb[0] + v * nc[0],
                    w * na[1] + u * nb[1] + v * nc[1],
                    w * na[2] + u * nb[2] + v * nc[2],
                ];
                normalize3(interpolated)
            }
            NormalSmoothingMethod::Garland => {
                // Area-weighted average with face normal
                let (na, nb, nc) = vertex_normals;
                let avg = [
                    (na[0] + nb[0] + nc[0]) / 3.0,
                    (na[1] + nb[1] + nc[1]) / 3.0,
                    (na[2] + nb[2] + nc[2]) / 3.0,
                ];
                // Blend with face normal based on crease angle
                let dot = dot3(normalize3(avg), face_normal);
                if dot.acos() < self.crease_angle {
                    normalize3(avg)
                } else {
                    face_normal
                }
            }
        }
    }

    /// Compute area-weighted vertex normals for a mesh.
    ///
    /// Returns vertex normals for each vertex index.
    pub fn compute_vertex_normals(vertices: &[[f64; 3]], indices: &[[usize; 3]]) -> Vec<[f64; 3]> {
        let n = vertices.len();
        let mut normals = vec![[0.0f64; 3]; n];

        for tri in indices {
            let a = vertices[tri[0]];
            let b = vertices[tri[1]];
            let c = vertices[tri[2]];
            let ab = sub3(b, a);
            let ac = sub3(c, a);
            let area_normal = cross3(ab, ac); // magnitude = 2 * area
            for &vi in tri {
                normals[vi][0] += area_normal[0];
                normals[vi][1] += area_normal[1];
                normals[vi][2] += area_normal[2];
            }
        }

        normals.iter_mut().for_each(|n| *n = normalize3(*n));
        normals
    }
}

// ---------------------------------------------------------------------------
// ContactCaching
// ---------------------------------------------------------------------------

/// LRU cache for body-pair contact manifolds.
///
/// Stores persistent manifolds indexed by (body_a, body_b) pair.
#[derive(Debug, Clone)]
pub struct ContactCaching {
    /// Cached manifolds by body pair.
    pub cache: std::collections::HashMap<(u32, u32), PersistentManifold>,
    /// Maximum cache size.
    pub max_size: usize,
    /// Frame counter for LRU eviction.
    pub frame: u64,
}

impl ContactCaching {
    /// Create a new contact cache.
    pub fn new(max_size: usize) -> Self {
        Self {
            cache: std::collections::HashMap::new(),
            max_size,
            frame: 0,
        }
    }

    /// Get or create a persistent manifold for a body pair.
    pub fn get_or_create(&mut self, body_a: u32, body_b: u32) -> &mut PersistentManifold {
        let key = (body_a.min(body_b), body_a.max(body_b));
        if !self.cache.contains_key(&key) {
            if self.cache.len() >= self.max_size {
                self.evict_oldest();
            }
            self.cache
                .insert(key, PersistentManifold::new(body_a, body_b, 0.01));
        }
        self.cache
            .entry(key)
            .or_insert_with(|| PersistentManifold::new(body_a, body_b, 0.01))
    }

    /// Evict the oldest (by last_frame) manifold.
    fn evict_oldest(&mut self) {
        let oldest_key = self
            .cache
            .iter()
            .min_by_key(|(_, m)| m.manifold.last_frame)
            .map(|(k, _)| *k);
        if let Some(k) = oldest_key {
            self.cache.remove(&k);
        }
    }

    /// Advance frame counter.
    pub fn advance_frame(&mut self) {
        self.frame += 1;
    }

    /// Remove stale manifolds (not updated in last N frames).
    pub fn purge_stale(&mut self, max_age: u64) {
        let frame = self.frame;
        self.cache
            .retain(|_, m| frame.saturating_sub(m.manifold.last_frame) <= max_age);
    }

    /// Number of cached manifolds.
    pub fn num_cached(&self) -> usize {
        self.cache.len()
    }
}

// ---------------------------------------------------------------------------
// ContactEventSystem
// ---------------------------------------------------------------------------

/// Contact event type.
#[derive(Debug, Clone, PartialEq)]
pub enum ContactEvent {
    /// Contact began (first contact point appeared).
    Begin(u32, u32),
    /// Contact ended (all contact points removed).
    End(u32, u32),
    /// Contact persisted (same bodies, still in contact).
    Persist(u32, u32),
}

/// Contact event listener trait.
pub trait ContactListener {
    /// Called when a contact event is emitted.
    fn on_contact_event(&mut self, event: ContactEvent);
}

/// Event system for contact begin/end notifications.
#[derive(Debug, Clone)]
pub struct ContactEventSystem {
    /// Current contact pairs (body_a, body_b).
    pub active_contacts: std::collections::HashSet<(u32, u32)>,
    /// Pending events this frame.
    pub pending_events: Vec<ContactEvent>,
}

impl ContactEventSystem {
    /// Create a new contact event system.
    pub fn new() -> Self {
        Self {
            active_contacts: std::collections::HashSet::new(),
            pending_events: Vec::new(),
        }
    }

    /// Update contact set and emit events.
    ///
    /// # Arguments
    /// * `new_contacts` – all currently active contact pairs this frame
    pub fn update(&mut self, new_contacts: &[(u32, u32)]) {
        let new_set: std::collections::HashSet<_> = new_contacts
            .iter()
            .map(|&(a, b)| (a.min(b), a.max(b)))
            .collect();

        self.pending_events.clear();

        // Emit End for pairs no longer in contact
        for &pair in &self.active_contacts {
            if !new_set.contains(&pair) {
                self.pending_events.push(ContactEvent::End(pair.0, pair.1));
            }
        }

        // Emit Begin for new pairs, Persist for existing
        for &pair in &new_set {
            if !self.active_contacts.contains(&pair) {
                self.pending_events
                    .push(ContactEvent::Begin(pair.0, pair.1));
            } else {
                self.pending_events
                    .push(ContactEvent::Persist(pair.0, pair.1));
            }
        }

        self.active_contacts = new_set;
    }

    /// Dispatch pending events to a listener.
    pub fn dispatch<L: ContactListener>(&self, listener: &mut L) {
        for event in &self.pending_events {
            listener.on_contact_event(event.clone());
        }
    }

    /// Number of active contact pairs.
    pub fn num_active(&self) -> usize {
        self.active_contacts.len()
    }

    /// Number of pending events.
    pub fn num_events(&self) -> usize {
        self.pending_events.len()
    }
}

impl Default for ContactEventSystem {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contact_point_new() {
        let cp = ContactPoint::new([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0], 0.1);
        assert!((cp.depth - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_contact_point_midpoint() {
        let cp = ContactPoint::new([0.0; 3], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1);
        let mid = cp.midpoint();
        assert!((mid[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clip_polygon_by_plane_keep_all() {
        let poly = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        // All points above plane y=0 (pointing +y)
        let clipped = clip_polygon_by_plane(&poly, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(clipped.len(), 3, "all points should be kept");
    }

    #[test]
    fn test_clip_polygon_by_plane_half() {
        let poly = vec![
            [0.0, -1.0, 0.0_f64],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        // Clip against y>=0
        let clipped = clip_polygon_by_plane(&poly, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            !clipped.is_empty(),
            "some points should remain after clipping"
        );
        for p in &clipped {
            assert!(
                p[1] >= -1e-10,
                "all remaining points should satisfy constraint"
            );
        }
    }

    #[test]
    fn test_sat_face_axis_separated() {
        let ref_pts = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let inc_pts = vec![[3.0_f64, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let sep = sat_face_axis([1.0, 0.0, 0.0], &ref_pts, &inc_pts);
        assert!(sep > 0.0, "should detect separation: {}", sep);
    }

    #[test]
    fn test_sat_face_axis_overlapping() {
        let ref_pts = vec![[0.0_f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let inc_pts = vec![[1.0_f64, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let sep = sat_face_axis([1.0, 0.0, 0.0], &ref_pts, &inc_pts);
        assert!(sep < 0.0, "should detect overlap: {}", sep);
    }

    #[test]
    fn test_find_incident_face() {
        let normals = vec![[0.0, 1.0, 0.0_f64], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]];
        let (idx, _) = find_incident_face(&normals, [0.0, 1.0, 0.0]);
        // Most anti-parallel: [0,-1,0] (dot = -1)
        assert_eq!(idx, 1, "most anti-parallel face should be found");
    }

    #[test]
    fn test_polyhedron_face_face_contact() {
        let pg = PolyhedronContact::new(4);
        // Two overlapping quads
        let ref_face = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let inc_face = vec![
            [0.2_f64, 0.2, -0.1],
            [0.8, 0.2, -0.1],
            [0.8, 0.8, -0.1],
            [0.2, 0.8, -0.1],
        ];
        let contacts = pg.face_face_contact(&ref_face, &inc_face, [0.0, 0.0, 1.0], 0.2);
        // The incident quad lies 0.1 below the reference plane, fully inside
        // the reference quad: all four corners are contacts of depth 0.1.
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        for (contact, corner) in contacts.iter().zip(inc_face.iter()) {
            assert_eq!(contact.normal, [0.0, 0.0, 1.0]);
            assert!((contact.depth - 0.1).abs() < 1e-12, "{}", contact.depth);
            assert_eq!(contact.pos_b, *corner);
            assert!(norm3(sub3(contact.pos_a, [corner[0], corner[1], 0.0])) < 1e-12);
        }
    }

    #[test]
    fn test_polyhedron_face_face_contact_clips_and_applies_the_margin() {
        let square = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let up = [0.0, 0.0, 1.0];
        let shifted = |dx: f64, z: f64| square.map(|p| [p[0] + dx, p[1], z]);

        // Half of the incident quad hangs over the reference quad: the clipped
        // polygon is [0.5, 1] x [0, 1].
        let pg = PolyhedronContact::new(8);
        let contacts = pg.face_face_contact(&square, &shifted(0.5, -0.1), up, 0.2);
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        for contact in &contacts {
            assert!(contact.pos_b[0] >= 0.5 - 1e-12 && contact.pos_b[0] <= 1.0 + 1e-12);
            assert!((contact.depth - 0.1).abs() < 1e-12);
        }

        // Above the reference plane but inside the margin: kept with a
        // negative depth; beyond the margin: dropped.
        let near = pg.face_face_contact(&square, &shifted(0.0, 0.1), up, 0.2);
        assert_eq!(near.len(), 4);
        assert!(near.iter().all(|c| (c.depth + 0.1).abs() < 1e-12));
        assert!(
            pg.face_face_contact(&square, &shifted(0.0, 0.3), up, 0.2)
                .is_empty()
        );

        // The number of contacts is capped, and empty input gives no contact.
        let capped =
            PolyhedronContact::new(2).face_face_contact(&square, &shifted(0.0, -0.1), up, 0.2);
        assert_eq!(capped.len(), 2);
        assert!(pg.face_face_contact(&[], &square, up, 0.2).is_empty());
        assert!(pg.face_face_contact(&square, &[], up, 0.2).is_empty());
    }

    #[test]
    fn test_polyhedron_edge_edge_closest() {
        let pg = PolyhedronContact::new(4);
        let pa = [0.0, 0.0, 0.0_f64];
        let da = [1.0, 0.0, 0.0];
        let pb = [0.5, 0.5, -1.0_f64];
        let db = [0.0, 0.0, 1.0];
        let (ca, cb, dist) = pg.edge_edge_closest(pa, da, pb, db);
        // Edge A runs along x at y = z = 0; edge B ends at (0.5, 0.5, 0): the
        // closest points are the middle of A and the end point of B.
        assert!(norm3(sub3(ca, [0.5, 0.0, 0.0])) < 1e-12, "{ca:?}");
        assert!(norm3(sub3(cb, [0.5, 0.5, 0.0])) < 1e-12, "{cb:?}");
        assert!((dist - 0.5).abs() < 1e-12, "distance {dist}");
    }

    /// Distance from `point` to the segment `start + t * dir`, `t` in `[0, 1]`.
    fn point_segment_distance(point: [f64; 3], start: [f64; 3], dir: [f64; 3]) -> f64 {
        let length_sq = dot3(dir, dir);
        let t = if length_sq > 0.0 {
            (dot3(sub3(point, start), dir) / length_sq).clamp(0.0, 1.0)
        } else {
            0.0
        };
        norm3(sub3(point, add3(start, scale3(dir, t))))
    }

    /// Exact distance between two segments: the minimum of a convex quadratic
    /// over the parameter square is either at its interior critical point or on
    /// the boundary, and the boundary consists of four point-segment problems.
    fn exact_segment_distance(pa: [f64; 3], da: [f64; 3], pb: [f64; 3], db: [f64; 3]) -> f64 {
        let mut best = point_segment_distance(pa, pb, db)
            .min(point_segment_distance(add3(pa, da), pb, db))
            .min(point_segment_distance(pb, pa, da))
            .min(point_segment_distance(add3(pb, db), pa, da));
        let r = sub3(pa, pb);
        let (a, b, e) = (dot3(da, da), dot3(da, db), dot3(db, db));
        let (c, f) = (dot3(da, r), dot3(db, r));
        let denom = a * e - b * b;
        if denom > 1e-12 {
            let s = (b * f - c * e) / denom;
            let t = (a * f - b * c) / denom;
            if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
                let gap = sub3(add3(pa, scale3(da, s)), add3(pb, scale3(db, t)));
                best = best.min(norm3(gap));
            }
        }
        best
    }

    #[test]
    fn test_polyhedron_edge_edge_closest_matches_the_exact_distance() {
        let pg = PolyhedronContact::new(1);
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut unit = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) as f64 / (1u64 << 31) as f64 * 2.0 - 1.0
        };
        let mut vec3 = |scale: f64| [unit() * scale, unit() * scale, unit() * scale];
        for case in 0..2000 {
            let pa = vec3(2.0);
            let pb = vec3(2.0);
            let da = vec3(2.0);
            // Every tenth pair is exactly parallel, every tenth pair has
            // perpendicular edge directions.
            let db = match case % 10 {
                0 => scale3(da, 0.7),
                5 => [da[1], -da[0], 0.0],
                _ => vec3(2.0),
            };
            let (ca, cb, dist) = pg.edge_edge_closest(pa, da, pb, db);
            let exact = exact_segment_distance(pa, da, pb, db);
            assert!(
                (dist - exact).abs() < 1e-9,
                "case {case}: distance {dist} vs exact {exact}"
            );
            assert!((norm3(sub3(ca, cb)) - dist).abs() < 1e-12);
            // The closest points lie on their segments.
            assert!(point_segment_distance(ca, pa, da) < 1e-9, "case {case}");
            assert!(point_segment_distance(cb, pb, db) < 1e-9, "case {case}");
        }
    }

    #[test]
    fn test_polyhedron_edge_edge_closest_degenerate_edges() {
        let pg = PolyhedronContact::new(1);
        let origin = [0.0, 0.0, 0.0_f64];
        let along_x = [2.0, 0.0, 0.0];
        let point_b = [1.0, 3.0, 0.0];
        // Edge B is a point: the closest point of A is its projection.
        let (ca, cb, dist) = pg.edge_edge_closest(origin, along_x, point_b, [0.0; 3]);
        assert!(norm3(sub3(ca, [1.0, 0.0, 0.0])) < 1e-12);
        assert_eq!(cb, point_b);
        assert!((dist - 3.0).abs() < 1e-12);
        // Edge A is a point.
        let (ca, cb, dist) = pg.edge_edge_closest(point_b, [0.0; 3], origin, along_x);
        assert_eq!(ca, point_b);
        assert!(norm3(sub3(cb, [1.0, 0.0, 0.0])) < 1e-12);
        assert!((dist - 3.0).abs() < 1e-12);
        // Both are points.
        let (ca, cb, dist) = pg.edge_edge_closest(origin, [0.0; 3], point_b, [0.0; 3]);
        assert_eq!((ca, cb), (origin, point_b));
        assert!((dist - norm3(point_b)).abs() < 1e-12);
    }

    #[test]
    fn test_box_box_manifold_separated() {
        let bbox = BoxBoxManifold::new(0.01);
        let id = [[1.0, 0.0, 0.0_f64], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let contacts = bbox.compute([0.0; 3], [1.0; 3], id, [5.0, 0.0, 0.0], [1.0; 3], id);
        assert!(
            contacts.is_empty(),
            "separated boxes should have no contacts"
        );
    }

    #[test]
    fn test_capsule_capsule_separated() {
        let contacts = CapsuleCapsuleManifold::compute(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
            0.5,
        );
        assert!(
            contacts.is_empty(),
            "separated capsules should have no contacts"
        );
    }

    #[test]
    fn test_capsule_capsule_intersecting() {
        let contacts = CapsuleCapsuleManifold::compute(
            [0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            0.5,
            [0.3, 0.0, 0.0],
            [0.3, 2.0, 0.0],
            0.5,
        );
        assert!(
            !contacts.is_empty(),
            "close capsules should generate contact"
        );
        assert!(contacts[0].depth > 0.0, "contact depth should be positive");
    }

    #[test]
    fn test_mesh_convex_sphere_contact() {
        let manifold_gen = MeshConvexManifold::new(8, 1e-3);
        let tri = ([0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let n = [0.0, 0.0, 1.0_f64];
        let sphere_center = [0.25, 0.25, 0.1]; // slightly above triangle
        let contact = manifold_gen.triangle_sphere_contact(tri, n, sphere_center, 0.5);
        assert!(contact.is_some(), "sphere above triangle should contact");
    }

    #[test]
    fn test_manifold_reduction_under_limit() {
        let pts: Vec<_> = (0..3)
            .map(|i| {
                ContactPoint::new(
                    [i as f64, 0.0, 0.0],
                    [i as f64, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    0.1,
                )
            })
            .collect();
        let reduced = ManifoldReduction::reduce(&pts);
        assert_eq!(reduced.len(), 3, "under 4 points should not reduce");
    }

    #[test]
    fn test_manifold_reduction_over_limit() {
        let pts: Vec<_> = (0..10)
            .map(|i| {
                let d = (i as f64) * 0.5;
                ContactPoint::new(
                    [i as f64, 0.0, 0.0],
                    [i as f64, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    d,
                )
            })
            .collect();
        let reduced = ManifoldReduction::reduce(&pts);
        assert!(
            reduced.len() <= MAX_MANIFOLD_POINTS,
            "reduced should have at most 4 points"
        );
    }

    #[test]
    fn test_persistent_manifold_update() {
        let mut pm = PersistentManifold::new(1, 2, 0.1);
        let pts = vec![ContactPoint::new([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0], 0.1)];
        pm.update(&pts, 1);
        assert_eq!(pm.num_points(), 1);
    }

    #[test]
    fn test_persistent_manifold_warm_start() {
        let mut pm = PersistentManifold::new(1, 2, 0.1);
        let mut cp = ContactPoint::new([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0], 0.1);
        cp.impulse = 5.0;
        pm.manifold.add_point(cp);
        pm.point_ages.push(0);
        let impulse = pm.warm_start_impulse(0);
        assert!((impulse - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_feature_contact_classification() {
        let face_normals = vec![[0.0, 1.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let ft = FeatureContact::classify_feature([0.0, 1.0, 0.0], &face_normals);
        assert_eq!(
            ft,
            FeatureType::Face,
            "aligned with face normal should be Face"
        );
    }

    #[test]
    fn test_contact_normal_smoothing_phong() {
        let smoother = ContactNormalSmoothing::new(NormalSmoothingMethod::Phong, 0.3);
        let face_n = [0.0, 0.0, 1.0_f64];
        let vn = ([0.0, 0.0, 1.0_f64], [0.0, 0.1, 0.995], [0.0, -0.1, 0.995]);
        let smooth_n = smoother.smooth_normal(face_n, vn, 0.5, 0.0);
        let len = norm3(smooth_n);
        assert!(
            (len - 1.0).abs() < 1e-6,
            "smoothed normal should be unit length"
        );
    }

    #[test]
    fn test_contact_normal_vertex_normals() {
        let verts = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let indices = vec![[0, 1, 2]];
        let normals = ContactNormalSmoothing::compute_vertex_normals(&verts, &indices);
        assert_eq!(normals.len(), 3);
        for n in &normals {
            let len = norm3(*n);
            assert!(
                (len - 1.0).abs() < 1e-6 || len < 1e-10,
                "vertex normals should be unit: {}",
                len
            );
        }
    }

    #[test]
    fn test_contact_caching_get_or_create() {
        let mut cache = ContactCaching::new(10);
        let m = cache.get_or_create(1, 2);
        assert_eq!(m.manifold.body_a, 1);
    }

    #[test]
    fn test_contact_caching_size_limit() {
        let mut cache = ContactCaching::new(2);
        cache.get_or_create(1, 2).manifold.last_frame = 1;
        cache.get_or_create(3, 4).manifold.last_frame = 2;
        cache.get_or_create(5, 6).manifold.last_frame = 3;
        assert!(cache.num_cached() <= 2, "cache should not exceed max_size");
    }

    #[test]
    fn test_contact_event_system_begin() {
        let mut evs = ContactEventSystem::new();
        evs.update(&[(1, 2)]);
        let begin = evs
            .pending_events
            .iter()
            .any(|e| matches!(e, ContactEvent::Begin(1, 2) | ContactEvent::Begin(2, 1)));
        assert!(begin, "should emit Begin event for new contact");
    }

    #[test]
    fn test_contact_event_system_end() {
        let mut evs = ContactEventSystem::new();
        evs.update(&[(1, 2)]);
        evs.update(&[]); // contact ended
        let end = evs
            .pending_events
            .iter()
            .any(|e| matches!(e, ContactEvent::End(_, _)));
        assert!(end, "should emit End event when contact is removed");
    }

    #[test]
    fn test_contact_event_system_persist() {
        let mut evs = ContactEventSystem::new();
        evs.update(&[(1, 2)]);
        evs.update(&[(1, 2)]); // still in contact
        let persist = evs
            .pending_events
            .iter()
            .any(|e| matches!(e, ContactEvent::Persist(_, _)));
        assert!(persist, "should emit Persist event for ongoing contact");
    }

    #[test]
    fn test_contact_manifold_add_points() {
        let mut mf = ContactManifoldEx::new(1, 2);
        for i in 0..5 {
            mf.add_point(ContactPoint::new(
                [i as f64, 0.0, 0.0],
                [0.0; 3],
                [0.0, 1.0, 0.0],
                0.1,
            ));
        }
        assert!(
            mf.num_points() <= MAX_MANIFOLD_POINTS,
            "manifold should not exceed max points"
        );
    }

    #[test]
    fn test_contact_manifold_prune() {
        let mut mf = ContactManifoldEx::new(1, 2);
        mf.add_point(ContactPoint::new([0.0; 3], [0.0; 3], [0.0, 1.0, 0.0], 0.1));
        mf.add_point(ContactPoint::new(
            [1.0, 0.0, 0.0],
            [0.0; 3],
            [0.0, 1.0, 0.0],
            -0.5,
        ));
        mf.prune_inactive();
        assert_eq!(mf.num_points(), 1, "negative depth point should be pruned");
    }

    // -- ManifoldReduction ---------------------------------------------------

    fn point_at(x: f64, y: f64, depth: f64) -> ContactPoint {
        ContactPoint::new([x, y, 0.0], [x, y, 0.0], [0.0, 0.0, 1.0], depth)
    }

    #[test]
    fn test_reduce_indices_identity_under_limit() {
        let pts: Vec<_> = (0..MAX_MANIFOLD_POINTS)
            .map(|i| point_at(i as f64, 0.0, 0.1))
            .collect();
        assert_eq!(ManifoldReduction::reduce_indices(&pts), vec![0, 1, 2, 3]);
        assert!(ManifoldReduction::reduce_indices(&[]).is_empty());
    }

    #[test]
    fn test_reduce_indices_distinct_for_coplanar_octagon() {
        // Eight coplanar points of equal depth: the reduction must keep four
        // distinct points that span a large quadrilateral.
        let pts: Vec<_> = (0..8)
            .map(|i| {
                let angle = std::f64::consts::FRAC_PI_4 * i as f64;
                point_at(angle.cos(), angle.sin(), 0.05)
            })
            .collect();
        let kept = ManifoldReduction::reduce_indices(&pts);
        assert_eq!(kept.len(), MAX_MANIFOLD_POINTS);
        let mut unique = kept.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), MAX_MANIFOLD_POINTS, "duplicates in {kept:?}");

        // Area of the kept quadrilateral (shoelace over the angular order).
        let mut by_angle: Vec<[f64; 3]> = kept.iter().map(|&i| pts[i].pos_a).collect();
        by_angle.sort_by(|a, b| {
            a[1].atan2(a[0])
                .partial_cmp(&b[1].atan2(b[0]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut twice_area = 0.0;
        for k in 0..by_angle.len() {
            let a = by_angle[k];
            let b = by_angle[(k + 1) % by_angle.len()];
            twice_area += a[0] * b[1] - b[0] * a[1];
        }
        // The best quadrilateral of a unit-circle octagon is a square: area 2.
        assert!(
            twice_area.abs() * 0.5 > 1.9,
            "kept points span too little area: {}",
            twice_area.abs() * 0.5
        );
    }

    #[test]
    fn test_reduce_indices_keeps_deepest_and_matches_reduce() {
        let mut pts: Vec<_> = (0..9)
            .map(|i| point_at((i % 3) as f64, (i / 3) as f64, 0.01))
            .collect();
        pts[4].depth = 0.5;
        let kept = ManifoldReduction::reduce_indices(&pts);
        assert_eq!(kept[0], 4, "deepest point is selected first");
        let reduced = ManifoldReduction::reduce(&pts);
        assert_eq!(reduced.len(), kept.len());
        for (point, &i) in reduced.iter().zip(kept.iter()) {
            assert_eq!(*point, pts[i]);
        }
    }

    #[test]
    fn test_reduce_indices_never_prefers_nan_candidates() {
        // The first point carries NaN data. As the first candidate of every
        // selection step it must not stick just because `x > NaN` is false.
        let mut pts: Vec<_> = (0..8)
            .map(|i| {
                let angle = std::f64::consts::FRAC_PI_4 * i as f64;
                point_at(angle.cos(), angle.sin(), 0.05)
            })
            .collect();
        pts[0].depth = f64::NAN;
        pts[0].pos_a = [f64::NAN, f64::NAN, 0.0];
        let kept = ManifoldReduction::reduce_indices(&pts);
        assert_eq!(kept.len(), MAX_MANIFOLD_POINTS);
        let mut unique = kept.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), MAX_MANIFOLD_POINTS, "duplicates in {kept:?}");
        assert!(!kept.contains(&0), "NaN point was kept: {kept:?}");
        for &i in &kept {
            assert!(pts[i].pos_a.iter().all(|c| c.is_finite()));
        }
    }

    // -- PersistentManifold lifetime ----------------------------------------

    #[test]
    fn test_persistent_manifold_unmatched_point_expires() {
        let mut pm = PersistentManifold::new(1, 2, 0.1);
        assert_eq!(pm.max_lifetime, 5);
        pm.update(&[point_at(0.0, 0.0, 0.1)], 1);
        assert_eq!(pm.num_points(), 1);
        assert_eq!(pm.point_ages, vec![0]);

        // The bodies separate: no new contacts arrive. The ghost point ages
        // and survives exactly `max_lifetime` further updates.
        for frame in 2..=6 {
            pm.update(&[], frame);
            assert_eq!(pm.num_points(), 1, "still alive at frame {frame}");
            assert_eq!(pm.point_ages, vec![(frame - 1) as usize]);
        }
        pm.update(&[], 7);
        assert_eq!(pm.num_points(), 0, "ghost contact must expire");
        assert!(pm.point_ages.is_empty());
        assert_eq!(pm.manifold.last_frame, 7);
    }

    #[test]
    fn test_persistent_manifold_rematched_point_resets_age_and_keeps_impulse() {
        let mut pm = PersistentManifold::new(1, 2, 0.1);
        pm.update(&[point_at(0.0, 0.0, 0.1)], 1);
        pm.manifold.points[0].impulse = 3.5;
        for frame in 2..=4 {
            pm.update(&[], frame);
        }
        assert_eq!(pm.point_ages, vec![3]);

        // The contact reappears close to where it was.
        pm.update(&[point_at(0.02, 0.0, 0.2)], 5);
        assert_eq!(pm.num_points(), 1);
        assert_eq!(pm.point_ages, vec![0], "re-matched point starts a new life");
        assert!((pm.warm_start_impulse(0) - 3.5).abs() < 1e-12);
        assert!((pm.manifold.points[0].depth - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_persistent_manifold_ages_follow_reduction() {
        let mut pm = PersistentManifold::new(1, 2, 0.05);
        // One contact at a corner of the future patch.
        pm.update(&[point_at(1.0, 1.0, 0.1)], 1);
        // Five unrelated new contacts arrive (the old one stays unmatched):
        // six candidates are reduced to four.
        let fresh = [
            point_at(-1.0, -1.0, 0.5),
            point_at(1.0, -1.0, 0.2),
            point_at(-1.0, 1.0, 0.2),
            point_at(0.0, 0.0, 0.2),
            point_at(0.1, 0.2, 0.2),
        ];
        pm.update(&fresh, 2);
        assert_eq!(pm.num_points(), MAX_MANIFOLD_POINTS);
        assert_eq!(pm.point_ages.len(), MAX_MANIFOLD_POINTS);
        for (point, &age) in pm.manifold.points.iter().zip(pm.point_ages.iter()) {
            let is_old_corner = point.pos_a == [1.0, 1.0, 0.0];
            assert_eq!(
                age,
                if is_old_corner { 1 } else { 0 },
                "wrong age for point at {:?}",
                point.pos_a
            );
        }
        assert!(
            pm.manifold
                .points
                .iter()
                .any(|p| p.pos_a == [1.0, 1.0, 0.0]),
            "the extreme old corner must survive the reduction"
        );
    }

    #[test]
    fn test_persistent_manifold_tolerates_missing_ages() {
        // Points added directly to the public manifold without an age entry
        // must not silently disappear on the next update.
        let mut pm = PersistentManifold::new(1, 2, 0.1);
        pm.manifold.add_point(point_at(0.0, 0.0, 0.1));
        pm.update(&[], 1);
        assert_eq!(pm.num_points(), 1);
        assert_eq!(pm.point_ages, vec![1]);
    }
}
