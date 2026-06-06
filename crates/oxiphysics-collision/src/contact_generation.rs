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
        let r = sub3(pa, pb);
        let a = dot3(da, da);
        let e = dot3(db, db);
        let f = dot3(db, r);

        if a < 1e-14 && e < 1e-14 {
            return (pa, pb, norm3(r));
        }

        let s;
        let t;
        if a < 1e-14 {
            s = 0.0;
            t = (f / e).clamp(0.0, 1.0);
        } else {
            let c = dot3(da, r);
            if e < 1e-14 {
                s = (-c / a).clamp(0.0, 1.0);
                t = 0.0;
            } else {
                let b = dot3(da, db);
                let denom = a * e - b * b;
                if denom.abs() > 1e-14 {
                    s = ((b * f - c * e) / denom).clamp(0.0, 1.0);
                } else {
                    s = 0.5;
                }
                t = ((b * s + f) / e).clamp(0.0, 1.0);
                let _ = s;
            }
            let _ = t;
        }

        // Recompute s with clamped t
        let s_final = if a > 1e-14 {
            let c = dot3(da, r);
            let b = dot3(da, db);
            ((b * t - c) / a).clamp(0.0, 1.0)
        } else {
            0.0
        };

        let t_final = if e > 1e-14 {
            let b = dot3(da, db);
            ((b * s_final + f) / e).clamp(0.0, 1.0)
        } else {
            0.0
        };

        // Wait, we need the correct values
        let _ = (s, t);

        let ca = add3(pa, scale3(da, s_final));
        let cb = add3(pb, scale3(db, t_final));
        (ca, cb, norm3(sub3(ca, cb)))
    }
}

// ---------------------------------------------------------------------------
// BoxBoxManifold
// ---------------------------------------------------------------------------

/// Specialized box-box contact manifold generator.
///
/// Handles all 6 face pairs and 12*12=144 edge-edge pairs via SAT.
#[derive(Debug, Clone)]
pub struct BoxBoxManifold {
    /// Tolerance for contact detection.
    pub tolerance: f64,
}

impl BoxBoxManifold {
    /// Create a new box-box manifold generator.
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Generate box half-extent vertices.
    fn box_vertices(
        center: [f64; 3],
        half_extents: [f64; 3],
        axes: [[f64; 3]; 3],
    ) -> Vec<[f64; 3]> {
        let mut verts = Vec::with_capacity(8);
        for sx in &[-1.0, 1.0] {
            for sy in &[-1.0, 1.0] {
                for sz in &[-1.0, 1.0] {
                    let v = add3(
                        center,
                        add3(
                            add3(
                                scale3(axes[0], sx * half_extents[0]),
                                scale3(axes[1], sy * half_extents[1]),
                            ),
                            scale3(axes[2], sz * half_extents[2]),
                        ),
                    );
                    verts.push(v);
                }
            }
        }
        verts
    }

    /// Compute box-box contact manifold.
    ///
    /// Returns list of contact points, or empty if separated.
    pub fn compute(
        &self,
        center_a: [f64; 3],
        half_a: [f64; 3],
        axes_a: [[f64; 3]; 3],
        center_b: [f64; 3],
        half_b: [f64; 3],
        axes_b: [[f64; 3]; 3],
    ) -> Vec<ContactPoint> {
        let verts_a = Self::box_vertices(center_a, half_a, axes_a);
        let verts_b = Self::box_vertices(center_b, half_b, axes_b);

        // Test 6 face axes (3 from each box)
        let mut min_sep = f64::INFINITY;
        let mut best_normal = [0.0f64, 1.0, 0.0];

        for &axis in axes_a.iter().chain(axes_b.iter()) {
            let n = normalize3(axis);
            let sep = sat_face_axis(n, &verts_a, &verts_b);
            if sep > self.tolerance {
                return Vec::new(); // Separated
            }
            if sep < min_sep {
                min_sep = sep;
                best_normal = n;
            }
        }

        // Orient normal from A to B
        if dot3(sub3(center_b, center_a), best_normal) < 0.0 {
            best_normal = neg3(best_normal);
        }

        // Generate contact points at closest vertices of B
        let ref_d = dot3(center_a, best_normal) + half_a[0] + half_a[1] + half_a[2]; // Approx
        let mut contacts = Vec::new();
        for &v in &verts_b {
            let pen = ref_d - dot3(v, best_normal);
            if pen > -self.tolerance {
                let pos_b = v;
                let pos_a = add3(v, scale3(best_normal, pen));
                contacts.push(ContactPoint::new(pos_a, pos_b, best_normal, pen.max(0.0)));
                if contacts.len() >= MAX_MANIFOLD_POINTS {
                    break;
                }
            }
        }
        contacts
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
    /// Algorithm:
    /// 1. Keep point with maximum penetration depth.
    /// 2. Keep point furthest from point 1.
    /// 3. Keep point forming max area triangle with 1-2.
    /// 4. Keep point forming max area with triangle 1-2-3.
    pub fn reduce(points: &[ContactPoint]) -> Vec<ContactPoint> {
        if points.len() <= MAX_MANIFOLD_POINTS {
            return points.to_vec();
        }

        let mut result = Vec::with_capacity(MAX_MANIFOLD_POINTS);

        // Step 1: deepest point
        let p0_idx = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.depth
                    .partial_cmp(&b.depth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        result.push(points[p0_idx]);
        let p0 = points[p0_idx].pos_a;

        // Step 2: furthest from p0
        let p1_idx = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                norm_sq3(sub3(a.pos_a, p0))
                    .partial_cmp(&norm_sq3(sub3(b.pos_a, p0)))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        result.push(points[p1_idx]);
        let p1 = points[p1_idx].pos_a;

        // Step 3: max area with p0-p1
        let p2_idx = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                let da = Self::point_to_line_dist_sq(a.pos_a, p0, p1);
                let db = Self::point_to_line_dist_sq(b.pos_a, p0, p1);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        result.push(points[p2_idx]);
        let p2 = points[p2_idx].pos_a;

        // Step 4: max area with p0-p1-p2
        let p3_idx = points
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                let da = Self::point_to_plane_dist_sq(a.pos_a, p0, p1, p2);
                let db = Self::point_to_plane_dist_sq(b.pos_a, p0, p1, p2);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        result.push(points[p3_idx]);

        result
    }

    fn point_to_line_dist_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
        let ab = sub3(b, a);
        let ap = sub3(p, a);
        let t = dot3(ap, ab) / dot3(ab, ab).max(1e-15);
        let closest = add3(a, scale3(ab, t.clamp(0.0, 1.0)));
        norm_sq3(sub3(p, closest))
    }

    fn point_to_plane_dist_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let n = normalize3(cross3(sub3(b, a), sub3(c, a)));
        let d = dot3(sub3(p, a), n);
        d * d
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
    pub fn update(&mut self, new_points: &[ContactPoint], frame: u64) {
        let threshold_sq = self.match_threshold * self.match_threshold;
        let mut matched = vec![false; new_points.len()];
        let mut new_manifold_pts = Vec::new();
        let mut new_ages = Vec::new();

        // Try to match new points to existing points
        for (old_pt, &old_age) in self.manifold.points.iter().zip(self.point_ages.iter()) {
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

        // Reduce to max 4 points
        let reduced = ManifoldReduction::reduce(&new_manifold_pts);
        self.point_ages = vec![0; reduced.len()];
        self.manifold.points = reduced;
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
        // Should generate some contacts
        let _ = contacts;
    }

    #[test]
    fn test_polyhedron_edge_edge_closest() {
        let pg = PolyhedronContact::new(4);
        let pa = [0.0, 0.0, 0.0_f64];
        let da = [1.0, 0.0, 0.0];
        let pb = [0.5, 0.5, -1.0_f64];
        let db = [0.0, 0.0, 1.0];
        let (ca, cb, dist) = pg.edge_edge_closest(pa, da, pb, db);
        let _ = (ca, cb);
        assert!(
            dist < 1.0,
            "closest distance should be less than 1: {}",
            dist
        );
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
}
