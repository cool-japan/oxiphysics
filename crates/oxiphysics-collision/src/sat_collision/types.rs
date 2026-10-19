//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Contact result from the arbitrary convex polyhedra SAT test.
#[derive(Clone, Debug)]
pub struct PolyhedraContact {
    /// Minimum-overlap separating axis (contact normal).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Feature type that produced the minimum overlap.
    pub feature: ContactFeature,
}
/// An oriented bounding box defined by a center, three orthogonal unit axes,
/// and per-axis half-extents.
#[derive(Clone, Debug)]
pub struct Obb {
    /// Center of the box in world space.
    pub center: [f64; 3],
    /// Three orthogonal unit axes (columns of the rotation matrix).
    pub axes: [[f64; 3]; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}
impl Obb {
    /// Construct an axis-aligned OBB from a center point and half-extents.
    pub fn axis_aligned(center: [f64; 3], half_extents: [f64; 3]) -> Self {
        Obb {
            center,
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            half_extents,
        }
    }
    /// Construct an OBB from an axis-aligned bounding box given by `min` and `max`.
    pub fn from_aabb(min: [f64; 3], max: [f64; 3]) -> Self {
        let center = [
            (min[0] + max[0]) * 0.5,
            (min[1] + max[1]) * 0.5,
            (min[2] + max[2]) * 0.5,
        ];
        let half_extents = [
            (max[0] - min[0]) * 0.5,
            (max[1] - min[1]) * 0.5,
            (max[2] - min[2]) * 0.5,
        ];
        Obb {
            center,
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            half_extents,
        }
    }
    /// Return all 8 corner vertices of the OBB in world space.
    pub fn vertices(&self) -> [[f64; 3]; 8] {
        let c = self.center;
        let a = self.axes;
        let h = self.half_extents;
        let mut verts = [[0.0f64; 3]; 8];
        for (i, s) in [
            [-1.0f64, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
        .iter()
        .enumerate()
        {
            verts[i] = add(
                add(add(c, scale(a[0], s[0] * h[0])), scale(a[1], s[1] * h[1])),
                scale(a[2], s[2] * h[2]),
            );
        }
        verts
    }
    /// Return the support point (furthest vertex) in direction `dir`.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let verts = self.vertices();
        let mut best = verts[0];
        let mut best_d = dot(best, dir);
        for &v in verts.iter().skip(1) {
            let d = dot(v, dir);
            if d > best_d {
                best_d = d;
                best = v;
            }
        }
        best
    }
    /// Project the OBB onto `axis` and return the (min, max) scalar interval.
    pub fn project_onto_axis(&self, axis: [f64; 3]) -> (f64, f64) {
        let center_proj = dot(self.center, axis);
        let radius = self.half_extents[0] * dot(self.axes[0], axis).abs()
            + self.half_extents[1] * dot(self.axes[1], axis).abs()
            + self.half_extents[2] * dot(self.axes[2], axis).abs();
        (center_proj - radius, center_proj + radius)
    }
}
/// Describes which geometric features are in contact.
#[derive(Clone, Debug, PartialEq)]
pub enum ContactFeature {
    /// Face-face contact.
    FaceFace { face_a: usize, face_b: usize },
    /// Face-edge contact.
    FaceEdge { face: usize, edge: usize },
    /// Edge-edge contact.
    EdgeEdge { edge_a: usize, edge_b: usize },
    /// Vertex-face contact.
    VertexFace { vertex: usize, face: usize },
}
/// Stateless namespace for OBB–OBB collision tests using SAT.
pub struct ObbCollision;
impl ObbCollision {
    /// Project both OBBs onto `axis`. Returns `None` if separated, or
    /// `Some(overlap)` if overlapping.
    pub fn test_overlap_on_axis(obb_a: &Obb, obb_b: &Obb, axis: [f64; 3]) -> Option<f64> {
        let len = length(axis);
        if len < 1e-12 {
            return None;
        }
        let norm_axis = scale(axis, 1.0 / len);
        let (min_a, max_a) = obb_a.project_onto_axis(norm_axis);
        let (min_b, max_b) = obb_b.project_onto_axis(norm_axis);
        let overlap = f64::min(max_a, max_b) - f64::max(min_a, min_b);
        if overlap < 0.0 { None } else { Some(overlap) }
    }
    /// Penetration of the two boxes along `axis`: the shortest translation of one
    /// box along the (normalised) axis after which the projections of the boxes
    /// no longer overlap. `None` if the axis separates the boxes (or has no
    /// direction).
    ///
    /// This is the quantity a minimum-translation search has to minimise. It
    /// equals the length of the common part of the two projected intervals only
    /// while neither interval contains the other: for a box that lies within the
    /// extent of the other along the axis the interval overlap is just the
    /// thickness of the inner box, whereas the translation needed to separate
    /// them is larger.
    pub fn penetration_on_axis(obb_a: &Obb, obb_b: &Obb, axis: [f64; 3]) -> Option<f64> {
        let len = length(axis);
        if len < 1e-12 {
            return None;
        }
        let norm_axis = scale(axis, 1.0 / len);
        let (min_a, max_a) = obb_a.project_onto_axis(norm_axis);
        let (min_b, max_b) = obb_b.project_onto_axis(norm_axis);
        let push = f64::min(max_a - min_b, max_b - min_a);
        if push < 0.0 { None } else { Some(push) }
    }
    /// Full OBB–OBB SAT test using 15 candidate axes.
    /// Returns `None` if separated, or `Some(SatContact)` with the minimum
    /// penetration axis: the axis (oriented from B toward A) with the shortest
    /// separating translation, see [`ObbCollision::penetration_on_axis`], and
    /// that translation as `depth`.
    pub fn obb_obb_sat(obb_a: &Obb, obb_b: &Obb) -> Option<SatContact> {
        let mut min_depth = f64::MAX;
        let mut best_normal = [0.0f64; 3];
        let mut best_feature = ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        };
        // A cross-product axis only wins when it beats the best face overlap by
        // more than rounding noise. Without the margin the cross product of two
        // parallel-face edges (equal to the face axis up to rounding) could win
        // and a flat face contact would be reported as a single edge contact.
        let extent = obb_a
            .half_extents
            .iter()
            .chain(obb_b.half_extents.iter())
            .fold(1.0_f64, |m, h| m.max(h.abs()));
        let edge_bias = 1e-6 * extent;
        let mut test = |axis: [f64; 3], feature: ContactFeature, bias: f64| -> bool {
            match Self::penetration_on_axis(obb_a, obb_b, axis) {
                None => false,
                Some(overlap) => {
                    if overlap < min_depth - bias {
                        min_depth = overlap;
                        let diff = sub(obb_a.center, obb_b.center);
                        let norm = normalize(axis);
                        best_normal = if dot(norm, diff) >= 0.0 {
                            norm
                        } else {
                            scale(norm, -1.0)
                        };
                        best_feature = feature;
                    }
                    true
                }
            }
        };
        for (i, &ax) in obb_a.axes.iter().enumerate() {
            if !test(
                ax,
                ContactFeature::FaceFace {
                    face_a: i * 2,
                    face_b: 0,
                },
                0.0,
            ) {
                return None;
            }
        }
        for (j, &bx) in obb_b.axes.iter().enumerate() {
            if !test(
                bx,
                ContactFeature::FaceFace {
                    face_a: 0,
                    face_b: j * 2,
                },
                0.0,
            ) {
                return None;
            }
        }
        for (i, &ai) in obb_a.axes.iter().enumerate() {
            for (j, &bj) in obb_b.axes.iter().enumerate() {
                let axis = cross(ai, bj);
                if length(axis) > 1e-6
                    && !test(
                        axis,
                        ContactFeature::EdgeEdge {
                            edge_a: i,
                            edge_b: j,
                        },
                        edge_bias,
                    )
                {
                    return None;
                }
            }
        }
        Some(SatContact {
            normal: best_normal,
            depth: min_depth,
            contact_type: best_feature,
        })
    }
    /// Fast boolean overlap test for two OBBs.
    pub fn box_box_fast(a: &Obb, b: &Obb) -> bool {
        Self::obb_obb_sat(a, b).is_some()
    }
}
/// Contact information produced by a SAT test.
#[derive(Clone, Debug)]
pub struct SatContact {
    /// Contact normal pointing from B towards A.
    pub normal: [f64; 3],
    /// Penetration depth (positive means overlap).
    pub depth: f64,
    /// The type of contact feature.
    pub contact_type: ContactFeature,
}
/// Generate contact points from SAT overlap information.
pub struct SatContactPointGenerator;
impl SatContactPointGenerator {
    /// Generate contact points from an OBB-OBB SAT result.
    ///
    /// `contact` is the result of [`ObbCollision::obb_obb_sat`] for the boxes
    /// `obb_a` and `obb_b` (normal from B toward A). The returned points lie
    /// half-way through the overlap of the two boxes:
    ///
    /// * **Face contact** (`FaceFace`, `FaceEdge`, `VertexFace`): the face of
    ///   the box whose face normal is best aligned with the contact normal is
    ///   the *reference face*; the most anti-parallel face of the other box is
    ///   the *incident face*. The incident quad is clipped against the four side
    ///   planes of the reference face (Sutherland-Hodgman, see
    ///   [`ContactManifoldBuilder::clip_polygon_by_plane`]) and the points of the
    ///   clipped polygon that lie below the reference plane are kept. Every
    ///   returned point is the middle of the segment between such an incident
    ///   point and its projection onto the reference plane, so the result is the
    ///   contact patch (typically four points, at most eight) halfway between
    ///   the two surfaces.
    /// * **Edge contact** (`EdgeEdge`): the two box edges (along the local axes
    ///   `edge_a` of A and `edge_b` of B) that face each other along the normal
    ///   are located and the middle of the segment between their closest points
    ///   is returned (one point).
    ///
    /// The contact normal is taken from `contact`, not recomputed, so the points
    /// are consistent with the axis the SAT test selected. The feature indices
    /// of a face contact are not used to pick the faces: the faces are derived
    /// from the normal, which also identifies the reference box when A and B
    /// both present a face along the axis.
    pub fn from_obb_obb(obb_a: &Obb, obb_b: &Obb, contact: &SatContact) -> Vec<[f64; 3]> {
        let normal = Self::unit_normal_toward_a(obb_a, obb_b, contact.normal);
        match &contact.contact_type {
            ContactFeature::EdgeEdge { edge_a, edge_b } => {
                Self::edge_edge_point(obb_a, obb_b, normal, edge_a % 3, edge_b % 3)
            }
            ContactFeature::FaceFace { .. }
            | ContactFeature::FaceEdge { .. }
            | ContactFeature::VertexFace { .. } => {
                let patch = Self::face_contact_patch(obb_a, obb_b, normal);
                if patch.is_empty() {
                    vec![Self::support_midpoint(obb_a, normal, contact.depth)]
                } else {
                    patch
                }
            }
        }
    }
    /// Unit contact normal pointing from B toward A; the direction between the
    /// centres (or +Y for coincident centres) when `normal` has no direction.
    fn unit_normal_toward_a(obb_a: &Obb, obb_b: &Obb, normal: [f64; 3]) -> [f64; 3] {
        let len = length(normal);
        if len > 1e-12 {
            return scale(normal, 1.0 / len);
        }
        let centres = sub(obb_a.center, obb_b.center);
        let centres_len = length(centres);
        if centres_len > 1e-12 {
            scale(centres, 1.0 / centres_len)
        } else {
            [0.0, 1.0, 0.0]
        }
    }
    /// The face of `obb` whose outward normal is best aligned with `outward`:
    /// `(local axis index, outward face normal, |cosine| of the alignment)`.
    fn best_face(obb: &Obb, outward: [f64; 3]) -> (usize, [f64; 3], f64) {
        let mut axis = 0;
        let mut alignment = -1.0;
        let mut sign = 1.0;
        for (k, ax) in obb.axes.iter().enumerate() {
            let d = dot(*ax, outward);
            if d.abs() > alignment {
                alignment = d.abs();
                axis = k;
                sign = if d >= 0.0 { 1.0 } else { -1.0 };
            }
        }
        (axis, scale(obb.axes[axis], sign), alignment)
    }
    /// Face-contact patch: clip the incident face against the reference face.
    /// Empty when nothing of the incident face lies over the reference face.
    fn face_contact_patch(obb_a: &Obb, obb_b: &Obb, normal: [f64; 3]) -> Vec<[f64; 3]> {
        // A face of A looks toward B along -normal, a face of B toward A along
        // +normal. The box that presents its face more squarely is the reference.
        let (axis_a, face_a, alignment_a) = Self::best_face(obb_a, scale(normal, -1.0));
        let (axis_b, face_b, alignment_b) = Self::best_face(obb_b, normal);
        let (reference, ref_axis, ref_normal, incident) = if alignment_a >= alignment_b - 1e-9 {
            (obb_a, axis_a, face_a, obb_b)
        } else {
            (obb_b, axis_b, face_b, obb_a)
        };
        let (inc_axis, inc_normal, _) = Self::best_face(incident, scale(ref_normal, -1.0));

        // Corners of the incident face, in cyclic order.
        let inc_u = (inc_axis + 1) % 3;
        let inc_v = (inc_axis + 2) % 3;
        let inc_centre = add(
            incident.center,
            scale(inc_normal, incident.half_extents[inc_axis]),
        );
        let du = scale(incident.axes[inc_u], incident.half_extents[inc_u]);
        let dv = scale(incident.axes[inc_v], incident.half_extents[inc_v]);
        let mut polygon = vec![
            add(add(inc_centre, du), dv),
            add(sub(inc_centre, du), dv),
            sub(sub(inc_centre, du), dv),
            sub(add(inc_centre, du), dv),
        ];

        // Clip against the four side planes of the reference face.
        for side in [(ref_axis + 1) % 3, (ref_axis + 2) % 3] {
            let axis = reference.axes[side];
            let offset = dot(reference.center, axis);
            let half = reference.half_extents[side];
            polygon = ContactManifoldBuilder::clip_polygon_by_plane(&polygon, axis, offset - half);
            polygon = ContactManifoldBuilder::clip_polygon_by_plane(
                &polygon,
                scale(axis, -1.0),
                -(offset + half),
            );
        }

        // Keep the points below the reference plane; each contact point is the
        // middle between the incident point and its projection onto the plane.
        let plane = dot(reference.center, ref_normal) + reference.half_extents[ref_axis];
        let largest_half = obb_a
            .half_extents
            .iter()
            .chain(obb_b.half_extents.iter())
            .fold(0.0_f64, |m, h| m.max(h.abs()));
        let slop = 1e-9 * (1.0 + largest_half);
        let mut points: Vec<[f64; 3]> = Vec::with_capacity(polygon.len());
        for p in polygon {
            let depth = plane - dot(p, ref_normal);
            if depth < -slop {
                continue;
            }
            let point = add(p, scale(ref_normal, 0.5 * depth.max(0.0)));
            let repeated = points
                .last()
                .is_some_and(|last| length(sub(*last, point)) <= slop);
            if !repeated {
                points.push(point);
            }
        }
        // The clip can repeat the first vertex at the end of the loop.
        if points.len() > 1
            && let (Some(first), Some(last)) = (points.first(), points.last())
            && length(sub(*first, *last)) <= slop
        {
            points.pop();
        }
        points
    }
    /// Middle between A's vertex that lies deepest toward B (the support point
    /// of A along `-normal`) and the point one `depth` further along `normal`.
    /// Used when a face contact leaves no clipped patch.
    fn support_midpoint(obb_a: &Obb, normal: [f64; 3], depth: f64) -> [f64; 3] {
        let deepest_a = obb_a.support(scale(normal, -1.0));
        add(deepest_a, scale(normal, 0.5 * depth))
    }
    /// The edge of `obb` parallel to local axis `axis` that lies farthest in
    /// direction `toward`, as `(start, direction)` with `direction` spanning the
    /// whole edge (length `2 * half_extent`).
    fn support_edge(obb: &Obb, axis: usize, toward: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let mut middle = obb.center;
        for (k, ax) in obb.axes.iter().enumerate() {
            if k == axis {
                continue;
            }
            let sign = if dot(*ax, toward) >= 0.0 { 1.0 } else { -1.0 };
            middle = add(middle, scale(*ax, sign * obb.half_extents[k]));
        }
        let half = scale(obb.axes[axis], obb.half_extents[axis]);
        (sub(middle, half), scale(half, 2.0))
    }
    /// Edge-edge contact point: the middle of the closest points of the real
    /// edges (A's along local axis `axis_a`, B's along `axis_b`) that face each
    /// other along `normal` (pointing from B toward A).
    fn edge_edge_point(
        obb_a: &Obb,
        obb_b: &Obb,
        normal: [f64; 3],
        axis_a: usize,
        axis_b: usize,
    ) -> Vec<[f64; 3]> {
        let (start_a, dir_a) = Self::support_edge(obb_a, axis_a, scale(normal, -1.0));
        let (start_b, dir_b) = Self::support_edge(obb_b, axis_b, normal);
        let (_, _, on_a, on_b) = closest_points_on_segments(start_a, dir_a, start_b, dir_b);
        vec![scale(add(on_a, on_b), 0.5)]
    }
}
/// SAT test between an OBB and a triangle.
pub struct ObbTriangleCollision;
impl ObbTriangleCollision {
    /// Test OBB vs triangle for intersection using SAT.
    ///
    /// Tests 3 OBB face axes, 1 triangle normal, and 9 edge-edge axes (3x3).
    pub fn test(obb: &Obb, tri: [[f64; 3]; 3]) -> Option<SatContact> {
        let mut min_depth = f64::MAX;
        let mut best_normal = [0.0f64; 3];
        let mut best_feature = ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        };
        let tri_edges = [
            sub(tri[1], tri[0]),
            sub(tri[2], tri[1]),
            sub(tri[0], tri[2]),
        ];
        let test_axis = |axis: [f64; 3]| -> Option<f64> {
            let len = length(axis);
            if len < 1e-12 {
                return Some(f64::MAX);
            }
            let norm = scale(axis, 1.0 / len);
            let (min_obb, max_obb) = obb.project_onto_axis(norm);
            let tri_proj: Vec<f64> = tri.iter().map(|&v| dot(v, norm)).collect();
            let min_tri = tri_proj.iter().cloned().fold(f64::MAX, f64::min);
            let max_tri = tri_proj.iter().cloned().fold(f64::MIN, f64::max);
            let overlap = f64::min(max_obb, max_tri) - f64::max(min_obb, min_tri);
            if overlap < 0.0 { None } else { Some(overlap) }
        };
        for (i, &ax) in obb.axes.iter().enumerate() {
            match test_axis(ax) {
                None => return None,
                Some(depth) if depth < min_depth => {
                    min_depth = depth;
                    best_normal = ax;
                    best_feature = ContactFeature::FaceFace {
                        face_a: i,
                        face_b: 0,
                    };
                }
                _ => {}
            }
        }
        let tri_normal = normalize(cross(tri_edges[0], tri_edges[2]));
        match test_axis(tri_normal) {
            None => return None,
            Some(depth) if depth < min_depth => {
                min_depth = depth;
                best_normal = tri_normal;
                best_feature = ContactFeature::FaceFace {
                    face_a: 0,
                    face_b: 1,
                };
            }
            _ => {}
        }
        for (i, &ai) in obb.axes.iter().enumerate() {
            for (j, &te) in tri_edges.iter().enumerate() {
                let axis = cross(ai, te);
                if length(axis) < 1e-9 {
                    continue;
                }
                match test_axis(axis) {
                    None => return None,
                    Some(depth) if depth < min_depth => {
                        min_depth = depth;
                        best_normal = normalize(axis);
                        best_feature = ContactFeature::EdgeEdge {
                            edge_a: i,
                            edge_b: j,
                        };
                    }
                    _ => {}
                }
            }
        }
        let tri_centroid = scale(add(add(tri[0], tri[1]), tri[2]), 1.0 / 3.0);
        let dir_to_tri = sub(tri_centroid, obb.center);
        if dot(best_normal, dir_to_tri) < 0.0 {
            best_normal = scale(best_normal, -1.0);
        }
        Some(SatContact {
            normal: best_normal,
            depth: min_depth,
            contact_type: best_feature,
        })
    }
}
/// A convex polytope stored as vertex, face, and edge lists.
#[derive(Clone, Debug)]
pub struct ConvexPolytope {
    /// Vertices in local space.
    pub vertices: Vec<[f64; 3]>,
    /// Faces as lists of vertex indices.
    pub faces: Vec<Vec<usize>>,
    /// Edges as pairs of vertex indices.
    pub edges: Vec<[usize; 2]>,
    /// Outward-pointing face normals (one per face).
    pub face_normals: Vec<[f64; 3]>,
}
impl ConvexPolytope {
    /// Decompose an OBB into a convex polytope (box with 6 faces and 12 edges).
    pub fn from_obb(obb: &Obb) -> Self {
        let verts = obb.vertices();
        let vertices: Vec<[f64; 3]> = verts.to_vec();
        let faces: Vec<Vec<usize>> = vec![
            vec![0, 2, 3, 1],
            vec![4, 5, 7, 6],
            vec![0, 1, 5, 4],
            vec![2, 6, 7, 3],
            vec![0, 4, 6, 2],
            vec![1, 3, 7, 5],
        ];
        let face_normals: Vec<[f64; 3]> = vec![
            scale(obb.axes[2], -1.0),
            obb.axes[2],
            scale(obb.axes[1], -1.0),
            obb.axes[1],
            scale(obb.axes[0], -1.0),
            obb.axes[0],
        ];
        let edges: Vec<[usize; 2]> = vec![
            [0, 1],
            [2, 3],
            [4, 5],
            [6, 7],
            [0, 2],
            [1, 3],
            [4, 6],
            [5, 7],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
        ];
        ConvexPolytope {
            vertices,
            faces,
            edges,
            face_normals,
        }
    }
    /// Return the furthest vertex in direction `dir` under transform `(translation, rotation)`.
    /// The rotation matrix columns are the local axes.
    pub fn support(&self, dir: [f64; 3], transform: ([f64; 3], [[f64; 3]; 3])) -> [f64; 3] {
        let (translation, axes) = transform;
        let local_dir = [dot(dir, axes[0]), dot(dir, axes[1]), dot(dir, axes[2])];
        let mut best_local = self.vertices[0];
        let mut best_d = dot(self.vertices[0], local_dir);
        for &v in self.vertices.iter().skip(1) {
            let d = dot(v, local_dir);
            if d > best_d {
                best_d = d;
                best_local = v;
            }
        }
        add(translation, rotate(axes, best_local))
    }
}
/// Stateless SAT collision test for general convex polytopes.
pub struct PolytopeCollision;
impl PolytopeCollision {
    /// SAT test between two convex polytopes each under their own transform.
    /// Tests all face normals of A, all face normals of B, and edge cross-products.
    pub fn sat_test(
        a: &ConvexPolytope,
        ta: ([f64; 3], [[f64; 3]; 3]),
        b: &ConvexPolytope,
        tb: ([f64; 3], [[f64; 3]; 3]),
    ) -> Option<SatContact> {
        let mut min_depth = f64::MAX;
        let mut best_normal = [0.0f64; 3];
        let mut best_feature = ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        };
        let project = |poly: &ConvexPolytope, t: ([f64; 3], [[f64; 3]; 3]), axis: [f64; 3]| {
            let (trans, axes) = t;
            let mut mn = f64::MAX;
            let mut mx = f64::MIN;
            for &lv in &poly.vertices {
                let wv = add(trans, rotate(axes, lv));
                let p = dot(wv, axis);
                if p < mn {
                    mn = p;
                }
                if p > mx {
                    mx = p;
                }
            }
            (mn, mx)
        };
        let mut test_axis = |axis: [f64; 3], feature: ContactFeature| -> bool {
            let len = length(axis);
            if len < 1e-9 {
                return true;
            }
            let norm = scale(axis, 1.0 / len);
            let (min_a, max_a) = project(a, ta, norm);
            let (min_b, max_b) = project(b, tb, norm);
            let overlap = f64::min(max_a, max_b) - f64::max(min_a, min_b);
            if overlap < 0.0 {
                return false;
            }
            if overlap < min_depth {
                min_depth = overlap;
                let diff = sub(ta.0, tb.0);
                best_normal = if dot(norm, diff) >= 0.0 {
                    norm
                } else {
                    scale(norm, -1.0)
                };
                best_feature = feature;
            }
            true
        };
        for (fi, &ln) in a.face_normals.iter().enumerate() {
            let wn = rotate(ta.1, ln);
            if !test_axis(
                wn,
                ContactFeature::FaceFace {
                    face_a: fi,
                    face_b: 0,
                },
            ) {
                return None;
            }
        }
        for (fi, &ln) in b.face_normals.iter().enumerate() {
            let wn = rotate(tb.1, ln);
            if !test_axis(
                wn,
                ContactFeature::FaceFace {
                    face_a: 0,
                    face_b: fi,
                },
            ) {
                return None;
            }
        }
        for (ei, &[ia0, ia1]) in a.edges.iter().enumerate() {
            let ea = normalize(sub(
                add(ta.0, rotate(ta.1, a.vertices[ia1])),
                add(ta.0, rotate(ta.1, a.vertices[ia0])),
            ));
            for (ej, &[ib0, ib1]) in b.edges.iter().enumerate() {
                let eb = normalize(sub(
                    add(tb.0, rotate(tb.1, b.vertices[ib1])),
                    add(tb.0, rotate(tb.1, b.vertices[ib0])),
                ));
                let axis = cross(ea, eb);
                if !test_axis(
                    axis,
                    ContactFeature::EdgeEdge {
                        edge_a: ei,
                        edge_b: ej,
                    },
                ) {
                    return None;
                }
            }
        }
        Some(SatContact {
            normal: best_normal,
            depth: min_depth,
            contact_type: best_feature,
        })
    }
}
/// A convex polyhedron for SAT tests, described by vertices and face normals.
#[derive(Clone, Debug)]
pub struct ConvexPolyhedron {
    /// Vertices in world space.
    pub vertices: Vec<[f64; 3]>,
    /// Face normals (outward-pointing, unit length).
    pub face_normals: Vec<[f64; 3]>,
    /// Edge directions (one per unique edge pair).
    pub edge_directions: Vec<[f64; 3]>,
}
impl ConvexPolyhedron {
    /// Create a convex polyhedron from vertices, face normals, and edge directions.
    pub fn new(
        vertices: Vec<[f64; 3]>,
        face_normals: Vec<[f64; 3]>,
        edge_directions: Vec<[f64; 3]>,
    ) -> Self {
        ConvexPolyhedron {
            vertices,
            face_normals,
            edge_directions,
        }
    }
    /// Support function: vertex furthest along `dir`.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        self.vertices
            .iter()
            .max_by(|a, b| {
                dot(**a, dir)
                    .partial_cmp(&dot(**b, dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or([0.0; 3])
    }
    /// Project all vertices onto `axis`, returning `(min, max)`.
    pub fn project(&self, axis: [f64; 3]) -> (f64, f64) {
        let mut mn = f64::INFINITY;
        let mut mx = f64::NEG_INFINITY;
        for &v in &self.vertices {
            let d = dot(v, axis);
            mn = mn.min(d);
            mx = mx.max(d);
        }
        (mn, mx)
    }
}
/// Cached separating axis from the previous frame.
///
/// Stores the best axis from the last SAT query so it can be checked first,
/// dramatically reducing the average number of axes tested when objects
/// move slowly relative to each other.
#[derive(Clone, Debug)]
pub struct SatAxisCache {
    /// The cached separating or minimum-penetration axis.
    pub axis: [f64; 3],
    /// The overlap (negative = separation) along the cached axis.
    pub overlap: f64,
    /// Whether the cache is valid for the current frame.
    pub valid: bool,
}
impl SatAxisCache {
    /// Create an empty (invalid) cache.
    pub fn new() -> Self {
        SatAxisCache {
            axis: [0.0, 1.0, 0.0],
            overlap: 0.0,
            valid: false,
        }
    }
    /// Invalidate the cache (e.g. when a body teleports or is removed).
    pub fn invalidate(&mut self) {
        self.valid = false;
    }
    /// Update the cache with a new axis and overlap.
    pub fn update(&mut self, axis: [f64; 3], overlap: f64) {
        self.axis = axis;
        self.overlap = overlap;
        self.valid = true;
    }
    /// Check the cached axis against two OBBs.
    ///
    /// Returns `Some(overlap)` if still overlapping in the cached direction,
    /// or `None` if the cache is invalid or the axis separates the OBBs.
    pub fn check(&self, obb_a: &Obb, obb_b: &Obb) -> Option<f64> {
        if !self.valid {
            return None;
        }
        ObbCollision::test_overlap_on_axis(obb_a, obb_b, self.axis)
    }
}
/// A flat-array OBB BVH (distinct from the recursive `ObbTree`).
#[derive(Clone, Debug, Default)]
pub struct ObbBvh {
    /// All nodes (root = index 0).
    pub nodes: Vec<ObbBvhNode>,
}
impl ObbBvh {
    /// Create an empty BVH.
    pub fn new() -> Self {
        ObbBvh { nodes: Vec::new() }
    }
    /// Build a balanced BVH from a slice of leaf OBBs.
    pub fn build(leaf_obbs: &[Obb]) -> Self {
        if leaf_obbs.is_empty() {
            return Self::new();
        }
        let mut bvh = ObbBvh { nodes: Vec::new() };
        let indices: Vec<usize> = (0..leaf_obbs.len()).collect();
        bvh.build_recursive(leaf_obbs, &indices);
        bvh
    }
    fn build_recursive(&mut self, obbs: &[Obb], indices: &[usize]) -> usize {
        let node_idx = self.nodes.len();
        if indices.len() == 1 {
            self.nodes.push(ObbBvhNode {
                bounds: obbs[indices[0]].clone(),
                left: usize::MAX,
                right: usize::MAX,
                geometry_index: Some(indices[0]),
            });
            return node_idx;
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for &i in indices {
            let verts = obbs[i].vertices();
            for v in &verts {
                for k in 0..3 {
                    mn[k] = mn[k].min(v[k]);
                    mx[k] = mx[k].max(v[k]);
                }
            }
        }
        let merged = Obb::from_aabb(mn, mx);
        let extents = [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]];
        let split_axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
            0
        } else if extents[1] >= extents[2] {
            1
        } else {
            2
        };
        let split_val = (mn[split_axis] + mx[split_axis]) * 0.5;
        let (mut left_idx, mut right_idx): (Vec<usize>, Vec<usize>) = indices
            .iter()
            .partition(|&&i| obbs[i].center[split_axis] <= split_val);
        if left_idx.is_empty() || right_idx.is_empty() {
            left_idx = indices[..indices.len() / 2].to_vec();
            right_idx = indices[indices.len() / 2..].to_vec();
        }
        self.nodes.push(ObbBvhNode {
            bounds: merged,
            left: usize::MAX,
            right: usize::MAX,
            geometry_index: None,
        });
        let left_child = self.build_recursive(obbs, &left_idx);
        let right_child = self.build_recursive(obbs, &right_idx);
        self.nodes[node_idx].left = left_child;
        self.nodes[node_idx].right = right_child;
        node_idx
    }
    /// Query all leaf geometry indices whose bounds overlap `query_obb`.
    pub fn query_overlapping(&self, query_obb: &Obb) -> Vec<usize> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut results = Vec::new();
        self.query_recursive(query_obb, 0, &mut results);
        results
    }
    fn query_recursive(&self, query: &Obb, node_idx: usize, results: &mut Vec<usize>) {
        let node = &self.nodes[node_idx];
        if ObbCollision::obb_obb_sat(&node.bounds, query).is_none() {
            return;
        }
        if node.is_leaf() {
            if let Some(g) = node.geometry_index {
                results.push(g);
            }
        } else {
            if node.left < self.nodes.len() {
                self.query_recursive(query, node.left, results);
            }
            if node.right < self.nodes.len() {
                self.query_recursive(query, node.right, results);
            }
        }
    }
    /// Number of leaf nodes.
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }
    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}
/// A capsule defined by two endpoint centers and a radius.
#[derive(Clone, Debug)]
pub struct Capsule {
    /// Center of the first end-cap.
    pub p0: [f64; 3],
    /// Center of the second end-cap.
    pub p1: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}
impl Capsule {
    /// Create a new capsule.
    pub fn new(p0: [f64; 3], p1: [f64; 3], radius: f64) -> Self {
        Self { p0, p1, radius }
    }
    /// Axis direction (p1 - p0), not normalised.
    pub fn axis(&self) -> [f64; 3] {
        sub(self.p1, self.p0)
    }
    /// Capsule half-length (excluding radius).
    pub fn half_length(&self) -> f64 {
        length(self.axis()) * 0.5
    }
    /// Closest point on the capsule axis segment to a query point.
    pub fn closest_axis_point(&self, q: [f64; 3]) -> [f64; 3] {
        let ax = self.axis();
        let ax_len_sq = dot(ax, ax);
        if ax_len_sq < 1e-24 {
            return self.p0;
        }
        let t = dot(sub(q, self.p0), ax) / ax_len_sq;
        let t_clamped = t.clamp(0.0, 1.0);
        add(self.p0, scale(ax, t_clamped))
    }
    /// Support function in direction `dir`.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        let d = dot(self.axis(), dir);
        let axis_pt = if d >= 0.0 { self.p1 } else { self.p0 };
        let dir_len = length(dir);
        if dir_len < 1e-12 {
            return axis_pt;
        }
        add(axis_pt, scale(dir, self.radius / dir_len))
    }
}
/// A bounding volume hierarchy of OBBs for triangle meshes.
#[derive(Debug)]
pub struct ObbTree {
    /// The OBB for this node.
    pub obb: Obb,
    /// Left child.
    pub left: Option<Box<ObbTree>>,
    /// Right child.
    pub right: Option<Box<ObbTree>>,
    /// Triangle indices stored at this leaf (empty for internal nodes).
    pub triangle_indices: Vec<usize>,
}
impl ObbTree {
    /// Build an OBB tree from a triangle mesh.
    pub fn build(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Self {
        let indices: Vec<usize> = (0..triangles.len()).collect();
        Self::build_recursive(vertices, triangles, &indices)
    }
    fn triangle_centroid(vertices: &[[f64; 3]], tri: [usize; 3]) -> [f64; 3] {
        let (a, b, c) = (vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
        scale(add(add(a, b), c), 1.0 / 3.0)
    }
    fn compute_obb_for_triangles(
        vertices: &[[f64; 3]],
        triangles: &[[usize; 3]],
        indices: &[usize],
    ) -> Obb {
        let mut pts: Vec<[f64; 3]> = Vec::new();
        for &ti in indices {
            for &vi in &triangles[ti] {
                pts.push(vertices[vi]);
            }
        }
        if pts.is_empty() {
            return Obb {
                center: [0.0; 3],
                axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                half_extents: [0.0; 3],
            };
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for &p in &pts {
            for k in 0..3 {
                if p[k] < mn[k] {
                    mn[k] = p[k];
                }
                if p[k] > mx[k] {
                    mx[k] = p[k];
                }
            }
        }
        Obb::from_aabb(mn, mx)
    }
    fn build_recursive(vertices: &[[f64; 3]], triangles: &[[usize; 3]], indices: &[usize]) -> Self {
        let obb = Self::compute_obb_for_triangles(vertices, triangles, indices);
        if indices.len() <= 1 {
            return ObbTree {
                obb,
                left: None,
                right: None,
                triangle_indices: indices.to_vec(),
            };
        }
        let longest_axis = {
            let he = obb.half_extents;
            if he[0] >= he[1] && he[0] >= he[2] {
                0
            } else if he[1] >= he[2] {
                1
            } else {
                2
            }
        };
        let mut sorted_indices = indices.to_vec();
        sorted_indices.sort_by(|&ia, &ib| {
            let ca = Self::triangle_centroid(vertices, triangles[ia]);
            let cb = Self::triangle_centroid(vertices, triangles[ib]);
            ca[longest_axis]
                .partial_cmp(&cb[longest_axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted_indices.len() / 2;
        let left_indices = &sorted_indices[..mid];
        let right_indices = &sorted_indices[mid..];
        let left = Box::new(Self::build_recursive(vertices, triangles, left_indices));
        let right = Box::new(Self::build_recursive(vertices, triangles, right_indices));
        ObbTree {
            obb,
            left: Some(left),
            right: Some(right),
            triangle_indices: Vec::new(),
        }
    }
    /// Query overlapping triangle pairs between two trees.
    pub fn query_overlap(&self, other: &ObbTree) -> Vec<(usize, usize)> {
        let mut result = Vec::new();
        self.query_recursive(other, &mut result);
        result
    }
    fn query_recursive(&self, other: &ObbTree, result: &mut Vec<(usize, usize)>) {
        if ObbCollision::box_box_fast(&self.obb, &other.obb) {
            match (&self.left, &self.right, &other.left, &other.right) {
                (None, None, None, None) => {
                    for &ta in &self.triangle_indices {
                        for &tb in &other.triangle_indices {
                            result.push((ta, tb));
                        }
                    }
                }
                _ => {
                    let self_children: Vec<&ObbTree> = [&self.left, &self.right]
                        .iter()
                        .filter_map(|c| c.as_deref())
                        .collect();
                    let other_children: Vec<&ObbTree> = [&other.left, &other.right]
                        .iter()
                        .filter_map(|c| c.as_deref())
                        .collect();
                    if self_children.is_empty() {
                        for oc in &other_children {
                            self.query_recursive(oc, result);
                        }
                    } else {
                        for sc in &self_children {
                            sc.query_recursive(other, result);
                        }
                    }
                }
            }
        }
    }
}
/// Builds contact manifolds via polygon clipping.
pub struct ContactManifoldBuilder;
impl ContactManifoldBuilder {
    /// Clip `polygon` by the half-space `dot(p, plane_normal) >= plane_d`.
    /// Returns the clipped polygon vertices.
    pub fn clip_polygon_by_plane(
        polygon: &[[f64; 3]],
        plane_normal: [f64; 3],
        plane_d: f64,
    ) -> Vec<[f64; 3]> {
        if polygon.is_empty() {
            return Vec::new();
        }
        let mut output = Vec::new();
        let n = polygon.len();
        for i in 0..n {
            let a = polygon[i];
            let b = polygon[(i + 1) % n];
            let da = dot(a, plane_normal) - plane_d;
            let db = dot(b, plane_normal) - plane_d;
            if da >= 0.0 {
                output.push(a);
            }
            if (da >= 0.0) != (db >= 0.0) {
                let t = da / (da - db);
                output.push(add(a, scale(sub(b, a), t)));
            }
        }
        output
    }
    /// Build a contact manifold by clipping the incident face polygon against
    /// the planes bounding the reference face (Sutherland-Hodgman algorithm).
    pub fn build_face_face(
        poly_a: &ConvexPolytope,
        face_a_idx: usize,
        ta: ([f64; 3], [[f64; 3]; 3]),
        poly_b: &ConvexPolytope,
        face_b_idx: usize,
        tb: ([f64; 3], [[f64; 3]; 3]),
    ) -> Vec<[f64; 3]> {
        let face_a = &poly_a.faces[face_a_idx];
        let ref_verts: Vec<[f64; 3]> = face_a
            .iter()
            .map(|&vi| add(ta.0, rotate(ta.1, poly_a.vertices[vi])))
            .collect();
        let face_b = &poly_b.faces[face_b_idx];
        let mut clipped: Vec<[f64; 3]> = face_b
            .iter()
            .map(|&vi| add(tb.0, rotate(tb.1, poly_b.vertices[vi])))
            .collect();
        let nref = ref_verts.len();
        let ref_normal = rotate(ta.1, poly_a.face_normals[face_a_idx]);
        for i in 0..nref {
            if clipped.is_empty() {
                break;
            }
            let edge_start = ref_verts[i];
            let edge_end = ref_verts[(i + 1) % nref];
            let edge_dir = sub(edge_end, edge_start);
            let side_normal = normalize(cross(ref_normal, edge_dir));
            let plane_d = dot(side_normal, edge_start);
            clipped = Self::clip_polygon_by_plane(&clipped, side_normal, plane_d);
        }
        let ref_plane_d = dot(ref_normal, ref_verts[0]);
        clipped.retain(|p| dot(*p, ref_normal) <= ref_plane_d + 1e-9);
        clipped
    }
}
/// SAT test between an OBB and a capsule.
pub struct ObbCapsuleCollision;
impl ObbCapsuleCollision {
    /// Relative distance below which the capsule axis counts as touching or
    /// crossing the box, so that the contact normal cannot come from the
    /// (undefined) direction of a zero distance.
    const TOUCH_TOLERANCE: f64 = 1e-12;
    /// Compute the closest point on the OBB surface to a query point.
    pub fn closest_point_on_obb(obb: &Obb, q: [f64; 3]) -> [f64; 3] {
        let d = sub(q, obb.center);
        let mut result = obb.center;
        for i in 0..3 {
            let dist = dot(d, obb.axes[i]);
            let clamped = dist.clamp(-obb.half_extents[i], obb.half_extents[i]);
            result = add(result, scale(obb.axes[i], clamped));
        }
        result
    }
    /// Closest pair of points between the segment `p0 -> p1` and the solid OBB.
    ///
    /// Returns `(point_on_segment, point_on_box)`.
    ///
    /// * When the segment touches or crosses the box the two points coincide:
    ///   they are the middle of the part of the segment inside the box, found by
    ///   clipping the segment against the three slabs of the box.
    /// * Otherwise, in the frame of the box, the squared distance from the point
    ///   `p(s) = p0 + s (p1 - p0)` to the box is `sum_k e_k(s)^2`, where `e_k`
    ///   is the excess of coordinate `k` beyond `+-half_extent_k`. Each `e_k` is
    ///   zero inside the slab and linear in `s` outside it, so the distance is a
    ///   convex piecewise quadratic of `s` whose pieces change at the (at most
    ///   six) parameters where a coordinate crosses a slab boundary. The minimum
    ///   of every piece is found in closed form, which makes the result exact;
    ///   no iteration or tolerance is involved.
    pub fn closest_points_segment_obb(
        obb: &Obb,
        p0: [f64; 3],
        p1: [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let to_local = |p: [f64; 3]| {
            let rel = sub(p, obb.center);
            [
                dot(rel, obb.axes[0]),
                dot(rel, obb.axes[1]),
                dot(rel, obb.axes[2]),
            ]
        };
        let half = obb.half_extents;
        let start = to_local(p0);
        let dir = sub(to_local(p1), start);
        let clamp_to_box = |p: [f64; 3]| {
            [
                p[0].max(-half[0]).min(half[0]),
                p[1].max(-half[1]).min(half[1]),
                p[2].max(-half[2]).min(half[2]),
            ]
        };
        let distance_sq = |s: f64| {
            let p = add(start, scale(dir, s));
            let nearest = clamp_to_box(p);
            let excess = sub(p, nearest);
            dot(excess, excess)
        };
        let to_world_pair = |s: f64| {
            let on_segment = add(p0, scale(sub(p1, p0), s));
            let nearest = clamp_to_box(add(start, scale(dir, s)));
            let mut on_box = obb.center;
            for (axis, coordinate) in obb.axes.iter().zip(nearest.iter()) {
                on_box = add(on_box, scale(*axis, *coordinate));
            }
            (on_segment, on_box)
        };

        // The part of the segment inside the box: clip it against the slabs.
        let mut inside_from = 0.0_f64;
        let mut inside_to = 1.0_f64;
        let mut within_slabs = true;
        for k in 0..3 {
            if dir[k] == 0.0 {
                within_slabs &= start[k].abs() <= half[k];
            } else {
                let a = (-half[k] - start[k]) / dir[k];
                let b = (half[k] - start[k]) / dir[k];
                inside_from = inside_from.max(a.min(b));
                inside_to = inside_to.min(a.max(b));
            }
        }
        if within_slabs && inside_from <= inside_to {
            return to_world_pair(0.5 * (inside_from + inside_to));
        }

        // Parameters at which a coordinate crosses a slab boundary.
        let mut cuts = vec![0.0_f64, 1.0];
        for k in 0..3 {
            if dir[k] != 0.0 {
                for bound in [-half[k], half[k]] {
                    let s = (bound - start[k]) / dir[k];
                    if s > 0.0 && s < 1.0 {
                        cuts.push(s);
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);

        let mut best_s = 0.0;
        let mut best = f64::INFINITY;
        for pair in cuts.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            if hi <= lo {
                continue;
            }
            // On (lo, hi) every coordinate is below, inside or above its slab,
            // so the squared distance is `quad * s^2 + 2 * lin * s + const`.
            let mid = 0.5 * (lo + hi);
            let mut quad = 0.0;
            let mut lin = 0.0;
            for k in 0..3 {
                let at_mid = start[k] + mid * dir[k];
                let offset = if at_mid > half[k] {
                    start[k] - half[k]
                } else if at_mid < -half[k] {
                    start[k] + half[k]
                } else {
                    continue;
                };
                quad += dir[k] * dir[k];
                lin += offset * dir[k];
            }
            let s = if quad > 0.0 {
                (-lin / quad).max(lo).min(hi)
            } else {
                mid
            };
            let value = distance_sq(s);
            if value < best {
                best = value;
                best_s = s;
            }
        }

        to_world_pair(best_s)
    }
    /// Minimum-translation axis and depth for a capsule whose axis crosses the
    /// box, as `(normal from the box toward the capsule, depth)`.
    ///
    /// The Minkowski difference of the box and the capsule axis is a convex
    /// polytope whose faces have the normals of the three box axes and of the
    /// cross products of each box axis with the capsule axis; the capsule radius
    /// then only inflates it. The shortest translation that separates the shapes
    /// is therefore found among those six axes: along a unit axis `L` it is
    /// `r_box + r_capsule - |c . L|`, with `r_box = sum_k h_k |u_k . L|`,
    /// `r_capsule = radius + |d . L| / 2` and `c` the offset between the
    /// capsule midpoint and the box centre.
    fn crossing_axis(obb: &Obb, capsule: &Capsule) -> ([f64; 3], f64) {
        let dir = capsule.axis();
        let dir_len = length(dir);
        let midpoint = scale(add(capsule.p0, capsule.p1), 0.5);
        let offset = sub(midpoint, obb.center);
        let mut best_overlap = f64::INFINITY;
        let mut best_axis = [0.0, 1.0, 0.0];
        let mut consider = |axis: [f64; 3]| {
            let box_radius: f64 = (0..3)
                .map(|k| obb.half_extents[k] * dot(obb.axes[k], axis).abs())
                .sum();
            let capsule_radius = capsule.radius + 0.5 * dot(dir, axis).abs();
            let separation = dot(offset, axis);
            let overlap = box_radius + capsule_radius - separation.abs();
            if overlap < best_overlap {
                best_overlap = overlap;
                best_axis = if separation >= 0.0 {
                    axis
                } else {
                    scale(axis, -1.0)
                };
            }
        };
        for axis in obb.axes {
            consider(axis);
        }
        for axis in obb.axes {
            let crossed = cross(axis, dir);
            let crossed_len = length(crossed);
            if crossed_len > 1e-9 * dir_len {
                consider(scale(crossed, 1.0 / crossed_len));
            }
        }
        (best_axis, best_overlap.max(0.0))
    }
    /// Test OBB vs Capsule for intersection.
    ///
    /// Returns `Some((depth, contact_point, normal))` if the capsule penetrates
    /// the box, `None` otherwise.
    ///
    /// * `normal` is the unit direction from the OBB toward the capsule; moving
    ///   the capsule by `depth` along it just separates the shapes.
    /// * `contact_point` lies in both shapes. When the capsule axis stays outside
    ///   the box it is half-way between the box point closest to the axis and
    ///   the deepest point of the capsule; when the axis crosses the box it is
    ///   the middle of the part of the axis inside the box.
    ///
    /// The capsule axis is a segment, so the shapes overlap exactly when the
    /// segment is closer to the box than the radius. The closest pair of the
    /// segment and the box is computed exactly
    /// ([`ObbCapsuleCollision::closest_points_segment_obb`]); when the axis stays
    /// outside the box, `depth = radius - distance` and the normal is the
    /// direction of that pair. When the axis touches or crosses the box the
    /// closest pair has no direction, and the minimum-translation axis is taken
    /// from the box axes and their cross products with the capsule axis.
    pub fn test(obb: &Obb, capsule: &Capsule) -> Option<(f64, [f64; 3], [f64; 3])> {
        let (on_axis, on_obb) = Self::closest_points_segment_obb(obb, capsule.p0, capsule.p1);
        let diff = sub(on_axis, on_obb);
        let dist = length(diff);
        if dist >= capsule.radius {
            return None;
        }
        let extent = obb.half_extents.iter().fold(0.0_f64, |m, h| m.max(h.abs()));
        if dist > Self::TOUCH_TOLERANCE * (1.0 + extent) {
            let normal = scale(diff, 1.0 / dist);
            let depth = capsule.radius - dist;
            let contact_point = sub(on_obb, scale(normal, 0.5 * depth));
            return Some((depth, contact_point, normal));
        }
        // The closest pair is the middle of the part of the axis inside the box:
        // a point of the capsule that lies in the box.
        let (normal, depth) = Self::crossing_axis(obb, capsule);
        Some((depth, on_axis, normal))
    }
}
/// SAT test between two arbitrary convex polyhedra.
pub struct ConvexPolyhedraSat;
impl ConvexPolyhedraSat {
    /// Test intersection between two convex polyhedra.
    ///
    /// Tests all face normals of both polyhedra plus all edge-cross-product axes.
    /// Returns `Some(contact)` if the polyhedra overlap, `None` if separated.
    pub fn test(a: &ConvexPolyhedron, b: &ConvexPolyhedron) -> Option<PolyhedraContact> {
        let mut min_depth = f64::INFINITY;
        let mut best_normal = [0.0_f64, 1.0, 0.0];
        let mut best_feature = ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        };
        for (i, &n) in a.face_normals.iter().enumerate() {
            let len = length(n);
            if len < 1e-12 {
                continue;
            }
            let axis = scale(n, 1.0 / len);
            let (a_min, a_max) = a.project(axis);
            let (b_min, b_max) = b.project(axis);
            let overlap = a_max.min(b_max) - a_min.max(b_min);
            if overlap < 0.0 {
                return None;
            }
            if overlap < min_depth {
                min_depth = overlap;
                best_normal = axis;
                best_feature = ContactFeature::FaceFace {
                    face_a: i,
                    face_b: 0,
                };
            }
        }
        for (j, &n) in b.face_normals.iter().enumerate() {
            let len = length(n);
            if len < 1e-12 {
                continue;
            }
            let axis = scale(n, 1.0 / len);
            let (a_min, a_max) = a.project(axis);
            let (b_min, b_max) = b.project(axis);
            let overlap = a_max.min(b_max) - a_min.max(b_min);
            if overlap < 0.0 {
                return None;
            }
            if overlap < min_depth {
                min_depth = overlap;
                best_normal = axis;
                best_feature = ContactFeature::FaceFace {
                    face_a: 0,
                    face_b: j,
                };
            }
        }
        for (ia, &ea) in a.edge_directions.iter().enumerate() {
            for (ib, &eb) in b.edge_directions.iter().enumerate() {
                let axis = cross(ea, eb);
                let len = length(axis);
                if len < 1e-12 {
                    continue;
                }
                let axis = scale(axis, 1.0 / len);
                let (a_min, a_max) = a.project(axis);
                let (b_min, b_max) = b.project(axis);
                let overlap = a_max.min(b_max) - a_min.max(b_min);
                if overlap < 0.0 {
                    return None;
                }
                if overlap < min_depth {
                    min_depth = overlap;
                    best_normal = axis;
                    best_feature = ContactFeature::EdgeEdge {
                        edge_a: ia,
                        edge_b: ib,
                    };
                }
            }
        }
        Some(PolyhedraContact {
            normal: best_normal,
            depth: min_depth,
            feature: best_feature,
        })
    }
    /// Returns `true` if the two convex polyhedra are separated.
    pub fn is_separated(a: &ConvexPolyhedron, b: &ConvexPolyhedron) -> bool {
        Self::test(a, b).is_none()
    }
}
/// A node in a flat-array OBB BVH.
#[derive(Clone, Debug)]
pub struct ObbBvhNode {
    /// Bounding OBB for this node.
    pub bounds: Obb,
    /// Left child index (`usize::MAX` = no child / leaf).
    pub left: usize,
    /// Right child index (`usize::MAX` = no child / leaf).
    pub right: usize,
    /// For leaf nodes: geometry index.
    pub geometry_index: Option<usize>,
}
impl ObbBvhNode {
    /// Returns `true` if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left == usize::MAX && self.right == usize::MAX
    }
}
