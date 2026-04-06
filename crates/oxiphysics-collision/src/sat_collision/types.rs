//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::if_same_then_else, clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// Contact result from the arbitrary convex polyhedra SAT test.
#[derive(Clone, Debug)]
#[allow(dead_code)]
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
    /// Full OBB–OBB SAT test using 15 candidate axes.
    /// Returns `None` if separated, or `Some(SatContact)` with the minimum
    /// penetration axis.
    pub fn obb_obb_sat(obb_a: &Obb, obb_b: &Obb) -> Option<SatContact> {
        let mut min_depth = f64::MAX;
        let mut best_normal = [0.0f64; 3];
        let mut best_feature = ContactFeature::FaceFace {
            face_a: 0,
            face_b: 0,
        };
        let mut test = |axis: [f64; 3], feature: ContactFeature| -> bool {
            match Self::test_overlap_on_axis(obb_a, obb_b, axis) {
                None => false,
                Some(overlap) => {
                    if overlap < min_depth {
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
        for i in 0..3 {
            if !test(
                obb_a.axes[i],
                ContactFeature::FaceFace {
                    face_a: i * 2,
                    face_b: 0,
                },
            ) {
                return None;
            }
        }
        for j in 0..3 {
            if !test(
                obb_b.axes[j],
                ContactFeature::FaceFace {
                    face_a: 0,
                    face_b: j * 2,
                },
            ) {
                return None;
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let axis = cross(obb_a.axes[i], obb_b.axes[j]);
                if length(axis) > 1e-6
                    && !test(
                        axis,
                        ContactFeature::EdgeEdge {
                            edge_a: i,
                            edge_b: j,
                        },
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
    /// Uses the minimum-overlap axis to determine contact type, then
    /// extracts the relevant vertex/edge/face points.
    pub fn from_obb_obb(obb_a: &Obb, obb_b: &Obb, contact: &SatContact) -> Vec<[f64; 3]> {
        match &contact.contact_type {
            ContactFeature::FaceFace { .. } => {
                let support_b = obb_b.support(scale(contact.normal, -1.0));
                let ref_d = dot(obb_a.center, contact.normal) + contact.depth * 0.5;
                let proj = dot(support_b, contact.normal);
                let adjusted = add(support_b, scale(contact.normal, ref_d - proj));
                vec![adjusted]
            }
            ContactFeature::EdgeEdge { edge_a, edge_b } => {
                let dir_a = obb_a.axes[edge_a % 3];
                let dir_b = obb_b.axes[edge_b % 3];
                let p_a = obb_a.support(scale(contact.normal, -1.0));
                let p_b = obb_b.support(contact.normal);
                let (cp, _) = edge_edge_contact(p_a, dir_a, p_b, dir_b);
                vec![cp]
            }
            ContactFeature::VertexFace { .. } | ContactFeature::FaceEdge { .. } => {
                let v = obb_b.support(scale(contact.normal, -1.0));
                vec![v]
            }
        }
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
        for i in 0..3 {
            match test_axis(obb.axes[i]) {
                None => return None,
                Some(depth) if depth < min_depth => {
                    min_depth = depth;
                    best_normal = obb.axes[i];
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
        for i in 0..3 {
            for j in 0..3 {
                let axis = cross(obb.axes[i], tri_edges[j]);
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
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
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn new() -> Self {
        SatAxisCache {
            axis: [0.0, 1.0, 0.0],
            overlap: 0.0,
            valid: false,
        }
    }
    /// Invalidate the cache (e.g. when a body teleports or is removed).
    #[allow(dead_code)]
    pub fn invalidate(&mut self) {
        self.valid = false;
    }
    /// Update the cache with a new axis and overlap.
    #[allow(dead_code)]
    pub fn update(&mut self, axis: [f64; 3], overlap: f64) {
        self.axis = axis;
        self.overlap = overlap;
        self.valid = true;
    }
    /// Check the cached axis against two OBBs.
    ///
    /// Returns `Some(overlap)` if still overlapping in the cached direction,
    /// or `None` if the cache is invalid or the axis separates the OBBs.
    #[allow(dead_code)]
    pub fn check(&self, obb_a: &Obb, obb_b: &Obb) -> Option<f64> {
        if !self.valid {
            return None;
        }
        ObbCollision::test_overlap_on_axis(obb_a, obb_b, self.axis)
    }
}
/// A flat-array OBB BVH (distinct from the recursive `ObbTree`).
#[derive(Clone, Debug, Default)]
#[allow(dead_code)]
pub struct ObbBvh {
    /// All nodes (root = index 0).
    pub nodes: Vec<ObbBvhNode>,
}
impl ObbBvh {
    /// Create an empty BVH.
    #[allow(dead_code)]
    pub fn new() -> Self {
        ObbBvh { nodes: Vec::new() }
    }
    /// Build a balanced BVH from a slice of leaf OBBs.
    #[allow(dead_code)]
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
        if left_idx.is_empty() {
            left_idx = indices[..indices.len() / 2].to_vec();
            right_idx = indices[indices.len() / 2..].to_vec();
        } else if right_idx.is_empty() {
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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }
    /// Total number of nodes.
    #[allow(dead_code)]
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
    /// Test OBB vs Capsule for intersection.
    ///
    /// Returns `Some(depth, contact_point, normal)` if intersecting.
    pub fn test(obb: &Obb, capsule: &Capsule) -> Option<(f64, [f64; 3], [f64; 3])> {
        let closest_on_axis = capsule.closest_axis_point(obb.center);
        let closest_on_obb = Self::closest_point_on_obb(obb, closest_on_axis);
        let diff = sub(closest_on_axis, closest_on_obb);
        let dist = length(diff);
        if dist < capsule.radius {
            let depth = capsule.radius - dist;
            let normal = if dist > 1e-10 {
                scale(diff, 1.0 / dist)
            } else {
                let center_diff = sub(closest_on_axis, obb.center);
                let cl = length(center_diff);
                if cl > 1e-10 {
                    scale(center_diff, 1.0 / cl)
                } else {
                    [0.0, 1.0, 0.0]
                }
            };
            let contact_point = add(closest_on_obb, scale(normal, depth * 0.5));
            Some((depth, contact_point, normal))
        } else {
            None
        }
    }
}
/// SAT test between two arbitrary convex polyhedra.
#[allow(dead_code)]
pub struct ConvexPolyhedraSat;
impl ConvexPolyhedraSat {
    /// Test intersection between two convex polyhedra.
    ///
    /// Tests all face normals of both polyhedra plus all edge-cross-product axes.
    /// Returns `Some(contact)` if the polyhedra overlap, `None` if separated.
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn is_separated(a: &ConvexPolyhedron, b: &ConvexPolyhedron) -> bool {
        Self::test(a, b).is_none()
    }
}
/// A node in a flat-array OBB BVH.
#[derive(Clone, Debug)]
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn is_leaf(&self) -> bool {
        self.left == usize::MAX && self.right == usize::MAX
    }
}
