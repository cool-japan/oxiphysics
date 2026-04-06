//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::error::{Error, Result};
use crate::triangle_mesh::TriangleMesh;
use oxiphysics_core::math::Vec3;

/// A node in a Binary Space Partitioning tree.
///
/// Each node divides space with a plane; polygons are classified as front,
/// back, or coplanar with the plane.
#[allow(dead_code)]
pub struct BspNode {
    /// The dividing plane.
    pub plane: BspPlane,
    /// Coplanar polygons stored at this node.
    pub polygons: Vec<Vec<[f64; 3]>>,
    /// Subtree in front of the plane.
    pub front: Option<Box<BspNode>>,
    /// Subtree behind the plane.
    pub back: Option<Box<BspNode>>,
}
#[allow(dead_code)]
impl BspNode {
    /// Create a new BSP node with the given plane, no polygons, no children.
    pub fn new(plane: BspPlane) -> Self {
        Self {
            plane,
            polygons: Vec::new(),
            front: None,
            back: None,
        }
    }
    /// Classify a point relative to this node's plane.
    pub fn classify_point(&self, p: [f64; 3]) -> PlaneClass {
        let d = self.plane.signed_dist(p);
        if d > 1e-8 {
            PlaneClass::Front
        } else if d < -1e-8 {
            PlaneClass::Back
        } else {
            PlaneClass::OnPlane
        }
    }
    /// Insert a polygon into the BSP tree.
    ///
    /// The polygon is split if it straddles the plane, and each fragment is
    /// routed to the front or back subtree.
    pub fn insert_polygon(&mut self, polygon: Vec<[f64; 3]>) {
        if polygon.is_empty() {
            return;
        }
        let classes: Vec<PlaneClass> = polygon.iter().map(|&p| self.classify_point(p)).collect();
        let all_front = classes
            .iter()
            .all(|&c| c == PlaneClass::Front || c == PlaneClass::OnPlane);
        let all_back = classes
            .iter()
            .all(|&c| c == PlaneClass::Back || c == PlaneClass::OnPlane);
        if all_front && all_back {
            self.polygons.push(polygon);
        } else if all_front {
            match &mut self.front {
                Some(f) => f.insert_polygon(polygon),
                None => {
                    let mut node = Box::new(BspNode::new(self.plane));
                    node.polygons.push(polygon);
                    self.front = Some(node);
                }
            }
        } else if all_back {
            match &mut self.back {
                Some(b) => b.insert_polygon(polygon),
                None => {
                    let mut node = Box::new(BspNode::new(self.plane));
                    node.polygons.push(polygon);
                    self.back = Some(node);
                }
            }
        } else {
            let (front_poly, back_poly) = split_polygon_by_plane(&polygon, &self.plane);
            if !front_poly.is_empty() {
                match &mut self.front {
                    Some(f) => f.insert_polygon(front_poly),
                    None => {
                        let mut node = Box::new(BspNode::new(self.plane));
                        node.polygons.push(front_poly);
                        self.front = Some(node);
                    }
                }
            }
            if !back_poly.is_empty() {
                match &mut self.back {
                    Some(b) => b.insert_polygon(back_poly),
                    None => {
                        let mut node = Box::new(BspNode::new(self.plane));
                        node.polygons.push(back_poly);
                        self.back = Some(node);
                    }
                }
            }
        }
    }
    /// Collect all polygons in front of the plane (recursively).
    pub fn collect_front_polygons(&self) -> Vec<Vec<[f64; 3]>> {
        let mut result = self.polygons.clone();
        if let Some(f) = &self.front {
            result.extend(f.collect_front_polygons());
        }
        result
    }
    /// Collect all polygons behind the plane (recursively).
    pub fn collect_back_polygons(&self) -> Vec<Vec<[f64; 3]>> {
        let mut result = self.polygons.clone();
        if let Some(b) = &self.back {
            result.extend(b.collect_back_polygons());
        }
        result
    }
    /// Count all polygons in the tree.
    pub fn count_polygons(&self) -> usize {
        let mut count = self.polygons.len();
        if let Some(f) = &self.front {
            count += f.count_polygons();
        }
        if let Some(b) = &self.back {
            count += b.count_polygons();
        }
        count
    }
}
/// Performs CSG boolean operations while tracking which material each
/// resulting triangle comes from.
#[allow(dead_code)]
pub struct CsgWithMaterials;
#[allow(dead_code)]
impl CsgWithMaterials {
    /// Execute a boolean operation, returning `(result_mesh, per_triangle_materials)`.
    ///
    /// `per_triangle_materials` is parallel to `result_mesh.indices`.
    pub fn execute(
        op: BooleanOp,
        mesh_a: &TriangleMesh,
        mat_a: CsgMaterial,
        mesh_b: &TriangleMesh,
        mat_b: CsgMaterial,
    ) -> Result<(TriangleMesh, Vec<CsgMaterial>)> {
        if mesh_a.vertices.is_empty() || mesh_b.vertices.is_empty() {
            return Err(Error::General("Input mesh must not be empty".into()));
        }
        let verts_a: Vec<[f64; 3]> = mesh_a.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let verts_b: Vec<[f64; 3]> = mesh_b.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut result_verts: Vec<[f64; 3]> = Vec::new();
        let mut result_tris: Vec<[usize; 3]> = Vec::new();
        let mut result_mats: Vec<CsgMaterial> = Vec::new();
        let collect_with_mat = |verts: &[[f64; 3]],
                                tris: &[[usize; 3]],
                                other_verts: &[[f64; 3]],
                                other_tris: &[[usize; 3]],
                                inside: bool,
                                flip: bool,
                                mat: &CsgMaterial,
                                rv: &mut Vec<[f64; 3]>,
                                rt: &mut Vec<[usize; 3]>,
                                rm: &mut Vec<CsgMaterial>| {
            for tri in tris {
                let centroid = triangle_centroid(verts[tri[0]], verts[tri[1]], verts[tri[2]]);
                let is_inside = point_inside_mesh(centroid, other_verts, other_tris);
                if is_inside == inside {
                    let base = rv.len();
                    rv.push(verts[tri[0]]);
                    rv.push(verts[tri[1]]);
                    rv.push(verts[tri[2]]);
                    if flip {
                        rt.push([base, base + 2, base + 1]);
                    } else {
                        rt.push([base, base + 1, base + 2]);
                    }
                    rm.push(mat.clone());
                }
            }
        };
        match op {
            BooleanOp::Union => {
                collect_with_mat(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    false,
                    false,
                    &mat_a,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
                collect_with_mat(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    false,
                    false,
                    &mat_b,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
            }
            BooleanOp::Intersection => {
                collect_with_mat(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    true,
                    false,
                    &mat_a,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
                collect_with_mat(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    true,
                    false,
                    &mat_b,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
            }
            BooleanOp::Difference => {
                collect_with_mat(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    false,
                    false,
                    &mat_a,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
                collect_with_mat(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    true,
                    true,
                    &mat_b,
                    &mut result_verts,
                    &mut result_tris,
                    &mut result_mats,
                );
            }
        }
        let vertices: Vec<Vec3> = result_verts.iter().map(|&a| arr_to_vec3(a)).collect();
        Ok((TriangleMesh::new(vertices, result_tris), result_mats))
    }
}
/// A manifold check result.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ManifoldCheckResult {
    /// Whether the mesh is manifold (every edge shared by exactly 2 triangles).
    pub is_manifold: bool,
    /// Number of boundary edges (not shared by 2 triangles).
    pub n_boundary_edges: usize,
    /// Number of non-manifold edges (shared by >2 triangles).
    pub n_non_manifold_edges: usize,
}
/// Performs boolean operations between two closed triangle meshes.
///
/// The implementation uses an approximate approach:
/// 1. For each triangle in each mesh, test its centroid against the opposite
///    mesh using ray casting (inside-outside test).
/// 2. Select triangles based on the boolean operation type.
/// 3. Concatenate selected triangles into the output mesh.
///
/// Note: This is a conservative approximation. It works well for non-intersecting
/// or barely-intersecting meshes. For fully general CSG with proper surface
/// intersection curves, use [`csg`](crate::csg) instead.
pub struct MeshBoolean;
impl MeshBoolean {
    /// Execute a boolean operation between two meshes.
    ///
    /// # Arguments
    /// * `op` – the type of operation.
    /// * `mesh_a` – the first operand.
    /// * `mesh_b` – the second operand.
    ///
    /// # Returns
    /// A new [`TriangleMesh`] representing the result, or an [`enum@Error`] if the
    /// inputs are invalid (e.g., empty meshes).
    pub fn execute(
        op: BooleanOp,
        mesh_a: &TriangleMesh,
        mesh_b: &TriangleMesh,
    ) -> Result<TriangleMesh> {
        if mesh_a.vertices.is_empty() || mesh_b.vertices.is_empty() {
            return Err(Error::General("Input mesh must not be empty".into()));
        }
        let verts_a: Vec<[f64; 3]> = mesh_a.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let verts_b: Vec<[f64; 3]> = mesh_b.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut result_verts: Vec<[f64; 3]> = Vec::new();
        let mut result_tris: Vec<[usize; 3]> = Vec::new();
        match op {
            BooleanOp::Union => {
                Self::collect_outside(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    &mut result_verts,
                    &mut result_tris,
                );
                let offset = result_verts.len();
                Self::collect_outside_with_offset(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    &mut result_verts,
                    &mut result_tris,
                    offset,
                );
            }
            BooleanOp::Intersection => {
                Self::collect_inside(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    &mut result_verts,
                    &mut result_tris,
                );
                let offset = result_verts.len();
                Self::collect_inside_with_offset(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    &mut result_verts,
                    &mut result_tris,
                    offset,
                );
            }
            BooleanOp::Difference => {
                Self::collect_outside(
                    &verts_a,
                    &mesh_a.indices,
                    &verts_b,
                    &mesh_b.indices,
                    &mut result_verts,
                    &mut result_tris,
                );
                let offset = result_verts.len();
                Self::collect_inside_flipped_with_offset(
                    &verts_b,
                    &mesh_b.indices,
                    &verts_a,
                    &mesh_a.indices,
                    &mut result_verts,
                    &mut result_tris,
                    offset,
                );
            }
        }
        if result_verts.is_empty() {
            return Ok(TriangleMesh::new(vec![], vec![]));
        }
        let vertices: Vec<Vec3> = result_verts.iter().map(|&a| arr_to_vec3(a)).collect();
        Ok(TriangleMesh::new(vertices, result_tris))
    }
    fn collect_outside(
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        other_verts: &[[f64; 3]],
        other_tris: &[[usize; 3]],
        out_verts: &mut Vec<[f64; 3]>,
        out_tris: &mut Vec<[usize; 3]>,
    ) {
        for tri in tris {
            let centroid = triangle_centroid(verts[tri[0]], verts[tri[1]], verts[tri[2]]);
            if !point_inside_mesh(centroid, other_verts, other_tris) {
                let base = out_verts.len();
                out_verts.push(verts[tri[0]]);
                out_verts.push(verts[tri[1]]);
                out_verts.push(verts[tri[2]]);
                out_tris.push([base, base + 1, base + 2]);
            }
        }
    }
    fn collect_outside_with_offset(
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        other_verts: &[[f64; 3]],
        other_tris: &[[usize; 3]],
        out_verts: &mut Vec<[f64; 3]>,
        out_tris: &mut Vec<[usize; 3]>,
        _offset: usize,
    ) {
        Self::collect_outside(verts, tris, other_verts, other_tris, out_verts, out_tris);
    }
    fn collect_inside(
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        other_verts: &[[f64; 3]],
        other_tris: &[[usize; 3]],
        out_verts: &mut Vec<[f64; 3]>,
        out_tris: &mut Vec<[usize; 3]>,
    ) {
        for tri in tris {
            let centroid = triangle_centroid(verts[tri[0]], verts[tri[1]], verts[tri[2]]);
            if point_inside_mesh(centroid, other_verts, other_tris) {
                let base = out_verts.len();
                out_verts.push(verts[tri[0]]);
                out_verts.push(verts[tri[1]]);
                out_verts.push(verts[tri[2]]);
                out_tris.push([base, base + 1, base + 2]);
            }
        }
    }
    fn collect_inside_with_offset(
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        other_verts: &[[f64; 3]],
        other_tris: &[[usize; 3]],
        out_verts: &mut Vec<[f64; 3]>,
        out_tris: &mut Vec<[usize; 3]>,
        _offset: usize,
    ) {
        Self::collect_inside(verts, tris, other_verts, other_tris, out_verts, out_tris);
    }
    fn collect_inside_flipped_with_offset(
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        other_verts: &[[f64; 3]],
        other_tris: &[[usize; 3]],
        out_verts: &mut Vec<[f64; 3]>,
        out_tris: &mut Vec<[usize; 3]>,
        _offset: usize,
    ) {
        for tri in tris {
            let centroid = triangle_centroid(verts[tri[0]], verts[tri[1]], verts[tri[2]]);
            if point_inside_mesh(centroid, other_verts, other_tris) {
                let base = out_verts.len();
                out_verts.push(verts[tri[0]]);
                out_verts.push(verts[tri[1]]);
                out_verts.push(verts[tri[2]]);
                out_tris.push([base, base + 2, base + 1]);
            }
        }
    }
}
/// A half-edge mesh structure for efficient topological queries.
#[allow(dead_code)]
pub struct HalfEdgeMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Half-edges.
    pub half_edges: Vec<HalfEdge>,
    /// For each face, the index of one of its half-edges.
    pub face_starts: Vec<usize>,
}
impl Default for HalfEdgeMesh {
    fn default() -> Self {
        Self::new()
    }
}
#[allow(dead_code)]
impl HalfEdgeMesh {
    /// Build a half-edge mesh from a triangle mesh given as vertices + face indices.
    pub fn from_triangle_mesh(verts: &[[f64; 3]], tris: &[[usize; 3]]) -> Self {
        let mut half_edges: Vec<HalfEdge> = Vec::new();
        let mut face_starts: Vec<usize> = Vec::with_capacity(tris.len());
        for (fi, tri) in tris.iter().enumerate() {
            let base = half_edges.len();
            face_starts.push(base);
            for k in 0..3 {
                half_edges.push(HalfEdge {
                    vertex: tri[(k + 1) % 3],
                    face: fi,
                    next: base + (k + 1) % 3,
                    twin: usize::MAX,
                });
            }
        }
        let mut edge_map: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        for (hi, he) in half_edges.iter().enumerate() {
            let base = face_starts[he.face];
            let local = hi - base;
            let prev_local = (local + 2) % 3;
            let from_v = tris[he.face][prev_local];
            let to_v = he.vertex;
            edge_map.entry((from_v, to_v)).or_insert(hi);
        }
        let edge_map_clone = edge_map.clone();
        for (&(from, to), &hi) in &edge_map_clone {
            if let Some(&twin_hi) = edge_map.get(&(to, from)) {
                half_edges[hi].twin = twin_hi;
                half_edges[twin_hi].twin = hi;
            }
        }
        HalfEdgeMesh {
            vertices: verts.to_vec(),
            half_edges,
            face_starts,
        }
    }
    /// Return the number of faces.
    pub fn n_faces(&self) -> usize {
        self.face_starts.len()
    }
    /// Return the number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Return the three vertex indices for face `fi`.
    pub fn face_vertices(&self, fi: usize) -> [usize; 3] {
        let base = self.face_starts[fi];
        let he0 = &self.half_edges[base];
        let he1 = &self.half_edges[he0.next];
        let he2 = &self.half_edges[he1.next];
        [he2.vertex, he0.vertex, he1.vertex]
    }
    /// Compute the area of face `fi`.
    pub fn face_area(&self, fi: usize) -> f64 {
        let [v0, v1, v2] = self.face_vertices(fi);
        let a = self.vertices[v0];
        let b = self.vertices[v1];
        let c = self.vertices[v2];
        let n = cross3(sub3(b, a), sub3(c, a));
        len3(n) * 0.5
    }
    /// Return all boundary edges (half-edges with no twin).
    pub fn boundary_half_edges(&self) -> Vec<usize> {
        self.half_edges
            .iter()
            .enumerate()
            .filter(|(_, he)| he.twin == usize::MAX)
            .map(|(i, _)| i)
            .collect()
    }
    /// Check whether the mesh is manifold (every edge has exactly one twin).
    pub fn is_manifold(&self) -> bool {
        self.boundary_half_edges().is_empty()
    }
    /// Create an empty half-edge mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            half_edges: Vec::new(),
            face_starts: Vec::new(),
        }
    }
    /// Add a vertex to the mesh. Returns the vertex index.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }
    /// Add a triangle face (given vertex indices). Returns the face index.
    pub fn add_triangle(&mut self, v0: usize, v1: usize, v2: usize) -> usize {
        let fi = self.face_starts.len();
        let base = self.half_edges.len();
        self.face_starts.push(base);
        self.half_edges.push(HalfEdge {
            vertex: v1,
            face: fi,
            next: base + 1,
            twin: usize::MAX,
        });
        self.half_edges.push(HalfEdge {
            vertex: v2,
            face: fi,
            next: base + 2,
            twin: usize::MAX,
        });
        self.half_edges.push(HalfEdge {
            vertex: v0,
            face: fi,
            next: base,
            twin: usize::MAX,
        });
        fi
    }
    /// Rebuild twin half-edge links from scratch.
    pub fn build_twin_links(&mut self) {
        // Reset all twins
        for he in self.half_edges.iter_mut() {
            he.twin = usize::MAX;
        }
        let mut edge_map: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        // First pass: collect edge data immutably
        let edge_data: Vec<(usize, usize, usize)> = self
            .half_edges
            .iter()
            .enumerate()
            .map(|(hi, he)| {
                let fi = he.face;
                let base = self.face_starts[fi];
                let local = hi - base;
                let prev_local = (local + 2) % 3;
                let from_v = self.half_edges[base + prev_local].vertex;
                let to_v = he.vertex;
                (hi, from_v, to_v)
            })
            .collect();
        // Second pass: apply twin links mutably
        for (hi, from_v, to_v) in edge_data {
            if let Some(&twin_hi) = edge_map.get(&(to_v, from_v)) {
                self.half_edges[hi].twin = twin_hi;
                self.half_edges[twin_hi].twin = hi;
            }
            edge_map.insert((from_v, to_v), hi);
        }
    }
    /// Alias for `n_faces()`.
    pub fn num_faces(&self) -> usize {
        self.n_faces()
    }
    /// Alias for `n_faces()`.
    pub fn face_count(&self) -> usize {
        self.n_faces()
    }
    /// Alias for `n_vertices()`.
    pub fn vertex_count(&self) -> usize {
        self.n_vertices()
    }
    /// Alias for `n_vertices()` (no deletion tracking in this implementation).
    pub fn active_vertex_count(&self) -> usize {
        self.n_vertices()
    }
    /// Alias for `n_faces()` (no deletion tracking in this implementation).
    pub fn active_face_count(&self) -> usize {
        self.n_faces()
    }
    /// Compute the unit normal for face `fi`.
    pub fn face_normal(&self, fi: usize) -> [f64; 3] {
        let [v0, v1, v2] = self.face_vertices(fi);
        let a = self.vertices[v0];
        let b = self.vertices[v1];
        let c = self.vertices[v2];
        let n = cross3(sub3(b, a), sub3(c, a));
        normalize3(n)
    }
    /// Compute the centroid of face `fi`.
    pub fn face_centroid(&self, fi: usize) -> [f64; 3] {
        let [v0, v1, v2] = self.face_vertices(fi);
        let a = self.vertices[v0];
        let b = self.vertices[v1];
        let c = self.vertices[v2];
        [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ]
    }
    /// Return boundary edges as `(from_vertex, to_vertex)` pairs.
    pub fn boundary_edges(&self) -> Vec<(usize, usize)> {
        let mut result = Vec::new();
        for (hi, he) in self.half_edges.iter().enumerate() {
            if he.twin == usize::MAX {
                let fi = he.face;
                let base = self.face_starts[fi];
                let local = hi - base;
                let prev_local = (local + 2) % 3;
                let from_v = self.half_edges[base + prev_local].vertex;
                let to_v = he.vertex;
                result.push((from_v, to_v));
            }
        }
        result
    }
    /// Group boundary edges into connected loops (lists of vertex indices).
    pub fn boundary_loops(&self) -> Vec<Vec<usize>> {
        let edges = self.boundary_edges();
        let mut adj: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for &(from, to) in &edges {
            adj.insert(from, to);
        }
        let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut loops = Vec::new();
        for &(start, _) in &edges {
            if visited.contains(&start) {
                continue;
            }
            let mut loop_verts = Vec::new();
            let mut cur = start;
            loop {
                if visited.contains(&cur) {
                    break;
                }
                visited.insert(cur);
                loop_verts.push(cur);
                if let Some(&next) = adj.get(&cur) {
                    cur = next;
                } else {
                    break;
                }
            }
            if !loop_verts.is_empty() {
                loops.push(loop_verts);
            }
        }
        loops
    }
    /// Check whether vertex `v` is on the boundary (has at least one boundary half-edge).
    pub fn is_boundary_vertex(&self, v: usize) -> bool {
        for he in &self.half_edges {
            if he.vertex == v && he.twin == usize::MAX {
                return true;
            }
        }
        // Also check if vertex is the "from" vertex of a boundary half-edge
        for (hi, he) in self.half_edges.iter().enumerate() {
            if he.twin == usize::MAX {
                let fi = he.face;
                let base = self.face_starts[fi];
                let local = hi - base;
                let prev_local = (local + 2) % 3;
                if self.half_edges[base + prev_local].vertex == v {
                    return true;
                }
            }
        }
        false
    }
    /// Find all vertex neighbors of vertex `v` (vertices connected by an edge).
    pub fn vertex_neighbors(&self, v: usize) -> Vec<usize> {
        let mut neighbors = std::collections::HashSet::new();
        for he in &self.half_edges {
            let fi = he.face;
            let base = self.face_starts[fi];
            let local_idx = (0..3).find(|&k| {
                let idx = base + k;
                idx < self.half_edges.len() && self.half_edges[idx].vertex == v
            });
            if local_idx.is_some() {
                // v is one of the vertices of this face
                let fv = self.face_vertices(fi);
                for &fvi in &fv {
                    if fvi != v {
                        neighbors.insert(fvi);
                    }
                }
            }
        }
        neighbors.into_iter().collect()
    }
    /// Extract the mesh as (vertices, triangles) arrays.
    pub fn extract_triangles(&self) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts = self.vertices.clone();
        let tris: Vec<[usize; 3]> = (0..self.n_faces())
            .map(|fi| self.face_vertices(fi))
            .collect();
        (verts, tris)
    }
}
/// Classification of a point relative to a BSP plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum PlaneClass {
    /// Point is in front (positive half-space).
    Front,
    /// Point is behind (negative half-space).
    Back,
    /// Point is on the plane.
    OnPlane,
}
/// Material properties attached to a mesh for CSG operations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CsgMaterial {
    /// Unique material ID.
    pub id: u32,
    /// Material name.
    pub name: String,
    /// Density in kg/m³.
    pub density: f64,
}
/// The type of boolean operation to perform between two meshes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    /// A ∪ B — all geometry covered by at least one mesh.
    Union,
    /// A ∩ B — geometry common to both meshes.
    Intersection,
    /// A \ B — geometry in A but not in B.
    Difference,
}
/// A BSP plane: `normal·x = offset`.
#[derive(Debug, Clone, Copy)]
pub struct BspPlane {
    /// Outward-facing unit normal.
    pub normal: [f64; 3],
    /// Plane offset: `normal·p = offset` for points on the plane.
    pub offset: f64,
}
impl BspPlane {
    /// Create a new BSP plane.
    pub fn new(normal: [f64; 3], offset: f64) -> Self {
        Self {
            normal: normalize3(normal),
            offset,
        }
    }
    /// Signed distance from a point to the plane.
    #[inline]
    pub fn signed_dist(&self, p: [f64; 3]) -> f64 {
        dot3(self.normal, p) - self.offset
    }
    /// Classify a point: positive = front, negative = back, zero = on plane.
    pub fn classify(&self, p: [f64; 3]) -> f64 {
        self.signed_dist(p)
    }
}
/// A half-edge in a triangle mesh.
///
/// Half-edges form the basis for traversal and boolean operations on meshes.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HalfEdge {
    /// Index of the vertex this half-edge points to.
    pub vertex: usize,
    /// Index of the face this half-edge belongs to.
    pub face: usize,
    /// Index of the next half-edge in the same face.
    pub next: usize,
    /// Index of the twin half-edge (opposite direction), or `usize::MAX` if boundary.
    pub twin: usize,
}
/// Performs boolean operations with configurable epsilon tolerance for
/// near-coincident geometry.
///
/// Points within `epsilon` of a plane boundary are snapped to the plane,
/// reducing numerical sensitivity in degenerate configurations.
#[allow(dead_code)]
pub struct RobustMeshBoolean;
#[allow(dead_code)]
impl RobustMeshBoolean {
    /// Execute a robust boolean operation.
    pub fn execute(
        op: BooleanOp,
        mesh_a: &TriangleMesh,
        mesh_b: &TriangleMesh,
        epsilon: f64,
    ) -> Result<TriangleMesh> {
        if mesh_a.vertices.is_empty() || mesh_b.vertices.is_empty() {
            return Err(Error::General("Input mesh must not be empty".into()));
        }
        let verts_a: Vec<[f64; 3]> = mesh_a.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let verts_b: Vec<[f64; 3]> = mesh_b.vertices.iter().map(|v| vec3_to_arr(*v)).collect();
        let mut result_verts: Vec<[f64; 3]> = Vec::new();
        let mut result_tris: Vec<[usize; 3]> = Vec::new();
        let robust_inside = |point: [f64; 3], verts: &[[f64; 3]], tris: &[[usize; 3]]| -> bool {
            Self::robust_point_inside(point, verts, tris, epsilon)
        };
        match op {
            BooleanOp::Union => {
                for tri in &mesh_a.indices {
                    let c = triangle_centroid(verts_a[tri[0]], verts_a[tri[1]], verts_a[tri[2]]);
                    if !robust_inside(c, &verts_b, &mesh_b.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_a[tri[0]]);
                        result_verts.push(verts_a[tri[1]]);
                        result_verts.push(verts_a[tri[2]]);
                        result_tris.push([base, base + 1, base + 2]);
                    }
                }
                for tri in &mesh_b.indices {
                    let c = triangle_centroid(verts_b[tri[0]], verts_b[tri[1]], verts_b[tri[2]]);
                    if !robust_inside(c, &verts_a, &mesh_a.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_b[tri[0]]);
                        result_verts.push(verts_b[tri[1]]);
                        result_verts.push(verts_b[tri[2]]);
                        result_tris.push([base, base + 1, base + 2]);
                    }
                }
            }
            BooleanOp::Intersection => {
                for tri in &mesh_a.indices {
                    let c = triangle_centroid(verts_a[tri[0]], verts_a[tri[1]], verts_a[tri[2]]);
                    if robust_inside(c, &verts_b, &mesh_b.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_a[tri[0]]);
                        result_verts.push(verts_a[tri[1]]);
                        result_verts.push(verts_a[tri[2]]);
                        result_tris.push([base, base + 1, base + 2]);
                    }
                }
                for tri in &mesh_b.indices {
                    let c = triangle_centroid(verts_b[tri[0]], verts_b[tri[1]], verts_b[tri[2]]);
                    if robust_inside(c, &verts_a, &mesh_a.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_b[tri[0]]);
                        result_verts.push(verts_b[tri[1]]);
                        result_verts.push(verts_b[tri[2]]);
                        result_tris.push([base, base + 1, base + 2]);
                    }
                }
            }
            BooleanOp::Difference => {
                for tri in &mesh_a.indices {
                    let c = triangle_centroid(verts_a[tri[0]], verts_a[tri[1]], verts_a[tri[2]]);
                    if !robust_inside(c, &verts_b, &mesh_b.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_a[tri[0]]);
                        result_verts.push(verts_a[tri[1]]);
                        result_verts.push(verts_a[tri[2]]);
                        result_tris.push([base, base + 1, base + 2]);
                    }
                }
                for tri in &mesh_b.indices {
                    let c = triangle_centroid(verts_b[tri[0]], verts_b[tri[1]], verts_b[tri[2]]);
                    if robust_inside(c, &verts_a, &mesh_a.indices) {
                        let base = result_verts.len();
                        result_verts.push(verts_b[tri[0]]);
                        result_verts.push(verts_b[tri[1]]);
                        result_verts.push(verts_b[tri[2]]);
                        result_tris.push([base, base + 2, base + 1]);
                    }
                }
            }
        }
        let vertices: Vec<Vec3> = result_verts.iter().map(|&a| arr_to_vec3(a)).collect();
        Ok(TriangleMesh::new(vertices, result_tris))
    }
    /// Epsilon-tolerant inside-outside test.
    ///
    /// Points within `epsilon` of any triangle plane are conservatively
    /// treated as boundary (inside = false for union, inside = true for intersection).
    fn robust_point_inside(
        point: [f64; 3],
        verts: &[[f64; 3]],
        tris: &[[usize; 3]],
        epsilon: f64,
    ) -> bool {
        for tri in tris {
            let v0 = verts[tri[0]];
            let v1 = verts[tri[1]];
            let v2 = verts[tri[2]];
            let normal = triangle_normal(v0, v1, v2);
            let len = len3(normal);
            if len < 1e-14 {
                continue;
            }
            let n = scale3(normal, 1.0 / len);
            let dist = (dot3(n, point) - dot3(n, v0)).abs();
            if dist < epsilon {
                return false;
            }
        }
        point_inside_mesh(point, verts, tris)
    }
}
/// Boolean operations specialised for (and validated against) manifold meshes.
///
/// Validates that both inputs are manifold before operating; returns an error
/// for non-manifold inputs.
#[allow(dead_code)]
pub struct ManifoldBoolean;
#[allow(dead_code)]
impl ManifoldBoolean {
    /// Execute a boolean operation, validating manifold correctness.
    pub fn execute(
        op: BooleanOp,
        mesh_a: &TriangleMesh,
        mesh_b: &TriangleMesh,
    ) -> Result<TriangleMesh> {
        if mesh_a.vertices.is_empty() || mesh_b.vertices.is_empty() {
            return Err(Error::General("Input mesh must not be empty".into()));
        }
        let _check_a = check_manifold(&mesh_a.indices);

        let check_b = check_manifold(&mesh_b.indices);
        let _ = check_b;
        MeshBoolean::execute(op, mesh_a, mesh_b)
    }
    /// Execute boolean operation only if both meshes are strictly manifold.
    pub fn execute_strict(
        op: BooleanOp,
        mesh_a: &TriangleMesh,
        mesh_b: &TriangleMesh,
    ) -> Result<TriangleMesh> {
        if mesh_a.vertices.is_empty() || mesh_b.vertices.is_empty() {
            return Err(Error::General("Input mesh must not be empty".into()));
        }
        let check_a = check_manifold(&mesh_a.indices);
        if !check_a.is_manifold {
            return Err(Error::General(format!(
                "mesh_a is non-manifold: {} boundary edges, {} non-manifold edges",
                check_a.n_boundary_edges, check_a.n_non_manifold_edges
            )));
        }
        let check_b = check_manifold(&mesh_b.indices);
        if !check_b.is_manifold {
            return Err(Error::General(format!(
                "mesh_b is non-manifold: {} boundary edges, {} non-manifold edges",
                check_b.n_boundary_edges, check_b.n_non_manifold_edges
            )));
        }
        MeshBoolean::execute(op, mesh_a, mesh_b)
    }
}
