//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::functions::{
    DEFAULT_CELL_SIZE, EPSILON, faces_share_vertex, v3_add, v3_cross, v3_max, v3_min, v3_norm,
    v3_scale, v3_sub, vertex_face_contact,
};

/// Lightweight AABB for deformable BVH nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeformAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl DeformAabb {
    /// Create a new AABB from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Create an empty (inverted) AABB.
    pub fn empty() -> Self {
        Self {
            min: [f64::MAX; 3],
            max: [f64::MIN; 3],
        }
    }
    /// Expand to include a point.
    pub fn expand_point(&mut self, p: [f64; 3]) {
        self.min = v3_min(self.min, p);
        self.max = v3_max(self.max, p);
    }
    /// Expand by a uniform margin.
    pub fn expand_margin(&mut self, m: f64) {
        self.min = v3_sub(self.min, [m; 3]);
        self.max = v3_add(self.max, [m; 3]);
    }
    /// Union of two AABBs.
    pub fn union(&self, other: &Self) -> Self {
        Self {
            min: v3_min(self.min, other.min),
            max: v3_max(self.max, other.max),
        }
    }
    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Test overlap with another AABB expanded by a margin on both sides.
    pub fn overlaps_with_margin(&self, other: &Self, margin: f64) -> bool {
        self.min[0] - margin <= other.max[0] + margin
            && self.max[0] + margin >= other.min[0] - margin
            && self.min[1] - margin <= other.max[1] + margin
            && self.max[1] + margin >= other.min[1] - margin
            && self.min[2] - margin <= other.max[2] + margin
            && self.max[2] + margin >= other.min[2] - margin
    }
    /// Surface area heuristic cost.
    pub fn surface_area(&self) -> f64 {
        let d = v3_sub(self.max, self.min);
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }
    /// Compute AABB for a triangle.
    pub fn from_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        Self {
            min: v3_min(v3_min(a, b), c),
            max: v3_max(v3_max(a, b), c),
        }
    }
    /// Centre of the AABB.
    pub fn centre(&self) -> [f64; 3] {
        v3_scale(v3_add(self.min, self.max), 0.5)
    }
}
/// A position-based contact constraint for deformable bodies.
#[derive(Debug, Clone)]
pub struct DeformConstraint {
    /// Vertex index in the mesh.
    pub vertex: usize,
    /// Target position (on the surface or corrected).
    pub target: [f64; 3],
    /// Contact normal.
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
}
/// Penalty force parameters.
#[derive(Debug, Clone, Copy)]
pub struct PenaltyParams {
    /// Stiffness coefficient.
    pub stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
    /// Friction coefficient.
    pub friction: f64,
}
/// Deformable triangle mesh with per-vertex positions and velocities.
#[derive(Debug, Clone)]
pub struct DeformableMesh {
    /// Current vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Previous frame vertex positions (for CCD).
    pub prev_positions: Vec<[f64; 3]>,
    /// Current vertex velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Triangle faces.
    pub faces: Vec<TriFace>,
    /// Per-vertex inverse mass (0 = pinned).
    pub inv_mass: Vec<f64>,
}
impl DeformableMesh {
    /// Create a new deformable mesh.
    pub fn new(positions: Vec<[f64; 3]>, faces: Vec<TriFace>, inv_mass: Vec<f64>) -> Self {
        let n = positions.len();
        Self {
            prev_positions: positions.clone(),
            velocities: vec![[0.0; 3]; n],
            positions,
            faces,
            inv_mass,
        }
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.positions.len()
    }
    /// Number of triangles.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Compute the face normal (not normalised) for face at `idx`.
    pub fn face_normal(&self, idx: usize) -> [f64; 3] {
        let f = &self.faces[idx];
        let ab = v3_sub(self.positions[f.v1], self.positions[f.v0]);
        let ac = v3_sub(self.positions[f.v2], self.positions[f.v0]);
        v3_cross(ab, ac)
    }
    /// Compute normalised face normal.
    pub fn face_normal_unit(&self, idx: usize) -> [f64; 3] {
        v3_norm(self.face_normal(idx))
    }
    /// AABB for a single triangle face.
    pub fn face_aabb(&self, idx: usize) -> DeformAabb {
        let f = &self.faces[idx];
        DeformAabb::from_triangle(
            self.positions[f.v0],
            self.positions[f.v1],
            self.positions[f.v2],
        )
    }
    /// Overall bounding box of the mesh.
    pub fn bounding_box(&self) -> DeformAabb {
        let mut aabb = DeformAabb::empty();
        for &p in &self.positions {
            aabb.expand_point(p);
        }
        aabb
    }
    /// Collect unique edges from all faces.
    pub fn collect_edges(&self) -> Vec<(usize, usize)> {
        let mut set = std::collections::HashSet::new();
        for f in &self.faces {
            for e in f.edges() {
                set.insert(e);
            }
        }
        set.into_iter().collect()
    }
    /// Store current positions as previous (call before time-step).
    pub fn snapshot_positions(&mut self) {
        self.prev_positions = self.positions.clone();
    }
}
/// Type of deformable contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeformContactType {
    /// Vertex penetrating a face.
    VertexFace,
    /// Two edges overlapping.
    EdgeEdge,
    /// Vertex vs rigid body.
    VertexRigid,
    /// Face vs rigid body.
    FaceRigid,
}
/// A rigid body represented as an infinite plane.
#[derive(Debug, Clone, Copy)]
pub struct RigidPlane {
    /// Plane normal (unit).
    pub normal: [f64; 3],
    /// Signed distance from origin.
    pub offset: f64,
}
/// A triangle face defined by three vertex indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TriFace {
    /// Index of vertex 0.
    pub v0: usize,
    /// Index of vertex 1.
    pub v1: usize,
    /// Index of vertex 2.
    pub v2: usize,
}
impl TriFace {
    /// Create a new triangle face.
    pub fn new(v0: usize, v1: usize, v2: usize) -> Self {
        Self { v0, v1, v2 }
    }
    /// Return sorted edge pairs for adjacency comparisons.
    pub fn edges(&self) -> [(usize, usize); 3] {
        [
            (self.v0.min(self.v1), self.v0.max(self.v1)),
            (self.v1.min(self.v2), self.v1.max(self.v2)),
            (self.v2.min(self.v0), self.v2.max(self.v0)),
        ]
    }
}
/// Spatial hash grid for accelerating deformable collision queries.
#[derive(Debug, Clone)]
pub struct SpatialHash {
    /// Cell size (uniform).
    pub cell_size: f64,
    /// Inverse cell size.
    pub(super) inv_cell: f64,
    /// Map from cell key to list of triangle indices.
    pub(super) cells: HashMap<CellKey, Vec<usize>>,
}
impl SpatialHash {
    /// Create a new spatial hash with given cell size.
    pub fn new(cell_size: f64) -> Self {
        let cs = if cell_size < EPSILON {
            DEFAULT_CELL_SIZE
        } else {
            cell_size
        };
        Self {
            cell_size: cs,
            inv_cell: 1.0 / cs,
            cells: HashMap::new(),
        }
    }
    /// Clear the spatial hash.
    pub fn clear(&mut self) {
        self.cells.clear();
    }
    /// Convert a position to a cell key.
    fn cell_key(&self, p: [f64; 3]) -> CellKey {
        CellKey {
            x: (p[0] * self.inv_cell).floor() as i64,
            y: (p[1] * self.inv_cell).floor() as i64,
            z: (p[2] * self.inv_cell).floor() as i64,
        }
    }
    /// Insert a triangle into the hash. Rasterises its AABB into cells.
    pub fn insert_triangle(&mut self, idx: usize, aabb: &DeformAabb) {
        let min_key = self.cell_key(aabb.min);
        let max_key = self.cell_key(aabb.max);
        for x in min_key.x..=max_key.x {
            for y in min_key.y..=max_key.y {
                for z in min_key.z..=max_key.z {
                    let key = CellKey { x, y, z };
                    self.cells.entry(key).or_default().push(idx);
                }
            }
        }
    }
    /// Build the spatial hash from a deformable mesh.
    pub fn build_from_mesh(&mut self, mesh: &DeformableMesh) {
        self.clear();
        for fi in 0..mesh.num_faces() {
            let aabb = mesh.face_aabb(fi);
            self.insert_triangle(fi, &aabb);
        }
    }
    /// Query all triangles overlapping a given AABB.
    pub fn query_aabb(&self, aabb: &DeformAabb) -> Vec<usize> {
        let mut result = Vec::new();
        let min_key = self.cell_key(aabb.min);
        let max_key = self.cell_key(aabb.max);
        let mut seen = std::collections::HashSet::new();
        for x in min_key.x..=max_key.x {
            for y in min_key.y..=max_key.y {
                for z in min_key.z..=max_key.z {
                    let key = CellKey { x, y, z };
                    if let Some(items) = self.cells.get(&key) {
                        for &i in items {
                            if seen.insert(i) {
                                result.push(i);
                            }
                        }
                    }
                }
            }
        }
        result
    }
    /// Query potential self-collision pairs using the hash.
    pub fn query_self_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = std::collections::HashSet::new();
        for tris in self.cells.values() {
            for (i, &a) in tris.iter().enumerate() {
                for &b in &tris[i + 1..] {
                    let pair = if a < b { (a, b) } else { (b, a) };
                    pairs.insert(pair);
                }
            }
        }
        pairs.into_iter().collect()
    }
    /// Detect self-collision using spatial hashing (alternative to BVH).
    pub fn detect_self_collision_hashed(
        &self,
        mesh: &DeformableMesh,
        config: &DeformCollisionConfig,
    ) -> Vec<DeformContact> {
        let mut contacts = Vec::new();
        if !config.self_collision {
            return contacts;
        }
        let pairs = self.query_self_pairs();
        for (fi_a, fi_b) in pairs {
            let fa = &mesh.faces[fi_a];
            let fb = &mesh.faces[fi_b];
            if config.adjacency_skip && faces_share_vertex(fa, fb) {
                continue;
            }
            let tb0 = mesh.positions[fb.v0];
            let tb1 = mesh.positions[fb.v1];
            let tb2 = mesh.positions[fb.v2];
            for &vi in &[fa.v0, fa.v1, fa.v2] {
                let p = mesh.positions[vi];
                if let Some(c) = vertex_face_contact(p, tb0, tb1, tb2, config.thickness, vi, fi_b) {
                    contacts.push(c);
                }
            }
        }
        contacts
    }
    /// Number of occupied cells.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }
}
/// A contact manifold collecting multiple contact points for a deformable pair.
#[derive(Debug, Clone)]
pub struct DeformManifold {
    /// Body A identifier.
    pub body_a: usize,
    /// Body B identifier.
    pub body_b: usize,
    /// Contact points in this manifold.
    pub contacts: Vec<DeformContact>,
    /// Maximum contacts to retain.
    pub max_contacts: usize,
}
impl DeformManifold {
    /// Create a new empty manifold.
    pub fn new(body_a: usize, body_b: usize, max_contacts: usize) -> Self {
        Self {
            body_a,
            body_b,
            contacts: Vec::new(),
            max_contacts,
        }
    }
    /// Add a contact, keeping only the deepest if over limit.
    pub fn add_contact(&mut self, contact: DeformContact) {
        self.contacts.push(contact);
        if self.contacts.len() > self.max_contacts {
            self.reduce();
        }
    }
    /// Reduce contacts to `max_contacts` by keeping the deepest ones.
    pub fn reduce(&mut self) {
        self.contacts.sort_by(|a, b| {
            b.depth
                .partial_cmp(&a.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.contacts.truncate(self.max_contacts);
    }
    /// Total number of contacts.
    pub fn len(&self) -> usize {
        self.contacts.len()
    }
    /// Whether the manifold has no contacts.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }
    /// Average contact normal.
    pub fn average_normal(&self) -> [f64; 3] {
        if self.contacts.is_empty() {
            return [0.0, 1.0, 0.0];
        }
        let mut sum = [0.0; 3];
        for c in &self.contacts {
            sum = v3_add(sum, c.normal);
        }
        v3_norm(sum)
    }
    /// Maximum penetration depth.
    pub fn max_depth(&self) -> f64 {
        self.contacts
            .iter()
            .map(|c| c.depth)
            .fold(0.0_f64, f64::max)
    }
}
/// Configuration for deformable collision detection.
#[derive(Debug, Clone)]
pub struct DeformCollisionConfig {
    /// Proximity thickness for vertex-face tests.
    pub thickness: f64,
    /// Enable edge-edge tests.
    pub edge_edge: bool,
    /// Enable self-collision.
    pub self_collision: bool,
    /// Adjacency skip count (skip contacts between vertices sharing faces).
    pub adjacency_skip: bool,
}
/// A rigid body represented as a sphere for deformable-rigid tests.
#[derive(Debug, Clone, Copy)]
pub struct RigidSphere {
    /// Centre of the sphere.
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
}
/// A rigid body represented as an axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub struct RigidBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
/// Cell key for spatial hashing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CellKey {
    /// X cell index.
    pub x: i64,
    /// Y cell index.
    pub y: i64,
    /// Z cell index.
    pub z: i64,
}
/// Result of a continuous collision query.
#[derive(Debug, Clone)]
pub struct DeformCcdResult {
    /// Time of impact in \[0, 1\].
    pub toi: f64,
    /// Contact normal at impact.
    pub normal: [f64; 3],
    /// Contact position at impact.
    pub position: [f64; 3],
    /// Vertex or feature index on mesh A.
    pub feature_a: usize,
    /// Face or feature index on mesh B.
    pub feature_b: usize,
}
/// A BVH node for deformable meshes. Stores either leaf face indices or children.
#[derive(Debug, Clone)]
pub struct DeformBvhNode {
    /// Bounding box of this node.
    pub aabb: DeformAabb,
    /// Left child index (usize::MAX if leaf).
    pub left: usize,
    /// Right child index (usize::MAX if leaf).
    pub right: usize,
    /// Face indices (non-empty only for leaves).
    pub faces: Vec<usize>,
}
/// A contact point produced by deformable collision detection.
#[derive(Debug, Clone)]
pub struct DeformContact {
    /// World-space contact position.
    pub position: [f64; 3],
    /// Contact normal (pointing from body B to body A).
    pub normal: [f64; 3],
    /// Penetration depth (positive when overlapping).
    pub depth: f64,
    /// Barycentric coordinates on the face (if applicable).
    pub bary: [f64; 3],
    /// Type of contact.
    pub contact_type: DeformContactType,
    /// Index of the vertex or edge on mesh A.
    pub feature_a: usize,
    /// Index of the face or edge on mesh B.
    pub feature_b: usize,
}
impl DeformContact {
    /// Create a new deformable contact.
    pub fn new(
        position: [f64; 3],
        normal: [f64; 3],
        depth: f64,
        bary: [f64; 3],
        contact_type: DeformContactType,
        feature_a: usize,
        feature_b: usize,
    ) -> Self {
        Self {
            position,
            normal,
            depth,
            bary,
            contact_type,
            feature_a,
            feature_b,
        }
    }
}
/// Constraint solver parameters.
#[derive(Debug, Clone, Copy)]
pub struct ConstraintParams {
    /// Number of solver iterations.
    pub iterations: usize,
    /// Relaxation factor (0..1].
    pub relaxation: f64,
    /// Friction coefficient.
    pub friction: f64,
}
/// BVH tree built over a deformable mesh's triangles.
#[derive(Debug, Clone)]
pub struct DeformBvh {
    /// Flat node storage.
    pub nodes: Vec<DeformBvhNode>,
}
impl DeformBvh {
    /// Build a BVH from a deformable mesh.
    pub fn build(mesh: &DeformableMesh) -> Self {
        let face_indices: Vec<usize> = (0..mesh.num_faces()).collect();
        let mut nodes = Vec::new();
        Self::build_recursive(mesh, &face_indices, &mut nodes);
        Self { nodes }
    }
    /// Recursive BVH build helper.
    fn build_recursive(
        mesh: &DeformableMesh,
        indices: &[usize],
        nodes: &mut Vec<DeformBvhNode>,
    ) -> usize {
        let mut aabb = DeformAabb::empty();
        for &i in indices {
            aabb = aabb.union(&mesh.face_aabb(i));
        }
        if indices.len() <= 2 {
            let idx = nodes.len();
            nodes.push(DeformBvhNode {
                aabb,
                left: usize::MAX,
                right: usize::MAX,
                faces: indices.to_vec(),
            });
            return idx;
        }
        let d = v3_sub(aabb.max, aabb.min);
        let axis = if d[0] >= d[1] && d[0] >= d[2] {
            0
        } else if d[1] >= d[2] {
            1
        } else {
            2
        };
        let mid_val = (aabb.min[axis] + aabb.max[axis]) * 0.5;
        let mut left_set = Vec::new();
        let mut right_set = Vec::new();
        for &i in indices {
            let c = mesh.face_aabb(i).centre();
            if c[axis] < mid_val {
                left_set.push(i);
            } else {
                right_set.push(i);
            }
        }
        if left_set.is_empty() || right_set.is_empty() {
            let half = indices.len() / 2;
            left_set = indices[..half].to_vec();
            right_set = indices[half..].to_vec();
        }
        let node_idx = nodes.len();
        nodes.push(DeformBvhNode {
            aabb,
            left: usize::MAX,
            right: usize::MAX,
            faces: Vec::new(),
        });
        let left = Self::build_recursive(mesh, &left_set, nodes);
        let right = Self::build_recursive(mesh, &right_set, nodes);
        nodes[node_idx].left = left;
        nodes[node_idx].right = right;
        node_idx
    }
    /// Refit the BVH after vertex positions change (bottom-up).
    pub fn refit(&mut self, mesh: &DeformableMesh) {
        if self.nodes.is_empty() {
            return;
        }
        self.refit_node(mesh, 0);
    }
    /// Recursive refit helper.
    fn refit_node(&mut self, mesh: &DeformableMesh, idx: usize) -> DeformAabb {
        if self.nodes[idx].left == usize::MAX {
            let mut aabb = DeformAabb::empty();
            for &fi in &self.nodes[idx].faces {
                aabb = aabb.union(&mesh.face_aabb(fi));
            }
            self.nodes[idx].aabb = aabb;
            aabb
        } else {
            let l = self.nodes[idx].left;
            let r = self.nodes[idx].right;
            let la = self.refit_node(mesh, l);
            let ra = self.refit_node(mesh, r);
            let aabb = la.union(&ra);
            self.nodes[idx].aabb = aabb;
            aabb
        }
    }
    /// Query all face pairs whose AABBs overlap between `self` and `other`.
    pub fn query_overlap_pairs(&self, other: &DeformBvh) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return pairs;
        }
        Self::overlap_recurse(&self.nodes, 0, &other.nodes, 0, &mut pairs);
        pairs
    }
    /// Query overlapping face pairs, expanding each AABB by `margin` before testing.
    pub fn query_overlap_pairs_margin(
        &self,
        other: &DeformBvh,
        margin: f64,
    ) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return pairs;
        }
        Self::overlap_recurse_margin(&self.nodes, 0, &other.nodes, 0, margin, &mut pairs);
        pairs
    }
    fn overlap_recurse_margin(
        a_nodes: &[DeformBvhNode],
        a_idx: usize,
        b_nodes: &[DeformBvhNode],
        b_idx: usize,
        margin: f64,
        pairs: &mut Vec<(usize, usize)>,
    ) {
        if !a_nodes[a_idx]
            .aabb
            .overlaps_with_margin(&b_nodes[b_idx].aabb, margin)
        {
            return;
        }
        let a_leaf = a_nodes[a_idx].left == usize::MAX;
        let b_leaf = b_nodes[b_idx].left == usize::MAX;
        if a_leaf && b_leaf {
            for &fi in &a_nodes[a_idx].faces {
                for &fj in &b_nodes[b_idx].faces {
                    pairs.push((fi, fj));
                }
            }
            return;
        }
        if a_leaf
            || (!b_leaf && a_nodes[a_idx].aabb.surface_area() <= b_nodes[b_idx].aabb.surface_area())
        {
            Self::overlap_recurse_margin(
                a_nodes,
                a_idx,
                b_nodes,
                b_nodes[b_idx].left,
                margin,
                pairs,
            );
            Self::overlap_recurse_margin(
                a_nodes,
                a_idx,
                b_nodes,
                b_nodes[b_idx].right,
                margin,
                pairs,
            );
        } else {
            Self::overlap_recurse_margin(
                a_nodes,
                a_nodes[a_idx].left,
                b_nodes,
                b_idx,
                margin,
                pairs,
            );
            Self::overlap_recurse_margin(
                a_nodes,
                a_nodes[a_idx].right,
                b_nodes,
                b_idx,
                margin,
                pairs,
            );
        }
    }
    /// Recursive overlap traversal between two BVH trees.
    fn overlap_recurse(
        a_nodes: &[DeformBvhNode],
        a_idx: usize,
        b_nodes: &[DeformBvhNode],
        b_idx: usize,
        pairs: &mut Vec<(usize, usize)>,
    ) {
        if !a_nodes[a_idx].aabb.overlaps(&b_nodes[b_idx].aabb) {
            return;
        }
        let a_leaf = a_nodes[a_idx].left == usize::MAX;
        let b_leaf = b_nodes[b_idx].left == usize::MAX;
        if a_leaf && b_leaf {
            for &fi in &a_nodes[a_idx].faces {
                for &fj in &b_nodes[b_idx].faces {
                    pairs.push((fi, fj));
                }
            }
            return;
        }
        if a_leaf
            || (!b_leaf && a_nodes[a_idx].aabb.surface_area() <= b_nodes[b_idx].aabb.surface_area())
        {
            Self::overlap_recurse(a_nodes, a_idx, b_nodes, b_nodes[b_idx].left, pairs);
            Self::overlap_recurse(a_nodes, a_idx, b_nodes, b_nodes[b_idx].right, pairs);
        } else {
            Self::overlap_recurse(a_nodes, a_nodes[a_idx].left, b_nodes, b_idx, pairs);
            Self::overlap_recurse(a_nodes, a_nodes[a_idx].right, b_nodes, b_idx, pairs);
        }
    }
    /// Self-overlap query (avoids duplicate/same-face pairs).
    pub fn query_self_overlap(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if self.nodes.is_empty() {
            return pairs;
        }
        self.self_overlap_recurse(0, 0, &mut pairs);
        pairs
    }
    /// Recursive self-overlap helper.
    fn self_overlap_recurse(&self, a_idx: usize, b_idx: usize, pairs: &mut Vec<(usize, usize)>) {
        if a_idx > b_idx {
            return;
        }
        if !self.nodes[a_idx].aabb.overlaps(&self.nodes[b_idx].aabb) {
            return;
        }
        let a_leaf = self.nodes[a_idx].left == usize::MAX;
        let b_leaf = self.nodes[b_idx].left == usize::MAX;
        if a_leaf && b_leaf {
            for &fi in &self.nodes[a_idx].faces {
                for &fj in &self.nodes[b_idx].faces {
                    if fi < fj {
                        pairs.push((fi, fj));
                    }
                }
            }
            return;
        }
        if a_idx == b_idx {
            let l = self.nodes[a_idx].left;
            let r = self.nodes[a_idx].right;
            if l != usize::MAX {
                self.self_overlap_recurse(l, l, pairs);
                self.self_overlap_recurse(l, r, pairs);
                self.self_overlap_recurse(r, r, pairs);
            }
            return;
        }
        if a_leaf {
            self.self_overlap_recurse(a_idx, self.nodes[b_idx].left, pairs);
            self.self_overlap_recurse(a_idx, self.nodes[b_idx].right, pairs);
        } else {
            self.self_overlap_recurse(self.nodes[a_idx].left, b_idx, pairs);
            self.self_overlap_recurse(self.nodes[a_idx].right, b_idx, pairs);
        }
    }
}
