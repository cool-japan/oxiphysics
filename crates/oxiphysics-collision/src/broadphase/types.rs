//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use crate::types::CollisionPair;
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Real, Vec3};

/// O(n^2) brute force broad phase for reference/debugging.
pub struct BruteForceBroadPhase;
/// Object type classification for pair-count histograms.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ObjectType {
    /// Static immovable object.
    Static,
    /// Dynamic simulated object.
    Dynamic,
    /// Kinematic (script-driven) object.
    Kinematic,
}
/// A node in the dynamic BVH tree.
#[derive(Debug, Clone)]
pub(super) enum BvhNodeData {
    /// Internal node with two children.
    Internal { left: usize, right: usize },
    /// Leaf node storing an object index.
    Leaf { object_index: usize },
}
/// A hierarchical scene graph that accelerates broadphase by skipping
/// subtrees whose root AABB does not overlap the query.
#[derive(Debug, Default)]
#[allow(dead_code)]
pub struct BroadphaseSceneGraph {
    pub(super) nodes: Vec<SceneGraphNode>,
}
#[allow(dead_code)]
impl BroadphaseSceneGraph {
    /// Create an empty scene graph.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }
    /// Add a node. Returns its index.
    pub fn add_node(&mut self, aabb: Aabb, parent: Option<usize>) -> usize {
        let idx = self.nodes.len();
        if let Some(p) = parent
            && p < self.nodes.len()
        {
            self.nodes[p].children.push(idx);
        }
        self.nodes.push(SceneGraphNode {
            aabb,
            parent,
            children: Vec::new(),
            active: true,
        });
        idx
    }
    /// Remove a node (marks inactive, does not compact).
    pub fn remove_node(&mut self, idx: usize) {
        if idx < self.nodes.len() {
            self.nodes[idx].active = false;
        }
    }
    /// Update the AABB for a node.
    pub fn update_aabb(&mut self, idx: usize, aabb: Aabb) {
        if idx < self.nodes.len() {
            self.nodes[idx].aabb = aabb;
        }
    }
    /// Number of active nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.active).count()
    }
    /// Returns the indices of children of `parent`.
    pub fn children_of(&self, parent: usize) -> Vec<usize> {
        if parent < self.nodes.len() {
            self.nodes[parent]
                .children
                .iter()
                .copied()
                .filter(|&c| self.nodes.get(c).is_some_and(|n| n.active))
                .collect()
        } else {
            Vec::new()
        }
    }
    /// Collect all active AABBs with their original node indices.
    fn active_aabbs(&self) -> (Vec<Aabb>, Vec<usize>) {
        let mut aabbs = Vec::new();
        let mut indices = Vec::new();
        for (i, node) in self.nodes.iter().enumerate() {
            if node.active {
                aabbs.push(node.aabb.clone());
                indices.push(i);
            }
        }
        (aabbs, indices)
    }
    /// Find all overlapping pairs using a BVH.
    pub fn find_pairs_bvh(&self) -> Vec<CollisionPair> {
        let (aabbs, node_indices) = self.active_aabbs();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        let raw_pairs = bvh.find_all_pairs(&aabbs);
        raw_pairs
            .iter()
            .map(|p| {
                let a = node_indices[p.a];
                let b = node_indices[p.b];
                if a < b {
                    CollisionPair::new(a, b)
                } else {
                    CollisionPair::new(b, a)
                }
            })
            .collect()
    }
    /// Query all active nodes overlapping the given AABB.
    pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
        let (aabbs, node_indices) = self.active_aabbs();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&aabbs);
        bvh.query(query).iter().map(|&i| node_indices[i]).collect()
    }
}
/// Instruments a broadphase query to collect timing information.
#[derive(Debug, Default)]
#[allow(dead_code)]
pub struct BroadphaseProfiler {
    /// Elapsed nanoseconds for the last query.
    pub(super) last_elapsed_ns: u64,
    /// Total number of profiled calls.
    pub(super) call_count: u64,
    /// Total elapsed nanoseconds across all calls.
    pub(super) total_elapsed_ns: u64,
}
#[allow(dead_code)]
impl BroadphaseProfiler {
    /// Create a new profiler.
    pub fn new() -> Self {
        Self::default()
    }
    /// Run SAP broadphase and record timing.
    pub fn profile_sap(&mut self, aabbs: &[Aabb], axis: usize) -> Vec<CollisionPair> {
        let start = std::time::Instant::now();
        let result = SweepAndPrune::new(axis).find_pairs(aabbs);
        let elapsed = start.elapsed().as_nanos() as u64;
        self.last_elapsed_ns = elapsed;
        self.total_elapsed_ns += elapsed;
        self.call_count += 1;
        result
    }
    /// Run BVH broadphase and record timing.
    pub fn profile_bvh(&mut self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        let start = std::time::Instant::now();
        let mut bvh = BvhBroadphase::new();
        bvh.build(aabbs);
        let result = bvh.find_all_pairs(aabbs);
        let elapsed = start.elapsed().as_nanos() as u64;
        self.last_elapsed_ns = elapsed;
        self.total_elapsed_ns += elapsed;
        self.call_count += 1;
        result
    }
    /// Elapsed nanoseconds for the most recent profiled call.
    pub fn last_elapsed_ns(&self) -> u64 {
        self.last_elapsed_ns
    }
    /// Total number of calls profiled.
    pub fn call_count(&self) -> u64 {
        self.call_count
    }
    /// Average nanoseconds per call.
    pub fn avg_elapsed_ns(&self) -> f64 {
        if self.call_count == 0 {
            0.0
        } else {
            self.total_elapsed_ns as f64 / self.call_count as f64
        }
    }
    /// Reset all counters.
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}
/// A view frustum defined by 6 planes.
///
/// Each plane is stored as `(normal, offset)` where `normal · x = offset`
/// defines the plane boundary.  The normal points inward (toward the
/// visible region).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Frustum {
    /// Six frustum planes: `(inward_normal, offset)`.
    pub planes: [(Vec3, Real); 6],
}
#[allow(dead_code)]
impl Frustum {
    /// Create a frustum from six plane definitions.
    pub fn new(planes: [(Vec3, Real); 6]) -> Self {
        Self { planes }
    }
    /// Create a symmetric frustum given camera parameters.
    ///
    /// Camera looks along -Z.  Near plane at z = -near, far plane at z = -far.
    /// Inward normals point toward the visible region.
    pub fn from_perspective(fov_y_rad: Real, aspect: Real, near: Real, far: Real) -> Self {
        let half_fov = fov_y_rad * 0.5;
        let tan_half = half_fov.tan();
        let near_plane = (Vec3::new(0.0, 0.0, -1.0), near);
        let far_plane = (Vec3::new(0.0, 0.0, 1.0), -far);
        let cos_fov = (1.0 / (1.0 + tan_half * tan_half)).sqrt();
        let sin_fov = tan_half * cos_fov;
        let top_plane = (Vec3::new(0.0, -cos_fov, -sin_fov), 0.0);
        let bottom_plane = (Vec3::new(0.0, cos_fov, -sin_fov), 0.0);
        let tan_half_x = tan_half * aspect;
        let cos_fov_x = (1.0 / (1.0 + tan_half_x * tan_half_x)).sqrt();
        let sin_fov_x = tan_half_x * cos_fov_x;
        let left_plane = (Vec3::new(cos_fov_x, 0.0, -sin_fov_x), 0.0);
        let right_plane = (Vec3::new(-cos_fov_x, 0.0, -sin_fov_x), 0.0);
        Self {
            planes: [
                near_plane,
                far_plane,
                top_plane,
                bottom_plane,
                left_plane,
                right_plane,
            ],
        }
    }
    /// Test whether an AABB is (potentially) inside the frustum.
    ///
    /// Returns `true` if the AABB is not fully outside any plane.
    pub fn contains_aabb(&self, aabb: &Aabb) -> bool {
        for &(ref normal, offset) in &self.planes {
            let px = if normal.x >= 0.0 {
                aabb.max.x
            } else {
                aabb.min.x
            };
            let py = if normal.y >= 0.0 {
                aabb.max.y
            } else {
                aabb.min.y
            };
            let pz = if normal.z >= 0.0 {
                aabb.max.z
            } else {
                aabb.min.z
            };
            let p = Vec3::new(px, py, pz);
            if p.dot(normal) < offset {
                return false;
            }
        }
        true
    }
}
/// Incrementally updated AABB tree (insertion/removal without full rebuild).
///
/// Each object gets a fattened AABB; the tree is only updated when the
/// object moves outside its fat AABB.
#[allow(dead_code)]
pub struct DynamicAabbTree {
    /// Object AABBs (tight).
    pub(super) tight_aabbs: Vec<Aabb>,
    /// Fattened AABBs (tight + margin).
    pub(super) fat_aabbs: Vec<Aabb>,
    /// Fattening margin.
    pub(super) margin: Real,
    /// Internal BVH for queries (rebuilt on demand).
    pub(super) bvh: BvhBroadphase,
    /// Whether the tree needs rebuilding.
    pub(super) dirty: bool,
}
#[allow(dead_code)]
impl DynamicAabbTree {
    /// Create a new dynamic tree with the given fattening margin.
    pub fn new(margin: Real) -> Self {
        Self {
            tight_aabbs: Vec::new(),
            fat_aabbs: Vec::new(),
            margin,
            bvh: BvhBroadphase::new(),
            dirty: true,
        }
    }
    /// Insert a new object. Returns its index.
    pub fn insert(&mut self, aabb: Aabb) -> usize {
        let idx = self.tight_aabbs.len();
        let fat = self.fatten(&aabb);
        self.tight_aabbs.push(aabb);
        self.fat_aabbs.push(fat);
        self.dirty = true;
        idx
    }
    /// Update the AABB for an existing object.
    ///
    /// Returns `true` if the internal tree was invalidated.
    pub fn update(&mut self, index: usize, aabb: Aabb) -> bool {
        self.tight_aabbs[index] = aabb.clone();
        if !self.fat_aabbs[index].contains_aabb(&aabb) {
            self.fat_aabbs[index] = self.fatten(&aabb);
            self.dirty = true;
            return true;
        }
        false
    }
    /// Ensure the internal BVH is up to date.
    pub fn rebuild_if_dirty(&mut self) {
        if self.dirty {
            self.bvh.build(&self.fat_aabbs);
            self.dirty = false;
        }
    }
    /// Query all overlapping pairs.
    pub fn find_pairs(&mut self) -> Vec<CollisionPair> {
        self.rebuild_if_dirty();
        let mut pairs = Vec::new();
        for i in 0..self.tight_aabbs.len() {
            let hits = self.bvh.query(&self.tight_aabbs[i]);
            for j in hits {
                if j > i && self.tight_aabbs[i].intersects(&self.tight_aabbs[j]) {
                    pairs.push(CollisionPair::new(i, j));
                }
            }
        }
        pairs
    }
    /// Query all objects overlapping a given AABB.
    pub fn query(&mut self, aabb: &Aabb) -> Vec<usize> {
        self.rebuild_if_dirty();
        self.bvh.query(aabb)
    }
    /// Number of objects.
    pub fn len(&self) -> usize {
        self.tight_aabbs.len()
    }
    /// Whether the tree is empty.
    pub fn is_empty(&self) -> bool {
        self.tight_aabbs.is_empty()
    }
    fn fatten(&self, aabb: &Aabb) -> Aabb {
        let m = Vec3::new(self.margin, self.margin, self.margin);
        Aabb::new(aabb.min - m, aabb.max + m)
    }
}
#[allow(dead_code)]
impl DynamicAabbTree {
    /// Mark the tree as dirty so it will be rebuilt on the next query.
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }
    /// Returns `true` if the tree needs rebuilding.
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }
    /// Return all pairs whose AABBs overlap when *inflated* by `inflation`.
    ///
    /// This is a proximity query: objects within `inflation` of each other
    /// are included even if their tight AABBs do not touch.
    pub fn find_proximity_pairs(&mut self, inflation: Real) -> Vec<CollisionPair> {
        self.rebuild_if_dirty();
        let n = self.tight_aabbs.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            let inflated = inflate_aabb(&self.tight_aabbs[i], inflation);
            let hits = self.bvh.query(&inflated);
            for j in hits {
                if j > i {
                    let inflated_j = inflate_aabb(&self.tight_aabbs[j], inflation);
                    if inflated.intersects(&inflated_j) {
                        pairs.push(CollisionPair::new(i, j));
                    }
                }
            }
        }
        pairs
    }
    /// Batch-insert many AABBs at once (more efficient than repeated `insert`).
    ///
    /// Appends all AABBs and marks the tree dirty; call `rebuild_if_dirty`
    /// (or any query method) to make the tree consistent.
    pub fn batch_insert(&mut self, aabbs: &[Aabb]) {
        for aabb in aabbs {
            let fat = self.fatten(aabb);
            self.tight_aabbs.push(aabb.clone());
            self.fat_aabbs.push(fat);
        }
        self.dirty = true;
    }
    /// Validate parent/child link consistency in the underlying BVH.
    ///
    /// Calls `rebuild_if_dirty` first, then delegates to BvhBroadphase::validate.
    pub fn validate(&mut self) -> bool {
        self.rebuild_if_dirty();
        self.bvh.validate()
    }
    /// Return tree statistics: (total objects, is_dirty).
    pub fn tree_info(&self) -> (usize, bool) {
        (self.tight_aabbs.len(), self.dirty)
    }
    /// Remove all objects and reset the tree.
    pub fn clear(&mut self) {
        self.tight_aabbs.clear();
        self.fat_aabbs.clear();
        self.bvh = BvhBroadphase::new();
        self.dirty = false;
    }
}
/// Hints for GPU-accelerated broad-phase.
///
/// These values can be used to configure a GPU compute shader for
/// uniform-grid or sort-based broadphase.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GpuBroadphaseHints {
    /// Number of objects in the scene.
    pub n_objects: usize,
    /// Scene AABB minimum corner.
    pub scene_min: [Real; 3],
    /// Scene AABB maximum corner.
    pub scene_max: [Real; 3],
    /// Average AABB half-extent across all objects.
    pub avg_half_extent: Real,
    /// Maximum AABB half-extent.
    pub max_half_extent: Real,
    /// Recommended number of GPU threads.
    pub recommended_threads: usize,
}
#[allow(dead_code)]
impl GpuBroadphaseHints {
    /// Compute GPU hints from a slice of AABBs.
    pub fn from_aabbs(aabbs: &[Aabb]) -> Self {
        if aabbs.is_empty() {
            return Self {
                n_objects: 0,
                scene_min: [0.0; 3],
                scene_max: [0.0; 3],
                avg_half_extent: 0.0,
                max_half_extent: 0.0,
                recommended_threads: 64,
            };
        }
        let mut scene_min = [Real::MAX; 3];
        let mut scene_max = [Real::MIN; 3];
        let mut sum_extent = 0.0_f64;
        let mut max_extent = 0.0_f64;
        for aabb in aabbs {
            for i in 0..3 {
                if aabb.min[i] < scene_min[i] {
                    scene_min[i] = aabb.min[i];
                }
                if aabb.max[i] > scene_max[i] {
                    scene_max[i] = aabb.max[i];
                }
            }
            let half = aabb.half_extents();
            let extent = half.x.max(half.y).max(half.z);
            sum_extent += extent;
            if extent > max_extent {
                max_extent = extent;
            }
        }
        let avg_half_extent = sum_extent / aabbs.len() as f64;
        let recommended_threads = aabbs.len().next_power_of_two().max(64);
        Self {
            n_objects: aabbs.len(),
            scene_min,
            scene_max,
            avg_half_extent,
            max_half_extent: max_extent,
            recommended_threads,
        }
    }
    /// Suggested uniform-grid cell size (2× average extent for good fill).
    pub fn suggested_cell_size(&self) -> Real {
        if self.n_objects == 0 {
            return 1.0;
        }
        (self.avg_half_extent * 2.0).max(1e-6)
    }
    /// Estimated number of non-empty grid cells.
    pub fn estimated_cells(&self) -> usize {
        if self.n_objects == 0 {
            return 0;
        }
        let cell = self.suggested_cell_size();
        let mut count = 1usize;
        for i in 0..3 {
            let extent = (self.scene_max[i] - self.scene_min[i]).max(0.0);
            count *= ((extent / cell).ceil() as usize).max(1);
        }
        count
    }
    /// Return flat array of AABB data suitable for GPU upload:
    /// `[min_x, min_y, min_z, max_x, max_y, max_z]` per object.
    pub fn to_flat_buffer(aabbs: &[Aabb]) -> Vec<f32> {
        let mut buf = Vec::with_capacity(aabbs.len() * 6);
        for aabb in aabbs {
            buf.push(aabb.min.x as f32);
            buf.push(aabb.min.y as f32);
            buf.push(aabb.min.z as f32);
            buf.push(aabb.max.x as f32);
            buf.push(aabb.max.y as f32);
            buf.push(aabb.max.z as f32);
        }
        buf
    }
}
/// Quality metrics for a static BVH.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct BvhQualityMetrics {
    /// SAH cost normalised by root surface area.
    pub sah_cost: f64,
    /// Maximum depth of any leaf.
    pub max_depth: usize,
    /// Average depth of leaves.
    pub avg_leaf_depth: f64,
    /// Number of leaf nodes.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
}
/// A node in the dynamic BVH tree.
#[derive(Debug, Clone)]
pub(super) struct BvhNode {
    pub(super) aabb: Aabb,
    pub(super) parent: Option<usize>,
    pub(super) data: BvhNodeData,
}
/// Statistics from a broad-phase query.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct BroadphaseStats {
    /// Number of objects.
    pub num_objects: usize,
    /// Number of candidate pairs found.
    pub num_pairs: usize,
    /// Number of AABB-AABB tests performed (estimated).
    pub num_tests: usize,
}
/// A node in the broadphase scene graph.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SceneGraphNode {
    /// World-space AABB for this node.
    pub aabb: Aabb,
    /// Optional parent node index.
    pub parent: Option<usize>,
    /// Child node indices.
    pub children: Vec<usize>,
    /// Whether this node is active (not removed).
    pub active: bool,
}
/// Histogram entry: how many collision pairs exist between two object types.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct PairCountHistogram {
    /// Static vs Static pairs.
    pub static_static: usize,
    /// Static vs Dynamic pairs.
    pub static_dynamic: usize,
    /// Static vs Kinematic pairs.
    pub static_kinematic: usize,
    /// Dynamic vs Dynamic pairs.
    pub dynamic_dynamic: usize,
    /// Dynamic vs Kinematic pairs.
    pub dynamic_kinematic: usize,
    /// Kinematic vs Kinematic pairs.
    pub kinematic_kinematic: usize,
}
impl PairCountHistogram {
    /// Total pair count across all categories.
    #[allow(dead_code)]
    pub fn total(&self) -> usize {
        self.static_static
            + self.static_dynamic
            + self.static_kinematic
            + self.dynamic_dynamic
            + self.dynamic_kinematic
            + self.kinematic_kinematic
    }
}
/// Top-down Bounding Volume Hierarchy for efficient broad-phase queries.
///
/// This implementation uses a static top-down build and is best used by
/// calling [`BvhBroadphase::build`] once per frame or after bulk updates.
/// For incremental updates, see [`crate::DynamicBvh`].
pub struct BvhBroadphase {
    pub(super) nodes: Vec<BvhNode>,
    pub(super) root: Option<usize>,
    /// Maps object index -> node index in the tree.
    pub(super) object_to_node: Vec<Option<usize>>,
}
impl BvhBroadphase {
    /// Create a new empty BVH.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: None,
            object_to_node: Vec::new(),
        }
    }
    /// Build the BVH from a set of AABBs using top-down construction.
    pub fn build(&mut self, aabbs: &[Aabb]) {
        self.nodes.clear();
        self.object_to_node.clear();
        self.object_to_node.resize(aabbs.len(), None);
        if aabbs.is_empty() {
            self.root = None;
            return;
        }
        let indices: Vec<usize> = (0..aabbs.len()).collect();
        self.root = Some(self.build_recursive(aabbs, &indices));
    }
    fn build_recursive(&mut self, aabbs: &[Aabb], indices: &[usize]) -> usize {
        if indices.len() == 1 {
            let idx = indices[0];
            let node_idx = self.nodes.len();
            self.nodes.push(BvhNode {
                aabb: aabbs[idx].clone(),
                parent: None,
                data: BvhNodeData::Leaf { object_index: idx },
            });
            self.object_to_node[idx] = Some(node_idx);
            return node_idx;
        }
        let mut combined = aabbs[indices[0]].clone();
        for &i in &indices[1..] {
            combined = combined.merge(&aabbs[i]);
        }
        let extent = combined.half_extents();
        let split_axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let mut sorted = indices.to_vec();
        sorted.sort_by(|&a, &b| {
            aabbs[a].center()[split_axis]
                .partial_cmp(&aabbs[b].center()[split_axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted.len() / 2;
        let left = self.build_recursive(aabbs, &sorted[..mid]);
        let right = self.build_recursive(aabbs, &sorted[mid..]);
        let node_idx = self.nodes.len();
        self.nodes.push(BvhNode {
            aabb: combined,
            parent: None,
            data: BvhNodeData::Internal { left, right },
        });
        self.nodes[left].parent = Some(node_idx);
        self.nodes[right].parent = Some(node_idx);
        node_idx
    }
    /// Query all pairs of overlapping AABBs.
    pub fn find_all_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        let mut pairs = Vec::new();
        if self.root.is_none() {
            return pairs;
        }
        let leaves: Vec<(usize, usize)> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(node_idx, node)| {
                if let BvhNodeData::Leaf { object_index } = node.data {
                    Some((node_idx, object_index))
                } else {
                    None
                }
            })
            .collect();
        for (i, &(_, obj_a)) in leaves.iter().enumerate() {
            for &(_, obj_b) in &leaves[(i + 1)..] {
                if aabbs[obj_a].intersects(&aabbs[obj_b]) {
                    pairs.push(CollisionPair::new(obj_a, obj_b));
                }
            }
        }
        pairs
    }
    /// Query all objects whose AABB overlaps the given AABB.
    pub fn query(&self, aabb: &Aabb) -> Vec<usize> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_recursive(root, aabb, &mut results);
        }
        results
    }
    fn query_recursive(&self, node_idx: usize, aabb: &Aabb, results: &mut Vec<usize>) {
        let node = &self.nodes[node_idx];
        if !node.aabb.intersects(aabb) {
            return;
        }
        match node.data {
            BvhNodeData::Leaf { object_index } => {
                results.push(object_index);
            }
            BvhNodeData::Internal { left, right } => {
                self.query_recursive(left, aabb, results);
                self.query_recursive(right, aabb, results);
            }
        }
    }
    /// Ray query: find all objects whose AABB is hit by the ray.
    pub fn ray_query(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Vec<usize> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.ray_query_recursive(root, ray_origin, ray_direction, max_toi, &mut results);
        }
        results
    }
    fn ray_query_recursive(
        &self,
        node_idx: usize,
        ray_origin: &Vec3,
        ray_direction: &Vec3,
        max_toi: Real,
        results: &mut Vec<usize>,
    ) {
        let node = &self.nodes[node_idx];
        if !ray_aabb_intersect(ray_origin, ray_direction, &node.aabb, max_toi) {
            return;
        }
        match node.data {
            BvhNodeData::Leaf { object_index } => {
                results.push(object_index);
            }
            BvhNodeData::Internal { left, right } => {
                self.ray_query_recursive(left, ray_origin, ray_direction, max_toi, results);
                self.ray_query_recursive(right, ray_origin, ray_direction, max_toi, results);
            }
        }
    }
}
#[allow(dead_code)]
impl BvhBroadphase {
    /// Compute quality metrics for this BVH.
    ///
    /// Returns `None` if the tree is empty.
    pub fn quality_metrics(&self) -> Option<BvhQualityMetrics> {
        let root = self.root?;
        let root_sa = self.nodes[root].aabb.surface_area();
        if root_sa <= 0.0 {
            return None;
        }
        let mut sah = 0.0_f64;
        let mut max_depth = 0usize;
        let mut total_leaf_depth = 0usize;
        let mut leaf_count = 0usize;
        let mut internal_count = 0usize;
        let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
        while let Some((idx, depth)) = stack.pop() {
            if depth > max_depth {
                max_depth = depth;
            }
            match self.nodes[idx].data {
                BvhNodeData::Leaf { .. } => {
                    leaf_count += 1;
                    total_leaf_depth += depth;
                }
                BvhNodeData::Internal { left, right } => {
                    internal_count += 1;
                    sah += self.nodes[idx].aabb.surface_area() / root_sa;
                    stack.push((left, depth + 1));
                    stack.push((right, depth + 1));
                }
            }
        }
        let avg_leaf_depth = if leaf_count > 0 {
            total_leaf_depth as f64 / leaf_count as f64
        } else {
            0.0
        };
        Some(BvhQualityMetrics {
            sah_cost: sah,
            max_depth,
            avg_leaf_depth,
            leaf_count,
            internal_count,
        })
    }
    /// Rebuild the BVH from `aabbs` if the SAH cost exceeds `threshold`.
    ///
    /// Returns `true` if a rebuild was triggered.
    pub fn rebuild_if_degraded(&mut self, aabbs: &[Aabb], threshold: f64) -> bool {
        let degraded = match self.quality_metrics() {
            Some(q) => q.sah_cost > threshold,
            None => false,
        };
        if degraded {
            self.build(aabbs);
        }
        degraded
    }
    /// SAH cost normalised by root surface area, or `0.0` for an empty tree.
    pub fn sah_cost(&self) -> f64 {
        self.quality_metrics().map(|q| q.sah_cost).unwrap_or(0.0)
    }
    /// Validate all parent/child links and AABB containment.
    ///
    /// Returns `true` if the tree is consistent.
    pub fn validate(&self) -> bool {
        match self.root {
            None => true,
            Some(root) => {
                if self.nodes[root].parent.is_some() {
                    return false;
                }
                self.validate_node(root)
            }
        }
    }
    fn validate_node(&self, idx: usize) -> bool {
        match self.nodes[idx].data {
            BvhNodeData::Leaf { .. } => true,
            BvhNodeData::Internal { left, right } => {
                if self.nodes[left].parent != Some(idx) {
                    return false;
                }
                if self.nodes[right].parent != Some(idx) {
                    return false;
                }
                if !self.nodes[idx].aabb.contains_aabb(&self.nodes[left].aabb) {
                    return false;
                }
                if !self.nodes[idx].aabb.contains_aabb(&self.nodes[right].aabb) {
                    return false;
                }
                self.validate_node(left) && self.validate_node(right)
            }
        }
    }
    /// Stack-based ray traversal against this BVH; returns hit object indices.
    ///
    /// Equivalent to `ray_query` but uses an explicit stack instead of recursion.
    pub fn ray_query_iterative(
        &self,
        ray_origin: &Vec3,
        ray_direction: &Vec3,
        max_toi: Real,
    ) -> Vec<usize> {
        let mut results = Vec::new();
        let Some(root) = self.root else {
            return results;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !ray_aabb_intersect(ray_origin, ray_direction, &node.aabb, max_toi) {
                continue;
            }
            match node.data {
                BvhNodeData::Leaf { object_index } => results.push(object_index),
                BvhNodeData::Internal { left, right } => {
                    stack.push(left);
                    stack.push(right);
                }
            }
        }
        results
    }
    /// Frustum culling against this BVH using a stack-based traversal.
    ///
    /// `planes`: six `(inward_normal, offset)` pairs.
    /// Returns object indices of leaves that are potentially inside the frustum.
    pub fn frustum_query_iterative(&self, planes: &[(Vec3, Real); 6]) -> Vec<usize> {
        let mut results = Vec::new();
        let Some(root) = self.root else {
            return results;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !aabb_in_frustum(&node.aabb, planes) {
                continue;
            }
            match node.data {
                BvhNodeData::Leaf { object_index } => results.push(object_index),
                BvhNodeData::Internal { left, right } => {
                    stack.push(left);
                    stack.push(right);
                }
            }
        }
        results
    }
    /// Batch insert a slice of AABBs and rebuild.
    ///
    /// Clears any existing tree and builds from the provided `aabbs`.
    pub fn batch_insert(&mut self, aabbs: &[Aabb]) {
        self.build(aabbs);
    }
}
/// Sweep and Prune (SAP) broad phase on a single axis.
pub struct SweepAndPrune {
    pub(super) axis: usize,
}
impl SweepAndPrune {
    /// Create a new SAP on the given axis (0=X, 1=Y, 2=Z).
    pub fn new(axis: usize) -> Self {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        Self { axis }
    }
    /// Create SAP on X axis.
    pub fn x_axis() -> Self {
        Self::new(0)
    }
}
/// Cache of collision pairs from the previous frame to speed up the
/// next broadphase query (warm starting).
///
/// On each frame:
/// 1. Filter cached pairs to those still overlapping.
/// 2. Merge with newly detected pairs.
/// 3. Update the cache with the merged set.
#[derive(Debug, Default)]
#[allow(dead_code)]
pub struct BroadphaseWarmstart {
    pub(super) cached_pairs: Vec<CollisionPair>,
}
#[allow(dead_code)]
impl BroadphaseWarmstart {
    /// Create a new empty warmstart cache.
    pub fn new() -> Self {
        Self::default()
    }
    /// Update the cache with the pairs from the current frame.
    pub fn update(&mut self, pairs: &[CollisionPair]) {
        self.cached_pairs = pairs.to_vec();
    }
    /// Number of pairs in the cache.
    pub fn cached_pair_count(&self) -> usize {
        self.cached_pairs.len()
    }
    /// Filter cached pairs, retaining only those whose AABBs still overlap.
    pub fn filter_still_overlapping(&self, aabbs: &[Aabb]) -> Vec<CollisionPair> {
        self.cached_pairs
            .iter()
            .filter(|p| {
                let max_idx = aabbs.len();
                p.a < max_idx && p.b < max_idx && aabbs[p.a].intersects(&aabbs[p.b])
            })
            .copied()
            .collect()
    }
    /// Merge cached pairs that are still overlapping with fresh new pairs.
    ///
    /// Deduplicates the output.
    pub fn merge_with_new(
        &self,
        new_pairs: &[CollisionPair],
        aabbs: &[Aabb],
    ) -> Vec<CollisionPair> {
        let mut merged: Vec<CollisionPair> = self.filter_still_overlapping(aabbs);
        for &p in new_pairs {
            if !merged.iter().any(|e| e.a == p.a && e.b == p.b) {
                merged.push(p);
            }
        }
        merged
    }
    /// Clear the cache.
    pub fn clear(&mut self) {
        self.cached_pairs.clear();
    }
    /// Run a full broadphase using SAP with warmstart.
    ///
    /// Returns the merged pair set (warm + fresh) and updates the cache.
    pub fn broadphase_with_warmstart(&mut self, aabbs: &[Aabb], axis: usize) -> Vec<CollisionPair> {
        let fresh = SweepAndPrune::new(axis).find_pairs(aabbs);
        let merged = self.merge_with_new(&fresh, aabbs);
        self.update(&merged);
        merged
    }
}
