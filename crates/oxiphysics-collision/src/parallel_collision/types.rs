//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock};

use super::functions::NodeIdx;
#[allow(unused_imports)]
use super::functions::*;

/// Thread-safe union-find for parallel island detection.
///
/// Uses path-compression and union-by-rank. Multiple threads may call
/// `find` concurrently; `union` takes an exclusive lock.
#[derive(Debug)]
pub struct ParallelUnionFind {
    pub(super) parent: Mutex<Vec<usize>>,
    pub(super) rank: Mutex<Vec<usize>>,
    pub(super) n: usize,
}
impl ParallelUnionFind {
    /// Create a new union-find for `n` elements.
    pub fn new(n: usize) -> Self {
        Self {
            parent: Mutex::new((0..n).collect()),
            rank: Mutex::new(vec![0; n]),
            n,
        }
    }
    /// Find the root of the set containing `x` with path compression.
    pub fn find(&self, x: usize) -> usize {
        let mut parent = self.parent.lock().unwrap_or_else(|e| e.into_inner());
        let mut root = x;
        while parent[root] != root {
            root = parent[root];
        }
        let mut node = x;
        while parent[node] != root {
            let next = parent[node];
            parent[node] = root;
            node = next;
        }
        root
    }
    /// Union the sets containing `x` and `y`.
    pub fn union(&self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        let mut parent = self.parent.lock().unwrap_or_else(|e| e.into_inner());
        let mut rank = self.rank.lock().unwrap_or_else(|e| e.into_inner());
        if rank[rx] < rank[ry] {
            parent[rx] = ry;
        } else if rank[rx] > rank[ry] {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx] += 1;
        }
    }
    /// Collect all islands as groups of element indices.
    pub fn islands(&self) -> Vec<Vec<usize>> {
        let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..self.n {
            let root = self.find(i);
            map.entry(root).or_default().push(i);
        }
        map.into_values().collect()
    }
    /// Number of distinct connected components.
    pub fn num_components(&self) -> usize {
        let mut roots = std::collections::HashSet::new();
        for i in 0..self.n {
            roots.insert(self.find(i));
        }
        roots.len()
    }
}
/// Axis-aligned bounding box using plain `[f64; 3]` components.
///
/// Designed for SIMD-friendly packing in structure-of-arrays layouts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParAabb {
    /// Minimum corner `[x_min, y_min, z_min]`.
    pub min: [f64; 3],
    /// Maximum corner `[x_max, y_max, z_max]`.
    pub max: [f64; 3],
}
impl ParAabb {
    /// Create a new AABB from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Create an AABB from a center point and half-extents.
    pub fn from_center_half_extents(center: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            min: sub3(center, half),
            max: add3(center, half),
        }
    }
    /// Merge two AABBs returning the smallest enclosing AABB.
    pub fn merge(&self, other: &Self) -> Self {
        Self {
            min: min3(self.min, other.min),
            max: max3(self.max, other.max),
        }
    }
    /// Return the center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        scale3(add3(self.min, self.max), 0.5)
    }
    /// Return the half-extents of the AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        scale3(sub3(self.max, self.min), 0.5)
    }
    /// Surface area of the AABB (for SAH heuristic).
    pub fn surface_area(&self) -> f64 {
        let d = sub3(self.max, self.min);
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }
    /// Test whether two AABBs overlap.
    ///
    /// Vectorized-friendly: checks all axes without branching.
    #[inline(always)]
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Expand the AABB by `margin` on all sides.
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }
    /// Check if this AABB contains a point.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
}
/// Result of a CCD query between two swept spheres.
#[derive(Debug, Clone, Copy)]
pub struct CcdHit {
    /// Time of impact in `[0, 1]`.
    pub toi: f64,
    /// Index of body A.
    pub body_a: u32,
    /// Index of body B.
    pub body_b: u32,
    /// Contact normal at the time of impact.
    pub normal: [f64; 3],
}
/// Descriptor for a GPU-mapped AABB buffer.
///
/// Holds layout information for a flat byte buffer that can be directly
/// uploaded to a compute shader. All fields are `u32` for alignment.
#[derive(Debug, Clone, Copy)]
pub struct GpuAabbBufferDesc {
    /// Byte offset of the min_x array.
    pub offset_min_x: u32,
    /// Byte offset of the min_y array.
    pub offset_min_y: u32,
    /// Byte offset of the min_z array.
    pub offset_min_z: u32,
    /// Byte offset of the max_x array.
    pub offset_max_x: u32,
    /// Byte offset of the max_y array.
    pub offset_max_y: u32,
    /// Byte offset of the max_z array.
    pub offset_max_z: u32,
    /// Number of AABBs in the buffer.
    pub count: u32,
    /// Stride between elements in bytes (typically 4 for f32).
    pub stride: u32,
}
impl GpuAabbBufferDesc {
    /// Compute descriptor for a tightly-packed f32 SoA buffer of `count` AABBs.
    pub fn new_f32_soa(count: u32) -> Self {
        let stride = 4u32;
        let array_bytes = count * stride;
        Self {
            offset_min_x: 0,
            offset_min_y: array_bytes,
            offset_min_z: array_bytes * 2,
            offset_max_x: array_bytes * 3,
            offset_max_y: array_bytes * 4,
            offset_max_z: array_bytes * 5,
            count,
            stride,
        }
    }
    /// Total byte size of the GPU buffer.
    pub fn total_bytes(&self) -> u32 {
        self.count * self.stride * 6
    }
    /// Serialize a `SoaAabbBuffer` to a flat `Vec`u8` for GPU upload.
    pub fn serialize_soa(buf: &SoaAabbBuffer) -> Vec<u8> {
        let n = buf.len();
        let mut out = Vec::with_capacity(n * 6 * 4);
        let write_slice = |dst: &mut Vec<u8>, src: &[f32]| {
            for &v in src {
                dst.extend_from_slice(&v.to_le_bytes());
            }
        };
        write_slice(&mut out, &buf.min_x);
        write_slice(&mut out, &buf.min_y);
        write_slice(&mut out, &buf.min_z);
        write_slice(&mut out, &buf.max_x);
        write_slice(&mut out, &buf.max_y);
        write_slice(&mut out, &buf.max_z);
        out
    }
}
#[derive(Debug, Clone)]
pub(super) struct AabbTreeInner {
    pub(super) nodes: Vec<ParAabbNode>,
    pub(super) root: NodeIdx,
    pub(super) free_list: Vec<NodeIdx>,
}
impl AabbTreeInner {
    fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: u32::MAX,
            free_list: Vec::new(),
        }
    }
    fn alloc_node(&mut self) -> NodeIdx {
        if let Some(idx) = self.free_list.pop() {
            idx
        } else {
            let idx = self.nodes.len() as NodeIdx;
            self.nodes.push(ParAabbNode {
                aabb: ParAabb::default(),
                left: u32::MAX,
                right: u32::MAX,
                handle: u32::MAX,
                height: 0,
            });
            idx
        }
    }
    fn free_node(&mut self, idx: NodeIdx) {
        self.free_list.push(idx);
    }
    fn fix_aabb(&mut self, idx: NodeIdx) {
        let left = self.nodes[idx as usize].left;
        let right = self.nodes[idx as usize].right;
        if left != u32::MAX && right != u32::MAX {
            let la = self.nodes[left as usize].aabb;
            let ra = self.nodes[right as usize].aabb;
            self.nodes[idx as usize].aabb = la.merge(&ra);
            let lh = self.nodes[left as usize].height;
            let rh = self.nodes[right as usize].height;
            self.nodes[idx as usize].height = 1 + lh.max(rh);
        }
    }
    /// Insert a leaf with the given AABB and body handle.
    fn insert(&mut self, aabb: ParAabb, handle: u32) -> NodeIdx {
        let leaf = self.alloc_node();
        self.nodes[leaf as usize].aabb = aabb;
        self.nodes[leaf as usize].handle = handle;
        self.nodes[leaf as usize].left = u32::MAX;
        self.nodes[leaf as usize].right = u32::MAX;
        self.nodes[leaf as usize].height = 0;
        if self.root == u32::MAX {
            self.root = leaf;
            return leaf;
        }
        let mut best = self.root;
        loop {
            let node = &self.nodes[best as usize];
            if node.is_leaf() {
                break;
            }
            let merged_area = node.aabb.merge(&aabb).surface_area();
            let left = node.left;
            let right = node.right;
            let la = self.nodes[left as usize].aabb.merge(&aabb).surface_area();
            let ra = self.nodes[right as usize].aabb.merge(&aabb).surface_area();
            if merged_area < la.min(ra) {
                break;
            }
            best = if la < ra { left } else { right };
        }
        let old_parent = self.alloc_node();
        let sib_aabb = self.nodes[best as usize].aabb;
        self.nodes[old_parent as usize].aabb = sib_aabb.merge(&aabb);
        self.nodes[old_parent as usize].left = best;
        self.nodes[old_parent as usize].right = leaf;
        let sib_h = self.nodes[best as usize].height;
        self.nodes[old_parent as usize].height = sib_h + 1;
        self.nodes[old_parent as usize].handle = u32::MAX;
        if self.root == best {
            self.root = old_parent;
        }
        leaf
    }
    /// Query all leaf handles whose AABB overlaps `query`.
    fn query_overlaps(&self, query: &ParAabb) -> Vec<u32> {
        let mut result = Vec::new();
        if self.root == u32::MAX {
            return result;
        }
        let mut stack = vec![self.root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx as usize];
            if !node.aabb.overlaps(query) {
                continue;
            }
            if node.is_leaf() {
                result.push(node.handle);
            } else {
                if node.left != u32::MAX {
                    stack.push(node.left);
                }
                if node.right != u32::MAX {
                    stack.push(node.right);
                }
            }
        }
        result
    }
}
/// Discriminant for parallel narrowphase shape kinds.
#[derive(Debug, Clone, Copy)]
pub enum ParShapeKind {
    /// Sphere shape.
    Sphere(ParSphere),
    /// Axis-aligned box shape.
    Box(ParBox),
}
/// State of a body for CCD sweep.
#[derive(Debug, Clone, Copy)]
pub struct CcdBodyState {
    /// Position at the start of the time step.
    pub pos0: [f64; 3],
    /// Position at the end of the time step.
    pub pos1: [f64; 3],
    /// Radius (for sphere-based CCD).
    pub radius: f64,
    /// Body index.
    pub index: u32,
}
impl CcdBodyState {
    /// Create a new CCD body state.
    pub fn new(pos0: [f64; 3], pos1: [f64; 3], radius: f64, index: u32) -> Self {
        Self {
            pos0,
            pos1,
            radius,
            index,
        }
    }
    /// Compute the swept AABB for CCD broadphase culling.
    pub fn swept_aabb(&self) -> ParAabb {
        let margin = [self.radius; 3];
        let mn0 = sub3(self.pos0, margin);
        let mx0 = add3(self.pos0, margin);
        let mn1 = sub3(self.pos1, margin);
        let mx1 = add3(self.pos1, margin);
        ParAabb {
            min: min3(mn0, mn1),
            max: max3(mx0, mx1),
        }
    }
}
/// GPU-uploadable BVH tree serialized as a flat array of `GpuBvhNode`.
#[derive(Debug, Clone)]
pub struct GpuBvhBuffer {
    /// Flat node array.
    pub nodes: Vec<GpuBvhNode>,
}
impl GpuBvhBuffer {
    /// Create an empty GPU BVH buffer.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }
    /// Build a simple linear BVH from a list of AABBs and handles.
    ///
    /// Uses median-split on the longest axis (SAH is deferred to GPU).
    pub fn build(aabbs: &[(ParAabb, u32)]) -> Self {
        let mut buf = Self::new();
        if aabbs.is_empty() {
            return buf;
        }
        Self::build_recursive(&mut buf.nodes, aabbs);
        buf
    }
    fn build_recursive(nodes: &mut Vec<GpuBvhNode>, aabbs: &[(ParAabb, u32)]) -> u32 {
        if aabbs.len() == 1 {
            let (aabb, handle) = &aabbs[0];
            let node = GpuBvhNode::leaf(
                [aabb.min[0] as f32, aabb.min[1] as f32, aabb.min[2] as f32],
                [aabb.max[0] as f32, aabb.max[1] as f32, aabb.max[2] as f32],
                *handle,
            );
            let idx = nodes.len() as u32;
            nodes.push(node);
            return idx;
        }
        let mut combined = aabbs[0].0;
        for (a, _) in aabbs.iter().skip(1) {
            combined = combined.merge(a);
        }
        let extent = sub3(combined.max, combined.min);
        let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
            0
        } else if extent[1] >= extent[2] {
            1
        } else {
            2
        };
        let mut sorted: Vec<_> = aabbs.to_vec();
        sorted.sort_by(|(a, _), (b, _)| {
            let ca = (a.min[axis] + a.max[axis]) * 0.5;
            let cb = (b.min[axis] + b.max[axis]) * 0.5;
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted.len() / 2;
        let my_idx = nodes.len() as u32;
        nodes.push(GpuBvhNode::internal([0.0; 3], [0.0; 3], 0, 0));
        let left = Self::build_recursive(nodes, &sorted[..mid]);
        let right = Self::build_recursive(nodes, &sorted[mid..]);
        nodes[my_idx as usize] = GpuBvhNode::internal(
            [
                combined.min[0] as f32,
                combined.min[1] as f32,
                combined.min[2] as f32,
            ],
            [
                combined.max[0] as f32,
                combined.max[1] as f32,
                combined.max[2] as f32,
            ],
            left,
            right,
        );
        my_idx
    }
    /// Number of nodes in the BVH.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Serialize to bytes for GPU upload.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.nodes.len() * 32);
        for node in &self.nodes {
            for &v in &node.min {
                out.extend_from_slice(&v.to_le_bytes());
            }
            out.extend_from_slice(&node.left_or_leaf.to_le_bytes());
            for &v in &node.max {
                out.extend_from_slice(&v.to_le_bytes());
            }
            out.extend_from_slice(&node.right_or_handle.to_le_bytes());
        }
        out
    }
}
/// A contact point produced by narrowphase.
#[derive(Debug, Clone, Copy)]
pub struct ParContact {
    /// Index of body A.
    pub body_a: u32,
    /// Index of body B.
    pub body_b: u32,
    /// Contact normal pointing from B toward A.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
    /// Contact point in world space.
    pub point: [f64; 3],
}
/// Frame-level parallel collision result.
#[derive(Debug, Clone)]
pub struct ParCollisionResult {
    /// All reduced contacts.
    pub contacts: Vec<ParContact>,
    /// Simulation islands (groups of touching bodies).
    pub islands: Vec<Vec<usize>>,
    /// CCD hits this frame (empty if CCD disabled).
    pub ccd_hits: Vec<CcdHit>,
}
/// A box shape (OBB) for narrowphase tests, represented by center + half-extents.
///
/// For simplicity the box is axis-aligned (AABB) in this narrowphase.
#[derive(Debug, Clone, Copy)]
pub struct ParBox {
    /// Center of the box in world space.
    pub center: [f64; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}
/// A sphere shape for narrowphase tests.
#[derive(Debug, Clone, Copy)]
pub struct ParSphere {
    /// Center of the sphere in world space.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}
/// A node in the parallel AABB tree.
#[derive(Debug, Clone)]
pub struct ParAabbNode {
    /// Tight AABB enclosing this node's subtree.
    pub aabb: ParAabb,
    /// Index of the left child (`NodeIdx::MAX` = leaf).
    pub left: NodeIdx,
    /// Index of the right child (`NodeIdx::MAX` = leaf).
    pub right: NodeIdx,
    /// Body handle if this is a leaf node.
    pub handle: u32,
    /// Height in the tree (0 = leaf).
    pub height: i32,
}
impl ParAabbNode {
    /// Returns `true` if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left == u32::MAX
    }
}
/// Structure-of-Arrays layout for AABB data.
///
/// Stores min/max for each axis in separate arrays, enabling SIMD-wide
/// loads. Suitable for vectorized broadphase queries.
#[derive(Debug, Clone)]
pub struct SoaAabbBuffer {
    /// x-component of min corners.
    pub min_x: Vec<f32>,
    /// y-component of min corners.
    pub min_y: Vec<f32>,
    /// z-component of min corners.
    pub min_z: Vec<f32>,
    /// x-component of max corners.
    pub max_x: Vec<f32>,
    /// y-component of max corners.
    pub max_y: Vec<f32>,
    /// z-component of max corners.
    pub max_z: Vec<f32>,
    /// Body handle associated with each AABB slot.
    pub handles: Vec<u32>,
}
impl SoaAabbBuffer {
    /// Create an empty SoA AABB buffer.
    pub fn new() -> Self {
        Self {
            min_x: Vec::new(),
            min_y: Vec::new(),
            min_z: Vec::new(),
            max_x: Vec::new(),
            max_y: Vec::new(),
            max_z: Vec::new(),
            handles: Vec::new(),
        }
    }
    /// Push an AABB into the buffer.
    pub fn push(&mut self, aabb: &ParAabb, handle: u32) {
        self.min_x.push(aabb.min[0] as f32);
        self.min_y.push(aabb.min[1] as f32);
        self.min_z.push(aabb.min[2] as f32);
        self.max_x.push(aabb.max[0] as f32);
        self.max_y.push(aabb.max[1] as f32);
        self.max_z.push(aabb.max[2] as f32);
        self.handles.push(handle);
    }
    /// Number of AABBs stored.
    pub fn len(&self) -> usize {
        self.handles.len()
    }
    /// Returns `true` if the buffer contains no AABBs.
    pub fn is_empty(&self) -> bool {
        self.handles.is_empty()
    }
    /// Query a single AABB against the buffer, returning matching indices.
    pub fn query(&self, aabb: &ParAabb) -> Vec<usize> {
        let qminx = aabb.min[0] as f32;
        let qminy = aabb.min[1] as f32;
        let qminz = aabb.min[2] as f32;
        let qmaxx = aabb.max[0] as f32;
        let qmaxy = aabb.max[1] as f32;
        let qmaxz = aabb.max[2] as f32;
        let mut result = Vec::new();
        let n = self.handles.len();
        for i in 0..n {
            let hit = (qminx <= self.max_x[i])
                & (qmaxx >= self.min_x[i])
                & (qminy <= self.max_y[i])
                & (qmaxy >= self.min_y[i])
                & (qminz <= self.max_z[i])
                & (qmaxz >= self.min_z[i]);
            if hit {
                result.push(i);
            }
        }
        result
    }
    /// Rebuild the buffer from a slice of (aabb, handle) pairs.
    pub fn rebuild(&mut self, entries: &[(ParAabb, u32)]) {
        self.min_x.clear();
        self.min_y.clear();
        self.min_z.clear();
        self.max_x.clear();
        self.max_y.clear();
        self.max_z.clear();
        self.handles.clear();
        for (aabb, handle) in entries {
            self.push(aabb, *handle);
        }
    }
}
/// Work-stealing broadphase queue.
///
/// Divides candidate pair generation into chunks that can be stolen by
/// worker threads when their local queues become empty.
#[derive(Debug)]
pub struct WorkStealingBroadphase {
    pub(super) aabbs: Arc<Vec<ParAabb>>,
    pub(super) chunk_size: usize,
    pub(super) pending: Mutex<Vec<BroadphaseWorkItem>>,
}
impl WorkStealingBroadphase {
    /// Create a new work-stealing broadphase for the given set of AABBs.
    pub fn new(aabbs: Vec<ParAabb>, chunk_size: usize) -> Self {
        let pairs = aabb_all_pairs_overlap(&aabbs);
        let pending = pairs
            .into_iter()
            .map(|(a, b)| BroadphaseWorkItem {
                body_a: a,
                body_b: b,
            })
            .collect();
        Self {
            aabbs: Arc::new(aabbs),
            chunk_size,
            pending: Mutex::new(pending),
        }
    }
    /// Steal up to `chunk_size` work items.
    pub fn steal(&self) -> Vec<BroadphaseWorkItem> {
        let mut q = self.pending.lock().unwrap_or_else(|e| e.into_inner());
        let n = self.chunk_size.min(q.len());
        q.drain(..n).collect()
    }
    /// Returns `true` when no more work items remain.
    pub fn is_done(&self) -> bool {
        self.pending
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .is_empty()
    }
    /// Total number of candidate pairs remaining.
    pub fn remaining(&self) -> usize {
        self.pending.lock().unwrap_or_else(|e| e.into_inner()).len()
    }
    /// Access the underlying AABB slice.
    pub fn aabbs(&self) -> &[ParAabb] {
        &self.aabbs
    }
}
/// Thread-safe AABB tree for parallel broadphase collision detection.
///
/// The tree is built once per frame and queried from multiple threads
/// simultaneously via `Arc<RwLock<>>`. Insertions/removals lock the writer,
/// queries only need a read lock.
#[derive(Debug)]
pub struct ParallelAabbTree {
    pub(super) inner: RwLock<AabbTreeInner>,
}
impl ParallelAabbTree {
    /// Create an empty parallel AABB tree.
    pub fn new() -> Self {
        Self {
            inner: RwLock::new(AabbTreeInner::new()),
        }
    }
    /// Insert a body into the tree.
    pub fn insert(&self, aabb: ParAabb, handle: u32) -> NodeIdx {
        self.inner
            .write()
            .unwrap_or_else(|e| e.into_inner())
            .insert(aabb, handle)
    }
    /// Query overlapping bodies for a given AABB.
    ///
    /// Multiple threads can call this simultaneously.
    pub fn query(&self, aabb: &ParAabb) -> Vec<u32> {
        self.inner
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .query_overlaps(aabb)
    }
    /// Number of nodes currently allocated.
    pub fn node_count(&self) -> usize {
        self.inner
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .nodes
            .len()
    }
}
/// A batch of shape pairs to be processed in parallel.
#[derive(Debug)]
pub struct NarrowphaseBatch {
    /// Pairs to process: each entry is (index_a, kind_a, index_b, kind_b).
    pub pairs: Vec<(u32, ParShapeKind, u32, ParShapeKind)>,
}
impl NarrowphaseBatch {
    /// Create an empty batch.
    pub fn new() -> Self {
        Self { pairs: Vec::new() }
    }
    /// Add a pair to the batch.
    pub fn push(&mut self, a_idx: u32, a: ParShapeKind, b_idx: u32, b: ParShapeKind) {
        self.pairs.push((a_idx, a, b_idx, b));
    }
    /// Process all pairs sequentially, returning all contacts.
    ///
    /// In a real engine this loop body would be dispatched across threads.
    pub fn process(&self) -> Vec<ParContact> {
        let mut contacts = Vec::new();
        for &(a_idx, ref a_shape, b_idx, ref b_shape) in &self.pairs {
            match (a_shape, b_shape) {
                (ParShapeKind::Sphere(sa), ParShapeKind::Sphere(sb)) => {
                    if let Some(c) = par_sphere_sphere(a_idx, b_idx, sa, sb) {
                        contacts.push(c);
                    }
                }
                (ParShapeKind::Sphere(sa), ParShapeKind::Box(bb)) => {
                    if let Some(c) = par_sphere_box(a_idx, b_idx, sa, bb) {
                        contacts.push(c);
                    }
                }
                (ParShapeKind::Box(ba), ParShapeKind::Sphere(sb)) => {
                    if let Some(mut c) = par_sphere_box(b_idx, a_idx, sb, ba) {
                        c.normal = scale3(c.normal, -1.0);
                        c.body_a = a_idx;
                        c.body_b = b_idx;
                        contacts.push(c);
                    }
                }
                (ParShapeKind::Box(_ba), ParShapeKind::Box(_bb)) => {
                    if let (ParShapeKind::Box(ba_), ParShapeKind::Box(bb_)) = (a_shape, b_shape)
                        && let Some(c) = par_box_box(a_idx, b_idx, ba_, bb_)
                    {
                        contacts.push(c);
                    }
                }
            }
        }
        contacts
    }
}
/// A GPU-mapped BVH node for use in compute shaders.
///
/// Packed as a 32-byte struct (two `\[f32; 3\]` + two `u32`).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct GpuBvhNode {
    /// AABB min corner (f32 for bandwidth efficiency).
    pub min: [f32; 3],
    /// Left child index, or `u32::MAX` for leaf.
    pub left_or_leaf: u32,
    /// AABB max corner.
    pub max: [f32; 3],
    /// Right child index, or body handle for leaf.
    pub right_or_handle: u32,
}
impl GpuBvhNode {
    /// Create a leaf node.
    pub fn leaf(min: [f32; 3], max: [f32; 3], handle: u32) -> Self {
        Self {
            min,
            max,
            left_or_leaf: u32::MAX,
            right_or_handle: handle,
        }
    }
    /// Create an internal node.
    pub fn internal(min: [f32; 3], max: [f32; 3], left: u32, right: u32) -> Self {
        Self {
            min,
            max,
            left_or_leaf: left,
            right_or_handle: right,
        }
    }
    /// Returns `true` if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left_or_leaf == u32::MAX
    }
}
/// Lock-free contact buffer using atomic index for concurrent writes.
///
/// Multiple threads can write contacts via `push_contact` using an atomic
/// slot counter. Each slot is protected by a `Mutex` for safe parallel access.
/// The buffer must be pre-allocated to capacity before parallel writes begin.
#[derive(Debug)]
pub struct LockFreeContactBuffer {
    pub(super) contacts: Vec<Mutex<Option<ParContact>>>,
    pub(super) write_idx: AtomicUsize,
    pub(super) capacity: usize,
}
impl LockFreeContactBuffer {
    /// Create a new lock-free contact buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        let mut contacts = Vec::with_capacity(capacity);
        for _ in 0..capacity {
            contacts.push(Mutex::new(None));
        }
        Self {
            contacts,
            write_idx: AtomicUsize::new(0),
            capacity,
        }
    }
    /// Push a contact, returning the slot index or `None` if full.
    pub fn push_contact(&self, contact: ParContact) -> Option<usize> {
        let idx = self.write_idx.fetch_add(1, Ordering::AcqRel);
        if idx < self.capacity {
            *self.contacts[idx].lock().unwrap_or_else(|e| e.into_inner()) = Some(contact);
            Some(idx)
        } else {
            self.write_idx.fetch_sub(1, Ordering::AcqRel);
            None
        }
    }
    /// Number of contacts written so far.
    pub fn len(&self) -> usize {
        self.write_idx.load(Ordering::Acquire).min(self.capacity)
    }
    /// Returns `true` if no contacts have been written.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Reset the buffer for reuse each frame.
    pub fn reset(&self) {
        let old = self.write_idx.swap(0, Ordering::AcqRel);
        for i in 0..old.min(self.capacity) {
            *self.contacts[i].lock().unwrap_or_else(|e| e.into_inner()) = None;
        }
    }
    /// Collect all written contacts into a `Vec`.
    pub fn collect(&self) -> Vec<ParContact> {
        let n = self.len();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            if let Some(c) = *self.contacts[i].lock().unwrap_or_else(|e| e.into_inner()) {
                out.push(c);
            }
        }
        out
    }
}
/// Configuration for the parallel collision pipeline.
#[derive(Debug, Clone)]
pub struct ParCollisionConfig {
    /// Maximum number of contacts per body pair after reduction.
    pub max_contacts_per_pair: usize,
    /// AABB expansion margin for predictive contacts (world units).
    pub aabb_margin: f64,
    /// Whether to run CCD sweep this frame.
    pub enable_ccd: bool,
}
/// A work item for the work-stealing broadphase scheduler.
#[derive(Debug, Clone)]
pub struct BroadphaseWorkItem {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
}
