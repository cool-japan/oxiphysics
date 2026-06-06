use super::functions::*;
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Broadphase collision detection kernel (CPU mock of a GPU dispatch).
///
/// Detects all overlapping AABB pairs from a flat list of `AabbGpu`.
#[derive(Debug, Default)]
pub struct BroadphaseGpuKernel {
    /// Inflation margin added to each AABB before testing.
    pub margin: f32,
}
impl BroadphaseGpuKernel {
    /// Create a new `BroadphaseGpuKernel` with the given inflation margin.
    pub fn new(margin: f32) -> Self {
        Self { margin }
    }
    /// Dispatch the broadphase over `aabbs`.
    ///
    /// Returns all overlapping pairs (canonical form, `a < b`).
    pub fn dispatch(&self, aabbs: &[AabbGpu]) -> Vec<CollisionPair> {
        let mut pairs = Vec::new();
        for (i, aabb_i) in aabbs.iter().enumerate() {
            let expanded_i = aabb_i.expanded(self.margin);
            for (j, aabb_j) in aabbs.iter().enumerate().skip(i + 1) {
                if expanded_i.overlaps(aabb_j) {
                    pairs.push(CollisionPair::new(i as u32, j as u32));
                }
            }
        }
        pairs
    }
    /// BVH-accelerated broadphase dispatch.
    pub fn dispatch_bvh(&self, aabbs: &[AabbGpu]) -> Vec<CollisionPair> {
        let mut builder = GpuBvhBuilder::new();
        let root = builder.build(aabbs);
        let mut pairs = Vec::new();
        for (i, aabb) in aabbs.iter().enumerate() {
            let query = aabb.expanded(self.margin);
            let hits = builder.query_overlaps(root, &query);
            for j in hits {
                if j as usize > i {
                    pairs.push(CollisionPair::new(i as u32, j));
                }
            }
        }
        pairs
    }
}
/// Axis-aligned bounding box for GPU collision passes.
///
/// Corners are stored as `[f32; 3]` to match typical GPU buffer layouts.
#[derive(Debug, Clone, PartialEq)]
pub struct AabbGpu {
    /// Component-wise minimum corner.
    pub min: [f32; 3],
    /// Component-wise maximum corner.
    pub max: [f32; 3],
}
impl AabbGpu {
    /// Construct an `AabbGpu` from explicit min/max corners.
    pub fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self { min, max }
    }
    /// Construct a degenerate `AabbGpu` that contains exactly one point.
    pub fn from_point(p: [f32; 3]) -> Self {
        Self { min: p, max: p }
    }
    /// Return the smallest `AabbGpu` that contains both `self` and `other`.
    pub fn merge(&self, other: &AabbGpu) -> AabbGpu {
        AabbGpu {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }
    /// Returns `true` if this box overlaps `other` (touching counts as overlap).
    pub fn overlaps(&self, other: &AabbGpu) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Returns `true` if the point `p` lies inside or on the surface.
    pub fn contains_point(&self, p: [f32; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
    /// Surface area of the AABB (sum of the six face areas).
    pub fn surface_area(&self) -> f32 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        2.0 * (dx * dy + dy * dz + dz * dx)
    }
    /// Volume of the AABB.
    pub fn volume(&self) -> f32 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        dx * dy * dz
    }
    /// Centre of the AABB.
    pub fn centre(&self) -> [f32; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> [f32; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }
    /// Expand the AABB by `margin` in every direction.
    pub fn expanded(&self, margin: f32) -> AabbGpu {
        AabbGpu {
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
}
/// Persistent 4-point contact manifold with warm-start impulse data.
///
/// Retains up to 4 contact points between a pair of bodies across frames,
/// selecting the deepest/most-separated subset for stability.
#[derive(Debug, Clone)]
pub struct PersistentManifoldGpu {
    /// Index of body A.
    pub body_a: u32,
    /// Index of body B.
    pub body_b: u32,
    /// Retained contact points (up to 4).
    pub points: [Option<ManifoldPoint>; 4],
    /// Number of active contact points.
    pub num_points: usize,
}
impl PersistentManifoldGpu {
    /// Create an empty manifold for the given body pair.
    pub fn new(body_a: u32, body_b: u32) -> Self {
        Self {
            body_a,
            body_b,
            points: [None, None, None, None],
            num_points: 0,
        }
    }
    /// Add or replace a contact point.
    ///
    /// If there are fewer than 4 points the new point is simply appended.
    /// Otherwise the shallowest existing point is replaced if it is shallower
    /// than the new one.
    pub fn add_contact(&mut self, new_pt: ManifoldPoint) {
        if self.num_points < 4 {
            self.points[self.num_points] = Some(new_pt);
            self.num_points += 1;
            return;
        }
        let mut shallowest_idx = 0;
        let mut shallowest_depth = f32::INFINITY;
        for (i, opt) in self.points.iter().enumerate() {
            if let Some(p) = opt
                && p.depth < shallowest_depth
            {
                shallowest_depth = p.depth;
                shallowest_idx = i;
            }
        }
        if new_pt.depth > shallowest_depth {
            self.points[shallowest_idx] = Some(new_pt);
        }
    }
    /// Remove stale contact points whose depth has become positive (separated).
    pub fn remove_stale(&mut self, threshold: f32) {
        for slot in self.points.iter_mut() {
            if let Some(p) = slot
                && p.depth < threshold
            {
                *slot = None;
            }
        }
        self.num_points = self.points.iter().filter(|s| s.is_some()).count();
    }
    /// Reset warm-start impulses to zero.
    pub fn reset_warm_start(&mut self) {
        for slot in self.points.iter_mut().flatten() {
            slot.warm_impulse_normal = 0.0;
            slot.warm_impulse_t1 = 0.0;
            slot.warm_impulse_t2 = 0.0;
        }
    }
    /// Return an iterator over active contact points.
    pub fn active_points(&self) -> impl Iterator<Item = &ManifoldPoint> {
        self.points.iter().flatten()
    }
}
/// Builds a flat BVH from a list of AABBs using Morton-code sorting.
///
/// The builder sorts primitives by their Morton code (computed from the centroid
/// of each AABB), then constructs the BVH bottom-up.
#[derive(Debug, Default)]
pub struct GpuBvhBuilder {
    pub(super) nodes: Vec<BvhNodeGpu>,
}
impl GpuBvhBuilder {
    /// Create a new, empty `GpuBvhBuilder`.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }
    /// Build a BVH from a slice of `AabbGpu` primitives.
    ///
    /// Returns the index of the root node in the returned flat node array.
    pub fn build(&mut self, aabbs: &[AabbGpu]) -> usize {
        self.nodes.clear();
        if aabbs.is_empty() {
            return 0;
        }
        let scene_min = [
            aabbs.iter().map(|a| a.min[0]).fold(f32::INFINITY, f32::min),
            aabbs.iter().map(|a| a.min[1]).fold(f32::INFINITY, f32::min),
            aabbs.iter().map(|a| a.min[2]).fold(f32::INFINITY, f32::min),
        ];
        let scene_max = [
            aabbs
                .iter()
                .map(|a| a.max[0])
                .fold(f32::NEG_INFINITY, f32::max),
            aabbs
                .iter()
                .map(|a| a.max[1])
                .fold(f32::NEG_INFINITY, f32::max),
            aabbs
                .iter()
                .map(|a| a.max[2])
                .fold(f32::NEG_INFINITY, f32::max),
        ];
        let extent = [
            (scene_max[0] - scene_min[0]).max(1e-6),
            (scene_max[1] - scene_min[1]).max(1e-6),
            (scene_max[2] - scene_min[2]).max(1e-6),
        ];
        let mut indexed: Vec<(u32, usize)> = aabbs
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let c = a.centre();
                let nx = (c[0] - scene_min[0]) / extent[0];
                let ny = (c[1] - scene_min[1]) / extent[1];
                let nz = (c[2] - scene_min[2]) / extent[2];
                (morton_code(nx, ny, nz), i)
            })
            .collect();
        indexed.sort_by_key(|&(code, _)| code);
        let leaf_base = 0usize;
        for &(_, prim_idx) in &indexed {
            self.nodes
                .push(BvhNodeGpu::leaf(aabbs[prim_idx].clone(), prim_idx as u32));
        }
        let mut current_level: Vec<usize> = (leaf_base..leaf_base + indexed.len()).collect();
        while current_level.len() > 1 {
            let mut next_level = Vec::new();
            let mut i = 0;
            while i < current_level.len() {
                if i + 1 < current_level.len() {
                    let l = current_level[i];
                    let r = current_level[i + 1];
                    let merged = self.nodes[l].aabb.merge(&self.nodes[r].aabb);
                    let internal = BvhNodeGpu::internal(merged, l as u32, r as u32);
                    let idx = self.nodes.len();
                    self.nodes.push(internal);
                    next_level.push(idx);
                    i += 2;
                } else {
                    next_level.push(current_level[i]);
                    i += 1;
                }
            }
            current_level = next_level;
        }
        current_level[0]
    }
    /// Return a reference to the built flat node array.
    pub fn nodes(&self) -> &[BvhNodeGpu] {
        &self.nodes
    }
    /// Query all leaf primitives whose AABB overlaps `query`.
    pub fn query_overlaps(&self, root: usize, query: &AabbGpu) -> Vec<u32> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut result = Vec::new();
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.aabb.overlaps(query) {
                continue;
            }
            if node.is_leaf {
                result.push(node.primitive_index);
            } else {
                stack.push(node.left_child as usize);
                stack.push(node.right_child as usize);
            }
        }
        result
    }
}
/// A single contact point in the manifold.
#[derive(Debug, Clone)]
pub struct ManifoldPoint {
    /// Contact point in world space.
    pub position: [f32; 3],
    /// Contact normal (pointing from B to A).
    pub normal: [f32; 3],
    /// Penetration depth.
    pub depth: f32,
    /// Accumulated impulse for warm starting (normal direction).
    pub warm_impulse_normal: f32,
    /// Accumulated impulse for warm starting (tangent 1).
    pub warm_impulse_t1: f32,
    /// Accumulated impulse for warm starting (tangent 2).
    pub warm_impulse_t2: f32,
}
impl ManifoldPoint {
    /// Create a new manifold point from a contact result.
    pub fn from_contact(c: &ContactResult) -> Self {
        Self {
            position: c.contact_point,
            normal: c.normal,
            depth: c.depth,
            warm_impulse_normal: 0.0,
            warm_impulse_t1: 0.0,
            warm_impulse_t2: 0.0,
        }
    }
}
/// Narrowphase collision detection kernel (CPU mock with simplified GJK/SAT).
///
/// For AABB vs AABB pairs the kernel uses SAT (Separating Axis Theorem) which
/// is exact.  For more general shapes a GJK iteration would be used.
#[derive(Debug, Default)]
pub struct NarrowphaseGpuKernel;
impl NarrowphaseGpuKernel {
    /// Create a new `NarrowphaseGpuKernel`.
    pub fn new() -> Self {
        Self
    }
    /// Test a single AABB pair using SAT; returns `Some(ContactResult)` on
    /// contact, `None` otherwise.
    pub fn test_aabb_pair(
        &self,
        prim_a: u32,
        aabb_a: &AabbGpu,
        prim_b: u32,
        aabb_b: &AabbGpu,
    ) -> Option<ContactResult> {
        let axes: [[f32; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut min_depth = f32::INFINITY;
        let mut min_axis = [1.0f32, 0.0, 0.0];
        for axis in &axes {
            let a_min = dot3f(aabb_a.min, *axis);
            let a_max = dot3f(aabb_a.max, *axis);
            let b_min = dot3f(aabb_b.min, *axis);
            let b_max = dot3f(aabb_b.max, *axis);
            let overlap = a_max.min(b_max) - a_min.max(b_min);
            if overlap <= 0.0 {
                return None;
            }
            if overlap < min_depth {
                min_depth = overlap;
                let ca = dot3f(aabb_a.centre(), *axis);
                let cb = dot3f(aabb_b.centre(), *axis);
                let sign = if ca > cb { 1.0 } else { -1.0 };
                min_axis = scale3f(*axis, sign);
            }
        }
        let ca = aabb_a.centre();
        let cb = aabb_b.centre();
        let contact_point = [
            (ca[0] + cb[0]) * 0.5,
            (ca[1] + cb[1]) * 0.5,
            (ca[2] + cb[2]) * 0.5,
        ];
        Some(ContactResult {
            prim_a,
            prim_b,
            contact_point,
            normal: min_axis,
            depth: min_depth,
        })
    }
    /// Run GJK for two AABB shapes (CPU mock).
    ///
    /// Returns `true` if the shapes intersect, along with the final simplex
    /// distance (0.0 on intersection).
    pub fn gjk_intersect(&self, a: &AabbGpu, b: &AabbGpu) -> (bool, f32) {
        let ca = a.centre();
        let cb = b.centre();
        let mut d = sub3f(ca, cb);
        if len_sq3f(d) < 1e-9 {
            d = [1.0, 0.0, 0.0];
        }
        let s0 = minkowski_support(a, b, d);
        let mut simplex = vec![s0];
        d = [-s0.v[0], -s0.v[1], -s0.v[2]];
        for _ in 0..32 {
            if len_sq3f(d) < 1e-12 {
                return (true, 0.0);
            }
            let a_sup = minkowski_support(a, b, d);
            if dot3f(a_sup.v, d) < 0.0 {
                return (false, len3f(d));
            }
            simplex.push(a_sup);
            let (intersecting, new_d) = update_simplex(&mut simplex);
            if intersecting {
                return (true, 0.0);
            }
            d = new_d;
        }
        (true, 0.0)
    }
    /// Dispatch narrowphase over a list of broadphase pairs and AABBs.
    pub fn dispatch(&self, pairs: &[CollisionPair], aabbs: &[AabbGpu]) -> Vec<ContactResult> {
        let mut contacts = Vec::new();
        for pair in pairs {
            let a = pair.a as usize;
            let b = pair.b as usize;
            if a < aabbs.len()
                && b < aabbs.len()
                && let Some(c) = self.test_aabb_pair(pair.a, &aabbs[a], pair.b, &aabbs[b])
            {
                contacts.push(c);
            }
        }
        contacts
    }
}
/// A single BVH node in GPU-friendly flat representation.
///
/// Leaf nodes have `is_leaf == true` and `primitive_index` set.
/// Internal nodes use `left_child` / `right_child` as indices into the node
/// array.
#[derive(Debug, Clone)]
pub struct BvhNodeGpu {
    /// Bounding box for this node.
    pub aabb: AabbGpu,
    /// Index of the left child node (internal node only).
    pub left_child: u32,
    /// Index of the right child node (internal node only).
    pub right_child: u32,
    /// Whether this node is a leaf.
    pub is_leaf: bool,
    /// Index of the primitive stored in this leaf.
    pub primitive_index: u32,
}
impl BvhNodeGpu {
    /// Create a new internal BVH node.
    pub fn internal(aabb: AabbGpu, left_child: u32, right_child: u32) -> Self {
        Self {
            aabb,
            left_child,
            right_child,
            is_leaf: false,
            primitive_index: u32::MAX,
        }
    }
    /// Create a new leaf BVH node.
    pub fn leaf(aabb: AabbGpu, primitive_index: u32) -> Self {
        Self {
            aabb,
            left_child: u32::MAX,
            right_child: u32::MAX,
            is_leaf: true,
            primitive_index,
        }
    }
}
/// Combined broadphase + narrowphase collision pipeline.
///
/// Holds the configuration for both stages and maintains a list of persistent
/// manifolds from the previous frame for warm starting.
#[derive(Debug)]
pub struct CollisionGpuPipeline {
    /// Broadphase kernel configuration.
    pub broadphase: BroadphaseGpuKernel,
    /// Narrowphase kernel.
    pub narrowphase: NarrowphaseGpuKernel,
    /// Persistent manifolds from the previous timestep.
    pub manifolds: Vec<PersistentManifoldGpu>,
    /// Contact list from the most recent dispatch.
    pub contacts: Vec<ContactResult>,
}
impl CollisionGpuPipeline {
    /// Create a new pipeline with default broadphase margin.
    pub fn new(margin: f32) -> Self {
        Self {
            broadphase: BroadphaseGpuKernel::new(margin),
            narrowphase: NarrowphaseGpuKernel::new(),
            manifolds: Vec::new(),
            contacts: Vec::new(),
        }
    }
    /// Run one full collision detection pass over `aabbs`.
    ///
    /// 1. Broadphase: find overlapping AABB pairs.
    /// 2. Narrowphase: compute contact points for each pair.
    /// 3. Update persistent manifolds.
    pub fn dispatch(&mut self, aabbs: &[AabbGpu]) {
        let pairs = self.broadphase.dispatch(aabbs);
        self.contacts = self.narrowphase.dispatch(&pairs, aabbs);
        self.update_manifolds();
    }
    /// Run the full pipeline using BVH-accelerated broadphase.
    pub fn dispatch_bvh(&mut self, aabbs: &[AabbGpu]) {
        let pairs = self.broadphase.dispatch_bvh(aabbs);
        self.contacts = self.narrowphase.dispatch(&pairs, aabbs);
        self.update_manifolds();
    }
    /// Update the persistent manifold list from the current contact set.
    fn update_manifolds(&mut self) {
        self.manifolds.retain(|m| m.num_points > 0);
        for contact in &self.contacts {
            let existing = self.manifolds.iter_mut().find(|m| {
                (m.body_a == contact.prim_a && m.body_b == contact.prim_b)
                    || (m.body_a == contact.prim_b && m.body_b == contact.prim_a)
            });
            let pt = ManifoldPoint::from_contact(contact);
            if let Some(manifold) = existing {
                manifold.add_contact(pt);
            } else {
                let mut new_m = PersistentManifoldGpu::new(contact.prim_a, contact.prim_b);
                new_m.add_contact(pt);
                self.manifolds.push(new_m);
            }
        }
    }
    /// Clear all manifolds and contact data.
    pub fn reset(&mut self) {
        self.manifolds.clear();
        self.contacts.clear();
    }
    /// Return the total number of contact points across all manifolds.
    pub fn total_contact_points(&self) -> usize {
        self.manifolds.iter().map(|m| m.num_points).sum()
    }
}
/// Result of a GJK distance query.
#[derive(Debug, Clone)]
pub struct GjkResult {
    /// Whether the shapes are intersecting (distance = 0).
    pub intersecting: bool,
    /// Minimum distance between the shapes (0 if intersecting).
    pub distance: f32,
    /// Closest point on shape A.
    pub closest_a: [f32; 3],
    /// Closest point on shape B.
    pub closest_b: [f32; 3],
    /// Separation axis (unit vector from B to A).
    pub axis: [f32; 3],
}
/// Full GPU collision pipeline: broadphase (SAP) + narrowphase (GJK) +
/// persistent contact cache with warmstarting.
///
/// Intended as the primary collision API. Configure once, call `step` each
/// physics tick to get a fresh contact list.
#[derive(Debug)]
pub struct GpuCollisionPipeline {
    /// SAP broadphase.
    pub broadphase: GpuBroadphase,
    /// GJK narrowphase.
    pub narrowphase: GpuNarrowphase,
    /// Persistent contact cache.
    pub contact_cache: GpuContactCache,
    /// AABB tree (rebuilt each frame or on-demand).
    pub aabb_tree: GpuAabbTree,
    /// Contact list from the most recent step.
    pub contacts: Vec<ContactResult>,
    /// Persistent manifolds.
    pub manifolds: Vec<PersistentManifoldGpu>,
    /// Cumulative stats since pipeline creation.
    pub cumulative_stats: CollisionKernelStats,
}
impl GpuCollisionPipeline {
    /// Create a new `GpuCollisionPipeline` with sensible defaults.
    pub fn new(margin: f32, max_cache_age: u32) -> Self {
        Self {
            broadphase: GpuBroadphase::new(0, margin),
            narrowphase: GpuNarrowphase::new(32),
            contact_cache: GpuContactCache::new(max_cache_age),
            aabb_tree: GpuAabbTree::new(),
            contacts: Vec::new(),
            manifolds: Vec::new(),
            cumulative_stats: CollisionKernelStats::new(),
        }
    }
    /// Run one full collision step using SAP broadphase + GJK narrowphase.
    ///
    /// Steps:
    /// 1. Rebuild the AABB tree.
    /// 2. Run SAP broadphase.
    /// 3. Run GJK narrowphase on each pair.
    /// 4. Refresh contact cache; apply warmstart data.
    /// 5. Update persistent manifolds.
    /// 6. Age out stale cache entries.
    pub fn step(&mut self, aabbs: &[AabbGpu]) {
        self.aabb_tree.build(aabbs);
        self.broadphase.sweep_axis = GpuBroadphase::choose_best_axis(aabbs);
        let pairs = self.broadphase.dispatch(aabbs);
        self.contacts = self.narrowphase.dispatch(&pairs, aabbs);
        for (pair, contact) in pairs.iter().zip(self.contacts.iter()) {
            self.contact_cache.insert(pair, contact);
        }
        self.update_manifolds();
        self.contact_cache.advance_frame();
        let mut step_stats = self.broadphase.stats.clone();
        step_stats.accumulate(&self.narrowphase.stats);
        self.cumulative_stats.accumulate(&step_stats);
    }
    fn update_manifolds(&mut self) {
        self.manifolds.retain(|m| m.num_points > 0);
        for contact in &self.contacts {
            let existing = self.manifolds.iter_mut().find(|m| {
                (m.body_a == contact.prim_a && m.body_b == contact.prim_b)
                    || (m.body_a == contact.prim_b && m.body_b == contact.prim_a)
            });
            let mut pt = ManifoldPoint::from_contact(contact);
            if let Some(entry) = self.contact_cache.find(contact.prim_a, contact.prim_b) {
                entry.apply_warm_start(&mut pt);
            }
            if let Some(manifold) = existing {
                manifold.add_contact(pt);
            } else {
                let mut new_m = PersistentManifoldGpu::new(contact.prim_a, contact.prim_b);
                new_m.add_contact(pt);
                self.manifolds.push(new_m);
            }
        }
    }
    /// Clear all contacts, manifolds, and the contact cache.
    pub fn reset(&mut self) {
        self.contacts.clear();
        self.manifolds.clear();
        self.contact_cache.clear();
    }
    /// Total active contact points across all manifolds.
    pub fn total_contact_points(&self) -> usize {
        self.manifolds.iter().map(|m| m.num_points).sum()
    }
    /// Get a snapshot of per-frame broadphase + narrowphase stats.
    pub fn frame_stats(&self) -> CollisionKernelStats {
        let mut s = self.broadphase.stats.clone();
        s.accumulate(&self.narrowphase.stats);
        s
    }
}
/// Broadphase entry used in Sort-and-Sweep.
#[derive(Debug, Clone)]
pub(super) struct SapEntry {
    /// Minimum along the sweep axis.
    pub(super) lo: f32,
    /// Maximum along the sweep axis.
    pub(super) hi: f32,
    /// Primitive index.
    pub(super) idx: u32,
}
/// A pair of primitive indices that may be in contact.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CollisionPair {
    /// Index of the first primitive.
    pub a: u32,
    /// Index of the second primitive.
    pub b: u32,
}
impl CollisionPair {
    /// Construct a canonical pair with `a <= b`.
    pub fn new(a: u32, b: u32) -> Self {
        if a <= b {
            Self { a, b }
        } else {
            Self { a: b, b: a }
        }
    }
}
/// Contact information produced by the narrowphase kernel.
#[derive(Debug, Clone)]
pub struct ContactResult {
    /// Index of the first primitive.
    pub prim_a: u32,
    /// Index of the second primitive.
    pub prim_b: u32,
    /// Contact point in world space (midpoint approximation).
    pub contact_point: [f32; 3],
    /// Contact normal pointing from B towards A.
    pub normal: [f32; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f32,
}
/// Support point difference used internally by GJK.
#[derive(Debug, Clone, Copy)]
pub(super) struct SupportPoint {
    /// The point in the Minkowski difference.
    pub(super) v: [f32; 3],
}
impl SupportPoint {
    pub(super) fn new(v: [f32; 3]) -> Self {
        Self { v }
    }
}
/// Statistics for a single collision kernel dispatch.
#[derive(Debug, Clone, Default)]
pub struct CollisionKernelStats {
    /// Total bytes transferred (read + write) during the dispatch.
    pub bytes_transferred: u64,
    /// Number of primitive pairs tested in the broadphase.
    pub broadphase_pairs_tested: u64,
    /// Number of pairs that passed the broadphase (potential contacts).
    pub broadphase_hits: u64,
    /// Number of narrowphase GJK/SAT queries executed.
    pub narrowphase_queries: u64,
    /// Number of narrowphase queries that resulted in a contact.
    pub narrowphase_contacts: u64,
    /// Elapsed wall-clock time for this dispatch in seconds.
    pub elapsed_secs: f64,
}
impl CollisionKernelStats {
    /// Create a zero-initialised stats record.
    pub fn new() -> Self {
        Self::default()
    }
    /// Broadphase hit-rate in `[0, 1]`.
    ///
    /// Returns `0.0` when no pairs were tested.
    pub fn broadphase_hit_rate(&self) -> f64 {
        if self.broadphase_pairs_tested == 0 {
            return 0.0_f64;
        }
        self.broadphase_hits as f64 / self.broadphase_pairs_tested as f64
    }
    /// Narrowphase contact-rate in `[0, 1]`.
    ///
    /// Returns `0.0` when no narrowphase queries were executed.
    pub fn narrowphase_contact_rate(&self) -> f64 {
        if self.narrowphase_queries == 0 {
            return 0.0_f64;
        }
        self.narrowphase_contacts as f64 / self.narrowphase_queries as f64
    }
    /// Estimated memory bandwidth in GB/s.
    ///
    /// Returns `0.0` when elapsed time is zero.
    pub fn bandwidth_gb_s(&self) -> f64 {
        if self.elapsed_secs <= 0.0_f64 {
            return 0.0_f64;
        }
        self.bytes_transferred as f64 / self.elapsed_secs / 1.0e9_f64
    }
    /// Pair throughput in pairs-per-second.
    ///
    /// Returns `0.0` when elapsed time is zero.
    pub fn pair_throughput(&self) -> f64 {
        if self.elapsed_secs <= 0.0_f64 {
            return 0.0_f64;
        }
        self.broadphase_pairs_tested as f64 / self.elapsed_secs
    }
    /// Accumulate another stats record into this one.
    pub fn accumulate(&mut self, other: &CollisionKernelStats) {
        self.bytes_transferred += other.bytes_transferred;
        self.broadphase_pairs_tested += other.broadphase_pairs_tested;
        self.broadphase_hits += other.broadphase_hits;
        self.narrowphase_queries += other.narrowphase_queries;
        self.narrowphase_contacts += other.narrowphase_contacts;
        self.elapsed_secs += other.elapsed_secs;
    }
}
/// A single entry in the contact cache for warmstarting.
#[derive(Debug, Clone)]
pub struct ContactCacheEntry {
    /// Body pair key (a, b) with a < b.
    pub key: (u32, u32),
    /// Accumulated normal impulse from the previous frame.
    pub accumulated_normal: f32,
    /// Accumulated tangent impulse (direction 1) from the previous frame.
    pub accumulated_tangent_1: f32,
    /// Accumulated tangent impulse (direction 2) from the previous frame.
    pub accumulated_tangent_2: f32,
    /// Contact point used for positional warmstart.
    pub contact_point: [f32; 3],
    /// Age in frames since this entry was last refreshed.
    pub age_frames: u32,
}
impl ContactCacheEntry {
    /// Create a new cache entry for the given pair and contact.
    pub fn new(pair: &CollisionPair, contact: &ContactResult) -> Self {
        let (a, b) = if pair.a <= pair.b {
            (pair.a, pair.b)
        } else {
            (pair.b, pair.a)
        };
        Self {
            key: (a, b),
            accumulated_normal: 0.0_f32,
            accumulated_tangent_1: 0.0_f32,
            accumulated_tangent_2: 0.0_f32,
            contact_point: contact.contact_point,
            age_frames: 0,
        }
    }
    /// Apply warmstart data to a manifold point.
    pub fn apply_warm_start(&self, pt: &mut ManifoldPoint) {
        pt.warm_impulse_normal = self.accumulated_normal;
        pt.warm_impulse_t1 = self.accumulated_tangent_1;
        pt.warm_impulse_t2 = self.accumulated_tangent_2;
    }
    /// Update from a resolved contact (call after solver).
    pub fn update(&mut self, normal_impulse: f32, t1: f32, t2: f32) {
        self.accumulated_normal = normal_impulse;
        self.accumulated_tangent_1 = t1;
        self.accumulated_tangent_2 = t2;
        self.age_frames = 0;
    }
}
/// Flat AABB tree (BVH) built with parallel bottom-up Morton-code construction.
///
/// The tree stores nodes in a flat `Vec`BvhNodeGpu` for cache-friendly GPU
/// traversal.  Construction uses a parallel (CPU-mock) bottom-up merge pass
/// after sorting leaf centres by Morton code.
#[derive(Debug, Default)]
pub struct GpuAabbTree {
    /// Flat array of BVH nodes.
    pub nodes: Vec<BvhNodeGpu>,
    /// Index of the root node.
    pub root: usize,
    /// Number of leaf primitives.
    pub num_primitives: usize,
    /// Morton codes computed during last build.
    pub morton_codes: Vec<u32>,
}
impl GpuAabbTree {
    /// Create a new empty `GpuAabbTree`.
    pub fn new() -> Self {
        Self::default()
    }
    /// Build the tree from a slice of `AabbGpu` primitives.
    ///
    /// Sorts primitives by Morton code, inserts leaf nodes, then performs
    /// a CPU-mock parallel bottom-up merge to form internal nodes.
    pub fn build(&mut self, aabbs: &[AabbGpu]) {
        self.nodes.clear();
        self.morton_codes.clear();
        self.num_primitives = aabbs.len();
        if aabbs.is_empty() {
            self.root = 0;
            return;
        }
        let mut scene_min = [f32::INFINITY; 3];
        let mut scene_max = [f32::NEG_INFINITY; 3];
        for a in aabbs {
            for d in 0..3 {
                scene_min[d] = scene_min[d].min(a.min[d]);
                scene_max[d] = scene_max[d].max(a.max[d]);
            }
        }
        let extent = [
            (scene_max[0] - scene_min[0]).max(1.0e-6_f32),
            (scene_max[1] - scene_min[1]).max(1.0e-6_f32),
            (scene_max[2] - scene_min[2]).max(1.0e-6_f32),
        ];
        let mut order: Vec<(u32, usize)> = aabbs
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let c = a.centre();
                let nx = (c[0] - scene_min[0]) / extent[0];
                let ny = (c[1] - scene_min[1]) / extent[1];
                let nz = (c[2] - scene_min[2]) / extent[2];
                let code = morton_code(nx, ny, nz);
                self.morton_codes.push(code);
                (code, i)
            })
            .collect();
        order.sort_by_key(|&(code, _)| code);
        for &(_, prim_idx) in &order {
            self.nodes
                .push(BvhNodeGpu::leaf(aabbs[prim_idx].clone(), prim_idx as u32));
        }
        let mut level: Vec<usize> = (0..aabbs.len()).collect();
        while level.len() > 1 {
            let mut next = Vec::with_capacity(level.len().div_ceil(2));
            let mut i = 0;
            while i < level.len() {
                if i + 1 < level.len() {
                    let l = level[i];
                    let r = level[i + 1];
                    let merged = self.nodes[l].aabb.merge(&self.nodes[r].aabb);
                    let parent = BvhNodeGpu::internal(merged, l as u32, r as u32);
                    let idx = self.nodes.len();
                    self.nodes.push(parent);
                    next.push(idx);
                    i += 2;
                } else {
                    next.push(level[i]);
                    i += 1;
                }
            }
            level = next;
        }
        self.root = level[0];
    }
    /// Query all leaf primitives whose AABB overlaps `query`.
    pub fn query(&self, query: &AabbGpu) -> Vec<u32> {
        if self.nodes.is_empty() {
            return Vec::new();
        }
        let mut result = Vec::new();
        let mut stack = vec![self.root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.aabb.overlaps(query) {
                continue;
            }
            if node.is_leaf {
                result.push(node.primitive_index);
            } else {
                stack.push(node.left_child as usize);
                stack.push(node.right_child as usize);
            }
        }
        result
    }
    /// Tree depth (longest path from root to leaf).
    pub fn depth(&self) -> usize {
        if self.nodes.is_empty() {
            return 0;
        }
        self.depth_from(self.root)
    }
    fn depth_from(&self, idx: usize) -> usize {
        let node = &self.nodes[idx];
        if node.is_leaf {
            1
        } else {
            let ld = self.depth_from(node.left_child as usize);
            let rd = self.depth_from(node.right_child as usize);
            1 + ld.max(rd)
        }
    }
    /// Count the number of leaf nodes.
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf).count()
    }
    /// Surface Area Heuristic cost of the tree.
    ///
    /// SAH cost = sum over all internal nodes of `surface_area(node)`.
    pub fn sah_cost(&self) -> f32 {
        self.nodes
            .iter()
            .filter(|n| !n.is_leaf)
            .map(|n| n.aabb.surface_area())
            .sum()
    }
}
/// GPU narrowphase kernel performing parallel GJK distance queries.
///
/// Each pair from the broadphase is resolved with a dedicated GJK iteration
/// to compute exact contact information.
#[derive(Debug, Default)]
pub struct GpuNarrowphase {
    /// Maximum GJK iterations before declaring intersection.
    pub max_iterations: usize,
    /// Stats from the most recent dispatch.
    pub stats: CollisionKernelStats,
}
impl GpuNarrowphase {
    /// Create a new `GpuNarrowphase` with a given iteration limit.
    pub fn new(max_iterations: usize) -> Self {
        Self {
            max_iterations: max_iterations.max(4),
            stats: CollisionKernelStats::new(),
        }
    }
    /// Run a GJK distance query for two `AabbGpu` shapes.
    pub fn gjk_distance(&self, a: &AabbGpu, b: &AabbGpu) -> GjkResult {
        let ca = a.centre();
        let cb = b.centre();
        let mut d = sub3f(ca, cb);
        if len_sq3f(d) < 1.0e-9_f32 {
            d = [1.0_f32, 0.0_f32, 0.0_f32];
        }
        let s0 = minkowski_support(a, b, d);
        let mut simplex: Vec<SupportPoint> = vec![s0];
        d = [-s0.v[0], -s0.v[1], -s0.v[2]];
        for _ in 0..self.max_iterations {
            if len_sq3f(d) < 1.0e-12_f32 {
                return GjkResult {
                    intersecting: true,
                    distance: 0.0_f32,
                    closest_a: ca,
                    closest_b: cb,
                    axis: norm3f(sub3f(ca, cb)),
                };
            }
            let sup = minkowski_support(a, b, d);
            if dot3f(sup.v, d) < 0.0_f32 {
                let dist = len3f(d);
                let axis = norm3f(d);
                return GjkResult {
                    intersecting: false,
                    distance: dist,
                    closest_a: aabb_support(a, axis),
                    closest_b: aabb_support(b, [-axis[0], -axis[1], -axis[2]]),
                    axis,
                };
            }
            simplex.push(sup);
            let (inter, new_d) = update_simplex(&mut simplex);
            if inter {
                return GjkResult {
                    intersecting: true,
                    distance: 0.0_f32,
                    closest_a: ca,
                    closest_b: cb,
                    axis: norm3f(sub3f(ca, cb)),
                };
            }
            d = new_d;
        }
        GjkResult {
            intersecting: true,
            distance: 0.0_f32,
            closest_a: ca,
            closest_b: cb,
            axis: norm3f(sub3f(ca, cb)),
        }
    }
    /// Compute penetration depth from an overlapping AABB pair using SAT.
    ///
    /// Returns `None` if the shapes do not overlap.
    pub fn sat_depth(&self, a: &AabbGpu, b: &AabbGpu) -> Option<(f32, [f32; 3])> {
        let mut min_depth = f32::INFINITY;
        let mut min_axis = [1.0_f32, 0.0_f32, 0.0_f32];
        for d in 0..3usize {
            let mut axis = [0.0_f32; 3];
            axis[d] = 1.0_f32;
            let a_min = a.min[d];
            let a_max = a.max[d];
            let b_min = b.min[d];
            let b_max = b.max[d];
            let overlap = a_max.min(b_max) - a_min.max(b_min);
            if overlap <= 0.0_f32 {
                return None;
            }
            if overlap < min_depth {
                min_depth = overlap;
                let sign = if (a.min[d] + a.max[d]) > (b.min[d] + b.max[d]) {
                    1.0_f32
                } else {
                    -1.0_f32
                };
                axis[d] = sign;
                min_axis = axis;
            }
        }
        Some((min_depth, min_axis))
    }
    /// Dispatch the narrowphase over a list of broadphase pairs.
    ///
    /// Returns full `ContactResult` list with GJK/SAT contact data.
    pub fn dispatch(&mut self, pairs: &[CollisionPair], aabbs: &[AabbGpu]) -> Vec<ContactResult> {
        let mut contacts = Vec::new();
        self.stats.narrowphase_queries += pairs.len() as u64;
        for pair in pairs {
            let ai = pair.a as usize;
            let bi = pair.b as usize;
            if ai >= aabbs.len() || bi >= aabbs.len() {
                continue;
            }
            let gjk = self.gjk_distance(&aabbs[ai], &aabbs[bi]);
            if gjk.intersecting {
                let (depth, normal) = self
                    .sat_depth(&aabbs[ai], &aabbs[bi])
                    .unwrap_or((0.0_f32, gjk.axis));
                let ca = aabbs[ai].centre();
                let cb = aabbs[bi].centre();
                let contact_point = [
                    (ca[0] + cb[0]) * 0.5_f32,
                    (ca[1] + cb[1]) * 0.5_f32,
                    (ca[2] + cb[2]) * 0.5_f32,
                ];
                contacts.push(ContactResult {
                    prim_a: pair.a,
                    prim_b: pair.b,
                    contact_point,
                    normal,
                    depth,
                });
                self.stats.narrowphase_contacts += 1;
            }
        }
        contacts
    }
}
/// Persistent contact cache for GPU broadphase/narrowphase warmstarting.
///
/// Maps body-pair keys to cached impulse data from the previous frame.
/// Entries are aged out after `max_age_frames` frames without a refresh.
#[derive(Debug, Default)]
pub struct GpuContactCache {
    /// Stored entries.
    pub entries: Vec<ContactCacheEntry>,
    /// Maximum age before an entry is evicted.
    pub max_age_frames: u32,
}
impl GpuContactCache {
    /// Create a new contact cache with the given maximum age.
    pub fn new(max_age_frames: u32) -> Self {
        Self {
            entries: Vec::new(),
            max_age_frames: max_age_frames.max(1),
        }
    }
    /// Look up a cached entry for the given pair key.
    pub fn find(&self, a: u32, b: u32) -> Option<&ContactCacheEntry> {
        let key = if a <= b { (a, b) } else { (b, a) };
        self.entries.iter().find(|e| e.key == key)
    }
    /// Look up a cached entry mutably.
    pub fn find_mut(&mut self, a: u32, b: u32) -> Option<&mut ContactCacheEntry> {
        let key = if a <= b { (a, b) } else { (b, a) };
        self.entries.iter_mut().find(|e| e.key == key)
    }
    /// Insert or refresh a contact cache entry.
    pub fn insert(&mut self, pair: &CollisionPair, contact: &ContactResult) {
        let a = pair.a;
        let b = pair.b;
        let key = if a <= b { (a, b) } else { (b, a) };
        if let Some(e) = self.entries.iter_mut().find(|e| e.key == key) {
            e.contact_point = contact.contact_point;
            e.age_frames = 0;
        } else {
            self.entries.push(ContactCacheEntry::new(pair, contact));
        }
    }
    /// Advance the cache by one frame: age entries and evict stale ones.
    pub fn advance_frame(&mut self) {
        for e in self.entries.iter_mut() {
            e.age_frames += 1;
        }
        let max_age = self.max_age_frames;
        self.entries.retain(|e| e.age_frames <= max_age);
    }
    /// Number of entries currently in the cache.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
/// GPU-accelerated broadphase using Sort-and-Sweep (SAP).
///
/// Sorts AABB projections along a chosen axis (default X), sweeps for
/// overlapping intervals, then validates the remaining axes and records
/// broadphase pairs.
#[derive(Debug, Default)]
pub struct GpuBroadphase {
    /// Axis along which to sort/sweep: 0=X, 1=Y, 2=Z.
    pub sweep_axis: usize,
    /// Inflation margin applied to each AABB before testing.
    pub margin: f32,
    /// Stats from the most recent dispatch.
    pub stats: CollisionKernelStats,
}
impl GpuBroadphase {
    /// Create a new `GpuBroadphase` with the given sweep axis and margin.
    pub fn new(sweep_axis: usize, margin: f32) -> Self {
        Self {
            sweep_axis: sweep_axis % 3,
            margin,
            stats: CollisionKernelStats::new(),
        }
    }
    /// Choose the axis with the widest spread for better SAP performance.
    ///
    /// Returns the axis index (0, 1, or 2).
    pub fn choose_best_axis(aabbs: &[AabbGpu]) -> usize {
        if aabbs.is_empty() {
            return 0;
        }
        let mut spread = [0.0_f32; 3];
        for (d, s) in spread.iter_mut().enumerate() {
            let lo = aabbs.iter().map(|a| a.min[d]).fold(f32::INFINITY, f32::min);
            let hi = aabbs
                .iter()
                .map(|a| a.max[d])
                .fold(f32::NEG_INFINITY, f32::max);
            *s = hi - lo;
        }
        if spread[0] >= spread[1] && spread[0] >= spread[2] {
            0
        } else if spread[1] >= spread[2] {
            1
        } else {
            2
        }
    }
    /// Run the SAP broadphase and return overlapping pairs.
    ///
    /// Complexity: O(n log n + k) where k is the number of pairs found.
    pub fn dispatch(&mut self, aabbs: &[AabbGpu]) -> Vec<CollisionPair> {
        let axis = self.sweep_axis;
        let _n = aabbs.len();
        let mut pairs = Vec::new();
        let mut entries: Vec<SapEntry> = aabbs
            .iter()
            .enumerate()
            .map(|(i, a)| SapEntry {
                lo: a.min[axis] - self.margin,
                hi: a.max[axis] + self.margin,
                idx: i as u32,
            })
            .collect();
        entries.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap_or(std::cmp::Ordering::Equal));
        let mut tests = 0u64;
        let mut hits = 0u64;
        for (i, ei) in entries.iter().enumerate() {
            for ej in entries.iter().skip(i + 1) {
                if ej.lo > ei.hi {
                    break;
                }
                tests += 1;
                let pi = ei.idx as usize;
                let pj = ej.idx as usize;
                let expanded_i = aabbs[pi].expanded(self.margin);
                if expanded_i.overlaps(&aabbs[pj]) {
                    hits += 1;
                    pairs.push(CollisionPair::new(ei.idx, ej.idx));
                }
            }
        }
        self.stats.broadphase_pairs_tested += tests;
        self.stats.broadphase_hits += hits;
        self.stats.bytes_transferred += std::mem::size_of_val(aabbs) as u64;
        pairs
    }
}
