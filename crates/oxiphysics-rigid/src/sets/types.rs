//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::body::RigidBody;
use crate::collider::Collider;
use oxiphysics_core::math::Vec3;
use oxiphysics_core::{BodyHandle, ColliderHandle, MassProperties};

/// A 3-D axis-aligned bounding box used for BVH construction and queries.
#[derive(Debug, Clone, PartialEq)]
pub struct Aabb3 {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl Aabb3 {
    /// Create a new AABB from `min` and `max` corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Create from a single point (degenerate AABB).
    pub fn from_point(p: [f64; 3]) -> Self {
        Self { min: p, max: p }
    }
    /// Create the smallest AABB enclosing a set of points.
    pub fn from_points(pts: &[[f64; 3]]) -> Option<Self> {
        if pts.is_empty() {
            return None;
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for p in &pts[1..] {
            for i in 0..3 {
                if p[i] < mn[i] {
                    mn[i] = p[i];
                }
                if p[i] > mx[i] {
                    mx[i] = p[i];
                }
            }
        }
        Some(Self { min: mn, max: mx })
    }
    /// Test whether `point` is inside this AABB (inclusive).
    pub fn contains(&self, point: [f64; 3]) -> bool {
        for ((p, mn), mx) in point.iter().zip(self.min.iter()).zip(self.max.iter()) {
            if p < mn || p > mx {
                return false;
            }
        }
        true
    }
    /// Test whether this AABB overlaps with `other`.
    pub fn overlaps(&self, other: &Aabb3) -> bool {
        for i in 0..3 {
            if self.min[i] > other.max[i] || self.max[i] < other.min[i] {
                return false;
            }
        }
        true
    }
    /// Surface area of this AABB.  Used as cost metric in SAH BVH.
    pub fn surface_area(&self) -> f64 {
        let lx = self.max[0] - self.min[0];
        let ly = self.max[1] - self.min[1];
        let lz = self.max[2] - self.min[2];
        2.0 * (lx * ly + lx * lz + ly * lz)
    }
    /// Volume of this AABB.
    pub fn volume(&self) -> f64 {
        (self.max[0] - self.min[0]) * (self.max[1] - self.min[1]) * (self.max[2] - self.min[2])
    }
    /// Centre of this AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Merge this AABB with `other`, producing the tightest enclosing AABB.
    pub fn merge(&self, other: &Aabb3) -> Aabb3 {
        let mut mn = [0.0f64; 3];
        let mut mx = [0.0f64; 3];
        for i in 0..3 {
            mn[i] = self.min[i].min(other.min[i]);
            mx[i] = self.max[i].max(other.max[i]);
        }
        Aabb3 { min: mn, max: mx }
    }
    /// Expand the AABB by `margin` in all directions.
    pub fn expand(&self, margin: f64) -> Aabb3 {
        Aabb3 {
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
    /// Compute the squared distance from `point` to the closest point on this
    /// AABB.  Returns 0 when the point is inside.
    pub fn sq_dist_to_point(&self, point: [f64; 3]) -> f64 {
        let mut sq = 0.0f64;
        for ((p, mn), mx) in point.iter().zip(self.min.iter()).zip(self.max.iter()) {
            if p < mn {
                let d = mn - p;
                sq += d * d;
            } else if p > mx {
                let d = p - mx;
                sq += d * d;
            }
        }
        sq
    }
}
/// Arena storage for colliders with generational handles.
#[derive(Debug, Default)]
pub struct ColliderSet {
    pub(super) colliders: Vec<Option<Collider>>,
    pub(super) generations: Vec<u32>,
    pub(super) free_list: Vec<u32>,
}
impl ColliderSet {
    /// Create a new empty collider set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a collider and return its handle.
    pub fn insert(&mut self, collider: Collider) -> ColliderHandle {
        if let Some(index) = self.free_list.pop() {
            let idx = index as usize;
            self.colliders[idx] = Some(collider);
            ColliderHandle::new(index, self.generations[idx])
        } else {
            let index = self.colliders.len() as u32;
            self.colliders.push(Some(collider));
            self.generations.push(0);
            ColliderHandle::new(index, 0)
        }
    }
    /// Get a reference to a collider by handle.
    pub fn get(&self, handle: ColliderHandle) -> Option<&Collider> {
        let idx = handle.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == handle.generation {
            self.colliders[idx].as_ref()
        } else {
            None
        }
    }
    /// Get a mutable reference to a collider by handle.
    pub fn get_mut(&mut self, handle: ColliderHandle) -> Option<&mut Collider> {
        let idx = handle.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == handle.generation {
            self.colliders[idx].as_mut()
        } else {
            None
        }
    }
    /// Remove a collider by handle.
    pub fn remove(&mut self, handle: ColliderHandle) -> Option<Collider> {
        let idx = handle.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == handle.generation {
            self.generations[idx] += 1;
            self.free_list.push(handle.index);
            self.colliders[idx].take()
        } else {
            None
        }
    }
    /// Number of active colliders.
    pub fn len(&self) -> usize {
        self.colliders.iter().filter(|c| c.is_some()).count()
    }
    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Insert a collider and return its handle together with the computed mass properties.
    pub fn insert_with_props(&mut self, collider: Collider) -> (ColliderHandle, MassProperties) {
        let props = collider.mass_properties();
        let handle = self.insert(collider);
        (handle, props)
    }
    /// Iterate over all active colliders.
    pub fn iter(&self) -> impl Iterator<Item = (ColliderHandle, &Collider)> {
        self.colliders.iter().enumerate().filter_map(move |(i, c)| {
            c.as_ref()
                .map(|col| (ColliderHandle::new(i as u32, self.generations[i]), col))
        })
    }
}
/// A complete BVH over a set of axis-aligned bounding boxes.
#[derive(Debug, Clone, Default)]
pub struct BvhTree {
    pub(super) nodes: Vec<BvhNode>,
    pub(super) root: Option<usize>,
}
impl BvhTree {
    /// Build a BVH from a list of `(id, aabb)` pairs.
    pub fn build(primitives: &[(usize, Aabb3)]) -> Self {
        let mut tree = BvhTree::default();
        if primitives.is_empty() {
            return tree;
        }
        let indices: Vec<usize> = (0..primitives.len()).collect();
        let root = tree.build_recursive(primitives, &indices);
        tree.root = Some(root);
        tree
    }
    fn build_recursive(&mut self, prims: &[(usize, Aabb3)], indices: &[usize]) -> usize {
        if indices.len() == 1 {
            let (id, aabb) = &prims[indices[0]];
            let node_idx = self.nodes.len();
            self.nodes.push(BvhNode::Leaf {
                aabb: aabb.clone(),
                id: *id,
            });
            return node_idx;
        }
        let combined = indices
            .iter()
            .map(|&i| prims[i].1.clone())
            .reduce(|a, b| a.merge(&b))
            .expect("operation should succeed");
        let extents = [
            combined.max[0] - combined.min[0],
            combined.max[1] - combined.min[1],
            combined.max[2] - combined.min[2],
        ];
        let axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
            0
        } else if extents[1] >= extents[2] {
            1
        } else {
            2
        };
        let mut sorted: Vec<usize> = indices.to_vec();
        sorted.sort_by(|&a, &b| {
            let ca = prims[a].1.center()[axis];
            let cb = prims[b].1.center()[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted.len() / 2;
        let (left_ids, right_ids) = sorted.split_at(mid);
        let left = self.build_recursive(prims, left_ids);
        let right = self.build_recursive(prims, right_ids);
        let node_idx = self.nodes.len();
        self.nodes.push(BvhNode::Internal {
            aabb: combined,
            left,
            right,
        });
        node_idx
    }
    /// Query the BVH for all primitives whose AABB overlaps `query`.
    pub fn query(&self, query: &Aabb3) -> Vec<usize> {
        let mut results = Vec::new();
        if let Some(root) = self.root {
            self.query_recursive(root, query, &mut results);
        }
        results
    }
    fn query_recursive(&self, node_idx: usize, query: &Aabb3, results: &mut Vec<usize>) {
        match &self.nodes[node_idx] {
            BvhNode::Leaf { aabb, id } => {
                if aabb.overlaps(query) {
                    results.push(*id);
                }
            }
            BvhNode::Internal { aabb, left, right } => {
                if aabb.overlaps(query) {
                    self.query_recursive(*left, query, results);
                    self.query_recursive(*right, query, results);
                }
            }
        }
    }
    /// Number of nodes in the tree.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Whether the tree is empty (no primitives).
    pub fn is_empty(&self) -> bool {
        self.root.is_none()
    }
    /// Rebuild the BVH from the current body set, using each body's AABB
    /// (centred at position, expanded by `margin`).
    pub fn rebuild_from_body_set(body_set: &RigidBodySet, margin: f64) -> Self {
        let primitives: Vec<(usize, Aabb3)> = body_set
            .iter()
            .map(|(h, b)| {
                let p = b.transform.position;
                let aabb = Aabb3::new(
                    [p.x - margin, p.y - margin, p.z - margin],
                    [p.x + margin, p.y + margin, p.z + margin],
                );
                (h.index as usize, aabb)
            })
            .collect();
        BvhTree::build(&primitives)
    }
}
/// Cache of recently seen body pairs (e.g. for warm-starting constraints).
///
/// Each entry maps a sorted pair `(min_idx, max_idx)` to a user-defined
/// `f64` value (e.g. accumulated impulse).
#[derive(Debug, Default, Clone)]
pub struct BodyPairCache {
    pub(super) data: std::collections::HashMap<(u32, u32), f64>,
}
impl BodyPairCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        Self::default()
    }
    fn key(a: oxiphysics_core::BodyHandle, b: oxiphysics_core::BodyHandle) -> (u32, u32) {
        let (lo, hi) = if a.index <= b.index {
            (a.index, b.index)
        } else {
            (b.index, a.index)
        };
        (lo, hi)
    }
    /// Insert or update a cached value for the body pair `(a, b)`.
    pub fn set(&mut self, a: oxiphysics_core::BodyHandle, b: oxiphysics_core::BodyHandle, v: f64) {
        self.data.insert(Self::key(a, b), v);
    }
    /// Retrieve the cached value for `(a, b)`, if present.
    pub fn get(
        &self,
        a: oxiphysics_core::BodyHandle,
        b: oxiphysics_core::BodyHandle,
    ) -> Option<f64> {
        self.data.get(&Self::key(a, b)).copied()
    }
    /// Remove an entry.
    pub fn remove(&mut self, a: oxiphysics_core::BodyHandle, b: oxiphysics_core::BodyHandle) {
        self.data.remove(&Self::key(a, b));
    }
    /// Remove all entries that reference a body handle no longer in `set`.
    pub fn evict_stale(&mut self, set: &RigidBodySet) {
        let valid: std::collections::HashSet<u32> = set.iter().map(|(h, _)| h.index).collect();
        self.data
            .retain(|(a, b), _| valid.contains(a) && valid.contains(b));
    }
    /// Clear the entire cache.
    pub fn clear(&mut self) {
        self.data.clear();
    }
    /// Number of cached pairs.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
/// Adjacency-list graph of bodies connected by constraints.
///
/// Used for island detection: connected components of the graph correspond to
/// independent simulation islands that can be solved in isolation.
#[derive(Debug, Default)]
pub struct ConstraintGraph {
    /// Adjacency list: `edges[i]` holds all body handles connected to body `i`.
    pub(super) edges: std::collections::HashMap<u32, Vec<u32>>,
}
impl ConstraintGraph {
    /// Create an empty graph.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a constraint edge between two body handles (undirected).
    pub fn add_constraint(&mut self, a: BodyHandle, b: BodyHandle) {
        self.edges.entry(a.index).or_default().push(b.index);
        self.edges.entry(b.index).or_default().push(a.index);
    }
    /// Remove all edges.
    pub fn clear(&mut self) {
        self.edges.clear();
    }
    /// Find all connected components (islands) using BFS.
    ///
    /// Returns a list of islands, where each island is a list of `BodyHandle`
    /// values (with generation 0 -- generation information is not stored in the
    /// graph).
    ///
    /// Bodies that appear in no constraints are returned as singleton islands.
    /// `all_handles` should be the set of all active body handles.
    pub fn find_islands(&self, all_handles: &[BodyHandle]) -> Vec<Vec<BodyHandle>> {
        let mut visited: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut islands: Vec<Vec<BodyHandle>> = Vec::new();
        for &handle in all_handles {
            if visited.contains(&handle.index) {
                continue;
            }
            let mut island = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(handle.index);
            visited.insert(handle.index);
            while let Some(idx) = queue.pop_front() {
                island.push(BodyHandle::new(idx, 0));
                if let Some(neighbors) = self.edges.get(&idx) {
                    for &nb in neighbors {
                        if !visited.contains(&nb) {
                            visited.insert(nb);
                            queue.push_back(nb);
                        }
                    }
                }
            }
            islands.push(island);
        }
        islands
    }
    /// Number of constraint edges (counting both directions).
    pub fn edge_count(&self) -> usize {
        self.edges.values().map(|v| v.len()).sum::<usize>() / 2
    }
    /// Number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.edges.len()
    }
    /// Degree of a specific node.
    pub fn degree(&self, handle: BodyHandle) -> usize {
        self.edges.get(&handle.index).map_or(0, |v| v.len())
    }
    /// Check if two bodies are directly connected.
    pub fn are_connected(&self, a: BodyHandle, b: BodyHandle) -> bool {
        self.edges
            .get(&a.index)
            .is_some_and(|neighbors| neighbors.contains(&b.index))
    }
}
/// A node in a top-down axis-aligned bounding-volume hierarchy (BVH).
///
/// The tree is built by recursively splitting along the longest axis at the
/// median position (median-cut SAH approximation).  Leaf nodes store a single
/// primitive ID.
#[derive(Debug, Clone)]
pub enum BvhNode {
    /// An internal node with a bounding volume and two children.
    Internal {
        /// Bounding volume of this subtree.
        aabb: Aabb3,
        /// Left child index into the flat node list.
        left: usize,
        /// Right child index into the flat node list.
        right: usize,
    },
    /// A leaf node storing a single primitive.
    Leaf {
        /// Bounding volume of this leaf.
        aabb: Aabb3,
        /// Primitive identifier (e.g. body slot index).
        id: usize,
    },
}
/// Arena storage for rigid bodies with generational handles.
#[derive(Debug, Default)]
pub struct RigidBodySet {
    pub(super) bodies: Vec<Option<RigidBody>>,
    pub(super) generations: Vec<u32>,
    pub(super) free_list: Vec<u32>,
}
impl RigidBodySet {
    /// Create a new empty body set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a body and return its handle.
    pub fn insert(&mut self, body: RigidBody) -> BodyHandle {
        if let Some(index) = self.free_list.pop() {
            let idx = index as usize;
            self.bodies[idx] = Some(body);
            BodyHandle::new(index, self.generations[idx])
        } else {
            let index = self.bodies.len() as u32;
            self.bodies.push(Some(body));
            self.generations.push(0);
            BodyHandle::new(index, 0)
        }
    }
    /// Get a reference to a body by handle.
    pub fn get(&self, handle: BodyHandle) -> Option<&RigidBody> {
        let idx = handle.index as usize;
        if idx < self.bodies.len() && self.generations[idx] == handle.generation {
            self.bodies[idx].as_ref()
        } else {
            None
        }
    }
    /// Get a mutable reference to a body by handle.
    pub fn get_mut(&mut self, handle: BodyHandle) -> Option<&mut RigidBody> {
        let idx = handle.index as usize;
        if idx < self.bodies.len() && self.generations[idx] == handle.generation {
            self.bodies[idx].as_mut()
        } else {
            None
        }
    }
    /// Remove a body by handle.
    pub fn remove(&mut self, handle: BodyHandle) -> Option<RigidBody> {
        let idx = handle.index as usize;
        if idx < self.bodies.len() && self.generations[idx] == handle.generation {
            self.generations[idx] += 1;
            self.free_list.push(handle.index);
            self.bodies[idx].take()
        } else {
            None
        }
    }
    /// Number of active bodies.
    pub fn len(&self) -> usize {
        self.bodies.iter().filter(|b| b.is_some()).count()
    }
    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Total capacity (including free slots).
    pub fn capacity(&self) -> usize {
        self.bodies.len()
    }
    /// Number of free slots available for reuse.
    pub fn free_count(&self) -> usize {
        self.free_list.len()
    }
    /// Iterate over all active bodies.
    pub fn iter(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.bodies.iter().enumerate().filter_map(move |(i, b)| {
            b.as_ref()
                .map(|body| (BodyHandle::new(i as u32, self.generations[i]), body))
        })
    }
    /// Iterate mutably over all active bodies.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (BodyHandle, &mut RigidBody)> {
        let generations = &self.generations;
        self.bodies
            .iter_mut()
            .enumerate()
            .filter_map(move |(i, b)| {
                b.as_mut()
                    .map(|body| (BodyHandle::new(i as u32, generations[i]), body))
            })
    }
    /// Iterate over bodies that are currently active (not sleeping, not static).
    pub fn active_bodies(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.state == crate::body::BodyState::Active)
    }
    /// Iterate over bodies that are currently sleeping.
    pub fn sleeping_bodies(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.state == crate::body::BodyState::Sleeping)
    }
    /// Iterate over dynamic bodies only.
    pub fn dynamic_bodies(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.body_type == crate::body::BodyType::Dynamic)
    }
    /// Iterate over static bodies only.
    pub fn static_bodies(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.body_type == crate::body::BodyType::Static)
    }
    /// Iterate over kinematic bodies only.
    pub fn kinematic_bodies(&self) -> impl Iterator<Item = (BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.body_type == crate::body::BodyType::Kinematic)
    }
    /// Collect all handles of currently active bodies.
    pub fn active_handles(&self) -> Vec<BodyHandle> {
        self.active_bodies().map(|(h, _)| h).collect()
    }
    /// Collect all valid handles.
    pub fn all_handles(&self) -> Vec<BodyHandle> {
        self.iter().map(|(h, _)| h).collect()
    }
    /// Find bodies whose mass is within the given range `[min_mass, max_mass]`.
    pub fn bodies_by_mass_range(
        &self,
        min_mass: f64,
        max_mass: f64,
    ) -> Vec<(BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| b.mass >= min_mass && b.mass <= max_mass)
            .collect()
    }
    /// Find the body with the largest mass.
    pub fn heaviest_body(&self) -> Option<(BodyHandle, &RigidBody)> {
        self.iter().max_by(|(_, a), (_, b)| {
            a.mass
                .partial_cmp(&b.mass)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
    /// Find the body with the smallest (positive) mass.
    pub fn lightest_dynamic_body(&self) -> Option<(BodyHandle, &RigidBody)> {
        self.dynamic_bodies()
            .filter(|(_, b)| b.mass > 0.0)
            .min_by(|(_, a), (_, b)| {
                a.mass
                    .partial_cmp(&b.mass)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }
    /// Find the body with the highest speed (linear velocity magnitude).
    pub fn fastest_body(&self) -> Option<(BodyHandle, &RigidBody)> {
        self.iter().max_by(|(_, a), (_, b)| {
            a.velocity
                .norm()
                .partial_cmp(&b.velocity.norm())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
    /// Find all bodies within `radius` of `point`.
    pub fn bodies_within_radius(
        &self,
        point: oxiphysics_core::math::Vec3,
        radius: f64,
    ) -> Vec<(BodyHandle, &RigidBody)> {
        let r2 = radius * radius;
        self.iter()
            .filter(|(_, b)| (b.transform.position - point).norm_squared() <= r2)
            .collect()
    }
    /// Find all bodies inside an AABB defined by `min_corner` and `max_corner`.
    pub fn bodies_in_aabb(
        &self,
        min_corner: oxiphysics_core::math::Vec3,
        max_corner: oxiphysics_core::math::Vec3,
    ) -> Vec<(BodyHandle, &RigidBody)> {
        self.iter()
            .filter(|(_, b)| {
                let p = b.transform.position;
                p.x >= min_corner.x
                    && p.x <= max_corner.x
                    && p.y >= min_corner.y
                    && p.y <= max_corner.y
                    && p.z >= min_corner.z
                    && p.z <= max_corner.z
            })
            .collect()
    }
    /// Find the nearest body to a given point.
    pub fn nearest_body(
        &self,
        point: oxiphysics_core::math::Vec3,
    ) -> Option<(BodyHandle, &RigidBody, f64)> {
        self.iter()
            .map(|(h, b)| {
                let dist = (b.transform.position - point).norm();
                (h, b, dist)
            })
            .min_by(|(_, _, da), (_, _, db)| {
                da.partial_cmp(db).unwrap_or(std::cmp::Ordering::Equal)
            })
    }
    /// Apply a function to all bodies matching a predicate.
    pub fn for_each_mut(
        &mut self,
        predicate: impl Fn(&RigidBody) -> bool,
        action: impl Fn(&mut RigidBody),
    ) {
        for (_, body) in self.iter_mut() {
            if predicate(body) {
                action(body);
            }
        }
    }
    /// Remove all bodies matching a predicate. Returns removed bodies and handles.
    pub fn remove_where(
        &mut self,
        predicate: impl Fn(&RigidBody) -> bool,
    ) -> Vec<(BodyHandle, RigidBody)> {
        let to_remove: Vec<BodyHandle> = self
            .iter()
            .filter(|(_, b)| predicate(b))
            .map(|(h, _)| h)
            .collect();
        let mut removed = Vec::new();
        for h in to_remove {
            if let Some(body) = self.remove(h) {
                removed.push((h, body));
            }
        }
        removed
    }
    /// Set velocity of all dynamic bodies to zero.
    pub fn freeze_all(&mut self) {
        use oxiphysics_core::math::Vec3;
        for (_, body) in self.iter_mut() {
            if body.body_type == crate::body::BodyType::Dynamic {
                body.velocity = Vec3::zeros();
                body.angular_velocity = Vec3::zeros();
            }
        }
    }
    /// Apply uniform gravity to all dynamic, non-sleeping bodies for time step `dt`.
    ///
    /// This adds `mass * g` to each body's force accumulator.
    pub fn apply_gravity(&mut self, g: [f64; 3], _dt: f64) {
        let gravity = Vec3::new(g[0], g[1], g[2]);
        for (_, body) in self.iter_mut() {
            if body.body_type != crate::body::BodyType::Static
                && body.state != crate::body::BodyState::Sleeping
            {
                body.force_accumulator += gravity * body.mass;
            }
        }
    }
    /// Integrate accumulated forces into velocities for all dynamic bodies.
    pub fn step_velocities(&mut self, dt: f64) {
        for (_, body) in self.iter_mut() {
            if body.body_type != crate::body::BodyType::Static
                && body.state != crate::body::BodyState::Sleeping
            {
                if body.inverse_mass > 0.0 {
                    body.velocity += body.force_accumulator * (body.inverse_mass * dt);
                }
                let inv_inertia_diag = Vec3::new(
                    body.world_inverse_inertia[(0, 0)],
                    body.world_inverse_inertia[(1, 1)],
                    body.world_inverse_inertia[(2, 2)],
                );
                body.angular_velocity += Vec3::new(
                    body.torque_accumulator.x * inv_inertia_diag.x * dt,
                    body.torque_accumulator.y * inv_inertia_diag.y * dt,
                    body.torque_accumulator.z * inv_inertia_diag.z * dt,
                );
            }
        }
    }
    /// Integrate velocities into positions for all dynamic bodies.
    pub fn step_positions(&mut self, dt: f64) {
        for (_, body) in self.iter_mut() {
            if body.body_type != crate::body::BodyType::Static
                && body.state != crate::body::BodyState::Sleeping
            {
                body.integrate_velocity(dt);
            }
        }
    }
    /// Zero all force and torque accumulators.
    pub fn clear_forces(&mut self) {
        for (_, body) in self.iter_mut() {
            body.force_accumulator = Vec3::zeros();
            body.torque_accumulator = Vec3::zeros();
        }
    }
    /// Attach a collider to a body and update the body's mass and inertia tensor
    /// from the collider's shape geometry and density.
    pub fn attach_collider_with_auto_inertia(
        &mut self,
        body_handle: BodyHandle,
        collider: &Collider,
    ) {
        if let Some(body) = self.get_mut(body_handle) {
            let props = collider.mass_properties();
            body.set_mass_properties(&props);
        }
    }
    /// Compute aggregate statistics about the body set.
    pub fn statistics(&self) -> BodySetStatistics {
        let mut stats = BodySetStatistics {
            total_count: self.len(),
            ..Default::default()
        };
        for (_, body) in self.iter() {
            match body.body_type {
                crate::body::BodyType::Dynamic => stats.dynamic_count += 1,
                crate::body::BodyType::Static => stats.static_count += 1,
                crate::body::BodyType::Kinematic => stats.kinematic_count += 1,
            }
            match body.state {
                crate::body::BodyState::Active => stats.active_count += 1,
                crate::body::BodyState::Sleeping => stats.sleeping_count += 1,
            }
            stats.total_mass += body.mass;
            let speed = body.velocity.norm();
            if speed > stats.max_speed {
                stats.max_speed = speed;
            }
            stats.total_kinetic_energy += 0.5 * body.mass * speed * speed;
        }
        if stats.dynamic_count > 0 {
            stats.avg_mass = stats.total_mass / stats.dynamic_count as f64;
        }
        stats
    }
    /// Export all body positions as a flat array `[x0, y0, z0, x1, y1, z1, ...]`.
    pub fn export_positions(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.len() * 3);
        for (_, body) in self.iter() {
            let p = body.transform.position;
            out.push(p.x);
            out.push(p.y);
            out.push(p.z);
        }
        out
    }
    /// Export all body velocities as a flat array.
    pub fn export_velocities(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.len() * 3);
        for (_, body) in self.iter() {
            let v = body.velocity;
            out.push(v.x);
            out.push(v.y);
            out.push(v.z);
        }
        out
    }
    /// Export all body masses as a flat array.
    pub fn export_masses(&self) -> Vec<f64> {
        self.iter().map(|(_, b)| b.mass).collect()
    }
    /// Import positions from a flat array, updating bodies in iteration order.
    pub fn import_positions(&mut self, data: &[f64]) {
        let mut idx = 0;
        for (_, body) in self.iter_mut() {
            if idx + 2 < data.len() {
                body.transform.position = Vec3::new(data[idx], data[idx + 1], data[idx + 2]);
                idx += 3;
            }
        }
    }
    /// Import velocities from a flat array, updating bodies in iteration order.
    pub fn import_velocities(&mut self, data: &[f64]) {
        let mut idx = 0;
        for (_, body) in self.iter_mut() {
            if idx + 2 < data.len() {
                body.velocity = Vec3::new(data[idx], data[idx + 1], data[idx + 2]);
                idx += 3;
            }
        }
    }
    /// Compute the total linear and angular momentum of all dynamic bodies.
    ///
    /// Returns `(linear_momentum, angular_momentum)` as `[f64; 3]` arrays.
    ///
    /// * Linear momentum:  `p = Σ m_i * v_i`
    /// * Angular momentum: `L = Σ m_i * (r_i × v_i)` where `r_i` is the
    ///   world-space position of body `i`.
    pub fn compute_total_momentum(&self) -> ([f64; 3], [f64; 3]) {
        let mut linear = [0.0f64; 3];
        let mut angular = [0.0f64; 3];
        for (_, body) in self.iter() {
            if body.inverse_mass <= 0.0 {
                continue;
            }
            let m = body.mass;
            let v = [body.velocity.x, body.velocity.y, body.velocity.z];
            let r = [
                body.transform.position.x,
                body.transform.position.y,
                body.transform.position.z,
            ];
            linear[0] += m * v[0];
            linear[1] += m * v[1];
            linear[2] += m * v[2];
            angular[0] += m * (r[1] * v[2] - r[2] * v[1]);
            angular[1] += m * (r[2] * v[0] - r[0] * v[2]);
            angular[2] += m * (r[0] * v[1] - r[1] * v[0]);
        }
        (linear, angular)
    }
    /// Compute the mass-weighted centre of mass of all dynamic bodies.
    ///
    /// Returns `Some([x, y, z])` or `None` if there are no dynamic bodies with
    /// positive mass.
    pub fn compute_center_of_mass(&self) -> Option<[f64; 3]> {
        let mut total_mass = 0.0f64;
        let mut com = [0.0f64; 3];
        for (_, body) in self.iter() {
            if body.inverse_mass <= 0.0 || body.mass <= 0.0 {
                continue;
            }
            let m = body.mass;
            let p = body.transform.position;
            com[0] += m * p.x;
            com[1] += m * p.y;
            com[2] += m * p.z;
            total_mass += m;
        }
        if total_mass < 1e-30 {
            return None;
        }
        Some([
            com[0] / total_mass,
            com[1] / total_mass,
            com[2] / total_mass,
        ])
    }
    /// Apply a radial explosion impulse to all dynamic bodies within `blast_radius`.
    ///
    /// The impulse magnitude applied to body `i` is:
    ///
    /// ```text
    /// J_i = peak_impulse * (1 - d_i / blast_radius)
    /// ```
    ///
    /// where `d_i` is the distance from `epicentre` to the body's centre of
    /// mass.  Bodies at the epicentre receive `peak_impulse`; bodies at or
    /// beyond `blast_radius` receive nothing.  The impulse direction is the
    /// unit vector from `epicentre` to each body.
    pub fn apply_explosion_impulse(
        &mut self,
        epicentre: [f64; 3],
        blast_radius: f64,
        peak_impulse: f64,
    ) {
        if blast_radius < 1e-30 {
            return;
        }
        let impulses: Vec<(BodyHandle, Vec3)> = self
            .iter()
            .filter_map(|(h, body)| {
                if body.inverse_mass <= 0.0 {
                    return None;
                }
                let p = body.transform.position;
                let dx = p.x - epicentre[0];
                let dy = p.y - epicentre[1];
                let dz = p.z - epicentre[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist >= blast_radius || dist < 1e-30 {
                    return None;
                }
                let falloff = 1.0 - dist / blast_radius;
                let mag = peak_impulse * falloff;
                let dir = Vec3::new(dx / dist, dy / dist, dz / dist);
                Some((h, dir * mag))
            })
            .collect();
        for (h, imp) in impulses {
            if let Some(body) = self.get_mut(h) {
                body.apply_impulse(imp);
            }
        }
    }
    /// Apply new transforms (position + velocity) to a batch of bodies in one
    /// pass — SIMD-friendly layout.
    ///
    /// `transforms` is a flat slice with layout:
    /// `[px, py, pz, vx, vy, vz, ...]` repeated for each handle in `handles`.
    ///
    /// Bodies for which the handle is invalid are silently skipped.  The slice
    /// length must be at least `6 * handles.len()`; extra elements are ignored.
    pub fn batch_update_transforms(&mut self, handles: &[BodyHandle], transforms: &[f64]) {
        let stride = 6usize;
        for (i, &h) in handles.iter().enumerate() {
            let base = i * stride;
            if base + 5 >= transforms.len() {
                break;
            }
            if let Some(body) = self.get_mut(h) {
                body.transform.position =
                    Vec3::new(transforms[base], transforms[base + 1], transforms[base + 2]);
                body.velocity = Vec3::new(
                    transforms[base + 3],
                    transforms[base + 4],
                    transforms[base + 5],
                );
            }
        }
    }
}
/// Aggregate statistics for a `RigidBodySet`.
#[derive(Debug, Clone, Default)]
pub struct BodySetStatistics {
    /// Total number of bodies.
    pub total_count: usize,
    /// Number of dynamic bodies.
    pub dynamic_count: usize,
    /// Number of static bodies.
    pub static_count: usize,
    /// Number of kinematic bodies.
    pub kinematic_count: usize,
    /// Number of active bodies.
    pub active_count: usize,
    /// Number of sleeping bodies.
    pub sleeping_count: usize,
    /// Total mass of all bodies.
    pub total_mass: f64,
    /// Average mass of dynamic bodies.
    pub avg_mass: f64,
    /// Maximum speed of any body.
    pub max_speed: f64,
    /// Total kinetic energy (linear only).
    pub total_kinetic_energy: f64,
}
/// Sub-step interpolation helper for render-time smoothing.
///
/// Stores the "previous" positions/velocities of all bodies so that
/// intermediate positions can be linearly interpolated using a sub-step
/// blend factor `alpha ∈ [0, 1]`.
#[derive(Debug, Default, Clone)]
pub struct BodyInterpolator {
    /// Map from body handle (index) to previous position.
    pub(super) prev_positions: std::collections::HashMap<u32, oxiphysics_core::math::Vec3>,
    /// Map from body handle (index) to previous velocity.
    pub(super) prev_velocities: std::collections::HashMap<u32, oxiphysics_core::math::Vec3>,
}
impl BodyInterpolator {
    /// Create a new empty interpolator.
    pub fn new() -> Self {
        Self::default()
    }
    /// Snapshot the current positions and velocities of all bodies as the
    /// "previous" state for next interpolation.
    pub fn snapshot_previous(&mut self, set: &RigidBodySet) {
        self.prev_positions.clear();
        self.prev_velocities.clear();
        for (h, body) in set.iter() {
            self.prev_positions.insert(h.index, body.transform.position);
            self.prev_velocities.insert(h.index, body.velocity);
        }
    }
    /// Interpolate the position of a body at sub-step blend factor `alpha`.
    ///
    /// `alpha = 0` → previous position; `alpha = 1` → current position.
    /// Falls back to the current position if no previous snapshot exists.
    pub fn interpolate_position(
        &self,
        handle: oxiphysics_core::BodyHandle,
        body: &crate::body::RigidBody,
        alpha: f64,
    ) -> oxiphysics_core::math::Vec3 {
        let curr = body.transform.position;
        if let Some(&prev) = self.prev_positions.get(&handle.index) {
            prev + (curr - prev) * alpha
        } else {
            curr
        }
    }
    /// Interpolate the velocity of a body at sub-step blend factor `alpha`.
    pub fn interpolate_velocity(
        &self,
        handle: oxiphysics_core::BodyHandle,
        body: &crate::body::RigidBody,
        alpha: f64,
    ) -> oxiphysics_core::math::Vec3 {
        let curr = body.velocity;
        if let Some(&prev) = self.prev_velocities.get(&handle.index) {
            prev + (curr - prev) * alpha
        } else {
            curr
        }
    }
    /// Clear all stored snapshots.
    pub fn clear(&mut self) {
        self.prev_positions.clear();
        self.prev_velocities.clear();
    }
    /// Number of bodies currently tracked.
    pub fn len(&self) -> usize {
        self.prev_positions.len()
    }
    /// Whether no bodies are tracked.
    pub fn is_empty(&self) -> bool {
        self.prev_positions.is_empty()
    }
}
/// Uniform-grid spatial hash for O(1) broad-phase neighbour queries.
///
/// Bodies are mapped to integer cells of side `cell_size`.  All bodies that
/// fall into the same cell — or an immediately adjacent cell for radius queries
/// — are returned as candidates.
#[derive(Debug, Clone)]
pub struct SpatialHash {
    /// Side-length of each cubic cell.
    pub(super) cell_size: f64,
    /// Map from (ix, iy, iz) hash to list of body slot indices.
    pub(super) cells: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
    /// All (id, pos) pairs inserted (used for radius queries).
    pub(super) entries: Vec<(usize, [f64; 3])>,
}
impl SpatialHash {
    /// Create a new spatial hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: std::collections::HashMap::new(),
            entries: Vec::new(),
        }
    }
    /// Insert a body slot index at world position `pos`.
    pub fn insert(&mut self, id: usize, pos: [f64; 3]) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push(id);
        self.entries.push((id, pos));
    }
    /// Remove all entries.
    pub fn clear(&mut self) {
        self.cells.clear();
        self.entries.clear();
    }
    /// Return all IDs in the cell that contains `pos`.
    pub fn query_cell(&self, pos: [f64; 3]) -> Vec<usize> {
        let key = self.cell_key(pos);
        self.cells.get(&key).cloned().unwrap_or_default()
    }
    /// Return all IDs within `radius` of `center`, using the hash for
    /// broad-phase filtering then exact distance check.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        let extra = (radius / self.cell_size).ceil() as i64 + 1;
        let base = self.cell_key(center);
        let mut candidates = std::collections::HashSet::new();
        for dx in -extra..=extra {
            for dy in -extra..=extra {
                for dz in -extra..=extra {
                    let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                    if let Some(ids) = self.cells.get(&key) {
                        for &id in ids {
                            candidates.insert(id);
                        }
                    }
                }
            }
        }
        candidates
            .into_iter()
            .filter(|&id| {
                if let Some(&(_, pos)) = self.entries.iter().find(|(eid, _)| *eid == id) {
                    let dx = pos[0] - center[0];
                    let dy = pos[1] - center[1];
                    let dz = pos[2] - center[2];
                    dx * dx + dy * dy + dz * dz <= r2
                } else {
                    false
                }
            })
            .collect()
    }
    fn cell_key(&self, pos: [f64; 3]) -> (i64, i64, i64) {
        (
            (pos[0] / self.cell_size).floor() as i64,
            (pos[1] / self.cell_size).floor() as i64,
            (pos[2] / self.cell_size).floor() as i64,
        )
    }
    /// Total number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether there are no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
