// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced spatial queries.
//!
//! Implements a dynamic AABB tree, overlapping-pair queries (brute-force,
//! sweep-and-prune, AABB tree), frustum culling, sphere/capsule/point/convex
//! sweep queries, contact pair filtering by layer mask, and a batched query
//! pipeline.

// ---------------------------------------------------------------------------
// AABB type
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Clone, Debug, PartialEq)]
pub struct QueryAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl QueryAabb {
    /// Construct from min and max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        QueryAabb { min, max }
    }

    /// Expand each side by `margin`.
    pub fn fatten(&self, margin: f64) -> Self {
        QueryAabb {
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

    /// Surface area of the AABB.
    pub fn surface_area(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Union of two AABBs.
    pub fn union(&self, other: &QueryAabb) -> QueryAabb {
        QueryAabb {
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

    /// Test if two AABBs overlap.
    pub fn overlaps(&self, other: &QueryAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Test if this AABB fully contains `other` (touching faces count as contained).
    pub fn contains_aabb(&self, other: &QueryAabb) -> bool {
        self.min
            .iter()
            .zip(other.min.iter())
            .all(|(outer, inner)| outer <= inner)
            && self
                .max
                .iter()
                .zip(other.max.iter())
                .all(|(outer, inner)| inner <= outer)
    }

    /// Test if this AABB contains a point.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Squared distance from a point to the AABB (0 if inside).
    pub fn point_dist_sq(&self, p: [f64; 3]) -> f64 {
        p.iter()
            .zip(self.min.iter())
            .zip(self.max.iter())
            .map(|((&pi, &mini), &maxi)| {
                if pi < mini {
                    (mini - pi) * (mini - pi)
                } else if pi > maxi {
                    (pi - maxi) * (pi - maxi)
                } else {
                    0.0
                }
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Helper functions (public)
// ---------------------------------------------------------------------------

/// Test if two AABBs overlap.
pub fn aabb_overlap(a: &QueryAabb, b: &QueryAabb) -> bool {
    a.overlaps(b)
}

/// Test if an AABB passes a frustum culling test (6 plane tests).
///
/// `planes` is a slice of (normal, d) pairs where plane equation = n·x + d = 0.
/// Returns true if the AABB is at least partially inside the frustum.
pub fn frustum_aabb_test(aabb: &QueryAabb, planes: &[([f64; 3], f64)]) -> bool {
    for &(n, d) in planes {
        // P-vertex: corner furthest in direction n
        let px = if n[0] > 0.0 { aabb.max[0] } else { aabb.min[0] };
        let py = if n[1] > 0.0 { aabb.max[1] } else { aabb.min[1] };
        let pz = if n[2] > 0.0 { aabb.max[2] } else { aabb.min[2] };
        let dist = n[0] * px + n[1] * py + n[2] * pz + d;
        if dist < 0.0 {
            return false;
        } // AABB fully outside this plane
    }
    true
}

/// Compute the swept AABB of an object moving from `aabb` by `motion`.
pub fn sweep_aabb_motion(aabb: &QueryAabb, motion: [f64; 3]) -> QueryAabb {
    let mut min = aabb.min;
    let mut max = aabb.max;
    for ((mn, mx), &m) in min.iter_mut().zip(max.iter_mut()).zip(motion.iter()) {
        if m < 0.0 {
            *mn += m;
        } else {
            *mx += m;
        }
    }
    QueryAabb { min, max }
}

/// Rebalance an AABB tree.
///
/// Insertions, removals and updates already keep the tree height-balanced with
/// AVL-style rotations. This additionally rebuilds the whole hierarchy from its
/// leaves with a median split (see [`AabbTree::rebuild`]), which restores a
/// tight, perfectly balanced tree after a long series of incremental changes.
/// Leaf handles returned by [`AabbTree::insert`] stay valid.
pub fn rebalance_tree(tree: &mut AabbTree) {
    tree.rebuild();
}

// ---------------------------------------------------------------------------
// Dynamic AABB tree node
// ---------------------------------------------------------------------------

/// A node in the dynamic AABB tree.
#[derive(Clone, Debug)]
pub struct DynamicNode {
    /// Fat AABB (enlarged by margin).
    pub fat_aabb: QueryAabb,
    /// Tight AABB (actual object bounds).
    pub aabb: QueryAabb,
    /// Parent node index, or None for root.
    pub parent: Option<usize>,
    /// Children node indices (None for leaves).
    pub children: [Option<usize>; 2],
    /// Height in the tree (0 = leaf).
    pub height: i32,
    /// User data (shape ID).
    pub user_data: usize,
    /// True if this is a leaf node.
    pub is_leaf: bool,
}

impl DynamicNode {
    /// Create a leaf node.
    pub fn leaf(aabb: QueryAabb, margin: f64, user_data: usize) -> Self {
        let fat_aabb = aabb.fatten(margin);
        DynamicNode {
            fat_aabb,
            aabb,
            parent: None,
            children: [None, None],
            height: 0,
            user_data,
            is_leaf: true,
        }
    }

    /// Create an internal node.
    pub fn internal(fat_aabb: QueryAabb) -> Self {
        DynamicNode {
            fat_aabb: fat_aabb.clone(),
            aabb: fat_aabb,
            parent: None,
            children: [None, None],
            height: 1,
            user_data: usize::MAX,
            is_leaf: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Dynamic AABB tree
// ---------------------------------------------------------------------------

/// Dynamic AABB tree for broad-phase collision detection.
///
/// Every leaf stores its tight AABB and a *fat* AABB enlarged by `margin`;
/// queries test the fat boxes. Internal nodes store the union of their
/// children's fat boxes. [`AabbTree::insert`], [`AabbTree::remove`] and
/// [`AabbTree::update`] refit every ancestor and keep the tree close to
/// height-balanced with AVL-style rotations; new leaves are placed with a
/// branch-and-bound surface-area-heuristic search.
pub struct AabbTree {
    /// Node storage.
    pub nodes: Vec<DynamicNode>,
    /// Root node index.
    pub root: Option<usize>,
    /// Free node list (indices of removed nodes).
    free_nodes: Vec<usize>,
    /// Fat AABB margin.
    pub margin: f64,
}

impl AabbTree {
    /// Construct an empty AABB tree with the given fat-AABB `margin`.
    ///
    /// A negative or NaN margin is treated as zero.
    pub fn new(margin: f64) -> Self {
        AabbTree {
            nodes: Vec::new(),
            root: None,
            free_nodes: Vec::new(),
            margin: margin.max(0.0),
        }
    }

    fn alloc_node(&mut self, node: DynamicNode) -> usize {
        if let Some(idx) = self.free_nodes.pop() {
            self.nodes[idx] = node;
            idx
        } else {
            self.nodes.push(node);
            self.nodes.len() - 1
        }
    }

    /// True if `idx` names a live leaf (not a free slot and not an internal node).
    fn is_live_leaf(&self, idx: usize) -> bool {
        self.nodes
            .get(idx)
            .is_some_and(|node| node.is_leaf && node.height == 0)
    }

    /// Return a node to the free list. Free slots are marked with `height == -1`.
    fn release_node(&mut self, idx: usize) {
        let node = &mut self.nodes[idx];
        node.parent = None;
        node.children = [None, None];
        node.height = -1;
        node.is_leaf = false;
        node.user_data = usize::MAX;
        self.free_nodes.push(idx);
    }

    /// Both children of an internal node, or `None` for a leaf or free slot.
    fn children_of(&self, idx: usize) -> Option<(usize, usize)> {
        let node = self.nodes.get(idx)?;
        if node.is_leaf {
            return None;
        }
        match node.children {
            [Some(left), Some(right)] => Some((left, right)),
            _ => None,
        }
    }

    /// Replace the child `old` of `parent` by `new`.
    fn replace_child(&mut self, parent: usize, old: usize, new: usize) {
        for slot in self.nodes[parent].children.iter_mut() {
            if *slot == Some(old) {
                *slot = Some(new);
            }
        }
    }

    /// Recompute the boxes and the height of an internal node from its children.
    fn refit_node(&mut self, idx: usize) {
        let Some((left, right)) = self.children_of(idx) else {
            return;
        };
        let fat = self.nodes[left].fat_aabb.union(&self.nodes[right].fat_aabb);
        let height = 1 + self.nodes[left].height.max(self.nodes[right].height);
        let node = &mut self.nodes[idx];
        node.aabb = fat.clone();
        node.fat_aabb = fat;
        node.height = height;
    }

    /// Rotate the child `up` of `top` above `top`, moving the shorter of the
    /// grandchildren below `top` (Box2D-style AVL rotation). Returns the node
    /// that now roots the subtree. Leaves are never moved.
    fn rotate_up(&mut self, top: usize, up: usize) -> usize {
        let Some((first, _)) = self.children_of(top) else {
            return top;
        };
        let up_slot = if first == up { 0 } else { 1 };
        let Some((grand_f, grand_g)) = self.children_of(up) else {
            return top;
        };

        // `up` takes the place of `top` below `top`'s old parent.
        let old_parent = self.nodes[top].parent;
        self.nodes[up].children[0] = Some(top);
        self.nodes[up].parent = old_parent;
        self.nodes[top].parent = Some(up);
        match old_parent {
            Some(p) => self.replace_child(p, top, up),
            None => self.root = Some(up),
        }

        // The taller grandchild stays with `up`, the shorter one moves to `top`.
        let (keep, give) = if self.nodes[grand_f].height > self.nodes[grand_g].height {
            (grand_f, grand_g)
        } else {
            (grand_g, grand_f)
        };
        self.nodes[up].children[1] = Some(keep);
        self.nodes[top].children[up_slot] = Some(give);
        self.nodes[give].parent = Some(top);
        self.refit_node(top);
        self.refit_node(up);
        up
    }

    /// Rotate if the subtrees below `idx` differ in height by more than one.
    /// Returns the node that now roots the subtree.
    ///
    /// A single rotation is applied per call, as in Box2D's dynamic tree. A
    /// leaf inserted next to a tall subtree can leave a residual imbalance
    /// below the rotated node; that costs query speed but never correctness.
    fn balance(&mut self, idx: usize) -> usize {
        let Some((left, right)) = self.children_of(idx) else {
            return idx;
        };
        let difference = self.nodes[right].height - self.nodes[left].height;
        if difference > 1 {
            self.rotate_up(idx, right)
        } else if difference < -1 {
            self.rotate_up(idx, left)
        } else {
            idx
        }
    }

    /// Walk from `start` to the root, rebalancing and refitting every ancestor.
    fn rebalance_and_refit_upwards(&mut self, start: usize) {
        let mut cursor = Some(start);
        while let Some(idx) = cursor {
            let idx = self.balance(idx);
            self.refit_node(idx);
            cursor = self.nodes[idx].parent;
        }
    }

    /// Branch-and-bound surface-area-heuristic search for the node below which a
    /// new leaf with fat box `leaf_fat` should be placed as a sibling.
    ///
    /// The cost of a candidate is the area of its union with the leaf plus the
    /// area growth ("inherited cost") it forces on every ancestor; subtrees that
    /// cannot beat the best candidate are pruned.
    fn find_best_sibling(&self, root: usize, leaf_fat: &QueryAabb) -> usize {
        let leaf_area = leaf_fat.surface_area();
        let mut best = root;
        let mut best_cost = self.nodes[root].fat_aabb.union(leaf_fat).surface_area();
        let mut stack: Vec<(usize, f64)> = vec![(root, 0.0)];
        while let Some((idx, inherited)) = stack.pop() {
            let node = &self.nodes[idx];
            let combined_area = node.fat_aabb.union(leaf_fat).surface_area();
            let cost = combined_area + inherited;
            if cost < best_cost {
                best_cost = cost;
                best = idx;
            }
            if let Some((left, right)) = self.children_of(idx) {
                let child_inherited = inherited + combined_area - node.fat_aabb.surface_area();
                if leaf_area + child_inherited < best_cost {
                    stack.push((left, child_inherited));
                    stack.push((right, child_inherited));
                }
            }
        }
        best
    }

    /// Link the (detached) leaf `leaf` into the tree and restore the ancestors'
    /// boxes, heights and balance.
    fn attach_leaf(&mut self, leaf: usize) {
        self.nodes[leaf].parent = None;
        let Some(root) = self.root else {
            self.root = Some(leaf);
            return;
        };
        let leaf_fat = self.nodes[leaf].fat_aabb.clone();
        let sibling = self.find_best_sibling(root, &leaf_fat);
        let old_parent = self.nodes[sibling].parent;
        let parent_fat = self.nodes[sibling].fat_aabb.union(&leaf_fat);
        let new_parent = self.alloc_node(DynamicNode::internal(parent_fat));
        self.nodes[new_parent].parent = old_parent;
        self.nodes[new_parent].children = [Some(sibling), Some(leaf)];
        self.nodes[new_parent].height = 1 + self.nodes[sibling].height.max(0);
        self.nodes[sibling].parent = Some(new_parent);
        self.nodes[leaf].parent = Some(new_parent);
        match old_parent {
            Some(parent) => self.replace_child(parent, sibling, new_parent),
            None => self.root = Some(new_parent),
        }
        self.rebalance_and_refit_upwards(new_parent);
    }

    /// Unlink the leaf `leaf` from the tree (its slot stays allocated) and
    /// restore the ancestors' boxes, heights and balance.
    fn detach_leaf(&mut self, leaf: usize) {
        if self.root == Some(leaf) {
            self.root = None;
            self.nodes[leaf].parent = None;
            return;
        }
        let Some(parent) = self.nodes[leaf].parent else {
            return;
        };
        let sibling = match self.nodes[parent].children {
            [Some(a), Some(b)] => {
                if a == leaf {
                    b
                } else {
                    a
                }
            }
            _ => return,
        };
        let grandparent = self.nodes[parent].parent;
        self.nodes[sibling].parent = grandparent;
        match grandparent {
            Some(g) => self.replace_child(g, parent, sibling),
            None => self.root = Some(sibling),
        }
        self.release_node(parent);
        self.nodes[leaf].parent = None;
        if let Some(g) = grandparent {
            self.rebalance_and_refit_upwards(g);
        }
    }

    /// Insert an AABB with user data, returns the leaf node index.
    ///
    /// The leaf is stored with a fat AABB (enlarged by the tree margin); all
    /// ancestors are refitted and the tree is rebalanced, so the new leaf is
    /// immediately reachable by [`AabbTree::query_aabb`] and
    /// [`AabbTree::raycast`].
    pub fn insert(&mut self, aabb: QueryAabb, user_data: usize) -> usize {
        let leaf_idx = self.alloc_node(DynamicNode::leaf(aabb, self.margin, user_data));
        self.attach_leaf(leaf_idx);
        leaf_idx
    }

    /// Remove a leaf node by index.
    ///
    /// An index that does not name a live leaf (out of range, already removed,
    /// or an internal node) is ignored.
    pub fn remove(&mut self, leaf_idx: usize) {
        if !self.is_live_leaf(leaf_idx) {
            return;
        }
        self.detach_leaf(leaf_idx);
        self.release_node(leaf_idx);
    }

    /// Update a leaf node's AABB.
    ///
    /// If the new AABB is still contained in the leaf's fat AABB only the tight
    /// box is replaced and the tree is untouched. Otherwise the leaf is
    /// re-inserted with a fresh fat AABB; the leaf index (the handle returned by
    /// [`AabbTree::insert`]) stays valid in both cases. An index that does not
    /// name a live leaf is ignored.
    pub fn update(&mut self, leaf_idx: usize, new_aabb: QueryAabb) {
        if !self.is_live_leaf(leaf_idx) {
            return;
        }
        if self.nodes[leaf_idx].fat_aabb.contains_aabb(&new_aabb) {
            self.nodes[leaf_idx].aabb = new_aabb;
            return;
        }
        self.detach_leaf(leaf_idx);
        let fat = new_aabb.fatten(self.margin);
        let leaf = &mut self.nodes[leaf_idx];
        leaf.aabb = new_aabb;
        leaf.fat_aabb = fat;
        self.attach_leaf(leaf_idx);
    }

    /// Height of the tree (0 for an empty tree or a single leaf).
    pub fn height(&self) -> i32 {
        self.root.map_or(0, |root| self.nodes[root].height)
    }

    /// Rebuild the whole hierarchy from the live leaves with a median split along
    /// the axis of largest centre spread.
    ///
    /// The result is a perfectly balanced tree; leaf indices are preserved,
    /// internal nodes are reallocated.
    pub fn rebuild(&mut self) {
        let mut leaves: Vec<usize> = Vec::new();
        for idx in 0..self.nodes.len() {
            if self.is_live_leaf(idx) {
                leaves.push(idx);
            } else if self.nodes[idx].height > 0 {
                self.release_node(idx);
            }
        }
        self.root = None;
        if leaves.is_empty() {
            return;
        }
        let root = self.build_subtree(&mut leaves);
        self.nodes[root].parent = None;
        self.root = Some(root);
    }

    /// Build a balanced subtree over `leaves` and return its root.
    fn build_subtree(&mut self, leaves: &mut [usize]) -> usize {
        if let [only] = *leaves {
            self.nodes[only].parent = None;
            return only;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &leaf in leaves.iter() {
            let centre = self.nodes[leaf].fat_aabb.center();
            for ((l, h), c) in lo.iter_mut().zip(hi.iter_mut()).zip(centre) {
                *l = l.min(c);
                *h = h.max(c);
            }
        }
        let axis = (0..3)
            .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
            .unwrap_or(0);
        let nodes = &self.nodes;
        leaves.sort_by(|&a, &b| {
            nodes[a].fat_aabb.center()[axis].total_cmp(&nodes[b].fat_aabb.center()[axis])
        });
        let (left_half, right_half) = leaves.split_at_mut(leaves.len() / 2);
        let left = self.build_subtree(left_half);
        let right = self.build_subtree(right_half);
        let fat = self.nodes[left].fat_aabb.union(&self.nodes[right].fat_aabb);
        let parent = self.alloc_node(DynamicNode::internal(fat));
        self.nodes[parent].children = [Some(left), Some(right)];
        self.nodes[parent].height = 1 + self.nodes[left].height.max(self.nodes[right].height);
        self.nodes[left].parent = Some(parent);
        self.nodes[right].parent = Some(parent);
        parent
    }

    /// Check the structural invariants of the tree.
    ///
    /// Verifies parent/child links, node heights, that every internal fat box
    /// contains its children's fat boxes, that every leaf's fat box contains its
    /// tight box, and that every allocated slot is either reachable from the root
    /// or on the free list. (Rotations keep the tree close to balanced but, as
    /// in other surface-area-heuristic dynamic trees, do not enforce a strict
    /// height difference of at most one.)
    pub fn validate(&self) -> bool {
        let Some(root) = self.root else {
            return self.free_nodes.len() == self.nodes.len();
        };
        if self
            .nodes
            .get(root)
            .is_none_or(|node| node.parent.is_some())
        {
            return false;
        }
        let mut reachable = 0usize;
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            reachable += 1;
            if reachable > self.nodes.len() {
                return false;
            }
            let node = &self.nodes[idx];
            if node.is_leaf {
                if node.height != 0
                    || node.children != [None, None]
                    || !node.fat_aabb.contains_aabb(&node.aabb)
                {
                    return false;
                }
                continue;
            }
            let Some((left, right)) = self.children_of(idx) else {
                return false;
            };
            let (Some(left_node), Some(right_node)) = (self.nodes.get(left), self.nodes.get(right))
            else {
                return false;
            };
            if left_node.parent != Some(idx) || right_node.parent != Some(idx) {
                return false;
            }
            if node.height != 1 + left_node.height.max(right_node.height) {
                return false;
            }
            if !node.fat_aabb.contains_aabb(&left_node.fat_aabb)
                || !node.fat_aabb.contains_aabb(&right_node.fat_aabb)
            {
                return false;
            }
            stack.push(left);
            stack.push(right);
        }
        reachable + self.free_nodes.len() == self.nodes.len()
    }

    /// Query all leaf nodes whose fat AABB overlaps `query_aabb`.
    pub fn query_aabb(&self, query_aabb: &QueryAabb) -> Vec<usize> {
        let mut result = Vec::new();
        let Some(root) = self.root else {
            return result;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.fat_aabb.overlaps(query_aabb) {
                continue;
            }
            if node.is_leaf {
                result.push(node.user_data);
            } else {
                for c in 0..2 {
                    if let Some(child) = node.children[c] {
                        stack.push(child);
                    }
                }
            }
        }
        result
    }

    /// Raycast: find all leaf nodes whose fat AABB is hit by ray (origin, dir, max_t).
    pub fn raycast(&self, origin: [f64; 3], dir: [f64; 3], max_t: f64) -> Vec<usize> {
        let mut result = Vec::new();
        let Some(root) = self.root else {
            return result;
        };
        let mut stack = vec![root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !ray_aabb_hit(&node.fat_aabb, origin, dir, max_t) {
                continue;
            }
            if node.is_leaf {
                result.push(node.user_data);
            } else {
                for c in 0..2 {
                    if let Some(child) = node.children[c] {
                        stack.push(child);
                    }
                }
            }
        }
        result
    }

    /// Number of live nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len().saturating_sub(self.free_nodes.len())
    }
}

fn ray_aabb_hit(aabb: &QueryAabb, origin: [f64; 3], dir: [f64; 3], max_t: f64) -> bool {
    let mut t_min = 0.0_f64;
    let mut t_max = max_t;
    for ((&o, &d), (&mn, &mx)) in origin
        .iter()
        .zip(dir.iter())
        .zip(aabb.min.iter().zip(aabb.max.iter()))
    {
        if d.abs() < 1e-14 {
            if o < mn || o > mx {
                return false;
            }
        } else {
            let inv = 1.0 / d;
            let t1 = (mn - o) * inv;
            let t2 = (mx - o) * inv;
            let ta = t1.min(t2);
            let tb = t1.max(t2);
            t_min = t_min.max(ta);
            t_max = t_max.min(tb);
            if t_min > t_max {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Pair query
// ---------------------------------------------------------------------------

/// An overlapping pair of shape IDs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShapePair {
    /// First shape ID.
    pub a: usize,
    /// Second shape ID.
    pub b: usize,
}

impl ShapePair {
    /// Create a new pair (a < b always).
    pub fn new(a: usize, b: usize) -> Self {
        if a < b {
            ShapePair { a, b }
        } else {
            ShapePair { a: b, b: a }
        }
    }
}

/// Overlapping pair detection using multiple algorithms.
pub struct PairQuery;

impl PairQuery {
    /// Brute-force O(n²) pair detection.
    pub fn brute_force(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let mut pairs = Vec::new();
        for (i, a) in aabbs.iter().enumerate() {
            for (j, b) in aabbs[i + 1..].iter().enumerate() {
                if a.overlaps(b) {
                    pairs.push(ShapePair::new(i, i + 1 + j));
                }
            }
        }
        pairs
    }

    /// Sweep-and-prune along X axis.
    pub fn sweep_and_prune(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let mut sorted: Vec<usize> = (0..aabbs.len()).collect();
        sorted.sort_unstable_by(|&a, &b| {
            aabbs[a].min[0]
                .partial_cmp(&aabbs[b].min[0])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut pairs = Vec::new();
        for (i, &a) in sorted.iter().enumerate() {
            for &b in &sorted[i + 1..] {
                if aabbs[b].min[0] > aabbs[a].max[0] {
                    break;
                }
                if aabbs[a].overlaps(&aabbs[b]) {
                    pairs.push(ShapePair::new(a, b));
                }
            }
        }
        pairs
    }

    /// AABB tree-based pair detection.
    ///
    /// The tree (whose boxes are fattened) only supplies candidates; a pair is
    /// reported when the tight AABBs overlap, so the result equals
    /// [`PairQuery::brute_force`].
    pub fn aabb_tree(aabbs: &[QueryAabb]) -> Vec<ShapePair> {
        let mut tree = AabbTree::new(0.01);
        for (i, aabb) in aabbs.iter().enumerate() {
            tree.insert(aabb.clone(), i);
        }
        let mut pairs = Vec::new();
        for (i, aabb) in aabbs.iter().enumerate() {
            let hits = tree.query_aabb(aabb);
            for h in hits {
                if h > i && aabbs.get(h).is_some_and(|other| other.overlaps(aabb)) {
                    pairs.push(ShapePair::new(i, h));
                }
            }
        }
        // Deduplicate
        pairs.sort_unstable_by_key(|p| (p.a, p.b));
        pairs.dedup();
        pairs
    }
}

// ---------------------------------------------------------------------------
// Frustum culling
// ---------------------------------------------------------------------------

/// Frustum defined by 6 planes.
pub struct FrustumCulling {
    /// Planes as (normal, d): n·x + d = 0.
    pub planes: Vec<([f64; 3], f64)>,
}

impl FrustumCulling {
    /// Construct from 6 planes.
    pub fn new(planes: Vec<([f64; 3], f64)>) -> Self {
        FrustumCulling { planes }
    }

    /// Cull AABBs: return indices of visible (not fully outside) objects.
    pub fn cull(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| frustum_aabb_test(aabb, &self.planes))
            .map(|(i, _)| i)
            .collect()
    }

    /// Test a single AABB.
    pub fn test(&self, aabb: &QueryAabb) -> bool {
        frustum_aabb_test(aabb, &self.planes)
    }
}

// ---------------------------------------------------------------------------
// Sphere query
// ---------------------------------------------------------------------------

/// Query for shapes overlapping a sphere.
#[derive(Clone, Debug)]
pub struct SphereQuery {
    /// Center.
    pub center: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl SphereQuery {
    /// Construct a sphere query.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        SphereQuery { center, radius }
    }

    /// Return shape indices overlapping this sphere (against AABB list), sorted by distance.
    pub fn query_sorted(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let r2 = self.radius * self.radius;
        let mut hits: Vec<(usize, f64)> = aabbs
            .iter()
            .enumerate()
            .filter_map(|(i, aabb)| {
                let d2 = aabb.point_dist_sq(self.center);
                if d2 <= r2 { Some((i, d2)) } else { None }
            })
            .collect();
        hits.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        hits.into_iter().map(|(i, _)| i).collect()
    }

    /// Sphere AABB for broad-phase.
    pub fn aabb(&self) -> QueryAabb {
        let r = self.radius;
        QueryAabb {
            min: [self.center[0] - r, self.center[1] - r, self.center[2] - r],
            max: [self.center[0] + r, self.center[1] + r, self.center[2] + r],
        }
    }
}

// ---------------------------------------------------------------------------
// Capsule query
// ---------------------------------------------------------------------------

/// Query for shapes along a capsule sweep path.
pub struct CapsuleQuery {
    /// Start of capsule axis.
    pub p0: [f64; 3],
    /// End of capsule axis.
    pub p1: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl CapsuleQuery {
    /// Construct a capsule query.
    pub fn new(p0: [f64; 3], p1: [f64; 3], radius: f64) -> Self {
        CapsuleQuery { p0, p1, radius }
    }

    /// AABB enclosing the capsule.
    pub fn aabb(&self) -> QueryAabb {
        let r = self.radius;
        QueryAabb {
            min: [
                self.p0[0].min(self.p1[0]) - r,
                self.p0[1].min(self.p1[1]) - r,
                self.p0[2].min(self.p1[2]) - r,
            ],
            max: [
                self.p0[0].max(self.p1[0]) + r,
                self.p0[1].max(self.p1[1]) + r,
                self.p0[2].max(self.p1[2]) + r,
            ],
        }
    }

    /// Return shapes hit by this capsule (broad-phase AABB test).
    pub fn query(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let caps_aabb = self.aabb();
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| aabb.overlaps(&caps_aabb))
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Point query
// ---------------------------------------------------------------------------

/// Query for the N shapes nearest to a point.
pub struct PointQuery {
    /// Query point.
    pub point: [f64; 3],
}

impl PointQuery {
    /// Construct a point query.
    pub fn new(point: [f64; 3]) -> Self {
        PointQuery { point }
    }

    /// Return the `n` shape indices nearest to the point, sorted by distance.
    pub fn nearest_n(&self, aabbs: &[QueryAabb], n: usize) -> Vec<usize> {
        let mut dists: Vec<(usize, f64)> = aabbs
            .iter()
            .enumerate()
            .map(|(i, aabb)| (i, aabb.point_dist_sq(self.point).sqrt()))
            .collect();
        dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        dists.into_iter().take(n).map(|(i, _)| i).collect()
    }

    /// Exact distance to AABB i.
    pub fn distance_to(&self, aabb: &QueryAabb) -> f64 {
        aabb.point_dist_sq(self.point).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Convex sweep
// ---------------------------------------------------------------------------

/// Convex sweep: find all shapes hit by a box swept along `motion`.
pub struct ConvexSweep {
    /// AABB of the convex shape at rest.
    pub shape_aabb: QueryAabb,
    /// Sweep motion vector.
    pub motion: [f64; 3],
}

impl ConvexSweep {
    /// Construct a convex sweep.
    pub fn new(shape_aabb: QueryAabb, motion: [f64; 3]) -> Self {
        ConvexSweep { shape_aabb, motion }
    }

    /// Swept AABB enclosing the entire motion path.
    pub fn swept_aabb(&self) -> QueryAabb {
        sweep_aabb_motion(&self.shape_aabb, self.motion)
    }

    /// Return shapes hit by the swept volume.
    pub fn query(&self, aabbs: &[QueryAabb]) -> Vec<usize> {
        let swept = self.swept_aabb();
        aabbs
            .iter()
            .enumerate()
            .filter(|(_, aabb)| aabb.overlaps(&swept))
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Contact pair filter
// ---------------------------------------------------------------------------

/// Contact pair filter using collision layer masks.
pub struct ContactPairFilter {
    /// Layer mask for each shape.
    pub layer_masks: Vec<u32>,
    /// Group for each shape.
    pub groups: Vec<u32>,
}

impl ContactPairFilter {
    /// Construct a contact pair filter.
    pub fn new(layer_masks: Vec<u32>, groups: Vec<u32>) -> Self {
        ContactPairFilter {
            layer_masks,
            groups,
        }
    }

    /// True if shapes `a` and `b` should collide (bitmask AND ≠ 0 and same group check).
    pub fn should_collide(&self, a: usize, b: usize) -> bool {
        let mask_a = self.layer_masks.get(a).copied().unwrap_or(0xFFFF_FFFF);
        let mask_b = self.layer_masks.get(b).copied().unwrap_or(0xFFFF_FFFF);
        (mask_a & mask_b) != 0
    }

    /// Filter a list of shape pairs.
    pub fn filter(&self, pairs: &[ShapePair]) -> Vec<ShapePair> {
        pairs
            .iter()
            .filter(|p| self.should_collide(p.a, p.b))
            .cloned()
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Query pipeline
// ---------------------------------------------------------------------------

/// A pending query (for batched processing).
#[derive(Clone, Debug)]
pub enum PendingQuery {
    /// AABB overlap query.
    AabbQuery(QueryAabb),
    /// Sphere overlap query.
    SphereQuery(SphereQuery),
    /// Point nearest-N query.
    PointQuery {
        /// Query point in 3D space.
        point: [f64; 3],
        /// Number of nearest neighbours to find.
        n: usize,
    },
}

/// Result of a pending query.
#[derive(Clone, Debug)]
pub struct QueryResult {
    /// Shape indices returned by the query.
    pub hits: Vec<usize>,
}

/// Multi-query pipeline with deferred (next-frame) result dispatch.
pub struct QueryPipeline {
    /// Registered AABBs.
    pub aabbs: Vec<QueryAabb>,
    /// Pending queries.
    pub pending: Vec<PendingQuery>,
    /// Results from previous frame.
    pub results: Vec<QueryResult>,
}

impl QueryPipeline {
    /// Construct an empty pipeline.
    pub fn new() -> Self {
        QueryPipeline {
            aabbs: Vec::new(),
            pending: Vec::new(),
            results: Vec::new(),
        }
    }

    /// Register a shape AABB.
    pub fn add_shape(&mut self, aabb: QueryAabb) -> usize {
        self.aabbs.push(aabb);
        self.aabbs.len() - 1
    }

    /// Submit a query for deferred processing.
    pub fn submit(&mut self, query: PendingQuery) {
        self.pending.push(query);
    }

    /// Process all pending queries and store results (simulated next-frame).
    pub fn flush(&mut self) {
        let aabbs = &self.aabbs;
        let mut new_results = Vec::with_capacity(self.pending.len());
        for query in self.pending.drain(..) {
            let hits = match query {
                PendingQuery::AabbQuery(aabb) => aabbs
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| a.overlaps(&aabb))
                    .map(|(i, _)| i)
                    .collect(),
                PendingQuery::SphereQuery(sq) => {
                    let sphere_aabb = sq.aabb();
                    aabbs
                        .iter()
                        .enumerate()
                        .filter(|(_, a)| a.overlaps(&sphere_aabb))
                        .map(|(i, _)| i)
                        .collect()
                }
                PendingQuery::PointQuery { point, n } => {
                    let pq = PointQuery::new(point);
                    pq.nearest_n(aabbs, n)
                }
            };
            new_results.push(QueryResult { hits });
        }
        self.results = new_results;
    }

    /// Number of shapes registered.
    pub fn n_shapes(&self) -> usize {
        self.aabbs.len()
    }
}

impl Default for QueryPipeline {
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

    fn unit_aabb() -> QueryAabb {
        QueryAabb::new([0.0; 3], [1.0; 3])
    }

    fn shifted_aabb(offset: f64) -> QueryAabb {
        QueryAabb::new([offset; 3], [offset + 1.0; 3])
    }

    #[test]
    fn test_aabb_overlaps_touching() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([1.0; 3], [2.0; 3]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_union() {
        let a = QueryAabb::new([0.0; 3], [1.0; 3]);
        let b = QueryAabb::new([1.0; 3], [2.0; 3]);
        let u = a.union(&b);
        assert_eq!(u.min, [0.0; 3]);
        assert_eq!(u.max, [2.0; 3]);
    }

    #[test]
    fn test_aabb_contains_point() {
        let a = unit_aabb();
        assert!(a.contains_point([0.5, 0.5, 0.5]));
        assert!(!a.contains_point([2.0, 0.5, 0.5]));
    }

    #[test]
    fn test_aabb_point_dist_sq_inside() {
        let a = unit_aabb();
        assert_eq!(a.point_dist_sq([0.5, 0.5, 0.5]), 0.0);
    }

    #[test]
    fn test_aabb_point_dist_sq_outside() {
        let a = unit_aabb();
        let d2 = a.point_dist_sq([2.0, 0.5, 0.5]);
        assert!((d2 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_fatten() {
        let a = unit_aabb();
        let fat = a.fatten(0.5);
        assert!((fat.min[0] - (-0.5)).abs() < 1e-12);
        assert!((fat.max[0] - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_tree_insert_and_query() {
        let mut tree = AabbTree::new(0.1);
        let idx = tree.insert(unit_aabb(), 42);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.contains(&42), "hits={:?}, leaf_idx={}", hits, idx);
    }

    #[test]
    fn test_aabb_tree_remove() {
        let mut tree = AabbTree::new(0.1);
        let leaf = tree.insert(unit_aabb(), 0);
        tree.remove(leaf);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.is_empty());
    }

    #[test]
    fn test_aabb_tree_multiple_inserts() {
        let mut tree = AabbTree::new(0.1);
        tree.insert(unit_aabb(), 0);
        tree.insert(shifted_aabb(5.0), 1);
        let hits = tree.query_aabb(&unit_aabb());
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
    }

    #[test]
    fn test_aabb_tree_raycast() {
        let mut tree = AabbTree::new(0.1);
        tree.insert(unit_aabb(), 7);
        let hits = tree.raycast([-5.0, 0.5, 0.5], [1.0, 0.0, 0.0], 20.0);
        assert!(hits.contains(&7));
    }

    #[test]
    fn test_pair_query_brute_force() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5), shifted_aabb(10.0)];
        let pairs = PairQuery::brute_force(&aabbs);
        // 0 and 1 overlap, 2 does not
        assert!(
            pairs
                .iter()
                .any(|p| (p.a == 0 && p.b == 1) || (p.a == 1 && p.b == 0))
        );
        assert!(!pairs.iter().any(|p| p.a == 2 || p.b == 2));
    }

    #[test]
    fn test_pair_query_sweep_and_prune() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5)];
        let pairs = PairQuery::sweep_and_prune(&aabbs);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_pair_query_aabb_tree() {
        let aabbs = vec![unit_aabb(), shifted_aabb(0.5), shifted_aabb(10.0)];
        let pairs = PairQuery::aabb_tree(&aabbs);
        assert!(pairs.iter().any(|p| p.a == 0 && p.b == 1));
    }

    #[test]
    fn test_frustum_culling_inside() {
        // Frustum: half-space x > -10
        let planes = vec![([1.0_f64, 0.0, 0.0], 10.0)];
        let fc = FrustumCulling::new(planes);
        let visible = fc.cull(&[unit_aabb()]);
        assert!(visible.contains(&0));
    }

    #[test]
    fn test_frustum_culling_outside() {
        // Frustum: x > 100
        let planes = vec![([1.0_f64, 0.0, 0.0], -200.0)];
        let fc = FrustumCulling::new(planes);
        let visible = fc.cull(&[unit_aabb()]);
        assert!(visible.is_empty());
    }

    #[test]
    fn test_sphere_query_hit() {
        let sq = SphereQuery::new([0.5, 0.5, 0.5], 0.1);
        let hits = sq.query_sorted(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_sphere_query_miss() {
        let sq = SphereQuery::new([10.0, 10.0, 10.0], 0.1);
        let hits = sq.query_sorted(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_capsule_query_hit() {
        let cq = CapsuleQuery::new([0.5, 0.5, -1.0], [0.5, 0.5, 2.0], 0.1);
        let hits = cq.query(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_capsule_query_miss() {
        let cq = CapsuleQuery::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0], 0.1);
        let hits = cq.query(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_point_query_nearest() {
        let aabbs = vec![unit_aabb(), shifted_aabb(5.0)];
        let pq = PointQuery::new([0.5, 0.5, 0.5]);
        let nearest = pq.nearest_n(&aabbs, 1);
        assert_eq!(nearest[0], 0);
    }

    #[test]
    fn test_point_query_distance_inside() {
        let aabb = unit_aabb();
        let pq = PointQuery::new([0.5, 0.5, 0.5]);
        assert_eq!(pq.distance_to(&aabb), 0.0);
    }

    #[test]
    fn test_convex_sweep_hit() {
        let cs = ConvexSweep::new(QueryAabb::new([-1.0; 3], [0.0; 3]), [2.0, 2.0, 2.0]);
        let hits = cs.query(&[unit_aabb()]);
        assert!(hits.contains(&0));
    }

    #[test]
    fn test_convex_sweep_miss() {
        let cs = ConvexSweep::new(QueryAabb::new([-10.0; 3], [-9.0; 3]), [-1.0, 0.0, 0.0]);
        let hits = cs.query(&[unit_aabb()]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_contact_pair_filter_passes() {
        let cf = ContactPairFilter::new(vec![0xF, 0xF], vec![0, 0]);
        assert!(cf.should_collide(0, 1));
    }

    #[test]
    fn test_contact_pair_filter_blocks() {
        let cf = ContactPairFilter::new(vec![0x1, 0x2], vec![0, 0]);
        assert!(!cf.should_collide(0, 1));
    }

    #[test]
    fn test_contact_pair_filter_filter() {
        let cf = ContactPairFilter::new(vec![0xF, 0xF, 0x1], vec![0; 3]);
        let pairs = vec![ShapePair::new(0, 1), ShapePair::new(0, 2)];
        let filtered = cf.filter(&pairs);
        // 0↔1: mask 0xF&0xF ≠ 0 → passes; 0↔2: mask 0xF&0x1 = 0x1 ≠ 0 → passes
        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_query_pipeline_flush() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        qp.submit(PendingQuery::AabbQuery(unit_aabb()));
        qp.flush();
        assert_eq!(qp.results.len(), 1);
        assert!(qp.results[0].hits.contains(&0));
    }

    #[test]
    fn test_query_pipeline_sphere_query() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        let sq = SphereQuery::new([0.5, 0.5, 0.5], 0.5);
        qp.submit(PendingQuery::SphereQuery(sq));
        qp.flush();
        assert!(!qp.results.is_empty());
    }

    #[test]
    fn test_query_pipeline_point_query() {
        let mut qp = QueryPipeline::new();
        qp.add_shape(unit_aabb());
        qp.add_shape(shifted_aabb(5.0));
        qp.submit(PendingQuery::PointQuery {
            point: [0.5, 0.5, 0.5],
            n: 1,
        });
        qp.flush();
        assert_eq!(qp.results[0].hits.len(), 1);
        assert_eq!(qp.results[0].hits[0], 0);
    }

    #[test]
    fn test_sweep_aabb_motion_positive() {
        let aabb = unit_aabb();
        let swept = sweep_aabb_motion(&aabb, [2.0, 0.0, 0.0]);
        assert!((swept.max[0] - 3.0).abs() < 1e-12);
        assert!((swept.min[0] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_sweep_aabb_motion_negative() {
        let aabb = unit_aabb();
        let swept = sweep_aabb_motion(&aabb, [-1.0, 0.0, 0.0]);
        assert!((swept.min[0] - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_overlap_function() {
        let a = unit_aabb();
        let b = shifted_aabb(0.5);
        assert!(aabb_overlap(&a, &b));
    }

    // -----------------------------------------------------------------------
    // AabbTree: ancestor refit on insert, fat-AABB containment on update
    // -----------------------------------------------------------------------

    /// Deterministic xorshift generator for the randomized tests.
    struct TestRng(u64);

    impl TestRng {
        fn new(seed: u64) -> Self {
            TestRng(seed.max(1))
        }

        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        /// Uniform in `[0, 1)`.
        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }

        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.unit()
        }

        fn index(&mut self, len: usize) -> usize {
            (self.next_u64() % len as u64) as usize
        }
    }

    fn random_aabb(rng: &mut TestRng, extent: f64, max_size: f64) -> QueryAabb {
        let min = [
            rng.range(-extent, extent),
            rng.range(-extent, extent),
            rng.range(-extent, extent),
        ];
        let size = [
            rng.range(0.05, max_size),
            rng.range(0.05, max_size),
            rng.range(0.05, max_size),
        ];
        QueryAabb::new(min, [min[0] + size[0], min[1] + size[1], min[2] + size[2]])
    }

    fn box_at_x(x: f64) -> QueryAabb {
        QueryAabb::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0])
    }

    /// Ids of all live leaves whose fat box overlaps `query`, sorted. This reads the
    /// leaf-level data only and is the reference for what a query must return.
    fn leaf_level_hits(tree: &AabbTree, query: &QueryAabb) -> Vec<usize> {
        let mut hits: Vec<usize> = tree
            .nodes
            .iter()
            .filter(|n| n.is_leaf && n.height == 0 && n.fat_aabb.overlaps(query))
            .map(|n| n.user_data)
            .collect();
        hits.sort_unstable();
        hits
    }

    fn sorted(mut v: Vec<usize>) -> Vec<usize> {
        v.sort_unstable();
        v
    }

    #[test]
    fn test_aabb_contains_aabb() {
        let outer = QueryAabb::new([0.0; 3], [4.0; 3]);
        assert!(outer.contains_aabb(&QueryAabb::new([1.0; 3], [2.0; 3])));
        assert!(outer.contains_aabb(&outer));
        assert!(outer.contains_aabb(&QueryAabb::new([0.0; 3], [4.0, 1.0, 1.0])));
        // Partial overlap on one axis is not containment.
        assert!(!outer.contains_aabb(&QueryAabb::new([1.0, 1.0, 1.0], [5.0, 2.0, 2.0])));
        assert!(!outer.contains_aabb(&QueryAabb::new([-0.5, 1.0, 1.0], [2.0, 2.0, 2.0])));
        assert!(!outer.contains_aabb(&QueryAabb::new([10.0; 3], [11.0; 3])));
        // Containment is not symmetric.
        assert!(!QueryAabb::new([1.0; 3], [2.0; 3]).contains_aabb(&outer));
    }

    #[test]
    fn test_aabb_tree_new_clamps_negative_margin() {
        assert_eq!(AabbTree::new(-1.0).margin, 0.0);
        assert_eq!(AabbTree::new(f64::NAN).margin, 0.0);
        assert_eq!(AabbTree::new(0.25).margin, 0.25);
    }

    #[test]
    fn test_aabb_tree_insert_refits_ancestors_for_far_apart_boxes() {
        // A=[0,1], B=[10,11], C=[100,101]: before the fix C was linked below B's
        // parent without refitting the root, so a query for C was pruned at the root.
        let mut tree = AabbTree::new(0.01);
        let boxes = [box_at_x(0.0), box_at_x(10.0), box_at_x(100.0)];
        for (id, aabb) in boxes.iter().enumerate() {
            tree.insert(aabb.clone(), id);
            assert!(tree.validate(), "tree invalid after inserting box {id}");
        }
        for (id, aabb) in boxes.iter().enumerate() {
            assert_eq!(sorted(tree.query_aabb(aabb)), vec![id], "box {id}");
        }
        let root = tree.root.expect("tree has a root");
        let root_box = &tree.nodes[root].fat_aabb;
        assert!(root_box.contains_aabb(&boxes[2].fatten(0.01)));
        assert!(root_box.contains_aabb(&boxes[0].fatten(0.01)));
    }

    #[test]
    fn test_aabb_tree_two_hundred_random_boxes_are_all_reachable() {
        let mut rng = TestRng::new(0x5EED_1234);
        let mut tree = AabbTree::new(0.05);
        let mut boxes = Vec::new();
        for id in 0..200 {
            let aabb = random_aabb(&mut rng, 50.0, 4.0);
            tree.insert(aabb.clone(), id);
            boxes.push(aabb);
            assert!(tree.validate(), "tree invalid after insert {id}");
        }
        for (id, aabb) in boxes.iter().enumerate() {
            let hits = tree.query_aabb(aabb);
            assert!(hits.contains(&id), "box {id} unreachable, hits={hits:?}");
        }
        // The hierarchy must return exactly what a scan of the leaves returns.
        for _ in 0..100 {
            let query = random_aabb(&mut rng, 55.0, 12.0);
            assert_eq!(
                sorted(tree.query_aabb(&query)),
                leaf_level_hits(&tree, &query)
            );
        }
    }

    #[test]
    fn test_aabb_tree_sorted_insertion_stays_balanced() {
        // Monotone insertion is the worst case for an unbalanced tree.
        let mut tree = AabbTree::new(0.1);
        for id in 0..1024 {
            tree.insert(box_at_x(id as f64 * 2.0), id);
        }
        assert!(tree.validate());
        assert!(tree.height() <= 16, "height {}", tree.height());
        assert_eq!(sorted(tree.query_aabb(&box_at_x(2.0 * 700.0))), vec![700]);
    }

    #[test]
    fn test_aabb_tree_update_inside_fat_aabb_leaves_tree_untouched() {
        let mut tree = AabbTree::new(0.1);
        let leaf = tree.insert(unit_aabb(), 5);
        tree.insert(shifted_aabb(5.0), 6);
        let node_count = tree.n_nodes();
        let fat_before = tree.nodes[leaf].fat_aabb.clone();
        let moved = QueryAabb::new([0.05; 3], [1.05; 3]);
        tree.update(leaf, moved.clone());
        assert_eq!(tree.n_nodes(), node_count);
        assert_eq!(tree.nodes[leaf].fat_aabb, fat_before);
        assert_eq!(tree.nodes[leaf].aabb, moved);
        assert!(tree.validate());
    }

    #[test]
    fn test_aabb_tree_update_partially_outside_fat_aabb_reinserts() {
        // The unit box has a fat AABB of [-0.1, 1.1]. Moving it to x in [0.6, 1.6]
        // still overlaps that fat box, but is not contained in it: the old
        // overlap test kept the stale fat AABB and hid the part beyond x = 1.1.
        let mut tree = AabbTree::new(0.1);
        let leaf = tree.insert(unit_aabb(), 3);
        tree.update(leaf, QueryAabb::new([0.6, 0.0, 0.0], [1.6, 1.0, 1.0]));
        let probe = QueryAabb::new([1.5, 0.0, 0.0], [1.55, 1.0, 1.0]);
        assert_eq!(tree.query_aabb(&probe), vec![3]);
        assert!(
            tree.nodes[leaf]
                .fat_aabb
                .contains_aabb(&tree.nodes[leaf].aabb)
        );
        assert!(tree.validate());
        // And the ray test sees the protruding part as well.
        assert_eq!(
            tree.raycast([1.55, 0.5, 5.0], [0.0, 0.0, -1.0], 10.0),
            vec![3]
        );
    }

    #[test]
    fn test_aabb_tree_update_keeps_leaf_handles_valid() {
        let mut tree = AabbTree::new(0.05);
        let handles: Vec<usize> = (0..20)
            .map(|i| tree.insert(box_at_x(i as f64 * 3.0), i))
            .collect();
        // Move every box far away: each update has to re-insert its leaf.
        for (i, &handle) in handles.iter().enumerate() {
            let new_box = box_at_x(1000.0 + i as f64 * 3.0);
            tree.update(handle, new_box.clone());
            assert!(tree.validate(), "invalid after updating {i}");
            assert!(tree.nodes[handle].is_leaf);
            assert_eq!(tree.nodes[handle].user_data, i);
            assert_eq!(tree.nodes[handle].aabb, new_box);
        }
        for (i, &handle) in handles.iter().enumerate() {
            let hits = tree.query_aabb(&box_at_x(1000.0 + i as f64 * 3.0));
            assert_eq!(hits, vec![i]);
            assert!(tree.query_aabb(&box_at_x(i as f64 * 3.0)).is_empty());
            // The handle obtained from `insert` still addresses the same leaf.
            tree.remove(handle);
            assert!(tree.validate(), "invalid after removing {i}");
        }
        assert!(tree.root.is_none());
        assert_eq!(tree.n_nodes(), 0);
    }

    #[test]
    fn test_aabb_tree_random_insert_remove_update_mix_matches_leaf_scan() {
        let mut rng = TestRng::new(0xC0FF_EE11);
        let mut tree = AabbTree::new(0.1);
        // (handle, tight box) of every live object; the id is the handle.
        let mut live: Vec<(usize, QueryAabb)> = Vec::new();
        for step in 0..600 {
            let choice = rng.unit();
            if live.is_empty() || choice < 0.45 {
                let aabb = random_aabb(&mut rng, 30.0, 3.0);
                let handle = tree.insert(aabb.clone(), usize::MAX - step);
                tree.nodes[handle].user_data = handle;
                live.push((handle, aabb));
            } else if choice < 0.65 {
                let k = rng.index(live.len());
                let (handle, _) = live.swap_remove(k);
                tree.remove(handle);
            } else {
                let k = rng.index(live.len());
                let (handle, old) = live[k].clone();
                let new_box = if rng.unit() < 0.5 {
                    // A small nudge that usually stays inside the fat AABB.
                    let d = rng.range(-0.05, 0.05);
                    QueryAabb::new(
                        [old.min[0] + d, old.min[1], old.min[2]],
                        [old.max[0] + d, old.max[1], old.max[2]],
                    )
                } else {
                    random_aabb(&mut rng, 30.0, 3.0)
                };
                tree.update(handle, new_box.clone());
                live[k].1 = new_box;
            }
            assert!(tree.validate(), "tree invalid at step {step}");
            assert_eq!(
                tree.n_nodes(),
                if live.is_empty() {
                    0
                } else {
                    2 * live.len() - 1
                }
            );
            if step % 20 == 0 {
                for (handle, tight) in &live {
                    assert!(
                        tree.query_aabb(tight).contains(handle),
                        "object {handle} lost at step {step}"
                    );
                }
                let query = random_aabb(&mut rng, 32.0, 10.0);
                assert_eq!(
                    sorted(tree.query_aabb(&query)),
                    leaf_level_hits(&tree, &query),
                    "query mismatch at step {step}"
                );
            }
        }
    }

    #[test]
    fn test_aabb_tree_remove_and_update_ignore_invalid_indices() {
        let mut tree = AabbTree::new(0.1);
        let a = tree.insert(unit_aabb(), 0);
        let b = tree.insert(shifted_aabb(5.0), 1);
        let c = tree.insert(shifted_aabb(10.0), 2);
        let root = tree.root.expect("tree has a root");
        assert!(!tree.nodes[root].is_leaf);
        // Internal node, out-of-range index.
        tree.remove(root);
        tree.remove(10_000);
        tree.update(root, unit_aabb());
        tree.update(10_000, unit_aabb());
        assert!(tree.validate());
        assert_eq!(tree.n_nodes(), 5);
        // Double removal must not corrupt the free list.
        tree.remove(b);
        tree.remove(b);
        tree.update(b, unit_aabb());
        assert!(tree.validate());
        assert_eq!(tree.n_nodes(), 3);
        let d = tree.insert(shifted_aabb(20.0), 3);
        let e = tree.insert(shifted_aabb(30.0), 4);
        assert!(tree.validate());
        let mut leaves = vec![a, c, d, e];
        leaves.sort_unstable();
        leaves.dedup();
        assert_eq!(leaves.len(), 4, "free list handed out a slot twice");
        assert_eq!(
            sorted(tree.query_aabb(&QueryAabb::new([-1.0; 3], [100.0; 3]))),
            vec![0, 2, 3, 4]
        );
    }

    #[test]
    fn test_aabb_tree_validate_rejects_corrupted_trees() {
        fn three_leaf_tree() -> AabbTree {
            let mut tree = AabbTree::new(0.1);
            tree.insert(box_at_x(0.0), 0);
            tree.insert(box_at_x(10.0), 1);
            tree.insert(box_at_x(100.0), 2);
            tree
        }
        assert!(three_leaf_tree().validate());

        // A root box that no longer covers a child.
        let mut tree = three_leaf_tree();
        let root = tree.root.expect("root");
        tree.nodes[root].fat_aabb = box_at_x(0.0);
        assert!(!tree.validate());

        // A stale height.
        let mut tree = three_leaf_tree();
        let root = tree.root.expect("root");
        tree.nodes[root].height += 1;
        assert!(!tree.validate());

        // A broken parent link.
        let mut tree = three_leaf_tree();
        let root = tree.root.expect("root");
        let (left, _) = tree.children_of(root).expect("children");
        tree.nodes[left].parent = None;
        assert!(!tree.validate());

        // A leaf whose tight box left its fat box.
        let mut tree = three_leaf_tree();
        let leaf = (0..tree.nodes.len())
            .find(|&i| tree.nodes[i].is_leaf)
            .expect("leaf");
        tree.nodes[leaf].aabb = box_at_x(500.0);
        assert!(!tree.validate());

        // A leaked slot (neither reachable nor free).
        let mut tree = three_leaf_tree();
        tree.nodes.push(DynamicNode::leaf(unit_aabb(), 0.1, 9));
        assert!(!tree.validate());
    }

    #[test]
    fn test_rebalance_tree_rebuilds_a_balanced_tree_and_keeps_handles() {
        let mut tree = AabbTree::new(0.05);
        let handles: Vec<usize> = (0..64)
            .map(|i| tree.insert(box_at_x(i as f64 * 5.0), i))
            .collect();
        // Remove every other box and re-add some elsewhere to churn the tree.
        for &handle in handles.iter().step_by(2) {
            tree.remove(handle);
        }
        let extra: Vec<usize> = (100..110)
            .map(|i| tree.insert(box_at_x(i as f64 * 5.0), i))
            .collect();
        assert!(tree.validate());
        rebalance_tree(&mut tree);
        assert!(tree.validate());
        // 32 + 10 = 42 leaves: a median-split tree has height ceil(log2(42)) = 6.
        assert_eq!(tree.height(), 6);
        assert_eq!(tree.n_nodes(), 2 * 42 - 1);
        for (i, &handle) in handles.iter().enumerate().filter(|(i, _)| i % 2 == 1) {
            assert_eq!(tree.nodes[handle].user_data, i);
            assert_eq!(tree.query_aabb(&box_at_x(i as f64 * 5.0)), vec![i]);
        }
        for (k, &handle) in extra.iter().enumerate() {
            assert_eq!(tree.nodes[handle].user_data, 100 + k);
        }
        // The rebuilt tree stays fully usable.
        let id = tree.insert(box_at_x(-40.0), 999);
        assert!(tree.validate());
        assert_eq!(tree.query_aabb(&box_at_x(-40.0)), vec![999]);
        tree.remove(id);
        assert!(tree.validate());
        // Rebuilding an empty tree and a single-leaf tree is a no-op.
        let mut empty = AabbTree::new(0.1);
        rebalance_tree(&mut empty);
        assert!(empty.validate() && empty.root.is_none());
        let mut single = AabbTree::new(0.1);
        let only = single.insert(unit_aabb(), 8);
        rebalance_tree(&mut single);
        assert_eq!(single.root, Some(only));
        assert!(single.validate());
    }

    #[test]
    fn test_pair_query_aabb_tree_matches_brute_force_on_random_boxes() {
        let mut rng = TestRng::new(0xABCD_0042);
        let boxes: Vec<QueryAabb> = (0..200).map(|_| random_aabb(&mut rng, 8.0, 3.0)).collect();
        let brute = PairQuery::brute_force(&boxes);
        assert!(brute.len() > 50, "scene too sparse: {} pairs", brute.len());
        assert_eq!(PairQuery::aabb_tree(&boxes), brute);
        let mut swept = PairQuery::sweep_and_prune(&boxes);
        swept.sort_unstable_by_key(|p| (p.a, p.b));
        assert_eq!(swept, brute);
    }

    #[test]
    fn test_pair_query_aabb_tree_ignores_fat_margin_only_overlaps() {
        // 0.005 apart: the tree's fat boxes (margin 0.01) overlap, the boxes do not.
        let boxes = vec![
            QueryAabb::new([0.0; 3], [1.0; 3]),
            QueryAabb::new([1.005, 0.0, 0.0], [2.0, 1.0, 1.0]),
        ];
        assert!(PairQuery::brute_force(&boxes).is_empty());
        assert!(PairQuery::aabb_tree(&boxes).is_empty());
    }

    #[test]
    fn test_aabb_tree_raycast_reaches_every_leaf_after_many_inserts() {
        let mut tree = AabbTree::new(0.1);
        for id in 0..50 {
            tree.insert(box_at_x(id as f64 * 3.0), id);
        }
        let along_x = tree.raycast([-5.0, 0.5, 0.5], [1.0, 0.0, 0.0], 1000.0);
        assert_eq!(sorted(along_x), (0..50).collect::<Vec<_>>());
        let along_y = tree.raycast([7.0 * 3.0 + 0.5, -5.0, 0.5], [0.0, 1.0, 0.0], 100.0);
        assert_eq!(along_y, vec![7]);
        let miss = tree.raycast([-5.0, 5.0, 0.5], [1.0, 0.0, 0.0], 1000.0);
        assert!(miss.is_empty());
    }
}
