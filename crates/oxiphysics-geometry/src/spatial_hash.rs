#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial hashing, loose octree, k-nearest neighbours, ray grid traversal,
//! and spatial statistics.

#![allow(dead_code)]

use std::collections::{BinaryHeap, HashMap};

// ─── Helper math ─────────────────────────────────────────────────────────────

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
}

// ─── SpatialHash3D ────────────────────────────────────────────────────────────

/// A uniform-grid spatial hash that maps 3-D positions to arbitrary data.
///
/// Items are bucketed by their cell coordinate
/// `(floor(x/cell_size), floor(y/cell_size), floor(z/cell_size))`.
pub struct SpatialHash3D<T: Clone> {
    /// Side length of each cubic cell.
    pub cell_size: f64,
    cells: HashMap<(i32, i32, i32), Vec<(usize, T)>>,
    /// Total number of items currently stored.
    pub n_items: usize,
}

impl<T: Clone> SpatialHash3D<T> {
    /// Creates a new empty spatial hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
            n_items: 0,
        }
    }

    /// Maps a single coordinate to its cell index via `floor(x / cell_size)`.
    pub fn cell_coord(&self, x: f64) -> i32 {
        (x / self.cell_size).floor() as i32
    }

    /// Returns the 3-D cell key for a world-space position.
    pub fn cell_key(&self, pos: [f64; 3]) -> (i32, i32, i32) {
        (
            self.cell_coord(pos[0]),
            self.cell_coord(pos[1]),
            self.cell_coord(pos[2]),
        )
    }

    /// Inserts `(id, data)` at `pos` and returns the cell key it was placed in.
    pub fn insert(&mut self, id: usize, pos: [f64; 3], data: T) -> (i32, i32, i32) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push((id, data));
        self.n_items += 1;
        key
    }

    /// Removes the entry with the given `id` from the cell that contains `pos`.
    pub fn remove(&mut self, id: usize, pos: [f64; 3]) {
        let key = self.cell_key(pos);
        if let Some(bucket) = self.cells.get_mut(&key) {
            let before = bucket.len();
            bucket.retain(|(item_id, _)| *item_id != id);
            let removed = before - bucket.len();
            self.n_items = self.n_items.saturating_sub(removed);
            if bucket.is_empty() {
                self.cells.remove(&key);
            }
        }
    }

    /// Returns all items within the cells overlapping a sphere of `radius`
    /// centred at `center`. Note: without stored positions, exact distance
    /// filtering is not possible; use `SpatialHashPos3D` for exact queries.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<(usize, &T)> {
        let r_cells = (radius / self.cell_size).ceil() as i32;
        let cx = self.cell_coord(center[0]);
        let cy = self.cell_coord(center[1]);
        let cz = self.cell_coord(center[2]);

        let mut result = Vec::new();
        for dx in -r_cells..=r_cells {
            for dy in -r_cells..=r_cells {
                for dz in -r_cells..=r_cells {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(bucket) = self.cells.get(&key) {
                        for (id, data) in bucket {
                            result.push((*id, data));
                        }
                    }
                }
            }
        }
        result
    }

    /// Returns all items stored in the given cell key.
    pub fn query_cell(&self, key: (i32, i32, i32)) -> Vec<(usize, &T)> {
        match self.cells.get(&key) {
            Some(bucket) => bucket.iter().map(|(id, data)| (*id, data)).collect(),
            None => Vec::new(),
        }
    }

    /// Removes all items.
    pub fn clear(&mut self) {
        self.cells.clear();
        self.n_items = 0;
    }

    /// Returns the total number of items stored.
    pub fn len(&self) -> usize {
        self.n_items
    }

    /// Returns `true` if no items are stored.
    pub fn is_empty(&self) -> bool {
        self.n_items == 0
    }

    /// Returns the number of non-empty cells.
    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }
}

// ─── SpatialHashPos3D ─────────────────────────────────────────────────────────

/// A uniform-grid spatial hash with position storage for exact radius queries.
pub struct SpatialHashPos3D<T: Clone> {
    /// Side length of each cubic cell.
    pub cell_size: f64,
    cells: HashMap<(i32, i32, i32), Vec<(usize, [f64; 3], T)>>,
    /// Total number of items currently stored.
    pub n_items: usize,
}

impl<T: Clone> SpatialHashPos3D<T> {
    /// Creates a new empty hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
            n_items: 0,
        }
    }

    /// Maps a coordinate to its cell index.
    pub fn cell_coord(&self, x: f64) -> i32 {
        (x / self.cell_size).floor() as i32
    }

    /// Returns the 3-D cell key for a world-space position.
    pub fn cell_key(&self, pos: [f64; 3]) -> (i32, i32, i32) {
        (
            self.cell_coord(pos[0]),
            self.cell_coord(pos[1]),
            self.cell_coord(pos[2]),
        )
    }

    /// Inserts item and returns its cell key.
    pub fn insert(&mut self, id: usize, pos: [f64; 3], data: T) -> (i32, i32, i32) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push((id, pos, data));
        self.n_items += 1;
        key
    }

    /// Removes the entry with the given id from the cell containing pos.
    pub fn remove(&mut self, id: usize, pos: [f64; 3]) {
        let key = self.cell_key(pos);
        if let Some(bucket) = self.cells.get_mut(&key) {
            let before = bucket.len();
            bucket.retain(|(item_id, _, _)| *item_id != id);
            let removed = before - bucket.len();
            self.n_items = self.n_items.saturating_sub(removed);
            if bucket.is_empty() {
                self.cells.remove(&key);
            }
        }
    }

    /// Returns all items within a sphere, with exact distance filtering.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<(usize, &T)> {
        let r_cells = (radius / self.cell_size).ceil() as i32;
        let cx = self.cell_coord(center[0]);
        let cy = self.cell_coord(center[1]);
        let cz = self.cell_coord(center[2]);
        let r2 = radius * radius;
        let mut result = Vec::new();

        for dx in -r_cells..=r_cells {
            for dy in -r_cells..=r_cells {
                for dz in -r_cells..=r_cells {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(bucket) = self.cells.get(&key) {
                        for (id, pos, data) in bucket {
                            let d2 = dist_sq(*pos, center);
                            if d2 <= r2 {
                                result.push((*id, data));
                            }
                        }
                    }
                }
            }
        }
        result
    }

    /// Returns all items within an axis-aligned box defined by `min` and `max`.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<(usize, &T)> {
        let cmin = [
            self.cell_coord(min[0]),
            self.cell_coord(min[1]),
            self.cell_coord(min[2]),
        ];
        let cmax = [
            self.cell_coord(max[0]),
            self.cell_coord(max[1]),
            self.cell_coord(max[2]),
        ];
        let mut result = Vec::new();
        for cx in cmin[0]..=cmax[0] {
            for cy in cmin[1]..=cmax[1] {
                for cz in cmin[2]..=cmax[2] {
                    if let Some(bucket) = self.cells.get(&(cx, cy, cz)) {
                        for (id, pos, data) in bucket {
                            if pos[0] >= min[0]
                                && pos[0] <= max[0]
                                && pos[1] >= min[1]
                                && pos[1] <= max[1]
                                && pos[2] >= min[2]
                                && pos[2] <= max[2]
                            {
                                result.push((*id, data));
                            }
                        }
                    }
                }
            }
        }
        result
    }

    /// Find the k nearest neighbours of `center`.
    ///
    /// Uses an expanding search radius starting at cell_size and doubling
    /// until at least k candidates are found, then filters to the true k
    /// nearest.
    pub fn k_nearest(&self, center: [f64; 3], k: usize) -> Vec<(usize, f64, &T)> {
        if k == 0 || self.n_items == 0 {
            return Vec::new();
        }

        // Collect all items with their distances
        let mut all: Vec<(usize, f64, &T)> = Vec::new();
        for bucket in self.cells.values() {
            for (id, pos, data) in bucket {
                let d2 = dist_sq(*pos, center);
                all.push((*id, d2, data));
            }
        }
        all.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        all.truncate(k);
        // Convert d2 to distance
        all.iter()
            .map(|(id, d2, data)| (*id, d2.sqrt(), *data))
            .collect()
    }

    /// Returns all items stored in the given cell.
    pub fn query_cell(&self, key: (i32, i32, i32)) -> Vec<(usize, &T)> {
        match self.cells.get(&key) {
            Some(bucket) => bucket.iter().map(|(id, _, data)| (*id, data)).collect(),
            None => Vec::new(),
        }
    }

    /// Removes all items.
    pub fn clear(&mut self) {
        self.cells.clear();
        self.n_items = 0;
    }

    /// Returns the total number of items stored.
    pub fn len(&self) -> usize {
        self.n_items
    }

    /// Returns `true` if no items are stored.
    pub fn is_empty(&self) -> bool {
        self.n_items == 0
    }

    /// Returns the number of non-empty cells.
    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }

    /// Compute statistics about the spatial hash distribution.
    pub fn statistics(&self) -> SpatialHashStats {
        let n_cells = self.cells.len();
        if n_cells == 0 {
            return SpatialHashStats {
                n_items: 0,
                n_occupied_cells: 0,
                min_bucket_size: 0,
                max_bucket_size: 0,
                avg_bucket_size: 0.0,
                load_factor: 0.0,
            };
        }
        let sizes: Vec<usize> = self.cells.values().map(|b| b.len()).collect();
        let min_b = *sizes.iter().min().unwrap_or(&0);
        let max_b = *sizes.iter().max().unwrap_or(&0);
        let sum: usize = sizes.iter().sum();
        SpatialHashStats {
            n_items: self.n_items,
            n_occupied_cells: n_cells,
            min_bucket_size: min_b,
            max_bucket_size: max_b,
            avg_bucket_size: sum as f64 / n_cells as f64,
            load_factor: sum as f64 / n_cells as f64,
        }
    }

    /// Update the position of an item (remove from old cell, insert into new).
    pub fn update_position(&mut self, id: usize, old_pos: [f64; 3], new_pos: [f64; 3], data: T) {
        self.remove(id, old_pos);
        self.insert(id, new_pos, data);
    }
}

/// Statistics about a spatial hash's distribution.
#[derive(Debug, Clone)]
pub struct SpatialHashStats {
    /// Total items stored.
    pub n_items: usize,
    /// Number of occupied cells.
    pub n_occupied_cells: usize,
    /// Smallest bucket.
    pub min_bucket_size: usize,
    /// Largest bucket.
    pub max_bucket_size: usize,
    /// Average items per occupied cell.
    pub avg_bucket_size: f64,
    /// Load factor (same as avg_bucket_size).
    pub load_factor: f64,
}

// ─── Ray traversal through uniform grid (3D DDA) ────────────────────────────

/// Result of a ray traversal step through the grid.
#[derive(Debug, Clone, Copy)]
pub struct RayGridStep {
    /// The cell key visited.
    pub cell: (i32, i32, i32),
    /// Parameter t at entry to this cell.
    pub t_enter: f64,
    /// Parameter t at exit from this cell.
    pub t_exit: f64,
}

/// Traverse a ray through a uniform grid using 3D DDA (Amanatides-Woo).
///
/// Returns the sequence of cells visited by the ray from `origin` in
/// `direction`, up to `max_t`.
pub fn ray_traverse_grid(
    origin: [f64; 3],
    direction: [f64; 3],
    cell_size: f64,
    max_t: f64,
    max_steps: usize,
) -> Vec<RayGridStep> {
    let mut result = Vec::new();
    if cell_size <= 0.0 {
        return result;
    }

    // Normalize direction
    let dir_len =
        (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
            .sqrt();
    if dir_len < 1e-15 {
        return result;
    }
    let dir = [
        direction[0] / dir_len,
        direction[1] / dir_len,
        direction[2] / dir_len,
    ];

    // Current cell
    let mut cell = [
        (origin[0] / cell_size).floor() as i32,
        (origin[1] / cell_size).floor() as i32,
        (origin[2] / cell_size).floor() as i32,
    ];

    // Step direction (+1 or -1) and t_delta (how much t to cross one cell)
    let mut step = [0i32; 3];
    let mut t_max = [0.0f64; 3]; // t at which ray crosses next cell boundary
    let mut t_delta = [f64::INFINITY; 3];

    for i in 0..3 {
        if dir[i] > 0.0 {
            step[i] = 1;
            let next_boundary = (cell[i] as f64 + 1.0) * cell_size;
            t_max[i] = (next_boundary - origin[i]) / dir[i];
            t_delta[i] = cell_size / dir[i];
        } else if dir[i] < 0.0 {
            step[i] = -1;
            let next_boundary = cell[i] as f64 * cell_size;
            t_max[i] = (next_boundary - origin[i]) / dir[i];
            t_delta[i] = cell_size / (-dir[i]);
        } else {
            step[i] = 0;
            t_max[i] = f64::INFINITY;
            t_delta[i] = f64::INFINITY;
        }
    }

    let mut t = 0.0;
    for _ in 0..max_steps {
        // Find the axis with smallest t_max
        let min_t = t_max[0].min(t_max[1]).min(t_max[2]);
        let t_exit = min_t.min(max_t);

        result.push(RayGridStep {
            cell: (cell[0], cell[1], cell[2]),
            t_enter: t,
            t_exit,
        });

        if min_t >= max_t {
            break;
        }

        // Advance along the axis with smallest t_max
        if t_max[0] <= t_max[1] && t_max[0] <= t_max[2] {
            cell[0] += step[0];
            t = t_max[0];
            t_max[0] += t_delta[0];
        } else if t_max[1] <= t_max[2] {
            cell[1] += step[1];
            t = t_max[1];
            t_max[1] += t_delta[1];
        } else {
            cell[2] += step[2];
            t = t_max[2];
            t_max[2] += t_delta[2];
        }
    }

    result
}

// ─── LooseOctree ──────────────────────────────────────────────────────────────

/// A loose octree node that stores items at every level.
pub struct LooseOctree<T: Clone> {
    /// World-space centre of this node.
    pub center: [f64; 3],
    /// Half-size of the *tight* cell (the loose cell is `2 * half_size`).
    pub half_size: f64,
    /// Remaining subdivision depth (0 = leaf).
    pub depth: u32,
    /// Items stored at this node (not yet pushed to children).
    pub items: Vec<(usize, [f64; 3], T)>,
    /// Eight children, present after the first subdivision.
    pub children: Option<Box<[LooseOctree<T>; 8]>>,
}

impl<T: Clone> LooseOctree<T> {
    /// Creates a new root node.
    pub fn new(center: [f64; 3], half_size: f64, max_depth: u32) -> Self {
        Self {
            center,
            half_size,
            depth: max_depth,
            items: Vec::new(),
            children: None,
        }
    }

    /// Returns the octant index (0..7) for `pos` relative to `self.center`.
    pub fn child_index(&self, pos: [f64; 3]) -> usize {
        let mut idx = 0usize;
        if pos[0] >= self.center[0] {
            idx |= 1;
        }
        if pos[1] >= self.center[1] {
            idx |= 2;
        }
        if pos[2] >= self.center[2] {
            idx |= 4;
        }
        idx
    }

    /// Builds child centres for octant `i`.
    fn child_center(&self, i: usize) -> [f64; 3] {
        let q = self.half_size * 0.5;
        [
            self.center[0] + if i & 1 != 0 { q } else { -q },
            self.center[1] + if i & 2 != 0 { q } else { -q },
            self.center[2] + if i & 4 != 0 { q } else { -q },
        ]
    }

    /// Subdivides this node, creating eight children.
    fn subdivide(&mut self) {
        if self.children.is_some() {
            return;
        }
        let child_depth = self.depth.saturating_sub(1);
        let child_half = self.half_size * 0.5;
        let children = std::array::from_fn(|i| {
            LooseOctree::new(self.child_center(i), child_half, child_depth)
        });
        self.children = Some(Box::new(children));
    }

    /// Inserts `(id, pos, data)` into the tree.
    pub fn insert(&mut self, id: usize, pos: [f64; 3], data: T, max_items: usize) {
        self.items.push((id, pos, data.clone()));

        if self.items.len() > max_items && self.depth > 0 {
            self.subdivide();
            let items: Vec<_> = self.items.drain(..).collect();
            for (iid, ipos, idata) in items {
                let ci = self.child_index(ipos);
                if let Some(children) = self.children.as_mut() {
                    children[ci].insert(iid, ipos, idata, max_items);
                }
            }
        }
    }

    /// Returns all items whose stored position is within `radius` of `center`.
    pub fn query_sphere(&self, center: [f64; 3], radius: f64) -> Vec<(usize, &T)> {
        let loose = self.half_size * 2.0;
        let mut closest2 = 0.0f64;
        for i in 0..3 {
            let lo = self.center[i] - loose;
            let hi = self.center[i] + loose;
            if center[i] < lo {
                closest2 += (center[i] - lo).powi(2);
            } else if center[i] > hi {
                closest2 += (center[i] - hi).powi(2);
            }
        }
        if closest2 > radius * radius {
            return Vec::new();
        }

        let r2 = radius * radius;
        let mut result = Vec::new();

        for (id, pos, data) in &self.items {
            let d2 = dist_sq(*pos, center);
            if d2 <= r2 {
                result.push((*id, data));
            }
        }

        if let Some(children) = &self.children {
            for child in children.iter() {
                result.extend(child.query_sphere(center, radius));
            }
        }

        result
    }

    /// Returns all items within an axis-aligned box.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<(usize, &T)> {
        // Quick rejection: does loose node box overlap query box?
        let loose = self.half_size * 2.0;
        for i in 0..3 {
            if self.center[i] - loose > max[i] || self.center[i] + loose < min[i] {
                return Vec::new();
            }
        }

        let mut result = Vec::new();
        for (id, pos, data) in &self.items {
            if pos[0] >= min[0]
                && pos[0] <= max[0]
                && pos[1] >= min[1]
                && pos[1] <= max[1]
                && pos[2] >= min[2]
                && pos[2] <= max[2]
            {
                result.push((*id, data));
            }
        }

        if let Some(children) = &self.children {
            for child in children.iter() {
                result.extend(child.query_aabb(min, max));
            }
        }

        result
    }

    /// Counts all items stored in this node and all descendants.
    pub fn item_count(&self) -> usize {
        let mut count = self.items.len();
        if let Some(children) = &self.children {
            for child in children.iter() {
                count += child.item_count();
            }
        }
        count
    }

    /// Returns the maximum depth reached in the tree.
    pub fn max_depth_reached(&self) -> u32 {
        if let Some(children) = &self.children {
            let child_max = children
                .iter()
                .map(|c| c.max_depth_reached())
                .max()
                .unwrap_or(0);
            child_max + 1
        } else {
            0
        }
    }

    /// Count the total number of nodes in the tree.
    pub fn node_count(&self) -> usize {
        let mut count = 1;
        if let Some(children) = &self.children {
            for child in children.iter() {
                count += child.node_count();
            }
        }
        count
    }

    /// K-nearest neighbours search via the octree.
    pub fn k_nearest(&self, center: [f64; 3], k: usize) -> Vec<(usize, f64, &T)> {
        // Use a max-heap to keep track of the k closest
        let mut heap: BinaryHeap<KnnEntry<&T>> = BinaryHeap::new();
        self.knn_recurse(center, k, &mut heap);

        let mut result: Vec<(usize, f64, &T)> =
            heap.into_iter().map(|e| (e.id, e.dist, e.data)).collect();
        result.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        result
    }

    fn knn_recurse<'a>(
        &'a self,
        center: [f64; 3],
        k: usize,
        heap: &mut BinaryHeap<KnnEntry<&'a T>>,
    ) {
        // Check items at this node
        for (id, pos, data) in &self.items {
            let d = dist_sq(*pos, center).sqrt();
            if heap.len() < k {
                heap.push(KnnEntry {
                    dist: d,
                    id: *id,
                    data,
                });
            } else if let Some(top) = heap.peek()
                && d < top.dist
            {
                heap.pop();
                heap.push(KnnEntry {
                    dist: d,
                    id: *id,
                    data,
                });
            }
        }

        // Recurse into children that could contain closer points
        if let Some(children) = &self.children {
            for child in children.iter() {
                // Minimum distance from center to child's loose box
                let loose = child.half_size * 2.0;
                let mut min_d2 = 0.0f64;
                for i in 0..3 {
                    let lo = child.center[i] - loose;
                    let hi = child.center[i] + loose;
                    if center[i] < lo {
                        min_d2 += (center[i] - lo).powi(2);
                    } else if center[i] > hi {
                        min_d2 += (center[i] - hi).powi(2);
                    }
                }
                let min_d = min_d2.sqrt();
                let should_search =
                    heap.len() < k || heap.peek().is_some_and(|top| min_d < top.dist);
                if should_search {
                    child.knn_recurse(center, k, heap);
                }
            }
        }
    }
}

/// Internal entry for KNN max-heap (ordered by distance descending).
struct KnnEntry<D> {
    dist: f64,
    id: usize,
    data: D,
}

impl<D> PartialEq for KnnEntry<D> {
    fn eq(&self, other: &Self) -> bool {
        self.dist == other.dist
    }
}

impl<D> Eq for KnnEntry<D> {}

impl<D> PartialOrd for KnnEntry<D> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<D> Ord for KnnEntry<D> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.dist
            .partial_cmp(&other.dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

// ─── Parallel Spatial Hash ────────────────────────────────────────────────────

/// A parallel-build spatial hash.
///
/// All entries are inserted in one batch call, making the build trivially
/// parallelisable (simulated here sequentially).  Queries are identical to
/// `SpatialHashPos3D`.
pub struct ParallelSpatialHash<T: Clone> {
    inner: SpatialHashPos3D<T>,
}

impl<T: Clone> ParallelSpatialHash<T> {
    /// Build from a slice of `(id, position, data)` tuples.
    pub fn build(entries: &[(usize, [f64; 3], T)], cell_size: f64) -> Self {
        let mut inner = SpatialHashPos3D::new(cell_size);
        for (id, pos, data) in entries {
            inner.insert(*id, *pos, data.clone());
        }
        Self { inner }
    }

    /// Query all items within `radius` of `center`.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<(usize, &T)> {
        self.inner.query_radius(center, radius)
    }

    /// Query all items within an AABB.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<(usize, &T)> {
        self.inner.query_aabb(min, max)
    }

    /// Number of items.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether the hash is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Number of occupied cells.
    pub fn n_cells(&self) -> usize {
        self.inner.n_cells()
    }
}

// ─── Spatial Hash Merging ─────────────────────────────────────────────────────

/// Merge two `SpatialHashPos3D` instances (same cell size) into a new one.
///
/// IDs are *not* re-numbered; callers must ensure they are disjoint.
pub fn merge_spatial_hashes<T: Clone>(
    a: &SpatialHashPos3D<T>,
    b: &SpatialHashPos3D<T>,
) -> SpatialHashPos3D<T> {
    assert!(
        (a.cell_size - b.cell_size).abs() < 1e-12,
        "merge requires identical cell sizes"
    );
    let mut merged = SpatialHashPos3D::new(a.cell_size);

    // Copy all entries from a
    for (key, bucket) in &a.cells {
        for (id, pos, data) in bucket {
            // Insert directly into merged map
            merged
                .cells
                .entry(*key)
                .or_default()
                .push((*id, *pos, data.clone()));
            merged.n_items += 1;
        }
    }

    // Copy all entries from b
    for (key, bucket) in &b.cells {
        for (id, pos, data) in bucket {
            merged
                .cells
                .entry(*key)
                .or_default()
                .push((*id, *pos, data.clone()));
            merged.n_items += 1;
        }
    }

    merged
}

// ─── Persistent Spatial Hash ──────────────────────────────────────────────────

/// A spatial hash that persists across frames.
///
/// Tracks the last known position of each object so that `insert_or_update`
/// correctly removes the object from its old cell before re-inserting.
pub struct PersistentSpatialHash<T: Clone> {
    /// Underlying hash with positions.
    hash: SpatialHashPos3D<T>,
    /// Tracks last position of each object: `id → last_position`.
    last_positions: HashMap<usize, [f64; 3]>,
    /// Current frame number.
    frame: u64,
}

impl<T: Clone> PersistentSpatialHash<T> {
    /// Create a new persistent hash with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            hash: SpatialHashPos3D::new(cell_size),
            last_positions: HashMap::new(),
            frame: 0,
        }
    }

    /// Insert a new object or update an existing one.
    ///
    /// If the object already exists, it is first removed from its old cell.
    pub fn insert_or_update(&mut self, id: usize, pos: [f64; 3], data: T) {
        if let Some(&old_pos) = self.last_positions.get(&id) {
            self.hash.remove(id, old_pos);
        }
        self.hash.insert(id, pos, data);
        self.last_positions.insert(id, pos);
    }

    /// Remove an object.
    pub fn remove(&mut self, id: usize) {
        if let Some(&pos) = self.last_positions.get(&id) {
            self.hash.remove(id, pos);
            self.last_positions.remove(&id);
        }
    }

    /// Query all items within `radius` of `center`.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<(usize, &T)> {
        self.hash.query_radius(center, radius)
    }

    /// Query all items in an AABB.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<(usize, &T)> {
        self.hash.query_aabb(min, max)
    }

    /// Number of objects.
    pub fn len(&self) -> usize {
        self.hash.len()
    }

    /// Whether the hash is empty.
    pub fn is_empty(&self) -> bool {
        self.hash.is_empty()
    }

    /// Advance the frame counter (call once per simulation step).
    pub fn advance_frame(&mut self) {
        self.frame += 1;
    }

    /// Current frame number.
    pub fn current_frame(&self) -> u64 {
        self.frame
    }

    /// Compute statistics.
    pub fn statistics(&self) -> SpatialHashStats {
        self.hash.statistics()
    }
}

// ─── GPU-Friendly Spatial Hash Layout ────────────────────────────────────────

/// A spatial hash with an interleaved, flat-array layout for GPU upload.
///
/// Data is stored as sorted arrays of cell keys and per-item records,
/// matching the layout expected by GPU compute shaders.
pub struct GpuSpatialHashLayout<T: Clone> {
    cell_size: f32,
    /// Sorted cell keys.
    cell_keys: Vec<(i32, i32, i32)>,
    /// Cell start/end indices into `items`.
    cell_offsets: Vec<(usize, usize)>,
    /// Items sorted by cell key: (id, position, data).
    items: Vec<(usize, [f32; 3], T)>,
}

impl<T: Clone> GpuSpatialHashLayout<T> {
    /// Build the GPU-friendly layout from a slice of `(id, position, data)`.
    pub fn build(entries: &[(usize, [f32; 3], T)], cell_size: f32) -> Self {
        // Bucket entries by cell key
        let mut buckets: HashMap<(i32, i32, i32), Vec<(usize, [f32; 3], T)>> = HashMap::new();
        for (id, pos, data) in entries {
            let key = (
                (pos[0] / cell_size).floor() as i32,
                (pos[1] / cell_size).floor() as i32,
                (pos[2] / cell_size).floor() as i32,
            );
            buckets
                .entry(key)
                .or_default()
                .push((*id, *pos, data.clone()));
        }

        // Sort cell keys for cache-friendly GPU traversal
        let mut cell_keys: Vec<(i32, i32, i32)> = buckets.keys().copied().collect();
        cell_keys.sort_unstable();

        let mut items = Vec::new();
        let mut cell_offsets = Vec::new();

        for &key in &cell_keys {
            let start = items.len();
            for entry in &buckets[&key] {
                items.push(entry.clone());
            }
            cell_offsets.push((start, items.len()));
        }

        Self {
            cell_size,
            cell_keys,
            cell_offsets,
            items,
        }
    }

    /// Total number of items.
    pub fn n_items(&self) -> usize {
        self.items.len()
    }

    /// Number of non-empty cells.
    pub fn n_cells(&self) -> usize {
        self.cell_keys.len()
    }

    /// Query items in a specific cell.
    pub fn query_cell_flat(&self, cx: i32, cy: i32, cz: i32) -> Vec<(usize, &T)> {
        let key = (cx, cy, cz);
        if let Ok(pos) = self.cell_keys.binary_search(&key) {
            let (start, end) = self.cell_offsets[pos];
            self.items[start..end]
                .iter()
                .map(|(id, _, data)| (*id, data))
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Export all items as a flat `f32` buffer for GPU upload.
    ///
    /// Format per item: `[x, y, z, cell_x_f32, cell_y_f32, cell_z_f32]`
    pub fn to_flat_f32_buffer(&self) -> Vec<f32> {
        let mut buf = Vec::with_capacity(self.items.len() * 6);
        for (item_idx, (_, pos, _)) in self.items.iter().enumerate() {
            // Find which cell this item belongs to
            let key = self
                .cell_keys
                .iter()
                .enumerate()
                .find(|(ci, _)| {
                    let (start, end) = self.cell_offsets[*ci];
                    item_idx >= start && item_idx < end
                })
                .map(|(_, &k)| k)
                .unwrap_or((0, 0, 0));
            buf.push(pos[0]);
            buf.push(pos[1]);
            buf.push(pos[2]);
            buf.push(key.0 as f32);
            buf.push(key.1 as f32);
            buf.push(key.2 as f32);
        }
        buf
    }

    /// Cell size used for this layout.
    pub fn cell_size(&self) -> f32 {
        self.cell_size
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── SpatialHash3D ──────────────────────────────────────────────────────

    #[test]
    fn spatial_hash_basic_insert_query() {
        let mut sh: SpatialHash3D<&str> = SpatialHash3D::new(1.0);
        sh.insert(0, [0.5, 0.5, 0.5], "a");
        assert_eq!(sh.len(), 1);
        assert_eq!(sh.n_cells(), 1);
        let items = sh.query_cell((0, 0, 0));
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].0, 0);
    }

    #[test]
    fn spatial_hash_remove() {
        let mut sh: SpatialHash3D<i32> = SpatialHash3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], 42);
        sh.insert(1, [0.5, 0.0, 0.0], 7);
        sh.remove(0, [0.0, 0.0, 0.0]);
        assert_eq!(sh.len(), 1);
    }

    #[test]
    fn spatial_hash_cell_key() {
        let sh: SpatialHash3D<()> = SpatialHash3D::new(2.0);
        assert_eq!(sh.cell_key([0.5, 0.5, 0.5]), (0, 0, 0));
        assert_eq!(sh.cell_key([2.5, 4.5, -0.5]), (1, 2, -1));
    }

    // ── SpatialHashPos3D ───────────────────────────────────────────────────

    #[test]
    fn spatial_hash_insert_query_radius() {
        let mut sh: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], "origin");
        sh.insert(1, [0.5, 0.0, 0.0], "near");
        sh.insert(2, [10.0, 0.0, 0.0], "far");

        let results = sh.query_radius([0.0, 0.0, 0.0], 1.0);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0), "origin should be found");
        assert!(ids.contains(&1), "near should be found");
        assert!(!ids.contains(&2), "far should NOT be found");
    }

    #[test]
    fn spatial_hash_clear_len_zero() {
        let mut sh: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], 42);
        sh.insert(1, [1.0, 1.0, 1.0], 7);
        assert_eq!(sh.len(), 2);
        sh.clear();
        assert_eq!(sh.len(), 0);
        assert!(sh.is_empty());
    }

    #[test]
    fn spatial_hash_query_radius_excludes_outside() {
        let mut sh: SpatialHashPos3D<u32> = SpatialHashPos3D::new(2.0);
        let radius = 1.5_f64;
        sh.insert(0, [0.0, 0.0, 0.0], 0);
        sh.insert(1, [radius + 0.1, 0.0, 0.0], 1);
        sh.insert(2, [radius - 0.1, 0.0, 0.0], 2);

        let results = sh.query_radius([0.0, 0.0, 0.0], radius);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&2));
        assert!(
            !ids.contains(&1),
            "item just outside radius must be excluded"
        );
    }

    #[test]
    fn spatial_hash_query_aabb() {
        let mut sh: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.5, 0.5, 0.5], "inside");
        sh.insert(1, [1.5, 0.5, 0.5], "inside2");
        sh.insert(2, [5.0, 5.0, 5.0], "outside");

        let results = sh.query_aabb([0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&1));
        assert!(!ids.contains(&2));
    }

    #[test]
    fn spatial_hash_k_nearest() {
        let mut sh: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], "a");
        sh.insert(1, [1.0, 0.0, 0.0], "b");
        sh.insert(2, [2.0, 0.0, 0.0], "c");
        sh.insert(3, [10.0, 0.0, 0.0], "d");

        let nearest = sh.k_nearest([0.0, 0.0, 0.0], 2);
        assert_eq!(nearest.len(), 2);
        assert_eq!(nearest[0].0, 0, "closest should be item 0");
        assert_eq!(nearest[1].0, 1, "second closest should be item 1");
    }

    #[test]
    fn spatial_hash_k_nearest_more_than_items() {
        let mut sh: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], 1);

        let nearest = sh.k_nearest([0.0, 0.0, 0.0], 5);
        assert_eq!(nearest.len(), 1, "should return all items when k > n_items");
    }

    #[test]
    fn spatial_hash_statistics() {
        let mut sh: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], 1);
        sh.insert(1, [0.5, 0.0, 0.0], 2);
        sh.insert(2, [5.0, 5.0, 5.0], 3);

        let stats = sh.statistics();
        assert_eq!(stats.n_items, 3);
        assert_eq!(stats.n_occupied_cells, 2);
        assert_eq!(stats.max_bucket_size, 2);
        assert_eq!(stats.min_bucket_size, 1);
    }

    #[test]
    fn spatial_hash_statistics_empty() {
        let sh: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        let stats = sh.statistics();
        assert_eq!(stats.n_items, 0);
        assert_eq!(stats.n_occupied_cells, 0);
    }

    #[test]
    fn spatial_hash_update_position() {
        let mut sh: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        sh.insert(0, [0.0, 0.0, 0.0], "a");
        sh.update_position(0, [0.0, 0.0, 0.0], [5.0, 5.0, 5.0], "a");
        assert_eq!(sh.len(), 1);

        // Should no longer be at origin
        let at_origin = sh.query_radius([0.0, 0.0, 0.0], 0.5);
        assert!(at_origin.is_empty(), "item should no longer be at origin");

        // Should be at new position
        let at_new = sh.query_radius([5.0, 5.0, 5.0], 0.5);
        assert_eq!(at_new.len(), 1);
    }

    // ── Ray traversal ──────────────────────────────────────────────────────

    #[test]
    fn ray_traverse_along_x_axis() {
        let steps = ray_traverse_grid([0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1.0, 3.0, 100);
        assert!(!steps.is_empty());
        // Should visit cells (0,0,0), (1,0,0), (2,0,0), (3,0,0)
        let cells: Vec<(i32, i32, i32)> = steps.iter().map(|s| s.cell).collect();
        assert!(cells.contains(&(0, 0, 0)));
        assert!(cells.contains(&(1, 0, 0)));
        assert!(cells.contains(&(2, 0, 0)));
    }

    #[test]
    fn ray_traverse_diagonal() {
        let steps = ray_traverse_grid([0.5, 0.5, 0.5], [1.0, 1.0, 0.0], 1.0, 3.0, 100);
        assert!(!steps.is_empty());
        // First cell should be (0,0,0)
        assert_eq!(steps[0].cell, (0, 0, 0));
    }

    #[test]
    fn ray_traverse_zero_direction() {
        let steps = ray_traverse_grid([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0, 100);
        assert!(steps.is_empty(), "zero direction should produce no steps");
    }

    #[test]
    fn ray_traverse_negative_direction() {
        let steps = ray_traverse_grid([2.5, 0.5, 0.5], [-1.0, 0.0, 0.0], 1.0, 3.0, 100);
        assert!(!steps.is_empty());
        let cells: Vec<(i32, i32, i32)> = steps.iter().map(|s| s.cell).collect();
        assert!(cells.contains(&(2, 0, 0)));
        assert!(cells.contains(&(1, 0, 0)));
        assert!(cells.contains(&(0, 0, 0)));
    }

    #[test]
    fn ray_traverse_max_steps_limit() {
        let steps = ray_traverse_grid([0.5, 0.5, 0.5], [1.0, 0.0, 0.0], 1.0, 1000.0, 5);
        assert!(steps.len() <= 5, "should respect max_steps limit");
    }

    // ── LooseOctree ────────────────────────────────────────────────────────

    #[test]
    fn loose_octree_insert_query_sphere() {
        let mut tree: LooseOctree<&str> = LooseOctree::new([0.0; 3], 16.0, 4);
        tree.insert(0, [1.0, 0.0, 0.0], "a", 2);
        tree.insert(1, [0.0, 1.0, 0.0], "b", 2);
        tree.insert(2, [100.0, 0.0, 0.0], "far", 2);

        let results = tree.query_sphere([0.0, 0.0, 0.0], 5.0);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0), "item 0 should be found");
        assert!(ids.contains(&1), "item 1 should be found");
        assert!(!ids.contains(&2), "far item must not appear");
    }

    #[test]
    fn loose_octree_item_count() {
        let mut tree: LooseOctree<i32> = LooseOctree::new([0.0; 3], 32.0, 3);
        for i in 0..5_usize {
            tree.insert(i, [i as f64, 0.0, 0.0], i as i32, 2);
        }
        assert!(
            tree.item_count() >= 5,
            "item_count should be at least 5, got {}",
            tree.item_count()
        );
    }

    #[test]
    fn loose_octree_query_aabb() {
        let mut tree: LooseOctree<&str> = LooseOctree::new([0.0; 3], 16.0, 4);
        tree.insert(0, [1.0, 1.0, 1.0], "inside", 2);
        tree.insert(1, [-1.0, -1.0, -1.0], "inside2", 2);
        tree.insert(2, [50.0, 50.0, 50.0], "outside", 2);

        let results = tree.query_aabb([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&1));
        assert!(!ids.contains(&2));
    }

    #[test]
    fn loose_octree_k_nearest() {
        let mut tree: LooseOctree<&str> = LooseOctree::new([0.0; 3], 32.0, 4);
        tree.insert(0, [0.0, 0.0, 0.0], "origin", 2);
        tree.insert(1, [1.0, 0.0, 0.0], "near", 2);
        tree.insert(2, [5.0, 0.0, 0.0], "mid", 2);
        tree.insert(3, [100.0, 0.0, 0.0], "far", 2);

        let nearest = tree.k_nearest([0.0, 0.0, 0.0], 2);
        assert_eq!(nearest.len(), 2);
        assert_eq!(nearest[0].0, 0, "closest should be item 0");
        assert_eq!(nearest[1].0, 1, "second closest should be item 1");
    }

    #[test]
    fn loose_octree_max_depth() {
        let mut tree: LooseOctree<i32> = LooseOctree::new([0.0; 3], 32.0, 4);
        // Insert enough items to trigger subdivision
        for i in 0..10_usize {
            tree.insert(i, [i as f64 * 0.1, 0.0, 0.0], i as i32, 2);
        }
        let depth = tree.max_depth_reached();
        assert!(depth > 0, "should have subdivided at least once");
    }

    #[test]
    fn loose_octree_node_count() {
        let tree: LooseOctree<i32> = LooseOctree::new([0.0; 3], 32.0, 4);
        assert_eq!(tree.node_count(), 1, "empty tree has one root node");
    }

    #[test]
    fn loose_octree_node_count_after_subdivision() {
        let mut tree: LooseOctree<i32> = LooseOctree::new([0.0; 3], 32.0, 4);
        for i in 0..10_usize {
            tree.insert(i, [i as f64, 0.0, 0.0], i as i32, 2);
        }
        assert!(tree.node_count() > 1, "should have created child nodes");
    }

    // ── Parallel spatial hash tests ────────────────────────────────────

    #[test]
    fn parallel_spatial_hash_build_query() {
        let entries: Vec<(usize, [f64; 3], i32)> = (0..8)
            .map(|i| (i, [i as f64, 0.0, 0.0], i as i32))
            .collect();
        let sh = ParallelSpatialHash::build(&entries, 2.0);
        let results = sh.query_radius([0.0, 0.0, 0.0], 2.5);
        let ids: Vec<usize> = results.iter().map(|(id, _)| *id).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&1));
        assert!(ids.contains(&2));
    }

    #[test]
    fn parallel_spatial_hash_len() {
        let entries: Vec<(usize, [f64; 3], &str)> =
            vec![(0, [0.0, 0.0, 0.0], "a"), (1, [1.0, 0.0, 0.0], "b")];
        let sh = ParallelSpatialHash::build(&entries, 1.0);
        assert_eq!(sh.len(), 2);
        assert!(!sh.is_empty());
    }

    #[test]
    fn parallel_spatial_hash_empty() {
        let entries: Vec<(usize, [f64; 3], i32)> = Vec::new();
        let sh = ParallelSpatialHash::build(&entries, 1.0);
        assert!(sh.is_empty());
        assert_eq!(sh.len(), 0);
    }

    // ── Spatial hash merging tests ─────────────────────────────────────

    #[test]
    fn spatial_hash_merge_basic() {
        let mut a: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        a.insert(0, [0.0, 0.0, 0.0], 10);
        a.insert(1, [1.0, 0.0, 0.0], 20);

        let mut b: SpatialHashPos3D<i32> = SpatialHashPos3D::new(1.0);
        b.insert(2, [5.0, 0.0, 0.0], 30);
        b.insert(3, [6.0, 0.0, 0.0], 40);

        let merged = merge_spatial_hashes(&a, &b);
        assert_eq!(merged.len(), 4);
    }

    #[test]
    fn spatial_hash_merge_query_all() {
        let mut a: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        a.insert(0, [0.0, 0.0, 0.0], "a");

        let mut b: SpatialHashPos3D<&str> = SpatialHashPos3D::new(1.0);
        b.insert(1, [10.0, 0.0, 0.0], "b");

        let merged = merge_spatial_hashes(&a, &b);
        let r0 = merged.query_radius([0.0, 0.0, 0.0], 0.5);
        let r1 = merged.query_radius([10.0, 0.0, 0.0], 0.5);
        assert_eq!(r0.len(), 1);
        assert_eq!(r1.len(), 1);
    }

    // ── Persistent spatial hash tests ─────────────────────────────────

    #[test]
    fn persistent_hash_update_positions() {
        let mut ph: PersistentSpatialHash<i32> = PersistentSpatialHash::new(1.0);
        ph.insert_or_update(0, [0.0, 0.0, 0.0], 100);
        ph.insert_or_update(1, [5.0, 0.0, 0.0], 200);

        // Query frame 0
        let r = ph.query_radius([0.0, 0.0, 0.0], 0.5);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0, 0);

        // Move item 0 to new position
        ph.insert_or_update(0, [5.5, 0.0, 0.0], 100);

        let r2 = ph.query_radius([0.0, 0.0, 0.0], 0.5);
        assert!(r2.is_empty(), "item 0 should have moved");
    }

    #[test]
    fn persistent_hash_remove() {
        let mut ph: PersistentSpatialHash<i32> = PersistentSpatialHash::new(1.0);
        ph.insert_or_update(0, [0.0, 0.0, 0.0], 1);
        ph.insert_or_update(1, [0.5, 0.0, 0.0], 2);
        ph.remove(0);
        assert_eq!(ph.len(), 1);
    }

    #[test]
    fn persistent_hash_frame_advance() {
        let mut ph: PersistentSpatialHash<i32> = PersistentSpatialHash::new(1.0);
        ph.insert_or_update(0, [0.0, 0.0, 0.0], 42);
        ph.advance_frame();
        assert_eq!(ph.current_frame(), 1);
        // Data persists across frames
        let r = ph.query_radius([0.0, 0.0, 0.0], 0.5);
        assert_eq!(r.len(), 1);
    }

    // ── GPU-friendly layout tests ──────────────────────────────────────

    #[test]
    fn gpu_layout_build_basic() {
        let entries: Vec<(usize, [f32; 3], u32)> = vec![
            (0, [0.0, 0.0, 0.0], 10),
            (1, [1.0, 0.0, 0.0], 20),
            (2, [2.0, 0.0, 0.0], 30),
        ];
        let layout = GpuSpatialHashLayout::build(&entries, 1.0);
        assert_eq!(layout.n_items(), 3);
    }

    #[test]
    fn gpu_layout_query_cell() {
        let entries: Vec<(usize, [f32; 3], u32)> = vec![
            (0, [0.5, 0.0, 0.0], 1),
            (1, [0.7, 0.0, 0.0], 2),
            (2, [5.0, 0.0, 0.0], 3),
        ];
        let layout = GpuSpatialHashLayout::build(&entries, 1.0);
        let cell_items = layout.query_cell_flat(0, 0, 0);
        assert_eq!(cell_items.len(), 2, "two items in cell (0,0,0)");
    }

    #[test]
    fn gpu_layout_flat_buffer_format() {
        let entries: Vec<(usize, [f32; 3], u32)> = vec![(0, [1.0, 2.0, 3.0], 99)];
        let layout = GpuSpatialHashLayout::build(&entries, 1.0);
        let buf = layout.to_flat_f32_buffer();
        // Format: [x, y, z, cell_x as f32, cell_y as f32, cell_z as f32]
        assert!(!buf.is_empty());
        assert_eq!(buf[0], 1.0_f32);
        assert_eq!(buf[1], 2.0_f32);
        assert_eq!(buf[2], 3.0_f32);
    }
}
