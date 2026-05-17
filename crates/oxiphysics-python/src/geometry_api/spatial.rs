// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial data structures: heightfield, spatial hash, AABB, BVH.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::shapes::{add3, dot3, len3, normalize3, scale3, sub3};

// ---------------------------------------------------------------------------
// PyHeightfield
// ---------------------------------------------------------------------------

/// A heightfield terrain represented as a regular grid.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyHeightfield {
    /// Number of columns (X axis).
    #[pyo3(get, set)]
    pub cols: usize,
    /// Number of rows (Z axis).
    #[pyo3(get, set)]
    pub rows: usize,
    /// Height values in row-major order (rows × cols entries).
    #[pyo3(get)]
    pub data: Vec<f64>,
    /// World-space scale per cell in X.
    #[pyo3(get, set)]
    pub scale_x: f64,
    /// World-space scale per cell in Z.
    #[pyo3(get, set)]
    pub scale_z: f64,
    /// Vertical scale factor applied to raw height values.
    #[pyo3(get, set)]
    pub scale_y: f64,
}

#[pymethods]
impl PyHeightfield {
    /// Create a new flat heightfield.
    #[new]
    pub fn new(cols: usize, rows: usize, scale_x: f64, scale_z: f64, scale_y: f64) -> Self {
        Self {
            cols,
            rows,
            data: vec![0.0; cols * rows],
            scale_x,
            scale_z,
            scale_y,
        }
    }

    /// Set height at grid cell (col, row).
    pub fn set_height(&mut self, col: usize, row: usize, h: f64) {
        if col < self.cols && row < self.rows {
            self.data[row * self.cols + col] = h;
        }
    }

    /// Bilinearly interpolated height at world-space (x, z).
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let gx = x / self.scale_x;
        let gz = z / self.scale_z;
        let col = gx.floor() as isize;
        let row = gz.floor() as isize;
        let fx = gx - col as f64;
        let fz = gz - row as f64;

        let sample = |c: isize, r: isize| -> f64 {
            let c = c.clamp(0, self.cols as isize - 1) as usize;
            let r = r.clamp(0, self.rows as isize - 1) as usize;
            self.data[r * self.cols + c] * self.scale_y
        };

        let h00 = sample(col, row);
        let h10 = sample(col + 1, row);
        let h01 = sample(col, row + 1);
        let h11 = sample(col + 1, row + 1);

        let h0 = h00 + (h10 - h00) * fx;
        let h1 = h01 + (h11 - h01) * fx;
        h0 + (h1 - h0) * fz
    }

    /// Central-difference surface normal at world-space (x, z).
    pub fn normal_at(&self, x: f64, z: f64) -> Vec<f64> {
        let eps = self.scale_x.min(self.scale_z) * 0.5;
        let dhdx = (self.height_at(x + eps, z) - self.height_at(x - eps, z)) / (2.0 * eps);
        let dhdz = (self.height_at(x, z + eps) - self.height_at(x, z - eps)) / (2.0 * eps);
        normalize3([-dhdx, 1.0, -dhdz]).to_vec()
    }

    /// Raycast against the heightfield. Returns the hit distance if found.
    pub fn raycast(&self, origin: Vec<f64>, direction: Vec<f64>) -> Option<f64> {
        if origin.len() < 3 || direction.len() < 3 {
            return None;
        }
        let origin_arr = [origin[0], origin[1], origin[2]];
        let direction_arr = [direction[0], direction[1], direction[2]];
        let dir = normalize3(direction_arr);
        let max_dist = (self.cols as f64 * self.scale_x + self.rows as f64 * self.scale_z) * 2.0;
        let steps = 1000usize;
        let step = max_dist / steps as f64;
        let mut t = 0.0;
        for _ in 0..steps {
            let p = [
                origin_arr[0] + dir[0] * t,
                origin_arr[1] + dir[1] * t,
                origin_arr[2] + dir[2] * t,
            ];
            let h = self.height_at(p[0], p[2]);
            if p[1] <= h {
                return Some(t);
            }
            t += step;
        }
        None
    }
}

// ---------------------------------------------------------------------------
// PySpatialHash
// ---------------------------------------------------------------------------

/// A spatial hash map for fast point proximity queries.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySpatialHash {
    /// Cell size determines resolution of the hash grid.
    #[pyo3(get, set)]
    pub cell_size: f64,
    /// Points indexed by ID.
    pub points: Vec<[f64; 3]>,
    /// Internal grid map: cell key -> list of point indices.
    grid: HashMap<(i64, i64, i64), Vec<usize>>,
}

#[pymethods]
impl PySpatialHash {
    /// Create a new spatial hash with the given cell size.
    #[new]
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            points: vec![],
            grid: HashMap::new(),
        }
    }

    /// Insert a point [x, y, z] and return its ID.
    pub fn insert(&mut self, p: Vec<f64>) -> usize {
        if p.len() < 3 {
            return 0;
        }
        let pt = [p[0], p[1], p[2]];
        let id = self.points.len();
        self.points.push(pt);
        let cell = self.cell_of(pt);
        self.grid.entry(cell).or_default().push(id);
        id
    }

    /// Remove a point by ID (marks its cell entry as removed).
    pub fn remove(&mut self, id: usize) {
        if id >= self.points.len() {
            return;
        }
        let p = self.points[id];
        let cell = self.cell_of(p);
        if let Some(list) = self.grid.get_mut(&cell) {
            list.retain(|&x| x != id);
        }
    }

    /// Query all points in the cell containing `p`.
    pub fn query(&self, p: Vec<f64>) -> Vec<usize> {
        if p.len() < 3 {
            return vec![];
        }
        let pt = [p[0], p[1], p[2]];
        let cell = self.cell_of(pt);
        self.grid.get(&cell).cloned().unwrap_or_default()
    }

    /// Query all points within a sphere of the given `radius` around `center`.
    pub fn sphere_query(&self, center: Vec<f64>, radius: f64) -> Vec<usize> {
        if center.len() < 3 {
            return vec![];
        }
        let center_arr = [center[0], center[1], center[2]];
        let r2 = radius * radius;
        let cells = (radius / self.cell_size).ceil() as i64 + 1;
        let base = self.cell_of(center_arr);
        let mut result = Vec::new();
        for dx in -cells..=cells {
            for dy in -cells..=cells {
                for dz in -cells..=cells {
                    let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                    if let Some(list) = self.grid.get(&key) {
                        for &id in list {
                            if id < self.points.len() {
                                let d = sub3(self.points[id], center_arr);
                                if dot3(d, d) <= r2 {
                                    result.push(id);
                                }
                            }
                        }
                    }
                }
            }
        }
        result
    }

    /// Get all points as a list of [x, y, z] lists.
    pub fn get_points(&self) -> Vec<Vec<f64>> {
        self.points.iter().map(|p| p.to_vec()).collect()
    }
}

impl PySpatialHash {
    fn cell_of(&self, p: [f64; 3]) -> (i64, i64, i64) {
        (
            (p[0] / self.cell_size).floor() as i64,
            (p[1] / self.cell_size).floor() as i64,
            (p[2] / self.cell_size).floor() as i64,
        )
    }
}

// ---------------------------------------------------------------------------
// PyAabb
// ---------------------------------------------------------------------------

/// An axis-aligned bounding box.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyAabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

#[pymethods]
impl PyAabb {
    /// Create a new AABB from min/max lists.
    #[new]
    pub fn new(min: Vec<f64>, max: Vec<f64>) -> Self {
        let mn = if min.len() >= 3 {
            [min[0], min[1], min[2]]
        } else {
            [0.0; 3]
        };
        let mx = if max.len() >= 3 {
            [max[0], max[1], max[2]]
        } else {
            [0.0; 3]
        };
        Self { min: mn, max: mx }
    }

    /// Check whether this AABB overlaps another.
    pub fn overlaps(&self, other: &PyAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Center of this AABB.
    pub fn center(&self) -> Vec<f64> {
        scale3(add3(self.min, self.max), 0.5).to_vec()
    }

    /// Half-extents of this AABB.
    pub fn half_extents(&self) -> Vec<f64> {
        scale3(sub3(self.max, self.min), 0.5).to_vec()
    }

    /// Surface area of this AABB.
    pub fn surface_area(&self) -> f64 {
        let e = sub3(self.max, self.min);
        2.0 * (e[0] * e[1] + e[1] * e[2] + e[2] * e[0])
    }

    /// Get min as list.
    pub fn get_min(&self) -> Vec<f64> {
        self.min.to_vec()
    }

    /// Get max as list.
    pub fn get_max(&self) -> Vec<f64> {
        self.max.to_vec()
    }
}

impl PyAabb {
    /// Internal constructor from fixed arrays.
    pub fn new_internal(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Internal center as fixed array.
    pub fn center_internal(&self) -> [f64; 3] {
        scale3(add3(self.min, self.max), 0.5)
    }
}

// ---------------------------------------------------------------------------
// PyBvhNode and PyBvh
// ---------------------------------------------------------------------------

/// A node in the BVH tree.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyBvhNode {
    /// Bounding box of this node.
    pub aabb: PyAabb,
    /// Left child index (or `usize::MAX` for leaf).
    #[pyo3(get)]
    pub left: usize,
    /// Right child index (or `usize::MAX` for leaf).
    #[pyo3(get)]
    pub right: usize,
    /// Leaf item index (meaningful only when `left == usize::MAX`).
    #[pyo3(get)]
    pub item: usize,
}

#[pymethods]
impl PyBvhNode {
    /// Get the AABB of this node.
    pub fn get_aabb(&self) -> PyAabb {
        self.aabb.clone()
    }
}

/// A BVH (Bounding Volume Hierarchy) over AABBs.
#[pyclass(skip_from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyBvh {
    /// All nodes in the BVH.
    pub nodes: Vec<PyBvhNode>,
    /// Root node index.
    #[pyo3(get)]
    pub root: usize,
    /// Original AABBs associated with leaf items.
    pub items: Vec<PyAabb>,
}

#[pymethods]
impl PyBvh {
    /// Build a BVH from a list of AABBs.
    #[staticmethod]
    pub fn build(aabbs: Vec<PyAabb>) -> Self {
        let n = aabbs.len();
        if n == 0 {
            return Self {
                nodes: vec![],
                root: 0,
                items: vec![],
            };
        }
        let mut bvh = Self {
            nodes: Vec::with_capacity(2 * n),
            root: 0,
            items: aabbs,
        };
        let indices: Vec<usize> = (0..n).collect();
        bvh.root = bvh.build_recursive(&indices);
        bvh
    }

    /// Query all leaf items whose AABB overlaps the query AABB.
    pub fn query_aabb(&self, query: &PyAabb) -> Vec<usize> {
        if self.nodes.is_empty() {
            return vec![];
        }
        let mut result = Vec::new();
        self.query_recursive(self.root, query, &mut result);
        result
    }

    /// Raycast against all leaf items; returns sorted (distance, item) pairs.
    pub fn raycast(&self, origin: Vec<f64>, direction: Vec<f64>) -> Vec<(f64, usize)> {
        if self.nodes.is_empty() || origin.len() < 3 || direction.len() < 3 {
            return vec![];
        }
        let origin_arr = [origin[0], origin[1], origin[2]];
        let direction_arr = [direction[0], direction[1], direction[2]];
        let dir = normalize3(direction_arr);
        let mut hits = Vec::new();
        self.raycast_recursive(self.root, origin_arr, dir, &mut hits);
        hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }

    /// Nearest-neighbor search: returns the index of the closest item's center.
    pub fn nearest_neighbor(&self, point: Vec<f64>) -> Option<usize> {
        if self.items.is_empty() || point.len() < 3 {
            return None;
        }
        let pt = [point[0], point[1], point[2]];
        let mut best_dist = f64::MAX;
        let mut best = 0usize;
        for (i, item) in self.items.iter().enumerate() {
            let d = len3(sub3(item.center_internal(), pt));
            if d < best_dist {
                best_dist = d;
                best = i;
            }
        }
        Some(best)
    }

    /// Return the number of nodes in the BVH.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get all items as list of AABBs.
    pub fn get_items(&self) -> Vec<PyAabb> {
        self.items.clone()
    }
}

impl PyBvh {
    fn build_recursive(&mut self, indices: &[usize]) -> usize {
        if indices.len() == 1 {
            let id = indices[0];
            let node = PyBvhNode {
                aabb: self.items[id].clone(),
                left: usize::MAX,
                right: usize::MAX,
                item: id,
            };
            let idx = self.nodes.len();
            self.nodes.push(node);
            return idx;
        }
        let combined = merge_aabbs(indices.iter().map(|&i| &self.items[i]));
        let axis = longest_axis(&combined);
        let mut sorted = indices.to_vec();
        sorted.sort_by(|&a, &b| {
            let ca = self.items[a].center_internal()[axis];
            let cb = self.items[b].center_internal()[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = sorted.len() / 2;
        let left = self.build_recursive(&sorted[..mid]);
        let right = self.build_recursive(&sorted[mid..]);
        let node = PyBvhNode {
            aabb: combined,
            left,
            right,
            item: usize::MAX,
        };
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }

    fn query_recursive(&self, node_idx: usize, query: &PyAabb, out: &mut Vec<usize>) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if !node.aabb.overlaps(query) {
            return;
        }
        if node.left == usize::MAX {
            out.push(node.item);
            return;
        }
        self.query_recursive(node.left, query, out);
        self.query_recursive(node.right, query, out);
    }

    fn raycast_recursive(
        &self,
        node_idx: usize,
        origin: [f64; 3],
        dir: [f64; 3],
        out: &mut Vec<(f64, usize)>,
    ) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if let Some(t) = ray_aabb_intersect(origin, dir, &node.aabb) {
            if node.left == usize::MAX {
                out.push((t, node.item));
                return;
            }
            self.raycast_recursive(node.left, origin, dir, out);
            self.raycast_recursive(node.right, origin, dir, out);
        }
    }
}

// ---------------------------------------------------------------------------
// BVH internal helpers
// ---------------------------------------------------------------------------

fn merge_aabbs<'a>(iter: impl Iterator<Item = &'a PyAabb>) -> PyAabb {
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for a in iter {
        for i in 0..3 {
            if a.min[i] < mn[i] {
                mn[i] = a.min[i];
            }
            if a.max[i] > mx[i] {
                mx[i] = a.max[i];
            }
        }
    }
    PyAabb { min: mn, max: mx }
}

fn longest_axis(aabb: &PyAabb) -> usize {
    let e = sub3(aabb.max, aabb.min);
    if e[0] >= e[1] && e[0] >= e[2] {
        0
    } else if e[1] >= e[2] {
        1
    } else {
        2
    }
}

fn ray_aabb_intersect(origin: [f64; 3], dir: [f64; 3], aabb: &PyAabb) -> Option<f64> {
    let mut tmin = f64::NEG_INFINITY;
    let mut tmax = f64::INFINITY;
    for i in 0..3 {
        if dir[i].abs() < 1e-15 {
            if origin[i] < aabb.min[i] || origin[i] > aabb.max[i] {
                return None;
            }
        } else {
            let inv = 1.0 / dir[i];
            let t1 = (aabb.min[i] - origin[i]) * inv;
            let t2 = (aabb.max[i] - origin[i]) * inv;
            tmin = tmin.max(t1.min(t2));
            tmax = tmax.min(t1.max(t2));
        }
    }
    if tmax >= tmin && tmax >= 0.0 {
        Some(tmin.max(0.0))
    } else {
        None
    }
}
