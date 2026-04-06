#![allow(clippy::needless_range_loop, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Point cloud file format I/O.
//!
//! Provides readers and writers for common point cloud formats:
//!
//! - **PLY** format (ASCII and binary little-endian)
//! - **PCD** format (Point Cloud Data, PCL-compatible)
//! - **XYZ** with optional normals and colors
//! - **LAS/LAZ** header parsing (LiDAR point cloud metadata)
//! - **E57** basic structure definition
//!
//! Also provides point cloud processing utilities:
//!
//! - Point cloud data structure (positions, normals, colors, intensity)
//! - Voxel grid downsampling
//! - Statistical outlier removal
//! - Bounding box computation
//! - Point cloud merging / concatenation
//! - Normal estimation from k-nearest neighbors
//! - KD-tree for nearest neighbor queries
//! - Point cloud to mesh (ball pivoting concept)
//! - Coordinate system transformations

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::io::{self, Write};

// ---------------------------------------------------------------------------
// Core Point Cloud Data Structure
// ---------------------------------------------------------------------------

/// A single point with optional attributes.
#[derive(Debug, Clone, Default)]
pub struct PointCloudPoint {
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Normal vector `[nx, ny, nz]` (unit length, if present).
    pub normal: Option<[f64; 3]>,
    /// RGB color `[r, g, b]` in `[0, 255]`.
    pub color: Option<[u8; 3]>,
    /// Scalar intensity value.
    pub intensity: Option<f64>,
    /// Classification label.
    pub classification: Option<u8>,
}

/// A point cloud: a collection of points with optional per-point attributes.
#[derive(Debug, Clone)]
pub struct PointCloud {
    /// Positions of all points.
    pub positions: Vec<[f64; 3]>,
    /// Per-point normals (if present, same length as `positions`).
    pub normals: Option<Vec<[f64; 3]>>,
    /// Per-point RGB colors (if present).
    pub colors: Option<Vec<[u8; 3]>>,
    /// Per-point intensity values (if present).
    pub intensities: Option<Vec<f64>>,
    /// Per-point classification labels (if present).
    pub classifications: Option<Vec<u8>>,
}

impl PointCloud {
    /// Create an empty point cloud.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            normals: None,
            colors: None,
            intensities: None,
            classifications: None,
        }
    }

    /// Create a point cloud from positions only.
    pub fn from_positions(positions: Vec<[f64; 3]>) -> Self {
        Self {
            positions,
            normals: None,
            colors: None,
            intensities: None,
            classifications: None,
        }
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Whether normals are present.
    pub fn has_normals(&self) -> bool {
        self.normals.is_some()
    }

    /// Whether colors are present.
    pub fn has_colors(&self) -> bool {
        self.colors.is_some()
    }

    /// Whether intensities are present.
    pub fn has_intensities(&self) -> bool {
        self.intensities.is_some()
    }

    /// Add a point with position only.
    pub fn add_point(&mut self, pos: [f64; 3]) {
        self.positions.push(pos);
        if let Some(ref mut n) = self.normals {
            n.push([0.0; 3]);
        }
        if let Some(ref mut c) = self.colors {
            c.push([255, 255, 255]);
        }
        if let Some(ref mut i) = self.intensities {
            i.push(0.0);
        }
        if let Some(ref mut cl) = self.classifications {
            cl.push(0);
        }
    }

    /// Add a point with all attributes.
    pub fn add_full_point(&mut self, point: PointCloudPoint) {
        self.positions.push(point.position);
        if let Some(ref mut n) = self.normals {
            n.push(point.normal.unwrap_or([0.0; 3]));
        } else if let Some(n) = point.normal {
            let mut norms = vec![[0.0_f64; 3]; self.positions.len() - 1];
            norms.push(n);
            self.normals = Some(norms);
        }
        if let Some(ref mut c) = self.colors {
            c.push(point.color.unwrap_or([255, 255, 255]));
        } else if let Some(c) = point.color {
            let mut cols = vec![[255_u8; 3]; self.positions.len() - 1];
            cols.push(c);
            self.colors = Some(cols);
        }
        if let Some(ref mut i) = self.intensities {
            i.push(point.intensity.unwrap_or(0.0));
        } else if let Some(intensity) = point.intensity {
            let mut ints = vec![0.0_f64; self.positions.len() - 1];
            ints.push(intensity);
            self.intensities = Some(ints);
        }
        if let Some(ref mut cl) = self.classifications {
            cl.push(point.classification.unwrap_or(0));
        }
    }

    /// Get a single point as a `PointCloudPoint`.
    pub fn get_point(&self, idx: usize) -> Option<PointCloudPoint> {
        if idx >= self.positions.len() {
            return None;
        }
        Some(PointCloudPoint {
            position: self.positions[idx],
            normal: self.normals.as_ref().map(|n| n[idx]),
            color: self.colors.as_ref().map(|c| c[idx]),
            intensity: self.intensities.as_ref().map(|i| i[idx]),
            classification: self.classifications.as_ref().map(|c| c[idx]),
        })
    }

    /// Reserve capacity for `n` additional points.
    pub fn reserve(&mut self, n: usize) {
        self.positions.reserve(n);
        if let Some(ref mut normals) = self.normals {
            normals.reserve(n);
        }
        if let Some(ref mut colors) = self.colors {
            colors.reserve(n);
        }
        if let Some(ref mut ints) = self.intensities {
            ints.reserve(n);
        }
    }
}

impl Default for PointCloud {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Bounding Box
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box in 3D.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl BoundingBox {
    /// Compute the bounding box of a point cloud.
    pub fn from_point_cloud(cloud: &PointCloud) -> Option<Self> {
        if cloud.is_empty() {
            return None;
        }
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for p in &cloud.positions {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        Some(Self { min, max })
    }

    /// Center of the bounding box.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Extents (half-widths) of the bounding box.
    pub fn extents(&self) -> [f64; 3] {
        [
            0.5 * (self.max[0] - self.min[0]),
            0.5 * (self.max[1] - self.min[1]),
            0.5 * (self.max[2] - self.min[2]),
        ]
    }

    /// Diagonal length.
    pub fn diagonal(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Volume of the bounding box.
    pub fn volume(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        dx * dy * dz
    }

    /// Check if a point is inside the bounding box.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Merge two bounding boxes.
    pub fn merge(&self, other: &BoundingBox) -> BoundingBox {
        BoundingBox {
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
}

// ---------------------------------------------------------------------------
// KD-Tree for Nearest Neighbor Queries
// ---------------------------------------------------------------------------

/// A node in a KD-tree for 3D points.
#[derive(Debug)]
enum KdNode {
    /// Leaf node containing point indices.
    Leaf {
        /// Indices into the original point cloud.
        indices: Vec<usize>,
    },
    /// Internal node splitting on a given axis.
    Split {
        /// Axis (0=x, 1=y, 2=z).
        axis: usize,
        /// Split value.
        split_value: f64,
        /// Left subtree (points with coord <= split_value).
        left: Box<KdNode>,
        /// Right subtree.
        right: Box<KdNode>,
    },
}

/// A KD-tree for efficient nearest neighbor queries in 3D.
#[derive(Debug)]
pub struct KdTree {
    /// Root node.
    root: Option<KdNode>,
    /// Reference to positions (copied).
    positions: Vec<[f64; 3]>,
}

/// Maximum leaf size for the KD-tree.
const KD_LEAF_SIZE: usize = 16;

impl KdTree {
    /// Build a KD-tree from a point cloud.
    pub fn build(cloud: &PointCloud) -> Self {
        let positions = cloud.positions.clone();
        let indices: Vec<usize> = (0..positions.len()).collect();
        let root = if indices.is_empty() {
            None
        } else {
            Some(Self::build_recursive(&positions, indices, 0))
        };
        Self { root, positions }
    }

    /// Build from raw positions.
    pub fn from_positions(positions: &[[f64; 3]]) -> Self {
        let pos = positions.to_vec();
        let indices: Vec<usize> = (0..pos.len()).collect();
        let root = if indices.is_empty() {
            None
        } else {
            Some(Self::build_recursive(&pos, indices, 0))
        };
        Self {
            root,
            positions: pos,
        }
    }

    fn build_recursive(positions: &[[f64; 3]], mut indices: Vec<usize>, depth: usize) -> KdNode {
        if indices.len() <= KD_LEAF_SIZE {
            return KdNode::Leaf { indices };
        }

        let axis = depth % 3;
        indices.sort_by(|&a, &b| {
            positions[a][axis]
                .partial_cmp(&positions[b][axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mid = indices.len() / 2;
        let split_value = positions[indices[mid]][axis];

        let left_indices = indices[..mid].to_vec();
        let right_indices = indices[mid..].to_vec();

        KdNode::Split {
            axis,
            split_value,
            left: Box::new(Self::build_recursive(positions, left_indices, depth + 1)),
            right: Box::new(Self::build_recursive(positions, right_indices, depth + 1)),
        }
    }

    /// Find the nearest neighbor to a query point.
    ///
    /// Returns `(index, squared_distance)`.
    pub fn nearest(&self, query: [f64; 3]) -> Option<(usize, f64)> {
        let root = self.root.as_ref()?;
        let mut best_idx = 0;
        let mut best_dist_sq = f64::MAX;
        Self::nearest_recursive(
            root,
            &self.positions,
            query,
            &mut best_idx,
            &mut best_dist_sq,
        );
        Some((best_idx, best_dist_sq))
    }

    fn nearest_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        best_idx: &mut usize,
        best_dist_sq: &mut f64,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if d < *best_dist_sq {
                        *best_dist_sq = d;
                        *best_idx = i;
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::nearest_recursive(first, positions, query, best_idx, best_dist_sq);
                // Check if we need to search the other subtree
                if diff * diff < *best_dist_sq {
                    Self::nearest_recursive(second, positions, query, best_idx, best_dist_sq);
                }
            }
        }
    }

    /// Find k nearest neighbors.
    ///
    /// Returns a vector of `(index, squared_distance)` sorted by distance.
    pub fn k_nearest(&self, query: [f64; 3], k: usize) -> Vec<(usize, f64)> {
        if self.root.is_none() || k == 0 {
            return Vec::new();
        }
        let mut results: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        Self::knn_recursive(
            self.root.as_ref().expect("value should be Some"),
            &self.positions,
            query,
            k,
            &mut results,
        );
        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        results.truncate(k);
        results
    }

    fn knn_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        k: usize,
        results: &mut Vec<(usize, f64)>,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if results.len() < k {
                        results.push((i, d));
                        results.sort_by(|a, b| {
                            a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                        });
                    } else if d < results[k - 1].1 {
                        results[k - 1] = (i, d);
                        results.sort_by(|a, b| {
                            a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                        });
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::knn_recursive(first, positions, query, k, results);
                let worst_dist = if results.len() < k {
                    f64::MAX
                } else {
                    results[results.len() - 1].1
                };
                if diff * diff < worst_dist {
                    Self::knn_recursive(second, positions, query, k, results);
                }
            }
        }
    }

    /// Find all points within a given radius (squared).
    pub fn radius_search(&self, query: [f64; 3], radius_sq: f64) -> Vec<(usize, f64)> {
        let mut results = Vec::new();
        if let Some(ref root) = self.root {
            Self::radius_recursive(root, &self.positions, query, radius_sq, &mut results);
        }
        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        results
    }

    fn radius_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        radius_sq: f64,
        results: &mut Vec<(usize, f64)>,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if d <= radius_sq {
                        results.push((i, d));
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::radius_recursive(first, positions, query, radius_sq, results);
                if diff * diff <= radius_sq {
                    Self::radius_recursive(second, positions, query, radius_sq, results);
                }
            }
        }
    }

    /// Number of points in the tree.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the tree is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

/// Squared distance between two 3D points.
fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Euclidean distance between two 3D points.
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    dist_sq(a, b).sqrt()
}

/// Normalize a 3D vector. Returns `[0; 3]` if length < eps.
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Cross product of two 3D vectors.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Voxel Grid Downsampling
// ---------------------------------------------------------------------------

/// Downsample a point cloud using a voxel grid filter.
///
/// Each voxel cell of side `voxel_size` keeps only the centroid of points
/// falling within it.
pub fn voxel_grid_downsample(cloud: &PointCloud, voxel_size: f64) -> PointCloud {
    if cloud.is_empty() || voxel_size <= 0.0 {
        return cloud.clone();
    }

    let inv_size = 1.0 / voxel_size;
    let mut voxels: std::collections::HashMap<(i64, i64, i64), (Vec<usize>, [f64; 3], usize)> =
        std::collections::HashMap::new();

    for (i, p) in cloud.positions.iter().enumerate() {
        let vx = (p[0] * inv_size).floor() as i64;
        let vy = (p[1] * inv_size).floor() as i64;
        let vz = (p[2] * inv_size).floor() as i64;
        let key = (vx, vy, vz);
        let entry = voxels
            .entry(key)
            .or_insert_with(|| (Vec::new(), [0.0; 3], 0));
        entry.0.push(i);
        entry.1[0] += p[0];
        entry.1[1] += p[1];
        entry.1[2] += p[2];
        entry.2 += 1;
    }

    let mut result = PointCloud::new();
    let has_normals = cloud.has_normals();
    let has_colors = cloud.has_colors();
    let has_intensities = cloud.has_intensities();

    if has_normals {
        result.normals = Some(Vec::new());
    }
    if has_colors {
        result.colors = Some(Vec::new());
    }
    if has_intensities {
        result.intensities = Some(Vec::new());
    }

    for (indices, sum, count) in voxels.values() {
        let n = *count as f64;
        let centroid = [sum[0] / n, sum[1] / n, sum[2] / n];
        result.positions.push(centroid);

        if has_normals {
            let normals = cloud.normals.as_ref().expect("value should be Some");
            let mut avg_n = [0.0; 3];
            for &i in indices {
                avg_n[0] += normals[i][0];
                avg_n[1] += normals[i][1];
                avg_n[2] += normals[i][2];
            }
            result
                .normals
                .as_mut()
                .expect("value should be Some")
                .push(normalize3(avg_n));
        }

        if has_colors {
            let colors = cloud.colors.as_ref().expect("value should be Some");
            let mut avg_c = [0.0_f64; 3];
            for &i in indices {
                avg_c[0] += colors[i][0] as f64;
                avg_c[1] += colors[i][1] as f64;
                avg_c[2] += colors[i][2] as f64;
            }
            result.colors.as_mut().expect("value should be Some").push([
                (avg_c[0] / n).round().clamp(0.0, 255.0) as u8,
                (avg_c[1] / n).round().clamp(0.0, 255.0) as u8,
                (avg_c[2] / n).round().clamp(0.0, 255.0) as u8,
            ]);
        }

        if has_intensities {
            let ints = cloud.intensities.as_ref().expect("value should be Some");
            let avg_i: f64 = indices.iter().map(|&i| ints[i]).sum::<f64>() / n;
            result
                .intensities
                .as_mut()
                .expect("value should be Some")
                .push(avg_i);
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Statistical Outlier Removal
// ---------------------------------------------------------------------------

/// Remove statistical outliers from a point cloud.
///
/// For each point, computes the mean distance to its `k` nearest neighbors.
/// Points whose mean distance is more than `std_dev_multiplier` standard
/// deviations above the global mean are removed.
pub fn statistical_outlier_removal(
    cloud: &PointCloud,
    k: usize,
    std_dev_multiplier: f64,
) -> PointCloud {
    if cloud.len() <= k {
        return cloud.clone();
    }

    let tree = KdTree::build(cloud);
    let mut mean_dists = Vec::with_capacity(cloud.len());

    for (i, pos) in cloud.positions.iter().enumerate() {
        let neighbors = tree.k_nearest(*pos, k + 1);
        let sum: f64 = neighbors
            .iter()
            .filter(|(idx, _)| *idx != i)
            .take(k)
            .map(|(_, d_sq)| d_sq.sqrt())
            .sum();
        mean_dists.push(sum / k as f64);
    }

    let global_mean: f64 = mean_dists.iter().sum::<f64>() / mean_dists.len() as f64;
    let variance: f64 = mean_dists
        .iter()
        .map(|d| (d - global_mean).powi(2))
        .sum::<f64>()
        / mean_dists.len() as f64;
    let std_dev = variance.sqrt();
    let threshold = global_mean + std_dev_multiplier * std_dev;

    let mut result = PointCloud::new();
    if cloud.has_normals() {
        result.normals = Some(Vec::new());
    }
    if cloud.has_colors() {
        result.colors = Some(Vec::new());
    }
    if cloud.has_intensities() {
        result.intensities = Some(Vec::new());
    }

    for (i, &md) in mean_dists.iter().enumerate() {
        if md <= threshold {
            result.positions.push(cloud.positions[i]);
            if let Some(ref normals) = cloud.normals {
                result
                    .normals
                    .as_mut()
                    .expect("value should be Some")
                    .push(normals[i]);
            }
            if let Some(ref colors) = cloud.colors {
                result
                    .colors
                    .as_mut()
                    .expect("value should be Some")
                    .push(colors[i]);
            }
            if let Some(ref ints) = cloud.intensities {
                result
                    .intensities
                    .as_mut()
                    .expect("value should be Some")
                    .push(ints[i]);
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Normal Estimation
// ---------------------------------------------------------------------------

/// Estimate normals for a point cloud using k-nearest neighbors.
///
/// For each point, finds `k` nearest neighbors, computes the covariance
/// matrix, and takes the eigenvector with smallest eigenvalue as the normal.
/// Normals are oriented consistently by checking alignment with the centroid
/// direction.
pub fn estimate_normals(cloud: &PointCloud, k: usize) -> Vec<[f64; 3]> {
    if cloud.is_empty() {
        return Vec::new();
    }

    let tree = KdTree::build(cloud);
    let mut normals = Vec::with_capacity(cloud.len());

    // Compute centroid for normal orientation
    let centroid = compute_centroid(cloud);

    for (i, pos) in cloud.positions.iter().enumerate() {
        let neighbors = tree.k_nearest(*pos, k + 1);
        let neighbor_positions: Vec<[f64; 3]> = neighbors
            .iter()
            .filter(|(idx, _)| *idx != i)
            .take(k)
            .map(|(idx, _)| cloud.positions[*idx])
            .collect();

        if neighbor_positions.len() < 3 {
            normals.push([0.0, 0.0, 1.0]);
            continue;
        }

        // Compute mean of neighbors
        let n_neigh = neighbor_positions.len() as f64;
        let mut mean = [0.0; 3];
        for np in &neighbor_positions {
            mean[0] += np[0];
            mean[1] += np[1];
            mean[2] += np[2];
        }
        mean[0] /= n_neigh;
        mean[1] /= n_neigh;
        mean[2] /= n_neigh;

        // Compute covariance matrix (3x3, symmetric)
        let mut cov = [[0.0_f64; 3]; 3];
        for np in &neighbor_positions {
            let d = [np[0] - mean[0], np[1] - mean[1], np[2] - mean[2]];
            for r in 0..3 {
                for c in 0..3 {
                    cov[r][c] += d[r] * d[c];
                }
            }
        }
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] /= n_neigh;
            }
        }

        // Find smallest eigenvector using power iteration on (max_eval*I - cov)
        let normal = smallest_eigenvector_3x3(cov);

        // Orient normal to point away from centroid
        let to_centroid = [
            centroid[0] - pos[0],
            centroid[1] - pos[1],
            centroid[2] - pos[2],
        ];
        let dot_val =
            normal[0] * to_centroid[0] + normal[1] * to_centroid[1] + normal[2] * to_centroid[2];
        if dot_val > 0.0 {
            normals.push([-normal[0], -normal[1], -normal[2]]);
        } else {
            normals.push(normal);
        }
    }

    normals
}

/// Compute the centroid of a point cloud.
pub fn compute_centroid(cloud: &PointCloud) -> [f64; 3] {
    if cloud.is_empty() {
        return [0.0; 3];
    }
    let n = cloud.len() as f64;
    let mut c = [0.0; 3];
    for p in &cloud.positions {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Approximate smallest eigenvector of a 3x3 symmetric matrix using
/// inverse power iteration.
fn smallest_eigenvector_3x3(m: [[f64; 3]; 3]) -> [f64; 3] {
    // Estimate largest eigenvalue via power iteration
    let mut v = [1.0, 0.0, 0.0];
    for _ in 0..20 {
        let mv = mat_vec_3x3(m, v);
        let len = (mv[0] * mv[0] + mv[1] * mv[1] + mv[2] * mv[2]).sqrt();
        if len < 1e-15 {
            return v;
        }
        v = [mv[0] / len, mv[1] / len, mv[2] / len];
    }
    let max_eval =
        v[0] * mat_vec_3x3(m, v)[0] + v[1] * mat_vec_3x3(m, v)[1] + v[2] * mat_vec_3x3(m, v)[2];

    // Shifted matrix: B = max_eval * I - M
    let shift = max_eval + 1.0; // ensure positive definite
    let b = [
        [shift - m[0][0], -m[0][1], -m[0][2]],
        [-m[1][0], shift - m[1][1], -m[1][2]],
        [-m[2][0], -m[2][1], shift - m[2][2]],
    ];

    // Power iteration on B gives largest eigenvector of B = smallest of M
    let mut w = [0.577, 0.577, 0.577]; // ~(1,1,1)/sqrt(3)
    for _ in 0..30 {
        let bw = mat_vec_3x3(b, w);
        let len = (bw[0] * bw[0] + bw[1] * bw[1] + bw[2] * bw[2]).sqrt();
        if len < 1e-15 {
            return w;
        }
        w = [bw[0] / len, bw[1] / len, bw[2] / len];
    }
    w
}

/// Multiply 3x3 matrix by 3-vector.
fn mat_vec_3x3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// Point Cloud Merging
// ---------------------------------------------------------------------------

/// Merge two point clouds into one.
pub fn merge_point_clouds(a: &PointCloud, b: &PointCloud) -> PointCloud {
    let mut result = PointCloud::new();
    result.positions = a
        .positions
        .iter()
        .chain(b.positions.iter())
        .copied()
        .collect();

    if a.has_normals() || b.has_normals() {
        let mut normals = Vec::with_capacity(a.len() + b.len());
        if let Some(ref na) = a.normals {
            normals.extend_from_slice(na);
        } else {
            normals.extend(std::iter::repeat_n([0.0_f64; 3], a.len()));
        }
        if let Some(ref nb) = b.normals {
            normals.extend_from_slice(nb);
        } else {
            normals.extend(std::iter::repeat_n([0.0_f64; 3], b.len()));
        }
        result.normals = Some(normals);
    }

    if a.has_colors() || b.has_colors() {
        let mut colors = Vec::with_capacity(a.len() + b.len());
        if let Some(ref ca) = a.colors {
            colors.extend_from_slice(ca);
        } else {
            colors.extend(std::iter::repeat_n([255_u8; 3], a.len()));
        }
        if let Some(ref cb) = b.colors {
            colors.extend_from_slice(cb);
        } else {
            colors.extend(std::iter::repeat_n([255_u8; 3], b.len()));
        }
        result.colors = Some(colors);
    }

    if a.has_intensities() || b.has_intensities() {
        let mut ints = Vec::with_capacity(a.len() + b.len());
        if let Some(ref ia) = a.intensities {
            ints.extend_from_slice(ia);
        } else {
            ints.extend(std::iter::repeat_n(0.0_f64, a.len()));
        }
        if let Some(ref ib) = b.intensities {
            ints.extend_from_slice(ib);
        } else {
            ints.extend(std::iter::repeat_n(0.0_f64, b.len()));
        }
        result.intensities = Some(ints);
    }

    result
}

// ---------------------------------------------------------------------------
// Coordinate System Transformations
// ---------------------------------------------------------------------------

/// Apply a rigid transformation to a point cloud.
///
/// `rotation` is a 3x3 matrix (row-major), `translation` is `[tx, ty, tz]`.
pub fn transform_point_cloud(
    cloud: &PointCloud,
    rotation: [[f64; 3]; 3],
    translation: [f64; 3],
) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        let rotated = mat_vec_3x3(rotation, *pos);
        *pos = [
            rotated[0] + translation[0],
            rotated[1] + translation[1],
            rotated[2] + translation[2],
        ];
    }
    if let Some(ref mut normals) = result.normals {
        for n in normals.iter_mut() {
            *n = normalize3(mat_vec_3x3(rotation, *n));
        }
    }
    result
}

/// Scale a point cloud uniformly.
pub fn scale_point_cloud(cloud: &PointCloud, scale: f64) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        pos[0] *= scale;
        pos[1] *= scale;
        pos[2] *= scale;
    }
    result
}

/// Translate a point cloud.
pub fn translate_point_cloud(cloud: &PointCloud, offset: [f64; 3]) -> PointCloud {
    let mut result = cloud.clone();
    for pos in &mut result.positions {
        pos[0] += offset[0];
        pos[1] += offset[1];
        pos[2] += offset[2];
    }
    result
}

/// Center a point cloud at the origin.
pub fn center_point_cloud(cloud: &PointCloud) -> PointCloud {
    let c = compute_centroid(cloud);
    translate_point_cloud(cloud, [-c[0], -c[1], -c[2]])
}

// ---------------------------------------------------------------------------
// PLY Format I/O
// ---------------------------------------------------------------------------

/// PLY file encoding type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlyEncoding {
    /// ASCII text format.
    Ascii,
    /// Binary little-endian format.
    BinaryLittleEndian,
}

/// Write a point cloud to PLY ASCII format.
pub fn write_ply_ascii(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    writeln!(writer, "ply")?;
    writeln!(writer, "format ascii 1.0")?;
    writeln!(writer, "element vertex {}", cloud.len())?;
    writeln!(writer, "property float x")?;
    writeln!(writer, "property float y")?;
    writeln!(writer, "property float z")?;
    if cloud.has_normals() {
        writeln!(writer, "property float nx")?;
        writeln!(writer, "property float ny")?;
        writeln!(writer, "property float nz")?;
    }
    if cloud.has_colors() {
        writeln!(writer, "property uchar red")?;
        writeln!(writer, "property uchar green")?;
        writeln!(writer, "property uchar blue")?;
    }
    if cloud.has_intensities() {
        writeln!(writer, "property float intensity")?;
    }
    writeln!(writer, "end_header")?;

    for i in 0..cloud.len() {
        let p = cloud.positions[i];
        write!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
        if let Some(ref normals) = cloud.normals {
            let n = normals[i];
            write!(writer, " {:.6} {:.6} {:.6}", n[0], n[1], n[2])?;
        }
        if let Some(ref colors) = cloud.colors {
            let c = colors[i];
            write!(writer, " {} {} {}", c[0], c[1], c[2])?;
        }
        if let Some(ref ints) = cloud.intensities {
            write!(writer, " {:.6}", ints[i])?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Write a point cloud to PLY binary little-endian format (to a byte buffer).
pub fn write_ply_binary_le(cloud: &PointCloud) -> Vec<u8> {
    let mut buf = Vec::new();

    let header = format!(
        "ply\nformat binary_little_endian 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nend_header\n",
        cloud.len()
    );
    buf.extend_from_slice(header.as_bytes());

    for p in &cloud.positions {
        buf.extend_from_slice(&(p[0] as f32).to_le_bytes());
        buf.extend_from_slice(&(p[1] as f32).to_le_bytes());
        buf.extend_from_slice(&(p[2] as f32).to_le_bytes());
    }
    buf
}

/// Parse a PLY ASCII string into a point cloud.
pub fn parse_ply_ascii(data: &str) -> Option<PointCloud> {
    let mut lines = data.lines();
    let first = lines.next()?;
    if first.trim() != "ply" {
        return None;
    }

    let mut vertex_count = 0usize;
    let mut has_normals = false;
    let mut has_colors = false;
    let mut _has_intensity = false;
    let mut in_header = true;

    while in_header {
        let line = lines.next()?;
        let line = line.trim();
        if line.starts_with("element vertex") {
            vertex_count = line.split_whitespace().last()?.parse().ok()?;
        } else if line.starts_with("property") {
            if line.contains(" nx") || line.contains(" ny") || line.contains(" nz") {
                has_normals = true;
            }
            if line.contains(" red") || line.contains(" green") || line.contains(" blue") {
                has_colors = true;
            }
            if line.contains(" intensity") {
                _has_intensity = true;
            }
        } else if line == "end_header" {
            in_header = false;
        }
    }

    let mut cloud = PointCloud::new();
    cloud.positions.reserve(vertex_count);
    if has_normals {
        cloud.normals = Some(Vec::with_capacity(vertex_count));
    }
    if has_colors {
        cloud.colors = Some(Vec::with_capacity(vertex_count));
    }

    for _ in 0..vertex_count {
        let line = lines.next()?;
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }
        let x: f64 = vals[0].parse().ok()?;
        let y: f64 = vals[1].parse().ok()?;
        let z: f64 = vals[2].parse().ok()?;
        cloud.positions.push([x, y, z]);

        let mut idx = 3;
        if has_normals && vals.len() >= idx + 3 {
            let nx: f64 = vals[idx].parse().ok()?;
            let ny: f64 = vals[idx + 1].parse().ok()?;
            let nz: f64 = vals[idx + 2].parse().ok()?;
            cloud
                .normals
                .as_mut()
                .expect("value should be Some")
                .push([nx, ny, nz]);
            idx += 3;
        }
        if has_colors && vals.len() >= idx + 3 {
            let r: u8 = vals[idx].parse().ok()?;
            let g: u8 = vals[idx + 1].parse().ok()?;
            let b: u8 = vals[idx + 2].parse().ok()?;
            cloud
                .colors
                .as_mut()
                .expect("value should be Some")
                .push([r, g, b]);
        }
    }

    Some(cloud)
}

// ---------------------------------------------------------------------------
// PCD Format I/O (PCL-compatible)
// ---------------------------------------------------------------------------

/// Write a point cloud to PCD ASCII format.
pub fn write_pcd_ascii(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    writeln!(writer, "# .PCD v0.7 - Point Cloud Data")?;
    writeln!(writer, "VERSION 0.7")?;
    writeln!(writer, "FIELDS x y z")?;
    writeln!(writer, "SIZE 4 4 4")?;
    writeln!(writer, "TYPE F F F")?;
    writeln!(writer, "COUNT 1 1 1")?;
    writeln!(writer, "WIDTH {}", cloud.len())?;
    writeln!(writer, "HEIGHT 1")?;
    writeln!(writer, "VIEWPOINT 0 0 0 1 0 0 0")?;
    writeln!(writer, "POINTS {}", cloud.len())?;
    writeln!(writer, "DATA ascii")?;

    for p in &cloud.positions {
        writeln!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
    }
    Ok(())
}

/// Parse a PCD ASCII string into a point cloud.
pub fn parse_pcd_ascii(data: &str) -> Option<PointCloud> {
    let mut lines = data.lines();
    let mut num_points = 0usize;
    let mut data_started = false;

    let mut cloud = PointCloud::new();

    for line in &mut lines {
        let line = line.trim();
        if line.starts_with("POINTS") {
            num_points = line.split_whitespace().last()?.parse().ok()?;
            cloud.positions.reserve(num_points);
        } else if line.starts_with("DATA") {
            data_started = true;
            break;
        }
    }

    if !data_started {
        return None;
    }

    for _ in 0..num_points {
        let line = lines.next()?;
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }
        let x: f64 = vals[0].parse().ok()?;
        let y: f64 = vals[1].parse().ok()?;
        let z: f64 = vals[2].parse().ok()?;
        cloud.positions.push([x, y, z]);
    }

    Some(cloud)
}

// ---------------------------------------------------------------------------
// XYZ Format I/O
// ---------------------------------------------------------------------------

/// Write a point cloud to XYZ format (with optional normals and colors).
pub fn write_xyz(cloud: &PointCloud, writer: &mut dyn Write) -> io::Result<()> {
    for i in 0..cloud.len() {
        let p = cloud.positions[i];
        write!(writer, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
        if let Some(ref normals) = cloud.normals {
            let n = normals[i];
            write!(writer, " {:.6} {:.6} {:.6}", n[0], n[1], n[2])?;
        }
        if let Some(ref colors) = cloud.colors {
            let c = colors[i];
            write!(writer, " {} {} {}", c[0], c[1], c[2])?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Parse an XYZ file (space-separated x y z \[nx ny nz\] \[r g b\]).
pub fn parse_xyz(data: &str) -> PointCloud {
    let mut cloud = PointCloud::new();
    let mut first_line_fields = 0;

    for (line_idx, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let vals: Vec<&str> = line.split_whitespace().collect();
        if vals.len() < 3 {
            continue;
        }

        if line_idx == 0 || first_line_fields == 0 {
            first_line_fields = vals.len();
            if first_line_fields >= 6 {
                cloud.normals = Some(Vec::new());
            }
            if first_line_fields >= 9 {
                cloud.colors = Some(Vec::new());
            }
        }

        let x: f64 = vals[0].parse().unwrap_or(0.0);
        let y: f64 = vals[1].parse().unwrap_or(0.0);
        let z: f64 = vals[2].parse().unwrap_or(0.0);
        cloud.positions.push([x, y, z]);

        if vals.len() >= 6
            && let Some(ref mut normals) = cloud.normals
        {
            let nx: f64 = vals[3].parse().unwrap_or(0.0);
            let ny: f64 = vals[4].parse().unwrap_or(0.0);
            let nz: f64 = vals[5].parse().unwrap_or(0.0);
            normals.push([nx, ny, nz]);
        }

        if vals.len() >= 9
            && let Some(ref mut colors) = cloud.colors
        {
            let r: u8 = vals[6].parse().unwrap_or(0);
            let g: u8 = vals[7].parse().unwrap_or(0);
            let b: u8 = vals[8].parse().unwrap_or(0);
            colors.push([r, g, b]);
        }
    }
    cloud
}

// ---------------------------------------------------------------------------
// LAS Header Parsing (LiDAR)
// ---------------------------------------------------------------------------

/// LAS file header structure (subset of fields).
#[derive(Debug, Clone)]
pub struct LasHeader {
    /// File signature (should be "LASF").
    pub file_signature: String,
    /// File source ID.
    pub file_source_id: u16,
    /// Version major.
    pub version_major: u8,
    /// Version minor.
    pub version_minor: u8,
    /// Point data record format.
    pub point_format: u8,
    /// Point data record length.
    pub point_record_length: u16,
    /// Number of point records.
    pub num_points: u64,
    /// Scale factors `[sx, sy, sz]`.
    pub scale: [f64; 3],
    /// Offset `[ox, oy, oz]`.
    pub offset: [f64; 3],
    /// Min bounds `[min_x, min_y, min_z]`.
    pub min_bounds: [f64; 3],
    /// Max bounds `[max_x, max_y, max_z]`.
    pub max_bounds: [f64; 3],
}

impl LasHeader {
    /// Create a mock LAS header for testing.
    pub fn mock(num_points: u64) -> Self {
        Self {
            file_signature: "LASF".to_string(),
            file_source_id: 0,
            version_major: 1,
            version_minor: 4,
            point_format: 0,
            point_record_length: 20,
            num_points,
            scale: [0.001, 0.001, 0.001],
            offset: [0.0, 0.0, 0.0],
            min_bounds: [0.0, 0.0, 0.0],
            max_bounds: [100.0, 100.0, 50.0],
        }
    }

    /// Parse a LAS header from raw bytes (minimal parsing).
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 227 {
            return None;
        }
        let sig = std::str::from_utf8(&data[0..4]).ok()?;
        if sig != "LASF" {
            return None;
        }
        let file_source_id = u16::from_le_bytes([data[4], data[5]]);
        let version_major = data[24];
        let version_minor = data[25];
        let point_format = data[104];
        let point_record_length = u16::from_le_bytes([data[105], data[106]]);

        // Number of point records (legacy 32-bit field at offset 107)
        let num_points_legacy =
            u32::from_le_bytes([data[107], data[108], data[109], data[110]]) as u64;

        let parse_f64 = |offset: usize| -> f64 {
            if offset + 8 > data.len() {
                return 0.0;
            }
            let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap_or([0; 8]);
            f64::from_le_bytes(bytes)
        };

        let scale = [parse_f64(131), parse_f64(139), parse_f64(147)];
        let offset = [parse_f64(155), parse_f64(163), parse_f64(171)];
        let max_x = parse_f64(179);
        let min_x = parse_f64(187);
        let max_y = parse_f64(195);
        let min_y = parse_f64(203);
        let max_z = parse_f64(211);
        let min_z = parse_f64(219);

        Some(Self {
            file_signature: sig.to_string(),
            file_source_id,
            version_major,
            version_minor,
            point_format,
            point_record_length,
            num_points: num_points_legacy,
            scale,
            offset,
            min_bounds: [min_x, min_y, min_z],
            max_bounds: [max_x, max_y, max_z],
        })
    }

    /// Bounding box volume.
    pub fn bounding_volume(&self) -> f64 {
        let dx = self.max_bounds[0] - self.min_bounds[0];
        let dy = self.max_bounds[1] - self.min_bounds[1];
        let dz = self.max_bounds[2] - self.min_bounds[2];
        dx.abs() * dy.abs() * dz.abs()
    }
}

// ---------------------------------------------------------------------------
// E57 Basic Structure
// ---------------------------------------------------------------------------

/// E57 scan header (simplified).
#[derive(Debug, Clone)]
pub struct E57ScanHeader {
    /// Scan name/description.
    pub name: String,
    /// Number of points in the scan.
    pub num_points: u64,
    /// Whether the scan has color data.
    pub has_color: bool,
    /// Whether the scan has intensity data.
    pub has_intensity: bool,
    /// Scanner position `[x, y, z]`.
    pub scanner_position: [f64; 3],
    /// Scanner orientation (quaternion `[x, y, z, w]`).
    pub scanner_orientation: [f64; 4],
}

impl E57ScanHeader {
    /// Create a mock E57 scan header for testing.
    pub fn mock(name: &str, num_points: u64) -> Self {
        Self {
            name: name.to_string(),
            num_points,
            has_color: true,
            has_intensity: true,
            scanner_position: [0.0; 3],
            scanner_orientation: [0.0, 0.0, 0.0, 1.0],
        }
    }
}

/// E57 file structure (simplified).
#[derive(Debug, Clone)]
pub struct E57File {
    /// Major version.
    pub version_major: u32,
    /// Minor version.
    pub version_minor: u32,
    /// Scans in the file.
    pub scans: Vec<E57ScanHeader>,
}

impl E57File {
    /// Create a mock E57 file for testing.
    pub fn mock() -> Self {
        Self {
            version_major: 1,
            version_minor: 0,
            scans: vec![E57ScanHeader::mock("scan_001", 1000)],
        }
    }

    /// Total number of points across all scans.
    pub fn total_points(&self) -> u64 {
        self.scans.iter().map(|s| s.num_points).sum()
    }

    /// Number of scans.
    pub fn num_scans(&self) -> usize {
        self.scans.len()
    }
}

// ---------------------------------------------------------------------------
// Ball Pivoting (Simplified Concept)
// ---------------------------------------------------------------------------

/// A triangle in a mesh (indices into a point cloud).
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// Vertex indices.
    pub indices: [usize; 3],
}

/// Simplified ball-pivoting mesh generation.
///
/// Given a point cloud with normals, attempts to create a triangle mesh by
/// pivoting a ball of radius `ball_radius` over the point cloud surface.
///
/// This is a conceptual implementation suitable for small point clouds.
pub fn ball_pivoting_mesh(cloud: &PointCloud, ball_radius: f64) -> Vec<Triangle> {
    if cloud.len() < 3 {
        return Vec::new();
    }

    let tree = KdTree::build(cloud);
    let mut triangles = Vec::new();
    let mut used = vec![false; cloud.len()];
    let search_radius_sq = (2.0 * ball_radius) * (2.0 * ball_radius);

    // Seed triangle: find 3 nearest points that form a valid triangle
    let neighbors = tree.k_nearest(cloud.positions[0], 20.min(cloud.len()));
    let mut found_seed = false;
    let mut seed = [0usize; 3];

    'outer: for i in 0..neighbors.len() {
        for j in (i + 1)..neighbors.len() {
            for k in (j + 1)..neighbors.len() {
                let a = neighbors[i].0;
                let b = neighbors[j].0;
                let c = neighbors[k].0;
                let ab = dist(cloud.positions[a], cloud.positions[b]);
                let bc = dist(cloud.positions[b], cloud.positions[c]);
                let ca = dist(cloud.positions[c], cloud.positions[a]);
                if ab < 2.0 * ball_radius && bc < 2.0 * ball_radius && ca < 2.0 * ball_radius {
                    seed = [a, b, c];
                    found_seed = true;
                    break 'outer;
                }
            }
        }
    }

    if !found_seed {
        return triangles;
    }

    triangles.push(Triangle { indices: seed });
    used[seed[0]] = true;
    used[seed[1]] = true;
    used[seed[2]] = true;

    // Simple expansion: for each edge of existing triangles, try to find a
    // third point that forms a valid triangle
    let max_triangles = cloud.len().min(5000);
    let mut edge_queue: Vec<(usize, usize)> =
        vec![(seed[0], seed[1]), (seed[1], seed[2]), (seed[2], seed[0])];

    while let Some((a, b)) = edge_queue.pop() {
        if triangles.len() >= max_triangles {
            break;
        }
        let mid = [
            0.5 * (cloud.positions[a][0] + cloud.positions[b][0]),
            0.5 * (cloud.positions[a][1] + cloud.positions[b][1]),
            0.5 * (cloud.positions[a][2] + cloud.positions[b][2]),
        ];
        let nearby = tree.radius_search(mid, search_radius_sq);
        for (c, _) in nearby {
            if c == a || c == b || used[c] {
                continue;
            }
            let ac = dist(cloud.positions[a], cloud.positions[c]);
            let bc = dist(cloud.positions[b], cloud.positions[c]);
            if ac < 2.0 * ball_radius && bc < 2.0 * ball_radius {
                triangles.push(Triangle { indices: [a, b, c] });
                used[c] = true;
                edge_queue.push((a, c));
                edge_queue.push((c, b));
                break;
            }
        }
    }

    triangles
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_cube_cloud(n_per_face: usize) -> PointCloud {
        let mut cloud = PointCloud::new();
        let step = 1.0 / n_per_face as f64;
        // Generate points on 6 faces of a unit cube
        for i in 0..n_per_face {
            for j in 0..n_per_face {
                let u = i as f64 * step;
                let v = j as f64 * step;
                cloud.positions.push([u, v, 0.0]); // bottom
                cloud.positions.push([u, v, 1.0]); // top
                cloud.positions.push([0.0, u, v]); // left
                cloud.positions.push([1.0, u, v]); // right
                cloud.positions.push([u, 0.0, v]); // front
                cloud.positions.push([u, 1.0, v]); // back
            }
        }
        cloud
    }

    fn make_simple_cloud() -> PointCloud {
        PointCloud::from_positions(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ])
    }

    #[test]
    fn test_point_cloud_new() {
        let cloud = PointCloud::new();
        assert!(cloud.is_empty());
        assert_eq!(cloud.len(), 0);
    }

    #[test]
    fn test_point_cloud_from_positions() {
        let cloud = make_simple_cloud();
        assert_eq!(cloud.len(), 8);
        assert!(!cloud.is_empty());
        assert!(!cloud.has_normals());
        assert!(!cloud.has_colors());
    }

    #[test]
    fn test_point_cloud_add_point() {
        let mut cloud = PointCloud::new();
        cloud.add_point([1.0, 2.0, 3.0]);
        assert_eq!(cloud.len(), 1);
        assert!((cloud.positions[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_cloud_get_point() {
        let cloud = make_simple_cloud();
        let pt = cloud.get_point(0).unwrap();
        assert!((pt.position[0]).abs() < 1e-10);
        assert!(cloud.get_point(100).is_none());
    }

    #[test]
    fn test_bounding_box_basic() {
        let cloud = make_simple_cloud();
        let bb = BoundingBox::from_point_cloud(&cloud).unwrap();
        assert!((bb.min[0]).abs() < 1e-10);
        assert!((bb.max[0] - 1.0).abs() < 1e-10);
        assert!((bb.max[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_box_center() {
        let cloud = make_simple_cloud();
        let bb = BoundingBox::from_point_cloud(&cloud).unwrap();
        let c = bb.center();
        assert!((c[0] - 0.5).abs() < 1e-10);
        assert!((c[1] - 0.5).abs() < 1e-10);
        assert!((c[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_box_volume() {
        let cloud = make_simple_cloud();
        let bb = BoundingBox::from_point_cloud(&cloud).unwrap();
        assert!((bb.volume() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_box_diagonal() {
        let cloud = make_simple_cloud();
        let bb = BoundingBox::from_point_cloud(&cloud).unwrap();
        let expected = (3.0_f64).sqrt();
        assert!((bb.diagonal() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_bounding_box_contains() {
        let cloud = make_simple_cloud();
        let bb = BoundingBox::from_point_cloud(&cloud).unwrap();
        assert!(bb.contains([0.5, 0.5, 0.5]));
        assert!(!bb.contains([2.0, 0.0, 0.0]));
    }

    #[test]
    fn test_bounding_box_merge() {
        let bb1 = BoundingBox {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let bb2 = BoundingBox {
            min: [-1.0, -1.0, -1.0],
            max: [0.5, 0.5, 0.5],
        };
        let merged = bb1.merge(&bb2);
        assert!((merged.min[0] + 1.0).abs() < 1e-10);
        assert!((merged.max[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kdtree_nearest() {
        let cloud = make_simple_cloud();
        let tree = KdTree::build(&cloud);
        let (idx, d_sq) = tree.nearest([0.1, 0.1, 0.1]).unwrap();
        assert_eq!(idx, 0); // nearest to origin
        assert!(d_sq < 0.1);
    }

    #[test]
    fn test_kdtree_k_nearest() {
        let cloud = make_simple_cloud();
        let tree = KdTree::build(&cloud);
        let results = tree.k_nearest([0.0, 0.0, 0.0], 3);
        assert_eq!(results.len(), 3);
        // First should be the origin itself
        assert_eq!(results[0].0, 0);
        assert!(results[0].1 < 1e-10);
    }

    #[test]
    fn test_kdtree_radius_search() {
        let cloud = make_simple_cloud();
        let tree = KdTree::build(&cloud);
        let results = tree.radius_search([0.0, 0.0, 0.0], 1.01); // radius^2 ~1
        // Should find the origin and 3 neighbors at distance 1
        assert!(!results.is_empty());
        assert!(results.len() <= 4);
    }

    #[test]
    fn test_kdtree_empty() {
        let cloud = PointCloud::new();
        let tree = KdTree::build(&cloud);
        assert!(tree.is_empty());
        assert!(tree.nearest([0.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn test_voxel_grid_downsample() {
        let cloud = make_cube_cloud(10); // 600 points
        assert!(cloud.len() > 100);
        let downsampled = voxel_grid_downsample(&cloud, 0.5);
        assert!(downsampled.len() < cloud.len());
        assert!(!downsampled.is_empty());
    }

    #[test]
    fn test_voxel_grid_preserves_bounds() {
        let cloud = make_simple_cloud();
        let bb_orig = BoundingBox::from_point_cloud(&cloud).unwrap();
        let downsampled = voxel_grid_downsample(&cloud, 0.5);
        let bb_down = BoundingBox::from_point_cloud(&downsampled).unwrap();
        // Downsampled bounds should be within original bounds (with some tolerance)
        assert!(bb_down.min[0] >= bb_orig.min[0] - 0.5);
        assert!(bb_down.max[0] <= bb_orig.max[0] + 0.5);
    }

    #[test]
    fn test_statistical_outlier_removal() {
        let mut cloud = make_simple_cloud();
        // Add an outlier far away
        cloud.positions.push([100.0, 100.0, 100.0]);
        let filtered = statistical_outlier_removal(&cloud, 3, 1.0);
        // The outlier should be removed
        assert!(filtered.len() <= cloud.len());
        // Check that the outlier is gone
        for p in &filtered.positions {
            assert!(p[0] < 50.0, "Outlier not removed: {p:?}");
        }
    }

    #[test]
    fn test_estimate_normals() {
        // Plane z=0
        let mut cloud = PointCloud::new();
        for i in 0..10 {
            for j in 0..10 {
                cloud.positions.push([i as f64, j as f64, 0.0]);
            }
        }
        let normals = estimate_normals(&cloud, 6);
        assert_eq!(normals.len(), 100);
        // Most normals should point roughly along z-axis
        let mut z_aligned = 0;
        for n in &normals {
            if n[2].abs() > 0.5 {
                z_aligned += 1;
            }
        }
        assert!(z_aligned > 50, "Only {z_aligned}/100 normals z-aligned");
    }

    #[test]
    fn test_compute_centroid() {
        let cloud = make_simple_cloud();
        let c = compute_centroid(&cloud);
        assert!((c[0] - 0.5).abs() < 1e-10);
        assert!((c[1] - 0.5).abs() < 1e-10);
        assert!((c[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_merge_point_clouds() {
        let a = PointCloud::from_positions(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let b = PointCloud::from_positions(vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let merged = merge_point_clouds(&a, &b);
        assert_eq!(merged.len(), 4);
        assert!((merged.positions[2][0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_point_cloud() {
        let cloud = PointCloud::from_positions(vec![[1.0, 0.0, 0.0]]);
        // 90-degree rotation around z-axis
        let rot = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let trans = [0.0, 0.0, 0.0];
        let result = transform_point_cloud(&cloud, rot, trans);
        assert!((result.positions[0][0]).abs() < 1e-10);
        assert!((result.positions[0][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_scale_point_cloud() {
        let cloud = PointCloud::from_positions(vec![[1.0, 2.0, 3.0]]);
        let scaled = scale_point_cloud(&cloud, 2.0);
        assert!((scaled.positions[0][0] - 2.0).abs() < 1e-10);
        assert!((scaled.positions[0][1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_translate_point_cloud() {
        let cloud = PointCloud::from_positions(vec![[1.0, 2.0, 3.0]]);
        let translated = translate_point_cloud(&cloud, [10.0, 20.0, 30.0]);
        assert!((translated.positions[0][0] - 11.0).abs() < 1e-10);
    }

    #[test]
    fn test_center_point_cloud() {
        let cloud = PointCloud::from_positions(vec![[2.0, 4.0, 6.0], [4.0, 6.0, 8.0]]);
        let centered = center_point_cloud(&cloud);
        let c = compute_centroid(&centered);
        assert!(c[0].abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
        assert!(c[2].abs() < 1e-10);
    }

    #[test]
    fn test_write_parse_ply_ascii() {
        let mut cloud = PointCloud::from_positions(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        cloud.normals = Some(vec![[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]);

        let mut buf = Vec::new();
        write_ply_ascii(&cloud, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let parsed = parse_ply_ascii(&text).unwrap();
        assert_eq!(parsed.len(), 2);
        assert!((parsed.positions[0][0] - 1.0).abs() < 0.001);
        assert!(parsed.has_normals());
    }

    #[test]
    fn test_write_ply_binary_le() {
        let cloud = PointCloud::from_positions(vec![[1.0, 2.0, 3.0]]);
        let buf = write_ply_binary_le(&cloud);
        // Check header starts with "ply"
        assert!(buf.starts_with(b"ply"));
        assert!(buf.len() > 12); // header + 3 floats
    }

    #[test]
    fn test_write_parse_pcd_ascii() {
        let cloud =
            PointCloud::from_positions(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let mut buf = Vec::new();
        write_pcd_ascii(&cloud, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let parsed = parse_pcd_ascii(&text).unwrap();
        assert_eq!(parsed.len(), 3);
        assert!((parsed.positions[1][0] - 4.0).abs() < 0.001);
    }

    #[test]
    fn test_write_parse_xyz() {
        let mut cloud = PointCloud::from_positions(vec![[1.0, 2.0, 3.0]]);
        cloud.normals = Some(vec![[0.0, 0.0, 1.0]]);
        cloud.colors = Some(vec![[255, 128, 0]]);

        let mut buf = Vec::new();
        write_xyz(&cloud, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        let parsed = parse_xyz(&text);
        assert_eq!(parsed.len(), 1);
        assert!((parsed.positions[0][0] - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_parse_xyz_basic() {
        let data = "1.0 2.0 3.0\n4.0 5.0 6.0\n";
        let cloud = parse_xyz(data);
        assert_eq!(cloud.len(), 2);
        assert!(!cloud.has_normals());
    }

    #[test]
    fn test_las_header_mock() {
        let header = LasHeader::mock(50000);
        assert_eq!(header.num_points, 50000);
        assert_eq!(header.file_signature, "LASF");
        assert!(header.bounding_volume() > 0.0);
    }

    #[test]
    fn test_e57_structure() {
        let e57 = E57File::mock();
        assert_eq!(e57.num_scans(), 1);
        assert_eq!(e57.total_points(), 1000);
    }

    #[test]
    fn test_ball_pivoting_mesh() {
        // Small planar cloud
        let mut cloud = PointCloud::new();
        for i in 0..5 {
            for j in 0..5 {
                cloud.positions.push([i as f64 * 0.5, j as f64 * 0.5, 0.0]);
            }
        }
        let triangles = ball_pivoting_mesh(&cloud, 0.6);
        assert!(!triangles.is_empty(), "Should produce some triangles");
        for tri in &triangles {
            assert!(tri.indices[0] < cloud.len());
            assert!(tri.indices[1] < cloud.len());
            assert!(tri.indices[2] < cloud.len());
        }
    }

    #[test]
    fn test_dist_sq_and_dist() {
        let a = [1.0, 0.0, 0.0];
        let b = [4.0, 0.0, 0.0];
        assert!((dist_sq(a, b) - 9.0).abs() < 1e-10);
        assert!((dist(a, b) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize3(v);
        assert!((n[0] - 0.6).abs() < 1e-10);
        assert!((n[1] - 0.8).abs() < 1e-10);

        let z = normalize3([0.0, 0.0, 0.0]);
        assert!(z[0].abs() < 1e-10);
    }

    #[test]
    fn test_cross3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kdtree_large() {
        let mut cloud = PointCloud::new();
        for i in 0..100 {
            cloud
                .positions
                .push([(i % 10) as f64, (i / 10) as f64, 0.0]);
        }
        let tree = KdTree::build(&cloud);
        assert_eq!(tree.len(), 100);

        let (idx, _d) = tree.nearest([4.5, 4.5, 0.0]).unwrap();
        // Should be one of the points near (4.5, 4.5)
        let p = cloud.positions[idx];
        assert!((p[0] - 4.5).abs() <= 1.0);
        assert!((p[1] - 4.5).abs() <= 1.0);
    }

    #[test]
    fn test_point_cloud_reserve() {
        let mut cloud = PointCloud::new();
        cloud.reserve(100);
        assert!(cloud.is_empty());
        cloud.add_point([1.0, 2.0, 3.0]);
        assert_eq!(cloud.len(), 1);
    }

    #[test]
    fn test_bounding_box_empty() {
        let cloud = PointCloud::new();
        assert!(BoundingBox::from_point_cloud(&cloud).is_none());
    }
}
