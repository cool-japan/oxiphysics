//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::HashMap;

/// Static utility methods for filtering point clouds.
pub struct PointCloudFilter;
impl PointCloudFilter {
    /// Downsamples by averaging points that fall into the same voxel.
    pub fn voxel_downsample(cloud: &PointCloud, voxel_size: f64) -> PointCloud {
        if cloud.points.is_empty() || voxel_size <= 0.0 {
            return PointCloud::new();
        }
        let inv = 1.0 / voxel_size;
        let mut voxels: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for (i, &p) in cloud.points.iter().enumerate() {
            let key = (
                (p[0] * inv).floor() as i64,
                (p[1] * inv).floor() as i64,
                (p[2] * inv).floor() as i64,
            );
            voxels.entry(key).or_default().push(i);
        }
        let mut out = PointCloud::new();
        for indices in voxels.values() {
            let n = indices.len() as f64;
            let mut sum = [0.0f64; 3];
            for &i in indices {
                let p = cloud.points[i];
                for k in 0..3 {
                    sum[k] += p[k];
                }
            }
            out.add_point([sum[0] / n, sum[1] / n, sum[2] / n]);
        }
        out
    }
    /// Removes points whose mean distance to their k nearest neighbors
    /// exceeds `mean + std_factor * std`.
    pub fn statistical_outlier_removal(
        cloud: &PointCloud,
        k: usize,
        std_factor: f64,
    ) -> PointCloud {
        if cloud.points.is_empty() {
            return PointCloud::new();
        }
        let tree = KdTree3D::build(&cloud.points);
        let mean_dists: Vec<f64> = cloud
            .points
            .iter()
            .map(|&p| {
                let neighbors = tree.k_nearest(p, k + 1);
                let dists: Vec<f64> = neighbors
                    .iter()
                    .filter(|&&(idx, _)| cloud.points[idx] != p || idx != 0)
                    .map(|&(_, d)| d.sqrt())
                    .collect();
                if dists.is_empty() {
                    0.0
                } else {
                    dists.iter().sum::<f64>() / dists.len() as f64
                }
            })
            .collect();
        let mean = mean_dists.iter().sum::<f64>() / mean_dists.len() as f64;
        let variance = mean_dists
            .iter()
            .map(|&d| (d - mean) * (d - mean))
            .sum::<f64>()
            / mean_dists.len() as f64;
        let std = variance.sqrt();
        let threshold = mean + std_factor * std;
        let mut out = PointCloud::new();
        for (i, &md) in mean_dists.iter().enumerate() {
            if md <= threshold {
                out.add_point(cloud.points[i]);
                if !cloud.normals.is_empty() {
                    out.normals.push(cloud.normals[i]);
                }
                if !cloud.colors.is_empty() {
                    out.colors.push(cloud.colors[i]);
                }
            }
        }
        out
    }
    /// Removes points that have fewer than `min_neighbors` within `radius`.
    pub fn radius_outlier_removal(
        cloud: &PointCloud,
        radius: f64,
        min_neighbors: usize,
    ) -> PointCloud {
        if cloud.points.is_empty() {
            return PointCloud::new();
        }
        let tree = KdTree3D::build(&cloud.points);
        let mut out = PointCloud::new();
        for (i, &p) in cloud.points.iter().enumerate() {
            let count = tree.range_search(p, radius).len();
            if count.saturating_sub(1) >= min_neighbors {
                out.add_point(p);
                if !cloud.normals.is_empty() {
                    out.normals.push(cloud.normals[i]);
                }
                if !cloud.colors.is_empty() {
                    out.colors.push(cloud.colors[i]);
                }
            }
        }
        out
    }
}
/// A 3-D k-d tree for efficient nearest-neighbor queries.
pub struct KdTree3D {
    pub(super) points: Vec<[f64; 3]>,
    pub(super) root: Option<Box<KdNode3D>>,
}
impl KdTree3D {
    /// Builds a k-d tree from a slice of points.
    pub fn build(points: &[[f64; 3]]) -> Self {
        let mut indices: Vec<usize> = (0..points.len()).collect();
        let root = Self::build_recursive(points, &mut indices, 0);
        Self {
            points: points.to_vec(),
            root,
        }
    }
    fn build_recursive(
        points: &[[f64; 3]],
        indices: &mut [usize],
        depth: usize,
    ) -> Option<Box<KdNode3D>> {
        if indices.is_empty() {
            return None;
        }
        let axis = depth % 3;
        indices.sort_unstable_by(|&a, &b| {
            points[a][axis]
                .partial_cmp(&points[b][axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = indices.len() / 2;
        let idx = indices[mid];
        let left = Self::build_recursive(points, &mut indices[..mid], depth + 1);
        let right = Self::build_recursive(points, &mut indices[mid + 1..], depth + 1);
        Some(Box::new(KdNode3D {
            point_idx: idx,
            split_axis: axis,
            left,
            right,
        }))
    }
    /// Returns `(index, squared_distance)` for the closest point to `query`.
    pub fn nearest_neighbor(&self, query: [f64; 3]) -> (usize, f64) {
        let mut best = (0, f64::MAX);
        if let Some(ref root) = self.root {
            Self::nn_search(&self.points, root, query, &mut best);
        }
        best
    }
    fn nn_search(points: &[[f64; 3]], node: &KdNode3D, query: [f64; 3], best: &mut (usize, f64)) {
        let p = points[node.point_idx];
        let d2 = dist2(p, query);
        if d2 < best.1 {
            *best = (node.point_idx, d2);
        }
        let axis = node.split_axis;
        let diff = query[axis] - p[axis];
        let (near, far) = if diff <= 0.0 {
            (node.left.as_deref(), node.right.as_deref())
        } else {
            (node.right.as_deref(), node.left.as_deref())
        };
        if let Some(near_node) = near {
            Self::nn_search(points, near_node, query, best);
        }
        if diff * diff < best.1
            && let Some(far_node) = far
        {
            Self::nn_search(points, far_node, query, best);
        }
    }
    /// Returns up to `k` nearest neighbors, sorted by ascending squared distance.
    pub fn k_nearest(&self, query: [f64; 3], k: usize) -> Vec<(usize, f64)> {
        if k == 0 || self.points.is_empty() {
            return Vec::new();
        }
        let mut heap: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        if let Some(ref root) = self.root {
            Self::knn_search(&self.points, root, query, k, &mut heap);
        }
        heap.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        heap
    }
    fn knn_search(
        points: &[[f64; 3]],
        node: &KdNode3D,
        query: [f64; 3],
        k: usize,
        heap: &mut Vec<(usize, f64)>,
    ) {
        let p = points[node.point_idx];
        let d2 = dist2(p, query);
        heap.push((node.point_idx, d2));
        if heap.len() > k {
            let worst_pos = heap
                .iter()
                .enumerate()
                .max_by(|a, b| {
                    a.1.1
                        .partial_cmp(&b.1.1)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            heap.swap_remove(worst_pos);
        }
        let worst_dist2 = heap.iter().map(|e| e.1).fold(0.0_f64, f64::max);
        let axis = node.split_axis;
        let diff = query[axis] - p[axis];
        let (near, far) = if diff <= 0.0 {
            (node.left.as_deref(), node.right.as_deref())
        } else {
            (node.right.as_deref(), node.left.as_deref())
        };
        if let Some(near_node) = near {
            Self::knn_search(points, near_node, query, k, heap);
        }
        let current_worst = heap.iter().map(|e| e.1).fold(0.0_f64, f64::max);
        if (heap.len() < k || diff * diff < current_worst)
            && let Some(far_node) = far
        {
            Self::knn_search(points, far_node, query, k, heap);
        }
        let _ = worst_dist2;
    }
    /// Returns indices of all points within `radius` of `center`.
    pub fn range_search(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        let r2 = radius * radius;
        if let Some(ref root) = self.root {
            Self::range_recursive(&self.points, root, center, r2, &mut result);
        }
        result
    }
    fn range_recursive(
        points: &[[f64; 3]],
        node: &KdNode3D,
        center: [f64; 3],
        r2: f64,
        result: &mut Vec<usize>,
    ) {
        let p = points[node.point_idx];
        if dist2(p, center) <= r2 {
            result.push(node.point_idx);
        }
        let axis = node.split_axis;
        let diff = center[axis] - p[axis];
        if (diff * diff <= r2 || diff <= 0.0)
            && let Some(ref left) = node.left
        {
            Self::range_recursive(points, left, center, r2, result);
        }
        if (diff * diff <= r2 || diff > 0.0)
            && let Some(ref right) = node.right
        {
            Self::range_recursive(points, right, center, r2, result);
        }
    }
}
/// Internal node for the 3-D k-d tree.
pub struct KdNode3D {
    /// Index into the original points slice.
    pub point_idx: usize,
    /// Axis on which this node splits (0 = X, 1 = Y, 2 = Z).
    pub split_axis: usize,
    /// Left subtree.
    pub left: Option<Box<KdNode3D>>,
    /// Right subtree.
    pub right: Option<Box<KdNode3D>>,
}
/// Estimates per-point normals via PCA of k-nearest-neighbor covariance.
pub struct NormalEstimation {
    /// Number of nearest neighbors to use for the local PCA.
    pub k_neighbors: usize,
}
impl NormalEstimation {
    /// Creates a new estimator with the given neighbor count.
    pub fn new(k_neighbors: usize) -> Self {
        Self { k_neighbors }
    }
    /// Estimates normals for every point in `cloud`.
    pub fn estimate(&self, cloud: &PointCloud) -> Vec<[f64; 3]> {
        let tree = KdTree3D::build(&cloud.points);
        cloud
            .points
            .iter()
            .map(|&p| {
                let neighbors = tree.k_nearest(p, self.k_neighbors);
                let pts: Vec<[f64; 3]> = neighbors
                    .iter()
                    .map(|&(idx, _)| cloud.points[idx])
                    .collect();
                if pts.len() < 3 {
                    return [0.0, 0.0, 1.0];
                }
                let cov = Self::covariance_matrix(&pts);
                Self::smallest_eigenvector_3x3_sym(&cov)
            })
            .collect()
    }
    /// Computes the 3×3 covariance matrix of the given points.
    pub fn covariance_matrix(pts: &[[f64; 3]]) -> [[f64; 3]; 3] {
        if pts.is_empty() {
            return [[0.0; 3]; 3];
        }
        let n = pts.len() as f64;
        let mut centroid = [0.0f64; 3];
        for &p in pts {
            for (c, p_i) in centroid.iter_mut().zip(p.iter()) {
                *c += p_i;
            }
        }
        for c in centroid.iter_mut() {
            *c /= n;
        }
        let mut cov = [[0.0f64; 3]; 3];
        for &p in pts {
            let d = sub(p, centroid);
            for (i, cov_row) in cov.iter_mut().enumerate() {
                for (j, cov_ij) in cov_row.iter_mut().enumerate() {
                    *cov_ij += d[i] * d[j];
                }
            }
        }
        cov
    }
    /// Returns the eigenvector corresponding to the smallest eigenvalue of a
    /// symmetric 3×3 matrix using shifted power iteration.
    pub fn smallest_eigenvector_3x3_sym(m: &[[f64; 3]; 3]) -> [f64; 3] {
        let mut v = [1.0f64, 1.0, 1.0];
        for _ in 0..20 {
            v = mat_vec_3(m, v);
            let len = length(v);
            if len < 1e-15 {
                break;
            }
            v = scale(v, 1.0 / len);
        }
        let lambda_max = dot(mat_vec_3(m, v), v);
        let shifted = [
            [lambda_max - m[0][0], -m[0][1], -m[0][2]],
            [-m[1][0], lambda_max - m[1][1], -m[1][2]],
            [-m[2][0], -m[2][1], lambda_max - m[2][2]],
        ];
        let mut u = [1.0f64, 0.0, 0.0];
        for _ in 0..50 {
            u = mat_vec_3(&shifted, u);
            let len = length(u);
            if len < 1e-15 {
                break;
            }
            u = scale(u, 1.0 / len);
        }
        normalize(u)
    }
}
/// Result of a RANSAC plane fitting.
pub struct RansacPlaneResult {
    /// Plane normal (unit vector).
    pub normal: [f64; 3],
    /// A point on the plane (centroid of inliers).
    pub point_on_plane: [f64; 3],
    /// Indices of inlier points.
    pub inliers: Vec<usize>,
    /// Number of inliers.
    pub n_inliers: usize,
}
/// Result of an ICP alignment.
pub struct IcpResult {
    /// 3×3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation vector.
    pub translation: [f64; 3],
    /// Mean squared distance after convergence.
    pub final_error: f64,
    /// Number of iterations performed.
    pub iterations: usize,
}
impl IcpResult {
    /// Applies the stored transform to `cloud` and returns the result.
    pub fn apply_to(&self, cloud: &PointCloud) -> PointCloud {
        let mut out = PointCloud::new();
        for &p in &cloud.points {
            let rp = mat_vec_3(&self.rotation, p);
            out.add_point(add(rp, self.translation));
        }
        if !cloud.normals.is_empty() {
            for &n in &cloud.normals {
                out.normals.push(mat_vec_3(&self.rotation, n));
            }
        }
        out.colors = cloud.colors.clone();
        out
    }
}
/// A collection of 3-D points with optional per-point normals and colors.
pub struct PointCloud {
    /// XYZ positions.
    pub points: Vec<[f64; 3]>,
    /// Per-point surface normals (may be empty).
    pub normals: Vec<[f64; 3]>,
    /// Per-point RGB colors in \[0, 1\] (may be empty).
    pub colors: Vec<[f32; 3]>,
}
impl PointCloud {
    /// Creates an empty point cloud.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            normals: Vec::new(),
            colors: Vec::new(),
        }
    }
    /// Appends a point.
    pub fn add_point(&mut self, p: [f64; 3]) {
        self.points.push(p);
    }
    /// Returns the number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Returns `true` if the cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Returns the axis-aligned bounding box as `(min, max)`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.points.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = self.points[0];
        let mut mx = self.points[0];
        for &p in &self.points {
            for i in 0..3 {
                if p[i] < mn[i] {
                    mn[i] = p[i];
                }
                if p[i] > mx[i] {
                    mx[i] = p[i];
                }
            }
        }
        (mn, mx)
    }
    /// Returns the centroid (mean position) of all points.
    pub fn centroid(&self) -> [f64; 3] {
        if self.points.is_empty() {
            return [0.0; 3];
        }
        let n = self.points.len() as f64;
        let mut sum = [0.0f64; 3];
        for &p in &self.points {
            for i in 0..3 {
                sum[i] += p[i];
            }
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
    /// Translates all points by `offset`.
    pub fn translate(&mut self, offset: [f64; 3]) {
        for p in &mut self.points {
            *p = add(*p, offset);
        }
    }
    /// Uniformly scales all points around the origin.
    pub fn scale_uniform(&mut self, s: f64) {
        for p in &mut self.points {
            *p = scale(*p, s);
        }
    }
    /// Builds a height-field point cloud from a regular `nx × ny` grid.
    ///
    /// Points are placed at `(i * dx, j * dx, height_fn(i * dx, j * dx))`.
    pub fn from_grid(nx: usize, ny: usize, height_fn: impl Fn(f64, f64) -> f64, dx: f64) -> Self {
        let mut cloud = Self::new();
        for i in 0..nx {
            for j in 0..ny {
                let x = i as f64 * dx;
                let y = j as f64 * dx;
                let z = height_fn(x, y);
                cloud.add_point([x, y, z]);
            }
        }
        cloud
    }
}
impl PointCloud {
    /// Constructs a PointCloud from a vec of points (expansion API).
    pub fn from_points(points: Vec<[f64; 3]>) -> Self {
        Self {
            points,
            normals: Vec::new(),
            colors: Vec::new(),
        }
    }
    /// Estimate per-point normals using PCA of the `k` nearest neighbors.
    /// Populates `self.normals`.
    pub fn estimate_normals_pca(&mut self, k: usize) {
        let tree = KdTree3D::build(&self.points);
        self.normals = self
            .points
            .iter()
            .map(|&p| {
                let neighbors = tree.k_nearest(p, k);
                let pts: Vec<[f64; 3]> =
                    neighbors.iter().map(|&(idx, _)| self.points[idx]).collect();
                if pts.len() < 3 {
                    return [0.0, 0.0, 1.0];
                }
                let cov = NormalEstimation::covariance_matrix(&pts);
                NormalEstimation::smallest_eigenvector_3x3_sym(&cov)
            })
            .collect();
    }
    /// Return a new downsampled cloud using a voxel grid of side `voxel_size`.
    pub fn voxel_downsample(&self, voxel_size: f64) -> Self {
        PointCloudFilter::voxel_downsample(self, voxel_size)
    }
    /// Return a new cloud with statistical outliers removed.
    pub fn statistical_outlier_removal(&self, k: usize, std_ratio: f64) -> Self {
        PointCloudFilter::statistical_outlier_removal(self, k, std_ratio)
    }
}
impl PointCloud {
    /// Farthest Point Sampling: select `k` well-spread points from this cloud.
    ///
    /// Returns a new PointCloud with the selected points (normals / colors
    /// are copied if present).
    pub fn farthest_point_sample(&self, k: usize) -> Self {
        let indices = farthest_point_sampling(&self.points, k, 0);
        let mut out = PointCloud::new();
        for &i in &indices {
            out.add_point(self.points[i]);
            if !self.normals.is_empty() {
                out.normals.push(self.normals[i]);
            }
            if !self.colors.is_empty() {
                out.colors.push(self.colors[i]);
            }
        }
        out
    }
    /// Fit a plane to the point cloud using RANSAC.
    ///
    /// Returns `None` if fewer than 3 points are present.
    pub fn fit_plane_ransac(&self, n_iter: usize, threshold: f64) -> Option<RansacPlaneResult> {
        ransac_fit_plane(&self.points, n_iter, threshold)
    }
    /// Compute the PCA-based OBB of this cloud.
    ///
    /// Returns `(principal_axes, half_extents, center)`.
    pub fn pca_obb(&self) -> ([[f64; 3]; 3], [f64; 3], [f64; 3]) {
        pca_obb(&self.points)
    }
    /// Compute PCA-based principal curvatures for each point.
    ///
    /// For each point a k-NN neighbourhood is gathered, the 3×3 covariance
    /// matrix is formed, and the two smallest eigenvalues (λ₁ ≤ λ₂) are
    /// returned as the approximate principal curvatures κ₁ and κ₂.
    ///
    /// Returns a `Vec<(f64, f64)>` of `(kappa1, kappa2)` for every point.
    pub fn compute_principal_curvatures(&self, k: usize) -> Vec<(f64, f64)> {
        if self.points.is_empty() {
            return Vec::new();
        }
        let tree = KdTree3D::build(&self.points);
        self.points
            .iter()
            .map(|&p| {
                let neighbors = tree.k_nearest(p, k + 1);
                let pts: Vec<[f64; 3]> =
                    neighbors.iter().map(|&(idx, _)| self.points[idx]).collect();
                if pts.len() < 3 {
                    return (0.0, 0.0);
                }
                let n = pts.len() as f64;
                let mut centroid = [0.0f64; 3];
                for &q in &pts {
                    for (c, q_i) in centroid.iter_mut().zip(q.iter()) {
                        *c += q_i;
                    }
                }
                for c in centroid.iter_mut() {
                    *c /= n;
                }
                let mut cov = [[0.0f64; 3]; 3];
                for &q in &pts {
                    let d = sub(q, centroid);
                    for (i, cov_row) in cov.iter_mut().enumerate() {
                        for (j, cov_ij) in cov_row.iter_mut().enumerate() {
                            *cov_ij += d[i] * d[j];
                        }
                    }
                }
                for cov_row in cov.iter_mut() {
                    for cov_ij in cov_row.iter_mut() {
                        *cov_ij /= n;
                    }
                }
                let eigenvalues = jacobi_eigenvalues_3x3(cov);
                let mut ev = eigenvalues;
                ev.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                (ev[0], ev[1])
            })
            .collect()
    }
    /// Align this cloud onto `target` using ICP (point-to-point) and
    /// return the registration result (rotation + translation + error).
    ///
    /// Calls the existing `IcpRegistration::align` helper.
    pub fn icp_register(&self, target: &PointCloud) -> IcpResult {
        IcpRegistration::align(self, target)
    }
}
/// Iterative Closest Point registration.
pub struct IcpRegistration {
    /// Maximum number of ICP iterations.
    pub max_iterations: usize,
    /// Convergence threshold on mean squared distance change.
    pub convergence_tol: f64,
}
impl IcpRegistration {
    /// Creates a new ICP solver.
    pub fn new(max_iterations: usize, convergence_tol: f64) -> Self {
        Self {
            max_iterations,
            convergence_tol,
        }
    }
    /// Aligns `source` onto `target` and returns the transform.
    pub fn align(source: &PointCloud, target: &PointCloud) -> IcpResult {
        let icp = Self::new(50, 1e-6);
        icp.run(source, target)
    }
    /// Run ICP alignment: returns [`IcpResult`] after at most `self.max_iterations` steps.
    pub fn run(&self, source: &PointCloud, target: &PointCloud) -> IcpResult {
        let mut rot = identity_3x3();
        let mut trans = [0.0f64; 3];
        let mut current: Vec<[f64; 3]> = source.points.clone();
        let target_tree = KdTree3D::build(&target.points);
        let mut prev_error = f64::MAX;
        let mut iters = 0;
        for iter in 0..self.max_iterations {
            iters = iter + 1;
            let pairs: Vec<(usize, usize)> = current
                .iter()
                .enumerate()
                .map(|(i, &p)| {
                    let (ti, _) = target_tree.nearest_neighbor(p);
                    (i, ti)
                })
                .collect();
            let error: f64 = pairs
                .iter()
                .map(|&(si, ti)| dist2(current[si], target.points[ti]))
                .sum::<f64>()
                / pairs.len() as f64;
            if (prev_error - error).abs() < self.convergence_tol {
                break;
            }
            prev_error = error;
            let src_centroid =
                mean_points(&pairs.iter().map(|&(si, _)| current[si]).collect::<Vec<_>>());
            let tgt_centroid = mean_points(
                &pairs
                    .iter()
                    .map(|&(_, ti)| target.points[ti])
                    .collect::<Vec<_>>(),
            );
            let mut h = [[0.0f64; 3]; 3];
            for &(si, ti) in &pairs {
                let s = sub(current[si], src_centroid);
                let t = sub(target.points[ti], tgt_centroid);
                for i in 0..3 {
                    for j in 0..3 {
                        h[i][j] += s[i] * t[j];
                    }
                }
            }
            let (u, _s_vals, v) = svd_3x3_sym_jacobi(h);
            let det = mat_det_3(&mat_mul_3(&v, &mat_transpose_3(&u)));
            let mut d = identity_3x3();
            if det < 0.0 {
                d[2][2] = -1.0;
            }
            let r_step = mat_mul_3(&v, &mat_mul_3(&d, &mat_transpose_3(&u)));
            let t_step = sub(tgt_centroid, mat_vec_3(&r_step, src_centroid));
            for p in &mut current {
                *p = add(mat_vec_3(&r_step, *p), t_step);
            }
            rot = mat_mul_3(&r_step, &rot);
            trans = add(mat_vec_3(&r_step, trans), t_step);
        }
        let target_tree2 = KdTree3D::build(&target.points);
        let final_error = if current.is_empty() {
            0.0
        } else {
            current
                .iter()
                .map(|&p| {
                    let (ti, _) = target_tree2.nearest_neighbor(p);
                    dist2(p, target.points[ti])
                })
                .sum::<f64>()
                / current.len() as f64
        };
        IcpResult {
            rotation: rot,
            translation: trans,
            final_error,
            iterations: iters,
        }
    }
    /// Returns all closest-point pairs `(source_idx, target_idx)`.
    pub fn closest_point_pairs(source: &PointCloud, target: &PointCloud) -> Vec<(usize, usize)> {
        let tree = KdTree3D::build(&target.points);
        source
            .points
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                let (ti, _) = tree.nearest_neighbor(p);
                (i, ti)
            })
            .collect()
    }
}
