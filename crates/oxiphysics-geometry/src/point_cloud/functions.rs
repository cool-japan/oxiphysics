//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    IcpRegistration, KdTree3D, NormalEstimation, PointCloud, PointCloudFilter, RansacPlaneResult,
};

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-15 {
        [0.0, 0.0, 1.0]
    } else {
        scale(a, 1.0 / len)
    }
}
pub(super) fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}
pub(super) fn mat_vec_3(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
pub(super) fn identity_3x3() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
pub(super) fn mat_transpose_3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}
pub(super) fn mat_mul_3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
pub(super) fn mat_det_3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
pub(super) fn mean_points(pts: &[[f64; 3]]) -> [f64; 3] {
    if pts.is_empty() {
        return [0.0; 3];
    }
    let n = pts.len() as f64;
    let mut s = [0.0f64; 3];
    for &p in pts {
        for i in 0..3 {
            s[i] += p[i];
        }
    }
    [s[0] / n, s[1] / n, s[2] / n]
}
/// One-sided Jacobi SVD for 3×3 matrices.
/// Returns (U, singular_values, V) such that A = U * diag(s) * V^T.
pub(super) fn svd_3x3_sym_jacobi(a: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3]) {
    let ata = mat_mul_3(&mat_transpose_3(&a), &a);
    let (v, eigenvalues) = jacobi_eigen_3x3(ata);
    let s = [
        eigenvalues[0].max(0.0).sqrt(),
        eigenvalues[1].max(0.0).sqrt(),
        eigenvalues[2].max(0.0).sqrt(),
    ];
    let mut u = [[0.0f64; 3]; 3];
    for j in 0..3 {
        let vj = [v[0][j], v[1][j], v[2][j]];
        let avj = mat_vec_3(&a, vj);
        let col = if s[j] > 1e-12 {
            scale(avj, 1.0 / s[j])
        } else {
            gram_schmidt_fallback(&u, j)
        };
        u[0][j] = col[0];
        u[1][j] = col[1];
        u[2][j] = col[2];
    }
    (u, s, v)
}
pub(super) fn gram_schmidt_fallback(u: &[[f64; 3]; 3], col: usize) -> [f64; 3] {
    let candidates: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for &cand in &candidates {
        let mut v = cand;
        for j in 0..col {
            let uj = [u[0][j], u[1][j], u[2][j]];
            let proj = dot(v, uj);
            v = sub(v, scale(uj, proj));
        }
        let len = length(v);
        if len > 1e-12 {
            return scale(v, 1.0 / len);
        }
    }
    [0.0, 0.0, 1.0]
}
/// Jacobi eigendecomposition of a symmetric 3×3 matrix.
/// Returns (eigenvectors_col_matrix, eigenvalues).
pub(super) fn jacobi_eigen_3x3(mut a: [[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3]) {
    let mut v = identity_3x3();
    for _ in 0..50 {
        let mut max_val = 0.0f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }
        let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            -1.0 / (-theta + (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for r in 0..3 {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];
            }
        }
        for r in 0..3 {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = c * vrp - s * vrq;
            v[r][q] = s * vrp + c * vrq;
        }
    }
    let eigenvalues = [a[0][0], a[1][1], a[2][2]];
    (v, eigenvalues)
}
/// Thin wrapper that returns only the eigenvalues from `jacobi_eigen_3x3`.
pub(super) fn jacobi_eigenvalues_3x3(a: [[f64; 3]; 3]) -> [f64; 3] {
    let (_vecs, vals) = jacobi_eigen_3x3(a);
    vals
}
/// Estimate per-point normals for `cloud` using k-nearest-neighbor PCA.
///
/// Delegates to [`NormalEstimation`] with the given neighbor count.
pub fn estimate_normals(cloud: &PointCloud, k_neighbors: usize) -> Vec<[f64; 3]> {
    NormalEstimation::new(k_neighbors).estimate(cloud)
}
/// Downsample `cloud` by averaging points in each voxel of side `voxel_size`.
///
/// Delegates to [`PointCloudFilter::voxel_downsample`].
pub fn voxel_downsample(cloud: &PointCloud, voxel_size: f64) -> PointCloud {
    PointCloudFilter::voxel_downsample(cloud, voxel_size)
}
/// Remove statistical outliers from `cloud`.
///
/// Points whose mean distance to their `k` nearest neighbors exceeds
/// `mean + std_ratio * std` are removed.
/// Delegates to [`PointCloudFilter::statistical_outlier_removal`].
pub fn statistical_outlier_removal(cloud: &PointCloud, k: usize, std_ratio: f64) -> PointCloud {
    PointCloudFilter::statistical_outlier_removal(cloud, k, std_ratio)
}
/// Compute the axis-aligned bounding box of `cloud`.
///
/// Returns `(min_corner, max_corner)` as `[f64;3]` arrays.
/// Returns `([0;3], [0;3])` for an empty cloud.
pub fn compute_bounding_box(cloud: &PointCloud) -> ([f64; 3], [f64; 3]) {
    cloud.bounding_box()
}
/// Align `source` onto `target` using point-to-point ICP.
///
/// Returns `(transform_4x4_row_major, final_rms_error)`.
/// The transform encodes rotation R and translation t as a 4×4 homogeneous
/// matrix stored row-major: `[R | t; 0 0 0 1]`.
pub fn icp_align(source: &PointCloud, target: &PointCloud, max_iter: usize) -> ([f64; 16], f64) {
    let icp = IcpRegistration::new(max_iter, 1e-8);
    let result = icp.run(source, target);
    let r = result.rotation;
    let t = result.translation;
    let m = [
        r[0][0], r[0][1], r[0][2], t[0], r[1][0], r[1][1], r[1][2], t[1], r[2][0], r[2][1],
        r[2][2], t[2], 0.0, 0.0, 0.0, 1.0,
    ];
    (m, result.final_error)
}
/// Point-to-point ICP alignment returning a 4×3 matrix (rotation 3 rows + translation as last row).
///
/// Layout: rows 0-2 are the 3×3 rotation matrix, row 3 is the translation vector.
pub fn icp_point_to_point(
    source: &[[f64; 3]],
    target: &[[f64; 3]],
    max_iter: usize,
) -> [[f64; 3]; 4] {
    let src_cloud = PointCloud::from_points(source.to_vec());
    let tgt_cloud = PointCloud::from_points(target.to_vec());
    let icp = IcpRegistration::new(max_iter, 1e-8);
    let result = icp.run(&src_cloud, &tgt_cloud);
    let r = result.rotation;
    [r[0], r[1], r[2], result.translation]
}
/// Compute per-point normals for a point slice using PCA of k nearest neighbors.
pub fn compute_point_cloud_normals(points: &[[f64; 3]], k: usize) -> Vec<[f64; 3]> {
    let cloud = PointCloud::from_points(points.to_vec());
    NormalEstimation::new(k).estimate(&cloud)
}
/// Simplified FPFH-like descriptor for point at `idx` within radius `r`.
///
/// Computes a 33-bin histogram of angular features between the point's normal
/// and its neighbors, returning a `Vec`f64` of length 33.
pub fn fpfh_feature(points: &[[f64; 3]], normals: &[[f64; 3]], idx: usize, r: f64) -> Vec<f64> {
    let tree = KdTree3D::build(points);
    let neighbors = tree.range_search(points[idx], r);
    let n_bins = 33usize;
    let mut hist = vec![0.0f64; n_bins];
    let ni = normals[idx];
    let pi = points[idx];
    for &j in &neighbors {
        if j == idx {
            continue;
        }
        let nj = normals[j];
        let pj = points[j];
        let d = sub(pj, pi);
        let d_len = length(d);
        if d_len < 1e-12 {
            continue;
        }
        let cos_alpha = dot(ni, nj).clamp(-1.0, 1.0);
        let d_norm = scale(d, 1.0 / d_len);
        let cos_phi = dot(ni, d_norm).clamp(-1.0, 1.0);
        let bin_a = ((cos_alpha + 1.0) / 2.0 * (n_bins / 3) as f64) as usize;
        let bin_b = ((cos_phi + 1.0) / 2.0 * (n_bins / 3) as f64) as usize;
        let bin_c = (dot(nj, d_norm).clamp(-1.0, 1.0) + 1.0) / 2.0 * (n_bins / 3) as f64;
        let bin_c = bin_c as usize;
        hist[(bin_a).min(n_bins / 3 - 1)] += 1.0;
        hist[(n_bins / 3 + bin_b).min(2 * n_bins / 3 - 1)] += 1.0;
        hist[(2 * n_bins / 3 + bin_c).min(n_bins - 1)] += 1.0;
    }
    let total: f64 = hist.iter().sum();
    if total > 0.0 {
        for v in &mut hist {
            *v /= total;
        }
    }
    hist
}
/// Selects `k` points from a cloud using Farthest Point Sampling (FPS).
///
/// Starts from the point at index `seed_idx`, then iteratively picks the
/// point farthest from the already-selected set (measured as minimum squared
/// distance to any selected point).
///
/// Returns indices into `points`.
pub fn farthest_point_sampling(points: &[[f64; 3]], k: usize, seed_idx: usize) -> Vec<usize> {
    let n = points.len();
    if n == 0 || k == 0 {
        return Vec::new();
    }
    let k = k.min(n);
    let mut selected = Vec::with_capacity(k);
    let mut min_dist = vec![f64::MAX; n];
    let seed = seed_idx % n;
    selected.push(seed);
    for i in 0..n {
        min_dist[i] = dist2_pts(points[i], points[seed]);
    }
    for _ in 1..k {
        let next = min_dist
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        selected.push(next);
        for i in 0..n {
            let d = dist2_pts(points[i], points[next]);
            if d < min_dist[i] {
                min_dist[i] = d;
            }
        }
    }
    selected
}
pub(super) fn dist2_pts(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}
/// Fit a plane to a point cloud using RANSAC.
///
/// Randomly samples 3 points per iteration, computes the plane normal, and
/// counts inliers within `threshold` distance to the plane.
///
/// # Parameters
/// - `points`: input point cloud.
/// - `n_iter`: number of RANSAC iterations (use ≥ 100 for noisy data).
/// - `threshold`: maximum point-to-plane distance to count as inlier.
///
/// Returns `None` if the cloud has fewer than 3 points.
pub fn ransac_fit_plane(
    points: &[[f64; 3]],
    n_iter: usize,
    threshold: f64,
) -> Option<RansacPlaneResult> {
    let n = points.len();
    if n < 3 {
        return None;
    }
    let mut best_inliers: Vec<usize> = Vec::new();
    let mut best_normal = [0.0f64; 3];
    let mut best_d = 0.0_f64;
    let step1 = 1usize;
    let step2 = n / 3 + 1;
    let step3 = 2 * n / 3 + 1;
    for iter in 0..n_iter {
        let i0 = (iter * step1) % n;
        let i1 = (iter * step2 + 1) % n;
        let i2 = (iter * step3 + 2) % n;
        if i0 == i1 || i1 == i2 || i0 == i2 {
            continue;
        }
        let p0 = points[i0];
        let p1 = points[i1];
        let p2 = points[i2];
        let e1 = sub(p1, p0);
        let e2 = sub(p2, p0);
        let normal_raw = cross(e1, e2);
        let normal = normalize(normal_raw);
        if length(normal) < 0.5 {
            continue;
        }
        let d = dot(normal, p0);
        let inliers: Vec<usize> = points
            .iter()
            .enumerate()
            .filter(|(_, p)| (dot(normal, **p) - d).abs() <= threshold)
            .map(|(i, _)| i)
            .collect();
        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_normal = normal;
            best_d = d;
        }
    }
    if best_inliers.is_empty() {
        return None;
    }
    let n_in = best_inliers.len() as f64;
    let mut centroid = [0.0f64; 3];
    for &i in &best_inliers {
        for k in 0..3 {
            centroid[k] += points[i][k];
        }
    }
    for k in 0..3 {
        centroid[k] /= n_in;
    }
    let n_inliers = best_inliers.len();
    let point_on_plane = sub(
        centroid,
        scale(best_normal, dot(best_normal, centroid) - best_d),
    );
    Some(RansacPlaneResult {
        normal: best_normal,
        point_on_plane,
        inliers: best_inliers,
        n_inliers,
    })
}
/// Axis-aligned bounding box extent: `(dimensions, center)`.
///
/// `dimensions\[i\]` is the length along axis i; `center` is the AABB center.
pub fn aabb_extent(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    if points.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = points[0];
    let mut mx = points[0];
    for &p in points.iter().skip(1) {
        for i in 0..3 {
            if p[i] < mn[i] {
                mn[i] = p[i];
            }
            if p[i] > mx[i] {
                mx[i] = p[i];
            }
        }
    }
    let dims = [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]];
    let center = [
        (mn[0] + mx[0]) * 0.5,
        (mn[1] + mx[1]) * 0.5,
        (mn[2] + mx[2]) * 0.5,
    ];
    (dims, center)
}
/// PCA-based Oriented Bounding Box (OBB) approximation.
///
/// Returns `(axes, half_extents, center)` where `axes` is a 3×3 column
/// matrix of the principal directions (sorted by decreasing variance).
pub fn pca_obb(points: &[[f64; 3]]) -> ([[f64; 3]; 3], [f64; 3], [f64; 3]) {
    if points.is_empty() {
        return (
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0; 3],
            [0.0; 3],
        );
    }
    let n = points.len() as f64;
    let mut centroid = [0.0f64; 3];
    for &p in points {
        for i in 0..3 {
            centroid[i] += p[i];
        }
    }
    for i in 0..3 {
        centroid[i] /= n;
    }
    let mut cov = [[0.0f64; 3]; 3];
    for &p in points {
        let d = sub(p, centroid);
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    let axis0 = power_iteration_largest(&cov);
    let deflated = deflate_matrix(&cov, &axis0);
    let axis1 = power_iteration_largest(&deflated);
    let axis2 = normalize(cross(axis0, axis1));
    let axes = [axis0, axis1, axis2];
    let mut min_proj = [f64::MAX; 3];
    let mut max_proj = [f64::MIN; 3];
    for &p in points {
        let d = sub(p, centroid);
        for k in 0..3 {
            let proj = dot(d, axes[k]);
            if proj < min_proj[k] {
                min_proj[k] = proj;
            }
            if proj > max_proj[k] {
                max_proj[k] = proj;
            }
        }
    }
    let half_extents = [
        (max_proj[0] - min_proj[0]) * 0.5,
        (max_proj[1] - min_proj[1]) * 0.5,
        (max_proj[2] - min_proj[2]) * 0.5,
    ];
    (axes, half_extents, centroid)
}
pub(super) fn power_iteration_largest(m: &[[f64; 3]; 3]) -> [f64; 3] {
    let mut v = [1.0f64, 0.0, 0.0];
    for _ in 0..30 {
        let mv = mat_vec_3_pca(m, v);
        let len = length(mv);
        if len < 1e-15 {
            break;
        }
        v = scale(mv, 1.0 / len);
    }
    normalize(v)
}
pub(super) fn deflate_matrix(m: &[[f64; 3]; 3], axis: &[f64; 3]) -> [[f64; 3]; 3] {
    let mv = mat_vec_3_pca(m, *axis);
    let lambda = dot(mv, *axis);
    let mut result = *m;
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] -= lambda * axis[i] * axis[j];
        }
    }
    result
}
pub(super) fn mat_vec_3_pca(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
