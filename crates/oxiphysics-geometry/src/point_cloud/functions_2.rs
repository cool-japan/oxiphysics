//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IcpRegistration;
    use crate::KdTree3D;
    use crate::PointCloud;
    use crate::PointCloudFilter;
    /// Helper: build a grid cloud and return it.
    fn flat_grid(nx: usize, ny: usize, dx: f64) -> PointCloud {
        PointCloud::from_grid(nx, ny, |_x, _y| 0.0, dx)
    }
    #[test]
    fn test_from_grid_point_count() {
        let cloud = flat_grid(5, 7, 1.0);
        assert_eq!(cloud.len(), 35);
    }
    #[test]
    fn test_centroid_symmetric_cloud() {
        let mut cloud = PointCloud::new();
        for &x in &[-1.0f64, 1.0] {
            for &y in &[-1.0f64, 1.0] {
                for &z in &[-1.0f64, 1.0] {
                    cloud.add_point([x, y, z]);
                }
            }
        }
        let c = cloud.centroid();
        assert!(c[0].abs() < 1e-10, "cx = {}", c[0]);
        assert!(c[1].abs() < 1e-10, "cy = {}", c[1]);
        assert!(c[2].abs() < 1e-10, "cz = {}", c[2]);
    }
    #[test]
    fn test_kdtree_nearest_single_point() {
        let pts = vec![[1.0f64, 2.0, 3.0]];
        let tree = KdTree3D::build(&pts);
        let (idx, d2) = tree.nearest_neighbor([1.0, 2.0, 3.0]);
        assert_eq!(idx, 0);
        assert!(d2 < 1e-12);
    }
    #[test]
    fn test_kdtree_k_nearest_returns_k() {
        let pts: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let tree = KdTree3D::build(&pts);
        let k = 5;
        let result = tree.k_nearest([0.0, 0.0, 0.0], k);
        assert_eq!(result.len(), k);
        assert!(result[0].1 < result[result.len() - 1].1 + 1e-12);
    }
    #[test]
    fn test_voxel_downsample_fewer_points() {
        let cloud = flat_grid(10, 10, 0.1);
        let downsampled = PointCloudFilter::voxel_downsample(&cloud, 1.0);
        assert!(
            downsampled.len() < cloud.len(),
            "expected fewer points after downsampling, got {} vs {}",
            downsampled.len(),
            cloud.len()
        );
    }
    #[test]
    fn test_bounding_box_contains_centroid() {
        let cloud = flat_grid(5, 5, 1.0);
        let (mn, mx) = cloud.bounding_box();
        let c = cloud.centroid();
        for i in 0..3 {
            assert!(
                mn[i] <= c[i] + 1e-9 && c[i] <= mx[i] + 1e-9,
                "axis {}: centroid {} not in [{}, {}]",
                i,
                c[i],
                mn[i],
                mx[i]
            );
        }
    }
    #[test]
    fn test_scale_uniform_doubles_distances() {
        let mut cloud = PointCloud::new();
        cloud.add_point([1.0, 0.0, 0.0]);
        cloud.add_point([3.0, 0.0, 0.0]);
        let d_before = dist2(cloud.points[0], cloud.points[1]).sqrt();
        cloud.scale_uniform(2.0);
        let d_after = dist2(cloud.points[0], cloud.points[1]).sqrt();
        assert!((d_after - 2.0 * d_before).abs() < 1e-10);
    }
    #[test]
    fn test_estimate_normals_sphere_surface() {
        let mut cloud = PointCloud::new();
        let n = 10usize;
        for i in 0..n {
            for j in 0..n {
                let theta = std::f64::consts::PI * i as f64 / (n - 1) as f64;
                let phi = 2.0 * std::f64::consts::PI * j as f64 / n as f64;
                let x = theta.sin() * phi.cos();
                let y = theta.sin() * phi.sin();
                let z = theta.cos();
                cloud.add_point([x, y, z]);
            }
        }
        let normals = estimate_normals(&cloud, 8);
        assert_eq!(normals.len(), cloud.len());
        let mut ok = 0usize;
        for (p, n) in cloud.points.iter().zip(normals.iter()) {
            let r = (p[0].powi(2) + p[1].powi(2) + p[2].powi(2)).sqrt();
            if r < 1e-6 {
                continue;
            }
            let dot = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / r;
            if dot.abs() > 0.5 {
                ok += 1;
            }
        }
        assert!(
            ok as f64 > cloud.len() as f64 * 0.7,
            "fewer than 70% of sphere normals are radial: {ok}/{}",
            cloud.len()
        );
    }
    #[test]
    fn test_voxel_downsample_free_fn() {
        let cloud = flat_grid(8, 8, 0.1);
        let down = voxel_downsample(&cloud, 0.5);
        assert!(down.len() < cloud.len(), "expected fewer points");
    }
    #[test]
    fn test_statistical_outlier_removal_keeps_inliers() {
        let cloud = flat_grid(5, 5, 1.0);
        let filtered = statistical_outlier_removal(&cloud, 4, 2.0);
        assert!(!filtered.is_empty(), "expected some points to remain");
    }
    #[test]
    fn test_compute_bounding_box_free_fn() {
        let cloud = flat_grid(5, 5, 1.0);
        let (mn, mx) = compute_bounding_box(&cloud);
        assert!((mn[0] - 0.0).abs() < 1e-9, "min x={}", mn[0]);
        assert!((mn[1] - 0.0).abs() < 1e-9, "min y={}", mn[1]);
        assert!((mx[0] - 4.0).abs() < 1e-9, "max x={}", mx[0]);
        assert!((mx[1] - 4.0).abs() < 1e-9, "max y={}", mx[1]);
    }
    #[test]
    fn test_icp_align_identity() {
        let cloud = flat_grid(4, 4, 1.0);
        let (m, err) = icp_align(&cloud, &cloud, 20);
        assert!(
            err < 1e-6,
            "ICP error on identity should be near 0, got {err}"
        );
        assert!((m[0] - 1.0).abs() < 1e-6, "m[0]={}", m[0]);
        assert!((m[5] - 1.0).abs() < 1e-6, "m[5]={}", m[5]);
        assert!((m[10] - 1.0).abs() < 1e-6, "m[10]={}", m[10]);
    }
    #[test]
    fn test_icp_align_translation() {
        let mut source = PointCloud::new();
        let mut target = PointCloud::new();
        for i in 0..8 {
            for j in 0..8 {
                target.add_point([i as f64, j as f64, 0.0]);
                source.add_point([i as f64 + 0.05, j as f64, 0.0]);
            }
        }
        let (_m, err) = icp_align(&source, &target, 100);
        assert!(
            err < 0.01,
            "ICP should reduce error after translation, err={err}"
        );
    }
    #[test]
    fn test_kdtree_range_search() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let tree = KdTree3D::build(&pts);
        let mut found = tree.range_search([0.0, 0.0, 0.0], 1.5);
        found.sort_unstable();
        assert_eq!(found, vec![0, 1], "range search result: {found:?}");
    }
    #[test]
    fn test_radius_outlier_removal() {
        let mut cloud = PointCloud::new();
        for i in 0..5 {
            for j in 0..5 {
                cloud.add_point([i as f64 * 0.1, j as f64 * 0.1, 0.0]);
            }
        }
        cloud.add_point([100.0, 100.0, 100.0]);
        let filtered = PointCloudFilter::radius_outlier_removal(&cloud, 0.3, 2);
        assert!(filtered.len() < cloud.len(), "outlier should be removed");
        for p in &filtered.points {
            assert!(p[0] < 50.0, "outlier survived: {:?}", p);
        }
    }
    #[test]
    fn test_normal_estimation_unit_vectors() {
        let cloud = flat_grid(5, 5, 1.0);
        let normals = estimate_normals(&cloud, 6);
        for (i, n) in normals.iter().enumerate() {
            let len = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();
            assert!((len - 1.0).abs() < 1e-10, "normal {i} not unit: len={len}");
        }
    }
    #[test]
    fn test_translate_moves_centroid() {
        let mut cloud = flat_grid(4, 4, 1.0);
        let before = cloud.centroid();
        cloud.translate([5.0, -3.0, 2.0]);
        let after = cloud.centroid();
        assert!((after[0] - (before[0] + 5.0)).abs() < 1e-9);
        assert!((after[1] - (before[1] - 3.0)).abs() < 1e-9);
        assert!((after[2] - (before[2] + 2.0)).abs() < 1e-9);
    }
    #[test]
    fn test_bounding_box_correct() {
        let pts = vec![[1.0, 2.0, 3.0], [-1.0, 0.0, 5.0], [0.0, 4.0, 1.0]];
        let cloud = PointCloud::from_points(pts);
        let (mn, mx) = cloud.bounding_box();
        assert!((mn[0] - (-1.0)).abs() < 1e-10);
        assert!((mn[1] - 0.0).abs() < 1e-10);
        assert!((mn[2] - 1.0).abs() < 1e-10);
        assert!((mx[0] - 1.0).abs() < 1e-10);
        assert!((mx[1] - 4.0).abs() < 1e-10);
        assert!((mx[2] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_centroid_expansion() {
        let pts = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let cloud = PointCloud::from_points(pts);
        let c = cloud.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - (2.0 / 3.0)).abs() < 1e-10);
        assert!(c[2].abs() < 1e-10);
    }
    #[test]
    fn test_downsample_reduces_count() {
        let cloud = flat_grid(10, 10, 0.1);
        let down = cloud.voxel_downsample(1.0);
        assert!(
            down.len() < cloud.len(),
            "expected fewer points after voxel downsampling: {} < {}",
            down.len(),
            cloud.len()
        );
    }
    #[test]
    fn test_outlier_removal_expansion() {
        let mut pts: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        pts.push([100.0, 100.0, 100.0]);
        let cloud = PointCloud::from_points(pts.clone());
        let filtered = cloud.statistical_outlier_removal(4, 1.5);
        assert!(filtered.len() < cloud.len(), "outlier should be removed");
    }
    #[test]
    fn test_icp_zero_transform_same_set() {
        let pts: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64, j as f64, 0.0]))
            .collect();
        let transform = icp_point_to_point(&pts, &pts, 20);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (transform[i][j] - expected).abs() < 1e-4,
                    "rotation[{i}][{j}] = {} expected {expected}",
                    transform[i][j]
                );
            }
        }
        for &v in &transform[3] {
            assert!(v.abs() < 1e-4, "translation should be near zero: {v}");
        }
    }
    #[test]
    fn test_compute_point_cloud_normals_unit() {
        let cloud = flat_grid(5, 5, 1.0);
        let normals = compute_point_cloud_normals(&cloud.points, 6);
        assert_eq!(normals.len(), cloud.len());
        for n in &normals {
            let len = (n[0].powi(2) + n[1].powi(2) + n[2].powi(2)).sqrt();
            assert!((len - 1.0).abs() < 1e-10, "normal not unit: len={len}");
        }
    }
    #[test]
    fn test_fpfh_feature_length_and_finite() {
        let cloud = flat_grid(5, 5, 1.0);
        let normals = compute_point_cloud_normals(&cloud.points, 6);
        let feat = fpfh_feature(&cloud.points, &normals, 12, 2.0);
        assert_eq!(feat.len(), 33, "FPFH descriptor should have 33 bins");
        for &v in &feat {
            assert!(v.is_finite(), "FPFH value must be finite: {v}");
        }
    }
    #[test]
    fn test_fps_returns_k_points() {
        let pts: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let selected = farthest_point_sampling(&pts, 5, 0);
        assert_eq!(selected.len(), 5, "FPS should return exactly k points");
    }
    #[test]
    fn test_fps_no_duplicates() {
        let pts: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let selected = farthest_point_sampling(&pts, 10, 0);
        let unique: std::collections::HashSet<usize> = selected.iter().copied().collect();
        assert_eq!(
            unique.len(),
            selected.len(),
            "FPS should return distinct indices"
        );
    }
    #[test]
    fn test_fps_valid_indices() {
        let pts: Vec<[f64; 3]> = (0..15).map(|i| [i as f64, (i as f64).sin(), 0.0]).collect();
        let selected = farthest_point_sampling(&pts, 7, 3);
        for &idx in &selected {
            assert!(
                idx < pts.len(),
                "FPS index {idx} out of range (n={})",
                pts.len()
            );
        }
    }
    #[test]
    fn test_fps_spread_for_line() {
        let pts: Vec<[f64; 3]> = (0..=10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let selected = farthest_point_sampling(&pts, 2, 0);
        assert_eq!(selected.len(), 2);
        let d = dist2_pts(pts[selected[0]], pts[selected[1]]).sqrt();
        assert!(
            d >= 9.0,
            "two FPS points on a line should be near the ends: d={d}"
        );
    }
    #[test]
    fn test_fps_point_cloud_method() {
        let cloud = flat_grid(8, 8, 1.0);
        let sampled = cloud.farthest_point_sample(10);
        assert_eq!(sampled.len(), 10);
    }
    #[test]
    fn test_fps_k_larger_than_n_clamped() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let selected = farthest_point_sampling(&pts, 100, 0);
        assert_eq!(selected.len(), 3, "FPS k > n should return n points");
    }
    #[test]
    fn test_ransac_plane_flat_grid() {
        let cloud = flat_grid(6, 6, 1.0);
        let result = cloud.fit_plane_ransac(50, 0.01);
        assert!(
            result.is_some(),
            "RANSAC should find a plane in a flat grid"
        );
        let r = result.unwrap();
        assert!(
            r.normal[2].abs() > 0.9,
            "plane normal should point ~z: {:?}",
            r.normal
        );
        assert!(
            r.n_inliers >= 36,
            "all flat-grid points should be inliers, got {}",
            r.n_inliers
        );
    }
    #[test]
    fn test_ransac_plane_normal_is_unit() {
        let pts: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64, j as f64, 0.0]))
            .collect();
        let result = ransac_fit_plane(&pts, 30, 0.01).expect("should find plane");
        let len = length(result.normal);
        assert!((len - 1.0).abs() < 1e-6, "plane normal not unit: len={len}");
    }
    #[test]
    fn test_ransac_plane_few_points_returns_none() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let result = ransac_fit_plane(&pts, 10, 0.1);
        assert!(result.is_none(), "fewer than 3 points should return None");
    }
    #[test]
    fn test_ransac_plane_inlier_distance() {
        let pts: Vec<[f64; 3]> = (0..6)
            .flat_map(|i| (0..6).map(move |j| [i as f64, j as f64, 0.0]))
            .collect();
        let threshold = 0.05;
        let result = ransac_fit_plane(&pts, 20, threshold).expect("should find plane");
        let n = result.normal;
        let d = dot(n, result.point_on_plane);
        for &i in &result.inliers {
            let dist = (dot(n, pts[i]) - d).abs();
            assert!(
                dist <= threshold + 1e-9,
                "inlier {i} at distance {dist} > threshold {threshold}"
            );
        }
    }
    #[test]
    fn test_aabb_extent_unit_cube() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let (dims, center) = aabb_extent(&pts);
        for i in 0..3 {
            assert!((dims[i] - 1.0).abs() < 1e-9, "dims[{i}]={}", dims[i]);
            assert!((center[i] - 0.5).abs() < 1e-9, "center[{i}]={}", center[i]);
        }
    }
    #[test]
    fn test_aabb_extent_empty() {
        let (dims, center) = aabb_extent(&[]);
        assert_eq!(dims, [0.0; 3]);
        assert_eq!(center, [0.0; 3]);
    }
    #[test]
    fn test_pca_obb_axes_orthonormal() {
        let cloud = flat_grid(5, 5, 1.0);
        let (axes, _half_extents, _center) = pca_obb(&cloud.points);
        for (k, ax) in axes.iter().enumerate() {
            let len = length(*ax);
            assert!((len - 1.0).abs() < 0.1, "axis {k} not unit: len={len}");
        }
    }
    #[test]
    fn test_pca_obb_contains_all_points() {
        let pts: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64, j as f64, 0.0]))
            .collect();
        let (axes, half_extents, center) = pca_obb(&pts);
        for &p in &pts {
            let d = sub(p, center);
            for k in 0..3 {
                let proj = dot(d, axes[k]).abs();
                assert!(
                    proj <= half_extents[k] + 1e-6,
                    "point {:?} projection {proj} > half_extent {} on axis {k}",
                    p,
                    half_extents[k]
                );
            }
        }
    }
    #[test]
    fn test_pca_obb_point_cloud_method() {
        let cloud = flat_grid(4, 4, 1.0);
        let (axes, he, _center) = cloud.pca_obb();
        assert_eq!(axes.len(), 3);
        assert_eq!(he.len(), 3);
        for &h in &he {
            assert!(
                h.is_finite() && h >= 0.0,
                "half_extent must be finite and non-negative: {h}"
            );
        }
    }
    #[test]
    fn test_voxel_downsample_single_voxel() {
        let pts = vec![[0.1, 0.1, 0.0], [0.2, 0.1, 0.0], [0.1, 0.2, 0.0]];
        let cloud = PointCloud::from_points(pts.clone());
        let down = PointCloudFilter::voxel_downsample(&cloud, 1.0);
        assert_eq!(down.len(), 1, "all points in one voxel → 1 output point");
        let avg_x = pts.iter().map(|p| p[0]).sum::<f64>() / 3.0;
        let avg_y = pts.iter().map(|p| p[1]).sum::<f64>() / 3.0;
        assert!((down.points[0][0] - avg_x).abs() < 1e-9);
        assert!((down.points[0][1] - avg_y).abs() < 1e-9);
    }
    #[test]
    fn test_statistical_outlier_removal_preserves_dense_cluster() {
        let mut pts: Vec<[f64; 3]> = (0..5)
            .flat_map(|i| (0..5).map(move |j| [i as f64 * 0.1, j as f64 * 0.1, 0.0]))
            .collect();
        pts.push([1000.0, 1000.0, 1000.0]);
        let cloud = PointCloud::from_points(pts);
        let filtered = PointCloudFilter::statistical_outlier_removal(&cloud, 5, 1.0);
        for p in &filtered.points {
            assert!(p[0] < 100.0, "outlier survived: {:?}", p);
        }
    }
    #[test]
    fn test_knn_sorted_ascending() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let tree = KdTree3D::build(&pts);
        let result = tree.k_nearest([0.0, 0.0, 0.0], 5);
        for w in result.windows(2) {
            assert!(w[0].1 <= w[1].1 + 1e-12, "k-nearest not sorted: {:?}", w);
        }
    }
    #[test]
    fn test_icp_result_apply_to_cloud() {
        let src = flat_grid(3, 3, 1.0);
        let tgt = flat_grid(3, 3, 1.0);
        let result = IcpRegistration::align(&src, &tgt);
        let transformed = result.apply_to(&src);
        assert_eq!(
            transformed.len(),
            src.len(),
            "apply_to should preserve point count"
        );
    }
    #[test]
    fn test_point_cloud_from_points() {
        let pts = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let cloud = PointCloud::from_points(pts.clone());
        assert_eq!(cloud.len(), 2);
        for (i, &p) in pts.iter().enumerate() {
            for k in 0..3 {
                assert!((cloud.points[i][k] - p[k]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_estimate_normals_pca_populated() {
        let mut cloud = flat_grid(5, 5, 1.0);
        cloud.estimate_normals_pca(6);
        assert_eq!(
            cloud.normals.len(),
            cloud.points.len(),
            "normals should match point count"
        );
    }
    #[test]
    fn test_point_cloud_is_empty() {
        let cloud = PointCloud::new();
        assert!(cloud.is_empty());
        let mut cloud2 = PointCloud::new();
        cloud2.add_point([1.0, 0.0, 0.0]);
        assert!(!cloud2.is_empty());
    }
    #[test]
    fn test_principal_curvatures_count_matches_points() {
        let cloud = flat_grid(5, 5, 1.0);
        let curvs = cloud.compute_principal_curvatures(6);
        assert_eq!(
            curvs.len(),
            cloud.len(),
            "curvature output length should match point count"
        );
    }
    #[test]
    fn test_principal_curvatures_finite_values() {
        let cloud = flat_grid(5, 5, 1.0);
        let curvs = cloud.compute_principal_curvatures(5);
        for (i, (k1, k2)) in curvs.iter().enumerate() {
            assert!(k1.is_finite(), "kappa1 at point {i} is not finite: {k1}");
            assert!(k2.is_finite(), "kappa2 at point {i} is not finite: {k2}");
        }
    }
    #[test]
    fn test_principal_curvatures_flat_grid_small() {
        let cloud = flat_grid(6, 6, 1.0);
        let curvs = cloud.compute_principal_curvatures(8);
        for (k1, _k2) in &curvs {
            assert!(
                *k1 >= -1e-6,
                "smallest eigenvalue of covariance should be non-negative: {k1}"
            );
        }
    }
    #[test]
    fn test_principal_curvatures_empty_cloud() {
        let cloud = PointCloud::new();
        let curvs = cloud.compute_principal_curvatures(5);
        assert!(
            curvs.is_empty(),
            "empty cloud should yield empty curvatures"
        );
    }
    #[test]
    fn test_voxel_downsample_method_reduces_count() {
        let cloud = flat_grid(10, 10, 0.1);
        let down = cloud.voxel_downsample(1.0);
        assert!(
            down.len() < cloud.len(),
            "voxel_downsample should reduce point count: {} < {}",
            down.len(),
            cloud.len()
        );
    }
    #[test]
    fn test_voxel_downsample_method_nonempty_result() {
        let cloud = flat_grid(4, 4, 0.5);
        let down = cloud.voxel_downsample(0.1);
        assert!(!down.is_empty(), "downsampled cloud should not be empty");
    }
    #[test]
    fn test_icp_register_same_cloud_low_error() {
        let cloud = flat_grid(4, 4, 1.0);
        let result = cloud.icp_register(&cloud);
        assert!(
            result.final_error < 1e-3,
            "ICP same-cloud error should be near 0, got {}",
            result.final_error
        );
    }
    #[test]
    fn test_icp_register_returns_valid_transform() {
        let cloud = flat_grid(4, 4, 1.0);
        let result = cloud.icp_register(&cloud);
        for k in 0..3 {
            assert!(
                (result.rotation[k][k] - 1.0).abs() < 0.1,
                "rotation[{k}][{k}] should be close to 1 for identity, got {}",
                result.rotation[k][k]
            );
        }
    }
    #[test]
    fn test_icp_register_apply_to_preserves_count() {
        let src = flat_grid(3, 3, 1.0);
        let tgt = flat_grid(3, 3, 1.0);
        let result = src.icp_register(&tgt);
        let transformed = result.apply_to(&src);
        assert_eq!(
            transformed.len(),
            src.len(),
            "apply_to should preserve point count"
        );
    }
}
