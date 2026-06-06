//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::point_cloud_io::functions::{cross3, dist, normalize3};
    use crate::point_cloud_io::*;
    fn make_cube_cloud(n_per_face: usize) -> PointCloud {
        let mut cloud = PointCloud::new();
        let step = 1.0 / n_per_face as f64;
        for i in 0..n_per_face {
            for j in 0..n_per_face {
                let u = i as f64 * step;
                let v = j as f64 * step;
                cloud.positions.push([u, v, 0.0]);
                cloud.positions.push([u, v, 1.0]);
                cloud.positions.push([0.0, u, v]);
                cloud.positions.push([1.0, u, v]);
                cloud.positions.push([u, 0.0, v]);
                cloud.positions.push([u, 1.0, v]);
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
        assert_eq!(idx, 0);
        assert!(d_sq < 0.1);
    }
    #[test]
    fn test_kdtree_k_nearest() {
        let cloud = make_simple_cloud();
        let tree = KdTree::build(&cloud);
        let results = tree.k_nearest([0.0, 0.0, 0.0], 3);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].0, 0);
        assert!(results[0].1 < 1e-10);
    }
    #[test]
    fn test_kdtree_radius_search() {
        let cloud = make_simple_cloud();
        let tree = KdTree::build(&cloud);
        let results = tree.radius_search([0.0, 0.0, 0.0], 1.01);
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
        let cloud = make_cube_cloud(10);
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
        assert!(bb_down.min[0] >= bb_orig.min[0] - 0.5);
        assert!(bb_down.max[0] <= bb_orig.max[0] + 0.5);
    }
    #[test]
    fn test_statistical_outlier_removal() {
        let mut cloud = make_simple_cloud();
        cloud.positions.push([100.0, 100.0, 100.0]);
        let filtered = statistical_outlier_removal(&cloud, 3, 1.0);
        assert!(filtered.len() <= cloud.len());
        for p in &filtered.positions {
            assert!(p[0] < 50.0, "Outlier not removed: {p:?}");
        }
    }
    #[test]
    fn test_estimate_normals() {
        let mut cloud = PointCloud::new();
        for i in 0..10 {
            for j in 0..10 {
                cloud.positions.push([i as f64, j as f64, 0.0]);
            }
        }
        let normals = estimate_normals(&cloud, 6);
        assert_eq!(normals.len(), 100);
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
        assert!(buf.starts_with(b"ply"));
        assert!(buf.len() > 12);
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
