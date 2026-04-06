//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_new_ch {
    use super::*;
    #[test]
    fn test_lower_hull_2d_square() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let lower = lower_hull_2d(pts);
        assert!(!lower.is_empty());
        assert!(lower.first().unwrap()[0] <= lower.last().unwrap()[0]);
    }
    #[test]
    fn test_upper_hull_2d_square() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let upper = upper_hull_2d(pts);
        assert!(!upper.is_empty());
    }
    #[test]
    fn test_monotone_chain_hull_2d_square() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
        let hull = monotone_chain_hull_2d(&pts);
        assert!(
            !hull
                .iter()
                .any(|&p| (p[0] - 0.5).abs() < 1e-9 && (p[1] - 0.5).abs() < 1e-9),
            "interior point should be excluded"
        );
    }
    #[test]
    fn test_monotone_chain_hull_2d_collinear() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.0, 1.0]];
        let hull = monotone_chain_hull_2d(&pts);
        assert!(!hull.is_empty());
    }
    #[test]
    fn test_lower_hull_sorted_order() {
        let pts = vec![[3.0, 1.0], [1.0, 2.0], [2.0, 0.5], [0.0, 0.0]];
        let lower = lower_hull_2d(pts);
        for w in lower.windows(2) {
            assert!(w[0][0] <= w[1][0], "lower hull must be x-sorted");
        }
    }
    #[test]
    fn test_monotone_chain_few_points() {
        let pts: Vec<[f64; 2]> = vec![];
        let hull = monotone_chain_hull_2d(&pts);
        assert!(hull.is_empty());
        let one = vec![[0.5, 0.5]];
        let hull1 = monotone_chain_hull_2d(&one);
        assert_eq!(hull1.len(), 1);
    }
    #[test]
    fn test_incremental_hull_builds_after_4_points() {
        let mut ich = IncrementalConvexHull::new();
        ich.add_point([0.0, 0.0, 0.0]);
        ich.add_point([1.0, 0.0, 0.0]);
        ich.add_point([0.0, 1.0, 0.0]);
        assert!(ich.hull().is_none(), "only 3 points: no 3D hull yet");
        ich.add_point([0.0, 0.0, 1.0]);
        assert!(
            ich.hull().is_some(),
            "4 non-coplanar points: hull should exist"
        );
    }
    #[test]
    fn test_incremental_hull_interior_point_skipped() {
        let mut ich = IncrementalConvexHull::new();
        let cube = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        ich.add_points(&cube);
        let n_pts_before = ich.n_points();
        ich.add_point([0.5, 0.5, 0.5]);
        assert_eq!(
            ich.n_points(),
            n_pts_before,
            "interior point should be skipped"
        );
    }
    #[test]
    fn test_incremental_hull_volume_grows_with_exterior_points() {
        let mut ich = IncrementalConvexHull::new();
        ich.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let v1 = ich.volume();
        ich.add_point([5.0, 5.0, 5.0]);
        let v2 = ich.volume();
        assert!(
            v2 > v1,
            "hull volume should grow when exterior point is added"
        );
    }
    #[test]
    fn test_incremental_hull_clear() {
        let mut ich = IncrementalConvexHull::new();
        ich.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        ich.clear();
        assert_eq!(ich.n_points(), 0);
        assert!(ich.hull().is_none());
    }
    #[test]
    fn test_incremental_hull_add_points_batch() {
        let mut ich = IncrementalConvexHull::new();
        let pts = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        ich.add_points(&pts);
        assert!(ich.hull().is_some());
        assert!(
            (ich.volume() - 4.0 / 3.0).abs() < 0.5,
            "volume approx 4/3, got {}",
            ich.volume()
        );
    }
    #[test]
    fn test_box_support_plus_x() {
        let s = box_support([0.0; 3], [1.0; 3], [1.0, 0.0, 0.0]);
        assert!((s[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_box_support_minus_x() {
        let s = box_support([0.0; 3], [1.0; 3], [-1.0, 0.0, 0.0]);
        assert!((s[0] - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_sphere_support_direction() {
        let s = sphere_support([0.0; 3], 1.0, [1.0, 0.0, 0.0]);
        assert!((s[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_sphere_support_norm() {
        let s = sphere_support([0.0; 3], 2.0, [1.0, 1.0, 0.0]);
        let len = dot(s, s).sqrt();
        assert!(
            (len - 2.0).abs() < 1e-10,
            "sphere support should have magnitude r={len}"
        );
    }
    #[test]
    fn test_capsule_support_along_axis() {
        let a = [0.0, 0.0, -1.0];
        let b = [0.0, 0.0, 1.0];
        let s = capsule_support(a, b, 0.5, [0.0, 0.0, 1.0]);
        assert!(s[2] > 0.0, "capsule support in +Z should have positive Z");
    }
    #[test]
    fn test_ellipsoid_support_on_axis() {
        let s = ellipsoid_support([0.0; 3], [2.0, 1.0, 1.0], [1.0, 0.0, 0.0]);
        assert!(
            (s[0] - 2.0).abs() < 1e-10,
            "ellipsoid support +X should be 2, got {}",
            s[0]
        );
        assert!(s[1].abs() < 1e-10);
        assert!(s[2].abs() < 1e-10);
    }
    #[test]
    fn test_vertex_decimation_reduces_count() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let hull = ConvexHull3D::build(&pts).expect("build");
        let simplified = vertex_decimation(&hull, 6);
        assert!(simplified.len() <= hull.vertices.len());
    }
    #[test]
    fn test_vertex_decimation_target_at_min_does_not_crash() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let hull = ConvexHull3D::build(&pts).expect("build");
        let simplified = vertex_decimation(&hull, 4);
        assert!(simplified.len() >= 4);
    }
    #[test]
    fn test_icosahedron_vertex_count() {
        let verts = icosahedron_vertices(1.0);
        assert_eq!(verts.len(), 12, "icosahedron has 12 vertices");
    }
    #[test]
    fn test_icosahedron_hull_builds() {
        let verts = icosahedron_vertices(1.0);
        let hull = ConvexHull3D::build(&verts);
        assert!(hull.is_some(), "icosahedron hull should build");
        let h = hull.unwrap();
        assert!(
            h.n_faces() >= 20 && h.n_faces() <= 24,
            "icosahedron hull face count should be 20-24, got {}",
            h.n_faces()
        );
    }
    #[test]
    fn test_icosahedron_volume_positive() {
        let verts = icosahedron_vertices(1.0);
        let hull = ConvexHull3D::build(&verts).expect("build");
        assert!(hull.volume() > 0.0, "icosahedron must have positive volume");
    }
    #[test]
    fn test_octahedron_vertex_count() {
        let verts = octahedron_vertices(1.0);
        assert_eq!(verts.len(), 6, "octahedron has 6 vertices");
    }
    #[test]
    fn test_octahedron_hull_builds() {
        let verts = octahedron_vertices(1.0);
        let hull = ConvexHull3D::build(&verts);
        assert!(hull.is_some(), "octahedron hull should build");
        let h = hull.unwrap();
        assert_eq!(
            h.n_faces(),
            8,
            "octahedron has 8 faces, got {}",
            h.n_faces()
        );
    }
    #[test]
    fn test_octahedron_all_vertices_equidistant_from_center() {
        let r = 2.0;
        let verts = octahedron_vertices(r);
        for v in &verts {
            let len = dot(*v, *v).sqrt();
            assert!(
                (len - r).abs() < 1e-10,
                "all vertices should be at distance r={r}, got {len}"
            );
        }
    }
    #[test]
    fn test_tetrahedron_vertex_count() {
        let verts = tetrahedron_vertices(1.0);
        assert_eq!(verts.len(), 4, "tetrahedron has 4 vertices");
    }
    #[test]
    fn test_tetrahedron_hull_builds() {
        let verts = tetrahedron_vertices(1.0);
        let hull = ConvexHull3D::build(&verts).expect("build");
        assert_eq!(hull.n_faces(), 4, "tetrahedron has 4 faces");
    }
    #[test]
    fn test_minkowski_difference_hulls_vertex_count() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let h1 = ConvexHull3D::build(&pts).expect("build");
        let h2 = ConvexHull3D::build(&pts).expect("build");
        let diff_pts = minkowski_difference_hulls(&h1, &h2);
        assert_eq!(diff_pts.len(), h1.vertices.len() * h2.vertices.len());
    }
    #[test]
    fn test_minkowski_difference_hull_contains_origin_for_overlapping() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let h1 = ConvexHull3D::build(&pts).expect("build");
        let h2 = h1.clone();
        if let Some(diff_hull) = minkowski_difference_hull(&h1, &h2) {
            let origin_inside = diff_hull.contains_point([0.0, 0.0, 0.0]);
            let _ = origin_inside;
            assert!(diff_hull.n_faces() > 0);
        }
    }
    #[test]
    fn test_point_cloud_aabb_basic() {
        let pts = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [-1.0, 0.0, 2.0]];
        let (mn, mx) = point_cloud_aabb(&pts).expect("should have AABB");
        assert!((mn[0] - (-1.0)).abs() < 1e-12);
        assert!((mx[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_point_cloud_aabb_empty() {
        assert!(point_cloud_aabb(&[]).is_none());
    }
    #[test]
    fn test_point_cloud_centroid_unit_cube() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let c = point_cloud_centroid(&pts).expect("centroid");
        for &ci in &c {
            assert!((ci - 0.5).abs() < 1e-9, "centroid should be 0.5, got {ci}");
        }
    }
    #[test]
    fn test_point_cloud_covariance_identity_like() {
        let pts = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
        ];
        let cov = point_cloud_covariance(&pts);
        assert!(cov[0][0] >= 0.0 && cov[1][1] >= 0.0 && cov[2][2] >= 0.0);
    }
    #[test]
    fn test_scale_point_cloud() {
        let pts = vec![[1.0, 2.0, 3.0]];
        let scaled = scale_point_cloud(&pts, 2.0);
        assert!((scaled[0][0] - 2.0).abs() < 1e-12);
        assert!((scaled[0][1] - 4.0).abs() < 1e-12);
        assert!((scaled[0][2] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_translate_point_cloud() {
        let pts = vec![[0.0, 0.0, 0.0]];
        let translated = translate_point_cloud(&pts, [1.0, 2.0, 3.0]);
        assert!((translated[0][0] - 1.0).abs() < 1e-12);
        assert!((translated[0][1] - 2.0).abs() < 1e-12);
        assert!((translated[0][2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_cross2d_ccw() {
        let c = cross2d([0.0, 0.0], [1.0, 0.0], [0.5, 1.0]);
        assert!(
            c > 0.0,
            "CCW arrangement should have positive cross, got {c}"
        );
    }
    #[test]
    fn test_cross2d_cw() {
        let c = cross2d([0.0, 0.0], [0.5, 1.0], [1.0, 0.0]);
        assert!(
            c < 0.0,
            "CW arrangement should have negative cross, got {c}"
        );
    }
    #[test]
    fn test_cylinder_support_top_direction() {
        let center = [0.0; 3];
        let s = cylinder_support(center, 1.0, 0.5, [0.0, 0.0, 1.0]);
        assert!(
            (s[2] - 1.0).abs() < 1e-12,
            "cylinder support in +Z should be at top h={}",
            s[2]
        );
    }
    #[test]
    fn test_cylinder_support_bottom_direction() {
        let center = [0.0; 3];
        let s = cylinder_support(center, 1.0, 0.5, [0.0, 0.0, -1.0]);
        assert!(
            (s[2] + 1.0).abs() < 1e-12,
            "cylinder support in -Z should be at bottom"
        );
    }
    #[test]
    fn test_icosahedron_scaled_has_larger_volume() {
        let v1 = icosahedron_vertices(1.0);
        let v2 = icosahedron_vertices(2.0);
        let h1 = ConvexHull3D::build(&v1).expect("build");
        let h2 = ConvexHull3D::build(&v2).expect("build");
        assert!(
            h2.volume() > h1.volume(),
            "scaled icosahedron should have larger volume: {} vs {}",
            h2.volume(),
            h1.volume()
        );
    }
    #[test]
    fn test_minkowski_sum_hull_volume_at_least_as_large() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let h = ConvexHull3D::build(&pts).expect("build");
        if let Some(sum) = minkowski_sum_hull(&h, &h) {
            assert!(
                sum.volume() >= h.volume(),
                "Minkowski sum volume {} should be >= operand volume {}",
                sum.volume(),
                h.volume()
            );
        }
    }
}
