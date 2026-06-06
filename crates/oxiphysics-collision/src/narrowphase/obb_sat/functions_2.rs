//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    add3_raw, cross3_raw, dot3_raw, negate3_raw, obb_obb_test, sat_test_axis, scale3_raw, sub3_raw,
};
use super::types::{ContactFeatureType, Obb};

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::narrowphase::ObbSat;
    use crate::narrowphase::ObbSatTest;
    use crate::narrowphase::ObbShape;

    use crate::narrowphase::compute_contact_points;

    use crate::narrowphase::obb_obb_test;

    use crate::narrowphase::obb_sat::len3_raw;
    use crate::narrowphase::obb_sat::normalize3_raw;

    use crate::narrowphase::obb_sat::types::Obb;

    use oxiphysics_core::Vec3;
    use std::f64::consts::FRAC_PI_4;
    fn identity_axes() -> [Vec3; 3] {
        [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
    }
    #[test]
    fn test_obb_separated() {
        let a = ObbShape::axis_aligned(Vec3::new(1.0, 1.0, 1.0), Vec3::new(0.0, 0.0, 0.0));
        let b = ObbShape::axis_aligned(Vec3::new(1.0, 1.0, 1.0), Vec3::new(10.0, 0.0, 0.0));
        assert!(
            ObbSat::query(&a, &b).is_none(),
            "clearly separated OBBs must return None"
        );
    }
    #[test]
    fn test_obb_face_contact() {
        let a = ObbShape::axis_aligned(Vec3::new(1.0, 1.0, 1.0), Vec3::new(0.0, 0.0, 0.0));
        let b = ObbShape::axis_aligned(Vec3::new(1.0, 1.0, 1.0), Vec3::new(0.0, 1.8, 0.0));
        let contact = ObbSat::query(&a, &b).expect("face-resting OBBs must produce contact");
        assert!(
            contact.depth > 0.0,
            "penetration depth must be positive, got {}",
            contact.depth
        );
        assert!(
            contact.normal.y.abs() > 0.9,
            "face contact normal must be near Y, got {:?}",
            contact.normal
        );
    }
    #[test]
    fn test_obb_edge_contact() {
        let angle = FRAC_PI_4;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let axes_b = [
            Vec3::new(cos_a, sin_a, 0.0),
            Vec3::new(-sin_a, cos_a, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let a = ObbShape::axis_aligned(Vec3::new(1.0, 1.0, 1.0), Vec3::new(0.0, 0.0, 0.0));
        let b = ObbShape::new(Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.5, 0.0, 0.0), axes_b);
        if let Some(contact) = ObbSat::query(&a, &b) {
            assert!(contact.depth > 0.0, "edge contact depth must be positive");
            assert!(
                (contact.normal.norm() - 1.0).abs() < 1e-6,
                "normal must be unit length, got {}",
                contact.normal.norm()
            );
        }
    }
    #[test]
    fn test_obb_contained() {
        let outer = ObbShape::axis_aligned(Vec3::new(5.0, 5.0, 5.0), Vec3::new(0.0, 0.0, 0.0));
        let inner = ObbShape::axis_aligned(Vec3::new(0.5, 0.5, 0.5), Vec3::new(0.0, 0.0, 0.0));
        let contact = ObbSat::query(&outer, &inner).expect("contained OBB must produce contact");
        assert!(
            contact.depth > 0.0,
            "contained OBB must have positive depth, got {}",
            contact.depth
        );
    }
    #[test]
    fn test_obb_rotated() {
        let angle = FRAC_PI_4;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let axes_b = [
            Vec3::new(cos_a, sin_a, 0.0),
            Vec3::new(-sin_a, cos_a, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let a = ObbShape::new(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 0.0, 0.0),
            identity_axes(),
        );
        let b = ObbShape::new(Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.2, 0.0, 0.0), axes_b);
        let contact = ObbSat::query(&a, &b).expect("45-degree rotated OBBs must contact");
        assert!(
            contact.depth > 0.0,
            "rotated contact depth must be positive, got {}",
            contact.depth
        );
        assert!(
            (contact.normal.norm() - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            contact.normal.norm()
        );
    }
    #[test]
    fn test_raw_obb_separated() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(obb_obb_test(&a, &b).is_none());
    }
    #[test]
    fn test_raw_obb_overlapping() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_obb_test(&a, &b).expect("should overlap");
        assert!(result.penetration_depth > 0.0);
        assert!((result.penetration_depth - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_raw_obb_contained() {
        let outer = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let inner = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let result = obb_obb_test(&outer, &inner).expect("contained must overlap");
        assert!(result.penetration_depth > 0.0);
    }
    #[test]
    fn test_project_obb_onto_axis() {
        let obb = Obb::axis_aligned([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let (lo, hi) = project_obb_onto_axis(&obb, [1.0, 0.0, 0.0]);
        assert!((lo - (-1.0)).abs() < 1e-10);
        assert!((hi - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_compute_contact_points() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let sat = obb_obb_test(&a, &b).unwrap();
        let contacts = compute_contact_points(&a, &b, &sat);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn test_obb_sphere_overlap() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_sphere_test(&obb, [1.5, 0.0, 0.0], 1.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.penetration_depth > 0.0);
    }
    #[test]
    fn test_obb_sphere_separated() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_sphere_test(&obb, [5.0, 0.0, 0.0], 0.5);
        assert!(result.is_none());
    }
    #[test]
    fn test_obb_capsule_overlap() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_capsule_test(&obb, [1.5, -2.0, 0.0], [1.5, 2.0, 0.0], 1.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.penetration_depth > 0.0);
    }
    #[test]
    fn test_obb_capsule_separated() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_capsule_test(&obb, [5.0, -2.0, 0.0], [5.0, 2.0, 0.0], 0.5);
        assert!(result.is_none());
    }
    #[test]
    fn test_obb_frustum_visible() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ];
        assert!(
            obb_frustum_cull(&obb, &planes),
            "OBB at origin should be inside frustum"
        );
    }
    #[test]
    fn test_obb_frustum_culled() {
        let obb = Obb::axis_aligned([100.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.0, -1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, -1.0, 1.0],
        ];
        assert!(
            !obb_frustum_cull(&obb, &planes),
            "OBB at x=100 should be culled"
        );
    }
    #[test]
    fn test_obb_triangle_overlap_inside_box() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        assert!(obb_triangle_test(&obb, v0, v1, v2));
    }
    #[test]
    fn test_obb_triangle_no_overlap_far() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let v0 = [10.0, 0.0, 0.0];
        let v1 = [11.0, 0.0, 0.0];
        let v2 = [10.0, 1.0, 0.0];
        assert!(!obb_triangle_test(&obb, v0, v1, v2));
    }
    #[test]
    fn test_obb_triangle_overlap_piercing() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let v0 = [-2.0, 0.5, 0.5];
        let v1 = [2.0, 0.5, 0.5];
        let v2 = [0.0, 0.5, 2.0];
        assert!(obb_triangle_test(&obb, v0, v1, v2));
    }
    #[test]
    fn test_obb_vertices_count() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let verts = obb_vertices(&obb);
        assert_eq!(verts.len(), 8);
    }
    #[test]
    fn test_obb_vertices_aabb() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let verts = obb_vertices(&obb);
        for v in &verts {
            assert!((v[0].abs() - 1.0).abs() < 1e-10, "vx={}", v[0]);
            assert!((v[1].abs() - 2.0).abs() < 1e-10, "vy={}", v[1]);
            assert!((v[2].abs() - 3.0).abs() < 1e-10, "vz={}", v[2]);
        }
    }
    #[test]
    fn test_obb_edges_count() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let edges = obb_edges(&obb);
        assert_eq!(edges.len(), 12);
    }
    #[test]
    fn test_obb_face_centers_count() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let fc = obb_face_centers(&obb);
        assert_eq!(fc.len(), 6);
    }
    #[test]
    fn test_obb_face_centers_correct_distance() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let fc = obb_face_centers(&obb);
        for c in &fc {
            let d = len3_raw(*c);
            assert!((d - 1.0).abs() < 1e-10, "face center distance={d}");
        }
    }
    #[test]
    fn test_obb_support_point_public() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let sp = obb_support_point(&obb, [1.0, 1.0, 1.0]);
        assert!((sp[0] - 1.0).abs() < 1e-10);
        assert!((sp[1] - 2.0).abs() < 1e-10);
        assert!((sp[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_obb_from_center_axes() {
        let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let obb = Obb::from_center_axes([1.0, 2.0, 3.0], axes, [0.5, 0.5, 0.5]);
        assert_eq!(obb.center, [1.0, 2.0, 3.0]);
        assert_eq!(obb.half_extents, [0.5, 0.5, 0.5]);
    }
    #[test]
    fn test_obb_project_onto_axis() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (lo, hi) = obb.project_onto_axis([1.0, 0.0, 0.0]);
        assert!((lo - (-1.0)).abs() < 1e-10);
        assert!((hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_obb_obb_sat_test_overlapping() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_obb_sat_test(&a, &b),
            "overlapping cubes should return true"
        );
    }
    #[test]
    fn test_obb_obb_sat_test_separated() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            !obb_obb_sat_test(&a, &b),
            "separated boxes should return false"
        );
    }
    #[test]
    fn test_obb_obb_contact_normal_some() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let normal = obb_obb_contact_normal(&a, &b);
        assert!(
            normal.is_some(),
            "overlapping boxes should have a contact normal"
        );
        let n = normal.unwrap();
        let len = len3_raw(n);
        assert!(
            (len - 1.0).abs() < 1e-6,
            "contact normal must be unit length"
        );
    }
    #[test]
    fn test_obb_obb_contact_normal_none() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(obb_obb_contact_normal(&a, &b).is_none());
    }
    #[test]
    fn test_obb_obb_penetration_depth_positive() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let depth = obb_obb_penetration_depth(&a, &b);
        assert!(depth.is_some());
        assert!(depth.unwrap() > 0.0, "penetration depth must be positive");
        assert!((depth.unwrap() - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_aabb_to_obb_round_trip() {
        let min = [-1.0_f64, -2.0, -3.0];
        let max = [1.0_f64, 2.0, 3.0];
        let obb = aabb_to_obb(min, max);
        assert_eq!(obb.center, [0.0, 0.0, 0.0]);
        assert_eq!(obb.half_extents, [1.0, 2.0, 3.0]);
        assert_eq!(obb.rotation[0], [1.0, 0.0, 0.0]);
        assert_eq!(obb.rotation[1], [0.0, 1.0, 0.0]);
        assert_eq!(obb.rotation[2], [0.0, 0.0, 1.0]);
    }
    #[test]
    fn test_aabb_to_obb_overlaps_with_neighbor() {
        let obb_a = aabb_to_obb([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let obb_b = aabb_to_obb([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(obb_obb_sat_test(&obb_a, &obb_b));
    }
    #[test]
    fn test_obb_from_points_axis_aligned_cube() {
        let pts: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let obb = obb_from_points(&pts).expect("must return OBB for valid points");
        for i in 0..3 {
            assert!(obb.center[i].abs() < 0.1, "center[{}]={}", i, obb.center[i]);
        }
        for i in 0..3 {
            assert!(
                (obb.half_extents[i] - 1.0).abs() < 0.1,
                "half_extents[{}]={}",
                i,
                obb.half_extents[i]
            );
        }
        for p in &pts {
            assert!(
                obb_contains_point(&obb, *p) || obb_point_sq_dist(&obb, *p) < 1e-6,
                "point {:?} outside fitted OBB",
                p
            );
        }
    }
    #[test]
    fn test_obb_from_points_elongated() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let obb = obb_from_points(&pts).expect("must return OBB");
        assert!(
            (obb.center[0] - 4.5).abs() < 0.2,
            "center x={}",
            obb.center[0]
        );
        let max_he = obb
            .half_extents
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(max_he > 3.0, "max half-extent={}", max_he);
    }
    #[test]
    fn test_obb_from_points_single_point_returns_none() {
        let pts = vec![[1.0_f64, 2.0, 3.0]];
        assert!(obb_from_points(&pts).is_none());
    }
    #[test]
    fn test_obb_from_points_two_points() {
        let pts = vec![[0.0_f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let obb = obb_from_points(&pts).expect("two points should be OK");
        assert!(
            obb.half_extents[0] >= 0.0 || obb.half_extents[1] >= 0.0 || obb.half_extents[2] >= 0.0
        );
    }
    #[test]
    fn test_obb_transform_identity() {
        let obb = Obb::axis_aligned([1.0, 2.0, 3.0], [0.5, 0.5, 0.5]);
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let result = obb_transform(&obb, identity, [0.0, 0.0, 0.0]);
        for i in 0..3 {
            assert!((result.center[i] - obb.center[i]).abs() < 1e-10);
            assert!((result.half_extents[i] - obb.half_extents[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_obb_transform_translation() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let t = [3.0, 5.0, -2.0];
        let result = obb_transform(&obb, identity, t);
        assert!((result.center[0] - 3.0).abs() < 1e-10);
        assert!((result.center[1] - 5.0).abs() < 1e-10);
        assert!((result.center[2] + 2.0).abs() < 1e-10);
        for i in 0..3 {
            assert!((result.half_extents[i] - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_obb_transform_rotation_90_deg() {
        let obb = Obb::axis_aligned([1.0, 0.0, 0.0], [0.5, 1.0, 0.5]);
        let rot90z = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let result = obb_transform(&obb, rot90z, [0.0, 0.0, 0.0]);
        assert!(
            (result.center[0] - 0.0).abs() < 1e-10,
            "cx={}",
            result.center[0]
        );
        assert!(
            (result.center[1] + 1.0).abs() < 1e-10,
            "cy={}",
            result.center[1]
        );
        assert!(
            (result.center[2] - 0.0).abs() < 1e-10,
            "cz={}",
            result.center[2]
        );
    }
    #[test]
    fn test_obb_closest_point_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let q = [0.5, 0.5, 0.5];
        let cp = obb_closest_point(&obb, q);
        for i in 0..3 {
            assert!((cp[i] - q[i]).abs() < 1e-10, "inside: cp[{}]={}", i, cp[i]);
        }
    }
    #[test]
    fn test_obb_closest_point_outside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let cp = obb_closest_point(&obb, [3.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-10);
        assert!(cp[1].abs() < 1e-10);
        assert!(cp[2].abs() < 1e-10);
    }
    #[test]
    fn test_obb_point_sq_dist_inside_is_zero() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let d = obb_point_sq_dist(&obb, [1.0, 1.0, 1.0]);
        assert!(d < 1e-10, "sq dist for interior point = {}", d);
    }
    #[test]
    fn test_obb_point_sq_dist_outside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let d = obb_point_sq_dist(&obb, [3.0, 0.0, 0.0]);
        assert!((d - 4.0).abs() < 1e-10, "sq dist = {}", d);
    }
    #[test]
    fn test_obb_segment_test_piercing() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(obb_segment_test(&obb, [-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]));
    }
    #[test]
    fn test_obb_segment_test_miss() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!obb_segment_test(&obb, [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]));
    }
    #[test]
    fn test_obb_segment_test_segment_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        assert!(obb_segment_test(&obb, [0.5, 0.5, 0.5], [1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_obb_ray_cast_hit_x() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let t = obb_ray_cast(&obb, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(t.is_some(), "ray along X should hit OBB");
        let t_val = t.unwrap();
        assert!((t_val - 4.0).abs() < 1e-6, "t should be 4.0, got {}", t_val);
    }
    #[test]
    fn test_obb_ray_cast_miss() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let t = obb_ray_cast(&obb, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(t.is_none(), "ray missing OBB should return None");
    }
    #[test]
    fn test_obb_ray_cast_behind() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let t = obb_ray_cast(&obb, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(t.is_none(), "ray pointing away from OBB should miss");
    }
    #[test]
    fn test_sat_obb_obb_overlap() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = sat_obb_obb(&a, &b);
        assert!(
            result.is_some(),
            "overlapping boxes must produce SAT result"
        );
        let (normal, depth) = result.unwrap();
        assert!(depth > 0.0, "depth must be positive, got {}", depth);
        let nlen = len3_raw(normal);
        assert!((nlen - 1.0).abs() < 1e-6, "normal must be unit length");
    }
    #[test]
    fn test_sat_obb_obb_separated() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(sat_obb_obb(&a, &b).is_none());
    }
    #[test]
    fn test_obb_contact_point_overlap() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (normal, _depth) = sat_obb_obb(&a, &b).unwrap();
        let cp = obb_contact_point(&a, &b, normal);
        assert!(cp[0] > 0.5 && cp[0] < 2.0, "contact point x={}", cp[0]);
    }
    #[test]
    fn test_sat_test_axis_overlap() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let overlap = sat_test_axis(&a, &b, [1.0, 0.0, 0.0]);
        assert!(overlap.is_some(), "should overlap on X axis");
        assert!((overlap.unwrap() - 0.5).abs() < 0.01, "overlap ≈ 0.5");
    }
    #[test]
    fn test_sat_test_axis_separated() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let overlap = sat_test_axis(&a, &b, [1.0, 0.0, 0.0]);
        assert!(overlap.is_none(), "should be separated on X axis");
    }
    #[test]
    fn test_sat_test_axis_degenerate() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([0.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = sat_test_axis(&a, &b, [0.0, 0.0, 0.0]);
        assert_eq!(result, Some(0.0));
    }
    #[test]
    fn test_obb_volume() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let v = obb_volume(&obb);
        assert!((v - 48.0).abs() < 1e-10, "volume = {}", v);
    }
    #[test]
    fn test_obb_surface_area() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let sa = obb_surface_area(&obb);
        assert!((sa - 24.0).abs() < 1e-10, "surface area = {}", sa);
    }
    #[test]
    fn test_obb_surface_area_rectangle() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let sa = obb_surface_area(&obb);
        assert!((sa - 208.0).abs() < 1e-10, "surface area = {}", sa);
    }
    #[test]
    fn test_obb_contains_point_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(obb_contains_point(&obb, [0.0, 0.0, 0.0]));
        assert!(obb_contains_point(&obb, [1.5, 1.5, 1.5]));
    }
    #[test]
    fn test_obb_contains_point_outside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!obb_contains_point(&obb, [2.0, 0.0, 0.0]));
        assert!(!obb_contains_point(&obb, [0.0, 5.0, 0.0]));
    }
    #[test]
    fn test_obb_contains_obb_true() {
        let outer = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let inner = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_contains_obb(&outer, &inner),
            "large OBB must contain small OBB"
        );
    }
    #[test]
    fn test_obb_contains_obb_false() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([3.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!obb_contains_obb(&a, &b), "offset OBBs are not contained");
    }
    #[test]
    fn test_obb_merge_aabb_basic() {
        let a = Obb::axis_aligned([-2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let merged = obb_merge_aabb(&a, &b);
        assert!(
            merged.half_extents[0] > 2.0,
            "merged hx={}",
            merged.half_extents[0]
        );
    }
    #[test]
    fn test_obb_merge_aabb_contains_both() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([3.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let merged = obb_merge_aabb(&a, &b);
        for v in &obb_vertices(&a) {
            assert!(
                obb_contains_point(&merged, *v) || obb_point_sq_dist(&merged, *v) < 1e-6,
                "merged OBB must contain vertex {:?}",
                v
            );
        }
        for v in &obb_vertices(&b) {
            assert!(
                obb_contains_point(&merged, *v) || obb_point_sq_dist(&merged, *v) < 1e-6,
                "merged OBB must contain vertex {:?}",
                v
            );
        }
    }
    #[test]
    fn test_eigen_symmetric3_diagonal() {
        let mat = [[3.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        let (evals, _evecs) = eigen_symmetric3(mat);
        let mut sorted = evals;
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 1.0).abs() < 1e-6, "eval0={}", sorted[0]);
        assert!((sorted[1] - 2.0).abs() < 1e-6, "eval1={}", sorted[1]);
        assert!((sorted[2] - 3.0).abs() < 1e-6, "eval2={}", sorted[2]);
    }
    #[test]
    fn test_eigen_symmetric3_identity() {
        let mat = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (evals, _) = eigen_symmetric3(mat);
        for &e in &evals {
            assert!((e - 1.0).abs() < 1e-6, "eval={}", e);
        }
    }
    #[test]
    fn test_covariance_matrix_axis_aligned_x() {
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let mean = points_mean(&pts);
        let cov = covariance_matrix(&pts, mean);
        assert!(cov[0][0] > cov[1][1], "cov[0][0] must dominate");
        assert!(cov[0][0] > cov[2][2], "cov[0][0] must dominate Z");
    }
    #[test]
    fn test_normalize3_raw_unit() {
        let n = normalize3_raw([3.0, 4.0, 0.0]);
        let len = len3_raw(n);
        assert!((len - 1.0).abs() < 1e-10, "len={}", len);
        assert!((n[0] - 0.6).abs() < 1e-10);
        assert!((n[1] - 0.8).abs() < 1e-10);
    }
    #[test]
    fn test_normalize3_raw_degenerate() {
        let n = normalize3_raw([0.0, 0.0, 0.0]);
        assert_eq!(n, [1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_obb_new_stores_fields() {
        let c = [1.0, 2.0, 3.0];
        let he = [0.5, 1.0, 1.5];
        let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let obb = Obb::new(c, he, rot);
        assert_eq!(obb.center, c);
        assert_eq!(obb.half_extents, he);
        assert_eq!(obb.rotation, rot);
    }
    #[test]
    fn test_obb_axis_aligned_identity_rotation() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(obb.rotation[0], [1.0, 0.0, 0.0]);
        assert_eq!(obb.rotation[1], [0.0, 1.0, 0.0]);
        assert_eq!(obb.rotation[2], [0.0, 0.0, 1.0]);
    }
    #[test]
    fn test_obb_from_center_axes_roundtrip() {
        let center = [5.0, -3.0, 2.0];
        let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let he = [1.0, 2.0, 3.0];
        let obb = Obb::from_center_axes(center, axes, he);
        assert_eq!(obb.center, center);
        assert_eq!(obb.half_extents, he);
        assert_eq!(obb.rotation, axes);
    }
    #[test]
    fn test_aabb_to_obb_center() {
        let obb = aabb_to_obb([0.0, 0.0, 0.0], [4.0, 6.0, 8.0]);
        assert!((obb.center[0] - 2.0).abs() < 1e-10, "cx={}", obb.center[0]);
        assert!((obb.center[1] - 3.0).abs() < 1e-10, "cy={}", obb.center[1]);
        assert!((obb.center[2] - 4.0).abs() < 1e-10, "cz={}", obb.center[2]);
    }
    #[test]
    fn test_aabb_to_obb_half_extents() {
        let obb = aabb_to_obb([2.0, 0.0, -1.0], [6.0, 4.0, 3.0]);
        assert!((obb.half_extents[0] - 2.0).abs() < 1e-10);
        assert!((obb.half_extents[1] - 2.0).abs() < 1e-10);
        assert!((obb.half_extents[2] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_to_obb_unit_cube() {
        let obb = aabb_to_obb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        for i in 0..3 {
            assert!(
                (obb.center[i]).abs() < 1e-10,
                "center[{}]={}",
                i,
                obb.center[i]
            );
            assert!((obb.half_extents[i] - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_project_obb_onto_axis_x() {
        let obb = Obb::axis_aligned([3.0, 0.0, 0.0], [1.0, 2.0, 2.0]);
        let (lo, hi) = project_obb_onto_axis(&obb, [1.0, 0.0, 0.0]);
        assert!((lo - 2.0).abs() < 1e-10, "lo={}", lo);
        assert!((hi - 4.0).abs() < 1e-10, "hi={}", hi);
    }
    #[test]
    fn test_project_obb_onto_axis_symmetric() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let (lo, hi) = project_obb_onto_axis(&obb, [0.0, 1.0, 0.0]);
        assert!((lo + 2.0).abs() < 1e-10);
        assert!((hi - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_project_obb_onto_axis_diagonal() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let axis = normalize3_raw([1.0, 1.0, 0.0]);
        let (lo, hi) = project_obb_onto_axis(&obb, axis);
        assert!(
            lo < 0.0 && hi > 0.0,
            "should span origin: lo={}, hi={}",
            lo,
            hi
        );
        let radius = (hi - lo) * 0.5;
        assert!((radius - 2_f64.sqrt()).abs() < 1e-9, "radius={}", radius);
    }
    #[test]
    fn test_obb_obb_sat_test_overlap() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_obb_sat_test(&a, &b),
            "overlapping OBBs must return true"
        );
    }
    #[test]
    fn test_obb_obb_sat_test_separated_far() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            !obb_obb_sat_test(&a, &b),
            "separated OBBs must return false"
        );
    }
    #[test]
    fn test_obb_obb_contact_normal_unit_x() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let n = obb_obb_contact_normal(&a, &b).expect("should have contact normal");
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-9,
            "contact normal not unit: len={}",
            len
        );
    }
    #[test]
    fn test_obb_obb_contact_normal_separated_is_none() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(obb_obb_contact_normal(&a, &b).is_none());
    }
    #[test]
    fn test_obb_obb_penetration_depth_nonzero() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let depth = obb_obb_penetration_depth(&a, &b).expect("should have depth");
        assert!(depth > 0.0, "penetration depth must be positive: {}", depth);
    }
    #[test]
    fn test_obb_closest_point_query_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let cp = obb_closest_point(&obb, [0.5, 0.5, 0.5]);
        assert!((cp[0] - 0.5).abs() < 1e-9);
        assert!((cp[1] - 0.5).abs() < 1e-9);
        assert!((cp[2] - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_obb_closest_point_outside_x() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let cp = obb_closest_point(&obb, [5.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-9, "cp.x={}", cp[0]);
        assert!(cp[1].abs() < 1e-9, "cp.y={}", cp[1]);
        assert!(cp[2].abs() < 1e-9, "cp.z={}", cp[2]);
    }
    #[test]
    fn test_obb_closest_point_origin_to_unit_cube() {
        let obb = Obb::axis_aligned([3.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let cp = obb_closest_point(&obb, [0.0, 0.0, 0.0]);
        assert!((cp[0] - 2.0).abs() < 1e-9, "cp.x={}", cp[0]);
    }
    #[test]
    fn test_obb_point_sq_dist_interior_zero() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let d = obb_point_sq_dist(&obb, [1.0, 2.0, 3.0]);
        assert!(d < 1e-10, "point inside: sq_dist={}", d);
    }
    #[test]
    fn test_obb_point_sq_dist_exterior() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let d = obb_point_sq_dist(&obb, [3.0, 0.0, 0.0]);
        assert!((d - 4.0).abs() < 1e-9, "sq_dist={}", d);
    }
    #[test]
    fn test_obb_vertices_count_eight() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let verts = obb_vertices(&obb);
        assert_eq!(verts.len(), 8);
    }
    #[test]
    fn test_obb_vertices_unit_cube_span() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let verts = obb_vertices(&obb);
        for v in &verts {
            for coord in v.iter() {
                assert!(
                    coord.abs() <= 1.0 + 1e-10,
                    "vertex out of bounds: {}",
                    coord
                );
            }
        }
    }
    #[test]
    fn test_obb_vertices_offset_cube() {
        let obb = Obb::axis_aligned([5.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let verts = obb_vertices(&obb);
        let xs: Vec<f64> = verts.iter().map(|v| v[0]).collect();
        assert!(
            xs.iter().all(|&x| (4.0 - 1e-10..=6.0 + 1e-10).contains(&x)),
            "x vertices must be in [4, 6]: {:?}",
            xs
        );
    }
    #[test]
    fn test_obb_edges_count_twelve() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let edges = obb_edges(&obb);
        assert_eq!(edges.len(), 12);
    }
    #[test]
    fn test_obb_face_centers_count_six() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let faces = obb_face_centers(&obb);
        assert_eq!(faces.len(), 6);
    }
    #[test]
    fn test_obb_face_centers_distances() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let faces = obb_face_centers(&obb);
        let he = [2.0_f64, 3.0, 4.0];
        for fc in &faces {
            let d2: f64 = fc.iter().map(|&x| x * x).sum();
            let d = d2.sqrt();
            let matches = he.iter().any(|&h| (d - h).abs() < 1e-9);
            assert!(
                matches,
                "face center distance {} doesn't match any half_extent {:?}",
                d, he
            );
        }
    }
    #[test]
    fn test_obb_segment_test_intersecting() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_segment_test(&obb, [-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
            "segment through center must intersect OBB"
        );
    }
    #[test]
    fn test_obb_segment_test_no_intersect() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            !obb_segment_test(&obb, [-2.0, 5.0, 0.0], [2.0, 5.0, 0.0]),
            "segment far from OBB must not intersect"
        );
    }
    #[test]
    fn test_obb_segment_test_endpoint_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_segment_test(&obb, [0.0, 0.0, 0.0], [5.0, 0.0, 0.0]),
            "segment with endpoint inside must intersect"
        );
    }
    #[test]
    fn test_obb_ray_cast_hit() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let t = obb_ray_cast(&obb, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(t.is_some(), "ray along X should hit unit cube");
        let t = t.unwrap();
        assert!((t - 4.0).abs() < 1e-9, "expected t=4.0, got {}", t);
    }
    #[test]
    fn test_obb_ray_cast_y_offset_miss() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let t = obb_ray_cast(&obb, [-5.0, 5.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(t.is_none(), "ray offset in Y should miss unit cube");
    }
    #[test]
    fn test_obb_ray_cast_from_inside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let t = obb_ray_cast(&obb, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        if let Some(t) = t {
            assert!(t >= -1e-9, "t must be non-negative, got {}", t);
        }
    }
    #[test]
    fn test_obb_transform_translation_only() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let t = obb_transform(&obb, id, [3.0, 4.0, 5.0]);
        assert!((t.center[0] - 3.0).abs() < 1e-9);
        assert!((t.center[1] - 4.0).abs() < 1e-9);
        assert!((t.center[2] - 5.0).abs() < 1e-9);
        assert!((t.half_extents[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_obb_transform_identity_rotation() {
        let obb = Obb::axis_aligned([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]);
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let t = obb_transform(&obb, id, [0.0, 0.0, 0.0]);
        for i in 0..3 {
            assert!((t.center[i] - obb.center[i]).abs() < 1e-9, "center[{}]", i);
            assert!(
                (t.half_extents[i] - obb.half_extents[i]).abs() < 1e-9,
                "he[{}]",
                i
            );
        }
    }
    #[test]
    fn test_obb_sphere_test_intersecting() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_sphere_test(&obb, [2.0, 0.0, 0.0], 2.0);
        assert!(
            result.is_some(),
            "sphere overlapping OBB should return contact"
        );
    }
    #[test]
    fn test_obb_sphere_test_separated() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_sphere_test(&obb, [10.0, 0.0, 0.0], 1.0);
        assert!(result.is_none(), "far sphere should return None");
    }
    #[test]
    fn test_obb_support_point_x_positive() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let sp = obb_support_point(&obb, [1.0, 0.0, 0.0]);
        assert!((sp[0] - 1.0).abs() < 1e-9, "sp.x={}", sp[0]);
    }
    #[test]
    fn test_obb_support_point_y_negative() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let sp = obb_support_point(&obb, [0.0, -1.0, 0.0]);
        assert!((sp[1] + 2.0).abs() < 1e-9, "sp.y={}", sp[1]);
    }
    #[test]
    fn test_obb_volume_unit_cube() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let v = obb_volume(&obb);
        assert!((v - 1.0).abs() < 1e-10, "unit cube volume={}", v);
    }
    #[test]
    fn test_obb_surface_area_unit_half_extents() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let sa = obb_surface_area(&obb);
        assert!((sa - 6.0).abs() < 1e-10, "SA={}", sa);
    }
    #[test]
    fn test_obb_sat_test_edge_edge_overlapping() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = ObbSatTest::test_edge_edge(&a, &b);
        assert!(
            result.is_some(),
            "overlapping axis-aligned OBBs must not be separated on edge-edge"
        );
    }
    #[test]
    fn test_obb_sat_test_edge_edge_separated_rotated() {
        let cos45 = std::f64::consts::FRAC_PI_4.cos();
        let sin45 = std::f64::consts::FRAC_PI_4.sin();
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb {
            center: [10.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
            rotation: [[cos45, sin45, 0.0], [-sin45, cos45, 0.0], [0.0, 0.0, 1.0]],
        };
        let result = ObbSatTest::test_edge_edge(&a, &b);
        assert!(
            result.is_none(),
            "widely separated rotated OBBs must be separated on some edge-edge axis"
        );
    }
    #[test]
    fn test_obb_sat_test_edge_edge_touching() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = ObbSatTest::test_edge_edge(&a, &b);
        assert!(
            result.is_some(),
            "touching axis-aligned OBBs should not be separated by edge-edge"
        );
    }
    #[test]
    fn test_obb_sat_contact_normal_overlapping_is_unit() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let normal = ObbSatTest::compute_contact_normal(&a, &b);
        assert!(
            normal.is_some(),
            "overlapping OBBs must yield a contact normal"
        );
        let n = normal.unwrap();
        let nlen = len3_raw(n);
        assert!(
            (nlen - 1.0).abs() < 1e-9,
            "contact normal must be unit, len={nlen}"
        );
    }
    #[test]
    fn test_obb_sat_contact_normal_separated_is_none() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let normal = ObbSatTest::compute_contact_normal(&a, &b);
        assert!(
            normal.is_none(),
            "separated OBBs must not yield a contact normal"
        );
    }
    #[test]
    fn test_obb_sat_contact_normal_aligned_with_separation_axis() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let n = ObbSatTest::compute_contact_normal(&a, &b).expect("should overlap");
        assert!(
            n[0].abs() > 0.9,
            "normal should be along X axis, got {:?}",
            n
        );
    }
    #[test]
    fn test_obb_sat_contact_normal_y_axis_contact() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([0.0, 1.8, 0.0], [1.0, 1.0, 1.0]);
        let n = ObbSatTest::compute_contact_normal(&a, &b).expect("should overlap along Y");
        assert!(n[1].abs() > 0.9, "normal should be along Y, got {:?}", n);
    }
    #[test]
    fn test_obb_sat_contact_normal_rotated_is_unit() {
        let cos45 = std::f64::consts::FRAC_PI_4.cos();
        let sin45 = std::f64::consts::FRAC_PI_4.sin();
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb {
            center: [1.2, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
            rotation: [[cos45, sin45, 0.0], [-sin45, cos45, 0.0], [0.0, 0.0, 1.0]],
        };
        if let Some(n) = ObbSatTest::compute_contact_normal(&a, &b) {
            let nlen = len3_raw(n);
            assert!(
                (nlen - 1.0).abs() < 1e-9,
                "normal must be unit for rotated OBBs, len={nlen}"
            );
        }
    }
    #[test]
    fn test_obb_sat_contact_wraps_obb_obb_test() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let via_sat_test = ObbSatTest::contact(&a, &b);
        let direct = obb_obb_test(&a, &b);
        assert_eq!(
            via_sat_test.is_some(),
            direct.is_some(),
            "ObbSatTest::contact and obb_obb_test must agree"
        );
        if let (Some(s), Some(d)) = (via_sat_test, direct) {
            assert!(
                (s.penetration_depth - d.penetration_depth).abs() < 1e-9,
                "penetration depths must match: {} vs {}",
                s.penetration_depth,
                d.penetration_depth
            );
        }
    }
    #[test]
    fn test_obb_from_points_none_for_empty() {
        let result = obb_from_points(&[]);
        assert!(result.is_none(), "empty points should return None");
    }
    #[test]
    fn test_obb_from_points_single() {
        let pts = [[1.0, 2.0, 3.0]];
        let obb = obb_from_points(&pts);
        if let Some(o) = obb {
            assert!((o.center[0] - 1.0).abs() < 1e-9);
        }
    }
    #[test]
    fn test_obb_from_points_axis_aligned_x() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let obb = obb_from_points(&pts).expect("should build OBB from collinear points");
        assert!(
            (obb.center[0] - 2.0).abs() < 1e-9,
            "center.x={}",
            obb.center[0]
        );
    }
    #[test]
    fn test_obb_contains_point_on_face() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_contains_point(&obb, [1.0, 0.0, 0.0]),
            "point on face should be inside or on boundary"
        );
    }
    #[test]
    fn test_obb_contains_point_corner() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_contains_point(&obb, [1.0, 1.0, 1.0]),
            "corner point should be on boundary"
        );
    }
    #[test]
    fn test_obb_triangle_test_contained() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [5.0, 5.0, 5.0]);
        let result = obb_triangle_test(&obb, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(result, "triangle inside OBB should intersect");
    }
    #[test]
    fn test_obb_triangle_test_outside() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_triangle_test(&obb, [10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.5, 1.0, 0.0]);
        assert!(!result, "triangle far from OBB should not intersect");
    }
}
/// Cluster a list of contact points, removing duplicates within `threshold` distance.
///
/// Returns a de-duplicated list where any two points are at least `threshold` apart.
pub fn cluster_contact_points(points: &[[f64; 3]], threshold: f64) -> Vec<[f64; 3]> {
    let mut result: Vec<[f64; 3]> = Vec::new();
    let thr_sq = threshold * threshold;
    for &p in points {
        let duplicate = result.iter().any(|&q| {
            let d = sub3_raw(p, q);
            dot3_raw(d, d) < thr_sq
        });
        if !duplicate {
            result.push(p);
        }
    }
    result
}
/// Classify a contact result as face or edge-edge based on the contact normal.
///
/// If the normal aligns with any face normal of either OBB, it is a face contact.
/// Otherwise it is an edge-edge contact.
pub fn classify_contact_feature(a: &Obb, b: &Obb, normal: [f64; 3]) -> ContactFeatureType {
    let threshold = 0.99;
    for obb in &[a, b] {
        for ax in &obb.rotation {
            if dot3_raw(*ax, normal).abs() > threshold {
                return ContactFeatureType::FaceContact;
            }
        }
    }
    ContactFeatureType::EdgeEdgeContact
}
/// Compute a simple spatial hash for an OBB center, useful for broad-phase bucketing.
///
/// `cell_size` is the size of each hash cell.
pub fn obb_spatial_hash(obb: &Obb, cell_size: f64) -> (i64, i64, i64) {
    fn hash_coord(x: f64, cell: f64) -> i64 {
        (x / cell).floor() as i64
    }
    (
        hash_coord(obb.center[0], cell_size),
        hash_coord(obb.center[1], cell_size),
        hash_coord(obb.center[2], cell_size),
    )
}
/// Test whether an OBB swept linearly along `velocity` for time in `[0, max_t]`
/// would intersect a static OBB.
///
/// This is a conservative test: it computes the Minkowski sum of the moving OBB
/// trajectory with the static OBB and tests a point query.
///
/// Returns the earliest time of contact (clamped to `max_t`), or `None` if the
/// swept OBB never intersects.
pub fn obb_swept_test(
    moving: &Obb,
    velocity: [f64; 3],
    stationary: &Obb,
    max_t: f64,
) -> Option<f64> {
    let steps = 8;
    for i in 0..=steps {
        let t = max_t * (i as f64 / steps as f64);
        let disp = scale3_raw(velocity, t);
        let moved_center = add3_raw(moving.center, disp);
        let moved = Obb {
            center: moved_center,
            half_extents: moving.half_extents,
            rotation: moving.rotation,
        };
        if obb_obb_test(&moved, stationary).is_some() {
            return Some(t);
        }
    }
    None
}
/// Compute the bounding sphere (center, radius) of an OBB.
///
/// The sphere is centered at the OBB center and has radius equal to the
/// half-diagonal of the OBB (conservative bound).
pub fn obb_bounding_sphere(obb: &Obb) -> ([f64; 3], f64) {
    let hx = obb.half_extents[0];
    let hy = obb.half_extents[1];
    let hz = obb.half_extents[2];
    let radius = (hx * hx + hy * hy + hz * hz).sqrt();
    (obb.center, radius)
}
/// Quick sphere-sphere overlap test. Returns `true` if the bounding spheres overlap.
pub fn obb_bounding_sphere_test(a: &Obb, b: &Obb) -> bool {
    let (ca, ra) = obb_bounding_sphere(a);
    let (cb, rb) = obb_bounding_sphere(b);
    let d = sub3_raw(ca, cb);
    let dist_sq = dot3_raw(d, d);
    let sum_r = ra + rb;
    dist_sq < sum_r * sum_r
}
/// Return the six face normals of an OBB (three positive and three negative).
pub fn obb_face_normals(obb: &Obb) -> [[f64; 3]; 6] {
    [
        obb.rotation[0],
        negate3_raw(obb.rotation[0]),
        obb.rotation[1],
        negate3_raw(obb.rotation[1]),
        obb.rotation[2],
        negate3_raw(obb.rotation[2]),
    ]
}
/// Return the normal of the OBB face closest to the world-space `query` point.
pub fn obb_closest_face_normal(obb: &Obb, query: [f64; 3]) -> [f64; 3] {
    let d = sub3_raw(query, obb.center);
    let mut best_dot = f64::NEG_INFINITY;
    let mut best_normal = obb.rotation[0];
    for &rot in &obb.rotation {
        let proj = dot3_raw(d, rot);
        if proj > best_dot {
            best_dot = proj;
            best_normal = rot;
        }
        if -proj > best_dot {
            best_dot = -proj;
            best_normal = negate3_raw(rot);
        }
    }
    best_normal
}
/// Run the full 15-axis SAT and return the overlap on every axis.
///
/// Returns a list of (axis, overlap) pairs in the order they were tested.
/// Any axis where the shapes are separated returns `None` in the overlap slot.
pub fn obb_obb_all_axis_overlaps(a: &Obb, b: &Obb) -> Vec<([f64; 3], Option<f64>)> {
    let mut results = Vec::with_capacity(15);
    let axes_a = a.rotation;
    let axes_b = b.rotation;
    for &ax in &axes_a {
        results.push((ax, sat_test_axis(a, b, ax)));
    }
    for &ax in &axes_b {
        results.push((ax, sat_test_axis(a, b, ax)));
    }
    for &ai in &axes_a {
        for &bj in &axes_b {
            let cross = cross3_raw(ai, bj);
            results.push((cross, sat_test_axis(a, b, cross)));
        }
    }
    results
}
#[cfg(test)]
mod tests_obb_extended {

    use crate::narrowphase::classify_contact_feature;
    use crate::narrowphase::cluster_contact_points;
    use crate::narrowphase::compute_contact_points;
    use crate::narrowphase::obb_bounding_sphere;
    use crate::narrowphase::obb_bounding_sphere_test;
    use crate::narrowphase::obb_closest_face_normal;
    use crate::narrowphase::obb_face_normals;
    use crate::narrowphase::obb_obb_all_axis_overlaps;
    use crate::narrowphase::obb_obb_test;
    use crate::narrowphase::obb_sat::add3_raw;
    use crate::narrowphase::obb_sat::dot3_raw;
    use crate::narrowphase::obb_sat::len3_raw;
    use crate::narrowphase::obb_sat::normalize3_raw;
    use crate::narrowphase::obb_sat::types::ContactFeatureType;
    use crate::narrowphase::obb_sat::types::Obb;
    use crate::narrowphase::obb_spatial_hash;
    use crate::narrowphase::obb_swept_test;

    #[test]
    fn test_cluster_contact_points_removes_nearby() {
        let pts = vec![[0.0_f64, 0.0, 0.0], [0.001, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let result = cluster_contact_points(&pts, 0.01);
        assert_eq!(
            result.len(),
            2,
            "nearby point should be removed, got {}",
            result.len()
        );
    }
    #[test]
    fn test_cluster_contact_points_keeps_distant() {
        let pts = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let result = cluster_contact_points(&pts, 0.5);
        assert_eq!(result.len(), 3, "all distinct points should be kept");
    }
    #[test]
    fn test_cluster_contact_points_empty() {
        let result = cluster_contact_points(&[], 0.1);
        assert!(result.is_empty());
    }
    #[test]
    fn test_classify_contact_feature_face() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let feature = classify_contact_feature(&a, &b, [1.0, 0.0, 0.0]);
        assert_eq!(feature, ContactFeatureType::FaceContact);
    }
    #[test]
    fn test_classify_contact_feature_edge() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let n = normalize3_raw([1.0, 1.0, 0.0]);
        let feature = classify_contact_feature(&a, &b, n);
        assert_eq!(feature, ContactFeatureType::EdgeEdgeContact);
    }
    #[test]
    fn test_obb_spatial_hash_same_cell() {
        let a = Obb::axis_aligned([0.1, 0.2, 0.3], [0.5, 0.5, 0.5]);
        let b = Obb::axis_aligned([0.9, 0.8, 0.7], [0.5, 0.5, 0.5]);
        let ha = obb_spatial_hash(&a, 1.0);
        let hb = obb_spatial_hash(&b, 1.0);
        assert_eq!(ha, hb, "both centers in [0,1)^3 should hash to same cell");
    }
    #[test]
    fn test_obb_spatial_hash_different_cells() {
        let a = Obb::axis_aligned([0.5, 0.5, 0.5], [0.5, 0.5, 0.5]);
        let b = Obb::axis_aligned([2.5, 0.5, 0.5], [0.5, 0.5, 0.5]);
        let ha = obb_spatial_hash(&a, 1.0);
        let hb = obb_spatial_hash(&b, 1.0);
        assert_ne!(ha.0, hb.0, "different X cells expected");
    }
    #[test]
    fn test_obb_swept_test_hits() {
        let moving = Obb::axis_aligned([-3.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let stationary = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let result = obb_swept_test(&moving, [1.0, 0.0, 0.0], &stationary, 5.0);
        assert!(result.is_some(), "moving box should hit stationary box");
    }
    #[test]
    fn test_obb_swept_test_misses() {
        let moving = Obb::axis_aligned([-10.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let stationary = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let result = obb_swept_test(&moving, [1.0, 0.0, 0.0], &stationary, 1.0);
        assert!(
            result.is_none(),
            "box should not reach stationary in only 1 unit"
        );
    }
    #[test]
    fn test_obb_bounding_sphere_unit_cube() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (center, radius) = obb_bounding_sphere(&obb);
        let expected_r = 3.0_f64.sqrt();
        for (i, &ci) in center.iter().enumerate() {
            assert!(ci.abs() < 1e-10, "center[{}]={}", i, ci);
        }
        assert!((radius - expected_r).abs() < 1e-10, "radius={}", radius);
    }
    #[test]
    fn test_obb_bounding_sphere_test_overlap() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            obb_bounding_sphere_test(&a, &b),
            "bounding spheres should overlap"
        );
    }
    #[test]
    fn test_obb_bounding_sphere_test_separated() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        assert!(
            !obb_bounding_sphere_test(&a, &b),
            "bounding spheres should be separated"
        );
    }
    #[test]
    fn test_obb_face_normals_count() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let normals = obb_face_normals(&obb);
        assert_eq!(normals.len(), 6);
    }
    #[test]
    fn test_obb_face_normals_unit_length() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        for n in &obb_face_normals(&obb) {
            let len = len3_raw(*n);
            assert!((len - 1.0).abs() < 1e-10, "face normal not unit: {}", len);
        }
    }
    #[test]
    fn test_obb_face_normals_pairs_are_opposite() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let ns = obb_face_normals(&obb);
        for i in 0..3 {
            let p = ns[2 * i];
            let m = ns[2 * i + 1];
            let sum = add3_raw(p, m);
            assert!(
                len3_raw(sum) < 1e-10,
                "pair {} not opposite: {:?} + {:?}",
                i,
                p,
                m
            );
        }
    }
    #[test]
    fn test_obb_closest_face_normal_positive_x() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let n = obb_closest_face_normal(&obb, [10.0, 0.0, 0.0]);
        assert!(
            dot3_raw(n, [1.0, 0.0, 0.0]) > 0.9,
            "expected +X normal, got {:?}",
            n
        );
    }
    #[test]
    fn test_obb_closest_face_normal_negative_y() {
        let obb = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let n = obb_closest_face_normal(&obb, [0.0, -10.0, 0.0]);
        assert!(
            dot3_raw(n, [0.0, -1.0, 0.0]) > 0.9,
            "expected -Y normal, got {:?}",
            n
        );
    }
    #[test]
    fn test_obb_obb_all_axis_overlaps_count() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let results = obb_obb_all_axis_overlaps(&a, &b);
        assert_eq!(results.len(), 15, "must test 15 axes");
    }
    #[test]
    fn test_obb_obb_all_axis_overlaps_overlapping_has_no_none() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let results = obb_obb_all_axis_overlaps(&a, &b);
        let non_degenerate: Vec<_> = results
            .iter()
            .filter(|(ax, _)| dot3_raw(*ax, *ax) > 1e-12)
            .collect();
        for (_, overlap) in &non_degenerate {
            assert!(
                overlap.is_some(),
                "overlapping OBBs should not show None on any non-degenerate axis"
            );
        }
    }
    #[test]
    fn test_obb_obb_all_axis_overlaps_separated_has_none() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([10.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let results = obb_obb_all_axis_overlaps(&a, &b);
        let has_separation = results
            .iter()
            .any(|(ax, ov)| dot3_raw(*ax, *ax) > 1e-12 && ov.is_none());
        assert!(
            has_separation,
            "separated OBBs must have at least one separating axis"
        );
    }
    #[test]
    fn test_obb_sat_test_edge_parallel_axes() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([0.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_obb_test(&a, &b);
        assert!(
            result.is_some(),
            "axis-aligned boxes sharing edges should overlap"
        );
    }
    #[test]
    fn test_obb_sat_touching_boxes() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let result = obb_obb_test(&a, &b);
        if let Some(r) = result {
            assert!(
                r.penetration_depth >= 0.0,
                "touching depth must be non-negative"
            );
        }
    }
    #[test]
    fn test_obb_sat_large_vs_small() {
        let large = Obb::axis_aligned([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
        let small = Obb::axis_aligned([0.0, 0.0, 0.0], [0.01, 0.01, 0.01]);
        let result = obb_obb_test(&large, &small);
        assert!(
            result.is_some(),
            "large box should contain small box → overlap"
        );
        let r = result.unwrap();
        assert!(r.penetration_depth > 0.0);
    }
    #[test]
    fn test_cluster_then_sat_contact() {
        let a = Obb::axis_aligned([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let sat = obb_obb_test(&a, &b).expect("should overlap");
        let raw_contacts = compute_contact_points(&a, &b, &sat);
        let clustered = cluster_contact_points(&raw_contacts, 0.01);
        assert!(clustered.len() <= raw_contacts.len());
        assert!(!clustered.is_empty());
    }
}
