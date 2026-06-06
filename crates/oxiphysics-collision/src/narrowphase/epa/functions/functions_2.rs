//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::super::types::{EpaConfig, EpaFaceRaw, EpaPolytope, EpaStats};
    use super::super::functions_impl::{
        compute_face_normal_raw, epa_barycentric_origin, epa_contact_normal_from_face,
        epa_filter_degenerate_faces, epa_has_converged, epa_newell_normal, epa_penetration,
        epa_penetration_early_exit, epa_penetration_with_config, epa_penetration_with_stats,
        epa_recompute_normals, epa_silhouette, epa_solve, epa_triangle_area, epa_validate_polytope,
        epa_with_fallback, find_horizon,
    };
    use super::super::functions_impl::{
        epa_cross3, epa_dot3, epa_len3, epa_negate3, epa_scale3, epa_sub3,
    };
    use crate::Epa;
    use crate::narrowphase::EpaSolver;
    use crate::narrowphase::epa::epa_add3;
    use crate::narrowphase::epa_penetration_capped;
    use crate::narrowphase::gjk::{Gjk, GjkResult};
    use oxiphysics_core::Transform;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_geometry::Sphere;
    #[test]
    fn test_epa_sphere_sphere() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        if let GjkResult::Intersecting(simplex) = result {
            let contact = Epa::penetration_depth(&s1, &t1, &s2, &t2, &simplex);
            assert!(contact.is_some());
            let c = contact.unwrap();
            assert!((c.depth - 1.0).abs() < 0.2, "depth = {}", c.depth);
        } else {
            panic!("Expected intersection");
        }
    }
    /// Two overlapping unit-radius spheres with centers 1 m apart.
    /// Expected penetration depth = r1 + r2 - dist = 1 + 1 - 1 = 1 m.
    /// The contact normal produced by EPA must be a unit vector and the depth
    /// must be non-negative.
    ///
    /// Note: EPA expands the Minkowski-difference polytope outward from the
    /// origin, so the returned normal points along the minimum-penetration
    /// axis.  With A at (0,0,0) and B at (1,0,0) the axis is X; EPA may
    /// return either +X or −X depending on which face it finds first.  We
    /// therefore only check magnitude, not sign.
    #[test]
    fn test_epa_sphere_sphere_overlap() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        let simplex = match result {
            GjkResult::Intersecting(s) => s,
            GjkResult::Separated { .. } => panic!("expected intersection"),
        };
        let contact = Epa::penetration_depth(&s1, &t1, &s2, &t2, &simplex)
            .expect("EPA must return a contact");
        assert!(
            (contact.depth - 1.0).abs() < 0.2,
            "expected depth ≈ 1.0, got {}",
            contact.depth
        );
        assert!(contact.depth >= 0.0, "depth must be non-negative");
        let nlen = contact.normal.norm();
        assert!(
            (nlen - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            nlen
        );
        assert!(
            contact.normal.x.abs() > 0.9,
            "normal should be aligned with X axis, got normal = {:?}",
            contact.normal
        );
    }
    fn make_tetrahedron() -> [[f64; 3]; 4] {
        [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ]
    }
    #[test]
    fn test_epa_polytope_from_simplex() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        assert_eq!(polytope.vertices.len(), 4);
        assert_eq!(polytope.faces.len(), 4);
        for face in &polytope.faces {
            assert!(face.distance >= 0.0, "face distance must be non-negative");
            let nlen = epa_len3(face.normal);
            assert!((nlen - 1.0).abs() < 1e-6, "normal must be unit length");
        }
    }
    #[test]
    fn test_epa_find_closest_face() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let idx = polytope.find_closest_face();
        let closest_dist = polytope.faces[idx].distance;
        for (i, face) in polytope.faces.iter().enumerate() {
            if i != idx {
                assert!(face.distance >= closest_dist - 1e-10);
            }
        }
    }
    #[test]
    fn test_epa_expand() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_faces = polytope.faces.len();
        let expanded = polytope.expand([3.0, 0.0, 0.0]);
        assert!(expanded, "expansion should succeed");
        assert!(polytope.vertices.len() == 5);
        assert!(polytope.faces.len() >= initial_faces);
    }
    #[test]
    fn test_find_horizon() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let horizon = find_horizon(&polytope.faces, &polytope.vertices, [3.0, 0.0, 0.0]);
        assert!(!horizon.is_empty(), "horizon must have edges");
    }
    #[test]
    fn test_epa_penetration_spheres() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(result.is_some(), "EPA must find penetration");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0, "depth must be positive");
        let nlen = epa_len3(pen.normal);
        assert!(
            (nlen - 1.0).abs() < 1e-4,
            "normal must be unit length, got {}",
            nlen
        );
    }
    #[test]
    fn test_epa_penetration_depth_value() {
        let simplex: [[f64; 3]; 4] = [
            [0.3, 0.3, 0.3],
            [0.3, -0.3, -0.3],
            [-0.3, 0.3, -0.3],
            [-0.3, -0.3, 0.3],
        ];
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        let pen = result.expect("EPA must converge");
        assert!(
            (pen.depth - 1.0).abs() < 0.2,
            "depth should be ~1.0, got {}",
            pen.depth
        );
    }
    #[test]
    fn test_epa_face_normals_outward() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        for face in &polytope.faces {
            let v = polytope.vertices[face.vertices[0]];
            let d = epa_dot3(face.normal, v);
            assert!(d >= -1e-10, "face normal must point outward, dot = {}", d);
        }
    }
    #[test]
    fn test_epa_zero_iter() {
        let simplex = make_tetrahedron();
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0], 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            0,
            1e-6,
        );
        assert!(
            result.is_some(),
            "should return best result even with 0 iterations"
        );
    }
    #[test]
    fn test_epa_config_default() {
        let cfg = EpaConfig::default();
        assert_eq!(cfg.max_iter, 64);
        assert!(cfg.tolerance > 0.0);
    }
    #[test]
    fn test_epa_config_high_precision() {
        let cfg = EpaConfig::high_precision();
        assert!(cfg.max_iter > 64);
        assert!(cfg.tolerance < 1e-6);
    }
    #[test]
    fn test_epa_has_converged() {
        assert!(
            epa_has_converged(1.0, 1.0 + 1e-8, 1e-6),
            "tiny change should converge"
        );
        assert!(
            !epa_has_converged(1.0, 1.1, 1e-6),
            "large change should not converge"
        );
    }
    #[test]
    fn test_epa_barycentric_origin_equilateral() {
        let a = [1.0, 0.0, 0.0];
        let b = [-0.5, 0.866, 0.0];
        let c = [-0.5, -0.866, 0.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-10, "weights must sum to 1");
        assert!(
            u > 0.0 && v > 0.0 && w > 0.0,
            "all weights must be positive"
        );
    }
    #[test]
    fn test_epa_barycentric_sum_to_one() {
        let a = [2.0, 0.0, 0.0];
        let b = [0.0, 2.0, 0.0];
        let c = [0.0, 0.0, 2.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-10, "u+v+w={}", u + v + w);
    }
    #[test]
    fn test_epa_contact_normal_from_face_unit_length() {
        let face = EpaFaceRaw {
            vertices: [0, 1, 2],
            normal: [1.0, 0.0, 0.0],
            distance: 1.5,
        };
        let n = epa_contact_normal_from_face(&face);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_epa_validate_polytope_valid() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        assert!(
            epa_validate_polytope(&polytope),
            "initial polytope must be valid"
        );
    }
    #[test]
    fn test_epa_validate_polytope_invalid_negative_distance() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        polytope.faces[0].distance = -1.0;
        assert!(
            !epa_validate_polytope(&polytope),
            "negative distance should fail validation"
        );
    }
    #[test]
    fn test_epa_solve_sphere_cso() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (normal, depth) = epa_solve(
            simplex,
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            64,
        );
        assert!(depth > 0.0, "penetration depth must be positive");
        let nlen = epa_len3(normal);
        assert!((nlen - 1.0).abs() < 1e-4, "normal must be unit length");
        assert!(
            (depth - 1.0).abs() < 0.2,
            "expected depth ~1.0, got {}",
            depth
        );
    }
    #[test]
    fn test_compute_face_normal_outward() {
        let v0 = [1.0, 0.0, 0.0_f64];
        let v1 = [0.0, 1.0, 0.0_f64];
        let v2 = [-1.0, 0.0, 0.0_f64];
        let n = compute_face_normal_raw(v0, v1, v2);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-6, "normal must be unit length");
        assert!(
            n[2].abs() > 0.9,
            "normal should be along Z for XY-plane triangle, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_solve_expand_adds_face() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_face_count = polytope.faces.len();
        let expanded = polytope.expand([3.0, 0.0, 0.0]);
        assert!(expanded, "expansion with a far point should succeed");
        assert!(
            polytope.faces.len() >= initial_face_count,
            "after expansion, face count should not decrease"
        );
    }
    #[test]
    fn test_epa_penetration_with_config() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let cfg = EpaConfig::high_precision();
        let result = epa_penetration_with_config(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            &cfg,
        );
        assert!(
            result.is_some(),
            "EPA with high-precision config should succeed"
        );
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_newell_normal_xy_triangle() {
        let verts = vec![[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(
            (len - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            len
        );
        assert!(
            n[2].abs() > 0.9,
            "normal must be along Z for XY triangle, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_newell_normal_degenerate() {
        let n = epa_newell_normal(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(n, [0.0, 1.0, 0.0]);
    }
    #[test]
    fn test_epa_newell_normal_collinear() {
        let verts = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(len > 0.5, "normal length={}", len);
    }
    #[test]
    fn test_epa_newell_normal_quad_xz_plane() {
        let verts = vec![
            [1.0_f64, 0.0, 1.0],
            [-1.0, 0.0, 1.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
        ];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-6, "normal unit len={}", len);
        assert!(
            n[1].abs() > 0.9,
            "expected Y normal for XZ-plane quad, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_recompute_normals_preserves_validity() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        epa_recompute_normals(&mut polytope);
        for face in &polytope.faces {
            let v = polytope.vertices[face.vertices[0]];
            let d = epa_dot3(face.normal, v);
            assert!(
                d >= -1e-6,
                "recomputed distance must be non-negative, got {}",
                d
            );
        }
    }
    #[test]
    fn test_epa_triangle_area_unit() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = epa_triangle_area(a, b, c);
        assert!((area - 0.5).abs() < 1e-10, "area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_equilateral() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [1.0, 3.0_f64.sqrt(), 0.0];
        let area = epa_triangle_area(a, b, c);
        let expected = 3.0_f64.sqrt();
        assert!(
            (area - expected).abs() < 1e-6,
            "area={}, expected={}",
            area,
            expected
        );
    }
    #[test]
    fn test_epa_triangle_area_degenerate() {
        let p = [1.0_f64, 2.0, 3.0];
        let area = epa_triangle_area(p, p, p);
        assert!(area < 1e-10, "degenerate triangle area={}", area);
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_tiny() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        polytope.vertices.push([0.01, 0.0, 0.0]);
        polytope.vertices.push([0.01, 1e-8, 0.0]);
        polytope.vertices.push([0.01, 0.0, 1e-8]);
        let n = polytope.vertices.len();
        polytope.faces.push(EpaFaceRaw {
            vertices: [n - 3, n - 2, n - 1],
            normal: [1.0, 0.0, 0.0],
            distance: 0.01,
        });
        let count_before = polytope.faces.len();
        epa_filter_degenerate_faces(&mut polytope, 1e-5);
        let count_after = polytope.faces.len();
        assert!(
            count_after < count_before,
            "degenerate face should be removed"
        );
    }
    #[test]
    fn test_epa_filter_degenerate_faces_keeps_valid() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let count_before = polytope.faces.len();
        epa_filter_degenerate_faces(&mut polytope, 1e-12);
        assert_eq!(
            polytope.faces.len(),
            count_before,
            "no faces should be removed with tiny threshold"
        );
    }
    #[test]
    fn test_epa_stats_converged() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(result.is_some(), "EPA with stats must succeed");
        assert!(stats.iterations > 0, "must have at least 1 iteration");
        assert!(stats.final_depth > 0.0, "final depth must be positive");
        assert!(stats.face_count >= 4, "must have at least 4 faces");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_stats_iterations_bounded() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let max_iter = 10;
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            max_iter,
            1e-6,
        );
        assert!(
            stats.iterations <= max_iter,
            "iterations={} > max={}",
            stats.iterations,
            max_iter
        );
        assert!(result.is_some() || !stats.converged);
    }
    #[test]
    fn test_epa_silhouette_nonempty() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let new_point = [3.0, 0.0, 0.0];
        let (horizon, visible) = epa_silhouette(&polytope, new_point);
        assert!(
            !visible.is_empty(),
            "some faces must be visible from far point"
        );
        assert!(!horizon.is_empty(), "horizon must have edges");
    }
    #[test]
    fn test_epa_silhouette_no_visible_from_inside() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let (_horizon, visible) = epa_silhouette(&polytope, [0.0, 0.0, 0.0]);
        assert!(visible.is_empty(), "origin should see no faces as visible");
    }
    #[test]
    fn test_epa_penetration_early_exit_converges() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_early_exit(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-5,
        );
        assert!(result.is_some(), "early-exit EPA must return a result");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0, "depth must be positive");
        let nlen = epa_len3(pen.normal);
        assert!((nlen - 1.0).abs() < 1e-4, "normal must be unit length");
    }
    #[test]
    fn test_epa_with_fallback_sphere_support() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let pen = epa_with_fallback(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(pen.depth >= 0.0, "fallback depth must be non-negative");
        let nlen = epa_len3(pen.normal);
        assert!(
            (nlen - 1.0).abs() < 1e-4,
            "fallback normal must be unit length"
        );
    }
    #[test]
    fn test_epa_with_fallback_returns_positive_depth() {
        let simplex: [[f64; 3]; 4] = [
            [0.3, 0.3, 0.3],
            [0.3, -0.3, -0.3],
            [-0.3, 0.3, -0.3],
            [-0.3, -0.3, 0.3],
        ];
        let pen = epa_with_fallback(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(pen.depth > 0.0, "fallback depth={}", pen.depth);
    }
    #[test]
    fn test_epa_penetration_capped_respects_max_faces() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_capped(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
            10,
        );
        assert!(result.is_some(), "capped EPA must return a result");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_penetration_capped_4_face_limit() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_capped(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
            4,
        );
        assert!(result.is_some());
    }
    #[test]
    fn test_epa_add3_basic() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let r = epa_add3(a, b);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_epa_dot3_orthogonal() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(
            (epa_dot3(a, b)).abs() < 1e-15,
            "orthogonal vectors must have dot=0"
        );
    }
    #[test]
    fn test_epa_dot3_parallel() {
        let a = [3.0_f64, 0.0, 0.0];
        let b = [5.0, 0.0, 0.0];
        assert!((epa_dot3(a, b) - 15.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_sub3_basic() {
        let a = [5.0_f64, 3.0, 1.0];
        let b = [1.0, 2.0, 3.0];
        let r = epa_sub3(a, b);
        assert_eq!(r, [4.0, 1.0, -2.0]);
    }
    #[test]
    fn test_epa_scale3_halving() {
        let v = [2.0_f64, 4.0, 8.0];
        let r = epa_scale3(v, 0.5);
        assert_eq!(r, [1.0, 2.0, 4.0]);
    }
    #[test]
    fn test_epa_negate3_roundtrip() {
        let v = [1.0_f64, -2.0, 3.0];
        let neg = epa_negate3(v);
        let back = epa_negate3(neg);
        assert_eq!(back, v);
    }
    #[test]
    fn test_epa_cross3_anticommutative() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let ab = epa_cross3(a, b);
        let ba = epa_cross3(b, a);
        for i in 0..3 {
            assert!(
                (ab[i] + ba[i]).abs() < 1e-15,
                "cross must be anti-commutative: idx={}",
                i
            );
        }
    }
    #[test]
    fn test_epa_cross3_x_cross_y_is_z() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = epa_cross3(a, b);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_len3_pythagorean() {
        let v = [3.0_f64, 4.0, 0.0];
        assert!((epa_len3(v) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_epa_len3_unit_axes() {
        assert!((epa_len3([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-15);
        assert!((epa_len3([0.0, 1.0, 0.0]) - 1.0).abs() < 1e-15);
        assert!((epa_len3([0.0, 0.0, 1.0]) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_add3_commutative() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [7.0, -1.0, 0.5];
        let ab = epa_add3(a, b);
        let ba = epa_add3(b, a);
        assert_eq!(ab, ba);
    }
    #[test]
    fn test_epa_config_default_values() {
        let cfg = EpaConfig::default();
        assert_eq!(cfg.max_iter, 64);
        assert!((cfg.tolerance - 1e-6).abs() < 1e-15);
        assert_eq!(cfg.max_faces, 256);
    }
    #[test]
    fn test_epa_config_high_precision_tighter() {
        let hp = EpaConfig::high_precision();
        assert!(hp.max_iter > 64, "high_precision must have more iterations");
        assert!(
            hp.tolerance < 1e-6,
            "high_precision must have tighter tolerance"
        );
    }
    #[test]
    fn test_epa_has_converged_yes() {
        assert!(
            epa_has_converged(1.0, 1.0 + 1e-8, 1e-6),
            "change of 1e-8 within tol=1e-6 must converge"
        );
    }
    #[test]
    fn test_epa_has_converged_no() {
        assert!(
            !epa_has_converged(0.0, 0.01, 1e-6),
            "change of 0.01 outside tol=1e-6 must not converge"
        );
    }
    #[test]
    fn test_epa_has_converged_exact() {
        assert!(
            epa_has_converged(5.0, 5.0, 1e-15),
            "zero change within any positive tol must converge"
        );
    }
    #[test]
    fn test_compute_face_normal_raw_xy_plane() {
        let n = compute_face_normal_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            n[2].abs() > 0.9,
            "normal z-component must dominate: n={:?}",
            n
        );
    }
    #[test]
    fn test_compute_face_normal_raw_unit_length() {
        let n = compute_face_normal_raw([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let len = epa_len3(n);
        assert!(
            (len - 1.0).abs() < 1e-9,
            "face normal should be unit length: len={}",
            len
        );
    }
    #[test]
    fn test_epa_barycentric_origin_centroid() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [-0.5, 0.866_025, 0.0];
        let c = [-0.5, -0.866_025, 0.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!(
            (u + v + w - 1.0).abs() < 1e-9,
            "bary coords must sum to 1: {}",
            u + v + w
        );
        assert!(
            u >= 0.0 && v >= 0.0 && w >= 0.0,
            "bary coords must be non-negative"
        );
    }
    #[test]
    fn test_epa_barycentric_origin_sum_one() {
        let a = [2.0_f64, 0.0, 0.0];
        let b = [0.0, 2.0, 0.0];
        let c = [0.0, 0.0, 2.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-9, "sum={}", u + v + w);
    }
    #[test]
    fn test_epa_triangle_area_right_triangle() {
        let area = epa_triangle_area([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]);
        assert!((area - 6.0).abs() < 1e-9, "area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_collapsed_point() {
        let area = epa_triangle_area([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(area < 1e-12, "degenerate triangle area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_unit_equilateral() {
        let area = epa_triangle_area(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3_f64.sqrt() / 2.0, 0.0],
        );
        assert!((area - 3_f64.sqrt() / 4.0).abs() < 1e-9, "area={}", area);
    }
    #[test]
    fn test_epa_polytope_from_gjk_simplex_face_count() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        assert_eq!(poly.faces.len(), 4, "tetrahedron must have 4 faces");
        assert_eq!(poly.vertices.len(), 4, "tetrahedron must have 4 vertices");
    }
    #[test]
    fn test_epa_polytope_find_closest_face_non_panics_empty() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let idx = poly.find_closest_face();
        assert_eq!(idx, 0);
    }
    #[test]
    fn test_epa_polytope_closest_face_empty_default() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let (idx, dist, normal) = poly.closest_face();
        assert_eq!(idx, 0);
        assert!((dist - 0.0).abs() < 1e-15);
        let _ = normal;
    }
    #[test]
    fn test_epa_polytope_from_gjk_simplex_normals_unit() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        for (i, f) in poly.faces.iter().enumerate() {
            let len = epa_len3(f.normal);
            assert!(
                (len - 1.0).abs() < 1e-9,
                "face {} normal not unit: len={}",
                i,
                len
            );
        }
    }
    #[test]
    fn test_epa_validate_polytope_tetra_valid() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        let valid = epa_validate_polytope(&poly);
        assert!(valid, "well-formed polytope must validate");
    }
    #[test]
    fn test_epa_validate_polytope_no_panic_empty() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let _valid = epa_validate_polytope(&poly);
    }
    #[test]
    fn test_epa_newell_normal_square_xy_plane() {
        let verts = vec![
            [1.0_f64, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(len > 1e-9, "Newell normal must be non-zero for XY square");
        let normalized = epa_scale3(n, 1.0 / len);
        assert!(
            normalized[2].abs() > 0.9,
            "XY square Newell normal must be along Z: {:?}",
            normalized
        );
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_nothing_if_all_valid() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        let before = poly.faces.len();
        epa_filter_degenerate_faces(&mut poly, 0.0);
        assert_eq!(poly.faces.len(), before, "no faces should be filtered");
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_all_if_threshold_large() {
        let simplex: [[f64; 3]; 4] = [
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
            [-0.1, -0.1, -0.1],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        epa_filter_degenerate_faces(&mut poly, 1000.0);
        assert_eq!(
            poly.faces.len(),
            0,
            "all faces should be filtered with large threshold"
        );
    }
    #[test]
    fn test_epa_recompute_normals_no_panic() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        epa_recompute_normals(&mut poly);
        for (i, f) in poly.faces.iter().enumerate() {
            let len = epa_len3(f.normal);
            assert!(
                (len - 1.0).abs() < 1e-9,
                "face {} recomputed normal not unit: len={}",
                i,
                len
            );
        }
    }
    #[test]
    fn test_epa_stats_default_zero() {
        let s = EpaStats::default();
        assert_eq!(s.iterations, 0);
        assert_eq!(s.face_count, 0);
        assert!((s.final_depth - 0.0).abs() < 1e-15);
        assert!(!s.converged);
    }
    #[test]
    fn test_epa_penetration_with_stats_converged_flag() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                epa_scale3(n, 1.0)
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(stats.iterations > 0, "must have at least one iteration");
        if result.is_some() {
            assert!(stats.final_depth >= 0.0, "final depth must be non-negative");
        }
    }
    #[test]
    fn test_epa_penetration_with_stats_depth_positive() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, _stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0], 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        if let Some(pen) = result {
            assert!(
                pen.depth > 0.0,
                "penetration depth must be positive: {}",
                pen.depth
            );
            let nlen = epa_len3(pen.normal);
            assert!(
                (nlen - 1.0).abs() < 1e-6,
                "penetration normal must be unit: len={}",
                nlen
            );
        }
    }
    #[test]
    fn test_epa_contact_normal_from_face_passthrough() {
        let face = EpaFaceRaw {
            vertices: [0, 1, 2],
            normal: [0.0, 1.0, 0.0],
            distance: 0.5,
        };
        let n = epa_contact_normal_from_face(&face);
        assert_eq!(n, [0.0, 1.0, 0.0]);
    }
    #[test]
    fn test_epa_solver_penetration_axis_unit_normal() {
        let simplex = make_tetrahedron();
        let solver = EpaSolver::from_simplex(&simplex);
        let (normal, depth) = solver.compute_penetration_axis();
        let nlen = epa_len3(normal);
        assert!(
            (nlen - 1.0).abs() < 1e-9,
            "penetration axis must be unit, len={nlen}"
        );
        assert!(depth >= 0.0, "depth must be non-negative, got {depth}");
    }
    #[test]
    fn test_epa_solver_penetration_axis_depth_positive_for_tetrahedron() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let solver = EpaSolver::from_simplex(&simplex);
        let (_normal, depth) = solver.compute_penetration_axis();
        assert!(
            depth > 0.0,
            "depth must be > 0 for a simplex enclosing origin, got {depth}"
        );
    }
    #[test]
    fn test_epa_solver_penetration_axis_empty_polytope_fallback() {
        let simplex = make_tetrahedron();
        let mut solver = EpaSolver::from_simplex(&simplex);
        solver.polytope = Some(EpaPolytope {
            vertices: vec![],
            faces: vec![],
        });
        let (normal, depth) = solver.compute_penetration_axis();
        assert_eq!(normal, [0.0, 1.0, 0.0]);
        assert!((depth - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_solver_refine_polytope_increases_vertex_count() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let before = solver.vertex_count();
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        let expanded = solver.refine_polytope(&mut support_fn);
        if expanded {
            assert!(
                solver.vertex_count() > before,
                "refine should add a vertex when expanding"
            );
        }
    }
    #[test]
    fn test_epa_solver_refine_polytope_increments_iterations() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        assert_eq!(solver.iterations(), 0);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            epa_scale3(dir, 1.0 / len)
        };
        solver.refine_polytope(&mut support_fn);
        assert_eq!(
            solver.iterations(),
            1,
            "one call must increment iteration counter"
        );
    }
    #[test]
    fn test_epa_solver_refine_polytope_false_on_convergence() {
        let simplex: [[f64; 3]; 4] = [
            [0.001, 0.001, 0.001],
            [0.001, -0.001, -0.001],
            [-0.001, 0.001, -0.001],
            [-0.001, -0.001, 0.001],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex).with_tolerance(1.0);
        let mut support_fn = |_dir: [f64; 3]| -> [f64; 3] { [0.001, 0.0, 0.0] };
        let result = solver.refine_polytope(&mut support_fn);
        assert!(
            !result,
            "large tolerance should cause immediate convergence"
        );
    }
    #[test]
    fn test_epa_solver_solve_returns_some() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        let result = solver.solve(&mut support_fn, 64);
        assert!(result.is_some(), "solver must return a penetration result");
    }
    #[test]
    fn test_epa_solver_solve_unit_normal() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        if let Some(pen) = solver.solve(&mut support_fn, 64) {
            let nlen = epa_len3(pen.normal);
            assert!(
                (nlen - 1.0).abs() < 1e-6,
                "penetration normal must be unit, len={nlen}"
            );
        }
    }
    #[test]
    fn test_epa_solver_face_count_grows() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let initial_faces = solver.face_count();
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        for _ in 0..5 {
            solver.refine_polytope(&mut support_fn);
        }
        assert!(
            solver.face_count() >= initial_faces,
            "face count should not decrease: {} < {}",
            solver.face_count(),
            initial_faces
        );
    }
    #[test]
    fn test_epa_solver_with_tolerance_setter() {
        let simplex = make_tetrahedron();
        let solver = EpaSolver::from_simplex(&simplex).with_tolerance(1e-9);
        let (normal, _depth) = solver.compute_penetration_axis();
        let nlen = epa_len3(normal);
        assert!((nlen - 1.0).abs() < 1e-9);
    }
}
