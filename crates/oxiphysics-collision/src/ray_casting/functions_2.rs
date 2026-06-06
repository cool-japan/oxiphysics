//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::ray_casting::types::{
        Aabb, Bvh, Capsule, ConvexMesh, Heightfield, Obb, Ray, RayDifferential, RayTree, Scene,
        Sphere,
    };
    use std::f64::consts::PI;
    #[test]
    fn test_ray_at_origin() {
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = r.at(0.0);
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_ray_at_t() {
        let r = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = r.at(3.0);
        assert!((p[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_ray_normalised_direction() {
        let r = Ray::new_normalised([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        let l = len(r.direction);
        assert!((l - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_sphere_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let sphere = Sphere::new([0.0, 0.0, 0.0], 1.0);
        let hit = ray_sphere(&ray, &sphere);
        assert!(hit.is_some(), "ray should hit sphere");
        let t = hit.unwrap().t;
        assert!((t - 4.0).abs() < 1e-10, "t should be 4.0, got {}", t);
    }
    #[test]
    fn test_ray_sphere_miss() {
        let ray = Ray::new([0.0, 5.0, -5.0], [0.0, 0.0, 1.0]);
        let sphere = Sphere::new([0.0, 0.0, 0.0], 1.0);
        assert!(ray_sphere(&ray, &sphere).is_none());
    }
    #[test]
    fn test_ray_sphere_normal_unit_length() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let sphere = Sphere::new([0.0, 0.0, 0.0], 1.0);
        let hit = ray_sphere(&ray, &sphere).unwrap();
        let l = len(hit.normal);
        assert!((l - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_ray_sphere_hit_point_on_surface() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let sphere = Sphere::new([0.0, 0.0, 0.0], 2.0);
        let hit = ray_sphere(&ray, &sphere).unwrap();
        let d = (len(sub(hit.point, sphere.centre)) - sphere.radius).abs();
        assert!(d < 1e-8, "hit point should lie on sphere surface");
    }
    #[test]
    fn test_ray_sphere_inside() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let sphere = Sphere::new([0.0, 0.0, 0.0], 2.0);
        let hit = ray_sphere(&ray, &sphere);
        assert!(hit.is_some());
    }
    #[test]
    fn test_ray_aabb_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let result = ray_aabb(&ray, &aabb);
        assert!(result.is_some());
    }
    #[test]
    fn test_ray_aabb_miss() {
        let ray = Ray::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let result = ray_aabb(&ray, &aabb);
        assert!(result.is_none());
    }
    #[test]
    fn test_ray_aabb_hit_entry_exit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb::new([-1.0, -1.0, 2.0], [1.0, 1.0, 4.0]);
        let (t_enter, t_exit) = ray_aabb(&ray, &aabb).unwrap();
        assert!((t_enter - 7.0).abs() < 1e-10, "entry t should be 7.0");
        assert!((t_exit - 9.0).abs() < 1e-10, "exit t should be 9.0");
    }
    #[test]
    fn test_aabb_surface_area() {
        let aabb = Aabb::new([0.0; 3], [1.0, 2.0, 3.0]);
        let sa = aabb.surface_area();
        assert!((sa - 22.0).abs() < 1e-10, "surface area should be 22");
    }
    #[test]
    fn test_ray_obb_identity_matches_aabb() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let aabb = Aabb::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let obb = Obb::from_aabb(&aabb);
        let hit = ray_obb(&ray, &obb);
        assert!(hit.is_some(), "ray should hit OBB aligned with AABB");
    }
    #[test]
    fn test_ray_obb_normal_unit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let obb = Obb::new(
            [0.0, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [1.0, 1.0, 1.0],
        );
        let hit = ray_obb(&ray, &obb).unwrap();
        let l = len(hit.normal);
        assert!((l - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_ray_triangle_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let hit = ray_triangle(
            &ray,
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0,
        );
        assert!(hit.is_some(), "ray should hit triangle");
        assert!((hit.unwrap().t - 5.0).abs() < 1e-8);
    }
    #[test]
    fn test_ray_triangle_miss_outside() {
        let ray = Ray::new([10.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let hit = ray_triangle(
            &ray,
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0,
        );
        assert!(hit.is_none());
    }
    #[test]
    fn test_ray_triangle_parallel_miss() {
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hit = ray_triangle(&ray, [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 0);
        let _ = hit;
    }
    #[test]
    fn test_ray_triangle_barycentric_coords() {
        let ray = Ray::new([0.333, 0.333, -5.0], [0.0, 0.0, 1.0]);
        let hit = ray_triangle(&ray, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0).unwrap();
        let u = hit.uv[0];
        let v = hit.uv[1];
        assert!(
            u >= 0.0 && v >= 0.0 && u + v <= 1.0 + 1e-8,
            "barycentric coords out of range"
        );
    }
    #[test]
    fn test_ray_capsule_cylinder_body_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cap = Capsule::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.5);
        let hit = ray_capsule(&ray, &cap);
        assert!(hit.is_some(), "ray should hit capsule");
    }
    #[test]
    fn test_ray_capsule_miss() {
        let ray = Ray::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cap = Capsule::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.5);
        let hit = ray_capsule(&ray, &cap);
        assert!(hit.is_none());
    }
    #[test]
    fn test_ray_capsule_normal_unit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let cap = Capsule::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.5);
        let hit = ray_capsule(&ray, &cap).unwrap();
        let l = len(hit.normal);
        assert!(
            (l - 1.0).abs() < 1e-7,
            "capsule normal should be unit length"
        );
    }
    #[test]
    fn test_ray_convex_mesh_tetrahedron() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let tris = vec![(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)];
        let mesh = ConvexMesh::new(verts);
        let ray = Ray::new([0.5, 0.3, -5.0], [0.0, 0.0, 1.0]);
        let hit = ray_convex_mesh(&ray, &mesh, &tris);
        assert!(hit.is_some(), "ray should hit tetrahedron");
    }
    #[test]
    fn test_ray_convex_mesh_miss() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let tris = vec![(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)];
        let mesh = ConvexMesh::new(verts);
        let ray = Ray::new([10.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(ray_convex_mesh(&ray, &mesh, &tris).is_none());
    }
    #[test]
    fn test_ray_heightfield_flat() {
        let heights = vec![0.0f64; 25];
        let hf = Heightfield::new(heights, 5, 5, 1.0);
        let ray = Ray::new([2.0, 5.0, 2.0], [0.0, -1.0, 0.0]);
        let hit = ray_heightfield(&ray, &hf, 100);
        assert!(
            hit.is_some(),
            "ray pointing down should hit flat heightfield"
        );
    }
    #[test]
    fn test_ray_heightfield_above_always_miss() {
        let heights = vec![-100.0f64; 25];
        let hf = Heightfield::new(heights, 5, 5, 1.0);
        let ray = Ray {
            origin: [2.0, 5.0, 2.0],
            direction: [0.0, 1.0, 0.0],
            t_min: 0.0,
            t_max: 1000.0,
        };
        let hit = ray_heightfield(&ray, &hf, 200);
        assert!(
            hit.is_none(),
            "ray going up should not hit terrain at y=-100"
        );
    }
    #[test]
    fn test_batch_ray_spheres_all_hit() {
        let spheres = vec![Sphere::new([0.0, 0.0, 0.0], 1.0)];
        let rays: Vec<Ray> = (0..4)
            .map(|i| Ray::new([i as f64 * 0.1, 0.0, -5.0], [0.0, 0.0, 1.0]))
            .collect();
        let result = batch_ray_spheres(&rays, &spheres);
        for h in &result.hits {
            assert!(h.is_some());
        }
    }
    #[test]
    fn test_batch_ray_triangles_hit_count() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let indices = vec![(0, 1, 2)];
        let rays = vec![
            Ray::new([0.3, 0.2, -5.0], [0.0, 0.0, 1.0]),
            Ray::new([5.0, 5.0, -5.0], [0.0, 0.0, 1.0]),
        ];
        let result = batch_ray_triangles(&rays, &verts, &indices);
        assert!(result.hits[0].is_some(), "first ray should hit");
        assert!(result.hits[1].is_none(), "second ray should miss");
    }
    #[test]
    fn test_bvh_hit() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let indices = vec![(0, 1, 2)];
        let bvh = Bvh::build(verts, indices);
        let ray = Ray::new([0.4, 0.2, -5.0], [0.0, 0.0, 1.0]);
        let hit = bvh.intersect(&ray);
        assert!(hit.is_some(), "BVH should hit the triangle");
    }
    #[test]
    fn test_bvh_miss() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let indices = vec![(0, 1, 2)];
        let bvh = Bvh::build(verts, indices);
        let ray = Ray::new([10.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(bvh.intersect(&ray).is_none());
    }
    #[test]
    fn test_bvh_multi_triangles() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [3.5, 1.0, 0.0],
        ];
        let indices = vec![(0, 1, 2), (3, 4, 5)];
        let bvh = Bvh::build(verts, indices);
        let r1 = Ray::new([0.4, 0.2, -5.0], [0.0, 0.0, 1.0]);
        let r2 = Ray::new([3.4, 0.2, -5.0], [0.0, 0.0, 1.0]);
        assert!(bvh.intersect(&r1).is_some());
        assert!(bvh.intersect(&r2).is_some());
    }
    #[test]
    fn test_pick_ray_centre_aligned() {
        let ray = pick_ray(
            400.0,
            300.0,
            800.0,
            600.0,
            PI / 3.0,
            0.1,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
        );
        assert!(
            ray.direction[2] < -0.9,
            "centre ray should point mostly in -Z"
        );
    }
    #[test]
    fn test_pick_ray_normalised() {
        let ray = pick_ray(
            100.0,
            200.0,
            800.0,
            600.0,
            PI / 3.0,
            0.1,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
        );
        let l = len(ray.direction);
        assert!((l - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_shadow_ray_occluded() {
        let occluder = Sphere::new([0.0, 0.0, 5.0], 1.0);
        let in_shadow = is_in_shadow([0.0, 0.0, 0.0], [0.0, 0.0, 10.0], &[occluder]);
        assert!(in_shadow, "sphere should cast shadow");
    }
    #[test]
    fn test_shadow_ray_unoccluded() {
        let occluder = Sphere::new([5.0, 0.0, 0.0], 1.0);
        let in_shadow = is_in_shadow([0.0, 0.0, 0.0], [0.0, 0.0, 10.0], &[occluder]);
        assert!(!in_shadow, "off-axis sphere should not block ray");
    }
    #[test]
    fn test_soft_shadow_fully_lit() {
        let factor = soft_shadow_factor([0.0, 0.0, 0.0], [0.0, 10.0, 0.0], 1.0, &[], 16, 42);
        assert!((factor - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_reflect_perpendicular() {
        let d = [0.0, -1.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let r = reflect_ray(d, n);
        assert!(
            (r[1] - 1.0).abs() < 1e-8,
            "reflection of downward ray off upward normal is upward"
        );
    }
    #[test]
    fn test_refract_no_tir() {
        let d = normalize([0.0, -1.0, 0.0]);
        let n = [0.0, 1.0, 0.0];
        let r = refract_ray(d, n, 1.0, 1.5);
        assert!(r.is_some(), "should not have TIR going from air to glass");
    }
    #[test]
    fn test_refract_tir() {
        let d = normalize([0.99, -0.14, 0.0]);
        let n = [0.0, 1.0, 0.0];
        let r = refract_ray(d, n, 1.5, 1.0);
        let _ = r;
    }
    #[test]
    fn test_schlick_zero_at_normal() {
        let r = schlick_reflectance(1.0, 1.0, 1.0);
        assert!(r.abs() < 1e-10, "no reflectance for identical media");
    }
    #[test]
    fn test_ray_point_closest_t() {
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let q = [3.0, 1.0, 0.0];
        let (t, p) = ray_closest_point(&ray, q);
        assert!((t - 3.0).abs() < 1e-10);
        assert!((p[0] - 3.0).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
    }
    #[test]
    fn test_ray_point_dist_sq_origin() {
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let d = ray_point_dist_sq(&ray, [0.0, 3.0, 0.0]);
        assert!((d - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_segment_point_distance_midpoint() {
        let (t, p, d) = segment_point_distance([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]);
        assert!((t - 0.5).abs() < 1e-10);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((d - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_ray_distance_parallel() {
        let r1 = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let r2 = Ray::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        let (_t1, _t2, d) = ray_ray_distance(&r1, &r2);
        assert!(
            (d - 1.0).abs() < 1e-8,
            "parallel rays distance = 1.0, got {}",
            d
        );
    }
    #[test]
    fn test_ray_differential_footprint_positive() {
        let r = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let rx = Ray::new([0.01, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let ry = Ray::new([0.0, 0.01, -5.0], [0.0, 0.0, 1.0]);
        let diff = RayDifferential::new(r, rx, ry);
        let fp = diff.footprint_at(5.0);
        assert!(fp >= 0.0, "footprint should be non-negative");
    }
    #[test]
    fn test_ray_tree_primary_hit() {
        let mut tree = RayTree::new(3);
        let spheres = vec![Sphere::new([0.0, 0.0, 0.0], 1.0)];
        let primary = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        tree.trace(primary, &spheres);
        assert!(!tree.nodes.is_empty());
        assert!(tree.nodes[0].hit.is_some(), "primary ray should hit sphere");
    }
    #[test]
    fn test_ray_tree_no_hit_empty_scene() {
        let mut tree = RayTree::new(3);
        let primary = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        tree.trace(primary, &[]);
        assert!(tree.nodes[0].hit.is_none());
    }
    #[test]
    fn test_scene_sphere_hit() {
        let mut scene = Scene::new();
        scene.add_sphere(Sphere::new([0.0, 0.0, 0.0], 1.0));
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(scene.cast(&ray).is_some());
    }
    #[test]
    fn test_scene_triangle_hit() {
        let mut scene = Scene::new();
        scene.add_triangle([-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]);
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(scene.cast(&ray).is_some());
    }
    #[test]
    fn test_scene_closest_primitive() {
        let mut scene = Scene::new();
        scene.add_sphere(Sphere::new([0.0, 0.0, 0.0], 0.5));
        scene.add_sphere(Sphere::new([0.0, 0.0, 5.0], 0.5));
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let hit = scene.cast(&ray).unwrap();
        assert!(hit.t < 9.0, "should hit near sphere first, got t={}", hit.t);
    }
}
