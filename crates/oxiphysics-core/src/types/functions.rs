//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use crate::math::{Mat3, Quat, Real, Vec3};

#[cfg(test)]
mod prop_tests {
    use crate::Aabb;
    use crate::Quat;
    use crate::Transform;
    use crate::Vec3;
    use proptest::prelude::*;
    fn coord_strategy() -> impl Strategy<Value = f64> {
        -100.0_f64..100.0_f64
    }
    fn vec3_strategy() -> impl Strategy<Value = Vec3> {
        (coord_strategy(), coord_strategy(), coord_strategy())
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }
    fn positive_coord_strategy() -> impl Strategy<Value = f64> {
        0.01_f64..100.0_f64
    }
    fn positive_vec3_strategy() -> impl Strategy<Value = Vec3> {
        (
            positive_coord_strategy(),
            positive_coord_strategy(),
            positive_coord_strategy(),
        )
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }
    fn transform_strategy() -> impl Strategy<Value = Transform> {
        (
            vec3_strategy(),
            coord_strategy(),
            coord_strategy(),
            coord_strategy(),
            coord_strategy(),
        )
            .prop_map(|(pos, qi, qj, qk, qw)| {
                let raw = nalgebra::Quaternion::new(qw, qi, qj, qk);
                let norm = raw.norm();
                let quat = if norm < 1e-10 {
                    Quat::identity()
                } else {
                    Quat::from_quaternion(nalgebra::Quaternion::new(
                        raw.w / norm,
                        raw.i / norm,
                        raw.j / norm,
                        raw.k / norm,
                    ))
                };
                Transform::new(pos, quat)
            })
    }
    proptest! {
        #[test] fn prop_transform_compose_associative(ta in transform_strategy(), tb in
        transform_strategy(), tc in transform_strategy(),) { let lhs = ta.compose(& tb)
        .compose(& tc); let rhs = ta.compose(& tb.compose(& tc)); let pos_diff = (lhs
        .position - rhs.position).norm(); prop_assert!(pos_diff < 1e-6,
        "positions differ by {}", pos_diff); } #[test] fn prop_transform_inverse(ta in
        transform_strategy()) { let inv = ta.inverse(); let composed = ta.compose(& inv);
        let pos_diff = composed.position.norm(); prop_assert!(pos_diff < 1e-6,
        "T * T^-1 position not zero: norm={}", pos_diff); let rot_diff = (composed
        .rotation * Quat::identity().inverse()).angle(); prop_assert!(rot_diff < 1e-6,
        "T * T^-1 rotation not identity: angle={}", rot_diff); } #[test] fn
        prop_aabb_merge_contains_both(min_a in vec3_strategy(), ext_a in
        positive_vec3_strategy(), min_b in vec3_strategy(), ext_b in
        positive_vec3_strategy(),) { let a = Aabb::new(min_a, min_a + ext_a); let b =
        Aabb::new(min_b, min_b + ext_b); let merged = a.merge(& b); prop_assert!(merged
        .contains_point(& a.min), "merged does not contain a.min"); prop_assert!(merged
        .contains_point(& a.max), "merged does not contain a.max"); prop_assert!(merged
        .contains_point(& b.min), "merged does not contain b.min"); prop_assert!(merged
        .contains_point(& b.max), "merged does not contain b.max"); } #[test] fn
        prop_aabb_surface_area_positive(min in vec3_strategy(), ext in
        positive_vec3_strategy(),) { let aabb = Aabb::new(min, min + ext); let sa = aabb
        .surface_area(); prop_assert!(sa > 0.0, "surface area not positive: {}", sa); }
        #[test] fn prop_vec3_dot_commutative(a in vec3_strategy(), b in vec3_strategy(),)
        { let ab = a.dot(& b); let ba = b.dot(& a); prop_assert!((ab - ba).abs() < 1e-6,
        "dot not commutative: {} vs {}", ab, ba); } #[test] fn
        prop_vec3_cross_anti_commutative(a in vec3_strategy(), b in vec3_strategy(),) {
        let ab = a.cross(& b); let ba = b.cross(& a); let diff = (ab + ba).norm();
        prop_assert!(diff < 1e-6, "cross not anti-commutative: diff norm={}", diff); }
    }
}
#[cfg(test)]
mod tests {
    use crate::Aabb;
    use crate::BodyHandle;
    use crate::ColliderHandle;
    use crate::DefaultConfigs;
    use crate::MassProperties;
    use crate::PhysicsConfig;
    use crate::PhysicsConfigBuilder;
    use crate::Quat;
    use crate::TimeStep;
    use crate::Transform;
    use crate::Vec3;
    #[test]
    fn test_aabb_intersection() {
        let a = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(0.5, 0.5, 0.5), Vec3::new(1.5, 1.5, 1.5));
        let c = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        assert!(a.intersects(&b));
        assert!(!a.intersects(&c));
    }
    #[test]
    fn test_aabb_merge() {
        let a = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        let merged = a.merge(&b);
        assert_eq!(merged.min, Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(merged.max, Vec3::new(3.0, 3.0, 3.0));
    }
    #[test]
    fn test_aabb_contains_point() {
        let aabb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        assert!(aabb.contains_point(&Vec3::new(0.5, 0.5, 0.5)));
        assert!(!aabb.contains_point(&Vec3::new(2.0, 0.5, 0.5)));
    }
    #[test]
    fn test_aabb_surface_area() {
        let aabb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 3.0));
        let sa = aabb.surface_area();
        assert!((sa - 22.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_compose() {
        let t1 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.0, 2.0, 0.0));
        let composed = t1.compose(&t2);
        let p = composed.transform_point(&Vec3::zeros());
        assert!((p.x - 1.0).abs() < 1e-10);
        assert!((p.y - 2.0).abs() < 1e-10);
        assert!(p.z.abs() < 1e-10);
    }
    #[test]
    fn test_transform_inverse() {
        let t = Transform::from_position(Vec3::new(1.0, 2.0, 3.0));
        let inv = t.inverse();
        let composed = t.compose(&inv);
        assert!(composed.position.norm() < 1e-10);
        let angle = (composed.rotation * Quat::identity().inverse()).angle();
        assert!(angle < 1e-10);
    }
    #[test]
    fn test_transform_point() {
        let t = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let p = t.transform_point(&Vec3::new(1.0, 0.0, 0.0));
        assert!((p.x - 11.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_vector() {
        let t = Transform::from_position(Vec3::new(10.0, 5.0, 3.0));
        let v = Vec3::new(1.0, 0.0, 0.0);
        let tv = t.transform_vector(&v);
        assert!((tv.x - 1.0).abs() < 1e-10);
        assert!(tv.y.abs() < 1e-10);
        assert!(tv.z.abs() < 1e-10);
    }
    #[test]
    fn test_transform_default() {
        let t = Transform::default();
        assert_eq!(t.position, Vec3::zeros());
    }
    #[test]
    fn test_mass_properties() {
        let mp = MassProperties::point_mass(2.0);
        assert!((mp.inverse_mass() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_body_handle() {
        let h = BodyHandle::new(0, 0);
        assert_eq!(h.index, 0);
        assert_eq!(h.generation, 0);
    }
    #[test]
    fn test_physics_config_default() {
        let cfg = PhysicsConfig::default();
        assert!((cfg.gravity.y + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_transform_to_matrix4_identity() {
        let t = Transform::default();
        let m = t.to_matrix4();
        assert!((m[0][0] - 1.0).abs() < 1e-10);
        assert!((m[1][1] - 1.0).abs() < 1e-10);
        assert!((m[2][2] - 1.0).abs() < 1e-10);
        assert!((m[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_euler_angles_identity() {
        let t = Transform::default();
        let (roll, pitch, yaw) = t.euler_angles();
        assert!(roll.abs() < 1e-10);
        assert!(pitch.abs() < 1e-10);
        assert!(yaw.abs() < 1e-10);
    }
    #[test]
    fn test_transform_from_axis_angle() {
        let t = Transform::from_axis_angle(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::PI / 2.0,
        );
        let v = t.transform_vector(&Vec3::new(1.0, 0.0, 0.0));
        assert!((v.y - 1.0).abs() < 1e-10, "v = {:?}", v);
    }
    #[test]
    fn test_transform_lerp() {
        let t1 = Transform::from_position(Vec3::zeros());
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let mid = t1.lerp(&t2, 0.5);
        assert!((mid.position.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_from_center_half_extents() {
        let aabb =
            Aabb::from_center_half_extents(Vec3::new(5.0, 5.0, 5.0), Vec3::new(1.0, 2.0, 3.0));
        assert!((aabb.min.x - 4.0).abs() < 1e-10);
        assert!((aabb.max.x - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_from_points() {
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, -2.0, -3.0),
        ];
        let aabb = Aabb::from_points(&points).unwrap();
        assert!((aabb.min.x + 1.0).abs() < 1e-10);
        assert!((aabb.max.z - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_from_points_empty() {
        let result = Aabb::from_points(&[]);
        assert!(result.is_none());
    }
    #[test]
    fn test_aabb_longest_axis() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 3.0, 2.0));
        assert_eq!(aabb.longest_axis(), 1);
    }
    #[test]
    fn test_aabb_diagonal() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(3.0, 4.0, 0.0));
        assert!((aabb.diagonal() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_physics_config_validate_valid() {
        let cfg = PhysicsConfig::default();
        assert!(cfg.validate().is_ok());
    }
    #[test]
    fn test_physics_config_validate_invalid() {
        let cfg = PhysicsConfig {
            solver_iterations: 0,
            ..PhysicsConfig::default()
        };
        assert!(cfg.validate().is_err());
    }
    #[test]
    fn test_timestep_validate() {
        let ts = TimeStep::new(0.01);
        assert!(ts.validate().is_ok());
        let ts_bad = TimeStep::new(-0.01);
        assert!(ts_bad.validate().is_err());
    }
    #[test]
    fn test_timestep_frequency() {
        let ts = TimeStep::new(0.01);
        assert!((ts.frequency() - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_builder_pattern() {
        let cfg = PhysicsConfigBuilder::new()
            .gravity(Vec3::new(0.0, -1.625, 0.0))
            .solver_iterations(16)
            .ccd_enabled(false)
            .build();
        assert!((cfg.gravity.y + 1.625).abs() < 1e-10);
        assert_eq!(cfg.solver_iterations, 16);
        assert!(!cfg.ccd_enabled);
    }
    #[test]
    fn test_builder_default() {
        let builder = PhysicsConfigBuilder::default();
        let cfg = builder.build();
        assert!((cfg.gravity.y + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_default_configs_rigid_body() {
        let cfg = DefaultConfigs::rigid_body();
        assert!(cfg.ccd_enabled);
    }
    #[test]
    fn test_default_configs_soft_body() {
        let cfg = DefaultConfigs::soft_body();
        assert_eq!(cfg.solver_iterations, 20);
    }
    #[test]
    fn test_default_configs_fluid() {
        let cfg = DefaultConfigs::fluid();
        assert_eq!(cfg.solver_iterations, 50);
        assert!(cfg.gravity.norm() < 1e-10);
    }
    #[test]
    fn test_default_configs_space() {
        let cfg = DefaultConfigs::space();
        assert!(cfg.gravity.norm() < 1e-10);
        assert!(cfg.ccd_enabled);
    }
    #[test]
    fn test_default_configs_moon() {
        let cfg = DefaultConfigs::moon();
        assert!((cfg.gravity.y + 1.625).abs() < 1e-10);
    }
    #[test]
    fn test_default_configs_mars() {
        let cfg = DefaultConfigs::mars();
        assert!((cfg.gravity.y + 3.72).abs() < 1e-10);
    }
    #[test]
    fn test_mass_properties_sphere() {
        let mp = MassProperties::sphere(10.0, 1.0);
        assert!((mp.mass - 10.0).abs() < 1e-10);
        let i = mp.local_inertia[(0, 0)];
        assert!((i - 4.0).abs() < 1e-10, "I_sphere = {i}");
    }
    #[test]
    fn test_mass_properties_cuboid() {
        let mp = MassProperties::cuboid(12.0, 1.0, 2.0, 3.0);
        assert!((mp.mass - 12.0).abs() < 1e-10);
        let ixx = mp.local_inertia[(0, 0)];
        assert!((ixx - 13.0).abs() < 1e-10, "Ixx = {ixx}");
    }
    #[test]
    fn test_mass_properties_cylinder() {
        let mp = MassProperties::cylinder(6.0, 1.0, 2.0);
        assert!((mp.mass - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_mass_properties_is_static() {
        let mp = MassProperties::point_mass(0.0);
        assert!(mp.is_static());
        let mp2 = MassProperties::point_mass(1.0);
        assert!(!mp2.is_static());
    }
    #[test]
    fn test_mass_properties_combine() {
        let mp1 = MassProperties::point_mass(2.0);
        let mp2 = MassProperties::point_mass(3.0);
        let combined = mp1.combine(&mp2);
        assert!((combined.mass - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_body_handle_null() {
        let h = BodyHandle::null();
        assert!(!h.is_valid());
        let h2 = BodyHandle::new(1, 0);
        assert!(h2.is_valid());
    }
    #[test]
    fn test_collider_handle_null() {
        let h = ColliderHandle::null();
        assert!(!h.is_valid());
    }
}
/// Multiply two quaternions `a * b` (Hamilton product), both stored as `[x,y,z,w]`.
pub(super) fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let (ax, ay, az, aw) = (a[0], a[1], a[2], a[3]);
    let (bx, by, bz, bw) = (b[0], b[1], b[2], b[3]);
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}
/// Rotate a vector `v` by quaternion `q` (`[x,y,z,w]`) using `q * v * q^-1`.
pub(super) fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let pv: [f64; 4] = [v[0], v[1], v[2], 0.0];
    let qc: [f64; 4] = [-q[0], -q[1], -q[2], q[3]];
    let tmp = quat_mul(q, pv);
    let res = quat_mul(tmp, qc);
    [res[0], res[1], res[2]]
}
#[cfg(test)]
mod types_expanded_tests {
    use crate::Aabb;
    use crate::BoundingBox3D;
    use crate::Capsule;
    use crate::MassProperties;
    use crate::Plane;
    use crate::Plane3D;
    use crate::Ray;
    use crate::Ray3D;
    use crate::Sphere;
    use crate::Transform;
    use crate::Transform3D;
    use crate::Triangle;
    use crate::Vec3;
    use std::f64::consts::PI;
    #[test]
    fn test_ray_at() {
        let ray = Ray::new(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let p = ray.at(3.0);
        assert!((p.x - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_from_to() {
        let ray = Ray::from_to(Vec3::zeros(), Vec3::new(5.0, 0.0, 0.0));
        assert!((ray.max_t - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_normalizes_direction() {
        let ray = Ray::new(Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0));
        assert!((ray.direction.norm() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_ray_hit() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let hit = aabb.ray_intersect(&ray);
        assert!(hit.is_some(), "ray should hit unit cube");
        let h = hit.unwrap();
        assert!((h.t_min - 4.0).abs() < 1e-10, "t_min={}", h.t_min);
    }
    #[test]
    fn test_aabb_ray_miss() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(-5.0, 5.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(aabb.ray_intersect(&ray).is_none());
    }
    #[test]
    fn test_aabb_ray_inside() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let hit = aabb.ray_intersect(&ray);
        assert!(hit.is_some());
    }
    #[test]
    fn test_aabb_intersection_overlapping() {
        let a = Aabb::new(Vec3::zeros(), Vec3::new(2.0, 2.0, 2.0));
        let b = Aabb::new(Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
        let inter = a.intersection(&b).unwrap();
        assert!((inter.min.x - 1.0).abs() < 1e-10);
        assert!((inter.max.x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_intersection_disjoint() {
        let a = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let b = Aabb::new(Vec3::new(2.0, 2.0, 2.0), Vec3::new(3.0, 3.0, 3.0));
        assert!(a.intersection(&b).is_none());
    }
    #[test]
    fn test_aabb_sq_distance_outside() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(2.0, 0.5, 0.5);
        assert!((aabb.sq_distance_to_point(&p) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_sq_distance_inside_is_zero() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(0.5, 0.5, 0.5);
        assert_eq!(aabb.sq_distance_to_point(&p), 0.0);
    }
    #[test]
    fn test_aabb_closest_point() {
        let aabb = Aabb::new(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        let p = Vec3::new(2.0, 0.5, 0.5);
        let cp = aabb.closest_point(&p);
        assert!((cp.x - 1.0).abs() < 1e-10);
        assert!((cp.y - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_oriented_box() {
        let t = Transform::default();
        let he = Vec3::new(1.0, 2.0, 3.0);
        let aabb = Aabb::oriented_box_aabb(he, &t);
        assert!((aabb.max.x - 1.0).abs() < 1e-10);
        assert!((aabb.max.y - 2.0).abs() < 1e-10);
        assert!((aabb.max.z - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_contains_point() {
        let s = Sphere::new(Vec3::zeros(), 2.0);
        assert!(s.contains_point(&Vec3::new(1.0, 0.0, 0.0)));
        assert!(!s.contains_point(&Vec3::new(3.0, 0.0, 0.0)));
    }
    #[test]
    fn test_sphere_overlaps_sphere() {
        let a = Sphere::new(Vec3::zeros(), 1.0);
        let b = Sphere::new(Vec3::new(1.5, 0.0, 0.0), 1.0);
        let c = Sphere::new(Vec3::new(5.0, 0.0, 0.0), 1.0);
        assert!(a.overlaps_sphere(&b));
        assert!(!a.overlaps_sphere(&c));
    }
    #[test]
    fn test_sphere_volume() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let expected = 4.0 / 3.0 * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_sphere_ray_intersect_hit() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let ray = Ray::new(Vec3::new(-5.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let t = s.ray_intersect(&ray);
        assert!(t.is_some(), "ray should hit sphere");
        assert!((t.unwrap() - 4.0).abs() < 1e-10, "t={}", t.unwrap());
    }
    #[test]
    fn test_sphere_ray_intersect_miss() {
        let s = Sphere::new(Vec3::zeros(), 1.0);
        let ray = Ray::new(Vec3::new(-5.0, 2.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(s.ray_intersect(&ray).is_none());
    }
    #[test]
    fn test_sphere_bounding_sphere() {
        let pts = vec![Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let s = Sphere::bounding_sphere(&pts).unwrap();
        assert!(s.radius >= 1.0 - 1e-9);
    }
    #[test]
    fn test_sphere_bounding_sphere_single_point() {
        let pts = vec![Vec3::new(3.0, 4.0, 0.0)];
        let s = Sphere::bounding_sphere(&pts).unwrap();
        assert_eq!(s.radius, 0.0);
    }
    #[test]
    fn test_capsule_contains_point_on_axis() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        assert!(c.contains_point(&Vec3::new(0.0, 0.0, 0.0)));
    }
    #[test]
    fn test_capsule_volume() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        let v = c.volume();
        assert!(v > 0.0, "volume should be positive: {v}");
    }
    #[test]
    fn test_capsule_aabb_upright() {
        let c = Capsule::upright(Vec3::zeros(), 1.0, 0.5);
        let aabb = c.aabb();
        assert!(aabb.max.y >= 1.0, "y max should be >= 1.0: {}", aabb.max.y);
    }
    #[test]
    fn test_plane_signed_distance() {
        let p = Plane::from_point_normal(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        let above = Vec3::new(0.0, 5.0, 0.0);
        let below = Vec3::new(0.0, -1.0, 0.0);
        assert!(p.signed_distance(&above) > 0.0);
        assert!(p.signed_distance(&below) < 0.0);
    }
    #[test]
    fn test_plane_project_point() {
        let p = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let projected = p.project_point(&Vec3::new(3.0, 5.0, 2.0));
        assert!(
            projected.y.abs() < 1e-10,
            "projected y should be 0: {}",
            projected.y
        );
    }
    #[test]
    fn test_plane_ray_intersect() {
        let plane = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let ray = Ray::new(Vec3::new(0.0, 5.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let t = plane.ray_intersect(&ray);
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-10, "t={}", t.unwrap());
    }
    #[test]
    fn test_plane_ray_parallel_miss() {
        let plane = Plane::from_point_normal(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0));
        let ray = Ray::new(Vec3::new(0.0, 5.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(plane.ray_intersect(&ray).is_none());
    }
    #[test]
    fn test_triangle_area() {
        let t = Triangle::new(
            Vec3::zeros(),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.0, 4.0, 0.0),
        );
        assert!((t.area() - 6.0).abs() < 1e-10, "area={}", t.area());
    }
    #[test]
    fn test_triangle_normal() {
        let t = Triangle::new(
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let n = t.normal();
        assert!((n.z - 1.0).abs() < 1e-10, "n={:?}", n);
    }
    #[test]
    fn test_triangle_ray_hit() {
        let tri = Triangle::new(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let ray = Ray::new(Vec3::new(0.0, 0.3, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let t = tri.ray_intersect(&ray);
        assert!(t.is_some(), "ray should hit triangle");
    }
    #[test]
    fn test_triangle_ray_miss() {
        let tri = Triangle::new(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let ray = Ray::new(Vec3::new(5.0, 5.0, -1.0), Vec3::new(0.0, 0.0, 1.0));
        assert!(tri.ray_intersect(&ray).is_none());
    }
    #[test]
    fn test_triangle_barycentric_centroid() {
        let tri = Triangle::new(
            Vec3::zeros(),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.0, 3.0, 0.0),
        );
        let centroid = tri.centroid();
        let (u, v, w) = tri.barycentric(&centroid);
        assert!((u - 1.0 / 3.0).abs() < 1e-10, "u={}", u);
        assert!((v - 1.0 / 3.0).abs() < 1e-10, "v={}", v);
        assert!((w - 1.0 / 3.0).abs() < 1e-10, "w={}", w);
    }
    #[test]
    fn test_mass_properties_cone() {
        let mp = MassProperties::cone(3.0, 1.0, 2.0);
        assert!((mp.mass - 3.0).abs() < 1e-10);
        assert!(mp.local_inertia[(0, 0)] > 0.0);
    }
    #[test]
    fn test_mass_properties_spherical_shell() {
        let mp = MassProperties::spherical_shell(1.0, 1.0);
        assert!((mp.local_inertia[(0, 0)] - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_mass_properties_capsule() {
        let mp = MassProperties::capsule(5.0, 0.5, 1.0);
        assert!((mp.mass - 5.0).abs() < 1e-10);
        assert!(mp.local_inertia[(0, 0)] > 0.0);
    }
    #[test]
    fn test_mass_properties_shifted_inertia() {
        let mp = MassProperties::sphere(1.0, 1.0);
        let shifted = mp.shifted_inertia(Vec3::new(1.0, 0.0, 0.0));
        assert!(
            (shifted[(0, 0)] - 0.4).abs() < 1e-10,
            "Ixx={}",
            shifted[(0, 0)]
        );
        assert!(
            (shifted[(1, 1)] - 1.4).abs() < 1e-10,
            "Iyy={}",
            shifted[(1, 1)]
        );
    }
    #[test]
    fn test_bounding_box3d_merge_many_empty() {
        let result = BoundingBox3D::merge_many(&[]);
        assert!(result.is_none(), "merge of empty list should be None");
    }
    #[test]
    fn test_bounding_box3d_merge_many_single() {
        let b = BoundingBox3D::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let merged = BoundingBox3D::merge_many(&[b]).unwrap();
        assert_eq!(merged.min, [1.0, 2.0, 3.0]);
        assert_eq!(merged.max, [4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_bounding_box3d_merge_many_two() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = BoundingBox3D::new([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0]);
        let merged = BoundingBox3D::merge_many(&[a, b]).unwrap();
        assert_eq!(merged.min, [-1.0, -1.0, -1.0]);
        assert_eq!(merged.max, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn test_bounding_box3d_intersection_overlap() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = BoundingBox3D::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let inter = a.intersection(&b).unwrap();
        assert_eq!(inter.min, [1.0, 1.0, 1.0]);
        assert_eq!(inter.max, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn test_bounding_box3d_intersection_no_overlap() {
        let a = BoundingBox3D::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = BoundingBox3D::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(
            a.intersection(&b).is_none(),
            "non-overlapping boxes should return None"
        );
    }
    #[test]
    fn test_ray3d_parameter_at_origin() {
        let ray = Ray3D::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = ray.parameter_at(0.0);
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_ray3d_parameter_at_positive_t() {
        let ray = Ray3D::new([1.0, 2.0, 3.0], [0.0, 0.0, 1.0]);
        let p = ray.parameter_at(5.0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_plane3d_project_point_on_plane() {
        let plane = Plane3D::new([0.0, 0.0, 1.0], 0.0);
        let proj = plane.project_point([1.0, 2.0, 0.0]);
        assert!((proj[0] - 1.0).abs() < 1e-10);
        assert!((proj[1] - 2.0).abs() < 1e-10);
        assert!(proj[2].abs() < 1e-10);
    }
    #[test]
    fn test_plane3d_project_point_above_plane() {
        let plane = Plane3D::new([0.0, 0.0, 1.0], 0.0);
        let proj = plane.project_point([1.0, 2.0, 3.0]);
        assert!((proj[0] - 1.0).abs() < 1e-10);
        assert!((proj[1] - 2.0).abs() < 1e-10);
        assert!(proj[2].abs() < 1e-10);
    }
    #[test]
    fn test_transform3d_identity_apply() {
        let t = Transform3D::identity();
        let p = t.apply([3.0, 4.0, 5.0]);
        assert!((p[0] - 3.0).abs() < 1e-10);
        assert!((p[1] - 4.0).abs() < 1e-10);
        assert!((p[2] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform3d_translation_apply() {
        let mut t = Transform3D::identity();
        t.position = [1.0, 2.0, 3.0];
        let p = t.apply([0.0, 0.0, 0.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform3d_inverse_cancels() {
        let mut t = Transform3D::identity();
        t.position = [1.0, 2.0, 3.0];
        let inv = t.inverse();
        let composed = t.compose(&inv);
        assert!(
            composed.position[0].abs() < 1e-9,
            "x={}",
            composed.position[0]
        );
        assert!(
            composed.position[1].abs() < 1e-9,
            "y={}",
            composed.position[1]
        );
        assert!(
            composed.position[2].abs() < 1e-9,
            "z={}",
            composed.position[2]
        );
    }
}
