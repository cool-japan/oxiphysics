//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Ray;

#[cfg(test)]
mod tests {

    use crate::Collider;
    use crate::RigidBody;

    use crate::world::PhysicsWorld;

    use crate::world::RigidWorld;
    use crate::world::WorldConfig;
    use oxiphysics_core::BodyHandle;
    #[test]
    fn test_add_remove_body_count() {
        let mut world = RigidWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        assert_eq!(world.body_count(), 0);
        let id0 = world.add_body([0.0; 3], [0.0; 3], 1.0);
        let id1 = world.add_body([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert_eq!(world.body_count(), 2);
        assert!(world.remove_body(id0));
        assert_eq!(world.body_count(), 1);
        assert!(!world.remove_body(id0));
        assert_eq!(world.body_count(), 1);
        assert!(world.remove_body(id1));
        assert_eq!(world.body_count(), 0);
    }
    #[test]
    fn test_step_free_fall_accelerates_downward() {
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 60.0;
        let mut world = RigidWorld::new(gravity, dt);
        let id = world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.step();
        let body = world.get_body(id).unwrap();
        let expected_vy = gravity[1] * dt;
        assert!(
            (body.velocity[1] - expected_vy).abs() < 1e-12,
            "vy = {}, expected {}",
            body.velocity[1],
            expected_vy
        );
        let expected_y = expected_vy * dt;
        assert!(
            (body.position[1] - expected_y).abs() < 1e-12,
            "y = {}, expected {}",
            body.position[1],
            expected_y
        );
    }
    #[test]
    fn test_static_body_stays_put() {
        let mut world = RigidWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        let id = world.add_body([0.0; 3], [0.0; 3], 0.0);
        world.step();
        world.step();
        world.step();
        let body = world.get_body(id).unwrap();
        assert_eq!(body.velocity, [0.0, 0.0, 0.0]);
        assert_eq!(body.position, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_kinetic_energy_two_unit_mass_bodies_speed_one() {
        let mut world = RigidWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        world.add_body([0.0; 3], [1.0, 0.0, 0.0], 1.0);
        world.add_body([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        let ke = world.kinetic_energy();
        assert!(
            (ke - 1.0).abs() < 1e-12,
            "kinetic energy = {}, expected 1.0",
            ke
        );
    }
    #[test]
    fn test_apply_impulse_changes_velocity() {
        let mut world = RigidWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let id = world.add_body([0.0; 3], [0.0; 3], 0.5);
        world.apply_impulse_to_body(id, [4.0, 0.0, 0.0]);
        let body = world.get_body(id).unwrap();
        assert!(
            (body.velocity[0] - 2.0).abs() < 1e-12,
            "vx = {}, expected 2.0",
            body.velocity[0]
        );
        assert!((body.velocity[1]).abs() < 1e-12);
        assert!((body.velocity[2]).abs() < 1e-12);
    }
    fn make_unit_body() -> RigidBody {
        let mut b = RigidBody::new(1.0);
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        b
    }
    #[test]
    fn test_physics_world_body_count() {
        let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        assert_eq!(world.body_count(), 0);
        let h0 = world.add_rigid_body(make_unit_body());
        let h1 = world.add_rigid_body(make_unit_body());
        assert_eq!(world.body_count(), 2);
        world.remove_rigid_body(h0);
        assert_eq!(world.body_count(), 1);
        world.remove_rigid_body(h1);
        assert_eq!(world.body_count(), 0);
    }
    #[test]
    fn test_physics_world_gravity_integration() {
        let dt = 1.0 / 60.0;
        let gravity = [0.0, -9.81, 0.0];
        let mut world = PhysicsWorld::new(gravity, dt);
        let mut body = make_unit_body();
        body.gravity_scale = 1.0;
        let h = world.add_rigid_body(body);
        world.step();
        let b = world.get_body(h).unwrap();
        let expected_vy = gravity[1] * dt;
        assert!(
            (b.velocity.y - expected_vy).abs() < 1e-10,
            "vy={} expected={}",
            b.velocity.y,
            expected_vy
        );
    }
    #[test]
    fn test_physics_world_impulse() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let h = world.add_rigid_body(make_unit_body());
        world.apply_impulse(h, [5.0, 0.0, 0.0]);
        let b = world.get_body(h).unwrap();
        assert!((b.velocity.x - 5.0).abs() < 1e-10, "vx={}", b.velocity.x);
    }
    #[test]
    fn test_physics_world_apply_force() {
        let dt = 0.1;
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        let h = world.add_rigid_body(make_unit_body());
        world.apply_force(h, [10.0, 0.0, 0.0]);
        world.step();
        let b = world.get_body(h).unwrap();
        assert!(
            (b.velocity.x - 1.0).abs() < 1e-9,
            "vx={} expected=1.0",
            b.velocity.x
        );
    }
    #[test]
    fn test_physics_world_get_body_invalid_handle() {
        let world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let invalid = BodyHandle::new(99, 99);
        assert!(world.get_body(invalid).is_none());
    }
    #[test]
    fn test_physics_world_active_body_count() {
        let mut world = PhysicsWorld::default_earth();
        world.add_rigid_body(make_unit_body());
        world.add_rigid_body(make_unit_body());
        assert_eq!(world.active_body_count(), 2);
    }
    #[test]
    fn test_physics_world_sleep_reduces_active_count() {
        let mut world = PhysicsWorld::default_earth();
        world.time_before_sleep = 0.1;
        world.linear_sleep_threshold = 100.0;
        world.angular_sleep_threshold = 100.0;
        let mut body = make_unit_body();
        body.gravity_scale = 0.0;
        world.add_rigid_body(body);
        let steps = ((world.time_before_sleep / world.dt) as usize) + 5;
        for _ in 0..steps {
            world.step();
        }
        assert_eq!(world.sleeping_body_count(), 1);
        assert_eq!(world.active_body_count(), 0);
    }
    #[test]
    fn test_physics_world_set_gravity() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        world.set_gravity([0.0, -9.81, 0.0]);
        assert!((world.gravity[1] + 9.81).abs() < 1e-12);
    }
    #[test]
    fn test_physics_world_step_advances_time() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let dt = world.dt;
        world.step();
        world.step();
        assert!((world.time - 2.0 * dt).abs() < 1e-12);
    }
    #[test]
    fn test_physics_world_static_body_unaffected() {
        let dt = 1.0 / 60.0;
        let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], dt);
        let body = RigidBody::new_static();
        let h = world.add_rigid_body(body);
        for _ in 0..10 {
            world.step();
        }
        let b = world.get_body(h).unwrap();
        assert!(b.transform.position.norm() < 1e-10, "static body moved");
    }
    #[test]
    fn test_physics_world_add_collider() {
        use oxiphysics_geometry::Sphere;
        use std::sync::Arc;
        let mut world = PhysicsWorld::default_earth();
        let h = world.add_rigid_body(make_unit_body());
        let shape: Arc<dyn oxiphysics_geometry::Shape> = Arc::new(Sphere::new(0.5));
        let collider = Collider::new(shape);
        let ch = world.add_collider(collider, Some(h));
        assert!(world.colliders.get(ch).is_some());
    }
    #[test]
    fn test_world_config_apply_and_get_roundtrip() {
        let mut world = PhysicsWorld::default_earth();
        let cfg = WorldConfig {
            gravity: [0.0, -1.62, 0.0],
            dt: 1.0 / 120.0,
            ..Default::default()
        };
        world.apply_config(&cfg);
        let back = world.get_config();
        assert!((back.gravity[1] + 1.62).abs() < 1e-12);
        assert!((back.dt - 1.0 / 120.0).abs() < 1e-12);
    }
    #[test]
    fn test_world_config_default_values() {
        let cfg = WorldConfig::default();
        assert!((cfg.gravity[1] + 9.81).abs() < 1e-10);
        assert!((cfg.dt - 1.0 / 60.0).abs() < 1e-10);
    }
    #[test]
    fn test_world_statistics_body_counts() {
        let mut world = PhysicsWorld::default_earth();
        world.add_rigid_body(make_unit_body());
        world.add_rigid_body(make_unit_body());
        world.add_rigid_body(RigidBody::new_static());
        let stats = world.statistics();
        assert_eq!(stats.total_bodies, 3);
        assert_eq!(stats.static_bodies, 1);
        assert_eq!(stats.sleeping_bodies, 0);
    }
    #[test]
    fn test_world_statistics_sim_time_advances() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        world.step();
        world.step();
        let stats = world.statistics();
        assert!((stats.sim_time - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_world_statistics_kinetic_energy_nonzero_after_gravity() {
        let mut world = PhysicsWorld::default_earth();
        let mut b = make_unit_body();
        b.gravity_scale = 1.0;
        world.add_rigid_body(b);
        world.step();
        let stats = world.statistics();
        assert!(
            stats.kinetic_energy > 0.0,
            "KE should be nonzero after gravity step"
        );
    }
    #[test]
    fn test_bodies_within_radius_finds_nearby() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        let h0 = world.add_rigid_body(make_unit_body());
        let mut b1 = make_unit_body();
        b1.transform.position = oxiphysics_core::math::Vec3::new(10.0, 0.0, 0.0);
        let _h1 = world.add_rigid_body(b1);
        let found = world.bodies_within_radius([0.0; 3], 1.0);
        assert!(found.contains(&h0));
        assert_eq!(found.len(), 1);
    }
    #[test]
    fn test_bodies_above_y_threshold() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        let _h0 = world.add_rigid_body(make_unit_body());
        let mut b1 = make_unit_body();
        b1.transform.position = oxiphysics_core::math::Vec3::new(0.0, 5.0, 0.0);
        let h1 = world.add_rigid_body(b1);
        let above = world.bodies_above_y(3.0);
        assert!(above.contains(&h1));
        assert_eq!(above.len(), 1);
    }
    #[test]
    fn test_bodies_exceeding_ke_after_impulse() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        let h = world.add_rigid_body(make_unit_body());
        world.apply_impulse(h, [100.0, 0.0, 0.0]);
        let high_ke = world.bodies_exceeding_ke(1.0);
        assert!(high_ke.contains(&h));
    }
    #[test]
    fn test_apply_force_batch() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        let h0 = world.add_rigid_body(make_unit_body());
        let h1 = world.add_rigid_body(make_unit_body());
        let handles = vec![h0, h1];
        world.apply_force_batch(&handles, [10.0, 0.0, 0.0]);
        world.step();
        let b0 = world.get_body(h0).unwrap();
        let b1 = world.get_body(h1).unwrap();
        assert!(b0.velocity.x > 0.0, "b0 vx={}", b0.velocity.x);
        assert!(b1.velocity.x > 0.0, "b1 vx={}", b1.velocity.x);
    }
    #[test]
    fn test_apply_impulse_batch() {
        let mut world = PhysicsWorld::new([0.0; 3], 0.1);
        let h0 = world.add_rigid_body(make_unit_body());
        let h1 = world.add_rigid_body(make_unit_body());
        let handles = vec![h0, h1];
        world.apply_impulse_batch(&handles, [5.0, 0.0, 0.0]);
        assert!((world.get_body(h0).unwrap().velocity.x - 5.0).abs() < 1e-10);
        assert!((world.get_body(h1).unwrap().velocity.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_wake_all_wakes_sleeping_bodies() {
        let mut world = PhysicsWorld::default_earth();
        world.time_before_sleep = 0.01;
        world.linear_sleep_threshold = 1000.0;
        world.angular_sleep_threshold = 1000.0;
        let mut b = make_unit_body();
        b.gravity_scale = 0.0;
        world.add_rigid_body(b);
        for _ in 0..100 {
            world.step();
        }
        assert_eq!(world.sleeping_body_count(), 1);
        world.wake_all();
        assert_eq!(world.sleeping_body_count(), 0);
    }
    #[test]
    fn test_set_gravity_scale_all() {
        let mut world = PhysicsWorld::default_earth();
        world.add_rigid_body(make_unit_body());
        world.add_rigid_body(make_unit_body());
        world.set_gravity_scale_all(crate::body::BodyType::Dynamic, 0.0);
        for (_, b) in world.bodies.iter() {
            assert!((b.gravity_scale).abs() < 1e-12);
        }
    }
    #[test]
    fn test_rigid_world_snapshot_roundtrip() {
        let mut world = RigidWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        world.add_body([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 1.0);
        world.add_body([4.0, 5.0, 6.0], [0.4, 0.5, 0.6], 0.5);
        world.step();
        let snap = world.to_snapshot();
        let restored = RigidWorld::from_snapshot(snap);
        assert_eq!(restored.body_count(), world.body_count());
        assert!((restored.time - world.time).abs() < 1e-12);
        assert_eq!(restored.gravity, world.gravity);
    }
    #[test]
    fn test_rigid_world_activate_deactivate() {
        let mut world = RigidWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        let id = world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.deactivate_body(id);
        let pos_before = world.get_body(id).unwrap().position;
        world.step();
        let pos_after = world.get_body(id).unwrap().position;
        assert_eq!(pos_before, pos_after, "deactivated body should not move");
        world.activate_body(id);
        world.step();
        let pos_after2 = world.get_body(id).unwrap().position;
        assert!(pos_after2[1] < pos_before[1], "activated body should fall");
    }
    #[test]
    fn test_rigid_world_active_body_count() {
        let mut world = RigidWorld::new([0.0; 3], 0.1);
        let id0 = world.add_body([0.0; 3], [0.0; 3], 1.0);
        let _id1 = world.add_body([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert_eq!(world.active_body_count(), 2);
        world.deactivate_body(id0);
        assert_eq!(world.active_body_count(), 1);
    }
    #[test]
    fn test_rigid_world_global_impulse() {
        let mut world = RigidWorld::new([0.0; 3], 0.1);
        let id = world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.apply_global_impulse([3.0, 0.0, 0.0]);
        let b = world.get_body(id).unwrap();
        assert!((b.velocity[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_rigid_world_most_energetic_body() {
        let mut world = RigidWorld::new([0.0; 3], 0.1);
        let id_slow = world.add_body([0.0; 3], [1.0, 0.0, 0.0], 1.0);
        let id_fast = world.add_body([5.0, 0.0, 0.0], [10.0, 0.0, 0.0], 1.0);
        let most = world.most_energetic_body().unwrap();
        assert_eq!(
            most.id, id_fast,
            "fast body should be most energetic (slow={id_slow})"
        );
    }
    #[test]
    fn test_rigid_world_step_n() {
        let mut world = RigidWorld::new([0.0, -9.81, 0.0], 0.1);
        world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.step_n(5);
        assert!((world.time - 0.5).abs() < 1e-10, "time={}", world.time);
    }
}
/// Helper: ray–sphere intersection.
///
/// Returns the smallest positive `t` along the ray that hits a sphere of
/// `radius` centred at `centre`, or `None` if there is no intersection.
pub(super) fn ray_sphere_intersect(ray: &Ray, centre: [f64; 3], radius: f64) -> Option<f64> {
    let oc = [
        ray.origin[0] - centre[0],
        ray.origin[1] - centre[1],
        ray.origin[2] - centre[2],
    ];
    let d = ray.direction;
    let a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let b = 2.0 * (oc[0] * d[0] + oc[1] * d[1] + oc[2] * d[2]);
    let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 || a < 1e-30 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t0 = (-b - sqrt_disc) / (2.0 * a);
    let t1 = (-b + sqrt_disc) / (2.0 * a);
    if t0 > 1e-6 {
        Some(t0)
    } else if t1 > 1e-6 {
        Some(t1)
    } else {
        None
    }
}
#[cfg(test)]
mod tests_world_new {

    use crate::RigidBody;
    use crate::world::AabbQuery;
    use crate::world::BodyConfig;
    use crate::world::GravityRegion;
    use crate::world::PhysicsWorld;
    use crate::world::Ray;
    use crate::world::RigidWorld;

    fn make_world() -> RigidWorld {
        RigidWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0)
    }
    #[test]
    fn test_raycast_hits_body_directly_ahead() {
        let mut world = make_world();
        world.add_body([5.0, 0.0, 0.0], [0.0; 3], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let result = world.raycast_nearest(&ray);
        assert!(result.is_some(), "should hit body at x=5");
        let r = result.unwrap();
        assert!(r.t > 0.0 && r.t < 10.0, "t = {}", r.t);
    }
    #[test]
    fn test_raycast_misses_body_to_side() {
        let mut world = make_world();
        world.add_body([0.0, 100.0, 0.0], [0.0; 3], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let result = world.raycast_nearest(&ray);
        assert!(result.is_none(), "should miss");
    }
    #[test]
    fn test_raycast_all_sorted_by_t() {
        let mut world = make_world();
        world.add_body([10.0, 0.0, 0.0], [0.0; 3], 1.0);
        world.add_body([3.0, 0.0, 0.0], [0.0; 3], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hits = world.raycast_all(&ray);
        if hits.len() >= 2 {
            assert!(
                hits[0].t <= hits[1].t,
                "hits should be sorted ascending by t"
            );
        }
    }
    #[test]
    fn test_raycast_empty_world() {
        let world = make_world();
        let ray = Ray::new([0.0; 3], [0.0, 1.0, 0.0]);
        assert!(world.raycast_nearest(&ray).is_none());
    }
    #[test]
    fn test_sphere_cast_wider_than_ray() {
        let mut world = make_world();
        world.add_body([5.0, 0.5, 0.0], [0.0; 3], 1.0);
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let hits = world.sphere_cast(&ray, 0.8);
        assert!(!hits.is_empty(), "sphere cast should hit the nearby body");
    }
    #[test]
    fn test_sphere_cast_empty_world() {
        let world = make_world();
        let ray = Ray::new([0.0; 3], [1.0, 0.0, 0.0]);
        assert!(world.sphere_cast(&ray, 1.0).is_empty());
    }
    #[test]
    fn test_aabb_overlap_query_finds_inside_bodies() {
        let mut world = make_world();
        let id_in = world.add_body([1.0, 1.0, 1.0], [0.0; 3], 1.0);
        let _id_out = world.add_body([100.0, 0.0, 0.0], [0.0; 3], 1.0);
        let query = AabbQuery::new([0.0; 3], [5.0; 3]);
        let found = world.aabb_overlap_query(&query);
        assert!(found.contains(&id_in), "should find body at [1,1,1]");
        assert_eq!(found.len(), 1);
    }
    #[test]
    fn test_aabb_query_contains_point() {
        let q = AabbQuery::new([0.0; 3], [2.0; 3]);
        assert!(q.contains_point([1.0, 1.0, 1.0]));
        assert!(!q.contains_point([3.0, 1.0, 1.0]));
    }
    #[test]
    fn test_aabb_overlaps() {
        let q1 = AabbQuery::new([0.0; 3], [2.0; 3]);
        let q2 = AabbQuery::new([1.0; 3], [3.0; 3]);
        let q3 = AabbQuery::new([5.0; 3], [6.0; 3]);
        assert!(q1.overlaps(&q2));
        assert!(!q1.overlaps(&q3));
    }
    #[test]
    fn test_aabb_bounding_sphere_query() {
        let mut world = make_world();
        let id = world.add_body([0.5, 0.5, 0.5], [0.0; 3], 1.0);
        let query = AabbQuery::new([-1.0; 3], [1.5; 3]);
        let found = world.aabb_bounding_sphere_query(&query);
        assert!(found.contains(&id));
    }
    #[test]
    fn test_add_bodies_from_configs_count() {
        let mut world = make_world();
        let configs = vec![
            BodyConfig::dynamic_at([0.0; 3]),
            BodyConfig::dynamic_at([1.0, 0.0, 0.0]),
            BodyConfig::static_at([2.0, 0.0, 0.0]),
        ];
        let ids = world.add_bodies_from_configs(&configs);
        assert_eq!(ids.len(), 3);
        assert_eq!(world.body_count(), 3);
    }
    #[test]
    fn test_add_bodies_from_configs_ids_unique() {
        let mut world = make_world();
        let configs = vec![
            BodyConfig::dynamic_at([0.0; 3]),
            BodyConfig::dynamic_at([1.0, 0.0, 0.0]),
        ];
        let ids = world.add_bodies_from_configs(&configs);
        assert_ne!(ids[0], ids[1], "ids should be unique");
    }
    #[test]
    fn test_body_config_static_has_zero_inv_mass() {
        let cfg = BodyConfig::static_at([0.0; 3]);
        assert_eq!(cfg.inv_mass, 0.0);
    }
    #[test]
    fn test_state_vector_roundtrip() {
        let mut world = make_world();
        world.add_body([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 1.0);
        world.add_body([4.0, 5.0, 6.0], [0.0; 3], 0.5);
        world.step();
        let sv = world.to_state_vector();
        let restored = RigidWorld::from_state_vector(&sv, world.dt);
        assert_eq!(restored.body_count(), world.body_count());
        assert!((restored.time - world.time).abs() < 1e-12);
        assert_eq!(restored.gravity, world.gravity);
    }
    #[test]
    fn test_state_vector_length() {
        let mut world = make_world();
        world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.add_body([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        let sv = world.to_state_vector();
        assert_eq!(
            sv.len(),
            2 * 9 + 4,
            "expected {} values, got {}",
            2 * 9 + 4,
            sv.len()
        );
    }
    #[test]
    fn test_state_vector_positions_preserved() {
        let mut world = make_world();
        world.add_body([7.0, 8.0, 9.0], [0.0; 3], 1.0);
        let sv = world.to_state_vector();
        let restored = RigidWorld::from_state_vector(&sv, world.dt);
        let b = restored.bodies.first().unwrap();
        assert!((b.position[0] - 7.0).abs() < 1e-12);
        assert!((b.position[1] - 8.0).abs() < 1e-12);
        assert!((b.position[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_gravity_region_inside_applies_gravity() {
        let region = GravityRegion::new([0.0; 3], 10.0, [0.0, 5.0, 0.0]);
        let g = region.gravity_at([0.0, 0.0, 0.0]);
        assert!(g.is_some());
        assert!((g.unwrap()[1] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_gravity_region_outside_returns_none() {
        let region = GravityRegion::new([0.0; 3], 1.0, [0.0, -9.81, 0.0]);
        let g = region.gravity_at([100.0, 0.0, 0.0]);
        assert!(g.is_none());
    }
    #[test]
    fn test_gravity_region_falloff_at_boundary() {
        let region = GravityRegion::new([0.0; 3], 10.0, [0.0, -10.0, 0.0]).with_falloff();
        let g_centre = region.gravity_at([0.0, 0.0, 0.0]).unwrap();
        assert!((g_centre[1] - (-10.0)).abs() < 1e-10);
        let g_half = region.gravity_at([5.0, 0.0, 0.0]).unwrap();
        assert!((g_half[1] - (-5.0)).abs() < 1e-10);
        assert!(region.gravity_at([10.1, 0.0, 0.0]).is_none());
    }
    #[test]
    fn test_apply_gravity_regions_affects_velocity() {
        let mut world = RigidWorld::new([0.0; 3], 0.1);
        world.add_body([0.0; 3], [0.0; 3], 1.0);
        let regions = vec![GravityRegion::new([0.0; 3], 10.0, [0.0, -9.81, 0.0])];
        world.apply_gravity_regions(&regions);
        let b = &world.bodies[0];
        assert!(
            b.velocity[1] < 0.0,
            "body should have gained downward velocity"
        );
    }
    #[test]
    fn test_rigid_world_stats_counts() {
        let mut world = make_world();
        world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.add_body([1.0, 0.0, 0.0], [0.0; 3], 0.0);
        let id = world.add_body([2.0, 0.0, 0.0], [0.0; 3], 1.0);
        world.deactivate_body(id);
        let stats = world.stats();
        assert_eq!(stats.total_bodies, 3);
        assert_eq!(stats.static_bodies, 1);
        assert_eq!(stats.active_bodies, 1);
        assert_eq!(stats.inactive_bodies, 1);
    }
    #[test]
    fn test_rigid_world_stats_sim_time() {
        let mut world = RigidWorld::new([0.0; 3], 0.1);
        world.add_body([0.0; 3], [0.0; 3], 1.0);
        world.step_n(3);
        let stats = world.stats();
        assert!(
            (stats.sim_time - 0.3).abs() < 1e-10,
            "sim_time = {}",
            stats.sim_time
        );
    }
    fn make_unit_body_pw() -> RigidBody {
        let mut b = RigidBody::new(1.0);
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        b
    }
    #[test]
    fn test_step_substep_advances_time() {
        let dt = 1.0 / 60.0;
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        world.add_rigid_body(make_unit_body_pw());
        world.step_substep(dt * 3.0);
        assert!(
            (world.time - dt * 3.0).abs() < 1e-9,
            "time after 3 sub-steps: {} expected {}",
            world.time,
            dt * 3.0
        );
    }
    #[test]
    fn test_step_substep_single_step_matches_step() {
        let dt = 1.0 / 60.0;
        let gravity = [0.0, -9.81, 0.0];
        let mut world_a = PhysicsWorld::new(gravity, dt);
        let ha = world_a.add_rigid_body(make_unit_body_pw());
        world_a.step();
        let mut world_b = PhysicsWorld::new(gravity, dt);
        let hb = world_b.add_rigid_body(make_unit_body_pw());
        world_b.step_substep(dt);
        let vy_a = world_a.get_body(ha).unwrap().velocity.y;
        let vy_b = world_b.get_body(hb).unwrap().velocity.y;
        assert!(
            (vy_a - vy_b).abs() < 1e-10,
            "single substep should match single step: {vy_a} vs {vy_b}"
        );
    }
    #[test]
    fn test_step_substep_zero_dt_no_advance() {
        let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
        world.add_rigid_body(make_unit_body_pw());
        world.step_substep(0.0);
        assert_eq!(world.time, 0.0, "zero total_dt should not advance time");
    }
    #[test]
    fn test_energy_error_zero_when_ke_equals_reference() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let h = world.add_rigid_body(make_unit_body_pw());
        world.apply_impulse(h, [2.0, 0.0, 0.0]);
        let ke = world.kinetic_energy();
        let err = world.compute_energy_error(ke);
        assert!(
            err.abs() < 1e-10,
            "error should be zero when reference = current KE: {err}"
        );
    }
    #[test]
    fn test_energy_error_nonzero_when_different() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let h = world.add_rigid_body(make_unit_body_pw());
        world.apply_impulse(h, [1.0, 0.0, 0.0]);
        let err = world.compute_energy_error(100.0);
        assert!(
            err > 0.0,
            "error should be positive when KE != reference: {err}"
        );
    }
    #[test]
    fn test_energy_error_finite_when_reference_zero() {
        let world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
        let err = world.compute_energy_error(0.0);
        assert!(err.is_finite(), "energy error should be finite: {err}");
    }
    #[test]
    fn test_wind_force_accelerates_body_in_wind_direction() {
        let dt = 0.1;
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        let h = world.add_rigid_body(make_unit_body_pw());
        world.apply_wind_force([10.0, 0.0, 0.0], 1.2);
        world.step();
        let vx = world.get_body(h).unwrap().velocity.x;
        assert!(
            vx > 0.0,
            "body should accelerate in wind direction: vx={vx}"
        );
    }
    #[test]
    fn test_wind_force_no_effect_on_static_body() {
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 0.1);
        let h = world.add_rigid_body(RigidBody::new_static());
        world.apply_wind_force([100.0, 0.0, 0.0], 1.2);
        world.step();
        let vx = world.get_body(h).unwrap().velocity.x;
        assert!(
            vx.abs() < 1e-12,
            "static body should not be affected by wind: vx={vx}"
        );
    }
    #[test]
    fn test_wind_force_zero_wind_no_acceleration() {
        let dt = 0.1;
        let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        let h = world.add_rigid_body(make_unit_body_pw());
        world.apply_wind_force([0.0, 0.0, 0.0], 1.2);
        world.step();
        let v = world.get_body(h).unwrap().velocity;
        assert!(
            v.x.abs() < 1e-12 && v.y.abs() < 1e-12 && v.z.abs() < 1e-12,
            "zero wind should produce no force"
        );
    }
    #[test]
    fn test_wind_force_stronger_wind_more_acceleration() {
        let dt = 0.1;
        let mut world_weak = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        let hw = world_weak.add_rigid_body(make_unit_body_pw());
        world_weak.apply_wind_force([5.0, 0.0, 0.0], 1.2);
        world_weak.step();
        let mut world_strong = PhysicsWorld::new([0.0, 0.0, 0.0], dt);
        let hs = world_strong.add_rigid_body(make_unit_body_pw());
        world_strong.apply_wind_force([20.0, 0.0, 0.0], 1.2);
        world_strong.step();
        let vx_weak = world_weak.get_body(hw).unwrap().velocity.x;
        let vx_strong = world_strong.get_body(hs).unwrap().velocity.x;
        assert!(
            vx_strong > vx_weak,
            "stronger wind should produce more acceleration: weak={vx_weak}, strong={vx_strong}"
        );
    }
}
