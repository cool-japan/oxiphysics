//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub(super) fn quaternion_multiply(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    [
        p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1],
        p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0],
        p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3],
        p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2],
    ]
}
pub(super) fn quaternion_normalize(q: [f64; 4]) -> [f64; 4] {
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / len, q[1] / len, q[2] / len, q[3] / len]
}

#[cfg(test)]
mod tests {
    use super::super::types::{
        WasmColliderShape, WasmJointHandle, WasmRaycastResult, WasmSceneSerializer,
    };
    use super::super::types_3::WasmWorld;
    use super::super::types_4::{
        WasmBodyType, WasmJoint, WasmJointType, WasmRigidBodyHandle, WasmSimulationConfig,
    };
    use super::*;
    fn make_world() -> WasmWorld {
        WasmWorld::earth()
    }
    #[test]
    fn test_world_create() {
        let world = make_world();
        assert_eq!(world.body_count(), 0);
        assert_eq!(world.collider_count(), 0);
        assert_eq!(world.joint_count(), 0);
    }
    #[test]
    fn test_add_dynamic_body() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        assert_eq!(world.body_count(), 1);
        let pos = world.get_position(h).unwrap();
        assert!((pos[1] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_add_static_body() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        let state = world.get_body_state(h).unwrap();
        assert_eq!(state.body_type, WasmBodyType::Static);
    }
    #[test]
    fn test_body_falls_under_gravity() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        let y0 = world.get_position(h).unwrap()[1];
        for _ in 0..60 {
            world.step();
        }
        let y1 = world.get_position(h).unwrap()[1];
        assert!(y1 < y0, "body should fall: y0={y0} y1={y1}");
    }
    #[test]
    fn test_apply_impulse() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.apply_impulse(h, [10.0, 0.0, 0.0]);
        let vel = world.get_velocity(h).unwrap();
        assert!((vel[0] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_apply_force() {
        let mut world = make_world();
        let h = world.add_dynamic_body(2.0, 0.0, 0.0, 0.0);
        world.apply_force(h, [20.0, 0.0, 0.0]);
        world.step();
        let vel = world.get_velocity(h).unwrap();
        assert!(vel[0] > 0.0);
    }
    #[test]
    fn test_set_position() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.set_position(h, [5.0, 5.0, 5.0]);
        let pos = world.get_position(h).unwrap();
        assert!((pos[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_remove_body() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        assert!(world.remove_body(h));
        assert_eq!(world.body_count(), 0);
    }
    #[test]
    fn test_sphere_collider() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        let ch = world.add_sphere_collider(h, 0.5);
        assert_eq!(world.collider_count(), 1);
        world.set_friction(ch, 0.3);
        world.set_restitution(ch, 0.5);
    }
    #[test]
    fn test_plane_collider_bounce() {
        let mut world = make_world();
        let ball = world.add_dynamic_body(1.0, 0.0, 2.0, 0.0);
        let sc = world.add_sphere_collider(ball, 0.5);
        world.set_restitution(sc, 0.8);
        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
        for _ in 0..200 {
            world.step();
        }
        let pos = world.get_position(ball).unwrap();
        assert!(pos[1] >= 0.0, "ball should be above or on the ground");
    }
    #[test]
    fn test_add_ball_joint() {
        let mut world = make_world();
        let a = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b = world.add_dynamic_body(1.0, 1.0, 0.0, 0.0);
        let jh = world.add_ball_joint(a, b, [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]);
        assert_eq!(world.joint_count(), 1);
        world.remove_joint(jh);
        assert_eq!(world.joint_count(), 0);
    }
    #[test]
    fn test_add_distance_joint() {
        let mut world = make_world();
        let a = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        let b = world.add_dynamic_body(1.0, 2.0, 0.0, 0.0);
        world.add_distance_joint(a, b, 2.0);
        assert_eq!(world.joint_count(), 1);
    }
    #[test]
    fn test_raycast_sphere_hit() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 1.0);
        let result = world.raycast([0.0, 0.0, -5.0], [0.0, 0.0, 1.0], 100.0);
        assert!(result.hit);
        assert!((result.toi - 4.0).abs() < 0.01);
    }
    #[test]
    fn test_raycast_miss() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.raycast([10.0, 10.0, -5.0], [0.0, 0.0, 1.0], 100.0);
        assert!(!result.hit);
    }
    #[test]
    fn test_raycast_miss_result() {
        let miss = WasmRaycastResult::miss();
        assert!(!miss.hit);
        assert!(miss.toi.is_infinite());
    }
    #[test]
    fn test_aabb_overlap_hit() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.aabb_overlap([-1.0; 3], [1.0; 3]);
        assert_eq!(result.count(), 1);
    }
    #[test]
    fn test_aabb_overlap_miss() {
        let mut world = make_world();
        let h = world.add_static_body(10.0, 0.0, 0.0);
        world.add_sphere_collider(h, 0.5);
        let result = world.aabb_overlap([-1.0; 3], [1.0; 3]);
        assert_eq!(result.count(), 0);
    }
    #[test]
    fn test_drain_events() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.wake_body(h);
        world.step();
        let _events = world.drain_events();
    }
    #[test]
    fn test_contact_events_on_collision() {
        let mut world = make_world();
        let ball = world.add_dynamic_body(1.0, 0.0, 1.0, 0.0);
        world.add_sphere_collider(ball, 0.5);
        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
        for _ in 0..100 {
            world.step();
        }
        let contacts = world.drain_contact_events();
        let _ = contacts;
    }
    #[test]
    fn test_reset() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 0.0, 0.0);
        world.step();
        world.reset();
        assert_eq!(world.body_count(), 0);
        assert!((world.time).abs() < 1e-15);
        assert_eq!(world.step_count, 0);
    }
    #[test]
    fn test_scene_serializer() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        world.step();
        let snapshot = WasmSceneSerializer::serialize(&world);
        assert_eq!(snapshot.bodies.len(), 1);
        let json = WasmSceneSerializer::to_json_string(&snapshot);
        assert!(json.contains("time"));
        assert!(json.contains("body_count"));
    }
    #[test]
    fn test_config_moon_gravity() {
        let cfg = WasmSimulationConfig::moon();
        assert!((cfg.gravity[1] + 1.625).abs() < 1e-10);
    }
    #[test]
    fn test_config_zero_gravity() {
        let cfg = WasmSimulationConfig::zero_gravity();
        assert!((cfg.gravity[1]).abs() < 1e-15);
    }
    #[test]
    fn test_step_with_dt() {
        let mut world = make_world();
        world.add_dynamic_body(1.0, 0.0, 5.0, 0.0);
        world.step_with_dt(0.1);
        assert!((world.config.dt - 1.0 / 60.0).abs() < 1e-10);
    }
    #[test]
    fn test_sensor_collider() {
        let mut world = make_world();
        let h = world.add_static_body(0.0, 0.0, 0.0);
        let ch = world.add_sphere_collider(h, 1.0);
        world.set_sensor(ch, true);
    }
    #[test]
    fn test_gravity_scale() {
        let mut world = make_world();
        let h = world.add_dynamic_body(1.0, 0.0, 10.0, 0.0);
        world.set_gravity_scale(h, 0.0);
        let y0 = world.get_position(h).unwrap()[1];
        for _ in 0..60 {
            world.step();
        }
        let y1 = world.get_position(h).unwrap()[1];
        assert!(
            (y1 - y0).abs() < 1e-3,
            "zero gravity scale: body should not fall"
        );
    }
    #[test]
    fn test_quaternion_multiply_identity() {
        let identity = [0.0, 0.0, 0.0, 1.0];
        let q = [0.1, 0.0, 0.0, 1.0];
        let result = quaternion_multiply(identity, q);
        assert!((result[0] - q[0]).abs() < 1e-10);
    }
    #[test]
    fn test_quaternion_normalize() {
        let q = [0.0, 0.0, 0.0, 2.0];
        let n = quaternion_normalize(q);
        assert!((n[3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_collider_shape_bounding_radius_sphere() {
        let s = WasmColliderShape::Sphere { radius: 2.5 };
        assert!((s.bounding_radius() - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_collider_shape_bounding_radius_box() {
        let b = WasmColliderShape::Box {
            half_extents: [1.0, 1.0, 1.0],
        };
        let r = b.bounding_radius();
        assert!((r - 3.0f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_wasm_joint_types() {
        let ha = WasmRigidBodyHandle(0);
        let hb = WasmRigidBodyHandle(1);
        let jh = WasmJointHandle(0);
        let j = WasmJoint::ball(jh, ha, hb, [0.0; 3], [0.0; 3]);
        assert_eq!(j.joint_type, WasmJointType::Ball);
    }
    #[test]
    fn test_multiple_bodies_simulation() {
        let mut world = make_world();
        let mut handles = Vec::new();
        for i in 0..5 {
            let h = world.add_dynamic_body(1.0, i as f64, 5.0, 0.0);
            world.add_sphere_collider(h, 0.3);
            handles.push(h);
        }
        let ground = world.add_static_body(0.0, 0.0, 0.0);
        world.add_plane_collider(ground, 0.0, 1.0, 0.0, 0.0);
        for _ in 0..60 {
            world.step();
        }
        for h in &handles {
            let pos = world.get_position(*h).unwrap();
            assert!(pos[1] >= -0.1, "body should not fall through ground");
        }
    }
}
