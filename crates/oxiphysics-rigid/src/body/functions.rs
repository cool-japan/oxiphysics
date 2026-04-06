//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::BodyState;
    use crate::BodyType;

    use crate::CenterOfMassTracker;

    use crate::ConstraintList;
    use crate::ConstraintType;
    use crate::ContactHistory;
    use crate::ForceAccumulator;

    use crate::PersistentForce;
    use crate::RigidBody;

    use crate::StateInterpolator;

    use nalgebra::Unit;
    use oxiphysics_core::Quat;
    use oxiphysics_core::Real;
    use oxiphysics_core::Steppable;
    use oxiphysics_core::TimeStep;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_geometry::Shape;
    #[test]
    fn test_rigid_body_step() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        body.step(&TimeStep::new(1.0));
        assert!((body.transform.position.x - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_static_body() {
        let body = RigidBody::new_static();
        assert_eq!(body.body_type, BodyType::Static);
        assert!((body.inverse_mass - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_apply_impulse() {
        let mut body = RigidBody::new(2.0);
        body.apply_impulse(Vec3::new(4.0, 0.0, 0.0));
        assert!((body.velocity.x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_gravity_integration() {
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        body.integrate_forces(1.0, &gravity);
        assert!((body.velocity.y + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_force_accumulator() {
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        body.apply_force(Vec3::new(10.0, 0.0, 0.0));
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        body.integrate_forces(0.1, &gravity);
        assert!((body.velocity.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sleep() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::zeros();
        body.angular_velocity = Vec3::zeros();
        for _ in 0..100 {
            body.check_sleep(0.01, 0.01, 0.01, 0.5);
        }
        assert_eq!(body.state, BodyState::Sleeping);
    }
    #[test]
    fn test_body_set_mass_from_box() {
        use oxiphysics_geometry::BoxShape;
        let mut body = RigidBody::new(1.0);
        let density = 500.0;
        let box_shape = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        body.set_mass_from_shape(&box_shape, density);
        let expected_mass = density * box_shape.volume();
        assert!(
            (body.mass - expected_mass).abs() < 1e-6,
            "mass={} expected={}",
            body.mass,
            expected_mass
        );
        assert!((body.inverse_mass - 1.0 / expected_mass).abs() < 1e-10);
    }
    #[test]
    fn test_body_inertia_nonzero() {
        use oxiphysics_geometry::Sphere;
        let mut body = RigidBody::new(10.0);
        let sphere = Sphere::new(0.5);
        body.set_inertia_from_shape(&sphere);
        assert!(body.world_inverse_inertia[(0, 0)] > 0.0);
        assert!(body.world_inverse_inertia[(1, 1)] > 0.0);
        assert!(body.world_inverse_inertia[(2, 2)] > 0.0);
    }
    #[test]
    fn test_body_sleep_timer_advances() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::zeros();
        body.angular_velocity = Vec3::zeros();
        let dt = 0.016;
        let lin_threshold = 0.1;
        let ang_threshold = 0.1;
        let time_before_sleep = 1.0;
        body.check_sleep(dt, lin_threshold, ang_threshold, time_before_sleep);
        assert!(body.sleep_timer > 0.0, "sleep_timer should have advanced");
        assert!((body.sleep_timer - dt).abs() < 1e-10);
        assert_eq!(body.state, BodyState::Active);
    }
    #[test]
    fn test_body_goes_to_sleep() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::zeros();
        body.angular_velocity = Vec3::zeros();
        let dt = 0.016_f64 as Real;
        let time_before_sleep = 0.5_f64 as Real;
        let steps = (time_before_sleep / dt).ceil() as u32 + 1;
        for _ in 0..steps {
            body.check_sleep(dt, 0.01, 0.01, time_before_sleep);
        }
        assert_eq!(body.state, BodyState::Sleeping);
    }
    #[test]
    fn test_body_wakes_up() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::zeros();
        body.angular_velocity = Vec3::zeros();
        let dt = 0.016_f64 as Real;
        let time_before_sleep = 0.5_f64 as Real;
        let steps = (time_before_sleep / dt).ceil() as u32 + 1;
        for _ in 0..steps {
            body.check_sleep(dt, 0.01, 0.01, time_before_sleep);
        }
        assert_eq!(body.state, BodyState::Sleeping);
        body.velocity = Vec3::new(3.0, 0.0, 0.0);
        body.wake_up();
        assert_eq!(body.state, BodyState::Active);
        assert!((body.velocity.x - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_kinematic_body_not_integrated() {
        let mut body = RigidBody::new_kinematic();
        body.velocity = Vec3::new(5.0, 0.0, 0.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let vel_before = body.velocity;
        body.integrate_forces(0.016, &gravity);
        assert!((body.velocity.x - vel_before.x).abs() < 1e-10);
        assert!((body.velocity.y - vel_before.y).abs() < 1e-10);
        assert!((body.velocity.z - vel_before.z).abs() < 1e-10);
    }
    #[test]
    fn test_static_body_immovable() {
        let mut body = RigidBody::new_static();
        body.velocity = Vec3::new(10.0, 0.0, 0.0);
        let pos_before = body.transform.position;
        body.integrate_velocity(0.016);
        assert!((body.transform.position.x - pos_before.x).abs() < 1e-10);
        assert!((body.transform.position.y - pos_before.y).abs() < 1e-10);
        assert!((body.transform.position.z - pos_before.z).abs() < 1e-10);
    }
    #[test]
    fn test_damping_reduces_velocity() {
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.1;
        body.gravity_scale = 0.0;
        body.velocity = Vec3::new(10.0, 0.0, 0.0);
        let dt = 0.016_f64 as Real;
        let initial_vx = body.velocity.x;
        let gravity = Vec3::zeros();
        body.integrate_forces(dt, &gravity);
        let expected_vx = initial_vx * (1.0 - 0.1_f64 as Real);
        assert!(
            (body.velocity.x - expected_vx).abs() < 1e-6,
            "velocity.x={} expected={}",
            body.velocity.x,
            expected_vx
        );
    }
    #[test]
    fn test_kinetic_energy() {
        let mut body = RigidBody::new(2.0);
        body.velocity = Vec3::new(3.0, 0.0, 0.0);
        let ke = body.kinetic_energy();
        assert!((ke - 9.0).abs() < 1e-10, "KE should be 9.0, got {ke}");
    }
    #[test]
    fn test_linear_momentum() {
        let mut body = RigidBody::new(5.0);
        body.velocity = Vec3::new(2.0, 0.0, 0.0);
        let p = body.linear_momentum();
        assert!((p.x - 10.0).abs() < 1e-10);
        assert!(p.y.abs() < 1e-10);
    }
    #[test]
    fn test_angular_momentum() {
        let mut body = RigidBody::new(1.0);
        body.angular_velocity = Vec3::new(0.0, 0.0, 5.0);
        let l = body.angular_momentum();
        assert!((l.z - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_velocity_at_point() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
        body.angular_velocity = Vec3::new(0.0, 0.0, 1.0);
        let point = body.transform.position + Vec3::new(1.0, 0.0, 0.0);
        let v = body.velocity_at_point(point);
        assert!((v.x - 1.0).abs() < 1e-10);
        assert!((v.y - 1.0).abs() < 1e-10);
        assert!(v.z.abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_roundtrip() {
        let mut body = RigidBody::new(3.0);
        body.velocity = Vec3::new(1.0, 2.0, 3.0);
        body.angular_velocity = Vec3::new(0.1, 0.2, 0.3);
        body.transform.position = Vec3::new(10.0, 20.0, 30.0);
        body.linear_damping = 0.05;
        body.angular_damping = 0.02;
        body.gravity_scale = 0.5;
        let snap = body.snapshot();
        assert_eq!(snap.body_type, "dynamic");
        assert_eq!(snap.state, "active");
        assert!((snap.mass - 3.0).abs() < 1e-10);
        let mut restored = RigidBody::new(1.0);
        restored.restore_from_snapshot(&snap);
        assert!((restored.mass - 3.0).abs() < 1e-10);
        assert!((restored.velocity.x - 1.0).abs() < 1e-10);
        assert!((restored.velocity.y - 2.0).abs() < 1e-10);
        assert!((restored.velocity.z - 3.0).abs() < 1e-10);
        assert!((restored.transform.position.x - 10.0).abs() < 1e-10);
        assert!((restored.linear_damping - 0.05).abs() < 1e-10);
    }
    #[test]
    fn test_merge_momentum_conservation() {
        let mut a = RigidBody::new(2.0);
        a.velocity = Vec3::new(3.0, 0.0, 0.0);
        let mut b = RigidBody::new(3.0);
        b.velocity = Vec3::new(-1.0, 0.0, 0.0);
        let merged = RigidBody::merge(&a, &b);
        assert!((merged.mass - 5.0).abs() < 1e-10);
        assert!((merged.velocity.x - 0.6).abs() < 1e-10);
    }
    #[test]
    fn test_merge_center_of_mass() {
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(2.0, 0.0, 0.0);
        let merged = RigidBody::merge(&a, &b);
        assert!((merged.transform.position.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_set_body_type_dynamic_to_kinematic() {
        let mut body = RigidBody::new(5.0);
        assert!(body.is_dynamic());
        body.set_body_type(BodyType::Kinematic);
        assert!(body.is_kinematic());
        assert!((body.inverse_mass - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_set_body_type_kinematic_to_dynamic() {
        let mut body = RigidBody::new_kinematic();
        body.mass = 5.0;
        body.set_body_type(BodyType::Dynamic);
        assert!(body.is_dynamic());
        assert!((body.inverse_mass - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_apply_angular_impulse() {
        let mut body = RigidBody::new(1.0);
        body.apply_angular_impulse(Vec3::new(0.0, 0.0, 2.0));
        assert!((body.angular_velocity.z - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_torque_accumulator() {
        let mut body = RigidBody::new(1.0);
        body.apply_torque(Vec3::new(0.0, 0.0, 5.0));
        body.apply_torque(Vec3::new(0.0, 0.0, 3.0));
        let t = body.accumulated_torque();
        assert!((t.z - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_clear_accumulators() {
        let mut body = RigidBody::new(1.0);
        body.apply_force(Vec3::new(10.0, 0.0, 0.0));
        body.apply_torque(Vec3::new(0.0, 5.0, 0.0));
        body.clear_accumulators();
        assert!(body.accumulated_force().norm() < 1e-10);
        assert!(body.accumulated_torque().norm() < 1e-10);
    }
    #[test]
    fn test_kinematic_target_motion() {
        let mut body = RigidBody::new_kinematic();
        let target_pos = Vec3::new(10.0, 0.0, 0.0);
        let target_rot = Quat::identity();
        body.set_kinematic_target(target_pos, target_rot);
        body.integrate_velocity(0.016);
        assert!((body.transform.position.x - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_kinematic_without_target_integrates() {
        let mut body = RigidBody::new_kinematic();
        body.velocity = Vec3::new(5.0, 0.0, 0.0);
        body.integrate_velocity(1.0);
        assert!((body.transform.position.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_angular_velocity_integration() {
        let mut body = RigidBody::new(1.0);
        body.angular_velocity = Vec3::new(0.0, 0.0, std::f64::consts::PI);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        body.gravity_scale = 0.0;
        body.integrate_angular_velocity(1.0);
        let angle = body.transform.rotation.angle();
        assert!(
            angle > 0.1,
            "Rotation angle should be non-trivial, got {angle}"
        );
    }
    #[test]
    fn test_force_at_point_produces_torque() {
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        body.gravity_scale = 0.0;
        let force = Vec3::new(0.0, 10.0, 0.0);
        let point = body.transform.position + Vec3::new(1.0, 0.0, 0.0);
        body.apply_force_at_point(force, point);
        assert!((body.torque_accumulator.z - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_step_count_increments() {
        let mut body = RigidBody::new(1.0);
        assert_eq!(body.step_count, 0);
        body.integrate_velocity(0.016);
        assert_eq!(body.step_count, 1);
        body.integrate_velocity(0.016);
        assert_eq!(body.step_count, 2);
    }
    #[test]
    fn test_world_aabb() {
        let mut body = RigidBody::new(1.0);
        body.transform.position = Vec3::new(5.0, 5.0, 5.0);
        let half = Vec3::new(1.0, 1.0, 1.0);
        let (min, max) = body.world_aabb(&half);
        assert!((min.x - 4.0).abs() < 1e-10);
        assert!((max.x - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_sleeping_body() {
        let mut body = RigidBody::new(1.0);
        body.state = BodyState::Sleeping;
        let snap = body.snapshot();
        assert_eq!(snap.state, "sleeping");
        let mut restored = RigidBody::new(1.0);
        restored.restore_from_snapshot(&snap);
        assert_eq!(restored.state, BodyState::Sleeping);
    }
    #[test]
    fn test_merge_zero_mass_bodies() {
        let a = RigidBody::new_static();
        let b = RigidBody::new_static();
        let merged = RigidBody::merge(&a, &b);
        assert_eq!(merged.body_type, BodyType::Static);
    }
    #[test]
    fn test_impulse_at_point_produces_angular_velocity() {
        let mut body = RigidBody::new(1.0);
        let impulse = Vec3::new(0.0, 1.0, 0.0);
        let point = body.transform.position + Vec3::new(1.0, 0.0, 0.0);
        body.apply_impulse_at_point(impulse, point);
        assert!((body.angular_velocity.z - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_is_kinematic() {
        let body = RigidBody::new_kinematic();
        assert!(body.is_kinematic());
        assert!(!body.is_dynamic());
        assert!(!body.is_static());
    }
    #[test]
    fn test_is_static() {
        let body = RigidBody::new_static();
        assert!(body.is_static());
        assert!(!body.is_dynamic());
        assert!(!body.is_kinematic());
    }
    #[test]
    fn test_force_accumulator_apply_force() {
        let mut acc = ForceAccumulator::new();
        acc.apply_force(Vec3::new(5.0, 0.0, 0.0));
        acc.apply_force(Vec3::new(3.0, 0.0, 0.0));
        assert!((acc.current_force.x - 8.0).abs() < 1e-10);
        assert_eq!(acc.total_force_applications, 2);
    }
    #[test]
    fn test_force_accumulator_peak_magnitude() {
        let mut acc = ForceAccumulator::new();
        acc.apply_force(Vec3::new(3.0, 4.0, 0.0));
        acc.apply_force(Vec3::new(1.0, 0.0, 0.0));
        assert!((acc.peak_force_magnitude - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_force_accumulator_flush() {
        let mut acc = ForceAccumulator::new();
        acc.apply_force(Vec3::new(10.0, 0.0, 0.0));
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        acc.flush(&mut body);
        assert!(acc.current_force.norm() < 1e-10);
        assert!((body.force_accumulator.x - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_force_accumulator_persistent_force() {
        let mut acc = ForceAccumulator::new();
        acc.add_persistent(PersistentForce::new("wind", Vec3::new(2.0, 0.0, 0.0)));
        assert_eq!(acc.persistent_count(), 1);
        assert_eq!(acc.enabled_count(), 1);
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        acc.flush(&mut body);
        assert!((body.force_accumulator.x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_force_accumulator_disable_persistent() {
        let mut acc = ForceAccumulator::new();
        acc.add_persistent(PersistentForce::new("thrust", Vec3::new(5.0, 0.0, 0.0)));
        acc.persistent[0].disable();
        assert_eq!(acc.enabled_count(), 0);
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        acc.flush(&mut body);
        assert!(body.force_accumulator.norm() < 1e-10);
    }
    #[test]
    fn test_force_accumulator_remove_persistent() {
        let mut acc = ForceAccumulator::new();
        acc.add_persistent(PersistentForce::new("a", Vec3::zeros()));
        acc.add_persistent(PersistentForce::new("b", Vec3::zeros()));
        acc.remove_persistent("a");
        assert_eq!(acc.persistent_count(), 1);
        assert_eq!(acc.persistent[0].label, "b");
    }
    #[test]
    fn test_persistent_force_at_point_produces_torque() {
        let mut acc = ForceAccumulator::new();
        let force = Vec3::new(0.0, 10.0, 0.0);
        let point = Vec3::new(1.0, 0.0, 0.0);
        acc.add_persistent(PersistentForce::at_point("lift", force, point));
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        acc.flush(&mut body);
        assert!((body.torque_accumulator.z - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_com_tracker_record_and_latest() {
        let mut tracker = CenterOfMassTracker::new(10);
        tracker.record(0.0, Vec3::new(0.0, 0.0, 0.0));
        tracker.record(1.0, Vec3::new(3.0, 4.0, 0.0));
        assert_eq!(tracker.sample_count(), 2);
        let latest = tracker.latest().unwrap();
        assert!((latest.x - 3.0).abs() < 1e-10);
        assert!((latest.y - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_com_tracker_total_path_length() {
        let mut tracker = CenterOfMassTracker::new(10);
        tracker.record(0.0, Vec3::new(0.0, 0.0, 0.0));
        tracker.record(1.0, Vec3::new(3.0, 0.0, 0.0));
        tracker.record(2.0, Vec3::new(3.0, 4.0, 0.0));
        let path = tracker.total_path_length();
        assert!(
            (path - 7.0).abs() < 1e-10,
            "path length should be 7, got {path}"
        );
    }
    #[test]
    fn test_com_tracker_average_speed() {
        let mut tracker = CenterOfMassTracker::new(10);
        tracker.record(0.0, Vec3::new(0.0, 0.0, 0.0));
        tracker.record(2.0, Vec3::new(6.0, 0.0, 0.0));
        let speed = tracker.average_speed();
        assert!(
            (speed - 3.0).abs() < 1e-10,
            "average speed should be 3, got {speed}"
        );
    }
    #[test]
    fn test_com_tracker_max_history() {
        let mut tracker = CenterOfMassTracker::new(3);
        for i in 0..6 {
            tracker.record(i as Real, Vec3::zeros());
        }
        assert_eq!(tracker.sample_count(), 3, "should be capped at 3");
    }
    #[test]
    fn test_com_tracker_clear() {
        let mut tracker = CenterOfMassTracker::new(10);
        tracker.record(0.0, Vec3::zeros());
        tracker.clear();
        assert_eq!(tracker.sample_count(), 0);
        assert!(tracker.latest().is_none());
    }
    #[test]
    fn test_com_tracker_no_samples() {
        let tracker = CenterOfMassTracker::new(10);
        assert!(tracker.latest().is_none());
        assert!((tracker.total_path_length() - 0.0).abs() < 1e-10);
        assert!((tracker.average_speed() - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_constraint_list_add_remove() {
        let mut list = ConstraintList::new();
        let id1 = list.add(ConstraintType::RotationLock, "no spin");
        let _id2 = list.add(ConstraintType::PositionLock, "no move");
        assert_eq!(list.len(), 2);
        let removed = list.remove(id1);
        assert!(removed);
        assert_eq!(list.len(), 1);
        assert!(!list.remove(id1));
    }
    #[test]
    fn test_constraint_list_get() {
        let mut list = ConstraintList::new();
        let id = list.add(ConstraintType::Fixed, "fixed");
        let c = list.get(id).unwrap();
        assert_eq!(c.constraint_type, ConstraintType::Fixed);
        assert!(c.active);
    }
    #[test]
    fn test_constraint_rotation_lock_zeroes_angular_velocity() {
        let mut list = ConstraintList::new();
        list.add(ConstraintType::RotationLock, "no spin");
        let mut body = RigidBody::new(1.0);
        body.angular_velocity = Vec3::new(5.0, 3.0, 1.0);
        list.apply_to(&mut body);
        assert!(body.angular_velocity.norm() < 1e-10);
    }
    #[test]
    fn test_constraint_position_lock_zeroes_velocity() {
        let mut list = ConstraintList::new();
        list.add(ConstraintType::PositionLock, "frozen");
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(3.0, 4.0, 5.0);
        list.apply_to(&mut body);
        assert!(body.velocity.norm() < 1e-10);
        body.angular_velocity = Vec3::new(1.0, 0.0, 0.0);
        list.apply_to(&mut body);
        assert!((body.angular_velocity.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_constraint_fixed_zeroes_all() {
        let mut list = ConstraintList::new();
        list.add(ConstraintType::Fixed, "full lock");
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(1.0, 2.0, 3.0);
        body.angular_velocity = Vec3::new(4.0, 5.0, 6.0);
        list.apply_to(&mut body);
        assert!(body.velocity.norm() < 1e-10);
        assert!(body.angular_velocity.norm() < 1e-10);
    }
    #[test]
    fn test_constraint_axis_lock_x() {
        let mut list = ConstraintList::new();
        list.add(ConstraintType::AxisLock(0), "lock x");
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(5.0, 3.0, 1.0);
        list.apply_to(&mut body);
        assert!(body.velocity.x.abs() < 1e-10);
        assert!((body.velocity.y - 3.0).abs() < 1e-10);
        assert!((body.velocity.z - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_constraint_deactivate() {
        let mut list = ConstraintList::new();
        let id = list.add(ConstraintType::Fixed, "lock");
        list.constraints[0].deactivate();
        let active = list.active_constraints();
        assert!(active.is_empty());
        list.constraints[0].activate();
        assert_eq!(list.active_constraints().len(), 1);
        let _ = id;
    }
    #[test]
    fn test_contact_history_record_and_latest() {
        let mut history = ContactHistory::new(10);
        history.record(0.1, 5, Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 3.0);
        history.record(0.2, 7, Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 5.0);
        assert_eq!(history.count(), 2);
        assert!((history.latest().unwrap().impulse - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_history_peak_impulse() {
        let mut history = ContactHistory::new(10);
        history.record(0.1, 0, Vec3::zeros(), Vec3::zeros(), 2.0);
        history.record(0.2, 1, Vec3::zeros(), Vec3::zeros(), 8.0);
        history.record(0.3, 2, Vec3::zeros(), Vec3::zeros(), 4.0);
        assert!((history.peak_impulse() - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_history_average_impulse() {
        let mut history = ContactHistory::new(10);
        history.record(0.0, 0, Vec3::zeros(), Vec3::zeros(), 4.0);
        history.record(0.0, 1, Vec3::zeros(), Vec3::zeros(), 6.0);
        assert!((history.average_impulse() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_history_contacts_with() {
        let mut history = ContactHistory::new(10);
        history.record(0.0, 3, Vec3::zeros(), Vec3::zeros(), 1.0);
        history.record(0.0, 7, Vec3::zeros(), Vec3::zeros(), 2.0);
        history.record(0.0, 3, Vec3::zeros(), Vec3::zeros(), 3.0);
        let with_3 = history.contacts_with(3);
        assert_eq!(with_3.len(), 2);
    }
    #[test]
    fn test_contact_history_max_entries() {
        let mut history = ContactHistory::new(3);
        for i in 0..6 {
            history.record(i as Real, i as u32, Vec3::zeros(), Vec3::zeros(), 1.0);
        }
        assert_eq!(history.count(), 3);
    }
    #[test]
    fn test_contact_history_clear() {
        let mut history = ContactHistory::new(10);
        history.record(0.0, 0, Vec3::zeros(), Vec3::zeros(), 1.0);
        history.clear();
        assert!(history.is_empty());
    }
    #[test]
    fn test_state_interpolator_push_and_sample_count() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.0, Vec3::zeros(), Quat::identity());
        interp.push(1.0, Vec3::new(1.0, 0.0, 0.0), Quat::identity());
        assert_eq!(interp.sample_count(), 2);
    }
    #[test]
    fn test_state_interpolator_interpolate_position_midpoint() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.0, Vec3::new(0.0, 0.0, 0.0), Quat::identity());
        interp.push(2.0, Vec3::new(4.0, 0.0, 0.0), Quat::identity());
        let pos = interp
            .interpolate_position(1.0)
            .expect("should interpolate");
        assert!(
            (pos.x - 2.0).abs() < 1e-10,
            "midpoint x should be 2.0, got {}",
            pos.x
        );
    }
    #[test]
    fn test_state_interpolator_position_at_boundary() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.0, Vec3::new(0.0, 0.0, 0.0), Quat::identity());
        interp.push(1.0, Vec3::new(10.0, 0.0, 0.0), Quat::identity());
        let p0 = interp.interpolate_position(0.0).unwrap();
        assert!(p0.x.abs() < 1e-10);
        let p1 = interp.interpolate_position(1.0).unwrap();
        assert!((p1.x - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_state_interpolator_out_of_range() {
        let mut interp = StateInterpolator::new(5);
        interp.push(1.0, Vec3::zeros(), Quat::identity());
        interp.push(2.0, Vec3::new(1.0, 0.0, 0.0), Quat::identity());
        assert!(interp.interpolate_position(0.0).is_none());
        assert!(interp.interpolate_position(3.0).is_none());
    }
    #[test]
    fn test_state_interpolator_too_few_samples() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.0, Vec3::zeros(), Quat::identity());
        assert!(interp.interpolate_position(0.0).is_none());
    }
    #[test]
    fn test_state_interpolator_time_range() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.5, Vec3::zeros(), Quat::identity());
        interp.push(1.5, Vec3::zeros(), Quat::identity());
        interp.push(3.0, Vec3::zeros(), Quat::identity());
        let (t_min, t_max) = interp.time_range().unwrap();
        assert!((t_min - 0.5).abs() < 1e-10);
        assert!((t_max - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_state_interpolator_push_from_body() {
        let body = RigidBody::new(1.0);
        let mut interp = StateInterpolator::new(5);
        interp.push_from_body(0.0, &body);
        assert_eq!(interp.sample_count(), 1);
    }
    #[test]
    fn test_state_interpolator_clear() {
        let mut interp = StateInterpolator::new(5);
        interp.push(0.0, Vec3::zeros(), Quat::identity());
        interp.clear();
        assert_eq!(interp.sample_count(), 0);
        assert!(interp.time_range().is_none());
    }
    #[test]
    fn test_state_interpolator_rotation() {
        let mut interp = StateInterpolator::new(5);
        let r0 = Quat::identity();
        let r1 = Quat::from_axis_angle(
            &Unit::new_normalize(Vec3::new(0.0, 1.0, 0.0)),
            std::f64::consts::PI,
        );
        interp.push(0.0, Vec3::zeros(), r0);
        interp.push(1.0, Vec3::zeros(), r1);
        let mid_rot = interp.interpolate_rotation(0.5);
        assert!(mid_rot.is_some(), "should produce a rotation at t=0.5");
    }
    #[test]
    fn test_state_interpolator_max_samples() {
        let mut interp = StateInterpolator::new(3);
        for i in 0..6 {
            interp.push(i as Real, Vec3::zeros(), Quat::identity());
        }
        assert_eq!(interp.sample_count(), 3);
    }
}
#[cfg(test)]
mod new_body_tests {

    use crate::AerodynamicParams;

    use crate::BuoyancyParams;

    use crate::CompoundPrimitive;
    use crate::CompoundRigidBody;

    use crate::MagnusEffect;

    use crate::RigidBody;
    use crate::RollingConstraint;

    use crate::VariableMassBody;

    use oxiphysics_core::Real;

    #[test]
    fn compound_body_total_mass_is_sum_of_parts() {
        let prims = vec![
            CompoundPrimitive::new_box([0.0, 0.0, 0.0], [0.5, 0.5, 0.5], 1000.0),
            CompoundPrimitive::new_sphere([1.0, 0.0, 0.0], 0.3, 800.0),
        ];
        let total = prims.iter().map(|p| p.mass).sum::<Real>();
        let compound = CompoundRigidBody::from_primitives(prims);
        assert!(
            (compound.total_mass - total).abs() < 1e-6,
            "total_mass={}",
            compound.total_mass
        );
    }
    #[test]
    fn compound_body_is_valid() {
        let prims = vec![CompoundPrimitive::new_box([0.0; 3], [1.0, 1.0, 1.0], 100.0)];
        let compound = CompoundRigidBody::from_primitives(prims);
        assert!(compound.is_valid());
    }
    #[test]
    fn compound_body_com_at_origin_symmetric() {
        let prims = vec![
            CompoundPrimitive::new_box([-1.0, 0.0, 0.0], [0.5, 0.5, 0.5], 1000.0),
            CompoundPrimitive::new_box([1.0, 0.0, 0.0], [0.5, 0.5, 0.5], 1000.0),
        ];
        let compound = CompoundRigidBody::from_primitives(prims);
        assert!(
            compound.local_com[0].abs() < 1e-10,
            "CoM x should be 0: {}",
            compound.local_com[0]
        );
    }
    #[test]
    fn compound_body_inertia_positive_diagonal() {
        let prims = vec![
            CompoundPrimitive::new_sphere([0.0, 0.0, 0.0], 0.5, 1000.0),
            CompoundPrimitive::new_box([1.0, 0.0, 0.0], [0.3, 0.3, 0.3], 500.0),
        ];
        let compound = CompoundRigidBody::from_primitives(prims);
        for &i in &compound.inertia_diag {
            assert!(i > 0.0, "inertia_diag component should be positive: {i}");
        }
    }
    #[test]
    fn compound_to_rigid_body_has_correct_mass() {
        let prims = vec![CompoundPrimitive::new_box(
            [0.0; 3],
            [0.5, 0.5, 0.5],
            1000.0,
        )];
        let compound = CompoundRigidBody::from_primitives(prims);
        let body = compound.to_rigid_body();
        assert!((body.mass - compound.total_mass).abs() < 1e-6);
    }
    #[test]
    fn compound_primitive_box_mass_correct() {
        let prim = CompoundPrimitive::new_box([0.0; 3], [0.5, 0.5, 0.5], 1000.0);
        assert!((prim.mass - 1000.0).abs() < 1e-6, "mass={}", prim.mass);
    }
    #[test]
    fn compound_cylinder_primitive_mass_nonzero() {
        let prim = CompoundPrimitive::new_cylinder([0.0; 3], 0.5, 1.0, 1000.0);
        assert!(prim.mass > 0.0, "mass={}", prim.mass);
    }
    #[test]
    fn compound_primitive_count() {
        let prims = vec![
            CompoundPrimitive::new_box([0.0; 3], [0.5; 3], 100.0),
            CompoundPrimitive::new_sphere([1.0, 0.0, 0.0], 0.3, 100.0),
            CompoundPrimitive::new_cylinder([0.0, 2.0, 0.0], 0.2, 0.5, 100.0),
        ];
        let compound = CompoundRigidBody::from_primitives(prims);
        assert_eq!(compound.primitive_count(), 3);
    }
    #[test]
    fn variable_mass_tsiolkovsky_positive() {
        let body = VariableMassBody::new(10.0, 90.0, 300.0, 1000.0, [0.0, 1.0, 0.0]);
        let dv = body.tsiolkovsky_delta_v(10.0);
        assert!(dv > 0.0, "delta-v should be positive: {dv}");
        assert!(dv > 5000.0, "should be large: {dv}");
    }
    #[test]
    fn variable_mass_remaining_propellant_decreases() {
        let mut body = VariableMassBody::new(10.0, 100.0, 300.0, 1000.0, [0.0, 1.0, 0.0]);
        body.throttle = 1.0;
        let before = body.remaining_propellant();
        body.step(0.1, &[0.0, -9.81, 0.0]);
        let after = body.remaining_propellant();
        assert!(
            after < before,
            "propellant should decrease: before={before}, after={after}"
        );
    }
    #[test]
    fn variable_mass_no_consumption_when_throttle_zero() {
        let mut body = VariableMassBody::new(10.0, 100.0, 300.0, 1000.0, [0.0, 1.0, 0.0]);
        body.throttle = 0.0;
        let before = body.remaining_propellant();
        body.step(1.0, &[0.0, -9.81, 0.0]);
        assert!((body.remaining_propellant() - before).abs() < 1e-10);
    }
    #[test]
    fn variable_mass_flow_rate_positive_at_throttle_one() {
        let body = VariableMassBody::new(10.0, 100.0, 300.0, 1000.0, [0.0, 1.0, 0.0]);
        let mut b2 = body;
        b2.throttle = 1.0;
        assert!(b2.mass_flow_rate() > 0.0);
    }
    #[test]
    fn variable_mass_flow_rate_zero_at_throttle_zero() {
        let mut body = VariableMassBody::new(10.0, 100.0, 300.0, 1000.0, [0.0, 1.0, 0.0]);
        body.throttle = 0.0;
        assert_eq!(body.mass_flow_rate(), 0.0);
    }
    #[test]
    fn variable_mass_current_mass_decreases_during_burn() {
        let mut body = VariableMassBody::new(10.0, 50.0, 250.0, 500.0, [0.0, 1.0, 0.0]);
        body.throttle = 1.0;
        let m0 = body.current_mass();
        for _ in 0..10 {
            body.step(0.1, &[0.0, 0.0, 0.0]);
        }
        assert!(body.current_mass() < m0, "mass should decrease during burn");
    }
    #[test]
    fn buoyancy_force_upward() {
        let params = BuoyancyParams::water();
        let [fx, fy, fz] = params.buoyant_force(0.001);
        assert!(fy > 0.0, "buoyancy should be upward: {fy}");
        assert_eq!(fx, 0.0);
        assert_eq!(fz, 0.0);
    }
    #[test]
    fn buoyancy_archimedes_principle() {
        let params = BuoyancyParams::water();
        let [_, fy, _] = params.buoyant_force(1e-3);
        assert!((fy - 1000.0 * 9.81 * 1e-3).abs() < 1e-4, "fy={fy}");
    }
    #[test]
    fn buoyancy_net_force_negative_when_heavy() {
        let params = BuoyancyParams::water();
        let net = params.net_vertical_force(10.0, 5e-4);
        assert!(net < 0.0, "heavy object sinks: net={net}");
    }
    #[test]
    fn buoyancy_net_force_positive_when_light() {
        let params = BuoyancyParams::water();
        let net = params.net_vertical_force(0.1, 1e-3);
        assert!(net > 0.0, "light object floats: net={net}");
    }
    #[test]
    fn buoyancy_apply_to_body() {
        let params = BuoyancyParams::water();
        let mut body = RigidBody::new(1.0);
        body.linear_damping = 0.0;
        body.gravity_scale = 0.0;
        params.apply_to_body(&mut body, 1e-3);
        assert!(
            body.force_accumulator.y > 0.0,
            "buoyancy should add upward force"
        );
    }
    #[test]
    fn drag_force_opposes_velocity() {
        let params = AerodynamicParams::standard();
        let v = [10.0, 0.0, 0.0];
        let [fx, _, _] = params.drag_force(v);
        assert!(fx < 0.0, "drag should oppose x-velocity: {fx}");
    }
    #[test]
    fn drag_force_zero_at_rest() {
        let params = AerodynamicParams::standard();
        let [fx, fy, fz] = params.drag_force([0.0, 0.0, 0.0]);
        assert_eq!(fx, 0.0);
        assert_eq!(fy, 0.0);
        assert_eq!(fz, 0.0);
    }
    #[test]
    fn drag_force_quadratic_in_speed() {
        let params = AerodynamicParams::standard();
        let [f1, _, _] = params.drag_force([10.0, 0.0, 0.0]);
        let [f2, _, _] = params.drag_force([20.0, 0.0, 0.0]);
        let ratio = f2 / f1;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "drag should be quadratic: ratio={ratio}"
        );
    }
    #[test]
    fn drag_dynamic_pressure() {
        let mut params = AerodynamicParams::standard();
        params.wind_velocity = [10.0, 0.0, 0.0];
        let q = params.dynamic_pressure();
        assert!((q - 61.25).abs() < 0.01, "q={q}");
    }
    #[test]
    fn magnus_force_perpendicular_to_velocity_and_spin() {
        let magnus = MagnusEffect::sphere(0.1, 1.225);
        let f = magnus.force([0.0, 0.0, 10.0], [10.0, 0.0, 0.0]);
        assert!(
            f[1].abs() > 0.0,
            "Magnus force should have y component: {:?}",
            f
        );
    }
    #[test]
    fn magnus_force_zero_no_spin() {
        let magnus = MagnusEffect::sphere(0.1, 1.225);
        let f = magnus.force([0.0, 0.0, 0.0], [10.0, 0.0, 0.0]);
        assert!(
            f.iter().all(|&x| x.abs() < 1e-14),
            "no spin → no Magnus: {f:?}"
        );
    }
    #[test]
    fn magnus_force_zero_no_velocity() {
        let magnus = MagnusEffect::sphere(0.1, 1.225);
        let f = magnus.force([10.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(
            f.iter().all(|&x| x.abs() < 1e-14),
            "no velocity → no Magnus: {f:?}"
        );
    }
    #[test]
    fn magnus_force_scales_with_fluid_density() {
        let m1 = MagnusEffect::sphere(0.1, 1.0);
        let m2 = MagnusEffect::sphere(0.1, 2.0);
        let f1 = m1.force([0.0, 0.0, 10.0], [10.0, 0.0, 0.0]);
        let f2 = m2.force([0.0, 0.0, 10.0], [10.0, 0.0, 0.0]);
        assert!(
            (f2[1] / f1[1] - 2.0).abs() < 1e-10,
            "force should scale with density"
        );
    }
    #[test]
    fn rolling_no_slip_after_enforce() {
        let rc = RollingConstraint::sphere(0.5);
        let mut body = RigidBody::new(1.0);
        body.velocity = oxiphysics_core::math::Vec3::new(2.0, 0.0, 0.0);
        body.angular_velocity = oxiphysics_core::math::Vec3::zeros();
        rc.enforce(&mut body);
        let slip = rc.slip_velocity(&body);
        let slip_speed = (slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]).sqrt();
        assert!(
            slip_speed < 2.0,
            "slip should decrease after enforcement: {slip_speed}"
        );
    }
    #[test]
    fn rolling_friction_torque_positive() {
        let rc = RollingConstraint::sphere(0.3);
        let tau = rc.rolling_friction_torque(1.0, 9.81);
        assert!(
            tau > 0.0,
            "rolling friction torque should be positive: {tau}"
        );
    }
    #[test]
    fn rolling_slip_zero_at_correct_omega() {
        let rc = RollingConstraint::sphere(1.0);
        let mut body = RigidBody::new(1.0);
        body.velocity = oxiphysics_core::math::Vec3::new(2.0, 0.0, 0.0);
        body.angular_velocity = oxiphysics_core::math::Vec3::new(0.0, 0.0, -2.0);
        let slip = rc.slip_velocity(&body);
        let slip_speed = (slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2]).sqrt();
        assert!(
            slip_speed < 0.01,
            "rolling without slipping: slip={slip_speed}"
        );
    }
    #[test]
    fn rolling_static_body_unchanged() {
        let rc = RollingConstraint::sphere(0.5);
        let mut body = RigidBody::new_static();
        let v_before = body.velocity;
        let w_before = body.angular_velocity;
        rc.enforce(&mut body);
        assert!((body.velocity - v_before).norm() < 1e-14);
        assert!((body.angular_velocity - w_before).norm() < 1e-14);
    }
}
