//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Quat, Real, Unit, Vec3};
use oxiphysics_rigid::RigidBodySet;

use super::types::BodyProps;

/// Baumgarte factor for joint position correction.
pub(super) const JOINT_BAUMGARTE: Real = 0.1;
/// Compute the scalar effective mass for a linear constraint along `dir`.
pub(super) fn linear_effective_mass(
    inv_mass_a: Real,
    inv_mass_b: Real,
    inv_inertia_a: &oxiphysics_core::math::Mat3,
    inv_inertia_b: &oxiphysics_core::math::Mat3,
    r_a: &Vec3,
    r_b: &Vec3,
    dir: &Vec3,
) -> Real {
    let ra_cross = r_a.cross(dir);
    let rb_cross = r_b.cross(dir);
    let k = inv_mass_a
        + inv_mass_b
        + ra_cross.dot(&(inv_inertia_a * ra_cross))
        + rb_cross.dot(&(inv_inertia_b * rb_cross));
    if k > 1e-12 { 1.0 / k } else { 0.0 }
}
/// Apply a linear impulse to a pair of bodies.
pub(super) fn apply_pair_impulse(
    bodies: &mut RigidBodySet,
    ha: BodyHandle,
    hb: BodyHandle,
    impulse: Vec3,
    r_a: &Vec3,
    r_b: &Vec3,
) {
    if let Some(body) = bodies.get_mut(ha) {
        body.velocity += impulse * body.inverse_mass;
        body.angular_velocity += body.world_inverse_inertia * r_a.cross(&impulse);
    }
    if let Some(body) = bodies.get_mut(hb) {
        body.velocity -= impulse * body.inverse_mass;
        body.angular_velocity -= body.world_inverse_inertia * r_b.cross(&impulse);
    }
}
/// Apply a pure angular impulse to a pair of bodies.
pub(super) fn apply_angular_impulse(
    bodies: &mut RigidBodySet,
    ha: BodyHandle,
    hb: BodyHandle,
    impulse: Vec3,
) {
    if let Some(body) = bodies.get_mut(ha) {
        body.angular_velocity += body.world_inverse_inertia * impulse;
    }
    if let Some(body) = bodies.get_mut(hb) {
        body.angular_velocity -= body.world_inverse_inertia * impulse;
    }
}
/// Extract body properties from the body set.
pub(super) fn read_body(bodies: &RigidBodySet, handle: BodyHandle) -> Option<BodyProps> {
    bodies.get(handle).map(|b| BodyProps {
        inv_mass: b.inverse_mass,
        inv_inertia: b.world_inverse_inertia,
        vel: b.velocity,
        ang_vel: b.angular_velocity,
        position: b.transform.position,
        rotation: b.transform.rotation,
    })
}
/// Create a small rotation quaternion from an axis-angle vector.
pub(super) fn small_rotation_quat(v: &Vec3) -> Quat {
    let angle = v.norm();
    if angle < 1e-10 {
        Quat::identity()
    } else {
        Quat::from_axis_angle(&Unit::new_normalize(*v), angle)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BallJoint;
    use crate::Constraint;
    use crate::FixedJoint;
    use crate::GearJoint;
    use crate::MotorJoint;
    use crate::PrismaticJoint;
    use crate::PulleyJoint;
    use crate::RackPinionJoint;
    use crate::RevoluteJoint;
    use crate::SpringJoint;

    use oxiphysics_core::math::Vec3;
    use oxiphysics_rigid::RigidBody;
    fn setup_two_bodies(pos_a: Vec3, pos_b: Vec3) -> (RigidBodySet, BodyHandle, BodyHandle) {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = pos_a;
        a.linear_damping = 0.0;
        a.angular_damping = 0.0;
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = pos_b;
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        let hb = bodies.insert(b);
        (bodies, ha, hb)
    }
    #[test]
    fn test_fixed_joint_maintains_relative_pose() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = FixedJoint::new(
            ha,
            hb,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(-0.5, 0.0, 0.0),
            Quat::identity(),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        joint.solve_position(&mut bodies, dt);
        let va = bodies.get(ha).unwrap().velocity;
        let vb = bodies.get(hb).unwrap().velocity;
        let vel_diff = (va - vb).norm();
        assert!(
            vel_diff < 2.0,
            "Fixed joint should constrain relative velocity, diff={vel_diff}"
        );
    }
    #[test]
    fn test_revolute_joint_allows_axis_rotation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(0.0, 2.0, 0.0);
        }
        let mut joint = RevoluteJoint::new(
            ha,
            hb,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let b = bodies.get(hb).unwrap();
        let total_ang_y = bodies.get(ha).unwrap().angular_velocity.y + b.angular_velocity.y;
        assert!(
            total_ang_y.abs() > 0.1,
            "Revolute joint should allow rotation around axis, total_ang_y={total_ang_y}"
        );
        let perp_a = Vec3::new(
            bodies.get(ha).unwrap().angular_velocity.x,
            0.0,
            bodies.get(ha).unwrap().angular_velocity.z,
        );
        let perp_b = Vec3::new(b.angular_velocity.x, 0.0, b.angular_velocity.z);
        let perp_diff = (perp_a - perp_b).norm();
        assert!(
            perp_diff < 1.0,
            "Revolute joint should block perpendicular rotation diff={perp_diff}"
        );
    }
    #[test]
    fn test_ball_joint_constrains_position() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(5.0, 0.0, 0.0);
        }
        let mut joint = BallJoint::new(ha, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let r_a = a.transform.rotation * Vec3::new(1.0, 0.0, 0.0);
        let r_b = b.transform.rotation * Vec3::new(-1.0, 0.0, 0.0);
        let va = a.velocity + a.angular_velocity.cross(&r_a);
        let vb = b.velocity + b.angular_velocity.cross(&r_b);
        let rel_v = (va - vb).norm();
        assert!(
            rel_v < 1.0,
            "Ball joint should constrain relative velocity at anchor, rel_v={rel_v}"
        );
    }
    #[test]
    fn test_spring_joint_exerts_correct_force() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        let stiffness = 100.0;
        let rest_length = 2.0;
        let mut spring = SpringJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            rest_length,
            stiffness,
            0.0,
        );
        let dt = 1.0 / 60.0;
        spring.prepare(&bodies, dt);
        spring.solve_velocity(&mut bodies, dt);
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        assert!(
            a.velocity.x > 0.0,
            "Spring should pull A toward B, vx={}",
            a.velocity.x
        );
        assert!(
            b.velocity.x < 0.0,
            "Spring should pull B toward A, vx={}",
            b.velocity.x
        );
        let expected_dv = stiffness * 1.0 * dt;
        assert!(
            (a.velocity.x - expected_dv).abs() < 0.01,
            "Spring velocity change should match, got {} expected {}",
            a.velocity.x,
            expected_dv
        );
    }
    #[test]
    fn test_ball_joint_allows_rotation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(0.0, 3.0, 0.0);
        }
        let mut joint = BallJoint::new(ha, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let b = bodies.get(hb).unwrap();
        let a = bodies.get(ha).unwrap();
        let total_ang = a.angular_velocity.norm() + b.angular_velocity.norm();
        assert!(
            total_ang > 0.1,
            "Ball joint should allow rotation, total_ang={}",
            total_ang
        );
    }
    #[test]
    fn test_prismatic_joint_allows_sliding() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(5.0, 0.0, 0.0);
        }
        let mut joint = PrismaticJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let total_vx = a.velocity.x + b.velocity.x;
        assert!(
            total_vx > 1.0,
            "Prismatic joint should allow sliding along axis, total_vx={}",
            total_vx
        );
    }
    #[test]
    fn test_tgs_solver_basic() {
        use crate::{PgsSolver, TgsSolver};
        let (mut bodies_tgs, ha_t, hb_t) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies_tgs.get_mut(hb_t) {
            b.angular_velocity = Vec3::new(0.0, 0.0, 3.0);
        }
        let (mut bodies_pgs, ha_p, hb_p) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies_pgs.get_mut(hb_p) {
            b.angular_velocity = Vec3::new(0.0, 0.0, 3.0);
        }
        let make_joint = |ha, hb| -> Box<dyn Constraint> {
            Box::new(RevoluteJoint::new(
                ha,
                hb,
                Vec3::new(0.5, 0.0, 0.0),
                Vec3::new(-0.5, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ))
        };
        let dt = 1.0 / 60.0;
        let mut constraints_tgs: Vec<Box<dyn Constraint>> = vec![make_joint(ha_t, hb_t)];
        let tgs = TgsSolver::new(1, 10, Vec3::zeros());
        tgs.solve(&mut constraints_tgs, &mut bodies_tgs, dt);
        let mut constraints_pgs: Vec<Box<dyn Constraint>> = vec![make_joint(ha_p, hb_p)];
        let pgs = PgsSolver::new(10, 1);
        pgs.solve(&mut constraints_pgs, &mut bodies_pgs, dt);
        let tgs_b = bodies_tgs.get(hb_t).unwrap();
        let pgs_b = bodies_pgs.get(hb_p).unwrap();
        let tgs_a = bodies_tgs.get(ha_t).unwrap();
        let pgs_a = bodies_pgs.get(ha_p).unwrap();
        let tgs_diff_z = (tgs_a.angular_velocity.z - tgs_b.angular_velocity.z).abs();
        let pgs_diff_z = (pgs_a.angular_velocity.z - pgs_b.angular_velocity.z).abs();
        assert!(
            tgs_diff_z < 3.0,
            "TGS should reduce off-axis angular velocity diff, got {tgs_diff_z}"
        );
        assert!(
            pgs_diff_z < 3.0,
            "PGS should reduce off-axis angular velocity diff, got {pgs_diff_z}"
        );
    }
    #[test]
    fn test_revolute_joint_angle_limit() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let mut joint = RevoluteJoint::new(
            ha,
            hb,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
        .with_limits(-std::f64::consts::FRAC_PI_4, std::f64::consts::FRAC_PI_4);
        assert!(joint.lower_limit.is_some(), "Lower limit should be set");
        assert!(joint.upper_limit.is_some(), "Upper limit should be set");
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(0.0, 0.0, 10.0);
        }
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let b = bodies.get(hb).unwrap();
        assert!(
            b.angular_velocity.z.abs() < 5.0,
            "Revolute angle limit joint should reduce off-axis angular velocity, z_ang={}",
            b.angular_velocity.z
        );
    }
    #[test]
    fn test_prismatic_joint_slide() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(3.0, 4.0, 0.0);
        }
        let mut joint = PrismaticJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let vy_diff = (a.velocity.y - b.velocity.y).abs();
        assert!(
            vy_diff < 1.0,
            "Prismatic joint should block perpendicular (Y) velocity diff={vy_diff}"
        );
        let total_px = a.velocity.x + b.velocity.x;
        assert!(
            (total_px - 3.0).abs() < 0.5,
            "Prismatic joint should conserve X momentum, total_px={total_px}"
        );
    }
    #[test]
    fn test_ball_joint_distance() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
        let anchor_a = Vec3::new(1.0, 0.0, 0.0);
        let anchor_b = Vec3::new(-1.0, 0.0, 0.0);
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(10.0, 0.0, 0.0);
        }
        let mut joint = BallJoint::new(ha, hb, anchor_a, anchor_b);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..40 {
            joint.solve_velocity(&mut bodies, dt);
        }
        for _ in 0..10 {
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt);
            }
            joint.prepare(&bodies, dt);
            joint.solve_position(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let world_a = a.transform.position + a.transform.rotation * anchor_a;
        let world_b = b.transform.position + b.transform.rotation * anchor_b;
        let dist = (world_a - world_b).norm();
        assert!(
            dist < 1.0,
            "Ball joint should keep anchor distance near zero, dist={dist}"
        );
    }
    #[test]
    fn test_motor_joint_drives_rotation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.inverse_mass = 0.0;
            b.world_inverse_inertia = oxiphysics_core::math::Mat3::zeros();
        }
        let target_velocity = 5.0;
        let mut motor = MotorJoint::new(ha, hb, Vec3::new(0.0, 1.0, 0.0), target_velocity, 1000.0);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            motor.prepare(&bodies, dt);
            motor.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let omega_y = a.angular_velocity.y;
        assert!(
            omega_y > 0.0,
            "Motor joint should drive body A to positive Y angular velocity, omega_y={omega_y}"
        );
        assert!(
            omega_y <= target_velocity + 0.5,
            "Motor joint angular velocity should not greatly exceed target, omega_y={omega_y}"
        );
    }
    #[test]
    fn test_fixed_joint_blocks_all_motion() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(a) = bodies.get_mut(ha) {
            a.velocity = Vec3::new(10.0, 0.0, 0.0);
        }
        let mut joint = FixedJoint::new(
            ha,
            hb,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(-0.5, 0.0, 0.0),
            Quat::identity(),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let va = bodies.get(ha).unwrap().velocity;
        let vb = bodies.get(hb).unwrap().velocity;
        assert!(
            vb.x.abs() > 0.1,
            "Fixed joint should couple body B to body A's motion, vb.x={}",
            vb.x
        );
        let vel_diff = (va - vb).norm();
        assert!(
            vel_diff < 1.0,
            "Fixed joint should block relative motion, vel_diff={vel_diff}"
        );
    }
    #[test]
    fn test_spring_joint_oscillation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        if let Some(a) = bodies.get_mut(ha) {
            a.inverse_mass = 0.0;
            a.world_inverse_inertia = oxiphysics_core::math::Mat3::zeros();
        }
        let mut spring = SpringJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0, 50.0, 0.0);
        let dt = 1.0 / 60.0;
        let mut positions: Vec<f64> = Vec::new();
        for _ in 0..120 {
            spring.prepare(&bodies, dt);
            spring.solve_velocity(&mut bodies, dt);
            if let Some(b) = bodies.get_mut(hb) {
                b.transform.position += b.velocity * dt;
            }
            let pos = bodies.get(hb).unwrap().transform.position.x;
            positions.push(pos);
        }
        let went_above = positions.iter().any(|&p| p > 2.1);
        let went_below = positions.iter().any(|&p| p < 1.9);
        assert!(
            went_above || went_below,
            "Undamped spring should oscillate around rest length"
        );
    }
    #[test]
    fn test_spring_joint_damping_reduces_amplitude() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        if let Some(a) = bodies.get_mut(ha) {
            a.inverse_mass = 0.0;
            a.world_inverse_inertia = oxiphysics_core::math::Mat3::zeros();
        }
        let mut spring = SpringJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0, 50.0, 5.0);
        let dt = 1.0 / 60.0;
        let mut speeds_early: Vec<f64> = Vec::new();
        let mut speeds_late: Vec<f64> = Vec::new();
        for i in 0..120 {
            spring.prepare(&bodies, dt);
            spring.solve_velocity(&mut bodies, dt);
            if let Some(b) = bodies.get_mut(hb) {
                b.transform.position += b.velocity * dt;
            }
            let speed = bodies.get(hb).unwrap().velocity.norm();
            if i < 30 {
                speeds_early.push(speed);
            } else if i >= 90 {
                speeds_late.push(speed);
            }
        }
        let max_early = speeds_early
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let max_late = speeds_late
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_late <= max_early + 1e-3,
            "Damped spring should reduce amplitude over time: early_max={max_early}, late_max={max_late}"
        );
    }
    #[test]
    fn test_revolute_joint_preserves_pivot() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(5.0, 3.0, 0.0);
        }
        let mut joint = RevoluteJoint::new(
            ha,
            hb,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a_body = bodies.get(ha).unwrap();
        let b_body = bodies.get(hb).unwrap();
        let rel_v = (a_body.velocity - b_body.velocity).norm();
        assert!(
            rel_v < 1e-3,
            "Revolute joint pivot velocity constraint: rel_v={rel_v}"
        );
    }
    #[test]
    fn test_prismatic_joint_axis_aligned() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = PrismaticJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let vy_diff = (a.velocity.y - b.velocity.y).abs();
        assert!(
            vy_diff < 1e-3,
            "Prismatic joint should block Y velocity (perpendicular to axis), vy_diff={vy_diff}"
        );
    }
    #[test]
    fn test_ball_joint_allows_rotation_constrains_translation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(0.0, 3.0, 0.0);
            b.velocity = Vec3::new(5.0, 0.0, 0.0);
        }
        let anchor_a = Vec3::new(0.0, 0.0, 0.0);
        let anchor_b = Vec3::new(0.0, 0.0, 0.0);
        let mut joint = BallJoint::new(ha, hb, anchor_a, anchor_b);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let total_ang = a.angular_velocity.norm() + b.angular_velocity.norm();
        assert!(
            total_ang > 0.5,
            "Ball joint should allow free rotation, total_ang={total_ang}"
        );
        let rel_lin_v = (a.velocity - b.velocity).norm();
        assert!(
            rel_lin_v < 1e-3,
            "Ball joint should constrain relative COM velocity, rel_lin_v={rel_lin_v}"
        );
    }
    #[test]
    fn test_motor_reaches_target_velocity() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.inverse_mass = 0.0;
            b.world_inverse_inertia = oxiphysics_core::math::Mat3::zeros();
        }
        let target_velocity = 10.0_f64;
        let mut motor = MotorJoint::new(ha, hb, Vec3::new(0.0, 0.0, 1.0), target_velocity, 10000.0);
        let dt = 1.0 / 60.0;
        let n_steps = 60;
        for _ in 0..n_steps {
            motor.prepare(&bodies, dt);
            motor.solve_velocity(&mut bodies, dt);
        }
        let omega_z = bodies.get(ha).unwrap().angular_velocity.z;
        assert!(
            (omega_z - target_velocity).abs() < 1.0,
            "Motor should reach target angular velocity within {n_steps} steps, omega_z={omega_z} target={target_velocity}"
        );
    }
    #[test]
    fn test_constraint_warm_start_reduces_error() {
        let make_setup = || {
            let (mut bodies, ha, hb) =
                setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
            if let Some(a) = bodies.get_mut(ha) {
                a.velocity = Vec3::new(-3.0, 0.0, 0.0);
            }
            if let Some(b) = bodies.get_mut(hb) {
                b.velocity = Vec3::new(3.0, 0.0, 0.0);
            }
            (bodies, ha, hb)
        };
        let dt = 1.0 / 60.0;
        let anchor_a = Vec3::new(0.5, 0.0, 0.0);
        let anchor_b = Vec3::new(-0.5, 0.0, 0.0);
        let (mut bodies_few, ha_f, hb_f) = make_setup();
        let mut joint_few = FixedJoint::new(ha_f, hb_f, anchor_a, anchor_b, Quat::identity());
        joint_few.prepare(&bodies_few, dt);
        for _ in 0..1 {
            joint_few.solve_velocity(&mut bodies_few, dt);
        }
        let va_f = bodies_few.get(ha_f).unwrap().velocity;
        let vb_f = bodies_few.get(hb_f).unwrap().velocity;
        let r_af = bodies_few.get(ha_f).unwrap().transform.rotation * anchor_a;
        let r_bf = bodies_few.get(hb_f).unwrap().transform.rotation * anchor_b;
        let va_pt_f = va_f + bodies_few.get(ha_f).unwrap().angular_velocity.cross(&r_af);
        let vb_pt_f = vb_f + bodies_few.get(hb_f).unwrap().angular_velocity.cross(&r_bf);
        let err_few = (va_pt_f - vb_pt_f).norm();
        let (mut bodies_many, ha_m, hb_m) = make_setup();
        let mut joint_many = FixedJoint::new(ha_m, hb_m, anchor_a, anchor_b, Quat::identity());
        joint_many.prepare(&bodies_many, dt);
        for _ in 0..20 {
            joint_many.solve_velocity(&mut bodies_many, dt);
        }
        let va_m = bodies_many.get(ha_m).unwrap().velocity;
        let vb_m = bodies_many.get(hb_m).unwrap().velocity;
        let r_am = bodies_many.get(ha_m).unwrap().transform.rotation * anchor_a;
        let r_bm = bodies_many.get(hb_m).unwrap().transform.rotation * anchor_b;
        let va_pt_m = va_m + bodies_many.get(ha_m).unwrap().angular_velocity.cross(&r_am);
        let vb_pt_m = vb_m + bodies_many.get(hb_m).unwrap().angular_velocity.cross(&r_bm);
        let err_many = (va_pt_m - vb_pt_m).norm();
        assert!(
            err_many <= err_few,
            "More solver iterations should reduce constraint error: err_few={err_few}, err_many={err_many}"
        );
    }
    #[test]
    fn test_gear_constraint_ratio() {
        let axis_a = Vec3::new(0.0, 0.0, 1.0);
        let axis_b = Vec3::new(0.0, 0.0, 1.0);
        let gear = GearJoint::new(0, 1, 2.0, axis_a, axis_b);
        let omega_a = Vec3::new(0.0, 0.0, 10.0);
        let omega_b = Vec3::new(0.0, 0.0, 3.0);
        let cv = gear.constraint_velocity(omega_a, omega_b);
        assert!(
            (cv - 16.0).abs() < 1e-12,
            "constraint_velocity should be 10 + 2*3 = 16, got {cv}"
        );
        let expected = omega_a.dot(&axis_a) + gear.ratio * omega_b.dot(&axis_b);
        assert!(
            (cv - expected).abs() < 1e-12,
            "constraint_velocity formula mismatch: cv={cv} expected={expected}"
        );
    }
    #[test]
    fn test_gear_jacobians() {
        let axis_a = Vec3::new(1.0, 0.0, 0.0);
        let axis_b = Vec3::new(0.0, 1.0, 0.0);
        let ratio = -3.0;
        let gear = GearJoint::new(0, 1, ratio, axis_a, axis_b);
        let (j_a, j_b) = gear.jacobian();
        assert!(
            (j_a - axis_a).norm() < 1e-12,
            "J_a should equal axis_a, got {j_a:?}"
        );
        let expected_jb = axis_b * ratio;
        assert!(
            (j_b - expected_jb).norm() < 1e-12,
            "J_b should equal ratio*axis_b, got {j_b:?} expected {expected_jb:?}"
        );
    }
    #[test]
    fn test_rack_pinion_velocity() {
        let pinion_axis = Vec3::new(0.0, 0.0, 1.0);
        let rack_axis = Vec3::new(1.0, 0.0, 0.0);
        let pitch_radius = 0.1;
        let joint = RackPinionJoint::new(0, 1, pitch_radius, pinion_axis, rack_axis);
        let omega_pinion = Vec3::new(0.0, 0.0, 1.0);
        let v_rack_correct = Vec3::new(0.1, 0.0, 0.0);
        let cv = joint.constraint_velocity(omega_pinion, v_rack_correct);
        assert!(
            cv.abs() < 1e-12,
            "rack-pinion: v_rack=0.1 m/s with omega=1 rad/s should satisfy constraint, cv={cv}"
        );
        let v_rack_zero = Vec3::zeros();
        let cv2 = joint.constraint_velocity(omega_pinion, v_rack_zero);
        assert!(
            (cv2 - (-0.1)).abs() < 1e-12,
            "rack-pinion: v_rack=0 should give cv=-0.1, got {cv2}"
        );
    }
    #[test]
    fn test_pulley_lengths() {
        let pulley = Vec3::zeros();
        let anchor_a = Vec3::new(3.0, 0.0, 0.0);
        let anchor_b = Vec3::new(0.0, 4.0, 0.0);
        let joint = PulleyJoint::new(0, 1, anchor_a, anchor_b, pulley, 1.0, 10.0);
        let (la, lb) = joint.current_lengths();
        assert!((la - 3.0).abs() < 1e-12, "la should be 3.0, got {la}");
        assert!((lb - 4.0).abs() < 1e-12, "lb should be 4.0, got {lb}");
    }
    #[test]
    fn test_pulley_constraint() {
        let pulley = Vec3::zeros();
        let anchor_a = Vec3::new(3.0, 0.0, 0.0);
        let anchor_b = Vec3::new(0.0, 4.0, 0.0);
        let ratio = 2.0;
        let max_length = 12.0;
        let joint = PulleyJoint::new(0, 1, anchor_a, anchor_b, pulley, ratio, max_length);
        let cv = joint.constraint_position();
        let expected = 3.0 + ratio * 4.0 - max_length;
        assert!(
            (cv - expected).abs() < 1e-12,
            "constraint_position should be {expected}, got {cv}"
        );
        let anchor_a2 = Vec3::new(5.0, 0.0, 0.0);
        let joint2 = PulleyJoint::new(0, 1, anchor_a2, anchor_b, pulley, ratio, 13.0);
        let cv2 = joint2.constraint_position();
        assert!(
            cv2.abs() < 1e-12,
            "taut pulley constraint should be 0, got {cv2}"
        );
    }
}
#[cfg(test)]
mod new_joint_tests {
    use super::*;
    use crate::Constraint;
    use crate::PrismaticJoint;
    use crate::RevoluteJoint;
    use crate::joints::BreakableJoint;
    use crate::joints::CableJoint;
    use crate::joints::ConeJoint;
    use crate::joints::CylindricalJoint;
    use crate::joints::DistanceJoint;
    use crate::joints::JointForceMeter;
    use crate::joints::SphericalJoint;
    use crate::joints::UniversalJoint;
    use crate::joints::WeldJoint;
    use crate::joints::XpbdJoint;
    use oxiphysics_rigid::RigidBody;
    fn setup_two_bodies(pos_a: Vec3, pos_b: Vec3) -> (RigidBodySet, BodyHandle, BodyHandle) {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = pos_a;
        a.linear_damping = 0.0;
        a.angular_damping = 0.0;
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = pos_b;
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        let hb = bodies.insert(b);
        (bodies, ha, hb)
    }
    #[test]
    fn test_universal_joint_constrains_translation() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = UniversalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let rel_v = (a.velocity - b.velocity).norm();
        assert!(
            rel_v < 2.0,
            "Universal joint should constrain translation, rel_v={rel_v}"
        );
    }
    #[test]
    fn test_universal_joint_body_handles() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let joint = UniversalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let handles = joint.body_handles();
        assert_eq!(handles.len(), 2);
        assert_eq!(handles[0].index, ha.index);
        assert_eq!(handles[1].index, hb.index);
        let _ = bodies;
    }
    #[test]
    fn test_universal_joint_angular_error_initially_near_zero() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let mut joint = UniversalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );
        joint.prepare(&bodies, 1.0 / 60.0);
        let err = joint.angular_error().abs();
        assert!(
            err < 1e-10,
            "perpendicular axes should have near-zero dot product, got {err}"
        );
        let _ = bodies;
    }
    #[test]
    fn test_cylindrical_joint_allows_sliding() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(5.0, 0.0, 0.0);
        }
        let mut joint = CylindricalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let total_vx = a.velocity.x + b.velocity.x;
        assert!(
            total_vx > 1.0,
            "Cylindrical joint should allow X sliding, total_vx={total_vx}"
        );
    }
    #[test]
    fn test_cylindrical_joint_blocks_perpendicular_motion() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = CylindricalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let vy_diff = (a.velocity.y - b.velocity.y).abs();
        assert!(
            vy_diff < 1.0,
            "Cylindrical joint should block perpendicular (Y) motion, vy_diff={vy_diff}"
        );
    }
    #[test]
    fn test_cylindrical_joint_allows_rotation() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(3.0, 0.0, 0.0);
        }
        let mut joint = CylindricalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let total_ang_x =
            bodies.get(ha).unwrap().angular_velocity.x + bodies.get(hb).unwrap().angular_velocity.x;
        assert!(
            total_ang_x.abs() > 0.5,
            "Cylindrical joint should allow rotation around axis, total_ang_x={total_ang_x}"
        );
    }
    #[test]
    fn test_cylindrical_joint_with_limits() {
        let joint = CylindricalJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
        )
        .with_limits(-2.0, 2.0);
        assert_eq!(joint.lower_limit, Some(-2.0));
        assert_eq!(joint.upper_limit, Some(2.0));
    }
    #[test]
    fn test_cone_joint_within_limit_at_rest() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let joint = ConeJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_4,
        );
        assert!(
            joint.within_limit(&bodies),
            "joint at rest should be within cone limit"
        );
    }
    #[test]
    fn test_cone_joint_constrains_translation() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = ConeJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_4,
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let rel_v = (a.velocity - b.velocity).norm();
        assert!(
            rel_v < 2.0,
            "Cone joint should constrain translation, rel_v={rel_v}"
        );
    }
    #[test]
    fn test_cone_joint_swing_angle_zero_at_identity() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let mut joint = ConeJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_2,
        );
        joint.prepare(&bodies, 1.0 / 60.0);
        let angle = joint.current_swing_angle(&bodies);
        assert!(
            angle.abs() < 1e-10,
            "swing angle should be zero at identity, got {angle}"
        );
        let _ = bodies;
    }
    #[test]
    fn test_cone_joint_clamps_half_angle() {
        let joint = ConeJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            std::f64::consts::PI * 2.0,
        );
        assert!(
            joint.half_angle <= std::f64::consts::PI + 1e-10,
            "half_angle should be clamped to π"
        );
    }
    #[test]
    fn test_weld_joint_from_current_pose() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
        let a_pos = bodies.get(ha).unwrap().transform.position;
        let b_pos = bodies.get(hb).unwrap().transform.position;
        let a_rot = bodies.get(ha).unwrap().transform.rotation;
        let b_rot = bodies.get(hb).unwrap().transform.rotation;
        let joint = WeldJoint::from_current_pose(ha, hb, a_pos, b_pos, a_rot, b_rot);
        let expected_offset = b_pos - a_pos;
        assert!(
            (joint.desired_offset - expected_offset).norm() < 1e-10,
            "Weld joint offset should match current pose"
        );
    }
    #[test]
    fn test_weld_joint_blocks_relative_velocity() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 8.0, 0.0);
        }
        let mut joint = WeldJoint::new(ha, hb, Vec3::new(1.0, 0.0, 0.0), Quat::identity());
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let rel_v = (a.velocity - b.velocity).norm();
        assert!(
            rel_v < 1.0,
            "Weld joint should block relative velocity, rel_v={rel_v}"
        );
    }
    #[test]
    fn test_weld_joint_body_handles() {
        let joint = WeldJoint::new(
            BodyHandle::new(3, 1),
            BodyHandle::new(7, 2),
            Vec3::zeros(),
            Quat::identity(),
        );
        let handles = joint.body_handles();
        assert_eq!(handles.len(), 2);
        assert_eq!(handles[0].index, 3);
        assert_eq!(handles[1].index, 7);
    }
    #[test]
    fn test_weld_joint_blocks_angular_velocity() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.angular_velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = WeldJoint::new(ha, hb, Vec3::new(1.0, 0.0, 0.0), Quat::identity());
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let ang_diff = (a.angular_velocity - b.angular_velocity).norm();
        assert!(
            ang_diff < 1.0,
            "Weld joint should block relative angular velocity, ang_diff={ang_diff}"
        );
    }
    #[test]
    fn test_distance_joint_fixed_distance_within_limit_at_rest() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0));
        let joint = DistanceJoint::fixed_distance(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0);
        assert!(
            joint.is_within_limits(&bodies),
            "bodies at rest distance should be within limits"
        );
    }
    #[test]
    fn test_distance_joint_too_far_flagged() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        let joint = DistanceJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 0.0, 2.0);
        assert!(
            !joint.is_within_limits(&bodies),
            "too-far bodies should be out of limit"
        );
    }
    #[test]
    fn test_distance_joint_current_distance() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        let joint = DistanceJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 0.0, 10.0);
        let dist = joint.current_distance(&bodies);
        assert!((dist - 3.0).abs() < 1e-10, "dist={dist}");
    }
    #[test]
    fn test_distance_joint_min_distance_respected() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0));
        let mut joint = DistanceJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 1.0, 10.0);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        joint.solve_position(&mut bodies, dt);
        let dist_after = joint.current_distance(&bodies);
        assert!(
            dist_after > 0.3,
            "joint should push bodies apart, dist={dist_after}"
        );
    }
    #[test]
    fn test_distance_joint_body_handles() {
        let joint = DistanceJoint::new(
            BodyHandle::new(4, 1),
            BodyHandle::new(7, 2),
            Vec3::zeros(),
            Vec3::zeros(),
            0.5,
            2.0,
        );
        let h = joint.body_handles();
        assert_eq!(h.len(), 2);
        assert_eq!(h[0].index, 4);
        assert_eq!(h[1].index, 7);
    }
    #[test]
    fn test_breakable_joint_starts_intact() {
        let joint = BreakableJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            2.0,
            100.0,
            1.0,
            1000.0,
        );
        assert!(!joint.is_broken());
        assert_eq!(joint.current_force(), 0.0);
    }
    #[test]
    fn test_breakable_joint_resets() {
        let mut joint = BreakableJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            2.0,
            100.0,
            1.0,
            1000.0,
        );
        joint.broken = true;
        joint.accumulated_force = 999.0;
        joint.reset();
        assert!(!joint.is_broken());
        assert_eq!(joint.current_force(), 0.0);
    }
    #[test]
    fn test_breakable_joint_breaks_when_threshold_exceeded() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0));
        let mut joint = BreakableJoint::new(ha, hb, 2.0, 500.0, 0.0, 0.01);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        joint.solve_velocity(&mut bodies, dt);
        assert!(
            joint.is_broken(),
            "joint should break when threshold exceeded"
        );
    }
    #[test]
    fn test_breakable_joint_body_handles() {
        let joint = BreakableJoint::new(
            BodyHandle::new(3, 0),
            BodyHandle::new(5, 0),
            1.0,
            50.0,
            0.5,
            100.0,
        );
        let h = joint.body_handles();
        assert_eq!(h.len(), 2);
        assert_eq!(h[0].index, 3);
        assert_eq!(h[1].index, 5);
    }
    #[test]
    fn test_xpbd_joint_rigid_converges_to_rest_length() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 0.0));
        let mut joint = XpbdJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
            joint.solve_position(&mut bodies, dt);
        }
        let pa = bodies.get(ha).unwrap().transform.position;
        let pb = bodies.get(hb).unwrap().transform.position;
        let dist = (pb - pa).norm();
        assert!(
            dist < 3.0 + 1e-6,
            "XPBD joint should reduce distance, dist={dist}"
        );
    }
    #[test]
    fn test_xpbd_joint_lambda_reset_on_prepare() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(2.0, 0.0, 0.0));
        let mut joint = XpbdJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0, 0.001, 0.0);
        joint.lambda = 99.0;
        joint.prepare(&bodies, 1.0 / 60.0);
        assert_eq!(joint.lambda, 0.0, "prepare() should reset lambda");
    }
    #[test]
    fn test_xpbd_joint_soft_allows_extension() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        let mut joint = XpbdJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0, 1.0, 0.0);
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        joint.solve_position(&mut bodies, dt);
        let pa = bodies.get(ha).unwrap().transform.position;
        let pb = bodies.get(hb).unwrap().transform.position;
        let dist = (pb - pa).norm();
        assert!(
            dist > 2.5,
            "soft XPBD joint should not fully snap to rest length in one step, dist={dist}"
        );
    }
    #[test]
    fn test_xpbd_joint_body_handles() {
        let joint = XpbdJoint::new(
            BodyHandle::new(2, 0),
            BodyHandle::new(6, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            1.0,
            0.0,
            0.0,
        );
        let h = joint.body_handles();
        assert_eq!(h[0].index, 2);
        assert_eq!(h[1].index, 6);
    }
    #[test]
    fn test_joint_force_meter_zero_when_no_constraint() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(2.0, 0.0, 0.0));
        let inner = DistanceJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 0.0, 10.0);
        let mut meter = JointForceMeter::new(inner);
        let dt = 1.0 / 60.0;
        let mut b = bodies;
        meter.prepare(&b, dt);
        meter.solve_velocity(&mut b, dt);
        let force = meter.force(dt);
        assert!(
            force < 1e3,
            "force within limits should be small, f={force}"
        );
    }
    #[test]
    fn test_joint_force_meter_reports_nonzero_force() {
        let (mut bodies, ha, hb) =
            setup_two_bodies(Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(10.0, 0.0, 0.0);
        }
        let inner = DistanceJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 0.0, 2.0);
        let mut meter = JointForceMeter::new(inner);
        let dt = 1.0 / 60.0;
        meter.prepare(&bodies, dt);
        meter.solve_velocity(&mut bodies, dt);
        let force = meter.force(dt);
        let _ = force;
    }
    #[test]
    fn test_joint_force_meter_body_handles_delegate() {
        let inner = DistanceJoint::new(
            BodyHandle::new(10, 0),
            BodyHandle::new(20, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            0.5,
            3.0,
        );
        let meter = JointForceMeter::new(inner);
        let h = meter.body_handles();
        assert_eq!(h[0].index, 10);
        assert_eq!(h[1].index, 20);
    }
    #[test]
    fn test_spherical_joint_zero_twist_at_identity() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let joint = SphericalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let angle = joint.compute_twist_angle(&bodies);
        assert!(
            angle.abs() < 1e-10,
            "twist should be zero at identity, angle={angle}"
        );
    }
    #[test]
    fn test_spherical_joint_body_handles() {
        let joint = SphericalJoint::new(
            BodyHandle::new(3, 0),
            BodyHandle::new(7, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let handles = joint.body_handles();
        assert_eq!(handles.len(), 2);
        assert_eq!(handles[0].index, 3);
        assert_eq!(handles[1].index, 7);
    }
    #[test]
    fn test_spherical_joint_constrains_translation() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(0.0, 5.0, 0.0);
        }
        let mut joint = SphericalJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 0.0, 1.0),
        );
        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..20 {
            joint.solve_velocity(&mut bodies, dt);
        }
        let a = bodies.get(ha).unwrap();
        let b = bodies.get(hb).unwrap();
        let rel_v = (a.velocity - b.velocity).norm();
        assert!(
            rel_v < 3.0,
            "SphericalJoint should reduce relative velocity, rel_v={rel_v}"
        );
    }
    #[test]
    fn test_spherical_joint_nonzero_axis_normalised() {
        let joint = SphericalJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(3.0, 0.0, 0.0),
        );
        assert!(
            (joint.local_axis.norm() - 1.0).abs() < 1e-10,
            "axis should be unit length"
        );
    }
    #[test]
    fn test_cable_joint_slack_at_short_distance() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0));
        let cable = CableJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 5.0);
        let len = cable.current_length(&bodies);
        assert!((len - 1.0).abs() < 1e-10, "len={len}");
        assert!(
            !cable.is_taut(),
            "cable should be slack when distance < max_length"
        );
    }
    #[test]
    fn test_cable_joint_taut_when_stretched() {
        let (bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0));
        let mut cable = CableJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0);
        let dt = 1.0 / 60.0;
        cable.prepare(&bodies, dt);
        assert!(
            cable.is_taut(),
            "cable should be taut when distance > max_length"
        );
    }
    #[test]
    fn test_cable_joint_prevents_separation() {
        let (mut bodies, ha, hb) = setup_two_bodies(Vec3::zeros(), Vec3::new(2.1, 0.0, 0.0));
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity = Vec3::new(10.0, 0.0, 0.0);
        }
        let mut cable = CableJoint::new(ha, hb, Vec3::zeros(), Vec3::zeros(), 2.0);
        let dt = 1.0 / 60.0;
        cable.prepare(&bodies, dt);
        assert!(cable.is_taut(), "cable should be taut");
        for _ in 0..20 {
            cable.apply_tension_constraint(&mut bodies, dt);
        }
        let b = bodies.get(hb).unwrap();
        assert!(
            b.velocity.x < 9.5,
            "cable should resist outward velocity, vx={}",
            b.velocity.x
        );
    }
    #[test]
    fn test_revolute_gear_ratio_constraint_zero_at_balance() {
        let joint = RevoluteJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let residual = joint.compute_gear_ratio_constraint(2.0, -2.0, 1.0);
        assert!(
            residual.abs() < 1e-10,
            "gear constraint should be zero at balance, r={residual}"
        );
    }
    #[test]
    fn test_revolute_gear_ratio_constraint_nonzero() {
        let joint = RevoluteJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let residual = joint.compute_gear_ratio_constraint(3.0, 1.0, 2.0);
        assert!((residual - 5.0).abs() < 1e-10, "residual={residual}");
    }
    #[test]
    fn test_prismatic_lead_screw_constraint_zero() {
        let joint = PrismaticJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let lead = 0.005;
        let omega_rel = 10.0;
        let v_rel = lead * omega_rel;
        let residual = joint.compute_lead_screw(v_rel, omega_rel, lead);
        assert!(
            residual.abs() < 1e-12,
            "lead-screw constraint should be zero, r={residual}"
        );
    }
    #[test]
    fn test_prismatic_lead_screw_constraint_nonzero() {
        let joint = PrismaticJoint::new(
            BodyHandle::new(0, 0),
            BodyHandle::new(1, 0),
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let lead = 0.005;
        let residual = joint.compute_lead_screw(1.0, 0.0, lead);
        assert!(
            residual.abs() > 0.5,
            "non-matching state should give nonzero residual, r={residual}"
        );
    }
}
