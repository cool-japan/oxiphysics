//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::Transform;
use oxiphysics_core::math::{Quat, Unit, Vec3};

use super::types::{
    AbaLinkData, ArticulatedChain2D, ArticulatedLink2D, ExternalForce, JointDof, RneaLink,
    SpatialInertia, SpatialVec,
};

/// Build the incremental `Transform` produced by a joint at coordinate `q`.
pub(super) fn joint_displacement(joint: &JointDof, q: f64) -> Transform {
    match joint {
        JointDof::Fixed | JointDof::Spherical | JointDof::Planar => Transform::default(),
        JointDof::Revolute(axis) => {
            let ax = Vec3::new(axis[0], axis[1], axis[2]);
            let norm = ax.norm();
            if norm < 1e-12 {
                return Transform::default();
            }
            let unit_ax = Unit::new_normalize(ax);
            let rotation = Quat::from_axis_angle(&unit_ax, q);
            Transform::new(Vec3::zeros(), rotation)
        }
        JointDof::Prismatic(axis) => {
            let ax = Vec3::new(axis[0], axis[1], axis[2]);
            let norm = ax.norm();
            if norm < 1e-12 {
                return Transform::default();
            }
            let unit_ax = ax / norm;
            Transform::new(unit_ax * q, Quat::identity())
        }
    }
}
/// Period of a simple pendulum: T = 2π√(L/g).
pub fn simple_pendulum_period(length: f64, g: f64) -> f64 {
    use std::f64::consts::TAU;
    TAU * (length / g).sqrt()
}
/// One RK4 step of a double pendulum.
///
/// Returns `(theta1_new, omega1_new, theta2_new, omega2_new)`.
pub fn double_pendulum_step(
    theta1: f64,
    omega1: f64,
    theta2: f64,
    omega2: f64,
    m1: f64,
    m2: f64,
    l1: f64,
    l2: f64,
    g: f64,
    dt: f64,
) -> (f64, f64, f64, f64) {
    let derivs = |th1: f64, om1: f64, th2: f64, om2: f64| -> (f64, f64, f64, f64) {
        let delta = th2 - th1;
        let sd = delta.sin();
        let cd = delta.cos();
        let denom1 = (m1 + m2) * l1 - m2 * l1 * cd * cd;
        let denom2 = (l2 / l1) * denom1;
        let alpha1 =
            (m2 * l1 * om1 * om1 * sd * cd + m2 * g * th2.sin() * cd + m2 * l2 * om2 * om2 * sd
                - (m1 + m2) * g * th1.sin())
                / denom1;
        let alpha2 = (-(m1 + m2) * l1 * om1 * om1 * sd * cd / l2 - (m1 + m2) * g * th2.sin() / l2
            + (m1 + m2) * g * th1.sin() * cd / l2
            - m2 * om2 * om2 * sd * cd)
            / (denom2 / l2);
        let _ = alpha2;
        let alpha2_correct = (-m2 * l2 * om2 * om2 * sd * cd - (m1 + m2) * g * th2.sin()
            + (m1 + m2) * g * th1.sin() * cd
            - (m1 + m2) * l1 * om1 * om1 * sd)
            / denom2;
        (om1, alpha1, om2, alpha2_correct)
    };
    let (k1_t1, k1_o1, k1_t2, k1_o2) = derivs(theta1, omega1, theta2, omega2);
    let (k2_t1, k2_o1, k2_t2, k2_o2) = derivs(
        theta1 + 0.5 * dt * k1_t1,
        omega1 + 0.5 * dt * k1_o1,
        theta2 + 0.5 * dt * k1_t2,
        omega2 + 0.5 * dt * k1_o2,
    );
    let (k3_t1, k3_o1, k3_t2, k3_o2) = derivs(
        theta1 + 0.5 * dt * k2_t1,
        omega1 + 0.5 * dt * k2_o1,
        theta2 + 0.5 * dt * k2_t2,
        omega2 + 0.5 * dt * k2_o2,
    );
    let (k4_t1, k4_o1, k4_t2, k4_o2) = derivs(
        theta1 + dt * k3_t1,
        omega1 + dt * k3_o1,
        theta2 + dt * k3_t2,
        omega2 + dt * k3_o2,
    );
    let theta1_new = theta1 + dt / 6.0 * (k1_t1 + 2.0 * k2_t1 + 2.0 * k3_t1 + k4_t1);
    let omega1_new = omega1 + dt / 6.0 * (k1_o1 + 2.0 * k2_o1 + 2.0 * k3_o1 + k4_o1);
    let theta2_new = theta2 + dt / 6.0 * (k1_t2 + 2.0 * k2_t2 + 2.0 * k3_t2 + k4_t2);
    let omega2_new = omega2 + dt / 6.0 * (k1_o2 + 2.0 * k2_o2 + 2.0 * k3_o2 + k4_o2);
    (theta1_new, omega1_new, theta2_new, omega2_new)
}
/// Return the world-space position of the last link's tip in a 2-D chain.
pub fn end_effector_position(
    chain: &ArticulatedChain2D,
    base_pos: [f64; 2],
    base_angle: f64,
) -> [f64; 2] {
    let fk = chain.forward_kinematics(base_pos, base_angle);
    if fk.is_empty() {
        base_pos
    } else {
        fk[fk.len() - 1].0
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::articulated::ArticulatedChain;

    use crate::articulated::JointType;
    use crate::articulated::LinkId;

    use std::f64::consts::FRAC_PI_2;
    fn zero_inertia() -> [[f64; 3]; 3] {
        [[0.0; 3]; 3]
    }
    fn build_revolute_arm(len: f64) -> (ArticulatedChain, LinkId, LinkId) {
        let mut chain = ArticulatedChain::new();
        let root = chain.add_root(1.0, zero_inertia());
        let child = chain.add_link(
            root,
            JointDof::Revolute([0.0, 0.0, 1.0]),
            Transform::from_position(Vec3::new(0.0, len, 0.0)),
            1.0,
            zero_inertia(),
        );
        (chain, root, child)
    }
    #[test]
    fn test_chain_two_links() {
        let (chain, _root, child) = build_revolute_arm(1.0);
        let fk = chain.compute_forward_kinematics();
        let child_tf = &fk[child.0 as usize].1;
        assert!(
            (child_tf.position.x).abs() < 1e-10,
            "x={}",
            child_tf.position.x
        );
        assert!(
            (child_tf.position.y - 1.0).abs() < 1e-10,
            "y={}",
            child_tf.position.y
        );
        assert!(
            (child_tf.position.z).abs() < 1e-10,
            "z={}",
            child_tf.position.z
        );
    }
    #[test]
    fn test_chain_revolute_90deg() {
        let (chain, _root, _child) = build_revolute_arm(1.0);
        let mut chain2 = ArticulatedChain::new();
        let root2 = chain2.add_root(1.0, zero_inertia());
        let _ = chain;
        chain2.set_joint_position(root2, FRAC_PI_2);
        let child2 = chain2.add_link(
            root2,
            JointDof::Revolute([0.0, 0.0, 1.0]),
            Transform::from_position(Vec3::new(0.0, 1.0, 0.0)),
            1.0,
            zero_inertia(),
        );
        let mut chain3 = ArticulatedChain::new();
        let root3 = chain3.add_root(1.0, zero_inertia());
        chain3.links[root3.0 as usize].joint_type = JointDof::Revolute([0.0, 0.0, 1.0]);
        chain3.set_joint_position(root3, FRAC_PI_2);
        let child3 = chain3.add_link(
            root3,
            JointDof::Fixed,
            Transform::from_position(Vec3::new(0.0, 1.0, 0.0)),
            1.0,
            zero_inertia(),
        );
        let fk3 = chain3.compute_forward_kinematics();
        let pos3 = fk3[child3.0 as usize].1.position;
        assert!(
            (pos3.x + 1.0).abs() < 1e-10,
            "expected x=-1, got x={}",
            pos3.x
        );
        assert!(pos3.y.abs() < 1e-10, "expected y=0, got y={}", pos3.y);
        assert!(pos3.z.abs() < 1e-10, "expected z=0, got z={}", pos3.z);
        let _ = child2;
    }
    #[test]
    fn test_chain_prismatic() {
        let mut chain = ArticulatedChain::new();
        let root = chain.add_root(1.0, zero_inertia());
        let child = chain.add_link(
            root,
            JointDof::Prismatic([0.0, 1.0, 0.0]),
            Transform::default(),
            1.0,
            zero_inertia(),
        );
        chain.set_joint_position(child, 2.5);
        let fk = chain.compute_forward_kinematics();
        let pos = fk[child.0 as usize].1.position;
        assert!(pos.x.abs() < 1e-10, "x={}", pos.x);
        assert!((pos.y - 2.5).abs() < 1e-10, "y={}", pos.y);
        assert!(pos.z.abs() < 1e-10, "z={}", pos.z);
    }
    #[test]
    fn test_chain_integrate() {
        let mut chain = ArticulatedChain::new();
        let root = chain.add_root(1.0, zero_inertia());
        let child = chain.add_link(
            root,
            JointDof::Revolute([0.0, 0.0, 1.0]),
            Transform::default(),
            1.0,
            zero_inertia(),
        );
        chain.set_joint_velocity(child, 1.0);
        chain.integrate(0.1);
        let q = chain.links[child.0 as usize].joint_position;
        assert!((q - 0.1).abs() < 1e-10, "q={}", q);
    }
    #[test]
    fn test_chain_joint_limits() {
        let mut chain = ArticulatedChain::new();
        let root = chain.add_root(1.0, zero_inertia());
        let child = chain.add_link(
            root,
            JointDof::Revolute([0.0, 0.0, 1.0]),
            Transform::default(),
            1.0,
            zero_inertia(),
        );
        chain.links[child.0 as usize].joint_limit = Some((0.0, 0.5));
        chain.set_joint_position(child, 0.45);
        chain.set_joint_velocity(child, 10.0);
        chain.apply_velocity_limits();
        chain.integrate(0.1);
        let q = chain.links[child.0 as usize].joint_position;
        assert!(q <= 0.5 + 1e-10, "q={} should be <= 0.5", q);
    }
    #[test]
    fn test_jacobian_single_revolute() {
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let joint_pos = Vec3::zeros();
        let ee_pos = Vec3::new(0.0, 1.0, 0.0);
        let col = ArticulatedChain::revolute_jacobian_column(axis, joint_pos, ee_pos);
        assert!((col[0] + 1.0).abs() < 1e-10, "linear.x={}", col[0]);
        assert!(col[1].abs() < 1e-10, "linear.y={}", col[1]);
        assert!(col[2].abs() < 1e-10, "linear.z={}", col[2]);
        assert!(col[3].abs() < 1e-10, "angular.x={}", col[3]);
        assert!(col[4].abs() < 1e-10, "angular.y={}", col[4]);
        assert!((col[5] - 1.0).abs() < 1e-10, "angular.z={}", col[5]);
    }
    #[test]
    fn test_2d_chain_new_empty() {
        let chain = ArticulatedChain2D::new();
        assert_eq!(chain.link_count(), 0);
    }
    #[test]
    fn test_2d_chain_add_link() {
        let mut chain = ArticulatedChain2D::new();
        let idx = chain.add_link(1.0, 1.0, None);
        assert_eq!(idx, 0);
        assert_eq!(chain.link_count(), 1);
    }
    #[test]
    fn test_2d_chain_add_two_links() {
        let mut chain = ArticulatedChain2D::new();
        let a = chain.add_link(1.0, 1.0, None);
        let b = chain.add_link(1.0, 1.0, Some(a));
        assert_eq!(chain.link_count(), 2);
        assert_eq!(chain.links[b].parent, Some(a));
    }
    #[test]
    fn test_2d_fk_single_link_zero_angle() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 2.0, None);
        let fk = chain.forward_kinematics([0.0, 0.0], 0.0);
        let (tip, angle) = fk[0];
        assert!((tip[0] - 2.0).abs() < 1e-10, "tip.x={}", tip[0]);
        assert!(tip[1].abs() < 1e-10, "tip.y={}", tip[1]);
        assert!(angle.abs() < 1e-10);
    }
    #[test]
    fn test_2d_fk_single_link_90deg() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_angle = FRAC_PI_2;
        let fk = chain.forward_kinematics([0.0, 0.0], 0.0);
        let (tip, _) = fk[0];
        assert!(tip[0].abs() < 1e-10, "tip.x={}", tip[0]);
        assert!((tip[1] - 1.0).abs() < 1e-10, "tip.y={}", tip[1]);
    }
    #[test]
    fn test_2d_fk_two_links_straight() {
        let mut chain = ArticulatedChain2D::new();
        let a = chain.add_link(1.0, 1.0, None);
        chain.add_link(1.0, 1.0, Some(a));
        let fk = chain.forward_kinematics([0.0, 0.0], 0.0);
        let (tip, _) = fk[1];
        assert!((tip[0] - 2.0).abs() < 1e-10, "tip.x={}", tip[0]);
        assert!(tip[1].abs() < 1e-10, "tip.y={}", tip[1]);
    }
    #[test]
    fn test_2d_fk_base_offset() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let fk = chain.forward_kinematics([3.0, 4.0], 0.0);
        let (tip, _) = fk[0];
        assert!((tip[0] - 4.0).abs() < 1e-10, "tip.x={}", tip[0]);
        assert!((tip[1] - 4.0).abs() < 1e-10, "tip.y={}", tip[1]);
    }
    #[test]
    fn test_2d_end_effector_empty_chain() {
        let chain = ArticulatedChain2D::new();
        let ee = end_effector_position(&chain, [1.0, 2.0], 0.0);
        assert_eq!(ee, [1.0, 2.0]);
    }
    #[test]
    fn test_2d_end_effector_single_link() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 3.0, None);
        let ee = end_effector_position(&chain, [0.0, 0.0], 0.0);
        assert!((ee[0] - 3.0).abs() < 1e-10);
        assert!(ee[1].abs() < 1e-10);
    }
    #[test]
    fn test_2d_jacobian_single_link() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let jac = chain.jacobian([0.0, 0.0], 0.0);
        assert_eq!(jac.len(), 1);
        assert!(jac[0][0].abs() < 1e-10, "dx/dθ={}", jac[0][0]);
        assert!((jac[0][1] - 1.0).abs() < 1e-10, "dy/dθ={}", jac[0][1]);
    }
    #[test]
    fn test_2d_jacobian_empty_chain() {
        let chain = ArticulatedChain2D::new();
        let jac = chain.jacobian([0.0, 0.0], 0.0);
        assert!(jac.is_empty());
    }
    #[test]
    fn test_2d_apply_torques_basic() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let inertia = 1.0 / 3.0;
        let torque = 1.0;
        chain.apply_torques(&[torque], 1.0);
        let expected_omega = torque / inertia;
        let expected_theta = expected_omega;
        assert!((chain.links[0].joint_velocity - expected_omega).abs() < 1e-10);
        assert!((chain.links[0].joint_angle - expected_theta).abs() < 1e-10);
    }
    #[test]
    fn test_2d_kinetic_energy_zero() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(2.0, 1.0, None);
        assert_eq!(chain.total_kinetic_energy(), 0.0);
    }
    #[test]
    fn test_2d_kinetic_energy_nonzero() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_velocity = 2.0;
        let ke = chain.total_kinetic_energy();
        let expected = 0.5 * (1.0 / 3.0) * 4.0;
        assert!((ke - expected).abs() < 1e-10, "ke={}", ke);
    }
    #[test]
    fn test_2d_potential_energy_horizontal() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let pe = chain.potential_energy([0.0, 0.0], 0.0, 9.81);
        assert!(pe.abs() < 1e-10, "pe={}", pe);
    }
    #[test]
    fn test_2d_potential_energy_vertical() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_angle = FRAC_PI_2;
        let pe = chain.potential_energy([0.0, 0.0], 0.0, 9.81);
        assert!((pe - 9.81 * 0.5).abs() < 1e-10, "pe={}", pe);
    }
    #[test]
    fn test_2d_inverse_dynamics_no_gravity_no_accel() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let torques = chain.inverse_dynamics(&[0.0], [0.0, 0.0]);
        assert_eq!(torques.len(), 1);
        assert!(torques[0].abs() < 1e-10);
    }
    #[test]
    fn test_2d_inverse_dynamics_inertia_only() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let torques = chain.inverse_dynamics(&[3.0], [0.0, 0.0]);
        assert_eq!(torques.len(), 1);
        let expected = (1.0 / 3.0) * 3.0;
        assert!((torques[0] - expected).abs() < 1e-10, "τ={}", torques[0]);
    }
    #[test]
    fn test_2d_default_chain() {
        let chain = ArticulatedChain2D::default();
        assert_eq!(chain.link_count(), 0);
    }
    #[test]
    fn test_simple_pendulum_period() {
        let t = simple_pendulum_period(1.0, 9.81);
        let expected = 2.0 * std::f64::consts::PI * (1.0_f64 / 9.81).sqrt();
        assert!((t - expected).abs() < 1e-10, "T={}", t);
    }
    #[test]
    fn test_simple_pendulum_longer_arm() {
        let t1 = simple_pendulum_period(1.0, 9.81);
        let t4 = simple_pendulum_period(4.0, 9.81);
        assert!((t4 - 2.0 * t1).abs() < 1e-10);
    }
    #[test]
    fn test_double_pendulum_step_stationary() {
        let (t1, o1, t2, o2) =
            double_pendulum_step(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 9.81, 0.001);
        let _ = (t1, o1, t2, o2);
    }
    #[test]
    fn test_double_pendulum_step_small_angle() {
        let (t1, o1, t2, o2) =
            double_pendulum_step(0.1, 0.0, 0.1, 0.0, 1.0, 1.0, 1.0, 1.0, 9.81, 0.001);
        assert!(t1.is_finite());
        assert!(o1.is_finite());
        assert!(t2.is_finite());
        assert!(o2.is_finite());
    }
    #[test]
    fn test_double_pendulum_energy_conservation_approx() {
        let theta1_0 = 0.1_f64;
        let theta2_0 = 0.1_f64;
        let m1 = 1.0_f64;
        let m2 = 1.0_f64;
        let l1 = 1.0_f64;
        let l2 = 1.0_f64;
        let g = 9.81_f64;
        let pe0 = -(m1 + m2) * g * l1 * theta1_0.cos() - m2 * g * l2 * theta2_0.cos();
        let mut th1 = theta1_0;
        let mut om1 = 0.0_f64;
        let mut th2 = theta2_0;
        let mut om2 = 0.0_f64;
        for _ in 0..100 {
            let (a, b, c, d) = double_pendulum_step(th1, om1, th2, om2, m1, m2, l1, l2, g, 0.001);
            th1 = a;
            om1 = b;
            th2 = c;
            om2 = d;
        }
        let pe1 = -(m1 + m2) * g * l1 * th1.cos() - m2 * g * l2 * th2.cos();
        let ke1 = 0.5 * m1 * l1 * l1 * om1 * om1
            + 0.5
                * m2
                * (l1 * l1 * om1 * om1
                    + l2 * l2 * om2 * om2
                    + 2.0 * l1 * l2 * om1 * om2 * (th2 - th1).cos());
        let e0 = pe0;
        let e1 = pe1 + ke1;
        let rel_err = (e1 - e0).abs() / e0.abs().max(1e-10);
        assert!(
            rel_err < 0.05,
            "energy not conserved: e0={e0}, e1={e1}, rel={rel_err}"
        );
    }
    #[test]
    fn test_2d_link_body_fields() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(3.0, 2.0, None);
        let link = &chain.links[0];
        assert_eq!(link.body.mass, 3.0);
        assert_eq!(link.body.length, 2.0);
        assert!((link.body.inertia - 4.0).abs() < 1e-10);
        assert_eq!(link.joint_type, JointType::Revolute);
    }
    #[test]
    fn test_2d_joint_types() {
        assert_ne!(JointType::Revolute, JointType::Prismatic);
        assert_ne!(JointType::Revolute, JointType::Fixed);
        assert_ne!(JointType::Prismatic, JointType::Fixed);
        assert_eq!(JointType::Fixed, JointType::Fixed);
    }
}
/// Compute the inverse dynamics (joint torques) for a single revolute joint
/// chain using the Recursive Newton-Euler Algorithm.
///
/// This is a simplified 1-DOF (revolute-z) version for `n` links.
///
/// Returns `Vec`f64` of joint torques.
pub fn rnea_revolute_chain(
    n: usize,
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    qdd: &[f64],
    gravity: [f64; 3],
) -> Vec<f64> {
    assert!(n > 0);
    let mut links: Vec<RneaLink> = vec![RneaLink::default(); n];
    let a_base = SpatialVec::new([0.0; 3], [-gravity[0], -gravity[1], -gravity[2]]);
    for i in 0..n {
        let q_i = q[i];
        let qd_i = qd[i];
        let qdd_i = qdd[i];
        let s = SpatialVec::new([0.0, 0.0, 1.0], [0.0; 3]);
        links[i].velocity = if i == 0 {
            s.scale(qd_i)
        } else {
            links[i - 1].velocity.add(&s.scale(qd_i))
        };
        let parent_acc = if i == 0 {
            a_base
        } else {
            links[i - 1].acceleration
        };
        let v_i = links[i].velocity;
        let coriolis = v_i.cross_motion(&s.scale(qd_i));
        links[i].acceleration = parent_acc.add(&s.scale(qdd_i)).add(&coriolis);
        let m = masses.get(i).copied().unwrap_or(1.0);
        let l = lengths.get(i).copied().unwrap_or(1.0);
        let com = [l / 2.0 * q_i.cos(), l / 2.0 * q_i.sin(), 0.0];
        let i_zz = m * l * l / 3.0;
        let inertia = SpatialInertia::new(m, com, [[0.0; 3], [0.0; 3], [0.0, 0.0, i_zz]]);
        links[i].force = inertia.multiply(&links[i].acceleration);
    }
    let mut torques = vec![0.0_f64; n];
    let s = SpatialVec::new([0.0, 0.0, 1.0], [0.0; 3]);
    for i in (0..n).rev() {
        let f = links[i].force;
        torques[i] = f.dot(&s);
    }
    torques
}
/// Compute the composite rigid-body inertia for a single link by combining
/// the link's own inertia with those of its children.
///
/// This is a scalar (1-D) approximation: `I_A_i = I_i + Σ_c I_A_c`.
pub fn composite_rigid_body_inertia(i: usize, inertias: &[f64], children: &[Vec<usize>]) -> f64 {
    let mut total = inertias.get(i).copied().unwrap_or(0.0);
    if let Some(kids) = children.get(i) {
        for &c in kids {
            total += composite_rigid_body_inertia(c, inertias, children);
        }
    }
    total
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn mat3_mul_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
#[cfg(test)]
mod tests_featherstone {

    use crate::articulated::JointImpedanceController;

    use crate::articulated::LoopClosureConstraint;
    use crate::articulated::PlanarJoint;
    use crate::articulated::ScrewJoint;
    use crate::articulated::SpatialInertia;
    use crate::articulated::SpatialVec;

    use crate::articulated::composite_rigid_body_inertia;

    use crate::articulated::rnea_revolute_chain;

    #[test]
    fn test_spatial_vec_dot() {
        let a = SpatialVec::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let b = SpatialVec::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        let d = a.dot(&b);
        assert!((d - 2.0).abs() < 1e-12, "dot = {d}");
    }
    #[test]
    fn test_spatial_vec_scale() {
        let v = SpatialVec::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let v2 = v.scale(2.0);
        assert!((v2.angular[0] - 2.0).abs() < 1e-12);
        assert!((v2.linear[2] - 12.0).abs() < 1e-12);
    }
    #[test]
    fn test_spatial_vec_add() {
        let a = SpatialVec::new([1.0; 3], [2.0; 3]);
        let b = SpatialVec::new([3.0; 3], [4.0; 3]);
        let c = a.add(&b);
        assert!((c.angular[0] - 4.0).abs() < 1e-12);
        assert!((c.linear[0] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_spatial_inertia_sphere() {
        let si = SpatialInertia::sphere(2.0, 0.5, [0.0; 3]);
        let expected = 2.0 / 5.0 * 2.0 * 0.25;
        assert!((si.i_rot[0][0] - expected).abs() < 1e-10);
    }
    #[test]
    fn test_rnea_single_link_gravity_torque() {
        let taus =
            rnea_revolute_chain(1, &[1.0], &[1.0], &[0.0], &[0.0], &[0.0], [0.0, -9.81, 0.0]);
        assert_eq!(taus.len(), 1);
        assert!(taus[0].is_finite(), "torque should be finite: {}", taus[0]);
    }
    #[test]
    fn test_composite_rigid_body_inertia_leaf() {
        let inertias = vec![5.0, 3.0, 2.0];
        let children: Vec<Vec<usize>> = vec![vec![1, 2], vec![], vec![]];
        let total = composite_rigid_body_inertia(0, &inertias, &children);
        assert!((total - 10.0).abs() < 1e-10, "total = {total}");
    }
    #[test]
    fn test_loop_closure_satisfied() {
        let mut c = LoopClosureConstraint::new(0, 1);
        c.violation = [0.0001, 0.0, 0.0];
        assert!(c.is_satisfied(0.001));
        c.violation = [0.01, 0.0, 0.0];
        assert!(!c.is_satisfied(0.001));
    }
    #[test]
    fn test_screw_joint_advance() {
        let mut sj = ScrewJoint::new(0.01);
        sj.advance(2.0 * std::f64::consts::PI);
        assert!(
            (sj.displacement - 0.01).abs() < 1e-10,
            "disp = {}",
            sj.displacement
        );
    }
    #[test]
    fn test_screw_joint_lead() {
        let sj = ScrewJoint::new(0.005);
        let lead = sj.lead();
        assert!((lead - 0.005 / (2.0 * std::f64::consts::PI)).abs() < 1e-12);
    }
    #[test]
    fn test_planar_joint_integrate() {
        let mut pj = PlanarJoint::new();
        pj.integrate(1.0, 2.0, 0.0, 1.0);
        assert!((pj.tx - 1.0).abs() < 1e-12);
        assert!((pj.ty - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_planar_joint_homogeneous_identity_at_zero() {
        let pj = PlanarJoint::new();
        let h = pj.to_homogeneous();
        assert!((h[0][0] - 1.0).abs() < 1e-12);
        assert!((h[1][0]).abs() < 1e-12);
        assert!((h[0][2]).abs() < 1e-12);
    }
    #[test]
    fn test_impedance_controller_torque() {
        let ctrl = JointImpedanceController::new(100.0, 10.0);
        let tau = ctrl.torque(0.1, 0.05);
        assert!((tau - 10.5).abs() < 1e-10, "tau = {tau}");
    }
    #[test]
    fn test_impedance_controller_torques_vector() {
        let ctrl = JointImpedanceController::new(100.0, 10.0);
        let taus = ctrl.torques(&[0.1, 0.2], &[0.0, 0.0]);
        assert!((taus[0] - 10.0).abs() < 1e-10);
        assert!((taus[1] - 20.0).abs() < 1e-10);
    }
}
/// Scalar (1-DOF per link) Articulated Body Algorithm for a revolute chain.
///
/// Given joint torques and the current state (q, q_dot), returns
/// joint accelerations q_ddot via the ABA.
///
/// Assumptions:
/// - All joints are revolute around the z-axis.
/// - The joint motion subspace is scalar `s = 1`.
/// - Inertias are scalar (moment about z).
///
/// # Arguments
/// * `n`       – number of links
/// * `inertias`– moment of inertia of each link about its joint axis (kg·m²)
/// * `masses`  – mass of each link (kg)
/// * `lengths` – link length (m), used for gravity torque
/// * `q`       – joint positions (rad)
/// * `qd`      – joint velocities (rad/s)
/// * `torques` – applied joint torques (N·m)
/// * `gravity` – gravitational acceleration vector [gx, gy, gz] (m/s²)
///
/// Returns `Vec`f64` of joint accelerations (rad/s²).
pub fn aba_revolute_chain(
    n: usize,
    inertias: &[f64],
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    torques: &[f64],
    gravity: [f64; 3],
) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut data: Vec<AbaLinkData> = vec![AbaLinkData::default(); n];
    let mut v_cumulative = 0.0_f64;
    for i in 0..n {
        let qd_i = qd.get(i).copied().unwrap_or(0.0);
        v_cumulative += qd_i;
        data[i].v = v_cumulative;
        let i_i = inertias.get(i).copied().unwrap_or(1.0);
        data[i].abi = i_i;
        let m_i = masses.get(i).copied().unwrap_or(1.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let abs_angle: f64 = q[..=i].iter().sum();
        let g_perp = -gravity[0] * abs_angle.sin() + gravity[1] * abs_angle.cos();
        data[i].bias = -m_i * g_perp * (l_i / 2.0);
        let _ = q_i;
    }
    for i in (0..(n - 1)).rev() {
        let child = i + 1;
        data[i].abi += data[child].abi;
        data[i].bias += data[child].bias;
    }
    let mut a_parent = 0.0_f64;
    let mut result = vec![0.0_f64; n];
    for i in 0..n {
        let tau_i = torques.get(i).copied().unwrap_or(0.0);
        let abi_i = data[i].abi.max(1e-30);
        let qdd_i = (tau_i - data[i].bias - a_parent) / abi_i;
        result[i] = qdd_i;
        a_parent += qdd_i;
    }
    result
}
/// Assemble the joint-space inertia matrix H (n×n) for a revolute chain.
///
/// Uses the Composite Rigid Body (CRB) algorithm:
/// H_ij = s_i · I_C_i · s_j  where s is the joint motion axis and I_C is
/// the composite inertia.
///
/// For a simple scalar revolute chain (all joints about z):
/// H_ij = Σ_{k = max(i,j)}^{n-1}  I_k  (contribution of links distal to both joints).
///
/// Returns the n×n symmetric matrix as a flat `Vec`f64` in row-major order.
pub fn joint_space_inertia_matrix(n: usize, inertias: &[f64]) -> Vec<f64> {
    let mut h = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let start = i.max(j);
            let val: f64 = (start..n)
                .map(|k| inertias.get(k).copied().unwrap_or(0.0))
                .sum();
            h[i * n + j] = val;
        }
    }
    h
}
/// Compute gravity compensation torques for a revolute chain.
///
/// The gravity compensation torque at joint i is:
/// τ_g_i = Σ_{j=i}^{n-1}  m_j * g_perp_j * r_{ij}
///
/// where `g_perp_j` is the component of gravity perpendicular to the
/// link j axis, and `r_{ij}` is the distance from joint i to the COM of link j.
pub fn gravity_compensation_torques(
    n: usize,
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    gravity: [f64; 3],
) -> Vec<f64> {
    let mut torques = vec![0.0_f64; n];
    for i in 0..n {
        let mut tau = 0.0;
        for j in i..n {
            let m_j = masses.get(j).copied().unwrap_or(1.0);
            let l_j = lengths.get(j).copied().unwrap_or(1.0);
            let abs_angle: f64 = q[..=j].iter().sum();
            let g_perp = -gravity[0] * abs_angle.sin() + gravity[1] * abs_angle.cos();
            let r: f64 = if j > i {
                lengths[i..j].iter().sum::<f64>() + l_j / 2.0
            } else {
                l_j / 2.0
            };
            tau -= m_j * g_perp * r;
        }
        torques[i] = tau;
    }
    torques
}
/// Compute forward dynamics (joint accelerations) with applied external forces.
///
/// External forces are projected onto the joint axis and added to the
/// applied torques before running the ABA.
pub fn forward_dynamics_with_external(
    n: usize,
    inertias: &[f64],
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    torques: &[f64],
    gravity: [f64; 3],
    external: &[ExternalForce],
) -> Vec<f64> {
    let mut tau_combined: Vec<f64> = torques
        .iter()
        .copied()
        .chain(std::iter::repeat(0.0))
        .take(n)
        .collect();
    for ext in external {
        if ext.link_idx < n {
            tau_combined[ext.link_idx] += ext.torque;
        }
    }
    aba_revolute_chain(n, inertias, masses, lengths, q, qd, &tau_combined, gravity)
}
/// Compute the 2×n analytical Jacobian of the end-effector position for a
/// planar revolute chain.
///
/// The end-effector is at the tip of the last link.
///
/// Returns `n` columns `\[∂x/∂θ_i, ∂y/∂θ_i\]`.
pub fn end_effector_jacobian_2d(
    n: usize,
    lengths: &[f64],
    q: &[f64],
    base: [f64; 2],
) -> Vec<[f64; 2]> {
    if n == 0 {
        return Vec::new();
    }
    let mut origins = vec![[0.0_f64; 2]; n + 1];
    origins[0] = base;
    let mut cum_angle = 0.0_f64;
    for i in 0..n {
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        cum_angle += q_i;
        origins[i + 1] = [
            origins[i][0] + l_i * cum_angle.cos(),
            origins[i][1] + l_i * cum_angle.sin(),
        ];
    }
    let ee = origins[n];
    (0..n)
        .map(|i| {
            let dx = ee[0] - origins[i][0];
            let dy = ee[1] - origins[i][1];
            [-dy, dx]
        })
        .collect()
}
/// Compute the 6×n spatial Jacobian for a planar revolute chain (all joints
/// revolute about z, links in the XY plane).
///
/// Each column is `\[0, 0, 1, -z_{ei} * 0, y_{ori→ee}, -x_{ori→ee}\]`
/// simplified to `\[linear(3), angular(3)\]`.
///
/// Returns `n` columns of `\[f64; 6\]`: `\[v_x, v_y, v_z, ω_x, ω_y, ω_z\]`.
pub fn spatial_jacobian_revolute_z(
    n: usize,
    lengths: &[f64],
    q: &[f64],
    base: [f64; 2],
) -> Vec<[f64; 6]> {
    if n == 0 {
        return Vec::new();
    }
    let mut origins = vec![[0.0_f64; 2]; n + 1];
    origins[0] = base;
    let mut cum_angle = 0.0_f64;
    for i in 0..n {
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        cum_angle += q_i;
        origins[i + 1] = [
            origins[i][0] + l_i * cum_angle.cos(),
            origins[i][1] + l_i * cum_angle.sin(),
        ];
    }
    let ee = origins[n];
    (0..n)
        .map(|i| {
            let dx = ee[0] - origins[i][0];
            let dy = ee[1] - origins[i][1];
            [-dy, dx, 0.0, 0.0, 0.0, 1.0]
        })
        .collect()
}
/// Damped least-squares (DLS) inverse kinematics step for a 2-D chain.
///
/// Computes Δq = J† Δx using  J† = Jᵀ (J Jᵀ + λ² I)⁻¹.
///
/// For a 2×n Jacobian acting on the 2-D position error,
/// returns the joint angle updates Δq.
///
/// * `jacobian` – n columns `\[∂x/∂θ_i, ∂y/∂θ_i\]`
/// * `dx`       – position error [Δx, Δy]
/// * `lambda`   – damping factor (prevents singularity)
pub fn dls_ik_step_2d(jacobian: &[[f64; 2]], dx: [f64; 2], lambda: f64) -> Vec<f64> {
    let n = jacobian.len();
    if n == 0 {
        return Vec::new();
    }
    let mut jjt = [[0.0_f64; 2]; 2];
    for col in jacobian {
        jjt[0][0] += col[0] * col[0];
        jjt[0][1] += col[0] * col[1];
        jjt[1][0] += col[1] * col[0];
        jjt[1][1] += col[1] * col[1];
    }
    let lam2 = lambda * lambda;
    jjt[0][0] += lam2;
    jjt[1][1] += lam2;
    let det = jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0];
    if det.abs() < 1e-30 {
        return vec![0.0; n];
    }
    let inv = [
        [jjt[1][1] / det, -jjt[0][1] / det],
        [-jjt[1][0] / det, jjt[0][0] / det],
    ];
    let alpha = [
        inv[0][0] * dx[0] + inv[0][1] * dx[1],
        inv[1][0] * dx[0] + inv[1][1] * dx[1],
    ];
    jacobian
        .iter()
        .map(|col| col[0] * alpha[0] + col[1] * alpha[1])
        .collect()
}
/// Compute the inverse dynamics for a planar revolute chain using the
/// Newton-Euler recursive algorithm in 3-D (revolute about z).
///
/// Returns joint torques.
pub fn rnea_3d_revolute_z(
    n: usize,
    masses: &[f64],
    inertias_zz: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    qdd: &[f64],
    gravity: [f64; 3],
) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut omega = vec![[0.0_f64; 3]; n];
    let mut alpha = vec![[0.0_f64; 3]; n];
    let mut a_com = vec![[0.0_f64; 3]; n];
    let a_base = [-gravity[0], -gravity[1], -gravity[2]];
    let mut cum_angle = 0.0_f64;
    for i in 0..n {
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let qd_i = qd.get(i).copied().unwrap_or(0.0);
        let qdd_i = qdd.get(i).copied().unwrap_or(0.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        cum_angle += q_i;
        let parent_omega = if i == 0 { [0.0; 3] } else { omega[i - 1] };
        let parent_alpha = if i == 0 { [0.0; 3] } else { alpha[i - 1] };
        let parent_a = if i == 0 { a_base } else { a_com[i - 1] };
        omega[i] = [parent_omega[0], parent_omega[1], parent_omega[2] + qd_i];
        let w_cross_qd_z = [
            omega[i][1] * qd_i - omega[i][2] * 0.0,
            omega[i][2] * 0.0 - omega[i][0] * qd_i,
            0.0,
        ];
        alpha[i] = [
            parent_alpha[0] + w_cross_qd_z[0],
            parent_alpha[1] + w_cross_qd_z[1],
            parent_alpha[2] + qdd_i,
        ];
        let r_com = [
            (l_i / 2.0) * cum_angle.cos(),
            (l_i / 2.0) * cum_angle.sin(),
            0.0,
        ];
        let ax_r = cross3w(alpha[i], r_com);
        let w_r = cross3w(omega[i], r_com);
        let w_w_r = cross3w(omega[i], w_r);
        a_com[i] = add3w(add3w(parent_a, ax_r), w_w_r);
    }
    let mut f = vec![[0.0_f64; 3]; n];
    let mut tau_net = vec![[0.0_f64; 3]; n];
    let mut torques = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let m_i = masses.get(i).copied().unwrap_or(1.0);
        let i_zz = inertias_zz.get(i).copied().unwrap_or(1.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        let abs_angle: f64 = q[..=i].iter().sum();
        let r_com = [
            (l_i / 2.0) * abs_angle.cos(),
            (l_i / 2.0) * abs_angle.sin(),
            0.0,
        ];
        let fi = scale3w(a_com[i], m_i);
        let fi_child = if i + 1 < n { f[i + 1] } else { [0.0; 3] };
        f[i] = add3w(fi, fi_child);
        let i_alpha = [0.0, 0.0, i_zz * alpha[i][2]];
        let w_times_iw = cross3w(omega[i], [0.0, 0.0, i_zz * omega[i][2]]);
        let r_cross_f = cross3w(r_com, fi);
        let tau_child = if i + 1 < n { tau_net[i + 1] } else { [0.0; 3] };
        tau_net[i] = [
            i_alpha[0] + w_times_iw[0] + r_cross_f[0] + tau_child[0],
            i_alpha[1] + w_times_iw[1] + r_cross_f[1] + tau_child[1],
            i_alpha[2] + w_times_iw[2] + r_cross_f[2] + tau_child[2],
        ];
        torques[i] = tau_net[i][2];
    }
    torques
}
#[inline]
pub(super) fn add3w(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3w(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn cross3w(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests_aba {

    use crate::articulated::ExternalForce;

    use crate::articulated::Wrench;
    use crate::articulated::aba_revolute_chain;

    use crate::articulated::dls_ik_step_2d;
    use crate::articulated::end_effector_jacobian_2d;
    use crate::articulated::forward_dynamics_with_external;
    use crate::articulated::gravity_compensation_torques;
    use crate::articulated::joint_space_inertia_matrix;
    use crate::articulated::rnea_3d_revolute_z;

    use crate::articulated::spatial_jacobian_revolute_z;
    #[test]
    fn test_aba_single_link_no_gravity() {
        let qdd = aba_revolute_chain(
            1,
            &[1.0],
            &[1.0],
            &[1.0],
            &[0.0],
            &[0.0],
            &[1.0],
            [0.0, 0.0, 0.0],
        );
        assert_eq!(qdd.len(), 1);
        assert!(qdd[0].is_finite(), "qdd must be finite: {}", qdd[0]);
        assert!(
            qdd[0] > 0.0,
            "acceleration should be positive for positive torque"
        );
    }
    #[test]
    fn test_aba_zero_torque_with_gravity_compensation() {
        let m = 1.0;
        let l = 1.0;
        let g_mag = 9.81;
        let q = [0.0_f64];
        let gravity = [0.0, -g_mag, 0.0];
        let gc = gravity_compensation_torques(1, &[m], &[l], &q, gravity);
        let qdd = aba_revolute_chain(
            1,
            &[m * l * l / 3.0],
            &[m],
            &[l],
            &q,
            &[0.0],
            &[gc[0]],
            gravity,
        );
        assert!(
            qdd[0].abs() < 1e-9,
            "qdd should be ~0 with gravity compensation: {}",
            qdd[0]
        );
    }
    #[test]
    fn test_aba_empty_chain() {
        let result = aba_revolute_chain(0, &[], &[], &[], &[], &[], &[], [0.0; 3]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_aba_two_links_finite() {
        let qdd = aba_revolute_chain(
            2,
            &[0.5, 0.2],
            &[1.0, 0.5],
            &[1.0, 0.8],
            &[0.1, -0.2],
            &[0.0; 2],
            &[0.0; 2],
            [0.0, -9.81, 0.0],
        );
        assert_eq!(qdd.len(), 2);
        assert!(qdd[0].is_finite());
        assert!(qdd[1].is_finite());
    }
    #[test]
    fn test_jsim_single_link() {
        let h = joint_space_inertia_matrix(1, &[5.0]);
        assert_eq!(h.len(), 1);
        assert!((h[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_jsim_two_links() {
        let h = joint_space_inertia_matrix(2, &[3.0, 2.0]);
        assert_eq!(h.len(), 4);
        assert!((h[0] - 5.0).abs() < 1e-10, "H[0][0]={}", h[0]);
        assert!((h[1] - 2.0).abs() < 1e-10, "H[0][1]={}", h[1]);
        assert!((h[2] - 2.0).abs() < 1e-10, "H[1][0]={}", h[2]);
        assert!((h[3] - 2.0).abs() < 1e-10, "H[1][1]={}", h[3]);
    }
    #[test]
    fn test_jsim_symmetric() {
        let h = joint_space_inertia_matrix(3, &[1.0, 2.0, 3.0]);
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (h[i * n + j] - h[j * n + i]).abs() < 1e-10,
                    "H not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_gravity_compensation_single_horizontal() {
        let gc = gravity_compensation_torques(1, &[1.0], &[1.0], &[0.0], [0.0, -9.81, 0.0]);
        assert_eq!(gc.len(), 1);
        assert!(
            gc[0].abs() > 0.0,
            "gravity torque should be nonzero: {}",
            gc[0]
        );
        assert!(gc[0].is_finite());
    }
    #[test]
    fn test_gravity_compensation_two_links() {
        let gc = gravity_compensation_torques(
            2,
            &[1.0, 1.0],
            &[1.0, 1.0],
            &[0.0, 0.0],
            [0.0, -9.81, 0.0],
        );
        assert_eq!(gc.len(), 2);
        assert!(
            gc[0].abs() > gc[1].abs(),
            "joint 0 should carry more gravity torque"
        );
    }
    #[test]
    fn test_gravity_compensation_zero_gravity() {
        let gc =
            gravity_compensation_torques(2, &[1.0, 1.0], &[1.0, 1.0], &[0.0, 0.0], [0.0, 0.0, 0.0]);
        for &t in &gc {
            assert_eq!(t, 0.0, "no gravity → no compensation torque");
        }
    }
    #[test]
    fn test_forward_dynamics_external_force_adds_to_torque() {
        let qdd_no_ext =
            aba_revolute_chain(1, &[1.0], &[1.0], &[1.0], &[0.0], &[0.0], &[0.0], [0.0; 3]);
        let ext = vec![ExternalForce {
            link_idx: 0,
            torque: 1.0,
        }];
        let qdd_with_ext = forward_dynamics_with_external(
            1,
            &[1.0],
            &[1.0],
            &[1.0],
            &[0.0],
            &[0.0],
            &[0.0],
            [0.0; 3],
            &ext,
        );
        assert!(
            qdd_with_ext[0] > qdd_no_ext[0],
            "external force should increase acceleration"
        );
    }
    #[test]
    fn test_forward_dynamics_external_invalid_idx() {
        let ext = vec![ExternalForce {
            link_idx: 99,
            torque: 100.0,
        }];
        let qdd = forward_dynamics_with_external(
            1,
            &[1.0],
            &[1.0],
            &[1.0],
            &[0.0],
            &[0.0],
            &[0.0],
            [0.0; 3],
            &ext,
        );
        assert_eq!(qdd.len(), 1);
        assert!(qdd[0].is_finite());
    }
    #[test]
    fn test_ee_jacobian_single_link_at_zero() {
        let jac = end_effector_jacobian_2d(1, &[1.0], &[0.0], [0.0, 0.0]);
        assert_eq!(jac.len(), 1);
        assert!(jac[0][0].abs() < 1e-10, "dx/dθ={}", jac[0][0]);
        assert!((jac[0][1] - 1.0).abs() < 1e-10, "dy/dθ={}", jac[0][1]);
    }
    #[test]
    fn test_ee_jacobian_two_links_straight() {
        let jac = end_effector_jacobian_2d(2, &[1.0, 1.0], &[0.0, 0.0], [0.0, 0.0]);
        assert_eq!(jac.len(), 2);
        assert!(jac[0][0].abs() < 1e-10);
        assert!((jac[0][1] - 2.0).abs() < 1e-10);
        assert!(jac[1][0].abs() < 1e-10);
        assert!((jac[1][1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ee_jacobian_empty() {
        let jac = end_effector_jacobian_2d(0, &[], &[], [0.0, 0.0]);
        assert!(jac.is_empty());
    }
    #[test]
    fn test_spatial_jacobian_single_link() {
        let jac = spatial_jacobian_revolute_z(1, &[1.0], &[0.0], [0.0, 0.0]);
        assert_eq!(jac.len(), 1);
        assert!(jac[0][0].abs() < 1e-10);
        assert!((jac[0][1] - 1.0).abs() < 1e-10);
        assert_eq!(jac[0][2], 0.0);
        assert_eq!(jac[0][3], 0.0);
        assert_eq!(jac[0][4], 0.0);
        assert_eq!(jac[0][5], 1.0);
    }
    #[test]
    fn test_dls_ik_single_joint_small_error() {
        let jac = vec![[0.0_f64, 1.0]];
        let dq = dls_ik_step_2d(&jac, [0.0, 0.1], 0.01);
        assert_eq!(dq.len(), 1);
        assert!(
            dq[0] > 0.0,
            "dq should be positive for positive y error: {}",
            dq[0]
        );
    }
    #[test]
    fn test_dls_ik_empty_jacobian() {
        let dq = dls_ik_step_2d(&[], [1.0, 0.0], 0.01);
        assert!(dq.is_empty());
    }
    #[test]
    fn test_dls_ik_reduces_error() {
        let jac = end_effector_jacobian_2d(2, &[1.0, 1.0], &[0.0, 0.0], [0.0, 0.0]);
        let dx = [0.0, 0.5];
        let dq = dls_ik_step_2d(&jac, dx, 0.01);
        assert_eq!(dq.len(), 2);
        for d in &dq {
            assert!(d.is_finite(), "dq must be finite: {d}");
        }
        let total: f64 = dq.iter().map(|x| x.abs()).sum();
        assert!(total > 0.0, "joint changes should be nonzero");
    }
    #[test]
    fn test_rnea_3d_single_link_gravity() {
        let taus = rnea_3d_revolute_z(
            1,
            &[1.0],
            &[1.0 / 3.0],
            &[1.0],
            &[0.0],
            &[0.0],
            &[0.0],
            [0.0, -9.81, 0.0],
        );
        assert_eq!(taus.len(), 1);
        assert!(taus[0].is_finite(), "torque must be finite: {}", taus[0]);
    }
    #[test]
    fn test_rnea_3d_no_gravity_no_motion() {
        let taus = rnea_3d_revolute_z(
            2,
            &[1.0, 0.5],
            &[0.33, 0.1],
            &[1.0, 0.8],
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        assert_eq!(taus.len(), 2);
        for &t in &taus {
            assert!(
                t.abs() < 1e-10,
                "zero gravity+motion → zero torques, got {t}"
            );
        }
    }
    #[test]
    fn test_rnea_3d_empty_chain() {
        let taus = rnea_3d_revolute_z(0, &[], &[], &[], &[], &[], &[], [0.0; 3]);
        assert!(taus.is_empty());
    }
    #[test]
    fn test_wrench_add() {
        let w1 = Wrench {
            force: [1.0, 2.0, 3.0],
            torque: [4.0, 5.0, 6.0],
        };
        let w2 = Wrench {
            force: [1.0, 1.0, 1.0],
            torque: [1.0, 1.0, 1.0],
        };
        let w = w1.add(&w2);
        assert!((w.force[0] - 2.0).abs() < 1e-10);
        assert!((w.torque[2] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_wrench_scale() {
        let w = Wrench {
            force: [2.0, 4.0, 6.0],
            torque: [1.0, 0.0, 0.0],
        };
        let ws = w.scale(0.5);
        assert!((ws.force[0] - 1.0).abs() < 1e-10);
        assert!((ws.force[1] - 2.0).abs() < 1e-10);
    }
}
/// Compute the full articulated body inertia (ABI) for a single-chain open
/// kinematic tree using the Featherstone ABA pass-1 backward sweep.
///
/// Returns the ABI at the root as a `Vec`f64` (flattened 6×6 matrix in
/// row-major order).  For a single-DOF revolute chain (z-axis) this reduces
/// to the scalar inertia at each joint.
///
/// # Arguments
/// * `n`       – number of links (indexed 0..n-1)
/// * `inertias` – scalar moment of inertia about joint axis for each link
/// * `masses`  – mass of each link
/// * `lengths` – link length
pub fn abi_backward_pass(n: usize, inertias: &[f64], masses: &[f64], lengths: &[f64]) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut abi = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let i_i = inertias.get(i).copied().unwrap_or(1.0);
        let m_i = masses.get(i).copied().unwrap_or(1.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        let abi_i = i_i + m_i * (l_i / 2.0) * (l_i / 2.0);
        abi[i] = abi_i + if i + 1 < n { abi[i + 1] } else { 0.0 };
    }
    abi
}
/// Compute the net generalised forces at each joint due to inertial and
/// velocity-product effects only (no external torques, no gravity).
///
/// Useful for isolating the Coriolis / centrifugal bias vector `C(q, q_dot)`.
///
/// Returns `Vec`f64` of bias forces (length `n`).
pub fn coriolis_bias_forces(
    n: usize,
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let mut omega_z = vec![0.0_f64; n];
    for i in 0..n {
        let qd_i = qd.get(i).copied().unwrap_or(0.0);
        omega_z[i] = if i == 0 { qd_i } else { omega_z[i - 1] + qd_i };
    }
    let mut bias = vec![0.0_f64; n];
    for i in 0..n {
        let mut tau = 0.0_f64;
        for j in i..n {
            let m_j = masses.get(j).copied().unwrap_or(1.0);
            let l_j = lengths.get(j).copied().unwrap_or(1.0);
            let abs_angle: f64 = q[..=j].iter().sum();
            let r_com = l_j / 2.0;
            let omega_sq = omega_z[j] * omega_z[j];
            let accel_perp = omega_sq * r_com;
            let dist_to_com: f64 = if j > i {
                lengths[i..j].iter().sum::<f64>() + l_j / 2.0
            } else {
                l_j / 2.0
            };
            let _ = abs_angle;
            tau += m_j * accel_perp * dist_to_com;
        }
        bias[i] = tau;
    }
    bias
}
/// Compute joint accelerations from torques using the joint-space equation:
///
/// H(q) * q_dd = τ - C(q, q_dot) - G(q)
///
/// Solves by computing H, C, G and inverting H (scalar chain → diagonal solve).
///
/// Returns joint accelerations `q_ddot`.
pub fn forward_dynamics_full(
    n: usize,
    inertias: &[f64],
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    torques: &[f64],
    gravity: [f64; 3],
) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let h = joint_space_inertia_matrix(n, inertias);
    let g_vec = gravity_compensation_torques(n, masses, lengths, q, gravity);
    let c_vec = coriolis_bias_forces(n, masses, lengths, q, qd);
    let diag_h: Vec<f64> = (0..n).map(|i| h[i * n + i]).collect();
    let mut qdd = vec![0.0_f64; n];
    for i in 0..n {
        let tau_i = torques.get(i).copied().unwrap_or(0.0);
        let tau_eff =
            tau_i - c_vec.get(i).copied().unwrap_or(0.0) - g_vec.get(i).copied().unwrap_or(0.0);
        qdd[i] = tau_eff / diag_h[i].max(1e-30);
    }
    qdd
}
/// Compute the total kinetic energy of a revolute chain.
///
/// T = 0.5 * q_dot^T H(q) q_dot
pub fn kinetic_energy(n: usize, inertias: &[f64], q: &[f64], qd: &[f64]) -> f64 {
    let h = joint_space_inertia_matrix(n, inertias);
    let _ = q;
    let mut ke = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let qd_i = qd.get(i).copied().unwrap_or(0.0);
            let qd_j = qd.get(j).copied().unwrap_or(0.0);
            ke += h[i * n + j] * qd_i * qd_j;
        }
    }
    0.5 * ke
}
/// Compute the total gravitational potential energy for a revolute planar chain.
///
/// PE = Σ_i m_i * g * h_com_i  where h_com_i is the y-coordinate of link i's CoM.
pub fn potential_energy_chain(
    n: usize,
    masses: &[f64],
    lengths: &[f64],
    q: &[f64],
    gravity_y: f64,
    base_pos: [f64; 2],
) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut cum_angle = 0.0_f64;
    let mut origins = vec![[0.0_f64; 2]; n + 1];
    origins[0] = base_pos;
    for i in 0..n {
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        cum_angle += q_i;
        origins[i + 1] = [
            origins[i][0] + l_i * cum_angle.cos(),
            origins[i][1] + l_i * cum_angle.sin(),
        ];
    }
    let mut pe = 0.0_f64;
    for i in 0..n {
        let m_i = masses.get(i).copied().unwrap_or(1.0);
        let com_y = (origins[i][1] + origins[i + 1][1]) / 2.0;
        pe += m_i * gravity_y * com_y;
    }
    pe
}
/// Compute per-link end-point linear velocities for a planar revolute chain.
///
/// v_tip_i = J_i * q_dot (2-D)
///
/// Returns `\[vx, vy\]` for each link's tip.
pub fn link_tip_velocities(
    n: usize,
    lengths: &[f64],
    q: &[f64],
    qd: &[f64],
    base_pos: [f64; 2],
) -> Vec<[f64; 2]> {
    if n == 0 {
        return Vec::new();
    }
    let mut cum_angle = 0.0_f64;
    let mut origins = vec![[0.0_f64; 2]; n + 1];
    origins[0] = base_pos;
    let mut angles = vec![0.0_f64; n];
    for i in 0..n {
        let q_i = q.get(i).copied().unwrap_or(0.0);
        let l_i = lengths.get(i).copied().unwrap_or(1.0);
        cum_angle += q_i;
        angles[i] = cum_angle;
        origins[i + 1] = [
            origins[i][0] + l_i * angles[i].cos(),
            origins[i][1] + l_i * angles[i].sin(),
        ];
    }
    let mut velocities = vec![[0.0_f64; 2]; n];
    for i in 0..n {
        let tip = origins[i + 1];
        let mut vx = 0.0_f64;
        let mut vy = 0.0_f64;
        for (j, orig_j) in origins.iter().enumerate().take(i + 1) {
            let dx = tip[0] - orig_j[0];
            let dy = tip[1] - orig_j[1];
            let qd_j = qd.get(j).copied().unwrap_or(0.0);
            vx += -dy * qd_j;
            vy += dx * qd_j;
        }
        velocities[i] = [vx, vy];
    }
    velocities
}
/// Apply hard joint-limit stops to a revolute chain.
///
/// For each link that has a joint limit and whose current position is outside
/// `\[q_min, q_max\]`, the joint is clamped and the velocity is set to zero
/// (inelastic collision with the stop).
///
/// Returns `true` if any limit was triggered.
pub fn enforce_joint_limits(links: &mut [ArticulatedLink2D], limits: &[(f64, f64)]) -> bool {
    let mut triggered = false;
    for (link, &(q_min, q_max)) in links.iter_mut().zip(limits.iter()) {
        if link.joint_angle < q_min {
            link.joint_angle = q_min;
            link.joint_velocity = link.joint_velocity.max(0.0);
            triggered = true;
        } else if link.joint_angle > q_max {
            link.joint_angle = q_max;
            link.joint_velocity = link.joint_velocity.min(0.0);
            triggered = true;
        }
    }
    triggered
}
/// Sample the reachable workspace of a planar revolute chain uniformly.
///
/// For `n_samples` random configurations, compute the end-effector position.
/// Returns a `Vec<\[f64; 2\]>` of sampled end-effector positions.
///
/// Uses a deterministic LCG for reproducible results.
pub fn sample_workspace(
    chain: &ArticulatedChain2D,
    base_pos: [f64; 2],
    n_samples: usize,
    q_mins: &[f64],
    q_maxs: &[f64],
) -> Vec<[f64; 2]> {
    let n = chain.links.len();
    if n == 0 || n_samples == 0 {
        return Vec::new();
    }
    let mut positions = Vec::with_capacity(n_samples);
    let mut lcg_state: u64 = 0x1234567890abcdef;
    let lcg_next = |s: &mut u64| -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*s >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..n_samples {
        let mut test_chain = chain.clone();
        for i in 0..n {
            let q_min = q_mins.get(i).copied().unwrap_or(-std::f64::consts::PI);
            let q_max = q_maxs.get(i).copied().unwrap_or(std::f64::consts::PI);
            let t = lcg_next(&mut lcg_state);
            test_chain.links[i].joint_angle = q_min + t * (q_max - q_min);
        }
        let fk = test_chain.forward_kinematics(base_pos, 0.0);
        if let Some(&(tip, _)) = fk.last() {
            positions.push(tip);
        }
    }
    positions
}
/// Cubic polynomial joint trajectory between two configurations.
///
/// Given start `q0` and end `q1` joint positions with zero initial and
/// final velocities, returns the position at normalised time `t ∈ \[0, 1\]`:
///
/// q(t) = q0 + (q1 - q0) * (3t² - 2t³)
pub fn cubic_spline_position(q0: f64, q1: f64, t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    q0 + (q1 - q0) * (3.0 * t * t - 2.0 * t * t * t)
}
/// Velocity of the cubic spline at normalised time `t`:
///
/// q_dot(t) = (q1 - q0) * 6t(1 - t) / T
///
/// where `T` is the total trajectory duration.
pub fn cubic_spline_velocity(q0: f64, q1: f64, t: f64, duration: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    (q1 - q0) * 6.0 * t * (1.0 - t) / duration.max(1e-12)
}
/// Sample a full multi-joint cubic spline trajectory.
///
/// Returns `n_steps` joint configurations `\[n_steps × n_joints\]` at evenly
/// spaced time points from 0 to `duration`.
pub fn sample_cubic_trajectory(
    q_start: &[f64],
    q_end: &[f64],
    _duration: f64,
    n_steps: usize,
) -> Vec<Vec<f64>> {
    assert_eq!(q_start.len(), q_end.len());
    let n_joints = q_start.len();
    (0..n_steps)
        .map(|step| {
            let t = if n_steps <= 1 {
                0.0
            } else {
                step as f64 / (n_steps - 1) as f64
            };
            (0..n_joints)
                .map(|j| cubic_spline_position(q_start[j], q_end[j], t))
                .collect()
        })
        .collect()
}
