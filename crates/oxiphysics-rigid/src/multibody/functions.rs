//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::FbMultibodySystem;

/// 3×3 matrix stored row-major.
pub(super) type Mat3 = [[f64; 3]; 3];
pub(super) fn mat3_identity() -> Mat3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
pub(super) fn mat3_mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
pub(super) fn mat3_vec_mul(m: Mat3, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
pub(super) fn mat3_transpose(m: Mat3) -> Mat3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
/// Rotation matrix from axis-angle representation.
pub(super) fn rot_from_axis_angle(axis: [f64; 3], angle: f64) -> Mat3 {
    let n = vec3_norm(axis);
    if n < 1e-12 {
        return mat3_identity();
    }
    let u = vec3_scale(axis, 1.0 / n);
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    let [ux, uy, uz] = u;
    [
        [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
        [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
        [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
    ]
}
/// Compute the spatial cross-product of two 6-vectors.
///
/// For spatial vectors `v1 = [omega; v]` and `v2 = [omega2; v2]`:
/// `v1 × v2 = [omega × omega2; omega × v2 + v × omega2]`
pub fn spatial_cross_product(v1: &[f64; 6], v2: &[f64; 6]) -> [f64; 6] {
    let w1 = [v1[0], v1[1], v1[2]];
    let u1 = [v1[3], v1[4], v1[5]];
    let w2 = [v2[0], v2[1], v2[2]];
    let u2 = [v2[3], v2[4], v2[5]];
    let top = vec3_cross(w1, w2);
    let bot_a = vec3_cross(w1, u2);
    let bot_b = vec3_cross(u1, w2);
    [
        top[0],
        top[1],
        top[2],
        bot_a[0] + bot_b[0],
        bot_a[1] + bot_b[1],
        bot_a[2] + bot_b[2],
    ]
}
/// Transform a 3×3 inertia tensor `I` by rotation matrix `R`:
/// `I' = R * I * R^T`.
///
/// Both `I` and `R` are passed as row-major flat 9-element arrays.
pub fn inertia_transform(i_mat: &[f64; 9], r_mat: &[f64; 9]) -> [f64; 9] {
    let i3 = [
        [i_mat[0], i_mat[1], i_mat[2]],
        [i_mat[3], i_mat[4], i_mat[5]],
        [i_mat[6], i_mat[7], i_mat[8]],
    ];
    let r3 = [
        [r_mat[0], r_mat[1], r_mat[2]],
        [r_mat[3], r_mat[4], r_mat[5]],
        [r_mat[6], r_mat[7], r_mat[8]],
    ];
    let rt = mat3_transpose(r3);
    let ri = mat3_mul(r3, i3);
    let result = mat3_mul(ri, rt);
    [
        result[0][0],
        result[0][1],
        result[0][2],
        result[1][0],
        result[1][1],
        result[1][2],
        result[2][0],
        result[2][1],
        result[2][2],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::articulated::JointDof;

    use crate::multibody::ArticulatedBody;

    use crate::multibody::FeatherstoneAba;
    use crate::multibody::InverseKinematics;
    use crate::multibody::Joint;
    use crate::multibody::JointDof as MbJointDof;
    use crate::multibody::JointSpaceInertia;
    use crate::multibody::JointType;
    use crate::multibody::KinematicChain;

    use crate::multibody::MultiBodyLink;
    use crate::multibody::MultibodyChain;
    use crate::multibody::MultibodyCollision;
    use crate::multibody::MultibodyKinematics;
    use crate::multibody::MultibodySystem;
    use crate::multibody::MultibodyTree;
    use crate::multibody::OperationalSpaceControl;
    use crate::multibody::RecursiveNewtonEuler;
    use crate::multibody::RigidLink;

    use std::f64::consts::PI;
    #[test]
    fn revolute_has_one_dof() {
        let j = JointDof::Revolute([0.0, 0.0, 1.0]);
        assert_eq!(j.num_dof(), 1);
    }
    #[test]
    fn fixed_has_zero_dof() {
        assert_eq!(JointDof::Fixed.num_dof(), 0);
    }
    #[test]
    fn spherical_has_three_dof() {
        assert_eq!(JointDof::Spherical.num_dof(), 3);
    }
    #[test]
    fn revolute_clamp_in_limits() {
        let j = JointDof::Revolute([0.0, 0.0, 1.0]);
        // Clamp to [-PI, PI] (default revolute limits)
        assert!(j.clamp(2.0 * PI) <= PI);
        assert!(j.clamp(-2.0 * PI) >= -PI);
        assert!((j.clamp(0.5) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn revolute_transform_identity_at_zero() {
        let j = JointDof::Revolute([0.0, 0.0, 1.0]);
        let (rot, trans) = j.transform(0.0);
        assert!((rot[0][0] - 1.0).abs() < 1e-10);
        assert_eq!(trans, [0.0; 3]);
    }
    #[test]
    fn prismatic_transform_translates() {
        let j = JointDof::Prismatic([1.0, 0.0, 0.0]);
        let (_, trans) = j.transform(0.5);
        assert!((trans[0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn rigid_link_kinetic_energy_zero_at_rest() {
        let link = RigidLink::new_box("l0", 1.0, [0.1, 0.2, 0.1]);
        assert_eq!(link.kinetic_energy([0.0; 3], [0.0; 3]), 0.0);
    }
    #[test]
    fn rigid_link_inertia_diagonal_positive() {
        let link = RigidLink::new_box("l0", 1.0, [0.1, 0.2, 0.1]);
        let m = link.inertia_matrix();
        assert!(m[0][0] > 0.0 && m[1][1] > 0.0 && m[2][2] > 0.0);
    }
    #[test]
    fn rigid_link_kinetic_energy_moving() {
        let link = RigidLink::new_box("l0", 2.0, [0.1, 0.1, 0.1]);
        let ke = link.kinetic_energy([1.0, 0.0, 0.0], [0.0; 3]);
        assert!((ke - 1.0).abs() < 1e-10, "ke={ke}");
    }
    #[test]
    fn chain_6dof_arm_has_6_dofs() {
        let chain = MultibodyChain::default_6dof_arm();
        assert_eq!(chain.total_dof(), 6);
    }
    #[test]
    fn chain_6dof_arm_has_6_links() {
        let chain = MultibodyChain::default_6dof_arm();
        assert_eq!(chain.num_links(), 6);
    }
    #[test]
    fn chain_total_mass_positive() {
        let chain = MultibodyChain::default_6dof_arm();
        assert!(chain.total_mass() > 0.0);
    }
    #[test]
    fn chain_set_joint_angles() {
        let mut chain = MultibodyChain::default_6dof_arm();
        chain.set_joint_angles(&[0.1, 0.2, 0.3, 0.1, 0.2, 0.1]);
        assert!((chain.q[0] - 0.1).abs() < 1e-10);
        assert!((chain.q[2] - 0.3).abs() < 1e-10);
    }
    #[test]
    fn humanoid_tree_has_many_nodes() {
        let t = MultibodyTree::humanoid();
        assert!(t.num_nodes() >= 10, "nodes: {}", t.num_nodes());
    }
    #[test]
    fn humanoid_tree_total_mass_positive() {
        let t = MultibodyTree::humanoid();
        assert!(t.total_mass() > 0.0);
    }
    #[test]
    fn tree_add_node_sets_parent() {
        let mut t = MultibodyTree::new();
        let link0 = RigidLink::new_box("root", 10.0, [0.1; 3]);
        let root_id = t.add_node(None, link0, MbJointDof::Fixed);
        let link1 = RigidLink::new_box("child", 5.0, [0.05; 3]);
        let child_id = t.add_node(
            Some(root_id),
            link1,
            MbJointDof::Revolute {
                axis: [0.0, 0.0, 1.0],
                limits: (-PI, PI),
            },
        );
        assert_eq!(t.nodes[child_id].parent, Some(root_id));
    }
    #[test]
    fn aba_forward_dynamics_length_matches_dof() {
        let chain = MultibodyChain::default_6dof_arm();
        let mut aba = FeatherstoneAba::new(chain, [0.0, -9.81, 0.0]);
        let torques = vec![0.0_f64; 6];
        let qdd = aba.forward_dynamics(&torques);
        assert_eq!(qdd.len(), 6);
    }
    #[test]
    fn aba_integrate_changes_positions() {
        let chain = MultibodyChain::default_6dof_arm();
        let q0: Vec<f64> = chain.q.clone();
        let mut aba = FeatherstoneAba::new(chain, [0.0, -9.81, 0.0]);
        let torques = vec![10.0_f64; 6];
        aba.integrate(&torques, 0.01);
        let changed = aba
            .chain
            .q
            .iter()
            .zip(q0.iter())
            .any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(changed, "joint positions should change after integration");
    }
    #[test]
    fn rnea_inverse_dynamics_length_matches_dof() {
        let chain = MultibodyChain::default_6dof_arm();
        let rnea = RecursiveNewtonEuler::new(chain, [0.0, -9.81, 0.0]);
        let tau = rnea.inverse_dynamics();
        assert_eq!(tau.len(), 6);
    }
    #[test]
    fn rnea_gravity_torque_nonzero() {
        let mut chain = MultibodyChain::default_6dof_arm();
        chain.q[1] = 0.5;
        let rnea = RecursiveNewtonEuler::new(chain, [0.0, -9.81, 0.0]);
        let tau = rnea.inverse_dynamics();
        let nonzero = tau.iter().any(|t| t.abs() > 1e-10);
        assert!(nonzero, "gravity torque should be nonzero in non-zero pose");
    }
    #[test]
    fn crba_diagonal_positive() {
        let chain = MultibodyChain::default_6dof_arm();
        let crba = JointSpaceInertia::new(chain);
        let diag = crba.diagonal();
        for (i, d) in diag.iter().enumerate() {
            assert!(*d > 0.0, "diagonal[{i}] should be positive: {d}");
        }
    }
    #[test]
    fn crba_h_diagonal_count() {
        let chain = MultibodyChain::default_6dof_arm();
        let n = chain.total_dof();
        let crba = JointSpaceInertia::new(chain);
        let diag = crba.diagonal();
        assert_eq!(diag.len(), n);
    }
    #[test]
    fn fk_produces_transforms_per_link() {
        let chain = MultibodyChain::default_6dof_arm();
        let kin = MultibodyKinematics::new(chain);
        let fk = kin.forward_kinematics();
        assert_eq!(fk.len(), 6);
    }
    #[test]
    fn fk_identity_at_zero_config() {
        let chain = MultibodyChain::default_6dof_arm();
        let kin = MultibodyKinematics::new(chain);
        let fk = kin.forward_kinematics();
        let rot = fk[0].rotation;
        assert!((rot[0][0] - 1.0).abs() < 1e-6, "R[0][0] should be 1");
    }
    #[test]
    fn ik_step_moves_toward_target() {
        let chain = MultibodyChain::default_6dof_arm();
        let mut kin = MultibodyKinematics::new(chain);
        let target = [0.5, 0.3, 0.2];
        let initial_ee = kin.end_effector_position();
        kin.inverse_kinematics_step(target, 0.1);
        let new_ee = kin.end_effector_position();
        let d_before = vec3_norm(vec3_sub(target, initial_ee));
        let d_after = vec3_norm(vec3_sub(target, new_ee));
        assert!(
            d_after <= d_before + 1e-6,
            "IK should reduce error: {d_before} → {d_after}"
        );
    }
    #[test]
    fn osc_task_force_direction() {
        let chain = MultibodyChain::default_6dof_arm();
        let kin = MultibodyKinematics::new(chain);
        let osc = OperationalSpaceControl::new(kin);
        let target = [1.0, 0.0, 0.0];
        let f = osc.task_force(target);
        assert!(vec3_dot(f, target) > 0.0 || vec3_norm(f) < 1e-10);
    }
    #[test]
    fn osc_null_space_length_matches_dof() {
        let chain = MultibodyChain::default_6dof_arm();
        let n = chain.total_dof();
        let kin = MultibodyKinematics::new(chain);
        let osc = OperationalSpaceControl::new(kin);
        let ns = osc.null_space_projection();
        assert_eq!(ns.len(), n);
    }
    #[test]
    fn collision_no_self_collision_spread_out() {
        let tree = MultibodyTree::humanoid();
        let n = tree.num_nodes();
        let radii = vec![0.05_f64; n];
        let col = MultibodyCollision::new(tree, radii);
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect();
        let results = col.detect(&positions);
        assert!(results.is_empty(), "no collision when spread apart");
    }
    #[test]
    fn collision_detects_overlap() {
        let mut tree = MultibodyTree::new();
        let l0 = RigidLink::new_box("l0", 1.0, [0.1; 3]);
        let l1 = RigidLink::new_box("l1", 1.0, [0.1; 3]);
        let l2 = RigidLink::new_box("l2", 1.0, [0.1; 3]);
        let r0 = tree.add_node(None, l0, MbJointDof::Fixed);
        tree.add_node(Some(r0), l1, MbJointDof::Fixed);
        tree.add_node(None, l2, MbJointDof::Fixed);
        let radii = vec![0.5_f64; 3];
        let col = MultibodyCollision::new(tree, radii);
        let positions = vec![[0.0; 3], [10.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let results = col.detect(&positions);
        assert!(
            !results.is_empty(),
            "should detect collision between 0 and 2"
        );
    }
    #[test]
    fn mat3_identity_mul() {
        let i = mat3_identity();
        let v = [1.0, 2.0, 3.0];
        let result = mat3_vec_mul(i, v);
        assert_eq!(result, v);
    }
    #[test]
    fn rot_from_axis_angle_90deg_z() {
        let r = rot_from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let v = mat3_vec_mul(r, [1.0, 0.0, 0.0]);
        assert!((v[0] - 0.0).abs() < 1e-6, "x: {}", v[0]);
        assert!((v[1] - 1.0).abs() < 1e-6, "y: {}", v[1]);
        assert!((v[2] - 0.0).abs() < 1e-6, "z: {}", v[2]);
    }
    #[test]
    fn joint_type_revolute_one_dof() {
        assert_eq!(JointType::Revolute.degrees_of_freedom(), 1);
    }
    #[test]
    fn joint_type_prismatic_one_dof() {
        assert_eq!(JointType::Prismatic.degrees_of_freedom(), 1);
    }
    #[test]
    fn joint_type_spherical_three_dof() {
        assert_eq!(JointType::Spherical.degrees_of_freedom(), 3);
    }
    #[test]
    fn joint_type_universal_two_dof() {
        assert_eq!(JointType::Universal.degrees_of_freedom(), 2);
    }
    #[test]
    fn joint_type_fixed_zero_dof() {
        assert_eq!(JointType::Fixed.degrees_of_freedom(), 0);
    }
    #[test]
    fn joint_type_planar_three_dof() {
        assert_eq!(JointType::Planar.degrees_of_freedom(), 3);
    }
    #[test]
    fn joint_type_cylindrical_two_dof() {
        assert_eq!(JointType::Cylindrical.degrees_of_freedom(), 2);
    }
    #[test]
    fn joint_type_has_rotation() {
        assert!(JointType::Revolute.has_rotation());
        assert!(JointType::Spherical.has_rotation());
        assert!(!JointType::Prismatic.has_rotation());
        assert!(!JointType::Fixed.has_rotation());
    }
    #[test]
    fn joint_type_has_translation() {
        assert!(JointType::Prismatic.has_translation());
        assert!(JointType::Planar.has_translation());
        assert!(!JointType::Revolute.has_translation());
        assert!(!JointType::Fixed.has_translation());
    }
    #[test]
    fn joint_new_revolute_properties() {
        let j = Joint::new_revolute(0, 1, -PI, PI);
        assert_eq!(j.joint_type, JointType::Revolute);
        assert_eq!(j.parent_body, 0);
        assert_eq!(j.child_body, 1);
        assert_eq!(j.num_dof(), 1);
    }
    #[test]
    fn joint_new_prismatic_properties() {
        let j = Joint::new_prismatic(1, 2, 0.0, 0.5);
        assert_eq!(j.joint_type, JointType::Prismatic);
        assert!((j.upper_limit - 0.5).abs() < 1e-12);
    }
    #[test]
    fn joint_fixed_zero_dof() {
        let j = Joint::new_fixed(0, 1);
        assert_eq!(j.num_dof(), 0);
    }
    #[test]
    fn joint_clamp_position() {
        let j = Joint::new_revolute(0, 1, -1.0, 1.0);
        assert_eq!(j.clamp_position(2.0), 1.0);
        assert_eq!(j.clamp_position(-2.0), -1.0);
        assert!((j.clamp_position(0.5) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn multibody_link_new_sphere_inertia() {
        let link = MultiBodyLink::new_sphere(1.0, 0.5);
        let [ixx, iyy, izz] = link.inertia_diag();
        let expected = 0.4 * 1.0 * 0.25;
        assert!((ixx - expected).abs() < 1e-10);
        assert!((iyy - expected).abs() < 1e-10);
        assert!((izz - expected).abs() < 1e-10);
    }
    #[test]
    fn multibody_link_new_box_positive_inertia() {
        let link = MultiBodyLink::new_box_link(2.0, 0.1, 0.2, 0.1);
        let [ixx, iyy, izz] = link.inertia_diag();
        assert!(ixx > 0.0 && iyy > 0.0 && izz > 0.0);
    }
    #[test]
    fn multibody_link_mass_stored() {
        let link = MultiBodyLink::new(3.0, 1.0, 2.0, 3.0, [0.0; 3]);
        assert!((link.mass - 3.0).abs() < 1e-12);
    }
    #[test]
    fn multibody_system_add_link_increases_ndof() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        assert_eq!(sys.n_dof(), 1);
    }
    #[test]
    fn multibody_system_n_links_correct() {
        let mut sys = MultibodySystem::new();
        sys.add_link(MultiBodyLink::new(1.0, 0.1, 0.1, 0.1, [0.0; 3]));
        sys.add_link(MultiBodyLink::new(2.0, 0.2, 0.2, 0.2, [0.0; 3]));
        assert_eq!(sys.n_links(), 2);
    }
    #[test]
    fn multibody_system_total_mass() {
        let mut sys = MultibodySystem::new();
        sys.add_link(MultiBodyLink::new(1.0, 0.1, 0.1, 0.1, [0.0; 3]));
        sys.add_link(MultiBodyLink::new(2.0, 0.2, 0.2, 0.2, [0.0; 3]));
        assert!((sys.total_mass() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn multibody_system_forward_kinematics_length() {
        let mut sys = MultibodySystem::new();
        sys.add_link(MultiBodyLink::new_sphere(1.0, 0.1));
        sys.add_link(MultiBodyLink::new_sphere(1.0, 0.1));
        let fk = sys.forward_kinematics();
        assert_eq!(fk.len(), 2);
    }
    #[test]
    fn multibody_system_jacobian_shape() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let j = sys.jacobian(0);
        assert_eq!(j.len(), 3);
        assert_eq!(j[0].len(), 1);
    }
    #[test]
    fn multibody_system_mass_matrix_diagonal_positive() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(2.0, 0.2);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let m = sys.mass_matrix();
        assert!(m[0][0] > 0.0, "mass matrix diagonal should be positive");
    }
    #[test]
    fn multibody_system_gravity_vector_nonzero() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let gvec = sys.gravity_vector();
        assert_eq!(gvec.len(), 1);
        assert!(gvec[0].abs() > 0.0, "gravity vector should be nonzero");
    }
    #[test]
    fn multibody_system_forward_dynamics_length() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let qdd = sys.forward_dynamics(&[0.0]);
        assert_eq!(qdd.len(), 1);
    }
    #[test]
    fn multibody_system_inverse_dynamics_consistency() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let tau_in = vec![5.0_f64];
        let qdd = sys.forward_dynamics(&tau_in);
        let tau_out = sys.inverse_dynamics(&sys.positions.clone(), &sys.velocities.clone(), &qdd);
        assert_eq!(tau_out.len(), 1);
        assert!(
            (tau_out[0] - tau_in[0]).abs() < 1e-6,
            "tau roundtrip: {:.4} != {:.4}",
            tau_out[0],
            tau_in[0]
        );
    }
    #[test]
    fn multibody_system_integrate_changes_state() {
        let mut sys = MultibodySystem::new();
        let mut link = MultiBodyLink::new_sphere(1.0, 0.1);
        link.parent_joint = Some(Joint::new_revolute(0, 1, -PI, PI));
        sys.add_link(link);
        let q0 = sys.positions[0];
        sys.integrate_euler(&[1.0], 0.01);
        assert!(
            (sys.positions[0] - q0).abs() > 1e-12,
            "position should change after integration"
        );
    }
    #[test]
    fn kinematic_chain_planar_n_links() {
        let chain = KinematicChain::planar(3, 0.5);
        assert_eq!(chain.n_links, 3);
        assert_eq!(chain.q.len(), 3);
    }
    #[test]
    fn kinematic_chain_dh_transform_identity_at_zero() {
        let chain = KinematicChain::planar(2, 0.0);
        let t = chain.dh_transform(0);
        assert!((t[0][0] - 1.0).abs() < 1e-10);
        assert!((t[1][1] - 1.0).abs() < 1e-10);
        assert!((t[2][2] - 1.0).abs() < 1e-10);
        assert!((t[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn kinematic_chain_end_effector_pose_4x4() {
        let chain = KinematicChain::planar(2, 0.5);
        let pose = chain.end_effector_pose();
        assert!((pose[3][3] - 1.0).abs() < 1e-10);
        assert!((pose[3][0]).abs() < 1e-10);
    }
    #[test]
    fn kinematic_chain_total_reach() {
        let chain = KinematicChain::planar(4, 0.25);
        assert!((chain.total_reach() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn kinematic_chain_jacobian_shape() {
        let chain = KinematicChain::planar(3, 0.5);
        let j = chain.jacobian_matrix();
        assert_eq!(j.len(), 3);
        assert_eq!(j[0].len(), 3);
    }
    #[test]
    fn kinematic_chain_end_effector_moves_with_q() {
        let mut chain = KinematicChain::planar(2, 1.0);
        let ee0 = chain.end_effector_position();
        chain.q[0] = PI / 2.0;
        let ee1 = chain.end_effector_position();
        let moved = (0..3).any(|k| (ee0[k] - ee1[k]).abs() > 1e-6);
        assert!(moved, "end effector should move when joint angle changes");
    }
    #[test]
    fn ik_jacobian_transpose_reachable_target() {
        let mut chain = KinematicChain::planar(3, 0.3);
        let target = [0.5, 0.3, 0.0];
        let _result = InverseKinematics::jacobian_transpose(target, &mut chain, 1e-3, 100);
    }
    #[test]
    fn ik_ccd_reachable_target() {
        let mut chain = KinematicChain::planar(3, 0.3);
        let target = [0.4, 0.2, 0.0];
        let _result = InverseKinematics::cyclic_coordinate_descent(target, &mut chain, 1e-3, 100);
    }
    #[test]
    fn ik_fabrik_basic() {
        let mut chain = KinematicChain::planar(3, 0.3);
        let target = [0.8, 0.0, 0.0];
        let _result = InverseKinematics::fabrik(target, &mut chain, 1e-2, 500);
        assert_eq!(
            chain.n_links, 3,
            "chain should still have 3 links after FABRIK"
        );
    }
    #[test]
    fn ik_fabrik_unreachable_stretches() {
        let mut chain = KinematicChain::planar(2, 0.5);
        let target = [2.0, 0.0, 0.0];
        let _result = InverseKinematics::fabrik(target, &mut chain, 1e-3, 50);
    }
    #[test]
    fn articulated_body_new_correct_counts() {
        let ab = ArticulatedBody::new(4);
        assert_eq!(ab.n_joints, 4);
        assert_eq!(ab.joint_angles.len(), 4);
        assert_eq!(ab.joint_velocities.len(), 4);
    }
    #[test]
    fn articulated_body_step_changes_angles() {
        let mut ab = ArticulatedBody::new(2);
        let initial = ab.joint_angles.clone();
        ab.step(&[1.0, 1.0], 0.01);
        let changed = ab
            .joint_angles
            .iter()
            .zip(initial.iter())
            .any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(changed, "joint angles should change after step");
    }
    #[test]
    fn articulated_body_kinetic_energy_zero_at_rest() {
        let ab = ArticulatedBody::new(3);
        assert!(ab.kinetic_energy().abs() < 1e-12);
    }
    #[test]
    fn articulated_body_kinetic_energy_positive_when_moving() {
        let mut ab = ArticulatedBody::new(2);
        ab.joint_velocities[0] = 1.0;
        assert!(ab.kinetic_energy() > 0.0);
    }
    #[test]
    fn articulated_body_potential_energy_zero_at_zero_angle() {
        let ab = ArticulatedBody::new(2);
        assert!(ab.potential_energy().abs() < 1e-12);
    }
    #[test]
    fn articulated_body_clamp_angles_wraps_to_pi() {
        let mut ab = ArticulatedBody::new(1);
        ab.joint_angles[0] = 4.0 * PI;
        ab.clamp_angles();
        assert!(ab.joint_angles[0].abs() < PI + 1e-10);
    }
    #[test]
    fn articulated_body_step_with_spring_restores() {
        let mut ab = ArticulatedBody::new(1);
        ab.joint_stiffness[0] = 10.0;
        ab.rest_angles[0] = 0.0;
        ab.joint_angles[0] = 1.0;
        ab.joint_damping[0] = 0.5;
        for _ in 0..100 {
            ab.step(&[0.0], 0.01);
        }
        assert!(
            ab.joint_angles[0].abs() < 1.0,
            "spring should pull angle toward 0"
        );
    }
}
/// Compute the centre-of-mass positions of all links in world space.
///
/// Returns one `[f64; 3]` per link.  The root link's CoM is placed at the
/// origin; each subsequent link is offset along its joint axis by
/// `joint_pos`.
pub fn forward_kinematics(sys: &FbMultibodySystem) -> Vec<[f64; 3]> {
    let n = sys.links.len();
    let mut positions = vec![[0.0_f64; 3]; n];
    for (i, link) in sys.links.iter().enumerate() {
        let parent_pos = if link.parent < 0 {
            [0.0_f64; 3]
        } else {
            positions[link.parent as usize]
        };
        let axis = link.joint_axis;
        let q = link.joint_pos;
        positions[i] = [
            parent_pos[0] + axis[0] * q,
            parent_pos[1] + axis[1] * q,
            parent_pos[2] + axis[2] * q,
        ];
    }
    positions
}
/// Compute the Jacobian column for a single link (the 3×n_dof partial
/// derivatives of `point` with respect to each joint variable).
///
/// Returns one `[f64; 3]` per degree of freedom (n_dof entries).
/// For revolute joints the column is `joint_axis × (point − joint_origin)`;
/// for prismatic joints it is simply `joint_axis`.
///
/// This simplified implementation treats every joint as prismatic (translation
/// along its axis), which is correct for a chain of prismatic joints.
pub fn jacobian_column(
    sys: &FbMultibodySystem,
    _link_idx: usize,
    _point: [f64; 3],
) -> Vec<[f64; 3]> {
    let mut cols = Vec::with_capacity(sys.n_dof);
    for link in &sys.links {
        if link.parent >= 0 {
            cols.push(link.joint_axis);
        }
    }
    cols
}
/// Compute the approximate joint-space mass matrix (n_dof × n_dof).
///
/// Uses a diagonal approximation: `M[i][i] = link.mass` for each DOF link.
/// Off-diagonal couplings are omitted for this simplified model.
pub fn mass_matrix(sys: &FbMultibodySystem) -> Vec<Vec<f64>> {
    let n = sys.n_dof;
    let mut m = vec![vec![0.0_f64; n]; n];
    let mut dof_idx = 0usize;
    for link in &sys.links {
        if link.parent >= 0 {
            m[dof_idx][dof_idx] = link.mass;
            dof_idx += 1;
        }
    }
    m
}
/// Compute the Coriolis/centrifugal vector **C(q, dq)·dq** (n_dof × 1).
///
/// In this simplified single-axis joint model each DOF is treated as independent,
/// so the Coriolis coupling between joints is neglected and each entry is zero.
/// A more complete implementation would account for cross-joint inertia coupling.
pub fn coriolis_vector(sys: &FbMultibodySystem) -> Vec<f64> {
    vec![0.0; sys.n_dof]
}
/// Compute the gravity generalised-force vector **G(q)** (n_dof × 1).
///
/// `G_i = -mass_i · (g · joint_axis_i)`
pub fn gravity_vector(sys: &FbMultibodySystem, g: [f64; 3]) -> Vec<f64> {
    let mut gv = Vec::with_capacity(sys.n_dof);
    for link in &sys.links {
        if link.parent >= 0 {
            let g_proj =
                g[0] * link.joint_axis[0] + g[1] * link.joint_axis[1] + g[2] * link.joint_axis[2];
            gv.push(-link.mass * g_proj);
        }
    }
    gv
}
/// Compute joint accelerations from applied torques/forces.
///
/// Solves `M · qdd = tau − C − G` for `qdd` using the diagonal mass matrix.
pub fn forward_dynamics(sys: &FbMultibodySystem, tau: &[f64], g: [f64; 3]) -> Vec<f64> {
    let c = coriolis_vector(sys);
    let gv = gravity_vector(sys, g);
    let mut qdd = Vec::with_capacity(sys.n_dof);
    let mut dof_idx = 0usize;
    for link in &sys.links {
        if link.parent >= 0 {
            let tau_i = tau.get(dof_idx).copied().unwrap_or(0.0);
            let rhs = tau_i - c[dof_idx] - gv[dof_idx];
            let m_ii = link.mass.max(1e-15);
            qdd.push(rhs / m_ii);
            dof_idx += 1;
        }
    }
    qdd
}
/// Integrate the multibody system by one time step `dt` (semi-implicit Euler).
pub fn integrate_multibody(sys: &mut FbMultibodySystem, tau: &[f64], g: [f64; 3], dt: f64) {
    let qdd = forward_dynamics(sys, tau, g);
    let mut dof_idx = 0usize;
    for link in sys.links.iter_mut() {
        if link.parent >= 0 {
            let acc = qdd.get(dof_idx).copied().unwrap_or(0.0);
            link.joint_vel += acc * dt;
            link.joint_pos += link.joint_vel * dt;
            dof_idx += 1;
        }
    }
}
/// Compute the total kinetic energy of the system (J).
///
/// `KE = Σ 0.5 · mass_i · joint_vel_i²`
pub fn kinetic_energy_multibody(sys: &FbMultibodySystem) -> f64 {
    sys.links
        .iter()
        .filter(|l| l.parent >= 0)
        .map(|l| 0.5 * l.mass * l.joint_vel * l.joint_vel)
        .sum()
}
#[cfg(test)]
mod simple_api_tests {
    use super::*;

    use crate::multibody::FbMultibodySystem;

    use crate::multibody::Link;

    use crate::multibody::coriolis_vector;
    use crate::multibody::forward_dynamics;
    use crate::multibody::gravity_vector;
    use crate::multibody::integrate_multibody;
    use crate::multibody::jacobian_column;
    use crate::multibody::kinetic_energy_multibody;
    use crate::multibody::mass_matrix;

    pub(super) const EPS: f64 = 1e-10;
    fn two_link_system() -> FbMultibodySystem {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [0.0, 0.0, 1.0], 0.0, 0.0));
        sys.add_link(Link::new(
            2.0,
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            0,
            [0.0, 0.0, 1.0],
            0.5,
            1.0,
        ));
        sys
    }
    #[test]
    fn link_new_stores_fields() {
        let l = Link::new(3.0, [1.0; 6], 0, [1.0, 0.0, 0.0], 0.5, 2.0);
        assert!((l.mass - 3.0).abs() < EPS);
        assert_eq!(l.parent, 0);
        assert!((l.joint_pos - 0.5).abs() < EPS);
        assert!((l.joint_vel - 2.0).abs() < EPS);
    }
    #[test]
    fn link_revolute_z_defaults() {
        let l = Link::revolute_z(5.0, 0.2, -1);
        assert_eq!(l.parent, -1);
        assert!((l.joint_axis[2] - 1.0).abs() < EPS);
        assert!((l.inertia[2] - 0.2).abs() < EPS);
    }
    #[test]
    fn system_new_empty() {
        let sys = FbMultibodySystem::new();
        assert_eq!(sys.num_links(), 0);
        assert_eq!(sys.n_dof, 0);
    }
    #[test]
    fn system_add_root_no_dof() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [0.0, 0.0, 1.0], 0.0, 0.0));
        assert_eq!(sys.n_dof, 0);
        assert_eq!(sys.num_links(), 1);
    }
    #[test]
    fn system_add_child_increments_dof() {
        let sys = two_link_system();
        assert_eq!(sys.n_dof, 1);
        assert_eq!(sys.num_links(), 2);
    }
    #[test]
    fn system_total_mass() {
        let sys = two_link_system();
        assert!((sys.total_mass() - 3.0).abs() < EPS);
    }
    #[test]
    fn system_default_same_as_new() {
        let s1 = FbMultibodySystem::new();
        let s2 = FbMultibodySystem::default();
        assert_eq!(s1.num_links(), s2.num_links());
    }
    #[test]
    fn fk_root_at_origin() {
        let sys = two_link_system();
        let pos = forward_kinematics(&sys);
        for c in pos[0] {
            assert!(c.abs() < EPS, "root should be at origin");
        }
    }
    #[test]
    fn fk_child_offset_by_joint_pos() {
        let sys = two_link_system();
        let pos = forward_kinematics(&sys);
        assert!((pos[1][2] - 0.5).abs() < EPS, "child z={}", pos[1][2]);
    }
    #[test]
    fn fk_three_link_chain() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [1.0, 0.0, 0.0], 0.0, 0.0));
        sys.add_link(Link::new(1.0, [0.0; 6], 0, [1.0, 0.0, 0.0], 1.0, 0.0));
        sys.add_link(Link::new(1.0, [0.0; 6], 1, [1.0, 0.0, 0.0], 1.0, 0.0));
        let pos = forward_kinematics(&sys);
        assert!((pos[2][0] - 2.0).abs() < EPS, "tip x={}", pos[2][0]);
    }
    #[test]
    fn fk_returns_correct_count() {
        let sys = two_link_system();
        assert_eq!(forward_kinematics(&sys).len(), 2);
    }
    #[test]
    fn jacobian_column_size_equals_ndof() {
        let sys = two_link_system();
        let jac = jacobian_column(&sys, 1, [0.0, 0.0, 0.5]);
        assert_eq!(jac.len(), sys.n_dof);
    }
    #[test]
    fn jacobian_column_contains_joint_axis() {
        let sys = two_link_system();
        let jac = jacobian_column(&sys, 1, [0.0, 0.0, 0.0]);
        assert!((jac[0][2] - 1.0).abs() < EPS, "expected z-axis");
    }
    #[test]
    fn mass_matrix_size() {
        let sys = two_link_system();
        let m = mass_matrix(&sys);
        assert_eq!(m.len(), sys.n_dof);
        assert_eq!(m[0].len(), sys.n_dof);
    }
    #[test]
    fn mass_matrix_diagonal_equals_link_mass() {
        let sys = two_link_system();
        let m = mass_matrix(&sys);
        assert!((m[0][0] - 2.0).abs() < EPS, "M[0][0]={}", m[0][0]);
    }
    #[test]
    fn mass_matrix_positive_definite() {
        let sys = two_link_system();
        let m = mass_matrix(&sys);
        for (i, row) in m.iter().enumerate() {
            assert!(row[i] > 0.0);
        }
    }
    #[test]
    fn coriolis_vector_zero_at_rest() {
        let mut sys = two_link_system();
        sys.links[1].joint_vel = 0.0;
        let c = coriolis_vector(&sys);
        assert!(c[0].abs() < EPS);
    }
    #[test]
    fn coriolis_vector_zero_simplified_model() {
        let sys = two_link_system();
        let c = coriolis_vector(&sys);
        assert!(c[0].abs() < EPS, "simplified model returns zero coriolis");
    }
    #[test]
    fn coriolis_vector_length_equals_ndof() {
        let sys = two_link_system();
        assert_eq!(coriolis_vector(&sys).len(), sys.n_dof);
    }
    #[test]
    fn gravity_vector_zero_perpendicular_to_axis() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [0.0, 0.0, 1.0], 0.0, 0.0));
        sys.add_link(Link::new(1.0, [0.0; 6], 0, [1.0, 0.0, 0.0], 0.0, 0.0));
        let gv = gravity_vector(&sys, [0.0, 0.0, -9.81]);
        assert!(gv[0].abs() < EPS, "no gravity component along x-axis");
    }
    #[test]
    fn gravity_vector_nonzero_parallel_to_axis() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [0.0, 0.0, 1.0], 0.0, 0.0));
        sys.add_link(Link::new(1.0, [0.0; 6], 0, [0.0, 0.0, 1.0], 0.0, 0.0));
        let gv = gravity_vector(&sys, [0.0, 0.0, -9.81]);
        assert!(gv[0].abs() > 0.0);
    }
    #[test]
    fn gravity_vector_length_equals_ndof() {
        let sys = two_link_system();
        assert_eq!(gravity_vector(&sys, [0.0, -9.81, 0.0]).len(), sys.n_dof);
    }
    #[test]
    fn forward_dynamics_length() {
        let sys = two_link_system();
        let qdd = forward_dynamics(&sys, &[0.0], [0.0, -9.81, 0.0]);
        assert_eq!(qdd.len(), sys.n_dof);
    }
    #[test]
    fn forward_dynamics_zero_tau_gravity_driven() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [0.0, 0.0, 1.0], 0.0, 0.0));
        sys.add_link(Link::new(2.0, [0.0; 6], 0, [0.0, 0.0, 1.0], 0.0, 0.0));
        let g = [0.0, 0.0, -9.81];
        let qdd = forward_dynamics(&sys, &[0.0], g);
        assert!((qdd[0] + 9.81).abs() < 1e-6, "qdd={}", qdd[0]);
    }
    #[test]
    fn forward_dynamics_with_torque_changes_accel() {
        let sys = two_link_system();
        let qdd_no_tau = forward_dynamics(&sys, &[0.0], [0.0, 0.0, 0.0]);
        let qdd_tau = forward_dynamics(&sys, &[10.0], [0.0, 0.0, 0.0]);
        assert!((qdd_tau[0] - qdd_no_tau[0]).abs() > EPS);
    }
    #[test]
    fn integrate_changes_position() {
        let mut sys = two_link_system();
        let q0 = sys.links[1].joint_pos;
        integrate_multibody(&mut sys, &[1.0], [0.0, 0.0, 0.0], 0.01);
        assert!((sys.links[1].joint_pos - q0).abs() > EPS);
    }
    #[test]
    fn integrate_changes_velocity() {
        let mut sys = two_link_system();
        let v0 = sys.links[1].joint_vel;
        integrate_multibody(&mut sys, &[1.0], [0.0, 0.0, 0.0], 0.01);
        let changed = (sys.links[1].joint_vel - v0).abs() > EPS;
        assert!(changed, "velocity should change under applied torque");
    }
    #[test]
    fn integrate_zero_tau_zero_gravity_constant_velocity() {
        let mut sys = FbMultibodySystem::new();
        sys.add_link(Link::new(1.0, [0.0; 6], -1, [1.0, 0.0, 0.0], 0.0, 0.0));
        sys.add_link(Link::new(1.0, [0.0; 6], 0, [1.0, 0.0, 0.0], 0.0, 2.0));
        let v0 = sys.links[1].joint_vel;
        integrate_multibody(&mut sys, &[0.0], [0.0, 0.0, 0.0], 0.01);
        assert!(
            (sys.links[1].joint_vel - v0).abs() < EPS,
            "no force → constant velocity"
        );
    }
    #[test]
    fn kinetic_energy_zero_at_rest() {
        let mut sys = two_link_system();
        sys.links[1].joint_vel = 0.0;
        assert!(kinetic_energy_multibody(&sys).abs() < EPS);
    }
    #[test]
    fn kinetic_energy_positive_when_moving() {
        let sys = two_link_system();
        assert!(kinetic_energy_multibody(&sys) > 0.0);
    }
    #[test]
    fn kinetic_energy_formula() {
        let sys = two_link_system();
        let ke = kinetic_energy_multibody(&sys);
        assert!((ke - 1.0).abs() < EPS, "KE=0.5*2*1=1.0, got {ke}");
    }
    #[test]
    fn kinetic_energy_increases_with_speed() {
        let mut sys = two_link_system();
        sys.links[1].joint_vel = 2.0;
        let ke2 = kinetic_energy_multibody(&sys);
        sys.links[1].joint_vel = 1.0;
        let ke1 = kinetic_energy_multibody(&sys);
        assert!(ke2 > ke1);
    }
}
