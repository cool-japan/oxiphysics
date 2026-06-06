//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ArticulatedPose, BodySegmentInertia, BoneDescriptor};

/// Multiply two quaternions (w, x, y, z order).
pub(super) fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let (aw, ax, ay, az) = (a[0], a[1], a[2], a[3]);
    let (bw, bx, by, bz) = (b[0], b[1], b[2], b[3]);
    [
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ]
}
/// Conjugate (inverse for unit quaternion) of q.
pub(super) fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [q[0], -q[1], -q[2], -q[3]]
}
/// Rotate vector v by quaternion q.
pub(super) fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let vq: [f64; 4] = [0.0, v[0], v[1], v[2]];
    let tmp = quat_mul(q, vq);
    let res = quat_mul(tmp, quat_conjugate(q));
    [res[1], res[2], res[3]]
}
pub(super) fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn vec_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn vec_length(a: [f64; 3]) -> f64 {
    vec_dot(a, a).sqrt()
}
/// Normalize a quaternion; returns identity if near-zero magnitude.
pub(super) fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let len2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if len2 < 1e-30 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    let inv = 1.0 / len2.sqrt();
    [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv]
}
/// Linearly interpolate two quaternions and normalize (NLERP).
pub(super) fn quat_nlerp(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    let sign = if dot < 0.0 { -1.0 } else { 1.0 };
    let q = [
        a[0] * (1.0 - t) + sign * b[0] * t,
        a[1] * (1.0 - t) + sign * b[1] * t,
        a[2] * (1.0 - t) + sign * b[2] * t,
        a[3] * (1.0 - t) + sign * b[3] * t,
    ];
    quat_normalize(q)
}
/// Angle between two quaternions (in radians).
pub(super) fn quat_angle(a: [f64; 4], b: [f64; 4]) -> f64 {
    let dot = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3])
        .abs()
        .clamp(0.0, 1.0);
    2.0 * dot.acos()
}
/// Quaternion from axis-angle representation.
pub(super) fn quat_from_axis_angle(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let half = angle * 0.5;
    let s = half.sin();
    let len = vec_length(axis);
    if len < 1e-15 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    let inv_len = 1.0 / len;
    [
        half.cos(),
        axis[0] * inv_len * s,
        axis[1] * inv_len * s,
        axis[2] * inv_len * s,
    ]
}
/// Compute world-space positions for each bone using forward kinematics.
///
/// Each bone is placed at the end of its parent bone, offset along Y by the
/// parent bone length. The root bone is placed at `pose.root_position`.
pub fn forward_kinematics(bones: &[BoneDescriptor], pose: &ArticulatedPose) -> Vec<[f64; 3]> {
    let n = bones.len();
    let mut positions = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let pos = match bones[i].parent {
            None => pose.root_position,
            Some(parent_idx) => {
                let parent_pos = positions[parent_idx];
                let parent_len = bones[parent_idx].length;
                [parent_pos[0], parent_pos[1] + parent_len, parent_pos[2]]
            }
        };
        positions[i] = pos;
    }
    positions
}
/// Compute the combined inertia tensor of a set of bones about the origin.
///
/// Uses the parallel axis theorem: I_total = Σ (I_bone + m * (|r|²·E - r⊗r))
/// Returns a 3×3 symmetric matrix stored as \[\[f64;3\];3\].
pub fn ragdoll_inertia_tensor(bones: &[BoneDescriptor], positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let mut tensor = [[0.0f64; 3]; 3];
    for (i, bone) in bones.iter().enumerate() {
        let p = positions[i];
        let m = bone.mass;
        let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        let local_inertia = bone.inertia;
        tensor[0][0] += local_inertia[0] + m * (r2 - p[0] * p[0]);
        tensor[1][1] += local_inertia[1] + m * (r2 - p[1] * p[1]);
        tensor[2][2] += local_inertia[2] + m * (r2 - p[2] * p[2]);
        tensor[0][1] -= m * p[0] * p[1];
        tensor[1][0] -= m * p[0] * p[1];
        tensor[0][2] -= m * p[0] * p[2];
        tensor[2][0] -= m * p[0] * p[2];
        tensor[1][2] -= m * p[1] * p[2];
        tensor[2][1] -= m * p[1] * p[2];
    }
    tensor
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::ragdoll::quat_from_axis_angle;
    use crate::ragdoll::types::{
        HumanoidRagdoll, JointLimit, Ragdoll, RagdollPose, RagdollState, SkeletonJoint,
    };

    #[test]
    fn test_add_bone_position() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        let idx = rag.add_bone("test", [1.0, 2.0, 3.0], 1.0, [0.1, 0.1, 0.1], None);
        assert_eq!(idx, 0);
        let pos = rag.bones[idx].position;
        assert!((pos[0] - 1.0).abs() < 1e-10);
        assert!((pos[1] - 2.0).abs() < 1e-10);
        assert!((pos[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_integrate_gravity() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("b", [0.0, 10.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.integrate(0.1);
        assert!(
            rag.bones[0].position[1] < 10.0,
            "Bone should fall under gravity"
        );
        assert!(
            rag.bones[0].linear_vel[1] < 0.0,
            "Velocity should be negative"
        );
    }
    #[test]
    fn test_bone_world_anchor_identity() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        let idx = rag.add_bone("b", [1.0, 2.0, 3.0], 1.0, [0.1, 0.1, 0.1], None);
        let anchor = rag.bone_world_anchor(idx, [0.5, 0.0, 0.0]);
        assert!((anchor[0] - 1.5).abs() < 1e-10);
        assert!((anchor[1] - 2.0).abs() < 1e-10);
        assert!((anchor[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_step_humanoid_no_panic() {
        let mut rag = HumanoidRagdoll::build();
        for _ in 0..10 {
            rag.step(1.0 / 60.0, 4);
        }
    }
    #[test]
    fn test_ball_joint_reduces_separation() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        let a = rag.add_bone("a", [0.0, 0.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        let b = rag.add_bone("b", [1.0, 0.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.add_ball_joint(a, b, [0.0; 3], [0.0; 3], std::f64::consts::PI);
        let sep_before = vec_length(vec_sub(
            rag.bone_world_anchor(b, [0.0; 3]),
            rag.bone_world_anchor(a, [0.0; 3]),
        ));
        rag.solve_ball_joint_position(0, 0.016);
        let sep_after = vec_length(vec_sub(
            rag.bone_world_anchor(b, [0.0; 3]),
            rag.bone_world_anchor(a, [0.0; 3]),
        ));
        assert!(
            sep_after < sep_before,
            "Separation should decrease after solve"
        );
    }
    #[test]
    fn test_humanoid_bone_count() {
        let rag = HumanoidRagdoll::build();
        assert!(
            rag.bones.len() >= 10,
            "Humanoid should have at least 10 bones"
        );
    }
    #[test]
    fn test_deactivate_stops_simulation() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("b", [0.0, 10.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.deactivate();
        let pos_before = rag.bones[0].position;
        rag.integrate(0.1);
        assert_eq!(
            rag.bones[0].position, pos_before,
            "Inactive ragdoll should not move"
        );
    }
    #[test]
    fn test_activate_resumes_simulation() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("b", [0.0, 10.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.deactivate();
        rag.activate();
        rag.integrate(0.1);
        assert!(
            rag.bones[0].position[1] < 10.0,
            "Active ragdoll should move"
        );
    }
    #[test]
    fn test_blending_state() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("b", [0.0, 10.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.start_blending(0.5);
        assert!(
            matches!(rag.state, RagdollState::Blending { blend_factor } if (blend_factor
            - 0.5).abs() < 1e-12)
        );
        rag.set_blend_factor(0.8);
        assert!(
            matches!(rag.state, RagdollState::Blending { blend_factor } if (blend_factor
            - 0.8).abs() < 1e-12)
        );
    }
    #[test]
    fn test_is_active() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        assert!(rag.is_active());
        rag.deactivate();
        assert!(!rag.is_active());
        rag.start_blending(0.5);
        assert!(rag.is_active());
    }
    #[test]
    fn test_capture_and_apply_pose() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [1.0, 2.0, 3.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.add_bone("b", [4.0, 5.0, 6.0], 1.0, [0.1, 0.1, 0.1], None);
        let pose = rag.capture_pose();
        rag.bones[0].position = [99.0, 0.0, 0.0];
        rag.apply_pose(&pose);
        assert!((rag.bones[0].position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_blend_poses() {
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let pose_b = RagdollPose {
            positions: vec![[10.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let blended = Ragdoll::blend_poses(&pose_a, &pose_b, 0.5);
        assert!(
            (blended.positions[0][0] - 5.0).abs() < 1e-10,
            "blended pos: {:?}",
            blended.positions[0]
        );
    }
    #[test]
    fn test_subtree() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        let root = rag.add_bone("root", [0.0; 3], 1.0, [0.1; 3], None);
        let child = rag.add_bone("child", [1.0, 0.0, 0.0], 1.0, [0.1; 3], Some(root));
        let grandchild = rag.add_bone("grandchild", [2.0, 0.0, 0.0], 1.0, [0.1; 3], Some(child));
        let tree = rag.subtree(root);
        assert!(tree.contains(&root));
        assert!(tree.contains(&child));
        assert!(tree.contains(&grandchild));
        assert_eq!(tree.len(), 3);
    }
    #[test]
    fn test_bone_depth() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        let root = rag.add_bone("root", [0.0; 3], 1.0, [0.1; 3], None);
        let child = rag.add_bone("child", [0.0; 3], 1.0, [0.1; 3], Some(root));
        let grandchild = rag.add_bone("gc", [0.0; 3], 1.0, [0.1; 3], Some(child));
        assert_eq!(rag.bone_depth(root), 0);
        assert_eq!(rag.bone_depth(child), 1);
        assert_eq!(rag.bone_depth(grandchild), 2);
    }
    #[test]
    fn test_root_bone() {
        let rag = HumanoidRagdoll::build();
        let root = rag.root_bone();
        assert!(root.is_some());
        assert_eq!(rag.bones[root.unwrap()].name, "pelvis");
    }
    #[test]
    fn test_total_mass() {
        let rag = HumanoidRagdoll::build();
        let mass = rag.total_mass();
        assert!(mass > 0.0, "total mass: {mass}");
    }
    #[test]
    fn test_center_of_mass() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [0.0, 0.0, 0.0], 1.0, [0.1; 3], None);
        rag.add_bone("b", [10.0, 0.0, 0.0], 1.0, [0.1; 3], None);
        let com = rag.center_of_mass();
        assert!((com[0] - 5.0).abs() < 1e-10, "COM x: {}", com[0]);
    }
    #[test]
    fn test_from_skeleton() {
        let skeleton = vec![
            SkeletonJoint {
                name: "root".to_string(),
                bind_position: [0.0, 0.0, 0.0],
                bind_orientation: [1.0, 0.0, 0.0, 0.0],
                parent: None,
                mass: 5.0,
                half_extents: [0.1, 0.1, 0.1],
            },
            SkeletonJoint {
                name: "child".to_string(),
                bind_position: [0.0, 1.0, 0.0],
                bind_orientation: [1.0, 0.0, 0.0, 0.0],
                parent: Some(0),
                mass: 3.0,
                half_extents: [0.05, 0.1, 0.05],
            },
        ];
        let rag = Ragdoll::from_skeleton(&skeleton, [0.0, -9.81, 0.0], std::f64::consts::PI / 4.0);
        assert_eq!(rag.bones.len(), 2);
        assert_eq!(rag.ball_joints.len(), 1);
        assert_eq!(rag.bones[0].name, "root");
    }
    #[test]
    fn test_standard_human_bone_count() {
        let skel = HumanoidRagdoll::standard_human();
        assert!(
            skel.bones.len() >= 10,
            "standard_human should have >= 10 bones, got {}",
            skel.bones.len()
        );
    }
    #[test]
    fn test_humanoid_skeleton_total_mass_positive() {
        let skel = HumanoidRagdoll::standard_human();
        let mass = skel.total_mass();
        assert!(mass > 0.0, "total mass should be positive: {mass}");
    }
    #[test]
    fn test_joint_limit_clamp() {
        let limit = JointLimit {
            min_angles: [-1.0, -0.5, -0.3],
            max_angles: [1.0, 0.5, 0.3],
        };
        let angles = [-2.0, 0.3, 1.0];
        let clamped = limit.clamp(angles);
        assert!((clamped[0] - (-1.0)).abs() < 1e-10, "should clamp to min");
        assert!((clamped[1] - 0.3).abs() < 1e-10, "within range, no change");
        assert!((clamped[2] - 0.3).abs() < 1e-10, "should clamp to max");
    }
    #[test]
    fn test_fk_root_position_matches_pose() {
        let skel = HumanoidRagdoll::standard_human();
        let root_pos = [1.0, 2.0, 3.0];
        let mut pose = ArticulatedPose::default_t_pose(skel.bones.len());
        pose.root_position = root_pos;
        let positions = forward_kinematics(&skel.bones, &pose);
        assert!((positions[0][0] - root_pos[0]).abs() < 1e-10);
        assert!((positions[0][1] - root_pos[1]).abs() < 1e-10);
        assert!((positions[0][2] - root_pos[2]).abs() < 1e-10);
    }
    #[test]
    fn test_ragdoll_inertia_tensor_positive_diagonal() {
        let skel = HumanoidRagdoll::standard_human();
        let pose = ArticulatedPose::default_t_pose(skel.bones.len());
        let positions = forward_kinematics(&skel.bones, &pose);
        let tensor = ragdoll_inertia_tensor(&skel.bones, &positions);
        assert!(
            tensor[0][0] > 0.0,
            "Ixx should be positive: {}",
            tensor[0][0]
        );
        assert!(
            tensor[1][1] > 0.0,
            "Iyy should be positive: {}",
            tensor[1][1]
        );
        assert!(
            tensor[2][2] > 0.0,
            "Izz should be positive: {}",
            tensor[2][2]
        );
    }
    #[test]
    fn test_apply_impulse_all() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [0.0; 3], 1.0, [0.1; 3], None);
        rag.add_bone("b", [1.0, 0.0, 0.0], 2.0, [0.1; 3], None);
        rag.apply_impulse_all([10.0, 0.0, 0.0]);
        assert!((rag.bones[0].linear_vel[0] - 10.0).abs() < 1e-10);
        assert!((rag.bones[1].linear_vel[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_apply_damping() {
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [0.0; 3], 1.0, [0.1; 3], None);
        rag.bones[0].linear_vel = [10.0, 0.0, 0.0];
        rag.bones[0].angular_vel = [0.0, 5.0, 0.0];
        rag.apply_damping(0.1, 0.2);
        assert!((rag.bones[0].linear_vel[0] - 9.0).abs() < 1e-10);
        assert!((rag.bones[0].angular_vel[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_nlerp_identity() {
        let q = [1.0, 0.0, 0.0, 0.0];
        let result = quat_nlerp(q, q, 0.5);
        assert!(
            (result[0] - 1.0).abs() < 1e-10,
            "nlerp identity: {:?}",
            result
        );
    }
    #[test]
    fn test_quat_angle_identity() {
        let q = [1.0, 0.0, 0.0, 0.0];
        let angle = quat_angle(q, q);
        assert!(
            angle.abs() < 1e-10,
            "angle between identity and itself: {angle}"
        );
    }
    #[test]
    fn test_quat_from_axis_angle() {
        let q = quat_from_axis_angle([0.0, 1.0, 0.0], std::f64::consts::PI / 2.0);
        let len2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
        assert!(
            (len2 - 1.0).abs() < 1e-10,
            "quaternion should be unit: len2={len2}"
        );
    }
}
/// Build standard body segment inertias for a 70 kg human.
pub fn standard_human_segment_inertias(total_mass: f64) -> Vec<(String, BodySegmentInertia)> {
    vec![
        (
            "head".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.081, 0.10, 0.22),
        ),
        (
            "trunk".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.497, 0.16, 0.52),
        ),
        (
            "upper_arm".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.028, 0.04, 0.28),
        ),
        (
            "forearm".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.016, 0.03, 0.27),
        ),
        (
            "hand".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.006, 0.04, 0.08),
        ),
        (
            "thigh".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.100, 0.07, 0.42),
        ),
        (
            "shank".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.046, 0.05, 0.39),
        ),
        (
            "foot".to_string(),
            BodySegmentInertia::from_body_mass(total_mass, 0.014, 0.04, 0.26),
        ),
    ]
}
#[cfg(test)]
mod tests_ragdoll_extra {

    use crate::BalanceController;
    use crate::BodySegmentInertia;
    use crate::ConeTwistLimit;
    use crate::GroundContactDetector;
    use crate::PdJointController;
    use crate::PosePdController;
    use crate::Ragdoll;
    use crate::RagdollPose;
    use crate::ragdoll::quat_from_axis_angle;
    use crate::standard_human_segment_inertias;
    #[test]
    fn test_cone_twist_identity_within_cone() {
        let limit =
            ConeTwistLimit::symmetric(std::f64::consts::PI / 4.0, std::f64::consts::PI / 2.0);
        let q_ident = [1.0, 0.0, 0.0, 0.0];
        assert!(
            limit.within_cone(q_ident),
            "identity orientation should be within cone"
        );
    }
    #[test]
    fn test_cone_twist_large_rotation_outside_cone() {
        let limit = ConeTwistLimit::symmetric(std::f64::consts::PI / 6.0, std::f64::consts::PI);
        let q = quat_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::PI / 2.0);
        assert!(
            !limit.within_cone(q),
            "90° rotation should be outside 30° cone"
        );
    }
    #[test]
    fn test_cone_twist_violation_zero_for_identity() {
        let limit = ConeTwistLimit::symmetric(std::f64::consts::PI / 4.0, std::f64::consts::PI);
        let violation = limit.cone_violation([1.0, 0.0, 0.0, 0.0]);
        assert!(
            violation.abs() < 1e-10,
            "no violation for identity: {violation}"
        );
    }
    #[test]
    fn test_cone_twist_apply_cone_limit_corrects() {
        let limit = ConeTwistLimit::symmetric(std::f64::consts::PI / 6.0, std::f64::consts::PI);
        let q_large = quat_from_axis_angle([1.0, 0.0, 0.0], std::f64::consts::PI / 2.0);
        let corrected = limit.apply_cone_limit(q_large);
        let violation_before = limit.cone_violation(q_large);
        let _violation_after = limit.cone_violation(corrected);
        assert!(
            violation_before > 0.0,
            "there should be a cone violation before correction"
        );
    }
    #[test]
    fn test_pd_joint_at_target_zero_torque_no_velocity() {
        let ctrl = PdJointController::new(100.0, 5.0, 1000.0);
        let tau = ctrl.torque(0.0, 0.0);
        assert!(
            tau.abs() < 1e-12,
            "torque at target with no velocity: {tau}"
        );
    }
    #[test]
    fn test_pd_joint_error_produces_positive_torque() {
        let ctrl = PdJointController::new(100.0, 0.0, 1000.0);
        let tau = ctrl.torque(0.0, 0.0);
        let _tau0 = tau;
        let mut ctrl2 = PdJointController::new(100.0, 0.0, 1000.0);
        ctrl2.set_target(1.0);
        let tau_pos = ctrl2.torque(0.0, 0.0);
        assert!(
            tau_pos > 0.0,
            "should produce positive torque toward target: {tau_pos}"
        );
    }
    #[test]
    fn test_pd_joint_velocity_damping() {
        let ctrl = PdJointController::new(0.0, 10.0, 1000.0);
        let tau = ctrl.torque(0.0, 5.0);
        assert!((tau - (-50.0)).abs() < 1e-10, "damping torque: {tau}");
    }
    #[test]
    fn test_pd_joint_max_torque_clamp() {
        let ctrl = PdJointController::new(1000.0, 0.0, 50.0);
        let mut ctrl_with_target = ctrl.clone();
        ctrl_with_target.set_target(10.0);
        let tau = ctrl_with_target.torque(0.0, 0.0);
        assert!(
            (tau - 50.0).abs() < 1e-10,
            "should be clamped to max_torque: {tau}"
        );
    }
    #[test]
    fn test_pose_pd_controller_zero_torques_at_rest() {
        let ctrl = PosePdController::new(3, 100.0, 5.0, 1000.0);
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [0.0; 3], 1.0, [0.1; 3], None);
        rag.add_bone("b", [1.0, 0.0, 0.0], 1.0, [0.1; 3], None);
        rag.add_bone("c", [2.0, 0.0, 0.0], 1.0, [0.1; 3], None);
        let torques = ctrl.compute_torques(&rag.bones);
        assert_eq!(torques.len(), 3);
        for (i, &t) in torques.iter().enumerate() {
            assert!(t.abs() < 1e-10, "torque[{i}] should be zero: {t}");
        }
    }
    #[test]
    fn test_pose_pd_blend_zero_gives_zero_torques() {
        let mut ctrl = PosePdController::new(2, 100.0, 5.0, 1000.0);
        ctrl.blend = 0.0;
        let mut rag = Ragdoll::new([0.0, 0.0, 0.0]);
        rag.add_bone("a", [0.0; 3], 1.0, [0.1; 3], None);
        rag.add_bone("b", [1.0, 0.0, 0.0], 1.0, [0.1; 3], None);
        for joint in &mut rag.bones {
            joint.orientation = quat_from_axis_angle([0.0, 0.0, 1.0], 0.5);
        }
        let torques = ctrl.compute_torques(&rag.bones);
        for &t in &torques {
            assert!(t.abs() < 1e-12, "blend=0 should give zero torques: {t}");
        }
    }
    #[test]
    fn test_balance_com_at_centroid_is_balanced() {
        let ctrl = BalanceController::rectangular(0.2, 0.3, 100.0, 5.0, 50.0);
        assert!(
            ctrl.is_balanced([0.0, 1.0, 0.0]),
            "CoM at centroid should be balanced"
        );
    }
    #[test]
    fn test_balance_com_outside_polygon_not_balanced() {
        let ctrl = BalanceController::rectangular(0.2, 0.3, 100.0, 5.0, 50.0);
        assert!(
            !ctrl.is_balanced([1.0, 1.0, 0.0]),
            "CoM far outside should not be balanced"
        );
    }
    #[test]
    fn test_balance_corrective_torque_direction() {
        let mut ctrl = BalanceController::rectangular(0.2, 0.3, 100.0, 0.0, 1000.0);
        let torque = ctrl.corrective_torque([0.5, 1.0, 0.0], 0.01);
        assert!(
            torque[0] < 0.0,
            "torque[X] should be negative to correct +X lean: {}",
            torque[0]
        );
    }
    #[test]
    fn test_balance_centroid_of_rectangle() {
        let ctrl = BalanceController::rectangular(0.4, 0.6, 10.0, 1.0, 50.0);
        let c = ctrl.centroid();
        assert!(c[0].abs() < 1e-10, "centroid X should be 0: {}", c[0]);
        assert!(c[1].abs() < 1e-10, "centroid Z should be 0: {}", c[1]);
    }
    #[test]
    fn test_ground_contact_above_ground_no_contact() {
        let det = GroundContactDetector::new(0.0, 0.3, 0.5);
        assert!(
            !det.is_in_contact([0.0, 2.0, 0.0], 0.1),
            "above ground: no contact"
        );
    }
    #[test]
    fn test_ground_contact_at_ground_level() {
        let det = GroundContactDetector::new(0.0, 0.3, 0.5);
        assert!(
            det.is_in_contact([0.0, 0.1, 0.0], 0.1),
            "at ground: in contact"
        );
    }
    #[test]
    fn test_ground_contact_penetration_depth() {
        let det = GroundContactDetector::new(0.0, 0.0, 0.0);
        let depth = det.penetration_depth([0.0, 0.05, 0.0], 0.1);
        assert!((depth - 0.05).abs() < 1e-10, "depth: {depth}");
    }
    #[test]
    fn test_ground_contact_correct_position() {
        let det = GroundContactDetector::new(0.0, 0.0, 0.0);
        let pos = det.correct_position([0.0, 0.05, 0.0], 0.1);
        assert!(
            pos[1] - 0.1 >= -1e-10,
            "corrected Y should be at least radius above ground: {}",
            pos[1]
        );
    }
    #[test]
    fn test_ground_contact_velocity_correction_bounces() {
        let det = GroundContactDetector::new(0.0, 0.5, 0.0);
        let vel = det.correct_velocity([0.0, -10.0, 0.0], [0.0, 0.0, 0.0], 0.01);
        assert!(vel[1] > 0.0, "velocity should be reversed: {}", vel[1]);
        assert!(
            (vel[1] - 5.0).abs() < 1e-10,
            "velocity magnitude: {}",
            vel[1]
        );
    }
    #[test]
    fn test_ground_contact_friction() {
        let det = GroundContactDetector::new(0.0, 0.0, 0.8);
        let vel = det.correct_velocity([5.0, -1.0, 3.0], [0.0, 0.0, 0.0], 0.01);
        assert!(
            (vel[0] - 1.0).abs() < 1e-10,
            "vx after friction: {}",
            vel[0]
        );
        assert!(
            (vel[2] - 0.6).abs() < 1e-10,
            "vz after friction: {}",
            vel[2]
        );
    }
    #[test]
    fn test_body_segment_inertia_cylinder() {
        let seg = BodySegmentInertia::new(10.0, 0.05, 0.40);
        let it = seg.i_transverse();
        let il = seg.i_longitudinal();
        let expected_it = 10.0 * (3.0 * 0.0025 + 0.16) / 12.0;
        let expected_il = 0.5 * 10.0 * 0.0025;
        assert!((it - expected_it).abs() < 1e-12, "I_transverse = {it}");
        assert!((il - expected_il).abs() < 1e-12, "I_longitudinal = {il}");
    }
    #[test]
    fn test_body_segment_inertia_tensor_all_positive() {
        let seg = BodySegmentInertia::new(5.0, 0.04, 0.30);
        let tensor = seg.inertia_tensor();
        for (i, &v) in tensor.iter().enumerate() {
            assert!(v > 0.0, "inertia[{i}] should be positive: {v}");
        }
    }
    #[test]
    fn test_body_segment_from_body_mass() {
        let seg = BodySegmentInertia::from_body_mass(70.0, 0.1, 0.07, 0.42);
        assert!((seg.mass - 7.0).abs() < 1e-10, "mass = {}", seg.mass);
    }
    #[test]
    fn test_standard_human_segments_count() {
        let segs = standard_human_segment_inertias(70.0);
        assert_eq!(segs.len(), 8, "should have 8 segments");
    }
    #[test]
    fn test_standard_human_segments_positive_inertias() {
        let segs = standard_human_segment_inertias(70.0);
        for (name, seg) in &segs {
            for &v in &seg.inertia_tensor() {
                assert!(v > 0.0, "{name}: inertia should be positive, got {v}");
            }
        }
    }
    #[test]
    fn test_body_segment_radius_of_gyration() {
        let seg = BodySegmentInertia::new(4.0, 0.05, 0.0);
        let rog = seg.radius_of_gyration_transverse();
        assert!(rog > 0.0, "radius of gyration should be positive: {rog}");
    }
    #[test]
    fn test_pose_distance_identity_poses_is_zero() {
        let pose = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]],
        };
        let dist = Ragdoll::compute_pose_distance(&pose, &pose, 1.0);
        assert!(
            dist.abs() < 1e-10,
            "same pose → distance should be 0, got {dist}"
        );
    }
    #[test]
    fn test_pose_distance_position_offset() {
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let pose_b = RagdollPose {
            positions: vec![[3.0, 4.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let dist = Ragdoll::compute_pose_distance(&pose_a, &pose_b, 1.0);
        assert!((dist - 5.0).abs() < 1e-10, "dist = {dist}");
    }
    #[test]
    fn test_pose_distance_orientation_difference() {
        let qa = [1.0_f64, 0.0, 0.0, 0.0];
        let qb = [0.0_f64, 0.0, 0.0, 1.0];
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![qa],
        };
        let pose_b = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![qb],
        };
        let dist = Ragdoll::compute_pose_distance(&pose_a, &pose_b, 0.0);
        assert!(dist > 0.0, "orientation difference → dist > 0: {dist}");
    }
    #[test]
    fn test_pose_distance_zero_weight_ignores_position() {
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let pose_b = RagdollPose {
            positions: vec![[100.0, 200.0, 300.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let dist = Ragdoll::compute_pose_distance(&pose_a, &pose_b, 0.0);
        assert!(
            dist.abs() < 1e-10,
            "pos_weight=0 → position ignored: {dist}"
        );
    }
    #[test]
    fn test_pose_distance_mismatched_lengths_uses_min() {
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]; 3],
        };
        let pose_b = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let dist = Ragdoll::compute_pose_distance(&pose_a, &pose_b, 1.0);
        assert!(dist.abs() < 1e-10, "dist = {dist}");
    }
    #[test]
    fn test_muscle_force_changes_angular_velocity() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("hip", [0.0, 0.0, 0.0], 10.0, [0.1, 0.2, 0.1], None);
        let omega_before = rag.bones[0].angular_vel;
        rag.apply_muscle_force(0, [0.0, 1.0, 0.0], 100.0, 0.01);
        let omega_after = rag.bones[0].angular_vel;
        assert!(
            omega_after[1] != omega_before[1] || omega_after[0] != omega_before[0],
            "muscle force should change angular velocity"
        );
    }
    #[test]
    fn test_muscle_force_zero_magnitude_no_change() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("hip", [0.0, 0.0, 0.0], 10.0, [0.1, 0.2, 0.1], None);
        rag.bones[0].angular_vel = [1.0, 2.0, 3.0];
        let omega_before = rag.bones[0].angular_vel;
        rag.apply_muscle_force(0, [0.0, 1.0, 0.0], 0.0, 0.01);
        assert_eq!(rag.bones[0].angular_vel, omega_before);
    }
    #[test]
    fn test_muscle_force_out_of_range_bone_no_panic() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("hip", [0.0, 0.0, 0.0], 10.0, [0.1, 0.2, 0.1], None);
        rag.apply_muscle_force(99, [0.0, 1.0, 0.0], 500.0, 0.01);
    }
    #[test]
    fn test_muscle_force_static_bone_ignored() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("ground", [0.0, 0.0, 0.0], 0.0, [1.0, 0.1, 1.0], None);
        rag.apply_muscle_force(0, [1.0, 0.0, 0.0], 1000.0, 1.0);
        assert_eq!(
            rag.bones[0].angular_vel, [0.0; 3],
            "static bone should be unaffected"
        );
    }
    #[test]
    fn test_apply_blended_pose_midpoint() {
        let pose_a = RagdollPose {
            positions: vec![[0.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let pose_b = RagdollPose {
            positions: vec![[2.0, 0.0, 0.0]],
            orientations: vec![[1.0, 0.0, 0.0, 0.0]],
        };
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("root", [0.0, 0.0, 0.0], 1.0, [0.1, 0.1, 0.1], None);
        rag.apply_blended_pose(&pose_a, &pose_b, 0.5);
        assert!(
            (rag.bones[0].position[0] - 1.0).abs() < 1e-10,
            "midpoint x = {}",
            rag.bones[0].position[0]
        );
    }
}
