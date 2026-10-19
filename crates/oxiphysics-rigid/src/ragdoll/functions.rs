//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ArticulatedPose, BodySegmentInertia, BoneDescriptor, BoneTransform};

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
/// Angular velocity (world frame) that rotates orientation `from` into `to`
/// within `dt` seconds along the shortest arc.
///
/// This is the logarithm map of `to ⊗ from⁻¹` divided by `dt`; a rotation by
/// `ω·dt` about the axis of `ω` (as applied by [`crate::Ragdoll::integrate`])
/// is recovered exactly.
pub(super) fn quat_delta_angular_velocity(from: [f64; 4], to: [f64; 4], dt: f64) -> [f64; 3] {
    let mut delta = quat_mul(to, quat_conjugate(from));
    if delta[0] < 0.0 {
        delta = [-delta[0], -delta[1], -delta[2], -delta[3]];
    }
    let axis_part = [delta[1], delta[2], delta[3]];
    let sin_half = vec_length(axis_part);
    if sin_half < 1e-12 {
        // Small-angle limit of angle / sin(angle / 2).
        return vec_scale(axis_part, 2.0 / dt);
    }
    let angle = 2.0 * sin_half.atan2(delta[0]);
    vec_scale(axis_part, angle / (sin_half * dt))
}
/// Quaternion (w, x, y, z) for Euler angles `[ax, ay, az]` in radians.
///
/// The rotations are composed in intrinsic X-Y-Z order, `q = qx(ax) ⊗ qy(ay) ⊗
/// qz(az)`, which as a rotation matrix is `Rx(ax) · Ry(ay) · Rz(az)`: a vector
/// is rotated about Z first, then Y, then X (all about fixed axes).
pub fn euler_xyz_to_quat(angles: [f64; 3]) -> [f64; 4] {
    let qx = quat_from_axis_angle([1.0, 0.0, 0.0], angles[0]);
    let qy = quat_from_axis_angle([0.0, 1.0, 0.0], angles[1]);
    let qz = quat_from_axis_angle([0.0, 0.0, 1.0], angles[2]);
    quat_normalize(quat_mul(quat_mul(qx, qy), qz))
}
/// Compute the world-space frame of every bone using forward kinematics.
///
/// Every bone extends along its local +Y axis from its proximal joint, and a
/// child attaches at the distal end of its parent.  A root bone sits at
/// `pose.root_position` with orientation
/// `normalize(pose.root_orientation) ⊗ euler(joint_angles[root])`; a child's
/// frame is
///
/// ```text
/// position_i    = position_parent + R_parent · [0, length_parent, 0]
/// orientation_i = orientation_parent ⊗ euler(joint_angles[i])
/// ```
///
/// with [`euler_xyz_to_quat`] for the joint rotation.  Bones may be listed in
/// any order (a parent may come after its child).  A bone whose parent index is
/// out of range, or that lies on a parent cycle, cannot be attached and is
/// treated as an additional root placed at `pose.root_position`.
pub fn forward_kinematics_transforms(
    bones: &[BoneDescriptor],
    pose: &ArticulatedPose,
) -> Vec<BoneTransform> {
    let n = bones.len();
    let root_orientation = quat_normalize(pose.root_orientation);
    let joint_rotation = |idx: usize| -> [f64; 4] {
        pose.joint_angles
            .get(idx)
            .map_or([1.0, 0.0, 0.0, 0.0], |&angles| euler_xyz_to_quat(angles))
    };
    let mut resolved: Vec<Option<BoneTransform>> = vec![None; n];
    for start in 0..n {
        if resolved[start].is_some() {
            continue;
        }
        // Climb through unresolved ancestors until a resolved bone, a root or
        // an invalid link ends the chain.
        let mut chain = vec![start];
        let mut anchor: Option<BoneTransform> = None;
        let mut cursor = start;
        while let Some(parent) = bones[cursor].parent {
            if parent >= n {
                break;
            }
            if let Some(done) = resolved[parent] {
                anchor = Some(done);
                break;
            }
            if chain.contains(&parent) {
                break;
            }
            chain.push(parent);
            cursor = parent;
        }
        // Place the chain from its top-most bone down to `start`.
        let mut above = anchor;
        for &idx in chain.iter().rev() {
            let local = joint_rotation(idx);
            let transform = match (above, bones[idx].parent) {
                (Some(parent_frame), Some(parent_idx)) => {
                    let along_parent = quat_rotate(
                        parent_frame.orientation,
                        [0.0, bones[parent_idx].length, 0.0],
                    );
                    BoneTransform {
                        position: vec_add(parent_frame.position, along_parent),
                        orientation: quat_normalize(quat_mul(parent_frame.orientation, local)),
                    }
                }
                _ => BoneTransform {
                    position: pose.root_position,
                    orientation: quat_normalize(quat_mul(root_orientation, local)),
                },
            };
            resolved[idx] = Some(transform);
            above = Some(transform);
        }
    }
    resolved
        .into_iter()
        .map(|frame| {
            frame.unwrap_or(BoneTransform {
                position: pose.root_position,
                orientation: root_orientation,
            })
        })
        .collect()
}
/// Compute world-space positions for each bone using forward kinematics.
///
/// Returns the proximal joint position of every bone as computed by
/// [`forward_kinematics_transforms`]: the root bone sits at
/// `pose.root_position`, and each child sits at the distal end of its parent,
/// `length_parent` along the parent's rotated +Y axis.  The root orientation
/// and the joint angles of `pose` rotate the whole chain; with all angles zero
/// and the identity root orientation the bones stack along +Y.
pub fn forward_kinematics(bones: &[BoneDescriptor], pose: &ArticulatedPose) -> Vec<[f64; 3]> {
    forward_kinematics_transforms(bones, pose)
        .into_iter()
        .map(|frame| frame.position)
        .collect()
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
        HumanoidRagdoll, HumanoidSkeleton, JointLimit, Ragdoll, RagdollPose, RagdollState,
        SkeletonJoint,
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
    fn test_bone(parent: Option<usize>, length: f64, mass: f64) -> BoneDescriptor {
        BoneDescriptor {
            name: "bone".to_string(),
            parent,
            length,
            mass,
            inertia: [0.1, 0.1, 0.1],
        }
    }
    fn assert_close(a: [f64; 3], b: [f64; 3], tol: f64, what: &str) {
        let diff = vec_length(vec_sub(a, b));
        assert!(diff < tol, "{what}: {a:?} vs {b:?}");
    }
    /// Largest component difference between two unit quaternions, modulo sign.
    fn quat_difference(a: [f64; 4], b: [f64; 4]) -> f64 {
        let direct = (0..4).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max);
        let flipped = (0..4).map(|k| (a[k] + b[k]).abs()).fold(0.0, f64::max);
        direct.min(flipped)
    }
    #[test]
    fn test_fk_root_rotation_about_z_moves_child_to_negative_x() {
        let bones = vec![test_bone(None, 1.0, 1.0), test_bone(Some(0), 0.5, 1.0)];
        let mut pose = ArticulatedPose::default_t_pose(2);
        pose.root_orientation = quat_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let positions = forward_kinematics(&bones, &pose);
        assert_close(positions[0], [0.0; 3], 1e-12, "root");
        assert_close(positions[1], [-1.0, 0.0, 0.0], 1e-12, "child");
        // The root position is only an offset.
        pose.root_position = [1.0, 2.0, 3.0];
        let shifted = forward_kinematics(&bones, &pose);
        assert_close(shifted[1], [0.0, 2.0, 3.0], 1e-12, "shifted child");
    }
    #[test]
    fn test_fk_normalises_the_root_orientation() {
        let bones = vec![test_bone(None, 2.0, 1.0), test_bone(Some(0), 1.0, 1.0)];
        let mut pose = ArticulatedPose::default_t_pose(2);
        pose.root_orientation = [3.0, 0.0, 0.0, 0.0];
        assert_close(
            forward_kinematics(&bones, &pose)[1],
            [0.0, 2.0, 0.0],
            1e-12,
            "scaled identity",
        );
        // A scaled quaternion for the same 90 degree turn about Z.
        let q = quat_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        pose.root_orientation = [4.0 * q[0], 4.0 * q[1], 4.0 * q[2], 4.0 * q[3]];
        assert_close(
            forward_kinematics(&bones, &pose)[1],
            [-2.0, 0.0, 0.0],
            1e-12,
            "scaled turn",
        );
    }
    #[test]
    fn test_fk_root_joint_angles_compose_with_the_root_orientation() {
        let bones = vec![test_bone(None, 1.0, 1.0), test_bone(Some(0), 1.0, 1.0)];
        let quarter = std::f64::consts::FRAC_PI_2;
        let mut by_joint = ArticulatedPose::default_t_pose(2);
        by_joint.joint_angles[0] = [0.0, 0.0, quarter];
        let mut by_root = ArticulatedPose::default_t_pose(2);
        by_root.root_orientation = quat_from_axis_angle([0.0, 0.0, 1.0], quarter);
        assert_close(
            forward_kinematics(&bones, &by_joint)[1],
            forward_kinematics(&bones, &by_root)[1],
            1e-12,
            "joint angle vs root orientation",
        );
        // Both together add up to a half turn.
        by_joint.root_orientation = quat_from_axis_angle([0.0, 0.0, 1.0], quarter);
        assert_close(
            forward_kinematics(&bones, &by_joint)[1],
            [0.0, -1.0, 0.0],
            1e-12,
            "half turn",
        );
    }
    #[test]
    fn test_fk_joint_angles_bend_the_chain() {
        let bones = vec![
            test_bone(None, 1.0, 1.0),
            test_bone(Some(0), 2.0, 1.0),
            test_bone(Some(1), 0.5, 1.0),
        ];
        let mut pose = ArticulatedPose::default_t_pose(3);
        // Bend the middle bone 90 degrees about X: its +Y axis becomes +Z.
        pose.joint_angles[1] = [std::f64::consts::FRAC_PI_2, 0.0, 0.0];
        let frames = forward_kinematics_transforms(&bones, &pose);
        assert_close(frames[1].position, [0.0, 1.0, 0.0], 1e-12, "middle joint");
        assert_close(frames[2].position, [0.0, 1.0, 2.0], 1e-12, "tip joint");
        // The tip inherits the bend: its local +Y points along world +Z.
        assert_close(
            quat_rotate(frames[2].orientation, [0.0, 1.0, 0.0]),
            [0.0, 0.0, 1.0],
            1e-12,
            "tip axis",
        );
        // An additional twist about the bone's own axis changes nothing but the frame.
        pose.joint_angles[1] = [std::f64::consts::FRAC_PI_2, 0.0, 0.0];
        pose.joint_angles[2] = [0.0, 0.0, 0.3];
        let twisted = forward_kinematics(&bones, &pose);
        assert_close(
            twisted[2],
            [0.0, 1.0, 2.0],
            1e-12,
            "tip joint after distal twist",
        );
    }
    #[test]
    fn test_euler_xyz_order_is_intrinsic_x_then_y_then_z() {
        let quarter = std::f64::consts::FRAC_PI_2;
        let q = euler_xyz_to_quat([quarter, quarter, 0.0]);
        // R = Rx(90) * Ry(90): x -> -z (about Y) -> +y (about X).
        assert_close(
            quat_rotate(q, [1.0, 0.0, 0.0]),
            [0.0, 1.0, 0.0],
            1e-12,
            "Rx*Ry applied to x",
        );
        let q_reversed = euler_xyz_to_quat([0.0, quarter, quarter]);
        // R = Ry(90) * Rz(90): x -> y (about Z) -> y (about Y).
        assert_close(
            quat_rotate(q_reversed, [1.0, 0.0, 0.0]),
            [0.0, 1.0, 0.0],
            1e-12,
            "Ry*Rz applied to x",
        );
        assert_close(
            quat_rotate(q_reversed, [0.0, 0.0, 1.0]),
            [1.0, 0.0, 0.0],
            1e-12,
            "Ry*Rz applied to z",
        );
        let identity = euler_xyz_to_quat([0.0; 3]);
        assert_eq!(identity, [1.0, 0.0, 0.0, 0.0]);
        let q_all = euler_xyz_to_quat([0.3, -0.7, 1.1]);
        let len2: f64 = q_all.iter().map(|c| c * c).sum();
        assert!((len2 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fk_accepts_children_listed_before_their_parents() {
        // Same skeleton as the root/child/grandchild chain, stored in reverse.
        let reversed = vec![
            test_bone(Some(1), 0.5, 1.0),
            test_bone(Some(2), 2.0, 1.0),
            test_bone(None, 1.0, 1.0),
        ];
        let forward = vec![
            test_bone(None, 1.0, 1.0),
            test_bone(Some(0), 2.0, 1.0),
            test_bone(Some(1), 0.5, 1.0),
        ];
        let mut pose_forward = ArticulatedPose::default_t_pose(3);
        pose_forward.root_position = [1.0, 2.0, 3.0];
        pose_forward.root_orientation = quat_from_axis_angle([0.0, 0.0, 1.0], 0.9);
        pose_forward.joint_angles[1] = [0.4, 0.0, 0.0];
        let mut pose_reversed = pose_forward.clone();
        pose_reversed.joint_angles = vec![[0.0; 3], [0.4, 0.0, 0.0], [0.0; 3]];
        let expected = forward_kinematics(&forward, &pose_forward);
        let got = forward_kinematics(&reversed, &pose_reversed);
        assert_close(got[2], expected[0], 1e-12, "root");
        assert_close(got[1], expected[1], 1e-12, "middle");
        assert_close(got[0], expected[2], 1e-12, "tip");
        assert_close(got[2], [1.0, 2.0, 3.0], 1e-12, "root position");
    }
    #[test]
    fn test_fk_detached_bones_are_placed_at_the_root_without_panicking() {
        let bones = vec![
            test_bone(None, 1.0, 1.0),
            test_bone(Some(99), 1.0, 1.0),
            test_bone(Some(2), 1.0, 1.0),
            test_bone(Some(3), 1.0, 1.0),
            test_bone(Some(2), 1.0, 1.0),
            test_bone(Some(5), 1.0, 1.0),
        ];
        let mut pose = ArticulatedPose::default_t_pose(6);
        pose.root_position = [0.5, 0.0, 0.0];
        let first = forward_kinematics(&bones, &pose);
        let second = forward_kinematics(&bones, &pose);
        assert_eq!(first, second, "resolution is deterministic");
        assert_eq!(first.len(), 6);
        assert!(first.iter().all(|p| p.iter().all(|c| c.is_finite())));
        // Dangling parent index: an extra root.
        assert_close(first[1], [0.5, 0.0, 0.0], 1e-12, "dangling parent");
        // The 2 <-> 3 cycle is broken once; its members and the bone hanging
        // off it stay finite and attached to each other.
        assert_close(
            vec_sub(first[4], first[2]),
            [0.0, 1.0, 0.0],
            1e-12,
            "child of a cycle member",
        );
        // A bone that is its own parent is a root.
        let looped = vec![test_bone(Some(0), 1.0, 1.0)];
        let single = forward_kinematics(&looped, &ArticulatedPose::default_t_pose(1));
        assert_close(single[0], [0.0; 3], 1e-12, "self parent");
    }
    #[test]
    fn test_fk_missing_joint_angles_default_to_zero() {
        let bones = vec![test_bone(None, 1.0, 1.0), test_bone(Some(0), 1.0, 1.0)];
        let mut pose = ArticulatedPose::default_t_pose(0);
        pose.root_position = [0.0, 1.0, 0.0];
        assert_close(
            forward_kinematics(&bones, &pose)[1],
            [0.0, 2.0, 0.0],
            1e-12,
            "no joint angles",
        );
        assert!(forward_kinematics(&[], &pose).is_empty());
    }
    #[test]
    fn test_humanoid_center_of_mass_follows_root_orientation_and_joint_angles() {
        let skeleton = HumanoidSkeleton {
            bones: vec![test_bone(None, 1.0, 2.0), test_bone(Some(0), 1.0, 1.0)],
            joints: Vec::new(),
        };
        // Bone midpoints at y = 0.5 and 1.5.
        let rest = ArticulatedPose::default_t_pose(2);
        assert_close(
            skeleton.center_of_mass(&rest),
            [0.0, 2.5 / 3.0, 0.0],
            1e-12,
            "rest",
        );
        // A 90 degree root turn about Z lays the stack along -X.
        let mut turned = ArticulatedPose::default_t_pose(2);
        turned.root_orientation =
            quat_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        assert_close(
            skeleton.center_of_mass(&turned),
            [-2.5 / 3.0, 0.0, 0.0],
            1e-12,
            "turned root",
        );
        // Folding the second bone 90 degrees about X moves its midpoint to (0, 1, 0.5).
        let mut bent = ArticulatedPose::default_t_pose(2);
        bent.joint_angles[1] = [std::f64::consts::FRAC_PI_2, 0.0, 0.0];
        assert_close(
            skeleton.center_of_mass(&bent),
            [0.0, (2.0 * 0.5 + 1.0) / 3.0, 0.5 / 3.0],
            1e-12,
            "bent joint",
        );
        // A standard human rotates rigidly with its root.
        let human = HumanoidRagdoll::standard_human();
        let mut pose = ArticulatedPose::default_t_pose(human.bones.len());
        let upright = human.center_of_mass(&pose);
        assert!(upright[1] > 0.0 && upright[0].abs() < 1e-12);
        pose.root_orientation = quat_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let sideways = human.center_of_mass(&pose);
        assert_close(
            sideways,
            [-upright[1], 0.0, 0.0],
            1e-9,
            "human centre of mass",
        );
    }
    fn hanging_chain(bones_below: usize) -> Ragdoll {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        // Static anchor bone (zero mass) at the origin.
        rag.add_bone("anchor", [0.0, 0.0, 0.0], 0.0, [0.05, 0.05, 0.05], None);
        for k in 0..bones_below {
            let centre_y = -0.5 - k as f64;
            rag.add_bone(
                "link",
                [0.0, centre_y, 0.0],
                1.0,
                [0.05, 0.5, 0.05],
                Some(k),
            );
            // Top of the new link meets the bottom of the previous one.
            let anchor_on_previous = if k == 0 { [0.0; 3] } else { [0.0, -0.5, 0.0] };
            rag.add_ball_joint(
                k,
                k + 1,
                anchor_on_previous,
                [0.0, 0.5, 0.0],
                std::f64::consts::PI,
            );
        }
        rag
    }
    fn max_anchor_error(rag: &Ragdoll) -> f64 {
        rag.ball_joints
            .iter()
            .map(|j| {
                let a = rag.bone_world_anchor(j.bone_a, j.local_anchor_a);
                let b = rag.bone_world_anchor(j.bone_b, j.local_anchor_b);
                vec_length(vec_sub(a, b))
            })
            .fold(0.0, f64::max)
    }
    #[test]
    fn test_hanging_chain_from_a_static_root_comes_to_rest() {
        let mut rag = hanging_chain(3);
        let dt = 1.0 / 60.0;
        let mut worst_speed = 0.0_f64;
        for step in 0..1000 {
            rag.step(dt, 20);
            if step >= 500 {
                for bone in &rag.bones {
                    worst_speed = worst_speed.max(vec_length(bone.linear_vel));
                }
            }
        }
        assert!(
            worst_speed < 1e-3,
            "hanging chain keeps moving at {worst_speed} m/s"
        );
        let anchor_error = max_anchor_error(&rag);
        assert!(anchor_error < 1e-4, "joint anchor error {anchor_error}");
        // Without the velocity update the speed grew as g*t, i.e. ~163 m/s here.
        assert_eq!(
            rag.bones[0].linear_vel, [0.0; 3],
            "static bones keep zero velocity"
        );
        for bone in rag.bones.iter().skip(1) {
            assert!(bone.angular_vel.iter().all(|w| w.abs() < 1e-9));
        }
    }
    #[test]
    fn test_free_fall_matches_semi_implicit_euler() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        rag.add_bone("ball", [0.0, 100.0, 0.0], 2.0, [0.1, 0.1, 0.1], None);
        let dt = 0.01;
        let steps = 100;
        for _ in 0..steps {
            rag.step(dt, 4);
        }
        let t = dt * steps as f64;
        assert!((rag.bones[0].linear_vel[1] + 9.81 * t).abs() < 1e-9);
        let expected_y = 100.0 - 9.81 * dt * dt * (steps * (steps + 1)) as f64 / 2.0;
        assert!((rag.bones[0].position[1] - expected_y).abs() < 1e-9);
    }
    #[test]
    fn test_unconstrained_bone_keeps_linear_and_angular_velocity() {
        let mut rag = Ragdoll::new([0.0; 3]);
        rag.add_bone("spinner", [1.0, 2.0, 3.0], 1.0, [0.1, 0.2, 0.3], None);
        rag.bones[0].linear_vel = [1.0, -2.0, 0.5];
        rag.bones[0].angular_vel = [0.5, -1.5, 2.0];
        let start_orientation = rag.bones[0].orientation;
        let dt = 1.0 / 120.0;
        for _ in 0..240 {
            rag.step(dt, 3);
        }
        assert_close(
            rag.bones[0].linear_vel,
            [1.0, -2.0, 0.5],
            1e-12,
            "linear velocity",
        );
        assert_close(
            rag.bones[0].angular_vel,
            [0.5, -1.5, 2.0],
            1e-9,
            "angular velocity",
        );
        assert_close(
            rag.bones[0].position,
            [3.0, -2.0, 4.0],
            1e-9,
            "position after 2 s",
        );
        // The orientation advanced by |omega| * t about the fixed world axis omega.
        let omega = [0.5, -1.5, 2.0];
        let speed = vec_length(omega);
        let expected = quat_normalize(quat_mul(
            quat_from_axis_angle(omega, speed * 2.0),
            start_orientation,
        ));
        assert!(quat_difference(rag.bones[0].orientation, expected) < 1e-9);
    }
    #[test]
    fn test_step_turns_constraint_corrections_into_velocity() {
        let mut rag = Ragdoll::new([0.0; 3]);
        let pivot = rag.add_bone("pivot", [0.0; 3], 0.0, [0.1; 3], None);
        let bob = rag.add_bone("bob", [0.0, 0.0, 0.0], 1.0, [0.1; 3], Some(pivot));
        rag.add_ball_joint(pivot, bob, [0.0; 3], [0.0; 3], std::f64::consts::PI);
        // Displace the bob; the joint pulls it back by the full offset in one pass.
        rag.bones[bob].position = [0.3, -0.4, 0.0];
        let dt = 0.1;
        rag.step(dt, 1);
        assert_close(
            rag.bones[bob].position,
            [0.0; 3],
            1e-12,
            "bob is projected onto the pivot",
        );
        assert_close(
            rag.bones[bob].linear_vel,
            [-3.0, 4.0, 0.0],
            1e-9,
            "velocity of the correction",
        );
        assert_eq!(rag.bones[pivot].linear_vel, [0.0; 3]);
    }
    #[test]
    fn test_step_derives_angular_velocity_from_orientation_projection() {
        let mut rag = Ragdoll::new([0.0; 3]);
        let parent = rag.add_bone("parent", [0.0; 3], 0.0, [0.1; 3], None);
        let child = rag.add_bone("child", [0.0; 3], 1.0, [0.1; 3], Some(parent));
        rag.add_ball_joint(parent, child, [0.0; 3], [0.0; 3], 0.5);
        let before = quat_from_axis_angle([0.0, 0.0, 1.0], 1.2);
        rag.bones[child].orientation = before;
        let dt = 1.0 / 60.0;
        rag.step(dt, 1);
        let after = rag.bones[child].orientation;
        assert!(
            quat_angle(after, before) > 0.1,
            "the swing limit must pull the child back"
        );
        let omega = rag.bones[child].angular_vel;
        assert!(
            omega[2] < 0.0,
            "the correction rotates the child back towards the parent: {omega:?}"
        );
        // omega * dt reproduces the rotation that took `before` to `after`.
        let replay = quat_normalize(quat_mul(
            quat_from_axis_angle(omega, vec_length(omega) * dt),
            before,
        ));
        assert!(
            quat_difference(replay, after) < 1e-9,
            "replayed orientation differs"
        );
    }
    #[test]
    fn test_inactive_ragdoll_stays_at_zero_velocity_while_constraints_run() {
        let mut rag = Ragdoll::new([0.0, -9.81, 0.0]);
        let pivot = rag.add_bone("pivot", [0.0; 3], 0.0, [0.1; 3], None);
        let bob = rag.add_bone("bob", [0.0; 3], 1.0, [0.1; 3], Some(pivot));
        rag.add_ball_joint(pivot, bob, [0.0; 3], [0.0; 3], std::f64::consts::PI);
        rag.deactivate();
        rag.bones[bob].position = [0.3, -0.4, 0.0];
        rag.step(0.1, 1);
        assert_eq!(rag.bones[bob].linear_vel, [0.0; 3]);
        assert_eq!(rag.bones[bob].angular_vel, [0.0; 3]);
    }
    #[test]
    fn test_step_with_zero_dt_does_not_produce_non_finite_velocities() {
        let mut rag = hanging_chain(2);
        rag.bones[1].linear_vel = [0.2, 0.0, -0.1];
        rag.step(0.0, 4);
        for bone in &rag.bones {
            assert!(
                bone.linear_vel
                    .iter()
                    .chain(&bone.angular_vel)
                    .all(|c| c.is_finite())
            );
        }
        assert_eq!(rag.bones[1].linear_vel, [0.2, 0.0, -0.1]);
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
