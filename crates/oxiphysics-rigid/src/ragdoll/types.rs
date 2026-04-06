//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

/// Skeleton joint definition for building ragdolls from skeletons.
#[derive(Debug, Clone)]
pub struct SkeletonJoint {
    /// Name of the joint/bone.
    pub name: String,
    /// Position in bind pose.
    pub bind_position: [f64; 3],
    /// Orientation in bind pose (w, x, y, z).
    pub bind_orientation: [f64; 4],
    /// Parent joint index, if any.
    pub parent: Option<usize>,
    /// Mass estimate for this body part.
    pub mass: f64,
    /// Half-extents for collision shape.
    pub half_extents: [f64; 3],
}
/// Activation state of a ragdoll.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RagdollState {
    /// Fully simulated.
    Active,
    /// Not simulated (kinematic / animation-driven).
    Inactive,
    /// Blending between animation and simulation.
    Blending {
        /// Blend factor 0..1 (0 = animation, 1 = physics).
        blend_factor: f64,
    },
}
/// A pose snapshot: position and orientation for each bone.
#[derive(Debug, Clone)]
pub struct RagdollPose {
    /// Position per bone.
    pub positions: Vec<[f64; 3]>,
    /// Orientation per bone.
    pub orientations: Vec<[f64; 4]>,
}
/// Balance controller that checks whether the center of mass (CoM) lies above
/// the support polygon and computes a corrective lean torque when needed.
///
/// Assumes a flat ground contact polygon defined by foot positions.
#[derive(Debug, Clone)]
pub struct BalanceController {
    /// Support polygon vertices (XZ plane, Y = ground level).
    pub support_polygon: Vec<[f64; 2]>,
    /// Proportional gain for lean correction torque.
    pub kp: f64,
    /// Derivative gain for lean correction.
    pub kd: f64,
    /// Maximum corrective torque (N·m).
    pub max_torque: f64,
    /// Previous CoM XZ error for derivative term.
    pub(super) prev_error: [f64; 2],
}
impl BalanceController {
    /// Create a new balance controller with a rectangular support polygon.
    ///
    /// `width` is left-right (X), `depth` is front-back (Z).
    pub fn rectangular(width: f64, depth: f64, kp: f64, kd: f64, max_torque: f64) -> Self {
        let hw = width * 0.5;
        let hd = depth * 0.5;
        Self {
            support_polygon: vec![[-hw, -hd], [hw, -hd], [hw, hd], [-hw, hd]],
            kp,
            kd,
            max_torque,
            prev_error: [0.0; 2],
        }
    }
    /// Check whether the XZ point `p` lies inside the support polygon using
    /// the ray-casting algorithm.
    pub fn point_in_polygon(&self, p: [f64; 2]) -> bool {
        let n = self.support_polygon.len();
        if n < 3 {
            return false;
        }
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let vi = self.support_polygon[i];
            let vj = self.support_polygon[j];
            if ((vi[1] > p[1]) != (vj[1] > p[1]))
                && (p[0] < (vj[0] - vi[0]) * (p[1] - vi[1]) / (vj[1] - vi[1]) + vi[0])
            {
                inside = !inside;
            }
            j = i;
        }
        inside
    }
    /// Compute the centroid of the support polygon.
    pub fn centroid(&self) -> [f64; 2] {
        let n = self.support_polygon.len() as f64;
        if n < 1.0 {
            return [0.0; 2];
        }
        let sum = self
            .support_polygon
            .iter()
            .fold([0.0f64; 2], |acc, v| [acc[0] + v[0], acc[1] + v[1]]);
        [sum[0] / n, sum[1] / n]
    }
    /// Compute the CoM error (XZ displacement from centroid).
    pub fn com_error(&self, com: [f64; 3]) -> [f64; 2] {
        let c = self.centroid();
        [com[0] - c[0], com[2] - c[1]]
    }
    /// Compute a corrective torque vector `[tx, 0, tz]` to restore balance.
    ///
    /// * `com` — current center of mass in world space.
    /// * `dt`  — time step.
    pub fn corrective_torque(&mut self, com: [f64; 3], dt: f64) -> [f64; 3] {
        let err = self.com_error(com);
        let d_err_x = if dt > 1e-15 {
            (err[0] - self.prev_error[0]) / dt
        } else {
            0.0
        };
        let d_err_z = if dt > 1e-15 {
            (err[1] - self.prev_error[1]) / dt
        } else {
            0.0
        };
        self.prev_error = err;
        let tx = -(self.kp * err[0] + self.kd * d_err_x).clamp(-self.max_torque, self.max_torque);
        let tz = -(self.kp * err[1] + self.kd * d_err_z).clamp(-self.max_torque, self.max_torque);
        [tx, 0.0, tz]
    }
    /// Whether the ragdoll is currently balanced (CoM inside support polygon).
    pub fn is_balanced(&self, com: [f64; 3]) -> bool {
        self.point_in_polygon([com[0], com[2]])
    }
}
/// Hinge (revolute) joint between two bones.
#[derive(Debug, Clone)]
pub struct HingeJoint {
    /// Index of the first bone.
    pub bone_a: usize,
    /// Index of the second bone.
    pub bone_b: usize,
    /// Anchor in bone A local space.
    pub local_anchor_a: [f64; 3],
    /// Anchor in bone B local space.
    pub local_anchor_b: [f64; 3],
    /// Hinge axis expressed in bone A local space.
    pub local_axis_a: [f64; 3],
    /// Minimum allowed angle (radians).
    pub min_angle: f64,
    /// Maximum allowed angle (radians).
    pub max_angle: f64,
}
/// Descriptor for a single bone in an articulated skeleton.
pub struct BoneDescriptor {
    /// Human-readable name.
    pub name: String,
    /// Index of the parent bone, if any.
    pub parent: Option<usize>,
    /// Length of the bone (m).
    pub length: f64,
    /// Mass of the bone segment (kg).
    pub mass: f64,
    /// Principal-axis inertia of the bone (kg·m²).
    pub inertia: [f64; 3],
}
/// Multi-joint PD controller that drives all bones in a ragdoll toward a
/// reference pose captured from an animation or pre-defined stance.
#[derive(Debug, Clone)]
pub struct PosePdController {
    /// One PD controller per bone.
    pub joints: Vec<PdJointController>,
    /// Blending weight: 0 = pure physics, 1 = full pose-following.
    pub blend: f64,
}
impl PosePdController {
    /// Create controllers for `n_bones` joints with uniform gains.
    pub fn new(n_bones: usize, kp: f64, kd: f64, max_torque: f64) -> Self {
        Self {
            joints: (0..n_bones)
                .map(|_| PdJointController::new(kp, kd, max_torque))
                .collect(),
            blend: 1.0,
        }
    }
    /// Set the target pose (reference quaternion per joint).
    pub fn set_pose(&mut self, pose: &RagdollPose) {
        for (i, joint) in self.joints.iter_mut().enumerate() {
            if i < pose.orientations.len() {
                let q = pose.orientations[i];
                let angle = 2.0
                    * (q[0] * q[3] + q[1] * q[2]).atan2(1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]));
                joint.set_target(angle);
            }
        }
    }
    /// Compute the torques for all joints given current bone orientations and
    /// angular velocities.
    ///
    /// Returns a vec of torques (N·m per joint).
    pub fn compute_torques(&self, bones: &[RagdollBone]) -> Vec<f64> {
        self.joints
            .iter()
            .zip(bones.iter())
            .map(|(ctrl, bone)| {
                let q = bone.orientation;
                let angle = 2.0
                    * (q[0] * q[3] + q[1] * q[2]).atan2(1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]));
                let ang_vel = bone.angular_vel[2];
                ctrl.torque(angle, ang_vel) * self.blend
            })
            .collect()
    }
}
/// An articulated pose described by joint angles and a root transform.
pub struct ArticulatedPose {
    /// Joint angles for each bone (Euler angles in radians).
    pub joint_angles: Vec<[f64; 3]>,
    /// Root position in world space.
    pub root_position: [f64; 3],
    /// Root orientation as quaternion (w, x, y, z).
    pub root_orientation: [f64; 4],
}
impl ArticulatedPose {
    /// Create a default T-pose (all joint angles zero).
    pub fn default_t_pose(n_bones: usize) -> Self {
        Self {
            joint_angles: vec![[0.0; 3]; n_bones],
            root_position: [0.0; 3],
            root_orientation: [1.0, 0.0, 0.0, 0.0],
        }
    }
}
/// Builder for a standard 15-bone humanoid ragdoll.
pub struct HumanoidRagdoll;
impl HumanoidRagdoll {
    /// Build a standard 15-bone humanoid ragdoll (Y-up coordinate system).
    #[allow(clippy::new_ret_no_self)]
    pub fn new() -> Ragdoll {
        let gravity = [0.0, -9.81, 0.0];
        let mut rag = Ragdoll::new(gravity);
        let pi = std::f64::consts::PI;
        let pelvis = rag.add_bone("pelvis", [0.0, 1.00, 0.0], 8.0, [0.12, 0.10, 0.08], None);
        let spine = rag.add_bone(
            "spine",
            [0.0, 1.25, 0.0],
            5.0,
            [0.10, 0.10, 0.07],
            Some(pelvis),
        );
        let chest = rag.add_bone(
            "chest",
            [0.0, 1.55, 0.0],
            7.0,
            [0.14, 0.12, 0.08],
            Some(spine),
        );
        let neck = rag.add_bone(
            "neck",
            [0.0, 1.78, 0.0],
            1.5,
            [0.04, 0.06, 0.04],
            Some(chest),
        );
        let _head = rag.add_bone(
            "head",
            [0.0, 1.95, 0.0],
            5.0,
            [0.10, 0.12, 0.10],
            Some(neck),
        );
        let uarm_l = rag.add_bone(
            "upper_arm_l",
            [-0.20, 1.55, 0.0],
            2.5,
            [0.04, 0.13, 0.04],
            Some(chest),
        );
        let larm_l = rag.add_bone(
            "lower_arm_l",
            [-0.20, 1.28, 0.0],
            1.5,
            [0.03, 0.12, 0.03],
            Some(uarm_l),
        );
        let uarm_r = rag.add_bone(
            "upper_arm_r",
            [0.20, 1.55, 0.0],
            2.5,
            [0.04, 0.13, 0.04],
            Some(chest),
        );
        let larm_r = rag.add_bone(
            "lower_arm_r",
            [0.20, 1.28, 0.0],
            1.5,
            [0.03, 0.12, 0.03],
            Some(uarm_r),
        );
        let uleg_l = rag.add_bone(
            "upper_leg_l",
            [-0.10, 0.70, 0.0],
            5.0,
            [0.05, 0.18, 0.05],
            Some(pelvis),
        );
        let lleg_l = rag.add_bone(
            "lower_leg_l",
            [-0.10, 0.35, 0.0],
            3.0,
            [0.04, 0.18, 0.04],
            Some(uleg_l),
        );
        let foot_l = rag.add_bone(
            "foot_l",
            [-0.10, 0.06, 0.0],
            1.0,
            [0.05, 0.04, 0.10],
            Some(lleg_l),
        );
        let uleg_r = rag.add_bone(
            "upper_leg_r",
            [0.10, 0.70, 0.0],
            5.0,
            [0.05, 0.18, 0.05],
            Some(pelvis),
        );
        let lleg_r = rag.add_bone(
            "lower_leg_r",
            [0.10, 0.35, 0.0],
            3.0,
            [0.04, 0.18, 0.04],
            Some(uleg_r),
        );
        let foot_r = rag.add_bone(
            "foot_r",
            [0.10, 0.06, 0.0],
            1.0,
            [0.05, 0.04, 0.10],
            Some(lleg_r),
        );
        rag.add_ball_joint(pelvis, spine, [0.0, 0.12, 0.0], [0.0, -0.12, 0.0], pi / 6.0);
        rag.add_ball_joint(spine, chest, [0.0, 0.12, 0.0], [0.0, -0.12, 0.0], pi / 8.0);
        rag.add_ball_joint(chest, neck, [0.0, 0.12, 0.0], [0.0, -0.06, 0.0], pi / 6.0);
        rag.add_ball_joint(neck, _head, [0.0, 0.06, 0.0], [0.0, -0.10, 0.0], pi / 4.0);
        rag.add_ball_joint(
            chest,
            uarm_l,
            [-0.14, 0.10, 0.0],
            [0.0, 0.13, 0.0],
            pi / 2.0,
        );
        rag.add_ball_joint(
            uarm_l,
            larm_l,
            [0.0, -0.13, 0.0],
            [0.0, 0.12, 0.0],
            pi / 2.0,
        );
        rag.add_ball_joint(chest, uarm_r, [0.14, 0.10, 0.0], [0.0, 0.13, 0.0], pi / 2.0);
        rag.add_ball_joint(
            uarm_r,
            larm_r,
            [0.0, -0.13, 0.0],
            [0.0, 0.12, 0.0],
            pi / 2.0,
        );
        rag.add_ball_joint(
            pelvis,
            uleg_l,
            [-0.08, -0.10, 0.0],
            [0.0, 0.18, 0.0],
            pi / 2.0,
        );
        rag.add_ball_joint(
            pelvis,
            uleg_r,
            [0.08, -0.10, 0.0],
            [0.0, 0.18, 0.0],
            pi / 2.0,
        );
        let axis_z = [0.0, 0.0, 1.0];
        rag.add_hinge_joint(
            uleg_l,
            lleg_l,
            [0.0, -0.18, 0.0],
            [0.0, 0.18, 0.0],
            axis_z,
            0.0,
            pi * 2.0 / 3.0,
        );
        rag.add_hinge_joint(
            lleg_l,
            foot_l,
            [0.0, -0.18, 0.0],
            [0.0, 0.04, 0.0],
            axis_z,
            -pi / 4.0,
            pi / 4.0,
        );
        rag.add_hinge_joint(
            uleg_r,
            lleg_r,
            [0.0, -0.18, 0.0],
            [0.0, 0.18, 0.0],
            axis_z,
            0.0,
            pi * 2.0 / 3.0,
        );
        rag.add_hinge_joint(
            lleg_r,
            foot_r,
            [0.0, -0.18, 0.0],
            [0.0, 0.04, 0.0],
            axis_z,
            -pi / 4.0,
            pi / 4.0,
        );
        rag
    }
}
impl HumanoidRagdoll {
    /// Build a standard ~15-bone human skeleton descriptor.
    pub fn standard_human() -> HumanoidSkeleton {
        let pi = std::f64::consts::PI;
        let bone = |name: &str, parent: Option<usize>, length: f64, mass: f64| -> BoneDescriptor {
            let r = length * 0.05;
            BoneDescriptor {
                name: name.to_string(),
                parent,
                length,
                mass,
                inertia: [
                    mass * (3.0 * r * r + length * length) / 12.0,
                    mass * (3.0 * r * r + length * length) / 12.0,
                    0.5 * mass * r * r,
                ],
            }
        };
        let sym_joint = |lim: f64| -> JointLimit {
            JointLimit {
                min_angles: [-lim, -lim, -lim],
                max_angles: [lim, lim, lim],
            }
        };
        let bones = vec![
            bone("pelvis", None, 0.20, 8.0),
            bone("spine", Some(0), 0.25, 5.0),
            bone("chest", Some(1), 0.25, 7.0),
            bone("neck", Some(2), 0.10, 1.5),
            bone("head", Some(3), 0.20, 5.0),
            bone("upper_arm_l", Some(2), 0.28, 2.5),
            bone("lower_arm_l", Some(5), 0.25, 1.5),
            bone("hand_l", Some(6), 0.08, 0.6),
            bone("upper_arm_r", Some(2), 0.28, 2.5),
            bone("lower_arm_r", Some(8), 0.25, 1.5),
            bone("hand_r", Some(9), 0.08, 0.6),
            bone("upper_leg_l", Some(0), 0.42, 9.0),
            bone("lower_leg_l", Some(11), 0.38, 4.5),
            bone("upper_leg_r", Some(0), 0.42, 9.0),
            bone("lower_leg_r", Some(13), 0.38, 4.5),
        ];
        let joints = vec![
            sym_joint(pi / 6.0),
            sym_joint(pi / 6.0),
            sym_joint(pi / 8.0),
            sym_joint(pi / 6.0),
            sym_joint(pi / 4.0),
            sym_joint(pi / 2.0),
            JointLimit {
                min_angles: [0.0, -pi / 2.0, 0.0],
                max_angles: [pi * 2.0 / 3.0, pi / 2.0, 0.0],
            },
            sym_joint(pi / 4.0),
            sym_joint(pi / 2.0),
            JointLimit {
                min_angles: [0.0, -pi / 2.0, 0.0],
                max_angles: [pi * 2.0 / 3.0, pi / 2.0, 0.0],
            },
            sym_joint(pi / 4.0),
            sym_joint(pi / 2.0),
            JointLimit {
                min_angles: [0.0, -pi / 6.0, 0.0],
                max_angles: [pi * 2.0 / 3.0, pi / 6.0, 0.0],
            },
            sym_joint(pi / 2.0),
            JointLimit {
                min_angles: [0.0, -pi / 6.0, 0.0],
                max_angles: [pi * 2.0 / 3.0, pi / 6.0, 0.0],
            },
        ];
        HumanoidSkeleton { bones, joints }
    }
}
/// Detects ground contact for ragdoll feet and computes contact response.
///
/// Models a flat ground plane at `y = ground_y`.
#[derive(Debug, Clone)]
pub struct GroundContactDetector {
    /// Y-coordinate of the ground plane.
    pub ground_y: f64,
    /// Contact restitution (bounciness) ∈ \[0, 1\].
    pub restitution: f64,
    /// Contact friction coefficient.
    pub friction: f64,
}
impl GroundContactDetector {
    /// Create a new ground contact detector.
    pub fn new(ground_y: f64, restitution: f64, friction: f64) -> Self {
        Self {
            ground_y,
            restitution,
            friction,
        }
    }
    /// Check whether the bone at `position` (with `radius`) is touching the ground.
    pub fn is_in_contact(&self, position: [f64; 3], radius: f64) -> bool {
        position[1] - radius <= self.ground_y
    }
    /// Compute the penetration depth of a bone below ground.
    pub fn penetration_depth(&self, position: [f64; 3], radius: f64) -> f64 {
        let depth = self.ground_y - (position[1] - radius);
        depth.max(0.0)
    }
    /// Apply a position correction to push the bone above ground.
    ///
    /// Returns the corrected position.
    pub fn correct_position(&self, mut position: [f64; 3], radius: f64) -> [f64; 3] {
        let depth = self.penetration_depth(position, radius);
        if depth > 0.0 {
            position[1] += depth;
        }
        position
    }
    /// Apply a velocity correction (contact response) for a bone hitting the ground.
    ///
    /// Reverses the normal (Y) component with restitution and applies friction to
    /// the tangential (XZ) components.
    ///
    /// Returns the corrected velocity.
    pub fn correct_velocity(
        &self,
        mut velocity: [f64; 3],
        position: [f64; 3],
        radius: f64,
    ) -> [f64; 3] {
        if !self.is_in_contact(position, radius) {
            return velocity;
        }
        if velocity[1] < 0.0 {
            velocity[1] = -velocity[1] * self.restitution;
        }
        let friction_factor = (1.0 - self.friction).max(0.0);
        velocity[0] *= friction_factor;
        velocity[2] *= friction_factor;
        velocity
    }
    /// Update the ground plane height (e.g. for terrain following).
    pub fn set_ground(&mut self, y: f64) {
        self.ground_y = y;
    }
}
/// Per-axis joint angle limits (min/max in radians).
pub struct JointLimit {
    /// Minimum allowed angles per axis (radians).
    pub min_angles: [f64; 3],
    /// Maximum allowed angles per axis (radians).
    pub max_angles: [f64; 3],
}
impl JointLimit {
    /// Clamp angles to \[min, max\] per axis.
    pub fn clamp(&self, angles: [f64; 3]) -> [f64; 3] {
        [
            angles[0].clamp(self.min_angles[0], self.max_angles[0]),
            angles[1].clamp(self.min_angles[1], self.max_angles[1]),
            angles[2].clamp(self.min_angles[2], self.max_angles[2]),
        ]
    }
}
/// Cone-twist joint limit that enforces a swing cone and axial twist range.
///
/// The swing angle is the angle between the current bone direction and the
/// rest direction.  The twist is the rotation about the bone axis.
#[derive(Debug, Clone)]
pub struct ConeTwistLimit {
    /// Half-angle of the swing cone (radians).
    pub cone_angle: f64,
    /// Minimum twist angle about the bone axis (radians).
    pub twist_min: f64,
    /// Maximum twist angle about the bone axis (radians).
    pub twist_max: f64,
    /// Rest orientation of the joint (quaternion w,x,y,z).
    pub rest_orientation: [f64; 4],
}
impl ConeTwistLimit {
    /// Create a symmetric cone-twist limit.
    pub fn symmetric(cone_angle: f64, twist_range: f64) -> Self {
        Self {
            cone_angle,
            twist_min: -twist_range,
            twist_max: twist_range,
            rest_orientation: [1.0, 0.0, 0.0, 0.0],
        }
    }
    /// Create a cone-twist limit with a custom rest orientation.
    pub fn with_rest(cone_angle: f64, twist_min: f64, twist_max: f64, rest: [f64; 4]) -> Self {
        Self {
            cone_angle,
            twist_min,
            twist_max,
            rest_orientation: rest,
        }
    }
    /// Check whether the given orientation `q` satisfies the cone constraint.
    ///
    /// Computes the relative rotation from rest and measures the swing angle.
    pub fn within_cone(&self, q: [f64; 4]) -> bool {
        let q_rel = quat_mul(quat_conjugate(self.rest_orientation), q);
        let w = q_rel[0].clamp(-1.0, 1.0);
        let swing = 2.0 * w.abs().acos();
        swing <= self.cone_angle
    }
    /// Check whether the twist (rotation about the primary axis) is within limits.
    ///
    /// Extracts the twist component about the Y-axis of the bone.
    pub fn within_twist(&self, q: [f64; 4]) -> bool {
        let q_rel = quat_mul(quat_conjugate(self.rest_orientation), q);
        let twist = 2.0
            * (q_rel[0] * q_rel[2] - q_rel[3] * q_rel[1])
                .atan2(1.0 - 2.0 * (q_rel[2] * q_rel[2] + q_rel[1] * q_rel[1]));
        twist >= self.twist_min && twist <= self.twist_max
    }
    /// Apply the cone limit by projecting the orientation back to the cone boundary
    /// if it violates the constraint.
    ///
    /// Returns the corrected orientation.
    pub fn apply_cone_limit(&self, q: [f64; 4]) -> [f64; 4] {
        if self.within_cone(q) {
            return q;
        }

        quat_nlerp(self.rest_orientation, q, 0.5)
    }
    /// Compute the angular error (rad) between the current orientation and the
    /// nearest point on the cone boundary.
    pub fn cone_violation(&self, q: [f64; 4]) -> f64 {
        let q_rel = quat_mul(quat_conjugate(self.rest_orientation), q);
        let w = q_rel[0].clamp(-1.0, 1.0);
        let swing = 2.0 * w.abs().acos();
        (swing - self.cone_angle).max(0.0)
    }
}
/// A single rigid body segment within a ragdoll.
#[derive(Debug, Clone)]
pub struct RagdollBone {
    /// Human-readable name (e.g. "pelvis", "upper_arm_l").
    pub name: String,
    /// Center-of-mass position in world space.
    pub position: [f64; 3],
    /// Orientation quaternion (w, x, y, z).
    pub orientation: [f64; 4],
    /// Linear velocity (world space).
    pub linear_vel: [f64; 3],
    /// Angular velocity (world space).
    pub angular_vel: [f64; 3],
    /// Total mass in kg.
    pub mass: f64,
    /// Cached inverse mass (0 for static/kinematic bones).
    pub inv_mass: f64,
    /// Diagonal principal-axis inertia in local frame.
    pub inertia: [f64; 3],
    /// Cached inverse inertia.
    pub inv_inertia: [f64; 3],
    /// Half-extents of the box approximation used for broad-phase collision.
    pub half_extents: [f64; 3],
    /// Index of the parent bone, if any.
    pub parent: Option<usize>,
    /// Indices of child bones.
    pub children: Vec<usize>,
}
impl RagdollBone {
    /// Compute box inertia tensor diagonal from mass and half-extents.
    fn box_inertia(mass: f64, he: [f64; 3]) -> [f64; 3] {
        let f = mass / 3.0;
        [
            f * (he[1] * he[1] + he[2] * he[2]),
            f * (he[0] * he[0] + he[2] * he[2]),
            f * (he[0] * he[0] + he[1] * he[1]),
        ]
    }
}
/// Articulated ragdoll composed of bones connected by joints.
#[derive(Debug, Clone)]
pub struct Ragdoll {
    /// All bones in this ragdoll.
    pub bones: Vec<RagdollBone>,
    /// Ball-and-socket joints.
    pub ball_joints: Vec<BallSocketJoint>,
    /// Hinge joints.
    pub hinge_joints: Vec<HingeJoint>,
    /// Gravitational acceleration vector.
    pub gravity: [f64; 3],
    /// Current activation state.
    pub state: RagdollState,
}
impl Ragdoll {
    /// Create an empty ragdoll with the given gravity vector.
    pub fn new(gravity: [f64; 3]) -> Self {
        Self {
            bones: Vec::new(),
            ball_joints: Vec::new(),
            hinge_joints: Vec::new(),
            gravity,
            state: RagdollState::Active,
        }
    }
    /// Add a bone and return its index.
    pub fn add_bone(
        &mut self,
        name: &str,
        pos: [f64; 3],
        mass: f64,
        half_extents: [f64; 3],
        parent: Option<usize>,
    ) -> usize {
        let idx = self.bones.len();
        let (inv_mass, inertia, inv_inertia) = if mass > 0.0 {
            let i = RagdollBone::box_inertia(mass, half_extents);
            let ii = [
                if i[0] > 0.0 { 1.0 / i[0] } else { 0.0 },
                if i[1] > 0.0 { 1.0 / i[1] } else { 0.0 },
                if i[2] > 0.0 { 1.0 / i[2] } else { 0.0 },
            ];
            (1.0 / mass, i, ii)
        } else {
            (0.0_f64, [0.0_f64; 3], [0.0_f64; 3])
        };
        if let Some(p) = parent
            && p < self.bones.len()
        {
            self.bones[p].children.push(idx);
        }
        self.bones.push(RagdollBone {
            name: name.to_string(),
            position: pos,
            orientation: [1.0, 0.0, 0.0, 0.0],
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            mass,
            inv_mass,
            inertia,
            inv_inertia,
            half_extents,
            parent,
            children: Vec::new(),
        });
        idx
    }
    /// Add a ball-and-socket joint and return its index.
    pub fn add_ball_joint(
        &mut self,
        bone_a: usize,
        bone_b: usize,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        swing_limit: f64,
    ) -> usize {
        let idx = self.ball_joints.len();
        self.ball_joints.push(BallSocketJoint {
            bone_a,
            bone_b,
            local_anchor_a: anchor_a,
            local_anchor_b: anchor_b,
            swing_limit,
            twist_limit: std::f64::consts::PI,
        });
        idx
    }
    /// Add a hinge joint and return its index.
    #[allow(clippy::too_many_arguments)]
    pub fn add_hinge_joint(
        &mut self,
        bone_a: usize,
        bone_b: usize,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        axis: [f64; 3],
        min_angle: f64,
        max_angle: f64,
    ) -> usize {
        let idx = self.hinge_joints.len();
        self.hinge_joints.push(HingeJoint {
            bone_a,
            bone_b,
            local_anchor_a: anchor_a,
            local_anchor_b: anchor_b,
            local_axis_a: axis,
            min_angle,
            max_angle,
        });
        idx
    }
    /// Return the world-space position of `local_anchor` on the given bone.
    pub fn bone_world_anchor(&self, bone_idx: usize, local_anchor: [f64; 3]) -> [f64; 3] {
        let bone = &self.bones[bone_idx];
        vec_add(bone.position, quat_rotate(bone.orientation, local_anchor))
    }
    /// Semi-implicit Euler integration for all dynamic bones.
    pub fn integrate(&mut self, dt: f64) {
        if self.state == RagdollState::Inactive {
            return;
        }
        let gravity = self.gravity;
        let blend = match self.state {
            RagdollState::Active => 1.0,
            RagdollState::Blending { blend_factor } => blend_factor,
            RagdollState::Inactive => return,
        };
        for bone in &mut self.bones {
            if bone.inv_mass == 0.0 {
                continue;
            }
            let g_scaled = vec_scale(gravity, blend);
            bone.linear_vel = vec_add(bone.linear_vel, vec_scale(g_scaled, dt));
            bone.position = vec_add(bone.position, vec_scale(bone.linear_vel, dt));
            let omega = bone.angular_vel;
            let omega_q: [f64; 4] = [0.0, omega[0], omega[1], omega[2]];
            let dq = quat_mul(omega_q, bone.orientation);
            let q = bone.orientation;
            bone.orientation = quat_normalize([
                q[0] + 0.5 * dq[0] * dt,
                q[1] + 0.5 * dq[1] * dt,
                q[2] + 0.5 * dq[2] * dt,
                q[3] + 0.5 * dq[3] * dt,
            ]);
        }
    }
    /// Position-level correction for a ball-and-socket joint.
    pub fn solve_ball_joint_position(&mut self, joint_idx: usize, _dt: f64) {
        let j = &self.ball_joints[joint_idx];
        let (ia, ib) = (j.bone_a, j.bone_b);
        let la = j.local_anchor_a;
        let lb = j.local_anchor_b;
        let wa = self.bone_world_anchor(ia, la);
        let wb = self.bone_world_anchor(ib, lb);
        let err = vec_sub(wb, wa);
        let inv_a = self.bones[ia].inv_mass;
        let inv_b = self.bones[ib].inv_mass;
        let total_inv = inv_a + inv_b;
        if total_inv == 0.0 {
            return;
        }
        let weight_a = inv_a / total_inv;
        let weight_b = inv_b / total_inv;
        self.bones[ia].position = vec_add(self.bones[ia].position, vec_scale(err, weight_a));
        self.bones[ib].position = vec_sub(self.bones[ib].position, vec_scale(err, weight_b));
    }
    /// Enforce swing limit on a ball-and-socket joint.
    ///
    /// If the relative orientation exceeds the swing limit, the bones
    /// are rotated back to the limit boundary.
    pub fn enforce_swing_limit(&mut self, joint_idx: usize) {
        let j = &self.ball_joints[joint_idx];
        let (ia, ib) = (j.bone_a, j.bone_b);
        let swing_limit = j.swing_limit;
        let qa = self.bones[ia].orientation;
        let qb = self.bones[ib].orientation;
        let q_rel = quat_mul(quat_conjugate(qa), qb);
        let angle = quat_angle([1.0, 0.0, 0.0, 0.0], q_rel);
        if angle > swing_limit && angle > 1e-10 {
            let t = swing_limit / angle;
            let clamped = quat_nlerp([1.0, 0.0, 0.0, 0.0], q_rel, t);
            let new_qb = quat_mul(qa, clamped);
            self.bones[ib].orientation = quat_normalize(new_qb);
        }
    }
    /// Enforce hinge angle limits.
    pub fn enforce_hinge_limits(&mut self, joint_idx: usize) {
        let j = &self.hinge_joints[joint_idx];
        let (ia, ib) = (j.bone_a, j.bone_b);
        let min_angle = j.min_angle;
        let max_angle = j.max_angle;
        let qa = self.bones[ia].orientation;
        let qb = self.bones[ib].orientation;
        let q_rel = quat_mul(quat_conjugate(qa), qb);
        let angle = 2.0 * q_rel[0].acos();
        if angle < min_angle {
            let clamped_q = quat_from_axis_angle(j.local_axis_a, min_angle);
            self.bones[ib].orientation = quat_normalize(quat_mul(qa, clamped_q));
        } else if angle > max_angle {
            let clamped_q = quat_from_axis_angle(j.local_axis_a, max_angle);
            self.bones[ib].orientation = quat_normalize(quat_mul(qa, clamped_q));
        }
    }
    /// Advance simulation by `dt` seconds with `n_iter` constraint iterations.
    pub fn step(&mut self, dt: f64, n_iter: usize) {
        self.integrate(dt);
        for _ in 0..n_iter {
            for ji in 0..self.ball_joints.len() {
                self.solve_ball_joint_position(ji, dt);
                self.enforce_swing_limit(ji);
            }
            for ji in 0..self.hinge_joints.len() {
                let j = &self.hinge_joints[ji];
                let (ia, ib) = (j.bone_a, j.bone_b);
                let la = j.local_anchor_a;
                let lb = j.local_anchor_b;
                let wa = self.bone_world_anchor(ia, la);
                let wb = self.bone_world_anchor(ib, lb);
                let err = vec_sub(wb, wa);
                let inv_a = self.bones[ia].inv_mass;
                let inv_b = self.bones[ib].inv_mass;
                let total_inv = inv_a + inv_b;
                if total_inv == 0.0 {
                    continue;
                }
                let weight_a = inv_a / total_inv;
                let weight_b = inv_b / total_inv;
                self.bones[ia].position =
                    vec_add(self.bones[ia].position, vec_scale(err, weight_a));
                self.bones[ib].position =
                    vec_sub(self.bones[ib].position, vec_scale(err, weight_b));
                self.enforce_hinge_limits(ji);
            }
        }
    }
    /// Activate the ragdoll for full simulation.
    pub fn activate(&mut self) {
        self.state = RagdollState::Active;
    }
    /// Deactivate the ragdoll (kinematic mode).
    pub fn deactivate(&mut self) {
        self.state = RagdollState::Inactive;
        for bone in &mut self.bones {
            bone.linear_vel = [0.0; 3];
            bone.angular_vel = [0.0; 3];
        }
    }
    /// Start blending from animation to physics.
    pub fn start_blending(&mut self, blend_factor: f64) {
        self.state = RagdollState::Blending {
            blend_factor: blend_factor.clamp(0.0, 1.0),
        };
    }
    /// Update the blend factor.
    pub fn set_blend_factor(&mut self, factor: f64) {
        if let RagdollState::Blending { .. } = self.state {
            self.state = RagdollState::Blending {
                blend_factor: factor.clamp(0.0, 1.0),
            };
        }
    }
    /// Whether the ragdoll is currently active.
    pub fn is_active(&self) -> bool {
        !matches!(self.state, RagdollState::Inactive)
    }
    /// Capture the current pose.
    pub fn capture_pose(&self) -> RagdollPose {
        RagdollPose {
            positions: self.bones.iter().map(|b| b.position).collect(),
            orientations: self.bones.iter().map(|b| b.orientation).collect(),
        }
    }
    /// Apply a pose to all bones, overriding positions and orientations.
    pub fn apply_pose(&mut self, pose: &RagdollPose) {
        for (i, bone) in self.bones.iter_mut().enumerate() {
            if i < pose.positions.len() {
                bone.position = pose.positions[i];
            }
            if i < pose.orientations.len() {
                bone.orientation = pose.orientations[i];
            }
        }
    }
    /// Blend between two poses.
    ///
    /// Returns a new pose that is `(1-t)*pose_a + t*pose_b`.
    pub fn blend_poses(pose_a: &RagdollPose, pose_b: &RagdollPose, t: f64) -> RagdollPose {
        let t = t.clamp(0.0, 1.0);
        let n = pose_a.positions.len().min(pose_b.positions.len());
        let mut result = RagdollPose {
            positions: Vec::with_capacity(n),
            orientations: Vec::with_capacity(n),
        };
        for i in 0..n {
            let pos = vec_add(
                vec_scale(pose_a.positions[i], 1.0 - t),
                vec_scale(pose_b.positions[i], t),
            );
            let ori = quat_nlerp(pose_a.orientations[i], pose_b.orientations[i], t);
            result.positions.push(pos);
            result.orientations.push(ori);
        }
        result
    }
    /// Compute the geodesic distance between two poses.
    ///
    /// Returns the sum of angular distances (radians) between corresponding
    /// bone orientations, plus the L2 norm of position differences scaled by
    /// `pos_weight`.
    pub fn compute_pose_distance(
        pose_a: &RagdollPose,
        pose_b: &RagdollPose,
        pos_weight: f64,
    ) -> f64 {
        let n = pose_a.positions.len().min(pose_b.positions.len());
        let mut total = 0.0_f64;
        for i in 0..n {
            let qa = pose_a.orientations[i];
            let qb = pose_b.orientations[i];
            let angle = quat_angle(qa, qb);
            total += angle;
            let dp = vec_sub(pose_a.positions[i], pose_b.positions[i]);
            total += pos_weight * vec_length(dp);
        }
        total
    }
    /// Alpha-blend between two ragdoll poses and apply to this ragdoll.
    ///
    /// `alpha = 0` results in `pose_a`; `alpha = 1` results in `pose_b`.
    /// Positions are linearly interpolated; orientations use NLERP.
    pub fn apply_blended_pose(&mut self, pose_a: &RagdollPose, pose_b: &RagdollPose, alpha: f64) {
        let blended = Self::blend_poses(pose_a, pose_b, alpha);
        self.apply_pose(&blended);
    }
    /// Apply a muscle torque to a bone along a specified joint axis.
    ///
    /// Computes a torque `τ = force_magnitude * axis` and adds it as an
    /// angular impulse `Δω = I⁻¹ · τ · dt` to the bone's angular velocity.
    ///
    /// - `bone_idx`: target bone index
    /// - `axis_world`: world-space axis (should be unit length)
    /// - `force_magnitude`: signed torque magnitude (N·m)
    /// - `dt`: time step (s)
    pub fn apply_muscle_force(
        &mut self,
        bone_idx: usize,
        axis_world: [f64; 3],
        force_magnitude: f64,
        dt: f64,
    ) {
        if bone_idx >= self.bones.len() {
            return;
        }
        let bone = &mut self.bones[bone_idx];
        if bone.inv_mass == 0.0 {
            return;
        }
        let torque = vec_scale(axis_world, force_magnitude);
        let ii = bone.inv_inertia;
        let delta_omega = [
            ii[0] * torque[0] * dt,
            ii[1] * torque[1] * dt,
            ii[2] * torque[2] * dt,
        ];
        bone.angular_vel = vec_add(bone.angular_vel, delta_omega);
    }
    /// Get all bones in a subtree rooted at `root_idx`.
    pub fn subtree(&self, root_idx: usize) -> Vec<usize> {
        let mut result = vec![root_idx];
        let mut queue = vec![root_idx];
        while let Some(idx) = queue.pop() {
            for &child in &self.bones[idx].children {
                result.push(child);
                queue.push(child);
            }
        }
        result
    }
    /// Get the depth of a bone in the hierarchy (root = 0).
    pub fn bone_depth(&self, bone_idx: usize) -> usize {
        let mut depth = 0;
        let mut current = bone_idx;
        while let Some(parent) = self.bones[current].parent {
            depth += 1;
            current = parent;
        }
        depth
    }
    /// Find the root bone (first bone with no parent).
    pub fn root_bone(&self) -> Option<usize> {
        self.bones.iter().position(|b| b.parent.is_none())
    }
    /// Total mass of the ragdoll.
    pub fn total_mass(&self) -> f64 {
        self.bones.iter().map(|b| b.mass).sum()
    }
    /// Center of mass of the ragdoll.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.bones.iter().map(|b| b.mass).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for bone in &self.bones {
            com = vec_add(com, vec_scale(bone.position, bone.mass));
        }
        vec_scale(com, 1.0 / total_mass)
    }
    /// Apply an impulse to all bones (e.g. explosion force).
    pub fn apply_impulse_all(&mut self, impulse: [f64; 3]) {
        for bone in &mut self.bones {
            if bone.inv_mass > 0.0 {
                bone.linear_vel = vec_add(bone.linear_vel, vec_scale(impulse, bone.inv_mass));
            }
        }
    }
    /// Apply linear damping to all bones.
    pub fn apply_damping(&mut self, linear_damping: f64, angular_damping: f64) {
        let lin_factor = (1.0 - linear_damping).clamp(0.0, 1.0);
        let ang_factor = (1.0 - angular_damping).clamp(0.0, 1.0);
        for bone in &mut self.bones {
            bone.linear_vel = vec_scale(bone.linear_vel, lin_factor);
            bone.angular_vel = vec_scale(bone.angular_vel, ang_factor);
        }
    }
    /// Build a ragdoll from a skeleton definition.
    ///
    /// Creates bones from skeleton joints and automatically adds
    /// ball-and-socket joints between parent-child pairs.
    pub fn from_skeleton(
        skeleton: &[SkeletonJoint],
        gravity: [f64; 3],
        default_swing_limit: f64,
    ) -> Self {
        let mut ragdoll = Ragdoll::new(gravity);
        for joint in skeleton {
            ragdoll.add_bone(
                &joint.name,
                joint.bind_position,
                joint.mass,
                joint.half_extents,
                joint.parent,
            );
        }
        for (i, joint) in skeleton.iter().enumerate() {
            if let Some(parent) = joint.parent {
                let anchor_in_parent = vec_sub(joint.bind_position, skeleton[parent].bind_position);
                ragdoll.add_ball_joint(parent, i, anchor_in_parent, [0.0; 3], default_swing_limit);
            }
        }
        ragdoll
    }
}
/// Ball-and-socket joint between two bones.
#[derive(Debug, Clone)]
pub struct BallSocketJoint {
    /// Index of the first bone.
    pub bone_a: usize,
    /// Index of the second bone.
    pub bone_b: usize,
    /// Anchor point in bone A local space.
    pub local_anchor_a: [f64; 3],
    /// Anchor point in bone B local space.
    pub local_anchor_b: [f64; 3],
    /// Maximum swing angle in radians.
    pub swing_limit: f64,
    /// Maximum twist angle in radians.
    pub twist_limit: f64,
}
/// Proportional-derivative controller for a single ragdoll joint angle.
///
/// Given a target angle (from a reference pose) and current angle/velocity,
/// computes the corrective torque to drive the bone toward the target.
#[derive(Debug, Clone)]
pub struct PdJointController {
    /// Proportional gain (N·m/rad).
    pub kp: f64,
    /// Derivative gain (N·m·s/rad).
    pub kd: f64,
    /// Maximum output torque magnitude (N·m).
    pub max_torque: f64,
    /// Target joint angle (rad).
    pub target_angle: f64,
}
impl PdJointController {
    /// Create a new PD joint controller.
    pub fn new(kp: f64, kd: f64, max_torque: f64) -> Self {
        Self {
            kp,
            kd,
            max_torque,
            target_angle: 0.0,
        }
    }
    /// Set the target angle.
    pub fn set_target(&mut self, angle: f64) {
        self.target_angle = angle;
    }
    /// Compute the control torque.
    ///
    /// * `angle` — current joint angle (rad)
    /// * `velocity` — current joint angular velocity (rad/s)
    pub fn torque(&self, angle: f64, velocity: f64) -> f64 {
        let error = self.target_angle - angle;
        let raw = self.kp * error - self.kd * velocity;
        raw.clamp(-self.max_torque, self.max_torque)
    }
    /// Compute the control torque from an angular error directly.
    pub fn torque_from_error(&self, error: f64, velocity: f64) -> f64 {
        let raw = self.kp * error - self.kd * velocity;
        raw.clamp(-self.max_torque, self.max_torque)
    }
    /// Torque required for a given target quaternion vs. current quaternion.
    ///
    /// Extracts the rotation angle difference and applies PD control.
    pub fn torque_from_quaternions(
        &self,
        q_current: [f64; 4],
        q_target: [f64; 4],
        angular_vel: f64,
    ) -> f64 {
        let error_angle = quat_angle(q_current, q_target);
        self.torque_from_error(error_angle, angular_vel)
    }
}
/// Inertia properties estimated from standard body segment geometry.
///
/// Body segments are modeled as solid cylinders.
/// I_xx = I_zz = m*(3r² + h²)/12,  I_yy = m*r²/2.
#[derive(Debug, Clone)]
pub struct BodySegmentInertia {
    /// Segment mass (kg).
    pub mass: f64,
    /// Cylinder radius (m).
    pub radius: f64,
    /// Cylinder height (m).
    pub height: f64,
}
impl BodySegmentInertia {
    /// Create a new body segment.
    pub fn new(mass: f64, radius: f64, height: f64) -> Self {
        Self {
            mass,
            radius,
            height,
        }
    }
    /// Estimate inertia for a standard body segment using body mass fraction.
    ///
    /// * `total_body_mass` — total body mass (kg)
    /// * `mass_fraction` — fraction of body mass for this segment
    /// * `radius` — approximate radius of the segment (m)
    /// * `height` — length/height of the segment (m)
    pub fn from_body_mass(
        total_body_mass: f64,
        mass_fraction: f64,
        radius: f64,
        height: f64,
    ) -> Self {
        Self::new(total_body_mass * mass_fraction, radius, height)
    }
    /// Principal inertia about the transverse axes (X and Z), in kg·m².
    pub fn i_transverse(&self) -> f64 {
        self.mass * (3.0 * self.radius * self.radius + self.height * self.height) / 12.0
    }
    /// Principal inertia about the longitudinal axis (Y), in kg·m².
    pub fn i_longitudinal(&self) -> f64 {
        0.5 * self.mass * self.radius * self.radius
    }
    /// Full principal inertia tensor as \[Ixx, Iyy, Izz\].
    pub fn inertia_tensor(&self) -> [f64; 3] {
        [
            self.i_transverse(),
            self.i_longitudinal(),
            self.i_transverse(),
        ]
    }
    /// Radius of gyration about the transverse axis.
    pub fn radius_of_gyration_transverse(&self) -> f64 {
        if self.mass < f64::EPSILON {
            return 0.0;
        }
        (self.i_transverse() / self.mass).sqrt()
    }
}
/// A humanoid skeleton consisting of bone descriptors and joint limits.
pub struct HumanoidSkeleton {
    /// All bones in the skeleton.
    pub bones: Vec<BoneDescriptor>,
    /// Joint limits per bone.
    pub joints: Vec<JointLimit>,
}
impl HumanoidSkeleton {
    /// Total mass of all bones.
    pub fn total_mass(&self) -> f64 {
        self.bones.iter().map(|b| b.mass).sum()
    }
    /// Center of mass given a pose (FK-based).
    ///
    /// Uses forward kinematics to compute world positions and weights by mass.
    pub fn center_of_mass(&self, pose: &ArticulatedPose) -> [f64; 3] {
        let positions = forward_kinematics(&self.bones, pose);
        let total: f64 = self.bones.iter().map(|b| b.mass).sum();
        if total < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for (i, bone) in self.bones.iter().enumerate() {
            let p = positions[i];
            com[0] += bone.mass * p[0];
            com[1] += bone.mass * p[1];
            com[2] += bone.mass * p[2];
        }
        [com[0] / total, com[1] / total, com[2] / total]
    }
}
