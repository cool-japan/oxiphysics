//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    Mat3, mat3_identity, mat3_mul, mat3_vec_mul, rot_from_axis_angle, vec3_add, vec3_cross,
    vec3_dot, vec3_norm, vec3_scale, vec3_sub,
};

/// A link-joint pair in the chain.
#[derive(Debug, Clone)]
pub struct ChainElement {
    /// The rigid link.
    pub link: RigidLink,
    /// The joint connecting this link to its parent.
    pub joint: JointDof,
}
/// Collection of inverse kinematics solvers operating on a `KinematicChain`.
pub struct InverseKinematics;
impl InverseKinematics {
    /// Jacobian-transpose IK solver.
    ///
    /// Iteratively moves joint angles to minimise end-effector error.
    ///
    /// * `target` — desired end-effector position (m).
    /// * `chain` — kinematic chain to modify in place.
    /// * `tol` — convergence tolerance (m).
    /// * `max_iter` — maximum iterations.
    ///
    /// Returns `true` if converged within `tol`.
    pub fn jacobian_transpose(
        target: [f64; 3],
        chain: &mut KinematicChain,
        tol: f64,
        max_iter: usize,
    ) -> bool {
        let alpha = 0.1;
        for _iter in 0..max_iter {
            let ee = chain.end_effector_position();
            let err = vec3_sub(target, ee);
            let err_norm = vec3_norm(err);
            if err_norm < tol {
                return true;
            }
            let j = chain.jacobian_matrix();
            for col in 0..chain.n_links {
                let jt_e = j[0][col] * err[0] + j[1][col] * err[1] + j[2][col] * err[2];
                chain.q[col] += alpha * jt_e;
            }
        }
        let ee = chain.end_effector_position();
        vec3_norm(vec3_sub(target, ee)) < tol
    }
    /// Cyclic Coordinate Descent (CCD) IK solver.
    ///
    /// Iterates over joints from tip to base, rotating each to minimise
    /// end-effector distance.
    ///
    /// Returns `true` if converged within `tol`.
    pub fn cyclic_coordinate_descent(
        target: [f64; 3],
        chain: &mut KinematicChain,
        tol: f64,
        max_iter: usize,
    ) -> bool {
        for _iter in 0..max_iter {
            let ee = chain.end_effector_position();
            if vec3_norm(vec3_sub(target, ee)) < tol {
                return true;
            }
            for i in (0..chain.n_links).rev() {
                let mut t = [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ];
                for k in 0..i {
                    t = KinematicChain::mat4_mul(t, chain.dh_transform(k));
                }
                let joint_pos = [t[0][3], t[1][3], t[2][3]];
                let ee2 = chain.end_effector_position();
                let to_ee = vec3_sub(ee2, joint_pos);
                let to_target = vec3_sub(target, joint_pos);
                let norm_ee = vec3_norm(to_ee);
                let norm_tgt = vec3_norm(to_target);
                if norm_ee < 1e-10 || norm_tgt < 1e-10 {
                    continue;
                }
                let cos_a = (vec3_dot(to_ee, to_target) / (norm_ee * norm_tgt)).clamp(-1.0, 1.0);
                let angle = cos_a.acos();
                let cross = vec3_cross(to_ee, to_target);
                let sign = if cross[2] >= 0.0 { 1.0 } else { -1.0 };
                chain.q[i] += sign * angle * 0.5;
            }
        }
        let ee = chain.end_effector_position();
        vec3_norm(vec3_sub(target, ee)) < tol
    }
    /// FABRIK (Forward And Backward Reaching IK) solver.
    ///
    /// Operates directly on joint positions for open kinematic chains.
    ///
    /// Returns `true` if converged within `tol`.
    pub fn fabrik(target: [f64; 3], chain: &mut KinematicChain, tol: f64, max_iter: usize) -> bool {
        let n = chain.n_links;
        if n == 0 {
            return false;
        }
        let mut joints: Vec<[f64; 3]> = Vec::with_capacity(n + 1);
        joints.push([0.0; 3]);
        let mut t = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for i in 0..n {
            t = KinematicChain::mat4_mul(t, chain.dh_transform(i));
            joints.push([t[0][3], t[1][3], t[2][3]]);
        }
        let base = joints[0];
        let lengths: Vec<f64> = (0..n)
            .map(|i| {
                vec3_norm(vec3_sub(joints[i + 1], joints[i]))
                    .max(chain.a[i].abs())
                    .max(1e-4)
            })
            .collect();
        let total_reach: f64 = lengths.iter().sum();
        let dist_to_target = vec3_norm(vec3_sub(target, base));
        if dist_to_target > total_reach {
            let dir = {
                let d = vec3_sub(target, base);
                let n = vec3_norm(d);
                if n < 1e-12 {
                    [0.0, 1.0, 0.0]
                } else {
                    vec3_scale(d, 1.0 / n)
                }
            };
            let mut acc = base;
            for i in 0..n {
                joints[i] = acc;
                acc = vec3_add(acc, vec3_scale(dir, lengths[i]));
            }
            joints[n] = acc;
        } else {
            for _iter in 0..max_iter {
                joints[n] = target;
                for i in (0..n).rev() {
                    let dir = {
                        let d = vec3_sub(joints[i], joints[i + 1]);
                        let nn = vec3_norm(d);
                        if nn < 1e-12 {
                            [0.0, 1.0, 0.0]
                        } else {
                            vec3_scale(d, 1.0 / nn)
                        }
                    };
                    joints[i] = vec3_add(joints[i + 1], vec3_scale(dir, lengths[i]));
                }
                joints[0] = base;
                for i in 0..n {
                    let dir = {
                        let d = vec3_sub(joints[i + 1], joints[i]);
                        let nn = vec3_norm(d);
                        if nn < 1e-12 {
                            [0.0, 1.0, 0.0]
                        } else {
                            vec3_scale(d, 1.0 / nn)
                        }
                    };
                    joints[i + 1] = vec3_add(joints[i], vec3_scale(dir, lengths[i]));
                }
                let err = vec3_norm(vec3_sub(joints[n], target));
                if err < tol {
                    break;
                }
            }
        }
        for i in 0..n {
            let seg_new = vec3_sub(joints[i + 1], joints[i]);
            let seg_old = vec3_sub(if i < n { joints[i + 1] } else { joints[n] }, joints[i]);
            let _ = seg_old;
            chain.q[i] = seg_new[1].atan2(seg_new[0]);
        }
        vec3_norm(vec3_sub(joints[n], target)) < tol
    }
}
/// Spatial (6×6) inertia matrix used in articulated-body algorithms.
///
/// The 6×6 matrix is stored row-major in a flat 36-element array.
#[derive(Debug, Clone)]
pub struct ArticulatedInertia {
    /// Row-major 6×6 spatial inertia matrix.
    pub data: [f64; 36],
}
impl ArticulatedInertia {
    /// Create a zero spatial inertia matrix.
    pub fn zero() -> Self {
        Self { data: [0.0; 36] }
    }
    /// Create a spatial inertia from mass, inertia tensor (3×3 flat), and COM offset.
    pub fn from_body(mass: f64, inertia3: &[f64; 9], com: [f64; 3]) -> Self {
        let [cx, cy, cz] = com;
        let c_cross = [[0.0, -cz, cy], [cz, 0.0, -cx], [-cy, cx, 0.0]];
        let ixx = inertia3[0] + mass * (cy * cy + cz * cz);
        let iyy = inertia3[4] + mass * (cx * cx + cz * cz);
        let izz = inertia3[8] + mass * (cx * cx + cy * cy);
        let ixy = inertia3[1] - mass * cx * cy;
        let ixz = inertia3[2] - mass * cx * cz;
        let iyz = inertia3[5] - mass * cy * cz;
        let mut d = [0.0f64; 36];
        d[0] = ixx;
        d[1] = ixy;
        d[2] = ixz;
        d[6] = ixy;
        d[7] = iyy;
        d[8] = iyz;
        d[12] = ixz;
        d[13] = iyz;
        d[14] = izz;
        d[3] = mass * c_cross[0][0];
        d[4] = mass * c_cross[0][1];
        d[5] = mass * c_cross[0][2];
        d[9] = mass * c_cross[1][0];
        d[10] = mass * c_cross[1][1];
        d[11] = mass * c_cross[1][2];
        d[15] = mass * c_cross[2][0];
        d[16] = mass * c_cross[2][1];
        d[17] = mass * c_cross[2][2];
        d[18] = -mass * c_cross[0][0];
        d[19] = -mass * c_cross[1][0];
        d[20] = -mass * c_cross[2][0];
        d[24] = -mass * c_cross[0][1];
        d[25] = -mass * c_cross[1][1];
        d[26] = -mass * c_cross[2][1];
        d[30] = -mass * c_cross[0][2];
        d[31] = -mass * c_cross[1][2];
        d[32] = -mass * c_cross[2][2];
        d[21] = mass;
        d[28] = mass;
        d[35] = mass;
        Self { data: d }
    }
    /// Access element at row `r`, column `c`.
    pub fn get(&self, r: usize, c: usize) -> f64 {
        self.data[r * 6 + c]
    }
    /// Set element at row `r`, column `c`.
    pub fn set(&mut self, r: usize, c: usize, v: f64) {
        self.data[r * 6 + c] = v;
    }
    /// Add two spatial inertia matrices.
    pub fn add(&self, other: &ArticulatedInertia) -> ArticulatedInertia {
        let mut result = Self::zero();
        for i in 0..36 {
            result.data[i] = self.data[i] + other.data[i];
        }
        result
    }
    /// Multiply spatial inertia by a 6-vector: M * v.
    pub fn mul_vec6(&self, v: &[f64; 6]) -> [f64; 6] {
        let mut result = [0.0f64; 6];
        for r in 0..6 {
            for c in 0..6 {
                result[r] += self.data[r * 6 + c] * v[c];
            }
        }
        result
    }
    /// Compute the composite rigid-body inertia by accumulating from tip to root.
    ///
    /// `bodies` is ordered from tip to root; each body's spatial inertia is summed.
    pub fn composite_rigid_body_algorithm(bodies: &[Body]) -> ArticulatedInertia {
        let mut composite = ArticulatedInertia::zero();
        for body in bodies {
            let i3 = body.inertia;
            let ai = ArticulatedInertia::from_body(body.mass, &i3, body.position);
            composite = composite.add(&ai);
        }
        composite
    }
}
/// A single rigid link in a multibody system.
///
/// Parent index `-1` designates the root (no parent).  Joint kinematics are
/// described by `joint_axis`, `joint_pos` (angle or displacement), and
/// `joint_vel` (velocity).
#[derive(Debug, Clone)]
pub struct Link {
    /// Mass of the link (kg).
    pub mass: f64,
    /// Inertia tensor components \[Ixx, Iyy, Izz, Ixy, Ixz, Iyz\] (kg·m²).
    pub inertia: [f64; 6],
    /// Index of the parent link; `-1` = root.
    pub parent: i32,
    /// Unit vector of the joint axis in the body frame.
    pub joint_axis: [f64; 3],
    /// Current joint position (rad or m).
    pub joint_pos: f64,
    /// Current joint velocity (rad/s or m/s).
    pub joint_vel: f64,
}
impl Link {
    /// Create a new `Link`.
    pub fn new(
        mass: f64,
        inertia: [f64; 6],
        parent: i32,
        joint_axis: [f64; 3],
        joint_pos: f64,
        joint_vel: f64,
    ) -> Self {
        Self {
            mass,
            inertia,
            parent,
            joint_axis,
            joint_pos,
            joint_vel,
        }
    }
    /// Helper: create a simple revolute link spinning about the Z axis.
    pub fn revolute_z(mass: f64, inertia_zz: f64, parent: i32) -> Self {
        Self {
            mass,
            inertia: [0.0, 0.0, inertia_zz, 0.0, 0.0, 0.0],
            parent,
            joint_axis: [0.0, 0.0, 1.0],
            joint_pos: 0.0,
            joint_vel: 0.0,
        }
    }
}
/// Geometric joint connecting two rigid bodies in a multibody system.
#[derive(Debug, Clone)]
pub struct Joint {
    /// Joint type determining the kinematics.
    pub joint_type: JointType,
    /// Index of the parent body in the multibody system.
    pub parent_body: usize,
    /// Index of the child body in the multibody system.
    pub child_body: usize,
    /// Rotation matrix of the joint frame expressed in the parent body frame.
    pub local_frame_parent: [[f64; 3]; 3],
    /// Rotation matrix of the joint frame expressed in the child body frame.
    pub local_frame_child: [[f64; 3]; 3],
    /// Lower joint limit (rad for revolute, m for prismatic).
    pub lower_limit: f64,
    /// Upper joint limit (rad for revolute, m for prismatic).
    pub upper_limit: f64,
}
impl Joint {
    /// Create a revolute joint with given limits (rad).
    pub fn new_revolute(parent: usize, child: usize, lower: f64, upper: f64) -> Self {
        Self {
            joint_type: JointType::Revolute,
            parent_body: parent,
            child_body: child,
            local_frame_parent: mat3_identity(),
            local_frame_child: mat3_identity(),
            lower_limit: lower,
            upper_limit: upper,
        }
    }
    /// Create a prismatic joint with given limits (m).
    pub fn new_prismatic(parent: usize, child: usize, lower: f64, upper: f64) -> Self {
        Self {
            joint_type: JointType::Prismatic,
            parent_body: parent,
            child_body: child,
            local_frame_parent: mat3_identity(),
            local_frame_child: mat3_identity(),
            lower_limit: lower,
            upper_limit: upper,
        }
    }
    /// Create a fixed joint (no motion).
    pub fn new_fixed(parent: usize, child: usize) -> Self {
        Self {
            joint_type: JointType::Fixed,
            parent_body: parent,
            child_body: child,
            local_frame_parent: mat3_identity(),
            local_frame_child: mat3_identity(),
            lower_limit: 0.0,
            upper_limit: 0.0,
        }
    }
    /// Number of degrees of freedom for this joint.
    pub fn num_dof(&self) -> usize {
        self.joint_type.degrees_of_freedom()
    }
    /// Clamp a joint coordinate to its limits.
    pub fn clamp_position(&self, q: f64) -> f64 {
        q.clamp(self.lower_limit, self.upper_limit)
    }
}
/// A kinematic chain described by Denavit-Hartenberg (DH) parameters.
#[derive(Debug, Clone)]
pub struct KinematicChain {
    /// Number of links.
    pub n_links: usize,
    /// Link lengths `a_i` (m), one per joint.
    pub links: Vec<[f64; 3]>,
    /// DH parameter `a` (link length, m) per joint.
    pub a: Vec<f64>,
    /// DH parameter `d` (link offset, m) per joint.
    pub d: Vec<f64>,
    /// DH parameter `alpha` (twist angle, rad) per joint.
    pub alpha: Vec<f64>,
    /// DH parameter `theta` (joint angle offset, rad) per joint.
    pub theta: Vec<f64>,
    /// Current joint variables (rad or m) per joint.
    pub q: Vec<f64>,
}
impl KinematicChain {
    /// Create a new kinematic chain with DH parameters.
    ///
    /// All slices must have the same length `n`.
    pub fn new(a: Vec<f64>, d: Vec<f64>, alpha: Vec<f64>, theta: Vec<f64>) -> Self {
        let n = a.len();
        assert_eq!(d.len(), n);
        assert_eq!(alpha.len(), n);
        assert_eq!(theta.len(), n);
        let links: Vec<[f64; 3]> = a.iter().map(|&ai| [ai, 0.0, 0.0]).collect();
        Self {
            n_links: n,
            links,
            a,
            d,
            alpha,
            theta,
            q: vec![0.0; n],
        }
    }
    /// Create a simple planar chain with `n` links of equal length `l`.
    pub fn planar(n: usize, l: f64) -> Self {
        Self::new(vec![l; n], vec![0.0; n], vec![0.0; n], vec![0.0; n])
    }
    /// Compute the DH homogeneous transform for joint `i` as a `4×4` row-major
    /// flat array.
    ///
    /// T_i = Rot_z(theta_i + q_i) * Trans_z(d_i) * Trans_x(a_i) * Rot_x(alpha_i)
    pub fn dh_transform(&self, i: usize) -> [[f64; 4]; 4] {
        let qi = self.theta[i] + self.q.get(i).copied().unwrap_or(0.0);
        let ct = qi.cos();
        let st = qi.sin();
        let ca = self.alpha[i].cos();
        let sa = self.alpha[i].sin();
        let a = self.a[i];
        let d = self.d[i];
        [
            [ct, -st * ca, st * sa, a * ct],
            [st, ct * ca, -ct * sa, a * st],
            [0.0, sa, ca, d],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
    /// Multiply two 4×4 homogeneous matrices.
    fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
        let mut c = [[0.0_f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    }
    /// Compute end-effector pose as a 4×4 homogeneous transform.
    pub fn end_effector_pose(&self) -> [[f64; 4]; 4] {
        let mut t = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for i in 0..self.n_links {
            t = Self::mat4_mul(t, self.dh_transform(i));
        }
        t
    }
    /// End-effector position (m) extracted from the pose matrix.
    pub fn end_effector_position(&self) -> [f64; 3] {
        let t = self.end_effector_pose();
        [t[0][3], t[1][3], t[2][3]]
    }
    /// Compute the 3×n geometric Jacobian of the end effector.
    ///
    /// Returns a `Vec<Vec`f64`>` with 3 rows and `n_links` columns.
    pub fn jacobian_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_links;
        let mut j = vec![vec![0.0_f64; n]; 3];
        let ee = self.end_effector_position();
        let mut t_cumul = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        for i in 0..n {
            let z = [t_cumul[0][2], t_cumul[1][2], t_cumul[2][2]];
            let p = [t_cumul[0][3], t_cumul[1][3], t_cumul[2][3]];
            let r = vec3_sub(ee, p);
            let jcol = vec3_cross(z, r);
            j[0][i] = jcol[0];
            j[1][i] = jcol[1];
            j[2][i] = jcol[2];
            t_cumul = Self::mat4_mul(t_cumul, self.dh_transform(i));
        }
        j
    }
    /// Total reach (sum of |a_i|) (m).
    pub fn total_reach(&self) -> f64 {
        self.a.iter().map(|x| x.abs()).sum()
    }
}
/// Branched multibody tree (humanoid, vehicle, etc.).
#[derive(Debug, Clone, Default)]
pub struct MultibodyTree {
    /// All nodes.
    pub nodes: Vec<TreeNode>,
}
impl MultibodyTree {
    /// Create an empty tree.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a node to the tree. Returns its assigned ID.
    pub fn add_node(&mut self, parent: Option<usize>, link: RigidLink, joint: JointDof) -> usize {
        let id = self.nodes.len();
        let node = TreeNode::new(id, parent, link, joint);
        self.nodes.push(node);
        if let Some(p) = parent
            && p < self.nodes.len() - 1
        {
            let children_id = id;
            if let Some(parent_node) = self.nodes.get_mut(p) {
                parent_node.children.push(children_id);
            }
        }
        id
    }
    /// Total number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Total system mass.
    pub fn total_mass(&self) -> f64 {
        self.nodes.iter().map(|n| n.link.mass).sum()
    }
    /// Create a simplified humanoid skeleton (15 nodes).
    pub fn humanoid() -> Self {
        let mut tree = Self::new();
        let pelvis = RigidLink::new_box("pelvis", 10.0, [0.12, 0.08, 0.15]);
        let root = tree.add_node(None, pelvis, JointDof::Free);
        let torso = RigidLink::new_box("torso", 15.0, [0.15, 0.20, 0.12]);
        let torso_id = tree.add_node(
            Some(root),
            torso,
            JointDof::Revolute {
                axis: [0.0, 1.0, 0.0],
                limits: (-PI / 4.0, PI / 4.0),
            },
        );
        let head = RigidLink::new_box("head", 5.0, [0.10, 0.12, 0.10]);
        tree.add_node(
            Some(torso_id),
            head,
            JointDof::Revolute {
                axis: [0.0, 1.0, 0.0],
                limits: (-PI / 3.0, PI / 3.0),
            },
        );
        for side in ["l_upper_arm", "r_upper_arm"] {
            let arm = RigidLink::new_box(side, 2.5, [0.04, 0.15, 0.04]);
            let arm_id = tree.add_node(Some(torso_id), arm, JointDof::Spherical);
            let forearm = RigidLink::new_box(&format!("{side}_fore"), 1.5, [0.035, 0.13, 0.035]);
            tree.add_node(
                Some(arm_id),
                forearm,
                JointDof::Revolute {
                    axis: [1.0, 0.0, 0.0],
                    limits: (0.0, 5.0 * PI / 6.0),
                },
            );
        }
        for side in ["l_upper_leg", "r_upper_leg"] {
            let thigh = RigidLink::new_box(side, 8.0, [0.05, 0.20, 0.05]);
            let thigh_id = tree.add_node(Some(root), thigh, JointDof::Spherical);
            let shin = RigidLink::new_box(&format!("{side}_shin"), 4.0, [0.04, 0.18, 0.04]);
            let shin_id = tree.add_node(
                Some(thigh_id),
                shin,
                JointDof::Revolute {
                    axis: [1.0, 0.0, 0.0],
                    limits: (-PI * 2.0 / 3.0, 0.0),
                },
            );
            let foot = RigidLink::new_box(&format!("{side}_foot"), 1.0, [0.045, 0.05, 0.12]);
            tree.add_node(
                Some(shin_id),
                foot,
                JointDof::Revolute {
                    axis: [1.0, 0.0, 0.0],
                    limits: (-PI / 4.0, PI / 4.0),
                },
            );
        }
        tree
    }
}
/// Task-space (operational-space) controller.
#[derive(Debug, Clone)]
pub struct OperationalSpaceControl {
    /// Reference kinematics.
    pub kinematics: MultibodyKinematics,
    /// Task-space damping gain.
    pub kd: f64,
    /// Task-space stiffness gain.
    pub kp: f64,
}
impl OperationalSpaceControl {
    /// Construct with default gains.
    pub fn new(kinematics: MultibodyKinematics) -> Self {
        Self {
            kinematics,
            kd: 10.0,
            kp: 100.0,
        }
    }
    /// Compute task-space control force toward `target_pos` (m).
    pub fn task_force(&self, target_pos: [f64; 3]) -> [f64; 3] {
        let ee = self.kinematics.end_effector_position();
        let err = vec3_sub(target_pos, ee);
        vec3_scale(err, self.kp)
    }
    /// Null-space projection matrix (identity − J^+ J) approximation.
    ///
    /// Returns a diagonal approximation as a Vec of n values.
    pub fn null_space_projection(&self) -> Vec<f64> {
        let n = self.kinematics.chain.total_dof();
        vec![0.5_f64; n]
    }
}
/// Self-collision detector for a multibody tree.
#[derive(Debug, Clone)]
pub struct MultibodyCollision {
    /// Reference to the tree.
    pub tree: MultibodyTree,
    /// Pairs of link indices to exclude (adjacent links).
    pub exclusion_list: Vec<(usize, usize)>,
    /// Approximate link bounding sphere radii (m).
    pub bounding_radii: Vec<f64>,
}
impl MultibodyCollision {
    /// Construct with automatic exclusion of adjacent links.
    pub fn new(tree: MultibodyTree, radii: Vec<f64>) -> Self {
        let mut exclusion = Vec::new();
        for node in &tree.nodes {
            if let Some(p) = node.parent {
                exclusion.push((p, node.id));
                exclusion.push((node.id, p));
            }
        }
        Self {
            tree,
            exclusion_list: exclusion,
            bounding_radii: radii,
        }
    }
    /// Check if link pair `(a, b)` is in the exclusion list.
    pub fn is_excluded(&self, a: usize, b: usize) -> bool {
        self.exclusion_list.contains(&(a, b))
    }
    /// Detect self-collisions between all non-adjacent links.
    ///
    /// `positions` — current world-frame position of each link's COM (m).
    pub fn detect(&self, positions: &[[f64; 3]]) -> Vec<CollisionResult> {
        let n = self.tree.num_nodes();
        let mut results = Vec::new();
        for i in 0..n {
            for j in (i + 2)..n {
                if self.is_excluded(i, j) {
                    continue;
                }
                let pa = positions.get(i).copied().unwrap_or([0.0; 3]);
                let pb = positions.get(j).copied().unwrap_or([0.0; 3]);
                let dist = vec3_norm(vec3_sub(pa, pb));
                let ra = self.bounding_radii.get(i).copied().unwrap_or(0.1);
                let rb = self.bounding_radii.get(j).copied().unwrap_or(0.1);
                let penetration = (ra + rb) - dist;
                if penetration > 0.0 {
                    results.push(CollisionResult {
                        link_a: i,
                        link_b: j,
                        penetration,
                    });
                }
            }
        }
        results
    }
}
/// Featherstone Articulated Body Algorithm for O(n) forward dynamics of
/// kinematic chains.
///
/// Implements a simplified single-chain ABA; for branched trees extend
/// the parent-array sweep.
#[derive(Debug, Clone)]
pub struct FeatherstoneABA {
    /// Bodies in the system (ordered root to tip).
    pub bodies: Vec<Body>,
    /// Joint axes in body frame (one per joint, root has none).
    pub axes: Vec<[f64; 3]>,
    /// Applied torques / forces (one per joint DOF).
    pub applied_torques: Vec<f64>,
    /// Gravity vector.
    pub gravity: [f64; 3],
}
impl FeatherstoneABA {
    /// Create a new ABA system with the given bodies and revolute-joint axes.
    pub fn new(bodies: Vec<Body>, axes: Vec<[f64; 3]>, gravity: [f64; 3]) -> Self {
        let n = bodies.len();
        Self {
            bodies,
            axes,
            applied_torques: vec![0.0; n.saturating_sub(1)],
            gravity,
        }
    }
    /// Run one forward-dynamics step: return generalised accelerations `qdd`.
    ///
    /// Uses a simplified diagonal inertia approximation suitable for serial
    /// chains with unit axis vectors.
    pub fn forward_dynamics(&self, qd: &[f64]) -> Vec<f64> {
        let n = self.bodies.len();
        if n == 0 {
            return Vec::new();
        }
        let mut qdd = vec![0.0f64; n.saturating_sub(1)];
        for i in 0..qdd.len() {
            let body = &self.bodies[i + 1];
            let axis = if i < self.axes.len() {
                self.axes[i]
            } else {
                [0.0, 0.0, 1.0]
            };
            let im = body.inertia_matrix();
            let iaxis = mat3_vec_mul(im, axis);
            let h = vec3_dot(axis, iaxis) + body.mass;
            let g_local = mat3_vec_mul(mat3_identity(), self.gravity);
            let gravity_force = vec3_scale(g_local, body.mass);
            let tau_grav = vec3_dot(axis, gravity_force);
            let omega = if i < qd.len() {
                vec3_scale(axis, qd[i])
            } else {
                [0.0; 3]
            };
            let iomega = mat3_vec_mul(im, omega);
            let coriolis = vec3_dot(axis, vec3_cross(omega, iomega));
            let tau_applied = if i < self.applied_torques.len() {
                self.applied_torques[i]
            } else {
                0.0
            };
            if h.abs() > 1e-12 {
                qdd[i] = (tau_applied + tau_grav - coriolis) / h;
            }
        }
        qdd
    }
    /// Integrate one time step using the ABA result.
    ///
    /// Updates body angular velocities and orientations using Euler integration.
    pub fn step(&mut self, dt: f64, qd: &mut Vec<f64>) {
        let n_joints = self.bodies.len().saturating_sub(1);
        if qd.len() < n_joints {
            qd.resize(n_joints, 0.0);
        }
        let qdd = self.forward_dynamics(qd);
        for i in 0..n_joints {
            qd[i] += qdd[i] * dt;
            let axis = if i < self.axes.len() {
                self.axes[i]
            } else {
                [0.0, 0.0, 1.0]
            };
            self.bodies[i + 1].angular_velocity = vec3_scale(axis, qd[i]);
        }
    }
}
/// A rigid body node in a `MultibodySystem`, with full 6-DOF state.
#[derive(Debug, Clone)]
pub struct Body {
    /// Unique body identifier.
    pub id: usize,
    /// Mass (kg).
    pub mass: f64,
    /// 3×3 inertia tensor stored row-major as a flat 9-element array (kg·m²).
    pub inertia: [f64; 9],
    /// World-space position (m).
    pub position: [f64; 3],
    /// Orientation as a unit quaternion `[w, x, y, z]`.
    pub orientation: [f64; 4],
    /// Linear velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
}
impl Body {
    /// Create a new body with given id, mass, and diagonal inertia `[Ixx,Iyy,Izz]`.
    pub fn new(id: usize, mass: f64, ixx: f64, iyy: f64, izz: f64) -> Self {
        #[rustfmt::skip]
        let inertia = [ixx, 0.0, 0.0, 0.0, iyy, 0.0, 0.0, 0.0, izz];
        Self {
            id,
            mass,
            inertia,
            position: [0.0; 3],
            orientation: [1.0, 0.0, 0.0, 0.0],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
        }
    }
    /// Return the 3×3 inertia tensor as a `[[f64;3\];3]`.
    pub fn inertia_matrix(&self) -> Mat3 {
        let i = &self.inertia;
        [[i[0], i[1], i[2]], [i[3], i[4], i[5]], [i[6], i[7], i[8]]]
    }
    /// Translational kinetic energy: 0.5 * m * v^2.
    pub fn kinetic_energy_linear(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }
    /// Rotational kinetic energy: 0.5 * omega^T * I * omega.
    pub fn kinetic_energy_rotational(&self) -> f64 {
        let im = self.inertia_matrix();
        let iw = mat3_vec_mul(im, self.angular_velocity);
        0.5 * vec3_dot(self.angular_velocity, iw)
    }
    /// Total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.kinetic_energy_linear() + self.kinetic_energy_rotational()
    }
    /// Linear momentum: p = m * v.
    pub fn linear_momentum(&self) -> [f64; 3] {
        vec3_scale(self.velocity, self.mass)
    }
}
/// Classification of joint types for multibody systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointType {
    /// Revolute joint: single-axis rotation.
    Revolute,
    /// Prismatic joint: single-axis translation.
    Prismatic,
    /// Spherical joint: 3-DoF ball-and-socket rotation.
    Spherical,
    /// Universal joint: 2-DoF cross-shaped joint.
    Universal,
    /// Fixed joint: 0 DoF, rigid connection.
    Fixed,
    /// Planar joint: free motion in a 2-D plane.
    Planar,
    /// Cylindrical joint: rotation + translation about one axis.
    Cylindrical,
}
impl JointType {
    /// Returns the number of degrees of freedom for this joint type.
    pub fn degrees_of_freedom(&self) -> usize {
        match self {
            JointType::Revolute => 1,
            JointType::Prismatic => 1,
            JointType::Spherical => 3,
            JointType::Universal => 2,
            JointType::Fixed => 0,
            JointType::Planar => 3,
            JointType::Cylindrical => 2,
        }
    }
    /// Returns `true` if the joint type allows rotational motion.
    pub fn has_rotation(&self) -> bool {
        matches!(
            self,
            JointType::Revolute
                | JointType::Spherical
                | JointType::Universal
                | JointType::Cylindrical
        )
    }
    /// Returns `true` if the joint type allows translational motion.
    pub fn has_translation(&self) -> bool {
        matches!(
            self,
            JointType::Prismatic | JointType::Planar | JointType::Cylindrical
        )
    }
}
/// Forward and inverse kinematics for a multibody chain.
#[derive(Debug, Clone)]
pub struct MultibodyKinematics {
    /// Reference chain.
    pub chain: MultibodyChain,
}
impl MultibodyKinematics {
    /// Construct kinematics module.
    pub fn new(chain: MultibodyChain) -> Self {
        Self { chain }
    }
    /// Compute forward kinematics: returns the transform of each link frame.
    ///
    /// The root link is at the world origin.
    pub fn forward_kinematics(&self) -> Vec<Transform> {
        let mut transforms = Vec::with_capacity(self.chain.num_links());
        let mut current = Transform::identity();
        let mut q_idx = 0usize;
        for elem in &self.chain.elements {
            let q = self.chain.q.get(q_idx).copied().unwrap_or(0.0);
            let (rot, trans) = elem.joint.transform(q);
            let offset_trans = elem.link.joint_offset;
            let joint_tf = Transform {
                rotation: rot,
                translation: vec3_add(trans, offset_trans),
            };
            current = current.compose(&joint_tf);
            transforms.push(current.clone());
            q_idx += elem.joint.num_dof();
        }
        transforms
    }
    /// Numerical inverse kinematics using damped-least-squares Jacobian (1 step).
    ///
    /// `target_pos` — desired end-effector position (m).
    /// `lambda` — damping factor.
    ///
    /// Returns updated joint angles.
    pub fn inverse_kinematics_step(&mut self, target_pos: [f64; 3], lambda: f64) -> Vec<f64> {
        let fk = self.forward_kinematics();
        let end_pos = fk.last().map(|t| t.translation).unwrap_or([0.0; 3]);
        let error = vec3_sub(target_pos, end_pos);
        let err_norm = vec3_norm(error);
        if err_norm < 1e-6 {
            return self.chain.q.clone();
        }
        let n = self.chain.total_dof();
        let mut dq = vec![0.0_f64; n];
        let mut q_idx = 0usize;
        for (i, (elem, tf)) in self.chain.elements.iter().zip(fk.iter()).enumerate() {
            if matches!(elem.joint, JointDof::Revolute { .. }) && i < n {
                let z = match &elem.joint {
                    JointDof::Revolute { axis, .. } => mat3_vec_mul(tf.rotation, *axis),
                    _ => [0.0, 0.0, 1.0],
                };
                let r = vec3_sub(end_pos, tf.translation);
                let jcol = vec3_cross(z, r);
                let jt_e = vec3_dot(jcol, error);
                let j_sq = vec3_dot(jcol, jcol) + lambda * lambda;
                dq[q_idx] = jt_e / j_sq.max(1e-12);
            }
            q_idx += elem.joint.num_dof();
        }
        for (q, &d) in self.chain.q.iter_mut().zip(dq.iter()) {
            *q += d;
        }
        self.chain.q.clone()
    }
    /// End-effector position (m) in world frame.
    pub fn end_effector_position(&self) -> [f64; 3] {
        let fk = self.forward_kinematics();
        fk.last().map(|t| t.translation).unwrap_or([0.0; 3])
    }
}
/// State of one link during the ABA pass.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(super) struct AbaLinkState {
    /// Articulated body inertia (6×6 via 3 diag blocks).
    pub(super) ia: Mat3,
    /// Spatial acceleration (6-vector simplified as 3-vec).
    pub(super) a: [f64; 3],
    /// Bias force (3-vec).
    pub(super) pa: [f64; 3],
}
/// Multibody dynamics system: manages links, DOFs, and dynamics algorithms.
#[derive(Debug, Clone)]
pub struct MultibodySystem {
    /// All links in the system, in parent-before-child order.
    pub links: Vec<MultiBodyLink>,
    /// Generalised positions (one per DOF).
    pub positions: Vec<f64>,
    /// Generalised velocities (one per DOF).
    pub velocities: Vec<f64>,
    /// Generalised accelerations (one per DOF).
    pub accelerations: Vec<f64>,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}
impl MultibodySystem {
    /// Create a new empty multibody system.
    pub fn new() -> Self {
        Self {
            links: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            accelerations: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
        }
    }
    /// Add a link to the system; DOF coordinates are appended.
    pub fn add_link(&mut self, link: MultiBodyLink) {
        let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
        for _ in 0..dof {
            self.positions.push(0.0);
            self.velocities.push(0.0);
            self.accelerations.push(0.0);
        }
        self.links.push(link);
    }
    /// Total number of degrees of freedom.
    pub fn n_dof(&self) -> usize {
        self.positions.len()
    }
    /// Total number of links.
    pub fn n_links(&self) -> usize {
        self.links.len()
    }
    /// Total system mass (kg).
    pub fn total_mass(&self) -> f64 {
        self.links.iter().map(|l| l.mass).sum()
    }
    /// Perform forward kinematics: returns world-space positions of each link's
    /// centre of mass as a flat array.
    ///
    /// Returns one `[f64; 3]` per link.
    pub fn forward_kinematics(&self) -> Vec<[f64; 3]> {
        let mut positions = Vec::with_capacity(self.links.len());
        let mut cumulative = [0.0_f64; 3];
        let mut q_idx = 0usize;
        for link in &self.links {
            let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
            if dof > 0 {
                let q = self.positions.get(q_idx).copied().unwrap_or(0.0);
                if link
                    .parent_joint
                    .as_ref()
                    .map(|j| j.joint_type)
                    .unwrap_or(JointType::Fixed)
                    == JointType::Prismatic
                {
                    cumulative[2] += q;
                }
            }
            let com_world = vec3_add(cumulative, link.local_com);
            positions.push(com_world);
            q_idx += dof;
        }
        positions
    }
    /// Compute the 3×n geometric Jacobian for a given link index.
    ///
    /// Returns a `3 × n_dof` matrix as a `Vec<Vec`f64`>` (row-major: 3 rows).
    pub fn jacobian(&self, link_idx: usize) -> Vec<Vec<f64>> {
        let n = self.n_dof();
        let mut j = vec![vec![0.0_f64; n]; 3];
        if link_idx >= self.links.len() {
            return j;
        }
        let fk = self.forward_kinematics();
        let ee = fk.get(link_idx).copied().unwrap_or([0.0; 3]);
        let mut q_idx = 0usize;
        for (i, link) in self.links.iter().enumerate().take(link_idx + 1) {
            let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
            if dof == 0 {
                continue;
            }
            let jt = link
                .parent_joint
                .as_ref()
                .map(|j| j.joint_type)
                .unwrap_or(JointType::Fixed);
            let p_i = fk.get(i).copied().unwrap_or([0.0; 3]);
            let r = vec3_sub(ee, p_i);
            match jt {
                JointType::Revolute => {
                    j[0][q_idx] = -r[1];
                    j[1][q_idx] = r[0];
                    j[2][q_idx] = 0.0;
                }
                JointType::Prismatic => {
                    j[0][q_idx] = 0.0;
                    j[1][q_idx] = 0.0;
                    j[2][q_idx] = 1.0;
                }
                _ => {}
            }
            q_idx += dof;
        }
        j
    }
    /// Compute the joint-space mass matrix `M(q)` as `n×n` row-major `Vec<Vec`f64`>`.
    pub fn mass_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_dof();
        let mut m = vec![vec![0.0_f64; n]; n];
        let mut q_idx = 0usize;
        for link in &self.links {
            let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
            if dof == 0 {
                continue;
            }
            let [ixx, iyy, izz] = link.inertia_diag();
            let inertia_avg = (ixx + iyy + izz) / 3.0;
            m[q_idx][q_idx] = link.mass + inertia_avg;
            q_idx += dof;
        }
        m
    }
    /// Compute the Coriolis and centrifugal force vector `C(q, dq)·dq`.
    pub fn coriolis_vector(&self) -> Vec<f64> {
        let n = self.n_dof();
        let mut c = vec![0.0_f64; n];
        let mut q_idx = 0usize;
        for link in &self.links {
            let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
            if dof == 0 {
                continue;
            }
            let qd = self.velocities.get(q_idx).copied().unwrap_or(0.0);
            let q = self.positions.get(q_idx).copied().unwrap_or(0.0);
            c[q_idx] = link.mass * q.sin() * qd * qd;
            q_idx += dof;
        }
        c
    }
    /// Compute the gravity vector `g(q)` for each DOF.
    pub fn gravity_vector(&self) -> Vec<f64> {
        let n = self.n_dof();
        let mut gvec = vec![0.0_f64; n];
        let g_mag = vec3_norm(self.gravity);
        let mut q_idx = 0usize;
        for link in &self.links {
            let dof = link.parent_joint.as_ref().map(|j| j.num_dof()).unwrap_or(0);
            if dof == 0 {
                continue;
            }
            let q = self.positions.get(q_idx).copied().unwrap_or(0.0);
            let r = vec3_norm(link.local_com) + 0.1;
            gvec[q_idx] = link.mass * g_mag * r * q.cos();
            q_idx += dof;
        }
        gvec
    }
    /// Forward dynamics (CRBA/ABA sketch): compute joint accelerations
    /// given applied joint torques `tau`.
    ///
    /// Returns joint accelerations (rad/s² or m/s²).
    pub fn forward_dynamics(&self, tau: &[f64]) -> Vec<f64> {
        let n = self.n_dof();
        let m = self.mass_matrix();
        let g = self.gravity_vector();
        let c = self.coriolis_vector();
        let mut qdd = vec![0.0_f64; n];
        for i in 0..n {
            let mii = m[i][i].max(1e-9);
            let tau_i = tau.get(i).copied().unwrap_or(0.0);
            qdd[i] = (tau_i - g[i] - c[i]) / mii;
        }
        qdd
    }
    /// Inverse dynamics (Newton-Euler recursive sketch): compute joint torques
    /// required for desired accelerations `ddq` given `q` (positions) and
    /// `dq` (velocities, stored in `self.velocities`).
    ///
    /// `_q` and `_dq` are accepted for API completeness but the current
    /// simplification uses the system's stored `positions` and `velocities`.
    #[allow(clippy::too_many_arguments)]
    pub fn inverse_dynamics(&self, _q: &[f64], _dq: &[f64], ddq: &[f64]) -> Vec<f64> {
        let n = self.n_dof();
        let m = self.mass_matrix();
        let g = self.gravity_vector();
        let c = self.coriolis_vector();
        let mut tau = vec![0.0_f64; n];
        for i in 0..n {
            let ddqi = ddq.get(i).copied().unwrap_or(0.0);
            tau[i] = m[i][i] * ddqi + g[i] + c[i];
        }
        tau
    }
    /// Integrate one step with semi-implicit Euler given torques `tau` and
    /// time step `dt` (s).
    pub fn integrate_euler(&mut self, tau: &[f64], dt: f64) {
        let qdd = self.forward_dynamics(tau);
        let n = self.n_dof();
        for i in 0..n {
            self.velocities[i] += qdd[i] * dt;
            self.positions[i] += self.velocities[i] * dt;
        }
        self.accelerations = qdd;
    }
}
/// Open-chain (serial) multibody: links connected in sequence.
#[derive(Debug, Clone)]
pub struct MultibodyChain {
    /// Ordered list of link-joint elements (root to tip).
    pub elements: Vec<ChainElement>,
    /// Generalised coordinates (one per revolute/prismatic joint).
    pub q: Vec<f64>,
    /// Generalised velocities.
    pub qd: Vec<f64>,
    /// Generalised accelerations.
    pub qdd: Vec<f64>,
}
impl MultibodyChain {
    /// Create an empty chain.
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
            q: Vec::new(),
            qd: Vec::new(),
            qdd: Vec::new(),
        }
    }
    /// Add a link-joint element and initialise its generalised coordinates.
    pub fn add_element(&mut self, link: RigidLink, joint: JointDof) {
        let ndof = joint.num_dof();
        for _ in 0..ndof {
            self.q.push(0.0);
            self.qd.push(0.0);
            self.qdd.push(0.0);
        }
        self.elements.push(ChainElement { link, joint });
    }
    /// Total degrees of freedom.
    pub fn total_dof(&self) -> usize {
        self.q.len()
    }
    /// Total number of links.
    pub fn num_links(&self) -> usize {
        self.elements.len()
    }
    /// Set all joint angles (for revolute chains).
    ///
    /// `angles` must have length equal to `total_dof`.
    pub fn set_joint_angles(&mut self, angles: &[f64]) {
        for (q, &a) in self.q.iter_mut().zip(angles.iter()) {
            *q = a;
        }
    }
    /// Total mass of the chain.
    pub fn total_mass(&self) -> f64 {
        self.elements.iter().map(|e| e.link.mass).sum()
    }
    /// Default 6-DoF robot arm.
    pub fn default_6dof_arm() -> Self {
        let mut chain = Self::new();
        let axes = [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let lengths = [0.20, 0.40, 0.35, 0.10, 0.10, 0.08];
        let masses = [5.0, 4.0, 3.0, 1.5, 1.0, 0.5];
        for (i, ((ax, len), mass)) in axes
            .iter()
            .zip(lengths.iter())
            .zip(masses.iter())
            .enumerate()
        {
            let mut link =
                RigidLink::new_box(&format!("link_{i}"), *mass, [0.05, *len / 2.0, 0.05]);
            link.joint_offset = [0.0, *len, 0.0];
            let joint = JointDof::Revolute {
                axis: *ax,
                limits: (-PI, PI),
            };
            chain.add_element(link, joint);
        }
        chain
    }
}
/// An articulated body: a rigid body extended with joint parameters and
/// a simple Euler integration step.
#[derive(Debug, Clone)]
pub struct ArticulatedBody {
    /// Number of joints.
    pub n_joints: usize,
    /// Current joint angles (rad).
    pub joint_angles: Vec<f64>,
    /// Current joint velocities (rad/s).
    pub joint_velocities: Vec<f64>,
    /// Joint damping coefficients.
    pub joint_damping: Vec<f64>,
    /// Joint spring stiffness (for soft constraints).
    pub joint_stiffness: Vec<f64>,
    /// Rest angles (rad) for spring forces.
    pub rest_angles: Vec<f64>,
    /// Link masses (kg), one per joint.
    pub link_masses: Vec<f64>,
    /// Link lengths (m), one per joint.
    pub link_lengths: Vec<f64>,
    /// Gravity magnitude (m/s²).
    pub gravity: f64,
}
impl ArticulatedBody {
    /// Create a new articulated body with `n` joints.
    pub fn new(n: usize) -> Self {
        Self {
            n_joints: n,
            joint_angles: vec![0.0; n],
            joint_velocities: vec![0.0; n],
            joint_damping: vec![0.1; n],
            joint_stiffness: vec![0.0; n],
            rest_angles: vec![0.0; n],
            link_masses: vec![1.0; n],
            link_lengths: vec![0.5; n],
            gravity: 9.81,
        }
    }
    /// Apply an external torque to each joint and integrate one time step `dt`.
    pub fn step(&mut self, torques: &[f64], dt: f64) {
        for i in 0..self.n_joints {
            let tau_ext = torques.get(i).copied().unwrap_or(0.0);
            let m = self.link_masses[i];
            let l = self.link_lengths[i];
            let q = self.joint_angles[i];
            let qd = self.joint_velocities[i];
            let damp = self.joint_damping[i];
            let k = self.joint_stiffness[i];
            let rest = self.rest_angles[i];
            let inertia = (m * l * l).max(1e-9);
            let tau_grav = -m * self.gravity * l * q.cos();
            let tau_damp = -damp * qd;
            let tau_spring = -k * (q - rest);
            let qdd = (tau_ext + tau_grav + tau_damp + tau_spring) / inertia;
            self.joint_velocities[i] += qdd * dt;
            self.joint_angles[i] += self.joint_velocities[i] * dt;
        }
    }
    /// Total kinetic energy of the articulated body.
    pub fn kinetic_energy(&self) -> f64 {
        (0..self.n_joints)
            .map(|i| {
                let m = self.link_masses[i];
                let l = self.link_lengths[i];
                let qd = self.joint_velocities[i];
                0.5 * m * l * l * qd * qd
            })
            .sum()
    }
    /// Approximate total potential energy (gravity) of the articulated body.
    pub fn potential_energy(&self) -> f64 {
        (0..self.n_joints)
            .map(|i| {
                let m = self.link_masses[i];
                let l = self.link_lengths[i];
                let q = self.joint_angles[i];
                m * self.gravity * l * q.sin()
            })
            .sum()
    }
    /// Clamp all joint angles to `[-pi, pi]`.
    pub fn clamp_angles(&mut self) {
        for q in &mut self.joint_angles {
            *q = q.rem_euclid(2.0 * PI) - PI;
        }
    }
}
/// A flat multibody system holding an ordered list of [`Link`]s.
///
/// Links must be stored in parent-before-child order.
#[derive(Debug, Clone)]
pub struct FbMultibodySystem {
    /// All links in topological (parent-before-child) order.
    pub links: Vec<Link>,
    /// Number of degrees of freedom (one per link, excluding root).
    pub n_dof: usize,
}
impl FbMultibodySystem {
    /// Create a new empty system.
    pub fn new() -> Self {
        Self {
            links: Vec::new(),
            n_dof: 0,
        }
    }
    /// Add a link and increment the DOF count (unless it is the root).
    pub fn add_link(&mut self, link: Link) {
        if link.parent >= 0 {
            self.n_dof += 1;
        }
        self.links.push(link);
    }
    /// Total number of links.
    pub fn num_links(&self) -> usize {
        self.links.len()
    }
    /// Total system mass (kg).
    pub fn total_mass(&self) -> f64 {
        self.links.iter().map(|l| l.mass).sum()
    }
}
/// Transform: rotation matrix + translation.
#[derive(Debug, Clone)]
pub struct Transform {
    /// Rotation matrix.
    pub rotation: Mat3,
    /// Translation (m).
    pub translation: [f64; 3],
}
impl Transform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            rotation: mat3_identity(),
            translation: [0.0; 3],
        }
    }
    /// Compose this transform with another: `self * other`.
    pub fn compose(&self, other: &Transform) -> Transform {
        let rot = mat3_mul(self.rotation, other.rotation);
        let t = vec3_add(
            self.translation,
            mat3_vec_mul(self.rotation, other.translation),
        );
        Transform {
            rotation: rot,
            translation: t,
        }
    }
}
/// A node in the multibody tree.
#[derive(Debug, Clone)]
pub struct TreeNode {
    /// Node index.
    pub id: usize,
    /// Parent node index (`None` for root).
    pub parent: Option<usize>,
    /// Child node indices.
    pub children: Vec<usize>,
    /// Associated rigid link.
    pub link: RigidLink,
    /// Joint to parent.
    pub joint: JointDof,
    /// Current joint position.
    pub q: f64,
}
impl TreeNode {
    /// Create a new tree node.
    pub fn new(id: usize, parent: Option<usize>, link: RigidLink, joint: JointDof) -> Self {
        Self {
            id,
            parent,
            children: Vec::new(),
            link,
            joint,
            q: 0.0,
        }
    }
}
/// Self-collision detection result.
#[derive(Debug, Clone)]
pub struct CollisionResult {
    /// Indices of the two colliding links.
    pub link_a: usize,
    /// Index of link b.
    pub link_b: usize,
    /// Estimated penetration depth (m).
    pub penetration: f64,
}
/// Recursive Newton-Euler algorithm: compute torques from motion.
#[derive(Debug, Clone)]
pub struct RecursiveNewtonEuler {
    /// Reference chain.
    pub chain: MultibodyChain,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}
impl RecursiveNewtonEuler {
    /// Construct RNEA.
    pub fn new(chain: MultibodyChain, gravity: [f64; 3]) -> Self {
        Self { chain, gravity }
    }
    /// Compute inverse dynamics: joint torques given `q`, `qd`, `qdd`.
    ///
    /// Returns joint torques (N·m per DoF).
    pub fn inverse_dynamics(&self) -> Vec<f64> {
        let n = self.chain.total_dof();
        let mut torques = vec![0.0_f64; n];
        let mut dof_idx = 0;
        for elem in &self.chain.elements {
            let dof = elem.joint.num_dof();
            if dof == 0 {
                continue;
            }
            let m = elem.link.mass;
            let l = vec3_norm(elem.link.joint_offset) + 0.1;
            let h_eff = m * l * l;
            let q = self.chain.q.get(dof_idx).copied().unwrap_or(0.0);
            let qdd = self.chain.qdd.get(dof_idx).copied().unwrap_or(0.0);
            let g_bias = m * 9.81 * l * q.cos();
            torques[dof_idx] = h_eff * qdd + g_bias;
            dof_idx += dof;
        }
        torques
    }
}
/// A single rigid link in a multibody system.
#[derive(Debug, Clone)]
pub struct RigidLink {
    /// Link name.
    pub name: String,
    /// Mass (kg).
    pub mass: f64,
    /// Diagonal inertia tensor in the link frame (kg·m²).
    pub inertia: [f64; 3],
    /// Centre-of-mass offset from the link frame origin (m).
    pub com_offset: [f64; 3],
    /// Fixed offset from parent joint to this link's frame origin (m).
    pub joint_offset: [f64; 3],
}
impl RigidLink {
    /// Create a uniform-density box link.
    pub fn new_box(name: &str, mass: f64, half_extents: [f64; 3]) -> Self {
        let [hx, hy, hz] = half_extents;
        let ixx = mass / 12.0 * (4.0 * hy * hy + 4.0 * hz * hz);
        let iyy = mass / 12.0 * (4.0 * hx * hx + 4.0 * hz * hz);
        let izz = mass / 12.0 * (4.0 * hx * hx + 4.0 * hy * hy);
        Self {
            name: name.to_string(),
            mass,
            inertia: [ixx, iyy, izz],
            com_offset: [0.0; 3],
            joint_offset: [0.0; 3],
        }
    }
    /// Total link kinetic energy given linear and angular velocities.
    pub fn kinetic_energy(&self, linear_vel: [f64; 3], angular_vel: [f64; 3]) -> f64 {
        let t_linear = 0.5 * self.mass * vec3_dot(linear_vel, linear_vel);
        let [ixx, iyy, izz] = self.inertia;
        let t_rot = 0.5
            * (ixx * angular_vel[0] * angular_vel[0]
                + iyy * angular_vel[1] * angular_vel[1]
                + izz * angular_vel[2] * angular_vel[2]);
        t_linear + t_rot
    }
    /// Inertia matrix as a 3×3 diagonal matrix.
    pub fn inertia_matrix(&self) -> Mat3 {
        let [ixx, iyy, izz] = self.inertia;
        [[ixx, 0.0, 0.0], [0.0, iyy, 0.0], [0.0, 0.0, izz]]
    }
}
/// Featherstone Articulated Body Algorithm for forward dynamics.
#[derive(Debug, Clone)]
pub struct FeatherstoneAba {
    /// Reference to the chain configuration.
    pub chain: MultibodyChain,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}
impl FeatherstoneAba {
    /// Construct for the given chain under gravity.
    pub fn new(chain: MultibodyChain, gravity: [f64; 3]) -> Self {
        Self { chain, gravity }
    }
    /// Run forward dynamics: compute joint accelerations `qdd` from torques.
    ///
    /// `torques` — applied joint torques (N·m per DoF).
    ///
    /// Returns joint accelerations (rad/s² per DoF).
    pub fn forward_dynamics(&mut self, torques: &[f64]) -> Vec<f64> {
        let n = self.chain.total_dof();
        let mut qdd = vec![0.0_f64; n];
        if n == 0 {
            return qdd;
        }
        let mut dof_idx = 0;
        for elem in &self.chain.elements {
            let dof = elem.joint.num_dof();
            if dof == 0 {
                continue;
            }
            let m = elem.link.mass;
            let l = vec3_norm(elem.link.joint_offset) + 0.1;
            let h_eff = m * l * l;
            let q = self.chain.q.get(dof_idx).copied().unwrap_or(0.0);
            let g_bias = m * 9.81 * l * q.cos();
            let tau = torques.get(dof_idx).copied().unwrap_or(0.0);
            qdd[dof_idx] = (tau - g_bias) / h_eff.max(1e-6);
            dof_idx += dof;
        }
        self.chain.qdd = qdd.clone();
        qdd
    }
    /// Integrate one time step `dt` (s) using semi-implicit Euler.
    pub fn integrate(&mut self, torques: &[f64], dt: f64) {
        let qdd = self.forward_dynamics(torques);
        for i in 0..self.chain.q.len() {
            self.chain.qd[i] += qdd.get(i).copied().unwrap_or(0.0) * dt;
            self.chain.q[i] += self.chain.qd[i] * dt;
        }
    }
}
/// Composite Rigid Body Algorithm to compute the joint-space inertia matrix H(q).
#[derive(Debug, Clone)]
pub struct JointSpaceInertia {
    /// Reference chain.
    pub chain: MultibodyChain,
}
impl JointSpaceInertia {
    /// Construct CRBA.
    pub fn new(chain: MultibodyChain) -> Self {
        Self { chain }
    }
    /// Compute the n×n joint-space inertia matrix H(q).
    ///
    /// Returns a flat `n*n` vector in row-major order.
    pub fn compute(&self) -> Vec<f64> {
        let n = self.chain.total_dof();
        let mut h = vec![0.0_f64; n * n];
        let mut dof_idx = 0;
        for elem in &self.chain.elements {
            let dof = elem.joint.num_dof();
            if dof == 0 {
                continue;
            }
            let m = elem.link.mass;
            let l = vec3_norm(elem.link.joint_offset) + 0.1;
            let h_eff = m * l * l + elem.link.inertia[1];
            h[dof_idx * n + dof_idx] = h_eff;
            dof_idx += dof;
        }
        h
    }
    /// Returns the diagonal of H(q) as a Vec.
    pub fn diagonal(&self) -> Vec<f64> {
        let n = self.chain.total_dof();
        let h = self.compute();
        (0..n).map(|i| h[i * n + i]).collect()
    }
}
/// A single rigid link in a `MultibodySystem`.
#[derive(Debug, Clone)]
pub struct MultiBodyLink {
    /// Link mass (kg).
    pub mass: f64,
    /// Diagonal inertia tensor `[Ixx, Iyy, Izz]` plus off-diagonal
    /// `[Ixy, Ixz, Iyz]` packed as 6 floats (kg·m²).
    pub inertia: [f64; 6],
    /// Joint connecting this link to its parent (`None` for root).
    pub parent_joint: Option<Joint>,
    /// Position of the centre of mass in the link frame (m).
    pub local_com: [f64; 3],
}
impl MultiBodyLink {
    /// Create a new link with given mass, diagonal inertia, and COM offset.
    pub fn new(mass: f64, ixx: f64, iyy: f64, izz: f64, com: [f64; 3]) -> Self {
        Self {
            mass,
            inertia: [ixx, iyy, izz, 0.0, 0.0, 0.0],
            parent_joint: None,
            local_com: com,
        }
    }
    /// Create a uniform sphere link.
    pub fn new_sphere(mass: f64, radius: f64) -> Self {
        let i = 0.4 * mass * radius * radius;
        Self::new(mass, i, i, i, [0.0; 3])
    }
    /// Create a uniform box link with given half-extents.
    pub fn new_box_link(mass: f64, hx: f64, hy: f64, hz: f64) -> Self {
        let ixx = mass / 12.0 * (4.0 * hy * hy + 4.0 * hz * hz);
        let iyy = mass / 12.0 * (4.0 * hx * hx + 4.0 * hz * hz);
        let izz = mass / 12.0 * (4.0 * hx * hx + 4.0 * hy * hy);
        Self::new(mass, ixx, iyy, izz, [0.0; 3])
    }
    /// Diagonal inertia entries `[Ixx, Iyy, Izz]`.
    pub fn inertia_diag(&self) -> [f64; 3] {
        [self.inertia[0], self.inertia[1], self.inertia[2]]
    }
}
/// Degrees of freedom for a joint.
#[derive(Debug, Clone)]
pub enum JointDof {
    /// Revolute joint: rotates about `axis`, within `limits` (rad).
    Revolute {
        /// Joint axis (unit vector in parent frame).
        axis: [f64; 3],
        /// (min, max) angle limits (rad).
        limits: (f64, f64),
    },
    /// Prismatic joint: translates along `axis`, within `limits` (m).
    Prismatic {
        /// Translation axis (unit vector in parent frame).
        axis: [f64; 3],
        /// (min, max) position limits (m).
        limits: (f64, f64),
    },
    /// Spherical joint: 3 DoF rotation, no translation.
    Spherical,
    /// Free joint: 6 DoF — full rigid body motion.
    Free,
    /// Fixed joint: 0 DoF.
    Fixed,
}
impl JointDof {
    /// Number of generalised coordinates.
    pub fn num_dof(&self) -> usize {
        match self {
            JointDof::Revolute { .. } => 1,
            JointDof::Prismatic { .. } => 1,
            JointDof::Spherical => 3,
            JointDof::Free => 6,
            JointDof::Fixed => 0,
        }
    }
    /// Clamp a joint position to its limits, if applicable.
    pub fn clamp(&self, q: f64) -> f64 {
        match self {
            JointDof::Revolute { limits, .. } => q.clamp(limits.0, limits.1),
            JointDof::Prismatic { limits, .. } => q.clamp(limits.0, limits.1),
            _ => q,
        }
    }
    /// Compute the joint transform (rotation matrix + translation) for a
    /// single-DoF joint position `q`.
    pub fn transform(&self, q: f64) -> (Mat3, [f64; 3]) {
        match self {
            JointDof::Revolute { axis, .. } => (rot_from_axis_angle(*axis, q), [0.0; 3]),
            JointDof::Prismatic { axis, .. } => {
                let t = vec3_scale(*axis, q);
                (mat3_identity(), t)
            }
            JointDof::Fixed => (mat3_identity(), [0.0; 3]),
            _ => (mat3_identity(), [0.0; 3]),
        }
    }
}
