//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    add3, add3w, cross3, dot3, joint_displacement, mat3_mul_vec, scale3, scale3w,
};
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};

/// Joint type for a 2-D articulated chain.
#[derive(Debug, Clone, PartialEq)]
pub enum JointType {
    /// Revolute (rotation) joint.
    Revolute,
    /// Prismatic (translation) joint.
    Prismatic,
    /// Fixed (rigid) joint.
    Fixed,
    /// Spherical (ball-and-socket) joint: 3 rotational DOF.
    Spherical,
    /// Planar joint: 2 translation + 1 rotation DOF.
    Planar,
    /// Universal (Cardan/Hooke) joint: 2 rotational DOF.
    Universal,
    /// Cylindrical joint: 1 rotation + 1 translation DOF.
    Cylindrical,
}
impl JointType {
    /// Number of degrees of freedom for this joint type.
    pub fn degrees_of_freedom(&self) -> usize {
        match self {
            JointType::Fixed => 0,
            JointType::Revolute | JointType::Prismatic => 1,
            JointType::Universal | JointType::Cylindrical => 2,
            JointType::Spherical | JointType::Planar => 3,
        }
    }
    /// Whether this joint allows rotation.
    pub fn has_rotation(&self) -> bool {
        matches!(
            self,
            JointType::Revolute
                | JointType::Spherical
                | JointType::Universal
                | JointType::Cylindrical
                | JointType::Planar
        )
    }
    /// Whether this joint allows translation.
    pub fn has_translation(&self) -> bool {
        matches!(
            self,
            JointType::Prismatic | JointType::Planar | JointType::Cylindrical
        )
    }
}
/// A screw joint combining rotation and translation along the same axis.
///
/// For every radian of rotation, the joint translates by `pitch / (2π)` metres.
#[derive(Debug, Clone)]
pub struct ScrewJoint {
    /// Screw pitch (m per full revolution = 2π rad).
    pub pitch: f64,
    /// Current angular position (rad).
    pub angle: f64,
    /// Current linear displacement (m).
    pub displacement: f64,
}
impl ScrewJoint {
    /// Create a screw joint with the given pitch.
    pub fn new(pitch: f64) -> Self {
        Self {
            pitch,
            angle: 0.0,
            displacement: 0.0,
        }
    }
    /// Advance the screw by `delta_angle` radians.
    pub fn advance(&mut self, delta_angle: f64) {
        self.angle += delta_angle;
        self.displacement = self.angle * self.pitch / (2.0 * std::f64::consts::PI);
    }
    /// Lead (displacement per radian).
    pub fn lead(&self) -> f64 {
        self.pitch / (2.0 * std::f64::consts::PI)
    }
}
/// A single link with joint information in a 2-D chain.
#[derive(Debug, Clone)]
pub struct ArticulatedLink2D {
    /// Body properties.
    pub body: LinkBody,
    /// Joint type.
    pub joint_type: JointType,
    /// Current joint angle \[rad\] (for revolute) or displacement \[m\] (for prismatic).
    pub joint_angle: f64,
    /// Current joint velocity \[rad/s or m/s\].
    pub joint_velocity: f64,
    /// Applied joint torque \[N·m\] or force \[N\].
    pub joint_torque: f64,
    /// Index of parent link, `None` for root.
    pub parent: Option<usize>,
}
/// Per-link external force (generalised force along the revolute z axis, N·m).
#[derive(Debug, Clone, Default)]
pub struct ExternalForce {
    /// Link index.
    pub link_idx: usize,
    /// Torque contribution along the joint axis (N·m).
    pub torque: f64,
}
/// A single link in an articulated chain.
#[derive(Debug, Clone)]
pub struct ArticulatedLink {
    /// Unique identifier for this link.
    pub id: LinkId,
    /// Parent link, or `None` for the root.
    pub parent: Option<LinkId>,
    /// Child links.
    pub children: Vec<LinkId>,
    /// Joint type connecting this link to its parent.
    pub joint_type: JointDof,
    /// Current joint position (generalised coordinate q).
    pub joint_position: f64,
    /// Current joint velocity (generalised velocity q_dot).
    pub joint_velocity: f64,
    /// Optional joint limits `(min, max)`.
    pub joint_limit: Option<(f64, f64)>,
    /// Pose of this link relative to the parent joint origin.
    pub local_transform: Transform,
    /// Mass of this link in kg.
    pub mass: f64,
    /// Local inertia tensor (3×3).
    pub inertia: [[f64; 3]; 3],
}
/// Degrees of freedom for a joint.
#[derive(Debug, Clone)]
pub enum JointDof {
    /// Fixed joint: no relative motion.
    Fixed,
    /// Revolute joint: rotation around the given axis.
    Revolute(
        /// Joint rotation axis in local (parent) frame.
        [f64; 3],
    ),
    /// Prismatic joint: translation along the given axis.
    Prismatic(
        /// Joint translation axis in local (parent) frame.
        [f64; 3],
    ),
    /// Spherical joint: free rotation (3 DOF).
    Spherical,
    /// Planar joint: motion in a plane (2 translation + 1 rotation DOF).
    Planar,
}
impl JointDof {
    /// Number of degrees of freedom.
    pub fn num_dof(&self) -> usize {
        match self {
            JointDof::Fixed => 0,
            JointDof::Revolute(_) | JointDof::Prismatic(_) => 1,
            JointDof::Spherical => 3,
            JointDof::Planar => 3,
        }
    }
    /// Clamp a scalar joint position within typical limits.
    ///
    /// For revolute joints the value is clamped to `[-π, π]` by default;
    /// for prismatic joints it is returned unchanged.
    pub fn clamp(&self, value: f64) -> f64 {
        match self {
            JointDof::Fixed => 0.0,
            JointDof::Revolute(_) => value.clamp(-std::f64::consts::PI, std::f64::consts::PI),
            JointDof::Prismatic(_) => value,
            JointDof::Spherical | JointDof::Planar => value,
        }
    }
    /// Compute the relative transform (rotation, translation) for a given
    /// scalar joint position `q`.
    ///
    /// Returns `(rotation_3x3, translation_3)`.
    pub fn transform(&self, q: f64) -> ([[f64; 3]; 3], [f64; 3]) {
        match self {
            JointDof::Fixed => (
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
            ),
            JointDof::Revolute(axis) => {
                let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
                if n < 1e-12 {
                    return (
                        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        [0.0; 3],
                    );
                }
                let u = [axis[0] / n, axis[1] / n, axis[2] / n];
                let c = q.cos();
                let s = q.sin();
                let t = 1.0 - c;
                let rot = [
                    [
                        t * u[0] * u[0] + c,
                        t * u[0] * u[1] - s * u[2],
                        t * u[0] * u[2] + s * u[1],
                    ],
                    [
                        t * u[0] * u[1] + s * u[2],
                        t * u[1] * u[1] + c,
                        t * u[1] * u[2] - s * u[0],
                    ],
                    [
                        t * u[0] * u[2] - s * u[1],
                        t * u[1] * u[2] + s * u[0],
                        t * u[2] * u[2] + c,
                    ],
                ];
                (rot, [0.0; 3])
            }
            JointDof::Prismatic(axis) => {
                let trans = [axis[0] * q, axis[1] * q, axis[2] * q];
                ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], trans)
            }
            JointDof::Spherical | JointDof::Planar => (
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
            ),
        }
    }
}
/// Body properties for a single link in a 2-D articulated chain.
#[derive(Debug, Clone)]
pub struct LinkBody {
    /// Mass of the link \[kg\].
    pub mass: f64,
    /// Length of the link \[m\].
    pub length: f64,
    /// Scalar moment of inertia (2-D) \[kg·m²\].
    pub inertia: f64,
    /// Position of the link origin in world space \[x, y\].
    pub position: [f64; 2],
    /// Orientation angle \[rad\].
    pub angle: f64,
    /// Angular velocity \[rad/s\].
    pub angular_velocity: f64,
    /// Linear velocity \[vx, vy\] \[m/s\].
    pub linear_velocity: [f64; 2],
}
/// Per-link intermediate data for the ABA forward pass.
#[derive(Debug, Clone, Default)]
pub struct AbaLinkData {
    /// Articulated-body inertia (scalar, 1-D approximation).
    pub abi: f64,
    /// Bias force (scalar).
    pub bias: f64,
    /// Intermediate velocity propagation term.
    pub v: f64,
    /// Computed joint acceleration.
    pub qdd: f64,
}
/// A 6×6 spatial inertia matrix stored as two 3×3 blocks for efficiency.
///
/// The matrix has the form:
/// ```text
/// Ia = [ I_c + m*c× * c×ᵀ   m*c× ]
///      [ m*c×ᵀ               m*I₃ ]
/// ```
/// where `c×` is the skew-symmetric cross-product matrix of the CoM offset.
#[derive(Debug, Clone, Copy, Default)]
pub struct SpatialInertia6 {
    /// Link mass \[kg\].
    pub mass: f64,
    /// Centre of mass in local frame \[m\].
    pub com: [f64; 3],
    /// Inertia tensor about CoM \[kg·m²\] (3×3, row-major).
    pub i_rot: [[f64; 3]; 3],
}
impl SpatialInertia6 {
    /// Create from physical parameters.
    pub fn new(mass: f64, com: [f64; 3], i_rot: [[f64; 3]; 3]) -> Self {
        Self { mass, com, i_rot }
    }
    /// Thin-rod inertia about one end, along x-axis.
    pub fn thin_rod_x(mass: f64, length: f64) -> Self {
        let i_perp = mass * length * length / 3.0;
        Self::new(
            mass,
            [length / 2.0, 0.0, 0.0],
            [[0.0, 0.0, 0.0], [0.0, i_perp, 0.0], [0.0, 0.0, i_perp]],
        )
    }
    /// Multiply the spatial inertia by a 6-D spatial acceleration vector.
    /// Convert to the `SpatialInertia` type used in the core RNEA code.
    fn to_core_spatial_inertia(self) -> SpatialInertia {
        SpatialInertia::new(self.mass, self.com, self.i_rot)
    }
    /// Multiply this inertia matrix by a spatial vector.
    ///
    /// Result is the net spatial force: `f = I_A * a`.
    pub fn mul_vec(&self, a: &SpatialVec) -> SpatialVec {
        self.to_core_spatial_inertia().multiply(a)
    }
}
/// Wrench: force and torque pair in 3-D.
#[derive(Debug, Clone, Copy, Default)]
pub struct Wrench {
    /// Force vector (N).
    pub force: [f64; 3],
    /// Torque vector (N·m).
    pub torque: [f64; 3],
}
impl Wrench {
    /// Create a wrench from torque and force.
    pub fn new(torque: [f64; 3], force: [f64; 3]) -> Self {
        Self { force, torque }
    }
    /// Create a zero wrench.
    pub fn zero() -> Self {
        Self::default()
    }
    /// Add two wrenches.
    pub fn add(&self, other: &Wrench) -> Wrench {
        Wrench {
            force: add3w(self.force, other.force),
            torque: add3w(self.torque, other.torque),
        }
    }
    /// Scale a wrench.
    pub fn scale(&self, s: f64) -> Wrench {
        Wrench {
            force: scale3w(self.force, s),
            torque: scale3w(self.torque, s),
        }
    }
}
/// A loop-closure constraint that forces two link endpoints to coincide.
#[derive(Debug, Clone)]
pub struct LoopClosureConstraint {
    /// Index of the first link in the loop.
    pub link_a: usize,
    /// Index of the second link in the loop.
    pub link_b: usize,
    /// Constraint violation (position error, m).
    pub violation: [f64; 3],
}
impl LoopClosureConstraint {
    /// Create a new constraint.
    pub fn new(link_a: usize, link_b: usize) -> Self {
        Self {
            link_a,
            link_b,
            violation: [0.0; 3],
        }
    }
    /// Check whether the constraint is satisfied within tolerance.
    pub fn is_satisfied(&self, tol: f64) -> bool {
        let v = &self.violation;
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() < tol
    }
}
/// Open kinematic chain of articulated rigid bodies.
#[derive(Debug, Clone)]
pub struct ArticulatedChain {
    /// All links, stored in insertion order.
    pub links: Vec<ArticulatedLink>,
    /// World-space transform of the root link.
    pub root_transform: Transform,
}
impl ArticulatedChain {
    /// Create an empty chain.
    pub fn new() -> Self {
        Self {
            links: Vec::new(),
            root_transform: Transform::default(),
        }
    }
    /// Add the root link (no parent, fixed to world).
    ///
    /// Returns the [`LinkId`] of the new root link.
    pub fn add_root(&mut self, mass: f64, inertia: [[f64; 3]; 3]) -> LinkId {
        let id = LinkId(self.links.len() as u32);
        self.links.push(ArticulatedLink {
            id,
            parent: None,
            children: Vec::new(),
            joint_type: JointDof::Fixed,
            joint_position: 0.0,
            joint_velocity: 0.0,
            joint_limit: None,
            local_transform: Transform::default(),
            mass,
            inertia,
        });
        id
    }
    /// Add a child link attached to `parent` via `joint`.
    ///
    /// `local_tf` is the resting pose of this link relative to the parent joint
    /// origin (applied *before* the joint displacement).
    ///
    /// Returns the [`LinkId`] of the new link.
    pub fn add_link(
        &mut self,
        parent: LinkId,
        joint: JointDof,
        local_tf: Transform,
        mass: f64,
        inertia: [[f64; 3]; 3],
    ) -> LinkId {
        let id = LinkId(self.links.len() as u32);
        self.links.push(ArticulatedLink {
            id,
            parent: Some(parent),
            children: Vec::new(),
            joint_type: joint,
            joint_position: 0.0,
            joint_velocity: 0.0,
            joint_limit: None,
            local_transform: local_tf,
            mass,
            inertia,
        });
        let parent_idx = parent.0 as usize;
        self.links[parent_idx].children.push(id);
        id
    }
    /// Set the joint position for `id`.
    pub fn set_joint_position(&mut self, id: LinkId, q: f64) {
        self.links[id.0 as usize].joint_position = q;
    }
    /// Set the joint velocity for `id`.
    pub fn set_joint_velocity(&mut self, id: LinkId, q_dot: f64) {
        self.links[id.0 as usize].joint_velocity = q_dot;
    }
    /// Compute the world-space transform of every link via forward kinematics.
    ///
    /// The walk is breadth-first from the root; each child accumulates its
    /// parent's world transform, the link's `local_transform`, and then the
    /// joint displacement.
    pub fn compute_forward_kinematics(&self) -> Vec<(LinkId, Transform)> {
        let n = self.links.len();
        let mut world_tfs: Vec<Option<Transform>> = vec![None; n];
        for link in &self.links {
            let idx = link.id.0 as usize;
            let parent_tf: Transform = match link.parent {
                None => self.root_transform.clone(),
                Some(pid) => world_tfs[pid.0 as usize]
                    .as_ref()
                    .expect("parent must be processed before child")
                    .clone(),
            };
            let after_local = parent_tf.compose(&link.local_transform);
            let joint_tf = joint_displacement(&link.joint_type, link.joint_position);
            world_tfs[idx] = Some(after_local.compose(&joint_tf));
        }
        world_tfs
            .into_iter()
            .enumerate()
            .map(|(i, tf)| (LinkId(i as u32), tf.expect("all links must be reachable")))
            .collect()
    }
    /// Compute one column of the geometric Jacobian for a revolute joint.
    ///
    /// Returns a 6-vector `[linear (3), angular (3)]`.
    ///
    /// * `joint_axis_world` – joint rotation axis expressed in world frame.
    /// * `joint_pos_world`  – position of the joint origin in world frame.
    /// * `ee_pos_world`     – position of the end-effector in world frame.
    pub fn revolute_jacobian_column(
        joint_axis_world: Vec3,
        joint_pos_world: Vec3,
        ee_pos_world: Vec3,
    ) -> [f64; 6] {
        let r = ee_pos_world - joint_pos_world;
        let linear = joint_axis_world.cross(&r);
        [
            linear.x,
            linear.y,
            linear.z,
            joint_axis_world.x,
            joint_axis_world.y,
            joint_axis_world.z,
        ]
    }
    /// Compute the 6×N geometric Jacobian for the end-effector at `ee_link`.
    ///
    /// Columns correspond to revolute joints on the path from root to
    /// `ee_link` (fixed/prismatic joints contribute a zero column; spherical
    /// and planar joints are treated as fixed for this single-DOF
    /// approximation).
    ///
    /// The result is a `Vec` of `[f64; 6]` columns, one per link on the
    /// ancestor path (including `ee_link`).
    pub fn compute_jacobian(&self, ee_link: LinkId) -> Vec<[f64; 6]> {
        let mut path: Vec<LinkId> = Vec::new();
        let mut cur = Some(ee_link);
        while let Some(id) = cur {
            path.push(id);
            cur = self.links[id.0 as usize].parent;
        }
        path.reverse();
        let world_tfs = self.compute_forward_kinematics();
        let ee_pos = world_tfs[ee_link.0 as usize].1.position;
        let mut columns: Vec<[f64; 6]> = Vec::with_capacity(path.len());
        for link_id in &path {
            let link = &self.links[link_id.0 as usize];
            let world_tf = &world_tfs[link_id.0 as usize].1;
            let col = match &link.joint_type {
                JointDof::Revolute(axis) => {
                    let axis_local = Vec3::new(axis[0], axis[1], axis[2]);
                    let axis_world = world_tf.transform_vector(&axis_local).normalize();
                    let joint_pos = world_tf.position;
                    Self::revolute_jacobian_column(axis_world, joint_pos, ee_pos)
                }
                JointDof::Prismatic(axis) => {
                    let axis_local = Vec3::new(axis[0], axis[1], axis[2]);
                    let axis_world = world_tf.transform_vector(&axis_local).normalize();
                    [axis_world.x, axis_world.y, axis_world.z, 0.0, 0.0, 0.0]
                }
                JointDof::Fixed | JointDof::Spherical | JointDof::Planar => [0.0; 6],
            };
            columns.push(col);
        }
        columns
    }
    /// Clamp each link's `joint_velocity` to its declared joint limits.
    pub fn apply_velocity_limits(&mut self) {
        for link in &mut self.links {
            if let Some((min, max)) = link.joint_limit {
                link.joint_velocity = link.joint_velocity.clamp(min, max);
            }
        }
    }
    /// Euler-integrate joint positions by `dt` and clamp to joint limits.
    pub fn integrate(&mut self, dt: Real) {
        for link in &mut self.links {
            link.joint_position += link.joint_velocity * dt;
            if let Some((min, max)) = link.joint_limit {
                link.joint_position = link.joint_position.clamp(min, max);
            }
        }
    }
}
/// A 6-D spatial vector `[ω; v]` (angular then linear) in Plücker coordinates.
#[derive(Debug, Clone, Copy, Default)]
pub struct SpatialVec {
    /// Angular (rotational) component.
    pub angular: [f64; 3],
    /// Linear (translational) component.
    pub linear: [f64; 3],
}
impl SpatialVec {
    /// Create from explicit components.
    pub fn new(angular: [f64; 3], linear: [f64; 3]) -> Self {
        Self { angular, linear }
    }
    /// Dot product of two spatial vectors: a · b = aω·bv + av·bω.
    pub fn dot(&self, other: &SpatialVec) -> f64 {
        dot3(self.angular, other.linear) + dot3(self.linear, other.angular)
    }
    /// Spatial cross product: f(v) × f(w) for spatial motion vectors.
    pub fn cross_motion(&self, other: &SpatialVec) -> SpatialVec {
        SpatialVec {
            angular: cross3(self.angular, other.angular),
            linear: [
                cross3(self.angular, other.linear)[0] + cross3(self.linear, other.angular)[0],
                cross3(self.angular, other.linear)[1] + cross3(self.linear, other.angular)[1],
                cross3(self.angular, other.linear)[2] + cross3(self.linear, other.angular)[2],
            ],
        }
    }
    /// Scale a spatial vector by a scalar.
    pub fn scale(&self, s: f64) -> SpatialVec {
        SpatialVec {
            angular: [
                self.angular[0] * s,
                self.angular[1] * s,
                self.angular[2] * s,
            ],
            linear: [self.linear[0] * s, self.linear[1] * s, self.linear[2] * s],
        }
    }
    /// Add two spatial vectors.
    pub fn add(&self, other: &SpatialVec) -> SpatialVec {
        SpatialVec {
            angular: [
                self.angular[0] + other.angular[0],
                self.angular[1] + other.angular[1],
                self.angular[2] + other.angular[2],
            ],
            linear: [
                self.linear[0] + other.linear[0],
                self.linear[1] + other.linear[1],
                self.linear[2] + other.linear[2],
            ],
        }
    }
}
/// Open 2-D kinematic chain of articulated rigid bodies.
#[derive(Debug, Clone)]
pub struct ArticulatedChain2D {
    /// All links in insertion order.
    pub links: Vec<ArticulatedLink2D>,
}
impl ArticulatedChain2D {
    /// Create an empty 2-D chain.
    pub fn new() -> Self {
        Self { links: Vec::new() }
    }
    /// Add a new revolute link to the chain.
    ///
    /// Returns the index of the newly added link.
    pub fn add_link(&mut self, mass: f64, length: f64, parent: Option<usize>) -> usize {
        let idx = self.links.len();
        self.links.push(ArticulatedLink2D {
            body: LinkBody {
                mass,
                length,
                inertia: mass * length * length / 3.0,
                position: [0.0; 2],
                angle: 0.0,
                angular_velocity: 0.0,
                linear_velocity: [0.0; 2],
            },
            joint_type: JointType::Revolute,
            joint_angle: 0.0,
            joint_velocity: 0.0,
            joint_torque: 0.0,
            parent,
        });
        idx
    }
    /// Return the number of links in the chain.
    pub fn link_count(&self) -> usize {
        self.links.len()
    }
    /// Forward kinematics: returns `(end_pos, global_angle)` for each link's tip.
    ///
    /// `base_pos` and `base_angle` define the chain root in world space.
    pub fn forward_kinematics(&self, base_pos: [f64; 2], base_angle: f64) -> Vec<([f64; 2], f64)> {
        let mut results: Vec<([f64; 2], f64)> = Vec::with_capacity(self.links.len());
        let mut origins: Vec<[f64; 2]> = Vec::with_capacity(self.links.len());
        let mut angles: Vec<f64> = Vec::with_capacity(self.links.len());
        for (i, link) in self.links.iter().enumerate() {
            let (ox, oy, parent_angle): (f64, f64, f64) = match link.parent {
                None => (base_pos[0], base_pos[1], base_angle),
                Some(p) => {
                    let tip = results[p].0;
                    let pa = angles[p];
                    (tip[0], tip[1], pa)
                }
            };
            let global_angle: f64 = parent_angle + link.joint_angle;
            let tip_x = ox + link.body.length * global_angle.cos();
            let tip_y = oy + link.body.length * global_angle.sin();
            results.push(([tip_x, tip_y], global_angle));
            origins.push([ox, oy]);
            angles.push(global_angle);
            let _ = i;
        }
        results
    }
    /// Geometric Jacobian: partial derivative of the end-effector position
    /// w.r.t. each joint angle (2 × n_joints, simplified).
    ///
    /// Returns one `[dx/dθ_i, dy/dθ_i]` per joint (link) for the last link's tip.
    pub fn jacobian(&self, base_pos: [f64; 2], base_angle: f64) -> Vec<[f64; 2]> {
        let n = self.links.len();
        if n == 0 {
            return Vec::new();
        }
        let fk = self.forward_kinematics(base_pos, base_angle);
        let ee = fk[n - 1].0;
        let mut joint_origins: Vec<[f64; 2]> = Vec::with_capacity(n);
        {
            let mut origins: Vec<[f64; 2]> = Vec::with_capacity(n);
            for (i, link) in self.links.iter().enumerate() {
                let origin = match link.parent {
                    None => base_pos,
                    Some(p) => origins[p],
                };
                let angle = fk[i].1;
                let tip = [
                    origin[0] + link.body.length * angle.cos(),
                    origin[1] + link.body.length * angle.sin(),
                ];
                joint_origins.push(origin);
                origins.push(tip);
            }
        }
        let cols: Vec<[f64; 2]> = joint_origins
            .iter()
            .map(|jo| {
                let dx = ee[0] - jo[0];
                let dy = ee[1] - jo[1];
                [-dy, dx]
            })
            .collect();
        cols
    }
    /// Newton-Euler inverse dynamics: given joint accelerations and gravity,
    /// returns the required joint torques.
    ///
    /// Uses a simplified scalar approach for a planar chain.
    pub fn inverse_dynamics(&self, joint_accelerations: &[f64], gravity: [f64; 2]) -> Vec<f64> {
        let n = self.links.len();
        let mut torques = vec![0.0_f64; n];
        for (i, torque) in torques.iter_mut().enumerate() {
            let mut tau = 0.0;
            for (j, &ja_j) in joint_accelerations.iter().enumerate().skip(i).take(n - i) {
                tau += self.links[j].body.inertia * ja_j;
                let angle_j = {
                    let mut a = 0.0;
                    for k in 0..=j {
                        a += self.links[k].joint_angle;
                    }
                    a
                };
                let g_perp = -gravity[0] * angle_j.sin() + gravity[1] * angle_j.cos();
                let r = {
                    let mut dist = 0.0;
                    for k in (i + 1)..j {
                        dist += self.links[k].body.length;
                    }
                    if j > i {
                        dist += self.links[j].body.length * 0.5;
                    } else {
                        dist = self.links[j].body.length * 0.5;
                    }
                    dist
                };
                tau += self.links[j].body.mass * g_perp * r;
            }
            *torque = tau;
        }
        torques
    }
    /// Integrate joint velocities and angles by `dt` using applied torques.
    ///
    /// Uses semi-implicit Euler: ω_new = ω + (τ/I)*dt, θ_new = θ + ω_new*dt.
    pub fn apply_torques(&mut self, torques: &[f64], dt: f64) {
        for (i, link) in self.links.iter_mut().enumerate() {
            if i < torques.len() {
                let alpha = torques[i] / link.body.inertia.max(1e-12);
                link.joint_velocity += alpha * dt;
                link.joint_angle += link.joint_velocity * dt;
            }
        }
    }
    /// Compute the total rotational kinetic energy of the chain.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.links
            .iter()
            .map(|l| 0.5 * l.body.inertia * l.joint_velocity * l.joint_velocity)
            .sum()
    }
    /// Compute the gravitational potential energy of all links.
    ///
    /// `gravity_magnitude` is |g| (positive), height is the y-coordinate of each link's COM.
    pub fn potential_energy(
        &self,
        base_pos: [f64; 2],
        base_angle: f64,
        gravity_magnitude: f64,
    ) -> f64 {
        let fk = self.forward_kinematics(base_pos, base_angle);
        let mut energy = 0.0;
        let mut origins: Vec<[f64; 2]> = Vec::with_capacity(self.links.len());
        {
            let mut tips: Vec<[f64; 2]> = Vec::with_capacity(self.links.len());
            for (i, link) in self.links.iter().enumerate() {
                let origin = match link.parent {
                    None => base_pos,
                    Some(p) => tips[p],
                };
                let tip = fk[i].0;
                let com_y = (origin[1] + tip[1]) * 0.5;
                energy += link.body.mass * gravity_magnitude * com_y;
                origins.push(origin);
                tips.push(tip);
            }
        }
        let _ = origins;
        energy
    }
}
/// Joint-space impedance controller: τ = K_d(q_d - q) + B_d(qd_d - qd).
#[derive(Debug, Clone)]
pub struct JointImpedanceController {
    /// Desired stiffness (Nm/rad).
    pub stiffness: f64,
    /// Desired damping (Nms/rad).
    pub damping: f64,
}
impl JointImpedanceController {
    /// Create a new impedance controller.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }
    /// Compute the joint torque for the given position and velocity errors.
    pub fn torque(&self, pos_error: f64, vel_error: f64) -> f64 {
        self.stiffness * pos_error + self.damping * vel_error
    }
    /// Compute torques for a whole joint-space vector.
    pub fn torques(&self, pos_errors: &[f64], vel_errors: &[f64]) -> Vec<f64> {
        pos_errors
            .iter()
            .zip(vel_errors.iter())
            .map(|(&pe, &ve)| self.torque(pe, ve))
            .collect()
    }
}
/// Opaque identifier for a link in an articulated chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LinkId(pub u32);
/// Per-link RNEA data.
#[derive(Debug, Clone, Default)]
pub struct RneaLink {
    /// Spatial velocity of the link.
    pub velocity: SpatialVec,
    /// Spatial acceleration of the link.
    pub acceleration: SpatialVec,
    /// Net spatial force on the link.
    pub force: SpatialVec,
}
/// Spatial rigid-body inertia matrix (6×6 symmetric) in body coordinates.
///
/// Stored in the compact form as per Featherstone:
/// ```text
/// I = [ I_c + m r×r×ᵀ  m r× ]
///     [ m r×ᵀ           m I₃ ]
/// ```
/// Here we store only mass `m`, CoM offset `c`, and rotational inertia `I_rot`.
#[derive(Debug, Clone, Copy)]
pub struct SpatialInertia {
    /// Body mass (kg).
    pub mass: f64,
    /// Centre of mass position in body frame.
    pub com: [f64; 3],
    /// 3×3 rotational inertia tensor about CoM.
    pub i_rot: [[f64; 3]; 3],
}
impl SpatialInertia {
    /// Create from mass, CoM, and rotational inertia.
    pub fn new(mass: f64, com: [f64; 3], i_rot: [[f64; 3]; 3]) -> Self {
        Self { mass, com, i_rot }
    }
    /// Create a uniform sphere inertia.
    pub fn sphere(mass: f64, radius: f64, com: [f64; 3]) -> Self {
        let i = 2.0 / 5.0 * mass * radius * radius;
        Self::new(mass, com, [[i, 0.0, 0.0], [0.0, i, 0.0], [0.0, 0.0, i]])
    }
    /// Create a box inertia at the origin.
    pub fn box_inertia(mass: f64, lx: f64, ly: f64, lz: f64) -> Self {
        let ixx = mass / 12.0 * (ly * ly + lz * lz);
        let iyy = mass / 12.0 * (lx * lx + lz * lz);
        let izz = mass / 12.0 * (lx * lx + ly * ly);
        Self::new(
            mass,
            [0.0; 3],
            [[ixx, 0.0, 0.0], [0.0, iyy, 0.0], [0.0, 0.0, izz]],
        )
    }
    /// Create a sphere inertia at the origin.
    pub fn sphere_inertia(mass: f64, radius: f64) -> Self {
        Self::sphere(mass, radius, [0.0; 3])
    }
    /// Multiply by a spatial acceleration vector: `I × a`.
    pub fn multiply(&self, a: &SpatialVec) -> SpatialVec {
        let alpha = a.angular;
        let a_lin = a.linear;
        let alpha_cross_c = cross3(alpha, self.com);
        let v_eff = add3(a_lin, alpha_cross_c);
        let lin = scale3(v_eff, self.mass);
        let i_alpha = mat3_mul_vec(self.i_rot, alpha);
        let c_cross_veff = cross3(self.com, v_eff);
        let angular = add3(i_alpha, scale3(c_cross_veff, self.mass));
        SpatialVec {
            angular,
            linear: lin,
        }
    }
}
/// A planar joint: 2 translational + 1 rotational DOF in the XY-plane.
#[derive(Debug, Clone, Default)]
pub struct PlanarJoint {
    /// Translation in x (m).
    pub tx: f64,
    /// Translation in y (m).
    pub ty: f64,
    /// Rotation about z (rad).
    pub rz: f64,
}
impl PlanarJoint {
    /// Create a zero-configuration planar joint.
    pub fn new() -> Self {
        Self::default()
    }
    /// Apply generalised velocities `[vx, vy, ωz]` for time step `dt`.
    pub fn integrate(&mut self, vx: f64, vy: f64, omega_z: f64, dt: f64) {
        self.tx += vx * dt;
        self.ty += vy * dt;
        self.rz += omega_z * dt;
    }
    /// Pose as a 3×3 homogeneous matrix `[[r00, r01, tx], [r10, r11, ty], [0, 0, 1]]`.
    pub fn to_homogeneous(&self) -> [[f64; 3]; 3] {
        let c = self.rz.cos();
        let s = self.rz.sin();
        [[c, -s, self.tx], [s, c, self.ty], [0.0, 0.0, 1.0]]
    }
}
