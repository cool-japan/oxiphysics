//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Quat, Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use super::functions::{JOINT_BAUMGARTE, apply_pair_impulse, read_body};

/// Ball joint (spherical joint): constrains two anchor points to coincide.
///
/// Allows free rotation but prevents relative translation at the joint point.
#[derive(Debug, Clone)]
pub struct BallJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl BallJoint {
    /// Create a new ball joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
}
/// Prismatic joint: allows translation along a single axis.
///
/// Constrains bodies so they can only slide along a specified axis.
/// Optionally enforces distance limits.
#[derive(Debug, Clone)]
pub struct PrismaticJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Sliding axis in body A's local space.
    pub local_axis: Vec3,
    /// Optional lower distance limit.
    pub lower_limit: Option<Real>,
    /// Optional upper distance limit.
    pub upper_limit: Option<Real>,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) world_axis: Vec3,
    pub(super) perp1: Vec3,
    pub(super) perp2: Vec3,
}
impl PrismaticJoint {
    /// Create a prismatic joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis: Vec3,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis: local_axis.normalize(),
            lower_limit: None,
            upper_limit: None,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            world_axis: Vec3::zeros(),
            perp1: Vec3::zeros(),
            perp2: Vec3::zeros(),
        }
    }
    /// Set distance limits.
    pub fn with_limits(mut self, lower: Real, upper: Real) -> Self {
        self.lower_limit = Some(lower);
        self.upper_limit = Some(upper);
        self
    }
    /// Compute the lead-screw (helix) constraint residual for a screw joint.
    ///
    /// A lead-screw joint couples linear translation along the slide axis to
    /// rotation around the same axis with pitch `lead` (metres per radian).
    ///
    /// The velocity-level constraint is:
    ///   `v_rel · axis  - lead * ω_rel · axis = 0`
    ///
    /// # Parameters
    /// * `v_rel_along_axis` – relative linear velocity projected onto the slide axis (m/s).
    /// * `omega_rel_along_axis` – relative angular velocity projected onto the slide axis (rad/s).
    /// * `lead` – thread lead in metres/radian (positive for right-handed thread).
    ///
    /// # Returns
    /// The constraint violation (velocity residual).
    pub fn compute_lead_screw(
        &self,
        v_rel_along_axis: Real,
        omega_rel_along_axis: Real,
        lead: Real,
    ) -> Real {
        v_rel_along_axis - lead * omega_rel_along_axis
    }
    /// Apply a lead-screw impulse to enforce the screw coupling.
    ///
    /// Computes and applies the impulse that drives the lead-screw constraint to zero.
    /// Returns the impulse magnitude.
    ///
    /// # Parameters
    /// * `bodies`  – the rigid body set.
    /// * `lead`    – thread lead in metres/radian (positive → right-handed).
    pub fn apply_lead_screw_impulse(&self, bodies: &mut RigidBodySet, lead: Real) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let axis = self.world_axis;
        let v_a = a.vel.dot(&axis);
        let v_b = b.vel.dot(&axis);
        let omega_a = a.ang_vel.dot(&axis);
        let omega_b = b.ang_vel.dot(&axis);
        let v_rel = v_a - v_b;
        let omega_rel = omega_a - omega_b;
        let c_vel = v_rel - lead * omega_rel;
        let k_lin = a.inv_mass + b.inv_mass;
        let k_ang =
            lead * lead * (axis.dot(&(a.inv_inertia * axis)) + axis.dot(&(b.inv_inertia * axis)));
        let k = k_lin + k_ang;
        if k < 1e-12 {
            return 0.0;
        }
        let lambda = -c_vel / k;
        let lin_impulse = axis * lambda;
        let ang_impulse = axis * (lead * lambda);
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.velocity += lin_impulse * body.inverse_mass;
            body.angular_velocity += body.world_inverse_inertia * ang_impulse;
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.velocity -= lin_impulse * body.inverse_mass;
            body.angular_velocity -= body.world_inverse_inertia * ang_impulse;
        }
        lambda
    }
}
/// Revolute joint: allows rotation around a single axis.
///
/// Constrains the anchor points to coincide and restricts rotation to a
/// single axis. Optionally enforces angle limits.
#[derive(Debug, Clone)]
pub struct RevoluteJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Rotation axis in body A's local space.
    pub local_axis: Vec3,
    /// Optional lower angle limit (radians).
    pub lower_limit: Option<Real>,
    /// Optional upper angle limit (radians).
    pub upper_limit: Option<Real>,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) world_axis: Vec3,
}
impl RevoluteJoint {
    /// Create a revolute joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis: Vec3,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis: local_axis.normalize(),
            lower_limit: None,
            upper_limit: None,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            world_axis: Vec3::zeros(),
        }
    }
    /// Set angle limits.
    pub fn with_limits(mut self, lower: Real, upper: Real) -> Self {
        self.lower_limit = Some(lower);
        self.upper_limit = Some(upper);
        self
    }
    /// Compute a gear-ratio velocity constraint between this joint and another revolute joint.
    ///
    /// The constraint is: `ω_a · axis_a + gear_ratio * ω_b_other · axis_b_other = 0`.
    ///
    /// Returns the constraint violation (velocity-level residual).  A value of
    /// zero means the gear constraint is satisfied.
    ///
    /// # Parameters
    /// * `omega_a` – angular velocity of this joint's body A along its world axis.
    /// * `omega_b_other` – angular velocity along the other joint's world axis.
    /// * `gear_ratio` – tooth ratio N_other / N_self.  Negative for external gears.
    pub fn compute_gear_ratio_constraint(
        &self,
        omega_a: Real,
        omega_b_other: Real,
        gear_ratio: Real,
    ) -> Real {
        omega_a + gear_ratio * omega_b_other
    }
    /// Compute and apply the gear-ratio impulse between two revolute joints.
    ///
    /// Solves for `λ` such that after applying it the gear constraint is satisfied
    /// to first order:  `ω_a_proj + gear_ratio * ω_b_proj = 0`.
    ///
    /// Returns the impulse magnitude `λ` applied.
    ///
    /// # Parameters
    /// * `bodies` – the rigid body set.
    /// * `other_body` – the handle of the body of the second revolute joint.
    /// * `other_world_axis` – the world-space axis of the second revolute joint.
    /// * `gear_ratio` – the gear ratio `N_other / N_self`.
    pub fn apply_gear_ratio_impulse(
        &self,
        bodies: &mut RigidBodySet,
        other_body: BodyHandle,
        other_world_axis: Vec3,
        gear_ratio: Real,
    ) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, other_body) {
            Some(p) => p,
            None => return 0.0,
        };
        let omega_a = a.ang_vel.dot(&self.world_axis);
        let omega_b = b.ang_vel.dot(&other_world_axis);
        let c_vel = omega_a + gear_ratio * omega_b;
        let ka = self.world_axis.dot(&(a.inv_inertia * self.world_axis));
        let kb = other_world_axis.dot(&(b.inv_inertia * other_world_axis));
        let k = ka + gear_ratio * gear_ratio * kb;
        if k < 1e-12 {
            return 0.0;
        }
        let lambda = -c_vel / k;
        let impulse_a = self.world_axis * lambda;
        let impulse_b = other_world_axis * (gear_ratio * lambda);
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.angular_velocity += body.world_inverse_inertia * impulse_a;
        }
        if let Some(body) = bodies.get_mut(other_body) {
            body.angular_velocity += body.world_inverse_inertia * impulse_b;
        }
        lambda
    }
}
/// Rack-and-pinion joint: converts rotation of a pinion gear into linear
/// translation of a rack.
///
/// The constraint is: v_rack · rack_axis - pitch_radius * (ω_pinion · pinion_axis) = 0
#[derive(Debug, Clone)]
pub struct RackPinionJoint {
    /// Handle of the rotating pinion body.
    pub body_pinion: u32,
    /// Handle of the translating rack body.
    pub body_rack: u32,
    /// Pitch radius of the pinion (m).
    pub pitch_radius: f64,
    /// Rotation axis of the pinion (local frame).
    pub pinion_axis: Vec3,
    /// Translation axis of the rack (local frame).
    pub rack_axis: Vec3,
    /// Accumulated constraint impulse.
    pub accumulated_lambda: f64,
}
impl RackPinionJoint {
    /// Create a new rack-and-pinion joint.
    pub fn new(
        body_pinion: u32,
        body_rack: u32,
        pitch_radius: f64,
        pinion_axis: Vec3,
        rack_axis: Vec3,
    ) -> Self {
        Self {
            body_pinion,
            body_rack,
            pitch_radius,
            pinion_axis,
            rack_axis,
            accumulated_lambda: 0.0,
        }
    }
    /// Velocity constraint: v_rack · rack_axis - pitch_radius * (ω_pinion · pinion_axis) = 0
    ///
    /// A positive pinion rotation (ω_pinion · pinion_axis > 0) drives positive rack
    /// translation (v_rack · rack_axis > 0).
    pub fn constraint_velocity(&self, omega_pinion: Vec3, v_rack: Vec3) -> f64 {
        v_rack.dot(&self.rack_axis) - self.pitch_radius * omega_pinion.dot(&self.pinion_axis)
    }
    /// Solve the velocity constraint and return the impulse magnitude.
    ///
    /// `inv_inertia` is the scalar inverse inertia of the pinion projected onto
    /// `pinion_axis`; `inv_mass` is the scalar inverse mass of the rack.
    pub fn solve(
        &mut self,
        omega_pinion: Vec3,
        v_rack: Vec3,
        inv_inertia: f64,
        inv_mass: f64,
    ) -> f64 {
        let c_vel = self.constraint_velocity(omega_pinion, v_rack);
        let j_pinion_mag_sq =
            self.pitch_radius * self.pitch_radius * self.pinion_axis.norm_squared();
        let j_rack_mag_sq = self.rack_axis.norm_squared();
        let k = j_pinion_mag_sq * inv_inertia + j_rack_mag_sq * inv_mass;
        if k < 1e-12 {
            return 0.0;
        }
        let lambda = -c_vel / k;
        self.accumulated_lambda += lambda;
        lambda
    }
}
/// Cone (angular) joint: allows rotation within a cone around a reference axis.
///
/// Constrains the angle between the joint axis and a reference direction to be
/// within `half_angle`.  This is used for shoulder joints, hip joints, and
/// other joints with angular swing limits.
///
/// The positional constraint is the same as a ball joint.  Only the swing
/// (cone) is constrained angularly; twist around the axis is unconstrained.
#[derive(Debug, Clone)]
pub struct ConeJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Reference axis in body A's local space.
    pub local_axis: Vec3,
    /// Half angle of the cone (radians).
    pub half_angle: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) world_axis: Vec3,
}
impl ConeJoint {
    /// Create a new cone joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis: Vec3,
        half_angle: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis: local_axis.normalize(),
            half_angle: half_angle.clamp(0.0, std::f64::consts::PI),
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            world_axis: Vec3::zeros(),
        }
    }
    /// Compute the current swing angle between body A's axis and body B's corresponding axis.
    ///
    /// Returns the angle in radians \[0, π\].
    pub fn current_swing_angle(&self, bodies: &RigidBodySet) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let axis_a = (a.rotation * self.local_axis).normalize();
        let axis_b = (b.rotation * self.local_axis).normalize();
        let cos_a = axis_a.dot(&axis_b).clamp(-1.0, 1.0);
        cos_a.acos()
    }
    /// Returns `true` if the joint is currently within the cone limit.
    pub fn within_limit(&self, bodies: &RigidBodySet) -> bool {
        self.current_swing_angle(bodies) <= self.half_angle + 1e-6
    }
}
/// Weld joint: rigidly locks two bodies together at an offset.
///
/// Unlike the `FixedJoint` which uses arbitrary local anchors, the weld joint
/// stores the *desired* relative pose (position offset + rotation) and enforces
/// it. Useful for dynamically welding two bodies at runtime.
#[derive(Debug, Clone)]
pub struct WeldJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Desired offset from body A's origin to body B's origin (in world space
    /// at the time of creation).
    pub desired_offset: Vec3,
    /// Desired relative rotation (B relative to A).
    pub desired_relative_rotation: Quat,
    /// Baumgarte correction factor.
    pub baumgarte: Real,
    pub(super) eff_mass_x: Real,
    pub(super) eff_mass_y: Real,
    pub(super) eff_mass_z: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl WeldJoint {
    /// Create a new weld joint that freezes the current relative pose.
    ///
    /// Call this when body A is at `pos_a` with rotation `rot_a` and body B
    /// is at `pos_b` with rotation `rot_b`.
    pub fn from_current_pose(
        body_a: BodyHandle,
        body_b: BodyHandle,
        pos_a: Vec3,
        pos_b: Vec3,
        rot_a: Quat,
        rot_b: Quat,
    ) -> Self {
        let desired_offset = pos_b - pos_a;
        let desired_relative_rotation = rot_a.inverse() * rot_b;
        Self {
            body_a,
            body_b,
            desired_offset,
            desired_relative_rotation,
            baumgarte: JOINT_BAUMGARTE,
            eff_mass_x: 0.0,
            eff_mass_y: 0.0,
            eff_mass_z: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Create a weld joint with explicit desired offset and relative rotation.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        desired_offset: Vec3,
        desired_relative_rotation: Quat,
    ) -> Self {
        Self {
            body_a,
            body_b,
            desired_offset,
            desired_relative_rotation,
            baumgarte: JOINT_BAUMGARTE,
            eff_mass_x: 0.0,
            eff_mass_y: 0.0,
            eff_mass_z: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
}
/// Read body properties needed for constraint solving.
pub(super) struct BodyProps {
    pub(super) inv_mass: Real,
    pub(super) inv_inertia: oxiphysics_core::math::Mat3,
    pub(super) vel: Vec3,
    pub(super) ang_vel: Vec3,
    pub(super) position: Vec3,
    pub(super) rotation: Quat,
}
/// Motor joint: drives a body to reach a target angular velocity around an axis.
///
/// Applies angular impulses to body A (and reaction to body B) to reach the
/// target angular velocity `target_velocity` along `local_axis` on body A.
/// The maximum torque is clamped by `max_torque`.
#[derive(Debug, Clone)]
pub struct MotorJoint {
    /// Handle of body A (the driven body).
    pub body_a: BodyHandle,
    /// Handle of body B (the reaction body, often static).
    pub body_b: BodyHandle,
    /// Motor axis in body A's local space.
    pub local_axis: Vec3,
    /// Target angular velocity (rad/s) along the axis.
    pub target_velocity: Real,
    /// Maximum torque the motor can exert (N*m).
    pub max_torque: Real,
    pub(super) world_axis: Vec3,
}
impl MotorJoint {
    /// Create a new motor joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_axis: Vec3,
        target_velocity: Real,
        max_torque: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_axis: local_axis.normalize(),
            target_velocity,
            max_torque,
            world_axis: Vec3::zeros(),
        }
    }
}
/// Fixed joint: locks the relative transform between two bodies.
///
/// Constrains both position and orientation so that the two anchor frames
/// maintain a constant relative pose.
#[derive(Debug, Clone)]
pub struct FixedJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Target relative rotation (B relative to A).
    pub target_rotation: Quat,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) eff_mass_x: Real,
    pub(super) eff_mass_y: Real,
    pub(super) eff_mass_z: Real,
}
impl FixedJoint {
    /// Create a fixed joint between two bodies at the given local anchors.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        target_rotation: Quat,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            target_rotation,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            eff_mass_x: 0.0,
            eff_mass_y: 0.0,
            eff_mass_z: 0.0,
        }
    }
}
/// Distance joint: constrains two bodies to stay within \[min_distance, max_distance\]
/// from each other (measured anchor-to-anchor).
///
/// When the distance is within the allowed range, no constraint is active.
/// When the distance drops below `min_distance` a repulsive impulse is applied;
/// when it exceeds `max_distance` an attractive impulse is applied.
#[derive(Debug, Clone)]
pub struct DistanceJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Minimum allowed distance between anchors (≥ 0).
    pub min_distance: Real,
    /// Maximum allowed distance between anchors (≥ min_distance).
    pub max_distance: Real,
    /// Accumulated impulse for warm-starting.
    pub accumulated_lambda: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) eff_mass: Real,
    pub(super) constraint_dir: Vec3,
}
impl DistanceJoint {
    /// Create a new distance joint with given local anchors and distance bounds.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        min_distance: Real,
        max_distance: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            min_distance: min_distance.max(0.0),
            max_distance: max_distance.max(min_distance),
            accumulated_lambda: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            eff_mass: 0.0,
            constraint_dir: Vec3::new(1.0, 0.0, 0.0),
        }
    }
    /// Create a distance joint that locks the current distance (min == max).
    pub fn fixed_distance(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        distance: Real,
    ) -> Self {
        Self::new(
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            distance,
            distance,
        )
    }
    /// Current anchor-to-anchor distance.
    pub fn current_distance(&self, bodies: &RigidBodySet) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let world_a = a.position + a.rotation * self.local_anchor_a;
        let world_b = b.position + b.rotation * self.local_anchor_b;
        (world_b - world_a).norm()
    }
    /// Whether the current distance is within the allowed range.
    pub fn is_within_limits(&self, bodies: &RigidBodySet) -> bool {
        let d = self.current_distance(bodies);
        d >= self.min_distance - 1e-6 && d <= self.max_distance + 1e-6
    }
}
/// Extended Position-Based Dynamics (XPBD) spring joint.
///
/// Implements a soft positional constraint with a compliance parameter `α`
/// (inverse stiffness, m/N) and damping `γ` (N·s/m).  At `α = 0` the joint
/// is rigid; increasing `α` makes it more compliant (spring-like).
///
/// The XPBD formulation directly projects positions, making it
/// substep-stable even for large compliance values.
#[derive(Debug, Clone)]
pub struct XpbdJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Rest distance between anchors.
    pub rest_length: Real,
    /// Compliance α (m/N). 0.0 = rigid.
    pub compliance: Real,
    /// Damping coefficient γ.
    pub damping: Real,
    /// Accumulated XPBD Lagrange multiplier (reset each substep).
    pub lambda: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl XpbdJoint {
    /// Create a new XPBD joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        rest_length: Real,
        compliance: Real,
        damping: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            rest_length,
            compliance,
            damping,
            lambda: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Reset the accumulated Lagrange multiplier at the start of a new substep.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}
/// Pulley joint: two bodies connected by an inextensible rope over a fixed pulley.
///
/// The position-level constraint is: la + ratio * lb <= max_length
/// where la = |anchor_a - pulley_center| and lb = |anchor_b - pulley_center|.
#[derive(Debug, Clone)]
pub struct PulleyJoint {
    /// Handle of body A.
    pub body_a: u32,
    /// Handle of body B.
    pub body_b: u32,
    /// Attachment point on body A (world space).
    pub anchor_a: Vec3,
    /// Attachment point on body B (world space).
    pub anchor_b: Vec3,
    /// Fixed pulley location (world space).
    pub pulley_center: Vec3,
    /// Block-and-tackle ratio (1.0 = simple pulley).
    pub ratio: f64,
    /// Total rope length (la + ratio * lb <= max_length).
    pub max_length: f64,
    /// Accumulated constraint impulse.
    pub accumulated_lambda: f64,
}
impl PulleyJoint {
    /// Create a new pulley joint.
    pub fn new(
        body_a: u32,
        body_b: u32,
        anchor_a: Vec3,
        anchor_b: Vec3,
        pulley_center: Vec3,
        ratio: f64,
        max_length: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            pulley_center,
            ratio,
            max_length,
            accumulated_lambda: 0.0,
        }
    }
    /// Compute current rope segment lengths.
    ///
    /// Returns `(la, lb)` where `la = |anchor_a - pulley_center|` and
    /// `lb = |anchor_b - pulley_center|`.
    pub fn current_lengths(&self) -> (f64, f64) {
        let la = (self.anchor_a - self.pulley_center).norm();
        let lb = (self.anchor_b - self.pulley_center).norm();
        (la, lb)
    }
    /// Position-level constraint value: la + ratio * lb - max_length.
    ///
    /// Negative means slack; zero or positive means the rope is taut (constraint active).
    pub fn constraint_position(&self) -> f64 {
        let (la, lb) = self.current_lengths();
        la + self.ratio * lb - self.max_length
    }
    /// Unit vectors pointing from each anchor toward the pulley center.
    ///
    /// Returns `(u_a, u_b)` where `u_a = (pulley - anchor_a) / la` and
    /// `u_b = (pulley - anchor_b) / lb`.  If a segment has zero length the
    /// zero vector is returned for that element.
    pub fn unit_vectors(&self) -> (Vec3, Vec3) {
        let da = self.pulley_center - self.anchor_a;
        let la = da.norm();
        let u_a = if la > 1e-12 { da / la } else { Vec3::zeros() };
        let db = self.pulley_center - self.anchor_b;
        let lb = db.norm();
        let u_b = if lb > 1e-12 { db / lb } else { Vec3::zeros() };
        (u_a, u_b)
    }
}
/// Gear joint: constrains angular velocity ratio between two bodies.
///
/// The constraint is: ω_a · axis_a + ratio * ω_b · axis_b = 0
/// For external gears, ratio is negative (e.g. -N_b/N_a).
/// For internal gears, ratio is positive.
#[derive(Debug, Clone)]
pub struct GearJoint {
    /// Handle of body A (driving gear).
    pub body_a: u32,
    /// Handle of body B (driven gear).
    pub body_b: u32,
    /// Tooth ratio N_b / N_a (negative for external gears).
    pub ratio: f64,
    /// Rotation axis on body A (local frame).
    pub axis_a: Vec3,
    /// Rotation axis on body B (local frame).
    pub axis_b: Vec3,
    /// Accumulated constraint impulse (for warm-starting).
    pub accumulated_lambda: f64,
}
impl GearJoint {
    /// Create a new gear joint.
    pub fn new(body_a: u32, body_b: u32, ratio: f64, axis_a: Vec3, axis_b: Vec3) -> Self {
        Self {
            body_a,
            body_b,
            ratio,
            axis_a,
            axis_b,
            accumulated_lambda: 0.0,
        }
    }
    /// Velocity-level constraint: ω_a · axis_a + ratio * ω_b · axis_b = 0
    pub fn constraint_velocity(&self, omega_a: Vec3, omega_b: Vec3) -> f64 {
        omega_a.dot(&self.axis_a) + self.ratio * omega_b.dot(&self.axis_b)
    }
    /// Jacobian rows for body_a and body_b angular parts.
    ///
    /// Returns (J_a, J_b) where the constraint is J_a · ω_a + J_b · ω_b = 0.
    pub fn jacobian(&self) -> (Vec3, Vec3) {
        (self.axis_a, self.axis_b * self.ratio)
    }
    /// Solve velocity constraint and return the impulse magnitude.
    ///
    /// `inv_inertia_a` and `inv_inertia_b` are scalar inverse inertias projected
    /// onto the respective constraint axes (i.e. axis^T * I^{-1} * axis).
    pub fn solve(
        &mut self,
        omega_a: Vec3,
        omega_b: Vec3,
        inv_inertia_a: f64,
        inv_inertia_b: f64,
    ) -> f64 {
        let c_vel = self.constraint_velocity(omega_a, omega_b);
        let (j_a, j_b) = self.jacobian();
        let k = j_a.norm_squared() * inv_inertia_a + j_b.norm_squared() * inv_inertia_b;
        if k < 1e-12 {
            return 0.0;
        }
        let lambda = -c_vel / k;
        self.accumulated_lambda += lambda;
        lambda
    }
}
/// Cable (inextensible rope) joint: prevents two bodies from moving further apart
/// than `max_length`, measured anchor-to-anchor.
///
/// A cable only pulls — it goes slack when the distance is below `max_length`.
/// The tension impulse is clamped to be non-negative (attractive only).
#[derive(Debug, Clone)]
pub struct CableJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Maximum cable length (the rope goes taut at this distance).
    pub max_length: Real,
    /// Accumulated constraint impulse (for warm-starting).
    pub accumulated_lambda: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) eff_mass: Real,
    pub(super) cable_dir: Vec3,
    pub(super) taut: bool,
}
impl CableJoint {
    /// Create a new cable joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        max_length: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            max_length: max_length.max(0.0),
            accumulated_lambda: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            eff_mass: 0.0,
            cable_dir: Vec3::new(1.0, 0.0, 0.0),
            taut: false,
        }
    }
    /// Returns `true` if the cable is currently taut (distance ≥ max_length).
    pub fn is_taut(&self) -> bool {
        self.taut
    }
    /// Current anchor-to-anchor distance.
    pub fn current_length(&self, bodies: &RigidBodySet) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let world_a = a.position + a.rotation * self.local_anchor_a;
        let world_b = b.position + b.rotation * self.local_anchor_b;
        (world_b - world_a).norm()
    }
    /// Apply the tension constraint for a single velocity-level iteration.
    ///
    /// `cable_dir` points from body A to body B.  A positive
    /// `rel_v = (v_b - v_a) · cable_dir` means the bodies are separating —
    /// the cable must resist by applying an attractive impulse along `–cable_dir`
    /// to body B (and along `+cable_dir` to body A via Newton's third law).
    pub fn apply_tension_constraint(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        if !self.taut {
            return;
        }
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let va = a.vel + a.ang_vel.cross(&self.r_a);
        let vb = b.vel + b.ang_vel.cross(&self.r_b);
        let rel_v = (vb - va).dot(&self.cable_dir);
        if rel_v <= 0.0 || self.eff_mass < 1e-12 {
            return;
        }
        let lambda = rel_v * self.eff_mass;
        let impulse = self.cable_dir * lambda;
        apply_pair_impulse(
            bodies,
            self.body_a,
            self.body_b,
            impulse,
            &self.r_a,
            &self.r_b,
        );
        self.accumulated_lambda += lambda;
    }
}
/// Universal joint (Hooke's joint): allows rotation around two perpendicular axes
/// while constraining a third.
///
/// Used in drive shafts and robotic arms to transmit rotation between non-aligned
/// shafts. Constrains translation (like a ball joint) and one rotational DOF.
///
/// The constraint allows rotation around `axis_a` and `axis_b` where
/// `axis_a ⊥ axis_b` in the design configuration.  Rotation around the cross
/// product `axis_a × axis_b` is constrained.
#[derive(Debug, Clone)]
pub struct UniversalJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// First free rotation axis in body A's local space.
    pub local_axis_a: Vec3,
    /// Second free rotation axis in body B's local space.
    pub local_axis_b: Vec3,
    pub(super) world_axis_a: Vec3,
    pub(super) world_axis_b: Vec3,
    pub(super) constrained_axis: Vec3,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl UniversalJoint {
    /// Create a new universal joint.
    ///
    /// `local_axis_a` and `local_axis_b` should be approximately perpendicular
    /// at the design configuration.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis_a: Vec3,
        local_axis_b: Vec3,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis_a: local_axis_a.normalize(),
            local_axis_b: local_axis_b.normalize(),
            world_axis_a: Vec3::zeros(),
            world_axis_b: Vec3::zeros(),
            constrained_axis: Vec3::zeros(),
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Current angular constraint violation: the dot product of the two free
    /// axes should be zero (they should remain perpendicular).
    ///
    /// Returns the dot product `axis_a · axis_b` (zero when satisfied).
    pub fn angular_error(&self) -> f64 {
        self.world_axis_a.dot(&self.world_axis_b)
    }
}
/// A wrapper around any joint that breaks when the constraint force exceeds
/// a threshold.
///
/// After breaking, the joint no longer applies any impulse.  Callers should
/// periodically check [`BreakableJoint::is_broken`] and remove the constraint
/// from their list if it is broken.
#[derive(Debug, Clone)]
pub struct BreakableJoint {
    /// Whether the joint has broken.
    pub broken: bool,
    /// Maximum allowed impulse per step before the joint breaks.
    pub break_threshold: Real,
    /// Total accumulated impulse magnitude since last reset.
    pub accumulated_force: Real,
    pub(super) body_a: BodyHandle,
    pub(super) body_b: BodyHandle,
    /// Underlying spring-type constraint (encoded as spring parameters).
    pub(super) stiffness: Real,
    pub(super) natural_length: Real,
    pub(super) damping: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl BreakableJoint {
    /// Create a breakable spring-like joint.
    ///
    /// `break_threshold` – impulse magnitude above which the joint snaps.
    /// `natural_length` – rest distance between body COMs.
    /// `stiffness` – spring stiffness (N/m).
    /// `damping` – spring damping (N·s/m).
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        natural_length: Real,
        stiffness: Real,
        damping: Real,
        break_threshold: Real,
    ) -> Self {
        Self {
            broken: false,
            break_threshold,
            accumulated_force: 0.0,
            body_a,
            body_b,
            stiffness,
            natural_length,
            damping,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Whether the joint has permanently broken.
    pub fn is_broken(&self) -> bool {
        self.broken
    }
    /// Reset the joint so it is active again (e.g. after a repair in game logic).
    pub fn reset(&mut self) {
        self.broken = false;
        self.accumulated_force = 0.0;
    }
    /// Current accumulated impulse (force proxy).
    pub fn current_force(&self) -> Real {
        self.accumulated_force
    }
}
/// Records the impulse applied by a constraint during the last velocity solve
/// step and exposes it as an approximation of constraint force / torque.
///
/// Wrap any joint with `JointForceMeter::new(inner)` to obtain per-step
/// force readout without modifying the underlying constraint.
#[derive(Debug, Clone)]
pub struct JointForceMeter<C: Constraint + Clone> {
    /// The wrapped constraint.
    pub inner: C,
    /// Accumulated impulse magnitude from the last `solve_velocity` call.
    pub last_impulse_magnitude: Real,
    /// Accumulated torque impulse magnitude from the last `solve_velocity` call.
    pub last_torque_magnitude: Real,
}
impl<C: Constraint + Clone> JointForceMeter<C> {
    /// Wrap an existing constraint.
    pub fn new(inner: C) -> Self {
        Self {
            inner,
            last_impulse_magnitude: 0.0,
            last_torque_magnitude: 0.0,
        }
    }
    /// Approximate constraint force from the last step given `dt`.
    pub fn force(&self, dt: Real) -> Real {
        if dt > 1e-12 {
            self.last_impulse_magnitude / dt
        } else {
            0.0
        }
    }
    /// Approximate constraint torque from the last step given `dt`.
    pub fn torque(&self, dt: Real) -> Real {
        if dt > 1e-12 {
            self.last_torque_magnitude / dt
        } else {
            0.0
        }
    }
}
/// A spherical joint with twist measurement.
///
/// Acts as a ball joint (constrains anchor coincidence) and additionally
/// exposes `compute_twist_angle` to measure the twist rotation about the
/// primary axis without enforcing any angular limit.
#[derive(Debug, Clone)]
pub struct SphericalJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Primary axis in body A's local frame (the "twist" axis).
    pub local_axis: Vec3,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl SphericalJoint {
    /// Create a new spherical joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis: Vec3,
    ) -> Self {
        let axis = if local_axis.norm() > 1e-10 {
            local_axis.normalize()
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis: axis,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Compute the twist angle of body B relative to body A about the primary axis.
    ///
    /// Returns the angle in radians in `[-π, π]`.  The twist is extracted from
    /// the relative quaternion by projecting out the swing component.
    ///
    /// Algorithm (swing-twist decomposition):
    ///   1. Compute `q_rel = q_a^{-1} * q_b`.
    ///   2. Extract the twist component whose axis is aligned with `local_axis`.
    ///   3. Return twice the arc-tangent of the projected vector part.
    pub fn compute_twist_angle(&self, bodies: &RigidBodySet) -> f64 {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let q_rel = a.rotation.inverse() * b.rotation;
        let axis = self.local_axis;
        let qv = Vec3::new(q_rel.i, q_rel.j, q_rel.k);
        let dot = qv.dot(&axis);
        let twist_vec = axis * dot;
        let twist_w = q_rel.w;
        let twist_sin_half = twist_vec.norm();
        let angle = 2.0 * twist_sin_half.atan2(twist_w.abs());
        let sign = if qv.dot(&axis) < 0.0 { -1.0 } else { 1.0 };
        angle * sign
    }
}
/// Spring-damper joint between two points on two bodies.
///
/// Exerts a spring force `F = -k * (distance - rest_length) - c * v_relative`
/// along the line connecting the two anchor points.
#[derive(Debug, Clone)]
pub struct SpringJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Rest length of the spring.
    pub rest_length: Real,
    /// Spring stiffness (N/m).
    pub stiffness: Real,
    /// Damping coefficient (N*s/m).
    pub damping: Real,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
    pub(super) direction: Vec3,
    pub(super) current_distance: Real,
}
impl SpringJoint {
    /// Create a new spring joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        rest_length: Real,
        stiffness: Real,
        damping: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            rest_length,
            stiffness,
            damping,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            direction: Vec3::zeros(),
            current_distance: 0.0,
        }
    }
}
/// Cylindrical joint: allows both rotation around and translation along a single axis.
///
/// Two degrees of freedom are left unconstrained: rotation around and translation
/// along the cylinder axis.  Four DOFs are constrained (two linear perpendicular
/// to the axis, two rotational perpendicular to the axis).
#[derive(Debug, Clone)]
pub struct CylindricalJoint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Local anchor on body A.
    pub local_anchor_a: Vec3,
    /// Local anchor on body B.
    pub local_anchor_b: Vec3,
    /// Cylinder axis in body A's local frame.
    pub local_axis: Vec3,
    /// Optional lower translation limit along the axis.
    pub lower_limit: Option<Real>,
    /// Optional upper translation limit along the axis.
    pub upper_limit: Option<Real>,
    pub(super) world_axis: Vec3,
    pub(super) perp1: Vec3,
    pub(super) perp2: Vec3,
    pub(super) r_a: Vec3,
    pub(super) r_b: Vec3,
}
impl CylindricalJoint {
    /// Create a new cylindrical joint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        local_anchor_a: Vec3,
        local_anchor_b: Vec3,
        local_axis: Vec3,
    ) -> Self {
        Self {
            body_a,
            body_b,
            local_anchor_a,
            local_anchor_b,
            local_axis: local_axis.normalize(),
            lower_limit: None,
            upper_limit: None,
            world_axis: Vec3::zeros(),
            perp1: Vec3::zeros(),
            perp2: Vec3::zeros(),
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
        }
    }
    /// Set translation limits along the axis.
    pub fn with_limits(mut self, lower: Real, upper: Real) -> Self {
        self.lower_limit = Some(lower);
        self.upper_limit = Some(upper);
        self
    }
    /// Current displacement of body B relative to body A along the cylinder axis.
    pub fn current_displacement(&self, bodies: &RigidBodySet) -> Real {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return 0.0,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return 0.0,
        };
        let delta = b.position - a.position;
        delta.dot(&self.world_axis)
    }
}
