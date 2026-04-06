//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Mat3, Real, Vec3};

use super::types::{BodyState, BodyType, KinematicTarget};

/// A rigid body with full dynamics state.
#[derive(Debug, Clone)]
pub struct RigidBody {
    /// The body's spatial transform.
    pub transform: Transform,
    /// Linear velocity.
    pub velocity: Vec3,
    /// Angular velocity.
    pub angular_velocity: Vec3,
    /// Mass in kg.
    pub mass: Real,
    /// Inverse mass (0 for static bodies).
    pub inverse_mass: Real,
    /// Local-space inertia tensor.
    pub local_inertia: Mat3,
    /// Inverse of world-space inertia tensor (recomputed each frame).
    pub world_inverse_inertia: Mat3,
    /// Accumulated force for the current step.
    pub force_accumulator: Vec3,
    /// Accumulated torque for the current step.
    pub torque_accumulator: Vec3,
    /// Body type.
    pub body_type: BodyType,
    /// Activity state.
    pub state: BodyState,
    /// Linear damping factor (0..1).
    pub linear_damping: Real,
    /// Angular damping factor (0..1).
    pub angular_damping: Real,
    /// Time spent below sleep thresholds.
    pub sleep_timer: Real,
    /// Gravity scale (1.0 = normal, 0.0 = no gravity).
    pub gravity_scale: Real,
    /// Optional kinematic target for scripted motion.
    pub kinematic_target: Option<KinematicTarget>,
    /// Number of integration steps this body has undergone.
    pub step_count: u64,
}
