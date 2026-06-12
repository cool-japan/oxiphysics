// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body dynamics simulation.
//!
//! Provides rigid body state, integration, force/torque accumulators,
//! body and collider sets, and semi-implicit Euler integration.
//!
//! # Top-level utilities
//!
//! In addition to the submodule re-exports, this crate exposes:
//!
//! * [`integrate_bodies`] — semi-implicit Euler integration over a
//!   [`RigidBodySet`].
//! * [`update_sleeping`] — velocity-threshold sleeping detection.
//! * [`apply_force_and_wake`] — apply a force to a body and wake it if needed.
//! * [`BodyActivationEvent`] — enum signalling wake/sleep transitions.
//!
//! ## Examples
//!
//! ### Gravity integration
//!
//! ```rust
//! use oxiphysics_rigid::{RigidBodySet, RigidBody, integrate_bodies};
//! use oxiphysics_core::math::Vec3;
//!
//! let mut set = RigidBodySet::new();
//! let h = set.insert({
//!     let mut b = RigidBody::new(1.0);
//!     b.linear_damping = 0.0;
//!     b
//! });
//! let gravity = Vec3::new(0.0, -9.81, 0.0);
//! // After 10 steps of dt=0.1 s, the body should have fallen
//! for _ in 0..10 {
//!     integrate_bodies(&mut set, 0.1, &gravity);
//! }
//! let body = set.get(h).expect("body must exist");
//! assert!(body.transform.position.y < 0.0, "body should have fallen");
//! ```
//!
//! ### Apply force wakes a sleeping body
//!
//! ```rust
//! use oxiphysics_rigid::{RigidBodySet, RigidBody, BodyState, apply_force_and_wake};
//! use oxiphysics_core::math::Vec3;
//!
//! let mut set = RigidBodySet::new();
//! let h = set.insert({
//!     let mut b = RigidBody::new(1.0);
//!     b.state = BodyState::Sleeping;
//!     b
//! });
//! apply_force_and_wake(&mut set, h, Vec3::new(10.0, 0.0, 0.0));
//! let body = set.get(h).expect("body must exist");
//! assert_eq!(body.state, BodyState::Active, "body should be awake after force");
//! ```
//!
//! ### Stationary body falls asleep
//!
//! ```rust
//! use oxiphysics_rigid::{RigidBodySet, RigidBody, BodyActivationEvent, update_sleeping};
//! use oxiphysics_core::math::Vec3;
//!
//! let mut set = RigidBodySet::new();
//! // Body at rest (zero velocity), will fall asleep after ~0.5s
//! set.insert(RigidBody::new(1.0));
//! let dt = 0.1_f64;
//! let mut fell_asleep = false;
//! for _ in 0..20 {
//!     let events = update_sleeping(&mut set, dt, 0.05, 0.05);
//!     if events.iter().any(|e| matches!(e, BodyActivationEvent::FellAsleep(_))) {
//!         fell_asleep = true;
//!     }
//! }
//! assert!(fell_asleep, "stationary body should fall asleep");
//! ```
#![warn(missing_docs)]

mod error;
pub use error::*;

mod body;
pub use body::*;

mod collider;
pub use collider::*;

mod sets;
pub use sets::*;

pub mod articulated;
pub mod sleeping;
pub use sleeping::{
    ActivityLevel, AdaptiveSleepThreshold, EnergyDecayTracker, EnergyTracker, FrozenBody,
    IslandSleepDecision, IslandSleepEvaluator, IslandSleepManager, PseudoSleepBody, SleepConfig,
    SleepDurationTracker, SleepEventType, SleepHistory, SleepHysteresis, SleepManager,
    SleepPredictor, SleepState, SleepStats, SleepThresholdAdapter, SleepWakePattern, SleepingBody,
    SleepingSystem, VelocityDamper, WakePropagator, kinetic_energy_rigid, potential_energy_rigid,
};
pub mod ccd;
pub mod impulse;
pub mod motors;
pub mod pipeline;
pub mod solver;
pub mod world;

pub mod ragdoll;
pub use ragdoll::{
    ArticulatedPose, BalanceController, BallSocketJoint, BodySegmentInertia, BoneDescriptor,
    ConeTwistLimit, GroundContactDetector, HingeJoint, HumanoidRagdoll, HumanoidSkeleton,
    JointLimit, PdJointController, PosePdController, Ragdoll, RagdollBone, RagdollPose,
    RagdollState, SkeletonJoint, forward_kinematics, ragdoll_inertia_tensor,
    standard_human_segment_inertias,
};

pub mod fluid_coupling;
pub use fluid_coupling::{
    AddedMass, AddedMassTensor, BuoyancyForce, FloatingBody, FlowRegime, FluidCoupling,
    FluidProperties, FroudeKrylovForce, FroudeNumber, HullResistance, HydrodynamicDrag,
    MorisonElement, OrientedBuoyancy, PropellerThrust, ReynoldsUtils, SloshingModel,
    SphereBuoyancy, SubmergedBodyCoupling, SubmergedVolumeFraction, VivParams,
    VortexInducedVibration, WaveForce, added_mass_cylinder, added_mass_sphere,
    drag_coefficient_cylinder, drag_coefficient_sphere, slamming_force,
};

// ── BodyActivationEvent ───────────────────────────────────────────────────────

/// Event emitted when a rigid body transitions between active and sleeping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BodyActivationEvent {
    /// The body was previously sleeping and has been woken up.
    Woke(oxiphysics_core::BodyHandle),
    /// The body was previously active and has fallen asleep.
    FellAsleep(oxiphysics_core::BodyHandle),
}

// ── integrate_bodies ──────────────────────────────────────────────────────────

/// Integrate all active dynamic bodies in `set` by one time step `dt` under
/// a constant `gravity` vector.
///
/// Uses semi-implicit Euler: velocity is updated first, then position.
/// Static and sleeping bodies are skipped.
pub fn integrate_bodies(
    set: &mut RigidBodySet,
    dt: oxiphysics_core::math::Real,
    gravity: &oxiphysics_core::math::Vec3,
) {
    for (_, body) in set.iter_mut() {
        body.integrate_forces(dt, gravity);
        body.integrate_velocity(dt);
    }
}

// ── update_sleeping ───────────────────────────────────────────────────────────

/// Check every dynamic body in `set` and transition it to/from sleeping based
/// on velocity thresholds.
///
/// A body enters the `Sleeping` state after its linear and angular speeds have
/// both been below the respective thresholds for `time_before_sleep` seconds.
/// The `dt` parameter is the duration of the last simulation step.
///
/// Returns a list of [`BodyActivationEvent`]s describing all state changes.
pub fn update_sleeping(
    set: &mut RigidBodySet,
    dt: oxiphysics_core::math::Real,
    linear_threshold: oxiphysics_core::math::Real,
    angular_threshold: oxiphysics_core::math::Real,
) -> Vec<BodyActivationEvent> {
    let time_before_sleep = 0.5; // seconds — could be made a parameter
    let mut events = Vec::new();

    for (handle, body) in set.iter_mut() {
        let was_sleeping = body.state == BodyState::Sleeping;
        body.check_sleep(dt, linear_threshold, angular_threshold, time_before_sleep);
        let is_sleeping = body.state == BodyState::Sleeping;

        if !was_sleeping && is_sleeping {
            events.push(BodyActivationEvent::FellAsleep(handle));
        } else if was_sleeping && !is_sleeping {
            events.push(BodyActivationEvent::Woke(handle));
        }
    }

    events
}

// ── apply_force_and_wake ──────────────────────────────────────────────────────

/// Apply `force` to the body identified by `handle` and wake it if it was
/// sleeping.
///
/// No-ops silently if `handle` is invalid or the body is not dynamic.
pub fn apply_force_and_wake(
    set: &mut RigidBodySet,
    handle: oxiphysics_core::BodyHandle,
    force: oxiphysics_core::math::Vec3,
) {
    if let Some(body) = set.get_mut(handle) {
        if body.state == BodyState::Sleeping {
            body.wake_up();
        }
        body.apply_force(force);
    }
}

// ── apply_torque_and_wake ─────────────────────────────────────────────────────

/// Apply `torque` to the body identified by `handle` and wake it if it was
/// sleeping.
///
/// No-ops silently if `handle` is invalid or the body is not dynamic.
pub fn apply_torque_and_wake(
    set: &mut RigidBodySet,
    handle: oxiphysics_core::BodyHandle,
    torque: oxiphysics_core::math::Vec3,
) {
    if let Some(body) = set.get_mut(handle) {
        if body.state == BodyState::Sleeping {
            body.wake_up();
        }
        body.apply_torque(torque);
    }
}

// ── apply_impulse_and_wake ───────────────────────────────────────────────────

/// Apply an instantaneous `impulse` to the body and wake it.
///
/// Unlike forces, impulses immediately change velocity.
pub fn apply_impulse_and_wake(
    set: &mut RigidBodySet,
    handle: oxiphysics_core::BodyHandle,
    impulse: oxiphysics_core::math::Vec3,
) {
    if let Some(body) = set.get_mut(handle) {
        if body.state == BodyState::Sleeping {
            body.wake_up();
        }
        body.apply_impulse(impulse);
    }
}

// ── compute_kinetic_energy ──────────────────────────────────────────────────

/// Compute the total kinetic energy (linear) of a body: 0.5 * m * v^2.
pub fn kinetic_energy_linear(body: &RigidBody) -> oxiphysics_core::math::Real {
    0.5 * body.mass * body.velocity.norm_squared()
}

/// Compute the rotational kinetic energy: 0.5 * omega^T * I * omega.
pub fn kinetic_energy_rotational(body: &RigidBody) -> oxiphysics_core::math::Real {
    let i_omega = body.local_inertia * body.angular_velocity;
    0.5 * body.angular_velocity.dot(&i_omega)
}

/// Compute total kinetic energy (linear + rotational).
pub fn kinetic_energy_total(body: &RigidBody) -> oxiphysics_core::math::Real {
    kinetic_energy_linear(body) + kinetic_energy_rotational(body)
}

// ── angular momentum ────────────────────────────────────────────────────────

/// Compute the angular momentum of a body: L = I * omega.
pub fn angular_momentum(body: &RigidBody) -> oxiphysics_core::math::Vec3 {
    body.local_inertia * body.angular_velocity
}

/// Compute the linear momentum of a body: p = m * v.
pub fn linear_momentum(body: &RigidBody) -> oxiphysics_core::math::Vec3 {
    body.velocity * body.mass
}

// ── mass properties from shapes ─────────────────────────────────────────────

/// Compute mass properties (mass and inertia tensor) for a solid sphere.
///
/// * `radius` - Sphere radius.
/// * `density` - Material density in kg/m^3.
///
/// Returns `(mass, inertia_diagonal)` where the inertia tensor is diagonal
/// with identical entries I = (2/5) * m * r^2.
pub fn sphere_mass_properties(
    radius: oxiphysics_core::math::Real,
    density: oxiphysics_core::math::Real,
) -> (oxiphysics_core::math::Real, oxiphysics_core::math::Real) {
    let volume = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
    let mass = density * volume;
    let inertia = (2.0 / 5.0) * mass * radius * radius;
    (mass, inertia)
}

/// Compute mass properties for a solid box.
///
/// * `half_extents` - Half-extents `[hx, hy, hz]` of the box.
/// * `density` - Material density in kg/m^3.
///
/// Returns `(mass, [Ixx, Iyy, Izz])`.
pub fn box_mass_properties(
    half_extents: [oxiphysics_core::math::Real; 3],
    density: oxiphysics_core::math::Real,
) -> (
    oxiphysics_core::math::Real,
    [oxiphysics_core::math::Real; 3],
) {
    let [hx, hy, hz] = half_extents;
    let sx = 2.0 * hx;
    let sy = 2.0 * hy;
    let sz = 2.0 * hz;
    let volume = sx * sy * sz;
    let mass = density * volume;
    let ixx = mass / 12.0 * (sy * sy + sz * sz);
    let iyy = mass / 12.0 * (sx * sx + sz * sz);
    let izz = mass / 12.0 * (sx * sx + sy * sy);
    (mass, [ixx, iyy, izz])
}

/// Compute mass properties for a solid cylinder aligned along the Y axis.
///
/// * `radius` - Cylinder radius.
/// * `height` - Cylinder height.
/// * `density` - Material density in kg/m^3.
///
/// Returns `(mass, [Ixx, Iyy, Izz])`.
pub fn cylinder_mass_properties(
    radius: oxiphysics_core::math::Real,
    height: oxiphysics_core::math::Real,
    density: oxiphysics_core::math::Real,
) -> (
    oxiphysics_core::math::Real,
    [oxiphysics_core::math::Real; 3],
) {
    let volume = std::f64::consts::PI * radius * radius * height;
    let mass = density * volume;
    let iyy = 0.5 * mass * radius * radius;
    let ixx = mass / 12.0 * (3.0 * radius * radius + height * height);
    let izz = ixx;
    (mass, [ixx, iyy, izz])
}

// ── wake_all ────────────────────────────────────────────────────────────────

/// Wake all sleeping bodies in the set.
pub fn wake_all(set: &mut RigidBodySet) {
    for (_, body) in set.iter_mut() {
        if body.state == BodyState::Sleeping {
            body.wake_up();
        }
    }
}

/// Count the number of sleeping bodies in the set.
pub fn count_sleeping(set: &RigidBodySet) -> usize {
    set.iter()
        .filter(|(_, b)| b.state == BodyState::Sleeping)
        .count()
}

/// Count the number of active dynamic bodies in the set.
pub fn count_active(set: &RigidBodySet) -> usize {
    set.iter()
        .filter(|(_, b)| b.is_dynamic() && b.state == BodyState::Active)
        .count()
}

/// Compute the center of mass of all dynamic bodies.
pub fn center_of_mass(set: &RigidBodySet) -> oxiphysics_core::math::Vec3 {
    let mut total_mass = 0.0;
    let mut weighted_pos = oxiphysics_core::math::Vec3::zeros();
    for (_, body) in set.iter() {
        if body.is_dynamic() {
            total_mass += body.mass;
            weighted_pos += body.transform.position * body.mass;
        }
    }
    if total_mass > 1e-10 {
        weighted_pos / total_mass
    } else {
        oxiphysics_core::math::Vec3::zeros()
    }
}

/// Compute the total linear momentum of all dynamic bodies.
pub fn total_linear_momentum(set: &RigidBodySet) -> oxiphysics_core::math::Vec3 {
    let mut total = oxiphysics_core::math::Vec3::zeros();
    for (_, body) in set.iter() {
        if body.is_dynamic() {
            total += body.velocity * body.mass;
        }
    }
    total
}

/// Compute total kinetic energy of all dynamic bodies.
pub fn total_kinetic_energy(set: &RigidBodySet) -> oxiphysics_core::math::Real {
    let mut total = 0.0;
    for (_, body) in set.iter() {
        if body.is_dynamic() {
            total += kinetic_energy_total(body);
        }
    }
    total
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod lib_tests {
    use super::*;
    use oxiphysics_core::math::Vec3;

    fn gravity() -> Vec3 {
        Vec3::new(0.0, -9.81, 0.0)
    }

    // ── insert / remove / handle validity ────────────────────────────────

    #[test]
    fn insert_and_get_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        assert!(set.get(h).is_some());
    }

    #[test]
    fn remove_invalidates_handle() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.remove(h);
        assert!(set.get(h).is_none());
    }

    #[test]
    fn handle_generation_prevents_stale_access() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.remove(h1);
        let h2 = set.insert(RigidBody::new(2.0));
        // Slot is reused but generation is bumped
        assert!(set.get(h1).is_none(), "stale handle should be invalid");
        assert!(set.get(h2).is_some(), "fresh handle should be valid");
    }

    #[test]
    fn iter_yields_all_active_bodies() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(2.0));
        set.insert(RigidBody::new(3.0));
        assert_eq!(set.iter().count(), 3);
    }

    // ── gravity integration ───────────────────────────────────────────────

    #[test]
    fn gravity_integration_increases_downward_velocity() {
        let mut set = RigidBodySet::new();
        let h = set.insert({
            let mut b = RigidBody::new(1.0);
            b.linear_damping = 0.0;
            b
        });
        integrate_bodies(&mut set, 1.0, &gravity());
        let vy = set.get(h).unwrap().velocity.y;
        assert!(vy < 0.0, "vy={vy} should be negative after gravity");
        assert!((vy + 9.81).abs() < 1e-6, "vy={vy} expected ≈ -9.81");
    }

    #[test]
    fn static_body_unaffected_by_gravity_integration() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new_static());
        let pos_before = set.get(h).unwrap().transform.position;
        integrate_bodies(&mut set, 1.0, &gravity());
        let pos_after = set.get(h).unwrap().transform.position;
        assert!(
            (pos_after - pos_before).norm() < 1e-12,
            "static body must not move"
        );
    }

    // ── sleeping threshold ────────────────────────────────────────────────

    #[test]
    fn sleeping_threshold_triggers_fell_asleep_event() {
        let mut set = RigidBodySet::new();
        let h = set.insert({
            let mut b = RigidBody::new(1.0);
            b.velocity = Vec3::zeros();
            b.angular_velocity = Vec3::zeros();
            b
        });
        // Drive the sleep timer above 0.5 s (time_before_sleep)
        let dt = 0.016;
        let steps = ((0.5_f64 / dt).ceil() as usize) + 5;
        let mut all_events: Vec<BodyActivationEvent> = Vec::new();
        for _ in 0..steps {
            let ev = update_sleeping(&mut set, dt, 0.1, 0.1);
            all_events.extend(ev);
        }
        let fell_asleep = all_events
            .iter()
            .any(|e| matches!(e, BodyActivationEvent::FellAsleep(_)));
        assert!(fell_asleep, "expected a FellAsleep event");
        assert_eq!(set.get(h).unwrap().state, BodyState::Sleeping);
    }

    #[test]
    fn fast_body_does_not_fall_asleep() {
        let mut set = RigidBodySet::new();
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.velocity = Vec3::new(10.0, 0.0, 0.0);
            b
        });
        let mut events: Vec<BodyActivationEvent> = Vec::new();
        for _ in 0..100 {
            let ev = update_sleeping(&mut set, 0.016, 0.1, 0.1);
            events.extend(ev);
        }
        assert!(
            events
                .iter()
                .all(|e| !matches!(e, BodyActivationEvent::FellAsleep(_))),
            "fast body should not fall asleep"
        );
    }

    // ── wake on force ─────────────────────────────────────────────────────

    #[test]
    fn apply_force_and_wake_wakes_sleeping_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b
        });
        apply_force_and_wake(&mut set, h, Vec3::new(1.0, 0.0, 0.0));
        let body = set.get(h).unwrap();
        assert_eq!(body.state, BodyState::Active, "body should be awake");
    }

    #[test]
    fn apply_force_and_wake_accumulates_force_on_active_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        apply_force_and_wake(&mut set, h, Vec3::new(5.0, 0.0, 0.0));
        let fx = set.get(h).unwrap().force_accumulator.x;
        assert!(
            (fx - 5.0).abs() < 1e-12,
            "force_accumulator.x={fx} expected 5.0"
        );
    }

    #[test]
    fn apply_force_and_wake_invalid_handle_no_panic() {
        let mut set = RigidBodySet::new();
        let bad_handle = oxiphysics_core::BodyHandle::new(99, 0);
        // Must not panic
        apply_force_and_wake(&mut set, bad_handle, Vec3::new(1.0, 0.0, 0.0));
    }

    // ── BodyActivationEvent ───────────────────────────────────────────────

    #[test]
    fn woke_event_emitted_when_sleeping_body_given_velocity() {
        let mut set = RigidBodySet::new();
        let h = set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b.sleep_timer = 1.0;
            b
        });
        // Manually give it a high velocity, then call update_sleeping
        if let Some(body) = set.get_mut(h) {
            body.velocity = Vec3::new(5.0, 0.0, 0.0);
            body.state = BodyState::Active; // simulate external wake
        }
        // Now trigger with zero velocity — should not emit Woke again
        if let Some(body) = set.get_mut(h) {
            body.velocity = Vec3::zeros();
        }
        let events = update_sleeping(&mut set, 0.016, 0.1, 0.1);
        // Body was Active and velocity is now zero → no Woke event expected
        assert!(
            events
                .iter()
                .all(|e| !matches!(e, BodyActivationEvent::Woke(_))),
            "no Woke event expected here"
        );
    }

    // ── torque and impulse wake tests ────────────────────────────────────

    #[test]
    fn apply_torque_and_wake_wakes_sleeping_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b
        });
        apply_torque_and_wake(&mut set, h, Vec3::new(0.0, 1.0, 0.0));
        assert_eq!(set.get(h).unwrap().state, BodyState::Active);
    }

    #[test]
    fn apply_impulse_and_wake_changes_velocity() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(2.0));
        apply_impulse_and_wake(&mut set, h, Vec3::new(4.0, 0.0, 0.0));
        // impulse / mass = 4 / 2 = 2
        let vx = set.get(h).unwrap().velocity.x;
        assert!((vx - 2.0).abs() < 1e-10, "vx={vx}");
    }

    // ── kinetic energy tests ─────────────────────────────────────────────

    #[test]
    fn kinetic_energy_linear_stationary_is_zero() {
        let body = RigidBody::new(5.0);
        assert!(kinetic_energy_linear(&body).abs() < 1e-10);
    }

    #[test]
    fn kinetic_energy_linear_moving_body() {
        let mut body = RigidBody::new(2.0);
        body.velocity = Vec3::new(3.0, 0.0, 0.0);
        // KE = 0.5 * 2 * 9 = 9
        assert!((kinetic_energy_linear(&body) - 9.0).abs() < 1e-10);
    }

    #[test]
    fn kinetic_energy_total_includes_rotation() {
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::new(1.0, 0.0, 0.0);
        body.angular_velocity = Vec3::new(0.0, 1.0, 0.0);
        let total = kinetic_energy_total(&body);
        let linear = kinetic_energy_linear(&body);
        let rotational = kinetic_energy_rotational(&body);
        assert!((total - linear - rotational).abs() < 1e-10);
        assert!(total > linear, "total KE should include rotational");
    }

    // ── momentum tests ───────────────────────────────────────────────────

    #[test]
    fn linear_momentum_correct() {
        let mut body = RigidBody::new(3.0);
        body.velocity = Vec3::new(2.0, 0.0, 0.0);
        let p = linear_momentum(&body);
        assert!((p.x - 6.0).abs() < 1e-10);
    }

    #[test]
    fn angular_momentum_zero_when_stationary() {
        let body = RigidBody::new(1.0);
        let l = angular_momentum(&body);
        assert!(l.norm() < 1e-10);
    }

    // ── mass properties tests ────────────────────────────────────────────

    #[test]
    fn sphere_mass_properties_unit_sphere() {
        let (mass, inertia) = sphere_mass_properties(1.0, 1000.0);
        let expected_volume = (4.0 / 3.0) * std::f64::consts::PI;
        let expected_mass = 1000.0 * expected_volume;
        assert!((mass - expected_mass).abs() < 1e-6, "mass={mass}");
        let expected_inertia = (2.0 / 5.0) * expected_mass;
        assert!(
            (inertia - expected_inertia).abs() < 1e-6,
            "inertia={inertia}"
        );
    }

    #[test]
    fn box_mass_properties_unit_cube() {
        let (mass, inertias) = box_mass_properties([0.5, 0.5, 0.5], 1000.0);
        // Volume = 1.0 m^3
        assert!((mass - 1000.0).abs() < 1e-6);
        // Ixx = Iyy = Izz = m/12 * (1^2 + 1^2) = 1000/6
        for i in inertias {
            assert!((i - 1000.0 / 6.0).abs() < 1e-6, "inertia={i}");
        }
    }

    #[test]
    fn cylinder_mass_properties_basic() {
        let (mass, inertias) = cylinder_mass_properties(1.0, 2.0, 1000.0);
        let expected_volume = std::f64::consts::PI * 1.0 * 1.0 * 2.0;
        let expected_mass = 1000.0 * expected_volume;
        assert!((mass - expected_mass).abs() < 1e-6);
        // Iyy = 0.5 * m * r^2
        let expected_iyy = 0.5 * expected_mass;
        assert!((inertias[1] - expected_iyy).abs() < 1e-6);
        // Ixx = Izz
        assert!((inertias[0] - inertias[2]).abs() < 1e-10);
    }

    // ── set utility tests ────────────────────────────────────────────────

    #[test]
    fn wake_all_wakes_sleeping_bodies() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b
        });
        let h2 = set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b
        });
        wake_all(&mut set);
        assert_eq!(set.get(h1).unwrap().state, BodyState::Active);
        assert_eq!(set.get(h2).unwrap().state, BodyState::Active);
    }

    #[test]
    fn count_sleeping_and_active() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.state = BodyState::Sleeping;
            b
        });
        set.insert(RigidBody::new(1.0));
        assert_eq!(count_sleeping(&set), 1);
        assert_eq!(count_active(&set), 2);
    }

    #[test]
    fn center_of_mass_two_equal_bodies() {
        let mut set = RigidBodySet::new();
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.transform.position = Vec3::new(0.0, 0.0, 0.0);
            b
        });
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.transform.position = Vec3::new(10.0, 0.0, 0.0);
            b
        });
        let com = center_of_mass(&set);
        assert!((com.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn total_linear_momentum_conservation() {
        let mut set = RigidBodySet::new();
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.velocity = Vec3::new(3.0, 0.0, 0.0);
            b
        });
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.velocity = Vec3::new(-3.0, 0.0, 0.0);
            b
        });
        let p = total_linear_momentum(&set);
        assert!(p.norm() < 1e-10, "momentum should be zero: {:?}", p);
    }

    #[test]
    fn total_kinetic_energy_positive() {
        let mut set = RigidBodySet::new();
        set.insert({
            let mut b = RigidBody::new(1.0);
            b.velocity = Vec3::new(5.0, 0.0, 0.0);
            b
        });
        let ke = total_kinetic_energy(&set);
        assert!(ke > 0.0, "KE should be positive: {ke}");
        assert!((ke - 12.5).abs() < 1e-10, "KE=0.5*1*25=12.5, got {ke}");
    }
}

pub mod aerospace_maneuver;
pub mod aerospace_rigid;
pub mod buoyancy_model;
pub mod cables;
pub mod cfd_coupling;
pub mod constraint_forces;
pub mod constraint_graph;
pub mod contact_dynamics;
pub mod contact_geometry;
pub mod continuum_robot;
pub mod crowd_simulation;
pub mod deformable_coupling;
pub mod deformable_rigid;
pub mod destruction_rigid;
pub mod fluid_dynamics;
pub mod fluid_structure;
pub mod fracture;
pub mod fracture_dynamics;
pub mod fracture_rigid;
pub mod granular;
pub mod granular_media;
pub mod granular_rigid;
pub mod gyroscope_rigid;
pub mod gyroscopic;
pub mod impact_dynamics;
pub mod impact_mechanics;
pub mod kinematics;
pub mod locomotion;
pub mod marine_dynamics;
pub mod marine_rigid;
pub mod mechanism_rigid;
pub mod multibody;
pub mod multibody_dynamics;
pub mod multiscale_rigid;
pub mod neural_control;
pub mod orbital_mechanics;
pub mod planetary_mechanics;
pub mod robotic_systems;
pub mod satellite_dynamics;
pub mod sensor_model;
pub mod sensor_simulation;
pub mod soft_constraints;
pub mod space_mechanics;
pub mod spatial_algebra;
pub mod swarm_dynamics;
pub mod swarm_physics;
pub mod swarm_robotics;
pub mod terrain_dynamics;
pub mod terrain_interaction;
pub mod thruster_dynamics;
pub mod train_dynamics;
pub mod train_rigid;
pub mod underwater_rigid;
pub mod underwater_vehicle;
pub mod wind_turbine;
