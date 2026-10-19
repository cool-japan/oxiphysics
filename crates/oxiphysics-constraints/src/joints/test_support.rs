// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Helpers shared by the joint regression tests.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_rigid::{RigidBody, RigidBodySet};

use crate::traits::Constraint;

/// Dynamic body with a diagonal (principal-axis) inertia tensor and no damping.
pub(super) fn body_with_inertia(mass: Real, principal_inertia: Vec3) -> RigidBody {
    let mut body = RigidBody::new(mass);
    body.local_inertia = Mat3::from_diagonal(&principal_inertia);
    body.update_world_inertia();
    body
}

/// Rotate a body about the world axis `axis` by `angle` and refresh its
/// world-space inertia.
pub(super) fn rotate_body(body: &mut RigidBody, axis: Vec3, angle: Real) {
    let unit_axis = oxiphysics_core::math::Unit::new_normalize(axis);
    body.transform.rotation = Quat::from_axis_angle(&unit_axis, angle) * body.transform.rotation;
    body.update_world_inertia();
}

/// Velocity of the anchor (`local_*` in the body frame) of each body, A minus B.
pub(super) fn anchor_relative_velocity(
    bodies: &RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    local_anchor_a: Vec3,
    local_anchor_b: Vec3,
) -> Vec3 {
    let a = bodies.get(handle_a).expect("body A exists");
    let b = bodies.get(handle_b).expect("body B exists");
    let r_a = a.transform.rotation * local_anchor_a;
    let r_b = b.transform.rotation * local_anchor_b;
    (a.velocity + a.angular_velocity.cross(&r_a)) - (b.velocity + b.angular_velocity.cross(&r_b))
}

/// Angle in `[0, π]` of the rotation of body B relative to body A.
pub(super) fn relative_rotation_angle(
    bodies: &RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
) -> Real {
    let a = bodies.get(handle_a).expect("body A exists");
    let b = bodies.get(handle_b).expect("body B exists");
    (a.transform.rotation.inverse() * b.transform.rotation).angle()
}

/// Angle in `[0, π]` of the world-space rotation of one body.
pub(super) fn rotation_angle(bodies: &RigidBodySet, handle: BodyHandle) -> Real {
    bodies
        .get(handle)
        .expect("body exists")
        .transform
        .rotation
        .angle()
}

/// Total linear momentum of all bodies.
pub(super) fn linear_momentum(bodies: &RigidBodySet) -> Vec3 {
    bodies
        .iter()
        .map(|(_, body)| body.velocity * body.mass)
        .fold(Vec3::zeros(), |sum, p| sum + p)
}

/// Total kinetic energy of all bodies.
pub(super) fn total_kinetic_energy(bodies: &RigidBodySet) -> Real {
    bodies.iter().map(|(_, body)| body.kinetic_energy()).sum()
}

/// One full solver step for a single constraint, in the usual order: integrate
/// forces, `prepare`, velocity iterations, position iterations, integrate
/// velocities.
pub(super) fn step_joint<C: Constraint>(
    joint: &mut C,
    bodies: &mut RigidBodySet,
    gravity: &Vec3,
    dt: Real,
    velocity_iterations: usize,
    position_iterations: usize,
) {
    for (_, body) in bodies.iter_mut() {
        body.integrate_forces(dt, gravity);
    }
    joint.prepare(bodies, dt);
    for _ in 0..velocity_iterations {
        joint.solve_velocity(bodies, dt);
    }
    for _ in 0..position_iterations {
        joint.solve_position(bodies, dt);
    }
    for (_, body) in bodies.iter_mut() {
        body.integrate_velocity(dt);
    }
}
