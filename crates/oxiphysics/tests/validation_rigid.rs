// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase 21.1 validation — rigid-body solver against analytical closed-form
//! reference problems.
//!
//! Three tests, each anchored to a closed-form physics result:
//!
//! * `test_restitution_ladder` — a ball is dropped from `h = 1 m` under
//!   gravity onto a static body with restitution `e = 0.5`. Rebound apex
//!   heights follow a geometric series `h_n = h_0 * e^(2n)`, so the first
//!   rebound reaches `0.25 m` (the second impact's drop height) and the
//!   next one `0.0625 m`. The baseline `rigid_bounce_height_2` locks the
//!   `0.25 m` value with an 8% relative tolerance.
//!
//! * `test_angular_momentum_conservation` — a torque-free box in zero
//!   gravity with an initial spin is integrated for 10 s. The magnitude
//!   of the world-frame angular momentum `‖R · I_local · ω‖` must be
//!   conserved within 1% (classical Euler result).
//!
//! * `test_stacked_box_equilibrium` — three unit boxes stacked along `+y`
//!   under Earth gravity are integrated for 3 s; the stack must reach a
//!   near-rest state with small penetration, small lateral drift, and
//!   small residual kinetic energy.
//!
//! Tests 1 and 3 run the full `PhysicsWorld` pipeline: contacts come from
//! the attached collider shapes and the collider restitution / friction
//! reach the sequential-impulse solver.

#[path = "regression_harness.rs"]
mod harness;

use oxiphysics::geometry::{BoxShape, Shape, Sphere};
use oxiphysics::rigid::world::PhysicsWorld;
use oxiphysics::rigid::{BodyType, Collider, RigidBody, RigidBodySet, integrate_bodies};
use oxiphysics_core::math::{Mat3, Vec3};
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Test 1 — restitution ladder (geometric bounce series)
// ---------------------------------------------------------------------------

/// Drop a ball from `h = 1 m` onto a static body with `e = 0.5` and record
/// the rebound apex heights. The analytical values are
/// `h_1 = 1.0 * 0.5^2 = 0.25 m` and `h_2 = 0.25 * 0.5^2 = 0.0625 m`.
///
/// Both colliders carry `e = 0.5` (the contact uses `max(e_a, e_b)`), the
/// ball is a radius-0.5 sphere and the ground a wide box whose top face is
/// the plane `y = -0.5`, so the ball rests with its centre at `y = 0`.
#[test]
fn test_restitution_ladder() {
    // dt = 1/240 s gives sub-millimetre bounce resolution without
    // inflating the test runtime.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 240.0);

    // Ball — mass 1, sphere collider r = 0.5. Damping disabled so kinetic
    // energy is only lost through the coefficient of restitution.
    let mut ball = RigidBody::new(1.0);
    ball.transform.position = Vec3::new(0.0, 1.0, 0.0);
    ball.linear_damping = 0.0;
    ball.angular_damping = 0.0;
    let ball_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
    ball.set_inertia_from_shape(ball_shape.as_ref());
    ball.update_world_inertia();
    let ball_h = world.add_rigid_body(ball);
    world.add_collider(
        Collider::new(ball_shape).with_restitution(0.5),
        Some(ball_h),
    );

    // Static ground: top face at y = -0.5, i.e. the ball touches it when its
    // centre reaches y = 0 (1 m below the release height).
    let mut ground = RigidBody::new_static();
    ground.transform.position = Vec3::new(0.0, -1.0, 0.0);
    let ground_h = world.add_rigid_body(ground);
    world.add_collider(
        Collider::new(Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)))).with_restitution(0.5),
        Some(ground_h),
    );

    // Step the world for 4 s (960 steps): long enough for several bounces.
    // An apex is the step at which the vertical velocity turns from upward to
    // downward; apexes lower than 1 cm are contact jitter, not bounces.
    let mut prev_vy = 0.0_f64;
    let mut peaks: Vec<f64> = Vec::new();
    for _ in 0..960 {
        world.step();
        let b = world
            .get_body(ball_h)
            .expect("ball handle inserted above must exist");
        let y = b.transform.position.y;
        let vy = b.velocity.y;
        if prev_vy > 0.0 && vy <= 0.0 && y > 0.01 {
            peaks.push(y);
        }
        prev_vy = vy;
    }

    assert!(
        peaks.len() >= 2,
        "expected at least two rebound apexes, got {}: {peaks:?}",
        peaks.len()
    );
    let h1 = peaks[0];
    let h2 = peaks[1];

    let baseline = harness::load_baseline("rigid_bounce_height_2")
        .expect("baseline `rigid_bounce_height_2` must exist");
    eprintln!(
        "[test_restitution_ladder] rebound apexes = {peaks:?} (analytical 0.25, 0.0625 m; baseline tol_rel = {:.1}%)",
        baseline.tolerance_rel * 100.0
    );
    assert_close!(h1, baseline);
    let expected_h2 = h1 * 0.25;
    assert!(
        (h2 - expected_h2).abs() < 0.15 * expected_h2,
        "second rebound {h2} must be e^2 = 0.25 of the first ({expected_h2})"
    );
}

// ---------------------------------------------------------------------------
// Test 2 — angular momentum conservation (torque-free spin)
// ---------------------------------------------------------------------------

/// A box with side lengths `0.2 × 0.4 × 0.3 m` and uniform unit mass is
/// given an initial angular velocity `ω = (1.0, 0.1, 0.0) rad/s` under
/// zero gravity and zero damping. Stepping only the integrator
/// (`integrate_bodies`) for 1000 × 0.01 s = 10 s must conserve the
/// world-frame angular momentum magnitude `‖R · I_local · ω‖` within 1%.
///
/// Inertia tensor for a solid box of full extents `(sx, sy, sz)` and mass
/// `m` has diagonal components `m/12 · (sy² + sz²)`, `m/12 · (sx² + sz²)`,
/// `m/12 · (sx² + sy²)`. For `(0.2, 0.4, 0.3)` at `m = 1`:
///
/// * `Ixx = (0.16 + 0.09) / 12 = 0.020833…`
/// * `Iyy = (0.04 + 0.09) / 12 = 0.010833…`
/// * `Izz = (0.04 + 0.16) / 12 = 0.016667…`
#[test]
fn test_angular_momentum_conservation() {
    let mut set = RigidBodySet::new();

    let mut body = RigidBody::new(1.0);
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body.local_inertia = Mat3::from_diagonal(&Vec3::new(
        0.020_833_333_333_333_332,
        0.010_833_333_333_333_334,
        0.016_666_666_666_666_666,
    ));
    body.update_world_inertia();
    body.angular_velocity = Vec3::new(1.0, 0.1, 0.0);
    let handle = set.insert(body);

    // Helper: world-frame angular momentum L = R · I_local · ω.
    let world_angular_momentum = |set: &RigidBodySet| -> Vec3 {
        let b = set.get(handle).expect("handle inserted above must exist");
        let rot = b.transform.rotation.to_rotation_matrix();
        let l_local = b.local_inertia * b.angular_velocity;
        rot.matrix() * l_local
    };

    let initial_l = world_angular_momentum(&set);
    let initial_mag = initial_l.norm();
    assert!(
        initial_mag > 0.0,
        "initial angular momentum must be non-zero"
    );

    let gravity = Vec3::zeros();
    let dt = 0.01_f64;
    let steps = 1000usize;

    let mut max_drift: f64 = 0.0;
    for _ in 0..steps {
        integrate_bodies(&mut set, dt, &gravity);
        let mag = world_angular_momentum(&set).norm();
        let drift = (mag - initial_mag).abs() / initial_mag;
        if drift > max_drift {
            max_drift = drift;
        }
    }

    let final_mag = world_angular_momentum(&set).norm();
    let final_drift = (final_mag - initial_mag).abs() / initial_mag;
    eprintln!(
        "[test_angular_momentum_conservation] steps = {steps}, dt = {dt} s, \
         |L_initial| = {initial_mag:.9}, |L_final| = {final_mag:.9}, \
         final_drift = {final_drift:.3e}, max_drift = {max_drift:.3e}",
    );

    assert!(
        max_drift < 0.01,
        "angular momentum drift exceeds 1% tolerance: \
         max_drift = {max_drift:.3e} ({:.3}%), final_drift = {final_drift:.3e} ({:.3}%), \
         |L_initial| = {initial_mag:.9}, |L_final| = {final_mag:.9}",
        max_drift * 100.0,
        final_drift * 100.0,
    );
}

// ---------------------------------------------------------------------------
// Test 3 — stacked box equilibrium
// ---------------------------------------------------------------------------

/// Three unit boxes stacked along `+y` under Earth gravity must settle to
/// near rest within 3 s. The test asserts:
///
/// * total downward penetration below the expected contact plane `< 2e-3 m`
///   (i.e. the stack has not sunk into the ground);
/// * lateral drift from the spawn column `< 5e-3 m` per body;
/// * per-body translational kinetic energy `< 1e-3 J` (i.e. the stack has
///   come to rest).
///
/// Each body carries a unit `BoxShape` collider (half-extent 0.5) and the
/// ground is a wide static box whose top face is the plane `y = 0`, so the
/// box-box manifolds are full 4-point face contacts.
#[test]
fn test_stacked_box_equilibrium() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 100.0);

    // Three unit-mass unit boxes resting at y = 0.5, 1.5, 2.5.
    let mut handles = Vec::with_capacity(3);
    let expected_rest_y = [0.5_f64, 1.5, 2.5];
    for &y in &expected_rest_y {
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(0.0, y + 0.5, 0.0); // drop from +0.5 m above rest
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        let shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
        b.set_inertia_from_shape(shape.as_ref());
        b.update_world_inertia();
        let h = world.add_rigid_body(b);
        world.add_collider(Collider::new(shape), Some(h));
        handles.push(h);
    }

    // Static ground box: top face at y = 0, aligned with expected_rest_y[0] = 0.5.
    let mut ground = RigidBody::new_static();
    ground.transform.position = Vec3::new(0.0, -0.5, 0.0);
    let ground_h = world.add_rigid_body(ground);
    world.add_collider(
        Collider::new(Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)))),
        Some(ground_h),
    );

    // Step the world for 3 s (300 steps).
    for _ in 0..300 {
        world.step();
    }

    // Sanity check every body against the analytical resting configuration.
    for (i, &h) in handles.iter().enumerate() {
        let b = world
            .get_body(h)
            .expect("stacked body handle inserted above must exist");
        assert_eq!(
            b.body_type,
            BodyType::Dynamic,
            "stacked body {i} was demoted to non-dynamic type"
        );
        let pos = b.transform.position;
        let vel = b.velocity;
        let ke = 0.5 * b.mass * vel.norm_squared();

        let penetration = (expected_rest_y[i] - pos.y).max(0.0);
        let lateral = (pos.x * pos.x + pos.z * pos.z).sqrt();

        eprintln!(
            "[test_stacked_box_equilibrium] body {i}: pos = ({:.6}, {:.6}, {:.6}), \
             vel = ({:.6}, {:.6}, {:.6}), penetration = {penetration:.3e}, \
             lateral = {lateral:.3e}, KE = {ke:.3e}",
            pos.x, pos.y, pos.z, vel.x, vel.y, vel.z,
        );

        assert!(
            penetration < 2e-3,
            "body {i} penetration {penetration:.3e} m exceeds 2e-3 m tolerance \
             (expected y ~ {:.3}, got y = {:.6})",
            expected_rest_y[i],
            pos.y,
        );
        assert!(
            lateral < 5e-3,
            "body {i} lateral drift {lateral:.3e} m exceeds 5e-3 m tolerance \
             (pos = ({:.6}, {:.6}, {:.6}))",
            pos.x,
            pos.y,
            pos.z,
        );
        assert!(
            ke < 1e-3,
            "body {i} kinetic energy {ke:.3e} J exceeds 1e-3 J tolerance \
             (vel = ({:.6}, {:.6}, {:.6}))",
            vel.x,
            vel.y,
            vel.z,
        );
    }
}
