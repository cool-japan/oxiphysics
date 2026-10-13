// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Damping of a [`RigidBody`] is a rate in 1/s and is zero by default.
//!
//! Driven through [`PhysicsWorld::step`], a default body must therefore follow
//! the closed-form projectile motion of semi-implicit Euler whatever the step
//! size (the former per-step factor of 0.99 was dt-dependent and saturated a
//! free fall near `99 g dt`, about 16 m/s at 60 Hz). A body with an explicit
//! damping rate approaches the *same* terminal velocity `-g / rate` at every
//! step size.

use oxiphysics_rigid::RigidBody;
use oxiphysics_rigid::world::PhysicsWorld;

const GRAVITY: f64 = 9.81;
const STEP_RATES: [u32; 4] = [30, 60, 240, 1000];

#[test]
fn default_body_free_falls_with_the_same_trajectory_at_every_step_rate() {
    for steps_per_second in STEP_RATES {
        let dt = 1.0 / f64::from(steps_per_second);
        let mut world = PhysicsWorld::new([0.0, -GRAVITY, 0.0], dt);
        let handle = world.add_rigid_body(RigidBody::new(1.0));
        for _ in 0..steps_per_second {
            world.step();
        }
        let body = world.get_body(handle).expect("the body is in the world");
        // After exactly one second: v = -g t, and semi-implicit Euler (velocity
        // first, then position) gives y = -g dt^2 n (n + 1) / 2.
        let n = f64::from(steps_per_second);
        let expected_speed = -GRAVITY;
        let expected_height = -GRAVITY * dt * dt * n * (n + 1.0) / 2.0;
        assert!(
            ((body.velocity.y - expected_speed) / expected_speed).abs() < 1e-9,
            "{steps_per_second} Hz: vy = {}, expected {expected_speed}",
            body.velocity.y
        );
        assert!(
            ((body.transform.position.y - expected_height) / expected_height).abs() < 1e-9,
            "{steps_per_second} Hz: y = {}, expected {expected_height}",
            body.transform.position.y
        );
        assert!(body.velocity.x == 0.0 && body.velocity.z == 0.0);
    }
}

#[test]
fn default_projectile_keeps_its_horizontal_speed() {
    // With the old default every step removed 1 % of the horizontal velocity:
    // 0.99^240 = 9 % of it survived one second at 240 Hz.
    for steps_per_second in STEP_RATES {
        let dt = 1.0 / f64::from(steps_per_second);
        let mut world = PhysicsWorld::new([0.0, -GRAVITY, 0.0], dt);
        let mut projectile = RigidBody::new(2.0);
        projectile.velocity.x = 12.0;
        projectile.velocity.y = 20.0;
        let handle = world.add_rigid_body(projectile);
        for _ in 0..steps_per_second {
            world.step();
        }
        let body = world.get_body(handle).expect("the body is in the world");
        assert!(
            (body.velocity.x - 12.0).abs() < 1e-12,
            "{steps_per_second} Hz: vx = {}",
            body.velocity.x
        );
        assert!(
            (body.velocity.y - (20.0 - GRAVITY)).abs() < 1e-9,
            "{steps_per_second} Hz: vy = {}",
            body.velocity.y
        );
    }
}

#[test]
fn damped_body_reaches_the_step_independent_terminal_velocity() {
    // dv/dt = -g - rate v settles at -g / rate. The implicit-Euler damping
    // factor 1 / (1 + dt rate) has exactly that fixed point for every dt.
    let rate = 2.0;
    let terminal = -GRAVITY / rate;
    for steps_per_second in STEP_RATES {
        let dt = 1.0 / f64::from(steps_per_second);
        let mut world = PhysicsWorld::new([0.0, -GRAVITY, 0.0], dt);
        let mut body = RigidBody::new(1.0);
        body.linear_damping = rate;
        let handle = world.add_rigid_body(body);
        // 20 s is 40 time constants (1 / rate = 0.5 s).
        for _ in 0..(20 * steps_per_second) {
            world.step();
        }
        let vy = world
            .get_body(handle)
            .expect("the body is in the world")
            .velocity
            .y;
        assert!(
            ((vy - terminal) / terminal).abs() < 1e-9,
            "{steps_per_second} Hz: vy = {vy}, terminal velocity {terminal}"
        );
    }
}
