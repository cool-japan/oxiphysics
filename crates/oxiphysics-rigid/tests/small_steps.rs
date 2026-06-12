// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Integration tests for the PhysX-5 / Macklin et al. 2019 "Small Steps"
//! sub-stepping solver (`PhysicsWorld::step_small_steps`).

use oxiphysics_rigid::world::PhysicsWorld;
use oxiphysics_rigid::{BodyType, RigidBody};

/// Spawn an undamped dynamic unit-mass sphere (radius 0.5) at `(x, y, z)` and
/// return its handle.
fn spawn_sphere(world: &mut PhysicsWorld, x: f64, y: f64, z: f64) -> oxiphysics_core::BodyHandle {
    let mut body = RigidBody::new(1.0);
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    let handle = world.add_rigid_body(body);
    let b = world.get_body_mut(handle).expect("body just added");
    b.transform.position.x = x;
    b.transform.position.y = y;
    b.transform.position.z = z;
    handle
}

/// Spawn a static sphere whose radius is `mass.cbrt() * 0.5`, centred at the
/// origin. The `.mass` field is set AFTER `new_static()` (which sets mass = 0)
/// because the broad/narrowphase radius formula reads `body.mass.cbrt()`.
fn spawn_static_sphere(world: &mut PhysicsWorld, mass: f64) -> oxiphysics_core::BodyHandle {
    let mut ground = RigidBody::new_static();
    ground.mass = mass;
    world.add_rigid_body(ground)
}

#[test]
fn test_small_steps_100_sphere_pile_stable() {
    // scale: REDUCED from a 100-sphere 10x10 grid to a single vertical column
    // of 5 spheres directly above the ground centre.
    //
    // Why: the ground is a *sphere* of radius 10. A flat 10x10 grid places its
    // outer spheres (x,z = +/-4.5) far out on the strongly curved upper surface
    // where the contact normal points steeply sideways, so those spheres roll
    // off the dome and fall past the ground centre (observed y < 0). That is an
    // inherent geometric instability of stacking a flat grid on a convex sphere,
    // NOT a solver defect (the on-axis settling, energy and bounce tests all
    // pass). The spec authorises this reduction ("IF the 100-grid won't stay
    // bounded ... REDUCE to a single vertical column of 5 spheres ... and
    // COMMENT why"). A vertical column sits on the dome apex where the normal is
    // vertical, so it stays bounded and exercises the same stacking solver.
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.solver_config.restitution = 0.1;

    // Static ground sphere: mass 8000 => radius 8000.cbrt()*0.5 = 20*0.5 = 10.0.
    spawn_static_sphere(&mut world, 8000.0);

    // Vertical column of 5 dynamic unit spheres (radius 0.5) on the dome apex.
    let mut handles = Vec::new();
    for &y in &[10.5_f64, 11.5, 12.5, 13.5, 14.5] {
        handles.push(spawn_sphere(&mut world, 0.0, y, 0.0));
    }
    assert_eq!(handles.len(), 5);

    for _ in 0..120 {
        world.step_small_steps();
    }

    // Stability: every dynamic body is finite and has not tunnelled below the
    // ground sphere centre.
    for h in &handles {
        let b = world.get_body(*h).expect("dynamic body present");
        assert!(b.body_type == BodyType::Dynamic);
        let p = &b.transform.position;
        assert!(p.x.is_finite(), "x not finite: {}", p.x);
        assert!(p.y.is_finite(), "y not finite: {}", p.y);
        assert!(p.z.is_finite(), "z not finite: {}", p.z);
        assert!(p.y > 0.0, "body tunnelled below ground centre: y = {}", p.y);
    }

    // Bounded penetration over the contacts cached at the final step.
    let max_depth = world
        .last_contacts
        .iter()
        .map(|c| c.depth)
        .fold(0.0_f64, f64::max);
    assert!(
        max_depth < 0.05,
        "max penetration too large: {}",
        max_depth
    );
}

#[test]
fn test_small_steps_energy_non_increasing() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.solver_config.restitution = 0.0;

    spawn_static_sphere(&mut world, 8000.0);

    // Vertical stack of 5 unit spheres directly above the ground centre.
    for &y in &[10.5_f64, 11.5, 12.5, 13.5, 14.5] {
        spawn_sphere(&mut world, 0.0, y, 0.0);
    }

    let mut ke_history = Vec::with_capacity(100);
    for _ in 0..100 {
        world.step_small_steps();
        ke_history.push(world.kinetic_energy());
    }

    // Gravity adds KE during free fall, so we only bound the peak and require
    // the system to have settled by the end of the run.
    let max_ke = ke_history.iter().cloned().fold(0.0_f64, f64::max);
    assert!(max_ke.is_finite(), "KE diverged (NaN/Inf): {}", max_ke);
    assert!(max_ke < 1000.0, "peak KE too high: {}", max_ke);

    let final_ke = *ke_history.last().expect("ran at least one frame");
    assert!(final_ke < 1.0, "final KE not small: {}", final_ke);

    // Strict non-increase only holds after landing: check the last 20 frames.
    // (Gravity injects KE during the free-fall phase earlier in the run.)
    let n = ke_history.len();
    for i in (n - 20)..n {
        let prev = ke_history[i - 1];
        let cur = ke_history[i];
        assert!(
            cur <= prev + 1e-3,
            "KE increased in settled tail at frame {}: {} -> {}",
            i,
            prev,
            cur
        );
    }
}

#[test]
fn test_small_steps_fast_body_no_energy_gain() {
    let mut world = PhysicsWorld::new([0.0, 0.0, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.solver_config.restitution = 0.0;

    // Static wall sphere: mass 64 => radius 64.cbrt()*0.5 = 4*0.5 = 2.0 at origin.
    spawn_static_sphere(&mut world, 64.0);

    // One fast dynamic unit sphere at x = -3, moving toward the wall at 50 m/s.
    let ball = spawn_sphere(&mut world, -3.0, 0.0, 0.0);
    world
        .get_body_mut(ball)
        .expect("ball present")
        .velocity
        .x = 50.0;

    // Initial KE = 0.5 * 1 * 50^2 = 1250 J.
    for _ in 0..60 {
        world.step_small_steps();
    }

    let final_ke = world.kinetic_energy();
    assert!(
        final_ke <= 1250.0 + 1e-3,
        "small-steps solver injected energy: {} > 1250",
        final_ke
    );

    let p = &world.get_body(ball).expect("ball present").transform.position;
    assert!(p.x.is_finite());
    assert!(p.y.is_finite());
    assert!(p.z.is_finite());
}

#[test]
fn test_small_steps_bouncing_ball_restitution() {
    let mut world = PhysicsWorld::new([0.0, -9.81, 0.0], 1.0 / 60.0);
    world.solver_config.substeps = 4;
    world.solver_config.restitution = 0.5;

    spawn_static_sphere(&mut world, 8000.0);

    // Resting contact is at y = 10.5 (ground top 10.0 + ball radius 0.5).
    // Drop from y0 = 12.5 => H = 2.0 m above rest, starting from rest.
    let ball = spawn_sphere(&mut world, 0.0, 12.5, 0.0);

    // Measurement: read y and vy each frame. Flag `engaged` once the ball
    // descends below 10.7 (contact established). After engagement, the first
    // time vy transitions from negative to positive marks the start of the
    // rebound; we then track the peak y until vy turns negative again. The peak
    // rebound height above rest is h_prime = peak_y - 10.5.
    let mut engaged = false;
    let mut prev_vy = 0.0_f64;
    let mut tracking = false;
    let mut peak_y = f64::NEG_INFINITY;
    let mut h_prime = f64::NAN;

    for _ in 0..300 {
        world.step_small_steps();
        let b = world.get_body(ball).expect("ball present");
        let y = b.transform.position.y;
        let vy = b.velocity.y;

        if y < 10.7 {
            engaged = true;
        }

        if engaged && !tracking && prev_vy < 0.0 && vy >= 0.0 {
            // Rebound started.
            tracking = true;
            peak_y = y;
        } else if tracking {
            if vy >= 0.0 {
                if y > peak_y {
                    peak_y = y;
                }
            } else if h_prime.is_nan() {
                // vy turned negative again => peak reached; record once.
                h_prime = peak_y - 10.5;
            }
        }

        prev_vy = vy;
    }

    // If the ball was still ascending at the end, use the tracked peak.
    if h_prime.is_nan() && tracking {
        h_prime = peak_y - 10.5;
    }

    assert!(h_prime.is_finite(), "no rebound detected");
    assert!(h_prime > 0.0, "ball did not rebound: h' = {}", h_prime);

    // Analytic energy restitution: h'/H = e^2 = 0.25. Allow a tolerant window
    // (numerical sub-stepping, gravity during contact, discrete sampling).
    let ratio = h_prime / 2.0;
    assert!(
        (0.15..0.40).contains(&ratio),
        "rebound ratio out of window: h'/H = {} (h' = {})",
        ratio,
        h_prime
    );
}
