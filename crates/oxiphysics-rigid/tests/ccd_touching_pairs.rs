// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for how the CCD pass treats touching and overlapping pairs.
//!
//! `check_sphere_sphere` / `compute_collision_time_sphere` used to return
//! `Some(0.0)` for every overlapping pair, even a separating one, and
//! `run_with_stats` advanced such a pair by a forced 1e-8 per pass until its
//! iteration cap ran out, so the rest of the step was never integrated.

use oxiphysics_rigid::ccd::{CcdPipeline, CcdSweep, SubstepIntegrator, SweptBody, TunnelTest};

fn sphere(x: f64, vx: f64, radius: f64, inv_mass: f64) -> SweptBody {
    SweptBody::new([x, 0.0, 0.0], [vx, 0.0, 0.0], radius, inv_mass)
}

// ── time of impact of overlapping spheres ───────────────────────────────────

#[test]
fn overlapping_spheres_impact_only_while_they_approach() {
    let at = |va: f64, vb: f64| {
        let a = sphere(0.0, va, 1.0, 1.0);
        let b = sphere(0.5, vb, 1.0, 1.0);
        let pipeline = TunnelTest::check_sphere_sphere(&a, &b, 1.0);
        let sweep = CcdSweep::compute_collision_time_sphere(
            a.pos, a.vel, a.radius, b.pos, b.vel, b.radius, 1.0,
        );
        assert_eq!(pipeline, sweep, "va {va}, vb {vb}");
        pipeline
    };
    // Approaching: impact now.
    assert_eq!(at(1.0, 0.0), Some(0.0));
    assert_eq!(at(0.0, -1.0), Some(0.0));
    assert_eq!(at(0.3, -0.3), Some(0.0));
    // Resting, co-moving and separating overlaps are not impacts.
    assert_eq!(at(0.0, 0.0), None);
    assert_eq!(at(2.0, 2.0), None);
    assert_eq!(at(-1.0, 0.0), None);
    assert_eq!(at(0.0, 1.0), None);
    assert_eq!(at(-0.3, 0.3), None);
}

#[test]
fn separated_spheres_keep_their_first_contact_time() {
    // Unchanged behaviour for pairs that do not overlap at t = 0.
    let a = sphere(-5.0, 1.0, 0.5, 1.0);
    let b = sphere(5.0, -1.0, 0.5, 1.0);
    let t = TunnelTest::check_sphere_sphere(&a, &b, 10.0).expect("head-on spheres collide");
    assert!((t - 4.5).abs() < 1e-12, "t = {t}");
    assert_eq!(TunnelTest::check_sphere_sphere(&a, &b, 4.0), None);
    // Exactly touching and closing: contact now (to rounding).
    let touching_a = sphere(-0.5, 1.0, 0.5, 1.0);
    let touching_b = sphere(0.5, -1.0, 0.5, 1.0);
    let t = TunnelTest::check_sphere_sphere(&touching_a, &touching_b, 1.0)
        .expect("touching, closing spheres collide");
    assert!(t.abs() < 1e-12, "t = {t}");
}

// ── run_with_stats ──────────────────────────────────────────────────────────

#[test]
fn resting_overlapping_pair_advances_the_full_step() {
    // Both spheres overlap by 0.5 and move together at 1 m/s.
    let mut bodies = vec![sphere(0.0, 1.0, 0.5, 1.0), sphere(0.5, 1.0, 0.5, 1.0)];
    let stats = CcdPipeline::new(2.0).run_with_stats(&mut bodies, &[0, 1], 0.5);
    assert_eq!(stats.events_found, 0);
    assert_eq!(stats.impulses_applied, 0);
    assert!((stats.time_covered - 2.0).abs() < 1e-12);
    assert!(
        (bodies[0].pos[0] - 2.0).abs() < 1e-12,
        "{:?}",
        bodies[0].pos
    );
    assert!(
        (bodies[1].pos[0] - 2.5).abs() < 1e-12,
        "{:?}",
        bodies[1].pos
    );
    assert_eq!(bodies[0].vel, [1.0, 0.0, 0.0]);
    assert_eq!(bodies[1].vel, [1.0, 0.0, 0.0]);

    // A completely stationary overlap as well.
    let mut still = vec![sphere(0.0, 0.0, 0.5, 1.0), sphere(0.5, 0.0, 0.5, 1.0)];
    let stats = CcdPipeline::new(1.0).run_with_stats(&mut still, &[0, 1], 0.5);
    assert_eq!(stats.events_found, 0);
    assert!((stats.time_covered - 1.0).abs() < 1e-12);
    assert_eq!(still[0].pos, [0.0; 3]);
}

#[test]
fn separating_overlapping_pair_advances_the_full_step() {
    let mut bodies = vec![sphere(0.0, -1.0, 0.5, 1.0), sphere(0.5, 2.0, 0.5, 1.0)];
    let stats = CcdPipeline::new(3.0).run_with_stats(&mut bodies, &[0, 1], 1.0);
    assert_eq!(stats.events_found, 0);
    assert!((stats.time_covered - 3.0).abs() < 1e-12);
    assert!((bodies[0].pos[0] + 3.0).abs() < 1e-12);
    assert!((bodies[1].pos[0] - 6.5).abs() < 1e-12);
}

#[test]
fn approaching_overlapping_pair_is_resolved_once_and_the_step_completes() {
    // Overlapping by 0.2 and closing at 2 m/s: an impact at t = 0.
    let mut bodies = vec![sphere(-0.4, 1.0, 0.5, 1.0), sphere(0.4, -1.0, 0.5, 1.0)];
    let stats = CcdPipeline::new(4.0).run_with_stats(&mut bodies, &[0, 1], 1.0);
    assert_eq!(stats.events_found, 1);
    assert_eq!(stats.impulses_applied, 1);
    assert!((stats.time_covered - 4.0).abs() < 1e-12);
    // Elastic exchange of velocities, then 4 s of separation.
    assert!(
        (bodies[0].vel[0] + 1.0).abs() < 1e-12,
        "{:?}",
        bodies[0].vel
    );
    assert!(
        (bodies[1].vel[0] - 1.0).abs() < 1e-12,
        "{:?}",
        bodies[1].vel
    );
    assert!((bodies[0].pos[0] - (-0.4 + 1e-8 - (4.0 - 1e-8))).abs() < 1e-9);
    assert!((bodies[1].pos[0] - (0.4 - 1e-8 + (4.0 - 1e-8))).abs() < 1e-9);
}

#[test]
fn forced_minimum_advance_never_exceeds_the_remaining_time() {
    // A step shorter than the forced 1e-8 advance of an impact at t = 0 must
    // still end exactly at its own end.
    let dt = 5e-9;
    let mut bodies = vec![sphere(-0.4, 1.0, 0.5, 1.0), sphere(0.4, -1.0, 0.5, 1.0)];
    let stats = CcdPipeline::new(dt).run_with_stats(&mut bodies, &[0, 1], 1.0);
    assert_eq!(stats.impulses_applied, 1);
    assert!(
        (stats.time_covered - dt).abs() < 1e-18,
        "{}",
        stats.time_covered
    );
    assert!(
        (bodies[0].pos[0] - (-0.4 + dt)).abs() < 1e-15,
        "{:?}",
        bodies[0].pos
    );
    assert!(
        (bodies[1].pos[0] - (0.4 - dt)).abs() < 1e-15,
        "{:?}",
        bodies[1].pos
    );
}

#[test]
fn pair_of_immovable_bodies_does_not_stall_the_step() {
    // Two kinematic bodies on a collision course cannot be resolved; they must
    // not consume the iteration budget and the whole step is integrated.
    let mut bodies = vec![sphere(-2.0, 1.0, 0.5, 0.0), sphere(2.0, -1.0, 0.5, 0.0)];
    let stats = CcdPipeline::new(10.0).run_with_stats(&mut bodies, &[0, 1], 1.0);
    assert_eq!(stats.events_found, 0);
    assert_eq!(stats.impulses_applied, 0);
    assert!((stats.time_covered - 10.0).abs() < 1e-12);
    assert!((bodies[0].pos[0] - 8.0).abs() < 1e-12);
    assert!((bodies[1].pos[0] + 8.0).abs() < 1e-12);
}

#[test]
fn immovable_pair_next_to_a_dynamic_collision_is_ignored() {
    // Slice indices 0/1: kinematic pair on a collision course (unresolvable);
    // indices 2/3: a dynamic head-on pair that must still be resolved.
    let mut bodies = vec![
        SweptBody::new([0.0, 20.0, 0.0], [1.0, 0.0, 0.0], 0.5, 0.0),
        SweptBody::new([4.0, 20.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 0.0),
        sphere(-4.5, 1.0, 0.5, 1.0),
        sphere(4.5, -1.0, 0.5, 1.0),
    ];
    let stats = CcdPipeline::new(10.0).run_with_stats(&mut bodies, &[0, 1, 2, 3], 1.0);
    assert_eq!(stats.impulses_applied, 1);
    assert!((stats.time_covered - 10.0).abs() < 1e-9);
    assert!((bodies[2].vel[0] + 1.0).abs() < 1e-12);
    assert!((bodies[3].vel[0] - 1.0).abs() < 1e-12);
    assert!((bodies[0].pos[0] - 10.0).abs() < 1e-12);
}

#[test]
fn step_is_fully_consumed_when_the_impact_cap_is_reached() {
    // A small sphere shuttles between two static spheres: every 0.2 s it hits
    // one of them, 100 impacts in 20 s, more than the 64 a pass resolves.
    // The remainder of the step must still be integrated.
    let wall_a = SweptBody::new([-1.2, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
    let wall_b = SweptBody::new([1.2, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
    let ball = sphere(0.0, 1.0, 0.1, 1.0);
    let mut bodies = vec![ball, wall_a, wall_b];
    let stats = CcdPipeline::new(20.0).run_with_stats(&mut bodies, &[0, 1, 2], 1.0);
    assert_eq!(stats.events_found, 64);
    assert_eq!(stats.impulses_applied, 64);
    assert!(
        (stats.time_covered - 20.0).abs() < 1e-9,
        "time covered {}",
        stats.time_covered
    );
    // 64 impacts take 0.1 + 63 * 0.2 = 12.7 s; the last one is on wall A at
    // x = -0.1 and sends the ball to +x, which then flies for 7.3 s.
    assert!((bodies[0].vel[0] - 1.0).abs() < 1e-9, "{:?}", bodies[0].vel);
    assert!((bodies[0].pos[0] - 7.2).abs() < 1e-6, "{:?}", bodies[0].pos);
    // The walls never moved.
    assert_eq!(bodies[1].pos, [-1.2, 0.0, 0.0]);
    assert_eq!(bodies[2].pos, [1.2, 0.0, 0.0]);
}

#[test]
fn response_skips_events_that_cannot_be_answered() {
    // resolve_event keeps its contract for the helper used by the pass.
    let integrator = SubstepIntegrator::new(1);
    let mut bodies = vec![sphere(-0.5, 1.0, 0.5, 0.0), sphere(0.5, -1.0, 0.5, 0.0)];
    let events = CcdPipeline::new(1.0).detect_tunneling(&bodies, &[0, 1]);
    assert_eq!(events.len(), 1);
    assert!(!integrator.resolve_event(&mut bodies, &events[0], 1.0));
    assert_eq!(bodies[0].vel, [1.0, 0.0, 0.0]);
    assert_eq!(bodies[1].vel, [-1.0, 0.0, 0.0]);
}
