// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`PhysicsWorld::step_with_accumulator`]: input
//! validation, the substep cap and the leftover-time bookkeeping.

use super::{Body, PhysicsWorld};
use crate::error::PhysicsError;
use crate::types::BodyHandle;

const FIXED_DT: f64 = 1.0 / 60.0;

fn world_with_falling_body() -> (PhysicsWorld, BodyHandle) {
    let mut world = PhysicsWorld::new();
    let mut body = Body::new_dynamic([0.0, 10.0, 0.0], 1.0);
    body.linear_damping = 0.0;
    let handle = world.add_body(body);
    (world, handle)
}

/// Name of the offending field of an `InvalidInput` error.
fn invalid_field(error: &PhysicsError) -> &str {
    match error {
        PhysicsError::InvalidInput { field, .. } => field.as_str(),
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

#[test]
fn zero_fixed_dt_is_rejected_instead_of_hanging() {
    let (mut world, handle) = world_with_falling_body();
    world.fixed_dt = 0.0;
    let error = world
        .step_with_accumulator(0.5)
        .expect_err("fixed_dt = 0 must be an error");
    assert_eq!(invalid_field(&error), "fixed_dt");
    assert_eq!(world.time, 0.0, "nothing may have been simulated");
    assert_eq!(world.accumulator, 0.0, "the accumulator must be untouched");
    let y = world.body(handle).map(|b| b.transform.position.y);
    assert_eq!(y, Some(10.0));
}

#[test]
fn negative_nan_and_infinite_fixed_dt_are_rejected() {
    for bad in [-FIXED_DT, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let (mut world, _) = world_with_falling_body();
        world.fixed_dt = bad;
        let error = world
            .step_with_accumulator(0.5)
            .expect_err("a non-positive or non-finite fixed_dt must be an error");
        assert_eq!(invalid_field(&error), "fixed_dt", "fixed_dt = {bad}");
        assert_eq!(world.accumulator, 0.0);
    }
}

#[test]
fn zero_max_substeps_is_rejected() {
    let (mut world, _) = world_with_falling_body();
    world.max_substeps = 0;
    let error = world
        .step_with_accumulator(0.5)
        .expect_err("max_substeps = 0 would freeze the world");
    assert_eq!(invalid_field(&error), "max_substeps");
    assert_eq!(world.accumulator, 0.0);
}

#[test]
fn nan_dt_is_rejected_and_leaves_the_accumulator_alone() {
    let (mut world, _) = world_with_falling_body();
    world.accumulator = 0.004;
    let error = world
        .step_with_accumulator(f64::NAN)
        .expect_err("NaN dt must be an error");
    assert_eq!(invalid_field(&error), "dt");
    assert_eq!(
        world.accumulator, 0.004,
        "the accumulator must be unchanged"
    );
    // The world is not poisoned: a later valid call still steps.
    let steps = world
        .step_with_accumulator(FIXED_DT)
        .expect("a finite dt after the rejected one must work");
    assert_eq!(steps, 1);
}

#[test]
fn infinite_and_negative_dt_are_rejected() {
    for bad in [f64::INFINITY, f64::NEG_INFINITY, -1e-3, -1.0] {
        let (mut world, _) = world_with_falling_body();
        world.accumulator = 0.001;
        let error = world
            .step_with_accumulator(bad)
            .expect_err("a negative or non-finite dt must be an error");
        assert_eq!(invalid_field(&error), "dt", "dt = {bad}");
        assert_eq!(world.accumulator, 0.001);
        assert_eq!(world.time, 0.0);
    }
}

#[test]
fn non_finite_accumulator_is_reported_not_propagated() {
    let (mut world, _) = world_with_falling_body();
    world.accumulator = f64::NAN;
    let error = world
        .step_with_accumulator(FIXED_DT)
        .expect_err("a poisoned accumulator must be an error");
    assert_eq!(invalid_field(&error), "accumulator");
}

#[test]
fn a_huge_frame_runs_at_most_max_substeps_and_keeps_one_step_of_backlog() {
    let (mut world, _) = world_with_falling_body();
    assert_eq!(world.max_substeps, 8, "documented default");
    let steps = world
        .step_with_accumulator(1.0e6)
        .expect("a huge but finite dt is valid input");
    assert_eq!(steps, world.max_substeps);
    assert!(
        (world.time - 8.0 * FIXED_DT).abs() < 1e-12,
        "time = {}",
        world.time
    );
    assert!(
        (world.accumulator - FIXED_DT).abs() < 1e-15,
        "the backlog must be clipped to one fixed step, got {}",
        world.accumulator
    );
    // The single retained step is paid out by the next call, then it is empty.
    let next = world
        .step_with_accumulator(0.0)
        .expect("zero dt is valid input");
    assert_eq!(next, 1);
    assert!(world.accumulator < FIXED_DT);
}

#[test]
fn max_substeps_is_configurable() {
    let (mut world, _) = world_with_falling_body();
    world.max_substeps = 3;
    let steps = world.step_with_accumulator(1.0).expect("valid input");
    assert_eq!(steps, 3);
    assert!((world.time - 3.0 * FIXED_DT).abs() < 1e-12);
}

#[test]
fn ordinary_frames_consume_whole_substeps_and_carry_the_remainder() {
    let (mut world, handle) = world_with_falling_body();
    let steps = world
        .step_with_accumulator(2.5 * FIXED_DT)
        .expect("valid input");
    assert_eq!(steps, 2);
    assert!((world.accumulator - 0.5 * FIXED_DT).abs() < 1e-12);
    assert!((world.time - 2.0 * FIXED_DT).abs() < 1e-12);

    // 0.5 + 0.6 = 1.1 fixed steps: exactly one more substep, 0.1 left over.
    let steps = world
        .step_with_accumulator(0.6 * FIXED_DT)
        .expect("valid input");
    assert_eq!(steps, 1);
    assert!((world.accumulator - 0.1 * FIXED_DT).abs() < 1e-12);

    // A frame shorter than one substep integrates nothing at all.
    let before = world.time;
    let steps = world
        .step_with_accumulator(0.1 * FIXED_DT)
        .expect("valid input");
    assert_eq!(steps, 0);
    assert_eq!(world.time, before);

    // Three substeps of semi-implicit Euler under g = 9.81 from rest.
    let y = world.body(handle).map(|b| b.transform.position.y);
    let expected = 10.0 - 9.81 * FIXED_DT * FIXED_DT * (1.0 + 2.0 + 3.0);
    assert!(
        y.is_some_and(|y| (y - expected).abs() < 1e-9),
        "y = {y:?}, expected {expected}"
    );
}
