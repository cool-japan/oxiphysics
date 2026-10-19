// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the lateral dynamics of the single-track motorcycle
//! simulator.
//!
//! The simulator used to have no lateral-velocity state at all (the slip angles
//! only saw the yaw rate) and never integrated the steering angle, so the steering
//! rate had no effect on the steering angle.  The bike now carries the lateral
//! velocity `v_y` with `m (v_y' + v r) = F_y,f cos(delta) + F_y,r` and
//! `I_z r' = a F_y,f cos(delta) - b F_y,r`, and the steering angle follows the
//! steering rate.

use oxiphysics_vehicle::motorcycle_dynamics::MotorcycleDynamicsSimulator;

const GRAVITY: f64 = 9.81;

/// A bike rolling freely at `speed` (m/s) whose roll motion is frozen and which
/// has no aerodynamic drag, so that only the lateral dynamics are exercised and
/// the speed stays constant.
fn frozen_bike(speed: f64) -> MotorcycleDynamicsSimulator {
    let mut sim = MotorcycleDynamicsSimulator::new();
    sim.state.speed = speed;
    sim.state.rear_wheel_omega = speed / sim.params.rear_wheel_radius;
    sim.state.front_wheel_omega = speed / sim.params.front_wheel_radius;
    sim.params.roll_inertia = 1.0e12;
    sim.fairing.frontal_area = 0.0;
    sim
}

fn assert_lateral_state_is_sane(sim: &MotorcycleDynamicsSimulator) {
    let s = &sim.state;
    for (name, value) in [
        ("lateral velocity", s.lateral_velocity),
        ("yaw rate", s.yaw_rate),
        ("steer angle", s.steer_angle),
        ("speed", s.speed),
        ("lean", s.lean_angle),
    ] {
        assert!(value.is_finite(), "{name} is not finite: {value}");
    }
}

#[test]
fn a_new_bike_has_no_lateral_velocity() {
    let sim = MotorcycleDynamicsSimulator::new();
    assert_eq!(sim.state.lateral_velocity, 0.0);
}

#[test]
fn straight_running_builds_up_no_lateral_motion() {
    let mut sim = frozen_bike(25.0);
    for _ in 0..400 {
        sim.step(0.005);
    }
    assert_eq!(sim.state.lateral_velocity, 0.0);
    assert_eq!(sim.state.yaw_rate, 0.0);
    assert_eq!(sim.state.steer_angle, 0.0);
}

// ---------------------------------------------------------------------------
// Steering
// ---------------------------------------------------------------------------

#[test]
fn steer_angle_follows_the_steer_rate() {
    let mut sim = frozen_bike(20.0);
    sim.state.steer_angle = 0.01;
    sim.state.steer_rate = 0.2;
    for _ in 0..25 {
        sim.step(0.004);
    }
    // 25 steps of 4 ms at 0.2 rad/s add 0.02 rad
    assert!(
        (sim.state.steer_angle - 0.03).abs() < 1e-12,
        "steer angle {}",
        sim.state.steer_angle
    );
    // and back again: 50 steps of 4 ms at -0.1 rad/s take off 0.02 rad
    sim.state.steer_rate = -0.1;
    for _ in 0..50 {
        sim.step(0.004);
    }
    assert!(
        (sim.state.steer_angle - 0.01).abs() < 1e-12,
        "steer angle {}",
        sim.state.steer_angle
    );
}

#[test]
fn constant_steer_angle_stays_constant_without_a_steer_rate() {
    let mut sim = frozen_bike(20.0);
    sim.state.steer_angle = 0.015;
    for _ in 0..200 {
        sim.step(0.005);
    }
    assert_eq!(sim.state.steer_angle, 0.015);
}

#[test]
fn a_steering_rate_turns_the_bike_towards_the_steered_side() {
    for sign in [1.0_f64, -1.0] {
        let mut sim = frozen_bike(20.0);
        sim.state.steer_rate = sign * 0.05;
        for _ in 0..100 {
            sim.step(0.005);
            assert_lateral_state_is_sane(&sim);
        }
        assert!(
            sim.state.steer_angle * sign > 0.0,
            "steer angle {}",
            sim.state.steer_angle
        );
        assert!(
            sim.state.yaw_rate * sign > 0.0,
            "steering {sign} gave yaw rate {}",
            sim.state.yaw_rate
        );
    }
}

#[test]
fn right_and_left_steering_are_mirror_images() {
    let run = |steer: f64| {
        let mut sim = frozen_bike(18.0);
        sim.state.steer_angle = steer;
        for _ in 0..300 {
            sim.step(0.005);
        }
        (sim.state.yaw_rate, sim.state.lateral_velocity)
    };
    let (yaw_right, v_y_right) = run(0.01);
    let (yaw_left, v_y_left) = run(-0.01);
    assert!(yaw_right > 0.0);
    assert!(
        (yaw_right + yaw_left).abs() < 1e-12,
        "{yaw_right} vs {yaw_left}"
    );
    assert!(
        (v_y_right + v_y_left).abs() < 1e-12,
        "{v_y_right} vs {v_y_left}"
    );
}

// ---------------------------------------------------------------------------
// Lateral velocity
// ---------------------------------------------------------------------------

#[test]
fn steady_cornering_matches_the_linear_bicycle_model() {
    let (speed, steer) = (20.0, 0.005);
    let mut sim = frozen_bike(speed);
    sim.state.steer_angle = steer;
    for _ in 0..600 {
        sim.step(0.005);
        assert_lateral_state_is_sane(&sim);
    }
    // Without aerodynamics the axle loads are the static ones, so the cornering
    // stiffness of each axle is proportional to its load: neutral steer.
    let wheelbase = sim.params.wheelbase;
    let a = sim.params.a;
    let b = wheelbase - a;
    let expected_yaw = speed * steer / wheelbase;
    assert!(
        (sim.state.yaw_rate - expected_yaw).abs() < 0.01 * expected_yaw,
        "yaw rate {} rad/s, bicycle model {expected_yaw} rad/s",
        sim.state.yaw_rate
    );
    // Force balance of the rear axle gives v_y = r (b - v^2 / (kappa g)) with the
    // slope kappa = B C D of the Magic Formula per unit load.
    let kappa = sim.rear_tire.lat.b * sim.rear_tire.lat.c * sim.rear_tire.lat.d;
    let expected_v_y = expected_yaw * (b - speed * speed / (kappa * GRAVITY));
    assert!(
        expected_v_y < 0.0,
        "at 20 m/s the tail slides to the outside of the turn"
    );
    assert!(
        (sim.state.lateral_velocity - expected_v_y).abs() < 0.03 * expected_v_y.abs(),
        "lateral velocity {} m/s, bicycle model {expected_v_y} m/s",
        sim.state.lateral_velocity
    );
    // the speed was not touched by the cornering
    assert!(
        (sim.state.speed - speed).abs() < 1e-6,
        "speed {}",
        sim.state.speed
    );
}

#[test]
fn a_lateral_velocity_disturbance_decays() {
    // A bike sliding sideways is pulled back by the slip it creates at both axles.
    let mut sim = frozen_bike(20.0);
    sim.state.lateral_velocity = 0.5;
    let mut peak_after_half_a_second: f64 = 0.0;
    for step in 0..400 {
        sim.step(0.005);
        assert_lateral_state_is_sane(&sim);
        if step >= 100 {
            peak_after_half_a_second =
                peak_after_half_a_second.max(sim.state.lateral_velocity.abs());
        }
    }
    assert!(
        peak_after_half_a_second < 0.05,
        "lateral velocity still reaches {peak_after_half_a_second} m/s"
    );
    assert!(
        sim.state.lateral_velocity.abs() < 1e-3 && sim.state.yaw_rate.abs() < 1e-3,
        "v_y {} m/s, r {} rad/s after 2 s",
        sim.state.lateral_velocity,
        sim.state.yaw_rate
    );
}

#[test]
fn a_bike_at_rest_has_no_lateral_motion() {
    let mut sim = frozen_bike(0.0);
    sim.state.lateral_velocity = 0.3;
    sim.state.yaw_rate = 0.4;
    sim.state.steer_angle = 0.2;
    sim.step(0.01);
    assert_eq!(sim.state.lateral_velocity, 0.0);
    assert_eq!(sim.state.yaw_rate, 0.0);
    assert_eq!(sim.state.speed, 0.0);
}

#[test]
fn a_crawling_bike_pivots_about_its_rear_contact_point() {
    let speed = 0.05;
    let mut sim = frozen_bike(speed);
    sim.state.steer_angle = 0.1;
    sim.step(0.01);
    let wheelbase = sim.params.wheelbase;
    let b = wheelbase - sim.params.a;
    let expected_yaw = speed * 0.1 / wheelbase;
    assert!(
        (sim.state.yaw_rate - expected_yaw).abs() < 1e-12,
        "yaw rate {}",
        sim.state.yaw_rate
    );
    assert!((sim.state.lateral_velocity - b * expected_yaw).abs() < 1e-12);
}

#[test]
fn a_zero_time_step_leaves_the_lateral_state_untouched() {
    let mut sim = frozen_bike(15.0);
    sim.state.lateral_velocity = 0.2;
    sim.state.yaw_rate = 0.1;
    sim.state.steer_angle = 0.02;
    sim.state.steer_rate = 0.3;
    sim.step(0.0);
    assert_eq!(sim.state.lateral_velocity, 0.2);
    assert_eq!(sim.state.yaw_rate, 0.1);
    assert_eq!(sim.state.steer_angle, 0.02);
}

// ---------------------------------------------------------------------------
// Roll
// ---------------------------------------------------------------------------

#[test]
fn the_lateral_tyre_force_rolls_the_bike_upright_from_the_first_step() {
    // In a right-hand turn the inertial force of the lateral acceleration acts to
    // the left of the centre of mass and rolls the bike to the left.  The roll
    // dynamics see the lateral acceleration of the centre of mass (the lateral
    // tyre force over the mass), so the roll rate reacts at once; the yaw rate is
    // still zero in the first step.
    for sign in [1.0_f64, -1.0] {
        let mut sim = MotorcycleDynamicsSimulator::new();
        let speed = 20.0;
        sim.state.speed = speed;
        sim.state.rear_wheel_omega = speed / sim.params.rear_wheel_radius;
        sim.state.front_wheel_omega = speed / sim.params.front_wheel_radius;
        sim.fairing.frontal_area = 0.0;
        sim.state.steer_angle = sign * 0.01;
        sim.step(0.005);
        assert_lateral_state_is_sane(&sim);
        assert!(
            sim.state.lean_rate * sign < 0.0,
            "steer {sign}: lean rate {} rad/s should point away from the turn",
            sim.state.lean_rate
        );
        // upright and without gravity moment the only roll torque is -m a_y h
        let torque = sim.state.lean_rate * sim.params.roll_inertia / 0.005;
        let lateral_accel = -torque / (sim.params.total_mass * sim.params.com_height);
        assert!(
            lateral_accel * sign > 0.0 && lateral_accel.abs() < 5.0,
            "implied lateral acceleration {lateral_accel} m/s^2"
        );
    }
}

// ---------------------------------------------------------------------------
// Robustness
// ---------------------------------------------------------------------------

/// Steering angle that gives roughly 3 m/s^2 of lateral acceleration at `speed`.
fn gentle_steer(speed: f64, wheelbase: f64) -> f64 {
    (3.0 * wheelbase / (speed * speed)).min(0.1)
}

#[test]
fn steering_is_stable_for_coarse_steps_at_every_speed() {
    // The tyre relaxation makes the response of a slow bike a lightly damped
    // oscillation (rate v / sigma), so give it ten seconds to settle.
    let duration = 10.0_f64;
    for speed in [0.5_f64, 1.0, 2.0, 3.5, 5.0, 10.0, 30.0, 60.0] {
        for dt in [0.005_f64, 0.01, 0.02, 0.05] {
            let mut sim = frozen_bike(speed);
            let steer = gentle_steer(speed, sim.params.wheelbase);
            sim.state.steer_angle = steer;
            let steps = (duration / dt).round() as usize;
            let one_second_before_the_end = steps - (1.0 / dt).round() as usize;
            let mut yaw_one_second_before_the_end = 0.0;
            for step in 0..steps {
                sim.step(dt);
                assert_lateral_state_is_sane(&sim);
                if step == one_second_before_the_end {
                    yaw_one_second_before_the_end = sim.state.yaw_rate;
                }
            }
            let expected_yaw = speed * steer / sim.params.wheelbase;
            assert!(
                (sim.state.yaw_rate - expected_yaw).abs() < 0.1 * expected_yaw,
                "v = {speed} m/s, dt = {dt} s: yaw rate {} rad/s, bicycle model {expected_yaw} rad/s",
                sim.state.yaw_rate
            );
            assert!(
                (sim.state.yaw_rate - yaw_one_second_before_the_end).abs() < 0.05 * expected_yaw,
                "v = {speed} m/s, dt = {dt} s: the response has not settled"
            );
            assert!(
                sim.state.lateral_velocity.abs() < 2.0,
                "v = {speed} m/s, dt = {dt} s: lateral velocity {}",
                sim.state.lateral_velocity
            );
        }
    }
}

#[test]
fn low_speed_steering_with_a_coarse_step_does_not_blow_up() {
    // 2 m/s and 50 ms steps: the lateral subsystem is stiffest here.
    let mut sim = frozen_bike(2.0);
    sim.state.steer_angle = 0.1;
    for _ in 0..300 {
        sim.step(0.05);
        assert_lateral_state_is_sane(&sim);
        assert!(
            sim.state.yaw_rate.abs() < 1.0,
            "yaw rate {}",
            sim.state.yaw_rate
        );
    }
    let expected_yaw = 2.0 * 0.1 / sim.params.wheelbase;
    assert!((sim.state.yaw_rate - expected_yaw).abs() < 0.1 * expected_yaw);
}
