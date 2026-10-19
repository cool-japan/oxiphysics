// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the lateral dynamics of the motorcycle model.

use super::*;

const MASS: f64 = 200.0;
const YAW_INERTIA: f64 = 20.0;
/// Distance from the centre of mass to the front axle (m).
const A: f64 = 0.65;
/// Distance from the centre of mass to the rear axle (m).
const B: f64 = 0.78;
const GRAVITY: f64 = 9.81;

const AT_REST: LateralState = LateralState {
    lateral_velocity: 0.0,
    yaw_rate: 0.0,
};

fn chassis() -> Chassis {
    Chassis {
        mass: MASS,
        yaw_inertia: YAW_INERTIA,
        a: A,
        b: B,
    }
}

/// Static axle loads (N): the front axle carries the share `b / L` of the weight.
fn static_loads() -> (f64, f64) {
    let wheelbase = A + B;
    (
        MASS * GRAVITY * B / wheelbase,
        MASS * GRAVITY * A / wheelbase,
    )
}

fn input(speed: f64, steer: f64) -> LateralInput {
    let (front_normal, rear_normal) = static_loads();
    LateralInput {
        speed,
        steer,
        lean: 0.0,
        front_normal,
        rear_normal,
        front_slip_ratio: 0.0,
        rear_slip_ratio: 0.0,
    }
}

fn fresh_tires() -> (PacejkaMotoTire, PacejkaMotoTire) {
    (PacejkaMotoTire::new(), PacejkaMotoTire::new())
}

// ---------------------------------------------------------------------------
// Sub-stepping
// ---------------------------------------------------------------------------

#[test]
fn cornering_stiffness_is_the_initial_slope_of_the_magic_formula() {
    let tire = PacejkaMotoTire::new();
    let normal = 1000.0;
    let alpha = 2.0e-3;
    let secant = tire.forces(alpha, 0.0, normal, 0.0).0 / alpha;
    let stiffness = cornering_stiffness(&tire, normal);
    assert!((stiffness - tire.lat.b * tire.lat.c * tire.lat.d * normal).abs() < 1e-9);
    assert!(
        (secant - stiffness).abs() < 2.0e-3 * stiffness,
        "secant {secant} N/rad, initial slope {stiffness} N/rad"
    );
    assert_eq!(cornering_stiffness(&tire, 0.0), 0.0);
    assert_eq!(cornering_stiffness(&tire, -5.0), 0.0);
}

#[test]
fn substep_count_follows_the_step_length() {
    let (front, rear) = fresh_tires();
    let count = |dt: f64| substep_count(dt, chassis(), input(20.0, 0.0), &front, &rear);
    assert_eq!(count(0.0), 1);
    assert_eq!(count(1.0e-4), 1);
    // natural frequency about 54 rad/s: a 10 ms step needs a couple of sub-steps,
    // a second about a hundred
    assert!((1..=3).contains(&count(0.01)), "{}", count(0.01));
    let one_second = count(1.0);
    assert!((90..=130).contains(&one_second), "{one_second}");
    assert!(count(0.05) >= count(0.01));
    assert_eq!(count(1.0e9), MAX_SUBSTEPS);
    assert_eq!(count(f64::NAN), MAX_SUBSTEPS);
}

#[test]
fn stiffer_tyres_need_more_substeps() {
    let (soft_front, soft_rear) = fresh_tires();
    let (mut stiff_front, mut stiff_rear) = fresh_tires();
    stiff_front.lat.b *= 4.0;
    stiff_rear.lat.b *= 4.0;
    let soft = substep_count(0.05, chassis(), input(20.0, 0.0), &soft_front, &soft_rear);
    let stiff = substep_count(0.05, chassis(), input(20.0, 0.0), &stiff_front, &stiff_rear);
    assert!(stiff > soft, "{stiff} sub-steps vs {soft}");
    // a longer relaxation length lowers the natural frequency
    let (mut long_front, mut long_rear) = fresh_tires();
    long_front.relaxation_length *= 4.0;
    long_rear.relaxation_length *= 4.0;
    let relaxed = substep_count(0.05, chassis(), input(20.0, 0.0), &long_front, &long_rear);
    assert!(relaxed < soft, "{relaxed} sub-steps vs {soft}");
}

// ---------------------------------------------------------------------------
// Degenerate steps
// ---------------------------------------------------------------------------

#[test]
fn a_step_of_zero_length_changes_nothing() {
    let (mut front, mut rear) = fresh_tires();
    front.filtered_slip_angle = 0.05;
    rear.filtered_slip_angle = -0.02;
    let state = LateralState {
        lateral_velocity: 0.3,
        yaw_rate: 0.2,
    };
    for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        let (next, accel) = advance(
            dt,
            chassis(),
            input(20.0, 0.03),
            state,
            &mut front,
            &mut rear,
        );
        assert_eq!(next.lateral_velocity, 0.3, "dt = {dt}");
        assert_eq!(next.yaw_rate, 0.2, "dt = {dt}");
        assert_eq!(front.filtered_slip_angle, 0.05, "dt = {dt}");
        assert_eq!(rear.filtered_slip_angle, -0.02, "dt = {dt}");
        assert!(accel.is_finite());
    }
}

#[test]
fn below_the_grip_speed_the_bike_follows_the_kinematic_motion() {
    let wheelbase = A + B;
    let moving = LateralState {
        lateral_velocity: 0.3,
        yaw_rate: 0.4,
    };
    // at rest nothing can rotate or slide, and the tyres lose their slip at once
    let (mut front, mut rear) = fresh_tires();
    front.filtered_slip_angle = 0.1;
    rear.filtered_slip_angle = -0.1;
    let (rest, accel) = advance(
        0.01,
        chassis(),
        input(0.0, 0.2),
        moving,
        &mut front,
        &mut rear,
    );
    assert_eq!(rest.yaw_rate, 0.0);
    assert_eq!(rest.lateral_velocity, 0.0);
    assert_eq!(accel, 0.0);
    assert_eq!(front.filtered_slip_angle, 0.0);
    assert_eq!(rear.filtered_slip_angle, 0.0);
    // crawling at the grip speed: r = v δ / L and v_y = b r
    let (speed, steer) = (GRIP_SPEED, 0.1);
    let (crawl, accel) = advance(
        0.01,
        chassis(),
        input(speed, steer),
        moving,
        &mut front,
        &mut rear,
    );
    assert!((crawl.yaw_rate - speed * steer / wheelbase).abs() < 1e-15);
    assert!((crawl.lateral_velocity - B * crawl.yaw_rate).abs() < 1e-15);
    assert!((accel - speed * crawl.yaw_rate).abs() < 1e-15);
    // just above it the slip-angle dynamics take over (and keep the lateral velocity)
    let (above, _) = advance(
        0.01,
        chassis(),
        input(2.0 * GRIP_SPEED, steer),
        moving,
        &mut front,
        &mut rear,
    );
    assert!(
        (above.lateral_velocity - B * above.yaw_rate).abs() > 1e-3,
        "the dynamic model must not snap to the kinematic state"
    );
}

#[test]
fn huge_steps_stay_finite_and_bounded() {
    let wheelbase = A + B;
    for dt in [0.5, 2.0, 10.0] {
        for speed in [0.5_f64, 5.0, 40.0] {
            // gentle steering: about 3 m/s^2 of lateral acceleration
            let steer = (3.0 * wheelbase / (speed * speed)).min(0.1);
            let (mut front, mut rear) = fresh_tires();
            let (next, accel) = advance(
                dt,
                chassis(),
                input(speed, steer),
                AT_REST,
                &mut front,
                &mut rear,
            );
            assert!(
                next.yaw_rate.is_finite() && next.lateral_velocity.is_finite() && accel.is_finite(),
                "dt = {dt} s, v = {speed} m/s: {next:?}"
            );
            // A slow bike answers a steering step with a lightly damped oscillation
            // around the neutral-steer yaw rate v * delta / L, whose amplitude is
            // about the size of the step.
            let expected_yaw = speed * steer / wheelbase;
            assert!(
                next.yaw_rate > -0.3 * expected_yaw && next.yaw_rate < 2.5 * expected_yaw,
                "dt = {dt} s, v = {speed} m/s: yaw rate {} rad/s vs {expected_yaw} rad/s",
                next.yaw_rate
            );
            assert!(
                next.lateral_velocity.abs() < 5.0,
                "dt = {dt} s, v = {speed} m/s: {next:?}"
            );
            // ten seconds are long enough for every response to settle
            if dt >= 10.0 {
                assert!(
                    (next.yaw_rate - expected_yaw).abs() < 0.1 * expected_yaw,
                    "dt = {dt} s, v = {speed} m/s: yaw rate {} rad/s vs {expected_yaw} rad/s",
                    next.yaw_rate
                );
            }
        }
    }
}

#[test]
fn steering_far_beyond_the_grip_limit_stays_within_what_the_tyres_can_give() {
    // 0.3 rad of steer at speed asks for far more lateral force than the tyres can
    // give.  Whatever the slip angles do, each tyre force is at most D * N, so the
    // lateral acceleration is at most D * g and the yaw acceleration at most
    // D * (a N_f + b N_r) / I_z.
    let (front_normal, rear_normal) = static_loads();
    let tire = PacejkaMotoTire::new();
    let max_accel = tire.lat.d * GRAVITY;
    let max_yaw_accel = tire.lat.d * (A * front_normal + B * rear_normal) / YAW_INERTIA;
    for speed in [5.0, 40.0] {
        let (mut front, mut rear) = fresh_tires();
        let dt = 2.0;
        let (next, accel) = advance(
            dt,
            chassis(),
            input(speed, 0.3),
            AT_REST,
            &mut front,
            &mut rear,
        );
        assert!(
            next.yaw_rate.is_finite() && next.lateral_velocity.is_finite() && accel.is_finite(),
            "v = {speed} m/s: {next:?}"
        );
        assert!(
            accel.abs() <= max_accel + 1e-9,
            "v = {speed} m/s: lateral acceleration {accel} m/s^2 exceeds {max_accel} m/s^2"
        );
        assert!(
            next.yaw_rate.abs() <= max_yaw_accel * dt + 1e-9,
            "v = {speed} m/s: yaw rate {} rad/s exceeds {}",
            next.yaw_rate,
            max_yaw_accel * dt
        );
    }
}

// ---------------------------------------------------------------------------
// Response
// ---------------------------------------------------------------------------

#[test]
fn a_right_steer_pushes_the_bike_to_the_right_and_yaws_it_right() {
    let (mut front, mut rear) = fresh_tires();
    let (next, accel) = advance(
        0.01,
        chassis(),
        input(20.0, 0.05),
        AT_REST,
        &mut front,
        &mut rear,
    );
    assert!(next.yaw_rate > 0.0, "yaw rate {}", next.yaw_rate);
    assert!(next.lateral_velocity > 0.0, "v_y {}", next.lateral_velocity);
    assert!(accel > 0.0, "lateral acceleration {accel}");
    assert!(front.filtered_slip_angle > 0.0);
    // the mirror image
    let (mut mirror_front, mut mirror_rear) = fresh_tires();
    let (mirror, mirror_accel) = advance(
        0.01,
        chassis(),
        input(20.0, -0.05),
        AT_REST,
        &mut mirror_front,
        &mut mirror_rear,
    );
    assert!((mirror.yaw_rate + next.yaw_rate).abs() < 1e-15);
    assert!((mirror.lateral_velocity + next.lateral_velocity).abs() < 1e-15);
    assert!((mirror_accel + accel).abs() < 1e-15);
}

#[test]
fn lateral_velocity_is_driven_by_the_lateral_force_balance() {
    // First sub-step from rest: only the front tyre has built up slip, so
    // m v_y' = F_f cos(delta) and I_z r' = a F_f cos(delta), whatever the step.
    let (mut front, mut rear) = fresh_tires();
    let dt = 1.0e-4;
    assert_eq!(
        substep_count(dt, chassis(), input(20.0, 0.05), &front, &rear),
        1
    );
    let steer = 0.05_f64;
    let (next, accel) = advance(
        dt,
        chassis(),
        input(20.0, steer),
        AT_REST,
        &mut front,
        &mut rear,
    );
    let front_force = front
        .forces(front.filtered_slip_angle, 0.0, static_loads().0, 0.0)
        .0;
    let rear_force = rear
        .forces(rear.filtered_slip_angle, 0.0, static_loads().1, 0.0)
        .0;
    // the rear tyre has seen no slip yet (a yaw rate and lateral velocity of zero)
    assert!(rear_force.abs() < 1e-12, "rear force {rear_force} N");
    let side_force = front_force * steer.cos();
    assert!((accel - side_force / MASS).abs() < 1e-9 * side_force.abs().max(1.0));
    assert!((next.lateral_velocity - dt * side_force / MASS).abs() < 1e-12);
    assert!((next.yaw_rate - dt * A * side_force / YAW_INERTIA).abs() < 1e-12);
}

#[test]
fn neutral_steer_steady_state_matches_the_linear_bicycle_model() {
    let (mut front, mut rear) = fresh_tires();
    let (speed, steer) = (20.0, 0.005);
    let wheelbase = A + B;
    let mut state = AT_REST;
    let mut accel = 0.0;
    for _ in 0..600 {
        let (next, mean_accel) = advance(
            0.005,
            chassis(),
            input(speed, steer),
            state,
            &mut front,
            &mut rear,
        );
        state = next;
        accel = mean_accel;
    }
    // Cornering stiffness proportional to the axle load: neutral steer, r = v δ / L.
    let expected_yaw = speed * steer / wheelbase;
    assert!(
        (state.yaw_rate - expected_yaw).abs() < 0.01 * expected_yaw,
        "yaw rate {} rad/s, bicycle model {expected_yaw} rad/s",
        state.yaw_rate
    );
    // The rear axle slips with α_r = v r / (κ g): v_y = b r − v α_r = r (b − v² / (κ g)).
    let kappa = front.lat.b * front.lat.c * front.lat.d;
    let expected_v_y = expected_yaw * (B - speed * speed / (kappa * GRAVITY));
    assert!(expected_v_y < 0.0, "the tail slides outwards at this speed");
    assert!(
        (state.lateral_velocity - expected_v_y).abs() < 0.03 * expected_v_y.abs(),
        "lateral velocity {} m/s, bicycle model {expected_v_y} m/s",
        state.lateral_velocity
    );
    // steady cornering: the lateral acceleration of the centre of mass is v r
    assert!(
        (accel - speed * state.yaw_rate).abs() < 0.02 * speed * state.yaw_rate.abs(),
        "lateral acceleration {accel} m/s^2, v r = {}",
        speed * state.yaw_rate
    );
}

#[test]
fn an_unstable_oversteering_bike_diverges_instead_of_being_damped() {
    // A rear axle much softer than the front one (a·C_f > b·C_r) is unstable
    // above its critical speed.  The explicit sub-steps must follow the growing
    // yaw rate rather than damp it.
    let (mut front, mut rear) = fresh_tires();
    rear.lat.b *= 0.15;
    let speed = 30.0;
    let mut state = LateralState {
        lateral_velocity: 0.0,
        yaw_rate: 0.01,
    };
    let mut peak: f64 = 0.0;
    for _ in 0..200 {
        let (next, _) = advance(
            0.005,
            chassis(),
            input(speed, 0.0),
            state,
            &mut front,
            &mut rear,
        );
        state = next;
        peak = peak.max(state.yaw_rate.abs());
        assert!(state.yaw_rate.is_finite() && state.lateral_velocity.is_finite());
    }
    assert!(
        peak > 0.1,
        "the initial 0.01 rad/s disturbance only reached {peak} rad/s"
    );
}
