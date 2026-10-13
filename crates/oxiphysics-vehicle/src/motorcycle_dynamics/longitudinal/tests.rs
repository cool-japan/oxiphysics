// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the longitudinal wheel/ground coupling of the motorcycle model.

use super::*;

const REAR_RADIUS: f64 = 0.315;
const FRONT_RADIUS: f64 = 0.305;
const REAR_INERTIA: f64 = 0.8;
const FRONT_INERTIA: f64 = 0.6;
const MASS: f64 = 200.0;

// ---------------------------------------------------------------------------
// Slip ratio
// ---------------------------------------------------------------------------

#[test]
fn slip_ratio_is_the_textbook_formula_above_the_floor() {
    let r = 0.31;
    // wheel rim 10 % faster than the ground at 20 m/s
    assert!((slip_ratio(22.0 / r, r, 20.0) - 0.1).abs() < 1e-12);
    // locked wheel
    assert!((slip_ratio(0.0, r, 20.0) + 1.0).abs() < 1e-12);
    // rolling without slip
    assert!(slip_ratio(20.0 / r, r, 20.0).abs() < 1e-12);
}

#[test]
fn slip_ratio_is_finite_and_continuous_at_standstill() {
    let r = 0.31;
    // a spinning wheel at rest has a finite slip ratio: rim speed over the floor
    let at_rest = slip_ratio(5.0 / r, r, 0.0);
    assert!((at_rest - 5.0 / SLIP_SPEED_FLOOR).abs() < 1e-12);
    assert_eq!(slip_ratio(0.0, r, 0.0), 0.0);
    // continuity across the floor
    let omega = 3.0 / r;
    let below = slip_ratio(omega, r, SLIP_SPEED_FLOOR - 1e-9);
    let above = slip_ratio(omega, r, SLIP_SPEED_FLOOR + 1e-9);
    assert!((below - above).abs() < 1e-6);
}

#[test]
fn slip_ratio_derivatives_match_finite_differences() {
    let r = 0.31;
    let h = 1e-6;
    for (omega, speed) in [
        (70.0, 20.0),
        (3.0, 0.2),
        (10.0, 0.9),
        (0.0, 5.0),
        (50.0, 1.5),
        (65.0, 12.0),
        (1.0, 0.0),
    ] {
        let (d_omega, d_speed) = slip_ratio_derivatives(omega, r, speed);
        let fd_omega =
            (slip_ratio(omega + h, r, speed) - slip_ratio(omega - h, r, speed)) / (2.0 * h);
        let fd_speed =
            (slip_ratio(omega, r, speed + h) - slip_ratio(omega, r, speed - h)) / (2.0 * h);
        assert!(
            (d_omega - fd_omega).abs() < 1e-6 * (1.0 + fd_omega.abs()),
            "d/d omega at ({omega}, {speed}): {d_omega} vs {fd_omega}"
        );
        assert!(
            (d_speed - fd_speed).abs() < 1e-6 * (1.0 + fd_speed.abs()),
            "d/d speed at ({omega}, {speed}): {d_speed} vs {fd_speed}"
        );
    }
}

// ---------------------------------------------------------------------------
// Tyre slope
// ---------------------------------------------------------------------------

#[test]
fn longitudinal_slope_is_the_tyre_stiffness_at_small_slip() {
    let tire = PacejkaMotoTire::new();
    let normal = 1000.0;
    // B * C * D * F_z of the default longitudinal Magic Formula coefficients
    let expected = 11.0 * 1.75 * 1.05 * normal;
    let slope = longitudinal_slope(&tire, 0.0, 0.0, normal, 0.0);
    assert!(
        (slope - expected).abs() < 1e-3 * expected,
        "slope {slope} vs stiffness {expected}"
    );
}

#[test]
fn longitudinal_slope_is_zero_beyond_the_force_peak() {
    let tire = PacejkaMotoTire::new();
    let normal = 1000.0;
    // a slip ratio of 0.6 is far past the peak, where the raw slope is negative
    let raw = (tire.forces(0.0, 0.6 + 1e-3, normal, 0.0).1
        - tire.forces(0.0, 0.6 - 1e-3, normal, 0.0).1)
        / 2e-3;
    assert!(raw < 0.0, "raw slope {raw}");
    assert_eq!(longitudinal_slope(&tire, 0.0, 0.6, normal, 0.0), 0.0);
}

#[test]
fn longitudinal_slope_shrinks_with_lateral_force() {
    // under cornering load the friction ellipse leaves less longitudinal grip
    let tire = PacejkaMotoTire::new();
    let straight = longitudinal_slope(&tire, 0.0, 0.02, 1000.0, 0.0);
    let cornering = longitudinal_slope(&tire, 0.1, 0.02, 1000.0, 0.0);
    assert!(straight > 0.0);
    assert!(cornering < straight);
}

#[test]
fn tyre_evaluation_reports_rolling_and_locked_forces() {
    let tire = PacejkaMotoTire::new();
    let normal = 1000.0;
    // rolling slightly faster than the ground: a driving force
    let rolling = TyreLongitudinal::evaluate(&tire, 0.0, normal, 0.0, (22.0 / 0.31, 0.31), 20.0);
    assert!(rolling.force > 0.0 && rolling.slope > 0.0);
    assert_eq!(rolling.force, tire.forces(0.0, 0.1, normal, 0.0).1);
    // the locked wheel at 20 m/s is at slip ratio -1: a large braking force
    assert_eq!(rolling.locked_force, tire.forces(0.0, -1.0, normal, 0.0).1);
    assert!(rolling.locked_force < -0.9 * normal);
    // at standstill a locked wheel does not slip at all
    let at_rest = TyreLongitudinal::evaluate(&tire, 0.0, normal, 0.0, (0.0, 0.31), 0.0);
    assert_eq!(at_rest.force, 0.0);
    assert_eq!(at_rest.locked_force, 0.0);
}

// ---------------------------------------------------------------------------
// advance
// ---------------------------------------------------------------------------

fn tyre(force: f64, slope: f64) -> TyreLongitudinal {
    TyreLongitudinal {
        force,
        slope,
        locked_force: force,
        locked_slope: slope,
    }
}

fn rear_wheel(omega: f64, tyre: TyreLongitudinal) -> WheelInput {
    WheelInput {
        radius: REAR_RADIUS,
        inertia: REAR_INERTIA,
        omega,
        drive_torque: 0.0,
        brake_torque: 0.0,
        tyre,
    }
}

fn front_wheel(omega: f64, tyre: TyreLongitudinal) -> WheelInput {
    WheelInput {
        radius: FRONT_RADIUS,
        inertia: FRONT_INERTIA,
        omega,
        drive_torque: 0.0,
        brake_torque: 0.0,
        tyre,
    }
}

fn chassis(speed: f64, drag: f64) -> ChassisInput {
    ChassisInput {
        mass: MASS,
        speed,
        drag,
    }
}

/// A tyre that exerts no force whatever the slip.
fn dead_tyre() -> TyreLongitudinal {
    tyre(0.0, 0.0)
}

#[test]
fn step_agrees_with_explicit_euler_for_tiny_steps() {
    let drag = 51.2;
    let rear = WheelInput {
        drive_torque: 224.0,
        ..rear_wheel(20.0 / REAR_RADIUS, tyre(400.0, 18_000.0))
    };
    let front = front_wheel(20.0 / FRONT_RADIUS, tyre(-30.0, 21_000.0));
    let dt = 1e-7;
    let (state, accel) = advance(dt, chassis(20.0, drag), rear, front);

    let f_v = (400.0 - 30.0 - drag) / MASS;
    let f_r = (224.0 - REAR_RADIUS * 400.0) / REAR_INERTIA;
    let f_f = FRONT_RADIUS * 30.0 / FRONT_INERTIA;
    assert!((accel - f_v).abs() < 1e-4 * f_v.abs());
    assert!((state.speed - (20.0 + dt * f_v)).abs() < 1e-10);
    assert!((state.rear_omega - (rear.omega + dt * f_r)).abs() < 1e-9);
    assert!((state.front_omega - (front.omega + dt * f_f)).abs() < 1e-9);
}

#[test]
fn drive_torque_reaches_the_bike_only_through_the_tyre_force() {
    // With a vanishing tyre force the wheel spins up and the bike does not move.
    let rear = WheelInput {
        drive_torque: 224.0,
        ..rear_wheel(0.0, dead_tyre())
    };
    let front = front_wheel(0.0, dead_tyre());
    let (state, accel) = advance(0.01, chassis(0.0, 0.0), rear, front);
    assert_eq!(state.speed, 0.0);
    assert_eq!(accel, 0.0);
    assert!((state.rear_omega - 0.01 * 224.0 / REAR_INERTIA).abs() < 1e-12);
    assert_eq!(state.front_omega, 0.0);
}

#[test]
fn tyre_reaction_slows_a_driven_wheel_and_spins_up_a_braked_one() {
    let omega = 20.0 / REAR_RADIUS;
    let front = front_wheel(20.0 / FRONT_RADIUS, dead_tyre());
    let dt = 0.01;
    // driving force: the wheel feels a retarding torque
    let driving = rear_wheel(omega, tyre(300.0, 0.0));
    let (state, _) = advance(dt, chassis(20.0, 0.0), driving, front);
    assert!((state.rear_omega - (omega - dt * REAR_RADIUS * 300.0 / REAR_INERTIA)).abs() < 1e-9);
    // braking force: the wheel is spun up by the ground
    let braking = rear_wheel(omega, tyre(-300.0, 0.0));
    let (state, _) = advance(dt, chassis(20.0, 0.0), braking, front);
    assert!((state.rear_omega - (omega + dt * REAR_RADIUS * 300.0 / REAR_INERTIA)).abs() < 1e-9);
}

#[test]
fn step_stays_bounded_for_huge_steps_and_stiff_tyres() {
    let rear = WheelInput {
        drive_torque: 224.0,
        brake_torque: 50.0,
        ..rear_wheel(1.5, tyre(400.0, 18_000.0))
    };
    let front = WheelInput {
        brake_torque: 50.0,
        ..front_wheel(1.6, tyre(-30.0, 21_000.0))
    };
    for dt in [0.01, 0.1, 1.0, 10.0] {
        let (state, accel) = advance(dt, chassis(0.5, 0.0), rear, front);
        assert!(state.speed.is_finite() && state.speed >= 0.0, "dt {dt}");
        assert!(
            state.rear_omega.is_finite() && state.rear_omega >= 0.0,
            "dt {dt}"
        );
        assert!(
            state.front_omega.is_finite() && state.front_omega >= 0.0,
            "dt {dt}"
        );
        assert!(accel.is_finite(), "dt {dt}");
        assert!(state.rear_omega * REAR_RADIUS < 1e3, "dt {dt}");
    }
}

#[test]
fn free_braked_wheel_loses_exactly_the_brake_impulse() {
    // a dead tyre leaves only the brake: dt * T / I = 1.25 rad/s per step
    let rear = WheelInput {
        brake_torque: 100.0,
        ..rear_wheel(50.0, dead_tyre())
    };
    let (state, _) = advance(
        0.01,
        chassis(20.0, 0.0),
        rear,
        front_wheel(65.0, dead_tyre()),
    );
    assert!((state.rear_omega - 48.75).abs() < 1e-12);
}

#[test]
fn brake_overshoot_stops_the_wheel_without_reversing_it() {
    // The wheel runs faster than the ground (driving force 500 N retards it by
    // 157 N m) and is braked on top: the sum would take it below zero.
    let rear = WheelInput {
        brake_torque: 100.0,
        ..rear_wheel(2.0, tyre(500.0, 0.0))
    };
    let (state, _) = advance(
        0.01,
        chassis(20.0, 0.0),
        rear,
        front_wheel(65.0, dead_tyre()),
    );
    assert_eq!(state.rear_omega, 0.0);
}

#[test]
fn brake_holds_a_wheel_against_the_drive_torque() {
    // Bike at rest, rear brake (200 N m) stronger than the drive torque (56 N m):
    // the wheel stays at zero and the bike does not creep.
    let locked = tyre(0.0, 18_000.0);
    let rear = WheelInput {
        drive_torque: 56.0,
        brake_torque: 200.0,
        ..rear_wheel(0.0, locked)
    };
    let front = WheelInput {
        brake_torque: 200.0,
        ..front_wheel(0.0, locked)
    };
    let (state, accel) = advance(0.01, chassis(0.0, 0.0), rear, front);
    assert_eq!(state.rear_omega, 0.0);
    assert_eq!(state.front_omega, 0.0);
    assert_eq!(state.speed, 0.0);
    assert_eq!(accel, 0.0);
}

#[test]
fn drive_torque_beyond_the_brake_spins_the_wheel_up() {
    let rear = WheelInput {
        drive_torque: 300.0,
        brake_torque: 200.0,
        ..rear_wheel(0.0, tyre(0.0, 0.0))
    };
    let (state, _) = advance(0.01, chassis(0.0, 0.0), rear, front_wheel(0.0, dead_tyre()));
    assert!((state.rear_omega - 0.01 * 100.0 / REAR_INERTIA).abs() < 1e-12);
}

#[test]
fn ground_reaction_beyond_the_brake_turns_a_slow_wheel() {
    // Bike moving at 10 m/s, locked tyre force -1000 N: the ground turns a locked
    // wheel forward with 315 N m.  A 100 N m brake cannot hold it ...
    let ground = TyreLongitudinal {
        force: -1000.0,
        slope: 0.0,
        locked_force: -1000.0,
        locked_slope: 0.0,
    };
    let weak = WheelInput {
        brake_torque: 100.0,
        ..rear_wheel(0.0, ground)
    };
    let (state, _) = advance(
        0.01,
        chassis(10.0, 0.0),
        weak,
        front_wheel(31.0, dead_tyre()),
    );
    let expected = 0.01 * (315.0 - 100.0) / REAR_INERTIA;
    assert!(
        (state.rear_omega - expected).abs() < 1e-9,
        "{}",
        state.rear_omega
    );
    // ... but a 400 N m brake can, and the bike then sees the locked tyre force.
    let strong = WheelInput {
        brake_torque: 400.0,
        ..rear_wheel(0.0, ground)
    };
    let (state, accel) = advance(
        0.01,
        chassis(10.0, 0.0),
        strong,
        front_wheel(31.0, dead_tyre()),
    );
    assert_eq!(state.rear_omega, 0.0);
    assert!((accel - (-1000.0 / MASS)).abs() < 1e-9, "accel {accel}");
}

#[test]
fn fast_wheel_is_not_held_even_by_a_strong_brake() {
    // 50 rad/s cannot be removed within one step by a 100 N m brake (1.25 rad/s)
    let rear = WheelInput {
        brake_torque: 100.0,
        ..rear_wheel(50.0, dead_tyre())
    };
    assert!(!rear.is_held(0.01));
    let slow = WheelInput { omega: 1.0, ..rear };
    assert!(slow.is_held(0.01));
    // no brake, no hold; no time, no hold
    let unbraked = WheelInput {
        brake_torque: 0.0,
        ..slow
    };
    assert!(!unbraked.is_held(0.01));
    assert!(!slow.is_held(0.0));
}

#[test]
fn braking_bike_comes_to_rest() {
    // dead tyres leave a drag of 0.01 N on a bike that creeps at 0.5 mm/s
    let rear = rear_wheel(0.0, dead_tyre());
    let front = front_wheel(0.0, dead_tyre());
    let (state, accel) = advance(0.01, chassis(5.0e-4, 0.01), rear, front);
    assert_eq!(state.speed, 0.0);
    assert!(accel < 0.0);
}

#[test]
fn accelerating_bike_that_is_still_slow_keeps_moving() {
    // A bike pulling away from rest must not be snapped back to zero.
    let rear = WheelInput {
        drive_torque: 224.0,
        ..rear_wheel(0.0, tyre(0.0, 18_000.0))
    };
    let front = front_wheel(0.0, tyre(0.0, 21_000.0));
    for dt in [1e-4, 1e-3, 1e-2] {
        let (state, accel) = advance(dt, chassis(0.0, 0.0), rear, front);
        assert!(state.speed > 0.0, "dt {dt}: speed {}", state.speed);
        assert!(accel > 0.0, "dt {dt}");
    }
}

#[test]
fn zero_time_step_changes_nothing() {
    let rear = WheelInput {
        drive_torque: 224.0,
        brake_torque: 10.0,
        ..rear_wheel(60.0, tyre(300.0, 18_000.0))
    };
    let front = front_wheel(62.0, tyre(-20.0, 21_000.0));
    let (state, accel) = advance(0.0, chassis(19.0, 51.0), rear, front);
    assert_eq!(state.speed, 19.0);
    assert_eq!(state.rear_omega, 60.0);
    assert_eq!(state.front_omega, 62.0);
    assert!((accel - (300.0 - 20.0 - 51.0) / MASS).abs() < 1e-12);
}
