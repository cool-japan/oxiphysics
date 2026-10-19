// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the longitudinal and yaw dynamics of the single-track
//! motorcycle simulator.
//!
//! The simulator used to add the drive and brake forces to the chassis on top of
//! the tyre forces they already produce, and the wheel reaction used `|F_x|`, so
//! the tyre reaction always slowed the wheel down: braking drove the wheels into
//! lock-up.  Wheels and bike are now coupled through the signed tyre force only.
//! The lateral slip angles and the yaw moment used swapped axle lever arms, which
//! made the yaw dynamics unstable.

use oxiphysics_vehicle::motorcycle_dynamics::MotorcycleDynamicsSimulator;

const DT: f64 = 0.01;
const AIR_DENSITY: f64 = 1.225;

/// A bike rolling freely at `speed` (m/s), upright, no inputs.
fn rolling_bike(speed: f64) -> MotorcycleDynamicsSimulator {
    let mut sim = MotorcycleDynamicsSimulator::new();
    sim.state.speed = speed;
    sim.state.rear_wheel_omega = speed / sim.params.rear_wheel_radius;
    sim.state.front_wheel_omega = speed / sim.params.front_wheel_radius;
    sim
}

/// Slip ratio `(omega r - v) / max(v, 1 m/s)` of the rear wheel.
fn rear_slip(sim: &MotorcycleDynamicsSimulator) -> f64 {
    let v = sim.state.speed;
    (sim.state.rear_wheel_omega * sim.params.rear_wheel_radius - v) / v.max(1.0)
}

/// Slip ratio `(omega r - v) / max(v, 1 m/s)` of the front wheel.
fn front_slip(sim: &MotorcycleDynamicsSimulator) -> f64 {
    let v = sim.state.speed;
    (sim.state.front_wheel_omega * sim.params.front_wheel_radius - v) / v.max(1.0)
}

/// Mass of the bike plus the translational equivalent of the wheel inertias,
/// which the bike has to spin up or down as well when it changes speed.
fn effective_mass(sim: &MotorcycleDynamicsSimulator) -> f64 {
    let p = &sim.params;
    p.total_mass
        + p.rear_wheel_inertia / (p.rear_wheel_radius * p.rear_wheel_radius)
        + p.front_wheel_inertia / (p.front_wheel_radius * p.front_wheel_radius)
}

/// Aerodynamic drag of the default sport-bike fairing at zero pitch (N).
fn drag(sim: &MotorcycleDynamicsSimulator, speed: f64) -> f64 {
    let f = &sim.fairing;
    0.5 * AIR_DENSITY * speed * speed * f.frontal_area * f.cd_zero
}

fn assert_state_is_sane(sim: &MotorcycleDynamicsSimulator) {
    let s = &sim.state;
    for (name, value) in [
        ("speed", s.speed),
        ("rear wheel", s.rear_wheel_omega),
        ("front wheel", s.front_wheel_omega),
        ("yaw rate", s.yaw_rate),
        ("lean", s.lean_angle),
        ("front load", s.front_normal),
        ("rear load", s.rear_normal),
    ] {
        assert!(value.is_finite(), "{name} is not finite: {value}");
    }
    assert!(s.speed >= 0.0, "negative speed {}", s.speed);
    assert!(s.rear_wheel_omega >= 0.0, "reversed rear wheel");
    assert!(s.front_wheel_omega >= 0.0, "reversed front wheel");
}

// ---------------------------------------------------------------------------
// Driving
// ---------------------------------------------------------------------------

#[test]
fn constant_drive_torque_accelerates_the_bike_from_rest() {
    let mut sim = MotorcycleDynamicsSimulator::new();
    sim.set_engine_torque(80.0);
    let mut previous = 0.0;
    for step in 0..500 {
        sim.step(DT);
        assert_state_is_sane(&sim);
        assert!(
            sim.state.speed >= previous,
            "speed fell from {previous} to {} at step {step}",
            sim.state.speed
        );
        previous = sim.state.speed;
    }
    // About 3 m/s^2 for five seconds.
    assert!(
        sim.state.speed > 10.0,
        "only {} m/s after 5 s of full drive torque",
        sim.state.speed
    );
    let slip = rear_slip(&sim);
    assert!(
        slip > 0.0 && slip < 0.2,
        "rear slip ratio {slip} did not settle below 0.2"
    );
}

#[test]
fn rear_slip_settles_below_the_peak_under_constant_drive_torque() {
    let mut sim = rolling_bike(20.0);
    sim.set_engine_torque(80.0);
    for step in 0..300 {
        sim.step(DT);
        assert_state_is_sane(&sim);
        if step >= 30 {
            let slip = rear_slip(&sim);
            assert!(
                slip > 0.0 && slip < 0.2,
                "rear slip ratio {slip} at step {step}"
            );
        }
    }
}

#[test]
fn drive_torque_is_not_counted_twice() {
    // The wheel torque T pushes the bike with T / r through the tyre only.  The
    // rotating wheels add I / r^2 to the mass that has to be accelerated.
    let mut sim = rolling_bike(20.0);
    sim.set_engine_torque(80.0);
    let wheel_torque = 80.0 * 42.0 / 15.0; // engine torque x sprocket ratio
    let t_end = 0.5;
    let steps = (t_end / DT) as usize;
    let v0 = sim.state.speed;
    for _ in 0..steps {
        sim.step(DT);
    }
    let measured = (sim.state.speed - v0) / t_end;
    let mean_speed = 0.5 * (v0 + sim.state.speed);
    let expected = (wheel_torque / sim.params.rear_wheel_radius - drag(&sim, mean_speed))
        / effective_mass(&sim);
    assert!(
        (measured - expected).abs() < 0.05 * expected,
        "acceleration {measured} m/s^2, expected {expected} m/s^2"
    );
}

#[test]
fn free_rolling_bike_only_loses_speed_to_drag() {
    let mut sim = rolling_bike(20.0);
    let t_end = 1.0;
    let steps = (t_end / DT) as usize;
    let mut previous = sim.state.speed;
    for _ in 0..steps {
        sim.step(DT);
        assert_state_is_sane(&sim);
        assert!(sim.state.speed <= previous);
        previous = sim.state.speed;
    }
    let measured = (20.0 - sim.state.speed) / t_end;
    let expected = drag(&sim, 20.0) / effective_mass(&sim);
    assert!(
        (measured - expected).abs() < 0.1 * expected,
        "deceleration {measured} m/s^2, expected {expected} m/s^2"
    );
    // The wheels keep rolling with the bike.
    assert!(rear_slip(&sim).abs() < 0.01);
    assert!(front_slip(&sim).abs() < 0.01);
}

#[test]
fn launch_does_not_depend_on_the_step_size() {
    let mut final_speeds = Vec::new();
    for dt in [0.001_f64, 0.005, 0.01, 0.02] {
        let mut sim = MotorcycleDynamicsSimulator::new();
        sim.set_engine_torque(80.0);
        let steps = (3.0 / dt).round() as usize;
        let mut previous = 0.0;
        for step in 0..steps {
            sim.step(dt);
            assert_state_is_sane(&sim);
            assert!(
                sim.state.speed >= previous,
                "dt {dt}: speed fell from {previous} to {} at step {step}",
                sim.state.speed
            );
            previous = sim.state.speed;
        }
        final_speeds.push(sim.state.speed);
    }
    let reference = final_speeds[0];
    for (index, speed) in final_speeds.iter().enumerate() {
        assert!(
            (speed - reference).abs() < 0.03 * reference,
            "step size #{index}: {speed} m/s vs {reference} m/s"
        );
    }
}

#[test]
fn zero_time_step_leaves_the_motion_untouched() {
    let mut sim = rolling_bike(15.0);
    sim.set_engine_torque(50.0);
    let before = sim.state.clone();
    sim.step(0.0);
    assert_eq!(sim.state.speed, before.speed);
    assert_eq!(sim.state.rear_wheel_omega, before.rear_wheel_omega);
    assert_eq!(sim.state.front_wheel_omega, before.front_wheel_omega);
    assert_state_is_sane(&sim);
}

// ---------------------------------------------------------------------------
// Braking
// ---------------------------------------------------------------------------

#[test]
fn moderate_braking_does_not_lock_the_wheels() {
    let mut sim = rolling_bike(20.0);
    sim.set_brakes(100.0, 40.0);
    for _ in 0..10 {
        sim.step(DT);
        assert_state_is_sane(&sim);
    }
    let v = sim.state.speed;
    assert!(v < 20.0 && v > 19.5, "speed after 0.1 s: {v}");
    let rim_front = sim.state.front_wheel_omega * sim.params.front_wheel_radius;
    let rim_rear = sim.state.rear_wheel_omega * sim.params.rear_wheel_radius;
    assert!(
        rim_front > 0.9 * v,
        "front wheel nearly locked: {rim_front} m/s"
    );
    assert!(
        rim_rear > 0.9 * v,
        "rear wheel nearly locked: {rim_rear} m/s"
    );
    assert!(front_slip(&sim) < 0.0, "braking slip is negative");
}

#[test]
fn moderate_braking_decelerates_with_the_brake_force() {
    let mut sim = rolling_bike(20.0);
    let (front_torque, rear_torque) = (100.0, 40.0);
    sim.set_brakes(front_torque, rear_torque);
    // let the slip build up, then measure over a second
    for _ in 0..30 {
        sim.step(DT);
    }
    let v_start = sim.state.speed;
    for _ in 0..100 {
        sim.step(DT);
        assert_state_is_sane(&sim);
    }
    let measured = (v_start - sim.state.speed) / 1.0;
    let mean_speed = 0.5 * (v_start + sim.state.speed);
    let brake_force =
        front_torque / sim.params.front_wheel_radius + rear_torque / sim.params.rear_wheel_radius;
    let expected = (brake_force + drag(&sim, mean_speed)) / effective_mass(&sim);
    assert!(
        (measured - expected).abs() < 0.05 * expected,
        "deceleration {measured} m/s^2, expected {expected} m/s^2"
    );
}

#[test]
fn heavy_front_braking_locks_the_front_wheel_but_the_bike_keeps_going() {
    // 800 N m asks for 2.6 kN at the ground, more than the front tyre can give
    let mut sim = rolling_bike(20.0);
    sim.set_brakes(800.0, 0.0);
    for _ in 0..50 {
        sim.step(DT);
        assert_state_is_sane(&sim);
    }
    assert_eq!(
        sim.state.front_wheel_omega, 0.0,
        "front wheel should be locked"
    );
    assert!(sim.state.speed < 18.0 && sim.state.speed > 0.0);
    // The unbraked rear wheel keeps rolling.
    assert!(rear_slip(&sim).abs() < 0.1);
}

#[test]
fn braking_to_a_standstill_never_reverses_anything() {
    let mut sim = rolling_bike(15.0);
    sim.set_brakes(300.0, 200.0);
    let mut stopped_at = None;
    for step in 0..1000 {
        sim.step(DT);
        assert_state_is_sane(&sim);
        if stopped_at.is_none() && sim.state.speed == 0.0 {
            stopped_at = Some(step);
        }
        if stopped_at.is_some() {
            assert_eq!(sim.state.speed, 0.0, "bike rolled again at step {step}");
            assert_eq!(sim.state.rear_wheel_omega, 0.0);
            assert_eq!(sim.state.front_wheel_omega, 0.0);
        }
    }
    let stopped_at = stopped_at.expect("the bike should have stopped within 10 s");
    // Full braking from 15 m/s needs less than about 3 s.
    assert!(stopped_at < 400, "stopped only after {stopped_at} steps");
}

#[test]
fn brakes_hold_a_stationary_bike_against_the_throttle_torque() {
    let mut sim = MotorcycleDynamicsSimulator::new();
    sim.set_brakes(200.0, 200.0);
    for _ in 0..100 {
        sim.step(DT);
    }
    assert_eq!(sim.state.speed, 0.0);
    assert_eq!(sim.state.rear_wheel_omega, 0.0);
    assert_eq!(sim.state.front_wheel_omega, 0.0);

    // A rear brake that is stronger than the drive torque keeps the wheel locked.
    sim.set_engine_torque(20.0);
    sim.set_brakes(200.0, 200.0);
    for _ in 0..100 {
        sim.step(DT);
    }
    assert_eq!(sim.state.speed, 0.0);
    assert_eq!(sim.state.rear_wheel_omega, 0.0);
}

// ---------------------------------------------------------------------------
// Yaw
// ---------------------------------------------------------------------------

/// A bike rolling at `speed` whose roll motion is frozen, so that only the yaw
/// dynamics are exercised.
fn frozen_lean_bike(speed: f64) -> MotorcycleDynamicsSimulator {
    let mut sim = rolling_bike(speed);
    sim.params.roll_inertia = 1.0e12;
    sim
}

#[test]
fn constant_steer_gives_the_bicycle_model_yaw_rate() {
    for steer in [0.01_f64, -0.01] {
        // coasting: a drive force on the rear tyre would eat into its lateral grip
        let mut sim = frozen_lean_bike(20.0);
        sim.state.steer_angle = steer;
        for _ in 0..400 {
            sim.step(0.005);
            assert_state_is_sane(&sim);
        }
        // Steady cornering of the linear single-track model with the lateral
        // velocity as a state: r = v * delta / (L + K_us * v^2), with the
        // understeer gradient K_us = m / L * (b / C_f - a / C_r) of the axle
        // loads the bike carries.  The drag of the coasting bike moves a little
        // load onto the front axle, which makes it slightly oversteering.
        let wheelbase = sim.params.wheelbase;
        let a = sim.params.a;
        let b = wheelbase - a;
        let tyre = &sim.front_tire.lat;
        let stiffness_per_load = tyre.b * tyre.c * tyre.d;
        let c_front = stiffness_per_load * sim.state.front_normal;
        let c_rear = stiffness_per_load * sim.state.rear_normal;
        let understeer_gradient = sim.params.total_mass / wheelbase * (b / c_front - a / c_rear);
        let v = sim.state.speed;
        let expected = v * steer / (wheelbase + understeer_gradient * v * v);
        let kinematic = v * steer / wheelbase;
        let yaw = sim.state.yaw_rate;
        assert!(
            (yaw - expected).abs() < 0.05 * expected.abs(),
            "steer {steer}: yaw rate {yaw}, bicycle model {expected}"
        );
        // and it stays close to the neutral-steer value v * delta / L
        assert!(
            (yaw - kinematic).abs() < 0.2 * kinematic.abs(),
            "steer {steer}: yaw rate {yaw}, kinematic {kinematic}"
        );
    }
}

#[test]
fn yaw_disturbance_decays() {
    // The tyre relaxation makes the response oscillate a little, but the
    // oscillation has to die out; with the lever arms of the axles swapped it
    // used to grow without bound.
    let mut sim = frozen_lean_bike(20.0);
    sim.state.yaw_rate = 0.3;
    let mut peak_after_a_tenth_of_a_second: f64 = 0.0;
    for step in 0..200 {
        sim.step(0.005);
        assert_state_is_sane(&sim);
        if step >= 20 {
            peak_after_a_tenth_of_a_second =
                peak_after_a_tenth_of_a_second.max(sim.state.yaw_rate.abs());
        }
    }
    assert!(
        peak_after_a_tenth_of_a_second < 0.05,
        "yaw rate still reaches {peak_after_a_tenth_of_a_second} rad/s"
    );
    assert!(
        sim.state.yaw_rate.abs() < 1e-3,
        "yaw rate after 1 s: {}",
        sim.state.yaw_rate
    );
}

#[test]
fn yaw_response_scales_with_the_steer_angle() {
    let response = |steer: f64| -> f64 {
        let mut sim = frozen_lean_bike(20.0);
        sim.state.steer_angle = steer;
        for _ in 0..200 {
            sim.step(0.005);
        }
        sim.state.yaw_rate
    };
    let small = response(0.005);
    let double = response(0.01);
    assert!(small > 0.0 && double > 0.0);
    assert!(
        (double / small - 2.0).abs() < 0.05,
        "doubling the steer angle gives {}x yaw rate",
        double / small
    );
}

// ---------------------------------------------------------------------------
// Robustness
// ---------------------------------------------------------------------------

/// Small deterministic xorshift generator for the randomised scenarios.
struct XorShift(u64);

impl XorShift {
    fn next_unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

#[test]
fn random_inputs_and_step_sizes_keep_the_state_physical() {
    let mut rng = XorShift(0xfeed_beef_1234_5678);
    for case in 0..24 {
        let dt = [0.005, 0.01, 0.02, 0.05][case % 4];
        let start_speed = rng.uniform(0.0, 60.0);
        let mut sim = rolling_bike(start_speed);
        // freeze the roll motion: the inverted-pendulum lean is not under test
        sim.params.roll_inertia = 1.0e12;
        let slip = rng.uniform(-0.3, 0.3);
        sim.state.rear_wheel_omega *= 1.0 + slip;
        sim.state.front_wheel_omega *= 1.0 + slip;
        let steps = (3.0 / dt) as usize;
        for step in 0..steps {
            if step.is_multiple_of(10) {
                let throttle = if rng.next_unit() < 0.5 {
                    rng.uniform(0.0, 150.0)
                } else {
                    0.0
                };
                sim.set_engine_torque(throttle);
                if rng.next_unit() < 0.5 {
                    sim.set_brakes(rng.uniform(0.0, 800.0), rng.uniform(0.0, 600.0));
                } else {
                    sim.set_brakes(0.0, 0.0);
                }
            }
            sim.step(dt);
            assert_state_is_sane(&sim);
        }
    }
}
