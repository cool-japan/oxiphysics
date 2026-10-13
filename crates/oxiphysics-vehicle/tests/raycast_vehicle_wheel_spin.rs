// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the signed wheel-spin integration of `RaycastVehicle`.
//!
//! Wheel spin used to be clamped to `>= 0`, which made reverse gear unable to
//! move the car, and wheels with index `>= 4` were never integrated, so on
//! 6/8-wheel vehicles they stayed locked and braked permanently.  The brake is
//! now a dissipative torque that can stop a wheel but never reverse it.

use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBody;
use oxiphysics_vehicle::drivetrain::DriveLayout;
use oxiphysics_vehicle::{
    LinearSuspension, LinearTire, RaycastVehicle, Wheel, WheelState, flat_ground_query,
};

/// Fine step of the settling phase and of the manoeuvres that pass through
/// standstill.  Below 1 m/s the slip ratio is normalised by 0.5 m/s, so the
/// tyre couples the chassis speed to the wheel like a damper of `k_x / 0.5`
/// Ns/m; the chassis side of that coupling is explicit and needs a step of a
/// few milliseconds at most (about 6 ms for a 1.5 t car on the default tyre).
/// 0.1 ms keeps it far inside that limit.  Tests of the game-rate steps start
/// from rolling speed instead.
const DT: f64 = 0.0001;
const GRAVITY: f64 = 9.81;
const WHEEL_RADIUS: f64 = 0.3;
/// Time the suspension needs to settle on its springs before a manoeuvre.
const SETTLE_SECONDS: f64 = 2.0;

struct Rig {
    vehicle: RaycastVehicle,
    chassis: RigidBody,
    /// Time step of [`Rig::step`]; [`DT`] unless a test switches it.
    dt: f64,
}

impl Rig {
    /// A 1500 kg chassis released just above its static ride height; gravity is
    /// on, so the springs carry the weight and the tyres see realistic loads.
    fn new(vehicle: RaycastVehicle) -> Self {
        let mut chassis = RigidBody::new(1500.0);
        chassis.transform.position = Vec3::new(0.0, 0.52, 0.0);
        chassis.linear_damping = 0.0;
        chassis.angular_damping = 0.0;
        Self {
            vehicle,
            chassis,
            dt: DT,
        }
    }

    fn rear_wheel_drive() -> Self {
        let mut vehicle = RaycastVehicle::new_default_4wheel();
        vehicle.drivetrain.layout = DriveLayout::RearWheelDrive;
        Self::new(vehicle)
    }

    /// Use the time step `dt` for every following [`Rig::step`].
    fn set_step(&mut self, dt: f64) {
        self.dt = dt;
    }

    /// Advance the vehicle and its chassis by `seconds` of simulated time.
    fn advance_seconds(&mut self, seconds: f64) {
        let steps = (seconds / self.dt).round() as usize;
        for _ in 0..steps {
            self.step();
        }
    }

    /// Let the car settle on its suspension with no driver input, then check
    /// that every wheel carries load and the chassis is at rest.  Every
    /// manoeuvre starts from here, so a test can never pass or fail because the
    /// car left the ground.
    fn settle(&mut self) {
        self.vehicle.set_inputs(0.0, 0.0, 0.0);
        self.advance_seconds(SETTLE_SECONDS);
        self.assert_grounded("after settling");
        assert!(
            self.chassis.velocity.norm() < 0.01,
            "the settled chassis must be at rest, v = {:?}",
            self.chassis.velocity
        );
        // The slow heave/pitch modes of the six-wheeler are still creeping by about
        // a micrometre per second after two seconds; a chattering wheel would be
        // off by whole rad/s.
        for wheel in 0..self.vehicle.num_wheels() {
            assert!(
                self.spin(wheel).abs() < 1e-4,
                "wheel {wheel} must not spin at rest: {}",
                self.spin(wheel)
            );
        }
    }

    /// Every wheel must be on the ground and loaded.
    fn assert_grounded(&self, context: &str) {
        assert!(
            self.vehicle.all_wheels_in_contact(),
            "{context}: a wheel left the ground (contacts: {:?}, chassis y = {})",
            self.vehicle
                .wheel_states
                .iter()
                .map(|w| w.is_in_contact)
                .collect::<Vec<_>>(),
            self.chassis.transform.position.y
        );
        for (i, state) in self.vehicle.wheel_states.iter().enumerate() {
            assert!(
                state.suspension_force > 0.0,
                "{context}: wheel {i} carries no load"
            );
        }
    }

    fn step(&mut self) {
        let gravity = Vec3::new(0.0, -GRAVITY, 0.0);
        let ground = flat_ground_query(0.0);
        let suspension = LinearSuspension;
        let tire = LinearTire::default();
        self.chassis.force_accumulator = Vec3::zeros();
        self.chassis.torque_accumulator = Vec3::zeros();
        let dt = self.dt;
        self.vehicle
            .update(dt, &mut self.chassis, &gravity, &ground, &suspension, &tire);
        self.chassis.integrate_forces(dt, &gravity);
        self.chassis.integrate_velocity(dt);
    }

    fn spin(&self, wheel: usize) -> f64 {
        self.vehicle.wheel_states[wheel].spin_velocity
    }
}

/// A default 4-wheel car with a third axle (wheels 4 and 5) behind the rear axle.
fn six_wheel_vehicle() -> RaycastVehicle {
    let mut vehicle = RaycastVehicle::new_default_4wheel();
    vehicle.drivetrain.layout = DriveLayout::RearWheelDrive;
    for x in [-0.75, 0.75] {
        vehicle.wheels.push(Wheel {
            connection_point: Vec3::new(x, -0.1, -2.5),
            ..Wheel::default()
        });
        vehicle.wheel_states.push(WheelState::default());
    }
    assert_eq!(vehicle.num_wheels(), 6);
    vehicle
}

#[test]
fn settled_rigs_carry_the_weight_on_every_wheel() {
    // Sanity check of the rig itself: with gravity on, the springs carry the
    // whole weight of the car, for the 4- and the 6-wheeler alike.
    for vehicle in [RaycastVehicle::new_default_4wheel(), six_wheel_vehicle()] {
        let wheels = vehicle.num_wheels();
        let mut rig = Rig::new(vehicle);
        rig.settle();
        let load: f64 = rig
            .vehicle
            .wheel_states
            .iter()
            .map(|state| state.suspension_force)
            .sum();
        let weight = 1500.0 * GRAVITY;
        assert!(
            (load - weight).abs() < 0.02 * weight,
            "{wheels} wheels carry {load} N, weight is {weight} N"
        );
    }
}

#[test]
fn reverse_gear_spins_the_driven_wheels_backwards_and_moves_the_car_backwards() {
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(-1);
    rig.vehicle.set_inputs(1.0, 0.0, 0.0);
    rig.advance_seconds(2.0);
    rig.assert_grounded("reversing");

    for wheel in [2, 3] {
        assert!(
            rig.spin(wheel) < -1.0,
            "driven wheel {wheel} must spin backwards, spin = {}",
            rig.spin(wheel)
        );
    }
    assert!(
        rig.chassis.velocity.z < -0.5,
        "reverse gear must move the chassis backwards, vz = {}",
        rig.chassis.velocity.z
    );
    assert!(
        rig.vehicle.state.speed < -0.5,
        "signed forward speed must be negative, speed = {}",
        rig.vehicle.state.speed
    );
}

#[test]
fn reverse_gear_is_the_mirror_image_of_first_gear_in_wheel_speed_direction() {
    let mut forward = Rig::rear_wheel_drive();
    forward.settle();
    forward.vehicle.drivetrain.gearbox.shift(1);
    forward.vehicle.set_inputs(1.0, 0.0, 0.0);
    forward.advance_seconds(1.0);

    let mut backward = Rig::rear_wheel_drive();
    backward.settle();
    backward.vehicle.drivetrain.gearbox.shift(-1);
    backward.vehicle.set_inputs(1.0, 0.0, 0.0);
    backward.advance_seconds(1.0);

    for wheel in [2, 3] {
        assert!(forward.spin(wheel) > 1.0, "{}", forward.spin(wheel));
        assert!(backward.spin(wheel) < -1.0, "{}", backward.spin(wheel));
    }
    // The undriven front wheels follow the ground in both directions.
    assert!(forward.spin(0) > 0.0 && backward.spin(0) < 0.0);
    assert!(forward.chassis.velocity.z > 0.0 && backward.chassis.velocity.z < 0.0);
}

#[test]
fn forward_gear_still_spins_forwards_and_accelerates() {
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(1);
    rig.vehicle.set_inputs(1.0, 0.0, 0.0);
    rig.advance_seconds(2.0);
    rig.assert_grounded("accelerating");

    assert!(rig.spin(2) > 1.0 && rig.spin(3) > 1.0);
    assert!(
        rig.chassis.velocity.z > 0.5,
        "vz = {}",
        rig.chassis.velocity.z
    );
}

#[test]
fn free_rolling_wheels_follow_a_backwards_moving_ground() {
    // Regression: with spin clamped to >= 0 the wheels stayed at rest while the
    // car rolled backwards, so a permanent slip ratio of +1 braked the car.
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(0);
    rig.chassis.velocity = Vec3::new(0.0, 0.0, -8.0);
    rig.vehicle.set_inputs(0.0, 0.0, 0.0);
    rig.advance_seconds(1.0);
    rig.assert_grounded("rolling backwards");

    let v = rig.chassis.velocity.z;
    assert!(v < -7.0, "a coasting car must not be braked: vz = {v}");
    for wheel in 0..4 {
        let expected = v / WHEEL_RADIUS;
        assert!(
            rig.spin(wheel) < 0.0 && (rig.spin(wheel) - expected).abs() < 0.1 * expected.abs(),
            "wheel {wheel}: spin {} should follow v / r = {expected}",
            rig.spin(wheel)
        );
    }
}

#[test]
fn braked_car_at_rest_does_not_creep_or_spin_the_wheels() {
    // The legacy drivetrain brake term is `-brake * signum(spin)`, and
    // `0.0_f64.signum()` is `+1.0`: a braked wheel at rest used to receive a
    // negative torque that the `>= 0` clamp then hid.  Without the clamp the
    // wheels must still stay put.
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(1);
    rig.vehicle.set_inputs(0.0, 1.0, 0.0);
    for _ in 0..(2.0 / DT) as usize {
        rig.step();
        for wheel in 0..4 {
            assert!(
                rig.spin(wheel).abs() < 1e-9,
                "braked wheel {wheel} must stay at rest, spin = {}",
                rig.spin(wheel)
            );
        }
    }
    assert!(
        rig.chassis.velocity.z.abs() < 1e-3,
        "the braked chassis must not creep: vz = {}",
        rig.chassis.velocity.z
    );
}

#[test]
fn brake_holds_the_wheels_against_engine_torque_in_both_directions() {
    // The brake (3000 N m per wheel) exceeds the drive torque available at idle
    // in first gear (about 580 N m per driven wheel), so the wheels stay locked.
    for gear in [1, -1] {
        let mut rig = Rig::rear_wheel_drive();
        rig.settle();
        rig.vehicle.drivetrain.gearbox.shift(gear);
        rig.vehicle.set_inputs(1.0, 1.0, 0.0);
        rig.advance_seconds(1.0);
        for wheel in 0..4 {
            assert!(
                rig.spin(wheel).abs() < 1e-9,
                "gear {gear}, wheel {wheel}: spin = {}",
                rig.spin(wheel)
            );
        }
        assert!(
            rig.chassis.velocity.z.abs() < 1e-3,
            "gear {gear}: vz = {}",
            rig.chassis.velocity.z
        );
    }
}

#[test]
fn light_brake_does_not_prevent_launch_when_drive_torque_exceeds_total_brake_torque() {
    // The brake acts on all four wheels but only the rear pair is driven.  At
    // idle the engine delivers about 1170 N m in total; 5 % brake is 150 N m on
    // each of the four wheels (600 N m in total), so the car must pull away.
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(1);
    rig.vehicle.set_inputs(1.0, 0.05, 0.0);
    rig.advance_seconds(2.0);
    rig.assert_grounded("launching against a light brake");
    assert!(rig.spin(2) > 1.0 && rig.spin(3) > 1.0);
    assert!(
        rig.chassis.velocity.z > 0.5,
        "vz = {}",
        rig.chassis.velocity.z
    );
}

#[test]
fn heavy_brake_holds_the_car_although_the_rear_wheels_alone_could_break_away() {
    // 15 % brake is 450 N m per wheel: below the 580 N m drive torque of one
    // rear wheel, but the locked front axle (7000 N of friction) holds the car,
    // and the total brake torque (1800 N m) exceeds the engine torque (1170 N m).
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(1);
    rig.vehicle.set_inputs(1.0, 0.15, 0.0);
    rig.advance_seconds(2.0);
    assert!(
        rig.chassis.velocity.z.abs() < 0.05,
        "vz = {}",
        rig.chassis.velocity.z
    );
    for wheel in [0, 1] {
        assert!(
            rig.spin(wheel).abs() < 1e-9,
            "front wheel {wheel} must stay locked"
        );
    }
}

#[test]
fn brake_stops_a_rolling_car_without_reversing_wheels_or_chassis() {
    let mut rig = Rig::rear_wheel_drive();
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(0);
    let speed = 15.0;
    rig.chassis.velocity = Vec3::new(0.0, 0.0, speed);
    for state in &mut rig.vehicle.wheel_states {
        state.spin_velocity = speed / WHEEL_RADIUS;
    }
    rig.vehicle.set_inputs(0.0, 1.0, 0.0);

    let mut previous_spin = [speed / WHEEL_RADIUS; 4];
    for _ in 0..(6.0 / DT) as usize {
        rig.step();
        for wheel in 0..4 {
            let spin = rig.spin(wheel);
            assert!(
                spin >= -1e-9,
                "the brake must never reverse wheel {wheel}: spin = {spin}"
            );
            // A forward-rolling wheel under braking only slows down.
            assert!(
                spin <= previous_spin[wheel] + 1e-9,
                "braking must not accelerate wheel {wheel}: {previous_spin:?} -> {spin}"
            );
            previous_spin[wheel] = spin;
        }
        assert!(
            rig.chassis.velocity.z > -0.25,
            "the chassis must not be driven backwards: vz = {}",
            rig.chassis.velocity.z
        );
    }
    assert!(
        rig.chassis.velocity.z.abs() < 0.25,
        "the car must come to a stop, vz = {}",
        rig.chassis.velocity.z
    );
    for wheel in 0..4 {
        assert!(rig.spin(wheel).abs() < 1e-9);
    }
}

#[test]
fn extra_axle_wheels_roll_with_the_ground_instead_of_locking() {
    // Regression: wheels with index >= 4 were never integrated, kept zero spin
    // and therefore produced a permanent slip ratio of -1 (locked-wheel braking).
    let mut rig = Rig::new(six_wheel_vehicle());
    rig.settle();
    let speed = 20.0;
    rig.chassis.velocity = Vec3::new(0.0, 0.0, speed);
    rig.vehicle.drivetrain.gearbox.shift(0);
    rig.vehicle.set_inputs(0.0, 0.0, 0.0);
    rig.advance_seconds(1.5);
    rig.assert_grounded("coasting on six wheels");

    // Spinning up six wheels from rest costs about 0.6 m/s and drag another 0.2;
    // two permanently locked wheels would have cost about 3 m/s.
    let v = rig.chassis.velocity.z;
    assert!(
        v > 18.5,
        "coasting must not be braked by locked wheels: vz = {v}"
    );
    let expected = v / WHEEL_RADIUS;
    for wheel in 0..6 {
        assert!(
            (rig.spin(wheel) - expected).abs() < 0.05 * expected,
            "wheel {wheel}: spin {} should be close to v / r = {expected}",
            rig.spin(wheel)
        );
    }
    for wheel in 4..6 {
        assert!(
            rig.vehicle.wheel_states[wheel].slip_ratio.abs() < 0.05,
            "extra-axle wheel {wheel} must roll without slip, slip ratio = {}",
            rig.vehicle.wheel_states[wheel].slip_ratio
        );
    }
}

#[test]
fn extra_axle_wheels_are_braked_and_not_driven() {
    let mut rig = Rig::new(six_wheel_vehicle());
    rig.settle();
    rig.vehicle.drivetrain.gearbox.shift(1);
    // Full throttle: only the legacy drivetrain wheels receive drive torque, the
    // extra axle follows the ground.
    rig.vehicle.set_inputs(1.0, 0.0, 0.0);
    rig.advance_seconds(2.0);
    rig.assert_grounded("accelerating on six wheels");
    let v = rig.chassis.velocity.z;
    assert!(v > 0.5, "the six-wheeler must accelerate, vz = {v}");
    for wheel in 4..6 {
        assert!(
            rig.spin(wheel) > 0.0 && rig.spin(wheel) < 1.2 * v / WHEEL_RADIUS,
            "undriven wheel {wheel} must roll at ground speed, spin = {} (v/r = {})",
            rig.spin(wheel),
            v / WHEEL_RADIUS
        );
    }

    // Now brake hard: the extra axle must be braked too.
    rig.vehicle.set_inputs(0.0, 1.0, 0.0);
    rig.advance_seconds(4.0);
    for wheel in 0..6 {
        assert!(
            rig.spin(wheel).abs() < 1e-9,
            "wheel {wheel} must be stopped by the brake, spin = {}",
            rig.spin(wheel)
        );
    }
    assert!(rig.chassis.velocity.z.abs() < 0.25);
}

/// Time steps of typical game loops: 60, 120 and 240 Hz.
const GAME_STEPS: [f64; 3] = [1.0 / 60.0, 1.0 / 120.0, 1.0 / 240.0];

/// Run `rig` at the step `dt` for `seconds`, calling `check` after every step
/// once `skip_seconds` of the manoeuvre have passed.
fn run_checked(rig: &mut Rig, dt: f64, seconds: f64, skip_seconds: f64, check: impl Fn(&Rig)) {
    rig.set_step(dt);
    let steps = (seconds / dt).round() as usize;
    for step in 0..steps {
        rig.step();
        if step as f64 * dt >= skip_seconds {
            check(rig);
        }
    }
}

#[test]
fn free_rolling_wheels_follow_the_ground_at_game_time_steps() {
    // The tyre/wheel coupling is stiff (about 800 / s at 20 m/s): explicit Euler
    // would flip the wheel spin by tens of rad/s every frame at these steps and
    // the tyre force would alternate between +mu*Fz and -mu*Fz.  The wheel must
    // stay on the rolling condition at every step.
    for dt in GAME_STEPS {
        let mut rig = Rig::rear_wheel_drive();
        rig.settle();
        rig.vehicle.drivetrain.gearbox.shift(0);
        let speed = 20.0;
        rig.chassis.velocity = Vec3::new(0.0, 0.0, speed);
        for state in &mut rig.vehicle.wheel_states {
            state.spin_velocity = speed / WHEEL_RADIUS;
        }
        rig.vehicle.set_inputs(0.0, 0.0, 0.0);
        run_checked(&mut rig, dt, 3.0, 0.5, |rig| {
            let v = rig.chassis.velocity.z;
            let expected = v / WHEEL_RADIUS;
            for wheel in 0..4 {
                assert!(
                    (rig.spin(wheel) - expected).abs() < 0.03 * expected,
                    "dt = {dt}: wheel {wheel} spin {} must follow v / r = {expected}",
                    rig.spin(wheel)
                );
            }
        });
        rig.assert_grounded("coasting at a game-rate step");
        assert!(
            rig.chassis.velocity.z > 19.5,
            "dt = {dt}: a coasting car must not be braked, vz = {}",
            rig.chassis.velocity.z
        );
    }
}

#[test]
fn extra_axle_wheels_follow_the_ground_at_game_time_steps() {
    // The 20 m/s coasting scenario of the six-wheel regression, at game steps.
    for dt in GAME_STEPS {
        let mut rig = Rig::new(six_wheel_vehicle());
        rig.settle();
        rig.vehicle.drivetrain.gearbox.shift(0);
        let speed = 20.0;
        rig.chassis.velocity = Vec3::new(0.0, 0.0, speed);
        rig.vehicle.set_inputs(0.0, 0.0, 0.0);
        // The wheels start at rest and must spin up to the ground speed.
        run_checked(&mut rig, dt, 3.0, 2.0, |rig| {
            let expected = rig.chassis.velocity.z / WHEEL_RADIUS;
            for wheel in 0..6 {
                assert!(
                    (rig.spin(wheel) - expected).abs() < 0.03 * expected,
                    "dt = {dt}: wheel {wheel} spin {} must follow v / r = {expected}",
                    rig.spin(wheel)
                );
            }
        });
        rig.assert_grounded("coasting on six wheels at a game-rate step");
        for wheel in 4..6 {
            assert!(rig.vehicle.wheel_states[wheel].slip_ratio.abs() < 0.03);
        }
    }
}

#[test]
fn braking_at_game_time_steps_stops_the_wheels_without_reversing_them() {
    for dt in GAME_STEPS {
        let mut rig = Rig::rear_wheel_drive();
        rig.settle();
        rig.vehicle.drivetrain.gearbox.shift(0);
        let speed = 15.0;
        rig.chassis.velocity = Vec3::new(0.0, 0.0, speed);
        for state in &mut rig.vehicle.wheel_states {
            state.spin_velocity = speed / WHEEL_RADIUS;
        }
        rig.vehicle.set_inputs(0.0, 1.0, 0.0);
        // One second of full braking: the car slows from 15 m/s to about 5 m/s and
        // is still rolling, so this stays clear of the standstill range.
        run_checked(&mut rig, dt, 1.0, 0.0, |rig| {
            for wheel in 0..4 {
                assert!(
                    rig.spin(wheel) >= 0.0,
                    "dt = {dt}: the brake reversed wheel {wheel}: {}",
                    rig.spin(wheel)
                );
            }
        });
        let v = rig.chassis.velocity.z;
        assert!(
            (4.5..6.5).contains(&v),
            "dt = {dt}: full braking for 1 s from 15 m/s must leave about 5 m/s, vz = {v}"
        );
        for wheel in 0..4 {
            assert!(
                rig.spin(wheel).abs() < 1e-9,
                "dt = {dt}: wheel {wheel} must be locked by the full brake, spin = {}",
                rig.spin(wheel)
            );
        }
    }
}
