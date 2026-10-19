// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the gauge-pressure bookkeeping of [`CornerBrake`] and
//! for the stateful ABS pressure modulation of [`BrakeSystem`].

use super::*;

const DT: f64 = 0.001;
const SPEED: f64 = 20.0;
const RADIUS: f64 = 0.3;
const AMBIENT: f64 = 25.0;
/// Pedal force that drives the pedal to its end stop (maximum pressure).
const FULL_PEDAL: f64 = 120.0;

fn rolling_wheels() -> [f64; 4] {
    [SPEED / RADIUS; 4]
}

/// A system with a fully pressed pedal, settled with rolling (unslipped) wheels.
fn settled_full_pedal() -> (BrakeSystem, [f64; 4]) {
    let mut bs = BrakeSystem::new_default();
    for _ in 0..200 {
        bs.step(FULL_PEDAL, rolling_wheels(), SPEED, RADIUS, AMBIENT, DT);
    }
    let base = [
        bs.corners[0].gauge_pressure,
        bs.corners[1].gauge_pressure,
        bs.corners[2].gauge_pressure,
        bs.corners[3].gauge_pressure,
    ];
    (bs, base)
}

// ---------------------------------------------------------------------------
// F4779: gauge pressure
// ---------------------------------------------------------------------------

#[test]
fn new_corner_brake_is_released() {
    let caliper = BrakeCaliper::new(4, 0.042, PadMaterial::semi_metallic());
    let rotor = BrakeRotor::new_vented(0.16, 0.088, 0.028, 8.0);
    let mut corner = CornerBrake::new(caliper, rotor);
    assert_eq!(corner.gauge_pressure, 0.0);
    assert_eq!(corner.torque(), 0.0);
    assert_eq!(corner.braking_force(0.3), 0.0);
    // Releasing heat generation: no pressure, no heat, rotor stays at ambient.
    corner.update_thermal(80.0, AMBIENT, 1.0);
    assert_eq!(corner.rotor.temperature, AMBIENT);
    assert_eq!(corner.caliper.pad_temperature, AMBIENT);
}

#[test]
fn idle_brake_system_has_no_drag_torque_and_no_heating() {
    // Before the fix the corners stored ATM + p as the applied pressure, so
    // the seal-force offset made the front corners drag with ~55 N·m even with
    // the pedal released, heating the rotors.
    let mut bs = BrakeSystem::new_default();
    for _ in 0..10_000 {
        bs.step(0.0, [50.0; 4], 15.0, 0.3, AMBIENT, 0.01);
    }
    for (i, corner) in bs.corners.iter().enumerate() {
        assert_eq!(corner.gauge_pressure, 0.0, "corner {i}");
        assert_eq!(corner.torque(), 0.0, "corner {i} drags");
        assert!(
            (corner.rotor.temperature - AMBIENT).abs() < 1e-12,
            "corner {i} rotor temperature {}",
            corner.rotor.temperature
        );
        assert!((corner.caliper.pad_temperature - AMBIENT).abs() < 1e-12);
    }
    assert_eq!(bs.energy_dissipated_j, 0.0);
    assert_eq!(bs.total_braking_power([50.0; 4]), 0.0);
    assert!((bs.max_rotor_temperature() - AMBIENT).abs() < 1e-12);
}

#[test]
fn corner_torque_agrees_with_the_torque_returned_by_step() {
    let mut bs = BrakeSystem::new_default();
    let mut torques = [0.0; 4];
    for _ in 0..300 {
        torques = bs.step(FULL_PEDAL, rolling_wheels(), SPEED, RADIUS, AMBIENT, DT);
    }
    for (i, torque) in torques.iter().enumerate() {
        assert!(*torque > 50.0, "corner {i} must brake: {torque}");
        // `step` reports the torque before the rotor heated up during the step,
        // so the two agree up to the tiny friction change of one millisecond.
        assert!(
            (bs.corners[i].torque() - torque).abs() < 1e-4 * torque,
            "corner {i}: CornerBrake::torque {} vs step {torque}",
            bs.corners[i].torque()
        );
    }
    let omega = rolling_wheels();
    let expected_power: f64 = torques
        .iter()
        .zip(omega.iter())
        .map(|(t, w)| t * w.abs())
        .sum();
    assert!((bs.total_braking_power(omega) - expected_power).abs() < 1e-4 * expected_power);
}

#[test]
fn corner_thermal_state_is_driven_by_the_gauge_pressure() {
    let caliper = BrakeCaliper::new(4, 0.042, PadMaterial::semi_metallic());
    let rotor = BrakeRotor::new_vented(0.16, 0.088, 0.028, 8.0);
    let mut corner = CornerBrake::new(caliper, rotor);
    corner.gauge_pressure = 1.0e6;
    let torque = corner.torque();
    assert!(torque > 100.0);
    corner.update_thermal(60.0, AMBIENT, 0.5);
    // Heat input equals torque times wheel speed; 60 % goes into the rotor.
    assert!(corner.rotor.temperature > AMBIENT + 0.1);
    assert!(corner.caliper.pad_temperature > AMBIENT);
}

// ---------------------------------------------------------------------------
// F4780: stateful ABS modulation
// ---------------------------------------------------------------------------

#[test]
fn abs_releases_on_every_step_while_the_wheel_over_slips() {
    // The old controller released once per hold_time (20 ms): 4 kPa out of
    // several hundred kPa.
    let mut abs = AbsBrakeController::new();
    let mut total = 0.0;
    for k in 0..10 {
        let dp = abs.regulate_wheel(0, 600_000.0, 0.5, DT);
        assert!(
            (dp + abs.release_rate * DT).abs() < 1e-9,
            "step {k}: dp = {dp}"
        );
        total += dp;
    }
    assert!((total + 40_000.0).abs() < 1e-6);
    assert_eq!(abs.phase[0], AbsPhase::Release);
    assert!(abs.active[0] && !abs.active[1]);
}

#[test]
fn abs_release_never_takes_more_than_the_applied_pressure() {
    let mut abs = AbsBrakeController::new();
    assert_eq!(abs.regulate_wheel(0, 1_000.0, 0.9, DT), -1_000.0);
    assert_eq!(abs.regulate_wheel(0, 0.0, 0.9, DT), 0.0);
}

#[test]
fn abs_enforces_a_minimum_release_dwell_before_building() {
    let mut abs = AbsBrakeController::new();
    assert!(abs.regulate_wheel(1, 500_000.0, 0.9, DT) < 0.0);
    // The wheel recovers immediately (slip 0), but the release lasts hold_time.
    let mut release_steps_after_recovery = 0;
    loop {
        let dp = abs.regulate_wheel(1, 400_000.0, 0.0, DT);
        if dp >= 0.0 {
            assert!(dp > 0.0, "build expected once the dwell has elapsed");
            break;
        }
        release_steps_after_recovery += 1;
        assert!(release_steps_after_recovery < 100, "never left release");
    }
    assert!(
        (18..=21).contains(&release_steps_after_recovery),
        "release continued for {release_steps_after_recovery} ms after recovery"
    );
    assert_eq!(abs.phase[1], AbsPhase::Build);
}

#[test]
fn abs_switches_from_build_back_to_release_without_waiting() {
    let mut abs = AbsBrakeController::new();
    abs.regulate_wheel(2, 500_000.0, 0.9, DT); // start releasing
    abs.hold_time = 0.0; // allow the switch to build at once
    assert!(abs.regulate_wheel(2, 400_000.0, 0.0, DT) > 0.0);
    assert_eq!(abs.phase[2], AbsPhase::Build);
    abs.hold_time = 1.0; // a long dwell must not delay a new release
    assert!(abs.regulate_wheel(2, 400_000.0, 0.9, DT) < 0.0);
    assert_eq!(abs.phase[2], AbsPhase::Release);
}

#[test]
fn abs_holds_pressure_inside_the_slip_band() {
    let mut abs = AbsBrakeController::new();
    abs.regulate_wheel(3, 500_000.0, 0.9, DT);
    abs.hold_time = 0.0;
    for _ in 0..5 {
        assert_eq!(abs.regulate_wheel(3, 300_000.0, 0.15, DT), 0.0);
    }
    assert_eq!(abs.phase[3], AbsPhase::Hold);
}

#[test]
fn abs_modulate_passes_the_driver_demand_through_when_not_regulating() {
    let mut abs = AbsBrakeController::new();
    for slip in [0.0, 0.05, 0.12, 0.15, 0.19] {
        let applied = abs.modulate(0, 750_000.0, slip, DT);
        assert_eq!(applied, 750_000.0, "slip {slip}");
        assert!(!abs.active[0]);
    }
}

#[test]
fn abs_modulated_pressure_accumulates_and_can_rebuild_from_zero() {
    let mut abs = AbsBrakeController::new();
    abs.hold_time = 0.0;
    let demand = 600_000.0;
    let mut previous = demand;
    for _ in 0..200 {
        let applied = abs.modulate(0, demand, 1.0, DT);
        assert!(applied <= previous, "pressure must fall monotonically");
        previous = applied;
    }
    // 200 ms at 4 MPa/s would remove 800 kPa: fully dumped, never negative.
    assert_eq!(previous, 0.0);
    assert!(abs.active[0]);
    // Rebuilding from zero used to be impossible (the build step was capped
    // at 10 % of the current pressure).
    let mut applied = 0.0;
    for _ in 0..10 {
        applied = abs.modulate(0, demand, 0.0, DT);
    }
    assert!(
        (applied - 10.0 * abs.build_rate * DT).abs() < 1e-6,
        "rebuilt pressure {applied}"
    );
    // Once the demand is restored the intervention ends.
    for _ in 0..1_000 {
        applied = abs.modulate(0, demand, 0.0, DT);
    }
    assert_eq!(applied, demand);
    assert!(!abs.active[0]);
    assert_eq!(abs.phase[0], AbsPhase::Hold);
}

#[test]
fn abs_modulated_pressure_never_exceeds_a_reduced_driver_demand() {
    let mut abs = AbsBrakeController::new();
    abs.modulate(0, 600_000.0, 0.9, DT);
    abs.hold_time = 0.0;
    let applied = abs.modulate(0, 100_000.0, 0.0, DT);
    assert!(applied <= 100_000.0, "applied {applied}");
    // The pedal is released altogether.
    let applied = abs.modulate(0, 0.0, 0.0, DT);
    assert_eq!(applied, 0.0);
    assert!(!abs.active[0], "no demand left to modulate");
}

#[test]
fn abs_reset_clears_the_modulation_state() {
    let mut abs = AbsBrakeController::new();
    for _ in 0..30 {
        abs.modulate(1, 600_000.0, 0.9, DT);
    }
    assert!(abs.is_active());
    assert!(abs.modulated_pressure[1] < 600_000.0);
    abs.reset();
    assert!(!abs.is_active());
    assert_eq!(abs.phase, [AbsPhase::Hold; 4]);
    assert_eq!(abs.modulated_pressure, [0.0; 4]);
    assert_eq!(abs.cycle_timer, [0.0; 4]);
}

#[test]
fn locked_wheels_lose_most_of_their_brake_pressure_within_100_ms() {
    let (mut bs, base) = settled_full_pedal();
    assert!(!bs.abs.is_active(), "rolling wheels must not trigger ABS");
    assert!(base[0] > 600_000.0 && base[2] > 400_000.0, "base {base:?}");
    assert_eq!(base[0], base[1]);
    assert_eq!(base[2], base[3]);

    let locked = [0.0; 4];
    let mut previous = base;
    for step in 0..100 {
        bs.step(FULL_PEDAL, locked, SPEED, RADIUS, AMBIENT, DT);
        for (i, (corner, before)) in bs.corners.iter().zip(previous.iter_mut()).enumerate() {
            let p = corner.gauge_pressure;
            assert!(
                p <= *before + 1e-9,
                "corner {i} pressure rose at step {step}: {before} -> {p}"
            );
            *before = p;
        }
    }
    assert!(bs.abs.is_active());
    for (i, (corner, demand)) in bs.corners.iter().zip(base.iter()).enumerate() {
        let p = corner.gauge_pressure;
        assert!(
            p < 0.5 * demand,
            "corner {i}: {p} Pa is not below half of the demand {demand} Pa"
        );
        // The torque follows the modulated pressure, not the pedal demand.
        let demand_torque = corner.caliper.braking_torque(*demand, &corner.rotor);
        assert!(corner.torque() < 0.5 * demand_torque);
    }
}

#[test]
fn abs_pressure_recovers_to_the_driver_demand_after_the_wheels_spin_up() {
    let (mut bs, base) = settled_full_pedal();
    for _ in 0..100 {
        bs.step(FULL_PEDAL, [0.0; 4], SPEED, RADIUS, AMBIENT, DT);
    }
    assert!(bs.abs.is_active());
    let dumped = bs.corners[0].gauge_pressure;
    for _ in 0..400 {
        bs.step(FULL_PEDAL, rolling_wheels(), SPEED, RADIUS, AMBIENT, DT);
    }
    assert!(
        !bs.abs.is_active(),
        "ABS must end once the demand is restored"
    );
    for (i, (corner, demand)) in bs.corners.iter().zip(base.iter()).enumerate() {
        assert!(
            (corner.gauge_pressure - demand).abs() < 1e-6,
            "corner {i}: {} vs demand {demand}",
            corner.gauge_pressure
        );
    }
    assert!(bs.corners[0].gauge_pressure > dumped);
}

#[test]
fn abs_holds_the_pressure_while_the_slip_stays_in_the_target_band() {
    let (mut bs, _) = settled_full_pedal();
    for _ in 0..40 {
        bs.step(FULL_PEDAL, [0.0; 4], SPEED, RADIUS, AMBIENT, DT);
    }
    // Slip 0.15 is the middle of the 0.10..0.20 band.
    let in_band = [0.85 * SPEED / RADIUS; 4];
    for _ in 0..30 {
        bs.step(FULL_PEDAL, in_band, SPEED, RADIUS, AMBIENT, DT);
    }
    let held = bs.corners[0].gauge_pressure;
    for _ in 0..50 {
        bs.step(FULL_PEDAL, in_band, SPEED, RADIUS, AMBIENT, DT);
        assert_eq!(bs.corners[0].gauge_pressure, held);
    }
    assert!(bs.abs.is_active());
    assert_eq!(bs.abs.phase[0], AbsPhase::Hold);
}

#[test]
fn abs_resets_cleanly_through_the_brake_system() {
    let (mut bs, _) = settled_full_pedal();
    for _ in 0..50 {
        bs.step(FULL_PEDAL, [0.0; 4], SPEED, RADIUS, AMBIENT, DT);
    }
    assert!(bs.abs.is_active());
    bs.abs.reset();
    assert!(!bs.abs.is_active());
    bs.step(FULL_PEDAL, rolling_wheels(), SPEED, RADIUS, AMBIENT, DT);
    assert!(
        bs.corners[0].gauge_pressure > 600_000.0,
        "after a reset the driver demand applies again"
    );
}

// ---------------------------------------------------------------------------
// Closed loop: ABS on a slippery surface
// ---------------------------------------------------------------------------

/// Friction coefficient of the rational slip curve with its peak `mu_peak`
/// at slip `0.15`.
fn friction(slip: f64, mu_peak: f64) -> f64 {
    let s_peak = 0.15;
    mu_peak * 2.0 * s_peak * slip / (s_peak * s_peak + slip * slip)
}

/// Brake from 20 m/s on a slippery road; returns the final vehicle speed and
/// the lowest front wheel speed seen (rad/s).
fn emergency_stop(abs_enabled: bool) -> (f64, f64) {
    let mu_peak = 0.15;
    let wheel_inertia = 1.2;
    let vehicle_mass = 1500.0;
    let normal_load = [4_100.0, 4_100.0, 3_300.0, 3_300.0];
    let mut bs = BrakeSystem::new_default();
    if !abs_enabled {
        bs.abs.release_rate = 0.0;
    }
    let mut speed = SPEED;
    let mut omega = [SPEED / RADIUS; 4];
    let mut min_front_omega = omega[0];
    for _ in 0..1_500 {
        let torques = bs.step(FULL_PEDAL, omega, speed, RADIUS, AMBIENT, DT);
        let mut longitudinal = 0.0;
        for i in 0..4 {
            let slip = AbsBrakeController::wheel_slip(speed, omega[i] * RADIUS);
            let tyre_force = friction(slip, mu_peak) * normal_load[i];
            longitudinal += tyre_force;
            let accel = (tyre_force * RADIUS - torques[i]) / wheel_inertia;
            omega[i] = (omega[i] + accel * DT).max(0.0);
        }
        speed = (speed - longitudinal / vehicle_mass * DT).max(0.0);
        min_front_omega = min_front_omega.min(omega[0].min(omega[1]));
    }
    (speed, min_front_omega)
}

#[test]
fn abs_prevents_wheel_lock_and_shortens_the_stop_on_a_slippery_road() {
    let (speed_without, min_omega_without) = emergency_stop(false);
    let (speed_with, min_omega_with) = emergency_stop(true);
    // Without any pressure modulation the front wheels lock.
    assert!(
        min_omega_without < 1.0,
        "front wheels should lock without ABS: {min_omega_without}"
    );
    // With ABS they keep rolling at a useful fraction of the vehicle speed.
    assert!(
        min_omega_with > 15.0,
        "ABS must keep the front wheels turning: {min_omega_with}"
    );
    // Regulating around peak friction decelerates the car more than sliding.
    assert!(
        speed_with < speed_without - 0.2,
        "speed after 1.5 s: ABS {speed_with}, locked {speed_without}"
    );
}
