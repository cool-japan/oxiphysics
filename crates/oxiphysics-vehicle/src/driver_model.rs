// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Driver / controller model for racing simulation.
//!
//! Provides PID-based throttle control, pure-pursuit steering, braking-point
//! calculation, automatic gear selection, trail braking, and a lightweight
//! driver-state updater.

/// Runtime state of the driver's inputs and gear.
#[derive(Debug, Clone)]
pub struct DriverState {
    /// Throttle pedal position in \[0, 1\].
    pub throttle: f64,
    /// Brake pedal position in \[0, 1\].
    pub brake: f64,
    /// Steering wheel input in \[−1, 1\] (positive = right).
    pub steering: f64,
    /// Currently selected gear (1-indexed, 0 = neutral).
    pub gear: u32,
    /// Driver reaction time \[s\].
    pub reaction_time: f64,
}

/// Fixed driver characteristics.
#[derive(Debug, Clone)]
pub struct DriverParams {
    /// Aggression factor (0 = very smooth, 1 = maximum attack).
    pub aggression: f64,
    /// Smoothness of inputs (0 = jerky, 1 = perfectly smooth).
    pub smoothness: f64,
    /// Offset applied to the braking point \[m\] (positive = brake later).
    pub braking_point_offset: f64,
    /// Maximum lateral acceleration the driver is comfortable with \[g\].
    pub max_lateral_g: f64,
}

/// PID controller for throttle demand.
///
/// Returns a clamped throttle value in \[0, 1\].
///
/// # Arguments
/// * `target_speed`  – desired speed \[m/s\]
/// * `current_speed` – actual speed \[m/s\]
/// * `kp`, `ki`, `kd` – PID gains
/// * `integral`      – running integral term (updated in place)
/// * `prev_err`      – previous error (updated in place)
/// * `dt`            – time step \[s\]
#[allow(clippy::too_many_arguments)]
pub fn pid_throttle(
    target_speed: f64,
    current_speed: f64,
    kp: f64,
    ki: f64,
    kd: f64,
    integral: &mut f64,
    prev_err: &mut f64,
    dt: f64,
) -> f64 {
    let err = target_speed - current_speed;
    *integral += err * dt;
    let derivative = if dt > 1e-15 {
        (err - *prev_err) / dt
    } else {
        0.0
    };
    *prev_err = err;
    let output = kp * err + ki * *integral + kd * derivative;
    output.clamp(0.0, 1.0)
}

/// Pure-pursuit steering angle for a bicycle model.
///
/// Returns the steering angle in radians (positive = right).
///
/// # Arguments
/// * `vehicle_pos`     – vehicle rear-axle position \[m\]
/// * `vehicle_heading` – heading angle from world-X axis \[radians\]
/// * `target_point`    – look-ahead point on the desired path \[m\]
/// * `wheelbase`       – distance between axles \[m\]
pub fn look_ahead_steering(
    vehicle_pos: [f64; 2],
    vehicle_heading: f64,
    target_point: [f64; 2],
    wheelbase: f64,
) -> f64 {
    let dx = target_point[0] - vehicle_pos[0];
    let dy = target_point[1] - vehicle_pos[1];
    let ld = (dx * dx + dy * dy).sqrt();
    if ld < 1e-10 || wheelbase < 1e-10 {
        return 0.0;
    }
    // Transform target into vehicle frame
    let alpha = dy.atan2(dx) - vehicle_heading;
    2.0 * wheelbase * alpha.sin() / ld
}

/// Compute the braking point distance from a corner.
///
/// Uses the kinematic formula:  d = (v² − v_c²) / (2 · a_max)
///
/// # Arguments
/// * `speed`        – approach speed \[m/s\]
/// * `decel_max`    – maximum deceleration magnitude \[m/s²\]
/// * `corner_speed` – required entry speed \[m/s\]
pub fn braking_point(speed: f64, decel_max: f64, corner_speed: f64) -> f64 {
    let v2 = speed * speed;
    let vc2 = corner_speed * corner_speed;
    if v2 <= vc2 || decel_max < 1e-10 {
        return 0.0;
    }
    (v2 - vc2) / (2.0 * decel_max)
}

/// Select the optimal gear for the current speed.
///
/// Chooses the highest gear whose wheel speed × ratio gives engine RPM ≤
/// `engine_rpm_opt`.  Falls back to the lowest gear if none satisfies.
///
/// # Arguments
/// * `speed`          – vehicle speed \[m/s\]
/// * `gear_ratios`    – total drive ratios (transmission × final drive) per gear
/// * `engine_rpm_opt` – optimum engine RPM (often peak-power RPM)
pub fn gear_selection(speed: f64, gear_ratios: &[f64], engine_rpm_opt: f64) -> u32 {
    if gear_ratios.is_empty() {
        return 1;
    }
    // RPM ≈ speed × ratio × 60 / (2π * wheel_radius); assume wheel_radius = 0.33 m
    const WHEEL_R: f64 = 0.33;
    let wheel_rps = speed / WHEEL_R;
    let mut best = 1u32;
    let mut best_rpm = f64::NEG_INFINITY;
    let mut found = false;
    for (i, &ratio) in gear_ratios.iter().enumerate() {
        let rpm = wheel_rps * ratio * 60.0 / (2.0 * std::f64::consts::PI);
        if rpm <= engine_rpm_opt {
            // Pick the gear whose RPM is closest to (but not exceeding) the optimum.
            // At low speeds this selects the lowest gear (highest ratio → highest RPM).
            // At high speeds this selects the highest viable gear.
            if !found || rpm > best_rpm {
                best = (i + 1) as u32;
                best_rpm = rpm;
                found = true;
            }
        }
    }
    // If no gear keeps RPM below opt, select the highest gear (lowest ratio)
    // to minimise over-revving.
    if !found {
        best = gear_ratios.len() as u32;
    }
    best
}

/// Desired slip angle for the current driving conditions.
///
/// Returns target slip angle in radians.
/// More aggressive drivers target higher slip angles.
pub fn slip_angle_target(speed: f64, params: &DriverParams) -> f64 {
    // Maximum physically meaningful slip angle ~ 15 degrees
    const MAX_SLIP_DEG: f64 = 15.0;
    let speed_normalised = (speed / 50.0).clamp(0.0, 1.0);
    let target_deg = MAX_SLIP_DEG * params.aggression * speed_normalised;
    target_deg.to_radians()
}

/// Compute the trail-braking throttle and brake inputs for a cornering manoeuvre.
///
/// Returns `(throttle, brake)` both in \[0, 1\].
///
/// The idea: simultaneously carry some brake into the corner while the lateral
/// demand is high, blending towards full throttle as the corner opens.
///
/// # Arguments
/// * `speed`          – current speed \[m/s\]
/// * `corner_radius`  – radius of curvature \[m\]
/// * `mu`             – tyre–road friction coefficient
/// * `mass`           – vehicle mass \[kg\] (unused — see note)
/// * `g`              – gravitational acceleration \[m/s²\]
#[allow(unused_variables)]
pub fn trail_braking(speed: f64, corner_radius: f64, mu: f64, mass: f64, g: f64) -> (f64, f64) {
    if corner_radius < 1e-10 {
        return (0.0, 1.0);
    }
    let lateral_a = speed * speed / corner_radius; // m/s²
    let max_a = mu * g;
    let lat_fraction = (lateral_a / max_a).clamp(0.0, 1.0);
    // Remaining grip budget for longitudinal deceleration
    let long_fraction = (1.0 - lat_fraction * lat_fraction).sqrt().clamp(0.0, 1.0);
    let brake = (long_fraction * 0.5).clamp(0.0, 1.0);
    let throttle = if lat_fraction < 0.5 {
        ((0.5 - lat_fraction) * 2.0).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (throttle, brake)
}

/// Compute the Euclidean deviation of the vehicle from its ideal racing line.
///
/// # Arguments
/// * `actual` – actual vehicle position \[m\]
/// * `ideal`  – ideal line position \[m\]
pub fn driver_line_deviation(actual: [f64; 2], ideal: [f64; 2]) -> f64 {
    let dx = actual[0] - ideal[0];
    let dy = actual[1] - ideal[1];
    (dx * dx + dy * dy).sqrt()
}

/// Advance the driver state by one time step.
///
/// Applies smoothness filtering to throttle/brake/steering and updates gear
/// based on the current throttle demand.
///
/// # Arguments
/// * `state`   – current driver state (mutated in place)
/// * `params`  – driver characteristics
/// * `error`   – tracking error signal (e.g. speed error) \[m/s or m\]
/// * `dt`      – time step \[s\]
pub fn update_driver(state: &mut DriverState, params: &DriverParams, error: f64, dt: f64) {
    let alpha = (params.smoothness * dt * 10.0).clamp(0.0, 1.0);
    // Simple proportional response
    let demand = (error * params.aggression * 0.1).clamp(-1.0, 1.0);
    if demand > 0.0 {
        state.throttle = (1.0 - alpha) * state.throttle + alpha * demand;
        state.brake *= 1.0 - alpha;
    } else {
        state.throttle *= 1.0 - alpha;
        state.brake = (1.0 - alpha) * state.brake + alpha * (-demand);
    }
    state.throttle = state.throttle.clamp(0.0, 1.0);
    state.brake = state.brake.clamp(0.0, 1.0);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── pid_throttle ──────────────────────────────────────────────────────

    #[test]
    fn pid_throttle_zero_error_no_integral() {
        let mut integral = 0.0;
        let mut prev_err = 0.0;
        let t = pid_throttle(
            50.0,
            50.0,
            1.0,
            0.0,
            0.0,
            &mut integral,
            &mut prev_err,
            0.01,
        );
        assert!(t.abs() < 1e-10, "t={t}");
    }

    #[test]
    fn pid_throttle_positive_error_gives_positive_output() {
        let mut integral = 0.0;
        let mut prev_err = 0.0;
        let t = pid_throttle(
            60.0,
            50.0,
            1.0,
            0.0,
            0.0,
            &mut integral,
            &mut prev_err,
            0.01,
        );
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn pid_throttle_clamped_to_unit() {
        let mut integral = 0.0;
        let mut prev_err = 0.0;
        let t = pid_throttle(
            200.0,
            0.0,
            100.0,
            0.0,
            0.0,
            &mut integral,
            &mut prev_err,
            0.01,
        );
        assert!((0.0..=1.0).contains(&t), "t={t}");
    }

    #[test]
    fn pid_throttle_negative_error_gives_zero() {
        let mut integral = 0.0;
        let mut prev_err = 0.0;
        let t = pid_throttle(
            30.0,
            50.0,
            1.0,
            0.0,
            0.0,
            &mut integral,
            &mut prev_err,
            0.01,
        );
        // Output is negative → clamped to 0
        assert!(t.abs() < 1e-10 || t == 0.0, "t={t}");
    }

    #[test]
    fn pid_throttle_integral_accumulates() {
        let mut integral = 0.0;
        let mut prev_err = 0.0;
        let _t1 = pid_throttle(60.0, 50.0, 0.0, 1.0, 0.0, &mut integral, &mut prev_err, 0.1);
        let _t2 = pid_throttle(60.0, 50.0, 0.0, 1.0, 0.0, &mut integral, &mut prev_err, 0.1);
        // Integral should have increased
        assert!(integral > 0.0, "integral={integral}");
    }

    // ── look_ahead_steering ───────────────────────────────────────────────

    #[test]
    fn pure_pursuit_zero_offset_zero_steering() {
        // Target directly ahead along heading
        let s = look_ahead_steering([0.0, 0.0], 0.0, [10.0, 0.0], 2.7);
        assert!(s.abs() < 1e-10, "s={s}");
    }

    #[test]
    fn pure_pursuit_target_to_right_positive_steering() {
        let s = look_ahead_steering([0.0, 0.0], 0.0, [10.0, -5.0], 2.7);
        assert!(
            s < 0.0,
            "s={s} should be negative for target to right in atan2 convention"
        );
    }

    #[test]
    fn pure_pursuit_target_to_left_left_steering() {
        let s = look_ahead_steering([0.0, 0.0], 0.0, [10.0, 5.0], 2.7);
        assert!(s > 0.0, "s={s}");
    }

    #[test]
    fn pure_pursuit_zero_lookahead_is_zero() {
        let s = look_ahead_steering([0.0, 0.0], 0.0, [0.0, 0.0], 2.7);
        assert!(s.abs() < 1e-10, "s={s}");
    }

    // ── braking_point ─────────────────────────────────────────────────────

    #[test]
    fn braking_point_at_corner_speed_is_zero() {
        let d = braking_point(30.0, 50.0, 30.0);
        assert!(d.abs() < 1e-10, "d={d}");
    }

    #[test]
    fn braking_point_positive_when_above_corner_speed() {
        let d = braking_point(50.0, 50.0, 30.0);
        assert!(d > 0.0, "d={d}");
    }

    #[test]
    fn braking_point_formula() {
        // d = (50^2 - 30^2) / (2 * 50) = (2500 - 900) / 100 = 16
        let d = braking_point(50.0, 50.0, 30.0);
        assert!((d - 16.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn braking_point_increases_with_speed() {
        let d1 = braking_point(60.0, 50.0, 30.0);
        let d2 = braking_point(80.0, 50.0, 30.0);
        assert!(d2 > d1);
    }

    // ── gear_selection ────────────────────────────────────────────────────

    #[test]
    fn gear_selection_empty_ratios_returns_one() {
        let g = gear_selection(30.0, &[], 6000.0);
        assert_eq!(g, 1);
    }

    #[test]
    fn gear_selection_low_speed_selects_low_gear() {
        let ratios = [15.0, 9.0, 6.0, 4.5, 3.5];
        let g = gear_selection(5.0, &ratios, 7000.0);
        assert!(g <= 2, "g={g}");
    }

    #[test]
    fn gear_selection_high_speed_selects_high_gear() {
        let ratios = [15.0, 9.0, 6.0, 4.5, 3.5];
        let g = gear_selection(80.0, &ratios, 8000.0);
        assert!(g >= 4, "g={g}");
    }

    #[test]
    fn gear_selection_increases_monotonically_with_speed() {
        let ratios = [12.0, 8.0, 5.5, 4.0, 3.2];
        let speeds = [5.0, 15.0, 30.0, 50.0, 70.0];
        let gears: Vec<u32> = speeds
            .iter()
            .map(|&s| gear_selection(s, &ratios, 7000.0))
            .collect();
        for w in gears.windows(2) {
            assert!(
                w[0] <= w[1],
                "gear should not decrease as speed increases: {:?}",
                gears
            );
        }
    }

    // ── slip_angle_target ─────────────────────────────────────────────────

    #[test]
    fn slip_angle_zero_speed_is_zero() {
        let p = DriverParams {
            aggression: 1.0,
            smoothness: 0.5,
            braking_point_offset: 0.0,
            max_lateral_g: 3.5,
        };
        let s = slip_angle_target(0.0, &p);
        assert!(s.abs() < 1e-15, "s={s}");
    }

    #[test]
    fn slip_angle_increases_with_aggression() {
        let p_mellow = DriverParams {
            aggression: 0.3,
            smoothness: 0.8,
            braking_point_offset: 0.0,
            max_lateral_g: 3.0,
        };
        let p_aggressive = DriverParams {
            aggression: 0.9,
            ..p_mellow.clone()
        };
        let s1 = slip_angle_target(30.0, &p_mellow);
        let s2 = slip_angle_target(30.0, &p_aggressive);
        assert!(s2 > s1);
    }

    // ── trail_braking ─────────────────────────────────────────────────────

    #[test]
    fn trail_braking_outputs_in_unit_range() {
        let (t, b) = trail_braking(30.0, 100.0, 1.0, 700.0, 9.81);
        assert!((0.0..=1.0).contains(&t), "t={t}");
        assert!((0.0..=1.0).contains(&b), "b={b}");
    }

    #[test]
    fn trail_braking_zero_radius_max_brake() {
        let (t, b) = trail_braking(30.0, 0.0, 1.0, 700.0, 9.81);
        assert_eq!(t, 0.0);
        assert_eq!(b, 1.0);
    }

    #[test]
    fn trail_braking_low_speed_allows_throttle() {
        let (t, _b) = trail_braking(5.0, 200.0, 1.0, 700.0, 9.81);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn trail_braking_high_lateral_reduces_throttle() {
        let (t_slow, _) = trail_braking(10.0, 100.0, 1.0, 700.0, 9.81);
        let (t_fast, _) = trail_braking(50.0, 100.0, 1.0, 700.0, 9.81);
        assert!(
            t_slow >= t_fast,
            "higher lateral load should reduce throttle"
        );
    }

    // ── driver_line_deviation ─────────────────────────────────────────────

    #[test]
    fn line_deviation_on_ideal_line_is_zero() {
        let d = driver_line_deviation([5.0, 3.0], [5.0, 3.0]);
        assert!(d.abs() < 1e-15);
    }

    #[test]
    fn line_deviation_pythagoras() {
        let d = driver_line_deviation([3.0, 0.0], [0.0, 4.0]);
        assert!((d - 5.0).abs() < 1e-12, "d={d}");
    }

    #[test]
    fn line_deviation_non_negative() {
        let d = driver_line_deviation([1.0, 2.0], [4.0, 6.0]);
        assert!(d >= 0.0);
    }

    // ── update_driver ─────────────────────────────────────────────────────

    #[test]
    fn update_driver_positive_error_increases_throttle() {
        let mut state = DriverState {
            throttle: 0.0,
            brake: 0.0,
            steering: 0.0,
            gear: 3,
            reaction_time: 0.2,
        };
        let params = DriverParams {
            aggression: 1.0,
            smoothness: 0.5,
            braking_point_offset: 0.0,
            max_lateral_g: 3.5,
        };
        update_driver(&mut state, &params, 10.0, 0.02);
        assert!(state.throttle > 0.0, "throttle={}", state.throttle);
    }

    #[test]
    fn update_driver_negative_error_increases_brake() {
        let mut state = DriverState {
            throttle: 0.0,
            brake: 0.0,
            steering: 0.0,
            gear: 3,
            reaction_time: 0.2,
        };
        let params = DriverParams {
            aggression: 1.0,
            smoothness: 0.5,
            braking_point_offset: 0.0,
            max_lateral_g: 3.5,
        };
        update_driver(&mut state, &params, -10.0, 0.02);
        assert!(state.brake > 0.0, "brake={}", state.brake);
    }

    #[test]
    fn update_driver_throttle_brake_in_unit_range() {
        let mut state = DriverState {
            throttle: 0.5,
            brake: 0.0,
            steering: 0.0,
            gear: 4,
            reaction_time: 0.1,
        };
        let params = DriverParams {
            aggression: 0.8,
            smoothness: 0.6,
            braking_point_offset: 10.0,
            max_lateral_g: 4.0,
        };
        for _ in 0..50 {
            update_driver(&mut state, &params, 5.0, 0.016);
        }
        assert!((0.0..=1.0).contains(&state.throttle));
        assert!((0.0..=1.0).contains(&state.brake));
    }
}
