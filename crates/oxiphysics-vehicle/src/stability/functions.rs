//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    AbsConfig, AbsState, ArcConfig, ArcOutput, EbdConfig, EscConfig, EscIntervention, EscState,
    EspLevel, HandlingBalance, TcsConfig, TcsState, TorqueVectoringOutput,
};

/// Compute the longitudinal slip ratio for a driven wheel.
///
/// `slip = (v_wheel - v_vehicle) / v_vehicle`, clamped to `[-1, 1]`.
///
/// Returns `0.0` when vehicle speed is near zero to avoid division by zero.
pub fn compute_slip_ratio(wheel_speed: f64, vehicle_speed: f64) -> f64 {
    if vehicle_speed.abs() < 1e-6 {
        return 0.0;
    }
    ((wheel_speed - vehicle_speed) / vehicle_speed).clamp(-1.0, 1.0)
}
/// Update TCS state and return the torque scaling factor `[0, 1]`.
///
/// When slip exceeds `config.slip_threshold`, the factor is multiplied by
/// `config.torque_reduction_rate` (floored at `min_torque_fraction`).
/// When slip is within threshold, the factor recovers toward 1.0.
pub fn tcs_update(state: &mut TcsState, slip: f64, config: &TcsConfig) -> f64 {
    state.slip_ratio = slip;
    if slip.abs() > config.slip_threshold {
        state.is_active = true;
        state.torque_reduction_factor = (state.torque_reduction_factor
            * config.torque_reduction_rate)
            .max(config.min_torque_fraction);
    } else {
        state.is_active = false;
        state.torque_reduction_factor = (state.torque_reduction_factor + 0.05).min(1.0);
    }
    state.torque_reduction_factor
}
/// Compute target yaw rate from the bicycle model (linear).
///
/// r = v / (L + Kus · v²) · δ
///
/// where L = `wheelbase`, Kus = `understeer_gradient`, δ = `steer_angle_deg` in radians.
pub fn target_yaw_rate(
    speed: f64,
    steer_angle_deg: f64,
    wheelbase: f64,
    understeer_gradient: f64,
) -> f64 {
    let delta = steer_angle_deg.to_radians();
    let denom = wheelbase + understeer_gradient * speed * speed;
    if denom.abs() < 1e-10 {
        return 0.0;
    }
    (speed / denom) * delta
}
/// Update ESC state and compute per-wheel intervention.
pub fn esc_update(
    state: &mut EscState,
    actual_yaw: f64,
    target_yaw: f64,
    sideslip_deg: f64,
    config: &EscConfig,
) -> EscIntervention {
    let yaw_error = actual_yaw - target_yaw;
    state.yaw_rate_error = yaw_error;
    let yaw_active = yaw_error.abs() > config.yaw_rate_threshold;
    let sideslip_active = sideslip_deg.abs() > config.sideslip_threshold;
    if !yaw_active && !sideslip_active {
        state.is_active = false;
        state.intervention_type = String::from("none");
        return EscIntervention::default();
    }
    state.is_active = true;
    let mut intervention = EscIntervention::default();
    if yaw_error > 0.0 {
        state.intervention_type = String::from("oversteer");
        let strength = (yaw_error.abs() * config.oversteer_gain).min(1.0);
        intervention.brake_front_left = strength;
        intervention.brake_rear_left = strength * 0.5;
        intervention.throttle_reduction = strength;
    } else {
        state.intervention_type = String::from("understeer");
        let strength = (yaw_error.abs() * config.understeer_gain).min(1.0);
        intervention.brake_front_right = strength;
        intervention.brake_rear_right = strength * 0.5;
        intervention.throttle_reduction = strength * 0.5;
    }
    intervention
}
/// Update ABS state for all four wheels and return modulated brake pressures.
///
/// For each wheel: if `|slip_ratio| > slip_target`, reduce pressure by
/// `pressure_decrease_rate * dt`; otherwise increase toward 1.0 by
/// `pressure_increase_rate * dt`.
pub fn abs_update(
    state: &mut AbsState,
    slip_ratios: &[f64; 4],
    dt: f64,
    config: &AbsConfig,
) -> [f64; 4] {
    let mut any_active = false;
    let mut pressures = [0.0f64; 4];
    for i in 0..4 {
        let slip = slip_ratios[i].abs();
        if slip > config.slip_target {
            any_active = true;
            state.pressure[i] = (state.pressure[i] - config.pressure_decrease_rate * dt).max(0.0);
        } else {
            state.pressure[i] = (state.pressure[i] + config.pressure_increase_rate * dt).min(1.0);
        }
        pressures[i] = state.pressure[i];
    }
    state.is_active = any_active;
    state.cycle_phase = if any_active {
        String::from("decrease")
    } else {
        String::from("increase")
    };
    pressures
}
/// Lateral (centripetal) acceleration for a vehicle in a turn.
///
/// `a_lat = v² / R`
pub fn lateral_acceleration(speed: f64, turn_radius: f64) -> f64 {
    if turn_radius.abs() < 1e-10 {
        return 0.0;
    }
    speed * speed / turn_radius
}
/// Lateral g-force threshold at which rollover occurs (in g).
///
/// `rollover_g = track_width / (2 · h_cg)`
pub fn rollover_threshold(track_width: f64, cg_height: f64) -> f64 {
    track_width / (2.0 * cg_height)
}
/// Static Stability Factor (SSF).
///
/// `SSF = track_width / (2 · h_cg)`
///
/// Higher SSF = more rollover resistance.
pub fn static_stability_factor(track_width: f64, cg_height: f64) -> f64 {
    track_width / (2.0 * cg_height)
}
/// Understeer gradient Kus (rad / (m/s²)) from vehicle and tire parameters.
///
/// `Kus = m · (b / (Cf · L) - a / (Cr · L))`
///
/// where `L = a + b`.
/// - Positive Kus → understeering vehicle.
/// - Negative Kus → oversteering vehicle.
pub fn understeer_gradient(m: f64, cf: f64, cr: f64, a: f64, b: f64) -> f64 {
    let l = a + b;
    m * (b / (cf * l) - a / (cr * l))
}
/// Critical speed for an oversteering vehicle (Kus < 0).
///
/// `v_cr = sqrt(-wheelbase / Kus)`
///
/// Returns `f64::INFINITY` when `Kus >= 0` (neutral or understeering).
pub fn critical_speed_oversteer(wheelbase: f64, kus: f64) -> f64 {
    if kus >= 0.0 {
        return f64::INFINITY;
    }
    (-wheelbase / kus).sqrt()
}
/// Lateral Load Transfer (LLT) distribution between axles during cornering.
///
/// Returns the rollover index for a given axle: 0 = no load transfer, 1 = rollover.
///
/// `LTD = (F_outer - F_inner) / (F_outer + F_inner)`
///
/// For a single axle: LTD = m * a_lat * h_cg / (track_width * normal_load)
pub fn lateral_load_transfer(
    lateral_accel: f64,
    mass: f64,
    cg_height: f64,
    track_width: f64,
    axle_normal_load: f64,
) -> f64 {
    if axle_normal_load < 1e-12 || track_width < 1e-12 {
        return 0.0;
    }
    let transfer = mass * lateral_accel * cg_height / track_width;
    (transfer / axle_normal_load).clamp(-1.0, 1.0)
}
/// Rollover index (RI) combining lateral load transfer from both axles.
///
/// `RI = max(|LTD_front|, |LTD_rear|)` in \[0, 1\]; 1 indicates impending rollover.
pub fn rollover_index(
    lateral_accel: f64,
    mass_front: f64,
    mass_rear: f64,
    cg_height: f64,
    track_width_front: f64,
    track_width_rear: f64,
    normal_load_front: f64,
    normal_load_rear: f64,
) -> f64 {
    let ltd_f = lateral_load_transfer(
        lateral_accel,
        mass_front,
        cg_height,
        track_width_front,
        normal_load_front,
    );
    let ltd_r = lateral_load_transfer(
        lateral_accel,
        mass_rear,
        cg_height,
        track_width_rear,
        normal_load_rear,
    );
    ltd_f.abs().max(ltd_r.abs()).min(1.0)
}
/// Compute active roll control torques to counteract body roll.
///
/// Returns the anti-roll bar torques for front and rear axles.
pub fn arc_compute(roll_angle: f64, roll_rate: f64, config: &ArcConfig) -> ArcOutput {
    let total = -(config.roll_gain * roll_angle + config.roll_rate_gain * roll_rate);
    let clamped = total.clamp(-config.max_torque, config.max_torque);
    ArcOutput {
        front_torque: clamped * config.front_fraction,
        rear_torque: clamped * (1.0 - config.front_fraction),
    }
}
/// Compute ESP intervention level from stability margins.
///
/// Returns the intervention level based on yaw rate error and sideslip.
pub fn esp_intervention_level(
    yaw_rate_error: f64,
    sideslip_deg: f64,
    yaw_threshold_light: f64,
    yaw_threshold_heavy: f64,
    sideslip_threshold_deg: f64,
) -> EspLevel {
    let yaw_abs = yaw_rate_error.abs();
    let ss_abs = sideslip_deg.abs();
    if yaw_abs >= yaw_threshold_heavy || ss_abs >= sideslip_threshold_deg * 1.5 {
        EspLevel::Heavy
    } else if yaw_abs >= yaw_threshold_light || ss_abs >= sideslip_threshold_deg {
        EspLevel::Light
    } else {
        EspLevel::None
    }
}
/// Compute torque vectoring adjustments to correct under/oversteer.
///
/// Positive `yaw_error` = oversteer (actual yaw > desired).
/// Uses a simple proportional strategy.
pub fn torque_vectoring(
    yaw_error: f64,
    base_torque: f64,
    tv_gain: f64,
    max_delta: f64,
) -> TorqueVectoringOutput {
    let delta = (tv_gain * yaw_error).clamp(-max_delta, max_delta);
    TorqueVectoringOutput {
        rear_left: base_torque * 0.5 + delta,
        rear_right: base_torque * 0.5 - delta,
    }
}
/// Slip angle (small-angle approximation) at an axle.
///
/// `alpha = atan(v_lateral / v_longitudinal)` ≈ `v_lat / v_lon` for small angles.
pub fn axle_slip_angle_rad(v_lateral: f64, v_longitudinal: f64) -> f64 {
    if v_longitudinal.abs() < 1e-6 {
        return 0.0;
    }
    (v_lateral / v_longitudinal).atan()
}
/// Compute the handling balance from front and rear slip angles.
///
/// `neutral_threshold_rad` defines the dead-band around zero.
pub fn handling_balance(
    alpha_front_rad: f64,
    alpha_rear_rad: f64,
    neutral_threshold_rad: f64,
) -> HandlingBalance {
    let delta = alpha_front_rad - alpha_rear_rad;
    if delta.abs() <= neutral_threshold_rad.abs() {
        HandlingBalance::Neutral
    } else if delta > 0.0 {
        HandlingBalance::Understeer
    } else {
        HandlingBalance::Oversteer
    }
}
/// Compute front and rear brake force fractions using EBD.
///
/// Returns `(front_fraction, rear_fraction)` summing to 1.0.
pub fn ebd_brake_split(decel_g: f64, config: &EbdConfig) -> (f64, f64) {
    let min_front = config.front_bias.min(config.max_front_fraction);
    let front = (config.front_bias + config.decel_bias_gain * decel_g.abs())
        .clamp(min_front, config.max_front_fraction);
    let rear = 1.0 - front;
    if rear < config.min_rear_fraction {
        let adjusted_front = 1.0 - config.min_rear_fraction;
        return (adjusted_front, config.min_rear_fraction);
    }
    (front, rear)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AbsController;
    use crate::ElectronicStabilityControl;
    use crate::TractionControl;
    use crate::stability::AbsEnhanced;
    use crate::stability::BicycleModel;
    use crate::stability::EmergencyBrakeAssist;
    use crate::stability::HillHold;
    use crate::stability::PitchRollState;
    use crate::stability::RolloverDetection;
    use crate::stability::RolloverPrevention;
    use crate::stability::StabilityController;
    use crate::stability::StabilityEnvelopeMonitor;
    use crate::stability::StabilityEnvelopeStatus;
    use crate::stability::YawRateObserver;
    #[test]
    fn test_slip_ratio_zero_vehicle_speed() {
        let slip = compute_slip_ratio(5.0, 0.0);
        assert_eq!(slip, 0.0, "zero vehicle speed should return 0");
    }
    #[test]
    fn test_slip_ratio_no_slip() {
        let slip = compute_slip_ratio(10.0, 10.0);
        assert!(slip.abs() < 1e-10, "equal speeds → slip = 0, got {slip}");
    }
    #[test]
    fn test_slip_ratio_spinning_wheel() {
        let slip = compute_slip_ratio(15.0, 10.0);
        assert!((slip - 0.5).abs() < 1e-10, "expected 0.5, got {slip}");
    }
    #[test]
    fn test_slip_ratio_clamped() {
        let slip = compute_slip_ratio(100.0, 1.0);
        assert_eq!(slip, 1.0, "slip should clamp to 1.0");
        let slip_neg = compute_slip_ratio(-100.0, 1.0);
        assert_eq!(slip_neg, -1.0, "negative slip should clamp to -1.0");
    }
    #[test]
    fn test_tcs_update_no_slip() {
        let config = TcsConfig::default();
        let mut state = TcsState::default();
        let factor = tcs_update(&mut state, 0.05, &config);
        assert!(!state.is_active);
        assert!(
            factor > 0.99,
            "factor should recover toward 1.0, got {factor}"
        );
    }
    #[test]
    fn test_tcs_update_slip() {
        let config = TcsConfig::default();
        let mut state = TcsState::default();
        let factor = tcs_update(&mut state, 0.5, &config);
        assert!(state.is_active);
        assert!(factor < 1.0, "factor should be reduced, got {factor}");
        assert!(
            factor >= config.min_torque_fraction,
            "factor should not go below min"
        );
    }
    #[test]
    fn test_tcs_min_torque_floor() {
        let config = TcsConfig::default();
        let mut state = TcsState {
            torque_reduction_factor: config.min_torque_fraction,
            ..Default::default()
        };
        let factor = tcs_update(&mut state, 0.9, &config);
        assert!(
            factor >= config.min_torque_fraction,
            "factor should not go below min_torque_fraction"
        );
    }
    #[test]
    fn test_abs_no_slip_no_reduction() {
        let abs = AbsController::new(0.2, 0.8);
        let torque = abs.apply(0.05, 1000.0);
        assert_eq!(torque, 1000.0);
    }
    #[test]
    fn test_abs_excessive_slip_reduces_torque() {
        let abs = AbsController::new(0.2, 0.8);
        let torque = abs.apply(0.35, 1000.0);
        assert!(
            (torque - 800.0).abs() < 1e-9,
            "expected 800.0, got {torque}"
        );
    }
    #[test]
    fn test_abs_is_active() {
        let abs = AbsController::new(0.2, 0.8);
        assert!(!abs.is_active(0.1));
        assert!(!abs.is_active(0.2));
        assert!(abs.is_active(0.21));
        assert!(abs.is_active(-0.25));
    }
    #[test]
    fn test_abs_update_reduces_pressure_on_lockup() {
        let config = AbsConfig::default();
        let mut state = AbsState::default();
        let slips = [0.5f64; 4];
        let pressures = abs_update(&mut state, &slips, 0.1, &config);
        assert!(state.is_active);
        for p in pressures.iter() {
            assert!(*p < 1.0, "pressure should be reduced, got {p}");
        }
    }
    #[test]
    fn test_abs_update_increases_pressure_when_normal() {
        let config = AbsConfig::default();
        let mut state = AbsState {
            pressure: vec![0.5, 0.5, 0.5, 0.5],
            ..Default::default()
        };
        let slips = [0.01f64; 4];
        let pressures = abs_update(&mut state, &slips, 0.1, &config);
        assert!(!state.is_active);
        for p in pressures.iter() {
            assert!(*p > 0.5, "pressure should increase, got {p}");
        }
    }
    #[test]
    fn test_abs_update_cycle_phase() {
        let config = AbsConfig::default();
        let mut state = AbsState::default();
        abs_update(&mut state, &[0.5, 0.5, 0.5, 0.5], 0.05, &config);
        assert_eq!(state.cycle_phase, "decrease");
        abs_update(&mut state, &[0.0, 0.0, 0.0, 0.0], 0.05, &config);
        assert_eq!(state.cycle_phase, "increase");
    }
    #[test]
    fn test_traction_no_spin_full_torque() {
        let tc = TractionControl::new(0.15, 0.7);
        let torque = tc.apply(0.05, 500.0);
        assert_eq!(torque, 500.0);
    }
    #[test]
    fn test_traction_wheel_spin_reduces_torque() {
        let tc = TractionControl::new(0.15, 0.7);
        let torque = tc.apply(0.3, 500.0);
        assert!(
            (torque - 350.0).abs() < 1e-9,
            "expected 350.0, got {torque}"
        );
    }
    #[test]
    fn test_esc_oversteer_correction() {
        let esc = ElectronicStabilityControl::new(0.05);
        let (fl, fr, rl, rr) = esc.correction_torques(0.5, 0.1, 1000.0);
        assert!(fl > 0.0);
        assert_eq!(fr, 0.0);
        assert!(rl > 0.0);
        assert_eq!(rr, 0.0);
    }
    #[test]
    fn test_esc_understeer_correction() {
        let esc = ElectronicStabilityControl::new(0.05);
        let (fl, fr, rl, rr) = esc.correction_torques(0.1, 0.5, 1000.0);
        assert_eq!(fl, 0.0);
        assert!(fr > 0.0);
        assert_eq!(rl, 0.0);
        assert!(rr > 0.0);
    }
    #[test]
    fn test_esc_within_threshold_no_correction() {
        let esc = ElectronicStabilityControl::new(0.1);
        let result = esc.correction_torques(0.15, 0.1, 1000.0);
        assert_eq!(result, (0.0, 0.0, 0.0, 0.0));
    }
    #[test]
    fn test_esc_update_no_intervention() {
        let config = EscConfig::default();
        let mut state = EscState::default();
        let intervention = esc_update(&mut state, 0.1, 0.1, 1.0, &config);
        assert!(!state.is_active);
        assert_eq!(intervention.brake_front_left, 0.0);
        assert_eq!(intervention.throttle_reduction, 0.0);
    }
    #[test]
    fn test_esc_update_oversteer() {
        let config = EscConfig::default();
        let mut state = EscState::default();
        let intervention = esc_update(&mut state, 1.0, 0.0, 0.0, &config);
        assert!(state.is_active);
        assert_eq!(state.intervention_type, "oversteer");
        assert!(intervention.brake_front_left > 0.0);
        assert_eq!(intervention.brake_front_right, 0.0);
    }
    #[test]
    fn test_esc_update_understeer() {
        let config = EscConfig::default();
        let mut state = EscState::default();
        let intervention = esc_update(&mut state, 0.0, 1.0, 0.0, &config);
        assert!(state.is_active);
        assert_eq!(state.intervention_type, "understeer");
        assert!(intervention.brake_front_right > 0.0);
        assert_eq!(intervention.brake_front_left, 0.0);
    }
    #[test]
    fn test_target_yaw_rate_zero_steer() {
        let r = target_yaw_rate(20.0, 0.0, 2.5, 0.003);
        assert_eq!(r, 0.0);
    }
    #[test]
    fn test_target_yaw_rate_positive() {
        let r = target_yaw_rate(10.0, 5.0, 2.5, 0.0);
        let expected = 10.0 / 2.5 * 5.0f64.to_radians();
        assert!((r - expected).abs() < 1e-9, "expected {expected}, got {r}");
    }
    #[test]
    fn test_lateral_acceleration() {
        let a = lateral_acceleration(20.0, 100.0);
        assert!((a - 4.0).abs() < 1e-10, "expected 4.0 m/s², got {a}");
    }
    #[test]
    fn test_lateral_acceleration_zero_radius() {
        let a = lateral_acceleration(20.0, 0.0);
        assert_eq!(a, 0.0);
    }
    #[test]
    fn test_rollover_threshold() {
        let thresh = rollover_threshold(1.5, 0.5);
        assert!((thresh - 1.5).abs() < 1e-10, "expected 1.5, got {thresh}");
    }
    #[test]
    fn test_static_stability_factor() {
        let ssf = static_stability_factor(1.6, 0.5);
        let expected = 1.6 / (2.0 * 0.5);
        assert!((ssf - expected).abs() < 1e-10);
    }
    #[test]
    fn test_understeer_gradient_neutral() {
        let m = 1500.0;
        let a = 1.2;
        let b = 1.3;
        let l = a + b;
        let cr = 50000.0;
        let cf = cr * b / a;
        let kus = understeer_gradient(m, cf, cr, a, b);
        let expected = m * (b / (cf * l) - a / (cr * l));
        assert!((kus - expected).abs() < 1e-6);
    }
    #[test]
    fn test_understeer_gradient_sign() {
        let kus_under = understeer_gradient(1500.0, 80000.0, 60000.0, 1.0, 1.5);
        let kus_over = understeer_gradient(1500.0, 60000.0, 80000.0, 1.5, 1.0);
        assert!(
            kus_under > kus_over,
            "understeer should have larger Kus than oversteer"
        );
    }
    #[test]
    fn test_critical_speed_oversteer() {
        let kus = -0.002;
        let v = critical_speed_oversteer(2.5, kus);
        let expected = (-2.5 / -0.002f64).sqrt();
        assert!((v - expected).abs() < 1e-9);
    }
    #[test]
    fn test_critical_speed_understeer_returns_infinity() {
        let v = critical_speed_oversteer(2.5, 0.003);
        assert!(v.is_infinite());
    }
    #[test]
    fn test_bicycle_model_new() {
        let m = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let expected_iz = 1500.0 * (1.2 * 1.2 + 1.4 * 1.4) / 3.0;
        assert!((m.iz - expected_iz).abs() < 1e-6);
    }
    #[test]
    fn test_bicycle_model_lateral_forces_straight() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let (fyf, fyr) = bm.lateral_forces(20.0, 0.0, 0.0, 0.0);
        assert!(fyf.abs() < 1e-10);
        assert!(fyr.abs() < 1e-10);
    }
    #[test]
    fn test_bicycle_model_lateral_forces_steer() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let delta = 5.0f64.to_radians();
        let (fyf, fyr) = bm.lateral_forces(20.0, 0.0, 0.0, delta);
        assert!(fyf > 0.0, "Fyf should be positive, got {fyf}");
        assert!(fyr.abs() < 1e-10, "Fyr should be ~0, got {fyr}");
    }
    #[test]
    fn test_bicycle_model_step_produces_yaw() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let delta = 5.0f64.to_radians();
        let (_beta_new, r_new) = bm.step(20.0, 0.0, 0.0, delta, 0.01);
        assert!(r_new.abs() > 0.0, "yaw rate should develop, got {r_new}");
    }
    #[test]
    fn test_bicycle_model_zero_speed() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let (fyf, fyr) = bm.lateral_forces(0.0, 0.1, 0.2, 0.1);
        assert_eq!(fyf, 0.0);
        assert_eq!(fyr, 0.0);
    }
    #[test]
    fn test_pitch_roll_state_integrate_zero() {
        let mut s = PitchRollState::default();
        s.integrate(0.0, 0.0, 0.1);
        assert_eq!(s.roll, 0.0);
        assert_eq!(s.pitch, 0.0);
    }
    #[test]
    fn test_pitch_roll_state_integrate_nonzero() {
        let mut s = PitchRollState::default();
        s.integrate(0.1, 0.05, 1.0);
        assert!((s.roll - 0.1).abs() < 1e-12);
        assert!((s.pitch - 0.05).abs() < 1e-12);
    }
    #[test]
    fn test_pitch_roll_state_accel_correct() {
        let mut s = PitchRollState {
            roll: 0.5,
            pitch: 0.3,
            ..Default::default()
        };
        s.correct_with_accel(0.0, 0.0, 0.0);
        assert!(
            s.roll.abs() < 1e-12,
            "roll should be corrected to 0, got {}",
            s.roll
        );
        assert!(s.pitch.abs() < 1e-12);
    }
    #[test]
    fn test_rollover_detection_no_risk() {
        let det = RolloverDetection::default_car();
        assert!(!det.is_at_risk(2.0, 0.1));
    }
    #[test]
    fn test_rollover_detection_lat_accel_risk() {
        let det = RolloverDetection::default_car();
        assert!(det.is_at_risk(10.0, 0.0));
    }
    #[test]
    fn test_rollover_detection_roll_angle_risk() {
        let det = RolloverDetection::default_car();
        assert!(det.is_at_risk(1.0, 0.5));
    }
    #[test]
    fn test_rollover_risk_level_zero_when_safe() {
        let det = RolloverDetection::default_car();
        let risk = det.risk_level(0.0, 0.0);
        assert_eq!(risk, 0.0);
    }
    #[test]
    fn test_rollover_risk_level_clamped_at_one() {
        let det = RolloverDetection::default_car();
        let risk = det.risk_level(100.0, 10.0);
        assert!(
            (risk - 1.0).abs() < 1e-12,
            "risk should clamp to 1.0, got {risk}"
        );
    }
    #[test]
    fn test_stability_controller_default_inactive() {
        let ctrl = StabilityController::new();
        assert!(!ctrl.any_active());
    }
    #[test]
    fn test_stability_controller_tcs_activates_on_slip() {
        let mut ctrl = StabilityController::new();
        let (drive, _, _) = ctrl.update(0.8, &[0.0; 4], 500.0, 0.0, 0.0, 0.0, 0.0, 0.01);
        assert!(ctrl.tcs_state.is_active, "TCS should be active");
        assert!(drive < 500.0, "drive torque should be reduced");
    }
    #[test]
    fn test_stability_controller_abs_activates_on_lockup() {
        let mut ctrl = StabilityController::new();
        let (_, pressures, _) = ctrl.update(0.0, &[0.5; 4], 0.0, 1000.0, 0.0, 0.0, 0.0, 0.1);
        assert!(ctrl.abs_state.is_active, "ABS should be active");
        for p in pressures.iter() {
            assert!(*p < 1.0, "brake pressure should be reduced");
        }
    }
    #[test]
    fn test_stability_controller_esc_activates_on_yaw_error() {
        let mut ctrl = StabilityController::new();
        let (_, _, esc) = ctrl.update(0.0, &[0.0; 4], 0.0, 0.0, 1.0, 0.0, 0.0, 0.01);
        assert!(ctrl.esc_state.is_active, "ESC should be active");
        assert!(
            esc.brake_front_left > 0.0,
            "oversteer: FL brake should apply"
        );
    }
    #[test]
    fn test_lateral_load_transfer_zero() {
        let ltd = lateral_load_transfer(0.0, 750.0, 0.5, 1.5, 7500.0);
        assert!(ltd.abs() < 1e-12, "zero lat accel → LTD = 0");
    }
    #[test]
    fn test_lateral_load_transfer_clamped() {
        let ltd = lateral_load_transfer(100.0, 1500.0, 0.5, 1.5, 100.0);
        assert!((ltd - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rollover_index_range() {
        let ri = rollover_index(5.0, 750.0, 750.0, 0.5, 1.5, 1.5, 7500.0, 7500.0);
        assert!(
            (0.0..=1.0).contains(&ri),
            "rollover index should be in [0,1]"
        );
    }
    #[test]
    fn test_arc_no_roll_zero_torque() {
        let config = ArcConfig::default();
        let out = arc_compute(0.0, 0.0, &config);
        assert!(out.front_torque.abs() < 1e-12);
        assert!(out.rear_torque.abs() < 1e-12);
    }
    #[test]
    fn test_arc_positive_roll_opposes_motion() {
        let config = ArcConfig::default();
        let out = arc_compute(0.05, 0.1, &config);
        assert!(out.front_torque < 0.0, "front torque should oppose roll");
        assert!(out.rear_torque < 0.0, "rear torque should oppose roll");
    }
    #[test]
    fn test_arc_torque_clamped() {
        let config = ArcConfig::default();
        let out = arc_compute(10.0, 10.0, &config);
        assert!(out.front_torque.abs() <= config.max_torque);
        assert!(out.rear_torque.abs() <= config.max_torque);
    }
    #[test]
    fn test_arc_front_rear_split() {
        let config = ArcConfig {
            front_fraction: 0.7,
            ..ArcConfig::default()
        };
        let out = arc_compute(0.1, 0.0, &config);
        let ratio = if out.rear_torque.abs() > 1e-12 {
            out.front_torque / out.rear_torque
        } else {
            1.0
        };
        assert!(
            (ratio - 0.7 / 0.3).abs() < 1e-6,
            "front/rear ratio should match fraction"
        );
    }
    #[test]
    fn test_yaw_observer_default_state() {
        let obs = YawRateObserver::new(1.0, 2.0);
        assert!((obs.beta_est).abs() < 1e-12);
        assert!((obs.r_est).abs() < 1e-12);
    }
    #[test]
    fn test_yaw_observer_tracks_measurement() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let mut obs = YawRateObserver::new(2.0, 5.0);
        for _ in 0..50 {
            obs.update(&bm, 20.0, 0.0, 0.3, 0.01);
        }
        assert!(
            (obs.r_est - 0.3).abs() < 0.2,
            "observer should roughly track r_meas=0.3, got {}",
            obs.r_est
        );
    }
    #[test]
    fn test_yaw_observer_zero_speed_no_update() {
        let bm = BicycleModel::new(1500.0, 80000.0, 70000.0, 1.2, 1.4);
        let mut obs = YawRateObserver::new(1.0, 1.0);
        obs.beta_est = 0.1;
        obs.r_est = 0.2;
        obs.update(&bm, 0.0, 0.0, 0.5, 0.01);
        assert!((obs.beta_est - 0.1).abs() < 1e-12);
        assert!((obs.r_est - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_esp_level_none() {
        let level = esp_intervention_level(0.01, 2.0, 0.05, 0.15, 8.0);
        assert_eq!(level, EspLevel::None);
    }
    #[test]
    fn test_esp_level_light() {
        let level = esp_intervention_level(0.08, 3.0, 0.05, 0.15, 8.0);
        assert_eq!(level, EspLevel::Light);
    }
    #[test]
    fn test_esp_level_heavy_yaw() {
        let level = esp_intervention_level(0.5, 3.0, 0.05, 0.15, 8.0);
        assert_eq!(level, EspLevel::Heavy);
    }
    #[test]
    fn test_esp_level_heavy_sideslip() {
        let level = esp_intervention_level(0.01, 14.0, 0.05, 0.15, 8.0);
        assert_eq!(level, EspLevel::Heavy);
    }
    #[test]
    fn test_torque_vectoring_zero_error() {
        let out = torque_vectoring(0.0, 200.0, 500.0, 50.0);
        assert!((out.rear_left - 100.0).abs() < 1e-10);
        assert!((out.rear_right - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_torque_vectoring_oversteer_shifts_right() {
        let out = torque_vectoring(0.05, 200.0, 500.0, 50.0);
        assert!(
            out.rear_left > out.rear_right,
            "oversteer: left wheel should have more torque, got L={} R={}",
            out.rear_left,
            out.rear_right
        );
    }
    #[test]
    fn test_torque_vectoring_clamped() {
        let out = torque_vectoring(100.0, 200.0, 500.0, 30.0);
        assert!((out.rear_left - (100.0 + 30.0)).abs() < 1e-10);
        assert!((out.rear_right - (100.0 - 30.0)).abs() < 1e-10);
    }
    #[test]
    fn test_stability_envelope_normal() {
        let mon = StabilityEnvelopeMonitor::default();
        let status = mon.evaluate(1.5, 0.1, 3.0);
        assert_eq!(status, StabilityEnvelopeStatus::Normal);
    }
    #[test]
    fn test_stability_envelope_warning_slip() {
        let mon = StabilityEnvelopeMonitor::default();
        let status = mon.evaluate(1.5, 0.25, 3.0);
        assert_eq!(status, StabilityEnvelopeStatus::Warning);
    }
    #[test]
    fn test_stability_envelope_critical_rollover() {
        let mon = StabilityEnvelopeMonitor::default();
        let status = mon.evaluate(0.5, 0.1, 3.0);
        assert_eq!(status, StabilityEnvelopeStatus::Critical);
    }
    #[test]
    fn test_stability_score_normal() {
        let mon = StabilityEnvelopeMonitor::default();
        let score = mon.stability_score(2.0, 0.05, 2.0);
        assert!(
            score > 0.5,
            "normal conditions should give score > 0.5, got {score}"
        );
    }
    #[test]
    fn test_stability_score_range() {
        let mon = StabilityEnvelopeMonitor::default();
        for (ssf, slip, ss) in [(2.0, 0.1, 5.0), (0.8, 0.5, 20.0), (1.5, 0.25, 10.0)] {
            let score = mon.stability_score(ssf, slip, ss);
            assert!((0.0..=1.0).contains(&score), "score out of [0,1]: {score}");
        }
    }
    #[test]
    fn test_handling_balance_neutral() {
        let bal = handling_balance(0.05, 0.04, 0.02);
        assert_eq!(bal, HandlingBalance::Neutral);
    }
    #[test]
    fn test_handling_balance_understeer() {
        let bal = handling_balance(0.15, 0.03, 0.02);
        assert_eq!(bal, HandlingBalance::Understeer);
    }
    #[test]
    fn test_handling_balance_oversteer() {
        let bal = handling_balance(0.02, 0.18, 0.02);
        assert_eq!(bal, HandlingBalance::Oversteer);
    }
    #[test]
    fn test_axle_slip_angle_zero_speed() {
        let alpha = axle_slip_angle_rad(2.0, 0.0);
        assert_eq!(alpha, 0.0, "zero longitudinal speed → slip angle = 0");
    }
    #[test]
    fn test_axle_slip_angle_positive_lateral() {
        let alpha = axle_slip_angle_rad(1.0, 10.0);
        assert!(alpha > 0.0, "positive lateral speed → positive slip angle");
        assert!(
            alpha < std::f64::consts::FRAC_PI_2,
            "should be less than 90°"
        );
    }
    #[test]
    fn test_ebd_sums_to_one() {
        let config = EbdConfig::default();
        let (front, rear) = ebd_brake_split(0.5, &config);
        assert!(
            (front + rear - 1.0).abs() < 1e-12,
            "front + rear must sum to 1.0"
        );
    }
    #[test]
    fn test_ebd_front_increases_with_decel() {
        let config = EbdConfig::default();
        let (f0, _) = ebd_brake_split(0.0, &config);
        let (f1, _) = ebd_brake_split(1.0, &config);
        assert!(f1 >= f0, "higher decel should shift more braking to front");
    }
    #[test]
    fn test_ebd_rear_minimum_enforced() {
        let config = EbdConfig::default();
        let (_, rear) = ebd_brake_split(10.0, &config);
        assert!(
            rear >= config.min_rear_fraction - 1e-12,
            "rear fraction must not drop below minimum"
        );
    }
    #[test]
    fn test_eba_not_active_slow_application() {
        let mut eba = EmergencyBrakeAssist::default_car();
        let _p1 = eba.update(0.0, 1.0);
        let p2 = eba.update(0.1, 1.0);
        assert!(!eba.is_active, "slow application should not trigger EBA");
        assert!((p2 - 0.1).abs() < 1e-9, "EBA inactive: pressure unchanged");
    }
    #[test]
    fn test_eba_activates_on_fast_application() {
        let mut eba = EmergencyBrakeAssist::default_car();
        let _p1 = eba.update(0.0, 0.01);
        let p2 = eba.update(0.8, 0.01);
        assert!(eba.is_active, "fast application should trigger EBA");
        assert!(p2 > 0.8, "EBA active: pressure should be amplified");
    }
    #[test]
    fn test_eba_deactivates_on_release() {
        let mut eba = EmergencyBrakeAssist::default_car();
        eba.is_active = true;
        let _p = eba.update(0.0, 0.1);
        assert!(!eba.is_active, "EBA should deactivate when pedal released");
    }
    #[test]
    fn test_hill_hold_engages_on_slope() {
        let mut hh = HillHold::default_car();
        let pressure = hh.update(5.0, 0.0, 0.0);
        assert!(hh.is_engaged, "should engage on 5° slope when stationary");
        assert!(
            pressure > 0.0,
            "engaged hill hold should provide brake pressure"
        );
    }
    #[test]
    fn test_hill_hold_not_engaged_on_flat() {
        let mut hh = HillHold::default_car();
        let pressure = hh.update(1.0, 0.0, 0.0);
        assert!(!hh.is_engaged, "should not engage on 1° (below threshold)");
        assert_eq!(pressure, 0.0);
    }
    #[test]
    fn test_hill_hold_releases_with_throttle() {
        let mut hh = HillHold::default_car();
        hh.update(5.0, 0.0, 0.0);
        assert!(hh.is_engaged);
        let pressure = hh.update(5.0, 0.0, 0.5);
        assert!(
            !hh.is_engaged,
            "hill hold should release when throttle applied"
        );
        assert_eq!(pressure, 0.0);
    }
    #[test]
    fn test_rollover_prevention_not_active_below_threshold() {
        let rp = RolloverPrevention::default_suv();
        assert!(!rp.is_active(5.0), "5 m/s² is below 7 m/s² threshold");
        let (thr_cut, brake) = rp.evaluate(5.0);
        assert_eq!(thr_cut, 0.0);
        assert_eq!(brake, 0.0);
    }
    #[test]
    fn test_rollover_prevention_active_above_threshold() {
        let rp = RolloverPrevention::default_suv();
        assert!(rp.is_active(9.0), "9 m/s² exceeds 7 m/s² threshold");
        let (thr_cut, brake) = rp.evaluate(9.0);
        assert!(thr_cut > 0.0, "throttle cut should be applied");
        assert!(brake > 0.0, "brake should be applied");
    }
    #[test]
    fn test_rollover_prevention_symmetric() {
        let rp = RolloverPrevention::default_suv();
        let (_, brake_pos) = rp.evaluate(9.0);
        let (_, brake_neg) = rp.evaluate(-9.0);
        assert!(
            (brake_pos - brake_neg).abs() < 1e-12,
            "symmetric for left/right corners"
        );
    }
    #[test]
    fn test_abs_enhanced_no_slip_full_pressure() {
        let mut abs = AbsEnhanced::default_car();
        let slips = [0.0f64; 4];
        let speeds = [100.0f64; 4];
        let pressures = abs.update(&slips, &speeds, 0.01);
        for p in pressures.iter() {
            assert!(*p >= 1.0 - 1e-9, "no slip → full pressure, got {p}");
        }
    }
    #[test]
    fn test_abs_enhanced_activates_on_slip() {
        let mut abs = AbsEnhanced::default_car();
        let slips = [0.5f64; 4];
        let speeds = [100.0f64; 4];
        let pressures = abs.update(&slips, &speeds, 0.01);
        assert!(abs.any_active(), "ABS should be active");
        for p in pressures.iter() {
            assert!(*p < 1.0, "pressure should be reduced, got {p}");
        }
    }
    #[test]
    fn test_abs_enhanced_activates_on_decel() {
        let mut abs = AbsEnhanced::default_car();
        abs.prev_wheel_speeds = [200.0; 4];
        let slips = [0.0f64; 4];
        let speeds = [50.0f64; 4];
        let pressures = abs.update(&slips, &speeds, 0.01);
        assert!(
            abs.any_active(),
            "high wheel deceleration should trigger ABS"
        );
        for p in pressures.iter() {
            assert!(
                *p < 1.0,
                "pressure should be reduced on high decel, got {p}"
            );
        }
    }
}
#[cfg(test)]
mod tests_stability_new {

    use crate::TractionControl;

    use crate::stability::AbsPressureModulator;

    use crate::stability::EscYawMoment;

    use crate::stability::RolloverDetection;

    #[test]
    fn test_tc_slip_ratio_zero_at_zero_vehicle_speed() {
        let slip = TractionControl::compute_slip_ratio(100.0, 0.3, 0.0);
        assert_eq!(slip, 0.0, "zero vehicle speed → slip = 0");
    }
    #[test]
    fn test_tc_slip_ratio_positive_when_spinning() {
        let slip = TractionControl::compute_slip_ratio(50.0, 0.3, 10.0);
        assert!((slip - 0.5).abs() < 1e-9, "slip should be 0.5, got {slip}");
    }
    #[test]
    fn test_tc_slip_ratio_clamped_to_minus_one() {
        let slip = TractionControl::compute_slip_ratio(0.0, 0.3, 30.0);
        assert!(
            (slip - (-1.0)).abs() < 1e-9,
            "locked wheel → slip = -1, got {slip}"
        );
    }
    #[test]
    fn test_abs_modulator_dumps_pressure_on_high_slip() {
        let mut abs = AbsPressureModulator::default_car();
        abs.pressure = 1.0;
        let p = abs.modulate_brake_pressure(0.30, 0.01);
        assert!(p < 1.0, "should dump pressure on high slip, got {p}");
    }
    #[test]
    fn test_abs_modulator_builds_pressure_on_low_slip() {
        let mut abs = AbsPressureModulator::default_car();
        abs.pressure = 0.5;
        let p = abs.modulate_brake_pressure(0.01, 0.1);
        assert!(p > 0.5, "should build pressure on low slip, got {p}");
    }
    #[test]
    fn test_abs_modulator_pressure_stays_bounded() {
        let mut abs = AbsPressureModulator::default_car();
        abs.pressure = 1.0;
        for _ in 0..100 {
            abs.modulate_brake_pressure(0.0, 0.01);
        }
        assert!(abs.pressure <= 1.0, "pressure must not exceed 1.0");
        for _ in 0..100 {
            abs.modulate_brake_pressure(0.5, 0.01);
        }
        assert!(abs.pressure >= 0.0, "pressure must not go below 0");
    }
    #[test]
    fn test_esc_yaw_moment_zero_when_no_error() {
        let esc = EscYawMoment::default_sedan();
        let mz = esc.compute_yaw_moment_demand(0.0, 0.0, 0.0);
        assert_eq!(mz, 0.0, "no error → zero yaw moment");
    }
    #[test]
    fn test_esc_yaw_moment_corrects_oversteer() {
        let esc = EscYawMoment::default_sedan();
        let mz = esc.compute_yaw_moment_demand(0.4, 0.2, 0.0);
        assert!(
            mz > 0.0,
            "oversteer correction should be positive, got {mz}"
        );
    }
    #[test]
    fn test_esc_yaw_moment_clamped_to_max() {
        let esc = EscYawMoment::default_sedan();
        let mz = esc.compute_yaw_moment_demand(10.0, 0.0, 0.0);
        assert!(
            (mz - esc.max_moment).abs() < 1e-6,
            "should clamp to max_moment: {mz}"
        );
    }
    #[test]
    fn test_ttl_infinity_when_not_building() {
        let rd = RolloverDetection::default_car();
        let ttl = rd.compute_ttl(2.0, 0.0, 0.0);
        assert!(ttl.is_infinite(), "zero accel rate → TTL = infinity: {ttl}");
    }
    #[test]
    fn test_ttl_finite_when_approaching_threshold() {
        let rd = RolloverDetection::default_car();
        let ttl = rd.compute_ttl(5.0, 2.0, 0.0);
        assert!(
            ttl.is_finite() && ttl > 0.0,
            "should return finite positive TTL: {ttl}"
        );
    }
    #[test]
    fn test_ttl_zero_when_over_threshold() {
        let rd = RolloverDetection::default_car();
        let ttl = rd.compute_ttl(rd.g_threshold + 1.0, 1.0, 0.0);
        assert_eq!(ttl, 0.0, "over threshold → TTL = 0");
    }
}
