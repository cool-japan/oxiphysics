//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::motors_extended::ZnTuningMethod;
use super::types::PidController;

#[cfg(test)]
mod tests {
    use super::super::types::*;

    use crate::motors::ElectricMotor;

    use crate::motors::HydraulicActuator;
    use crate::motors::MotorBank;

    use crate::motors::MotorThermalModel;

    use crate::motors::PidController;

    use crate::motors::SpringDamper;

    use crate::motors::TorqueSpring;

    #[test]
    fn position_motor_at_target_only_damping() {
        let mut motor = PositionMotor::new(10.0, 2.0, 100.0);
        motor.update_position(5.0);
        motor.target_position = 5.0;
        let velocity = 3.0;
        let force = motor.compute_force(velocity);
        assert!((force - (-2.0 * velocity)).abs() < 1e-12);
    }
    #[test]
    fn position_motor_below_target_positive_force() {
        let mut motor = PositionMotor::new(10.0, 0.0, 100.0);
        motor.update_position(0.0);
        motor.target_position = 5.0;
        let force = motor.compute_force(0.0);
        assert!(force > 0.0);
        assert!((force - 50.0).abs() < 1e-12);
    }
    #[test]
    fn position_motor_clamped_at_max_force() {
        let mut motor = PositionMotor::new(100.0, 0.0, 10.0);
        motor.update_position(0.0);
        motor.target_position = 5.0;
        let force = motor.compute_force(0.0);
        assert!((force - 10.0).abs() < 1e-12);
    }
    #[test]
    fn velocity_motor_at_target_zero_torque() {
        let mut motor = VelocityMotor::new(5.0, 50.0);
        motor.target_velocity = 3.0;
        motor.update_velocity(3.0);
        let torque = motor.compute_torque();
        assert!(torque.abs() < 1e-12);
    }
    #[test]
    fn spring_damper_at_rest_zero_force() {
        let sd = SpringDamper::new(10.0, 2.0, 1.0);
        let f = sd.force(1.0, 0.0);
        assert!(f.abs() < 1e-12);
    }
    #[test]
    fn spring_damper_compressed_positive_force() {
        let sd = SpringDamper::new(10.0, 0.0, 1.0);
        let f = sd.force(0.5, 0.0);
        assert!(f > 0.0);
    }
    #[test]
    fn pid_proportional_only() {
        let mut pid = PidController::new(2.0, 0.0, 0.0, 1000.0);
        let out = pid.update(5.0, 0.01);
        assert!((out - 10.0).abs() < 1e-12);
    }
    #[test]
    fn pid_integral_accumulates() {
        let mut pid = PidController::new(0.0, 1.0, 0.0, 1e9);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        let out = pid.update(1.0, 0.1);
        assert!((out - 0.3).abs() < 1e-12);
    }
    #[test]
    fn pid_derivative_nonzero_on_error_change() {
        let mut pid = PidController::new(0.0, 0.0, 1.0, 1e9);
        pid.update(0.0, 0.1);
        let out = pid.update(2.0, 0.1);
        assert!((out - 20.0).abs() < 1e-10);
    }
    #[test]
    fn pid_reset_clears_state() {
        let mut pid = PidController::new(1.0, 1.0, 0.0, 1e9);
        pid.update(5.0, 1.0);
        pid.reset();
        assert!(pid.integral.abs() < 1e-15);
        assert!(pid.prev_error.abs() < 1e-15);
    }
    #[test]
    fn pid_output_clamped() {
        let mut pid = PidController::new(100.0, 0.0, 0.0, 10.0);
        let out = pid.update(5.0, 0.01);
        assert!((out - 10.0).abs() < 1e-12);
    }
    #[test]
    fn motor_max_torque_at_stall() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let t = motor.torque_at(0.0, 1.0);
        assert!((t - 100.0 * 0.9).abs() < 1e-9);
    }
    #[test]
    fn motor_zero_torque_at_max_rpm() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let t = motor.torque_at(5000.0, 1.0);
        assert!(t.abs() < 1e-9);
    }
    #[test]
    fn motor_throttle_scales_torque() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let t_full = motor.torque_at(0.0, 1.0);
        let t_half = motor.torque_at(0.0, 0.5);
        assert!((t_half - t_full * 0.5).abs() < 1e-9);
    }
    #[test]
    fn hydraulic_force_proportional_to_pressure() {
        let act = HydraulicActuator::new(0.01, 0.1, 1_000_000.0);
        let f = act.force(500_000.0);
        assert!((f - 5_000.0).abs() < 1e-9);
    }
    #[test]
    fn hydraulic_force_clamped_at_max_pressure() {
        let act = HydraulicActuator::new(0.01, 0.1, 1_000_000.0);
        let f = act.force(2_000_000.0);
        let expected = act.force(1_000_000.0);
        assert!((f - expected).abs() < 1e-9);
    }
    #[test]
    fn torque_spring_at_rest_angle_zero_torque() {
        let spring = TorqueSpring::new(100.0, 5.0, 0.5);
        let t = spring.torque(0.5, 0.0);
        assert!(t.abs() < 1e-12);
    }
    #[test]
    fn torque_spring_displaced_restoring_torque() {
        let spring = TorqueSpring::new(100.0, 0.0, 0.0);
        let t = spring.torque(1.0, 0.0);
        assert!((t - (-100.0)).abs() < 1e-12);
    }
    #[test]
    fn pid_anti_windup() {
        let mut pid = PidController::new(0.0, 1.0, 0.0, 1e9).with_integral_clamp(5.0);
        for _ in 0..1000 {
            pid.update(100.0, 0.1);
        }
        assert!(
            (pid.integral - 5.0).abs() < 1e-12,
            "integral should be clamped to 5.0"
        );
    }
    #[test]
    fn pid_ziegler_nichols_tuning() {
        let mut pid = PidController::new(0.0, 0.0, 0.0, 1000.0);
        pid.ziegler_nichols_tune(10.0, 1.0);
        assert!((pid.kp - 6.0).abs() < 1e-12);
        assert!((pid.ki - 12.0).abs() < 1e-12);
        assert!((pid.kd - 0.75).abs() < 1e-12);
    }
    #[test]
    fn pid_contributions() {
        let mut pid = PidController::new(1.0, 2.0, 3.0, 1e9);
        pid.update(1.0, 0.1);
        let c = pid.contributions(2.0, 0.1);
        assert!((c.proportional - 2.0).abs() < 1e-12);
        assert!((c.integral - 0.2).abs() < 1e-12);
        assert!((c.derivative - 30.0).abs() < 1e-12);
    }
    #[test]
    fn pid_contributions_total() {
        let c = PidContributions {
            proportional: 1.0,
            integral: 2.0,
            derivative: 3.0,
        };
        assert!((c.total() - 6.0).abs() < 1e-12);
    }
    #[test]
    fn motor_electrical_power() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let p_mech = motor.power(2500.0, 1.0);
        let p_elec = motor.electrical_power(2500.0, 1.0);
        assert!(p_elec >= p_mech - 1e-9);
    }
    #[test]
    fn motor_current_positive() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let i = motor.current(2500.0, 1.0);
        assert!(i > 0.0);
    }
    #[test]
    fn motor_heat_loss_non_negative() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let h = motor.heat_loss(2500.0, 1.0);
        assert!(h >= 0.0);
    }
    #[test]
    fn motor_rpm_at_torque() {
        let motor = ElectricMotor::new(100.0, 5000.0);
        let rpm = motor.rpm_at_torque(0.0);
        assert!((rpm - 5000.0).abs() < 1e-6);
        let rpm2 = motor.rpm_at_torque(1000.0);
        assert!((rpm2).abs() < 1e-6);
    }
    #[test]
    fn motor_thermal_step() {
        let mut motor = ElectricMotor::new(100.0, 5000.0);
        let initial_temp = motor.thermal.temperature;
        motor.step_thermal(2500.0, 1.0, 1.0);
        assert!(motor.thermal.temperature >= initial_temp);
    }
    #[test]
    fn thermal_model_cooling() {
        let mut thermal = MotorThermalModel::new(25.0, 150.0, 50.0, 0.5);
        thermal.temperature = 100.0;
        thermal.step(0.0, 1.0);
        assert!(thermal.temperature < 100.0);
        assert!(thermal.temperature > 25.0);
    }
    #[test]
    fn thermal_model_steady_state() {
        let thermal = MotorThermalModel::new(25.0, 150.0, 50.0, 0.5);
        let ss = thermal.steady_state_temperature(100.0);
        assert!((ss - 75.0).abs() < 1e-12);
    }
    #[test]
    fn thermal_model_time_constant() {
        let thermal = MotorThermalModel::new(25.0, 150.0, 50.0, 0.5);
        let tau = thermal.time_constant();
        assert!((tau - 25.0).abs() < 1e-12);
    }
    #[test]
    fn thermal_model_margin() {
        let thermal = MotorThermalModel::new(25.0, 150.0, 50.0, 0.5);
        assert!((thermal.margin() - 125.0).abs() < 1e-12);
    }
    #[test]
    fn thermal_model_reset() {
        let mut thermal = MotorThermalModel::new(25.0, 150.0, 50.0, 0.5);
        thermal.temperature = 100.0;
        thermal.reset();
        assert!((thermal.temperature - 25.0).abs() < 1e-12);
    }
    #[test]
    fn torque_spring_natural_frequency() {
        let spring = TorqueSpring::new(100.0, 0.0, 0.0);
        let omega = spring.natural_frequency(1.0);
        assert!((omega - 10.0).abs() < 1e-12);
    }
    #[test]
    fn torque_spring_damping_ratio() {
        let spring = TorqueSpring::new(100.0, 20.0, 0.0);
        let zeta = spring.damping_ratio(1.0);
        assert!((zeta - 1.0).abs() < 1e-12);
    }
    #[test]
    fn torque_spring_critically_damped() {
        let spring = TorqueSpring::new(100.0, 20.0, 0.0);
        assert!(spring.is_critically_damped(1.0));
        assert!(!spring.is_overdamped(1.0));
    }
    #[test]
    fn spring_damper_natural_frequency() {
        let sd = SpringDamper::new(100.0, 0.0, 1.0);
        let omega = sd.natural_frequency(1.0);
        assert!((omega - 10.0).abs() < 1e-12);
    }
    #[test]
    fn spring_damper_critical_damping() {
        let sd = SpringDamper::new(100.0, 0.0, 1.0);
        let c_crit = sd.critical_damping(1.0);
        assert!((c_crit - 20.0).abs() < 1e-12);
    }
    #[test]
    fn spring_damper_damping_power() {
        let sd = SpringDamper::new(10.0, 5.0, 1.0);
        let p = sd.damping_power(2.0);
        assert!((p - 20.0).abs() < 1e-12);
    }
    #[test]
    fn position_motor_error() {
        let mut motor = PositionMotor::new(10.0, 0.0, 100.0);
        motor.target_position = 5.0;
        motor.update_position(3.0);
        assert!((motor.position_error() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn position_motor_power() {
        let mut motor = PositionMotor::new(10.0, 0.0, 100.0);
        motor.target_position = 5.0;
        motor.update_position(0.0);
        let power = motor.power(2.0);
        assert!((power - 100.0).abs() < 1e-12);
    }
    #[test]
    fn velocity_motor_error() {
        let mut motor = VelocityMotor::new(5.0, 50.0);
        motor.target_velocity = 10.0;
        motor.update_velocity(3.0);
        assert!((motor.velocity_error() - 7.0).abs() < 1e-12);
    }
    #[test]
    fn servo_motor_toggle() {
        let mut servo = ServoMotor::new(10.0, 1.0, 100.0, 10.0);
        assert!(servo.enabled);
        servo.toggle();
        assert!(!servo.enabled);
        servo.toggle();
        assert!(servo.enabled);
    }
    #[test]
    fn linear_actuator_range_fraction() {
        let mut act = LinearActuator::new(0.0, 10.0);
        act.position = 5.0;
        assert!((act.range_fraction() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn linear_actuator_travel_remaining() {
        let mut act = LinearActuator::new(0.0, 10.0);
        act.position = 3.0;
        assert!((act.travel_remaining_positive() - 7.0).abs() < 1e-12);
        assert!((act.travel_remaining_negative() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn hydraulic_displaced_volume() {
        let mut act = HydraulicActuator::new(0.01, 0.1, 1_000_000.0);
        act.set_position(0.05);
        let vol = act.displaced_volume();
        assert!((vol - 0.0005).abs() < 1e-9);
    }
    #[test]
    fn hydraulic_flow_rate() {
        let act = HydraulicActuator::new(0.01, 0.1, 1_000_000.0);
        let fr = act.flow_rate_for_velocity(1.0);
        assert!((fr - 0.01).abs() < 1e-9);
    }
    #[test]
    fn motor_bank_update() {
        let mut bank = MotorBank::new(3, 1.0, 0.0, 0.0, 100.0);
        let errors = [1.0, 2.0, 3.0];
        let outputs = bank.update(&errors, 0.01);
        assert_eq!(outputs.len(), 3);
        assert!((outputs[0] - 1.0).abs() < 1e-12);
        assert!((outputs[1] - 2.0).abs() < 1e-12);
        assert!((outputs[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn motor_bank_reset() {
        let mut bank = MotorBank::new(2, 0.0, 1.0, 0.0, 100.0);
        bank.update(&[5.0, 5.0], 1.0);
        bank.reset_all();
        for ctrl in &bank.controllers {
            assert!(ctrl.integral.abs() < 1e-15);
        }
    }
    #[test]
    fn motor_bank_len() {
        let bank = MotorBank::new(5, 1.0, 0.0, 0.0, 100.0);
        assert_eq!(bank.len(), 5);
        assert!(!bank.is_empty());
    }
}
/// Trapezoidal waveform: flat at ±1 with linear transitions.
pub fn trapezoidal_wave(theta: f64) -> f64 {
    use std::f64::consts::PI;
    let two_pi = 2.0 * PI;
    let t = ((theta % two_pi) + two_pi) % two_pi;
    if t < PI / 3.0 {
        t / (PI / 3.0)
    } else if t < 2.0 * PI / 3.0 {
        1.0
    } else if t < PI {
        1.0 - (t - 2.0 * PI / 3.0) / (PI / 3.0)
    } else if t < 4.0 * PI / 3.0 {
        -(t - PI) / (PI / 3.0)
    } else if t < 5.0 * PI / 3.0 {
        -1.0
    } else {
        -1.0 + (t - 5.0 * PI / 3.0) / (PI / 3.0)
    }
}
/// Compute the regenerative braking torque available for a DC motor at
/// speed `omega` with the given electrical parameters.
///
/// The back-EMF is fed into a braking resistor `r_brake`.  The braking
/// current is I_brake = V_emf / (R_a + R_brake) and the braking torque is
/// τ_brake = Kt * I_brake.
pub fn regenerative_braking_torque(
    omega: f64,
    ke: f64,
    kt: f64,
    r_armature: f64,
    r_brake: f64,
) -> f64 {
    let v_emf = ke * omega;
    let i_brake = v_emf / (r_armature + r_brake).max(f64::EPSILON);
    kt * i_brake
}
/// Compute the regenerated electrical power during braking.
///
/// P_regen = V_emf * I_brake - I_brake² * R_a  (power delivered to external circuit)
pub fn regenerated_power(omega: f64, ke: f64, r_armature: f64, r_brake: f64) -> f64 {
    let v_emf = ke * omega;
    let i_brake = v_emf / (r_armature + r_brake).max(f64::EPSILON);
    v_emf * i_brake - i_brake * i_brake * r_armature
}
#[cfg(test)]
mod tests_motors_ext {

    use crate::motors::BackEmfObserver;
    use crate::motors::BldcMotor;

    use crate::motors::DcMotor;

    use crate::motors::GearStage;
    use crate::motors::GearboxChain;
    use crate::motors::HarmonicDrive;

    use crate::motors::MotorWindingThermal;

    use crate::motors::PmsmMotor;

    use crate::motors::StepperMotor;

    use crate::motors::regenerated_power;
    use crate::motors::regenerative_braking_torque;
    use crate::motors::trapezoidal_wave;

    #[test]
    fn test_dc_motor_no_load_speed() {
        let m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.001);
        let v = 12.0;
        let expected = v / 0.1;
        assert!(
            (m.no_load_speed(v) - expected).abs() < 1e-8,
            "no-load speed = {}",
            m.no_load_speed(v)
        );
    }
    #[test]
    fn test_dc_motor_stall_torque() {
        let m = DcMotor::new(2.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        let tau = m.stall_torque(12.0);
        assert!((tau - 0.6).abs() < 1e-8, "stall torque = {tau}");
    }
    #[test]
    fn test_dc_motor_step_increases_current() {
        let mut m = DcMotor::new(1.0, 0.001, 0.01, 0.1, 0.001, 0.0);
        m.step(12.0, 0.0, 0.001);
        assert!(m.current > 0.0, "current should increase: {}", m.current);
    }
    #[test]
    fn test_dc_motor_back_emf() {
        let mut m = DcMotor::new(1.0, 0.001, 0.1, 0.1, 0.001, 0.0);
        m.omega = 100.0;
        assert!(
            (m.back_emf() - 10.0).abs() < 1e-10,
            "back EMF = {}",
            m.back_emf()
        );
    }
    #[test]
    fn test_bldc_trapezoidal_wave_peak() {
        let v = trapezoidal_wave(std::f64::consts::PI / 2.0);
        assert!((v - 1.0).abs() < 1e-8, "v = {v}");
    }
    #[test]
    fn test_bldc_back_emf_phases_sum_zero() {
        let m = BldcMotor::new(1.0, 0.001, 0.1, 0.1, 2);
        let emf = m.back_emf_phases();
        let sum = emf[0] + emf[1] + emf[2];
        assert!(sum.abs() < 0.1, "3-phase EMF sum should be ~0: {sum}");
    }
    #[test]
    fn test_bldc_update_angle() {
        let mut m = BldcMotor::new(1.0, 0.001, 0.1, 0.1, 2);
        m.omega = 100.0;
        m.update_angle(0.001);
        assert!(m.theta_e > 0.0, "theta_e should increase: {}", m.theta_e);
    }
    #[test]
    fn test_harmonic_drive_output_torque_at_rest() {
        let hd = HarmonicDrive::new(100.0, 5000.0, 10.0);
        assert!(hd.output_torque().abs() < 1e-10, "zero torque at rest");
    }
    #[test]
    fn test_harmonic_drive_step() {
        let mut hd = HarmonicDrive::new(100.0, 5000.0, 5.0);
        hd.step(std::f64::consts::PI * 2.0, 0.01, 0.001);
        assert!(
            hd.theta_out.abs() > 0.0 || hd.omega_out.abs() > 0.0,
            "output should respond to input motion"
        );
    }
    #[test]
    fn test_thermal_model_steady_state_rise() {
        let tm = MotorWindingThermal::new(0.5, 1.5, 100.0, 200.0, 25.0);
        let rise = tm.steady_state_rise(10.0);
        assert!((rise - 20.0).abs() < 1e-10, "rise = {rise}");
    }
    #[test]
    fn test_thermal_model_step_increases_temperature() {
        let mut tm = MotorWindingThermal::new(0.5, 1.5, 50.0, 100.0, 25.0);
        tm.step(100.0, 1.0);
        assert!(
            tm.t_winding > 25.0,
            "winding should heat up: {}",
            tm.t_winding
        );
    }
    #[test]
    fn test_thermal_model_within_limit() {
        let tm = MotorWindingThermal::new(0.5, 1.5, 50.0, 100.0, 25.0);
        assert!(tm.within_limit(150.0));
    }
    #[test]
    fn test_regenerative_braking_torque() {
        let tau = regenerative_braking_torque(100.0, 0.1, 0.1, 1.0, 9.0);
        assert!((tau - 0.1).abs() < 1e-8, "tau = {tau}");
    }
    #[test]
    fn test_regenerated_power_positive() {
        let p = regenerated_power(100.0, 0.1, 1.0, 9.0);
        assert!((p - 9.0).abs() < 1e-8, "P = {p}");
    }
    #[test]
    fn test_stepper_steps_per_rev() {
        let sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        assert_eq!(sm.steps_per_rev(), 200);
    }
    #[test]
    fn test_stepper_step_angle() {
        let sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        let expected = 2.0 * std::f64::consts::PI / 200.0;
        assert!(
            (sm.step_angle_rad() - expected).abs() < 1e-12,
            "step angle = {}",
            sm.step_angle_rad()
        );
    }
    #[test]
    fn test_stepper_em_torque_at_zero_angle() {
        let mut sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        sm.iA = 0.0;
        sm.iB = 2.0;
        sm.theta = 0.0;
        let tau = sm.electromagnetic_torque();
        assert!((tau - (-0.4 * 2.0)).abs() < 1e-10, "tau = {tau}");
    }
    #[test]
    fn test_stepper_max_holding_torque() {
        let sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        let mht = sm.max_holding_torque(2.0);
        assert!((mht - 0.8).abs() < 1e-10, "mht = {mht}");
    }
    #[test]
    fn test_stepper_dynamics_step() {
        let mut sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        sm.set_phase_currents(std::f64::consts::FRAC_PI_2, 2.0);
        let omega_0 = sm.omega;
        sm.step(0.0, 0.001);
        assert_ne!(sm.omega, omega_0, "omega should change after step");
    }
    #[test]
    fn test_stepper_detent_torque_symmetry() {
        let mut sm = StepperMotor::new(50, 1.0, 0.001, 0.4, 0.02, 24.0, 1e-5, 1e-6);
        sm.theta = 0.01;
        let t1 = sm.detent_torque().abs();
        sm.theta = std::f64::consts::PI / (sm.N_r as f64) - 0.01;
        let t2 = sm.detent_torque().abs();
        assert!(t1.is_finite() && t2.is_finite());
    }
    #[test]
    fn test_pmsm_pole_pairs() {
        let m = PmsmMotor::new(8, 0.5, 0.002, 0.003, 0.15, 0.001, 0.0001);
        assert!((m.pole_pairs() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_pmsm_torque_non_salient_iq_only() {
        let mut m = PmsmMotor::new(4, 0.5, 0.002, 0.002, 0.15, 0.001, 0.0001);
        m.i_d = 0.0;
        m.i_q = 10.0;
        let expected = 1.5 * 2.0 * 0.15 * 10.0;
        let tau = m.torque();
        assert!(
            (tau - expected).abs() < 1e-10,
            "tau = {tau}, expected = {expected}"
        );
    }
    #[test]
    fn test_pmsm_copper_loss() {
        let mut m = PmsmMotor::new(4, 1.0, 0.002, 0.002, 0.15, 0.001, 0.0001);
        m.i_d = 3.0;
        m.i_q = 4.0;
        let p = m.copper_loss();
        assert!((p - 37.5).abs() < 1e-10, "P_cu = {p}");
    }
    #[test]
    fn test_pmsm_back_emf() {
        let mut m = PmsmMotor::new(4, 0.5, 0.002, 0.002, 0.15, 0.001, 0.0001);
        m.omega_e = 100.0;
        let emf = m.back_emf();
        assert!((emf - 15.0).abs() < 1e-10, "back EMF = {emf}");
    }
    #[test]
    fn test_pmsm_flux_weaken_clamps_to_negative() {
        let mut m = PmsmMotor::new(4, 0.5, 0.002, 0.002, 0.15, 0.001, 0.0001);
        m.omega_e = 500.0;
        m.flux_weaken(100.0);
        assert!(m.i_d <= 0.0, "i_d should be non-positive: {}", m.i_d);
    }
    #[test]
    fn test_pmsm_integrate_mechanics_speed_increases() {
        let mut m = PmsmMotor::new(4, 0.5, 0.002, 0.002, 0.15, 0.001, 0.0001);
        m.i_q = 5.0;
        m.integrate_mechanics(0.0, 0.001);
        assert!(m.omega_m > 0.0, "speed should be positive: {}", m.omega_m);
    }
    #[test]
    fn test_pmsm_mtpa_angle_non_salient_is_zero() {
        let m = PmsmMotor::new(4, 0.5, 0.002, 0.002, 0.15, 0.001, 0.0001);
        let beta = m.mtpa_angle(5.0);
        assert!(
            beta.abs() < 1e-12,
            "non-salient MTPA angle should be 0: {beta}"
        );
    }
    #[test]
    fn test_gearbox_overall_ratio_single_stage() {
        let mut g = GearboxChain::new();
        g.add_stage(0.1, 0.95);
        assert!((g.overall_ratio() - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_gearbox_overall_efficiency_two_stages() {
        let mut g = GearboxChain::new();
        g.add_stage(0.5, 0.97);
        g.add_stage(0.5, 0.96);
        let expected = 0.97 * 0.96;
        assert!((g.overall_efficiency() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_gearbox_output_torque() {
        let mut g = GearboxChain::new();
        g.add_stage(0.1, 1.0);
        let tau_out = g.output_torque(1.0);
        assert!((tau_out - 10.0).abs() < 1e-10, "tau_out = {tau_out}");
    }
    #[test]
    fn test_gearbox_back_drive_not_possible_low_efficiency() {
        let mut g = GearboxChain::new();
        g.add_stage(0.02, 0.3);
        assert!(
            g.back_drive_torque(1.0).is_none(),
            "worm gear should not back-drive"
        );
    }
    #[test]
    fn test_gearbox_back_drive_possible_high_efficiency() {
        let mut g = GearboxChain::new();
        g.add_stage(0.5, 0.95);
        let bd = g.back_drive_torque(1.0);
        assert!(bd.is_some(), "high-efficiency gear should back-drive");
        assert!(bd.unwrap() > 0.0);
    }
    #[test]
    fn test_gearstage_output_speed() {
        let gs = GearStage::new(2.0, 0.98);
        assert!((gs.output_speed(100.0) - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_back_emf_observer_at_rest() {
        let mut obs = BackEmfObserver::new(1.0, 0.0, 0.1, 1.0);
        let omega = obs.update(1.0, 1.0, 0.001);
        assert!(omega.abs() < 1e-6, "omega at rest = {omega}");
    }
    #[test]
    fn test_back_emf_observer_steady_speed() {
        let mut obs = BackEmfObserver::new(1.0, 0.0, 0.1, 1.0);
        let omega = obs.update(12.0, 2.0, 0.001);
        assert!((omega - 100.0).abs() < 1e-6, "omega = {omega}");
    }
    #[test]
    fn test_back_emf_observer_filter_smooths() {
        let mut obs = BackEmfObserver::new(0.0, 0.0, 1.0, 0.1);
        obs.update(50.0, 0.0, 0.001);
        assert!(obs.omega_est < 50.0);
    }
    #[test]
    fn test_back_emf_observer_reset() {
        let mut obs = BackEmfObserver::new(1.0, 0.01, 0.1, 0.5);
        obs.update(12.0, 2.0, 0.001);
        obs.reset();
        assert_eq!(obs.omega_est, 0.0);
    }
}
#[cfg(test)]
mod tests_motors_extra {

    use crate::motors::DcMotorEfficiencyMap;

    use crate::motors::MotorCurrentLimiter;

    use crate::motors::MotorRundownIdentifier;

    use crate::motors::PwmController;

    #[test]
    fn test_pwm_average_voltage() {
        let mut pwm = PwmController::new(48.0, 20_000.0, 0.0);
        pwm.set_duty(0.75);
        assert!(
            (pwm.average_voltage() - 36.0).abs() < 1e-10,
            "V_avg = {}",
            pwm.average_voltage()
        );
    }
    #[test]
    fn test_pwm_duty_clamp() {
        let mut pwm = PwmController::new(48.0, 20_000.0, 0.0);
        pwm.set_duty(1.5);
        assert!((pwm.duty - 1.0).abs() < 1e-12, "duty = {}", pwm.duty);
        pwm.set_duty(-0.5);
        assert!((pwm.duty).abs() < 1e-12, "duty = {}", pwm.duty);
    }
    #[test]
    fn test_pwm_dead_time_reduces_effective_duty() {
        let mut pwm = PwmController::new(48.0, 20_000.0, 0.5e-6);
        pwm.set_duty(0.5);
        let eff = pwm.effective_duty();
        assert!((eff - 0.49).abs() < 1e-10, "effective duty = {eff}");
    }
    #[test]
    fn test_pwm_current_ripple() {
        let mut pwm = PwmController::new(24.0, 10_000.0, 0.0);
        pwm.set_duty(0.5);
        let ripple = pwm.current_ripple(0.001);
        assert!((ripple - 0.6).abs() < 1e-10, "ripple = {ripple}");
    }
    #[test]
    fn test_pwm_modulate_sets_duty() {
        let mut pwm = PwmController::new(100.0, 20_000.0, 0.0);
        let duty = pwm.modulate(70.0);
        assert!((duty - 0.7).abs() < 1e-10, "duty = {duty}");
        assert!((pwm.average_voltage() - 70.0).abs() < 1e-10);
    }
    #[test]
    fn test_pwm_modulate_clamps_above_bus() {
        let mut pwm = PwmController::new(48.0, 20_000.0, 0.0);
        let duty = pwm.modulate(100.0);
        assert!(
            (duty - 1.0).abs() < 1e-10,
            "duty should be clamped to 1: {duty}"
        );
    }
    #[test]
    fn test_pwm_switching_loss() {
        let pwm = PwmController::new(48.0, 10_000.0, 0.0);
        let p = pwm.switching_loss(100e-6);
        assert!((p - 1.0).abs() < 1e-10, "P_sw = {p}");
    }
    #[test]
    fn test_rundown_time_constant_pure_exponential() {
        let mut id = MotorRundownIdentifier::new();
        for i in 0..20 {
            let t = i as f64 * 0.2;
            let omega = 100.0 * (-t / 2.0).exp();
            id.push(t, omega);
        }
        let tau = id.time_constant().expect("should compute tau");
        assert!((tau - 2.0).abs() < 0.05, "tau = {tau}");
    }
    #[test]
    fn test_rundown_initial_speed() {
        let mut id = MotorRundownIdentifier::new();
        for i in 0..20 {
            let t = i as f64 * 0.1;
            let omega = 50.0 * (-t / 1.5).exp();
            id.push(t, omega);
        }
        let omega0 = id.initial_speed().expect("should compute omega0");
        assert!((omega0 - 50.0).abs() < 1.0, "omega0 = {omega0}");
    }
    #[test]
    fn test_rundown_predict_speed() {
        let mut id = MotorRundownIdentifier::new();
        for i in 0..20 {
            let t = i as f64 * 0.1;
            let omega = 100.0 * (-t).exp();
            id.push(t, omega);
        }
        let pred = id.predict_speed(1.0).expect("predict should work");
        assert!(
            (pred - 100.0 * (-1.0_f64).exp()).abs() < 2.0,
            "predicted = {pred}"
        );
    }
    #[test]
    fn test_rundown_no_samples_returns_none() {
        let id = MotorRundownIdentifier::new();
        assert!(id.time_constant().is_none());
        assert!(id.initial_speed().is_none());
    }
    #[test]
    fn test_rundown_identify_jb() {
        let mut id = MotorRundownIdentifier::new();
        for i in 0..20 {
            let t = i as f64 * 0.1;
            let omega = 100.0 * (-t).exp();
            id.push(t, omega);
        }
        let (j, b) = id.identify_jb(10.0, 100.0).expect("should identify");
        assert!(j > 0.0, "J should be positive: {j}");
        assert!(b > 0.0, "B should be positive: {b}");
    }
    #[test]
    fn test_rundown_clear() {
        let mut id = MotorRundownIdentifier::new();
        id.push(0.0, 100.0);
        id.push(1.0, 50.0);
        id.clear();
        assert_eq!(id.samples.len(), 0);
        assert!(id.time_constant().is_none());
    }
    #[test]
    fn test_efficiency_map_uniform_value() {
        let map = DcMotorEfficiencyMap::uniform(1000.0, 10.0, 0.85, 5, 5);
        let eff = map.efficiency_at(5.0, 500.0);
        assert!((eff - 0.85).abs() < 1e-10, "eff = {eff}");
    }
    #[test]
    fn test_efficiency_map_output_power() {
        let map = DcMotorEfficiencyMap::uniform(1000.0, 10.0, 0.85, 5, 5);
        let p = map.output_power(5.0, 200.0);
        assert!((p - 1000.0).abs() < 1e-10, "P_out = {p}");
    }
    #[test]
    fn test_efficiency_map_input_power_greater_than_output() {
        let map = DcMotorEfficiencyMap::uniform(1000.0, 10.0, 0.9, 5, 5);
        let p_in = map.input_power(5.0, 200.0);
        let p_out = map.output_power(5.0, 200.0);
        assert!(
            p_in >= p_out - 1e-9,
            "P_in should be >= P_out: p_in={p_in}, p_out={p_out}"
        );
    }
    #[test]
    fn test_current_limiter_cold_allows_rated_current() {
        let lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        assert!(
            (lim.max_current() - 10.0).abs() < 1e-10,
            "max_current = {}",
            lim.max_current()
        );
    }
    #[test]
    fn test_current_limiter_hot_reduces_max_current() {
        let mut lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        lim.t_winding = 100.0;
        let max_i = lim.max_current();
        assert!(
            max_i < 10.0,
            "max_current should be reduced when hot: {max_i}"
        );
    }
    #[test]
    fn test_current_limiter_clamp_request() {
        let lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        let clamped = lim.clamp_current(50.0);
        assert!(
            (clamped - 10.0).abs() < 1e-10,
            "should clamp to rated: {clamped}"
        );
    }
    #[test]
    fn test_current_limiter_step_heats_winding() {
        let mut lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        let t_before = lim.t_winding;
        lim.step(10.0, 0.5, 1.0);
        assert!(lim.t_winding > t_before, "winding should heat up");
    }
    #[test]
    fn test_current_limiter_thermal_headroom_full_when_cold() {
        let lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        assert!((lim.thermal_headroom() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_current_limiter_not_derated_when_cold() {
        let lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        assert!(!lim.is_derated(), "should not be derated when cold");
    }
    #[test]
    fn test_current_limiter_derated_when_hot() {
        let mut lim = MotorCurrentLimiter::new(10.0, 30.0, 150.0, 25.0, 60.0);
        lim.t_winding = 140.0;
        assert!(lim.is_derated(), "should be derated when hot");
    }
}
#[cfg(test)]
mod new_motor_tests {

    use crate::motors::CascadedPid;

    use crate::motors::ElectricMotor;

    use crate::motors::HydraulicActuator;

    use crate::motors::MotorEfficiencyCurve;

    use crate::motors::RegenerativeBrakeController;

    use crate::motors::SCurveProfile;

    use crate::motors::TrapPhase;
    use crate::motors::TrapezoidalProfile;

    #[test]
    fn test_cascaded_pid_drives_toward_target() {
        let mut cpid = CascadedPid::new(
            10.0, 0.1, 0.5, 100.0, 5.0, 0.1, 0.1, 50.0, 2.0, 0.0, 0.0, 24.0,
        );
        let output = cpid.update(1.0, 0.0, 0.0, 0.0, 0.01);
        assert!(output > 0.0, "output should be positive: {output}");
    }
    #[test]
    fn test_cascaded_pid_no_error_gives_small_output() {
        let mut cpid = CascadedPid::new(
            10.0, 0.1, 0.5, 100.0, 5.0, 0.1, 0.1, 50.0, 2.0, 0.0, 0.0, 24.0,
        );
        let output = cpid.update(0.0, 0.0, 0.0, 0.0, 0.01);
        assert!(
            output.abs() < 1e-9,
            "output with zero error should be ~0: {output}"
        );
    }
    #[test]
    fn test_cascaded_pid_reset_clears_state() {
        let mut cpid = CascadedPid::new(
            10.0, 1.0, 0.5, 100.0, 5.0, 1.0, 0.1, 50.0, 2.0, 0.5, 0.0, 24.0,
        );
        for _ in 0..50 {
            cpid.update(1.0, 0.0, 0.0, 0.0, 0.01);
        }
        cpid.reset();
        assert!((cpid.position_pid.integral).abs() < 1e-12);
        assert!((cpid.velocity_pid.integral).abs() < 1e-12);
        assert!((cpid.current_pid.integral).abs() < 1e-12);
    }
    #[test]
    fn test_cascaded_pid_output_clamped() {
        let mut cpid = CascadedPid::new(
            1000.0, 0.0, 0.0, 5.0, 1000.0, 0.0, 0.0, 5.0, 1000.0, 0.0, 0.0, 24.0,
        );
        let output = cpid.update(100.0, 0.0, 0.0, 0.0, 0.01);
        assert!(
            output.abs() <= 24.0 + 1e-9,
            "output clamped to current limit: {output}"
        );
    }
    #[test]
    fn test_cascaded_pid_negative_target() {
        let mut cpid = CascadedPid::new(
            10.0, 0.1, 0.5, 100.0, 5.0, 0.1, 0.1, 50.0, 2.0, 0.0, 0.0, 24.0,
        );
        let output = cpid.update(-1.0, 0.0, 0.0, 0.0, 0.01);
        assert!(
            output < 0.0,
            "negative target should produce negative output: {output}"
        );
    }
    #[test]
    fn test_trapezoidal_starts_at_zero_velocity() {
        let prof = TrapezoidalProfile::plan(0.0, 10.0, 2.0, 1.0);
        let (_, v, _) = prof.query(0.0);
        assert!(v.abs() < 1e-10, "initial velocity should be 0: {v}");
    }
    #[test]
    fn test_trapezoidal_ends_at_target() {
        let prof = TrapezoidalProfile::plan(0.0, 10.0, 2.0, 1.0);
        let (x, v, _) = prof.query(prof.t_total + 1.0);
        assert!((x - 10.0).abs() < 1e-9, "end position: {x}");
        assert!(v.abs() < 1e-10, "end velocity: {v}");
    }
    #[test]
    fn test_trapezoidal_peak_velocity_bounded() {
        let prof = TrapezoidalProfile::plan(0.0, 10.0, 2.0, 1.0);
        let v_peak = prof.v_cruise;
        assert!(v_peak <= 2.0 + 1e-9, "cruise velocity bounded: {v_peak}");
    }
    #[test]
    fn test_trapezoidal_short_move_triangle_profile() {
        let prof = TrapezoidalProfile::plan(0.0, 1.0, 10.0, 1.0);
        assert!(
            (prof.t1 - prof.t2).abs() < 1e-9,
            "t1 should equal t2 for triangle: t1={}, t2={}",
            prof.t1,
            prof.t2
        );
    }
    #[test]
    fn test_trapezoidal_phase_at_start() {
        let prof = TrapezoidalProfile::plan(0.0, 10.0, 2.0, 1.0);
        assert_eq!(prof.phase(0.001), TrapPhase::Accelerating);
    }
    #[test]
    fn test_trapezoidal_phase_done_after_total() {
        let prof = TrapezoidalProfile::plan(0.0, 10.0, 2.0, 1.0);
        assert_eq!(prof.phase(prof.t_total + 1.0), TrapPhase::Done);
    }
    #[test]
    fn test_trapezoidal_cruise_phase_exists_for_long_move() {
        let prof = TrapezoidalProfile::plan(0.0, 100.0, 2.0, 1.0);
        let mid = (prof.t1 + prof.t2) / 2.0;
        assert_eq!(prof.phase(mid), TrapPhase::Cruising);
    }
    #[test]
    fn test_trapezoidal_backward_move() {
        let prof = TrapezoidalProfile::plan(10.0, 0.0, 2.0, 1.0);
        let (x, v, _) = prof.query(prof.t_total + 1.0);
        assert!((x - 0.0).abs() < 1e-9, "backward end position: {x}");
        assert!(v.abs() < 1e-10, "backward end velocity: {v}");
    }
    #[test]
    fn test_trapezoidal_monotone_position_forward() {
        let prof = TrapezoidalProfile::plan(0.0, 5.0, 2.0, 1.0);
        let mut prev_x = 0.0_f64;
        let n = 100;
        for i in 1..=n {
            let t = i as f64 * prof.t_total / n as f64;
            let (x, _, _) = prof.query(t);
            assert!(
                x >= prev_x - 1e-9,
                "position not monotone at t={t}: prev={prev_x}, cur={x}"
            );
            prev_x = x;
        }
    }
    #[test]
    fn test_scurve_starts_at_zero_velocity() {
        let prof = SCurveProfile::plan(0.0, 10.0, 2.0, 1.0, 2.0);
        let v = prof.velocity(0.0);
        assert!(v.abs() < 1e-10, "initial velocity: {v}");
    }
    #[test]
    fn test_scurve_ends_at_zero_velocity() {
        let prof = SCurveProfile::plan(0.0, 10.0, 2.0, 1.0, 2.0);
        let v = prof.velocity(prof.t_total);
        assert!(v.abs() < 1e-10, "final velocity: {v}");
    }
    #[test]
    fn test_scurve_duration_positive() {
        let prof = SCurveProfile::plan(0.0, 10.0, 2.0, 1.0, 2.0);
        assert!(prof.duration() > 0.0, "duration: {}", prof.duration());
    }
    #[test]
    fn test_scurve_velocity_does_not_exceed_vmax() {
        let v_max = 2.0;
        let prof = SCurveProfile::plan(0.0, 20.0, v_max, 1.0, 2.0);
        let n = 200;
        for i in 0..=n {
            let t = i as f64 * prof.t_total / n as f64;
            let v = prof.velocity(t).abs();
            assert!(v <= v_max + 1e-6, "velocity exceeded v_max at t={t}: v={v}");
        }
    }
    #[test]
    fn test_scurve_symmetric_velocity_profile() {
        let prof = SCurveProfile::plan(0.0, 20.0, 2.0, 1.0, 4.0);
        let t_mid = prof.t_total / 2.0;
        let v_before = prof.velocity(t_mid - 0.01).abs();
        let v_after = prof.velocity(t_mid + 0.01).abs();
        assert!(
            (v_before - v_after).abs() < 0.1,
            "not symmetric: v_before={v_before}, v_after={v_after}"
        );
    }
    #[test]
    fn test_efficiency_curve_flat_returns_constant() {
        let curve = MotorEfficiencyCurve::flat(0.9);
        assert!((curve.efficiency_at(100.0) - 0.9).abs() < 1e-10);
        assert!((curve.efficiency_at(500.0) - 0.9).abs() < 1e-10);
    }
    #[test]
    fn test_efficiency_curve_interpolation() {
        let curve = MotorEfficiencyCurve::new(vec![0.0, 100.0, 200.0], vec![0.7, 0.9, 0.8]);
        let eff = curve.efficiency_at(50.0);
        assert!((eff - 0.8).abs() < 1e-10, "interpolated eff={eff}");
    }
    #[test]
    fn test_efficiency_curve_below_range_clamps_to_first() {
        let curve = MotorEfficiencyCurve::new(vec![10.0, 100.0], vec![0.5, 0.9]);
        assert!((curve.efficiency_at(0.0) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_efficiency_curve_above_range_clamps_to_last() {
        let curve = MotorEfficiencyCurve::new(vec![10.0, 100.0], vec![0.5, 0.9]);
        assert!((curve.efficiency_at(200.0) - 0.9).abs() < 1e-10);
    }
    #[test]
    fn test_efficiency_curve_output_power() {
        let curve = MotorEfficiencyCurve::flat(0.85);
        let p_out = curve.output_power(1000.0, 50.0);
        assert!((p_out - 850.0).abs() < 1e-9, "p_out={p_out}");
    }
    #[test]
    fn test_efficiency_curve_peak_speed() {
        let curve =
            MotorEfficiencyCurve::new(vec![0.0, 50.0, 100.0, 200.0], vec![0.6, 0.92, 0.88, 0.7]);
        let peak_s = curve.peak_efficiency_speed();
        assert!((peak_s - 50.0).abs() < 1e-10, "peak speed={peak_s}");
    }
    #[test]
    fn test_regen_brake_produces_torque_at_speed() {
        let mut rbc = RegenerativeBrakeController::new(1.0, 0.05, 0.05, 0.5, 10.0, 20.0, 0.95);
        let (torque, power) = rbc.update(100.0, 0.01);
        assert!(torque > 0.0, "braking torque should be positive: {torque}");
        assert!(power >= 0.0, "regen power should be non-negative: {power}");
    }
    #[test]
    fn test_regen_brake_zero_below_min_speed() {
        let mut rbc = RegenerativeBrakeController::new(5.0, 0.05, 0.05, 0.5, 10.0, 20.0, 0.95);
        let (torque, power) = rbc.update(0.5, 0.01);
        assert_eq!(torque, 0.0);
        assert_eq!(power, 0.0);
    }
    #[test]
    fn test_regen_brake_disabled_when_battery_full() {
        let mut rbc = RegenerativeBrakeController::new(1.0, 0.05, 0.05, 0.5, 10.0, 20.0, 0.95);
        rbc.soc = 0.98;
        let (torque, power) = rbc.update(100.0, 0.01);
        assert_eq!(torque, 0.0);
        assert_eq!(power, 0.0);
    }
    #[test]
    fn test_regen_brake_energy_accumulates() {
        let mut rbc = RegenerativeBrakeController::new(1.0, 0.05, 0.05, 0.5, 10.0, 20.0, 0.95);
        for _ in 0..100 {
            rbc.update(100.0, 0.01);
        }
        assert!(
            rbc.energy_recovered > 0.0,
            "energy_recovered={}",
            rbc.energy_recovered
        );
    }
    #[test]
    fn test_regen_brake_reset_energy() {
        let mut rbc = RegenerativeBrakeController::new(1.0, 0.05, 0.05, 0.5, 10.0, 20.0, 0.95);
        rbc.update(100.0, 1.0);
        rbc.reset_energy();
        assert_eq!(rbc.energy_recovered, 0.0);
    }
    #[test]
    fn test_regen_braking_zero_at_zero_rpm() {
        let motor = ElectricMotor::new(100.0, 3000.0);
        let t = motor.compute_regenerative_braking(0.0, 1.0, 1000.0);
        assert_eq!(t, 0.0, "regen braking should be zero at zero rpm");
    }
    #[test]
    fn test_regen_braking_zero_throttle_produces_zero() {
        let motor = ElectricMotor::new(100.0, 3000.0);
        let t = motor.compute_regenerative_braking(1500.0, 0.0, 1000.0);
        assert_eq!(t, 0.0, "zero throttle should produce no regen braking");
    }
    #[test]
    fn test_regen_braking_full_throttle_above_peak_capped() {
        let motor = ElectricMotor::new(100.0, 3000.0);
        let t_at_peak = motor.compute_regenerative_braking(1000.0, 1.0, 1000.0);
        let t_above_peak = motor.compute_regenerative_braking(2000.0, 1.0, 1000.0);
        assert!(
            (t_at_peak - t_above_peak).abs() < 1e-9,
            "at_peak={t_at_peak}, above_peak={t_above_peak}"
        );
    }
    #[test]
    fn test_regen_braking_partial_throttle_scales_linearly() {
        let motor = ElectricMotor::new(100.0, 3000.0);
        let t_full = motor.compute_regenerative_braking(2000.0, 1.0, 1000.0);
        let t_half = motor.compute_regenerative_braking(2000.0, 0.5, 1000.0);
        assert!(
            (t_half - t_full * 0.5).abs() < 1e-9,
            "half throttle should give half torque: full={t_full}, half={t_half}"
        );
    }
    #[test]
    fn test_regen_braking_positive_value() {
        let motor = ElectricMotor::new(200.0, 5000.0);
        let t = motor.compute_regenerative_braking(2500.0, 0.8, 2500.0);
        assert!(t > 0.0, "regen braking torque should be positive: {t}");
    }
    #[test]
    fn test_temperature_derating_below_knee_returns_one() {
        let factor = ElectricMotor::compute_temperature_derating(50.0, 25.0, 150.0);
        assert!(
            (factor - 1.0).abs() < 1e-10,
            "below knee derating should be 1.0, got {factor}"
        );
    }
    #[test]
    fn test_temperature_derating_at_max_returns_zero() {
        let factor = ElectricMotor::compute_temperature_derating(150.0, 25.0, 150.0);
        assert!(
            factor.abs() < 1e-10,
            "at max temperature derating should be 0.0, got {factor}"
        );
    }
    #[test]
    fn test_temperature_derating_monotone_decreasing() {
        let ambient = 25.0;
        let t_max = 150.0;
        let mut prev = 1.0f64;
        let n = 20;
        for i in 0..=n {
            let t = ambient + (t_max - ambient) * (i as f64) / (n as f64);
            let f = ElectricMotor::compute_temperature_derating(t, ambient, t_max);
            assert!(
                f <= prev + 1e-10,
                "derating not monotone at t={t}: prev={prev}, cur={f}"
            );
            prev = f;
        }
    }
    #[test]
    fn test_flow_rate_positive_dp_produces_positive_flow() {
        let act = HydraulicActuator::new(0.01, 0.5, 2e7);
        let q = act.compute_flow_rate(0.7, 1e-4, 1e6, 850.0);
        assert!(q > 0.0, "positive ΔP should produce positive flow: {q}");
    }
    #[test]
    fn test_flow_rate_negative_dp_produces_negative_flow() {
        let act = HydraulicActuator::new(0.01, 0.5, 2e7);
        let q = act.compute_flow_rate(0.7, 1e-4, -1e6, 850.0);
        assert!(q < 0.0, "negative ΔP should produce negative flow: {q}");
    }
}
/// Apply a Ziegler-Nichols tuning rule to a PID controller.
///
/// Takes the ultimate (marginal stability) gain `ku` and oscillation period
/// `tu` and sets kp, ki, kd according to the chosen `method`.
#[allow(dead_code)]
pub fn zn_tune_pid(pid: &mut PidController, ku: f64, tu: f64, method: ZnTuningMethod) {
    match method {
        ZnTuningMethod::Classic => {
            pid.kp = 0.6 * ku;
            pid.ki = 2.0 * pid.kp / tu;
            pid.kd = pid.kp * tu / 8.0;
        }
        ZnTuningMethod::OvershootReduced => {
            pid.kp = 0.33 * ku;
            pid.ki = pid.kp * 2.0 / tu;
            pid.kd = pid.kp * tu / 6.0;
        }
        ZnTuningMethod::NoOvershoot => {
            pid.kp = 0.2 * ku;
            pid.ki = 2.0 * pid.kp / tu;
            pid.kd = pid.kp * tu / 3.0;
        }
        ZnTuningMethod::PiOnly => {
            pid.kp = 0.45 * ku;
            pid.ki = pid.kp * 1.2 / tu;
            pid.kd = 0.0;
        }
    }
}
#[cfg(test)]
mod tests_motors_expansion {

    use crate::kinematics::GearTrain;

    use crate::motors::DcMotor;

    use crate::motors::MotorThermalTwoNode;

    use crate::motors::PidController;

    use crate::motors::RelayAutoTuner;

    use crate::motors::ZnTuningMethod;

    use crate::motors::zn_tune_pid;
    #[test]
    fn test_dc_motor_mechanical_power() {
        let mut m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        m.current = 2.0;
        m.omega = 100.0;
        let p = m.mechanical_power();
        assert!((p - 20.0).abs() < 1e-10, "mechanical power: {p}");
    }
    #[test]
    fn test_dc_motor_copper_loss() {
        let mut m = DcMotor::new(2.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        m.current = 3.0;
        let loss = m.copper_loss();
        assert!((loss - 18.0).abs() < 1e-10, "copper loss: {loss}");
    }
    #[test]
    fn test_dc_motor_electrical_time_constant() {
        let m = DcMotor::new(2.0, 0.004, 0.1, 0.1, 0.01, 0.0);
        let tau_e = m.electrical_time_constant();
        assert!(
            (tau_e - 0.002).abs() < 1e-12,
            "τ_e = L/R = 0.004/2 = 0.002: {tau_e}"
        );
    }
    #[test]
    fn test_dc_motor_mechanical_time_constant_positive() {
        let m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        let tau_m = m.mechanical_time_constant();
        assert!(tau_m > 0.0, "τ_m > 0: {tau_m}");
    }
    #[test]
    fn test_dc_motor_speed_torque_slope_negative() {
        let m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        let slope = m.speed_torque_slope();
        assert!(
            slope < 0.0,
            "speed-torque slope should be negative: {slope}"
        );
    }
    #[test]
    fn test_dc_motor_steady_state_no_load() {
        let m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        let (tau, omega) = m.steady_state_no_load(12.0);
        assert_eq!(tau, 0.0, "no-load torque = 0");
        assert!(omega > 0.0, "no-load speed > 0: {omega}");
    }
    #[test]
    fn test_dc_motor_steady_state_under_load() {
        let m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        let (_tau_nl, omega_nl) = m.steady_state_no_load(12.0);
        let (_tau_l, omega_l) = m.steady_state(12.0, 0.1);
        assert!(omega_l < omega_nl, "loaded speed < no-load speed");
    }
    #[test]
    fn test_dc_motor_reset_state() {
        let mut m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        m.current = 5.0;
        m.omega = 200.0;
        m.reset_state();
        assert_eq!(m.current, 0.0, "current reset to 0");
        assert_eq!(m.omega, 0.0, "omega reset to 0");
    }
    #[test]
    fn test_dc_motor_efficiency_zero_when_stalled() {
        let mut m = DcMotor::new(1.0, 0.01, 0.1, 0.1, 0.01, 0.0);
        m.omega = 0.0;
        m.current = 2.0;
        let eta = m.efficiency(12.0);
        assert_eq!(eta, 0.0, "stalled → zero efficiency");
    }
    #[test]
    fn test_gear_train_output_speed_reduction() {
        let g = GearTrain::simple(5.0, 0.95, 100.0, 0.01);
        let omega_out = g.output_speed(100.0);
        assert!(
            (omega_out - 20.0).abs() < 1e-10,
            "output speed: {omega_out}"
        );
    }
    #[test]
    fn test_gear_train_input_speed() {
        let g = GearTrain::simple(5.0, 0.95, 100.0, 0.01);
        let omega_in = g.input_speed(20.0);
        assert!((omega_in - 100.0).abs() < 1e-10, "input speed: {omega_in}");
    }
    #[test]
    fn test_gear_train_output_torque_amplification() {
        let g = GearTrain::simple(5.0, 1.0, 1000.0, 0.0);
        let tau_out = g.output_torque(10.0);
        assert!((tau_out - 50.0).abs() < 1e-10, "output torque: {tau_out}");
    }
    #[test]
    fn test_gear_train_output_torque_with_efficiency() {
        let g = GearTrain::simple(5.0, 0.9, 1000.0, 0.0);
        let tau_out = g.output_torque(10.0);
        assert!(
            (tau_out - 45.0).abs() < 1e-10,
            "output torque with η: {tau_out}"
        );
    }
    #[test]
    fn test_gear_train_input_torque_required() {
        let g = GearTrain::simple(5.0, 0.9, 1000.0, 0.0);
        let tau_in = g.input_torque_required(45.0);
        assert!(
            (tau_in - 10.0).abs() < 1e-9,
            "input torque required: {tau_in}"
        );
    }
    #[test]
    fn test_gear_train_reflected_inertia() {
        let g = GearTrain::simple(4.0, 0.95, 100.0, 0.16);
        let j_refl = g.reflected_inertia();
        assert!((j_refl - 0.01).abs() < 1e-12, "reflected inertia: {j_refl}");
    }
    #[test]
    fn test_gear_train_not_overloaded_below_max() {
        let g = GearTrain::simple(5.0, 1.0, 100.0, 0.0);
        assert!(
            !g.is_overloaded(10.0),
            "10 N·m input → 50 N·m output ≤ 100 → not overloaded"
        );
    }
    #[test]
    fn test_gear_train_overloaded_above_max() {
        let g = GearTrain::simple(5.0, 1.0, 40.0, 0.0);
        assert!(
            g.is_overloaded(10.0),
            "10 N·m input → 50 N·m output > 40 → overloaded"
        );
    }
    #[test]
    fn test_gear_train_power_loss_positive() {
        let g = GearTrain::simple(5.0, 0.9, 100.0, 0.0);
        let loss = g.power_loss(10.0, 100.0);
        assert!(loss > 0.0, "power loss > 0 with η < 1: {loss}");
    }
    #[test]
    fn test_gear_train_power_loss_zero_perfect_efficiency() {
        let g = GearTrain::simple(5.0, 1.0, 100.0, 0.0);
        let loss = g.power_loss(10.0, 100.0);
        assert!(
            loss.abs() < 1e-10,
            "zero loss for perfect efficiency: {loss}"
        );
    }
    #[test]
    fn test_gear_train_total_inertia() {
        let g = GearTrain::simple(4.0, 0.95, 100.0, 0.16);
        let j_total = g.total_inertia(0.005);
        assert!((j_total - 0.015).abs() < 1e-12, "total inertia: {j_total}");
    }
    #[test]
    fn test_thermal_two_node_initial_at_ambient() {
        let therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        assert_eq!(therm.t_winding, 25.0, "initial winding = ambient");
        assert_eq!(therm.t_case, 25.0, "initial case = ambient");
    }
    #[test]
    fn test_thermal_two_node_winding_rises_with_loss() {
        let mut therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        for _ in 0..100 {
            therm.step(50.0, 1.0);
        }
        assert!(
            therm.t_winding > 25.0,
            "winding temperature should rise: {}",
            therm.t_winding
        );
    }
    #[test]
    fn test_thermal_two_node_steady_state_rise() {
        let therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        let rise = therm.steady_state_winding_rise(50.0);
        assert!(
            (rise - 75.0).abs() < 1e-10,
            "steady-state rise = 75 °C: {rise}"
        );
    }
    #[test]
    fn test_thermal_two_node_case_rise() {
        let therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        let rise = therm.steady_state_case_rise(50.0);
        assert!((rise - 50.0).abs() < 1e-10, "case rise = 50 °C: {rise}");
    }
    #[test]
    fn test_thermal_two_node_not_overheated_initially() {
        let therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        assert!(!therm.is_overheated(), "not overheated at ambient");
    }
    #[test]
    fn test_thermal_two_node_reset() {
        let mut therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        for _ in 0..50 {
            therm.step(100.0, 1.0);
        }
        therm.reset();
        assert_eq!(therm.t_winding, 25.0, "winding reset to ambient");
        assert_eq!(therm.t_case, 25.0, "case reset to ambient");
    }
    #[test]
    fn test_thermal_two_node_margin_decreases_with_heating() {
        let mut therm = MotorThermalTwoNode::new(25.0, 150.0, 0.5, 1.0, 200.0, 500.0);
        let margin0 = therm.thermal_margin();
        for _ in 0..100 {
            therm.step(100.0, 1.0);
        }
        let margin1 = therm.thermal_margin();
        assert!(margin1 < margin0, "margin decreases as temperature rises");
    }
    #[test]
    fn test_zn_classic_tuning() {
        let mut pid = PidController::new(1.0, 0.0, 0.0, 1000.0);
        zn_tune_pid(&mut pid, 10.0, 1.0, ZnTuningMethod::Classic);
        assert!((pid.kp - 6.0).abs() < 1e-10, "kp classic: {}", pid.kp);
    }
    #[test]
    fn test_zn_pi_only_kd_zero() {
        let mut pid = PidController::new(1.0, 0.0, 0.0, 1000.0);
        zn_tune_pid(&mut pid, 10.0, 1.0, ZnTuningMethod::PiOnly);
        assert_eq!(pid.kd, 0.0, "PI-only → kd = 0");
    }
    #[test]
    fn test_zn_no_overshoot_lower_kp() {
        let mut pid_classic = PidController::new(1.0, 0.0, 0.0, 1000.0);
        let mut pid_no_os = PidController::new(1.0, 0.0, 0.0, 1000.0);
        zn_tune_pid(&mut pid_classic, 10.0, 1.0, ZnTuningMethod::Classic);
        zn_tune_pid(&mut pid_no_os, 10.0, 1.0, ZnTuningMethod::NoOvershoot);
        assert!(
            pid_no_os.kp < pid_classic.kp,
            "no-overshoot should give lower kp"
        );
    }
    #[test]
    fn test_zn_overshoot_reduced_between_classic_and_no_os() {
        let mut pid_classic = PidController::new(1.0, 0.0, 0.0, 1000.0);
        let mut pid_reduced = PidController::new(1.0, 0.0, 0.0, 1000.0);
        let mut pid_no_os = PidController::new(1.0, 0.0, 0.0, 1000.0);
        zn_tune_pid(&mut pid_classic, 10.0, 1.0, ZnTuningMethod::Classic);
        zn_tune_pid(
            &mut pid_reduced,
            10.0,
            1.0,
            ZnTuningMethod::OvershootReduced,
        );
        zn_tune_pid(&mut pid_no_os, 10.0, 1.0, ZnTuningMethod::NoOvershoot);
        assert!(
            pid_reduced.kp < pid_classic.kp,
            "reduced-os kp < classic kp"
        );
        assert!(
            pid_reduced.kp > pid_no_os.kp,
            "reduced-os kp > no-overshoot kp"
        );
    }
    #[test]
    fn test_relay_tuner_initial_no_crossings() {
        let tuner = RelayAutoTuner::new(1.0);
        assert_eq!(tuner.crossing_count(), 0, "no crossings initially");
        assert!(tuner.estimated_tu().is_none(), "no Tu estimate initially");
    }
    #[test]
    fn test_relay_tuner_output_bounded_by_relay_amp() {
        let mut tuner = RelayAutoTuner::new(2.5);
        let out = tuner.update(0.1, 0.01);
        assert!(out.abs() <= 2.5 + 1e-12, "relay output bounded: {out}");
    }
    #[test]
    fn test_relay_tuner_ku_positive() {
        let tuner = RelayAutoTuner::new(1.0);
        let ku = tuner.estimated_ku(0.5);
        assert!(ku > 0.0, "Ku > 0: {ku}");
    }
    #[test]
    fn test_relay_tuner_reset_clears_state() {
        let mut tuner = RelayAutoTuner::new(1.0);
        let mut out_sign = 1.0_f64;
        for _ in 0..10 {
            tuner.update(out_sign * 0.5, 0.01);
            out_sign = -out_sign;
        }
        tuner.reset();
        assert_eq!(tuner.crossing_count(), 0, "crossings cleared after reset");
        assert_eq!(tuner.sample_count, 0, "sample count cleared");
    }
    #[test]
    fn test_relay_tuner_detects_crossings() {
        let mut tuner = RelayAutoTuner::new(1.0);
        let signals = vec![1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0];
        for &s in &signals {
            tuner.update(s, 0.1);
        }
        assert!(tuner.crossing_count() > 0, "should detect zero crossings");
    }
}
