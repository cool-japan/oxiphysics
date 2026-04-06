//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::types::*;

    #[test]
    fn motor_pid_pure_proportional() {
        let mut pid = MotorPid::new(2.0, 0.0, 0.0);
        let out = pid.update(3.0, 0.1);
        assert!((out - 6.0).abs() < 1e-10, "P output should be kp*error = 6");
    }
    #[test]
    fn motor_pid_integral_accumulates() {
        let mut pid = MotorPid::new(0.0, 1.0, 0.0);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        let out = pid.update(1.0, 0.1);
        assert!((pid.integral - 0.3).abs() < 1e-10, "Integral should be 0.3");
        assert!((out - 0.3).abs() < 1e-10, "I output should equal integral");
    }
    #[test]
    fn motor_pid_integral_clamped() {
        let mut pid = MotorPid::new(0.0, 1.0, 0.0).with_integral_clamp(0.15);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        assert!(
            (pid.integral - 0.15).abs() < 1e-10,
            "Integral should be clamped to 0.15"
        );
    }
    #[test]
    fn motor_pid_zero_dt_returns_zero() {
        let mut pid = MotorPid::new(100.0, 100.0, 100.0);
        let out = pid.update(1.0, 0.0);
        assert_eq!(out, 0.0);
    }
    #[test]
    fn motor_pid_reset() {
        let mut pid = MotorPid::new(1.0, 1.0, 1.0);
        pid.update(5.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.prev_error, 0.0);
    }
    #[test]
    fn motor_pid_accessors() {
        let mut pid = MotorPid::new(1.0, 1.0, 0.0);
        pid.update(5.0, 0.1);
        assert!(pid.integral_value() > 0.0);
        assert!((pid.prev_error_value() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn motor_limits_force_clamp() {
        let limits = MotorLimits::force_only(10.0);
        assert_eq!(limits.clamp_force(15.0), 10.0);
        assert_eq!(limits.clamp_force(-15.0), -10.0);
        assert_eq!(limits.clamp_force(5.0), 5.0);
    }
    #[test]
    fn motor_limits_velocity_clamp() {
        let limits = MotorLimits::force_and_velocity(10.0, 5.0);
        assert_eq!(limits.clamp_velocity(10.0), 5.0);
        assert_eq!(limits.clamp_velocity(-10.0), -5.0);
    }
    #[test]
    fn motor_limits_velocity_no_clamp() {
        let limits = MotorLimits::force_only(10.0);
        assert_eq!(limits.clamp_velocity(100.0), 100.0);
    }
    #[test]
    fn motor_limits_position_range() {
        let limits = MotorLimits::force_only(10.0).with_position_range(-1.0, 1.0);
        assert!(limits.position_in_range(0.0));
        assert!(limits.position_in_range(-1.0));
        assert!(limits.position_in_range(1.0));
        assert!(!limits.position_in_range(-1.5));
        assert!(!limits.position_in_range(1.5));
    }
    #[test]
    fn motor_limits_rate_limiting() {
        let mut limits = MotorLimits::force_only(100.0);
        limits.max_force_rate = Some(10.0);
        let result = limits.rate_limit_force(20.0, 0.0, 1.0);
        assert!(
            (result - 10.0).abs() < 1e-10,
            "Rate should be limited to 10 N/s"
        );
    }
    #[test]
    fn motor_friction_zero_velocity() {
        let friction = MotorFriction::coulomb_viscous(1.0, 0.1);
        assert_eq!(friction.friction_force(0.0), 0.0);
    }
    #[test]
    fn motor_friction_coulomb_direction() {
        let friction = MotorFriction::coulomb_viscous(1.0, 0.0);
        let f_pos = friction.friction_force(1.0);
        let f_neg = friction.friction_force(-1.0);
        assert!(f_pos > 0.0, "Positive velocity → positive friction");
        assert!(f_neg < 0.0, "Negative velocity → negative friction");
    }
    #[test]
    fn motor_friction_viscous_proportional() {
        let friction = MotorFriction::coulomb_viscous(0.0, 2.0);
        let f = friction.friction_force(3.0);
        assert!((f - 6.0).abs() < 1e-10);
    }
    #[test]
    fn motor_friction_stribeck() {
        let friction = MotorFriction::stribeck(1.0, 0.1, 2.0, 0.01);
        let f_low = friction.friction_force(0.001).abs();
        let f_high = friction.friction_force(1.0).abs();
        assert!(f_low > 0.0);
        assert!(f_high > 0.0);
    }
    #[test]
    fn motor_friction_with_threshold() {
        let friction = MotorFriction::coulomb_viscous(1.0, 0.1);
        let f = friction.friction_force_with_threshold(0.001, 0.01);
        assert_eq!(f, 0.0, "Below threshold should give zero friction");
    }
    #[test]
    fn linear_motor_velocity_mode_reaches_target() {
        let mut motor = LinearMotorConstraint::velocity_mode(5.0, 1000.0);
        let dt = 0.01;
        let mut vel = 0.0;
        let mut pos = 0.0;
        for _ in 0..500 {
            let force = motor.step(dt, vel, pos);
            vel += force * dt;
            pos += vel * dt;
        }
        assert!(
            (vel - 5.0).abs() < 0.5,
            "Linear motor should reach target velocity 5 m/s, got {vel}"
        );
    }
    #[test]
    fn linear_motor_force_clamped() {
        let mut motor = LinearMotorConstraint::velocity_mode(100.0, 10.0);
        let force = motor.step(0.01, 0.0, 0.0);
        assert!(
            force.abs() <= 10.0,
            "Force must be within max_force limit: {force}"
        );
    }
    #[test]
    fn linear_motor_position_mode_converges() {
        let mut motor = LinearMotorConstraint::position_mode(1.0, 1000.0);
        let dt = 0.01;
        let mut vel = 0.0;
        let mut pos = 0.0;
        for _ in 0..2000 {
            let force = motor.step(dt, vel, pos);
            let damping = -5.0 * vel;
            vel += (force + damping) * dt;
            pos += vel * dt;
        }
        assert!(
            (pos - 1.0).abs() < 0.1,
            "Position motor should converge to 1.0 m, got {pos}"
        );
    }
    #[test]
    fn linear_motor_combined_mode() {
        let mut motor = LinearMotorConstraint::combined_mode(1.0, 0.0, 1000.0);
        let force = motor.step(0.01, 0.0, 0.0);
        assert!(force > 0.0, "Should produce force toward target");
    }
    #[test]
    fn linear_motor_with_friction() {
        let mut motor = LinearMotorConstraint::velocity_mode(5.0, 1000.0);
        let friction = MotorFriction::coulomb_viscous(0.5, 0.1);
        let force = motor.step_with_friction(0.01, 1.0, 0.0, &friction);
        assert!(force.abs() <= 1000.0);
    }
    #[test]
    fn angular_motor_reaches_target_velocity() {
        let mut motor = AngularMotorConstraint::new(10.0, 1000.0);
        let dt = 0.005;
        let mut ang_vel = 0.0f64;
        let inertia = 1.0f64;
        for _ in 0..1000 {
            let torque = motor.step(dt, ang_vel, 0.0);
            ang_vel += torque / inertia * dt;
        }
        assert!(
            (ang_vel - 10.0).abs() < 1.0,
            "Angular motor should reach 10 rad/s, got {ang_vel}"
        );
    }
    #[test]
    fn angular_motor_torque_clamped() {
        let mut motor = AngularMotorConstraint::new(1000.0, 5.0);
        let torque = motor.step(0.01, 0.0, 0.0);
        assert!(
            torque.abs() <= 5.0,
            "Torque must not exceed max_torque: {torque}"
        );
    }
    #[test]
    fn angular_motor_negative_target() {
        let mut motor = AngularMotorConstraint::new(-5.0, 1000.0);
        let dt = 0.005;
        let mut ang_vel = 0.0f64;
        for _ in 0..1000 {
            let torque = motor.step(dt, ang_vel, 0.0);
            ang_vel += torque * dt;
        }
        assert!(
            ang_vel < -4.0,
            "Motor should spin backward to -5 rad/s, got {ang_vel}"
        );
    }
    #[test]
    fn angular_motor_set_target() {
        let mut motor = AngularMotorConstraint::new(5.0, 100.0);
        motor.set_target(10.0);
        assert_eq!(motor.target_angular_velocity, 10.0);
    }
    #[test]
    fn servo_converges_to_target_position() {
        let mut servo = ServoMotorConstraint::with_pd_gains(1.0, 1000.0, 200.0, 20.0);
        let dt = 0.01;
        let mut pos = 0.0f64;
        let mut vel = 0.0f64;
        let mass = 1.0f64;
        for _ in 0..1000 {
            let torque = servo.step(dt, vel, pos);
            vel += torque / mass * dt;
            pos += vel * dt;
        }
        assert!(
            (pos - 1.0).abs() < 0.05,
            "Servo should converge to 1.0, got {pos}"
        );
    }
    #[test]
    fn servo_output_clamped() {
        let mut servo = ServoMotorConstraint::new(100.0, 5.0, 1000.0, 0.0, 0.0);
        let out = servo.step(0.01, 0.0, 0.0);
        assert!(
            out.abs() <= 5.0,
            "Servo output must be within max_torque: {out}"
        );
    }
    #[test]
    fn servo_set_target_updates() {
        let mut servo = ServoMotorConstraint::with_pd_gains(0.0, 100.0, 10.0, 1.0);
        servo.set_target(2.0);
        assert_eq!(servo.target_position, 2.0);
    }
    #[test]
    fn servo_reset_clears_pid() {
        let mut servo = ServoMotorConstraint::new(1.0, 100.0, 1.0, 1.0, 0.0);
        servo.step(0.1, 0.0, 0.0);
        servo.reset();
        assert_eq!(servo.pid.integral, 0.0);
        assert_eq!(servo.pid.prev_error, 0.0);
    }
    #[test]
    fn impedance_control_force_at_target() {
        let motor = ImpedanceControlMotor::new(100.0, 10.0, 1000.0);
        let force = motor.step(0.0, 0.0);
        assert_eq!(force, 0.0, "At target → zero force");
    }
    #[test]
    fn impedance_control_stiffness() {
        let mut motor = ImpedanceControlMotor::new(100.0, 0.0, 1000.0);
        motor.set_desired(1.0, 0.0);
        let force = motor.step(0.0, 0.0);
        assert!((force - 100.0).abs() < 1e-10, "F = K * x_error = 100");
    }
    #[test]
    fn impedance_control_damping() {
        let mut motor = ImpedanceControlMotor::new(0.0, 10.0, 1000.0);
        motor.set_desired(0.0, 5.0);
        let force = motor.step(0.0, 0.0);
        assert!((force - 50.0).abs() < 1e-10, "F = D * v_error = 50");
    }
    #[test]
    fn impedance_control_clamped() {
        let mut motor = ImpedanceControlMotor::new(1e6, 0.0, 100.0);
        motor.set_desired(1.0, 0.0);
        let force = motor.step(0.0, 0.0);
        assert!(
            (force - 100.0).abs() < 1e-10,
            "Force should be clamped to 100"
        );
    }
    #[test]
    fn impedance_magnitude() {
        let motor = ImpedanceControlMotor::new(100.0, 10.0, 1000.0);
        let z = motor.impedance_magnitude(10.0);
        let expected = (100.0f64 + 100.0).sqrt();
        assert!((z - expected).abs() < 1e-10);
    }
    #[test]
    fn motor_solver_batch_step() {
        let mut solver = MotorSolver::new();
        let linear = MotorConstraint::Linear(LinearMotorConstraint::velocity_mode(1.0, 100.0));
        let angular = MotorConstraint::Angular(AngularMotorConstraint::new(2.0, 100.0));
        solver.add(linear);
        solver.add(angular);
        let outputs = solver.step(0.01, &[0.0, 0.0], &[0.0, 0.0]);
        assert_eq!(outputs.len(), 2);
        assert!(outputs[0] != 0.0, "Linear motor should produce force");
        assert!(outputs[1] != 0.0, "Angular motor should produce torque");
    }
    #[test]
    fn motor_solver_empty() {
        let mut solver = MotorSolver::new();
        let outputs = solver.step(0.01, &[], &[]);
        assert!(outputs.is_empty());
    }
    #[test]
    fn motor_solver_add_returns_index() {
        let mut solver = MotorSolver::new();
        let m1 = MotorConstraint::Angular(AngularMotorConstraint::new(0.0, 10.0));
        let m2 = MotorConstraint::Angular(AngularMotorConstraint::new(0.0, 10.0));
        let i1 = solver.add(m1);
        let i2 = solver.add(m2);
        assert_eq!(i1, 0);
        assert_eq!(i2, 1);
    }
    #[test]
    fn motor_solver_remove() {
        let mut solver = MotorSolver::new();
        solver.add(MotorConstraint::Angular(AngularMotorConstraint::new(
            0.0, 10.0,
        )));
        solver.add(MotorConstraint::Angular(AngularMotorConstraint::new(
            0.0, 10.0,
        )));
        assert_eq!(solver.len(), 2);
        solver.remove(0);
        assert_eq!(solver.len(), 1);
    }
    #[test]
    fn motor_solver_is_empty() {
        let solver = MotorSolver::new();
        assert!(solver.is_empty());
    }
    #[test]
    fn motor_solver_energy_expenditure() {
        let solver = MotorSolver::new();
        let forces = vec![10.0, -5.0];
        let velocities = vec![2.0, 3.0];
        let energy = solver.energy_expenditure(&forces, &velocities, 0.01);
        assert!((energy - 0.35).abs() < 1e-10);
    }
    #[test]
    fn motor_solver_impedance_variant() {
        let mut solver = MotorSolver::new();
        let imp = ImpedanceControlMotor::new(100.0, 10.0, 1000.0);
        solver.add(MotorConstraint::Impedance(imp));
        let outputs = solver.step(0.01, &[0.0], &[0.0]);
        assert_eq!(outputs.len(), 1);
    }
    #[test]
    fn motor_thermal_copper_loss() {
        let thermal = MotorThermal::new(100.0, 0.5, 20.0, 2.0, 150.0);
        let p = thermal.copper_loss(3.0);
        assert!((p - 18.0).abs() < 1e-10, "Copper loss = {p}");
    }
    #[test]
    fn motor_thermal_iron_loss() {
        let thermal = MotorThermal::new(100.0, 0.5, 20.0, 2.0, 150.0);
        let p = thermal.iron_loss(100.0, 0.001, 0.01);
        assert!((p - 11.0).abs() < 1e-8, "Iron loss = {p}");
    }
    #[test]
    fn motor_thermal_step_heats_up() {
        let mut thermal = MotorThermal::new(100.0, 0.5, 20.0, 2.0, 150.0);
        let t_before = thermal.temperature;
        thermal.step(100.0, 1.0);
        assert!(thermal.temperature > t_before, "Motor should heat up");
    }
    #[test]
    fn motor_thermal_steady_state() {
        let thermal = MotorThermal::new(100.0, 0.5, 20.0, 2.0, 150.0);
        let t_ss = thermal.steady_state_temperature(50.0);
        assert!((t_ss - 45.0).abs() < 1e-10, "T_ss = {t_ss}");
    }
    #[test]
    fn motor_thermal_time_constant() {
        let thermal = MotorThermal::new(200.0, 0.5, 20.0, 1.0, 150.0);
        let tau = thermal.time_constant();
        assert!((tau - 100.0).abs() < 1e-10, "tau = {tau}");
    }
    #[test]
    fn motor_thermal_overheated_flag() {
        let mut thermal = MotorThermal::new(10.0, 0.01, 20.0, 1.0, 50.0);
        thermal.temperature = 100.0;
        assert!(thermal.is_overheated(), "Should detect overheat");
        thermal.reset();
        assert!(!thermal.is_overheated(), "After reset: not overheated");
    }
    #[test]
    fn motor_thermal_reset_to_ambient() {
        let mut thermal = MotorThermal::new(100.0, 0.5, 25.0, 1.0, 150.0);
        thermal.temperature = 80.0;
        thermal.reset();
        assert!((thermal.temperature - 25.0).abs() < 1e-10);
    }
    #[test]
    fn drive_cycle_empty_returns_zero() {
        let dc = DriveCycle::new();
        assert_eq!(dc.velocity_at(1.0), 0.0);
    }
    #[test]
    fn drive_cycle_single_point() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(0.0, 5.0);
        assert_eq!(dc.velocity_at(0.0), 5.0);
        assert_eq!(dc.velocity_at(10.0), 5.0);
    }
    #[test]
    fn drive_cycle_linear_interpolation() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(0.0, 0.0);
        dc.add_setpoint(10.0, 10.0);
        let v = dc.velocity_at(5.0);
        assert!((v - 5.0).abs() < 1e-10, "Interpolated v = {v}");
    }
    #[test]
    fn drive_cycle_duration() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(0.0, 0.0);
        dc.add_setpoint(5.0, 1.0);
        dc.add_setpoint(10.0, 0.0);
        assert!((dc.duration() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn drive_cycle_peak_velocity() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(0.0, 0.0);
        dc.add_setpoint(2.0, 5.0);
        dc.add_setpoint(4.0, -8.0);
        dc.add_setpoint(6.0, 0.0);
        assert!(
            (dc.peak_velocity() - 8.0).abs() < 1e-10,
            "Peak = {}",
            dc.peak_velocity()
        );
    }
    #[test]
    fn drive_cycle_rms_velocity_constant() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(0.0, 3.0);
        dc.add_setpoint(10.0, 3.0);
        let rms = dc.rms_velocity(0.1);
        assert!((rms - 3.0).abs() < 0.01, "RMS of constant = {rms}");
    }
    #[test]
    fn drive_cycle_setpoints_sorted() {
        let mut dc = DriveCycle::new();
        dc.add_setpoint(5.0, 2.0);
        dc.add_setpoint(1.0, 0.0);
        dc.add_setpoint(3.0, 1.0);
        assert!(dc.setpoints[0].0 < dc.setpoints[1].0, "Should be sorted");
    }
    #[test]
    fn regen_braking_zero_below_min_speed() {
        let mut rb = RegenerativeBraking::new(0.9, 100.0, 1.0);
        let t = rb.compute_regen_torque(50.0, 0.5, 0.01);
        assert_eq!(t, 0.0, "Below min speed: no regen");
    }
    #[test]
    fn regen_braking_torque_negative() {
        let mut rb = RegenerativeBraking::new(0.9, 100.0, 0.1);
        let t = rb.compute_regen_torque(50.0, 10.0, 0.01);
        assert!(t < 0.0, "Regen torque should oppose motion: {t}");
    }
    #[test]
    fn regen_braking_energy_accumulates() {
        let mut rb = RegenerativeBraking::new(0.9, 100.0, 0.1);
        rb.compute_regen_torque(50.0, 10.0, 1.0);
        assert!(rb.energy_recovered > 0.0, "Energy should accumulate");
    }
    #[test]
    fn regen_braking_limited_by_max_torque() {
        let mut rb = RegenerativeBraking::new(1.0, 10.0, 0.1);
        let t = rb.compute_regen_torque(1000.0, 5.0, 0.01);
        assert!(t.abs() <= 10.0, "Regen torque capped at max: {t}");
    }
    #[test]
    fn regen_braking_reset_energy() {
        let mut rb = RegenerativeBraking::new(0.9, 100.0, 0.1);
        rb.compute_regen_torque(50.0, 10.0, 1.0);
        rb.reset_energy();
        assert_eq!(rb.energy_recovered, 0.0);
    }
    #[test]
    fn regen_braking_power() {
        let rb = RegenerativeBraking::new(0.8, 100.0, 0.0);
        let p = rb.regen_power(50.0, 10.0);
        assert!((p - 400.0).abs() < 1e-10, "Regen power = {p}");
    }
    #[test]
    fn fault_detector_no_fault() {
        let mut fd = MotorFaultDetector::new(150.0, 10.0, 0.1, 100.0);
        let fault = fd.check(25.0, 1.0, 0.001, 50.0);
        assert_eq!(*fault, MotorFault::None);
    }
    #[test]
    fn fault_detector_over_temperature() {
        let mut fd = MotorFaultDetector::new(100.0, 10.0, 0.1, 100.0);
        let fault = fd.check(150.0, 1.0, 0.001, 50.0);
        assert_eq!(*fault, MotorFault::OverTemperature);
    }
    #[test]
    fn fault_detector_over_current() {
        let mut fd = MotorFaultDetector::new(150.0, 5.0, 0.1, 100.0);
        let fault = fd.check(25.0, 10.0, 0.001, 50.0);
        assert_eq!(*fault, MotorFault::OverCurrent);
    }
    #[test]
    fn fault_detector_over_speed() {
        let mut fd = MotorFaultDetector::new(150.0, 10.0, 0.1, 50.0);
        let fault = fd.check(25.0, 1.0, 0.001, 100.0);
        assert_eq!(*fault, MotorFault::OverSpeed);
    }
    #[test]
    fn fault_detector_position_error() {
        let mut fd = MotorFaultDetector::new(150.0, 10.0, 0.05, 100.0);
        let fault = fd.check(25.0, 1.0, 0.1, 50.0);
        assert_eq!(*fault, MotorFault::PositionError);
    }
    #[test]
    fn fault_detector_stall() {
        let mut fd = MotorFaultDetector::new(150.0, 10.0, 0.1, 100.0);
        fd.stall_speed_threshold = 0.5;
        fd.stall_current_threshold = 0.5;
        let fault = fd.check(25.0, 9.0, 0.001, 0.1);
        assert_eq!(*fault, MotorFault::Stall);
    }
    #[test]
    fn fault_detector_has_fault_and_clear() {
        let mut fd = MotorFaultDetector::new(100.0, 5.0, 0.1, 50.0);
        fd.check(200.0, 1.0, 0.0, 10.0);
        assert!(fd.has_fault());
        fd.clear();
        assert!(!fd.has_fault());
    }
    #[test]
    fn multi_axis_step_produces_forces() {
        let mut coord = MultiAxisCoordinator::new(3, 100.0, 1000.0);
        coord.set_targets(&[1.0, 2.0, 3.0]);
        let forces = coord.step(0.01, &[0.0, 0.0, 0.0]);
        assert_eq!(forces.len(), 3);
        for &f in &forces {
            assert!(f > 0.0, "Should produce positive force toward target: {f}");
        }
    }
    #[test]
    fn multi_axis_force_clamped() {
        let mut coord = MultiAxisCoordinator::new(2, 1e6, 10.0);
        coord.set_targets(&[100.0, 100.0]);
        let forces = coord.step(0.01, &[0.0, 0.0]);
        for &f in &forces {
            assert!(f.abs() <= 10.0, "Force exceeds limit: {f}");
        }
    }
    #[test]
    fn multi_axis_zero_error_zero_force() {
        let mut coord = MultiAxisCoordinator::new(2, 100.0, 1000.0);
        coord.set_targets(&[1.0, 2.0]);
        let forces = coord.step(0.01, &[1.0, 2.0]);
        for &f in &forces {
            assert!(f.abs() < 1e-10, "At target: force should be zero: {f}");
        }
    }
    #[test]
    fn multi_axis_sync_error() {
        let mut coord = MultiAxisCoordinator::new(2, 100.0, 1000.0);
        coord.set_targets(&[1.0, 1.0]);
        let err = coord.sync_error(&[0.5, 0.8]);
        assert!((err - 0.3).abs() < 1e-10, "Sync error = {err}");
    }
    #[test]
    fn multi_axis_reset_all_pids() {
        let mut coord = MultiAxisCoordinator::new(2, 100.0, 1000.0);
        coord.set_targets(&[1.0, 1.0]);
        coord.step(0.01, &[0.0, 0.0]);
        coord.reset_all();
        for pid in &coord.pids {
            assert_eq!(pid.integral, 0.0);
        }
    }
    #[test]
    fn multi_axis_coupling_effect() {
        let mut coord = MultiAxisCoordinator::new(2, 100.0, 1000.0);
        coord.set_targets(&[1.0, 0.0]);
        coord.set_coupling(1, 0, 10.0);
        let forces = coord.step(0.01, &[0.0, 0.0]);
        assert!(
            forces[1].abs() > 0.0,
            "Coupling should create force on axis 1"
        );
    }
}
#[cfg(test)]
mod extended_motor_tests {

    use crate::motor_constraints::ActiveCompliance;
    use crate::motor_constraints::BallScrew;
    use crate::motor_constraints::CablePulleySystem;
    use crate::motor_constraints::HarmonicDrive;
    use crate::motor_constraints::ImpedanceController;
    use crate::motor_constraints::RackAndPinionMotor;
    use crate::motor_constraints::SeriesElasticActuator;
    use crate::motor_constraints::WormGear;
    #[test]
    fn harmonic_drive_output_torque_amplified() {
        let hd = HarmonicDrive::new(100.0, 10_000.0, 0.01);
        let t_out = hd.output_torque(1.0);
        assert!(
            (t_out - 100.0).abs() < 1e-10,
            "Torque should be amplified by ratio: {t_out}"
        );
    }
    #[test]
    fn harmonic_drive_reflected_inertia() {
        let hd = HarmonicDrive::new(100.0, 10_000.0, 0.1);
        let reflected = hd.reflected_inertia();
        assert!(
            (reflected - 1e-5).abs() < 1e-15,
            "reflected inertia: {reflected}"
        );
    }
    #[test]
    fn harmonic_drive_input_angle() {
        let mut hd = HarmonicDrive::new(50.0, 0.0, 0.01);
        hd.output_angle = 0.1;
        assert!((hd.input_angle() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn harmonic_drive_step_moves_output() {
        let mut hd = HarmonicDrive::new(100.0, 0.0, 0.01);
        hd.step(1.0, 0.01);
        assert!(hd.output_velocity != 0.0, "output velocity should change");
    }
    #[test]
    fn harmonic_drive_friction_opposes_motion() {
        let mut hd = HarmonicDrive::new(100.0, 0.0, 0.01);
        hd = hd.with_friction(5.0);
        hd.output_velocity = 1.0;
        let t = hd.output_torque(0.0);
        assert!(
            t < 0.0,
            "Friction should produce negative torque when moving forward"
        );
    }
    #[test]
    fn worm_gear_ratio() {
        let wg = WormGear::new(2, 50, 15.0, 0.1);
        assert!((wg.gear_ratio() - 25.0).abs() < 1e-10);
    }
    #[test]
    fn worm_gear_driving_efficiency_positive() {
        let wg = WormGear::new(2, 50, 15.0, 0.1);
        let eff = wg.driving_efficiency();
        assert!(eff > 0.0 && eff <= 1.0, "driving efficiency: {eff}");
    }
    #[test]
    fn worm_gear_self_locking_at_low_lead_angle() {
        let wg = WormGear::new(1, 50, 2.0, 0.15);
        assert!(
            wg.is_self_locking(),
            "Should be self-locking at low lead angle"
        );
    }
    #[test]
    fn worm_gear_not_self_locking_at_high_lead_angle() {
        let wg = WormGear::new(4, 50, 45.0, 0.05);
        assert!(!wg.is_self_locking(), "High lead angle → not self-locking");
    }
    #[test]
    fn worm_gear_output_torque_accounts_efficiency() {
        let wg = WormGear::new(2, 50, 15.0, 0.1);
        let t_ideal = 1.0 * wg.gear_ratio();
        let t_actual = wg.output_torque(1.0);
        assert!(
            t_actual < t_ideal,
            "Actual torque should be less than ideal due to friction"
        );
    }
    #[test]
    fn ball_screw_linear_velocity() {
        let bs = BallScrew::new(0.005, 0.95, 5000.0);
        let v = bs.linear_velocity(2.0 * std::f64::consts::PI);
        assert!((v - 0.005).abs() < 1e-6, "1 rev/s → 5 mm/s: {v}");
    }
    #[test]
    fn ball_screw_axial_force_clamped() {
        let bs = BallScrew::new(0.005, 0.95, 1000.0);
        let f = bs.axial_force(1e6);
        assert!((f - 1000.0).abs() < 1e-6, "Should be clamped to max: {f}");
    }
    #[test]
    fn ball_screw_round_trip_torque_force() {
        let bs = BallScrew::new(0.01, 0.90, 10_000.0);
        let force = 500.0;
        let torque = bs.required_torque(force);
        let f_back = bs.axial_force(torque);
        assert!(
            f_back <= force + 0.01,
            "Round trip: f_back={f_back}, original={force}"
        );
    }
    #[test]
    fn ball_screw_step_advances_position() {
        let mut bs = BallScrew::new(0.01, 0.95, 10_000.0);
        let omega = 10.0;
        bs.step(omega, 1.0);
        let expected = bs.linear_velocity(omega);
        assert!(
            bs.position.abs() > 0.0,
            "Position should advance: {}",
            bs.position
        );
        let _ = expected;
    }
    #[test]
    fn rack_pinion_motor_rack_force() {
        let rp = RackAndPinionMotor::new(0.02, 0.9, 1000.0);
        let f = rp.rack_force(10.0);
        assert!((f - 450.0).abs() < 1e-6, "rack force: {f}");
    }
    #[test]
    fn rack_pinion_motor_required_torque() {
        let rp = RackAndPinionMotor::new(0.02, 0.9, 1000.0);
        let t = rp.required_torque(450.0);
        assert!((t - 10.0).abs() < 1e-6, "required torque: {t}");
    }
    #[test]
    fn rack_pinion_motor_step() {
        let mut rp = RackAndPinionMotor::new(0.05, 0.95, 5000.0);
        rp.step(10.0, 0.1);
        assert!(
            (rp.position - 0.05).abs() < 1e-10,
            "position: {}",
            rp.position
        );
    }
    #[test]
    fn cable_pulley_tension_zero_at_rest() {
        let cp = CablePulleySystem::new(0.05, 5000.0, 1.0);
        assert_eq!(cp.tension(), 0.0);
    }
    #[test]
    fn cable_pulley_tension_proportional_to_elongation() {
        let mut cp = CablePulleySystem::new(0.05, 5000.0, 1.0);
        cp.elongation = 0.01;
        let t = cp.tension();
        assert!((t - 50.0).abs() < 1e-6, "tension: {t}");
    }
    #[test]
    fn cable_pulley_elastic_energy() {
        let mut cp = CablePulleySystem::new(0.05, 2000.0, 1.0);
        cp.elongation = 0.02;
        let e = cp.elastic_energy();
        assert!((e - 0.4).abs() < 1e-10, "energy: {e}");
    }
    #[test]
    fn cable_pulley_update_elongation() {
        let mut cp = CablePulleySystem::new(0.05, 5000.0, 1.0);
        cp.update_elongation(0.05, 0.01);
        let expected = (0.01_f64 - 0.05 * 0.05).max(0.0);
        assert!((cp.elongation - expected).abs() < 1e-12);
    }
    #[test]
    fn sea_output_force_proportional_to_deflection() {
        let mut sea = SeriesElasticActuator::new(1000.0, 1.0);
        sea.motor_position = 0.01;
        sea.load_position = 0.0;
        assert!((sea.output_force() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn sea_force_error() {
        let mut sea = SeriesElasticActuator::new(1000.0, 1.0);
        sea.target_force = 10.0;
        sea.motor_position = 0.005;
        assert!((sea.force_error() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn sea_desired_motor_velocity() {
        let mut sea = SeriesElasticActuator::new(1000.0, 0.001);
        sea.target_force = 10.0;
        let v = sea.desired_motor_velocity();
        assert!((v - 0.01).abs() < 1e-12);
    }
    #[test]
    fn sea_elastic_energy() {
        let mut sea = SeriesElasticActuator::new(500.0, 1.0);
        sea.motor_position = 0.02;
        assert!((sea.elastic_energy() - 0.1).abs() < 1e-10);
    }
    #[test]
    fn active_compliance_zero_error_zero_force() {
        let ac = ActiveCompliance::new(1000.0, 50.0, 500.0);
        let f = ac.compliance_force(0.0, 0.0);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn active_compliance_restoring_force() {
        let ac = ActiveCompliance::new(1000.0, 0.0, 5000.0);
        let f = ac.compliance_force(0.01, 0.0);
        assert!((f - (-10.0)).abs() < 1e-10);
    }
    #[test]
    fn active_compliance_damping_ratio() {
        let ac = ActiveCompliance::new(1000.0, 20.0, 1000.0);
        let ratio = ac.damping_ratio(1.0);
        assert!((ratio - 0.3162).abs() < 0.001, "damping ratio: {ratio}");
    }
    #[test]
    fn active_compliance_max_force_clamped() {
        let ac = ActiveCompliance::new(1e6, 0.0, 100.0);
        let f = ac.compliance_force(1.0, 0.0);
        assert!((f.abs() - 100.0).abs() < 1e-10, "Should be clamped: {f}");
    }
    #[test]
    fn active_compliance_impedance_magnitude_increases_with_freq() {
        let ac = ActiveCompliance::new(100.0, 10.0, 1000.0);
        let z_low = ac.impedance_magnitude(1.0);
        let z_high = ac.impedance_magnitude(100.0);
        assert!(z_high > z_low, "higher ω → higher impedance magnitude");
    }
    #[test]
    fn impedance_controller_natural_frequency() {
        let ic = ImpedanceController::new(1.0, 100.0, 10.0, 1000.0);
        let wn = ic.natural_frequency();
        assert!((wn - 10.0).abs() < 1e-10);
    }
    #[test]
    fn impedance_controller_damping_ratio() {
        let ic = ImpedanceController::new(1.0, 100.0, 20.0, 1000.0);
        let ratio = ic.damping_ratio();
        assert!((ratio - 1.0).abs() < 1e-10, "ratio: {ratio}");
    }
    #[test]
    fn impedance_controller_compute_force_clamped() {
        let mut ic = ImpedanceController::new(1.0, 1e8, 1e6, 50.0);
        ic.set_desired(1.0, 0.0, 0.0);
        let f = ic.compute_force(0.0, 0.0, 0.0);
        assert!((f.abs() - 50.0).abs() < 1e-10, "Should be clamped: {f}");
    }
    #[test]
    fn impedance_controller_external_force_compensation() {
        let mut ic = ImpedanceController::new(1.0, 0.0, 0.0, 1000.0);
        ic.set_desired(0.0, 0.0, 0.0);
        let f = ic.compute_force(0.0, 0.0, 10.0);
        assert!((f - (-10.0)).abs() < 1e-10, "Force: {f}");
    }
    #[test]
    fn impedance_controller_critical_damping() {
        let ic = ImpedanceController::new(2.0, 200.0, 40.0, 1000.0);
        let dc = ic.critical_damping();
        assert!((dc - 40.0).abs() < 1e-8, "critical damping: {dc}");
    }
}
