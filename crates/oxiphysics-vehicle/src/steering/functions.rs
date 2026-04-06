//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::AckermannSteering;
    use crate::steering::AckermannGeometry;

    use crate::steering::FourWheelSteering;

    use crate::steering::FullFourWheelSteering;

    use crate::steering::RearWheelSteering;
    use crate::steering::SpeedDependentRatio;
    use crate::steering::SteerByWire;

    use crate::steering::SteerLimiter;
    use crate::steering::SteeringColumn;

    use crate::steering::SteeringFeedback;
    use crate::steering::SteeringFeel;

    use crate::steering::TorqueSteerCompensation;

    #[test]
    fn test_ackermann_straight() {
        let steering = AckermannSteering::default();
        let (left, right) = steering.compute_wheel_angles(0.0);
        assert!((left).abs() < 1e-10);
        assert!((right).abs() < 1e-10);
    }
    #[test]
    fn test_ackermann_inner_greater_than_outer() {
        let steering = AckermannSteering::new(2.5, 1.5, 0.6);
        let (left, right) = steering.compute_wheel_angles(0.5);
        assert!(left > 0.0);
        assert!(right > 0.0);
        assert!(
            right > left,
            "inner(right)={right} should > outer(left)={left}"
        );
    }
    #[test]
    fn test_ackermann_symmetry() {
        let steering = AckermannSteering::new(2.5, 1.5, 0.6);
        let (left_r, right_r) = steering.compute_wheel_angles(0.5);
        let (left_l, right_l) = steering.compute_wheel_angles(-0.5);
        assert!((left_r - (-right_l)).abs() < 1e-10);
        assert!((right_r - (-left_l)).abs() < 1e-10);
    }
    #[test]
    fn test_turning_radius() {
        let steering = AckermannSteering::new(2.5, 1.5, 0.6);
        let r = steering.turning_radius(1.0);
        assert!(r.is_some());
        let radius = r.unwrap();
        assert!(radius > 0.0 && radius < 100.0, "radius={radius}");
    }
    #[test]
    fn test_turning_radius_straight() {
        let steering = AckermannSteering::default();
        assert!(steering.turning_radius(0.0).is_none());
    }
    #[test]
    fn test_parallel_steering_factor_zero() {
        let steering = AckermannSteering::new(2.5, 1.5, 0.6).with_ackermann_factor(0.0);
        let (left, right) = steering.compute_wheel_angles(0.5);
        assert!((left - right).abs() < 1e-10, "left={left}, right={right}");
    }
    #[test]
    fn test_ackermann_steering() {
        let wheelbase = 2.5_f64;
        let track_width = 1.5_f64;
        let steering = AckermannSteering::new(wheelbase, track_width, 0.6);
        let (left_angle, right_angle) = steering.compute_wheel_angles(0.5);
        assert!(left_angle > 0.0 && right_angle > 0.0);
        assert!(
            right_angle > left_angle,
            "inner(right)={right_angle:.4} should > outer(left)={left_angle:.4}"
        );
        let outer_angle = left_angle;
        let inner_angle = right_angle;
        let cot_outer = 1.0 / outer_angle.tan();
        let cot_inner = 1.0 / inner_angle.tan();
        let ackermann_diff = cot_outer - cot_inner;
        let expected_diff = track_width / wheelbase;
        assert!(
            (ackermann_diff - expected_diff).abs() < 1e-6,
            "cot(outer)-cot(inner)={ackermann_diff:.6}, expected W/L={expected_diff:.6}"
        );
    }
    #[test]
    fn test_geometry_inner_outer() {
        let g = AckermannGeometry::new(2.5, 1.5);
        let (inner, outer) = g.steering_angles(0.3);
        assert!(inner > 0.0 && outer > 0.0);
        assert!(inner > outer, "inner={inner:.4} should > outer={outer:.4}");
    }
    #[test]
    fn test_geometry_turning_radius_zero_is_infinity() {
        let g = AckermannGeometry::new(2.5, 1.5);
        let r = g.turning_radius(0.0);
        assert!(r.is_infinite(), "expected infinity, got {r}");
    }
    #[test]
    fn test_geometry_max_steer_angle_positive() {
        let g = AckermannGeometry::new(2.5, 1.5);
        let max_a = g.max_steering_angle();
        assert!(max_a > 0.0 && max_a < std::f64::consts::FRAC_PI_2);
    }
    #[test]
    fn test_servo_rate_limiting() {
        let sbw = SteerByWire::new(1.0, 0.5, 10.0);
        let next = sbw.update(0.0, 1.0);
        assert!((next - 0.5).abs() < 1e-10, "expected 0.5 rad, got {next}");
    }
    #[test]
    fn test_servo_reaches_target() {
        let sbw = SteerByWire::new(0.1, 1.0, 10.0);
        let next = sbw.update(0.0, 1.0);
        assert!((next - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_torque_feedback_clamped() {
        let sbw = SteerByWire::new(0.0, 1.0, 5.0);
        let t = sbw.apply_torque_feedback(4.0, 4.0);
        assert!((t - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_speed_sensitive_rear_steering() {
        let fws = FourWheelSteering::new(0.3, 30.0);
        let front = 0.2_f64;
        let rear_low = fws.rear_angle(front, 0.0);
        let rear_high = fws.rear_angle(front, 60.0);
        assert!(
            rear_low < 0.0,
            "expected counter-steer at low speed, got {rear_low}"
        );
        assert!(
            rear_high > 0.0,
            "expected in-phase at high speed, got {rear_high}"
        );
    }
    #[test]
    fn test_dead_zone_suppresses_small_inputs() {
        let mut col = SteeringColumn::new(0.1, 0.0, 0.0);
        let out = col.process_input(0.05, 0.0);
        assert!(
            out.abs() < 1e-10,
            "small input inside dead zone should be zero, got {out}"
        );
    }
    #[test]
    fn test_dead_zone_passes_large_inputs() {
        let mut col = SteeringColumn::new(0.1, 0.0, 0.0);
        let out = col.process_input(0.5, 0.0);
        assert!(
            out.abs() > 0.0,
            "large input outside dead zone should be non-zero"
        );
    }
    #[test]
    fn test_steer_limiter_at_zero_speed() {
        let lim = SteerLimiter::new(0.6, 0.2, 50.0);
        let out = lim.limit(0.6, 0.0);
        assert!((out - 0.6).abs() < 1e-10);
    }
    #[test]
    fn test_steer_limiter_at_max_speed() {
        let lim = SteerLimiter::new(0.6, 0.2, 50.0);
        let out = lim.limit(0.6, 50.0);
        assert!((out - 0.2).abs() < 1e-10, "expected 0.2, got {out}");
    }
    #[test]
    fn test_aligning_torque_sign() {
        let feel = SteeringFeel::new(0.05, 0.03, 0.02);
        let t = feel.aligning_torque(1000.0, 0.05);
        assert!(
            t < 0.0,
            "aligning torque should be negative (restoring), got {t}"
        );
    }
    #[test]
    fn test_rear_wheel_steering_angle() {
        let rws = RearWheelSteering::new(0.4, 1.0);
        let angle = rws.compute_rear_angle(0.5);
        assert!((angle - 0.2).abs() < 1e-10, "expected 0.2, got {angle}");
    }
    #[test]
    fn test_rear_wheel_steering_clamped() {
        let rws = RearWheelSteering::new(0.4, 2.0);
        let angle = rws.compute_rear_angle(1.0);
        assert!((angle - 0.4).abs() < 1e-10, "Should be clamped to max");
    }
    #[test]
    fn test_rear_wheel_steering_zero_input() {
        let rws = RearWheelSteering::new(0.4, 1.0);
        let angle = rws.compute_rear_angle(0.0);
        assert!(angle.abs() < 1e-10);
    }
    #[test]
    fn test_rear_wheel_turning_radius() {
        let rws = RearWheelSteering::new(0.4, 1.0);
        let r = rws.turning_radius(1.0, 2.5);
        assert!(r.is_some());
        assert!(r.unwrap() > 0.0);
    }
    #[test]
    fn test_rear_wheel_turning_radius_straight() {
        let rws = RearWheelSteering::new(0.4, 1.0);
        assert!(rws.turning_radius(0.0, 2.5).is_none());
    }
    #[test]
    fn test_speed_ratio_at_zero_speed() {
        let sdr = SpeedDependentRatio::new(12.0, 16.0, 30.0);
        assert!((sdr.ratio_at_speed(0.0) - 12.0).abs() < 1e-10);
    }
    #[test]
    fn test_speed_ratio_at_transition_speed() {
        let sdr = SpeedDependentRatio::new(12.0, 16.0, 30.0);
        assert!((sdr.ratio_at_speed(30.0) - 16.0).abs() < 1e-10);
    }
    #[test]
    fn test_speed_ratio_at_half_speed() {
        let sdr = SpeedDependentRatio::new(12.0, 16.0, 30.0);
        let r = sdr.ratio_at_speed(15.0);
        assert!((r - 14.0).abs() < 1e-10, "expected 14.0, got {r}");
    }
    #[test]
    fn test_speed_ratio_apply() {
        let sdr = SpeedDependentRatio::new(12.0, 16.0, 30.0);
        let road = sdr.apply(120.0_f64.to_radians(), 0.0);
        let expected = 120.0_f64.to_radians() / 12.0;
        assert!((road - expected).abs() < 1e-6);
    }
    #[test]
    fn test_steering_feedback_centering() {
        let fb = SteeringFeedback::new(0.0, 0.0, 100.0, 500.0);
        let force = fb.compute_force(0.0, 0.1, 0.0);
        assert!((force - (-10.0)).abs() < 1e-6, "expected -10, got {force}");
    }
    #[test]
    fn test_steering_feedback_clamped() {
        let fb = SteeringFeedback::new(1.0, 0.0, 0.0, 50.0);
        let force = fb.compute_force(1000.0, 0.0, 0.0);
        assert!((force.abs() - 50.0).abs() < 1e-6, "Should be clamped to 50");
    }
    #[test]
    fn test_steering_feedback_damping() {
        let fb = SteeringFeedback::new(0.0, 10.0, 0.0, 500.0);
        let force = fb.compute_force(0.0, 0.0, 1.0);
        assert!((force - (-10.0)).abs() < 1e-6, "Damping force");
    }
    #[test]
    fn test_torque_steer_correction_zero_torque() {
        let ts = TorqueSteerCompensation::new(0.001, 0.05, 1.0);
        assert!(ts.correction(0.0).abs() < 1e-10);
    }
    #[test]
    fn test_torque_steer_correction_positive_torque() {
        let ts = TorqueSteerCompensation::new(0.001, 0.05, 1.0);
        let c = ts.correction(100.0);
        assert!((c - (-0.05)).abs() < 1e-10, "expected -0.05, got {c}");
    }
    #[test]
    fn test_torque_steer_apply() {
        let ts = TorqueSteerCompensation::new(0.0001, 0.05, 1.0);
        let base = 0.1;
        let result = ts.apply(base, 50.0);
        let expected = base + ts.correction(50.0);
        assert!((result - expected).abs() < 1e-10);
    }
    #[test]
    fn test_full_4ws_straight() {
        let front = AckermannSteering::default();
        let rear = FourWheelSteering::new(0.3, 30.0);
        let ratio = SpeedDependentRatio::new(14.0, 18.0, 40.0);
        let fws = FullFourWheelSteering::new(front, rear, ratio);
        let (fl, fr, rl, rr) = fws.compute_all_angles(0.0, 10.0);
        assert!(fl.abs() < 1e-10);
        assert!(fr.abs() < 1e-10);
        assert!(rl.abs() < 1e-10);
        assert!(rr.abs() < 1e-10);
    }
    #[test]
    fn test_full_4ws_rear_angles_equal() {
        let front = AckermannSteering::default();
        let rear = FourWheelSteering::new(0.3, 30.0);
        let ratio = SpeedDependentRatio::new(14.0, 18.0, 40.0);
        let fws = FullFourWheelSteering::new(front, rear, ratio);
        let (_fl, _fr, rl, rr) = fws.compute_all_angles(0.5, 10.0);
        assert!((rl - rr).abs() < 1e-10, "Rear angles should be equal");
    }
    #[test]
    fn test_full_4ws_nonzero_turn() {
        let front = AckermannSteering::default();
        let rear = FourWheelSteering::new(0.3, 30.0);
        let ratio = SpeedDependentRatio::new(14.0, 18.0, 40.0);
        let fws = FullFourWheelSteering::new(front, rear, ratio);
        let (fl, fr, _rl, _rr) = fws.compute_all_angles(0.5, 10.0);
        assert!(fl.abs() > 0.0 || fr.abs() > 0.0);
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::steering::DynamicSteeringRatio;

    use crate::steering::RackAndPinion;
    use crate::steering::SpeedDependentRatio;

    use crate::steering::SteeringColumnCompliance;
    use crate::steering::SteeringFeel;
    use crate::steering::SteeringFeelEnhanced;
    use crate::steering::UndersteerAnalyzer;
    #[test]
    fn test_rack_travel_linear() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let travel = rp.rack_travel(1.0);
        assert!((travel - 0.02).abs() < 1e-10, "expected 0.02, got {travel}");
    }
    #[test]
    fn test_rack_travel_clamped() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let travel = rp.rack_travel(100.0);
        assert!((travel - 0.08).abs() < 1e-10, "should be clamped to max");
    }
    #[test]
    fn test_steer_angle_from_travel_small_angle() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let travel = 0.01;
        let angle = rp.steer_angle_from_travel(travel);
        let expected = (travel / 0.15_f64).atan();
        assert!((angle - expected).abs() < 1e-10);
    }
    #[test]
    fn test_steer_angle_end_to_end() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let angle = rp.steer_angle(1.0);
        let direct = rp.steer_angle_from_travel(rp.rack_travel(1.0));
        assert!((angle - direct).abs() < 1e-10);
    }
    #[test]
    fn test_mechanical_advantage_positive() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let ma = rp.mechanical_advantage(0.0);
        assert!(ma > 0.0, "mechanical advantage should be positive");
    }
    #[test]
    fn test_required_pinion_torque_proportional() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        let t1 = rp.required_pinion_torque(100.0);
        let t2 = rp.required_pinion_torque(200.0);
        assert!(
            (t2 - 2.0 * t1).abs() < 1e-10,
            "torque should be proportional to force"
        );
    }
    #[test]
    fn test_rack_deflection_zero_force() {
        let rp = RackAndPinion::new(0.02, 0.15, 0.08, 1e-7);
        assert!(rp.rack_deflection(0.0).abs() < 1e-12);
    }
    #[test]
    fn test_understeer_gradient_sign_understeer() {
        let a = UndersteerAnalyzer::new(2.5, 80_000.0, 100_000.0, 1500.0, 0.6);
        let k = a.understeer_gradient();
        assert!(k > 0.0, "expected understeer (K>0), got K={k}");
    }
    #[test]
    fn test_understeer_gradient_sign_oversteer() {
        let a = UndersteerAnalyzer::new(2.5, 120_000.0, 60_000.0, 1500.0, 0.35);
        let k = a.understeer_gradient();
        assert!(k < 0.0, "expected oversteer (K<0), got K={k}");
    }
    #[test]
    fn test_critical_speed_only_for_oversteer() {
        let us = UndersteerAnalyzer::new(2.5, 80_000.0, 100_000.0, 1500.0, 0.6);
        assert!(
            us.critical_speed().is_none(),
            "understeer has no critical speed"
        );
        let os = UndersteerAnalyzer::new(2.5, 120_000.0, 60_000.0, 1500.0, 0.35);
        assert!(
            os.critical_speed().is_some(),
            "oversteer should have critical speed"
        );
    }
    #[test]
    fn test_characteristic_speed_only_for_understeer() {
        let us = UndersteerAnalyzer::new(2.5, 80_000.0, 100_000.0, 1500.0, 0.6);
        assert!(us.characteristic_speed().is_some());
        let os = UndersteerAnalyzer::new(2.5, 120_000.0, 60_000.0, 1500.0, 0.35);
        assert!(os.characteristic_speed().is_none());
    }
    #[test]
    fn test_classify_understeer() {
        let a = UndersteerAnalyzer::new(2.5, 80_000.0, 100_000.0, 1500.0, 0.6);
        assert_eq!(a.classify(), "understeer");
    }
    #[test]
    fn test_classify_oversteer() {
        let a = UndersteerAnalyzer::new(2.5, 120_000.0, 60_000.0, 1500.0, 0.35);
        assert_eq!(a.classify(), "oversteer");
    }
    #[test]
    fn test_predicted_steer_angle_increases_with_speed_for_understeer() {
        let a = UndersteerAnalyzer::new(2.5, 80_000.0, 100_000.0, 1500.0, 0.6);
        let angle_low = a.predicted_steer_angle(10.0, 50.0);
        let angle_high = a.predicted_steer_angle(30.0, 50.0);
        assert!(
            angle_high > angle_low,
            "understeer: more steer at higher speed: low={angle_low:.4}, high={angle_high:.4}"
        );
    }
    #[test]
    fn test_column_twist_converges_to_zero_no_torque() {
        let mut col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        for _ in 0..10 {
            col.step(0.0, 0.0, 0.0, 0.01);
        }
        assert!(col.twist_angle().abs() < 1e-6, "twist should decay to zero");
    }
    #[test]
    fn test_column_twist_nonzero_under_torque() {
        let mut col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        col.step(50.0, 0.0, 0.0, 0.01);
        assert!(
            col.twist_angle().abs() > 0.0,
            "non-zero torque should produce twist"
        );
    }
    #[test]
    fn test_column_reset() {
        let mut col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        col.step(100.0, 0.0, 0.0, 0.01);
        col.reset();
        assert!(col.twist_angle().abs() < 1e-12, "reset should zero twist");
    }
    #[test]
    fn test_column_natural_frequency_positive() {
        let col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        let f = col.natural_frequency();
        assert!(f > 0.0, "natural frequency must be positive");
    }
    #[test]
    fn test_dynamic_ratio_reduces_at_large_angle() {
        let base = SpeedDependentRatio::new(14.0, 18.0, 40.0);
        let dsr = DynamicSteeringRatio::new(base, 4.0, 2.0);
        let ratio_small = dsr.ratio_at(10.0, 0.1);
        let ratio_large = dsr.ratio_at(10.0, 2.0);
        assert!(
            ratio_large < ratio_small,
            "ratio should decrease at large angles: small={ratio_small:.2}, large={ratio_large:.2}"
        );
    }
    #[test]
    fn test_dynamic_ratio_minimum_is_one() {
        let base = SpeedDependentRatio::new(2.0, 2.0, 40.0);
        let dsr = DynamicSteeringRatio::new(base, 100.0, 0.01);
        let ratio = dsr.ratio_at(0.0, 100.0);
        assert!(ratio >= 1.0, "ratio must not go below 1.0, got {ratio}");
    }
    #[test]
    fn test_enhanced_feel_zero_speed_gives_zero() {
        let base = SteeringFeel::new(0.05, 0.03, 0.02);
        let feel = SteeringFeelEnhanced::new(base, 0.05, 0.2, 30.0);
        let torque = feel.aligning_torque(1000.0, 0.05, 0.1, 0.0);
        assert!(
            torque.abs() < 1e-9,
            "zero speed should produce zero torque, got {torque}"
        );
    }
    #[test]
    fn test_enhanced_feel_increases_with_speed() {
        let base = SteeringFeel::new(0.05, 0.03, 0.02);
        let feel = SteeringFeelEnhanced::new(base, 0.05, 0.2, 30.0);
        let t_low = feel.aligning_torque(1000.0, 0.05, 0.3, 5.0).abs();
        let t_high = feel.aligning_torque(1000.0, 0.05, 0.3, 25.0).abs();
        assert!(
            t_high > t_low,
            "feel should increase with speed: low={t_low:.4}, high={t_high:.4}"
        );
    }
    #[test]
    fn test_enhanced_feel_on_centre_is_lighter() {
        let base = SteeringFeel::new(0.05, 0.03, 0.02);
        let feel = SteeringFeelEnhanced::new(base, 0.3, 0.1, 30.0);
        let t_centre = feel.aligning_torque(1000.0, 0.05, 0.0, 20.0).abs();
        let t_off = feel.aligning_torque(1000.0, 0.05, 0.5, 20.0).abs();
        assert!(
            t_centre < t_off,
            "on-centre feel should be lighter: centre={t_centre:.4}, off={t_off:.4}"
        );
    }
}
#[cfg(test)]
mod tests_steering_extended {

    use crate::steering::FrontRearSteeringCorrelation;
    use crate::steering::ModelPredictiveSteering;
    use crate::steering::OversteerUndersteerDetector;

    use crate::steering::SteerByWire;
    use crate::steering::SteerByWireSimulator;
    use crate::steering::SteeringColumnCompliance;
    use crate::steering::SteeringFeel;

    #[test]
    fn test_mpc_bicycle_step_straight() {
        let mpc = ModelPredictiveSteering::new(1.0, 10, 2.5, 0.5, 1.0, 0.5, 0.1);
        let (x, y, psi) = mpc.bicycle_step(0.0, 0.0, 0.0, 10.0, 0.0, 0.1);
        assert!((x - 1.0).abs() < 1e-10, "x should advance: {x}");
        assert!(y.abs() < 1e-10, "y should not change for straight: {y}");
        assert!(
            psi.abs() < 1e-10,
            "psi should not change for straight: {psi}"
        );
    }
    #[test]
    fn test_mpc_bicycle_step_turning() {
        let mpc = ModelPredictiveSteering::new(1.0, 10, 2.5, 0.5, 1.0, 0.5, 0.1);
        let (_, _, psi) = mpc.bicycle_step(0.0, 0.0, 0.0, 10.0, 0.1, 0.1);
        assert!(psi.abs() > 0.0, "turning should change heading: {psi}");
    }
    #[test]
    fn test_mpc_compute_cost_zero_error() {
        let mpc = ModelPredictiveSteering::new(1.0, 5, 2.5, 0.5, 1.0, 0.5, 0.1);
        let cost = mpc.compute_cost(0.0, 0.0, 0.0, 5.0, 0.0, 0.0);
        assert!(cost >= 0.0);
    }
    #[test]
    fn test_mpc_optimal_steer_corrects_cte() {
        let mpc = ModelPredictiveSteering::new(5.0, 50, 2.5, 0.5, 200.0, 0.0, 0.0);
        let steer = mpc.optimal_steer(0.0, 0.5, 0.0, 5.0, 0.0, 21);
        assert!(
            steer <= 0.0,
            "Should steer left or straight to correct positive CTE: {steer}"
        );
    }
    #[test]
    fn test_mpc_optimal_steer_within_limits() {
        let mpc = ModelPredictiveSteering::new(1.0, 5, 2.5, 0.3, 1.0, 0.5, 0.1);
        let steer = mpc.optimal_steer(0.0, 1.0, 0.0, 5.0, 0.0, 9);
        assert!(
            steer.abs() <= 0.3,
            "Optimal steer must stay within limits: {steer}"
        );
    }
    #[test]
    fn test_understeer_detector_expected_yaw_rate() {
        let det = OversteerUndersteerDetector::new(2.5, 0.01);
        let yaw = det.expected_yaw_rate(10.0, 0.1);
        assert!(yaw > 0.0, "positive steer → positive yaw: {yaw}");
    }
    #[test]
    fn test_understeer_detector_neutral() {
        let det = OversteerUndersteerDetector::new(2.5, 0.05);
        let expected = det.expected_yaw_rate(10.0, 0.1);
        let class = det.classify(expected, 10.0, 0.1);
        assert_eq!(class, "neutral");
    }
    #[test]
    fn test_understeer_detector_understeer_classification() {
        let det = OversteerUndersteerDetector::new(2.5, 0.01);
        let expected = det.expected_yaw_rate(10.0, 0.1);
        let class = det.classify(expected * 0.1, 10.0, 0.1);
        assert_eq!(class, "understeer");
    }
    #[test]
    fn test_understeer_detector_oversteer_classification() {
        let det = OversteerUndersteerDetector::new(2.5, 0.01);
        let expected = det.expected_yaw_rate(10.0, 0.1);
        let class = det.classify(expected * 2.0, 10.0, 0.1);
        assert_eq!(class, "oversteer");
    }
    #[test]
    fn test_understeer_oversteer_index_neutral() {
        let det = OversteerUndersteerDetector::new(2.5, 0.01);
        let expected = det.expected_yaw_rate(10.0, 0.05);
        let idx = det.oversteer_index(expected, 10.0, 0.05);
        assert!(idx.abs() < 1e-10, "Measured = expected → index ≈ 0: {idx}");
    }
    #[test]
    fn test_sbw_simulator_step_returns_angle() {
        let servo = SteerByWire::new(0.1, 1.0, 20.0);
        let feel = SteeringFeel::new(0.05, 0.03, 0.02);
        let col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        let mut sim = SteerByWireSimulator::new(servo, feel, col);
        let (angle, _feedback) = sim.step(0.1, 500.0, 0.05, 0.01);
        assert!(
            angle.abs() > 0.0 || angle.abs() < 1.0,
            "angle should be valid: {angle}"
        );
    }
    #[test]
    fn test_sbw_simulator_column_torque() {
        let servo = SteerByWire::new(0.0, 1.0, 20.0);
        let feel = SteeringFeel::new(0.05, 0.03, 0.02);
        let col = SteeringColumnCompliance::new(5000.0, 50.0, 14.0);
        let sim = SteerByWireSimulator::new(servo, feel, col);
        let torque = sim.column_torque(1000.0, 0.05);
        assert!(
            torque.abs() <= 20.0,
            "column torque must not exceed torque_limit: {torque}"
        );
    }
    #[test]
    fn test_correlation_expected_rear_low_speed() {
        let corr = FrontRearSteeringCorrelation::new(-0.3, 0.2, 30.0, 0.01);
        let rear = corr.expected_rear(0.1, 0.0);
        assert!((rear - (-0.03)).abs() < 1e-10, "expected -0.03, got {rear}");
    }
    #[test]
    fn test_correlation_expected_rear_high_speed() {
        let corr = FrontRearSteeringCorrelation::new(-0.3, 0.2, 30.0, 0.01);
        let rear = corr.expected_rear(0.1, 30.0);
        assert!((rear - 0.02).abs() < 1e-10, "expected 0.02, got {rear}");
    }
    #[test]
    fn test_correlation_is_correlated_within_tolerance() {
        let corr = FrontRearSteeringCorrelation::new(0.0, 0.3, 30.0, 0.05);
        let front = 0.1;
        let speed = 30.0;
        let expected_rear = corr.expected_rear(front, speed);
        assert!(corr.is_correlated(front, expected_rear, speed));
    }
    #[test]
    fn test_correlation_is_not_correlated_outside_tolerance() {
        let corr = FrontRearSteeringCorrelation::new(0.0, 0.3, 30.0, 0.01);
        let front = 0.1;
        let speed = 30.0;
        assert!(!corr.is_correlated(front, 0.5, speed));
    }
    #[test]
    fn test_correlation_error_zero_when_perfect() {
        let corr = FrontRearSteeringCorrelation::new(0.0, 0.3, 30.0, 0.01);
        let front = 0.2;
        let speed = 30.0;
        let rear = corr.expected_rear(front, speed);
        let err = corr.correlation_error(front, rear, speed);
        assert!(err.abs() < 1e-12, "error should be zero: {err}");
    }
}
