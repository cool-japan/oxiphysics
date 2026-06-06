//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod extra_tests {
    use super::super::*;
    use crate::traits::types::AngularMotorData;
    use crate::traits::types::AngularMotorMode;
    use crate::traits::types::ConeLimitData;
    use crate::traits::types::ConstraintResidualTracker;
    use crate::traits::types::DistanceLimitData;
    use crate::traits::types::ErrorNorm;
    use crate::traits::types::LagrangeMultiplier;
    use crate::traits::types::SpringConstraintData;
    use crate::traits::types::XpbdParams;
    #[test]
    fn error_norm_l1() {
        let r = vec![1.0, -2.0, 3.0];
        assert!((compute_error_norm(&r, ErrorNorm::L1) - 6.0).abs() < 1e-12);
    }
    #[test]
    fn error_norm_l2() {
        let r = vec![3.0, 4.0];
        assert!((compute_error_norm(&r, ErrorNorm::L2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn error_norm_linf() {
        let r = vec![1.0, -5.0, 3.0];
        assert!((compute_error_norm(&r, ErrorNorm::LInfinity) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn error_norm_empty() {
        assert_eq!(compute_error_norm(&[], ErrorNorm::L2), 0.0);
    }
    #[test]
    fn error_norm_weighted_l2_falls_back_to_l2() {
        let r = vec![3.0, 4.0];
        assert!((compute_error_norm(&r, ErrorNorm::WeightedL2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn weighted_l2_norm_correct() {
        let r = vec![2.0, 3.0];
        let w = vec![1.0, 2.0];
        let expected = 40.0_f64.sqrt();
        assert!((compute_weighted_l2_norm(&r, &w) - expected).abs() < 1e-10);
    }
    #[test]
    fn normalize_residuals_unit_norm() {
        let r = vec![3.0, 4.0];
        let n = normalize_residuals(&r);
        let l2: f64 = n.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((l2 - 1.0).abs() < 1e-10);
    }
    #[test]
    fn normalize_residuals_zero_stays_zero() {
        let r = vec![0.0, 0.0];
        let n = normalize_residuals(&r);
        assert!(n.iter().all(|x| x.abs() < 1e-15));
    }
    #[test]
    fn lagrange_bilateral_unbounded() {
        let mut lm = LagrangeMultiplier::bilateral();
        let actual = lm.apply_delta(2.0);
        assert!((actual - 2.0).abs() < 1e-12);
        assert!((lm.lambda - 2.0).abs() < 1e-12);
    }
    #[test]
    fn lagrange_unilateral_clamps_below_zero() {
        let mut lm = LagrangeMultiplier::unilateral();
        let actual = lm.apply_delta(-1.0);
        assert!(
            actual.abs() < 1e-12,
            "unilateral cannot go negative: actual={actual}"
        );
        assert!(lm.lambda.abs() < 1e-12);
    }
    #[test]
    fn lagrange_reset() {
        let mut lm = LagrangeMultiplier::bilateral();
        lm.apply_delta(5.0);
        lm.reset();
        assert_eq!(lm.lambda, 0.0);
    }
    #[test]
    fn lagrange_warm_start() {
        let mut lm = LagrangeMultiplier::bilateral();
        lm.warm_start(10.0, 0.8);
        assert!((lm.lambda - 8.0).abs() < 1e-12);
    }
    #[test]
    fn lagrange_with_bounds_clamping() {
        let mut lm = LagrangeMultiplier::with_bounds(-1.0, 1.0);
        lm.apply_delta(5.0);
        assert!(
            (lm.lambda - 1.0).abs() < 1e-12,
            "clamped to upper bound: {}",
            lm.lambda
        );
    }
    #[test]
    fn lagrange_default_is_bilateral() {
        let lm = LagrangeMultiplier::default();
        assert_eq!(lm.lower_bound, f64::NEG_INFINITY);
        assert_eq!(lm.upper_bound, f64::INFINITY);
    }
    #[test]
    fn lagrange_apply_delta_returns_actual_change() {
        let mut lm = LagrangeMultiplier::unilateral();
        let d1 = lm.apply_delta(-3.0);
        assert_eq!(d1, 0.0);
        let d2 = lm.apply_delta(2.0);
        assert!((d2 - 2.0).abs() < 1e-12, "d2={d2}");
    }
    #[test]
    fn xpbd_alpha_tilde_correct() {
        let p = XpbdParams::new(1e-4, 0.01);
        let at = p.alpha_tilde(0.01);
        assert!((at - 1.0).abs() < 1e-10);
    }
    #[test]
    fn xpbd_beta_tilde_correct() {
        let p = XpbdParams::new(0.0, 0.1);
        let bt = p.beta_tilde(0.01);
        assert!((bt - 10.0).abs() < 1e-10);
    }
    #[test]
    fn xpbd_delta_lambda_zero_compliance() {
        let p = XpbdParams::rigid();
        let dl = p.delta_lambda(0.1, 0.0, 2.0, 0.0, 0.01);
        assert!((dl + 0.05).abs() < 1e-10);
    }
    #[test]
    fn xpbd_delta_lambda_soft() {
        let p = XpbdParams::new(1e-4, 0.0);
        let dt = 0.01;
        let at = p.alpha_tilde(dt);
        let dl = p.delta_lambda(0.1, 0.0, 2.0, 0.0, dt);
        let expected = -(0.1 + at * 0.0) / (2.0 + at);
        assert!((dl - expected).abs() < 1e-10);
    }
    #[test]
    fn xpbd_delta_lambda_zero_dt() {
        let p = XpbdParams::new(1e-4, 0.01);
        let dl = p.delta_lambda(1.0, 0.0, 1.0, 0.0, 0.0);
        assert!((dl + 1.0).abs() < 1e-10);
    }
    #[test]
    fn xpbd_delta_lambda_degenerate_mass() {
        let p = XpbdParams::rigid();
        let dl = p.delta_lambda(1.0, 0.0, 0.0, 0.0, 0.01);
        assert!(dl.abs() < 1e-10, "zero mass → zero delta");
    }
    #[test]
    fn xpbd_params_default_is_rigid() {
        let p = XpbdParams::default();
        assert_eq!(p.compliance, 0.0);
        assert_eq!(p.damping, 0.0);
    }
    #[test]
    fn spring_force_at_rest() {
        let s = SpringConstraintData::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 1.0);
        assert!((s.spring_force(1.0)).abs() < 1e-12);
    }
    #[test]
    fn spring_force_stretched() {
        let s = SpringConstraintData::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 0.0);
        assert!((s.spring_force(1.5) - 50.0).abs() < 1e-10);
    }
    #[test]
    fn spring_force_compressed() {
        let s = SpringConstraintData::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 100.0, 0.0);
        assert!((s.spring_force(0.5) + 50.0).abs() < 1e-10);
    }
    #[test]
    fn spring_active_flag() {
        let s = SpringConstraintData::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 10.0, 0.0);
        assert!(s.active);
    }
    #[test]
    fn distance_limit_within_bounds() {
        let dl = DistanceLimitData::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 2.0);
        assert!(
            (dl.violation(1.0)).abs() < 1e-12,
            "within bounds → no violation"
        );
    }
    #[test]
    fn distance_limit_below_lower() {
        let dl = DistanceLimitData::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 2.0);
        let v = dl.violation(0.2);
        assert!(v > 0.0, "below lower → positive violation: {v}");
    }
    #[test]
    fn distance_limit_above_upper() {
        let dl = DistanceLimitData::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, 2.0);
        let v = dl.violation(3.0);
        assert!(v < 0.0, "above upper → negative violation: {v}");
    }
    #[test]
    fn distance_limit_lower_active() {
        let dl = DistanceLimitData::new([0.0; 3], [0.0, 0.0, 0.0], 1.0, 5.0);
        assert!(dl.lower_active(0.5));
        assert!(!dl.lower_active(2.0));
    }
    #[test]
    fn distance_limit_upper_active() {
        let dl = DistanceLimitData::new([0.0; 3], [0.0, 0.0, 0.0], 0.0, 3.0);
        assert!(dl.upper_active(4.0));
        assert!(!dl.upper_active(2.0));
    }
    #[test]
    fn angular_motor_off_produces_zero_torque() {
        let mut m = AngularMotorData::new(100.0);
        let t = m.compute_torque(0.0, 0.0, 0.01);
        assert_eq!(t, 0.0);
    }
    #[test]
    fn angular_motor_velocity_mode_drives_toward_target() {
        let mut m = AngularMotorData::new(1000.0);
        m.mode = AngularMotorMode::Velocity;
        m.target_velocity = 5.0;
        let t = m.compute_torque(0.0, 0.0, 0.01);
        assert!(t > 0.0, "motor should apply positive torque: {t}");
    }
    #[test]
    fn angular_motor_position_mode_drives_toward_target() {
        let mut m = AngularMotorData::new(1000.0);
        m.mode = AngularMotorMode::Position;
        m.target_angle = std::f64::consts::PI / 4.0;
        let t = m.compute_torque(0.0, 0.0, 0.01);
        assert!(
            t > 0.0,
            "motor should apply positive torque toward target: {t}"
        );
    }
    #[test]
    fn angular_motor_torque_clamped_to_max() {
        let mut m = AngularMotorData::new(1.0);
        m.mode = AngularMotorMode::Position;
        m.target_angle = 100.0;
        m.kp = 10000.0;
        let t = m.compute_torque(0.0, 0.0, 0.01);
        assert!((t - 1.0).abs() < 1e-12, "torque clamped to max: {t}");
    }
    #[test]
    fn angular_motor_at_target() {
        let m = AngularMotorData::new(100.0);
        assert!(
            m.at_target(0.0, 1e-3),
            "at target angle 0 with tolerance 1e-3"
        );
        assert!(!m.at_target(1.0, 1e-3), "not at target angle 1.0");
    }
    #[test]
    fn angular_motor_default_mode_is_off() {
        let m = AngularMotorData::default();
        assert_eq!(m.mode, AngularMotorMode::Off);
    }
    #[test]
    fn cone_limit_aligned_axis_zero_angle() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], 0.5);
        let angle = cl.current_angle([0.0, 1.0, 0.0]);
        assert!(angle.abs() < 1e-10, "aligned axis → zero angle: {angle}");
    }
    #[test]
    fn cone_limit_perpendicular_axis_pi_half() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], 1.0);
        let angle = cl.current_angle([1.0, 0.0, 0.0]);
        assert!(
            (angle - std::f64::consts::PI / 2.0).abs() < 1e-10,
            "angle={angle}"
        );
    }
    #[test]
    fn cone_limit_not_violated_within() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], std::f64::consts::PI / 4.0);
        assert!(
            !cl.is_violated([0.0, 1.0, 0.0]),
            "aligned axis should not violate"
        );
    }
    #[test]
    fn cone_limit_violated_outside() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], 0.1);
        assert!(
            cl.is_violated([1.0, 0.0, 0.0]),
            "perpendicular axis should violate 0.1 rad cone"
        );
    }
    #[test]
    fn cone_limit_residual_zero_when_within() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], std::f64::consts::PI / 4.0);
        assert_eq!(cl.residual([0.0, 1.0, 0.0]), 0.0);
    }
    #[test]
    fn cone_limit_residual_positive_when_violated() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], 0.1);
        let r = cl.residual([1.0, 0.0, 0.0]);
        assert!(r > 0.0, "residual should be positive when violated: {r}");
    }
    #[test]
    fn cone_limit_update_active() {
        let mut cl = ConeLimitData::new([0.0, 1.0, 0.0], 0.1);
        cl.update_active([1.0, 0.0, 0.0]);
        assert!(cl.active, "should be active when violated");
    }
    #[test]
    fn cone_limit_degenerate_zero_axis() {
        let cl = ConeLimitData::new([0.0, 1.0, 0.0], 0.5);
        let angle = cl.current_angle([0.0, 0.0, 0.0]);
        assert_eq!(angle, 0.0, "degenerate zero axis should return 0");
    }
    #[test]
    fn tracker_initial_zero() {
        let t = ConstraintResidualTracker::new(3, ErrorNorm::L2, 1e-6);
        assert_eq!(t.error_norm(), 0.0);
        assert_eq!(t.iteration_count, 0);
    }
    #[test]
    fn tracker_update_and_norm() {
        let mut t = ConstraintResidualTracker::new(2, ErrorNorm::L2, 1e-6);
        t.update(&[3.0, 4.0]);
        assert!((t.error_norm() - 5.0).abs() < 1e-10);
        assert_eq!(t.iteration_count, 1);
    }
    #[test]
    fn tracker_converged_after_small_residual() {
        let mut t = ConstraintResidualTracker::new(2, ErrorNorm::L2, 1e-4);
        t.update(&[1e-5, 1e-5]);
        assert!(t.has_converged());
    }
    #[test]
    fn tracker_not_converged_with_large_residual() {
        let mut t = ConstraintResidualTracker::new(2, ErrorNorm::L2, 1e-6);
        t.update(&[1.0, 2.0]);
        assert!(!t.has_converged());
    }
    #[test]
    fn tracker_reset_clears_state() {
        let mut t = ConstraintResidualTracker::new(2, ErrorNorm::L2, 1e-6);
        t.update(&[5.0, 5.0]);
        t.reset();
        assert_eq!(t.error_norm(), 0.0);
        assert_eq!(t.iteration_count, 0);
    }
    #[test]
    fn angular_motor_mode_default_is_off() {
        let mode = AngularMotorMode::default();
        assert_eq!(mode, AngularMotorMode::Off);
    }
    #[test]
    fn angular_motor_mode_eq() {
        assert_eq!(AngularMotorMode::Velocity, AngularMotorMode::Velocity);
        assert_ne!(AngularMotorMode::Velocity, AngularMotorMode::Position);
    }
    #[test]
    fn error_norm_l1_single_value() {
        assert!((compute_error_norm(&[-3.0], ErrorNorm::L1) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn error_norm_linf_single_value() {
        assert!((compute_error_norm(&[7.0], ErrorNorm::LInfinity) - 7.0).abs() < 1e-12);
    }
    #[test]
    fn lagrange_bilateral_accumulates() {
        let mut lm = LagrangeMultiplier::bilateral();
        lm.apply_delta(1.0);
        lm.apply_delta(2.0);
        assert!((lm.lambda - 3.0).abs() < 1e-12);
    }
    #[test]
    fn distance_limit_violation_at_lower_boundary() {
        let dl = DistanceLimitData::new([0.0; 3], [0.0; 3], 1.0, 3.0);
        assert!((dl.violation(1.0)).abs() < 1e-12);
    }
}
