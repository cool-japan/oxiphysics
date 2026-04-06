//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

// Types are re-exported via mod.rs; import in test modules via use super::*.

#[cfg(test)]
mod tests {

    use crate::Differential;
    use crate::DriveLayout;

    use crate::EngineCurve;
    use crate::Gearbox;

    use crate::drivetrain::BsfcMap;
    use crate::drivetrain::Clutch;

    use crate::drivetrain::DifferentialMode;
    use crate::drivetrain::DriveshaftCompliance;

    use crate::drivetrain::DrivetrainLegacy;
    use crate::drivetrain::Engine;
    use crate::drivetrain::EngineBraking;
    use crate::drivetrain::EngineCurveLegacy;
    use crate::drivetrain::EngineTorqueCurve;
    use crate::drivetrain::GearboxLegacy;
    use crate::drivetrain::GearboxManual;

    use crate::drivetrain::LimitedSlipDifferential;

    use crate::drivetrain::NvhDriveline;
    use crate::drivetrain::OpenDifferential;

    use crate::drivetrain::PolynomialTorqueCurve;

    use crate::drivetrain::TorqueConverter;

    use crate::drivetrain::TorsenDifferential;
    use oxiphysics_core::Real;
    #[test]
    fn test_engine_curve_torque_interpolates() {
        let curve =
            EngineCurve::new(vec![1000.0, 3000.0, 5000.0], vec![100.0, 300.0, 200.0]).unwrap();
        let t = curve.torque_at_rpm(2000.0);
        assert!(
            (t - 200.0).abs() < 1e-9,
            "interpolated torque {t}, expected 200"
        );
    }
    #[test]
    fn test_engine_available_torque_scales_with_throttle() {
        let mut engine = Engine::typical_4cylinder();
        engine.rpm = 3000.0;
        engine.throttle = 0.5;
        let half = engine.available_torque();
        engine.throttle = 1.0;
        let full = engine.available_torque();
        assert!((full - 2.0 * half).abs() < 1e-9, "full={full}, half={half}");
    }
    #[test]
    fn test_gearbox_total_ratio_changes_on_shift() {
        let mut gb = Gearbox::typical_6speed();
        gb.current_gear = 1;
        let r1 = gb.total_ratio();
        gb.current_gear = 2;
        let r2 = gb.total_ratio();
        assert!(
            r1 > r2,
            "1st gear ratio {r1} should exceed 2nd gear ratio {r2}"
        );
    }
    #[test]
    fn test_differential_split_torque_open_50_50() {
        let diff = Differential::open();
        let (l, r) = diff.split_torque(200.0, 0.0);
        assert!((l - 100.0).abs() < 1e-9, "left={l}");
        assert!((r - 100.0).abs() < 1e-9, "right={r}");
    }
    #[test]
    fn test_clutch_transmitted_torque_clamped_to_max() {
        let mut clutch = Clutch::new(300.0);
        clutch.engage();
        let t = clutch.transmitted_torque(200.0);
        assert!((t - 200.0).abs() < 1e-9, "t={t}");
        let t_capped = clutch.transmitted_torque(500.0);
        assert!((t_capped - 300.0).abs() < 1e-9, "t_capped={t_capped}");
    }
    #[test]
    fn test_engine_update_rpm_changes_with_load() {
        let mut engine = Engine::typical_4cylinder();
        engine.rpm = 2000.0;
        engine.throttle = 1.0;
        let initial = engine.rpm;
        engine.update(0.0, 0.1);
        assert!(
            engine.rpm > initial,
            "rpm should increase without load: {}",
            engine.rpm
        );
    }
    #[test]
    fn test_gearbox_shift_up_gear_increases() {
        let mut gb = Gearbox::typical_6speed();
        gb.current_gear = 1;
        gb.shift_up();
        assert_eq!(gb.current_gear, 2, "gear should be 2 after shift_up from 1");
    }
    #[test]
    fn test_engine_torque_lookup() {
        let engine = EngineCurveLegacy::new(vec![
            (1000.0, 100.0),
            (3000.0, 300.0),
            (5000.0, 400.0),
            (7000.0, 350.0),
        ]);
        let t = engine.torque_at_rpm(3000.0);
        assert!((t - 300.0).abs() < 1e-10, "t={t}");
        let t = engine.torque_at_rpm(2000.0);
        assert!((t - 200.0).abs() < 1e-10, "t={t}");
        let t = engine.torque_at_rpm(500.0);
        assert!((t - 100.0).abs() < 1e-10, "t={t}");
    }
    #[test]
    fn test_flat_engine_curve() {
        let engine = EngineCurveLegacy::flat(300.0, 800.0, 7000.0);
        let t = engine.torque_at_rpm(4000.0);
        assert!((t - 300.0).abs() < 1e-10);
    }
    #[test]
    fn test_gearbox_ratio() {
        let mut gearbox = GearboxLegacy::new(vec![3.5, 2.5, 1.8, 1.3, 1.0], 3.2, 3.7);
        gearbox.shift(1);
        let r = gearbox.total_ratio();
        assert!((r - 3.5 * 3.7).abs() < 1e-10);
        gearbox.shift(3);
        let r = gearbox.total_ratio();
        assert!((r - 1.8 * 3.7).abs() < 1e-10);
        gearbox.shift(0);
        assert!((gearbox.total_ratio()).abs() < 1e-10);
        gearbox.shift(-1);
        let r = gearbox.total_ratio();
        assert!((r - (-3.2 * 3.7)).abs() < 1e-10);
    }
    #[test]
    fn test_gearbox_shift_up_down() {
        let mut gearbox = GearboxLegacy::default();
        gearbox.shift(1);
        assert_eq!(gearbox.current_gear, 1);
        gearbox.shift_up();
        assert_eq!(gearbox.current_gear, 2);
        gearbox.shift_down();
        assert_eq!(gearbox.current_gear, 1);
        gearbox.shift_down();
        assert_eq!(gearbox.current_gear, 0);
        gearbox.shift_down();
        assert_eq!(gearbox.current_gear, -1);
        gearbox.shift_down();
        assert_eq!(gearbox.current_gear, -1);
    }
    #[test]
    fn test_differential_open() {
        let diff = DifferentialMode::Open;
        let (l, r) = diff.split_torque(100.0, 10.0, 20.0);
        assert!((l - 50.0).abs() < 1e-10);
        assert!((r - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_differential_limited_slip() {
        let diff = DifferentialMode::LimitedSlip { torque_bias: 3.0 };
        let (l, r) = diff.split_torque(100.0, 5.0, 20.0);
        assert!(l > r, "l={l}, r={r}");
        assert!((l + r - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_drivetrain_compute_torques() {
        let mut dt = DrivetrainLegacy::new();
        dt.layout = DriveLayout::RearWheelDrive;
        dt.gearbox.shift(1);
        let wheel_speeds = [0.0, 0.0, 10.0, 10.0];
        let torques = dt.compute_wheel_torques(1.0, 0.0, &wheel_speeds);
        assert!((torques[0]).abs() < 1e-10);
        assert!((torques[1]).abs() < 1e-10);
        assert!(torques[2] > 0.0, "rl={}", torques[2]);
        assert!(torques[3] > 0.0, "rr={}", torques[3]);
    }
    #[test]
    fn test_engine_torque_curve() {
        let engine = EngineCurveLegacy::default();
        let torque_idle = engine.torque_at_rpm(engine.idle_rpm);
        let torque_redline = engine.torque_at_rpm(engine.redline_rpm);
        let mut peak_torque: Real = 0.0;
        let mut peak_rpm: Real = 0.0;
        let steps = 1000;
        for i in 0..=steps {
            let rpm = engine.idle_rpm
                + (engine.redline_rpm - engine.idle_rpm) * i as Real / steps as Real;
            let t = engine.torque_at_rpm(rpm);
            if t > peak_torque {
                peak_torque = t;
                peak_rpm = rpm;
            }
        }
        assert!(peak_torque > torque_idle);
        assert!(peak_torque > torque_redline);
        let rpm_range = engine.redline_rpm - engine.idle_rpm;
        let peak_fraction = (peak_rpm - engine.idle_rpm) / rpm_range;
        assert!(peak_fraction > 0.1 && peak_fraction < 0.9);
    }
    #[test]
    fn test_engine_curve_idle() {
        let curve = EngineTorqueCurve::typical_na();
        let torque = curve.torque_at_rpm(curve.idle_rpm);
        assert!(
            torque > 0.0,
            "torque at idle should be positive, got {torque}"
        );
    }
    #[test]
    fn test_engine_curve_interpolation() {
        let curve = EngineTorqueCurve::typical_na();
        let midpoint_rpm = (curve.rpm_points[0] + curve.rpm_points[1]) * 0.5;
        let expected = (curve.torque_points[0] + curve.torque_points[1]) * 0.5;
        let got = curve.torque_at_rpm(midpoint_rpm);
        assert!((got - expected).abs() < 0.1);
    }
    #[test]
    fn test_engine_peak_torque_rpm() {
        let curve = EngineTorqueCurve::typical_na();
        let peak_rpm = curve.peak_torque_rpm();
        assert!((peak_rpm - 3000.0).abs() < 1.0);
    }
    #[test]
    fn test_clutch_disengaged() {
        let clutch = Clutch::new(500.0);
        let t = clutch.transmitted_torque(300.0);
        assert!(
            t.abs() < 1e-10,
            "disengaged clutch should transmit 0 N·m, got {t}"
        );
    }
    #[test]
    fn test_clutch_engaged() {
        let mut clutch = Clutch::new(500.0);
        clutch.engage();
        let t = clutch.transmitted_torque(300.0);
        assert!((t - 300.0).abs() < 1e-10, "got {t}");
        let t_capped = clutch.transmitted_torque(700.0);
        assert!((t_capped - 500.0).abs() < 1e-10, "got {t_capped}");
    }
    #[test]
    fn test_gearbox_manual_total_ratio() {
        let gb = GearboxManual::six_speed();
        let expected = 3.82 * 3.73;
        let got = gb.total_ratio();
        assert!((got - expected).abs() < 1e-9);
    }
    #[test]
    fn test_gearbox_manual_shift_up_down() {
        let mut gb = GearboxManual::six_speed();
        assert_eq!(gb.current_gear, 1);
        let shifted = gb.shift_up();
        assert!(shifted);
        assert_eq!(gb.current_gear, 2);
        gb.shift_up();
        assert_eq!(gb.current_gear, 3);
        let shifted = gb.shift_down();
        assert!(shifted);
        assert_eq!(gb.current_gear, 2);
        gb.shift_to(6);
        let at_max = gb.shift_up();
        assert!(!at_max);
        assert_eq!(gb.current_gear, 6);
    }
    #[test]
    fn test_gearbox_manual_wheel_rpm() {
        let gb = GearboxManual::six_speed();
        let engine_rpm = 3000.0;
        let total = gb.total_ratio();
        let expected_wheel_rpm = engine_rpm / total;
        let got = gb.wheel_rpm(engine_rpm);
        assert!((got - expected_wheel_rpm).abs() < 1e-9);
    }
    #[test]
    fn test_open_diff_equal_split() {
        let diff = OpenDifferential::new();
        let (left, right) = diff.torque_split(200.0);
        assert!((left - 100.0).abs() < 1e-10);
        assert!((right - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_lsd_biases_torque() {
        let mut lsd = LimitedSlipDifferential::new(50.0);
        lsd.left_wheel_speed = 5.0;
        lsd.right_wheel_speed = 15.0;
        let (left, right) = lsd.torque_split(200.0);
        assert!(
            left > right,
            "LSD should bias to slower (left) wheel, left={left:.2}, right={right:.2}"
        );
        assert!((left + right - 200.0).abs() < 1e-6);
    }
    #[test]
    fn test_poly_torque_curve_positive_at_mid_rpm() {
        let curve = PolynomialTorqueCurve::typical_na();
        let t = curve.torque_at_rpm(3000.0);
        assert!(t > 0.0, "torque at 3000 RPM should be positive, got {t}");
    }
    #[test]
    fn test_poly_torque_curve_clamped() {
        let curve = PolynomialTorqueCurve::typical_na();
        let t_below = curve.torque_at_rpm(0.0);
        let t_min = curve.torque_at_rpm(curve.min_rpm);
        assert!(
            (t_below - t_min).abs() < 1e-9,
            "below min_rpm should be clamped"
        );
    }
    #[test]
    fn test_poly_torque_curve_peak_rpm_in_range() {
        let curve = PolynomialTorqueCurve::typical_na();
        let peak = curve.peak_torque_rpm();
        assert!(
            peak >= curve.min_rpm && peak <= curve.max_rpm,
            "peak RPM {peak} should be in [{}, {}]",
            curve.min_rpm,
            curve.max_rpm
        );
    }
    #[test]
    fn test_torque_converter_lockup() {
        let tc = TorqueConverter::new(250.0, 0.65);
        let out = tc.output_torque(200.0, 1.0);
        assert!(
            (out - 200.0).abs() < 1e-9,
            "locked-up TC passes torque 1:1, got {out}"
        );
    }
    #[test]
    fn test_torque_converter_stall_multiplies_torque() {
        let tc = TorqueConverter::new(250.0, 0.65);
        let out = tc.output_torque(200.0, 0.0);
        assert!(out > 200.0, "at stall TC should multiply torque, got {out}");
    }
    #[test]
    fn test_torque_converter_capacity_positive() {
        let tc = TorqueConverter::new(250.0, 0.65);
        let cap = tc.capacity_factor(2000.0);
        assert!(cap > 0.0, "capacity factor should be positive, got {cap}");
    }
    #[test]
    fn test_engine_braking_zero_at_idle() {
        let eb = EngineBraking::new(20.0, 800.0);
        let t = eb.braking_torque(800.0);
        assert!(t.abs() < 1e-9, "no engine braking at idle, got {t}");
    }
    #[test]
    fn test_engine_braking_increases_with_rpm() {
        let eb = EngineBraking::new(20.0, 800.0);
        let t1 = eb.braking_torque(3000.0);
        let t2 = eb.braking_torque(5000.0);
        assert!(
            t2 > t1,
            "engine braking should increase with RPM, t1={t1}, t2={t2}"
        );
    }
    #[test]
    fn test_engine_braking_is_positive() {
        let eb = EngineBraking::new(20.0, 800.0);
        let t = eb.braking_torque(4000.0);
        assert!(
            t >= 0.0,
            "engine braking torque should be non-negative, got {t}"
        );
    }
    #[test]
    fn test_driveshaft_no_twist_at_equal_speeds() {
        let ds = DriveshaftCompliance::new(5000.0, 50.0);
        let t = ds.transmitted_torque(100.0, 100.0);
        assert!(t.abs() < 1e-9, "no twist angle → no torque, got {t}");
    }
    #[test]
    fn test_driveshaft_torque_positive_when_engine_faster() {
        let ds = DriveshaftCompliance::new(5000.0, 50.0);
        let t = ds.transmitted_torque(200.0, 100.0);
        assert!(
            t > 0.0,
            "engine faster → positive transmitted torque, got {t}"
        );
    }
    #[test]
    fn test_driveshaft_natural_frequency_positive() {
        let ds = DriveshaftCompliance::new(5000.0, 50.0);
        let freq = ds.natural_frequency(0.3);
        assert!(
            freq > 0.0,
            "natural frequency should be positive, got {freq}"
        );
    }
    #[test]
    fn test_torsen_equal_speeds_is_50_50() {
        let torsen = TorsenDifferential::new(4.0);
        let (l, r) = torsen.torque_split(200.0, 50.0, 50.0);
        assert!((l - 100.0).abs() < 1e-6, "equal speeds → 50/50, l={l}");
        assert!((r - 100.0).abs() < 1e-6, "equal speeds → 50/50, r={r}");
    }
    #[test]
    fn test_torsen_biases_toward_slower_wheel() {
        let torsen = TorsenDifferential::new(4.0);
        let (l, r) = torsen.torque_split(200.0, 20.0, 80.0);
        assert!(
            l > r,
            "Torsen should bias to slower (left) wheel, l={l:.2}, r={r:.2}"
        );
        assert!((l + r - 200.0).abs() < 1e-6, "total torque conserved");
    }
    #[test]
    fn test_bsfc_map_in_range() {
        let map = BsfcMap::default_gasoline();
        let bsfc = map.lookup(3000.0, 0.7);
        assert!(
            bsfc > 200.0 && bsfc < 500.0,
            "BSFC should be in realistic range, got {bsfc}"
        );
    }
    #[test]
    fn test_bsfc_low_load_higher_than_optimal() {
        let map = BsfcMap::default_gasoline();
        let bsfc_light = map.lookup(3000.0, 0.1);
        let bsfc_opt = map.lookup(3000.0, 0.65);
        assert!(
            bsfc_light > bsfc_opt,
            "light load should have higher BSFC, light={bsfc_light}, opt={bsfc_opt}"
        );
    }
    #[test]
    fn test_bsfc_fuel_flow_positive() {
        let map = BsfcMap::default_gasoline();
        let flow = map.fuel_flow_g_per_s(3000.0, 0.7, 150.0);
        assert!(flow > 0.0, "fuel flow should be positive, got {flow}");
    }
    #[test]
    fn test_nvh_resonance_detected() {
        let nvh = NvhDriveline::new(50.0, 5.0);
        let at_res = nvh.is_resonant(50.0, 2.0);
        assert!(
            at_res,
            "excitation at natural freq within bandwidth should be resonant"
        );
    }
    #[test]
    fn test_nvh_no_resonance_far_away() {
        let nvh = NvhDriveline::new(50.0, 5.0);
        let far = nvh.is_resonant(100.0, 2.0);
        assert!(!far, "excitation far from natural freq should not resonate");
    }
    #[test]
    fn test_nvh_amplitude_positive() {
        let nvh = NvhDriveline::new(50.0, 5.0);
        let amp = nvh.vibration_amplitude(50.0);
        assert!(
            amp > 0.0,
            "amplitude at resonance should be positive, got {amp}"
        );
    }
}
#[cfg(test)]
mod tests_drivetrain_new {

    use crate::Differential;
    use crate::Drivetrain;
    use crate::Gearbox;

    use crate::drivetrain::Engine;
    use crate::drivetrain::LaunchControl;

    #[test]
    fn test_lsd_equal_speeds_gives_50_50() {
        let diff = Differential::new(0.5);
        let (l, r) = diff.compute_torque_split_limited_slip(200.0, 10.0, 10.0);
        assert!(
            (l + r - 200.0).abs() < 1e-6,
            "total torque conserved: l={l}, r={r}"
        );
        assert!((l - r).abs() < 1e-3, "equal speeds → symmetric split");
    }
    #[test]
    fn test_lsd_biases_to_slower_wheel() {
        let diff = Differential::new(1.0);
        let (l, r) = diff.compute_torque_split_limited_slip(200.0, 5.0, 15.0);
        assert!(
            l > r,
            "slower (left) wheel should get more torque: l={l}, r={r}"
        );
        assert!((l + r - 200.0).abs() < 1e-6, "total torque conserved");
    }
    #[test]
    fn test_lsd_open_diff_remains_50_50() {
        let diff = Differential::open();
        let (l, r) = diff.compute_torque_split_limited_slip(100.0, 20.0, 5.0);
        assert!(
            (l - 50.0).abs() < 1e-6 && (r - 50.0).abs() < 1e-6,
            "open diff → 50/50: l={l}, r={r}"
        );
    }
    #[test]
    fn test_launch_control_limits_throttle_on_spin() {
        let mut lc = LaunchControl::default_racing();
        let out = lc.compute_launch_control(1.0, 4500.0, 0.3, 1.0, 0.01);
        assert!(out < 1.0, "excessive slip should reduce throttle: {out}");
    }
    #[test]
    fn test_launch_control_ramps_after_launch_speed() {
        let mut lc = LaunchControl::default_racing();
        lc.throttle_output = 0.5;
        let out = lc.compute_launch_control(1.0, 3000.0, 0.05, 5.0, 0.1);
        assert!(out > 0.5, "should ramp throttle after launch speed: {out}");
    }
    #[test]
    fn test_launch_control_output_bounded() {
        let mut lc = LaunchControl::default_racing();
        lc.throttle_output = 0.0;
        for _ in 0..200 {
            let out = lc.compute_launch_control(1.0, 4500.0, 0.05, 0.5, 0.01);
            assert!((0.0..=1.0).contains(&out), "throttle out of bounds: {out}");
        }
    }
    #[test]
    fn test_exhaust_temp_increases_with_throttle() {
        let mut eng = Engine::typical_4cylinder();
        eng.rpm = 3000.0;
        eng.throttle = 0.2;
        let t_low = eng.compute_exhaust_temperature();
        eng.throttle = 0.9;
        let t_high = eng.compute_exhaust_temperature();
        assert!(
            t_high > t_low,
            "higher throttle → higher EGT: {t_low} vs {t_high}"
        );
    }
    #[test]
    fn test_exhaust_temp_in_realistic_range() {
        let mut eng = Engine::typical_4cylinder();
        eng.throttle = 1.0;
        eng.rpm = 6000.0;
        let t = eng.compute_exhaust_temperature();
        assert!(t > 600.0 && t < 1300.0, "EGT should be 600–1300 K, got {t}");
    }
    #[test]
    fn test_exhaust_temp_minimum_at_idle() {
        let mut eng = Engine::typical_4cylinder();
        eng.throttle = 0.0;
        eng.rpm = eng.idle_rpm;
        let t = eng.compute_exhaust_temperature();
        let t_max = {
            eng.throttle = 1.0;
            eng.rpm = eng.max_rpm;
            eng.compute_exhaust_temperature()
        };
        assert!(
            t < t_max,
            "idle should have lower EGT than WOT: idle={t}, wot={t_max}"
        );
    }
    #[test]
    fn test_powertrain_inertia_exceeds_engine_inertia() {
        let dt = Drivetrain::new(Engine::typical_4cylinder(), Gearbox::typical_6speed(), 0.32);
        let i_total = dt.compute_powertrain_inertia(1.5);
        assert!(
            i_total > dt.engine.inertia,
            "total inertia should exceed engine alone: engine={}, total={i_total}",
            dt.engine.inertia
        );
    }
    #[test]
    fn test_powertrain_inertia_higher_in_lower_gear() {
        let mut dt1 = Drivetrain::new(Engine::typical_4cylinder(), Gearbox::typical_6speed(), 0.32);
        dt1.gearbox.current_gear = 1;
        let i1 = dt1.compute_powertrain_inertia(1.5);
        let mut dt6 = Drivetrain::new(Engine::typical_4cylinder(), Gearbox::typical_6speed(), 0.32);
        dt6.gearbox.current_gear = 6;
        let i6 = dt6.compute_powertrain_inertia(1.5);
        assert!(
            i1 < i6,
            "lower gear (higher ratio) → less wheel inertia referred back → lower total: gear1={i1}, gear6={i6}"
        );
    }
}
#[cfg(test)]
mod tests_drivetrain_extended {

    use crate::drivetrain::AutomaticTransmission;
    use crate::drivetrain::AwdCenterDifferential;
    use crate::drivetrain::AwdCenterMode;
    use crate::drivetrain::CvtTransmission;
    use crate::drivetrain::DriveshaftNonlinear;

    use crate::drivetrain::LockupState;
    use crate::drivetrain::PlanetaryGearSet;
    use crate::drivetrain::ShiftMode;
    use crate::drivetrain::TorqueConverterLockup;

    #[test]
    fn test_planetary_basic_ratio() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        assert!(
            (pg.basic_ratio() - (-3.0)).abs() < 1e-9,
            "basic ratio = -3.0, got {}",
            pg.basic_ratio()
        );
    }
    #[test]
    fn test_planetary_sun_to_ring_ratio() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        assert!((pg.ratio_sun_to_ring_fixed_carrier() - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_planetary_sun_to_carrier_ratio() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        let expected = 1.0 + 30.0 / 90.0;
        assert!((pg.ratio_sun_to_carrier_fixed_ring() - expected).abs() < 1e-9);
    }
    #[test]
    fn test_planetary_ring_speed_willis() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        let ring = pg.ring_speed(100.0, 0.0);
        assert!(
            (ring - 300.0).abs() < 1e-9,
            "ring_speed={ring}, expected 300"
        );
    }
    #[test]
    fn test_planetary_ring_torque_from_sun() {
        let pg = PlanetaryGearSet::new(30, 90, 0.97);
        let t = pg.ring_torque_from_sun(100.0);
        assert!((t - 291.0).abs() < 1e-9, "ring torque={t}");
    }
    #[test]
    fn test_planetary_min_planet_count_divides_sum() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        let n = pg.min_planet_count();
        assert!(
            (pg.z_sun + pg.z_ring).is_multiple_of(n),
            "sum must divide evenly by n={n}"
        );
    }
    #[test]
    fn test_planetary_sun_speed_round_trip() {
        let pg = PlanetaryGearSet::new(30, 90, 1.0);
        let sun_orig = 200.0;
        let carrier = 50.0;
        let ring = pg.ring_speed(sun_orig, carrier);
        let sun_back = pg.sun_speed(ring, carrier);
        assert!(
            (sun_back - sun_orig).abs() < 1e-9,
            "round-trip sun speed: {sun_back}"
        );
    }
    #[test]
    fn test_tc_lockup_initial_state_open() {
        let tc = TorqueConverterLockup::typical_passenger_car();
        assert_eq!(tc.state, LockupState::Open);
        assert!(!tc.is_locked());
    }
    #[test]
    fn test_tc_lockup_transitions_to_engaging() {
        let mut tc = TorqueConverterLockup::typical_passenger_car();
        let sr = tc.update(1000.0, 900.0, 0.01);
        assert!(sr > 0.0);
        assert_ne!(tc.state, LockupState::Open, "should transition from Open");
    }
    #[test]
    fn test_tc_lockup_fully_locks_after_engage_time() {
        let mut tc = TorqueConverterLockup::typical_passenger_car();
        tc.engage_time = 0.1;
        for _ in 0..50 {
            tc.update(1000.0, 900.0, 0.01);
        }
        assert!(tc.is_locked(), "should be locked after engage_time");
    }
    #[test]
    fn test_tc_lockup_reopens_on_low_sr() {
        let mut tc = TorqueConverterLockup::typical_passenger_car();
        tc.engage_time = 0.1;
        for _ in 0..50 {
            tc.update(1000.0, 900.0, 0.01);
        }
        assert!(tc.is_locked());
        tc.update(1000.0, 500.0, 0.1);
        assert_eq!(tc.state, LockupState::Open, "should open when SR drops");
    }
    #[test]
    fn test_tc_effective_output_open_mode() {
        let tc = TorqueConverterLockup::typical_passenger_car();
        let out = tc.effective_output_torque(100.0, 0.0);
        assert!(
            out > 100.0,
            "at stall converter should amplify torque: {out}"
        );
    }
    #[test]
    fn test_at_initial_gear_is_1() {
        let at = AutomaticTransmission::zf_6hp();
        assert_eq!(at.current_gear, 1);
    }
    #[test]
    fn test_at_total_ratio_first_gear() {
        let at = AutomaticTransmission::zf_6hp();
        let expected = 4.17 * 3.46;
        assert!(
            (at.total_ratio() - expected).abs() < 1e-6,
            "1st gear ratio: {}",
            at.total_ratio()
        );
    }
    #[test]
    fn test_at_output_torque_positive() {
        let at = AutomaticTransmission::zf_6hp();
        let t = at.output_torque(100.0);
        assert!(t > 0.0, "output torque must be positive: {t}");
    }
    #[test]
    fn test_at_upshift_at_high_rpm() {
        let mut at = AutomaticTransmission::zf_6hp();
        at.shift_mode = ShiftMode::Economy;
        let shifted = at.step_shift_strategy(2500.0, 0.5, 0.1);
        assert!(shifted, "should request upshift at high RPM");
    }
    #[test]
    fn test_at_no_shift_below_threshold() {
        let mut at = AutomaticTransmission::zf_6hp();
        at.shift_mode = ShiftMode::Economy;
        let shifted = at.step_shift_strategy(800.0, 0.3, 0.1);
        assert!(!shifted, "should not shift at low RPM");
    }
    #[test]
    fn test_at_manual_upshift() {
        let mut at = AutomaticTransmission::zf_6hp();
        at.shift_mode = ShiftMode::Manual;
        at.manual_upshift();
        assert!(at.is_shifting(), "manual upshift should start shift");
    }
    #[test]
    fn test_at_manual_downshift() {
        let mut at = AutomaticTransmission::zf_6hp();
        at.current_gear = 3;
        at.shift_mode = ShiftMode::Manual;
        at.manual_downshift();
        assert!(at.is_shifting());
        for _ in 0..10 {
            at.step_shift_strategy(3000.0, 0.5, 0.05);
        }
        assert_eq!(at.current_gear, 2, "should have downshifted to gear 2");
    }
    #[test]
    fn test_at_torque_reduction_during_shift() {
        let mut at = AutomaticTransmission::zf_6hp();
        at.step_shift_strategy(2500.0, 0.5, 0.1);
        let factor = at.torque_reduction_factor();
        assert!(factor < 1.0, "torque reduction during shift: {factor}");
    }
    #[test]
    fn test_at_sport_mode_holds_gear_longer() {
        let mut at_eco = AutomaticTransmission::zf_6hp();
        at_eco.shift_mode = ShiftMode::Economy;
        let shifted_eco = at_eco.step_shift_strategy(3000.0, 0.5, 0.01);
        let mut at_sport = AutomaticTransmission::zf_6hp();
        at_sport.shift_mode = ShiftMode::Sport;
        let shifted_sport = at_sport.step_shift_strategy(3000.0, 0.5, 0.01);
        assert!(!shifted_sport, "sport mode should hold gear at 3000 RPM");
        let _ = shifted_eco;
    }
    #[test]
    fn test_cvt_initial_ratio_is_max() {
        let cvt = CvtTransmission::jatco_jf011e();
        assert!(
            (cvt.current_ratio - cvt.ratio_max).abs() < 1e-9,
            "starts in low gear"
        );
    }
    #[test]
    fn test_cvt_total_ratio_positive() {
        let cvt = CvtTransmission::jatco_jf011e();
        assert!(cvt.total_ratio() > 0.0);
    }
    #[test]
    fn test_cvt_output_torque_increases_with_ratio() {
        let mut cvt = CvtTransmission::jatco_jf011e();
        cvt.current_ratio = cvt.ratio_max;
        let t_low = cvt.output_torque(100.0);
        cvt.current_ratio = cvt.ratio_min;
        let t_high_gear = cvt.output_torque(100.0);
        assert!(
            t_low > t_high_gear,
            "lower ratio gear = less torque multiplication"
        );
    }
    #[test]
    fn test_cvt_step_ratio_moves_toward_target() {
        let mut cvt = CvtTransmission::jatco_jf011e();
        cvt.current_ratio = cvt.ratio_max;
        let new_r = cvt.step_ratio(cvt.ratio_min, 1.0);
        assert!(
            new_r < cvt.ratio_max,
            "ratio should move toward target: {new_r}"
        );
    }
    #[test]
    fn test_cvt_ideal_ratio_clamped() {
        let cvt = CvtTransmission::jatco_jf011e();
        let r = cvt.ideal_ratio_for_wheel_speed(200.0);
        assert!(
            r >= cvt.ratio_min && r <= cvt.ratio_max,
            "ratio should be clamped: {r}"
        );
    }
    #[test]
    fn test_cvt_ratio_coverage_above_1() {
        let cvt = CvtTransmission::jatco_jf011e();
        assert!(cvt.ratio_coverage() > 1.0, "coverage must exceed 1");
    }
    #[test]
    fn test_cvt_zero_wheel_speed_returns_ratio_max() {
        let cvt = CvtTransmission::jatco_jf011e();
        let r = cvt.ideal_ratio_for_wheel_speed(0.0);
        assert!((r - cvt.ratio_max).abs() < 1e-9);
    }
    #[test]
    fn test_nonlinear_driveshaft_no_torque_in_lash_band() {
        let ds = DriveshaftNonlinear::new(5000.0, 50.0, 0.05, 1000.0);
        let t = ds.transmitted_torque(0.0, 0.0);
        assert!(t.abs() < 1e-9, "no torque in lash band: {t}");
    }
    #[test]
    fn test_nonlinear_driveshaft_torque_outside_lash() {
        let mut ds = DriveshaftNonlinear::new(5000.0, 50.0, 0.01, 1000.0);
        ds.deflection = 0.1;
        let t = ds.transmitted_torque(100.0, 0.0);
        assert!(
            t > 0.0,
            "deflection outside lash should produce torque: {t}"
        );
    }
    #[test]
    fn test_nonlinear_driveshaft_torque_limited() {
        let mut ds = DriveshaftNonlinear::new(50_000.0, 500.0, 0.0, 500.0);
        ds.deflection = 1.0;
        let t = ds.transmitted_torque(200.0, 0.0);
        assert!(t.abs() <= 500.0, "torque must be clamped to limit: {t}");
    }
    #[test]
    fn test_nonlinear_driveshaft_step_advances_deflection() {
        let mut ds = DriveshaftNonlinear::new(5000.0, 50.0, 0.01, 1000.0);
        ds.step(100.0, 50.0, 0.1);
        assert!(
            ds.deflection > 0.0,
            "deflection should increase: {}",
            ds.deflection
        );
    }
    #[test]
    fn test_awd_fixed_split_respects_fraction() {
        let awd = AwdCenterDifferential::fixed_40_60();
        let (front, rear) = awd.split_torque(1000.0);
        assert!((front - 400.0).abs() < 1e-9, "front={front}");
        assert!((rear - 600.0).abs() < 1e-9, "rear={rear}");
    }
    #[test]
    fn test_awd_torque_conserved_viscous() {
        let awd = AwdCenterDifferential::torsen_centre();
        let (front, rear) = awd.split_torque(500.0);
        assert!(
            (front + rear - 500.0).abs() < 1e-9,
            "total torque must be conserved"
        );
    }
    #[test]
    fn test_awd_active_mode_follows_command() {
        let awd = AwdCenterDifferential {
            mode: AwdCenterMode::Active,
            commanded_front: 0.7,
            front_speed: 0.0,
            rear_speed: 0.0,
        };
        let (front, rear) = awd.split_torque(100.0);
        assert!((front - 70.0).abs() < 1e-9, "front={front}");
        assert!((rear - 30.0).abs() < 1e-9, "rear={rear}");
    }
    #[test]
    fn test_awd_viscous_biases_on_speed_diff() {
        let mut awd = AwdCenterDifferential::torsen_centre();
        awd.front_speed = 50.0;
        awd.rear_speed = 10.0;
        let (front_slip, rear_slip) = awd.split_torque(1000.0);
        awd.front_speed = 10.0;
        awd.rear_speed = 10.0;
        let (front_equal, _rear_equal) = awd.split_torque(1000.0);
        assert!(
            front_slip < front_equal,
            "front torque should reduce when front slips"
        );
        assert!(front_slip + rear_slip - 1000.0 < 1e-6, "total conserved");
    }
}
