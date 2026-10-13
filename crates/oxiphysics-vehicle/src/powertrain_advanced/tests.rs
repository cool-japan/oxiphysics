// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the advanced powertrain models.

use super::*;

// CombustionEngine tests

#[test]
fn test_engine_create_default_values() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    assert_eq!(eng.cylinders, 4);
    assert_eq!(eng.displacement, 2.0e-3);
    assert!((eng.compression_ratio - 10.5).abs() < 1e-9);
}

#[test]
fn test_volumetric_efficiency_at_full_throttle() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let ve = eng.volumetric_efficiency(3500.0, 1.0);
    assert!(ve > 0.5 && ve <= 1.0);
}

#[test]
fn test_volumetric_efficiency_decreases_at_high_rpm() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let ve_mid = eng.volumetric_efficiency(3500.0, 1.0);
    let ve_high = eng.volumetric_efficiency(7000.0, 1.0);
    assert!(ve_mid > ve_high);
}

#[test]
fn test_bmep_positive() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let bmep = eng.bmep(200.0);
    assert!(bmep > 0.0);
}

#[test]
fn test_wiebe_fraction_zero_before_start() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let x = eng.wiebe_burned_fraction(-20.0);
    assert_eq!(x, 0.0);
}

#[test]
fn test_wiebe_fraction_increases_monotonically() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let x1 = eng.wiebe_burned_fraction(0.0);
    let x2 = eng.wiebe_burned_fraction(20.0);
    let x3 = eng.wiebe_burned_fraction(50.0);
    assert!(x1 <= x2 && x2 <= x3);
}

#[test]
fn test_engine_torque_at_throttle_zero() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let t = eng.torque(3500.0, 0.0);
    assert_eq!(t, 0.0);
}

#[test]
fn test_engine_power_positive_at_full_throttle() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let p = eng.power(3000.0, 1.0);
    assert!(p > 0.0);
}

#[test]
fn test_bsfc_positive() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let bsfc = eng.bsfc(2000.0, 0.8);
    assert!(bsfc > 0.0 && bsfc < 1000.0);
}

#[test]
fn test_fuel_flow_rate_positive() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let ffr = eng.fuel_flow_rate(3000.0, 1.0);
    assert!(ffr > 0.0);
}

#[test]
fn test_knock_probability_clipped_to_unity() {
    let eng = CombustionEngine::new(2.0e-3, 4, 14.0); // high CR
    let kp = eng.knock_probability(7000.0, 1.0);
    assert!((0.0..=1.0).contains(&kp));
}

#[test]
fn test_engine_step_increases_crank_angle() {
    let mut eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    eng.set_throttle(0.5);
    eng.rpm = 3000.0;
    let before = eng.crank_angle;
    eng.step(0.01, 50.0);
    assert!(eng.crank_angle != before || eng.rpm > eng.idle_rpm);
}

// FuelInjection tests

#[test]
fn test_injection_port_create() {
    let fi = FuelInjection::new_port_injection(4);
    assert_eq!(fi.spec.count, 4);
    assert!(fi.spec.rail_pressure > 0.0);
}

#[test]
fn test_injection_gdi_higher_pressure() {
    let port = FuelInjection::new_port_injection(4);
    let gdi = FuelInjection::new_gdi(4);
    assert!(gdi.spec.rail_pressure > port.spec.rail_pressure);
}

#[test]
fn test_spray_penetration_zero_at_zero_time() {
    let fi = FuelInjection::new_gdi(4);
    let pen = fi.spray_penetration(0.0, 1.2);
    assert_eq!(pen, 0.0);
}

#[test]
fn test_spray_penetration_increases_with_time() {
    let fi = FuelInjection::new_gdi(4);
    let p1 = fi.spray_penetration(1e-4, 1.2);
    let p2 = fi.spray_penetration(1e-3, 1.2);
    assert!(p2 > p1);
}

#[test]
fn test_compute_duration_positive() {
    let fi = FuelInjection::new_port_injection(4);
    let dur = fi.compute_duration(5e-5);
    assert!(dur > 0.0);
}

#[test]
fn test_set_lambda_sets_fields() {
    let mut fi = FuelInjection::new_port_injection(4);
    fi.set_lambda_and_compute(0.03, 1.0);
    assert!((fi.lambda - 1.0).abs() < 1e-9);
    assert!(fi.fuel_mass_cycle > 0.0);
}

// ExhaustSystem tests

#[test]
fn test_exhaust_backpressure_increases_with_rpm() {
    let exh = ExhaustSystem::new();
    let bp1 = exh.back_pressure(1000.0);
    let bp2 = exh.back_pressure(6000.0);
    assert!(bp2 > bp1);
}

#[test]
fn test_exhaust_temp_increases_with_load() {
    let exh = ExhaustSystem::new();
    let t1 = exh.exhaust_temp(3000.0, 0.3);
    let t2 = exh.exhaust_temp(3000.0, 0.9);
    assert!(t2 > t1);
}

// TurboCharger tests

#[test]
fn test_turbo_create_default() {
    let tc = TurboCharger::new();
    assert_eq!(tc.turbo_rpm, 0.0);
    assert_eq!(tc.boost_pressure, 0.0);
}

#[test]
fn test_turbo_compressor_pr_ge_one() {
    let tc = TurboCharger::new();
    let pr = tc.compressor_pr(0.10, 80_000.0);
    assert!(pr >= 1.0);
}

#[test]
fn test_turbo_efficiency_in_range() {
    let tc = TurboCharger::new();
    let eff = tc.compressor_efficiency(0.15, 100_000.0);
    assert!((0.0..=1.0).contains(&eff));
}

#[test]
fn test_turbo_step_builds_boost() {
    let mut tc = TurboCharger::new();
    tc.turbo_rpm = 80_000.0;
    for _ in 0..100 {
        tc.step(0.01, 50_000.0);
    }
    assert!(tc.boost_pressure >= 0.0);
}

#[test]
fn test_wastegate_opens_above_setpoint() {
    let mut tc = TurboCharger::new();
    tc.boost_pressure = tc.wastegate_setpoint + 10_000.0;
    tc.step(0.001, 1_000.0);
    assert!(tc.wastegate_dc > 0.0);
}

// ElectricMotor tests

#[test]
fn test_motor_available_torque_in_constant_region() {
    let m = ElectricMotor::new_pmsm(250.0, 150_000.0, 10_000.0);
    let t = m.available_torque(2000.0); // below base RPM
    assert!((t - 250.0).abs() < 1e-9);
}

#[test]
fn test_motor_torque_decreases_above_base_rpm() {
    let m = ElectricMotor::new_pmsm(250.0, 150_000.0, 10_000.0);
    let t_base = m.available_torque(m.base_rpm);
    let t_high = m.available_torque(m.max_rpm);
    assert!(t_high <= t_base);
}

#[test]
fn test_motor_efficiency_in_range() {
    let m = ElectricMotor::new_pmsm(250.0, 150_000.0, 10_000.0);
    let eff = m.efficiency(100.0, 3000.0);
    assert!((0.0..=1.0).contains(&eff));
}

#[test]
fn test_motor_derating_factor_at_normal_temp() {
    let m = ElectricMotor::new_pmsm(250.0, 150_000.0, 10_000.0);
    assert!((m.derating_factor() - 1.0).abs() < 1e-9);
}

// Battery tests

#[test]
fn test_battery_ocv_increases_with_soc() {
    let bat = Battery::new_nmc(60.0, 400.0);
    let v1 = bat.ocv(0.2);
    let v2 = bat.ocv(0.8);
    assert!(v2 > v1);
}

#[test]
fn test_battery_soc_decreases_on_discharge() {
    let mut bat = Battery::new_nmc(60.0, 400.0);
    let soc0 = bat.soc;
    bat.step(100.0, 1.0); // 100A discharge for 1s
    assert!(bat.soc < soc0);
}

#[test]
fn test_battery_terminal_voltage_drops_under_load() {
    let bat = Battery::new_nmc(60.0, 400.0);
    let v0 = bat.terminal_voltage(0.0);
    let v_load = bat.terminal_voltage(200.0);
    assert!(v_load < v0);
}

#[test]
fn test_battery_soc_clamped() {
    let mut bat = Battery::new_nmc(10.0, 400.0);
    bat.soc = 0.01;
    bat.step(10_000.0, 1.0);
    assert!(bat.soc >= 0.0);
}

// HybridPowertrain tests

#[test]
fn test_hybrid_determine_mode_regen_on_brake() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mot = ElectricMotor::new_pmsm(200.0, 100_000.0, 8000.0);
    let bat = Battery::new_nmc(40.0, 350.0);
    let hyb = HybridPowertrain::new(eng, mot, bat);
    let mode = hyb.determine_mode(0.0, true);
    assert_eq!(mode, HybridMode::RegenBraking);
}

#[test]
fn test_hybrid_power_split_electric_only() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mot = ElectricMotor::new_pmsm(200.0, 100_000.0, 8000.0);
    let bat = Battery::new_nmc(40.0, 350.0);
    let mut hyb = HybridPowertrain::new(eng, mot, bat);
    hyb.mode = HybridMode::ElectricOnly;
    let (eng_p, _mot_p) = hyb.power_split(50_000.0);
    assert_eq!(eng_p, 0.0);
}

#[test]
fn test_hybrid_step_increments_fuel() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mot = ElectricMotor::new_pmsm(200.0, 100_000.0, 8000.0);
    let bat = Battery::new_nmc(40.0, 350.0);
    let mut hyb = HybridPowertrain::new(eng, mot, bat);
    hyb.mode = HybridMode::EngineOnly;
    hyb.engine.rpm = 3000.0;
    hyb.ems = EmsStrategy::Cdcs;
    hyb.battery.soc = 0.5; // below the target, so the engine carries the demand
    hyb.step(1.0, 30_000.0, false);
    assert_eq!(hyb.mode, HybridMode::EngineOnly);
    let expected = hyb.engine.fuel_flow_for_power(3000.0, 30_000.0) * 1.0;
    assert!(hyb.total_fuel_kg > 0.0);
    assert!((hyb.total_fuel_kg - expected).abs() < 1e-15);
}

// CVT tests

#[test]
fn test_cvt_output_torque_ratio_multiplied() {
    let cvt = CvtTransmission::new();
    let out = cvt.output_torque(100.0);
    assert!(out > 100.0); // ratio > 1 multiplies torque (minus losses)
}

#[test]
fn test_cvt_step_changes_ratio() {
    let mut cvt = CvtTransmission::new();
    cvt.ratio = 2.0;
    // rate_limit = 0.5/s, so after 1s we move 0.5 toward target 0.5
    // new ratio = 2.0 - 0.5 = 1.5
    cvt.step(1.0, 0.5);
    assert!(cvt.ratio < 2.0, "ratio should decrease toward target");
    assert!(cvt.ratio >= cvt.ratio_min, "ratio stays above minimum");
}

#[test]
fn test_cvt_ratio_stays_within_bounds() {
    let mut cvt = CvtTransmission::new();
    cvt.step(10.0, 0.0);
    assert!(cvt.ratio >= cvt.ratio_min);
    cvt.step(10.0, 100.0);
    assert!(cvt.ratio <= cvt.ratio_max);
}

#[test]
fn test_cvt_efficiency_in_range() {
    let cvt = CvtTransmission::new();
    let eff = cvt.efficiency(100.0, 1.5);
    assert!((0.0..=1.0).contains(&eff));
}

// DCT tests

#[test]
fn test_dct_create_7speed() {
    let dct = DualClutchTransmission::new_7speed();
    assert_eq!(dct.ratios.len(), 7);
    assert_eq!(dct.current_gear, 1);
}

#[test]
fn test_dct_upshift_increments_gear() {
    let mut dct = DualClutchTransmission::new_7speed();
    dct.upshift();
    // Simulate shift
    for _ in 0..20 {
        dct.step(0.01);
    }
    assert_eq!(dct.current_gear, 2);
}

#[test]
fn test_dct_downshift_decrements_gear() {
    let mut dct = DualClutchTransmission::new_7speed();
    // First go to gear 2
    dct.upshift();
    for _ in 0..20 {
        dct.step(0.01);
    }
    // Now downshift
    dct.downshift();
    for _ in 0..20 {
        dct.step(0.01);
    }
    assert_eq!(dct.current_gear, 1);
}

#[test]
fn test_dct_total_ratio_positive() {
    let dct = DualClutchTransmission::new_7speed();
    assert!(dct.total_ratio() > 0.0);
}

#[test]
fn test_dct_output_torque_positive() {
    let dct = DualClutchTransmission::new_7speed();
    let t = dct.output_torque(200.0);
    assert!(t > 0.0);
}

// ThrottleByWire tests

#[test]
fn test_tbw_default_throttle_zero() {
    let tbw = ThrottleByWire::new();
    assert_eq!(tbw.throttle_pos, 0.0);
}

#[test]
fn test_tbw_pedal_increases_throttle() {
    let mut tbw = ThrottleByWire::new();
    tbw.process_pedal(0.5);
    assert!(tbw.target_throttle > 0.0);
}

#[test]
fn test_tbw_step_converges_to_target() {
    let mut tbw = ThrottleByWire::new();
    tbw.process_pedal(1.0);
    for _ in 0..100 {
        tbw.step(0.1);
    }
    assert!((tbw.throttle_pos - tbw.target_throttle).abs() < 0.01);
}

#[test]
fn test_tbw_fault_limits_throttle() {
    let mut tbw = ThrottleByWire::new();
    tbw.trigger_fault();
    tbw.process_pedal(1.0);
    assert!(tbw.target_throttle <= tbw.limp_home_limit);
}

#[test]
fn test_tbw_tcs_limit() {
    let mut tbw = ThrottleByWire::new();
    tbw.set_tcs_limit(0.3);
    tbw.process_pedal(1.0);
    assert!(tbw.target_throttle <= 0.3 + 1e-9);
}

#[test]
fn test_tbw_clear_fault() {
    let mut tbw = ThrottleByWire::new();
    tbw.trigger_fault();
    tbw.clear_fault();
    tbw.process_pedal(0.8);
    for _ in 0..100 {
        tbw.step(0.1);
    }
    assert!(tbw.throttle_pos > tbw.limp_home_limit);
}

#[test]
fn test_intercooler_reduces_temp() {
    let ic = Intercooler::new();
    let outlet = ic.outlet_temp(180.0);
    assert!((25.0..180.0).contains(&outlet));
}

#[test]
fn test_pedal_map_sport_aggressive() {
    let sport = PedalMap::sport();
    let t = sport.evaluate(0.2);
    let linear = PedalMap::linear();
    let t_lin = linear.evaluate(0.2);
    assert!(t > t_lin);
}

#[test]
fn test_wiebe_heat_release_nonnegative() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let hrr = eng.heat_release_rate(10.0, 1000.0);
    assert!(hrr >= 0.0);
}

#[test]
fn test_engine_set_throttle_clamps() {
    let mut eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    eng.set_throttle(2.0);
    assert!(eng.throttle <= 1.0);
    eng.set_throttle(-1.0);
    assert!(eng.throttle >= 0.0);
}

#[test]
fn test_hybrid_fuel_economy_nonnegative() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mot = ElectricMotor::new_pmsm(200.0, 100_000.0, 8000.0);
    let bat = Battery::new_nmc(40.0, 350.0);
    let hyb = HybridPowertrain::new(eng, mot, bat);
    let fe = hyb.fuel_economy(20_000.0, 15.0);
    assert!(fe >= 0.0);
}

#[test]
fn test_dct_recommended_gear_in_range() {
    let dct = DualClutchTransmission::new_7speed();
    let g = dct.recommended_gear(20.0, 0.5);
    assert!(g >= 1 && g <= dct.ratios.len());
}

// ─── Fuel flow for a requested power ────────────────────────────────────────

#[test]
fn test_fuel_flow_for_power_is_power_times_bsfc() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let power_w = 60_000.0;
    let relative_load = power_w / eng.peak_power;
    let expected = power_w / 1000.0 * eng.bsfc(3000.0, relative_load) / 3_600_000.0;
    let got = eng.fuel_flow_for_power(3000.0, power_w);
    assert!((got - expected).abs() < 1e-15, "{got} vs {expected}");
    // 60 kW at ~277 g/kWh is about 4.6 g/s
    assert!((got - 4.6e-3).abs() < 0.1e-3, "fuel flow {got} kg/s");
}

#[test]
fn test_fuel_flow_for_power_zero_negative_and_nan() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    assert_eq!(eng.fuel_flow_for_power(3000.0, 0.0), 0.0);
    assert_eq!(eng.fuel_flow_for_power(3000.0, -5_000.0), 0.0);
    assert_eq!(eng.fuel_flow_for_power(3000.0, f64::NAN), 0.0);
}

#[test]
fn test_fuel_flow_for_power_increases_with_power() {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mut previous = 0.0;
    for step in 1..=14 {
        let flow = eng.fuel_flow_for_power(3000.0, step as f64 * 10_000.0);
        assert!(flow > previous, "flow not increasing at step {step}");
        previous = flow;
    }
}

// ─── Hybrid powertrain: battery current sign and generator credit ───────────

fn hybrid_fixture() -> HybridPowertrain {
    let eng = CombustionEngine::new(2.0e-3, 4, 10.5);
    let mot = ElectricMotor::new_pmsm(200.0, 100_000.0, 8000.0);
    let bat = Battery::new_nmc(40.0, 350.0);
    HybridPowertrain::new(eng, mot, bat)
}

/// SOC lost by the 40 Ah / 350 V fixture battery when `power_w` flows out of it
/// for `dt` seconds.
fn soc_drop(power_w: f64, dt: f64) -> f64 {
    power_w / 350.0 * dt / 3600.0 / 40.0
}

#[test]
fn test_hybrid_electric_only_step_lowers_soc() {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Cdcs; // SOC 0.8 above the 0.6 target -> electric only
    let soc0 = hyb.battery.soc;
    hyb.step(10.0, 20_000.0, false);
    assert_eq!(hyb.mode, HybridMode::ElectricOnly);
    let drop = soc0 - hyb.battery.soc;
    assert!(drop > 0.0, "motoring must discharge the battery");
    assert!(
        (drop - soc_drop(20_000.0, 10.0)).abs() < 1e-12,
        "drop {drop}"
    );
}

#[test]
fn test_hybrid_regen_braking_step_raises_soc() {
    let mut hyb = hybrid_fixture();
    let soc0 = hyb.battery.soc;
    hyb.step(2.0, 30_000.0, true);
    assert_eq!(hyb.mode, HybridMode::RegenBraking);
    let gain = hyb.battery.soc - soc0;
    assert!(gain > 0.0, "regenerative braking must charge the battery");
    assert!(
        (gain - soc_drop(30_000.0, 2.0)).abs() < 1e-12,
        "gain {gain}"
    );
    assert!((hyb.regen_energy_j - 60_000.0).abs() < 1e-9);
}

#[test]
fn test_hybrid_rule_based_parallel_step_discharges_battery() {
    let mut hyb = hybrid_fixture();
    let soc0 = hyb.battery.soc;
    hyb.step(5.0, 50_000.0, false);
    assert_eq!(hyb.mode, HybridMode::Parallel);
    // the motor carries 40 % of the demand
    let drop = soc0 - hyb.battery.soc;
    assert!(
        (drop - soc_drop(20_000.0, 5.0)).abs() < 1e-12,
        "drop {drop}"
    );
    assert_eq!(hyb.regen_energy_j, 0.0);
}

#[test]
fn test_hybrid_engine_only_step_leaves_battery_alone() {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Cdcs;
    hyb.battery.soc = 0.5;
    hyb.step(10.0, 40_000.0, false);
    assert_eq!(hyb.mode, HybridMode::EngineOnly);
    assert_eq!(hyb.battery.soc, 0.5);
}

#[test]
fn test_hybrid_charging_mode_credits_engine_power_to_battery() {
    let mut hyb = hybrid_fixture();
    hyb.mode = HybridMode::Charging;
    let soc0 = hyb.battery.soc;
    hyb.advance(10.0, 20_000.0);
    let gain = hyb.battery.soc - soc0;
    assert!(gain > 0.0, "charging must raise the SOC");
    assert!(
        (gain - soc_drop(20_000.0, 10.0)).abs() < 1e-12,
        "gain {gain}"
    );
    assert!(hyb.total_fuel_kg > 0.0, "the generator engine burns fuel");
    assert_eq!(
        hyb.regen_energy_j, 0.0,
        "engine charging is not regeneration"
    );
}

#[test]
fn test_hybrid_series_mode_battery_sees_only_the_difference() {
    // demand below the engine rating: the generator covers the motor entirely
    let mut hyb = hybrid_fixture();
    hyb.mode = HybridMode::Series;
    let soc0 = hyb.battery.soc;
    hyb.advance(10.0, 100_000.0);
    assert!((hyb.battery.soc - soc0).abs() < 1e-15);

    // demand above the engine rating: the battery makes up the rest
    let mut hyb = hybrid_fixture();
    hyb.mode = HybridMode::Series;
    let soc0 = hyb.battery.soc;
    hyb.advance(10.0, 200_000.0);
    let drop = soc0 - hyb.battery.soc;
    let missing_w = 200_000.0 - hyb.engine.peak_power;
    assert!(
        (drop - soc_drop(missing_w, 10.0)).abs() < 1e-12,
        "drop {drop}"
    );
}

#[test]
fn test_hybrid_battery_power_per_mode() {
    let mut hyb = hybrid_fixture();
    for mode in [
        HybridMode::ElectricOnly,
        HybridMode::EngineOnly,
        HybridMode::Parallel,
        HybridMode::RegenBraking,
    ] {
        hyb.mode = mode;
        assert_eq!(hyb.battery_power(30_000.0, 12_000.0), 12_000.0, "{mode:?}");
        assert_eq!(hyb.battery_power(0.0, -7_000.0), -7_000.0, "{mode:?}");
    }
    hyb.mode = HybridMode::Series;
    assert_eq!(hyb.battery_power(30_000.0, 50_000.0), 20_000.0);
    hyb.mode = HybridMode::Charging;
    assert_eq!(hyb.battery_power(30_000.0, 0.0), -30_000.0);
}

#[test]
fn test_hybrid_power_split_conserves_energy_in_every_mode() {
    // Lossless bookkeeping: the engine power plus the power drawn from the battery
    // is the mechanical power at the wheels.  The wheels receive the demand while
    // driving, the (negative) recuperated power while braking and nothing at all
    // while the car stands and charges.  The Series generator feeds the motor
    // directly, so it is counted once, through the battery difference.
    for ems in [EmsStrategy::RuleBased, EmsStrategy::Cdcs, EmsStrategy::Ecms] {
        for soc in [0.3, 0.6, 0.85] {
            for demand_w in [5_000.0, 30_000.0, 80_000.0, 95_000.0] {
                for (mode, wheel_power_w) in [
                    (HybridMode::ElectricOnly, demand_w),
                    (HybridMode::EngineOnly, demand_w),
                    (HybridMode::Parallel, demand_w),
                    (HybridMode::Series, demand_w),
                    (HybridMode::RegenBraking, -demand_w),
                    (HybridMode::Charging, 0.0),
                ] {
                    let mut hyb = hybrid_fixture();
                    hyb.ems = ems;
                    hyb.engine.rpm = 3000.0;
                    hyb.battery.soc = soc;
                    hyb.mode = mode;
                    let (engine_w, motor_w) = hyb.power_split(demand_w);
                    let delivered_w = engine_w + hyb.battery_power(engine_w, motor_w);
                    assert!(
                        (delivered_w - wheel_power_w).abs() < 1e-6,
                        "{ems:?} {mode:?} soc {soc} demand {demand_w}: engine {engine_w} W, \
                         motor {motor_w} W give {delivered_w} W instead of {wheel_power_w} W"
                    );
                }
            }
        }
    }
}

#[test]
fn test_hybrid_regen_counter_ignores_engine_charging_while_driving() {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Ecms;
    hyb.engine.rpm = 3000.0;
    hyb.battery.soc = 0.3;
    hyb.step(1.0, 30_000.0, false);
    assert!(hyb.battery.soc > 0.3, "the engine recharges the pack");
    assert_eq!(hyb.regen_energy_j, 0.0);
}

// ─── Hybrid powertrain: motor speed follows the engine ──────────────────────

#[test]
fn test_hybrid_step_sets_motor_speed_from_engine_speed() {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Cdcs;
    hyb.engine.rpm = 3000.0;
    hyb.motor_ratio = 2.0;
    hyb.step(1.0, 50_000.0, false);
    assert_eq!(hyb.mode, HybridMode::ElectricOnly);
    assert!(
        (hyb.motor.rpm - 6000.0).abs() < 1e-9,
        "rpm {}",
        hyb.motor.rpm
    );
    let omega = 6000.0 * std::f64::consts::PI / 30.0;
    assert!((hyb.motor.torque_cmd - 50_000.0 / omega).abs() < 1e-9);
}

#[test]
fn test_hybrid_motor_heats_while_it_runs() {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Cdcs;
    hyb.engine.rpm = 3000.0;
    assert_eq!(hyb.motor.winding_temp, 25.0);
    hyb.step(1.0, 50_000.0, false);
    assert!(hyb.motor.winding_temp > 25.0, "{}", hyb.motor.winding_temp);
}

#[test]
fn test_hybrid_default_motor_ratio_is_one() {
    assert_eq!(hybrid_fixture().motor_ratio, 1.0);
}

// ─── ECMS ───────────────────────────────────────────────────────────────────

fn ecms_fixture(soc: f64) -> HybridPowertrain {
    let mut hyb = hybrid_fixture();
    hyb.ems = EmsStrategy::Ecms;
    hyb.engine.rpm = 3000.0;
    hyb.battery.soc = soc;
    hyb
}

#[test]
fn test_ecms_low_soc_selects_less_motor_power_than_high_soc() {
    for demand_w in [30_000.0, 60_000.0, 100_000.0, 150_000.0] {
        let (_, motor_high) = ecms_fixture(0.85).ecms_split(demand_w);
        let (_, motor_low) = ecms_fixture(0.30).ecms_split(demand_w);
        assert!(
            motor_low < motor_high,
            "demand {demand_w}: motor {motor_low} W at SOC 0.3 vs {motor_high} W at SOC 0.85"
        );
    }
}

#[test]
fn test_ecms_motor_power_is_monotone_in_soc() {
    for demand_w in [40_000.0, 100_000.0, 150_000.0] {
        let mut previous = f64::NEG_INFINITY;
        for step in 0..=12 {
            let soc = 0.25 + 0.05 * step as f64;
            let (_, motor) = ecms_fixture(soc).ecms_split(demand_w);
            assert!(
                motor >= previous,
                "demand {demand_w}, SOC {soc}: motor {motor} W after {previous} W"
            );
            previous = motor;
        }
    }
}

#[test]
fn test_ecms_reads_the_equivalence_factor() {
    let mut cheap = ecms_fixture(0.6);
    cheap.lambda_ecms = 1.0;
    let mut dear = ecms_fixture(0.6);
    dear.lambda_ecms = 6.0;
    let (engine_cheap, motor_cheap) = cheap.ecms_split(60_000.0);
    let (engine_dear, motor_dear) = dear.ecms_split(60_000.0);
    assert!(
        motor_cheap > 0.0 && engine_cheap == 0.0,
        "cheap electricity: {engine_cheap}/{motor_cheap}"
    );
    assert!(
        motor_dear < 0.0 && engine_dear > 60_000.0,
        "dear electricity: {engine_dear}/{motor_dear}"
    );
}

#[test]
fn test_ecms_reads_the_soc_target() {
    // Moving the target up makes the same SOC look low: electricity gets dearer.
    let mut relaxed = ecms_fixture(0.6);
    relaxed.soc_target = 0.3;
    let mut strict = ecms_fixture(0.6);
    strict.soc_target = 0.9;
    let (_, motor_relaxed) = relaxed.ecms_split(100_000.0);
    let (_, motor_strict) = strict.ecms_split(100_000.0);
    assert!(motor_strict < motor_relaxed);
}

#[test]
fn test_ecms_split_meets_demand_and_respects_limits() {
    let motor_max = 100_000.0;
    let engine_max = 147_000.0;
    for soc in [0.1, 0.2, 0.25, 0.5, 0.6, 0.85, 0.9, 0.95] {
        let hyb = ecms_fixture(soc);
        let discharge_allowed = soc > 0.2;
        let charge_allowed = soc < 0.9;
        let deliverable = engine_max + if discharge_allowed { motor_max } else { 0.0 };
        for demand_w in [
            0.0, 5_000.0, 40_000.0, 90_000.0, 140_000.0, 200_000.0, 240_000.0, 400_000.0,
        ] {
            let (engine, motor) = hyb.ecms_split(demand_w);
            let context = format!("SOC {soc}, demand {demand_w}: engine {engine}, motor {motor}");
            assert!((0.0..=engine_max).contains(&engine), "{context}");
            assert!(motor.abs() <= motor_max, "{context}");
            assert!(discharge_allowed || motor <= 0.0, "{context}");
            assert!(charge_allowed || motor >= 0.0, "{context}");
            if demand_w <= deliverable {
                assert!((engine + motor - demand_w).abs() < 1e-6, "{context}");
            } else {
                // beyond what the powertrain can deliver both run flat out
                assert_eq!(engine, engine_max, "{context}");
                assert_eq!(
                    motor,
                    if discharge_allowed { motor_max } else { 0.0 },
                    "{context}"
                );
            }
        }
    }
}

#[test]
fn test_ecms_runs_both_flat_out_beyond_the_combined_rating() {
    let (engine, motor) = ecms_fixture(0.6).ecms_split(400_000.0);
    assert_eq!(engine, 147_000.0);
    assert_eq!(motor, 100_000.0);
}

#[test]
fn test_ecms_cannot_discharge_an_empty_battery() {
    let (engine, motor) = ecms_fixture(0.15).ecms_split(150_000.0);
    assert_eq!(engine, 147_000.0);
    assert_eq!(motor, 0.0);
}

#[test]
fn test_ecms_overheated_motor_contributes_nothing() {
    let mut hyb = ecms_fixture(0.85);
    hyb.motor.winding_temp = 150.0;
    assert_eq!(hyb.motor.derating_factor(), 0.0);
    assert_eq!(hyb.ecms_split(60_000.0), (60_000.0, 0.0));
}

#[test]
fn test_ecms_choice_is_the_minimum_of_its_cost() {
    for soc in [0.3, 0.45, 0.6, 0.75, 0.85] {
        let hyb = ecms_fixture(soc);
        let s = hyb.lambda_ecms * (1.0 + (hyb.soc_target - soc));
        let cost = |demand_w: f64, motor_w: f64| -> f64 {
            let engine_w = (demand_w - motor_w).clamp(0.0, hyb.engine.peak_power);
            hyb.engine.fuel_flow_for_power(hyb.engine.rpm, engine_w) + s * motor_w / 43.0e6
        };
        for demand_w in [30_000.0, 80_000.0, 120_000.0, 150_000.0] {
            let (_, chosen) = hyb.ecms_split(demand_w);
            let lo = (-100_000.0_f64).max(demand_w - hyb.engine.peak_power);
            let hi = 100_000.0_f64.min(demand_w);
            let fine_best = (0..=4000)
                .map(|i| cost(demand_w, lo + (hi - lo) * i as f64 / 4000.0))
                .fold(f64::INFINITY, f64::min);
            let gap = cost(demand_w, chosen) - fine_best;
            assert!(
                (-1e-12..5e-6).contains(&gap),
                "SOC {soc}, demand {demand_w}: chosen cost exceeds the optimum by {gap} kg/s"
            );
        }
    }
}

#[test]
fn test_ecms_negative_demand_recuperates_within_limits() {
    assert_eq!(ecms_fixture(0.5).ecms_split(-30_000.0), (0.0, -30_000.0));
    assert_eq!(ecms_fixture(0.5).ecms_split(-250_000.0), (0.0, -100_000.0));
    // a full battery cannot take the energy back
    assert_eq!(ecms_fixture(0.95).ecms_split(-30_000.0), (0.0, 0.0));
    assert_eq!(ecms_fixture(0.5).ecms_split(f64::NAN), (0.0, 0.0));
}

#[test]
fn test_ecms_determine_mode_follows_the_split() {
    assert_eq!(
        ecms_fixture(0.85).determine_mode(30_000.0, false),
        HybridMode::ElectricOnly
    );
    assert_eq!(
        ecms_fixture(0.6).determine_mode(150_000.0, false),
        HybridMode::Parallel
    );
    assert_eq!(
        ecms_fixture(0.15).determine_mode(150_000.0, false),
        HybridMode::EngineOnly
    );
    assert_eq!(
        ecms_fixture(0.6).determine_mode(-20_000.0, false),
        HybridMode::RegenBraking
    );
    assert_eq!(
        ecms_fixture(0.6).determine_mode(20_000.0, true),
        HybridMode::RegenBraking
    );
}

#[test]
fn test_ecms_power_split_follows_ecms_split_in_propulsion_modes() {
    let mut hyb = ecms_fixture(0.6);
    let expected = hyb.ecms_split(150_000.0);
    for mode in [
        HybridMode::ElectricOnly,
        HybridMode::EngineOnly,
        HybridMode::Parallel,
    ] {
        hyb.mode = mode;
        assert_eq!(hyb.power_split(150_000.0), expected, "{mode:?}");
    }
    hyb.mode = HybridMode::RegenBraking;
    assert_eq!(hyb.power_split(40_000.0), (0.0, -40_000.0));
}

#[test]
fn test_ecms_step_at_high_soc_drives_electrically() {
    let mut hyb = ecms_fixture(0.85);
    hyb.step(5.0, 30_000.0, false);
    assert_eq!(hyb.mode, HybridMode::ElectricOnly);
    assert_eq!(hyb.total_fuel_kg, 0.0);
    assert!((0.85 - hyb.battery.soc - soc_drop(30_000.0, 5.0)).abs() < 1e-12);
    assert_eq!(hyb.engine.throttle, 0.0);
}

#[test]
fn test_ecms_step_at_low_soc_runs_the_engine_and_recharges() {
    let mut hyb = ecms_fixture(0.30);
    hyb.step(5.0, 30_000.0, false);
    assert_eq!(hyb.mode, HybridMode::Parallel);
    assert!(hyb.total_fuel_kg > 0.0);
    assert!(hyb.battery.soc > 0.30, "SOC {}", hyb.battery.soc);
    assert!(
        hyb.engine.throttle > 0.5,
        "throttle {}",
        hyb.engine.throttle
    );
}
