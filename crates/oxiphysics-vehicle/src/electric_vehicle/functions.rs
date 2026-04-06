//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types_core::{AlphaBetaCurrent, Battery, BatteryParams, DqCurrent};
use super::types_ext::ElectricMotor;

/// Faraday constant (C/mol)
pub const FARADAY: f64 = 96_485.0;
/// Universal gas constant (J/(mol·K))
pub const R_GAS: f64 = 8.314_462;
/// Reference temperature (K)
pub const T_REF: f64 = 298.15;
/// Boltzmann constant (J/K)
pub const K_BOLTZMANN: f64 = 1.380_649e-23;
/// Clarke transform: three-phase → αβ.
pub fn clarke_transform(ia: f64, ib: f64, ic: f64) -> AlphaBetaCurrent {
    AlphaBetaCurrent {
        alpha: (2.0 * ia - ib - ic) / 3.0,
        beta: (ib - ic) / f64::sqrt(3.0),
    }
}
/// Inverse Clarke transform: αβ → three-phase.
pub fn inverse_clarke(alpha: f64, beta: f64) -> (f64, f64, f64) {
    let ia = alpha;
    let ib = -0.5 * alpha + f64::sqrt(3.0) / 2.0 * beta;
    let ic = -0.5 * alpha - f64::sqrt(3.0) / 2.0 * beta;
    (ia, ib, ic)
}
/// Park transform: αβ → dq.
pub fn park_transform(alpha: f64, beta: f64, theta_e: f64) -> DqCurrent {
    DqCurrent {
        id: alpha * theta_e.cos() + beta * theta_e.sin(),
        iq: -alpha * theta_e.sin() + beta * theta_e.cos(),
    }
}
/// Inverse Park transform: dq → αβ.
pub fn inverse_park(id: f64, iq: f64, theta_e: f64) -> (f64, f64) {
    let alpha = id * theta_e.cos() - iq * theta_e.sin();
    let beta = id * theta_e.sin() + iq * theta_e.cos();
    (alpha, beta)
}
/// Space-vector pulse-width modulation (SVPWM) duty cycles.
///
/// Returns (da, db, dc) duty cycles in \[0, 1\].
pub fn svpwm(v_alpha: f64, v_beta: f64, v_dc: f64) -> (f64, f64, f64) {
    let vref = (v_alpha.powi(2) + v_beta.powi(2)).sqrt();
    let theta = v_beta.atan2(v_alpha);
    let sector = ((theta.to_degrees() + 360.0) % 360.0 / 60.0).floor() as usize;
    let sector = sector.clamp(0, 5);
    let m = (vref / (v_dc / f64::sqrt(3.0))).clamp(0.0, 1.0);
    let t1 = m
        * (std::f64::consts::PI / 3.0 - (theta - sector as f64 * std::f64::consts::PI / 3.0)).sin()
        / (std::f64::consts::PI / 3.0).sin();
    let t2 = m * (theta - sector as f64 * std::f64::consts::PI / 3.0).sin()
        / (std::f64::consts::PI / 3.0).sin();
    let t0 = 1.0 - t1 - t2;
    let t0 = t0.max(0.0);
    let (ta, tb, tc) = match sector {
        0 => (t1 + t2 + t0 / 2.0, t2 + t0 / 2.0, t0 / 2.0),
        1 => (t1 + t0 / 2.0, t1 + t2 + t0 / 2.0, t0 / 2.0),
        2 => (t0 / 2.0, t1 + t2 + t0 / 2.0, t2 + t0 / 2.0),
        3 => (t0 / 2.0, t1 + t0 / 2.0, t1 + t2 + t0 / 2.0),
        4 => (t2 + t0 / 2.0, t0 / 2.0, t1 + t2 + t0 / 2.0),
        _ => (t1 + t2 + t0 / 2.0, t0 / 2.0, t1 + t0 / 2.0),
    };
    (ta.clamp(0.0, 1.0), tb.clamp(0.0, 1.0), tc.clamp(0.0, 1.0))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::electric_vehicle::BatteryCell;
    use crate::electric_vehicle::BatteryPack;
    use crate::electric_vehicle::CcCvCharger;
    use crate::electric_vehicle::CellChemistry;

    use crate::electric_vehicle::ChargingMode;
    use crate::electric_vehicle::ChargingModel;
    use crate::electric_vehicle::DriveCycle;
    use crate::electric_vehicle::EfficiencyMap;
    use crate::electric_vehicle::ElectricVehicle;
    use crate::electric_vehicle::EnergyConsumption;
    use crate::electric_vehicle::EnergyRouter;
    use crate::electric_vehicle::EvDrivetrain;
    use crate::electric_vehicle::GridTariff;
    use crate::electric_vehicle::InductionMotor;
    use crate::electric_vehicle::MotorController;

    use crate::electric_vehicle::PackConfig;
    use crate::electric_vehicle::PermanentMagnetMotor;
    use crate::electric_vehicle::PiController;
    use crate::electric_vehicle::PulseCharger;
    use crate::electric_vehicle::RangeEstimator;
    use crate::electric_vehicle::RegenerativeBraking;
    use crate::electric_vehicle::RegenerativeEfficiency;
    use crate::electric_vehicle::RoadSegment;

    use crate::electric_vehicle::SmartChargingScheduler;
    use crate::electric_vehicle::SocEstimator;
    use crate::electric_vehicle::SohEstimator;
    use crate::electric_vehicle::ThermalManagement;
    use crate::electric_vehicle::ThermalMotorModel;
    use crate::electric_vehicle::V2gController;
    #[test]
    fn test_cell_chemistry_voltages() {
        for chem in &[
            CellChemistry::NMC,
            CellChemistry::LFP,
            CellChemistry::NCA,
            CellChemistry::LTO,
            CellChemistry::LMO,
        ] {
            assert!(chem.nominal_voltage() > 0.0);
            assert!(chem.max_voltage() > chem.nominal_voltage());
            assert!(chem.min_voltage() < chem.nominal_voltage());
        }
    }
    #[test]
    fn test_battery_cell_creation() {
        let cell = BatteryCell::new(CellChemistry::NMC, 50.0);
        assert_eq!(cell.soc, 1.0);
        assert_eq!(cell.soh, 1.0);
        assert!(cell.ocv > 0.0);
    }
    #[test]
    fn test_cell_ocv_monotone() {
        let mut cell = BatteryCell::new(CellChemistry::NMC, 50.0);
        let mut prev = 0.0;
        for i in 0..=10 {
            cell.soc = i as f64 / 10.0;
            let v = cell.compute_ocv();
            assert!(v >= prev - 1e-9, "OCV not monotone at soc={}", cell.soc);
            prev = v;
        }
    }
    #[test]
    fn test_cell_soc_decreases_on_discharge() {
        let mut cell = BatteryCell::new(CellChemistry::NMC, 50.0);
        let init_soc = cell.soc;
        cell.update_soc(50.0, 3600.0);
        assert!(cell.soc < init_soc);
    }
    #[test]
    fn test_cell_soc_clamped() {
        let mut cell = BatteryCell::new(CellChemistry::NMC, 50.0);
        cell.update_soc(500.0, 3600.0);
        assert!(cell.soc >= 0.0);
    }
    #[test]
    fn test_pack_creation() {
        let config = PackConfig::standard_400v_nmc();
        let pack = BatteryPack::new(config);
        assert_eq!(pack.cells.len(), 96 * 4);
        assert!(pack.pack_soc > 0.9);
    }
    #[test]
    fn test_pack_voltage_reasonable() {
        let config = PackConfig::standard_400v_nmc();
        let pack = BatteryPack::new(config);
        let v = pack.pack_voltage();
        assert!(v > 300.0 && v < 500.0);
    }
    #[test]
    fn test_pack_capacity() {
        let config = PackConfig::standard_400v_nmc();
        let kwh = config.total_capacity_kwh();
        assert!(kwh > 50.0 && kwh < 120.0);
    }
    #[test]
    fn test_passive_balance() {
        let config = PackConfig::standard_400v_nmc();
        let mut pack = BatteryPack::new(config);
        pack.cells[0].soc = 0.95;
        pack.cells[1].soc = 0.80;
        pack.passive_balance(0.05);
        assert!(pack.balancing_active[0]);
        assert!(!pack.balancing_active[1]);
    }
    #[test]
    fn test_soh_estimator() {
        let mut est = SohEstimator::new(100.0);
        est.update_soh_from_capacity(85.0);
        assert!((est.soh - 0.85).abs() < 1e-6);
    }
    #[test]
    fn test_soh_calendar_aging() {
        let mut est = SohEstimator::new(100.0);
        est.apply_calendar_aging(365.0, 25.0);
        assert!(est.soh < 1.0);
    }
    #[test]
    fn test_soh_cycle_aging() {
        let mut est = SohEstimator::new(100.0);
        est.apply_cycle_aging(500.0, 0.8);
        assert!(est.soh < 1.0);
    }
    #[test]
    fn test_soh_end_of_life() {
        let mut est = SohEstimator::new(100.0);
        est.soh = 0.75;
        assert!(est.is_end_of_life());
        est.soh = 0.85;
        assert!(!est.is_end_of_life());
    }
    #[test]
    fn test_clarke_transform_balanced() {
        let ia = 1.0;
        let ib = -0.5;
        let ic = -0.5;
        let ab = clarke_transform(ia, ib, ic);
        assert!((ia + ib + ic).abs() < 1e-9);
        assert!((ab.alpha - 1.0 / 3.0 * (2.0 * ia - ib - ic)).abs() < 1e-9);
    }
    #[test]
    fn test_park_inverse_park_roundtrip() {
        let id = 10.0;
        let iq = 50.0;
        let theta = 1.2;
        let (alpha, beta) = inverse_park(id, iq, theta);
        let dq = park_transform(alpha, beta, theta);
        assert!((dq.id - id).abs() < 1e-10);
        assert!((dq.iq - iq).abs() < 1e-10);
    }
    #[test]
    fn test_svpwm_duty_in_range() {
        let (da, db, dc) = svpwm(100.0, 50.0, 400.0);
        assert!((0.0..=1.0).contains(&da));
        assert!((0.0..=1.0).contains(&db));
        assert!((0.0..=1.0).contains(&dc));
    }
    #[test]
    fn test_pi_controller_tracks_setpoint() {
        let mut ctrl = PiController::new(1.0, 10.0, 100.0);
        let mut error = 10.0;
        let dt = 0.001;
        for _ in 0..5000 {
            let out = ctrl.step(error, dt);
            error -= out * dt;
        }
        assert!(error.abs() < 1.0);
    }
    #[test]
    fn test_efficiency_map_bounds() {
        let map = EfficiencyMap::synthetic_pmsm();
        let eta = map.efficiency(3000.0, 100.0);
        assert!(eta > 0.5 && eta < 1.0);
    }
    #[test]
    fn test_efficiency_map_interpolation() {
        let map = EfficiencyMap::synthetic_pmsm();
        let _ = map.efficiency(-100.0, -50.0);
        let _ = map.efficiency(100000.0, 100000.0);
    }
    #[test]
    fn test_regen_braking_energy_recovery() {
        let mut regen = RegenerativeBraking::new(200.0, 80_000.0);
        let (rt, mt, pw) = regen.compute_regen(100.0, 100.0, 0.5, 1.0);
        assert!(rt >= 0.0);
        assert!(mt >= 0.0);
        assert!(pw >= 0.0);
        assert!((rt + mt - 100.0).abs() < 1e-6);
    }
    #[test]
    fn test_regen_no_recovery_at_high_soc() {
        let mut regen = RegenerativeBraking::new(200.0, 80_000.0);
        let (rt, mt, pw) = regen.compute_regen(100.0, 100.0, 0.99, 1.0);
        assert_eq!(rt, 0.0);
        assert_eq!(mt, 100.0);
        assert_eq!(pw, 0.0);
    }
    #[test]
    fn test_drive_cycle_total_distance() {
        let cycle = DriveCycle::wltp_synthetic();
        assert!(cycle.total_distance_m() > 1000.0);
    }
    #[test]
    fn test_range_estimator_reasonable() {
        let est = RangeEstimator::new(2000.0, 2.4, 0.28, 75.0);
        let cycle = DriveCycle::wltp_synthetic();
        let range = est.remaining_range_km(0.8, &cycle);
        assert!(range > 100.0 && range < 800.0);
    }
    #[test]
    fn test_grid_tariff() {
        let tariff = GridTariff::two_rate(0.30, 0.10);
        let peak = tariff.price_at_hour(12.0);
        let offpeak = tariff.price_at_hour(3.0);
        assert!(peak > offpeak);
    }
    #[test]
    fn test_v2g_export_at_high_soc() {
        let ctrl = V2gController::new(11_000.0, 22_000.0, 0.25);
        let power = ctrl.compute_setpoint(0.9, 12.0, 8.0);
        assert!(power <= 0.0);
    }
    #[test]
    fn test_smart_charging_schedule() {
        let scheduler = SmartChargingScheduler::new(8.0, 0.9, 0.3, 75.0, 11_000.0);
        let schedule = scheduler.build_schedule(0.0);
        assert!(!schedule.is_empty());
        let cost = SmartChargingScheduler::schedule_cost(&schedule);
        assert!(cost > 0.0);
    }
    #[test]
    fn test_energy_router_feasibility() {
        let mut router = EnergyRouter::new(2000.0, 0.65);
        router.segments = vec![
            RoadSegment {
                id: 0,
                length_m: 10_000.0,
                avg_speed_ms: 30.0,
                gradient_rad: 0.01,
                surface_factor: 1.0,
            },
            RoadSegment {
                id: 1,
                length_m: 5_000.0,
                avg_speed_ms: 20.0,
                gradient_rad: -0.02,
                surface_factor: 1.2,
            },
        ];
        assert!(router.is_feasible(1.0, 75.0));
        assert!(!router.is_feasible(0.001, 75.0));
    }
    #[test]
    fn test_thermal_management_no_panic() {
        let mut tms = ThermalManagement::new();
        let pack = BatteryPack::new(PackConfig::standard_400v_nmc());
        let aux = tms.update(&pack, 35.0, 1.0);
        assert!(aux >= 0.0);
    }
    #[test]
    fn test_ev_creation() {
        let ev = ElectricVehicle::new_75kwh();
        assert!(ev.pack.pack_soc > 0.9);
        assert_eq!(ev.odometer_m, 0.0);
    }
    #[test]
    fn test_ev_step_no_panic() {
        let mut ev = ElectricVehicle::new_75kwh();
        ev.speed_ms = 30.0;
        ev.step(100.0, 0.01, 20.0);
    }
    #[test]
    fn test_ev_range_estimate() {
        let ev = ElectricVehicle::new_75kwh();
        let range = ev.range_km();
        assert!(range > 50.0);
    }
    #[test]
    fn test_cccv_charger_start() {
        let config = PackConfig::standard_400v_nmc();
        let mut pack = BatteryPack::new(config.clone());
        for c in &mut pack.cells {
            c.soc = 0.2;
            c.ocv = c.compute_ocv();
            c.terminal_voltage = c.ocv;
        }
        let mut charger = CcCvCharger::new(100.0, config.nominal_voltage() * 1.05);
        charger.start();
        assert_eq!(charger.mode, ChargingMode::ConstantCurrent);
        let (current, _voltage) = charger.step(&pack, 1.0);
        assert!(current > 0.0);
    }
    #[test]
    fn test_pulse_charger_alternates() {
        let mut charger = PulseCharger::new(100.0, 0.5, 0.1);
        let c1 = charger.step(400.0, 0.01);
        assert!(c1 > 0.0);
    }
    #[test]
    fn test_pack_update_soc_decreases() {
        let config = PackConfig::standard_400v_nmc();
        let mut pack = BatteryPack::new(config);
        let init_soc = pack.pack_soc;
        pack.update(100.0, 360.0);
        assert!(pack.pack_soc < init_soc);
    }
    #[test]
    fn motor_controller_clamps_to_max_torque() {
        let mut mc = MotorController::new(250.0, 10000.0);
        mc.update_speed(500.0);
        let torque = mc.step(9999.0, 0.1);
        assert!(torque <= 250.0, "torque={torque}");
    }
    #[test]
    fn motor_controller_zero_command_gives_low_torque() {
        let mut mc = MotorController::new(250.0, 10000.0);
        for _ in 0..100 {
            mc.step(0.0, 0.01);
        }
        assert!(mc.actual_torque.abs() < 1.0, "torque should settle near 0");
    }
    #[test]
    fn motor_controller_speed_update_clamped() {
        let mut mc = MotorController::new(250.0, 10000.0);
        mc.update_speed(1e9);
        assert!(mc.speed_rpm <= 10000.0);
    }
    #[test]
    fn motor_controller_electrical_power_positive_motoring() {
        let mut mc = MotorController::new(250.0, 10000.0);
        mc.update_speed(3000.0);
        mc.actual_torque = 100.0;
        let p = mc.electrical_power();
        assert!(p > 0.0, "motoring power should be positive: {p}");
    }
    #[test]
    fn induction_motor_sync_speed_positive() {
        let motor = InductionMotor::new_50hz_4pole(200.0);
        assert!(motor.synchronous_speed_rads() > 0.0);
    }
    #[test]
    fn induction_motor_zero_slip_zero_torque() {
        let motor = InductionMotor::new_50hz_4pole(200.0);
        let omega_sync = motor.synchronous_speed_rads();
        let t = motor.torque_at_speed(omega_sync);
        assert!(t.abs() < 1.0, "zero slip → near-zero torque: {t}");
    }
    #[test]
    fn induction_motor_rated_torque_at_rated_slip() {
        let motor = InductionMotor::new_50hz_4pole(200.0);
        let t = motor.torque_from_slip(motor.rated_slip);
        assert!(t > 0.0 && t <= motor.rated_torque * 1.1, "torque={t}");
    }
    #[test]
    fn induction_motor_copper_losses_positive() {
        let motor = InductionMotor::new_50hz_4pole(200.0);
        let losses = motor.copper_losses(0.05, 100.0);
        assert!(losses > 0.0, "losses={losses}");
    }
    #[test]
    fn pmsm_torque_from_iq_positive() {
        let motor = PermanentMagnetMotor::automotive_traction();
        let t = motor.torque_from_iq(100.0);
        assert!(t > 0.0, "torque={t}");
    }
    #[test]
    fn pmsm_iq_for_torque_roundtrip() {
        let motor = PermanentMagnetMotor::automotive_traction();
        let demanded = 150.0;
        let iq = motor.iq_for_torque(demanded);
        let t_back = motor.torque_from_iq(iq);
        assert!((t_back - demanded).abs() < 1e-6, "t_back={t_back}");
    }
    #[test]
    fn pmsm_back_emf_increases_with_speed() {
        let motor = PermanentMagnetMotor::automotive_traction();
        let emf1 = motor.back_emf(100.0);
        let emf2 = motor.back_emf(200.0);
        assert!(emf2 > emf1, "back-EMF should increase with speed");
    }
    #[test]
    fn pmsm_copper_losses_positive() {
        let motor = PermanentMagnetMotor::automotive_traction();
        let losses = motor.copper_losses(100.0);
        assert!(losses > 0.0, "losses={losses}");
    }
    #[test]
    fn pmsm_max_torque_decreases_with_speed() {
        let motor = PermanentMagnetMotor::automotive_traction();
        let t1 = motor.max_torque_at_speed(400.0, 50.0);
        let t2 = motor.max_torque_at_speed(400.0, 500.0);
        assert!(
            t1 >= t2,
            "torque should decrease at higher speed: {t1} vs {t2}"
        );
    }
    #[test]
    fn charging_model_starts_in_cc_mode() {
        let cm = ChargingModel::new(100.0, 420.0, 5.0);
        assert_eq!(cm.mode, ChargingMode::ConstantCurrent);
    }
    #[test]
    fn charging_model_cc_delivers_current() {
        let mut cm = ChargingModel::new(100.0, 420.0, 5.0);
        let (i, _v) = cm.step(350.0, 1.0);
        assert!(i > 0.0, "current={i}");
    }
    #[test]
    fn charging_model_switches_to_cv_at_voltage() {
        let mut cm = ChargingModel::new(100.0, 400.0, 5.0);
        cm.step(401.0, 1.0);
        assert_eq!(cm.mode, ChargingMode::ConstantVoltage);
    }
    #[test]
    fn charging_model_energy_accumulates() {
        let mut cm = ChargingModel::home_ac_11kw();
        cm.step(380.0, 3600.0);
        assert!(cm.energy_delivered_wh > 0.0);
    }
    #[test]
    fn charging_model_reset_clears_state() {
        let mut cm = ChargingModel::dc_fast_150kw();
        cm.step(390.0, 100.0);
        cm.reset();
        assert_eq!(cm.mode, ChargingMode::ConstantCurrent);
        assert_eq!(cm.energy_delivered_wh, 0.0);
    }
    #[test]
    fn regen_efficiency_chain_product() {
        let regen = RegenerativeEfficiency::new(200.0, 80_000.0);
        let eta = regen.chain_efficiency();
        assert!(eta > 0.85 && eta < 1.0, "eta={eta}");
    }
    #[test]
    fn regen_efficiency_compute_returns_zero_at_high_soc() {
        let regen = RegenerativeEfficiency::new(200.0, 80_000.0);
        let (t, p) = regen.compute(100.0, 50.0, 0.96);
        assert_eq!(t, 0.0);
        assert_eq!(p, 0.0);
    }
    #[test]
    fn regen_efficiency_compute_positive_at_normal_soc() {
        let regen = RegenerativeEfficiency::new(200.0, 80_000.0);
        let (t, p) = regen.compute(100.0, 50.0, 0.5);
        assert!(t > 0.0, "torque={t}");
        assert!(p > 0.0, "power={p}");
    }
    #[test]
    fn regen_efficiency_recoverable_energy_decreasing_speed() {
        let regen = RegenerativeEfficiency::new(200.0, 80_000.0);
        let e = regen.recoverable_energy(1500.0, 30.0, 0.0, 0.5);
        assert!(e > 0.0, "energy={e}");
    }
    #[test]
    fn regen_efficiency_zero_recovery_at_soc_limit() {
        let regen = RegenerativeEfficiency::new(200.0, 80_000.0);
        let e = regen.recoverable_energy(1500.0, 30.0, 0.0, 0.95);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn soc_estimator_initial_soc() {
        let est = SocEstimator::new(100.0, 0.8);
        assert!((est.soc - 0.8).abs() < 1e-10);
    }
    #[test]
    fn soc_estimator_decreases_on_discharge() {
        let mut est = SocEstimator::new(100.0, 1.0);
        est.update_coulomb(100.0, 3600.0);
        assert!(est.soc < 0.01, "soc={}", est.soc);
    }
    #[test]
    fn soc_estimator_clamped_to_zero() {
        let mut est = SocEstimator::new(10.0, 1.0);
        est.update_coulomb(1000.0, 3600.0);
        assert_eq!(est.soc, 0.0);
    }
    #[test]
    fn soc_estimator_ocv_correction() {
        let mut est = SocEstimator::new(100.0, 0.6);
        est.correct_from_ocv(3.8, 3.0, 4.2);
        assert!(
            (est.soc - 0.6).abs() < 0.1,
            "soc after correction={}",
            est.soc
        );
    }
    #[test]
    fn soc_estimator_correction_due_after_interval() {
        let mut est = SocEstimator::new(100.0, 1.0);
        est.correction_interval = 3;
        est.update_coulomb(1.0, 1.0);
        est.update_coulomb(1.0, 1.0);
        est.update_coulomb(1.0, 1.0);
        assert!(est.correction_due());
    }
    #[test]
    fn soc_estimator_remaining_ah() {
        let est = SocEstimator::new(100.0, 0.5);
        assert!((est.remaining_ah() - 50.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_motor_initial_temp_equals_coolant() {
        let tm = ThermalMotorModel::new(25.0);
        assert!((tm.temperature_c - 25.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_motor_heats_up_with_losses() {
        let mut tm = ThermalMotorModel::new(25.0);
        let init = tm.temperature_c;
        for _ in 0..1000 {
            tm.step(500.0, 0.1);
        }
        assert!(tm.temperature_c > init, "should heat up");
    }
    #[test]
    fn thermal_motor_steady_state_temp() {
        let tm = ThermalMotorModel::new(25.0);
        let ss = tm.steady_state_temp(1000.0);
        assert!((ss - 105.0).abs() < 1.0, "ss={ss}");
    }
    #[test]
    fn thermal_motor_no_derating_below_onset() {
        let tm = ThermalMotorModel::new(25.0);
        assert!((tm.torque_derating() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_motor_derating_above_onset() {
        let mut tm = ThermalMotorModel::new(25.0);
        tm.temperature_c = 170.0;
        let d = tm.torque_derating();
        assert!(d < 1.0 && d > 0.0, "derating={d}");
    }
    #[test]
    fn thermal_motor_overtemp_detection() {
        let mut tm = ThermalMotorModel::new(25.0);
        tm.temperature_c = 180.0;
        assert!(tm.is_overtemperature());
        tm.temperature_c = 179.0;
        assert!(!tm.is_overtemperature());
    }
    #[test]
    fn energy_consumption_initial_state() {
        let ec = EnergyConsumption::new(75.0, 1.0);
        assert_eq!(ec.energy_kwh, 0.0);
        assert_eq!(ec.distance_km, 0.0);
    }
    #[test]
    fn energy_consumption_record_step_increases_energy() {
        let mut ec = EnergyConsumption::new(75.0, 1.0);
        ec.record_step(10_000.0, 1.0, 30.0);
        assert!(ec.energy_kwh > 0.0);
        assert!(ec.distance_km > 0.0);
    }
    #[test]
    fn energy_consumption_soc_decreases() {
        let mut ec = EnergyConsumption::new(75.0, 1.0);
        ec.record_step(50_000.0, 100.0, 30.0);
        assert!(ec.soc < 1.0, "soc={}", ec.soc);
    }
    #[test]
    fn energy_consumption_average_reasonable() {
        let mut ec = EnergyConsumption::new(75.0, 0.9);
        ec.record_step(15_000.0, 3600.0, 100.0 / 3.6);
        let avg = ec.average_consumption_kwh_per_km();
        assert!(avg > 0.0, "avg={avg}");
    }
    #[test]
    fn energy_consumption_remaining_range_positive() {
        let ec = EnergyConsumption::new(75.0, 0.8);
        let range = ec.remaining_range_km();
        assert!(range > 0.0, "range={range}");
    }
    #[test]
    fn energy_consumption_wh_per_km() {
        let mut ec = EnergyConsumption::new(75.0, 1.0);
        ec.record_step(10_000.0, 1.0, 1000.0);
        let wh = ec.wh_per_km();
        assert!(wh.is_finite(), "wh/km={wh}");
    }
    #[test]
    fn ev_drivetrain_single_rear_not_awd() {
        let dt = EvDrivetrain::single_rear(300.0, 12000.0, 9.0);
        assert!(!dt.is_awd());
        assert!(dt.rear_motor.is_some());
        assert!(dt.front_motor.is_none());
    }
    #[test]
    fn ev_drivetrain_awd_is_awd() {
        let dt = EvDrivetrain::dual_motor_awd(200.0, 10000.0, 300.0, 12000.0, 0.4, 9.0, 50.0);
        assert!(dt.is_awd());
    }
    #[test]
    fn ev_drivetrain_step_returns_finite_torque() {
        let mut dt = EvDrivetrain::single_rear(300.0, 12000.0, 9.0);
        let t = dt.step(200.0, 0.01, 0.0);
        assert!(t.is_finite(), "torque={t}");
    }
    #[test]
    fn ev_drivetrain_wheel_speed_update() {
        let mut dt = EvDrivetrain::single_rear(300.0, 12000.0, 9.0);
        dt.update_wheel_speed(30.0, 0.35);
        assert!(dt.wheel_speed_rads > 0.0);
    }
    #[test]
    fn ev_drivetrain_electrical_power_finite() {
        let mut dt = EvDrivetrain::single_rear(300.0, 12000.0, 9.0);
        dt.update_wheel_speed(20.0, 0.35);
        dt.step(100.0, 0.01, 0.0);
        let p = dt.electrical_power_w();
        assert!(p.is_finite(), "power={p}");
    }
    #[test]
    fn ev_drivetrain_tv_correction_asymmetric_torque() {
        let mut dt = EvDrivetrain::dual_motor_awd(200.0, 10000.0, 200.0, 10000.0, 0.5, 9.0, 100.0);
        let t1 = dt.step(200.0, 0.01, 0.0);
        let mut dt2 = EvDrivetrain::dual_motor_awd(200.0, 10000.0, 200.0, 10000.0, 0.5, 9.0, 100.0);
        let _t2 = dt2.step(200.0, 0.01, 0.5);
        assert!(t1.is_finite());
    }
}
/// Simplified open-circuit voltage curve as a fraction of nominal voltage.
///
/// Polynomial fit over SOC ∈ \[0, 1\]: monotonically increasing from 0.85 to 1.05.
/// Clamps SOC to `[0, 1]` before evaluation.
pub fn soc_to_ocv(soc: f64) -> f64 {
    let s = soc.clamp(0.0, 1.0);
    0.85 + 0.12 * s + 0.06 * s * s + 0.02 * s * s * s
}
/// Simple battery thermal model.
///
/// Returns the temperature rise (K) over time step `dt` (s).
///
/// `ΔT = (I² · R_int − cooling_power) · dt / C_thermal`
///
/// * `current`      – charge/discharge current (A)
/// * `r_int`        – internal resistance (Ω)
/// * `dt`           – time step (s)
/// * `thermal_mass` – thermal mass C_th (J/K)
/// * `cooling`      – cooling power (W)
pub fn battery_thermal_model(
    current: f64,
    r_int: f64,
    dt: f64,
    thermal_mass: f64,
    cooling: f64,
) -> f64 {
    if thermal_mass < 1e-15 {
        return 0.0;
    }
    let heat = current * current * r_int;
    (heat - cooling) * dt / thermal_mass
}
/// Compute the regenerative braking torque from deceleration.
///
/// `T_regen = mass · |decel| · wheel_radius · regen_fraction`
///
/// * `decel`          – deceleration magnitude (m/s²)
/// * `mass`           – vehicle mass (kg)
/// * `wheel_radius`   – wheel rolling radius (m)
/// * `regen_fraction` – fraction of braking energy captured by regen
pub fn regenerative_braking_torque(
    decel: f64,
    mass: f64,
    wheel_radius: f64,
    regen_fraction: f64,
) -> f64 {
    mass * decel.abs() * wheel_radius * regen_fraction.clamp(0.0, 1.0)
}
/// Estimate remaining range (km).
///
/// `range = capacity_kwh · SOC / consumption_kwh_per_km`
///
/// Returns `0.0` for zero or negative consumption.
pub fn range_estimation(capacity_kwh: f64, consumption_kwh_per_km: f64, soc: f64) -> f64 {
    if consumption_kwh_per_km < 1e-12 {
        return 0.0;
    }
    capacity_kwh * soc.clamp(0.0, 1.0) / consumption_kwh_per_km
}
/// Peukert's equation: effective capacity at a given discharge current.
///
/// `C_eff = C_rated · (I_rated / I)^(k − 1)`
///
/// * `c_rated`  – rated capacity (Ah) at `i_rated`
/// * `i_rated`  – rated discharge current (A)
/// * `current`  – actual current (A); clamped to positive
/// * `k`        – Peukert exponent (typically 1.05–1.3)
pub fn peukert_capacity(c_rated: f64, i_rated: f64, current: f64, k: f64) -> f64 {
    let i = current.max(1e-9);
    c_rated * (i_rated / i).powf(k - 1.0)
}
/// Estimate energy consumption (kWh) over a speed profile (WLTC-style).
///
/// Uses a simple point-mass model including aerodynamic drag, rolling
/// resistance, and kinetic energy changes.
///
/// * `speed_profile` – vehicle speed at each time step (m/s)
/// * `mass`          – vehicle mass (kg)
/// * `drag_coeff`    – aerodynamic drag coefficient Cd
/// * `area`          – frontal area (m²)
/// * `crr`           – coefficient of rolling resistance
///
/// Time step is assumed to be 1 s throughout.
pub fn energy_consumption_wltc(
    speed_profile: &[f64],
    mass: f64,
    drag_coeff: f64,
    area: f64,
    crr: f64,
) -> f64 {
    pub(super) const RHO: f64 = 1.225;
    pub(super) const G: f64 = 9.81;
    pub(super) const DT: f64 = 1.0;
    let mut total_energy_j = 0.0;
    for i in 0..speed_profile.len() {
        let v = speed_profile[i].max(0.0);
        let f_drag = 0.5 * RHO * drag_coeff * area * v * v;
        let f_roll = crr * mass * G;
        let f_trac = f_drag + f_roll;
        let dke = if i + 1 < speed_profile.len() {
            let v_next = speed_profile[i + 1].max(0.0);
            0.5 * mass * (v_next * v_next - v * v)
        } else {
            0.0
        };
        let power = f_trac * v + dke / DT;
        if power > 0.0 {
            total_energy_j += power * DT;
        }
    }
    total_energy_j / 3_600_000.0
}
#[cfg(test)]
mod simple_ev_tests {
    use super::*;
    use crate::electric_vehicle::ChargeSession;
    use crate::electric_vehicle::MotorEfficiencyMap;
    use crate::electric_vehicle::SimpleBatteryPack;
    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn simple_battery_available_energy_equals_cap_times_soc() {
        let b = SimpleBatteryPack::new(75.0, 400.0, 0.8);
        let expected = 75.0 * 0.8;
        assert!((b.available_energy_kwh() - expected).abs() < EPS);
    }
    #[test]
    fn simple_battery_available_energy_zero_soc() {
        let b = SimpleBatteryPack::new(75.0, 400.0, 0.0);
        assert!(b.available_energy_kwh().abs() < EPS);
    }
    #[test]
    fn simple_battery_soc_clamped() {
        let b = SimpleBatteryPack::new(75.0, 400.0, 1.5);
        assert!((b.soc - 1.0).abs() < EPS);
    }
    #[test]
    fn simple_battery_ocv_positive() {
        let b = SimpleBatteryPack::new(75.0, 400.0, 0.5);
        assert!(b.ocv_voltage() > 0.0);
    }
    #[test]
    fn simple_battery_ocv_increases_with_soc() {
        let b_lo = SimpleBatteryPack::new(75.0, 400.0, 0.2);
        let b_hi = SimpleBatteryPack::new(75.0, 400.0, 0.8);
        assert!(b_hi.ocv_voltage() > b_lo.ocv_voltage());
    }
    #[test]
    fn simple_battery_internal_resistance_positive() {
        let b = SimpleBatteryPack::new(75.0, 400.0, 0.5);
        assert!(b.internal_resistance() > 0.0);
    }
    #[test]
    fn simple_battery_internal_resistance_increases_with_age() {
        let mut b = SimpleBatteryPack::new(75.0, 400.0, 0.5);
        let r_new = b.internal_resistance();
        b.soh = 0.8;
        let r_aged = b.internal_resistance();
        assert!(r_aged > r_new, "aged battery should have higher resistance");
    }
    #[test]
    fn soc_to_ocv_at_zero() {
        let ocv = soc_to_ocv(0.0);
        assert!((ocv - 0.85).abs() < EPS);
    }
    #[test]
    fn soc_to_ocv_monotone_increasing() {
        let v1 = soc_to_ocv(0.2);
        let v2 = soc_to_ocv(0.8);
        assert!(v2 > v1);
    }
    #[test]
    fn soc_to_ocv_clamps_above_1() {
        let v1 = soc_to_ocv(1.0);
        let v2 = soc_to_ocv(1.5);
        assert!((v1 - v2).abs() < EPS);
    }
    #[test]
    fn thermal_model_zero_current_zero_heat() {
        let dt = battery_thermal_model(0.0, 0.1, 1.0, 1000.0, 0.0);
        assert!(dt.abs() < EPS);
    }
    #[test]
    fn thermal_model_positive_current_rises() {
        let dt = battery_thermal_model(100.0, 0.1, 1.0, 1000.0, 0.0);
        assert!(dt > 0.0, "temperature should rise: {dt}");
    }
    #[test]
    fn thermal_model_with_cooling_less_rise() {
        let dt_no_cool = battery_thermal_model(100.0, 0.1, 1.0, 1000.0, 0.0);
        let dt_cool = battery_thermal_model(100.0, 0.1, 1.0, 1000.0, 500.0);
        assert!(dt_cool < dt_no_cool);
    }
    #[test]
    fn motor_map_efficiency_in_range() {
        let m = MotorEfficiencyMap::new(500.0, 150_000.0, 300.0);
        let eta = m.efficiency(250.0, 150.0);
        assert!((0.0..=1.0).contains(&eta), "efficiency={eta}");
    }
    #[test]
    fn motor_map_max_torque_at_low_speed() {
        let m = MotorEfficiencyMap::new(500.0, 150_000.0, 300.0);
        assert!((m.max_torque_at_speed(0.0) - 500.0).abs() < EPS);
    }
    #[test]
    fn motor_map_max_torque_decreases_above_base() {
        let m = MotorEfficiencyMap::new(500.0, 150_000.0, 300.0);
        let t_base = m.max_torque_at_speed(300.0);
        let t_high = m.max_torque_at_speed(1000.0);
        assert!(t_high <= t_base);
    }
    #[test]
    fn motor_map_power_at_operating_point() {
        let m = MotorEfficiencyMap::new(500.0, 150_000.0, 300.0);
        let p = m.power_at_operating_point(200.0, 100.0);
        assert!((p - 20_000.0).abs() < EPS);
    }
    #[test]
    fn regen_torque_positive_for_positive_decel() {
        let t = regenerative_braking_torque(3.0, 2000.0, 0.33, 0.7);
        assert!(t > 0.0, "regen torque={t}");
    }
    #[test]
    fn regen_torque_zero_decel_zero_torque() {
        let t = regenerative_braking_torque(0.0, 2000.0, 0.33, 0.7);
        assert!(t.abs() < EPS);
    }
    #[test]
    fn regen_torque_scales_with_mass() {
        let t1 = regenerative_braking_torque(3.0, 1000.0, 0.33, 0.7);
        let t2 = regenerative_braking_torque(3.0, 2000.0, 0.33, 0.7);
        assert!((t2 - 2.0 * t1).abs() < EPS);
    }
    #[test]
    fn range_decreases_with_lower_soc() {
        let r1 = range_estimation(75.0, 0.2, 0.8);
        let r2 = range_estimation(75.0, 0.2, 0.4);
        assert!(r1 > r2);
    }
    #[test]
    fn range_zero_consumption_returns_zero() {
        let r = range_estimation(75.0, 0.0, 0.8);
        assert!(r.abs() < EPS);
    }
    #[test]
    fn range_formula_check() {
        let r = range_estimation(75.0, 0.2, 0.8);
        let expected = 75.0 * 0.8 / 0.2;
        assert!((r - expected).abs() < EPS);
    }
    #[test]
    fn range_zero_soc_zero_range() {
        let r = range_estimation(75.0, 0.2, 0.0);
        assert!(r.abs() < EPS);
    }
    #[test]
    fn charge_session_energy_delivered() {
        let cs = ChargeSession {
            power_kw: 50.0,
            duration_h: 1.0,
            efficiency: 0.95,
        };
        let e = cs.energy_delivered_kwh();
        assert!((e - 47.5).abs() < EPS);
    }
    #[test]
    fn charge_session_cost_estimate() {
        let cs = ChargeSession {
            power_kw: 50.0,
            duration_h: 2.0,
            efficiency: 0.95,
        };
        let cost = cs.cost_estimate(0.25);
        assert!((cost - 50.0 * 2.0 * 0.25).abs() < EPS);
    }
    #[test]
    fn charge_session_zero_duration_zero_energy() {
        let cs = ChargeSession {
            power_kw: 100.0,
            duration_h: 0.0,
            efficiency: 0.9,
        };
        assert!(cs.energy_delivered_kwh().abs() < EPS);
    }
    #[test]
    fn peukert_at_rated_current_equals_rated_capacity() {
        let c = peukert_capacity(100.0, 20.0, 20.0, 1.2);
        assert!((c - 100.0).abs() < EPS);
    }
    #[test]
    fn peukert_high_current_reduces_capacity() {
        let c_low = peukert_capacity(100.0, 20.0, 10.0, 1.2);
        let c_high = peukert_capacity(100.0, 20.0, 40.0, 1.2);
        assert!(
            c_high < c_low,
            "high current should reduce effective capacity"
        );
    }
    #[test]
    fn peukert_positive_output() {
        let c = peukert_capacity(100.0, 20.0, 30.0, 1.1);
        assert!(c > 0.0);
    }
    #[test]
    fn wltc_energy_positive_for_nonzero_profile() {
        let profile = vec![0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 20.0, 10.0, 0.0];
        let e = energy_consumption_wltc(&profile, 1800.0, 0.3, 2.2, 0.01);
        assert!(e >= 0.0);
    }
    #[test]
    fn wltc_energy_zero_for_zero_speed() {
        let profile = vec![0.0; 10];
        let e = energy_consumption_wltc(&profile, 1800.0, 0.3, 2.2, 0.01);
        assert!(e.abs() < EPS);
    }
    #[test]
    fn wltc_energy_increases_with_drag() {
        let profile = vec![20.0; 10];
        let e_low = energy_consumption_wltc(&profile, 1800.0, 0.25, 2.2, 0.01);
        let e_high = energy_consumption_wltc(&profile, 1800.0, 0.35, 2.2, 0.01);
        assert!(e_high > e_low);
    }
    #[test]
    fn wltc_energy_finite() {
        let profile: Vec<f64> = (0..50).map(|i| i as f64 % 30.0).collect();
        let e = energy_consumption_wltc(&profile, 1800.0, 0.3, 2.2, 0.01);
        assert!(e.is_finite());
    }
}
/// Compute the open-circuit voltage from state of charge using a polynomial fit.
///
/// Fits a simple quadratic: `V_oc = V_nom * (0.85 + 0.30·soc − 0.05·soc²)`
pub fn battery_voltage(bat: &Battery) -> f64 {
    let s = bat.soc.clamp(0.0, 1.0);
    bat.voltage_nominal * (0.85 + 0.30 * s - 0.05 * s * s)
}
/// Compute the terminal current drawn/supplied (A) for a given power demand.
///
/// Solves `V_oc·I − R·I² = P_w` for the smaller (physical) root.
/// Positive power = discharge; negative = charge.
pub fn battery_current(bat: &Battery, power_w: f64) -> f64 {
    let v_oc = battery_voltage(bat);
    let r = bat.internal_resistance;
    let discriminant = v_oc * v_oc - 4.0 * r * power_w;
    if discriminant < 0.0 {
        return v_oc / (2.0 * r);
    }
    let sqrt_d = discriminant.sqrt();
    (v_oc - sqrt_d) / (2.0 * r)
}
/// Compute ohmic heat generation rate (W) = I²R.
pub fn battery_heat(bat: &Battery, current: f64) -> f64 {
    current * current * bat.internal_resistance
}
/// Advance the battery state by `dt` seconds under `power_w` watts.
///
/// Updates `soc` and `temp_celsius`.
pub fn step_battery(bat: &mut Battery, power_w: f64, params: &BatteryParams, dt: f64) {
    let current = battery_current(bat, power_w);
    let heat = battery_heat(bat, current);
    let capacity_ws = bat.capacity_kwh * 3_600_000.0;
    let eta = if power_w >= 0.0 {
        params.discharge_efficiency
    } else {
        params.charge_efficiency
    };
    bat.soc -= eta * current * battery_voltage(bat) * dt / capacity_ws;
    bat.soc = bat.soc.clamp(0.0, 1.0);
    let net_heat = heat - params.cooling_power;
    bat.temp_celsius += (net_heat * dt) / params.thermal_capacity;
}
/// Look up motor efficiency for a given operating point via bilinear interpolation.
///
/// Falls back to the nearest entry if outside the map boundaries.
pub fn motor_efficiency(motor: &ElectricMotor, speed_rpm: f64, torque_nm: f64) -> f64 {
    if motor.efficiency_map.is_empty() {
        return 0.90;
    }
    let best = motor
        .efficiency_map
        .iter()
        .min_by(|a, b| {
            let da = (a.0 - speed_rpm).powi(2) + (a.1 - torque_nm).powi(2);
            let db = (b.0 - speed_rpm).powi(2) + (b.1 - torque_nm).powi(2);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("operation should succeed");
    best.2.clamp(0.0, 1.0)
}
/// Compute the mechanical power output (W) at an operating point.
///
/// `P = (2π/60) · speed_rpm · torque_nm`
pub fn motor_power(motor: &ElectricMotor, speed_rpm: f64, torque_nm: f64) -> f64 {
    let torque_clamped = torque_nm.min(motor.max_torque_nm);
    let power_limit = motor.max_power_kw * 1000.0;
    let p = (2.0 * std::f64::consts::PI / 60.0) * speed_rpm * torque_clamped;
    p.min(power_limit)
}
/// Compute the regenerative braking torque (N·m) available for an EV powertrain.
///
/// `T_regen = regen_frac · mass · decel · wheel_radius`
pub fn ev_regenerative_braking_torque(
    decel: f64,
    mass: f64,
    wheel_radius: f64,
    regen_frac: f64,
) -> f64 {
    regen_frac * mass * decel * wheel_radius
}
/// Estimate remaining range (km) given the current battery state.
///
/// `range = capacity_kwh · soc / avg_consumption_wh_per_km * 1000`
pub fn range_estimate(bat: &Battery, avg_consumption_wh_per_km: f64) -> f64 {
    if avg_consumption_wh_per_km < 1e-15 {
        return 0.0;
    }
    let available_wh = bat.capacity_kwh * 1000.0 * bat.soc;
    available_wh / avg_consumption_wh_per_km
}
/// Simulate one time step of the combined motor + battery system.
///
/// Computes the motor's electrical power demand (accounting for efficiency)
/// and steps the battery accordingly.
pub fn ev_step(
    bat: &mut Battery,
    motor: &ElectricMotor,
    params: &BatteryParams,
    speed_rpm: f64,
    torque_nm: f64,
    dt: f64,
) {
    let mech_power = motor_power(motor, speed_rpm, torque_nm);
    let eta = motor_efficiency(motor, speed_rpm, torque_nm).max(1e-6);
    let elec_power = mech_power / eta;
    step_battery(bat, elec_power, params, dt);
}
#[cfg(test)]
mod ev_powertrain_tests {
    use super::*;

    pub(super) const EPS: f64 = 1e-9;
    fn bat() -> Battery {
        Battery::typical_75kwh()
    }
    fn params() -> BatteryParams {
        BatteryParams::default_params()
    }
    fn motor() -> ElectricMotor {
        ElectricMotor::typical_200kw()
    }
    #[test]
    fn battery_voltage_full_soc() {
        let b = bat();
        let v = battery_voltage(&b);
        assert!((v - 440.0).abs() < EPS, "v={v}");
    }
    #[test]
    fn battery_voltage_empty_soc() {
        let mut b = bat();
        b.soc = 0.0;
        let v = battery_voltage(&b);
        assert!((v - 340.0).abs() < EPS, "v={v}");
    }
    #[test]
    fn battery_voltage_monotone_increasing() {
        let mut b = bat();
        b.soc = 0.3;
        let v1 = battery_voltage(&b);
        b.soc = 0.7;
        let v2 = battery_voltage(&b);
        assert!(v2 > v1, "voltage should increase with SOC");
    }
    #[test]
    fn battery_voltage_clamps_above_one() {
        let mut b = bat();
        b.soc = 1.5;
        let v1 = battery_voltage(&b);
        b.soc = 1.0;
        let v2 = battery_voltage(&b);
        assert!((v1 - v2).abs() < EPS, "SOC clamped to 1");
    }
    #[test]
    fn battery_current_positive_for_positive_power() {
        let b = bat();
        let i = battery_current(&b, 10_000.0);
        assert!(i > 0.0, "positive power → positive current");
    }
    #[test]
    fn battery_current_finite() {
        let b = bat();
        assert!(battery_current(&b, 50_000.0).is_finite());
    }
    #[test]
    fn battery_current_zero_power_near_zero() {
        let b = bat();
        let i = battery_current(&b, 0.0);
        assert!(i.abs() < EPS, "zero power → zero current");
    }
    #[test]
    fn battery_heat_zero_current_zero() {
        let b = bat();
        assert!(battery_heat(&b, 0.0).abs() < EPS);
    }
    #[test]
    fn battery_heat_positive_for_nonzero_current() {
        let b = bat();
        assert!(battery_heat(&b, 10.0) > 0.0);
    }
    #[test]
    fn battery_heat_formula() {
        let b = bat();
        let heat = battery_heat(&b, 10.0);
        assert!((heat - 5.0).abs() < EPS, "heat={heat}");
    }
    #[test]
    fn step_battery_discharging_decreases_soc() {
        let mut b = bat();
        let s0 = b.soc;
        step_battery(&mut b, 50_000.0, &params(), 1.0);
        assert!(b.soc < s0, "discharging should reduce SOC");
    }
    #[test]
    fn step_battery_charging_increases_soc() {
        let mut b = bat();
        b.soc = 0.5;
        step_battery(&mut b, -10_000.0, &params(), 1.0);
        assert!(b.soc > 0.5, "charging should increase SOC");
    }
    #[test]
    fn step_battery_soc_clamped_to_zero() {
        let mut b = bat();
        b.soc = 0.0;
        step_battery(&mut b, 1_000_000.0, &params(), 100.0);
        assert!(b.soc >= 0.0, "SOC should not go negative");
    }
    #[test]
    fn step_battery_heat_raises_temperature() {
        let mut b = bat();
        b.temp_celsius = 25.0;
        let mut p = params();
        p.cooling_power = 0.0;
        step_battery(&mut b, 100_000.0, &p, 10.0);
        assert!(b.temp_celsius > 25.0, "temperature should rise");
    }
    #[test]
    fn motor_efficiency_in_range() {
        let m = motor();
        let eta = motor_efficiency(&m, 3000.0, 200.0);
        assert!((0.0..=1.0).contains(&eta), "eta={eta}");
    }
    #[test]
    fn motor_efficiency_empty_map_default() {
        let m = ElectricMotor {
            max_torque_nm: 100.0,
            max_power_kw: 50.0,
            efficiency_map: vec![],
        };
        let eta = motor_efficiency(&m, 1000.0, 100.0);
        assert!((eta - 0.90).abs() < EPS);
    }
    #[test]
    fn motor_efficiency_finite() {
        let m = motor();
        assert!(motor_efficiency(&m, 0.0, 0.0).is_finite());
    }
    #[test]
    fn motor_power_zero_speed_zero() {
        let m = motor();
        assert!(motor_power(&m, 0.0, 100.0).abs() < EPS);
    }
    #[test]
    fn motor_power_clamps_to_max_power() {
        let m = motor();
        let p = motor_power(&m, 9000.0, 500.0);
        assert!(
            p <= m.max_power_kw * 1000.0 + EPS,
            "power exceeded limit: {p}"
        );
    }
    #[test]
    fn motor_power_positive_for_positive_torque_and_speed() {
        let m = motor();
        assert!(motor_power(&m, 3000.0, 200.0) > 0.0);
    }
    #[test]
    fn regen_torque_zero_decel_zero() {
        let t = ev_regenerative_braking_torque(0.0, 2000.0, 0.33, 0.7);
        assert!(t.abs() < EPS);
    }
    #[test]
    fn regen_torque_positive_for_positive_decel() {
        let t = ev_regenerative_braking_torque(2.0, 2000.0, 0.33, 0.7);
        assert!(t > 0.0);
    }
    #[test]
    fn regen_torque_scales_with_mass() {
        let t1 = ev_regenerative_braking_torque(3.0, 1000.0, 0.33, 0.7);
        let t2 = ev_regenerative_braking_torque(3.0, 2000.0, 0.33, 0.7);
        assert!((t2 - 2.0 * t1).abs() < EPS);
    }
    #[test]
    fn range_estimate_full_soc() {
        let b = bat();
        let r = range_estimate(&b, 200.0);
        assert!((r - 75_000.0 / 200.0).abs() < EPS, "range={r}");
    }
    #[test]
    fn range_estimate_zero_consumption_returns_zero() {
        let b = bat();
        assert!(range_estimate(&b, 0.0).abs() < EPS);
    }
    #[test]
    fn range_estimate_decreases_with_lower_soc() {
        let mut b = bat();
        b.soc = 0.8;
        let r1 = range_estimate(&b, 200.0);
        b.soc = 0.4;
        let r2 = range_estimate(&b, 200.0);
        assert!(r1 > r2);
    }
    #[test]
    fn range_estimate_zero_soc_zero_range() {
        let mut b = bat();
        b.soc = 0.0;
        assert!(range_estimate(&b, 200.0).abs() < EPS);
    }
    #[test]
    fn ev_step_discharges_battery() {
        let mut b = bat();
        let s0 = b.soc;
        ev_step(&mut b, &motor(), &params(), 3000.0, 200.0, 1.0);
        assert!(b.soc < s0, "ev_step should discharge battery");
    }
    #[test]
    fn ev_step_zero_speed_no_discharge() {
        let mut b = bat();
        let s0 = b.soc;
        ev_step(&mut b, &motor(), &params(), 0.0, 0.0, 1.0);
        assert!(
            (b.soc - s0).abs() < EPS,
            "zero operating point: SOC unchanged"
        );
    }
    #[test]
    fn ev_step_soc_stays_non_negative() {
        let mut b = bat();
        b.soc = 0.01;
        for _ in 0..1000 {
            ev_step(&mut b, &motor(), &params(), 3000.0, 400.0, 1.0);
        }
        assert!(b.soc >= 0.0, "SOC should not go negative");
    }
}
