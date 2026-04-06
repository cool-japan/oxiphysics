//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::fuel_systems::CombustionModel;
    use crate::fuel_systems::EmissionsModel;
    use crate::fuel_systems::EnergyRecovery;
    use crate::fuel_systems::ExhaustSystem;
    use crate::fuel_systems::FuelInjector;
    use crate::fuel_systems::FuelProperties;
    use crate::fuel_systems::FuelPump;
    use crate::fuel_systems::FuelTank;
    use crate::fuel_systems::FuelTankCoM;
    use crate::fuel_systems::HybridBattery;
    use crate::fuel_systems::HybridEnergyManagement;
    use crate::fuel_systems::KersSystem;
    use crate::fuel_systems::RegenerativeBraking;
    use crate::fuel_systems::TersSystem;
    use std::f64::consts::PI;
    #[test]
    fn fuel_tank_initial_mass_correct() {
        let tank = FuelTank::new_gasoline(60.0, 0.5);
        let expected = 60.0 * 0.5 * 0.745;
        assert!(
            (tank.mass_kg() - expected).abs() < 1e-6,
            "mass_kg={}, expected={expected}",
            tank.mass_kg()
        );
    }
    #[test]
    fn fuel_tank_fill_fraction_clamped() {
        let tank = FuelTank::new_gasoline(60.0, 1.5);
        assert_eq!(tank.fill_fraction(), 1.0);
    }
    #[test]
    fn fuel_tank_consume_reduces_level() {
        let mut tank = FuelTank::new_gasoline(60.0, 1.0);
        let m0 = tank.mass_kg();
        let consumed = tank.consume(5.0);
        assert!((consumed - 5.0).abs() < 1e-6);
        assert!(tank.mass_kg() < m0, "level should decrease after consuming");
    }
    #[test]
    fn fuel_tank_consume_clamped_to_available() {
        let mut tank = FuelTank::new_gasoline(10.0, 0.1);
        let available = tank.mass_kg();
        let consumed = tank.consume(100.0);
        assert!(
            (consumed - available).abs() < 1e-6,
            "consumed={consumed}, available={available}"
        );
    }
    #[test]
    fn fuel_tank_refuel_clamps_to_capacity() {
        let mut tank = FuelTank::new_gasoline(60.0, 0.5);
        tank.refuel(100.0);
        assert!((tank.level_l - 60.0).abs() < 1e-6);
    }
    #[test]
    fn fuel_tank_is_empty_after_full_consume() {
        let mut tank = FuelTank::new_diesel(20.0, 1.0);
        tank.consume(1000.0);
        assert!(tank.is_empty());
    }
    #[test]
    fn fuel_tank_slosh_angle_changes_with_lateral_accel() {
        let mut tank = FuelTank::new_gasoline(50.0, 0.8);
        let angle0 = tank.slosh_angle;
        for _ in 0..50 {
            tank.integrate_slosh(5.0, 0.02);
        }
        assert!(
            tank.slosh_angle != angle0,
            "slosh angle should change under lateral accel"
        );
    }
    #[test]
    fn fuel_tank_slosh_angle_bounded() {
        let mut tank = FuelTank::new_gasoline(50.0, 0.8);
        for _ in 0..1000 {
            tank.integrate_slosh(50.0, 0.02);
        }
        assert!(
            tank.slosh_angle.abs() <= PI / 3.0 + 1e-9,
            "slosh angle={}",
            tank.slosh_angle
        );
    }
    #[test]
    fn fuel_injector_gdi_pulse_mass_positive() {
        let inj = FuelInjector::new_gdi();
        assert!(
            inj.pulse_mass_mg() > 0.0,
            "GDI pulse mass should be positive"
        );
    }
    #[test]
    fn fuel_injector_pfi_pulse_mass_positive() {
        let inj = FuelInjector::new_pfi();
        assert!(
            inj.pulse_mass_mg() > 0.0,
            "PFI pulse mass should be positive"
        );
    }
    #[test]
    fn fuel_injector_smd_positive() {
        let inj = FuelInjector::new_gdi();
        let smd = inj.sauter_mean_diameter_um();
        assert!(
            smd > 0.0 && smd < 500.0,
            "SMD={smd} µm should be in realistic range"
        );
    }
    #[test]
    fn fuel_injector_fire_updates_injected_mass() {
        let mut inj = FuelInjector::new_gdi();
        inj.fire();
        assert!(inj.injected_mass_mg > 0.0);
    }
    #[test]
    fn fuel_injector_reset_cycle_clears_mass() {
        let mut inj = FuelInjector::new_gdi();
        inj.fire();
        inj.reset_cycle();
        assert_eq!(inj.injected_mass_mg, 0.0);
    }
    #[test]
    fn fuel_injector_spray_cone_angle_reasonable() {
        let inj = FuelInjector::new_gdi();
        let angle = inj.spray_cone_angle_deg();
        assert!(angle > 10.0 && angle < 45.0, "angle={angle}");
    }
    #[test]
    fn fuel_injector_gdi_higher_pressure_than_pfi() {
        let gdi = FuelInjector::new_gdi();
        let pfi = FuelInjector::new_pfi();
        assert!(gdi.rail_pressure_bar > pfi.rail_pressure_bar);
    }
    #[test]
    fn combustion_model_heat_release_positive() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        assert!(model.heat_release_j() > 0.0);
    }
    #[test]
    fn combustion_model_lambda_stoich() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        let lam = model.lambda();
        assert!(
            (lam - 1.0).abs() < 1e-6,
            "lambda={lam} should be 1.0 at stoich"
        );
    }
    #[test]
    fn combustion_model_wiebe_fraction_zero_before_start() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        assert_eq!(model.wiebe_fraction(-20.0), 0.0);
    }
    #[test]
    fn combustion_model_wiebe_fraction_approaches_one_at_end() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        let end_angle = model.wiebe_start_deg + model.wiebe_duration_deg;
        let x = model.wiebe_fraction(end_angle);
        assert!(x > 0.99, "wiebe fraction at end should approach 1, got {x}");
    }
    #[test]
    fn combustion_model_nox_hot_engine() {
        let mut model = CombustionModel::new_gasoline_si(5000.0, 50.0);
        model.air_mass_mg = 50.0 * 14.7;
        let nox = model.nox_mg_per_cycle();
        assert!(nox >= 0.0, "NOx should be non-negative");
    }
    #[test]
    fn combustion_model_co_increases_rich() {
        let mut model_rich = CombustionModel::new_gasoline_si(3000.0, 20.0);
        model_rich.air_mass_mg = 20.0 * 14.7 * 0.85;
        let co_rich = model_rich.co_mg_per_cycle();
        let model_stoich = CombustionModel::new_gasoline_si(3000.0, 20.0);
        let co_stoich = model_stoich.co_mg_per_cycle();
        assert!(co_rich > co_stoich, "CO should be higher for rich mixture");
    }
    #[test]
    fn combustion_model_peak_temperature_positive() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        assert!(model.peak_temperature_k() > 300.0);
    }
    #[test]
    fn combustion_model_heat_release_rate_zero_outside_window() {
        let model = CombustionModel::new_gasoline_si(3000.0, 20.0);
        assert_eq!(model.heat_release_rate(-90.0), 0.0);
        assert_eq!(model.heat_release_rate(100.0), 0.0);
    }
    #[test]
    fn exhaust_system_pipe_outlet_temp_less_than_manifold() {
        let ex = ExhaustSystem::new_passenger_car();
        assert!(
            ex.pipe_outlet_temp_k() < ex.manifold_temp_k,
            "pipe outlet should be cooler than manifold"
        );
    }
    #[test]
    fn exhaust_system_backpressure_increases_with_flow() {
        let ex = ExhaustSystem::new_passenger_car();
        let bp1 = ex.backpressure_pa(0.03);
        let bp2 = ex.backpressure_pa(0.05);
        assert!(bp2 > bp1, "backpressure should increase with flow");
    }
    #[test]
    fn exhaust_system_cold_catalyst_no_conversion() {
        let mut ex = ExhaustSystem::new_passenger_car();
        ex.catalyst.temperature_k = 300.0;
        let pass_through = ex.emission_pass_through(1.0);
        assert!(
            (pass_through - 1.0).abs() < 1e-6,
            "cold catalyst should pass all emissions"
        );
    }
    #[test]
    fn exhaust_system_hot_catalyst_converts_emissions() {
        let mut ex = ExhaustSystem::new_passenger_car();
        ex.catalyst.temperature_k = 800.0;
        let pass_through = ex.emission_pass_through(1.0);
        assert!(
            pass_through < 0.5,
            "hot catalyst should convert most emissions"
        );
    }
    #[test]
    fn exhaust_system_update_warms_catalyst() {
        let mut ex = ExhaustSystem::new_passenger_car();
        let t0 = ex.catalyst.temperature_k;
        for _ in 0..100 {
            ex.update(1.0);
        }
        assert!(
            ex.catalyst.temperature_k > t0 || ex.catalyst.temperature_k > 400.0,
            "catalyst should warm up"
        );
    }
    #[test]
    fn exhaust_system_catalyst_lit_off_when_hot() {
        let mut ex = ExhaustSystem::new_passenger_car();
        ex.catalyst.temperature_k = 700.0;
        assert!(ex.catalyst_lit_off());
    }
    #[test]
    fn hybrid_battery_initial_soc_correct() {
        let bat = HybridBattery::new_nmc(5.0, 0.8);
        assert!((bat.soc - 0.8).abs() < 1e-6);
    }
    #[test]
    fn hybrid_battery_ocv_increases_with_soc() {
        let mut bat = HybridBattery::new_nmc(5.0, 0.1);
        let v1 = bat.open_circuit_voltage();
        bat.soc = 0.9;
        let v2 = bat.open_circuit_voltage();
        assert!(v2 > v1, "OCV should increase with SOC");
    }
    #[test]
    fn hybrid_battery_charge_increases_soc() {
        let mut bat = HybridBattery::new_nmc(5.0, 0.5);
        let soc0 = bat.soc;
        bat.charge(10.0, 60.0);
        assert!(bat.soc > soc0, "SOC should increase when charging");
    }
    #[test]
    fn hybrid_battery_discharge_decreases_soc() {
        let mut bat = HybridBattery::new_nmc(5.0, 0.8);
        let soc0 = bat.soc;
        bat.discharge(10.0, 60.0);
        assert!(bat.soc < soc0, "SOC should decrease when discharging");
    }
    #[test]
    fn hybrid_battery_soc_clamped_to_one_when_overcharged() {
        let mut bat = HybridBattery::new_nmc(1.0, 0.99);
        bat.charge(1e6, 1.0);
        assert!(bat.soc <= 1.0, "SOC must not exceed 1.0");
    }
    #[test]
    fn hybrid_battery_terminal_voltage_less_than_ocv_under_load() {
        let bat = HybridBattery::new_nmc(5.0, 0.8);
        let v_term = bat.terminal_voltage(5.0);
        let v_oc = bat.open_circuit_voltage();
        assert!(
            v_term < v_oc,
            "terminal voltage should be less than OCV under load"
        );
    }
    #[test]
    fn hybrid_battery_soh_starts_at_one() {
        let bat = HybridBattery::new_nmc(5.0, 0.5);
        assert!((bat.state_of_health - 1.0).abs() < 1e-6);
    }
    #[test]
    fn hybrid_battery_thermal_update_changes_temperature() {
        let mut bat = HybridBattery::new_nmc(5.0, 0.5);
        let t0 = bat.temperature_c;
        bat.update_thermal(500.0, 25.0, 10.0);
        assert!(
            bat.temperature_c != t0,
            "temperature should change after thermal update"
        );
    }
    #[test]
    fn kers_absorb_increases_stored_energy() {
        let mut kers = KersSystem::new(0.25, 60_000.0 * PI / 30.0);
        kers.omega_rad_s = 1000.0;
        let e0 = kers.stored_energy_j;
        kers.absorb_energy(50.0, 0.1);
        assert!(
            kers.stored_energy_j >= e0,
            "stored energy should increase or stay same"
        );
    }
    #[test]
    fn kers_release_reduces_speed() {
        let mut kers = KersSystem::new(0.25, 60_000.0 * PI / 30.0);
        kers.omega_rad_s = 5000.0;
        kers.stored_energy_j = 0.5 * 0.25 * 5000.0 * 5000.0;
        let omega0 = kers.omega_rad_s;
        kers.release_energy(100.0, 0.1);
        assert!(kers.omega_rad_s < omega0, "flywheel should slow down");
    }
    #[test]
    fn kers_soc_between_zero_and_one() {
        let kers = KersSystem::new(0.25, 60_000.0 * PI / 30.0);
        let soc = kers.state_of_charge();
        assert!((0.0..=1.0).contains(&soc), "SOC={soc}");
    }
    #[test]
    fn ters_carnot_efficiency_between_zero_and_one() {
        let ters = TersSystem::new_teg(1000.0);
        let eta = ters.carnot_efficiency();
        assert!((0.0..=1.0).contains(&eta), "eta={eta}");
    }
    #[test]
    fn ters_recoverable_power_positive() {
        let ters = TersSystem::new_teg(1000.0);
        let p = ters.recoverable_power_w(0.05);
        assert!(p >= 0.0, "recoverable power should be non-negative: {p}");
    }
    #[test]
    fn ters_teg_max_power_positive() {
        let ters = TersSystem::new_teg(1000.0);
        let p = ters.teg_max_power_w();
        assert!(p > 0.0, "TEG max power should be positive: {p}");
    }
    #[test]
    fn regen_braking_zero_below_min_speed() {
        let regen = RegenerativeBraking::new(200.0, 0.6, 0.85);
        let torque = regen.regen_torque(500.0, 0.5);
        assert_eq!(torque, 0.0);
    }
    #[test]
    fn regen_braking_capped_at_max() {
        let regen = RegenerativeBraking::new(200.0, 0.6, 0.85);
        let torque = regen.regen_torque(10_000.0, 50.0);
        assert!(
            torque <= 200.0 + 1e-9,
            "regen torque should be capped: {torque}"
        );
    }
    #[test]
    fn energy_recovery_step_does_not_panic() {
        let mut er = EnergyRecovery::new();
        er.step(200.0, 50.0, 0.04, 0.01);
    }
    #[test]
    fn energy_recovery_total_stored_energy_positive() {
        let er = EnergyRecovery::new();
        assert!(er.total_stored_energy_j() > 0.0);
    }
    #[test]
    fn fuel_properties_gasoline_stoich_afr() {
        let fp = FuelProperties::gasoline();
        assert!((fp.stoich_afr - 14.7).abs() < 0.1);
    }
    #[test]
    fn fuel_pump_gdi_flow_at_zero_back_pressure() {
        let pump = FuelPump::new_gdi_hp();
        let flow = pump.flow_rate_ml_per_min(0.0_f64);
        assert!(
            flow > 0.0_f64,
            "flow at zero back-pressure must be positive: {flow}"
        );
    }
    #[test]
    fn fuel_pump_gdi_zero_flow_at_max_pressure() {
        let pump = FuelPump::new_gdi_hp();
        let flow = pump.flow_rate_ml_per_min(pump.max_pressure_bar);
        assert!(
            flow < 1.0_f64,
            "flow should be near zero at max pressure: {flow}"
        );
    }
    #[test]
    fn fuel_pump_flow_decreases_with_back_pressure() {
        let pump = FuelPump::new_gdi_hp();
        let q1 = pump.flow_rate_ml_per_min(50.0_f64);
        let q2 = pump.flow_rate_ml_per_min(150.0_f64);
        assert!(
            q1 > q2,
            "flow should decrease with back-pressure: q1={q1} q2={q2}"
        );
    }
    #[test]
    fn fuel_pump_hydraulic_power_positive() {
        let pump = FuelPump::new_gdi_hp();
        let p = pump.hydraulic_power_w(100.0_f64);
        assert!(p > 0.0_f64, "hydraulic power must be positive: {p}");
    }
    #[test]
    fn fuel_pump_electric_power_positive() {
        let pump = FuelPump::new_transfer();
        assert!(pump.electric_power_w() > 0.0_f64);
    }
    #[test]
    fn fuel_pump_overall_efficiency_in_range() {
        let pump = FuelPump::new_gdi_hp();
        let eff = pump.overall_efficiency(100.0_f64);
        assert!((0.0_f64..=1.0_f64).contains(&eff), "efficiency={eff}");
    }
    #[test]
    fn emissions_model_firing_rate_positive() {
        let em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2500.0_f64);
        assert!(em.firing_rate_hz() > 0.0_f64);
    }
    #[test]
    fn emissions_model_nox_nonnegative() {
        let em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2800.0_f64);
        assert!(em.nox_raw_g_s() >= 0.0_f64);
    }
    #[test]
    fn emissions_model_nox_higher_at_hotter_temp() {
        let em_cool =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 1800.0_f64);
        let em_hot =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2800.0_f64);
        assert!(
            em_hot.nox_raw_g_s() >= em_cool.nox_raw_g_s(),
            "NOx should increase with peak temperature"
        );
    }
    #[test]
    fn emissions_model_co_positive_rich() {
        let mut em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2500.0_f64);
        em.lambda = 0.85_f64;
        assert!(em.co_raw_g_s() > 0.0_f64);
    }
    #[test]
    fn emissions_model_tailpipe_nox_leq_raw() {
        let em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2800.0_f64);
        assert!(em.nox_tailpipe_g_s() <= em.nox_raw_g_s() + 1.0e-15_f64);
    }
    #[test]
    fn emissions_model_pm_increases_rich() {
        let em_lean =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.1_f64, 2200.0_f64);
        let mut em_rich =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 0.85_f64, 2200.0_f64);
        em_rich.lambda = 0.85_f64;
        assert!(
            em_rich.pm_raw_g_s() > em_lean.pm_raw_g_s(),
            "PM should be higher for rich mixtures"
        );
    }
    #[test]
    fn emissions_model_fuel_flow_positive() {
        let em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2500.0_f64);
        assert!(em.fuel_flow_g_s() > 0.0_f64);
    }
    #[test]
    fn emissions_model_bsfc_realistic() {
        let em =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2500.0_f64);
        let bsfc = em.bsfc_g_per_kwh();
        assert!(bsfc > 150.0_f64 && bsfc < 500.0_f64, "BSFC={bsfc} g/kWh");
    }
    #[test]
    fn emissions_model_egr_reduces_nox() {
        let mut em_no_egr =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2800.0_f64);
        em_no_egr.egr_fraction = 0.0_f64;
        let mut em_egr =
            EmissionsModel::new_gasoline_na(2.0_f64, 4, 3000.0_f64, 20.0_f64, 1.0_f64, 2800.0_f64);
        em_egr.egr_fraction = 0.2_f64;
        assert!(
            em_egr.nox_raw_g_s() < em_no_egr.nox_raw_g_s(),
            "EGR should reduce NOx"
        );
    }
    #[test]
    fn hem_power_split_sums_to_demand() {
        let hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.6_f64);
        let demand = 60_000.0_f64;
        let (ice, em) = hem.power_split(demand);
        let total = ice + em;
        assert!(
            (total - demand).abs() <= demand.abs() + 1.0_f64,
            "Power split sum={total} demand={demand}"
        );
    }
    #[test]
    fn hem_ev_mode_available_high_soc() {
        let mut hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.75_f64);
        hem.soc = 0.75_f64;
        assert!(hem.ev_mode_available());
    }
    #[test]
    fn hem_ev_mode_unavailable_low_soc() {
        let mut hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.2_f64);
        hem.soc = 0.2_f64;
        assert!(!hem.ev_mode_available());
    }
    #[test]
    fn hem_update_soc_charging_increases_soc() {
        let mut hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.5_f64);
        let soc0 = hem.soc;
        hem.update_soc(-10_000.0_f64, 10.0_f64);
        assert!(hem.soc > soc0, "SOC should increase when generating");
    }
    #[test]
    fn hem_update_soc_discharging_decreases_soc() {
        let mut hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.5_f64);
        let soc0 = hem.soc;
        hem.update_soc(10_000.0_f64, 10.0_f64);
        assert!(hem.soc < soc0, "SOC should decrease when motoring");
    }
    #[test]
    fn hem_fuel_consumption_positive_under_load() {
        let hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.6_f64);
        let fc = hem.fuel_consumption_g_s(50_000.0_f64);
        assert!(
            fc > 0.0_f64,
            "Fuel consumption should be positive under load: {fc}"
        );
    }
    #[test]
    fn hem_charge_mode_when_soc_low() {
        let mut hem = HybridEnergyManagement::new(100_000.0_f64, 50_000.0_f64, 8.8_f64, 0.2_f64);
        hem.soc = 0.2_f64;
        let (ice, em) = hem.power_split(30_000.0_f64);
        assert!(ice > 0.0_f64, "ICE must run in charge mode");
        assert!(
            em < 0.0_f64 || ice > 30_000.0_f64,
            "EM should be generating or ICE should over-deliver: ice={ice} em={em}"
        );
    }
    #[test]
    fn tank_com_full_level_fuel_mass() {
        let tank = FuelTankCoM::new(
            [1.5_f64, 0.0_f64, 0.3_f64],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            1.0_f64,
        );
        let expected = 60.0_f64 * 0.745_f64;
        assert!((tank.fuel_mass_kg() - expected).abs() < 1.0e-6_f64);
    }
    #[test]
    fn tank_com_empty_fill_fraction_zero() {
        let mut tank = FuelTankCoM::new(
            [0.0_f64; 3],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            0.0_f64,
        );
        tank.fill_l = 0.0_f64;
        assert_eq!(tank.fill_fraction(), 0.0_f64);
    }
    #[test]
    fn tank_com_fuel_com_lower_when_less_fuel() {
        let tank_full = FuelTankCoM::new(
            [0.0_f64, 0.0_f64, 0.3_f64],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            1.0_f64,
        );
        let tank_half = FuelTankCoM::new(
            [0.0_f64, 0.0_f64, 0.3_f64],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            0.3_f64,
        );
        assert!(
            tank_half.fuel_com()[2] < tank_full.fuel_com()[2],
            "CoM should be lower with less fuel"
        );
    }
    #[test]
    fn tank_com_consume_reduces_fill() {
        let mut tank = FuelTankCoM::new(
            [0.0_f64; 3],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            1.0_f64,
        );
        let fill0 = tank.fill_l;
        tank.consume(5.0_f64);
        assert!(tank.fill_l < fill0, "fill should decrease after consume");
    }
    #[test]
    fn tank_com_combined_com_is_finite() {
        let tank = FuelTankCoM::new(
            [1.5_f64, 0.0_f64, 0.35_f64],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            0.5_f64,
        );
        let com = tank.combined_com(3.0_f64);
        assert!(com[0].is_finite() && com[1].is_finite() && com[2].is_finite());
    }
    #[test]
    fn tank_com_slosh_shifts_lateral_com() {
        let mut tank = FuelTankCoM::new(
            [0.0_f64, 0.0_f64, 0.3_f64],
            0.2_f64,
            0.4_f64,
            0.6_f64,
            0.745_f64,
            60.0_f64,
            0.7_f64,
        );
        let y0 = tank.fuel_com()[1];
        tank.slosh_lateral_factor = 0.1_f64;
        let y1 = tank.fuel_com()[1];
        assert!(y1 != y0, "Lateral CoM should shift with sloshing");
    }
}
