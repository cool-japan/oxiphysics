//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute frictional heat power from longitudinal slip.
///
/// Q_long = mu * Fn * |omega * r - v_x|
pub fn heat_from_longitudinal_slip(
    friction: f64,
    normal_force: f64,
    wheel_speed: f64,
    vehicle_speed: f64,
) -> f64 {
    friction * normal_force * (wheel_speed - vehicle_speed).abs()
}
/// Compute frictional heat power from lateral slip.
///
/// Q_lat = mu * Fn * |v_y|
pub fn heat_from_lateral_slip(friction: f64, normal_force: f64, lateral_velocity: f64) -> f64 {
    friction * normal_force * lateral_velocity.abs()
}
/// Total frictional heat from combined longitudinal and lateral slip.
pub fn heat_from_combined_slip(
    friction: f64,
    normal_force: f64,
    long_slip_vel: f64,
    lat_slip_vel: f64,
) -> f64 {
    let v = (long_slip_vel * long_slip_vel + lat_slip_vel * lat_slip_vel).sqrt();
    friction * normal_force * v
}
/// Compute heat generation power from traction force and slip velocity.
///
/// Physical model: Q = F_traction × v_slip
///
/// This is the canonical formula for tire heating from slip.
pub fn tread_slip_heat(traction_force_n: f64, slip_velocity_m_s: f64) -> f64 {
    traction_force_n.abs() * slip_velocity_m_s.abs()
}
/// Compute slip velocity (m/s) from wheel speed, vehicle speed, and tire radius.
///
/// `v_slip = |omega * r - v|`
pub fn slip_velocity(wheel_omega_rad_s: f64, tire_radius_m: f64, vehicle_speed_m_s: f64) -> f64 {
    (wheel_omega_rad_s * tire_radius_m - vehicle_speed_m_s).abs()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ConvectiveCooling;
    use crate::InnerLinerTemperature;
    use crate::RubberDegradation;
    use crate::TemperatureGrip;
    use crate::ThermalDiffusion1D;
    use crate::TireBeltModel;
    use crate::TireCompound;
    use crate::TireContactHeat;
    use crate::TireGroundHeatExchange;
    use crate::TireHeatBalance;
    use crate::TireInflationPressure;
    use crate::TireLayerModel;
    use crate::TireLayerParams;
    use crate::TireOperatingWindow;
    use crate::TireTemperatureMap;
    use crate::TireTemperaturePredictor;
    use crate::TireThermalModel;
    use crate::TireThermalParams;
    use crate::TireWarmup;
    use crate::VulcanizationDegradation;
    fn medium_params() -> TireThermalParams {
        TireCompound::Medium.default_params()
    }
    #[test]
    fn test_grip_at_optimal_temperature() {
        let p = medium_params();
        let max_grip = p.max_grip;
        let t_opt = p.optimal_temperature;
        let mut model = TireThermalModel::new(p);
        model.state.surface_temp = t_opt;
        let grip = model.grip_coefficient();
        assert!(
            (grip - max_grip).abs() < 1e-10,
            "grip at T_opt should equal max_grip, got {grip}"
        );
    }
    #[test]
    fn test_grip_at_ambient_less_than_max() {
        let p = medium_params();
        let max_grip = p.max_grip;
        let model = TireThermalModel::new(p);
        let grip = model.grip_coefficient();
        assert!(
            grip < max_grip,
            "cold grip should be less than max_grip, got {grip} vs {max_grip}"
        );
    }
    #[test]
    fn test_step_increases_temperature_with_heat() {
        let p = medium_params();
        let initial_temp = p.ambient_temperature;
        let mut model = TireThermalModel::new(p);
        model.compute_heat_generation(5000.0, 2.0, 1.4);
        model.step(0.1);
        assert!(
            model.state.surface_temp > initial_temp,
            "surface temp should rise with heat generation"
        );
    }
    #[test]
    fn test_step_cools_toward_ambient_no_heat() {
        let p = medium_params();
        let t_amb = p.ambient_temperature;
        let mut model = TireThermalModel::new(p);
        model.state.surface_temp = t_amb + 60.0;
        model.state.core_temp = t_amb + 60.0;
        model.state.heat_generation_rate = 0.0;
        let temp_before = model.state.surface_temp;
        model.step(1.0);
        assert!(
            model.state.surface_temp < temp_before,
            "surface should cool when no heat is generated"
        );
    }
    #[test]
    fn test_soft_has_lower_t_opt_than_hard() {
        let soft = TireCompound::Soft.default_params();
        let hard = TireCompound::Hard.default_params();
        assert!(
            soft.optimal_temperature < hard.optimal_temperature,
            "Soft T_opt ({}) should be < Hard T_opt ({})",
            soft.optimal_temperature,
            hard.optimal_temperature
        );
    }
    #[test]
    fn test_temperature_map_rotate_changes_contact_sector() {
        let mut map = TireTemperatureMap::new(36, 293.15);
        let initial = map.contact_sector;
        map.rotate(10.0, 1.0);
        assert_ne!(
            map.contact_sector, initial,
            "contact sector should change after rotation"
        );
    }
    #[test]
    fn test_simulate_warmup_returns_finite_time() {
        let p = TireCompound::Soft.default_params();
        let model = TireThermalModel::new(p);
        let mut warmup = TireWarmup::new(model);
        let t = warmup.simulate_warmup(4000.0, 2.0, 0.5, 2000);
        assert!(
            t.is_finite(),
            "warmup should complete for reasonable inputs, got {t}"
        );
    }
    #[test]
    fn test_three_layer_initial_at_ambient() {
        let p = TireLayerParams::default_medium();
        let t_amb = p.ambient_temperature;
        let model = TireLayerModel::new(p);
        for &t in &model.state.temps {
            assert!(
                (t - t_amb).abs() < 1e-10,
                "all layers should start at ambient"
            );
        }
    }
    #[test]
    fn test_three_layer_tread_heats_first() {
        let p = TireLayerParams::default_medium();
        let mut model = TireLayerModel::new(p);
        model.compute_heat_generation(5000.0, 2.0, 1.3);
        model.step(0.1);
        assert!(
            model.state.temps[0] > model.state.temps[1],
            "tread should heat faster than carcass"
        );
        assert!(
            model.state.temps[1] >= model.state.temps[2],
            "carcass should be >= sidewall"
        );
    }
    #[test]
    fn test_three_layer_heat_propagates() {
        let p = TireLayerParams {
            layer_mass: [3.0, 4.0, 2.5],
            layer_cp: [1100.0, 1000.0, 950.0],
            inter_layer_conductivity: [5.0, 5.0],
            contact_area: 1.0,
            ambient_temperature: 293.15,
            tread_convection: 0.0,
            sidewall_convection: 0.0,
            sidewall_area: 0.0,
        };
        let t_amb = p.ambient_temperature;
        let mut model = TireLayerModel::new(p);
        model.compute_heat_generation(5000.0, 3.0, 1.4);
        for _ in 0..5000 {
            model.step(0.01);
        }
        assert!(
            model.state.temps[2] > t_amb + 0.01,
            "sidewall should eventually warm up, got {}",
            model.state.temps[2]
        );
    }
    #[test]
    fn test_three_layer_cooling() {
        let p = TireLayerParams::default_medium();
        let t_amb = p.ambient_temperature;
        let mut model = TireLayerModel::new(p);
        model.state.temps = [t_amb + 50.0, t_amb + 50.0, t_amb + 50.0];
        model.state.heat_generation_rate = 0.0;
        let initial_mean = model.mean_temperature();
        for _ in 0..100 {
            model.step(0.1);
        }
        assert!(
            model.mean_temperature() < initial_mean,
            "should cool without heat input"
        );
    }
    #[test]
    fn test_three_layer_gradient() {
        let p = TireLayerParams::default_medium();
        let t_amb = p.ambient_temperature;
        let mut model = TireLayerModel::new(p);
        model.state.temps = [t_amb + 60.0, t_amb + 30.0, t_amb + 10.0];
        assert!(
            (model.tread_sidewall_gradient() - 50.0).abs() < 1e-10,
            "gradient should be 50 K"
        );
    }
    #[test]
    fn test_three_layer_combined_heat_generation() {
        let p = TireLayerParams::default_medium();
        let mut model = TireLayerModel::new(p);
        model.compute_heat_generation_combined(5000.0, 1.0, 1.0, 1.0);
        let expected = 5000.0 * (2.0_f64).sqrt();
        assert!(
            (model.state.heat_generation_rate - expected).abs() < 1e-6,
            "combined heat should use total slip speed"
        );
    }
    #[test]
    fn test_temperature_grip_peak() {
        let g = TemperatureGrip {
            mu_peak: 1.5,
            mu_cold: 0.8,
            t_optimal: 370.0,
            sigma: 20.0,
        };
        let mu = g.grip_at(370.0);
        assert!((mu - 1.5).abs() < 1e-10, "at T_opt grip should be peak");
    }
    #[test]
    fn test_temperature_grip_far_from_optimal() {
        let g = TemperatureGrip {
            mu_peak: 1.5,
            mu_cold: 0.8,
            t_optimal: 370.0,
            sigma: 20.0,
        };
        let mu = g.grip_at(200.0);
        assert!(
            (mu - 0.8).abs() < 0.01,
            "far from optimal should approach cold grip, got {mu}"
        );
    }
    #[test]
    fn test_temperature_grip_derivative_zero_at_peak() {
        let g = TemperatureGrip {
            mu_peak: 1.5,
            mu_cold: 0.8,
            t_optimal: 370.0,
            sigma: 20.0,
        };
        let d = g.grip_derivative(370.0);
        assert!(d.abs() < 1e-10, "derivative at T_opt should be zero");
    }
    #[test]
    fn test_temperature_grip_in_window() {
        let g = TemperatureGrip {
            mu_peak: 1.5,
            mu_cold: 0.8,
            t_optimal: 370.0,
            sigma: 20.0,
        };
        assert!(g.in_window(375.0));
        assert!(!g.in_window(400.0));
    }
    #[test]
    fn test_convective_cooling_zero_speed() {
        let c = ConvectiveCooling {
            h_base: 20.0,
            area: 0.05,
            h_speed: 10.0,
        };
        let h = c.effective_h(0.0);
        assert!((h - 20.0).abs() < 1e-10);
    }
    #[test]
    fn test_convective_cooling_with_speed() {
        let c = ConvectiveCooling {
            h_base: 20.0,
            area: 0.05,
            h_speed: 10.0,
        };
        let h = c.effective_h(25.0);
        assert!((h - 70.0).abs() < 1e-10);
    }
    #[test]
    fn test_convective_cooling_heat_flux() {
        let c = ConvectiveCooling {
            h_base: 20.0,
            area: 0.05,
            h_speed: 0.0,
        };
        let flux = c.heat_flux(350.0, 300.0, 0.0);
        assert!((flux - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_convective_cooling_cools_surface() {
        let c = ConvectiveCooling {
            h_base: 20.0,
            area: 0.05,
            h_speed: 0.0,
        };
        let t_new = c.cool(350.0, 300.0, 0.0, 1.0, 500.0);
        assert!(t_new < 350.0, "surface should cool");
        assert!(t_new > 300.0, "should not cool below ambient");
    }
    #[test]
    fn test_diffusion_uniform_stays_uniform() {
        let mut d = ThermalDiffusion1D::new(10, 1e-4, 0.001, 350.0);
        d.step(0.001);
        for &t in &d.temps {
            assert!((t - 350.0).abs() < 1e-10, "uniform should stay uniform");
        }
    }
    #[test]
    fn test_diffusion_hot_spot_spreads() {
        let mut d = ThermalDiffusion1D::new(10, 1e-4, 0.001, 300.0);
        d.temps[5] = 400.0;
        let initial_max = d.max_temperature();
        for _ in 0..100 {
            d.step(0.001);
        }
        assert!(
            d.max_temperature() < initial_max,
            "hot spot should spread out"
        );
        assert!(
            d.temps[4] > 300.0 && d.temps[6] > 300.0,
            "neighbours should warm up"
        );
    }
    #[test]
    fn test_diffusion_add_heat() {
        let mut d = ThermalDiffusion1D::new(5, 1e-4, 0.001, 300.0);
        d.add_heat_at_node(2, 100.0, 10.0);
        assert!((d.temps[2] - 310.0).abs() < 1e-10);
    }
    #[test]
    fn test_circumferential_diffusion() {
        let mut map = TireTemperatureMap::new(36, 300.0);
        map.temps[0] = 400.0;
        let initial_var = map.temperature_variance();
        for _ in 0..200 {
            map.diffuse_circumferential(0.5, 0.01, 1.0);
        }
        assert!(
            map.temperature_variance() < initial_var,
            "variance should decrease as heat diffuses"
        );
    }
    #[test]
    fn test_temperature_map_variance_uniform() {
        let map = TireTemperatureMap::new(10, 350.0);
        assert!(
            map.temperature_variance() < 1e-10,
            "uniform should have zero variance"
        );
    }
    #[test]
    fn test_heat_from_longitudinal_slip() {
        let q = heat_from_longitudinal_slip(1.0, 5000.0, 12.0, 10.0);
        assert!((q - 10000.0).abs() < 1e-10);
    }
    #[test]
    fn test_heat_from_lateral_slip() {
        let q = heat_from_lateral_slip(1.0, 5000.0, 2.0);
        assert!((q - 10000.0).abs() < 1e-10);
    }
    #[test]
    fn test_heat_from_combined_slip() {
        let q = heat_from_combined_slip(1.0, 5000.0, 3.0, 4.0);
        assert!((q - 25000.0).abs() < 1e-10);
    }
    #[test]
    fn test_heat_zero_when_no_slip() {
        assert!(heat_from_longitudinal_slip(1.0, 5000.0, 10.0, 10.0).abs() < 1e-10);
        assert!(heat_from_lateral_slip(1.0, 5000.0, 0.0).abs() < 1e-10);
        assert!(heat_from_combined_slip(1.0, 5000.0, 0.0, 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_inflation_pressure_cold_matches_reference() {
        let m = TireInflationPressure::road_tire_cold();
        let p = m.pressure_at(m.t_cold);
        assert!(
            (p - m.p_cold).abs() < 1e-6,
            "at t_cold pressure should equal p_cold"
        );
    }
    #[test]
    fn test_inflation_pressure_rises_with_temperature() {
        let m = TireInflationPressure::road_tire_cold();
        let p_hot = m.pressure_at(373.15);
        assert!(p_hot > m.p_cold, "pressure should rise with temperature");
    }
    #[test]
    fn test_inflation_pressure_rise_positive() {
        let m = TireInflationPressure::road_tire_cold();
        let rise = m.pressure_rise(373.15);
        assert!(rise > 0.0, "pressure rise should be positive when hot");
    }
    #[test]
    fn test_inflation_pressure_over_inflation_check() {
        let m = TireInflationPressure::road_tire_cold();
        let is_over = m.is_over_inflated(403.15, 250_000.0);
        assert!(
            is_over,
            "tire should be over-inflated at 130°C above 2.5 bar limit"
        );
    }
    #[test]
    fn test_heat_balance_net_stored() {
        let mut b = TireHeatBalance::new();
        b.record_step(1000.0, 400.0, 1.0);
        assert!((b.net_heat_stored() - 600.0).abs() < 1e-9);
    }
    #[test]
    fn test_heat_balance_cooling_efficiency() {
        let mut b = TireHeatBalance::new();
        b.record_step(1000.0, 250.0, 1.0);
        assert!((b.cooling_efficiency() - 0.25).abs() < 1e-9);
    }
    #[test]
    fn test_heat_balance_peak_rate_tracked() {
        let mut b = TireHeatBalance::new();
        b.record_step(500.0, 100.0, 0.5);
        b.record_step(1500.0, 200.0, 0.5);
        b.record_step(800.0, 100.0, 0.5);
        assert!(
            (b.peak_heat_rate - 1500.0).abs() < 1e-9,
            "peak should be 1500 W"
        );
    }
    #[test]
    fn test_heat_balance_zero_consumed_gives_zero_efficiency() {
        let b = TireHeatBalance::new();
        assert_eq!(b.cooling_efficiency(), 0.0);
    }
    #[test]
    fn test_rubber_degradation_zero_at_start() {
        let d = RubberDegradation::racing_compound();
        assert!((d.degradation).abs() < 1e-12);
        assert!((d.residual_strength() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rubber_degradation_increases_above_t_ref() {
        let mut d = RubberDegradation::racing_compound();
        let t_hot = d.t_ref + 50.0;
        for _ in 0..100 {
            d.step(t_hot, 1.0);
        }
        assert!(
            d.degradation > 0.0,
            "degradation should increase above t_ref"
        );
    }
    #[test]
    fn test_rubber_degradation_clamped_at_one() {
        let mut d = RubberDegradation::racing_compound();
        d.degradation = 0.999;
        d.step(d.t_ref + 1000.0, 1000.0);
        assert!(d.degradation <= 1.0, "degradation clamped at 1.0");
    }
    #[test]
    fn test_rubber_degradation_is_thermally_degraded() {
        let mut d = RubberDegradation::racing_compound();
        d.degradation = 0.7;
        assert!(d.is_thermally_degraded(0.5), "0.7 ≥ 0.5 threshold");
        assert!(!d.is_thermally_degraded(0.8), "0.7 < 0.8 threshold");
    }
    #[test]
    fn test_ground_heat_conduction_positive_when_tread_hotter() {
        let m = TireGroundHeatExchange::new(50.0, 0.6, 300.0);
        let flux = m.conductive_flux(380.0);
        assert!(flux > 0.0, "tread hotter than road → heat flows to road");
    }
    #[test]
    fn test_ground_heat_friction_split() {
        let m = TireGroundHeatExchange::new(50.0, 0.6, 300.0);
        let tire_heat = m.tire_friction_heat(1000.0);
        let road_heat = m.road_friction_heat(1000.0);
        assert!((tire_heat - 600.0).abs() < 1e-9, "60% to tire");
        assert!((road_heat - 400.0).abs() < 1e-9, "40% to road");
        assert!(
            (tire_heat + road_heat - 1000.0).abs() < 1e-9,
            "total conserved"
        );
    }
    #[test]
    fn test_ground_heat_net_positive_when_friction_dominates() {
        let m = TireGroundHeatExchange::new(1.0, 0.9, 300.0);
        let net = m.net_tread_heat(310.0, 10000.0);
        assert!(net > 0.0, "strong friction should heat tread");
    }
    #[test]
    fn test_operating_window_in_optimal() {
        let w = TireOperatingWindow::soft_compound();
        let mid = (w.t_min + w.t_max) / 2.0;
        assert!(w.in_optimal_window(mid), "midpoint should be in window");
    }
    #[test]
    fn test_operating_window_below_minimum() {
        let w = TireOperatingWindow::soft_compound();
        assert!(
            !w.in_optimal_window(w.t_min - 5.0),
            "below t_min is out of window"
        );
    }
    #[test]
    fn test_operating_window_critical_threshold() {
        let w = TireOperatingWindow::hard_compound();
        assert!(
            w.is_critical(w.t_critical + 1.0),
            "above t_critical is critical"
        );
        assert!(
            !w.is_critical(w.t_critical - 1.0),
            "below t_critical is not critical"
        );
    }
    #[test]
    fn test_operating_window_warm_up_margin_negative_when_cold() {
        let w = TireOperatingWindow::soft_compound();
        let margin = w.warm_up_margin(w.t_min - 10.0);
        assert!(margin < 0.0, "cold tire has negative margin");
    }
    #[test]
    fn test_operating_window_normalised_position_in_range() {
        let w = TireOperatingWindow::soft_compound();
        let pos = w.normalised_position((w.t_min + w.t_max) / 2.0);
        assert!((pos - 0.5).abs() < 1e-9, "midpoint should be 0.5");
    }
    #[test]
    fn test_predictor_heating_increases_temp() {
        let t_next =
            TireTemperaturePredictor::predict_surface_temp(350.0, 2000.0, 500.0, 9000.0, 1.0);
        assert!(
            t_next > 350.0,
            "net positive heat → temperature should rise"
        );
    }
    #[test]
    fn test_predictor_cooling_decreases_temp() {
        let t_next =
            TireTemperaturePredictor::predict_surface_temp(380.0, 0.0, 1500.0, 9000.0, 1.0);
        assert!(t_next < 380.0, "no heat in, cooling → temperature drops");
    }
    #[test]
    fn test_predictor_time_to_target_zero_when_already_there() {
        let t = TireTemperaturePredictor::time_to_target(373.0, 373.0, 0.0, 0.0, 9000.0);
        assert!(t < 1e-6, "already at target → time = 0");
    }
    #[test]
    fn test_predictor_time_to_target_finite_with_net_heat() {
        let t = TireTemperaturePredictor::time_to_target(300.0, 373.0, 5000.0, 1000.0, 9000.0);
        assert!(t.is_finite() && t > 0.0, "should converge, got {t}");
        assert!((t - 164.25).abs() < 1e-6, "expected 164.25 s, got {t}");
    }
    #[test]
    fn test_predictor_time_to_target_infinite_when_unreachable() {
        let t = TireTemperaturePredictor::time_to_target(300.0, 400.0, 500.0, 1000.0, 9000.0);
        assert!(
            t.is_infinite(),
            "can't heat with net cooling → infinite time"
        );
    }
    #[test]
    fn test_tire_belt_model_initial_at_ambient() {
        let m = TireBeltModel::medium_racing();
        let t_amb = m.t_ambient;
        assert!((m.tread_temp - t_amb).abs() < 1e-10);
        assert!((m.belt_temp - t_amb).abs() < 1e-10);
        assert!((m.carcass_temp - t_amb).abs() < 1e-10);
    }
    #[test]
    fn test_tire_belt_model_tread_heats_first() {
        let mut m = TireBeltModel::medium_racing();
        for _ in 0..100 {
            m.step(20000.0, 0.01);
        }
        assert!(m.tread_temp > m.belt_temp, "tread should be hottest");
        assert!(m.belt_temp > m.carcass_temp, "belt hotter than carcass");
    }
    #[test]
    fn test_tire_belt_model_cools_without_heat() {
        let mut m = TireBeltModel::medium_racing();
        m.tread_temp = 400.0;
        m.belt_temp = 380.0;
        m.carcass_temp = 360.0;
        let mean_before = m.mean_temp();
        for _ in 0..200 {
            m.step(0.0, 0.1);
        }
        assert!(m.mean_temp() < mean_before, "should cool toward ambient");
    }
    #[test]
    fn test_tire_belt_model_gradient() {
        let mut m = TireBeltModel::medium_racing();
        for _ in 0..200 {
            m.step(15000.0, 0.01);
        }
        let grad = m.gradient();
        assert!(
            grad > 0.0,
            "tread should be hotter than carcass under heating"
        );
    }
    #[test]
    fn test_vulcanization_starts_undercured() {
        let v = VulcanizationDegradation::racing_slick();
        assert_eq!(v.cure_state, 0.0, "starts uncured");
        assert_eq!(v.degradation_factor(), 0.0, "no degradation at zero cure");
    }
    #[test]
    fn test_vulcanization_grip_under_cure_is_low() {
        let v = VulcanizationDegradation::racing_slick();
        let grip = v.grip_multiplier();
        assert_eq!(grip, 0.0, "uncured rubber has zero grip");
    }
    #[test]
    fn test_vulcanization_optimal_gives_full_grip() {
        let mut v = VulcanizationDegradation::racing_slick();
        v.cure_state = 1.0;
        let grip = v.grip_multiplier();
        assert!((grip - 1.0).abs() < 1e-9, "at optimal cure: full grip");
    }
    #[test]
    fn test_vulcanization_over_cure_degrades_grip() {
        let mut v = VulcanizationDegradation::racing_slick();
        v.cure_state = 5.0;
        let grip = v.grip_multiplier();
        assert!(grip < 1.0, "over-cured rubber should lose grip, got {grip}");
    }
    #[test]
    fn test_vulcanization_rate_increases_with_temperature() {
        let v = VulcanizationDegradation::racing_slick();
        let rate_cool = v.cure_rate(350.0);
        let rate_hot = v.cure_rate(450.0);
        assert!(
            rate_hot > rate_cool,
            "higher temperature → faster cure rate"
        );
    }
    #[test]
    fn test_tread_slip_heat_proportional() {
        let q = tread_slip_heat(1000.0, 2.0);
        assert!((q - 2000.0).abs() < 1e-9, "Q = F * v = 2000 W, got {q}");
    }
    #[test]
    fn test_tread_slip_heat_zero_when_no_slip() {
        let q = tread_slip_heat(1000.0, 0.0);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_tread_slip_heat_absolute_value() {
        let q = tread_slip_heat(-500.0, -3.0);
        assert!((q - 1500.0).abs() < 1e-9, "heat must be positive, got {q}");
    }
    #[test]
    fn test_slip_velocity_zero_when_no_slip() {
        let v_slip = slip_velocity(10.0, 0.3, 3.0);
        assert!(v_slip.abs() < 1e-9, "no slip → zero slip velocity");
    }
    #[test]
    fn test_slip_velocity_positive_when_spinning() {
        let v_slip = slip_velocity(15.0, 0.3, 3.0);
        assert!((v_slip - 1.5).abs() < 1e-9, "slip = 1.5 m/s, got {v_slip}");
    }
    #[test]
    fn test_inner_liner_starts_at_ambient() {
        let liner = InnerLinerTemperature::default_tpms();
        assert!((liner.temp - liner.t_ambient).abs() < 1e-9);
    }
    #[test]
    fn test_inner_liner_tracks_carcass() {
        let mut liner = InnerLinerTemperature::default_tpms();
        for _ in 0..1000 {
            liner.update(400.0, 0.1);
        }
        assert!(
            liner.temp > 350.0,
            "liner should track hot carcass, got {}",
            liner.temp
        );
    }
    #[test]
    fn test_inner_liner_pressure_rises_with_temp() {
        let mut liner = InnerLinerTemperature::default_tpms();
        liner.temp = 373.15;
        let p = liner.tpms_pressure(220_000.0, 293.15);
        assert!(
            p > 220_000.0,
            "hot liner → higher pressure than cold, got {p}"
        );
    }
    #[test]
    fn test_inner_liner_pressure_warning() {
        let mut liner = InnerLinerTemperature::default_tpms();
        liner.temp = 500.0;
        let warn = liner.pressure_warning(220_000.0, 293.15, 0.0, 300_000.0);
        assert!(warn, "over-heated liner should trigger pressure warning");
    }
    #[test]
    fn test_tire_contact_heat_split_sums_to_total() {
        let m = TireContactHeat::dry_asphalt();
        let total = 8000.0;
        let tire_heat = m.tire_slip_heat(total);
        let road_heat = m.road_slip_heat(total);
        assert!(
            (tire_heat + road_heat - total).abs() < 1e-9,
            "heat split must sum to total"
        );
    }
    #[test]
    fn test_tire_contact_heat_tire_gets_majority() {
        let m = TireContactHeat::dry_asphalt();
        let tire_heat = m.tire_slip_heat(1000.0);
        let road_heat = m.road_slip_heat(1000.0);
        assert!(
            tire_heat > road_heat,
            "tire should get more heat than road (fraction=0.6)"
        );
    }
    #[test]
    fn test_tire_contact_heat_net_positive_with_high_friction() {
        let m = TireContactHeat::dry_asphalt();
        let net = m.net_tread_heat(320.0, 50000.0);
        assert!(net > 0.0, "high friction → net positive heat to tread");
    }
    #[test]
    fn test_tire_contact_heat_net_negative_cool_road() {
        let m = TireContactHeat::new(0.6, 280.0, 200.0);
        let net = m.net_tread_heat(400.0, 0.0);
        assert!(
            net < 0.0,
            "cool road and no friction → net heat negative (tread cools)"
        );
    }
}
#[cfg(test)]
mod tests_tire_thermal_new {

    use crate::TireCompound;
    use crate::TireLayerModel;
    use crate::TireLayerParams;
    use crate::TireThermalModel;

    fn default_medium_model() -> TireThermalModel {
        TireThermalModel::new(TireCompound::Medium.default_params())
    }
    fn default_layer_model() -> TireLayerModel {
        TireLayerModel::new(TireLayerParams::default_medium())
    }
    #[test]
    fn test_cornering_heat_zero_when_no_slip() {
        let mut model = default_medium_model();
        model.compute_heat_generation_cornering(3000.0, 0.0, 0.9);
        assert_eq!(
            model.state.heat_generation_rate, 0.0,
            "no lateral slip → no heat"
        );
    }
    #[test]
    fn test_cornering_heat_positive_with_slip() {
        let mut model = default_medium_model();
        model.compute_heat_generation_cornering(4000.0, 2.0, 0.9);
        assert!(
            model.state.heat_generation_rate > 0.0,
            "lateral force + slip → positive heat: {}",
            model.state.heat_generation_rate
        );
    }
    #[test]
    fn test_cornering_heat_scales_with_force() {
        let mut m1 = default_medium_model();
        let mut m2 = default_medium_model();
        m1.compute_heat_generation_cornering(2000.0, 1.5, 0.9);
        m2.compute_heat_generation_cornering(4000.0, 1.5, 0.9);
        assert!(
            m2.state.heat_generation_rate > m1.state.heat_generation_rate,
            "higher force → more heat"
        );
    }
    #[test]
    fn test_peak_grip_temperature_equals_optimal() {
        let model = default_medium_model();
        let (t_opt, _) = model.compute_peak_grip_temperature();
        assert!(
            (t_opt - model.params.optimal_temperature).abs() < 1e-9,
            "peak grip T should match optimal_temperature: {t_opt}"
        );
    }
    #[test]
    fn test_peak_grip_window_positive() {
        let model = default_medium_model();
        let (_, window) = model.compute_peak_grip_temperature();
        assert!(window > 0.0, "window should be positive: {window}");
    }
    #[test]
    fn test_blister_risk_zero_at_cold_tire() {
        let model = default_medium_model();
        let risk = model.compute_blister_risk(model.params.optimal_temperature + 40.0, 8000.0);
        assert_eq!(risk, 0.0, "cold tire → blister risk = 0");
    }
    #[test]
    fn test_blister_risk_increases_above_limit() {
        let mut model = default_medium_model();
        model.state.surface_temp = model.params.optimal_temperature + 100.0;
        model.state.heat_generation_rate = 9000.0;
        let risk = model.compute_blister_risk(model.params.optimal_temperature + 40.0, 8000.0);
        assert!(
            risk > 0.0,
            "hot tire with high flux → blister risk > 0: {risk}"
        );
    }
    #[test]
    fn test_blister_risk_bounded_to_one() {
        let mut model = default_medium_model();
        model.state.surface_temp = 900.0;
        model.state.heat_generation_rate = 1_000_000.0;
        let risk = model.compute_blister_risk(350.0, 1.0);
        assert!(risk <= 1.0, "blister risk must not exceed 1.0: {risk}");
    }
    #[test]
    fn test_layer_model_cornering_heat_positive() {
        let mut model = default_layer_model();
        model.compute_heat_generation_cornering(5000.0, 1.5, 0.85);
        assert!(
            model.state.heat_generation_rate > 0.0,
            "layer model: lateral heat should be positive"
        );
    }
    #[test]
    fn test_layer_model_blister_risk_zero_cold() {
        let model = default_layer_model();
        let t_limit = model.params.ambient_temperature + 80.0 + 40.0;
        let risk = model.compute_blister_risk(t_limit, 8000.0);
        assert_eq!(risk, 0.0, "cold layer model → blister risk = 0");
    }
    #[test]
    fn test_layer_model_peak_grip_temperature_above_ambient() {
        let model = default_layer_model();
        let (t_opt, _) = model.compute_peak_grip_temperature();
        assert!(
            t_opt > model.params.ambient_temperature,
            "peak grip T should be above ambient: {t_opt} vs {}",
            model.params.ambient_temperature
        );
    }
}
