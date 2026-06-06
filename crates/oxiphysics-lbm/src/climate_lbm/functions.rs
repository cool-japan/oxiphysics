//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Stefan-Boltzmann constant (W m⁻² K⁻⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_4e-8;
/// Solar constant (W m⁻²).
pub const SOLAR_CONSTANT: f64 = 1361.0;
/// Pre-industrial atmospheric CO₂ (ppm).
pub const CO2_PREINDUSTRIAL: f64 = 280.0;
/// Current atmospheric CO₂ baseline used in examples (ppm).
pub const CO2_CURRENT: f64 = 420.0;
/// Earth mean radius (m).
pub const EARTH_RADIUS: f64 = 6.371e6;
/// Specific heat capacity of sea water (J kg⁻¹ K⁻¹).
pub const CP_SEAWATER: f64 = 3993.0;
/// Density of sea water (kg m⁻³).
pub const RHO_SEAWATER: f64 = 1025.0;
/// Latent heat of ice fusion (J kg⁻¹).
pub const LATENT_HEAT_ICE: f64 = 3.34e5;
/// D2Q9 velocity set weights.
pub(super) const W9: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D2Q9 velocity vectors (ex, ey).
pub(super) const E9X: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
pub(super) const E9Y: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// Convert temperature from Celsius to Kelvin.
pub fn celsius_to_kelvin(t_c: f64) -> f64 {
    t_c + 273.15
}
/// Convert temperature from Kelvin to Celsius.
pub fn kelvin_to_celsius(t_k: f64) -> f64 {
    t_k - 273.15
}
/// Blackbody emission (W m⁻²) at temperature T (K).
pub fn blackbody_emission(t_k: f64) -> f64 {
    STEFAN_BOLTZMANN * t_k.powi(4)
}
/// Effective radiating temperature for a given OLR (K).
pub fn effective_radiating_temperature(olr: f64) -> f64 {
    (olr / STEFAN_BOLTZMANN).powf(0.25)
}
/// Atmospheric transmission window fraction (0–1).
///
/// Simplified model: window = 0.4 − 0.02 · ln(CO₂ / 280).
pub fn transmission_window(co2_ppm: f64) -> f64 {
    (0.4 - 0.02 * (co2_ppm / CO2_PREINDUSTRIAL).ln()).clamp(0.0, 1.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::climate_lbm::types::*;
    #[test]
    fn test_ebm_absorbed_shortwave() {
        let ebm = EnergyBalanceModel::new();
        let asw = ebm.absorbed_shortwave();
        assert!((asw - 1361.0 / 4.0 * 0.70).abs() < 0.1);
    }
    #[test]
    fn test_ebm_outgoing_longwave() {
        let ebm = EnergyBalanceModel::new();
        let olr = ebm.outgoing_longwave();
        assert!((olr - (203.3 + 2.09 * 14.0)).abs() < 1e-10);
    }
    #[test]
    fn test_ebm_net_toa_flux() {
        let ebm = EnergyBalanceModel::new();
        let flux = ebm.net_toa_flux();
        assert!(flux.is_finite());
    }
    #[test]
    fn test_ebm_equilibrium_temperature() {
        let ebm = EnergyBalanceModel::new();
        let t_eq = ebm.equilibrium_temperature(0.30, 0.0);
        assert!((t_eq - 14.0).abs() < 5.0);
    }
    #[test]
    fn test_ebm_step_converges() {
        let mut ebm = EnergyBalanceModel::new();
        ebm.temperature = 0.0;
        for _ in 0..1000 {
            ebm.step(0.1, 0.0);
        }
        let t_eq = ebm.equilibrium_temperature(0.30, 0.0);
        assert!((ebm.temperature - t_eq).abs() < 0.5);
    }
    #[test]
    fn test_ebm_run_to_equilibrium() {
        let mut ebm = EnergyBalanceModel::new();
        ebm.temperature = 0.0;
        let steps = ebm.run_to_equilibrium(0.05, 1e-4, 50000);
        assert!(steps < 50000);
    }
    #[test]
    fn test_co2_radiative_forcing_zero_at_ref() {
        let model = Co2FeedbackModel::new();
        assert!(model.radiative_forcing().abs() < 1e-10);
    }
    #[test]
    fn test_co2_radiative_forcing_positive_above_ref() {
        let mut model = Co2FeedbackModel::new();
        model.co2_ppm = 560.0;
        assert!(model.radiative_forcing() > 0.0);
    }
    #[test]
    fn test_co2_doubling_forcing() {
        let f = Co2FeedbackModel::co2_doubling_forcing();
        assert!((f - 3.71).abs() < 0.1);
    }
    #[test]
    fn test_co2_step_increases_ppm() {
        let mut model = Co2FeedbackModel::new();
        let ppm0 = model.co2_ppm;
        model.step(1.0);
        assert!(model.co2_ppm > ppm0);
    }
    #[test]
    fn test_co2_doubling_year_none_initially() {
        let model = Co2FeedbackModel::new();
        assert!(model.doubling_year().is_none());
    }
    #[test]
    fn test_ice_fraction_warm() {
        let mut ice = IceAlbedoFeedback::new();
        ice.temperature = 20.0;
        assert!(ice.ice_fraction() < 0.01);
    }
    #[test]
    fn test_ice_fraction_cold() {
        let mut ice = IceAlbedoFeedback::new();
        ice.temperature = -20.0;
        assert!(ice.ice_fraction() > 0.99);
    }
    #[test]
    fn test_effective_albedo_range() {
        let ice = IceAlbedoFeedback::new();
        let a = ice.effective_albedo();
        assert!((0.0..=1.0).contains(&a));
    }
    #[test]
    fn test_ice_albedo_feedback_positive() {
        let ice = IceAlbedoFeedback::new();
        assert!(ice.feedback_parameter() > 0.0);
    }
    #[test]
    fn test_ocean_heat_capacity() {
        let ocean = OceanHeatUptake::new();
        let c = ocean.heat_capacity();
        assert!((c - RHO_SEAWATER * CP_SEAWATER * 50.0).abs() < 1.0);
    }
    #[test]
    fn test_ocean_air_sea_flux_equilibrium() {
        let ocean = OceanHeatUptake::new();
        assert!(ocean.air_sea_flux().abs() < 1e-10);
    }
    #[test]
    fn test_ocean_step_warms_when_forced() {
        let mut ocean = OceanHeatUptake::new();
        ocean.atm_temperature = 20.0;
        let t0 = ocean.temperature;
        ocean.step(86400.0);
        assert!(ocean.temperature > t0);
    }
    #[test]
    fn test_radiative_forcing_total() {
        let mut rf = RadiativeForcing::zero();
        rf.co2 = 2.0;
        rf.aerosol_direct = -0.5;
        assert!((rf.total() - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_radiative_forcing_from_co2_ratio() {
        let f = RadiativeForcing::from_co2_ratio(2.0);
        assert!((f - 5.35 * 2_f64.ln()).abs() < 1e-10);
    }
    #[test]
    fn test_radiative_forcing_anthropogenic() {
        let mut rf = RadiativeForcing::zero();
        rf.co2 = 3.0;
        rf.solar = 0.1;
        assert!((rf.anthropogenic() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_ghg_co2_forcing_zero_at_ref() {
        let ghg = GreenhouseGasEffect::new();
        assert!(ghg.co2_forcing().abs() < 1e-10);
    }
    #[test]
    fn test_ghg_total_forcing_positive_with_elevated_co2() {
        let mut ghg = GreenhouseGasEffect::new();
        ghg.co2_ppm = 560.0;
        assert!(ghg.total_forcing() > 0.0);
    }
    #[test]
    fn test_ghg_effective_emissivity_range() {
        let ghg = GreenhouseGasEffect::new();
        let eps = ghg.effective_emissivity();
        assert!((0.0..=1.0).contains(&eps));
    }
    #[test]
    fn test_ghg_greenhouse_enhancement_positive() {
        let ghg = GreenhouseGasEffect::new();
        assert!(ghg.greenhouse_enhancement() > 0.0);
    }
    #[test]
    fn test_climate_sensitivity_ecs_positive() {
        let cs = ClimateSensitivity::new();
        assert!(cs.ecs() > 0.0);
    }
    #[test]
    fn test_climate_sensitivity_ecs_range() {
        let cs = ClimateSensitivity::new();
        let ecs = cs.ecs();
        assert!(ecs > 1.5 && ecs < 6.0);
    }
    #[test]
    fn test_tcr_less_than_ecs() {
        let cs = ClimateSensitivity::new();
        assert!(cs.tcr() < cs.ecs());
    }
    #[test]
    fn test_delta_t_eq_proportional_to_forcing() {
        let cs = ClimateSensitivity::new();
        let dt1 = cs.delta_t_eq(1.0);
        let dt2 = cs.delta_t_eq(2.0);
        assert!((dt2 / dt1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_sea_level_rate_positive_for_warming() {
        let sl = SeaLevelRiseModel::new();
        assert!(sl.rate(1.0) > 0.0);
    }
    #[test]
    fn test_sea_level_step_increases() {
        let mut sl = SeaLevelRiseModel::new();
        sl.step(1.0, 2.0);
        assert!(sl.sea_level > 0.0);
    }
    #[test]
    fn test_sea_level_ice_fraction_range() {
        let sl = SeaLevelRiseModel::new();
        let f = sl.ice_fraction();
        assert!((0.0..=1.0).contains(&f));
    }
    #[test]
    fn test_carbon_total_conserved_approx() {
        let mut cc = CarbonCycle::new();
        let total0 = cc.total_carbon();
        cc.step(1.0);
        let total1 = cc.total_carbon();
        assert!((total1 - total0 - cc.emission_rate).abs() < 5.0);
    }
    #[test]
    fn test_carbon_co2_ppm() {
        let cc = CarbonCycle::new();
        assert!((cc.co2_ppm() - 590.0 / 2.12).abs() < 0.1);
    }
    #[test]
    fn test_carbon_npp_increases_with_co2() {
        let mut cc = CarbonCycle::new();
        let npp0 = cc.npp();
        cc.atmosphere *= 2.0;
        assert!(cc.npp() > npp0);
    }
    #[test]
    fn test_lbm_grid_initial_temperature() {
        let grid = ClimateLbmGrid::new(8, 4, 15.0);
        assert!((grid.global_mean_temperature() - 15.0).abs() < 1.0);
    }
    #[test]
    fn test_lbm_grid_stream_conserves_total() {
        let mut grid = ClimateLbmGrid::new(8, 4, 15.0);
        let t0: f64 = grid
            .f
            .iter()
            .flat_map(|r| r.iter())
            .flat_map(|c| c.iter())
            .sum();
        grid.stream();
        let t1: f64 = grid
            .f
            .iter()
            .flat_map(|r| r.iter())
            .flat_map(|c| c.iter())
            .sum();
        assert!((t0 - t1).abs() < 1e-8);
    }
    #[test]
    fn test_lbm_grid_collide_no_nan() {
        let mut grid = ClimateLbmGrid::new(4, 4, 14.0);
        grid.collide();
        for row in &grid.f {
            for cell in row {
                for &v in cell.iter() {
                    assert!(v.is_finite());
                }
            }
        }
    }
    #[test]
    fn test_lbm_grid_step_finite() {
        let mut grid = ClimateLbmGrid::new(4, 4, 14.0);
        grid.step(0.01, 3.7);
        assert!(grid.global_mean_temperature().is_finite());
    }
    #[test]
    fn test_climate_simulation_snapshot() {
        let mut sim = ClimateSimulation::new(4, 4, 1850.0);
        sim.step(1.0);
        let snap = sim.snapshot();
        assert!(snap.co2_ppm > CO2_PREINDUSTRIAL);
        assert!(snap.year > 1850.0);
    }
    #[test]
    fn test_climate_simulation_sea_level_rises() {
        let mut sim = ClimateSimulation::new(4, 4, 1850.0);
        let delta_t = 4.0;
        sim.sea_level.step(10.0, delta_t);
        assert!(sim.sea_level.sea_level > 0.0);
    }
    #[test]
    fn test_celsius_kelvin_roundtrip() {
        let t_c = 25.0;
        let t_k = celsius_to_kelvin(t_c);
        assert!((kelvin_to_celsius(t_k) - t_c).abs() < 1e-10);
    }
    #[test]
    fn test_blackbody_emission_sun() {
        let sigma_t4 = blackbody_emission(5778.0);
        assert!((sigma_t4 - 5.67e-8 * 5778_f64.powi(4)).abs() / sigma_t4 < 1e-4);
    }
    #[test]
    fn test_effective_radiating_temperature_roundtrip() {
        let t = 255.0;
        let olr = blackbody_emission(t);
        let t_back = effective_radiating_temperature(olr);
        assert!((t_back - t).abs() < 1e-6);
    }
    #[test]
    fn test_transmission_window_decreases_with_co2() {
        let w0 = transmission_window(280.0);
        let w1 = transmission_window(560.0);
        assert!(w1 < w0);
    }
    #[test]
    fn test_enso_step_finite() {
        let mut enso = EnsoProxy::new();
        enso.step(1.0);
        assert!(enso.nino34().is_finite());
    }
    #[test]
    fn test_enso_phase_name() {
        let mut enso = EnsoProxy::new();
        enso.sst_east = 1.0;
        assert_eq!(enso.phase_name(), "El Nino");
        enso.sst_east = -1.0;
        assert_eq!(enso.phase_name(), "La Nina");
        enso.sst_east = 0.0;
        assert_eq!(enso.phase_name(), "Neutral");
    }
    #[test]
    fn test_milankovitch_insolation_positive() {
        let m = MilankovitchForcing::present_day();
        let ins = m.insolation_65n_summer();
        assert!(ins > 0.0);
    }
    #[test]
    fn test_milankovitch_insolation_toa_finite() {
        let m = MilankovitchForcing::present_day();
        let ins = m.insolation_toa(0.0, 172.0);
        assert!(ins.is_finite() && ins >= 0.0);
    }
    #[test]
    fn test_permafrost_stefan_depth_nonnegative() {
        let pf = PermafrostModel::new();
        assert!(pf.stefan_depth() >= 0.0);
    }
    #[test]
    fn test_permafrost_carbon_release_increases_with_temperature() {
        let mut pf = PermafrostModel::new();
        let cr0 = pf.carbon_release_rate();
        pf.update_temperature(0.0);
        let cr1 = pf.carbon_release_rate();
        assert!(cr1 >= cr0);
    }
    #[test]
    fn test_zonal_mean_profile_from_grid() {
        let grid = ClimateLbmGrid::new(8, 4, 14.0);
        let profile = ZonalMeanProfile::from_grid(&grid);
        assert_eq!(profile.latitudes.len(), 4);
        for &t in &profile.temperature {
            assert!(t.is_finite());
        }
    }
    #[test]
    fn test_hadley_cell_width_positive() {
        let hc = HadleyCell::new();
        assert!(hc.cell_width_radians() > 0.0);
    }
    #[test]
    fn test_hadley_upwelling_positive() {
        let hc = HadleyCell::new();
        assert!(hc.upwelling_velocity() > 0.0);
    }
    #[test]
    fn test_monsoon_active() {
        let m = MonsoonIndex::new();
        assert!(m.is_monsoon_active());
    }
    #[test]
    fn test_monsoon_precipitation_positive_when_active() {
        let m = MonsoonIndex::new();
        assert!(m.precipitation_proxy() > 0.0);
    }
    #[test]
    fn test_cloud_feedback_forcing_finite() {
        let cf = CloudFeedback::new();
        let f = cf.net_cloud_forcing();
        assert!(f.is_finite());
    }
    #[test]
    fn test_low_cloud_decreases_with_warming() {
        let mut cf = CloudFeedback::new();
        let lc0 = cf.low_cloud_fraction();
        cf.temperature = 20.0;
        let lc1 = cf.low_cloud_fraction();
        assert!(lc1 < lc0);
    }
    #[test]
    fn test_arctic_amplification_factor_finite() {
        let grid_ref = ClimateLbmGrid::new(8, 12, 14.0);
        let mut grid_now = ClimateLbmGrid::new(8, 12, 14.0);
        for j in 10..12 {
            for i in 0..8 {
                grid_now.temp[j][i] += 4.0;
            }
        }
        for j in 0..10 {
            for i in 0..8 {
                grid_now.temp[j][i] += 1.0;
            }
        }
        let p_ref = ZonalMeanProfile::from_grid(&grid_ref);
        let p_now = ZonalMeanProfile::from_grid(&grid_now);
        let aa = ArcticAmplification::compute(&p_now, &p_ref);
        assert!(aa.amplification_factor.is_finite());
        assert!(aa.is_strong());
    }
}
