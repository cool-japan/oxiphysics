//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::aerospace::AblativeHeatShield;

    use crate::aerospace::CfrpLaminate;
    use crate::aerospace::CfrpPly;
    use crate::aerospace::CmcMaterial;

    use crate::aerospace::CrashAbsorber;
    use crate::aerospace::DamageToleranceFlaw;
    use crate::aerospace::FatigueSpectrum;

    use crate::aerospace::HoneycombSandwich;

    use crate::aerospace::SuperalloyMaterial;
    use crate::aerospace::ThermalProtection;
    use crate::aerospace::TitaniumAlloy;

    #[test]
    fn cfrp_ply_transformed_stiffness_0deg() {
        let ply = CfrpPly::im7_8552(0.0, 0.125);
        let q = ply.transformed_stiffness();
        assert!(q[0] > 100.0, "Q11bar at 0° should be large");
    }
    #[test]
    fn cfrp_laminate_thickness() {
        let lam = CfrpLaminate::quasi_isotropic(0.125);
        assert!((lam.thickness - 8.0 * 0.125).abs() < 1e-10);
    }
    #[test]
    fn cfrp_a_matrix_positive_diagonal() {
        let lam = CfrpLaminate::quasi_isotropic(0.125);
        let [a11, a22, _a12, a66] = lam.a_matrix();
        assert!(a11 > 0.0);
        assert!(a22 > 0.0);
        assert!(a66 > 0.0);
    }
    #[test]
    fn cfrp_effective_ex_positive() {
        let lam = CfrpLaminate::quasi_isotropic(0.125);
        assert!(lam.effective_ex() > 0.0);
    }
    #[test]
    fn cfrp_first_ply_failure_no_load() {
        let lam = CfrpLaminate::quasi_isotropic(0.125);
        let rf = lam.first_ply_failure_index(0.0, 0.0, 0.0);
        assert!(rf.is_infinite() || rf > 1.0);
    }
    #[test]
    fn cmc_sic_sic_designation() {
        let m = CmcMaterial::sic_sic();
        assert_eq!(m.designation, "SiC/SiC");
    }
    #[test]
    fn cmc_critical_crack_length_positive() {
        let m = CmcMaterial::sic_sic();
        let ac = m.critical_crack_length(200.0, 1.0);
        assert!(ac > 0.0);
    }
    #[test]
    fn cmc_thermal_shock_resistance_positive() {
        let m = CmcMaterial::sic_sic();
        let r = m.thermal_shock_resistance(0.25);
        assert!(r > 0.0);
    }
    #[test]
    fn cmc_c_sic_higher_max_temp() {
        let sic = CmcMaterial::sic_sic();
        let csic = CmcMaterial::c_sic();
        assert!(csic.max_temperature > sic.max_temperature);
    }
    #[test]
    fn superalloy_in718_properties() {
        let m = SuperalloyMaterial::in718();
        assert_eq!(m.designation, "IN718");
        assert!(m.sigma_02_650 > 0.0);
    }
    #[test]
    fn superalloy_oxidation_mass_gain_sqrt_law() {
        let m = SuperalloyMaterial::in718();
        let mg1 = m.oxidation_mass_gain(100.0);
        let mg4 = m.oxidation_mass_gain(400.0);
        assert!((mg4 / mg1 - 2.0).abs() < 0.01);
    }
    #[test]
    fn superalloy_creep_estimate_positive() {
        let m = SuperalloyMaterial::in718();
        let life = m.creep_rupture_estimate(923.15, 500.0);
        assert!(life > 0.0);
    }
    #[test]
    fn ti64_crack_growth_below_threshold() {
        let ti = TitaniumAlloy::ti64();
        let rate = ti.crack_growth_rate(1.0);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn ti64_crack_growth_above_threshold() {
        let ti = TitaniumAlloy::ti64();
        let rate = ti.crack_growth_rate(10.0);
        assert!(rate > 0.0);
    }
    #[test]
    fn ti64_fatigue_life_positive() {
        let ti = TitaniumAlloy::ti64();
        let n = ti.fatigue_life_cycles(0.001, 0.01, 100.0);
        assert!(n > 0.0);
    }
    #[test]
    fn ti64_scc_critical_stress() {
        let ti = TitaniumAlloy::ti64();
        let sigma = ti.scc_critical_stress(0.001, 60.0);
        assert!(sigma > 0.0);
    }
    #[test]
    fn avcoat_recession_rate_positive() {
        let a = AblativeHeatShield::avcoat();
        let rate = a.recession_rate(1e6);
        assert!(rate > 0.0);
    }
    #[test]
    fn avcoat_mass_loss_increases_with_time() {
        let a = AblativeHeatShield::avcoat();
        let m1 = a.mass_loss(1e6, 1.0);
        let m2 = a.mass_loss(1e6, 2.0);
        assert!(m2 > m1);
    }
    #[test]
    fn pica_conductivity_increases_with_temp() {
        let p = ThermalProtection::pica();
        let k_low = p.conductivity_at(100.0);
        let k_high = p.conductivity_at(1000.0);
        assert!(k_high > k_low);
    }
    #[test]
    fn pica_radiation_flux_positive() {
        let p = ThermalProtection::pica();
        let q = p.radiation_heat_flux(1500.0);
        assert!(q > 0.0);
    }
    #[test]
    fn pica_required_thickness_positive() {
        let p = ThermalProtection::pica();
        let t = p.required_thickness(100_000.0, 200.0);
        assert!(t > 0.0);
    }
    #[test]
    fn honeycomb_total_thickness() {
        let h = HoneycombSandwich::aluminum_honeycomb(6.0, 50.0, 2.0, 20.0);
        assert!((h.total_thickness() - 24.0).abs() < 0.01);
    }
    #[test]
    fn honeycomb_flexural_rigidity_positive() {
        let h = HoneycombSandwich::aluminum_honeycomb(6.0, 50.0, 2.0, 20.0);
        assert!(h.flexural_rigidity() > 0.0);
    }
    #[test]
    fn honeycomb_wrinkling_stress_positive() {
        let h = HoneycombSandwich::aluminum_honeycomb(6.0, 50.0, 2.0, 20.0);
        assert!(h.wrinkling_stress() > 0.0);
    }
    #[test]
    fn crash_absorber_peak_force() {
        let c = CrashAbsorber::aluminum_foam(2.0, 10000.0, 200.0);
        assert!((c.peak_force() - 2.0 * 10000.0 * 0.8).abs() < 0.01);
    }
    #[test]
    fn crash_absorber_sea_positive() {
        let c = CrashAbsorber::aluminum_foam(2.0, 10000.0, 200.0);
        assert!(c.specific_energy_absorption() > 0.0);
    }
    #[test]
    fn crash_absorber_stroke_for_energy() {
        let c = CrashAbsorber::aluminum_foam(2.0, 10000.0, 200.0);
        let stroke = c.stroke_for_energy(1000.0);
        assert!(stroke > 0.0);
    }
    #[test]
    fn fatigue_spectrum_total_flights() {
        let s = FatigueSpectrum::medium_transport();
        let expected = 120_000.0 * 0.5;
        assert!((s.total_flights() - expected).abs() < 0.01);
    }
    #[test]
    fn fatigue_spectrum_exceedances_at_1g() {
        let s = FatigueSpectrum::medium_transport();
        let e = s.exceedances_at(1.0);
        assert!(e > 0.0);
    }
    #[test]
    fn fatigue_spectrum_miners_rule() {
        let s = FatigueSpectrum::medium_transport();
        let pairs = vec![(1000.0, 2000.0), (500.0, 10000.0)];
        let damage = s.miners_rule_damage(&pairs);
        assert!((damage - 0.55).abs() < 0.01);
    }
    #[test]
    fn damage_tolerance_requires_repair() {
        let flaw = DamageToleranceFlaw::new(0.001, 0.02, 5000.0);
        assert!(!flaw.requires_repair(0.005));
        assert!(flaw.requires_repair(0.015));
    }
    #[test]
    fn damage_tolerance_is_critical() {
        let flaw = DamageToleranceFlaw::new(0.001, 0.02, 5000.0);
        assert!(!flaw.is_critical(0.01));
        assert!(flaw.is_critical(0.025));
    }
    #[test]
    fn damage_tolerance_life_factor() {
        let flaw = DamageToleranceFlaw::new(0.001, 0.02, 5000.0);
        assert!((flaw.life_factor() - 20.0).abs() < 0.01);
    }
    #[test]
    fn damage_tolerance_inspections() {
        let flaw = DamageToleranceFlaw::new(0.001, 0.02, 5000.0);
        let n = flaw.inspections_before_critical(1e-6);
        assert!(n > 0.0);
    }
}
#[cfg(test)]
mod aerospace_extended_tests {

    use crate::aerospace::AblationCharModel;

    use crate::aerospace::AdhesiveJoint;
    use crate::aerospace::AerogelInsulation;
    use crate::aerospace::CarbonCarbonComposite;
    use crate::aerospace::CeramicMatrixCompositeExtended;

    use crate::aerospace::CompositeLaminateClt;

    use crate::aerospace::FoamCoreSandwich;

    use crate::aerospace::NickelSuperalloyExtended;
    use crate::aerospace::SpaceEnvironmentEffects;

    use crate::aerospace::TitaniumAlloyExtended;
    #[test]
    fn ablation_char_fraction_zero_below_onset() {
        let m = AblationCharModel::pica_char();
        assert_eq!(m.char_fraction(300.0), 0.0);
    }
    #[test]
    fn ablation_char_fraction_one_above_end() {
        let m = AblationCharModel::pica_char();
        assert_eq!(m.char_fraction(800.0), 1.0);
    }
    #[test]
    fn ablation_effective_density_between_bounds() {
        let m = AblationCharModel::pica_char();
        let rho = m.effective_density(0.5);
        assert!(rho > m.rho_char && rho < m.rho_virgin);
    }
    #[test]
    fn ablation_recession_velocity_positive() {
        let m = AblationCharModel::carbon_phenolic();
        let v = m.recession_velocity(1.0e7);
        assert!(v > 0.0);
    }
    #[test]
    fn ablation_pyrolysis_gas_flux_positive() {
        let m = AblationCharModel::pica_char();
        let flux = m.pyrolysis_gas_flux(1.0e6);
        assert!(flux > 0.0);
    }
    #[test]
    fn ablation_backface_cools_with_char() {
        let m = AblationCharModel::pica_char();
        let t_back = m.backface_temperature(1.0e5, 0.01, 1600.0);
        assert!(t_back < 1600.0);
    }
    #[test]
    fn clt_cross_ply_thickness() {
        let lam = CompositeLaminateClt::cross_ply(4, 0.125, 161.0, 11.4, 5.17, 0.32);
        assert!((lam.total_thickness_mm() - 8.0 * 0.125).abs() < 1.0e-10);
    }
    #[test]
    fn clt_a11_positive() {
        let lam = CompositeLaminateClt::cross_ply(2, 0.125, 161.0, 11.4, 5.17, 0.32);
        assert!(lam.a11() > 0.0);
    }
    #[test]
    fn clt_a22_positive() {
        let lam = CompositeLaminateClt::cross_ply(2, 0.125, 161.0, 11.4, 5.17, 0.32);
        assert!(lam.a22() > 0.0);
    }
    #[test]
    fn clt_interlaminar_shear_positive() {
        let lam = CompositeLaminateClt::cross_ply(4, 0.125, 161.0, 11.4, 5.17, 0.32);
        let tau = lam.interlaminar_shear_stress(10.0);
        assert!(tau > 0.0);
    }
    #[test]
    fn clt_ye_delamination_rf_no_stress() {
        let lam = CompositeLaminateClt::cross_ply(2, 0.125, 161.0, 11.4, 5.17, 0.32);
        let rf = lam.ye_delamination_rf(0.0, 0.0, 0.0, 50.0, 80.0, 80.0);
        assert!(rf.is_infinite());
    }
    #[test]
    fn clt_ye_delamination_rf_safe() {
        let lam = CompositeLaminateClt::cross_ply(2, 0.125, 161.0, 11.4, 5.17, 0.32);
        let rf = lam.ye_delamination_rf(10.0, 5.0, 5.0, 50.0, 80.0, 80.0);
        assert!(rf > 0.0 && rf.is_finite());
    }
    #[test]
    fn ti64_e_decreases_with_temperature() {
        let ti = TitaniumAlloyExtended::ti64_elevated();
        let e_rt = ti.youngs_modulus_at(20.0);
        let e_hot = ti.youngs_modulus_at(300.0);
        assert!(e_hot < e_rt);
    }
    #[test]
    fn ti64_yield_decreases_with_temperature() {
        let ti = TitaniumAlloyExtended::ti64_elevated();
        let sy_rt = ti.yield_strength_at(20.0);
        let sy_hot = ti.yield_strength_at(300.0);
        assert!(sy_hot < sy_rt);
    }
    #[test]
    fn ti64_creep_rate_positive() {
        let ti = TitaniumAlloyExtended::ti64_elevated();
        let rate = ti.creep_rate(500.0, 773.15);
        assert!(rate > 0.0);
    }
    #[test]
    fn ti6246_higher_strength_than_ti64() {
        let ti64 = TitaniumAlloyExtended::ti64_elevated();
        let ti6246 = TitaniumAlloyExtended::ti6246();
        assert!(ti6246.sigma_y_rt > ti64.sigma_y_rt);
    }
    #[test]
    fn ti64_oxidation_parabolic() {
        let ti = TitaniumAlloyExtended::ti64_elevated();
        let mg1 = ti.oxidation_gain_600(100.0);
        let mg4 = ti.oxidation_gain_600(400.0);
        assert!((mg4 / mg1 - 2.0).abs() < 0.01);
    }
    #[test]
    fn cmsx4_designation() {
        let m = NickelSuperalloyExtended::cmsx4();
        assert_eq!(m.designation, "CMSX-4");
    }
    #[test]
    fn superalloy_oxidation_gain_positive() {
        let m = NickelSuperalloyExtended::cmsx4();
        let mg = m.oxidation_gain(1000.0);
        assert!(mg > 0.0);
    }
    #[test]
    fn superalloy_gamma_prime_strengthening_positive() {
        let mut m = NickelSuperalloyExtended::cmsx4();
        m.gamma_prime_vf = 0.40;
        let delta_sigma = m.gamma_prime_strengthening(200.0, 80.0);
        assert!(delta_sigma > 0.0, "delta_sigma={delta_sigma}");
    }
    #[test]
    fn superalloy_tbc_fraction_between_0_1() {
        let m = NickelSuperalloyExtended::rene142();
        let frac = m.tbc_temperature_fraction(2.0, 15.0, 200.0);
        assert!(frac > 0.0 && frac < 1.0);
    }
    #[test]
    fn cc_effective_conductivity_positive() {
        let cc = CarbonCarbonComposite::cc_2d_woven();
        assert!(cc.effective_conductivity() > 0.0);
    }
    #[test]
    fn cc_oxidation_zero_below_onset() {
        let cc = CarbonCarbonComposite::cc_2d_woven();
        let rec = cc.oxidation_recession_um(400.0, 100.0);
        assert_eq!(rec, 0.0);
    }
    #[test]
    fn cc_oxidation_positive_above_onset() {
        let cc = CarbonCarbonComposite::cc_2d_woven();
        let rec = cc.oxidation_recession_um(600.0, 10.0);
        assert!(rec > 0.0);
    }
    #[test]
    fn cc_strength_retention_high_temp() {
        let cc = CarbonCarbonComposite::cc_3d_nozzle();
        let r = cc.strength_retention(2000.0);
        assert!(r > 0.9);
    }
    #[test]
    fn foam_sandwich_total_thickness() {
        let panel = FoamCoreSandwich::rohacell51_cfrp(1.0, 20.0);
        assert!((panel.total_thickness() - 22.0).abs() < 0.01);
    }
    #[test]
    fn foam_sandwich_bending_stiffness_positive() {
        let panel = FoamCoreSandwich::rohacell51_cfrp(2.0, 30.0);
        assert!(panel.bending_stiffness() > 0.0);
    }
    #[test]
    fn foam_sandwich_max_load_shear_positive() {
        let panel = FoamCoreSandwich::airex_c70_gfrp(2.0, 25.0);
        let q = panel.max_load_shear(500.0);
        assert!(q > 0.0);
    }
    #[test]
    fn foam_sandwich_wrinkling_stress_positive() {
        let panel = FoamCoreSandwich::rohacell51_cfrp(1.5, 20.0);
        assert!(panel.wrinkling_stress_foam() > 0.0);
    }
    #[test]
    fn cmc_ext_modulus_positive() {
        let cmc = CeramicMatrixCompositeExtended::hi_nicalon_sic_sic();
        assert!(cmc.composite_modulus() > 0.0);
    }
    #[test]
    fn cmc_ext_pullout_energy_positive() {
        let cmc = CeramicMatrixCompositeExtended::hi_nicalon_sic_sic();
        assert!(cmc.pullout_energy() > 0.0);
    }
    #[test]
    fn cmc_ext_ultimate_strength_positive() {
        let cmc = CeramicMatrixCompositeExtended::tyranno_sa();
        assert!(cmc.ultimate_strength() > 0.0);
    }
    #[test]
    fn cmc_ext_thermal_shock_positive() {
        let cmc = CeramicMatrixCompositeExtended::hi_nicalon_sic_sic();
        let r = cmc.thermal_shock_r_prime(15.0, 4.0e-6);
        assert!(r > 0.0);
    }
    #[test]
    fn adhesive_rf_shear_safe() {
        let joint = AdhesiveJoint::fm300(0.2, 25.0);
        let rf = joint.rf_shear(20.0);
        assert!(rf > 1.0);
    }
    #[test]
    fn adhesive_volkersen_positive() {
        let joint = AdhesiveJoint::fm300(0.2, 25.0);
        let tau = joint.volkersen_peak_shear(70_000.0, 2.0, 100.0);
        assert!(tau > 0.0);
    }
    #[test]
    fn adhesive_t_peel_energy_positive() {
        let joint = AdhesiveJoint::ec3448(0.25, 20.0);
        let g = joint.t_peel_energy(5.0);
        assert!(g > 0.0);
    }
    #[test]
    fn adhesive_climbing_drum_g1_positive() {
        let joint = AdhesiveJoint::fm300(0.2, 30.0);
        let g = joint.climbing_drum_g1(0.5, 0.025);
        assert!(g > 0.0);
    }
    #[test]
    fn space_kapton_passes_outgassing() {
        let m = SpaceEnvironmentEffects::kapton_hn();
        assert!(m.passes_outgassing());
    }
    #[test]
    fn space_darkening_increases_absorptance() {
        let m = SpaceEnvironmentEffects::kapton_hn();
        let a_initial = m.alpha_initial;
        let a_dark = m.radiation_darkened_absorptance(1.0e6);
        assert!(a_dark >= a_initial);
    }
    #[test]
    fn space_ao_erosion_positive_in_leo() {
        let m = SpaceEnvironmentEffects::kapton_hn();
        let depth = m.ao_erosion_depth_um(5.0, 1.4);
        assert!(depth > 0.0);
    }
    #[test]
    fn space_alpha_over_epsilon_positive() {
        let m = SpaceEnvironmentEffects::solar_black_paint();
        assert!(m.alpha_over_epsilon() > 0.0);
    }
    #[test]
    fn aerogel_knudsen_ambient_less_than_1() {
        let a = AerogelInsulation::silica_aerogel();
        let kn = a.knudsen_number(101_325.0);
        assert!(kn > 0.0);
    }
    #[test]
    fn aerogel_effective_conductivity_positive() {
        let a = AerogelInsulation::silica_aerogel();
        let k = a.effective_conductivity(101_325.0);
        assert!(k > 0.0);
    }
    #[test]
    fn aerogel_vacuum_lower_than_ambient() {
        let a = AerogelInsulation::silica_aerogel();
        let k_vac = a.vacuum_conductivity();
        let k_amb = a.effective_conductivity(101_325.0);
        assert!(k_vac < k_amb);
    }
    #[test]
    fn aerogel_pyrogel_higher_density_than_silica() {
        let silica = AerogelInsulation::silica_aerogel();
        let pyrogel = AerogelInsulation::pyrogel_xt();
        assert!(pyrogel.density > silica.density);
    }
}
