//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Compute Terzaghi bearing capacity factors (Nc, Nq, Ng) for angle φ (radians).
pub fn terzaghi_bearing_factors(phi: f64) -> (f64, f64, f64) {
    if phi < 1e-6 {
        return (5.14, 1.0, 0.0);
    }
    let nq = (PI * phi.tan()).exp() * ((PI / 4.0 + phi / 2.0).tan()).powi(2);
    let nc = (nq - 1.0) / phi.tan();
    let ng = (nq - 1.0) * (1.4 * phi).tan();
    (nc, nq, ng)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::construction::Admixture;
    use crate::construction::AdmixtureType;
    use crate::construction::AggregateGrading;
    use crate::construction::AngleSection;
    use crate::construction::AsphaltMixture;
    use crate::construction::ChannelSection;
    use crate::construction::ConcreteMaterial;
    use crate::construction::ConcreteMixDesign;
    use crate::construction::CoverRequirement;
    use crate::construction::CrossLaminatedTimber;
    use crate::construction::EurocodeLoad;
    use crate::construction::ExposureClass;
    use crate::construction::FoundationSoil;
    use crate::construction::Geogrid;
    use crate::construction::GeosyntheticReinforcement;
    use crate::construction::GeotextileFilter;
    use crate::construction::GlulamSection;
    use crate::construction::HssSection;
    use crate::construction::LaminatedVeneerLumber;
    use crate::construction::MasonryPrism;
    use crate::construction::MasonryShearWall;
    use crate::construction::MasonryUnit;
    use crate::construction::MoistureCorrectionTimber;
    use crate::construction::PileFoundation;
    use crate::construction::PrestressedConcrete;
    use crate::construction::ReinforcedConcrete;
    use crate::construction::RetainingWall;
    use crate::construction::SteelRebarBond;
    use crate::construction::SteelSection;
    use crate::construction::StirrupDesign;
    use crate::construction::StructuralLoad;
    use crate::construction::TerzaghiSettlement;
    use crate::construction::TimberMaterial;
    #[test]
    fn concrete_ec_aci() {
        let c = ConcreteMaterial::new(28.0);
        let expected_ec = 4700.0 * 28.0_f64.sqrt();
        assert!((c.ec - expected_ec).abs() < 1.0);
    }
    #[test]
    fn concrete_ft_is_tenth_fc() {
        let c = ConcreteMaterial::new(30.0);
        assert!((c.ft - 3.0).abs() < 0.01);
    }
    #[test]
    fn concrete_shear_modulus() {
        let c = ConcreteMaterial::new(28.0);
        let g = c.shear_modulus();
        let expected = c.ec / (2.0 * 1.2);
        assert!((g - expected).abs() < 1.0);
    }
    #[test]
    fn concrete_modulus_of_rupture() {
        let c = ConcreteMaterial::new(28.0);
        let fr = c.modulus_of_rupture();
        assert!((fr - 0.62 * 28.0_f64.sqrt()).abs() < 0.01);
    }
    #[test]
    fn concrete_ultimate_compressive_strain() {
        let c = ConcreteMaterial::new(28.0);
        assert!((c.ultimate_compressive_strain() - 0.003).abs() < 1e-10);
    }
    #[test]
    fn rc_moment_capacity_positive() {
        let c = ConcreteMaterial::new(28.0);
        let rc = ReinforcedConcrete::new(c, 300.0, 500.0, 1500.0, 420.0, 0.0, 60.0);
        assert!(rc.moment_capacity() > 0.0);
    }
    #[test]
    fn rc_design_moment_lt_nominal() {
        let c = ConcreteMaterial::new(28.0);
        let rc = ReinforcedConcrete::new(c, 300.0, 500.0, 1500.0, 420.0, 0.0, 60.0);
        assert!(rc.design_moment_capacity() < rc.moment_capacity());
    }
    #[test]
    fn rc_reinforcement_ratio() {
        let c = ConcreteMaterial::new(28.0);
        let rc = ReinforcedConcrete::new(c, 300.0, 500.0, 600.0, 420.0, 0.0, 60.0);
        let rho = rc.reinforcement_ratio();
        assert!((rho - 600.0 / 150000.0).abs() < 1e-10);
    }
    #[test]
    fn rc_shear_capacity_positive() {
        let c = ConcreteMaterial::new(28.0);
        let rc = ReinforcedConcrete::new(c, 300.0, 500.0, 1500.0, 420.0, 0.0, 60.0);
        assert!(rc.shear_capacity() > 0.0);
    }
    #[test]
    fn rc_rho_min() {
        let c = ConcreteMaterial::new(28.0);
        let rc = ReinforcedConcrete::new(c, 300.0, 500.0, 1500.0, 420.0, 0.0, 60.0);
        let rmin = rc.rho_min();
        assert!(rmin > 0.0 && rmin < 0.01);
    }
    #[test]
    fn steel_i_section_area() {
        let s = SteelSection::i_section(150.0, 10.0, 300.0, 8.0, 250.0);
        assert!((s.area - (2.0 * 150.0 * 10.0 + 280.0 * 8.0)).abs() < 1.0);
    }
    #[test]
    fn steel_plastic_moment_positive() {
        let s = SteelSection::i_section(200.0, 12.0, 400.0, 10.0, 250.0);
        assert!(s.plastic_moment() > 0.0);
    }
    #[test]
    fn steel_rx_ry_positive() {
        let s = SteelSection::i_section(200.0, 12.0, 400.0, 10.0, 250.0);
        assert!(s.rx() > 0.0);
        assert!(s.ry() > 0.0);
    }
    #[test]
    fn steel_buckling_stress_decreases_with_slenderness() {
        let s = SteelSection::i_section(200.0, 12.0, 400.0, 10.0, 250.0);
        let fe_50 = s.elastic_buckling_stress(50.0);
        let fe_100 = s.elastic_buckling_stress(100.0);
        assert!(fe_50 > fe_100);
    }
    #[test]
    fn timber_adjusted_fb() {
        let t = TimberMaterial::douglas_fir_ss();
        let fb_adj = t.adjusted_fb(1.0, 1.0);
        assert!((fb_adj - t.mor).abs() < 0.01);
    }
    #[test]
    fn timber_adjusted_fb_with_duration() {
        let t = TimberMaterial::douglas_fir_ss();
        let fb_adj = t.adjusted_fb(1.25, 0.85);
        assert!((fb_adj - t.mor * 1.25 * 0.85).abs() < 0.01);
    }
    #[test]
    fn timber_shear_modulus_lr() {
        let t = TimberMaterial::douglas_fir_ss();
        assert!((t.shear_modulus_lr() - t.e_longitudinal / 16.0).abs() < 0.01);
    }
    #[test]
    fn masonry_em_is_700_fm() {
        let m = MasonryUnit::clay_brick(20.0);
        assert!((m.em - 14000.0).abs() < 1.0);
    }
    #[test]
    fn masonry_shear_modulus() {
        let m = MasonryUnit::clay_brick(20.0);
        assert!((m.shear_modulus() - 0.4 * m.em).abs() < 1.0);
    }
    #[test]
    fn asphalt_air_voids_criterion() {
        let a = AsphaltMixture::superpave_12_5();
        assert!(a.meets_air_voids_criterion());
    }
    #[test]
    fn asphalt_bulk_density() {
        let a = AsphaltMixture::superpave_12_5();
        let gmb = a.bulk_density(2.5);
        assert!(gmb > 0.0 && gmb < 2.5);
    }
    #[test]
    fn foundation_bearing_capacity_positive() {
        let soil = FoundationSoil::medium_clay();
        let qu = soil.ultimate_bearing_capacity(1.5, 1.0);
        assert!(qu > 0.0);
    }
    #[test]
    fn foundation_allowable_lt_ultimate() {
        let soil = FoundationSoil::medium_clay();
        let qu = soil.ultimate_bearing_capacity(1.5, 1.0);
        let qa = soil.allowable_bearing_capacity(1.5, 1.0);
        assert!((qu / 3.0 - qa).abs() < 0.01);
    }
    #[test]
    fn foundation_immediate_settlement_positive() {
        let soil = FoundationSoil::medium_clay();
        let si = soil.immediate_settlement(100.0, 1.5, 5.0, 0.45);
        assert!(si > 0.0);
    }
    #[test]
    fn terzaghi_factors_cohesive() {
        let (nc, nq, ng) = terzaghi_bearing_factors(0.0);
        assert!((nc - 5.14).abs() < 0.01);
        assert!((nq - 1.0).abs() < 0.01);
        assert!(ng.abs() < 0.01);
    }
    #[test]
    fn retaining_wall_rankine_ka_phi30() {
        let w = RetainingWall::new(4.0, 18.0, 30.0);
        let ka = w.rankine_ka();
        assert!((ka - 1.0 / 3.0).abs() < 0.01);
    }
    #[test]
    fn retaining_wall_kp_gt_ka() {
        let w = RetainingWall::new(4.0, 18.0, 30.0);
        assert!(w.rankine_kp() > w.rankine_ka());
    }
    #[test]
    fn retaining_wall_active_thrust_positive() {
        let w = RetainingWall::new(4.0, 18.0, 30.0);
        assert!(w.active_thrust_rankine() > 0.0);
    }
    #[test]
    fn rebar_development_length_positive() {
        let r = SteelRebarBond::new(20.0, 420.0, 28.0, 40.0);
        assert!(r.development_length() > 0.0);
    }
    #[test]
    fn rebar_hook_development_min_criteria() {
        let r = SteelRebarBond::new(20.0, 420.0, 28.0, 40.0);
        let ldh = r.hook_development_length();
        assert!(ldh >= 8.0 * 20.0);
        assert!(ldh >= 150.0);
    }
    #[test]
    fn structural_load_lrfd_combo1() {
        let l = StructuralLoad::new(10.0, 5.0, 2.0, 3.0, 1.0);
        assert!((l.lrfd_combo1() - 14.0).abs() < 0.01);
    }
    #[test]
    fn structural_load_lrfd_combo2() {
        let l = StructuralLoad::new(10.0, 5.0, 2.0, 3.0, 1.0);
        let expected = 12.0 + 8.0 + 0.5;
        assert!((l.lrfd_combo2() - expected).abs() < 0.01);
    }
    #[test]
    fn structural_load_governing_lrfd() {
        let l = StructuralLoad::new(10.0, 5.0, 2.0, 3.0, 1.0);
        let gov = l.governing_lrfd();
        assert!(gov >= l.lrfd_combo1());
        assert!(gov >= l.lrfd_combo2());
    }
    #[test]
    fn structural_load_asd_combo() {
        let l = StructuralLoad::new(10.0, 5.0, 0.0, 0.0, 0.0);
        assert!((l.asd_combo_dl() - 15.0).abs() < 0.01);
    }
    #[test]
    fn mix_design_water_content() {
        let m = ConcreteMixDesign::new(0.5, 350.0, 720.0, 1080.0, 2.0);
        assert!((m.water - 175.0).abs() < 0.01);
    }
    #[test]
    fn mix_design_abrams_strength() {
        let m = ConcreteMixDesign::c30_normal_weight();
        let fc = m.abrams_strength();
        assert!(fc > 20.0 && fc < 60.0, "fc={fc}");
    }
    #[test]
    fn mix_design_aggregate_cement_ratio() {
        let m = ConcreteMixDesign::c30_normal_weight();
        let acr = m.aggregate_cement_ratio();
        assert!(acr > 4.0 && acr < 8.0, "acr={acr}");
    }
    #[test]
    fn mix_design_meets_wc_limit() {
        let m = ConcreteMixDesign::c30_normal_weight();
        assert!(m.meets_durability_wc_limit(0.55));
        assert!(!m.meets_durability_wc_limit(0.45));
    }
    #[test]
    fn mix_design_paste_volume() {
        let m = ConcreteMixDesign::c30_normal_weight();
        let pv = m.paste_volume();
        assert!(pv > 0.1 && pv < 0.5, "pv={pv}");
    }
    #[test]
    fn mix_design_high_performance_slump() {
        let m = ConcreteMixDesign::high_performance();
        assert!(m.slump > 100.0, "slump={}", m.slump);
    }
    #[test]
    fn admixture_adjusted_water() {
        let a = Admixture::new(AdmixtureType::WaterReducer, 0.3, 15.0, 0.0);
        let w = a.adjusted_water(175.0);
        assert!((w - 175.0 * 0.85).abs() < 0.01);
    }
    #[test]
    fn admixture_adjusted_wc_ratio() {
        let a = Admixture::new(AdmixtureType::Superplasticiser, 1.5, 25.0, 0.0);
        let wc = a.adjusted_wc_ratio(0.5);
        assert!((wc - 0.375).abs() < 0.001);
    }
    #[test]
    fn admixture_is_accelerator() {
        let a = Admixture::new(AdmixtureType::Accelerator, 2.0, 0.0, -30.0);
        assert!(a.is_accelerator());
        let b = Admixture::new(AdmixtureType::Retarder, 0.5, 0.0, 60.0);
        assert!(!b.is_accelerator());
    }
    #[test]
    fn aggregate_fineness_modulus() {
        let g = AggregateGrading::astm_c33_fine();
        let fm = g.fineness_modulus();
        assert!(fm > 2.0 && fm < 4.0, "fm={fm}");
    }
    #[test]
    fn aggregate_max_size() {
        let g = AggregateGrading::astm_c33_coarse();
        let ms = g.maximum_size();
        assert!(ms > 0.0, "ms={ms}");
    }
    #[test]
    fn aggregate_combined_fm() {
        let m = ConcreteMixDesign::c30_normal_weight();
        let fm = m.combined_fineness_modulus(2.8, 7.0);
        assert!(fm.is_finite(), "fm={fm}");
    }
    #[test]
    fn prestressed_initial_force() {
        let ps = PrestressedConcrete::new(
            200_000.0, 1.0e10, 300.0, 200.0, 1200.0, 1860.0, 1200.0, 150.0, 40.0, 12_000.0,
        );
        let pi = ps.initial_prestress_force();
        assert!(pi > 0.0, "pi={pi}");
    }
    #[test]
    fn prestressed_effective_force_after_loss() {
        let ps = PrestressedConcrete::new(
            200_000.0, 1.0e10, 300.0, 200.0, 1200.0, 1860.0, 1200.0, 150.0, 40.0, 12_000.0,
        );
        let pe = ps.effective_prestress_force(0.2);
        assert!(pe < ps.initial_prestress_force(), "pe={pe}");
    }
    #[test]
    fn prestressed_cracking_moment_positive() {
        let ps = PrestressedConcrete::new(
            200_000.0, 1.0e10, 300.0, 200.0, 1200.0, 1860.0, 1200.0, 150.0, 40.0, 12_000.0,
        );
        let pe = ps.effective_prestress_force(0.2);
        let mcr = ps.cracking_moment(pe);
        assert!(mcr > 0.0, "mcr={mcr}");
    }
    #[test]
    fn prestressed_losses_positive() {
        let ps = PrestressedConcrete::new(
            200_000.0, 1.0e10, 300.0, 200.0, 1200.0, 1860.0, 1200.0, 150.0, 40.0, 12_000.0,
        );
        let losses = ps.total_losses();
        assert!(losses > 0.0, "losses={losses}");
    }
    #[test]
    fn stirrup_phi_vn_positive() {
        let s = StirrupDesign::new(300.0, 500.0, 28.0, 420.0, 71.0, 2, 150.0);
        assert!(s.phi_vn() > 0.0);
    }
    #[test]
    fn stirrup_vn_gt_vc() {
        let s = StirrupDesign::new(300.0, 500.0, 28.0, 420.0, 71.0, 2, 150.0);
        assert!(s.vn() > s.vc());
    }
    #[test]
    fn stirrup_max_spacing() {
        let s = StirrupDesign::new(300.0, 500.0, 28.0, 420.0, 71.0, 2, 150.0);
        assert!(s.max_spacing() <= 600.0);
        assert!(s.max_spacing() <= s.d / 2.0 + 0.1);
    }
    #[test]
    fn cover_interior_small_bar() {
        let c = CoverRequirement::minimum_cover(&ExposureClass::Interior, 12.0);
        assert!((c - 20.0).abs() < 0.01);
    }
    #[test]
    fn cover_severe_exceeds_interior() {
        let c_int = CoverRequirement::minimum_cover(&ExposureClass::Interior, 20.0);
        let c_sev = CoverRequirement::minimum_cover(&ExposureClass::Severe, 20.0);
        assert!(c_sev > c_int);
    }
    #[test]
    fn cover_effective_depth() {
        let d = CoverRequirement::effective_depth(600.0, 40.0, 20.0);
        assert!((d - 550.0).abs() < 0.01);
    }
    #[test]
    fn glulam_area() {
        let g = GlulamSection::df_24f_v4();
        assert!((g.area() - 275.0 * 570.0).abs() < 1.0);
    }
    #[test]
    fn glulam_allowable_moment_positive() {
        let g = GlulamSection::df_24f_v4();
        let m = g.allowable_moment(1.0, 0.9);
        assert!(m > 0.0, "m={m}");
    }
    #[test]
    fn glulam_midspan_deflection_positive() {
        let g = GlulamSection::df_24f_v4();
        let delta = g.midspan_deflection(2.0, 9_000.0);
        assert!(delta > 0.0, "delta={delta}");
    }
    #[test]
    fn glulam_allowable_shear_positive() {
        let g = GlulamSection::df_24f_v4();
        assert!(g.allowable_shear() > 0.0);
    }
    #[test]
    fn clt_effective_stiffness_positive() {
        let clt = CrossLaminatedTimber::five_layer(1000.0, 175.0);
        let ei = clt.effective_bending_stiffness();
        assert!(ei > 0.0, "ei={ei}");
    }
    #[test]
    fn clt_axial_stiffness_positive() {
        let clt = CrossLaminatedTimber::five_layer(1000.0, 175.0);
        let ea = clt.effective_axial_stiffness();
        assert!(ea > 0.0, "ea={ea}");
    }
    #[test]
    fn lvl_euler_buckling_decreases_with_length() {
        let lvl = LaminatedVeneerLumber::microllam_1_9e(89.0, 302.0);
        let p1 = lvl.euler_buckling_load(3000.0);
        let p2 = lvl.euler_buckling_load(6000.0);
        assert!(p1 > p2, "p1={p1}, p2={p2}");
    }
    #[test]
    fn lvl_d_to_b_ratio() {
        let lvl = LaminatedVeneerLumber::microllam_1_9e(89.0, 302.0);
        let ratio = lvl.d_to_b_ratio();
        assert!(ratio > 1.0 && ratio < 10.0, "ratio={ratio}");
    }
    #[test]
    fn moisture_cm_e_dry() {
        let m = MoistureCorrectionTimber::douglas_fir();
        let cm = m.cm_factor_e();
        assert!(cm > 0.8 && cm <= 1.0, "cm={cm}");
    }
    #[test]
    fn moisture_dimensional_change_positive() {
        let m = MoistureCorrectionTimber::douglas_fir();
        let delta = m.dimensional_change(150.0);
        assert!(delta >= 0.0, "delta={delta}");
    }
    #[test]
    fn masonry_prism_ht_correction() {
        let p = MasonryPrism::new(30.0, "S", 2.0, 18.0);
        let f = p.ht_correction_factor();
        assert!((f - 0.75).abs() < 0.01, "f={f}");
    }
    #[test]
    fn masonry_prism_em_is_900_fm() {
        let p = MasonryPrism::new(25.0, "N", 5.0, 20.0);
        let em = p.em();
        let fm = p.fm_corrected();
        assert!((em - 900.0 * fm).abs() < 1.0);
    }
    #[test]
    fn masonry_shear_wall_capacity_positive() {
        let w = MasonryShearWall::new(3000.0, 190.0, 3200.0, 12.0, 400.0, 420.0, 0.2);
        let v = w.in_plane_shear_capacity();
        assert!(v > 0.0, "v={v}");
    }
    #[test]
    fn masonry_shear_wall_aspect_ratio() {
        let w = MasonryShearWall::new(3000.0, 190.0, 3200.0, 12.0, 400.0, 420.0, 0.2);
        let ar = w.aspect_ratio();
        assert!((ar - 3200.0 / 3000.0).abs() < 0.001);
    }
    #[test]
    fn settlement_nc_primary() {
        let s = TerzaghiSettlement::normally_consolidated(0.35, 1.0, 5.0, 50.0, 2.0);
        let sc = s.primary_settlement(50.0);
        assert!(sc > 0.0, "sc={sc}");
    }
    #[test]
    fn settlement_oc_smaller_than_nc() {
        let nc = TerzaghiSettlement::normally_consolidated(0.35, 1.0, 5.0, 50.0, 2.0);
        let oc = TerzaghiSettlement::over_consolidated(0.35, 0.07, 1.0, 5.0, 50.0, 200.0, 2.0);
        let sc_nc = nc.primary_settlement(30.0);
        let sc_oc = oc.primary_settlement(30.0);
        assert!(sc_oc < sc_nc, "sc_oc={sc_oc}, sc_nc={sc_nc}");
    }
    #[test]
    fn settlement_consolidation_time_positive() {
        let s = TerzaghiSettlement::normally_consolidated(0.35, 1.0, 5.0, 50.0, 2.0);
        let t = s.consolidation_time(0.9);
        assert!(t > 0.0, "t={t}");
    }
    #[test]
    fn settlement_ocr_normally_consolidated() {
        let s = TerzaghiSettlement::normally_consolidated(0.35, 1.0, 5.0, 50.0, 2.0);
        assert!((s.ocr() - 1.0).abs() < 0.001);
    }
    #[test]
    fn geotextile_retention_criterion() {
        let g = GeotextileFilter::woven_road_drainage();
        assert!(g.retention_criterion_met(1.0));
        assert!(!g.retention_criterion_met(0.5));
    }
    #[test]
    fn geotextile_filter_ratio() {
        let g = GeotextileFilter::woven_road_drainage();
        let fr = g.filter_ratio();
        assert!(fr.is_finite() && fr > 0.0, "fr={fr}");
    }
    #[test]
    fn geogrid_allowable_strength() {
        let g = Geogrid::bx1100();
        let ta = g.allowable_strength(1.2, 1.3);
        assert!(ta > 0.0 && ta < g.ltds, "ta={ta}");
    }
    #[test]
    fn geogrid_ux1500hs_tult_md() {
        let g = Geogrid::ux1500hs();
        assert!(g.tult_md > 100.0, "tult_md={}", g.tult_md);
    }
    #[test]
    fn geogrid_interaction_coefficient_is_real_ci() {
        // Ci = tan(δ) / tan(φ_soil); must vary with φ_soil and NOT be the old
        // constant 0.8 (the previous body was 0.8·tan(φ)/tan(φ) ≡ 0.8).
        let g = Geogrid::bx1100(); // δ = 30°
        let ci_30 = g.interaction_coefficient(30.0);
        let ci_36 = g.interaction_coefficient(36.0);
        let ci_25 = g.interaction_coefficient(25.0);

        // δ = φ_soil = 30° ⇒ Ci = 1.0 (full friction mobilised).
        assert!((ci_30 - 1.0).abs() < 1e-9, "ci_30={ci_30}");
        // Stiffer soil (φ > δ) ⇒ Ci < 1; weaker soil (φ < δ) ⇒ Ci > 1.
        assert!(ci_36 < ci_30, "ci_36={ci_36} ci_30={ci_30}");
        assert!(ci_25 > ci_30, "ci_25={ci_25} ci_30={ci_30}");
        // Genuinely input-dependent: the three values differ from each other.
        assert!((ci_36 - ci_25).abs() > 0.05, "ci_36={ci_36} ci_25={ci_25}");

        // Closed-form check against tan(δ)/tan(φ).
        let expected = (30f64).to_radians().tan() / (36f64).to_radians().tan();
        assert!(
            (ci_36 - expected).abs() < 1e-9,
            "ci_36={ci_36} exp={expected}"
        );

        // Must NOT be pinned at the fabricated constant 0.8 for every input.
        assert!(
            (ci_36 - 0.8).abs() > 0.01 || (ci_25 - 0.8).abs() > 0.01,
            "Ci collapsed to the old constant 0.8"
        );

        // Different geogrid (different δ) ⇒ different Ci at the same soil angle.
        let g2 = Geogrid::ux1500hs(); // δ = 32°
        let ci2_30 = g2.interaction_coefficient(30.0);
        assert!(
            (ci2_30 - ci_30).abs() > 0.01,
            "ci2_30={ci2_30} ci_30={ci_30}"
        );

        // Degenerate soil angle is handled without NaN/∞.
        let ci_zero = g.interaction_coefficient(0.0);
        assert!(ci_zero.is_finite(), "ci_zero={ci_zero}");
        assert_eq!(ci_zero, 0.0, "ci_zero={ci_zero}");
        assert!(
            g.interaction_coefficient(-5.0).is_finite(),
            "negative φ_soil must stay finite"
        );
    }
    #[test]
    fn geosynthetic_horizontal_stress_increases_with_depth() {
        let r = GeosyntheticReinforcement::new("geogrid", 10.0, 0.6, 34.0, 20.0, 0.9);
        let s1 = r.horizontal_stress(3.0);
        let s2 = r.horizontal_stress(6.0);
        assert!(s2 > s1, "s2={s2}, s1={s1}");
    }
    #[test]
    fn geosynthetic_tension_fos_infinite_at_surface() {
        let r = GeosyntheticReinforcement::new("geogrid", 10.0, 0.6, 34.0, 20.0, 0.9);
        let fos = r.tension_fos(0.0);
        assert!(fos.is_infinite() || fos > 100.0, "fos={fos}");
    }
    #[test]
    fn eurocode_uls_610_dominates_dead_live() {
        let e = EurocodeLoad::new(10.0, 5.0, 2.0, 1.0, 0.5, 0.0);
        let uls = e.uls_combo_610();
        assert!(uls > 10.0 * 1.35, "uls={uls}");
    }
    #[test]
    fn eurocode_governing_uls_geq_each_combo() {
        let e = EurocodeLoad::new(10.0, 5.0, 2.0, 1.0, 0.5, 0.0);
        let gov = e.governing_uls();
        assert!(gov >= e.uls_combo_610());
        assert!(gov >= e.uls_combo_610a());
        assert!(gov >= e.uls_combo_610b());
    }
    #[test]
    fn eurocode_sls_characteristic_less_than_uls() {
        let e = EurocodeLoad::new(10.0, 5.0, 2.0, 1.0, 0.5, 0.0);
        assert!(e.sls_characteristic() < e.uls_combo_610());
    }
    #[test]
    fn eurocode_quasi_permanent_less_than_frequent() {
        let e = EurocodeLoad::new(10.0, 5.0, 2.0, 1.0, 0.5, 0.0);
        assert!(e.sls_quasi_permanent() <= e.sls_frequent());
    }
    #[test]
    fn hss_square_area() {
        let h = HssSection::square(100.0, 6.0, 350.0);
        let expected = 100.0 * 100.0 - 88.0 * 88.0;
        assert!((h.area() - expected).abs() < 1.0);
    }
    #[test]
    fn hss_rectangular_ix_gt_iy() {
        let h = HssSection::rectangular(100.0, 200.0, 6.0, 350.0);
        assert!(h.ix() > h.iy(), "ix={}, iy={}", h.ix(), h.iy());
    }
    #[test]
    fn hss_compact_flange_check() {
        let h = HssSection::square(100.0, 6.0, 350.0);
        let compact = h.flange_is_compact();
        let _ = compact;
        assert!(h.compact_flange_limit() > 0.0);
    }
    #[test]
    fn hss_plastic_moment_positive() {
        let h = HssSection::rectangular(150.0, 250.0, 8.0, 350.0);
        assert!(h.plastic_moment() > 0.0);
    }
    #[test]
    fn angle_area_positive() {
        let a = AngleSection::equal_leg(100.0, 10.0, 250.0);
        assert!(a.area() > 0.0, "area={}", a.area());
    }
    #[test]
    fn angle_centroid_between_t_and_l() {
        let a = AngleSection::equal_leg(100.0, 10.0, 250.0);
        let cx = a.centroid();
        assert!(cx > 10.0 && cx < 100.0, "cx={cx}");
    }
    #[test]
    fn channel_area_positive() {
        let c = ChannelSection::new(230.0, 67.0, 11.0, 7.0, 250.0);
        assert!(c.area() > 0.0);
    }
    #[test]
    fn channel_moment_capacity_positive() {
        let c = ChannelSection::new(230.0, 67.0, 11.0, 7.0, 250.0);
        assert!(c.moment_capacity() > 0.0);
    }
    #[test]
    fn channel_shear_capacity_positive() {
        let c = ChannelSection::new(230.0, 67.0, 11.0, 7.0, 250.0);
        assert!(c.shear_capacity() > 0.0);
    }
    #[test]
    fn pile_ultimate_capacity_positive() {
        let p = PileFoundation::concrete_pile(0.5, 15.0, 40.0, 800.0);
        assert!(p.ultimate_capacity() > 0.0);
    }
    #[test]
    fn pile_allowable_lt_ultimate() {
        let p = PileFoundation::concrete_pile(0.5, 15.0, 40.0, 800.0);
        assert!(p.allowable_capacity() < p.ultimate_capacity());
    }
    #[test]
    fn pile_elastic_compression_positive() {
        let p = PileFoundation::concrete_pile(0.5, 15.0, 40.0, 800.0);
        let delta = p.elastic_compression(500.0);
        assert!(delta > 0.0, "delta={delta}");
    }
    #[test]
    fn pile_skin_friction_plus_end_bearing_is_ultimate() {
        let p = PileFoundation::concrete_pile(0.5, 15.0, 40.0, 800.0);
        let total = p.skin_friction() + p.end_bearing();
        assert!((total - p.ultimate_capacity()).abs() < 0.01);
    }
}
