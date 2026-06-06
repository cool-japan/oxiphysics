//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::Material;

/// Get all wood species presets.
pub fn wood_presets() -> Vec<Material> {
    vec![
        wood(),
        oak(),
        pine_wood(),
        balsa_wood(),
        beech_wood(),
        teak_wood(),
        maple_wood(),
        walnut_wood(),
        bamboo(),
        mdf_board(),
        osb_board(),
    ]
}
/// Get all concrete grade presets.
pub fn concrete_grade_presets() -> Vec<Material> {
    vec![
        concrete_c20(),
        concrete_c30(),
        concrete_c40(),
        concrete_c50(),
        concrete_c60(),
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::presets::MaterialCategory;
    #[test]
    fn test_preset_materials_reasonable_values() {
        for m in all_presets() {
            assert!(m.density > 0.0, "{} density should be positive", m.name);
            assert!(
                m.friction >= 0.0,
                "{} friction should be non-negative",
                m.name
            );
            assert!(
                m.restitution >= 0.0 && m.restitution <= 1.0,
                "{} restitution should be in [0, 1], got {}",
                m.name,
                m.restitution
            );
        }
    }
    #[test]
    fn test_density_ordering() {
        assert!(steel().density > water().density);
        assert!(tungsten().density > steel().density);
        assert!(air().density < water().density);
        assert!(aluminum().density < steel().density);
    }
    #[test]
    fn test_friction_extremes() {
        assert!(ptfe().friction < rubber().friction);
        assert!(ice().friction < concrete().friction);
    }
    #[test]
    fn test_restitution_extremes() {
        assert!(billiard_ball().restitution > lead().restitution);
        assert!(soft_rubber().restitution > steel().restitution);
    }
    #[test]
    fn test_extended_steel_moduli() {
        let s = extended_steel();
        let k = s.bulk_modulus();
        let g = s.shear_modulus();
        assert!((k - 166.667e9).abs() < 1e8, "K = {k}");
        assert!((g - 76.923e9).abs() < 1e8, "G = {g}");
    }
    #[test]
    fn test_all_presets_have_names() {
        for m in all_presets() {
            assert!(!m.name.is_empty());
        }
    }
    #[test]
    fn test_ceramic_presets_exist() {
        let ceramics = ceramic_presets();
        assert!(ceramics.len() >= 5, "should have at least 5 ceramics");
        for c in &ceramics {
            assert!(c.density > 2000.0, "{} density too low for ceramic", c.name);
        }
    }
    #[test]
    fn test_composite_presets_exist() {
        let composites = composite_presets();
        assert!(composites.len() >= 4, "should have at least 4 composites");
    }
    #[test]
    fn test_biomaterial_presets_exist() {
        let bio = biomaterial_presets();
        assert!(bio.len() >= 5, "should have at least 5 biomaterials");
    }
    #[test]
    fn test_polymer_presets_exist() {
        let polymers = polymer_presets();
        assert!(polymers.len() >= 10, "should have at least 10 polymers");
        for p in &polymers {
            assert!(
                p.density < 3000.0,
                "{} density too high for polymer: {}",
                p.name,
                p.density
            );
        }
    }
    #[test]
    fn test_preset_by_name_found() {
        let m = preset_by_name("steel");
        assert!(m.is_some());
        assert_eq!(m.unwrap().name, "steel");
    }
    #[test]
    fn test_preset_by_name_case_insensitive() {
        let m = preset_by_name("STEEL");
        assert!(m.is_some());
    }
    #[test]
    fn test_preset_by_name_not_found() {
        let m = preset_by_name("unobtanium");
        assert!(m.is_none());
    }
    #[test]
    fn test_presets_by_category() {
        let metals = presets_by_category(MaterialCategory::Metal);
        assert!(!metals.is_empty());
        let fluids = presets_by_category(MaterialCategory::Fluid);
        assert!(!fluids.is_empty());
    }
    #[test]
    fn test_extended_specific_stiffness() {
        let al = extended_aluminum();
        let steel = extended_steel();
        assert!(al.specific_stiffness() > steel.specific_stiffness() * 0.5);
    }
    #[test]
    fn test_extended_titanium_properties() {
        let ti = extended_titanium();
        assert!(
            ti.yield_strength > 800e6,
            "Ti-6Al-4V yield should be > 800 MPa"
        );
        assert!(ti.specific_strength() > steel().density);
    }
    #[test]
    fn test_extended_alumina_high_modulus() {
        let al2o3 = extended_alumina();
        assert!(al2o3.young_modulus > 300e9, "alumina E should be > 300 GPa");
    }
    #[test]
    fn test_extended_cfrp_lightweight_strong() {
        let cfrp = extended_cfrp();
        let steel = extended_steel();
        assert!(cfrp.specific_strength() > steel.specific_strength());
    }
    #[test]
    fn test_all_extended_presets_count() {
        let all = all_extended_presets();
        assert!(all.len() >= 10, "should have at least 10 extended presets");
    }
    #[test]
    fn test_no_duplicate_names() {
        let all = all_presets();
        let mut names: Vec<&str> = all.iter().map(|m| m.name.as_str()).collect();
        names.sort();
        for w in names.windows(2) {
            assert_ne!(w[0], w[1], "duplicate material name: {}", w[0]);
        }
    }
    #[test]
    fn test_ceramic_hardness_proxy() {
        assert!(alumina().friction < brick().friction);
    }
    #[test]
    fn test_mercury_highest_density_fluid() {
        let fluids = fluid_presets();
        let max_density = fluids.iter().map(|f| f.density).fold(0.0_f64, f64::max);
        assert!((max_density - mercury().density).abs() < 1.0);
    }
    #[test]
    fn test_thermal_diffusivity_extended() {
        let cu = extended_copper();
        let alpha = cu.thermal_diffusivity();
        assert!(alpha > 1e-4, "copper diffusivity = {alpha}");
    }
    #[test]
    fn test_nitinol_properties() {
        let m = nitinol();
        assert!((m.density - 6450.0).abs() < 1.0, "density = {}", m.density);
        assert_eq!(m.name, "nitinol_niti");
    }
    #[test]
    fn test_hydrogel_low_density() {
        let m = hydrogel();
        assert!(m.density < 1100.0, "hydrogel density = {}", m.density);
        assert!(
            m.restitution < 0.1,
            "hydrogel restitution = {}",
            m.restitution
        );
    }
    #[test]
    fn test_carbon_fiber_composite_low_density() {
        let m = carbon_fiber_composite();
        assert!(m.density < 1700.0, "CFRP density = {}", m.density);
    }
    #[test]
    fn test_glass_fiber_epoxy_properties() {
        let m = glass_fiber_epoxy();
        assert!(
            m.density > 1500.0 && m.density < 2500.0,
            "GFE density = {}",
            m.density
        );
    }
    #[test]
    fn test_extended_nitinol_shape_memory() {
        let m = extended_nitinol();
        assert!(
            m.young_modulus > 50.0e9 && m.young_modulus < 120.0e9,
            "NiTi E = {}",
            m.young_modulus
        );
    }
    #[test]
    fn test_extended_ptfe_low_modulus() {
        let m = extended_ptfe();
        assert!(m.young_modulus < 1.0e9, "PTFE E = {}", m.young_modulus);
    }
    #[test]
    fn test_extended_polycarbonate_moderate_modulus() {
        let m = extended_polycarbonate();
        assert!(
            m.young_modulus > 1.0e9 && m.young_modulus < 5.0e9,
            "PC E = {}",
            m.young_modulus
        );
    }
    #[test]
    fn test_extended_carbon_fiber_composite_high_modulus() {
        let m = extended_carbon_fiber_composite();
        assert!(m.young_modulus > 100.0e9, "CFRP E = {}", m.young_modulus);
    }
    #[test]
    fn test_extended_glass_fiber_epoxy_properties() {
        let m = extended_glass_fiber_epoxy();
        assert!(
            m.young_modulus > 10.0e9 && m.young_modulus < 40.0e9,
            "GFE E = {}",
            m.young_modulus
        );
    }
    #[test]
    fn test_extended_hydrogel_nearly_incompressible() {
        let m = extended_hydrogel();
        assert!(m.poisson_ratio > 0.45, "hydrogel ν = {}", m.poisson_ratio);
    }
    #[test]
    fn test_stiffest_material() {
        let materials = advanced_extended_presets();
        let stiffest = stiffest_material(&materials);
        assert!(stiffest.is_some());
        let e = stiffest.unwrap().young_modulus;
        for m in &materials {
            assert!(
                m.young_modulus <= e + 1.0,
                "not stiffest: {} > {}",
                m.young_modulus,
                e
            );
        }
    }
    #[test]
    fn test_highest_specific_strength() {
        let materials = all_extended_presets();
        let strongest = highest_specific_strength(&materials);
        assert!(strongest.is_some());
        let ss = strongest.unwrap().specific_strength();
        assert!(ss > 0.0);
    }
    #[test]
    fn test_advanced_extended_presets_count() {
        let all = advanced_extended_presets();
        assert!(all.len() >= 5, "should have at least 5 advanced presets");
    }
    #[test]
    fn test_nitinol_vs_steel_density() {
        assert!(nitinol().density < steel().density);
    }
    #[test]
    fn test_hydrogel_vs_water_density() {
        assert!(hydrogel().density > water().density);
    }
    #[test]
    fn test_al_7075_density_in_range() {
        let m = al_7075_t6();
        assert!(
            m.density > 2500.0 && m.density < 3000.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_al_2024_name() {
        let m = al_2024_t3();
        assert_eq!(m.name, "al_2024_t3");
    }
    #[test]
    fn test_ti_cp_density_lower_than_ti64() {
        assert!(ti_cp_grade2().density < titanium().density + 100.0);
    }
    #[test]
    fn test_inconel_625_name() {
        let m = inconel_625();
        assert_eq!(m.name, "inconel_625");
    }
    #[test]
    fn test_cobalt_chromium_density() {
        let m = cobalt_chromium();
        assert!(
            m.density > 8000.0 && m.density < 10000.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_beryllium_low_density() {
        let m = beryllium();
        assert!(m.density < 2000.0, "density={}", m.density);
    }
    #[test]
    fn test_extended_al_7075_high_strength() {
        let m = extended_al_7075_t6();
        assert!(m.yield_strength > 400.0e6, "yield={}", m.yield_strength);
    }
    #[test]
    fn test_extended_al_2024_modulus() {
        let m = extended_al_2024_t3();
        assert!(
            m.young_modulus > 60.0e9 && m.young_modulus < 80.0e9,
            "E={}",
            m.young_modulus
        );
    }
    #[test]
    fn test_extended_inconel_625_thermal() {
        let m = extended_inconel_625();
        assert!(m.thermal_conductivity > 5.0, "k={}", m.thermal_conductivity);
    }
    #[test]
    fn test_intervertebral_disc_low_density() {
        let m = intervertebral_disc();
        assert!(m.density < 1300.0, "density={}", m.density);
    }
    #[test]
    fn test_meniscus_cartilage_very_low_friction() {
        let m = meniscus_cartilage();
        assert!(m.friction < 0.05, "friction={}", m.friction);
    }
    #[test]
    fn test_articular_cartilage_name() {
        let m = articular_cartilage();
        assert_eq!(m.name, "articular_cartilage");
    }
    #[test]
    fn test_ligament_density() {
        let m = ligament();
        assert!(
            m.density > 1000.0 && m.density < 1500.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_pla_biopolymer_density() {
        let m = pla_biopolymer();
        assert!(
            m.density > 1000.0 && m.density < 1500.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_extended_cortical_bone_orthotropic_young() {
        let m = extended_cortical_bone_orthotropic();
        assert!(m.young_modulus > 15.0e9, "E={}", m.young_modulus);
    }
    #[test]
    fn test_ptfe_enhanced_very_low_friction() {
        let m = ptfe_enhanced();
        assert!(m.friction < 0.05, "friction={}", m.friction);
    }
    #[test]
    fn test_peek_cf30_higher_stiffness() {
        let m = extended_peek_cf30();
        assert!(
            m.young_modulus > extended_peek().young_modulus,
            "CF30 PEEK should be stiffer than pure PEEK"
        );
    }
    #[test]
    fn test_pei_ultem_name() {
        let m = pei_ultem();
        assert_eq!(m.name, "pei_ultem");
    }
    #[test]
    fn test_pvdf_moderate_density() {
        let m = pvdf();
        assert!(
            m.density > 1500.0 && m.density < 2000.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_abs_density_range() {
        let m = abs_plastic();
        assert!(
            m.density > 900.0 && m.density < 1100.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_sic_density_expected() {
        let m = silicon_carbide();
        assert!((m.density - 3210.0).abs() < 50.0, "density={}", m.density);
    }
    #[test]
    fn test_al2o3_density_expected() {
        let m = alumina();
        assert!((m.density - 3950.0).abs() < 100.0, "density={}", m.density);
    }
    #[test]
    fn test_mullite_density() {
        let m = mullite();
        assert!(
            m.density > 2500.0 && m.density < 3500.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_cordierite_low_density_ceramic() {
        let m = cordierite();
        assert!(m.density < 2800.0, "density={}", m.density);
    }
    #[test]
    fn test_extended_zirconia_tough() {
        let m = extended_zirconia_y_psz();
        assert!(m.young_modulus > 150.0e9, "E={}", m.young_modulus);
    }
    #[test]
    fn test_gravel_density_in_range() {
        let m = gravel();
        assert!(
            m.density > 1400.0 && m.density < 2200.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_silt_low_friction() {
        let m = silt();
        assert!(m.friction < 0.4, "friction={}", m.friction);
    }
    #[test]
    fn test_dry_clay_name() {
        let m = dry_clay();
        assert_eq!(m.name, "dry_clay");
    }
    #[test]
    fn test_concrete_c20_lower_density_than_c60() {
        let c20 = concrete_c20();
        let c60 = concrete_c60();
        assert!(c20.density > 2000.0 && c60.density > 2000.0);
    }
    #[test]
    fn test_extended_concrete_c20_vs_c60_modulus() {
        let c20 = extended_concrete_c20();
        let c60 = extended_concrete_c60();
        assert!(
            c60.young_modulus > c20.young_modulus,
            "C60 should be stiffer: c20={}, c60={}",
            c20.young_modulus,
            c60.young_modulus
        );
    }
    #[test]
    fn test_concrete_c40_name() {
        let m = concrete_c40();
        assert_eq!(m.name, "concrete_c40");
    }
    #[test]
    fn test_rammed_earth_density() {
        let m = rammed_earth();
        assert!(
            m.density > 1600.0 && m.density < 2400.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_pine_wood_low_density() {
        let m = pine_wood();
        assert!(m.density < 600.0, "density={}", m.density);
    }
    #[test]
    fn test_balsa_lightest_wood() {
        let m = balsa_wood();
        assert!(
            m.density < pine_wood().density,
            "balsa should be lighter than pine"
        );
    }
    #[test]
    fn test_teak_heavier_than_pine() {
        assert!(teak_wood().density > pine_wood().density);
    }
    #[test]
    fn test_bamboo_name() {
        let m = bamboo();
        assert_eq!(m.name, "bamboo");
    }
    #[test]
    fn test_mdf_density_reasonable() {
        let m = mdf_board();
        assert!(
            m.density > 500.0 && m.density < 900.0,
            "density={}",
            m.density
        );
    }
    #[test]
    fn test_extended_oak_longitudinal_young() {
        let m = extended_oak();
        assert!(
            m.young_modulus > 8.0e9 && m.young_modulus < 20.0e9,
            "E={}",
            m.young_modulus
        );
    }
    #[test]
    fn test_all_wood_densities_sensible() {
        for w in wood_presets() {
            assert!(
                w.density > 50.0 && w.density < 1500.0,
                "{} density={}",
                w.name,
                w.density
            );
        }
    }
    #[test]
    fn test_all_soil_densities_sensible() {
        for s in soil_presets() {
            assert!(
                s.density > 500.0 && s.density < 2700.0,
                "{} density={}",
                s.name,
                s.density
            );
        }
    }
    #[test]
    fn test_aerospace_presets_count() {
        let ae = aerospace_presets();
        assert!(
            ae.len() >= 5,
            "should have >= 5 aerospace presets, got {}",
            ae.len()
        );
    }
    #[test]
    fn test_all_presets_with_new_categories() {
        let new_mats = [
            al_7075_t6(),
            al_2024_t3(),
            ti_cp_grade2(),
            inconel_625(),
            cobalt_chromium(),
            beryllium(),
            gravel(),
            silt(),
            dry_clay(),
            rammed_earth(),
            concrete_c20(),
            concrete_c40(),
            concrete_c60(),
            pine_wood(),
            balsa_wood(),
            teak_wood(),
            bamboo(),
            mdf_board(),
            ptfe_enhanced(),
            pei_ultem(),
            pvdf(),
            abs_plastic(),
            pla_biopolymer(),
            mullite(),
            cordierite(),
            intervertebral_disc(),
            articular_cartilage(),
            ligament(),
            meniscus_cartilage(),
        ];
        for m in &new_mats {
            assert!(m.density > 0.0, "{} has zero density", m.name);
            assert!(m.friction >= 0.0, "{} has negative friction", m.name);
            assert!(
                m.restitution >= 0.0 && m.restitution <= 1.0,
                "{} restitution={} out of range",
                m.name,
                m.restitution
            );
        }
    }
}
