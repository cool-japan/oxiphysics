//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::Material;

use super::types::{ExtendedMaterial, MaterialCategory};

/// Steel (structural): high density, moderate friction, low restitution.
pub fn steel() -> Material {
    Material::new("steel", 7800.0, 0.6, 0.3)
}
/// Stainless steel (304): similar to structural steel, slightly lower density.
pub fn stainless_steel() -> Material {
    Material::new("stainless_steel", 7900.0, 0.55, 0.25)
}
/// Cast iron: very high density, moderate friction, low restitution.
pub fn cast_iron() -> Material {
    Material::new("cast_iron", 7200.0, 0.55, 0.2)
}
/// Aluminum (6061 alloy): low density, moderate friction, moderate restitution.
pub fn aluminum() -> Material {
    Material::new("aluminum", 2700.0, 0.47, 0.3)
}
/// Copper: high density, moderate friction, low restitution.
pub fn copper() -> Material {
    Material::new("copper", 8960.0, 0.5, 0.25)
}
/// Brass (Cu-Zn alloy): moderate-high density, moderate friction.
pub fn brass() -> Material {
    Material::new("brass", 8500.0, 0.5, 0.28)
}
/// Titanium (Ti-6Al-4V): moderate density, high strength.
pub fn titanium() -> Material {
    Material::new("titanium", 4430.0, 0.5, 0.32)
}
/// Lead: very high density, low friction, very low restitution.
pub fn lead() -> Material {
    Material::new("lead", 11340.0, 0.4, 0.05)
}
/// Tungsten: highest density common metal, low friction.
pub fn tungsten() -> Material {
    Material::new("tungsten", 19300.0, 0.35, 0.15)
}
/// Magnesium: very low density, moderate friction.
pub fn magnesium() -> Material {
    Material::new("magnesium", 1740.0, 0.45, 0.3)
}
/// Nickel: high density, corrosion resistant.
pub fn nickel() -> Material {
    Material::new("nickel", 8908.0, 0.53, 0.28)
}
/// Zinc: moderate density, used in galvanizing.
pub fn zinc() -> Material {
    Material::new("zinc", 7133.0, 0.45, 0.2)
}
/// Silver: high density, excellent conductivity.
pub fn silver() -> Material {
    Material::new("silver", 10490.0, 0.4, 0.35)
}
/// Gold: very high density, soft, low friction.
pub fn gold() -> Material {
    Material::new("gold", 19320.0, 0.35, 0.2)
}
/// Inconel 718: nickel-based superalloy.
pub fn inconel() -> Material {
    Material::new("inconel_718", 8190.0, 0.5, 0.25)
}
/// Alumina (Al2O3): hard, brittle, high-temperature ceramic.
pub fn alumina() -> Material {
    Material::new("alumina", 3950.0, 0.35, 0.5)
}
/// Silicon carbide (SiC): extremely hard, excellent thermal conductivity.
pub fn silicon_carbide() -> Material {
    Material::new("silicon_carbide", 3210.0, 0.3, 0.45)
}
/// Silicon nitride (Si3N4): high strength, thermal shock resistant.
pub fn silicon_nitride() -> Material {
    Material::new("silicon_nitride", 3200.0, 0.32, 0.48)
}
/// Zirconia (ZrO2): high toughness ceramic, used in dental implants.
pub fn zirconia() -> Material {
    Material::new("zirconia", 5680.0, 0.35, 0.42)
}
/// Boron carbide (B4C): extremely hard, lightweight armour material.
pub fn boron_carbide() -> Material {
    Material::new("boron_carbide", 2520.0, 0.25, 0.4)
}
/// Porcelain: traditional ceramic, moderate properties.
pub fn porcelain() -> Material {
    Material::new("porcelain", 2400.0, 0.4, 0.45)
}
/// Fused silica (SiO2): optical-grade glass ceramic.
pub fn fused_silica() -> Material {
    Material::new("fused_silica", 2200.0, 0.35, 0.5)
}
/// Carbon fibre reinforced polymer (CFRP), unidirectional layup.
pub fn cfrp() -> Material {
    Material::new("cfrp", 1600.0, 0.3, 0.4)
}
/// Glass fibre reinforced polymer (GFRP / fibreglass).
pub fn gfrp() -> Material {
    Material::new("gfrp", 2000.0, 0.35, 0.45)
}
/// Kevlar (aramid) fibre composite.
pub fn kevlar_composite() -> Material {
    Material::new("kevlar_composite", 1440.0, 0.35, 0.4)
}
/// Ceramic matrix composite (CMC): SiC/SiC.
pub fn ceramic_matrix_composite() -> Material {
    Material::new("cmc_sic_sic", 2800.0, 0.3, 0.35)
}
/// Metal matrix composite (MMC): Al-SiC.
pub fn metal_matrix_composite() -> Material {
    Material::new("mmc_al_sic", 2900.0, 0.4, 0.35)
}
/// Plywood (laminated wood composite).
pub fn plywood() -> Material {
    Material::new("plywood", 640.0, 0.45, 0.38)
}
/// Rubber: low density, very high friction, high restitution.
pub fn rubber() -> Material {
    Material::new("rubber", 1100.0, 0.9, 0.8)
}
/// Natural rubber (soft): very high restitution.
pub fn soft_rubber() -> Material {
    Material::new("soft_rubber", 950.0, 0.95, 0.92)
}
/// HDPE (high-density polyethylene): low friction, moderate restitution.
pub fn hdpe() -> Material {
    Material::new("hdpe", 960.0, 0.25, 0.45)
}
/// PTFE (Teflon): extremely low friction.
pub fn ptfe() -> Material {
    Material::new("ptfe", 2200.0, 0.04, 0.35)
}
/// Nylon (PA 6): moderate density, moderate friction.
pub fn nylon() -> Material {
    Material::new("nylon", 1150.0, 0.4, 0.5)
}
/// Polycarbonate: transparent, moderate density.
pub fn polycarbonate() -> Material {
    Material::new("polycarbonate", 1200.0, 0.45, 0.55)
}
/// PVC (rigid): moderate density, moderate friction.
pub fn pvc() -> Material {
    Material::new("pvc", 1380.0, 0.5, 0.4)
}
/// Acrylic (PMMA): similar to glass in optical terms.
pub fn acrylic() -> Material {
    Material::new("acrylic", 1190.0, 0.5, 0.5)
}
/// Polyethylene terephthalate (PET): common plastic bottle material.
pub fn pet() -> Material {
    Material::new("pet", 1380.0, 0.4, 0.45)
}
/// Polypropylene (PP): lightweight, chemical resistant.
pub fn polypropylene() -> Material {
    Material::new("polypropylene", 905.0, 0.3, 0.4)
}
/// Polystyrene (PS): rigid, brittle, lightweight.
pub fn polystyrene() -> Material {
    Material::new("polystyrene", 1050.0, 0.4, 0.35)
}
/// Polyurethane foam: lightweight, energy absorbing.
pub fn polyurethane_foam() -> Material {
    Material::new("polyurethane_foam", 30.0, 0.6, 0.1)
}
/// Epoxy resin: thermoset, used as matrix in composites.
pub fn epoxy() -> Material {
    Material::new("epoxy", 1200.0, 0.5, 0.35)
}
/// PEEK (polyetheretherketone): high-performance engineering polymer.
pub fn peek() -> Material {
    Material::new("peek", 1300.0, 0.4, 0.45)
}
/// Bone (cortical): hard, brittle, high density for biological material.
pub fn cortical_bone() -> Material {
    Material::new("cortical_bone", 1900.0, 0.4, 0.3)
}
/// Cancellous (trabecular/spongy) bone: porous, lower density.
pub fn cancellous_bone() -> Material {
    Material::new("cancellous_bone", 800.0, 0.35, 0.2)
}
/// Cartilage: soft, low friction (synovial joint).
pub fn cartilage() -> Material {
    Material::new("cartilage", 1200.0, 0.01, 0.4)
}
/// Silicone (medical grade): biocompatible elastomer.
pub fn silicone_medical() -> Material {
    Material::new("silicone_medical", 1100.0, 0.6, 0.75)
}
/// Hydroxyapatite: mineral component of bone, used in implants.
pub fn hydroxyapatite() -> Material {
    Material::new("hydroxyapatite", 3160.0, 0.3, 0.35)
}
/// Dental enamel: hardest biological tissue.
pub fn dental_enamel() -> Material {
    Material::new("dental_enamel", 2950.0, 0.35, 0.4)
}
/// Collagen (tendon/ligament): flexible biological material.
pub fn collagen() -> Material {
    Material::new("collagen", 1300.0, 0.2, 0.3)
}
/// Bioglass (45S5): bioactive glass for bone repair.
pub fn bioglass() -> Material {
    Material::new("bioglass_45s5", 2700.0, 0.35, 0.4)
}
/// UHMWPE (ultra-high molecular weight polyethylene): used in joint replacements.
pub fn uhmwpe() -> Material {
    Material::new("uhmwpe", 930.0, 0.1, 0.4)
}
/// Wood (generic): low density, moderate friction.
pub fn wood() -> Material {
    Material::new("wood", 600.0, 0.5, 0.4)
}
/// Oak (hardwood): denser than generic wood, lower restitution.
pub fn oak() -> Material {
    Material::new("oak", 750.0, 0.55, 0.35)
}
/// Glass: moderate density, low friction, moderate restitution.
pub fn glass() -> Material {
    Material::new("glass", 2500.0, 0.4, 0.5)
}
/// Safety glass (tempered).
pub fn tempered_glass() -> Material {
    Material::new("tempered_glass", 2500.0, 0.4, 0.5)
}
/// Concrete (normal weight).
pub fn concrete() -> Material {
    Material::new("concrete", 2400.0, 0.65, 0.2)
}
/// Reinforced concrete.
pub fn reinforced_concrete() -> Material {
    Material::new("reinforced_concrete", 2500.0, 0.65, 0.2)
}
/// Brick (fired clay).
pub fn brick() -> Material {
    Material::new("brick", 1900.0, 0.7, 0.15)
}
/// Granite (natural stone).
pub fn granite() -> Material {
    Material::new("granite", 2700.0, 0.75, 0.2)
}
/// Marble.
pub fn marble() -> Material {
    Material::new("marble", 2700.0, 0.45, 0.3)
}
/// Asphalt (road surface).
pub fn asphalt() -> Material {
    Material::new("asphalt", 2300.0, 0.9, 0.1)
}
/// Dry sand.
pub fn dry_sand() -> Material {
    Material::new("dry_sand", 1600.0, 0.6, 0.1)
}
/// Wet clay.
pub fn wet_clay() -> Material {
    Material::new("wet_clay", 1800.0, 0.3, 0.05)
}
/// Water: standard density (1000 kg/m^3).
pub fn water() -> Material {
    Material::new("water", 1000.0, 0.0, 0.0)
}
/// Seawater (3.5% salinity).
pub fn seawater() -> Material {
    Material::new("seawater", 1025.0, 0.0, 0.0)
}
/// Air (at STP).
pub fn air() -> Material {
    Material::new("air", 1.225, 0.0, 0.0)
}
/// Ice (0 degrees C).
pub fn ice() -> Material {
    Material::new("ice", 917.0, 0.03, 0.7)
}
/// Glycerin.
pub fn glycerin() -> Material {
    Material::new("glycerin", 1260.0, 0.0, 0.0)
}
/// Ethanol.
pub fn ethanol() -> Material {
    Material::new("ethanol", 789.0, 0.0, 0.0)
}
/// Mercury: liquid metal, very high density.
pub fn mercury() -> Material {
    Material::new("mercury", 13534.0, 0.0, 0.0)
}
/// Tennis ball.
pub fn tennis_ball() -> Material {
    Material::new("tennis_ball", 580.0, 0.6, 0.75)
}
/// Golf ball.
pub fn golf_ball() -> Material {
    Material::new("golf_ball", 1100.0, 0.35, 0.8)
}
/// Basketball.
pub fn basketball() -> Material {
    Material::new("basketball", 580.0, 0.8, 0.72)
}
/// Billiard ball (phenolic resin).
pub fn billiard_ball() -> Material {
    Material::new("billiard_ball", 1550.0, 0.2, 0.95)
}
/// Bowling ball (polyester).
pub fn bowling_ball() -> Material {
    Material::new("bowling_ball", 1600.0, 0.3, 0.65)
}
/// Extended preset: structural steel (S355).
pub fn extended_steel() -> ExtendedMaterial {
    ExtendedMaterial::new(
        steel(),
        200.0e9,
        0.3,
        355.0e6,
        510.0e6,
        50.0,
        500.0,
        12.0e-6,
    )
}
/// Extended preset: aluminum 6061-T6.
pub fn extended_aluminum() -> ExtendedMaterial {
    ExtendedMaterial::new(
        aluminum(),
        69.0e9,
        0.33,
        275.0e6,
        310.0e6,
        167.0,
        896.0,
        23.6e-6,
    )
}
/// Extended preset: HDPE polymer.
pub fn extended_hdpe() -> ExtendedMaterial {
    ExtendedMaterial::new(hdpe(), 1.0e9, 0.44, 25.0e6, 32.0e6, 0.46, 1900.0, 150.0e-6)
}
/// Extended preset: concrete C30.
pub fn extended_concrete() -> ExtendedMaterial {
    ExtendedMaterial::new(concrete(), 30.0e9, 0.2, 30.0e6, 3.0e6, 1.7, 880.0, 10.0e-6)
}
/// Extended preset: copper.
pub fn extended_copper() -> ExtendedMaterial {
    ExtendedMaterial::new(
        copper(),
        120.0e9,
        0.34,
        70.0e6,
        220.0e6,
        401.0,
        385.0,
        17.0e-6,
    )
}
/// Extended preset: titanium Ti-6Al-4V.
pub fn extended_titanium() -> ExtendedMaterial {
    ExtendedMaterial::new(
        titanium(),
        114.0e9,
        0.33,
        880.0e6,
        950.0e6,
        6.7,
        526.0,
        8.6e-6,
    )
}
/// Extended preset: alumina (Al2O3) ceramic.
pub fn extended_alumina() -> ExtendedMaterial {
    ExtendedMaterial::new(
        alumina(),
        380.0e9,
        0.22,
        300.0e6,
        300.0e6,
        30.0,
        880.0,
        8.0e-6,
    )
}
/// Extended preset: silicon carbide (SiC) ceramic.
pub fn extended_silicon_carbide() -> ExtendedMaterial {
    ExtendedMaterial::new(
        silicon_carbide(),
        410.0e9,
        0.14,
        400.0e6,
        400.0e6,
        120.0,
        750.0,
        4.0e-6,
    )
}
/// Extended preset: CFRP composite (quasi-isotropic).
pub fn extended_cfrp() -> ExtendedMaterial {
    ExtendedMaterial::new(cfrp(), 70.0e9, 0.3, 600.0e6, 800.0e6, 5.0, 1000.0, 2.0e-6)
}
/// Extended preset: PEEK polymer.
pub fn extended_peek() -> ExtendedMaterial {
    ExtendedMaterial::new(peek(), 3.6e9, 0.4, 100.0e6, 100.0e6, 0.25, 2180.0, 47.0e-6)
}
/// Extended preset: cortical bone.
pub fn extended_cortical_bone() -> ExtendedMaterial {
    ExtendedMaterial::new(
        cortical_bone(),
        17.0e9,
        0.3,
        130.0e6,
        150.0e6,
        0.32,
        1260.0,
        11.0e-6,
    )
}
/// Get all metal presets.
pub fn metal_presets() -> Vec<Material> {
    vec![
        steel(),
        stainless_steel(),
        cast_iron(),
        aluminum(),
        copper(),
        brass(),
        titanium(),
        lead(),
        tungsten(),
        magnesium(),
        nickel(),
        zinc(),
        silver(),
        gold(),
        inconel(),
    ]
}
/// Get all ceramic presets.
pub fn ceramic_presets() -> Vec<Material> {
    vec![
        alumina(),
        silicon_carbide(),
        silicon_nitride(),
        zirconia(),
        boron_carbide(),
        porcelain(),
        fused_silica(),
    ]
}
/// Get all composite presets.
pub fn composite_presets() -> Vec<Material> {
    vec![
        cfrp(),
        gfrp(),
        kevlar_composite(),
        ceramic_matrix_composite(),
        metal_matrix_composite(),
        plywood(),
    ]
}
/// Get all polymer presets.
pub fn polymer_presets() -> Vec<Material> {
    vec![
        rubber(),
        soft_rubber(),
        hdpe(),
        ptfe(),
        nylon(),
        polycarbonate(),
        pvc(),
        acrylic(),
        pet(),
        polypropylene(),
        polystyrene(),
        polyurethane_foam(),
        epoxy(),
        peek(),
    ]
}
/// Get all biomaterial presets.
pub fn biomaterial_presets() -> Vec<Material> {
    vec![
        cortical_bone(),
        cancellous_bone(),
        cartilage(),
        silicone_medical(),
        hydroxyapatite(),
        dental_enamel(),
        collagen(),
        bioglass(),
        uhmwpe(),
    ]
}
/// Get all building material presets.
pub fn building_presets() -> Vec<Material> {
    vec![
        wood(),
        oak(),
        glass(),
        tempered_glass(),
        concrete(),
        reinforced_concrete(),
        brick(),
        granite(),
        marble(),
        asphalt(),
        dry_sand(),
        wet_clay(),
    ]
}
/// Get all fluid presets.
pub fn fluid_presets() -> Vec<Material> {
    vec![
        water(),
        seawater(),
        air(),
        ice(),
        glycerin(),
        ethanol(),
        mercury(),
    ]
}
/// Get all sports material presets.
pub fn sports_presets() -> Vec<Material> {
    vec![
        tennis_ball(),
        golf_ball(),
        basketball(),
        billiard_ball(),
        bowling_ball(),
    ]
}
/// Get a list of all standard preset materials.
pub fn all_presets() -> Vec<Material> {
    let mut all = Vec::new();
    all.extend(metal_presets());
    all.extend(ceramic_presets());
    all.extend(composite_presets());
    all.extend(polymer_presets());
    all.extend(biomaterial_presets());
    all.extend(building_presets());
    all.extend(fluid_presets());
    all.extend(sports_presets());
    all
}
/// Get all presets in a given category.
pub fn presets_by_category(category: MaterialCategory) -> Vec<Material> {
    match category {
        MaterialCategory::Metal => metal_presets(),
        MaterialCategory::Polymer => polymer_presets(),
        MaterialCategory::Ceramic => ceramic_presets(),
        MaterialCategory::Composite => composite_presets(),
        MaterialCategory::Biomaterial => biomaterial_presets(),
        MaterialCategory::Building => building_presets(),
        MaterialCategory::Fluid => fluid_presets(),
        MaterialCategory::Sports => sports_presets(),
    }
}
/// Look up a preset material by name (case-insensitive).
///
/// Returns `None` if no match is found.
pub fn preset_by_name(name: &str) -> Option<Material> {
    let lower = name.to_lowercase();
    all_presets()
        .into_iter()
        .find(|m| m.name.to_lowercase() == lower)
}
/// Get all extended material presets.
pub fn all_extended_presets() -> Vec<ExtendedMaterial> {
    vec![
        extended_steel(),
        extended_aluminum(),
        extended_hdpe(),
        extended_concrete(),
        extended_copper(),
        extended_titanium(),
        extended_alumina(),
        extended_silicon_carbide(),
        extended_cfrp(),
        extended_peek(),
        extended_cortical_bone(),
    ]
}
/// Shape memory alloy Nitinol (NiTi) – austenite phase.
pub fn nitinol() -> Material {
    Material::new("nitinol_niti", 6450.0, 0.5, 0.35)
}
/// Hastelloy C-276: corrosion-resistant nickel alloy.
pub fn hastelloy() -> Material {
    Material::new("hastelloy_c276", 8890.0, 0.5, 0.25)
}
/// Maraging steel (18Ni-300): ultra-high strength.
pub fn maraging_steel() -> Material {
    Material::new("maraging_steel_300", 8000.0, 0.58, 0.25)
}
/// Hydrogel (crosslinked polymer in water): very low density, very low restitution.
pub fn hydrogel() -> Material {
    Material::new("hydrogel", 1050.0, 0.05, 0.05)
}
/// Polyimide (PI): high-temperature polymer.
pub fn polyimide() -> Material {
    Material::new("polyimide", 1430.0, 0.4, 0.4)
}
/// Carbon fiber composite (aerospace grade, high-modulus).
pub fn carbon_fiber_composite() -> Material {
    Material::new("carbon_fiber_composite", 1550.0, 0.28, 0.38)
}
/// Glass fiber epoxy composite (hand lay-up).
pub fn glass_fiber_epoxy() -> Material {
    Material::new("glass_fiber_epoxy", 1900.0, 0.35, 0.40)
}
/// Trabecular bone (porous cancellous bone, medium density).
pub fn trabecular_bone() -> Material {
    Material::new("trabecular_bone", 500.0, 0.3, 0.15)
}
/// PEEK biomedical grade (implants).
pub fn peek_biomedical() -> Material {
    Material::new("peek_biomedical", 1310.0, 0.38, 0.42)
}
/// Extended preset: Nitinol (NiTi) shape memory alloy.
pub fn extended_nitinol() -> ExtendedMaterial {
    ExtendedMaterial::new(
        nitinol(),
        83.0e9,
        0.33,
        195.0e6,
        895.0e6,
        18.0,
        490.0,
        11.0e-6,
    )
}
/// Extended preset: PTFE (Teflon).
pub fn extended_ptfe() -> ExtendedMaterial {
    ExtendedMaterial::new(ptfe(), 0.5e9, 0.46, 23.0e6, 31.0e6, 0.25, 1000.0, 135.0e-6)
}
/// Extended preset: polycarbonate.
pub fn extended_polycarbonate() -> ExtendedMaterial {
    ExtendedMaterial::new(
        polycarbonate(),
        2.4e9,
        0.37,
        60.0e6,
        65.0e6,
        0.2,
        1200.0,
        68.0e-6,
    )
}
/// Extended preset: carbon fiber composite (high modulus).
pub fn extended_carbon_fiber_composite() -> ExtendedMaterial {
    ExtendedMaterial::new(
        carbon_fiber_composite(),
        120.0e9,
        0.27,
        800.0e6,
        1000.0e6,
        7.0,
        900.0,
        1.5e-6,
    )
}
/// Extended preset: glass fiber epoxy.
pub fn extended_glass_fiber_epoxy() -> ExtendedMaterial {
    ExtendedMaterial::new(
        glass_fiber_epoxy(),
        25.0e9,
        0.27,
        350.0e6,
        450.0e6,
        0.35,
        1000.0,
        15.0e-6,
    )
}
/// Extended preset: hydrogel.
pub fn extended_hydrogel() -> ExtendedMaterial {
    ExtendedMaterial::new(
        hydrogel(),
        0.001e9,
        0.49,
        0.01e6,
        0.02e6,
        0.6,
        4000.0,
        300.0e-6,
    )
}
/// Extended preset: Nitinol shape memory alloy (extended).
pub fn extended_niti() -> ExtendedMaterial {
    extended_nitinol()
}
/// Extended preset: cortical bone (updated).
pub fn extended_cortical_bone_v2() -> ExtendedMaterial {
    ExtendedMaterial::new(
        cortical_bone(),
        18.0e9,
        0.30,
        150.0e6,
        170.0e6,
        0.33,
        1300.0,
        10.0e-6,
    )
}
/// Compare two materials by density (lighter first).
#[allow(dead_code)]
pub fn compare_by_density(a: &Material, b: &Material) -> std::cmp::Ordering {
    a.density
        .partial_cmp(&b.density)
        .unwrap_or(std::cmp::Ordering::Equal)
}
/// Find the material with the highest Young's modulus from a list.
#[allow(dead_code)]
pub fn stiffest_material(materials: &[ExtendedMaterial]) -> Option<&ExtendedMaterial> {
    materials.iter().max_by(|a, b| {
        a.young_modulus
            .partial_cmp(&b.young_modulus)
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}
/// Find the material with the highest specific strength (σ_y / ρ).
#[allow(dead_code)]
pub fn highest_specific_strength(materials: &[ExtendedMaterial]) -> Option<&ExtendedMaterial> {
    materials.iter().max_by(|a, b| {
        a.specific_strength()
            .partial_cmp(&b.specific_strength())
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}
/// Get all new advanced presets as extended materials.
pub fn advanced_extended_presets() -> Vec<ExtendedMaterial> {
    vec![
        extended_nitinol(),
        extended_ptfe(),
        extended_polycarbonate(),
        extended_carbon_fiber_composite(),
        extended_glass_fiber_epoxy(),
        extended_hydrogel(),
    ]
}
/// Aluminum 7075-T6: high-strength aerospace alloy.
///
/// Used in aircraft structures; higher strength than 6061, slightly lower ductility.
pub fn al_7075_t6() -> Material {
    Material::new("al_7075_t6", 2810.0, 0.47, 0.28)
}
/// Aluminum 2024-T3: fatigue-resistant aerospace alloy.
///
/// Classic aircraft alloy used in fuselage skin and wing structures.
pub fn al_2024_t3() -> Material {
    Material::new("al_2024_t3", 2780.0, 0.47, 0.3)
}
/// Commercially pure titanium Grade 2 (cp-Ti): biocompatible, corrosion resistant.
///
/// Lower strength than Ti-6Al-4V but excellent formability and biocompatibility.
pub fn ti_cp_grade2() -> Material {
    Material::new("ti_cp_grade2", 4510.0, 0.48, 0.3)
}
/// Inconel 625: nickel-chromium-molybdenum superalloy for extreme environments.
///
/// Used in gas turbines, marine applications, and chemical processing.
pub fn inconel_625() -> Material {
    Material::new("inconel_625", 8440.0, 0.5, 0.25)
}
/// Cobalt-chromium alloy (CoCr): high wear resistance, used in dental and orthopaedic implants.
pub fn cobalt_chromium() -> Material {
    Material::new("cobalt_chromium", 8300.0, 0.45, 0.22)
}
/// Beryllium: extremely lightweight structural metal with very high specific stiffness.
///
/// Density 1850 kg/m³, used in aerospace mirrors and precision instruments.
pub fn beryllium() -> Material {
    Material::new("beryllium", 1850.0, 0.4, 0.25)
}
/// Molybdenum: refractory metal for high-temperature structural applications.
pub fn molybdenum() -> Material {
    Material::new("molybdenum", 10220.0, 0.45, 0.2)
}
/// Rhenium: ultra-high melting point metal used in turbine alloys.
pub fn rhenium() -> Material {
    Material::new("rhenium", 21020.0, 0.4, 0.15)
}
/// Extended preset: Aluminum 7075-T6 aerospace alloy.
pub fn extended_al_7075_t6() -> ExtendedMaterial {
    ExtendedMaterial::new(
        al_7075_t6(),
        71.7e9,
        0.33,
        503.0e6,
        572.0e6,
        130.0,
        960.0,
        23.6e-6,
    )
}
/// Extended preset: Aluminum 2024-T3 aerospace alloy.
pub fn extended_al_2024_t3() -> ExtendedMaterial {
    ExtendedMaterial::new(
        al_2024_t3(),
        73.1e9,
        0.33,
        345.0e6,
        483.0e6,
        121.0,
        875.0,
        23.2e-6,
    )
}
/// Extended preset: Inconel 625.
pub fn extended_inconel_625() -> ExtendedMaterial {
    ExtendedMaterial::new(
        inconel_625(),
        207.0e9,
        0.31,
        517.0e6,
        930.0e6,
        12.8,
        410.0,
        13.0e-6,
    )
}
/// Extended preset: Cobalt-Chromium alloy.
pub fn extended_cobalt_chromium() -> ExtendedMaterial {
    ExtendedMaterial::new(
        cobalt_chromium(),
        210.0e9,
        0.3,
        500.0e6,
        900.0e6,
        14.7,
        450.0,
        12.0e-6,
    )
}
/// Extended preset: Beryllium (aerospace grade).
pub fn extended_beryllium() -> ExtendedMaterial {
    ExtendedMaterial::new(
        beryllium(),
        287.0e9,
        0.07,
        207.0e6,
        380.0e6,
        210.0,
        1825.0,
        11.4e-6,
    )
}
/// Get all aerospace alloy presets.
pub fn aerospace_presets() -> Vec<Material> {
    vec![
        al_7075_t6(),
        al_2024_t3(),
        ti_cp_grade2(),
        inconel_625(),
        cobalt_chromium(),
        beryllium(),
        molybdenum(),
        rhenium(),
    ]
}
/// Intervertebral disc (nucleus pulposus): nearly incompressible, viscoelastic.
pub fn intervertebral_disc() -> Material {
    Material::new("intervertebral_disc", 1020.0, 0.02, 0.3)
}
/// Articular cartilage: thin layer covering joint surfaces, very low friction.
pub fn articular_cartilage() -> Material {
    Material::new("articular_cartilage", 1100.0, 0.005, 0.35)
}
/// Meniscus cartilage: fibrocartilage in knee joint.
pub fn meniscus_cartilage() -> Material {
    Material::new("meniscus_cartilage", 1200.0, 0.01, 0.3)
}
/// Ligament: dense connective tissue connecting bone to bone.
pub fn ligament() -> Material {
    Material::new("ligament", 1200.0, 0.3, 0.2)
}
/// Tendon: fibrous connective tissue connecting muscle to bone.
pub fn tendon() -> Material {
    Material::new("tendon", 1165.0, 0.25, 0.15)
}
/// Skin (dermis layer): viscoelastic biological tissue.
pub fn skin_dermis() -> Material {
    Material::new("skin_dermis", 1090.0, 0.5, 0.2)
}
/// PLA (polylactic acid) biopolymer: biodegradable polymer for implants and 3-D printing.
pub fn pla_biopolymer() -> Material {
    Material::new("pla_biopolymer", 1240.0, 0.45, 0.4)
}
/// Extended preset: cortical bone with orthotropic representative values.
pub fn extended_cortical_bone_orthotropic() -> ExtendedMaterial {
    ExtendedMaterial::new(
        cortical_bone(),
        20.0e9,
        0.3,
        140.0e6,
        160.0e6,
        0.34,
        1260.0,
        11.0e-6,
    )
}
/// PTFE (Teflon) enhanced grade: lowest friction engineering polymer.
pub fn ptfe_enhanced() -> Material {
    Material::new("ptfe_enhanced", 2180.0, 0.03, 0.3)
}
/// PEI (Ultem): high-strength, high-temperature amorphous thermoplastic.
pub fn pei_ultem() -> Material {
    Material::new("pei_ultem", 1270.0, 0.42, 0.4)
}
/// PVDF (polyvinylidene fluoride): piezoelectric polymer with good chemical resistance.
pub fn pvdf() -> Material {
    Material::new("pvdf", 1780.0, 0.4, 0.38)
}
/// ABS (acrylonitrile butadiene styrene): common engineering thermoplastic.
pub fn abs_plastic() -> Material {
    Material::new("abs", 1050.0, 0.5, 0.45)
}
/// PA12 (polyamide 12): flexible nylon variant used in SLS 3-D printing.
pub fn pa12_nylon() -> Material {
    Material::new("pa12_nylon", 1010.0, 0.38, 0.42)
}
/// Extended preset: PEEK-CF30 (30% carbon fibre filled PEEK).
pub fn extended_peek_cf30() -> ExtendedMaterial {
    ExtendedMaterial::new(peek(), 15.0e9, 0.38, 220.0e6, 230.0e6, 0.9, 1100.0, 20.0e-6)
}
/// Mullite (3Al₂O₃·2SiO₂): refractory ceramic, excellent thermal stability.
pub fn mullite() -> Material {
    Material::new("mullite", 2800.0, 0.3, 0.4)
}
/// Cordierite (Mg₂Al₄Si₅O₁₈): very low thermal expansion ceramic.
pub fn cordierite() -> Material {
    Material::new("cordierite", 2510.0, 0.3, 0.35)
}
/// Yttria-partially-stabilised zirconia (Y-PSZ): tough ceramic used in thermal barrier coatings.
pub fn zirconia_y_psz() -> Material {
    Material::new("zirconia_y_psz", 5900.0, 0.35, 0.38)
}
/// Extended preset: Y-PSZ zirconia.
pub fn extended_zirconia_y_psz() -> ExtendedMaterial {
    ExtendedMaterial::new(
        zirconia_y_psz(),
        200.0e9,
        0.31,
        900.0e6,
        1100.0e6,
        2.0,
        460.0,
        10.5e-6,
    )
}
/// Concrete grade C20 (characteristic strength 20 MPa).
pub fn concrete_c20() -> Material {
    Material::new("concrete_c20", 2350.0, 0.65, 0.15)
}
/// Concrete grade C30 (characteristic strength 30 MPa).
pub fn concrete_c30() -> Material {
    Material::new("concrete_c30", 2380.0, 0.65, 0.18)
}
/// Concrete grade C40 (characteristic strength 40 MPa).
pub fn concrete_c40() -> Material {
    Material::new("concrete_c40", 2400.0, 0.65, 0.2)
}
/// Concrete grade C50 (characteristic strength 50 MPa).
pub fn concrete_c50() -> Material {
    Material::new("concrete_c50", 2420.0, 0.65, 0.2)
}
/// Concrete grade C60 (high-performance concrete, 60 MPa).
pub fn concrete_c60() -> Material {
    Material::new("concrete_c60", 2450.0, 0.65, 0.22)
}
/// Extended preset: Concrete C20.
pub fn extended_concrete_c20() -> ExtendedMaterial {
    ExtendedMaterial::new(
        concrete_c20(),
        28.0e9,
        0.2,
        20.0e6,
        2.2e6,
        1.5,
        850.0,
        10.0e-6,
    )
}
/// Extended preset: Concrete C40.
pub fn extended_concrete_c40() -> ExtendedMaterial {
    ExtendedMaterial::new(
        concrete_c40(),
        35.0e9,
        0.2,
        40.0e6,
        3.8e6,
        1.7,
        870.0,
        10.0e-6,
    )
}
/// Extended preset: Concrete C60 high-performance.
pub fn extended_concrete_c60() -> ExtendedMaterial {
    ExtendedMaterial::new(
        concrete_c60(),
        42.0e9,
        0.2,
        60.0e6,
        5.0e6,
        1.9,
        900.0,
        10.0e-6,
    )
}
/// Gravel: coarse granular soil with high internal friction angle.
pub fn gravel() -> Material {
    Material::new("gravel", 1700.0, 0.7, 0.08)
}
/// Silt: fine-grained soil between sand and clay in particle size.
pub fn silt() -> Material {
    Material::new("silt", 1850.0, 0.25, 0.05)
}
/// Dry clay: fine-grained plastic soil in dry state.
pub fn dry_clay() -> Material {
    Material::new("dry_clay", 1650.0, 0.35, 0.06)
}
/// Saturated clay: clay with high water content, very cohesive.
pub fn saturated_clay() -> Material {
    Material::new("saturated_clay", 2000.0, 0.15, 0.04)
}
/// Dense sand: compacted granular soil.
pub fn dense_sand() -> Material {
    Material::new("dense_sand", 1900.0, 0.65, 0.1)
}
/// Loose sand: poorly compacted granular soil.
pub fn loose_sand() -> Material {
    Material::new("loose_sand", 1400.0, 0.5, 0.08)
}
/// Peat: organic soil with high compressibility.
pub fn peat() -> Material {
    Material::new("peat", 800.0, 0.2, 0.05)
}
/// Rammed earth: compacted earth used in construction.
pub fn rammed_earth() -> Material {
    Material::new("rammed_earth", 2000.0, 0.55, 0.12)
}
/// Get all soil presets.
pub fn soil_presets() -> Vec<Material> {
    vec![
        dry_sand(),
        wet_clay(),
        gravel(),
        silt(),
        dry_clay(),
        saturated_clay(),
        dense_sand(),
        loose_sand(),
        peat(),
        rammed_earth(),
    ]
}
/// Pine wood (Scots pine): common softwood structural material.
pub fn pine_wood() -> Material {
    Material::new("pine_wood", 530.0, 0.5, 0.38)
}
/// Balsa wood: extremely lightweight wood used in model aircraft and sandwich panels.
pub fn balsa_wood() -> Material {
    Material::new("balsa_wood", 160.0, 0.45, 0.35)
}
/// Beech wood: European hardwood used in furniture and flooring.
pub fn beech_wood() -> Material {
    Material::new("beech_wood", 720.0, 0.55, 0.36)
}
/// Teak wood: dense tropical hardwood with natural oil resistance.
pub fn teak_wood() -> Material {
    Material::new("teak_wood", 800.0, 0.5, 0.32)
}
/// Maple wood: hard, dense North American hardwood.
pub fn maple_wood() -> Material {
    Material::new("maple_wood", 705.0, 0.55, 0.35)
}
/// Walnut wood: premium hardwood for furniture and gunstocks.
pub fn walnut_wood() -> Material {
    Material::new("walnut_wood", 640.0, 0.5, 0.33)
}
/// Bamboo: rapidly renewable grass with mechanical properties rivalling hardwood.
pub fn bamboo() -> Material {
    Material::new("bamboo", 700.0, 0.45, 0.36)
}
/// MDF board (medium-density fibreboard): engineered wood product.
pub fn mdf_board() -> Material {
    Material::new("mdf_board", 750.0, 0.5, 0.3)
}
/// Oriented strand board (OSB): structural engineered wood panel.
pub fn osb_board() -> Material {
    Material::new("osb_board", 640.0, 0.5, 0.28)
}
/// Extended preset: Oak (longitudinal direction).
pub fn extended_oak() -> ExtendedMaterial {
    ExtendedMaterial::new(oak(), 12.0e9, 0.3, 55.0e6, 90.0e6, 0.17, 1700.0, 5.0e-6)
}
/// Extended preset: Pine wood (longitudinal direction).
pub fn extended_pine_wood() -> ExtendedMaterial {
    ExtendedMaterial::new(
        pine_wood(),
        10.0e9,
        0.3,
        40.0e6,
        70.0e6,
        0.14,
        1600.0,
        4.5e-6,
    )
}
