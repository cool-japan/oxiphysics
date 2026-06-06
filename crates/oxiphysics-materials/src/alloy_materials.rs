// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Metallic alloy material models: composition, mechanical properties,
//! strengthening mechanisms, phase diagrams, thermal and corrosion properties.

// ─── AlloyComposition ────────────────────────────────────────────────────────

/// Chemical composition of a metallic alloy.
#[derive(Debug, Clone)]
pub struct AlloyComposition {
    /// List of (element symbol, weight fraction) pairs.
    pub elements: Vec<(String, f64)>,
    /// Bulk density in kg/m³.
    pub density: f64,
    /// (solidus, liquidus) melting range in K.
    pub melting_range: (f64, f64),
}

impl AlloyComposition {
    /// Construct a new `AlloyComposition`.
    pub fn new(elements: Vec<(String, f64)>, density: f64, melting_range: (f64, f64)) -> Self {
        Self {
            elements,
            density,
            melting_range,
        }
    }

    /// Return `true` if all weight fractions sum to approximately 1.
    pub fn validate(&self) -> bool {
        let sum: f64 = self.elements.iter().map(|(_, f)| f).sum();
        (sum - 1.0).abs() < 1e-6
    }

    /// Compute mixture density using the rule of mixtures (weight-fraction average
    /// of inverse densities — Reuss bound).
    ///
    /// `densities` is a list of `(symbol, density_kg_m3)` pairs.
    pub fn density_mixture(&self, densities: &[(String, f64)]) -> f64 {
        let mut inv_sum = 0.0;
        let mut frac_sum = 0.0;
        for (sym, frac) in &self.elements {
            if let Some((_, d)) = densities.iter().find(|(s, _)| s == sym)
                && *d > 0.0
            {
                inv_sum += frac / d;
                frac_sum += frac;
            }
        }
        if inv_sum < 1e-30 || frac_sum < 1e-12 {
            self.density
        } else {
            frac_sum / inv_sum
        }
    }
}

// ─── AlloySeries ─────────────────────────────────────────────────────────────

/// Classification of common alloy families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlloySeries {
    /// Aluminium 1xxx series (≥99% Al).
    Series1xxx,
    /// Aluminium 2xxx series (Cu-based).
    Series2xxx,
    /// Aluminium 3xxx series (Mn-based).
    Series3xxx,
    /// Aluminium 4xxx series (Si-based).
    Series4xxx,
    /// Aluminium 5xxx series (Mg-based).
    Series5xxx,
    /// Aluminium 6xxx series (Mg+Si).
    Series6xxx,
    /// Aluminium 7xxx series (Zn-based).
    Series7xxx,
    /// Aluminium 8xxx series (other elements).
    Series8xxx,
    /// Aluminium 9xxx series (reserved).
    Series9xxx,
    /// Austenitic stainless steel 304.
    StainlessSteel304,
    /// Austenitic stainless steel 316.
    StainlessSteel316,
    /// Nickel superalloy Inconel 718.
    Inconel718,
    /// Titanium alloy Ti-6Al-4V.
    TitaniumTi6Al4V,
    /// Alloy steel AISI 4140.
    Tool4140,
    /// Maraging steel 300 grade.
    Maraging300,
    /// Nickel-based alloy Hastelloy C-276.
    Hastelloy,
}

// ─── AlloyMechanicalProps ─────────────────────────────────────────────────────

/// Mechanical properties of a metallic alloy.
#[derive(Debug, Clone)]
pub struct AlloyMechanicalProps {
    /// 0.2% proof stress / yield strength in MPa.
    pub yield_strength: f64,
    /// Ultimate tensile strength in MPa.
    pub uts: f64,
    /// Elongation at fracture in % (gauge length = 50 mm).
    pub elongation: f64,
    /// Vickers hardness HV.
    pub hardness_hv: f64,
    /// Young's modulus in GPa.
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Plane-strain fracture toughness K_IC in MPa·√m.
    pub fracture_toughness: f64,
}

impl AlloyMechanicalProps {
    /// Construct `AlloyMechanicalProps` from individual values.
    pub fn new(
        yield_strength: f64,
        uts: f64,
        elongation: f64,
        hardness_hv: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
        fracture_toughness: f64,
    ) -> Self {
        Self {
            yield_strength,
            uts,
            elongation,
            hardness_hv,
            youngs_modulus,
            poisson_ratio,
            fracture_toughness,
        }
    }

    /// Compute safety factor = yield_strength / applied_stress.
    pub fn safety_factor(&self, applied_stress: f64) -> f64 {
        if applied_stress.abs() < 1e-15 {
            f64::INFINITY
        } else {
            self.yield_strength / applied_stress
        }
    }

    /// Return `true` if the alloy is considered brittle (elongation < 5 %).
    pub fn is_brittle(&self) -> bool {
        self.elongation < 5.0
    }
}

// ─── Hall-Petch equation ──────────────────────────────────────────────────────

/// Hall-Petch grain-boundary strengthening model.
pub struct HallPetch;

impl HallPetch {
    /// Compute yield strength:  σ_y = σ_0 + k / √(grain_size).
    ///
    /// `grain_size` is in metres; `k` is the Hall-Petch slope in MPa·√m.
    pub fn yield_strength(sigma_0: f64, k: f64, grain_size: f64) -> f64 {
        sigma_0 + k / grain_size.max(1e-30).sqrt()
    }

    /// Invert Hall-Petch to find required grain size for a target yield stress.
    pub fn grain_size_from_yield(sigma_y: f64, sigma_0: f64, k: f64) -> f64 {
        let diff = sigma_y - sigma_0;
        if diff.abs() < 1e-15 {
            return f64::INFINITY;
        }
        (k / diff).powi(2)
    }
}

// ─── Strengthening mechanisms ─────────────────────────────────────────────────

/// Collection of metallic strengthening mechanism models.
pub struct Strengthening;

impl Strengthening {
    /// Solid-solution strengthening increment: Δσ_ss = k_ss · c^(2/3).
    ///
    /// `c` is solute concentration (mol fraction); `k_ss` is a material constant (MPa).
    pub fn solid_solution_strengthening(c: f64, k_ss: f64) -> f64 {
        k_ss * c.max(0.0).powf(2.0 / 3.0)
    }

    /// Orowan precipitation-hardening increment.
    ///
    /// Δσ_ph ≈ 0.81 · G · b / (2π · λ · √(1-ν)) where ν ≈ 0.3 is absorbed.
    ///
    /// Here we use the simplified Orowan–Ashby form:
    /// Δσ ≈ 0.13 · G · b / particle_spacing
    pub fn precipitation_hardening(particle_spacing: f64, shear_modulus: f64, burgers: f64) -> f64 {
        0.13 * shear_modulus * burgers / particle_spacing.max(1e-30)
    }

    /// Hollomon work-hardening (power-law): σ = k · ε_p^n.
    ///
    /// `eps_p` is plastic strain, `k` is strength coefficient (MPa), `n` is strain-hardening exponent.
    pub fn work_hardening(eps_p: f64, k: f64, n: f64) -> f64 {
        k * eps_p.max(0.0).powf(n)
    }

    /// Combine strengthening contributions by simple summation.
    ///
    /// `ss` = solid-solution, `ph` = precipitation, `wh` = work-hardening,
    /// `gb` = grain-boundary (Hall-Petch) increment.
    pub fn combined_strengthening(ss: f64, ph: f64, wh: f64, gb: f64) -> f64 {
        ss + ph + wh + gb
    }
}

// ─── Binary phase diagram ─────────────────────────────────────────────────────

/// Phase label in a binary alloy phase diagram.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinaryPhase {
    /// Fully liquid.
    Liquid,
    /// Single-phase alpha solid solution.
    Alpha,
    /// Single-phase beta solid solution.
    Beta,
    /// Two-phase α + β region.
    AlphaPlusBeta,
    /// Eutectic mixture.
    Eutectic,
}

/// Simple binary phase diagram helper.
pub struct PhaseDiagram;

impl PhaseDiagram {
    /// Determine the equilibrium phase for composition `x` (mol fraction of B)
    /// and temperature `temp` (K) in a simple eutectic binary system.
    ///
    /// `eutectic` = `(x_eutectic, T_eutectic)`.
    /// `liquidus` = `[T_melt_A, T_melt_B]` (melting points of pure A and B).
    pub fn phase_at_composition(
        x: f64,
        temp: f64,
        eutectic: (f64, f64),
        liquidus: [f64; 2],
    ) -> BinaryPhase {
        let (x_e, t_e) = eutectic;
        let t_a = liquidus[0];
        let t_b = liquidus[1];

        // Linear liquidus on A-rich side: T_liq_a(x) = T_A - (T_A - T_e) * x / x_e
        let t_liq_a = if x_e > 1e-12 {
            t_a - (t_a - t_e) * x / x_e
        } else {
            t_a
        };

        // Linear liquidus on B-rich side: T_liq_b(x) = T_B - (T_B - T_e) * (1-x) / (1-x_e)
        let t_liq_b = if (1.0 - x_e) > 1e-12 {
            t_b - (t_b - t_e) * (1.0 - x) / (1.0 - x_e)
        } else {
            t_b
        };

        let t_liquidus = if x <= x_e { t_liq_a } else { t_liq_b };

        if temp > t_liquidus {
            BinaryPhase::Liquid
        } else if (temp - t_e).abs() < 5.0 && (x - x_e).abs() < 0.02 {
            BinaryPhase::Eutectic
        } else if temp < t_e {
            BinaryPhase::AlphaPlusBeta
        } else if x < x_e {
            BinaryPhase::Alpha
        } else {
            BinaryPhase::Beta
        }
    }
}

// ─── Thermal properties ───────────────────────────────────────────────────────

/// Compute alloy thermal conductivity by linear rule of mixtures.
///
/// `k1`, `k2` are the component conductivities (W/m·K); `x` is the mole fraction of component 2.
pub fn alloy_thermal_conductivity(k1: f64, k2: f64, x: f64) -> f64 {
    (1.0 - x) * k1 + x * k2
}

/// Compute alloy specific heat capacity by mass-weighted mixture rule.
///
/// Each tuple is `(mass_fraction, cp)` in J/(kg·K).
pub fn alloy_specific_heat(cp_components: &[(f64, f64)]) -> f64 {
    cp_components.iter().map(|(f, cp)| f * cp).sum()
}

// ─── Corrosion ────────────────────────────────────────────────────────────────

/// Compute the Pitting Resistance Equivalent Number (PREN) for stainless steels.
///
/// PREN = %Cr + 3.3·%Mo + 16·%N
pub fn pitting_resistance_equivalent(cr: f64, mo: f64, n_pct: f64) -> f64 {
    cr + 3.3 * mo + 16.0 * n_pct
}

/// Assess galvanic corrosion risk from the electrode-potential difference.
///
/// Returns a static string describing the risk level.
pub fn galvanic_corrosion_risk(e1: f64, e2: f64) -> &'static str {
    let delta = (e1 - e2).abs();
    if delta < 0.1 {
        "negligible"
    } else if delta < 0.25 {
        "low"
    } else if delta < 0.5 {
        "moderate"
    } else {
        "high"
    }
}

// ─── AluminumAlloyT6 ──────────────────────────────────────────────────────────

/// Typical T6-temper mechanical properties for common aluminium alloys.
pub struct AluminumAlloyT6;

impl AluminumAlloyT6 {
    /// Return typical T6 mechanical properties for a given aluminium alloy series.
    pub fn properties(series: AlloySeries) -> AlloyMechanicalProps {
        match series {
            AlloySeries::Series2xxx => {
                AlloyMechanicalProps::new(324.0, 469.0, 10.0, 130.0, 73.1, 0.33, 37.0)
            }
            AlloySeries::Series6xxx => {
                AlloyMechanicalProps::new(276.0, 310.0, 12.0, 95.0, 68.9, 0.33, 29.0)
            }
            AlloySeries::Series7xxx => {
                AlloyMechanicalProps::new(503.0, 572.0, 11.0, 150.0, 71.7, 0.33, 29.0)
            }
            _ => AlloyMechanicalProps::new(100.0, 150.0, 8.0, 50.0, 69.0, 0.33, 20.0),
        }
    }
}

// ─── NickelSuperalloy ─────────────────────────────────────────────────────────

/// Creep and oxidation model for nickel-based superalloys.
pub struct NickelSuperalloy;

impl NickelSuperalloy {
    /// Compute steady-state creep rate using the Norton power law:
    /// ε̇ = A · σ^n · exp(-Q / (R·T))
    ///
    /// - `sigma`: stress in MPa
    /// - `temp`: temperature in K
    /// - `a_coeff`: pre-exponential factor (s⁻¹·MPa⁻ⁿ)
    /// - `n_exp`: stress exponent
    /// - `q_activation`: activation energy in J/mol
    pub fn creep_rate(sigma: f64, temp: f64, a_coeff: f64, n_exp: f64, q_activation: f64) -> f64 {
        const R: f64 = 8.314; // J/(mol·K)
        a_coeff * sigma.powf(n_exp) * (-q_activation / (R * temp.max(1e-3))).exp()
    }

    /// Estimate parabolic oxidation mass gain: Δm² = k_p · t.
    ///
    /// Returns Δm (mg/cm²) at time `t` (s) using parabolic rate constant `k_p`.
    pub fn oxidation_mass_gain(k_p: f64, t: f64) -> f64 {
        (k_p * t.max(0.0)).sqrt()
    }
}

// ─── WeldabilityIndex ─────────────────────────────────────────────────────────

/// Carbon equivalent for steel weldability assessment.
pub struct WeldabilityIndex;

impl WeldabilityIndex {
    /// Compute the International Institute of Welding (IIW) carbon equivalent:
    ///
    /// CE = C + Mn/6 + (Cr+Mo+V)/5 + (Ni+Cu)/15
    ///
    /// All values in weight %.
    pub fn carbon_equivalent_iiw(
        c: f64,
        mn: f64,
        cr: f64,
        mo: f64,
        v: f64,
        ni: f64,
        cu: f64,
    ) -> f64 {
        c + mn / 6.0 + (cr + mo + v) / 5.0 + (ni + cu) / 15.0
    }

    /// Classify weldability from the IIW carbon equivalent.
    pub fn weldability_class(ce: f64) -> &'static str {
        if ce < 0.35 {
            "excellent"
        } else if ce < 0.45 {
            "good"
        } else if ce < 0.60 {
            "fair"
        } else {
            "poor"
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn stainless_304() -> AlloyComposition {
        AlloyComposition::new(
            vec![
                ("Fe".to_string(), 0.69),
                ("Cr".to_string(), 0.19),
                ("Ni".to_string(), 0.10),
                ("Mn".to_string(), 0.02),
            ],
            7900.0,
            (1673.0, 1723.0),
        )
    }

    fn al6061_t6() -> AlloyMechanicalProps {
        AlloyMechanicalProps::new(276.0, 310.0, 12.0, 95.0, 68.9, 0.33, 29.0)
    }

    #[test]
    fn test_alloy_composition_validate_valid() {
        let alloy = stainless_304();
        assert!(alloy.validate());
    }

    #[test]
    fn test_alloy_composition_validate_invalid() {
        let alloy = AlloyComposition::new(
            vec![("Fe".to_string(), 0.5), ("Cr".to_string(), 0.3)],
            7900.0,
            (1673.0, 1723.0),
        );
        assert!(!alloy.validate());
    }

    #[test]
    fn test_density_mixture_rule_of_mixtures() {
        let alloy = AlloyComposition::new(
            vec![("A".to_string(), 0.5), ("B".to_string(), 0.5)],
            0.0,
            (1000.0, 1200.0),
        );
        let densities = vec![("A".to_string(), 2000.0), ("B".to_string(), 4000.0)];
        let d = alloy.density_mixture(&densities);
        // Reuss bound = 1 / (0.5/2000 + 0.5/4000) = 2666.67
        assert!((d - 2666.67).abs() < 1.0, "d = {:.6}", d);
    }

    #[test]
    fn test_density_mixture_missing_element() {
        let alloy = AlloyComposition::new(vec![("X".to_string(), 1.0)], 7000.0, (1000.0, 1100.0));
        let densities = vec![("Fe".to_string(), 7874.0)];
        // Falls back to alloy.density
        let d = alloy.density_mixture(&densities);
        assert!((d - 7000.0).abs() < 1e-6);
    }

    #[test]
    fn test_safety_factor_normal() {
        let props = al6061_t6();
        let sf = props.safety_factor(138.0);
        assert!((sf - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_safety_factor_zero_stress() {
        let props = al6061_t6();
        assert!(props.safety_factor(0.0).is_infinite());
    }

    #[test]
    fn test_is_brittle_ductile() {
        let props = al6061_t6();
        assert!(!props.is_brittle());
    }

    #[test]
    fn test_is_brittle_true() {
        let props = AlloyMechanicalProps::new(600.0, 700.0, 2.0, 700.0, 210.0, 0.28, 50.0);
        assert!(props.is_brittle());
    }

    #[test]
    fn test_hall_petch_yield_strength() {
        // σ_0 = 50 MPa, k = 0.5 MPa·mm^0.5, grain_size = 0.25 mm² → k/√d = 1.0
        let sy = HallPetch::yield_strength(50.0, 0.5, 0.25);
        assert!((sy - 51.0).abs() < 1e-9, "sy = {:.6}", sy);
    }

    #[test]
    fn test_hall_petch_grain_size_inversion() {
        let sigma_0 = 50.0;
        let k = 0.5;
        let d0 = 1.0e-6;
        let sy = HallPetch::yield_strength(sigma_0, k, d0);
        let d_back = HallPetch::grain_size_from_yield(sy, sigma_0, k);
        assert!((d_back - d0).abs() < 1e-18, "d_back = {:.6e}", d_back);
    }

    #[test]
    fn test_hall_petch_grain_size_same_yield() {
        let d = HallPetch::grain_size_from_yield(50.0, 50.0, 0.5);
        assert!(d.is_infinite());
    }

    #[test]
    fn test_solid_solution_strengthening_zero() {
        assert!((Strengthening::solid_solution_strengthening(0.0, 100.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_solid_solution_strengthening_positive() {
        let ds = Strengthening::solid_solution_strengthening(0.01, 500.0);
        assert!(ds > 0.0);
    }

    #[test]
    fn test_precipitation_hardening_positive() {
        let ds = Strengthening::precipitation_hardening(1e-7, 26e3, 2.86e-10);
        assert!(ds > 0.0);
    }

    #[test]
    fn test_work_hardening_hollomon() {
        // σ = 500 · 0.2^0.2
        let expected = 500.0 * 0.2_f64.powf(0.2);
        let result = Strengthening::work_hardening(0.2, 500.0, 0.2);
        assert!((result - expected).abs() < 1e-9);
    }

    #[test]
    fn test_combined_strengthening() {
        let total = Strengthening::combined_strengthening(50.0, 30.0, 20.0, 10.0);
        assert!((total - 110.0).abs() < 1e-10);
    }

    #[test]
    fn test_phase_diagram_liquid() {
        // Above liquidus → Liquid
        let phase = PhaseDiagram::phase_at_composition(0.3, 1500.0, (0.5, 800.0), [1200.0, 1400.0]);
        assert_eq!(phase, BinaryPhase::Liquid);
    }

    #[test]
    fn test_phase_diagram_alpha() {
        let phase = PhaseDiagram::phase_at_composition(0.1, 900.0, (0.5, 600.0), [1200.0, 1100.0]);
        assert_eq!(phase, BinaryPhase::Alpha);
    }

    #[test]
    fn test_phase_diagram_beta() {
        let phase = PhaseDiagram::phase_at_composition(0.8, 900.0, (0.5, 600.0), [1200.0, 1100.0]);
        assert_eq!(phase, BinaryPhase::Beta);
    }

    #[test]
    fn test_phase_diagram_two_phase() {
        let phase = PhaseDiagram::phase_at_composition(0.3, 500.0, (0.5, 600.0), [1200.0, 1100.0]);
        assert_eq!(phase, BinaryPhase::AlphaPlusBeta);
    }

    #[test]
    fn test_alloy_thermal_conductivity_pure_components() {
        let k = alloy_thermal_conductivity(15.0, 400.0, 0.0);
        assert!((k - 15.0).abs() < 1e-10);
        let k2 = alloy_thermal_conductivity(15.0, 400.0, 1.0);
        assert!((k2 - 400.0).abs() < 1e-10);
    }

    #[test]
    fn test_alloy_thermal_conductivity_midpoint() {
        let k = alloy_thermal_conductivity(10.0, 20.0, 0.5);
        assert!((k - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_alloy_specific_heat_single() {
        let cp = alloy_specific_heat(&[(1.0, 500.0)]);
        assert!((cp - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_alloy_specific_heat_mixture() {
        let cp = alloy_specific_heat(&[(0.7, 500.0), (0.3, 900.0)]);
        assert!((cp - 620.0).abs() < 1e-10);
    }

    #[test]
    fn test_pren_calculation() {
        // 316L: Cr≈17, Mo≈2.5, N≈0.03
        let pren = pitting_resistance_equivalent(17.0, 2.5, 0.03);
        let expected = 17.0 + 3.3 * 2.5 + 16.0 * 0.03;
        assert!((pren - expected).abs() < 1e-10);
    }

    #[test]
    fn test_galvanic_corrosion_negligible() {
        assert_eq!(galvanic_corrosion_risk(0.0, 0.05), "negligible");
    }

    #[test]
    fn test_galvanic_corrosion_high() {
        assert_eq!(galvanic_corrosion_risk(0.0, 1.0), "high");
    }

    #[test]
    fn test_galvanic_corrosion_low() {
        assert_eq!(galvanic_corrosion_risk(0.1, 0.3), "low");
    }

    #[test]
    fn test_aluminum_t6_6xxx() {
        let props = AluminumAlloyT6::properties(AlloySeries::Series6xxx);
        assert!((props.yield_strength - 276.0).abs() < 1e-6);
    }

    #[test]
    fn test_aluminum_t6_7xxx_high_strength() {
        let p7 = AluminumAlloyT6::properties(AlloySeries::Series7xxx);
        let p6 = AluminumAlloyT6::properties(AlloySeries::Series6xxx);
        assert!(p7.yield_strength > p6.yield_strength);
    }

    #[test]
    fn test_nickel_creep_rate_positive() {
        let rate = NickelSuperalloy::creep_rate(200.0, 1073.0, 1e-15, 4.0, 290_000.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_nickel_creep_rate_increases_with_stress() {
        let r1 = NickelSuperalloy::creep_rate(100.0, 1073.0, 1e-15, 4.0, 290_000.0);
        let r2 = NickelSuperalloy::creep_rate(200.0, 1073.0, 1e-15, 4.0, 290_000.0);
        assert!(r2 > r1);
    }

    #[test]
    fn test_oxidation_mass_gain_zero_time() {
        let dm = NickelSuperalloy::oxidation_mass_gain(1e-12, 0.0);
        assert!((dm - 0.0).abs() < 1e-20);
    }

    #[test]
    fn test_oxidation_mass_gain_positive() {
        let dm = NickelSuperalloy::oxidation_mass_gain(1e-12, 3600.0);
        assert!(dm > 0.0);
    }

    #[test]
    fn test_carbon_equivalent_low_alloy() {
        let ce = WeldabilityIndex::carbon_equivalent_iiw(0.15, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!((ce - 0.3167).abs() < 0.001, "ce = {:.6}", ce);
    }

    #[test]
    fn test_weldability_class_excellent() {
        assert_eq!(WeldabilityIndex::weldability_class(0.30), "excellent");
    }

    #[test]
    fn test_weldability_class_poor() {
        assert_eq!(WeldabilityIndex::weldability_class(0.65), "poor");
    }

    #[test]
    fn test_alloy_series_debug() {
        let s = format!("{:?}", AlloySeries::Inconel718);
        assert!(s.contains("Inconel718"));
    }

    #[test]
    fn test_binary_phase_debug() {
        let s = format!("{:?}", BinaryPhase::AlphaPlusBeta);
        assert!(s.contains("Alpha"));
    }
}
