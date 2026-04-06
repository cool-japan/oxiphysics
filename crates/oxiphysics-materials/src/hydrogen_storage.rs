// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hydrogen storage material models.
//!
//! Covers compressed gas, liquid hydrogen, metal hydrides, chemical hydrides,
//! and metal-organic frameworks (MOFs).  All functions use SI units unless
//! stated otherwise; pressures are given in bar and temperatures in Kelvin.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

/// Universal gas constant (J / mol·K).
const R_GAS: f64 = 8.314;
/// Molar mass of H₂ (kg / mol).
const M_H2: f64 = 2.016e-3;
/// Lower heating value of H₂ (MJ / kg).
const LHV_H2_MJ_KG: f64 = 119.96;
/// Van der Waals constant *a* for H₂ (Pa·m⁶/mol²).
const VDW_A_H2: f64 = 0.02476e0; // 0.02476 L²·atm/mol² → converted to Pa·m⁶/mol²
/// Van der Waals constant *b* for H₂ (m³/mol).
const VDW_B_H2: f64 = 26.61e-6; // 26.61 cm³/mol

/// Supported hydrogen storage technology types.
#[derive(Clone, Debug, PartialEq)]
pub enum HydrogenStorageType {
    /// Compressed gaseous hydrogen (CGH₂).
    CompressedGas,
    /// Liquid hydrogen (LH₂, cryogenic ~20 K).
    LiquidH2,
    /// Reversible metal hydride (e.g. LaNi₅H₆, MgH₂).
    MetalHydride,
    /// Chemical / irreversible hydride (e.g. NaBH₄, NH₃BH₃).
    ChemicalHydride,
    /// Metal-organic framework physisorption.
    MOF,
}

// ---------------------------------------------------------------------------
// Metal Hydride
// ---------------------------------------------------------------------------

/// Reversible metal hydride model using van 't Hoff thermodynamics.
///
/// Equilibrium plateau pressure is described by
/// `ln(P_eq / P_ref) = ΔH/(R·T) − ΔS/R`
/// where P_ref = 1 bar.
#[derive(Clone, Debug)]
pub struct MetalHydride {
    /// Common name of the host metal or alloy (e.g. `"LaNi5"`).
    pub metal_name: String,
    /// Maximum gravimetric capacity (wt%).
    pub max_capacity_wt_percent: f64,
    /// Desorption onset temperature at 1 bar (K).
    pub desorption_temperature: f64,
    /// Van 't Hoff enthalpy ΔH (J / mol H₂), typically negative for absorption.
    pub van_hoff_enthalpy: f64,
    /// Van 't Hoff entropy ΔS (J / mol·K).
    pub van_hoff_entropy: f64,
}

impl MetalHydride {
    /// Construct a new metal hydride descriptor.
    pub fn new(
        metal_name: impl Into<String>,
        max_capacity_wt_percent: f64,
        desorption_temperature: f64,
        van_hoff_enthalpy: f64,
        van_hoff_entropy: f64,
    ) -> Self {
        Self {
            metal_name: metal_name.into(),
            max_capacity_wt_percent,
            desorption_temperature,
            van_hoff_enthalpy,
            van_hoff_entropy,
        }
    }

    /// LaNi₅ preset: ~1.4 wt%, moderate temperature, well-known system.
    pub fn lani5() -> Self {
        Self::new("LaNi5", 1.4, 310.0, -30_500.0, -108.0)
    }

    /// MgH₂ preset: ~7.6 wt%, high temperature (requires ~300 °C).
    pub fn mgh2() -> Self {
        Self::new("MgH2", 7.6, 573.0, -74_500.0, -135.0)
    }

    /// Equilibrium plateau pressure (bar) at temperature `t_k` (K) using
    /// the van 't Hoff relation.
    ///
    /// `P_eq = exp(ΔH / (R·T) − ΔS / R)` \[bar\]
    pub fn equilibrium_pressure(&self, t_k: f64) -> f64 {
        // ln(P_eq / 1 bar) = ΔH/(R·T) − ΔS/R
        let ln_p = self.van_hoff_enthalpy / (R_GAS * t_k) - self.van_hoff_entropy / R_GAS;
        ln_p.exp()
    }
}

// ---------------------------------------------------------------------------
// Compressed-Gas Pressure Vessel
// ---------------------------------------------------------------------------

/// High-pressure hydrogen vessel (compressed gas storage).
#[derive(Clone, Debug)]
pub struct PressureVesselH2 {
    /// Internal volume (litres).
    pub volume_l: f64,
    /// Maximum allowable pressure (bar).
    pub max_pressure_bar: f64,
    /// Gas temperature (K).
    pub temperature: f64,
}

impl PressureVesselH2 {
    /// Construct a new pressure vessel with the given parameters.
    pub fn new(volume_l: f64, max_pressure_bar: f64, temperature: f64) -> Self {
        Self {
            volume_l,
            max_pressure_bar,
            temperature,
        }
    }

    /// Typical 700-bar Type-IV composite vessel (~120 L, 25 °C).
    pub fn type_iv_700bar() -> Self {
        Self::new(120.0, 700.0, 298.15)
    }

    /// Mass of stored H₂ (kg) using the ideal gas law as approximation.
    ///
    /// m = P·V·M / (R·T), where P is in Pa and V is in m³.
    pub fn mass_h2(&self) -> f64 {
        let p_pa = self.max_pressure_bar * 1e5;
        let v_m3 = self.volume_l * 1e-3;
        p_pa * v_m3 * M_H2 / (R_GAS * self.temperature)
    }

    /// Gravimetric energy density (MJ / kg of H₂) assuming LHV.
    pub fn energy_density(&self) -> f64 {
        // The energy stored per kg of H₂ is simply the lower heating value.
        LHV_H2_MJ_KG
    }

    /// Total stored energy (MJ).
    pub fn stored_energy_mj(&self) -> f64 {
        self.mass_h2() * LHV_H2_MJ_KG
    }
}

// ---------------------------------------------------------------------------
// MOF Adsorption
// ---------------------------------------------------------------------------

/// Metal-organic framework (MOF) hydrogen adsorption model.
///
/// Uses a simplified Langmuir-type isotherm combined with the
/// Clausius-Clapeyron equation for temperature dependence.
#[derive(Clone, Debug)]
pub struct MofAdsorption {
    /// BET surface area (m² / g).
    pub surface_area_m2g: f64,
    /// Total pore volume (cm³ / g).
    pub pore_volume: f64,
    /// Isosteric heat of adsorption |Q_st| (J / mol), positive value.
    pub isosteric_heat: f64,
    /// Saturation capacity at reference conditions (wt%).
    pub saturation_capacity_wt: f64,
    /// Reference temperature for the Langmuir constant (K).
    pub t_ref: f64,
    /// Langmuir constant at reference temperature (1 / bar).
    pub k_lang_ref: f64,
}

impl MofAdsorption {
    /// Construct a new MOF adsorption model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface_area_m2g: f64,
        pore_volume: f64,
        isosteric_heat: f64,
        saturation_capacity_wt: f64,
        t_ref: f64,
        k_lang_ref: f64,
    ) -> Self {
        Self {
            surface_area_m2g,
            pore_volume,
            isosteric_heat,
            saturation_capacity_wt,
            t_ref,
            k_lang_ref,
        }
    }

    /// MOF-5 preset (IRMOF-1): well-studied benchmark material.
    pub fn mof5() -> Self {
        // MOF-5: ~3800 m²/g, ~1.55 cm³/g, Q_st ~4–5 kJ/mol, ~7 wt% at 77 K / 40 bar
        Self::new(3800.0, 1.55, 4500.0, 7.0, 77.0, 0.05)
    }

    /// Langmuir constant at temperature `t_k` (K) via Clausius-Clapeyron.
    fn k_lang(&self, t_k: f64) -> f64 {
        self.k_lang_ref * (self.isosteric_heat / R_GAS * (1.0 / t_k - 1.0 / self.t_ref)).exp()
    }

    /// Hydrogen adsorption capacity (wt%) at temperature `t_k` (K) and
    /// pressure `p_bar` (bar) using a Langmuir isotherm.
    pub fn capacity_at_condition(&self, t_k: f64, p_bar: f64) -> f64 {
        let k = self.k_lang(t_k);
        self.saturation_capacity_wt * k * p_bar / (1.0 + k * p_bar)
    }
}

// ---------------------------------------------------------------------------
// Hydrogen Safety
// ---------------------------------------------------------------------------

/// Flammability and safety data for hydrogen.
#[derive(Clone, Debug)]
pub struct HydrogenSafety {
    /// Lower flammability limit in air (volume fraction, 0–1).
    pub lower_flammability_limit: f64,
    /// Upper flammability limit in air (volume fraction, 0–1).
    pub upper_flammability_limit: f64,
    /// Auto-ignition temperature (K).
    pub autoignition_temp: f64,
}

impl HydrogenSafety {
    /// Standard H₂ safety limits per IEC/ISO standards.
    ///
    /// LFL = 4 vol%, UFL = 75 vol%, auto-ignition ≈ 773 K (500 °C).
    pub fn standard() -> Self {
        Self {
            lower_flammability_limit: 0.04,
            upper_flammability_limit: 0.75,
            autoignition_temp: 773.0,
        }
    }

    /// Returns `true` when the volume fraction `concentration` (0–1) falls
    /// within the flammable range `[LFL, UFL]`.
    pub fn is_flammable_mixture(&self, concentration: f64) -> bool {
        concentration >= self.lower_flammability_limit
            && concentration <= self.upper_flammability_limit
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// H₂ density (kg / m³) at pressure `p_bar` (bar) and temperature `t_k` (K)
/// using the van der Waals equation of state.
///
/// The van der Waals equation `(P + a/V²)(V − b) = RT` is solved via
/// bisection on the molar volume V (m³/mol), then ρ = M/V.
pub fn h2_density_kg_m3(p_bar: f64, t_k: f64) -> f64 {
    let p_pa = p_bar * 1e5;
    // a and b for H₂ in SI units (Pa·m⁶/mol², m³/mol)
    // a = 0.02476 L²·atm/mol² = 0.02476e-3 m⁶ * 101325 Pa / mol²
    let a = 0.02476e-3 * 101_325.0; // Pa·m⁶/mol²
    let b = VDW_B_H2; // m³/mol

    // f(V) = (P + a/V²)(V − b) − R·T
    let f = |v: f64| -> f64 { (p_pa + a / (v * v)) * (v - b) - R_GAS * t_k };

    // Bracket: lower bound just above co-volume b, upper bound = 10× ideal-gas volume
    let v_lo = b * 1.001_f64;
    let v_hi = 10.0 * R_GAS * t_k / p_pa.max(1.0);

    // If the bracket does not straddle the root (can happen at very high P),
    // fall back to the ideal-gas result.
    if f(v_lo) * f(v_hi) > 0.0 {
        return M_H2 * p_pa / (R_GAS * t_k);
    }

    // Bisection — 80 iterations gives sub-ULP accuracy.
    let mut lo = v_lo;
    let mut hi = v_hi;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(lo) * f(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    M_H2 / (0.5 * (lo + hi))
}

/// Gravimetric energy density (MJ / kg) for a system storing hydrogen at the
/// given gravimetric capacity `wt_percent` (wt%).
///
/// The formula accounts only for the H₂ fraction of the system mass:
/// `E_grav = (wt% / 100) * LHV_H2`.
pub fn gravimetric_energy_density_mj_kg(wt_percent: f64) -> f64 {
    (wt_percent / 100.0) * LHV_H2_MJ_KG
}

/// Volumetric storage capacity (kg H₂ / m³) given the density of the storage
/// medium `density` (kg / m³) and the gravimetric capacity `wt_percent` (wt%).
pub fn volumetric_capacity(density: f64, wt_percent: f64) -> f64 {
    density * wt_percent / 100.0
}

/// Critical wrinkling wavelength for a thin sheet under biaxial tension —
/// re-exported here for completeness; used when sizing inflatable vessels.
///
/// `λ = 2π * sqrt(D / T)` where D is bending stiffness and T is tension.
pub fn vessel_wrinkling_wavelength(bending_stiffness: f64, tension: f64) -> f64 {
    2.0 * PI * (bending_stiffness / tension).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── van der Waals density ────────────────────────────────────────────────

    #[test]
    fn test_h2_density_ambient() {
        // At 1 bar / 298 K ideal-gas ≈ 0.0822 kg/m³; vdW should be very close.
        let rho = h2_density_kg_m3(1.0, 298.15);
        assert!(rho > 0.07 && rho < 0.10, "ρ(1 bar,298K)={rho}");
    }

    #[test]
    fn test_h2_density_700bar() {
        // At 700 bar / 298 K compressed H₂ (vdW) is in the range 30–100 kg/m³.
        let rho = h2_density_kg_m3(700.0, 298.15);
        assert!(rho > 30.0 && rho < 100.0, "ρ(700 bar,298K)={rho}");
    }

    #[test]
    fn test_h2_density_increases_with_pressure() {
        let rho_low = h2_density_kg_m3(50.0, 298.15);
        let rho_high = h2_density_kg_m3(200.0, 298.15);
        assert!(rho_high > rho_low, "Density must increase with pressure");
    }

    #[test]
    fn test_h2_density_decreases_with_temperature() {
        let rho_cold = h2_density_kg_m3(100.0, 200.0);
        let rho_hot = h2_density_kg_m3(100.0, 400.0);
        assert!(rho_cold > rho_hot, "Density must decrease with temperature");
    }

    #[test]
    fn test_h2_density_positive() {
        let rho = h2_density_kg_m3(350.0, 333.0);
        assert!(rho > 0.0, "Density must be positive");
    }

    // ── van 't Hoff equilibrium pressure ────────────────────────────────────

    #[test]
    fn test_lani5_equilibrium_pressure_room_temp() {
        let mh = MetalHydride::lani5();
        // LaNi₅ has plateau pressure ~2–3 bar at ~313 K.
        let p = mh.equilibrium_pressure(313.0);
        assert!(p > 0.5 && p < 20.0, "LaNi5 P_eq at 313 K = {p} bar");
    }

    #[test]
    fn test_mgh2_equilibrium_pressure_high_temp() {
        let mh = MetalHydride::mgh2();
        // MgH₂ needs ~573 K to reach ~1 bar desorption pressure.
        let p = mh.equilibrium_pressure(573.0);
        assert!(p > 0.1 && p < 50.0, "MgH2 P_eq at 573 K = {p} bar");
    }

    #[test]
    fn test_equilibrium_pressure_increases_with_temperature() {
        let mh = MetalHydride::lani5();
        let p_low = mh.equilibrium_pressure(300.0);
        let p_high = mh.equilibrium_pressure(350.0);
        assert!(
            p_high > p_low,
            "Higher T must give higher P_eq: {p_low} vs {p_high}"
        );
    }

    #[test]
    fn test_equilibrium_pressure_positive() {
        let mh = MetalHydride::lani5();
        assert!(mh.equilibrium_pressure(400.0) > 0.0);
    }

    #[test]
    fn test_metal_hydride_capacity_positive() {
        let mh = MetalHydride::mgh2();
        assert!(mh.max_capacity_wt_percent > 0.0);
    }

    #[test]
    fn test_metal_hydride_new_roundtrip() {
        let mh = MetalHydride::new("TestMH", 3.5, 400.0, -50_000.0, -120.0);
        assert_eq!(mh.metal_name, "TestMH");
        assert!((mh.max_capacity_wt_percent - 3.5).abs() < EPS);
    }

    // ── Pressure vessel ──────────────────────────────────────────────────────

    #[test]
    fn test_pressure_vessel_mass_positive() {
        let pv = PressureVesselH2::type_iv_700bar();
        let m = pv.mass_h2();
        assert!(m > 0.0, "Mass must be positive, got {m}");
    }

    #[test]
    fn test_pressure_vessel_mass_scales_with_volume() {
        let pv1 = PressureVesselH2::new(100.0, 700.0, 298.15);
        let pv2 = PressureVesselH2::new(200.0, 700.0, 298.15);
        let ratio = pv2.mass_h2() / pv1.mass_h2();
        assert!(
            (ratio - 2.0).abs() < 1e-9,
            "Mass must scale linearly with volume, ratio={ratio}"
        );
    }

    #[test]
    fn test_pressure_vessel_mass_scales_with_pressure() {
        let pv1 = PressureVesselH2::new(100.0, 350.0, 298.15);
        let pv2 = PressureVesselH2::new(100.0, 700.0, 298.15);
        let ratio = pv2.mass_h2() / pv1.mass_h2();
        assert!(
            (ratio - 2.0).abs() < 1e-9,
            "Mass must scale linearly with pressure, ratio={ratio}"
        );
    }

    #[test]
    fn test_pressure_vessel_energy_density_lhv() {
        let pv = PressureVesselH2::new(50.0, 350.0, 298.15);
        assert!((pv.energy_density() - LHV_H2_MJ_KG).abs() < EPS);
    }

    #[test]
    fn test_pressure_vessel_stored_energy() {
        let pv = PressureVesselH2::new(100.0, 350.0, 298.15);
        let expected = pv.mass_h2() * LHV_H2_MJ_KG;
        assert!((pv.stored_energy_mj() - expected).abs() < EPS);
    }

    #[test]
    fn test_pressure_vessel_type_iv_reasonable_mass() {
        // 700 bar, 120 L → should store ~4–8 kg H₂ (ideal gas approximation).
        let pv = PressureVesselH2::type_iv_700bar();
        let m = pv.mass_h2();
        assert!(
            m > 3.0 && m < 10.0,
            "Type-IV 700-bar vessel should store ~5 kg H₂, got {m}"
        );
    }

    // ── MOF adsorption ───────────────────────────────────────────────────────

    #[test]
    fn test_mof_capacity_positive() {
        let mof = MofAdsorption::mof5();
        let c = mof.capacity_at_condition(77.0, 40.0);
        assert!(c > 0.0, "MOF capacity must be positive, got {c}");
    }

    #[test]
    fn test_mof_capacity_below_saturation() {
        let mof = MofAdsorption::mof5();
        let c = mof.capacity_at_condition(77.0, 40.0);
        assert!(
            c <= mof.saturation_capacity_wt,
            "Capacity {c} must not exceed saturation {}",
            mof.saturation_capacity_wt
        );
    }

    #[test]
    fn test_mof_capacity_increases_with_pressure() {
        let mof = MofAdsorption::mof5();
        let c_low = mof.capacity_at_condition(77.0, 1.0);
        let c_high = mof.capacity_at_condition(77.0, 40.0);
        assert!(c_high > c_low, "Higher pressure must give higher capacity");
    }

    #[test]
    fn test_mof_capacity_decreases_with_temperature() {
        let mof = MofAdsorption::mof5();
        let c_cold = mof.capacity_at_condition(77.0, 10.0);
        let c_warm = mof.capacity_at_condition(298.0, 10.0);
        assert!(c_cold > c_warm, "Lower T must give higher physisorption");
    }

    #[test]
    fn test_mof_capacity_zero_pressure() {
        let mof = MofAdsorption::mof5();
        let c = mof.capacity_at_condition(77.0, 0.0);
        assert!(c.abs() < EPS, "Zero pressure must give zero capacity");
    }

    // ── Hydrogen safety ──────────────────────────────────────────────────────

    #[test]
    fn test_safety_lfl_not_flammable_below() {
        let safety = HydrogenSafety::standard();
        // 3% is below LFL (4%) → not flammable
        assert!(!safety.is_flammable_mixture(0.03));
    }

    #[test]
    fn test_safety_at_lfl_flammable() {
        let safety = HydrogenSafety::standard();
        assert!(safety.is_flammable_mixture(0.04));
    }

    #[test]
    fn test_safety_midrange_flammable() {
        let safety = HydrogenSafety::standard();
        assert!(safety.is_flammable_mixture(0.20));
    }

    #[test]
    fn test_safety_at_ufl_flammable() {
        let safety = HydrogenSafety::standard();
        assert!(safety.is_flammable_mixture(0.75));
    }

    #[test]
    fn test_safety_above_ufl_not_flammable() {
        let safety = HydrogenSafety::standard();
        // 80% is above UFL (75%) → rich mixture, not flammable
        assert!(!safety.is_flammable_mixture(0.80));
    }

    #[test]
    fn test_safety_zero_concentration_not_flammable() {
        let safety = HydrogenSafety::standard();
        assert!(!safety.is_flammable_mixture(0.0));
    }

    #[test]
    fn test_safety_autoignition_temp() {
        let safety = HydrogenSafety::standard();
        assert!((safety.autoignition_temp - 773.0).abs() < EPS);
    }

    #[test]
    fn test_safety_lfl_4_percent() {
        let safety = HydrogenSafety::standard();
        assert!((safety.lower_flammability_limit - 0.04).abs() < EPS);
    }

    #[test]
    fn test_safety_ufl_75_percent() {
        let safety = HydrogenSafety::standard();
        assert!((safety.upper_flammability_limit - 0.75).abs() < EPS);
    }

    // ── Free functions ───────────────────────────────────────────────────────

    #[test]
    fn test_gravimetric_energy_density_7wt() {
        // 7 wt% → (7/100) * 119.96 ≈ 8.397 MJ/kg
        let e = gravimetric_energy_density_mj_kg(7.0);
        assert!((e - 0.07 * LHV_H2_MJ_KG).abs() < 1e-9, "Got {e}");
    }

    #[test]
    fn test_gravimetric_energy_density_scales_linearly() {
        let e1 = gravimetric_energy_density_mj_kg(2.0);
        let e2 = gravimetric_energy_density_mj_kg(4.0);
        assert!((e2 - 2.0 * e1).abs() < EPS, "Must scale linearly");
    }

    #[test]
    fn test_volumetric_capacity() {
        // density=1000 kg/m³, 5 wt% → 50 kg H₂/m³
        let vc = volumetric_capacity(1000.0, 5.0);
        assert!((vc - 50.0).abs() < EPS, "Got {vc}");
    }

    #[test]
    fn test_volumetric_capacity_zero_wt() {
        let vc = volumetric_capacity(800.0, 0.0);
        assert!(vc.abs() < EPS);
    }

    #[test]
    fn test_vessel_wrinkling_wavelength() {
        let lam = vessel_wrinkling_wavelength(1.0, 1.0);
        let expected = 2.0 * PI;
        assert!((lam - expected).abs() < 1e-10, "Got {lam}");
    }

    // ── HydrogenStorageType enum ─────────────────────────────────────────────

    #[test]
    fn test_storage_type_eq() {
        assert_eq!(
            HydrogenStorageType::CompressedGas,
            HydrogenStorageType::CompressedGas
        );
        assert_ne!(HydrogenStorageType::LiquidH2, HydrogenStorageType::MOF);
    }

    #[test]
    fn test_storage_type_clone() {
        let t = HydrogenStorageType::MetalHydride;
        let t2 = t.clone();
        assert_eq!(t, t2);
    }

    // ── Edge / regression cases ──────────────────────────────────────────────

    #[test]
    fn test_mof_new_fields_accessible() {
        let mof = MofAdsorption::new(2000.0, 1.0, 5000.0, 5.0, 77.0, 0.03);
        assert!((mof.surface_area_m2g - 2000.0).abs() < EPS);
        assert!((mof.pore_volume - 1.0).abs() < EPS);
    }

    #[test]
    fn test_h2_density_linearity_low_pressure() {
        // At very low pressures vdW ≈ ideal gas: ρ ∝ P (within ~5%).
        let rho_1 = h2_density_kg_m3(1.0, 300.0);
        let rho_2 = h2_density_kg_m3(2.0, 300.0);
        let ratio = rho_2 / rho_1;
        assert!(
            (ratio - 2.0).abs() < 0.15,
            "Low-P density should be ~linear in P: ratio={ratio}"
        );
    }
}

// ---------------------------------------------------------------------------
// Extended hydrogen storage functions (van der Waals, Sievert, Clausius-Clapeyron,
// Langmuir, BEP, HydrideStorage, energy density, storage efficiency, safe pressure)
// ---------------------------------------------------------------------------

/// Molar volume of H₂ \[m³/mol\] via the van der Waals equation of state.
///
/// Solves `(P + a/V²)(V − b) = RT` by bisection.
///
/// * `p`   – pressure \[Pa\]
/// * `t`   – temperature \[K\]
/// * `a`   – van der Waals constant a \[Pa·m⁶/mol²\]
/// * `b`   – van der Waals constant b \[m³/mol\]
pub fn van_der_waals_h2(p: f64, t: f64, a: f64, b: f64) -> f64 {
    // f(V) = (P + a/V²)(V − b) − RT
    let f = |v: f64| -> f64 { (p + a / (v * v)) * (v - b) - R_GAS * t };

    let v_lo = b * 1.001_f64;
    let v_hi = 10.0 * R_GAS * t / p.max(1.0);

    // Fall back to ideal-gas volume if bracket has no sign change.
    if f(v_lo) * f(v_hi) > 0.0 {
        return R_GAS * t / p;
    }

    let mut lo = v_lo;
    let mut hi = v_hi;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(lo) * f(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Sievert's law: hydrogen concentration in a metal \[mol/m³ or wt-frac equivalent\].
///
/// c = K_S · √P_H₂
///
/// * `k_s`   – Sievert constant
/// * `p_h2`  – hydrogen partial pressure \[Pa\]
pub fn sievert_law(k_s: f64, p_h2: f64) -> f64 {
    k_s * p_h2.sqrt()
}

/// Clausius-Clapeyron pressure ratio for hydrogen desorption.
///
/// P₂/P₁ = exp(−ΔH/R · (1/T₂ − 1/T₁))
///
/// * `delta_h` – enthalpy of desorption \[J/mol\] (positive for endothermic)
/// * `t1`      – reference temperature \[K\]
/// * `t2`      – target temperature \[K\]
///
/// Returns P₂/P₁.
pub fn clausius_clapeyron_h2(delta_h: f64, t1: f64, t2: f64) -> f64 {
    (-delta_h / R_GAS * (1.0 / t2 - 1.0 / t1)).exp()
}

/// Langmuir adsorption isotherm: fractional surface coverage.
///
/// θ = k · P / (1 + k · P)
///
/// * `k`        – Langmuir equilibrium constant \[1/Pa\]
/// * `pressure` – pressure \[Pa\]
///
/// Returns θ ∈ \[0, 1).
pub fn langmuir_hydrogen(k: f64, pressure: f64) -> f64 {
    k * pressure / (1.0 + k * pressure)
}

/// Brønsted-Evans-Polanyi (BEP) activation energy \[J/mol\].
///
/// E_a = E_a0 + α · ΔE
///
/// * `barrier`  – intrinsic barrier E_a0 \[J/mol\]
/// * `delta_e`  – reaction energy ΔE \[J/mol\]
/// * `alpha`    – BEP slope (0 < α < 1, typically 0.5)
pub fn bep_relationship(barrier: f64, delta_e: f64, alpha: f64) -> f64 {
    barrier + alpha * delta_e
}

/// A simple metal-hydride hydrogen storage tank model.
///
/// Encapsulates capacity, operating temperature and pressure, and an
/// Arrhenius activation energy for desorption kinetics.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HydrideStorage {
    /// Maximum gravimetric capacity \[wt%\].
    pub capacity_wt_pct: f64,
    /// Operating temperature \[K\].
    pub temperature: f64,
    /// Operating pressure \[Pa\].
    pub pressure: f64,
    /// Activation energy for desorption \[J/mol\].
    pub activation_energy: f64,
}

impl HydrideStorage {
    /// Creates a new `HydrideStorage` tank.
    ///
    /// `activation_energy` defaults to 50 kJ/mol if not specified; use
    /// [`HydrideStorage::with_activation_energy`] for full control.
    pub fn new(cap: f64, t: f64, p: f64) -> Self {
        Self {
            capacity_wt_pct: cap,
            temperature: t,
            pressure: p,
            activation_energy: 50_000.0,
        }
    }

    /// Creates a `HydrideStorage` with an explicit activation energy.
    pub fn with_activation_energy(cap: f64, t: f64, p: f64, ea: f64) -> Self {
        Self {
            capacity_wt_pct: cap,
            temperature: t,
            pressure: p,
            activation_energy: ea,
        }
    }

    /// Desorption rate \[wt%/s\] via an Arrhenius expression.
    ///
    /// rate = cap · A · exp(−E_a / (R · T))
    ///
    /// where A = 1e6 s⁻¹ is a pre-exponential factor.
    pub fn desorption_rate(&self) -> f64 {
        let a_pre = 1.0e6_f64;
        self.capacity_wt_pct * a_pre * (-self.activation_energy / (R_GAS * self.temperature)).exp()
    }

    /// Gravimetric density \[kWh_H₂/kg_H₂\] — numerically equal to LHV.
    ///
    /// Returns wt% fraction times the lower heating value of H₂.
    pub fn gravimetric_density(&self) -> f64 {
        (self.capacity_wt_pct / 100.0) * LHV_H2_MJ_KG / 3.6 // MJ/kg → kWh/kg
    }

    /// Time to charge from 0 to `target_wt_pct` at the current desorption rate \[s\].
    ///
    /// Approximates charging as the inverse of the desorption rate scaled to
    /// target capacity.
    pub fn charging_time(&self, target_wt_pct: f64) -> f64 {
        let rate = self.desorption_rate();
        if rate <= 0.0 {
            return f64::INFINITY;
        }
        target_wt_pct / rate
    }
}

/// Volumetric energy density of compressed H₂ \[kWh/L\].
///
/// * `pressure_mpa`  – storage pressure \[MPa\]
/// * `temperature_k` – gas temperature \[K\]
pub fn hydrogen_energy_density_volumetric(pressure_mpa: f64, temperature_k: f64) -> f64 {
    // ideal-gas density: ρ = P*M/(R*T) [kg/m³]
    let p_pa = pressure_mpa * 1.0e6;
    let rho = p_pa * M_H2 / (R_GAS * temperature_k);
    // energy: rho [kg/m³] * LHV [MJ/kg] / 3.6 → kWh/m³, then /1000 → kWh/L
    rho * LHV_H2_MJ_KG / 3.6 / 1000.0
}

/// Specific energy of a hydrogen storage system \[kWh/kg_system\].
///
/// * `useful_energy`      – usable H₂ energy \[kWh\]
/// * `total_system_mass`  – total mass of system including vessel \[kg\]
pub fn storage_efficiency(useful_energy: f64, total_system_mass: f64) -> f64 {
    useful_energy / total_system_mass
}

/// Safe maximum pressure for a cylindrical pressure vessel \[Pa\] from hoop stress.
///
/// P = 2 · t · UTS / D
///
/// * `vessel_diameter` – inner diameter \[m\]
/// * `wall_thickness`  – wall thickness \[m\]
/// * `uts`             – ultimate tensile strength of wall material \[Pa\]
pub fn safe_pressure_limit(vessel_diameter: f64, wall_thickness: f64, uts: f64) -> f64 {
    2.0 * wall_thickness * uts / vessel_diameter
}

// ---------------------------------------------------------------------------
// Extended tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod extended_tests {

    use crate::additive_manufacturing_materials::R_GAS;
    use crate::hydrogen_storage::HydrideStorage;
    use crate::hydrogen_storage::bep_relationship;
    use crate::hydrogen_storage::clausius_clapeyron_h2;
    use crate::hydrogen_storage::hydrogen_energy_density_volumetric;
    use crate::hydrogen_storage::langmuir_hydrogen;
    use crate::hydrogen_storage::safe_pressure_limit;
    use crate::hydrogen_storage::sievert_law;
    use crate::hydrogen_storage::storage_efficiency;
    use crate::hydrogen_storage::van_der_waals_h2;

    const EPS: f64 = 1e-10;
    /// van der Waals constants for H₂ (SI units).
    const A_VDW: f64 = 0.02476e-3 * 101_325.0;
    const B_VDW: f64 = 26.61e-6;

    // ── van_der_waals_h2 ──────────────────────────────────────────────────

    #[test]
    fn test_vdw_positive_volume() {
        let v = van_der_waals_h2(1e5, 300.0, A_VDW, B_VDW);
        assert!(v > 0.0, "Molar volume must be positive, got {v}");
    }

    #[test]
    fn test_vdw_volume_decreases_with_pressure() {
        let v_low = van_der_waals_h2(1e5, 300.0, A_VDW, B_VDW);
        let v_high = van_der_waals_h2(1e6, 300.0, A_VDW, B_VDW);
        assert!(v_low > v_high, "Higher pressure → smaller molar volume");
    }

    #[test]
    fn test_vdw_approaches_ideal_at_low_p() {
        // At 1 bar / 300 K: V_ideal = RT/P ≈ 0.02494 m³/mol.
        // vdW volume differs from ideal by a few percent due to a and b corrections.
        let v = van_der_waals_h2(1e5, 300.0, A_VDW, B_VDW);
        let v_ideal = R_GAS * 300.0 / 1e5;
        // Allow up to 5% deviation at 1 bar.
        assert!(
            (v - v_ideal).abs() / v_ideal < 0.05,
            "vdW should ≈ ideal at low P: vdW={v}, ideal={v_ideal}"
        );
    }

    #[test]
    fn test_vdw_finite() {
        let v = van_der_waals_h2(7e7, 300.0, A_VDW, B_VDW);
        assert!(v.is_finite() && v > 0.0);
    }

    // ── sievert_law ───────────────────────────────────────────────────────

    #[test]
    fn test_sievert_positive() {
        let c = sievert_law(1e-3, 1e5);
        assert!(c > 0.0, "Sievert concentration must be positive");
    }

    #[test]
    fn test_sievert_sqrt_dependence() {
        // c ∝ √P: doubling P should give √2 × c.
        let c1 = sievert_law(1e-3, 1e5);
        let c2 = sievert_law(1e-3, 4e5);
        let ratio = c2 / c1;
        assert!(
            (ratio - 2.0).abs() < 1e-9,
            "Sievert: c should scale as √P, ratio={ratio}"
        );
    }

    #[test]
    fn test_sievert_zero_pressure() {
        let c = sievert_law(1e-3, 0.0);
        assert!(c.abs() < EPS, "Zero pressure → zero concentration");
    }

    #[test]
    fn test_sievert_increases_with_pressure() {
        let c1 = sievert_law(1e-3, 1e5);
        let c2 = sievert_law(1e-3, 9e5);
        assert!(c2 > c1);
    }

    // ── clausius_clapeyron_h2 ─────────────────────────────────────────────

    #[test]
    fn test_clausius_clapeyron_positive_ratio() {
        // Positive ΔH (endothermic desorption), T2 > T1 → P2/P1 > 1.
        let ratio = clausius_clapeyron_h2(30_000.0, 300.0, 350.0);
        assert!(
            ratio > 1.0,
            "Higher T → higher equilibrium P for endothermic: ratio={ratio}"
        );
    }

    #[test]
    fn test_clausius_clapeyron_t1_eq_t2() {
        let ratio = clausius_clapeyron_h2(30_000.0, 300.0, 300.0);
        assert!((ratio - 1.0).abs() < 1e-9, "Equal temps → ratio=1");
    }

    #[test]
    fn test_clausius_clapeyron_finite() {
        let ratio = clausius_clapeyron_h2(74_500.0, 573.0, 700.0);
        assert!(ratio.is_finite() && ratio > 0.0);
    }

    // ── langmuir_hydrogen ─────────────────────────────────────────────────

    #[test]
    fn test_langmuir_in_zero_one() {
        let theta = langmuir_hydrogen(1e-5, 1e5);
        assert!(
            (0.0..1.0).contains(&theta),
            "Langmuir theta must be in [0,1): {theta}"
        );
    }

    #[test]
    fn test_langmuir_increases_with_pressure() {
        let t1 = langmuir_hydrogen(1e-5, 1e5);
        let t2 = langmuir_hydrogen(1e-5, 1e6);
        assert!(t2 > t1, "Higher pressure → higher Langmuir coverage");
    }

    #[test]
    fn test_langmuir_zero_pressure() {
        let theta = langmuir_hydrogen(1e-5, 0.0);
        assert!(theta.abs() < EPS, "Zero pressure → zero Langmuir coverage");
    }

    #[test]
    fn test_langmuir_approaches_one_high_pressure() {
        let theta = langmuir_hydrogen(1.0, 1e9);
        assert!(theta > 0.999, "Very high P → θ → 1, got {theta}");
    }

    // ── bep_relationship ──────────────────────────────────────────────────

    #[test]
    fn test_bep_positive_barrier() {
        let ea = bep_relationship(50_000.0, 10_000.0, 0.5);
        assert!(ea > 0.0);
    }

    #[test]
    fn test_bep_formula() {
        let ea = bep_relationship(50_000.0, 20_000.0, 0.3);
        let expected = 50_000.0 + 0.3 * 20_000.0;
        assert!((ea - expected).abs() < EPS);
    }

    // ── HydrideStorage ────────────────────────────────────────────────────

    #[test]
    fn test_hydride_storage_desorption_rate_positive() {
        let hs = HydrideStorage::new(6.0, 573.0, 1e5);
        assert!(
            hs.desorption_rate() > 0.0,
            "Desorption rate must be positive"
        );
    }

    #[test]
    fn test_hydride_storage_desorption_rate_increases_with_temperature() {
        let hs_cold = HydrideStorage::new(6.0, 400.0, 1e5);
        let hs_hot = HydrideStorage::new(6.0, 600.0, 1e5);
        assert!(
            hs_hot.desorption_rate() > hs_cold.desorption_rate(),
            "Higher T → faster desorption"
        );
    }

    #[test]
    fn test_hydride_storage_gravimetric_density_positive() {
        let hs = HydrideStorage::new(6.0, 573.0, 1e5);
        assert!(hs.gravimetric_density() > 0.0);
    }

    #[test]
    fn test_hydride_storage_charging_time_positive() {
        let hs = HydrideStorage::new(6.0, 573.0, 1e5);
        let t = hs.charging_time(5.0);
        assert!(
            t > 0.0 && t.is_finite(),
            "Charging time must be positive and finite"
        );
    }

    #[test]
    fn test_hydride_storage_new_fields() {
        let hs = HydrideStorage::new(7.6, 573.0, 1e5);
        assert!((hs.capacity_wt_pct - 7.6).abs() < EPS);
        assert!((hs.temperature - 573.0).abs() < EPS);
    }

    // ── hydrogen_energy_density_volumetric ────────────────────────────────

    #[test]
    fn test_energy_density_positive() {
        let e = hydrogen_energy_density_volumetric(70.0, 298.15);
        assert!(e > 0.0, "Energy density must be positive");
    }

    #[test]
    fn test_energy_density_increases_with_pressure() {
        let e1 = hydrogen_energy_density_volumetric(35.0, 298.15);
        let e2 = hydrogen_energy_density_volumetric(70.0, 298.15);
        assert!(
            e2 > e1,
            "Higher pressure → higher volumetric energy density"
        );
    }

    // ── storage_efficiency ────────────────────────────────────────────────

    #[test]
    fn test_storage_efficiency_positive() {
        let eff = storage_efficiency(10.0, 100.0);
        assert!(eff > 0.0);
    }

    #[test]
    fn test_storage_efficiency_formula() {
        let eff = storage_efficiency(15.0, 75.0);
        assert!((eff - 0.2).abs() < 1e-12);
    }

    // ── safe_pressure_limit ───────────────────────────────────────────────

    #[test]
    fn test_safe_pressure_positive() {
        let p = safe_pressure_limit(0.3, 0.01, 500e6);
        assert!(p > 0.0, "Safe pressure must be positive");
    }

    #[test]
    fn test_safe_pressure_increases_with_wall_thickness() {
        let p_thin = safe_pressure_limit(0.3, 0.005, 500e6);
        let p_thick = safe_pressure_limit(0.3, 0.015, 500e6);
        assert!(p_thick > p_thin, "Thicker wall → higher safe pressure");
    }

    #[test]
    fn test_safe_pressure_formula() {
        let d = 0.4_f64;
        let t = 0.01_f64;
        let uts = 800e6_f64;
        let expected = 2.0 * t * uts / d;
        let got = safe_pressure_limit(d, t, uts);
        assert!((got - expected).abs() < 1e-6);
    }

    #[test]
    fn test_safe_pressure_increases_with_uts() {
        let p1 = safe_pressure_limit(0.3, 0.01, 400e6);
        let p2 = safe_pressure_limit(0.3, 0.01, 800e6);
        assert!(p2 > p1, "Higher UTS → higher safe pressure");
    }
}
