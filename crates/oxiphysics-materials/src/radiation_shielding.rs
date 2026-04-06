// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Radiation shielding and nuclear material models.
//!
//! Implements attenuation calculations, half-value layer, buildup factors,
//! dose rate from kerma, and radioactivation analysis for gamma and neutron
//! radiation in shielding materials.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::LN_2;

// ─────────────────────────────────────────────────────────────────────────────
// Core struct
// ─────────────────────────────────────────────────────────────────────────────

/// A material used for radiation shielding.
#[derive(Debug, Clone)]
pub struct RadiationMaterial {
    /// Name of the material (e.g. "Lead", "Concrete").
    pub name: String,
    /// Atomic number Z of the primary constituent element.
    pub atomic_number: u32,
    /// Mass density \[kg/m³\].
    pub density: f64,
    /// Mass attenuation coefficient μ/ρ \[m²/kg\] at the relevant photon energy.
    pub mass_attenuation: f64,
    /// Thermal neutron absorption cross-section \[barn = 1e-28 m²\].
    pub neutron_cross_section: f64,
}

impl RadiationMaterial {
    /// Create a new radiation shielding material.
    ///
    /// # Arguments
    /// * `name` – Human-readable label.
    /// * `atomic_number` – Atomic number Z.
    /// * `density` – Mass density \[kg/m³\].
    /// * `mass_attenuation` – μ/ρ at the relevant energy \[m²/kg\].
    /// * `neutron_cross_section` – σ_a in barn for thermal neutrons.
    pub fn new(
        name: impl Into<String>,
        atomic_number: u32,
        density: f64,
        mass_attenuation: f64,
        neutron_cross_section: f64,
    ) -> Self {
        Self {
            name: name.into(),
            atomic_number,
            density,
            mass_attenuation,
            neutron_cross_section,
        }
    }

    /// Preset: lead (Pb).  μ/ρ at 1 MeV ≈ 7.1e-3 m²/kg.
    pub fn lead() -> Self {
        Self::new("Lead", 82, 11_340.0, 7.1e-3, 0.171)
    }

    /// Preset: ordinary concrete.  μ/ρ at 1 MeV ≈ 6.5e-3 m²/kg.
    pub fn concrete() -> Self {
        Self::new("Concrete", 14, 2_300.0, 6.5e-3, 0.17)
    }

    /// Preset: water.  μ/ρ at 1 MeV ≈ 6.7e-3 m²/kg.
    pub fn water() -> Self {
        Self::new("Water", 1, 1_000.0, 6.7e-3, 0.333)
    }

    /// Preset: polyethylene (high-density neutron moderator).
    pub fn polyethylene() -> Self {
        Self::new("Polyethylene", 6, 950.0, 2.0e-3, 0.333)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Attenuation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the linear attenuation coefficient μ = ρ · (μ/ρ) \[1/m\].
///
/// # Arguments
/// * `density` – Material density \[kg/m³\].
/// * `mass_attenuation` – Mass attenuation coefficient μ/ρ \[m²/kg\].
pub fn linear_attenuation(density: f64, mass_attenuation: f64) -> f64 {
    density * mass_attenuation
}

/// Compute the transmitted intensity fraction after thickness `x` \[m\]:
/// I/I₀ = exp(-μ · x), neglecting buildup.
///
/// # Arguments
/// * `mu` – Linear attenuation coefficient \[1/m\].
/// * `thickness` – Shield thickness \[m\].
pub fn narrow_beam_transmission(mu: f64, thickness: f64) -> f64 {
    (-mu * thickness).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Half-Value Layer
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the half-value layer HVL = ln(2) / μ \[m\].
///
/// The HVL is the thickness of material that attenuates a narrow beam of
/// radiation to half its original intensity.
///
/// # Arguments
/// * `mu` – Linear attenuation coefficient \[1/m\].
pub fn half_value_layer(mu: f64) -> f64 {
    if mu.abs() < 1e-300 {
        f64::INFINITY
    } else {
        LN_2 / mu
    }
}

/// Compute the tenth-value layer TVL = ln(10) / μ \[m\].
///
/// # Arguments
/// * `mu` – Linear attenuation coefficient \[1/m\].
pub fn tenth_value_layer(mu: f64) -> f64 {
    if mu.abs() < 1e-300 {
        f64::INFINITY
    } else {
        10.0_f64.ln() / mu
    }
}

/// Determine the number of HVLs needed to achieve a given transmission fraction.
///
/// # Arguments
/// * `transmission` – Desired I/I₀ (0 < transmission ≤ 1).
pub fn hvls_for_transmission(transmission: f64) -> f64 {
    assert!(transmission > 0.0 && transmission <= 1.0);
    -transmission.log2()
}

// ─────────────────────────────────────────────────────────────────────────────
// Buildup Factor (Taylor two-term approximation)
// ─────────────────────────────────────────────────────────────────────────────

/// Taylor two-term buildup factor B(μx) = A·exp(α₁·μx) + (1-A)·exp(α₂·μx).
///
/// Accounts for scattered photons that still deposit dose inside the shield.
///
/// # Arguments
/// * `mu_x` – Dimensionless penetration depth (μ · x).
/// * `cap_a` – Taylor coefficient A (material and energy dependent).
/// * `alpha1` – Taylor exponent α₁.
/// * `alpha2` – Taylor exponent α₂.
pub fn buildup_factor(mu_x: f64, cap_a: f64, alpha1: f64, alpha2: f64) -> f64 {
    cap_a * (alpha1 * mu_x).exp() + (1.0 - cap_a) * (alpha2 * mu_x).exp()
}

/// Berger buildup factor B(μx) = 1 + C·(μx)·exp(D·μx).
///
/// Alternative simple form used for low-Z materials.
///
/// # Arguments
/// * `mu_x` – Dimensionless penetration depth (μ · x).
/// * `c_coeff` – Berger coefficient C.
/// * `d_coeff` – Berger coefficient D.
pub fn berger_buildup_factor(mu_x: f64, c_coeff: f64, d_coeff: f64) -> f64 {
    1.0 + c_coeff * mu_x * (d_coeff * mu_x).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Dose Rate
// ─────────────────────────────────────────────────────────────────────────────

/// Compute air-kerma dose rate \[Gy/s\] from photon fluence rate and energy.
///
/// K̇ = Φ̇ · E · (μ_en/ρ)_air
///
/// # Arguments
/// * `fluence_rate` – Photon fluence rate \[photons/m²/s\].
/// * `energy_j` – Photon energy \[J\].
/// * `mu_en_over_rho_air` – Mass energy-absorption coefficient of air \[m²/kg\]
///   (≈ 2.7e-3 m²/kg at 1 MeV).
pub fn dose_rate(fluence_rate: f64, energy_j: f64, mu_en_over_rho_air: f64) -> f64 {
    fluence_rate * energy_j * mu_en_over_rho_air
}

/// Compute dose rate behind a slab shield \[Gy/s\], including buildup.
///
/// D̊ = Φ̇₀ · E · (μ_en/ρ)_air · B(μx) · exp(-μx)
///
/// # Arguments
/// * `fluence_rate_0` – Source-side photon fluence rate \[photons/m²/s\].
/// * `energy_j` – Photon energy \[J\].
/// * `mu_en_over_rho_air` – (μ_en/ρ) for air \[m²/kg\].
/// * `mu` – Linear attenuation coefficient of shielding material \[1/m\].
/// * `thickness` – Shield thickness \[m\].
/// * `buildup` – Pre-computed buildup factor B(μx).
pub fn shielded_dose_rate(
    fluence_rate_0: f64,
    energy_j: f64,
    mu_en_over_rho_air: f64,
    mu: f64,
    thickness: f64,
    buildup: f64,
) -> f64 {
    let unshielded = dose_rate(fluence_rate_0, energy_j, mu_en_over_rho_air);
    unshielded * buildup * (-mu * thickness).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Radioactivation
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a radioactivation analysis.
#[derive(Debug, Clone, Copy)]
pub struct ActivationResult {
    /// Saturation activity \[Bq\].
    pub saturation_activity: f64,
    /// Activity at end of irradiation \[Bq\].
    pub activity_at_eoi: f64,
    /// Activity after cooling time t_cool \[Bq\].
    pub activity_after_cooling: f64,
    /// Decay constant λ \[1/s\].
    pub decay_constant: f64,
}

/// Radioactivation under constant neutron flux (single-nuclide model).
///
/// Solves A(t) = N·σ·Φ·(1 − exp(−λ·t_irr)) · exp(−λ·t_cool)
/// where N = (m/M)·N_A is the number of target atoms.
///
/// # Arguments
/// * `mass_kg` – Mass of target nuclide \[kg\].
/// * `molar_mass` – Molar mass \[kg/mol\].
/// * `cross_section_m2` – Microscopic activation cross-section \[m²\].
/// * `neutron_flux` – Thermal neutron flux \[n/m²/s\].
/// * `half_life_s` – Half-life of the product nuclide \[s\].
/// * `irradiation_time_s` – Duration of irradiation \[s\].
/// * `cooling_time_s` – Decay time after end of irradiation \[s\].
pub fn activation_analysis(
    mass_kg: f64,
    molar_mass: f64,
    cross_section_m2: f64,
    neutron_flux: f64,
    half_life_s: f64,
    irradiation_time_s: f64,
    cooling_time_s: f64,
) -> ActivationResult {
    const AVOGADRO: f64 = 6.022_140_76e23;
    let n_atoms = (mass_kg / molar_mass) * AVOGADRO;
    let lambda = LN_2 / half_life_s;
    let saturation_activity = n_atoms * cross_section_m2 * neutron_flux * lambda;
    // Activity at end of irradiation:
    let activity_at_eoi = saturation_activity * (1.0 - (-lambda * irradiation_time_s).exp());
    // After cooling:
    let activity_after_cooling = activity_at_eoi * (-lambda * cooling_time_s).exp();
    ActivationResult {
        saturation_activity,
        activity_at_eoi,
        activity_after_cooling,
        decay_constant: lambda,
    }
}

/// Compute the effective dose rate \[Sv/s\] from activity using a dose
/// conversion factor (DCF) \[Sv/Bq\].
///
/// # Arguments
/// * `activity_bq` – Activity of the source \[Bq\].
/// * `dcf` – Dose conversion factor \[Sv/Bq\].
pub fn effective_dose_rate(activity_bq: f64, dcf: f64) -> f64 {
    activity_bq * dcf
}

// ─────────────────────────────────────────────────────────────────────────────
// Neutron shielding helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute neutron removal cross-section contribution Σ_r = N · σ_r
/// where N = (ρ/M)·N_A is the atomic number density.
///
/// # Arguments
/// * `density` – Material density \[kg/m³\].
/// * `molar_mass` – Molar mass \[kg/mol\].
/// * `removal_cross_section_m2` – Fast-neutron removal cross-section \[m²\].
pub fn macroscopic_removal_cross_section(
    density: f64,
    molar_mass: f64,
    removal_cross_section_m2: f64,
) -> f64 {
    const AVOGADRO: f64 = 6.022_140_76e23;
    let number_density = (density / molar_mass) * AVOGADRO;
    number_density * removal_cross_section_m2
}

/// Compute the mean free path for neutrons: λ = 1 / Σ_t \[m\].
///
/// # Arguments
/// * `macroscopic_cross_section` – Total macroscopic cross-section Σ_t \[1/m\].
pub fn mean_free_path(macroscopic_cross_section: f64) -> f64 {
    1.0 / macroscopic_cross_section
}

// ─────────────────────────────────────────────────────────────────────────────
// Shield design helper
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the required shield thickness to achieve a target dose-rate reduction
/// factor R (= D̊_0 / D̊_target), neglecting buildup.
///
/// x = ln(R) / μ
///
/// # Arguments
/// * `reduction_factor` – D̊_0 / D̊_target (must be > 1).
/// * `mu` – Linear attenuation coefficient \[1/m\].
pub fn required_thickness(reduction_factor: f64, mu: f64) -> f64 {
    assert!(reduction_factor > 1.0);
    reduction_factor.ln() / mu
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Linear attenuation: μ = ρ * (μ/ρ)
    #[test]
    fn test_linear_attenuation_lead() {
        let mat = RadiationMaterial::lead();
        let mu = linear_attenuation(mat.density, mat.mass_attenuation);
        // Lead: 11340 * 7.1e-3 ≈ 80.5 m⁻¹
        assert!(
            (mu - 80.514).abs() < 0.5,
            "Lead μ should be ~80.5, got {mu}"
        );
    }

    // 2. μ = 0 → infinite HVL
    #[test]
    fn test_hvl_zero_mu() {
        assert_eq!(half_value_layer(0.0), f64::INFINITY);
    }

    // 3. HVL × μ = ln(2)
    #[test]
    fn test_hvl_identity() {
        let mu = 50.0;
        let hvl = half_value_layer(mu);
        assert!((hvl * mu - LN_2).abs() < EPS, "HVL * μ must equal ln2");
    }

    // 4. Narrow-beam transmission at x = HVL should be 0.5
    #[test]
    fn test_narrow_beam_half_at_hvl() {
        let mu = 30.0;
        let hvl = half_value_layer(mu);
        let t = narrow_beam_transmission(mu, hvl);
        assert!(
            (t - 0.5).abs() < EPS,
            "Transmission at HVL must be 0.5, got {t}"
        );
    }

    // 5. Narrow-beam transmission at x = 0 must be 1.0
    #[test]
    fn test_narrow_beam_zero_thickness() {
        let t = narrow_beam_transmission(100.0, 0.0);
        assert!((t - 1.0).abs() < EPS);
    }

    // 6. Transmission decreases monotonically with thickness.
    #[test]
    fn test_transmission_monotone() {
        let mu = 20.0;
        let t1 = narrow_beam_transmission(mu, 0.01);
        let t2 = narrow_beam_transmission(mu, 0.05);
        let t3 = narrow_beam_transmission(mu, 0.1);
        assert!(
            t1 > t2 && t2 > t3,
            "Transmission must decrease with thickness"
        );
    }

    // 7. Tenth-value layer: TVL * μ = ln(10)
    #[test]
    fn test_tvl_identity() {
        let mu = 80.0;
        let tvl = tenth_value_layer(mu);
        assert!(
            (tvl * mu - 10.0_f64.ln()).abs() < EPS,
            "TVL * μ must equal ln(10)"
        );
    }

    // 8. TVL = HVL * log2(10)
    #[test]
    fn test_tvl_vs_hvl_ratio() {
        let mu = 45.0;
        let hvl = half_value_layer(mu);
        let tvl = tenth_value_layer(mu);
        let ratio = tvl / hvl;
        let expected = 10.0_f64.log2(); // ≈ 3.3219
        assert!(
            (ratio - expected).abs() < EPS,
            "TVL/HVL must equal log2(10), got {ratio}"
        );
    }

    // 9. Taylor buildup at μx = 0 equals 1.0 (no penetration → no scatter).
    #[test]
    fn test_taylor_buildup_at_zero() {
        let b = buildup_factor(0.0, 0.5, -0.1, 0.05);
        assert!((b - 1.0).abs() < EPS, "Buildup at μx=0 must be 1, got {b}");
    }

    // 10. Berger buildup at μx = 0 equals 1.0.
    #[test]
    fn test_berger_buildup_at_zero() {
        let b = berger_buildup_factor(0.0, 2.0, 0.3);
        assert!(
            (b - 1.0).abs() < EPS,
            "Berger buildup at μx=0 must be 1, got {b}"
        );
    }

    // 11. Berger buildup is always ≥ 1 for positive μx.
    #[test]
    fn test_berger_buildup_gte_one() {
        for &mu_x in &[0.5, 1.0, 2.0, 5.0, 10.0] {
            let b = berger_buildup_factor(mu_x, 1.5, 0.2);
            assert!(b >= 1.0, "Berger buildup must be ≥ 1, got {b} at μx={mu_x}");
        }
    }

    // 12. Dose rate scales linearly with fluence rate.
    #[test]
    fn test_dose_rate_linear_fluence() {
        let e_j = 1.6e-13; // 1 MeV in J
        let mu_en = 2.7e-3;
        let d1 = dose_rate(1e6, e_j, mu_en);
        let d2 = dose_rate(2e6, e_j, mu_en);
        assert!(
            (d2 / d1 - 2.0).abs() < EPS,
            "Dose rate must scale linearly with fluence"
        );
    }

    // 13. Dose rate is zero for zero fluence.
    #[test]
    fn test_dose_rate_zero_fluence() {
        let d = dose_rate(0.0, 1.6e-13, 2.7e-3);
        assert_eq!(d, 0.0);
    }

    // 14. Shielded dose rate at zero thickness equals unshielded (buildup=1).
    #[test]
    fn test_shielded_dose_no_shield() {
        let d_unshielded = dose_rate(1e8, 1.6e-13, 2.7e-3);
        let d_shielded = shielded_dose_rate(1e8, 1.6e-13, 2.7e-3, 80.0, 0.0, 1.0);
        assert!((d_shielded - d_unshielded).abs() < EPS * d_unshielded);
    }

    // 15. Shielded dose rate decreases with increasing thickness.
    #[test]
    fn test_shielded_dose_decreases_with_thickness() {
        let (phi0, e, mu_en, mu) = (1e10, 1.6e-13, 2.7e-3, 80.0);
        let d1 = shielded_dose_rate(phi0, e, mu_en, mu, 0.01, 1.0);
        let d2 = shielded_dose_rate(phi0, e, mu_en, mu, 0.05, 1.0);
        assert!(d1 > d2, "Dose rate must decrease with shield thickness");
    }

    // 16. Activation analysis: saturation activity scales linearly with flux.
    #[test]
    fn test_activation_saturation_scales_with_flux() {
        let r1 = activation_analysis(1e-3, 0.059, 1e-28, 1e14, 3.7e3, 1e6, 0.0);
        let r2 = activation_analysis(1e-3, 0.059, 1e-28, 2e14, 3.7e3, 1e6, 0.0);
        let ratio = r2.saturation_activity / r1.saturation_activity;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "Saturation activity must scale with flux, ratio={ratio}"
        );
    }

    // 17. Activity at end of irradiation → saturation for very long irradiation.
    #[test]
    fn test_activation_approaches_saturation() {
        let half_life = 600.0; // 10 min
        // 100 half-lives of irradiation → should be ≈ saturation
        let result =
            activation_analysis(1e-3, 0.059, 1e-28, 1e14, half_life, 100.0 * half_life, 0.0);
        let ratio = result.activity_at_eoi / result.saturation_activity;
        assert!(
            (ratio - 1.0).abs() < 1e-4,
            "Long irradiation should approach saturation, ratio={ratio}"
        );
    }

    // 18. Activity decays to zero after many half-lives of cooling.
    #[test]
    fn test_activation_decay_after_cooling() {
        let half_life = 600.0;
        let result = activation_analysis(
            1e-3,
            0.059,
            1e-28,
            1e14,
            half_life,
            10.0 * half_life,
            100.0 * half_life,
        );
        assert!(
            result.activity_after_cooling < result.activity_at_eoi * 1e-25,
            "Activity should be negligible after 100 half-lives"
        );
    }

    // 19. Decay constant λ = ln(2) / t½
    #[test]
    fn test_decay_constant() {
        let half_life = 1234.5;
        let result = activation_analysis(1e-3, 0.059, 1e-28, 1e14, half_life, 1.0, 0.0);
        let expected_lambda = LN_2 / half_life;
        assert!(
            (result.decay_constant - expected_lambda).abs() < EPS,
            "Decay constant mismatch: {} vs {}",
            result.decay_constant,
            expected_lambda
        );
    }

    // 20. Required thickness: round-trip with transmission.
    #[test]
    fn test_required_thickness_round_trip() {
        let mu = 80.0;
        let target_reduction = 1000.0;
        let x = required_thickness(target_reduction, mu);
        let actual_reduction = 1.0 / narrow_beam_transmission(mu, x);
        assert!(
            (actual_reduction - target_reduction).abs() < 1e-6,
            "Round-trip failed: {actual_reduction} vs {target_reduction}"
        );
    }

    // 21. HVLs for transmission: 1 HVL → 0.5 transmission.
    #[test]
    fn test_hvls_for_half_transmission() {
        let n = hvls_for_transmission(0.5);
        assert!(
            (n - 1.0).abs() < EPS,
            "Need exactly 1 HVL for 50% transmission, got {n}"
        );
    }

    // 22. HVLs for transmission: 1/4 transmission → 2 HVLs.
    #[test]
    fn test_hvls_for_quarter_transmission() {
        let n = hvls_for_transmission(0.25);
        assert!(
            (n - 2.0).abs() < EPS,
            "Need 2 HVLs for 25% transmission, got {n}"
        );
    }

    // 23. Macroscopic removal cross-section is positive for valid inputs.
    #[test]
    fn test_macroscopic_removal_cross_section_positive() {
        // Water: ρ=1000, M≈0.018 kg/mol, σ_r≈3e-30 m²
        let sigma = macroscopic_removal_cross_section(1000.0, 0.018, 3e-30);
        assert!(
            sigma > 0.0,
            "Macroscopic cross-section must be positive, got {sigma}"
        );
    }

    // 24. Mean free path is reciprocal of macroscopic cross-section.
    #[test]
    fn test_mean_free_path_reciprocal() {
        let sigma_t = 100.0; // 1/m
        let mfp = mean_free_path(sigma_t);
        assert!(
            (mfp - 0.01).abs() < EPS,
            "MFP must be 1/Σ_t = 0.01 m, got {mfp}"
        );
    }

    // 25. Water preset has density 1000 kg/m³.
    #[test]
    fn test_water_preset_density() {
        let water = RadiationMaterial::water();
        assert_eq!(water.density, 1000.0);
    }

    // 26. Concrete preset has lower density than lead.
    #[test]
    fn test_density_ordering() {
        let lead = RadiationMaterial::lead();
        let concrete = RadiationMaterial::concrete();
        assert!(
            lead.density > concrete.density,
            "Lead must be denser than concrete"
        );
    }

    // 27. Effective dose rate scales linearly with activity.
    #[test]
    fn test_effective_dose_rate_linear() {
        let dcf = 1e-15; // Sv/Bq
        let d1 = effective_dose_rate(1e6, dcf);
        let d2 = effective_dose_rate(3e6, dcf);
        assert!(
            (d2 / d1 - 3.0).abs() < EPS,
            "Effective dose rate must scale with activity"
        );
    }

    // 28. Required thickness increases with higher reduction factors.
    #[test]
    fn test_required_thickness_monotone() {
        let mu = 50.0;
        let x10 = required_thickness(10.0, mu);
        let x100 = required_thickness(100.0, mu);
        let x1000 = required_thickness(1000.0, mu);
        assert!(
            x10 < x100 && x100 < x1000,
            "Required thickness must increase with reduction factor"
        );
    }

    // 29. Taylor buildup is always ≥ 1 for positive α₁ and α₂.
    #[test]
    fn test_taylor_buildup_gte_one_positive_alphas() {
        // With both exponents positive the buildup only grows with depth.
        let (cap_a, a1, a2) = (0.4, 0.1, 0.05);
        for &mu_x in &[0.0, 0.5, 1.0, 2.0, 5.0] {
            let b = buildup_factor(mu_x, cap_a, a1, a2);
            assert!(b >= 1.0, "Buildup must be ≥ 1 at μx={mu_x}, got {b}");
        }
    }

    // 30. Polyethylene has lower density than concrete.
    #[test]
    fn test_polyethylene_density() {
        let pe = RadiationMaterial::polyethylene();
        let concrete = RadiationMaterial::concrete();
        assert!(
            pe.density < concrete.density,
            "Polyethylene must be lighter than concrete"
        );
    }
}
