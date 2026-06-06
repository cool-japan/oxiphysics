// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polymer mechanical models: worm-like chain, freely-jointed chain, viscoelasticity,
//! Prony-series relaxation, glass transition, Flory-Huggins mixing, rubber elasticity,
//! crystallinity effects, and crazing.
//!
//! All functions use SI units unless otherwise stated.

// ---------------------------------------------------------------------------
// PolymerChain
// ---------------------------------------------------------------------------

/// Geometric and flexibility parameters of a single polymer chain.
#[derive(Debug, Clone)]
pub struct PolymerChain {
    /// Number of Kuhn segments (or equivalent rigid units).
    pub n_segments: usize,
    /// Kuhn / statistical segment length (m).
    pub segment_length: f64,
    /// Persistence length (m) — half of Kuhn length for a worm-like chain.
    pub persistence_length: f64,
    /// Contour length L_c (m) = n_segments × segment_length.
    pub contour_length: f64,
}

impl PolymerChain {
    /// Construct a [`PolymerChain`] and derive `contour_length` automatically.
    pub fn new(n_segments: usize, segment_length: f64, persistence_length: f64) -> Self {
        Self {
            n_segments,
            segment_length,
            persistence_length,
            contour_length: n_segments as f64 * segment_length,
        }
    }

    /// End-to-end distance root-mean-square for an ideal (Gaussian) chain (m).
    pub fn rms_end_to_end(&self) -> f64 {
        (self.n_segments as f64).sqrt() * self.segment_length
    }
}

// ---------------------------------------------------------------------------
// ViscoelasticPolymer
// ---------------------------------------------------------------------------

/// Linear viscoelastic polymer described by a Prony series and basic moduli.
#[derive(Debug, Clone)]
pub struct ViscoelasticPolymer {
    /// Instantaneous (glassy) elastic modulus E₀ (Pa).
    pub elastic_modulus: f64,
    /// Zero-shear viscosity η₀ (Pa·s).
    pub viscosity: f64,
    /// Prony relaxation times τᵢ (s).
    pub relaxation_times: Vec<f64>,
    /// Prony partial moduli Eᵢ (Pa) — must match length of `relaxation_times`.
    pub moduli: Vec<f64>,
}

impl ViscoelasticPolymer {
    /// Create a new [`ViscoelasticPolymer`].
    pub fn new(
        elastic_modulus: f64,
        viscosity: f64,
        relaxation_times: Vec<f64>,
        moduli: Vec<f64>,
    ) -> Self {
        Self {
            elastic_modulus,
            viscosity,
            relaxation_times,
            moduli,
        }
    }
}

// ---------------------------------------------------------------------------
// Worm-like chain (WLC) model
// ---------------------------------------------------------------------------

/// Marko-Siggia worm-like chain force–extension relation (interpolation formula).
///
/// `F = (k_B T / L_p) · [ 1/(4(1−x)²) − 1/4 + x ]`
///
/// where x = extension / L_c.
///
/// # Arguments
/// * `extension` – end-to-end extension r (m), must be < `contour`.
/// * `contour`   – contour length L_c (m).
/// * `persistence` – persistence length L_p (m).
/// * `kb_temp`   – thermal energy k_B T (J).
///
/// Returns the restoring force (N).  Returns a large value as `extension → contour`.
pub fn worm_like_chain_force(extension: f64, contour: f64, persistence: f64, kb_temp: f64) -> f64 {
    let x = (extension / contour).min(1.0 - 1e-9);
    let x = x.max(0.0);
    let prefactor = kb_temp / persistence;
    prefactor * (0.25 / (1.0 - x).powi(2) - 0.25 + x)
}

// ---------------------------------------------------------------------------
// Freely-jointed chain (FJC) model
// ---------------------------------------------------------------------------

/// Freely-jointed chain (Langevin / inverse-Langevin) force–extension.
///
/// Uses the Cohen (1991) Padé approximation to the inverse Langevin function:
///
/// `β = x(3 − x²)/(1 − x²)`   (Padé)
///
/// `F = (k_B T / b) · β(r / N b)`
///
/// where b is the Kuhn length and N is the number of segments.
///
/// # Arguments
/// * `extension` – end-to-end distance r (m).
/// * `contour`   – contour length L_c = N b (m).
/// * `kb_temp`   – thermal energy k_B T (J).
///
/// Returns force (N).
pub fn freely_jointed_chain_force(extension: f64, contour: f64, kb_temp: f64) -> f64 {
    // b = segment length derived from contour, assume n_segments = 1 for this call.
    // The caller should pass contour = N * b; b appears via kb_temp / b.
    // We treat b = contour / N ≈ segment scale.  For a single-segment call:
    let x = (extension / contour).clamp(-1.0 + 1e-9, 1.0 - 1e-9);
    // Padé approximation to inverse Langevin
    let beta = x * (3.0 - x * x) / (1.0 - x * x);
    // Force = (k_B T / b_kuhn) * beta, b_kuhn = contour (single segment)
    (kb_temp / contour) * beta
}

// ---------------------------------------------------------------------------
// Entropic spring constant
// ---------------------------------------------------------------------------

/// Effective spring constant of a Gaussian polymer chain (entropic spring).
///
/// `k = 3 k_B T / (N b²)`  — valid for small extensions r << √N b.
///
/// # Arguments
/// * `chain`   – polymer chain parameters.
/// * `kb_temp` – thermal energy k_B T (J).
///
/// Returns spring constant (N/m).
pub fn entropic_spring_constant(chain: &PolymerChain, kb_temp: f64) -> f64 {
    let n = chain.n_segments as f64;
    let b = chain.segment_length;
    3.0 * kb_temp / (n * b * b)
}

// ---------------------------------------------------------------------------
// Prony series stress relaxation
// ---------------------------------------------------------------------------

/// Time-domain relaxation modulus from a generalised Maxwell (Prony series) model.
///
/// `E(t) = Σᵢ Eᵢ · exp(−t / τᵢ)`
///
/// The equilibrium modulus is excluded (spring-pot limit); add it separately if needed.
///
/// # Arguments
/// * `polymer` – viscoelastic material with Prony parameters.
/// * `t`       – time (s).
///
/// Returns E(t) (Pa).
pub fn prony_series_relaxation(polymer: &ViscoelasticPolymer, t: f64) -> f64 {
    polymer
        .moduli
        .iter()
        .zip(polymer.relaxation_times.iter())
        .map(|(&e_i, &tau_i)| e_i * (-t / tau_i).exp())
        .sum()
}

// ---------------------------------------------------------------------------
// Glass-transition temperature (Fox / Flory equation for MW dependence)
// ---------------------------------------------------------------------------

/// Fox-Flory molecular-weight dependence of the glass-transition temperature.
///
/// `T_g = T_g∞ − K_g / M_w`
///
/// # Arguments
/// * `molecular_weight` – number-average molecular weight M_w (g/mol).
/// * `tg_inf`          – bulk T_g at infinite MW (K).
/// * `kg`              – Fox-Flory constant (K·g/mol), typically 1×10⁵ to 2×10⁶.
///
/// Returns T_g (K).
pub fn glass_transition_temp(molecular_weight: f64, tg_inf: f64, kg: f64) -> f64 {
    tg_inf - kg / molecular_weight
}

// ---------------------------------------------------------------------------
// Flory-Huggins interaction parameter χ
// ---------------------------------------------------------------------------

/// Estimate the Flory-Huggins interaction parameter χ from solubility parameters.
///
/// `χ = (V_r / k_B T) · (δ_polymer − δ_solvent)²`
///
/// # Arguments
/// * `solubility_polymer`  – solubility parameter of polymer δ_p (Pa^0.5).
/// * `solubility_solvent`  – solubility parameter of solvent δ_s (Pa^0.5).
/// * `molar_vol`           – molar volume of solvent (m³/mol).
/// * `kb_temp`             – thermal energy k_B T per molecule (J).
///
/// Returns dimensionless χ.
pub fn flory_huggins_chi(
    solubility_polymer: f64,
    solubility_solvent: f64,
    molar_vol: f64,
    kb_temp: f64,
) -> f64 {
    let delta_sq = (solubility_polymer - solubility_solvent).powi(2);
    molar_vol * delta_sq / kb_temp
}

// ---------------------------------------------------------------------------
// Crystallinity effect on modulus (Rule of mixtures / Takayanagi)
// ---------------------------------------------------------------------------

/// Simple rule-of-mixtures estimate for the modulus of a semi-crystalline polymer.
///
/// `E = (1 − φ) · E_a + φ · E_c`
///
/// where φ is the degree of crystallinity.
///
/// # Arguments
/// * `e_amorphous`    – modulus of the amorphous phase (Pa).
/// * `e_crystalline`  – modulus of the crystalline phase (Pa).
/// * `crystallinity`  – volume fraction of crystalline phase φ (0–1).
///
/// Returns effective modulus (Pa).
pub fn crystallinity_effect_on_modulus(
    e_amorphous: f64,
    e_crystalline: f64,
    crystallinity: f64,
) -> f64 {
    let phi = crystallinity.clamp(0.0, 1.0);
    (1.0 - phi) * e_amorphous + phi * e_crystalline
}

// ---------------------------------------------------------------------------
// Crazing stress (Dugdale / empirical)
// ---------------------------------------------------------------------------

/// Estimate crazing initiation stress using a power-law rate-dependent model.
///
/// `σ_craze = E · (ε̇ · τ)^(1/n)`
///
/// where n is the Eyring flow exponent (≈ viscosity / elastic_modulus) and
/// τ = viscosity / elastic_modulus is the Maxwell relaxation time.
///
/// # Arguments
/// * `polymer`     – viscoelastic material parameters.
/// * `strain_rate` – applied strain rate (1/s).
///
/// Returns crazing stress (Pa).
pub fn crazing_stress(polymer: &ViscoelasticPolymer, strain_rate: f64) -> f64 {
    // Maxwell relaxation time τ = η / E
    let tau = polymer.viscosity / polymer.elastic_modulus;
    // Empirical exponent n = 10 (typical for glassy polymers)
    let n = 10.0_f64;
    polymer.elastic_modulus * (strain_rate * tau).powf(1.0 / n)
}

// ---------------------------------------------------------------------------
// Rubber elastic modulus
// ---------------------------------------------------------------------------

/// Affine network model rubber elastic shear modulus.
///
/// `G = n_c · k_B · T`
///
/// where n_c is the number density of network chains (m⁻³).
///
/// # Arguments
/// * `n_chains` – network chain density n_c (chains/m³).
/// * `kb_temp`  – thermal energy k_B T (J).
///
/// Returns shear modulus G (Pa).
pub fn rubber_elastic_modulus(n_chains: f64, kb_temp: f64) -> f64 {
    n_chains * kb_temp
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;
    const KB_ROOM: f64 = 4.11e-21; // k_B * 300 K ≈ 4.14e-21 J

    // --- PolymerChain ---

    #[test]
    fn test_polymer_chain_contour() {
        let c = PolymerChain::new(100, 1e-9, 0.5e-9);
        assert!((c.contour_length - 100.0e-9).abs() < EPS, "contour = N*b");
    }

    #[test]
    fn test_polymer_chain_rms() {
        let c = PolymerChain::new(100, 1e-9, 0.5e-9);
        let expected = 10.0e-9;
        assert!(
            (c.rms_end_to_end() - expected).abs() < 1e-20,
            "rms = √N * b"
        );
    }

    #[test]
    fn test_polymer_chain_n_zero_contour() {
        let c = PolymerChain::new(0, 1e-9, 0.5e-9);
        assert_eq!(c.contour_length, 0.0);
    }

    // --- worm_like_chain_force ---

    #[test]
    fn test_wlc_force_positive() {
        let f = worm_like_chain_force(20e-9, 30e-9, 50e-9, KB_ROOM);
        assert!(f > 0.0, "WLC force must be positive");
    }

    #[test]
    fn test_wlc_force_increases_with_extension() {
        let f1 = worm_like_chain_force(10e-9, 30e-9, 50e-9, KB_ROOM);
        let f2 = worm_like_chain_force(25e-9, 30e-9, 50e-9, KB_ROOM);
        assert!(f2 > f1, "WLC force increases with extension");
    }

    #[test]
    fn test_wlc_force_near_zero_extension() {
        // At very small extension x << 1, force ~ k_BT/Lp * x (linear)
        let f = worm_like_chain_force(1e-12, 30e-9, 50e-9, KB_ROOM);
        assert!(f >= 0.0, "WLC force non-negative");
    }

    #[test]
    fn test_wlc_force_scales_with_temp() {
        let f1 = worm_like_chain_force(20e-9, 30e-9, 50e-9, KB_ROOM);
        let f2 = worm_like_chain_force(20e-9, 30e-9, 50e-9, 2.0 * KB_ROOM);
        assert!((f2 - 2.0 * f1).abs() < 1e-25, "WLC linear in kBT");
    }

    // --- freely_jointed_chain_force ---

    #[test]
    fn test_fjc_force_positive_extension() {
        let f = freely_jointed_chain_force(5e-9, 30e-9, KB_ROOM);
        assert!(f > 0.0, "FJC force positive for positive extension");
    }

    #[test]
    fn test_fjc_force_negative_extension() {
        let f = freely_jointed_chain_force(-5e-9, 30e-9, KB_ROOM);
        assert!(f < 0.0, "FJC force negative for compression");
    }

    #[test]
    fn test_fjc_force_zero_at_zero() {
        let f = freely_jointed_chain_force(0.0, 30e-9, KB_ROOM);
        assert!(f.abs() < EPS, "FJC force zero at zero extension");
    }

    #[test]
    fn test_fjc_force_increases() {
        let f1 = freely_jointed_chain_force(5e-9, 30e-9, KB_ROOM);
        let f2 = freely_jointed_chain_force(15e-9, 30e-9, KB_ROOM);
        assert!(f2 > f1, "FJC force increases with extension");
    }

    // --- entropic_spring_constant ---

    #[test]
    fn test_entropic_spring_positive() {
        let c = PolymerChain::new(100, 1e-9, 0.5e-9);
        let k = entropic_spring_constant(&c, KB_ROOM);
        assert!(k > 0.0, "spring constant must be positive");
    }

    #[test]
    fn test_entropic_spring_decreases_with_n() {
        let c1 = PolymerChain::new(10, 1e-9, 0.5e-9);
        let c2 = PolymerChain::new(100, 1e-9, 0.5e-9);
        let k1 = entropic_spring_constant(&c1, KB_ROOM);
        let k2 = entropic_spring_constant(&c2, KB_ROOM);
        assert!(k1 > k2, "more segments → softer spring");
    }

    #[test]
    fn test_entropic_spring_formula() {
        let c = PolymerChain::new(50, 2e-9, 1e-9);
        let k = entropic_spring_constant(&c, KB_ROOM);
        let expected = 3.0 * KB_ROOM / (50.0 * 4.0e-18);
        assert!((k - expected).abs() < EPS, "spring formula");
    }

    // --- prony_series_relaxation ---

    #[test]
    fn test_prony_at_t0() {
        let p = ViscoelasticPolymer::new(1e9, 1e6, vec![1.0, 10.0], vec![0.5e9, 0.5e9]);
        let e = prony_series_relaxation(&p, 0.0);
        assert!((e - 1.0e9).abs() < 1.0, "E(0) = sum of Prony moduli");
    }

    #[test]
    fn test_prony_relaxes_to_zero() {
        let p = ViscoelasticPolymer::new(1e9, 1e6, vec![0.001], vec![1e9]);
        let e = prony_series_relaxation(&p, 1000.0);
        assert!(e < 1.0, "E(t→∞) → 0 for single-arm Prony");
    }

    #[test]
    fn test_prony_monotone_decreasing() {
        let p = ViscoelasticPolymer::new(1e9, 1e6, vec![1.0, 5.0], vec![0.4e9, 0.6e9]);
        let e1 = prony_series_relaxation(&p, 1.0);
        let e2 = prony_series_relaxation(&p, 10.0);
        assert!(e2 < e1, "relaxation must decrease with time");
    }

    #[test]
    fn test_prony_empty_series() {
        let p = ViscoelasticPolymer::new(1e9, 1e6, vec![], vec![]);
        let e = prony_series_relaxation(&p, 5.0);
        assert_eq!(e, 0.0, "empty Prony series → zero");
    }

    // --- glass_transition_temp ---

    #[test]
    fn test_tg_approaches_tg_inf() {
        // Very high MW → T_g ≈ T_g∞
        let tg = glass_transition_temp(1e7, 373.0, 1e5);
        assert!((tg - 373.0).abs() < 0.1, "high MW → T_g ≈ T_g∞");
    }

    #[test]
    fn test_tg_increases_with_mw() {
        let tg1 = glass_transition_temp(10_000.0, 373.0, 1e5);
        let tg2 = glass_transition_temp(100_000.0, 373.0, 1e5);
        assert!(tg2 > tg1, "T_g increases with MW");
    }

    #[test]
    fn test_tg_known_formula() {
        let tg = glass_transition_temp(50_000.0, 400.0, 2e6);
        let expected = 400.0 - 2e6 / 50_000.0;
        assert!((tg - expected).abs() < EPS, "Fox-Flory formula");
    }

    // --- flory_huggins_chi ---

    #[test]
    fn test_chi_zero_delta_delta() {
        // Identical solubility parameters → χ = 0 (ideal mixing)
        let chi = flory_huggins_chi(20_000.0, 20_000.0, 1e-4, KB_ROOM);
        assert!(chi.abs() < EPS, "identical δ → χ = 0");
    }

    #[test]
    fn test_chi_positive_for_different_params() {
        let chi = flory_huggins_chi(20_000.0, 18_000.0, 1e-4, KB_ROOM);
        assert!(chi > 0.0, "χ > 0 for different solubility parameters");
    }

    #[test]
    fn test_chi_scales_with_molar_vol() {
        let chi1 = flory_huggins_chi(20_000.0, 18_000.0, 1e-4, KB_ROOM);
        let chi2 = flory_huggins_chi(20_000.0, 18_000.0, 2e-4, KB_ROOM);
        assert!((chi2 - 2.0 * chi1).abs() < EPS, "χ linear in V_r");
    }

    // --- crystallinity_effect_on_modulus ---

    #[test]
    fn test_crystallinity_zero_gives_amorphous() {
        let e = crystallinity_effect_on_modulus(1e8, 5e9, 0.0);
        assert!((e - 1e8).abs() < EPS, "0% crystallinity → E_amorphous");
    }

    #[test]
    fn test_crystallinity_one_gives_crystalline() {
        let e = crystallinity_effect_on_modulus(1e8, 5e9, 1.0);
        assert!((e - 5e9).abs() < EPS, "100% crystallinity → E_crystalline");
    }

    #[test]
    fn test_crystallinity_intermediate() {
        let e = crystallinity_effect_on_modulus(1e8, 5e9, 0.5);
        let expected = 0.5 * 1e8 + 0.5 * 5e9;
        assert!((e - expected).abs() < EPS, "50% crystallinity");
    }

    #[test]
    fn test_crystallinity_increases_modulus() {
        let e1 = crystallinity_effect_on_modulus(1e8, 5e9, 0.3);
        let e2 = crystallinity_effect_on_modulus(1e8, 5e9, 0.7);
        assert!(e2 > e1, "more crystalline → stiffer");
    }

    // --- crazing_stress ---

    #[test]
    fn test_crazing_stress_positive() {
        let p = ViscoelasticPolymer::new(3e9, 1e12, vec![1.0], vec![1e9]);
        let s = crazing_stress(&p, 1.0);
        assert!(s > 0.0, "crazing stress must be positive");
    }

    #[test]
    fn test_crazing_stress_increases_with_rate() {
        let p = ViscoelasticPolymer::new(3e9, 1e12, vec![1.0], vec![1e9]);
        let s1 = crazing_stress(&p, 1.0);
        let s2 = crazing_stress(&p, 100.0);
        assert!(s2 > s1, "higher strain rate → higher crazing stress");
    }

    // --- rubber_elastic_modulus ---

    #[test]
    fn test_rubber_modulus_positive() {
        let g = rubber_elastic_modulus(1e26, KB_ROOM);
        assert!(g > 0.0, "rubber modulus must be positive");
    }

    #[test]
    fn test_rubber_modulus_formula() {
        let n = 1e26_f64;
        let g = rubber_elastic_modulus(n, KB_ROOM);
        assert!((g - n * KB_ROOM).abs() < EPS, "G = n k_B T");
    }

    #[test]
    fn test_rubber_modulus_doubles_with_density() {
        let g1 = rubber_elastic_modulus(1e26, KB_ROOM);
        let g2 = rubber_elastic_modulus(2e26, KB_ROOM);
        assert!((g2 - 2.0 * g1).abs() < EPS, "G linear in n");
    }

    #[test]
    fn test_rubber_modulus_mpa_range() {
        // Typical crosslinked rubber: n ~ 1e26 /m³, G ~ 0.4 MPa
        let g = rubber_elastic_modulus(1e26, KB_ROOM);
        // k_BT at 300K ≈ 4.14e-21 J → G ≈ 0.4 MPa
        assert!(g > 1e4 && g < 1e7, "G in reasonable rubber range");
    }
}
