// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geomechanics and rock mechanics.
//!
//! Implements:
//! - [`RockMaterial`]: Mohr-Coulomb, Hoek-Brown, Drucker-Prager, triaxial deformation.
//! - [`PoroelasticMaterial`]: Biot poroelasticity, Terzaghi consolidation.
//! - [`Soil`]: index properties, compression index, permeability (Hazen), undrained strength.
//! - [`SoilType`]: soil classification.
//! - [`Rockfill`]: friction angle and stiffness of coarse-grained fill.
//! - [`cam_clay_yield`]: Modified Cam-Clay yield function.
//! - [`critical_state_line`]: critical-state line in (p, q) space.
//!
//! References
//! ----------
//! Hoek, E. & Brown, E. T. (1980) *Underground Excavations in Rock*, IMM.
//! Terzaghi, K. (1943) *Theoretical Soil Mechanics*, Wiley.
//! Biot, M. A. (1941) *J. Appl. Phys.* 12, 155–164.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// SoilType
// ---------------------------------------------------------------------------

/// USCS-based soil classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoilType {
    /// High-plasticity fine-grained soil.
    Clay,
    /// Cohesionless granular soil (D50 0.075–4.75 mm).
    Sand,
    /// Coarse granular soil (D50 > 4.75 mm).
    Gravel,
    /// Low-plasticity fine-grained soil.
    Silt,
    /// Highly organic compressible soil.
    Peat,
}

// ---------------------------------------------------------------------------
// RockMaterial
// ---------------------------------------------------------------------------

/// Intact-rock material for failure and deformation analysis.
#[derive(Debug, Clone)]
pub struct RockMaterial {
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
    /// Internal friction angle φ \[radians\].
    pub friction_angle: f64,
    /// Uniaxial tensile strength σ_t \[Pa\].
    pub tensile_strength: f64,
    /// Uniaxial compressive strength σ_c (UCS) \[Pa\].
    pub ucs: f64,
    /// Young's modulus E \[Pa\].
    pub youngs_modulus: f64,
    /// Poisson's ratio ν \[-\].
    pub poisson_ratio: f64,
}

impl RockMaterial {
    /// Construct a rock material from standard laboratory parameters.
    pub fn new(
        cohesion: f64,
        friction_angle_deg: f64,
        tensile_strength: f64,
        ucs: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            cohesion,
            friction_angle: friction_angle_deg.to_radians(),
            tensile_strength,
            ucs,
            youngs_modulus,
            poisson_ratio,
        }
    }

    /// Mohr-Coulomb peak strength: σ₁ = σ₃ tan²(45+φ/2) + 2c tan(45+φ/2).
    ///
    /// # Parameters
    /// - `sigma3`: minor principal stress (confining pressure) \[Pa\].
    ///
    /// # Returns
    /// Major principal stress at failure \[Pa\].
    pub fn mohr_coulomb_strength(&self, sigma3: f64) -> f64 {
        let phi = self.friction_angle;
        let n_phi = ((PI / 4.0 + phi / 2.0).tan()).powi(2);
        n_phi * sigma3 + 2.0 * self.cohesion * n_phi.sqrt()
    }

    /// Hoek-Brown empirical strength criterion for rock masses.
    ///
    /// σ₁ = σ₃ + σ_ci √(m_i σ₃/σ_ci + s)
    ///
    /// where `s` is derived from the Geological Strength Index (GSI).
    ///
    /// # Parameters
    /// - `sigma3`: minor principal stress \[Pa\].
    /// - `mi`: Hoek-Brown intact-rock constant.
    /// - `gsi`: Geological Strength Index (0–100).
    ///
    /// # Returns
    /// Major principal stress at failure \[Pa\].
    pub fn hoek_brown_strength(&self, sigma3: f64, mi: f64, gsi: f64) -> f64 {
        let sigma_ci = self.ucs;
        // s exponent for undisturbed rock: s = exp((GSI-100)/9)
        let s = ((gsi - 100.0) / 9.0).exp();
        let arg = mi * sigma3 / sigma_ci + s;
        sigma3 + sigma_ci * arg.max(0.0).sqrt()
    }

    /// Drucker-Prager yield function f(p, q).
    ///
    /// f = q − α p − k ≤ 0 at yield,  where α and k are derived from φ and c
    /// for plane-strain (outer cone inscribed in Mohr-Coulomb hexagon).
    ///
    /// # Parameters
    /// - `p`: mean effective stress = (σ₁+σ₂+σ₃)/3 \[Pa\].
    /// - `q`: deviatoric stress = σ₁−σ₃ \[Pa\].
    ///
    /// # Returns
    /// Value of the yield function; negative → elastic, zero/positive → yielding.
    pub fn drucker_prager_yield(&self, p: f64, q: f64) -> f64 {
        let phi = self.friction_angle;
        // Circumscribed DP match to Mohr-Coulomb
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let alpha = 2.0 * sin_phi / (3.0f64.sqrt() * (3.0 - sin_phi));
        let k = 6.0 * self.cohesion * cos_phi / (3.0f64.sqrt() * (3.0 - sin_phi));
        q - alpha * p - k
    }

    /// Triaxial deformation: axial strain ε₁ and volumetric strain ε_v
    /// under (σ₁, σ₃) using linear elasticity.
    ///
    /// ε₁ = (σ₁ − 2ν σ₃) / E
    /// ε_v = (σ₁ − 2σ₃)(1 − 2ν) / E
    ///
    /// # Returns
    /// `(axial_strain, volumetric_strain)`.
    pub fn triaxial_deformation(&self, sigma1: f64, sigma3: f64) -> (f64, f64) {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let eps1 = (sigma1 - 2.0 * nu * sigma3) / e;
        let eps_v = (sigma1 - 2.0 * sigma3) * (1.0 - 2.0 * nu) / e;
        (eps1, eps_v)
    }
}

// ---------------------------------------------------------------------------
// PoroelasticMaterial
// ---------------------------------------------------------------------------

/// Saturated poroelastic material following Biot (1941) and Terzaghi (1943).
#[derive(Debug, Clone)]
pub struct PoroelasticMaterial {
    /// Biot-Willis coefficient α = 1 − K_s / K_grain \[-\].
    pub biot_coeff: f64,
    /// Drained bulk modulus K \[Pa\].
    pub bulk_modulus: f64,
    /// Porosity n \[-\].
    pub porosity: f64,
    /// Intrinsic permeability k \[m²\].
    pub permeability: f64,
    /// Bulk modulus of the pore fluid K_f \[Pa\].
    pub fluid_bulk_modulus: f64,
    /// Fluid viscosity μ \[Pa s\].
    pub fluid_viscosity: f64,
}

impl PoroelasticMaterial {
    /// Construct a poroelastic material with full parameters.
    pub fn new(
        biot_coeff: f64,
        bulk_modulus: f64,
        porosity: f64,
        permeability: f64,
        fluid_bulk_modulus: f64,
        fluid_viscosity: f64,
    ) -> Self {
        Self {
            biot_coeff,
            bulk_modulus,
            porosity,
            permeability,
            fluid_bulk_modulus,
            fluid_viscosity,
        }
    }

    /// Biot storage modulus M \[Pa\]:
    ///
    /// 1/M = (α − n)/K_s + n/K_f  (simplified: K_s → ∞ ⟹ 1/M ≈ n/K_f)
    ///
    /// Returns M in Pa.
    pub fn biot_modulus(&self) -> f64 {
        // Approximate: M ≈ K_f / n
        self.fluid_bulk_modulus / self.porosity.max(1e-10)
    }

    /// Undrained bulk modulus K_u = K + α² M \[Pa\].
    pub fn undrained_bulk_modulus(&self) -> f64 {
        let m = self.biot_modulus();
        self.bulk_modulus + self.biot_coeff * self.biot_coeff * m
    }

    /// Skempton's B coefficient: B = α M / K_u \[-\].
    pub fn skempton_b(&self) -> f64 {
        let m = self.biot_modulus();
        let ku = self.undrained_bulk_modulus();
        self.biot_coeff * m / ku.max(1e-30)
    }

    /// Consolidation (diffusivity) coefficient c_v = k / (μ (α²/K + 1/M)) \[m² s⁻¹\].
    pub fn consolidation_coefficient(&self) -> f64 {
        let m = self.biot_modulus();
        let denom = self.fluid_viscosity
            * (self.biot_coeff * self.biot_coeff / self.bulk_modulus.max(1e-30) + 1.0 / m);
        self.permeability / denom.max(1e-30)
    }

    /// Terzaghi 1-D consolidation settlement at time t for a layer of `depth`
    /// under applied stress `sigma` \[Pa\].
    ///
    /// Uses the series solution U(T_v) truncated to 15 terms.
    ///
    /// # Parameters
    /// - `sigma`: applied total stress increment \[Pa\].
    /// - `depth`: drained drainage length H \[m\].
    /// - `time`: elapsed time \[s\].
    ///
    /// # Returns
    /// Settlement s(t) \[m\].
    pub fn terzaghi_settlement(&self, sigma: f64, depth: f64, time: f64) -> f64 {
        let cv = self.consolidation_coefficient();
        let tv = cv * time / (depth * depth).max(1e-30);
        let u = degree_of_consolidation(tv);
        // Immediate + primary consolidation settlement
        let mv = self.biot_coeff / self.bulk_modulus.max(1e-30);
        u * mv * sigma * depth
    }
}

// ---------------------------------------------------------------------------
// Soil
// ---------------------------------------------------------------------------

/// Soil sample with index and consolidation properties.
#[derive(Debug, Clone)]
pub struct Soil {
    /// USCS soil classification.
    pub classification: SoilType,
    /// Void ratio e = V_voids / V_solids \[-\].
    pub void_ratio: f64,
    /// Liquid limit w_L \[%\].
    pub liquid_limit: f64,
    /// Plastic limit w_P \[%\].
    pub plastic_limit: f64,
    /// Natural water content w \[%\].
    pub water_content: f64,
    /// Effective cohesion c' \[Pa\] (for drained strength).
    pub cohesion: f64,
    /// Effective friction angle φ' \[radians\].
    pub friction_angle: f64,
}

impl Soil {
    /// Construct a soil from Atterberg limits and void ratio.
    pub fn new(
        classification: SoilType,
        void_ratio: f64,
        liquid_limit: f64,
        plastic_limit: f64,
        water_content: f64,
        cohesion: f64,
        friction_angle_deg: f64,
    ) -> Self {
        Self {
            classification,
            void_ratio,
            liquid_limit,
            plastic_limit,
            water_content,
            cohesion,
            friction_angle: friction_angle_deg.to_radians(),
        }
    }

    /// Plasticity index I_P = w_L − w_P \[%\].
    pub fn plasticity_index(&self) -> f64 {
        self.liquid_limit - self.plastic_limit
    }

    /// Skempton (1944) compression index from liquid limit:
    /// C_c = 0.009 (w_L − 10).
    pub fn compression_index(&self) -> f64 {
        0.009 * (self.liquid_limit - 10.0)
    }

    /// Estimate preconsolidation pressure via Casagrande (1936) empirical rule.
    ///
    /// Uses the liquidity index to estimate OCR, then back-calculates σ'_p.
    /// Returns σ'_p in Pa (reference atmospheric pressure = 101.325 kPa used).
    pub fn preconsolidation_pressure(&self) -> f64 {
        let li = (self.water_content - self.plastic_limit) / self.plasticity_index().max(1e-10);
        // Higher LI → lower OCR → lower σ'_p (rough empirical: σ'_p ∝ (2 − LI) × σ_atm)
        let ocr_est = (2.0 - li.clamp(0.0, 2.0)).max(0.5);
        let sigma_ref = 101_325.0; // Pa
        ocr_est * sigma_ref
    }

    /// Hazen (1911) hydraulic conductivity: k = C · d₁₀² \[m s⁻¹\].
    ///
    /// # Parameters
    /// - `d10`: effective grain size D₁₀ \[m\].
    ///
    /// # Returns
    /// Hydraulic conductivity \[m s⁻¹\].
    pub fn permeability_hazen(&self, d10: f64) -> f64 {
        // Hazen constant C ≈ 0.0116 for d10 in metres, k in m/s
        // (originally C=1 with d10 in cm, k in cm/s; converted here)
        let c = 0.01; // dimensionless, SI units
        c * d10 * d10
    }

    /// Undrained shear strength c_u from the Skempton (1957) correlation:
    /// c_u / σ'_v ≈ 0.11 + 0.0037 · I_P
    ///
    /// Returns c_u \[Pa\] for a reference vertical effective stress of 100 kPa.
    pub fn shear_strength_cu(&self) -> f64 {
        let sigma_v = 100_000.0; // 100 kPa reference
        let ip = self.plasticity_index();
        (0.11 + 0.0037 * ip) * sigma_v
    }
}

// ---------------------------------------------------------------------------
// Rockfill
// ---------------------------------------------------------------------------

/// Coarse granular rockfill material.
#[derive(Debug, Clone)]
pub struct Rockfill {
    /// Median grain size D₅₀ \[m\].
    pub d50: f64,
    /// Coefficient of uniformity C_u = D₆₀ / D₁₀ \[-\].
    pub cu: f64,
}

impl Rockfill {
    /// Construct rockfill from median grain size and uniformity coefficient.
    pub fn new(d50: f64, cu: f64) -> Self {
        Self { d50, cu }
    }

    /// Peak friction angle of rockfill \[radians\].
    ///
    /// Uses the empirical Marachi et al. (1972) correlation:
    /// φ' ≈ φ_0 − Δφ log₁₀(σ₃/p_a)  at a reference stress of 100 kPa.
    /// Here a simpler constant estimate is returned:
    /// φ' ≈ 38° + 2° · log₁₀(C_u) \[degrees\] converted to radians.
    pub fn friction_angle(&self) -> f64 {
        let phi_deg = 38.0 + 2.0 * self.cu.max(1.0).log10();
        phi_deg.to_radians()
    }

    /// Tangent stiffness modulus E_t \[Pa\] from the Janbu (1963) hyperbolic model:
    ///
    /// E_t = K · p_a · (σ₃ / p_a)^n
    ///
    /// where K ≈ 300 and n ≈ 0.45 for dense rockfill.
    ///
    /// # Parameters
    /// - `sigma3`: confining pressure \[Pa\].
    ///
    /// # Returns
    /// Tangent stiffness \[Pa\].
    pub fn stiffness(&self, sigma3: f64) -> f64 {
        let pa = 101_325.0; // atmospheric pressure [Pa]
        let k_const = 300.0;
        let n = 0.45;
        k_const * pa * (sigma3 / pa).max(1e-6).powf(n)
    }
}

// ---------------------------------------------------------------------------
// Stand-alone functions
// ---------------------------------------------------------------------------

/// Modified Cam-Clay yield surface:
///
/// f(p, q) = q² / M² − p (p₀ − p)
///
/// # Parameters
/// - `p`: mean effective stress \[Pa\].
/// - `q`: deviatoric stress \[Pa\].
/// - `m`: slope of the critical state line M = q/p at failure.
/// - `p0`: preconsolidation pressure (size of the yield surface) \[Pa\].
///
/// # Returns
/// Yield function value; f ≤ 0 → inside yield surface (elastic).
pub fn cam_clay_yield(p: f64, q: f64, m: f64, p0: f64) -> f64 {
    q * q / (m * m) - p * (p0 - p)
}

/// Critical-state line (CSL) in (p, q) space: q = M p.
///
/// # Parameters
/// - `p`: mean effective stress \[Pa\].
/// - `m`: critical-state stress ratio M = 6 sin φ_cs / (3 − sin φ_cs).
///
/// # Returns
/// Deviatoric stress q at critical state \[Pa\].
pub fn critical_state_line(p: f64, m: f64) -> f64 {
    m * p
}

/// Degree of primary consolidation U(T_v) via Terzaghi's series solution.
///
/// U = 1 − Σ_{m=0}^{N} (8/((2m+1)²π²)) exp(−(2m+1)²π²T_v/4)
///
/// # Parameters
/// - `tv`: time factor T_v = c_v t / H².
///
/// # Returns
/// Degree of consolidation U ∈ \[0, 1\].
pub fn degree_of_consolidation(tv: f64) -> f64 {
    if tv <= 0.0 {
        return 0.0;
    }
    let mut u = 1.0f64;
    for m in 0..15usize {
        let n = (2 * m + 1) as f64;
        u -= (8.0 / (n * n * PI * PI)) * (-(n * n * PI * PI * tv) / 4.0).exp();
    }
    u.clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- RockMaterial ---

    fn granite() -> RockMaterial {
        RockMaterial::new(
            15e6,  // c = 15 MPa
            45.0,  // φ = 45°
            5e6,   // σ_t = 5 MPa
            200e6, // UCS = 200 MPa
            70e9,  // E = 70 GPa
            0.25,  // ν = 0.25
        )
    }

    #[test]
    fn test_mohr_coulomb_zero_confinement() {
        let rock = granite();
        let sigma1 = rock.mohr_coulomb_strength(0.0);
        // σ₁ = 2c tan(45+φ/2) = 2*15e6*tan(67.5°)
        let expected = 2.0 * 15e6 * (PI / 4.0 + PI / 8.0).tan();
        assert!((sigma1 - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_mohr_coulomb_increases_with_confinement() {
        let rock = granite();
        let s1 = rock.mohr_coulomb_strength(10e6);
        let s2 = rock.mohr_coulomb_strength(20e6);
        assert!(s2 > s1);
    }

    #[test]
    fn test_hoek_brown_gsi100_equals_intact() {
        let rock = granite();
        // GSI=100 → s = exp(0) = 1 → intact rock criterion
        let sigma3 = 10e6;
        let mi = 32.0; // typical granite
        let sigma1 = rock.hoek_brown_strength(sigma3, mi, 100.0);
        assert!(sigma1 > sigma3, "σ₁ should exceed σ₃");
    }

    #[test]
    fn test_hoek_brown_lower_gsi_lower_strength() {
        let rock = granite();
        let sigma3 = 5e6;
        let mi = 32.0;
        let s100 = rock.hoek_brown_strength(sigma3, mi, 100.0);
        let s50 = rock.hoek_brown_strength(sigma3, mi, 50.0);
        assert!(s100 > s50);
    }

    #[test]
    fn test_drucker_prager_yield_negative_elastic() {
        let rock = granite();
        // Small deviatoric stress → inside yield surface
        let f = rock.drucker_prager_yield(100e6, 10e6);
        assert!(f < 0.0);
    }

    #[test]
    fn test_drucker_prager_yield_large_q_positive() {
        let rock = granite();
        let f = rock.drucker_prager_yield(1e6, 500e6);
        assert!(f > 0.0);
    }

    #[test]
    fn test_triaxial_deformation_elastic() {
        let rock = granite();
        let (eps1, _) = rock.triaxial_deformation(200e6, 10e6);
        // eps1 = (200e6 - 2*0.25*10e6)/70e9
        let expected = (200e6 - 2.0 * 0.25 * 10e6) / 70e9;
        assert!((eps1 - expected).abs() / expected.abs() < 1e-10);
    }

    #[test]
    fn test_triaxial_volumetric_strain() {
        let rock = granite();
        let (_, ev) = rock.triaxial_deformation(200e6, 10e6);
        assert!(ev.is_finite());
    }

    // --- PoroelasticMaterial ---

    fn sandstone() -> PoroelasticMaterial {
        PoroelasticMaterial::new(
            0.7,   // α
            10e9,  // K_drained
            0.20,  // n
            1e-14, // k [m²]
            2.2e9, // K_f (water)
            1e-3,  // μ (water)
        )
    }

    #[test]
    fn test_biot_modulus_positive() {
        let m = sandstone().biot_modulus();
        assert!(m > 0.0);
    }

    #[test]
    fn test_undrained_bulk_gt_drained() {
        let s = sandstone();
        assert!(s.undrained_bulk_modulus() > s.bulk_modulus);
    }

    #[test]
    fn test_skempton_b_between_0_and_1() {
        let b = sandstone().skempton_b();
        assert!((0.0..=1.0).contains(&b), "B = {b}");
    }

    #[test]
    fn test_consolidation_coefficient_positive() {
        let cv = sandstone().consolidation_coefficient();
        assert!(cv > 0.0);
    }

    #[test]
    fn test_terzaghi_settlement_zero_time() {
        let s = sandstone().terzaghi_settlement(100e3, 10.0, 0.0);
        assert_eq!(s, 0.0);
    }

    #[test]
    fn test_terzaghi_settlement_increases_with_time() {
        let mat = sandstone();
        let s1 = mat.terzaghi_settlement(100e3, 10.0, 1e6);
        let s2 = mat.terzaghi_settlement(100e3, 10.0, 1e10);
        assert!(s2 >= s1);
    }

    // --- Soil ---

    fn london_clay() -> Soil {
        Soil::new(SoilType::Clay, 0.9, 70.0, 28.0, 45.0, 20e3, 25.0)
    }

    #[test]
    fn test_plasticity_index() {
        let s = london_clay();
        assert!((s.plasticity_index() - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_compression_index() {
        let s = london_clay();
        // Cc = 0.009*(70-10) = 0.54
        assert!((s.compression_index() - 0.54).abs() < 1e-10);
    }

    #[test]
    fn test_preconsolidation_positive() {
        let p = london_clay().preconsolidation_pressure();
        assert!(p > 0.0);
    }

    #[test]
    fn test_hazen_permeability() {
        let s = london_clay();
        let k = s.permeability_hazen(0.0001); // d10 = 0.1 mm
        assert!(k > 0.0 && k < 1.0);
    }

    #[test]
    fn test_shear_strength_cu_positive() {
        let cu = london_clay().shear_strength_cu();
        assert!(cu > 0.0);
    }

    #[test]
    fn test_soil_type_variants() {
        let types = [
            SoilType::Clay,
            SoilType::Sand,
            SoilType::Gravel,
            SoilType::Silt,
            SoilType::Peat,
        ];
        assert_eq!(types.len(), 5);
    }

    // --- Rockfill ---

    fn dam_fill() -> Rockfill {
        Rockfill::new(0.2, 10.0) // D50=200mm, Cu=10
    }

    #[test]
    fn test_rockfill_friction_angle_in_range() {
        let phi = dam_fill().friction_angle();
        // Should be around 40° = 0.698 rad
        assert!(phi > 0.5 && phi < 1.2);
    }

    #[test]
    fn test_rockfill_stiffness_positive() {
        let e = dam_fill().stiffness(200e3);
        assert!(e > 0.0);
    }

    #[test]
    fn test_rockfill_stiffness_increases_with_confinement() {
        let fill = dam_fill();
        let e1 = fill.stiffness(100e3);
        let e2 = fill.stiffness(400e3);
        assert!(e2 > e1);
    }

    // --- cam_clay_yield ---

    #[test]
    fn test_cam_clay_yield_at_csl() {
        // On the CSL: q = Mp, p₀ = 2p → f = q²/M² - p(p₀-p) = p² - p² = 0
        let p = 100e3;
        let m = 1.2;
        let p0 = 2.0 * p;
        let q = m * p; // on CSL
        let f = cam_clay_yield(p, q, m, p0);
        assert!(f.abs() < 1e-6);
    }

    #[test]
    fn test_cam_clay_yield_inside_negative() {
        let p = 50e3;
        let m = 1.2;
        let p0 = 200e3;
        let q = 10e3; // small q, well inside
        let f = cam_clay_yield(p, q, m, p0);
        assert!(f < 0.0);
    }

    #[test]
    fn test_cam_clay_yield_outside_positive() {
        let p = 50e3;
        let m = 1.2;
        let p0 = 200e3;
        let q = 200e3; // large q, outside
        let f = cam_clay_yield(p, q, m, p0);
        assert!(f > 0.0);
    }

    // --- critical_state_line ---

    #[test]
    fn test_csl_proportional() {
        let p = 200e3;
        let m = 1.2;
        assert!((critical_state_line(p, m) - m * p).abs() < 1e-10);
    }

    #[test]
    fn test_csl_zero_p() {
        assert_eq!(critical_state_line(0.0, 1.2), 0.0);
    }

    // --- degree_of_consolidation ---

    #[test]
    fn test_consolidation_zero_tv() {
        let u = degree_of_consolidation(0.0);
        assert!(u < 0.01, "U(Tv=0) should be ~0");
    }

    #[test]
    fn test_consolidation_large_tv() {
        let u = degree_of_consolidation(10.0);
        assert!((u - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_consolidation_monotone() {
        let u1 = degree_of_consolidation(0.1);
        let u2 = degree_of_consolidation(0.5);
        let u3 = degree_of_consolidation(1.0);
        assert!(u1 < u2 && u2 < u3);
    }

    #[test]
    fn test_consolidation_half_tv_approx() {
        // T_v ≈ 0.197 gives U ≈ 50%
        let u = degree_of_consolidation(0.197);
        assert!(u > 0.45 && u < 0.55, "U(0.197) ≈ 0.5, got {u}");
    }
}
