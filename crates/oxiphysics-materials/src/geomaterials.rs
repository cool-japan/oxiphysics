// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geomechanical and geo-material constitutive models.
//!
//! Covers:
//! - Soil constitutive models: linear elastic, Mohr-Coulomb, Drucker-Prager
//! - Rock failure criteria: Hoek-Brown, Griffith
//! - Consolidation: Terzaghi 1D, Biot 3D poroelastic
//! - Liquefaction potential: SPT-based, cyclic stress ratio
//! - Permeability models: Darcy, Kozeny-Carman
//! - Compression index (Cc, Cs), preconsolidation pressure
//! - Effective stress principle
//! - Undrained shear strength
//! - Settlement calculation

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Stress tensor utilities (Voigt notation: [σxx, σyy, σzz, τxy, τyz, τzx])
// ---------------------------------------------------------------------------

/// Compute the mean stress (hydrostatic pressure) from a stress vector (Voigt notation).
///
/// p = (σxx + σyy + σzz) / 3
pub fn mean_stress(sigma: &[f64; 6]) -> f64 {
    (sigma[0] + sigma[1] + sigma[2]) / 3.0
}

/// Compute the deviatoric stress vector from a Voigt stress vector.
pub fn deviatoric_stress(sigma: &[f64; 6]) -> [f64; 6] {
    let p = mean_stress(sigma);
    [
        sigma[0] - p,
        sigma[1] - p,
        sigma[2] - p,
        sigma[3],
        sigma[4],
        sigma[5],
    ]
}

/// Compute the von Mises equivalent stress (J2 invariant based).
///
/// q = sqrt(3/2 * s:s)  where s is the deviatoric stress tensor.
pub fn von_mises_stress(sigma: &[f64; 6]) -> f64 {
    let s = deviatoric_stress(sigma);
    let j2 =
        0.5 * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    (3.0 * j2).sqrt()
}

/// Compute principal stresses from a 2D Mohr's circle (σ1 ≥ σ2).
///
/// For 2D plane stress: `sigma` = `[σxx, σyy, τxy]`.
/// Returns `(σ1, σ2, θp)` where θp is the principal angle in radians.
pub fn principal_stresses_2d(sxx: f64, syy: f64, txy: f64) -> (f64, f64, f64) {
    let avg = (sxx + syy) / 2.0;
    let r = (((sxx - syy) / 2.0).powi(2) + txy.powi(2)).sqrt();
    let sigma1 = avg + r;
    let sigma2 = avg - r;
    let theta_p = 0.5 * txy.atan2((sxx - syy) / 2.0);
    (sigma1, sigma2, theta_p)
}

// ---------------------------------------------------------------------------
// Linear elastic isotropic model
// ---------------------------------------------------------------------------

/// Linear isotropic elastic material for soils.
#[derive(Debug, Clone)]
pub struct LinearElasticSoil {
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio ν \[-\].
    pub poisson_ratio: f64,
}

impl LinearElasticSoil {
    /// Create a new linear elastic soil.
    pub fn new(young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
        }
    }

    /// Bulk modulus K = E / (3(1-2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Compute elastic strain from stress increment (Voigt notation).
    ///
    /// ε = C⁻¹ σ  (isotropic compliance).
    pub fn strain_from_stress(&self, sigma: &[f64; 6]) -> [f64; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let g = self.shear_modulus();
        [
            (sigma[0] - nu * (sigma[1] + sigma[2])) / e,
            (sigma[1] - nu * (sigma[0] + sigma[2])) / e,
            (sigma[2] - nu * (sigma[0] + sigma[1])) / e,
            sigma[3] / (2.0 * g),
            sigma[4] / (2.0 * g),
            sigma[5] / (2.0 * g),
        ]
    }

    /// Compute stress from strain increment (Voigt notation).
    pub fn stress_from_strain(&self, eps: &[f64; 6]) -> [f64; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let g = self.shear_modulus();
        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let _ev = eps[0] + eps[1] + eps[2];
        [
            factor * ((1.0 - nu) * eps[0] + nu * (eps[1] + eps[2])),
            factor * ((1.0 - nu) * eps[1] + nu * (eps[0] + eps[2])),
            factor * ((1.0 - nu) * eps[2] + nu * (eps[0] + eps[1])),
            2.0 * g * eps[3],
            2.0 * g * eps[4],
            2.0 * g * eps[5],
        ]
    }

    /// Vertical settlement under applied vertical stress `delta_sigma_v` over depth `h`.
    ///
    /// δ = Δσ_v · h / E_oed   where E_oed = E(1-ν)/((1+ν)(1-2ν)).
    pub fn settlement_1d(&self, delta_sigma_v: f64, h: f64) -> f64 {
        let nu = self.poisson_ratio;
        let e_oed = self.young_modulus * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        delta_sigma_v * h / e_oed
    }
}

// ---------------------------------------------------------------------------
// Mohr-Coulomb failure criterion
// ---------------------------------------------------------------------------

/// Mohr-Coulomb failure criterion for soils and rocks.
///
/// Failure envelope: τ = c + σ_n tan(φ)
#[derive(Debug, Clone)]
pub struct MohrCoulomb {
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
    /// Friction angle φ \[radians\].
    pub friction_angle: f64,
}

impl MohrCoulomb {
    /// Create a Mohr-Coulomb model with cohesion `c` and friction angle `phi` (radians).
    pub fn new(cohesion: f64, friction_angle: f64) -> Self {
        Self {
            cohesion,
            friction_angle,
        }
    }

    /// Shear strength at a given normal stress.
    pub fn shear_strength(&self, sigma_n: f64) -> f64 {
        self.cohesion + sigma_n * self.friction_angle.tan()
    }

    /// Check if a stress state (normal stress σ_n, shear stress τ) has failed.
    ///
    /// Returns `true` if `τ ≥ c + σ_n tan(φ)`.
    pub fn is_failed(&self, sigma_n: f64, tau: f64) -> bool {
        tau >= self.shear_strength(sigma_n)
    }

    /// Evaluate the yield function F for principal stresses σ1 ≥ σ3.
    ///
    /// F = σ1 - σ3 * N_φ - 2c√N_φ  where N_φ = (1 + sin φ)/(1 - sin φ).
    pub fn yield_function(&self, sigma1: f64, sigma3: f64) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let n_phi = (1.0 + sin_phi) / (1.0 - sin_phi);
        sigma1 - sigma3 * n_phi - 2.0 * self.cohesion * n_phi.sqrt()
    }

    /// Compute the uniaxial compressive strength (UCS).
    ///
    /// UCS = 2c cos(φ) / (1 - sin(φ)).
    pub fn ucs(&self) -> f64 {
        let phi = self.friction_angle;
        2.0 * self.cohesion * phi.cos() / (1.0 - phi.sin())
    }

    /// Compute the uniaxial tensile strength.
    ///
    /// σ_t = 2c cos(φ) / (1 + sin(φ)).
    pub fn tensile_strength(&self) -> f64 {
        let phi = self.friction_angle;
        2.0 * self.cohesion * phi.cos() / (1.0 + phi.sin())
    }

    /// Compute Drucker-Prager equivalent parameters (inner cone matching Mohr-Coulomb).
    ///
    /// Returns `(alpha, k)` for the DP criterion F = α·I1 + √J2 - k.
    pub fn to_drucker_prager_inner(&self) -> (f64, f64) {
        let sin_phi = self.friction_angle.sin();
        let alpha = 2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi));
        let k =
            6.0 * self.cohesion * self.friction_angle.cos() / (3.0_f64.sqrt() * (3.0 - sin_phi));
        (alpha, k)
    }

    /// Compute the plane-strain lateral earth pressure coefficient K0 = ν/(1-ν) for elastic.
    pub fn at_rest_pressure_coefficient(poisson_ratio: f64) -> f64 {
        poisson_ratio / (1.0 - poisson_ratio)
    }

    /// Rankine's active earth pressure coefficient Ka.
    ///
    /// Ka = tan²(45° - φ/2).
    pub fn active_pressure_coefficient(&self) -> f64 {
        let angle = PI / 4.0 - self.friction_angle / 2.0;
        angle.tan().powi(2)
    }

    /// Rankine's passive earth pressure coefficient Kp.
    ///
    /// Kp = tan²(45° + φ/2).
    pub fn passive_pressure_coefficient(&self) -> f64 {
        let angle = PI / 4.0 + self.friction_angle / 2.0;
        angle.tan().powi(2)
    }
}

// ---------------------------------------------------------------------------
// Drucker-Prager model
// ---------------------------------------------------------------------------

/// Drucker-Prager elasto-plastic constitutive model.
///
/// Yield surface: F = α·I1 + √J2 - k ≤ 0
/// where I1 = σ_kk (first stress invariant) and J2 is the second deviatoric invariant.
#[derive(Debug, Clone)]
pub struct DruckerPragerModel {
    /// Cohesion-like parameter k \[Pa\].
    pub k: f64,
    /// Pressure sensitivity parameter α \[-\].
    pub alpha: f64,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio \[-\].
    pub poisson_ratio: f64,
}

impl DruckerPragerModel {
    /// Create a Drucker-Prager model.
    pub fn new(k: f64, alpha: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            k,
            alpha,
            young_modulus,
            poisson_ratio,
        }
    }

    /// Create a Drucker-Prager model matching Mohr-Coulomb (inner cone).
    pub fn from_mohr_coulomb(mc: &MohrCoulomb, young_modulus: f64, poisson_ratio: f64) -> Self {
        let (alpha, k) = mc.to_drucker_prager_inner();
        Self::new(k, alpha, young_modulus, poisson_ratio)
    }

    /// Evaluate the yield function F for a Voigt stress state.
    ///
    /// F = α·I1 + √J2 - k
    pub fn yield_function(&self, sigma: &[f64; 6]) -> f64 {
        let i1 = sigma[0] + sigma[1] + sigma[2];
        let s = deviatoric_stress(sigma);
        let j2 = 0.5 * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2])
            + s[3] * s[3]
            + s[4] * s[4]
            + s[5] * s[5];
        self.alpha * i1 + j2.sqrt() - self.k
    }

    /// Check if the stress state is beyond yield.
    pub fn is_yielded(&self, sigma: &[f64; 6]) -> bool {
        self.yield_function(sigma) > 0.0
    }

    /// Compute the plastic flow direction (gradient of yield surface w.r.t. stress).
    ///
    /// ∂F/∂σ = α·δ + (1/(2√J2)) s
    pub fn flow_direction(&self, sigma: &[f64; 6]) -> [f64; 6] {
        let s = deviatoric_stress(sigma);
        let j2 = 0.5 * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2])
            + s[3] * s[3]
            + s[4] * s[4]
            + s[5] * s[5];
        let sqrt_j2 = j2.sqrt().max(1e-20);
        let factor = 1.0 / (2.0 * sqrt_j2);
        [
            self.alpha + factor * s[0],
            self.alpha + factor * s[1],
            self.alpha + factor * s[2],
            factor * s[3],
            factor * s[4],
            factor * s[5],
        ]
    }

    /// Return-mapping projection: project trial stress back onto yield surface.
    ///
    /// Uses the radial return algorithm for the Drucker-Prager cone.
    pub fn return_mapping(&self, sigma_trial: &[f64; 6]) -> [f64; 6] {
        if !self.is_yielded(sigma_trial) {
            return *sigma_trial;
        }
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let k_bulk = e / (3.0 * (1.0 - 2.0 * nu));
        let g = e / (2.0 * (1.0 + nu));

        let i1_trial = sigma_trial[0] + sigma_trial[1] + sigma_trial[2];
        let s_trial = deviatoric_stress(sigma_trial);
        let j2_trial = 0.5
            * (s_trial[0] * s_trial[0] + s_trial[1] * s_trial[1] + s_trial[2] * s_trial[2])
            + s_trial[3] * s_trial[3]
            + s_trial[4] * s_trial[4]
            + s_trial[5] * s_trial[5];
        let sqrt_j2 = j2_trial.sqrt().max(1e-20);

        // Plastic multiplier: Δγ = F_trial / (G + 9α²K)
        let f_trial = self.alpha * i1_trial + sqrt_j2 - self.k;
        let d_gamma = f_trial / (g + 9.0 * self.alpha * self.alpha * k_bulk);

        // Check for apex return (deviatoric stress would go to zero or flip sign)
        let scale = 1.0 - g * d_gamma / sqrt_j2;
        if scale <= 0.0 {
            // Project to apex: J2 = 0, p = k / (3α)
            let p_apex = self.k / (3.0 * self.alpha.max(1e-30));
            return [p_apex, p_apex, p_apex, 0.0, 0.0, 0.0];
        }

        // Corrected stress
        let p_new = i1_trial / 3.0 - 3.0 * self.alpha * k_bulk * d_gamma;
        [
            p_new + scale * s_trial[0],
            p_new + scale * s_trial[1],
            p_new + scale * s_trial[2],
            scale * s_trial[3],
            scale * s_trial[4],
            scale * s_trial[5],
        ]
    }
}

// ---------------------------------------------------------------------------
// Hoek-Brown rock failure criterion
// ---------------------------------------------------------------------------

/// Hoek-Brown failure criterion for jointed rock masses.
///
/// σ1 = σ3 + σci · (mb · σ3/σci + s)^a
#[derive(Debug, Clone)]
pub struct HoekBrown {
    /// Uniaxial compressive strength of intact rock σci \[Pa\].
    pub sigma_ci: f64,
    /// Hoek-Brown constant for the rock mass mb \[-\].
    pub mb: f64,
    /// Rock mass parameter s \[-\] (0 for completely jointed, 1 for intact).
    pub s: f64,
    /// Rock mass parameter a \[-\] (typically 0.5 for intact).
    pub a: f64,
}

impl HoekBrown {
    /// Create a Hoek-Brown model.
    pub fn new(sigma_ci: f64, mb: f64, s: f64, a: f64) -> Self {
        Self { sigma_ci, mb, s, a }
    }

    /// Create for intact rock (s=1, a=0.5, mb from mi).
    pub fn intact_rock(sigma_ci: f64, mi: f64) -> Self {
        Self::new(sigma_ci, mi, 1.0, 0.5)
    }

    /// Compute σ1 at failure given confining stress σ3.
    pub fn sigma1_at_failure(&self, sigma3: f64) -> f64 {
        let term = self.mb * sigma3 / self.sigma_ci + self.s;
        if term < 0.0 {
            sigma3
        } else {
            sigma3 + self.sigma_ci * term.powf(self.a)
        }
    }

    /// Uniaxial compressive strength of the rock mass.
    ///
    /// UCS_m = σci · s^a.
    pub fn ucs_mass(&self) -> f64 {
        self.sigma_ci * self.s.powf(self.a)
    }

    /// Uniaxial tensile strength.
    ///
    /// σt = σci/2 · (mb - sqrt(mb² + 4s)).
    pub fn tensile_strength(&self) -> f64 {
        let disc = self.mb * self.mb + 4.0 * self.s;
        self.sigma_ci / 2.0 * (self.mb - disc.sqrt())
    }

    /// Geological Strength Index (GSI) from mb, s (approximate inverse).
    ///
    /// Uses the relation: mb = mi · exp((GSI-100)/28), s = exp((GSI-100)/9).
    /// Solves via s: GSI ≈ 9 · ln(s) + 100.
    pub fn gsi_from_s(&self) -> f64 {
        9.0 * self.s.ln() + 100.0
    }

    /// Mohr-Coulomb equivalent friction angle at confining stress σ3 \[radians\].
    ///
    /// Linearizes the Hoek-Brown envelope at the given confining stress.
    pub fn equivalent_friction_angle(&self, sigma3: f64) -> f64 {
        let h = self.sigma_ci
            * (self.mb * sigma3 / self.sigma_ci + self.s).powf(self.a - 1.0)
            * self.mb
            * self.a;
        let sin_phi = (h - 1.0) / (h + 1.0);
        sin_phi.asin().max(0.0)
    }

    /// Equivalent cohesion c at a given confining stress σ3 \[Pa\].
    pub fn equivalent_cohesion(&self, sigma3: f64) -> f64 {
        let phi = self.equivalent_friction_angle(sigma3);
        let sigma1 = self.sigma1_at_failure(sigma3);
        let tau = (sigma1 - sigma3) / 2.0;
        let sigma_n = (sigma1 + sigma3) / 2.0;
        (tau - sigma_n * phi.tan()).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// Griffith failure criterion for brittle materials
// ---------------------------------------------------------------------------

/// Griffith crack theory failure criterion.
///
/// For tension: σ_t = T0 (tensile strength).
/// For compression-tension: (σ1 - σ3)² + 8T0(σ1 + σ3) = 0.
#[derive(Debug, Clone)]
pub struct GriffithCriterion {
    /// Tensile strength T0 \[Pa\].
    pub tensile_strength: f64,
}

impl GriffithCriterion {
    /// Create a Griffith criterion model.
    pub fn new(tensile_strength: f64) -> Self {
        Self { tensile_strength }
    }

    /// Evaluate the Griffith yield function for principal stresses (σ1 ≥ σ3).
    ///
    /// Returns positive value if failed.
    pub fn yield_function(&self, sigma1: f64, sigma3: f64) -> f64 {
        let t0 = self.tensile_strength;
        if sigma1 + 3.0 * sigma3 < 0.0 {
            // Tensile regime: simple tensile cut-off (failure when σ3 ≤ -T0).
            -(sigma3 + t0)
        } else {
            // Compressive-dominated regime: Griffith parabolic criterion.
            (sigma1 - sigma3).powi(2) - 8.0 * t0 * (sigma1 + sigma3)
        }
    }

    /// Check if the stress state is at failure.
    pub fn is_failed(&self, sigma1: f64, sigma3: f64) -> bool {
        self.yield_function(sigma1, sigma3) >= 0.0
    }

    /// Uniaxial compressive strength from Griffith: UCS = 8 T0.
    pub fn ucs(&self) -> f64 {
        8.0 * self.tensile_strength
    }

    /// Critical stress intensity factor approximation K_Ic ≈ T0 √(π a) for penny-shaped crack.
    pub fn stress_intensity_factor(&self, crack_half_length: f64) -> f64 {
        self.tensile_strength * (PI * crack_half_length).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Effective stress principle (Terzaghi)
// ---------------------------------------------------------------------------

/// Compute effective stress from total stress and pore water pressure.
///
/// σ' = σ - u  (Terzaghi effective stress principle)
pub fn effective_stress(total_stress: f64, pore_pressure: f64) -> f64 {
    total_stress - pore_pressure
}

/// Compute pore water pressure at depth z below water table.
///
/// u = γ_w · z_w  where z_w is depth below phreatic surface.
pub fn hydrostatic_pore_pressure(depth_below_water_table: f64, unit_weight_water: f64) -> f64 {
    unit_weight_water * depth_below_water_table
}

/// Compute vertical total stress from overburden.
///
/// σ_v = Σ γ_i · h_i
pub fn overburden_stress(layers: &[(f64, f64)]) -> f64 {
    // Each layer: (unit_weight, thickness)
    layers.iter().map(|(gamma, h)| gamma * h).sum()
}

// ---------------------------------------------------------------------------
// Terzaghi 1D consolidation
// ---------------------------------------------------------------------------

/// Terzaghi 1D consolidation model for fully saturated soils.
#[derive(Debug, Clone)]
pub struct TerzaghiConsolidation {
    /// Coefficient of consolidation cv \[m²/s\].
    pub cv: f64,
    /// Drainage path length H \[m\] (half-thickness for double drainage).
    pub drainage_path: f64,
    /// Initial excess pore pressure u0 \[Pa\].
    pub initial_excess_pore_pressure: f64,
}

impl TerzaghiConsolidation {
    /// Create a Terzaghi consolidation model.
    pub fn new(cv: f64, drainage_path: f64, initial_excess_pore_pressure: f64) -> Self {
        Self {
            cv,
            drainage_path,
            initial_excess_pore_pressure,
        }
    }

    /// Time factor T_v = cv · t / H².
    pub fn time_factor(&self, t: f64) -> f64 {
        self.cv * t / (self.drainage_path * self.drainage_path)
    }

    /// Average degree of consolidation U(T_v) (Terzaghi series solution).
    ///
    /// U ≈ 1 - Σ (2/M²) exp(-M² T_v)  for M = π(2m+1)/2.
    pub fn average_degree_of_consolidation(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let tv = self.time_factor(t);
        let mut u = 0.0f64;
        for m in 0..=50 {
            let big_m = PI * (2 * m + 1) as f64 / 2.0;
            u += (2.0 / (big_m * big_m)) * (-big_m * big_m * tv).exp();
        }
        1.0 - u
    }

    /// Approximate average degree of consolidation using Terzaghi's approximation.
    ///
    /// For U < 0.6: U ≈ sqrt(4*Tv/π)
    /// For U ≥ 0.6: U ≈ 1 - 10^(-(Tv+0.085)/0.933) (approximate)
    pub fn average_doc_approx(&self, t: f64) -> f64 {
        let tv = self.time_factor(t);
        if tv <= 0.217 {
            (4.0 * tv / PI).sqrt()
        } else {
            1.0 - (-1.781 * tv + 0.933_f64.ln()).exp().min(1.0)
        }
    }

    /// Excess pore pressure distribution at depth z and time t.
    ///
    /// u(z,t) = u0 · Σ (4/((2m+1)π)) sin(M·z/H) exp(-M²·T_v)
    pub fn excess_pore_pressure(&self, z: f64, t: f64) -> f64 {
        let tv = self.time_factor(t);
        let h = self.drainage_path;
        let mut u = 0.0f64;
        for m in 0..=50 {
            let big_m = PI * (2 * m + 1) as f64 / 2.0;
            u += (4.0 / ((2 * m + 1) as f64 * PI))
                * (big_m * z / h).sin()
                * (-big_m * big_m * tv).exp();
        }
        self.initial_excess_pore_pressure * u
    }

    /// Settlement at time t \[m\] for a layer of thickness `H_total`.
    ///
    /// s(t) = U(t) · s_final
    pub fn settlement_at_time(&self, t: f64, final_settlement: f64) -> f64 {
        self.average_degree_of_consolidation(t) * final_settlement
    }
}

// ---------------------------------------------------------------------------
// Biot 3D poroelastic model
// ---------------------------------------------------------------------------

/// Biot 3D poroelastic parameters.
#[derive(Debug, Clone)]
pub struct BiotPoroelastic {
    /// Drained bulk modulus K_d \[Pa\].
    pub bulk_modulus_drained: f64,
    /// Shear modulus G \[Pa\] (same drained/undrained for incompressible fluid).
    pub shear_modulus: f64,
    /// Biot coefficient α \[-\] (= 1 - K_d/K_s).
    pub biot_coefficient: f64,
    /// Biot modulus M \[Pa\] (= 1/(n/K_f + (α-n)/K_s)).
    pub biot_modulus: f64,
    /// Permeability k \[m²\].
    pub permeability: f64,
    /// Fluid viscosity μ \[Pa·s\].
    pub fluid_viscosity: f64,
}

impl BiotPoroelastic {
    /// Create a Biot poroelastic model.
    pub fn new(
        bulk_modulus_drained: f64,
        shear_modulus: f64,
        biot_coefficient: f64,
        biot_modulus: f64,
        permeability: f64,
        fluid_viscosity: f64,
    ) -> Self {
        Self {
            bulk_modulus_drained,
            shear_modulus,
            biot_coefficient,
            biot_modulus,
            permeability,
            fluid_viscosity,
        }
    }

    /// Undrained bulk modulus K_u = K_d + α² M.
    pub fn bulk_modulus_undrained(&self) -> f64 {
        self.bulk_modulus_drained
            + self.biot_coefficient * self.biot_coefficient * self.biot_modulus
    }

    /// Skempton's pore pressure coefficient B.
    ///
    /// B = α M / K_u.
    pub fn skempton_b(&self) -> f64 {
        let ku = self.bulk_modulus_undrained();
        if ku < 1e-20 {
            0.0
        } else {
            self.biot_coefficient * self.biot_modulus / ku
        }
    }

    /// Hydraulic diffusivity c = k / (μ · S), where S = 1/M is specific storage.
    pub fn hydraulic_diffusivity(&self) -> f64 {
        let s = 1.0 / self.biot_modulus;
        self.permeability / (self.fluid_viscosity * s)
    }

    /// Pore pressure induced by isotropic total stress increment Δσ (undrained).
    ///
    /// Δu = B · Δσ_mean.
    pub fn undrained_pore_pressure_response(&self, delta_sigma_mean: f64) -> f64 {
        self.skempton_b() * delta_sigma_mean
    }

    /// Compute the effective stress from total stress Voigt vector and pore pressure.
    pub fn effective_stress_voigt(&self, sigma: &[f64; 6], pore_pressure: f64) -> [f64; 6] {
        let alpha = self.biot_coefficient;
        [
            sigma[0] - alpha * pore_pressure,
            sigma[1] - alpha * pore_pressure,
            sigma[2] - alpha * pore_pressure,
            sigma[3],
            sigma[4],
            sigma[5],
        ]
    }
}

// ---------------------------------------------------------------------------
// Liquefaction potential
// ---------------------------------------------------------------------------

/// SPT-based liquefaction potential assessment (Youd & Idriss, 2001 simplified procedure).
#[derive(Debug, Clone)]
pub struct LiquefactionAssessment {
    /// Peak ground acceleration amax \[g\].
    pub amax: f64,
    /// Earthquake magnitude Mw \[-\].
    pub magnitude: f64,
    /// Total vertical stress σ_v \[kPa\].
    pub total_vertical_stress: f64,
    /// Effective vertical stress σ'_v \[kPa\].
    pub effective_vertical_stress: f64,
    /// SPT blow count N60 (corrected to 60% energy) \[-\].
    pub n60: f64,
    /// Fines content FC \[%\].
    pub fines_content: f64,
}

impl LiquefactionAssessment {
    /// Create a liquefaction assessment.
    pub fn new(
        amax: f64,
        magnitude: f64,
        total_vertical_stress: f64,
        effective_vertical_stress: f64,
        n60: f64,
        fines_content: f64,
    ) -> Self {
        Self {
            amax,
            magnitude,
            total_vertical_stress,
            effective_vertical_stress,
            n60,
            fines_content,
        }
    }

    /// Stress reduction coefficient rd (simplified Seed & Idriss).
    ///
    /// rd ≈ 1.0 - 0.00765z for z ≤ 9.15 m, 1.174 - 0.0267z for 9.15 < z ≤ 23 m.
    pub fn stress_reduction_coefficient(depth_m: f64) -> f64 {
        if depth_m <= 9.15 {
            1.0 - 0.00765 * depth_m
        } else if depth_m <= 23.0 {
            1.174 - 0.0267 * depth_m
        } else {
            0.564
        }
    }

    /// Cyclic Stress Ratio CSR_7.5 induced by the earthquake.
    ///
    /// CSR = 0.65 · (σ_v/σ'_v) · amax · rd
    pub fn cyclic_stress_ratio(&self, depth_m: f64) -> f64 {
        let rd = Self::stress_reduction_coefficient(depth_m);
        0.65 * (self.total_vertical_stress / self.effective_vertical_stress) * self.amax * rd
    }

    /// Magnitude scaling factor MSF (Idriss, 1999).
    ///
    /// MSF = 10^2.24 / Mw^2.56.
    pub fn magnitude_scaling_factor(&self) -> f64 {
        10.0_f64.powf(2.24) / self.magnitude.powf(2.56)
    }

    /// Normalize SPT blow count to σ'_v = 100 kPa: (N1)60.
    ///
    /// CN = min(1.7, (Pa/σ'_v)^0.5)  where Pa = 100 kPa = 100 000 Pa.
    pub fn normalized_n1_60(&self) -> f64 {
        let cn = (100_000.0 / self.effective_vertical_stress).sqrt().min(1.7);
        cn * self.n60
    }

    /// Fines-corrected (N1)60cs.
    pub fn corrected_n1_60cs(&self) -> f64 {
        let n1_60 = self.normalized_n1_60();
        let delta_n = if self.fines_content < 5.0 {
            0.0
        } else if self.fines_content < 35.0 {
            -17.0 + 0.76 * (1600.0 / (self.fines_content + 9.0) + 17.0_f64.powi(2)).sqrt()
        } else {
            5.0
        };
        n1_60 + delta_n
    }

    /// Cyclic Resistance Ratio CRR_7.5 from SPT (Robertson & Wride, 1998).
    ///
    /// For (N1)60cs < 30:
    /// CRR = 1/(34 - N1_60cs) + N1_60cs/135 + 50/(10·N1_60cs+45)² - 1/200
    pub fn cyclic_resistance_ratio_7_5(&self) -> f64 {
        let n = self.corrected_n1_60cs();
        if n >= 30.0 {
            return f64::INFINITY; // Non-liquefiable
        }
        1.0 / (34.0 - n) + n / 135.0 + 50.0 / (10.0 * n + 45.0).powi(2) - 1.0 / 200.0
    }

    /// Factor of safety against liquefaction FL.
    ///
    /// FL = CRR_7.5 · MSF / CSR
    pub fn factor_of_safety(&self, depth_m: f64) -> f64 {
        let crr = self.cyclic_resistance_ratio_7_5();
        if crr.is_infinite() {
            return f64::INFINITY;
        }
        let csr = self.cyclic_stress_ratio(depth_m);
        let msf = self.magnitude_scaling_factor();
        crr * msf / csr
    }

    /// Check if liquefaction is predicted (FL < 1.0).
    pub fn is_liquefiable(&self, depth_m: f64) -> bool {
        self.factor_of_safety(depth_m) < 1.0
    }
}

// ---------------------------------------------------------------------------
// Darcy's law and permeability
// ---------------------------------------------------------------------------

/// Darcy flow through a porous medium.
///
/// q = -k/μ · ∇p  (Darcy velocity vector \[m/s\]).
pub fn darcy_flow(permeability: f64, viscosity: f64, pressure_gradient: f64) -> f64 {
    -permeability / viscosity * pressure_gradient
}

/// Kozeny-Carman permeability model.
///
/// k = (e³/(1+e)) · Cs · d²  where Cs is the Kozeny-Carman constant (≈ 1/180 for spheres),
/// d is particle diameter, and e is void ratio.
pub fn kozeny_carman_permeability(void_ratio: f64, particle_diameter: f64) -> f64 {
    let cs = 1.0 / 180.0;
    let e = void_ratio;
    cs * particle_diameter * particle_diameter * e * e * e / (1.0 + e)
}

/// Hydraulic conductivity K \[m/s\] from intrinsic permeability.
///
/// K = k · γ_w / μ  where γ_w = ρ_w · g.
pub fn hydraulic_conductivity(
    permeability: f64,
    fluid_density: f64,
    gravity: f64,
    viscosity: f64,
) -> f64 {
    permeability * fluid_density * gravity / viscosity
}

// ---------------------------------------------------------------------------
// Compression index and consolidation settlement
// ---------------------------------------------------------------------------

/// Soil compressibility parameters.
#[derive(Debug, Clone)]
pub struct SoilCompressibility {
    /// Compression index Cc \[-\] (slope of e-log(σ') in virgin compression).
    pub compression_index: f64,
    /// Recompression index Cs \[-\] (slope in recompression).
    pub recompression_index: f64,
    /// Initial void ratio e0 \[-\].
    pub initial_void_ratio: f64,
    /// Preconsolidation pressure σ'_p \[Pa\].
    pub preconsolidation_pressure: f64,
}

impl SoilCompressibility {
    /// Create soil compressibility parameters.
    pub fn new(
        compression_index: f64,
        recompression_index: f64,
        initial_void_ratio: f64,
        preconsolidation_pressure: f64,
    ) -> Self {
        Self {
            compression_index,
            recompression_index,
            initial_void_ratio,
            preconsolidation_pressure,
        }
    }

    /// Estimate Cc from Skempton's correlation: Cc ≈ 0.009(LL - 10).
    pub fn cc_from_liquid_limit(liquid_limit_percent: f64) -> f64 {
        0.009 * (liquid_limit_percent - 10.0).max(0.0)
    }

    /// Estimate Cc from Terzaghi's correlation: Cc ≈ 0.54(e0 - 0.35).
    pub fn cc_from_void_ratio(initial_void_ratio: f64) -> f64 {
        (0.54 * (initial_void_ratio - 0.35)).max(0.0)
    }

    /// Compute consolidation settlement for a layer of thickness H \[m\].
    ///
    /// Normally consolidated: δ = Cc·H/(1+e0) · log10(σ'f/σ'0)
    /// Overconsolidated (σ'f < σ'p): δ = Cs·H/(1+e0) · log10(σ'f/σ'0)
    /// Overconsolidated (σ'f > σ'p): sum of both branches.
    pub fn consolidation_settlement(
        &self,
        layer_thickness: f64,
        initial_effective_stress: f64,
        final_effective_stress: f64,
    ) -> f64 {
        let h = layer_thickness;
        let e0 = self.initial_void_ratio;
        let sigma0 = initial_effective_stress;
        let sigmap = self.preconsolidation_pressure;
        let sigmaf = final_effective_stress;

        if sigmaf <= sigmap {
            // Overconsolidated (recompression only)
            self.recompression_index * h / (1.0 + e0) * (sigmaf / sigma0).log10()
        } else if sigma0 >= sigmap {
            // Normally consolidated
            self.compression_index * h / (1.0 + e0) * (sigmaf / sigma0).log10()
        } else {
            // Crosses preconsolidation: two branches
            let s1 = self.recompression_index * h / (1.0 + e0) * (sigmap / sigma0).log10();
            let s2 = self.compression_index * h / (1.0 + e0) * (sigmaf / sigmap).log10();
            s1 + s2
        }
    }

    /// Overconsolidation ratio OCR = σ'_p / σ'_0.
    pub fn ocr(&self, initial_effective_stress: f64) -> f64 {
        self.preconsolidation_pressure / initial_effective_stress
    }

    /// Coefficient of volume change mv \[1/Pa\].
    ///
    /// mv ≈ Cc / (2.303 · (1 + e0) · σ'_avg).
    pub fn mv(&self, effective_stress: f64) -> f64 {
        self.compression_index / (2.303 * (1.0 + self.initial_void_ratio) * effective_stress)
    }
}

// ---------------------------------------------------------------------------
// Undrained shear strength
// ---------------------------------------------------------------------------

/// Undrained shear strength correlations and models.
#[derive(Debug, Clone)]
pub struct UndrainedShearStrength;

impl UndrainedShearStrength {
    /// Skempton's correlation: Su ≈ 0.11 + 0.0037·PI for NC clay.
    ///
    /// `pi` is the plasticity index in percent.
    pub fn skempton_nc(effective_overburden: f64, plasticity_index: f64) -> f64 {
        effective_overburden * (0.11 + 0.0037 * plasticity_index)
    }

    /// Su from vane shear test (corrected for anisotropy).
    ///
    /// Su_corrected = μ · Su_vane  where μ = 1.05 - 0.05·(PI/10) (Bjerrum).
    pub fn vane_corrected(su_vane: f64, plasticity_index: f64) -> f64 {
        let mu = (1.05 - 0.005 * plasticity_index).clamp(0.5, 1.0);
        mu * su_vane
    }

    /// Normalize undrained shear strength ratio Su/σ'_v for OCR soils.
    ///
    /// SHANSEP: (Su/σ'v) = S · OCR^m  where S ≈ 0.22 and m ≈ 0.80.
    pub fn shansep(effective_vertical_stress: f64, ocr: f64) -> f64 {
        let s = 0.22;
        let m = 0.80;
        effective_vertical_stress * s * ocr.powf(m)
    }

    /// Compute Su from triaxial test (undrained): Su = (σ1 - σ3)/2.
    pub fn from_triaxial(sigma1: f64, sigma3: f64) -> f64 {
        (sigma1 - sigma3) / 2.0
    }
}

// ---------------------------------------------------------------------------
// Settlement calculation
// ---------------------------------------------------------------------------

/// Immediate (elastic) settlement under a uniformly loaded flexible circular footing.
///
/// δ_i = q · B · (1-ν²) / E · I_s  where I_s = π/4 for rigid circular footing.
pub fn immediate_settlement_circular(
    applied_pressure: f64,
    diameter: f64,
    young_modulus: f64,
    poisson_ratio: f64,
) -> f64 {
    let b = diameter;
    let i_s = PI / 4.0;
    applied_pressure * b * (1.0 - poisson_ratio * poisson_ratio) / young_modulus * i_s
}

/// Elastic settlement of a rectangular footing (Steinbrenner's method).
///
/// δ = q · B · (1-ν²) / E · I_w  where I_w is the influence factor.
pub fn immediate_settlement_rectangular(
    applied_pressure: f64,
    width: f64,
    length: f64,
    young_modulus: f64,
    poisson_ratio: f64,
) -> f64 {
    // Simplified influence factor for L/B ratio
    let m = length / width;
    let i_w = 0.5 * (m * m + 1.0).sqrt() / m * (m + (m * m + 1.0).sqrt()).ln() / PI;
    applied_pressure * width * (1.0 - poisson_ratio * poisson_ratio) / young_modulus * i_w * 4.0
}

/// Secondary (creep) compression settlement.
///
/// δ_s = Cα · H / (1 + ep) · log10(t2/t1)
///
/// where Cα is the secondary compression index.
pub fn secondary_compression_settlement(
    secondary_compression_index: f64,
    layer_thickness: f64,
    void_ratio_at_end_of_primary: f64,
    t1: f64,
    t2: f64,
) -> f64 {
    if t2 <= t1 {
        return 0.0;
    }
    secondary_compression_index * layer_thickness / (1.0 + void_ratio_at_end_of_primary)
        * (t2 / t1).log10()
}

// ---------------------------------------------------------------------------
// Bearing capacity
// ---------------------------------------------------------------------------

/// Terzaghi bearing capacity for strip footing.
///
/// qu = c·Nc + q·Nq + 0.5·γ·B·Nγ
pub fn terzaghi_bearing_capacity_strip(
    cohesion: f64,
    surcharge: f64,
    unit_weight: f64,
    width: f64,
    friction_angle: f64,
) -> f64 {
    let phi = friction_angle;
    let nq = (PI * phi.tan()).exp() * ((PI / 4.0 + phi / 2.0).tan()).powi(2);
    let nc = if phi < 1e-9 {
        PI + 2.0
    } else {
        (nq - 1.0) / phi.tan()
    };
    let ng = 2.0 * (nq + 1.0) * phi.tan();
    cohesion * nc + surcharge * nq + 0.5 * unit_weight * width * ng
}

/// Meyerhof bearing capacity factors.
///
/// Returns `(Nc, Nq, Nγ)` for friction angle φ \[radians\].
pub fn meyerhof_bearing_factors(friction_angle: f64) -> (f64, f64, f64) {
    let phi = friction_angle;
    let nq = (PI * phi.tan()).exp() * ((PI / 4.0 + phi / 2.0).tan()).powi(2);
    let nc = if phi < 1e-9 {
        PI + 2.0
    } else {
        (nq - 1.0) / phi.tan()
    };
    let ng = (nq - 1.0) * (1.4 * phi).tan();
    (nc, nq, ng)
}

// ---------------------------------------------------------------------------
// Stress path utilities
// ---------------------------------------------------------------------------

/// Stress path in p-q space for triaxial loading.
///
/// p = (σ1 + 2σ3)/3  (mean total stress),  q = σ1 - σ3 (deviatoric).
pub fn triaxial_stress_path(sigma1: f64, sigma3: f64) -> (f64, f64) {
    let p = (sigma1 + 2.0 * sigma3) / 3.0;
    let q = sigma1 - sigma3;
    (p, q)
}

/// Compute the K0 (at-rest) stress for a given overburden.
pub fn k0_stress(unit_weight: f64, depth: f64, k0: f64) -> (f64, f64) {
    let sigma_v = unit_weight * depth;
    let sigma_h = k0 * sigma_v;
    (sigma_v, sigma_h)
}

// ---------------------------------------------------------------------------
// Rock Quality Designation (RQD)
// ---------------------------------------------------------------------------

/// Estimate RQD from volumetric joint count Jv.
///
/// RQD ≈ 115 - 3.3·Jv  (Palmström, 1982), clamped to \[0, 100\].
pub fn rqd_from_jv(jv: f64) -> f64 {
    (115.0 - 3.3 * jv).clamp(0.0, 100.0)
}

/// Estimate RQD from the frequency of discontinuities per unit length λ.
///
/// RQD = 100 · e^(-0.1λ) · (0.1λ + 1)  (Priest & Hudson, 1976).
pub fn rqd_from_frequency(lambda: f64) -> f64 {
    100.0 * (-0.1 * lambda).exp() * (0.1 * lambda + 1.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // -- Stress utilities --

    #[test]
    fn test_mean_stress_hydrostatic() {
        let sigma = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        assert!((mean_stress(&sigma) - 100.0).abs() < EPS);
    }

    #[test]
    fn test_deviatoric_stress_hydrostatic_is_zero() {
        let sigma = [200.0, 200.0, 200.0, 0.0, 0.0, 0.0];
        let s = deviatoric_stress(&sigma);
        for &v in &s[0..3] {
            assert!(v.abs() < EPS, "deviatoric not zero: {}", v);
        }
    }

    #[test]
    fn test_von_mises_uniaxial() {
        // Pure uniaxial: σxx = 300, rest zero → q = 300
        let sigma = [300.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let q = von_mises_stress(&sigma);
        assert!((q - 300.0).abs() < 1e-6, "von Mises = {}", q);
    }

    #[test]
    fn test_principal_stresses_2d_pure_shear() {
        // Pure shear τ = 100 → σ1 = 100, σ2 = -100
        let (s1, s2, _theta) = principal_stresses_2d(0.0, 0.0, 100.0);
        assert!((s1 - 100.0).abs() < EPS);
        assert!((s2 + 100.0).abs() < EPS);
    }

    #[test]
    fn test_principal_stresses_2d_uniaxial() {
        let (s1, s2, _) = principal_stresses_2d(200.0, 0.0, 0.0);
        assert!((s1 - 200.0).abs() < EPS);
        assert!(s2.abs() < EPS);
    }

    // -- Linear elastic soil --

    #[test]
    fn test_linear_elastic_bulk_modulus() {
        let soil = LinearElasticSoil::new(1e6, 0.3);
        let k = soil.bulk_modulus();
        // E / (3(1-2*0.3)) = 1e6 / 1.2
        assert!((k - 1e6 / 1.2).abs() < 1.0);
    }

    #[test]
    fn test_linear_elastic_strain_stress_roundtrip() {
        let soil = LinearElasticSoil::new(2e6, 0.25);
        let eps = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma = soil.stress_from_strain(&eps);
        let eps2 = soil.strain_from_stress(&sigma);
        for k in 0..6 {
            assert!(
                (eps[k] - eps2[k]).abs() < 1e-10,
                "roundtrip failed at {}: {} vs {}",
                k,
                eps[k],
                eps2[k]
            );
        }
    }

    #[test]
    fn test_linear_elastic_settlement_positive() {
        let soil = LinearElasticSoil::new(5e6, 0.3);
        let s = soil.settlement_1d(50e3, 2.0);
        assert!(s > 0.0, "settlement = {}", s);
    }

    // -- Mohr-Coulomb --

    #[test]
    fn test_mohr_coulomb_shear_strength() {
        let mc = MohrCoulomb::new(10e3, PI / 6.0); // c=10kPa, φ=30°
        let tau = mc.shear_strength(100e3);
        // τ = 10000 + 100000 * tan(30°) ≈ 10000 + 57735
        assert!((tau - (10000.0 + 100000.0 * (PI / 6.0).tan())).abs() < 1e-4);
    }

    #[test]
    fn test_mohr_coulomb_yield_function_below_yield() {
        let mc = MohrCoulomb::new(20e3, 30.0_f64.to_radians());
        // Hydrostatic compression should not yield
        let f = mc.yield_function(100e3, 100e3);
        assert!(f <= 0.0, "yield function = {}", f);
    }

    #[test]
    fn test_mohr_coulomb_ka_kp_product() {
        let mc = MohrCoulomb::new(0.0, 30.0_f64.to_radians());
        let ka = mc.active_pressure_coefficient();
        let kp = mc.passive_pressure_coefficient();
        // Ka * Kp should be 1.0 for φ > 0
        assert!((ka * kp - 1.0).abs() < 1e-9, "Ka*Kp = {}", ka * kp);
    }

    #[test]
    fn test_mohr_coulomb_ucs_positive() {
        let mc = MohrCoulomb::new(25e3, 35.0_f64.to_radians());
        assert!(mc.ucs() > 0.0);
    }

    #[test]
    fn test_mohr_coulomb_to_drucker_prager() {
        let mc = MohrCoulomb::new(10e3, 30.0_f64.to_radians());
        let (alpha, k) = mc.to_drucker_prager_inner();
        assert!(alpha > 0.0, "alpha = {}", alpha);
        assert!(k > 0.0, "k = {}", k);
    }

    // -- Drucker-Prager --

    #[test]
    fn test_drucker_prager_yield_hydrostatic() {
        let mc = MohrCoulomb::new(20e3, 30.0_f64.to_radians());
        let dp = DruckerPragerModel::from_mohr_coulomb(&mc, 10e6, 0.3);
        let sigma = [100e3, 100e3, 100e3, 0.0, 0.0, 0.0];
        let f = dp.yield_function(&sigma);
        // Hydrostatic: deviatoric = 0, so F = alpha*I1 - k
        let expected = dp.alpha * 300e3 - dp.k;
        assert!((f - expected).abs() < 1e-4, "F = {}", f);
    }

    #[test]
    fn test_drucker_prager_return_mapping_stress_on_yield() {
        let dp = DruckerPragerModel::new(50e3, 0.1, 10e6, 0.3);
        // Trial stress well beyond yield
        let sigma_trial = [500e3, 200e3, 200e3, 50e3, 0.0, 0.0];
        let sigma_ret = dp.return_mapping(&sigma_trial);
        let f = dp.yield_function(&sigma_ret);
        assert!(f <= 1e-3, "after return mapping F = {}", f);
    }

    #[test]
    fn test_drucker_prager_elastic_no_change() {
        let dp = DruckerPragerModel::new(1e9, 0.0, 10e6, 0.3);
        // Very large k → no yielding
        let sigma = [10e3, 5e3, 5e3, 1e3, 0.0, 0.0];
        let sigma_ret = dp.return_mapping(&sigma);
        for k in 0..6 {
            assert!((sigma[k] - sigma_ret[k]).abs() < 1e-6);
        }
    }

    // -- Hoek-Brown --

    #[test]
    fn test_hoek_brown_uniaxial_compression() {
        let hb = HoekBrown::intact_rock(100e6, 7.0);
        // At σ3 = 0: σ1 = σci * s^a = 100e6 * 1.0^0.5 = 100e6
        let sigma1 = hb.sigma1_at_failure(0.0);
        assert!((sigma1 - 100e6).abs() < 1.0, "sigma1 = {}", sigma1);
    }

    #[test]
    fn test_hoek_brown_ucs_mass() {
        let hb = HoekBrown::new(100e6, 5.0, 0.25, 0.5);
        let ucs = hb.ucs_mass();
        assert!((ucs - 100e6 * 0.5).abs() < 1.0);
    }

    #[test]
    fn test_hoek_brown_tensile_strength_negative() {
        let hb = HoekBrown::intact_rock(50e6, 7.0);
        let st = hb.tensile_strength();
        // Should be negative (tension)
        assert!(st < 0.0, "tensile = {}", st);
    }

    // -- Griffith --

    #[test]
    fn test_griffith_tensile_failure() {
        let g = GriffithCriterion::new(10e6);
        // σ3 = -T0 → tensile failure
        assert!(g.is_failed(5e6, -10e6));
    }

    #[test]
    fn test_griffith_ucs() {
        let g = GriffithCriterion::new(10e6);
        assert!((g.ucs() - 80e6).abs() < 1.0);
    }

    #[test]
    fn test_griffith_no_failure_below_threshold() {
        let g = GriffithCriterion::new(5e6);
        // Highly compressive stress: should not fail in tension
        assert!(!g.is_failed(100e6, 50e6));
    }

    // -- Effective stress --

    #[test]
    fn test_effective_stress_basic() {
        let sigma_eff = effective_stress(200e3, 80e3);
        assert!((sigma_eff - 120e3).abs() < EPS);
    }

    #[test]
    fn test_hydrostatic_pore_pressure() {
        let u = hydrostatic_pore_pressure(10.0, 9810.0);
        assert!((u - 98100.0).abs() < EPS);
    }

    #[test]
    fn test_overburden_stress() {
        let layers = vec![(18000.0, 2.0), (20000.0, 3.0)];
        let sigma_v = overburden_stress(&layers);
        assert!((sigma_v - (36000.0 + 60000.0)).abs() < EPS);
    }

    // -- Terzaghi consolidation --

    #[test]
    fn test_terzaghi_doc_at_t_zero() {
        let c = TerzaghiConsolidation::new(1e-7, 5.0, 100e3);
        let u = c.average_degree_of_consolidation(0.0);
        assert!(u.abs() < 1e-8, "U(0) = {}", u);
    }

    #[test]
    fn test_terzaghi_doc_increases_with_time() {
        let c = TerzaghiConsolidation::new(1e-7, 3.0, 100e3);
        let u1 = c.average_degree_of_consolidation(1e7);
        let u2 = c.average_degree_of_consolidation(5e7);
        assert!(u2 > u1, "U should increase with time");
    }

    #[test]
    fn test_terzaghi_time_factor() {
        let c = TerzaghiConsolidation::new(1e-8, 2.0, 50e3);
        let tv = c.time_factor(1e8);
        assert!((tv - 1e-8 * 1e8 / 4.0).abs() < EPS);
    }

    #[test]
    fn test_terzaghi_settlement_monotone() {
        let c = TerzaghiConsolidation::new(1e-7, 5.0, 100e3);
        let s1 = c.settlement_at_time(1e6, 0.1);
        let s2 = c.settlement_at_time(1e8, 0.1);
        assert!(s2 >= s1);
    }

    // -- Biot poroelastic --

    #[test]
    fn test_biot_skempton_b_range() {
        let biot = BiotPoroelastic::new(1e9, 3e8, 0.8, 5e9, 1e-12, 1e-3);
        let b = biot.skempton_b();
        assert!((0.0..=1.0).contains(&b), "B = {}", b);
    }

    #[test]
    fn test_biot_undrained_bulk_ge_drained() {
        let biot = BiotPoroelastic::new(1e9, 3e8, 0.7, 2e9, 1e-12, 1e-3);
        let ku = biot.bulk_modulus_undrained();
        assert!(ku >= biot.bulk_modulus_drained, "Ku < Kd");
    }

    #[test]
    fn test_biot_effective_stress() {
        let biot = BiotPoroelastic::new(1e9, 3e8, 1.0, 2e9, 1e-12, 1e-3);
        let sigma = [200e3, 200e3, 200e3, 0.0, 0.0, 0.0];
        let sigma_eff = biot.effective_stress_voigt(&sigma, 50e3);
        assert!((sigma_eff[0] - 150e3).abs() < EPS);
    }

    // -- Liquefaction --

    #[test]
    fn test_liquefaction_csr_positive() {
        let liq = LiquefactionAssessment::new(0.3, 7.0, 150e3, 80e3, 10.0, 5.0);
        let csr = liq.cyclic_stress_ratio(5.0);
        assert!(csr > 0.0);
    }

    #[test]
    fn test_liquefaction_msf_gt_one_for_small_mw() {
        let liq = LiquefactionAssessment::new(0.3, 6.0, 150e3, 80e3, 10.0, 5.0);
        let msf = liq.magnitude_scaling_factor();
        assert!(msf > 1.0, "MSF = {}", msf);
    }

    #[test]
    fn test_liquefaction_high_n_not_liquefiable() {
        let liq = LiquefactionAssessment::new(0.2, 7.5, 120e3, 80e3, 35.0, 2.0);
        assert!(!liq.is_liquefiable(5.0));
    }

    // -- Kozeny-Carman --

    #[test]
    fn test_kozeny_carman_increases_with_void_ratio() {
        let k1 = kozeny_carman_permeability(0.5, 0.1e-3);
        let k2 = kozeny_carman_permeability(1.0, 0.1e-3);
        assert!(k2 > k1);
    }

    #[test]
    fn test_darcy_flow_direction() {
        // Positive gradient → negative flow (upstream to downstream)
        let q = darcy_flow(1e-12, 1e-3, 1000.0);
        assert!(q < 0.0);
    }

    // -- Compression index --

    #[test]
    fn test_cc_from_liquid_limit() {
        let cc = SoilCompressibility::cc_from_liquid_limit(50.0);
        assert!((cc - 0.009 * 40.0).abs() < EPS);
    }

    #[test]
    fn test_consolidation_settlement_nc() {
        let sc = SoilCompressibility::new(0.3, 0.05, 0.8, 50e3);
        // NC case: σ0=100kPa, σf=200kPa, σp=50kPa < σ0 → fully NC
        let s = sc.consolidation_settlement(5.0, 100e3, 200e3);
        assert!(s > 0.0, "settlement = {}", s);
    }

    #[test]
    fn test_consolidation_settlement_oc() {
        let sc = SoilCompressibility::new(0.3, 0.05, 0.8, 200e3);
        // OC case: σ0=50kPa, σf=100kPa, σp=200kPa > σf → recompression only
        let s = sc.consolidation_settlement(5.0, 50e3, 100e3);
        assert!(s > 0.0);
    }

    #[test]
    fn test_consolidation_settlement_cross_preconsolidation() {
        let sc = SoilCompressibility::new(0.3, 0.05, 0.8, 100e3);
        // Crosses: σ0=50kPa, σp=100kPa, σf=200kPa
        let s = sc.consolidation_settlement(5.0, 50e3, 200e3);
        let s_nc = sc.consolidation_settlement(5.0, 50e3, 200e3);
        assert!(s > 0.0 && s_nc > 0.0);
    }

    #[test]
    fn test_ocr_calculation() {
        let sc = SoilCompressibility::new(0.3, 0.05, 0.8, 200e3);
        let ocr = sc.ocr(100e3);
        assert!((ocr - 2.0).abs() < EPS);
    }

    // -- Undrained shear strength --

    #[test]
    fn test_shansep_increases_with_ocr() {
        let su1 = UndrainedShearStrength::shansep(100e3, 1.0);
        let su2 = UndrainedShearStrength::shansep(100e3, 2.0);
        assert!(su2 > su1);
    }

    #[test]
    fn test_su_from_triaxial() {
        let su = UndrainedShearStrength::from_triaxial(300e3, 100e3);
        assert!((su - 100e3).abs() < EPS);
    }

    // -- Bearing capacity --

    #[test]
    fn test_terzaghi_bearing_capacity_cohesive() {
        // φ=0: qu = c·Nc (Nc = π+2 for φ=0)
        let qu = terzaghi_bearing_capacity_strip(50e3, 0.0, 0.0, 1.0, 0.0);
        let expected = 50e3 * (PI + 2.0);
        assert!((qu - expected).abs() < 1e-2, "qu = {}", qu);
    }

    #[test]
    fn test_meyerhof_factors_phi30() {
        let (nc, nq, ng) = meyerhof_bearing_factors(30.0_f64.to_radians());
        assert!(nc > 0.0 && nq > 0.0 && ng > 0.0);
        assert!(nq > 1.0, "Nq should be > 1 for φ > 0");
    }

    // -- Stress path --

    #[test]
    fn test_triaxial_stress_path() {
        let (p, q) = triaxial_stress_path(300e3, 100e3);
        assert!((p - (300e3 + 200e3) / 3.0).abs() < EPS);
        assert!((q - 200e3).abs() < EPS);
    }

    #[test]
    fn test_k0_stress() {
        let (sv, sh) = k0_stress(18000.0, 5.0, 0.5);
        assert!((sv - 90000.0).abs() < EPS);
        assert!((sh - 45000.0).abs() < EPS);
    }

    // -- RQD --

    #[test]
    fn test_rqd_from_jv_clamped() {
        assert!((rqd_from_jv(0.0) - 115.0_f64.min(100.0)).abs() < EPS);
        assert!((rqd_from_jv(50.0)).abs() < EPS);
    }

    #[test]
    fn test_rqd_from_frequency_positive() {
        let rqd = rqd_from_frequency(10.0);
        assert!(rqd > 0.0 && rqd <= 100.0);
    }

    // -- Secondary compression --

    #[test]
    fn test_secondary_compression_increases_with_time() {
        let s1 = secondary_compression_settlement(0.02, 3.0, 1.2, 1e6, 1e7);
        let s2 = secondary_compression_settlement(0.02, 3.0, 1.2, 1e6, 1e8);
        assert!(s2 > s1);
    }

    #[test]
    fn test_secondary_compression_zero_for_backward_time() {
        let s = secondary_compression_settlement(0.02, 3.0, 1.2, 1e7, 1e6);
        assert_eq!(s, 0.0);
    }
}
