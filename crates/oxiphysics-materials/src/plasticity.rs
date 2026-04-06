// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Plasticity models for metals, soils, and general inelastic materials.
//!
//! Provides:
//! - Isotropic hardening (linear and nonlinear)
//! - Kinematic hardening (Prager model)
//! - Combined isotropic-kinematic hardening (Chaboche)
//! - J2 return mapping algorithm (radial return)
//! - Drucker-Prager yield surface (geomaterials)
//! - Modified Cam-Clay (critical state soil mechanics)
//! - Plastic dissipation computation
//! - Limit load (collapse load factor)

// ─── Isotropic Hardening ──────────────────────────────────────────────────────

/// Linear isotropic hardening law.
///
/// σ_y(ε_p) = σ_y0 + H * ε_p
///
/// where σ_y0 is the initial yield stress, H is the hardening modulus,
/// and ε_p is the equivalent plastic strain.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct IsotropicHardening {
    /// Initial yield stress σ_y0 (Pa).
    pub sigma_y0: f64,
    /// Isotropic hardening modulus H (Pa).
    pub hardening_modulus: f64,
}

#[allow(dead_code)]
impl IsotropicHardening {
    /// Create a new linear isotropic hardening model.
    pub fn new(sigma_y0: f64, hardening_modulus: f64) -> Self {
        Self {
            sigma_y0,
            hardening_modulus,
        }
    }

    /// Current yield stress at equivalent plastic strain ε_p.
    pub fn yield_stress(&self, eps_p: f64) -> f64 {
        self.sigma_y0 + self.hardening_modulus * eps_p
    }

    /// Derivative of yield stress with respect to ε_p (= H for linear hardening).
    pub fn d_yield_stress(&self) -> f64 {
        self.hardening_modulus
    }

    /// Elastic perfectly-plastic (H = 0).
    pub fn elastic_perfectly_plastic(sigma_y0: f64) -> Self {
        Self::new(sigma_y0, 0.0)
    }

    /// Mild steel default (σ_y0 = 250 MPa, H = 2 GPa).
    pub fn mild_steel() -> Self {
        Self::new(250.0e6, 2.0e9)
    }

    /// Aluminium alloy (σ_y0 = 270 MPa, H = 1.5 GPa).
    pub fn aluminium_alloy() -> Self {
        Self::new(270.0e6, 1.5e9)
    }
}

// ─── Nonlinear Isotropic Hardening (Voce) ─────────────────────────────────────

/// Voce nonlinear isotropic (saturation) hardening.
///
/// σ_y(ε_p) = σ_∞ - (σ_∞ - σ_y0) * exp(-b * ε_p)
///
/// Stress saturates to σ_∞ as ε_p → ∞.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct VoceHardening {
    /// Initial yield stress (Pa).
    pub sigma_y0: f64,
    /// Saturation yield stress (Pa).
    pub sigma_inf: f64,
    /// Saturation rate b (dimensionless).
    pub b: f64,
}

#[allow(dead_code)]
impl VoceHardening {
    /// Create a Voce hardening model.
    pub fn new(sigma_y0: f64, sigma_inf: f64, b: f64) -> Self {
        Self {
            sigma_y0,
            sigma_inf,
            b,
        }
    }

    /// Yield stress at equivalent plastic strain ε_p.
    pub fn yield_stress(&self, eps_p: f64) -> f64 {
        self.sigma_inf - (self.sigma_inf - self.sigma_y0) * (-self.b * eps_p).exp()
    }

    /// Derivative dσ_y/dε_p.
    pub fn d_yield_stress(&self, eps_p: f64) -> f64 {
        self.b * (self.sigma_inf - self.sigma_y0) * (-self.b * eps_p).exp()
    }

    /// Copper-like parameters (σ_y0 = 50 MPa, σ_∞ = 300 MPa, b = 10).
    pub fn copper_like() -> Self {
        Self::new(50.0e6, 300.0e6, 10.0)
    }
}

// ─── Kinematic Hardening (Prager) ─────────────────────────────────────────────

/// Linear kinematic (Prager) hardening model.
///
/// Backstress evolution: dα = C * dε_p
///
/// where C is the kinematic hardening modulus and α is the backstress.
/// The yield function becomes: f(σ - α) - σ_y0 = 0.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PragerKinematic {
    /// Initial yield stress σ_y0 (Pa).
    pub sigma_y0: f64,
    /// Kinematic hardening modulus C (Pa).
    pub c_kinematic: f64,
}

#[allow(dead_code)]
impl PragerKinematic {
    /// Create a Prager kinematic hardening model.
    pub fn new(sigma_y0: f64, c_kinematic: f64) -> Self {
        Self {
            sigma_y0,
            c_kinematic,
        }
    }

    /// Backstress increment Δα (Voigt notation: \[xx, yy, zz, xy, yz, xz\]).
    ///
    /// Δα = C * Δε_p
    pub fn backstress_increment(&self, delta_ep: &[f64; 6]) -> [f64; 6] {
        let c = self.c_kinematic;
        [
            c * delta_ep[0],
            c * delta_ep[1],
            c * delta_ep[2],
            c * delta_ep[3],
            c * delta_ep[4],
            c * delta_ep[5],
        ]
    }

    /// Effective stress: σ_eff = σ - α (shifted stress).
    pub fn effective_stress(stress: &[f64; 6], backstress: &[f64; 6]) -> [f64; 6] {
        [
            stress[0] - backstress[0],
            stress[1] - backstress[1],
            stress[2] - backstress[2],
            stress[3] - backstress[3],
            stress[4] - backstress[4],
            stress[5] - backstress[5],
        ]
    }

    /// Von Mises norm of a stress/backstress tensor (Voigt notation).
    pub fn von_mises_norm(s: &[f64; 6]) -> f64 {
        let mean = (s[0] + s[1] + s[2]) / 3.0;
        let dev = [s[0] - mean, s[1] - mean, s[2] - mean];
        let j2 = 0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
            + s[3].powi(2)
            + s[4].powi(2)
            + s[5].powi(2);
        (3.0 * j2).sqrt() // = sqrt(3 * J2) = von Mises
    }

    /// Yield function f = ||σ - α||_vM - σ_y0.
    pub fn yield_function(&self, stress: &[f64; 6], backstress: &[f64; 6]) -> f64 {
        let s_eff = Self::effective_stress(stress, backstress);
        Self::von_mises_norm(&s_eff) - self.sigma_y0
    }
}

// ─── Chaboche Combined Hardening ──────────────────────────────────────────────

/// Chaboche combined isotropic-kinematic hardening model.
///
/// Isotropic part: R(p) = Q*(1 - exp(-b*p)) where p = ∫|dε_p|
/// Kinematic part: dα = C * dε_p - γ * α * dp  (Armstrong-Frederick rule)
///
/// Yield surface: f = ||σ - α||_vM - (σ_y0 + R) = 0
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Chaboche {
    /// Initial yield stress (Pa).
    pub sigma_y0: f64,
    /// Isotropic saturation stress Q (Pa).
    pub q: f64,
    /// Isotropic saturation rate b.
    pub b_iso: f64,
    /// Kinematic hardening modulus C (Pa).
    pub c_kin: f64,
    /// Kinematic recovery parameter γ (dimensionless).
    pub gamma_kin: f64,
}

#[allow(dead_code)]
impl Chaboche {
    /// Create a Chaboche model.
    pub fn new(sigma_y0: f64, q: f64, b_iso: f64, c_kin: f64, gamma_kin: f64) -> Self {
        Self {
            sigma_y0,
            q,
            b_iso,
            c_kin,
            gamma_kin,
        }
    }

    /// Default steel parameters (representative values).
    pub fn steel_default() -> Self {
        Self::new(
            250.0e6, // σ_y0
            80.0e6,  // Q
            15.0,    // b
            50.0e9,  // C
            500.0,   // γ
        )
    }

    /// Isotropic hardening stress R(p) = Q*(1 - exp(-b*p)).
    pub fn isotropic_stress(&self, p: f64) -> f64 {
        self.q * (1.0 - (-self.b_iso * p).exp())
    }

    /// Current yield stress: σ_y0 + R(p).
    pub fn yield_stress(&self, p: f64) -> f64 {
        self.sigma_y0 + self.isotropic_stress(p)
    }

    /// Armstrong-Frederick backstress increment.
    ///
    /// Δα = C * n̂ * Δp - γ * α * Δp
    /// where n̂ is the flow direction, Δp the plastic strain increment.
    pub fn backstress_increment(
        &self,
        flow_dir: &[f64; 6],
        alpha: &[f64; 6],
        delta_p: f64,
    ) -> [f64; 6] {
        let mut da = [0.0f64; 6];
        for i in 0..6 {
            da[i] = (self.c_kin * flow_dir[i] - self.gamma_kin * alpha[i]) * delta_p;
        }
        da
    }

    /// Accumulated backstress norm ||α|| = sqrt(2/3 α:α).
    pub fn backstress_norm(alpha: &[f64; 6]) -> f64 {
        let sum = alpha[0].powi(2)
            + alpha[1].powi(2)
            + alpha[2].powi(2)
            + 2.0 * (alpha[3].powi(2) + alpha[4].powi(2) + alpha[5].powi(2));
        (2.0 / 3.0 * sum).sqrt()
    }
}

// ─── J2 Return Mapping ────────────────────────────────────────────────────────

/// J2 (von Mises) return mapping algorithm with isotropic hardening.
///
/// Given a trial elastic stress, compute the corrected plastic stress
/// via the radial return algorithm.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct J2ReturnMapping {
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
    /// Initial yield stress (Pa).
    pub sigma_y0: f64,
    /// Isotropic hardening modulus H (Pa).
    pub hardening_modulus: f64,
}

#[allow(dead_code)]
impl J2ReturnMapping {
    /// Create a J2 return mapping integrator.
    pub fn new(
        shear_modulus: f64,
        bulk_modulus: f64,
        sigma_y0: f64,
        hardening_modulus: f64,
    ) -> Self {
        Self {
            shear_modulus,
            bulk_modulus,
            sigma_y0,
            hardening_modulus,
        }
    }

    /// Create from Young's modulus and Poisson's ratio.
    #[allow(non_snake_case)]
    pub fn from_young_poisson(E: f64, nu: f64, sigma_y0: f64, h: f64) -> Self {
        let g = E / (2.0 * (1.0 + nu));
        let k = E / (3.0 * (1.0 - 2.0 * nu));
        Self::new(g, k, sigma_y0, h)
    }

    /// Von Mises stress (√(3 J₂)) from full Voigt stress \[xx,yy,zz,xy,yz,xz\].
    pub fn von_mises(stress: &[f64; 6]) -> f64 {
        let mean = (stress[0] + stress[1] + stress[2]) / 3.0;
        let dev = [stress[0] - mean, stress[1] - mean, stress[2] - mean];
        let j2 = 0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
            + stress[3].powi(2)
            + stress[4].powi(2)
            + stress[5].powi(2);
        (3.0 * j2).sqrt()
    }

    /// Radial return mapping.
    ///
    /// Returns (updated_stress, delta_eps_p, delta_eps_p_equivalent).
    pub fn return_map(&self, trial_stress: &[f64; 6], eps_p_bar: f64) -> ([f64; 6], [f64; 6], f64) {
        let sigma_y = self.sigma_y0 + self.hardening_modulus * eps_p_bar;
        let vm = Self::von_mises(trial_stress);

        if vm <= sigma_y + f64::EPSILON * sigma_y {
            // Elastic — no plastic correction
            return (*trial_stress, [0.0; 6], 0.0);
        }

        // Plastic consistency parameter Δγ
        let delta_gamma = (vm - sigma_y) / (3.0 * self.shear_modulus + self.hardening_modulus);

        // Deviatoric trial stress
        let p = (trial_stress[0] + trial_stress[1] + trial_stress[2]) / 3.0;
        let dev_trial = [
            trial_stress[0] - p,
            trial_stress[1] - p,
            trial_stress[2] - p,
            trial_stress[3],
            trial_stress[4],
            trial_stress[5],
        ];

        // Scale factor for deviatoric correction
        let scale = 1.0 - 3.0 * self.shear_modulus * delta_gamma / vm;

        // Updated stress (pressure unchanged)
        let updated_stress = [
            p + dev_trial[0] * scale,
            p + dev_trial[1] * scale,
            p + dev_trial[2] * scale,
            dev_trial[3] * scale,
            dev_trial[4] * scale,
            dev_trial[5] * scale,
        ];

        // Plastic strain increment (Δε_p = Δγ * n̂ where n̂ = dev/||dev||_vm)
        let n_hat = [
            dev_trial[0] / vm * 3.0 / 2.0,
            dev_trial[1] / vm * 3.0 / 2.0,
            dev_trial[2] / vm * 3.0 / 2.0,
            dev_trial[3] / vm * 3.0 / 2.0,
            dev_trial[4] / vm * 3.0 / 2.0,
            dev_trial[5] / vm * 3.0 / 2.0,
        ];
        let delta_eps_p = [
            delta_gamma * n_hat[0],
            delta_gamma * n_hat[1],
            delta_gamma * n_hat[2],
            delta_gamma * n_hat[3],
            delta_gamma * n_hat[4],
            delta_gamma * n_hat[5],
        ];

        (updated_stress, delta_eps_p, delta_gamma)
    }

    /// Consistent tangent modulus (elastic-plastic, uniaxial).
    ///
    /// C_ep = 2G * \[1 - 3G*Δγ/q_trial - 3G/(3G+H)*q_trial/(q_trial)\]
    /// Simplified scalar: E_ep = E*H/(E+H).
    pub fn consistent_tangent_uniaxial(&self) -> f64 {
        let e = 9.0 * self.bulk_modulus * self.shear_modulus
            / (3.0 * self.bulk_modulus + self.shear_modulus);
        e * self.hardening_modulus / (e + self.hardening_modulus)
    }
}

// ─── Drucker-Prager ───────────────────────────────────────────────────────────

/// Drucker-Prager yield surface for geomaterials (soil, rock, concrete).
///
/// Yield function: F = √J₂ + α·I₁ - k = 0
///
/// where:
///   α = 2 sin φ / (√3 (3 - sin φ))   (outer cone, Mohr-Coulomb match in compression)
///   k = 6 c cos φ / (√3 (3 - sin φ))
///   I₁ = tr(σ), J₂ = ½ dev:dev
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct DruckerPragerModel {
    /// Internal friction angle φ (radians).
    pub friction_angle: f64,
    /// Cohesion c (Pa).
    pub cohesion: f64,
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
}

#[allow(dead_code)]
impl DruckerPragerModel {
    /// Create a Drucker-Prager model.
    pub fn new(friction_angle: f64, cohesion: f64, shear_modulus: f64, bulk_modulus: f64) -> Self {
        Self {
            friction_angle,
            cohesion,
            shear_modulus,
            bulk_modulus,
        }
    }

    /// Create from φ (degrees), c (Pa), Young's modulus (Pa), Poisson's ratio.
    #[allow(non_snake_case, clippy::too_many_arguments)]
    pub fn from_mohr_coulomb_degrees(phi_deg: f64, cohesion: f64, E: f64, nu: f64) -> Self {
        let phi = phi_deg.to_radians();
        let g = E / (2.0 * (1.0 + nu));
        let k = E / (3.0 * (1.0 - 2.0 * nu));
        Self::new(phi, cohesion, g, k)
    }

    /// DP α parameter (friction / pressure coupling).
    pub fn alpha_dp(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi))
    }

    /// DP k parameter (cohesion term).
    pub fn k_dp(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let cos_phi = self.friction_angle.cos();
        6.0 * self.cohesion * cos_phi / (3.0_f64.sqrt() * (3.0 - sin_phi))
    }

    /// Second deviatoric stress invariant J₂.
    fn j2(stress: &[f64; 6]) -> f64 {
        let mean = (stress[0] + stress[1] + stress[2]) / 3.0;
        let d = [stress[0] - mean, stress[1] - mean, stress[2] - mean];
        0.5 * (d[0].powi(2) + d[1].powi(2) + d[2].powi(2))
            + stress[3].powi(2)
            + stress[4].powi(2)
            + stress[5].powi(2)
    }

    /// Yield function value F. F > 0 indicates yielding.
    pub fn yield_function(&self, stress: &[f64; 6]) -> f64 {
        let i1 = stress[0] + stress[1] + stress[2];
        let sqrt_j2 = Self::j2(stress).sqrt();
        sqrt_j2 + self.alpha_dp() * i1 - self.k_dp()
    }

    /// Check if stress state is elastic (F ≤ 0).
    pub fn is_elastic(&self, stress: &[f64; 6]) -> bool {
        self.yield_function(stress) <= 0.0
    }

    /// Uniaxial tensile strength (from cohesion and friction angle).
    ///
    /// σ_t = 2c cos φ / (1 + sin φ)
    pub fn uniaxial_tensile_strength(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let cos_phi = self.friction_angle.cos();
        2.0 * self.cohesion * cos_phi / (1.0 + sin_phi)
    }

    /// Uniaxial compressive strength.
    ///
    /// σ_c = 2c cos φ / (1 - sin φ)
    pub fn uniaxial_compressive_strength(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let cos_phi = self.friction_angle.cos();
        2.0 * self.cohesion * cos_phi / (1.0 - sin_phi)
    }
}

// ─── Modified Cam-Clay ────────────────────────────────────────────────────────

/// Modified Cam-Clay (MCC) critical state soil mechanics model.
///
/// Yield surface: f = q² / M² + p'(p' - p'_c) = 0
///
/// where:
///   p' = I₁/3 = mean effective stress
///   q = √(3 J₂) = deviatoric stress (von Mises)
///   M = critical state stress ratio q/p' at failure
///   p'_c = pre-consolidation pressure (hardening variable)
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct ModifiedCamClay {
    /// Critical state stress ratio M = q/p at critical state.
    pub m: f64,
    /// Swelling line slope κ in e-ln(p) space.
    pub kappa: f64,
    /// Normal consolidation line slope λ in e-ln(p) space.
    pub lambda: f64,
    /// Initial void ratio e₀.
    pub e0: f64,
    /// Poisson's ratio ν (for elastic stiffness).
    pub poisson_ratio: f64,
    /// Initial pre-consolidation pressure p'_c0 (Pa).
    pub pc0: f64,
}

#[allow(dead_code)]
impl ModifiedCamClay {
    /// Create a Modified Cam-Clay model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(m: f64, kappa: f64, lambda: f64, e0: f64, poisson_ratio: f64, pc0: f64) -> Self {
        Self {
            m,
            kappa,
            lambda,
            e0,
            poisson_ratio,
            pc0,
        }
    }

    /// Soft clay default parameters.
    pub fn soft_clay() -> Self {
        Self::new(0.9, 0.05, 0.2, 1.5, 0.3, 100.0e3)
    }

    /// Elastic bulk modulus K = (1+e₀) p' / κ.
    pub fn bulk_modulus(&self, p_prime: f64) -> f64 {
        (1.0 + self.e0) * p_prime / self.kappa
    }

    /// Elastic shear modulus from K and Poisson's ratio.
    pub fn shear_modulus(&self, p_prime: f64) -> f64 {
        let k = self.bulk_modulus(p_prime);
        3.0 * k * (1.0 - 2.0 * self.poisson_ratio) / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Mean effective stress from stress vector (Voigt: \[xx,yy,zz,xy,yz,xz\]).
    pub fn mean_stress(stress: &[f64; 6]) -> f64 {
        (stress[0] + stress[1] + stress[2]) / 3.0
    }

    /// Deviatoric stress q = √(3 J₂).
    pub fn deviatoric_stress(stress: &[f64; 6]) -> f64 {
        let p = Self::mean_stress(stress);
        let dev = [stress[0] - p, stress[1] - p, stress[2] - p];
        let j2 = 0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
            + stress[3].powi(2)
            + stress[4].powi(2)
            + stress[5].powi(2);
        (3.0 * j2).sqrt()
    }

    /// Yield function f = q²/M² + p'(p' - p'_c).
    pub fn yield_function(&self, stress: &[f64; 6], pc: f64) -> f64 {
        let p = Self::mean_stress(stress);
        let q = Self::deviatoric_stress(stress);
        q * q / (self.m * self.m) + p * (p - pc)
    }

    /// Critical state (CSL) pressure at given void ratio.
    ///
    /// p_cs = p_c / 2 (at yield surface apex for MCC)
    pub fn critical_state_pressure(&self, pc: f64) -> f64 {
        pc / 2.0
    }

    /// Hardening: evolution of p'_c with plastic volumetric strain Δε_v_p.
    ///
    /// p'_c^new = p'_c * exp((1+e₀) * Δε_v_p / (λ - κ))
    pub fn update_preconsolidation_pressure(&self, pc: f64, delta_eps_v_p: f64) -> f64 {
        pc * ((1.0 + self.e0) * delta_eps_v_p / (self.lambda - self.kappa)).exp()
    }

    /// Compression index Cc ≈ λ * ln(10) (conventional oedometer measure).
    pub fn compression_index(&self) -> f64 {
        self.lambda * 10.0_f64.ln()
    }
}

// ─── Plastic Dissipation ──────────────────────────────────────────────────────

/// Compute plastic dissipation D_p = σ : dε_p.
///
/// Both stress and plastic strain increment in Voigt notation \[xx,yy,zz,xy,yz,xz\].
/// Shear components have factor 2 for engineering notation (γ = 2 ε).
#[allow(dead_code)]
pub fn plastic_dissipation(stress: &[f64; 6], delta_eps_p: &[f64; 6]) -> f64 {
    // Normal components
    let normal =
        stress[0] * delta_eps_p[0] + stress[1] * delta_eps_p[1] + stress[2] * delta_eps_p[2];
    // Shear components (factor 2 for full tensor contraction)
    let shear = 2.0
        * (stress[3] * delta_eps_p[3] + stress[4] * delta_eps_p[4] + stress[5] * delta_eps_p[5]);
    normal + shear
}

/// Cumulative plastic dissipation over a loading history.
///
/// Returns ∑ σ_n : Δε_p_n summed over the provided increments.
#[allow(dead_code)]
pub fn cumulative_dissipation(
    stresses: &[[f64; 6]],
    plastic_strain_increments: &[[f64; 6]],
) -> f64 {
    stresses
        .iter()
        .zip(plastic_strain_increments.iter())
        .map(|(s, dep)| plastic_dissipation(s, dep))
        .sum()
}

// ─── Limit Load (Collapse Load Factor) ────────────────────────────────────────

/// Limit load factor estimation using the lower bound theorem.
///
/// For a structure carrying load F with yield stress σ_y and elastic stress σ_e,
/// the lower bound collapse load factor is:
///
///   f_L = σ_y / max(σ_vonMises)
///
/// Returns the fraction of reference load at which yielding first occurs.
#[allow(dead_code)]
pub fn limit_load_factor_lower_bound(
    elastic_stresses_voigt: &[[f64; 6]],
    yield_stress: f64,
) -> f64 {
    if elastic_stresses_voigt.is_empty() || yield_stress <= 0.0 {
        return 0.0;
    }
    let max_vm = elastic_stresses_voigt
        .iter()
        .map(von_mises_stress_voigt)
        .fold(0.0_f64, f64::max);
    if max_vm < f64::EPSILON {
        return f64::INFINITY;
    }
    yield_stress / max_vm
}

/// Von Mises stress from Voigt notation \[xx, yy, zz, xy, yz, xz\].
#[allow(dead_code)]
pub fn von_mises_stress_voigt(s: &[f64; 6]) -> f64 {
    let (sxx, syy, szz, sxy, syz, sxz) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    let term1 = (sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2);
    let term2 = 6.0 * (sxy.powi(2) + syz.powi(2) + sxz.powi(2));
    ((term1 + term2) / 2.0).sqrt()
}

// ─── Power-Law Creep-Plasticity Interaction ───────────────────────────────────

/// Power-law plasticity (Ramberg-Osgood type): ε_total = σ/E + (σ/σ_ref)^n / E.
///
/// Used to represent combined elastic-plastic behavior.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct RambergOsgood {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Yield reference stress σ_ref (Pa).
    pub sigma_ref: f64,
    /// Hardening exponent n (n ≥ 1, n=1 linear).
    pub n: f64,
    /// Offset coefficient α (often 3/7).
    pub alpha: f64,
}

#[allow(dead_code)]
impl RambergOsgood {
    /// Create a Ramberg-Osgood model.
    pub fn new(young_modulus: f64, sigma_ref: f64, n: f64, alpha: f64) -> Self {
        Self {
            young_modulus,
            sigma_ref,
            n,
            alpha,
        }
    }

    /// Total strain at given uniaxial stress σ:
    ///
    /// ε = σ/E + α * (σ/σ_ref)^n * σ_ref / E
    pub fn total_strain(&self, sigma: f64) -> f64 {
        sigma / self.young_modulus
            + self.alpha * (sigma / self.sigma_ref).powf(self.n) * self.sigma_ref
                / self.young_modulus
    }

    /// Plastic strain = total strain - elastic strain.
    pub fn plastic_strain(&self, sigma: f64) -> f64 {
        self.total_strain(sigma) - sigma / self.young_modulus
    }

    /// Secant modulus E_s = σ / ε_total.
    pub fn secant_modulus(&self, sigma: f64) -> f64 {
        let eps = self.total_strain(sigma);
        if eps.abs() < f64::EPSILON {
            return self.young_modulus;
        }
        sigma / eps
    }

    /// Tangent modulus E_t = dσ/dε (inverse of dε/dσ).
    pub fn tangent_modulus(&self, sigma: f64) -> f64 {
        let d_eps_d_sigma = 1.0 / self.young_modulus
            + self.alpha * self.n / self.sigma_ref * (sigma / self.sigma_ref).powf(self.n - 1.0);
        if d_eps_d_sigma.abs() < f64::EPSILON {
            return self.young_modulus;
        }
        1.0 / d_eps_d_sigma
    }
}

// ─── J2 Consistent Tangent Modulus ────────────────────────────────────────────

/// Full 6×6 algorithmic (consistent) tangent modulus for J2 plasticity
/// with isotropic hardening, in Voigt notation.
///
/// The consistent tangent ensures quadratic convergence in Newton-Raphson FEM.
/// It is defined as:
///
///   C_ep = C_e - (6G²)/(3G+H) * (n⊗n) / σ_y_trial²
///          - 2G * Δγ/q_tr * (I_dev - n⊗n)
///
/// where n = s_trial / q_tr is the unit deviatoric stress direction,
/// Δγ is the plastic consistency parameter, G is the shear modulus,
/// H is the hardening modulus, and I_dev is the deviatoric projector.
///
/// For the elastic step, returns the elastic stiffness C_e.
#[allow(dead_code)]
pub struct J2ConsistentTangent {
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
    /// Isotropic hardening modulus H (Pa).
    pub hardening_modulus: f64,
}

#[allow(dead_code)]
impl J2ConsistentTangent {
    /// Create a new J2 consistent tangent calculator.
    pub fn new(shear_modulus: f64, bulk_modulus: f64, hardening_modulus: f64) -> Self {
        Self {
            shear_modulus,
            bulk_modulus,
            hardening_modulus,
        }
    }

    /// Create from Young's modulus and Poisson's ratio.
    #[allow(non_snake_case)]
    pub fn from_young_poisson(E: f64, nu: f64, hardening_modulus: f64) -> Self {
        let g = E / (2.0 * (1.0 + nu));
        let k = E / (3.0 * (1.0 - 2.0 * nu));
        Self::new(g, k, hardening_modulus)
    }

    /// Elastic 6×6 stiffness in Voigt notation (isotropic).
    ///
    /// Voigt order: \[σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy\].
    pub fn elastic_stiffness(&self) -> [f64; 36] {
        let g = self.shear_modulus;
        let k = self.bulk_modulus;
        let lam = k - 2.0 / 3.0 * g; // λ = K - 2/3 G
        let mut c = [0.0_f64; 36];
        // Normal-normal coupling
        for i in 0..3 {
            for j in 0..3 {
                c[i * 6 + j] = lam;
            }
            c[i * 6 + i] += 2.0 * g;
        }
        // Shear diagonal
        c[3 * 6 + 3] = g;
        c[4 * 6 + 4] = g;
        c[5 * 6 + 5] = g;
        c
    }

    /// Compute the consistent (algorithmic) tangent modulus.
    ///
    /// Requires the trial stress, yield stress at start of increment, and
    /// whether plastic flow occurred.
    ///
    /// # Arguments
    /// * `trial_stress` — trial Voigt stress \[xx,yy,zz,yz,xz,xy\]
    /// * `sigma_y` — current yield stress (before increment)
    /// * `delta_gamma` — plastic consistency parameter (0 if elastic)
    ///
    /// # Returns
    /// 6×6 consistent tangent (flat row-major).
    pub fn compute_consistent_tangent(
        &self,
        trial_stress: &[f64; 6],
        delta_gamma: f64,
    ) -> [f64; 36] {
        let g = self.shear_modulus;
        let h = self.hardening_modulus;

        // Von Mises norm of trial deviatoric stress
        let p = (trial_stress[0] + trial_stress[1] + trial_stress[2]) / 3.0;
        let dev = [
            trial_stress[0] - p,
            trial_stress[1] - p,
            trial_stress[2] - p,
            trial_stress[3],
            trial_stress[4],
            trial_stress[5],
        ];
        let q_tr = {
            let j2 = 0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
                + dev[3].powi(2)
                + dev[4].powi(2)
                + dev[5].powi(2);
            (3.0 * j2).sqrt()
        };

        // If elastic (no plastic flow), return elastic stiffness
        if delta_gamma < f64::EPSILON || q_tr < f64::EPSILON {
            return self.elastic_stiffness();
        }

        // Flow direction n = sqrt(3/2) * s / q_tr  (unit deviatoric direction)
        let scale = if q_tr > 1e-30 { 1.0 / q_tr } else { 0.0 };
        let n: [f64; 6] = [
            dev[0] * scale,
            dev[1] * scale,
            dev[2] * scale,
            dev[3] * scale,
            dev[4] * scale,
            dev[5] * scale,
        ];

        // Start from elastic stiffness
        let mut c = self.elastic_stiffness();

        // Correction 1: -(6G²)/(3G+H) * (n⊗n)
        let coeff1 = 6.0 * g * g / (3.0 * g + h);
        for i in 0..6 {
            for j in 0..6 {
                c[i * 6 + j] -= coeff1 * n[i] * n[j];
            }
        }

        // Correction 2: -2G * Δγ/q_tr * (I_dev - n⊗n)
        // I_dev in Voigt: diagonal [2/3, 2/3, 2/3, 1, 1, 1] - off-diag [-1/3]
        let coeff2 = 2.0 * g * delta_gamma / q_tr;
        // I_dev diagonal components
        let i_dev_diag = [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0];
        for i in 0..3 {
            for j in 0..3 {
                let i_dev_ij = if i == j { i_dev_diag[i] } else { -1.0 / 3.0 };
                c[i * 6 + j] -= coeff2 * (i_dev_ij - n[i] * n[j]);
            }
        }
        for k in 3..6 {
            c[k * 6 + k] -= coeff2 * (i_dev_diag[k] - n[k] * n[k]);
        }

        c
    }
}

// ─── Drucker-Prager Apex Return Mapping ───────────────────────────────────────

/// Apex return mapping for the Drucker-Prager yield surface.
///
/// When the trial stress is at or above the apex (hydrostatic tensile side),
/// the stress must be mapped to the apex point, which is the intersection of
/// the yield surface with the hydrostatic axis.
///
/// Apex stress: p_apex = k / (3α),  q_apex = 0
///
/// where α and k are the DP parameters.
#[allow(dead_code)]
pub struct DruckerPragerApexReturn {
    /// DP friction parameter α.
    pub alpha: f64,
    /// DP cohesion parameter k (Pa).
    pub k_dp: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
}

#[allow(dead_code)]
impl DruckerPragerApexReturn {
    /// Create from DP parameters.
    pub fn new(alpha: f64, k_dp: f64, bulk_modulus: f64, shear_modulus: f64) -> Self {
        Self {
            alpha,
            k_dp,
            bulk_modulus,
            shear_modulus,
        }
    }

    /// Create from DruckerPragerModel.
    pub fn from_dp_model(dp: &DruckerPragerModel) -> Self {
        Self::new(dp.alpha_dp(), dp.k_dp(), dp.bulk_modulus, dp.shear_modulus)
    }

    /// Mean apex pressure p_apex = k / (3α).
    pub fn apex_pressure(&self) -> f64 {
        if self.alpha.abs() < f64::EPSILON {
            return 0.0;
        }
        self.k_dp / (3.0 * self.alpha)
    }

    /// Check whether the trial stress state lies above the apex.
    ///
    /// The trial mean stress is compared to p_apex. If p < p_apex (hydrostatic
    /// tension above apex), the apex return is required.
    pub fn needs_apex_return(&self, trial_stress: &[f64; 6]) -> bool {
        let p = (trial_stress[0] + trial_stress[1] + trial_stress[2]) / 3.0;
        p < self.apex_pressure()
    }

    /// Perform apex return: map trial stress to the apex.
    ///
    /// Returns the corrected stress state (hydrostatic, q = 0).
    pub fn compute_apex_return(&self, _trial_stress: &[f64; 6]) -> [f64; 6] {
        let p_apex = self.apex_pressure();
        [p_apex, p_apex, p_apex, 0.0, 0.0, 0.0]
    }

    /// Compute the plastic multiplier needed for apex return.
    ///
    /// From consistency: f(σ_apex) = 0 → Δλ = (√J2_tr + α * I1_tr - k) / (3αK + G)
    pub fn compute_apex_multiplier(&self, trial_stress: &[f64; 6]) -> f64 {
        let i1 = trial_stress[0] + trial_stress[1] + trial_stress[2];
        let p = i1 / 3.0;
        let dev = [
            trial_stress[0] - p,
            trial_stress[1] - p,
            trial_stress[2] - p,
        ];
        let j2 = 0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
            + trial_stress[3].powi(2)
            + trial_stress[4].powi(2)
            + trial_stress[5].powi(2);
        let sqrt_j2 = j2.sqrt();
        let f_tr = sqrt_j2 + self.alpha * i1 - self.k_dp;
        let denom = 3.0 * self.alpha * self.bulk_modulus + self.shear_modulus;
        if denom.abs() < f64::EPSILON {
            return 0.0;
        }
        f_tr / denom
    }
}

// ─── Cam-Clay Preconsolidation Evolution ──────────────────────────────────────

/// Cam-Clay preconsolidation pressure evolution model.
///
/// Models the expansion of the yield surface as the soil is loaded beyond
/// its current consolidation state. The preconsolidation pressure p'_c
/// hardening is governed by plastic volumetric strain.
///
/// Includes:
/// - Elastic wall (swelling line)
/// - Normal consolidation line (NCL)
/// - Critical state line (CSL)
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct CamClayPreconsolidation {
    /// Compression index Cc (slope of NCL in e-log₁₀p space).
    pub cc: f64,
    /// Swelling index Cs (slope of elastic line in e-log₁₀p space).
    pub cs: f64,
    /// Initial void ratio e0.
    pub e0: f64,
    /// Initial preconsolidation pressure p'_c0 (Pa).
    pub pc0: f64,
}

#[allow(dead_code)]
impl CamClayPreconsolidation {
    /// Create a new Cam-Clay preconsolidation model.
    pub fn new(cc: f64, cs: f64, e0: f64, pc0: f64) -> Self {
        Self { cc, cs, e0, pc0 }
    }

    /// Lightly overconsolidated clay default (Cc=0.5, Cs=0.05).
    pub fn soft_clay() -> Self {
        Self::new(0.5, 0.05, 1.2, 100.0e3)
    }

    /// Lambda parameter: λ = Cc / ln(10).
    pub fn lambda(&self) -> f64 {
        self.cc / 10.0_f64.ln()
    }

    /// Kappa parameter: κ = Cs / ln(10).
    pub fn kappa(&self) -> f64 {
        self.cs / 10.0_f64.ln()
    }

    /// Update preconsolidation pressure after plastic volumetric strain increment.
    ///
    /// p'_c_new = p'_c * exp((1 + e0) * Δε_v_p / (λ - κ))
    pub fn compute_preconsolidation(&self, pc_old: f64, delta_eps_v_p: f64) -> f64 {
        let lambda = self.lambda();
        let kappa = self.kappa();
        if (lambda - kappa).abs() < f64::EPSILON {
            return pc_old;
        }
        pc_old * ((1.0 + self.e0) * delta_eps_v_p / (lambda - kappa)).exp()
    }

    /// Void ratio on the normal consolidation line at pressure p.
    ///
    /// e_NCL(p) = e_N - Cc * log10(p / p_ref) where e_N = e0 at p = pc0.
    pub fn void_ratio_ncl(&self, p: f64) -> f64 {
        self.e0 - self.cc * (p / self.pc0).log10()
    }

    /// Void ratio on the swelling line at pressure p, given current state (pc, e_curr).
    ///
    /// e_sw(p) = e_curr - Cs * log10(p / p0)
    pub fn void_ratio_swelling(&self, p: f64, p0: f64, e_curr: f64) -> f64 {
        e_curr - self.cs * (p / p0).log10()
    }

    /// Overconsolidation ratio OCR = p'_c / p'.
    pub fn ocr(&self, p_prime: f64, pc: f64) -> f64 {
        if p_prime.abs() < f64::EPSILON {
            return 1.0;
        }
        pc / p_prime
    }

    /// Critical state mean stress for the current preconsolidation.
    ///
    /// p'_cs = p'_c / 2 (for Modified Cam-Clay).
    pub fn critical_state_pressure(&self, pc: f64) -> f64 {
        pc / 2.0
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── IsotropicHardening tests ───────────────────────────────────────────────

    #[test]
    fn test_isotropic_hardening_at_zero_strain() {
        let h = IsotropicHardening::mild_steel();
        let sy = h.yield_stress(0.0);
        assert!((sy - 250.0e6).abs() < 1.0, "σ_y(0) = {sy}");
    }

    #[test]
    fn test_isotropic_hardening_increases_with_strain() {
        let h = IsotropicHardening::mild_steel();
        let sy0 = h.yield_stress(0.0);
        let sy1 = h.yield_stress(0.01);
        assert!(
            sy1 > sy0,
            "yield stress should increase with plastic strain"
        );
    }

    #[test]
    fn test_isotropic_hardening_linear_rate() {
        let h = IsotropicHardening::new(200.0e6, 1.0e9);
        let eps_p = 0.005;
        let sy = h.yield_stress(eps_p);
        let expected = 200.0e6 + 1.0e9 * eps_p;
        assert!(
            (sy - expected).abs() < 1.0,
            "linear hardening: {sy} vs {expected}"
        );
    }

    #[test]
    fn test_isotropic_hardening_derivative_is_h() {
        let h = IsotropicHardening::new(200.0e6, 5.0e9);
        assert!((h.d_yield_stress() - 5.0e9).abs() < 1.0);
    }

    #[test]
    fn test_elastic_perfectly_plastic_constant_yield() {
        let h = IsotropicHardening::elastic_perfectly_plastic(300.0e6);
        assert!((h.yield_stress(0.0) - h.yield_stress(1.0)).abs() < 1.0);
    }

    // ── VoceHardening tests ────────────────────────────────────────────────────

    #[test]
    fn test_voce_at_zero_strain() {
        let v = VoceHardening::copper_like();
        let sy = v.yield_stress(0.0);
        assert!((sy - 50.0e6).abs() < 1.0, "Voce σ_y(0) = {sy}");
    }

    #[test]
    fn test_voce_saturates_at_large_strain() {
        let v = VoceHardening::copper_like();
        let sy = v.yield_stress(100.0);
        assert!(
            (sy - v.sigma_inf).abs() / v.sigma_inf < 0.001,
            "Voce should saturate: {sy}"
        );
    }

    #[test]
    fn test_voce_derivative_positive_at_small_strain() {
        let v = VoceHardening::new(100.0e6, 400.0e6, 5.0);
        assert!(v.d_yield_stress(0.01) > 0.0, "dσ_y/dε_p should be positive");
    }

    #[test]
    fn test_voce_derivative_decreases_with_strain() {
        let v = VoceHardening::new(100.0e6, 400.0e6, 5.0);
        let d1 = v.d_yield_stress(0.0);
        let d2 = v.d_yield_stress(0.1);
        assert!(d2 < d1, "hardening rate should decrease with strain");
    }

    // ── PragerKinematic tests ──────────────────────────────────────────────────

    #[test]
    fn test_prager_backstress_increment_direction() {
        let pk = PragerKinematic::new(250.0e6, 10.0e9);
        let dep = [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.0];
        let da = pk.backstress_increment(&dep);
        assert!((da[0] - 10.0e9 * 0.001).abs() < 1.0, "Δα_xx");
    }

    #[test]
    fn test_prager_von_mises_norm_uniaxial() {
        // Uniaxial stress: σ_xx = 300 MPa
        let s = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = PragerKinematic::von_mises_norm(&s);
        assert!((vm - 300.0e6).abs() < 1.0, "VM norm for uniaxial: {vm}");
    }

    #[test]
    fn test_prager_yield_function_below_yield() {
        let pk = PragerKinematic::new(300.0e6, 10.0e9);
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let backstress = [0.0; 6];
        let f = pk.yield_function(&stress, &backstress);
        assert!(f < 0.0, "Below yield: f = {f}");
    }

    #[test]
    fn test_prager_effective_stress_with_backstress() {
        let stress = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let backstress = [50.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let s_eff = PragerKinematic::effective_stress(&stress, &backstress);
        assert!((s_eff[0] - 250.0e6).abs() < 1.0, "σ_eff_xx = {}", s_eff[0]);
    }

    // ── Chaboche tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_chaboche_isotropic_stress_zero_at_zero_strain() {
        let c = Chaboche::steel_default();
        assert!(c.isotropic_stress(0.0).abs() < 1.0);
    }

    #[test]
    fn test_chaboche_isotropic_stress_saturates() {
        let c = Chaboche::steel_default();
        let r_large = c.isotropic_stress(1000.0);
        assert!(
            (r_large - c.q).abs() / c.q < 0.001,
            "Chaboche R should saturate: {r_large}"
        );
    }

    #[test]
    fn test_chaboche_yield_stress_increases_with_p() {
        let c = Chaboche::steel_default();
        let sy0 = c.yield_stress(0.0);
        let sy1 = c.yield_stress(0.01);
        assert!(sy1 > sy0, "yield stress should increase with p");
    }

    #[test]
    fn test_chaboche_backstress_norm_positive() {
        let alpha = [50.0e6, -20.0e6, -30.0e6, 5.0e6, 0.0, 0.0];
        let norm = Chaboche::backstress_norm(&alpha);
        assert!(norm > 0.0, "backstress norm should be positive: {norm}");
    }

    #[test]
    fn test_chaboche_backstress_increment_sign() {
        let c = Chaboche::steel_default();
        let flow_dir = [1.0, -0.5, -0.5, 0.0, 0.0, 0.0];
        let alpha = [0.0; 6];
        let da = c.backstress_increment(&flow_dir, &alpha, 0.001);
        // With α=0, Δα = C * n̂ * Δp > 0 for positive flow direction
        assert!(
            da[0] > 0.0,
            "backstress increment xx should be positive: {}",
            da[0]
        );
    }

    // ── J2ReturnMapping tests ──────────────────────────────────────────────────

    #[test]
    fn test_j2_return_map_elastic() {
        let rm = J2ReturnMapping::from_young_poisson(200.0e9, 0.3, 250.0e6, 2.0e9);
        let trial = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (stress, _dep, delta_gamma) = rm.return_map(&trial, 0.0);
        assert!(
            delta_gamma < f64::EPSILON,
            "Should be elastic: Δγ={delta_gamma}"
        );
        assert!(
            (stress[0] - 100.0e6).abs() < 1.0,
            "stress unchanged: {}",
            stress[0]
        );
    }

    #[test]
    fn test_j2_return_map_plastic() {
        let rm = J2ReturnMapping::from_young_poisson(200.0e9, 0.3, 250.0e6, 2.0e9);
        let trial = [500.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (stress, _dep, delta_gamma) = rm.return_map(&trial, 0.0);
        assert!(delta_gamma > 0.0, "Should yield: Δγ={delta_gamma}");
        // After return, von Mises should equal updated yield stress
        let vm = J2ReturnMapping::von_mises(&stress);
        let sy_updated = 250.0e6 + 2.0e9 * delta_gamma;
        assert!(
            (vm - sy_updated).abs() / sy_updated < 0.001,
            "VM after return: {vm} vs σ_y: {sy_updated}"
        );
    }

    #[test]
    fn test_j2_von_mises_uniaxial() {
        let s = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = J2ReturnMapping::von_mises(&s);
        assert!((vm - 300.0e6).abs() < 1.0, "VM uniaxial: {vm}");
    }

    #[test]
    fn test_j2_consistent_tangent_less_than_elastic() {
        let rm = J2ReturnMapping::from_young_poisson(200.0e9, 0.3, 250.0e6, 2.0e9);
        let c_ep = rm.consistent_tangent_uniaxial();
        let e =
            9.0 * rm.bulk_modulus * rm.shear_modulus / (3.0 * rm.bulk_modulus + rm.shear_modulus);
        assert!(c_ep < e, "C_ep should be < E");
        assert!(c_ep > 0.0, "C_ep should be positive");
    }

    // ── DruckerPragerModel tests ───────────────────────────────────────────────

    #[test]
    fn test_drucker_prager_cohesion_only_yield_stress() {
        // φ = 0: purely cohesive, α_dp = 0, k_dp = 2c/√3
        let dp = DruckerPragerModel::from_mohr_coulomb_degrees(0.0, 1.0e6, 30.0e9, 0.3);
        let k = dp.k_dp();
        let expected = 2.0 * 1.0e6 / 3.0_f64.sqrt();
        assert!(
            (k - expected).abs() / expected < 0.01,
            "k_dp = {k} vs {expected}"
        );
    }

    #[test]
    fn test_drucker_prager_hydrostatic_elastic() {
        let dp = DruckerPragerModel::from_mohr_coulomb_degrees(30.0, 1.0e6, 30.0e9, 0.3);
        // Large hydrostatic compression: I1 << 0, J2 = 0 → F = α*I1 - k < 0
        let stress = [-10.0e6, -10.0e6, -10.0e6, 0.0, 0.0, 0.0];
        assert!(
            dp.is_elastic(&stress),
            "Deep hydrostatic compression should be elastic"
        );
    }

    #[test]
    fn test_drucker_prager_tensile_strength_positive() {
        let dp = DruckerPragerModel::from_mohr_coulomb_degrees(30.0, 500.0e3, 10.0e9, 0.3);
        let st = dp.uniaxial_tensile_strength();
        assert!(st > 0.0, "tensile strength should be positive: {st}");
    }

    #[test]
    fn test_drucker_prager_compressive_strength_gt_tensile() {
        let dp = DruckerPragerModel::from_mohr_coulomb_degrees(30.0, 500.0e3, 10.0e9, 0.3);
        let st = dp.uniaxial_tensile_strength();
        let sc = dp.uniaxial_compressive_strength();
        assert!(sc > st, "compressive > tensile for φ>0: sc={sc}, st={st}");
    }

    // ── ModifiedCamClay tests ──────────────────────────────────────────────────

    #[test]
    fn test_mcc_bulk_modulus_increases_with_pressure() {
        let mcc = ModifiedCamClay::soft_clay();
        let k1 = mcc.bulk_modulus(50.0e3);
        let k2 = mcc.bulk_modulus(100.0e3);
        assert!(k2 > k1, "K should increase with p'");
    }

    #[test]
    fn test_mcc_yield_function_on_yield_surface() {
        let mcc = ModifiedCamClay::soft_clay();
        let pc = mcc.pc0;
        // At the right apex: q=0, p'=p'_c → f = 0/M² + p'_c*(p'_c - p'_c) = 0
        // Stress [pc, pc, pc, 0, 0, 0] → mean p' = pc, q = 0
        let stress_apex = [pc, pc, pc, 0.0, 0.0, 0.0];
        let f = mcc.yield_function(&stress_apex, pc);
        assert!(
            f.abs() < 1.0,
            "MCC yield f at right apex should be ~0, got {f}"
        );
    }

    #[test]
    fn test_mcc_hardening_increases_pc_under_compression() {
        let mcc = ModifiedCamClay::soft_clay();
        let pc_new = mcc.update_preconsolidation_pressure(mcc.pc0, 0.01);
        assert!(
            pc_new > mcc.pc0,
            "pc should increase under plastic volumetric compression"
        );
    }

    #[test]
    fn test_mcc_compression_index_positive() {
        let mcc = ModifiedCamClay::soft_clay();
        assert!(mcc.compression_index() > 0.0);
    }

    #[test]
    fn test_mcc_critical_state_pressure_half_pc() {
        let mcc = ModifiedCamClay::soft_clay();
        let p_cs = mcc.critical_state_pressure(mcc.pc0);
        assert!((p_cs - mcc.pc0 / 2.0).abs() < 1e-6);
    }

    // ── Plastic dissipation tests ──────────────────────────────────────────────

    #[test]
    fn test_plastic_dissipation_positive() {
        let stress = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let dep = [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.0];
        let d = plastic_dissipation(&stress, &dep);
        assert!(d > 0.0, "plastic dissipation should be positive: {d}");
    }

    #[test]
    fn test_plastic_dissipation_zero_for_zero_strain() {
        let stress = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let dep = [0.0; 6];
        let d = plastic_dissipation(&stress, &dep);
        assert!(d.abs() < 1e-10, "zero strain → zero dissipation: {d}");
    }

    #[test]
    fn test_cumulative_dissipation() {
        let stresses = vec![
            [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64],
            [310.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64],
        ];
        let deps = vec![
            [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.0_f64],
            [0.001, -0.0005, -0.0005, 0.0, 0.0, 0.0_f64],
        ];
        let d_total = cumulative_dissipation(&stresses, &deps);
        assert!(d_total > 0.0, "cumulative dissipation should be positive");
    }

    // ── Limit load tests ───────────────────────────────────────────────────────

    #[test]
    fn test_limit_load_factor_below_yield() {
        // Single element with σ_xx = 100 MPa, yield = 250 MPa → factor = 2.5
        let stresses = vec![[100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64]];
        let f = limit_load_factor_lower_bound(&stresses, 250.0e6);
        assert!((f - 2.5).abs() < 0.01, "limit load factor: {f}");
    }

    #[test]
    fn test_limit_load_factor_multiple_elements() {
        let stresses = vec![
            [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64],
            [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64], // governing
        ];
        let f = limit_load_factor_lower_bound(&stresses, 250.0e6);
        assert!((f - 1.25).abs() < 0.01, "limit load factor: {f}");
    }

    #[test]
    fn test_limit_load_empty_returns_zero() {
        let f = limit_load_factor_lower_bound(&[], 250.0e6);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_von_mises_voigt_uniaxial() {
        let s = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress_voigt(&s);
        assert!((vm - 200.0e6).abs() < 1.0, "VM uniaxial: {vm}");
    }

    #[test]
    fn test_von_mises_voigt_hydrostatic_zero() {
        // Hydrostatic stress → J2 = 0 → VM = 0
        let s = [100.0e6, 100.0e6, 100.0e6, 0.0, 0.0, 0.0];
        let vm = von_mises_stress_voigt(&s);
        assert!(vm.abs() < 1.0, "VM hydrostatic should be ~0: {vm}");
    }

    // ── RambergOsgood tests ────────────────────────────────────────────────────

    #[test]
    fn test_ramberg_osgood_elastic_at_small_stress() {
        let ro = RambergOsgood::new(200.0e9, 400.0e6, 5.0, 3.0 / 7.0);
        let sigma = 1.0e6; // very small
        let eps = ro.total_strain(sigma);
        let eps_elastic = sigma / ro.young_modulus;
        assert!(
            (eps - eps_elastic).abs() / eps_elastic < 0.01,
            "small stress should be nearly elastic"
        );
    }

    #[test]
    fn test_ramberg_osgood_plastic_strain_zero_at_zero_stress() {
        let ro = RambergOsgood::new(200.0e9, 400.0e6, 5.0, 3.0 / 7.0);
        assert!(ro.plastic_strain(0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_ramberg_osgood_tangent_modulus_less_than_elastic() {
        let ro = RambergOsgood::new(200.0e9, 400.0e6, 5.0, 3.0 / 7.0);
        let et = ro.tangent_modulus(ro.sigma_ref);
        assert!(
            et < ro.young_modulus,
            "tangent should be < E at reference stress"
        );
        assert!(et > 0.0, "tangent should be positive");
    }

    #[test]
    fn test_ramberg_osgood_secant_less_than_elastic() {
        let ro = RambergOsgood::new(200.0e9, 400.0e6, 5.0, 3.0 / 7.0);
        let es = ro.secant_modulus(ro.sigma_ref);
        assert!(
            es < ro.young_modulus,
            "secant modulus at ref stress should be < E"
        );
    }

    // ── J2ConsistentTangent tests ──────────────────────────────────────────────

    /// Elastic step: consistent tangent equals elastic stiffness.
    #[test]
    fn test_j2_consistent_tangent_elastic_step() {
        let ct = J2ConsistentTangent::from_young_poisson(200.0e9, 0.3, 2.0e9);
        let trial = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let c = ct.compute_consistent_tangent(&trial, 0.0);
        let ce = ct.elastic_stiffness();
        for i in 0..36 {
            assert!(
                (c[i] - ce[i]).abs() < 1.0,
                "Elastic tangent mismatch at {i}: {} vs {}",
                c[i],
                ce[i]
            );
        }
    }

    /// Elastic stiffness: C11 = K + 4G/3 for isotropic.
    #[test]
    fn test_j2_consistent_tangent_c11() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ct = J2ConsistentTangent::from_young_poisson(e, nu, 2.0e9);
        let c = ct.elastic_stiffness();
        let k = e / (3.0 * (1.0 - 2.0 * nu));
        let g = e / (2.0 * (1.0 + nu));
        let c11 = k + 4.0 / 3.0 * g;
        assert!((c[0] - c11).abs() / c11 < 1e-8, "C11: {} vs {}", c[0], c11);
    }

    /// Elastic stiffness: C12 = K - 2G/3 for isotropic.
    #[test]
    fn test_j2_consistent_tangent_c12() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ct = J2ConsistentTangent::from_young_poisson(e, nu, 2.0e9);
        let c = ct.elastic_stiffness();
        let k = e / (3.0 * (1.0 - 2.0 * nu));
        let g = e / (2.0 * (1.0 + nu));
        let c12 = k - 2.0 / 3.0 * g;
        assert!(
            (c[1] - c12).abs() / c12.abs() < 1e-8,
            "C12: {} vs {}",
            c[1],
            c12
        );
    }

    /// Plastic step: consistent tangent differs from elastic stiffness.
    #[test]
    fn test_j2_consistent_tangent_plastic_step_differs() {
        let ct = J2ConsistentTangent::from_young_poisson(200.0e9, 0.3, 2.0e9);
        // Uniaxial trial stress well beyond yield (so delta_gamma > 0)
        let trial = [500.0e6, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let delta_gamma = 1.0e-4;
        let c_ep = ct.compute_consistent_tangent(&trial, delta_gamma);
        let c_e = ct.elastic_stiffness();
        let diff: f64 = c_ep
            .iter()
            .zip(c_e.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(
            diff > 1.0,
            "Plastic tangent should differ from elastic: total diff = {diff}"
        );
    }

    // ── DruckerPragerApexReturn tests ─────────────────────────────────────────

    /// Apex pressure: p_apex = k / (3α).
    #[test]
    fn test_dp_apex_pressure() {
        let phi = 30.0_f64.to_radians();
        let cohesion = 50.0e3_f64;
        let e = 50.0e6_f64;
        let nu = 0.3_f64;
        let dp = DruckerPragerModel::new(
            phi,
            cohesion,
            e / (2.0 * (1.0 + nu)),
            e / (3.0 * (1.0 - 2.0 * nu)),
        );
        let apr = DruckerPragerApexReturn::from_dp_model(&dp);
        let p_apex = apr.apex_pressure();
        let expected = dp.k_dp() / (3.0 * dp.alpha_dp());
        assert!(
            (p_apex - expected).abs() / expected < 1e-10,
            "Apex pressure: {} vs {}",
            p_apex,
            expected
        );
    }

    /// Apex return stress has zero deviatoric part.
    #[test]
    fn test_dp_apex_return_zero_deviatoric() {
        let apr = DruckerPragerApexReturn::new(0.2, 100.0e3, 30.0e6, 15.0e6);
        let trial = [-200.0e3, -100.0e3, -50.0e3, 10.0e3, 0.0, 0.0];
        let s_apex = apr.compute_apex_return(&trial);
        // All three normal components should equal p_apex
        assert!(
            (s_apex[0] - s_apex[1]).abs() < f64::EPSILON,
            "Apex return should be hydrostatic"
        );
        assert!(
            (s_apex[1] - s_apex[2]).abs() < f64::EPSILON,
            "Apex return should be hydrostatic"
        );
        // All shear components zero
        assert_eq!(s_apex[3], 0.0);
    }

    // ── CamClayPreconsolidation tests ─────────────────────────────────────────

    /// Preconsolidation increases with compressive volumetric plastic strain.
    #[test]
    fn test_cam_clay_preconsolidation_increases_under_compression() {
        let cc = CamClayPreconsolidation::soft_clay();
        let pc_new = cc.compute_preconsolidation(cc.pc0, 0.01);
        assert!(
            pc_new > cc.pc0,
            "Preconsolidation should increase under compression"
        );
    }

    /// Zero plastic strain → unchanged preconsolidation.
    #[test]
    fn test_cam_clay_preconsolidation_zero_strain_unchanged() {
        let cc = CamClayPreconsolidation::soft_clay();
        let pc_new = cc.compute_preconsolidation(cc.pc0, 0.0);
        assert!(
            (pc_new - cc.pc0).abs() / cc.pc0 < 1e-12,
            "Zero plastic strain: pc should not change: {} vs {}",
            pc_new,
            cc.pc0
        );
    }

    /// OCR = 1 when p' equals pc (normally consolidated).
    #[test]
    fn test_cam_clay_ocr_normally_consolidated() {
        let cc = CamClayPreconsolidation::soft_clay();
        let ocr = cc.ocr(cc.pc0, cc.pc0);
        assert!(
            (ocr - 1.0).abs() < f64::EPSILON,
            "OCR should be 1 at p' = pc: {ocr}"
        );
    }

    /// Critical state pressure = pc / 2.
    #[test]
    fn test_cam_clay_critical_state_pressure() {
        let cc = CamClayPreconsolidation::soft_clay();
        let pcs = cc.critical_state_pressure(cc.pc0);
        assert!(
            (pcs - cc.pc0 / 2.0).abs() < f64::EPSILON,
            "p'_cs should be pc/2: {pcs}"
        );
    }
}
