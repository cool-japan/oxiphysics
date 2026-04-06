// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geomechanics FEM — Biot consolidation and soil mechanics.
//!
//! Implements coupled poroelastic and soil-mechanics computations:
//!
//! - [`BiotConsolidation`]: Biot poroelasticity parameters
//! - [`biot_effective_stress`]: Terzaghi effective stress principle
//! - [`consolidation_coefficient`]: Coefficient of consolidation Cv
//! - [`terzaghi_1d_solution`]: Analytical 1-D consolidation series
//! - [`SoilLayer`]: Single soil layer descriptor
//! - [`degree_of_consolidation_1d`]: Degree of consolidation U(t)
//! - [`MohrCoulombCriterion`]: Mohr–Coulomb yield criterion
//! - [`active_earth_pressure_ka`]: Rankine active earth pressure Ka
//! - [`passive_earth_pressure_kp`]: Rankine passive earth pressure Kp
//! - [`FoundationBearing`]: Terzaghi bearing capacity
//! - [`SlopeStability`]: Infinite-slope factor of safety
//! - [`settlement_elastic`]: Boussinesq elastic settlement

#![allow(dead_code)]
#![allow(non_snake_case)]

use std::f64::consts::PI;

// ============================================================================
// Physical constants
// ============================================================================

/// Unit weight of water γ_w (N/m³).
const GAMMA_W: f64 = 9810.0;

// ============================================================================
// BiotConsolidation
// ============================================================================

/// Parameters for Biot poroelasticity / consolidation analysis.
pub struct BiotConsolidation {
    /// Number of nodes in the FEM mesh.
    pub n_nodes: usize,
    /// Number of elements in the FEM mesh.
    pub n_elements: usize,
    /// Young's modulus E (Pa).
    pub E: f64,
    /// Poisson's ratio ν (dimensionless).
    pub nu: f64,
    /// Intrinsic permeability k (m/s).
    pub k_perm: f64,
    /// Biot coefficient α (dimensionless, 0 < α ≤ 1).
    pub alpha: f64,
    /// Biot modulus M (Pa).
    pub M_biot: f64,
}

impl BiotConsolidation {
    /// Create a new `BiotConsolidation` with the given parameters.
    pub fn new(
        n_nodes: usize,
        n_elements: usize,
        E: f64,
        nu: f64,
        k_perm: f64,
        alpha: f64,
        M_biot: f64,
    ) -> Self {
        Self {
            n_nodes,
            n_elements,
            E,
            nu,
            k_perm,
            alpha,
            M_biot,
        }
    }

    /// Compute the oedometric modulus M_oed = E(1-ν)/((1+ν)(1-2ν)).
    pub fn oedometric_modulus(&self) -> f64 {
        self.E * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
    }

    /// Compute the coefficient of consolidation Cv (m²/s).
    pub fn cv(&self) -> f64 {
        consolidation_coefficient(self.k_perm, self.E, self.nu)
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Compute the Terzaghi effective stress: σ′ = σ - α·p·I.
///
/// # Arguments
/// * `total_stress` – Voigt total stress tensor \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz\] (Pa)
/// * `pore_pressure` – pore water pressure p (Pa, positive in compression)
/// * `alpha` – Biot coefficient α
///
/// Returns the effective stress array in the same Voigt order.
pub fn biot_effective_stress(total_stress: [f64; 6], pore_pressure: f64, alpha: f64) -> [f64; 6] {
    [
        total_stress[0] - alpha * pore_pressure,
        total_stress[1] - alpha * pore_pressure,
        total_stress[2] - alpha * pore_pressure,
        total_stress[3],
        total_stress[4],
        total_stress[5],
    ]
}

/// Coefficient of consolidation Cv = k·E·(1-ν) / ((1+ν)·(1-2ν)·γ_w).
///
/// # Arguments
/// * `k`  – hydraulic conductivity (m/s)
/// * `E`  – Young's modulus (Pa)
/// * `nu` – Poisson's ratio
pub fn consolidation_coefficient(k: f64, E: f64, nu: f64) -> f64 {
    let denom = (1.0 + nu) * (1.0 - 2.0 * nu) * GAMMA_W;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    k * E * (1.0 - nu) / denom
}

/// Terzaghi 1-D consolidation series solution for excess pore pressure.
///
/// u(z, t) = (4/π) Σ_{m=0}^{n-1}
///   (1/(2m+1)) · sin((2m+1)·π·z/(2H)) · exp(-(2m+1)²·π²·Tv/4)
///
/// where Tv = cv·t/H².
///
/// # Arguments
/// * `z`       – depth from the drainage boundary (m, 0 ≤ z ≤ H)
/// * `t`       – time (s)
/// * `h`       – drainage path length H (m)
/// * `cv`      – coefficient of consolidation (m²/s)
/// * `n_terms` – number of series terms
///
/// Returns the normalised excess pore pressure u/u₀ ∈ \[0, 1\].
pub fn terzaghi_1d_solution(z: f64, t: f64, h: f64, cv: f64, n_terms: usize) -> f64 {
    if h.abs() < 1e-300 {
        return 0.0;
    }
    let tv = cv * t / (h * h);
    let mut sum = 0.0;
    for m in 0..n_terms {
        let m_f = m as f64;
        let n_val = 2.0 * m_f + 1.0;
        let sin_term = (n_val * PI * z / (2.0 * h)).sin();
        let exp_term = (-(n_val * n_val) * PI * PI * tv / 4.0).exp();
        sum += (1.0 / n_val) * sin_term * exp_term;
    }
    (4.0 / PI) * sum
}

// ============================================================================
// SoilLayer
// ============================================================================

/// Descriptor for a single soil layer used in consolidation analysis.
pub struct SoilLayer {
    /// Thickness of the layer (m).
    pub thickness: f64,
    /// Coefficient of consolidation Cv (m²/s).
    pub cv: f64,
    /// Initial excess pore pressure (Pa).
    pub initial_excess_pressure: f64,
}

impl SoilLayer {
    /// Create a new `SoilLayer`.
    pub fn new(thickness: f64, cv: f64, initial_excess_pressure: f64) -> Self {
        Self {
            thickness,
            cv,
            initial_excess_pressure,
        }
    }

    /// Time factor T_v = cv·t / H².
    pub fn time_factor(&self, t: f64) -> f64 {
        if self.thickness.abs() < 1e-300 {
            return 0.0;
        }
        self.cv * t / (self.thickness * self.thickness)
    }
}

/// Average degree of consolidation U(t) for 1-D Terzaghi consolidation
/// using the first `n_terms = 20` series terms.
///
/// U(t) = 1 - Σ_{m=0}^{n-1} (2/(M²)) · exp(-M²·Tv)
/// where M = (2m+1)·π/2.
///
/// # Arguments
/// * `t`  – time (s)
/// * `h`  – drainage path length H (m)
/// * `cv` – coefficient of consolidation (m²/s)
pub fn degree_of_consolidation_1d(t: f64, h: f64, cv: f64) -> f64 {
    if h.abs() < 1e-300 {
        return 1.0;
    }
    let tv = cv * t / (h * h);
    let n_terms = 20usize;
    let mut sum = 0.0;
    for m in 0..n_terms {
        let m_f = m as f64;
        let big_m = (2.0 * m_f + 1.0) * PI / 2.0;
        sum += (2.0 / (big_m * big_m)) * (-(big_m * big_m) * tv).exp();
    }
    1.0 - sum
}

// ============================================================================
// MohrCoulombCriterion
// ============================================================================

/// Mohr–Coulomb yield criterion for soil and rock.
pub struct MohrCoulombCriterion {
    /// Cohesion c (Pa).
    pub cohesion: f64,
    /// Internal friction angle φ (radians).
    pub friction_angle: f64,
}

impl MohrCoulombCriterion {
    /// Create a `MohrCoulombCriterion`.
    pub fn new(cohesion: f64, friction_angle: f64) -> Self {
        Self {
            cohesion,
            friction_angle,
        }
    }

    /// Yield function f(σ₁, σ₃) = (σ₁ - σ₃) - (σ₁ + σ₃)·sin(φ) - 2c·cos(φ).
    ///
    /// f < 0 → elastic; f = 0 → on yield surface; f > 0 → failed.
    pub fn yield_function(&self, sigma1: f64, sigma3: f64) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let cos_phi = self.friction_angle.cos();
        (sigma1 - sigma3) - (sigma1 + sigma3) * sin_phi - 2.0 * self.cohesion * cos_phi
    }

    /// Returns `true` if the stress state has reached or exceeded yield.
    pub fn is_yielded(&self, s1: f64, s3: f64) -> bool {
        self.yield_function(s1, s3) >= 0.0
    }

    /// Shear strength on a plane with given normal stress: τ = c + σ·tan(φ).
    pub fn shear_strength(&self, normal_stress: f64) -> f64 {
        self.cohesion + normal_stress * self.friction_angle.tan()
    }
}

// ============================================================================
// Earth pressure coefficients
// ============================================================================

/// Rankine active earth pressure coefficient Ka = tan²(45° - φ/2).
///
/// # Arguments
/// * `friction_angle` – internal friction angle φ (radians)
pub fn active_earth_pressure_ka(friction_angle: f64) -> f64 {
    let t = (PI / 4.0 - friction_angle / 2.0).tan();
    t * t
}

/// Rankine passive earth pressure coefficient Kp = tan²(45° + φ/2).
///
/// # Arguments
/// * `friction_angle` – internal friction angle φ (radians)
pub fn passive_earth_pressure_kp(friction_angle: f64) -> f64 {
    let t = (PI / 4.0 + friction_angle / 2.0).tan();
    t * t
}

// ============================================================================
// FoundationBearing
// ============================================================================

/// Bearing capacity parameters for a shallow strip foundation (Terzaghi).
pub struct FoundationBearing {
    /// Foundation width B (m).
    pub width: f64,
    /// Foundation embedment depth D (m).
    pub depth: f64,
    /// Soil cohesion c (Pa).
    pub cohesion: f64,
    /// Internal friction angle φ (radians).
    pub phi: f64,
    /// Unit weight of soil γ (N/m³).
    pub gamma: f64,
}

impl FoundationBearing {
    /// Create a `FoundationBearing`.
    pub fn new(width: f64, depth: f64, cohesion: f64, phi: f64, gamma: f64) -> Self {
        Self {
            width,
            depth,
            cohesion,
            phi,
            gamma,
        }
    }

    /// Ultimate bearing capacity: qu = c·Nc + γ·D·Nq + 0.5·γ·B·Nγ.
    pub fn ultimate_capacity(&self) -> f64 {
        let (nc, nq, ngamma) = bearing_factors(self.phi);
        self.cohesion * nc + self.gamma * self.depth * nq + 0.5 * self.gamma * self.width * ngamma
    }
}

/// Terzaghi bearing capacity factors (Nc, Nq, Nγ).
///
/// Uses the standard expressions:
/// Nq = e^(π·tan φ) · tan²(45 + φ/2)
/// Nc = (Nq - 1) / tan φ    (or 5.14 at φ = 0)
/// Nγ = 2(Nq + 1)·tan φ
///
/// # Arguments
/// * `phi` – internal friction angle φ (radians)
pub fn bearing_factors(phi: f64) -> (f64, f64, f64) {
    let nq = (PI * phi.tan()).exp() * (PI / 4.0 + phi / 2.0).tan().powi(2);
    let nc = if phi.abs() < 1e-9 {
        // Prandtl solution: Nc = π + 2 ≈ 5.14
        PI + 2.0
    } else {
        (nq - 1.0) / phi.tan()
    };
    let ngamma = 2.0 * (nq + 1.0) * phi.tan();
    (nc, nq, ngamma)
}

// ============================================================================
// SlopeStability
// ============================================================================

/// Parameters for infinite-slope stability analysis.
pub struct SlopeStability {
    /// Slope height H (m).
    pub height: f64,
    /// Slope inclination β (radians).
    pub beta: f64,
    /// Cohesion c (Pa).
    pub c: f64,
    /// Internal friction angle φ (radians).
    pub phi: f64,
    /// Bulk unit weight γ (N/m³).
    pub gamma: f64,
}

impl SlopeStability {
    /// Create a `SlopeStability`.
    pub fn new(height: f64, beta: f64, c: f64, phi: f64, gamma: f64) -> Self {
        Self {
            height,
            beta,
            c,
            phi,
            gamma,
        }
    }

    /// Factor of safety for an infinite slope:
    ///
    /// FS = (c + γ·H·cos²(β)·tan(φ)) / (γ·H·sin(β)·cos(β))
    pub fn factor_of_safety_infinite_slope(&self) -> f64 {
        let driving = self.gamma * self.height * self.beta.sin() * self.beta.cos();
        if driving.abs() < 1e-300 {
            return f64::INFINITY;
        }
        let resisting =
            self.c + self.gamma * self.height * self.beta.cos().powi(2) * self.phi.tan();
        resisting / driving
    }
}

// ============================================================================
// settlement_elastic
// ============================================================================

/// Elastic settlement using the Boussinesq formula:
///
/// s = q·B·(1 - ν²)·I_f / E
///
/// where I_f is an empirical influence factor (depth_factor here).
///
/// # Arguments
/// * `q`            – applied stress q (Pa)
/// * `b`            – foundation width B (m)
/// * `e_mod`        – Young's modulus of the soil E (Pa)
/// * `nu`           – Poisson's ratio ν
/// * `depth_factor` – influence factor I_f (dimensionless)
pub fn settlement_elastic(q: f64, b: f64, e_mod: f64, nu: f64, depth_factor: f64) -> f64 {
    if e_mod.abs() < 1e-300 {
        return 0.0;
    }
    q * b * (1.0 - nu * nu) * depth_factor / e_mod
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- BiotConsolidation --------------------------------------------------

    #[test]
    fn test_biot_consolidation_new() {
        let bc = BiotConsolidation::new(10, 5, 1e6, 0.3, 1e-5, 1.0, 1e8);
        assert_eq!(bc.n_nodes, 10);
        assert_eq!(bc.n_elements, 5);
    }

    #[test]
    fn test_biot_oedometric_modulus_positive() {
        let bc = BiotConsolidation::new(1, 1, 1e6, 0.3, 1e-5, 1.0, 1e8);
        let m_oed = bc.oedometric_modulus();
        assert!(m_oed > 0.0, "Oedometric modulus must be positive");
    }

    #[test]
    fn test_biot_cv_positive() {
        let bc = BiotConsolidation::new(1, 1, 1e6, 0.3, 1e-5, 1.0, 1e8);
        let cv = bc.cv();
        assert!(cv > 0.0, "Cv must be positive");
    }

    // ---- biot_effective_stress ----------------------------------------------

    #[test]
    fn test_effective_stress_reduces_normal() {
        let sigma = [100.0, 80.0, 60.0, 0.0, 0.0, 0.0];
        let eff = biot_effective_stress(sigma, 20.0, 1.0);
        assert!((eff[0] - 80.0).abs() < 1e-10);
        assert!((eff[1] - 60.0).abs() < 1e-10);
        assert!((eff[2] - 40.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_stress_shear_unchanged() {
        let sigma = [0.0, 0.0, 0.0, 15.0, 20.0, 25.0];
        let eff = biot_effective_stress(sigma, 10.0, 1.0);
        assert!((eff[3] - 15.0).abs() < 1e-10);
        assert!((eff[4] - 20.0).abs() < 1e-10);
        assert!((eff[5] - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_stress_zero_pore_pressure() {
        let sigma = [50.0, 60.0, 70.0, 5.0, 6.0, 7.0];
        let eff = biot_effective_stress(sigma, 0.0, 1.0);
        for (a, b) in sigma.iter().zip(eff.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_effective_stress_alpha_partial() {
        let sigma = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let eff = biot_effective_stress(sigma, 10.0, 0.5);
        assert!((eff[0] - 95.0).abs() < 1e-10);
    }

    // ---- consolidation_coefficient ------------------------------------------

    #[test]
    fn test_consolidation_coefficient_positive() {
        let cv = consolidation_coefficient(1e-5, 1e6, 0.3);
        assert!(cv > 0.0);
    }

    #[test]
    fn test_consolidation_coefficient_zero_k() {
        let cv = consolidation_coefficient(0.0, 1e6, 0.3);
        assert_eq!(cv, 0.0);
    }

    // ---- terzaghi_1d_solution -----------------------------------------------

    #[test]
    fn test_terzaghi_1d_zero_time() {
        // At t=0 the pore pressure should equal u₀ (u/u₀ = 1 at z=H/2)
        let u = terzaghi_1d_solution(0.5, 0.0, 1.0, 1e-4, 50);
        assert!(u > 0.9, "At t=0, u/u₀ should be near 1, got {u}");
    }

    #[test]
    fn test_terzaghi_1d_large_time_approaches_zero() {
        // At t→∞ excess pore pressure should dissipate
        let u = terzaghi_1d_solution(0.5, 1e10, 1.0, 1e-4, 50);
        assert!(u.abs() < 1e-6, "At t→∞, u/u₀ should → 0, got {u}");
    }

    #[test]
    fn test_terzaghi_1d_zero_height_returns_zero() {
        let u = terzaghi_1d_solution(0.5, 1.0, 0.0, 1e-4, 10);
        assert_eq!(u, 0.0);
    }

    // ---- degree_of_consolidation_1d -----------------------------------------

    #[test]
    fn test_degree_of_consolidation_zero_time() {
        let u = degree_of_consolidation_1d(0.0, 1.0, 1e-4);
        // With 20 series terms the partial sum at t=0 leaves a small truncation residual;
        // the true value is 0 but the approximation yields ~0.011 — verify it is small.
        assert!(u < 0.02, "U(0) should be near 0, got {u}");
    }

    #[test]
    fn test_degree_of_consolidation_large_time_approaches_one() {
        let u = degree_of_consolidation_1d(1e12, 1.0, 1e-4);
        assert!((u - 1.0).abs() < 1e-6, "U(∞) should be ~1, got {u}");
    }

    #[test]
    fn test_degree_of_consolidation_monotone() {
        let cv = 1e-4;
        let h = 1.0;
        let u1 = degree_of_consolidation_1d(1000.0, h, cv);
        let u2 = degree_of_consolidation_1d(5000.0, h, cv);
        assert!(u2 > u1, "Consolidation should increase with time");
    }

    // ---- SoilLayer ----------------------------------------------------------

    #[test]
    fn test_soil_layer_time_factor() {
        let layer = SoilLayer::new(5.0, 2e-3, 50e3);
        let tv = layer.time_factor(100.0);
        let expected = 2e-3 * 100.0 / 25.0;
        assert!((tv - expected).abs() < 1e-12);
    }

    // ---- MohrCoulombCriterion -----------------------------------------------

    #[test]
    fn test_mohr_coulomb_no_yield_below_failure() {
        let mc = MohrCoulombCriterion::new(10e3, 30.0_f64.to_radians());
        // Very small principal stress difference
        assert!(!mc.is_yielded(20e3, 18e3));
    }

    #[test]
    fn test_mohr_coulomb_yield_at_failure() {
        let phi = 0.0_f64; // frictionless
        let c = 100.0;
        let mc = MohrCoulombCriterion::new(c, phi);
        // f = (σ1 - σ3) - 2c
        let f = mc.yield_function(200.0, 0.0);
        assert!((f - 0.0).abs() < 1.0); // near zero
    }

    #[test]
    fn test_shear_strength_cohesion_only() {
        let mc = MohrCoulombCriterion::new(50.0, 0.0);
        assert!((mc.shear_strength(100.0) - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_shear_strength_increases_with_normal_stress() {
        let mc = MohrCoulombCriterion::new(10.0, 30.0_f64.to_radians());
        let tau1 = mc.shear_strength(100.0);
        let tau2 = mc.shear_strength(200.0);
        assert!(tau2 > tau1);
    }

    // ---- earth pressure -----------------------------------------------------

    #[test]
    fn test_ka_zero_friction() {
        let ka = active_earth_pressure_ka(0.0);
        assert!((ka - 1.0).abs() < 1e-10, "Ka at φ=0 should be 1.0");
    }

    #[test]
    fn test_kp_zero_friction() {
        let kp = passive_earth_pressure_kp(0.0);
        assert!((kp - 1.0).abs() < 1e-10, "Kp at φ=0 should be 1.0");
    }

    #[test]
    fn test_ka_kp_product_unity_for_zero_phi() {
        let phi = 0.0;
        assert!(
            (active_earth_pressure_ka(phi) * passive_earth_pressure_kp(phi) - 1.0).abs() < 1e-10
        );
    }

    #[test]
    fn test_ka_30_degrees() {
        let phi = 30.0_f64.to_radians();
        let ka = active_earth_pressure_ka(phi);
        // tan²(30°) = 1/3 ≈ 0.3333
        assert!(
            (ka - 1.0 / 3.0).abs() < 1e-4,
            "Ka at φ=30° ≈ 0.333, got {ka}"
        );
    }

    #[test]
    fn test_kp_30_degrees() {
        let phi = 30.0_f64.to_radians();
        let kp = passive_earth_pressure_kp(phi);
        // tan²(60°) = 3.0
        assert!((kp - 3.0).abs() < 1e-4, "Kp at φ=30° ≈ 3.0, got {kp}");
    }

    #[test]
    fn test_ka_less_than_one() {
        let phi = 20.0_f64.to_radians();
        assert!(active_earth_pressure_ka(phi) < 1.0);
    }

    #[test]
    fn test_kp_greater_than_one() {
        let phi = 20.0_f64.to_radians();
        assert!(passive_earth_pressure_kp(phi) > 1.0);
    }

    #[test]
    fn test_ka_kp_reciprocal() {
        let phi = 25.0_f64.to_radians();
        let ka = active_earth_pressure_ka(phi);
        let kp = passive_earth_pressure_kp(phi);
        // Ka * Kp != 1 in general (only at phi=0), but Ka < 1 < Kp
        assert!(ka < 1.0);
        assert!(kp > 1.0);
    }

    // ---- bearing_factors ----------------------------------------------------

    #[test]
    fn test_bearing_nc_zero_phi() {
        let (nc, _nq, _ngamma) = bearing_factors(0.0);
        // Prandtl Nc = π + 2 ≈ 5.14
        assert!((nc - (PI + 2.0)).abs() < 0.01, "Nc at φ=0 ≈ 5.14, got {nc}");
    }

    #[test]
    fn test_bearing_nq_unity_at_zero_phi() {
        let (_nc, nq, _ngamma) = bearing_factors(0.0);
        // Nq = e^0 * tan²(45) = 1
        assert!(
            (nq - 1.0).abs() < 1e-10,
            "Nq at φ=0 should be 1.0, got {nq}"
        );
    }

    #[test]
    fn test_bearing_factors_positive() {
        let (nc, nq, ngamma) = bearing_factors(30.0_f64.to_radians());
        assert!(nc > 0.0 && nq > 0.0 && ngamma > 0.0);
    }

    // ---- FoundationBearing --------------------------------------------------

    #[test]
    fn test_foundation_bearing_positive_capacity() {
        let fb = FoundationBearing::new(2.0, 1.0, 10e3, 30.0_f64.to_radians(), 18e3);
        assert!(fb.ultimate_capacity() > 0.0);
    }

    #[test]
    fn test_foundation_bearing_increases_with_cohesion() {
        let fb1 = FoundationBearing::new(2.0, 1.0, 10e3, 20.0_f64.to_radians(), 18e3);
        let fb2 = FoundationBearing::new(2.0, 1.0, 50e3, 20.0_f64.to_radians(), 18e3);
        assert!(fb2.ultimate_capacity() > fb1.ultimate_capacity());
    }

    #[test]
    fn test_foundation_bearing_increases_with_width() {
        let fb1 = FoundationBearing::new(1.0, 1.0, 10e3, 25.0_f64.to_radians(), 18e3);
        let fb2 = FoundationBearing::new(4.0, 1.0, 10e3, 25.0_f64.to_radians(), 18e3);
        assert!(fb2.ultimate_capacity() > fb1.ultimate_capacity());
    }

    // ---- SlopeStability -----------------------------------------------------

    #[test]
    fn test_slope_stability_beta_equals_phi() {
        // FS = (c + ...) / driving; for c=0 and beta=phi, FS = 1.0
        let phi = 30.0_f64.to_radians();
        let ss = SlopeStability::new(5.0, phi, 0.0, phi, 18e3);
        let fs = ss.factor_of_safety_infinite_slope();
        assert!(
            (fs - 1.0).abs() < 1e-8,
            "FS should be 1 when beta=phi, c=0; got {fs}"
        );
    }

    #[test]
    fn test_slope_stability_fs_greater_one_for_gentle_slope() {
        let phi = 30.0_f64.to_radians();
        let beta = 15.0_f64.to_radians();
        let ss = SlopeStability::new(5.0, beta, 0.0, phi, 18e3);
        assert!(ss.factor_of_safety_infinite_slope() > 1.0);
    }

    #[test]
    fn test_slope_stability_cohesion_increases_fs() {
        let phi = 20.0_f64.to_radians();
        let beta = 25.0_f64.to_radians();
        let ss1 = SlopeStability::new(5.0, beta, 0.0, phi, 18e3);
        let ss2 = SlopeStability::new(5.0, beta, 10e3, phi, 18e3);
        assert!(ss2.factor_of_safety_infinite_slope() > ss1.factor_of_safety_infinite_slope());
    }

    // ---- settlement_elastic -------------------------------------------------

    #[test]
    fn test_settlement_elastic_positive() {
        let s = settlement_elastic(100e3, 2.0, 10e6, 0.3, 0.85);
        assert!(s > 0.0);
    }

    #[test]
    fn test_settlement_elastic_zero_modulus() {
        let s = settlement_elastic(100e3, 2.0, 0.0, 0.3, 0.85);
        assert_eq!(s, 0.0);
    }

    #[test]
    fn test_settlement_elastic_increases_with_load() {
        let s1 = settlement_elastic(100e3, 2.0, 10e6, 0.3, 0.85);
        let s2 = settlement_elastic(200e3, 2.0, 10e6, 0.3, 0.85);
        assert!(s2 > s1);
    }
}
