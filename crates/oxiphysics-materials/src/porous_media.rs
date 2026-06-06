// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Porous and cellular material models.
//!
//! Provides constitutive relations for porous solids, including Darcy/Forchheimer
//! flow through porous media, the Kozeny-Carman permeability model, effective
//! thermal conductivity, Gibson-Ashby cellular solid stiffness, and Plateau border
//! geometry for open-cell foams.

/// Porous medium descriptor.
#[derive(Debug, Clone, PartialEq)]
pub struct PorousMedia {
    /// Porosity φ – volume fraction of void space (0–1).
    pub porosity: f64,
    /// Intrinsic permeability κ \[m²\].
    pub permeability: f64,
    /// Tortuosity τ of the pore channels (≥ 1).
    pub tortuosity: f64,
    /// Dynamic viscosity μ of the saturating fluid \[Pa·s\].
    pub viscosity: f64,
    /// Fluid density ρ \[kg/m³\] (used in Forchheimer correction).
    pub fluid_density: f64,
}

impl PorousMedia {
    /// Create a new `PorousMedia`.
    pub fn new(
        porosity: f64,
        permeability: f64,
        tortuosity: f64,
        viscosity: f64,
        fluid_density: f64,
    ) -> Self {
        Self {
            porosity,
            permeability,
            tortuosity,
            viscosity,
            fluid_density,
        }
    }
}

/// Darcy's law — superficial velocity \[m/s\].
///
/// Returns the Darcy flux q = −(κ/μ) · ∇P in 1-D (scalar form).
///
/// # Arguments
/// * `medium` – porous medium properties
/// * `grad_p` – pressure gradient ∂P/∂x \[Pa/m\]  (negative for flow in +x)
pub fn darcy_flow(medium: &PorousMedia, grad_p: f64) -> f64 {
    -(medium.permeability / medium.viscosity) * grad_p
}

/// Forchheimer correction — superficial velocity including inertial term.
///
/// Solves the quadratic `β·ρ·q² + μ/κ·q + ∇P = 0` for q.
/// The Forchheimer coefficient β \[1/m\] accounts for inertial (high-Re) effects.
///
/// For a given pressure gradient only the physically meaningful (Darcy-like)
/// root is returned.
///
/// # Arguments
/// * `medium` – porous medium properties
/// * `grad_p` – pressure gradient \[Pa/m\]
/// * `beta`   – Forchheimer inertial resistance coefficient \[1/m\]
pub fn forchheimer_correction(medium: &PorousMedia, grad_p: f64, beta: f64) -> f64 {
    // β ρ q² + (μ/κ) q + ∇P = 0   →   q = [ -μ/κ ± sqrt((μ/κ)² - 4β ρ ∇P) ] / (2β ρ)
    //
    // When grad_p < 0 (pressure drives flow in +x), ∇P < 0, discriminant > 0.
    let a = beta * medium.fluid_density;
    let b_coeff = medium.viscosity / medium.permeability;
    let c = grad_p;

    if a.abs() < 1e-30 {
        // Degenerate: no inertial term → Darcy.
        return darcy_flow(medium, grad_p);
    }

    let discriminant = b_coeff * b_coeff - 4.0 * a * c;
    if discriminant < 0.0 {
        // No real solution – return Darcy approximation as fallback.
        return darcy_flow(medium, grad_p);
    }

    // Choose the root with the same sign as the Darcy solution.
    let sqrt_d = discriminant.sqrt();
    let q1 = (-b_coeff + sqrt_d) / (2.0 * a);
    let q2 = (-b_coeff - sqrt_d) / (2.0 * a);

    // The physically correct root has smaller magnitude (closer to Darcy).
    if q1.abs() < q2.abs() { q1 } else { q2 }
}

/// Kozeny-Carman permeability model.
///
/// Returns the intrinsic permeability κ \[m²\] for a packed-bed of
/// approximately spherical particles:
///
/// κ = φ³ · d_p² / (180 · (1−φ)²)
///
/// # Arguments
/// * `porosity`       – void fraction φ (0–1)
/// * `particle_diam`  – effective particle diameter d_p \[m\]
pub fn kozeny_carman(porosity: f64, particle_diam: f64) -> f64 {
    let phi = porosity;
    let omp = 1.0 - phi;
    phi.powi(3) * particle_diam.powi(2) / (180.0 * omp * omp)
}

/// Effective thermal conductivity — parallel and series mixture models.
///
/// Returns `(k_parallel, k_series)` where:
/// - `k_parallel` is the upper (Voigt / parallel) bound.
/// - `k_series`   is the lower (Reuss / series) bound.
///
/// # Arguments
/// * `k_solid`  – thermal conductivity of the solid phase \[W/(m·K)\]
/// * `k_fluid`  – thermal conductivity of the fluid/void phase \[W/(m·K)\]
/// * `porosity` – void fraction φ (0–1)
pub fn effective_thermal_conductivity(k_solid: f64, k_fluid: f64, porosity: f64) -> (f64, f64) {
    let phi = porosity;
    let k_par = (1.0 - phi) * k_solid + phi * k_fluid;
    let k_ser = 1.0 / ((1.0 - phi) / k_solid + phi / k_fluid);
    (k_par, k_ser)
}

/// Gibson-Ashby relative stiffness for cellular (open or closed-cell) solids.
///
/// Relative stiffness: E*/E_s = C · (ρ*/ρ_s)ⁿ
///
/// For open-cell foams: C ≈ 1.0, n = 2.
/// For closed-cell foams additional membrane stretching contributes (n ≈ 2, C adjusted).
///
/// # Arguments
/// * `relative_density` – ρ*/ρ_s (ratio of foam to solid density)
/// * `c_coeff`          – Gibson-Ashby proportionality constant (≈ 1.0)
/// * `exponent`         – power-law exponent (2 for open-cell, ~1.5–2 for closed-cell)
pub fn gibson_ashby_stiffness(relative_density: f64, c_coeff: f64, exponent: f64) -> f64 {
    c_coeff * relative_density.powf(exponent)
}

/// Open-cell Plateau border cross-sectional area \[m²\].
///
/// The Plateau border is the liquid-filled channel at junctions of three soap
/// films in an open-cell foam.  Its cross-sectional area is:
///
/// A_pb = (√3 − π/2) · r²
///
/// where r is the border radius of curvature \[m\].
///
/// # Arguments
/// * `radius` – Plateau border radius \[m\]
pub fn open_cell_plateau(radius: f64) -> f64 {
    (3.0_f64.sqrt() - std::f64::consts::PI / 2.0) * radius * radius
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    // ------------------------------------------------------------------
    // PorousMedia construction
    // ------------------------------------------------------------------

    #[test]
    fn test_porous_media_new() {
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        assert_eq!(m.porosity, 0.4);
        assert_eq!(m.permeability, 1e-12);
        assert_eq!(m.tortuosity, 1.5);
        assert_eq!(m.viscosity, 1e-3);
        assert_eq!(m.fluid_density, 1000.0);
    }

    // ------------------------------------------------------------------
    // Darcy flow
    // ------------------------------------------------------------------

    #[test]
    fn test_darcy_flow_positive_direction() {
        // Negative pressure gradient → positive flux (flow in +x).
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let q = darcy_flow(&m, -1000.0);
        assert!(
            q > 0.0,
            "Darcy flux should be positive for negative grad_p, got {q}"
        );
    }

    #[test]
    fn test_darcy_flow_zero_gradient() {
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let q = darcy_flow(&m, 0.0);
        assert!(q.abs() < EPS, "Zero pressure gradient → zero flux, got {q}");
    }

    #[test]
    fn test_darcy_flow_formula() {
        let m = PorousMedia::new(0.4, 2e-12, 1.5, 2e-3, 1000.0);
        let grad_p = -500.0;
        let q = darcy_flow(&m, grad_p);
        let expected = -(m.permeability / m.viscosity) * grad_p;
        assert!((q - expected).abs() < EPS * expected.abs().max(1.0));
    }

    #[test]
    fn test_darcy_flow_scales_with_permeability() {
        let m1 = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let m2 = PorousMedia::new(0.4, 2e-12, 1.5, 1e-3, 1000.0);
        let grad_p = -200.0;
        let q1 = darcy_flow(&m1, grad_p);
        let q2 = darcy_flow(&m2, grad_p);
        assert!(
            (q2 - 2.0 * q1).abs() < EPS,
            "Doubling κ should double Darcy flux"
        );
    }

    #[test]
    fn test_darcy_flow_scales_with_viscosity() {
        let m1 = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let m2 = PorousMedia::new(0.4, 1e-12, 1.5, 2e-3, 1000.0);
        let grad_p = -200.0;
        let q1 = darcy_flow(&m1, grad_p);
        let q2 = darcy_flow(&m2, grad_p);
        assert!(
            (q2 - q1 / 2.0).abs() < EPS,
            "Doubling μ should halve Darcy flux"
        );
    }

    // ------------------------------------------------------------------
    // Forchheimer correction
    // ------------------------------------------------------------------

    #[test]
    fn test_forchheimer_degenerate_matches_darcy() {
        // With β = 0 Forchheimer must equal Darcy.
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let grad_p = -1000.0;
        let qd = darcy_flow(&m, grad_p);
        let qf = forchheimer_correction(&m, grad_p, 0.0);
        assert!(
            (qf - qd).abs() < 1e-20,
            "β=0 Forchheimer should equal Darcy: qf={qf}, qd={qd}"
        );
    }

    #[test]
    fn test_forchheimer_reduces_velocity() {
        // Inertial resistance always reduces velocity relative to Darcy.
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let grad_p = -1e6;
        let qd = darcy_flow(&m, grad_p);
        let qf = forchheimer_correction(&m, grad_p, 1e8);
        assert!(
            qf.abs() <= qd.abs() + 1e-30,
            "Forchheimer velocity ({qf}) should not exceed Darcy ({qd})"
        );
    }

    #[test]
    fn test_forchheimer_sign_consistency() {
        let m = PorousMedia::new(0.4, 1e-12, 1.5, 1e-3, 1000.0);
        let grad_p = -5000.0;
        let qf = forchheimer_correction(&m, grad_p, 1e6);
        assert!(
            qf > 0.0,
            "Forchheimer flux should be positive for negative grad_p, got {qf}"
        );
    }

    // ------------------------------------------------------------------
    // Kozeny-Carman
    // ------------------------------------------------------------------

    #[test]
    fn test_kozeny_carman_positive() {
        let k = kozeny_carman(0.4, 1e-4);
        assert!(
            k > 0.0,
            "Kozeny-Carman permeability should be positive, got {k}"
        );
    }

    #[test]
    fn test_kozeny_carman_formula() {
        let phi = 0.35;
        let dp = 2e-4;
        let k = kozeny_carman(phi, dp);
        let expected = phi.powi(3) * dp.powi(2) / (180.0 * (1.0 - phi).powi(2));
        assert!(
            (k - expected).abs() < EPS * expected,
            "K-C formula mismatch: k={k}, expected={expected}"
        );
    }

    #[test]
    fn test_kozeny_carman_larger_particle_larger_perm() {
        let k1 = kozeny_carman(0.4, 1e-4);
        let k2 = kozeny_carman(0.4, 2e-4);
        assert!(k2 > k1, "Larger particles should give larger permeability");
    }

    #[test]
    fn test_kozeny_carman_higher_porosity_larger_perm() {
        let k1 = kozeny_carman(0.3, 1e-4);
        let k2 = kozeny_carman(0.5, 1e-4);
        assert!(k2 > k1, "Higher porosity should give larger permeability");
    }

    #[test]
    fn test_kozeny_carman_quadratic_dp() {
        // Permeability scales with dp².
        let k1 = kozeny_carman(0.4, 1e-4);
        let k2 = kozeny_carman(0.4, 2e-4);
        let ratio = k2 / k1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "K-C should scale as d_p², ratio={ratio}"
        );
    }

    // ------------------------------------------------------------------
    // Effective thermal conductivity
    // ------------------------------------------------------------------

    #[test]
    fn test_k_eff_parallel_series_ordering() {
        let (kp, ks) = effective_thermal_conductivity(1.0, 0.025, 0.4);
        assert!(kp >= ks, "Parallel k ({kp}) should be >= series k ({ks})");
    }

    #[test]
    fn test_k_eff_pure_solid() {
        // φ = 0 → both bounds equal k_solid.
        let (kp, ks) = effective_thermal_conductivity(1.0, 0.025, 0.0);
        assert!((kp - 1.0).abs() < EPS);
        assert!((ks - 1.0).abs() < EPS);
    }

    #[test]
    fn test_k_eff_pure_fluid() {
        // φ = 1 → both bounds equal k_fluid.
        let (kp, ks) = effective_thermal_conductivity(1.0, 0.025, 1.0);
        assert!((kp - 0.025).abs() < EPS);
        assert!((ks - 0.025).abs() < EPS);
    }

    #[test]
    fn test_k_eff_same_phases() {
        // When k_solid = k_fluid, both models return that value.
        let (kp, ks) = effective_thermal_conductivity(0.5, 0.5, 0.4);
        assert!((kp - 0.5).abs() < EPS);
        assert!((ks - 0.5).abs() < EPS);
    }

    #[test]
    fn test_k_eff_parallel_formula() {
        let ks_val = 50.0;
        let kf_val = 0.1;
        let phi = 0.3;
        let (kp, _) = effective_thermal_conductivity(ks_val, kf_val, phi);
        let expected = (1.0 - phi) * ks_val + phi * kf_val;
        assert!((kp - expected).abs() < EPS * expected);
    }

    // ------------------------------------------------------------------
    // Gibson-Ashby stiffness
    // ------------------------------------------------------------------

    #[test]
    fn test_gibson_ashby_open_cell() {
        // Open-cell: C=1, n=2.  E*/E_s = (ρ*/ρ_s)².
        let rho_rel = 0.1;
        let rel_stiff = gibson_ashby_stiffness(rho_rel, 1.0, 2.0);
        assert!(
            (rel_stiff - 0.01).abs() < EPS,
            "GA stiffness should be {:.3}, got {rel_stiff}",
            0.01_f64
        );
    }

    #[test]
    fn test_gibson_ashby_scales_with_density() {
        let s1 = gibson_ashby_stiffness(0.1, 1.0, 2.0);
        let s2 = gibson_ashby_stiffness(0.2, 1.0, 2.0);
        assert!(
            s2 > s1,
            "Higher relative density should give higher stiffness"
        );
    }

    #[test]
    fn test_gibson_ashby_exponent_effect() {
        let s1 = gibson_ashby_stiffness(0.15, 1.0, 1.5);
        let s2 = gibson_ashby_stiffness(0.15, 1.0, 2.0);
        // For ρ* < 1 (foam), higher exponent gives lower relative stiffness.
        assert!(
            s2 < s1,
            "For ρ*/ρ_s < 1, higher exponent should give lower relative stiffness"
        );
    }

    #[test]
    fn test_gibson_ashby_unit_density() {
        // At ρ*/ρ_s = 1.0, E*/E_s = C (regardless of exponent).
        let c = 1.3;
        let rel_stiff = gibson_ashby_stiffness(1.0, c, 2.5);
        assert!(
            (rel_stiff - c).abs() < EPS,
            "At unit relative density, stiffness should equal C={c}"
        );
    }

    // ------------------------------------------------------------------
    // Plateau border geometry
    // ------------------------------------------------------------------

    #[test]
    fn test_plateau_positive() {
        let a = open_cell_plateau(0.01);
        assert!(a > 0.0, "Plateau border area should be positive, got {a}");
    }

    #[test]
    fn test_plateau_scales_quadratically() {
        let a1 = open_cell_plateau(0.01);
        let a2 = open_cell_plateau(0.02);
        let ratio = a2 / a1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "Plateau area should scale as r², ratio={ratio}"
        );
    }

    #[test]
    fn test_plateau_formula() {
        let r = 0.005;
        let a = open_cell_plateau(r);
        let expected = (3.0_f64.sqrt() - std::f64::consts::PI / 2.0) * r * r;
        assert!((a - expected).abs() < EPS * expected.abs().max(1e-20));
    }

    #[test]
    fn test_plateau_constant_positive() {
        // The geometric constant (√3 - π/2) should be positive.
        let constant = 3.0_f64.sqrt() - std::f64::consts::PI / 2.0;
        assert!(
            constant > 0.0,
            "Plateau geometric constant should be positive, got {constant}"
        );
    }
}
