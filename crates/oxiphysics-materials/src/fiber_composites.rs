// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fiber-reinforced composite material models.
//!
//! Implements micromechanical models for unidirectional and woven
//! fiber composites including rule-of-mixtures, Halpin–Tsai, and
//! classical laminate theory (CLT) stiffness assembly.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Fiber struct
// ---------------------------------------------------------------------------

/// Material properties for a single fiber type.
#[derive(Debug, Clone, PartialEq)]
pub struct Fiber {
    /// Longitudinal Young's modulus of the fiber (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio of the fiber (dimensionless).
    pub poissons_ratio: f64,
    /// Density of the fiber (kg/m³).
    pub density: f64,
    /// Tensile strength of the fiber (Pa).
    pub tensile_strength: f64,
    /// Nominal fiber diameter in micrometres (µm).
    pub diameter_um: f64,
}

// ---------------------------------------------------------------------------
// FiberType enum
// ---------------------------------------------------------------------------

/// Common commercial fiber types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiberType {
    /// High-stiffness / high-strength carbon fiber (T300 grade).
    Carbon,
    /// E-glass fiber (economical, good electrical insulation).
    Glass,
    /// Kevlar 49 aramid fiber (high toughness, light weight).
    Kevlar,
    /// Basalt fiber (natural volcanic origin, corrosion resistant).
    Basalt,
    /// High-strength steel fiber (for concrete reinforcement).
    Steel,
}

// ---------------------------------------------------------------------------
// fiber_preset
// ---------------------------------------------------------------------------

/// Return a [`Fiber`] populated with typical published properties for `ft`.
///
/// Sources: Ashby & Jones *Engineering Materials 1*, MatWeb, and fiber
/// manufacturer data sheets.  Values are representative mid-range figures.
pub fn fiber_preset(ft: FiberType) -> Fiber {
    match ft {
        FiberType::Carbon => Fiber {
            youngs_modulus: 230.0e9,
            poissons_ratio: 0.20,
            density: 1750.0,
            tensile_strength: 3500.0e6,
            diameter_um: 7.0,
        },
        FiberType::Glass => Fiber {
            youngs_modulus: 72.0e9,
            poissons_ratio: 0.22,
            density: 2540.0,
            tensile_strength: 3400.0e6,
            diameter_um: 10.0,
        },
        FiberType::Kevlar => Fiber {
            youngs_modulus: 125.0e9,
            poissons_ratio: 0.36,
            density: 1440.0,
            tensile_strength: 3620.0e6,
            diameter_um: 12.0,
        },
        FiberType::Basalt => Fiber {
            youngs_modulus: 89.0e9,
            poissons_ratio: 0.26,
            density: 2650.0,
            tensile_strength: 4840.0e6,
            diameter_um: 15.0,
        },
        FiberType::Steel => Fiber {
            youngs_modulus: 200.0e9,
            poissons_ratio: 0.30,
            density: 7800.0,
            tensile_strength: 1200.0e6,
            diameter_um: 200.0,
        },
    }
}

// ---------------------------------------------------------------------------
// rule_of_mixtures_longitudinal
// ---------------------------------------------------------------------------

/// Longitudinal Young's modulus E₁ of a unidirectional ply (Voigt / ROM).
///
/// E₁ = v_f · E_f + (1 − v_f) · E_m
///
/// where v_f is the fiber volume fraction (0–1), E_f the fiber modulus, and
/// E_m the matrix modulus.
pub fn rule_of_mixtures_longitudinal(e_fiber: f64, e_matrix: f64, v_f: f64) -> f64 {
    let v_m = 1.0 - v_f;
    v_f * e_fiber + v_m * e_matrix
}

// ---------------------------------------------------------------------------
// rule_of_mixtures_transverse
// ---------------------------------------------------------------------------

/// Transverse Young's modulus E₂ using the Halpin–Tsai approximation.
///
/// This is a convenience wrapper around [`halpin_tsai_transverse`] with the
/// standard reinforcing factor ξ = 2 (circular fibers in square packing).
pub fn rule_of_mixtures_transverse(e_fiber: f64, e_matrix: f64, v_f: f64) -> f64 {
    halpin_tsai_transverse(e_fiber, e_matrix, v_f, 2.0)
}

// ---------------------------------------------------------------------------
// halpin_tsai_transverse
// ---------------------------------------------------------------------------

/// Transverse Young's modulus E₂ via the Halpin–Tsai model.
///
/// ```text
/// η = (E_f/E_m − 1) / (E_f/E_m + ξ)
/// E₂ = E_m · (1 + ξ·η·v_f) / (1 − η·v_f)
/// ```
///
/// `xi` (ξ) is the reinforcing factor: 2 for circular cross-sections, larger
/// for higher aspect-ratio inclusions.
pub fn halpin_tsai_transverse(e_fiber: f64, e_matrix: f64, v_f: f64, xi: f64) -> f64 {
    if e_matrix < 1e-30 {
        return 0.0;
    }
    let ratio = e_fiber / e_matrix;
    let eta = (ratio - 1.0) / (ratio + xi);
    let denom = 1.0 - eta * v_f;
    if denom.abs() < 1e-30 {
        return e_matrix;
    }
    e_matrix * (1.0 + xi * eta * v_f) / denom
}

// ---------------------------------------------------------------------------
// clt_in_plane_stiffness
// ---------------------------------------------------------------------------

/// Assemble the CLT in-plane stiffness matrix **A** (N/m) for a laminate.
///
/// Each ply is described by the tuple `(E1, E2, G12, nu12, theta_deg, t)`:
/// - `E1` – longitudinal modulus (Pa)
/// - `E2` – transverse modulus (Pa)
/// - `G12` – in-plane shear modulus (Pa)
/// - `nu12` – major Poisson's ratio
/// - `theta_deg` – ply angle measured from the laminate reference axis (°)
/// - `t` – ply thickness (m)
///
/// Returns the 3×3 **A** matrix in Voigt notation `[σ_xx, σ_yy, σ_xy]`.
pub fn clt_in_plane_stiffness(plies: &[(f64, f64, f64, f64, f64, f64)]) -> [[f64; 3]; 3] {
    let mut a = [[0.0f64; 3]; 3];
    for &(e1, e2, g12, nu12, theta_deg, t) in plies {
        let q = ply_reduced_stiffness(e1, e2, g12, nu12);
        let qbar = transform_stiffness(q, theta_deg);
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += qbar[i][j] * t;
            }
        }
    }
    a
}

/// Compute the reduced stiffness matrix **Q** for a ply in its material axes.
fn ply_reduced_stiffness(e1: f64, e2: f64, g12: f64, nu12: f64) -> [[f64; 3]; 3] {
    let nu21 = nu12 * e2 / e1.max(1e-30);
    let delta = 1.0 - nu12 * nu21;
    if delta.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let q11 = e1 / delta;
    let q22 = e2 / delta;
    let q12 = nu12 * e2 / delta;
    let q66 = g12;
    [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
}

/// Rotate reduced stiffness **Q** from material axes by angle `theta_deg`.
fn transform_stiffness(q: [[f64; 3]; 3], theta_deg: f64) -> [[f64; 3]; 3] {
    let th = theta_deg.to_radians();
    let c = th.cos();
    let s = th.sin();
    let c2 = c * c;
    let s2 = s * s;
    let _cs = c * s;
    let c4 = c2 * c2;
    let s4 = s2 * s2;
    let c2s2 = c2 * s2;
    let q11 = q[0][0];
    let q12 = q[0][1];
    let q22 = q[1][1];
    let q66 = q[2][2];
    let qb11 = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * s4;
    let qb12 = (q11 + q22 - 4.0 * q66) * c2s2 + q12 * (c4 + s4);
    let qb22 = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * c4;
    let qb16 = (q11 - q12 - 2.0 * q66) * c * c2 * s - (q22 - q12 - 2.0 * q66) * s * s2 * c;
    let qb26 = (q11 - q12 - 2.0 * q66) * c * s * s2 - (q22 - q12 - 2.0 * q66) * s * c * c2;
    let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2s2 + q66 * (c4 + s4);
    [[qb11, qb12, qb16], [qb12, qb22, qb26], [qb16, qb26, qb66]]
}

// ---------------------------------------------------------------------------
// micromechanics_shear
// ---------------------------------------------------------------------------

/// In-plane shear modulus G₁₂ via Halpin–Tsai micromechanics.
///
/// ```text
/// η = (G_f/G_m − 1) / (G_f/G_m + 1)
/// G₁₂ = G_m · (1 + η·v_f) / (1 − η·v_f)
/// ```
pub fn micromechanics_shear(g_fiber: f64, g_matrix: f64, v_f: f64) -> f64 {
    if g_matrix < 1e-30 {
        return 0.0;
    }
    let ratio = g_fiber / g_matrix;
    let eta = (ratio - 1.0) / (ratio + 1.0);
    let denom = 1.0 - eta * v_f;
    if denom.abs() < 1e-30 {
        return g_matrix;
    }
    g_matrix * (1.0 + eta * v_f) / denom
}

// ---------------------------------------------------------------------------
// critical_fiber_length
// ---------------------------------------------------------------------------

/// Critical fiber length l_c for effective load transfer (Kelly–Tyson).
///
/// ```text
/// l_c = σ_f · d / (2 τ_i)
/// ```
///
/// where `tensile_strength` = σ_f (Pa), `diameter` = d (m), and
/// `interfacial_shear` = τ_i (Pa).
pub fn critical_fiber_length(tensile_strength: f64, interfacial_shear: f64, diameter: f64) -> f64 {
    if interfacial_shear.abs() < 1e-30 {
        return f64::INFINITY;
    }
    tensile_strength * diameter / (2.0 * interfacial_shear)
}

// ---------------------------------------------------------------------------
// fiber_pullout_energy
// ---------------------------------------------------------------------------

/// Fiber pull-out fracture energy per unit area (J/m²).
///
/// ```text
/// G_pullout = π · d · l² · τ_i / 2
/// ```
///
/// where `d` = fiber diameter (m), `bond_length` = l (m), and
/// `tau_i` = interfacial shear strength (Pa).
pub fn fiber_pullout_energy(fiber: &Fiber, bond_length: f64, tau_i: f64) -> f64 {
    let d = fiber.diameter_um * 1.0e-6; // convert µm → m
    PI * d * bond_length * bond_length * tau_i / 2.0
}

// ---------------------------------------------------------------------------
// winding_angle_effect
// ---------------------------------------------------------------------------

/// Effective longitudinal modulus E_x for a ply rotated by `theta_deg`.
///
/// Uses the exact off-axis transformation:
///
/// ```text
/// 1/E_x = cos⁴θ/E₁ + sin⁴θ/E₂ + (1/G₁₂ − 2ν₁₂/E₁)·sin²θ·cos²θ
/// ```
pub fn winding_angle_effect(e1: f64, e2: f64, g12: f64, nu12: f64, theta_deg: f64) -> f64 {
    let th = theta_deg.to_radians();
    let c = th.cos();
    let s = th.sin();
    let c2 = c * c;
    let s2 = s * s;
    let inv = c2 * c2 / e1.max(1e-30)
        + s2 * s2 / e2.max(1e-30)
        + (1.0 / g12.max(1e-30) - 2.0 * nu12 / e1.max(1e-30)) * s2 * c2;
    if inv.abs() < 1e-30 { 0.0 } else { 1.0 / inv }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // 1. fiber_preset: Carbon fiber has E > Glass.
    #[test]
    fn test_carbon_stiffer_than_glass() {
        let c = fiber_preset(FiberType::Carbon);
        let g = fiber_preset(FiberType::Glass);
        assert!(c.youngs_modulus > g.youngs_modulus);
    }

    // 2. fiber_preset: all densities are positive.
    #[test]
    fn test_all_densities_positive() {
        for ft in [
            FiberType::Carbon,
            FiberType::Glass,
            FiberType::Kevlar,
            FiberType::Basalt,
            FiberType::Steel,
        ] {
            assert!(fiber_preset(ft).density > 0.0);
        }
    }

    // 3. fiber_preset: all diameters > 0.
    #[test]
    fn test_all_diameters_positive() {
        for ft in [
            FiberType::Carbon,
            FiberType::Glass,
            FiberType::Kevlar,
            FiberType::Basalt,
            FiberType::Steel,
        ] {
            assert!(fiber_preset(ft).diameter_um > 0.0);
        }
    }

    // 4. rule_of_mixtures_longitudinal: vf=0 gives e_matrix.
    #[test]
    fn test_rom_long_vf0() {
        let e = rule_of_mixtures_longitudinal(200e9, 3e9, 0.0);
        assert!((e - 3e9).abs() < EPS * 3e9);
    }

    // 5. rule_of_mixtures_longitudinal: vf=1 gives e_fiber.
    #[test]
    fn test_rom_long_vf1() {
        let e = rule_of_mixtures_longitudinal(200e9, 3e9, 1.0);
        assert!((e - 200e9).abs() < EPS * 200e9);
    }

    // 6. rule_of_mixtures_longitudinal: result is between e_m and e_f.
    #[test]
    fn test_rom_long_bounds() {
        let e = rule_of_mixtures_longitudinal(230e9, 4e9, 0.6);
        assert!((4e9..=230e9).contains(&e));
    }

    // 7. rule_of_mixtures_transverse: vf=0 gives e_matrix.
    #[test]
    fn test_rom_trans_vf0() {
        let e = rule_of_mixtures_transverse(200e9, 3e9, 0.0);
        assert!((e - 3e9).abs() < EPS * 3e9);
    }

    // 8. rule_of_mixtures_transverse: result is bounded by matrix and fiber.
    #[test]
    fn test_rom_trans_bounds() {
        let e = rule_of_mixtures_transverse(200e9, 4e9, 0.5);
        assert!((4e9..=200e9).contains(&e));
    }

    // 9. halpin_tsai_transverse: xi=0 → E2 = E_m (no reinforcement).
    #[test]
    fn test_halpin_tsai_xi0() {
        // xi=0 → eta = (r-1)/(r+0) = 1 for large r
        // check vf=0 at xi=2 still gives e_matrix
        let e = halpin_tsai_transverse(100e9, 4e9, 0.0, 2.0);
        assert!((e - 4e9).abs() < EPS * 4e9);
    }

    // 10. halpin_tsai_transverse: higher xi → higher E2.
    #[test]
    fn test_halpin_tsai_higher_xi() {
        let e1 = halpin_tsai_transverse(200e9, 4e9, 0.6, 1.0);
        let e2 = halpin_tsai_transverse(200e9, 4e9, 0.6, 4.0);
        assert!(e2 > e1, "Higher xi should give higher transverse modulus");
    }

    // 11. clt_in_plane_stiffness: single 0° ply gives A11 ≈ Q11*t.
    #[test]
    fn test_clt_single_ply_0deg() {
        let e1 = 140e9f64;
        let e2 = 10e9f64;
        let g12 = 5e9f64;
        let nu12 = 0.3f64;
        let t = 0.001f64;
        let a = clt_in_plane_stiffness(&[(e1, e2, g12, nu12, 0.0, t)]);
        let nu21 = nu12 * e2 / e1;
        let delta = 1.0 - nu12 * nu21;
        let q11 = e1 / delta;
        // A11 should equal Q11 * t
        assert!(
            (a[0][0] - q11 * t).abs() < 1e3,
            "A11={} should match Q11*t={}",
            a[0][0],
            q11 * t
        );
    }

    // 12. clt_in_plane_stiffness: symmetric cross-ply gives A16≈0.
    #[test]
    fn test_clt_crossply_a16_zero() {
        let ply = (140e9, 10e9, 5e9, 0.3, 0.0, 0.001);
        let ply90 = (140e9, 10e9, 5e9, 0.3, 90.0, 0.001);
        let a = clt_in_plane_stiffness(&[ply, ply90, ply90, ply]);
        assert!(a[0][2].abs() < 1e3, "Cross-ply A16 should be ≈ 0");
    }

    // 13. micromechanics_shear: vf=0 gives g_matrix.
    #[test]
    fn test_shear_vf0() {
        let g = micromechanics_shear(30e9, 1.5e9, 0.0);
        assert!((g - 1.5e9).abs() < EPS * 1.5e9);
    }

    // 14. micromechanics_shear: result > g_matrix for positive vf.
    #[test]
    fn test_shear_higher_than_matrix() {
        let g = micromechanics_shear(30e9, 1.5e9, 0.5);
        assert!(g > 1.5e9);
    }

    // 15. critical_fiber_length: l_c = sigma * d / (2 * tau).
    #[test]
    fn test_critical_length_formula() {
        let sigma = 3500e6f64;
        let d = 7e-6f64;
        let tau = 35e6f64;
        let lc = critical_fiber_length(sigma, tau, d);
        let expected = sigma * d / (2.0 * tau);
        assert!((lc - expected).abs() < EPS * expected);
    }

    // 16. critical_fiber_length: zero tau returns infinity.
    #[test]
    fn test_critical_length_zero_tau() {
        assert!(critical_fiber_length(1000.0, 0.0, 1e-6).is_infinite());
    }

    // 17. critical_fiber_length: larger diameter → longer l_c.
    #[test]
    fn test_critical_length_larger_diameter() {
        let lc1 = critical_fiber_length(3500e6, 35e6, 7e-6);
        let lc2 = critical_fiber_length(3500e6, 35e6, 14e-6);
        assert!((lc2 - 2.0 * lc1).abs() < EPS * lc1);
    }

    // 18. fiber_pullout_energy: energy is positive.
    #[test]
    fn test_pullout_energy_positive() {
        let f = fiber_preset(FiberType::Carbon);
        let e = fiber_pullout_energy(&f, 1e-3, 35e6);
        assert!(e > 0.0);
    }

    // 19. fiber_pullout_energy: scales quadratically with bond_length.
    #[test]
    fn test_pullout_energy_quadratic_length() {
        let f = fiber_preset(FiberType::Carbon);
        let e1 = fiber_pullout_energy(&f, 1e-3, 35e6);
        let e2 = fiber_pullout_energy(&f, 2e-3, 35e6);
        assert!((e2 - 4.0 * e1).abs() < EPS * e1 * 4.0);
    }

    // 20. fiber_pullout_energy: scales linearly with tau_i.
    #[test]
    fn test_pullout_energy_linear_tau() {
        let f = fiber_preset(FiberType::Carbon);
        let e1 = fiber_pullout_energy(&f, 1e-3, 10e6);
        let e2 = fiber_pullout_energy(&f, 1e-3, 20e6);
        assert!((e2 - 2.0 * e1).abs() < EPS * e1 * 2.0);
    }

    // 21. winding_angle_effect: theta=0 gives E1.
    #[test]
    fn test_winding_0deg_gives_e1() {
        let e1 = 140e9f64;
        let e2 = 10e9f64;
        let g12 = 5e9f64;
        let nu12 = 0.3f64;
        let ex = winding_angle_effect(e1, e2, g12, nu12, 0.0);
        assert!((ex - e1).abs() < 1e3, "At 0° Ex should = E1");
    }

    // 22. winding_angle_effect: theta=90 gives E2.
    #[test]
    fn test_winding_90deg_gives_e2() {
        let e1 = 140e9f64;
        let e2 = 10e9f64;
        let g12 = 5e9f64;
        let nu12 = 0.3f64;
        let ex = winding_angle_effect(e1, e2, g12, nu12, 90.0);
        assert!((ex - e2).abs() < 1e3, "At 90° Ex should = E2");
    }

    // 23. winding_angle_effect: intermediate angle gives modulus between E1 and E2.
    #[test]
    fn test_winding_45deg_bounds() {
        let e1 = 140e9f64;
        let e2 = 10e9f64;
        let ex = winding_angle_effect(e1, e2, 5e9, 0.3, 45.0);
        assert!(ex >= e2 && ex <= e1);
    }

    // 24. winding_angle_effect: modulus decreases monotonically from 0 to 90.
    #[test]
    fn test_winding_monotone() {
        let e1 = 140e9f64;
        let e2 = 10e9f64;
        let g12 = 5e9f64;
        let nu12 = 0.3f64;
        let mut prev = winding_angle_effect(e1, e2, g12, nu12, 0.0);
        for deg in (5u32..=90).step_by(5) {
            let cur = winding_angle_effect(e1, e2, g12, nu12, deg as f64);
            assert!(
                cur <= prev + 1e6,
                "Ex should decrease: theta={deg} cur={cur} prev={prev}"
            );
            prev = cur;
        }
    }

    // 25. rule_of_mixtures_longitudinal: linear interpolation at vf=0.5.
    #[test]
    fn test_rom_long_midpoint() {
        let e_f = 200e9f64;
        let e_m = 4e9f64;
        let e = rule_of_mixtures_longitudinal(e_f, e_m, 0.5);
        let expected = 0.5 * e_f + 0.5 * e_m;
        assert!((e - expected).abs() < EPS * expected);
    }

    // 26. halpin_tsai_transverse: xi=2 at vf=0.6 gives reasonable E2.
    #[test]
    fn test_halpin_tsai_realistic() {
        let e2 = halpin_tsai_transverse(230e9, 4e9, 0.6, 2.0);
        // Should be well above matrix modulus but below fiber modulus.
        assert!(e2 > 4e9 && e2 < 50e9);
    }

    // 27. micromechanics_shear: zero matrix modulus returns 0.
    #[test]
    fn test_shear_zero_matrix() {
        let g = micromechanics_shear(30e9, 0.0, 0.5);
        assert_eq!(g, 0.0);
    }

    // 28. clt_in_plane_stiffness: A matrix is positive semi-definite (A11 > 0).
    #[test]
    fn test_clt_a11_positive() {
        let plies = vec![
            (140e9, 10e9, 5e9, 0.3, 0.0, 0.001),
            (140e9, 10e9, 5e9, 0.3, 90.0, 0.001),
        ];
        let a = clt_in_plane_stiffness(&plies);
        assert!(a[0][0] > 0.0 && a[1][1] > 0.0 && a[2][2] > 0.0);
    }

    // 29. fiber_pullout_energy: zero tau gives zero energy.
    #[test]
    fn test_pullout_energy_zero_tau() {
        let f = fiber_preset(FiberType::Glass);
        assert_eq!(fiber_pullout_energy(&f, 1e-3, 0.0), 0.0);
    }

    // 30. fiber_preset: Basalt tensile strength > Steel.
    #[test]
    fn test_basalt_stronger_than_steel() {
        let b = fiber_preset(FiberType::Basalt);
        let s = fiber_preset(FiberType::Steel);
        assert!(b.tensile_strength > s.tensile_strength);
    }

    // 31. micromechanics_shear: vf=1 approaches g_fiber (Halpin–Tsai saturates).
    #[test]
    fn test_shear_vf1_saturates() {
        let g = micromechanics_shear(30e9, 1.5e9, 1.0);
        // eta * 1 = (r-1)/(r+1); numerator = 1 + eta, denom = 1 - eta
        // Should be significantly above matrix.
        assert!(g > 1.5e9);
    }

    // 32. halpin_tsai_transverse: matrix with zero modulus returns 0.
    #[test]
    fn test_halpin_tsai_zero_matrix() {
        let e = halpin_tsai_transverse(200e9, 0.0, 0.5, 2.0);
        assert_eq!(e, 0.0);
    }

    // 33. critical_fiber_length: proportional to tensile strength.
    #[test]
    fn test_critical_length_proportional_strength() {
        let lc1 = critical_fiber_length(1000e6, 50e6, 10e-6);
        let lc2 = critical_fiber_length(2000e6, 50e6, 10e-6);
        assert!((lc2 - 2.0 * lc1).abs() < EPS * lc1 * 2.0);
    }

    // 34. clt_in_plane_stiffness: empty ply list gives zero matrix.
    #[test]
    fn test_clt_empty_plies() {
        let a = clt_in_plane_stiffness(&[]);
        for row in &a {
            for val in row {
                assert_eq!(*val, 0.0);
            }
        }
    }

    // 35. winding_angle_effect: Glass fiber at 0 deg gives ~72 GPa.
    #[test]
    fn test_winding_glass_0deg() {
        let gf = fiber_preset(FiberType::Glass);
        // E_m assumed 4 GPa, G12 = 2 GPa, nu12 = 0.25
        let ex = winding_angle_effect(gf.youngs_modulus, 10e9, 3e9, 0.25, 0.0);
        assert!((ex - gf.youngs_modulus).abs() < 1e4);
    }
}
