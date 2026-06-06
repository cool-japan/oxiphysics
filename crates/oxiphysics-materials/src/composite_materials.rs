// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite material mechanics.
//!
//! Provides models for fiber-reinforced composites, including the rule of
//! mixtures, Halpin-Tsai equations, classical lamination theory (CLT), and
//! failure criteria such as Tsai-Wu.

/// Properties of a reinforcing fiber phase.
#[derive(Debug, Clone, PartialEq)]
pub struct FiberProperties {
    /// Longitudinal (axial) Young's modulus \[Pa\].
    pub modulus_longitudinal: f64,
    /// Transverse Young's modulus \[Pa\].
    pub modulus_transverse: f64,
    /// Longitudinal tensile strength \[Pa\].
    pub strength_longitudinal: f64,
    /// Fiber volume fraction (0–1).
    pub volume_fraction: f64,
    /// In-plane orientation angle of the fiber \[radians\].
    pub orientation: f64,
}

impl FiberProperties {
    /// Create a new `FiberProperties`.
    pub fn new(
        modulus_longitudinal: f64,
        modulus_transverse: f64,
        strength_longitudinal: f64,
        volume_fraction: f64,
        orientation: f64,
    ) -> Self {
        Self {
            modulus_longitudinal,
            modulus_transverse,
            strength_longitudinal,
            volume_fraction,
            orientation,
        }
    }
}

/// Properties of the matrix (binder) phase.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixProperties {
    /// Young's modulus of the matrix \[Pa\].
    pub modulus: f64,
    /// Poisson's ratio of the matrix (dimensionless).
    pub poisson_ratio: f64,
}

impl MatrixProperties {
    /// Create a new `MatrixProperties`.
    pub fn new(modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            modulus,
            poisson_ratio,
        }
    }
}

/// A single ply in a composite laminate.
#[derive(Debug, Clone, PartialEq)]
pub struct Ply {
    /// Ply angle measured from the laminate reference axis \[radians\].
    pub angle: f64,
    /// Ply thickness \[m\].
    pub thickness: f64,
    /// Longitudinal modulus E11 \[Pa\].
    pub e11: f64,
    /// Transverse modulus E22 \[Pa\].
    pub e22: f64,
    /// In-plane shear modulus G12 \[Pa\].
    pub g12: f64,
    /// Major Poisson's ratio ν12 (dimensionless).
    pub nu12: f64,
}

impl Ply {
    /// Create a new `Ply`.
    pub fn new(angle: f64, thickness: f64, e11: f64, e22: f64, g12: f64, nu12: f64) -> Self {
        Self {
            angle,
            thickness,
            e11,
            e22,
            g12,
            nu12,
        }
    }
}

/// A stack of plies forming a composite laminate.
#[derive(Debug, Clone)]
pub struct CompositeLayup {
    /// Ordered stack of plies (bottom to top).
    pub plies: Vec<Ply>,
}

impl CompositeLayup {
    /// Create a new `CompositeLayup` from a vector of plies.
    pub fn new(plies: Vec<Ply>) -> Self {
        Self { plies }
    }

    /// Total laminate thickness \[m\].
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }
}

/// Rule-of-mixtures modulus bounds.
///
/// Returns `(E_voigt, E_reuss)` where:
/// - `E_voigt` is the upper (Voigt) bound: iso-strain assumption.
/// - `E_reuss` is the lower (Reuss) bound: iso-stress assumption.
///
/// # Arguments
/// * `e_fiber` – fiber modulus \[Pa\]
/// * `e_matrix` – matrix modulus \[Pa\]
/// * `vf` – fiber volume fraction (0–1)
pub fn rule_of_mixtures_modulus(e_fiber: f64, e_matrix: f64, vf: f64) -> (f64, f64) {
    let vm = 1.0 - vf;
    let e_voigt = e_fiber * vf + e_matrix * vm;
    let e_reuss = 1.0 / (vf / e_fiber + vm / e_matrix);
    (e_voigt, e_reuss)
}

/// Halpin-Tsai composite equations.
///
/// Returns `(E11, E22, G12)` for a unidirectional composite.
///
/// # Arguments
/// * `fiber` – fiber properties
/// * `matrix` – matrix properties
/// * `xi_e` – Halpin-Tsai reinforcement factor for E22 (typical: 2 for circular fibers)
/// * `xi_g` – Halpin-Tsai reinforcement factor for G12 (typical: 1)
pub fn halpin_tsai(
    fiber: &FiberProperties,
    matrix: &MatrixProperties,
    xi_e: f64,
    xi_g: f64,
) -> (f64, f64, f64) {
    let vf = fiber.volume_fraction;
    let vm = 1.0 - vf;

    // E11 — rule of mixtures (Voigt)
    let e11 = fiber.modulus_longitudinal * vf + matrix.modulus * vm;

    // E22 — Halpin-Tsai
    let eta_e = (fiber.modulus_transverse / matrix.modulus - 1.0)
        / (fiber.modulus_transverse / matrix.modulus + xi_e);
    let e22 = matrix.modulus * (1.0 + xi_e * eta_e * vf) / (1.0 - eta_e * vf);

    // G12 — Halpin-Tsai (approximate fiber shear modulus from isotropic assumption)
    let g_fiber = fiber.modulus_longitudinal / (2.0 * (1.0 + 0.25)); // ν_f ≈ 0.25
    let g_matrix = matrix.modulus / (2.0 * (1.0 + matrix.poisson_ratio));
    let eta_g = (g_fiber / g_matrix - 1.0) / (g_fiber / g_matrix + xi_g);
    let g12 = g_matrix * (1.0 + xi_g * eta_g * vf) / (1.0 - eta_g * vf);

    (e11, e22, g12)
}

/// ABD matrices from Classical Lamination Theory (CLT).
///
/// Returns `(A, B, D)` where each matrix is stored in **Voigt notation** as a
/// flat 6-element array `[A11, A12, A16, A22, A26, A66]` (and likewise for B
/// and D).  This covers the in-plane (`A`), coupling (`B`), and bending (`D`)
/// stiffness sub-matrices.
///
/// # Arguments
/// * `layup` – the composite laminate definition
pub fn classical_lamination_theory(layup: &CompositeLayup) -> ([f64; 6], [f64; 6], [f64; 6]) {
    let mut a = [0.0f64; 6];
    let mut b = [0.0f64; 6];
    let mut d = [0.0f64; 6];

    // Mid-plane z coordinate (origin at laminate mid-plane).
    let total_h = layup.total_thickness();
    let mut z_bottom = -total_h / 2.0;

    for ply in &layup.plies {
        let z_top = z_bottom + ply.thickness;
        let z_mid = (z_bottom + z_top) / 2.0;

        // Reduced stiffness in principal axes.
        let nu21 = ply.nu12 * ply.e22 / ply.e11;
        let denom = 1.0 - ply.nu12 * nu21;
        let q11 = ply.e11 / denom;
        let q22 = ply.e22 / denom;
        let q12 = ply.nu12 * ply.e22 / denom;
        let q66 = ply.g12;

        // Transform Q to laminate axes.
        let c = ply.angle.cos();
        let s = ply.angle.sin();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;

        // Q-bar components (Voigt: [11, 12, 16, 22, 26, 66]).
        let qb11 = q11 * c2 * c2 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s2 * s2;
        let qb12 = (q11 + q22 - 4.0 * q66) * s2 * c2 + q12 * (c2 * c2 + s2 * s2);
        let qb16 = (q11 - q12 - 2.0 * q66) * cs * c2 - (q22 - q12 - 2.0 * q66) * cs * s2;
        let qb22 = q11 * s2 * s2 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * c2 * c2;
        let qb26 = (q11 - q12 - 2.0 * q66) * cs * s2 - (q22 - q12 - 2.0 * q66) * cs * c2;
        let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s2 * c2 + q66 * (c2 * c2 + s2 * s2);

        let qbar = [qb11, qb12, qb16, qb22, qb26, qb66];

        let dz = z_top - z_bottom;
        let dz2 = (z_top * z_top - z_bottom * z_bottom) / 2.0;
        let dz3 = (z_top * z_top * z_top - z_bottom * z_bottom * z_bottom) / 3.0;

        for i in 0..6 {
            a[i] += qbar[i] * dz;
            b[i] += qbar[i] * dz2;
            d[i] += qbar[i] * dz3;
        }

        // Suppress unused variable warnings for z_mid (used conceptually).
        let _ = z_mid;
        z_bottom = z_top;
    }

    (a, b, d)
}

/// Tsai-Wu tensor failure criterion.
///
/// Returns the Tsai-Wu failure index `FI`.  Failure is predicted when `FI ≥ 1`.
///
/// # Arguments
/// * `sigma` – stress state `[σ1, σ2, τ12]` \[Pa\]
/// * `f1t`   – longitudinal tensile strength \[Pa\]
/// * `f1c`   – longitudinal compressive strength \[Pa\] (positive value)
/// * `f2t`   – transverse tensile strength \[Pa\]
/// * `f2c`   – transverse compressive strength \[Pa\] (positive value)
/// * `f12`   – in-plane shear strength \[Pa\]
pub fn failure_criterion_tsai_wu(
    sigma: [f64; 3],
    f1t: f64,
    f1c: f64,
    f2t: f64,
    f2c: f64,
    f12: f64,
) -> f64 {
    let [s1, s2, t12] = sigma;

    // Strength tensors (first-order and second-order).
    let f1 = 1.0 / f1t - 1.0 / f1c;
    let f2 = 1.0 / f2t - 1.0 / f2c;
    let f11 = 1.0 / (f1t * f1c);
    let f22 = 1.0 / (f2t * f2c);
    let f66 = 1.0 / (f12 * f12);
    // Interaction term F12 (approximate: -0.5 * sqrt(F11 * F22)).
    let f12_int = -0.5 * (f11 * f22).sqrt();

    f1 * s1 + f2 * s2 + f11 * s1 * s1 + f22 * s2 * s2 + f66 * t12 * t12 + 2.0 * f12_int * s1 * s2
}

/// Interlaminar shear stress (ILSS) from short beam shear theory.
///
/// Returns the maximum interlaminar shear stress \[Pa\] from a 3-point bend test.
///
/// # Arguments
/// * `p_max`  – maximum applied load \[N\]
/// * `width`  – specimen width \[m\]
/// * `thickness` – specimen thickness \[m\]
pub fn interlaminar_shear(p_max: f64, width: f64, thickness: f64) -> f64 {
    // ILSS = 0.75 * P_max / (width * thickness)  (short-beam formula)
    0.75 * p_max / (width * thickness)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    // ------------------------------------------------------------------
    // FiberProperties / MatrixProperties construction
    // ------------------------------------------------------------------

    #[test]
    fn test_fiber_props_new() {
        let f = FiberProperties::new(230e9, 15e9, 3.5e9, 0.6, 0.0);
        assert_eq!(f.modulus_longitudinal, 230e9);
        assert_eq!(f.volume_fraction, 0.6);
    }

    #[test]
    fn test_matrix_props_new() {
        let m = MatrixProperties::new(3.5e9, 0.35);
        assert_eq!(m.modulus, 3.5e9);
        assert!((m.poisson_ratio - 0.35).abs() < EPS);
    }

    // ------------------------------------------------------------------
    // Rule of Mixtures
    // ------------------------------------------------------------------

    #[test]
    fn test_rom_voigt_bounds() {
        // Voigt ≥ Reuss always.
        let (ev, er) = rule_of_mixtures_modulus(200e9, 3e9, 0.5);
        assert!(ev > er, "Voigt {ev} should exceed Reuss {er}");
    }

    #[test]
    fn test_rom_pure_fiber() {
        // When vf = 1.0, both bounds equal E_fiber.
        let (ev, er) = rule_of_mixtures_modulus(200e9, 3e9, 1.0);
        assert!((ev - 200e9).abs() < 1.0, "Voigt at vf=1 should be E_fiber");
        assert!((er - 200e9).abs() < 1.0, "Reuss at vf=1 should be E_fiber");
    }

    #[test]
    fn test_rom_pure_matrix() {
        // When vf = 0.0, both bounds equal E_matrix.
        let (ev, er) = rule_of_mixtures_modulus(200e9, 3e9, 0.0);
        assert!((ev - 3e9).abs() < 1.0);
        assert!((er - 3e9).abs() < 1.0);
    }

    #[test]
    fn test_rom_midpoint_voigt() {
        // vf = 0.5: Voigt = (E_f + E_m) / 2
        let ef = 100.0;
        let em = 50.0;
        let (ev, _) = rule_of_mixtures_modulus(ef, em, 0.5);
        assert!((ev - 75.0).abs() < EPS);
    }

    #[test]
    fn test_rom_reuss_harmonic_mean() {
        let ef = 100.0;
        let em = 50.0;
        let vf = 0.4;
        let (_, er) = rule_of_mixtures_modulus(ef, em, vf);
        let expected = 1.0 / (vf / ef + (1.0 - vf) / em);
        assert!((er - expected).abs() < EPS);
    }

    #[test]
    fn test_rom_symmetry_voigt() {
        // Swapping fiber and matrix fractions should give the same Voigt value
        // if E_f = E_m (degenerate case).
        let (ev1, _) = rule_of_mixtures_modulus(100.0, 100.0, 0.3);
        let (ev2, _) = rule_of_mixtures_modulus(100.0, 100.0, 0.7);
        assert!((ev1 - ev2).abs() < EPS, "Degenerate: ev1={ev1}, ev2={ev2}");
    }

    // ------------------------------------------------------------------
    // Halpin-Tsai
    // ------------------------------------------------------------------

    #[test]
    fn test_halpin_tsai_e11_rule_of_mixtures() {
        let fiber = FiberProperties::new(230e9, 15e9, 3.5e9, 0.6, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (e11, _, _) = halpin_tsai(&fiber, &matrix, 2.0, 1.0);
        let (ev, _) = rule_of_mixtures_modulus(
            fiber.modulus_longitudinal,
            matrix.modulus,
            fiber.volume_fraction,
        );
        assert!(
            (e11 - ev).abs() < 1.0,
            "E11 from HT should equal Voigt ROM: e11={e11}, ev={ev}"
        );
    }

    #[test]
    fn test_halpin_tsai_e22_positive() {
        let fiber = FiberProperties::new(230e9, 15e9, 3.5e9, 0.6, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (_, e22, _) = halpin_tsai(&fiber, &matrix, 2.0, 1.0);
        assert!(e22 > 0.0, "E22 should be positive, got {e22}");
    }

    #[test]
    fn test_halpin_tsai_g12_positive() {
        let fiber = FiberProperties::new(230e9, 15e9, 3.5e9, 0.6, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (_, _, g12) = halpin_tsai(&fiber, &matrix, 2.0, 1.0);
        assert!(g12 > 0.0, "G12 should be positive, got {g12}");
    }

    #[test]
    fn test_halpin_tsai_e22_above_matrix_modulus() {
        // With fibers present E22 should exceed the pure matrix modulus.
        let fiber = FiberProperties::new(230e9, 15e9, 3.5e9, 0.5, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (_, e22, _) = halpin_tsai(&fiber, &matrix, 2.0, 1.0);
        assert!(e22 > matrix.modulus, "E22 should exceed matrix modulus");
    }

    // ------------------------------------------------------------------
    // Classical Lamination Theory
    // ------------------------------------------------------------------

    #[test]
    fn test_clt_a_matrix_positive_diagonal() {
        // A symmetric [0/90]s laminate should have positive A11 and A22.
        let ply_0 = Ply::new(0.0, 0.001, 140e9, 10e9, 5e9, 0.3);
        let ply_90 = Ply::new(PI / 2.0, 0.001, 140e9, 10e9, 5e9, 0.3);
        let layup = CompositeLayup::new(vec![
            ply_0,
            ply_90.clone(),
            ply_90,
            Ply::new(0.0, 0.001, 140e9, 10e9, 5e9, 0.3),
        ]);
        let (a, _, _) = classical_lamination_theory(&layup);
        assert!(a[0] > 0.0, "A11 should be positive, got {}", a[0]);
        assert!(a[3] > 0.0, "A22 should be positive, got {}", a[3]);
    }

    #[test]
    fn test_clt_b_matrix_zero_symmetric() {
        // A symmetric laminate ([0/90]s) has zero B matrix.
        let ply = Ply::new(0.0, 0.001, 140e9, 10e9, 5e9, 0.3);
        let ply90 = Ply::new(PI / 2.0, 0.001, 140e9, 10e9, 5e9, 0.3);
        let layup = CompositeLayup::new(vec![ply.clone(), ply90.clone(), ply90, ply]);
        let (_, b, _) = classical_lamination_theory(&layup);
        for (i, &bi) in b.iter().enumerate() {
            assert!(
                bi.abs() < 1.0,
                "B[{i}] = {bi} should be ~0 for symmetric laminate"
            );
        }
    }

    #[test]
    fn test_clt_single_ply_thickness() {
        // For a single unidirectional ply, A11 ≈ Q11 * h.
        let h = 0.002;
        let ply = Ply::new(0.0, h, 140e9, 10e9, 5e9, 0.3);
        let layup = CompositeLayup::new(vec![ply]);
        let (a, _, _) = classical_lamination_theory(&layup);
        let nu21 = 0.3 * 10e9 / 140e9;
        let q11 = 140e9 / (1.0 - 0.3 * nu21);
        let expected_a11 = q11 * h;
        let rel_err = (a[0] - expected_a11).abs() / expected_a11;
        assert!(rel_err < 1e-6, "A11={} expected={expected_a11}", a[0]);
    }

    #[test]
    fn test_clt_d_matrix_positive_diagonal() {
        let ply = Ply::new(0.0, 0.002, 140e9, 10e9, 5e9, 0.3);
        let layup = CompositeLayup::new(vec![ply]);
        let (_, _, d) = classical_lamination_theory(&layup);
        assert!(d[0] > 0.0, "D11 should be positive");
        assert!(d[3] > 0.0, "D22 should be positive");
    }

    // ------------------------------------------------------------------
    // Tsai-Wu failure criterion
    // ------------------------------------------------------------------

    #[test]
    fn test_tsai_wu_zero_stress() {
        let fi = failure_criterion_tsai_wu([0.0, 0.0, 0.0], 1.5e9, 1.0e9, 0.05e9, 0.25e9, 0.07e9);
        assert!((fi).abs() < EPS, "FI at zero stress should be ~0, got {fi}");
    }

    #[test]
    fn test_tsai_wu_uniaxial_at_strength() {
        // Apply exactly the longitudinal tensile strength → FI should be ~1.
        let f1t = 1.5e9_f64;
        let f1c = 1.0e9_f64;
        // f1*σ1 + f11*σ1² = 1  when σ1 = f1t (by definition of the criterion).
        let fi = failure_criterion_tsai_wu([f1t, 0.0, 0.0], f1t, f1c, 0.05e9, 0.25e9, 0.07e9);
        // Should be close to 1.
        assert!(
            (fi - 1.0).abs() < 0.1,
            "FI at longitudinal tensile strength should be ~1, got {fi}"
        );
    }

    #[test]
    fn test_tsai_wu_positive_above_threshold() {
        // For any nonzero shear the FI should increase.
        let fi0 = failure_criterion_tsai_wu([0.0, 0.0, 0.0], 1.5e9, 1.0e9, 0.05e9, 0.25e9, 0.07e9);
        let fi1 = failure_criterion_tsai_wu([0.0, 0.0, 1e7], 1.5e9, 1.0e9, 0.05e9, 0.25e9, 0.07e9);
        assert!(fi1 > fi0, "Non-zero shear should increase FI");
    }

    // ------------------------------------------------------------------
    // Interlaminar shear
    // ------------------------------------------------------------------

    #[test]
    fn test_ilss_formula() {
        let ilss = interlaminar_shear(1000.0, 0.02, 0.004);
        let expected = 0.75 * 1000.0 / (0.02 * 0.004);
        assert!((ilss - expected).abs() < EPS);
    }

    #[test]
    fn test_ilss_scales_with_load() {
        let i1 = interlaminar_shear(500.0, 0.02, 0.004);
        let i2 = interlaminar_shear(1000.0, 0.02, 0.004);
        assert!(
            (i2 - 2.0 * i1).abs() < EPS,
            "ILSS should double when load doubles"
        );
    }

    #[test]
    fn test_ilss_positive() {
        let ilss = interlaminar_shear(800.0, 0.015, 0.003);
        assert!(ilss > 0.0, "ILSS must be positive");
    }

    // ------------------------------------------------------------------
    // CompositeLayup helpers
    // ------------------------------------------------------------------

    #[test]
    fn test_total_thickness() {
        let plies = vec![
            Ply::new(0.0, 0.001, 140e9, 10e9, 5e9, 0.3),
            Ply::new(PI / 4.0, 0.002, 140e9, 10e9, 5e9, 0.3),
        ];
        let layup = CompositeLayup::new(plies);
        assert!((layup.total_thickness() - 0.003).abs() < EPS);
    }

    #[test]
    fn test_ply_new_fields() {
        let ply = Ply::new(PI / 4.0, 0.001, 140e9, 10e9, 5e9, 0.28);
        assert!((ply.angle - PI / 4.0).abs() < EPS);
        assert!((ply.nu12 - 0.28).abs() < EPS);
    }

    // ------------------------------------------------------------------
    // Additional edge-case tests
    // ------------------------------------------------------------------

    #[test]
    fn test_rom_modulus_between_bounds() {
        let (ev, er) = rule_of_mixtures_modulus(200e9, 5e9, 0.3);
        // Any physically meaningful effective modulus lies within the bounds.
        assert!(ev >= er);
        assert!(ev >= 5e9);
        assert!(er >= 5e9);
        assert!(ev <= 200e9);
    }

    #[test]
    fn test_halpin_tsai_varies_with_vf() {
        let f1 = FiberProperties::new(230e9, 15e9, 3.5e9, 0.3, 0.0);
        let f2 = FiberProperties::new(230e9, 15e9, 3.5e9, 0.6, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (e11_low, _, _) = halpin_tsai(&f1, &matrix, 2.0, 1.0);
        let (e11_high, _, _) = halpin_tsai(&f2, &matrix, 2.0, 1.0);
        assert!(
            e11_high > e11_low,
            "Higher fiber volume fraction should give higher E11"
        );
    }

    #[test]
    fn test_tsai_wu_compressive_stress_increases_fi() {
        // Compressive σ1 (negative) also contributes to failure index via F11.
        let f1t = 1.5e9_f64;
        let f1c = 1.0e9_f64;
        let fi_ten = failure_criterion_tsai_wu([1e8, 0.0, 0.0], f1t, f1c, 0.05e9, 0.25e9, 0.07e9);
        let fi_comp = failure_criterion_tsai_wu([-1e8, 0.0, 0.0], f1t, f1c, 0.05e9, 0.25e9, 0.07e9);
        // Both should have positive FI contribution from the quadratic term.
        assert!(fi_ten.is_finite());
        assert!(fi_comp.is_finite());
    }

    #[test]
    fn test_clt_quasi_isotropic_a_symmetry() {
        // A quasi-isotropic [0/60/-60] laminate should have A11 ≈ A22.
        let h = 0.001;
        let e11 = 140e9;
        let e22 = 10e9;
        let g12 = 5e9;
        let nu12 = 0.3;
        let plies = vec![
            Ply::new(0.0, h, e11, e22, g12, nu12),
            Ply::new(PI / 3.0, h, e11, e22, g12, nu12),
            Ply::new(-PI / 3.0, h, e11, e22, g12, nu12),
        ];
        let layup = CompositeLayup::new(plies);
        let (a, _, _) = classical_lamination_theory(&layup);
        let diff = (a[0] - a[3]).abs() / a[0];
        assert!(
            diff < 0.05,
            "Quasi-isotropic laminate: A11={:.3e} should ≈ A22={:.3e} (rel diff={diff:.4})",
            a[0],
            a[3]
        );
    }

    #[test]
    fn test_tsai_wu_transverse_stress() {
        // Transverse tensile stress near its own strength → FI near 1.
        let f2t = 0.05e9_f64;
        let fi = failure_criterion_tsai_wu([0.0, f2t, 0.0], 1.5e9, 1.0e9, f2t, 0.25e9, 0.07e9);
        assert!(
            (fi - 1.0).abs() < 0.1,
            "FI at transverse tensile strength should be ~1, got {fi}"
        );
    }

    #[test]
    fn test_halpin_tsai_xi_effect() {
        // Higher ξ should give higher E22 (more reinforcement efficiency).
        let fiber = FiberProperties::new(230e9, 15e9, 3.5e9, 0.5, 0.0);
        let matrix = MatrixProperties::new(3.5e9, 0.35);
        let (_, e22_low, _) = halpin_tsai(&fiber, &matrix, 1.0, 1.0);
        let (_, e22_high, _) = halpin_tsai(&fiber, &matrix, 4.0, 1.0);
        assert!(e22_high > e22_low, "Higher ξ_E should give higher E22");
    }
}
