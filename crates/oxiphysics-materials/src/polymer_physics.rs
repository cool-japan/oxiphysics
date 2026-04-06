// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polymer physics — chain models and rubber elasticity.
//!
//! Provides the Freely Jointed Chain (FJC), Worm-Like Chain (WLC),
//! Rouse/Zimm relaxation times, rubber elasticity (neo-Hookean), and
//! Flory-Huggins mixing thermodynamics.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 1. Langevin function and inverse
// ---------------------------------------------------------------------------

/// Langevin function L(x) = coth(x) − 1/x.
///
/// Returns 0 for x ≈ 0 (using the linear approximation x/3) to avoid
/// division by zero.
pub fn langevin_function(x: f64) -> f64 {
    if x.abs() < 1e-6 {
        x / 3.0
    } else {
        x.cosh() / x.sinh() - 1.0 / x
    }
}

/// Padé approximant for the inverse Langevin function.
///
/// Approximation: L⁻¹(y) ≈ y·(3 − y²) / (1 − y²), valid for |y| < 1.
pub fn inverse_langevin_approx(y: f64) -> f64 {
    let y2 = y * y;
    y * (3.0 - y2) / (1.0 - y2)
}

// ---------------------------------------------------------------------------
// 2. Freely Jointed Chain (FJC)
// ---------------------------------------------------------------------------

/// Freely Jointed Chain (FJC) model.
///
/// Describes a polymer as `n_segments` rigid links of equal length
/// `segment_length` joined by freely-rotating joints.
#[derive(Debug, Clone)]
pub struct FreelyJointedChain {
    /// Number of Kuhn segments.
    pub n_segments: usize,
    /// Kuhn segment length (m).
    pub segment_length: f64,
}

impl FreelyJointedChain {
    /// Create a new FJC model with `n` segments of length `l`.
    pub fn new(n: usize, l: f64) -> Self {
        Self {
            n_segments: n,
            segment_length: l,
        }
    }

    /// Contour (maximum extended) length L = n·b.
    pub fn contour_length(&self) -> f64 {
        self.n_segments as f64 * self.segment_length
    }

    /// Root-mean-square end-to-end distance ⟨r²⟩^½ = b·√n.
    pub fn end_to_end_rms(&self) -> f64 {
        self.segment_length * (self.n_segments as f64).sqrt()
    }

    /// Force required to extend the chain to fractional extension `x = r/L`.
    ///
    /// Uses the inverse Langevin approximation: f = (kT / b) · L⁻¹(x).
    ///
    /// # Arguments
    /// * `x`  – fractional extension r/L in (0, 1).
    /// * `kt` – thermal energy k_B T (J).
    pub fn force_extension(&self, x: f64, kt: f64) -> f64 {
        let x_clamped = x.clamp(1e-6, 1.0 - 1e-6);
        (kt / self.segment_length) * inverse_langevin_approx(x_clamped)
    }

    /// Langevin function evaluated at `x`.
    pub fn langevin(x: f64) -> f64 {
        langevin_function(x)
    }
}

// ---------------------------------------------------------------------------
// 3. Worm-Like Chain (WLC)
// ---------------------------------------------------------------------------

/// Worm-Like Chain (WLC) model — semiflexible polymer.
///
/// Described by persistence length `lp` and contour length `lc`.
#[derive(Debug, Clone)]
pub struct WormLikeChain {
    /// Persistence length (m).
    pub persistence_length: f64,
    /// Contour length (m).
    pub contour_length: f64,
}

impl WormLikeChain {
    /// Create a new WLC model.
    pub fn new(lp: f64, lc: f64) -> Self {
        Self {
            persistence_length: lp,
            contour_length: lc,
        }
    }

    /// RMS end-to-end distance for the WLC: ⟨r²⟩^½ = √(2·lp·lc).
    ///
    /// This is the Gaussian limit; always less than contour length.
    pub fn end_to_end_rms(&self) -> f64 {
        (2.0 * self.persistence_length * self.contour_length).sqrt()
    }

    /// Force at fractional extension `r/lc` via the Marko-Siggia interpolation.
    ///
    /// f = (kT / lp) · \[ 1/(4(1−x)²) − 1/4 + x \]  where x = extension/lc.
    ///
    /// # Arguments
    /// * `extension` – end-to-end distance (m), must be < contour_length.
    /// * `kt`        – thermal energy k_B T (J).
    pub fn force_extension(&self, extension: f64, kt: f64) -> f64 {
        let lc = self.contour_length;
        let x = (extension / lc).clamp(0.0, 1.0 - 1e-6);
        let factor = kt / self.persistence_length;
        factor * (0.25 / (1.0 - x).powi(2) - 0.25 + x)
    }
}

// ---------------------------------------------------------------------------
// 4. Rouse and Zimm relaxation times
// ---------------------------------------------------------------------------

/// Rouse relaxation time τ_R = n² · b² · η / (3π² · kT).
///
/// # Arguments
/// * `n`   – number of Rouse beads.
/// * `eta` – solvent viscosity (Pa·s).
/// * `b`   – bead spacing / Kuhn length (m).
/// * `kt`  – thermal energy k_B T (J).
pub fn rouse_relaxation_time(n: usize, eta: f64, b: f64, kt: f64) -> f64 {
    let n_f = n as f64;
    n_f * n_f * b * b * eta / (3.0 * PI * PI * kt)
}

/// Zimm relaxation time τ_Z ≈ η · Rg³ / kT.
///
/// # Arguments
/// * `n`   – number of beads (used to scale Rg if desired; here Rg is explicit).
/// * `eta` – solvent viscosity (Pa·s).
/// * `rg`  – radius of gyration (m).
/// * `kt`  – thermal energy k_B T (J).
pub fn zimm_relaxation_time(_n: usize, eta: f64, rg: f64, kt: f64) -> f64 {
    eta * rg.powi(3) / kt
}

// ---------------------------------------------------------------------------
// 5. Rubber elasticity (neo-Hookean)
// ---------------------------------------------------------------------------

/// Neo-Hookean rubber elasticity model.
///
/// Models a crosslinked elastomer network in terms of chain density,
/// thermal energy, and crosslink density.
#[derive(Debug, Clone)]
pub struct RubberElasticity {
    /// Number of network chains per unit volume (m⁻³).
    pub n_chains: f64,
    /// Thermal energy k_B T (J).
    pub kt: f64,
    /// Crosslink density ν (m⁻³); used in bulk modulus estimate.
    pub crosslink_density: f64,
}

impl RubberElasticity {
    /// Create a new rubber elasticity model.
    pub fn new(n_chains: f64, kt: f64, crosslink_density: f64) -> Self {
        Self {
            n_chains,
            kt,
            crosslink_density,
        }
    }

    /// Shear modulus G = n·kT  (Pa).
    pub fn shear_modulus(&self) -> f64 {
        self.n_chains * self.kt
    }

    /// Bulk modulus K ≈ (2/3)·G + ν·kT (compressibility contribution).
    pub fn bulk_modulus(&self) -> f64 {
        let g = self.shear_modulus();
        2.0 / 3.0 * g + self.crosslink_density * self.kt
    }

    /// Strain energy density W = (G/2)·(λ₁² + λ₂² + λ₃² − 3)  for uniaxial
    /// stretch λ (incompressible: λ₂ = λ₃ = 1/√λ).
    ///
    /// # Arguments
    /// * `lambda` – uniaxial stretch ratio (≥ 1 for extension).
    pub fn strain_energy(&self, lambda: f64) -> f64 {
        let g = self.shear_modulus();
        // I1 = λ² + 2/λ for incompressible uniaxial deformation.
        let i1 = lambda * lambda + 2.0 / lambda;
        g / 2.0 * (i1 - 3.0)
    }

    /// Nominal (1st Piola-Kirchhoff) stress for uniaxial stretch λ.
    ///
    /// σ = G·(λ − 1/λ²).
    pub fn stress_stretch(lambda: f64) -> f64 {
        // G factor removed — returns normalised stress / G.
        lambda - 1.0 / (lambda * lambda)
    }
}

// ---------------------------------------------------------------------------
// 6. Flory-Huggins mixing thermodynamics
// ---------------------------------------------------------------------------

/// Flory-Huggins free energy of mixing per lattice site (units of kT).
///
/// ΔF/(nkT) = φ·ln(φ)/N₁ + (1−φ)·ln(1−φ)/N₂ + χ·φ·(1−φ)
///
/// # Arguments
/// * `phi`  – volume fraction of component 1.
/// * `n`    – degree of polymerisation of component 1 (use 1 for small molecule).
/// * `chi`  – Flory-Huggins interaction parameter.
pub fn flory_huggins_free_energy(phi: f64, n: usize, chi: f64) -> f64 {
    let phi = phi.clamp(1e-10, 1.0 - 1e-10);
    let n_f = n as f64;
    phi * phi.ln() / n_f + (1.0 - phi) * (1.0 - phi).ln() + chi * phi * (1.0 - phi)
}

/// Spinodal compositions for a symmetric blend (N₁ = N₂ = N).
///
/// Spinodal: d²ΔF/dφ² = 0  ⟹  1/(N·φ) + 1/(N·(1−φ)) − 2χ = 0.
/// Returns (φ_low, φ_high) spinodal points.
///
/// # Arguments
/// * `chi` – Flory-Huggins parameter.
/// * `n1`  – degree of polymerisation of component 1.
/// * `n2`  – degree of polymerisation of component 2.
pub fn spinodal_composition(chi: f64, n1: usize, n2: usize) -> (f64, f64) {
    // Numerically scan for spinodal points (d²F/dφ² = 0).
    let n1_f = n1 as f64;
    let n2_f = n2 as f64;
    let d2f = |phi: f64| -> f64 {
        let phi = phi.clamp(1e-10, 1.0 - 1e-10);
        1.0 / (n1_f * phi) + 1.0 / (n2_f * (1.0 - phi)) - 2.0 * chi
    };
    // Search in (0, 0.5] and [0.5, 1).
    let find_root = |a: f64, b: f64| -> f64 {
        let mut lo = a;
        let mut hi = b;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if d2f(lo) * d2f(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    };

    let chi_crit = chi_critical(n1, n2);
    if chi <= chi_crit {
        // No spinodal – return symmetric midpoint.
        (0.5, 0.5)
    } else {
        let phi_low = find_root(1e-6, 0.5 - 1e-6);
        let phi_high = find_root(0.5 + 1e-6, 1.0 - 1e-6);
        (phi_low, phi_high)
    }
}

/// Critical Flory-Huggins parameter χ_c = (1/√N₁ + 1/√N₂)² / 2.
///
/// Phase separation occurs when χ > χ_c.
pub fn chi_critical(n1: usize, n2: usize) -> f64 {
    let term = 1.0 / (n1 as f64).sqrt() + 1.0 / (n2 as f64).sqrt();
    term * term / 2.0
}

// ---------------------------------------------------------------------------
// 7. Additional polymer formulae
// ---------------------------------------------------------------------------

/// Ideal-chain radius of gyration Rg = b·√(N/6).
///
/// # Arguments
/// * `n` – number of segments.
/// * `b` – segment length (m).
pub fn radius_of_gyration_ideal(n: usize, b: f64) -> f64 {
    b * ((n as f64) / 6.0).sqrt()
}

/// Mark-Houwink intrinsic viscosity \[η\] = K · M^α.
///
/// # Arguments
/// * `k`     – Mark-Houwink K constant (mL/g).
/// * `alpha` – Mark-Houwink exponent.
/// * `mw`    – molecular weight (g/mol).
pub fn polymer_viscosity_mark_houwink(k: f64, alpha: f64, mw: f64) -> f64 {
    k * mw.powf(alpha)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── FJC ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_fjc_contour_length() {
        let fjc = FreelyJointedChain::new(100, 0.38e-9);
        let expected = 100.0 * 0.38e-9;
        assert!((fjc.contour_length() - expected).abs() < EPS);
    }

    #[test]
    fn test_fjc_end_to_end_rms_scales_sqrt_n() {
        // r_rms = b * sqrt(n), so ratio r_rms / b = sqrt(n)
        let b = 0.38e-9_f64;
        let fjc100 = FreelyJointedChain::new(100, b);
        let fjc400 = FreelyJointedChain::new(400, b);
        let ratio = fjc400.end_to_end_rms() / fjc100.end_to_end_rms();
        assert!((ratio - 2.0).abs() < 1e-10, "ratio={ratio}");
    }

    #[test]
    fn test_fjc_rms_smaller_than_contour() {
        let fjc = FreelyJointedChain::new(1000, 1e-9);
        assert!(fjc.end_to_end_rms() < fjc.contour_length());
    }

    #[test]
    fn test_fjc_force_extension_positive() {
        let fjc = FreelyJointedChain::new(100, 1e-9);
        let f = fjc.force_extension(0.5, 4.1e-21);
        assert!(f > 0.0, "force should be positive, got {f}");
    }

    #[test]
    fn test_fjc_force_extension_increases() {
        let fjc = FreelyJointedChain::new(100, 1e-9);
        let kt = 4.1e-21_f64;
        let f1 = fjc.force_extension(0.3, kt);
        let f2 = fjc.force_extension(0.8, kt);
        assert!(f2 > f1, "force should increase with extension");
    }

    // ── Langevin function ────────────────────────────────────────────────────

    #[test]
    fn test_langevin_zero() {
        // L(0) = 0 via L'Hopital
        assert!(langevin_function(0.0).abs() < 1e-6);
    }

    #[test]
    fn test_langevin_approaches_one_for_large_x() {
        // L(x) → 1 as x → ∞; at x=100 the value is very close to 1.
        let val = langevin_function(100.0);
        assert!(
            val > 0.98 && val <= 1.0,
            "L(100)={val} should be in (0.98, 1]"
        );
    }

    #[test]
    fn test_langevin_positive_for_positive_x() {
        assert!(langevin_function(1.0) > 0.0);
    }

    #[test]
    fn test_langevin_odd_function() {
        let val_pos = langevin_function(2.0);
        let val_neg = langevin_function(-2.0);
        assert!((val_pos + val_neg).abs() < EPS, "L should be odd");
    }

    #[test]
    fn test_langevin_less_than_one() {
        assert!(langevin_function(5.0) < 1.0);
    }

    // ── Inverse Langevin ────────────────────────────────────────────────────

    #[test]
    fn test_inverse_langevin_near_zero() {
        // L⁻¹(0) = 0
        assert!(inverse_langevin_approx(0.0).abs() < EPS);
    }

    #[test]
    fn test_inverse_langevin_increases() {
        let v1 = inverse_langevin_approx(0.3);
        let v2 = inverse_langevin_approx(0.6);
        assert!(v2 > v1);
    }

    #[test]
    fn test_inverse_langevin_large_for_y_near_one() {
        let v = inverse_langevin_approx(0.99);
        assert!(v > 10.0, "inverse Langevin near 1 should be large, got {v}");
    }

    // ── WLC ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_wlc_end_to_end_rms_less_than_contour() {
        let wlc = WormLikeChain::new(50e-9, 1000e-9);
        assert!(wlc.end_to_end_rms() < wlc.contour_length);
    }

    #[test]
    fn test_wlc_end_to_end_rms_formula() {
        let lp = 50e-9_f64;
        let lc = 1000e-9_f64;
        let wlc = WormLikeChain::new(lp, lc);
        let expected = (2.0 * lp * lc).sqrt();
        assert!((wlc.end_to_end_rms() - expected).abs() < EPS);
    }

    #[test]
    fn test_wlc_force_extension_positive() {
        let wlc = WormLikeChain::new(50e-9, 1000e-9);
        let f = wlc.force_extension(500e-9, 4.1e-21);
        assert!(f > 0.0, "WLC force should be positive, got {f}");
    }

    #[test]
    fn test_wlc_force_extension_increases() {
        let wlc = WormLikeChain::new(50e-9, 1000e-9);
        let kt = 4.1e-21_f64;
        let f1 = wlc.force_extension(200e-9, kt);
        let f2 = wlc.force_extension(900e-9, kt);
        assert!(f2 > f1, "WLC force should increase with extension");
    }

    // ── Relaxation times ────────────────────────────────────────────────────

    #[test]
    fn test_rouse_relaxation_scales_n_squared() {
        let eta = 1e-3_f64;
        let b = 1e-9_f64;
        let kt = 4.1e-21_f64;
        let t100 = rouse_relaxation_time(100, eta, b, kt);
        let t200 = rouse_relaxation_time(200, eta, b, kt);
        let ratio = t200 / t100;
        assert!((ratio - 4.0).abs() < 1e-9, "ratio={ratio}");
    }

    #[test]
    fn test_rouse_relaxation_positive() {
        let t = rouse_relaxation_time(50, 1e-3, 1e-9, 4.1e-21);
        assert!(t > 0.0);
    }

    #[test]
    fn test_zimm_relaxation_positive() {
        let t = zimm_relaxation_time(100, 1e-3, 10e-9, 4.1e-21);
        assert!(t > 0.0);
    }

    #[test]
    fn test_zimm_relaxation_scales_rg_cubed() {
        let eta = 1e-3_f64;
        let kt = 4.1e-21_f64;
        let t1 = zimm_relaxation_time(100, eta, 10e-9, kt);
        let t2 = zimm_relaxation_time(100, eta, 20e-9, kt);
        let ratio = t2 / t1;
        assert!((ratio - 8.0).abs() < 1e-9, "ratio={ratio}");
    }

    // ── Rubber elasticity ────────────────────────────────────────────────────

    #[test]
    fn test_rubber_shear_modulus_equals_n_kt() {
        let n = 1e25_f64;
        let kt = 4.1e-21_f64;
        let re = RubberElasticity::new(n, kt, 1e22);
        let g = re.shear_modulus();
        assert!((g - n * kt).abs() < 1e-6, "G={g}, expected n*kT={}", n * kt);
    }

    #[test]
    fn test_rubber_bulk_modulus_positive() {
        let re = RubberElasticity::new(1e25, 4.1e-21, 1e22);
        assert!(re.bulk_modulus() > 0.0);
    }

    #[test]
    fn test_rubber_strain_energy_zero_at_lambda_one() {
        let re = RubberElasticity::new(1e25, 4.1e-21, 1e22);
        let w = re.strain_energy(1.0);
        assert!(w.abs() < 1e-10, "W at λ=1 should be 0, got {w}");
    }

    #[test]
    fn test_rubber_strain_energy_positive_for_stretch() {
        let re = RubberElasticity::new(1e25, 4.1e-21, 1e22);
        assert!(re.strain_energy(2.0) > 0.0);
    }

    #[test]
    fn test_rubber_stress_stretch_formula() {
        // At λ=1: stress = 1 - 1 = 0
        assert!(RubberElasticity::stress_stretch(1.0).abs() < EPS);
        // At λ=2: stress = 2 - 0.25 = 1.75
        let expected = 2.0 - 1.0 / 4.0;
        assert!((RubberElasticity::stress_stretch(2.0) - expected).abs() < EPS);
    }

    // ── Flory-Huggins ───────────────────────────────────────────────────────

    #[test]
    fn test_flory_huggins_symmetric_at_half() {
        // Symmetric blend (n1 = n2): free energy should be the same at phi and 1-phi.
        let fh_half = flory_huggins_free_energy(0.5, 100, 0.0);
        let fh_half2 = flory_huggins_free_energy(0.5, 100, 0.0);
        assert!((fh_half - fh_half2).abs() < EPS);
    }

    #[test]
    fn test_flory_huggins_chi_zero_is_entropy_only() {
        // With chi=0 the interaction term vanishes
        let phi = 0.3_f64;
        let n = 50_usize;
        let expected = phi * phi.ln() / n as f64 + (1.0 - phi) * (1.0 - phi).ln();
        let val = flory_huggins_free_energy(phi, n, 0.0);
        assert!((val - expected).abs() < EPS);
    }

    #[test]
    fn test_flory_huggins_negative_for_mixing_favorable() {
        // chi=0, purely entropic mixing is always negative (mixing lowers free energy)
        assert!(flory_huggins_free_energy(0.3, 1, 0.0) < 0.0);
    }

    // ── chi_critical ────────────────────────────────────────────────────────

    #[test]
    fn test_chi_critical_decreases_with_n() {
        let chi_n10 = chi_critical(10, 10);
        let chi_n100 = chi_critical(100, 100);
        assert!(chi_n100 < chi_n10, "chi_c should decrease with N");
    }

    #[test]
    fn test_chi_critical_symmetric_blend() {
        // For symmetric blend: chi_c = (2/sqrt(N))^2 / 2 = 2/N
        let n = 100_usize;
        let expected = 2.0 / n as f64;
        let got = chi_critical(n, n);
        assert!(
            (got - expected).abs() < 1e-10,
            "chi_c={got}, expected={expected}"
        );
    }

    #[test]
    fn test_spinodal_below_critical_chi_returns_midpoint() {
        let chi_c = chi_critical(100, 100);
        let (lo, hi) = spinodal_composition(chi_c * 0.5, 100, 100);
        assert!((lo - 0.5).abs() < 1e-6);
        assert!((hi - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_spinodal_above_critical_chi_gives_two_points() {
        let chi_c = chi_critical(100, 100);
        let (lo, hi) = spinodal_composition(chi_c * 2.0, 100, 100);
        assert!(lo < 0.5, "low spinodal should be < 0.5, got {lo}");
        assert!(hi > 0.5, "high spinodal should be > 0.5, got {hi}");
    }

    // ── radius_of_gyration_ideal ────────────────────────────────────────────

    #[test]
    fn test_radius_of_gyration_formula() {
        let n = 100_usize;
        let b = 1e-9_f64;
        let expected = b * ((n as f64) / 6.0).sqrt();
        let got = radius_of_gyration_ideal(n, b);
        assert!((got - expected).abs() < EPS);
    }

    #[test]
    fn test_radius_of_gyration_positive() {
        assert!(radius_of_gyration_ideal(100, 1e-9) > 0.0);
    }

    #[test]
    fn test_radius_of_gyration_scales_sqrt_n() {
        let b = 1e-9_f64;
        let r100 = radius_of_gyration_ideal(100, b);
        let r400 = radius_of_gyration_ideal(400, b);
        let ratio = r400 / r100;
        assert!((ratio - 2.0).abs() < 1e-10, "ratio={ratio}");
    }

    // ── Mark-Houwink ────────────────────────────────────────────────────────

    #[test]
    fn test_mark_houwink_formula() {
        let k = 1e-4_f64;
        let alpha = 0.7_f64;
        let mw = 1e6_f64;
        let expected = k * mw.powf(alpha);
        let got = polymer_viscosity_mark_houwink(k, alpha, mw);
        assert!((got - expected).abs() < EPS);
    }

    #[test]
    fn test_mark_houwink_increases_with_mw() {
        let v1 = polymer_viscosity_mark_houwink(1e-4, 0.7, 1e5);
        let v2 = polymer_viscosity_mark_houwink(1e-4, 0.7, 1e6);
        assert!(v2 > v1);
    }

    // ── FJC langevin static method ───────────────────────────────────────────

    #[test]
    fn test_fjc_langevin_static() {
        let v = FreelyJointedChain::langevin(1.0);
        assert!(v > 0.0 && v < 1.0);
    }

    // ── edge / boundary ─────────────────────────────────────────────────────

    #[test]
    fn test_fjc_single_segment() {
        let fjc = FreelyJointedChain::new(1, 1e-9);
        assert!((fjc.contour_length() - 1e-9).abs() < EPS);
        assert!((fjc.end_to_end_rms() - 1e-9).abs() < EPS);
    }

    #[test]
    fn test_wlc_rod_limit_lp_much_greater_than_lc() {
        // When lp >> lc the chain behaves as a rod: r_rms ≈ sqrt(2*lp*lc)
        let lp = 1.0_f64; // 1 m persistence length
        let lc = 1e-6_f64; // very short
        let wlc = WormLikeChain::new(lp, lc);
        let rms = wlc.end_to_end_rms();
        let expected = (2.0 * lp * lc).sqrt();
        assert!((rms - expected).abs() < EPS);
    }

    #[test]
    fn test_rubber_stress_zero_at_lambda_one_normalised() {
        // Normalised stress σ/G = λ − 1/λ² = 0 at λ=1
        assert!(RubberElasticity::stress_stretch(1.0).abs() < EPS);
    }
}
