// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Composite material failure criteria and progressive damage models.
//!
//! Provides:
//! - [`MaxStressCriterion`] — longitudinal/transverse/shear failure indices for UD plies
//! - [`TsaiHillCriterion`] — Tsai-Hill failure surface, load factor, safety margin
//! - [`TsaiWuCriterion`] — Tsai-Wu tensor criterion, interaction coefficient, biaxial loading
//! - [`PuckCriterion`] — fiber fracture (FF) and inter-fiber fracture (IFF) modes
//! - [`ProgressiveDamage`] — stiffness degradation, ply-by-ply failure, last-ply failure
//! - [`DelamCriterion`] — quadratic delamination, ERR (G_I/G_II/G_III), mixed-mode B-K law

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Stress state
// ---------------------------------------------------------------------------

/// In-plane stress state of a unidirectional ply in principal material axes.
///
/// Axes: 1 = fiber direction, 2 = transverse in-plane, 6 = in-plane shear.
#[derive(Debug, Clone, Copy)]
pub struct PlaneStress {
    /// Normal stress in fiber direction σ₁ (Pa).
    pub sigma1: f64,
    /// Normal stress transverse to fiber σ₂ (Pa).
    pub sigma2: f64,
    /// In-plane shear stress τ₁₂ (Pa).
    pub tau12: f64,
}

impl PlaneStress {
    /// Creates a new plane stress state.
    pub fn new(sigma1: f64, sigma2: f64, tau12: f64) -> Self {
        Self {
            sigma1,
            sigma2,
            tau12,
        }
    }

    /// Returns a zero stress state.
    pub fn zero() -> Self {
        Self {
            sigma1: 0.0,
            sigma2: 0.0,
            tau12: 0.0,
        }
    }

    /// Scales all stress components by a factor.
    pub fn scaled(self, factor: f64) -> Self {
        Self {
            sigma1: self.sigma1 * factor,
            sigma2: self.sigma2 * factor,
            tau12: self.tau12 * factor,
        }
    }
}

/// Interlaminar stress components at a ply interface.
#[derive(Debug, Clone, Copy)]
pub struct InterlaminaStress {
    /// Through-thickness normal stress σ_z (Pa).
    pub sigma_z: f64,
    /// Interlaminar shear stress τ_xz (Pa).
    pub tau_xz: f64,
    /// Interlaminar shear stress τ_yz (Pa).
    pub tau_yz: f64,
}

impl InterlaminaStress {
    /// Creates a new interlaminar stress state.
    pub fn new(sigma_z: f64, tau_xz: f64, tau_yz: f64) -> Self {
        Self {
            sigma_z,
            tau_xz,
            tau_yz,
        }
    }
}

// ---------------------------------------------------------------------------
// Ply strength properties
// ---------------------------------------------------------------------------

/// Strength properties of a unidirectional composite ply.
#[derive(Debug, Clone)]
pub struct PlyStrength {
    /// Longitudinal tensile strength F₁ᵗ (Pa).
    pub f1t: f64,
    /// Longitudinal compressive strength F₁c (Pa, positive value).
    pub f1c: f64,
    /// Transverse tensile strength F₂ᵗ (Pa).
    pub f2t: f64,
    /// Transverse compressive strength F₂c (Pa, positive value).
    pub f2c: f64,
    /// In-plane shear strength F₁₂ (Pa).
    pub f12: f64,
    /// Through-thickness tensile strength F_zt (Pa, for delamination).
    pub fzt: f64,
    /// Through-thickness shear strength F_zs (Pa, for delamination).
    pub fzs: f64,
}

impl PlyStrength {
    /// Returns typical strength values for a carbon/epoxy UD ply (T300/914).
    pub fn carbon_epoxy_t300() -> Self {
        Self {
            f1t: 1500.0e6,
            f1c: 900.0e6,
            f2t: 50.0e6,
            f2c: 200.0e6,
            f12: 70.0e6,
            fzt: 50.0e6,
            fzs: 100.0e6,
        }
    }

    /// Returns typical strength values for glass/epoxy UD ply (E-glass/epoxy).
    pub fn glass_epoxy() -> Self {
        Self {
            f1t: 780.0e6,
            f1c: 500.0e6,
            f2t: 28.0e6,
            f2c: 130.0e6,
            f12: 50.0e6,
            fzt: 28.0e6,
            fzs: 80.0e6,
        }
    }
}

// ---------------------------------------------------------------------------
// MaxStressCriterion
// ---------------------------------------------------------------------------

/// Maximum stress failure criterion for a unidirectional ply.
///
/// Failure occurs when any stress component exceeds its corresponding
/// allowable strength.  The failure index for each mode is stress/strength;
/// the ply fails when any index ≥ 1.
#[derive(Debug, Clone)]
pub struct MaxStressCriterion {
    /// Ply strength properties.
    pub strength: PlyStrength,
}

impl MaxStressCriterion {
    /// Constructs a new MaxStressCriterion with the given strength.
    pub fn new(strength: PlyStrength) -> Self {
        Self { strength }
    }

    /// Computes the longitudinal failure index (tension or compression).
    pub fn fi_longitudinal(&self, stress: PlaneStress) -> f64 {
        if stress.sigma1 >= 0.0 {
            stress.sigma1 / self.strength.f1t
        } else {
            stress.sigma1.abs() / self.strength.f1c
        }
    }

    /// Computes the transverse failure index (tension or compression).
    pub fn fi_transverse(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            stress.sigma2 / self.strength.f2t
        } else {
            stress.sigma2.abs() / self.strength.f2c
        }
    }

    /// Computes the shear failure index.
    pub fn fi_shear(&self, stress: PlaneStress) -> f64 {
        stress.tau12.abs() / self.strength.f12
    }

    /// Returns all three failure indices as \[FI_1, FI_2, FI_12\].
    pub fn failure_indices(&self, stress: PlaneStress) -> [f64; 3] {
        [
            self.fi_longitudinal(stress),
            self.fi_transverse(stress),
            self.fi_shear(stress),
        ]
    }

    /// Returns the maximum failure index (governing mode).
    pub fn max_fi(&self, stress: PlaneStress) -> f64 {
        let fis = self.failure_indices(stress);
        fis.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Returns `true` if any failure index ≥ 1 (first-ply failure).
    pub fn has_failed(&self, stress: PlaneStress) -> bool {
        self.max_fi(stress) >= 1.0
    }

    /// Computes the load factor — the factor by which the stress must be
    /// multiplied to reach first failure (inverse of max FI).
    pub fn load_factor(&self, stress: PlaneStress) -> f64 {
        let max = self.max_fi(stress);
        if max <= 0.0 { f64::INFINITY } else { 1.0 / max }
    }
}

// ---------------------------------------------------------------------------
// TsaiHillCriterion
// ---------------------------------------------------------------------------

/// Tsai-Hill failure criterion for a unidirectional ply.
///
/// The quadratic failure surface:
/// (σ₁/F₁)² − σ₁σ₂/F₁² + (σ₂/F₂)² + (τ₁₂/F₁₂)² = 1
///
/// where F₁ and F₂ depend on the sign of σ₁ and σ₂.
#[derive(Debug, Clone)]
pub struct TsaiHillCriterion {
    /// Ply strength properties.
    pub strength: PlyStrength,
}

impl TsaiHillCriterion {
    /// Constructs a TsaiHillCriterion.
    pub fn new(strength: PlyStrength) -> Self {
        Self { strength }
    }

    /// Computes the Tsai-Hill failure index H² (dimensionless).
    ///
    /// Failure occurs when H² ≥ 1.
    pub fn failure_index_sq(&self, stress: PlaneStress) -> f64 {
        let f1 = if stress.sigma1 >= 0.0 {
            self.strength.f1t
        } else {
            self.strength.f1c
        };
        let f2 = if stress.sigma2 >= 0.0 {
            self.strength.f2t
        } else {
            self.strength.f2c
        };
        let s1 = stress.sigma1 / f1;
        let s2 = stress.sigma2 / f2;
        let s12 = stress.tau12 / self.strength.f12;
        // The interaction term is σ₁σ₂ / F₁² (not σ₁σ₂ / (F₁F₂)).
        s1 * s1 - stress.sigma1 * stress.sigma2 / (f1 * f1) + s2 * s2 + s12 * s12
    }

    /// Returns the Tsai-Hill failure index H (= √H²).
    pub fn failure_index(&self, stress: PlaneStress) -> f64 {
        self.failure_index_sq(stress).sqrt()
    }

    /// Returns `true` if Tsai-Hill criterion is satisfied (failure H² ≥ 1).
    pub fn has_failed(&self, stress: PlaneStress) -> bool {
        self.failure_index_sq(stress) >= 1.0
    }

    /// Computes the safety margin R = 1/H (R > 1 means no failure).
    pub fn safety_margin(&self, stress: PlaneStress) -> f64 {
        let h = self.failure_index(stress);
        if h <= 0.0 { f64::INFINITY } else { 1.0 / h }
    }

    /// Computes the load factor RF such that stress × RF is at the failure surface.
    ///
    /// From H²(RF·σ) = 1:  RF = 1 / √H²(σ)
    pub fn load_factor(&self, stress: PlaneStress) -> f64 {
        let h2 = self.failure_index_sq(stress);
        if h2 <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / h2.sqrt()
        }
    }

    /// Finds the biaxial strength envelope at a given stress ratio σ₂/σ₁ = k.
    ///
    /// Returns σ₁ > 0 (fibre-direction tension) at failure along the ray
    /// σ₂ = k·σ₁, so F₁ = F₁ᵗ and F₂ is F₂ᵗ for k ≥ 0 (σ₂ in tension) or F₂c
    /// for k < 0 (σ₂ in compression).
    pub fn biaxial_strength(&self, k: f64) -> f64 {
        // H²(σ₁, k·σ₁, 0) = 1
        // A σ₁² = 1  →  σ₁ = 1/√A
        let f1 = self.strength.f1t;
        let f2 = if k >= 0.0 {
            self.strength.f2t
        } else {
            self.strength.f2c
        };
        let a = 1.0 / (f1 * f1) - k / (f1 * f1) + k * k / (f2 * f2);
        if a <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / a.sqrt()
        }
    }
}

// ---------------------------------------------------------------------------
// TsaiWuCriterion
// ---------------------------------------------------------------------------

/// Tsai-Wu tensor failure criterion for unidirectional composites.
///
/// F_i σ_i + F_ij σ_i σ_j = 1  (Einstein summation, i,j = 1,2,6)
///
/// where the strength tensors are determined from uniaxial and shear tests.
#[derive(Debug, Clone)]
pub struct TsaiWuCriterion {
    /// Ply strength properties.
    pub strength: PlyStrength,
    /// Biaxial interaction coefficient F₁₂* (default −0.5 based on Tsai-Wu normalization).
    pub f12_star: f64,
}

impl TsaiWuCriterion {
    /// Constructs a TsaiWuCriterion with the given strength and interaction coefficient.
    pub fn new(strength: PlyStrength, f12_star: f64) -> Self {
        Self { strength, f12_star }
    }

    /// Returns the first-order tensor components \[F₁, F₂\].
    pub fn linear_terms(&self) -> [f64; 2] {
        let f1 = 1.0 / self.strength.f1t - 1.0 / self.strength.f1c;
        let f2 = 1.0 / self.strength.f2t - 1.0 / self.strength.f2c;
        [f1, f2]
    }

    /// Returns the second-order tensor components \[F₁₁, F₂₂, F₆₆, F₁₂\].
    pub fn quadratic_terms(&self) -> [f64; 4] {
        let f11 = 1.0 / (self.strength.f1t * self.strength.f1c);
        let f22 = 1.0 / (self.strength.f2t * self.strength.f2c);
        let f66 = 1.0 / (self.strength.f12 * self.strength.f12);
        let f12 = self.f12_star
            / (self.strength.f1t * self.strength.f1c * self.strength.f2t * self.strength.f2c)
                .sqrt();
        [f11, f22, f66, f12]
    }

    /// Computes the Tsai-Wu failure index (scalar value; failure when ≥ 1).
    pub fn failure_index(&self, stress: PlaneStress) -> f64 {
        let [f1, f2] = self.linear_terms();
        let [f11, f22, f66, f12] = self.quadratic_terms();
        let s1 = stress.sigma1;
        let s2 = stress.sigma2;
        let s12 = stress.tau12;
        f1 * s1 + f2 * s2 + f11 * s1 * s1 + f22 * s2 * s2 + f66 * s12 * s12 + 2.0 * f12 * s1 * s2
    }

    /// Returns `true` if the Tsai-Wu criterion is violated.
    pub fn has_failed(&self, stress: PlaneStress) -> bool {
        self.failure_index(stress) >= 1.0
    }

    /// Computes the load factor RF using the quadratic formula.
    ///
    /// RF is the positive root of: H(RF·σ) = 1.
    pub fn load_factor(&self, stress: PlaneStress) -> f64 {
        let [f1, f2] = self.linear_terms();
        let [f11, f22, f66, f12] = self.quadratic_terms();
        let s1 = stress.sigma1;
        let s2 = stress.sigma2;
        let s12 = stress.tau12;
        // A RF² + B RF − 1 = 0
        let a = f11 * s1 * s1 + f22 * s2 * s2 + f66 * s12 * s12 + 2.0 * f12 * s1 * s2;
        let b = f1 * s1 + f2 * s2;
        if a.abs() < 1e-40 {
            if b.abs() < 1e-40 {
                return f64::INFINITY;
            }
            return 1.0 / b;
        }
        let disc = b * b + 4.0 * a;
        if disc < 0.0 {
            return f64::INFINITY;
        }
        (-b + disc.sqrt()) / (2.0 * a)
    }

    /// Checks the stability condition |F₁₂| < √(F₁₁ F₂₂).
    pub fn is_stable(&self) -> bool {
        let [f11, f22, _, f12] = self.quadratic_terms();
        f12.abs() < (f11 * f22).sqrt()
    }
}

// ---------------------------------------------------------------------------
// PuckCriterion
// ---------------------------------------------------------------------------

/// Puck fiber fracture (FF) mode indicator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PuckFfMode {
    /// No fiber fracture.
    None,
    /// Fiber fracture in tension (σ₁ > 0).
    Tension,
    /// Fiber fracture in compression (σ₁ < 0) — kinking / splitting.
    Compression,
}

/// Puck inter-fiber fracture (IFF) mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PuckIffMode {
    /// No IFF.
    None,
    /// Mode A — transverse tension fracture plane ψ = 0°.
    ModeA,
    /// Mode B — in-plane shear dominated, fracture angle ψ ≈ 0°.
    ModeB,
    /// Mode C — transverse compression + shear, oblique fracture plane.
    ModeC,
}

/// Puck action-plane failure criterion for unidirectional composites.
///
/// Separates fiber fracture (FF) from inter-fiber fracture (IFF) and
/// identifies the fracture plane orientation.
#[derive(Debug, Clone)]
pub struct PuckCriterion {
    /// Ply strength properties.
    pub strength: PlyStrength,
    /// Inclination parameter p₊₁₂ for IFF Mode A (≈ 0.30 for C/E).
    pub p_plus_12: f64,
    /// Inclination parameter p₋₁₂ for IFF Modes B/C (≈ 0.25 for C/E).
    pub p_minus_12: f64,
    /// Fiber modulus parallel E₁ (Pa).
    pub e1: f64,
    /// Fiber Poisson's ratio ν₁₂.
    pub nu12: f64,
}

impl PuckCriterion {
    /// Constructs a PuckCriterion with given material parameters.
    pub fn new(strength: PlyStrength, p_plus_12: f64, p_minus_12: f64, e1: f64, nu12: f64) -> Self {
        Self {
            strength,
            p_plus_12,
            p_minus_12,
            e1,
            nu12,
        }
    }

    /// Returns default Puck parameters for carbon/epoxy T300/914.
    pub fn carbon_epoxy_default() -> Self {
        Self::new(PlyStrength::carbon_epoxy_t300(), 0.30, 0.25, 130.0e9, 0.28)
    }

    /// Computes the fiber fracture (FF) exposure f_E,FF.
    ///
    /// f_E,FF < 1 → no FF;  f_E,FF ≥ 1 → fiber fracture.
    pub fn ff_exposure(&self, stress: PlaneStress) -> f64 {
        if stress.sigma1 >= 0.0 {
            // tension
            (stress.sigma1 / self.strength.f1t).abs()
        } else {
            // compression — kink band criterion
            (stress.sigma1.abs() / self.strength.f1c).abs()
        }
    }

    /// Returns the FF mode.
    pub fn ff_mode(&self, stress: PlaneStress) -> PuckFfMode {
        let fe = self.ff_exposure(stress);
        if fe < 1.0 {
            PuckFfMode::None
        } else if stress.sigma1 >= 0.0 {
            PuckFfMode::Tension
        } else {
            PuckFfMode::Compression
        }
    }

    /// Action-plane fracture parameters R⊥⊥ᴬ, p⊥⊥⁻ and τ₂₁c that govern the
    /// transverse-compression modes (B and C) of the Puck IFF criterion.
    ///
    /// With S₂₁ the in-plane shear strength, F₂c the transverse compressive
    /// strength and p⁻ = `p_minus_12`:
    ///
    /// - R⊥⊥ᴬ = S₂₁/(2p⁻) · (√(1 + 2p⁻F₂c/S₂₁) − 1)  (= F₂c / (1 + √(1 + 2p⁻F₂c/S₂₁)))
    /// - p⊥⊥⁻ = p⁻ · R⊥⊥ᴬ / S₂₁
    /// - τ₂₁c = S₂₁ · √(1 + 2p⊥⊥⁻)
    ///
    /// These satisfy F₂c = 2(1 + p⊥⊥⁻)·R⊥⊥ᴬ.
    pub fn action_plane_parameters(&self) -> PuckActionPlaneParameters {
        let s21 = self.strength.f12;
        let f2c = self.strength.f2c;
        let p = self.p_minus_12;
        // Rationalised form of S21/(2p)(sqrt(1+x)-1), finite for p -> 0.
        let r_perp_perp_a = f2c / (1.0 + (1.0 + 2.0 * p * f2c / s21).sqrt());
        let p_perp_perp = p * r_perp_perp_a / s21;
        let tau21_c = s21 * (1.0 + 2.0 * p_perp_perp).sqrt();
        PuckActionPlaneParameters {
            r_perp_perp_a,
            p_perp_perp,
            tau21_c,
        }
    }

    /// Returns `true` when σ₂ < 0 lies in the Mode B domain
    /// (|σ₂/τ₂₁| ≤ R⊥⊥ᴬ/τ₂₁c) rather than the Mode C domain.
    fn is_mode_b_domain(&self, stress: PlaneStress) -> bool {
        let prm = self.action_plane_parameters();
        stress.sigma2.abs() * prm.tau21_c <= prm.r_perp_perp_a * stress.tau12.abs()
    }

    /// Computes the IFF exposure for Mode A (σ₂ ≥ 0, transverse tension).
    ///
    /// f_E,IFF = √\[(τ₂₁/S₂₁)² + (1 − p⁺F₂ᵗ/S₂₁)²(σ₂/F₂ᵗ)²\] + p⁺σ₂/S₂₁
    ///
    /// The plane-stress (in-plane) form of Puck's Mode A on the action plane
    /// at θ = 0°.  It equals 1 for pure transverse tension σ₂ = F₂ᵗ and for
    /// pure shear τ₂₁ = S₂₁.  Returns 0 for σ₂ < 0.
    pub fn iff_exposure_mode_a(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 < 0.0 {
            return 0.0;
        }
        let s2 = stress.sigma2;
        let s12 = stress.tau12;
        let s21 = self.strength.f12;
        let f2t = self.strength.f2t;
        let tension_weight = 1.0 - self.p_plus_12 * f2t / s21;
        let sq = ((s12 / s21).powi(2) + (tension_weight * s2 / f2t).powi(2)).sqrt();
        sq + self.p_plus_12 * s2 / s21
    }

    /// Computes the IFF exposure for Mode B (σ₂ < 0, shear dominated).
    ///
    /// f_E,IFF = (1/S₂₁) · (√\[τ₂₁² + (p⁻σ₂)²\] + p⁻σ₂)
    ///
    /// Valid in the Mode B domain |σ₂/τ₂₁| ≤ R⊥⊥ᴬ/τ₂₁c; returns 0 for σ₂ ≥ 0.
    /// It is identically 0 at τ₂₁ = 0, which is why pure transverse compression
    /// is handled by [`PuckCriterion::iff_exposure_mode_c`].
    pub fn iff_exposure_mode_b(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            return 0.0;
        }
        let s2 = stress.sigma2; // negative
        let s12 = stress.tau12;
        let s21 = self.strength.f12;
        let sq = (s12 * s12 + (self.p_minus_12 * s2).powi(2)).sqrt();
        (sq + self.p_minus_12 * s2) / s21
    }

    /// Computes the IFF exposure for Mode C (σ₂ < 0, compression dominated).
    ///
    /// f_E,IFF = \[(τ₂₁/(2(1 + p⊥⊥⁻)S₂₁))² + (σ₂/F₂c)²\] · F₂c/(−σ₂)
    ///
    /// Valid in the Mode C domain |τ₂₁/σ₂| ≤ τ₂₁c/R⊥⊥ᴬ; returns 0 for σ₂ ≥ 0.
    /// For pure transverse compression it reduces to |σ₂|/F₂c.
    pub fn iff_exposure_mode_c(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            return 0.0;
        }
        let prm = self.action_plane_parameters();
        let s2 = stress.sigma2; // negative
        let shear_term = stress.tau12 / (2.0 * (1.0 + prm.p_perp_perp) * self.strength.f12);
        let comp_term = s2 / self.strength.f2c;
        (shear_term * shear_term + comp_term * comp_term) * self.strength.f2c / (-s2)
    }

    /// Computes the IFF exposure for the transverse-compression modes (σ₂ < 0).
    ///
    /// Selects Mode B when |σ₂/τ₂₁| ≤ R⊥⊥ᴬ/τ₂₁c and Mode C otherwise, so pure
    /// transverse compression (τ₂₁ = 0) is governed by Mode C and reaches 1 at
    /// σ₂ = −F₂c.  Returns 0 for σ₂ ≥ 0.
    pub fn iff_exposure_mode_bc(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            return 0.0;
        }
        if self.is_mode_b_domain(stress) {
            self.iff_exposure_mode_b(stress)
        } else {
            self.iff_exposure_mode_c(stress)
        }
    }

    /// Returns the maximum IFF exposure considering all modes.
    pub fn iff_exposure(&self, stress: PlaneStress) -> f64 {
        let fa = self.iff_exposure_mode_a(stress);
        let fbc = self.iff_exposure_mode_bc(stress);
        fa.max(fbc)
    }

    /// Identifies the IFF mode at the given stress state.
    pub fn iff_mode(&self, stress: PlaneStress) -> PuckIffMode {
        if self.iff_exposure(stress) < 1.0 {
            return PuckIffMode::None;
        }
        if stress.sigma2 >= 0.0 {
            PuckIffMode::ModeA
        } else if self.is_mode_b_domain(stress) {
            PuckIffMode::ModeB
        } else {
            PuckIffMode::ModeC
        }
    }

    /// Returns `true` if either FF or IFF criterion is violated.
    pub fn has_failed(&self, stress: PlaneStress) -> bool {
        self.ff_exposure(stress) >= 1.0 || self.iff_exposure(stress) >= 1.0
    }

    /// Fracture plane angle θ_fp (degrees) of the IFF action plane.
    ///
    /// θ_fp = 0° in Modes A and B.  In Mode C (Puck & Schürmann 2002)
    ///
    /// cos θ_fp = √\[((R⊥⊥ᴬ/S₂₁)²(τ₂₁/σ₂)² + 1) / (2(1 + p⊥⊥⁻))\]
    ///
    /// which is the plane that maximises the action-plane exposure; it equals
    /// arccos(1/√(2(1 + p⊥⊥⁻))) ≈ 51° for pure transverse compression and falls
    /// to 0° at the Mode B/C transition.
    pub fn fracture_plane_angle(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 || self.is_mode_b_domain(stress) {
            return 0.0;
        }
        let prm = self.action_plane_parameters();
        let shear_ratio = stress.tau12 / stress.sigma2;
        let weight = (prm.r_perp_perp_a / self.strength.f12).powi(2) * shear_ratio * shear_ratio;
        let cos_sq = (weight + 1.0) / (2.0 * (1.0 + prm.p_perp_perp));
        cos_sq.sqrt().clamp(0.0, 1.0).acos().to_degrees()
    }
}

/// Action-plane fracture parameters of the Puck IFF criterion for σ₂ < 0.
///
/// Obtained from [`PuckCriterion::action_plane_parameters`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PuckActionPlaneParameters {
    /// Fracture resistance R⊥⊥ᴬ of the action plane against shear due to
    /// transverse compression (Pa).
    pub r_perp_perp_a: f64,
    /// Inclination parameter p⊥⊥⁻ of the (⊥⊥) fracture envelope (dimensionless).
    pub p_perp_perp: f64,
    /// Shear stress τ₂₁c at the transition between Mode B and Mode C (Pa).
    pub tau21_c: f64,
}

// ---------------------------------------------------------------------------
// ProgressiveDamage
// ---------------------------------------------------------------------------

/// Stiffness degradation rule applied after ply failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegradationRule {
    /// Total degradation — the ply stiffness is knocked down to 1 % after
    /// fibre fracture and to 10 % after inter-fibre fracture.
    Total,
    /// Partial degradation — the ply stiffness is knocked down to 50 % after
    /// fibre fracture and to 70 % after inter-fibre fracture.
    Partial,
    /// No degradation — track failure but keep original stiffness.
    None,
}

/// State of a single ply in a laminate.
#[derive(Debug, Clone)]
pub struct PlyState {
    /// Ply angle (degrees from laminate x-axis).
    pub angle_deg: f64,
    /// Ply thickness (m).
    pub thickness: f64,
    /// Current (degraded) in-plane stiffness matrix Q \[Q11, Q22, Q12, Q66\] (Pa),
    /// always equal to `q_virgin · degradation`.
    pub q_matrix: [f64; 4],
    /// Virgin (undamaged) in-plane stiffness matrix \[Q11, Q22, Q12, Q66\] (Pa).
    ///
    /// Ply stresses are evaluated with this matrix when they are compared with
    /// the strengths of the intact material.
    pub q_virgin: [f64; 4],
    /// Whether fiber fracture has occurred.
    pub ff_failed: bool,
    /// Whether inter-fiber fracture has occurred.
    pub iff_failed: bool,
    /// Cumulative stiffness knock-down factor (0..1, 1 = no degradation).
    pub degradation: f64,
}

impl PlyState {
    /// Constructs a PlyState from elastic constants.
    ///
    /// E₁, E₂ (Pa), ν₁₂, G₁₂ (Pa).
    pub fn new(angle_deg: f64, thickness: f64, e1: f64, e2: f64, nu12: f64, g12: f64) -> Self {
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        let q11 = e1 / denom;
        let q22 = e2 / denom;
        let q12 = nu12 * e2 / denom;
        let q66 = g12;
        let q = [q11, q22, q12, q66];
        Self {
            angle_deg,
            thickness,
            q_matrix: q,
            q_virgin: q,
            ff_failed: false,
            iff_failed: false,
            degradation: 1.0,
        }
    }

    /// Applies stiffness degradation according to rule and failure mode.
    ///
    /// The knock-down is applied to the virgin stiffness exactly once
    /// (`q_matrix = q_virgin · degradation`), so repeated calls are idempotent
    /// and damage is irreversible: the cumulative factor never increases.
    pub fn apply_degradation(&mut self, rule: DegradationRule, ff: bool, iff: bool) {
        self.ff_failed |= ff;
        self.iff_failed |= iff;
        let factor = match rule {
            DegradationRule::Total => {
                if ff {
                    0.01
                } else if iff {
                    0.1
                } else {
                    1.0
                }
            }
            DegradationRule::Partial => {
                if ff {
                    0.5
                } else if iff {
                    0.7
                } else {
                    1.0
                }
            }
            DegradationRule::None => 1.0,
        };
        self.degradation = self.degradation.min(factor);
        for (q, q0) in self.q_matrix.iter_mut().zip(self.q_virgin.iter()) {
            *q = q0 * self.degradation;
        }
    }

    /// Returns the effective Q₁₁ after degradation.
    pub fn q11(&self) -> f64 {
        self.q_matrix[0]
    }

    /// Returns the effective Q₂₂ after degradation.
    pub fn q22(&self) -> f64 {
        self.q_matrix[1]
    }

    /// Returns the effective Q₁₂ after degradation.
    pub fn q12(&self) -> f64 {
        self.q_matrix[2]
    }

    /// Returns the effective Q₆₆ after degradation.
    pub fn q66(&self) -> f64 {
        self.q_matrix[3]
    }
}

/// Progressive damage analysis of a composite laminate.
///
/// Performs ply-by-ply failure analysis with the maximum-stress criterion,
/// applies stiffness degradation, and tracks first-ply failure (FPF)
/// and last-ply failure (LPF) loads.
#[derive(Debug, Clone)]
pub struct ProgressiveDamage {
    /// Ordered list of ply states.
    pub plies: Vec<PlyState>,
    /// Degradation rule applied when a ply fails.
    pub rule: DegradationRule,
    /// Strength per ply.
    pub strengths: Vec<PlyStrength>,
    /// Load factor λ of the first ply failure found by the last
    /// [`ProgressiveDamage::progressive_failure_sweep`] (0 = no sweep has been
    /// run yet, or the last sweep activated no failure mode).
    pub fpf_load: f64,
    /// Load factor λ of the last failure event of the last sweep, after which
    /// no further failure mode can be activated (0 = no sweep has been run
    /// yet, or the last sweep activated no failure mode).
    pub lpf_load: f64,
}

impl ProgressiveDamage {
    /// Constructs a ProgressiveDamage model from ply states, strengths, and rule.
    pub fn new(plies: Vec<PlyState>, strengths: Vec<PlyStrength>, rule: DegradationRule) -> Self {
        Self {
            plies,
            rule,
            strengths,
            fpf_load: 0.0,
            lpf_load: 0.0,
        }
    }

    /// Returns the number of intact (unfailed) plies.
    pub fn n_intact(&self) -> usize {
        self.plies
            .iter()
            .filter(|p| !p.ff_failed && !p.iff_failed)
            .count()
    }

    /// Returns `true` if all plies have failed (last-ply failure).
    pub fn all_failed(&self) -> bool {
        self.plies.iter().all(|p| p.ff_failed || p.iff_failed)
    }

    /// Computes the laminate ABD A-matrix (membrane stiffness, 3×3) using CLT.
    ///
    /// Each ply contributes its current (degraded) stiffness `q_matrix`, which
    /// already carries the knock-down factor, times its thickness.
    ///
    /// Returns \[A11, A22, A12, A16, A26, A66\] (Pa·m).
    pub fn a_matrix(&self) -> [f64; 6] {
        let mut a = [0.0_f64; 6];
        for ply in &self.plies {
            let t = ply.thickness;
            let q = transformed_q(&ply.q_matrix, ply.angle_deg.to_radians());
            a[0] += q[0] * t;
            a[1] += q[1] * t;
            a[2] += q[2] * t;
            a[3] += q[3] * t;
            a[4] += q[4] * t;
            a[5] += q[5] * t;
        }
        a
    }

    /// The full symmetric 3×3 membrane stiffness matrix
    /// \[\[A11, A12, A16\], \[A12, A22, A26\], \[A16, A26, A66\]\] (Pa·m).
    pub fn a_matrix_full(&self) -> [[f64; 3]; 3] {
        let [a11, a22, a12, a16, a26, a66] = self.a_matrix();
        [[a11, a12, a16], [a12, a22, a26], [a16, a26, a66]]
    }

    /// Event-driven progressive failure analysis under the proportional
    /// membrane load `λ · applied_load`, with `applied_load = [N₁, N₂, N₆]`
    /// (force resultants, Pa·m) in laminate axes.
    ///
    /// The laminate response is linear between failure events.  At every stage
    /// the mid-plane strains per unit load follow from the full coupled
    /// system `A ε = N` (including A₁₂, A₁₆ and A₂₆), the ply stresses per unit
    /// load follow from the virgin ply stiffness, and the load factor at which
    /// each still-active failure mode of the maximum-stress criterion is
    /// reached is `strength / |σ(λ = 1)|`.  The smallest of these (never
    /// smaller than the current load level, so cascades fail at the current
    /// level) is the next event: every ply and mode reaching its strength
    /// there fails simultaneously, the ply stiffness is degraded according to
    /// [`ProgressiveDamage::rule`], `A` is re-assembled and the next event is
    /// searched, until no failure mode can be activated any more.
    ///
    /// Fibre fracture (σ₁) ends the ply; inter-fibre fracture (σ₂ or τ₁₂)
    /// leaves the fibre mode active.  Only the membrane stiffness `A` is used
    /// (bending–extension coupling `B` is neglected, which is exact for
    /// symmetric laminates).  The analysis starts from the current damage
    /// state of the plies and updates them in place.
    ///
    /// Returns the load factor at first-ply failure (`f64::INFINITY` when no
    /// failure mode is ever activated, e.g. for zero load) and sets
    /// [`ProgressiveDamage::fpf_load`] and [`ProgressiveDamage::lpf_load`].
    ///
    /// # Errors
    /// * [`crate::Error::DimensionMismatch`] if there is not exactly one
    ///   strength entry per ply.
    /// * [`crate::Error::General`] if a load component or a membrane stiffness
    ///   entry is not finite, or a strength is not finite.
    /// * [`crate::Error::ParameterOutOfRange`] if one of the in-plane strengths
    ///   F₁ᵗ, F₁c, F₂ᵗ, F₂c, F₁₂ of a ply is not strictly positive (a zero
    ///   strength would turn the failure index infinite and the mode would be
    ///   silently treated as unloaded).
    /// * [`crate::Error::Singularity`] if the membrane stiffness matrix `A` is
    ///   singular (e.g. an empty laminate).
    pub fn progressive_failure_sweep(
        &mut self,
        applied_load: [f64; 3],
    ) -> Result<f64, crate::Error> {
        if self.strengths.len() != self.plies.len() {
            return Err(crate::Error::DimensionMismatch {
                expected: format!("{} ply strengths (one per ply)", self.plies.len()),
                actual: format!("{} ply strengths", self.strengths.len()),
            });
        }
        for &component in &applied_load {
            crate::require_finite("applied_load", component)?;
        }
        for strength in &self.strengths {
            let in_plane = [
                ("f1t", strength.f1t),
                ("f1c", strength.f1c),
                ("f2t", strength.f2t),
                ("f2c", strength.f2c),
                ("f12", strength.f12),
            ];
            for (name, value) in in_plane {
                crate::require_finite(name, value)?;
                crate::require_positive(name, value)?;
            }
        }

        let criteria: Vec<MaxStressCriterion> = self
            .strengths
            .iter()
            .cloned()
            .map(MaxStressCriterion::new)
            .collect();
        let rule = self.rule;
        let mut lambda = 0.0_f64;
        let mut first_failure = f64::INFINITY;
        let mut any_failure = false;

        loop {
            let membrane_stiffness = self.a_matrix_full();
            if membrane_stiffness
                .iter()
                .flatten()
                .any(|entry| !entry.is_finite())
            {
                return Err(crate::Error::General(
                    "laminate membrane stiffness matrix A is not finite \
                     (check the ply thicknesses and elastic constants)"
                        .to_string(),
                ));
            }
            let strains_per_unit_load =
                solve_3x3(membrane_stiffness, applied_load).ok_or_else(|| {
                    crate::Error::singular(
                        "progressive_failure_sweep",
                        "laminate membrane stiffness matrix A is singular",
                    )
                })?;

            let mut lambda_min = f64::INFINITY;
            let mut loads: Vec<PlyFailureLoads> = Vec::with_capacity(self.plies.len());
            for (ply, criterion) in self.plies.iter().zip(criteria.iter()) {
                let ply_loads = PlyFailureLoads::evaluate(ply, criterion, strains_per_unit_load);
                lambda_min = lambda_min.min(ply_loads.ff).min(ply_loads.iff);
                loads.push(ply_loads);
            }
            if !lambda_min.is_finite() {
                break;
            }

            let lambda_event = lambda.max(lambda_min);
            let threshold = lambda_event * (1.0 + SIMULTANEOUS_FAILURE_TOLERANCE);
            for (ply, ply_loads) in self.plies.iter_mut().zip(loads.iter()) {
                let ff = ply_loads.ff <= threshold;
                let iff = ply_loads.iff <= threshold;
                if ff || iff {
                    ply.apply_degradation(rule, ff, iff);
                }
            }
            if !any_failure {
                first_failure = lambda_event;
                any_failure = true;
            }
            lambda = lambda_event;
        }

        // Both load factors describe this sweep only: a sweep that finds no
        // active failure mode must not leave the values of an earlier sweep.
        self.fpf_load = if any_failure { first_failure } else { 0.0 };
        self.lpf_load = if any_failure { lambda } else { 0.0 };
        Ok(first_failure)
    }
}

/// Relative tolerance within which failure events are treated as simultaneous.
const SIMULTANEOUS_FAILURE_TOLERANCE: f64 = 1.0e-9;

/// Load factors at which the fibre and inter-fibre failure modes of one ply are
/// reached, for the current damage state and a proportional load.
struct PlyFailureLoads {
    /// Load factor of fibre fracture (σ₁ against F₁ᵗ/F₁c); `INFINITY` if inactive.
    ff: f64,
    /// Load factor of inter-fibre fracture (σ₂ against F₂ᵗ/F₂c or τ₁₂ against
    /// F₁₂); `INFINITY` if inactive.
    iff: f64,
}

impl PlyFailureLoads {
    /// Evaluates the still-active failure modes of `ply` for the laminate
    /// mid-plane strains `strains` (per unit load, laminate axes).
    fn evaluate(ply: &PlyState, criterion: &MaxStressCriterion, strains: [f64; 3]) -> Self {
        if ply.ff_failed {
            // Fibre fracture exhausts the ply.
            return Self {
                ff: f64::INFINITY,
                iff: f64::INFINITY,
            };
        }
        let stress = ply_stress_from_laminate_strain(ply, strains);
        let [fi_fibre, fi_transverse, fi_shear] = criterion.failure_indices(stress);
        let iff = if ply.iff_failed {
            f64::INFINITY
        } else {
            load_factor_from_index(fi_transverse.max(fi_shear))
        };
        Self {
            ff: load_factor_from_index(fi_fibre),
            iff,
        }
    }
}

/// Load factor `1 / FI` for a failure index per unit load; `INFINITY` when the
/// mode is unloaded.
fn load_factor_from_index(failure_index: f64) -> f64 {
    if failure_index > 0.0 && failure_index.is_finite() {
        1.0 / failure_index
    } else {
        f64::INFINITY
    }
}

/// Ply-axis stress of the intact material for laminate mid-plane strains
/// `[ε_x, ε_y, γ_xy]`: strains are rotated into the ply axes and multiplied by
/// the virgin ply stiffness.
fn ply_stress_from_laminate_strain(ply: &PlyState, strains: [f64; 3]) -> PlaneStress {
    let (s, c) = ply.angle_deg.to_radians().sin_cos();
    let [eps_x, eps_y, gamma_xy] = strains;
    let eps1 = c * c * eps_x + s * s * eps_y + c * s * gamma_xy;
    let eps2 = s * s * eps_x + c * c * eps_y - c * s * gamma_xy;
    let gamma12 = -2.0 * c * s * eps_x + 2.0 * c * s * eps_y + (c * c - s * s) * gamma_xy;
    let [q11, q22, q12, q66] = ply.q_virgin;
    PlaneStress::new(
        q11 * eps1 + q12 * eps2,
        q12 * eps1 + q22 * eps2,
        q66 * gamma12,
    )
}

/// Solves the 3×3 system `A x = b` by Gaussian elimination with partial
/// pivoting.  Returns `None` when `A` is singular relative to its largest
/// entry, or when `A` or `b` has a non-finite entry.
fn solve_3x3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    if a.iter().flatten().chain(b.iter()).any(|v| !v.is_finite()) {
        return None;
    }
    let scale = a.iter().flatten().fold(0.0_f64, |acc, v| acc.max(v.abs()));
    if scale <= 0.0 {
        return None;
    }
    let pivot_tolerance = 1.0e-14 * scale;
    let mut m = [[0.0_f64; 4]; 3];
    for (row, (a_row, rhs)) in m.iter_mut().zip(a.iter().zip(b.iter())) {
        row[..3].copy_from_slice(a_row);
        row[3] = *rhs;
    }
    for col in 0..3 {
        let mut pivot_row = col;
        for r in (col + 1)..3 {
            if m[r][col].abs() > m[pivot_row][col].abs() {
                pivot_row = r;
            }
        }
        if m[pivot_row][col].abs() <= pivot_tolerance {
            return None;
        }
        m.swap(col, pivot_row);
        let pivot_values = m[col];
        for row in m.iter_mut().skip(col + 1) {
            let factor = row[col] / pivot_values[col];
            for (value, pivot_value) in row.iter_mut().zip(pivot_values.iter()).skip(col) {
                *value -= factor * pivot_value;
            }
        }
    }
    let mut x = [0.0_f64; 3];
    for i in (0..3).rev() {
        let mut acc = m[i][3];
        for j in (i + 1)..3 {
            acc -= m[i][j] * x[j];
        }
        x[i] = acc / m[i][i];
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// DelamCriterion
// ---------------------------------------------------------------------------

/// Energy release rate components for delamination.
#[derive(Debug, Clone, Copy)]
pub struct EnergyReleaseRate {
    /// Mode I (opening) ERR G_I (J/m²).
    pub g1: f64,
    /// Mode II (sliding shear) ERR G_II (J/m²).
    pub g2: f64,
    /// Mode III (tearing shear) ERR G_III (J/m²).
    pub g3: f64,
}

impl EnergyReleaseRate {
    /// Creates a new EnergyReleaseRate.
    pub fn new(g1: f64, g2: f64, g3: f64) -> Self {
        Self { g1, g2, g3 }
    }

    /// Returns the total ERR G = G_I + G_II + G_III.
    pub fn total(&self) -> f64 {
        self.g1 + self.g2 + self.g3
    }

    /// Returns the mode mix ratio G_I / G_total.
    pub fn mode_i_ratio(&self) -> f64 {
        let gt = self.total();
        if gt < 1e-20 { 0.0 } else { self.g1 / gt }
    }

    /// Returns the mode mix ratio G_II / G_total.
    pub fn mode_ii_ratio(&self) -> f64 {
        let gt = self.total();
        if gt < 1e-20 { 0.0 } else { self.g2 / gt }
    }
}

/// Delamination fracture criterion for composite laminates.
///
/// Implements:
/// 1. Quadratic stress criterion (interlaminar stresses).
/// 2. Linear fracture mechanics criterion (G vs G_c).
/// 3. Benzeggagh-Kenane (B-K) mixed-mode criterion.
#[derive(Debug, Clone)]
pub struct DelamCriterion {
    /// Interlaminar tensile strength Z_t (Pa).
    pub z_t: f64,
    /// Interlaminar shear strength S_xz (Pa).
    pub s_xz: f64,
    /// Interlaminar shear strength S_yz (Pa).
    pub s_yz: f64,
    /// Mode I critical ERR G_Ic (J/m²).
    pub g1c: f64,
    /// Mode II critical ERR G_IIc (J/m²).
    pub g2c: f64,
    /// Mode III critical ERR G_IIIc (J/m²).
    pub g3c: f64,
    /// Benzeggagh-Kenane exponent η (C/E: ≈ 1.75).
    pub bk_eta: f64,
}

impl DelamCriterion {
    /// Constructs a DelamCriterion.
    pub fn new(z_t: f64, s_xz: f64, s_yz: f64, g1c: f64, g2c: f64, g3c: f64, bk_eta: f64) -> Self {
        Self {
            z_t,
            s_xz,
            s_yz,
            g1c,
            g2c,
            g3c,
            bk_eta,
        }
    }

    /// Returns default values for carbon/epoxy (T300/914).
    pub fn carbon_epoxy_default() -> Self {
        Self::new(50.0e6, 100.0e6, 100.0e6, 200.0, 400.0, 400.0, 1.75)
    }

    /// Computes the quadratic interlaminar stress criterion.
    ///
    /// (⟨σ_z⟩/Z_t)² + (τ_xz/S_xz)² + (τ_yz/S_yz)² ≥ 1 → delamination
    ///
    /// where ⟨·⟩ is the Macaulay bracket (only tensile σ_z contributes).
    pub fn quadratic_stress_index(&self, stress: InterlaminaStress) -> f64 {
        let sz = stress.sigma_z.max(0.0); // Macaulay bracket
        (sz / self.z_t).powi(2)
            + (stress.tau_xz / self.s_xz).powi(2)
            + (stress.tau_yz / self.s_yz).powi(2)
    }

    /// Returns `true` if the quadratic stress criterion is violated.
    pub fn stress_delamination(&self, stress: InterlaminaStress) -> bool {
        self.quadratic_stress_index(stress) >= 1.0
    }

    /// Computes the linear ERR criterion.
    ///
    /// G_I/G_Ic + G_II/G_IIc + G_III/G_IIIc ≥ 1 → delamination
    pub fn linear_err_index(&self, err: EnergyReleaseRate) -> f64 {
        err.g1 / self.g1c + err.g2 / self.g2c + err.g3 / self.g3c
    }

    /// Returns `true` if the linear ERR criterion is violated.
    pub fn err_delamination(&self, err: EnergyReleaseRate) -> bool {
        self.linear_err_index(err) >= 1.0
    }

    /// Computes the mixed-mode B-K criterion (Benzeggagh–Kenane 1996).
    ///
    /// G_Ic + (G_IIc − G_Ic) * (G_shear/G_total)^η = G_c(mm)
    /// Delamination when G_total ≥ G_c(mm).
    pub fn bk_criterion(&self, err: EnergyReleaseRate) -> f64 {
        let gt = err.total();
        if gt < 1e-20 {
            return 0.0;
        }
        let g_shear = err.g2 + err.g3;
        let mix = g_shear / gt;
        let gc_mm = self.g1c + (self.g2c - self.g1c) * mix.powf(self.bk_eta);
        gt / gc_mm
    }

    /// Returns `true` if B-K criterion is violated.
    pub fn bk_delamination(&self, err: EnergyReleaseRate) -> bool {
        self.bk_criterion(err) >= 1.0
    }

    /// Computes the mode I stress intensity factor K_I from σ_z and half-crack a (m).
    ///
    /// K_I = σ_z √(π a)
    pub fn stress_intensity_mode1(&self, sigma_z: f64, half_crack: f64) -> f64 {
        sigma_z * (PI * half_crack).sqrt()
    }

    /// Returns the critical G_I from mode I fracture toughness K_Ic
    /// using G = K²/E.
    ///
    /// G_Ic = K_Ic² / E_z  where E_z is the through-thickness modulus (Pa).
    pub fn g1c_from_k1c(k1c: f64, ez: f64) -> f64 {
        k1c * k1c / ez
    }

    /// Computes the crack driving force (J-integral proxy) from ERR (J/m²).
    pub fn j_integral_proxy(&self, err: EnergyReleaseRate) -> f64 {
        err.total()
    }

    /// Estimates the delamination growth rate da/dN (m/cycle) using a Paris law.
    ///
    /// da/dN = C (ΔG_equiv)^m
    ///
    /// where ΔG_equiv = ΔG_I + ΔG_II * (G_IIc/G_Ic)^α  (approximate).
    pub fn paris_growth_rate(
        &self,
        delta_g1: f64,
        delta_g2: f64,
        paris_c: f64,
        paris_m: f64,
    ) -> f64 {
        let alpha = 0.5;
        let g_equiv = delta_g1 + delta_g2 * (self.g2c / self.g1c).powf(alpha);
        paris_c * g_equiv.powf(paris_m)
    }
}

// ---------------------------------------------------------------------------
// Helper: transformed stiffness matrix (CLT)
// ---------------------------------------------------------------------------

/// Transforms the ply stiffness matrix Q to the laminate coordinate system.
///
/// Returns \[Q̄₁₁, Q̄₂₂, Q̄₁₂, Q̄₁₆, Q̄₂₆, Q̄₆₆\].
fn transformed_q(q: &[f64; 4], theta: f64) -> [f64; 6] {
    let c = theta.cos();
    let s = theta.sin();
    let c2 = c * c;
    let s2 = s * s;
    let c4 = c2 * c2;
    let s4 = s2 * s2;
    let cs = c * s;
    let cs2 = cs * cs;
    let [q11, q22, q12, q66] = *q;
    let q_bar11 = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * cs2 + q22 * s4;
    let q_bar22 = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * cs2 + q22 * c4;
    let q_bar12 = (q11 + q22 - 4.0 * q66) * cs2 + q12 * (c4 + s4);
    let q_bar16 = (q11 - q12 - 2.0 * q66) * c2 * cs - (q22 - q12 - 2.0 * q66) * s2 * cs;
    let q_bar26 = (q11 - q12 - 2.0 * q66) * s2 * cs - (q22 - q12 - 2.0 * q66) * c2 * cs;
    let q_bar66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * cs2 + q66 * (c4 + s4);
    [q_bar11, q_bar22, q_bar12, q_bar16, q_bar26, q_bar66]
}

/// Computes the reserve factor (inverse of failure index) for a given criterion value.
///
/// Returns `None` if the criterion value is zero.
pub fn reserve_factor(criterion_value: f64) -> Option<f64> {
    if criterion_value <= 0.0 {
        None
    } else {
        Some(1.0 / criterion_value)
    }
}

/// Returns the Tsai-Hill failure envelope points for σ₁ sweeping from −F₁c to F₁t.
///
/// Each entry is (σ₁, σ₂) on the tensile (σ₂ ≥ 0) branch of the failure locus
/// at τ₁₂ = 0, so that F₂ = F₂ᵗ applies and `failure_index_sq` equals 1 at
/// every returned point.
pub fn tsai_hill_envelope(strength: &PlyStrength, n_points: usize) -> Vec<(f64, f64)> {
    let mut pts = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let sigma1 =
            -strength.f1c + (strength.f1t + strength.f1c) * i as f64 / (n_points - 1).max(1) as f64;
        // Solve: (s1/F1)² - s1*s2/F1² + (s2/F2)² = 1 for s2
        let f1 = if sigma1 >= 0.0 {
            strength.f1t
        } else {
            strength.f1c
        };
        let f2 = strength.f2t;
        let a = 1.0 / (f2 * f2);
        let b = -sigma1 / (f1 * f1);
        let c_coeff = (sigma1 / f1).powi(2) - 1.0;
        let disc = b * b - 4.0 * a * c_coeff;
        if disc >= 0.0 {
            let s2 = (-b + disc.sqrt()) / (2.0 * a);
            pts.push((sigma1, s2));
        }
    }
    pts
}

#[cfg(test)]
mod tests;
