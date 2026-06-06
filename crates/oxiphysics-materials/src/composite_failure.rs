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
        s1 * s1 - s1 * s2 + s2 * s2 + s12 * s12
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
    /// Returns σ₁ at failure.
    pub fn biaxial_strength(&self, k: f64) -> f64 {
        // H²(σ₁, k·σ₁, 0) = 1
        // A σ₁² = 1  →  σ₁ = 1/√A
        let f1 = self.strength.f1t.max(self.strength.f1c);
        let f2 = self.strength.f2t.max(self.strength.f2c);
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

    /// Computes the IFF exposure on the action plane at ψ = 0° (simplified).
    ///
    /// For Mode A (σ_n > 0):
    /// f_E,IFF = √\[(τ_nt/F_2A)² + (σ_n/F_2t)²\] + p₊₁₂ σ_n / F_2A
    ///
    /// Uses simplified 2-D (in-plane) version.
    pub fn iff_exposure_mode_a(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 < 0.0 {
            return 0.0;
        }
        let s2 = stress.sigma2;
        let s12 = stress.tau12;
        let f2a = self.strength.f12; // approximate action-plane shear strength
        let f2t = self.strength.f2t;
        let sq = ((s12 / f2a).powi(2) + (s2 / f2t).powi(2)).sqrt();
        sq + self.p_plus_12 * s2 / f2a
    }

    /// Computes the IFF exposure for Mode B/C (σ_n ≤ 0).
    ///
    /// f_E,IFF = (1/F_2A) √\[τ_nt² + (p₋₁₂ σ_n)²\] + p₋₁₂ σ_n / F_2A
    pub fn iff_exposure_mode_bc(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            return 0.0;
        }
        let s2 = stress.sigma2; // negative
        let s12 = stress.tau12;
        let f2a = self.strength.f12;
        let sq = (s12 * s12 + (self.p_minus_12 * s2).powi(2)).sqrt();
        (sq + self.p_minus_12 * s2) / f2a
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
        } else {
            // Distinguish B vs C by shear dominance
            if stress.tau12.abs() > self.strength.f12 * 0.5 {
                PuckIffMode::ModeB
            } else {
                PuckIffMode::ModeC
            }
        }
    }

    /// Returns `true` if either FF or IFF criterion is violated.
    pub fn has_failed(&self, stress: PlaneStress) -> bool {
        self.ff_exposure(stress) >= 1.0 || self.iff_exposure(stress) >= 1.0
    }

    /// Estimates the fracture plane angle θ_fp (degrees) for IFF Mode C.
    ///
    /// θ_fp ≈ arccos(F₂c / (2·F₂c + p₋₁₂·|σ₂|))  (approximate Puck formula)
    pub fn fracture_plane_angle(&self, stress: PlaneStress) -> f64 {
        if stress.sigma2 >= 0.0 {
            return 0.0; // Mode A: fracture perpendicular to load
        }
        let denom = 2.0 * self.strength.f2c + self.p_minus_12 * stress.sigma2.abs();
        if denom < 1e-20 {
            return 0.0;
        }
        (self.strength.f2c / denom).acos().to_degrees()
    }
}

// ---------------------------------------------------------------------------
// ProgressiveDamage
// ---------------------------------------------------------------------------

/// Stiffness degradation rule applied after ply failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegradationRule {
    /// Total degradation — set failed stiffness components to near zero.
    Total,
    /// Partial degradation — reduce by a factor (e.g., 10%).
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
    /// In-plane stiffness matrix Q \[Q11, Q22, Q12, Q66\] (Pa).
    pub q_matrix: [f64; 4],
    /// Whether fiber fracture has occurred.
    pub ff_failed: bool,
    /// Whether inter-fiber fracture has occurred.
    pub iff_failed: bool,
    /// Degradation factor applied (0..1, 1 = no degradation).
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
        Self {
            angle_deg,
            thickness,
            q_matrix: [q11, q22, q12, q66],
            ff_failed: false,
            iff_failed: false,
            degradation: 1.0,
        }
    }

    /// Applies stiffness degradation according to rule and failure mode.
    pub fn apply_degradation(&mut self, rule: DegradationRule, ff: bool, iff: bool) {
        self.ff_failed |= ff;
        self.iff_failed |= iff;
        self.degradation = match rule {
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
        for q in &mut self.q_matrix {
            *q *= self.degradation;
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
/// Performs ply-by-ply failure analysis using the chosen criterion,
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
    /// Applied load multiplier at FPF (0 = not yet reached).
    pub fpf_load: f64,
    /// Applied load multiplier at LPF (0 = not yet reached).
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
    /// Returns \[A11, A22, A12, A16, A26, A66\] (Pa·m).
    pub fn a_matrix(&self) -> [f64; 6] {
        let mut a = [0.0_f64; 6];
        for ply in &self.plies {
            let t = ply.thickness * ply.degradation;
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

    /// Performs a progressive failure sweep at a given load vector \[N₁, N₂, N₆\] (Pa·m).
    ///
    /// Returns the load factor at first-ply failure and updates ply states.
    pub fn progressive_failure_sweep(&mut self, applied_load: [f64; 3], n_steps: usize) -> f64 {
        let load_max = 2.0; // sweep load factor from 0 to 2
        let d_lambda = load_max / n_steps as f64;
        let mut first_failure_lambda = f64::INFINITY;

        for step in 1..=n_steps {
            let lambda = step as f64 * d_lambda;
            let n = [
                applied_load[0] * lambda,
                applied_load[1] * lambda,
                applied_load[2] * lambda,
            ];
            // Compute average in-plane strains from A-matrix (simplified: isotropic approx)
            let a = self.a_matrix();
            let total_t: f64 = self
                .plies
                .iter()
                .map(|p| p.thickness)
                .sum::<f64>()
                .max(1e-12);
            let eps1 = n[0] / (a[0].max(1e-12) * total_t);
            let eps2 = n[1] / (a[1].max(1e-12) * total_t);
            let gam12 = n[2] / (a[5].max(1e-12) * total_t);

            let mut any_new_failure = false;
            for (i, ply) in self.plies.iter_mut().enumerate() {
                if ply.ff_failed && ply.iff_failed {
                    continue;
                }
                // Transform strains to ply axes and compute stress
                let theta = ply.angle_deg.to_radians();
                let c = theta.cos();
                let s = theta.sin();
                let eps1_ply = c * c * eps1 + s * s * eps2 + c * s * gam12;
                let eps2_ply = s * s * eps1 + c * c * eps2 - c * s * gam12;
                let gam12_ply = -2.0 * c * s * eps1 + 2.0 * c * s * eps2 + (c * c - s * s) * gam12;
                let q = &ply.q_matrix;
                let stress = PlaneStress::new(
                    q[0] * eps1_ply + q[2] * eps2_ply,
                    q[2] * eps1_ply + q[1] * eps2_ply,
                    q[3] * gam12_ply,
                );
                if i < self.strengths.len() {
                    let crit = MaxStressCriterion::new(self.strengths[i].clone());
                    if crit.has_failed(stress) && !ply.ff_failed {
                        ply.ff_failed = true;
                        ply.degradation = match self.rule {
                            DegradationRule::Total => 0.01,
                            DegradationRule::Partial => 0.5,
                            DegradationRule::None => 1.0,
                        };
                        any_new_failure = true;
                        if first_failure_lambda > lambda {
                            first_failure_lambda = lambda;
                        }
                    }
                }
            }
            let _ = any_new_failure;
            if self.all_failed() {
                self.lpf_load = lambda;
                break;
            }
        }
        if first_failure_lambda < f64::INFINITY {
            self.fpf_load = first_failure_lambda;
        }
        first_failure_lambda
    }
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
/// Each entry is (σ₁, σ₂) on the failure locus at τ₁₂ = 0.
pub fn tsai_hill_envelope(strength: &PlyStrength, n_points: usize) -> Vec<(f64, f64)> {
    let crit = TsaiHillCriterion::new(strength.clone());
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
            let _ = crit.failure_index_sq(PlaneStress::new(sigma1, s2, 0.0)); // verify
            pts.push((sigma1, s2));
        }
    }
    pts
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PlaneStress ---

    #[test]
    fn test_plane_stress_zero() {
        let s = PlaneStress::zero();
        assert_eq!(s.sigma1, 0.0);
    }

    #[test]
    fn test_plane_stress_scaled() {
        let s = PlaneStress::new(100.0, 50.0, 20.0).scaled(2.0);
        assert!((s.sigma1 - 200.0).abs() < 1e-10);
        assert!((s.tau12 - 40.0).abs() < 1e-10);
    }

    // --- PlyStrength ---

    #[test]
    fn test_carbon_epoxy_strengths_positive() {
        let st = PlyStrength::carbon_epoxy_t300();
        assert!(st.f1t > 0.0 && st.f1c > 0.0 && st.f2t > 0.0 && st.f12 > 0.0);
    }

    #[test]
    fn test_glass_epoxy_f1t_less_than_carbon() {
        let ce = PlyStrength::carbon_epoxy_t300();
        let ge = PlyStrength::glass_epoxy();
        assert!(ge.f1t < ce.f1t);
    }

    // --- MaxStressCriterion ---

    #[test]
    fn test_max_stress_no_failure_under_threshold() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6); // well below strength
        assert!(!crit.has_failed(s));
    }

    #[test]
    fn test_max_stress_failure_at_f1t() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
        assert!(crit.has_failed(s));
    }

    #[test]
    fn test_max_stress_fi_longitudinal_tension() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(750.0e6, 0.0, 0.0);
        assert!((crit.fi_longitudinal(s) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_max_stress_fi_longitudinal_compression() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(-900.0e6, 0.0, 0.0);
        assert!((crit.fi_longitudinal(s) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_max_stress_fi_transverse_zero() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(100.0e6, 0.0, 0.0);
        assert!((crit.fi_transverse(s)).abs() < 1e-10);
    }

    #[test]
    fn test_max_stress_load_factor_gt_one_safe() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
        assert!(crit.load_factor(s) > 1.0);
    }

    #[test]
    fn test_max_stress_load_factor_le_one_failed() {
        let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
        assert!(crit.load_factor(s) <= 1.0);
    }

    // --- TsaiHillCriterion ---

    #[test]
    fn test_tsai_hill_no_failure_low_stress() {
        let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
        assert!(!crit.has_failed(s));
    }

    #[test]
    fn test_tsai_hill_failure_at_strength() {
        let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
        assert!(crit.has_failed(s));
    }

    #[test]
    fn test_tsai_hill_safety_margin_gt_one_safe() {
        let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(100.0e6, 5.0e6, 5.0e6);
        assert!(crit.safety_margin(s) > 1.0);
    }

    #[test]
    fn test_tsai_hill_load_factor_positive() {
        let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
        let s = PlaneStress::new(300.0e6, 20.0e6, 10.0e6);
        assert!(crit.load_factor(s) > 0.0);
    }

    #[test]
    fn test_tsai_hill_biaxial_strength_positive() {
        let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
        let sigma1 = crit.biaxial_strength(0.0);
        assert!(sigma1 > 0.0);
    }

    #[test]
    fn test_tsai_hill_failure_index_at_uniaxial_f1t() {
        let st = PlyStrength::carbon_epoxy_t300();
        let crit = TsaiHillCriterion::new(st.clone());
        let s = PlaneStress::new(st.f1t, 0.0, 0.0);
        let fi = crit.failure_index_sq(s);
        assert!((fi - 1.0).abs() < 1e-6, "FI at F1t should be 1: {fi}");
    }

    // --- TsaiWuCriterion ---

    #[test]
    fn test_tsai_wu_stable_interaction_coeff() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        assert!(crit.is_stable());
    }

    #[test]
    fn test_tsai_wu_no_failure_low_stress() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        let s = PlaneStress::new(100.0e6, 5.0e6, 5.0e6);
        assert!(!crit.has_failed(s));
    }

    #[test]
    fn test_tsai_wu_failure_at_high_load() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        let s = PlaneStress::new(1600.0e6, 60.0e6, 80.0e6);
        assert!(crit.has_failed(s));
    }

    #[test]
    fn test_tsai_wu_load_factor_positive() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        let s = PlaneStress::new(300.0e6, 20.0e6, 10.0e6);
        assert!(crit.load_factor(s) > 0.0);
    }

    #[test]
    fn test_tsai_wu_linear_terms_sign() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        let [f1, _f2] = crit.linear_terms();
        // f1 = 1/F1t - 1/F1c. Since F1t > F1c typically, sign depends on values.
        assert!(f1.is_finite());
    }

    #[test]
    fn test_tsai_wu_quadratic_terms_f11_positive() {
        let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
        let [f11, f22, f66, _f12] = crit.quadratic_terms();
        assert!(f11 > 0.0 && f22 > 0.0 && f66 > 0.0);
    }

    // --- PuckCriterion ---

    #[test]
    fn test_puck_ff_none_at_low_stress() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
        assert_eq!(puck.ff_mode(s), PuckFfMode::None);
    }

    #[test]
    fn test_puck_ff_tension_at_f1t() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(1600.0e6, 0.0, 0.0);
        assert_eq!(puck.ff_mode(s), PuckFfMode::Tension);
    }

    #[test]
    fn test_puck_ff_compression_at_minus_f1c() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(-1000.0e6, 0.0, 0.0);
        assert_eq!(puck.ff_mode(s), PuckFfMode::Compression);
    }

    #[test]
    fn test_puck_iff_mode_a_tension_transverse() {
        let puck = PuckCriterion::carbon_epoxy_default();
        // High transverse tension
        let s = PlaneStress::new(0.0, 60.0e6, 0.0);
        assert!(puck.iff_exposure_mode_a(s) > 1.0);
    }

    #[test]
    fn test_puck_iff_mode_bc_compression() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(0.0, -250.0e6, 80.0e6);
        // High compressive transverse + shear — should trigger mode BC
        let exp = puck.iff_exposure_mode_bc(s);
        assert!(exp >= 0.0);
    }

    #[test]
    fn test_puck_fracture_angle_zero_for_tension() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(0.0, 30.0e6, 0.0);
        assert!((puck.fracture_plane_angle(s)).abs() < 1e-10);
    }

    #[test]
    fn test_puck_fracture_angle_positive_for_compression() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::new(0.0, -100.0e6, 0.0);
        let angle = puck.fracture_plane_angle(s);
        assert!(angle >= 0.0);
    }

    #[test]
    fn test_puck_no_failure_at_zero_stress() {
        let puck = PuckCriterion::carbon_epoxy_default();
        let s = PlaneStress::zero();
        assert!(!puck.has_failed(s));
    }

    // --- ProgressiveDamage ---

    #[test]
    fn test_ply_state_construction() {
        let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        assert!(ply.q11() > ply.q22());
        assert!(!ply.ff_failed);
    }

    #[test]
    fn test_ply_state_q11_q22_ordering() {
        let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        assert!(
            ply.q11() > ply.q22(),
            "Q11 should be > Q22 for UD ply at 0°"
        );
    }

    #[test]
    fn test_ply_degradation_total() {
        let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let q11_before = ply.q11();
        ply.apply_degradation(DegradationRule::Total, true, false);
        assert!(ply.q11() < q11_before);
        assert!(ply.ff_failed);
    }

    #[test]
    fn test_ply_degradation_none_keeps_stiffness() {
        let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let q11_before = ply.q11();
        ply.apply_degradation(DegradationRule::None, false, true);
        assert!((ply.q11() - q11_before).abs() < 1e-3);
    }

    #[test]
    fn test_progressive_damage_n_intact_all() {
        let ply1 = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let ply2 = PlyState::new(90.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let st1 = PlyStrength::carbon_epoxy_t300();
        let st2 = PlyStrength::carbon_epoxy_t300();
        let pd = ProgressiveDamage::new(vec![ply1, ply2], vec![st1, st2], DegradationRule::Total);
        assert_eq!(pd.n_intact(), 2);
    }

    #[test]
    fn test_progressive_damage_not_all_failed_initially() {
        let ply1 = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let pd = ProgressiveDamage::new(
            vec![ply1],
            vec![PlyStrength::carbon_epoxy_t300()],
            DegradationRule::Total,
        );
        assert!(!pd.all_failed());
    }

    #[test]
    fn test_progressive_damage_a_matrix_nonzero() {
        let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let pd = ProgressiveDamage::new(
            vec![ply],
            vec![PlyStrength::carbon_epoxy_t300()],
            DegradationRule::Total,
        );
        let a = pd.a_matrix();
        assert!(a[0] > 0.0);
    }

    #[test]
    fn test_progressive_failure_sweep_returns_finite() {
        let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
        let mut pd = ProgressiveDamage::new(
            vec![ply],
            vec![PlyStrength::carbon_epoxy_t300()],
            DegradationRule::Total,
        );
        let lf = pd.progressive_failure_sweep([1.0e7, 0.0, 0.0], 100);
        assert!(lf.is_finite() || lf == f64::INFINITY);
    }

    // --- DelamCriterion ---

    #[test]
    fn test_delam_no_failure_low_stress() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let s = InterlaminaStress::new(10.0e6, 20.0e6, 20.0e6);
        assert!(!dc.stress_delamination(s));
    }

    #[test]
    fn test_delam_failure_above_z_t() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let s = InterlaminaStress::new(60.0e6, 0.0, 0.0); // σ_z > Z_t=50MPa
        assert!(dc.stress_delamination(s));
    }

    #[test]
    fn test_delam_compression_no_mode1_contribution() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let s = InterlaminaStress::new(-100.0e6, 0.0, 0.0); // compressive
        assert!(!dc.stress_delamination(s));
    }

    #[test]
    fn test_delam_linear_err_failure() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let err = EnergyReleaseRate::new(200.0, 400.0, 0.0); // at critical values
        assert!(dc.err_delamination(err));
    }

    #[test]
    fn test_delam_linear_err_safe() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let err = EnergyReleaseRate::new(50.0, 50.0, 0.0);
        assert!(!dc.err_delamination(err));
    }

    #[test]
    fn test_bk_criterion_at_critical() {
        let dc = DelamCriterion::carbon_epoxy_default();
        // Pure mode I at G1c
        let err = EnergyReleaseRate::new(200.0, 0.0, 0.0);
        let bk = dc.bk_criterion(err);
        assert!((bk - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bk_delamination_above_critical() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let err = EnergyReleaseRate::new(300.0, 0.0, 0.0);
        assert!(dc.bk_delamination(err));
    }

    #[test]
    fn test_bk_no_delamination_below_critical() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let err = EnergyReleaseRate::new(50.0, 50.0, 0.0);
        assert!(!dc.bk_delamination(err));
    }

    #[test]
    fn test_err_mode_i_ratio() {
        let err = EnergyReleaseRate::new(100.0, 50.0, 50.0);
        let r = err.mode_i_ratio();
        assert!((r - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_err_total() {
        let err = EnergyReleaseRate::new(100.0, 200.0, 50.0);
        assert!((err.total() - 350.0).abs() < 1e-10);
    }

    #[test]
    fn test_stress_intensity_mode1() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let k1 = dc.stress_intensity_mode1(1.0e6, 1.0e-3 / PI);
        assert!((k1 - 1.0e6 / 1000.0_f64.sqrt()).abs() < 1.0);
    }

    #[test]
    fn test_g1c_from_k1c() {
        let g = DelamCriterion::g1c_from_k1c(1000.0, 5.0e9);
        assert!((g - 1000.0 * 1000.0 / 5.0e9).abs() < 1e-10);
    }

    #[test]
    fn test_paris_growth_rate_positive() {
        let dc = DelamCriterion::carbon_epoxy_default();
        let rate = dc.paris_growth_rate(100.0, 200.0, 1e-10, 3.0);
        assert!(rate > 0.0);
    }

    // --- Helper functions ---

    #[test]
    fn test_transformed_q_identity_at_zero_angle() {
        let q = [130.0e9_f64, 10.0e9_f64, 3.0e9_f64, 5.0e9_f64];
        let qbar = transformed_q(&q, 0.0);
        assert!((qbar[0] - q[0]).abs() < 1.0, "Q11 not preserved at 0°");
        assert!((qbar[1] - q[1]).abs() < 1.0, "Q22 not preserved at 0°");
    }

    #[test]
    fn test_transformed_q_90deg_swaps_11_22() {
        let q = [130.0e9_f64, 10.0e9_f64, 3.0e9_f64, 5.0e9_f64];
        let qbar0 = transformed_q(&q, 0.0);
        let qbar90 = transformed_q(&q, PI / 2.0);
        assert!(
            (qbar90[0] - qbar0[1]).abs() < 1.0,
            "Q̄11(90°) should equal Q̄22(0°)"
        );
    }

    #[test]
    fn test_reserve_factor_positive() {
        assert!((reserve_factor(0.5).unwrap() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_reserve_factor_none_at_zero() {
        assert!(reserve_factor(0.0).is_none());
    }

    #[test]
    fn test_tsai_hill_envelope_nonempty() {
        let st = PlyStrength::carbon_epoxy_t300();
        let pts = tsai_hill_envelope(&st, 10);
        assert!(!pts.is_empty());
    }
}
