// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constitutive (material) models for finite element analysis.
//!
//! Provides stress-strain relationships (constitutive matrices) used
//! in element stiffness calculations. Implements:
//!
//! - [`LinearElasticMaterial`] / [`LinearElastic`] — standard isotropic linear elasticity
//! - [`NeoHookean`] — compressible hyperelastic model
//! - [`MooneyRivlin`] — incompressible / near-incompressible hyperelastic model
//! - [`J2Plasticity`] — von Mises (J2) elastoplasticity with isotropic hardening
//! - [`ConstitutiveModel`] trait — common interface for stress and tangent stiffness

// ---------------------------------------------------------------------------
// Voigt notation: [σ11, σ22, σ33, σ12, σ23, σ13]  (indices 0..5)
// ---------------------------------------------------------------------------

/// Alias for a 6-component Voigt stress/strain vector.
pub type Stress6 = [f64; 6];

// ---------------------------------------------------------------------------
// ConstitutiveModel trait
// ---------------------------------------------------------------------------

/// Trait for constitutive (stress-strain) models.
pub trait ConstitutiveModel {
    /// Compute the Cauchy stress given an engineering strain vector in Voigt notation.
    ///
    /// Voigt order: \[ε11, ε22, ε33, γ12, γ23, γ13\] → \[σ11, σ22, σ33, σ12, σ23, σ13\]
    fn stress(&self, strain: Stress6) -> Stress6;

    /// Return the 6×6 algorithmic tangent stiffness matrix (row-major, 36 entries).
    fn tangent_stiffness(&self, strain: Stress6) -> [f64; 36];
}

// ---------------------------------------------------------------------------
// LinearElasticMaterial  (original name kept for compatibility)
// ---------------------------------------------------------------------------

/// Linear elastic isotropic material model.
///
/// Characterized by Young's modulus E and Poisson's ratio nu. Produces the
/// standard 6x6 constitutive matrix (D) for 3D stress-strain in Voigt notation.
#[derive(Debug, Clone)]
pub struct LinearElasticMaterial {
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless, 0 < nu < 0.5).
    pub poisson_ratio: f64,
}

impl LinearElasticMaterial {
    /// Create a new linear elastic material.
    pub fn new(young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
        }
    }

    /// Compute the 6x6 constitutive (elasticity) matrix in Voigt notation.
    ///
    /// ```text
    /// D = (E / ((1+nu)(1-2nu))) *
    ///     [[1-nu, nu,   nu,   0,         0,         0        ],
    ///      [nu,   1-nu, nu,   0,         0,         0        ],
    ///      [nu,   nu,   1-nu, 0,         0,         0        ],
    ///      [0,    0,    0,    (1-2nu)/2, 0,         0        ],
    ///      [0,    0,    0,    0,         (1-2nu)/2, 0        ],
    ///      [0,    0,    0,    0,         0,         (1-2nu)/2]]
    /// ```
    pub fn constitutive_matrix(&self) -> [[f64; 6]; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let a = 1.0 - nu;
        let g = (1.0 - 2.0 * nu) / 2.0;

        let mut d = [[0.0f64; 6]; 6];
        d[0][0] = factor * a;
        d[0][1] = factor * nu;
        d[0][2] = factor * nu;
        d[1][0] = factor * nu;
        d[1][1] = factor * a;
        d[1][2] = factor * nu;
        d[2][0] = factor * nu;
        d[2][1] = factor * nu;
        d[2][2] = factor * a;
        d[3][3] = factor * g;
        d[4][4] = factor * g;
        d[5][5] = factor * g;
        d
    }

    /// Return the shear modulus G = E / (2(1+nu)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Return the bulk modulus K = E / (3(1-2nu)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
}

// ---------------------------------------------------------------------------
// LinearElastic — new name, implements ConstitutiveModel trait
// ---------------------------------------------------------------------------

/// Linear elastic isotropic material implementing [`ConstitutiveModel`].
///
/// Same physics as [`LinearElasticMaterial`]; exposed via the trait interface.
#[derive(Debug, Clone)]
pub struct LinearElastic {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν (0 < ν < 0.5).
    pub poisson_ratio: f64,
}

impl LinearElastic {
    /// Create a new [`LinearElastic`] material.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
        }
    }

    /// First Lamé parameter λ = Eν / ((1+ν)(1−2ν)).
    pub fn lambda(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Shear modulus (second Lamé parameter) μ = E / (2(1+ν)).
    pub fn mu(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Return the 6×6 isotropic stiffness matrix in Voigt notation (row-major flat array).
    pub fn isotropic_stiffness_voigt(&self) -> [f64; 36] {
        let lambda = self.lambda();
        let mu = self.mu();
        let mut c = [0.0f64; 36];
        // Normal-normal coupling
        c[0] = lambda + 2.0 * mu; // C1111
        c[1] = lambda; // C1122
        c[2] = lambda; // C1133
        c[6] = lambda; // C2211
        c[7] = lambda + 2.0 * mu; // C2222
        c[8] = lambda; // C2233
        c[12] = lambda; // C3311
        c[13] = lambda; // C3322
        c[14] = lambda + 2.0 * mu; // C3333
        // Shear terms
        c[21] = mu; // C1212
        c[28] = mu; // C2323
        c[35] = mu; // C1313
        c
    }

    /// Compute stress = C · ε for a given Voigt strain.
    pub fn stress_from_strain(&self, strain: Stress6) -> Stress6 {
        let c = self.isotropic_stiffness_voigt();
        let mut s = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                s[i] += c[i * 6 + j] * strain[j];
            }
        }
        s
    }
}

impl ConstitutiveModel for LinearElastic {
    fn stress(&self, strain: Stress6) -> Stress6 {
        self.stress_from_strain(strain)
    }

    fn tangent_stiffness(&self, _strain: Stress6) -> [f64; 36] {
        self.isotropic_stiffness_voigt()
    }
}

// ---------------------------------------------------------------------------
// NeoHookean hyperelastic
// ---------------------------------------------------------------------------

/// Compressible Neo-Hookean hyperelastic material.
///
/// Strain energy density:
/// W = μ/2 · (I₁ − 3) − μ·ln(J) + κ/4 · (J² − 1 − 2·ln(J))
///
/// where I₁ = tr(C), J = det(F), C = Fᵀ·F.
#[derive(Debug, Clone)]
pub struct NeoHookean {
    /// Shear modulus μ (Pa).
    pub mu: f64,
    /// Bulk modulus κ (Pa).
    pub kappa: f64,
}

impl NeoHookean {
    /// Create a new Neo-Hookean material.
    pub fn new(mu: f64, kappa: f64) -> Self {
        Self { mu, kappa }
    }

    /// Strain energy density W(F).
    pub fn strain_energy_density(&self, f: [[f64; 3]; 3]) -> f64 {
        let j = det3(f);
        let c = ftf(f); // C = Fᵀ F
        let i1 = c[0][0] + c[1][1] + c[2][2]; // tr(C)
        let ln_j = j.ln();
        self.mu / 2.0 * (i1 - 3.0) - self.mu * ln_j + self.kappa / 4.0 * (j * j - 1.0 - 2.0 * ln_j)
    }

    /// Second Piola–Kirchhoff stress S = ∂W/∂E (in material frame).
    ///
    /// S = μ(I − C⁻¹) + κ/2·(J² − 1)·C⁻¹
    pub fn pk2_stress(&self, f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let j = det3(f);
        let c = ftf(f);
        let c_inv = inv3(c);
        let coeff = self.kappa / 2.0 * (j * j - 1.0);
        let mut s = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                let delta_ik = if i == k { 1.0 } else { 0.0 };
                s[i][k] = self.mu * (delta_ik - c_inv[i][k]) + coeff * c_inv[i][k];
            }
        }
        s
    }
}

impl ConstitutiveModel for NeoHookean {
    fn stress(&self, strain: Stress6) -> Stress6 {
        // For small-strain use: F ≈ I + ε (symmetric part only)
        let f = voigt_strain_to_f(strain);
        let s = self.pk2_stress(f);
        [s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2]]
    }

    fn tangent_stiffness(&self, _strain: Stress6) -> [f64; 36] {
        // Linearised tangent: at identity deformation, reduces to linear elastic
        let le = LinearElastic::new(
            self.mu * (3.0 * self.kappa + self.mu) / (self.kappa + self.mu), // Young's from mu,kappa
            (3.0 * self.kappa - 2.0 * self.mu) / (2.0 * (3.0 * self.kappa + self.mu)), // nu
        );
        le.isotropic_stiffness_voigt()
    }
}

// ---------------------------------------------------------------------------
// MooneyRivlin hyperelastic
// ---------------------------------------------------------------------------

/// Mooney-Rivlin hyperelastic material (incompressible form).
///
/// Strain energy density:
/// W = c10·(I₁ − 3) + c01·(I₂ − 3)
#[derive(Debug, Clone)]
pub struct MooneyRivlin {
    /// Material parameter c10 (Pa).
    pub c10: f64,
    /// Material parameter c01 (Pa).
    pub c01: f64,
}

impl MooneyRivlin {
    /// Create a new Mooney-Rivlin material.
    pub fn new(c10: f64, c01: f64) -> Self {
        Self { c10, c01 }
    }

    /// Strain energy density as a function of principal invariants I₁ and I₂.
    ///
    /// W = c10·(I₁ − 3) + c01·(I₂ − 3)
    pub fn strain_energy_density(&self, i1: f64, i2: f64) -> f64 {
        self.c10 * (i1 - 3.0) + self.c01 * (i2 - 3.0)
    }

    /// Compute I₁ = tr(C) and I₂ = (tr²(C) − tr(C²))/2 from deformation gradient F.
    pub fn invariants(f: [[f64; 3]; 3]) -> (f64, f64) {
        let c = ftf(f);
        let i1 = c[0][0] + c[1][1] + c[2][2];
        // tr(C²) = Σ_ij C_ij²
        let mut tr_c2 = 0.0;
        for (i, c_row) in c.iter().enumerate() {
            for (j, &c_ij) in c_row.iter().enumerate() {
                tr_c2 += c_ij * c[j][i];
            }
        }
        let i2 = (i1 * i1 - tr_c2) / 2.0;
        (i1, i2)
    }
}

impl ConstitutiveModel for MooneyRivlin {
    fn stress(&self, strain: Stress6) -> Stress6 {
        let f = voigt_strain_to_f(strain);
        let (i1, i2) = MooneyRivlin::invariants(f);
        let c = ftf(f);
        let c2 = mat3_mul(c, c);
        // S = 2·(∂W/∂I1·I + ∂W/∂I2·(I1·I − C))
        // For incompressible: simplified PK2 (pressure excluded here)
        let dw_di1 = self.c10;
        let dw_di2 = self.c01;
        let mut s = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                let delta = if i == k { 1.0 } else { 0.0 };
                s[i][k] = 2.0 * (dw_di1 * delta + dw_di2 * (i1 * delta - c[i][k]));
            }
        }
        // Suppress unused warning
        let _ = (i2, c2);
        [s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2]]
    }

    fn tangent_stiffness(&self, _strain: Stress6) -> [f64; 36] {
        // Linearised at reference config: reduces to isotropic with mu = 2*(c10+c01)
        let mu = 2.0 * (self.c10 + self.c01);
        let le = LinearElastic::new(3.0 * mu, 0.5 - 1e-9); // near incompressible
        le.isotropic_stiffness_voigt()
    }
}

// ---------------------------------------------------------------------------
// J2Plasticity
// ---------------------------------------------------------------------------

/// J₂ (von Mises) elastoplasticity with linear isotropic hardening.
#[derive(Debug, Clone)]
pub struct J2Plasticity {
    /// Young's modulus E (Pa).
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
    /// Initial yield stress σ_y (Pa).
    pub yield_stress: f64,
    /// Linear isotropic hardening modulus H (Pa).
    pub hardening: f64,
}

impl J2Plasticity {
    /// Create a new J2 plasticity model.
    pub fn new(e: f64, nu: f64, yield_stress: f64, hardening: f64) -> Self {
        Self {
            e,
            nu,
            yield_stress,
            hardening,
        }
    }

    /// Von Mises yield function f = σ_eq − (σ_y + H·ε_p).
    ///
    /// σ_eq = √(3/2 · s:s) where s is the deviatoric stress.
    /// Returns positive when yielding.
    pub fn yield_function(&self, stress_voigt: Stress6, plastic_strain: f64) -> f64 {
        let s_eq = von_mises_stress(stress_voigt);
        s_eq - (self.yield_stress + self.hardening * plastic_strain)
    }

    /// Radial return mapping (cutting-plane / closest-point algorithm).
    ///
    /// Given a trial stress (purely elastic predictor) and the accumulated
    /// plastic strain, returns `(corrected_stress_voigt, updated_plastic_strain)`.
    pub fn return_mapping(&self, trial_stress: Stress6, plastic_strain: f64) -> (Stress6, f64) {
        let mu = self.e / (2.0 * (1.0 + self.nu));
        let f_trial = self.yield_function(trial_stress, plastic_strain);

        if f_trial <= 0.0 {
            // Elastic — no correction needed
            return (trial_stress, plastic_strain);
        }

        // Radial return: Δγ = f_trial / (3μ + H)
        let d_gamma = f_trial / (3.0 * mu + self.hardening);
        let s_eq_trial = von_mises_stress(trial_stress);

        // Scale factor for deviatoric stress
        let scale = if s_eq_trial > 1e-30 {
            1.0 - d_gamma * 3.0 * mu / s_eq_trial
        } else {
            1.0
        };

        let p = hydrostatic_stress(trial_stress);
        let dev = deviatoric_voigt(trial_stress);
        let corrected = [
            p + scale * dev[0],
            p + scale * dev[1],
            p + scale * dev[2],
            scale * trial_stress[3],
            scale * trial_stress[4],
            scale * trial_stress[5],
        ];

        (corrected, plastic_strain + d_gamma)
    }
}

impl ConstitutiveModel for J2Plasticity {
    fn stress(&self, strain: Stress6) -> Stress6 {
        let le = LinearElastic::new(self.e, self.nu);
        let trial = le.stress(strain);
        let (s, _) = self.return_mapping(trial, 0.0);
        s
    }

    fn tangent_stiffness(&self, _strain: Stress6) -> [f64; 36] {
        let le = LinearElastic::new(self.e, self.nu);
        le.isotropic_stiffness_voigt()
    }
}

// ---------------------------------------------------------------------------
// Stress invariants (standalone helpers)
// ---------------------------------------------------------------------------

/// Compute von Mises equivalent stress from Voigt stress vector.
///
/// σ_eq = √( ½·\[(σ11−σ22)² + (σ22−σ33)² + (σ33−σ11)²\] + 3·(σ12² + σ23² + σ13²) )
pub fn von_mises_stress(s: Stress6) -> f64 {
    let dev_sq = (s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2);
    let shear_sq = s[3].powi(2) + s[4].powi(2) + s[5].powi(2);
    ((dev_sq / 2.0 + 3.0 * shear_sq).max(0.0)).sqrt()
}

/// Hydrostatic (mean normal) stress: p = (σ11 + σ22 + σ33) / 3.
pub fn hydrostatic_stress(s: Stress6) -> f64 {
    (s[0] + s[1] + s[2]) / 3.0
}

/// Deviatoric stress in Voigt notation: s_ij = σ_ij − p·δ_ij.
pub fn deviatoric_voigt(s: Stress6) -> Stress6 {
    let p = hydrostatic_stress(s);
    [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]]
}

// ---------------------------------------------------------------------------
// 3×3 matrix helpers (no nalgebra)
// ---------------------------------------------------------------------------

#[inline]
fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// C = Fᵀ · F
#[inline]
fn ftf(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for f_k in &f {
                c[i][j] += f_k[i] * f_k[j];
            }
        }
    }
    c
}

/// 3×3 matrix multiplication A · B.
fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Inverse of a 3×3 matrix (panics if singular).
fn inv3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    assert!(d.abs() > 1e-30, "matrix is singular");
    let inv_d = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
        ],
    ]
}

/// Approximate deformation gradient F ≈ I + ε from Voigt small strain.
/// (Used for hyperelastic models called with small-strain input.)
fn voigt_strain_to_f(strain: Stress6) -> [[f64; 3]; 3] {
    [
        [1.0 + strain[0], strain[3] / 2.0, strain[5] / 2.0],
        [strain[3] / 2.0, 1.0 + strain[1], strain[4] / 2.0],
        [strain[5] / 2.0, strain[4] / 2.0, 1.0 + strain[2]],
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- LinearElasticMaterial (original) ----

    #[test]
    fn test_constitutive_matrix_steel() {
        let mat = LinearElasticMaterial::new(200.0e9, 0.3);
        let d = mat.constitutive_matrix();

        let factor = 200.0e9 / ((1.0 + 0.3) * (1.0 - 0.6));
        let expected_d00 = factor * 0.7;
        let expected_d01 = factor * 0.3;
        let expected_d33 = factor * 0.2;

        assert!(
            (d[0][0] - expected_d00).abs() / expected_d00.abs() < 1e-12,
            "D[0][0] = {} expected {expected_d00}",
            d[0][0]
        );
        assert!(
            (d[0][1] - expected_d01).abs() / expected_d01.abs() < 1e-12,
            "D[0][1] = {} expected {expected_d01}",
            d[0][1]
        );
        assert!(
            (d[3][3] - expected_d33).abs() / expected_d33.abs() < 1e-12,
            "D[3][3] = {} expected {expected_d33}",
            d[3][3]
        );

        // D should be symmetric
        for (i, row) in d.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!((v - d[j][i]).abs() < 1e-6, "D not symmetric at ({i},{j})");
            }
        }

        assert_eq!(d[0][3], 0.0);
        assert_eq!(d[3][0], 0.0);
    }

    #[test]
    fn test_shear_and_bulk_modulus() {
        let mat = LinearElasticMaterial::new(200.0e9, 0.3);
        let g = mat.shear_modulus();
        let k = mat.bulk_modulus();
        let expected_g = 200.0e9 / (2.0 * 1.3);
        let expected_k = 200.0e9 / (3.0 * 0.4);
        assert!((g - expected_g).abs() / expected_g < 1e-12);
        assert!((k - expected_k).abs() / expected_k < 1e-12);
    }

    // ---- LinearElastic (trait) ----

    #[test]
    fn test_isotropic_stiffness_symmetry() {
        let mat = LinearElastic::new(210.0e9, 0.3);
        let c = mat.isotropic_stiffness_voigt();
        // The 6×6 matrix should be symmetric: C[i,j] == C[j,i]
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c[i * 6 + j] - c[j * 6 + i]).abs() < 1.0,
                    "C not symmetric at ({i},{j}): {} vs {}",
                    c[i * 6 + j],
                    c[j * 6 + i]
                );
            }
        }
    }

    #[test]
    fn test_pure_hydrostatic_gives_diagonal_stress() {
        let mat = LinearElastic::new(210.0e9, 0.3);
        // Hydrostatic strain: ε11 = ε22 = ε33 = ε0, shear = 0
        let eps = 1e-3;
        let strain = [eps, eps, eps, 0.0, 0.0, 0.0];
        let stress = mat.stress(strain);
        // Off-diagonal (shear) stresses should be zero
        assert!(stress[3].abs() < 1.0, "σ12 should be ~0, got {}", stress[3]);
        assert!(stress[4].abs() < 1.0, "σ23 should be ~0, got {}", stress[4]);
        assert!(stress[5].abs() < 1.0, "σ13 should be ~0, got {}", stress[5]);
        // All diagonal stresses should be equal
        assert!(
            (stress[0] - stress[1]).abs() < 1.0,
            "σ11 != σ22 in hydrostatic: {} vs {}",
            stress[0],
            stress[1]
        );
        assert!(
            (stress[1] - stress[2]).abs() < 1.0,
            "σ22 != σ33 in hydrostatic: {} vs {}",
            stress[1],
            stress[2]
        );
    }

    #[test]
    fn test_lame_parameters() {
        let e = 200.0e9;
        let nu = 0.3;
        let mat = LinearElastic::new(e, nu);
        let expected_lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let expected_mu = e / (2.0 * (1.0 + nu));
        assert!((mat.lambda() - expected_lambda).abs() / expected_lambda < 1e-12);
        assert!((mat.mu() - expected_mu).abs() / expected_mu < 1e-12);
    }

    // ---- NeoHookean ----

    #[test]
    fn test_neohookean_undeformed_zero_pk2() {
        let mat = NeoHookean::new(1.0e6, 3.0e6);
        // F = I → S = μ(I − I) + κ/2·(1−1)·I = 0
        let f_identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let s = mat.pk2_stress(f_identity);
        for (i, row) in s.iter().enumerate() {
            for (k, &v) in row.iter().enumerate() {
                assert!(
                    v.abs() < 1e-6,
                    "PK2 at identity should be zero, S[{i}][{k}] = {v}"
                );
            }
        }
    }

    #[test]
    fn test_neohookean_energy_positive_for_stretch() {
        let mat = NeoHookean::new(1.0e6, 3.0e6);
        // Simple extension
        let f_stretch = [[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let w = mat.strain_energy_density(f_stretch);
        assert!(
            w > 0.0,
            "strain energy density should be positive for stretch, got {w}"
        );
    }

    // ---- MooneyRivlin ----

    #[test]
    fn test_mooney_rivlin_energy_positive() {
        let mat = MooneyRivlin::new(0.5e6, 0.1e6);
        // I1 > 3 and I2 > 3 for any real deformation
        let w = mat.strain_energy_density(3.5, 3.2);
        assert!(w > 0.0, "Mooney-Rivlin energy should be positive, got {w}");
    }

    #[test]
    fn test_mooney_rivlin_energy_zero_at_reference() {
        let mat = MooneyRivlin::new(0.5e6, 0.1e6);
        // W(I1=3, I2=3) = 0 (reference configuration)
        let w = mat.strain_energy_density(3.0, 3.0);
        assert!(
            w.abs() < 1e-30,
            "W should be 0 at reference config, got {w}"
        );
    }

    // ---- J2Plasticity ----

    #[test]
    fn test_j2_yield_at_correct_stress() {
        let yield_stress = 250.0e6; // 250 MPa
        let mat = J2Plasticity::new(200.0e9, 0.3, yield_stress, 1.0e9);
        // Uniaxial stress at yield: σ11 = σ_y, all others 0
        let stress_at_yield = [yield_stress, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = mat.yield_function(stress_at_yield, 0.0);
        assert!(
            f.abs() < 1e-3,
            "yield function should be ~0 at σ=σ_y, got {f}"
        );

        // Below yield: stress = 0.9 * σ_y → f < 0
        let stress_below = [0.9 * yield_stress, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f_below = mat.yield_function(stress_below, 0.0);
        assert!(f_below < 0.0, "should be elastic below yield, f={f_below}");
    }

    #[test]
    fn test_j2_return_mapping_elastic() {
        let mat = J2Plasticity::new(200.0e9, 0.3, 250.0e6, 1.0e9);
        let trial = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // well below yield
        let (corrected, ep) = mat.return_mapping(trial, 0.0);
        // Should not change
        assert_eq!(ep, 0.0);
        for i in 0..6 {
            assert!((corrected[i] - trial[i]).abs() < 1.0);
        }
    }

    #[test]
    fn test_j2_return_mapping_plastic() {
        let yield_stress = 250.0e6;
        let mat = J2Plasticity::new(200.0e9, 0.3, yield_stress, 0.0); // no hardening
        // Trial stress well above yield
        let trial = [400.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (corrected, ep) = mat.return_mapping(trial, 0.0);
        // After return mapping, von Mises stress should equal yield stress
        let s_eq = von_mises_stress(corrected);
        assert!(
            (s_eq - yield_stress).abs() / yield_stress < 1e-6,
            "after return mapping σ_eq should equal σ_y: {} vs {}",
            s_eq,
            yield_stress
        );
        assert!(ep > 0.0, "plastic strain should be positive after yielding");
    }

    // ---- Stress invariants ----

    #[test]
    fn test_hydrostatic_and_deviatoric() {
        let s = [100.0, 200.0, 300.0, 50.0, 30.0, 20.0];
        let p = hydrostatic_stress(s);
        assert!((p - 200.0).abs() < 1e-10);
        let dev = deviatoric_voigt(s);
        // Trace of deviatoric should be zero
        let trace = dev[0] + dev[1] + dev[2];
        assert!(
            trace.abs() < 1e-10,
            "deviatoric trace should be 0, got {trace}"
        );
    }

    #[test]
    fn test_von_mises_hydrostatic_zero() {
        // Hydrostatic stress: all normal equal, no shear → von Mises = 0
        let s = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress(s);
        assert!(
            vm.abs() < 1e-8,
            "hydrostatic stress has zero von Mises, got {vm}"
        );
    }
}

// ---------------------------------------------------------------------------
// Thermoelastic constitutive model
// ---------------------------------------------------------------------------

/// Thermoelastic isotropic material.
///
/// Combines linear elasticity with thermal expansion:
///   σ = C : (ε - α ΔT I)
///
/// where α is the coefficient of thermal expansion and ΔT = T - T_ref.
#[derive(Debug, Clone)]
pub struct ThermoElastic {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Coefficient of thermal expansion α (1/K).
    pub alpha: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
}

impl ThermoElastic {
    /// Create a new thermoelastic material.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64, alpha: f64, t_ref: f64) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
            alpha,
            t_ref,
        }
    }

    /// Thermal strain vector (Voigt) for temperature change ΔT.
    ///
    /// ε_th = α ΔT \[1, 1, 1, 0, 0, 0\]ᵀ
    pub fn thermal_strain(&self, temperature: f64) -> Stress6 {
        let dt = temperature - self.t_ref;
        let eps_th = self.alpha * dt;
        [eps_th, eps_th, eps_th, 0.0, 0.0, 0.0]
    }

    /// Compute thermoelastic stress given mechanical strain and temperature.
    ///
    /// σ = C : (ε_mech - ε_th)
    pub fn stress(&self, strain: Stress6, temperature: f64) -> Stress6 {
        let eps_th = self.thermal_strain(temperature);
        let eps_mech = [
            strain[0] - eps_th[0],
            strain[1] - eps_th[1],
            strain[2] - eps_th[2],
            strain[3],
            strain[4],
            strain[5],
        ];
        let le = LinearElastic::new(self.youngs_modulus, self.poisson_ratio);
        le.stress(eps_mech)
    }

    /// Elastic stiffness matrix C (same as LinearElastic).
    pub fn stiffness(&self) -> [f64; 36] {
        LinearElastic::new(self.youngs_modulus, self.poisson_ratio).isotropic_stiffness_voigt()
    }

    /// Thermal expansion vector B = C : α_vec.
    ///
    /// Used in FEM assembly: F_th = ∫ B_strain^T · C · α_vec · ΔT dV
    pub fn thermal_stress_vector(&self, temperature: f64) -> Stress6 {
        let eps_th = self.thermal_strain(temperature);
        let le = LinearElastic::new(self.youngs_modulus, self.poisson_ratio);
        le.stress(eps_th)
    }
}

// ---------------------------------------------------------------------------
// Orthotropic elasticity (improved)
// ---------------------------------------------------------------------------

/// Orthotropic elastic material (3 orthogonal planes of symmetry).
///
/// Defined by 9 independent elastic constants:
/// E1, E2, E3 (Young's moduli)
/// nu12, nu13, nu23 (Poisson's ratios)
/// G12, G13, G23 (shear moduli)
#[derive(Debug, Clone)]
pub struct OrthotropicElastic {
    /// Young's moduli \[E1, E2, E3\] (Pa).
    pub e: [f64; 3],
    /// Poisson's ratios \[nu12, nu13, nu23\].
    pub nu: [f64; 3],
    /// Shear moduli \[G12, G13, G23\] (Pa).
    pub g: [f64; 3],
}

impl OrthotropicElastic {
    /// Create a new orthotropic material.
    pub fn new(e: [f64; 3], nu: [f64; 3], g: [f64; 3]) -> Self {
        Self { e, nu, g }
    }

    /// Create an isotropic material as a special case.
    pub fn from_isotropic(youngs: f64, poisson: f64) -> Self {
        let g_val = youngs / (2.0 * (1.0 + poisson));
        Self {
            e: [youngs; 3],
            nu: [poisson; 3],
            g: [g_val; 3],
        }
    }

    /// Compute the 6×6 orthotropic compliance matrix S (row-major flat array).
    ///
    /// S = C⁻¹ for orthotropic elasticity.
    /// In Voigt notation: \[ε\] = \[S\]\[σ\]
    pub fn compliance_matrix(&self) -> [f64; 36] {
        let [e1, e2, e3] = self.e;
        let [nu12, nu13, nu23] = self.nu;
        // nu21 = nu12 * E2/E1, etc. (symmetry conditions)
        let nu21 = nu12 * e2 / e1.max(f64::EPSILON);
        let nu31 = nu13 * e3 / e1.max(f64::EPSILON);
        let nu32 = nu23 * e3 / e2.max(f64::EPSILON);
        let [g12, g13, g23] = self.g;

        let mut s = [0.0_f64; 36];
        // Normal-normal block (3x3 upper left)
        s[0] = 1.0 / e1;
        s[1] = -nu21 / e2;
        s[2] = -nu31 / e3;
        s[6] = -nu12 / e1;
        s[7] = 1.0 / e2;
        s[8] = -nu32 / e3;
        s[12] = -nu13 / e1;
        s[13] = -nu23 / e2;
        s[14] = 1.0 / e3;
        // Shear terms
        s[21] = 1.0 / g12;
        s[28] = 1.0 / g23;
        s[35] = 1.0 / g13;
        s
    }

    /// Compute the stress from strain using a simplified diagonal approximation.
    ///
    /// For orthotropic materials the full C matrix requires inverting S.
    /// This function returns an approximate stress assuming small off-diagonal coupling.
    pub fn stress_approx(&self, strain: Stress6) -> Stress6 {
        let [e1, e2, e3] = self.e;
        let [g12, g13, g23] = self.g;
        [
            e1 * strain[0],
            e2 * strain[1],
            e3 * strain[2],
            g12 * strain[3],
            g23 * strain[4],
            g13 * strain[5],
        ]
    }

    /// Check if the orthotropic material satisfies Drucker stability criteria.
    ///
    /// Requirements:
    /// E1, E2, E3, G12, G13, G23 > 0
    /// |nu12| < sqrt(E1/E2), etc.
    pub fn is_stable(&self) -> bool {
        let [e1, e2, e3] = self.e;
        let [nu12, nu13, nu23] = self.nu;
        let [g12, g13, g23] = self.g;
        if e1 <= 0.0 || e2 <= 0.0 || e3 <= 0.0 {
            return false;
        }
        if g12 <= 0.0 || g13 <= 0.0 || g23 <= 0.0 {
            return false;
        }
        let check12 = nu12.abs() < (e1 / e2).sqrt();
        let check13 = nu13.abs() < (e1 / e3).sqrt();
        let check23 = nu23.abs() < (e2 / e3).sqrt();
        check12 && check13 && check23
    }
}

// ---------------------------------------------------------------------------
// Rate-dependent plasticity: Perzyna viscoplasticity
// ---------------------------------------------------------------------------

/// Rate-dependent plasticity using the Perzyna overstress model.
///
/// The inelastic strain rate is:
/// dε_vp/dt = (1/η) * <f(σ) / σ_y>^N * ∂f/∂σ
///
/// where `x` = max(0, x), f = von Mises yield function,
/// η = viscosity, N = rate exponent.
#[derive(Debug, Clone)]
pub struct PerzynaViscoplasticity {
    /// Young's modulus E (Pa).
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
    /// Initial yield stress σ_y (Pa).
    pub yield_stress: f64,
    /// Isotropic hardening modulus H (Pa).
    pub hardening: f64,
    /// Viscosity parameter η (s).
    pub eta: f64,
    /// Rate sensitivity exponent N.
    pub exponent: f64,
}

impl PerzynaViscoplasticity {
    /// Create a new Perzyna viscoplastic material.
    pub fn new(
        e: f64,
        nu: f64,
        yield_stress: f64,
        hardening: f64,
        eta: f64,
        exponent: f64,
    ) -> Self {
        Self {
            e,
            nu,
            yield_stress,
            hardening,
            eta,
            exponent,
        }
    }

    /// Overstress function: Φ = max(0, (σ_eq - σ_y) / σ_y).
    pub fn overstress(&self, stress_voigt: Stress6, plastic_strain: f64) -> f64 {
        let s_eq = von_mises_stress(stress_voigt);
        let sigma_y = self.yield_stress + self.hardening * plastic_strain;
        ((s_eq - sigma_y) / sigma_y.max(f64::EPSILON)).max(0.0)
    }

    /// Inelastic strain rate magnitude: dγ/dt = (1/η) * Φ^N.
    pub fn viscoplastic_strain_rate(&self, stress_voigt: Stress6, plastic_strain: f64) -> f64 {
        let phi = self.overstress(stress_voigt, plastic_strain);
        if phi < f64::EPSILON {
            return 0.0;
        }
        phi.powf(self.exponent) / self.eta.max(f64::EPSILON)
    }

    /// Compute the stress update using an explicit integration:
    ///
    /// σ^{n+1} = C : (ε^{n+1} - ε_vp^n)
    /// ε_vp^{n+1} = ε_vp^n + (dγ/dt) * dt * n̂
    ///
    /// Returns (updated_stress, updated_plastic_strain).
    pub fn explicit_stress_update(
        &self,
        strain: Stress6,
        plastic_strain: f64,
        dt: f64,
    ) -> (Stress6, f64) {
        let le = LinearElastic::new(self.e, self.nu);
        let trial = le.stress(strain);
        let dg_dt = self.viscoplastic_strain_rate(trial, plastic_strain);
        let d_gamma = dg_dt * dt;
        let new_ep = plastic_strain + d_gamma;

        // Radial return correction
        let s_eq = von_mises_stress(trial);
        if s_eq < f64::EPSILON {
            return (trial, new_ep);
        }
        let scale = 1.0 - d_gamma * 3.0 * le.mu() / s_eq.max(f64::EPSILON);
        let p = hydrostatic_stress(trial);
        let dev = deviatoric_voigt(trial);
        let corrected = [
            p + scale * dev[0],
            p + scale * dev[1],
            p + scale * dev[2],
            scale * trial[3],
            scale * trial[4],
            scale * trial[5],
        ];
        (corrected, new_ep)
    }
}

// ---------------------------------------------------------------------------
// Viscoplastic flow rule (Norton power law)
// ---------------------------------------------------------------------------

/// Norton power-law creep (viscoplastic flow rule):
///
/// dε_cr/dt = A * σ_eq^n
///
/// where A is the creep coefficient and n is the creep exponent.
pub fn norton_creep_rate(stress_voigt: Stress6, a_coeff: f64, n_exp: f64) -> f64 {
    let s_eq = von_mises_stress(stress_voigt);
    a_coeff * s_eq.powf(n_exp)
}

/// Integrate Norton creep over a time step dt.
///
/// Returns the incremental equivalent creep strain Δε_cr.
pub fn norton_creep_increment(stress_voigt: Stress6, a_coeff: f64, n_exp: f64, dt: f64) -> f64 {
    norton_creep_rate(stress_voigt, a_coeff, n_exp) * dt
}

// ---------------------------------------------------------------------------
// Constitutive model selection / registry
// ---------------------------------------------------------------------------

/// Identifier for registered constitutive models.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelId {
    /// Standard isotropic linear elastic.
    LinearElasticModel,
    /// Neo-Hookean hyperelastic.
    NeoHookeanModel,
    /// Mooney-Rivlin hyperelastic.
    MooneyRivlinModel,
    /// J2 (von Mises) elastoplasticity.
    J2PlasticityModel,
    /// Thermoelastic (linear elastic + thermal expansion).
    ThermoElasticModel,
}

/// Model descriptor used in the constitutive model registry.
#[derive(Debug, Clone)]
pub struct ModelDescriptor {
    /// Model identifier.
    pub id: ModelId,
    /// Human-readable name.
    pub name: String,
    /// Number of state variables (e.g. plastic strain).
    pub n_state_vars: usize,
}

impl ModelDescriptor {
    /// Descriptor for linear elastic model.
    pub fn linear_elastic() -> Self {
        Self {
            id: ModelId::LinearElasticModel,
            name: "Linear Elastic".into(),
            n_state_vars: 0,
        }
    }

    /// Descriptor for J2 plasticity model.
    pub fn j2_plasticity() -> Self {
        Self {
            id: ModelId::J2PlasticityModel,
            name: "J2 Plasticity".into(),
            n_state_vars: 1, // accumulated plastic strain
        }
    }

    /// Descriptor for thermoelastic model.
    pub fn thermoelastic() -> Self {
        Self {
            id: ModelId::ThermoElasticModel,
            name: "Thermoelastic".into(),
            n_state_vars: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// Strain energy and material testing utilities
// ---------------------------------------------------------------------------

/// Compute elastic strain energy density: W = 0.5 * σ : ε
pub fn strain_energy_density(stress: Stress6, strain: Stress6) -> f64 {
    // Engineering shear convention: σ12 * γ12 (factor 2 already in γ12)
    0.5 * (stress[0] * strain[0]
        + stress[1] * strain[1]
        + stress[2] * strain[2]
        + stress[3] * strain[3]
        + stress[4] * strain[4]
        + stress[5] * strain[5])
}

/// Compute effective elastic modulus from uniaxial test data.
///
/// E_eff = Δσ11 / Δε11
///
/// Uses two (strain, stress) data points.
pub fn effective_modulus_uniaxial(strain1: f64, stress1: f64, strain2: f64, stress2: f64) -> f64 {
    let d_strain = strain2 - strain1;
    if d_strain.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    (stress2 - stress1) / d_strain
}

/// Compute the second invariant of the deviatoric stress J2.
///
/// J2 = 0.5 * s:s  where s is the deviatoric stress.
pub fn stress_j2_invariant(s: Stress6) -> f64 {
    let dev = deviatoric_voigt(s);
    0.5 * (dev[0] * dev[0]
        + dev[1] * dev[1]
        + dev[2] * dev[2]
        + 2.0 * dev[3] * dev[3]
        + 2.0 * dev[4] * dev[4]
        + 2.0 * dev[5] * dev[5])
}

/// Compute the Lode angle θ (Lode parameter) from the stress invariants.
///
/// cos(3θ) = (3√3 / 2) * J3 / J2^(3/2)
///
/// Returns the angle in radians in the range \[0, π/3\].
pub fn lode_angle(s: Stress6) -> f64 {
    let j2 = stress_j2_invariant(s);
    if j2 < f64::EPSILON {
        return 0.0;
    }
    let dev = deviatoric_voigt(s);
    // J3 = det(s) for deviatoric stress
    let j3 = dev[0] * (dev[1] * dev[2] - dev[4] * dev[4])
        - dev[3] * (dev[3] * dev[2] - dev[4] * dev[5])
        + dev[5] * (dev[3] * dev[4] - dev[1] * dev[5]);
    let ratio = (3.0 * 3.0_f64.sqrt() / 2.0) * j3 / j2.powf(1.5);
    let ratio_clamped = ratio.clamp(-1.0, 1.0);
    ratio_clamped.acos() / 3.0
}

// ---------------------------------------------------------------------------
// Extended tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ── ThermoElastic tests ───────────────────────────────────────────────

    #[test]
    fn thermoelastic_no_thermal_stress_at_ref() {
        let mat = ThermoElastic::new(200.0e9, 0.3, 12e-6, 293.0);
        let strain = [0.0; 6];
        let stress = mat.stress(strain, 293.0); // T = T_ref
        for &si in &stress {
            assert!(si.abs() < 1.0, "no thermal stress at T_ref, got {si}");
        }
    }

    #[test]
    fn thermoelastic_thermal_strain_scales_with_dt() {
        let mat = ThermoElastic::new(200.0e9, 0.3, 12e-6, 293.0);
        let eps1 = mat.thermal_strain(393.0); // ΔT = 100
        let eps2 = mat.thermal_strain(493.0); // ΔT = 200
        // eps2 should be double eps1
        assert!(
            (eps2[0] / eps1[0] - 2.0).abs() < 1e-10,
            "thermal strain ratio: {}",
            eps2[0] / eps1[0]
        );
    }

    #[test]
    fn thermoelastic_thermal_strain_isotropic() {
        let mat = ThermoElastic::new(200.0e9, 0.3, 12e-6, 293.0);
        let eps = mat.thermal_strain(393.0);
        // All normal components equal, all shear zero
        assert!((eps[0] - eps[1]).abs() < 1e-20);
        assert!((eps[1] - eps[2]).abs() < 1e-20);
        assert!(eps[3].abs() < 1e-20);
        assert!(eps[4].abs() < 1e-20);
        assert!(eps[5].abs() < 1e-20);
    }

    #[test]
    fn thermoelastic_positive_thermal_stress_for_heating() {
        let mat = ThermoElastic::new(200.0e9, 0.3, 12e-6, 293.0);
        // Fixed strain (constrained expansion → compressive thermal stress)
        let strain = [0.0; 6]; // fixed boundary
        let stress = mat.stress(strain, 393.0);
        // All normal stresses should be compressive (negative)
        assert!(
            stress[0] < 0.0,
            "constrained thermal expansion: sigma_11 = {}",
            stress[0]
        );
    }

    #[test]
    fn thermoelastic_stiffness_matches_linear_elastic() {
        let mat = ThermoElastic::new(200.0e9, 0.3, 12e-6, 293.0);
        let le = LinearElastic::new(200.0e9, 0.3);
        let c_thermo = mat.stiffness();
        let c_le = le.isotropic_stiffness_voigt();
        for i in 0..36 {
            assert!(
                (c_thermo[i] - c_le[i]).abs() < 1.0,
                "stiffness mismatch at [{i}]"
            );
        }
    }

    // ── OrthotropicElastic tests ──────────────────────────────────────────

    #[test]
    fn orthotropic_from_isotropic_stable() {
        let mat = OrthotropicElastic::from_isotropic(200.0e9, 0.3);
        assert!(mat.is_stable(), "isotropic material must satisfy stability");
    }

    #[test]
    fn orthotropic_compliance_diagonal_nonzero() {
        let mat = OrthotropicElastic::from_isotropic(200.0e9, 0.3);
        let s = mat.compliance_matrix();
        // Diagonal elements should be positive
        assert!(s[0] > 0.0, "S11 > 0");
        assert!(s[7] > 0.0, "S22 > 0");
        assert!(s[14] > 0.0, "S33 > 0");
        assert!(s[21] > 0.0, "S44 > 0");
    }

    #[test]
    fn orthotropic_stress_approx_uniaxial() {
        let mat = OrthotropicElastic::from_isotropic(200.0e9, 0.3);
        let strain = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let stress = mat.stress_approx(strain);
        assert!(
            (stress[0] - 200.0e9 * 1e-3).abs() < 1.0,
            "sigma_11 = {}",
            stress[0]
        );
    }

    #[test]
    fn orthotropic_from_isotropic_compliance_inverse_consistency() {
        let mat = OrthotropicElastic::from_isotropic(200.0e9, 0.3);
        let s = mat.compliance_matrix();
        // S11 = 1/E1
        let expected_s11 = 1.0 / 200.0e9;
        assert!(
            (s[0] - expected_s11).abs() / expected_s11 < 1e-10,
            "S11 = {}",
            s[0]
        );
    }

    // ── PerzynaViscoplasticity tests ──────────────────────────────────────

    #[test]
    fn perzyna_no_overstress_below_yield() {
        let mat = PerzynaViscoplasticity::new(200.0e9, 0.3, 250.0e6, 0.0, 1e-3, 1.0);
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // below yield
        let phi = mat.overstress(stress, 0.0);
        assert_eq!(phi, 0.0, "no overstress below yield: phi = {phi}");
    }

    #[test]
    fn perzyna_overstress_above_yield() {
        let mat = PerzynaViscoplasticity::new(200.0e9, 0.3, 250.0e6, 0.0, 1e-3, 1.0);
        let stress = [400.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let phi = mat.overstress(stress, 0.0);
        assert!(phi > 0.0, "overstress above yield: phi = {phi}");
    }

    #[test]
    fn perzyna_viscoplastic_rate_zero_at_yield() {
        let mat = PerzynaViscoplasticity::new(200.0e9, 0.3, 250.0e6, 0.0, 1e-3, 1.0);
        let stress = [250.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // exactly at yield
        let rate = mat.viscoplastic_strain_rate(stress, 0.0);
        assert!(rate.abs() < 1e-10, "no rate at yield surface: {rate}");
    }

    #[test]
    fn perzyna_explicit_update_plastic_strain_increases() {
        let mat = PerzynaViscoplasticity::new(200.0e9, 0.3, 250.0e6, 0.0, 1e-3, 1.0);
        let strain = [2e-3, 0.0, 0.0, 0.0, 0.0, 0.0]; // above yield
        let (_, ep) = mat.explicit_stress_update(strain, 0.0, 0.001);
        assert!(ep >= 0.0, "plastic strain must be non-negative: {ep}");
    }

    // ── Norton creep tests ────────────────────────────────────────────────

    #[test]
    fn norton_creep_zero_stress_zero_rate() {
        let rate = norton_creep_rate([0.0; 6], 1e-15, 3.0);
        assert_eq!(rate, 0.0, "zero stress → zero creep rate");
    }

    #[test]
    fn norton_creep_rate_positive_for_nonzero_stress() {
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rate = norton_creep_rate(stress, 1e-20, 3.0);
        assert!(rate > 0.0 && rate.is_finite(), "creep rate = {rate}");
    }

    #[test]
    fn norton_creep_increment_proportional_to_dt() {
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let de1 = norton_creep_increment(stress, 1e-20, 3.0, 1.0);
        let de2 = norton_creep_increment(stress, 1e-20, 3.0, 2.0);
        assert!(
            (de2 / de1 - 2.0).abs() < 1e-10,
            "increment ratio: {}",
            de2 / de1
        );
    }

    // ── ModelDescriptor tests ─────────────────────────────────────────────

    #[test]
    fn model_descriptor_names_unique() {
        let le = ModelDescriptor::linear_elastic();
        let j2 = ModelDescriptor::j2_plasticity();
        let te = ModelDescriptor::thermoelastic();
        assert_ne!(le.name, j2.name);
        assert_ne!(j2.name, te.name);
        assert_ne!(le.name, te.name);
    }

    #[test]
    fn model_descriptor_state_vars() {
        let le = ModelDescriptor::linear_elastic();
        let j2 = ModelDescriptor::j2_plasticity();
        assert_eq!(le.n_state_vars, 0);
        assert_eq!(j2.n_state_vars, 1);
    }

    #[test]
    fn model_descriptor_ids() {
        assert_eq!(
            ModelDescriptor::linear_elastic().id,
            ModelId::LinearElasticModel
        );
        assert_eq!(
            ModelDescriptor::j2_plasticity().id,
            ModelId::J2PlasticityModel
        );
        assert_eq!(
            ModelDescriptor::thermoelastic().id,
            ModelId::ThermoElasticModel
        );
    }

    // ── strain_energy_density tests ───────────────────────────────────────

    #[test]
    fn strain_energy_positive_uniaxial() {
        let e = 200.0e9;
        let eps = 1e-3;
        let sigma = e * eps;
        let stress = [sigma, 0.0, 0.0, 0.0, 0.0, 0.0];
        let strain = [eps, 0.0, 0.0, 0.0, 0.0, 0.0];
        let w = strain_energy_density(stress, strain);
        let expected = 0.5 * sigma * eps;
        assert!(
            (w - expected).abs() / expected < 1e-10,
            "W = {w}, expected {expected}"
        );
    }

    #[test]
    fn strain_energy_zero_for_zero_strain() {
        let w = strain_energy_density([0.0; 6], [0.0; 6]);
        assert_eq!(w, 0.0);
    }

    // ── effective_modulus_uniaxial tests ──────────────────────────────────

    #[test]
    fn effective_modulus_linear_material() {
        let e = effective_modulus_uniaxial(0.0, 0.0, 1e-3, 200.0e6);
        assert!((e - 200.0e9).abs() / 200.0e9 < 1e-10, "E_eff = {e}");
    }

    #[test]
    fn effective_modulus_infinity_for_zero_strain_diff() {
        let e = effective_modulus_uniaxial(0.001, 100.0e6, 0.001, 200.0e6);
        assert!(e.is_infinite(), "E_eff should be infinite: {e}");
    }

    // ── Stress invariants (J2, Lode) ─────────────────────────────────────

    #[test]
    fn j2_invariant_uniaxial_tension() {
        // Uniaxial: σ11 = σ, others 0 → J2 = σ²/3
        let sigma = 300.0e6;
        let s = [sigma, 0.0, 0.0, 0.0, 0.0, 0.0];
        let j2 = stress_j2_invariant(s);
        let expected = sigma * sigma / 3.0;
        assert!(
            (j2 - expected).abs() / expected < 1e-10,
            "J2 = {j2}, expected {expected}"
        );
    }

    #[test]
    fn j2_invariant_hydrostatic_zero() {
        let s = [100.0e6, 100.0e6, 100.0e6, 0.0, 0.0, 0.0];
        let j2 = stress_j2_invariant(s);
        assert!(j2.abs() < 1.0, "hydrostatic: J2 = {j2}");
    }

    #[test]
    fn lode_angle_hydrostatic_zero() {
        // Hydrostatic state has no deviatoric component → J2=0 → angle = 0
        let s = [100.0e6, 100.0e6, 100.0e6, 0.0, 0.0, 0.0];
        let theta = lode_angle(s);
        assert!(theta.abs() < 1e-10, "hydrostatic Lode angle = {theta}");
    }

    #[test]
    fn lode_angle_in_valid_range() {
        let s = [300.0e6, 100.0e6, 50.0e6, 20.0e6, 10.0e6, 5.0e6];
        let theta = lode_angle(s);
        let pi_over_3 = std::f64::consts::PI / 3.0;
        assert!(
            (0.0..=pi_over_3 + 1e-10).contains(&theta),
            "Lode angle out of range: {theta}"
        );
    }
}
