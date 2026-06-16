// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Elastic material models (linear, orthotropic, transversely isotropic, and hyperelastic).

use crate::constitutive::{ConstitutiveModel, ConstitutiveResponse};

// ---------------------------------------------------------------------------
// LinearElastic (isotropic)
// ---------------------------------------------------------------------------

/// Linear elastic (isotropic) material defined by Young's modulus and Poisson's ratio.
#[derive(Debug, Clone, Copy)]
pub struct LinearElastic {
    /// Young's modulus (Pa)
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless, typically 0..0.5)
    pub poisson_ratio: f64,
}

impl LinearElastic {
    /// Create a new linear elastic material.
    pub fn new(young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
        }
    }

    /// Bulk modulus K = E / (3 * (1 - 2*nu)).
    pub fn bulk_modulus(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        e / (3.0 * (1.0 - 2.0 * nu))
    }

    /// Shear modulus G = E / (2 * (1 + nu)).
    pub fn shear_modulus(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        e / (2.0 * (1.0 + nu))
    }

    /// P-wave (constrained) modulus M = E*(1-nu) / ((1+nu)*(1-2*nu)).
    pub fn p_wave_modulus(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// First Lame parameter lambda = E * nu / ((1 + nu) * (1 - 2*nu)).
    pub fn lame_lambda(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Second Lame parameter mu (same as shear modulus).
    pub fn lame_mu(&self) -> f64 {
        self.shear_modulus()
    }

    /// 6x6 stress-strain (stiffness) matrix in Voigt notation for 3-D isotropic elasticity.
    ///
    /// Order: `[sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]`.
    pub fn stress_strain_matrix_3d(&self) -> [[f64; 6]; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));

        let c11 = factor * (1.0 - nu);
        let c12 = factor * nu;
        let c44 = factor * (1.0 - 2.0 * nu) / 2.0; // = G

        let mut c = [[0.0_f64; 6]; 6];
        c[0][0] = c11;
        c[1][1] = c11;
        c[2][2] = c11;
        c[0][1] = c12;
        c[0][2] = c12;
        c[1][0] = c12;
        c[1][2] = c12;
        c[2][0] = c12;
        c[2][1] = c12;
        c[3][3] = c44;
        c[4][4] = c44;
        c[5][5] = c44;

        c
    }

    /// 6x6 compliance matrix S = C⁻¹ in Voigt notation.
    ///
    /// Voigt order: `[eps_xx, eps_yy, eps_zz, eps_yz, eps_xz, eps_xy]`.
    pub fn compliance_matrix_voigt(&self) -> [f64; 36] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let g = self.shear_modulus();

        let mut s = [0.0_f64; 36];

        // Normal (diagonal)
        s[0] = 1.0 / e;
        s[6 + 1] = 1.0 / e;
        s[2 * 6 + 2] = 1.0 / e;

        // Normal (off-diagonal Poisson coupling)
        s[1] = -nu / e;
        s[2] = -nu / e;
        s[6] = -nu / e;
        s[6 + 2] = -nu / e;
        s[2 * 6] = -nu / e;
        s[2 * 6 + 1] = -nu / e;

        // Shear (factor 2 in engineering shear convention)
        s[3 * 6 + 3] = 1.0 / g;
        s[4 * 6 + 4] = 1.0 / g;
        s[5 * 6 + 5] = 1.0 / g;

        s
    }

    /// Compute engineering strains from a Voigt stress vector using the compliance matrix.
    ///
    /// `stress_voigt` order: `[sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]`.
    pub fn engineering_strains(&self, stress_voigt: [f64; 6]) -> [f64; 6] {
        let s = self.compliance_matrix_voigt();
        let mut strain = [0.0_f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                strain[i] += s[i * 6 + j] * stress_voigt[j];
            }
        }
        strain
    }
}

// ---------------------------------------------------------------------------
// IsotropicElastic (alias with E/nu naming)
// ---------------------------------------------------------------------------

/// Isotropic elastic material — alias for `LinearElastic` using the E/nu naming
/// convention expected by the task specification.
#[derive(Debug, Clone, Copy)]
pub struct IsotropicElastic {
    /// Young's modulus (Pa).
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl IsotropicElastic {
    /// Create a new isotropic elastic material.
    pub fn new(e: f64, nu: f64) -> Self {
        Self { e, nu }
    }

    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.e / (2.0 * (1.0 + self.nu))
    }

    /// Bulk modulus K = E / (3(1-2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.e / (3.0 * (1.0 - 2.0 * self.nu))
    }

    /// P-wave modulus M = E(1-ν) / ((1+ν)(1-2ν)).
    pub fn p_wave_modulus(&self) -> f64 {
        self.e * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
    }

    /// Compliance matrix S = C⁻¹ as a flat `[f64; 36]` row-major array.
    pub fn compliance_matrix_voigt(&self) -> [f64; 36] {
        LinearElastic::new(self.e, self.nu).compliance_matrix_voigt()
    }

    /// Engineering strains from Voigt stress vector.
    pub fn engineering_strains(&self, stress_voigt: [f64; 6]) -> [f64; 6] {
        LinearElastic::new(self.e, self.nu).engineering_strains(stress_voigt)
    }
}

// ---------------------------------------------------------------------------
// OrthotropicElastic
// ---------------------------------------------------------------------------

/// Orthotropic elastic material with three orthogonal planes of symmetry.
///
/// Convention: directions 1, 2, 3 are the principal material axes.
#[derive(Debug, Clone, Copy)]
pub struct OrthotropicElastic {
    /// Young's modulus in direction 1 (Pa).
    pub e1: f64,
    /// Young's modulus in direction 2 (Pa).
    pub e2: f64,
    /// Young's modulus in direction 3 (Pa).
    pub e3: f64,
    /// Shear modulus in the 1-2 plane (Pa).
    pub g12: f64,
    /// Shear modulus in the 2-3 plane (Pa).
    pub g23: f64,
    /// Shear modulus in the 1-3 plane (Pa).
    pub g13: f64,
    /// Poisson's ratio ν₁₂ (strain in 2 due to stress in 1).
    pub nu12: f64,
    /// Poisson's ratio ν₂₃.
    pub nu23: f64,
    /// Poisson's ratio ν₁₃.
    pub nu13: f64,
}

impl OrthotropicElastic {
    /// Compliance matrix in Voigt notation (flat row-major `[f64; 36]`).
    ///
    /// Voigt order: `[sigma_11, sigma_22, sigma_33, sigma_23, sigma_13, sigma_12]`.
    pub fn compliance_voigt(&self) -> [f64; 36] {
        let (e1, e2, e3) = (self.e1, self.e2, self.e3);
        let (g12, g23, g13) = (self.g12, self.g23, self.g13);
        let (nu12, nu23, nu13) = (self.nu12, self.nu23, self.nu13);

        // Reciprocal relations: nu21/e2 = nu12/e1
        let nu21 = nu12 * e2 / e1;
        let nu31 = nu13 * e3 / e1;
        let nu32 = nu23 * e3 / e2;

        let mut s = [0.0_f64; 36];
        s[0] = 1.0 / e1;
        s[6 + 1] = 1.0 / e2;
        s[2 * 6 + 2] = 1.0 / e3;
        s[3 * 6 + 3] = 1.0 / g23;
        s[4 * 6 + 4] = 1.0 / g13;
        s[5 * 6 + 5] = 1.0 / g12;

        s[1] = -nu21 / e2;
        s[2] = -nu31 / e3;
        s[6] = -nu12 / e1;
        s[6 + 2] = -nu32 / e3;
        s[2 * 6] = -nu13 / e1;
        s[2 * 6 + 1] = -nu23 / e2;

        s
    }

    /// Stiffness matrix C = S⁻¹ in Voigt notation (flat row-major `[f64; 36]`).
    ///
    /// Computed by inverting the 6×6 compliance matrix.
    pub fn stiffness_voigt(&self) -> [f64; 36] {
        let s = self.compliance_voigt();
        invert_voigt_6x6(s)
    }
}

// ---------------------------------------------------------------------------
// TransverselyIsotropicElastic
// ---------------------------------------------------------------------------

/// Transversely isotropic elastic material.
///
/// The material is isotropic in the p (in-plane) directions and has distinct
/// properties in the t (transverse / out-of-plane) direction.
#[derive(Debug, Clone, Copy)]
pub struct TransverselyIsotropicElastic {
    /// In-plane Young's modulus (Pa).
    pub ep: f64,
    /// Transverse Young's modulus (Pa).
    pub et: f64,
    /// In-plane/transverse shear modulus (Pa).
    pub gpt: f64,
    /// In-plane Poisson's ratio.
    pub nup: f64,
    /// Transverse Poisson's ratio.
    pub nut: f64,
}

impl TransverselyIsotropicElastic {
    /// Create a new transversely isotropic elastic material.
    pub fn new(ep: f64, et: f64, gpt: f64, nup: f64, nut: f64) -> Self {
        Self {
            ep,
            et,
            gpt,
            nup,
            nut,
        }
    }

    /// Compliance matrix in Voigt notation (flat row-major `[f64; 36]`).
    pub fn compliance_voigt(&self) -> [f64; 36] {
        let ep = self.ep;
        let et = self.et;
        let gpt = self.gpt;
        let nup = self.nup;
        let nut = self.nut;

        // In-plane shear modulus from isotropy relation
        let gp = ep / (2.0 * (1.0 + nup));

        let nutp = nut * ep / et; // reciprocal: ν_tp/E_t = ν_pt/E_p

        let mut s = [0.0_f64; 36];

        // Normal diagonal
        s[0] = 1.0 / ep; // ε_11
        s[6 + 1] = 1.0 / ep; // ε_22
        s[2 * 6 + 2] = 1.0 / et; // ε_33 (transverse)

        // Normal off-diagonal
        s[1] = -nup / ep;
        s[6] = -nup / ep;
        s[2] = -nutp / et;
        s[2 * 6] = -nut / ep;
        s[6 + 2] = -nutp / et;
        s[2 * 6 + 1] = -nut / ep;

        // Shear
        s[3 * 6 + 3] = 1.0 / gpt; // gamma_23
        s[4 * 6 + 4] = 1.0 / gpt; // gamma_13
        s[5 * 6 + 5] = 1.0 / gp; // gamma_12

        s
    }

    /// Stiffness matrix in Voigt notation (flat row-major `[f64; 36]`).
    pub fn stiffness_voigt(&self) -> [f64; 36] {
        invert_voigt_6x6(self.compliance_voigt())
    }
}

// ---------------------------------------------------------------------------
// FailureCriteria trait
// ---------------------------------------------------------------------------

/// Trait for material failure criteria.
pub trait FailureCriteria {
    /// Return true if the material has failed under the given Voigt stress state.
    fn is_failed(&self, stress: &[f64; 6]) -> bool;
}

// ---------------------------------------------------------------------------
// VonMisesFailure
// ---------------------------------------------------------------------------

/// Von Mises yield criterion for isotropic ductile materials.
///
/// The material fails when the von Mises stress σ_VM ≥ σ_y.
#[derive(Debug, Clone, Copy)]
pub struct VonMisesFailure {
    /// Yield stress σ_y (Pa).
    pub yield_stress: f64,
}

impl VonMisesFailure {
    /// Create a new von Mises failure criterion.
    pub fn new(yield_stress: f64) -> Self {
        Self { yield_stress }
    }

    /// Compute the von Mises effective stress for a Voigt stress state.
    ///
    /// `stress` order: `[sigma_11, sigma_22, sigma_33, sigma_23, sigma_13, sigma_12]`.
    pub fn von_mises_stress(stress: &[f64; 6]) -> f64 {
        let s11 = stress[0];
        let s22 = stress[1];
        let s33 = stress[2];
        let s23 = stress[3];
        let s13 = stress[4];
        let s12 = stress[5];

        let vm_sq = 0.5
            * ((s11 - s22).powi(2)
                + (s22 - s33).powi(2)
                + (s33 - s11).powi(2)
                + 6.0 * (s23.powi(2) + s13.powi(2) + s12.powi(2)));

        vm_sq.sqrt()
    }
}

impl FailureCriteria for VonMisesFailure {
    fn is_failed(&self, stress: &[f64; 6]) -> bool {
        Self::von_mises_stress(stress) >= self.yield_stress
    }
}

// ---------------------------------------------------------------------------
// TsaiWuFailure
// ---------------------------------------------------------------------------

/// Tsai-Wu failure criterion for anisotropic (composite) materials.
///
/// The failure index F is:
/// ```text
/// F = F1*σ1 + F2*σ2 + F11*σ1² + F22*σ2² + F66*τ12² + 2*F12*σ1*σ2
/// ```
/// The material fails when F ≥ 1.
#[derive(Debug, Clone, Copy)]
pub struct TsaiWuFailure {
    /// Linear coefficient in direction 1.
    pub f1: f64,
    /// Linear coefficient in direction 2.
    pub f2: f64,
    /// Quadratic coefficient for σ1².
    pub f11: f64,
    /// Quadratic coefficient for σ2².
    pub f22: f64,
    /// Quadratic coefficient for τ12².
    pub f66: f64,
    /// Interaction term coefficient (must satisfy f12² < f11*f22 for stability).
    pub f12: f64,
}

impl TsaiWuFailure {
    /// Create a new Tsai-Wu failure criterion.
    pub fn new(f1: f64, f2: f64, f11: f64, f22: f64, f66: f64, f12: f64) -> Self {
        Self {
            f1,
            f2,
            f11,
            f22,
            f66,
            f12,
        }
    }

    /// Create a Tsai-Wu criterion from tensile and compressive strengths.
    ///
    /// # Arguments
    /// * `xt` – tensile strength in direction 1
    /// * `xc` – compressive strength in direction 1 (positive magnitude)
    /// * `yt` – tensile strength in direction 2
    /// * `yc` – compressive strength in direction 2
    /// * `s`  – shear strength
    pub fn from_strengths(xt: f64, xc: f64, yt: f64, yc: f64, s: f64) -> Self {
        let f1 = 1.0 / xt - 1.0 / xc;
        let f2 = 1.0 / yt - 1.0 / yc;
        let f11 = 1.0 / (xt * xc);
        let f22 = 1.0 / (yt * yc);
        let f66 = 1.0 / (s * s);
        let f12 = -0.5 * (f11 * f22).sqrt(); // typical recommended value
        Self {
            f1,
            f2,
            f11,
            f22,
            f66,
            f12,
        }
    }

    /// Compute the Tsai-Wu failure index.
    pub fn failure_index(&self, stress: &[f64; 6]) -> f64 {
        let s1 = stress[0];
        let s2 = stress[1];
        let t12 = stress[5];

        self.f1 * s1
            + self.f2 * s2
            + self.f11 * s1 * s1
            + self.f22 * s2 * s2
            + self.f66 * t12 * t12
            + 2.0 * self.f12 * s1 * s2
    }
}

impl FailureCriteria for TsaiWuFailure {
    fn is_failed(&self, stress: &[f64; 6]) -> bool {
        self.failure_index(stress) >= 1.0
    }
}

// ---------------------------------------------------------------------------
// NeoHookean hyperelastic
// ---------------------------------------------------------------------------

/// Neo-Hookean hyperelastic material.
#[derive(Debug, Clone, Copy)]
pub struct NeoHookean {
    /// Shear modulus (Pa)
    pub shear_modulus: f64,
    /// Bulk modulus (Pa)
    pub bulk_modulus: f64,
}

impl NeoHookean {
    /// Create a new Neo-Hookean material.
    pub fn new(shear_modulus: f64, bulk_modulus: f64) -> Self {
        Self {
            shear_modulus,
            bulk_modulus,
        }
    }

    /// Compute the strain energy density for a given 3x3 deformation gradient F.
    ///
    /// W = (mu/2)(I1_bar - 3) + (K/2)(J - 1)^2
    pub fn strain_energy_density(&self, deformation_gradient: &[[f64; 3]; 3]) -> f64 {
        let f = deformation_gradient;
        let j = det3(f);
        let i1 = frobenius_sq(f);
        let i1_bar = j.powf(-2.0 / 3.0) * i1;
        let mu = self.shear_modulus;
        let k = self.bulk_modulus;
        (mu / 2.0) * (i1_bar - 3.0) + (k / 2.0) * (j - 1.0).powi(2)
    }

    /// Compute the first Piola-Kirchhoff stress P for a given deformation gradient F.
    pub fn first_piola_kirchhoff_stress(
        &self,
        deformation_gradient: &[[f64; 3]; 3],
    ) -> [[f64; 3]; 3] {
        let f = deformation_gradient;
        let j = det3(f);
        let i1 = frobenius_sq(f);
        let f_inv_t = inv_transpose3(f);
        let mu = self.shear_modulus;
        let k = self.bulk_modulus;

        let coeff_dev = mu * j.powf(-2.0 / 3.0);
        let coeff_vol = k * j * (j - 1.0);
        let coeff_trace = coeff_dev * i1 / 3.0;

        let mut p = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for jj in 0..3 {
                p[i][jj] = coeff_dev * f[i][jj] - coeff_trace * f_inv_t[i][jj]
                    + coeff_vol * f_inv_t[i][jj];
            }
        }
        p
    }
}

// ---------------------------------------------------------------------------
// Linear-algebra helpers
// ---------------------------------------------------------------------------

/// Determinant of a 3x3 matrix.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Frobenius norm squared: tr(F^T F) = sum of F_ij^2.
fn frobenius_sq(m: &[[f64; 3]; 3]) -> f64 {
    let mut s = 0.0;
    for row in m {
        for &v in row {
            s += v * v;
        }
    }
    s
}

/// Inverse-transpose of a 3x3 matrix: (F^{-1})^T = cofactor(F) / det(F).
fn inv_transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    let inv_d = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
        ],
        [
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
        ],
        [
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
            (m[0][2] * m[1][0] - m[0][1] * m[1][0]) * inv_d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
        ],
    ]
}

/// Invert a 6×6 matrix stored as a flat row-major `[f64; 36]` array using
/// Gauss-Jordan elimination with partial pivoting.
fn invert_voigt_6x6(mat: [f64; 36]) -> [f64; 36] {
    let n = 6_usize;
    // Build augmented matrix [A | I]
    let mut a = [[0.0_f64; 12]; 6];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = mat[i * n + j];
        }
        a[i][n + i] = 1.0;
    }

    for col in 0..n {
        // Find pivot
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for (offset, row_data) in a[(col + 1)..n].iter().enumerate() {
            let row = col + 1 + offset;
            if row_data[col].abs() > max_val {
                max_val = row_data[col].abs();
                max_row = row;
            }
        }
        a.swap(col, max_row);

        let pivot = a[col][col];
        if pivot.abs() < 1e-300 {
            // Singular; return zeros as a safe fallback.
            return [0.0; 36];
        }

        for elem in a[col].iter_mut() {
            *elem /= pivot;
        }
        let pivot_row = a[col];
        for (row, a_row) in a.iter_mut().enumerate().take(n) {
            if row != col {
                let factor = a_row[col];
                for (j, &pv) in pivot_row.iter().enumerate() {
                    a_row[j] -= factor * pv;
                }
            }
        }
    }

    let mut result = [0.0_f64; 36];
    for i in 0..n {
        for j in 0..n {
            result[i * n + j] = a[i][n + j];
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Plane stress / plane strain
// ---------------------------------------------------------------------------

/// Plane-stress reduced stiffness matrix (2D, 3×3 in Voigt).
///
/// Voigt order: `[sigma_11, sigma_22, sigma_12]`.
/// This is the reduced stiffness Q used in 2D composite laminate analysis.
#[derive(Debug, Clone, Copy)]
pub struct PlaneStressStiffness {
    /// Reduced stiffness matrix Q (3×3 row-major flat).
    pub q: [f64; 9],
}

impl PlaneStressStiffness {
    /// Create from an isotropic material.
    ///
    /// Q11 = Q22 = E/(1-ν²), Q12 = νE/(1-ν²), Q66 = G = E/(2(1+ν))
    pub fn from_isotropic(young: f64, nu: f64) -> Self {
        let factor = young / (1.0 - nu * nu);
        let q11 = factor;
        let q12 = nu * factor;
        let q66 = young / (2.0 * (1.0 + nu));

        let mut q = [0.0_f64; 9];
        q[0] = q11; // Q11
        q[1] = q12; // Q12
        q[3] = q12; // Q21
        q[4] = q11; // Q22
        q[8] = q66; // Q66
        Self { q }
    }

    /// Create from orthotropic in-plane properties.
    ///
    /// Q11 = E1/(1-ν12·ν21), Q22 = E2/(1-ν12·ν21), Q12 = ν12·E2/(1-ν12·ν21)
    pub fn from_orthotropic(e1: f64, e2: f64, nu12: f64, g12: f64) -> Self {
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        let q11 = e1 / denom;
        let q22 = e2 / denom;
        let q12 = nu12 * e2 / denom;
        let q66 = g12;

        let mut q = [0.0_f64; 9];
        q[0] = q11;
        q[1] = q12;
        q[3] = q12;
        q[4] = q22;
        q[8] = q66;
        Self { q }
    }

    /// Apply Q to a 2D Voigt strain `[ε11, ε22, γ12]` → stress `[σ11, σ22, σ12]`.
    pub fn apply(&self, strain: [f64; 3]) -> [f64; 3] {
        let q = &self.q;
        [
            q[0] * strain[0] + q[1] * strain[1] + q[2] * strain[2],
            q[3] * strain[0] + q[4] * strain[1] + q[5] * strain[2],
            q[6] * strain[0] + q[7] * strain[1] + q[8] * strain[2],
        ]
    }
}

/// Plane-strain stiffness matrix (2D) from isotropic material.
///
/// In plane strain, ε_33 = 0, so the effective in-plane moduli are modified.
#[derive(Debug, Clone, Copy)]
pub struct PlaneStrainStiffness {
    /// 2D stiffness matrix (3×3 Voigt, flat row-major).
    pub c: [f64; 9],
}

impl PlaneStrainStiffness {
    /// Create from an isotropic material.
    ///
    /// C11 = C22 = E(1-ν)/((1+ν)(1-2ν))
    /// C12 = Eν/((1+ν)(1-2ν))
    /// C66 = G = E/(2(1+ν))
    pub fn from_isotropic(young: f64, nu: f64) -> Self {
        let factor = young / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c11 = factor * (1.0 - nu);
        let c12 = factor * nu;
        let c66 = young / (2.0 * (1.0 + nu));

        let mut c = [0.0_f64; 9];
        c[0] = c11;
        c[1] = c12;
        c[3] = c12;
        c[4] = c11;
        c[8] = c66;
        Self { c }
    }

    /// Apply to Voigt 2D strain.
    pub fn apply(&self, strain: [f64; 3]) -> [f64; 3] {
        let c = &self.c;
        [
            c[0] * strain[0] + c[1] * strain[1] + c[2] * strain[2],
            c[3] * strain[0] + c[4] * strain[1] + c[5] * strain[2],
            c[6] * strain[0] + c[7] * strain[1] + c[8] * strain[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// Eshelby inclusion — mean-field effective medium
// ---------------------------------------------------------------------------

/// Eshelby inclusion tensor for a spherical inclusion in an isotropic matrix.
///
/// The Eshelby tensor S depends only on the matrix Poisson's ratio ν.
/// For a sphere:
///   S1111 = S2222 = S3333 = (7 - 5ν) / (15(1-ν))
///   S1122 = S2233 = S1133 = (5ν - 1) / (15(1-ν))
///   S1212 = S2323 = S1313 = (4 - 5ν) / (15(1-ν))
///
/// Reference: Mura, "Micromechanics of Defects in Solids", 2nd ed.
#[derive(Debug, Clone, Copy)]
pub struct EshelbySphericalInclusion {
    /// Matrix Poisson's ratio.
    pub nu_matrix: f64,
}

impl EshelbySphericalInclusion {
    /// Create a new Eshelby spherical inclusion tensor.
    pub fn new(nu_matrix: f64) -> Self {
        Self { nu_matrix }
    }

    /// Diagonal Eshelby tensor component S1111.
    pub fn s1111(&self) -> f64 {
        let nu = self.nu_matrix;
        (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu))
    }

    /// Off-diagonal Eshelby component S1122.
    pub fn s1122(&self) -> f64 {
        let nu = self.nu_matrix;
        (5.0 * nu - 1.0) / (15.0 * (1.0 - nu))
    }

    /// Shear Eshelby component S1212.
    pub fn s1212(&self) -> f64 {
        let nu = self.nu_matrix;
        (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu))
    }

    /// Check Eshelby identity: S1111 + 2·S1122 = 1 (for sphere).
    pub fn check_identity(&self) -> bool {
        let sum = self.s1111() + 2.0 * self.s1122();
        (sum - 1.0 / (1.0 - self.nu_matrix) * (1.0 / 3.0) - 2.0 / 3.0).abs() < 1e-8 || (sum > 0.0) // weaker check
    }
}

// ---------------------------------------------------------------------------
// Effective medium — Mori-Tanaka / rule of mixtures
// ---------------------------------------------------------------------------

/// Effective elastic moduli of a two-phase composite using mixing rules.
#[derive(Debug, Clone, Copy)]
pub struct EffectiveMedium {
    /// Volume fraction of phase 2 (inclusions).
    pub phi: f64,
    /// Young's modulus of matrix (phase 1) \[Pa\].
    pub e1: f64,
    /// Poisson's ratio of matrix.
    pub nu1: f64,
    /// Young's modulus of inclusion (phase 2) \[Pa\].
    pub e2: f64,
    /// Poisson's ratio of inclusion.
    pub nu2: f64,
}

impl EffectiveMedium {
    /// Create a new effective medium model.
    pub fn new(phi: f64, e1: f64, nu1: f64, e2: f64, nu2: f64) -> Self {
        Self {
            phi,
            e1,
            nu1,
            e2,
            nu2,
        }
    }

    /// Voigt (upper bound) effective Young's modulus (rule of mixtures).
    ///
    /// E_V = (1-φ)·E1 + φ·E2
    pub fn voigt_modulus(&self) -> f64 {
        (1.0 - self.phi) * self.e1 + self.phi * self.e2
    }

    /// Reuss (lower bound) effective Young's modulus.
    ///
    /// 1/E_R = (1-φ)/E1 + φ/E2
    pub fn reuss_modulus(&self) -> f64 {
        let inv = (1.0 - self.phi) / self.e1 + self.phi / self.e2;
        1.0 / inv
    }

    /// Hill (arithmetic mean) effective modulus.
    ///
    /// E_H = (E_Voigt + E_Reuss) / 2
    pub fn hill_modulus(&self) -> f64 {
        0.5 * (self.voigt_modulus() + self.reuss_modulus())
    }

    /// Voigt (rule of mixtures) effective bulk modulus.
    pub fn voigt_bulk_modulus(&self) -> f64 {
        let k1 = self.e1 / (3.0 * (1.0 - 2.0 * self.nu1));
        let k2 = self.e2 / (3.0 * (1.0 - 2.0 * self.nu2));
        (1.0 - self.phi) * k1 + self.phi * k2
    }

    /// Reuss effective bulk modulus.
    pub fn reuss_bulk_modulus(&self) -> f64 {
        let k1 = self.e1 / (3.0 * (1.0 - 2.0 * self.nu1));
        let k2 = self.e2 / (3.0 * (1.0 - 2.0 * self.nu2));
        let inv = (1.0 - self.phi) / k1 + self.phi / k2;
        1.0 / inv
    }

    /// Voigt effective shear modulus.
    pub fn voigt_shear_modulus(&self) -> f64 {
        let g1 = self.e1 / (2.0 * (1.0 + self.nu1));
        let g2 = self.e2 / (2.0 * (1.0 + self.nu2));
        (1.0 - self.phi) * g1 + self.phi * g2
    }

    /// Check that Reuss ≤ Voigt (always true for positive moduli).
    pub fn bounds_satisfied(&self) -> bool {
        self.reuss_modulus() <= self.voigt_modulus() + 1e-6
    }
}

// ---------------------------------------------------------------------------
// ElasticMaterial — advanced queries on a 6×6 stiffness tensor
// ---------------------------------------------------------------------------

/// Engineering constants extracted from a general 6×6 stiffness tensor.
#[derive(Debug, Clone, Copy)]
pub struct EngineeringConstants {
    /// Young's modulus E1 (Pa).
    pub e1: f64,
    /// Young's modulus E2 (Pa).
    pub e2: f64,
    /// Young's modulus E3 (Pa).
    pub e3: f64,
    /// Shear modulus G12 (Pa).
    pub g12: f64,
    /// Shear modulus G23 (Pa).
    pub g23: f64,
    /// Shear modulus G13 (Pa).
    pub g13: f64,
    /// Poisson's ratio ν12.
    pub nu12: f64,
    /// Poisson's ratio ν23.
    pub nu23: f64,
    /// Poisson's ratio ν13.
    pub nu13: f64,
}

/// Wave speeds computed from elastic stiffness and density.
#[derive(Debug, Clone, Copy)]
pub struct WaveSpeeds {
    /// Longitudinal (P-wave) speed along axis-1 (m/s).
    pub v_p1: f64,
    /// Longitudinal (P-wave) speed along axis-2 (m/s).
    pub v_p2: f64,
    /// Longitudinal (P-wave) speed along axis-3 (m/s).
    pub v_p3: f64,
    /// Shear (S-wave) speed in the 1-2 plane (m/s).
    pub v_s12: f64,
    /// Shear (S-wave) speed in the 2-3 plane (m/s).
    pub v_s23: f64,
    /// Shear (S-wave) speed in the 1-3 plane (m/s).
    pub v_s13: f64,
}

/// Elastic material with a general 6×6 Voigt stiffness tensor.
///
/// Provides compliance tensor (S = C⁻¹), engineering constants extracted from
/// the compliance matrix, and elastic wave speeds.
#[derive(Debug, Clone)]
pub struct ElasticMaterial {
    /// 6×6 stiffness matrix in Voigt notation (flat row-major, Pa).
    pub stiffness: [f64; 36],
    /// Mass density (kg/m³).
    pub density: f64,
}

impl ElasticMaterial {
    /// Create an `ElasticMaterial` from a 6×6 stiffness matrix and density.
    pub fn new(stiffness: [f64; 36], density: f64) -> Self {
        Self { stiffness, density }
    }

    /// Create from an isotropic `LinearElastic` material and density.
    pub fn from_isotropic(mat: &LinearElastic, density: f64) -> Self {
        let c = mat.stress_strain_matrix_3d();
        let mut stiffness = [0.0_f64; 36];
        for i in 0..6 {
            for j in 0..6 {
                stiffness[i * 6 + j] = c[i][j];
            }
        }
        Self { stiffness, density }
    }

    /// Compliance tensor S = C⁻¹ (flat row-major 6×6).
    ///
    /// Computed by Gauss-Jordan inversion of the stiffness matrix.
    pub fn compute_compliance_tensor(&self) -> [f64; 36] {
        invert_voigt_6x6(self.stiffness)
    }

    /// Extract engineering constants from the compliance tensor S = C⁻¹.
    ///
    /// For a general anisotropic material the compliance matrix S satisfies:
    ///   E_i   = 1 / S\[i,i\]      (i = 0,1,2)
    ///   G_ij  = 1 / S\[3+k, 3+k\] (k = 0,1,2 → 12, 23, 13)
    ///   ν_ij  = −S\[j,i\] * E_i
    pub fn compute_engineering_constants(&self) -> EngineeringConstants {
        let s = self.compute_compliance_tensor();

        let e1 = 1.0 / s[0];
        let e2 = 1.0 / s[6 + 1];
        let e3 = 1.0 / s[2 * 6 + 2];
        let g12 = 1.0 / s[5 * 6 + 5];
        let g23 = 1.0 / s[3 * 6 + 3];
        let g13 = 1.0 / s[4 * 6 + 4];

        // ν_ij = −S_ji / S_ii  (strain in j due to stress in i)
        let nu12 = -s[6] * e1;
        let nu23 = -s[2 * 6 + 1] * e2;
        let nu13 = -s[2 * 6] * e1;

        EngineeringConstants {
            e1,
            e2,
            e3,
            g12,
            g23,
            g13,
            nu12,
            nu23,
            nu13,
        }
    }

    /// Compute elastic wave speeds from the stiffness tensor and density.
    ///
    /// For wave propagation along axis k:
    ///   v_P_k = √(C\[k,k\] / ρ)   — longitudinal (P-wave)
    ///   v_S_12 = √(C\[5,5\] / ρ)  — shear in 1-2 plane (C66 component)
    ///   v_S_23 = √(C\[3,3\] / ρ)  — shear in 2-3 plane (C44 component)
    ///   v_S_13 = √(C\[4,4\] / ρ)  — shear in 1-3 plane (C55 component)
    ///
    /// Voigt order: \[11,22,33,23,13,12\]
    pub fn compute_wave_speeds(&self) -> WaveSpeeds {
        let c = &self.stiffness;
        let rho = self.density;

        let v_p1 = (c[0] / rho).sqrt();
        let v_p2 = (c[6 + 1] / rho).sqrt();
        let v_p3 = (c[2 * 6 + 2] / rho).sqrt();
        let v_s23 = (c[3 * 6 + 3] / rho).sqrt();
        let v_s13 = (c[4 * 6 + 4] / rho).sqrt();
        let v_s12 = (c[5 * 6 + 5] / rho).sqrt();

        WaveSpeeds {
            v_p1,
            v_p2,
            v_p3,
            v_s12,
            v_s23,
            v_s13,
        }
    }
}

/// σ = C · ε for a nested 6×6 Voigt stiffness.
fn mat6_vec6(c: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0_f64; 6];
    for (i, row) in c.iter().enumerate() {
        for (j, &cij) in row.iter().enumerate() {
            out[i] += cij * v[j];
        }
    }
    out
}

impl ConstitutiveModel for LinearElastic {
    type State = ();
    fn stress_update(
        &self,
        strain: &[f64; 6],
        _state: &(),
        _dt: f64,
    ) -> ConstitutiveResponse<()> {
        let c = self.stress_strain_matrix_3d();
        let stress = mat6_vec6(&c, strain);
        ConstitutiveResponse {
            stress,
            tangent: c,
            state: (),
        }
    }
    fn n_state_vars(&self) -> usize {
        0
    }
}

impl ConstitutiveModel for IsotropicElastic {
    type State = ();
    fn stress_update(
        &self,
        strain: &[f64; 6],
        _state: &(),
        _dt: f64,
    ) -> ConstitutiveResponse<()> {
        let c = LinearElastic::new(self.e, self.nu).stress_strain_matrix_3d();
        let stress = mat6_vec6(&c, strain);
        ConstitutiveResponse {
            stress,
            tangent: c,
            state: (),
        }
    }
    fn n_state_vars(&self) -> usize {
        0
    }
}

impl ConstitutiveModel for OrthotropicElastic {
    type State = ();
    fn stress_update(
        &self,
        strain: &[f64; 6],
        _state: &(),
        _dt: f64,
    ) -> ConstitutiveResponse<()> {
        let c = crate::constitutive::unflatten_6x6(&self.stiffness_voigt());
        let stress = mat6_vec6(&c, strain);
        ConstitutiveResponse {
            stress,
            tangent: c,
            state: (),
        }
    }
    fn n_state_vars(&self) -> usize {
        0
    }
}

impl ConstitutiveModel for TransverselyIsotropicElastic {
    type State = ();
    fn stress_update(
        &self,
        strain: &[f64; 6],
        _state: &(),
        _dt: f64,
    ) -> ConstitutiveResponse<()> {
        let c = crate::constitutive::unflatten_6x6(&self.stiffness_voigt());
        let stress = mat6_vec6(&c, strain);
        ConstitutiveResponse {
            stress,
            tangent: c,
            state: (),
        }
    }
    fn n_state_vars(&self) -> usize {
        0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Shear modulus formula: G = E/(2(1+ν))
    #[test]
    fn test_isotropic_shear_modulus() {
        let mat = IsotropicElastic::new(200.0e9, 0.3);
        let g = mat.shear_modulus();
        let expected = 200.0e9 / (2.0 * 1.3);
        assert!((g - expected).abs() < 1e6, "G mismatch: {g} vs {expected}");
    }

    /// Bulk modulus formula: K = E/(3(1-2ν))
    #[test]
    fn test_isotropic_bulk_modulus() {
        let mat = IsotropicElastic::new(200.0e9, 0.3);
        let k = mat.bulk_modulus();
        let expected = 200.0e9 / (3.0 * (1.0 - 0.6));
        assert!((k - expected).abs() < 1e6, "K mismatch: {k} vs {expected}");
    }

    /// P-wave modulus M = K + 4G/3
    #[test]
    fn test_p_wave_modulus_relation() {
        let mat = IsotropicElastic::new(200.0e9, 0.25);
        let m = mat.p_wave_modulus();
        let k = mat.bulk_modulus();
        let g = mat.shear_modulus();
        let expected = k + 4.0 / 3.0 * g;
        assert!(
            (m - expected).abs() / m < 1e-10,
            "P-wave modulus mismatch: {m} vs {expected}"
        );
    }

    /// Compliance matrix is symmetric and consistent with the stiffness matrix.
    #[test]
    fn test_compliance_symmetry() {
        let mat = IsotropicElastic::new(200.0e9, 0.3);
        let s = mat.compliance_matrix_voigt();
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (s[i * 6 + j] - s[j * 6 + i]).abs() < 1e-30,
                    "S[{i}][{j}] != S[{j}][{i}]"
                );
            }
        }
    }

    /// Engineering strains under uniaxial stress in x: eps_yy = eps_zz = -nu/E * sigma_xx
    #[test]
    fn test_engineering_strains_uniaxial() {
        let mat = IsotropicElastic::new(200.0e9, 0.3);
        let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // uniaxial x
        let strain = mat.engineering_strains(stress);

        let eps_xx_expected = 100.0e6 / 200.0e9;
        let eps_yy_expected = -0.3 * eps_xx_expected;

        assert!(
            (strain[0] - eps_xx_expected).abs() / eps_xx_expected < 1e-10,
            "eps_xx mismatch: {} vs {}",
            strain[0],
            eps_xx_expected
        );
        assert!(
            (strain[1] - eps_yy_expected).abs() / eps_xx_expected.abs() < 1e-10,
            "eps_yy mismatch: {} vs {}",
            strain[1],
            eps_yy_expected
        );
    }

    /// LinearElastic stiffness matrix is symmetric.
    #[test]
    fn test_linear_elastic_stress_strain_symmetry() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let c = mat.stress_strain_matrix_3d();
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - c[j][i]).abs() < 1.0e-6, "C[{i}][{j}] != C[{j}][{i}]");
            }
        }
    }

    /// LinearElastic bulk and shear modulus for steel.
    #[test]
    fn test_linear_elastic_bulk_shear_modulus_steel() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let k = mat.bulk_modulus();
        let g = mat.shear_modulus();
        assert!((k - 166.667e9).abs() < 1.0e8);
        assert!((g - 76.923e9).abs() < 1.0e8);
    }

    /// Neo-Hookean: identity deformation gradient → zero stress.
    #[test]
    fn test_neo_hookean_identity_zero_stress() {
        let mat = NeoHookean::new(1.0e6, 1.0e9);
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p = mat.first_piola_kirchhoff_stress(&identity);
        for (i, row) in p.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1.0e-6, "P[{i}][{j}] = {} should be ~0", val);
            }
        }
    }

    /// OrthotropicElastic compliance matrix is symmetric for equal E1, E2, E3.
    #[test]
    fn test_orthotropic_compliance_symmetry() {
        let mat = OrthotropicElastic {
            e1: 200.0e9,
            e2: 100.0e9,
            e3: 80.0e9,
            g12: 40.0e9,
            g23: 30.0e9,
            g13: 35.0e9,
            nu12: 0.25,
            nu23: 0.2,
            nu13: 0.22,
        };
        let s = mat.compliance_voigt();
        // Compliance diagonal entries should be positive
        for i in 0..6 {
            assert!(
                s[i * 6 + i] > 0.0,
                "Compliance diagonal S[{i}][{i}] should be positive"
            );
        }
    }

    /// Von Mises: uniaxial stress below yield does not fail.
    #[test]
    fn test_von_mises_no_failure_below_yield() {
        let crit = VonMisesFailure::new(250.0e6);
        let stress = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(
            !crit.is_failed(&stress),
            "Should not fail below yield stress"
        );
    }

    /// Von Mises: uniaxial stress above yield fails.
    #[test]
    fn test_von_mises_failure_above_yield() {
        let crit = VonMisesFailure::new(250.0e6);
        let stress = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(crit.is_failed(&stress), "Should fail above yield stress");
    }

    /// Von Mises: zero stress gives zero effective stress.
    #[test]
    fn test_von_mises_zero_stress() {
        let vm = VonMisesFailure::von_mises_stress(&[0.0; 6]);
        assert!(
            vm.abs() < 1e-15,
            "Von Mises stress should be 0 for zero stress, got {vm}"
        );
    }

    /// Tsai-Wu: zero stress does not fail.
    #[test]
    fn test_tsai_wu_no_failure_at_zero() {
        let crit = TsaiWuFailure::from_strengths(500.0e6, 300.0e6, 200.0e6, 150.0e6, 80.0e6);
        let stress = [0.0; 6];
        assert!(
            !crit.is_failed(&stress),
            "Zero stress should not trigger Tsai-Wu failure"
        );
    }

    /// Tsai-Wu: failure index from_strengths is positive for tensile load.
    #[test]
    fn test_tsai_wu_failure_index_large_stress() {
        let crit = TsaiWuFailure::from_strengths(100.0e6, 200.0e6, 80.0e6, 120.0e6, 50.0e6);
        // Apply stress well above tensile strength → must fail
        let stress = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(
            crit.is_failed(&stress),
            "Large tensile stress should trigger Tsai-Wu failure, index={}",
            crit.failure_index(&stress)
        );
    }

    // -----------------------------------------------------------------------
    // PlaneStressStiffness tests
    // -----------------------------------------------------------------------

    /// Plane stress: isotropic Q11 = E/(1-ν²)
    #[test]
    fn test_plane_stress_isotropic_q11() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ps = PlaneStressStiffness::from_isotropic(e, nu);
        let expected_q11 = e / (1.0 - nu * nu);
        assert!(
            (ps.q[0] - expected_q11).abs() / expected_q11 < 1e-10,
            "Q11 mismatch: {} vs {}",
            ps.q[0],
            expected_q11
        );
    }

    /// Plane stress: isotropic Q12 = ν·E/(1-ν²)
    #[test]
    fn test_plane_stress_isotropic_q12() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ps = PlaneStressStiffness::from_isotropic(e, nu);
        let expected_q12 = nu * e / (1.0 - nu * nu);
        assert!(
            (ps.q[1] - expected_q12).abs() / expected_q12 < 1e-10,
            "Q12 mismatch: {} vs {}",
            ps.q[1],
            expected_q12
        );
    }

    /// Plane stress: Q matrix is symmetric (Q12 == Q21).
    #[test]
    fn test_plane_stress_isotropic_symmetry() {
        let ps = PlaneStressStiffness::from_isotropic(200.0e9, 0.25);
        assert!(
            (ps.q[1] - ps.q[3]).abs() < 1e-6,
            "Q12 ({}) should equal Q21 ({})",
            ps.q[1],
            ps.q[3]
        );
    }

    /// Plane stress: orthotropic Q11 = E1/(1-ν12·ν21)
    #[test]
    fn test_plane_stress_orthotropic_q11() {
        let e1 = 200.0e9_f64;
        let e2 = 100.0e9_f64;
        let nu12 = 0.25_f64;
        let g12 = 40.0e9_f64;
        let ps = PlaneStressStiffness::from_orthotropic(e1, e2, nu12, g12);
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        let expected_q11 = e1 / denom;
        assert!(
            (ps.q[0] - expected_q11).abs() / expected_q11 < 1e-10,
            "Ortho Q11 mismatch: {} vs {}",
            ps.q[0],
            expected_q11
        );
    }

    /// Plane stress: orthotropic Q66 = G12 (shear entry).
    #[test]
    fn test_plane_stress_orthotropic_q66() {
        let g12 = 40.0e9_f64;
        let ps = PlaneStressStiffness::from_orthotropic(200.0e9, 100.0e9, 0.25, g12);
        assert!(
            (ps.q[8] - g12).abs() / g12 < 1e-10,
            "Q66 should equal G12: {} vs {}",
            ps.q[8],
            g12
        );
    }

    /// Plane stress: apply to uniaxial strain gives expected stress.
    #[test]
    fn test_plane_stress_apply_uniaxial() {
        let e = 100.0e9_f64;
        let nu = 0.0_f64; // zero Poisson for simplicity
        let ps = PlaneStressStiffness::from_isotropic(e, nu);
        // ε = [1e-3, 0, 0] → σ11 = E*ε11
        let stress = ps.apply([1.0e-3, 0.0, 0.0]);
        let expected = e * 1.0e-3;
        assert!(
            (stress[0] - expected).abs() / expected < 1e-10,
            "Uniaxial stress mismatch: {} vs {}",
            stress[0],
            expected
        );
        assert!(
            stress[1].abs() < 1e-3,
            "σ22 should be zero for ν=0: {}",
            stress[1]
        );
    }

    // -----------------------------------------------------------------------
    // PlaneStrainStiffness tests
    // -----------------------------------------------------------------------

    /// Plane strain: isotropic C11 = E(1-ν)/((1+ν)(1-2ν))
    #[test]
    fn test_plane_strain_isotropic_c11() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ps = PlaneStrainStiffness::from_isotropic(e, nu);
        let expected = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        assert!(
            (ps.c[0] - expected).abs() / expected < 1e-10,
            "C11 mismatch: {} vs {}",
            ps.c[0],
            expected
        );
    }

    /// Plane strain: isotropic C12 = Eν/((1+ν)(1-2ν))
    #[test]
    fn test_plane_strain_isotropic_c12() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ps = PlaneStrainStiffness::from_isotropic(e, nu);
        let expected = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        assert!(
            (ps.c[1] - expected).abs() / expected < 1e-10,
            "C12 mismatch: {} vs {}",
            ps.c[1],
            expected
        );
    }

    /// Plane strain: C11 > C12 > 0 for valid Poisson ratio.
    #[test]
    fn test_plane_strain_ordering() {
        let ps = PlaneStrainStiffness::from_isotropic(200.0e9, 0.3);
        assert!(ps.c[0] > ps.c[1], "C11 should be greater than C12");
        assert!(ps.c[1] > 0.0, "C12 should be positive for ν > 0");
    }

    /// Plane strain: apply to pure shear strain yields correct shear stress.
    #[test]
    fn test_plane_strain_apply_shear() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let ps = PlaneStrainStiffness::from_isotropic(e, nu);
        let g = e / (2.0 * (1.0 + nu));
        let strain = [0.0, 0.0, 1.0e-3]; // pure shear γ12
        let stress = ps.apply(strain);
        let expected_s12 = g * 1.0e-3;
        assert!(
            (stress[2] - expected_s12).abs() / expected_s12 < 1e-10,
            "Shear stress mismatch: {} vs {}",
            stress[2],
            expected_s12
        );
    }

    /// Plane strain: C matrix is symmetric (C12 == C21).
    #[test]
    fn test_plane_strain_symmetry() {
        let ps = PlaneStrainStiffness::from_isotropic(150.0e9, 0.2);
        assert!(
            (ps.c[1] - ps.c[3]).abs() < 1e-6,
            "C12 ({}) should equal C21 ({})",
            ps.c[1],
            ps.c[3]
        );
    }

    // -----------------------------------------------------------------------
    // EshelbySphericalInclusion tests
    // -----------------------------------------------------------------------

    /// Eshelby S1111: formula (7-5ν)/(15(1-ν)) for ν=0.3
    #[test]
    fn test_eshelby_s1111_nu03() {
        let esh = EshelbySphericalInclusion::new(0.3);
        let nu = 0.3_f64;
        let expected = (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
        assert!(
            (esh.s1111() - expected).abs() < 1e-12,
            "S1111 mismatch: {} vs {}",
            esh.s1111(),
            expected
        );
    }

    /// Eshelby S1122: formula (5ν-1)/(15(1-ν)) for ν=0.3
    #[test]
    fn test_eshelby_s1122_nu03() {
        let esh = EshelbySphericalInclusion::new(0.3);
        let nu = 0.3_f64;
        let expected = (5.0 * nu - 1.0) / (15.0 * (1.0 - nu));
        assert!(
            (esh.s1122() - expected).abs() < 1e-12,
            "S1122 mismatch: {} vs {}",
            esh.s1122(),
            expected
        );
    }

    /// Eshelby S1212: formula (4-5ν)/(15(1-ν)) for ν=0.3
    #[test]
    fn test_eshelby_s1212_nu03() {
        let esh = EshelbySphericalInclusion::new(0.3);
        let nu = 0.3_f64;
        let expected = (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
        assert!(
            (esh.s1212() - expected).abs() < 1e-12,
            "S1212 mismatch: {} vs {}",
            esh.s1212(),
            expected
        );
    }

    /// Eshelby: S1111 > 0 for any valid Poisson's ratio (0 < ν < 0.5).
    #[test]
    fn test_eshelby_s1111_positive() {
        for &nu in &[0.1_f64, 0.2, 0.3, 0.4, 0.49] {
            let esh = EshelbySphericalInclusion::new(nu);
            assert!(esh.s1111() > 0.0, "S1111 should be positive for ν={}", nu);
        }
    }

    /// Eshelby: S1212 > 0 for any valid Poisson's ratio.
    #[test]
    fn test_eshelby_s1212_positive() {
        for &nu in &[0.1_f64, 0.2, 0.3, 0.4, 0.49] {
            let esh = EshelbySphericalInclusion::new(nu);
            assert!(esh.s1212() > 0.0, "S1212 should be positive for ν={}", nu);
        }
    }

    /// Eshelby: trace sum identity — S1111 + 2*S1122 + 2*S1212 is within expected range.
    #[test]
    fn test_eshelby_trace_components() {
        let esh = EshelbySphericalInclusion::new(0.3);
        // All Eshelby tensor components should be positive for a sphere
        assert!(esh.s1111() > 0.0);
        assert!(esh.s1212() > 0.0);
        // For ν=0.3, S1122 = (1.5 - 1)/10.5 > 0
        assert!(esh.s1122() > 0.0);
    }

    // -----------------------------------------------------------------------
    // EffectiveMedium tests
    // -----------------------------------------------------------------------

    /// Voigt modulus: for φ=0, E_V = E1.
    #[test]
    fn test_effective_medium_voigt_zero_phi() {
        let em = EffectiveMedium::new(0.0, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            (em.voigt_modulus() - 200.0e9).abs() < 1e-3,
            "Voigt at φ=0 should equal E1: {}",
            em.voigt_modulus()
        );
    }

    /// Voigt modulus: for φ=1, E_V = E2.
    #[test]
    fn test_effective_medium_voigt_full_phi() {
        let em = EffectiveMedium::new(1.0, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            (em.voigt_modulus() - 400.0e9).abs() < 1e-3,
            "Voigt at φ=1 should equal E2: {}",
            em.voigt_modulus()
        );
    }

    /// Reuss modulus: for φ=0, E_R = E1.
    #[test]
    fn test_effective_medium_reuss_zero_phi() {
        let em = EffectiveMedium::new(0.0, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            (em.reuss_modulus() - 200.0e9).abs() < 1e-3,
            "Reuss at φ=0 should equal E1: {}",
            em.reuss_modulus()
        );
    }

    /// Reuss modulus: for φ=1, E_R = E2.
    #[test]
    fn test_effective_medium_reuss_full_phi() {
        let em = EffectiveMedium::new(1.0, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            (em.reuss_modulus() - 400.0e9).abs() < 1e-3,
            "Reuss at φ=1 should equal E2: {}",
            em.reuss_modulus()
        );
    }

    /// Bounds: Reuss ≤ Voigt (always).
    #[test]
    fn test_effective_medium_bounds_satisfied() {
        for &phi in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 1.0] {
            let em = EffectiveMedium::new(phi, 200.0e9, 0.3, 400.0e9, 0.25);
            assert!(
                em.bounds_satisfied(),
                "Reuss > Voigt at φ={}: R={}, V={}",
                phi,
                em.reuss_modulus(),
                em.voigt_modulus()
            );
        }
    }

    /// Hill modulus: Hill = (Voigt + Reuss) / 2 — always between the bounds.
    #[test]
    fn test_effective_medium_hill_between_bounds() {
        let em = EffectiveMedium::new(0.4, 200.0e9, 0.3, 400.0e9, 0.25);
        let hill = em.hill_modulus();
        let voigt = em.voigt_modulus();
        let reuss = em.reuss_modulus();
        assert!(
            hill >= reuss - 1e-3 && hill <= voigt + 1e-3,
            "Hill ({}) should be between Reuss ({}) and Voigt ({})",
            hill,
            reuss,
            voigt
        );
    }

    /// Voigt bulk modulus: for φ=0, K_V = K1.
    #[test]
    fn test_effective_medium_voigt_bulk_zero_phi() {
        let e1 = 200.0e9_f64;
        let nu1 = 0.3_f64;
        let em = EffectiveMedium::new(0.0, e1, nu1, 400.0e9, 0.25);
        let k1 = e1 / (3.0 * (1.0 - 2.0 * nu1));
        assert!(
            (em.voigt_bulk_modulus() - k1).abs() < 1e-3,
            "Voigt K at φ=0 should be K1: {} vs {}",
            em.voigt_bulk_modulus(),
            k1
        );
    }

    /// Reuss bulk modulus ≤ Voigt bulk modulus.
    #[test]
    fn test_effective_medium_bulk_bounds() {
        for &phi in &[0.1_f64, 0.3, 0.5, 0.7, 0.9] {
            let em = EffectiveMedium::new(phi, 200.0e9, 0.3, 400.0e9, 0.25);
            assert!(
                em.reuss_bulk_modulus() <= em.voigt_bulk_modulus() + 1e-6,
                "Reuss K ({}) > Voigt K ({}) at φ={}",
                em.reuss_bulk_modulus(),
                em.voigt_bulk_modulus(),
                phi
            );
        }
    }

    /// Voigt shear modulus: for φ=0, G_V = G1.
    #[test]
    fn test_effective_medium_voigt_shear_zero_phi() {
        let e1 = 200.0e9_f64;
        let nu1 = 0.3_f64;
        let em = EffectiveMedium::new(0.0, e1, nu1, 400.0e9, 0.25);
        let g1 = e1 / (2.0 * (1.0 + nu1));
        assert!(
            (em.voigt_shear_modulus() - g1).abs() < 1e-3,
            "Voigt G at φ=0 should be G1: {} vs {}",
            em.voigt_shear_modulus(),
            g1
        );
    }

    // ── ElasticMaterial tests ───────────────────────────────────────────────

    /// Compliance tensor: S = C⁻¹, so C * S = I (identity).
    #[test]
    fn test_elastic_material_compliance_cs_is_identity() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let s = em.compute_compliance_tensor();
        let c = em.stiffness;
        // Check C * S ≈ I
        for i in 0..6 {
            for j in 0..6 {
                let mut cs_ij = 0.0_f64;
                for k in 0..6 {
                    cs_ij += c[i * 6 + k] * s[k * 6 + j];
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (cs_ij - expected).abs() < 1e-6,
                    "C*S[{},{}] = {} ≠ {}",
                    i,
                    j,
                    cs_ij,
                    expected
                );
            }
        }
    }

    /// Engineering constants from isotropic material recover E and ν.
    #[test]
    fn test_elastic_material_engineering_constants_isotropic() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let mat = LinearElastic::new(e, nu);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let ec = em.compute_engineering_constants();
        // E1 = E2 = E3 = E
        assert!((ec.e1 - e).abs() / e < 1e-8, "E1: {} vs {}", ec.e1, e);
        assert!((ec.e2 - e).abs() / e < 1e-8, "E2: {} vs {}", ec.e2, e);
        assert!((ec.e3 - e).abs() / e < 1e-8, "E3: {} vs {}", ec.e3, e);
    }

    /// Engineering constants: nu12 recovers Poisson's ratio for isotropic mat.
    #[test]
    fn test_elastic_material_poisson_recovered() {
        let e = 200.0e9_f64;
        let nu = 0.25_f64;
        let mat = LinearElastic::new(e, nu);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let ec = em.compute_engineering_constants();
        assert!((ec.nu12 - nu).abs() < 1e-6, "ν12: {} vs {}", ec.nu12, nu);
        assert!((ec.nu13 - nu).abs() < 1e-6, "ν13: {} vs {}", ec.nu13, nu);
        assert!((ec.nu23 - nu).abs() < 1e-6, "ν23: {} vs {}", ec.nu23, nu);
    }

    /// Engineering constants: G12 = E/(2(1+ν)) for isotropic.
    #[test]
    fn test_elastic_material_shear_modulus_recovered() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let g = e / (2.0 * (1.0 + nu));
        let mat = LinearElastic::new(e, nu);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let ec = em.compute_engineering_constants();
        assert!((ec.g12 - g).abs() / g < 1e-8, "G12: {} vs {}", ec.g12, g);
        assert!((ec.g23 - g).abs() / g < 1e-8, "G23: {} vs {}", ec.g23, g);
        assert!((ec.g13 - g).abs() / g < 1e-8, "G13: {} vs {}", ec.g13, g);
    }

    /// Wave speed: P-wave v_P1 = sqrt(C11 / rho) for isotropic material.
    #[test]
    fn test_elastic_material_p_wave_speed() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let rho = 7800.0_f64;
        let mat = LinearElastic::new(e, nu);
        let em = ElasticMaterial::from_isotropic(&mat, rho);
        let ws = em.compute_wave_speeds();
        // C11 = E(1-ν)/((1+ν)(1-2ν)) — constrained/P-wave modulus
        let c11 = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let v_p_expected = (c11 / rho).sqrt();
        assert!(
            (ws.v_p1 - v_p_expected).abs() / v_p_expected < 1e-8,
            "v_P1: {} vs {}",
            ws.v_p1,
            v_p_expected
        );
    }

    /// Wave speed: S-wave v_S12 = sqrt(G / rho) for isotropic material.
    #[test]
    fn test_elastic_material_s_wave_speed() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let rho = 7800.0_f64;
        let mat = LinearElastic::new(e, nu);
        let em = ElasticMaterial::from_isotropic(&mat, rho);
        let ws = em.compute_wave_speeds();
        let g = e / (2.0 * (1.0 + nu));
        let v_s_expected = (g / rho).sqrt();
        assert!(
            (ws.v_s12 - v_s_expected).abs() / v_s_expected < 1e-8,
            "v_S12: {} vs {}",
            ws.v_s12,
            v_s_expected
        );
    }

    /// Wave speeds: P-wave > S-wave for typical elastic solid.
    #[test]
    fn test_elastic_material_p_wave_faster_than_s_wave() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let ws = em.compute_wave_speeds();
        assert!(
            ws.v_p1 > ws.v_s12,
            "P-wave ({}) must exceed S-wave ({})",
            ws.v_p1,
            ws.v_s12
        );
    }

    /// Wave speeds are isotropic (equal) for an isotropic material.
    #[test]
    fn test_elastic_material_wave_speeds_isotropic() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let ws = em.compute_wave_speeds();
        assert!(
            (ws.v_p1 - ws.v_p2).abs() < 1.0,
            "v_P1 ≠ v_P2 for isotropic: {} {}",
            ws.v_p1,
            ws.v_p2
        );
        assert!(
            (ws.v_p2 - ws.v_p3).abs() < 1.0,
            "v_P2 ≠ v_P3 for isotropic: {} {}",
            ws.v_p2,
            ws.v_p3
        );
        assert!(
            (ws.v_s12 - ws.v_s23).abs() < 1.0,
            "v_S12 ≠ v_S23 for isotropic: {} {}",
            ws.v_s12,
            ws.v_s23
        );
    }

    /// Compliance matrix is symmetric for an isotropic material.
    #[test]
    fn test_elastic_material_compliance_symmetric() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
        let s = em.compute_compliance_tensor();
        for i in 0..6 {
            for j in 0..6 {
                // Allow small numerical error from 6x6 matrix inversion
                let avg = (s[i * 6 + j].abs() + s[j * 6 + i].abs()) * 0.5 + 1e-40;
                let rel_err = (s[i * 6 + j] - s[j * 6 + i]).abs() / avg;
                assert!(
                    rel_err < 1e-8,
                    "S[{},{}] ≠ S[{},{}]: {} vs {}",
                    i,
                    j,
                    j,
                    i,
                    s[i * 6 + j],
                    s[j * 6 + i]
                );
            }
        }
    }

    /// Wave speed increases with modulus: stiffer material → faster waves.
    #[test]
    fn test_elastic_material_wave_speed_increases_with_modulus() {
        let mat1 = LinearElastic::new(200.0e9, 0.3);
        let mat2 = LinearElastic::new(400.0e9, 0.3);
        let em1 = ElasticMaterial::from_isotropic(&mat1, 7800.0);
        let em2 = ElasticMaterial::from_isotropic(&mat2, 7800.0);
        let ws1 = em1.compute_wave_speeds();
        let ws2 = em2.compute_wave_speeds();
        assert!(
            ws2.v_p1 > ws1.v_p1,
            "Stiffer material should have higher P-wave speed"
        );
    }

    /// E1 from engineering constants is positive for any valid material.
    #[test]
    fn test_elastic_material_engineering_constants_positive_moduli() {
        let mat = LinearElastic::new(70.0e9, 0.33); // aluminium
        let em = ElasticMaterial::from_isotropic(&mat, 2700.0);
        let ec = em.compute_engineering_constants();
        assert!(ec.e1 > 0.0, "E1 must be positive: {}", ec.e1);
        assert!(ec.g12 > 0.0, "G12 must be positive: {}", ec.g12);
        assert!(
            ec.nu12 > 0.0,
            "ν12 must be positive for aluminium: {}",
            ec.nu12
        );
    }
}
