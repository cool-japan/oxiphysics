#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biological growth mechanics: multiplicative decomposition, residual stress,
//! morphoelastic ODE, surface growth, and growth-driven buckling.
//!
//! Implements the continuum theory of finite growth (Rodriguez–Hoger–McCulloch,
//! 1994) using the multiplicative decomposition F = Fe · Fg, where Fg is the
//! inelastic growth tensor and Fe is the elastic accommodation tensor.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// Types
// ─────────────────────────────────────────────────────────────────────────────

/// A 3×3 matrix stored in row-major order (indices \[row*3 + col\]).
pub type Mat3 = [f64; 9];

/// A symmetric 6-component Voigt stress vector \[σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃\].
pub type Stress6 = [f64; 6];

// ─────────────────────────────────────────────────────────────────────────────
// Mat3 helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Return the identity 3×3 matrix.
pub fn mat3_identity() -> Mat3 {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}

/// Multiply two 3×3 matrices (row-major).
pub fn mat3_mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}

/// Compute the determinant of a 3×3 matrix.
pub fn mat3_det(m: Mat3) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}

/// Compute the inverse of a 3×3 matrix.  Panics if the matrix is singular
/// (|det| < 1e-300).
pub fn mat3_inv(m: Mat3) -> Mat3 {
    let det = mat3_det(m);
    assert!(
        det.abs() > 1e-300,
        "Cannot invert singular matrix (det = {det})"
    );
    let inv_det = 1.0 / det;
    [
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ]
}

/// Transpose a 3×3 matrix.
pub fn mat3_transpose(m: Mat3) -> Mat3 {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

/// Scale a matrix by a scalar.
pub fn mat3_scale(m: Mat3, s: f64) -> Mat3 {
    let mut out = m;
    for v in &mut out {
        *v *= s;
    }
    out
}

/// Add two matrices element-wise.
pub fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = a;
    for (i, bv) in b.iter().enumerate() {
        c[i] += bv;
    }
    c
}

/// Trace of a 3×3 matrix.
pub fn mat3_trace(m: Mat3) -> f64 {
    m[0] + m[4] + m[8]
}

// ─────────────────────────────────────────────────────────────────────────────
// GrowthTensor
// ─────────────────────────────────────────────────────────────────────────────

/// A growth tensor Fg encoding the reference-to-grown configuration map.
///
/// Stores a 3×3 matrix representing isotropic, transversely isotropic, or
/// fully anisotropic growth.
#[derive(Debug, Clone, Copy)]
pub struct GrowthTensor {
    /// The 3×3 growth tensor matrix (row-major).
    pub fg: Mat3,
}

impl GrowthTensor {
    /// Create an identity (no-growth) growth tensor.
    pub fn identity() -> Self {
        Self {
            fg: mat3_identity(),
        }
    }

    /// Create an isotropic volumetric growth tensor Fg = γ · I.
    ///
    /// # Arguments
    /// * `gamma` – Linear growth stretch ratio (γ > 1 for growth, 0 < γ < 1 for atrophy).
    pub fn isotropic(gamma: f64) -> Self {
        Self {
            fg: mat3_scale(mat3_identity(), gamma),
        }
    }

    /// Create a transversely isotropic growth tensor with different growth
    /// in the fibre direction and laterally.
    ///
    /// Fg = γ_f · (n ⊗ n) + γ_l · (I - n ⊗ n)
    ///
    /// # Arguments
    /// * `fibre_dir` – Unit fibre direction \[3\].
    /// * `gamma_fibre` – Stretch ratio along fibre direction.
    /// * `gamma_lateral` – Stretch ratio lateral to fibre.
    pub fn transversely_isotropic(
        fibre_dir: [f64; 3],
        gamma_fibre: f64,
        gamma_lateral: f64,
    ) -> Self {
        let n = fibre_dir;
        // n ⊗ n
        let nn = [
            n[0] * n[0],
            n[0] * n[1],
            n[0] * n[2],
            n[1] * n[0],
            n[1] * n[1],
            n[1] * n[2],
            n[2] * n[0],
            n[2] * n[1],
            n[2] * n[2],
        ];
        // I - n ⊗ n
        let i_minus_nn = mat3_add(mat3_identity(), mat3_scale(nn, -1.0));
        let term_f = mat3_scale(nn, gamma_fibre);
        let term_l = mat3_scale(i_minus_nn, gamma_lateral);
        Self {
            fg: mat3_add(term_f, term_l),
        }
    }

    /// Return the volumetric growth ratio Jg = det(Fg).
    pub fn volume_ratio(&self) -> f64 {
        mat3_det(self.fg)
    }

    /// Return the inverse of the growth tensor.
    pub fn inverse(&self) -> Mat3 {
        mat3_inv(self.fg)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Multiplicative decomposition
// ─────────────────────────────────────────────────────────────────────────────

/// Result of the multiplicative decomposition F = Fe · Fg.
#[derive(Debug, Clone, Copy)]
pub struct DecompositionResult {
    /// Total deformation gradient F.
    pub f_total: Mat3,
    /// Growth tensor Fg.
    pub fg: Mat3,
    /// Elastic deformation gradient Fe = F · Fg⁻¹.
    pub fe: Mat3,
    /// Elastic volume ratio Je = det(Fe).
    pub je: f64,
    /// Growth volume ratio Jg = det(Fg).
    pub jg: f64,
    /// Total volume ratio J = det(F).
    pub j_total: f64,
}

/// Perform the multiplicative decomposition F = Fe · Fg.
///
/// Given the total deformation gradient F and growth tensor Fg, extracts the
/// elastic part Fe = F · Fg⁻¹.
///
/// # Arguments
/// * `f_total` – Total deformation gradient (3×3 row-major).
/// * `growth` – Growth tensor Fg.
pub fn multiplicative_decomposition(f_total: Mat3, growth: &GrowthTensor) -> DecompositionResult {
    let fg_inv = growth.inverse();
    let fe = mat3_mul(f_total, fg_inv);
    let je = mat3_det(fe);
    let jg = growth.volume_ratio();
    let j_total = mat3_det(f_total);
    DecompositionResult {
        f_total,
        fg: growth.fg,
        fe,
        je,
        jg,
        j_total,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Residual stress from constrained growth (Neo-Hookean)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Cauchy stress tensor \[Pa\] from constrained growth using a
/// compressible Neo-Hookean elastic law:
///
/// σ = (μ/J_e) · (Be − I) + (λ·ln(Je)/Je) · I
///
/// where Be = Fe · FeT is the elastic left Cauchy-Green tensor.
///
/// # Arguments
/// * `fe` – Elastic deformation gradient (3×3 row-major).
/// * `mu` – Shear modulus \[Pa\].
/// * `lambda` – First Lamé parameter \[Pa\].
///
/// Returns the 3×3 Cauchy stress tensor (row-major).
pub fn growth_residual_stress(fe: Mat3, mu: f64, lambda: f64) -> Mat3 {
    let fe_t = mat3_transpose(fe);
    let be = mat3_mul(fe, fe_t); // Fe · FeT
    let je = mat3_det(fe);
    let je_safe = if je.abs() < 1e-300 { 1e-300 } else { je };
    let inv_je = 1.0 / je_safe;
    let ln_je = je_safe.ln();
    let identity = mat3_identity();
    // σ = (μ/Je)*(Be - I) + (λ*ln(Je)/Je)*I
    let be_minus_i = mat3_add(be, mat3_scale(identity, -1.0));
    let term1 = mat3_scale(be_minus_i, mu * inv_je);
    let term2 = mat3_scale(identity, lambda * ln_je * inv_je);
    mat3_add(term1, term2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Morphoelastic ODE for growth tensor evolution
// ─────────────────────────────────────────────────────────────────────────────

/// Growth stimulus function: a simple stress-driven growth rate.
///
/// Returns dγ/dt = k_g * max(0, σ_mean - σ_threshold) for isotropic growth.
///
/// # Arguments
/// * `sigma_mean` – Mean Cauchy stress (hydrostatic) \[Pa\].
/// * `sigma_threshold` – Growth threshold stress \[Pa\].
/// * `growth_rate` – Rate coefficient k_g \[1/(Pa·s)\].
pub fn growth_stimulus(sigma_mean: f64, sigma_threshold: f64, growth_rate: f64) -> f64 {
    growth_rate * (sigma_mean - sigma_threshold).max(0.0)
}

/// Advance the isotropic growth tensor by one explicit Euler step.
///
/// dγ/dt = growth_stimulus(σ_mean, σ_thr, k_g)
/// γ(t + dt) = γ(t) + dγ/dt * dt
///
/// Returns the updated growth stretch ratio γ.
///
/// # Arguments
/// * `gamma_current` – Current isotropic growth stretch γ.
/// * `sigma_cauchy` – Current Cauchy stress tensor (3×3 row-major).
/// * `sigma_threshold` – Growth threshold stress \[Pa\].
/// * `growth_rate` – Rate coefficient k_g \[1/(Pa·s)\].
/// * `dt` – Time step \[s\].
pub fn morphoelastic_evolution(
    gamma_current: f64,
    sigma_cauchy: Mat3,
    sigma_threshold: f64,
    growth_rate: f64,
    dt: f64,
) -> f64 {
    let sigma_mean = mat3_trace(sigma_cauchy) / 3.0;
    let dgamma_dt = growth_stimulus(sigma_mean, sigma_threshold, growth_rate);
    gamma_current + dgamma_dt * dt
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface area growth
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the surface area growth ratio for a thin shell undergoing
/// in-plane growth with stretch ratios λ₁ and λ₂.
///
/// A_grown / A_rest = λ₁ · λ₂
///
/// # Arguments
/// * `lambda1` – Principal in-plane growth stretch in direction 1.
/// * `lambda2` – Principal in-plane growth stretch in direction 2.
pub fn area_growth(lambda1: f64, lambda2: f64) -> f64 {
    lambda1 * lambda2
}

/// Compute the grown area of a surface patch.
///
/// # Arguments
/// * `rest_area` – Reference area \[m²\].
/// * `lambda1` – In-plane stretch ratio 1.
/// * `lambda2` – In-plane stretch ratio 2.
pub fn grown_area(rest_area: f64, lambda1: f64, lambda2: f64) -> f64 {
    rest_area * area_growth(lambda1, lambda2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Growth-driven buckling
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the critical surface growth ratio for wrinkling / buckling of a
/// thin film on a soft substrate (Huang et al. model).
///
/// γ_cr ≈ 1 + (3/4) * (E_s / E_f)^(2/3)
///
/// where E_f is the film modulus and E_s is the substrate modulus.
///
/// # Arguments
/// * `e_film` – Film (surface layer) Young's modulus \[Pa\].
/// * `e_substrate` – Substrate Young's modulus \[Pa\].
pub fn growth_driven_buckling(e_film: f64, e_substrate: f64) -> f64 {
    1.0 + 0.75 * (e_substrate / e_film).powf(2.0 / 3.0)
}

/// Estimate the critical axial growth for Euler buckling of a constrained rod.
///
/// A rod of length L with bending stiffness EI fixed at both ends buckles when
/// the compressive force exceeds the Euler critical load:
/// P_cr = 4π² · EI / L²
///
/// The corresponding critical growth strain: ε_cr = P_cr / (EA)
///
/// # Arguments
/// * `length` – Rod length \[m\].
/// * `ei` – Bending stiffness EI \[N·m²\].
/// * `ea` – Axial stiffness EA \[N\].
pub fn euler_buckling_growth_strain(length: f64, ei: f64, ea: f64) -> f64 {
    let p_cr = 4.0 * std::f64::consts::PI * std::f64::consts::PI * ei / (length * length);
    p_cr / ea
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Identity growth tensor has volume ratio 1.
    #[test]
    fn test_growth_tensor_identity_volume() {
        let gt = GrowthTensor::identity();
        assert!(
            (gt.volume_ratio() - 1.0).abs() < EPS,
            "Identity growth must have Jg = 1"
        );
    }

    // 2. Isotropic growth: Jg = γ³.
    #[test]
    fn test_isotropic_growth_volume() {
        let gamma = 1.3;
        let gt = GrowthTensor::isotropic(gamma);
        let expected = gamma * gamma * gamma;
        let actual = gt.volume_ratio();
        assert!(
            (actual - expected).abs() < EPS,
            "Isotropic growth Jg must be γ³ = {expected}, got {actual}"
        );
    }

    // 3. Multiplicative decomposition: identity growth → Fe = F.
    #[test]
    fn test_decomposition_identity_growth() {
        let f = [1.2, 0.1, 0.0, 0.0, 1.1, 0.05, 0.0, 0.0, 0.9_f64];
        let gt = GrowthTensor::identity();
        let result = multiplicative_decomposition(f, &gt);
        for i in 0..9 {
            assert!(
                (result.fe[i] - f[i]).abs() < EPS,
                "Fe should equal F when Fg = I, mismatch at index {i}"
            );
        }
    }

    // 4. Decomposition: J_total = Je * Jg (multiplicativity of determinants).
    #[test]
    fn test_decomposition_volume_multiplicativity() {
        let f = [1.5, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 0.8_f64];
        let gt = GrowthTensor::isotropic(1.1);
        let result = multiplicative_decomposition(f, &gt);
        let product = result.je * result.jg;
        assert!(
            (product - result.j_total).abs() < 1e-8,
            "J = Je * Jg must hold, got {product} vs {}",
            result.j_total
        );
    }

    // 5. Residual stress from identity elastic deformation is zero.
    #[test]
    fn test_residual_stress_identity_fe() {
        let fe = mat3_identity();
        let sigma = growth_residual_stress(fe, 1e5, 2e5);
        // For Fe = I: Be = I, Je = 1, so σ = (μ/1)*(I-I) + 0*I = 0.
        for (i, s) in sigma.iter().enumerate() {
            assert!(
                s.abs() < EPS,
                "Residual stress must be zero for Fe=I at index {i}, got {s}"
            );
        }
    }

    // 6. Residual stress is non-zero for elastic deformation (stretch).
    #[test]
    fn test_residual_stress_nonzero_for_stretch() {
        let fe = [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let sigma = growth_residual_stress(fe, 1e6, 2e6);
        // At least some stress component must be nonzero.
        let sigma_norm: f64 = sigma.iter().map(|s| s * s).sum::<f64>().sqrt();
        assert!(
            sigma_norm > 0.0,
            "Stretched configuration must have non-zero stress"
        );
    }

    // 7. Residual stress [0,0] (xx) is positive for tensile elastic stretch in x.
    #[test]
    fn test_residual_stress_sign_tension() {
        // Stretch x by factor 2, keep y and z at 1.
        let fe = [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let sigma = growth_residual_stress(fe, 1e6, 1e6);
        assert!(
            sigma[0] > 0.0,
            "Tensile x-stretch must give positive σ_xx, got {}",
            sigma[0]
        );
    }

    // 8. Growth stimulus is zero below threshold.
    #[test]
    fn test_growth_stimulus_below_threshold() {
        let stimulus = growth_stimulus(0.5e6, 1e6, 1e-10);
        assert_eq!(stimulus, 0.0, "Stimulus below threshold must be zero");
    }

    // 9. Growth stimulus is positive above threshold.
    #[test]
    fn test_growth_stimulus_above_threshold() {
        let stimulus = growth_stimulus(2e6, 1e6, 1e-10);
        assert!(stimulus > 0.0, "Stimulus above threshold must be positive");
    }

    // 10. Growth stimulus scales linearly with k_g.
    #[test]
    fn test_growth_stimulus_scales_with_rate() {
        let s1 = growth_stimulus(2e6, 1e6, 1e-10);
        let s2 = growth_stimulus(2e6, 1e6, 2e-10);
        assert!((s2 / s1 - 2.0).abs() < EPS, "Stimulus must scale with k_g");
    }

    // 11. Morphoelastic evolution with zero stimulus leaves γ unchanged.
    #[test]
    fn test_morphoelastic_no_change_below_threshold() {
        let sigma = mat3_scale(mat3_identity(), 0.5e6); // σ_mean = 0.5 MPa
        let gamma_new = morphoelastic_evolution(1.0, sigma, 1e6, 1e-10, 1.0);
        assert!(
            (gamma_new - 1.0).abs() < EPS,
            "γ must remain unchanged when σ < threshold"
        );
    }

    // 12. Morphoelastic evolution grows γ when σ > threshold.
    #[test]
    fn test_morphoelastic_grows() {
        let sigma = mat3_scale(mat3_identity(), 2e6); // σ_mean = 2 MPa
        let gamma_new = morphoelastic_evolution(1.0, sigma, 1e6, 1e-10, 1.0);
        assert!(gamma_new > 1.0, "γ must increase when σ > threshold");
    }

    // 13. Area growth ratio: λ₁ = λ₂ = 1 → ratio = 1.
    #[test]
    fn test_area_growth_identity() {
        let ratio = area_growth(1.0, 1.0);
        assert!(
            (ratio - 1.0).abs() < EPS,
            "Area growth with unit stretches must be 1"
        );
    }

    // 14. Area growth is multiplicative.
    #[test]
    fn test_area_growth_multiplicative() {
        let l1 = 1.5;
        let l2 = 2.0;
        let ratio = area_growth(l1, l2);
        assert!(
            (ratio - l1 * l2).abs() < EPS,
            "Area growth must equal λ₁·λ₂"
        );
    }

    // 15. Grown area = rest_area * area_growth_ratio.
    #[test]
    fn test_grown_area_formula() {
        let a0 = 4.0;
        let l1 = 1.3;
        let l2 = 1.7;
        let a = grown_area(a0, l1, l2);
        assert!((a - a0 * l1 * l2).abs() < EPS, "Grown area mismatch: {a}");
    }

    // 16. Buckling threshold is > 1 for E_s < E_f.
    #[test]
    fn test_buckling_threshold_greater_than_one() {
        let gamma_cr = growth_driven_buckling(1e9, 1e3); // very stiff film, soft substrate
        assert!(
            gamma_cr > 1.0,
            "Buckling threshold must be > 1, got {gamma_cr}"
        );
    }

    // 17. Buckling threshold decreases as substrate gets softer (lower E_s/E_f).
    #[test]
    fn test_buckling_threshold_decreases_softer_substrate() {
        // γ_cr = 1 + 0.75*(E_s/E_f)^(2/3): smaller E_s → smaller γ_cr.
        let ef = 1e9;
        let g1 = growth_driven_buckling(ef, 1e4); // stiffer substrate
        let g2 = growth_driven_buckling(ef, 1e3); // softer substrate
        assert!(
            g2 < g1,
            "Softer substrate (lower E_s) should give lower buckling threshold"
        );
    }

    // 18. Buckling threshold = 1.75 for E_s/E_f = 1.
    #[test]
    fn test_buckling_threshold_equal_moduli() {
        let gamma_cr = growth_driven_buckling(1e6, 1e6);
        assert!(
            (gamma_cr - 1.75).abs() < EPS,
            "Equal moduli: γ_cr must be 1.75, got {gamma_cr}"
        );
    }

    // 19. Transversely isotropic growth: volume ratio computed correctly.
    #[test]
    fn test_transversely_isotropic_volume() {
        // Fibre along x: γ_f in x, γ_l in y and z.
        let gt = GrowthTensor::transversely_isotropic([1.0, 0.0, 0.0], 2.0, 1.5);
        let jg = gt.volume_ratio();
        // Expected Jg = γ_f * γ_l^2 = 2.0 * 2.25 = 4.5
        let expected = 2.0 * 1.5 * 1.5;
        assert!(
            (jg - expected).abs() < 1e-8,
            "Transversely isotropic Jg mismatch: {jg} vs {expected}"
        );
    }

    // 20. mat3_identity determinant is 1.
    #[test]
    fn test_mat3_identity_det() {
        let det = mat3_det(mat3_identity());
        assert!(
            (det - 1.0).abs() < EPS,
            "Det of identity must be 1, got {det}"
        );
    }

    // 21. mat3_mul identity: A * I = A.
    #[test]
    fn test_mat3_mul_identity() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0_f64];
        let id = mat3_identity();
        let result = mat3_mul(a, id);
        for i in 0..9 {
            assert!(
                (result[i] - a[i]).abs() < EPS,
                "A*I must equal A at index {i}"
            );
        }
    }

    // 22. mat3_inv: A * A⁻¹ = I.
    #[test]
    fn test_mat3_inv_round_trip() {
        let a = [2.0, 1.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0, 4.0_f64];
        let a_inv = mat3_inv(a);
        let product = mat3_mul(a, a_inv);
        let identity = mat3_identity();
        for i in 0..9 {
            assert!(
                (product[i] - identity[i]).abs() < 1e-10,
                "A*A⁻¹ must equal I at index {i}, got {}",
                product[i]
            );
        }
    }

    // 23. Euler buckling strain: longer rod buckles at lower strain.
    #[test]
    fn test_euler_buckling_longer_rod_lower_strain() {
        let ei = 10.0;
        let ea = 1e4;
        let eps_short = euler_buckling_growth_strain(1.0, ei, ea);
        let eps_long = euler_buckling_growth_strain(2.0, ei, ea);
        assert!(
            eps_long < eps_short,
            "Longer rod must buckle at lower growth strain"
        );
    }

    // 24. Euler buckling strain is positive.
    #[test]
    fn test_euler_buckling_strain_positive() {
        let eps = euler_buckling_growth_strain(0.5, 5.0, 1000.0);
        assert!(eps > 0.0, "Buckling strain must be positive, got {eps}");
    }

    // 25. Decomposition of diagonal stretch: Fe is diagonal.
    #[test]
    fn test_decomposition_diagonal_stretch() {
        // F = diag(2, 3, 4), Fg = diag(2, 3, 4) → Fe = I
        let f = [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0_f64];
        let fg_mat = f;
        let gt = GrowthTensor { fg: fg_mat };
        let result = multiplicative_decomposition(f, &gt);
        let identity = mat3_identity();
        for i in 0..9 {
            assert!(
                (result.fe[i] - identity[i]).abs() < 1e-10,
                "Fe must be identity when F = Fg, mismatch at {i}: {}",
                result.fe[i]
            );
        }
    }
}
