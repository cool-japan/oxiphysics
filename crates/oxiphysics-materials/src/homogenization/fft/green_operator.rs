// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Periodic reference Green operator Γ⁰ in Fourier space (Moulinec–Suquet 1998).
//!
//! The isotropic reference Green operator for a homogeneous medium with Lamé
//! moduli (λ₀, μ₀) relates a Fourier-space polarization stress τ̂(ξ) to a strain
//! correction ε̃(ξ):
//!
//! ```text
//! ε̃_ij(ξ) = Γ⁰_{ijkl}(ξ) τ̂_kl(ξ)
//! ```
//!
//! With the normalized wave vector ξ̂ = ξ/|ξ|, the operator is
//!
//! ```text
//! Γ⁰_{ijkl}(ξ) =
//!     (1/(4μ₀)) (δ_ik ξ̂_j ξ̂_l + δ_il ξ̂_j ξ̂_k + δ_jk ξ̂_i ξ̂_l + δ_jl ξ̂_i ξ̂_k)
//!   - ((λ₀+μ₀)/(μ₀(λ₀+2μ₀))) ξ̂_i ξ̂_j ξ̂_k ξ̂_l
//! ```
//!
//! and `Γ⁰(0) = 0`.
//!
//! # Voigt conventions
//!
//! Tensors are passed as 6-component Voigt vectors with the ordering used by
//! [`crate::elastic::LinearElastic::stress_strain_matrix_3d`]:
//!
//! | Voigt | tensor (i,j) |
//! |-------|--------------|
//! | 0     | (0,0)        |
//! | 1     | (1,1)        |
//! | 2     | (2,2)        |
//! | 3     | (1,2)        |
//! | 4     | (0,2)        |
//! | 5     | (0,1)        |
//!
//! Stress is stored as the plain tensor components (σ₂₃ at index 3, etc.).
//! Strain uses the **engineering** convention: the shear entries (indices 3,4,5)
//! store γ = 2ε, matching the compliance/stiffness matrices of `LinearElastic`.
//! Applying Γ⁰ to a stress therefore yields an engineering-strain Voigt vector,
//! so that the full fixed-point loop `ε → σ = C:ε → −Γ⁰:σ` stays self-consistent
//! in the engineering convention.

/// Map a Voigt index (0..6) to the pair of 0-based tensor indices `(i, j)`.
#[inline]
pub(crate) const fn voigt_to_pair(v: usize) -> (usize, usize) {
    match v {
        0 => (0, 0),
        1 => (1, 1),
        2 => (2, 2),
        3 => (1, 2),
        4 => (0, 2),
        _ => (0, 1),
    }
}

/// Expand a stress Voigt vector into a full symmetric 3×3 tensor.
#[inline]
fn stress_voigt_to_full(tau: [f64; 6]) -> [[f64; 3]; 3] {
    let mut t = [[0.0_f64; 3]; 3];
    t[0][0] = tau[0];
    t[1][1] = tau[1];
    t[2][2] = tau[2];
    t[1][2] = tau[3];
    t[2][1] = tau[3];
    t[0][2] = tau[4];
    t[2][0] = tau[4];
    t[0][1] = tau[5];
    t[1][0] = tau[5];
    t
}

/// Collapse a full symmetric 3×3 strain tensor into an engineering-strain Voigt
/// vector. The off-diagonal (shear) components are doubled (γ = 2ε).
#[inline]
fn strain_full_to_voigt(e: [[f64; 3]; 3]) -> [f64; 6] {
    [
        e[0][0],
        e[1][1],
        e[2][2],
        e[1][2] + e[2][1],
        e[0][2] + e[2][0],
        e[0][1] + e[1][0],
    ]
}

/// Apply the periodic reference Green operator Γ⁰ to a polarization stress.
///
/// Returns the resulting strain as an **engineering-strain** Voigt vector
/// (`ε̃_ij = Γ⁰_{ijkl} τ_kl`).  For `ξ = 0` the result is the zero vector,
/// reflecting `Γ⁰(0) = 0`.
///
/// # Arguments
/// * `xi`        — wave vector ξ (not necessarily normalized; this routine
///   normalizes internally).
/// * `tau_voigt` — polarization stress in Voigt form (plain tensor components).
/// * `lambda0`   — reference Lamé first parameter λ₀.
/// * `mu0`       — reference shear modulus μ₀.
pub fn apply_gamma0(xi: [f64; 3], tau_voigt: [f64; 6], lambda0: f64, mu0: f64) -> [f64; 6] {
    let xi_norm_sq = xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2];
    if xi_norm_sq <= 0.0 {
        return [0.0; 6];
    }
    let inv_norm = 1.0 / xi_norm_sq.sqrt();
    let n = [xi[0] * inv_norm, xi[1] * inv_norm, xi[2] * inv_norm];

    let tau = stress_voigt_to_full(tau_voigt);

    // Precompute τ · n  (the vector  t_i = Σ_l τ_il n_l).
    let mut tn = [0.0_f64; 3];
    for (i, tni) in tn.iter_mut().enumerate() {
        *tni = tau[i][0] * n[0] + tau[i][1] * n[1] + tau[i][2] * n[2];
    }
    // Scalar  n · τ · n.
    let ntn = n[0] * tn[0] + n[1] * tn[1] + n[2] * tn[2];

    let inv_2mu = 1.0 / (2.0 * mu0);
    let longitudinal = (lambda0 + mu0) / (mu0 * (lambda0 + 2.0 * mu0));

    // ε̃_ij = Γ⁰_{ijkl} τ_kl
    //
    // The four-delta shear term of Γ⁰ carries the prefactor 1/(4μ₀):
    //   (1/(4μ₀))(δ_ik n_j n_l + δ_il n_j n_k + δ_jk n_i n_l + δ_jl n_i n_k) τ_kl.
    // Contracting against the symmetric τ collapses the four contributions to
    //   2 n_i (τn)_j + 2 n_j (τn)_i,
    // so the prefactor on [ n_i (τn)_j + n_j (τn)_i ] becomes 2·(1/(4μ₀)) = 1/(2μ₀).
    // The longitudinal correction subtracts (λ₀+μ₀)/(μ₀(λ₀+2μ₀)) · n_i n_j (nτn).
    let mut eps = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let shear = n[i] * tn[j] + n[j] * tn[i];
            eps[i][j] = inv_2mu * shear - longitudinal * n[i] * n[j] * ntn;
        }
    }

    strain_full_to_voigt(eps)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// Γ⁰ is even in ξ: Γ⁰(ξ) = Γ⁰(−ξ).
    #[test]
    fn test_gamma0_symmetry() {
        let lambda0 = 1.5;
        let mu0 = 0.8;
        let cases: [([f64; 3], [f64; 6]); 4] = [
            ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            ([0.3, -1.2, 0.7], [0.5, -0.2, 0.9, 0.4, -0.1, 0.3]),
            ([2.0, 2.0, -1.0], [1.0, 1.0, 1.0, 0.5, 0.5, 0.5]),
            ([0.0, 1.0, -3.0], [-0.7, 0.2, 0.4, 0.1, 0.6, -0.3]),
        ];
        for (xi, tau) in cases {
            let pos = apply_gamma0(xi, tau, lambda0, mu0);
            let neg = apply_gamma0([-xi[0], -xi[1], -xi[2]], tau, lambda0, mu0);
            for c in 0..6 {
                assert!(
                    approx(pos[c], neg[c], 1e-12),
                    "Γ⁰ not even at xi={xi:?} comp {c}: {} vs {}",
                    pos[c],
                    neg[c]
                );
            }
        }
    }

    /// For ξ=[1,0,0], τ=[1,0,0,0,0,0] (pure longitudinal stress along x),
    /// ε₁₁ = 1/(λ₀+2μ₀).
    #[test]
    fn test_gamma0_unit_longitudinal() {
        let lambda0 = 1.5;
        let mu0 = 0.8;
        let eps = apply_gamma0(
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            lambda0,
            mu0,
        );
        let expected = 1.0 / (lambda0 + 2.0 * mu0);
        assert!(
            approx(eps[0], expected, 1e-12),
            "ε₁₁ = {} expected {expected}",
            eps[0]
        );
        // Pure longitudinal loading along x produces no transverse normal strain
        // from Γ⁰ (the n_i n_j structure couples only to x).
        assert!(approx(eps[1], 0.0, 1e-12), "ε₂₂ = {} should be 0", eps[1]);
        assert!(approx(eps[2], 0.0, 1e-12), "ε₃₃ = {} should be 0", eps[2]);
        for (c, val) in eps[3..].iter().enumerate().map(|(i, v)| (i + 3, v)) {
            assert!(approx(*val, 0.0, 1e-12), "shear {c} = {} should be 0", val);
        }
    }

    /// Γ⁰(0) = 0.
    #[test]
    fn test_gamma0_zero_frequency() {
        let eps = apply_gamma0([0.0, 0.0, 0.0], [1.0, 2.0, 3.0, 0.4, 0.5, 0.6], 1.0, 1.0);
        for (c, val) in eps.iter().enumerate() {
            assert!(val.abs() < 1e-300, "Γ⁰(0) comp {c} nonzero: {}", val);
        }
    }

    /// Pure shear stress τ₁₂ with ξ along x gives the expected shear strain.
    ///
    /// For ξ̂ = e_x and τ = τ₁₂ (e_x⊗e_y + e_y⊗e_x) with τ₁₂ = 1:
    /// (τn)=[0,1,0], nτn=0, so ε₁₂ = (1/(2μ₀))·(n₁·(τn)₁ + ...) — by the
    /// contraction ε₁₂ = ε₂₁ = 1/(2μ₀), hence engineering γ₁₂ = ε₁₂+ε₂₁ = 1/μ₀.
    #[test]
    fn test_gamma0_shear_x_direction() {
        let mu0 = 0.8;
        let lambda0 = 1.5;
        let eps = apply_gamma0(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            lambda0,
            mu0,
        );
        let expected_gamma = 1.0 / mu0;
        assert!(
            approx(eps[5], expected_gamma, 1e-12),
            "γ₁₂ = {} expected {expected_gamma}",
            eps[5]
        );
    }

    /// Consistency check: applying Γ⁰ to a stress generated by the reference
    /// medium itself, σ = C⁰ : ε, should recover the deviatoric/compatible part.
    /// For a longitudinal strain ε = ε₁₁ e_x⊗e_x with ξ along x, Γ⁰ : (C⁰:ε)
    /// must return exactly ε (the field is already compatible & periodic).
    #[test]
    fn test_gamma0_reference_medium_recovers_strain() {
        let lambda0 = 1.5;
        let mu0 = 0.8;
        // ε = e_x⊗e_x  (ε₁₁ = 1)
        // σ = C⁰:ε  →  σ₁₁ = λ₀+2μ₀, σ₂₂ = σ₃₃ = λ₀.
        let sxx = lambda0 + 2.0 * mu0;
        let syy = lambda0;
        let tau = [sxx, syy, syy, 0.0, 0.0, 0.0];
        let eps = apply_gamma0([1.0, 0.0, 0.0], tau, lambda0, mu0);
        // Γ⁰ : σ should give back ε₁₁ = 1 and zero transverse for ξ ∥ x.
        assert!(approx(eps[0], 1.0, 1e-12), "ε₁₁ = {} expected 1", eps[0]);
        assert!(approx(eps[1], 0.0, 1e-12), "ε₂₂ = {} expected 0", eps[1]);
        assert!(approx(eps[2], 0.0, 1e-12), "ε₃₃ = {} expected 0", eps[2]);
    }
}
