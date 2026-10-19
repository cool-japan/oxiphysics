// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Overlap integrals for the minimal s/p valence basis shared by the PM3 and
//! SCC-DFTB engines.
//!
//! [`BasisFunction`] carries no Cartesian direction, so the semi-empirical engines
//! tag every basis function with an [`OrbType`] (s, px, py, pz) and evaluate the
//! overlap of the contracted Cartesian Gaussians here.  The contraction
//! coefficients stored in [`BasisFunction`] already include the primitive
//! normalisation constants (see [`BasisFunction::sto3g_s`]), so no extra
//! normalisation is applied to the primitives.

use crate::quantum_chemistry::BasisFunction;
use std::f64::consts::PI;

/// Orbital type of a single basis-function slot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OrbType {
    /// s-type orbital.
    S,
    /// px-type orbital.
    Px,
    /// py-type orbital.
    Py,
    /// pz-type orbital.
    Pz,
}

impl OrbType {
    /// Cartesian axis (0 = x, 1 = y, 2 = z) of a p orbital, `None` for s.
    pub(crate) fn axis(self) -> Option<usize> {
        match self {
            OrbType::S => None,
            OrbType::Px => Some(0),
            OrbType::Py => Some(1),
            OrbType::Pz => Some(2),
        }
    }

    /// The three p orbital types in x, y, z order.
    pub(crate) const P_SHELL: [OrbType; 3] = [OrbType::Px, OrbType::Py, OrbType::Pz];
}

/// Overlap of two contracted Cartesian Gaussians with s or p angular parts.
///
/// For primitives with exponents α (on centre A) and β (on centre B) the
/// Gaussian product theorem gives P = (αA + βB)/γ with γ = α + β and
///
/// * ⟨s|s⟩ = K
/// * ⟨s_A|p_k,B⟩ = K·(P − B)_k
/// * ⟨p_k,A|s_B⟩ = K·(P − A)_k
/// * ⟨p_k,A|p_l,B⟩ = K·[(P − A)_k (P − B)_l + δ_kl/(2γ)]
///
/// where K = (π/γ)^{3/2} exp(−αβ|A − B|²/γ).
pub(crate) fn contracted_overlap(
    bf_i: &BasisFunction,
    ot_i: OrbType,
    bf_j: &BasisFunction,
    ot_j: OrbType,
) -> f64 {
    let a = bf_i.center;
    let b = bf_j.center;
    let ab2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
    let mut s = 0.0_f64;
    for (&alpha, &ci) in bf_i.exponents.iter().zip(bf_i.coefficients.iter()) {
        for (&beta, &cj) in bf_j.exponents.iter().zip(bf_j.coefficients.iter()) {
            let gamma = alpha + beta;
            let k = (PI / gamma).powf(1.5) * (-(alpha * beta / gamma) * ab2).exp();
            let p = [
                (alpha * a[0] + beta * b[0]) / gamma,
                (alpha * a[1] + beta * b[1]) / gamma,
                (alpha * a[2] + beta * b[2]) / gamma,
            ];
            let angular = match (ot_i.axis(), ot_j.axis()) {
                (None, None) => 1.0,
                (None, Some(l)) => p[l] - b[l],
                (Some(k_ax), None) => p[k_ax] - a[k_ax],
                (Some(k_ax), Some(l)) => {
                    let mut v = (p[k_ax] - a[k_ax]) * (p[l] - b[l]);
                    if k_ax == l {
                        v += 1.0 / (2.0 * gamma);
                    }
                    v
                }
            };
            s += ci * cj * k * angular;
        }
    }
    s
}

/// Build the N×N overlap matrix of the tagged s/p basis.
///
/// The diagonal is renormalised to exactly one (S′_ij = S_ij / √(S_ii S_jj)) to
/// remove the residual rounding of the tabulated contraction coefficients.
pub(crate) fn build_overlap(basis: &[BasisFunction], orb_types: &[OrbType]) -> Vec<Vec<f64>> {
    let n = basis.len().min(orb_types.len());
    let mut s = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in i..n {
            let v = contracted_overlap(&basis[i], orb_types[i], &basis[j], orb_types[j]);
            s[i][j] = v;
            s[j][i] = v;
        }
    }
    let norms: Vec<f64> = (0..n)
        .map(|i| if s[i][i] > 0.0 { s[i][i].sqrt() } else { 1.0 })
        .collect();
    for (i, row) in s.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v /= norms[i] * norms[j];
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum_chemistry::CartesianAxis;

    #[test]
    fn normalised_s_and_p_self_overlap_is_one() {
        let s = BasisFunction::sto3g_s([0.1, -0.2, 0.3], 1.24);
        let p = BasisFunction::sto3g_p([0.1, -0.2, 0.3], 2.25, CartesianAxis::X);
        assert!((contracted_overlap(&s, OrbType::S, &s, OrbType::S) - 1.0).abs() < 1e-6);
        for ot in OrbType::P_SHELL {
            let v = contracted_overlap(&p, ot, &p, ot);
            assert!((v - 1.0).abs() < 1e-6, "p self-overlap {v}");
        }
    }

    #[test]
    fn same_centre_s_p_and_px_py_are_orthogonal() {
        let s = BasisFunction::sto3g_s([0.0; 3], 1.24);
        let p = BasisFunction::sto3g_p([0.0; 3], 2.25, CartesianAxis::X);
        assert!(contracted_overlap(&s, OrbType::S, &p, OrbType::Px).abs() < 1e-14);
        assert!(contracted_overlap(&p, OrbType::Px, &p, OrbType::Py).abs() < 1e-14);
    }

    #[test]
    fn s_p_overlap_is_antisymmetric_in_direction() {
        // An s orbital on +x overlaps a px orbital at the origin positively,
        // one on −x negatively.
        let p = BasisFunction::sto3g_p([0.0; 3], 2.25, CartesianAxis::X);
        let s_plus = BasisFunction::sto3g_s([1.5, 0.0, 0.0], 1.24);
        let s_minus = BasisFunction::sto3g_s([-1.5, 0.0, 0.0], 1.24);
        let v_plus = contracted_overlap(&p, OrbType::Px, &s_plus, OrbType::S);
        let v_minus = contracted_overlap(&p, OrbType::Px, &s_minus, OrbType::S);
        assert!(v_plus > 0.0);
        assert!((v_plus + v_minus).abs() < 1e-14);
        // and ⟨s|p⟩ = ⟨p|s⟩
        let v_swapped = contracted_overlap(&s_plus, OrbType::S, &p, OrbType::Px);
        assert!((v_plus - v_swapped).abs() < 1e-14);
    }

    #[test]
    fn overlap_matrix_is_positive_definite_for_water_basis() {
        let o = [0.0, 0.0, 0.0];
        let h1 = [0.0, 1.43, 1.11];
        let h2 = [0.0, -1.43, 1.11];
        let mut basis = vec![BasisFunction::sto3g_s(o, 2.6)];
        let mut types = vec![OrbType::S];
        for (ot, axis) in OrbType::P_SHELL.into_iter().zip(CartesianAxis::ALL) {
            basis.push(BasisFunction::sto3g_p(o, 2.6, axis));
            types.push(ot);
        }
        basis.push(BasisFunction::sto3g_s(h1, 1.24));
        types.push(OrbType::S);
        basis.push(BasisFunction::sto3g_s(h2, 1.24));
        types.push(OrbType::S);
        let s = build_overlap(&basis, &types);
        assert!(oxiphysics_core::numerical_methods::cholesky(&s).is_some());
    }
}
