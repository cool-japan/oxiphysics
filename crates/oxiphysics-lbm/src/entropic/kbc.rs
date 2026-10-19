// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Natural-moment decomposition and entropic stabiliser of the
//! Karlin–Bösch–Chikatamarla (KBC) collision for D2Q9.
//!
//! The populations are represented by the nine natural moments
//! `M_pq = Σ_i f_i c_ix^p c_iy^q` (`p, q ∈ {0, 1, 2}`), grouped as
//!
//! - **k** (kinematic, conserved): `ρ = M_00`, `j_x = M_10`, `j_y = M_01`;
//! - **s** (shear): the deviatoric stress `N = M_20 − M_02` and `Π_xy = M_11`;
//! - **h** (higher order): the trace `T = M_20 + M_02`, the third-order
//!   moments `Q_xyy = M_12`, `Q_yxx = M_21` and the fourth-order moment
//!   `A = M_22`.
//!
//! The non-equilibrium deviation `Δ = f − f_eq` has no k part (the
//! equilibrium carries the cell's density and momentum), so `Δ = Δs + Δh`.
//! The KBC collision is
//!
//! ```text
//! f' = f − β (2 Δs + γ Δh),
//! γ  = 1/β − (2 − 1/β) ⟨Δs|Δh⟩ / ⟨Δh|Δh⟩,   ⟨x|y⟩ = Σ_i x_i y_i / f_eq,i,
//! ```
//!
//! with `β = 1/(2τ)`, `τ = 3ν + 1/2`.  The shear moments are therefore relaxed
//! with the BGK frequency `2β = 1/τ` of the requested viscosity, while `γ`
//! minimises the (quadratic) entropy of the post-collision state over the
//! relaxation of the higher-order moments.  Because the trace `T` belongs to
//! the h part, the bulk viscosity is set by `γ` (for `⟨Δs|Δh⟩ = 0`,
//! `βγ = 1` and the bulk mode is relaxed straight to equilibrium).
//!
//! Reference: Karlin, Bösch & Chikatamarla (2014), Phys. Rev. E 90, 031302(R).

use super::C;

/// Relative magnitude (with respect to the populations) below which the
/// higher-order deviation `Δh` is indistinguishable from floating-point
/// rounding.
///
/// The stabiliser ratio `⟨Δs|Δh⟩ / ⟨Δh|Δh⟩` is invariant under a rescaling of
/// `Δh`, so it does not vanish when `Δh` does: a `Δh` made of rounding noise
/// would inject a spurious correction of the size of `Δs`.  Below
/// `⟨Δh|Δh⟩ ≤ (KBC_H_NOISE_FLOOR)² ρ` the h part is treated as zero and
/// `γ = 1/β` (the value for orthogonal s and h parts) is used.
pub const KBC_H_NOISE_FLOOR: f64 = 1e-12;

/// Shear (s) part of a deviation `delta` in the natural-moment basis.
///
/// The result carries exactly the moments `N = M_20 − M_02` and
/// `Π_xy = M_11` of `delta` and has every other natural moment (`ρ`, `j`,
/// `T`, `Q_xyy`, `Q_yxx`, `A`) equal to zero:
///
/// ```text
/// s(0,0) = 0,   s(±1,0) = N/4,   s(0,±1) = −N/4,   s(σx,σy) = σx σy Π_xy / 4.
/// ```
pub fn kbc_shear_part(delta: &[f64; 9]) -> [f64; 9] {
    let mut n_xy = 0.0;
    let mut pi_xy = 0.0;
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        n_xy += delta[i] * (cx * cx - cy * cy);
        pi_xy += delta[i] * cx * cy;
    }
    let mut s = [0.0; 9];
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        s[i] = 0.25 * (n_xy * (cx * cx - cy * cy) + pi_xy * cx * cy);
    }
    s
}

/// Split a non-equilibrium deviation `Δ = f − f_eq` into its shear part `Δs`
/// and its higher-order part `Δh = Δ − Δs` (trace, third- and fourth-order
/// moments, plus the rounding-level remainder of the conserved moments).
///
/// Returns `(Δs, Δh)` with `Δs + Δh = Δ`.
pub fn kbc_split(delta: &[f64; 9]) -> ([f64; 9], [f64; 9]) {
    let ds = kbc_shear_part(delta);
    let mut dh = [0.0; 9];
    for i in 0..9 {
        dh[i] = delta[i] - ds[i];
    }
    (ds, dh)
}

/// Entropic scalar product `⟨x|y⟩ = Σ_i x_i y_i / f_eq,i` (the Hessian of the
/// H-function at the equilibrium).
///
/// Returns `None` if an equilibrium population is not strictly positive or
/// the result is not finite.
pub fn entropic_scalar_product(x: &[f64; 9], y: &[f64; 9], f_eq: &[f64; 9]) -> Option<f64> {
    let mut acc = 0.0;
    for i in 0..9 {
        if f_eq[i] <= 0.0 {
            return None;
        }
        acc += x[i] * y[i] / f_eq[i];
    }
    acc.is_finite().then_some(acc)
}

/// KBC stabiliser `γ = 1/β − (2 − 1/β) ⟨Δs|Δh⟩ / ⟨Δh|Δh⟩`.
///
/// It is the exact minimiser over `γ` of the quadratic entropy estimate
/// `⟨Δ'|Δ'⟩` of the post-collision deviation
/// `Δ' = (1 − 2β) Δs + (1 − βγ) Δh`.  Falls back to `γ = 1/β` when `Δh` is at
/// the rounding level (see [`KBC_H_NOISE_FLOOR`]), when the scalar product is
/// undefined (non-positive equilibrium) or when the result is not finite.
pub fn kbc_stabiliser(beta: f64, ds: &[f64; 9], dh: &[f64; 9], f_eq: &[f64; 9]) -> f64 {
    let inv_beta = 1.0 / beta;
    let rho: f64 = f_eq.iter().sum();
    let floor = KBC_H_NOISE_FLOOR * KBC_H_NOISE_FLOOR * rho.abs();
    let (Some(sh), Some(hh)) = (
        entropic_scalar_product(ds, dh, f_eq),
        entropic_scalar_product(dh, dh, f_eq),
    ) else {
        return inv_beta;
    };
    if hh <= floor {
        return inv_beta;
    }
    let gamma = inv_beta - (2.0 - inv_beta) * sh / hh;
    if gamma.is_finite() { gamma } else { inv_beta }
}

/// One KBC collision of the populations `f` towards the equilibrium `f_eq`
/// (which must carry the density and momentum of `f`) with relaxation
/// parameter `β = 1/(2τ)`.
///
/// Returns the post-collision populations `f − β (2 Δs + γ Δh)` and the
/// stabiliser `γ` that was used.  Density and momentum are conserved exactly
/// up to rounding, because `Δs` has no conserved moments by construction and
/// `Δh = Δ − Δs` only carries the rounding-level k part of `Δ`.
pub fn kbc_collide_populations(f: &[f64; 9], f_eq: &[f64; 9], beta: f64) -> ([f64; 9], f64) {
    let mut delta = [0.0; 9];
    for i in 0..9 {
        delta[i] = f[i] - f_eq[i];
    }
    let (ds, dh) = kbc_split(&delta);
    let gamma = kbc_stabiliser(beta, &ds, &dh, f_eq);
    let mut post = [0.0; 9];
    for i in 0..9 {
        post[i] = f[i] - beta * (2.0 * ds[i] + gamma * dh[i]);
    }
    (post, gamma)
}
