// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Karlin-Bösch-Chikatamarla (KBC) decomposition of the D3Q19 non-equilibrium populations.
//!
//! The KBC collision `f* = f − β (2 Δs + γ Δh)` (see
//! [`KbcCollision`](super::KbcCollision)) needs the non-equilibrium part `Δf = f − feq` split
//! into
//!
//! - the *shear* part `Δs_i = w_i / (2 cs⁴) (c_i c_i) : Π_dev`, the second-order Hermite
//!   reconstruction of the traceless non-equilibrium stress `Π_dev = Π − tr(Π) I / 3`,
//!   `Π = Σ_i Δf_i c_i c_i`, and
//! - the *higher-order* part `Δh = Δf − Δs`, which keeps the trace of `Π` and every moment of
//!   order three and higher.
//!
//! `Δs` and `Δh` carry no mass and no momentum, and the second moment of `Δs` is exactly
//! `Π_dev`. The stabiliser `γ` extremises the quadratic entropy estimate of the post-collision
//! state in the entropic scalar product `⟨x|y⟩ = Σ_i x_i y_i / feq_i`.
//!
//! Reference: Karlin, Bösch, Chikatamarla, Phys. Rev. E 90, 031302(R) (2014).

use crate::lattice::{CS2, D3Q19_VELOCITIES, D3Q19_WEIGHTS};

use super::functions::compute_pi1_tensor;

/// Shear part `Δs` of the non-equilibrium populations for the non-equilibrium stress
/// `pi = Σ_i (f_i − feq_i) c_i c_i`, stored as `[xx, yy, zz, xy, xz, yz]`.
///
/// Only the traceless part of `pi` enters:
/// `Δs_i = w_i / (2 cs⁴) Σ_ab c_ia c_ib Π_dev,ab`.
pub fn kbc_shear_part_d3q19(pi: &[f64; 6]) -> [f64; 19] {
    let mean_normal = (pi[0] + pi[1] + pi[2]) / 3.0;
    let dev = [
        pi[0] - mean_normal,
        pi[1] - mean_normal,
        pi[2] - mean_normal,
        pi[3],
        pi[4],
        pi[5],
    ];
    let cs4 = CS2 * CS2;
    std::array::from_fn(|i| {
        let c = D3Q19_VELOCITIES[i];
        let (cx, cy, cz) = (f64::from(c[0]), f64::from(c[1]), f64::from(c[2]));
        let contraction = cx * cx * dev[0]
            + cy * cy * dev[1]
            + cz * cz * dev[2]
            + 2.0 * cx * cy * dev[3]
            + 2.0 * cx * cz * dev[4]
            + 2.0 * cy * cz * dev[5];
        D3Q19_WEIGHTS[i] / (2.0 * cs4) * contraction
    })
}

/// Split `f − feq` into the shear part `Δs` and the higher-order part `Δh = f − feq − Δs`.
///
/// Returns `(Δs, Δh)`. Both have zero mass and zero momentum whenever `f` and `feq` share them.
pub fn kbc_decompose_d3q19(f: &[f64; 19], feq: &[f64; 19]) -> ([f64; 19], [f64; 19]) {
    let delta_s = kbc_shear_part_d3q19(&compute_pi1_tensor(f, feq));
    let delta_h = std::array::from_fn(|i| f[i] - feq[i] - delta_s[i]);
    (delta_s, delta_h)
}

/// Entropic scalar product `⟨x|y⟩ = Σ_i x_i y_i / feq_i`.
pub fn kbc_entropic_product(x: &[f64; 19], y: &[f64; 19], feq: &[f64; 19]) -> f64 {
    x.iter()
        .zip(y.iter())
        .zip(feq.iter())
        .map(|((&xi, &yi), &feqi)| xi * yi / feqi)
        .sum()
}

/// Relative size below which the higher-order part `Δh` is treated as zero.
///
/// The ratio `⟨Δs|Δh⟩ / ⟨Δh|Δh⟩` is invariant under a rescaling of `Δh`, so it does not vanish
/// with `Δh`: a `Δh` made of round-off noise (for instance for a state that is an equilibrium
/// plus a pure shear part) would inject a spurious correction of the size of `Δs`. `Δh` is
/// treated as zero when `⟨Δh|Δh⟩ ≤ KBC_H_NOISE_FLOOR² Σ_i feq_i`, i.e. when it is below
/// `10⁻¹²` of the populations.
pub const KBC_H_NOISE_FLOOR: f64 = 1e-12;

/// KBC entropic stabiliser `γ = 1/β − (2 − 1/β) ⟨Δs|Δh⟩ / ⟨Δh|Δh⟩`.
///
/// It minimises the quadratic entropy estimate `⟨f* − feq|f* − feq⟩` of
/// `f* = f − β (2 Δs + γ Δh)` over `γ`. `γ = 2` (the BGK value) is returned when `Δh` vanishes
/// (below [`KBC_H_NOISE_FLOOR`], the value of `γ` is then immaterial), when `β` is not
/// positive, when the equilibrium populations are not all positive (the entropic product is
/// undefined) and when the result would not be finite.
pub fn kbc_stabilizer(beta: f64, delta_s: &[f64; 19], delta_h: &[f64; 19], feq: &[f64; 19]) -> f64 {
    const BGK_GAMMA: f64 = 2.0;
    if beta.is_nan() || beta <= 0.0 || feq.iter().any(|&x| x.is_nan() || x <= 0.0) {
        return BGK_GAMMA;
    }
    let hh = kbc_entropic_product(delta_h, delta_h, feq);
    let floor = KBC_H_NOISE_FLOOR * KBC_H_NOISE_FLOOR * feq.iter().sum::<f64>();
    if hh.is_nan() || hh <= floor {
        return BGK_GAMMA;
    }
    let sh = kbc_entropic_product(delta_s, delta_h, feq);
    let inv_beta = 1.0 / beta;
    let gamma = inv_beta - (2.0 - inv_beta) * sh / hh;
    if gamma.is_finite() { gamma } else { BGK_GAMMA }
}
