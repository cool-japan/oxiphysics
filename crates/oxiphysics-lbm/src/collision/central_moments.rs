// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact central-moment and cumulant transforms for the D3Q19 lattice.
//!
//! The D3Q19 velocity set has no corner velocities, so every monomial `x^a y^b z^c` with
//! `a, b, c <= 2` and all three exponents non-zero vanishes on the lattice nodes. The remaining
//! 19 monomials (exponent triples with at most two non-zero entries) are linearly independent on
//! the 19 nodes, so the 19 raw moments
//!
//! `m_abc = sum_i f_i c_x^a c_y^b c_z^c`
//!
//! carry exactly the information of the 19 populations. Shifting the lattice velocities by a
//! velocity `u` maps these monomials onto each other (a shifted monomial only expands into
//! monomials with a subset of the active axes), so the central moments
//!
//! `kappa_abc = sum_i f_i (c_x - u_x)^a (c_y - u_y)^b (c_z - u_z)^c`
//!
//! are an invertible, unit-triangular image of the raw moments for **every** `u`.
//!
//! This module provides the exact transforms between populations, raw moments, central moments
//! and cumulants. The ordering of the 19 entries is the one used by
//! [`compute_raw_moments_d3q19`](super::compute_raw_moments_d3q19) and
//! [`compute_central_moments_cumulant`](super::compute_central_moments_cumulant), see
//! [`MOMENT_EXPONENTS`].

/// Exponent triples `(a, b, c)` of the 19 moments `m_abc`, in storage order.
///
/// Order 0: `000`; order 1: `100, 010, 001`; order 2: `200, 020, 002, 110, 101, 011`;
/// order 3: `210, 201, 120, 021, 102, 012`; order 4: `220, 202, 022`.
pub const MOMENT_EXPONENTS: [[usize; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [2, 1, 0],
    [2, 0, 1],
    [1, 2, 0],
    [0, 2, 1],
    [1, 0, 2],
    [0, 1, 2],
    [2, 2, 0],
    [2, 0, 2],
    [0, 2, 2],
];

/// Position of the exponent triple `(a, b, c)` inside [`MOMENT_EXPONENTS`]
/// (`usize::MAX` for the eight triples that do not occur on D3Q19).
const MOMENT_INDEX: [[[usize; 3]; 3]; 3] = build_moment_index();

const fn build_moment_index() -> [[[usize; 3]; 3]; 3] {
    let mut table = [[[usize::MAX; 3]; 3]; 3];
    let mut k = 0;
    while k < 19 {
        let e = MOMENT_EXPONENTS[k];
        table[e[0]][e[1]][e[2]] = k;
        k += 1;
    }
    table
}

/// Binomial coefficients `C(n, k)` for `n, k <= 2`.
const BINOMIAL: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 2.0, 1.0]];

/// Translate the monomial moments by `shift`:
///
/// `out_alpha = sum_{beta <= alpha} prod_a C(alpha_a, beta_a) shift_a^(alpha_a - beta_a) src_beta`.
///
/// With `shift = -u` this maps raw moments to central moments about `u`; with `shift = +u` it
/// maps central moments back to raw moments. Because the active axes of `beta` are a subset of
/// those of `alpha`, every `beta` is one of the 19 stored moments.
fn translate_moments(src: &[f64; 19], shift: [f64; 3]) -> [f64; 19] {
    let mut out = [0.0_f64; 19];
    for (out_k, alpha) in out.iter_mut().zip(MOMENT_EXPONENTS.iter()) {
        let mut acc = 0.0_f64;
        for bx in 0..=alpha[0] {
            for by in 0..=alpha[1] {
                for bz in 0..=alpha[2] {
                    let coefficient = BINOMIAL[alpha[0]][bx]
                        * BINOMIAL[alpha[1]][by]
                        * BINOMIAL[alpha[2]][bz]
                        * shift[0].powi((alpha[0] - bx) as i32)
                        * shift[1].powi((alpha[1] - by) as i32)
                        * shift[2].powi((alpha[2] - bz) as i32);
                    acc += coefficient * src[MOMENT_INDEX[bx][by][bz]];
                }
            }
        }
        *out_k = acc;
    }
    out
}

/// Convert the 19 raw moments `m_abc` to central moments about the velocity `u`.
pub fn raw_to_central_moments_d3q19(m: &[f64; 19], u: [f64; 3]) -> [f64; 19] {
    translate_moments(m, [-u[0], -u[1], -u[2]])
}

/// Convert the 19 central moments about the velocity `u` back to raw moments `m_abc`.
///
/// Exact inverse of [`raw_to_central_moments_d3q19`] for the same `u`.
pub fn central_to_raw_moments_d3q19(kappa: &[f64; 19], u: [f64; 3]) -> [f64; 19] {
    translate_moments(kappa, u)
}

/// Reconstruct the 19 populations from their 19 raw moments.
///
/// This is the exact inverse of
/// [`compute_raw_moments_d3q19`](super::compute_raw_moments_d3q19):
/// `f = M^{-1} m`. The inversion is closed-form. The `xy`, `xz` and `yz` edge groups of four
/// populations are determined by the four moments `m_ab0` (`a, b in {1, 2}`) of their plane
/// (a 2x2 sign transform), the face populations follow from the even and odd moments along their
/// axis after subtracting the edge contributions, and the rest population from the density.
pub fn raw_moments_to_populations_d3q19(m: &[f64; 19]) -> [f64; 19] {
    let [
        m000,
        m100,
        m010,
        m001,
        m200,
        m020,
        m002,
        m110,
        m101,
        m011,
        m210,
        m201,
        m120,
        m021,
        m102,
        m012,
        m220,
        m202,
        m022,
    ] = *m;
    let mut f = [0.0_f64; 19];

    // Rest population.
    f[0] = m000 - (m200 + m020 + m002) + (m220 + m202 + m022);

    // Face populations: even and odd parts along each axis.
    let x_even = m200 - m220 - m202;
    let x_odd = m100 - m120 - m102;
    f[1] = 0.5 * (x_even + x_odd);
    f[2] = 0.5 * (x_even - x_odd);
    let y_even = m020 - m220 - m022;
    let y_odd = m010 - m210 - m012;
    f[3] = 0.5 * (y_even + y_odd);
    f[4] = 0.5 * (y_even - y_odd);
    let z_even = m002 - m202 - m022;
    let z_odd = m001 - m201 - m021;
    f[5] = 0.5 * (z_even + z_odd);
    f[6] = 0.5 * (z_even - z_odd);

    // Edge populations: velocities (sa, sb) inside each coordinate plane.
    // xy plane: 7 = (+,+), 8 = (-,+), 9 = (+,-), 10 = (-,-).
    f[7] = 0.25 * (m220 + m120 + m210 + m110);
    f[8] = 0.25 * (m220 - m120 + m210 - m110);
    f[9] = 0.25 * (m220 + m120 - m210 - m110);
    f[10] = 0.25 * (m220 - m120 - m210 + m110);
    // xz plane: 11 = (x+, z+), 12 = (x-, z+), 13 = (x+, z-), 14 = (x-, z-).
    f[11] = 0.25 * (m202 + m102 + m201 + m101);
    f[12] = 0.25 * (m202 - m102 + m201 - m101);
    f[13] = 0.25 * (m202 + m102 - m201 - m101);
    f[14] = 0.25 * (m202 - m102 - m201 + m101);
    // yz plane: 15 = (y+, z+), 16 = (y-, z+), 17 = (y+, z-), 18 = (y-, z-).
    f[15] = 0.25 * (m022 + m012 + m021 + m011);
    f[16] = 0.25 * (m022 - m012 + m021 - m011);
    f[17] = 0.25 * (m022 + m012 - m021 - m011);
    f[18] = 0.25 * (m022 - m012 - m021 + m011);
    f
}

/// Reconstruct the 19 populations from their central moments about the velocity `u`.
///
/// This is the exact inverse of
/// [`compute_central_moments_cumulant`](super::compute_central_moments_cumulant) for the same
/// `u`: the central moments are shifted back to raw moments (binomial translation) and the
/// populations are recovered with [`raw_moments_to_populations_d3q19`]. The density is the
/// zeroth central moment `kappa[0]`.
pub fn central_moments_to_populations_d3q19(kappa: &[f64; 19], u: [f64; 3]) -> [f64; 19] {
    raw_moments_to_populations_d3q19(&central_to_raw_moments_d3q19(kappa, u))
}

/// Density-weighted fourth-order cumulant correction: `(a b + 2 c^2) / rho`.
#[inline]
fn pair_product_correction(a: f64, b: f64, c: f64, inv_rho: f64) -> f64 {
    (a * b + 2.0 * c * c) * inv_rho
}

/// Inverse density, or zero for a vanishing density (the quadratic cumulant corrections then
/// drop out consistently in both transform directions).
#[inline]
fn inverse_density(rho: f64) -> f64 {
    if rho.abs() > f64::MIN_POSITIVE {
        1.0 / rho
    } else {
        0.0
    }
}

/// Convert central moments to (density-weighted) cumulants.
///
/// The cumulants are those of the normalised distribution `f / rho`, multiplied by `rho`. For
/// central moments about the mean velocity (`kappa_100 = kappa_010 = kappa_001 = 0`), the
/// cumulants of order 0 to 3 equal the central moments, and the three fourth-order cumulants of
/// the D3Q19 moment set are
///
/// `C_220 = kappa_220 - (kappa_200 kappa_020 + 2 kappa_110^2) / rho`
///
/// and its `202` and `022` counterparts. The returned array uses the ordering of
/// [`MOMENT_EXPONENTS`]. The density is `kappa[0]`.
pub fn central_moments_to_cumulants_d3q19(kappa: &[f64; 19]) -> [f64; 19] {
    let inv_rho = inverse_density(kappa[0]);
    let mut c = *kappa;
    c[16] = kappa[16] - pair_product_correction(kappa[4], kappa[5], kappa[7], inv_rho);
    c[17] = kappa[17] - pair_product_correction(kappa[4], kappa[6], kappa[8], inv_rho);
    c[18] = kappa[18] - pair_product_correction(kappa[5], kappa[6], kappa[9], inv_rho);
    c
}

/// Convert (density-weighted) cumulants back to central moments.
///
/// Exact inverse of [`central_moments_to_cumulants_d3q19`]: the second-order cumulants equal the
/// second-order central moments, which restore the fourth-order central moments.
pub fn cumulants_to_central_moments_d3q19(c: &[f64; 19]) -> [f64; 19] {
    let inv_rho = inverse_density(c[0]);
    let mut kappa = *c;
    kappa[16] = c[16] + pair_product_correction(c[4], c[5], c[7], inv_rho);
    kappa[17] = c[17] + pair_product_correction(c[4], c[6], c[8], inv_rho);
    kappa[18] = c[18] + pair_product_correction(c[5], c[6], c[9], inv_rho);
    kappa
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::{
        compute_central_moments_cumulant, compute_equilibrium_d3q19, compute_macroscopic,
        compute_raw_moments_d3q19,
    };
    use crate::lattice::D3Q19_VELOCITIES;

    /// Deterministic pseudo-random populations around a moving equilibrium.
    fn random_populations(seed: u64, amplitude: f64) -> [f64; 19] {
        let mut state = seed;
        let feq = compute_equilibrium_d3q19(1.2, [0.04, -0.03, 0.02]);
        std::array::from_fn(|i| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
            feq[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
        })
    }

    fn max_abs_diff(a: &[f64; 19], b: &[f64; 19]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn moment_exponents_have_at_most_two_active_axes_and_are_indexed() {
        for (k, e) in MOMENT_EXPONENTS.iter().enumerate() {
            let active = e.iter().filter(|&&p| p > 0).count();
            assert!(active <= 2, "moment {k} = {e:?} has three active axes");
            assert_eq!(MOMENT_INDEX[e[0]][e[1]][e[2]], k);
        }
        let missing = MOMENT_INDEX
            .iter()
            .flatten()
            .flatten()
            .filter(|&&idx| idx == usize::MAX)
            .count();
        assert_eq!(
            missing, 8,
            "exactly the eight all-axes-active triples are absent"
        );
    }

    #[test]
    fn raw_moment_ordering_matches_compute_raw_moments() {
        let f = random_populations(3, 0.1);
        let m = compute_raw_moments_d3q19(&f);
        for (k, e) in MOMENT_EXPONENTS.iter().enumerate() {
            let direct: f64 = f
                .iter()
                .zip(D3Q19_VELOCITIES.iter())
                .map(|(&fi, c)| {
                    fi * f64::from(c[0]).powi(e[0] as i32)
                        * f64::from(c[1]).powi(e[1] as i32)
                        * f64::from(c[2]).powi(e[2] as i32)
                })
                .sum();
            assert!(
                (m[k] - direct).abs() < 1e-14,
                "moment {k}: {} vs {direct}",
                m[k]
            );
        }
    }

    #[test]
    fn raw_moments_round_trip_to_populations_exactly() {
        for seed in 1..=10 {
            let f = random_populations(seed, 0.3);
            let back = raw_moments_to_populations_d3q19(&compute_raw_moments_d3q19(&f));
            assert!(max_abs_diff(&f, &back) < 1e-14, "seed {seed}");
        }
        // Single-population unit vectors check every column of the inverse.
        for i in 0..19 {
            let mut f = [0.0_f64; 19];
            f[i] = 1.0;
            let back = raw_moments_to_populations_d3q19(&compute_raw_moments_d3q19(&f));
            assert!(max_abs_diff(&f, &back) < 1e-15, "unit population {i}");
        }
    }

    #[test]
    fn central_moments_match_the_direct_definition_and_round_trip() {
        let velocities = [
            [0.0, 0.0, 0.0],
            [0.05, -0.02, 0.03],
            [-0.3, 0.2, 0.1],
            [0.7, 0.7, -0.4],
        ];
        for (seed, u) in velocities.iter().enumerate() {
            let f = random_populations(seed as u64 + 1, 0.2);
            let kappa_direct = compute_central_moments_cumulant(&f, *u);
            let kappa = raw_to_central_moments_d3q19(&compute_raw_moments_d3q19(&f), *u);
            assert!(
                max_abs_diff(&kappa, &kappa_direct) < 1e-13,
                "binomial shift differs from the direct central moments for u = {u:?}"
            );
            let back = central_moments_to_populations_d3q19(&kappa_direct, *u);
            assert!(max_abs_diff(&f, &back) < 1e-13, "round trip for u = {u:?}");
            let raw = central_to_raw_moments_d3q19(&kappa_direct, *u);
            assert!(max_abs_diff(&raw, &compute_raw_moments_d3q19(&f)) < 1e-13);
        }
    }

    #[test]
    fn cumulants_round_trip_and_vanish_at_rest_equilibrium() {
        for seed in 1..=6 {
            let f = random_populations(seed, 0.25);
            let (_, u) = compute_macroscopic(&f);
            let kappa = compute_central_moments_cumulant(&f, u);
            let back =
                cumulants_to_central_moments_d3q19(&central_moments_to_cumulants_d3q19(&kappa));
            assert!(max_abs_diff(&kappa, &back) < 1e-14, "seed {seed}");
        }
        // The D3Q19 equilibrium at rest has the Maxwellian cumulants up to fourth order:
        // C_2 = rho / 3 on the diagonal and zero for every other cumulant of order >= 1.
        let rho = 1.4;
        let feq = compute_equilibrium_d3q19(rho, [0.0; 3]);
        let c =
            central_moments_to_cumulants_d3q19(&compute_central_moments_cumulant(&feq, [0.0; 3]));
        assert!((c[0] - rho).abs() < 1e-14);
        for k in [4, 5, 6] {
            assert!((c[k] - rho / 3.0).abs() < 1e-14, "C[{k}] = {}", c[k]);
        }
        for k in [1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] {
            assert!(c[k].abs() < 1e-14, "C[{k}] = {}", c[k]);
        }
    }
}
