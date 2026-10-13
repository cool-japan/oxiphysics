// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eshelby tensor of a spheroidal inclusion in an isotropic matrix.
//!
//! For an ellipsoidal inclusion with semi-axes `a₁, a₂, a₃` in an infinite
//! isotropic matrix of Poisson's ratio `ν`, Mura (1987, eq. 11.17) expresses the
//! Eshelby tensor through the shape integrals
//! `I_i = 2π a₁a₂a₃ ∫₀^∞ ds / ((a_i² + s) Δ(s))` and
//! `I_ij = 2π a₁a₂a₃ ∫₀^∞ ds / ((a_i² + s)(a_j² + s) Δ(s))` with
//! `Δ(s) = √((a₁² + s)(a₂² + s)(a₃² + s))`.  For a spheroid (`a₁ = a₂`, aspect
//! ratio `α = a₃/a₁`) every shape integral follows from the two scalars
//!
//! * `g(α)`, with `I₁ = I₂ = 2π g` and `I₃ = 4π − 4π g`, and
//! * `h(α) = (3g − 2)/(α² − 1)`, with `I₁₃ = I₂₃ = 2π h`,
//!
//! through the identities `I₁₁ = (4π − I₁₃)/4` and
//! `α² I₃₃ = (4π − 2α² I₁₃)/3`.  `g` is the familiar closed form
//! `α (α√(α²−1) − acosh α)/(α²−1)^(3/2)` (prolate, `α > 1`) or
//! `α (acos α − α√(1−α²))/(1−α²)^(3/2)` (oblate, `α < 1`).  Close to the sphere
//! (`|α² − 1| < 1/4`) those expressions suffer catastrophic cancellation, so
//! `g` and `h` are evaluated from their convergent Taylor series in
//! `x = α² − 1` instead, which also yields the exact spherical tensor at
//! `α = 1`.
//!
//! The tensor is returned in Mandel notation (see [`super::mandel`]) in the frame
//! where the axis of revolution of the spheroid is direction 1.

use super::mandel::Mat6;
use crate::error::{Error, Result};
use std::f64::consts::PI;

/// Aspect ratios above this value (and their reciprocals below it) are
/// numerically indistinguishable from an infinite cylinder (respectively an
/// infinitely thin disc), so they are clamped to the limit.
pub const MAX_ASPECT_RATIO: f64 = 1.0e8;

/// `|α² − 1|` below which the Taylor series of `g` and `h` replaces the
/// closed forms.
const SERIES_RADIUS: f64 = 0.25;

/// Number of Taylor coefficients (orders `x⁰ … x^SERIES_ORDER`); at
/// `|x| < 1/4` the truncation error is below `1e-24`.
const SERIES_ORDER: usize = 40;

/// Evaluate a polynomial with the given coefficients (lowest order first).
fn horner(coefficients: &[f64], x: f64) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |acc, &coefficient| acc * x + coefficient)
}

/// Shape scalars `(g, h)` of a spheroid with aspect ratio `alpha`
/// (symmetry-axis semi-length over equatorial semi-length).
fn spheroid_shape_scalars(alpha: f64) -> (f64, f64) {
    // x = α² − 1 without cancellation near α = 1.
    let x = (alpha - 1.0) * (alpha + 1.0);
    if x.abs() < SERIES_RADIUS {
        // g = α Σ d_k x^k with d_k = 2 (−1)^k C(2k,k) / (4^k (2k+3)); the
        // binomial series of √(1+x) has coefficients s_j = C(1/2, j), and
        // h = (3g − 2)/x = 3 Σ_{n≥1} p_n x^(n−1) with p_n = Σ_k d_k s_(n−k).
        let mut d = [0.0_f64; SERIES_ORDER + 1];
        let mut s = [0.0_f64; SERIES_ORDER + 1];
        let mut central = 1.0_f64;
        let mut binomial_half = 1.0_f64;
        let mut sign = 1.0_f64;
        for (k, (d_k, s_k)) in d.iter_mut().zip(s.iter_mut()).enumerate() {
            *d_k = 2.0 * sign * central / (2 * k + 3) as f64;
            *s_k = binomial_half;
            central *= (2 * k + 1) as f64 / (2 * k + 2) as f64;
            binomial_half *= (0.5 - k as f64) / (k + 1) as f64;
            sign = -sign;
        }
        let g = alpha * horner(&d, x);
        let mut p = [0.0_f64; SERIES_ORDER + 1];
        for (n, p_n) in p.iter_mut().enumerate().skip(1) {
            *p_n = d
                .iter()
                .take(n + 1)
                .zip(s[..=n].iter().rev())
                .map(|(d_k, s_nk)| d_k * s_nk)
                .sum();
        }
        let h = 3.0 * horner(&p[1..], x);
        return (g, h);
    }
    let g = if alpha > 1.0 {
        let beta = x.sqrt();
        alpha * (alpha * beta - beta.asinh()) / (beta * beta * beta)
    } else {
        let beta = ((1.0 - alpha) * (1.0 + alpha)).sqrt();
        alpha * (alpha.acos() - alpha * beta) / (beta * beta * beta)
    };
    (g, (3.0 * g - 2.0) / x)
}

/// Eshelby tensor of a spheroidal inclusion in an isotropic matrix, in Mandel
/// notation (`[11, 22, 33, 23, 13, 12]`).
///
/// * `nu_matrix` — Poisson's ratio of the matrix, `-1 < ν ≤ 1/2`.
/// * `aspect_ratio` — semi-axis along the axis of revolution divided by the
///   equatorial semi-axis: `α > 1` prolate (fibres, needles), `α = 1` sphere,
///   `α < 1` oblate (platelets, discs).  Values beyond
///   [`MAX_ASPECT_RATIO`] (or its reciprocal) are clamped to the cylinder
///   (thin-disc) limit; `f64::INFINITY` is therefore a valid continuous-fibre
///   input.
///
/// The axis of revolution is direction 1, so the matrix has the structure
/// `S₁₁₁₁ = M[0][0]`, `S₁₁₂₂ = S₁₁₃₃ = M[0][1]`, `S₂₂₁₁ = S₃₃₁₁ = M[1][0]`,
/// `S₂₂₂₂ = S₃₃₃₃ = M[1][1]`, `S₂₂₃₃ = S₃₃₂₂ = M[1][2]` and the shear
/// diagonals `2 S₂₃₂₃ = M[3][3]`, `2 S₁₃₁₃ = M[4][4]`,
/// `2 S₁₂₁₂ = M[5][5]`.  Limits: the sphere gives
/// `S₁₁₁₁ = (7 − 5ν)/(15(1 − ν))`; an infinite fibre gives
/// `S₁₁₁₁ = 0`, `S₂₂₂₂ = (5 − 4ν)/(8(1 − ν))`, `S₂₂₁₁ = ν/(2(1 − ν))`,
/// `S₁₂₁₂ = 1/4`; an infinitely thin disc gives `S₁₁₁₁ = 1`,
/// `S₁₁₂₂ = ν/(1 − ν)`.  For every aspect ratio
/// `S_iikk = (1 + ν)/(1 − ν)` (sum over all `i, k`).
///
/// # Errors
/// [`Error::ParameterOutOfRange`] if `nu_matrix` is outside `(-1, 1/2]` or
/// `aspect_ratio` is not a positive number.
pub fn spheroid_eshelby_mandel(nu_matrix: f64, aspect_ratio: f64) -> Result<Mat6> {
    if !(nu_matrix.is_finite() && nu_matrix > -1.0 && nu_matrix <= 0.5) {
        return Err(Error::ParameterOutOfRange {
            param: "matrix_poisson_ratio",
            value: nu_matrix,
            min: -1.0,
            max: 0.5,
        });
    }
    if aspect_ratio.is_nan() || aspect_ratio <= 0.0 {
        return Err(Error::ParameterOutOfRange {
            param: "aspect_ratio",
            value: aspect_ratio,
            min: 0.0,
            max: f64::INFINITY,
        });
    }
    let alpha = aspect_ratio.clamp(1.0 / MAX_ASPECT_RATIO, MAX_ASPECT_RATIO);
    let (g, h) = spheroid_shape_scalars(alpha);
    let nu = nu_matrix;
    let c0 = 1.0 / (8.0 * PI * (1.0 - nu));
    let c1 = (1.0 - 2.0 * nu) * c0;
    let alpha_sq = alpha * alpha;
    // Shape integrals in the frame with the axis of revolution along direction 3
    // (Mura's convention): I1 = I2, I13 = I23, I11 = I22 = I12.
    let i1 = 2.0 * PI * g;
    let i3 = 4.0 * PI - 2.0 * i1;
    let i13 = 2.0 * PI * h;
    let i11 = (4.0 * PI - i13) / 4.0;
    let alpha_sq_i33 = (4.0 * PI - 2.0 * alpha_sq * i13) / 3.0;
    // Eshelby components (equatorial = e, axial = a).
    let s_eeee = 3.0 * c0 * i11 + c1 * i1;
    let s_eeee_cross = c0 * i11 - c1 * i1;
    let s_eeaa = c0 * alpha_sq * i13 - c1 * i1;
    let s_aaee = c0 * i13 - c1 * i3;
    let s_aaaa = 3.0 * c0 * alpha_sq_i33 + c1 * i3;
    let s_eeee_shear = c0 * i11 + c1 * i1;
    let s_ea_shear = 0.5 * c0 * (1.0 + alpha_sq) * i13 + 0.5 * c1 * (i1 + i3);
    // Relabel so that the axis of revolution is direction 1: 3 -> 1, 1 -> 2, 2 -> 3.
    let mut m = [[0.0; 6]; 6];
    m[0][0] = s_aaaa;
    m[0][1] = s_aaee;
    m[0][2] = s_aaee;
    m[1][0] = s_eeaa;
    m[2][0] = s_eeaa;
    m[1][1] = s_eeee;
    m[2][2] = s_eeee;
    m[1][2] = s_eeee_cross;
    m[2][1] = s_eeee_cross;
    m[3][3] = 2.0 * s_eeee_shear;
    m[4][4] = 2.0 * s_ea_shear;
    m[5][5] = 2.0 * s_ea_shear;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elastic::EshelbySphericalInclusion;

    fn trace_sum(m: &Mat6) -> f64 {
        m.iter().take(3).flat_map(|row| row.iter().take(3)).sum()
    }

    #[test]
    fn sphere_matches_the_spherical_eshelby_tensor() {
        for &nu in &[-0.2, 0.0, 0.2, 0.3, 0.35, 0.45, 0.5] {
            let m = spheroid_eshelby_mandel(nu, 1.0).unwrap();
            let sphere = EshelbySphericalInclusion::new(nu);
            for (i, row) in m.iter().enumerate().take(3) {
                assert!(
                    (row[i] - sphere.s1111()).abs() < 1e-14,
                    "S_iiii at ν = {nu}"
                );
                for (j, &entry) in row.iter().enumerate().take(3) {
                    if i != j {
                        assert!((entry - sphere.s1122()).abs() < 1e-14, "S_iijj at ν = {nu}");
                    }
                }
                assert!(
                    (m[i + 3][i + 3] - 2.0 * sphere.s1212()).abs() < 1e-14,
                    "S_ijij at ν = {nu}"
                );
            }
        }
    }

    #[test]
    fn infinite_cylinder_limit() {
        for &nu in &[0.2, 0.3, 0.35, 0.45] {
            let m = spheroid_eshelby_mandel(nu, f64::INFINITY).unwrap();
            let d = 1.0 - nu;
            assert!(m[0][0].abs() < 1e-12, "S1111 = {}", m[0][0]);
            assert!(m[0][1].abs() < 1e-12 && m[0][2].abs() < 1e-12);
            assert!((m[1][1] - (5.0 - 4.0 * nu) / (8.0 * d)).abs() < 1e-12);
            assert!((m[1][2] - (4.0 * nu - 1.0) / (8.0 * d)).abs() < 1e-12);
            assert!(
                (m[1][0] - nu / (2.0 * d)).abs() < 1e-12,
                "S2211 = {}",
                m[1][0]
            );
            assert!((m[2][0] - nu / (2.0 * d)).abs() < 1e-12);
            assert!((m[3][3] - 2.0 * (3.0 - 4.0 * nu) / (8.0 * d)).abs() < 1e-12);
            assert!((m[4][4] - 0.5).abs() < 1e-12 && (m[5][5] - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn slender_fibre_converges_to_the_cylinder() {
        let nu = 0.35;
        let m = spheroid_eshelby_mandel(nu, 1.0e4).unwrap();
        assert!(m[0][0].abs() < 1e-6, "S1111 = {}", m[0][0]);
        assert!((m[1][1] - (5.0 - 4.0 * nu) / (8.0 * (1.0 - nu))).abs() < 1e-6);
        assert!((m[1][0] - nu / (2.0 * (1.0 - nu))).abs() < 1e-6);
        assert!((m[5][5] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn thin_disc_limit() {
        for &nu in &[0.2, 0.35, 0.45] {
            let m = spheroid_eshelby_mandel(nu, 1.0e-8).unwrap();
            assert!((m[0][0] - 1.0).abs() < 1e-6, "S_nnnn = {}", m[0][0]);
            assert!(
                (m[0][1] - nu / (1.0 - nu)).abs() < 1e-6,
                "S_nnpp = {}",
                m[0][1]
            );
            assert!(m[1][0].abs() < 1e-6 && m[1][1].abs() < 1e-6 && m[1][2].abs() < 1e-6);
            assert!((m[4][4] - 1.0).abs() < 1e-6, "2 S_npnp = {}", m[4][4]);
            assert!(m[3][3].abs() < 1e-6);
        }
    }

    #[test]
    fn trace_identity_holds_for_every_aspect_ratio() {
        for &nu in &[0.0, 0.2, 0.35, 0.45] {
            let expected = (1.0 + nu) / (1.0 - nu);
            for &alpha in &[
                1e-6, 1e-3, 0.05, 0.3, 0.7, 0.9, 0.99, 0.999999, 1.0, 1.000001, 1.001, 1.1, 1.3,
                1.5, 3.0, 10.0, 100.0, 1.0e3, 1.0e6,
            ] {
                let m = spheroid_eshelby_mandel(nu, alpha).unwrap();
                assert!(
                    (trace_sum(&m) - expected).abs() < 1e-10,
                    "S_iikk = {} (expected {expected}) at ν = {nu}, α = {alpha}",
                    trace_sum(&m)
                );
            }
        }
    }

    #[test]
    fn transverse_isotropy_structure() {
        for &alpha in &[0.2, 0.9, 1.0, 1.7, 40.0] {
            let m = spheroid_eshelby_mandel(0.3, alpha).unwrap();
            assert_eq!(m[0][1], m[0][2]);
            assert_eq!(m[1][0], m[2][0]);
            assert_eq!(m[1][1], m[2][2]);
            assert_eq!(m[1][2], m[2][1]);
            assert_eq!(m[4][4], m[5][5]);
            // In-plane isotropy: S2323 = (S2222 - S2233)/2.
            assert!((m[3][3] - (m[1][1] - m[1][2])).abs() < 1e-13);
            // No coupling between normal and shear components.
            for row in m.iter().take(3) {
                assert!(row.iter().skip(3).all(|&v| v == 0.0));
            }
            for row in m.iter().skip(3) {
                assert!(row.iter().take(3).all(|&v| v == 0.0));
            }
        }
    }

    #[test]
    fn smooth_across_the_series_to_closed_form_switch() {
        // |alpha^2 - 1| = 1/4 on both the prolate and the oblate side.
        for &alpha0 in &[1.25_f64.sqrt(), 0.75_f64.sqrt()] {
            let below = spheroid_eshelby_mandel(0.35, alpha0 * (1.0 - 1e-9)).unwrap();
            let above = spheroid_eshelby_mandel(0.35, alpha0 * (1.0 + 1e-9)).unwrap();
            for (row_b, row_a) in below.iter().zip(above.iter()) {
                for (b, a) in row_b.iter().zip(row_a.iter()) {
                    assert!(
                        (a - b).abs() < 1e-8,
                        "jump {} across the switch",
                        (a - b).abs()
                    );
                }
            }
        }
    }

    #[test]
    fn continuity_through_the_sphere() {
        let sphere = spheroid_eshelby_mandel(0.35, 1.0).unwrap();
        for &alpha in &[1.0 - 1e-9, 1.0 + 1e-9, 1.0 + 1e-6, 1.0 - 1e-6] {
            let m = spheroid_eshelby_mandel(0.35, alpha).unwrap();
            for (row_m, row_s) in m.iter().zip(sphere.iter()) {
                for (a, b) in row_m.iter().zip(row_s.iter()) {
                    assert!((a - b).abs() < 1e-5, "α = {alpha}: {a} vs {b}");
                }
            }
        }
    }

    /// Reference values from numerical integration of Mura's shape integrals
    /// (adaptive quadrature at 40-digit precision), ν = 0.35.
    #[test]
    fn matches_numerically_integrated_reference_values() {
        let cases = [
            (1.0001, REF_ALPHA_1_0001),
            (1.01, REF_ALPHA_1_01),
            (1.1, REF_ALPHA_1_1),
            (1.2, REF_ALPHA_1_2),
            (0.99, REF_ALPHA_0_99),
            (0.9, REF_ALPHA_0_9),
            (0.8, REF_ALPHA_0_8),
            (1000.0, REF_ALPHA_1000),
            (10.0, REF_ALPHA_10),
            (3.0, REF_ALPHA_3),
            (1.5, REF_ALPHA_1_5),
            (0.5, REF_ALPHA_0_5),
            (0.1, REF_ALPHA_0_1),
        ];
        for (alpha, reference) in cases {
            let m = spheroid_eshelby_mandel(0.35, alpha).unwrap();
            let got = [
                m[0][0], m[0][1], m[1][0], m[1][1], m[1][2], m[3][3], m[4][4],
            ];
            for (g, r) in got.iter().zip(reference.iter()) {
                assert!((g - r).abs() < 1e-11, "α = {alpha}: {g} vs reference {r}");
            }
        }
    }

    // Order: S1111, S1122, S2211, S2222, S2233, 2 S2323, 2 S1313 (axis of revolution = 1).
    // Close to the sphere (series evaluation of the shape scalars) and just beyond it.
    const REF_ALPHA_1_0001: [f64; 7] = [
        0.538429012410212,
        0.076911649421162,
        0.076933186498167,
        0.538477800901191,
        0.076924395435912,
        0.461553405465278,
        0.461530990161065,
    ];
    const REF_ALPHA_1_01: [f64; 7] = [
        0.535222958939551,
        0.075790833214236,
        0.077930914801885,
        0.540075036095901,
        0.077053274879741,
        0.463021761216161,
        0.460802603826243,
    ];
    const REF_ALPHA_1_1: [f64; 7] = [
        0.507310958498447,
        0.066486463668688,
        0.086717943963304,
        0.553514527078978,
        0.078087124501344,
        0.475427402577634,
        0.455116292382981,
    ];
    const REF_ALPHA_1_2: [f64; 7] = [
        0.478741175445330,
        0.057761891315085,
        0.095889621294271,
        0.566449474088670,
        0.078989964040847,
        0.487459510047823,
        0.450450183164754,
    ];
    const REF_ALPHA_0_99: [f64; 7] = [
        0.541728542953403,
        0.078076713741422,
        0.075908938913651,
        0.536822105552147,
        0.076789508777617,
        0.460032596774530,
        0.462297324583886,
    ];
    const REF_ALPHA_0_9: [f64; 7] = [
        0.572455198843568,
        0.089510613792781,
        0.066500550523333,
        0.520803006070171,
        0.075419768653468,
        0.445383237416703,
        0.470277775799694,
    ];
    const REF_ALPHA_0_8: [f64; 7] = [
        0.609560059618798,
        0.104831090216430,
        0.055475432760443,
        0.499903507971528,
        0.073471477703738,
        0.426432030267791,
        0.482100447085178,
    ];
    const REF_ALPHA_1000: [f64; 7] = [
        0.000015986946476,
        -0.000001138679373,
        0.269224299045057,
        0.692306642207815,
        0.076923742414802,
        0.615382899793014,
        0.499993914422441,
    ];
    const REF_ALPHA_10: [f64; 7] = [
        0.044197783492040,
        -0.001032785279012,
        0.251813234477418,
        0.687230585026115,
        0.078351612490997,
        0.608878972535119,
        0.486231036960306,
    ];
    const REF_ALPHA_3: [f64; 7] = [
        0.211159449283273,
        0.007310873682203,
        0.188737844216653,
        0.655466037511640,
        0.081367058409406,
        0.574098979102234,
        0.451904748295607,
    ];
    const REF_ALPHA_1_5: [f64; 7] = [
        0.406149471067366,
        0.038867548099681,
        0.119920985615612,
        0.595950841334775,
        0.080647427877788,
        0.515303413456987,
        0.443322716403049,
    ];
    const REF_ALPHA_0_5: [f64; 7] = [
        0.741763950415026,
        0.176595241299771,
        0.020010397691541,
        0.407784256747508,
        0.063189667515206,
        0.344594589232302,
        0.549036473428829,
    ];
    const REF_ALPHA_0_1: [f64; 7] = [
        0.955582549986731,
        0.416120858323966,
        -0.009913365025412,
        0.131908122643237,
        0.022554647526382,
        0.109353475116855,
        0.835623864804970,
    ];

    #[test]
    fn rejects_invalid_arguments() {
        assert!(spheroid_eshelby_mandel(0.5, 2.0).is_ok());
        assert!(spheroid_eshelby_mandel(0.51, 2.0).is_err());
        assert!(spheroid_eshelby_mandel(-1.0, 2.0).is_err());
        assert!(spheroid_eshelby_mandel(f64::NAN, 2.0).is_err());
        assert!(spheroid_eshelby_mandel(0.3, 0.0).is_err());
        assert!(spheroid_eshelby_mandel(0.3, -2.0).is_err());
        assert!(spheroid_eshelby_mandel(0.3, f64::NAN).is_err());
    }

    #[test]
    fn infinite_aspect_ratio_equals_clamped_limit() {
        let inf = spheroid_eshelby_mandel(0.3, f64::INFINITY).unwrap();
        let clamped = spheroid_eshelby_mandel(0.3, MAX_ASPECT_RATIO).unwrap();
        assert_eq!(inf, clamped);
        let beyond = spheroid_eshelby_mandel(0.3, 10.0 * MAX_ASPECT_RATIO).unwrap();
        assert_eq!(beyond, clamped);
    }
}
