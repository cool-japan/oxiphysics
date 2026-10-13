// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Benveniste Mori-Tanaka homogenization of aligned spheroidal inclusions.
//!
//! For an isotropic matrix `C_m`, aligned inclusions of stiffness `C_i` and
//! volume fraction `f`, and the Eshelby tensor `S` of the inclusion shape, the
//! Mori-Tanaka strain-concentration tensor is
//! `A = [I + (1 − f) S C_m⁻¹ (C_i − C_m)]⁻¹` and the effective stiffness
//! `C = C_m + f (C_i − C_m) A` (Benveniste, 1987; equivalent to
//! `A_dil [(1 − f) I + f A_dil]⁻¹` with `A_dil = [I + S C_m⁻¹ (C_i − C_m)]⁻¹`).
//! Everything is carried out in Mandel notation (see [`super::mandel`]), where
//! the tensor double contraction is the ordinary 6×6 matrix product.
//!
//! All stiffness matrices in this module describe materials that are
//! transversely isotropic about **direction 1**.

use super::mandel::{
    Mat6, add6, identity6, inverse6, is_positive_definite, mul6, scale6, sub6, symmetrize6,
};
use crate::error::{Error, Result};

/// Effective engineering constants of a material that is transversely
/// isotropic about direction 1 (the filler's axis of revolution).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransverselyIsotropicModuli {
    /// Young's modulus along the axis of revolution E_a \[Pa\].
    pub e_axis: f64,
    /// Young's modulus perpendicular to the axis of revolution E_p \[Pa\].
    pub e_perp: f64,
    /// Shear modulus G_ap of the planes that contain the axis of revolution \[Pa\].
    pub g_axial: f64,
    /// Shear modulus G_pp in the plane perpendicular to the axis of revolution \[Pa\].
    pub g_perp: f64,
    /// Poisson's ratio ν_ap = −ε_perp/ε_axis for stress along the axis of revolution.
    pub nu_axis_perp: f64,
    /// Poisson's ratio ν_pp = −ε₃/ε₂ for stress along direction 2 in the
    /// plane perpendicular to the axis of revolution.
    pub nu_perp_perp: f64,
}

fn require_positive(param: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(Error::ParameterOutOfRange {
            param,
            value,
            min: 0.0,
            max: f64::INFINITY,
        })
    }
}

/// Mandel stiffness of a transversely isotropic material (axis of revolution =
/// direction 1) from engineering constants.
///
/// * `e_axis`, `e_perp` — Young's moduli along and perpendicular to the axis \[Pa\].
/// * `nu` — Poisson's ratio.  It is the transverse-plane ratio ν_pp and, for
///   loading along the *stiffer* principal direction, the contraction ratio of
///   the other direction (the reciprocal ratio follows from the symmetry of the
///   compliance, `ν_ap/E_a = ν_pa/E_p`, i.e. `S₁₂ = −ν / max(E_a, E_p)`).  For
///   an isotropic material (`E_a = E_p`) this is the usual Poisson's ratio.
/// * `g_axial` — shear modulus of the planes containing the axis \[Pa\].  The
///   shear modulus of the isotropy plane is `E_p / (2 (1 + ν))`.
///
/// # Errors
/// [`Error::ParameterOutOfRange`] for non-positive moduli, `ν` outside
/// `(-1, 1)`, or constants that do not give a positive-definite stiffness
/// (thermodynamically inadmissible combinations).
pub fn transversely_isotropic_stiffness(
    e_axis: f64,
    e_perp: f64,
    nu: f64,
    g_axial: f64,
) -> Result<Mat6> {
    require_positive("e_axis", e_axis)?;
    require_positive("e_perp", e_perp)?;
    if !(nu.is_finite() && nu > -1.0 && nu < 1.0) {
        return Err(Error::ParameterOutOfRange {
            param: "poisson_ratio",
            value: nu,
            min: -1.0,
            max: 1.0,
        });
    }
    require_positive("g_axial", g_axial)?;
    let e_max = e_axis.max(e_perp);
    let mut s = [[0.0; 6]; 6];
    s[0][0] = 1.0 / e_axis;
    s[1][1] = 1.0 / e_perp;
    s[2][2] = 1.0 / e_perp;
    for &(i, j) in &[(0, 1), (1, 0), (0, 2), (2, 0)] {
        s[i][j] = -nu / e_max;
    }
    s[1][2] = -nu / e_perp;
    s[2][1] = -nu / e_perp;
    s[3][3] = (1.0 + nu) / e_perp;
    s[4][4] = 0.5 / g_axial;
    s[5][5] = 0.5 / g_axial;
    if !is_positive_definite(&s) {
        return Err(Error::ParameterOutOfRange {
            param: "poisson_ratio",
            value: nu,
            min: -1.0,
            max: 1.0,
        });
    }
    inverse6(&s).ok_or_else(|| Error::Singularity {
        context: "transversely isotropic compliance",
        detail: "the compliance matrix built from the engineering constants is singular"
            .to_string(),
    })
}

/// Benveniste Mori-Tanaka effective stiffness (Mandel notation).
///
/// * `c_matrix` — stiffness of the matrix.
/// * `c_inclusion` — stiffness of the aligned inclusions.
/// * `eshelby` — Eshelby tensor of the common inclusion shape in the matrix,
///   expressed in the same frame (see
///   [`super::eshelby_spheroid::spheroid_eshelby_mandel`]).
/// * `volume_fraction` — inclusion volume fraction `f ∈ [0, 1]`.
///
/// The result is `C_m` for `f = 0` and `C_i` for `f = 1`, and is symmetrised
/// to remove round-off asymmetry.
///
/// # Errors
/// [`Error::ParameterOutOfRange`] for `f` outside `[0, 1]`;
/// [`Error::Singularity`] if the matrix stiffness or the concentration
/// operator cannot be inverted.
pub fn mori_tanaka_stiffness(
    c_matrix: &Mat6,
    c_inclusion: &Mat6,
    eshelby: &Mat6,
    volume_fraction: f64,
) -> Result<Mat6> {
    if !(0.0..=1.0).contains(&volume_fraction) {
        return Err(Error::ParameterOutOfRange {
            param: "volume_fraction",
            value: volume_fraction,
            min: 0.0,
            max: 1.0,
        });
    }
    let singular = |what: &str| Error::Singularity {
        context: "Mori-Tanaka homogenization",
        detail: format!("{what} is singular"),
    };
    let c_matrix_inv = inverse6(c_matrix).ok_or_else(|| singular("the matrix stiffness"))?;
    let contrast = sub6(c_inclusion, c_matrix);
    let coupling = mul6(eshelby, &mul6(&c_matrix_inv, &contrast));
    let operator = add6(&identity6(), &scale6(&coupling, 1.0 - volume_fraction));
    let concentration =
        inverse6(&operator).ok_or_else(|| singular("the strain-concentration operator"))?;
    let correction = mul6(&contrast, &concentration);
    Ok(symmetrize6(&add6(
        c_matrix,
        &scale6(&correction, volume_fraction),
    )))
}

/// Engineering constants of a transversely isotropic Mandel stiffness (axis of
/// revolution = direction 1), obtained from its compliance.  Round-off
/// differences between the two equivalent transverse directions are averaged.
///
/// # Errors
/// [`Error::Singularity`] if the stiffness cannot be inverted.
pub fn transversely_isotropic_moduli(stiffness: &Mat6) -> Result<TransverselyIsotropicModuli> {
    let s = inverse6(stiffness).ok_or_else(|| Error::Singularity {
        context: "transversely isotropic moduli",
        detail: "the effective stiffness is singular".to_string(),
    })?;
    let s_axis = s[0][0];
    let s_perp = 0.5 * (s[1][1] + s[2][2]);
    let s_cross = 0.25 * (s[0][1] + s[1][0] + s[0][2] + s[2][0]);
    let s_plane = 0.5 * (s[1][2] + s[2][1]);
    Ok(TransverselyIsotropicModuli {
        e_axis: 1.0 / s_axis,
        e_perp: 1.0 / s_perp,
        g_axial: 1.0 / (s[4][4] + s[5][5]),
        g_perp: 0.5 / s[3][3],
        nu_axis_perp: -s_cross / s_axis,
        nu_perp_perp: -s_plane / s_perp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nanocomposites::eshelby_spheroid::spheroid_eshelby_mandel;
    use crate::nanocomposites::mandel::isotropic_stiffness;

    const E_M: f64 = 3.5e9;
    const NU_M: f64 = 0.35;
    const E_F: f64 = 1.0e12;
    const NU_F: f64 = 0.2;

    fn matrix() -> Mat6 {
        isotropic_stiffness(E_M, NU_M)
    }

    fn isotropic_filler() -> Mat6 {
        isotropic_stiffness(E_F, NU_F)
    }

    fn isotropic_filler_ti() -> Mat6 {
        transversely_isotropic_stiffness(E_F, E_F, NU_F, E_F / (2.0 * (1.0 + NU_F))).unwrap()
    }

    fn bulk_and_shear(c: &Mat6) -> (f64, f64) {
        let k = c
            .iter()
            .take(3)
            .flat_map(|row| row.iter().take(3))
            .sum::<f64>()
            / 9.0;
        let g = (c[3][3] + c[4][4] + c[5][5]) / 6.0;
        (k, g)
    }

    fn relative_error(value: f64, reference: f64) -> f64 {
        (value - reference).abs() / reference.abs()
    }

    fn assert_close(label: &str, value: f64, reference: f64, tol: f64) {
        assert!(
            relative_error(value, reference) < tol,
            "{label}: {value} vs reference {reference}"
        );
    }

    #[test]
    fn spherical_inclusions_reproduce_the_closed_form_hashin_shtrikman_moduli() {
        // The Mori-Tanaka estimate for isotropic spheres coincides with the
        // Hashin-Shtrikman bulk/shear expressions of the inclusion phase.
        let km = E_M / (3.0 * (1.0 - 2.0 * NU_M));
        let gm = E_M / (2.0 * (1.0 + NU_M));
        let kf = E_F / (3.0 * (1.0 - 2.0 * NU_F));
        let gf = E_F / (2.0 * (1.0 + NU_F));
        let s = spheroid_eshelby_mandel(NU_M, 1.0).unwrap();
        for &f in &[0.05, 0.2, 0.5, 0.8] {
            let c = mori_tanaka_stiffness(&matrix(), &isotropic_filler(), &s, f).unwrap();
            let (k, g) = bulk_and_shear(&c);
            let k_hs = km + f / (1.0 / (kf - km) + 3.0 * (1.0 - f) / (3.0 * km + 4.0 * gm));
            let g_hs = gm
                + f / (1.0 / (gf - gm)
                    + 6.0 * (1.0 - f) * (km + 2.0 * gm) / (5.0 * gm * (3.0 * km + 4.0 * gm)));
            assert_close(&format!("K at f = {f}"), k, k_hs, 1e-12);
            assert_close(&format!("G at f = {f}"), g, g_hs, 1e-12);
            // Spheres leave the composite isotropic.
            let m = transversely_isotropic_moduli(&c).unwrap();
            assert_close("E_a vs E_p", m.e_axis, m.e_perp, 1e-12);
            assert_close("G_ap vs G_pp", m.g_axial, m.g_perp, 1e-12);
        }
    }

    #[test]
    fn volume_fraction_limits_recover_the_constituents() {
        let cm = matrix();
        let cf = transversely_isotropic_stiffness(230.0e9, 15.0e9, 0.2, 15.0e9).unwrap();
        for &alpha in &[1.0, 20.0, 500.0, 0.05] {
            let s = spheroid_eshelby_mandel(NU_M, alpha).unwrap();
            let c0 = mori_tanaka_stiffness(&cm, &cf, &s, 0.0).unwrap();
            let c1 = mori_tanaka_stiffness(&cm, &cf, &s, 1.0).unwrap();
            for (row_c0, row_m) in c0.iter().zip(cm.iter()) {
                for (a, b) in row_c0.iter().zip(row_m.iter()) {
                    assert!((a - b).abs() <= 1e-13 * 1.0e10, "f = 0: {a} vs {b}");
                }
            }
            for (row_c1, row_f) in c1.iter().zip(cf.iter()) {
                for (a, b) in row_c1.iter().zip(row_f.iter()) {
                    assert!((a - b).abs() <= 1e-12 * 1.0e12, "f = 1: {a} vs {b}");
                }
            }
        }
    }

    #[test]
    fn continuous_fibres_satisfy_hills_exact_relations() {
        // For aligned continuous fibres the plane-strain bulk modulus equals the
        // composite-cylinder value and the axial modulus follows Hill's relation
        //   E_a = f E_f + (1-f) E_m + 4 (nu_f - nu_m)^2 f (1-f) / (f/k_m + (1-f)/k_f + 1/G_m)
        // with the plane-strain bulk moduli k = K + G/3 of the constituents.
        let gm = E_M / (2.0 * (1.0 + NU_M));
        let gf = E_F / (2.0 * (1.0 + NU_F));
        let km = E_M / (3.0 * (1.0 - 2.0 * NU_M)) + gm / 3.0;
        let kf = E_F / (3.0 * (1.0 - 2.0 * NU_F)) + gf / 3.0;
        let s = spheroid_eshelby_mandel(NU_M, f64::INFINITY).unwrap();
        for &f in &[0.05, 0.3, 0.6] {
            let c = mori_tanaka_stiffness(&matrix(), &isotropic_filler(), &s, f).unwrap();
            let m = transversely_isotropic_moduli(&c).unwrap();
            let k_cca = km + f / (1.0 / (kf - km) + (1.0 - f) / (km + gm));
            let k_eff = 0.5 * (c[1][1] + c[1][2]);
            assert_close(&format!("plane-strain k at f = {f}"), k_eff, k_cca, 1e-10);
            let e_hill = f * E_F
                + (1.0 - f) * E_M
                + 4.0 * (NU_F - NU_M).powi(2) * f * (1.0 - f)
                    / (f / km + (1.0 - f) / kf + 1.0 / gm);
            assert_close(&format!("E_a at f = {f}"), m.e_axis, e_hill, 1e-9);
        }
    }

    // Reference engineering constants (E_a, E_p, nu_ap, nu_pp, G_ap, G_pp) from an independent
    // implementation: full 4th-order tensor algebra with an Eshelby tensor obtained by
    // 40-digit numerical integration of Mura's shape integrals.  Matrix: E = 3.5 GPa, nu = 0.35.
    struct Reference {
        label: &'static str,
        e_axis: f64,
        e_perp: f64,
        nu_axis_perp: f64,
        nu_perp_perp: f64,
        g_axial: f64,
        g_perp: f64,
    }

    fn check_against(reference: &Reference, m: &TransverselyIsotropicModuli) {
        let l = reference.label;
        assert_close(&format!("{l}: E_a"), m.e_axis, reference.e_axis, 1e-9);
        assert_close(&format!("{l}: E_p"), m.e_perp, reference.e_perp, 1e-9);
        assert_close(
            &format!("{l}: nu_ap"),
            m.nu_axis_perp,
            reference.nu_axis_perp,
            1e-9,
        );
        assert_close(
            &format!("{l}: nu_pp"),
            m.nu_perp_perp,
            reference.nu_perp_perp,
            1e-9,
        );
        assert_close(&format!("{l}: G_ap"), m.g_axial, reference.g_axial, 1e-9);
        assert_close(&format!("{l}: G_pp"), m.g_perp, reference.g_perp, 1e-9);
    }

    fn effective(cf: &Mat6, alpha: f64, f: f64) -> TransverselyIsotropicModuli {
        let s = spheroid_eshelby_mandel(NU_M, alpha).unwrap();
        let c = mori_tanaka_stiffness(&matrix(), cf, &s, f).unwrap();
        transversely_isotropic_moduli(&c).unwrap()
    }

    #[test]
    fn isotropic_filler_matches_independent_reference_values() {
        let cf = isotropic_filler_ti();
        let cases = [
            (
                100.0,
                0.05,
                Reference {
                    label: "alpha 100, f 0.05",
                    e_axis: 42249265665.269,
                    e_perp: 4260790429.5468,
                    nu_axis_perp: 0.34287141917260,
                    nu_perp_perp: 0.51457488387143,
                    g_axial: 1431957505.0198,
                    g_perp: 1406596159.3975,
                },
            ),
            (
                10.0,
                0.05,
                Reference {
                    label: "alpha 10, f 0.05",
                    e_axis: 7288841755.1660,
                    e_perp: 4023474056.4560,
                    nu_axis_perp: 0.34981620496965,
                    nu_perp_perp: 0.42904727318682,
                    g_axial: 1435670914.0580,
                    g_perp: 1407747011.5749,
                },
            ),
            (
                1.0,
                0.05,
                Reference {
                    label: "sphere, f 0.05",
                    e_axis: 3882576058.6431,
                    e_perp: 3882576058.6431,
                    nu_axis_perp: 0.34524418120763,
                    nu_perp_perp: 0.34524418120763,
                    g_axial: 1443074838.3382,
                    g_perp: 1443074838.3382,
                },
            ),
            (
                1000.0,
                0.01,
                Reference {
                    label: "alpha 1000, f 0.01",
                    e_axis: 13419208914.607,
                    e_perp: 3907778105.4184,
                    nu_axis_perp: 0.34807020796117,
                    nu_perp_perp: 0.48306657329818,
                    g_axial: 1322320343.5185,
                    g_perp: 1317465505.5194,
                },
            ),
            (
                0.01,
                0.02,
                Reference {
                    label: "oblate alpha 0.01, f 0.02",
                    e_axis: 4295082226.9394,
                    e_perp: 7460243731.7955,
                    nu_axis_perp: 0.23263289096669,
                    nu_perp_perp: 0.25001034805574,
                    g_axial: 1323191088.8859,
                    g_perp: 2984072789.2369,
                },
            ),
            (
                0.1,
                0.05,
                Reference {
                    label: "oblate alpha 0.1, f 0.05",
                    e_axis: 3973999027.4942,
                    e_perp: 4943114746.7554,
                    nu_axis_perp: 0.30468967817225,
                    nu_perp_perp: 0.29944858444959,
                    g_axial: 1377623503.2186,
                    g_perp: 1902004744.9007,
                },
            ),
        ];
        for (alpha, f, reference) in cases {
            check_against(&reference, &effective(&cf, alpha, f));
        }
    }

    #[test]
    fn anisotropic_fillers_match_independent_reference_values() {
        // Carbon-fibre-like filler (stiff along the axis), prolate alpha = 50, f = 0.3.
        let fibre = transversely_isotropic_stiffness(230.0e9, 15.0e9, 0.2, 15.0e9).unwrap();
        check_against(
            &Reference {
                label: "anisotropic fibre",
                e_axis: 62202830144.875,
                e_perp: 5446135641.9319,
                nu_axis_perp: 0.30790696537860,
                nu_perp_perp: 0.46547988436441,
                g_axial: 2172645476.3118,
                g_perp: 1858140701.9087,
            },
            &effective(&fibre, 50.0, 0.3),
        );
        // Graphite-like platelet (stiff in the plane, soft through the thickness),
        // oblate alpha = 1/200, f = 0.03.
        let platelet = transversely_isotropic_stiffness(36.0e9, 1.0e12, 0.165, 4.0e9).unwrap();
        check_against(
            &Reference {
                label: "anisotropic platelet",
                e_axis: 4812204437.7196,
                e_perp: 13371219335.181,
                nu_axis_perp: 0.15377112668225,
                nu_perp_perp: 0.2038357944038,
                g_axial: 1323304707.3809,
                g_perp: 5553589367.1459,
            },
            &effective(&platelet, 1.0 / 200.0, 0.03),
        );
    }

    #[test]
    fn effective_stiffness_lies_between_the_voigt_and_reuss_bounds() {
        // C_Reuss <= C_MT <= C_Voigt in the sense of quadratic forms.
        let cm = matrix();
        let fibre = transversely_isotropic_stiffness(230.0e9, 15.0e9, 0.2, 15.0e9).unwrap();
        for filler in [isotropic_filler_ti(), fibre] {
            let filler_compliance = inverse6(&filler).unwrap();
            let matrix_compliance = inverse6(&cm).unwrap();
            for &alpha in &[1.0, 8.0, 100.0, 0.05] {
                let s = spheroid_eshelby_mandel(NU_M, alpha).unwrap();
                for &f in &[0.02, 0.1, 0.4, 0.8] {
                    let c = mori_tanaka_stiffness(&cm, &filler, &s, f).unwrap();
                    let voigt = add6(&scale6(&cm, 1.0 - f), &scale6(&filler, f));
                    let reuss = inverse6(&add6(
                        &scale6(&matrix_compliance, 1.0 - f),
                        &scale6(&filler_compliance, f),
                    ))
                    .unwrap();
                    let margin = scale6(&identity6(), 1.0e-12 * 1.0e12);
                    let upper = add6(&sub6(&voigt, &c), &margin);
                    let lower = add6(&sub6(&c, &reuss), &margin);
                    assert!(
                        is_positive_definite(&upper),
                        "C_MT above Voigt at α = {alpha}, f = {f}"
                    );
                    assert!(
                        is_positive_definite(&lower),
                        "C_MT below Reuss at α = {alpha}, f = {f}"
                    );
                }
            }
        }
    }

    #[test]
    fn stiffness_builder_round_trips_the_engineering_constants() {
        // Stiff axis: nu is the contraction ratio for loading along the axis.
        let c = transversely_isotropic_stiffness(230.0e9, 15.0e9, 0.2, 15.0e9).unwrap();
        let m = transversely_isotropic_moduli(&c).unwrap();
        assert_close("E_a", m.e_axis, 230.0e9, 1e-12);
        assert_close("E_p", m.e_perp, 15.0e9, 1e-12);
        assert_close("G_ap", m.g_axial, 15.0e9, 1e-12);
        assert_close("G_pp", m.g_perp, 15.0e9 / 2.4, 1e-12);
        assert_close("nu_ap", m.nu_axis_perp, 0.2, 1e-12);
        assert_close("nu_pp", m.nu_perp_perp, 0.2, 1e-12);
        // Stiff plane: nu is the contraction ratio for loading in the plane, so the
        // axial ratio follows from the compliance symmetry nu_ap/E_a = nu_pa/E_p.
        let c = transversely_isotropic_stiffness(36.0e9, 1.0e12, 0.165, 4.0e9).unwrap();
        let m = transversely_isotropic_moduli(&c).unwrap();
        assert_close("E_a", m.e_axis, 36.0e9, 1e-12);
        assert_close("E_p", m.e_perp, 1.0e12, 1e-12);
        assert_close("nu_ap", m.nu_axis_perp, 0.165 * 36.0e9 / 1.0e12, 1e-12);
        assert_close("G_pp", m.g_perp, 1.0e12 / (2.0 * 1.165), 1e-12);
        // Isotropic constants give an isotropic stiffness.
        let iso =
            transversely_isotropic_stiffness(E_M, E_M, NU_M, E_M / (2.0 * (1.0 + NU_M))).unwrap();
        for (row_a, row_b) in iso.iter().zip(matrix().iter()) {
            for (a, b) in row_a.iter().zip(row_b.iter()) {
                assert!((a - b).abs() <= 1e-12 * 1.0e10);
            }
        }
    }

    #[test]
    fn invalid_constants_and_arguments_are_rejected() {
        assert!(transversely_isotropic_stiffness(-1.0, 1.0e9, 0.3, 1.0e9).is_err());
        assert!(transversely_isotropic_stiffness(1.0e9, 0.0, 0.3, 1.0e9).is_err());
        assert!(transversely_isotropic_stiffness(1.0e9, 1.0e9, 0.3, f64::NAN).is_err());
        assert!(transversely_isotropic_stiffness(1.0e9, 1.0e9, 1.0, 1.0e9).is_err());
        // nu = 0.6 is outside the positive-definite range of an isotropic material.
        let err = transversely_isotropic_stiffness(1.0e9, 1.0e9, 0.6, 1.0e9 / 3.2).unwrap_err();
        assert!(matches!(err, Error::ParameterOutOfRange { .. }), "{err}");
        let s = spheroid_eshelby_mandel(NU_M, 5.0).unwrap();
        assert!(mori_tanaka_stiffness(&matrix(), &isotropic_filler(), &s, -0.1).is_err());
        assert!(mori_tanaka_stiffness(&matrix(), &isotropic_filler(), &s, 1.1).is_err());
        assert!(mori_tanaka_stiffness(&matrix(), &isotropic_filler(), &s, f64::NAN).is_err());
        let zero = [[0.0; 6]; 6];
        assert!(matches!(
            mori_tanaka_stiffness(&zero, &isotropic_filler(), &s, 0.1),
            Err(Error::Singularity { .. })
        ));
        assert!(matches!(
            transversely_isotropic_moduli(&zero),
            Err(Error::Singularity { .. })
        ));
    }
}
