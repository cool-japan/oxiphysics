//! Regression tests for the stress-intensity-factor solutions: weight-function
//! quadrature, three-point-bend, crack at a hole and surface crack.

use super::functions::*;
use std::f64::consts::PI;

fn relative_error(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

const M1: f64 = 0.6147;
const M2: f64 = -0.1244;

// ---------------------------------------------------------------------------
// Weight-function SIF
// ---------------------------------------------------------------------------

/// Closed form of the uniform-stress weight-function SIF.
fn uniform_closed_form(stress: f64, a: f64, m1: f64, m2: f64, m3: f64) -> f64 {
    stress * 4.0 * a.sqrt() / (2.0 * PI).sqrt() * (1.0 + m1 / 2.0 + m2 / 3.0 + m3 / 4.0)
}

#[test]
fn weight_function_uniform_stress_matches_closed_form_for_any_panel_count() {
    let (stress, a) = (100.0e6, 0.01);
    for (m1, m2, m3) in [(M1, M2, 0.0), (0.6147, -0.1244, 0.05), (0.0, 0.0, 0.0)] {
        let exact = uniform_closed_form(stress, a, m1, m2, m3);
        for n_points in [0, 1, 2, 7, 100, 1000] {
            let k = sif_weight_function(&|_| stress, a, m1, m2, m3, n_points);
            assert!(
                relative_error(k, exact) < 1e-12,
                "n_points = {n_points}: K = {k}, exact = {exact}"
            );
        }
    }
}

#[test]
fn weight_function_endpoint_singularity_no_longer_dominates() {
    // The clamped trapezoid rule added an end-point term of ~4 sigma (about 20 times the
    // true K ~ 0.2 sigma for a = 0.01 m), set by the 1e-10 clamp rather than by n_points.
    let (stress, a) = (100.0e6, 0.01);
    let k = sif_weight_function(&|_| stress, a, M1, M2, 0.0, 1000);
    let edge_crack = 1.12 * stress * (PI * a).sqrt();
    assert!(
        (k - edge_crack).abs() < 0.02 * edge_crack,
        "K = {k}, 1.12 sigma sqrt(pi a) = {edge_crack}"
    );
    assert!(
        k < 0.25 * stress,
        "K = {k} must be ~0.2 sigma, not ~4 sigma"
    );
}

#[test]
fn weight_function_gauss_legendre_rule_integrates_u_powers_exactly() {
    // sigma(x) = sigma0 (1 - x/a)^k = sigma0 u^(2k): the integrand
    // u^(2k) (1 + M1 u + M2 u^2 + M3 u^3) has degree 2k + 3 <= 31 for k <= 14,
    // so one 16-point panel must reproduce  1/(2k+1) + M1/(2k+2) + M2/(2k+3) + M3/(2k+4).
    let (stress, a) = (1.0_f64, 0.04_f64);
    let (m1, m2, m3) = (0.6147_f64, -0.1244_f64, 0.05_f64);
    let prefactor = 4.0 * a.sqrt() / (2.0 * PI).sqrt();
    for k in 0..=14_i32 {
        let kf = f64::from(k);
        let exact = prefactor
            * (1.0 / (2.0 * kf + 1.0)
                + m1 / (2.0 * kf + 2.0)
                + m2 / (2.0 * kf + 3.0)
                + m3 / (2.0 * kf + 4.0));
        let profile = |x: f64| stress * (1.0 - x / a).powi(k);
        let value = sif_weight_function(&profile, a, m1, m2, m3, 1);
        assert!(
            relative_error(value, exact) < 1e-13,
            "k = {k}: {value} vs {exact}"
        );
    }
}

#[test]
fn weight_function_converges_for_a_rapidly_varying_profile() {
    // sigma(x) = sigma0 exp(-3 x/a), M = (0.6147, -0.1244, 0.05), a = 0.01 m; the reference
    // is a 40-digit adaptive quadrature of the substituted integral.
    let (stress, a) = (100.0e6, 0.01);
    let reference = 4_737_515.824_316_793;
    let profile = |x: f64| stress * (-3.0 * x / a).exp();
    let coarse = sif_weight_function(&profile, a, 0.6147, -0.1244, 0.05, 1);
    let fine = sif_weight_function(&profile, a, 0.6147, -0.1244, 0.05, 8);
    assert!(
        relative_error(coarse, reference) < 1e-9,
        "one panel: {coarse}"
    );
    assert!(
        relative_error(fine, reference) < 1e-13,
        "eight panels: {fine}"
    );
}

#[test]
fn weight_function_degenerate_crack_lengths_give_zero() {
    for a in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        assert_eq!(sif_weight_function(&|_| 1.0e6, a, M1, M2, 0.0, 10), 0.0);
    }
}

// ---------------------------------------------------------------------------
// Three-point bend (SENB)
// ---------------------------------------------------------------------------

#[test]
fn senb_sif_uses_w_to_the_three_halves() {
    // ASTM E399 SE(B): P = 10 kN, S = 0.2 m, B = 25 mm, W = 50 mm, a = 25 mm.
    let k = sif_three_point_bend(10.0e3, 0.2, 0.025, 0.05, 0.025);
    // f(0.5) = 2.6625 and K = P S / (B W^1.5) f = 19.05 MPa sqrt(m).
    let f_alpha = k * 0.025 * 0.05_f64.powf(1.5) / (10.0e3 * 0.2);
    assert!((f_alpha - 2.6625).abs() < 1e-9, "f(0.5) = {f_alpha}");
    assert!(
        (k - 19.05e6).abs() < 0.001 * 19.05e6,
        "K = {k}, expected 19.05e6 Pa sqrt(m)"
    );
}

#[test]
fn senb_sif_scales_as_w_to_the_minus_three_halves_at_fixed_alpha() {
    // Geometrically scaling W and a by 2 at fixed P, S, B lowers K by 2^1.5
    // (the old W^0.5 denominator gave 2^0.5).
    let k1 = sif_three_point_bend(8.0e3, 0.16, 0.02, 0.04, 0.016);
    let k2 = sif_three_point_bend(8.0e3, 0.16, 0.02, 0.08, 0.032);
    assert!(
        relative_error(k1 / k2, 2.0_f64.powf(1.5)) < 1e-12,
        "K1/K2 = {}",
        k1 / k2
    );
}

#[test]
fn senb_critical_crack_length_round_trips_for_several_widths() {
    for width in [0.01, 0.025, 0.05, 0.1] {
        let (force, span, thickness) = (3.0e3, 4.0 * width, 0.5 * width);
        for fraction in [0.1, 0.3, 0.5, 0.8] {
            let a_target = fraction * width;
            let k_ic = sif_three_point_bend(force, span, thickness, width, a_target);
            let a_crit = senb_critical_crack_length(force, span, thickness, width, k_ic, 200);
            let a_crit = a_crit.expect("the toughness of a known crack must bracket a root");
            assert!(
                (a_crit - a_target).abs() < 1e-9 * width,
                "W = {width}, a/W = {fraction}: a_crit = {a_crit}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Crack emanating from a circular hole
// ---------------------------------------------------------------------------

#[test]
fn hole_crack_factor_matches_the_bowie_limits() {
    let f_short = hole_crack_geometry_factor(1.0e-6);
    assert!((f_short - 3.365).abs() < 0.002, "F(a/R -> 0) = {f_short}");
    let f_unit = hole_crack_geometry_factor(1.0);
    assert!((f_unit - 1.444).abs() < 0.002, "F(a/R = 1) = {f_unit}");
    let f_long = hole_crack_geometry_factor(1.0e6);
    assert!((f_long - 1.0).abs() < 1e-3, "F(a/R -> inf) = {f_long}");
}

#[test]
fn hole_crack_factor_decreases_monotonically_from_the_stress_concentration() {
    let mut previous = f64::INFINITY;
    for exponent in -6..=6 {
        let lambda = 10.0_f64.powi(exponent);
        let factor = hole_crack_geometry_factor(lambda);
        assert!(
            factor < previous && factor > 1.0 - 1e-5,
            "a/R = {lambda}: F = {factor}"
        );
        previous = factor;
    }
    assert!(hole_crack_geometry_factor(f64::NAN).is_nan());
    // Negative ratios are clamped to the a/R -> 0 limit.
    assert_eq!(
        hole_crack_geometry_factor(-1.0),
        hole_crack_geometry_factor(0.0)
    );
}

#[test]
fn hole_crack_sif_is_about_3_4_times_the_plain_crack_for_short_cracks() {
    // The old polynomial 1 + 0.1215 L - 0.0598 L^2 (clamped to >= 1) ignored K_t = 3.
    let stress = 100.0e6;
    let radius = 1.0;
    let crack = 1.0e-5;
    let k_hole = sif_hole_edge_crack(stress, crack, radius);
    let k_plain = stress * (PI * crack).sqrt();
    let ratio = k_hole / k_plain;
    assert!((ratio - 3.364).abs() < 0.005, "K_hole / K_plain = {ratio}");
}

#[test]
fn hole_crack_sif_matches_the_geometry_factor_and_long_crack_limit() {
    let stress = 100.0e6;
    let k = sif_hole_edge_crack(stress, 0.001, 0.01); // a/R = 0.1
    let expected = stress * (PI * 0.001_f64).sqrt() * hole_crack_geometry_factor(0.1);
    assert!(relative_error(k, expected) < 1e-14);
    assert!(relative_error(k, 15_766_406.547_102_813) < 1e-12, "K = {k}");
    // A vanishing hole leaves a central crack of half-length a: K = sigma sqrt(pi a).
    let k_no_hole = sif_hole_edge_crack(stress, 0.002, 0.0);
    assert!(relative_error(k_no_hole, stress * (PI * 0.002_f64).sqrt()) < 1e-14);
    // Very long cracks approach the central-crack value as well.
    let k_long = sif_hole_edge_crack(stress, 1.0, 1.0e-6);
    assert!(relative_error(k_long, stress * PI.sqrt()) < 1e-5);
}

#[test]
fn hole_crack_sif_degenerate_inputs() {
    assert_eq!(sif_hole_edge_crack(100.0e6, 0.0, 0.01), 0.0);
    assert_eq!(sif_hole_edge_crack(100.0e6, -0.1, 0.01), 0.0);
    assert!(sif_hole_edge_crack(100.0e6, 0.01, -0.01).is_nan());
    assert!(sif_hole_edge_crack(100.0e6, 0.01, f64::NAN).is_nan());
}

// ---------------------------------------------------------------------------
// Surface crack
// ---------------------------------------------------------------------------

#[test]
fn surface_crack_shape_factor_is_2_464_for_a_semicircle() {
    let q = surface_crack_shape_factor(0.005, 0.005);
    assert!((q - 2.464).abs() < 1e-3, "Q(a = c) = {q}");
    // Q = (pi/2)^2 = 2.467 for a circle: the 1.464 coefficient reproduces it to 0.1 %.
    assert!(relative_error(q, (PI / 2.0).powi(2)) < 2e-3);
    // The old 4.593 coefficient gave 5.593.
    assert!(q < 3.0);
}

#[test]
fn surface_crack_shape_factor_is_symmetric_under_inversion_of_the_aspect_ratio() {
    let shallow = surface_crack_shape_factor(0.005, 0.01); // a/c = 0.5
    let deep = surface_crack_shape_factor(0.01, 0.005); // a/c = 2
    assert!(relative_error(shallow, 1.0 + 1.464 * 0.5_f64.powf(1.65)) < 1e-14);
    assert_eq!(shallow, deep);
    assert!(surface_crack_shape_factor(1.0e-6, 1.0) < 1.0001);
}

#[test]
fn surface_crack_sif_applies_the_deepest_point_correction() {
    let stress = 100.0e6;
    let a = 0.005;
    let k = sif_surface_crack(stress, a, a);
    let expected = 1.04 * stress * (PI * a / 2.464).sqrt();
    assert!(
        relative_error(k, expected) < 1e-12,
        "K = {k}, expected {expected}"
    );
    assert!(relative_error(k, 8_303_724.390_934_898) < 1e-9);
    // The old Q = 5.593 and missing F_s gave a SIF 0.64 times this value.
    let old_value = stress * (PI * a / 5.593).sqrt();
    assert!(k / old_value > 1.5, "ratio = {}", k / old_value);
}

#[test]
fn surface_crack_deepest_point_factor_follows_newman_raju() {
    assert!(relative_error(surface_crack_deepest_point_factor(0.005, 0.005), 1.04) < 1e-14);
    assert!(relative_error(surface_crack_deepest_point_factor(0.005, 0.01), 1.085) < 1e-14);
    assert!(relative_error(surface_crack_deepest_point_factor(0.002, 0.01), 1.112) < 1e-14);
    // a/c = 2: M1 = sqrt(0.5) (1 + 0.02) and f_phi(pi/2) = sqrt(0.5).
    let deep = surface_crack_deepest_point_factor(0.01, 0.005);
    assert!(relative_error(deep, 0.5 * 1.02) < 1e-14, "F_s = {deep}");
}

#[test]
fn surface_crack_sif_is_continuous_across_a_equals_c() {
    let stress = 100.0e6;
    let c = 0.01;
    let below = sif_surface_crack(stress, c * (1.0 - 1e-9), c);
    let at = sif_surface_crack(stress, c, c);
    let above = sif_surface_crack(stress, c * (1.0 + 1e-9), c);
    assert!(relative_error(below, at) < 1e-6, "{below} vs {at}");
    assert!(relative_error(above, at) < 1e-6, "{above} vs {at}");
}

#[test]
fn surface_crack_sif_grows_with_depth_and_degenerate_inputs_give_zero() {
    let k_shallow = sif_surface_crack(100.0e6, 0.005, 0.02);
    let k_deeper = sif_surface_crack(100.0e6, 0.008, 0.02);
    assert!(k_deeper > k_shallow);
    assert_eq!(sif_surface_crack(100.0e6, 0.0, 0.01), 0.0);
    assert_eq!(sif_surface_crack(100.0e6, 0.01, 0.0), 0.0);
    assert_eq!(sif_surface_crack(100.0e6, -0.01, 0.01), 0.0);
}
