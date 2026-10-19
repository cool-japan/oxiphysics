// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the composite failure criteria and progressive damage models.

use super::*;

// --- PlaneStress ---

#[test]
fn test_plane_stress_zero() {
    let s = PlaneStress::zero();
    assert_eq!(s.sigma1, 0.0);
}

#[test]
fn test_plane_stress_scaled() {
    let s = PlaneStress::new(100.0, 50.0, 20.0).scaled(2.0);
    assert!((s.sigma1 - 200.0).abs() < 1e-10);
    assert!((s.tau12 - 40.0).abs() < 1e-10);
}

// --- PlyStrength ---

#[test]
fn test_carbon_epoxy_strengths_positive() {
    let st = PlyStrength::carbon_epoxy_t300();
    assert!(st.f1t > 0.0 && st.f1c > 0.0 && st.f2t > 0.0 && st.f12 > 0.0);
}

#[test]
fn test_glass_epoxy_f1t_less_than_carbon() {
    let ce = PlyStrength::carbon_epoxy_t300();
    let ge = PlyStrength::glass_epoxy();
    assert!(ge.f1t < ce.f1t);
}

// --- MaxStressCriterion ---

#[test]
fn test_max_stress_no_failure_under_threshold() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6); // well below strength
    assert!(!crit.has_failed(s));
}

#[test]
fn test_max_stress_failure_at_f1t() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
    assert!(crit.has_failed(s));
}

#[test]
fn test_max_stress_fi_longitudinal_tension() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(750.0e6, 0.0, 0.0);
    assert!((crit.fi_longitudinal(s) - 0.5).abs() < 1e-6);
}

#[test]
fn test_max_stress_fi_longitudinal_compression() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(-900.0e6, 0.0, 0.0);
    assert!((crit.fi_longitudinal(s) - 1.0).abs() < 1e-6);
}

#[test]
fn test_max_stress_fi_transverse_zero() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(100.0e6, 0.0, 0.0);
    assert!((crit.fi_transverse(s)).abs() < 1e-10);
}

#[test]
fn test_max_stress_load_factor_gt_one_safe() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
    assert!(crit.load_factor(s) > 1.0);
}

#[test]
fn test_max_stress_load_factor_le_one_failed() {
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
    assert!(crit.load_factor(s) <= 1.0);
}

// --- TsaiHillCriterion ---

#[test]
fn test_tsai_hill_no_failure_low_stress() {
    let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
    assert!(!crit.has_failed(s));
}

#[test]
fn test_tsai_hill_failure_at_strength() {
    let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(1500.0e6, 0.0, 0.0);
    assert!(crit.has_failed(s));
}

#[test]
fn test_tsai_hill_safety_margin_gt_one_safe() {
    let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(100.0e6, 5.0e6, 5.0e6);
    assert!(crit.safety_margin(s) > 1.0);
}

#[test]
fn test_tsai_hill_load_factor_positive() {
    let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s = PlaneStress::new(300.0e6, 20.0e6, 10.0e6);
    assert!(crit.load_factor(s) > 0.0);
}

#[test]
fn test_tsai_hill_biaxial_strength_positive() {
    let crit = TsaiHillCriterion::new(PlyStrength::carbon_epoxy_t300());
    let sigma1 = crit.biaxial_strength(0.0);
    assert!(sigma1 > 0.0);
}

#[test]
fn test_tsai_hill_failure_index_at_uniaxial_f1t() {
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let s = PlaneStress::new(st.f1t, 0.0, 0.0);
    let fi = crit.failure_index_sq(s);
    assert!((fi - 1.0).abs() < 1e-6, "FI at F1t should be 1: {fi}");
}

#[test]
fn test_tsai_hill_interaction_term_uses_f1_squared() {
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let (s1, s2) = (600.0e6, 20.0e6);
    let h2 = crit.failure_index_sq(PlaneStress::new(s1, s2, 0.0));
    let expected = (s1 / st.f1t).powi(2) - s1 * s2 / (st.f1t * st.f1t) + (s2 / st.f2t).powi(2);
    assert!(
        (h2 - expected).abs() < 1e-12,
        "H^2 = {h2}, expected {expected}"
    );
    // With the interaction term sigma1*sigma2/(F1*F2) the index would be
    // (0.4)^2 - 0.4*0.4 + (0.4)^2 = 0.16 instead of 0.3147.
    assert!((h2 - 0.314_666_666_7).abs() < 1e-9);
}

#[test]
fn test_tsai_hill_interaction_term_sign_in_compression() {
    // sigma1 < 0 and sigma2 > 0 -> the interaction term raises H^2.
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let (s1, s2) = (-450.0e6, 25.0e6);
    let h2 = crit.failure_index_sq(PlaneStress::new(s1, s2, 0.0));
    let expected = (s1 / st.f1c).powi(2) - s1 * s2 / (st.f1c * st.f1c) + (s2 / st.f2t).powi(2);
    assert!(
        (h2 - expected).abs() < 1e-12,
        "H^2 = {h2}, expected {expected}"
    );
    let without_interaction = (s1 / st.f1c).powi(2) + (s2 / st.f2t).powi(2);
    assert!(h2 > without_interaction);
}

#[test]
fn test_tsai_hill_envelope_points_lie_on_failure_surface() {
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let pts = tsai_hill_envelope(&st, 50);
    assert_eq!(pts.len(), 50);
    for &(s1, s2) in &pts {
        assert!(s2 >= 0.0, "tensile branch expected, got sigma2 = {s2}");
        let h2 = crit.failure_index_sq(PlaneStress::new(s1, s2, 0.0));
        assert!(
            (h2 - 1.0).abs() < 1e-9,
            "envelope point ({s1}, {s2}) has H^2 = {h2}"
        );
    }
}

#[test]
fn test_tsai_hill_biaxial_strength_uniaxial_limit_is_f1t() {
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let s1 = crit.biaxial_strength(0.0);
    assert!((s1 - st.f1t).abs() <= 1e-12 * st.f1t, "sigma1 = {s1}");
}

#[test]
fn test_tsai_hill_biaxial_strength_lies_on_failure_surface() {
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st);
    for &k in &[-0.5, -0.1, -0.02, 0.0, 0.02, 0.1, 0.5] {
        let s1 = crit.biaxial_strength(k);
        assert!(s1.is_finite() && s1 > 0.0, "k = {k}: sigma1 = {s1}");
        let h2 = crit.failure_index_sq(PlaneStress::new(s1, k * s1, 0.0));
        assert!((h2 - 1.0).abs() < 1e-9, "k = {k}: H^2 = {h2}");
    }
}

#[test]
fn test_tsai_hill_biaxial_strength_uses_compressive_f2_for_negative_ratio() {
    // For k < 0 the transverse stress is compressive, so F2c (> F2t) bounds it.
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = TsaiHillCriterion::new(st.clone());
    let k = -0.2;
    let s1 = crit.biaxial_strength(k);
    let a = (1.0 - k) / (st.f1t * st.f1t) + k * k / (st.f2c * st.f2c);
    assert!((s1 - 1.0 / a.sqrt()).abs() <= 1e-12 * s1);
}

// --- TsaiWuCriterion ---

#[test]
fn test_tsai_wu_stable_interaction_coeff() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    assert!(crit.is_stable());
}

#[test]
fn test_tsai_wu_no_failure_low_stress() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    let s = PlaneStress::new(100.0e6, 5.0e6, 5.0e6);
    assert!(!crit.has_failed(s));
}

#[test]
fn test_tsai_wu_failure_at_high_load() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    let s = PlaneStress::new(1600.0e6, 60.0e6, 80.0e6);
    assert!(crit.has_failed(s));
}

#[test]
fn test_tsai_wu_load_factor_positive() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    let s = PlaneStress::new(300.0e6, 20.0e6, 10.0e6);
    assert!(crit.load_factor(s) > 0.0);
}

#[test]
fn test_tsai_wu_linear_terms_sign() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    let [f1, _f2] = crit.linear_terms();
    // f1 = 1/F1t - 1/F1c. Since F1t > F1c typically, sign depends on values.
    assert!(f1.is_finite());
}

#[test]
fn test_tsai_wu_quadratic_terms_f11_positive() {
    let crit = TsaiWuCriterion::new(PlyStrength::carbon_epoxy_t300(), -0.5);
    let [f11, f22, f66, _f12] = crit.quadratic_terms();
    assert!(f11 > 0.0 && f22 > 0.0 && f66 > 0.0);
}

// --- PuckCriterion ---

#[test]
fn test_puck_ff_none_at_low_stress() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(100.0e6, 10.0e6, 5.0e6);
    assert_eq!(puck.ff_mode(s), PuckFfMode::None);
}

#[test]
fn test_puck_ff_tension_at_f1t() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(1600.0e6, 0.0, 0.0);
    assert_eq!(puck.ff_mode(s), PuckFfMode::Tension);
}

#[test]
fn test_puck_ff_compression_at_minus_f1c() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(-1000.0e6, 0.0, 0.0);
    assert_eq!(puck.ff_mode(s), PuckFfMode::Compression);
}

#[test]
fn test_puck_iff_mode_a_tension_transverse() {
    let puck = PuckCriterion::carbon_epoxy_default();
    // High transverse tension
    let s = PlaneStress::new(0.0, 60.0e6, 0.0);
    assert!(puck.iff_exposure_mode_a(s) > 1.0);
}

#[test]
fn test_puck_iff_mode_bc_compression() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(0.0, -250.0e6, 80.0e6);
    // High compressive transverse + shear: |sigma2 tau_c| = 2.18e16 exceeds
    // R_perp_perp_A |tau12| = 6.25e15, so the compression-dominated Mode C
    // governs and the exposure (hand computed) is 1.409 639 827 8.
    let exp = puck.iff_exposure_mode_bc(s);
    assert!(
        (exp - 1.409_639_827_801_86).abs() < 1e-9,
        "exposure = {exp}"
    );
    assert_eq!(puck.iff_mode(s), PuckIffMode::ModeC);
    assert!(puck.has_failed(s));
}

#[test]
fn test_puck_fracture_angle_zero_for_tension() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(0.0, 30.0e6, 0.0);
    assert!((puck.fracture_plane_angle(s)).abs() < 1e-10);
}

#[test]
fn test_puck_fracture_angle_positive_for_compression() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::new(0.0, -100.0e6, 0.0);
    let angle = puck.fracture_plane_angle(s);
    // Pure transverse compression: acos(1 / sqrt(2 (1 + p_perp_perp))) = 51.303 deg
    // for the default carbon/epoxy parameters.
    assert!((angle - 51.303_356_35).abs() < 1e-6, "angle = {angle}");
    // With shear added (sigma2 = -250 MPa, tau12 = 80 MPa) the plane tilts back.
    let tilted = puck.fracture_plane_angle(PlaneStress::new(0.0, -250.0e6, 80.0e6));
    assert!(
        (tilted - 48.400_456_66).abs() < 1e-6,
        "tilted angle = {tilted}"
    );
}

#[test]
fn test_puck_no_failure_at_zero_stress() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let s = PlaneStress::zero();
    assert!(!puck.has_failed(s));
}

/// Puck's physically based action-plane exposure for sigma_n <= 0 with
/// p_perp_perp / R_perp_perp_A = p_perp_par / R_perp_par = kappa:
/// f = kappa sigma_n + sqrt(kappa^2 sigma_n^2 + (tau_nt/R_perp_perp_A)^2
///     + (tau_n1/R_perp_par)^2) on the plane at angle theta.
fn puck_plane_exposure(puck: &PuckCriterion, sigma2: f64, tau21: f64, theta: f64) -> f64 {
    let prm = puck.action_plane_parameters();
    let s21 = puck.strength.f12;
    let kappa = puck.p_minus_12 / s21;
    let sigma_n = sigma2 * theta.cos().powi(2);
    let tau_nt = -sigma2 * theta.sin() * theta.cos();
    let tau_n1 = tau21 * theta.cos();
    let q = (tau_nt / prm.r_perp_perp_a).powi(2) + (tau_n1 / s21).powi(2);
    kappa * sigma_n + (kappa * kappa * sigma_n * sigma_n + q).sqrt()
}

/// Returns (theta, exposure) of the most exposed action plane.
fn puck_most_exposed_plane(puck: &PuckCriterion, sigma2: f64, tau21: f64) -> (f64, f64) {
    let n = 20_000;
    let mut best = (0.0, f64::NEG_INFINITY);
    for i in 0..=n {
        let theta = std::f64::consts::FRAC_PI_2 * i as f64 / n as f64;
        let f = puck_plane_exposure(puck, sigma2, tau21, theta);
        if f > best.1 {
            best = (theta, f);
        }
    }
    best
}

#[test]
fn test_puck_pure_transverse_compression_reaches_unity_at_f2c() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let f2c = puck.strength.f2c;
    let at_strength = PlaneStress::new(0.0, -f2c, 0.0);
    assert!((puck.iff_exposure(at_strength) - 1.0).abs() < 1e-12);
    assert!((puck.iff_exposure_mode_c(at_strength) - 1.0).abs() < 1e-12);
    // Exposure grows linearly with |sigma2| and is NOT identically zero.
    let ten_times = PlaneStress::new(0.0, -10.0 * f2c, 0.0);
    assert!((puck.iff_exposure(ten_times) - 10.0).abs() < 1e-9);
    assert!(puck.has_failed(ten_times));
    // Below strength: no failure; above: Mode C.
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, -0.99 * f2c, 0.0)),
        PuckIffMode::None
    );
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, -1.01 * f2c, 0.0)),
        PuckIffMode::ModeC
    );
}

#[test]
fn test_puck_mode_a_reaches_unity_for_tension_and_shear() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let f2t = puck.strength.f2t;
    let s21 = puck.strength.f12;
    let tension = puck.iff_exposure_mode_a(PlaneStress::new(0.0, f2t, 0.0));
    assert!(
        (tension - 1.0).abs() < 1e-12,
        "tension exposure = {tension}"
    );
    for tau in [s21, -s21] {
        let shear = puck.iff_exposure(PlaneStress::new(0.0, 0.0, tau));
        assert!((shear - 1.0).abs() < 1e-12, "shear exposure = {shear}");
    }
    // Slightly below / above the transverse tensile strength.
    assert!(puck.iff_exposure(PlaneStress::new(0.0, 0.99 * f2t, 0.0)) < 1.0);
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, 1.01 * f2t, 0.0)),
        PuckIffMode::ModeA
    );
}

#[test]
fn test_puck_action_plane_parameters_relations() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let prm = puck.action_plane_parameters();
    let s21 = puck.strength.f12;
    let f2c = puck.strength.f2c;
    let p = puck.p_minus_12;
    let direct = s21 / (2.0 * p) * ((1.0 + 2.0 * p * f2c / s21).sqrt() - 1.0);
    assert!((prm.r_perp_perp_a - direct).abs() <= 1e-12 * direct);
    assert!((prm.p_perp_perp - p * prm.r_perp_perp_a / s21).abs() < 1e-15);
    assert!(
        (f2c - 2.0 * (1.0 + prm.p_perp_perp) * prm.r_perp_perp_a).abs() <= 1e-9 * f2c,
        "F2c = 2 (1 + p) R_perp_perp_A must hold"
    );
    assert!((prm.tau21_c - s21 * (1.0 + 2.0 * prm.p_perp_perp).sqrt()).abs() <= 1e-12 * s21);
    // Frictionless limit p- -> 0: R_perp_perp_A -> F2c / 2.
    let mut frictionless = puck.clone();
    frictionless.p_minus_12 = 0.0;
    let prm0 = frictionless.action_plane_parameters();
    assert!((prm0.r_perp_perp_a - 0.5 * f2c).abs() <= 1e-12 * f2c);
    assert_eq!(prm0.p_perp_perp, 0.0);
    assert!((prm0.tau21_c - s21).abs() <= 1e-12 * s21);
}

#[test]
fn test_puck_modes_b_and_c_agree_at_transition() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let prm = puck.action_plane_parameters();
    let s = PlaneStress::new(0.0, -prm.r_perp_perp_a, prm.tau21_c);
    let fb = puck.iff_exposure_mode_b(s);
    let fc = puck.iff_exposure_mode_c(s);
    assert!((fb - 1.0).abs() < 1e-12, "Mode B at transition = {fb}");
    assert!((fc - 1.0).abs() < 1e-12, "Mode C at transition = {fc}");
    // Both formulas are continuous across the switch for a swept shear.
    let tau_switch = prm.tau21_c;
    let below = PlaneStress::new(0.0, -prm.r_perp_perp_a, tau_switch * (1.0 - 1e-9));
    let above = PlaneStress::new(0.0, -prm.r_perp_perp_a, tau_switch * (1.0 + 1e-9));
    let f_below = puck.iff_exposure_mode_bc(below);
    let f_above = puck.iff_exposure_mode_bc(above);
    assert!((f_below - f_above).abs() < 1e-6, "{f_below} vs {f_above}");
}

#[test]
fn test_puck_iff_mode_selection_b_vs_c() {
    let puck = PuckCriterion::carbon_epoxy_default();
    // Shear dominated, moderate compression: Mode B (|sigma2/tau| < R/tau_c = 0.894).
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, -40.0e6, 80.0e6)),
        PuckIffMode::ModeB
    );
    // Compression dominated: Mode C.
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, -195.0e6, 40.0e6)),
        PuckIffMode::ModeC
    );
    // Tension: Mode A.
    assert_eq!(
        puck.iff_mode(PlaneStress::new(0.0, 60.0e6, 10.0e6)),
        PuckIffMode::ModeA
    );
}

#[test]
fn test_puck_closed_forms_match_most_exposed_action_plane() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let states = [
        (-100.0e6, 0.0),
        (-200.0e6, 0.0),
        (-150.0e6, 40.0e6),
        (-120.0e6, 60.0e6),
        (-60.0e6, 60.0e6),
        (-40.0e6, 60.0e6),
        (-30.0e6, 70.0e6),
        (-10.0e6, 35.0e6),
    ];
    for &(sigma2, tau21) in &states {
        let stress = PlaneStress::new(0.0, sigma2, tau21);
        let closed = puck.iff_exposure_mode_bc(stress);
        let (theta, brute) = puck_most_exposed_plane(&puck, sigma2, tau21);
        assert!(
            (closed - brute).abs() <= 1e-6 * brute,
            "({sigma2}, {tau21}): closed form {closed} vs action-plane maximum {brute}"
        );
        let angle = puck.fracture_plane_angle(stress);
        assert!(
            (angle - theta.to_degrees()).abs() < 0.02,
            "({sigma2}, {tau21}): theta_fp {angle} deg vs most exposed plane {} deg",
            theta.to_degrees()
        );
    }
}

#[test]
fn test_puck_fracture_plane_angle_pure_compression() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let prm = puck.action_plane_parameters();
    let angle = puck.fracture_plane_angle(PlaneStress::new(0.0, -150.0e6, 0.0));
    let expected = (1.0 / (2.0 * (1.0 + prm.p_perp_perp)).sqrt())
        .acos()
        .to_degrees();
    assert!((angle - expected).abs() < 1e-9, "angle = {angle}");
    // Puck's fracture angle for transverse compression is about 53 +/- 3 deg.
    assert!(angle > 48.0 && angle < 56.0, "angle = {angle}");
}

#[test]
fn test_puck_fracture_plane_angle_falls_to_zero_at_mode_b_transition() {
    let puck = PuckCriterion::carbon_epoxy_default();
    let prm = puck.action_plane_parameters();
    let sigma2: f64 = -100.0e6;
    // Shear at which |sigma2/tau| = R_perp_perp_A / tau_c (Mode B/C boundary).
    let tau_boundary = sigma2.abs() * prm.tau21_c / prm.r_perp_perp_a;
    let mut previous = f64::INFINITY;
    for fraction in [0.0, 0.3, 0.6, 0.9, 0.999] {
        let angle =
            puck.fracture_plane_angle(PlaneStress::new(0.0, sigma2, fraction * tau_boundary));
        assert!(angle < previous, "angle must decrease with shear: {angle}");
        previous = angle;
    }
    let at_boundary =
        puck.fracture_plane_angle(PlaneStress::new(0.0, sigma2, tau_boundary * (1.0 + 1e-9)));
    assert!(at_boundary < 0.05, "angle at transition = {at_boundary}");
    // Mode B and tension: fracture plane at 0 deg.
    assert_eq!(
        puck.fracture_plane_angle(PlaneStress::new(0.0, sigma2, 2.0 * tau_boundary)),
        0.0
    );
    assert_eq!(
        puck.fracture_plane_angle(PlaneStress::new(0.0, 40.0e6, 10.0e6)),
        0.0
    );
}

// --- ProgressiveDamage ---

#[test]
fn test_ply_state_construction() {
    let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    assert!(ply.q11() > ply.q22());
    assert!(!ply.ff_failed);
}

#[test]
fn test_ply_state_q11_q22_ordering() {
    let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    assert!(
        ply.q11() > ply.q22(),
        "Q11 should be > Q22 for UD ply at 0°"
    );
}

#[test]
fn test_ply_degradation_total() {
    let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let q11_before = ply.q11();
    ply.apply_degradation(DegradationRule::Total, true, false);
    assert!(ply.q11() < q11_before);
    assert!(ply.ff_failed);
}

#[test]
fn test_ply_degradation_none_keeps_stiffness() {
    let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let q11_before = ply.q11();
    ply.apply_degradation(DegradationRule::None, false, true);
    assert!((ply.q11() - q11_before).abs() < 1e-3);
}

#[test]
fn test_progressive_damage_n_intact_all() {
    let ply1 = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let ply2 = PlyState::new(90.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let st1 = PlyStrength::carbon_epoxy_t300();
    let st2 = PlyStrength::carbon_epoxy_t300();
    let pd = ProgressiveDamage::new(vec![ply1, ply2], vec![st1, st2], DegradationRule::Total);
    assert_eq!(pd.n_intact(), 2);
}

#[test]
fn test_progressive_damage_not_all_failed_initially() {
    let ply1 = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let pd = ProgressiveDamage::new(
        vec![ply1],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    assert!(!pd.all_failed());
}

#[test]
fn test_progressive_damage_a_matrix_nonzero() {
    let ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let pd = ProgressiveDamage::new(
        vec![ply],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    let a = pd.a_matrix();
    assert!(a[0] > 0.0);
}

// --- Ply degradation (idempotent knock-down) ---

#[test]
fn test_ply_degradation_is_idempotent_and_applied_once() {
    let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let q_virgin = ply.q_virgin;
    ply.apply_degradation(DegradationRule::Total, true, false);
    ply.apply_degradation(DegradationRule::Total, true, false);
    assert!((ply.q11() - 0.01 * q_virgin[0]).abs() <= 1e-12 * q_virgin[0]);
    assert!((ply.q66() - 0.01 * q_virgin[3]).abs() <= 1e-12 * q_virgin[3]);
    assert!((ply.degradation - 0.01).abs() < 1e-15);
    // A single 0 degree ply: A11 = 0.01 * Q11 * t (not 0.01^2).
    let pd = ProgressiveDamage::new(
        vec![ply],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    let a11 = pd.a_matrix()[0];
    assert!(
        (a11 - 0.01 * q_virgin[0] * 0.125e-3).abs() <= 1e-9 * a11,
        "A11 = {a11}"
    );
}

#[test]
fn test_ply_degradation_never_heals() {
    let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    ply.apply_degradation(DegradationRule::Total, true, false);
    let q11 = ply.q11();
    // A later matrix-only or empty call must not restore stiffness.
    ply.apply_degradation(DegradationRule::Total, false, true);
    assert_eq!(ply.q11(), q11);
    ply.apply_degradation(DegradationRule::Total, false, false);
    assert_eq!(ply.q11(), q11);
    assert!(ply.ff_failed && ply.iff_failed);
    assert!((ply.degradation - 0.01).abs() < 1e-15);
}

#[test]
fn test_ply_degradation_matrix_failure_then_fibre_failure_compounds_to_minimum() {
    let mut ply = PlyState::new(0.0, 0.125e-3, 130.0e9, 10.0e9, 0.28, 5.0e9);
    ply.apply_degradation(DegradationRule::Partial, false, true);
    assert!((ply.degradation - 0.7).abs() < 1e-15);
    ply.apply_degradation(DegradationRule::Partial, true, false);
    assert!((ply.degradation - 0.5).abs() < 1e-15);
    assert!((ply.q22() - 0.5 * ply.q_virgin[1]).abs() <= 1e-9 * ply.q_virgin[1]);
}

// --- Progressive failure sweep ---

const PLY_T: f64 = 0.125e-3;

fn t300_ply(angle_deg: f64) -> PlyState {
    PlyState::new(angle_deg, PLY_T, 130.0e9, 10.0e9, 0.28, 5.0e9)
}

fn t300_laminate(angles: &[f64], rule: DegradationRule) -> ProgressiveDamage {
    let plies: Vec<PlyState> = angles.iter().map(|&a| t300_ply(a)).collect();
    let strengths = vec![PlyStrength::carbon_epoxy_t300(); angles.len()];
    ProgressiveDamage::new(plies, strengths, rule)
}

/// Independent evaluation of the reduced stiffness \[Q11, Q22, Q12, Q66\].
fn reduced_stiffness() -> [f64; 4] {
    let (e1, e2, nu12, g12) = (130.0e9, 10.0e9, 0.28, 5.0e9);
    let nu21 = nu12 * e2 / e1;
    let d = 1.0 - nu12 * nu21;
    [e1 / d, e2 / d, nu12 * e2 / d, g12]
}

/// Stress transformation from laminate axes to ply axes.
fn rotate_stress_to_ply(angle_deg: f64, sx: f64, sy: f64, txy: f64) -> PlaneStress {
    let (s, c) = angle_deg.to_radians().sin_cos();
    PlaneStress::new(
        c * c * sx + s * s * sy + 2.0 * c * s * txy,
        s * s * sx + c * c * sy - 2.0 * c * s * txy,
        -c * s * sx + c * s * sy + (c * c - s * s) * txy,
    )
}

fn assert_close(actual: f64, expected: f64, rel: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= rel * expected.abs(),
        "{what}: {actual} vs expected {expected}"
    );
}

#[test]
fn test_sweep_single_ply_fibre_failure_load_factor_is_exact() {
    let st = PlyStrength::carbon_epoxy_t300();
    // N1 = Xt * t  ->  lambda_FPF = 1 exactly (not capped, not scaled by 1/t).
    let mut pd = t300_laminate(&[0.0], DegradationRule::Total);
    let lf = pd
        .progressive_failure_sweep([st.f1t * PLY_T, 0.0, 0.0])
        .unwrap();
    assert_close(lf, 1.0, 1e-9, "lambda_FPF at N1 = Xt t");
    assert_close(pd.fpf_load, 1.0, 1e-9, "fpf_load");
    assert!(pd.plies[0].ff_failed);
    // A tenth of that load fails at lambda = 10 (the old sweep returned INFINITY).
    let mut pd = t300_laminate(&[0.0], DegradationRule::Total);
    let lf = pd
        .progressive_failure_sweep([0.1 * st.f1t * PLY_T, 0.0, 0.0])
        .unwrap();
    assert_close(lf, 10.0, 1e-9, "lambda_FPF at N1 = 0.1 Xt t");
    // Compression: N1 = -Xc * t -> lambda = 1.
    let mut pd = t300_laminate(&[0.0], DegradationRule::Total);
    let lf = pd
        .progressive_failure_sweep([-st.f1c * PLY_T, 0.0, 0.0])
        .unwrap();
    assert_close(lf, 1.0, 1e-9, "lambda_FPF at N1 = -Xc t");
}

#[test]
fn test_sweep_single_off_axis_ply_matches_stress_rotation() {
    // For a single ply the laminate stress is N / t in every layup, so the
    // ply-axis stress follows from the pure stress transformation.  This
    // exercises the full coupled A (A16, A26 != 0 off axis).
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = MaxStressCriterion::new(st);
    let loads = [
        [120.0e3, 15.0e3, 8.0e3],
        [-60.0e3, 40.0e3, -12.0e3],
        [30.0e3, -20.0e3, 25.0e3],
    ];
    for &angle in &[0.0, 30.0, 45.0, 60.0, 90.0, -30.0, -75.0] {
        for n in &loads {
            let stress = rotate_stress_to_ply(angle, n[0] / PLY_T, n[1] / PLY_T, n[2] / PLY_T);
            let expected = crit.load_factor(stress);
            let mut pd = t300_laminate(&[angle], DegradationRule::None);
            let lf = pd.progressive_failure_sweep(*n).unwrap();
            assert_close(lf, expected, 1e-9, &format!("angle {angle}, load {n:?}"));
        }
    }
}

#[test]
fn test_sweep_cross_ply_poisson_coupling_matches_clt_solution() {
    // [0/90]s under N_x: closed-form 2x2 solution including A12.
    let [q11, q22, q12, _q66] = reduced_stiffness();
    let a11 = 2.0 * PLY_T * (q11 + q22);
    let a12 = 4.0 * PLY_T * q12;
    let det = a11 * a11 - a12 * a12;
    let nx = 100.0e3;
    let eps_x = nx * a11 / det;
    let eps_y = -nx * a12 / det;
    let st = PlyStrength::carbon_epoxy_t300();
    let crit = MaxStressCriterion::new(st);
    let s0 = PlaneStress::new(q11 * eps_x + q12 * eps_y, q12 * eps_x + q22 * eps_y, 0.0);
    let s90 = PlaneStress::new(q11 * eps_y + q12 * eps_x, q12 * eps_y + q22 * eps_x, 0.0);
    let expected = crit.load_factor(s0).min(crit.load_factor(s90));
    // The transverse (90 degree) plies govern the first failure (matrix cracking).
    assert!(crit.load_factor(s90) < crit.load_factor(s0));
    let mut pd = t300_laminate(&[0.0, 90.0, 90.0, 0.0], DegradationRule::Total);
    let lf = pd.progressive_failure_sweep([nx, 0.0, 0.0]).unwrap();
    assert_close(lf, expected, 1e-9, "cross-ply lambda_FPF");
    assert_close(pd.fpf_load, expected, 1e-9, "fpf_load");
    assert!(pd.plies[1].iff_failed && pd.plies[2].iff_failed);
}

#[test]
fn test_sweep_second_stage_uses_degraded_stiffness() {
    // Degrade the 90 degree plies by hand (matrix cracking, factor 0.1), then
    // sweep: the next failure event follows from the re-assembled A.
    let [q11, q22, q12, _q66] = reduced_stiffness();
    let d = 0.1;
    let a11 = 2.0 * PLY_T * (q11 + d * q22);
    let a22 = 2.0 * PLY_T * (q22 + d * q11);
    let a12 = 2.0 * PLY_T * (1.0 + d) * q12;
    let det = a11 * a22 - a12 * a12;
    let nx = 100.0e3;
    let eps_x = nx * a22 / det;
    let eps_y = -nx * a12 / det;
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let s0 = PlaneStress::new(q11 * eps_x + q12 * eps_y, q12 * eps_x + q22 * eps_y, 0.0);
    let s90 = PlaneStress::new(q11 * eps_y + q12 * eps_x, q12 * eps_y + q22 * eps_x, 0.0);
    let fi0 = crit.failure_indices(s0);
    let fi90 = crit.failure_indices(s90);
    // Active modes: both modes of the 0 degree plies, fibre mode of the 90 degree plies.
    let expected = [fi0[0], fi0[1].max(fi0[2]), fi90[0]]
        .into_iter()
        .filter(|fi| *fi > 0.0)
        .map(|fi| 1.0 / fi)
        .fold(f64::INFINITY, f64::min);

    let mut pd = t300_laminate(&[0.0, 90.0, 90.0, 0.0], DegradationRule::Total);
    pd.plies[1].apply_degradation(DegradationRule::Total, false, true);
    pd.plies[2].apply_degradation(DegradationRule::Total, false, true);
    let lf = pd.progressive_failure_sweep([nx, 0.0, 0.0]).unwrap();
    assert_close(lf, expected, 1e-9, "second-stage lambda");
    // The undamaged first-stage value is lower than the second stage.
    let mut fresh = t300_laminate(&[0.0, 90.0, 90.0, 0.0], DegradationRule::Total);
    let first = fresh.progressive_failure_sweep([nx, 0.0, 0.0]).unwrap();
    assert!(lf > 1.5 * first, "stage 2 {lf} vs stage 1 {first}");
    // A full run continues past the first event up to (at least) stage 2.
    assert!(fresh.lpf_load >= lf * (1.0 - 1e-9));
    assert!(fresh.lpf_load >= fresh.fpf_load);
}

#[test]
fn test_sweep_matches_independent_event_driven_reference() {
    // First- and last-failure load factors of the full progressive analysis
    // (DegradationRule::Total), taken from an independent offline
    // implementation of the same event-driven model.
    let cases: Vec<(Vec<f64>, [f64; 3], f64, f64)> = vec![
        // Cross-ply: matrix cracking of the 90 plies, fibre failure of the
        // 0 plies, finally fibre failure of the (degraded) 90 plies.
        (
            vec![0.0, 90.0, 90.0, 0.0],
            [100.0e3, 0.0, 0.0],
            1.766_990_291_262_14,
            20.112_814_285_714_3,
        ),
        // Balanced angle ply: the shear tau12 = sigma_x / 2 gives exactly 70/60.
        (
            vec![45.0, -45.0, -45.0, 45.0],
            [60.0e3, 0.0, 0.0],
            7.0 / 6.0,
            1.370_481_927_710_843_8,
        ),
        (
            vec![0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0],
            [80.0e3, 10.0e3, 5.0e3],
            3.461_240_166_124_5,
            5.876_231_231_366_68,
        ),
        // Coupled shear/biaxial load on [30/-30]s: collapse cascades at one load level.
        (
            vec![30.0, -30.0, -30.0, 30.0],
            [50.0e3, 20.0e3, 10.0e3],
            3.621_952_791_430_49,
            3.621_952_791_430_49,
        ),
    ];
    for (angles, load, fpf, lpf) in cases {
        let mut pd = t300_laminate(&angles, DegradationRule::Total);
        let first = pd.progressive_failure_sweep(load).unwrap();
        assert_close(first, fpf, 1e-9, &format!("{angles:?} FPF"));
        assert_close(pd.fpf_load, fpf, 1e-9, &format!("{angles:?} fpf_load"));
        assert_close(pd.lpf_load, lpf, 1e-9, &format!("{angles:?} lpf_load"));
        assert!(pd.all_failed() && pd.n_intact() == 0);
        assert!(pd.plies.iter().all(|p| p.ff_failed));
    }
}

#[test]
fn test_sweep_partial_degradation_matches_independent_reference() {
    // Same event-driven model with DegradationRule::Partial (matrix failure
    // keeps 70 % of the stiffness, fibre failure 50 %, factors compound to the
    // minimum); reference values from the independent offline implementation.
    let cases: Vec<(Vec<f64>, [f64; 3], f64, f64)> = vec![
        (
            vec![0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0],
            [80.0e3, 10.0e3, 5.0e3],
            3.461_240_166_124_498,
            17.502_861_857_197_98,
        ),
        (
            vec![0.0, 90.0, 90.0, 0.0],
            [100.0e3, 0.0, 0.0],
            1.766_990_291_262_136,
            92.420_228_571_428_47,
        ),
        (
            vec![30.0, -30.0, -30.0, 30.0],
            [50.0e3, 20.0e3, 10.0e3],
            3.621_952_791_430_491,
            11.671_590_324_860_787,
        ),
    ];
    for (angles, load, fpf, lpf) in cases {
        let mut pd = t300_laminate(&angles, DegradationRule::Partial);
        let first = pd.progressive_failure_sweep(load).unwrap();
        assert_close(first, fpf, 1e-9, &format!("{angles:?} FPF"));
        assert_close(pd.fpf_load, fpf, 1e-9, &format!("{angles:?} fpf_load"));
        assert_close(pd.lpf_load, lpf, 1e-9, &format!("{angles:?} lpf_load"));
        // Every ply ends fibre-failed with the compounded 50 % knock-down.
        assert!(pd.plies.iter().all(|p| p.ff_failed));
        assert!(pd.plies.iter().all(|p| (p.degradation - 0.5).abs() < 1e-15));
    }
}

#[test]
fn test_sweep_without_degradation_visits_every_mode_in_load_order() {
    // With DegradationRule::None the stiffness never changes, so each mode's
    // load factor is fixed by the initial state and the sweep must report the
    // smallest as FPF and the largest *occurring* event as LPF.
    let [q11, q22, q12, _q66] = reduced_stiffness();
    let a11 = 2.0 * PLY_T * (q11 + q22);
    let a12 = 4.0 * PLY_T * q12;
    let det = a11 * a11 - a12 * a12;
    let nx = 100.0e3;
    let eps_x = nx * a11 / det;
    let eps_y = -nx * a12 / det;
    let crit = MaxStressCriterion::new(PlyStrength::carbon_epoxy_t300());
    let plies = [
        PlaneStress::new(q11 * eps_x + q12 * eps_y, q12 * eps_x + q22 * eps_y, 0.0),
        PlaneStress::new(q11 * eps_y + q12 * eps_x, q12 * eps_y + q22 * eps_x, 0.0),
    ];
    let mut first = f64::INFINITY;
    let mut last = 0.0_f64;
    for stress in plies {
        let [fi_f, fi_t, fi_s] = crit.failure_indices(stress);
        let lam_f = if fi_f > 0.0 {
            1.0 / fi_f
        } else {
            f64::INFINITY
        };
        let lam_i = if fi_t.max(fi_s) > 0.0 {
            1.0 / fi_t.max(fi_s)
        } else {
            f64::INFINITY
        };
        first = first.min(lam_f).min(lam_i);
        // Fibre fracture ends the ply; otherwise the matrix mode precedes it.
        let ply_last = if lam_f.is_finite() { lam_f } else { lam_i };
        last = last.max(ply_last);
    }
    let mut pd = t300_laminate(&[0.0, 90.0, 90.0, 0.0], DegradationRule::None);
    let lf = pd.progressive_failure_sweep([nx, 0.0, 0.0]).unwrap();
    assert_close(lf, first, 1e-9, "FPF");
    assert_close(pd.fpf_load, first, 1e-9, "fpf_load");
    assert_close(pd.lpf_load, last, 1e-9, "lpf_load");
    assert!(pd.lpf_load > pd.fpf_load);
    assert!(pd.plies.iter().all(|p| p.ff_failed));
    assert_eq!(pd.n_intact(), 0);
    assert!(pd.all_failed());
    // No stiffness change with the None rule.
    assert!(pd.plies.iter().all(|p| p.degradation == 1.0));
}

#[test]
fn test_sweep_failure_load_scales_with_thickness() {
    // Doubling every ply thickness halves the strains for the same N, so the
    // failure load factor doubles (A has units Pa*m; the old N/(A t) scaling
    // made lambda grow with t^2).
    let n = [80.0e3, 10.0e3, 5.0e3];
    let build = |scale: f64| {
        let plies: Vec<PlyState> = [0.0, 45.0, -45.0, 90.0]
            .iter()
            .map(|&a| PlyState::new(a, PLY_T * scale, 130.0e9, 10.0e9, 0.28, 5.0e9))
            .collect();
        ProgressiveDamage::new(
            plies,
            vec![PlyStrength::carbon_epoxy_t300(); 4],
            DegradationRule::Total,
        )
    };
    let l1 = build(1.0).progressive_failure_sweep(n).unwrap();
    let l2 = build(2.0).progressive_failure_sweep(n).unwrap();
    assert!(l1.is_finite() && l1 > 0.0);
    assert_close(l2, 2.0 * l1, 1e-9, "lambda(2t) = 2 lambda(t)");
}

#[test]
fn test_sweep_symmetric_plies_fail_simultaneously() {
    let mut pd = t300_laminate(&[45.0, -45.0, -45.0, 45.0], DegradationRule::Total);
    let lf = pd.progressive_failure_sweep([60.0e3, 0.0, 0.0]).unwrap();
    // sigma_x = 60e3 / 0.5e-3 = 120 MPa and tau12 = sigma_x / 2 = 60 MPa in every
    // +-45 ply (stiffness independent), so IFF occurs at 70 / 60.
    assert_close(lf, 7.0 / 6.0, 1e-9, "lambda_FPF of the [+-45]s laminate");
    // The +45 and -45 plies carry equal stress magnitudes and fail together.
    assert!(pd.plies.iter().all(|p| p.iff_failed));
    assert_close(pd.fpf_load, lf, 1e-12, "fpf_load");
    assert!(pd.lpf_load >= pd.fpf_load);
}

#[test]
fn test_sweep_zero_load_never_fails() {
    let mut pd = t300_laminate(&[0.0, 90.0], DegradationRule::Total);
    let lf = pd.progressive_failure_sweep([0.0, 0.0, 0.0]).unwrap();
    assert_eq!(lf, f64::INFINITY);
    assert_eq!(pd.fpf_load, 0.0);
    assert_eq!(pd.lpf_load, 0.0);
    assert_eq!(pd.n_intact(), 2);
}

#[test]
fn test_sweep_load_factors_describe_the_last_sweep_only() {
    let mut pd = t300_laminate(&[0.0, 90.0, 90.0, 0.0], DegradationRule::Total);
    let load = [100.0e3, 0.0, 0.0];
    let first = pd.progressive_failure_sweep(load).unwrap();
    assert!(first.is_finite() && first > 0.0);
    assert_eq!(pd.fpf_load, first);
    assert!(pd.lpf_load > pd.fpf_load);
    assert!(pd.plies.iter().all(|p| p.ff_failed));
    // Every ply is exhausted: a second sweep finds no active failure mode and
    // must not report the load factors of the first one.
    let again = pd.progressive_failure_sweep(load).unwrap();
    assert_eq!(again, f64::INFINITY);
    assert_eq!(pd.fpf_load, 0.0);
    assert_eq!(pd.lpf_load, 0.0);
}

#[test]
fn test_sweep_reports_singular_and_invalid_input() {
    // Empty laminate: A = 0 is singular.
    let mut empty = ProgressiveDamage::new(vec![], vec![], DegradationRule::Total);
    assert!(matches!(
        empty.progressive_failure_sweep([1.0, 0.0, 0.0]),
        Err(crate::Error::Singularity { .. })
    ));
    // One strength entry missing.
    let mut pd = ProgressiveDamage::new(
        vec![t300_ply(0.0), t300_ply(90.0)],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    assert!(matches!(
        pd.progressive_failure_sweep([1.0e5, 0.0, 0.0]),
        Err(crate::Error::DimensionMismatch { .. })
    ));
    // Non-finite load.
    let mut pd = t300_laminate(&[0.0], DegradationRule::Total);
    assert!(pd.progressive_failure_sweep([f64::NAN, 0.0, 0.0]).is_err());
    assert!(
        pd.progressive_failure_sweep([0.0, f64::INFINITY, 0.0])
            .is_err()
    );
}

/// The T300 strengths with the `which`-th in-plane strength
/// (F1t, F1c, F2t, F2c, F12) replaced by `value`.
fn t300_strength_with(which: usize, value: f64) -> PlyStrength {
    let mut strength = PlyStrength::carbon_epoxy_t300();
    match which {
        0 => strength.f1t = value,
        1 => strength.f1c = value,
        2 => strength.f2t = value,
        3 => strength.f2c = value,
        _ => strength.f12 = value,
    }
    strength
}

#[test]
fn test_sweep_rejects_invalid_ply_strengths() {
    // A zero (or negative / NaN / infinite) strength would make the failure
    // index non-finite, which used to be read as "mode unloaded" and would let
    // the ply carry any load without ever failing.
    let names = ["f1t", "f1c", "f2t", "f2c", "f12"];
    for (which, name) in names.iter().enumerate() {
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut pd = ProgressiveDamage::new(
                vec![t300_ply(0.0)],
                vec![t300_strength_with(which, value)],
                DegradationRule::Total,
            );
            let result = pd.progressive_failure_sweep([1.0e5, 0.0, 0.0]);
            if value.is_finite() {
                // Zero / negative: out of the admissible range, named by field.
                match result {
                    Err(crate::Error::ParameterOutOfRange { param, .. }) => {
                        assert_eq!(param, *name, "wrong parameter reported for {value}");
                    }
                    other => panic!("{name} = {value} must be out of range, got {other:?}"),
                }
            } else {
                assert!(
                    matches!(result, Err(crate::Error::General(_))),
                    "{name} = {value} must be reported as non-finite, got {result:?}"
                );
            }
            assert!(
                !pd.plies[0].ff_failed && !pd.plies[0].iff_failed,
                "{name} = {value}: a rejected analysis must leave the plies untouched"
            );
            assert_eq!(pd.fpf_load, 0.0);
            assert_eq!(pd.lpf_load, 0.0);
        }
    }
    // Valid strengths are accepted.
    let mut pd = ProgressiveDamage::new(
        vec![t300_ply(0.0)],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    assert!(pd.progressive_failure_sweep([1.0e5, 0.0, 0.0]).is_ok());
}

#[test]
fn test_sweep_rejects_non_finite_membrane_stiffness() {
    // A NaN ply thickness makes A non-finite; this must be reported rather
    // than silently producing "no failure" (lambda = infinity).
    let bad_ply = PlyState::new(0.0, f64::NAN, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let mut pd = ProgressiveDamage::new(
        vec![bad_ply, t300_ply(90.0)],
        vec![PlyStrength::carbon_epoxy_t300(); 2],
        DegradationRule::Total,
    );
    assert!(matches!(
        pd.progressive_failure_sweep([1.0e5, 0.0, 0.0]),
        Err(crate::Error::General(_))
    ));
    assert_eq!(pd.n_intact(), 2);
    // Zero thickness everywhere: A = 0 is singular, not non-finite.
    let flat = PlyState::new(0.0, 0.0, 130.0e9, 10.0e9, 0.28, 5.0e9);
    let mut pd = ProgressiveDamage::new(
        vec![flat],
        vec![PlyStrength::carbon_epoxy_t300()],
        DegradationRule::Total,
    );
    assert!(matches!(
        pd.progressive_failure_sweep([1.0e5, 0.0, 0.0]),
        Err(crate::Error::Singularity { .. })
    ));
}

#[test]
fn test_solve_3x3_solves_and_detects_singular() {
    let a = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]];
    let x_true = [1.0, -2.0, 3.0];
    let mut b = [0.0; 3];
    for (bi, row) in b.iter_mut().zip(a.iter()) {
        *bi = row
            .iter()
            .zip(x_true.iter())
            .map(|(aij, xj)| aij * xj)
            .sum();
    }
    let x = solve_3x3(a, b).unwrap();
    for (i, (xi, expected)) in x.iter().zip(x_true.iter()).enumerate() {
        assert!((xi - expected).abs() < 1e-12, "x[{i}] = {xi}");
    }
    // Rank-deficient: third row = first row + second row.
    let singular = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]];
    assert!(solve_3x3(singular, [1.0, 2.0, 3.0]).is_none());
    assert!(solve_3x3([[0.0; 3]; 3], [1.0, 0.0, 0.0]).is_none());
    // Row pivoting is required here (zero leading entry).
    let pivoting = [[0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
    let x = solve_3x3(pivoting, [5.0, 3.0, 7.0]).unwrap();
    let check = [2.0 * x[1] + x[2], x[0] + x[1], 3.0 * x[0] + x[2]];
    assert!((check[0] - 5.0).abs() < 1e-12);
    assert!((check[1] - 3.0).abs() < 1e-12);
    assert!((check[2] - 7.0).abs() < 1e-12);
    // Non-finite entries are rejected instead of propagating NaN.
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(solve_3x3(identity, [1.0, 2.0, 3.0]).is_some());
    let mut nan_matrix = identity;
    nan_matrix[1][1] = f64::NAN;
    assert!(solve_3x3(nan_matrix, [1.0, 2.0, 3.0]).is_none());
    let mut inf_matrix = identity;
    inf_matrix[0][2] = f64::INFINITY;
    assert!(solve_3x3(inf_matrix, [1.0, 2.0, 3.0]).is_none());
    assert!(solve_3x3(identity, [1.0, f64::NAN, 3.0]).is_none());
}

// --- DelamCriterion ---

#[test]
fn test_delam_no_failure_low_stress() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let s = InterlaminaStress::new(10.0e6, 20.0e6, 20.0e6);
    assert!(!dc.stress_delamination(s));
}

#[test]
fn test_delam_failure_above_z_t() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let s = InterlaminaStress::new(60.0e6, 0.0, 0.0); // σ_z > Z_t=50MPa
    assert!(dc.stress_delamination(s));
}

#[test]
fn test_delam_compression_no_mode1_contribution() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let s = InterlaminaStress::new(-100.0e6, 0.0, 0.0); // compressive
    assert!(!dc.stress_delamination(s));
}

#[test]
fn test_delam_linear_err_failure() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let err = EnergyReleaseRate::new(200.0, 400.0, 0.0); // at critical values
    assert!(dc.err_delamination(err));
}

#[test]
fn test_delam_linear_err_safe() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let err = EnergyReleaseRate::new(50.0, 50.0, 0.0);
    assert!(!dc.err_delamination(err));
}

#[test]
fn test_bk_criterion_at_critical() {
    let dc = DelamCriterion::carbon_epoxy_default();
    // Pure mode I at G1c
    let err = EnergyReleaseRate::new(200.0, 0.0, 0.0);
    let bk = dc.bk_criterion(err);
    assert!((bk - 1.0).abs() < 1e-6);
}

#[test]
fn test_bk_delamination_above_critical() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let err = EnergyReleaseRate::new(300.0, 0.0, 0.0);
    assert!(dc.bk_delamination(err));
}

#[test]
fn test_bk_no_delamination_below_critical() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let err = EnergyReleaseRate::new(50.0, 50.0, 0.0);
    assert!(!dc.bk_delamination(err));
}

#[test]
fn test_err_mode_i_ratio() {
    let err = EnergyReleaseRate::new(100.0, 50.0, 50.0);
    let r = err.mode_i_ratio();
    assert!((r - 0.5).abs() < 1e-10);
}

#[test]
fn test_err_total() {
    let err = EnergyReleaseRate::new(100.0, 200.0, 50.0);
    assert!((err.total() - 350.0).abs() < 1e-10);
}

#[test]
fn test_stress_intensity_mode1() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let k1 = dc.stress_intensity_mode1(1.0e6, 1.0e-3 / PI);
    assert!((k1 - 1.0e6 / 1000.0_f64.sqrt()).abs() < 1.0);
}

#[test]
fn test_g1c_from_k1c() {
    let g = DelamCriterion::g1c_from_k1c(1000.0, 5.0e9);
    assert!((g - 1000.0 * 1000.0 / 5.0e9).abs() < 1e-10);
}

#[test]
fn test_paris_growth_rate_positive() {
    let dc = DelamCriterion::carbon_epoxy_default();
    let rate = dc.paris_growth_rate(100.0, 200.0, 1e-10, 3.0);
    assert!(rate > 0.0);
}

// --- Helper functions ---

#[test]
fn test_transformed_q_identity_at_zero_angle() {
    let q = [130.0e9_f64, 10.0e9_f64, 3.0e9_f64, 5.0e9_f64];
    let qbar = transformed_q(&q, 0.0);
    assert!((qbar[0] - q[0]).abs() < 1.0, "Q11 not preserved at 0°");
    assert!((qbar[1] - q[1]).abs() < 1.0, "Q22 not preserved at 0°");
}

#[test]
fn test_transformed_q_90deg_swaps_11_22() {
    let q = [130.0e9_f64, 10.0e9_f64, 3.0e9_f64, 5.0e9_f64];
    let qbar0 = transformed_q(&q, 0.0);
    let qbar90 = transformed_q(&q, PI / 2.0);
    assert!(
        (qbar90[0] - qbar0[1]).abs() < 1.0,
        "Q̄11(90°) should equal Q̄22(0°)"
    );
}

#[test]
fn test_reserve_factor_positive() {
    assert!((reserve_factor(0.5).unwrap() - 2.0).abs() < 1e-10);
}

#[test]
fn test_reserve_factor_none_at_zero() {
    assert!(reserve_factor(0.0).is_none());
}

#[test]
fn test_tsai_hill_envelope_nonempty() {
    let st = PlyStrength::carbon_epoxy_t300();
    let pts = tsai_hill_envelope(&st, 10);
    assert!(!pts.is_empty());
}
