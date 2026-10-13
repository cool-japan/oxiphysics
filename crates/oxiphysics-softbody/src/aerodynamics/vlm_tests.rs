// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`VortexLatticeMethod`].
//!
//! The expected values come from independent theory rather than from the
//! solver itself:
//!
//! * the closed-form downwash of one horseshoe vortex (Biot-Savart integrated
//!   analytically for a straight bound segment and two semi-infinite legs),
//! * the Glauert Fourier-series lifting-line solution of a rectangular wing,
//! * the lifting-line result for an elliptic wing, `2π·AR / (AR + 2)`,
//!   which a high-aspect-ratio rectangular wing must approach.

use super::helpers::gauss_solve_opt;
use super::types::VortexLatticeMethod;
use std::f64::consts::PI;

const DENSITY: f64 = 1.225;
const SPEED: f64 = 30.0;

/// Lift coefficient `L / (q S)` of `vlm` at [`SPEED`] and [`DENSITY`].
fn lift_coefficient(vlm: &VortexLatticeMethod) -> f64 {
    let lift = vlm
        .total_lift(SPEED, DENSITY)
        .expect("a regular wing must solve");
    lift / (0.5 * DENSITY * SPEED * SPEED * vlm.span * vlm.chord)
}

/// Lift-curve slope `dCL/dα` (1/rad) of a rectangular wing, measured at one
/// degree angle of attack.  The VLM lift is exactly proportional to `sin α`.
fn lift_slope(n_span: usize, n_chord: usize, span: f64, chord: f64) -> f64 {
    let alpha = 1.0_f64.to_radians();
    let vlm = VortexLatticeMethod::new(n_span, n_chord, span, chord, alpha);
    lift_coefficient(&vlm) / alpha.sin()
}

/// Glauert lifting-line lift slope (1/rad) of a rectangular wing with
/// two-dimensional section slope `2π`.
///
/// Solves `Σ A_n sin(nθ) · (4b / (a0 c) + n / sin θ) = 1` at `N` stations over
/// the half span for the odd Fourier coefficients and returns `π · AR · A_1`.
fn glauert_rectangular_lift_slope(aspect_ratio: f64) -> f64 {
    let a0 = 2.0 * PI;
    let stations = 40_usize;
    let mu = 4.0 * aspect_ratio / a0;
    let mut matrix = vec![0.0_f64; stations * stations];
    for k in 1..=stations {
        let theta = k as f64 * PI / (2.0 * stations as f64);
        for m in 0..stations {
            let order = (2 * m + 1) as f64;
            matrix[(k - 1) * stations + m] = (order * theta).sin() * (mu + order / theta.sin());
        }
    }
    let rhs = vec![1.0_f64; stations];
    let coefficients =
        gauss_solve_opt(stations, &mut matrix, &rhs).expect("the Glauert system is regular");
    PI * aspect_ratio * coefficients[0]
}

#[test]
fn collocation_points_sit_half_a_panel_behind_the_bound_vortex() {
    let (n_span, n_chord, span, chord) = (5_usize, 4_usize, 6.0_f64, 2.0_f64);
    let vlm = VortexLatticeMethod::new(n_span, n_chord, span, chord, 0.05);
    let dx = span / n_span as f64;
    let dy = chord / n_chord as f64;
    assert_eq!(vlm.n_panels(), n_span * n_chord);
    assert_eq!(vlm.panel_centers.len(), n_span * n_chord);
    assert_eq!(vlm.bound_vortex_pts.len(), n_span * n_chord);
    for i in 0..n_span {
        for j in 0..n_chord {
            let k = i * n_chord + j;
            let cp = vlm.panel_centers[k];
            let (left, right) = vlm.bound_vortex_pts[k];
            // Bound vortex at the panel quarter chord, collocation at 3/4.
            assert!(
                (left[1] - (j as f64 + 0.25) * dy).abs() < 1e-12,
                "panel {k}: bound vortex must sit at 1/4 chord"
            );
            assert!(
                (cp[1] - left[1] - 0.5 * dy).abs() < 1e-12,
                "panel {k}: collocation point must be 0.5*dy behind the bound vortex, \
                 got {}",
                cp[1] - left[1]
            );
            // Mid-span of the bound segment, in the wing plane.
            assert!((cp[0] - 0.5 * (left[0] + right[0])).abs() < 1e-12);
            assert!((right[0] - left[0] - dx).abs() < 1e-12);
            assert_eq!(cp[2], 0.0);
            assert_eq!(left[2], 0.0);
            assert_eq!(right[1], left[1], "bound vortex must run along X");
        }
    }
}

#[test]
fn single_horseshoe_matches_the_closed_form_downwash() {
    // One panel: bound vortex at 1/4 chord, collocation point at 3/4 chord,
    // i.e. a distance d = chord/2 behind the bound vortex at mid-span.
    // Downwash per unit circulation for a bound segment of half-length s and
    // two semi-infinite trailing legs, evaluated at that point:
    //   w = 1/(4π) · [ 2s / (d·√(s²+d²)) + (2/s) · (1 + d/√(s²+d²)) ].
    let (span, chord, speed) = (2.0_f64, 1.0_f64, 10.0_f64);
    let alpha = 3.0_f64.to_radians();
    let vlm = VortexLatticeMethod::new(1, 1, span, chord, alpha);
    let s = 0.5 * span;
    let d = 0.5 * chord;
    let root = (s * s + d * d).sqrt();
    let w_unit = (2.0 * s / (d * root) + (2.0 / s) * (1.0 + d / root)) / (4.0 * PI);
    // No penetration at the collocation point: V sin(α) - Γ · w = 0.
    let expected_gamma = speed * alpha.sin() / w_unit;
    let gammas = vlm
        .solve_gammas(speed, DENSITY)
        .expect("one regular panel must solve");
    assert_eq!(gammas.len(), 1);
    let rel = (gammas[0] - expected_gamma).abs() / expected_gamma;
    assert!(
        rel < 1e-4,
        "Γ = {} but the closed form gives {expected_gamma} (rel err {rel:e})",
        gammas[0]
    );
    // Kutta-Joukowski: L = ρ V Γ b.
    let lift = vlm.total_lift(speed, DENSITY).expect("lift");
    let expected_lift = DENSITY * speed * expected_gamma * span;
    assert!(((lift - expected_lift) / expected_lift).abs() < 1e-4);
}

#[test]
fn lift_changes_sign_with_angle_of_attack() {
    let alpha = 4.0_f64.to_radians();
    let up = VortexLatticeMethod::new(8, 3, 6.0, 1.0, alpha)
        .total_lift(SPEED, DENSITY)
        .expect("lift at +alpha");
    let down = VortexLatticeMethod::new(8, 3, 6.0, 1.0, -alpha)
        .total_lift(SPEED, DENSITY)
        .expect("lift at -alpha");
    let zero = VortexLatticeMethod::new(8, 3, 6.0, 1.0, 0.0)
        .total_lift(SPEED, DENSITY)
        .expect("lift at zero alpha");
    assert!(up > 0.0, "positive angle of attack must lift, got {up}");
    assert!(
        down < 0.0,
        "negative angle of attack must push down, got {down}"
    );
    assert!(
        ((up + down) / up).abs() < 1e-9,
        "the flat-plate lift must be odd in alpha: {up} vs {down}"
    );
    assert!(
        zero.abs() < 1e-9 * up,
        "no lift at zero incidence, got {zero}"
    );
}

#[test]
fn circulation_is_positive_for_positive_incidence() {
    let vlm = VortexLatticeMethod::new(6, 2, 5.0, 1.0, 3.0_f64.to_radians());
    let gammas = vlm.solve_gammas(SPEED, DENSITY).expect("solve");
    assert_eq!(gammas.len(), 12);
    assert!(
        gammas.iter().all(|&g| g > 0.0),
        "every panel of a lifting flat plate carries positive circulation: {gammas:?}"
    );
}

#[test]
fn lift_scales_with_sin_alpha_and_speed_squared() {
    let a1 = 1.0_f64.to_radians();
    let a2 = 2.0_f64.to_radians();
    let l1 = VortexLatticeMethod::new(8, 2, 6.0, 1.0, a1)
        .total_lift(SPEED, DENSITY)
        .expect("lift");
    let l2 = VortexLatticeMethod::new(8, 2, 6.0, 1.0, a2)
        .total_lift(SPEED, DENSITY)
        .expect("lift");
    // The no-penetration right-hand side is -V sin(alpha): lift is linear in
    // sin(alpha), so doubling the angle scales it by sin(2a)/sin(a) = 2 cos a.
    let expected = a2.sin() / a1.sin();
    assert!(
        ((l2 / l1) - expected).abs() < 1e-9,
        "lift ratio {} should equal sin(2a)/sin(a) = {expected}",
        l2 / l1
    );
    // Lift is proportional to V² (Γ ∝ V and L = ρ V Γ b).
    let vlm = VortexLatticeMethod::new(8, 2, 6.0, 1.0, a2);
    let slow = vlm.total_lift(10.0, DENSITY).expect("lift");
    let fast = vlm.total_lift(20.0, DENSITY).expect("lift");
    assert!(((fast / slow) - 4.0).abs() < 1e-9);
}

#[test]
fn lift_does_not_depend_on_the_chordwise_discretisation() {
    // The old implementation summed only the first n_span circulations, which
    // is wrong whenever n_chord > 1.
    let alpha = 3.0_f64.to_radians();
    let one = VortexLatticeMethod::new(16, 1, 8.0, 1.0, alpha)
        .total_lift(SPEED, DENSITY)
        .expect("lift n_chord=1");
    let four = VortexLatticeMethod::new(16, 4, 8.0, 1.0, alpha)
        .total_lift(SPEED, DENSITY)
        .expect("lift n_chord=4");
    let rel = (four - one).abs() / one;
    assert!(
        rel < 0.03,
        "n_chord=1 gives {one}, n_chord=4 gives {four} (relative difference {rel})"
    );
    // The total must use every panel: summing the circulations of all panels
    // of the 4-chord wing reproduces the lift.
    let vlm = VortexLatticeMethod::new(16, 4, 8.0, 1.0, alpha);
    let gammas = vlm.solve_gammas(SPEED, DENSITY).expect("gammas");
    let dx = 8.0 / 16.0;
    let by_hand = DENSITY * SPEED * gammas.iter().sum::<f64>() * dx;
    assert!(((four - by_hand) / four).abs() < 1e-12);
}

#[test]
fn spanwise_circulation_is_symmetric_and_peaks_at_mid_span() {
    let vlm = VortexLatticeMethod::new(16, 1, 8.0, 1.0, 2.0_f64.to_radians());
    let gammas = vlm.solve_gammas(SPEED, DENSITY).expect("solve");
    let n = gammas.len();
    for i in 0..n / 2 {
        let mirror = n - 1 - i;
        assert!(
            (gammas[i] - gammas[mirror]).abs() < 1e-9 * gammas[i].abs(),
            "Γ[{i}] = {} differs from its mirror Γ[{mirror}] = {}",
            gammas[i],
            gammas[mirror]
        );
    }
    for i in 0..n / 2 - 1 {
        assert!(
            gammas[i] < gammas[i + 1],
            "circulation must grow monotonically towards mid-span: {gammas:?}"
        );
    }
    let tip = gammas[0];
    let centre = gammas[n / 2];
    assert!(
        tip < 0.7 * centre,
        "the tip circulation must drop well below the root value: tip {tip}, centre {centre}"
    );
}

#[test]
fn rectangular_wing_lift_slope_matches_glauert_lifting_line() {
    // AR = 8 rectangular wing.  The converged VLM slope is about 4.65-4.75
    // per radian; Glauert lifting-line theory gives 4.84.  (The elliptic-wing
    // value 2π·AR/(AR+2) = 5.03 does not apply to a rectangular planform.)
    let aspect_ratio = 8.0;
    let reference = glauert_rectangular_lift_slope(aspect_ratio);
    assert!(
        (reference - 4.8377).abs() < 5e-3,
        "Glauert reference solution changed: {reference}"
    );
    let slope = lift_slope(32, 4, aspect_ratio, 1.0);
    let rel = (slope - reference).abs() / reference;
    assert!(
        rel < 0.05,
        "VLM slope {slope} vs lifting-line slope {reference} (relative difference {rel})"
    );
    let elliptic_ideal = 2.0 * PI * aspect_ratio / (aspect_ratio + 2.0);
    assert!(
        slope < elliptic_ideal && slope > 0.88 * elliptic_ideal,
        "a rectangular wing is less efficient than the elliptic ideal {elliptic_ideal}, got {slope}"
    );
}

#[test]
fn high_aspect_ratio_wing_approaches_the_elliptic_lift_slope() {
    // For large AR the planform matters little and the lift slope tends to
    // 2π·AR/(AR+2).
    let aspect_ratio = 50.0;
    let ideal = 2.0 * PI * aspect_ratio / (aspect_ratio + 2.0);
    let slope = lift_slope(100, 2, aspect_ratio, 1.0);
    let rel = (slope - ideal).abs() / ideal;
    assert!(
        rel < 0.05,
        "AR = 50 lift slope {slope} vs 2π·AR/(AR+2) = {ideal} (relative difference {rel})"
    );
}

#[test]
fn lift_slopes_match_an_independent_reference_implementation() {
    // Reference values from a separately coded vortex lattice (numpy, finite
    // straight-segment Biot-Savart, LAPACK solve) with the same panelling:
    // bound vortex at 1/4 chord, collocation point at 3/4 chord, trailing legs
    // 100 spans long, freestream tilted by the angle of attack.
    let cases = [
        (8, 1, 8.0, 1.0, 4.870_554_223_401_475),
        (16, 2, 8.0, 1.0, 4.746_406_130_210_446),
        (32, 4, 8.0, 1.0, 4.671_252_018_520_893),
        (6, 3, 4.0, 1.0, 3.992_617_074_997_988),
        (10, 2, 3.0, 1.5, 2.660_484_940_999_510_6),
    ];
    for (n_span, n_chord, span, chord, expected) in cases {
        let slope = lift_slope(n_span, n_chord, span, chord);
        assert!(
            ((slope - expected) / expected).abs() < 1e-8,
            "{n_span} x {n_chord} panels, span {span}, chord {chord}: slope {slope} vs {expected}"
        );
    }
}

#[test]
fn lift_slope_decreases_with_decreasing_aspect_ratio() {
    let slope_8 = lift_slope(24, 2, 8.0, 1.0);
    let slope_4 = lift_slope(24, 2, 4.0, 1.0);
    let slope_2 = lift_slope(24, 2, 2.0, 1.0);
    assert!(
        slope_2 < slope_4 && slope_4 < slope_8,
        "finite-wing lift slopes must grow with aspect ratio: {slope_2}, {slope_4}, {slope_8}"
    );
}

#[test]
fn invalid_inputs_are_reported_as_errors() {
    let vlm = VortexLatticeMethod::new(4, 2, 4.0, 1.0, 0.05);
    assert!(vlm.solve_gammas(f64::NAN, DENSITY).is_err());
    assert!(vlm.solve_gammas(f64::INFINITY, DENSITY).is_err());
    assert!(vlm.solve_gammas(-1.0, DENSITY).is_err());
    assert!(vlm.total_lift(SPEED, f64::NAN).is_err());
    let flat = VortexLatticeMethod::new(4, 2, 0.0, 1.0, 0.05);
    assert!(flat.solve_gammas(SPEED, DENSITY).is_err(), "zero span");
    let no_chord = VortexLatticeMethod::new(4, 2, 4.0, 0.0, 0.05);
    assert!(no_chord.solve_gammas(SPEED, DENSITY).is_err(), "zero chord");
    let bad_alpha = VortexLatticeMethod::new(4, 2, 4.0, 1.0, f64::NAN);
    assert!(bad_alpha.total_lift(SPEED, DENSITY).is_err(), "NaN alpha");
    // Zero freestream speed is a valid state with zero circulation.
    let still = vlm.solve_gammas(0.0, DENSITY).expect("rest is valid");
    assert!(still.iter().all(|&g| g == 0.0));
}

#[test]
fn inconsistent_panel_arrays_are_reported() {
    let mut vlm = VortexLatticeMethod::new(4, 2, 4.0, 1.0, 0.05);
    vlm.panel_centers.pop();
    let err = vlm
        .solve_gammas(SPEED, DENSITY)
        .expect_err("a missing collocation point must be rejected");
    assert!(err.to_string().contains("do not match"), "got: {err}");
}

#[test]
fn singular_influence_matrix_is_reported_instead_of_returning_zeros() {
    // Two identical panels make two identical rows in the AIC matrix.
    let mut vlm = VortexLatticeMethod::new(2, 1, 2.0, 1.0, 0.05);
    vlm.panel_centers[1] = vlm.panel_centers[0];
    vlm.bound_vortex_pts[1] = vlm.bound_vortex_pts[0];
    let err = vlm
        .solve_gammas(SPEED, DENSITY)
        .expect_err("duplicate panels give a singular system");
    assert!(err.to_string().contains("singular"), "got: {err}");
    assert!(vlm.total_lift(SPEED, DENSITY).is_err());
}
