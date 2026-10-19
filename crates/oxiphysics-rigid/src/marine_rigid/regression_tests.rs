// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for `marine_rigid`: the JONSWAP spectrum must honour the
//! requested significant wave height.

use super::*;
use crate::marine_dynamics;

/// Trapezoidal zeroth moment of the spectrum over `[0.2 w_p, 6 w_p]`.
fn zeroth_moment(omega_peak: f64, hs: f64, gamma: f64) -> f64 {
    let n = 10_000;
    let lo = 0.2 * omega_peak;
    let hi = 6.0 * omega_peak;
    let d_omega = (hi - lo) / n as f64;
    let mut m0 = 0.0;
    for i in 0..n {
        let a = lo + i as f64 * d_omega;
        let b = a + d_omega;
        m0 += 0.5
            * d_omega
            * (jonswap_spectrum(a, omega_peak, hs, gamma, 9.81)
                + jonswap_spectrum(b, omega_peak, hs, gamma, 9.81));
    }
    m0
}

#[test]
fn jonswap_spectrum_integrates_to_the_requested_significant_wave_height() {
    let omega_peak = 0.8;
    for &gamma in &[1.0, 3.3, 5.0] {
        for &hs in &[0.5, 2.0, 10.0] {
            let m0 = zeroth_moment(omega_peak, hs, gamma);
            let hs_from_spectrum = significant_wave_height(m0);
            assert!(
                (hs_from_spectrum - hs).abs() < 0.03 * hs,
                "gamma = {gamma}, Hs = {hs}: spectrum integrates to Hs = {hs_from_spectrum}"
            );
        }
    }
}

#[test]
fn jonswap_spectrum_scales_with_hs_squared() {
    let omega_peak = 0.7;
    for &omega in &[0.5, omega_peak, 1.1, 2.0] {
        let s1 = jonswap_spectrum(omega, omega_peak, 1.0, 3.3, 9.81);
        let s2 = jonswap_spectrum(omega, omega_peak, 2.0, 3.3, 9.81);
        let s3 = jonswap_spectrum(omega, omega_peak, 3.0, 3.3, 9.81);
        assert!(s1 > 0.0);
        assert!(
            (s2 / s1 - 4.0).abs() < 1e-12,
            "omega = {omega}: {}",
            s2 / s1
        );
        assert!(
            (s3 / s1 - 9.0).abs() < 1e-12,
            "omega = {omega}: {}",
            s3 / s1
        );
    }
}

#[test]
fn jonswap_spectrum_with_unit_gamma_is_pierson_moskowitz() {
    let omega_peak = 0.9;
    let tp = 2.0 * PI / omega_peak;
    for &omega in &[0.4, 0.9, 1.5, 3.0] {
        let jonswap = jonswap_spectrum(omega, omega_peak, 2.5, 1.0, 9.81);
        let pm = marine_dynamics::pierson_moskowitz(omega, 2.5, tp);
        assert!(
            (jonswap - pm).abs() <= 1e-12 * pm.max(1e-300),
            "omega = {omega}: {jonswap} vs {pm}"
        );
    }
}

#[test]
fn jonswap_peak_value_matches_the_closed_form() {
    // At omega = omega_p: S = A_gamma * (5/16) Hs^2 e^(-5/4) / omega_p * gamma.
    let (omega_peak, hs, gamma) = (0.6, 4.0, 3.3_f64);
    let a_gamma = 1.0 - 0.287 * gamma.ln();
    let expected = a_gamma * 5.0 / 16.0 * hs * hs * (-1.25_f64).exp() / omega_peak * gamma;
    let s = jonswap_spectrum(omega_peak, omega_peak, hs, gamma, 9.81);
    assert!((s - expected).abs() < 1e-12 * expected, "{s} vs {expected}");
    // Peak enhancement raises the peak above the unenhanced spectrum.
    let flat = jonswap_spectrum(omega_peak, omega_peak, hs, 1.0, 9.81);
    assert!(s > flat);
}

#[test]
fn jonswap_spectrum_ignores_gravity_and_handles_degenerate_inputs() {
    let reference = jonswap_spectrum(1.0, 0.8, 2.0, 3.3, 9.81);
    assert!(reference > 0.0);
    assert_eq!(reference, jonswap_spectrum(1.0, 0.8, 2.0, 3.3, 1.62));
    assert_eq!(jonswap_spectrum(0.0, 0.8, 2.0, 3.3, 9.81), 0.0);
    assert_eq!(jonswap_spectrum(-1.0, 0.8, 2.0, 3.3, 9.81), 0.0);
    assert_eq!(jonswap_spectrum(1.0, 0.0, 2.0, 3.3, 9.81), 0.0);
    assert_eq!(jonswap_spectrum(1.0, -0.8, 2.0, 3.3, 9.81), 0.0);
    assert_eq!(jonswap_spectrum(1.0, 0.8, 0.0, 3.3, 9.81), 0.0);
}

// ── mooring line ────────────────────────────────────────────────────────

#[test]
fn mooring_line_matches_the_inextensible_catenary_closed_form() {
    let line = MooringLine::new(300.0, 1000.0, 250.0, 100.0);
    let w = line.weight_per_length;
    let h = line.catenary_horizontal_tension();
    let a = h / w;
    assert!(h > 0.0 && h.is_finite());
    // Suspended length L_s = sqrt(z (z + 2a)); X = L - L_s + a asinh(L_s / a).
    let suspended = (100.0 * (100.0 + 2.0 * a)).sqrt();
    let span = 300.0 - suspended + a * (suspended / a).asinh();
    assert!((span - 250.0).abs() < 1e-8, "span closure {span}");
    assert!(suspended < 300.0, "part of the line lies on the seabed");
    // z = a (cosh(x_s / a) - 1) on the suspended part.
    let x_s = a * (suspended / a).asinh();
    assert!((a * ((x_s / a).cosh() - 1.0) - 100.0).abs() < 1e-6);
    // The fairlead carries only the suspended weight: T_V = w L_s < w L.
    let fairlead = line.fairlead_tension();
    let expected = h.hypot(w * suspended);
    assert!((fairlead - expected).abs() < 1e-6 * expected);
    assert!(fairlead < h.hypot(w * 300.0));
    assert!((fairlead - w * (a + 100.0)).abs() < 1e-6 * fairlead);
}

#[test]
fn mooring_tension_grows_as_the_fairlead_moves_away() {
    let mut previous = 0.0;
    for span in [210.0, 230.0, 250.0, 265.0, 270.0] {
        let line = MooringLine::new(300.0, 1000.0, span, 100.0);
        let h = line.catenary_horizontal_tension();
        assert!(h > previous, "span {span}: {h} <= {previous}");
        assert!(line.fairlead_tension() >= h);
        previous = h;
    }
}

#[test]
fn mooring_profile_runs_from_the_anchor_to_the_fairlead() {
    for &(length, w, span, height) in &[
        (300.0, 1000.0, 250.0, 100.0),
        (300.0, 1000.0, 280.0, 100.0),
        (300.0, 1000.0, 150.0, 100.0),
    ] {
        let line = MooringLine::new(length, w, span, height);
        let profile = line.catenary_profile(41);
        assert_eq!(profile.len(), 41);
        let first = profile[0];
        let last = profile[40];
        assert!(first[0].abs() < 1e-9 && first[1].abs() < 1e-9, "{first:?}");
        assert!(
            (last[0] - span).abs() < 1e-6 * span && (last[1] - height).abs() < 1e-6 * height,
            "span {span}: last point {last:?}"
        );
        for pair in profile.windows(2) {
            assert!(pair[1][0] >= pair[0][0] - 1e-9, "x must not decrease");
            assert!(pair[1][1] >= pair[0][1] - 1e-9, "z must not decrease");
            assert!(pair[1][1] >= -1e-9);
        }
    }
}

#[test]
fn mooring_line_degenerate_inputs_do_not_panic_or_produce_nan() {
    let slack = MooringLine::new(300.0, 1000.0, 150.0, 100.0);
    assert_eq!(slack.catenary_horizontal_tension(), 0.0);
    // A slack line still hangs from the fairlead by its own weight.
    assert!((slack.fairlead_tension() - 1000.0 * 100.0).abs() < 1e-6 * 1.0e5);

    let broken = MooringLine::new(0.0, 1000.0, 100.0, 10.0);
    assert_eq!(broken.catenary_horizontal_tension(), 0.0);
    assert_eq!(broken.fairlead_tension(), 0.0);
    let straight = broken.catenary_profile(5);
    assert_eq!(straight.len(), 5);
    assert_eq!(straight[0], [0.0, 0.0]);
    assert_eq!(straight[4], [100.0, 10.0]);

    // Line shorter than the anchor-fairlead distance: infinite tension.
    let too_short = MooringLine::new(100.0, 1000.0, 100.0, 50.0);
    assert!(too_short.catenary_horizontal_tension().is_infinite());
    assert!(too_short.fairlead_tension().is_infinite());
    let profile = too_short.catenary_profile(3);
    assert!(profile.iter().all(|p| p[0].is_finite() && p[1].is_finite()));

    // A fairlead below the anchor is treated as being level with it.
    let below = MooringLine::new(300.0, 1000.0, 150.0, -5.0);
    assert_eq!(below.catenary_horizontal_tension(), 0.0);
    assert_eq!(below.catenary_profile(0).len(), 2);
}
