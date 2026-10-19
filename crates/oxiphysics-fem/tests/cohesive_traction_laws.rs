// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the PPR and Xu-Needleman cohesive laws.
//!
//! * `PprParams` used an ad hoc polynomial whose work of separation was
//!   `1.6·10⁵·Γ_n` for `α = 4, λ = 0.1` (and `2.67·Γ_n` for `α = 2, λ = 0.5`)
//!   and whose peak was not at `λ·δ_n`.  The true PPR mode-I law integrates to
//!   exactly `Γ_n`.
//! * `XuNeedlemanParams::normal_traction` carried a spurious `(1 + Δn/δn)`
//!   factor, so the traction did not vanish at zero opening (and was not the
//!   derivative of the documented potential).

use oxiphysics_fem::Error;
use oxiphysics_fem::cohesive::{PprParams, XuNeedlemanParams};

fn ppr(gamma_n: f64, alpha: f64, lambda: f64, delta_n: f64) -> PprParams {
    PprParams::new(
        gamma_n,
        0.5 * gamma_n,
        alpha,
        alpha,
        delta_n,
        delta_n,
        lambda,
        lambda,
    )
}

/// Rectangle-rule independent check of the work of separation.
fn simpson_energy(params: &PprParams, intervals: usize) -> f64 {
    let h = params.delta_n / intervals as f64;
    let mut sum = params.normal_traction(0.0) + params.normal_traction(params.delta_n);
    for i in 1..intervals {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * params.normal_traction(i as f64 * h);
    }
    sum * h / 3.0
}

#[test]
fn work_of_separation_equals_the_fracture_energy() {
    // (α, λ_n): m = α(α−1)λ²/(1−αλ²) takes values 0.125, 1, 1, 0.08, 1.85, 6, 49.
    let cases = [
        (4.0, 0.1),
        (2.0, 0.5),
        (5.0, 0.2),
        (1.5, 0.3),
        (3.0, 0.4),
        (10.0, 0.2),
        (2.0, 0.7),
    ];
    for (alpha, lambda) in cases {
        let params = ppr(500.0, alpha, lambda, 1.0e-4);
        let energy = params.mode_i_energy_numerical();
        // The singular (1−s)^(α−1) end for α < 2 limits the quadrature accuracy.
        let tolerance = if alpha < 2.0 { 1.0e-4 } else { 1.0e-9 };
        assert!(
            (energy - 500.0).abs() <= tolerance * 500.0,
            "α = {alpha}, λ = {lambda}: ∫T dΔ = {energy}, Γ_n = 500"
        );
    }
}

#[test]
fn work_of_separation_is_confirmed_by_an_independent_quadrature() {
    let params = ppr(80.0, 4.0, 0.1, 2.0e-4);
    let energy = simpson_energy(&params, 20_000);
    assert!(
        (energy - 80.0).abs() < 1.0e-6 * 80.0,
        "Simpson energy {energy}"
    );
}

#[test]
fn peak_traction_sits_at_lambda_times_delta_n() {
    for (alpha, lambda) in [(4.0, 0.1), (2.0, 0.5), (5.0, 0.2), (3.0, 0.35)] {
        let params = ppr(300.0, alpha, lambda, 2.0e-4);
        let peak_at = params.delta_n_peak();
        assert!((peak_at - lambda * 2.0e-4).abs() < 1e-18);
        let sigma_max = params.sigma_max();
        let at_peak = params.normal_traction(peak_at);
        assert!(
            (at_peak - sigma_max).abs() <= 1e-12 * sigma_max,
            "α = {alpha}: T(δ_peak) = {at_peak}, σ_max = {sigma_max}"
        );
        // No sample over the opening range exceeds the peak, and the largest one
        // is the sample closest to the peak.
        let samples = 4000;
        let mut best = (0.0_f64, 0.0_f64);
        for i in 1..samples {
            let d = params.delta_n * i as f64 / samples as f64;
            let t = params.normal_traction(d);
            if t > best.1 {
                best = (d, t);
            }
        }
        assert!(best.1 <= sigma_max * (1.0 + 1e-12), "α = {alpha}");
        assert!(
            (best.0 - peak_at).abs() <= 2.0 * params.delta_n / samples as f64,
            "α = {alpha}: scanned maximum at {}, expected {peak_at}",
            best.0
        );
    }
}

/// `α = 2`, `λ = 0.5` gives `m = 1` and `T = 6·Γ_n·s(1 − s)/δ_n`.
#[test]
fn alpha_two_lambda_half_is_a_parabola() {
    let (gamma, delta) = (200.0, 1.0e-3);
    let params = ppr(gamma, 2.0, 0.5, delta);
    assert!((params.sigma_max() - 1.5 * gamma / delta).abs() < 1e-9 * gamma / delta);
    for s in [0.1, 0.25, 0.5, 0.8] {
        let expected = 6.0 * gamma / delta * s * (1.0 - s);
        let got = params.normal_traction(s * delta);
        assert!(
            (got - expected).abs() <= 1e-12 * expected,
            "s = {s}: {got} vs {expected}"
        );
    }
}

/// `α = 5`, `λ = 0.2` gives `m = 1` and `T = 30·Γ_n·s(1 − s)⁴/δ_n`.
#[test]
fn alpha_five_lambda_one_fifth_closed_form() {
    let (gamma, delta) = (50.0, 5.0e-4);
    let params = ppr(gamma, 5.0, 0.2, delta);
    assert!((params.normal_shape_exponent().expect("valid") - 1.0).abs() < 1e-14);
    for s in [0.05, 0.2, 0.5, 0.9] {
        let expected = 30.0 * gamma / delta * s * (1.0 - s).powi(4);
        let got = params.normal_traction(s * delta);
        assert!(
            (got - expected).abs() <= 1e-12 * expected,
            "s = {s}: {got} vs {expected}"
        );
    }
    assert!((params.sigma_max() - 2.4576 * gamma / delta).abs() < 1e-9 * gamma / delta);
}

#[test]
fn traction_vanishes_at_both_ends_and_outside_the_range() {
    let params = ppr(500.0, 4.0, 0.1, 1.0e-4);
    assert_eq!(params.normal_traction(0.0), 0.0);
    assert_eq!(params.normal_traction(-1.0e-5), 0.0);
    assert_eq!(params.normal_traction(1.0e-4), 0.0);
    assert_eq!(params.normal_traction(2.0e-4), 0.0);
    for i in 1..100 {
        let t = params.normal_traction(1.0e-4 * i as f64 / 100.0);
        assert!(t > 0.0 && t.is_finite(), "T({i}/100 δ) = {t}");
    }
}

#[test]
fn parameters_outside_the_ppr_domain_are_reported_not_invented() {
    // λ ≥ 1/√α makes 1 − αλ² non-positive: no finite positive exponent m.
    let beyond = ppr(500.0, 4.0, 0.5, 1.0e-4);
    assert!(matches!(beyond.validate(), Err(Error::InvalidInput(_))));
    assert!(beyond.normal_shape_exponent().is_none());
    assert!(beyond.sigma_max().is_nan());
    assert!(beyond.normal_traction(0.5e-4).is_nan());
    assert!(beyond.mode_i_energy_numerical().is_nan());

    for (gamma_n, alpha, delta_n, lambda) in [
        (-1.0, 4.0, 1e-4, 0.1),
        (500.0, 1.0, 1e-4, 0.1),
        (500.0, 0.5, 1e-4, 0.1),
        (500.0, 4.0, 0.0, 0.1),
        (500.0, 4.0, 1e-4, 0.0),
        (500.0, 4.0, 1e-4, -0.1),
        (f64::NAN, 4.0, 1e-4, 0.1),
        (500.0, f64::INFINITY, 1e-4, 0.1),
    ] {
        assert!(
            PprParams::try_new(gamma_n, 250.0, alpha, 4.0, delta_n, 1e-4, lambda, 0.1).is_err(),
            "({gamma_n}, {alpha}, {delta_n}, {lambda}) must be rejected"
        );
    }
    // The tangential parameters are validated as well.
    assert!(PprParams::try_new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.6).is_err());
    assert!(PprParams::try_new(500.0, 250.0, 4.0, 1.0, 1e-4, 1e-4, 0.1, 0.1).is_err());
    assert!(PprParams::try_new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.2).is_ok());
}

// ---------------------------------------------------------------------------
// Xu-Needleman
// ---------------------------------------------------------------------------

fn xu(phi_n: f64, phi_t: f64) -> XuNeedlemanParams {
    XuNeedlemanParams::new(phi_n, phi_t, 1.0e-4, 2.0e-4)
}

#[test]
fn xu_needleman_normal_traction_is_zero_at_zero_opening() {
    for params in [xu(500.0, 500.0), xu(500.0, 250.0), xu(500.0, 900.0)] {
        assert_eq!(params.normal_traction(0.0, 0.0), 0.0);
        // Also under pure tangential opening: Tn carries the factor Δn/δn.
        assert_eq!(params.normal_traction(0.0, 1.0e-4), 0.0);
        assert_eq!(params.tangential_traction(0.0, 0.0), 0.0);
        assert_eq!(params.potential(0.0, 0.0), 0.0);
    }
}

#[test]
fn xu_needleman_normal_peak_is_sigma_max_at_delta_n() {
    let params = xu(500.0, 250.0);
    let sigma_max = params.sigma_max();
    let at_characteristic_length = params.normal_traction(params.delta_n, 0.0);
    assert!(
        (at_characteristic_length - sigma_max).abs() <= 1e-12 * sigma_max,
        "T_n(δ_n) = {at_characteristic_length}, σ_max = {sigma_max}"
    );
    let samples = 4000;
    let mut best = (0.0_f64, 0.0_f64);
    for i in 0..=samples {
        let d = 6.0 * params.delta_n * i as f64 / samples as f64;
        let t = params.normal_traction(d, 0.0);
        if t > best.1 {
            best = (d, t);
        }
    }
    assert!(best.1 <= sigma_max * (1.0 + 1e-12));
    assert!((best.0 - params.delta_n).abs() <= 6.0 * params.delta_n / samples as f64);
}

#[test]
fn xu_needleman_tangential_peak_is_tau_max_at_delta_t_over_root_two() {
    let params = xu(500.0, 250.0);
    let at_peak = params.tangential_traction(0.0, params.delta_t / std::f64::consts::SQRT_2);
    assert!((at_peak - params.tau_max()).abs() <= 1e-12 * params.tau_max());
    // The traction is odd in Δt and decays for large tangential opening.
    assert!(params.tangential_traction(0.0, 5.0 * params.delta_t) < 1e-9 * params.tau_max());
}

/// `T_n = ∂Φ/∂Δn`, `T_t = ∂Φ/∂Δt` by central differences, including a
/// `φ_t ≠ φ_n` set (the documented potential only covered `φ_t = φ_n`).
#[test]
fn xu_needleman_tractions_are_the_derivatives_of_the_potential() {
    let points = [
        (0.3e-4, 0.0),
        (1.0e-4, 0.7e-4),
        (2.5e-4, 1.5e-4),
        (0.5e-4, 3.0e-4),
        (4.0e-4, 0.2e-4),
    ];
    for params in [xu(500.0, 500.0), xu(500.0, 250.0), xu(400.0, 900.0)] {
        for &(dn, dt) in &points {
            let hn = 1.0e-4 * params.delta_n;
            let ht = 1.0e-4 * params.delta_t;
            let numeric_n =
                (params.potential(dn + hn, dt) - params.potential(dn - hn, dt)) / (2.0 * hn);
            let numeric_t =
                (params.potential(dn, dt + ht) - params.potential(dn, dt - ht)) / (2.0 * ht);
            let analytic_n = params.normal_traction(dn, dt);
            let analytic_t = params.tangential_traction(dn, dt);
            let scale_n = params.sigma_max();
            let scale_t = params.tau_max();
            assert!(
                (numeric_n - analytic_n).abs() <= 1e-6 * scale_n,
                "φ=({}, {}), Δ=({dn:e}, {dt:e}): dΦ/dΔn = {numeric_n}, T_n = {analytic_n}",
                params.phi_n,
                params.phi_t
            );
            assert!(
                (numeric_t - analytic_t).abs() <= 1e-6 * scale_t,
                "φ=({}, {}), Δ=({dn:e}, {dt:e}): dΦ/dΔt = {numeric_t}, T_t = {analytic_t}",
                params.phi_n,
                params.phi_t
            );
        }
    }
}

#[test]
fn xu_needleman_work_of_separation_limits() {
    let params = xu(500.0, 250.0);
    // Complete normal separation (Δt = 0) costs φ_n; complete tangential
    // separation (Δn = 0) costs φ_t.
    let normal_only = params.potential(40.0 * params.delta_n, 0.0);
    assert!((normal_only - 500.0).abs() < 1e-9, "{normal_only}");
    let tangential_only = params.potential(0.0, 40.0 * params.delta_t);
    assert!((tangential_only - 250.0).abs() < 1e-9, "{tangential_only}");
    // For φ_t = φ_n the documented single-term potential is recovered.
    let equal = xu(500.0, 500.0);
    let (dn, dt) = (1.3e-4, 0.8e-4);
    let documented = 500.0
        * (1.0
            - (1.0 + dn / 1.0e-4)
                * (-dn / 1.0e-4_f64).exp()
                * (-(dt * dt) / (2.0e-4 * 2.0e-4_f64)).exp());
    assert!((equal.potential(dn, dt) - documented).abs() < 1e-9);
}
