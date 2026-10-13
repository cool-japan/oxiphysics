// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Combined-slip friction-circle regression tests for the tire force models.
//!
//! Every model must keep the resultant of its lateral and longitudinal force
//! below `mu * Fz`, while pure-slip forces stay exactly those of the underlying
//! Magic Formula / linear law.

use oxiphysics_vehicle::tire::{TransientBrushModel, TransientBrushState};
use oxiphysics_vehicle::{LinearTire, PacejkaCoeffs, PacejkaTire, TireModel};
use proptest::prelude::*;

/// Relative slack allowed on `hypot(Fy, Fx) <= mu * Fz`.
const CIRCLE_TOL: f64 = 1e-9;

fn resultant(forces: (f64, f64)) -> f64 {
    forces.0.hypot(forces.1)
}

/// Independent evaluation of the Magic Formula
/// `D * Fz * mu * sin(C * atan(B x - E (B x - atan(B x))))`.
fn magic_formula(coeffs: &PacejkaCoeffs, x: f64, normal: f64, mu: f64) -> f64 {
    let bx = coeffs.b * x;
    let inner = bx - coeffs.e * (bx - bx.atan());
    coeffs.d * normal * mu * (coeffs.c * inner.atan()).sin()
}

#[test]
fn pacejka_old_overshoot_case_is_limited() {
    // Regression: tan(alpha) = kappa = 0.3 with the default weight 0.5
    // produced a resultant of about 1.10 * mu * Fz.
    let tire = PacejkaTire::new();
    let (normal, mu) = (5000.0, 1.0);
    let alpha = 0.3_f64.atan();
    let forces = tire.compute_forces(alpha, 0.3, normal, mu);
    assert!(
        resultant(forces) <= mu * normal * (1.0 + CIRCLE_TOL),
        "resultant {} exceeds mu*Fz = {}",
        resultant(forces),
        mu * normal
    );
    // Both components stay active (the tire is not simply switched off).
    assert!(forces.0 > 0.0 && forces.1 > 0.0, "{forces:?}");
}

#[test]
fn pacejka_resultant_never_exceeds_friction_circle_on_grid() {
    let tire = PacejkaTire::new();
    let third_pi = std::f64::consts::FRAC_PI_3;
    for &normal in &[500.0, 4000.0, 15000.0] {
        for &mu in &[0.3, 1.0, 1.8] {
            for ia in 0..=120 {
                let alpha = -third_pi + 2.0 * third_pi * f64::from(ia) / 120.0;
                for ik in 0..=80 {
                    let kappa = -1.0 + 2.0 * f64::from(ik) / 80.0;
                    let f = tire.compute_forces(alpha, kappa, normal, mu);
                    assert!(
                        resultant(f) <= mu * normal * (1.0 + CIRCLE_TOL),
                        "alpha={alpha}, kappa={kappa}: |F|={} > mu*Fz={}",
                        resultant(f),
                        mu * normal
                    );
                }
            }
        }
    }
}

#[test]
fn pacejka_pure_slip_forces_are_unchanged() {
    let tire = PacejkaTire::new();
    let (normal, mu) = (5000.0, 1.1);
    for i in -40..=40 {
        let alpha = f64::from(i) * 0.02;
        let (fy, fx) = tire.compute_forces(alpha, 0.0, normal, mu);
        let expected = magic_formula(&tire.lateral, alpha.tan(), normal, mu);
        assert!(
            (fy - expected).abs() <= 1e-9 * expected.abs().max(1.0),
            "pure lateral slip alpha={alpha}: {fy} vs {expected}"
        );
        assert_eq!(fx, 0.0, "no longitudinal force at zero slip ratio");
    }
    for i in -40..=40 {
        let kappa = f64::from(i) * 0.025;
        let (fy, fx) = tire.compute_forces(0.0, kappa, normal, mu);
        let expected = magic_formula(&tire.longitudinal, kappa, normal, mu);
        assert!(
            (fx - expected).abs() <= 1e-9 * expected.abs().max(1.0),
            "pure longitudinal slip kappa={kappa}: {fx} vs {expected}"
        );
        assert_eq!(fy, 0.0, "no lateral force at zero slip angle");
    }
}

#[test]
fn pacejka_combined_slip_reduces_each_axis() {
    let tire = PacejkaTire::new();
    let (normal, mu) = (5000.0, 1.0);
    let (fy_pure, _) = tire.compute_forces(0.1, 0.0, normal, mu);
    let (_, fx_pure) = tire.compute_forces(0.0, 0.1, normal, mu);
    let (fy_comb, fx_comb) = tire.compute_forces(0.1, 0.1, normal, mu);
    assert!(
        fy_comb < fy_pure,
        "slip ratio must reduce lateral force: {fy_comb} vs {fy_pure}"
    );
    assert!(
        fx_comb < fx_pure,
        "slip angle must reduce longitudinal force: {fx_comb} vs {fx_pure}"
    );
    assert!(fy_comb > 0.0 && fx_comb > 0.0);
}

#[test]
fn pacejka_full_coupling_matches_the_magic_formula_weighting_functions() {
    // With the coupling weight at 1 and the friction ellipse not binding, the
    // forces are the pure-slip forces times G_yk = cos(atan(10 k)) and
    // G_xa = cos(atan(12 tan a)).
    let mut tire = PacejkaTire::new();
    tire.combined_slip_weight = 1.0;
    let (normal, mu) = (5000.0_f64, 1.0_f64);
    let (alpha, kappa) = (0.04_f64, 0.03_f64);
    let (fy, fx) = tire.compute_forces(alpha, kappa, normal, mu);
    let fy0 = magic_formula(&tire.lateral, alpha.tan(), normal, mu);
    let fx0 = magic_formula(&tire.longitudinal, kappa, normal, mu);
    let g_y_kappa = (10.0 * kappa).atan().cos();
    let g_x_alpha = (12.0 * alpha.tan()).atan().cos();
    assert!(
        (fy - fy0 * g_y_kappa).abs() < 1e-9 * normal,
        "Fy = {fy}, expected {}",
        fy0 * g_y_kappa
    );
    assert!(
        (fx - fx0 * g_x_alpha).abs() < 1e-9 * normal,
        "Fx = {fx}, expected {}",
        fx0 * g_x_alpha
    );
    // The ellipse indeed does not bind here, so the weighting is what we see.
    let utilisation =
        (fy / (tire.lateral.d * mu * normal)).hypot(fx / (tire.longitudinal.d * mu * normal));
    assert!(utilisation < 0.95, "utilisation {utilisation}");
    // Half coupling interpolates linearly between the pure and the fully weighted force.
    tire.combined_slip_weight = 0.5;
    let (fy_half, fx_half) = tire.compute_forces(alpha, kappa, normal, mu);
    assert!((fy_half - fy0 * (0.5 + 0.5 * g_y_kappa)).abs() < 1e-9 * normal);
    assert!((fx_half - fx0 * (0.5 + 0.5 * g_x_alpha)).abs() < 1e-9 * normal);
}

#[test]
fn pacejka_forces_are_odd_in_the_slips() {
    let tire = PacejkaTire::new();
    let (fy, fx) = tire.compute_forces(0.12, 0.2, 4000.0, 1.0);
    let (fy_m, fx_m) = tire.compute_forces(-0.12, -0.2, 4000.0, 1.0);
    assert!((fy + fy_m).abs() < 1e-9 && (fx + fx_m).abs() < 1e-9);
    let (fy_a, fx_a) = tire.compute_forces(-0.12, 0.2, 4000.0, 1.0);
    assert!((fy_a + fy).abs() < 1e-9 && (fx_a - fx).abs() < 1e-9);
}

#[test]
fn pacejka_ellipse_uses_peak_factors_of_each_axis() {
    // A tire whose lateral peak factor is 1.3 may exceed mu*Fz in pure lateral
    // slip; the limit must not clip that, but must bound mixed slip.
    let tire = PacejkaTire::with_coefficients(
        PacejkaCoeffs {
            b: 10.0,
            c: 1.9,
            d: 1.3,
            e: 0.97,
        },
        PacejkaCoeffs {
            b: 12.0,
            c: 2.3,
            d: 0.9,
            e: 0.97,
        },
        1.0,
    );
    let (normal, mu) = (5000.0, 1.0);
    let (fy_pure, _) = tire.compute_forces(0.2, 0.0, normal, mu);
    let expected = magic_formula(&tire.lateral, 0.2_f64.tan(), normal, mu);
    assert!((fy_pure - expected).abs() <= 1e-9 * expected);
    assert!(
        fy_pure > mu * normal,
        "pure lateral peak must reach D*mu*Fz"
    );
    for ia in 0..=30 {
        let alpha = -0.5 + f64::from(ia) / 30.0;
        for ik in 0..=30 {
            let kappa = -0.5 + f64::from(ik) / 30.0;
            let (fy, fx) = tire.compute_forces(alpha, kappa, normal, mu);
            let utilisation = (fy / (1.3 * mu * normal)).hypot(fx / (0.9 * mu * normal));
            assert!(
                utilisation <= 1.0 + CIRCLE_TOL,
                "alpha={alpha}, kappa={kappa}: ellipse utilisation {utilisation}"
            );
        }
    }
}

#[test]
fn pacejka_stiffness_factors_control_the_coupling() {
    // A larger B_yk makes the lateral force collapse faster with slip ratio.
    let soft = PacejkaTire::new().with_combined_slip_stiffness(12.0, 5.0);
    let stiff = PacejkaTire::new().with_combined_slip_stiffness(12.0, 20.0);
    let (normal, mu) = (5000.0, 1.0);
    let (fy_soft, _) = soft.compute_forces(0.05, 0.05, normal, mu);
    let (fy_stiff, _) = stiff.compute_forces(0.05, 0.05, normal, mu);
    assert!(
        fy_stiff < fy_soft,
        "larger B_yk must reduce Fy more: {fy_stiff} vs {fy_soft}"
    );
    assert_eq!(PacejkaTire::new().b_x_alpha, 12.0);
    assert_eq!(PacejkaTire::new().b_y_kappa, 10.0);
}

#[test]
fn pacejka_without_grip_returns_zero_forces() {
    let tire = PacejkaTire::new();
    assert_eq!(tire.compute_forces(0.1, 0.1, -100.0, 1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.1, 0.1, 100.0, 0.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.1, 0.1, 100.0, -1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.1, 0.1, f64::NAN, 1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.1, 0.1, 100.0, f64::NAN), (0.0, 0.0));
}

#[test]
fn pacejka_weight_zero_still_respects_friction_circle() {
    let mut tire = PacejkaTire::new();
    tire.combined_slip_weight = 0.0;
    let (normal, mu) = (5000.0, 1.0);
    for ia in 0..=40 {
        let alpha = -0.6 + 1.2 * f64::from(ia) / 40.0;
        for ik in 0..=40 {
            let kappa = -0.6 + 1.2 * f64::from(ik) / 40.0;
            let f = tire.compute_forces(alpha, kappa, normal, mu);
            assert!(resultant(f) <= mu * normal * (1.0 + CIRCLE_TOL));
        }
    }
}

#[test]
fn linear_tire_resultant_is_limited_to_friction_circle() {
    let tire = LinearTire::new(40000.0, 60000.0);
    let (normal, mu) = (3000.0, 0.9);
    // Both axes would exceed mu*Fz on their own when clamped independently.
    let (fy, fx) = tire.compute_forces(0.06, 0.04, normal, mu);
    assert!(resultant((fy, fx)) <= mu * normal * (1.0 + CIRCLE_TOL));
    // Radial scaling preserves the force direction.
    assert!(
        (fy / fx - (40000.0 * 0.06) / (60000.0 * 0.04)).abs() < 1e-12,
        "direction changed: fy={fy}, fx={fx}"
    );
    assert!((resultant((fy, fx)) - mu * normal).abs() < 1e-6);
}

#[test]
fn linear_tire_unsaturated_forces_are_unchanged_and_guard_grip() {
    let tire = LinearTire::new(40000.0, 60000.0);
    let (fy, fx) = tire.compute_forces(0.01, -0.02, 5000.0, 1.0);
    assert!((fy - 400.0).abs() < 1e-9 && (fx + 1200.0).abs() < 1e-9);
    // Negative or NaN load/friction must not panic and produce no force.
    assert_eq!(tire.compute_forces(0.5, 0.5, -10.0, 1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.5, 0.5, 10.0, -1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.5, 0.5, f64::NAN, 1.0), (0.0, 0.0));
    assert_eq!(tire.compute_forces(0.5, 0.5, 10.0, f64::NAN), (0.0, 0.0));
}

#[test]
fn transient_brush_forces_respect_friction_circle() {
    let brush = TransientBrushModel::new(50000.0, 80000.0, 0.3, 0.15);
    let state = TransientBrushState {
        slip_lat: 0.2,
        slip_lon: 0.15,
    };
    let (normal, mu) = (4000.0, 1.0);
    let (fy, fx) = state.forces(&brush, normal, mu);
    assert!(resultant((fy, fx)) <= mu * normal * (1.0 + CIRCLE_TOL));
    assert!((resultant((fy, fx)) - mu * normal).abs() < 1e-6);
    assert!(fy > 0.0 && fx > 0.0);
    // Small slips stay linear.
    let linear = TransientBrushState {
        slip_lat: 0.01,
        slip_lon: 0.005,
    };
    let (fy, fx) = linear.forces(&brush, normal, mu);
    assert!((fy - 500.0).abs() < 1e-9 && (fx - 400.0).abs() < 1e-9);
}

#[test]
fn transient_brush_forces_without_grip_do_not_panic() {
    let brush = TransientBrushModel::new(50000.0, 80000.0, 0.3, 0.15);
    let state = TransientBrushState {
        slip_lat: 0.2,
        slip_lon: 0.1,
    };
    assert_eq!(state.forces(&brush, 0.0, 1.0), (0.0, 0.0));
    assert_eq!(state.forces(&brush, -500.0, 1.0), (0.0, 0.0));
    assert_eq!(state.forces(&brush, 500.0, -1.0), (0.0, 0.0));
    assert_eq!(state.forces(&brush, f64::NAN, 1.0), (0.0, 0.0));
}

proptest! {
    #[test]
    fn pacejka_resultant_within_friction_circle(
        alpha in -std::f64::consts::FRAC_PI_3..std::f64::consts::FRAC_PI_3,
        kappa in -1.0f64..1.0,
        normal in 100.0f64..20000.0,
        mu in 0.2f64..2.0,
    ) {
        let tire = PacejkaTire::new();
        let f = tire.compute_forces(alpha, kappa, normal, mu);
        prop_assert!(resultant(f) <= mu * normal * (1.0 + CIRCLE_TOL));
    }

    #[test]
    fn linear_resultant_within_friction_circle(
        alpha in -std::f64::consts::FRAC_PI_3..std::f64::consts::FRAC_PI_3,
        kappa in -1.0f64..1.0,
        normal in 100.0f64..20000.0,
        mu in 0.2f64..2.0,
    ) {
        let tire = LinearTire::default();
        let f = tire.compute_forces(alpha, kappa, normal, mu);
        prop_assert!(resultant(f) <= mu * normal * (1.0 + CIRCLE_TOL));
    }

    #[test]
    fn transient_brush_resultant_within_friction_circle(
        slip_lat in -1.0f64..1.0,
        slip_lon in -1.0f64..1.0,
        normal in 100.0f64..20000.0,
        mu in 0.2f64..2.0,
    ) {
        let brush = TransientBrushModel::new(50000.0, 80000.0, 0.3, 0.15);
        let state = TransientBrushState { slip_lat, slip_lon };
        let f = state.forces(&brush, normal, mu);
        prop_assert!(resultant(f) <= mu * normal * (1.0 + CIRCLE_TOL));
    }

    #[test]
    fn pacejka_pure_lateral_slip_matches_magic_formula(
        alpha in -std::f64::consts::FRAC_PI_3..std::f64::consts::FRAC_PI_3,
        normal in 100.0f64..20000.0,
        mu in 0.2f64..2.0,
    ) {
        let tire = PacejkaTire::new();
        let (fy, fx) = tire.compute_forces(alpha, 0.0, normal, mu);
        let expected = magic_formula(&tire.lateral, alpha.tan(), normal, mu);
        prop_assert!((fy - expected).abs() <= 1e-9 * normal * mu);
        prop_assert_eq!(fx, 0.0);
    }

    #[test]
    fn pacejka_pure_longitudinal_slip_matches_magic_formula(
        kappa in -1.0f64..1.0,
        normal in 100.0f64..20000.0,
        mu in 0.2f64..2.0,
    ) {
        let tire = PacejkaTire::new();
        let (fy, fx) = tire.compute_forces(0.0, kappa, normal, mu);
        let expected = magic_formula(&tire.longitudinal, kappa, normal, mu);
        prop_assert!((fx - expected).abs() <= 1e-9 * normal * mu);
        prop_assert_eq!(fy, 0.0);
    }
}
