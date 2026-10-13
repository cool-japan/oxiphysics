// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the 1-D Euler Roe flux.
//!
//! The previous implementation used wave strengths that ignored the energy
//! jump and `γ`, so a pure pressure jump produced no dissipation at all and
//! `Σ α_k K_k` did not even reproduce the mass jump.  These tests check the
//! decomposition algebraically (`Σ α_k K_k = ΔU` and Roe's conservation
//! property `Σ λ_k α_k K_k = ΔF`) and the flux against exactly solvable
//! situations (isolated contacts, stationary shocks, supersonic upwinding,
//! sonic rarefactions).
//!
//! The Sod interface is deliberately *not* compared with
//! [`hllc_flux_euler_1d`] at a tight tolerance: at `x/t = 0` the exact contact
//! moves at `u* = 0.927` while the Roe linearisation places it at `ũ = 0`, and
//! HLLC with Davis wave speeds is itself 9-27 % away from the exact Godunov
//! flux, so the two correct solvers legitimately differ by about 10 % there.

use super::*;

const GAMMA: f64 = 1.4;

/// Riemann problems `(name, left, right)` with primitive states `[ρ, u, p]`.
const PROBLEMS: [(&str, [f64; 3], [f64; 3]); 6] = [
    ("sod", [1.0, 0.0, 1.0], [0.125, 0.0, 0.1]),
    ("pressure_jump", [1.0, 0.0, 1.0], [1.0, 0.0, 0.5]),
    ("toro_123", [1.0, -2.0, 0.4], [1.0, 2.0, 0.4]),
    ("toro_3", [1.0, 0.0, 1000.0], [1.0, 0.0, 0.01]),
    (
        "toro_4",
        [5.99924, 19.5975, 460.894],
        [5.99242, -6.19633, 46.0950],
    ),
    ("toro_5", [1.0, -19.59745, 1000.0], [1.0, -19.59745, 0.01]),
];

fn conserved(state: [f64; 3]) -> [f64; 3] {
    let [rho, u, p] = state;
    [rho, rho * u, p / (GAMMA - 1.0) + 0.5 * rho * u * u]
}

fn physical_flux(state: [f64; 3]) -> [f64; 3] {
    let [rho, u, p] = state;
    let e = p / (GAMMA - 1.0) + 0.5 * rho * u * u;
    [rho * u, rho * u * u + p, (e + p) * u]
}

fn max_abs(v: [f64; 3]) -> f64 {
    v.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}

fn roe(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    roe_flux_euler_1d(
        left[0], left[1], left[2], right[0], right[1], right[2], GAMMA,
    )
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `Σ_k weight_k · α_k · K_k` for the decomposition `dec`.
fn wave_sum(dec: &RoeDecomposition, weight: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0_f64; 3];
    for (i, o) in out.iter_mut().enumerate() {
        *o = (0..3)
            .map(|k| weight[k] * dec.alpha[k] * dec.eigenvectors[k][i])
            .sum();
    }
    out
}

#[test]
fn wave_strengths_reproduce_the_conserved_jump() {
    for (name, left, right) in PROBLEMS {
        let dec = roe_decomposition(left, right, GAMMA)
            .unwrap_or_else(|| panic!("{name}: Roe average must be valid"));
        let du = sub(conserved(right), conserved(left));
        let rebuilt = wave_sum(&dec, [1.0; 3]);
        let tol = 1e-12 * (1.0 + max_abs(du));
        for i in 0..3 {
            assert!(
                (rebuilt[i] - du[i]).abs() <= tol,
                "{name}: component {i}: sum alpha_k K_k = {} but dU = {}",
                rebuilt[i],
                du[i]
            );
        }
        // Eigenvalues are ordered u-a < u < u+a with a > 0.
        assert!(dec.sound_speed > 0.0, "{name}");
        assert!(dec.lambda[0] < dec.lambda[1] && dec.lambda[1] < dec.lambda[2]);
    }
}

/// Roe's property U: the linearisation is conservative, i.e.
/// `Σ λ_k α_k K_k = F_R - F_L` exactly (this needs the energy jump and `γ` in
/// the wave strengths).
#[test]
fn roe_linearisation_is_conservative() {
    for (name, left, right) in PROBLEMS {
        let dec = roe_decomposition(left, right, GAMMA)
            .unwrap_or_else(|| panic!("{name}: Roe average must be valid"));
        let df = sub(physical_flux(right), physical_flux(left));
        let rebuilt = wave_sum(&dec, dec.lambda);
        let tol = 1e-11 * (1.0 + max_abs(df));
        for i in 0..3 {
            assert!(
                (rebuilt[i] - df[i]).abs() <= tol,
                "{name}: component {i}: sum lambda_k alpha_k K_k = {} but dF = {}",
                rebuilt[i],
                df[i]
            );
        }
    }
}

/// A pure pressure jump (equal density, zero velocity) has `Δρ = Δ(ρu) = 0`
/// and used to get zero dissipation.  The mass flux must be the closed form
/// `-Δp / (2 ã)` (all numerical dissipation, since the central average of the
/// physical mass fluxes is zero).
#[test]
fn pressure_only_jump_is_dissipated() {
    let left = [1.0, 0.0, 1.0];
    let right = [1.0, 0.0, 0.5];
    let flux = roe(left, right);

    // ρ = 1, u = 0 on both sides: H = p/(γ-1) + p and ã² = (γ-1) · mean(H).
    let enthalpy = |p: f64| p / (GAMMA - 1.0) + p;
    let a_roe = ((GAMMA - 1.0) * 0.5 * (enthalpy(1.0) + enthalpy(0.5))).sqrt();
    let delta_p = 0.5 - 1.0;
    let expected_mass = -delta_p / (2.0 * a_roe);

    assert!(
        flux[0] > 0.1,
        "mass must flow from high to low pressure: {}",
        flux[0]
    );
    assert!(
        (flux[0] - expected_mass).abs() < 1e-12,
        "mass flux {} vs closed form {expected_mass}",
        flux[0]
    );
    // The momentum flux of the symmetric acoustic pair is the mean pressure.
    assert!((flux[1] - 0.75).abs() < 1e-12, "momentum flux {}", flux[1]);
}

/// An isolated contact is resolved exactly: the interface sees the upwind
/// physical state.
#[test]
fn isolated_contact_is_upwinded_exactly() {
    for u in [1.0, 0.25, -0.25, -1.0] {
        let left = [1.0, u, 1.0];
        let right = [0.5, u, 1.0];
        let flux = roe(left, right);
        let expected = if u > 0.0 {
            physical_flux(left)
        } else {
            physical_flux(right)
        };
        let tol = 1e-12 * (1.0 + max_abs(expected));
        for i in 0..3 {
            assert!(
                (flux[i] - expected[i]).abs() <= tol,
                "u={u} component {i}: {} vs {}",
                flux[i],
                expected[i]
            );
        }
    }
}

/// Normal-shock relations for a stationary shock in front of which the flow
/// has Mach number `mach` (`ρ₁ = p₁ = 1`).
fn stationary_shock(mach: f64) -> ([f64; 3], [f64; 3]) {
    let (rho1, p1) = (1.0, 1.0);
    let a1 = (GAMMA * p1 / rho1).sqrt();
    let u1 = mach * a1;
    let m2 = mach * mach;
    let p2 = p1 * (1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (m2 - 1.0));
    let rho2 = rho1 * ((GAMMA + 1.0) * m2) / ((GAMMA - 1.0) * m2 + 2.0);
    let u2 = u1 * rho1 / rho2;
    ([rho1, u1, p1], [rho2, u2, p2])
}

/// Roe is exact for an isolated stationary shock: all three fluxes are
/// continuous and the numerical flux equals them.  Harten's entropy fix adds
/// a bounded amount of dissipation to the `λ = 0` wave (documented smearing).
#[test]
fn stationary_shock_is_exact_without_entropy_fix() {
    let (left, right) = stationary_shock(2.0);
    let fl = physical_flux(left);
    let fr = physical_flux(right);
    for i in 0..3 {
        assert!(
            (fl[i] - fr[i]).abs() < 1e-12 * (1.0 + fl[i].abs()),
            "Rankine-Hugoniot violated in component {i}"
        );
    }

    let exact = roe_flux_euler_1d_with_entropy_fix(left, right, GAMMA, 0.0);
    for i in 0..3 {
        assert!(
            (exact[i] - fl[i]).abs() < 1e-12 * (1.0 + fl[i].abs()),
            "component {i}: {} vs {}",
            exact[i],
            fl[i]
        );
    }

    let fixed = roe(left, right);
    let mut max_dev = 0.0_f64;
    for i in 0..3 {
        max_dev = max_dev.max((fixed[i] - fl[i]).abs() / fl[i].abs());
    }
    assert!(max_dev > 1e-4, "the entropy fix must act on |lambda| = 0");
    assert!(max_dev < 0.05, "entropy-fix smearing too large: {max_dev}");
}

/// Supersonic flow (all `λ_k` of one sign, outside the entropy-fix band) is
/// fully upwinded.
#[test]
fn supersonic_flow_is_fully_upwinded() {
    let left = [1.0, 3.0, 1.0];
    let right = [0.8, 3.2, 0.9];
    let flux = roe(left, right);
    let expected = physical_flux(left);
    for i in 0..3 {
        assert!(
            (flux[i] - expected[i]).abs() < 1e-12 * (1.0 + expected[i].abs()),
            "rightward component {i}: {} vs {}",
            flux[i],
            expected[i]
        );
    }

    let left = [1.0, -3.2, 1.0];
    let right = [0.8, -3.0, 0.9];
    let flux = roe(left, right);
    let expected = physical_flux(right);
    for i in 0..3 {
        assert!(
            (flux[i] - expected[i]).abs() < 1e-12 * (1.0 + expected[i].abs()),
            "leftward component {i}: {} vs {}",
            flux[i],
            expected[i]
        );
    }
}

/// Consistency: equal states return the physical flux.
#[test]
fn uniform_state_returns_the_physical_flux() {
    for state in [[1.0, 0.1, 1.0], [0.3, -2.0, 4.0], [2.5, 0.0, 0.7]] {
        let flux = roe(state, state);
        let expected = physical_flux(state);
        for i in 0..3 {
            assert!(
                (flux[i] - expected[i]).abs() < 1e-13 * (1.0 + expected[i].abs()),
                "state {state:?} component {i}"
            );
        }
    }
}

/// Mirroring the problem (swap sides, flip `u`) mirrors the flux: mass and
/// energy flux change sign, momentum flux does not.
#[test]
fn flux_is_mirror_symmetric() {
    for (name, left, right) in PROBLEMS {
        let mirrored_left = [right[0], -right[1], right[2]];
        let mirrored_right = [left[0], -left[1], left[2]];
        let f = roe(left, right);
        let g = roe(mirrored_left, mirrored_right);
        let tol = 1e-11 * (1.0 + max_abs(f));
        assert!((f[0] + g[0]).abs() <= tol, "{name}: mass flux");
        assert!((f[1] - g[1]).abs() <= tol, "{name}: momentum flux");
        assert!((f[2] + g[2]).abs() <= tol, "{name}: energy flux");
    }
}

/// Weak waves: every consistent approximate Riemann solver converges to the
/// linearised flux, so Roe and HLLC agree to a fraction of a percent.
#[test]
fn weak_waves_agree_with_hllc() {
    let cases = [
        ([1.0, 0.1, 1.0], [0.98, 0.12, 0.97]),
        ([1.0, 0.3, 1.0], [0.99, 0.31, 0.985]),
        ([1.0, 0.0, 1.0], [0.99, 0.0, 0.99]),
    ];
    for (left, right) in cases {
        let r = roe(left, right);
        let h = hllc_flux_euler_1d(
            left[0], left[1], left[2], right[0], right[1], right[2], GAMMA,
        );
        let scale = max_abs(h);
        for i in 0..3 {
            assert!(
                (r[i] - h[i]).abs() < 0.01 * scale,
                "{left:?}/{right:?} component {i}: roe {} vs hllc {}",
                r[i],
                h[i]
            );
        }
    }
}

/// At the Sod interface the exact Godunov flux is the flux of the left star
/// state (the rarefaction tail moves at -0.070 < 0 < u* = 0.927).  The Roe
/// mass flux must be within 2 % of it (observed: 1.2 %).
#[test]
fn sod_mass_flux_matches_exact_godunov_flux() {
    // Toro (2009) table 4.1: p* = 0.30313, u* = 0.92745.
    let (p_star, u_star) = (0.30313, 0.92745);
    let rho_star_left = 1.0 * (p_star / 1.0_f64).powf(1.0 / GAMMA);
    let exact_mass = rho_star_left * u_star;

    let flux = roe([1.0, 0.0, 1.0], [0.125, 0.0, 0.1]);
    assert!(flux.iter().all(|f| f.is_finite()));
    assert!(
        (flux[0] - exact_mass).abs() < 0.02 * exact_mass,
        "mass flux {} vs exact {exact_mass}",
        flux[0]
    );
    // The same interface through HLLC is finite and flows the same way.
    let hllc = hllc_flux_euler_1d(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, GAMMA);
    assert!(hllc[0] > 0.0 && flux[0] > 0.0);
}

/// Left-facing isentropic rarefaction that starts subsonic
/// (`u_L = u_over_a · a_L`) and ends supersonic (`u_R = mach_right · a_R`),
/// together with the exact sonic-point state that the interface `x/t = 0`
/// sees inside the fan.
fn transonic_rarefaction(u_over_a: f64, mach_right: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let (rho_l, p_l) = (1.0_f64, 1.0_f64);
    let a_l = (GAMMA * p_l / rho_l).sqrt();
    let u_l = u_over_a * a_l;
    // The Riemann invariant u + 2a/(γ-1) is constant across a 1-rarefaction.
    let invariant = u_l + 2.0 * a_l / (GAMMA - 1.0);
    let a_r = invariant / (mach_right + 2.0 / (GAMMA - 1.0));
    let right = [
        rho_l * (a_r / a_l).powf(2.0 / (GAMMA - 1.0)),
        mach_right * a_r,
        p_l * (a_r / a_l).powf(2.0 * GAMMA / (GAMMA - 1.0)),
    ];
    let a_s = (GAMMA - 1.0) / (GAMMA + 1.0) * invariant;
    let sonic = [
        rho_l * (a_s / a_l).powf(2.0 / (GAMMA - 1.0)),
        a_s,
        p_l * (a_s / a_l).powf(2.0 * GAMMA / (GAMMA - 1.0)),
    ];
    ([rho_l, u_l, p_l], right, sonic)
}

/// White-box check of Harten's fix: for `|λ₁| < δ = 0.1 ã` the first
/// acoustic wave contributes `(λ₁² + δ²)/(2δ)` instead of `|λ₁|`, the contact
/// is never touched, and the third wave is unchanged when it is supersonic.
#[test]
fn entropy_fix_replaces_small_acoustic_eigenvalues() {
    let (left, right, _) = transonic_rarefaction(0.9, 1.1);
    let dec = roe_decomposition(left, right, GAMMA).expect("valid Roe average");
    let delta = ROE_ENTROPY_FIX_FRACTION * dec.sound_speed;
    assert!(
        dec.lambda[0].abs() < delta,
        "test setup: lambda_1 = {} must lie inside the fix band {delta}",
        dec.lambda[0]
    );
    assert!(dec.lambda[2].abs() > delta);

    let fixed = roe_flux_euler_1d_with_entropy_fix(left, right, GAMMA, ROE_ENTROPY_FIX_FRACTION);
    let plain = roe_flux_euler_1d_with_entropy_fix(left, right, GAMMA, 0.0);
    let lambda_hat = (dec.lambda[0] * dec.lambda[0] + delta * delta) / (2.0 * delta);
    for i in 0..3 {
        let expected_shift =
            -0.5 * (lambda_hat - dec.lambda[0].abs()) * dec.alpha[0] * dec.eigenvectors[0][i];
        assert!(
            ((fixed[i] - plain[i]) - expected_shift).abs() < 1e-12 * (1.0 + plain[i].abs()),
            "component {i}: shift {} vs expected {expected_shift}",
            fixed[i] - plain[i]
        );
    }
    assert!(max_abs(sub(fixed, plain)) > 1e-4);

    // A contact with ũ = 0 (λ₂ = 0) must not be smoothed.
    let a = roe([1.0, 0.0, 1.0], [1.0, 0.0, 0.5]);
    let b = roe_flux_euler_1d_with_entropy_fix([1.0, 0.0, 1.0], [1.0, 0.0, 0.5], GAMMA, 0.0);
    for i in 0..3 {
        assert!((a[i] - b[i]).abs() < 1e-15, "contact smoothed in {i}");
    }
}

/// At a sonic rarefaction the unfixed Roe flux contains an expansion shock;
/// the fix pulls the flux onto the exact sonic-point flux.
#[test]
fn entropy_fix_improves_transonic_rarefaction_flux() {
    let (left, right, sonic) = transonic_rarefaction(0.9, 1.1);
    let exact = physical_flux(sonic);
    let scale = max_abs(exact);
    let error = |flux: [f64; 3]| max_abs(sub(flux, exact)) / scale;

    let fixed = error(roe(left, right));
    let plain = error(roe_flux_euler_1d_with_entropy_fix(left, right, GAMMA, 0.0));
    assert!(fixed < 1e-3, "fixed flux error {fixed}");
    assert!(
        fixed < plain / 5.0,
        "entropy fix must reduce the error: fixed {fixed} vs plain {plain}"
    );
}

/// A Roe average that is not a valid state (here `ã² < 0` from negative
/// pressures) falls back to the HLL flux instead of dividing by zero.
#[test]
fn invalid_roe_average_falls_back_to_hll() {
    let left = [1.0, 0.0, -0.5];
    let right = [1.0, 0.0, -0.5];
    assert!(roe_decomposition(left, right, GAMMA).is_none());
    let flux = roe(left, right);
    let hll = hll_flux_euler_1d(1.0, 0.0, -0.5, 1.0, 0.0, -0.5, GAMMA);
    assert!(flux.iter().all(|f| f.is_finite()));
    for i in 0..3 {
        assert_eq!(flux[i], hll[i]);
    }
}
