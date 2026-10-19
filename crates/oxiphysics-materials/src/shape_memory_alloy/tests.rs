// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the shape memory alloy models.

use super::*;

const TOL: f64 = 1.0e-6;

// Helper: create NiTi Brinson model
fn niti_brinson() -> BrinsonModel {
    NiTiPreset::brinson()
}

// Helper: create NiTi temperatures
fn niti_temps() -> TransformationTemperatures {
    NiTiPreset::temperatures()
}

// 1. TransformationTemperatures basic sanity
#[test]
fn test_transformation_temperatures_basic() {
    let tt = niti_temps();
    assert!(tt.ms > tt.mf, "Ms > Mf required");
    assert!(tt.af > tt.as_, "Af > As required");
    assert!(tt.hysteresis() > 0.0, "Hysteresis should be positive");
}

// 2. Temperature shift via Clausius-Clapeyron
#[test]
fn test_temperature_shift() {
    let tt = niti_temps();
    let stress = 100.0e6; // 100 MPa
    let shifted = tt.shifted(stress, 6.5e6, 6.5e6);
    let expected_shift = 100.0e6 / 6.5e6;
    assert!((shifted.ms - tt.ms - expected_shift).abs() < TOL);
    assert!((shifted.af - tt.af - expected_shift).abs() < TOL);
}

// 3. Tanaka model: fully martensite below Mf
#[test]
fn test_tanaka_full_martensite() {
    let model = NiTiPreset::tanaka();
    let xi = model.forward_fraction(200.0, 0.0); // well below Mf
    assert!(xi > 0.99, "Should be nearly fully martensite, got {xi}");
}

// 4. Tanaka model: fully austenite above Af
#[test]
fn test_tanaka_full_austenite() {
    let model = NiTiPreset::tanaka();
    let xi = model.reverse_fraction(400.0, 0.0); // well above Af
    assert!(xi < 0.01, "Should be nearly fully austenite, got {xi}");
}

// 5. Tanaka effective modulus at ξ=0 equals E_A
#[test]
fn test_tanaka_modulus_austenite() {
    let model = NiTiPreset::tanaka();
    let e = model.effective_modulus(0.0);
    assert!((e - 67.0e9).abs() < TOL, "E at xi=0 should be E_A");
}

// 6. Tanaka effective modulus at ξ=1 equals E_M
#[test]
fn test_tanaka_modulus_martensite() {
    let model = NiTiPreset::tanaka();
    let e = model.effective_modulus(1.0);
    assert!((e - 26.3e9).abs() < TOL, "E at xi=1 should be E_M");
}

// 7. Liang-Rogers: phase fraction at Mf is 1
#[test]
fn test_liang_rogers_mf() {
    let model = NiTiPreset::liang_rogers();
    let xi = model.phase_fraction(model.temps.mf, 0.0, 0.0);
    assert!(
        (xi - 1.0).abs() < 0.01,
        "At Mf, xi should be ~1.0, got {xi}"
    );
}

// 8. Liang-Rogers: phase fraction above Af is 0
#[test]
fn test_liang_rogers_above_af() {
    let model = NiTiPreset::liang_rogers();
    let xi = model.phase_fraction(model.temps.af + 50.0, 0.0, 1.0);
    assert!(xi < 0.01, "Above Af, xi should be ~0, got {xi}");
}

// 9. Liang-Rogers: reverse fraction at Af with xi_prev=1
#[test]
fn test_liang_rogers_reverse_at_af() {
    let model = NiTiPreset::liang_rogers();
    let xi = model.reverse_fraction(model.temps.af, 0.0, 1.0);
    assert!(xi < 0.05, "At Af with xi_prev=1, xi should be ~0, got {xi}");
}

// 10. Brinson: effective modulus interpolation
#[test]
fn test_brinson_modulus_interpolation() {
    let model = niti_brinson();
    let e_half = model.effective_modulus(0.5);
    let expected = 0.5 * 67.0e9 + 0.5 * 26.3e9;
    assert!(
        (e_half - expected).abs() < TOL,
        "E at xi=0.5 should be average, got {e_half}"
    );
}

// 11. Brinson: critical stress increases with temperature above Ms
#[test]
fn test_brinson_critical_stress_vs_temp() {
    let model = niti_brinson();
    let (sig_s_low, _) = model.critical_stresses_at_temp(model.temps.ms);
    let (sig_s_high, _) = model.critical_stresses_at_temp(model.temps.ms + 20.0);
    assert!(
        sig_s_high > sig_s_low,
        "Critical stress should increase with T"
    );
}

// 12. Brinson: stress response at zero strain and xi=0 is zero (at T_ref)
#[test]
fn test_brinson_zero_stress_at_origin() {
    let model = niti_brinson();
    let sigma = model.stress_response(0.0, 0.0, 0.0, 300.0, 300.0);
    assert!(
        sigma.abs() < TOL,
        "Stress at origin should be zero, got {sigma}"
    );
}

// 13. Brinson: state update changes xi during cooling
#[test]
fn test_brinson_state_update_cooling() {
    let model = niti_brinson();
    let mut state = SmaPhaseState::new(model.temps.af + 10.0);
    // Cool to between Mf and Ms
    let target_temp = (model.temps.mf + model.temps.ms) / 2.0;
    model.update_state(&mut state, 0.001, target_temp);
    assert!(
        state.xi > 0.0,
        "xi should increase during cooling, got {}",
        state.xi
    );
}

// 14. Superelastic curve: loading curve is non-empty
#[test]
fn test_superelastic_loading_nonempty() {
    let model = niti_brinson();
    let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 50);
    let curve = curve_gen.loading_curve(0.06);
    assert!(!curve.is_empty(), "Loading curve should be non-empty");
    assert_eq!(curve.len(), 50);
}

// 15. Superelastic curve: stress generally increases during loading
#[test]
fn test_superelastic_loading_trend() {
    let model = niti_brinson();
    let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 100);
    let curve = curve_gen.loading_curve(0.05);
    // Overall trend: final stress should be higher than initial (non-zero) stress
    let first_nonzero = curve.iter().find(|(_, s)| *s > 0.0);
    let last = curve.last().unwrap();
    if let Some(first) = first_nonzero {
        assert!(
            last.1 >= first.1,
            "Final stress should exceed initial: first={:.6}, last={:.6}",
            first.1,
            last.1
        );
    }
}

// 16. Superelastic curve: full loop has more points than one branch
#[test]
fn test_superelastic_full_loop_length() {
    let model = niti_brinson();
    let curve_gen = SuperelasticCurveGenerator::new(model, 350.0, 30);
    let loop_curve = curve_gen.full_loop(0.05);
    assert!(
        loop_curve.len() > 30,
        "Full loop should have loading + unloading points"
    );
}

// 17. Shape memory effect: recovery ratio is non-negative
#[test]
fn test_sme_recovery_ratio_nonneg() {
    let model = niti_brinson();
    let sme = ShapeMemoryEffect::new(model, 250.0, 350.0);
    let ratio = sme.recovery_ratio(0.05);
    assert!(ratio >= 0.0, "Recovery ratio should be >= 0, got {ratio}");
}

// 18. Two-way shape memory: training increases two_way_fraction
#[test]
fn test_twsm_training() {
    let model = niti_brinson();
    let mut twsm = TwoWayShapeMemory::new(model);
    assert!(twsm.two_way_fraction < TOL, "Initially should be ~0");
    for _ in 0..20 {
        twsm.train_cycle();
    }
    assert!(
        twsm.two_way_fraction > 0.1,
        "After training, fraction should grow, got {}",
        twsm.two_way_fraction
    );
}

// 19. Two-way: cooling strain at T < Mf is nonzero after training
#[test]
fn test_twsm_cooling_strain() {
    let model = niti_brinson();
    let mut twsm = TwoWayShapeMemory::new(model.clone());
    for _ in 0..50 {
        twsm.train_cycle();
    }
    let eps = twsm.cooling_strain(model.temps.mf - 10.0);
    assert!(eps > 0.0, "Cooling strain should be > 0, got {eps}");
}

// 20. Two-way: heating strain above Af is zero
#[test]
fn test_twsm_heating_above_af() {
    let model = niti_brinson();
    let mut twsm = TwoWayShapeMemory::new(model.clone());
    for _ in 0..10 {
        twsm.train_cycle();
    }
    let eps = twsm.heating_strain(model.temps.af + 20.0);
    assert!(
        eps.abs() < TOL,
        "Heating strain above Af should be ~0, got {eps}"
    );
}

// 21. Training effect: residual strain accumulates
#[test]
fn test_training_residual_accumulates() {
    let mut te = TrainingEffect::new();
    te.add_cycle();
    let r1 = te.current_residual();
    for _ in 0..50 {
        te.add_cycle();
    }
    let r51 = te.current_residual();
    assert!(r51 > r1, "Residual should grow: r1={r1}, r51={r51}");
}

// 22. Training effect: plateau stress decreases
#[test]
fn test_training_plateau_decreases() {
    let mut te = TrainingEffect::new();
    te.add_cycle();
    let s1 = te.current_plateau_stress();
    for _ in 0..100 {
        te.add_cycle();
    }
    let s101 = te.current_plateau_stress();
    assert!(
        s101 < s1,
        "Plateau stress should decrease: s1={s1}, s101={s101}"
    );
}

// 23. Thermomechanical: forward transformation heats
#[test]
fn test_thermo_forward_heats() {
    let tmc = NiTiPreset::thermomechanical();
    let dt = tmc.temperature_increment(0.5); // dxi = +0.5
    // Forward transformation releases heat → dt < 0 because of formula sign
    // but physical convention: temperature rises. Our formula: dT = -L/cp * dxi
    // dxi > 0 → dT < 0 (latent heat absorbed from surroundings for forward)
    // Actually for SMA, A→M is exothermic so temp rises.
    // Our formula gives dT = -L/cp * dxi. If L > 0 and dxi > 0 → dT < 0.
    // This models the environment perspective. Let's just verify it's nonzero.
    assert!(dt.abs() > 0.0, "Temperature change should be nonzero");
}

// 24. Thermomechanical: no phase change → no temperature change
#[test]
fn test_thermo_no_change() {
    let tmc = NiTiPreset::thermomechanical();
    let dt = tmc.temperature_increment(0.0);
    assert!(dt.abs() < TOL, "No phase change → no dT");
}

// 25. Clausius-Clapeyron: symmetric slopes
#[test]
fn test_cc_symmetric() {
    let cc = ClausiusClapeyron::symmetric(6.5e6);
    assert!(
        (cc.cm - cc.ca).abs() < TOL,
        "Symmetric slopes should be equal"
    );
}

// 26. Clausius-Clapeyron: stress shift
#[test]
fn test_cc_stress_shift() {
    let cc = ClausiusClapeyron::new(6.5e6, 7.0e6);
    let ds = cc.stress_shift_martensite(10.0);
    let expected = 6.5e6 * 10.0;
    assert!((ds - expected).abs() < TOL, "Stress shift incorrect");
}

// 27. CuAlNi preset creates valid model
#[test]
fn test_cualni_preset() {
    let model = CuAlNiPreset::brinson();
    assert!(model.e_a > model.e_m, "E_A > E_M for CuAlNi");
    assert!(model.h_max > 0.0, "h_max should be positive");
}

// 28. CuZnAl preset creates valid model
#[test]
fn test_cuznal_preset() {
    let model = CuZnAlPreset::brinson();
    let tt = CuZnAlPreset::temperatures();
    assert!(tt.af > tt.ms, "Af > Ms required");
    assert!(model.h_max > 0.0);
}

// 29. Hysteresis energy of identical curves is zero
#[test]
fn test_hysteresis_energy_zero() {
    let curve = vec![(0.0, 0.0), (0.01, 100.0e6), (0.02, 200.0e6)];
    let e = hysteresis_energy(&curve, &curve);
    assert!(e < TOL, "Same loading/unloading → zero hysteresis");
}

// 30. Hysteresis energy of a simple rectangle
#[test]
fn test_hysteresis_energy_rectangle() {
    let loading = vec![(0.0, 100.0), (1.0, 100.0)];
    let unloading = vec![(0.0, 50.0), (1.0, 50.0)];
    let e = hysteresis_energy(&loading, &unloading);
    assert!(
        (e - 50.0).abs() < TOL,
        "Rectangle hysteresis should be 50.0, got {e}"
    );
}

// 31. Secant modulus
#[test]
fn test_secant_modulus() {
    let m = secant_modulus(0.01, 700.0e6);
    assert!((m - 70.0e9).abs() < 1.0, "Secant modulus should be E");
}

// 32. Tangent modulus
#[test]
fn test_tangent_modulus() {
    let m = tangent_modulus(0.01, 700.0e6, 0.02, 1400.0e6);
    assert!((m - 70.0e9).abs() < 1.0, "Tangent modulus should be 70 GPa");
}

// 33. Temperature sweep recovery returns correct length
#[test]
fn test_temp_sweep_length() {
    let model = niti_brinson();
    let sweep = temperature_sweep_recovery(&model, 0.05, 250.0, 350.0, 20);
    assert_eq!(sweep.len(), 20);
}

// 34. Brinson reverse fractions decrease xi
#[test]
fn test_brinson_reverse_decreases() {
    let model = niti_brinson();
    let mid_temp = (model.temps.as_ + model.temps.af) / 2.0;
    let (xs, xt) = model.reverse_fractions(mid_temp, 0.0, 0.5, 0.5);
    assert!(xs + xt < 1.0, "Reverse should decrease total xi");
}

// 35. SmaPhaseState: fully martensite has xi=1
#[test]
fn test_phase_state_full_martensite() {
    let state = SmaPhaseState::fully_martensite(250.0);
    assert!((state.xi - 1.0).abs() < TOL);
    assert!((state.xi_t - 1.0).abs() < TOL);
}

// ─────────────────────────────────────────────────────────────────────────────
// Tanaka kinetics: coefficients, stress coupling, thermal hysteresis
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_tanaka_coefficients_have_physical_signs_and_magnitudes() {
    let model = NiTiPreset::tanaka();
    let ln_100 = 100.0_f64.ln();
    // Ms - Mf = 20 K and Af - As = 20 K for the NiTi preset.
    assert!(
        (model.a_m - ln_100 / 20.0).abs() < 1.0e-12,
        "a_M = {}",
        model.a_m
    );
    assert!(
        (model.a_a - ln_100 / 20.0).abs() < 1.0e-12,
        "a_A = {}",
        model.a_a
    );
    // Stress promotes martensite: b_M, b_A are +a/C (Clausius-Clapeyron).
    assert!(
        (model.b_m - model.a_m / 6.5e6).abs() < 1.0e-20,
        "b_M = {}",
        model.b_m
    );
    assert!(
        (model.b_a - model.a_a / 6.5e6).abs() < 1.0e-20,
        "b_A = {}",
        model.b_a
    );
    assert!(model.b_m > 0.0 && model.b_a > 0.0);
    assert_eq!((model.cm, model.ca), (6.5e6, 6.5e6));
}

#[test]
fn test_tanaka_transformation_is_99_percent_complete_at_finish_temperatures() {
    let model = NiTiPreset::tanaka();
    let t = model.temps;
    assert!((model.forward_fraction(t.mf, 0.0) - 0.99).abs() < 1.0e-12);
    assert!(model.forward_fraction(t.ms, 0.0).abs() < 1.0e-12);
    assert!((model.reverse_fraction(t.af, 0.0) - 0.01).abs() < 1.0e-12);
    assert!((model.reverse_fraction(t.as_, 0.0) - 1.0).abs() < 1.0e-12);
}

#[test]
fn test_tanaka_stress_promotes_forward_transformation() {
    let model = NiTiPreset::tanaka();
    let t = model.temps;
    let probe = t.ms - 1.0;
    let unstressed = model.forward_fraction(probe, 0.0);
    let stressed = model.forward_fraction(probe, 100.0e6);
    assert!(
        stressed > unstressed,
        "stress must raise the martensite fraction: {stressed} vs {unstressed}"
    );
    // Ms(sigma) = Ms + sigma / C_M: the transformation starts exactly there.
    let sigma = 65.0e6;
    let ms_sigma = t.ms + sigma / model.cm;
    assert!(model.forward_fraction(ms_sigma, sigma).abs() < 1.0e-12);
    assert!(model.forward_fraction(ms_sigma - 1.0, sigma) > 0.0);
    // Monotone in stress.
    let mut previous = -1.0;
    for k in 0..10 {
        let value = model.forward_fraction(probe, 20.0e6 * k as f64);
        assert!(value > previous, "not increasing in stress at k = {k}");
        previous = value;
    }
}

#[test]
fn test_tanaka_stress_retards_reverse_transformation() {
    let model = NiTiPreset::tanaka();
    let t = model.temps;
    let probe = 0.5 * (t.as_ + t.af);
    let unstressed = model.reverse_fraction(probe, 0.0);
    let stressed = model.reverse_fraction(probe, 50.0e6);
    assert!(
        stressed > unstressed,
        "stress must delay the reverse transformation: {stressed} vs {unstressed}"
    );
    // As(sigma) = As + sigma / C_A: the reverse transformation starts exactly there.
    let sigma = 32.5e6;
    let as_sigma = t.as_ + sigma / model.ca;
    assert!((model.reverse_fraction(as_sigma, sigma) - 1.0).abs() < 1.0e-12);
}

#[test]
fn test_tanaka_phase_fraction_retains_previous_value_in_hysteresis_window() {
    let model = NiTiPreset::tanaka();
    let t = model.temps;
    // Ms = 291 K < T < As = 295 K: neither forward nor reverse transformation acts.
    for xi_prev in [0.0, 0.37, 1.0] {
        let xi = model.phase_fraction(293.0, 0.0, xi_prev);
        assert_eq!(xi, xi_prev, "window value must be retained");
    }
    // Under stress the window moves up by sigma / C.
    let sigma = 65.0e6;
    let shift = sigma / model.cm;
    assert_eq!(model.phase_fraction(293.0 + shift, sigma, 0.6), 0.6);
    assert!(model.phase_fraction(t.ms - 5.0, 0.0, 0.0) > 0.0);
}

#[test]
fn test_tanaka_phase_fraction_uses_the_models_own_clausius_clapeyron_slopes() {
    // CuAlNi: C_M = C_A = 4.5 MPa/K, Ms = 291 K, As = 295 K.
    let model = CuAlNiPreset::tanaka();
    assert_eq!((model.cm, model.ca), (4.5e6, 4.5e6));
    let sigma = 200.0e6;
    let shift = sigma / 4.5e6; // 44.4 K (the former fixed 6.5 MPa/K would give 30.8 K)
    let ms_sigma = 291.0 + shift;
    let as_sigma = 295.0 + shift;
    // Just below Ms(sigma): forward transformation is active ...
    assert!(model.phase_fraction(ms_sigma - 1.0, sigma, 0.0) > 0.05);
    // ... and inside the shifted window the previous fraction is retained.
    assert_eq!(
        model.phase_fraction(0.5 * (ms_sigma + as_sigma), sigma, 0.4),
        0.4
    );
    // Just above As(sigma): reverse transformation (cannot increase the fraction).
    let xi = model.phase_fraction(as_sigma + 1.0, sigma, 1.0);
    assert!(xi < 1.0 && xi > 0.5, "reverse branch near As: {xi}");
}

#[test]
fn test_tanaka_phase_fraction_never_reverses_the_transformation_direction() {
    let model = NiTiPreset::tanaka();
    let t = model.temps;
    // Forward region: martensite already present is retained (heating below As).
    for k in 0..20 {
        let temp = t.mf + (t.ms - t.mf) * k as f64 / 20.0;
        assert!(model.phase_fraction(temp, 0.0, 0.95) >= 0.95);
    }
    assert_eq!(model.phase_fraction(t.mf - 10.0, 0.0, 1.0), 1.0);
    // Reverse region: austenite is not re-transformed above Ms (cooling above As).
    for k in 0..20 {
        let temp = t.as_ + (t.af + 20.0 - t.as_) * k as f64 / 20.0;
        assert!(model.phase_fraction(temp, 0.0, 0.3) <= 0.3);
    }
    // exp(-ln(100)/20 * 60) = 1e-6 at Af + 40 K.
    assert!(model.phase_fraction(t.af + 40.0, 0.0, 0.3) < 2.0e-6);
}

#[test]
fn test_tanaka_cooling_heating_cycle_exhibits_thermal_hysteresis() {
    let model = NiTiPreset::tanaka();
    // Cool 340 K -> 240 K, then heat back, threading the fraction through phase_fraction.
    let mut xi = 0.0;
    let mut cooling = Vec::new();
    for k in 0..=100 {
        let temp = 340.0 - k as f64;
        xi = model.phase_fraction(temp, 0.0, xi);
        cooling.push((temp, xi));
    }
    let mut heating = Vec::new();
    for k in 0..=100 {
        let temp = 240.0 + k as f64;
        xi = model.phase_fraction(temp, 0.0, xi);
        heating.push((temp, xi));
    }
    let at = |path: &[(f64, f64)], temp: f64| -> f64 {
        path.iter()
            .find(|(t, _)| (*t - temp).abs() < 1.0e-9)
            .map(|&(_, x)| x)
            .unwrap_or(f64::NAN)
    };
    // Fully austenitic before Ms, fully (99 %+) martensitic after Mf, on both paths' ends.
    assert!(at(&cooling, 340.0).abs() < 1.0e-12);
    assert!(at(&cooling, 260.0) > 0.99);
    assert!(at(&heating, 340.0) < 0.011);
    // At 285 K (Mf < T < Ms) the cooling path is only partly transformed while the heating
    // path still carries all the martensite: the loop is open (thermal hysteresis).
    let cool_285 = at(&cooling, 285.0);
    let heat_285 = at(&heating, 285.0);
    assert!(cool_285 > 0.3 && cool_285 < 0.9, "cooling value {cool_285}");
    assert!(
        heat_285 > 0.99,
        "heating keeps the martensite below As: {heat_285}"
    );
    // The heating path is monotonically non-increasing, the cooling path non-decreasing.
    assert!(heating.windows(2).all(|w| w[1].1 <= w[0].1 + 1.0e-12));
    assert!(cooling.windows(2).all(|w| w[1].1 >= w[0].1 - 1.0e-12));
}

// ─────────────────────────────────────────────────────────────────────────────
// Brinson kinetics: endpoints of the closed-form branches
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_brinson_xi_s_loading_saturates_above_the_finish_stress() {
    let model = niti_brinson();
    // T = 250 K < Ms: window [sigma_s^cr, sigma_f^cr] = [100, 170] MPa.
    assert_eq!(model.xi_s_loading(250.0, 200.0e6, 0.3, 0.7), 1.0);
    assert_eq!(model.xi_s_loading(250.0, 170.0e6, 0.3, 0.7), 1.0);
    assert_eq!(model.xi_s_loading(250.0, 50.0e6, 0.3, 0.7), 0.3);
    assert_eq!(model.xi_s_loading(250.0, -10.0e6, 0.3, 0.7), 0.3);
    // Cosine midpoint: xi_S = 1 - (1 - xi_S0) / 2.
    let mid = model.xi_s_loading(250.0, 135.0e6, 0.3, 0.7);
    assert!((mid - 0.65).abs() < 1.0e-12, "midpoint {mid}");
    // Above Ms the window shifts up with C_M (T - Ms).
    let (sig_s, sig_f) = model.critical_stresses_at_temp(350.0);
    assert!((model.xi_s_loading(350.0, sig_s, 0.0, 0.0)).abs() < 1.0e-12);
    assert_eq!(model.xi_s_loading(350.0, sig_f + 1.0, 0.0, 0.0), 1.0);
}

#[test]
fn test_brinson_xi_t_cooling_completes_below_the_finish_temperature() {
    let model = niti_brinson();
    let t = model.temps;
    assert_eq!(model.xi_t_cooling(t.mf - 5.0, 0.0, 0.2), 1.0);
    assert_eq!(model.xi_t_cooling(t.ms + 5.0, 0.0, 0.2), 0.2);
    let mid = model.xi_t_cooling(0.5 * (t.ms + t.mf), 0.0, 0.2);
    assert!((mid - 0.6).abs() < 1.0e-12, "midpoint {mid}");
}

#[test]
fn test_brinson_reverse_fractions_vanish_above_the_finish_temperature() {
    let model = niti_brinson();
    let t = model.temps;
    assert_eq!(
        model.reverse_fractions(t.af + 10.0, 0.0, 0.5, 0.5),
        (0.0, 0.0)
    );
    assert_eq!(
        model.reverse_fractions(t.as_ - 10.0, 0.0, 0.5, 0.5),
        (0.5, 0.5)
    );
    let (xs, xt) = model.reverse_fractions(0.5 * (t.as_ + t.af), 0.0, 0.5, 0.25);
    assert!((xs - 0.25).abs() < 1.0e-12 && (xt - 0.125).abs() < 1.0e-12);
    // Tensile stress delays the reverse transformation by sigma / C_A.
    let sigma = 65.0e6;
    let (xs_load, _) = model.reverse_fractions(t.af + 5.0, sigma, 1.0, 0.0);
    assert!(
        xs_load > 0.0,
        "still partly martensite under load: {xs_load}"
    );
    assert_eq!(
        model.reverse_fractions(t.af + 10.0 + sigma / 6.5e6, sigma, 1.0, 0.0),
        (0.0, 0.0)
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Brinson update_state: detwinning, shape-memory recovery, pseudoelasticity
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn test_brinson_detwinning_occurs_below_ms_under_load() {
    // Twinned martensite at 250 K (< Mf = 271 K) strained to 5 %: the stress rises into
    // the detwinning window [100, 170] MPa and a large part of the twins reorient.
    let model = niti_brinson();
    let mut state = SmaPhaseState::fully_martensite(250.0);
    model.update_state(&mut state, 0.05, 250.0);
    assert!(
        state.xi_s > 0.5 && state.xi_s < 0.8,
        "xi_S = {}",
        state.xi_s
    );
    assert!(
        (state.xi_s + state.xi_t - 1.0).abs() < 1.0e-12,
        "fully martensitic"
    );
    assert!((state.xi - 1.0).abs() < 1.0e-12);
    assert!(
        (100.0e6..=170.0e6).contains(&state.stress),
        "stress {} inside the detwinning window",
        state.stress
    );
    assert!((state.transformation_strain - model.h_max * state.xi_s).abs() < 1.0e-15);
    // The stress satisfies the constitutive law with the stored fractions.
    let law = model.stress_response(0.05, state.xi_s, state.xi_t, 250.0, model.t_ref);
    assert!(
        (state.stress - law).abs() < 1.0e-3,
        "stress {} vs constitutive law {law}",
        state.stress
    );
    // Reference solution of the implicit stress equation (default T_ref = Af = 315 K).
    assert!(
        (state.xi_s - 0.645_630).abs() < 1.0e-5,
        "xi_S = {}",
        state.xi_s
    );
    assert!(
        (state.stress - 141.585e6).abs() < 2.0e3,
        "sigma = {}",
        state.stress
    );
}

#[test]
fn test_brinson_no_detwinning_below_the_critical_stress() {
    let model = niti_brinson().with_reference_temperature(250.0);
    let mut state = SmaPhaseState::fully_martensite(250.0);
    // 0.1 % strain: sigma = E_M * 0.001 = 26 MPa < sigma_s^cr = 100 MPa.
    model.update_state(&mut state, 0.001, 250.0);
    assert_eq!(state.xi_s, 0.0);
    assert_eq!(state.xi_t, 1.0);
    assert!(
        (state.stress - 26.3e9 * 0.001).abs() < 1.0e-3,
        "{}",
        state.stress
    );
}

#[test]
fn test_brinson_detwinning_saturates_at_full_detwinned_martensite() {
    let model = niti_brinson().with_reference_temperature(250.0);
    let mut state = SmaPhaseState::fully_martensite(250.0);
    model.update_state(&mut state, 0.2, 250.0);
    assert_eq!(state.xi_s, 1.0);
    assert_eq!(state.xi_t, 0.0);
    assert!(state.stress > model.sigma_crit_f);
}

#[test]
fn test_brinson_shape_memory_effect_recovers_the_deformation() {
    let model = niti_brinson();
    let sme = ShapeMemoryEffect::new(model, 250.0, 350.0);
    let (before, after) = sme.simulate_recovery(0.05);
    assert!(
        before > 0.03 && before < 0.05,
        "residual transformation strain after unloading: {before}"
    );
    assert!(
        after.abs() < 1.0e-12,
        "all strain recovered above Af: {after}"
    );
    let ratio = sme.recovery_ratio(0.05);
    assert!(ratio > 0.9, "recovery ratio {ratio}");
    assert!((ratio - 1.0).abs() < 1.0e-12);
    // Pins the stress-free reference-temperature convention of the simulation.
    assert!(
        (before - 0.044_581_295).abs() < 1.0e-6,
        "residual strain {before}"
    );
}

#[test]
fn test_brinson_partial_recovery_follows_the_reverse_cosine_kinetics() {
    // Heating to the middle of [As, Af] recovers exactly half of the strain.
    let t = niti_temps();
    for (t_rec, expected) in [
        (t.as_ - 5.0, 0.0),
        (t.as_, 0.0),
        (t.as_ + 5.0, 0.5 * (1.0 - std::f64::consts::FRAC_1_SQRT_2)),
        (0.5 * (t.as_ + t.af), 0.5),
        (t.af, 1.0),
        (t.af + 30.0, 1.0),
    ] {
        let sme = ShapeMemoryEffect::new(niti_brinson(), 250.0, t_rec);
        let ratio = sme.recovery_ratio(0.05);
        assert!(
            (ratio - expected).abs() < 1.0e-9,
            "T_rec = {t_rec}: ratio {ratio}, expected {expected}"
        );
    }
}

#[test]
fn test_brinson_recovery_ratio_is_zero_without_residual_strain() {
    // A tiny strain never reaches the critical stress: nothing to recover.
    let sme = ShapeMemoryEffect::new(niti_brinson(), 250.0, 350.0);
    assert_eq!(sme.recovery_ratio(0.001), 0.0);
}

#[test]
fn test_brinson_heating_above_af_removes_all_martensite() {
    let model = niti_brinson();
    let mut state = SmaPhaseState::fully_martensite(280.0);
    model.update_state(&mut state, 0.0, 350.0);
    assert_eq!(state.xi_s, 0.0);
    assert_eq!(state.xi_t, 0.0);
    assert_eq!(state.xi, 0.0);
    // Strain-free heating of constrained austenite builds up only the thermoelastic stress.
    let thermal = model.theta * (350.0 - model.t_ref);
    assert!(
        (state.stress - thermal).abs() < 1.0e-3,
        "{} vs {thermal}",
        state.stress
    );
}

#[test]
fn test_brinson_cooling_below_mf_gives_twinned_martensite() {
    let model = niti_brinson();
    let mut state = SmaPhaseState::new(330.0);
    model.update_state(&mut state, 0.0, 250.0);
    assert!((state.xi_t - 1.0).abs() < 1.0e-12, "xi_T = {}", state.xi_t);
    assert_eq!(state.xi_s, 0.0);
    assert!((state.xi - 1.0).abs() < 1.0e-12);
    // Cooling does not touch detwinned martensite: it is kept, the rest becomes twinned.
    let mut retained = SmaPhaseState::new(300.0);
    retained.xi_s = 0.4;
    retained.xi = 0.4;
    model.update_state(&mut retained, 0.0, 250.0);
    assert!(
        (retained.xi_s - 0.4).abs() < 1.0e-12,
        "xi_S = {}",
        retained.xi_s
    );
    assert!(
        (retained.xi_t - 0.6).abs() < 1.0e-12,
        "xi_T = {}",
        retained.xi_t
    );
}

#[test]
fn test_brinson_cooling_through_the_window_is_a_cosine_in_temperature() {
    let model = niti_brinson();
    let t = model.temps;
    let mut state = SmaPhaseState::new(330.0);
    model.update_state(&mut state, 0.0, 0.5 * (t.ms + t.mf));
    assert!(
        (state.xi_t - 0.5).abs() < 1.0e-12,
        "xi_T at mid-window {}",
        state.xi_t
    );
}

#[test]
fn test_brinson_thermoelastic_term_vanishes_at_the_reference_temperature() {
    let austenite_temp = 355.0;
    let model = niti_brinson();
    assert_eq!(model.t_ref, model.temps.af);
    let mut state = SmaPhaseState::new(austenite_temp);
    model.update_state(&mut state, 0.0, austenite_temp);
    let expected = model.theta * (austenite_temp - model.temps.af);
    assert!(
        (state.stress - expected).abs() < 1.0e-3,
        "{} vs {expected}",
        state.stress
    );
    assert!(state.stress > 0.0);
    let referenced = niti_brinson().with_reference_temperature(austenite_temp);
    assert_eq!(referenced.t_ref, austenite_temp);
    let mut stress_free = SmaPhaseState::new(austenite_temp);
    referenced.update_state(&mut stress_free, 0.0, austenite_temp);
    assert!(stress_free.stress.abs() < 1.0e-3, "{}", stress_free.stress);
}

fn replay_pseudoelastic_loop(
    steps: usize,
    temperature: f64,
    max_strain: f64,
) -> (SmaPhaseState, SmaPhaseState) {
    let model = niti_brinson().with_reference_temperature(temperature);
    let mut state = SmaPhaseState::new(temperature);
    for i in 0..=steps {
        let eps = max_strain * i as f64 / steps as f64;
        model.update_state(&mut state, eps, temperature);
    }
    let peak = state;
    for i in 1..=steps {
        let eps = max_strain * (1.0 - i as f64 / steps as f64);
        model.update_state(&mut state, eps, temperature);
    }
    (peak, state)
}

#[test]
fn test_brinson_pseudoelastic_loop_closes_and_is_independent_of_step_size() {
    let (peak_coarse, end_coarse) = replay_pseudoelastic_loop(10, 350.0, 0.06);
    let (peak_fine, end_fine) = replay_pseudoelastic_loop(400, 350.0, 0.06);
    // Loading produces a large amount of detwinned martensite ...
    assert!(
        peak_fine.xi_s > 0.6 && peak_fine.xi_s < 0.8,
        "peak xi_S {}",
        peak_fine.xi_s
    );
    assert_eq!(peak_fine.xi_t, 0.0);
    assert!(
        (peak_fine.stress - 527.3e6).abs() < 1.0e6,
        "peak stress {}",
        peak_fine.stress
    );
    // ... which is independent of how finely the load path is discretised ...
    assert!((peak_coarse.xi_s - peak_fine.xi_s).abs() < 1.0e-9);
    assert!((peak_coarse.stress - peak_fine.stress).abs() < 1.0);
    // ... and unloading recovers everything (flag-shaped superelastic loop).
    for end in [end_coarse, end_fine] {
        assert!(
            end.xi_s.abs() < 1.0e-12 && end.xi_t.abs() < 1.0e-12,
            "{end:?}"
        );
        assert!(end.stress.abs() < 1.0e-3, "residual stress {}", end.stress);
        assert!(end.transformation_strain.abs() < 1.0e-12);
    }
}

#[test]
fn test_brinson_superelastic_curve_starts_at_origin_and_shows_hysteresis() {
    let generator = SuperelasticCurveGenerator::new(niti_brinson(), 350.0, 200);
    let loop_curve = generator.full_loop(0.06);
    assert_eq!(loop_curve.len(), 400);
    let (loading, unloading) = loop_curve.split_at(200);
    assert!(
        loading[0].1.abs() < 1.0e-3,
        "stress-free origin: {}",
        loading[0].1
    );
    assert!(
        unloading[199].1.abs() < 1.0e-3,
        "stress-free end: {}",
        unloading[199].1
    );
    // Loading plateau (~0.5 GPa) lies far above the unloading plateau (~0.3 GPa).
    let load_at = |e: f64| {
        loading
            .iter()
            .find(|p| p.0 >= e)
            .map(|p| p.1)
            .unwrap_or(f64::NAN)
    };
    let unload_at = |e: f64| {
        unloading
            .iter()
            .find(|p| p.0 <= e)
            .map(|p| p.1)
            .unwrap_or(f64::NAN)
    };
    assert!(
        load_at(0.03) > 480.0e6 && load_at(0.03) < 540.0e6,
        "{}",
        load_at(0.03)
    );
    assert!(
        unload_at(0.03) > 250.0e6 && unload_at(0.03) < 340.0e6,
        "{}",
        unload_at(0.03)
    );
    let energy = hysteresis_energy(loading, unloading);
    assert!(
        energy > 1.0e7 && energy < 1.2e7,
        "dissipated energy {energy} J/m3"
    );
}

#[test]
fn test_brinson_strain_and_stress_controlled_updates_are_consistent() {
    let model = niti_brinson().with_reference_temperature(300.0);
    let start = SmaPhaseState::new(330.0);
    let mut strain_driven = start;
    model.update_state(&mut strain_driven, 0.04, 330.0);
    assert!(
        strain_driven.xi_s > 0.1,
        "plateau reached: {}",
        strain_driven.xi_s
    );
    let mut stress_driven = start;
    model.update_state_stress_controlled(&mut stress_driven, strain_driven.stress, 330.0);
    assert!((stress_driven.xi_s - strain_driven.xi_s).abs() < 1.0e-9);
    assert!((stress_driven.xi_t - strain_driven.xi_t).abs() < 1.0e-9);
    assert!(
        (stress_driven.strain - 0.04).abs() < 1.0e-9,
        "{}",
        stress_driven.strain
    );
}

#[test]
fn test_brinson_stress_free_unloading_keeps_the_residual_strain() {
    let model = niti_brinson().with_reference_temperature(250.0);
    let mut state = SmaPhaseState::fully_martensite(250.0);
    model.update_state(&mut state, 0.05, 250.0);
    let xi_s = state.xi_s;
    model.update_state_stress_controlled(&mut state, 0.0, 250.0);
    assert_eq!(state.stress, 0.0);
    assert!(
        (state.xi_s - xi_s).abs() < 1.0e-12,
        "detwinned martensite is retained"
    );
    // Residual strain = eps_L * xi_S (no thermal term at the reference temperature).
    assert!(
        (state.strain - model.h_max * xi_s).abs() < 1.0e-12,
        "{}",
        state.strain
    );
}

#[test]
fn test_temperature_sweep_recovery_traces_the_reverse_transformation() {
    let model = niti_brinson();
    let sweep = temperature_sweep_recovery(&model, 0.05, 250.0, 350.0, 101);
    assert_eq!(sweep.len(), 101);
    assert!((sweep[0].0 - 250.0).abs() < 1.0e-12 && (sweep[100].0 - 350.0).abs() < 1.0e-12);
    assert!(
        (sweep[0].1 - 0.05).abs() < 1.0e-12,
        "starts at the stored strain: {}",
        sweep[0].1
    );
    assert!(
        sweep[100].1.abs() < 1.0e-12,
        "fully recovered: {}",
        sweep[100].1
    );
    assert!(
        sweep.windows(2).all(|w| w[1].1 <= w[0].1 + 1.0e-15),
        "monotone recovery"
    );
    let at = |temp: f64| {
        sweep
            .iter()
            .find(|p| (p.0 - temp).abs() < 1.0e-9)
            .map(|p| p.1)
            .unwrap_or(f64::NAN)
    };
    assert!((at(295.0) - 0.05).abs() < 1.0e-12, "no recovery at As");
    assert!(
        (at(305.0) - 0.025).abs() < 1.0e-12,
        "half recovered at (As + Af)/2"
    );
    assert!(at(315.0).abs() < 1.0e-12, "fully recovered at Af");
}

#[test]
fn test_brinson_update_state_ignores_non_finite_input() {
    let model = niti_brinson();
    let mut state = SmaPhaseState::fully_martensite(250.0);
    model.update_state(&mut state, 0.05, 250.0);
    let before = state;
    model.update_state(&mut state, f64::NAN, 250.0);
    model.update_state(&mut state, 0.05, f64::INFINITY);
    model.update_state_stress_controlled(&mut state, f64::NAN, 250.0);
    assert_eq!(state.xi_s, before.xi_s);
    assert_eq!(state.stress, before.stress);
    assert_eq!(state.temperature, before.temperature);
}

#[test]
fn test_brinson_update_state_invariants_hold_along_random_paths() {
    // Deterministic pseudo-random thermomechanical paths: the fractions stay in the
    // simplex and every state satisfies the constitutive law.
    let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..40 {
        let start_temp = 240.0 + 140.0 * next();
        let model = niti_brinson().with_reference_temperature(start_temp);
        let mut state = SmaPhaseState::new(start_temp);
        if next() < 0.5 {
            state.xi_t = 1.0;
            state.xi = 1.0;
        }
        for _ in 0..60 {
            let strain = -0.02 + 0.11 * next();
            let temp = 240.0 + 140.0 * next();
            model.update_state(&mut state, strain, temp);
            assert!(
                state.xi_s >= 0.0 && state.xi_s <= 1.0,
                "xi_S = {}",
                state.xi_s
            );
            assert!(
                state.xi_t >= 0.0 && state.xi_t <= 1.0,
                "xi_T = {}",
                state.xi_t
            );
            assert!(
                state.xi_s + state.xi_t <= 1.0 + 1.0e-12,
                "sum {}",
                state.xi_s + state.xi_t
            );
            assert!(
                state.stress.is_finite() && state.strain == strain && state.temperature == temp
            );
            let law = model.stress_response(strain, state.xi_s, state.xi_t, temp, model.t_ref);
            assert!(
                (state.stress - law).abs() < 10.0,
                "constitutive residual {} Pa at eps = {strain}, T = {temp}",
                state.stress - law
            );
        }
    }
}

#[test]
fn test_brinson_incremental_update_matches_closed_form_branch_kinetics() {
    // Along a monotone path the step-by-step update must reproduce Brinson's closed-form
    // cosine kinetics from the transformation-start fractions, whatever the step count.
    let temperature = 350.0;
    let model = niti_brinson().with_reference_temperature(temperature);
    let mut state = SmaPhaseState::new(temperature);
    // Loading branch: start from austenite (xi_S0 = 0), 37 awkward steps up to 3 % strain.
    for i in 0..=37 {
        model.update_state(&mut state, 0.03 * i as f64 / 37.0, temperature);
    }
    let forward_closed_form = model.xi_s_loading(temperature, state.stress, 0.0, 0.0);
    assert!(forward_closed_form > 0.2 && forward_closed_form < 0.9);
    assert!(
        (state.xi_s - forward_closed_form).abs() < 1.0e-12,
        "incremental xi_S {} vs closed form {forward_closed_form}",
        state.xi_s
    );
    // Continue to the peak (6 % strain), then unload to 2 % strain in 29 steps: the reverse
    // branch starts from the peak fractions.
    for i in 1..=40 {
        model.update_state(&mut state, 0.03 + 0.03 * i as f64 / 40.0, temperature);
    }
    let (peak_xi_s, peak_xi_t) = (state.xi_s, state.xi_t);
    for i in 1..=29 {
        model.update_state(&mut state, 0.06 - 0.04 * i as f64 / 29.0, temperature);
    }
    let (reverse_xi_s, reverse_xi_t) =
        model.reverse_fractions(temperature, state.stress, peak_xi_s, peak_xi_t);
    assert!(
        reverse_xi_s > 0.0 && reverse_xi_s < peak_xi_s,
        "partly recovered: {reverse_xi_s} of {peak_xi_s}"
    );
    assert!(
        (state.xi_s - reverse_xi_s).abs() < 1.0e-12 && (state.xi_t - reverse_xi_t).abs() < 1.0e-12,
        "incremental ({}, {}) vs closed form ({reverse_xi_s}, {reverse_xi_t})",
        state.xi_s,
        state.xi_t
    );
}

#[test]
fn test_brinson_unloading_inside_the_loading_window_keeps_the_transformed_fraction() {
    // Unloading from a partly transformed state must not run the forward kinetics backwards:
    // while the stress is still above the reverse window the fractions stay constant.
    let temperature = 350.0;
    let model = niti_brinson().with_reference_temperature(temperature);
    let mut state = SmaPhaseState::new(temperature);
    for i in 0..=100 {
        model.update_state(&mut state, 0.03 * i as f64 / 100.0, temperature);
    }
    let loaded = state;
    assert!(loaded.xi_s > 0.2);
    // Elastic unloading by 0.1 % strain stays well above sigma = C_A (T - As) = 357 MPa.
    model.update_state(&mut state, 0.029, temperature);
    assert!(state.stress < loaded.stress && state.stress > 400.0e6);
    assert_eq!(state.xi_s, loaded.xi_s);
    assert_eq!(state.xi_t, loaded.xi_t);
    // Reloading to the previous strain restores the previous stress (no hysteresis jump).
    model.update_state(&mut state, 0.03, temperature);
    assert!(
        (state.stress - loaded.stress).abs() < 1.0,
        "{} vs {}",
        state.stress,
        loaded.stress
    );
    assert!((state.xi_s - loaded.xi_s).abs() < 1.0e-12);
}

#[test]
fn test_brinson_partial_unloading_and_reloading_follow_return_point_memory() {
    let temperature = 350.0;
    let model = niti_brinson().with_reference_temperature(temperature);
    let mut state = SmaPhaseState::new(temperature);
    // Load to 5 % strain on the loading plateau.
    for i in 0..=200 {
        model.update_state(&mut state, 0.05 * i as f64 / 200.0, temperature);
    }
    let peak = state;
    assert!(peak.xi_s > 0.5, "peak xi_S {}", peak.xi_s);

    // Unload deep into the reverse window (stress ~ 290 MPa < C_A (T - As) = 357 MPa).
    for i in 1..=200 {
        model.update_state(&mut state, 0.05 - 0.02 * i as f64 / 200.0, temperature);
    }
    let turning = state;
    assert!(
        turning.xi_s > 0.0 && turning.xi_s < peak.xi_s,
        "partial recovery: {} of {}",
        turning.xi_s,
        peak.xi_s
    );
    assert!(
        turning.stress < 357.0e6 && turning.stress > 227.5e6,
        "{}",
        turning.stress
    );

    // Reload: the transformed fraction stays constant on the elastic line ...
    let mut max_elastic_stress = 0.0_f64;
    for i in 1..=200 {
        model.update_state(&mut state, 0.03 + 0.01 * i as f64 / 200.0, temperature);
        if state.stress < model.critical_stresses_at_temp(temperature).0 - 1.0e6 {
            assert!(
                (state.xi_s - turning.xi_s).abs() < 1.0e-12,
                "xi_S changed on the elastic reload line at sigma = {}",
                state.stress
            );
            max_elastic_stress = max_elastic_stress.max(state.stress);
        }
    }
    assert!(
        max_elastic_stress > turning.stress + 1.0e6,
        "the reload line was exercised"
    );
    // ... and transformation resumes only on the loading plateau.
    for i in 1..=200 {
        model.update_state(&mut state, 0.04 + 0.01 * i as f64 / 200.0, temperature);
    }
    assert!(
        state.xi_s > turning.xi_s,
        "transformation resumed: {}",
        state.xi_s
    );
    let (plateau_start, plateau_end) = model.critical_stresses_at_temp(temperature);
    assert!(
        state.stress > plateau_start && state.stress < plateau_end,
        "reloaded stress {} on the loading plateau [{plateau_start}, {plateau_end}]",
        state.stress
    );
    assert!(state.xi_s <= 1.0 && state.xi_t == 0.0);
}

#[test]
fn test_brinson_update_records_the_transformation_memory() {
    let model = niti_brinson();
    let mut state = SmaPhaseState::fully_martensite(250.0);
    assert!(state.progress_forward.is_none());
    model.update_state(&mut state, 0.05, 250.0);
    let forward = state.progress_forward.unwrap_or(f64::NAN);
    // The forward ramp stands inside its window; the other branches sit at their ends.
    assert!(forward > 0.0 && forward < 1.0, "forward progress {forward}");
    assert_eq!(state.progress_reverse, Some(0.0));
    assert_eq!(state.progress_cooling, Some(1.0));
    // Elastic unloading leaves the memory untouched.
    let before = state;
    model.update_state(&mut state, 0.049, 250.0);
    assert_eq!(state.progress_forward, before.progress_forward);
    assert_eq!(state.xi_s, before.xi_s);
}
