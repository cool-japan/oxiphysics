// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Langmuir equilibrium constant derived from the
//! adsorption energy.
//!
//! `LangmuirAdsorption::equilibrium_constant_from_thermo` used to evaluate
//! `exp(+E_ads / k_B T)` for an adsorption energy `E_ads < 0`, so stronger
//! binding gave an exponentially *smaller* constant.

use super::*;

/// One electronvolt in joules.
const EV_TO_J: f64 = 1.602_176_634e-19;

/// Mass of an N₂ molecule (kg).
const N2_MASS_KG: f64 = 28.0134 * 1.660_539_066_60e-27;

/// Thermodynamic equilibrium constant of N₂ at `temperature_k` for the binding
/// energy `energy_ev` (eV, negative for exothermic adsorption).
fn equilibrium_constant(energy_ev: f64, temperature_k: f64) -> f64 {
    LangmuirAdsorption::equilibrium_constant_from_thermo(
        energy_ev * EV_TO_J,
        temperature_k,
        N2_MASS_KG,
    )
}

/// `λ³ / k_B T` with `λ = h / sqrt(2π m k_B T)`: the constant at zero binding energy.
fn inverse_translational_partition_function(temperature_k: f64) -> f64 {
    let kbt = KB * temperature_k;
    let thermal_wavelength = H_PLANCK / (2.0 * PI * N2_MASS_KG * kbt).sqrt();
    thermal_wavelength.powi(3) / kbt
}

#[test]
fn stronger_binding_gives_a_larger_equilibrium_constant() {
    let strong = equilibrium_constant(-0.5, 300.0);
    let weak = equilibrium_constant(-0.1, 300.0);
    assert!(weak > 0.0 && weak.is_finite(), "K(-0.1 eV) = {weak}");
    assert!(strong.is_finite(), "K(-0.5 eV) = {strong}");
    assert!(
        strong > weak,
        "K(-0.5 eV) = {strong} must exceed K(-0.1 eV) = {weak}"
    );
}

#[test]
fn equilibrium_constant_is_monotone_in_the_binding_energy() {
    let mut previous = 0.0;
    // Binding energy deepening from -0.05 eV to -1.0 eV.
    for step in 1..=20 {
        let energy_ev = -0.05 * step as f64;
        let k = equilibrium_constant(energy_ev, 300.0);
        assert!(
            k > previous,
            "K must grow as the binding deepens: E = {energy_ev} eV gives {k} after {previous}"
        );
        previous = k;
    }
}

#[test]
fn equilibrium_constant_ratio_is_the_boltzmann_factor_of_the_energy_difference() {
    let (e_strong_ev, e_weak_ev) = (-0.5, -0.1);
    let temperature_k = 300.0;
    let k_strong = equilibrium_constant(e_strong_ev, temperature_k);
    let k_weak = equilibrium_constant(e_weak_ev, temperature_k);
    // ln(K1/K2) = (E2 - E1) / k_B T = 0.4 eV / k_B T, about 15.5.
    let expected = (e_weak_ev - e_strong_ev) * EV_TO_J / (KB * temperature_k);
    let actual = (k_strong / k_weak).ln();
    assert!(expected > 15.0 && expected < 16.0, "expected = {expected}");
    assert!(
        (actual - expected).abs() <= 1.0e-9 * expected.abs().max(1.0),
        "ln(K1/K2) = {actual}, expected {expected}"
    );
}

#[test]
fn zero_binding_energy_gives_the_inverse_translational_partition_function() {
    for temperature_k in [77.0, 300.0, 1000.0] {
        let k = equilibrium_constant(0.0, temperature_k);
        let expected = inverse_translational_partition_function(temperature_k);
        assert!(
            (k - expected).abs() <= 1.0e-12 * expected,
            "T = {temperature_k} K: K = {k:e}, lambda^3/kT = {expected:e}"
        );
    }
}

#[test]
fn endothermic_adsorption_is_weaker_than_the_unbound_reference() {
    let reference = equilibrium_constant(0.0, 300.0);
    let endothermic = equilibrium_constant(0.2, 300.0);
    assert!(
        endothermic < reference,
        "a repulsive site must bind less than E = 0: {endothermic} vs {reference}"
    );
}

#[test]
fn equilibrium_constant_obeys_the_van_t_hoff_temperature_law() {
    // K = λ³ / (k_B T) exp(-E / k_B T) ∝ T^(-5/2) exp(-E / k_B T), so for E < 0
    //   ln(K(T1) / K(T2)) = (5/2) ln(T2 / T1) + (-E / k_B) (1/T1 - 1/T2).
    let energy_ev = -0.3;
    let (t1, t2) = (250.0, 350.0);
    let actual = (equilibrium_constant(energy_ev, t1) / equilibrium_constant(energy_ev, t2)).ln();
    let expected = 2.5 * (t2 / t1).ln() + (-energy_ev * EV_TO_J / KB) * (1.0 / t1 - 1.0 / t2);
    assert!(
        (actual - expected).abs() <= 1.0e-9 * expected.abs().max(1.0),
        "ln(K1/K2) = {actual}, expected {expected}"
    );
    // Exothermic binding weakens when the surface is heated.
    assert!(
        equilibrium_constant(energy_ev, 200.0) > equilibrium_constant(energy_ev, 300.0)
            && equilibrium_constant(energy_ev, 300.0) > equilibrium_constant(energy_ev, 400.0)
    );
}

#[test]
fn coverage_from_the_thermodynamic_constant_rises_with_binding_strength() {
    let pressure_pa = 1.0e5;
    let coverage = |energy_ev: f64| {
        LangmuirAdsorption::new(equilibrium_constant(energy_ev, 300.0), 100.0, 1.0e5)
            .langmuir_coverage(pressure_pa)
    };
    let strong = coverage(-0.5);
    let weak = coverage(-0.1);
    // K(-0.5 eV) is about 4e-4 Pa^-1 and K(-0.1 eV) about 8e-11 Pa^-1.
    assert!(strong > 0.9, "strongly bound surface: coverage {strong}");
    assert!(weak < 1.0e-3, "weakly bound surface: coverage {weak}");
}

#[test]
fn extreme_binding_energy_saturates_instead_of_overflowing() {
    // 50 eV is about 1900 k_B T at 300 K: exp(1900) overflows an f64.
    let k = equilibrium_constant(-50.0, 300.0);
    assert!(k.is_finite(), "K must stay finite, got {k}");
    assert!(k > 1.0e100, "K = {k:e}");
    let model = LangmuirAdsorption::new(k, 100.0, 1.0e5);
    let theta = model.langmuir_coverage(1.0e5);
    assert_eq!(theta, 1.0, "saturated surface is fully covered");
}

#[test]
fn langmuir_coverage_saturates_when_the_product_k_p_overflows() {
    let model = LangmuirAdsorption::new(f64::MAX, 100.0, 1.0e5);
    assert_eq!(model.langmuir_coverage(1.0e5), model.theta_max);
    // Zero pressure still gives an empty surface.
    assert_eq!(model.langmuir_coverage(0.0), 0.0);
}

#[test]
fn equilibrium_constant_degenerate_inputs() {
    assert_eq!(equilibrium_constant(-0.5, 0.0), 0.0);
    assert_eq!(
        LangmuirAdsorption::equilibrium_constant_from_thermo(-EV_TO_J, 300.0, 0.0),
        0.0
    );
    assert!(equilibrium_constant(f64::NAN, 300.0).is_nan());
}

#[test]
fn arrhenius_factor_is_the_boltzmann_factor_of_the_barrier() {
    let temperature_k = 400.0;
    let kbt = KB * temperature_k;
    let barrier_j = 0.7 * EV_TO_J;
    let factor = arrhenius_factor(barrier_j, temperature_k);
    let expected = (-barrier_j / kbt).exp();
    assert!((factor - expected).abs() <= 1.0e-12 * expected);
    // A negative barrier is an energy gain: the factor exceeds one.
    assert!(arrhenius_factor(-barrier_j, temperature_k) > 1.0);
    // The two are exact reciprocals.
    let product = factor * arrhenius_factor(-barrier_j, temperature_k);
    assert!((product - 1.0).abs() < 1.0e-12, "product = {product}");
}

#[test]
fn arrhenius_factor_saturates_at_both_extremes() {
    let temperature_k = 300.0;
    let kbt = KB * temperature_k;
    // exponent = -barrier / kbt below -700: flushed to zero.
    assert_eq!(arrhenius_factor(800.0 * kbt, temperature_k), 0.0);
    // exponent above 700: finite instead of infinity.
    assert_eq!(arrhenius_factor(-800.0 * kbt, temperature_k), f64::MAX);
    // Just inside the range the exponential is evaluated normally.
    let inside = arrhenius_factor(-600.0 * kbt, temperature_k);
    assert!(
        inside.is_finite() && inside > 1.0e260,
        "inside = {inside:e}"
    );
    // Zero temperature has no thermal activation at all.
    assert_eq!(arrhenius_factor(-1.0e-19, 0.0), 0.0);
}
