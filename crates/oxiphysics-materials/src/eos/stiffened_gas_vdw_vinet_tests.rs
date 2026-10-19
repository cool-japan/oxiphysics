//! Regression tests for the stiffened-gas, van der Waals and Vinet equations of state.
//!
//! Reference values quoted as literals were obtained independently of the
//! implementation (closed forms evaluated by hand and 40-digit polynomial root
//! finding for the van der Waals cubic).

use super::functions::EquationOfState;
use super::types::{StiffenedGasEos, VanDerWaalsEos, VinetEos};

fn relative_error(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

// ---------------------------------------------------------------------------
// Stiffened gas
// ---------------------------------------------------------------------------

#[test]
fn stiffened_gas_water_is_at_one_bar_and_1628_m_per_s_at_reference_density() {
    let water = StiffenedGasEos::water();
    assert_eq!(water.reference_pressure, 1.0e5);
    let p = water.pressure(1000.0);
    assert!(relative_error(p, 1.0e5) < 1e-12, "P(1000) = {p}");
    let c = water.sound_speed(1000.0);
    let expected = (6.59_f64 * (1.0e5 + 4.02e8) / 1000.0).sqrt();
    assert!((c - expected).abs() < 1e-9, "c = {c}, expected {expected}");
    assert!((c - 1628.0).abs() < 0.5, "c(1000) = {c}");
}

#[test]
fn stiffened_gas_old_hardcoded_energy_branch_gave_6_gpa_for_water() {
    // pressure() used to evaluate (gamma-1) rho (cv T) - p_inf with cv T = 1.2e6,
    // which is ~6.3 GPa for water at its reference density.
    let p = StiffenedGasEos::water().pressure(1000.0);
    assert!(
        p < 1.0e6,
        "reference-density pressure must be ~1 bar, got {p}"
    );
}

#[test]
fn stiffened_gas_air_is_standard_atmosphere_with_340_m_per_s() {
    let air = StiffenedGasEos::air();
    assert_eq!(air.reference_pressure, 101_325.0);
    assert!(relative_error(air.pressure(1.225), 101_325.0) < 1e-12);
    let c = air.sound_speed(1.225);
    assert!((c - 340.29).abs() < 0.01, "c = {c}");
}

#[test]
fn stiffened_gas_new_defaults_to_one_atmosphere_reference_pressure() {
    let eos = StiffenedGasEos::new(2.0, 1.0e6, 500.0, 0.0);
    assert_eq!(
        eos.reference_pressure,
        StiffenedGasEos::DEFAULT_REFERENCE_PRESSURE
    );
    assert!(relative_error(eos.pressure(500.0), 101_325.0) < 1e-12);
}

#[test]
fn stiffened_gas_with_reference_pressure_moves_the_isentrope() {
    let eos = StiffenedGasEos::water().with_reference_pressure(5.0e6);
    assert!(relative_error(eos.pressure(1000.0), 5.0e6) < 1e-12);
    let c = eos.sound_speed(1000.0);
    let expected = (6.59_f64 * (5.0e6 + 4.02e8) / 1000.0).sqrt();
    assert!((c - expected).abs() < 1e-9);
}

#[test]
fn stiffened_gas_density_pressure_round_trip_is_exact_to_1e_12() {
    for eos in [StiffenedGasEos::water(), StiffenedGasEos::air()] {
        for factor in [0.5, 0.9, 1.0, 1.1, 2.0, 5.0] {
            let rho = eos.reference_density * factor;
            let rho_back = eos.density_from_pressure(eos.pressure(rho));
            assert!(
                relative_error(rho_back, rho) < 1e-12,
                "gamma = {}: {rho} -> {rho_back}",
                eos.gamma
            );
        }
    }
}

#[test]
fn stiffened_gas_round_trip_near_the_cavitation_limit_is_limited_by_conditioning() {
    // At rho = 0.1 rho_ref the water pressure is -p_inf + ~100 Pa: forming P + p_inf
    // cancels ~7 digits, so the round trip is good to ~1e-9, not to machine precision.
    let water = StiffenedGasEos::water();
    for factor in [0.1, 0.2, 0.3] {
        let rho = water.reference_density * factor;
        let rho_back = water.density_from_pressure(water.pressure(rho));
        assert!(relative_error(rho_back, rho) < 1e-9, "{rho} -> {rho_back}");
    }
}

#[test]
fn stiffened_gas_isentrope_is_monotone_in_pressure_and_sound_speed() {
    let water = StiffenedGasEos::water();
    let mut previous_p = f64::NEG_INFINITY;
    let mut previous_c = 0.0;
    for i in 1..=60 {
        let rho = 100.0 * f64::from(i);
        let p = water.pressure(rho);
        let c = water.sound_speed(rho);
        assert!(p > previous_p && c > previous_c, "rho = {rho}");
        previous_p = p;
        previous_c = c;
    }
}

#[test]
fn stiffened_gas_energy_form_has_gamma_times_p_inf_offset() {
    let water = StiffenedGasEos::water();
    let p_at_reference_energy = water.pressure_from_energy(1000.0, water.reference_energy);
    let expected = -6.59 * 4.02e8;
    assert!(relative_error(p_at_reference_energy, expected) < 1e-14);
    // The old form subtracted p_inf instead of gamma * p_inf: (gamma - 1) p_inf = 2.2 GPa apart.
    assert!((p_at_reference_energy - (-4.02e8)).abs() > 2.0e9);
    // p_inf = 0 is the ideal gas: p = (gamma - 1) rho e.
    let air = StiffenedGasEos::air();
    let p = air.pressure_from_energy(1.225, 2.5e5);
    assert!(relative_error(p, 0.4 * 1.225 * 2.5e5) < 1e-14);
}

#[test]
fn stiffened_gas_energy_form_matches_isentrope_and_sound_speed_identity() {
    // On the isentrope e - e0 = (P + gamma p_inf) / ((gamma - 1) rho); the energy
    // form must return the isentrope pressure and the thermodynamic identity
    // c^2 = (gamma - 1) [ (e - e0) + P / rho ] must equal gamma (P + p_inf) / rho.
    for eos in [StiffenedGasEos::water(), StiffenedGasEos::air()] {
        for factor in [0.8, 1.0, 1.3, 2.0] {
            let rho = eos.reference_density * factor;
            let p = eos.pressure(rho);
            let e = eos.reference_energy + (p + eos.gamma * eos.p_inf) / ((eos.gamma - 1.0) * rho);
            let p_energy = eos.pressure_from_energy(rho, e);
            assert!(
                (p_energy - p).abs() <= 1e-9 * (p.abs() + eos.p_inf),
                "energy form {p_energy} vs isentrope {p}"
            );
            let c_squared = (eos.gamma - 1.0) * ((e - eos.reference_energy) + p / rho);
            let c = eos.sound_speed(rho);
            assert!(
                relative_error(c * c, c_squared) < 1e-12,
                "c^2 = {} vs {c_squared}",
                c * c
            );
        }
    }
}

#[test]
fn stiffened_gas_zero_density_limit_is_minus_p_inf() {
    let water = StiffenedGasEos::water();
    assert!(relative_error(water.pressure(0.0), -4.02e8) < 1e-14);
    assert!(relative_error(water.pressure(-3.0), -4.02e8) < 1e-14);
    assert_eq!(water.sound_speed(0.0), 0.0);
    assert_eq!(water.density_from_pressure(-4.02e8), 0.0);
    assert_eq!(water.density_from_pressure(-9.0e8), 0.0);
    assert!(water.density_from_pressure(f64::NAN).is_nan());
}

// ---------------------------------------------------------------------------
// Vinet
// ---------------------------------------------------------------------------

#[test]
fn vinet_sound_speed_at_reference_density_is_sqrt_k0_over_rho0() {
    let copper = VinetEos::copper();
    let rho0 = 1.0 / copper.v0;
    let c = copper.sound_speed(rho0);
    let expected = (136.7e9_f64 / 8960.0).sqrt();
    assert!(relative_error(c, expected) < 1e-9, "c = {c}");
    assert!((c - 3906.0).abs() < 0.5);
    assert!(relative_error(copper.bulk_modulus(copper.v0), 136.7e9) < 1e-12);
}

#[test]
fn vinet_sound_speed_is_not_identically_zero_and_has_golden_values() {
    // Golden values from the closed form K = K0 x^-2 exp(eta(1-x)) [2 - x + eta x (1-x)].
    let copper = VinetEos::copper();
    let rho0 = 1.0 / copper.v0;
    for (ratio, expected) in [
        (0.9, 3017.1819434599315),
        (1.05, 4316.300766966462),
        (1.2, 5454.830239840407),
        (1.5, 7442.268002696007),
    ] {
        let c = copper.sound_speed(ratio * rho0);
        assert!(
            relative_error(c, expected) < 1e-9,
            "rho/rho0 = {ratio}: c = {c}"
        );
    }
}

#[test]
fn vinet_sound_speed_matches_central_difference_of_the_pressure() {
    for eos in [VinetEos::copper(), VinetEos::diamond()] {
        let rho0 = 1.0 / eos.v0;
        for ratio in [0.9, 1.0, 1.05, 1.2, 1.5] {
            let rho = ratio * rho0;
            let h = 1e-5 * rho;
            let slope = (eos.pressure(rho + h) - eos.pressure(rho - h)) / (2.0 * h);
            let c = eos.sound_speed(rho);
            assert!(
                relative_error(c * c, slope) < 1e-5,
                "rho/rho0 = {ratio}: c^2 = {} vs dP/drho = {slope}",
                c * c
            );
        }
    }
}

#[test]
fn vinet_sound_speed_increases_monotonically_with_compression() {
    for eos in [VinetEos::copper(), VinetEos::diamond()] {
        let rho0 = 1.0 / eos.v0;
        let mut previous = 0.0;
        for i in 0..=100 {
            let rho = rho0 * (1.0 + 0.005 * f64::from(i));
            let c = eos.sound_speed(rho);
            assert!(
                c > previous,
                "rho/rho0 = {}: c = {c} after {previous}",
                rho / rho0
            );
            previous = c;
        }
    }
}

#[test]
fn vinet_sound_speed_degenerate_densities_are_zero_not_nan() {
    let copper = VinetEos::copper();
    assert_eq!(copper.sound_speed(0.0), 0.0);
    assert_eq!(copper.sound_speed(-1.0), 0.0);
    // Far into expansion the modulus is negative; the sound speed is clamped to 0.
    assert_eq!(copper.sound_speed(1.0e-3), 0.0);
}

// ---------------------------------------------------------------------------
// Van der Waals
// ---------------------------------------------------------------------------

/// Critical point of the model: `T_c = 8a/(27 R b)`, `P_c = a/(27 b^2)`.
fn critical_state(eos: &VanDerWaalsEos) -> (f64, f64) {
    (
        8.0 * eos.a / (27.0 * eos.r_gas * eos.b),
        eos.a / (27.0 * eos.b * eos.b),
    )
}

#[test]
fn vdw_co2_at_350_k_and_one_bar_matches_polynomial_root_reference() {
    let co2 = VanDerWaalsEos::co2();
    let rho = co2.density_at_temperature(1.0e5, 350.0);
    assert!(
        relative_error(rho, 1.516403796615695) < 1e-12,
        "rho = {rho}"
    );
    // Ideal gas gives 1.512 kg/m^3; the finding quotes ~1.52.
    assert!((rho - 1.52).abs() < 0.01);
    let via_trait = co2.with_temperature(350.0).density_from_pressure(1.0e5);
    assert_eq!(rho, via_trait);
}

#[test]
fn vdw_density_from_pressure_no_longer_saturates_at_the_upper_bracket() {
    // The bisection over [1e-5, 1e4] kg/m^3 used to return ~10000 for any target.
    let co2 = VanDerWaalsEos::co2();
    let rho = co2.density_from_pressure(1.0e5);
    assert!(
        relative_error(rho, 1.7714877921987426) < 1e-12,
        "rho = {rho}"
    );
    assert!(rho < 10.0);
    let rho_high = co2.density_from_pressure(1.0e8);
    assert!(
        relative_error(rho_high, 815.5306450976467) < 1e-12,
        "rho = {rho_high}"
    );
    assert!(rho_high < co2.molar_mass / co2.b);
}

#[test]
fn vdw_pressure_round_trip_holds_across_the_supercritical_isotherm() {
    let co2 = VanDerWaalsEos::co2().with_temperature(350.0);
    for pressure in [1.0e3, 1.0e5, 1.0e6, 5.0e6, 2.0e7, 1.0e8, 5.0e8] {
        let rho = co2.density_from_pressure(pressure);
        let back = co2.pressure(rho);
        assert!(
            relative_error(back, pressure) < 1e-9,
            "P = {pressure}: rho = {rho}, back = {back}"
        );
    }
    // The isotherm is monotone above T_c, so the density must rise with pressure.
    let mut previous = 0.0;
    for exponent in 3..=9 {
        let rho = co2.density_from_pressure(10.0_f64.powi(exponent));
        assert!(rho > previous);
        previous = rho;
    }
}

#[test]
fn vdw_returns_the_vapour_root_inside_the_three_root_window() {
    // (T, P, vapour density): P lies between the spinodal pressures, so the cubic
    // has three real roots; the other two densities are ~348/628 kg/m^3 (T = 250 K),
    // ~347/532 (280 K) and ~333/423 (300 K).
    let co2 = VanDerWaalsEos::co2();
    let critical_density = co2.molar_mass / (3.0 * co2.b);
    for (temperature, pressure, expected) in [
        (250.0, 2.0e6, 49.681850225606944),
        (280.0, 5.0e6, 147.26570596193957),
        (300.0, 7.0e6, 270.03494256221506),
    ] {
        let rho = co2.density_at_temperature(pressure, temperature);
        assert!(
            relative_error(rho, expected) < 1e-10,
            "T = {temperature}: rho = {rho}, expected {expected}"
        );
        assert!(
            rho < critical_density,
            "must be the gas-like root, got {rho}"
        );
        let back = co2.pressure_at_temperature(rho, temperature);
        assert!(relative_error(back, pressure) < 1e-9);
        // No root at a lower density: below the vapour root the isotherm lies strictly
        // under the target pressure.
        let mut lower = rho * (1.0 - 1e-6);
        while lower > 1e-6 {
            assert!(co2.pressure_at_temperature(lower, temperature) < pressure);
            lower *= 0.7;
        }
    }
    // The trait uses the configured temperature (default 300 K).
    assert_eq!(co2.temperature, 300.0);
    let rho_trait = co2.density_from_pressure(7.0e6);
    assert!(relative_error(rho_trait, 270.03494256221506) < 1e-10);
}

#[test]
fn vdw_temperature_field_drives_pressure_and_density() {
    let default = VanDerWaalsEos::co2();
    assert_eq!(default.temperature, VanDerWaalsEos::DEFAULT_TEMPERATURE);
    let hot = default.with_temperature(350.0);
    assert_eq!(hot.temperature, 350.0);
    let rho = 25.0;
    assert_eq!(hot.pressure(rho), hot.pressure_at_temperature(rho, 350.0));
    assert_eq!(
        default.pressure(rho),
        default.pressure_at_temperature(rho, 300.0)
    );
    assert!(hot.pressure(rho) > default.pressure(rho));
    assert!(hot.density_from_pressure(1.0e5) < default.density_from_pressure(1.0e5));
}

#[test]
fn vdw_pressure_is_infinite_at_and_beyond_the_covolume_density() {
    let co2 = VanDerWaalsEos::co2();
    let pole = co2.molar_mass / co2.b;
    assert_eq!(co2.pressure(pole), f64::INFINITY);
    assert_eq!(co2.pressure(1.2 * pole), f64::INFINITY);
    assert_eq!(co2.sound_speed(1.2 * pole), f64::INFINITY);
    let just_below = co2.pressure(0.99 * pole);
    assert!(
        just_below.is_finite() && just_below > 1.0e8,
        "P = {just_below}"
    );
    assert_eq!(co2.pressure(0.0), 0.0);
    assert_eq!(co2.pressure(-1.0), 0.0);
}

#[test]
fn vdw_non_positive_pressure_has_no_finite_vapour_volume() {
    let co2 = VanDerWaalsEos::co2();
    assert_eq!(co2.density_from_pressure(0.0), 0.0);
    assert_eq!(co2.density_from_pressure(-1.0e5), 0.0);
    assert!(co2.vapour_molar_volume(0.0, 300.0).is_infinite());
    assert!(co2.density_at_temperature(f64::NAN, 300.0).is_nan());
    assert!(co2.density_at_temperature(1.0e5, -1.0).is_nan());
}

#[test]
fn vdw_ideal_gas_limit_without_attraction_or_covolume() {
    let ideal = VanDerWaalsEos::new(0.0, 0.0, 0.028);
    let rho = ideal.density_at_temperature(1.0e5, 300.0);
    let expected = 1.0e5 * 0.028 / (8.314 * 300.0);
    assert!(relative_error(rho, expected) < 1e-12, "rho = {rho}");
}

#[test]
fn vdw_critical_point_has_the_triple_root_three_b() {
    let co2 = VanDerWaalsEos::co2();
    let (t_c, p_c) = critical_state(&co2);
    assert!((t_c - 303.9).abs() < 0.1);
    let v = co2.vapour_molar_volume(p_c, t_c);
    assert!(
        relative_error(v, 3.0 * co2.b) < 1e-4,
        "V = {v}, 3b = {}",
        3.0 * co2.b
    );
}

#[test]
fn vdw_sound_speed_is_the_isothermal_slope_at_the_configured_temperature() {
    let co2 = VanDerWaalsEos::co2().with_temperature(350.0);
    for rho in [0.5, 5.0, 50.0, 200.0] {
        let h = 1e-6 * rho;
        let slope = (co2.pressure(rho + h) - co2.pressure(rho - h)) / (2.0 * h);
        let c = co2.sound_speed(rho);
        assert!(
            relative_error(c * c, slope) < 1e-6,
            "rho = {rho}: c^2 = {} vs {slope}",
            c * c
        );
    }
    // Ideal-gas isothermal limit sqrt(R T / M) as rho -> 0, and it grows with T.
    let c_dilute = co2.sound_speed(1.0e-6);
    assert!(relative_error(c_dilute, (8.314_f64 * 350.0 / 0.044).sqrt()) < 1e-3);
    assert!(c_dilute > VanDerWaalsEos::co2().sound_speed(1.0e-6));
}

#[test]
fn vdw_sound_speed_vanishes_inside_the_spinodal_region() {
    // At 250 K the spinodal densities are M/V for V in (8.38e-5, 2.35e-4) m^3/mol,
    // i.e. rho in (187, 525) kg/m^3, where (dP/drho)_T < 0.
    let co2 = VanDerWaalsEos::co2().with_temperature(250.0);
    assert_eq!(co2.sound_speed(350.0), 0.0);
    assert!(co2.pressure_density_derivative_at_temperature(350.0, 250.0) < 0.0);
    assert!(co2.sound_speed(50.0) > 100.0);
}
