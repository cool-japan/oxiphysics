// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Bolton LCL temperature (the logarithm is taken of
//! the ratio `T / T_d` in kelvin) and the exact ISA density altitude (metres
//! throughout, no `120 ft/C` rule).

use super::*;

/// A one-level profile whose surface values are the arguments.
fn profile_with_surface(pressure: f64, temperature: f64, dew_point: f64) -> SoundingProfile {
    let mut profile = SoundingProfile::new("TEST");
    profile.add_level(SoundingLevel {
        pressure,
        height: 0.0,
        temperature,
        dew_point,
        wind_direction: 0.0,
        wind_speed: 0.0,
    });
    profile
}

// ── Bolton LCL ───────────────────────────────────────────────────────────

#[test]
fn bolton_lcl_temperature_matches_reference_value() {
    // T = 30 C, Td = 20 C: Bolton (1980) eq. 15 gives 290.815 K.
    let t_lcl = bolton_lcl_temperature(303.15, 293.15);
    assert!((t_lcl - 290.815).abs() < 0.01, "t_lcl={t_lcl}");
}

#[test]
fn bolton_lcl_temperature_of_saturated_air_is_the_air_temperature() {
    let t_lcl = bolton_lcl_temperature(298.15, 298.15);
    assert!((t_lcl - 298.15).abs() < 1e-9, "t_lcl={t_lcl}");
}

#[test]
fn bolton_lcl_temperature_clamps_dew_point_above_temperature() {
    let clamped = bolton_lcl_temperature(290.0, 295.0);
    assert!((clamped - 290.0).abs() < 1e-9, "clamped={clamped}");
}

#[test]
fn bolton_lcl_temperature_is_nan_for_non_physical_input() {
    assert!(bolton_lcl_temperature(f64::NAN, 280.0).is_nan());
    assert!(bolton_lcl_temperature(290.0, f64::INFINITY).is_nan());
    // A dew point at or below 56 K makes the formula singular.
    assert!(bolton_lcl_temperature(290.0, 56.0).is_nan());
    assert!(bolton_lcl_temperature(-5.0, -10.0).is_nan());
}

#[test]
fn lcl_pressure_for_warm_moist_surface() {
    // The old Celsius-difference logarithm gave ~221 hPa here; Bolton gives 864.5.
    let p_lcl = profile_with_surface(1000.0, 30.0, 20.0)
        .compute_lcl_pressure()
        .expect("physical surface values");
    assert!(p_lcl > 860.0 && p_lcl < 870.0, "p_lcl={p_lcl}");
}

#[test]
fn lcl_pressure_for_cooler_surface() {
    let p_lcl = profile_with_surface(1000.0, 20.0, 10.0)
        .compute_lcl_pressure()
        .expect("physical surface values");
    assert!(p_lcl > 857.0 && p_lcl < 866.0, "p_lcl={p_lcl}");
}

#[test]
fn lcl_pressure_of_saturated_surface_is_the_surface_pressure() {
    // T == Td used to evaluate ln(0) = -inf.
    let p_lcl = profile_with_surface(1000.0, 25.0, 25.0)
        .compute_lcl_pressure()
        .expect("saturated surface");
    assert!((p_lcl - 1000.0).abs() < 1e-9, "p_lcl={p_lcl}");
}

#[test]
fn lcl_pressure_clamps_dew_point_above_temperature() {
    let p_lcl = profile_with_surface(950.0, 20.0, 25.0)
        .compute_lcl_pressure()
        .expect("clamped to saturation");
    assert!((p_lcl - 950.0).abs() < 1e-9, "p_lcl={p_lcl}");
}

#[test]
fn lcl_pressure_decreases_with_dew_point_depression() {
    let moist = profile_with_surface(1000.0, 30.0, 25.0)
        .compute_lcl_pressure()
        .expect("moist");
    let medium = profile_with_surface(1000.0, 30.0, 20.0)
        .compute_lcl_pressure()
        .expect("medium");
    let dry = profile_with_surface(1000.0, 30.0, 10.0)
        .compute_lcl_pressure()
        .expect("dry");
    assert!(moist < 1000.0);
    assert!(medium < moist, "medium={medium} moist={moist}");
    assert!(dry < medium, "dry={dry} medium={medium}");
}

#[test]
fn lcl_height_agrees_with_the_espy_rule_of_thumb() {
    // The LCL rises ~125 m per degree of dew-point depression: ~1250 m for 10 K.
    let p_lcl = profile_with_surface(1000.0, 30.0, 20.0)
        .compute_lcl_pressure()
        .expect("physical surface values");
    let height = pressure_altitude(p_lcl, 1000.0);
    assert!(height > 1100.0 && height < 1400.0, "height={height}");
}

#[test]
fn lcl_pressure_is_none_without_a_physical_surface() {
    assert!(
        SoundingProfile::new("EMPTY")
            .compute_lcl_pressure()
            .is_none()
    );
    // -230 C is 43 K: below the 56 K singularity of the Bolton formula.
    assert!(
        profile_with_surface(1000.0, -220.0, -230.0)
            .compute_lcl_pressure()
            .is_none()
    );
}

#[test]
fn equivalent_potential_temperature_uses_the_kelvin_ratio_lcl() {
    let (t, td, p) = (300.0_f64, 295.0_f64, 1000.0_f64);
    let w = mixing_ratio(td - 273.15, p);
    let theta = potential_temperature(t, p, 1000.0);
    let t_lcl = 1.0 / (1.0 / (td - 56.0) + (t / td).ln() / 800.0) + 56.0;
    let expected = theta
        * (3.376 / t_lcl - 0.00254)
            .exp()
            .powf(w * 1000.0 * (1.0 + 0.81e-3 * w * 1000.0));
    let theta_e = equivalent_potential_temperature(t, td, p);
    assert!(
        (theta_e - expected).abs() < 1e-9,
        "theta_e={theta_e} expected={expected}"
    );
    // Moist air has the larger equivalent potential temperature.
    assert!(theta_e > equivalent_potential_temperature(t, 285.0, p));
    // Saturated air (T == Td) is finite, not ln(0).
    assert!(equivalent_potential_temperature(t, t, p).is_finite());
}

// ── Density altitude ─────────────────────────────────────────────────────

#[test]
fn density_altitude_at_isa_sea_level_is_zero() {
    assert!(density_altitude(0.0, 15.0).abs() < 1e-9);
}

#[test]
fn density_altitude_equals_pressure_altitude_at_isa_temperature() {
    for pressure_alt in [500.0, 1000.0, 2000.0, 3000.0] {
        let isa_temp = 15.0 - 0.0065 * pressure_alt;
        let da = density_altitude(pressure_alt, isa_temp);
        assert!((da - pressure_alt).abs() < 0.5, "pa={pressure_alt} da={da}");
    }
}

#[test]
fn density_altitude_correction_is_per_degree_in_metres_not_feet() {
    // 120 ft/C = 36.6 m/C; the exact ISA relation gives ~34.5 m/C here. The
    // old code added 120 *metres* per degree.
    let da = density_altitude(1000.0, 30.0);
    let dev = 30.0 - (15.0 - 0.0065 * 1000.0);
    let per_degree = (da - 1000.0) / dev;
    assert!(
        per_degree > 30.0 && per_degree < 40.0,
        "per_degree={per_degree}"
    );
}

#[test]
fn density_altitude_increases_with_temperature_and_can_go_negative() {
    let cold = density_altitude(1000.0, -10.0);
    let warm = density_altitude(1000.0, 10.0);
    let hot = density_altitude(1000.0, 40.0);
    assert!(cold < warm && warm < hot, "{cold} {warm} {hot}");
    // Very cold air at sea level is denser than the standard atmosphere.
    let arctic = density_altitude(0.0, -40.0);
    assert!(arctic < -2000.0 && arctic > -2500.0, "arctic={arctic}");
}

#[test]
fn density_altitude_hot_high_airfield() {
    // Pressure altitude 1600 m on a 30 C day: exact ISA value ~2469.6 m.
    let da = density_altitude(1600.0, 30.0);
    assert!(da > 2460.0 && da < 2480.0, "da={da}");
}

#[test]
fn density_altitude_is_nan_for_non_physical_input() {
    assert!(density_altitude(1000.0, -300.0).is_nan());
    // Above the ISA zero-pressure altitude (~44.3 km) the base is negative.
    assert!(density_altitude(50_000.0, 15.0).is_nan());
    assert!(density_altitude(f64::NAN, 15.0).is_nan());
    assert!(density_altitude(1000.0, f64::NAN).is_nan());
}

// ── Sounding text header ─────────────────────────────────────────────────

#[test]
fn sounding_text_header_without_an_elevation_value_does_not_panic() {
    // Six whitespace-separated fields: `parts[6]` used to index past the end.
    let text = "# SOUNDING 72451\n# LAT 40.78 LON -73.97 ELEV\n1000.0 10.0 20.0 15.0 180 10\n";
    let profile = SoundingProfile::from_text(text).expect("one level is present");
    assert!((profile.latitude - 40.78).abs() < 1e-12);
    assert!((profile.longitude + 73.97).abs() < 1e-12);
    assert_eq!(profile.elevation, 0.0);
    assert_eq!(profile.levels.len(), 1);
    // A header cut even shorter keeps the defaults.
    let text = "# SOUNDING 72451\n# LAT 40.78\n1000.0 10.0 20.0 15.0 180 10\n";
    let profile = SoundingProfile::from_text(text).expect("one level is present");
    assert_eq!(profile.latitude, 40.78);
    assert_eq!(profile.longitude, 0.0);
}

#[test]
fn sounding_text_round_trips_the_header() {
    let mut profile = SoundingProfile::new("72451");
    profile.latitude = 40.7812;
    profile.longitude = -73.9712;
    profile.elevation = 10.5;
    profile.add_level(SoundingLevel {
        pressure: 1000.0,
        height: 10.0,
        temperature: 20.0,
        dew_point: 15.0,
        wind_direction: 180.0,
        wind_speed: 10.0,
    });
    let parsed = SoundingProfile::from_text(&profile.to_text()).expect("round trip");
    assert!((parsed.latitude - 40.7812).abs() < 1e-9);
    assert!((parsed.longitude + 73.9712).abs() < 1e-9);
    assert!((parsed.elevation - 10.5).abs() < 1e-9);
}
