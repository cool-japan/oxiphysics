//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_new {
    use super::*;
    use std::f64::consts::PI;
    #[test]
    fn test_sound_speed_ideal_gas_air() {
        let c = sound_speed_ideal_gas(1.4, 287.0, 293.0);
        assert!(c > 340.0 && c < 350.0, "air sound speed ~343 m/s, got {c}");
    }
    #[test]
    fn test_sound_speed_ideal_gas_increases_with_temp() {
        let c_cold = sound_speed_ideal_gas(1.4, 287.0, 250.0);
        let c_warm = sound_speed_ideal_gas(1.4, 287.0, 350.0);
        assert!(c_warm > c_cold, "sound speed increases with temperature");
    }
    #[test]
    fn test_sound_speed_liquid_water() {
        let c = sound_speed_liquid(2.18e9, 998.0);
        assert!(
            c > 1400.0 && c < 1600.0,
            "water sound speed ~1480 m/s, got {c}"
        );
    }
    #[test]
    fn test_sound_speed_air_temperature_at_0c() {
        let c = sound_speed_air_temperature(0.0);
        assert!((c - 331.3).abs() < 0.01, "c(0°C) = 331.3 m/s, got {c}");
    }
    #[test]
    fn test_sound_speed_air_temperature_at_20c() {
        let c = sound_speed_air_temperature(20.0);
        assert!(c > 342.0 && c < 345.0, "c(20°C) ≈ 343 m/s, got {c}");
    }
    #[test]
    fn test_group_velocity_linear_dispersion() {
        let c = 343.0_f64;
        let omega_fn = |k: f64| c * k;
        let vg = group_velocity(10.0, 0.001, omega_fn);
        assert!(
            (vg - c).abs() < 0.001,
            "group velocity = c for linear dispersion, got {vg}"
        );
    }
    #[test]
    fn test_pressure_from_intensity_round_trip() {
        let z = 413.0_f64;
        let i = 1.0e-3;
        let p_peak = pressure_from_intensity(i, z);
        let i_back = p_peak * p_peak / (2.0 * z);
        assert!(
            (i_back - i).abs() / i < 1e-10,
            "round-trip I→p→I: {i_back} vs {i}"
        );
    }
    #[test]
    fn test_rms_pressure_is_peak_over_sqrt2() {
        let p_peak = 10.0;
        let p_rms = rms_pressure(p_peak);
        assert!(
            (p_rms - p_peak / 2.0_f64.sqrt()).abs() < 1e-12,
            "p_rms = p_peak/√2"
        );
    }
    #[test]
    fn test_particle_velocity_positive() {
        let u = particle_velocity(1.0, 413.0);
        assert!(u > 0.0, "particle velocity should be positive");
    }
    #[test]
    fn test_particle_velocity_zero_impedance() {
        let u = particle_velocity(1.0, 0.0);
        assert_eq!(u, 0.0, "zero impedance returns zero");
    }
    #[test]
    fn test_standing_wave_closed_closed_at_node() {
        let k = PI / 2.0;
        let p = standing_wave_pressure_closed_closed(1.0, k, 1.0);
        assert!(p.abs() < 1e-10, "pressure node at λ/4, got {p}");
    }
    #[test]
    fn test_standing_wave_closed_closed_antinode() {
        let p = standing_wave_pressure_closed_closed(1.0, PI / 2.0, 0.0);
        assert!((p - 2.0).abs() < 1e-12, "pressure antinode at x=0, got {p}");
    }
    #[test]
    fn test_standing_wave_open_closed_at_open_end() {
        let p = standing_wave_pressure_open_closed(1.0, PI / 4.0, 0.0);
        assert!(p.abs() < 1e-12, "pressure node at open end, got {p}");
    }
    #[test]
    fn test_closed_closed_resonance_fundamental() {
        let f = closed_closed_resonance(0.5, 343.0, 1);
        assert!(
            (f - 343.0).abs() < 1e-6,
            "fundamental f_1 = 343 Hz, got {f}"
        );
    }
    #[test]
    fn test_closed_closed_resonance_harmonics() {
        let f1 = closed_closed_resonance(0.5, 343.0, 1);
        let f2 = closed_closed_resonance(0.5, 343.0, 2);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-10,
            "second harmonic = 2 × fundamental"
        );
    }
    #[test]
    fn test_resonator_q_factor() {
        let q = resonator_q_factor(1000.0, 10.0);
        assert!((q - 100.0).abs() < 1e-10, "Q = f0/Δf = 100, got {q}");
    }
    #[test]
    fn test_resonance_decay_time() {
        let tau = resonance_decay_time(100.0, 1000.0);
        let expected = 100.0 / (PI * 1000.0);
        assert!((tau - expected).abs() < 1e-10, "τ = Q/(π·f0), got {tau}");
    }
    #[test]
    fn test_acoustic_energy_density_positive() {
        let w = acoustic_energy_density(1.0, 1.2, 343.0);
        assert!(w > 0.0, "energy density should be positive, got {w}");
    }
    #[test]
    fn test_acoustic_energy_density_scales_with_pressure_squared() {
        let w1 = acoustic_energy_density(1.0, 1.2, 343.0);
        let w2 = acoustic_energy_density(2.0, 1.2, 343.0);
        assert!((w2 / w1 - 4.0).abs() < 1e-10, "w ∝ p², ratio = {}", w2 / w1);
    }
    #[test]
    fn test_ula_array_factor_steering_direction_maximum() {
        let theta_steer = 0.0;
        let af = ula_array_factor(
            8,
            0.5 * 343.0 / 1000.0,
            1000.0,
            343.0,
            theta_steer,
            theta_steer,
        );
        assert!(
            (af - 1.0).abs() < 1e-10,
            "AF should be 1 at steering direction, got {af}"
        );
    }
    #[test]
    fn test_ula_array_factor_off_steering_less_than_one() {
        let theta_steer = 0.0;
        let theta_off = 0.3;
        let af = ula_array_factor(
            8,
            0.3 * 343.0 / 1000.0,
            1000.0,
            343.0,
            theta_off,
            theta_steer,
        );
        assert!(af <= 1.0 + 1e-10, "AF ≤ 1 off-axis, got {af}");
    }
    #[test]
    fn test_tube_impedance_open_end_at_resonance() {
        let l = 0.25_f64;
        let c = 343.0_f64;
        let f_quarter = c / (4.0 * l);
        let z_in = tube_input_impedance_open_end_imag(413.0, f_quarter, c, l);
        assert!(
            z_in.abs() < 1.0,
            "tube impedance ~0 at quarter-wave resonance, got {z_in}"
        );
    }
    #[test]
    fn test_tube_impedance_closed_end_at_resonance() {
        let l = 0.25_f64;
        let c = 343.0_f64;
        let f_q = c / (4.0 * l);
        let z_in = tube_input_impedance_closed_end_imag(413.0, f_q, c, l);
        assert!(
            z_in.abs() < 1.0,
            "closed-end tube Z_in ~0 at λ/4 frequency, got {z_in}"
        );
    }
    #[test]
    fn test_power_law_attenuation_positive() {
        let att = power_law_attenuation(0.5, 5.0, 1.0);
        assert!(att > 0.0, "attenuation should be positive, got {att}");
    }
    #[test]
    fn test_power_law_attenuation_frequency_dependence() {
        let att1 = power_law_attenuation(0.5, 2.0, 1.0);
        let att2 = power_law_attenuation(0.5, 4.0, 1.0);
        assert!(
            (att2 / att1 - 2.0).abs() < 1e-10,
            "linear frequency dependence"
        );
    }
    #[test]
    fn test_pulse_echo_attenuation_double_path() {
        let depth = 2.0_f64;
        let expected = 2.0 * power_law_attenuation(0.5, 5.0, 1.0) * depth;
        let round_trip = pulse_echo_attenuation_db(0.5, 5.0, 1.0, depth);
        assert!(
            (round_trip - expected).abs() < 1e-10,
            "round-trip formula, got {round_trip}"
        );
    }
    #[test]
    fn test_water_absorption_positive_at_10khz() {
        let wa = WaterAbsorption::seawater();
        let att = wa.absorption_db_per_km(10.0);
        assert!(
            att > 0.0,
            "seawater absorption should be positive, got {att}"
        );
    }
    #[test]
    fn test_water_absorption_increases_with_frequency() {
        let wa = WaterAbsorption::seawater();
        let att_low = wa.absorption_db_per_km(1.0);
        let att_high = wa.absorption_db_per_km(100.0);
        assert!(
            att_high > att_low,
            "absorption increases with frequency: {att_low} vs {att_high}"
        );
    }
    #[test]
    fn test_water_absorption_freshwater_less_than_seawater_at_10khz() {
        let fw = WaterAbsorption::freshwater();
        let sw = WaterAbsorption::seawater();
        let att_fw = fw.absorption_db_per_km(10.0);
        let att_sw = sw.absorption_db_per_km(10.0);
        assert!(
            att_fw < att_sw,
            "freshwater absorption < seawater at 10 kHz: {att_fw} vs {att_sw}"
        );
    }
    #[test]
    fn test_piston_radiation_resistance_low_ka() {
        let c = 343.0;
        let r = 0.001;
        let f = 100.0;
        let ka = 2.0 * PI * f * r / c;
        let r_normalised = piston_radiation_resistance_normalised(f, r, c);
        let expected = ka * ka / 2.0;
        assert!(
            (r_normalised - expected).abs() < 1e-6,
            "Low ka approximation: {r_normalised} vs {expected}"
        );
    }
    #[test]
    fn test_piston_radiation_resistance_high_ka_approaches_unity() {
        let r_n = piston_radiation_resistance_normalised(1e6, 1.0, 343.0);
        assert!((r_n - 1.0).abs() < 0.01, "High-ka piston R → 1, got {r_n}");
    }
    #[test]
    fn test_piston_directivity_index_positive_for_large_ka() {
        let di = piston_directivity_index_db(100_000.0, 0.01, 343.0);
        assert!(di > 0.0, "DI should be positive for ka >> 1, got {di}");
    }
    #[test]
    fn test_rubber_density_in_range() {
        let r = AcousticMaterial::rubber();
        assert!(
            r.density > 1000.0 && r.density < 1500.0,
            "rubber density ~1200 kg/m³, got {}",
            r.density
        );
    }
    #[test]
    fn test_aluminium_longitudinal_velocity() {
        let al = AcousticMaterial::aluminium();
        let c = al.longitudinal_velocity();
        assert!(c > 5500.0 && c < 7000.0, "Al c_L ~6300 m/s, got {c}");
    }
    #[test]
    fn test_glass_acoustic_impedance_large() {
        let gl = AcousticMaterial::glass();
        let z = gl.acoustic_impedance();
        assert!(z > 1e7, "glass impedance >> air, got {z}");
    }
    #[test]
    fn test_sound_speed_method_equals_longitudinal_velocity() {
        let m = AcousticMaterial::water();
        assert!((m.sound_speed() - m.longitudinal_velocity()).abs() < 1e-10);
    }
    #[test]
    fn test_relative_impedance_water_much_greater_than_one() {
        let m = AcousticMaterial::water();
        let rel = m.relative_impedance();
        assert!(rel > 1000.0, "Z_water / Z_air >> 1, got {rel}");
    }
}
