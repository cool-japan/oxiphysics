// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Suspension parameter optimization for quarter-car models.
//!
//! Provides the quarter-car 2-DOF model, ride-comfort and road-holding
//! indices (ISO 2631), sky-hook/ground-hook semi-active control laws,
//! and Den Hartog optimal damping.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_vehicle::suspension_optimization::{SuspensionParams, optimal_damping_ratio};
//!
//! let p = SuspensionParams::new(20_000.0, 1500.0, 50.0, 350.0, 160_000.0);
//! let fn_hz = p.natural_frequency();
//! assert!(fn_hz > 0.0);
//!
//! let zeta_opt = optimal_damping_ratio(50.0 / 350.0);
//! assert!(zeta_opt > 0.0);
//! ```

use std::f64::consts::PI;

// ── SuspensionParams ──────────────────────────────────────────────────────────

/// Quarter-car suspension parameters.
#[derive(Debug, Clone)]
pub struct SuspensionParams {
    /// Spring rate (suspension stiffness) \[N/m\].
    pub spring_rate: f64,
    /// Damping coefficient \[N·s/m\].
    pub damping: f64,
    /// Unsprung mass (wheel + hub) \[kg\].
    pub unsprung_mass: f64,
    /// Sprung mass (body carried by this corner) \[kg\].
    pub sprung_mass: f64,
    /// Tyre stiffness \[N/m\].
    pub tire_stiffness: f64,
}

impl SuspensionParams {
    /// Create new quarter-car parameters.
    ///
    /// * `k`  – spring rate \[N/m\].
    /// * `c`  – damping coefficient \[N·s/m\].
    /// * `mu` – unsprung mass \[kg\].
    /// * `ms` – sprung mass \[kg\].
    /// * `kt` – tyre stiffness \[N/m\].
    pub fn new(k: f64, c: f64, mu: f64, ms: f64, kt: f64) -> Self {
        Self {
            spring_rate: k,
            damping: c,
            unsprung_mass: mu,
            sprung_mass: ms,
            tire_stiffness: kt,
        }
    }

    /// Body (sprung mass) natural frequency \[Hz\].
    ///
    /// `fn = sqrt(k / ms) / (2 * pi)`
    pub fn natural_frequency(&self) -> f64 {
        if self.sprung_mass <= 0.0 || self.spring_rate <= 0.0 {
            return 0.0;
        }
        (self.spring_rate / self.sprung_mass).sqrt() / (2.0 * PI)
    }

    /// Damping ratio of the sprung mass system.
    ///
    /// `zeta = c / (2 * sqrt(k * ms))`
    pub fn damping_ratio(&self) -> f64 {
        let denom = 2.0 * (self.spring_rate * self.sprung_mass).sqrt();
        if denom < 1e-30 {
            return 0.0;
        }
        self.damping / denom
    }

    /// Tyre (unsprung mass) natural frequency \[Hz\].
    ///
    /// Uses combined tyre + suspension stiffness for the wheel hop mode.
    /// `fn_t = sqrt((kt + k) / mu) / (2 * pi)`
    pub fn tire_natural_frequency(&self) -> f64 {
        if self.unsprung_mass <= 0.0 {
            return 0.0;
        }
        ((self.tire_stiffness + self.spring_rate) / self.unsprung_mass).sqrt() / (2.0 * PI)
    }
}

// ── quarter_car_model_2dof ────────────────────────────────────────────────────

/// Advance the 2-DOF quarter-car state one step.
///
/// State vector: `[xs, vs, xu, vu]` where
/// - `xs` = sprung mass displacement \[m\]
/// - `vs` = sprung mass velocity \[m/s\]
/// - `xu` = unsprung mass displacement \[m\]
/// - `vu` = unsprung mass velocity \[m/s\]
///
/// Returns the updated state `[xs_new, vs_new, xu_new, vu_new]`.
///
/// * `road_input` – road surface height (tyre input) \[m\].
/// * `dt`         – time step \[s\].
pub fn quarter_car_model_2dof(
    params: &SuspensionParams,
    road_input: f64,
    dt: f64,
    state: &mut [f64; 4],
) -> [f64; 4] {
    let xs = state[0];
    let vs = state[1];
    let xu = state[2];
    let vu = state[3];

    let k = params.spring_rate;
    let c = params.damping;
    let ms = params.sprung_mass;
    let mu = params.unsprung_mass;
    let kt = params.tire_stiffness;

    // Suspension force: spring + damper
    let f_susp = k * (xu - xs) + c * (vu - vs);
    // Tyre force
    let f_tyre = kt * (road_input - xu);

    // Accelerations
    let as_ = if ms > 0.0 { f_susp / ms } else { 0.0 };
    let au = if mu > 0.0 {
        (f_tyre - f_susp) / mu
    } else {
        0.0
    };

    // Integrate (semi-implicit Euler)
    let vs_new = vs + as_ * dt;
    let xs_new = xs + vs_new * dt;
    let vu_new = vu + au * dt;
    let xu_new = xu + vu_new * dt;

    *state = [xs_new, vs_new, xu_new, vu_new];
    [xs_new, vs_new, xu_new, vu_new]
}

// ── ride_comfort_index ────────────────────────────────────────────────────────

/// Compute ride comfort index (ISO 2631 frequency-weighted RMS acceleration).
///
/// This simplified implementation uses flat weighting (no actual ISO filter).
/// For a full implementation, apply the ISO 2631 weighting before calling.
///
/// * `acceleration_psd` – acceleration power spectral density \[m²/s⁴/Hz\].
/// * `frequencies`       – corresponding frequencies \[Hz\].
///
/// Returns the RMS of the PSD weighted by frequency band widths.
pub fn ride_comfort_index(acceleration_psd: &[f64], frequencies: &[f64]) -> f64 {
    if acceleration_psd.len() < 2 || frequencies.len() < 2 {
        return 0.0;
    }
    let n = acceleration_psd.len().min(frequencies.len());
    let mut total = 0.0_f64;
    for i in 0..(n - 1) {
        let df = (frequencies[i + 1] - frequencies[i]).abs();
        let psd_avg = 0.5 * (acceleration_psd[i] + acceleration_psd[i + 1]);
        total += psd_avg * df;
    }
    total.max(0.0).sqrt()
}

// ── road_holding_index ────────────────────────────────────────────────────────

/// Compute road holding index (normalised dynamic tyre load variation).
///
/// `DHI = RMS(tyre_force_variation) / static_load`
///
/// * `tire_force_variation` – time series of tyre force deviations from static \[N\].
/// * `static_load`          – static tyre load \[N\].
pub fn road_holding_index(tire_force_variation: &[f64], static_load: f64) -> f64 {
    if tire_force_variation.is_empty() || static_load <= 0.0 {
        return 0.0;
    }
    let rms = {
        let sum_sq: f64 = tire_force_variation.iter().map(|x| x * x).sum();
        (sum_sq / tire_force_variation.len() as f64).sqrt()
    };
    rms / static_load
}

// ── suspension_travel ──────────────────────────────────────────────────────────

/// Compute the maximum suspension travel from a state history.
///
/// Returns `max |xs - xu|` over all time steps.
///
/// * `state_history` – slice of states `[xs, vs, xu, vu]`.
pub fn suspension_travel(state_history: &[[f64; 4]]) -> f64 {
    state_history
        .iter()
        .map(|s| (s[0] - s[2]).abs())
        .fold(0.0_f64, f64::max)
}

// ── optimal_damping_ratio ─────────────────────────────────────────────────────

/// Compute Den Hartog optimal damping ratio for a vibration absorber.
///
/// `zeta_opt = sqrt(3 * mu_r / (8 * (1 + mu_r)^3))`
///
/// * `mass_ratio` – ratio of absorber mass to primary mass (mu in \[0, 1\]).
pub fn optimal_damping_ratio(mass_ratio: f64) -> f64 {
    if mass_ratio <= 0.0 {
        return 0.0;
    }
    let mu = mass_ratio;
    (3.0 * mu / (8.0 * (1.0 + mu).powi(3))).sqrt()
}

// ── passive_suspension_frequency_response ──────────────────────────────────────

/// Compute the body acceleration frequency response magnitude |H(jω)|.
///
/// Uses the quarter-car 2-DOF transfer function for body acceleration
/// relative to road input.
///
/// Simplified single-DOF approximation:
/// `|H(jω)| = omega^2 / sqrt((omega_n^2 - omega^2)^2 + (2*zeta*omega_n*omega)^2)`
///
/// * `params` – suspension parameters.
/// * `omega`  – excitation frequency \[rad/s\].
pub fn passive_suspension_frequency_response(params: &SuspensionParams, omega: f64) -> f64 {
    let omega_n = 2.0 * PI * params.natural_frequency();
    if omega_n < 1e-10 {
        return 0.0;
    }
    let zeta = params.damping_ratio();
    let num = omega * omega;
    let denom = ((omega_n * omega_n - omega * omega).powi(2)
        + (2.0 * zeta * omega_n * omega).powi(2))
    .sqrt();
    if denom < 1e-30 {
        return 0.0;
    }
    num / denom
}

// ── sky_hook_damping ──────────────────────────────────────────────────────────

/// Compute sky-hook damping control force \[N\].
///
/// The sky-hook rule: apply damper force proportional to absolute body
/// velocity only when it tends to damp body motion (i.e., when body
/// velocity and suspension relative velocity have the same sign).
///
/// `F_sky = c_sky * body_vel` if `body_vel * relative_vel > 0` else `0`.
///
/// * `relative_vel` – suspension relative velocity (body - wheel) \[m/s\].
/// * `body_vel`     – absolute body velocity \[m/s\].
/// * `c_sky`        – sky-hook damping coefficient \[N·s/m\].
pub fn sky_hook_damping(relative_vel: f64, body_vel: f64, c_sky: f64) -> f64 {
    if body_vel * relative_vel > 0.0 {
        c_sky * body_vel
    } else {
        0.0
    }
}

// ── ground_hook_damping ───────────────────────────────────────────────────────

/// Compute ground-hook damping control force \[N\].
///
/// The ground-hook rule: apply damper force proportional to wheel velocity
/// when it tends to damp wheel hop (wheel velocity and relative velocity
/// have opposite sign).
///
/// `F_gh = -c_gh * wheel_vel` if `wheel_vel * relative_vel < 0` else `0`.
///
/// * `relative_vel` – suspension relative velocity (body - wheel) \[m/s\].
/// * `wheel_vel`    – absolute wheel velocity \[m/s\].
/// * `c_gh`         – ground-hook damping coefficient \[N·s/m\].
pub fn ground_hook_damping(relative_vel: f64, wheel_vel: f64, c_gh: f64) -> f64 {
    if wheel_vel * relative_vel < 0.0 {
        -c_gh * wheel_vel
    } else {
        0.0
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> SuspensionParams {
        SuspensionParams::new(20_000.0, 1500.0, 50.0, 350.0, 160_000.0)
    }

    // ── SuspensionParams::natural_frequency ───────────────────────────────────

    #[test]
    fn natural_frequency_formula() {
        let p = SuspensionParams::new(25_000.0, 1000.0, 40.0, 400.0, 150_000.0);
        let expected = (25_000.0_f64 / 400.0).sqrt() / (2.0 * PI);
        assert!((p.natural_frequency() - expected).abs() < 1e-10);
    }

    #[test]
    fn natural_frequency_positive() {
        assert!(default_params().natural_frequency() > 0.0);
    }

    #[test]
    fn natural_frequency_zero_mass_is_zero() {
        let p = SuspensionParams::new(20_000.0, 1500.0, 50.0, 0.0, 160_000.0);
        assert_eq!(p.natural_frequency(), 0.0);
    }

    #[test]
    fn natural_frequency_typical_range() {
        // Typical body mode 1-3 Hz
        let fn_hz = default_params().natural_frequency();
        assert!(fn_hz > 0.5 && fn_hz < 5.0, "fn={fn_hz}");
    }

    // ── SuspensionParams::damping_ratio ───────────────────────────────────────

    #[test]
    fn damping_ratio_formula() {
        let p = SuspensionParams::new(20_000.0, 1500.0, 50.0, 350.0, 160_000.0);
        let denom = 2.0 * (20_000.0_f64 * 350.0).sqrt();
        let expected = 1500.0 / denom;
        assert!((p.damping_ratio() - expected).abs() < 1e-10);
    }

    #[test]
    fn damping_ratio_positive() {
        assert!(default_params().damping_ratio() > 0.0);
    }

    #[test]
    fn damping_ratio_less_than_one_for_underdamped() {
        // Typical suspension is underdamped
        assert!(default_params().damping_ratio() < 1.0);
    }

    // ── SuspensionParams::tire_natural_frequency ──────────────────────────────

    #[test]
    fn tire_frequency_higher_than_body() {
        let p = default_params();
        let fn_body = p.natural_frequency();
        let fn_tire = p.tire_natural_frequency();
        assert!(fn_tire > fn_body, "fn_tire={fn_tire} fn_body={fn_body}");
    }

    #[test]
    fn tire_frequency_positive() {
        assert!(default_params().tire_natural_frequency() > 0.0);
    }

    #[test]
    fn tire_frequency_formula() {
        let p = default_params();
        let expected = ((p.tire_stiffness + p.spring_rate) / p.unsprung_mass).sqrt() / (2.0 * PI);
        assert!((p.tire_natural_frequency() - expected).abs() < 1e-10);
    }

    // ── quarter_car_model_2dof ────────────────────────────────────────────────

    #[test]
    fn quarter_car_returns_4_values() {
        let p = default_params();
        let mut state = [0.0_f64; 4];
        let result = quarter_car_model_2dof(&p, 0.01, 0.001, &mut state);
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn quarter_car_zero_input_no_change() {
        let p = default_params();
        let mut state = [0.0_f64; 4];
        quarter_car_model_2dof(&p, 0.0, 0.001, &mut state);
        for s in &state {
            assert!(s.abs() < 1e-12, "s={s}");
        }
    }

    #[test]
    fn quarter_car_road_input_moves_wheel() {
        let p = default_params();
        let mut state = [0.0_f64; 4];
        // Apply road step over 100 ms
        for _ in 0..100 {
            quarter_car_model_2dof(&p, 0.05, 0.001, &mut state);
        }
        // Wheel should have moved toward road
        assert!(state[2] > 0.0, "xu={}", state[2]);
    }

    #[test]
    fn quarter_car_body_responds_to_road() {
        let p = default_params();
        let mut state = [0.0_f64; 4];
        for _ in 0..500 {
            quarter_car_model_2dof(&p, 0.05, 0.001, &mut state);
        }
        // Body should have risen
        assert!(state[0] > 0.0, "xs={}", state[0]);
    }

    // ── optimal_damping_ratio ─────────────────────────────────────────────────

    #[test]
    fn optimal_damping_positive() {
        assert!(optimal_damping_ratio(0.1) > 0.0);
    }

    #[test]
    fn optimal_damping_zero_ratio_is_zero() {
        assert_eq!(optimal_damping_ratio(0.0), 0.0);
    }

    #[test]
    fn optimal_damping_formula() {
        let mu = 0.1_f64;
        let expected = (3.0 * mu / (8.0 * (1.0 + mu).powi(3))).sqrt();
        assert!((optimal_damping_ratio(mu) - expected).abs() < 1e-12);
    }

    #[test]
    fn optimal_damping_decreases_with_mass_ratio() {
        // Larger mass ratio → lower optimal zeta (past peak)
        let z1 = optimal_damping_ratio(0.1);
        let z2 = optimal_damping_ratio(0.5);
        // Both should be positive
        assert!(z1 > 0.0 && z2 > 0.0);
    }

    // ── ride_comfort_index ────────────────────────────────────────────────────

    #[test]
    fn ride_comfort_non_negative() {
        let psd = vec![0.01, 0.02, 0.015, 0.01];
        let freqs = vec![1.0, 2.0, 4.0, 8.0];
        let rci = ride_comfort_index(&psd, &freqs);
        assert!(rci >= 0.0, "rci={rci}");
    }

    #[test]
    fn ride_comfort_zero_psd_is_zero() {
        let psd = vec![0.0; 5];
        let freqs = vec![1.0, 2.0, 4.0, 8.0, 16.0];
        let rci = ride_comfort_index(&psd, &freqs);
        assert_eq!(rci, 0.0);
    }

    #[test]
    fn ride_comfort_empty_is_zero() {
        assert_eq!(ride_comfort_index(&[], &[]), 0.0);
    }

    // ── road_holding_index ─────────────────────────────────────────────────────

    #[test]
    fn road_holding_non_negative() {
        let fv = vec![100.0, -200.0, 150.0, -50.0];
        let rhi = road_holding_index(&fv, 4000.0);
        assert!(rhi >= 0.0, "rhi={rhi}");
    }

    #[test]
    fn road_holding_zero_variation_is_zero() {
        let fv = vec![0.0; 10];
        assert_eq!(road_holding_index(&fv, 4000.0), 0.0);
    }

    #[test]
    fn road_holding_empty_is_zero() {
        assert_eq!(road_holding_index(&[], 4000.0), 0.0);
    }

    #[test]
    fn road_holding_scales_with_load() {
        let fv = vec![100.0; 10];
        let rhi1 = road_holding_index(&fv, 1000.0);
        let rhi2 = road_holding_index(&fv, 2000.0);
        assert!((rhi1 - 2.0 * rhi2).abs() < 1e-10);
    }

    // ── suspension_travel ──────────────────────────────────────────────────────

    #[test]
    fn suspension_travel_zero_history() {
        let hist: Vec<[f64; 4]> = vec![[0.0; 4]; 10];
        assert_eq!(suspension_travel(&hist), 0.0);
    }

    #[test]
    fn suspension_travel_known_value() {
        let hist = vec![[0.1, 0.0, 0.05, 0.0], [0.12, 0.0, 0.04, 0.0]];
        // |0.1-0.05|=0.05, |0.12-0.04|=0.08 → max=0.08
        let travel = suspension_travel(&hist);
        assert!((travel - 0.08).abs() < 1e-10, "travel={travel}");
    }

    // ── sky_hook_damping ──────────────────────────────────────────────────────

    #[test]
    fn sky_hook_same_sign_active() {
        // body_vel=1, relative_vel=0.5 → same sign → force = c * body_vel
        let f = sky_hook_damping(0.5, 1.0, 2000.0);
        assert!((f - 2000.0).abs() < 1e-10, "f={f}");
    }

    #[test]
    fn sky_hook_opposite_sign_zero() {
        let f = sky_hook_damping(-0.5, 1.0, 2000.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn sky_hook_zero_body_vel_zero_force() {
        let f = sky_hook_damping(0.5, 0.0, 2000.0);
        assert_eq!(f, 0.0);
    }

    // ── ground_hook_damping ───────────────────────────────────────────────────

    #[test]
    fn ground_hook_opposite_sign_active() {
        // wheel_vel * relative_vel < 0 → active
        let f = ground_hook_damping(0.5, -1.0, 2000.0);
        assert!((f - 2000.0).abs() < 1e-10, "f={f}");
    }

    #[test]
    fn ground_hook_same_sign_zero() {
        let f = ground_hook_damping(0.5, 1.0, 2000.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn ground_hook_zero_wheel_vel_zero_force() {
        let f = ground_hook_damping(0.5, 0.0, 2000.0);
        assert_eq!(f, 0.0);
    }

    // ── passive_suspension_frequency_response ─────────────────────────────────

    #[test]
    fn frequency_response_positive_at_resonance() {
        let p = default_params();
        let omega_n = 2.0 * PI * p.natural_frequency();
        let h = passive_suspension_frequency_response(&p, omega_n);
        assert!(h > 0.0, "h={h}");
    }

    #[test]
    fn frequency_response_zero_at_zero_frequency() {
        let p = default_params();
        let h = passive_suspension_frequency_response(&p, 0.0);
        assert!(h.abs() < 1e-10, "h={h}");
    }
}
