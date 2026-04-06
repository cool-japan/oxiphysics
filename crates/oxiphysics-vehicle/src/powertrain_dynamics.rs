// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Powertrain dynamics and driveline simulation.
//!
//! Provides combustion engine torque/power/BSFC models, transmission
//! gear calculations, cycle efficiency formulae, fuel injection maps,
//! and clutch/gear-shift dynamics.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_vehicle::powertrain_dynamics::{CombustionEngine, otto_cycle_efficiency};
//!
//! let engine = CombustionEngine::new(2.0, 150.0, 300.0, 800.0, 6500.0);
//! let torque = engine.torque_at_rpm(3000.0);
//! assert!(torque > 0.0);
//!
//! let eta = otto_cycle_efficiency(10.0, 1.4);
//! assert!(eta > 0.0 && eta < 1.0);
//! ```

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ── helpers ───────────────────────────────────────────────────────────────────

/// Clamp `v` to `[lo, hi]`.
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// Bilinear interpolation in a 2-D table.
///
/// `xs` and `ys` must be sorted ascending.
/// `table[i][j]` corresponds to `(xs[i], ys[j])`.
fn bilinear_interp(xs: &[f64], ys: &[f64], table: &[Vec<f64>], x: f64, y: f64) -> f64 {
    if xs.is_empty() || ys.is_empty() || table.is_empty() {
        return 0.0;
    }
    // Find row index
    let ix = {
        let mut i = 0;
        while i + 1 < xs.len() && x > xs[i + 1] {
            i += 1;
        }
        i.min(xs.len().saturating_sub(2))
    };
    // Find column index
    let iy = {
        let mut j = 0;
        while j + 1 < ys.len() && y > ys[j + 1] {
            j += 1;
        }
        j.min(ys.len().saturating_sub(2))
    };
    let dx = if (xs[ix + 1] - xs[ix]).abs() > 1e-30 {
        (x - xs[ix]) / (xs[ix + 1] - xs[ix])
    } else {
        0.0
    };
    let dy = if (ys[iy + 1] - ys[iy]).abs() > 1e-30 {
        (y - ys[iy]) / (ys[iy + 1] - ys[iy])
    } else {
        0.0
    };
    let dx = clamp(dx, 0.0, 1.0);
    let dy = clamp(dy, 0.0, 1.0);
    let r0 =
        table[ix][iy] * (1.0 - dy) + table[ix].get(iy + 1).copied().unwrap_or(table[ix][iy]) * dy;
    let r1 = table[ix + 1].get(iy).copied().unwrap_or(table[ix][iy]) * (1.0 - dy)
        + table[ix + 1].get(iy + 1).copied().unwrap_or(table[ix][iy]) * dy;
    r0 * (1.0 - dx) + r1 * dx
}

// ── CombustionEngine ──────────────────────────────────────────────────────────

/// Simplified internal combustion engine model.
#[derive(Debug, Clone)]
pub struct CombustionEngine {
    /// Engine displacement \[litres\].
    pub displacement_l: f64,
    /// Maximum power \[kW\].
    pub max_power_kw: f64,
    /// Maximum torque \[N·m\].
    pub max_torque_nm: f64,
    /// Idle speed \[RPM\].
    pub idle_rpm: f64,
    /// Redline speed \[RPM\].
    pub redline_rpm: f64,
}

impl CombustionEngine {
    /// Create a new combustion engine.
    ///
    /// * `disp`   – displacement \[litres\].
    /// * `max_p`  – peak power \[kW\].
    /// * `max_t`  – peak torque \[N·m\].
    /// * `idle`   – idle RPM.
    /// * `red`    – redline RPM.
    pub fn new(disp: f64, max_p: f64, max_t: f64, idle: f64, red: f64) -> Self {
        Self {
            displacement_l: disp,
            max_power_kw: max_p,
            max_torque_nm: max_t,
            idle_rpm: idle,
            redline_rpm: red,
        }
    }

    /// Torque at `rpm` \[N·m\].
    ///
    /// Uses a parabolic curve peaking at mid-range: returns 0 at idle/redline,
    /// max torque at the midpoint of the operating range.
    pub fn torque_at_rpm(&self, rpm: f64) -> f64 {
        if rpm <= self.idle_rpm || rpm >= self.redline_rpm {
            return 0.0;
        }
        // Normalised position in operating range [0, 1]
        let x = (rpm - self.idle_rpm) / (self.redline_rpm - self.idle_rpm);
        // Parabolic: peak at x = 0.5
        let torque = self.max_torque_nm * 4.0 * x * (1.0 - x);
        torque.max(0.0)
    }

    /// Power at `rpm` \[kW\].
    ///
    /// `P = T * omega`, capped at `max_power_kw`.
    pub fn power_at_rpm(&self, rpm: f64) -> f64 {
        let torque = self.torque_at_rpm(rpm);
        let omega = rpm * 2.0 * PI / 60.0;
        let power_kw = torque * omega / 1000.0;
        power_kw.min(self.max_power_kw)
    }

    /// Brake-specific fuel consumption at `rpm` and `load` \[g/kWh\].
    ///
    /// `load` is a dimensionless load factor in \[0, 1\].
    /// Returns a simplified BSFC model: lowest near peak efficiency zone.
    pub fn bsfc_at_rpm(&self, rpm: f64, load: f64) -> f64 {
        // Minimum BSFC at mid-range RPM and full load
        let rpm_norm = (rpm - self.idle_rpm) / (self.redline_rpm - self.idle_rpm);
        let rpm_norm = clamp(rpm_norm, 0.0, 1.0);
        // Efficiency degrades at low/high RPM and low load
        let rpm_factor = 1.0 + 0.5 * (rpm_norm - 0.5).powi(2) * 4.0;
        let load_factor = if load > 0.05 {
            1.0 + 0.3 * (1.0 - load)
        } else {
            3.0 // Very high BSFC at near-idle
        };
        230.0 * rpm_factor * load_factor // 230 g/kWh baseline
    }
}

// ── Cycle efficiency formulae ─────────────────────────────────────────────────

/// Compute Otto cycle thermal efficiency.
///
/// `eta = 1 - 1 / r^(gamma - 1)`
///
/// * `compression_ratio` – volumetric compression ratio.
/// * `gamma`             – specific heat ratio (≈ 1.4 for air).
pub fn otto_cycle_efficiency(compression_ratio: f64, gamma: f64) -> f64 {
    if compression_ratio <= 1.0 {
        return 0.0;
    }
    1.0 - 1.0 / compression_ratio.powf(gamma - 1.0)
}

/// Compute Diesel cycle thermal efficiency.
///
/// `eta = 1 - (rc^gamma - 1) / (r^(gamma-1) * gamma * (rc - 1))`
///
/// * `r`     – compression ratio.
/// * `rc`    – cutoff ratio (V_after_combustion / V_TDC).
/// * `gamma` – specific heat ratio.
pub fn diesel_cycle_efficiency(r: f64, rc: f64, gamma: f64) -> f64 {
    if r <= 1.0 || rc <= 1.0 {
        return 0.0;
    }
    let numerator = rc.powf(gamma) - 1.0;
    let denominator = r.powf(gamma - 1.0) * gamma * (rc - 1.0);
    if denominator < 1e-30 {
        return 0.0;
    }
    1.0 - numerator / denominator
}

// ── Transmission ──────────────────────────────────────────────────────────────

/// Multi-speed gearbox model.
#[derive(Debug, Clone)]
pub struct Transmission {
    /// Individual gear ratios (index 0 = 1st gear).
    pub gear_ratios: Vec<f64>,
    /// Final drive ratio.
    pub final_drive: f64,
    /// Transmission mechanical efficiency (0, 1].
    pub efficiency: f64,
}

impl Transmission {
    /// Create a new transmission.
    ///
    /// * `ratios`      – gear ratios, index 0 = 1st gear.
    /// * `final_drive` – final drive ratio.
    /// * `eta`         – efficiency (typically 0.95-0.98).
    pub fn new(ratios: Vec<f64>, final_drive: f64, eta: f64) -> Self {
        Self {
            gear_ratios: ratios,
            final_drive,
            efficiency: eta,
        }
    }

    /// Compute wheel torque \[N·m\] for given engine torque and gear.
    ///
    /// `T_wheel = engine_torque * gear_ratio * final_drive * efficiency`
    ///
    /// Returns 0 if `gear` is out of range.
    pub fn wheel_torque(&self, engine_torque: f64, gear: usize) -> f64 {
        if gear == 0 || gear > self.gear_ratios.len() {
            return 0.0;
        }
        engine_torque * self.gear_ratios[gear - 1] * self.final_drive * self.efficiency
    }

    /// Compute vehicle speed \[m/s\] from engine RPM and wheel radius.
    ///
    /// * `engine_rpm` – engine speed \[RPM\].
    /// * `gear`       – gear number (1-based).
    /// * `wheel_r`    – wheel rolling radius \[m\].
    pub fn wheel_speed(&self, engine_rpm: f64, gear: usize, wheel_r: f64) -> f64 {
        if gear == 0 || gear > self.gear_ratios.len() || wheel_r <= 0.0 {
            return 0.0;
        }
        let total_ratio = self.gear_ratios[gear - 1] * self.final_drive;
        if total_ratio < 1e-10 {
            return 0.0;
        }
        let wheel_rpm = engine_rpm / total_ratio;
        wheel_rpm * 2.0 * PI * wheel_r / 60.0
    }

    /// Return the number of gears.
    pub fn gear_count(&self) -> usize {
        self.gear_ratios.len()
    }
}

// ── clutch_torque_capacity ────────────────────────────────────────────────────

/// Compute dry-clutch torque capacity \[N·m\].
///
/// `T = mu_c * F_normal * r_mean`
///
/// * `mu_c`     – clutch friction coefficient.
/// * `f_normal` – normal (clamping) force \[N\].
/// * `r_mean`   – mean friction radius \[m\].
pub fn clutch_torque_capacity(mu_c: f64, f_normal: f64, r_mean: f64) -> f64 {
    mu_c * f_normal * r_mean
}

// ── gear_shift_time_loss ──────────────────────────────────────────────────────

/// Estimate acceleration lost during a gear change \[m/s²·s = Δv loss, m/s\].
///
/// During `shift_time`, the drivetrain is disconnected so no tractive force.
/// Velocity loss = `engine_torque / (mass * speed) * shift_time` (simplified).
///
/// * `shift_time`     – duration of the gear change \[s\].
/// * `engine_torque`  – torque at the wheels just before shift \[N·m\].
/// * `mass`           – vehicle mass \[kg\].
/// * `speed`          – vehicle speed at time of shift \[m/s\].
pub fn gear_shift_time_loss(shift_time: f64, engine_torque: f64, mass: f64, speed: f64) -> f64 {
    if mass <= 0.0 || speed <= 0.0 {
        return 0.0;
    }
    (engine_torque / (mass * speed)) * shift_time
}

// ── engine_braking_torque ────────────────────────────────────────────────────

/// Estimate engine braking torque \[N·m\] during deceleration.
///
/// Simplified model: pumping + friction losses at zero throttle.
/// `T_brake ≈ displacement * rpm * (1 - throttle) * 0.003`
///
/// * `displacement` – engine displacement \[litres\].
/// * `rpm`          – engine speed \[RPM\].
/// * `throttle`     – throttle position \[0, 1\].
pub fn engine_braking_torque(displacement: f64, rpm: f64, throttle: f64) -> f64 {
    let throttle = clamp(throttle, 0.0, 1.0);
    (displacement * rpm * (1.0 - throttle) * 0.003).abs()
}

// ── fuel_injection_map ────────────────────────────────────────────────────────

/// Look up fuel injection quantity \[mg/stroke\] from a 2-D map.
///
/// Performs bilinear interpolation in the fuel map.
///
/// * `rpm`             – engine speed axis breakpoints \[RPM\].
/// * `load`            – engine load axis breakpoints \[0, 1\] or torque \[N·m\].
/// * `injection_table` – 2-D fuel map, indexed `[rpm_idx][load_idx]`.
pub fn fuel_injection_map(rpm: f64, load: f64, injection_table: &[Vec<f64>]) -> f64 {
    let n_rpm = injection_table.len();
    if n_rpm == 0 {
        return 0.0;
    }
    let n_load = injection_table[0].len();
    if n_load == 0 {
        return 0.0;
    }
    // Default axis: rpm [1000, 2000, ..., 7000], load [0.0, 0.5, 1.0]
    let rpm_axis: Vec<f64> = (0..n_rpm)
        .map(|i| 1000.0 + i as f64 * (6000.0 / (n_rpm as f64 - 1.0).max(1.0)))
        .collect();
    let load_axis: Vec<f64> = (0..n_load)
        .map(|j| j as f64 / (n_load as f64 - 1.0).max(1.0))
        .collect();
    bilinear_interp(&rpm_axis, &load_axis, injection_table, rpm, load)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_engine() -> CombustionEngine {
        CombustionEngine::new(2.0, 150.0, 300.0, 800.0, 6500.0)
    }

    fn default_transmission() -> Transmission {
        Transmission::new(vec![3.5, 2.1, 1.4, 1.0, 0.8], 3.7, 0.97)
    }

    // ── otto_cycle_efficiency ─────────────────────────────────────────────────

    #[test]
    fn otto_efficiency_formula() {
        let r = 10.0_f64;
        let gamma = 1.4_f64;
        let expected = 1.0 - 1.0 / r.powf(gamma - 1.0);
        assert!((otto_cycle_efficiency(r, gamma) - expected).abs() < 1e-12);
    }

    #[test]
    fn otto_efficiency_in_range() {
        let eta = otto_cycle_efficiency(10.0, 1.4);
        assert!(eta > 0.0 && eta < 1.0, "eta={eta}");
    }

    #[test]
    fn otto_efficiency_cr_1_is_zero() {
        assert_eq!(otto_cycle_efficiency(1.0, 1.4), 0.0);
    }

    #[test]
    fn otto_efficiency_increases_with_cr() {
        let e1 = otto_cycle_efficiency(8.0, 1.4);
        let e2 = otto_cycle_efficiency(12.0, 1.4);
        assert!(e2 > e1);
    }

    #[test]
    fn otto_efficiency_never_exceeds_one() {
        for cr in [5.0, 10.0, 15.0, 20.0] {
            assert!(otto_cycle_efficiency(cr, 1.4) < 1.0);
        }
    }

    // ── diesel_cycle_efficiency ────────────────────────────────────────────────

    #[test]
    fn diesel_efficiency_in_range() {
        let eta = diesel_cycle_efficiency(18.0, 2.0, 1.4);
        assert!(eta > 0.0 && eta < 1.0, "eta={eta}");
    }

    #[test]
    fn diesel_efficiency_cr_1_is_zero() {
        assert_eq!(diesel_cycle_efficiency(1.0, 2.0, 1.4), 0.0);
    }

    #[test]
    fn diesel_higher_cr_higher_efficiency() {
        let e1 = diesel_cycle_efficiency(15.0, 2.0, 1.4);
        let e2 = diesel_cycle_efficiency(20.0, 2.0, 1.4);
        assert!(e2 > e1, "e1={e1} e2={e2}");
    }

    // ── CombustionEngine ──────────────────────────────────────────────────────

    #[test]
    fn engine_torque_at_idle_is_zero() {
        let e = default_engine();
        assert_eq!(e.torque_at_rpm(800.0), 0.0);
    }

    #[test]
    fn engine_torque_at_redline_is_zero() {
        let e = default_engine();
        assert_eq!(e.torque_at_rpm(6500.0), 0.0);
    }

    #[test]
    fn engine_torque_positive_mid_range() {
        let e = default_engine();
        let t = e.torque_at_rpm(3650.0); // midpoint
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn engine_torque_peak_at_midpoint() {
        let e = default_engine();
        let mid = (e.idle_rpm + e.redline_rpm) / 2.0;
        let t_mid = e.torque_at_rpm(mid);
        assert!((t_mid - e.max_torque_nm).abs() < 1.0, "t_mid={t_mid}");
    }

    #[test]
    fn engine_power_at_rpm_positive() {
        let e = default_engine();
        let p = e.power_at_rpm(4000.0);
        assert!(p > 0.0, "p={p}");
    }

    #[test]
    fn engine_power_continuous_with_torque() {
        // Power should be continuous between nearby RPM values
        let e = default_engine();
        let p1 = e.power_at_rpm(3000.0);
        let p2 = e.power_at_rpm(3001.0);
        assert!((p1 - p2).abs() < 1.0, "p1={p1} p2={p2}");
    }

    #[test]
    fn engine_bsfc_positive() {
        let e = default_engine();
        let bsfc = e.bsfc_at_rpm(3000.0, 0.8);
        assert!(bsfc > 0.0, "bsfc={bsfc}");
    }

    #[test]
    fn engine_bsfc_higher_at_light_load() {
        let e = default_engine();
        let bsfc_full = e.bsfc_at_rpm(3000.0, 1.0);
        let bsfc_part = e.bsfc_at_rpm(3000.0, 0.3);
        assert!(
            bsfc_part > bsfc_full,
            "bsfc_part={bsfc_part} bsfc_full={bsfc_full}"
        );
    }

    // ── Transmission ──────────────────────────────────────────────────────────

    #[test]
    fn transmission_gear_count_correct() {
        let t = default_transmission();
        assert_eq!(t.gear_count(), 5);
    }

    #[test]
    fn transmission_wheel_torque_first_gear() {
        let t = default_transmission();
        // gear 1: ratio=3.5, final=3.7, eta=0.97
        let tw = t.wheel_torque(200.0, 1);
        let expected = 200.0 * 3.5 * 3.7 * 0.97;
        assert!((tw - expected).abs() < 1e-6 * expected, "tw={tw}");
    }

    #[test]
    fn transmission_wheel_torque_scales_with_engine_torque() {
        let t = default_transmission();
        let tw1 = t.wheel_torque(100.0, 3);
        let tw2 = t.wheel_torque(200.0, 3);
        assert!((tw2 - 2.0 * tw1).abs() < 1e-10, "tw1={tw1} tw2={tw2}");
    }

    #[test]
    fn transmission_wheel_torque_invalid_gear_zero() {
        let t = default_transmission();
        assert_eq!(t.wheel_torque(200.0, 0), 0.0);
        assert_eq!(t.wheel_torque(200.0, 99), 0.0);
    }

    #[test]
    fn transmission_wheel_speed_formula() {
        let t = default_transmission();
        let gear = 3;
        let total_ratio = t.gear_ratios[gear - 1] * t.final_drive;
        let wheel_r = 0.33;
        let rpm = 3000.0;
        let expected = (rpm / total_ratio) * 2.0 * PI * wheel_r / 60.0;
        let speed = t.wheel_speed(rpm, gear, wheel_r);
        assert!((speed - expected).abs() < 1e-10, "speed={speed}");
    }

    #[test]
    fn transmission_wheel_speed_positive() {
        let t = default_transmission();
        let speed = t.wheel_speed(3000.0, 3, 0.33);
        assert!(speed > 0.0, "speed={speed}");
    }

    #[test]
    fn transmission_higher_gear_higher_speed() {
        let t = default_transmission();
        let s1 = t.wheel_speed(3000.0, 1, 0.33);
        let s4 = t.wheel_speed(3000.0, 4, 0.33);
        assert!(s4 > s1, "s1={s1} s4={s4}");
    }

    // ── clutch_torque_capacity ─────────────────────────────────────────────────

    #[test]
    fn clutch_torque_formula() {
        let t = clutch_torque_capacity(0.35, 8000.0, 0.10);
        let expected = 0.35 * 8000.0 * 0.10;
        assert!((t - expected).abs() < 1e-10, "t={t}");
    }

    #[test]
    fn clutch_torque_positive() {
        assert!(clutch_torque_capacity(0.35, 8000.0, 0.10) > 0.0);
    }

    #[test]
    fn clutch_torque_scales_linearly() {
        let t1 = clutch_torque_capacity(0.3, 5000.0, 0.08);
        let t2 = clutch_torque_capacity(0.3, 10000.0, 0.08);
        assert!((t2 - 2.0 * t1).abs() < 1e-10);
    }

    // ── engine_braking_torque ─────────────────────────────────────────────────

    #[test]
    fn engine_braking_positive_at_zero_throttle() {
        let t = engine_braking_torque(2.0, 3000.0, 0.0);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn engine_braking_zero_at_full_throttle() {
        let t = engine_braking_torque(2.0, 3000.0, 1.0);
        assert!(t.abs() < 1e-10, "t={t}");
    }

    #[test]
    fn engine_braking_increases_with_rpm() {
        let t1 = engine_braking_torque(2.0, 2000.0, 0.0);
        let t2 = engine_braking_torque(2.0, 4000.0, 0.0);
        assert!(t2 > t1, "t1={t1} t2={t2}");
    }

    // ── fuel_injection_map ─────────────────────────────────────────────────────

    #[test]
    fn fuel_injection_map_positive() {
        let table: Vec<Vec<f64>> = vec![
            vec![10.0, 20.0, 30.0],
            vec![12.0, 22.0, 32.0],
            vec![14.0, 24.0, 34.0],
        ];
        let qty = fuel_injection_map(2000.0, 0.5, &table);
        assert!(qty > 0.0, "qty={qty}");
    }

    #[test]
    fn fuel_injection_map_empty_is_zero() {
        assert_eq!(fuel_injection_map(2000.0, 0.5, &[]), 0.0);
    }

    #[test]
    fn fuel_injection_map_higher_load_higher_fuel() {
        let table: Vec<Vec<f64>> = vec![
            vec![10.0, 20.0, 30.0],
            vec![12.0, 22.0, 32.0],
            vec![14.0, 24.0, 34.0],
        ];
        let q1 = fuel_injection_map(2000.0, 0.0, &table);
        let q2 = fuel_injection_map(2000.0, 1.0, &table);
        assert!(q2 > q1, "q1={q1} q2={q2}");
    }

    // ── gear_shift_time_loss ───────────────────────────────────────────────────

    #[test]
    fn gear_shift_time_loss_positive() {
        let loss = gear_shift_time_loss(0.25, 2000.0, 1500.0, 20.0);
        assert!(loss > 0.0, "loss={loss}");
    }

    #[test]
    fn gear_shift_zero_mass_is_zero() {
        assert_eq!(gear_shift_time_loss(0.25, 2000.0, 0.0, 20.0), 0.0);
    }

    #[test]
    fn gear_shift_zero_speed_is_zero() {
        assert_eq!(gear_shift_time_loss(0.25, 2000.0, 1500.0, 0.0), 0.0);
    }
}
