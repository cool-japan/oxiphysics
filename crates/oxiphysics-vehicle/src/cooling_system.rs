// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle cooling system thermal model.
//!
//! Implements NTU-effectiveness heat exchanger theory, thermostat valve
//! control, lumped-capacitance coolant ODE, and fan curve models.
//!
//! # Overview
//!
//! - [`CoolingLoop`] — coolant loop with flow rate and heat exchangers.
//! - [`radiator_effectiveness`] — NTU method: ε = 1 − exp(−NTU(1+Cr)).
//! - [`heat_rejection_rate`] — Q = ε · C_min · (T_hot − T_cold).
//! - [`thermostat_valve`] — temperature-dependent bypass flow fraction.
//! - [`coolant_temperature_ode`] — lumped capacitance thermal ODE step.
//! - [`fan_curve`] — fan pressure rise vs volumetric flow rate.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// CoolingLoop
// ─────────────────────────────────────────────────────────────────────────────

/// A vehicle cooling loop combining a coolant circuit and one or more heat
/// exchangers (radiator, charge air cooler, etc.).
#[derive(Debug, Clone)]
pub struct CoolingLoop {
    /// Total coolant flow rate through the loop \[kg/s\].
    pub flow_rate: f64,
    /// Specific heat capacity of the coolant \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Thermal capacity (mass × Cp) of the lumped coolant volume \[J/K\].
    pub thermal_capacity: f64,
    /// Ambient temperature \[K\].
    pub ambient_temp: f64,
    /// Current coolant temperature \[K\].
    pub coolant_temp: f64,
    /// Number of transfer units (NTU) of the main radiator.
    pub ntu: f64,
    /// Capacity rate ratio of the radiator Cr = C_min / C_max.
    pub capacity_rate_ratio: f64,
    /// Thermostat opening temperature \[K\].
    pub thermostat_open_temp: f64,
    /// Thermostat fully-open temperature \[K\] (linear ramp from open to full).
    pub thermostat_full_temp: f64,
}

impl CoolingLoop {
    /// Construct a cooling loop with given parameters.
    pub fn new(
        flow_rate: f64,
        specific_heat: f64,
        thermal_capacity: f64,
        ambient_temp: f64,
        coolant_temp: f64,
        ntu: f64,
        capacity_rate_ratio: f64,
        thermostat_open_temp: f64,
        thermostat_full_temp: f64,
    ) -> Self {
        Self {
            flow_rate,
            specific_heat,
            thermal_capacity,
            ambient_temp,
            coolant_temp,
            ntu,
            capacity_rate_ratio,
            thermostat_open_temp,
            thermostat_full_temp,
        }
    }

    /// Construct a typical water-glycol passenger car cooling loop.
    pub fn passenger_car() -> Self {
        Self::new(
            0.5,      // flow_rate [kg/s]
            3_800.0,  // specific_heat [J/(kg·K)] — 50/50 water-glycol
            15_000.0, // thermal_capacity [J/K]
            298.15,   // ambient_temp [K] = 25°C
            368.15,   // coolant_temp [K] = 95°C
            3.0,      // NTU
            0.5,      // Cr
            353.15,   // thermostat opens at 80°C
            368.15,   // fully open at 95°C
        )
    }

    /// Advance coolant temperature by one time step `dt` \[s\] given engine heat
    /// input `q_engine` \[W\].
    ///
    /// Internally calls [`coolant_temperature_ode`].
    pub fn step(&mut self, q_engine: f64, dt: f64) {
        let valve = thermostat_valve(
            self.coolant_temp,
            self.thermostat_open_temp,
            self.thermostat_full_temp,
        );
        let effective_ntu = self.ntu * valve;
        let eps = radiator_effectiveness(effective_ntu, self.capacity_rate_ratio);
        let c_min = self.flow_rate * self.specific_heat;
        let q_rad = heat_rejection_rate(eps, c_min, self.coolant_temp, self.ambient_temp);
        let dtemp = coolant_temperature_ode(self.thermal_capacity, q_engine, q_rad, dt);
        self.coolant_temp += dtemp;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NTU / effectiveness
// ─────────────────────────────────────────────────────────────────────────────

/// Compute radiator (heat exchanger) effectiveness using the NTU method.
///
/// For a single-pass counter-flow heat exchanger:
/// ε = (1 − exp(−NTU·(1+Cr))) / (1 + Cr)
///
/// For the cross-flow (mixed–unmixed) and other configurations a similar
/// formula is used.  This implementation uses the counter-flow formula which
/// is a common simplification for automotive radiators.
///
/// # Arguments
///
/// * `ntu` — number of transfer units (NTU = UA/C_min ≥ 0).
/// * `cr` — capacity rate ratio C_min/C_max ∈ \[0, 1\].
///
/// Returns effectiveness ε ∈ \[0, 1\].
pub fn radiator_effectiveness(ntu: f64, cr: f64) -> f64 {
    if ntu <= 0.0 {
        return 0.0;
    }
    let cr = cr.clamp(0.0, 1.0);
    if cr < 1e-9 {
        // Condensing/evaporating case: ε = 1 − exp(−NTU)
        return 1.0 - (-ntu).exp();
    }
    if (cr - 1.0).abs() < 1e-9 {
        // Balanced flow: ε = NTU / (1 + NTU)
        return ntu / (1.0 + ntu);
    }
    // General counter-flow formula
    let numerator = 1.0 - (-(ntu * (1.0 - cr))).exp();
    let denominator = 1.0 - cr * (-(ntu * (1.0 - cr))).exp();
    (numerator / denominator).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Heat rejection rate
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the heat rejection rate of a heat exchanger.
///
/// Q = ε · C_min · (T_hot − T_cold)
///
/// # Arguments
///
/// * `effectiveness` — ε ∈ \[0, 1\] from [`radiator_effectiveness`].
/// * `c_min` — minimum capacity rate C_min = (ṁ·cp)_min \[W/K\].
/// * `t_hot` — inlet temperature of the hot fluid \[K or °C\].
/// * `t_cold` — inlet temperature of the cold fluid \[K or °C\].
///
/// Returns heat rejection rate \[W\].
pub fn heat_rejection_rate(effectiveness: f64, c_min: f64, t_hot: f64, t_cold: f64) -> f64 {
    effectiveness * c_min * (t_hot - t_cold).max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Thermostat valve
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the thermostat valve open fraction \[0, 1\] as a function of coolant
/// temperature.
///
/// - Below `t_open`: valve is fully closed (flow bypasses radiator).
/// - Between `t_open` and `t_full`: linear ramp from 0 to 1.
/// - Above `t_full`: valve is fully open.
///
/// # Arguments
///
/// * `t_coolant` — current coolant temperature \[K or °C\].
/// * `t_open` — temperature at which the thermostat starts opening \[K or °C\].
/// * `t_full` — temperature at which the thermostat is fully open \[K or °C\].
///
/// Returns valve position ∈ \[0, 1\].
pub fn thermostat_valve(t_coolant: f64, t_open: f64, t_full: f64) -> f64 {
    if t_coolant <= t_open {
        return 0.0;
    }
    if t_coolant >= t_full {
        return 1.0;
    }
    let range = t_full - t_open;
    if range < 1e-9 {
        return 1.0;
    }
    (t_coolant - t_open) / range
}

// ─────────────────────────────────────────────────────────────────────────────
// Coolant temperature ODE (lumped capacitance)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the temperature change of the coolant over one Euler step.
///
/// Lumped capacitance model:
/// dT/dt = (Q_engine − Q_radiator) / C_thermal
///
/// # Arguments
///
/// * `thermal_capacity` — coolant thermal mass C = m·cp \[J/K\].
/// * `q_engine` — heat input from engine \[W\].
/// * `q_radiator` — heat removed by radiator \[W\].
/// * `dt` — time step \[s\].
///
/// Returns ΔT \[K\] to be added to the current temperature.
pub fn coolant_temperature_ode(
    thermal_capacity: f64,
    q_engine: f64,
    q_radiator: f64,
    dt: f64,
) -> f64 {
    if thermal_capacity < 1e-12 {
        return 0.0;
    }
    let dqdt = q_engine - q_radiator;
    dqdt / thermal_capacity * dt
}

// ─────────────────────────────────────────────────────────────────────────────
// Fan curve
// ─────────────────────────────────────────────────────────────────────────────

/// Fan operating point descriptor.
#[derive(Debug, Clone, Copy)]
pub struct FanOperatingPoint {
    /// Volumetric flow rate \[m³/s\].
    pub flow_rate: f64,
    /// Static pressure rise \[Pa\].
    pub pressure_rise: f64,
    /// Power consumed by the fan \[W\].
    pub power: f64,
}

/// Compute the fan static pressure rise and power for a given flow rate.
///
/// Uses a simplified quadratic fan curve:
/// ΔP = ΔP_max · (1 − (Q/Q_max)²)
/// P_fan = ΔP · Q / η
///
/// # Arguments
///
/// * `flow_rate` — volumetric flow rate \[m³/s\] at the desired operating point.
/// * `dp_max` — maximum (stall) pressure rise \[Pa\] at zero flow.
/// * `q_max` — free-delivery flow rate \[m³/s\] at zero pressure.
/// * `efficiency` — fan total-to-static efficiency ∈ (0, 1].
///
/// Returns a [`FanOperatingPoint`].
pub fn fan_curve(flow_rate: f64, dp_max: f64, q_max: f64, efficiency: f64) -> FanOperatingPoint {
    let q = flow_rate.clamp(0.0, q_max);
    let ratio = if q_max > 1e-12 { q / q_max } else { 0.0 };
    let pressure_rise = dp_max * (1.0 - ratio * ratio).max(0.0);
    let eta = efficiency.clamp(1e-6, 1.0);
    let power = pressure_rise * q / eta;
    FanOperatingPoint {
        flow_rate: q,
        pressure_rise,
        power,
    }
}

/// Return the flow rate at which the fan operating point intersects a system
/// curve modelled as ΔP_sys = R · Q².
///
/// Solves: ΔP_max · (1 − (Q/Q_max)²) = R · Q²
/// → Q² · (ΔP_max/Q_max² + R) = ΔP_max
/// → Q = sqrt(ΔP_max / (ΔP_max/Q_max² + R))
///
/// # Arguments
///
/// * `dp_max` — stall pressure \[Pa\].
/// * `q_max` — free-delivery flow \[m³/s\].
/// * `system_resistance` — system curve coefficient R \[Pa/(m³/s)²\].
///
/// Returns equilibrium flow rate \[m³/s\].
pub fn fan_system_operating_point(dp_max: f64, q_max: f64, system_resistance: f64) -> f64 {
    let denom = dp_max / (q_max * q_max).max(1e-18) + system_resistance;
    if denom < 1e-12 {
        return 0.0;
    }
    (dp_max / denom).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── CoolingLoop ──────────────────────────────────────────────────────

    #[test]
    fn passenger_car_loop_construction() {
        let loop_ = CoolingLoop::passenger_car();
        assert!(loop_.coolant_temp > loop_.ambient_temp);
        assert!(loop_.ntu > 0.0);
    }

    #[test]
    fn cooling_loop_step_reduces_temp_when_no_engine_heat() {
        let mut loop_ = CoolingLoop::passenger_car();
        let initial_temp = loop_.coolant_temp;
        // No engine heat → radiator should cool the coolant
        loop_.step(0.0, 1.0);
        assert!(
            loop_.coolant_temp <= initial_temp,
            "temp should not rise with no engine heat: {} vs {}",
            loop_.coolant_temp,
            initial_temp
        );
    }

    #[test]
    fn cooling_loop_step_raises_temp_with_high_engine_heat() {
        let mut loop_ = CoolingLoop::passenger_car();
        // Set coolant below thermostat opening so valve = 0 (no cooling)
        loop_.coolant_temp = loop_.thermostat_open_temp - 10.0;
        let t0 = loop_.coolant_temp;
        loop_.step(50_000.0, 1.0); // 50 kW
        assert!(
            loop_.coolant_temp > t0,
            "high engine heat should raise temp"
        );
    }

    #[test]
    fn cooling_loop_fields_accessible() {
        let loop_ = CoolingLoop::passenger_car();
        assert!(loop_.flow_rate > 0.0);
        assert!(loop_.specific_heat > 0.0);
    }

    // ── radiator_effectiveness ───────────────────────────────────────────

    #[test]
    fn effectiveness_zero_ntu_is_zero() {
        assert_eq!(radiator_effectiveness(0.0, 0.5), 0.0);
    }

    #[test]
    fn effectiveness_large_ntu_approaches_one() {
        let eps = radiator_effectiveness(20.0, 0.1);
        assert!(
            eps > 0.95,
            "high NTU should give effectiveness near 1: {eps}"
        );
    }

    #[test]
    fn effectiveness_balanced_flow_ntu1() {
        // Cr=1: ε = NTU/(1+NTU) = 1/2 for NTU=1
        let eps = radiator_effectiveness(1.0, 1.0);
        assert!((eps - 0.5).abs() < 1e-9, "balanced flow: ε={eps}");
    }

    #[test]
    fn effectiveness_zero_cr_analytical() {
        // Cr→0: ε = 1 − exp(−NTU)
        let ntu = 2.0;
        let eps = radiator_effectiveness(ntu, 0.0);
        let expected = 1.0 - (-ntu).exp();
        assert!((eps - expected).abs() < 1e-9);
    }

    #[test]
    fn effectiveness_clamped_to_unit_interval() {
        let eps = radiator_effectiveness(100.0, 0.5);
        assert!((0.0..=1.0).contains(&eps), "ε={eps}");
    }

    #[test]
    fn effectiveness_increases_with_ntu() {
        let e1 = radiator_effectiveness(1.0, 0.5);
        let e2 = radiator_effectiveness(2.0, 0.5);
        let e3 = radiator_effectiveness(4.0, 0.5);
        assert!(e1 < e2, "ε should increase with NTU");
        assert!(e2 < e3);
    }

    // ── heat_rejection_rate ──────────────────────────────────────────────

    #[test]
    fn heat_rejection_zero_effectiveness() {
        let q = heat_rejection_rate(0.0, 1000.0, 400.0, 300.0);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn heat_rejection_perfect_exchanger() {
        // ε=1, C_min=100, ΔT=50 → Q=5000
        let q = heat_rejection_rate(1.0, 100.0, 350.0, 300.0);
        assert!((q - 5000.0).abs() < 1e-9);
    }

    #[test]
    fn heat_rejection_negative_delta_t_is_zero() {
        // Cold > Hot → no heat transfer
        let q = heat_rejection_rate(0.8, 1000.0, 290.0, 310.0);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn heat_rejection_scales_with_c_min() {
        let q1 = heat_rejection_rate(0.5, 100.0, 400.0, 300.0);
        let q2 = heat_rejection_rate(0.5, 200.0, 400.0, 300.0);
        assert!((q2 / q1 - 2.0).abs() < 1e-9);
    }

    // ── thermostat_valve ─────────────────────────────────────────────────

    #[test]
    fn thermostat_below_open_temp_is_closed() {
        let v = thermostat_valve(350.0, 360.0, 375.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn thermostat_above_full_temp_is_open() {
        let v = thermostat_valve(380.0, 360.0, 375.0);
        assert_eq!(v, 1.0);
    }

    #[test]
    fn thermostat_midpoint_is_half() {
        let v = thermostat_valve(367.5, 360.0, 375.0);
        assert!((v - 0.5).abs() < 1e-9, "valve={v}");
    }

    #[test]
    fn thermostat_at_open_temp_is_zero() {
        let v = thermostat_valve(360.0, 360.0, 375.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn thermostat_at_full_temp_is_one() {
        let v = thermostat_valve(375.0, 360.0, 375.0);
        assert_eq!(v, 1.0);
    }

    #[test]
    fn thermostat_equal_temps_above_open_is_one() {
        // Degenerate case: t_coolant > t_open but range == 0 → fully open
        let v = thermostat_valve(371.0, 370.0, 370.0);
        assert_eq!(v, 1.0);
    }

    // ── coolant_temperature_ode ──────────────────────────────────────────

    #[test]
    fn ode_zero_heat_exchange_no_change() {
        let dt = coolant_temperature_ode(10_000.0, 5_000.0, 5_000.0, 1.0);
        assert!(dt.abs() < 1e-12);
    }

    #[test]
    fn ode_positive_net_heat_raises_temp() {
        let dt = coolant_temperature_ode(10_000.0, 8_000.0, 3_000.0, 1.0);
        // (8000-3000)/10000 * 1 = 0.5 K
        assert!((dt - 0.5).abs() < 1e-9);
    }

    #[test]
    fn ode_negative_net_heat_lowers_temp() {
        let dt = coolant_temperature_ode(10_000.0, 1_000.0, 6_000.0, 1.0);
        assert!(dt < 0.0);
    }

    #[test]
    fn ode_scales_with_dt() {
        let dt1 = coolant_temperature_ode(10_000.0, 5_000.0, 0.0, 1.0);
        let dt2 = coolant_temperature_ode(10_000.0, 5_000.0, 0.0, 2.0);
        assert!((dt2 / dt1 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn ode_zero_capacity_no_change() {
        let dt = coolant_temperature_ode(0.0, 5_000.0, 0.0, 1.0);
        assert_eq!(dt, 0.0);
    }

    // ── fan_curve ────────────────────────────────────────────────────────

    #[test]
    fn fan_zero_flow_max_pressure() {
        let pt = fan_curve(0.0, 500.0, 1.0, 0.7);
        assert!((pt.pressure_rise - 500.0).abs() < 1e-9);
    }

    #[test]
    fn fan_free_delivery_zero_pressure() {
        let pt = fan_curve(1.0, 500.0, 1.0, 0.7);
        assert!(pt.pressure_rise.abs() < 1e-9);
    }

    #[test]
    fn fan_midpoint_pressure() {
        // At Q = Q_max/2: ΔP = ΔP_max * (1 - 0.25) = 0.75 * ΔP_max
        let pt = fan_curve(0.5, 400.0, 1.0, 0.8);
        assert!(
            (pt.pressure_rise - 300.0).abs() < 1e-9,
            "ΔP={}",
            pt.pressure_rise
        );
    }

    #[test]
    fn fan_power_positive_nonzero_flow() {
        let pt = fan_curve(0.5, 400.0, 1.0, 0.8);
        assert!(pt.power > 0.0, "fan power should be positive");
    }

    #[test]
    fn fan_clamps_to_q_max() {
        let pt = fan_curve(2.0, 500.0, 1.0, 0.7);
        assert!(pt.flow_rate <= 1.0);
    }

    // ── fan_system_operating_point ───────────────────────────────────────

    #[test]
    fn fan_system_op_point_positive() {
        let q = fan_system_operating_point(500.0, 1.0, 100.0);
        assert!(q > 0.0 && q < 1.0, "operating point q={q}");
    }

    #[test]
    fn fan_system_op_point_zero_resistance() {
        // Zero resistance → Q = Q_max
        let q = fan_system_operating_point(500.0, 1.0, 0.0);
        assert!((q - 1.0).abs() < 1e-6, "q={q}");
    }

    #[test]
    fn fan_system_op_point_high_resistance() {
        let q_lo = fan_system_operating_point(500.0, 1.0, 1000.0);
        let q_hi = fan_system_operating_point(500.0, 1.0, 10.0);
        assert!(q_lo < q_hi, "higher resistance → lower flow");
    }

    #[test]
    fn cooling_loop_step_thermostat_closed_no_cooling() {
        let mut loop_ = CoolingLoop::passenger_car();
        // Force temp below thermostat: valve = 0 → no radiator heat
        loop_.coolant_temp = loop_.thermostat_open_temp - 20.0;
        loop_.thermal_capacity = 1_000_000.0; // large → tiny change
        let t0 = loop_.coolant_temp;
        loop_.step(0.0, 1.0);
        // With valve=0 and no engine heat, temp should be exactly unchanged
        assert!((loop_.coolant_temp - t0).abs() < 1e-9);
    }
}
