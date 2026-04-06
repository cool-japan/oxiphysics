// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! EV charging station and grid interaction model.
//!
//! Implements charging level definitions, battery CC-CV charging curves,
//! grid demand management, charging schedulers, and V2G (vehicle-to-grid)
//! capabilities.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_vehicle::charging_station::{
//!     ChargingStation, ChargingLevel, cc_cv_charging_time, range_from_soc,
//! };
//!
//! let station = ChargingStation::new(ChargingLevel::Level2, 4);
//! assert_eq!(station.n_ports, 4);
//!
//! let t = cc_cv_charging_time(75.0, 0.2, 0.8, 1.0);
//! assert!(t > 0.0);
//!
//! let range = range_from_soc(0.8, 6.0, 75.0);
//! assert!((range - 360.0).abs() < 1e-6);
//! ```

#![allow(dead_code)]

// ── ChargingLevel ─────────────────────────────────────────────────────────────

/// EV charging level / standard defining the nominal power.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChargingLevel {
    /// Level 1 AC charging (~1.4 kW, 120 V / 12 A).
    Level1,
    /// Level 2 AC charging (~7.2 kW, 240 V / 30 A).
    Level2,
    /// DC fast charging at 50 kW.
    Dcfc50,
    /// DC fast charging at 150 kW.
    Dcfc150,
    /// DC fast charging at 350 kW (ultra-rapid).
    Dcfc350,
}

impl ChargingLevel {
    /// Nominal power \[kW\] for this charging level.
    pub fn nominal_power_kw(self) -> f64 {
        match self {
            ChargingLevel::Level1 => 1.4,
            ChargingLevel::Level2 => 7.2,
            ChargingLevel::Dcfc50 => 50.0,
            ChargingLevel::Dcfc150 => 150.0,
            ChargingLevel::Dcfc350 => 350.0,
        }
    }

    /// Whether this level is a DC fast charger.
    pub fn is_dc_fast(self) -> bool {
        matches!(
            self,
            ChargingLevel::Dcfc50 | ChargingLevel::Dcfc150 | ChargingLevel::Dcfc350
        )
    }
}

// ── ChargingStation ───────────────────────────────────────────────────────────

/// A physical EV charging station with multiple ports.
#[derive(Debug, Clone)]
pub struct ChargingStation {
    /// Charging level (defines nominal power per port).
    pub level: ChargingLevel,
    /// Total number of charging ports.
    pub n_ports: usize,
    /// Number of currently available (unoccupied) ports.
    pub available_ports: usize,
    /// Actual maximum station power \[kW\] (may be lower than nominal × ports).
    pub power_kw: f64,
}

impl ChargingStation {
    /// Create a new station with the given level and port count.
    ///
    /// The total power is set to `level.nominal_power_kw() * n_ports`.
    pub fn new(level: ChargingLevel, n_ports: usize) -> Self {
        let power_kw = level.nominal_power_kw() * n_ports as f64;
        Self {
            level,
            n_ports,
            available_ports: n_ports,
            power_kw,
        }
    }

    /// Attempt to charge a vehicle battery for time step `dt` \[s\].
    ///
    /// Returns the new SoC after charging.  No port is actually reserved by
    /// this method – use in conjunction with `available_ports` tracking.
    ///
    /// * `vehicle_soc`      – current state of charge \[0, 1\].
    /// * `battery_capacity` – battery capacity \[kWh\].
    /// * `dt`               – time step \[s\].
    pub fn charge(&self, vehicle_soc: f64, battery_capacity: f64, dt: f64) -> f64 {
        if battery_capacity <= 0.0 || dt <= 0.0 {
            return vehicle_soc;
        }
        let power_per_port = self.power_kw / self.n_ports.max(1) as f64;
        let energy_kwh = power_per_port * (dt / 3600.0);
        let delta_soc = energy_kwh / battery_capacity;
        (vehicle_soc + delta_soc).min(1.0)
    }

    /// Reserve one port.  Returns `true` if a port was available.
    pub fn occupy_port(&mut self) -> bool {
        if self.available_ports > 0 {
            self.available_ports -= 1;
            true
        } else {
            false
        }
    }

    /// Release one port.
    pub fn release_port(&mut self) {
        if self.available_ports < self.n_ports {
            self.available_ports += 1;
        }
    }

    /// Whether at least one port is free.
    pub fn has_available_port(&self) -> bool {
        self.available_ports > 0
    }
}

// ── BatteryChargingModel ──────────────────────────────────────────────────────

/// CC-CV battery charging model.
///
/// During the constant-current (CC) phase power is limited by `max_power_kw`.
/// When SoC reaches `cv_start_soc` the model transitions to constant-voltage
/// (CV) and the accepted power tapers linearly to zero at SoC = 1.
#[derive(Debug, Clone)]
pub struct BatteryChargingModel {
    /// Total battery capacity \[kWh\].
    pub capacity_kwh: f64,
    /// Current state of charge \[0, 1\].
    pub soc: f64,
    /// Maximum charge power in CC phase \[kW\].
    pub max_power_kw: f64,
    /// SoC threshold at which CV phase begins (typically 0.8).
    pub cv_start_soc: f64,
    /// Round-trip charging efficiency (0–1).
    pub efficiency: f64,
}

impl BatteryChargingModel {
    /// Create a new battery model.
    pub fn new(capacity_kwh: f64, soc: f64, max_power_kw: f64) -> Self {
        Self {
            capacity_kwh,
            soc: soc.clamp(0.0, 1.0),
            max_power_kw,
            cv_start_soc: 0.8,
            efficiency: 0.95,
        }
    }

    /// Compute the accepted power \[kW\] at the current SoC.
    ///
    /// CC phase: full power.  CV phase: tapers linearly to zero at SoC = 1.
    pub fn accepted_power_kw(&self) -> f64 {
        if self.soc >= 1.0 {
            return 0.0;
        }
        if self.soc < self.cv_start_soc {
            self.max_power_kw
        } else {
            let cv_fraction = (1.0 - self.soc) / (1.0 - self.cv_start_soc);
            self.max_power_kw * cv_fraction
        }
    }

    /// Step the model for `dt` \[s\] given an available `power_kw` \[kW\].
    ///
    /// Returns the energy actually charged \[kWh\].
    pub fn step(&mut self, power_kw: f64, dt: f64) -> f64 {
        if self.capacity_kwh <= 0.0 || dt <= 0.0 {
            return 0.0;
        }
        let accepted = power_kw.min(self.accepted_power_kw());
        let energy_kwh = accepted * (dt / 3600.0) * self.efficiency;
        let delta_soc = energy_kwh / self.capacity_kwh;
        self.soc = (self.soc + delta_soc).min(1.0);
        energy_kwh
    }

    /// Whether the battery is considered fully charged (SoC ≥ 99%).
    pub fn is_full(&self) -> bool {
        self.soc >= 0.99
    }
}

// ── GridConnection ────────────────────────────────────────────────────────────

/// Electrical grid connection and demand management.
#[derive(Debug, Clone)]
pub struct GridConnection {
    /// Maximum subscribed demand \[kW\].
    pub max_demand_kw: f64,
    /// Current grid load \[kW\].
    pub current_load: f64,
    /// Demand charge rate \[$/kW\] per billing period.
    pub demand_charge_rate: f64,
    /// Peak-shaving threshold \[kW\].  Charging above this incurs extra cost.
    pub peak_shaving_threshold: f64,
}

impl GridConnection {
    /// Create a new grid connection.
    pub fn new(max_demand_kw: f64, demand_charge_rate: f64, peak_shaving_threshold: f64) -> Self {
        Self {
            max_demand_kw,
            current_load: 0.0,
            demand_charge_rate,
            peak_shaving_threshold,
        }
    }

    /// Available charging headroom \[kW\] below the peak-shaving threshold.
    pub fn available_charging_power(&self) -> f64 {
        (self.peak_shaving_threshold - self.current_load).max(0.0)
    }

    /// Update the current load and schedule charging accordingly.
    ///
    /// Returns the allowed charging power \[kW\] (limited by headroom and max demand).
    pub fn schedule_charging(&self, requested_kw: f64) -> f64 {
        let headroom = (self.max_demand_kw - self.current_load).max(0.0);
        let peak_headroom = self.available_charging_power();
        requested_kw.min(headroom).min(peak_headroom)
    }

    /// Compute the demand charge \[$\] for the billing period.
    ///
    /// * `peak_demand_kw` – recorded peak demand during the period \[kW\].
    pub fn demand_charge(&self, peak_demand_kw: f64) -> f64 {
        let excess = (peak_demand_kw - self.peak_shaving_threshold).max(0.0);
        excess * self.demand_charge_rate
    }
}

// ── ChargingScheduler ─────────────────────────────────────────────────────────

/// A charging session descriptor.
#[derive(Debug, Clone)]
pub struct ChargingSession {
    /// Vehicle arrival time \[s since epoch or simulation start\].
    pub arrival_time: f64,
    /// Vehicle departure time \[s\].
    pub departure_time: f64,
    /// Target SoC at departure \[0, 1\].
    pub soc_target: f64,
    /// Current SoC at arrival \[0, 1\].
    pub soc_init: f64,
    /// Battery capacity \[kWh\].
    pub capacity_kwh: f64,
}

/// Schedules multiple EV charging sessions to minimise peak demand.
#[derive(Debug, Clone)]
pub struct ChargingScheduler {
    /// Pending charging sessions.
    pub sessions: Vec<ChargingSession>,
    /// Grid connection for demand management.
    pub grid: GridConnection,
}

impl ChargingScheduler {
    /// Create a new empty scheduler.
    pub fn new(grid: GridConnection) -> Self {
        Self {
            sessions: Vec::new(),
            grid,
        }
    }

    /// Add a charging session to the queue.
    pub fn add_session(&mut self, session: ChargingSession) {
        self.sessions.push(session);
    }

    /// Simple greedy scheduler: sort sessions by departure time (earliest
    /// deadline first) and assign power up to the grid's available headroom.
    ///
    /// Returns a list of `(session_index, assigned_power_kw)` pairs.
    pub fn optimize_schedule(&self) -> Vec<(usize, f64)> {
        let mut order: Vec<usize> = (0..self.sessions.len()).collect();
        order.sort_by(|&a, &b| {
            self.sessions[a]
                .departure_time
                .partial_cmp(&self.sessions[b].departure_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut result = Vec::new();
        let mut remaining_headroom = self.grid.max_demand_kw - self.grid.current_load;

        for idx in order {
            let session = &self.sessions[idx];
            let energy_needed =
                (session.soc_target - session.soc_init).max(0.0) * session.capacity_kwh;
            let time_available = (session.departure_time - session.arrival_time).max(0.0);
            let power_needed = if time_available > 0.0 {
                energy_needed / (time_available / 3600.0)
            } else {
                0.0
            };

            let assigned = power_needed.min(remaining_headroom).max(0.0);
            remaining_headroom -= assigned;
            result.push((idx, assigned));
        }

        result
    }

    /// Total energy demand \[kWh\] across all pending sessions.
    pub fn total_energy_demand(&self) -> f64 {
        self.sessions
            .iter()
            .map(|s| (s.soc_target - s.soc_init).max(0.0) * s.capacity_kwh)
            .sum()
    }
}

// ── V2gCapability ─────────────────────────────────────────────────────────────

/// Vehicle-to-grid (V2G) discharge capability.
#[derive(Debug, Clone)]
pub struct V2gCapability {
    /// Maximum discharge power \[kW\].
    pub discharge_power_kw: f64,
    /// Minimum SoC at which discharge is permitted \[0, 1\].
    pub min_soc: f64,
    /// Discharge efficiency (0–1).
    pub efficiency: f64,
}

impl V2gCapability {
    /// Create a new V2G capability.
    pub fn new(discharge_power_kw: f64, min_soc: f64) -> Self {
        Self {
            discharge_power_kw,
            min_soc,
            efficiency: 0.93,
        }
    }

    /// Whether the battery can discharge at the current SoC.
    pub fn can_discharge(&self, soc: f64) -> bool {
        soc > self.min_soc
    }

    /// Compute the energy discharged to the grid \[kWh\] in `dt` \[s\].
    ///
    /// Returns the energy delivered (accounting for efficiency).
    /// Does not modify the battery SoC internally — that must be done by the
    /// caller.
    pub fn discharge(&self, dt: f64) -> f64 {
        let energy_from_battery = self.discharge_power_kw * (dt / 3600.0);
        energy_from_battery * self.efficiency
    }

    /// Compute the SoC delta when discharging for `dt` \[s\] from a battery of
    /// given capacity \[kWh\].
    pub fn delta_soc(&self, capacity_kwh: f64, dt: f64) -> f64 {
        if capacity_kwh <= 0.0 {
            return 0.0;
        }
        self.discharge_power_kw * (dt / 3600.0) / capacity_kwh
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Estimate the total CC-CV charging time \[s\].
///
/// * `capacity`   – battery capacity \[kWh\].
/// * `soc_init`   – initial SoC \[0, 1\].
/// * `soc_target` – target SoC \[0, 1\].
/// * `c_rate`     – charge rate in C (e.g., 1.0 = 1C charges full in 1 h).
/// * CV phase starts at SoC = 0.8 and the power tapers logarithmically.
///
/// Returns total time in seconds.
pub fn cc_cv_charging_time(capacity: f64, soc_init: f64, soc_target: f64, c_rate: f64) -> f64 {
    if soc_init >= soc_target || capacity <= 0.0 || c_rate <= 0.0 {
        return 0.0;
    }
    let cv_start = 0.8_f64;
    let power_cc = c_rate * capacity; // kW (when C=1 and capacity in kWh)
    let target = soc_target.min(1.0);

    // CC phase
    let soc_cc_end = cv_start.min(target);
    let t_cc = if soc_init < soc_cc_end {
        let energy_cc = (soc_cc_end - soc_init) * capacity;
        energy_cc / power_cc * 3600.0 // s
    } else {
        0.0
    };

    // CV phase (tapered power, approximated as ln taper)
    let soc_cv_start = soc_init.max(cv_start);
    let t_cv = if target > cv_start && soc_cv_start < target {
        // Use integral of 1/((1-soc)/(1-cv_start)) dSoC from soc_cv_start to target
        // Power in CV: P(soc) = P_cc * (1 - soc) / (1 - cv_start)
        // dt = capacity * dSoC / P(soc)
        // integral dt = capacity * (1 - cv_start) / P_cc
        //               * integral dSoC / (1-soc) from a to b
        //             = capacity*(1-cv_start)/P_cc * ln((1-a)/(1-b))
        let a = soc_cv_start;
        let b = target.min(0.9999); // avoid ln(0)
        let ratio = (1.0 - a) / (1.0 - b);
        if ratio <= 1.0 {
            0.0
        } else {
            capacity * (1.0 - cv_start) / power_cc * ratio.ln() * 3600.0
        }
    } else {
        0.0
    };

    t_cc + t_cv
}

/// Compute the achievable driving range \[km\] from a given SoC.
///
/// * `soc`                   – current SoC \[0, 1\].
/// * `efficiency_km_per_kwh` – vehicle energy efficiency \[km/kWh\].
/// * `capacity`              – battery capacity \[kWh\].
pub fn range_from_soc(soc: f64, efficiency_km_per_kwh: f64, capacity: f64) -> f64 {
    soc.max(0.0) * capacity * efficiency_km_per_kwh
}

/// Compute the SoC required to cover a given range \[km\].
///
/// Returns a value in \[0, 1\], clamped.
pub fn soc_for_range(range_km: f64, efficiency_km_per_kwh: f64, capacity: f64) -> f64 {
    if efficiency_km_per_kwh <= 0.0 || capacity <= 0.0 {
        return 1.0;
    }
    (range_km / (efficiency_km_per_kwh * capacity)).min(1.0)
}

/// Compute the cost of a charging session \[$\].
///
/// * `energy_kwh`      – energy charged \[kWh\].
/// * `price_per_kwh`   – electricity price \[$/kWh\].
pub fn charging_cost(energy_kwh: f64, price_per_kwh: f64) -> f64 {
    energy_kwh * price_per_kwh
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ChargingLevel ────────────────────────────────────────────────────────

    #[test]
    fn level1_power_is_1_4_kw() {
        assert!((ChargingLevel::Level1.nominal_power_kw() - 1.4).abs() < 1e-10);
    }

    #[test]
    fn level2_power_is_7_2_kw() {
        assert!((ChargingLevel::Level2.nominal_power_kw() - 7.2).abs() < 1e-10);
    }

    #[test]
    fn dcfc50_power_is_50_kw() {
        assert!((ChargingLevel::Dcfc50.nominal_power_kw() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn dcfc150_power_is_150_kw() {
        assert!((ChargingLevel::Dcfc150.nominal_power_kw() - 150.0).abs() < 1e-10);
    }

    #[test]
    fn dcfc350_power_is_350_kw() {
        assert!((ChargingLevel::Dcfc350.nominal_power_kw() - 350.0).abs() < 1e-10);
    }

    #[test]
    fn level1_is_not_dc_fast() {
        assert!(!ChargingLevel::Level1.is_dc_fast());
    }

    #[test]
    fn dcfc50_is_dc_fast() {
        assert!(ChargingLevel::Dcfc50.is_dc_fast());
    }

    // ── ChargingStation ──────────────────────────────────────────────────────

    #[test]
    fn station_new_sets_fields() {
        let s = ChargingStation::new(ChargingLevel::Level2, 4);
        assert_eq!(s.n_ports, 4);
        assert_eq!(s.available_ports, 4);
        assert!((s.power_kw - 7.2 * 4.0).abs() < 1e-10);
    }

    #[test]
    fn station_charge_increases_soc() {
        let s = ChargingStation::new(ChargingLevel::Level2, 1);
        let new_soc = s.charge(0.5, 75.0, 3600.0);
        assert!(new_soc > 0.5, "new_soc={new_soc}");
    }

    #[test]
    fn station_charge_capped_at_one() {
        let s = ChargingStation::new(ChargingLevel::Dcfc350, 1);
        // Very long dt → SoC capped at 1.0
        let new_soc = s.charge(0.9, 10.0, 1_000_000.0);
        assert!((new_soc - 1.0).abs() < 1e-10, "new_soc={new_soc}");
    }

    #[test]
    fn station_occupy_and_release_ports() {
        let mut s = ChargingStation::new(ChargingLevel::Level2, 2);
        assert!(s.occupy_port());
        assert!(s.occupy_port());
        assert!(!s.occupy_port()); // no ports left
        s.release_port();
        assert!(s.has_available_port());
    }

    #[test]
    fn station_has_available_port_initially() {
        let s = ChargingStation::new(ChargingLevel::Level2, 3);
        assert!(s.has_available_port());
    }

    // ── BatteryChargingModel ─────────────────────────────────────────────────

    #[test]
    fn battery_cc_phase_full_power() {
        let b = BatteryChargingModel::new(75.0, 0.5, 50.0);
        assert!((b.accepted_power_kw() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn battery_cv_phase_tapered_power() {
        let b = BatteryChargingModel::new(75.0, 0.9, 50.0);
        let p = b.accepted_power_kw();
        // at soc=0.9, fraction = (1-0.9)/(1-0.8) = 0.5 → power = 25
        assert!((p - 25.0).abs() < 1e-10, "p={p}");
    }

    #[test]
    fn battery_full_soc_zero_power() {
        let b = BatteryChargingModel::new(75.0, 1.0, 50.0);
        assert_eq!(b.accepted_power_kw(), 0.0);
    }

    #[test]
    fn battery_step_increases_soc() {
        let mut b = BatteryChargingModel::new(75.0, 0.5, 50.0);
        let energy = b.step(50.0, 3600.0);
        assert!(b.soc > 0.5, "soc={}", b.soc);
        assert!(energy > 0.0, "energy={energy}");
    }

    #[test]
    fn battery_step_does_not_exceed_full() {
        let mut b = BatteryChargingModel::new(10.0, 0.99, 100.0);
        b.step(100.0, 3600.0);
        assert!(b.soc <= 1.0, "soc={}", b.soc);
    }

    #[test]
    fn battery_is_full_when_soc_high() {
        let mut b = BatteryChargingModel::new(75.0, 0.99, 50.0);
        b.soc = 0.99;
        assert!(b.is_full());
    }

    #[test]
    fn battery_is_not_full_when_soc_low() {
        let b = BatteryChargingModel::new(75.0, 0.5, 50.0);
        assert!(!b.is_full());
    }

    // ── GridConnection ────────────────────────────────────────────────────────

    #[test]
    fn grid_available_charging_power_below_threshold() {
        let mut g = GridConnection::new(500.0, 10.0, 200.0);
        g.current_load = 150.0;
        assert!((g.available_charging_power() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn grid_available_power_zero_when_at_threshold() {
        let mut g = GridConnection::new(500.0, 10.0, 200.0);
        g.current_load = 200.0;
        assert_eq!(g.available_charging_power(), 0.0);
    }

    #[test]
    fn grid_schedule_charging_limits_to_headroom() {
        let mut g = GridConnection::new(100.0, 5.0, 80.0);
        g.current_load = 70.0;
        let allowed = g.schedule_charging(50.0);
        assert!(allowed <= 10.0 + 1e-10, "allowed={allowed}");
    }

    #[test]
    fn grid_demand_charge_zero_below_threshold() {
        let g = GridConnection::new(500.0, 10.0, 200.0);
        assert_eq!(g.demand_charge(199.0), 0.0);
    }

    #[test]
    fn grid_demand_charge_positive_above_threshold() {
        let g = GridConnection::new(500.0, 10.0, 200.0);
        let charge = g.demand_charge(250.0);
        assert!((charge - 500.0).abs() < 1e-10, "charge={charge}"); // 50 kW excess × $10/kW
    }

    // ── ChargingScheduler ─────────────────────────────────────────────────────

    #[test]
    fn scheduler_total_energy_demand_correct() {
        let grid = GridConnection::new(500.0, 10.0, 400.0);
        let mut sched = ChargingScheduler::new(grid);
        sched.add_session(ChargingSession {
            arrival_time: 0.0,
            departure_time: 3600.0,
            soc_target: 0.8,
            soc_init: 0.2,
            capacity_kwh: 75.0,
        });
        // energy = (0.8 - 0.2) * 75 = 45 kWh
        assert!((sched.total_energy_demand() - 45.0).abs() < 1e-10);
    }

    #[test]
    fn scheduler_optimize_returns_assignment_per_session() {
        let grid = GridConnection::new(500.0, 10.0, 400.0);
        let mut sched = ChargingScheduler::new(grid);
        sched.add_session(ChargingSession {
            arrival_time: 0.0,
            departure_time: 7200.0,
            soc_target: 0.9,
            soc_init: 0.3,
            capacity_kwh: 60.0,
        });
        sched.add_session(ChargingSession {
            arrival_time: 0.0,
            departure_time: 3600.0,
            soc_target: 0.8,
            soc_init: 0.5,
            capacity_kwh: 80.0,
        });
        let plan = sched.optimize_schedule();
        assert_eq!(plan.len(), 2);
        for (_, power) in &plan {
            assert!(*power >= 0.0, "power={power}");
        }
    }

    #[test]
    fn scheduler_empty_sessions_zero_energy() {
        let grid = GridConnection::new(500.0, 10.0, 400.0);
        let sched = ChargingScheduler::new(grid);
        assert_eq!(sched.total_energy_demand(), 0.0);
        assert!(sched.optimize_schedule().is_empty());
    }

    // ── V2gCapability ─────────────────────────────────────────────────────────

    #[test]
    fn v2g_can_discharge_above_min_soc() {
        let v2g = V2gCapability::new(10.0, 0.2);
        assert!(v2g.can_discharge(0.5));
    }

    #[test]
    fn v2g_cannot_discharge_at_min_soc() {
        let v2g = V2gCapability::new(10.0, 0.2);
        assert!(!v2g.can_discharge(0.2));
    }

    #[test]
    fn v2g_discharge_energy_positive() {
        let v2g = V2gCapability::new(10.0, 0.2);
        let energy = v2g.discharge(3600.0);
        // 10 kW × 1 h × 0.93 efficiency = 9.3 kWh
        assert!((energy - 9.3).abs() < 1e-6, "energy={energy}");
    }

    #[test]
    fn v2g_discharge_zero_dt() {
        let v2g = V2gCapability::new(10.0, 0.2);
        assert_eq!(v2g.discharge(0.0), 0.0);
    }

    #[test]
    fn v2g_delta_soc_correct() {
        let v2g = V2gCapability::new(10.0, 0.2);
        // delta_soc = 10 * (3600/3600) / 100 = 0.1
        let ds = v2g.delta_soc(100.0, 3600.0);
        assert!((ds - 0.1).abs() < 1e-10, "ds={ds}");
    }

    #[test]
    fn v2g_delta_soc_zero_capacity() {
        let v2g = V2gCapability::new(10.0, 0.2);
        assert_eq!(v2g.delta_soc(0.0, 3600.0), 0.0);
    }

    // ── cc_cv_charging_time ───────────────────────────────────────────────────

    #[test]
    fn cc_cv_time_positive_for_valid_inputs() {
        let t = cc_cv_charging_time(75.0, 0.2, 0.8, 1.0);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn cc_cv_time_zero_when_already_at_target() {
        let t = cc_cv_charging_time(75.0, 0.8, 0.8, 1.0);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn cc_cv_time_zero_for_zero_c_rate() {
        let t = cc_cv_charging_time(75.0, 0.2, 0.8, 0.0);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn cc_cv_time_less_for_higher_c_rate() {
        let t1 = cc_cv_charging_time(75.0, 0.2, 0.8, 0.5);
        let t2 = cc_cv_charging_time(75.0, 0.2, 0.8, 1.0);
        assert!(t2 < t1, "t1={t1} t2={t2}");
    }

    #[test]
    fn cc_cv_time_includes_cv_phase() {
        // From 0.2 to 0.9 → must include CV phase
        let t_no_cv = cc_cv_charging_time(75.0, 0.2, 0.8, 1.0);
        let t_with_cv = cc_cv_charging_time(75.0, 0.2, 0.9, 1.0);
        assert!(
            t_with_cv > t_no_cv,
            "t_no_cv={t_no_cv} t_with_cv={t_with_cv}"
        );
    }

    // ── range_from_soc ────────────────────────────────────────────────────────

    #[test]
    fn range_from_soc_full_battery() {
        let range = range_from_soc(1.0, 6.0, 75.0);
        assert!((range - 450.0).abs() < 1e-6, "range={range}");
    }

    #[test]
    fn range_from_soc_half_battery() {
        let range = range_from_soc(0.5, 6.0, 75.0);
        assert!((range - 225.0).abs() < 1e-6, "range={range}");
    }

    #[test]
    fn range_from_soc_known_value() {
        // soc=0.8, eff=6 km/kWh, cap=75 kWh → 360 km
        let range = range_from_soc(0.8, 6.0, 75.0);
        assert!((range - 360.0).abs() < 1e-6, "range={range}");
    }

    #[test]
    fn range_from_soc_zero_soc() {
        let range = range_from_soc(0.0, 6.0, 75.0);
        assert_eq!(range, 0.0);
    }

    // ── soc_for_range ─────────────────────────────────────────────────────────

    #[test]
    fn soc_for_range_basic() {
        // 300 km, 6 km/kWh, 75 kWh → 300/(6*75) = 0.6667
        let soc = soc_for_range(300.0, 6.0, 75.0);
        assert!((soc - 300.0 / 450.0).abs() < 1e-10, "soc={soc}");
    }

    #[test]
    fn soc_for_range_exceeds_full_is_clamped() {
        let soc = soc_for_range(10_000.0, 6.0, 75.0);
        assert!(soc <= 1.0, "soc={soc}");
    }

    // ── charging_cost ─────────────────────────────────────────────────────────

    #[test]
    fn charging_cost_basic() {
        let cost = charging_cost(45.0, 0.25);
        assert!((cost - 11.25).abs() < 1e-10, "cost={cost}");
    }

    #[test]
    fn charging_cost_zero_energy() {
        assert_eq!(charging_cost(0.0, 0.30), 0.0);
    }
}
