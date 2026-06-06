// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle energy recovery systems.
//!
//! Covers:
//! - [`BatteryModel`] — electrochemical state-of-charge model with OCV table
//! - [`RegenerativeBraking`] — regen torque blending and recoverable power
//! - [`KersSystem`] — flywheel-based kinetic energy recovery
//! - [`FuelCell`] — hydrogen fuel cell power module
//! - [`HybridPowerTrain`] — parallel hybrid power split and fuel consumption
//! - [`willans_line_efficiency`] — analytical engine efficiency model
//! - [`optimal_operating_line`] — torque-vs-RPM OOL for minimum fuel

// ---------------------------------------------------------------------------
// BatteryModel
// ---------------------------------------------------------------------------

/// Electrochemical battery model with state-of-charge tracking.
///
/// Uses an open-circuit-voltage (OCV) look-up table and a lumped internal
/// resistance to compute terminal voltage and energy throughput.
#[derive(Debug, Clone)]
pub struct BatteryModel {
    /// Nominal battery capacity (kWh).
    pub capacity_kwh: f64,
    /// Current state of charge, in the range \[0.0, 1.0\].
    pub soc: f64,
    /// Lumped internal resistance (Ω).
    pub internal_resistance: f64,
    /// OCV table: `Vec<(soc, voltage)>` sorted ascending by SOC.
    pub v_oc_table: Vec<(f64, f64)>,
}

impl BatteryModel {
    /// Create a default lithium-ion battery model.
    ///
    /// # Arguments
    /// * `capacity_kwh` — nominal capacity (kWh).
    /// * `initial_soc` — initial state of charge \[0.0, 1.0\].
    pub fn new(capacity_kwh: f64, initial_soc: f64) -> Self {
        // Typical NMC OCV curve
        let v_oc_table = vec![
            (0.0, 3.0),
            (0.1, 3.4),
            (0.2, 3.55),
            (0.3, 3.65),
            (0.4, 3.72),
            (0.5, 3.78),
            (0.6, 3.85),
            (0.7, 3.92),
            (0.8, 4.0),
            (0.9, 4.1),
            (1.0, 4.2),
        ];
        BatteryModel {
            capacity_kwh,
            soc: initial_soc.clamp(0.0, 1.0),
            internal_resistance: 0.05,
            v_oc_table,
        }
    }

    /// Create a battery model with a custom OCV table.
    ///
    /// # Arguments
    /// * `capacity_kwh` — nominal capacity (kWh).
    /// * `initial_soc` — initial state of charge \[0.0, 1.0\].
    /// * `internal_resistance` — lumped DC resistance (Ω).
    /// * `v_oc_table` — `(soc, voltage)` pairs sorted by SOC ascending.
    pub fn with_table(
        capacity_kwh: f64,
        initial_soc: f64,
        internal_resistance: f64,
        v_oc_table: Vec<(f64, f64)>,
    ) -> Self {
        BatteryModel {
            capacity_kwh,
            soc: initial_soc.clamp(0.0, 1.0),
            internal_resistance,
            v_oc_table,
        }
    }

    /// Interpolate the open-circuit voltage at a given SOC.
    ///
    /// Clamps `soc` to the table bounds and performs linear interpolation.
    pub fn interpolate_voc(&self, soc: f64) -> f64 {
        let table = &self.v_oc_table;
        if table.is_empty() {
            return 0.0;
        }
        let soc = soc.clamp(
            table.first().expect("collection should not be empty").0,
            table.last().expect("collection should not be empty").0,
        );
        for i in 1..table.len() {
            let (s0, v0) = table[i - 1];
            let (s1, v1) = table[i];
            if soc <= s1 {
                let t = (soc - s0) / (s1 - s0).max(1e-12);
                return v0 + t * (v1 - v0);
            }
        }
        table.last().expect("collection should not be empty").1
    }

    /// Compute the terminal voltage (V).
    ///
    /// Returns `V_oc(soc)` (no current correction since current is not tracked).
    pub fn voltage(&self) -> f64 {
        self.interpolate_voc(self.soc)
    }

    /// Charge the battery with `power_w` watts for `dt` seconds.
    ///
    /// Updates SOC according to coulomb counting.
    /// Clamps SOC to \[0.0, 1.0\].
    pub fn charge(&mut self, power_w: f64, dt: f64) {
        let capacity_wh = self.capacity_kwh * 1000.0;
        let delta_soc = (power_w * dt / 3600.0) / capacity_wh;
        self.soc = (self.soc + delta_soc).clamp(0.0, 1.0);
    }

    /// Discharge the battery delivering up to `power_w` watts for `dt` seconds.
    ///
    /// Returns the actual power delivered (W), which may be less than `power_w`
    /// if the battery is near empty.
    pub fn discharge(&mut self, power_w: f64, dt: f64) -> f64 {
        let capacity_wh = self.capacity_kwh * 1000.0;
        let max_energy_j = self.soc * capacity_wh * 3600.0;
        let requested_j = power_w * dt;
        let actual_j = requested_j.min(max_energy_j);
        let actual_power = actual_j / dt.max(1e-9);
        let delta_soc = actual_j / (capacity_wh * 3600.0);
        self.soc = (self.soc - delta_soc).clamp(0.0, 1.0);
        actual_power
    }

    /// Remaining energy in the battery (Wh).
    pub fn remaining_energy_wh(&self) -> f64 {
        self.soc * self.capacity_kwh * 1000.0
    }
}

// ---------------------------------------------------------------------------
// RegenerativeBraking
// ---------------------------------------------------------------------------

/// Regenerative braking module.
///
/// Blends regenerative and friction braking torques while respecting
/// maximum regen torque and power limits.
#[derive(Debug, Clone)]
pub struct RegenerativeBraking {
    /// Maximum regenerative braking torque (N·m).
    pub max_regen_torque: f64,
    /// Motor/generator mechanical-to-electrical efficiency \[0, 1\].
    pub efficiency: f64,
    /// Maximum regenerative power (kW).
    pub max_power_kw: f64,
}

impl RegenerativeBraking {
    /// Create a new regenerative braking module.
    ///
    /// # Arguments
    /// * `max_regen_torque` — peak regen torque (N·m).
    /// * `efficiency` — energy conversion efficiency \[0, 1\].
    /// * `max_power_kw` — peak electrical power limit (kW).
    pub fn new(max_regen_torque: f64, efficiency: f64, max_power_kw: f64) -> Self {
        RegenerativeBraking {
            max_regen_torque,
            efficiency: efficiency.clamp(0.0, 1.0),
            max_power_kw,
        }
    }

    /// Compute the regen and friction torque split for a braking event.
    ///
    /// # Arguments
    /// * `wheel_speed_rads` — current wheel angular velocity (rad/s).
    /// * `brake_demand` — normalised brake demand \[0.0, 1.0\].
    ///
    /// Returns `(regen_torque_nm, friction_torque_nm)`.
    pub fn regen_torque(&self, wheel_speed_rads: f64, brake_demand: f64) -> (f64, f64) {
        let demand = brake_demand.clamp(0.0, 1.0);
        if wheel_speed_rads < 0.1 {
            // Too slow for regen — all friction
            return (0.0, demand * self.max_regen_torque);
        }
        // Max regen limited by power
        let max_regen_by_power = if wheel_speed_rads > 1e-6 {
            (self.max_power_kw * 1000.0) / wheel_speed_rads
        } else {
            self.max_regen_torque
        };
        let available_regen = self.max_regen_torque.min(max_regen_by_power);
        let total_braking = demand * self.max_regen_torque;
        let regen = total_braking.min(available_regen);
        let friction = (total_braking - regen).max(0.0);
        (regen, friction)
    }

    /// Compute power recoverable from braking.
    ///
    /// # Arguments
    /// * `vehicle_speed` — vehicle velocity (m/s).
    /// * `decel_force` — braking force (N).
    ///
    /// Returns recoverable electrical power (W).
    pub fn recoverable_power(&self, vehicle_speed: f64, decel_force: f64) -> f64 {
        let mechanical_power = vehicle_speed * decel_force.abs();
        let available = mechanical_power.min(self.max_power_kw * 1000.0);
        available * self.efficiency
    }
}

// ---------------------------------------------------------------------------
// KersSystem
// ---------------------------------------------------------------------------

/// Flywheel-based Kinetic Energy Recovery System (KERS).
///
/// Stores kinetic energy in a spinning flywheel and can deploy it on demand.
#[derive(Debug, Clone)]
pub struct KersSystem {
    /// Flywheel moment of inertia (kg·m²).
    pub flywheel_inertia: f64,
    /// Current flywheel angular speed (rad/s).
    pub flywheel_speed: f64,
    /// Maximum deployment/storage power (kW).
    pub max_power_kw: f64,
    /// Round-trip energy conversion efficiency \[0, 1\].
    pub efficiency: f64,
}

impl KersSystem {
    /// Create a new KERS module.
    ///
    /// # Arguments
    /// * `flywheel_inertia` — moment of inertia of the flywheel (kg·m²).
    /// * `max_power_kw` — peak power in/out (kW).
    /// * `efficiency` — one-way conversion efficiency \[0, 1\].
    pub fn new(flywheel_inertia: f64, max_power_kw: f64, efficiency: f64) -> Self {
        KersSystem {
            flywheel_inertia,
            flywheel_speed: 0.0,
            max_power_kw,
            efficiency: efficiency.clamp(0.0, 1.0),
        }
    }

    /// Store power in the flywheel for `dt` seconds.
    ///
    /// # Arguments
    /// * `power` — input power (W).
    /// * `dt` — time step (s).
    pub fn store(&mut self, power: f64, dt: f64) {
        let stored_energy = power.min(self.max_power_kw * 1000.0) * self.efficiency * dt;
        // E = 0.5 * I * omega^2  =>  omega = sqrt(2E/I)
        let current_e = 0.5 * self.flywheel_inertia * self.flywheel_speed * self.flywheel_speed;
        let new_e = (current_e + stored_energy).max(0.0);
        self.flywheel_speed = (2.0 * new_e / self.flywheel_inertia.max(1e-12)).sqrt();
    }

    /// Deploy energy from the flywheel.
    ///
    /// # Arguments
    /// * `power_demand` — requested power output (W).
    /// * `dt` — time step (s).
    ///
    /// Returns actual power delivered (W).
    pub fn deploy(&mut self, power_demand: f64, dt: f64) -> f64 {
        let available_e = self.energy_stored();
        let requested_e = power_demand.min(self.max_power_kw * 1000.0) * dt;
        let actual_e = requested_e.min(available_e);
        let _actual_power = actual_e / dt.max(1e-9);
        let delivered_e = actual_e * self.efficiency;
        let new_e = (available_e - actual_e).max(0.0);
        self.flywheel_speed = (2.0 * new_e / self.flywheel_inertia.max(1e-12)).sqrt();
        delivered_e / dt.max(1e-9)
    }

    /// Energy currently stored in the flywheel (J).
    pub fn energy_stored(&self) -> f64 {
        0.5 * self.flywheel_inertia * self.flywheel_speed * self.flywheel_speed
    }

    /// Equivalent stored energy in Wh.
    pub fn energy_stored_wh(&self) -> f64 {
        self.energy_stored() / 3600.0
    }
}

// ---------------------------------------------------------------------------
// FuelCell
// ---------------------------------------------------------------------------

/// Hydrogen fuel cell power module.
///
/// Models power output as a function of hydrogen consumption and a
/// piece-wise-linear efficiency curve.
#[derive(Debug, Clone)]
pub struct FuelCell {
    /// Peak electrical power output (kW).
    pub max_power_kw: f64,
    /// Efficiency table: `Vec<(power_fraction, efficiency)>`, fractions in \[0,1\].
    pub efficiency_table: Vec<(f64, f64)>,
    /// Remaining hydrogen (kg).
    pub hydrogen_kg: f64,
}

impl FuelCell {
    /// Gravimetric energy density of hydrogen (J/kg).
    pub const H2_LHV_J_KG: f64 = 120.0e6;

    /// Create a new fuel cell model.
    ///
    /// # Arguments
    /// * `max_power_kw` — rated peak power (kW).
    /// * `hydrogen_kg` — initial hydrogen fuel load (kg).
    pub fn new(max_power_kw: f64, hydrogen_kg: f64) -> Self {
        let efficiency_table = vec![
            (0.0, 0.0),
            (0.1, 0.55),
            (0.2, 0.60),
            (0.4, 0.58),
            (0.6, 0.55),
            (0.8, 0.52),
            (1.0, 0.48),
        ];
        FuelCell {
            max_power_kw,
            efficiency_table,
            hydrogen_kg,
        }
    }

    /// Interpolate fuel cell efficiency at a given power fraction.
    ///
    /// # Arguments
    /// * `power_fraction` — fraction of `max_power_kw` demanded, in \[0, 1\].
    pub fn efficiency_at(&self, power_fraction: f64) -> f64 {
        let table = &self.efficiency_table;
        if table.is_empty() {
            return 0.5;
        }
        let pf = power_fraction.clamp(0.0, 1.0);
        for i in 1..table.len() {
            let (p0, e0) = table[i - 1];
            let (p1, e1) = table[i];
            if pf <= p1 {
                let t = (pf - p0) / (p1 - p0).max(1e-12);
                return e0 + t * (e1 - e0);
            }
        }
        table.last().map(|(_, e)| *e).unwrap_or(0.5)
    }

    /// Request power output from the fuel cell.
    ///
    /// Consumes hydrogen proportional to the chemical energy demanded.
    /// Returns actual electrical power output (W).
    ///
    /// # Arguments
    /// * `demand_kw` — desired electrical power (kW).
    /// * `dt` — time step (s).
    pub fn power_output(&mut self, demand_kw: f64, dt: f64) -> f64 {
        let demand_kw = demand_kw.min(self.max_power_kw).max(0.0);
        let power_fraction = demand_kw / self.max_power_kw.max(1e-12);
        let eta = self.efficiency_at(power_fraction);
        if eta < 1e-9 || self.hydrogen_kg <= 0.0 {
            return 0.0;
        }
        let electrical_w = demand_kw * 1000.0;
        let chemical_w = electrical_w / eta;
        let h2_consumed = chemical_w * dt / Self::H2_LHV_J_KG;
        let actual_h2 = h2_consumed.min(self.hydrogen_kg);
        self.hydrogen_kg -= actual_h2;

        actual_h2 * Self::H2_LHV_J_KG * eta / dt.max(1e-9)
    }

    /// Remaining hydrogen range estimate.
    ///
    /// # Arguments
    /// * `average_power_kw` — average power consumption (kW).
    ///
    /// Returns remaining operation time (s).
    pub fn remaining_time(&self, average_power_kw: f64) -> f64 {
        if average_power_kw < 1e-9 {
            return f64::INFINITY;
        }
        let eta = self.efficiency_at(average_power_kw / self.max_power_kw.max(1e-12));
        let chemical_w = average_power_kw * 1000.0 / eta.max(1e-9);
        self.hydrogen_kg * Self::H2_LHV_J_KG / chemical_w
    }
}

// ---------------------------------------------------------------------------
// HybridPowerTrain
// ---------------------------------------------------------------------------

/// Parallel hybrid powertrain combining ICE, battery, and regen braking.
#[derive(Debug, Clone)]
pub struct HybridPowerTrain {
    /// Peak ICE power (kW).
    pub ice_power_kw: f64,
    /// On-board battery.
    pub battery: BatteryModel,
    /// Regenerative braking module.
    pub regen: RegenerativeBraking,
}

impl HybridPowerTrain {
    /// Create a new hybrid powertrain.
    ///
    /// # Arguments
    /// * `ice_power_kw` — rated ICE peak power (kW).
    /// * `battery` — battery model.
    /// * `regen` — regen braking module.
    pub fn new(ice_power_kw: f64, battery: BatteryModel, regen: RegenerativeBraking) -> Self {
        HybridPowerTrain {
            ice_power_kw,
            battery,
            regen,
        }
    }

    /// Compute the power split between ICE, battery, and regen.
    ///
    /// Uses a rule-based energy management strategy:
    /// - If braking (`brake_force > 0`), route regen first.
    /// - Otherwise, prioritise battery if SOC > 0.3, else use ICE.
    ///
    /// # Arguments
    /// * `demand_kw` — total traction power demand (kW).
    /// * `speed_ms` — vehicle speed (m/s).
    /// * `brake_force` — braking force (N). Zero for driving, positive for braking.
    ///
    /// Returns `(ice_power_kw, battery_power_kw, regen_power_kw)`.
    pub fn split_power(&self, demand_kw: f64, speed_ms: f64, brake_force: f64) -> (f64, f64, f64) {
        if brake_force > 0.0 {
            let recoverable = self.regen.recoverable_power(speed_ms, brake_force) / 1000.0;
            return (0.0, 0.0, recoverable);
        }
        let demand = demand_kw.max(0.0);
        if self.battery.soc > 0.3 {
            let batt_fraction = ((self.battery.soc - 0.3) / 0.7).clamp(0.0, 1.0);
            let batt_max = (self.ice_power_kw * 0.5) * batt_fraction;
            let battery_power = demand.min(batt_max);
            let ice_power = (demand - battery_power).clamp(0.0, self.ice_power_kw);
            (ice_power, battery_power, 0.0)
        } else {
            let ice_power = demand.min(self.ice_power_kw);
            (ice_power, 0.0, 0.0)
        }
    }

    /// Estimate fuel consumption in L/100 km.
    ///
    /// Uses a simplified Willans-line model scaled to ICE power output.
    ///
    /// # Arguments
    /// * `power_kw` — effective ICE power (kW).
    /// * `speed_ms` — vehicle speed (m/s).
    pub fn fuel_consumption_l_per_100km(&self, power_kw: f64, speed_ms: f64) -> f64 {
        if speed_ms < 0.1 {
            return 0.0;
        }
        let bmep_fraction = (power_kw / self.ice_power_kw.max(1e-9)).clamp(0.0, 1.0);
        // Approximate BSFC (g/kWh) using Willans line
        let bsfc_g_kwh = 250.0 + 100.0 * (1.0 - bmep_fraction).powi(2);
        // Power to fuel flow: q_g_s = power_kw * bsfc_g_kwh / 3600
        let q_g_s = power_kw * bsfc_g_kwh / 3600.0;
        // Gasoline density ≈ 0.74 kg/L
        let q_l_s = q_g_s / 740.0;
        // L/100km = (L/s) / (m/s) * 100_000
        q_l_s / speed_ms * 100_000.0
    }
}

// ---------------------------------------------------------------------------
// Standalone functions
// ---------------------------------------------------------------------------

/// Willans-line engine efficiency model.
///
/// Computes indicated thermal efficiency from brake and friction mean
/// effective pressures using the Willans approximation:
///
/// η = BMEP / (BMEP + FMEP)
///
/// # Arguments
/// * `bmep` — brake mean effective pressure (Pa or relative).
/// * `fmep` — friction mean effective pressure (Pa or relative).
///
/// Returns efficiency in \[0, 1).
pub fn willans_line_efficiency(bmep: f64, fmep: f64) -> f64 {
    let denom = (bmep + fmep).abs();
    if denom < 1e-12 {
        return 0.0;
    }
    (bmep / denom).clamp(0.0, 1.0)
}

/// Compute the Optimal Operating Line (OOL) — torque vs. RPM pairs that
/// minimise fuel consumption for a given power demand.
///
/// Uses the maximum-torque-at-minimum-speed heuristic. Points are distributed
/// from idle RPM to `max_rpm`.
///
/// # Arguments
/// * `max_torque` — maximum engine torque (N·m).
/// * `max_rpm` — maximum engine speed (RPM).
/// * `n_points` — number of OOL points to generate.
///
/// Returns a `Vec<(torque_nm, rpm)>` sorted ascending by RPM.
pub fn optimal_operating_line(max_torque: f64, max_rpm: f64, n_points: usize) -> Vec<(f64, f64)> {
    if n_points == 0 {
        return Vec::new();
    }
    let idle_rpm = 800.0_f64.min(max_rpm * 0.1);
    let mut points = Vec::with_capacity(n_points);
    for k in 0..n_points {
        let t = k as f64 / (n_points - 1).max(1) as f64;
        let rpm = idle_rpm + t * (max_rpm - idle_rpm);
        // OOL: use peak torque at low speeds, then drop proportionally
        let torque_fraction = if t < 0.5 { 1.0 } else { 1.0 - (t - 0.5) * 0.8 };
        let torque = max_torque * torque_fraction.clamp(0.0, 1.0);
        points.push((torque, rpm));
    }
    points
}

/// Compute charging C-rate from power and capacity.
///
/// C-rate = power (kW) / capacity (kWh).
///
/// # Arguments
/// * `power_kw` — charging power (kW).
/// * `capacity_kwh` — battery capacity (kWh).
pub fn c_rate(power_kw: f64, capacity_kwh: f64) -> f64 {
    if capacity_kwh < 1e-12 {
        return 0.0;
    }
    power_kw / capacity_kwh
}

/// Compute the specific energy recovered per braking event.
///
/// # Arguments
/// * `vehicle_mass_kg` — vehicle mass (kg).
/// * `speed_initial_ms` — initial speed before braking (m/s).
/// * `speed_final_ms` — final speed after braking (m/s).
/// * `efficiency` — regen efficiency \[0, 1\].
///
/// Returns recovered energy (J).
pub fn braking_energy_recovered(
    vehicle_mass_kg: f64,
    speed_initial_ms: f64,
    speed_final_ms: f64,
    efficiency: f64,
) -> f64 {
    let ke_diff = 0.5
        * vehicle_mass_kg
        * (speed_initial_ms * speed_initial_ms - speed_final_ms * speed_final_ms);
    (ke_diff * efficiency.clamp(0.0, 1.0)).max(0.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── BatteryModel ──────────────────────────────────────────────────────

    #[test]
    fn battery_initial_soc_clamped() {
        let b = BatteryModel::new(60.0, 1.5);
        assert!((b.soc - 1.0).abs() < 1e-12);
    }

    #[test]
    fn battery_initial_soc_negative_clamped() {
        let b = BatteryModel::new(60.0, -0.1);
        assert!(b.soc.abs() < 1e-12);
    }

    #[test]
    fn battery_soc_starts_at_given_value() {
        let b = BatteryModel::new(60.0, 0.8);
        assert!((b.soc - 0.8).abs() < 1e-10);
    }

    #[test]
    fn battery_interpolate_voc_at_soc_0() {
        let b = BatteryModel::new(60.0, 0.0);
        let v = b.interpolate_voc(0.0);
        assert!((v - 3.0).abs() < 1e-6, "v={v}");
    }

    #[test]
    fn battery_interpolate_voc_at_soc_1() {
        let b = BatteryModel::new(60.0, 1.0);
        let v = b.interpolate_voc(1.0);
        assert!((v - 4.2).abs() < 1e-6, "v={v}");
    }

    #[test]
    fn battery_interpolate_voc_midpoint() {
        let b = BatteryModel::new(60.0, 0.5);
        let v = b.interpolate_voc(0.5);
        assert!(v > 3.0 && v < 4.2, "v={v}");
    }

    #[test]
    fn battery_voltage_uses_current_soc() {
        let b = BatteryModel::new(60.0, 0.5);
        let v = b.voltage();
        let expected = b.interpolate_voc(0.5);
        assert!((v - expected).abs() < 1e-10);
    }

    #[test]
    fn battery_charge_increases_soc() {
        let mut b = BatteryModel::new(60.0, 0.5);
        let soc_before = b.soc;
        b.charge(10000.0, 1.0); // 10 kW for 1 s
        assert!(b.soc > soc_before, "SOC should increase after charging");
    }

    #[test]
    fn battery_charge_clamps_at_full() {
        let mut b = BatteryModel::new(1.0, 0.99);
        b.charge(1_000_000.0, 100.0);
        assert!((b.soc - 1.0).abs() < 1e-10, "SOC should not exceed 1.0");
    }

    #[test]
    fn battery_discharge_decreases_soc() {
        let mut b = BatteryModel::new(60.0, 0.5);
        let soc_before = b.soc;
        b.discharge(10000.0, 1.0);
        assert!(b.soc < soc_before, "SOC should decrease after discharge");
    }

    #[test]
    fn battery_discharge_clamps_at_empty() {
        let mut b = BatteryModel::new(1.0, 0.01);
        b.discharge(1_000_000.0, 100.0);
        assert!(b.soc >= 0.0, "SOC must not go below 0");
    }

    #[test]
    fn battery_discharge_returns_bounded_power() {
        let mut b = BatteryModel::new(0.001, 0.5); // very small battery
        let delivered = b.discharge(1_000_000.0, 1.0);
        assert!(delivered <= 1_000_000.0, "power cannot exceed demand");
    }

    #[test]
    fn battery_remaining_energy_proportional_to_soc() {
        let b = BatteryModel::new(60.0, 0.5);
        assert!((b.remaining_energy_wh() - 30_000.0).abs() < 1e-6);
    }

    #[test]
    fn battery_interpolate_voc_empty_table() {
        let b = BatteryModel::with_table(10.0, 0.5, 0.1, vec![]);
        let v = b.interpolate_voc(0.5);
        assert!((v).abs() < 1e-10);
    }

    // ── RegenerativeBraking ───────────────────────────────────────────────

    #[test]
    fn regen_torque_zero_demand() {
        let r = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let (rt, ft) = r.regen_torque(100.0, 0.0);
        assert!(rt.abs() < 1e-10);
        assert!(ft.abs() < 1e-10);
    }

    #[test]
    fn regen_torque_slow_speed_all_friction() {
        let r = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let (rt, _ft) = r.regen_torque(0.05, 1.0);
        assert!(rt.abs() < 1e-10, "no regen below min speed");
    }

    #[test]
    fn regen_torque_sum_equals_total_demand() {
        let r = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let demand = 0.7;
        let total = demand * r.max_regen_torque;
        let (rt, ft) = r.regen_torque(50.0, demand);
        assert!(
            (rt + ft - total).abs() < 1e-6,
            "rt+ft should equal total demand"
        );
    }

    #[test]
    fn regen_torque_respects_power_limit() {
        let r = RegenerativeBraking::new(10000.0, 0.9, 1.0); // very low power limit
        let (rt, _) = r.regen_torque(10.0, 1.0);
        // At 10 rad/s and 1 kW limit: max torque = 1000/10 = 100 Nm
        assert!(
            rt <= 1000.0 / 10.0 + 1e-6,
            "regen torque={rt} exceeds power limit"
        );
    }

    #[test]
    fn recoverable_power_proportional_to_speed() {
        let r = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let p1 = r.recoverable_power(10.0, 1000.0);
        let p2 = r.recoverable_power(20.0, 1000.0);
        assert!(p2 >= p1, "more speed = more recoverable power");
    }

    #[test]
    fn recoverable_power_limited_by_max() {
        let r = RegenerativeBraking::new(200.0, 0.9, 1.0); // 1 kW max
        let p = r.recoverable_power(100.0, 10000.0);
        assert!(p <= 1000.0 * 0.9 + 1e-6, "power={p} should be capped");
    }

    // ── KersSystem ────────────────────────────────────────────────────────

    #[test]
    fn kers_initial_energy_zero() {
        let k = KersSystem::new(0.5, 60.0, 0.9);
        assert!(k.energy_stored().abs() < 1e-12);
    }

    #[test]
    fn kers_store_increases_energy() {
        let mut k = KersSystem::new(0.5, 60.0, 0.9);
        k.store(10000.0, 1.0);
        assert!(k.energy_stored() > 0.0);
    }

    #[test]
    fn kers_deploy_decreases_energy() {
        let mut k = KersSystem::new(0.5, 60.0, 0.9);
        k.store(10000.0, 5.0);
        let e_before = k.energy_stored();
        k.deploy(5000.0, 1.0);
        assert!(k.energy_stored() < e_before);
    }

    #[test]
    fn kers_deploy_limited_by_stored_energy() {
        let mut k = KersSystem::new(0.5, 60.0, 0.9);
        k.store(100.0, 0.01); // store small amount
        let power = k.deploy(1_000_000.0, 1.0);
        assert!(power <= 1_000_000.0, "cannot deliver more than stored");
    }

    #[test]
    fn kers_energy_in_wh() {
        let mut k = KersSystem::new(1.0, 100.0, 1.0);
        k.store(3600.0, 1.0); // 3600 J = 1 Wh
        assert!(k.energy_stored_wh() > 0.0);
    }

    #[test]
    fn kers_flywheel_speed_nonnegative() {
        let mut k = KersSystem::new(1.0, 100.0, 0.9);
        k.store(50000.0, 2.0);
        k.deploy(100_000.0, 10.0);
        assert!(k.flywheel_speed >= 0.0);
    }

    // ── FuelCell ──────────────────────────────────────────────────────────

    #[test]
    fn fuel_cell_efficiency_at_zero_is_zero() {
        let fc = FuelCell::new(100.0, 5.0);
        let eta = fc.efficiency_at(0.0);
        assert!(eta.abs() < 1e-10);
    }

    #[test]
    fn fuel_cell_efficiency_at_full_load() {
        let fc = FuelCell::new(100.0, 5.0);
        let eta = fc.efficiency_at(1.0);
        assert!(eta > 0.0 && eta <= 1.0, "eta={eta}");
    }

    #[test]
    fn fuel_cell_power_output_reduces_hydrogen() {
        let mut fc = FuelCell::new(100.0, 5.0);
        let h2_before = fc.hydrogen_kg;
        fc.power_output(50.0, 1.0);
        assert!(fc.hydrogen_kg < h2_before, "H2 should be consumed");
    }

    #[test]
    fn fuel_cell_zero_demand_no_h2_consumed() {
        let mut fc = FuelCell::new(100.0, 5.0);
        let h2_before = fc.hydrogen_kg;
        fc.power_output(0.0, 1.0);
        assert!((fc.hydrogen_kg - h2_before).abs() < 1e-12);
    }

    #[test]
    fn fuel_cell_no_power_when_empty() {
        let mut fc = FuelCell::new(100.0, 0.0);
        let p = fc.power_output(50.0, 1.0);
        assert!(p.abs() < 1e-10, "no power without H2");
    }

    #[test]
    fn fuel_cell_remaining_time_infinite_at_zero_demand() {
        let fc = FuelCell::new(100.0, 5.0);
        let t = fc.remaining_time(0.0);
        assert!(t.is_infinite());
    }

    #[test]
    fn fuel_cell_remaining_time_positive() {
        let fc = FuelCell::new(100.0, 5.0);
        let t = fc.remaining_time(50.0);
        assert!(t > 0.0, "remaining time should be positive");
    }

    // ── HybridPowerTrain ──────────────────────────────────────────────────

    #[test]
    fn hybrid_split_braking_gives_regen() {
        let batt = BatteryModel::new(60.0, 0.5);
        let regen = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let pt = HybridPowerTrain::new(150.0, batt, regen);
        let (ice, battery, regen_p) = pt.split_power(0.0, 20.0, 5000.0);
        assert!(regen_p > 0.0, "regen should be nonzero during braking");
        assert!(ice.abs() < 1e-10);
        assert!(battery.abs() < 1e-10);
    }

    #[test]
    fn hybrid_split_high_soc_uses_battery() {
        let batt = BatteryModel::new(60.0, 0.9); // high SOC
        let regen = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let pt = HybridPowerTrain::new(150.0, batt, regen);
        let (ice, battery, _) = pt.split_power(50.0, 20.0, 0.0);
        assert!(battery > 0.0 || ice > 0.0);
    }

    #[test]
    fn hybrid_split_low_soc_uses_ice() {
        let batt = BatteryModel::new(60.0, 0.1); // low SOC
        let regen = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let pt = HybridPowerTrain::new(150.0, batt, regen);
        let (ice, battery, _) = pt.split_power(80.0, 20.0, 0.0);
        assert!(ice > 0.0, "ICE should be used when SOC is low");
        assert!(
            battery.abs() < 1e-10,
            "battery should not be used at low SOC"
        );
    }

    #[test]
    fn hybrid_fuel_consumption_zero_speed() {
        let batt = BatteryModel::new(60.0, 0.5);
        let regen = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let pt = HybridPowerTrain::new(150.0, batt, regen);
        let fc = pt.fuel_consumption_l_per_100km(50.0, 0.0);
        assert!(fc.abs() < 1e-10);
    }

    #[test]
    fn hybrid_fuel_consumption_positive_at_speed() {
        let batt = BatteryModel::new(60.0, 0.5);
        let regen = RegenerativeBraking::new(200.0, 0.9, 50.0);
        let pt = HybridPowerTrain::new(150.0, batt, regen);
        let fc = pt.fuel_consumption_l_per_100km(50.0, 30.0);
        assert!(fc > 0.0, "fuel consumption should be positive");
    }

    // ── willans_line_efficiency ───────────────────────────────────────────

    #[test]
    fn willans_efficiency_zero_bmep() {
        let eta = willans_line_efficiency(0.0, 10.0);
        assert!(eta.abs() < 1e-10);
    }

    #[test]
    fn willans_efficiency_equal_pressures() {
        let eta = willans_line_efficiency(10.0, 10.0);
        assert!((eta - 0.5).abs() < 1e-10);
    }

    #[test]
    fn willans_efficiency_high_bmep() {
        let eta = willans_line_efficiency(90.0, 10.0);
        assert!(eta > 0.5, "eta={eta}");
        assert!(eta <= 1.0, "eta={eta} must not exceed 1");
    }

    #[test]
    fn willans_efficiency_zero_denom() {
        let eta = willans_line_efficiency(0.0, 0.0);
        assert!(eta.abs() < 1e-10);
    }

    // ── optimal_operating_line ────────────────────────────────────────────

    #[test]
    fn ool_returns_n_points() {
        let ool = optimal_operating_line(400.0, 6000.0, 10);
        assert_eq!(ool.len(), 10);
    }

    #[test]
    fn ool_zero_points_returns_empty() {
        let ool = optimal_operating_line(400.0, 6000.0, 0);
        assert!(ool.is_empty());
    }

    #[test]
    fn ool_rpm_ascending() {
        let ool = optimal_operating_line(400.0, 6000.0, 8);
        for i in 1..ool.len() {
            assert!(ool[i].1 >= ool[i - 1].1, "RPM should be ascending");
        }
    }

    #[test]
    fn ool_torque_nonnegative() {
        let ool = optimal_operating_line(400.0, 6000.0, 10);
        for (t, rpm) in &ool {
            assert!(*t >= 0.0, "torque={t} at rpm={rpm} should be nonneg");
        }
    }

    #[test]
    fn ool_torque_bounded_by_max() {
        let max_torque = 400.0;
        let ool = optimal_operating_line(max_torque, 6000.0, 10);
        for (t, _) in &ool {
            assert!(*t <= max_torque + 1e-10, "torque={t} exceeds max");
        }
    }

    // ── standalone helpers ────────────────────────────────────────────────

    #[test]
    fn c_rate_basic() {
        assert!((c_rate(30.0, 60.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn c_rate_zero_capacity() {
        assert!(c_rate(10.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn braking_energy_recovered_simple() {
        // 1000 kg, 20 m/s -> 10 m/s, eta=1.0
        let e = braking_energy_recovered(1000.0, 20.0, 10.0, 1.0);
        let expected = 0.5 * 1000.0 * (400.0 - 100.0);
        assert!((e - expected).abs() < 1e-6, "e={e} expected={expected}");
    }

    #[test]
    fn braking_energy_recovered_zero_speed_diff() {
        let e = braking_energy_recovered(1000.0, 20.0, 20.0, 0.9);
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn braking_energy_recovered_efficiency_scales() {
        let e1 = braking_energy_recovered(1000.0, 20.0, 0.0, 1.0);
        let e2 = braking_energy_recovered(1000.0, 20.0, 0.0, 0.5);
        assert!((e1 - 2.0 * e2).abs() < 1e-6, "e1={e1} e2={e2}");
    }
}
