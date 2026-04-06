#![allow(clippy::if_same_then_else)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle thermal management system.
//!
//! Provides lumped-parameter thermal models for heat exchangers, cooling
//! circuits, battery packs, cabin climate control, and radiators.

// ── ThermalNode ───────────────────────────────────────────────────────────────

/// A lumped thermal mass node with heat input, output, and first-order dynamics.
#[derive(Debug, Clone)]
pub struct ThermalNode {
    /// Current temperature of the node (°C or K — consistent with inputs).
    pub temperature: f64,
    /// Thermal mass C_p (J/K).
    pub thermal_mass: f64,
    /// Heat flowing into the node per second (W).
    pub heat_input: f64,
    /// Heat flowing out of the node per second (W).
    pub heat_output: f64,
}

impl ThermalNode {
    /// Create a new thermal node.
    ///
    /// * `temperature` — initial temperature.
    /// * `thermal_mass` — thermal capacitance in J/K.
    pub fn new(temperature: f64, thermal_mass: f64) -> Self {
        Self {
            temperature,
            thermal_mass,
            heat_input: 0.0,
            heat_output: 0.0,
        }
    }

    /// Advance the node temperature by `dt` seconds.
    ///
    /// Applies first-order Euler integration:
    /// `T += (Q_in - Q_out) / C_p * dt`.
    ///
    /// Returns the updated temperature.
    pub fn step(&mut self, dt: f64) -> f64 {
        if self.thermal_mass > 0.0 {
            let net_heat = self.heat_input - self.heat_output;
            self.temperature += net_heat / self.thermal_mass * dt;
        }
        self.temperature
    }
}

// ── HeatExchanger ─────────────────────────────────────────────────────────────

/// Effectiveness-NTU heat exchanger model.
#[derive(Debug, Clone)]
pub struct HeatExchanger {
    /// Heat exchanger effectiveness ε ∈ \[0, 1\].
    pub effectiveness: f64,
    /// Minimum heat capacity rate C_min = min(ṁ_h·cp_h, ṁ_c·cp_c) (W/K).
    pub c_min: f64,
    /// Maximum heat capacity rate C_max (W/K).
    pub c_max: f64,
    /// Inlet hot-side temperature (°C or K).
    pub hot_temp: f64,
    /// Inlet cold-side temperature (°C or K).
    pub cold_temp: f64,
}

impl HeatExchanger {
    /// Create a new heat exchanger.
    pub fn new(effectiveness: f64, c_min: f64, c_max: f64, hot_temp: f64, cold_temp: f64) -> Self {
        Self {
            effectiveness: effectiveness.clamp(0.0, 1.0),
            c_min,
            c_max,
            hot_temp,
            cold_temp,
        }
    }

    /// Compute heat transfer rate Q (W) using the effectiveness-NTU method:
    ///
    /// `Q = ε · C_min · (T_hot - T_cold)`.
    pub fn heat_transfer(&self) -> f64 {
        self.effectiveness * self.c_min * (self.hot_temp - self.cold_temp)
    }

    /// Compute effectiveness from the number of transfer units (NTU) and
    /// capacity ratio `C_r = C_min / C_max` for a counter-flow arrangement.
    ///
    /// For `C_r < 1`:
    /// `ε = (1 - exp(-NTU*(1-Cr))) / (1 - Cr*exp(-NTU*(1-Cr)))`
    ///
    /// For `C_r = 1`:
    /// `ε = NTU / (1 + NTU)`
    pub fn ntu_method(&mut self, ua: f64) -> f64 {
        let ntu = if self.c_min > 0.0 {
            ua / self.c_min
        } else {
            0.0
        };
        let cr = if self.c_max > 0.0 {
            self.c_min / self.c_max
        } else {
            1.0
        };

        let eps = if (cr - 1.0).abs() < 1e-6 {
            ntu / (1.0 + ntu)
        } else {
            let exp_term = (-ntu * (1.0 - cr)).exp();
            (1.0 - exp_term) / (1.0 - cr * exp_term)
        };
        self.effectiveness = eps.clamp(0.0, 1.0);
        self.effectiveness
    }
}

// ── CoolingCircuit ────────────────────────────────────────────────────────────

/// A closed cooling loop with multiple thermal nodes and a circulating pump.
#[derive(Debug, Clone)]
pub struct CoolingCircuit {
    /// Thermal nodes in the circuit (e.g., engine, radiator, heater core).
    pub nodes: Vec<ThermalNode>,
    /// Pump mechanical power (W).
    pub pump_power: f64,
    /// Coolant volumetric flow rate (L/s).
    pub flow_rate: f64,
}

impl CoolingCircuit {
    /// Create a cooling circuit with the given nodes and pump parameters.
    pub fn new(nodes: Vec<ThermalNode>, pump_power: f64, flow_rate: f64) -> Self {
        Self {
            nodes,
            pump_power,
            flow_rate,
        }
    }

    /// Advance all nodes by `dt` seconds.
    ///
    /// Pump heat (pump_power) is added as heat input to the first node.
    pub fn step(&mut self, dt: f64) {
        if let Some(first) = self.nodes.first_mut() {
            first.heat_input += self.pump_power;
        }
        for node in &mut self.nodes {
            node.step(dt);
        }
        // Restore pump_power contribution (avoid accumulation across steps)
        if let Some(first) = self.nodes.first_mut() {
            first.heat_input -= self.pump_power;
        }
    }

    /// Return the mean coolant temperature across all nodes (°C or K).
    pub fn mean_temperature(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        self.nodes.iter().map(|n| n.temperature).sum::<f64>() / self.nodes.len() as f64
    }
}

// ── BatteryThermal ────────────────────────────────────────────────────────────

/// Lumped thermal model for a battery pack.
#[derive(Debug, Clone)]
pub struct BatteryThermal {
    /// Number of cells in the pack.
    pub n_cells: usize,
    /// Per-cell thermal mass (J/K).
    pub cell_thermal_mass: f64,
    /// Per-cell DC internal resistance (Ω).
    pub internal_resistance: f64,
    /// Current pack temperature (°C or K).
    pub temperature: f64,
}

impl BatteryThermal {
    /// Create a new battery thermal model.
    pub fn new(n_cells: usize, cell_thermal_mass: f64, internal_resistance: f64) -> Self {
        Self {
            n_cells,
            cell_thermal_mass,
            internal_resistance,
            temperature: 25.0,
        }
    }

    /// Compute total pack heat generation from ohmic losses (W).
    ///
    /// `Q = n_cells · I² · R_cell`
    pub fn generate_heat(&self, current: f64) -> f64 {
        self.n_cells as f64 * current * current * self.internal_resistance
    }

    /// Estimate required cooling power (W) to keep pack at or below `t_max`.
    ///
    /// A simple proportional model: if `T > T_max`, return the excess times
    /// a cooling coefficient (100 W/K), otherwise 0.
    pub fn cooling_requirement(&self, t_max: f64) -> f64 {
        let excess = self.temperature - t_max;
        if excess > 0.0 { excess * 100.0 } else { 0.0 }
    }

    /// Total pack thermal mass (J/K).
    pub fn total_thermal_mass(&self) -> f64 {
        self.n_cells as f64 * self.cell_thermal_mass
    }
}

// ── ClimateControl ────────────────────────────────────────────────────────────

/// Vehicle cabin climate control with HVAC system.
#[derive(Debug, Clone)]
pub struct ClimateControl {
    /// Cabin thermal mass (J/K).
    pub cabin_thermal_mass: f64,
    /// HVAC system power limit (W).
    pub hvac_power: f64,
    /// Target cabin temperature setpoint (°C).
    pub target_temp: f64,
    /// Coefficient of performance for heating.
    pub cop_heating: f64,
    /// Coefficient of performance for cooling.
    pub cop_cooling: f64,
    /// Current cabin temperature (°C).
    pub cabin_temp: f64,
}

impl ClimateControl {
    /// Create a new climate control model.
    pub fn new(
        cabin_thermal_mass: f64,
        hvac_power: f64,
        target_temp: f64,
        cop_heating: f64,
        cop_cooling: f64,
    ) -> Self {
        Self {
            cabin_thermal_mass,
            hvac_power,
            target_temp,
            cop_heating,
            cop_cooling,
            cabin_temp: target_temp,
        }
    }

    /// Advance the cabin temperature model by `dt` seconds.
    ///
    /// Heat exchange with the ambient through a simple conductance, plus HVAC
    /// active conditioning.
    ///
    /// Returns the updated cabin temperature.
    pub fn step(&mut self, dt: f64, ambient_temp: f64) -> f64 {
        // Simple infiltration / conduction (UA_cabin ~ 50 W/K)
        let ua_cabin = 50.0;
        let q_infiltration = ua_cabin * (ambient_temp - self.cabin_temp);

        // HVAC action
        let error = self.target_temp - self.cabin_temp;
        let q_hvac = if error > 0.0 {
            // Heating needed
            (self.hvac_power * self.cop_heating).min(error * self.cabin_thermal_mass / dt)
        } else if error < 0.0 {
            // Cooling needed: negative heat into cabin
            -(self.hvac_power * self.cop_cooling).min((-error) * self.cabin_thermal_mass / dt)
        } else {
            0.0
        };

        let q_net = q_infiltration + q_hvac;
        self.cabin_temp += q_net / self.cabin_thermal_mass * dt;
        self.cabin_temp
    }

    /// Estimated HVAC power draw from the vehicle electrical system (W).
    pub fn electrical_power_draw(&self) -> f64 {
        let error = self.target_temp - self.cabin_temp;
        if error.abs() < 0.1 {
            0.0
        } else if error > 0.0 {
            self.hvac_power
        } else {
            self.hvac_power
        }
    }
}

// ── RadiatorModel ─────────────────────────────────────────────────────────────

/// Automotive radiator heat rejection model.
#[derive(Debug, Clone)]
pub struct RadiatorModel {
    /// Frontal area of the radiator (m²).
    pub area: f64,
    /// Convective heat transfer coefficient at reference airflow (W/m²K).
    pub h_conv: f64,
    /// Overall fin efficiency η_f ∈ \[0, 1\].
    pub fin_efficiency: f64,
}

impl RadiatorModel {
    /// Create a new radiator model.
    pub fn new(area: f64, h_conv: f64, fin_efficiency: f64) -> Self {
        Self {
            area,
            h_conv,
            fin_efficiency: fin_efficiency.clamp(0.0, 1.0),
        }
    }

    /// Compute heat rejection (W) to the airstream.
    ///
    /// `Q = η_f · h · A · (T_coolant - T_air) · (1 + k_airflow * airflow)`
    ///
    /// where `k_airflow = 0.05 (s/m)` models increased convection at higher speeds.
    pub fn heat_rejection(&self, t_coolant: f64, t_air: f64, airflow: f64) -> f64 {
        let k_airflow = 0.05_f64;
        let delta_t = t_coolant - t_air;
        let airflow_factor = 1.0 + k_airflow * airflow.max(0.0);
        self.fin_efficiency * self.h_conv * self.area * delta_t * airflow_factor
    }

    /// Required coolant temperature to achieve a target heat rejection (W).
    pub fn required_coolant_temp(&self, q_target: f64, t_air: f64, airflow: f64) -> f64 {
        let k_airflow = 0.05_f64;
        let airflow_factor = 1.0 + k_airflow * airflow.max(0.0);
        let denom = self.fin_efficiency * self.h_conv * self.area * airflow_factor;
        if denom > 1e-12 {
            t_air + q_target / denom
        } else {
            t_air
        }
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Compute the convective heat transfer coefficient using the Wilson
/// (Dittus-Boelter style) correlation:
///
/// `h = 0.023 · Re^0.8 · Pr^0.4 · k / D`
///
/// * `re` — Reynolds number
/// * `pr` — Prandtl number
/// * `d` — hydraulic diameter (m)
/// * `k` — fluid thermal conductivity (W/m·K)
pub fn wilson_heat_transfer_coeff(re: f64, pr: f64, d: f64, k: f64) -> f64 {
    if d > 0.0 {
        0.023 * re.powf(0.8) * pr.powf(0.4) * k / d
    } else {
        0.0
    }
}

/// Compute convective heat transfer using Newton's law of cooling:
///
/// `Q = h · A · ΔT` (W).
pub fn newton_cooling(h: f64, a: f64, delta_t: f64) -> f64 {
    h * a * delta_t
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ThermalNode ───────────────────────────────────────────────────────

    #[test]
    fn thermal_node_step_heats_up() {
        let mut node = ThermalNode::new(20.0, 1000.0);
        node.heat_input = 1000.0; // 1000 W in
        node.heat_output = 0.0;
        let t = node.step(1.0);
        assert!((t - 21.0).abs() < 1e-6, "T={t}");
    }

    #[test]
    fn thermal_node_step_cools_down() {
        let mut node = ThermalNode::new(80.0, 2000.0);
        node.heat_output = 2000.0;
        let t = node.step(1.0);
        assert!((t - 79.0).abs() < 1e-6, "T={t}");
    }

    #[test]
    fn thermal_node_zero_thermal_mass_unchanged() {
        let mut node = ThermalNode::new(25.0, 0.0);
        node.heat_input = 1000.0;
        let t = node.step(1.0);
        assert!((t - 25.0).abs() < 1e-10);
    }

    #[test]
    fn thermal_node_equilibrium_no_change() {
        let mut node = ThermalNode::new(60.0, 500.0);
        node.heat_input = 200.0;
        node.heat_output = 200.0;
        let t = node.step(1.0);
        assert!((t - 60.0).abs() < 1e-10);
    }

    #[test]
    fn thermal_node_small_dt() {
        let mut node = ThermalNode::new(20.0, 1000.0);
        node.heat_input = 1000.0;
        let t = node.step(0.001);
        assert!((t - 20.001).abs() < 1e-6, "T={t}");
    }

    // ── HeatExchanger ─────────────────────────────────────────────────────

    #[test]
    fn heat_exchanger_basic_heat_transfer() {
        let hx = HeatExchanger::new(0.8, 500.0, 1000.0, 90.0, 20.0);
        let q = hx.heat_transfer();
        // Q = 0.8 * 500 * (90-20) = 28000
        assert!((q - 28000.0).abs() < 1e-6, "Q={q}");
    }

    #[test]
    fn heat_exchanger_zero_delta_t_no_transfer() {
        let hx = HeatExchanger::new(0.9, 500.0, 1000.0, 50.0, 50.0);
        assert!(hx.heat_transfer().abs() < 1e-10);
    }

    #[test]
    fn heat_exchanger_ntu_method_cr_zero() {
        let mut hx = HeatExchanger::new(0.5, 100.0, 1e12, 80.0, 20.0);
        // With C_r → 0: ε = 1 - exp(-NTU)
        let eps = hx.ntu_method(100.0); // NTU = 1.0
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((eps - expected).abs() < 1e-6, "eps={eps}");
    }

    #[test]
    fn heat_exchanger_ntu_method_cr_one() {
        let mut hx = HeatExchanger::new(0.5, 100.0, 100.0, 80.0, 20.0);
        // C_r = 1: ε = NTU / (1 + NTU)
        let eps = hx.ntu_method(200.0); // NTU = 2.0
        let expected = 2.0 / 3.0;
        assert!((eps - expected).abs() < 1e-6, "eps={eps}");
    }

    #[test]
    fn heat_exchanger_effectiveness_clamped() {
        let mut hx = HeatExchanger::new(0.5, 100.0, 100.0, 80.0, 20.0);
        let eps = hx.ntu_method(1e6);
        assert!((0.0..=1.0).contains(&eps));
    }

    #[test]
    fn heat_exchanger_ntu_higher_ua_higher_eps() {
        let mut hx = HeatExchanger::new(0.5, 100.0, 200.0, 80.0, 20.0);
        let eps1 = hx.ntu_method(100.0);
        let eps2 = hx.ntu_method(300.0);
        assert!(eps2 > eps1, "higher UA should give higher effectiveness");
    }

    // ── CoolingCircuit ────────────────────────────────────────────────────

    #[test]
    fn cooling_circuit_mean_temp() {
        let nodes = vec![
            ThermalNode::new(40.0, 1000.0),
            ThermalNode::new(60.0, 1000.0),
        ];
        let circuit = CoolingCircuit::new(nodes, 100.0, 1.0);
        assert!((circuit.mean_temperature() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn cooling_circuit_empty_mean_zero() {
        let circuit = CoolingCircuit::new(vec![], 0.0, 0.0);
        assert!(circuit.mean_temperature().abs() < 1e-10);
    }

    #[test]
    fn cooling_circuit_step_increases_first_node_temp() {
        let mut circuit = CoolingCircuit::new(
            vec![ThermalNode::new(20.0, 1000.0)],
            1000.0, // 1000 W pump
            1.0,
        );
        circuit.step(1.0);
        // First node should have received 1000 W for 1 s → +1 K
        assert!(circuit.nodes[0].temperature > 20.0);
    }

    #[test]
    fn cooling_circuit_pump_power_not_accumulated() {
        let mut circuit = CoolingCircuit::new(vec![ThermalNode::new(20.0, 1000.0)], 1000.0, 1.0);
        circuit.step(1.0);
        // After step, heat_input should be restored to original value (0.0)
        assert!((circuit.nodes[0].heat_input - 0.0).abs() < 1e-10);
    }

    // ── BatteryThermal ────────────────────────────────────────────────────

    #[test]
    fn battery_heat_generation_ohmic() {
        let bat = BatteryThermal::new(100, 500.0, 0.005);
        // Q = 100 * 10^2 * 0.005 = 50 W
        let q = bat.generate_heat(10.0);
        assert!((q - 50.0).abs() < 1e-6, "Q={q}");
    }

    #[test]
    fn battery_heat_generation_zero_current() {
        let bat = BatteryThermal::new(100, 500.0, 0.005);
        assert!(bat.generate_heat(0.0).abs() < 1e-10);
    }

    #[test]
    fn battery_cooling_requirement_above_max() {
        let mut bat = BatteryThermal::new(96, 500.0, 0.005);
        bat.temperature = 45.0;
        let q = bat.cooling_requirement(40.0);
        assert!((q - 500.0).abs() < 1e-6, "Q={q}");
    }

    #[test]
    fn battery_cooling_requirement_below_max() {
        let bat = BatteryThermal::new(96, 500.0, 0.005);
        // temperature = 25.0, t_max = 40 → no cooling needed
        let q = bat.cooling_requirement(40.0);
        assert!(q.abs() < 1e-10);
    }

    #[test]
    fn battery_total_thermal_mass() {
        let bat = BatteryThermal::new(100, 500.0, 0.005);
        assert!((bat.total_thermal_mass() - 50000.0).abs() < 1e-6);
    }

    #[test]
    fn battery_heat_quadratic_current() {
        let bat = BatteryThermal::new(1, 500.0, 1.0);
        let q1 = bat.generate_heat(1.0);
        let q2 = bat.generate_heat(2.0);
        assert!((q2 - 4.0 * q1).abs() < 1e-10);
    }

    // ── ClimateControl ────────────────────────────────────────────────────

    #[test]
    fn climate_control_at_target_no_hvac() {
        let mut cc = ClimateControl::new(50000.0, 5000.0, 22.0, 3.5, 3.0);
        cc.cabin_temp = 22.0;
        // Step with ambient = 22°C — no hvac needed, cabin should stay ~22
        let t = cc.step(1.0, 22.0);
        assert!((t - 22.0).abs() < 0.1, "T={t}");
    }

    #[test]
    fn climate_control_heating_cold_ambient() {
        let mut cc = ClimateControl::new(50000.0, 5000.0, 22.0, 3.5, 3.0);
        cc.cabin_temp = 5.0;
        let t_initial = cc.cabin_temp;
        let t = cc.step(60.0, -10.0);
        // With heating, cabin should be warmer than initial
        assert!(t >= t_initial, "T={t} should be >= {t_initial}");
    }

    #[test]
    fn climate_control_cooling_hot_ambient() {
        let mut cc = ClimateControl::new(50000.0, 5000.0, 22.0, 3.5, 3.0);
        cc.cabin_temp = 40.0;
        let t = cc.step(1.0, 35.0);
        // Cabin should drop toward target
        assert!(t < 40.0, "T={t} should drop from 40");
    }

    #[test]
    fn climate_control_electrical_power_nonzero_when_off_target() {
        let mut cc = ClimateControl::new(50000.0, 5000.0, 22.0, 3.5, 3.0);
        cc.cabin_temp = 35.0;
        assert!(cc.electrical_power_draw() > 0.0);
    }

    #[test]
    fn climate_control_electrical_power_zero_at_target() {
        let mut cc = ClimateControl::new(50000.0, 5000.0, 22.0, 3.5, 3.0);
        cc.cabin_temp = 22.0;
        assert!(cc.electrical_power_draw().abs() < 1e-10);
    }

    // ── RadiatorModel ─────────────────────────────────────────────────────

    #[test]
    fn radiator_heat_rejection_basic() {
        let rad = RadiatorModel::new(0.5, 200.0, 0.9);
        // Q = 0.9 * 200 * 0.5 * (90-20) * (1+0.05*10) = 90 * 70 * 1.5 = 9450
        let q = rad.heat_rejection(90.0, 20.0, 10.0);
        assert!((q - 9450.0).abs() < 1e-6, "Q={q}");
    }

    #[test]
    fn radiator_zero_delta_t_zero_rejection() {
        let rad = RadiatorModel::new(0.5, 200.0, 0.9);
        assert!(rad.heat_rejection(50.0, 50.0, 5.0).abs() < 1e-10);
    }

    #[test]
    fn radiator_higher_airflow_more_rejection() {
        let rad = RadiatorModel::new(0.5, 200.0, 0.9);
        let q1 = rad.heat_rejection(90.0, 20.0, 0.0);
        let q2 = rad.heat_rejection(90.0, 20.0, 30.0);
        assert!(q2 > q1, "more airflow should increase rejection");
    }

    #[test]
    fn radiator_negative_airflow_treated_as_zero() {
        let rad = RadiatorModel::new(0.5, 200.0, 0.9);
        let q_neg = rad.heat_rejection(90.0, 20.0, -5.0);
        let q_zero = rad.heat_rejection(90.0, 20.0, 0.0);
        assert!((q_neg - q_zero).abs() < 1e-10);
    }

    #[test]
    fn radiator_required_coolant_temp_consistent() {
        let rad = RadiatorModel::new(0.5, 200.0, 0.9);
        let q_target = 5000.0;
        let t_coolant = rad.required_coolant_temp(q_target, 20.0, 0.0);
        let q_actual = rad.heat_rejection(t_coolant, 20.0, 0.0);
        assert!((q_actual - q_target).abs() < 1e-3, "q={q_actual}");
    }

    // ── wilson_heat_transfer_coeff ────────────────────────────────────────

    #[test]
    fn wilson_coeff_turbulent_water() {
        // Re=10000, Pr=7 (water), D=0.01 m, k=0.6 W/mK
        let h = wilson_heat_transfer_coeff(10000.0, 7.0, 0.01, 0.6);
        assert!(h > 0.0, "h={h}");
        // Expected ~ 0.023 * 10000^0.8 * 7^0.4 * 0.6 / 0.01 ≈ several kW/m²K
        assert!(
            h > 1000.0,
            "h={h} should be substantial for turbulent water"
        );
    }

    #[test]
    fn wilson_coeff_zero_diameter_zero() {
        let h = wilson_heat_transfer_coeff(10000.0, 7.0, 0.0, 0.6);
        assert!(h.abs() < 1e-10);
    }

    #[test]
    fn wilson_coeff_increases_with_reynolds() {
        let h1 = wilson_heat_transfer_coeff(10000.0, 7.0, 0.01, 0.6);
        let h2 = wilson_heat_transfer_coeff(20000.0, 7.0, 0.01, 0.6);
        assert!(h2 > h1);
    }

    #[test]
    fn wilson_coeff_air() {
        // Re=50000, Pr=0.71 (air), D=0.02 m, k=0.026 W/mK
        let h = wilson_heat_transfer_coeff(50000.0, 0.71, 0.02, 0.026);
        assert!(h > 0.0, "h={h}");
    }

    // ── newton_cooling ────────────────────────────────────────────────────

    #[test]
    fn newton_cooling_basic() {
        // Q = 10 * 2 * 30 = 600 W
        let q = newton_cooling(10.0, 2.0, 30.0);
        assert!((q - 600.0).abs() < 1e-10);
    }

    #[test]
    fn newton_cooling_zero_delta_t() {
        assert!(newton_cooling(100.0, 1.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn newton_cooling_negative_delta_t() {
        // Cold surface — negative Q means heat absorbed
        let q = newton_cooling(10.0, 1.0, -5.0);
        assert!((q + 50.0).abs() < 1e-10);
    }

    #[test]
    fn newton_cooling_proportional_to_area() {
        let q1 = newton_cooling(10.0, 1.0, 20.0);
        let q2 = newton_cooling(10.0, 3.0, 20.0);
        assert!((q2 - 3.0 * q1).abs() < 1e-10);
    }
}
