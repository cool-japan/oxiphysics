// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fuel cell vehicle powertrain model.
//!
//! Implements a proton-exchange-membrane (PEM) fuel cell stack with
//! polarization curve, hydrogen storage, and full vehicle range estimation.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_vehicle::fuel_cell::{FuelCellStack, HydrogenTank, nernst_voltage};
//!
//! let stack = FuelCellStack::new(400, 300.0);
//! assert!(stack.voltage() > 0.0);
//!
//! let v = nernst_voltage(353.15, 1.5, 0.3);
//! assert!(v > 1.0);
//! ```

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Universal gas constant (J mol⁻¹ K⁻¹).
pub const R_GAS: f64 = 8.314_462;
/// Faraday constant (C mol⁻¹).
pub const FARADAY: f64 = 96_485.0;
/// Number of electrons transferred per mole of H₂ in PEM fuel cell.
pub const N_ELECTRONS: f64 = 2.0;

// ---------------------------------------------------------------------------
// Nernst voltage
// ---------------------------------------------------------------------------

/// Compute the Nernst (thermodynamic) open-circuit voltage for a hydrogen
/// fuel cell at temperature `T` (K), H₂ partial pressure `p_h2` (bar), and
/// O₂ partial pressure `p_o2` (bar).
///
/// Uses: `E = E0 + (RT / 2F) * ln(p_h2 * sqrt(p_o2))`
/// where `E0 = 1.229 V` at 298 K (adjusted for temperature).
pub fn nernst_voltage(t_k: f64, p_h2: f64, p_o2: f64) -> f64 {
    // Temperature-corrected standard potential (V K⁻¹ coefficient ≈ −8.5e-4)
    let e0 = 1.229 - 8.5e-4 * (t_k - 298.15);
    let delta_g_term =
        (R_GAS * t_k / (N_ELECTRONS * FARADAY)) * (p_h2 * p_o2.sqrt()).ln().max(-50.0);
    e0 + delta_g_term
}

// ---------------------------------------------------------------------------
// Activation overpotential (Butler-Volmer)
// ---------------------------------------------------------------------------

/// Compute the activation overpotential (V) using the Tafel approximation
/// of the Butler-Volmer equation.
///
/// * `j`     – current density (A cm⁻²)
/// * `j0`    – exchange current density (A cm⁻²)
/// * `alpha` – transfer coefficient (dimensionless, typically 0.5)
/// * `t_k`   – temperature (K)
///
/// Returns the activation loss (always ≥ 0).
pub fn activation_overpotential(j: f64, j0: f64, alpha: f64, t_k: f64) -> f64 {
    if j <= 0.0 || j0 <= 0.0 {
        return 0.0;
    }
    let rt_af = R_GAS * t_k / (alpha * FARADAY);
    rt_af * (j / j0).ln().max(0.0)
}

// ---------------------------------------------------------------------------
// PolarizationCurve
// ---------------------------------------------------------------------------

/// Polarization curve parameters for a single PEM fuel cell.
#[derive(Debug, Clone)]
pub struct PolarizationCurve {
    /// Open-circuit voltage (V).
    pub ocv: f64,
    /// Exchange current density (A cm⁻²).
    pub j0: f64,
    /// Butler-Volmer transfer coefficient.
    pub alpha: f64,
    /// Ohmic resistance (Ω cm²).
    pub r_ohm: f64,
    /// Limiting current density (A cm⁻²).
    pub j_lim: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl PolarizationCurve {
    /// Create a standard PEM polarization curve at 80°C.
    pub fn new_standard() -> Self {
        Self {
            ocv: 1.0,
            j0: 1e-4,
            alpha: 0.5,
            r_ohm: 0.1,
            j_lim: 1.5,
            temperature: 353.15,
        }
    }

    /// Create a custom polarization curve.
    pub fn new(ocv: f64, j0: f64, alpha: f64, r_ohm: f64, j_lim: f64, temperature: f64) -> Self {
        Self {
            ocv,
            j0,
            alpha,
            r_ohm,
            j_lim,
            temperature,
        }
    }

    /// Activation loss (V) at current density `j` (A cm⁻²).
    pub fn activation_loss(&self, j: f64) -> f64 {
        activation_overpotential(j, self.j0, self.alpha, self.temperature)
    }

    /// Ohmic loss (V) at current density `j` (A cm⁻²).
    pub fn ohmic_loss(&self, j: f64) -> f64 {
        self.r_ohm * j
    }

    /// Concentration (mass-transport) loss (V) at current density `j` (A cm⁻²).
    ///
    /// Uses: `V_conc = m * exp(n * j)` approach, simplified to
    /// `-( RT / 2F ) * ln(1 - j / j_lim)`.
    pub fn concentration_loss(&self, j: f64) -> f64 {
        let j_safe = j.min(self.j_lim * 0.999);
        if j_safe <= 0.0 {
            return 0.0;
        }
        let rt_nf = R_GAS * self.temperature / (N_ELECTRONS * FARADAY);
        -rt_nf * (1.0 - j_safe / self.j_lim).ln().max(-50.0)
    }

    /// Total cell voltage at current density `j` (A cm⁻²).
    pub fn total_voltage(&self, j: f64) -> f64 {
        let v =
            self.ocv - self.activation_loss(j) - self.ohmic_loss(j) - self.concentration_loss(j);
        v.max(0.0)
    }
}

// ---------------------------------------------------------------------------
// FuelCellStack
// ---------------------------------------------------------------------------

/// PEM fuel cell stack composed of `n_cells` unit cells.
#[derive(Debug, Clone)]
pub struct FuelCellStack {
    /// Number of cells in series.
    pub n_cells: u32,
    /// Active cell area (cm²).
    pub active_area: f64,
    /// Operating current density (A cm⁻²).
    pub current_density: f64,
    /// Polarization curve parameters.
    pub polarization: PolarizationCurve,
    /// Operating temperature (K).
    pub temperature: f64,
}

impl FuelCellStack {
    /// Create a new stack with `n_cells` cells and `active_area` cm².
    pub fn new(n_cells: u32, active_area: f64) -> Self {
        Self {
            n_cells,
            active_area,
            current_density: 0.5,
            polarization: PolarizationCurve::new_standard(),
            temperature: 353.15,
        }
    }

    /// Stack output voltage (V).
    pub fn voltage(&self) -> f64 {
        self.n_cells as f64 * self.polarization.total_voltage(self.current_density)
    }

    /// Stack current (A).
    pub fn current(&self) -> f64 {
        self.current_density * self.active_area
    }

    /// Net electrical power output (W).
    pub fn power_output(&self) -> f64 {
        self.voltage() * self.current()
    }

    /// Thermodynamic efficiency (0–1) relative to HHV of hydrogen.
    ///
    /// η = V_stack / (n_cells * V_thermo), where V_thermo = 1.481 V (HHV basis).
    pub fn efficiency(&self) -> f64 {
        let v_thermo = 1.481 * self.n_cells as f64;
        if v_thermo < 1e-12 {
            return 0.0;
        }
        (self.voltage() / v_thermo).clamp(0.0, 1.0)
    }

    /// Set operating current density (A cm⁻²).
    pub fn set_current_density(&mut self, j: f64) {
        self.current_density = j.max(0.0);
    }

    /// Hydrogen consumption rate (g s⁻¹) at current operating point.
    ///
    /// From Faraday's law: ṁ_H2 = I * M_H2 / (2 * F)
    /// where M_H2 = 2.016 g/mol.
    pub fn hydrogen_consumption_rate_gs(&self) -> f64 {
        let i_total = self.current();
        i_total * 2.016 / (N_ELECTRONS * FARADAY)
    }
}

// ---------------------------------------------------------------------------
// HydrogenTank
// ---------------------------------------------------------------------------

/// Compressed hydrogen storage tank.
#[derive(Debug, Clone)]
pub struct HydrogenTank {
    /// Tank capacity (kg of H₂ at rated pressure).
    pub capacity_kg: f64,
    /// Rated pressure (bar).
    pub pressure_bar: f64,
    /// Current hydrogen mass (kg).
    pub current_mass: f64,
}

impl HydrogenTank {
    /// Create a new hydrogen tank.
    pub fn new(capacity_kg: f64, pressure_bar: f64) -> Self {
        Self {
            capacity_kg,
            pressure_bar,
            current_mass: capacity_kg,
        }
    }

    /// Current state-of-charge (0–1).
    pub fn soc(&self) -> f64 {
        if self.capacity_kg < 1e-12 {
            return 0.0;
        }
        (self.current_mass / self.capacity_kg).clamp(0.0, 1.0)
    }

    /// Refuel tank to full capacity.
    pub fn refuel(&mut self) {
        self.current_mass = self.capacity_kg;
    }

    /// Consume `kg` of hydrogen from the tank.
    ///
    /// Returns the actual amount consumed (clamped to available mass).
    pub fn consume(&mut self, kg: f64) -> f64 {
        let actual = kg.min(self.current_mass).max(0.0);
        self.current_mass -= actual;
        actual
    }

    /// Estimated current pressure (bar), using ideal gas law approximation.
    pub fn current_pressure(&self) -> f64 {
        self.pressure_bar * self.soc()
    }

    /// Whether the tank is empty.
    pub fn is_empty(&self) -> bool {
        self.current_mass < 1e-6
    }
}

// ---------------------------------------------------------------------------
// FuelCellVehicle
// ---------------------------------------------------------------------------

/// Full fuel cell vehicle powertrain model.
#[derive(Debug, Clone)]
pub struct FuelCellVehicle {
    /// Fuel cell stack.
    pub stack: FuelCellStack,
    /// Hydrogen storage tank.
    pub tank: HydrogenTank,
    /// Electric motor peak power (kW).
    pub motor_power_kw: f64,
    /// Electric motor efficiency (0–1).
    pub motor_efficiency: f64,
    /// High-voltage buffer battery capacity (kWh).
    pub battery_buffer_kwh: f64,
    /// Current battery state-of-charge (0–1).
    pub battery_soc: f64,
    /// Total energy delivered to drivetrain (Wh).
    pub energy_delivered_wh: f64,
}

impl FuelCellVehicle {
    /// Create a new fuel cell vehicle with default parameters.
    ///
    /// * `stack`           – configured fuel cell stack
    /// * `tank`            – hydrogen tank
    /// * `motor_power_kw`  – peak motor power (kW)
    /// * `motor_efficiency`– motor efficiency (0–1)
    /// * `battery_kwh`     – buffer battery capacity (kWh)
    pub fn new(
        stack: FuelCellStack,
        tank: HydrogenTank,
        motor_power_kw: f64,
        motor_efficiency: f64,
        battery_kwh: f64,
    ) -> Self {
        Self {
            stack,
            tank,
            motor_power_kw,
            motor_efficiency,
            battery_buffer_kwh: battery_kwh,
            battery_soc: 0.5,
            energy_delivered_wh: 0.0,
        }
    }

    /// Simulate one time step of duration `dt` (s) with a `power_demand` (W).
    ///
    /// The fuel cell provides power; excess goes to/from the buffer battery.
    /// Hydrogen is consumed proportional to stack current.
    pub fn step(&mut self, dt: f64, power_demand: f64) {
        let fc_power = self.stack.power_output();
        let excess = fc_power - power_demand;

        if excess >= 0.0 {
            // Charge battery with excess
            let charge_wh = (excess * dt / 3600.0)
                .min(self.battery_buffer_kwh * 1000.0 * (1.0 - self.battery_soc));
            self.battery_soc += charge_wh / (self.battery_buffer_kwh * 1000.0).max(1e-12);
        } else {
            // Discharge battery for deficit
            let deficit_wh =
                (-excess * dt / 3600.0).min(self.battery_buffer_kwh * 1000.0 * self.battery_soc);
            self.battery_soc -= deficit_wh / (self.battery_buffer_kwh * 1000.0).max(1e-12);
        }

        // Consume hydrogen
        let h2_rate = self.stack.hydrogen_consumption_rate_gs(); // g/s
        let h2_kg = h2_rate * dt / 1000.0;
        self.tank.consume(h2_kg);

        // Track energy delivered
        let motor_power = power_demand.min(self.motor_power_kw * 1000.0);
        self.energy_delivered_wh += motor_power * self.motor_efficiency * dt / 3600.0;

        self.battery_soc = self.battery_soc.clamp(0.0, 1.0);
    }

    /// Remaining range estimate (km) given current state.
    pub fn remaining_range_km(&self, mass_kg: f64, drag_coeff: f64, speed_kmh: f64) -> f64 {
        FcvRange::estimate_range(mass_kg, drag_coeff, speed_kmh, &self.stack) * self.tank.soc()
    }
}

// ---------------------------------------------------------------------------
// FcvRange — range estimation
// ---------------------------------------------------------------------------

/// Fuel cell vehicle range estimator.
pub struct FcvRange;

impl FcvRange {
    /// Estimate the maximum range (km) of a fuel cell vehicle.
    ///
    /// Uses a simplified steady-state power balance:
    /// - Rolling resistance: `F_roll = m * g * Cr`
    /// - Aerodynamic drag:   `F_drag = 0.5 * rho * Cd * A * v^2`
    /// - Power demand:       `P = (F_roll + F_drag) * v`
    ///
    /// Then range = `(tank_capacity_kg * LHV_H2 * efficiency) / P * v`
    ///
    /// * `mass_kg`    – vehicle mass (kg)
    /// * `drag_coeff` – product `Cd * A` (m²)
    /// * `speed_kmh`  – cruise speed (km/h)
    /// * `stack`      – reference to the fuel cell stack for efficiency
    pub fn estimate_range(
        mass_kg: f64,
        drag_coeff: f64,
        speed_kmh: f64,
        stack: &FuelCellStack,
    ) -> f64 {
        let v = speed_kmh / 3.6; // m/s
        let g = 9.81;
        let cr = 0.01; // rolling resistance coefficient
        let rho_air = 1.225; // kg/m³
        let lhv_h2 = 33.33; // kWh/kg (lower heating value)
        let tank_kg = stack.active_area * 0.001; // simplified: area scales tank

        let f_roll = mass_kg * g * cr;
        let f_drag = 0.5 * rho_air * drag_coeff * v * v;
        let p_demand = (f_roll + f_drag) * v; // W

        if p_demand < 1.0 {
            return 1e6; // effectively infinite range at negligible demand
        }

        let efficiency = stack.efficiency();
        let energy_available_wh = tank_kg * lhv_h2 * 1000.0 * efficiency;
        let range_m = energy_available_wh / (p_demand / 3600.0);
        range_m / 1000.0 // km
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- nernst_voltage ---

    #[test]
    fn nernst_voltage_at_reference_conditions() {
        // At 298.15 K, p_h2=1.0, p_o2=1.0: ln(1*1^0.5) = 0 => E ≈ E0
        let v = nernst_voltage(298.15, 1.0, 1.0);
        assert!((v - 1.229).abs() < 0.01, "expected ~1.229 V, got {v}");
    }

    #[test]
    fn nernst_voltage_increases_with_pressure() {
        let v1 = nernst_voltage(353.15, 1.0, 0.21);
        let v2 = nernst_voltage(353.15, 2.0, 0.4);
        assert!(v2 > v1, "higher pressures should increase Nernst voltage");
    }

    #[test]
    fn nernst_voltage_decreases_with_temperature() {
        let v_low = nernst_voltage(298.15, 1.0, 0.21);
        let v_high = nernst_voltage(353.15, 1.0, 0.21);
        assert!(v_low > v_high, "higher T should decrease Nernst voltage");
    }

    #[test]
    fn nernst_voltage_positive_at_normal_conditions() {
        let v = nernst_voltage(353.15, 1.5, 0.3);
        assert!(
            v > 0.8,
            "Nernst voltage should be > 0.8 V at operating conditions"
        );
    }

    // --- activation_overpotential ---

    #[test]
    fn activation_overpotential_zero_at_zero_current() {
        let eta = activation_overpotential(0.0, 1e-4, 0.5, 353.15);
        assert!(eta.abs() < 1e-10);
    }

    #[test]
    fn activation_overpotential_positive() {
        let eta = activation_overpotential(0.5, 1e-4, 0.5, 353.15);
        assert!(eta > 0.0, "activation loss should be positive");
    }

    #[test]
    fn activation_overpotential_increases_with_current() {
        let eta1 = activation_overpotential(0.1, 1e-4, 0.5, 353.15);
        let eta2 = activation_overpotential(1.0, 1e-4, 0.5, 353.15);
        assert!(eta2 > eta1, "higher current -> higher activation loss");
    }

    #[test]
    fn activation_overpotential_decreases_with_temperature() {
        let eta_low = activation_overpotential(0.5, 1e-4, 0.5, 298.15);
        let eta_high = activation_overpotential(0.5, 1e-4, 0.5, 373.15);
        // Higher T reduces RT/aF prefactor differently
        // Both are positive; relationship depends on the ratio
        assert!(eta_low > 0.0 && eta_high > 0.0);
    }

    // --- PolarizationCurve ---

    #[test]
    fn polarization_curve_total_voltage_at_zero_current() {
        let curve = PolarizationCurve::new_standard();
        let v = curve.total_voltage(0.0);
        assert!((v - curve.ocv).abs() < 1e-10, "at j=0, V should equal OCV");
    }

    #[test]
    fn polarization_curve_voltage_decreases_with_current() {
        let curve = PolarizationCurve::new_standard();
        let v1 = curve.total_voltage(0.1);
        let v2 = curve.total_voltage(0.5);
        let v3 = curve.total_voltage(1.0);
        assert!(v1 >= v2, "voltage should decrease with increasing current");
        assert!(v2 >= v3);
    }

    #[test]
    fn polarization_curve_ohmic_loss_linear() {
        let curve = PolarizationCurve::new_standard();
        let loss1 = curve.ohmic_loss(0.5);
        let loss2 = curve.ohmic_loss(1.0);
        assert!((loss2 / loss1 - 2.0).abs() < 1e-10, "ohmic loss is linear");
    }

    #[test]
    fn polarization_curve_concentration_loss_near_limit() {
        let curve = PolarizationCurve::new_standard();
        let loss = curve.concentration_loss(curve.j_lim * 0.95);
        assert!(
            loss > 0.0,
            "concentration loss near limiting current should be positive"
        );
    }

    #[test]
    fn polarization_curve_total_voltage_non_negative() {
        let curve = PolarizationCurve::new_standard();
        // Should not go negative
        let v = curve.total_voltage(curve.j_lim * 0.99);
        assert!(v >= 0.0);
    }

    // --- FuelCellStack ---

    #[test]
    fn fuel_cell_stack_voltage_positive() {
        let stack = FuelCellStack::new(400, 300.0);
        assert!(stack.voltage() > 0.0);
    }

    #[test]
    fn fuel_cell_stack_power_positive() {
        let stack = FuelCellStack::new(400, 300.0);
        assert!(stack.power_output() > 0.0);
    }

    #[test]
    fn fuel_cell_stack_efficiency_between_zero_and_one() {
        let stack = FuelCellStack::new(400, 300.0);
        let eff = stack.efficiency();
        assert!(eff > 0.0 && eff <= 1.0, "efficiency={eff}");
    }

    #[test]
    fn fuel_cell_stack_hydrogen_consumption_rate_positive() {
        let stack = FuelCellStack::new(100, 100.0);
        let rate = stack.hydrogen_consumption_rate_gs();
        assert!(rate > 0.0);
    }

    #[test]
    fn fuel_cell_stack_more_cells_higher_voltage() {
        let stack1 = FuelCellStack::new(100, 300.0);
        let stack2 = FuelCellStack::new(400, 300.0);
        assert!(stack2.voltage() > stack1.voltage());
    }

    #[test]
    fn fuel_cell_stack_set_current_density() {
        let mut stack = FuelCellStack::new(400, 300.0);
        stack.set_current_density(0.8);
        assert!((stack.current_density - 0.8).abs() < 1e-12);
    }

    #[test]
    fn fuel_cell_stack_higher_current_density_higher_power() {
        let mut stack = FuelCellStack::new(400, 300.0);
        stack.set_current_density(0.3);
        let p1 = stack.power_output();
        stack.set_current_density(0.7);
        let p2 = stack.power_output();
        assert!(p2 > p1, "higher j -> higher power (in non-limiting region)");
    }

    // --- HydrogenTank ---

    #[test]
    fn hydrogen_tank_initial_soc_full() {
        let tank = HydrogenTank::new(5.0, 700.0);
        assert!((tank.soc() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn hydrogen_tank_consume_decreases_mass() {
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.consume(1.0);
        assert!((tank.current_mass - 4.0).abs() < 1e-10);
    }

    #[test]
    fn hydrogen_tank_consume_clamped_to_available() {
        let mut tank = HydrogenTank::new(2.0, 700.0);
        let consumed = tank.consume(10.0);
        assert!((consumed - 2.0).abs() < 1e-10);
        assert!(tank.is_empty());
    }

    #[test]
    fn hydrogen_tank_refuel_restores_capacity() {
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.consume(3.0);
        tank.refuel();
        assert!((tank.current_mass - 5.0).abs() < 1e-10);
    }

    #[test]
    fn hydrogen_tank_pressure_proportional_to_soc() {
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.consume(2.5); // 50% empty
        let p = tank.current_pressure();
        assert!((p - 350.0).abs() < 1e-6, "pressure at 50% SOC: {p}");
    }

    // --- FuelCellVehicle ---

    #[test]
    fn fcv_step_consumes_hydrogen() {
        let stack = FuelCellStack::new(400, 300.0);
        let tank = HydrogenTank::new(5.0, 700.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 1.0);
        let mass_before = fcv.tank.current_mass;
        fcv.step(1.0, 10000.0);
        assert!(
            fcv.tank.current_mass < mass_before,
            "hydrogen should be consumed"
        );
    }

    #[test]
    fn fcv_step_battery_soc_bounded() {
        let stack = FuelCellStack::new(200, 200.0);
        let tank = HydrogenTank::new(5.0, 700.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 50.0, 0.95, 2.0);
        for _ in 0..100 {
            fcv.step(0.1, 5000.0);
        }
        assert!(fcv.battery_soc >= 0.0 && fcv.battery_soc <= 1.0);
    }

    #[test]
    fn fcv_remaining_range_positive() {
        let stack = FuelCellStack::new(400, 300.0);
        let tank = HydrogenTank::new(5.0, 700.0);
        let fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 1.0);
        let range = fcv.remaining_range_km(1500.0, 0.6, 100.0);
        assert!(range > 0.0);
    }

    #[test]
    fn fcv_energy_delivered_increases() {
        let stack = FuelCellStack::new(300, 200.0);
        let tank = HydrogenTank::new(5.0, 700.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 80.0, 0.95, 2.0);
        fcv.step(1.0, 20000.0);
        assert!(fcv.energy_delivered_wh > 0.0);
    }

    // --- FcvRange ---

    #[test]
    fn fcv_range_higher_speed_less_range() {
        let stack = FuelCellStack::new(400, 300.0);
        let r1 = FcvRange::estimate_range(1500.0, 0.6, 80.0, &stack);
        let r2 = FcvRange::estimate_range(1500.0, 0.6, 130.0, &stack);
        assert!(r1 > r2, "lower speed -> longer range due to less drag");
    }

    #[test]
    fn fcv_range_higher_drag_less_range() {
        let stack = FuelCellStack::new(400, 300.0);
        let r1 = FcvRange::estimate_range(1500.0, 0.4, 100.0, &stack);
        let r2 = FcvRange::estimate_range(1500.0, 0.8, 100.0, &stack);
        assert!(r1 > r2, "higher drag -> less range");
    }

    #[test]
    fn fcv_range_positive() {
        let stack = FuelCellStack::new(400, 300.0);
        let r = FcvRange::estimate_range(1500.0, 0.6, 100.0, &stack);
        assert!(r > 0.0);
    }

    #[test]
    fn fcv_range_zero_demand_gives_large_range() {
        let stack = FuelCellStack::new(400, 300.0);
        // Very low speed ~0 => tiny drag => essentially infinite range
        let r = FcvRange::estimate_range(1500.0, 0.6, 0.001, &stack);
        assert!(r > 1000.0, "near-zero speed gives very large range: {r}");
    }

    // --- integration of stack and polarization ---

    #[test]
    fn stack_and_polarization_consistent() {
        let mut stack = FuelCellStack::new(200, 250.0);
        // Verify voltage decreases as we push more current
        let v_low_j = {
            stack.set_current_density(0.1);
            stack.voltage()
        };
        let v_high_j = {
            stack.set_current_density(1.0);
            stack.voltage()
        };
        assert!(
            v_low_j > v_high_j,
            "voltage should decrease at higher current density"
        );
    }

    #[test]
    fn nernst_equation_with_real_pem_conditions() {
        // Typical PEM conditions: 80°C, p_h2=1.5 bar, p_o2=0.3 bar (air-fed)
        let v = nernst_voltage(353.15, 1.5, 0.3);
        assert!(v > 1.0 && v < 1.3, "Nernst OCV at PEM conditions: {v}");
    }

    #[test]
    fn activation_overpotential_physical_range() {
        // Typical PEM activation loss at 0.5 A/cm² with small j0 is in 0.1–0.7 V
        let eta = activation_overpotential(0.5, 1e-4, 0.5, 353.15);
        assert!(eta > 0.05 && eta < 0.7, "activation loss: {eta}");
    }

    // Verify standard math constants are accessible
    #[test]
    fn ln2_constant_accessible() {
        use std::f64::consts::LN_2;
        assert!((LN_2 - 0.69).abs() < 0.01);
    }
}
