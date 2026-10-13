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

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Universal gas constant (J mol⁻¹ K⁻¹).
pub const R_GAS: f64 = 8.314_462;
/// Faraday constant (C mol⁻¹).
pub const FARADAY: f64 = 96_485.0;
/// Number of electrons transferred per mole of H₂ in PEM fuel cell.
pub const N_ELECTRONS: f64 = 2.0;
/// Molar mass of molecular hydrogen (g mol⁻¹).
pub const M_H2_G_PER_MOL: f64 = 2.016;
/// Lower heating value of hydrogen (kWh kg⁻¹).
pub const LHV_H2_KWH_PER_KG: f64 = 33.33;
/// Thermoneutral cell voltage on a higher-heating-value basis (V).
pub const V_THERMONEUTRAL_HHV: f64 = 1.481;
/// Thermoneutral cell voltage on a lower-heating-value basis (V).
pub const V_THERMONEUTRAL_LHV: f64 = 1.254;

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
        self.efficiency_for_thermoneutral_voltage(V_THERMONEUTRAL_HHV)
    }

    /// Thermodynamic efficiency (0–1) relative to the lower heating value of
    /// hydrogen.
    ///
    /// η_LHV = V_stack / (n_cells * 1.254 V).  This is the efficiency that has
    /// to be combined with the LHV of the fuel (33.33 kWh/kg) when converting a
    /// hydrogen mass into usable electrical energy.
    pub fn efficiency_lhv(&self) -> f64 {
        self.efficiency_for_thermoneutral_voltage(V_THERMONEUTRAL_LHV)
    }

    /// Stack voltage divided by the thermoneutral voltage of all cells in series.
    fn efficiency_for_thermoneutral_voltage(&self, v_thermoneutral: f64) -> f64 {
        let v_thermo = v_thermoneutral * self.n_cells as f64;
        if v_thermo < 1e-12 {
            return 0.0;
        }
        (self.voltage() / v_thermo).clamp(0.0, 1.0)
    }

    /// Set operating current density (A cm⁻²).
    pub fn set_current_density(&mut self, j: f64) {
        self.current_density = j.max(0.0);
    }

    /// Hydrogen consumption rate (g s⁻¹) of the whole stack at the current
    /// operating point.
    ///
    /// From Faraday's law: every cell of the series stack carries the same
    /// current `I = j · A_cell` and consumes `I / (2 F)` mol s⁻¹ of H₂, so
    /// the stack consumes
    /// ṁ_H2 = n_cells * I * M_H2 / (2 * F)
    /// where M_H2 = 2.016 g/mol.
    pub fn hydrogen_consumption_rate_gs(&self) -> f64 {
        self.n_cells as f64 * self.current() * M_H2_G_PER_MOL / (N_ELECTRONS * FARADAY)
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
    /// Hydrogen is consumed at the Faraday rate of the whole stack
    /// ([`FuelCellStack::hydrogen_consumption_rate_gs`]).
    ///
    /// A stack cannot generate electricity without fuel: with an empty tank the
    /// stack output is zero and the demand is served from the buffer battery
    /// alone; when the remaining hydrogen lasts for only part of the step the
    /// stack output is scaled by that fraction.  The energy credited to the
    /// drivetrain never exceeds what the stack and battery actually supplied.
    pub fn step(&mut self, dt: f64, power_demand: f64) {
        // Hydrogen the stack would burn during this step and the part of it the
        // tank can still provide.
        let h2_needed_kg = self.stack.hydrogen_consumption_rate_gs() * dt / 1000.0;
        let h2_available_kg = self.tank.current_mass.max(0.0);
        let supply_fraction = if h2_available_kg <= 0.0 {
            0.0
        } else if h2_needed_kg <= h2_available_kg {
            1.0
        } else {
            h2_available_kg / h2_needed_kg
        };

        let fc_power = self.stack.power_output() * supply_fraction;
        let excess = fc_power - power_demand;
        let battery_capacity_wh = self.battery_buffer_kwh * 1000.0;

        // Power that actually reaches the drivetrain bus during this step (W).
        let mut supplied_power = power_demand;
        if excess >= 0.0 {
            // Charge battery with excess
            let charge_wh =
                (excess * dt / 3600.0).min(battery_capacity_wh * (1.0 - self.battery_soc));
            self.battery_soc += charge_wh / battery_capacity_wh.max(1e-12);
        } else {
            // Discharge battery for deficit
            let deficit_wh = (-excess * dt / 3600.0).min(battery_capacity_wh * self.battery_soc);
            self.battery_soc -= deficit_wh / battery_capacity_wh.max(1e-12);
            if dt > 0.0 {
                supplied_power = fc_power + deficit_wh * 3600.0 / dt;
            }
        }

        // Consume hydrogen (the tank clamps to what it still holds)
        self.tank.consume(h2_needed_kg);

        // Track energy delivered
        let motor_power = power_demand
            .min(self.motor_power_kw * 1000.0)
            .min(supplied_power);
        self.energy_delivered_wh += motor_power * self.motor_efficiency * dt / 3600.0;

        self.battery_soc = self.battery_soc.clamp(0.0, 1.0);
    }

    /// Remaining range estimate (km) given current state.
    ///
    /// Uses the hydrogen mass that is actually left in the tank
    /// ([`HydrogenTank::current_mass`]); see [`FcvRange::estimate_range`].
    pub fn remaining_range_km(&self, mass_kg: f64, drag_coeff: f64, speed_kmh: f64) -> f64 {
        FcvRange::estimate_range(
            mass_kg,
            drag_coeff,
            speed_kmh,
            self.tank.current_mass,
            &self.stack,
        )
    }
}

// ---------------------------------------------------------------------------
// FcvRange — range estimation
// ---------------------------------------------------------------------------

/// Fuel cell vehicle range estimator.
pub struct FcvRange;

impl FcvRange {
    /// Estimate the range (km) of a fuel cell vehicle cruising at a constant
    /// speed on a given hydrogen load.
    ///
    /// Steady-state road-load balance:
    /// - Rolling resistance: `F_roll = m * g * Cr`
    /// - Aerodynamic drag:   `F_drag = 0.5 * rho * Cd * A * v^2`
    ///
    /// The usable energy is the lower heating value of the hydrogen in the
    /// tank times the stack efficiency on the same (LHV) basis,
    /// `E = m_H2 * LHV_H2 * eta_LHV` with
    /// `eta_LHV = V_stack / (n_cells * 1.254 V)`
    /// ([`FuelCellStack::efficiency_lhv`]).  Spending `E` against the road
    /// load, `E = (F_roll + F_drag) * range`, gives
    ///
    /// `range = E / (F_roll + F_drag)`
    ///
    /// which is the distance travelled (not the driving time).
    ///
    /// * `mass_kg`    – vehicle mass (kg)
    /// * `drag_coeff` – product `Cd * A` (m²)
    /// * `speed_kmh`  – cruise speed (km/h)
    /// * `tank_kg_h2` – hydrogen mass available for the trip (kg)
    /// * `stack`      – fuel cell stack providing the operating-point efficiency
    ///
    /// Returns `0.0` when no usable energy is available (empty tank or a stack
    /// that produces no voltage) and `f64::INFINITY` when the vehicle meets no
    /// resisting force at all.
    pub fn estimate_range(
        mass_kg: f64,
        drag_coeff: f64,
        speed_kmh: f64,
        tank_kg_h2: f64,
        stack: &FuelCellStack,
    ) -> f64 {
        const GRAVITY: f64 = 9.81; // m/s²
        const ROLLING_RESISTANCE: f64 = 0.01; // dimensionless
        const AIR_DENSITY: f64 = 1.225; // kg/m³
        const JOULES_PER_KWH: f64 = 3.6e6;

        let v = speed_kmh / 3.6; // m/s
        let f_roll = mass_kg.max(0.0) * GRAVITY * ROLLING_RESISTANCE;
        let f_drag = 0.5 * AIR_DENSITY * drag_coeff.max(0.0) * v * v;
        let f_total = f_roll + f_drag; // N

        let energy_j =
            tank_kg_h2.max(0.0) * LHV_H2_KWH_PER_KG * JOULES_PER_KWH * stack.efficiency_lhv();
        if energy_j <= 0.0 {
            return 0.0;
        }
        if f_total <= 0.0 {
            return f64::INFINITY;
        }
        energy_j / f_total / 1000.0 // km
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
        let r1 = FcvRange::estimate_range(1500.0, 0.6, 80.0, 5.0, &stack);
        let r2 = FcvRange::estimate_range(1500.0, 0.6, 130.0, 5.0, &stack);
        assert!(r1 > r2, "lower speed -> longer range due to less drag");
    }

    #[test]
    fn fcv_range_higher_drag_less_range() {
        let stack = FuelCellStack::new(400, 300.0);
        let r1 = FcvRange::estimate_range(1500.0, 0.4, 100.0, 5.0, &stack);
        let r2 = FcvRange::estimate_range(1500.0, 0.8, 100.0, 5.0, &stack);
        assert!(r1 > r2, "higher drag -> less range");
    }

    #[test]
    fn fcv_range_positive() {
        let stack = FuelCellStack::new(400, 300.0);
        let r = FcvRange::estimate_range(1500.0, 0.6, 100.0, 5.0, &stack);
        assert!(r > 0.0);
    }

    #[test]
    fn fcv_range_at_vanishing_speed_is_energy_over_rolling_resistance() {
        let stack = FuelCellStack::new(400, 300.0);
        let tank_kg = 5.0;
        // Without aerodynamic drag only rolling resistance remains, so the
        // range is a finite, speed-independent distance (no sentinel value).
        let r = FcvRange::estimate_range(1500.0, 0.6, 0.001, tank_kg, &stack);
        let energy_j = tank_kg * 33.33 * 3.6e6 * stack.voltage() / (1.254 * 400.0);
        let expected_km = energy_j / (1500.0 * 9.81 * 0.01) / 1000.0;
        assert!(r.is_finite());
        assert!(
            (r - expected_km).abs() <= 1e-6 * expected_km,
            "range {r} km vs {expected_km} km"
        );
        assert!(r > 1000.0, "near-zero speed gives very large range: {r}");
    }

    #[test]
    fn fcv_range_equals_energy_over_road_load() {
        let stack = FuelCellStack::new(400, 300.0);
        let (mass_kg, cd_a, speed_kmh, tank_kg) = (1800.0_f64, 0.65_f64, 100.0_f64, 5.6_f64);
        let v = speed_kmh / 3.6;
        let road_load_n = mass_kg * 9.81 * 0.01 + 0.5 * 1.225 * cd_a * v * v;
        // LHV based stack efficiency: V_stack / (n_cells * 1.254 V)
        let eta_lhv = stack.voltage() / (1.254 * 400.0);
        let energy_j = tank_kg * 33.33 * 3.6e6 * eta_lhv;
        let expected_km = energy_j / road_load_n / 1000.0;
        let got = FcvRange::estimate_range(mass_kg, cd_a, speed_kmh, tank_kg, &stack);
        assert!(
            (got - expected_km).abs() <= 1e-9 * expected_km,
            "range {got} km vs closed form {expected_km} km"
        );
        // A 5.6 kg tank of a passenger FCV is good for a few hundred km.
        assert!(got > 300.0 && got < 700.0, "implausible range {got} km");
    }

    #[test]
    fn fcv_range_scales_linearly_with_tank_mass() {
        let stack = FuelCellStack::new(400, 300.0);
        let r_small = FcvRange::estimate_range(1500.0, 0.6, 100.0, 2.5, &stack);
        let r_large = FcvRange::estimate_range(1500.0, 0.6, 100.0, 5.0, &stack);
        assert!(
            (r_large / r_small - 2.0).abs() < 1e-12,
            "doubling the hydrogen must double the range"
        );
    }

    #[test]
    fn fcv_range_cruise_at_80_beats_130_kmh() {
        let stack = FuelCellStack::new(400, 300.0);
        let r_80 = FcvRange::estimate_range(1500.0, 0.6, 80.0, 5.0, &stack);
        let r_130 = FcvRange::estimate_range(1500.0, 0.6, 130.0, 5.0, &stack);
        assert!(r_80 > r_130);
    }

    #[test]
    fn fcv_range_is_zero_for_empty_tank() {
        let stack = FuelCellStack::new(400, 300.0);
        assert_eq!(
            FcvRange::estimate_range(1500.0, 0.6, 100.0, 0.0, &stack),
            0.0
        );
        assert_eq!(
            FcvRange::estimate_range(1500.0, 0.6, 100.0, -1.0, &stack),
            0.0
        );
    }

    #[test]
    fn fcv_range_is_infinite_without_any_resisting_force() {
        let stack = FuelCellStack::new(400, 300.0);
        let r = FcvRange::estimate_range(0.0, 0.0, 100.0, 5.0, &stack);
        assert!(r.is_infinite() && r > 0.0);
    }

    #[test]
    fn fcv_remaining_range_tracks_hydrogen_left_in_tank() {
        let stack = FuelCellStack::new(400, 300.0);
        let tank = HydrogenTank::new(5.0, 700.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 1.0);
        let full = fcv.remaining_range_km(1500.0, 0.6, 100.0);
        let expected_full = FcvRange::estimate_range(1500.0, 0.6, 100.0, 5.0, &fcv.stack);
        assert!((full - expected_full).abs() <= 1e-9 * expected_full);
        fcv.tank.consume(2.5);
        let half = fcv.remaining_range_km(1500.0, 0.6, 100.0);
        assert!(
            (half / full - 0.5).abs() < 1e-12,
            "half a tank must give half the range"
        );
    }

    // --- Faraday consumption of the whole stack ---

    #[test]
    fn hydrogen_consumption_matches_faraday_for_series_stack() {
        // 300 cells of 300 cm^2 at 0.5 A/cm^2: every cell carries 150 A.
        let stack = FuelCellStack::new(300, 300.0);
        let rate = stack.hydrogen_consumption_rate_gs();
        let expected = 300.0 * 150.0 * 2.016 / (2.0 * 96_485.0);
        assert!((rate - expected).abs() < 1e-12, "rate {rate} vs {expected}");
        assert!((rate - 0.470).abs() < 0.005, "rate {rate} g/s");
    }

    #[test]
    fn hydrogen_consumption_scales_with_cell_count() {
        let single = FuelCellStack::new(150, 300.0).hydrogen_consumption_rate_gs();
        let double = FuelCellStack::new(300, 300.0).hydrogen_consumption_rate_gs();
        assert!((double / single - 2.0).abs() < 1e-12);
    }

    #[test]
    fn hydrogen_consumption_is_consistent_with_stack_efficiency() {
        // The chemical power of the consumed hydrogen (HHV basis) must equal
        // the electrical output divided by the stack efficiency.
        let stack = FuelCellStack::new(250, 280.0);
        let mol_per_s = stack.hydrogen_consumption_rate_gs() / M_H2_G_PER_MOL;
        let chemical_power_w = mol_per_s * N_ELECTRONS * FARADAY * V_THERMONEUTRAL_HHV;
        let ratio = stack.power_output() / chemical_power_w;
        assert!((ratio - stack.efficiency()).abs() < 1e-12);
    }

    #[test]
    fn efficiency_lhv_is_hhv_efficiency_scaled_by_thermoneutral_voltages() {
        let stack = FuelCellStack::new(400, 300.0);
        let ratio = stack.efficiency_lhv() / stack.efficiency();
        assert!((ratio - 1.481 / 1.254).abs() < 1e-12);
    }

    #[test]
    fn fcv_step_drains_tank_at_faraday_rate() {
        let stack = FuelCellStack::new(300, 300.0);
        let rate_gs = stack.hydrogen_consumption_rate_gs();
        let tank = HydrogenTank::new(5.0, 700.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 1.0);
        fcv.step(10.0, 20_000.0);
        let consumed_kg = 5.0 - fcv.tank.current_mass;
        assert!(
            (consumed_kg - rate_gs * 10.0 / 1000.0).abs() < 1e-12,
            "consumed {consumed_kg} kg"
        );
    }

    #[test]
    fn fcv_step_with_empty_tank_produces_no_stack_power() {
        let stack = FuelCellStack::new(300, 300.0);
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.consume(5.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 2.0);
        fcv.battery_soc = 0.5;
        let demand_w = 10_000.0;
        fcv.step(10.0, demand_w);
        // The whole demand is served from the buffer battery.
        let expected_soc = 0.5 - demand_w * 10.0 / 3600.0 / 2000.0;
        assert!(
            (fcv.battery_soc - expected_soc).abs() < 1e-12,
            "soc {} vs {expected_soc}",
            fcv.battery_soc
        );
        assert_eq!(fcv.tank.current_mass, 0.0);
        let expected_wh = demand_w * 0.95 * 10.0 / 3600.0;
        assert!((fcv.energy_delivered_wh - expected_wh).abs() < 1e-9);
    }

    #[test]
    fn fcv_step_with_empty_tank_and_flat_battery_delivers_nothing() {
        let stack = FuelCellStack::new(300, 300.0);
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.consume(5.0);
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 2.0);
        fcv.battery_soc = 0.0;
        fcv.step(10.0, 10_000.0);
        assert_eq!(fcv.battery_soc, 0.0);
        assert_eq!(fcv.energy_delivered_wh, 0.0);
    }

    #[test]
    fn fcv_step_scales_stack_output_when_tank_runs_dry_mid_step() {
        let stack = FuelCellStack::new(300, 300.0);
        let rate_gs = stack.hydrogen_consumption_rate_gs();
        let stack_power_w = stack.power_output();
        let dt = 10.0;
        let needed_kg = rate_gs * dt / 1000.0;
        let mut tank = HydrogenTank::new(5.0, 700.0);
        tank.current_mass = needed_kg / 2.0; // only half of the step's hydrogen is left
        let mut fcv = FuelCellVehicle::new(stack, tank, 100.0, 0.95, 2.0);
        fcv.battery_soc = 0.5;
        // With no demand the stack output goes straight into the battery.
        fcv.step(dt, 0.0);
        let charged_wh = stack_power_w * 0.5 * dt / 3600.0;
        let expected_soc = 0.5 + charged_wh / 2000.0;
        assert!(
            (fcv.battery_soc - expected_soc).abs() < 1e-12,
            "soc {} vs {expected_soc}",
            fcv.battery_soc
        );
        assert!(fcv.tank.current_mass.abs() < 1e-15);
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
