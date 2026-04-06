// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Energy management constraint system for physics simulations.
//!
//! Provides constraints for energy budgets (kinetic/potential/thermal),
//! energy conservation enforcement, power flow, battery models (SOC/SOH),
//! thermal management, energy harvesting, regenerative braking, fuel cells,
//! renewable energy storage, demand response, Hamiltonian constraints,
//! symplectic integrator constraints, energy dissipation bounds, and
//! Lyapunov stability constraints.

// ── Helper functions ──────────────────────────────────────────────────────────

/// Dot product of two 3-vectors.
#[allow(dead_code)]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length squared of a 3-vector.
#[allow(dead_code)]
fn len2_3(v: [f64; 3]) -> f64 {
    v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
}

/// Length of a 3-vector.
#[allow(dead_code)]
fn len3(v: [f64; 3]) -> f64 {
    len2_3(v).sqrt()
}

// ── Energy Budget Constraints ─────────────────────────────────────────────────

/// Energy budget descriptor holding kinetic, potential, and thermal energy budgets.
#[derive(Debug, Clone)]
pub struct EnergyBudget {
    /// Maximum allowed kinetic energy (Joules).
    pub max_kinetic: f64,
    /// Maximum allowed potential energy (Joules).
    pub max_potential: f64,
    /// Maximum allowed thermal energy (Joules).
    pub max_thermal: f64,
    /// Minimum total energy floor (Joules).
    pub min_total: f64,
}

impl EnergyBudget {
    /// Create a new `EnergyBudget`.
    pub fn new(max_kinetic: f64, max_potential: f64, max_thermal: f64, min_total: f64) -> Self {
        Self {
            max_kinetic,
            max_potential,
            max_thermal,
            min_total,
        }
    }

    /// Check whether a given energy state satisfies budget constraints.
    ///
    /// Returns `true` if all budgets are satisfied.
    pub fn is_satisfied(&self, kinetic: f64, potential: f64, thermal: f64) -> bool {
        kinetic <= self.max_kinetic
            && potential <= self.max_potential
            && thermal <= self.max_thermal
            && (kinetic + potential + thermal) >= self.min_total
    }

    /// Compute the total energy violation (sum of excess above each limit, Joules).
    pub fn violation(&self, kinetic: f64, potential: f64, thermal: f64) -> f64 {
        let ke_v = (kinetic - self.max_kinetic).max(0.0);
        let pe_v = (potential - self.max_potential).max(0.0);
        let te_v = (thermal - self.max_thermal).max(0.0);
        ke_v + pe_v + te_v
    }
}

/// Compute kinetic energy of a rigid body.
///
/// KE = 0.5 * m * |v|² + 0.5 * ω^T * I * ω
///
/// * `mass`    – body mass (kg).
/// * `vel`     – linear velocity (m/s).
/// * `omega`   – angular velocity (rad/s).
/// * `inertia` – 3×3 inertia tensor (kg·m²), row-major.
#[allow(dead_code)]
pub fn kinetic_energy(mass: f64, vel: [f64; 3], omega: [f64; 3], inertia: [[f64; 3]; 3]) -> f64 {
    let ke_lin = 0.5 * mass * len2_3(vel);
    // rotational: 0.5 * ω^T I ω
    let iw = [
        inertia[0][0] * omega[0] + inertia[0][1] * omega[1] + inertia[0][2] * omega[2],
        inertia[1][0] * omega[0] + inertia[1][1] * omega[1] + inertia[1][2] * omega[2],
        inertia[2][0] * omega[0] + inertia[2][1] * omega[1] + inertia[2][2] * omega[2],
    ];
    let ke_rot = 0.5 * dot3(omega, iw);
    ke_lin + ke_rot
}

/// Compute gravitational potential energy.
///
/// PE = m * g * h  (where h = position\[2\] in default up-axis convention)
#[allow(dead_code)]
pub fn gravitational_potential_energy(mass: f64, position: [f64; 3], gravity: f64) -> f64 {
    mass * gravity * position[2]
}

/// Compute elastic (spring) potential energy.
///
/// PE_spring = 0.5 * k * x²
#[allow(dead_code)]
pub fn spring_potential_energy(stiffness: f64, displacement: f64) -> f64 {
    0.5 * stiffness * displacement * displacement
}

// ── Energy Conservation Enforcement ──────────────────────────────────────────

/// Energy conservation constraint enforcer.
///
/// Tracks the reference total energy and computes a correction impulse
/// whenever the current energy drifts beyond a tolerance.
#[derive(Debug, Clone)]
pub struct EnergyConservationConstraint {
    /// Reference (initial) total mechanical energy (Joules).
    pub reference_energy: f64,
    /// Allowed fractional energy drift before correction is applied.
    pub tolerance: f64,
    /// Baumgarte-style gain for energy correction.
    pub gain: f64,
}

impl EnergyConservationConstraint {
    /// Create a new `EnergyConservationConstraint`.
    pub fn new(reference_energy: f64, tolerance: f64, gain: f64) -> Self {
        Self {
            reference_energy,
            tolerance,
            gain,
        }
    }

    /// Compute the fractional energy drift: (E_current - E_ref) / E_ref.
    pub fn fractional_drift(&self, current_energy: f64) -> f64 {
        if self.reference_energy.abs() < 1e-15 {
            current_energy
        } else {
            (current_energy - self.reference_energy) / self.reference_energy
        }
    }

    /// Compute the energy correction impulse magnitude.
    ///
    /// Returns a signed scalar that scales correction forces to bring energy
    /// back toward the reference value.
    pub fn correction_impulse(&self, current_energy: f64, dt: f64) -> f64 {
        let drift = self.fractional_drift(current_energy);
        if drift.abs() > self.tolerance {
            -self.gain * drift / dt
        } else {
            0.0
        }
    }
}

// ── Power Flow Constraints ─────────────────────────────────────────────────

/// Power flow constraint describing max/min power transfer between two nodes.
#[derive(Debug, Clone)]
pub struct PowerFlowConstraint {
    /// Maximum power flow (Watts). Positive = forward.
    pub max_power: f64,
    /// Minimum power flow (Watts). Negative = reverse.
    pub min_power: f64,
    /// Line resistance (Ohm or equivalent).
    pub resistance: f64,
}

impl PowerFlowConstraint {
    /// Create a new `PowerFlowConstraint`.
    pub fn new(max_power: f64, min_power: f64, resistance: f64) -> Self {
        Self {
            max_power,
            min_power,
            resistance,
        }
    }

    /// Clamp a requested power value to the constraint bounds.
    pub fn clamp_power(&self, requested: f64) -> f64 {
        requested.max(self.min_power).min(self.max_power)
    }

    /// Compute resistive loss for a given power flow (Watts).
    ///
    /// Loss ≈ R * (P/V)² but simplified here as loss = R * P²  (unit-normalized).
    pub fn resistive_loss(&self, power: f64) -> f64 {
        self.resistance * power * power
    }

    /// Check whether the requested power violates either limit.
    pub fn is_violated(&self, power: f64) -> bool {
        power > self.max_power || power < self.min_power
    }
}

// ── Battery Model Constraints (SOC / SOH) ─────────────────────────────────

/// Battery state-of-charge and state-of-health constraint model.
#[derive(Debug, Clone)]
pub struct BatteryConstraint {
    /// Current state of charge \[0, 1\] (dimensionless).
    pub soc: f64,
    /// State of health \[0, 1\] (dimensionless, 1 = new).
    pub soh: f64,
    /// Nominal capacity (Joules or Wh).
    pub capacity: f64,
    /// Minimum SOC before discharge is blocked.
    pub soc_min: f64,
    /// Maximum SOC before charge is blocked.
    pub soc_max: f64,
    /// Maximum charge rate (W).
    pub max_charge_rate: f64,
    /// Maximum discharge rate (W).
    pub max_discharge_rate: f64,
    /// Internal resistance (Ohm).
    pub internal_resistance: f64,
}

impl BatteryConstraint {
    /// Create a new `BatteryConstraint` with default health.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        capacity: f64,
        soc_min: f64,
        soc_max: f64,
        max_charge_rate: f64,
        max_discharge_rate: f64,
        internal_resistance: f64,
        initial_soc: f64,
    ) -> Self {
        Self {
            soc: initial_soc.max(soc_min).min(soc_max),
            soh: 1.0,
            capacity,
            soc_min,
            soc_max,
            max_charge_rate,
            max_discharge_rate,
            internal_resistance,
        }
    }

    /// Available energy for discharge (Joules).
    pub fn available_energy(&self) -> f64 {
        (self.soc - self.soc_min) * self.capacity * self.soh
    }

    /// Headroom for charging (Joules).
    pub fn charge_headroom(&self) -> f64 {
        (self.soc_max - self.soc) * self.capacity * self.soh
    }

    /// Update SOC given a net power flow (positive = charging) over `dt` seconds.
    ///
    /// Returns the actual power delivered (after clamping).
    pub fn update(&mut self, requested_power: f64, dt: f64) -> f64 {
        let clamped = if requested_power > 0.0 {
            requested_power.min(self.max_charge_rate)
        } else {
            requested_power.max(-self.max_discharge_rate)
        };
        let delta_energy = clamped * dt;
        let new_energy = self.soc * self.capacity * self.soh + delta_energy;
        let new_soc = (new_energy / (self.capacity * self.soh))
            .max(self.soc_min)
            .min(self.soc_max);
        self.soc = new_soc;
        clamped
    }

    /// Terminal voltage approximation: V_oc - I * R_int  (simplified linear model).
    ///
    /// * `open_circuit_voltage` – open-circuit voltage (V).
    /// * `current`              – current draw (A, positive = discharge).
    pub fn terminal_voltage(&self, open_circuit_voltage: f64, current: f64) -> f64 {
        open_circuit_voltage - current * self.internal_resistance
    }

    /// Check whether the battery can supply the requested power.
    pub fn can_discharge(&self, power: f64) -> bool {
        power <= self.max_discharge_rate && self.soc > self.soc_min
    }

    /// Check whether the battery can absorb the requested charge power.
    pub fn can_charge(&self, power: f64) -> bool {
        power <= self.max_charge_rate && self.soc < self.soc_max
    }
}

// ── Thermal Management Constraints ────────────────────────────────────────

/// Thermal constraint for temperature-bounded simulation components.
#[derive(Debug, Clone)]
pub struct ThermalConstraint {
    /// Current temperature (K).
    pub temperature: f64,
    /// Minimum operating temperature (K).
    pub temp_min: f64,
    /// Maximum operating temperature (K).
    pub temp_max: f64,
    /// Thermal capacitance (J/K).
    pub thermal_capacitance: f64,
    /// Thermal conductance to ambient (W/K).
    pub thermal_conductance: f64,
    /// Ambient temperature (K).
    pub ambient_temperature: f64,
}

impl ThermalConstraint {
    /// Create a new `ThermalConstraint`.
    pub fn new(
        initial_temp: f64,
        temp_min: f64,
        temp_max: f64,
        thermal_capacitance: f64,
        thermal_conductance: f64,
        ambient_temperature: f64,
    ) -> Self {
        Self {
            temperature: initial_temp,
            temp_min,
            temp_max,
            thermal_capacitance,
            thermal_conductance,
            ambient_temperature,
        }
    }

    /// Compute natural cooling/heating power (W) from ambient.
    pub fn ambient_exchange_power(&self) -> f64 {
        self.thermal_conductance * (self.ambient_temperature - self.temperature)
    }

    /// Update temperature given heat input power (W) over `dt` seconds.
    pub fn update(&mut self, heat_input: f64, dt: f64) -> f64 {
        let exchange = self.ambient_exchange_power();
        let net_power = heat_input + exchange;
        let delta_t = net_power * dt / self.thermal_capacitance;
        self.temperature = (self.temperature + delta_t)
            .max(self.temp_min)
            .min(self.temp_max);
        self.temperature
    }

    /// Check whether the component is in thermal safe zone.
    pub fn is_in_limits(&self) -> bool {
        self.temperature >= self.temp_min && self.temperature <= self.temp_max
    }

    /// Thermal violation (excess above max or below min, Kelvin).
    pub fn violation(&self) -> f64 {
        let hi = (self.temperature - self.temp_max).max(0.0);
        let lo = (self.temp_min - self.temperature).max(0.0);
        hi + lo
    }
}

// ── Energy Harvesting Constraints ──────────────────────────────────────────

/// Energy harvesting constraint capturing ambient-source power availability.
#[derive(Debug, Clone)]
pub struct EnergyHarvestingConstraint {
    /// Available harvested power (W).
    pub available_power: f64,
    /// Harvesting efficiency \[0, 1\].
    pub efficiency: f64,
    /// Maximum storage buffer (J).
    pub buffer_capacity: f64,
    /// Current buffer fill level (J).
    pub buffer_level: f64,
}

impl EnergyHarvestingConstraint {
    /// Create a new `EnergyHarvestingConstraint`.
    pub fn new(
        available_power: f64,
        efficiency: f64,
        buffer_capacity: f64,
        initial_buffer: f64,
    ) -> Self {
        Self {
            available_power,
            efficiency,
            buffer_capacity,
            buffer_level: initial_buffer.min(buffer_capacity).max(0.0),
        }
    }

    /// Harvested power actually captured (after efficiency).
    pub fn harvested_power(&self) -> f64 {
        self.available_power * self.efficiency
    }

    /// Update buffer over `dt` seconds, given consumed power.
    ///
    /// Returns net energy added to buffer (J).
    pub fn update(&mut self, consumed_power: f64, dt: f64) -> f64 {
        let harvested = self.harvested_power() * dt;
        let consumed = consumed_power * dt;
        let delta = harvested - consumed;
        let new_level = (self.buffer_level + delta)
            .max(0.0)
            .min(self.buffer_capacity);
        let actual_delta = new_level - self.buffer_level;
        self.buffer_level = new_level;
        actual_delta
    }

    /// Fraction of buffer currently filled \[0, 1\].
    pub fn fill_fraction(&self) -> f64 {
        if self.buffer_capacity < 1e-15 {
            0.0
        } else {
            self.buffer_level / self.buffer_capacity
        }
    }
}

// ── Regenerative Braking Constraint ───────────────────────────────────────

/// Regenerative braking constraint model.
///
/// Converts kinetic energy during deceleration into recoverable electrical energy.
#[derive(Debug, Clone)]
pub struct RegenerativeBrakingConstraint {
    /// Regeneration efficiency \[0, 1\].
    pub regen_efficiency: f64,
    /// Maximum regenerative power (W).
    pub max_regen_power: f64,
    /// Minimum vehicle speed for regen to be active (m/s).
    pub min_speed: f64,
}

impl RegenerativeBrakingConstraint {
    /// Create a new `RegenerativeBrakingConstraint`.
    pub fn new(regen_efficiency: f64, max_regen_power: f64, min_speed: f64) -> Self {
        Self {
            regen_efficiency,
            max_regen_power,
            min_speed,
        }
    }

    /// Compute recoverable braking power for a given vehicle speed and
    /// deceleration force (N).
    ///
    /// Returns recovered power (W), always non-negative.
    pub fn recovered_power(&self, speed: f64, braking_force: f64) -> f64 {
        if speed < self.min_speed {
            return 0.0;
        }
        let raw_power = (speed * braking_force).abs();
        let regen = raw_power * self.regen_efficiency;
        regen.min(self.max_regen_power)
    }

    /// Compute the effective braking force after subtracting regen from total braking.
    ///
    /// Returns the mechanical (friction) braking force required.
    pub fn friction_braking_force(&self, total_braking_force: f64, speed: f64) -> f64 {
        if speed < self.min_speed {
            return total_braking_force;
        }
        let regen_force =
            (self.max_regen_power / speed.max(1e-6) / self.regen_efficiency.max(1e-6))
                .min(total_braking_force.abs());
        (total_braking_force.abs() - regen_force).max(0.0)
            * if total_braking_force < 0.0 { -1.0 } else { 1.0 }
    }
}

// ── Fuel Cell Constraints ─────────────────────────────────────────────────

/// Fuel cell constraint model including polarization losses.
#[derive(Debug, Clone)]
pub struct FuelCellConstraint {
    /// Nominal open-circuit voltage (V).
    pub open_circuit_voltage: f64,
    /// Maximum power output (W).
    pub max_power: f64,
    /// Activation loss coefficient (V).
    pub activation_loss: f64,
    /// Ohmic resistance (Ohm).
    pub ohmic_resistance: f64,
    /// Mass transport loss coefficient (V).
    pub concentration_loss: f64,
    /// Current hydrogen supply rate (mol/s).
    pub hydrogen_supply_rate: f64,
    /// Faraday constant (C/mol).
    pub faraday_constant: f64,
}

impl FuelCellConstraint {
    /// Create a new `FuelCellConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        open_circuit_voltage: f64,
        max_power: f64,
        activation_loss: f64,
        ohmic_resistance: f64,
        concentration_loss: f64,
        hydrogen_supply_rate: f64,
        faraday_constant: f64,
    ) -> Self {
        Self {
            open_circuit_voltage,
            max_power,
            activation_loss,
            ohmic_resistance,
            concentration_loss,
            hydrogen_supply_rate,
            faraday_constant,
        }
    }

    /// Compute cell voltage at a given current density (A/m²).
    ///
    /// V = V_oc - V_act - V_ohm - V_conc
    pub fn cell_voltage(&self, current_density: f64) -> f64 {
        let v_act = self.activation_loss * (current_density / 0.01 + 1.0).ln();
        let v_ohm = self.ohmic_resistance * current_density;
        let v_conc = self.concentration_loss * (current_density / 1.0).max(0.0);
        (self.open_circuit_voltage - v_act - v_ohm - v_conc).max(0.0)
    }

    /// Maximum hydrogen consumption rate at full power (mol/s).
    pub fn max_hydrogen_consumption(&self) -> f64 {
        // 2F per mole H2 in PEM fuel cell
        self.max_power / (2.0 * self.faraday_constant * self.open_circuit_voltage)
    }

    /// Check whether hydrogen supply is sufficient for the requested power.
    pub fn has_sufficient_fuel(&self, requested_power: f64) -> bool {
        let required = requested_power / (2.0 * self.faraday_constant * self.open_circuit_voltage);
        required <= self.hydrogen_supply_rate
    }
}

// ── Renewable Energy Storage Constraint ───────────────────────────────────

/// Renewable energy storage constraint (e.g., pumped hydro, compressed air).
#[derive(Debug, Clone)]
pub struct RenewableStorageConstraint {
    /// Total storage capacity (J).
    pub capacity: f64,
    /// Current stored energy (J).
    pub stored_energy: f64,
    /// Charge efficiency \[0, 1\].
    pub charge_efficiency: f64,
    /// Discharge efficiency \[0, 1\].
    pub discharge_efficiency: f64,
    /// Self-discharge rate per second \[0, 1\].
    pub self_discharge_rate: f64,
    /// Maximum charge power (W).
    pub max_charge_power: f64,
    /// Maximum discharge power (W).
    pub max_discharge_power: f64,
}

impl RenewableStorageConstraint {
    /// Create a new `RenewableStorageConstraint`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        capacity: f64,
        initial_energy: f64,
        charge_efficiency: f64,
        discharge_efficiency: f64,
        self_discharge_rate: f64,
        max_charge_power: f64,
        max_discharge_power: f64,
    ) -> Self {
        Self {
            capacity,
            stored_energy: initial_energy.min(capacity).max(0.0),
            charge_efficiency,
            discharge_efficiency,
            self_discharge_rate,
            max_charge_power,
            max_discharge_power,
        }
    }

    /// State of charge \[0, 1\].
    pub fn soc(&self) -> f64 {
        if self.capacity < 1e-15 {
            0.0
        } else {
            self.stored_energy / self.capacity
        }
    }

    /// Update stored energy given net power (positive = charge, negative = discharge)
    /// and time step `dt`.
    ///
    /// Returns actual power flow (W) after applying limits and efficiency.
    pub fn update(&mut self, requested_power: f64, dt: f64) -> f64 {
        // Apply self-discharge
        self.stored_energy *= 1.0 - self.self_discharge_rate * dt;
        self.stored_energy = self.stored_energy.max(0.0);

        let clamped = if requested_power >= 0.0 {
            requested_power.min(self.max_charge_power)
        } else {
            requested_power.max(-self.max_discharge_power)
        };

        let delta = if clamped >= 0.0 {
            clamped * self.charge_efficiency * dt
        } else {
            clamped / self.discharge_efficiency.max(1e-6) * dt
        };

        let new_stored = (self.stored_energy + delta).max(0.0).min(self.capacity);
        self.stored_energy = new_stored;
        clamped
    }
}

// ── Demand Response Constraint ─────────────────────────────────────────────

/// Demand response constraint for grid-interactive energy management.
#[derive(Debug, Clone)]
pub struct DemandResponseConstraint {
    /// Baseline power demand (W).
    pub baseline_demand: f64,
    /// Current curtailment fraction \[0, 1\].
    pub curtailment_fraction: f64,
    /// Maximum allowed curtailment fraction.
    pub max_curtailment: f64,
    /// Minimum response time (s).
    pub min_response_time: f64,
    /// Grid signal: 1.0 = normal, < 1.0 = curtail, > 1.0 = increase.
    pub grid_signal: f64,
}

impl DemandResponseConstraint {
    /// Create a new `DemandResponseConstraint`.
    pub fn new(baseline_demand: f64, max_curtailment: f64, min_response_time: f64) -> Self {
        Self {
            baseline_demand,
            curtailment_fraction: 0.0,
            max_curtailment,
            min_response_time,
            grid_signal: 1.0,
        }
    }

    /// Compute the allowed power demand after curtailment.
    pub fn allowed_demand(&self) -> f64 {
        self.baseline_demand * self.grid_signal * (1.0 - self.curtailment_fraction)
    }

    /// Apply a demand response event, updating curtailment fraction.
    ///
    /// `target_curtailment` is the desired fraction \[0, max_curtailment\].
    pub fn apply_event(&mut self, target_curtailment: f64) {
        self.curtailment_fraction = target_curtailment.max(0.0).min(self.max_curtailment);
    }

    /// Compute load shedding (W) relative to baseline.
    pub fn load_shedding(&self) -> f64 {
        self.baseline_demand - self.allowed_demand()
    }
}

// ── Hamiltonian Constraints ────────────────────────────────────────────────

/// Hamiltonian constraint ensuring the total mechanical energy equals a target value.
///
/// H = KE + PE = H_target
#[derive(Debug, Clone)]
pub struct HamiltonianConstraint {
    /// Target Hamiltonian value (J).
    pub target_hamiltonian: f64,
    /// Allowed deviation tolerance (J).
    pub tolerance: f64,
    /// Accumulated Hamiltonian violation (J·s).
    pub accumulated_violation: f64,
}

impl HamiltonianConstraint {
    /// Create a new `HamiltonianConstraint`.
    pub fn new(target_hamiltonian: f64, tolerance: f64) -> Self {
        Self {
            target_hamiltonian,
            tolerance,
            accumulated_violation: 0.0,
        }
    }

    /// Evaluate Hamiltonian error: H_current - H_target.
    pub fn hamiltonian_error(&self, ke: f64, pe: f64) -> f64 {
        (ke + pe) - self.target_hamiltonian
    }

    /// Check whether the Hamiltonian constraint is satisfied.
    pub fn is_satisfied(&self, ke: f64, pe: f64) -> bool {
        self.hamiltonian_error(ke, pe).abs() <= self.tolerance
    }

    /// Accumulate the constraint violation over a timestep `dt`.
    pub fn accumulate(&mut self, ke: f64, pe: f64, dt: f64) {
        self.accumulated_violation += self.hamiltonian_error(ke, pe).abs() * dt;
    }

    /// Compute correction velocity scale to restore the Hamiltonian.
    ///
    /// Returns a multiplicative scale applied to the velocity field.
    pub fn velocity_correction_scale(&self, ke: f64, pe: f64) -> f64 {
        let error = self.hamiltonian_error(ke, pe);
        if ke < 1e-15 || error.abs() < self.tolerance {
            return 1.0;
        }
        // Scale KE to match target: KE_corrected = KE + delta_H => v_scale = sqrt(KE_corrected / KE)
        let ke_target = (self.target_hamiltonian - pe).max(0.0);
        (ke_target / ke).sqrt()
    }
}

// ── Symplectic Integrator Constraints ─────────────────────────────────────

/// Symplectic integration constraint for Hamiltonian systems.
///
/// Enforces symplecticity preservation in numerical integration.
#[derive(Debug, Clone)]
pub struct SymplecticConstraint {
    /// Symplectic tolerance for volume preservation error.
    pub symplectic_tolerance: f64,
    /// Integration order (1 = Euler, 2 = leapfrog, 4 = Forest-Ruth).
    pub integration_order: usize,
    /// Accumulated symplectic error metric.
    pub accumulated_error: f64,
}

impl SymplecticConstraint {
    /// Create a new `SymplecticConstraint`.
    pub fn new(symplectic_tolerance: f64, integration_order: usize) -> Self {
        Self {
            symplectic_tolerance,
            integration_order,
            accumulated_error: 0.0,
        }
    }

    /// Compute the leapfrog (Störmer-Verlet) position update.
    ///
    /// q(t + dt) = q(t) + dt * p(t+dt/2) / m
    pub fn leapfrog_position_update(
        &self,
        position: [f64; 3],
        half_step_momentum: [f64; 3],
        mass: f64,
        dt: f64,
    ) -> [f64; 3] {
        let inv_m = if mass.abs() > 1e-15 { 1.0 / mass } else { 0.0 };
        [
            position[0] + dt * half_step_momentum[0] * inv_m,
            position[1] + dt * half_step_momentum[1] * inv_m,
            position[2] + dt * half_step_momentum[2] * inv_m,
        ]
    }

    /// Compute the leapfrog half-step momentum update.
    ///
    /// p(t + dt/2) = p(t - dt/2) + dt * F(q(t))
    pub fn leapfrog_momentum_update(
        &self,
        momentum: [f64; 3],
        force: [f64; 3],
        dt: f64,
    ) -> [f64; 3] {
        [
            momentum[0] + dt * force[0],
            momentum[1] + dt * force[1],
            momentum[2] + dt * force[2],
        ]
    }

    /// Estimate the symplectic area preservation error for a 2D phase space.
    ///
    /// For an ideal symplectic integrator, the phase space area is conserved.
    /// This computes | |J| - 1 | where J is the Jacobian determinant approximated
    /// from incremental position/momentum changes.
    pub fn phase_area_error(&mut self, dq: f64, dp: f64, dq_old: f64, dp_old: f64) -> f64 {
        // Approximate Jacobian as 2x2: J = [[dq/dq0, dq/dp0], [dp/dq0, dp/dp0]]
        // Simplified: use ratio of areas
        let area_new = dq * dp;
        let area_old = dq_old * dp_old;
        let error = if area_old.abs() > 1e-15 {
            ((area_new / area_old).abs() - 1.0).abs()
        } else {
            0.0
        };
        self.accumulated_error += error;
        error
    }

    /// Forest-Ruth coefficients for 4th-order symplectic integration.
    ///
    /// Returns (theta, xi) coefficients.
    pub fn forest_ruth_coefficients() -> (f64, f64) {
        let theta = 1.0 / (2.0 - 2.0_f64.powf(1.0 / 3.0));
        let xi = 1.0 - 2.0 * theta;
        (theta, xi)
    }
}

// ── Energy Dissipation Bounds ─────────────────────────────────────────────

/// Energy dissipation constraint bounding allowable dissipation rates.
#[derive(Debug, Clone)]
pub struct EnergyDissipationConstraint {
    /// Maximum allowed dissipation rate (W).
    pub max_dissipation_rate: f64,
    /// Minimum required dissipation rate (W) — for active cooling requirements.
    pub min_dissipation_rate: f64,
    /// Accumulated total dissipation (J).
    pub total_dissipation: f64,
    /// Dissipation budget limit (J).
    pub dissipation_budget: f64,
}

impl EnergyDissipationConstraint {
    /// Create a new `EnergyDissipationConstraint`.
    pub fn new(
        max_dissipation_rate: f64,
        min_dissipation_rate: f64,
        dissipation_budget: f64,
    ) -> Self {
        Self {
            max_dissipation_rate,
            min_dissipation_rate,
            total_dissipation: 0.0,
            dissipation_budget,
        }
    }

    /// Compute viscous dissipation: D = c * |v|²
    ///
    /// * `damping_coeff` – viscous damping coefficient (N·s/m).
    /// * `velocity`      – body velocity (m/s).
    pub fn viscous_dissipation(damping_coeff: f64, velocity: [f64; 3]) -> f64 {
        damping_coeff * len2_3(velocity)
    }

    /// Clamp dissipation rate to allowed bounds.
    pub fn clamp_dissipation(&self, rate: f64) -> f64 {
        rate.max(self.min_dissipation_rate)
            .min(self.max_dissipation_rate)
    }

    /// Accumulate dissipation over time step `dt`.
    ///
    /// Returns `true` if the budget has not been exceeded.
    pub fn accumulate(&mut self, rate: f64, dt: f64) -> bool {
        self.total_dissipation += rate.abs() * dt;
        self.total_dissipation <= self.dissipation_budget
    }

    /// Remaining dissipation budget (J).
    pub fn remaining_budget(&self) -> f64 {
        (self.dissipation_budget - self.total_dissipation).max(0.0)
    }
}

// ── Lyapunov Stability Constraints ────────────────────────────────────────

/// Lyapunov stability constraint based on an energy-like Lyapunov function.
///
/// Enforces V̇ ≤ -α * V to guarantee asymptotic stability.
#[derive(Debug, Clone)]
pub struct LyapunovConstraint {
    /// Decay rate α (s⁻¹), must be positive.
    pub decay_rate: f64,
    /// Current Lyapunov function value V (must be non-negative).
    pub lyapunov_value: f64,
    /// Maximum allowed Lyapunov value before control action is required.
    pub lyapunov_threshold: f64,
    /// Quadratic weight matrix (3×3) for V = x^T P x.
    pub weight_matrix: [[f64; 3]; 3],
}

impl LyapunovConstraint {
    /// Create a new `LyapunovConstraint` with an identity weight matrix.
    pub fn new(decay_rate: f64, lyapunov_threshold: f64) -> Self {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        Self {
            decay_rate,
            lyapunov_value: 0.0,
            lyapunov_threshold,
            weight_matrix: identity,
        }
    }

    /// Evaluate the quadratic Lyapunov function V = x^T P x.
    pub fn evaluate(&self, state: [f64; 3]) -> f64 {
        let px = [
            self.weight_matrix[0][0] * state[0]
                + self.weight_matrix[0][1] * state[1]
                + self.weight_matrix[0][2] * state[2],
            self.weight_matrix[1][0] * state[0]
                + self.weight_matrix[1][1] * state[1]
                + self.weight_matrix[1][2] * state[2],
            self.weight_matrix[2][0] * state[0]
                + self.weight_matrix[2][1] * state[1]
                + self.weight_matrix[2][2] * state[2],
        ];
        dot3(state, px)
    }

    /// Compute the required V̇ upper bound: V̇ ≤ -α * V.
    pub fn lyapunov_dot_bound(&self) -> f64 {
        -self.decay_rate * self.lyapunov_value
    }

    /// Check whether the current Lyapunov derivative satisfies the stability condition.
    pub fn is_stable(&self, v_dot: f64) -> bool {
        v_dot <= self.lyapunov_dot_bound()
    }

    /// Compute a Lyapunov-based control correction.
    ///
    /// Returns a correction vector to push V̇ below the bound.
    pub fn correction(&self, state: [f64; 3], v_dot_current: f64) -> [f64; 3] {
        let excess = (v_dot_current - self.lyapunov_dot_bound()).max(0.0);
        let norm = len3(state).max(1e-15);
        let scale = -excess / (norm * norm);
        [state[0] * scale, state[1] * scale, state[2] * scale]
    }
}

// ── Composite Energy Manager ───────────────────────────────────────────────

/// A composite energy manager that coordinates multiple energy constraints.
#[derive(Debug, Clone)]
pub struct EnergyManager {
    /// Battery subsystem.
    pub battery: BatteryConstraint,
    /// Thermal subsystem.
    pub thermal: ThermalConstraint,
    /// Energy dissipation tracker.
    pub dissipation: EnergyDissipationConstraint,
    /// Energy conservation enforcer.
    pub conservation: EnergyConservationConstraint,
    /// Simulation time (s).
    pub time: f64,
}

impl EnergyManager {
    /// Create a new `EnergyManager` with default subsystems.
    pub fn new(
        battery: BatteryConstraint,
        thermal: ThermalConstraint,
        dissipation: EnergyDissipationConstraint,
        reference_energy: f64,
    ) -> Self {
        Self {
            battery,
            thermal,
            dissipation,
            conservation: EnergyConservationConstraint::new(reference_energy, 0.01, 1.0),
            time: 0.0,
        }
    }

    /// Advance the energy manager by one timestep `dt`.
    ///
    /// * `electrical_power` – net electrical demand (W, positive = consuming).
    /// * `heat_input`       – net heat input to thermal system (W).
    /// * `dissipation_rate` – current energy dissipation rate (W).
    /// * `ke`               – current kinetic energy (J).
    /// * `pe`               – current potential energy (J).
    ///
    /// Returns `true` if all constraints are satisfied.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        electrical_power: f64,
        heat_input: f64,
        dissipation_rate: f64,
        ke: f64,
        pe: f64,
        dt: f64,
    ) -> bool {
        self.time += dt;
        self.battery.update(-electrical_power, dt);
        self.thermal.update(heat_input, dt);
        let _budget_ok = self.dissipation.accumulate(dissipation_rate, dt);
        let _correction = self.conservation.correction_impulse(ke + pe, dt);

        self.battery.soc > self.battery.soc_min
            && self.thermal.is_in_limits()
            && self.dissipation.total_dissipation <= self.dissipation.dissipation_budget
    }
}

// ── Solar Irradiance Model ─────────────────────────────────────────────────

/// Simple solar irradiance model for energy harvesting applications.
#[derive(Debug, Clone)]
pub struct SolarIrradianceModel {
    /// Peak irradiance (W/m²).
    pub peak_irradiance: f64,
    /// Panel area (m²).
    pub panel_area: f64,
    /// Panel efficiency \[0, 1\].
    pub panel_efficiency: f64,
    /// Atmospheric attenuation factor \[0, 1\].
    pub atmospheric_attenuation: f64,
}

impl SolarIrradianceModel {
    /// Create a new `SolarIrradianceModel`.
    pub fn new(
        peak_irradiance: f64,
        panel_area: f64,
        panel_efficiency: f64,
        atmospheric_attenuation: f64,
    ) -> Self {
        Self {
            peak_irradiance,
            panel_area,
            panel_efficiency,
            atmospheric_attenuation,
        }
    }

    /// Compute available solar power at a given solar elevation angle (radians).
    pub fn available_power(&self, elevation_angle: f64) -> f64 {
        if elevation_angle <= 0.0 {
            return 0.0;
        }
        let irradiance =
            self.peak_irradiance * elevation_angle.sin() * self.atmospheric_attenuation;
        irradiance * self.panel_area * self.panel_efficiency
    }

    /// Compute daily energy yield (J) assuming sinusoidal solar profile.
    ///
    /// * `day_length_hours` – number of daylight hours.
    pub fn daily_energy_yield(&self, day_length_hours: f64) -> f64 {
        // Average factor of (2/π) for sinusoidal elevation over the day
        let avg_factor = 2.0 / std::f64::consts::PI;
        let day_length_s = day_length_hours * 3600.0;
        self.peak_irradiance
            * avg_factor
            * self.atmospheric_attenuation
            * self.panel_area
            * self.panel_efficiency
            * day_length_s
    }
}

// ── Wind Energy Model ──────────────────────────────────────────────────────

/// Wind energy harvesting model using Betz law.
#[derive(Debug, Clone)]
pub struct WindEnergyModel {
    /// Air density (kg/m³).
    pub air_density: f64,
    /// Rotor swept area (m²).
    pub rotor_area: f64,
    /// Power coefficient Cp \[0, 0.593\] (Betz limit).
    pub power_coefficient: f64,
    /// Cut-in wind speed (m/s).
    pub cut_in_speed: f64,
    /// Cut-out wind speed (m/s).
    pub cut_out_speed: f64,
    /// Rated wind speed (m/s).
    pub rated_speed: f64,
    /// Rated power output (W).
    pub rated_power: f64,
}

impl WindEnergyModel {
    /// Create a new `WindEnergyModel`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        air_density: f64,
        rotor_area: f64,
        power_coefficient: f64,
        cut_in_speed: f64,
        cut_out_speed: f64,
        rated_speed: f64,
        rated_power: f64,
    ) -> Self {
        Self {
            air_density,
            rotor_area,
            power_coefficient: power_coefficient.min(0.593),
            cut_in_speed,
            cut_out_speed,
            rated_speed,
            rated_power,
        }
    }

    /// Compute wind power output (W) for a given wind speed (m/s).
    pub fn power_output(&self, wind_speed: f64) -> f64 {
        if wind_speed < self.cut_in_speed || wind_speed > self.cut_out_speed {
            return 0.0;
        }
        if wind_speed >= self.rated_speed {
            return self.rated_power;
        }
        let raw =
            0.5 * self.air_density * self.rotor_area * self.power_coefficient * wind_speed.powi(3);
        raw.min(self.rated_power)
    }

    /// Betz limit theoretical maximum power (W) for given wind speed.
    pub fn betz_limit_power(&self, wind_speed: f64) -> f64 {
        0.5 * self.air_density * self.rotor_area * 0.593 * wind_speed.powi(3)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── EnergyBudget ─────────────────────────────────────────────────────────

    #[test]
    fn test_energy_budget_satisfied() {
        let budget = EnergyBudget::new(100.0, 200.0, 50.0, 0.0);
        assert!(budget.is_satisfied(50.0, 100.0, 25.0));
    }

    #[test]
    fn test_energy_budget_violated() {
        let budget = EnergyBudget::new(100.0, 200.0, 50.0, 0.0);
        assert!(!budget.is_satisfied(150.0, 100.0, 25.0));
    }

    #[test]
    fn test_energy_budget_violation_magnitude() {
        let budget = EnergyBudget::new(100.0, 200.0, 50.0, 0.0);
        let v = budget.violation(120.0, 100.0, 60.0);
        // KE excess = 20, thermal excess = 10
        assert!((v - 30.0).abs() < 1e-10, "v={v}");
    }

    // ── Kinetic / Potential Energy ────────────────────────────────────────────

    #[test]
    fn test_kinetic_energy_linear_only() {
        let inertia = [[0.0; 3]; 3];
        let ke = kinetic_energy(2.0, [3.0, 0.0, 0.0], [0.0; 3], inertia);
        // 0.5 * 2 * 9 = 9
        assert!((ke - 9.0).abs() < 1e-10, "ke={ke}");
    }

    #[test]
    fn test_gravitational_pe() {
        let pe = gravitational_potential_energy(10.0, [0.0, 0.0, 5.0], 9.81);
        assert!((pe - 490.5).abs() < 1e-6, "pe={pe}");
    }

    #[test]
    fn test_spring_pe() {
        let pe = spring_potential_energy(200.0, 0.1);
        assert!((pe - 1.0).abs() < 1e-10, "pe={pe}");
    }

    // ── EnergyConservationConstraint ─────────────────────────────────────────

    #[test]
    fn test_conservation_fractional_drift_zero() {
        let c = EnergyConservationConstraint::new(1000.0, 0.01, 1.0);
        let drift = c.fractional_drift(1000.0);
        assert!(drift.abs() < 1e-10);
    }

    #[test]
    fn test_conservation_correction_within_tolerance() {
        let c = EnergyConservationConstraint::new(1000.0, 0.01, 1.0);
        let corr = c.correction_impulse(1005.0, 0.01);
        // drift = 0.005 < 0.01 → no correction
        assert!(corr.abs() < 1e-10, "corr={corr}");
    }

    #[test]
    fn test_conservation_correction_outside_tolerance() {
        let c = EnergyConservationConstraint::new(1000.0, 0.001, 1.0);
        let corr = c.correction_impulse(1050.0, 0.01);
        assert!(corr.abs() > 0.0, "corr={corr}");
    }

    // ── PowerFlowConstraint ───────────────────────────────────────────────────

    #[test]
    fn test_power_flow_clamp() {
        let pf = PowerFlowConstraint::new(1000.0, -500.0, 0.01);
        assert!((pf.clamp_power(2000.0) - 1000.0).abs() < 1e-10);
        assert!((pf.clamp_power(-1000.0) - (-500.0)).abs() < 1e-10);
    }

    #[test]
    fn test_power_flow_resistive_loss() {
        let pf = PowerFlowConstraint::new(1000.0, -500.0, 0.1);
        let loss = pf.resistive_loss(100.0);
        assert!((loss - 1000.0).abs() < 1e-10, "loss={loss}");
    }

    // ── BatteryConstraint ─────────────────────────────────────────────────────

    #[test]
    fn test_battery_available_energy() {
        let bat = BatteryConstraint::new(10000.0, 0.1, 0.9, 5000.0, 5000.0, 0.05, 0.5);
        // (0.5 - 0.1) * 10000 * 1.0 = 4000
        let avail = bat.available_energy();
        assert!((avail - 4000.0).abs() < 1e-6, "avail={avail}");
    }

    #[test]
    fn test_battery_update_charge() {
        let mut bat = BatteryConstraint::new(10000.0, 0.1, 0.9, 5000.0, 5000.0, 0.05, 0.5);
        bat.update(1000.0, 1.0); // 1000 W for 1 s = 1000 J
        // new_energy = 5000 + 1000 = 6000; soc = 6000/10000 = 0.6
        assert!((bat.soc - 0.6).abs() < 1e-6, "soc={}", bat.soc);
    }

    #[test]
    fn test_battery_cannot_discharge_below_min() {
        let bat = BatteryConstraint::new(10000.0, 0.5, 0.9, 5000.0, 5000.0, 0.05, 0.5);
        assert!(!bat.can_discharge(100.0));
    }

    // ── ThermalConstraint ─────────────────────────────────────────────────────

    #[test]
    fn test_thermal_in_limits() {
        let tc = ThermalConstraint::new(300.0, 250.0, 400.0, 1000.0, 10.0, 293.0);
        assert!(tc.is_in_limits());
    }

    #[test]
    fn test_thermal_violation_over() {
        let tc = ThermalConstraint::new(420.0, 250.0, 400.0, 1000.0, 10.0, 293.0);
        assert!((tc.violation() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_update_cooling() {
        let mut tc = ThermalConstraint::new(400.0, 250.0, 500.0, 1000.0, 100.0, 300.0);
        // Exchange power = 100 * (300 - 400) = -10000 W; delta_T = -10000 * 1 / 1000 = -10
        let new_t = tc.update(0.0, 1.0);
        assert!((new_t - 390.0).abs() < 1e-6, "new_t={new_t}");
    }

    // ── EnergyHarvestingConstraint ────────────────────────────────────────────

    #[test]
    fn test_harvesting_harvested_power() {
        let eh = EnergyHarvestingConstraint::new(1000.0, 0.2, 5000.0, 1000.0);
        assert!((eh.harvested_power() - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_harvesting_fill_fraction() {
        let eh = EnergyHarvestingConstraint::new(1000.0, 0.2, 5000.0, 2500.0);
        assert!((eh.fill_fraction() - 0.5).abs() < 1e-10);
    }

    // ── RegenerativeBrakingConstraint ─────────────────────────────────────────

    #[test]
    fn test_regen_braking_recovered_power() {
        let rb = RegenerativeBrakingConstraint::new(0.8, 50000.0, 1.0);
        let p = rb.recovered_power(20.0, 1000.0);
        // raw = 20000, regen = 16000
        assert!((p - 16000.0).abs() < 1e-6, "p={p}");
    }

    #[test]
    fn test_regen_braking_below_min_speed() {
        let rb = RegenerativeBrakingConstraint::new(0.8, 50000.0, 5.0);
        let p = rb.recovered_power(2.0, 1000.0);
        assert!(p.abs() < 1e-10);
    }

    // ── FuelCellConstraint ────────────────────────────────────────────────────

    #[test]
    fn test_fuel_cell_voltage_zero_current() {
        let fc = FuelCellConstraint::new(1.2, 10000.0, 0.05, 0.01, 0.02, 0.01, 96485.0);
        let v = fc.cell_voltage(0.0);
        // v_act = 0.05 * ln(1) = 0, v_ohm = 0, v_conc = 0 => V_oc
        assert!((v - 1.2).abs() < 1e-6, "v={v}");
    }

    #[test]
    fn test_fuel_cell_has_fuel() {
        let fc = FuelCellConstraint::new(1.2, 10000.0, 0.05, 0.01, 0.02, 1.0, 96485.0);
        assert!(fc.has_sufficient_fuel(1.0)); // trivially small power
    }

    // ── RenewableStorageConstraint ────────────────────────────────────────────

    #[test]
    fn test_renewable_storage_soc() {
        let rs = RenewableStorageConstraint::new(10000.0, 5000.0, 0.9, 0.9, 0.0, 2000.0, 2000.0);
        assert!((rs.soc() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_renewable_storage_update_charge() {
        let mut rs =
            RenewableStorageConstraint::new(10000.0, 5000.0, 0.9, 0.9, 0.0, 2000.0, 2000.0);
        rs.update(1000.0, 1.0);
        // 5000 + 1000 * 0.9 = 5900 => soc = 0.59
        assert!((rs.soc() - 0.59).abs() < 1e-6, "soc={}", rs.soc());
    }

    // ── DemandResponseConstraint ──────────────────────────────────────────────

    #[test]
    fn test_demand_response_allowed_demand() {
        let mut dr = DemandResponseConstraint::new(10000.0, 0.5, 1.0);
        dr.grid_signal = 1.0;
        dr.apply_event(0.3);
        let allowed = dr.allowed_demand();
        assert!((allowed - 7000.0).abs() < 1e-6, "allowed={allowed}");
    }

    #[test]
    fn test_demand_response_load_shedding() {
        let mut dr = DemandResponseConstraint::new(10000.0, 0.5, 1.0);
        dr.apply_event(0.2);
        let shed = dr.load_shedding();
        assert!((shed - 2000.0).abs() < 1e-6, "shed={shed}");
    }

    // ── HamiltonianConstraint ─────────────────────────────────────────────────

    #[test]
    fn test_hamiltonian_error_zero() {
        let hc = HamiltonianConstraint::new(100.0, 1.0);
        let err = hc.hamiltonian_error(60.0, 40.0);
        assert!(err.abs() < 1e-10);
    }

    #[test]
    fn test_hamiltonian_satisfied() {
        let hc = HamiltonianConstraint::new(100.0, 1.0);
        assert!(hc.is_satisfied(60.0, 40.5));
    }

    #[test]
    fn test_hamiltonian_velocity_scale_correct() {
        let hc = HamiltonianConstraint::new(100.0, 0.1);
        // KE=50, PE=40 => H=90, target=100 => scale = sqrt(60/50) = sqrt(1.2)
        let scale = hc.velocity_correction_scale(50.0, 40.0);
        let expected = (60.0_f64 / 50.0).sqrt();
        assert!((scale - expected).abs() < 1e-8, "scale={scale}");
    }

    // ── SymplecticConstraint ──────────────────────────────────────────────────

    #[test]
    fn test_symplectic_leapfrog_position() {
        let sc = SymplecticConstraint::new(1e-6, 2);
        let p = sc.leapfrog_position_update([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 2.0, 0.01);
        // dt * p / m = 0.01 * 10 / 2 = 0.05
        assert!((p[0] - 0.05).abs() < 1e-10, "p[0]={}", p[0]);
    }

    #[test]
    fn test_symplectic_leapfrog_momentum() {
        let sc = SymplecticConstraint::new(1e-6, 2);
        let p = sc.leapfrog_momentum_update([0.0; 3], [5.0, 0.0, 0.0], 0.01);
        assert!((p[0] - 0.05).abs() < 1e-10, "p[0]={}", p[0]);
    }

    #[test]
    fn test_symplectic_forest_ruth_coefficients() {
        let (theta, xi) = SymplecticConstraint::forest_ruth_coefficients();
        // 2*theta + xi = 1
        assert!((2.0 * theta + xi - 1.0).abs() < 1e-10);
    }

    // ── EnergyDissipationConstraint ───────────────────────────────────────────

    #[test]
    fn test_dissipation_viscous() {
        let d = EnergyDissipationConstraint::viscous_dissipation(2.0, [3.0, 4.0, 0.0]);
        // 2.0 * 25.0 = 50
        assert!((d - 50.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn test_dissipation_accumulate_within_budget() {
        let mut dc = EnergyDissipationConstraint::new(1000.0, 0.0, 5000.0);
        let ok = dc.accumulate(100.0, 1.0);
        assert!(ok);
        assert!((dc.remaining_budget() - 4900.0).abs() < 1e-6);
    }

    // ── LyapunovConstraint ────────────────────────────────────────────────────

    #[test]
    fn test_lyapunov_evaluate_identity() {
        let lc = LyapunovConstraint::new(1.0, 100.0);
        let v = lc.evaluate([3.0, 4.0, 0.0]);
        // x^T I x = 9 + 16 = 25
        assert!((v - 25.0).abs() < 1e-10, "v={v}");
    }

    #[test]
    fn test_lyapunov_is_stable() {
        let mut lc = LyapunovConstraint::new(1.0, 100.0);
        lc.lyapunov_value = 10.0;
        // bound = -1 * 10 = -10; v_dot = -15 < -10 => stable
        assert!(lc.is_stable(-15.0));
        assert!(!lc.is_stable(-5.0));
    }

    // ── SolarIrradianceModel ──────────────────────────────────────────────────

    #[test]
    fn test_solar_power_at_zenith() {
        let sol = SolarIrradianceModel::new(1000.0, 10.0, 0.2, 0.9);
        // elevation = pi/2 => sin = 1.0
        let p = sol.available_power(std::f64::consts::FRAC_PI_2);
        assert!((p - 1800.0).abs() < 1e-4, "p={p}");
    }

    #[test]
    fn test_solar_power_below_horizon() {
        let sol = SolarIrradianceModel::new(1000.0, 10.0, 0.2, 0.9);
        assert!(sol.available_power(-0.1).abs() < 1e-10);
    }

    // ── WindEnergyModel ───────────────────────────────────────────────────────

    #[test]
    fn test_wind_energy_below_cut_in() {
        let wm = WindEnergyModel::new(1.225, 50.0, 0.4, 3.0, 25.0, 12.0, 2_000_000.0);
        assert!(wm.power_output(2.0).abs() < 1e-10);
    }

    #[test]
    fn test_wind_energy_rated() {
        let wm = WindEnergyModel::new(1.225, 50.0, 0.4, 3.0, 25.0, 12.0, 2_000_000.0);
        let p = wm.power_output(15.0);
        assert!((p - 2_000_000.0).abs() < 1e-4, "p={p}");
    }

    #[test]
    fn test_wind_betz_limit() {
        let wm = WindEnergyModel::new(1.225, 50.0, 0.4, 3.0, 25.0, 12.0, 2_000_000.0);
        let betz = wm.betz_limit_power(10.0);
        let theoretical = 0.5 * 1.225 * 50.0 * 0.593 * 1000.0;
        assert!(
            (betz - theoretical).abs() < 1.0,
            "betz={betz} expected={theoretical}"
        );
    }

    // ── EnergyManager integration ──────────────────────────────────────────────

    #[test]
    fn test_energy_manager_step() {
        let bat = BatteryConstraint::new(100_000.0, 0.1, 0.9, 10000.0, 10000.0, 0.05, 0.5);
        let tc = ThermalConstraint::new(300.0, 250.0, 450.0, 5000.0, 20.0, 293.0);
        let dc = EnergyDissipationConstraint::new(5000.0, 0.0, 1_000_000.0);
        let mut em = EnergyManager::new(bat, tc, dc, 50000.0);
        let ok = em.step(500.0, 100.0, 200.0, 30000.0, 20000.0, 0.01);
        assert!(ok);
    }
}
