//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::VecDeque;

use super::functions::soc_to_ocv;

/// Smart charging scheduler using LP-inspired greedy algorithm.
#[derive(Debug, Clone)]
pub struct SmartChargingScheduler {
    /// Departure time (hours from now)
    pub departure_hour: f64,
    /// Target SOC at departure
    pub target_soc: f64,
    /// Current SOC
    pub current_soc: f64,
    /// Pack capacity (kWh)
    pub pack_capacity_kwh: f64,
    /// Max charge power (W)
    pub max_power_w: f64,
    /// Grid tariff
    pub tariff: GridTariff,
}
impl SmartChargingScheduler {
    /// Create a scheduler.
    pub fn new(
        departure_hour: f64,
        target_soc: f64,
        current_soc: f64,
        pack_capacity_kwh: f64,
        max_power_w: f64,
    ) -> Self {
        Self {
            departure_hour,
            target_soc,
            current_soc,
            pack_capacity_kwh,
            max_power_w,
            tariff: GridTariff::two_rate(0.30, 0.10),
        }
    }
    /// Build an optimal charging schedule (cheapest hours first).
    pub fn build_schedule(&self, current_hour: f64) -> Vec<ChargingSlot> {
        let energy_needed_kwh =
            (self.target_soc - self.current_soc).max(0.0) * self.pack_capacity_kwh;
        let hours_available = (self.departure_hour - current_hour).max(0.0);
        let slot_duration = 0.5_f64;
        let n_slots = (hours_available / slot_duration) as usize;
        let mut slots: Vec<ChargingSlot> = (0..n_slots)
            .map(|i| {
                let start = current_hour + i as f64 * slot_duration;
                let price = self.tariff.price_at_hour(start);
                ChargingSlot {
                    start_hour: start,
                    end_hour: start + slot_duration,
                    power_w: self.max_power_w,
                    price,
                }
            })
            .collect();
        slots.sort_by(|a, b| {
            a.price
                .partial_cmp(&b.price)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut remaining_kwh = energy_needed_kwh;
        let mut schedule = Vec::new();
        for slot in slots {
            if remaining_kwh <= 0.0 {
                break;
            }
            let energy_this_slot = slot.power_w * slot_duration / 1000.0;
            let fraction = (remaining_kwh / energy_this_slot).min(1.0);
            let mut chosen_slot = slot.clone();
            chosen_slot.power_w *= fraction;
            remaining_kwh -= energy_this_slot * fraction;
            schedule.push(chosen_slot);
        }
        schedule.sort_by(|a, b| {
            a.start_hour
                .partial_cmp(&b.start_hour)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        schedule
    }
    /// Total cost of a schedule (€).
    pub fn schedule_cost(schedule: &[ChargingSlot]) -> f64 {
        schedule
            .iter()
            .map(|s| s.power_w * (s.end_hour - s.start_hour) / 1000.0 * s.price)
            .sum()
    }
}
/// Tracks energy consumption and estimates range from battery capacity.
#[derive(Debug, Clone)]
pub struct EnergyConsumption {
    /// Total energy consumed (kWh).
    pub energy_kwh: f64,
    /// Total distance travelled (km).
    pub distance_km: f64,
    /// Battery usable capacity (kWh).
    pub battery_capacity_kwh: f64,
    /// Current SOC (0–1).
    pub soc: f64,
    /// Rolling average consumption (kWh/km) — updated each step.
    pub rolling_consumption_kwh_per_km: f64,
    /// Smoothing factor for the rolling average (0–1, lower = more smoothing).
    pub smoothing: f64,
}
impl EnergyConsumption {
    /// Create an energy consumption tracker.
    ///
    /// * `battery_capacity_kwh` — usable battery energy (kWh).
    /// * `initial_soc`          — starting state of charge (0–1).
    pub fn new(battery_capacity_kwh: f64, initial_soc: f64) -> Self {
        Self {
            energy_kwh: 0.0,
            distance_km: 0.0,
            battery_capacity_kwh,
            soc: initial_soc.clamp(0.0, 1.0),
            rolling_consumption_kwh_per_km: 0.2,
            smoothing: 0.05,
        }
    }
    /// Record a step: `power_w` (W) consumed for `dt` seconds at `speed_ms` (m/s).
    pub fn record_step(&mut self, power_w: f64, dt: f64, speed_ms: f64) {
        let energy_step = power_w * dt / 3_600_000.0;
        let dist_step = speed_ms * dt / 1000.0;
        self.energy_kwh += energy_step.max(0.0);
        self.distance_km += dist_step;
        let soc_delta = energy_step / self.battery_capacity_kwh.max(1e-6);
        self.soc = (self.soc - soc_delta).clamp(0.0, 1.0);
        if dist_step > 1e-6 {
            let inst = energy_step / dist_step;
            self.rolling_consumption_kwh_per_km = (1.0 - self.smoothing)
                * self.rolling_consumption_kwh_per_km
                + self.smoothing * inst;
        }
    }
    /// Average consumption over the trip (kWh/km), or 0 if no distance.
    pub fn average_consumption_kwh_per_km(&self) -> f64 {
        if self.distance_km < 1e-6 {
            return 0.0;
        }
        self.energy_kwh / self.distance_km
    }
    /// Estimated remaining range (km) based on rolling consumption and SOC.
    pub fn remaining_range_km(&self) -> f64 {
        let remaining_kwh = self.soc * self.battery_capacity_kwh;
        remaining_kwh / self.rolling_consumption_kwh_per_km.max(1e-6)
    }
    /// Wh per km (same as `average_consumption` × 1000).
    pub fn wh_per_km(&self) -> f64 {
        self.average_consumption_kwh_per_km() * 1000.0
    }
}
/// Configuration of a battery pack.
#[derive(Debug, Clone)]
pub struct PackConfig {
    /// Number of cells in series
    pub cells_series: usize,
    /// Number of cells in parallel
    pub cells_parallel: usize,
    /// Chemistry
    pub chemistry: CellChemistry,
    /// Nominal cell capacity (Ah)
    pub cell_capacity_ah: f64,
    /// Thermal resistance (°C/W) of pack cooling
    pub thermal_resistance: f64,
    /// Coolant temperature (°C)
    pub coolant_temp: f64,
}
impl PackConfig {
    /// Create a standard 400V NMC pack.
    pub fn standard_400v_nmc() -> Self {
        Self {
            cells_series: 96,
            cells_parallel: 4,
            chemistry: CellChemistry::NMC,
            cell_capacity_ah: 50.0,
            thermal_resistance: 0.05,
            coolant_temp: 20.0,
        }
    }
    /// Total nominal pack voltage (V).
    pub fn nominal_voltage(&self) -> f64 {
        self.cells_series as f64 * self.chemistry.nominal_voltage()
    }
    /// Total capacity (kWh).
    pub fn total_capacity_kwh(&self) -> f64 {
        let ah = self.cell_capacity_ah * self.cells_parallel as f64;
        let v = self.nominal_voltage();
        ah * v / 1000.0
    }
}
/// PI current controller for FOC.
#[derive(Debug, Clone)]
pub struct PiController {
    /// Proportional gain
    pub kp: f64,
    /// Integral gain
    pub ki: f64,
    /// Integral accumulator
    integrator: f64,
    /// Anti-windup limit
    pub limit: f64,
}
impl PiController {
    /// Create a PI controller.
    pub fn new(kp: f64, ki: f64, limit: f64) -> Self {
        Self {
            kp,
            ki,
            integrator: 0.0,
            limit,
        }
    }
    /// Step the controller. Returns control output.
    pub fn step(&mut self, error: f64, dt_s: f64) -> f64 {
        self.integrator = (self.integrator + error * dt_s).clamp(-self.limit, self.limit);
        let out = self.kp * error + self.ki * self.integrator;
        out.clamp(-self.limit, self.limit)
    }
    /// Reset integrator.
    pub fn reset(&mut self) {
        self.integrator = 0.0;
    }
}
/// State-of-charge estimator using Coulomb counting with periodic open-circuit
/// voltage (OCV) correction.
#[derive(Debug, Clone)]
pub struct SocEstimator {
    /// Current SOC estimate (0–1).
    pub soc: f64,
    /// Nominal capacity (Ah).
    pub capacity_ah: f64,
    /// Coulomb counting accumulated charge (Ah).
    pub accumulated_ah: f64,
    /// OCV correction gain: weight given to OCV-based estimate vs Coulomb
    /// counter on each correction call (0–1).
    pub ocv_correction_gain: f64,
    /// Number of Coulomb-counting steps between OCV corrections.
    pub correction_interval: u32,
    /// Step counter for OCV correction scheduling.
    pub step_counter: u32,
}
impl SocEstimator {
    /// Create a SOC estimator starting at `initial_soc`.
    ///
    /// * `capacity_ah`     — cell/pack capacity in Ah.
    /// * `initial_soc`     — starting SOC (0–1).
    pub fn new(capacity_ah: f64, initial_soc: f64) -> Self {
        Self {
            soc: initial_soc.clamp(0.0, 1.0),
            capacity_ah,
            accumulated_ah: 0.0,
            ocv_correction_gain: 0.1,
            correction_interval: 100,
            step_counter: 0,
        }
    }
    /// Update SOC from current (A) for time `dt` (s).
    ///
    /// Positive current = discharge; negative = charge.
    pub fn update_coulomb(&mut self, current_a: f64, dt: f64) {
        let delta_ah = current_a * dt / 3600.0;
        self.accumulated_ah += delta_ah;
        self.soc -= delta_ah / self.capacity_ah.max(1e-6);
        self.soc = self.soc.clamp(0.0, 1.0);
        self.step_counter += 1;
    }
    /// Apply an OCV-based correction using the measured open-circuit voltage
    /// `v_ocv` (V) and a linear OCV–SOC model with `v_min` and `v_max`.
    pub fn correct_from_ocv(&mut self, v_ocv: f64, v_min: f64, v_max: f64) {
        let soc_ocv = ((v_ocv - v_min) / (v_max - v_min).max(1e-6)).clamp(0.0, 1.0);
        self.soc = (1.0 - self.ocv_correction_gain) * self.soc + self.ocv_correction_gain * soc_ocv;
        self.step_counter = 0;
    }
    /// Returns `true` if an OCV correction is due.
    pub fn correction_due(&self) -> bool {
        self.step_counter >= self.correction_interval
    }
    /// Reset the estimator to a given SOC (e.g., after a full charge).
    pub fn reset(&mut self, soc: f64) {
        self.soc = soc.clamp(0.0, 1.0);
        self.accumulated_ah = 0.0;
        self.step_counter = 0;
    }
    /// Remaining charge (Ah) available.
    pub fn remaining_ah(&self) -> f64 {
        self.soc * self.capacity_ah
    }
}
/// Range estimator using energy consumption model.
#[derive(Debug, Clone)]
pub struct RangeEstimator {
    /// Vehicle mass (kg)
    pub mass_kg: f64,
    /// Frontal area (m²)
    pub frontal_area: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Rolling resistance coefficient
    pub crr: f64,
    /// Motor efficiency map
    pub eta_map: EfficiencyMap,
    /// Pack capacity (kWh)
    pub pack_capacity_kwh: f64,
    /// Auxiliary load (W)
    pub aux_load_w: f64,
}
impl RangeEstimator {
    /// Air density (kg/m³) at standard conditions.
    const RHO_AIR: f64 = 1.225;
    /// Create a range estimator.
    pub fn new(mass_kg: f64, frontal_area: f64, cd: f64, pack_capacity_kwh: f64) -> Self {
        Self {
            mass_kg,
            frontal_area,
            cd,
            crr: 0.01,
            eta_map: EfficiencyMap::synthetic_pmsm(),
            pack_capacity_kwh,
            aux_load_w: 500.0,
        }
    }
    /// Estimate specific energy consumption (Wh/km) for a drive cycle.
    pub fn specific_consumption_wh_per_km(&self, cycle: &DriveCycle, _soc: f64) -> f64 {
        let mut energy_wh = 0.0;
        let n = cycle.speed.len();
        for i in 1..n {
            let v = cycle.speed[i];
            let dt = 1.0;
            let a = (cycle.speed[i] - cycle.speed[i.saturating_sub(1)]) / dt;
            let f_drag = 0.5 * Self::RHO_AIR * self.cd * self.frontal_area * v.powi(2);
            let f_roll = self.mass_kg * 9.81 * self.crr;
            let f_accel = self.mass_kg * a;
            let f_total = f_drag + f_roll + f_accel;
            let p_mech = f_total * v;
            let rpm = v * 60.0 / (2.0 * std::f64::consts::PI * 0.3);
            let torque = if v > 0.1 {
                p_mech / (rpm * std::f64::consts::TAU / 60.0).max(1.0)
            } else {
                0.0
            };
            let eta = self.eta_map.efficiency(rpm.abs(), torque.abs());
            let p_elec = if p_mech > 0.0 {
                p_mech / eta.max(0.1)
            } else {
                p_mech * eta
            };
            let p_total = p_elec + self.aux_load_w;
            energy_wh += p_total * dt / 3600.0;
        }
        let dist_km = cycle.total_distance_m() / 1000.0;
        if dist_km > 0.0 {
            energy_wh / dist_km
        } else {
            200.0
        }
    }
    /// Estimate remaining range (km) from current SOC.
    pub fn remaining_range_km(&self, current_soc: f64, cycle: &DriveCycle) -> f64 {
        let wh_per_km = self.specific_consumption_wh_per_km(cycle, current_soc);
        let available_kwh = self.pack_capacity_kwh * current_soc;
        available_kwh * 1000.0 / wh_per_km.max(1.0)
    }
}
/// A simple lithium-ion battery model.
#[derive(Debug, Clone)]
pub struct Battery {
    /// Total capacity (kWh).
    pub capacity_kwh: f64,
    /// State of charge (0–1).
    pub soc: f64,
    /// Nominal open-circuit voltage (V).
    pub voltage_nominal: f64,
    /// Internal resistance (Ω).
    pub internal_resistance: f64,
    /// Battery temperature (°C).
    pub temp_celsius: f64,
}
impl Battery {
    /// Create a typical 75 kWh EV battery.
    pub fn typical_75kwh() -> Self {
        Self {
            capacity_kwh: 75.0,
            soc: 1.0,
            voltage_nominal: 400.0,
            internal_resistance: 0.05,
            temp_celsius: 25.0,
        }
    }
}
/// Models the efficiency of converting deceleration kinetic energy into
/// battery charge (regenerative braking).
#[derive(Debug, Clone)]
pub struct RegenerativeEfficiency {
    /// Motor efficiency during generation (0–1).
    pub motor_eta: f64,
    /// Inverter efficiency during regen (0–1).
    pub inverter_eta: f64,
    /// Battery charge efficiency (0–1).
    pub battery_charge_eta: f64,
    /// Maximum regen torque (N·m).
    pub max_regen_torque: f64,
    /// Maximum regen power (W).
    pub max_regen_power: f64,
    /// SOC above which regen is reduced or disabled (0–1).
    pub soc_limit: f64,
}
impl RegenerativeEfficiency {
    /// Create a typical EV regen efficiency model.
    ///
    /// * `max_torque` — maximum regen torque (N·m).
    /// * `max_power`  — maximum regen power (W).
    pub fn new(max_torque: f64, max_power: f64) -> Self {
        Self {
            motor_eta: 0.94,
            inverter_eta: 0.97,
            battery_charge_eta: 0.98,
            max_regen_torque: max_torque,
            max_regen_power: max_power,
            soc_limit: 0.95,
        }
    }
    /// End-to-end regen chain efficiency.
    pub fn chain_efficiency(&self) -> f64 {
        self.motor_eta * self.inverter_eta * self.battery_charge_eta
    }
    /// Regen torque (N·m) and battery charge rate (W) for a given braking
    /// demand `brake_torque` (N·m, positive) at wheel speed `omega` (rad/s)
    /// and current SOC.
    ///
    /// Returns `(regen_torque_nm, charge_power_w)`.
    pub fn compute(&self, brake_torque: f64, omega: f64, soc: f64) -> (f64, f64) {
        if soc >= self.soc_limit || omega < 0.5 {
            return (0.0, 0.0);
        }
        let soc_factor = ((self.soc_limit - soc) / (self.soc_limit - 0.0)).clamp(0.0, 1.0);
        let regen = brake_torque.min(self.max_regen_torque) * soc_factor;
        let power_mech = regen * omega;
        let charge_power = (power_mech * self.chain_efficiency()).min(self.max_regen_power);
        (regen, charge_power)
    }
    /// Recoverable energy (J) from a deceleration event: vehicle of `mass` kg
    /// slowing from `v_initial` to `v_final` (m/s).
    pub fn recoverable_energy(&self, mass: f64, v_initial: f64, v_final: f64, soc: f64) -> f64 {
        if soc >= self.soc_limit {
            return 0.0;
        }
        let ke_delta = 0.5 * mass * (v_initial * v_initial - v_final * v_final).max(0.0);
        ke_delta * self.chain_efficiency()
    }
}
/// Three-phase current in αβ (Clarke) frame.
#[derive(Debug, Clone, Copy)]
pub struct AlphaBetaCurrent {
    /// α component (A)
    pub alpha: f64,
    /// β component (A)
    pub beta: f64,
}
/// State of health estimator using incremental capacity analysis.
#[derive(Debug, Clone)]
pub struct SohEstimator {
    /// Voltage samples (V)
    voltage_buf: VecDeque<f64>,
    /// Capacity samples (Ah)
    capacity_buf: VecDeque<f64>,
    /// Reference nominal capacity (Ah)
    nominal_capacity: f64,
    /// Current SOH estimate
    pub soh: f64,
    /// Calendar aging coefficient (per day)
    pub calendar_aging_coeff: f64,
    /// Cycle aging coefficient (per full cycle)
    pub cycle_aging_coeff: f64,
}
impl SohEstimator {
    /// Create a new SOH estimator.
    pub fn new(nominal_capacity_ah: f64) -> Self {
        Self {
            voltage_buf: VecDeque::with_capacity(1000),
            capacity_buf: VecDeque::with_capacity(1000),
            nominal_capacity: nominal_capacity_ah,
            soh: 1.0,
            calendar_aging_coeff: 2e-5,
            cycle_aging_coeff: 1e-4,
        }
    }
    /// Add a measurement sample for ICA (Incremental Capacity Analysis).
    pub fn add_sample(&mut self, voltage: f64, capacity_ah: f64) {
        if self.voltage_buf.len() >= 1000 {
            self.voltage_buf.pop_front();
            self.capacity_buf.pop_front();
        }
        self.voltage_buf.push_back(voltage);
        self.capacity_buf.push_back(capacity_ah);
    }
    /// Estimate capacity from charge curve via coulomb counting.
    pub fn estimate_capacity_from_charge(&self) -> f64 {
        if self.capacity_buf.len() < 2 {
            return self.nominal_capacity;
        }
        *self.capacity_buf.back().expect("deque should not be empty")
            - *self
                .capacity_buf
                .front()
                .expect("deque should not be empty")
    }
    /// Update SOH from measured capacity.
    pub fn update_soh_from_capacity(&mut self, measured_capacity_ah: f64) {
        self.soh = (measured_capacity_ah / self.nominal_capacity).clamp(0.0, 1.0);
    }
    /// Apply calendar aging model (days elapsed, temperature °C).
    pub fn apply_calendar_aging(&mut self, days: f64, temperature_c: f64) {
        let arrhenius = ((temperature_c - 25.0) / 10.0 * 0.693).exp();
        let degradation = self.calendar_aging_coeff * days * arrhenius;
        self.soh = (self.soh - degradation).max(0.0);
    }
    /// Apply cycle aging model (cycle count, depth of discharge fraction).
    pub fn apply_cycle_aging(&mut self, cycle_count: f64, dod: f64) {
        let dod_factor = dod.powi(2);
        let degradation = self.cycle_aging_coeff * cycle_count * dod_factor;
        self.soh = (self.soh - degradation).max(0.0);
    }
    /// End-of-life threshold check (SOH < 0.8 is typically EoL).
    pub fn is_end_of_life(&self) -> bool {
        self.soh < 0.8
    }
    /// Compute remaining useful life (cycles) at current degradation rate.
    pub fn remaining_useful_life_cycles(&self, cycles_so_far: f64, dod: f64) -> f64 {
        let dod_factor = dod.powi(2);
        let rate_per_cycle = self.cycle_aging_coeff * dod_factor;
        if rate_per_cycle <= 0.0 {
            return f64::INFINITY;
        }
        let remaining_degradation = self.soh - 0.8;
        if remaining_degradation <= 0.0 {
            return 0.0;
        }
        remaining_degradation / rate_per_cycle + cycles_so_far
    }
}
/// PMSM (Permanent Magnet Synchronous Motor) model: torque constant, back-EMF,
/// and copper losses.
#[derive(Debug, Clone)]
pub struct PermanentMagnetMotor {
    /// Torque constant Kt (N·m/A).
    pub kt: f64,
    /// Back-EMF constant Ke (V·s/rad).
    pub ke: f64,
    /// Phase resistance Rs (Ω).
    pub rs: f64,
    /// Phase inductance Ls (H).
    pub ls: f64,
    /// Number of pole pairs.
    pub pole_pairs: u32,
    /// Peak current (A).
    pub peak_current: f64,
    /// Rotor flux linkage λ (Wb).
    pub lambda: f64,
}
impl PermanentMagnetMotor {
    /// Create a typical automotive traction PMSM (e.g., 150 kW peak).
    pub fn automotive_traction() -> Self {
        Self {
            kt: 1.2,
            ke: 1.2,
            rs: 0.015,
            ls: 0.0003,
            pole_pairs: 4,
            peak_current: 600.0,
            lambda: 0.1,
        }
    }
    /// Torque produced for a given q-axis current Iq (A).
    pub fn torque_from_iq(&self, iq: f64) -> f64 {
        1.5 * self.pole_pairs as f64 * self.lambda * iq
    }
    /// Back-EMF voltage (V) at mechanical speed `omega_mech` (rad/s).
    pub fn back_emf(&self, omega_mech: f64) -> f64 {
        self.ke * omega_mech * self.pole_pairs as f64
    }
    /// Copper losses (W) for phase current `i` (A).
    pub fn copper_losses(&self, i: f64) -> f64 {
        3.0 * self.rs * i * i
    }
    /// Required Iq (A) for a demanded torque (N·m).
    pub fn iq_for_torque(&self, torque_nm: f64) -> f64 {
        let denom = 1.5 * self.pole_pairs as f64 * self.lambda;
        if denom.abs() < 1e-12 {
            0.0
        } else {
            torque_nm / denom
        }
    }
    /// Terminal voltage (V) including back-EMF and resistive drop.
    pub fn terminal_voltage(&self, iq: f64, omega_mech: f64) -> f64 {
        let omega_e = omega_mech * self.pole_pairs as f64;

        self.rs * iq + omega_e * self.ls * 0.0 + omega_e * self.lambda
    }
    /// Maximum torque at a given dc-bus voltage `vdc` (V) and speed (rad/s).
    pub fn max_torque_at_speed(&self, vdc: f64, omega_mech: f64) -> f64 {
        let v_bemf = self.back_emf(omega_mech);
        let v_available = (vdc - v_bemf).max(0.0);
        let i_max = (v_available / self.rs.max(1e-6)).min(self.peak_current);
        self.torque_from_iq(i_max)
    }
}
/// Current in dq (Park) frame.
#[derive(Debug, Clone, Copy)]
pub struct DqCurrent {
    /// d-axis current (A) — flux component
    pub id: f64,
    /// q-axis current (A) — torque component
    pub iq: f64,
}
/// Thermal model for a motor: tracks winding temperature using a lumped
/// first-order thermal model (R·C network).
#[derive(Debug, Clone)]
pub struct ThermalMotorModel {
    /// Thermal resistance from winding to coolant (K/W).
    pub r_thermal: f64,
    /// Thermal capacitance of the winding mass (J/K).
    pub c_thermal: f64,
    /// Current winding temperature (°C).
    pub temperature_c: f64,
    /// Coolant/ambient temperature (°C).
    pub coolant_temp_c: f64,
    /// Maximum safe winding temperature (°C).
    pub max_temp_c: f64,
    /// Derating factor applied to torque when temperature exceeds 80 % of max.
    pub derate_factor: f64,
}
impl ThermalMotorModel {
    /// Create a typical EV traction motor thermal model.
    ///
    /// * `coolant_temp_c` — steady-state coolant temperature (°C).
    pub fn new(coolant_temp_c: f64) -> Self {
        Self {
            r_thermal: 0.08,
            c_thermal: 1500.0,
            temperature_c: coolant_temp_c,
            coolant_temp_c,
            max_temp_c: 180.0,
            derate_factor: 0.7,
        }
    }
    /// Step the thermal model: apply `heat_input_w` (W) for `dt` seconds.
    ///
    /// Returns updated winding temperature (°C).
    pub fn step(&mut self, heat_input_w: f64, dt: f64) -> f64 {
        let q_in = heat_input_w;
        let q_out = (self.temperature_c - self.coolant_temp_c) / self.r_thermal;
        let d_temp = (q_in - q_out) * dt / self.c_thermal;
        self.temperature_c += d_temp;
        self.temperature_c
    }
    /// Torque derating factor (0–1) based on current temperature.
    ///
    /// Returns 1.0 below 80 % of max temp, linearly reducing to `derate_factor`
    /// at max temp.
    pub fn torque_derating(&self) -> f64 {
        let onset = self.max_temp_c * 0.8;
        if self.temperature_c <= onset {
            1.0
        } else {
            let excess = (self.temperature_c - onset) / (self.max_temp_c - onset).max(1.0);
            1.0 - excess * (1.0 - self.derate_factor)
        }
    }
    /// Returns `true` if winding temperature is at or above the safety limit.
    pub fn is_overtemperature(&self) -> bool {
        self.temperature_c >= self.max_temp_c
    }
    /// Steady-state temperature (°C) at a constant heat input (W).
    pub fn steady_state_temp(&self, heat_w: f64) -> f64 {
        self.coolant_temp_c + heat_w * self.r_thermal
    }
}
/// Pulse charger with depolarization pulses.
#[derive(Debug, Clone)]
pub struct PulseCharger {
    /// Charge pulse current (A)
    pub pulse_current: f64,
    /// Rest current (0 = off)
    pub rest_current: f64,
    /// Pulse duration (s)
    pub pulse_duration: f64,
    /// Rest duration (s)
    pub rest_duration: f64,
    /// Phase timer (s)
    timer: f64,
    /// Current phase
    in_pulse: bool,
    /// Energy delivered (Wh)
    pub energy_in: f64,
}
impl PulseCharger {
    /// Create a pulse charger.
    pub fn new(pulse_current: f64, pulse_duration_s: f64, rest_duration_s: f64) -> Self {
        Self {
            pulse_current,
            rest_current: 0.0,
            pulse_duration: pulse_duration_s,
            rest_duration: rest_duration_s,
            timer: 0.0,
            in_pulse: true,
            energy_in: 0.0,
        }
    }
    /// Step the pulse charger. Returns current (A).
    pub fn step(&mut self, voltage: f64, dt_s: f64) -> f64 {
        self.timer += dt_s;
        if self.in_pulse && self.timer >= self.pulse_duration {
            self.in_pulse = false;
            self.timer = 0.0;
        } else if !self.in_pulse && self.timer >= self.rest_duration {
            self.in_pulse = true;
            self.timer = 0.0;
        }
        let current = if self.in_pulse {
            self.pulse_current
        } else {
            self.rest_current
        };
        self.energy_in += current * voltage * dt_s / 3600.0;
        current
    }
}
/// Parameters governing battery step behaviour.
#[derive(Debug, Clone)]
pub struct BatteryParams {
    /// Coulombic efficiency during charging (0–1).
    pub charge_efficiency: f64,
    /// Discharge efficiency (0–1).
    pub discharge_efficiency: f64,
    /// Thermal capacity of the pack (J/K).
    pub thermal_capacity: f64,
    /// Cooling power (W).
    pub cooling_power: f64,
}
impl BatteryParams {
    /// Default battery parameters.
    pub fn default_params() -> Self {
        Self {
            charge_efficiency: 0.95,
            discharge_efficiency: 0.97,
            thermal_capacity: 50_000.0,
            cooling_power: 2_000.0,
        }
    }
}
/// Driving cycle sample.
#[derive(Debug, Clone)]
pub struct DriveCycle {
    /// Time samples (s)
    pub time: Vec<f64>,
    /// Speed profile (m/s)
    pub speed: Vec<f64>,
}
impl DriveCycle {
    /// Create a synthetic WLTP-like drive cycle.
    pub fn wltp_synthetic() -> Self {
        let mut time = Vec::new();
        let mut speed = Vec::new();
        let mut t = 0.0_f64;
        let mut v = 0.0_f64;
        let phases = [
            (589.0, 46.5 / 3.6),
            (433.0, 76.6 / 3.6),
            (455.0, 97.4 / 3.6),
            (323.0, 131.3 / 3.6),
        ];
        for (duration, v_target) in phases.iter() {
            let n = (*duration / 1.0) as usize;
            for _ in 0..n {
                let dv = (v_target - v) * 0.05;
                v += dv;
                time.push(t);
                speed.push(v);
                t += 1.0;
            }
        }
        Self { time, speed }
    }
    /// Compute total distance (m).
    pub fn total_distance_m(&self) -> f64 {
        self.speed.iter().sum::<f64>() * 1.0
    }
}
/// A charging session slot.
#[derive(Debug, Clone)]
pub struct ChargingSlot {
    /// Start time (hours from now)
    pub start_hour: f64,
    /// End time (hours from now)
    pub end_hour: f64,
    /// Power level (W)
    pub power_w: f64,
    /// Price at this slot (€/kWh)
    pub price: f64,
}
/// Single or dual-motor EV drivetrain with torque vectoring.
#[derive(Debug, Clone)]
pub struct EvDrivetrain {
    /// Front motor controller (if present).
    pub front_motor: Option<MotorController>,
    /// Rear motor controller (if present).
    pub rear_motor: Option<MotorController>,
    /// Front axle torque ratio (0 = rear only, 1 = front only, 0.5 = equal).
    pub front_torque_ratio: f64,
    /// Final drive ratio (motor speed / wheel speed).
    pub final_drive_ratio: f64,
    /// Drivetrain mechanical efficiency (0–1).
    pub drivetrain_eta: f64,
    /// Current wheel speed (rad/s).
    pub wheel_speed_rads: f64,
    /// Total torque delivered to wheels this step (N·m).
    pub wheel_torque: f64,
    /// Torque vectoring gain: extra yaw correction torque per deg/s yaw error (N·m / (rad/s)).
    pub tv_gain: f64,
}
impl EvDrivetrain {
    /// Create a single rear-motor drivetrain.
    ///
    /// * `max_torque`       — peak motor torque (N·m).
    /// * `max_speed_rpm`    — peak motor speed (rpm).
    /// * `final_drive`      — gear/differential ratio.
    pub fn single_rear(max_torque: f64, max_speed_rpm: f64, final_drive: f64) -> Self {
        Self {
            front_motor: None,
            rear_motor: Some(MotorController::new(max_torque, max_speed_rpm)),
            front_torque_ratio: 0.0,
            final_drive_ratio: final_drive,
            drivetrain_eta: 0.95,
            wheel_speed_rads: 0.0,
            wheel_torque: 0.0,
            tv_gain: 0.0,
        }
    }
    /// Create a dual-motor (AWD) drivetrain with torque vectoring.
    ///
    /// * `front_torque` / `rear_torque` — peak motor torques (N·m).
    /// * `front_ratio`                  — fraction of total torque to front axle.
    /// * `final_drive`                  — combined drive ratio.
    /// * `tv_gain`                      — torque vectoring gain (N·m/(rad/s)).
    pub fn dual_motor_awd(
        front_torque: f64,
        front_speed_rpm: f64,
        rear_torque: f64,
        rear_speed_rpm: f64,
        front_ratio: f64,
        final_drive: f64,
        tv_gain: f64,
    ) -> Self {
        Self {
            front_motor: Some(MotorController::new(front_torque, front_speed_rpm)),
            rear_motor: Some(MotorController::new(rear_torque, rear_speed_rpm)),
            front_torque_ratio: front_ratio.clamp(0.0, 1.0),
            final_drive_ratio: final_drive,
            drivetrain_eta: 0.93,
            wheel_speed_rads: 0.0,
            wheel_torque: 0.0,
            tv_gain,
        }
    }
    /// Step the drivetrain: apply `driver_torque_demand` (N·m at wheels) for `dt` (s).
    ///
    /// `yaw_rate_error` (rad/s) is used for torque vectoring (0 if not needed).
    ///
    /// Returns total torque at wheels (N·m).
    pub fn step(&mut self, driver_torque_demand: f64, dt: f64, yaw_rate_error: f64) -> f64 {
        let motor_speed_rpm =
            self.wheel_speed_rads * self.final_drive_ratio * 60.0 / (2.0 * std::f64::consts::PI);
        let tv_correction = self.tv_gain * yaw_rate_error;
        let front_demand = driver_torque_demand * self.front_torque_ratio + tv_correction * 0.5;
        let rear_demand =
            driver_torque_demand * (1.0 - self.front_torque_ratio) - tv_correction * 0.5;
        let front_motor_demand = front_demand / self.final_drive_ratio.max(1e-6);
        let rear_motor_demand = rear_demand / self.final_drive_ratio.max(1e-6);
        let mut total_motor_torque = 0.0;
        if let Some(ref mut f) = self.front_motor {
            f.update_speed(motor_speed_rpm);
            total_motor_torque += f.step(front_motor_demand, dt);
        }
        if let Some(ref mut r) = self.rear_motor {
            r.update_speed(motor_speed_rpm);
            total_motor_torque += r.step(rear_motor_demand, dt);
        }
        self.wheel_torque = total_motor_torque * self.final_drive_ratio * self.drivetrain_eta;
        self.wheel_torque
    }
    /// Update wheel speed from vehicle speed and wheel radius.
    pub fn update_wheel_speed(&mut self, vehicle_speed_ms: f64, wheel_radius: f64) {
        self.wheel_speed_rads = vehicle_speed_ms / wheel_radius.max(1e-3);
    }
    /// Total electrical power draw (W).
    pub fn electrical_power_w(&self) -> f64 {
        let mut total = 0.0;
        if let Some(ref f) = self.front_motor {
            total += f.electrical_power();
        }
        if let Some(ref r) = self.rear_motor {
            total += r.electrical_power();
        }
        total
    }
    /// Returns `true` if both front and rear motors are present (AWD).
    pub fn is_awd(&self) -> bool {
        self.front_motor.is_some() && self.rear_motor.is_some()
    }
}
/// CC/CV charger state machine.
#[derive(Debug, Clone)]
pub struct CcCvCharger {
    /// Target charge current (A) for CC phase
    pub cc_current: f64,
    /// Target voltage (V) for CV phase
    pub cv_voltage: f64,
    /// Current at which to switch from CV to idle (A)
    pub termination_current: f64,
    /// Current charging mode
    pub mode: ChargingMode,
    /// Accumulated charge energy (Wh)
    pub energy_in: f64,
    /// Accumulated charge time (s)
    pub charge_time: f64,
}
impl CcCvCharger {
    /// Create a standard CC/CV charger.
    pub fn new(cc_current: f64, cv_voltage: f64) -> Self {
        Self {
            cc_current,
            cv_voltage,
            termination_current: cc_current * 0.05,
            mode: ChargingMode::Idle,
            energy_in: 0.0,
            charge_time: 0.0,
        }
    }
    /// Start charging.
    pub fn start(&mut self) {
        self.mode = ChargingMode::ConstantCurrent;
    }
    /// Compute charger output given pack state. Returns (current_A, voltage_V).
    pub fn step(&mut self, pack: &BatteryPack, dt_s: f64) -> (f64, f64) {
        let pack_v = pack.pack_voltage();
        match self.mode {
            ChargingMode::ConstantCurrent => {
                let i = self.cc_current;
                if pack_v >= self.cv_voltage {
                    self.mode = ChargingMode::ConstantVoltage;
                }
                let p = i * pack_v;
                self.energy_in += p * dt_s / 3600.0;
                self.charge_time += dt_s;
                (i, pack_v)
            }
            ChargingMode::ConstantVoltage => {
                let r_pack = pack.cells.iter().map(|c| c.r_internal).sum::<f64>()
                    / pack.config.cells_parallel as f64;
                let i = (self.cv_voltage - pack_v) / r_pack.max(0.001);
                let i = i.clamp(0.0, self.cc_current);
                if i < self.termination_current {
                    self.mode = ChargingMode::Idle;
                }
                let p = i * pack_v;
                self.energy_in += p * dt_s / 3600.0;
                self.charge_time += dt_s;
                (i, self.cv_voltage)
            }
            ChargingMode::Idle => (0.0, pack_v),
            _ => (0.0, pack_v),
        }
    }
    /// Charging efficiency (energy out / energy in).
    pub fn efficiency(&self, energy_stored_wh: f64) -> f64 {
        if self.energy_in <= 0.0 {
            return 0.0;
        }
        energy_stored_wh / self.energy_in
    }
}
/// Road segment for routing.
#[derive(Debug, Clone)]
pub struct RoadSegment {
    /// Segment index
    pub id: usize,
    /// Length (m)
    pub length_m: f64,
    /// Average speed (m/s)
    pub avg_speed_ms: f64,
    /// Gradient (positive = uphill)
    pub gradient_rad: f64,
    /// Road surface rolling resistance multiplier
    pub surface_factor: f64,
}
/// Induction motor model: torque-speed curve from slip frequency.
#[derive(Debug, Clone)]
pub struct InductionMotor {
    /// Number of pole pairs.
    pub pole_pairs: u32,
    /// Supply frequency (Hz).
    pub supply_frequency: f64,
    /// Rated slip at full torque (dimensionless, 0–1).
    pub rated_slip: f64,
    /// Rated torque (N·m).
    pub rated_torque: f64,
    /// Rotor resistance (Ω).
    pub rotor_resistance: f64,
    /// Stator resistance (Ω).
    pub stator_resistance: f64,
    /// Magnetising inductance (H).
    pub lm: f64,
    /// Rotor leakage inductance (H).
    pub lr: f64,
}
impl InductionMotor {
    /// Create a standard 4-pole 50 Hz induction motor.
    ///
    /// * `rated_torque` — rated torque at full load (N·m).
    pub fn new_50hz_4pole(rated_torque: f64) -> Self {
        Self {
            pole_pairs: 2,
            supply_frequency: 50.0,
            rated_slip: 0.05,
            rated_torque,
            rotor_resistance: 0.2,
            stator_resistance: 0.3,
            lm: 0.1,
            lr: 0.005,
        }
    }
    /// Synchronous speed (rad/s mechanical).
    pub fn synchronous_speed_rads(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.supply_frequency / self.pole_pairs as f64
    }
    /// Slip frequency for a given mechanical speed (rad/s).
    pub fn slip_frequency(&self, omega_mech: f64) -> f64 {
        let omega_sync = self.synchronous_speed_rads();
        (omega_sync - omega_mech) / omega_sync.max(1e-9)
    }
    /// Torque (N·m) using the Kloss formula for given slip `s` (−1 to 1).
    ///
    /// Returns positive torque for motoring (s > 0) and negative for generating (s < 0).
    pub fn torque_from_slip(&self, slip: f64) -> f64 {
        if slip.abs() < 1e-9 {
            return 0.0;
        }
        let ratio = slip / self.rated_slip;
        let linear = self.rated_torque * ratio;
        self.rated_torque
            * (ratio * 2.0 / (1.0 + ratio * ratio))
                .max(-self.rated_torque)
                .min(self.rated_torque)
            + linear * 0.0
    }
    /// Torque at a given mechanical speed (rad/s).
    pub fn torque_at_speed(&self, omega_mech: f64) -> f64 {
        let s = self.slip_frequency(omega_mech);
        self.torque_from_slip(s)
    }
    /// Copper losses (W) at given slip.
    pub fn copper_losses(&self, slip: f64, current_a: f64) -> f64 {
        (self.stator_resistance + self.rotor_resistance * slip * slip) * current_a * current_a
    }
}
/// Regenerative braking controller.
#[derive(Debug, Clone)]
pub struct RegenerativeBraking {
    /// Max regenerative torque (Nm)
    pub max_regen_torque: f64,
    /// Max charging power during regen (W)
    pub max_regen_power: f64,
    /// Blending factor (0 = mechanical only, 1 = full regen)
    pub blend_factor: f64,
    /// Energy recovered (Wh)
    pub energy_recovered: f64,
    /// Motor efficiency map
    pub efficiency_map: EfficiencyMap,
}
impl RegenerativeBraking {
    /// Create a regen controller.
    pub fn new(max_regen_torque: f64, max_regen_power: f64) -> Self {
        Self {
            max_regen_torque,
            max_regen_power,
            blend_factor: 0.7,
            energy_recovered: 0.0,
            efficiency_map: EfficiencyMap::synthetic_pmsm(),
        }
    }
    /// Compute regen torque given braking demand (Nm) and motor speed (rad/s).
    ///
    /// Returns (regen_torque_Nm, mechanical_brake_torque_Nm, recovered_power_W).
    pub fn compute_regen(
        &mut self,
        brake_demand_nm: f64,
        omega_rad_s: f64,
        soc: f64,
        dt_s: f64,
    ) -> (f64, f64, f64) {
        if omega_rad_s < 1.0 || soc >= 0.98 {
            return (0.0, brake_demand_nm, 0.0);
        }
        let rpm = omega_rad_s * 60.0 / std::f64::consts::TAU;
        let power_limit = (self.max_regen_power * (1.0 - soc / 0.98)).min(self.max_regen_power);
        let torque_from_power = power_limit / omega_rad_s;
        let regen_torque = brake_demand_nm
            .min(self.max_regen_torque)
            .min(torque_from_power)
            * self.blend_factor;
        let regen_torque = regen_torque.max(0.0);
        let mech_torque = brake_demand_nm - regen_torque;
        let eta = self.efficiency_map.efficiency(rpm, regen_torque);
        let recovered_power = regen_torque * omega_rad_s * eta;
        self.energy_recovered += recovered_power * dt_s / 3600.0;
        (regen_torque, mech_torque, recovered_power)
    }
}

/// Charging mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChargingMode {
    /// Constant current phase
    ConstantCurrent,
    /// Constant voltage phase
    ConstantVoltage,
    /// Pulse charging
    Pulse,
    /// Trickle charge
    Trickle,
    /// Idle / not charging
    Idle,
}
/// Chemistry type of a battery cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellChemistry {
    /// Lithium nickel manganese cobalt oxide
    NMC,
    /// Lithium iron phosphate
    LFP,
    /// Lithium nickel cobalt aluminium oxide
    NCA,
    /// Lithium titanate oxide
    LTO,
    /// Lithium manganese oxide
    LMO,
}
impl CellChemistry {
    /// Nominal cell voltage (V).
    pub fn nominal_voltage(&self) -> f64 {
        match self {
            CellChemistry::NMC => 3.6,
            CellChemistry::LFP => 3.2,
            CellChemistry::NCA => 3.6,
            CellChemistry::LTO => 2.3,
            CellChemistry::LMO => 3.7,
        }
    }
    /// Maximum charge voltage (V).
    pub fn max_voltage(&self) -> f64 {
        match self {
            CellChemistry::NMC => 4.2,
            CellChemistry::LFP => 3.65,
            CellChemistry::NCA => 4.2,
            CellChemistry::LTO => 2.85,
            CellChemistry::LMO => 4.2,
        }
    }
    /// Minimum discharge voltage (V).
    pub fn min_voltage(&self) -> f64 {
        match self {
            CellChemistry::NMC => 2.5,
            CellChemistry::LFP => 2.5,
            CellChemistry::NCA => 2.5,
            CellChemistry::LTO => 1.5,
            CellChemistry::LMO => 3.0,
        }
    }
    /// Typical thermal runaway onset temperature (°C).
    pub fn thermal_runaway_temp(&self) -> f64 {
        match self {
            CellChemistry::NMC => 180.0,
            CellChemistry::LFP => 270.0,
            CellChemistry::NCA => 150.0,
            CellChemistry::LTO => 300.0,
            CellChemistry::LMO => 250.0,
        }
    }
}
/// Motor controller: maps a torque command to actual delivered torque.
#[derive(Debug, Clone)]
pub struct MotorController {
    /// Maximum commandable torque (N·m).
    pub max_torque: f64,
    /// Maximum motor speed (rpm).
    pub max_speed_rpm: f64,
    /// Efficiency map used to derate/correct torque.
    pub efficiency_map: EfficiencyMap,
    /// Current speed estimate (rpm).
    pub speed_rpm: f64,
    /// Last delivered torque (N·m).
    pub actual_torque: f64,
    /// Torque filter time constant (s) — first-order lag to simulate actuator lag.
    pub torque_filter_tc: f64,
}
impl MotorController {
    /// Create a motor controller with a synthetic PMSM efficiency map.
    ///
    /// * `max_torque`   — peak torque (N·m).
    /// * `max_speed`    — peak speed (rpm).
    pub fn new(max_torque: f64, max_speed_rpm: f64) -> Self {
        Self {
            max_torque,
            max_speed_rpm,
            efficiency_map: EfficiencyMap::synthetic_pmsm(),
            speed_rpm: 0.0,
            actual_torque: 0.0,
            torque_filter_tc: 0.02,
        }
    }
    /// Step the controller: apply `commanded_torque` (N·m) for `dt` seconds.
    ///
    /// Returns actual torque delivered this step (N·m).
    pub fn step(&mut self, commanded_torque: f64, dt: f64) -> f64 {
        let clamped = commanded_torque.clamp(-self.max_torque, self.max_torque);
        let base_speed = self.max_speed_rpm * 0.3;
        let torque_limit = if self.speed_rpm > base_speed {
            clamped * (base_speed / self.speed_rpm.max(1.0)).min(1.0)
        } else {
            clamped
        };
        let eta = self
            .efficiency_map
            .efficiency(self.speed_rpm, torque_limit.abs());
        let target = torque_limit * eta;
        let alpha = dt / (self.torque_filter_tc + dt);
        self.actual_torque += alpha * (target - self.actual_torque);
        self.actual_torque
    }
    /// Update the internal speed estimate (rpm).
    pub fn update_speed(&mut self, speed_rpm: f64) {
        self.speed_rpm = speed_rpm.max(0.0).min(self.max_speed_rpm);
    }
    /// Electrical power draw for the delivered torque (W).
    ///
    /// Positive = consuming power from battery.
    pub fn electrical_power(&self) -> f64 {
        let omega = self.speed_rpm * (2.0 * std::f64::consts::PI / 60.0);
        let eta = self
            .efficiency_map
            .efficiency(self.speed_rpm, self.actual_torque.abs())
            .max(0.01);
        if self.actual_torque >= 0.0 {
            self.actual_torque * omega / eta
        } else {
            self.actual_torque * omega * eta
        }
    }
}
/// Single battery cell state.
#[derive(Debug, Clone)]
pub struct BatteryCell {
    /// State of charge \[0,1\]
    pub soc: f64,
    /// State of health \[0,1\] (capacity relative to nominal)
    pub soh: f64,
    /// Cell temperature (°C)
    pub temperature: f64,
    /// Internal resistance (Ω)
    pub r_internal: f64,
    /// Capacity (Ah)
    pub capacity_ah: f64,
    /// Cumulative charge throughput (Ah)
    pub charge_throughput: f64,
    /// Chemistry
    pub chemistry: CellChemistry,
    /// Open-circuit voltage (V)
    pub ocv: f64,
    /// Terminal voltage (V)
    pub terminal_voltage: f64,
}
impl BatteryCell {
    /// Create a new cell with given chemistry and capacity.
    pub fn new(chemistry: CellChemistry, capacity_ah: f64) -> Self {
        let ocv = chemistry.nominal_voltage();
        Self {
            soc: 1.0,
            soh: 1.0,
            temperature: 25.0,
            r_internal: 0.002,
            capacity_ah,
            charge_throughput: 0.0,
            chemistry,
            ocv,
            terminal_voltage: ocv,
        }
    }
    /// Compute open-circuit voltage from SOC using a polynomial model.
    ///
    /// Uses a 5th-order Bernstein polynomial fit typical for NMC/LFP cells.
    pub fn compute_ocv(&self) -> f64 {
        let s = self.soc.clamp(0.0, 1.0);
        let v_min = self.chemistry.min_voltage();
        let v_max = self.chemistry.max_voltage();
        let x = (s - 0.5) * 8.0;
        let sigmoid = 1.0 / (1.0 + (-x).exp());
        v_min + (v_max - v_min) * sigmoid
    }
    /// Update cell SOC given current (A) and time step (s).
    ///
    /// Positive current = discharge, negative = charge.
    pub fn update_soc(&mut self, current_a: f64, dt_s: f64) {
        let dq = current_a * dt_s / 3600.0;
        let effective_capacity = self.capacity_ah * self.soh;
        self.soc -= dq / effective_capacity;
        self.soc = self.soc.clamp(0.0, 1.0);
        self.charge_throughput += current_a.abs() * dt_s / 3600.0;
        self.ocv = self.compute_ocv();
        self.terminal_voltage = self.ocv - current_a * self.r_internal;
    }
    /// Estimate internal resistance as function of temperature and SOC.
    pub fn update_r_internal(&mut self) {
        let temp_factor = (-(self.temperature - 25.0) * 0.015).exp();
        let soc_factor = 1.0 + 0.3 * (1.0 - self.soc).powi(2);
        self.r_internal = 0.002 * temp_factor * soc_factor;
    }
    /// Compute heat generation rate (W).
    pub fn heat_generation_rate(&self, current_a: f64) -> f64 {
        current_a.powi(2) * self.r_internal
    }
}
/// 2D lookup table for motor efficiency.
#[derive(Debug, Clone)]
pub struct EfficiencyMap {
    /// Speed breakpoints (rpm)
    pub speed_rpm: Vec<f64>,
    /// Torque breakpoints (Nm)
    pub torque_nm: Vec<f64>,
    /// Efficiency table \[speed\]\[torque\]
    pub eta: Vec<Vec<f64>>,
}
impl EfficiencyMap {
    /// Create a synthetic efficiency map for a typical EV motor.
    pub fn synthetic_pmsm() -> Self {
        let speeds: Vec<f64> = (0..=20).map(|i| i as f64 * 500.0).collect();
        let torques: Vec<f64> = (0..=20).map(|i| i as f64 * 10.0).collect();
        let mut eta = vec![vec![0.0f64; 21]; 21];
        for (si, &s) in speeds.iter().enumerate() {
            for (ti, &t) in torques.iter().enumerate() {
                let p = s * t * std::f64::consts::PI / 30.0;
                let s_norm = s / 5000.0;
                let t_norm = t / 150.0;
                let eff = if p < 100.0 {
                    0.7
                } else {
                    0.96 - 0.04 * (s_norm - 0.5).powi(2) - 0.03 * (t_norm - 0.4).powi(2)
                };
                eta[si][ti] = eff.clamp(0.5, 0.98);
            }
        }
        Self {
            speed_rpm: speeds,
            torque_nm: torques,
            eta,
        }
    }
    /// Bilinear interpolation to get efficiency at (speed, torque).
    pub fn efficiency(&self, speed_rpm: f64, torque_nm: f64) -> f64 {
        let s = speed_rpm.clamp(
            *self
                .speed_rpm
                .first()
                .expect("collection should not be empty"),
            *self
                .speed_rpm
                .last()
                .expect("collection should not be empty"),
        );
        let t = torque_nm.clamp(
            *self
                .torque_nm
                .first()
                .expect("collection should not be empty"),
            *self
                .torque_nm
                .last()
                .expect("collection should not be empty"),
        );
        let si = self
            .speed_rpm
            .partition_point(|&x| x < s)
            .saturating_sub(1)
            .min(self.speed_rpm.len() - 2);
        let ti = self
            .torque_nm
            .partition_point(|&x| x < t)
            .saturating_sub(1)
            .min(self.torque_nm.len() - 2);
        let ds = (s - self.speed_rpm[si]) / (self.speed_rpm[si + 1] - self.speed_rpm[si] + 1e-12);
        let dt = (t - self.torque_nm[ti]) / (self.torque_nm[ti + 1] - self.torque_nm[ti] + 1e-12);
        let e00 = self.eta[si][ti];
        let e10 = self.eta[si + 1][ti];
        let e01 = self.eta[si][ti + 1];
        let e11 = self.eta[si + 1][ti + 1];
        (e00 * (1.0 - ds) * (1.0 - dt)
            + e10 * ds * (1.0 - dt)
            + e01 * (1.0 - ds) * dt
            + e11 * ds * dt)
            .clamp(0.0, 1.0)
    }
}
/// Grid tariff schedule (time-of-use pricing).
#[derive(Debug, Clone)]
pub struct GridTariff {
    /// Hour-of-day breakpoints \[0..24\]
    pub hours: Vec<f64>,
    /// Price (€/kWh) for each interval
    pub price: Vec<f64>,
}
impl GridTariff {
    /// Create a simple two-rate tariff (peak/off-peak).
    pub fn two_rate(peak_price: f64, offpeak_price: f64) -> Self {
        Self {
            hours: vec![0.0, 7.0, 23.0, 24.0],
            price: vec![offpeak_price, peak_price, offpeak_price],
        }
    }
    /// Get price for given hour of day.
    pub fn price_at_hour(&self, hour: f64) -> f64 {
        let h = hour % 24.0;
        for i in 0..self.hours.len().saturating_sub(1) {
            if h >= self.hours[i] && h < self.hours[i + 1] {
                return self.price[i];
            }
        }
        *self.price.last().unwrap_or(&0.1)
    }
}
/// Lightweight battery pack model for range and energy calculations.
///
/// Tracks capacity, nominal voltage, state-of-charge (SOC), state-of-health
/// Lightweight battery pack model for range and energy calculations.
#[derive(Debug, Clone)]
pub struct SimpleBatteryPack {
    /// Total pack capacity (kWh).
    pub capacity_kwh: f64,
    /// Nominal pack voltage (V).
    pub voltage_nominal: f64,
    /// State of charge in `[0, 1]`.
    pub soc: f64,
    /// State of health in `[0, 1]` (1 = new).
    pub soh: f64,
    /// Pack temperature (°C).
    pub temperature: f64,
}
impl SimpleBatteryPack {
    /// Construct a new pack.
    ///
    /// * `cap`     – total capacity (kWh)
    /// * `voltage` – nominal voltage (V)
    /// * `soc`     – initial state of charge in `[0, 1]`
    pub fn new(cap: f64, voltage: f64, soc: f64) -> Self {
        Self {
            capacity_kwh: cap,
            voltage_nominal: voltage,
            soc: soc.clamp(0.0, 1.0),
            soh: 1.0,
            temperature: 25.0,
        }
    }
    /// Available energy (kWh) = capacity × SOH × SOC.
    pub fn available_energy_kwh(&self) -> f64 {
        self.capacity_kwh * self.soh * self.soc
    }
    /// Open-circuit voltage estimated from SOC using `soc_to_ocv`.
    pub fn ocv_voltage(&self) -> f64 {
        soc_to_ocv(self.soc) * self.voltage_nominal
    }
    /// Internal resistance estimate (Ω) — increases with age and low temperature.
    ///
    /// A simple model: R_int = R_0 * (2 - SOH) * (1 + 0.01 * (25 - T))
    pub fn internal_resistance(&self) -> f64 {
        let r0 = 0.1;
        let age_factor = 2.0 - self.soh;
        let temp_factor = 1.0 + 0.01 * (25.0 - self.temperature);
        r0 * age_factor * temp_factor.max(0.1)
    }
}
/// Battery management system state.
#[derive(Debug, Clone)]
pub struct BatteryPack {
    /// Individual cells (flattened, row = series, col = parallel)
    pub cells: Vec<BatteryCell>,
    /// Configuration
    pub config: PackConfig,
    /// Pack-level SOC (energy-weighted)
    pub pack_soc: f64,
    /// Pack-level SOH
    pub pack_soh: f64,
    /// Max allowed continuous discharge current (A)
    pub max_discharge_current: f64,
    /// Max allowed continuous charge current (A)
    pub max_charge_current: f64,
    /// Current balancing state
    pub balancing_active: Vec<bool>,
    /// Total cycle count
    pub cycle_count: f64,
}
impl BatteryPack {
    /// Create a new battery pack from configuration.
    pub fn new(config: PackConfig) -> Self {
        let n = config.cells_series * config.cells_parallel;
        let cells: Vec<BatteryCell> = (0..n)
            .map(|_| BatteryCell::new(config.chemistry, config.cell_capacity_ah))
            .collect();
        let balancing_active = vec![false; n];
        Self {
            cells,
            pack_soc: 1.0,
            pack_soh: 1.0,
            max_discharge_current: 300.0,
            max_charge_current: 150.0,
            balancing_active,
            cycle_count: 0.0,
            config,
        }
    }
    /// Compute pack-level SOC as mean of cell SOCs.
    pub fn compute_pack_soc(&self) -> f64 {
        self.cells.iter().map(|c| c.soc).sum::<f64>() / self.cells.len() as f64
    }
    /// Get minimum SOC among all cells.
    pub fn min_cell_soc(&self) -> f64 {
        self.cells.iter().map(|c| c.soc).fold(f64::MAX, f64::min)
    }
    /// Get maximum SOC among all cells.
    pub fn max_cell_soc(&self) -> f64 {
        self.cells.iter().map(|c| c.soc).fold(f64::MIN, f64::max)
    }
    /// Get maximum cell temperature (°C).
    pub fn max_temperature(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| c.temperature)
            .fold(f64::MIN, f64::max)
    }
    /// Get minimum terminal voltage in series string (V).
    pub fn pack_voltage(&self) -> f64 {
        let s = self.config.cells_series;
        let p = self.config.cells_parallel;
        let mut total = 0.0;
        for col in 0..p {
            let mut string_v = 0.0;
            for row in 0..s {
                string_v += self.cells[row * p + col].terminal_voltage;
            }
            total += string_v;
        }
        total / p as f64
    }
    /// Perform passive cell balancing: bleed excess charge from high-SOC cells.
    ///
    /// Balances cells within `tolerance` SOC units of the minimum SOC.
    pub fn passive_balance(&mut self, tolerance: f64) {
        let min_soc = self.min_cell_soc();
        for (i, cell) in self.cells.iter_mut().enumerate() {
            self.balancing_active[i] = cell.soc > min_soc + tolerance;
        }
    }
    /// Apply balancing currents (small discharge on high-SOC cells).
    pub fn apply_balancing(&mut self, dt_s: f64) {
        let balance_current = 0.5;
        for (i, cell) in self.cells.iter_mut().enumerate() {
            if self.balancing_active[i] {
                cell.update_soc(balance_current, dt_s);
            }
        }
    }
    /// Update all cells given pack current (A) and time step (s).
    pub fn update(&mut self, pack_current_a: f64, dt_s: f64) {
        let cell_current = pack_current_a / self.config.cells_parallel as f64;
        for cell in &mut self.cells {
            cell.update_r_internal();
            cell.update_soc(cell_current, dt_s);
            let q_gen = cell.heat_generation_rate(cell_current);
            let dt_temp = (q_gen
                - (cell.temperature - self.config.coolant_temp) / self.config.thermal_resistance)
                * dt_s
                / 1000.0;
            cell.temperature += dt_temp;
        }
        self.pack_soc = self.compute_pack_soc();
        self.cycle_count += pack_current_a.abs() * dt_s
            / 3600.0
            / (self.config.cell_capacity_ah * self.config.cells_parallel as f64 * 2.0);
    }
}
