//! Extended types for electric vehicle simulation

#![allow(clippy::too_many_arguments)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{clarke_transform, inverse_park, park_transform, svpwm};
use super::types_core::{
    BatteryPack, CcCvCharger, CellChemistry, ChargingMode, DriveCycle, GridTariff, PackConfig,
    PiController, RangeEstimator, RegenerativeBraking, RoadSegment, SohEstimator,
};

/// Field-Oriented Controller for a PMSM (Permanent Magnet Synchronous Motor).
#[derive(Debug, Clone)]
pub struct FocController {
    /// d-axis PI controller
    pub id_ctrl: PiController,
    /// q-axis PI controller
    pub iq_ctrl: PiController,
    /// Motor electrical poles pairs
    pub pole_pairs: u32,
    /// Motor phase inductance (H)
    pub ls: f64,
    /// Motor phase resistance (Ω)
    pub rs: f64,
    /// Motor flux linkage (Wb)
    pub lambda: f64,
    /// Electrical angle (rad)
    pub theta_e: f64,
    /// Electrical angular velocity (rad/s)
    pub omega_e: f64,
    /// DC bus voltage (V)
    pub v_dc: f64,
}
impl FocController {
    /// Create a new FOC controller.
    pub fn new(pole_pairs: u32, ls: f64, rs: f64, lambda: f64, v_dc: f64) -> Self {
        let bandwidth = 1000.0;
        let kp = bandwidth * ls;
        let ki = bandwidth * rs;
        Self {
            id_ctrl: PiController::new(kp, ki, v_dc / f64::sqrt(3.0)),
            iq_ctrl: PiController::new(kp, ki, v_dc / f64::sqrt(3.0)),
            pole_pairs,
            ls,
            rs,
            lambda,
            theta_e: 0.0,
            omega_e: 0.0,
            v_dc,
        }
    }
    /// Update electrical angle from mechanical speed (rad/s).
    pub fn update_angle(&mut self, omega_mech: f64, dt_s: f64) {
        self.omega_e = omega_mech * self.pole_pairs as f64;
        self.theta_e = (self.theta_e + self.omega_e * dt_s) % (2.0 * std::f64::consts::PI);
    }
    /// Compute torque from q-axis current.
    pub fn torque_from_iq(&self, iq: f64) -> f64 {
        1.5 * self.pole_pairs as f64 * self.lambda * iq
    }
    /// Maximum torque per ampere (MTPA) — optimal id for given torque demand.
    ///
    /// For a surface-mount PMSM (no reluctance), id* = 0.
    pub fn mtpa_id(&self, _torque_demand: f64) -> f64 {
        0.0
    }
    /// Flux-weakening id to extend speed range above base speed.
    pub fn flux_weakening_id(&self, omega_mech: f64, iq: f64) -> f64 {
        let v_max = self.v_dc / f64::sqrt(3.0);
        let vd = self.rs * 0.0 - self.omega_e * self.ls * iq;
        let vq_max = (v_max.powi(2) - vd.powi(2)).sqrt().max(0.0);
        let vq_limit = vq_max;
        let vq_mtpa = self.rs * iq + self.omega_e * self.lambda;
        if vq_mtpa > vq_limit && omega_mech > 0.0 {
            -(self.lambda - vq_limit / (self.omega_e.max(1.0) * self.ls))
        } else {
            0.0
        }
    }
    /// Compute SVPWM duty cycles from current measurements.
    ///
    /// Returns (da, db, dc, torque_Nm).
    pub fn step(
        &mut self,
        ia: f64,
        ib: f64,
        ic: f64,
        id_ref: f64,
        iq_ref: f64,
        omega_mech: f64,
        dt_s: f64,
    ) -> (f64, f64, f64, f64) {
        self.update_angle(omega_mech, dt_s);
        let ab = clarke_transform(ia, ib, ic);
        let dq = park_transform(ab.alpha, ab.beta, self.theta_e);
        let ff_d = -self.omega_e * self.ls * dq.iq;
        let ff_q = self.omega_e * (self.ls * dq.id + self.lambda);
        let vd = self.id_ctrl.step(id_ref - dq.id, dt_s) + ff_d;
        let vq = self.iq_ctrl.step(iq_ref - dq.iq, dt_s) + ff_q;
        let (v_alpha, v_beta) = inverse_park(vd, vq, self.theta_e);
        let (da, db, dc) = svpwm(v_alpha, v_beta, self.v_dc);
        let torque = self.torque_from_iq(dq.iq);
        (da, db, dc, torque)
    }
}
/// CC-CV (Constant Current / Constant Voltage) charging model for a battery
/// pack.  Handles the two-phase charging profile and termination detection.
#[derive(Debug, Clone)]
pub struct ChargingModel {
    /// Constant current phase current (A).
    pub cc_current: f64,
    /// CV phase target voltage (V).
    pub cv_voltage: f64,
    /// Current charging mode.
    pub mode: ChargingMode,
    /// Minimum current to terminate charging (A).
    pub termination_current: f64,
    /// Total energy delivered (Wh).
    pub energy_delivered_wh: f64,
    /// Total time elapsed (s).
    pub time_s: f64,
    /// Whether charging is complete.
    pub complete: bool,
}
impl ChargingModel {
    /// Create a CC-CV charger with given parameters.
    ///
    /// * `cc_current`   — constant current phase amperage (A).
    /// * `cv_voltage`   — CV phase target voltage (V).
    /// * `termination`  — current at which charging terminates (A).
    pub fn new(cc_current: f64, cv_voltage: f64, termination_current: f64) -> Self {
        Self {
            cc_current,
            cv_voltage,
            mode: ChargingMode::ConstantCurrent,
            termination_current,
            energy_delivered_wh: 0.0,
            time_s: 0.0,
            complete: false,
        }
    }
    /// Standard 11 kW AC home charger model (32 A, 400 V).
    pub fn home_ac_11kw() -> Self {
        Self::new(32.0, 420.0, 3.0)
    }
    /// DC fast charger model (150 kW @ 400 V).
    pub fn dc_fast_150kw() -> Self {
        Self::new(375.0, 420.0, 10.0)
    }
    /// Step the charger for `dt` seconds given current pack voltage `v_pack` (V).
    ///
    /// Returns `(current_A, voltage_V)` applied this step.
    pub fn step(&mut self, v_pack: f64, dt: f64) -> (f64, f64) {
        if self.complete {
            return (0.0, v_pack);
        }
        self.time_s += dt;
        let (current, voltage) = match self.mode {
            ChargingMode::ConstantCurrent => {
                if v_pack >= self.cv_voltage {
                    self.mode = ChargingMode::ConstantVoltage;
                    (self.cc_current * 0.5, self.cv_voltage)
                } else {
                    (self.cc_current, v_pack)
                }
            }
            ChargingMode::ConstantVoltage => {
                let i = self.cc_current
                    * ((self.cv_voltage - v_pack) / (self.cv_voltage * 0.1 + 1.0)).clamp(0.0, 1.0);
                if i < self.termination_current {
                    self.complete = true;
                    (0.0, v_pack)
                } else {
                    (i, self.cv_voltage)
                }
            }
            _ => (0.0, v_pack),
        };
        self.energy_delivered_wh += current * voltage * dt / 3600.0;
        (current, voltage)
    }
    /// Reset the charger to start a new charging session.
    pub fn reset(&mut self) {
        self.mode = ChargingMode::ConstantCurrent;
        self.energy_delivered_wh = 0.0;
        self.time_s = 0.0;
        self.complete = false;
    }
}
/// Full electric vehicle system integrating all subsystems.
#[derive(Debug, Clone)]
pub struct ElectricVehicle {
    /// Battery pack
    pub pack: BatteryPack,
    /// FOC motor controller
    pub motor: FocController,
    /// Regen braking controller
    pub regen: RegenerativeBraking,
    /// CC/CV charger
    pub charger: CcCvCharger,
    /// V2G controller
    pub v2g: V2gController,
    /// Thermal management
    pub tms: ThermalManagement,
    /// SOH estimator
    pub soh_estimator: SohEstimator,
    /// Range estimator
    pub range_estimator: RangeEstimator,
    /// Vehicle speed (m/s)
    pub speed_ms: f64,
    /// Motor speed (rad/s)
    pub omega_motor: f64,
    /// Gear ratio
    pub gear_ratio: f64,
    /// Odometer (m)
    pub odometer_m: f64,
    /// Total energy consumed (kWh)
    pub energy_consumed_kwh: f64,
}
impl ElectricVehicle {
    /// Create a representative 75 kWh EV.
    pub fn new_75kwh() -> Self {
        let config = PackConfig {
            cells_series: 96,
            cells_parallel: 4,
            chemistry: CellChemistry::NMC,
            cell_capacity_ah: 60.0,
            thermal_resistance: 0.05,
            coolant_temp: 20.0,
        };
        let pack = BatteryPack::new(config.clone());
        let motor = FocController::new(4, 0.3e-3, 0.02, 0.1, 400.0);
        let regen = RegenerativeBraking::new(200.0, 80_000.0);
        let charger = CcCvCharger::new(200.0, config.nominal_voltage() * 1.05);
        let v2g = V2gController::new(11_000.0, 22_000.0, 0.25);
        let tms = ThermalManagement::new();
        let soh_estimator =
            SohEstimator::new(config.cell_capacity_ah * config.cells_parallel as f64);
        let range_estimator = RangeEstimator::new(2000.0, 2.4, 0.28, 75.0);
        Self {
            pack,
            motor,
            regen,
            charger,
            v2g,
            tms,
            soh_estimator,
            range_estimator,
            speed_ms: 0.0,
            omega_motor: 0.0,
            gear_ratio: 9.0,
            odometer_m: 0.0,
            energy_consumed_kwh: 0.0,
        }
    }
    /// Simulate one time step given driver torque demand (Nm) and dt (s).
    pub fn step(&mut self, torque_demand_nm: f64, dt_s: f64, t_ambient: f64) {
        let wheel_radius = 0.33;
        self.omega_motor = self.speed_ms / wheel_radius * self.gear_ratio;
        let iq_ref = (torque_demand_nm / (1.5 * self.motor.pole_pairs as f64 * self.motor.lambda))
            .clamp(-400.0, 400.0);
        let id_ref = self.motor.flux_weakening_id(self.omega_motor, iq_ref);
        let (_da, _db, _dc, torque_out) = self.motor.step(
            iq_ref * 0.99,
            iq_ref * (-0.5),
            iq_ref * (-0.5),
            id_ref,
            iq_ref,
            self.omega_motor,
            dt_s,
        );
        let p_mech = torque_out * self.omega_motor;
        let rpm = self.omega_motor * 60.0 / std::f64::consts::TAU;
        let eta = self
            .regen
            .efficiency_map
            .efficiency(rpm.abs(), torque_out.abs());
        let p_elec = if p_mech > 0.0 {
            p_mech / eta.max(0.1)
        } else {
            p_mech * eta
        };
        let pack_current = p_elec / self.pack.pack_voltage().max(1.0);
        self.pack.update(pack_current, dt_s);
        self.tms.update(&self.pack, t_ambient, dt_s);
        self.odometer_m += self.speed_ms * dt_s;
        if p_elec > 0.0 {
            self.energy_consumed_kwh += p_elec * dt_s / 3_600_000.0;
        }
    }
    /// Estimate current range (km).
    pub fn range_km(&self) -> f64 {
        let cycle = DriveCycle::wltp_synthetic();
        self.range_estimator
            .remaining_range_km(self.pack.pack_soc, &cycle)
    }
}
/// Simple charging session model.
#[derive(Debug, Clone)]
pub struct ChargeSession {
    /// Charging power (kW).
    pub power_kw: f64,
    /// Duration of the charging session (h).
    pub duration_h: f64,
    /// Charger-to-battery efficiency in `[0, 1]`.
    pub efficiency: f64,
}
impl ChargeSession {
    /// Energy delivered to the battery (kWh).
    ///
    /// `E = P · t · η`
    pub fn energy_delivered_kwh(&self) -> f64 {
        self.power_kw * self.duration_h * self.efficiency.clamp(0.0, 1.0)
    }
    /// Estimated cost of the charging session.
    ///
    /// `cost = P · t · price_per_kwh` (cost is on grid energy, before efficiency)
    pub fn cost_estimate(&self, price_per_kwh: f64) -> f64 {
        self.power_kw * self.duration_h * price_per_kwh
    }
}
/// A simple electric motor model with an efficiency map.
#[derive(Debug, Clone)]
pub struct ElectricMotor {
    /// Peak torque (N·m).
    pub max_torque_nm: f64,
    /// Peak power (kW).
    pub max_power_kw: f64,
    /// Efficiency map: `Vec<(speed_rpm, torque_nm, efficiency)>`.
    pub efficiency_map: Vec<(f64, f64, f64)>,
}
impl ElectricMotor {
    /// Create a simple 200 kW motor.
    pub fn typical_200kw() -> Self {
        Self {
            max_torque_nm: 450.0,
            max_power_kw: 200.0,
            efficiency_map: vec![
                (0.0, 0.0, 0.85),
                (1000.0, 200.0, 0.90),
                (3000.0, 300.0, 0.95),
                (6000.0, 200.0, 0.93),
                (9000.0, 100.0, 0.88),
            ],
        }
    }
}
/// Vehicle-to-Grid (V2G) controller.
#[derive(Debug, Clone)]
pub struct V2gController {
    /// Max power export to grid (W)
    pub max_export_w: f64,
    /// Max power import from grid (W)
    pub max_import_w: f64,
    /// Target SOC to maintain for vehicle usage
    pub reserve_soc: f64,
    /// Departure SOC target
    pub departure_soc: f64,
    /// V2G enabled flag
    pub v2g_enabled: bool,
    /// Grid tariff
    pub tariff: GridTariff,
    /// Total energy exported (kWh)
    pub energy_exported_kwh: f64,
    /// Total energy imported (kWh)
    pub energy_imported_kwh: f64,
    /// Revenue earned (€)
    pub revenue: f64,
}
impl V2gController {
    /// Create a V2G controller.
    pub fn new(max_export_w: f64, max_import_w: f64, reserve_soc: f64) -> Self {
        Self {
            max_export_w,
            max_import_w,
            reserve_soc,
            departure_soc: 0.8,
            v2g_enabled: true,
            tariff: GridTariff::two_rate(0.30, 0.10),
            energy_exported_kwh: 0.0,
            energy_imported_kwh: 0.0,
            revenue: 0.0,
        }
    }
    /// Compute power setpoint (W). Positive = charge, negative = discharge (export).
    pub fn compute_setpoint(
        &self,
        current_soc: f64,
        hour_of_day: f64,
        hours_to_departure: f64,
    ) -> f64 {
        if !self.v2g_enabled {
            return if current_soc < self.departure_soc {
                self.max_import_w
            } else {
                0.0
            };
        }
        let price = self.tariff.price_at_hour(hour_of_day);
        let mean_price = (self.tariff.price[0] + self.tariff.price[1]) / 2.0;
        if current_soc > self.reserve_soc + 0.1 && price > mean_price * 1.2 {
            -self.max_export_w * (current_soc - self.reserve_soc).min(0.3) / 0.3
        } else if current_soc < self.departure_soc && price < mean_price * 0.8 {
            let time_factor = if hours_to_departure > 0.0 {
                (self.departure_soc - current_soc) / hours_to_departure.max(0.1) * 100.0
            } else {
                self.max_import_w
            };
            time_factor.min(self.max_import_w)
        } else {
            0.0
        }
    }
    /// Update energy accounting.
    pub fn update_accounting(&mut self, power_w: f64, dt_s: f64, hour_of_day: f64) {
        let energy_kwh = power_w.abs() * dt_s / 3_600_000.0;
        let price = self.tariff.price_at_hour(hour_of_day);
        if power_w < 0.0 {
            self.energy_exported_kwh += energy_kwh;
            self.revenue += energy_kwh * price;
        } else if power_w > 0.0 {
            self.energy_imported_kwh += energy_kwh;
            self.revenue -= energy_kwh * price;
        }
    }
}
/// Simple energy-optimal route planner.
#[derive(Debug, Clone)]
pub struct EnergyRouter {
    /// Segments (pre-connected graph in order)
    pub segments: Vec<RoadSegment>,
    /// Vehicle mass (kg)
    pub mass_kg: f64,
    /// Drag area (cd * A)
    pub cd_area: f64,
    /// Drivetrain efficiency
    pub eta_drive: f64,
    /// Regen efficiency
    pub eta_regen: f64,
}
impl EnergyRouter {
    /// Air density.
    const RHO: f64 = 1.225;
    /// Gravity.
    const G: f64 = 9.81;
    /// Create an energy router.
    pub fn new(mass_kg: f64, cd_area: f64) -> Self {
        Self {
            segments: Vec::new(),
            mass_kg,
            cd_area,
            eta_drive: 0.92,
            eta_regen: 0.75,
        }
    }
    /// Compute energy for a single segment (Wh).
    pub fn segment_energy_wh(&self, seg: &RoadSegment) -> f64 {
        let v = seg.avg_speed_ms;
        let f_drag = 0.5 * Self::RHO * self.cd_area * v.powi(2);
        let f_roll = self.mass_kg * Self::G * 0.01 * seg.surface_factor;
        let f_grade = self.mass_kg * Self::G * seg.gradient_rad.sin();
        let f_total = f_drag + f_roll + f_grade;
        let p_mech = f_total * v;
        let p_elec = if p_mech > 0.0 {
            p_mech / self.eta_drive
        } else {
            p_mech * self.eta_regen
        };
        let t = seg.length_m / v.max(0.1);
        p_elec * t / 3600.0
    }
    /// Total route energy (Wh).
    pub fn total_energy_wh(&self) -> f64 {
        self.segments
            .iter()
            .map(|s| self.segment_energy_wh(s))
            .sum()
    }
    /// Estimate if route is feasible given pack SOC and capacity.
    pub fn is_feasible(&self, soc: f64, pack_kwh: f64) -> bool {
        let needed_kwh = self.total_energy_wh() / 1000.0;
        needed_kwh <= soc * pack_kwh
    }
    /// Greedy alternative route selector (returns segment subset with minimum energy).
    pub fn minimum_energy_route<'a>(
        &'a self,
        alt_segments: &'a [Vec<RoadSegment>],
    ) -> &'a [RoadSegment] {
        let mut best_idx = 0;
        let mut best_energy = f64::MAX;
        for (i, path) in alt_segments.iter().enumerate() {
            let e: f64 = path.iter().map(|s| self.segment_energy_wh(s)).sum();
            if e < best_energy {
                best_energy = e;
                best_idx = i;
            }
        }
        &alt_segments[best_idx]
    }
}
/// Simple electric motor efficiency map.
///
/// Models peak torque (constant up to `base_speed`, then field-weakening),
/// maximum power, and a uniform efficiency look-up.
#[derive(Debug, Clone)]
pub struct MotorEfficiencyMap {
    /// Maximum torque at stall / low speed (N·m).
    pub max_torque: f64,
    /// Peak power (W).
    pub max_power: f64,
    /// Corner speed where field-weakening begins (rad/s).
    pub base_speed: f64,
    /// Efficiency map stored row-by-row (torque × speed grid).
    pub efficiency_map: Vec<Vec<f64>>,
}
impl MotorEfficiencyMap {
    /// Create a motor map with uniform efficiency of 0.90.
    ///
    /// * `max_torque`  – stall torque (N·m)
    /// * `max_power`   – peak power (W)
    /// * `base_speed`  – corner speed (rad/s)
    pub fn new(max_torque: f64, max_power: f64, base_speed: f64) -> Self {
        let efficiency_map = vec![vec![0.90; 5]; 5];
        Self {
            max_torque,
            max_power,
            base_speed,
            efficiency_map,
        }
    }
    /// Look up efficiency at a given operating point.
    ///
    /// Returns the nearest cell in the efficiency map.  Falls back to 0.90
    /// for out-of-range or empty maps.
    pub fn efficiency(&self, torque: f64, speed: f64) -> f64 {
        if self.efficiency_map.is_empty() || self.efficiency_map[0].is_empty() {
            return 0.90;
        }
        let rows = self.efficiency_map.len();
        let cols = self.efficiency_map[0].len();
        let row = ((torque / self.max_torque.max(1e-9)) * (rows - 1) as f64)
            .round()
            .clamp(0.0, (rows - 1) as f64) as usize;
        let col = ((speed / self.max_power.max(1e-9)) * (cols - 1) as f64)
            .round()
            .clamp(0.0, (cols - 1) as f64) as usize;
        self.efficiency_map[row][col]
    }
    /// Maximum torque available at a given speed (N·m).
    ///
    /// Constant torque up to `base_speed`; above that, torque decreases as
    /// `P_max / speed`.
    pub fn max_torque_at_speed(&self, speed: f64) -> f64 {
        if speed < 1e-9 {
            return self.max_torque;
        }
        if speed <= self.base_speed {
            self.max_torque
        } else {
            (self.max_power / speed).min(self.max_torque)
        }
    }
    /// Mechanical power at the given operating point (W).
    pub fn power_at_operating_point(&self, torque: f64, speed: f64) -> f64 {
        torque * speed
    }
}
/// Battery thermal management system.
#[derive(Debug, Clone)]
pub struct ThermalManagement {
    /// Coolant flow rate (L/min)
    pub flow_rate_lpm: f64,
    /// Coolant specific heat (J/(kg·K))
    pub cp_coolant: f64,
    /// Coolant density (kg/L)
    pub rho_coolant: f64,
    /// Chiller COP
    pub chiller_cop: f64,
    /// Heater power (W)
    pub heater_power_w: f64,
    /// Target min temperature (°C)
    pub t_min: f64,
    /// Target max temperature (°C)
    pub t_max: f64,
    /// Current coolant temperature (°C)
    pub t_coolant: f64,
    /// Compressor power (W)
    pub compressor_power_w: f64,
}
impl ThermalManagement {
    /// Create a standard TMS.
    pub fn new() -> Self {
        Self {
            flow_rate_lpm: 10.0,
            cp_coolant: 3900.0,
            rho_coolant: 1.07,
            chiller_cop: 3.0,
            heater_power_w: 4000.0,
            t_min: 15.0,
            t_max: 35.0,
            t_coolant: 20.0,
            compressor_power_w: 0.0,
        }
    }
    /// Update coolant temperature given battery pack and ambient conditions.
    ///
    /// Returns auxiliary power consumed (W).
    pub fn update(&mut self, pack: &BatteryPack, t_ambient: f64, dt_s: f64) -> f64 {
        let t_pack_avg =
            pack.cells.iter().map(|c| c.temperature).sum::<f64>() / pack.cells.len() as f64;
        let q_gen: f64 = pack.cells.iter().map(|c| c.heat_generation_rate(1.0)).sum();
        let mut aux_power = 0.0;
        if t_pack_avg > self.t_max {
            let q_cool = (t_pack_avg - self.t_max) * 5000.0;
            self.compressor_power_w = q_cool / self.chiller_cop;
            aux_power += self.compressor_power_w;
            self.t_coolant = t_ambient - 5.0;
        } else if t_pack_avg < self.t_min {
            aux_power += self.heater_power_w;
            self.t_coolant = t_pack_avg + 5.0;
            self.compressor_power_w = 0.0;
        } else {
            self.compressor_power_w = 0.0;
        }
        let m_dot = self.flow_rate_lpm * self.rho_coolant / 60.0;
        let _dt_coolant = (q_gen - m_dot * self.cp_coolant * (self.t_coolant - t_ambient)) * dt_s
            / (m_dot * self.cp_coolant + 1.0);
        aux_power
    }
    /// Check if pre-conditioning is required before charging.
    pub fn needs_preconditioning(&self, pack: &BatteryPack) -> bool {
        let t_pack_avg =
            pack.cells.iter().map(|c| c.temperature).sum::<f64>() / pack.cells.len() as f64;
        t_pack_avg < self.t_min || t_pack_avg > self.t_max
    }
}
