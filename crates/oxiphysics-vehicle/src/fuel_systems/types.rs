//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use std::f64::consts::PI;

/// SOC-based rule strategy for a parallel hybrid electric vehicle (P-HEV).
///
/// Implements a charge-sustaining strategy with power-split between the
/// internal combustion engine (ICE) and the electric motor (EM).
#[derive(Debug, Clone)]
pub struct HybridEnergyManagement {
    /// Target state of charge (SOC) set-point \[0.0, 1.0\].
    pub soc_target: f64,
    /// SOC lower threshold: below this, force charge mode.
    pub soc_low: f64,
    /// SOC upper threshold: above this, use electric-only if possible.
    pub soc_high: f64,
    /// Maximum ICE power (W).
    pub ice_max_power_w: f64,
    /// Maximum EM power (W, motoring and generating).
    pub em_max_power_w: f64,
    /// Battery state of charge \[0.0, 1.0\].
    pub soc: f64,
    /// Battery capacity (kWh).
    pub battery_capacity_kwh: f64,
    /// ICE fuel consumption slope (g/kWh) — for fuel use estimation.
    pub ice_bsfc_g_per_kwh: f64,
}
impl HybridEnergyManagement {
    /// Construct a typical P-HEV energy manager.
    pub fn new(
        ice_max_power_w: f64,
        em_max_power_w: f64,
        battery_capacity_kwh: f64,
        initial_soc: f64,
    ) -> Self {
        Self {
            soc_target: 0.6_f64,
            soc_low: 0.3_f64,
            soc_high: 0.8_f64,
            ice_max_power_w,
            em_max_power_w,
            soc: initial_soc.clamp(0.0_f64, 1.0_f64),
            battery_capacity_kwh,
            ice_bsfc_g_per_kwh: 240.0_f64,
        }
    }
    /// Compute the power split between ICE and EM for a given driver demand.
    ///
    /// Returns `(ice_power_w, em_power_w)`.  Positive EM power = motoring;
    /// negative EM power = generating (charging).
    pub fn power_split(&self, demanded_power_w: f64) -> (f64, f64) {
        if self.soc < self.soc_low {
            let ice_power = (demanded_power_w + 0.1_f64 * self.ice_max_power_w)
                .clamp(0.0_f64, self.ice_max_power_w);
            let em_power = demanded_power_w - ice_power;
            return (ice_power, em_power.max(-self.em_max_power_w));
        }
        if self.soc > self.soc_high {
            let em_power = demanded_power_w.min(self.em_max_power_w);
            let ice_power = (demanded_power_w - em_power).max(0.0_f64);
            return (ice_power, em_power);
        }
        let ice_power = demanded_power_w.clamp(0.0_f64, self.ice_max_power_w);
        let em_power = demanded_power_w - ice_power;
        let em_clamped = em_power.clamp(-self.em_max_power_w, self.em_max_power_w);
        (ice_power, em_clamped)
    }
    /// Update SOC given EM power `em_power_w` and time step `dt` (s).
    ///
    /// Negative `em_power_w` means the motor is generating (charging).
    pub fn update_soc(&mut self, em_power_w: f64, dt: f64) {
        let energy_kwh = -em_power_w * dt / 3.6e6_f64;
        let d_soc = energy_kwh / self.battery_capacity_kwh.max(1.0e-9_f64);
        self.soc = (self.soc + d_soc).clamp(0.0_f64, 1.0_f64);
    }
    /// Estimate instantaneous fuel consumption (g/s) from ICE power.
    pub fn fuel_consumption_g_s(&self, ice_power_w: f64) -> f64 {
        let ice_kw = ice_power_w / 1000.0_f64;
        ice_kw * self.ice_bsfc_g_per_kwh / 3600.0_f64
    }
    /// True if electric-only operation is feasible.
    pub fn ev_mode_available(&self) -> bool {
        self.soc > self.soc_low
    }
}
/// Stoichiometric combustion parameters for a hydrocarbon fuel.
#[derive(Debug, Clone)]
pub struct FuelProperties {
    /// Stoichiometric air-fuel ratio (kg air / kg fuel).
    pub stoich_afr: f64,
    /// Lower heating value (MJ/kg).
    pub lhv_mj_per_kg: f64,
    /// Carbon-to-hydrogen mass ratio (C/H).
    pub ch_ratio: f64,
}
impl FuelProperties {
    /// Standard gasoline (iso-octane basis).
    pub fn gasoline() -> Self {
        Self {
            stoich_afr: 14.7,
            lhv_mj_per_kg: 44.0,
            ch_ratio: 5.25,
        }
    }
    /// Standard diesel (cetane basis).
    pub fn diesel() -> Self {
        Self {
            stoich_afr: 14.5,
            lhv_mj_per_kg: 42.5,
            ch_ratio: 5.83,
        }
    }
    /// Ethanol (E100).
    pub fn ethanol() -> Self {
        Self {
            stoich_afr: 9.0,
            lhv_mj_per_kg: 26.8,
            ch_ratio: 3.0,
        }
    }
}
/// Kinetic Energy Recovery System (KERS) — flywheel based.
#[derive(Debug, Clone)]
pub struct KersSystem {
    /// Flywheel moment of inertia (kg·m²).
    pub inertia_kg_m2: f64,
    /// Current flywheel angular velocity (rad/s).
    pub omega_rad_s: f64,
    /// Maximum flywheel speed (rad/s).
    pub max_omega_rad_s: f64,
    /// Gear ratio flywheel ↔ drivetrain.
    pub gear_ratio: f64,
    /// CVT efficiency.
    pub cvt_efficiency: f64,
    /// Stored energy (J).
    pub stored_energy_j: f64,
}
impl KersSystem {
    /// Construct a typical flywheel KERS.
    pub fn new(inertia_kg_m2: f64, max_omega_rad_s: f64) -> Self {
        Self {
            inertia_kg_m2,
            omega_rad_s: 0.0,
            max_omega_rad_s,
            gear_ratio: 5.0,
            cvt_efficiency: 0.9,
            stored_energy_j: 0.0,
        }
    }
    /// Store kinetic energy from braking torque `torque_nm` over `dt` (s).
    pub fn absorb_energy(&mut self, torque_nm: f64, dt: f64) {
        let power_w = torque_nm * self.omega_rad_s;
        let energy_in = power_w * dt * self.cvt_efficiency;
        let delta_omega =
            (torque_nm * dt / self.inertia_kg_m2).min(self.max_omega_rad_s - self.omega_rad_s);
        self.omega_rad_s = (self.omega_rad_s + delta_omega.max(0.0)).min(self.max_omega_rad_s);
        self.stored_energy_j = 0.5 * self.inertia_kg_m2 * self.omega_rad_s * self.omega_rad_s;
        let _ = energy_in;
    }
    /// Release energy from flywheel as a torque to the drivetrain for `dt` (s).
    ///
    /// Returns the torque actually delivered (N·m).
    pub fn release_energy(&mut self, requested_torque_nm: f64, dt: f64) -> f64 {
        if self.omega_rad_s < 1.0 {
            return 0.0;
        }
        let max_torque = self.inertia_kg_m2 * self.omega_rad_s / dt;
        let torque = requested_torque_nm.min(max_torque);
        let delta_omega = torque * dt / self.inertia_kg_m2;
        self.omega_rad_s = (self.omega_rad_s - delta_omega).max(0.0);
        self.stored_energy_j = 0.5 * self.inertia_kg_m2 * self.omega_rad_s * self.omega_rad_s;
        torque * self.cvt_efficiency
    }
    /// Fraction of maximum stored energy currently in flywheel.
    pub fn state_of_charge(&self) -> f64 {
        let max_e = 0.5 * self.inertia_kg_m2 * self.max_omega_rad_s * self.max_omega_rad_s;
        if max_e < 1e-15 {
            return 0.0;
        }
        (self.stored_energy_j / max_e).clamp(0.0, 1.0)
    }
}
/// Fuel injector: direct injection timing, pulse width, flow rate and atomization.
///
/// Models a solenoid injector for port or direct injection engines.
#[derive(Debug, Clone)]
pub struct FuelInjector {
    /// Maximum volumetric flow rate (mL/min) at full open.
    pub max_flow_ml_per_min: f64,
    /// Injector opening delay (ms).
    pub opening_delay_ms: f64,
    /// Injector closing delay (ms).
    pub closing_delay_ms: f64,
    /// Fuel rail pressure (bar).
    pub rail_pressure_bar: f64,
    /// Current injection pulse width (ms).
    pub pulse_width_ms: f64,
    /// Injection start angle (°CA before TDC).
    pub start_angle_btdc: f64,
    /// Injector nozzle diameter (mm).
    pub nozzle_diameter_mm: f64,
    /// Number of nozzle holes.
    pub n_holes: u32,
    /// Total fuel injected in the current cycle (mg).
    pub injected_mass_mg: f64,
}
impl FuelInjector {
    /// Construct a direct-injection gasoline injector with default parameters.
    pub fn new_gdi() -> Self {
        Self {
            max_flow_ml_per_min: 200.0,
            opening_delay_ms: 0.3,
            closing_delay_ms: 0.2,
            rail_pressure_bar: 200.0,
            pulse_width_ms: 2.5,
            start_angle_btdc: 300.0,
            nozzle_diameter_mm: 0.2,
            n_holes: 6,
            injected_mass_mg: 0.0,
        }
    }
    /// Construct a port fuel injection (PFI) injector.
    pub fn new_pfi() -> Self {
        Self {
            max_flow_ml_per_min: 300.0,
            opening_delay_ms: 0.5,
            closing_delay_ms: 0.3,
            rail_pressure_bar: 4.0,
            pulse_width_ms: 3.5,
            start_angle_btdc: 180.0,
            nozzle_diameter_mm: 0.4,
            n_holes: 2,
            injected_mass_mg: 0.0,
        }
    }
    /// Effective flow rate (mL/min) at the given rail pressure.
    ///
    /// Flow scales as √(rail_pressure / reference_pressure).
    pub fn effective_flow_rate_ml_per_min(&self) -> f64 {
        let ref_pressure = if self.rail_pressure_bar > 100.0 {
            200.0
        } else {
            4.0
        };
        self.max_flow_ml_per_min * (self.rail_pressure_bar / ref_pressure).sqrt()
    }
    /// Compute the fuel mass injected per pulse (mg) for gasoline density.
    ///
    /// Assumes fuel density 0.745 g/mL.
    pub fn pulse_mass_mg(&self) -> f64 {
        let density_g_per_ml = 0.745;
        let flow_ml_per_ms = self.effective_flow_rate_ml_per_min() / 60_000.0;
        let effective_pw =
            (self.pulse_width_ms - self.opening_delay_ms - self.closing_delay_ms).max(0.0);
        flow_ml_per_ms * effective_pw * density_g_per_ml * 1000.0
    }
    /// Sauter Mean Diameter (SMD) of the fuel spray droplets (µm).
    ///
    /// Uses the Hiroyasu-Arai correlation:
    /// SMD = C · d^0.25 · ΔP^{-0.25} · ρ_f^{0.25}
    /// where C is an empirical constant.
    pub fn sauter_mean_diameter_um(&self) -> f64 {
        let c = 2.33e-3;
        let d_m = self.nozzle_diameter_mm * 1e-3;
        let delta_p = self.rail_pressure_bar * 1e5;
        let rho_f: f64 = 745.0;
        let smd_m = c * d_m.powf(0.25) * delta_p.powf(-0.25) * rho_f.powf(0.25);
        smd_m * 1e6
    }
    /// Spray cone half-angle (degrees) — increases with rail pressure.
    pub fn spray_cone_angle_deg(&self) -> f64 {
        15.0 + 10.0 * (self.rail_pressure_bar / 200.0).sqrt()
    }
    /// Simulate one injection event: update `injected_mass_mg`.
    pub fn fire(&mut self) {
        self.injected_mass_mg = self.pulse_mass_mg();
    }
    /// Reset the injected mass counter.
    pub fn reset_cycle(&mut self) {
        self.injected_mass_mg = 0.0;
    }
}
/// Regenerative braking controller.
#[derive(Debug, Clone)]
pub struct RegenerativeBraking {
    /// Maximum regenerative braking torque (N·m).
    pub max_regen_torque_nm: f64,
    /// Blending threshold — fraction of total braking that is regenerative.
    pub regen_fraction: f64,
    /// Regeneration efficiency (motor/generator efficiency).
    pub regen_efficiency: f64,
    /// Minimum wheel speed for regeneration (rad/s).
    pub min_wheel_speed_rad_s: f64,
    /// Current regenerated power (W).
    pub regen_power_w: f64,
}
impl RegenerativeBraking {
    /// Construct a regenerative braking controller.
    pub fn new(max_regen_torque_nm: f64, regen_fraction: f64, regen_efficiency: f64) -> Self {
        Self {
            max_regen_torque_nm,
            regen_fraction,
            regen_efficiency,
            min_wheel_speed_rad_s: 2.0,
            regen_power_w: 0.0,
        }
    }
    /// Compute regenerative torque for a requested total braking torque.
    ///
    /// Returns the regenerative torque component (N·m).
    pub fn regen_torque(&self, total_brake_torque_nm: f64, wheel_speed_rad_s: f64) -> f64 {
        if wheel_speed_rad_s < self.min_wheel_speed_rad_s {
            return 0.0;
        }
        (total_brake_torque_nm * self.regen_fraction).min(self.max_regen_torque_nm)
    }
    /// Compute recovered power (W) from regen torque at given wheel speed.
    pub fn recovered_power_w(&self, total_brake_torque_nm: f64, wheel_speed_rad_s: f64) -> f64 {
        let tau = self.regen_torque(total_brake_torque_nm, wheel_speed_rad_s);
        tau * wheel_speed_rad_s * self.regen_efficiency
    }
    /// Update internal regen power state.
    pub fn update(&mut self, total_brake_torque_nm: f64, wheel_speed_rad_s: f64) {
        self.regen_power_w = self.recovered_power_w(total_brake_torque_nm, wheel_speed_rad_s);
    }
}
/// Three-way catalytic converter efficiency model.
#[derive(Debug, Clone)]
pub struct Catalyst {
    /// Catalyst bed temperature (K).
    pub temperature_k: f64,
    /// Light-off temperature (K) — conversion drops to ~0 below this.
    pub light_off_k: f64,
    /// Peak conversion efficiency (at stoich, hot catalyst).
    pub peak_efficiency: f64,
    /// Lambda window half-width for high efficiency.
    pub lambda_window: f64,
    /// Thermal time constant (s).
    pub thermal_tau_s: f64,
    /// Catalyst thermal mass × heat capacity (J/K).
    pub thermal_capacitance: f64,
}
impl Catalyst {
    /// Construct a standard TWC.
    pub fn new_twc() -> Self {
        Self {
            temperature_k: 300.0,
            light_off_k: 573.0,
            peak_efficiency: 0.98,
            lambda_window: 0.02,
            thermal_tau_s: 30.0,
            thermal_capacitance: 500.0,
        }
    }
    /// Conversion efficiency for a given lambda and current temperature.
    pub fn conversion_efficiency(&self, lambda: f64) -> f64 {
        if self.temperature_k < self.light_off_k {
            return 0.0;
        }
        let temp_factor = ((self.temperature_k - self.light_off_k) / 100.0).min(1.0);
        let lambda_factor = 1.0 - ((lambda - 1.0).abs() / self.lambda_window).min(1.0);
        self.peak_efficiency * temp_factor * lambda_factor
    }
    /// Update catalyst temperature given exhaust temperature `t_exhaust_k` and `dt` (s).
    pub fn update_temperature(&mut self, t_exhaust_k: f64, dt: f64) {
        let tau = self.thermal_tau_s;
        self.temperature_k += (t_exhaust_k - self.temperature_k) / tau * dt;
    }
}
/// Hybrid traction battery with SOC tracking, thermal model and aging.
#[derive(Debug, Clone)]
pub struct HybridBattery {
    /// Nominal capacity (kWh).
    pub capacity_kwh: f64,
    /// State of charge \[0.0, 1.0\].
    pub soc: f64,
    /// Internal resistance (Ω).
    pub internal_resistance_ohm: f64,
    /// Nominal cell voltage (V).
    pub nominal_voltage_v: f64,
    /// Battery temperature (°C).
    pub temperature_c: f64,
    /// Thermal resistance to ambient (°C/W).
    pub thermal_resistance: f64,
    /// Thermal capacitance (J/°C).
    pub thermal_capacitance: f64,
    /// State of health \[0.0, 1.0\] — decreases with aging.
    pub state_of_health: f64,
    /// Total energy throughput (kWh) — used for aging estimate.
    pub energy_throughput_kwh: f64,
    /// Cycle life at 100% DoD before 20% capacity loss.
    pub cycle_life: f64,
}
impl HybridBattery {
    /// Construct a typical NMC hybrid battery pack.
    ///
    /// # Arguments
    /// * `capacity_kwh` — nominal capacity.
    /// * `initial_soc` — starting state of charge \[0.0, 1.0\].
    pub fn new_nmc(capacity_kwh: f64, initial_soc: f64) -> Self {
        Self {
            capacity_kwh,
            soc: initial_soc.clamp(0.0, 1.0),
            internal_resistance_ohm: 0.1,
            nominal_voltage_v: 400.0,
            temperature_c: 25.0,
            thermal_resistance: 0.05,
            thermal_capacitance: 10_000.0,
            state_of_health: 1.0,
            energy_throughput_kwh: 0.0,
            cycle_life: 2000.0,
        }
    }
    /// Open-circuit voltage (V) as a function of SOC.
    ///
    /// Uses a piecewise linear approximation of an NMC OCV curve.
    pub fn open_circuit_voltage(&self) -> f64 {
        let table = [
            (0.0_f64, 3.0_f64),
            (0.1, 3.4),
            (0.2, 3.55),
            (0.4, 3.72),
            (0.6, 3.85),
            (0.8, 4.0),
            (0.9, 4.1),
            (1.0, 4.2),
        ];
        let soc = self.soc.clamp(0.0, 1.0);
        for i in 0..table.len() - 1 {
            let (s0, v0) = table[i];
            let (s1, v1) = table[i + 1];
            if soc <= s1 {
                let t = (soc - s0) / (s1 - s0);
                return v0 + t * (v1 - v0);
            }
        }
        table.last().expect("collection should not be empty").1
    }
    /// Available energy (kWh) given current SOC and SoH.
    pub fn available_energy_kwh(&self) -> f64 {
        self.capacity_kwh * self.soc * self.state_of_health
    }
    /// Charge the battery at power `power_kw` for `dt` seconds.
    ///
    /// Returns actual power accepted (limited by SOC ceiling).
    pub fn charge(&mut self, power_kw: f64, dt: f64) -> f64 {
        let energy_kwh = power_kw * dt / 3600.0;
        let head_room = self.capacity_kwh * (1.0 - self.soc) * self.state_of_health;
        let accepted_kwh = energy_kwh.min(head_room).max(0.0);
        self.soc = (self.soc + accepted_kwh / (self.capacity_kwh * self.state_of_health + 1e-15))
            .clamp(0.0, 1.0);
        self.energy_throughput_kwh += accepted_kwh;
        self._update_aging();
        accepted_kwh * 3600.0 / dt
    }
    /// Discharge the battery at power `power_kw` for `dt` seconds.
    ///
    /// Returns actual power delivered (limited by available energy).
    pub fn discharge(&mut self, power_kw: f64, dt: f64) -> f64 {
        let energy_kwh = power_kw * dt / 3600.0;
        let available = self.available_energy_kwh();
        let delivered_kwh = energy_kwh.min(available).max(0.0);
        self.soc = (self.soc - delivered_kwh / (self.capacity_kwh * self.state_of_health + 1e-15))
            .clamp(0.0, 1.0);
        self.energy_throughput_kwh += delivered_kwh;
        self._update_aging();
        delivered_kwh * 3600.0 / dt
    }
    /// Update thermal model for one time step `dt` (s) at heat generation `q_w` (W).
    pub fn update_thermal(&mut self, q_w: f64, ambient_c: f64, dt: f64) {
        let q_dissipated = (self.temperature_c - ambient_c) / self.thermal_resistance;
        let delta_temp = (q_w - q_dissipated) / self.thermal_capacitance * dt;
        self.temperature_c += delta_temp;
    }
    /// Compute ohmic heat generation for a given current (A).
    pub fn ohmic_heat_w(&self, current_a: f64) -> f64 {
        current_a * current_a * self.internal_resistance_ohm
    }
    /// Update state-of-health based on energy throughput.
    fn _update_aging(&mut self) {
        let full_cycle_energy = self.capacity_kwh;
        let n_equiv_cycles = self.energy_throughput_kwh / (full_cycle_energy + 1e-15);
        let degradation = 0.2 * (n_equiv_cycles / self.cycle_life).min(1.0);
        self.state_of_health = (1.0 - degradation).max(0.0);
    }
    /// Terminal voltage under load (V).
    ///
    /// V_terminal = V_oc - I × R_int.  Current derived from power ÷ V_oc.
    pub fn terminal_voltage(&self, power_kw: f64) -> f64 {
        let v_oc = self.open_circuit_voltage();
        if v_oc < 1e-3 {
            return 0.0;
        }
        let current_a = power_kw * 1000.0 / v_oc;
        v_oc - current_a * self.internal_resistance_ohm
    }
}
/// Extended emissions model for NOx, CO, HC and particulate matter.
///
/// Extends the per-cycle values in [`CombustionModel`] to mass-flow rates
/// and Euro-standard compliance checking.
#[derive(Debug, Clone)]
pub struct EmissionsModel {
    /// Engine displacement volume (L).
    pub displacement_l: f64,
    /// Number of cylinders.
    pub n_cylinders: u32,
    /// Engine speed (rpm).
    pub engine_rpm: f64,
    /// Fuel mass per cycle per cylinder (mg).
    pub fuel_mass_mg: f64,
    /// Lambda (equivalence ratio λ = AFR / stoich_AFR).
    pub lambda: f64,
    /// Peak in-cylinder temperature (K).
    pub peak_temp_k: f64,
    /// EGR fraction (0–1) — reduces NOx.
    pub egr_fraction: f64,
    /// Catalyst conversion efficiency (0–1).
    pub catalyst_efficiency: f64,
}
impl EmissionsModel {
    /// Construct an emissions model for a naturally-aspirated gasoline engine.
    pub fn new_gasoline_na(
        displacement_l: f64,
        n_cylinders: u32,
        engine_rpm: f64,
        fuel_mass_mg: f64,
        lambda: f64,
        peak_temp_k: f64,
    ) -> Self {
        Self {
            displacement_l,
            n_cylinders,
            engine_rpm,
            fuel_mass_mg,
            lambda,
            peak_temp_k,
            egr_fraction: 0.0_f64,
            catalyst_efficiency: 0.97_f64,
        }
    }
    /// Engine firing rate (cycles/s) for a 4-stroke engine.
    pub fn firing_rate_hz(&self) -> f64 {
        self.engine_rpm / 60.0_f64 * self.n_cylinders as f64 / 2.0_f64
    }
    /// Raw (pre-catalyst) NOx mass flow (g/s) using the Zeldovich mechanism.
    ///
    /// `NOx_raw ∝ exp(-Ea / (R·T_peak)) × (1 - EGR)²`
    pub fn nox_raw_g_s(&self) -> f64 {
        if self.peak_temp_k < 1500.0_f64 {
            return 0.0_f64;
        }
        let ea_over_r = 38_000.0_f64;
        let c_nox = 8.0e-5_f64;
        let egr_factor = (1.0_f64 - self.egr_fraction).powi(2);
        let per_cycle = self.fuel_mass_mg * c_nox * (-ea_over_r / self.peak_temp_k).exp();
        per_cycle * self.firing_rate_hz() * 1.0e-6_f64 * egr_factor
    }
    /// Post-catalyst NOx mass flow (g/s).
    pub fn nox_tailpipe_g_s(&self) -> f64 {
        self.nox_raw_g_s() * (1.0_f64 - self.catalyst_efficiency).max(0.0_f64)
    }
    /// CO mass flow (g/s) — rich-burn driven.
    pub fn co_raw_g_s(&self) -> f64 {
        let co_fraction =
            (1.0_f64 - self.lambda.clamp(0.7_f64, 1.5_f64)).max(0.0_f64) * 0.15_f64 + 0.002_f64;
        self.fuel_mass_mg * co_fraction * self.firing_rate_hz() * 1.0e-6_f64
    }
    /// HC mass flow (g/s).
    pub fn hc_raw_g_s(&self) -> f64 {
        let hc_frac = 0.002_f64 + 0.01_f64 * (self.lambda - 1.0_f64).powi(2);
        self.fuel_mass_mg * hc_frac * self.firing_rate_hz() * 1.0e-6_f64
    }
    /// Particulate matter (PM) mass flow (g/s).
    ///
    /// PM is primarily relevant for GDI and diesel; approximated as a
    /// rich-mixture and temperature function.
    pub fn pm_raw_g_s(&self) -> f64 {
        if self.lambda >= 1.0_f64 {
            return 1.0e-7_f64 * self.fuel_mass_mg * self.firing_rate_hz() * 1.0e-6_f64;
        }
        let soot_factor = (1.0_f64 - self.lambda).powi(2) * 5.0e-4_f64;
        self.fuel_mass_mg * soot_factor * self.firing_rate_hz() * 1.0e-6_f64
    }
    /// Fuel consumption rate (g/s) — total fuel flow.
    pub fn fuel_flow_g_s(&self) -> f64 {
        self.fuel_mass_mg * self.firing_rate_hz() * 1.0e-3_f64
    }
    /// Brake-specific fuel consumption proxy (g/kWh).
    ///
    /// Estimated from fuel flow divided by a rough indicated power.
    /// `BSFC ≈ ṁ_fuel / (η_i × LHV × ṁ_fuel)` — simplified to constant-η estimate.
    pub fn bsfc_g_per_kwh(&self) -> f64 {
        let lhv_kj_per_g = 44.0_f64;
        let eta_indicated = 0.40_f64;
        3600.0_f64 / (lhv_kj_per_g * eta_indicated)
    }
}
/// Electric fuel pump model with pressure-flow curve and power consumption.
#[derive(Debug, Clone)]
pub struct FuelPump {
    /// Maximum delivery pressure (bar).
    pub max_pressure_bar: f64,
    /// Maximum volumetric flow rate at zero back-pressure (mL/min).
    pub max_flow_ml_per_min: f64,
    /// Pump volumetric efficiency (0–1).
    pub volumetric_efficiency: f64,
    /// Motor efficiency (0–1).
    pub motor_efficiency: f64,
    /// Supply voltage (V).
    pub supply_voltage_v: f64,
    /// Current draw at nominal operating point (A).
    pub nominal_current_a: f64,
}
impl FuelPump {
    /// Construct a typical GDI high-pressure pump.
    pub fn new_gdi_hp() -> Self {
        Self {
            max_pressure_bar: 250.0_f64,
            max_flow_ml_per_min: 400.0_f64,
            volumetric_efficiency: 0.92_f64,
            motor_efficiency: 0.85_f64,
            supply_voltage_v: 12.0_f64,
            nominal_current_a: 8.0_f64,
        }
    }
    /// Construct a low-pressure transfer pump (port injection / feed).
    pub fn new_transfer() -> Self {
        Self {
            max_pressure_bar: 6.0_f64,
            max_flow_ml_per_min: 1200.0_f64,
            volumetric_efficiency: 0.88_f64,
            motor_efficiency: 0.80_f64,
            supply_voltage_v: 12.0_f64,
            nominal_current_a: 4.5_f64,
        }
    }
    /// Flow rate (mL/min) at a given back-pressure `back_pressure_bar`.
    ///
    /// Linear pump curve: `Q = Q_max · (1 − back_pressure / max_pressure)`.
    pub fn flow_rate_ml_per_min(&self, back_pressure_bar: f64) -> f64 {
        let ratio = (back_pressure_bar / self.max_pressure_bar).clamp(0.0_f64, 1.0_f64);
        self.max_flow_ml_per_min * (1.0_f64 - ratio) * self.volumetric_efficiency
    }
    /// Hydraulic power delivered to the fuel (W).
    pub fn hydraulic_power_w(&self, back_pressure_bar: f64) -> f64 {
        let q_m3_s = self.flow_rate_ml_per_min(back_pressure_bar) * 1.0e-6_f64 / 60.0_f64;
        let dp_pa = back_pressure_bar * 1.0e5_f64;
        q_m3_s * dp_pa
    }
    /// Electric power consumed by the pump motor (W).
    pub fn electric_power_w(&self) -> f64 {
        self.supply_voltage_v * self.nominal_current_a
    }
    /// Overall pump efficiency at a given operating point.
    pub fn overall_efficiency(&self, back_pressure_bar: f64) -> f64 {
        let p_hyd = self.hydraulic_power_w(back_pressure_bar);
        let p_elec = self.electric_power_w();
        if p_elec < 1.0e-9_f64 {
            return 0.0_f64;
        }
        (p_hyd / p_elec).clamp(0.0_f64, 1.0_f64)
    }
}
/// Complete exhaust system: pipes, muffler, catalyst.
#[derive(Debug, Clone)]
pub struct ExhaustSystem {
    /// Catalyst module.
    pub catalyst: Catalyst,
    /// Exhaust pipe length (m).
    pub pipe_length_m: f64,
    /// Exhaust pipe diameter (m).
    pub pipe_diameter_m: f64,
    /// Exhaust gas temperature at manifold outlet (K).
    pub manifold_temp_k: f64,
    /// Ambient temperature (K).
    pub ambient_temp_k: f64,
    /// Backpressure at rated mass flow (Pa).
    pub rated_backpressure_pa: f64,
    /// Muffler sound attenuation (dB).
    pub muffler_attenuation_db: f64,
}
impl ExhaustSystem {
    /// Construct a standard passenger car exhaust system.
    pub fn new_passenger_car() -> Self {
        Self {
            catalyst: Catalyst::new_twc(),
            pipe_length_m: 2.5,
            pipe_diameter_m: 0.06,
            manifold_temp_k: 1100.0,
            ambient_temp_k: 293.0,
            rated_backpressure_pa: 8_000.0,
            muffler_attenuation_db: 20.0,
        }
    }
    /// Exhaust gas temperature after the pipe (accounts for cooling).
    ///
    /// Uses a simple exponential cooling model.
    pub fn pipe_outlet_temp_k(&self) -> f64 {
        let t_in = self.manifold_temp_k;
        let t_amb = self.ambient_temp_k;
        let cooling_per_m = 30.0;
        (t_in - t_amb) * (-cooling_per_m * self.pipe_length_m / (t_in - t_amb + 1.0)).exp() + t_amb
    }
    /// Dynamic backpressure (Pa) at given mass flow rate (kg/s).
    ///
    /// Backpressure scales quadratically with flow velocity.
    pub fn backpressure_pa(&self, mass_flow_kg_s: f64) -> f64 {
        let rated_flow = 0.05;
        if rated_flow < 1e-15 {
            return 0.0;
        }
        let ratio = mass_flow_kg_s / rated_flow;
        self.rated_backpressure_pa * ratio * ratio
    }
    /// Post-catalyst emission factor \[0.0, 1.0\] for a given lambda.
    ///
    /// 1.0 = all emissions pass through; 0.0 = fully converted.
    pub fn emission_pass_through(&self, lambda: f64) -> f64 {
        1.0 - self.catalyst.conversion_efficiency(lambda)
    }
    /// Update the full thermal model for one time step `dt`.
    pub fn update(&mut self, dt: f64) {
        let t_outlet = self.pipe_outlet_temp_k();
        self.catalyst.update_temperature(t_outlet, dt);
    }
    /// Whether the catalyst has reached light-off temperature.
    pub fn catalyst_lit_off(&self) -> bool {
        self.catalyst.temperature_k >= self.catalyst.light_off_k
    }
}
/// Fuel tank with sloshing dynamics modelled as a pendulum.
///
/// The sloshing model represents the centre-of-gravity shift due to fuel
/// motion as a 1-D damped pendulum excited by lateral acceleration.
#[derive(Debug, Clone)]
pub struct FuelTank {
    /// Maximum tank capacity (litres).
    pub capacity_l: f64,
    /// Current fuel level (litres).
    pub level_l: f64,
    /// Fuel density (kg/L).
    pub density_kg_l: f64,
    /// Tank width (m) — characteristic dimension for sloshing.
    pub width_m: f64,
    /// Tank height from floor (m).
    pub height_m: f64,
    /// Sloshing pendulum angle (rad).
    pub slosh_angle: f64,
    /// Sloshing pendulum angular velocity (rad/s).
    pub slosh_ang_vel: f64,
    /// Sloshing damping coefficient (dimensionless).
    pub slosh_damping: f64,
}
impl FuelTank {
    /// Construct a fuel tank.
    ///
    /// # Arguments
    /// * `capacity_l` — maximum capacity in litres.
    /// * `level_fraction` — initial fill fraction \[0.0, 1.0\].
    /// * `density_kg_l` — fuel density (0.745 gasoline, 0.832 diesel).
    /// * `width_m` — tank width for sloshing pendulum length.
    pub fn new(capacity_l: f64, level_fraction: f64, density_kg_l: f64, width_m: f64) -> Self {
        Self {
            capacity_l,
            level_l: capacity_l * level_fraction.clamp(0.0, 1.0),
            density_kg_l,
            width_m,
            height_m: 0.15,
            slosh_angle: 0.0,
            slosh_ang_vel: 0.0,
            slosh_damping: 0.2,
        }
    }
    /// Construct a standard gasoline tank.
    pub fn new_gasoline(capacity_l: f64, level_fraction: f64) -> Self {
        Self::new(capacity_l, level_fraction, 0.745, 0.4)
    }
    /// Construct a standard diesel tank.
    pub fn new_diesel(capacity_l: f64, level_fraction: f64) -> Self {
        Self::new(capacity_l, level_fraction, 0.832, 0.4)
    }
    /// Current fuel mass (kg).
    pub fn mass_kg(&self) -> f64 {
        self.level_l * self.density_kg_l
    }
    /// Current fill fraction \[0.0, 1.0\].
    pub fn fill_fraction(&self) -> f64 {
        if self.capacity_l > 0.0 {
            (self.level_l / self.capacity_l).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }
    /// Consume `mass_kg` of fuel.
    ///
    /// Returns the amount actually consumed (limited by available fuel).
    pub fn consume(&mut self, mass_kg: f64) -> f64 {
        let available_kg = self.mass_kg();
        let consumed_kg = mass_kg.clamp(0.0, available_kg);
        let consumed_l = consumed_kg / self.density_kg_l;
        self.level_l = (self.level_l - consumed_l).max(0.0);
        consumed_kg
    }
    /// Refuel by adding `volume_l` litres.
    ///
    /// Clamps to capacity.
    pub fn refuel(&mut self, volume_l: f64) {
        self.level_l = (self.level_l + volume_l).min(self.capacity_l);
    }
    /// Integrate the sloshing pendulum for one time step `dt`.
    ///
    /// * `lateral_accel` — lateral acceleration of the vehicle (m/s²).
    ///
    /// The pendulum length is proportional to the fill fraction × tank width.
    pub fn integrate_slosh(&mut self, lateral_accel: f64, dt: f64) {
        let g = 9.81;
        let fill = self.fill_fraction();
        if fill < 1e-3 {
            return;
        }
        let l = 0.5 * self.width_m * fill;
        let omega_n_sq = g / l;
        let omega_n = omega_n_sq.sqrt();
        let theta_ddot = -omega_n_sq * self.slosh_angle
            - 2.0 * self.slosh_damping * omega_n * self.slosh_ang_vel
            - lateral_accel / l;
        self.slosh_ang_vel += theta_ddot * dt;
        self.slosh_angle += self.slosh_ang_vel * dt;
        self.slosh_angle = self.slosh_angle.clamp(-PI / 3.0, PI / 3.0);
    }
    /// Centre-of-gravity lateral offset due to sloshing (m).
    pub fn slosh_cg_offset(&self) -> f64 {
        0.5 * self.width_m * self.slosh_angle.sin() * self.fill_fraction()
    }
    /// Whether the tank is empty.
    pub fn is_empty(&self) -> bool {
        self.level_l < 1e-6
    }
}
/// Unified energy recovery system combining KERS, TERS and regenerative braking.
#[derive(Debug, Clone)]
pub struct EnergyRecovery {
    /// Flywheel KERS.
    pub kers: KersSystem,
    /// Thermal recovery TERS.
    pub ters: TersSystem,
    /// Regenerative braking.
    pub regen_braking: RegenerativeBraking,
    /// Battery to store recovered electrical energy.
    pub battery: HybridBattery,
}
impl EnergyRecovery {
    /// Construct a combined energy recovery system.
    pub fn new() -> Self {
        Self {
            kers: KersSystem::new(0.25, 60_000.0 * PI / 30.0),
            ters: TersSystem::new_teg(1000.0),
            regen_braking: RegenerativeBraking::new(200.0, 0.6, 0.85),
            battery: HybridBattery::new_nmc(1.5, 0.5),
        }
    }
    /// Total available stored energy (J).
    pub fn total_stored_energy_j(&self) -> f64 {
        self.kers.stored_energy_j + self.battery.available_energy_kwh() * 3.6e6
    }
    /// Step all recovery systems for `dt` seconds.
    ///
    /// * `brake_torque_nm` — wheel braking torque.
    /// * `wheel_speed_rad_s` — wheel angular speed.
    /// * `exhaust_mdot_kg_s` — exhaust mass flow rate.
    pub fn step(
        &mut self,
        brake_torque_nm: f64,
        wheel_speed_rad_s: f64,
        exhaust_mdot_kg_s: f64,
        dt: f64,
    ) {
        self.regen_braking
            .update(brake_torque_nm, wheel_speed_rad_s);
        let regen_w = self.regen_braking.regen_power_w;
        self.battery.charge(regen_w / 1000.0, dt);
        let ters_w = self.ters.recoverable_power_w(exhaust_mdot_kg_s);
        self.battery.charge(ters_w / 1000.0, dt);
        let heat_w = self
            .battery
            .ohmic_heat_w(regen_w / self.battery.nominal_voltage_v);
        self.battery.update_thermal(heat_w, 25.0, dt);
    }
}
/// Combustion model: heat release, stoichiometry and emissions.
#[derive(Debug, Clone)]
pub struct CombustionModel {
    /// Fuel properties.
    pub fuel: FuelProperties,
    /// Fuel mass per cycle (mg).
    pub fuel_mass_mg: f64,
    /// Air mass per cycle (mg).
    pub air_mass_mg: f64,
    /// Engine speed (rpm).
    pub engine_rpm: f64,
    /// Combustion efficiency \[0.0, 1.0\].
    pub combustion_efficiency: f64,
    /// Wiebe function exponent m (shape parameter).
    pub wiebe_m: f64,
    /// Wiebe function combustion duration (°CA).
    pub wiebe_duration_deg: f64,
    /// Wiebe function start of combustion (°CA aTDC).
    pub wiebe_start_deg: f64,
}
impl CombustionModel {
    /// Construct a combustion model for gasoline SI engines.
    pub fn new_gasoline_si(engine_rpm: f64, fuel_mass_mg: f64) -> Self {
        let fuel = FuelProperties::gasoline();
        let air_mass_mg = fuel_mass_mg * fuel.stoich_afr;
        Self {
            fuel,
            fuel_mass_mg,
            air_mass_mg,
            engine_rpm,
            combustion_efficiency: 0.95,
            wiebe_m: 2.0,
            wiebe_duration_deg: 60.0,
            wiebe_start_deg: -10.0,
        }
    }
    /// Equivalence ratio λ = (actual AFR) / (stoich AFR).
    pub fn lambda(&self) -> f64 {
        if self.fuel_mass_mg < 1e-15 {
            return f64::INFINITY;
        }
        let actual_afr = self.air_mass_mg / self.fuel_mass_mg;
        actual_afr / self.fuel.stoich_afr
    }
    /// Total heat released per cycle (J).
    pub fn heat_release_j(&self) -> f64 {
        let fuel_kg = self.fuel_mass_mg * 1e-6;
        let lhv_j_per_kg = self.fuel.lhv_mj_per_kg * 1e6;
        fuel_kg * lhv_j_per_kg * self.combustion_efficiency
    }
    /// Wiebe heat-release fraction at crank angle `theta_ca` (°CA aTDC).
    ///
    /// x_b(θ) = 1 − exp(−a ((θ − θ_s) / Δθ)^(m+1))
    /// with a = 5.0 (Wiebe efficiency factor).
    pub fn wiebe_fraction(&self, theta_ca: f64) -> f64 {
        let a = 5.0;
        let xi = (theta_ca - self.wiebe_start_deg) / self.wiebe_duration_deg;
        if xi <= 0.0 {
            return 0.0;
        }
        if xi >= 1.0 {
            return 1.0;
        }
        1.0 - (-a * xi.powf(self.wiebe_m + 1.0)).exp()
    }
    /// Instantaneous heat release rate (J/°CA) at crank angle `theta`.
    pub fn heat_release_rate(&self, theta_ca: f64) -> f64 {
        let a = 5.0;
        let xi = (theta_ca - self.wiebe_start_deg) / self.wiebe_duration_deg;
        if xi <= 0.0 || xi >= 1.0 {
            return 0.0;
        }
        let q_total = self.heat_release_j();
        q_total * a * (self.wiebe_m + 1.0) / self.wiebe_duration_deg
            * xi.powf(self.wiebe_m)
            * (-a * xi.powf(self.wiebe_m + 1.0)).exp()
    }
    /// Estimate peak cylinder temperature (K) for NOx model.
    ///
    /// Simple adiabatic flame temperature approximation.
    pub fn peak_temperature_k(&self) -> f64 {
        let t_intake = 350.0;
        let q = self.heat_release_j();
        let m_charge = (self.fuel_mass_mg + self.air_mass_mg) * 1e-6;
        let cp = 1000.0;
        if m_charge < 1e-15 {
            return t_intake;
        }
        t_intake + q / (m_charge * cp)
    }
    /// NOx emission (mg/cycle) using a simplified Zeldovich rate model.
    ///
    /// NOx ∝ exp(−Ea / (R T_peak)) for T_peak > threshold.
    pub fn nox_mg_per_cycle(&self) -> f64 {
        let t_peak = self.peak_temperature_k();
        let t_threshold = 1800.0;
        if t_peak < t_threshold {
            return 0.0;
        }
        let ea_over_r = 38_000.0;
        let c_nox = 1.5e-3;
        let fuel_mg = self.fuel_mass_mg;
        fuel_mg * c_nox * (-ea_over_r / t_peak).exp()
    }
    /// CO emission (mg/cycle): inverse function of lambda (richer = more CO).
    pub fn co_mg_per_cycle(&self) -> f64 {
        let lam = self.lambda().clamp(0.7, 1.5);
        let fuel_mg = self.fuel_mass_mg;
        let co_fraction = (1.0 - lam).max(0.0) * 0.15 + 0.002;
        fuel_mg * co_fraction
    }
    /// HC emission (mg/cycle): increases with lean misfire and rich combustion.
    pub fn hc_mg_per_cycle(&self) -> f64 {
        let lam = self.lambda().clamp(0.7, 2.0);
        let fuel_mg = self.fuel_mass_mg;
        let hc_frac = 0.002 + 0.01 * (lam - 1.0).powi(2);
        fuel_mg * hc_frac * (1.0 - self.combustion_efficiency)
    }
}
/// Fuel tank centre-of-mass model that accounts for varying fill level.
///
/// The tank is modelled as a rectangular box; CoM shifts vertically as
/// fuel level changes and laterally due to sloshing.
#[derive(Debug, Clone)]
pub struct FuelTankCoM {
    /// Tank body reference CoM in the vehicle frame (m).
    pub tank_com_ref: [f64; 3],
    /// Tank internal height (m).
    pub tank_internal_height_m: f64,
    /// Tank width (m).
    pub tank_width_m: f64,
    /// Tank length (m).
    pub tank_length_m: f64,
    /// Fuel density (kg/L).
    pub density_kg_l: f64,
    /// Maximum capacity (L).
    pub capacity_l: f64,
    /// Current fill level (L).
    pub fill_l: f64,
    /// Sloshing lateral offset factor (dimensionless, typically ±0.2).
    pub slosh_lateral_factor: f64,
}
impl FuelTankCoM {
    /// Construct a [`FuelTankCoM`] model.
    pub fn new(
        tank_com_ref: [f64; 3],
        tank_internal_height_m: f64,
        tank_width_m: f64,
        tank_length_m: f64,
        density_kg_l: f64,
        capacity_l: f64,
        fill_fraction: f64,
    ) -> Self {
        Self {
            tank_com_ref,
            tank_internal_height_m,
            tank_width_m,
            tank_length_m,
            density_kg_l,
            capacity_l,
            fill_l: capacity_l * fill_fraction.clamp(0.0_f64, 1.0_f64),
            slosh_lateral_factor: 0.0_f64,
        }
    }
    /// Fill fraction \[0, 1\].
    pub fn fill_fraction(&self) -> f64 {
        if self.capacity_l < 1.0e-9_f64 {
            return 0.0_f64;
        }
        (self.fill_l / self.capacity_l).clamp(0.0_f64, 1.0_f64)
    }
    /// Fuel mass (kg).
    pub fn fuel_mass_kg(&self) -> f64 {
        self.fill_l * self.density_kg_l
    }
    /// Centre of mass of the fuel in the vehicle frame (m).
    ///
    /// Vertical CoM of the fuel column = tank bottom + fill_height/2.
    /// Lateral CoM shifted by sloshing.
    pub fn fuel_com(&self) -> [f64; 3] {
        let fill_height = self.tank_internal_height_m * self.fill_fraction();
        let bottom_z = self.tank_com_ref[2] - 0.5_f64 * self.tank_internal_height_m;
        let fuel_z = bottom_z + 0.5_f64 * fill_height;
        let fuel_y = self.tank_com_ref[1]
            + self.tank_width_m * self.slosh_lateral_factor * self.fill_fraction();
        [self.tank_com_ref[0], fuel_y, fuel_z]
    }
    /// Combined (tank structure + fuel) CoM given tank structure mass `tank_mass_kg`.
    pub fn combined_com(&self, tank_mass_kg: f64) -> [f64; 3] {
        let m_fuel = self.fuel_mass_kg();
        let m_total = tank_mass_kg + m_fuel;
        if m_total < 1.0e-9_f64 {
            return self.tank_com_ref;
        }
        let f_com = self.fuel_com();
        [
            (tank_mass_kg * self.tank_com_ref[0] + m_fuel * f_com[0]) / m_total,
            (tank_mass_kg * self.tank_com_ref[1] + m_fuel * f_com[1]) / m_total,
            (tank_mass_kg * self.tank_com_ref[2] + m_fuel * f_com[2]) / m_total,
        ]
    }
    /// Consume `mass_kg` of fuel and update fill level.
    pub fn consume(&mut self, mass_kg: f64) {
        let consumed_l = (mass_kg / self.density_kg_l.max(1.0e-9_f64)).clamp(0.0_f64, self.fill_l);
        self.fill_l -= consumed_l;
    }
}
/// Thermal Energy Recovery System (TERS) — Rankine cycle / thermoelectric.
#[derive(Debug, Clone)]
pub struct TersSystem {
    /// Exhaust temperature (K).
    pub exhaust_temp_k: f64,
    /// Cold side temperature (K).
    pub cold_temp_k: f64,
    /// Heat exchanger effectiveness \[0.0, 1.0\].
    pub hx_effectiveness: f64,
    /// Maximum recoverable power (W).
    pub max_power_w: f64,
    /// Thermoelectric module Seebeck coefficient (V/K) — for TEG model.
    pub seebeck_coefficient: f64,
    /// Module internal resistance (Ω).
    pub module_resistance_ohm: f64,
}
impl TersSystem {
    /// Construct a thermoelectric generator (TEG) module.
    pub fn new_teg(max_power_w: f64) -> Self {
        Self {
            exhaust_temp_k: 800.0,
            cold_temp_k: 300.0,
            hx_effectiveness: 0.7,
            max_power_w,
            seebeck_coefficient: 1.5e-3,
            module_resistance_ohm: 0.5,
        }
    }
    /// Carnot efficiency limit for current hot/cold temperatures.
    pub fn carnot_efficiency(&self) -> f64 {
        if self.exhaust_temp_k <= self.cold_temp_k {
            return 0.0;
        }
        1.0 - self.cold_temp_k / self.exhaust_temp_k
    }
    /// Recoverable power (W) given mass flow `mdot_kg_s` and exhaust cp.
    pub fn recoverable_power_w(&self, mdot_kg_s: f64) -> f64 {
        let cp = 1100.0;
        let q_available =
            mdot_kg_s * cp * (self.exhaust_temp_k - self.cold_temp_k) * self.hx_effectiveness;
        let eta = self.carnot_efficiency() * 0.3;
        (q_available * eta).min(self.max_power_w).max(0.0)
    }
    /// TEG open-circuit voltage (V).
    pub fn teg_open_circuit_voltage(&self) -> f64 {
        self.seebeck_coefficient * (self.exhaust_temp_k - self.cold_temp_k)
    }
    /// TEG maximum power point power (W).
    pub fn teg_max_power_w(&self) -> f64 {
        let v_oc = self.teg_open_circuit_voltage();
        v_oc * v_oc / (4.0 * self.module_resistance_ohm)
    }
}
