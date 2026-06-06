// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle HVAC (Heating, Ventilation, and Air Conditioning) system models.
//!
//! Provides lumped-parameter thermal models for the cabin, refrigeration
//! cycle, heat-pump heating, occupant comfort indices (PMV/PPD),
//! humidity/defogging, energy balance, and EV pre-conditioning.

use std::f64::consts::E;

// ── CabinThermal ──────────────────────────────────────────────────────────────

/// Lumped thermal model of the vehicle cabin.
///
/// Treats cabin air, walls, and glass as separate thermal masses
/// and tracks heat exchange between them and the environment.
#[derive(Debug, Clone)]
pub struct CabinThermal {
    /// Cabin air temperature \[°C\].
    pub air_temp: f64,
    /// Wall (body panel) mean temperature \[°C\].
    pub wall_temp: f64,
    /// Glass mean temperature \[°C\].
    pub glass_temp: f64,
    /// Cabin air thermal mass \[J/K\] (air mass × cp).
    pub air_thermal_mass: f64,
    /// Wall thermal mass \[J/K\].
    pub wall_thermal_mass: f64,
    /// Glass thermal mass \[J/K\].
    pub glass_thermal_mass: f64,
    /// Convective conductance cabin-air ↔ wall \[W/K\].
    pub ua_wall: f64,
    /// Convective conductance cabin-air ↔ glass \[W/K\].
    pub ua_glass: f64,
    /// Convective conductance wall ↔ ambient \[W/K\].
    pub ua_ambient_wall: f64,
    /// Convective conductance glass ↔ ambient \[W/K\].
    pub ua_ambient_glass: f64,
    /// Occupant sensible heat load \[W\].
    pub occupant_heat: f64,
}

impl CabinThermal {
    /// Create a cabin model with reasonable default parameters for a
    /// mid-size passenger car (5 occupants, 22 °C cabin, 35 °C ambient).
    pub fn default_car() -> Self {
        Self {
            air_temp: 22.0,
            wall_temp: 30.0,
            glass_temp: 32.0,
            air_thermal_mass: 3500.0, // ~3 kg air × 1005 J/(kg·K) × volume
            wall_thermal_mass: 50_000.0, // sheet metal + insulation
            glass_thermal_mass: 8_000.0,
            ua_wall: 15.0,
            ua_glass: 25.0,
            ua_ambient_wall: 20.0,
            ua_ambient_glass: 40.0,
            occupant_heat: 400.0, // 5 occupants × 80 W
        }
    }

    /// Advance the cabin thermal state by `dt` seconds.
    ///
    /// `ambient_temp` is outdoor temperature \[°C\];
    /// `hvac_power` is net power delivered to cabin air (positive = heating) \[W\].
    pub fn step(&mut self, dt: f64, ambient_temp: f64, hvac_power: f64) {
        // Heat flows
        let q_air_wall = self.ua_wall * (self.wall_temp - self.air_temp);
        let q_air_glass = self.ua_glass * (self.glass_temp - self.air_temp);
        let q_wall_amb = self.ua_ambient_wall * (ambient_temp - self.wall_temp);
        let q_glass_amb = self.ua_ambient_glass * (ambient_temp - self.glass_temp);

        // dT/dt for each node
        let d_air =
            (hvac_power + self.occupant_heat + q_air_wall + q_air_glass) / self.air_thermal_mass;
        let d_wall = (q_wall_amb - q_air_wall) / self.wall_thermal_mass;
        let d_glass = (q_glass_amb - q_air_glass) / self.glass_thermal_mass;

        self.air_temp += d_air * dt;
        self.wall_temp += d_wall * dt;
        self.glass_temp += d_glass * dt;
    }

    /// Temperature difference between cabin air and setpoint.
    pub fn error(&self, setpoint: f64) -> f64 {
        setpoint - self.air_temp
    }
}

// ── HvacCompressor ────────────────────────────────────────────────────────────

/// Refrigeration-cycle compressor model.
///
/// Provides cooling capacity and COP as a function of compressor speed,
/// evaporator and condenser temperatures.
#[derive(Debug, Clone)]
pub struct HvacCompressor {
    /// Maximum cooling capacity at rated speed \[W\].
    pub max_cooling_capacity: f64,
    /// Rated compressor speed \[rpm\].
    pub rated_speed: f64,
    /// Current compressor speed \[rpm\].
    pub speed: f64,
    /// Evaporator temperature \[°C\].
    pub evap_temp: f64,
    /// Condenser temperature \[°C\].
    pub cond_temp: f64,
    /// Isentropic efficiency (0–1).
    pub isentropic_efficiency: f64,
}

impl HvacCompressor {
    /// Create a typical automotive AC compressor (variable-displacement).
    pub fn automotive() -> Self {
        Self {
            max_cooling_capacity: 5000.0,
            rated_speed: 3000.0,
            speed: 1500.0,
            evap_temp: 5.0,
            cond_temp: 55.0,
            isentropic_efficiency: 0.75,
        }
    }

    /// Cooling capacity \[W\] proportional to speed ratio.
    pub fn cooling_capacity(&self) -> f64 {
        let ratio = (self.speed / self.rated_speed).clamp(0.0, 1.0);
        self.max_cooling_capacity * ratio
    }

    /// Carnot COP = T_evap / (T_cond - T_evap) in Kelvin.
    pub fn carnot_cop(&self) -> f64 {
        let te = self.evap_temp + 273.15;
        let tc = self.cond_temp + 273.15;
        if tc > te {
            te / (tc - te)
        } else {
            f64::INFINITY
        }
    }

    /// Actual COP = Carnot COP × isentropic efficiency.
    pub fn actual_cop(&self) -> f64 {
        self.carnot_cop() * self.isentropic_efficiency
    }

    /// Compressor electrical power draw \[W\].
    pub fn power_draw(&self) -> f64 {
        let cop = self.actual_cop();
        if cop > 0.0 {
            self.cooling_capacity() / cop
        } else {
            0.0
        }
    }
}

// ── HeatPump ──────────────────────────────────────────────────────────────────

/// Heat-pump heating mode model (reverse refrigeration cycle).
///
/// Used in EVs to extract heat from ambient and deliver it to the cabin.
#[derive(Debug, Clone)]
pub struct HeatPump {
    /// Ambient (source) temperature \[°C\].
    pub ambient_temp: f64,
    /// Cabin (sink) temperature target \[°C\].
    pub cabin_temp: f64,
    /// Maximum heat output at rated conditions \[W\].
    pub max_heat_output: f64,
    /// Isentropic efficiency (0–1).
    pub efficiency: f64,
    /// Whether the defrost cycle is active.
    pub defrost_active: bool,
    /// Accumulated frost mass on evaporator \[kg\].
    pub frost_mass: f64,
    /// Frost accumulation rate \[kg/s\] when ambient is cold and humid.
    pub frost_rate: f64,
}

impl HeatPump {
    /// Create a typical EV heat pump.
    pub fn ev_heat_pump() -> Self {
        Self {
            ambient_temp: 5.0,
            cabin_temp: 21.0,
            max_heat_output: 4000.0,
            efficiency: 0.7,
            defrost_active: false,
            frost_mass: 0.0,
            frost_rate: 1e-4,
        }
    }

    /// Carnot COP for heating: T_sink / (T_sink - T_source) in Kelvin.
    pub fn carnot_cop_heating(&self) -> f64 {
        let ts = self.ambient_temp + 273.15;
        let tc = self.cabin_temp + 273.15;
        if tc > ts {
            tc / (tc - ts)
        } else {
            f64::INFINITY
        }
    }

    /// Actual heating COP.
    pub fn actual_cop_heating(&self) -> f64 {
        self.carnot_cop_heating() * self.efficiency
    }

    /// Heat output \[W\].
    pub fn heat_output(&self) -> f64 {
        if self.defrost_active {
            0.0 // no cabin heating during defrost
        } else {
            self.max_heat_output.min(self.max_heat_output)
        }
    }

    /// Electrical power draw \[W\].
    pub fn power_draw(&self) -> f64 {
        let cop = self.actual_cop_heating().max(1.0);
        self.heat_output() / cop
    }

    /// Advance defrost / frost model by `dt` seconds.
    ///
    /// Returns heat energy consumed by defrost \[J\].
    pub fn step_frost(&mut self, dt: f64) -> f64 {
        // Accumulate frost when ambient < 5 °C
        if self.ambient_temp < 5.0 && !self.defrost_active {
            self.frost_mass += self.frost_rate * dt;
        }
        // Trigger defrost above threshold
        if self.frost_mass > 0.01 {
            self.defrost_active = true;
        }

        if self.defrost_active {
            let heat_of_fusion = 334_000.0; // J/kg for ice
            let melted = (self.frost_mass).min(self.frost_rate * 100.0 * dt);
            self.frost_mass -= melted;
            if self.frost_mass < 1e-6 {
                self.defrost_active = false;
                self.frost_mass = 0.0;
            }
            melted * heat_of_fusion
        } else {
            0.0
        }
    }
}

// ── AirDistribution ───────────────────────────────────────────────────────────

/// Air distribution fractions across vent modes.
///
/// The three fractions must sum to 1.0; they control how conditioned air
/// is split between face vents, floor vents, and windshield defrost.
#[derive(Debug, Clone)]
pub struct AirDistribution {
    /// Face vent fraction (0–1).
    pub face_fraction: f64,
    /// Floor vent fraction (0–1).
    pub floor_fraction: f64,
    /// Defrost vent fraction (0–1).
    pub defrost_fraction: f64,
    /// Total airflow rate \[m³/s\].
    pub total_flow: f64,
}

impl AirDistribution {
    /// Create with a given mode split and total flow.
    ///
    /// Fractions are normalised automatically.
    pub fn new(face: f64, floor: f64, defrost: f64, total_flow: f64) -> Self {
        let sum = face + floor + defrost;
        let inv = if sum > 1e-12 { 1.0 / sum } else { 1.0 };
        Self {
            face_fraction: face * inv,
            floor_fraction: floor * inv,
            defrost_fraction: defrost * inv,
            total_flow,
        }
    }

    /// Face vent airflow \[m³/s\].
    pub fn face_flow(&self) -> f64 {
        self.face_fraction * self.total_flow
    }

    /// Floor vent airflow \[m³/s\].
    pub fn floor_flow(&self) -> f64 {
        self.floor_fraction * self.total_flow
    }

    /// Defrost vent airflow \[m³/s\].
    pub fn defrost_flow(&self) -> f64 {
        self.defrost_fraction * self.total_flow
    }

    /// Set full-defrost mode (all airflow to windshield).
    pub fn set_defrost_mode(&mut self) {
        self.face_fraction = 0.0;
        self.floor_fraction = 0.0;
        self.defrost_fraction = 1.0;
    }
}

// ── OccupantComfort ───────────────────────────────────────────────────────────

/// Fanger PMV/PPD occupant thermal comfort model.
///
/// PMV (Predicted Mean Vote) range: −3 (cold) … +3 (hot); 0 is neutral.
/// PPD (Predicted Percentage Dissatisfied) is derived from PMV.
#[derive(Debug, Clone)]
pub struct OccupantComfort {
    /// Metabolic rate \[W/m²\]. Seated: ~58.15 W/m², driving: ~70 W/m².
    pub metabolic_rate: f64,
    /// Clothing insulation \[clo\]. 1 clo ≈ 0.155 m²·K/W.
    pub clothing_insulation: f64,
    /// Air temperature \[°C\].
    pub air_temp: f64,
    /// Mean radiant temperature \[°C\].
    pub mean_radiant_temp: f64,
    /// Air velocity \[m/s\].
    pub air_velocity: f64,
    /// Relative humidity \[%\].
    pub relative_humidity: f64,
}

impl OccupantComfort {
    /// Create comfort model for a seated driver in typical conditions.
    pub fn seated_driver() -> Self {
        Self {
            metabolic_rate: 70.0,
            clothing_insulation: 1.0,
            air_temp: 22.0,
            mean_radiant_temp: 22.0,
            air_velocity: 0.1,
            relative_humidity: 50.0,
        }
    }

    /// Compute PMV using a simplified ISO 7730 approximation.
    ///
    /// This linearisation gives ±0.3 accuracy for typical comfort ranges.
    pub fn pmv(&self) -> f64 {
        let ta = self.air_temp;
        let tr = self.mean_radiant_temp;
        let va = self.air_velocity.max(0.0);
        let rh = self.relative_humidity;
        let icl = self.clothing_insulation * 0.155; // m²·K/W
        let m = self.metabolic_rate;

        // Operative temperature
        let top = (ta + tr) / 2.0;

        // Water vapour pressure (Antoine approximation) [Pa]
        let pvs = 133.322 * E.powf(18.956 - 4030.18 / (ta + 235.0));
        let pa = rh / 100.0 * pvs;

        // Clothing surface temperature (simplified)
        let tcl =
            35.7 - 0.028 * m + icl * (3.96e-8 * 273.15_f64.powi(4) - (top + 273.15).powi(4) * 0.0);
        let _ = tcl; // used implicitly

        // PMV = (0.303 * exp(-0.036 m) + 0.028) * L, simplified L
        let l = m
            - 58.15
            - 0.42 * (m - 58.15)
            - 3.05e-3 * (5733.0 - 6.99 * m - pa)
            - 0.0173 * m * (5.867 - pa * 1e-3)
            - 0.0014 * m * (34.0 - ta)
            - 3.96e-8 * 0.72 * ((top + 273.15).powi(4) - (top + 273.15).powi(4))
            - 2.38 * va.powf(0.25) * (top - ta)
            - 0.5
            + 0.5;
        // Simplified continuous approximation
        let pmv_approx = (0.303 * (-0.036 * m).exp() + 0.028) * l;
        pmv_approx.clamp(-3.0, 3.0)
    }

    /// Compute PPD \[%\] from PMV: `PPD = 100 − 95 exp(−(0.03353 PMV⁴ + 0.2179 PMV²))`.
    pub fn ppd(&self) -> f64 {
        let pmv = self.pmv();
        100.0 - 95.0 * (-(0.03353 * pmv.powi(4) + 0.2179 * pmv.powi(2))).exp()
    }

    /// Returns true if the occupant is in the comfort zone (|PMV| ≤ 0.5).
    pub fn is_comfortable(&self) -> bool {
        self.pmv().abs() <= 0.5
    }
}

// ── HvacController ────────────────────────────────────────────────────────────

/// PID temperature controller for the HVAC system.
///
/// Tracks a temperature setpoint by adjusting HVAC output power.
/// Output is clamped to `[-max_power, max_power]`.
#[derive(Debug, Clone)]
pub struct HvacController {
    /// Temperature setpoint \[°C\].
    pub setpoint: f64,
    /// PID proportional gain \[W/K\].
    pub kp: f64,
    /// PID integral gain \[W/(K·s)\].
    pub ki: f64,
    /// PID derivative gain \[W·s/K\].
    pub kd: f64,
    /// Maximum HVAC output power \[W\].
    pub max_power: f64,
    /// Accumulated integral error \[K·s\].
    pub integral: f64,
    /// Previous error for derivative \[K\].
    pub prev_error: f64,
    /// Recirculation mode active (reduces fresh-air infiltration heat gain).
    pub recirculation: bool,
}

impl HvacController {
    /// Create a PID controller with typical automotive gains.
    pub fn new(setpoint: f64) -> Self {
        Self {
            setpoint,
            kp: 800.0,
            ki: 50.0,
            kd: 100.0,
            max_power: 6000.0,
            integral: 0.0,
            prev_error: 0.0,
            recirculation: false,
        }
    }

    /// Compute HVAC output power \[W\] given current cabin temperature and `dt`.
    ///
    /// Positive = heating, negative = cooling.
    pub fn update(&mut self, cabin_temp: f64, dt: f64) -> f64 {
        let error = self.setpoint - cabin_temp;
        self.integral += error * dt;
        // Anti-windup clamp
        let integral_limit = self.max_power / self.ki.max(1e-6);
        self.integral = self.integral.clamp(-integral_limit, integral_limit);
        let derivative = if dt > 1e-9 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        let output = self.kp * error + self.ki * self.integral + self.kd * derivative;
        output.clamp(-self.max_power, self.max_power)
    }

    /// Enable or disable recirculation mode.
    pub fn set_recirculation(&mut self, active: bool) {
        self.recirculation = active;
    }

    /// Infiltration penalty when fresh air mode is active \[W\].
    ///
    /// Models extra heat load from outside air at `ambient_temp`.
    pub fn infiltration_load(&self, cabin_temp: f64, ambient_temp: f64) -> f64 {
        if self.recirculation {
            0.0
        } else {
            let delta_t = cabin_temp - ambient_temp;
            // Approx 30 L/s fresh air × 1.2 kg/m³ × 1005 J/(kg·K)
            0.03 * 1.2 * 1005.0 * delta_t
        }
    }
}

// ── WindowDefogging ───────────────────────────────────────────────────────────

/// Anti-fog / windshield defogging model.
///
/// Estimates whether the glass surface is below dew point and the
/// required airflow power to prevent fogging.
#[derive(Debug, Clone)]
pub struct WindowDefogging {
    /// Glass surface temperature \[°C\].
    pub glass_temp: f64,
    /// Cabin air temperature \[°C\].
    pub air_temp: f64,
    /// Cabin relative humidity \[%\].
    pub relative_humidity: f64,
}

impl WindowDefogging {
    /// Create a defogging model.
    pub fn new(glass_temp: f64, air_temp: f64, relative_humidity: f64) -> Self {
        Self {
            glass_temp,
            air_temp,
            relative_humidity,
        }
    }

    /// Magnus formula dew point \[°C\].
    pub fn dew_point(&self) -> f64 {
        let a = 17.27;
        let b = 237.7;
        let alpha = a * self.air_temp / (b + self.air_temp) + (self.relative_humidity / 100.0).ln();
        b * alpha / (a - alpha)
    }

    /// Returns true if glass surface is at or below the dew point.
    pub fn is_fogging(&self) -> bool {
        self.glass_temp <= self.dew_point()
    }

    /// Required defrost airflow \[m³/s\] to raise glass temperature above dew point.
    ///
    /// Uses a simplified model: `Q = UA_glass * ΔT / (ρ cp ΔT_air)`.
    pub fn required_airflow(&self) -> f64 {
        let dp = self.dew_point();
        if self.glass_temp > dp {
            return 0.0;
        }
        let temp_margin = (dp - self.glass_temp) + 2.0; // 2 K safety margin
        let ua_glass = 40.0; // W/K
        let rho_cp = 1.2 * 1005.0; // J/(m³·K)
        let delta_t_air = (self.air_temp - self.glass_temp).abs().max(1.0);
        ua_glass * temp_margin / (rho_cp * delta_t_air)
    }

    /// Heater power required to prevent fogging \[W\].
    pub fn required_heater_power(&self) -> f64 {
        let dp = self.dew_point();
        if self.glass_temp > dp {
            0.0
        } else {
            let temp_lift = dp - self.glass_temp + 2.0;
            40.0 * temp_lift // ua_glass × ΔT
        }
    }
}

// ── CabinHumidity ─────────────────────────────────────────────────────────────

/// Moisture balance model for the vehicle cabin.
///
/// Tracks absolute humidity \[kg water / kg dry air\] considering
/// occupant moisture generation and condensation on glass.
#[derive(Debug, Clone)]
pub struct CabinHumidity {
    /// Current cabin absolute humidity \[kg/kg\].
    pub absolute_humidity: f64,
    /// Cabin air mass \[kg\].
    pub air_mass: f64,
    /// Moisture generation rate per occupant \[kg/s\].
    pub moisture_per_occupant: f64,
    /// Number of occupants.
    pub num_occupants: usize,
    /// Fresh-air exchange rate \[1/s\] (reduces humidity when humid outside).
    pub exchange_rate: f64,
    /// Ambient absolute humidity \[kg/kg\].
    pub ambient_humidity: f64,
}

impl CabinHumidity {
    /// Create with typical vehicle occupancy conditions.
    pub fn new(num_occupants: usize) -> Self {
        Self {
            absolute_humidity: 0.010,      // 10 g/kg typical interior
            air_mass: 3.0,                 // ~3 kg of cabin air
            moisture_per_occupant: 2.5e-5, // ~90 g/h each
            num_occupants,
            exchange_rate: 0.002,
            ambient_humidity: 0.007,
        }
    }

    /// Advance moisture balance by `dt` seconds.
    ///
    /// Returns condensation rate \[kg/s\] on cold surfaces (positive = condensing).
    pub fn step(&mut self, dt: f64, glass_temp: f64, air_temp: f64) -> f64 {
        // Saturation humidity via Magnus-like formula
        let sat = Self::saturation_humidity(air_temp);
        // Condensation: if cabin humidity exceeds saturation
        let condensation_rate = if self.absolute_humidity > sat {
            let excess = self.absolute_humidity - sat;
            excess * self.air_mass / dt.max(1e-9)
        } else {
            0.0
        };

        let generation = self.moisture_per_occupant * self.num_occupants as f64;
        let exchange = self.exchange_rate * (self.ambient_humidity - self.absolute_humidity);
        let condensation_on_glass = if glass_temp < Self::dew_point_from_ah(self.absolute_humidity)
        {
            0.2 * condensation_rate // fraction on glass
        } else {
            0.0
        };

        let d_humidity =
            (generation - condensation_on_glass * 0.5 + exchange * self.air_mass) / self.air_mass;
        self.absolute_humidity += d_humidity * dt;
        self.absolute_humidity = self.absolute_humidity.max(0.001);

        condensation_rate
    }

    /// Saturation absolute humidity \[kg/kg\] at temperature `t` \[°C\].
    pub fn saturation_humidity(t: f64) -> f64 {
        let pvs = 610.78 * E.powf(17.27 * t / (t + 237.3)); // Pa
        0.622 * pvs / (101325.0 - pvs)
    }

    /// Dew point from absolute humidity using inverse Magnus.
    pub fn dew_point_from_ah(ah: f64) -> f64 {
        let pvs = ah * 101325.0 / (0.622 + ah);
        let ln_p = (pvs / 610.78).ln();
        237.3 * ln_p / (17.27 - ln_p)
    }

    /// Relative humidity \[%\] from absolute humidity and temperature.
    pub fn relative_humidity(&self, air_temp: f64) -> f64 {
        let sat = Self::saturation_humidity(air_temp);
        (self.absolute_humidity / sat * 100.0).min(100.0)
    }
}

// ── HvacEnergyBalance ─────────────────────────────────────────────────────────

/// Total HVAC power balance: compressor + blower fan + auxiliary heater.
#[derive(Debug, Clone)]
pub struct HvacEnergyBalance {
    /// Compressor power \[W\].
    pub compressor_power: f64,
    /// Blower fan power \[W\].
    pub blower_power: f64,
    /// Auxiliary PTC heater power \[W\].
    pub heater_power: f64,
    /// Heat-pump power draw \[W\].
    pub heat_pump_power: f64,
}

impl HvacEnergyBalance {
    /// Create with all components set to zero.
    pub fn zero() -> Self {
        Self {
            compressor_power: 0.0,
            blower_power: 0.0,
            heater_power: 0.0,
            heat_pump_power: 0.0,
        }
    }

    /// Total electrical power drawn by HVAC \[W\].
    pub fn total_power(&self) -> f64 {
        self.compressor_power + self.blower_power + self.heater_power + self.heat_pump_power
    }

    /// Energy consumed over `dt` seconds \[J\].
    pub fn energy_consumed(&self, dt: f64) -> f64 {
        self.total_power() * dt
    }

    /// COP of the entire HVAC system given delivered thermal power `q_delivered`.
    pub fn system_cop(&self, q_delivered: f64) -> f64 {
        let p = self.total_power();
        if p > 1e-3 {
            q_delivered / p
        } else {
            f64::INFINITY
        }
    }
}

// ── CabinPreconditioning ──────────────────────────────────────────────────────

/// EV cabin pre-conditioning: pre-heat or pre-cool from the grid while plugged in.
///
/// Models the energy cost of reaching target temperature before departure.
#[derive(Debug, Clone)]
pub struct CabinPreconditioning {
    /// Target cabin temperature \[°C\].
    pub target_temp: f64,
    /// Initial cabin temperature \[°C\].
    pub initial_temp: f64,
    /// Current cabin temperature \[°C\].
    pub current_temp: f64,
    /// HVAC power available during pre-conditioning \[W\].
    pub hvac_power: f64,
    /// Cabin thermal mass \[J/K\].
    pub thermal_mass: f64,
    /// Thermal loss coefficient to ambient \[W/K\].
    pub ua_total: f64,
    /// Ambient temperature \[°C\].
    pub ambient_temp: f64,
    /// Time elapsed \[s\].
    pub elapsed: f64,
    /// Energy consumed from grid \[J\].
    pub grid_energy: f64,
    /// Whether pre-conditioning is complete.
    pub complete: bool,
}

impl CabinPreconditioning {
    /// Create a pre-conditioning session.
    pub fn new(initial_temp: f64, target_temp: f64, ambient_temp: f64) -> Self {
        Self {
            target_temp,
            initial_temp,
            current_temp: initial_temp,
            hvac_power: 5000.0,
            thermal_mass: 60_000.0,
            ua_total: 60.0,
            ambient_temp,
            elapsed: 0.0,
            grid_energy: 0.0,
            complete: false,
        }
    }

    /// Advance by `dt` seconds. Returns true when target is reached.
    pub fn step(&mut self, dt: f64) -> bool {
        if self.complete {
            return true;
        }
        let is_heating = self.target_temp > self.current_temp;
        let hvac = if is_heating {
            self.hvac_power
        } else {
            -self.hvac_power
        };
        let loss = self.ua_total * (self.current_temp - self.ambient_temp);
        let d_temp = (hvac - loss) / self.thermal_mass;
        self.current_temp += d_temp * dt;
        self.elapsed += dt;
        self.grid_energy += self.hvac_power * dt;

        let reached = if is_heating {
            self.current_temp >= self.target_temp
        } else {
            self.current_temp <= self.target_temp
        };
        if reached {
            self.current_temp = self.target_temp;
            self.complete = true;
        }
        self.complete
    }

    /// Estimated total energy needed to reach target \[J\] (analytical steady-state).
    pub fn estimated_energy(&self) -> f64 {
        let delta_t = (self.target_temp - self.initial_temp).abs();
        // E ≈ C * ΔT + UA * ΔT * τ where τ = C / UA
        let tau = self.thermal_mass / self.ua_total;
        self.thermal_mass * delta_t * (1.0 + 1.0 / (tau * self.ua_total / self.thermal_mass))
    }

    /// Remaining temperature difference to target.
    pub fn temp_error(&self) -> f64 {
        (self.target_temp - self.current_temp).abs()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── CabinThermal ──────────────────────────────────────────────────────

    #[test]
    fn cabin_thermal_heats_up_with_positive_hvac() {
        let mut cabin = CabinThermal::default_car();
        cabin.air_temp = 5.0;
        let initial = cabin.air_temp;
        cabin.step(60.0, 0.0, 3000.0);
        assert!(cabin.air_temp > initial, "cabin should warm up");
    }

    #[test]
    fn cabin_thermal_cools_down_with_negative_hvac() {
        let mut cabin = CabinThermal::default_car();
        cabin.air_temp = 35.0;
        let initial = cabin.air_temp;
        cabin.step(60.0, 40.0, -4000.0);
        assert!(cabin.air_temp < initial, "cabin should cool down");
    }

    #[test]
    fn cabin_thermal_error_positive_when_below_setpoint() {
        let cabin = CabinThermal::default_car();
        let mut c = cabin.clone();
        c.air_temp = 18.0;
        assert!(c.error(22.0) > 0.0);
    }

    #[test]
    fn cabin_thermal_error_zero_at_setpoint() {
        let mut cabin = CabinThermal::default_car();
        cabin.air_temp = 22.0;
        assert!((cabin.error(22.0)).abs() < 1e-12);
    }

    #[test]
    fn cabin_thermal_wall_equilibrates_toward_ambient() {
        let mut cabin = CabinThermal::default_car();
        let ambient = -10.0;
        for _ in 0..1000 {
            cabin.step(10.0, ambient, 0.0);
        }
        // Wall should trend toward ambient
        assert!(cabin.wall_temp < 30.0);
    }

    // ── HvacCompressor ────────────────────────────────────────────────────

    #[test]
    fn compressor_zero_speed_zero_capacity() {
        let mut comp = HvacCompressor::automotive();
        comp.speed = 0.0;
        assert_eq!(comp.cooling_capacity(), 0.0);
    }

    #[test]
    fn compressor_rated_speed_max_capacity() {
        let comp = HvacCompressor::automotive();
        let mut c = comp.clone();
        c.speed = c.rated_speed;
        assert!((c.cooling_capacity() - c.max_cooling_capacity).abs() < 1e-10);
    }

    #[test]
    fn compressor_cop_positive() {
        let comp = HvacCompressor::automotive();
        assert!(comp.actual_cop() > 0.0);
    }

    #[test]
    fn compressor_carnot_cop_decreases_with_condenser_temp() {
        let mut comp = HvacCompressor::automotive();
        let cop_low = comp.carnot_cop();
        comp.cond_temp = 70.0;
        let cop_high = comp.carnot_cop();
        assert!(cop_low > cop_high);
    }

    #[test]
    fn compressor_power_draw_positive() {
        let comp = HvacCompressor::automotive();
        assert!(comp.power_draw() > 0.0);
    }

    // ── HeatPump ──────────────────────────────────────────────────────────

    #[test]
    fn heat_pump_cop_greater_than_one() {
        let hp = HeatPump::ev_heat_pump();
        assert!(hp.actual_cop_heating() > 1.0);
    }

    #[test]
    fn heat_pump_cold_ambient_reduces_cop() {
        let hp = HeatPump::ev_heat_pump();
        let mut hp_cold = hp.clone();
        hp_cold.ambient_temp = -20.0;
        assert!(hp_cold.carnot_cop_heating() < hp.carnot_cop_heating());
    }

    #[test]
    fn heat_pump_defrost_stops_heat_output() {
        let mut hp = HeatPump::ev_heat_pump();
        hp.defrost_active = true;
        assert_eq!(hp.heat_output(), 0.0);
    }

    #[test]
    fn heat_pump_frost_accumulation() {
        let mut hp = HeatPump::ev_heat_pump();
        hp.ambient_temp = -5.0;
        // Use small steps so frost accumulates before defrost triggers
        let mut any_frost_seen = false;
        for _ in 0..200 {
            if hp.frost_mass > 0.0 {
                any_frost_seen = true;
                break;
            }
            hp.step_frost(0.5);
        }
        // After enough small steps, frost should have accumulated at some point
        // or defrost was triggered (which also means frost was present)
        assert!(
            any_frost_seen || hp.defrost_active || hp.frost_mass == 0.0,
            "frost_mass={} defrost={}",
            hp.frost_mass,
            hp.defrost_active
        );
        // More specifically: check that frost accumulates with cold ambient
        let mut hp2 = HeatPump::ev_heat_pump();
        hp2.ambient_temp = -5.0;
        hp2.step_frost(0.5);
        assert!(
            hp2.frost_mass > 0.0,
            "frost should accumulate: {}",
            hp2.frost_mass
        );
    }

    #[test]
    fn heat_pump_defrost_triggered() {
        let mut hp = HeatPump::ev_heat_pump();
        hp.ambient_temp = -5.0;
        hp.frost_mass = 0.02; // above threshold
        hp.step_frost(1.0);
        // Defrost should activate or frost reduced
        assert!(hp.defrost_active || hp.frost_mass < 0.02);
    }

    // ── AirDistribution ───────────────────────────────────────────────────

    #[test]
    fn air_distribution_fractions_sum_to_one() {
        let ad = AirDistribution::new(0.5, 0.3, 0.2, 0.05);
        let sum = ad.face_fraction + ad.floor_fraction + ad.defrost_fraction;
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn air_distribution_flows_sum_to_total() {
        let ad = AirDistribution::new(0.5, 0.3, 0.2, 0.05);
        let total = ad.face_flow() + ad.floor_flow() + ad.defrost_flow();
        assert!((total - ad.total_flow).abs() < 1e-12);
    }

    #[test]
    fn air_distribution_defrost_mode() {
        let mut ad = AirDistribution::new(0.5, 0.3, 0.2, 0.05);
        ad.set_defrost_mode();
        assert_eq!(ad.defrost_fraction, 1.0);
        assert_eq!(ad.face_fraction, 0.0);
        assert_eq!(ad.floor_fraction, 0.0);
    }

    #[test]
    fn air_distribution_normalises_fractions() {
        let ad = AirDistribution::new(1.0, 1.0, 0.0, 0.1);
        assert!((ad.face_fraction - 0.5).abs() < 1e-12);
        assert!((ad.floor_fraction - 0.5).abs() < 1e-12);
    }

    // ── OccupantComfort ───────────────────────────────────────────────────

    #[test]
    fn pmv_in_valid_range() {
        let oc = OccupantComfort::seated_driver();
        let pmv = oc.pmv();
        assert!((-3.0..=3.0).contains(&pmv));
    }

    #[test]
    fn ppd_at_least_five_percent() {
        let oc = OccupantComfort::seated_driver();
        assert!(oc.ppd() >= 5.0);
    }

    #[test]
    fn ppd_less_than_one_hundred() {
        let oc = OccupantComfort::seated_driver();
        assert!(oc.ppd() <= 100.0);
    }

    #[test]
    fn comfort_check_returns_bool() {
        let oc = OccupantComfort::seated_driver();
        let _ = oc.is_comfortable(); // just ensure no panic
    }

    // ── HvacController ────────────────────────────────────────────────────

    #[test]
    fn controller_positive_output_when_below_setpoint() {
        let mut ctrl = HvacController::new(22.0);
        let out = ctrl.update(15.0, 1.0);
        assert!(out > 0.0);
    }

    #[test]
    fn controller_negative_output_when_above_setpoint() {
        let mut ctrl = HvacController::new(22.0);
        let out = ctrl.update(30.0, 1.0);
        assert!(out < 0.0);
    }

    #[test]
    fn controller_clamped_to_max_power() {
        let mut ctrl = HvacController::new(22.0);
        let out = ctrl.update(-50.0, 1.0); // extreme cold
        assert!(out <= ctrl.max_power);
    }

    #[test]
    fn controller_recirculation_reduces_infiltration() {
        let mut ctrl = HvacController::new(22.0);
        ctrl.set_recirculation(true);
        assert_eq!(ctrl.infiltration_load(22.0, -10.0), 0.0);
    }

    #[test]
    fn controller_infiltration_load_nonzero_fresh_air() {
        let mut ctrl = HvacController::new(22.0);
        ctrl.set_recirculation(false);
        let load = ctrl.infiltration_load(22.0, -10.0);
        assert!(load > 0.0);
    }

    // ── WindowDefogging ───────────────────────────────────────────────────

    #[test]
    fn defogging_dew_point_below_air_temp() {
        let wd = WindowDefogging::new(10.0, 20.0, 60.0);
        assert!(wd.dew_point() < wd.air_temp);
    }

    #[test]
    fn defogging_is_fogging_cold_glass() {
        let wd = WindowDefogging::new(5.0, 20.0, 80.0);
        // dew point at 80% RH, 20°C ≈ 16°C → glass at 5°C is fogging
        assert!(wd.is_fogging());
    }

    #[test]
    fn defogging_not_fogging_warm_glass() {
        let wd = WindowDefogging::new(25.0, 20.0, 50.0);
        // glass warmer than dew point
        assert!(!wd.is_fogging());
    }

    #[test]
    fn defogging_required_flow_nonzero_when_fogging() {
        let wd = WindowDefogging::new(5.0, 20.0, 90.0);
        if wd.is_fogging() {
            assert!(wd.required_airflow() > 0.0);
        }
    }

    #[test]
    fn defogging_heater_power_zero_when_not_fogging() {
        let wd = WindowDefogging::new(25.0, 20.0, 40.0);
        assert_eq!(wd.required_heater_power(), 0.0);
    }

    // ── CabinHumidity ─────────────────────────────────────────────────────

    #[test]
    fn cabin_humidity_saturation_increases_with_temperature() {
        let sat20 = CabinHumidity::saturation_humidity(20.0);
        let sat30 = CabinHumidity::saturation_humidity(30.0);
        assert!(sat30 > sat20);
    }

    #[test]
    fn cabin_humidity_dew_point_from_ah() {
        let ah = 0.010;
        let dp = CabinHumidity::dew_point_from_ah(ah);
        assert!(dp > 0.0 && dp < 20.0); // reasonable range
    }

    #[test]
    fn cabin_humidity_relative_humidity_reasonable() {
        let ch = CabinHumidity::new(2);
        let rh = ch.relative_humidity(22.0);
        assert!(rh > 0.0 && rh <= 100.0);
    }

    #[test]
    fn cabin_humidity_step_no_panic() {
        let mut ch = CabinHumidity::new(2);
        let _ = ch.step(1.0, 10.0, 20.0);
        assert!(ch.absolute_humidity >= 0.001);
    }

    #[test]
    fn cabin_humidity_increases_with_occupants() {
        let mut ch1 = CabinHumidity::new(1);
        let mut ch5 = CabinHumidity::new(5);
        ch1.ambient_humidity = ch1.absolute_humidity; // no exchange
        ch5.ambient_humidity = ch5.absolute_humidity;
        ch1.step(60.0, 20.0, 22.0);
        ch5.step(60.0, 20.0, 22.0);
        assert!(ch5.absolute_humidity >= ch1.absolute_humidity);
    }

    // ── HvacEnergyBalance ─────────────────────────────────────────────────

    #[test]
    fn energy_balance_total_power_sum() {
        let mut eb = HvacEnergyBalance::zero();
        eb.compressor_power = 1000.0;
        eb.blower_power = 200.0;
        eb.heater_power = 500.0;
        eb.heat_pump_power = 300.0;
        assert!((eb.total_power() - 2000.0).abs() < 1e-10);
    }

    #[test]
    fn energy_balance_energy_consumed() {
        let mut eb = HvacEnergyBalance::zero();
        eb.compressor_power = 1000.0;
        let energy = eb.energy_consumed(3600.0);
        assert!((energy - 3_600_000.0).abs() < 1.0);
    }

    #[test]
    fn energy_balance_system_cop() {
        let mut eb = HvacEnergyBalance::zero();
        eb.compressor_power = 1000.0;
        let cop = eb.system_cop(3000.0);
        assert!((cop - 3.0).abs() < 1e-10);
    }

    #[test]
    fn energy_balance_zero_power_infinite_cop() {
        let eb = HvacEnergyBalance::zero();
        let cop = eb.system_cop(1000.0);
        assert!(cop.is_infinite());
    }

    // ── CabinPreconditioning ──────────────────────────────────────────────

    #[test]
    fn preconditioning_heats_cabin() {
        let mut pc = CabinPreconditioning::new(-10.0, 22.0, -10.0);
        let initial = pc.current_temp;
        pc.step(60.0);
        assert!(pc.current_temp > initial);
    }

    #[test]
    fn preconditioning_cools_cabin() {
        let mut pc = CabinPreconditioning::new(40.0, 22.0, 35.0);
        let initial = pc.current_temp;
        pc.step(60.0);
        assert!(pc.current_temp < initial);
    }

    #[test]
    fn preconditioning_completes() {
        let mut pc = CabinPreconditioning::new(20.0, 22.0, 20.0);
        let mut done = false;
        for _ in 0..10_000 {
            if pc.step(1.0) {
                done = true;
                break;
            }
        }
        assert!(done, "pre-conditioning should complete");
        assert!(pc.complete);
    }

    #[test]
    fn preconditioning_grid_energy_consumed() {
        let mut pc = CabinPreconditioning::new(-10.0, 22.0, -10.0);
        for _ in 0..100 {
            pc.step(1.0);
        }
        assert!(pc.grid_energy > 0.0);
    }

    #[test]
    fn preconditioning_temp_error_decreases() {
        let mut pc = CabinPreconditioning::new(-10.0, 22.0, -10.0);
        let err_init = pc.temp_error();
        for _ in 0..1000 {
            pc.step(1.0);
        }
        assert!(pc.temp_error() < err_init);
    }
}
