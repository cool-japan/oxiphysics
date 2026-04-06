// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Brake system simulation: hydraulic circuits, ABS, thermal models, and regenerative braking.
//!
//! Provides comprehensive brake system physics including:
//! - Hydraulic pressure distribution across axles
//! - Caliper piston and pad friction models
//! - Thermal fade and rotor stress analysis
//! - ABS wheel-slip regulation
//! - Regenerative braking integration
//! - Brake-by-wire actuation
//! - Pedal feel simulation

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Atmospheric pressure in Pa.
const ATMOSPHERIC_PRESSURE: f64 = 101_325.0;

/// Steel specific heat capacity in J/(kg·K).
const STEEL_SPECIFIC_HEAT: f64 = 490.0;

/// Cast iron specific heat capacity in J/(kg·K).
const CAST_IRON_SPECIFIC_HEAT: f64 = 460.0;

/// Minimum brake line pressure considered active (Pa).
const MIN_ACTIVE_PRESSURE: f64 = 50_000.0;

/// Standard brake fluid density in kg/m³.
const BRAKE_FLUID_DENSITY: f64 = 1_050.0;

/// Speed of sound in brake fluid (bulk modulus approximation) in m/s.
const BRAKE_FLUID_BULK_MODULUS: f64 = 1.8e9;

/// Gravity acceleration in m/s².
const G: f64 = 9.81;

// ---------------------------------------------------------------------------
// Brake Pad Material
// ---------------------------------------------------------------------------

/// Brake pad friction material properties.
#[derive(Debug, Clone, Copy)]
pub struct PadMaterial {
    /// Cold friction coefficient (at ambient temperature).
    pub friction_cold: f64,
    /// Hot friction coefficient (at peak operating temperature).
    pub friction_hot: f64,
    /// Temperature at which hot coefficient applies (°C).
    pub peak_temp: f64,
    /// Fade onset temperature (°C) — friction starts dropping above this.
    pub fade_onset_temp: f64,
    /// Minimum friction coefficient at extreme heat.
    pub friction_fade_floor: f64,
    /// Pad thermal conductivity W/(m·K).
    pub conductivity: f64,
    /// Pad density kg/m³.
    pub density: f64,
    /// Pad specific heat J/(kg·K).
    pub specific_heat: f64,
    /// Wear coefficient (m³/(N·m)).
    pub wear_coefficient: f64,
}

impl PadMaterial {
    /// Standard organic/non-asbestos (NAO) compound.
    pub fn organic() -> Self {
        Self {
            friction_cold: 0.38,
            friction_hot: 0.42,
            peak_temp: 300.0,
            fade_onset_temp: 400.0,
            friction_fade_floor: 0.20,
            conductivity: 1.0,
            density: 2_200.0,
            specific_heat: 900.0,
            wear_coefficient: 3e-14,
        }
    }

    /// Semi-metallic compound — higher heat tolerance.
    pub fn semi_metallic() -> Self {
        Self {
            friction_cold: 0.35,
            friction_hot: 0.44,
            peak_temp: 450.0,
            fade_onset_temp: 600.0,
            friction_fade_floor: 0.22,
            conductivity: 4.0,
            density: 2_800.0,
            specific_heat: 700.0,
            wear_coefficient: 2e-14,
        }
    }

    /// Carbon-ceramic compound for high-performance applications.
    pub fn carbon_ceramic() -> Self {
        Self {
            friction_cold: 0.32,
            friction_hot: 0.48,
            peak_temp: 700.0,
            fade_onset_temp: 900.0,
            friction_fade_floor: 0.28,
            conductivity: 10.0,
            density: 2_000.0,
            specific_heat: 800.0,
            wear_coefficient: 5e-15,
        }
    }

    /// Compute friction coefficient as a function of temperature.
    pub fn friction_at_temp(&self, temp_c: f64) -> f64 {
        if temp_c <= self.peak_temp {
            // Linear ramp from cold to hot
            let t = (temp_c / self.peak_temp).clamp(0.0, 1.0);
            self.friction_cold + t * (self.friction_hot - self.friction_cold)
        } else if temp_c <= self.fade_onset_temp {
            self.friction_hot
        } else {
            // Thermal fade — exponential decay
            let excess = temp_c - self.fade_onset_temp;
            let range = self.friction_hot - self.friction_fade_floor;
            let decay = (-excess / 300.0).exp();
            self.friction_fade_floor + range * decay
        }
    }
}

// ---------------------------------------------------------------------------
// Brake Rotor
// ---------------------------------------------------------------------------

/// Brake rotor (disc) physical model.
#[derive(Debug, Clone)]
pub struct BrakeRotor {
    /// Rotor outer radius in m.
    pub outer_radius: f64,
    /// Rotor inner radius in m.
    pub inner_radius: f64,
    /// Rotor thickness in m.
    pub thickness: f64,
    /// Rotor mass in kg.
    pub mass: f64,
    /// Rotor material specific heat J/(kg·K).
    pub specific_heat: f64,
    /// Rotor thermal conductivity W/(m·K).
    pub conductivity: f64,
    /// Emissivity for radiation cooling (0–1).
    pub emissivity: f64,
    /// Current rotor temperature in °C.
    pub temperature: f64,
    /// Convective heat transfer coefficient W/(m²·K).
    pub h_conv: f64,
    /// Total accumulated wear depth in m.
    pub wear_depth: f64,
    /// Vented disc flag.
    pub vented: bool,
}

impl BrakeRotor {
    /// Create a new rotor with default cast-iron properties.
    pub fn new(outer_r: f64, inner_r: f64, thickness: f64, mass: f64) -> Self {
        Self {
            outer_radius: outer_r,
            inner_radius: inner_r,
            thickness,
            mass,
            specific_heat: CAST_IRON_SPECIFIC_HEAT,
            conductivity: 50.0,
            emissivity: 0.85,
            temperature: 25.0,
            h_conv: 60.0,
            wear_depth: 0.0,
            vented: false,
        }
    }

    /// Create a vented rotor (higher convective coefficient).
    pub fn new_vented(outer_r: f64, inner_r: f64, thickness: f64, mass: f64) -> Self {
        let mut r = Self::new(outer_r, inner_r, thickness, mass);
        r.h_conv = 120.0;
        r.vented = true;
        r
    }

    /// Effective braking radius (mean of inner and outer contact patch).
    pub fn effective_radius(&self) -> f64 {
        (self.outer_radius + self.inner_radius) * 0.5
    }

    /// Rotor swept area in m² (both sides).
    pub fn swept_area(&self) -> f64 {
        PI * (self.outer_radius * self.outer_radius - self.inner_radius * self.inner_radius)
    }

    /// Total rotor surface area for cooling (top + bottom + edge).
    pub fn cooling_area(&self) -> f64 {
        let face = 2.0 * self.swept_area();
        let edge = 2.0 * PI * (self.outer_radius + self.inner_radius) * self.thickness;
        face + edge
    }

    /// Update rotor temperature given heat input and time step.
    ///
    /// Uses lumped capacitance model with convection and radiation.
    pub fn update_temperature(&mut self, heat_in_watts: f64, ambient_temp: f64, dt: f64) {
        let thermal_mass = self.mass * self.specific_heat; // J/K
        let temp_k = self.temperature + 273.15;
        let ambient_k = ambient_temp + 273.15;

        // Convective cooling
        let q_conv = self.h_conv * self.cooling_area() * (self.temperature - ambient_temp);

        // Radiative cooling (Stefan-Boltzmann)
        const SIGMA: f64 = 5.67e-8;
        let q_rad =
            self.emissivity * SIGMA * self.cooling_area() * (temp_k.powi(4) - ambient_k.powi(4));

        let net_heat = heat_in_watts - q_conv - q_rad;
        self.temperature += net_heat * dt / thermal_mass;
        self.temperature = self.temperature.max(ambient_temp);
    }

    /// Compute hoop stress in rotor due to thermal gradient (simplified).
    ///
    /// Returns stress in Pa.
    pub fn thermal_hoop_stress(&self, delta_t: f64) -> f64 {
        // σ_θ ≈ E·α·ΔT / (1 - ν) for a disc with uniform ΔT
        let e = 165e9_f64; // Young's modulus cast iron, Pa
        let alpha = 11e-6_f64; // thermal expansion coefficient, 1/K
        let nu = 0.26_f64; // Poisson's ratio
        e * alpha * delta_t / (1.0 - nu)
    }

    /// Effective braking torque per unit clamping force (N·m per N).
    pub fn torque_per_clamp(&self) -> f64 {
        self.effective_radius()
    }

    /// Accumulate wear on the rotor given a wear rate.
    pub fn apply_wear(&mut self, energy_j: f64, wear_coeff: f64) {
        // Archard wear: V = k · E / hardness (simplified)
        let hardness = 1_500e6_f64; // Pa, approximate for cast iron
        let volume_worn = wear_coeff * energy_j / hardness;
        let area = self.swept_area();
        if area > 1e-10 {
            self.wear_depth += volume_worn / area;
        }
    }
}

// ---------------------------------------------------------------------------
// Brake Caliper
// ---------------------------------------------------------------------------

/// Hydraulic brake caliper piston model.
#[derive(Debug, Clone)]
pub struct BrakeCaliper {
    /// Number of pistons.
    pub num_pistons: u32,
    /// Piston diameter in m.
    pub piston_diameter: f64,
    /// Seal stiffness (spring-back) in N/m.
    pub seal_stiffness: f64,
    /// Piston displacement at zero pressure in m (running clearance).
    pub running_clearance: f64,
    /// Brake pad material.
    pub pad_material: PadMaterial,
    /// Current pad temperature in °C.
    pub pad_temperature: f64,
    /// Remaining pad thickness in m.
    pub pad_thickness: f64,
    /// Initial pad thickness in m.
    pub initial_pad_thickness: f64,
}

impl BrakeCaliper {
    /// Create a new caliper with given number of pistons and piston diameter.
    pub fn new(num_pistons: u32, piston_diameter: f64, pad_material: PadMaterial) -> Self {
        Self {
            num_pistons,
            piston_diameter,
            seal_stiffness: 50_000.0,
            running_clearance: 0.0002,
            pad_material,
            pad_temperature: 25.0,
            pad_thickness: 0.012,
            initial_pad_thickness: 0.012,
        }
    }

    /// Total piston area in m².
    pub fn piston_area(&self) -> f64 {
        let single = PI * 0.25 * self.piston_diameter * self.piston_diameter;
        single * self.num_pistons as f64
    }

    /// Clamping force given hydraulic pressure in Pa.
    pub fn clamping_force(&self, pressure: f64) -> f64 {
        let f_hydraulic = self.piston_area() * pressure;
        let f_seal = self.seal_stiffness * self.running_clearance;
        (f_hydraulic - f_seal).max(0.0)
    }

    /// Friction force on rotor face (one side) given hydraulic pressure and rotor temp.
    pub fn friction_force(&self, pressure: f64, rotor_temp: f64) -> f64 {
        let mu = self.pad_material.friction_at_temp(rotor_temp);
        self.clamping_force(pressure) * mu
    }

    /// Braking torque about wheel axis given pressure, rotor, and effective radius.
    pub fn braking_torque(&self, pressure: f64, rotor: &BrakeRotor) -> f64 {
        let ff = self.friction_force(pressure, rotor.temperature);
        ff * rotor.effective_radius() * 2.0 // two friction surfaces (disc both sides)
    }

    /// Heat generated at pad-rotor interface per second in W.
    pub fn heat_generation_rate(
        &self,
        pressure: f64,
        rotor: &BrakeRotor,
        wheel_speed_rad_s: f64,
    ) -> f64 {
        let torque = self.braking_torque(pressure, rotor);
        torque * wheel_speed_rad_s.abs()
    }

    /// Update pad temperature (lumped model).
    pub fn update_pad_temperature(&mut self, heat_w: f64, ambient: f64, dt: f64) {
        let pad_mass = self.pad_material.density
            * PI
            * 0.25
            * self.piston_diameter
            * self.piston_diameter
            * self.pad_thickness;
        let thermal_mass = pad_mass * self.pad_material.specific_heat;
        let q_conv = 30.0
            * PI
            * self.piston_diameter
            * self.pad_thickness
            * (self.pad_temperature - ambient);
        self.pad_temperature += (heat_w - q_conv) * dt / thermal_mass.max(1e-6);
        self.pad_temperature = self.pad_temperature.max(ambient);
    }

    /// Fraction of pad remaining (1.0 = new, 0.0 = worn out).
    pub fn pad_life_fraction(&self) -> f64 {
        (self.pad_thickness / self.initial_pad_thickness).clamp(0.0, 1.0)
    }

    /// Apply pad wear given contact energy.
    pub fn apply_pad_wear(&mut self, contact_force: f64, sliding_distance: f64) {
        let wear_vol = self.pad_material.wear_coefficient * contact_force * sliding_distance;
        let area = self.piston_area();
        if area > 1e-10 {
            self.pad_thickness = (self.pad_thickness - wear_vol / area).max(0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Hydraulic Brake Circuit
// ---------------------------------------------------------------------------

/// Hydraulic circuit model for a single brake circuit (front or rear).
#[derive(Debug, Clone)]
pub struct HydraulicCircuit {
    /// Master cylinder bore diameter in m.
    pub master_cylinder_bore: f64,
    /// Total line volume in m³.
    pub line_volume: f64,
    /// Brake fluid bulk modulus in Pa.
    pub fluid_bulk_modulus: f64,
    /// Line compliance (volume change per unit pressure) in m³/Pa.
    pub compliance: f64,
    /// Current line pressure in Pa.
    pub pressure: f64,
    /// Vapor-lock temperature threshold in °C (boiling point of fluid).
    pub vapor_lock_temp: f64,
    /// Current fluid temperature in °C.
    pub fluid_temperature: f64,
}

impl HydraulicCircuit {
    /// Create a new hydraulic circuit.
    pub fn new(master_bore: f64, line_volume_m3: f64) -> Self {
        Self {
            master_cylinder_bore: master_bore,
            line_volume: line_volume_m3,
            fluid_bulk_modulus: BRAKE_FLUID_BULK_MODULUS,
            compliance: line_volume_m3 / BRAKE_FLUID_BULK_MODULUS,
            pressure: ATMOSPHERIC_PRESSURE,
            vapor_lock_temp: 230.0,
            fluid_temperature: 25.0,
        }
    }

    /// Update line pressure given piston displacement input (m).
    ///
    /// Returns new gauge pressure in Pa.
    pub fn update_pressure(&mut self, piston_displacement: f64) -> f64 {
        let mc_area = PI * 0.25 * self.master_cylinder_bore * self.master_cylinder_bore;
        let delta_vol = mc_area * piston_displacement;
        let delta_pressure = delta_vol / self.compliance;
        self.pressure = (ATMOSPHERIC_PRESSURE + delta_pressure).max(ATMOSPHERIC_PRESSURE);
        self.pressure
    }

    /// Compute master cylinder displacement required for target pressure.
    pub fn displacement_for_pressure(&self, target_gauge_pressure: f64) -> f64 {
        let mc_area = PI * 0.25 * self.master_cylinder_bore * self.master_cylinder_bore;
        let delta_vol = target_gauge_pressure * self.compliance;
        delta_vol / mc_area
    }

    /// Check if vapor lock condition (fluid boiling) is occurring.
    pub fn is_vapor_locked(&self) -> bool {
        self.fluid_temperature >= self.vapor_lock_temp
    }

    /// Effective pressure considering vapor lock (zero if boiling).
    pub fn effective_pressure(&self) -> f64 {
        if self.is_vapor_locked() {
            ATMOSPHERIC_PRESSURE
        } else {
            self.pressure
        }
    }

    /// Gauge pressure (above atmospheric).
    pub fn gauge_pressure(&self) -> f64 {
        (self.pressure - ATMOSPHERIC_PRESSURE).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// Brake Bias Controller
// ---------------------------------------------------------------------------

/// Front/rear brake bias control.
#[derive(Debug, Clone)]
pub struct BrakeBiasController {
    /// Bias ratio: fraction of total pressure sent to front (0.0–1.0).
    pub front_bias: f64,
    /// Whether the bias is adjustable (e.g., brake balance bar).
    pub adjustable: bool,
    /// Minimum allowed front bias.
    pub min_bias: f64,
    /// Maximum allowed front bias.
    pub max_bias: f64,
}

impl BrakeBiasController {
    /// Create a new bias controller with given initial front bias fraction.
    pub fn new(front_bias: f64) -> Self {
        Self {
            front_bias: front_bias.clamp(0.3, 0.85),
            adjustable: true,
            min_bias: 0.3,
            max_bias: 0.85,
        }
    }

    /// Front brake pressure fraction.
    pub fn front_fraction(&self) -> f64 {
        self.front_bias
    }

    /// Rear brake pressure fraction.
    pub fn rear_fraction(&self) -> f64 {
        1.0 - self.front_bias
    }

    /// Adjust bias by a delta (positive = more front).
    pub fn adjust(&mut self, delta: f64) {
        if self.adjustable {
            self.front_bias = (self.front_bias + delta).clamp(self.min_bias, self.max_bias);
        }
    }

    /// Compute front and rear pressures from total master cylinder pressure.
    pub fn split_pressure(&self, total_pressure: f64) -> (f64, f64) {
        let front = total_pressure * self.front_bias;
        let rear = total_pressure * (1.0 - self.front_bias);
        (front, rear)
    }

    /// Optimal bias for a given deceleration and weight transfer.
    ///
    /// `decel_g` – deceleration in units of g.
    /// `front_weight_static` – static front axle weight fraction.
    /// `h_cg` – centre of gravity height in m.
    /// `wheelbase` – vehicle wheelbase in m.
    pub fn optimal_bias(decel_g: f64, front_weight_static: f64, h_cg: f64, wheelbase: f64) -> f64 {
        // Weight transfer fraction to front during braking
        let wt_transfer = decel_g * h_cg / wheelbase;
        let front_dynamic = front_weight_static + wt_transfer;
        front_dynamic.clamp(0.30, 0.85)
    }
}

// ---------------------------------------------------------------------------
// ABS Controller
// ---------------------------------------------------------------------------

/// ABS (Anti-lock Braking System) wheel-slip regulation controller.
#[derive(Debug, Clone)]
pub struct AbsBrakeController {
    /// Target wheel slip ratio (typical 0.1–0.2 for peak friction).
    pub target_slip: f64,
    /// Slip tolerance band (pressure release starts above target + band).
    pub slip_tolerance: f64,
    /// Current ABS active flag per wheel (4-wheel assumed).
    pub active: [bool; 4],
    /// Modulation pressure reduction rate in Pa/s.
    pub release_rate: f64,
    /// Modulation pressure build rate in Pa/s.
    pub build_rate: f64,
    /// Minimum cycle hold time in seconds.
    pub hold_time: f64,
    /// Time since last transition for each wheel.
    pub cycle_timer: [f64; 4],
}

impl AbsBrakeController {
    /// Create a default ABS controller.
    pub fn new() -> Self {
        Self {
            target_slip: 0.15,
            slip_tolerance: 0.05,
            active: [false; 4],
            release_rate: 4_000_000.0, // 4 MPa/s
            build_rate: 2_000_000.0,   // 2 MPa/s
            hold_time: 0.02,           // 20 ms hold
            cycle_timer: [0.0; 4],
        }
    }

    /// Compute wheel slip ratio from vehicle speed and wheel speed.
    ///
    /// `v_vehicle` – vehicle speed at wheel contact (m/s).
    /// `v_wheel` – wheel peripheral speed (m/s).
    pub fn wheel_slip(v_vehicle: f64, v_wheel: f64) -> f64 {
        if v_vehicle.abs() < 0.1 {
            return 0.0;
        }
        (v_vehicle - v_wheel) / v_vehicle.abs()
    }

    /// Regulate brake pressure for one wheel.
    ///
    /// Returns pressure modifier in Pa (negative = reduce, positive = increase).
    pub fn regulate_wheel(
        &mut self,
        wheel_idx: usize,
        current_pressure: f64,
        slip_ratio: f64,
        dt: f64,
    ) -> f64 {
        self.cycle_timer[wheel_idx] += dt;

        let over_threshold = slip_ratio > self.target_slip + self.slip_tolerance;
        let under_threshold = slip_ratio < self.target_slip - self.slip_tolerance;

        if over_threshold && self.cycle_timer[wheel_idx] >= self.hold_time {
            // Release phase
            self.active[wheel_idx] = true;
            self.cycle_timer[wheel_idx] = 0.0;
            -self.release_rate * dt
        } else if under_threshold && self.active[wheel_idx] {
            // Build phase
            let max_build = current_pressure * 0.1; // cap at 10% of current per step
            (self.build_rate * dt).min(max_build)
        } else {
            // Hold
            0.0
        }
    }

    /// Reset ABS state (e.g., when vehicle stops).
    pub fn reset(&mut self) {
        self.active = [false; 4];
        self.cycle_timer = [0.0; 4];
    }

    /// Is ABS currently intervening on any wheel?
    pub fn is_active(&self) -> bool {
        self.active.iter().any(|&a| a)
    }
}

impl Default for AbsBrakeController {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Regenerative Braking
// ---------------------------------------------------------------------------

/// Regenerative braking integration model.
#[derive(Debug, Clone)]
pub struct RegenerativeBraking {
    /// Maximum regenerative torque available per driven wheel (N·m).
    pub max_regen_torque: f64,
    /// Generator efficiency (0.0–1.0).
    pub efficiency: f64,
    /// Battery state of charge (0.0–1.0).
    pub soc: f64,
    /// Maximum SOC above which regen is disabled.
    pub max_soc: f64,
    /// Maximum recoverable power in W.
    pub max_power_w: f64,
    /// Blending threshold — if friction braking exceeds this, blend with regen.
    pub blend_threshold_torque: f64,
    /// Current regenerated energy in J.
    pub energy_recovered_j: f64,
}

impl RegenerativeBraking {
    /// Create a new regenerative braking model.
    pub fn new(max_torque: f64, efficiency: f64, max_power_w: f64) -> Self {
        Self {
            max_regen_torque: max_torque,
            efficiency,
            soc: 0.5,
            max_soc: 0.95,
            max_power_w,
            blend_threshold_torque: max_torque * 0.5,
            energy_recovered_j: 0.0,
        }
    }

    /// Available regenerative torque given current conditions.
    ///
    /// Reduces regen near max SOC and zero speed.
    pub fn available_torque(&self, wheel_speed_rad_s: f64) -> f64 {
        if self.soc >= self.max_soc {
            return 0.0;
        }
        if wheel_speed_rad_s.abs() < 0.5 {
            return 0.0;
        }
        // Power limit
        let power_limited = self.max_power_w / wheel_speed_rad_s.abs();
        // SOC derating
        let soc_factor = 1.0 - (self.soc / self.max_soc).powi(4);
        self.max_regen_torque.min(power_limited) * soc_factor
    }

    /// Apply regenerative braking, returns actual regen torque applied.
    pub fn apply(&mut self, requested_torque: f64, wheel_speed_rad_s: f64, dt: f64) -> f64 {
        let available = self.available_torque(wheel_speed_rad_s);
        let actual = requested_torque.min(available).max(0.0);
        let power = actual * wheel_speed_rad_s.abs();
        let energy = power * dt * self.efficiency;
        self.energy_recovered_j += energy;
        // Update SOC (assuming 100 Wh battery capacity for normalization)
        let battery_capacity = 3_600_000.0; // 1 kWh in J as default
        self.soc = (self.soc + energy / battery_capacity).min(self.max_soc);
        actual
    }

    /// Blend friction and regenerative braking to meet total demanded torque.
    ///
    /// Returns `(friction_torque, regen_torque)`.
    pub fn blend_torques(&self, demanded: f64, wheel_speed: f64) -> (f64, f64) {
        let max_regen = self.available_torque(wheel_speed);
        let regen = demanded.min(max_regen).max(0.0);
        let friction = (demanded - regen).max(0.0);
        (friction, regen)
    }
}

// ---------------------------------------------------------------------------
// Brake Pedal Simulator
// ---------------------------------------------------------------------------

/// Brake pedal feel / simulator model.
#[derive(Debug, Clone)]
pub struct BrakePedalSimulator {
    /// Pedal ratio (mechanical advantage).
    pub pedal_ratio: f64,
    /// Pedal travel range in m.
    pub max_travel: f64,
    /// Dead-band (free play) in m.
    pub free_play: f64,
    /// Non-linear pedal stiffness coefficients \[a₀, a₁, a₂\] (Pa/m polynomial).
    pub stiffness_poly: [f64; 3],
    /// Pedal return spring stiffness N/m.
    pub return_spring: f64,
    /// Pedal damping N·s/m.
    pub damping: f64,
    /// Current pedal position in m (0 = released).
    pub position: f64,
    /// Current pedal velocity in m/s.
    pub velocity: f64,
}

impl BrakePedalSimulator {
    /// Create a default pedal simulator.
    pub fn new() -> Self {
        Self {
            pedal_ratio: 4.0,
            max_travel: 0.15,
            free_play: 0.005,
            stiffness_poly: [0.0, 5_000_000.0, 15_000_000.0],
            return_spring: 80.0,
            damping: 5.0,
            position: 0.0,
            velocity: 0.0,
        }
    }

    /// Master cylinder force given pedal force applied by driver.
    pub fn mc_force(&self, driver_force: f64) -> f64 {
        driver_force * self.pedal_ratio
    }

    /// Pedal pressure output (Pa) given pedal travel in m.
    pub fn pressure_from_travel(&self, travel: f64) -> f64 {
        let t = (travel - self.free_play).max(0.0);
        let p =
            self.stiffness_poly[0] + self.stiffness_poly[1] * t + self.stiffness_poly[2] * t * t;
        p.max(0.0)
    }

    /// Driver-perceived pedal feedback force (N) for given travel.
    pub fn feedback_force(&self, travel: f64) -> f64 {
        let p = self.pressure_from_travel(travel);
        let mc_area = PI * 0.25 * 0.022_f64 * 0.022; // 22 mm MC bore
        p * mc_area / self.pedal_ratio
    }

    /// Simulate pedal dynamics given driver force and time step.
    ///
    /// Returns current pedal position (m).
    pub fn step(&mut self, driver_force_n: f64, dt: f64) -> f64 {
        let travel = self.position;
        let restore = self.pressure_from_travel(travel) * PI * 0.25 * 0.022 * 0.022
            + self.return_spring * travel
            + self.damping * self.velocity;
        let mc_force = self.mc_force(driver_force_n);
        let net = mc_force - restore;
        // Assume pedal + linkage effective mass of 1.5 kg
        let mass = 1.5_f64;
        let accel = net / mass;
        self.velocity += accel * dt;
        self.position = (self.position + self.velocity * dt).clamp(0.0, self.max_travel);
        if self.position <= 0.0 && self.velocity < 0.0 {
            self.velocity = 0.0;
        }
        self.position
    }
}

impl Default for BrakePedalSimulator {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Brake-by-Wire
// ---------------------------------------------------------------------------

/// Brake-by-wire (BBW) actuation model.
#[derive(Debug, Clone)]
pub struct BrakeByWire {
    /// Target pressure per wheel in Pa (set by control unit).
    pub target_pressure: [f64; 4],
    /// Actual pressure per wheel in Pa.
    pub actual_pressure: [f64; 4],
    /// Maximum pressure rate of change in Pa/s.
    pub pressure_rate: f64,
    /// Latency of the electronic actuator in seconds.
    pub latency: f64,
    /// Accumulated latency buffer for each wheel.
    latency_buffer: [f64; 4],
    /// Redundancy active flag (fallback to hydraulic).
    pub hydraulic_fallback: bool,
}

impl BrakeByWire {
    /// Create a new brake-by-wire system.
    pub fn new(max_pressure: f64) -> Self {
        let _ = max_pressure;
        Self {
            target_pressure: [ATMOSPHERIC_PRESSURE; 4],
            actual_pressure: [ATMOSPHERIC_PRESSURE; 4],
            pressure_rate: 20_000_000.0, // 20 MPa/s
            latency: 0.005,
            latency_buffer: [0.0; 4],
            hydraulic_fallback: false,
        }
    }

    /// Set target pressure for a wheel.
    pub fn set_target(&mut self, wheel: usize, pressure: f64) {
        if wheel < 4 {
            self.target_pressure[wheel] = pressure.max(ATMOSPHERIC_PRESSURE);
        }
    }

    /// Update all wheels' actual pressures towards targets.
    pub fn update(&mut self, dt: f64) {
        for i in 0..4 {
            self.latency_buffer[i] += dt;
            if self.latency_buffer[i] < self.latency {
                continue;
            }
            let diff = self.target_pressure[i] - self.actual_pressure[i];
            let max_change = self.pressure_rate * dt;
            self.actual_pressure[i] += diff.clamp(-max_change, max_change);
        }
    }

    /// Activate hydraulic fallback mode.
    pub fn activate_fallback(&mut self) {
        self.hydraulic_fallback = true;
    }
}

// ---------------------------------------------------------------------------
// Complete Brake System
// ---------------------------------------------------------------------------

/// Per-corner brake assembly.
#[derive(Debug, Clone)]
pub struct CornerBrake {
    /// Caliper for this corner.
    pub caliper: BrakeCaliper,
    /// Rotor for this corner.
    pub rotor: BrakeRotor,
    /// Current applied pressure in Pa.
    pub pressure: f64,
}

impl CornerBrake {
    /// Create a corner brake assembly.
    pub fn new(caliper: BrakeCaliper, rotor: BrakeRotor) -> Self {
        Self {
            caliper,
            rotor,
            pressure: ATMOSPHERIC_PRESSURE,
        }
    }

    /// Braking torque at this corner (N·m).
    pub fn torque(&self) -> f64 {
        self.caliper.braking_torque(self.pressure, &self.rotor)
    }

    /// Braking force at this corner (N) given wheel radius.
    pub fn braking_force(&self, wheel_radius: f64) -> f64 {
        if wheel_radius < 1e-6 {
            return 0.0;
        }
        self.torque() / wheel_radius
    }

    /// Update thermal state given wheel speed and time step.
    pub fn update_thermal(&mut self, wheel_speed_rad_s: f64, ambient: f64, dt: f64) {
        let heat = self
            .caliper
            .heat_generation_rate(self.pressure, &self.rotor, wheel_speed_rad_s);
        // Distribute heat 60/40 rotor/pad
        self.rotor.update_temperature(heat * 0.6, ambient, dt);
        self.caliper.update_pad_temperature(heat * 0.4, ambient, dt);
    }
}

/// Full 4-wheel vehicle brake system.
#[derive(Debug, Clone)]
pub struct BrakeSystem {
    /// Four corner assemblies (FL, FR, RL, RR).
    pub corners: [CornerBrake; 4],
    /// Front hydraulic circuit.
    pub front_circuit: HydraulicCircuit,
    /// Rear hydraulic circuit.
    pub rear_circuit: HydraulicCircuit,
    /// Brake bias controller.
    pub bias: BrakeBiasController,
    /// ABS controller.
    pub abs: AbsBrakeController,
    /// Regenerative braking (optional).
    pub regen: Option<RegenerativeBraking>,
    /// Pedal simulator.
    pub pedal: BrakePedalSimulator,
    /// Brake-by-wire (optional).
    pub bbw: Option<BrakeByWire>,
    /// Total energy dissipated by brakes in J.
    pub energy_dissipated_j: f64,
}

impl BrakeSystem {
    /// Create a default 4-wheel hydraulic brake system.
    pub fn new_default() -> Self {
        let pad = PadMaterial::semi_metallic();
        let mk_corner = |piston_d: f64, rotor_r: f64| {
            let cal = BrakeCaliper::new(4, piston_d, pad);
            let rot = BrakeRotor::new_vented(rotor_r, rotor_r * 0.55, 0.028, 8.0);
            CornerBrake::new(cal, rot)
        };

        Self {
            corners: [
                mk_corner(0.042, 0.160), // FL
                mk_corner(0.042, 0.160), // FR
                mk_corner(0.034, 0.130), // RL
                mk_corner(0.034, 0.130), // RR
            ],
            front_circuit: HydraulicCircuit::new(0.022, 80e-6),
            rear_circuit: HydraulicCircuit::new(0.019, 60e-6),
            bias: BrakeBiasController::new(0.60),
            abs: AbsBrakeController::new(),
            regen: None,
            pedal: BrakePedalSimulator::new(),
            bbw: None,
            energy_dissipated_j: 0.0,
        }
    }

    /// Apply driver brake input and simulate one timestep.
    ///
    /// `pedal_force` – force applied to brake pedal in N.
    /// `wheel_speeds` – angular velocity of each wheel (rad/s) \[FL, FR, RL, RR\].
    /// `vehicle_speed` – forward speed at wheel contact patches (m/s).
    /// `wheel_radius` – wheel rolling radius (m).
    /// `ambient_temp` – ambient temperature °C.
    /// `dt` – simulation time step in seconds.
    ///
    /// Returns braking torques \[FL, FR, RL, RR\] in N·m.
    pub fn step(
        &mut self,
        pedal_force: f64,
        wheel_speeds: [f64; 4],
        vehicle_speed: f64,
        wheel_radius: f64,
        ambient_temp: f64,
        dt: f64,
    ) -> [f64; 4] {
        // Step 1: pedal to master cylinder pressure
        let pedal_pos = self.pedal.step(pedal_force, dt);
        let total_pressure = self.pedal.pressure_from_travel(pedal_pos);

        // Step 2: split via bias
        let (p_front, p_rear) = self.bias.split_pressure(total_pressure);

        // Step 3: update circuit pressures
        self.front_circuit.pressure = ATMOSPHERIC_PRESSURE + p_front;
        self.rear_circuit.pressure = ATMOSPHERIC_PRESSURE + p_rear;

        let base_pressures = [p_front, p_front, p_rear, p_rear];

        // Step 4: ABS modulation
        let mut pressures = base_pressures;
        for i in 0..4 {
            let v_wheel = wheel_speeds[i] * wheel_radius;
            let slip = AbsBrakeController::wheel_slip(vehicle_speed, v_wheel);
            let dp = self
                .abs
                .regulate_wheel(i, ATMOSPHERIC_PRESSURE + pressures[i], slip, dt);
            pressures[i] = (pressures[i] + dp).max(0.0);
        }

        // Step 5: apply to corners
        let mut torques = [0.0_f64; 4];
        for i in 0..4 {
            self.corners[i].pressure = ATMOSPHERIC_PRESSURE + pressures[i];
            torques[i] = self.corners[i]
                .caliper
                .braking_torque(pressures[i], &self.corners[i].rotor);
            self.corners[i].update_thermal(wheel_speeds[i], ambient_temp, dt);
            let heat = torques[i] * wheel_speeds[i].abs();
            self.energy_dissipated_j += heat * dt;
        }

        torques
    }

    /// Enable regenerative braking integration.
    pub fn enable_regen(&mut self, regen: RegenerativeBraking) {
        self.regen = Some(regen);
    }

    /// Enable brake-by-wire.
    pub fn enable_bbw(&mut self, bbw: BrakeByWire) {
        self.bbw = Some(bbw);
    }

    /// Overall brake temperature check (returns max rotor temp °C).
    pub fn max_rotor_temperature(&self) -> f64 {
        self.corners
            .iter()
            .map(|c| c.rotor.temperature)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum pad life fraction across all corners.
    pub fn min_pad_life(&self) -> f64 {
        self.corners
            .iter()
            .map(|c| c.caliper.pad_life_fraction())
            .fold(f64::INFINITY, f64::min)
    }

    /// Total braking power dissipated currently in W.
    pub fn total_braking_power(&self, wheel_speeds: [f64; 4]) -> f64 {
        self.corners
            .iter()
            .zip(wheel_speeds.iter())
            .map(|(c, &ws)| {
                c.caliper
                    .braking_torque(c.pressure - ATMOSPHERIC_PRESSURE, &c.rotor)
                    * ws.abs()
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Brake Force Distribution
// ---------------------------------------------------------------------------

/// Electronic Brake Force Distribution (EBD) system.
#[derive(Debug, Clone)]
pub struct EbdSystem {
    /// Maximum longitudinal deceleration in m/s² before rear bias reduction.
    pub rear_threshold_decel: f64,
    /// Ramp rate for rear pressure reduction (Pa per m/s²).
    pub rear_ramp: f64,
    /// Enabled flag.
    pub enabled: bool,
}

impl EbdSystem {
    /// Create default EBD system.
    pub fn new() -> Self {
        Self {
            rear_threshold_decel: 5.0,
            rear_ramp: 200_000.0,
            enabled: true,
        }
    }

    /// Compute rear pressure multiplier based on deceleration.
    pub fn rear_pressure_factor(&self, deceleration_m_s2: f64) -> f64 {
        if !self.enabled {
            return 1.0;
        }
        if deceleration_m_s2 <= self.rear_threshold_decel {
            return 1.0;
        }
        let excess = deceleration_m_s2 - self.rear_threshold_decel;
        let reduction = self.rear_ramp * excess;
        let baseline = 6_000_000.0_f64; // assume 6 MPa baseline
        (1.0 - reduction / baseline).clamp(0.4, 1.0)
    }
}

impl Default for EbdSystem {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Brake Fade Index
// ---------------------------------------------------------------------------

/// Compute a 0–1 fade index (0 = no fade, 1 = complete fade) from rotor temperature.
pub fn fade_index(rotor_temp_c: f64, pad_material: &PadMaterial) -> f64 {
    if rotor_temp_c <= pad_material.fade_onset_temp {
        return 0.0;
    }
    let excess = rotor_temp_c - pad_material.fade_onset_temp;
    (excess / 400.0).clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Brake Thermal Stress
// ---------------------------------------------------------------------------

/// Compute peak thermal stress in a disc brake rotor.
///
/// Returns stress in Pa.
pub fn rotor_thermal_stress(rotor: &BrakeRotor, ambient_temp: f64) -> f64 {
    let delta_t = (rotor.temperature - ambient_temp).max(0.0);
    rotor.thermal_hoop_stress(delta_t)
}

// ---------------------------------------------------------------------------
// Stopping Distance Estimate
// ---------------------------------------------------------------------------

/// Estimate braking stopping distance.
///
/// `v0` – initial speed in m/s.
/// `decel` – mean deceleration in m/s².
/// Returns stopping distance in m.
pub fn stopping_distance(v0: f64, decel: f64) -> f64 {
    if decel < 1e-6 {
        return f64::INFINITY;
    }
    v0 * v0 / (2.0 * decel)
}

/// Required deceleration for a given speed and stopping distance.
pub fn required_deceleration(v0: f64, distance: f64) -> f64 {
    if distance < 1e-6 {
        return f64::INFINITY;
    }
    v0 * v0 / (2.0 * distance)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Pad material tests
    #[test]
    fn test_organic_pad_cold_friction() {
        let pad = PadMaterial::organic();
        let mu = pad.friction_at_temp(25.0);
        // At 25 °C the model linearly interpolates between friction_cold (0.38)
        // and friction_hot (0.42) over the range [0, peak_temp=300], so mu ≈ 0.383.
        assert!((mu - 0.38).abs() < 0.01);
    }

    #[test]
    fn test_organic_pad_hot_friction() {
        let pad = PadMaterial::organic();
        let mu = pad.friction_at_temp(300.0);
        assert!((mu - 0.42).abs() < 1e-3);
    }

    #[test]
    fn test_pad_fade_high_temp() {
        let pad = PadMaterial::organic();
        let mu_at_800 = pad.friction_at_temp(800.0);
        let mu_at_300 = pad.friction_at_temp(300.0);
        assert!(
            mu_at_800 < mu_at_300,
            "No fade detected at extreme temperature"
        );
    }

    #[test]
    fn test_pad_fade_floor_respected() {
        let pad = PadMaterial::organic();
        let mu = pad.friction_at_temp(2000.0);
        assert!(mu >= pad.friction_fade_floor - 1e-9);
    }

    #[test]
    fn test_semi_metallic_higher_peak_than_organic() {
        let organic = PadMaterial::organic();
        let semi = PadMaterial::semi_metallic();
        assert!(semi.fade_onset_temp > organic.fade_onset_temp);
    }

    // Rotor tests
    #[test]
    fn test_rotor_effective_radius() {
        let r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        let eff = r.effective_radius();
        assert!((eff - 0.124).abs() < 1e-9);
    }

    #[test]
    fn test_rotor_swept_area_positive() {
        let r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        assert!(r.swept_area() > 0.0);
    }

    #[test]
    fn test_rotor_cooling_area_positive() {
        let r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        assert!(r.cooling_area() > r.swept_area());
    }

    #[test]
    fn test_rotor_temperature_rises_under_heat() {
        let mut r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        let init = r.temperature;
        r.update_temperature(50_000.0, 25.0, 0.1);
        assert!(r.temperature > init);
    }

    #[test]
    fn test_rotor_cools_to_ambient_without_heat() {
        let mut r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        r.temperature = 400.0;
        for _ in 0..10000 {
            r.update_temperature(0.0, 25.0, 0.1);
        }
        assert!(r.temperature < 100.0);
    }

    #[test]
    fn test_rotor_thermal_stress_proportional() {
        let r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        let s1 = r.thermal_hoop_stress(100.0);
        let s2 = r.thermal_hoop_stress(200.0);
        assert!((s2 / s1 - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_rotor_wear_accumulates() {
        let mut r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        r.apply_wear(1e6, 1e-14);
        assert!(r.wear_depth > 0.0);
    }

    // Caliper tests
    #[test]
    fn test_caliper_piston_area_positive() {
        let cal = BrakeCaliper::new(4, 0.042, PadMaterial::organic());
        assert!(cal.piston_area() > 0.0);
    }

    #[test]
    fn test_caliper_clamping_increases_with_pressure() {
        let cal = BrakeCaliper::new(4, 0.042, PadMaterial::organic());
        let f1 = cal.clamping_force(1_000_000.0);
        let f2 = cal.clamping_force(2_000_000.0);
        assert!(f2 > f1);
    }

    #[test]
    fn test_caliper_braking_torque_positive() {
        let cal = BrakeCaliper::new(4, 0.042, PadMaterial::organic());
        let rot = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        let t = cal.braking_torque(1_000_000.0, &rot);
        assert!(t > 0.0);
    }

    #[test]
    fn test_caliper_pad_life_starts_full() {
        let cal = BrakeCaliper::new(4, 0.042, PadMaterial::organic());
        assert!((cal.pad_life_fraction() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_caliper_pad_wear_reduces_thickness() {
        let mut cal = BrakeCaliper::new(4, 0.042, PadMaterial::organic());
        cal.apply_pad_wear(10_000.0, 1000.0);
        assert!(cal.pad_thickness < cal.initial_pad_thickness);
    }

    // Hydraulic circuit tests
    #[test]
    fn test_hydraulic_pressure_rises_with_displacement() {
        let mut hc = HydraulicCircuit::new(0.022, 80e-6);
        let p = hc.update_pressure(0.01);
        assert!(p > ATMOSPHERIC_PRESSURE);
    }

    #[test]
    fn test_hydraulic_effective_pressure_normal() {
        let hc = HydraulicCircuit::new(0.022, 80e-6);
        assert!((hc.effective_pressure() - hc.pressure).abs() < 1e-6);
    }

    #[test]
    fn test_hydraulic_vapor_lock() {
        let mut hc = HydraulicCircuit::new(0.022, 80e-6);
        hc.fluid_temperature = 250.0;
        hc.pressure = 5_000_000.0;
        assert_eq!(hc.effective_pressure(), ATMOSPHERIC_PRESSURE);
    }

    // Bias controller tests
    #[test]
    fn test_bias_fractions_sum_to_one() {
        let b = BrakeBiasController::new(0.6);
        assert!((b.front_fraction() + b.rear_fraction() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_bias_pressure_split_consistent() {
        let b = BrakeBiasController::new(0.65);
        let (pf, pr) = b.split_pressure(6_000_000.0);
        assert!((pf + pr - 6_000_000.0).abs() < 1e-3);
    }

    #[test]
    fn test_bias_adjust_clamped() {
        let mut b = BrakeBiasController::new(0.6);
        b.adjust(1.0); // try to go way over
        assert!(b.front_bias <= b.max_bias);
    }

    // ABS tests
    #[test]
    fn test_abs_wheel_slip_zero_at_zero_speed() {
        let slip = AbsBrakeController::wheel_slip(0.0, 0.0);
        assert_eq!(slip, 0.0);
    }

    #[test]
    fn test_abs_wheel_slip_full_lock() {
        let slip = AbsBrakeController::wheel_slip(30.0, 0.0);
        assert!((slip - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_abs_no_intervention_below_threshold() {
        let mut abs = AbsBrakeController::new();
        let dp = abs.regulate_wheel(0, 2_000_000.0, 0.05, 0.01);
        assert_eq!(dp, 0.0); // hold
    }

    #[test]
    fn test_abs_activates_over_threshold() {
        let mut abs = AbsBrakeController::new();
        // Trigger hold time first
        abs.cycle_timer[0] = 0.1;
        let dp = abs.regulate_wheel(0, 2_000_000.0, 0.30, 0.01);
        assert!(dp < 0.0, "ABS should release pressure when slip is high");
    }

    // Regen braking tests
    #[test]
    fn test_regen_no_torque_at_zero_speed() {
        let regen = RegenerativeBraking::new(500.0, 0.9, 50_000.0);
        assert_eq!(regen.available_torque(0.0), 0.0);
    }

    #[test]
    fn test_regen_no_torque_at_max_soc() {
        let mut regen = RegenerativeBraking::new(500.0, 0.9, 50_000.0);
        regen.soc = 0.95;
        assert_eq!(regen.available_torque(100.0), 0.0);
    }

    #[test]
    fn test_regen_positive_torque_normal_conditions() {
        let regen = RegenerativeBraking::new(500.0, 0.9, 50_000.0);
        let t = regen.available_torque(50.0);
        assert!(t > 0.0);
    }

    #[test]
    fn test_regen_blend_conserves_total() {
        let regen = RegenerativeBraking::new(500.0, 0.9, 50_000.0);
        let demanded = 300.0;
        let (f, r) = regen.blend_torques(demanded, 50.0);
        assert!((f + r - demanded).abs() < 1e-6);
    }

    // Pedal simulator tests
    #[test]
    fn test_pedal_zero_force_zero_travel() {
        let mut p = BrakePedalSimulator::new();
        let pos = p.step(0.0, 0.01);
        assert_eq!(pos, 0.0);
    }

    #[test]
    fn test_pedal_positive_force_increases_position() {
        let mut p = BrakePedalSimulator::new();
        for _ in 0..100 {
            p.step(200.0, 0.01);
        }
        assert!(p.position > 0.0);
    }

    #[test]
    fn test_pedal_pressure_from_travel_zero_at_free_play() {
        let p = BrakePedalSimulator::new();
        let pressure = p.pressure_from_travel(p.free_play);
        assert_eq!(pressure, 0.0);
    }

    #[test]
    fn test_pedal_pressure_increases_with_travel() {
        let p = BrakePedalSimulator::new();
        let p1 = p.pressure_from_travel(0.02);
        let p2 = p.pressure_from_travel(0.05);
        assert!(p2 > p1);
    }

    // Brake system integration tests
    #[test]
    fn test_brake_system_torques_positive_under_braking() {
        let mut bs = BrakeSystem::new_default();
        let torques = bs.step(400.0, [100.0; 4], 30.0, 0.32, 25.0, 0.01);
        assert!(torques.iter().all(|&t| t >= 0.0));
    }

    #[test]
    fn test_brake_system_no_torque_no_pedal() {
        let mut bs = BrakeSystem::new_default();
        let torques = bs.step(0.0, [100.0; 4], 30.0, 0.32, 25.0, 0.01);
        assert!(torques.iter().all(|&t| t < 1.0));
    }

    #[test]
    fn test_brake_system_max_rotor_temp_ambient_initially() {
        let bs = BrakeSystem::new_default();
        assert!(bs.max_rotor_temperature() < 30.0);
    }

    #[test]
    fn test_brake_system_pad_life_starts_full() {
        let bs = BrakeSystem::new_default();
        assert!((bs.min_pad_life() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_brake_system_energy_accumulates() {
        let mut bs = BrakeSystem::new_default();
        for _ in 0..100 {
            bs.step(400.0, [80.0; 4], 25.0, 0.32, 25.0, 0.01);
        }
        assert!(bs.energy_dissipated_j > 0.0);
    }

    // EBD tests
    #[test]
    fn test_ebd_no_reduction_below_threshold() {
        let ebd = EbdSystem::new();
        assert!((ebd.rear_pressure_factor(3.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_ebd_reduction_above_threshold() {
        let ebd = EbdSystem::new();
        let factor = ebd.rear_pressure_factor(20.0);
        assert!(factor < 1.0);
    }

    #[test]
    fn test_ebd_disabled_returns_one() {
        let mut ebd = EbdSystem::new();
        ebd.enabled = false;
        assert!((ebd.rear_pressure_factor(30.0) - 1.0).abs() < 1e-9);
    }

    // Helper function tests
    #[test]
    fn test_stopping_distance_basic() {
        let d = stopping_distance(30.0, 9.0);
        assert!((d - 50.0).abs() < 0.1);
    }

    #[test]
    fn test_required_deceleration_basic() {
        let a = required_deceleration(30.0, 50.0);
        assert!((a - 9.0).abs() < 0.1);
    }

    #[test]
    fn test_fade_index_zero_below_onset() {
        let pad = PadMaterial::organic();
        assert_eq!(fade_index(300.0, &pad), 0.0);
    }

    #[test]
    fn test_fade_index_positive_above_onset() {
        let pad = PadMaterial::organic();
        let fi = fade_index(600.0, &pad);
        assert!(fi > 0.0);
    }

    #[test]
    fn test_rotor_thermal_stress_from_function() {
        let mut r = BrakeRotor::new(0.16, 0.088, 0.028, 8.0);
        r.temperature = 300.0;
        let stress = rotor_thermal_stress(&r, 25.0);
        assert!(stress > 0.0);
    }

    #[test]
    fn test_bbw_pressure_follows_target() {
        let mut bbw = BrakeByWire::new(10_000_000.0);
        bbw.set_target(0, 3_000_000.0);
        for _ in 0..100 {
            bbw.update(0.01);
        }
        assert!(bbw.actual_pressure[0] > ATMOSPHERIC_PRESSURE);
    }

    #[test]
    fn test_optimal_bias_increases_with_decel() {
        let b1 = BrakeBiasController::optimal_bias(0.5, 0.55, 0.45, 2.6);
        let b2 = BrakeBiasController::optimal_bias(1.0, 0.55, 0.45, 2.6);
        assert!(b2 >= b1);
    }

    #[test]
    fn test_carbon_ceramic_high_fade_onset() {
        let cc = PadMaterial::carbon_ceramic();
        assert!(cc.fade_onset_temp > 800.0);
    }
}
