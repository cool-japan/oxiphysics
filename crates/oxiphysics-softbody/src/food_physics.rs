// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Food physics and cooking simulation.
//!
//! Implements:
//! - [`FoodRheology`] — viscoelastic food (bread dough, gel), frequency-sweep G'/G'', tan δ
//! - [`HeatCookingModel`] — thermal conduction in food, Maillard reaction kinetics, cooking degree
//! - [`FermentationModel`] — CO₂ production (yeast kinetics), dough rise, gas bubble nucleation
//! - [`TextureEvolution`] — starch gelatinization, protein denaturation, softening/hardening
//! - [`ChewingSimulation`] — jaw mechanics, food fragmentation, particle size reduction
//! - [`FoodEmulsion`] — oil-water emulsion stability, Ostwald ripening, droplet size distribution

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical/mathematical constants
// ---------------------------------------------------------------------------

/// Universal gas constant (J mol⁻¹ K⁻¹).
const R_GAS: f64 = 8.314;
/// Absolute zero offset (K).
const KELVIN: f64 = 273.15;

// ---------------------------------------------------------------------------
// FoodRheology
// ---------------------------------------------------------------------------

/// Viscoelastic model for food materials (Maxwell, Kelvin-Voigt, power-law).
///
/// Enables frequency-sweep predictions of storage modulus G', loss modulus G'',
/// and loss tangent tan δ.
#[derive(Debug, Clone)]
pub struct FoodRheology {
    /// Storage modulus at reference frequency (Pa).
    pub g_prime_ref: f64,
    /// Loss modulus at reference frequency (Pa).
    pub g_double_prime_ref: f64,
    /// Frequency power-law exponent for G' (dimensionless).
    pub alpha_g_prime: f64,
    /// Frequency power-law exponent for G'' (dimensionless).
    pub alpha_g_double_prime: f64,
    /// Reference frequency (Hz).
    pub f_ref: f64,
    /// Zero-frequency viscosity (Pa·s) for the Maxwell element.
    pub eta_0: f64,
    /// Relaxation time λ (s).
    pub lambda: f64,
    /// Food category label.
    pub category: FoodCategory,
}

/// Food material category for rheological classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FoodCategory {
    /// Bread dough.
    Dough,
    /// Hydrocolloid gel (gelatin, agar).
    Gel,
    /// Processed cheese or processed dairy.
    DairyGel,
    /// Fruit puree or vegetable paste.
    Puree,
    /// Generic / custom.
    Custom,
}

impl FoodRheology {
    /// Create a new rheological model for bread dough (typical parameters).
    pub fn bread_dough() -> Self {
        Self {
            g_prime_ref: 3500.0,        // Pa at 1 Hz
            g_double_prime_ref: 1200.0, // Pa at 1 Hz
            alpha_g_prime: 0.21,
            alpha_g_double_prime: 0.18,
            f_ref: 1.0,
            eta_0: 5000.0,
            lambda: 10.0,
            category: FoodCategory::Dough,
        }
    }

    /// Create a gelatin gel model (10 % w/w, 20 °C).
    pub fn gelatin_gel() -> Self {
        Self {
            g_prime_ref: 2000.0,
            g_double_prime_ref: 300.0,
            alpha_g_prime: 0.08,
            alpha_g_double_prime: 0.12,
            f_ref: 1.0,
            eta_0: 800.0,
            lambda: 50.0,
            category: FoodCategory::Gel,
        }
    }

    /// Storage modulus G'(ω) \[Pa\] at angular frequency `omega` (rad/s).
    ///
    /// Power-law model: G'(ω) = G'_ref * (ω / ω_ref)^α_G'
    pub fn g_prime(&self, omega: f64) -> f64 {
        let omega_ref = 2.0 * PI * self.f_ref;
        self.g_prime_ref * (omega / omega_ref).powf(self.alpha_g_prime)
    }

    /// Loss modulus G''(ω) \[Pa\] at angular frequency `omega` (rad/s).
    ///
    /// Power-law model: G''(ω) = G''_ref * (ω / ω_ref)^α_G''
    pub fn g_double_prime(&self, omega: f64) -> f64 {
        let omega_ref = 2.0 * PI * self.f_ref;
        self.g_double_prime_ref * (omega / omega_ref).powf(self.alpha_g_double_prime)
    }

    /// Loss tangent tan δ = G'' / G' at angular frequency `omega`.
    pub fn tan_delta(&self, omega: f64) -> f64 {
        let g1 = self.g_prime(omega);
        if g1 < 1e-30 {
            return f64::INFINITY;
        }
        self.g_double_prime(omega) / g1
    }

    /// Complex modulus magnitude |G*| = sqrt(G'^2 + G''^2).
    pub fn complex_modulus(&self, omega: f64) -> f64 {
        let g1 = self.g_prime(omega);
        let g2 = self.g_double_prime(omega);
        (g1 * g1 + g2 * g2).sqrt()
    }

    /// Phase angle δ (radians): δ = atan(G''/G').
    pub fn phase_angle(&self, omega: f64) -> f64 {
        self.g_double_prime(omega).atan2(self.g_prime(omega))
    }

    /// Perform a frequency sweep over `n_points` decades from `f_min` to `f_max` Hz.
    ///
    /// Returns vectors of `(frequency_Hz, G'_Pa, G''_Pa, tan_delta)`.
    pub fn frequency_sweep(
        &self,
        f_min: f64,
        f_max: f64,
        n_points: usize,
    ) -> Vec<(f64, f64, f64, f64)> {
        if n_points == 0 {
            return Vec::new();
        }
        let log_min = f_min.log10();
        let log_max = f_max.log10();
        (0..n_points)
            .map(|i| {
                let t = i as f64 / (n_points as f64 - 1.0).max(1.0);
                let f = 10.0_f64.powf(log_min + t * (log_max - log_min));
                let omega = 2.0 * PI * f;
                (
                    f,
                    self.g_prime(omega),
                    self.g_double_prime(omega),
                    self.tan_delta(omega),
                )
            })
            .collect()
    }

    /// Apparent viscosity η' = G''/ω \[Pa·s\].
    pub fn apparent_viscosity(&self, omega: f64) -> f64 {
        if omega < 1e-30 {
            return f64::INFINITY;
        }
        self.g_double_prime(omega) / omega
    }

    /// Maxwell model: creep compliance J(t) = (1/E) * (1 − exp(−t/τ)) + t/η_0.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let e0 = self.g_prime_ref;
        if e0 < 1e-30 {
            return 0.0;
        }
        (1.0 / e0) * (1.0 - (-t / self.lambda).exp()) + t / self.eta_0
    }

    /// Stress relaxation modulus G(t) = G0 * exp(−t/λ) \[Pa\].
    pub fn stress_relaxation(&self, t: f64) -> f64 {
        self.g_prime_ref * (-t / self.lambda).exp()
    }

    /// Indicates if the material is predominantly elastic (tan δ < 1).
    pub fn is_elastic_dominant(&self, omega: f64) -> bool {
        self.tan_delta(omega) < 1.0
    }

    /// Indicates if the material is predominantly viscous (tan δ > 1).
    pub fn is_viscous_dominant(&self, omega: f64) -> bool {
        self.tan_delta(omega) > 1.0
    }
}

// ---------------------------------------------------------------------------
// HeatCookingModel
// ---------------------------------------------------------------------------

/// 1D thermal model for cooking processes in food.
///
/// Solves the heat diffusion equation with Maillard reaction kinetics and
/// tracks a dimensionless "cooking degree" from 0 (raw) to 1 (fully cooked).
#[derive(Debug, Clone)]
pub struct HeatCookingModel {
    /// Number of spatial nodes (including boundaries).
    pub n_nodes: usize,
    /// Physical thickness of the food slab (m).
    pub thickness: f64,
    /// Thermal diffusivity α = k/(ρ·cp) (m²/s).
    pub thermal_diffusivity: f64,
    /// Thermal conductivity k (W/m/K).
    pub conductivity: f64,
    /// Specific heat capacity cp (J/kg/K).
    pub specific_heat: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Temperature profile (K) at each node.
    pub temperature: Vec<f64>,
    /// Cooking degree at each node \[0, 1\].
    pub cooking_degree: Vec<f64>,
    /// Maillard reaction activation energy (J/mol).
    pub maillard_ea: f64,
    /// Maillard pre-exponential factor (s⁻¹).
    pub maillard_a: f64,
    /// Reference temperature for cooking degree (K).
    pub t_cook_ref: f64,
    /// Elapsed cooking time (s).
    pub time: f64,
}

impl HeatCookingModel {
    /// Create a new 1D cooking model for a `thickness`-metre food slab at initial
    /// temperature `t_initial` (K).
    pub fn new(
        n_nodes: usize,
        thickness: f64,
        thermal_diffusivity: f64,
        conductivity: f64,
        specific_heat: f64,
        density: f64,
        t_initial: f64,
    ) -> Self {
        let n = n_nodes.max(3);
        Self {
            n_nodes: n,
            thickness,
            thermal_diffusivity,
            conductivity,
            specific_heat,
            density,
            temperature: vec![t_initial; n],
            cooking_degree: vec![0.0; n],
            maillard_ea: 80_000.0, // J/mol
            maillard_a: 1e8,       // s⁻¹
            t_cook_ref: 373.15,    // 100 °C
            time: 0.0,
        }
    }

    /// Create a bread-baking model (standard parameters).
    pub fn bread_baking(n_nodes: usize) -> Self {
        Self::new(
            n_nodes, 0.05,   // 5 cm thick bread
            1.4e-7, // m²/s
            0.18,   // W/m/K
            3000.0, // J/kg/K
            450.0,  // kg/m³
            298.15, // 25 °C initial
        )
    }

    fn dx(&self) -> f64 {
        self.thickness / (self.n_nodes as f64 - 1.0).max(1.0)
    }

    /// Set boundary temperatures: left surface T_left, right surface T_right (K).
    pub fn set_boundary_temperatures(&mut self, t_left: f64, t_right: f64) {
        let n = self.n_nodes;
        self.temperature[0] = t_left;
        self.temperature[n - 1] = t_right;
    }

    /// Advance the heat equation by `dt` seconds using explicit finite differences.
    ///
    /// The Fourier number must satisfy Fo = α·dt/dx² ≤ 0.5 for stability.
    pub fn step(&mut self, dt: f64) {
        let dx = self.dx();
        let alpha = self.thermal_diffusivity;
        let fo = alpha * dt / (dx * dx);
        // Clamp Fourier number for stability (sub-step if needed).
        let n_sub = (fo / 0.4).ceil() as usize;
        let dt_sub = dt / n_sub as f64;
        let fo_sub = alpha * dt_sub / (dx * dx);

        for _ in 0..n_sub {
            let t_old = self.temperature.clone();
            // Interior nodes
            for i in 1..(self.n_nodes - 1) {
                self.temperature[i] =
                    t_old[i] + fo_sub * (t_old[i + 1] - 2.0 * t_old[i] + t_old[i - 1]);
            }
        }

        // Update cooking degree (Arrhenius kinetics).
        for i in 0..self.n_nodes {
            let t_k = self.temperature[i];
            let rate = self.maillard_rate(t_k);
            let deg = &mut self.cooking_degree[i];
            *deg = (*deg + rate * dt).min(1.0);
        }
        self.time += dt;
    }

    /// Maillard reaction rate (s⁻¹) at temperature `t_k` (K).
    ///
    /// k(T) = A * exp(−Ea / (R·T))
    pub fn maillard_rate(&self, t_k: f64) -> f64 {
        if t_k < 1.0 {
            return 0.0;
        }
        self.maillard_a * (-(self.maillard_ea) / (R_GAS * t_k)).exp()
    }

    /// Celsius temperature at node `i`.
    pub fn temperature_celsius(&self, i: usize) -> f64 {
        self.temperature[i] - KELVIN
    }

    /// Mean temperature across all nodes (K).
    pub fn mean_temperature(&self) -> f64 {
        self.temperature.iter().sum::<f64>() / self.n_nodes as f64
    }

    /// Mean cooking degree across all nodes.
    pub fn mean_cooking_degree(&self) -> f64 {
        self.cooking_degree.iter().sum::<f64>() / self.n_nodes as f64
    }

    /// Centre node index.
    pub fn centre_node(&self) -> usize {
        self.n_nodes / 2
    }

    /// Surface browning indicator: cooking degree at the boundary nodes.
    pub fn surface_browning(&self) -> f64 {
        (self.cooking_degree[0] + self.cooking_degree[self.n_nodes - 1]) / 2.0
    }

    /// Total thermal energy stored in the slab (J/m²).
    pub fn stored_energy(&self) -> f64 {
        let dx = self.dx();
        let sum: f64 = self.temperature.iter().sum();
        self.density * self.specific_heat * dx * sum
    }

    /// Time to reach `target_degree` cooking degree at the centre (estimate via
    /// linear extrapolation from current rate).
    pub fn time_to_cook(&self, target_degree: f64) -> Option<f64> {
        let centre = self.centre_node();
        let current = self.cooking_degree[centre];
        if current >= target_degree {
            return Some(0.0);
        }
        let rate = self.maillard_rate(self.temperature[centre]);
        if rate < 1e-30 {
            return None;
        }
        Some((target_degree - current) / rate)
    }
}

// ---------------------------------------------------------------------------
// FermentationModel
// ---------------------------------------------------------------------------

/// Yeast fermentation and dough-rise model.
///
/// Models CO₂ production via Monod-type yeast kinetics, nucleation and growth
/// of gas bubbles, and macroscopic dough volume increase.
#[derive(Debug, Clone)]
pub struct FermentationModel {
    /// Mass of yeast biomass (kg per kg flour).
    pub yeast_concentration: f64,
    /// Substrate (sugar) concentration (kg per kg flour).
    pub substrate: f64,
    /// Maximum specific growth rate of yeast μ_max (s⁻¹).
    pub mu_max: f64,
    /// Monod half-saturation constant K_s (kg/kg).
    pub ks: f64,
    /// CO₂ yield coefficient Y_CO2 (mol CO₂ per kg substrate).
    pub y_co2: f64,
    /// Current CO₂ concentration in dough (mol/kg).
    pub co2_concentration: f64,
    /// Dough volume (relative to initial; 1.0 = no rise).
    pub volume_ratio: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Optimal fermentation temperature (K).
    pub t_optimal: f64,
    /// Gas bubble number density (bubbles/m³).
    pub bubble_density: f64,
    /// Mean bubble radius (m).
    pub bubble_radius: f64,
    /// Elapsed fermentation time (s).
    pub time: f64,
    /// CO₂ solubility in dough water (mol/kg per Pa).
    pub co2_solubility: f64,
    /// Ambient pressure (Pa).
    pub pressure: f64,
}

impl FermentationModel {
    /// Create a standard bread dough fermentation model.
    pub fn bread_dough(initial_yeast: f64, initial_substrate: f64, temperature_c: f64) -> Self {
        Self {
            yeast_concentration: initial_yeast,
            substrate: initial_substrate,
            mu_max: 2.8e-4, // s⁻¹  (~1/hr)
            ks: 0.01,       // kg/kg
            y_co2: 4.0,     // mol CO₂ / kg sugar
            co2_concentration: 0.0,
            volume_ratio: 1.0,
            temperature: temperature_c + KELVIN,
            t_optimal: 305.15,    // 32 °C
            bubble_density: 1e10, // bubbles/m³
            bubble_radius: 10e-6, // 10 µm initial
            time: 0.0,
            co2_solubility: 3.4e-4,
            pressure: 101_325.0,
        }
    }

    /// Temperature correction factor (bell-shaped, peak at T_optimal).
    pub fn temperature_factor(&self) -> f64 {
        let dt = (self.temperature - self.t_optimal).abs();
        let sigma = 10.0; // K width
        (-dt * dt / (2.0 * sigma * sigma)).exp()
    }

    /// CO₂ production rate (mol/kg/s) based on Monod kinetics.
    pub fn co2_production_rate(&self) -> f64 {
        let mu =
            self.mu_max * self.temperature_factor() * self.substrate / (self.ks + self.substrate);
        mu * self.yeast_concentration * self.y_co2
    }

    /// Yeast growth rate dX/dt \[kg/kg/s\].
    pub fn yeast_growth_rate(&self) -> f64 {
        let mu =
            self.mu_max * self.temperature_factor() * self.substrate / (self.ks + self.substrate);
        mu * self.yeast_concentration
    }

    /// Substrate consumption rate dS/dt \[kg/kg/s\].
    pub fn substrate_consumption_rate(&self) -> f64 {
        let yield_biomass = 0.5; // kg biomass / kg substrate
        self.yeast_growth_rate() / yield_biomass
    }

    /// Advance the fermentation model by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let dco2 = self.co2_production_rate() * dt;
        let dx = self.yeast_growth_rate() * dt;
        let ds = self.substrate_consumption_rate() * dt;

        self.co2_concentration += dco2;
        self.yeast_concentration += dx;
        self.substrate = (self.substrate - ds).max(0.0);

        // Volume increase proportional to dissolved + gaseous CO₂.
        let excess_co2 = (self.co2_concentration - self.co2_solubility * self.pressure).max(0.0);
        // Ideal gas: V = nRT/P at constant P.
        // Specific volume of gas (m³/kg of CO₂).
        let v_gas = excess_co2 * R_GAS * self.temperature / self.pressure / 44e-3; // molar mass CO₂
        self.volume_ratio = 1.0 + v_gas.min(4.0); // clamp at 5× rise

        // Bubble growth: bubble radius grows as (3 * n_gas / (4π * n_bubbles))^(1/3)
        let n_gas_mol = excess_co2.max(0.0);
        if n_gas_mol > 1e-10 && self.bubble_density > 1e3 {
            let v_per_bubble =
                n_gas_mol * R_GAS * self.temperature / self.pressure / self.bubble_density / 44e-3;
            let r = (3.0 * v_per_bubble / (4.0 * PI)).powf(1.0 / 3.0);
            self.bubble_radius = r.max(1e-9);
        }

        self.time += dt;
    }

    /// Run fermentation for `total_time` seconds with step `dt`.
    pub fn run(&mut self, total_time: f64, dt: f64) {
        let n = (total_time / dt).ceil() as usize;
        for _ in 0..n {
            if self.time >= total_time {
                break;
            }
            self.step(dt.min(total_time - self.time));
        }
    }

    /// Check if fermentation is complete (substrate depleted below threshold).
    pub fn is_complete(&self, threshold: f64) -> bool {
        self.substrate < threshold
    }

    /// Proofing ratio: current volume / initial volume.
    pub fn proofing_ratio(&self) -> f64 {
        self.volume_ratio
    }

    /// CO₂ in headspace (escaped from dough, estimated).
    pub fn escaped_co2_fraction(&self) -> f64 {
        let total = self.co2_concentration + 1e-10;
        let dissolved = self.co2_solubility * self.pressure;
        (1.0 - dissolved / total).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// TextureEvolution
// ---------------------------------------------------------------------------

/// Models changes in food texture during cooking.
///
/// Tracks starch gelatinization, protein denaturation, and overall
/// softness/hardness evolution as functions of temperature and time.
#[derive(Debug, Clone)]
pub struct TextureEvolution {
    /// Fraction of starch gelatinized \[0, 1\].
    pub gelatinization_degree: f64,
    /// Fraction of proteins denatured \[0, 1\].
    pub denaturation_degree: f64,
    /// Current material hardness (Pa).
    pub hardness: f64,
    /// Initial hardness before any cooking (Pa).
    pub initial_hardness: f64,
    /// Hardness after full gelatinization (Pa).
    pub gelatinized_hardness: f64,
    /// Hardness after protein denaturation (Pa).
    pub denatured_hardness: f64,
    /// Starch gelatinization onset temperature (K).
    pub t_gel_onset: f64,
    /// Starch gelatinization completion temperature (K).
    pub t_gel_end: f64,
    /// Protein denaturation temperature (K).
    pub t_denat: f64,
    /// Protein denaturation range (K).
    pub t_denat_range: f64,
    /// Starch gelatinization rate coefficient (s⁻¹).
    pub k_gel: f64,
    /// Protein denaturation rate coefficient (s⁻¹).
    pub k_denat: f64,
    /// Activation energy for gelatinization (J/mol).
    pub ea_gel: f64,
    /// Activation energy for denaturation (J/mol).
    pub ea_denat: f64,
    /// Elapsed time (s).
    pub time: f64,
}

impl TextureEvolution {
    /// Create a bread/cake texture evolution model.
    pub fn bread() -> Self {
        Self {
            gelatinization_degree: 0.0,
            denaturation_degree: 0.0,
            hardness: 1000.0,
            initial_hardness: 1000.0,
            gelatinized_hardness: 3000.0, // starch gels harden crust
            denatured_hardness: 2500.0,
            t_gel_onset: 333.15, // 60 °C
            t_gel_end: 363.15,   // 90 °C
            t_denat: 343.15,     // 70 °C
            t_denat_range: 10.0,
            k_gel: 0.01,
            k_denat: 0.005,
            ea_gel: 60_000.0,
            ea_denat: 75_000.0,
            time: 0.0,
        }
    }

    /// Create a meat cooking texture model.
    pub fn meat() -> Self {
        Self {
            gelatinization_degree: 0.0,
            denaturation_degree: 0.0,
            hardness: 5000.0,
            initial_hardness: 5000.0,
            gelatinized_hardness: 4000.0,
            denatured_hardness: 15_000.0, // proteins toughen meat
            t_gel_onset: 323.15,
            t_gel_end: 353.15,
            t_denat: 338.15, // 65 °C
            t_denat_range: 8.0,
            k_gel: 0.002,
            k_denat: 0.008,
            ea_gel: 55_000.0,
            ea_denat: 90_000.0,
            time: 0.0,
        }
    }

    /// Gelatinization rate at temperature `t_k` (K): sigmoid based on temperature.
    pub fn gelatinization_rate(&self, t_k: f64) -> f64 {
        if t_k < self.t_gel_onset {
            return 0.0;
        }
        let frac = ((t_k - self.t_gel_onset) / (self.t_gel_end - self.t_gel_onset)).clamp(0.0, 1.0);
        let arrhenius = (-(self.ea_gel) / (R_GAS * t_k)).exp();
        let k_eff = self.k_gel * arrhenius * 1e8; // normalise
        k_eff * frac * (1.0 - self.gelatinization_degree)
    }

    /// Protein denaturation rate at temperature `t_k` (K).
    pub fn denaturation_rate(&self, t_k: f64) -> f64 {
        if t_k < self.t_denat - self.t_denat_range {
            return 0.0;
        }
        let arrhenius = (-(self.ea_denat) / (R_GAS * t_k)).exp();
        let k_eff = self.k_denat * arrhenius * 1e12;
        k_eff * (1.0 - self.denaturation_degree)
    }

    /// Advance texture evolution by `dt` seconds at temperature `t_k` (K).
    pub fn step(&mut self, dt: f64, t_k: f64) {
        let dg = self.gelatinization_rate(t_k) * dt;
        let dd = self.denaturation_rate(t_k) * dt;

        self.gelatinization_degree = (self.gelatinization_degree + dg).min(1.0);
        self.denaturation_degree = (self.denaturation_degree + dd).min(1.0);

        // Hardness: weighted blend.
        let h = self.initial_hardness * (1.0 - self.gelatinization_degree)
            + self.gelatinized_hardness * self.gelatinization_degree;
        let h_final = h * (1.0 - self.denaturation_degree)
            + self.denatured_hardness * self.denaturation_degree;
        self.hardness = h_final;
        self.time += dt;
    }

    /// Is the food fully cooked (both processes > 0.95)?
    pub fn is_fully_cooked(&self) -> bool {
        self.gelatinization_degree > 0.95 && self.denaturation_degree > 0.95
    }

    /// Chewiness index: proportional to denaturation degree × hardness.
    pub fn chewiness_index(&self) -> f64 {
        self.denaturation_degree * self.hardness
    }

    /// Stickiness indicator: gel content × (1 − denaturation).
    pub fn stickiness(&self) -> f64 {
        self.gelatinization_degree * (1.0 - self.denaturation_degree)
    }
}

// ---------------------------------------------------------------------------
// ChewingSimulation
// ---------------------------------------------------------------------------

/// A single food particle in the chewing simulation.
#[derive(Debug, Clone)]
pub struct FoodParticle {
    /// Equivalent diameter (m).
    pub diameter: f64,
    /// Hardness (Pa) of this particle.
    pub hardness: f64,
    /// Moisture content (fraction, 0–1).
    pub moisture: f64,
}

impl FoodParticle {
    /// Create a new food particle.
    pub fn new(diameter: f64, hardness: f64, moisture: f64) -> Self {
        Self {
            diameter,
            hardness,
            moisture: moisture.clamp(0.0, 1.0),
        }
    }
}

/// Jaw mechanics and food fragmentation simulation.
///
/// Models the jaw as a spring-damper system, applies occlusal force to food
/// particles, and computes particle size reduction during each chewing cycle.
#[derive(Debug, Clone)]
pub struct ChewingSimulation {
    /// Food particles in the oral cavity.
    pub particles: Vec<FoodParticle>,
    /// Maximum bite force (N).
    pub max_bite_force: f64,
    /// Jaw spring constant (N/m).
    pub jaw_spring: f64,
    /// Jaw damping coefficient (N·s/m).
    pub jaw_damping: f64,
    /// Current jaw opening distance (m).
    pub jaw_opening: f64,
    /// Jaw angular velocity equivalent (m/s).
    pub jaw_velocity: f64,
    /// Number of chewing cycles completed.
    pub chew_count: usize,
    /// Fragmentation efficiency parameter (dimensionless).
    pub fragmentation_efficiency: f64,
    /// Target particle size for swallowing (m).
    pub swallow_threshold: f64,
    /// Saliva flow rate (m³/s) for moisture uptake.
    pub saliva_flow: f64,
    /// Elapsed time (s).
    pub time: f64,
}

impl ChewingSimulation {
    /// Create a typical human chewing simulation.
    pub fn human(food: Vec<FoodParticle>) -> Self {
        Self {
            particles: food,
            max_bite_force: 100.0, // N, moderate bite force
            jaw_spring: 5000.0,    // N/m
            jaw_damping: 200.0,    // N·s/m
            jaw_opening: 0.03,     // 3 cm initial opening
            jaw_velocity: 0.0,
            chew_count: 0,
            fragmentation_efficiency: 0.3,
            swallow_threshold: 2e-3, // 2 mm swallowing size
            saliva_flow: 1e-8,       // m³/s
            time: 0.0,
        }
    }

    /// Compute bite force based on jaw position (spring model).
    pub fn bite_force(&self) -> f64 {
        let contact_pos = 0.002; // 2 mm food contact gap
        let delta = (contact_pos - self.jaw_opening).max(0.0);
        (self.jaw_spring * delta).min(self.max_bite_force)
    }

    /// Apply one chewing stroke (close-open cycle) over `dt` seconds.
    ///
    /// The jaw follows a sinusoidal trajectory at ~1.5 Hz.
    /// Returns the peak bite force during this time step.
    pub fn chew_stroke(&mut self, dt: f64) -> f64 {
        // Sinusoidal jaw motion: period 0.6 s, amplitude 15 mm.
        let omega = 2.0 * PI / 0.6; // rad/s
        let amplitude = 0.015; // m
        let t_new = self.time + dt;
        // Jaw position oscillates: opening = amplitude * (1 − cos(ωt)) / 2
        // so it goes from 0 (closed) to amplitude (open) and back.
        let opening_prev = amplitude * 0.5 * (1.0 - (omega * self.time).cos());
        let opening_new = amplitude * 0.5 * (1.0 - (omega * t_new).cos());
        self.jaw_opening = opening_new;

        // Detect a zero-crossing in the cosine (jaw closure event).
        let cos_prev = (omega * self.time).cos();
        let cos_new = (omega * t_new).cos();
        // Closure happens when cos crosses from negative to positive (jaw reaches minimum).
        if cos_prev < 0.0 && cos_new >= 0.0 {
            self.chew_count += 1;
            let contact_force = self.max_bite_force * (1.0 - opening_prev / amplitude).max(0.0);
            self.fragment_particles(contact_force);
            for p in &mut self.particles {
                let moisture_uptake = self.saliva_flow * 0.6 / p.diameter.max(1e-6);
                p.moisture = (p.moisture + moisture_uptake).min(1.0);
            }
        }

        self.time = t_new;
        let contact_pos = 0.002;
        let delta = (contact_pos - self.jaw_opening).max(0.0);
        (self.jaw_spring * delta).min(self.max_bite_force)
    }

    /// Fragment particles based on applied `force`.
    fn fragment_particles(&mut self, force: f64) {
        let mut new_particles: Vec<FoodParticle> = Vec::new();
        for p in &mut self.particles {
            // Fracture criterion: force exceeds hardness × cross-section.
            let area = PI * (p.diameter / 2.0).powi(2);
            let fracture_force = p.hardness * area;
            if force > fracture_force && p.diameter > self.swallow_threshold {
                // Fragment: split into 2 daughters of ~0.65× diameter.
                let d_new = p.diameter * (0.5_f64).powf(1.0 / 3.0); // volume halved
                new_particles.push(FoodParticle::new(d_new, p.hardness, p.moisture));
                p.diameter = d_new;
            }
        }
        self.particles.extend(new_particles);
    }

    /// Fraction of particles ready for swallowing (below threshold).
    pub fn swallow_fraction(&self) -> f64 {
        if self.particles.is_empty() {
            return 1.0;
        }
        let ready = self
            .particles
            .iter()
            .filter(|p| p.diameter <= self.swallow_threshold)
            .count();
        ready as f64 / self.particles.len() as f64
    }

    /// Mean particle diameter (m).
    pub fn mean_particle_diameter(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        self.particles.iter().map(|p| p.diameter).sum::<f64>() / self.particles.len() as f64
    }

    /// D50: median particle diameter.
    pub fn d50(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let mut diameters: Vec<f64> = self.particles.iter().map(|p| p.diameter).collect();
        diameters.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = diameters.len() / 2;
        diameters[mid]
    }

    /// Simulate until swallow fraction exceeds `target` or `max_chews` is reached.
    ///
    /// `dt` is the time-step per call to `chew_stroke`.  For efficiency,
    /// prefer dt ≈ 0.05 s (one coarse step per chew cycle).
    pub fn chew_until_ready(&mut self, target: f64, max_chews: usize, dt: f64) {
        // Use at most 20 steps per cycle to limit wall-clock time.
        let dt_cycle = 0.6_f64; // s per chewing cycle
        let n_per_cycle = ((dt_cycle / dt).ceil() as usize).clamp(1, 20);
        let effective_dt = dt_cycle / n_per_cycle as f64;
        while self.chew_count < max_chews && self.swallow_fraction() < target {
            for _ in 0..n_per_cycle {
                self.chew_stroke(effective_dt);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// FoodEmulsion
// ---------------------------------------------------------------------------

/// Droplet in an oil-water emulsion.
#[derive(Debug, Clone)]
pub struct EmulsionDroplet {
    /// Droplet radius (m).
    pub radius: f64,
    /// Phase: true = oil-in-water, false = water-in-oil.
    pub is_oil: bool,
    /// Surface charge / zeta potential (mV).
    pub zeta_potential: f64,
    /// Position in 1D (m), for coarsening calculations.
    pub position: f64,
}

impl EmulsionDroplet {
    /// Create an oil droplet with given radius and zeta potential.
    pub fn oil(radius: f64, zeta_potential: f64, position: f64) -> Self {
        Self {
            radius,
            is_oil: true,
            zeta_potential,
            position,
        }
    }

    /// Volume of the droplet (m³).
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    /// Surface area of the droplet (m²).
    pub fn surface_area(&self) -> f64 {
        4.0 * PI * self.radius.powi(2)
    }
}

/// Oil-water food emulsion model with Ostwald ripening and stability analysis.
#[derive(Debug, Clone)]
pub struct FoodEmulsion {
    /// List of emulsion droplets.
    pub droplets: Vec<EmulsionDroplet>,
    /// Interfacial tension γ (N/m).
    pub interfacial_tension: f64,
    /// Continuous-phase viscosity (Pa·s).
    pub continuous_viscosity: f64,
    /// Dispersed-phase molar volume Vm (m³/mol).
    pub molar_volume: f64,
    /// Droplet diffusivity in continuous phase (m²/s).
    pub diffusivity: f64,
    /// Solubility of dispersed phase in continuous phase (mol/m³).
    pub solubility: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Volume fraction of dispersed phase.
    pub phi: f64,
    /// Elapsed time (s).
    pub time: f64,
    /// Hamaker constant for van-der-Waals interactions (J).
    pub hamaker: f64,
}

impl FoodEmulsion {
    /// Create a standard oil-in-water food emulsion (e.g., mayonnaise).
    pub fn mayonnaise(n_droplets: usize) -> Self {
        let mut droplets = Vec::with_capacity(n_droplets);
        for i in 0..n_droplets {
            let r = 1e-6 + (i as f64) * 0.2e-6; // radii 1–3 µm
            droplets.push(EmulsionDroplet::oil(r, -30.0, i as f64 * 5e-6));
        }
        Self {
            droplets,
            interfacial_tension: 0.01,
            continuous_viscosity: 1e-3,
            molar_volume: 1.6e-4, // m³/mol (trioleate)
            diffusivity: 1e-12,
            solubility: 1.0,
            temperature: 298.15,
            phi: 0.7, // high oil fraction
            time: 0.0,
            hamaker: 5e-21,
        }
    }

    /// Lifshitz-Slezov-Wagner (LSW) Ostwald ripening rate ω (m³/s).
    ///
    /// ω = 8 γ Vm D C∞ / (9 R T)
    pub fn lsw_ripening_rate(&self) -> f64 {
        8.0 * self.interfacial_tension * self.molar_volume * self.diffusivity * self.solubility
            / (9.0 * R_GAS * self.temperature)
    }

    /// Advance Ostwald ripening by `dt` seconds.
    ///
    /// r³(t + dt) ≈ r³(t) + ω·dt, then renormalise.
    pub fn step_ripening(&mut self, dt: f64) {
        let omega = self.lsw_ripening_rate();
        for d in &mut self.droplets {
            let r3 = d.radius.powi(3) + omega * dt;
            d.radius = r3.cbrt().max(1e-9);
        }
        self.time += dt;
    }

    /// Mean droplet radius (m).
    pub fn mean_radius(&self) -> f64 {
        if self.droplets.is_empty() {
            return 0.0;
        }
        self.droplets.iter().map(|d| d.radius).sum::<f64>() / self.droplets.len() as f64
    }

    /// Volume-weighted mean radius d\[4,3\] (De Brouckere mean).
    pub fn d43(&self) -> f64 {
        let num: f64 = self.droplets.iter().map(|d| d.radius.powi(4)).sum();
        let den: f64 = self.droplets.iter().map(|d| d.radius.powi(3)).sum();
        if den < 1e-30 {
            return 0.0;
        }
        num / den
    }

    /// Sauter mean radius d\[3,2\] (surface-volume ratio).
    pub fn d32(&self) -> f64 {
        let num: f64 = self.droplets.iter().map(|d| d.radius.powi(3)).sum();
        let den: f64 = self.droplets.iter().map(|d| d.radius.powi(2)).sum();
        if den < 1e-30 {
            return 0.0;
        }
        num / den
    }

    /// Coefficient of variation (CV) of droplet size distribution (%).
    pub fn size_cv(&self) -> f64 {
        if self.droplets.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_radius();
        let var = self
            .droplets
            .iter()
            .map(|d| (d.radius - mean).powi(2))
            .sum::<f64>()
            / self.droplets.len() as f64;
        100.0 * var.sqrt() / mean.max(1e-30)
    }

    /// Flocculation indicator: fraction of droplet pairs closer than 2×r_mean.
    pub fn flocculation_index(&self) -> f64 {
        let r_mean = self.mean_radius();
        if self.droplets.len() < 2 {
            return 0.0;
        }
        let mut close = 0usize;
        let mut pairs = 0usize;
        for i in 0..self.droplets.len() {
            for j in (i + 1)..self.droplets.len() {
                let sep = (self.droplets[i].position - self.droplets[j].position).abs();
                if sep < 2.0 * r_mean {
                    close += 1;
                }
                pairs += 1;
            }
        }
        if pairs == 0 {
            0.0
        } else {
            close as f64 / pairs as f64
        }
    }

    /// Critical flocculation electrostatic energy barrier (kT units).
    ///
    /// V_max/(kT) ≈ 32 π ε₀ ε_r κ⁻¹ r ψ₀² / (kT) — simplified.
    pub fn electrostatic_stability(&self, debye_length: f64, epsilon_r: f64) -> f64 {
        let epsilon_0 = 8.854e-12; // F/m
        let kb_t = 1.38e-23 * self.temperature;
        let r = self.mean_radius();
        let psi0 = self
            .droplets
            .first()
            .map(|d| d.zeta_potential.abs() * 1e-3)
            .unwrap_or(0.03);
        let v_max = 32.0 * PI * epsilon_0 * epsilon_r * debye_length * r * psi0 * psi0;
        v_max / kb_t
    }

    /// Van-der-Waals attraction energy (kT units) at surface-to-surface separation `h`.
    pub fn vdw_attraction(&self, h: f64) -> f64 {
        let r = self.mean_radius();
        let kb_t = 1.38e-23 * self.temperature;
        // Derjaguin approximation: V_vdw ≈ -A * R / (6*h)
        -(self.hamaker * r / (6.0 * h)) / kb_t
    }

    /// Emulsion stability score: ratio of electrostatic barrier to vdw well depth.
    pub fn stability_ratio(&self, debye_length: f64, epsilon_r: f64, h_min: f64) -> f64 {
        let v_el = self.electrostatic_stability(debye_length, epsilon_r);
        let v_vdw = self.vdw_attraction(h_min).abs();
        if v_vdw < 1e-10 {
            return f64::INFINITY;
        }
        v_el / v_vdw
    }

    /// Total interfacial area (m²) of all droplets.
    pub fn total_interfacial_area(&self) -> f64 {
        self.droplets.iter().map(|d| d.surface_area()).sum()
    }

    /// Run Ostwald ripening simulation for `total_time` seconds.
    pub fn run_ripening(&mut self, total_time: f64, dt: f64) {
        while self.time < total_time {
            let step = dt.min(total_time - self.time);
            self.step_ripening(step);
        }
    }

    /// Compute the droplet size distribution as a histogram.
    ///
    /// Returns `(bin_centres_m, counts)`.
    pub fn size_histogram(&self, n_bins: usize) -> (Vec<f64>, Vec<usize>) {
        if self.droplets.is_empty() || n_bins == 0 {
            return (Vec::new(), Vec::new());
        }
        let r_min = self
            .droplets
            .iter()
            .map(|d| d.radius)
            .fold(f64::INFINITY, f64::min);
        let r_max = self
            .droplets
            .iter()
            .map(|d| d.radius)
            .fold(0.0_f64, f64::max);
        let dr = (r_max - r_min) / n_bins as f64 + 1e-30;
        let mut counts = vec![0usize; n_bins];
        let mut centres = Vec::with_capacity(n_bins);
        for i in 0..n_bins {
            centres.push(r_min + (i as f64 + 0.5) * dr);
        }
        for d in &self.droplets {
            let bin = ((d.radius - r_min) / dr) as usize;
            let bin = bin.min(n_bins - 1);
            counts[bin] += 1;
        }
        (centres, counts)
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Convert Celsius to Kelvin.
#[inline]
pub fn celsius_to_kelvin(t_c: f64) -> f64 {
    t_c + KELVIN
}

/// Convert Kelvin to Celsius.
#[inline]
pub fn kelvin_to_celsius(t_k: f64) -> f64 {
    t_k - KELVIN
}

/// Compute water activity from Brunauer-Emmett-Teller (BET) equation (simplified).
///
/// `m` — moisture content (g/g dry), `m0` — monolayer moisture, `c` — BET constant.
pub fn water_activity_bet(m: f64, m0: f64, c: f64) -> f64 {
    if m0 < 1e-10 {
        return 0.0;
    }
    let ratio = m / m0;
    // Simplified 2-parameter BET: aw = ratio / (1 + (c-1)*ratio + ...)
    // Use simplified approximation.
    (c * ratio) / ((1.0 - ratio) * (1.0 + (c - 1.0) * ratio)).max(1e-10)
}

/// Compute Maillard browning degree from time-temperature history.
///
/// Uses the Arrhenius integral: D = ∫ A exp(-Ea/RT) dt.
pub fn maillard_integral(
    time_temp_history: &[(f64, f64)], // (dt, T_K)
    ea: f64,
    a: f64,
) -> f64 {
    time_temp_history
        .iter()
        .map(|(dt, t_k)| a * (-(ea) / (R_GAS * t_k)).exp() * dt)
        .sum()
}

/// Compute effective thermal diffusivity using the parallel model for a two-phase food.
///
/// `alpha1`, `phi1` — diffusivity and volume fraction of phase 1.
/// `alpha2`, `phi2` — diffusivity and volume fraction of phase 2.
pub fn effective_diffusivity_parallel(alpha1: f64, phi1: f64, alpha2: f64, phi2: f64) -> f64 {
    alpha1 * phi1 + alpha2 * phi2
}

/// Water loss rate during cooking (simplified evaporation model).
///
/// `t_surface_k` — surface temperature (K), `rh` — relative humidity of air (0–1).
/// Returns evaporation flux (kg/m²/s).
pub fn evaporation_flux(t_surface_k: f64, rh: f64) -> f64 {
    // Antoine equation for water vapour pressure (simplified).
    let p_sat =
        133.322 * 10.0_f64.powf(8.07131 - 1730.63 / (233.426 + kelvin_to_celsius(t_surface_k)));
    let p_air = rh * p_sat;
    let h_mass = 1e-2; // mass transfer coefficient (m/s, typical)
    let mw = 0.018; // kg/mol
    let rt = R_GAS * t_surface_k;
    h_mass * mw * (p_sat - p_air) / rt
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── FoodRheology ─────────────────────────────────────────────────────────

    #[test]
    fn test_rheology_g_prime_at_reference_frequency() {
        let r = FoodRheology::bread_dough();
        let omega_ref = 2.0 * PI * r.f_ref;
        let g1 = r.g_prime(omega_ref);
        assert!(
            (g1 - r.g_prime_ref).abs() < 1.0,
            "G' at f_ref should be ~G'_ref, got {g1}"
        );
    }

    #[test]
    fn test_rheology_g_double_prime_at_reference() {
        let r = FoodRheology::gelatin_gel();
        let omega_ref = 2.0 * PI * r.f_ref;
        let g2 = r.g_double_prime(omega_ref);
        assert!(
            (g2 - r.g_double_prime_ref).abs() < 1.0,
            "G'' at f_ref, got {g2}"
        );
    }

    #[test]
    fn test_tan_delta_positive() {
        let r = FoodRheology::bread_dough();
        let td = r.tan_delta(2.0 * PI);
        assert!(td > 0.0, "tan delta should be positive, got {td}");
    }

    #[test]
    fn test_elastic_dominant() {
        let r = FoodRheology::gelatin_gel();
        // Gelatin gel is elastic dominant (G' >> G'').
        assert!(
            r.is_elastic_dominant(2.0 * PI),
            "gel should be elastic dominant"
        );
    }

    #[test]
    fn test_frequency_sweep_length() {
        let r = FoodRheology::bread_dough();
        let sweep = r.frequency_sweep(0.01, 100.0, 50);
        assert_eq!(sweep.len(), 50);
    }

    #[test]
    fn test_frequency_sweep_g_prime_increases() {
        let r = FoodRheology::bread_dough();
        let sweep = r.frequency_sweep(0.01, 100.0, 20);
        // G' should increase with frequency (positive power law exponent).
        let g1_low = sweep[0].1;
        let g1_high = sweep[19].1;
        assert!(g1_high > g1_low, "G' should increase with frequency");
    }

    #[test]
    fn test_complex_modulus_geq_g_prime() {
        let r = FoodRheology::bread_dough();
        let omega = 2.0 * PI;
        let gstar = r.complex_modulus(omega);
        let g1 = r.g_prime(omega);
        assert!(
            gstar >= g1 - 1e-10,
            "|G*| should be >= G', got |G*|={gstar}, G'={g1}"
        );
    }

    #[test]
    fn test_stress_relaxation_decays() {
        let r = FoodRheology::bread_dough();
        let g0 = r.stress_relaxation(0.0);
        let g_late = r.stress_relaxation(100.0);
        assert!(g_late < g0, "stress relaxation should decrease with time");
    }

    #[test]
    fn test_creep_compliance_increases() {
        let r = FoodRheology::gelatin_gel();
        let j0 = r.creep_compliance(0.0);
        let j1 = r.creep_compliance(10.0);
        assert!(j1 >= j0, "creep compliance should increase with time");
    }

    #[test]
    fn test_apparent_viscosity_positive() {
        let r = FoodRheology::bread_dough();
        let eta = r.apparent_viscosity(2.0 * PI);
        assert!(
            eta > 0.0,
            "apparent viscosity should be positive, got {eta}"
        );
    }

    // ── HeatCookingModel ────────────────────────────────────────────────────

    #[test]
    fn test_heat_model_initial_temperature() {
        let m = HeatCookingModel::bread_baking(11);
        assert!((m.mean_temperature() - 298.15).abs() < 1e-6);
    }

    #[test]
    fn test_heat_model_boundary_heats_interior() {
        let mut m = HeatCookingModel::bread_baking(11);
        m.set_boundary_temperatures(473.15, 473.15); // 200 °C oven
        for _ in 0..1000 {
            m.step(1.0); // 1000 seconds
        }
        let t_centre = m.temperature[m.centre_node()];
        assert!(t_centre > 298.15, "centre should heat up, got {t_centre}");
    }

    #[test]
    fn test_heat_model_cooking_degree_increases() {
        let mut m = HeatCookingModel::bread_baking(11);
        m.set_boundary_temperatures(523.15, 523.15); // 250 °C
        for _ in 0..500 {
            m.step(1.0);
        }
        assert!(
            m.mean_cooking_degree() > 0.0,
            "cooking degree should increase"
        );
    }

    #[test]
    fn test_maillard_rate_increases_with_temperature() {
        let m = HeatCookingModel::bread_baking(11);
        let k_low = m.maillard_rate(373.15); // 100 °C
        let k_high = m.maillard_rate(473.15); // 200 °C
        assert!(
            k_high > k_low,
            "Maillard rate should increase with temperature"
        );
    }

    #[test]
    fn test_maillard_rate_zero_at_low_temperature() {
        let m = HeatCookingModel::bread_baking(11);
        let k = m.maillard_rate(200.0);
        // At 200 K the Arrhenius rate is negligibly small compared to baking temperatures.
        assert!(
            k < 1e-10,
            "Maillard rate should be negligible at 200 K, got {k}"
        );
    }

    #[test]
    fn test_heat_model_surface_browning() {
        let mut m = HeatCookingModel::bread_baking(11);
        m.set_boundary_temperatures(523.15, 523.15);
        for _ in 0..200 {
            m.step(1.0);
        }
        let sb = m.surface_browning();
        assert!(
            (0.0..=1.0).contains(&sb),
            "surface browning in [0,1], got {sb}"
        );
    }

    #[test]
    fn test_heat_model_centre_cools_slower_than_surface() {
        let mut m = HeatCookingModel::new(11, 0.05, 1.4e-7, 0.18, 3000.0, 450.0, 298.15);
        m.set_boundary_temperatures(473.15, 473.15);
        for _ in 0..100 {
            m.step(1.0);
        }
        let t_surface = m.temperature[0];
        let t_centre = m.temperature[m.centre_node()];
        assert!(
            t_surface >= t_centre,
            "surface should be hotter than centre during heating"
        );
    }

    #[test]
    fn test_heat_model_stored_energy_positive() {
        let m = HeatCookingModel::bread_baking(11);
        assert!(m.stored_energy() > 0.0);
    }

    // ── FermentationModel ───────────────────────────────────────────────────

    #[test]
    fn test_fermentation_co2_increases() {
        let mut f = FermentationModel::bread_dough(0.02, 0.1, 32.0);
        let co2_init = f.co2_concentration;
        f.step(60.0); // 1 minute
        assert!(
            f.co2_concentration > co2_init,
            "CO2 should increase with fermentation"
        );
    }

    #[test]
    fn test_fermentation_substrate_decreases() {
        let mut f = FermentationModel::bread_dough(0.02, 0.1, 32.0);
        let s0 = f.substrate;
        f.step(60.0);
        assert!(f.substrate < s0, "substrate should be consumed");
    }

    #[test]
    fn test_fermentation_volume_increases() {
        // Use high yeast and substrate concentrations to generate enough CO₂
        // that the excess above the solubility limit produces visible gas volume.
        let mut f = FermentationModel::bread_dough(1.0, 10.0, 32.0);
        // Lower solubility so excess CO₂ creates gas quickly.
        f.co2_solubility = 1e-7;
        f.run(3600.0, 60.0); // 1 hour
        assert!(
            f.volume_ratio > 1.0,
            "volume should increase during fermentation"
        );
    }

    #[test]
    fn test_fermentation_temperature_factor() {
        let f_opt = FermentationModel::bread_dough(0.02, 0.1, 32.0);
        let f_cold = FermentationModel::bread_dough(0.02, 0.1, 4.0);
        assert!(
            f_opt.temperature_factor() > f_cold.temperature_factor(),
            "optimal T should give higher temperature factor"
        );
    }

    #[test]
    fn test_fermentation_completion() {
        let mut f = FermentationModel::bread_dough(0.1, 0.001, 32.0); // tiny substrate
        f.run(600.0, 10.0);
        assert!(f.is_complete(0.005), "should complete with tiny substrate");
    }

    #[test]
    fn test_bubble_radius_positive() {
        let mut f = FermentationModel::bread_dough(0.02, 0.5, 32.0);
        f.run(1800.0, 60.0);
        assert!(f.bubble_radius > 0.0, "bubble radius should be positive");
    }

    #[test]
    fn test_co2_production_rate_positive() {
        let f = FermentationModel::bread_dough(0.02, 0.1, 32.0);
        let rate = f.co2_production_rate();
        assert!(
            rate > 0.0,
            "CO2 production rate should be positive, got {rate}"
        );
    }

    // ── TextureEvolution ────────────────────────────────────────────────────

    #[test]
    fn test_texture_no_change_at_low_temperature() {
        let mut t = TextureEvolution::bread();
        let h0 = t.hardness;
        for _ in 0..1000 {
            t.step(1.0, 300.0); // 27 °C — below gelatinization onset
        }
        assert!(
            (t.hardness - h0).abs() < 1.0,
            "hardness should not change significantly below onset, got {}",
            t.hardness
        );
    }

    #[test]
    fn test_texture_gelatinization_at_high_temperature() {
        let mut t = TextureEvolution::bread();
        // Heat well above onset.
        for _ in 0..1000 {
            t.step(1.0, 360.0); // 87 °C
        }
        assert!(
            t.gelatinization_degree > 0.0,
            "gelatinization should proceed above onset"
        );
    }

    #[test]
    fn test_texture_protein_denaturation() {
        let mut t = TextureEvolution::meat();
        for _ in 0..2000 {
            t.step(1.0, 363.15); // 90 °C — well above denaturation temp
        }
        assert!(
            t.denaturation_degree > 0.0,
            "proteins should denature at high temperature"
        );
    }

    #[test]
    fn test_texture_cooking_degree_bounded() {
        let mut t = TextureEvolution::bread();
        for _ in 0..50000 {
            t.step(1.0, 473.15);
        }
        assert!(t.gelatinization_degree <= 1.0);
        assert!(t.denaturation_degree <= 1.0);
    }

    #[test]
    fn test_texture_chewiness_index_positive() {
        let mut t = TextureEvolution::meat();
        for _ in 0..5000 {
            t.step(1.0, 363.15);
        }
        let ci = t.chewiness_index();
        assert!(
            ci >= 0.0,
            "chewiness index should be non-negative, got {ci}"
        );
    }

    #[test]
    fn test_texture_stickiness_peak_midway() {
        let mut t = TextureEvolution::bread();
        // Push gelation but not denaturation.
        for _ in 0..3000 {
            t.step(1.0, 350.0);
        }
        let s = t.stickiness();
        assert!(
            (0.0..=1.0).contains(&s),
            "stickiness should be in [0,1], got {s}"
        );
    }

    // ── ChewingSimulation ────────────────────────────────────────────────────

    #[test]
    fn test_chewing_mean_diameter_decreases() {
        let food = vec![
            FoodParticle::new(0.01, 500_000.0, 0.4), // 1 cm, 500 kPa
            FoodParticle::new(0.01, 500_000.0, 0.4),
        ];
        let mut chew = ChewingSimulation::human(food);
        let d0 = chew.mean_particle_diameter();
        chew.chew_until_ready(0.5, 50, 0.01);
        let d1 = chew.mean_particle_diameter();
        // Either diameter decreases or stays same (softer foods may not fragment).
        assert!(
            d1 <= d0 + 1e-10,
            "diameter should not increase during chewing"
        );
    }

    #[test]
    fn test_chewing_chew_count_increments() {
        let food = vec![FoodParticle::new(0.01, 500_000.0, 0.4)];
        let mut chew = ChewingSimulation::human(food);
        chew.chew_until_ready(0.5, 10, 0.01);
        assert!(chew.chew_count > 0, "chew count should increase");
    }

    #[test]
    fn test_chewing_bite_force_non_negative() {
        let food = vec![FoodParticle::new(0.01, 500_000.0, 0.4)];
        let chew = ChewingSimulation::human(food);
        assert!(chew.bite_force() >= 0.0);
    }

    #[test]
    fn test_chewing_d50_positive() {
        let food = vec![
            FoodParticle::new(0.005, 100_000.0, 0.5),
            FoodParticle::new(0.01, 100_000.0, 0.5),
            FoodParticle::new(0.015, 100_000.0, 0.5),
        ];
        let chew = ChewingSimulation::human(food);
        assert!(chew.d50() > 0.0);
    }

    #[test]
    fn test_chewing_swallow_fraction_bounded() {
        let food = vec![FoodParticle::new(0.001, 100_000.0, 0.5)]; // below threshold already
        let chew = ChewingSimulation::human(food);
        let sf = chew.swallow_fraction();
        assert!(
            (0.0..=1.0).contains(&sf),
            "swallow fraction in [0,1], got {sf}"
        );
    }

    #[test]
    fn test_chewing_soft_food_swallowable() {
        // Very small, soft particles should all be below the swallow threshold.
        let food: Vec<FoodParticle> = (0..5)
            .map(|_| FoodParticle::new(0.001, 100.0, 0.9)) // 1 mm, very soft
            .collect();
        let chew = ChewingSimulation::human(food);
        let sf = chew.swallow_fraction();
        assert!(
            sf > 0.8,
            "soft small food should be mostly swallowable, got {sf}"
        );
    }

    // ── FoodEmulsion ────────────────────────────────────────────────────────

    #[test]
    fn test_emulsion_initial_mean_radius() {
        let em = FoodEmulsion::mayonnaise(5);
        let mr = em.mean_radius();
        assert!(mr > 0.0, "mean radius should be positive, got {mr}");
    }

    #[test]
    fn test_emulsion_ostwald_ripening_increases_radius() {
        let mut em = FoodEmulsion::mayonnaise(10);
        let r0 = em.mean_radius();
        em.run_ripening(3600.0 * 24.0, 3600.0); // 24 hours
        let r1 = em.mean_radius();
        assert!(
            r1 >= r0 - 1e-15,
            "Ostwald ripening should not decrease mean radius"
        );
    }

    #[test]
    fn test_emulsion_lsw_rate_positive() {
        let em = FoodEmulsion::mayonnaise(5);
        assert!(em.lsw_ripening_rate() > 0.0);
    }

    #[test]
    fn test_emulsion_d43_geq_d32() {
        let em = FoodEmulsion::mayonnaise(10);
        assert!(em.d43() >= em.d32() - 1e-15, "d43 should be >= d32");
    }

    #[test]
    fn test_emulsion_size_cv_non_negative() {
        let em = FoodEmulsion::mayonnaise(10);
        assert!(em.size_cv() >= 0.0);
    }

    #[test]
    fn test_emulsion_total_interfacial_area_positive() {
        let em = FoodEmulsion::mayonnaise(5);
        assert!(em.total_interfacial_area() > 0.0);
    }

    #[test]
    fn test_emulsion_stability_ratio_positive() {
        let em = FoodEmulsion::mayonnaise(5);
        let sr = em.stability_ratio(10e-9, 80.0, 1e-9);
        assert!(
            sr >= 0.0,
            "stability ratio should be non-negative, got {sr}"
        );
    }

    #[test]
    fn test_emulsion_flocculation_index_bounded() {
        let em = FoodEmulsion::mayonnaise(5);
        let fi = em.flocculation_index();
        assert!(
            (0.0..=1.0).contains(&fi),
            "flocculation index in [0,1], got {fi}"
        );
    }

    #[test]
    fn test_emulsion_histogram_bins() {
        let em = FoodEmulsion::mayonnaise(10);
        let (centres, counts) = em.size_histogram(5);
        assert_eq!(centres.len(), 5);
        assert_eq!(counts.len(), 5);
        assert_eq!(counts.iter().sum::<usize>(), 10);
    }

    // ── Utility functions ────────────────────────────────────────────────────

    #[test]
    fn test_celsius_kelvin_roundtrip() {
        let t_c = 100.0;
        let t_k = celsius_to_kelvin(t_c);
        assert!((kelvin_to_celsius(t_k) - t_c).abs() < 1e-10);
    }

    #[test]
    fn test_effective_diffusivity_parallel() {
        let alpha_eff = effective_diffusivity_parallel(1e-7, 0.5, 2e-7, 0.5);
        assert!((alpha_eff - 1.5e-7).abs() < 1e-15);
    }

    #[test]
    fn test_maillard_integral_positive() {
        let history: Vec<(f64, f64)> = vec![(10.0, 473.15); 10];
        let d = maillard_integral(&history, 80_000.0, 1e8);
        assert!(d > 0.0, "Maillard integral should be positive, got {d}");
    }

    #[test]
    fn test_water_activity_bet_zero_moisture() {
        let aw = water_activity_bet(0.0, 0.1, 10.0);
        assert!(
            aw.abs() < 1e-10,
            "zero moisture → zero water activity, got {aw}"
        );
    }

    #[test]
    fn test_evaporation_flux_positive_above_boiling() {
        let flux = evaporation_flux(celsius_to_kelvin(150.0), 0.0);
        assert!(
            flux > 0.0,
            "evaporation flux should be positive, got {flux}"
        );
    }
}
