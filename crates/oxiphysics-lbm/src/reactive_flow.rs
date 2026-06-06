// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reactive flow and combustion modeling for LBM.
//!
//! Implements one-step chemistry, species transport, flame properties,
//! ignition models, and a 2D reactive flow grid.

/// Universal gas constant (J/mol/K).
const R_U: f64 = 8.314;

// ─────────────────────────────────────────────────────────────
// Species
// ─────────────────────────────────────────────────────────────

/// Thermodynamic and transport properties of a chemical species.
pub struct Species {
    /// Species name (e.g. "CH4").
    pub name: String,
    /// Molecular weight (g/mol).
    pub molecular_weight: f64,
    /// Specific heat at constant pressure (J/kg/K).
    pub specific_heat_cp: f64,
    /// Standard enthalpy of formation (J/kg).
    pub formation_enthalpy: f64,
    /// Lewis number Le = α/D (thermal-to-mass diffusivity ratio).
    pub lewis_number: f64,
}

impl Species {
    /// Create a new species.
    pub fn new(
        name: &str,
        molecular_weight: f64,
        specific_heat_cp: f64,
        formation_enthalpy: f64,
        lewis_number: f64,
    ) -> Self {
        Self {
            name: name.to_string(),
            molecular_weight,
            specific_heat_cp,
            formation_enthalpy,
            lewis_number,
        }
    }
}

// ─────────────────────────────────────────────────────────────
// ChemicalReaction (one-step Arrhenius)
// ─────────────────────────────────────────────────────────────

/// Simplified one-step chemical reaction with Arrhenius kinetics.
pub struct ChemicalReaction {
    /// Stoichiometric coefficient for fuel.
    pub fuel_stoich: f64,
    /// Stoichiometric coefficient for oxidizer.
    pub oxidizer_stoich: f64,
    /// Stoichiometric coefficient for products.
    pub product_stoich: f64,
    /// Activation energy Eₐ (J/mol).
    pub activation_energy: f64,
    /// Pre-exponential factor A.
    pub pre_exponential: f64,
    /// Heat release Q (J/kg_fuel).
    pub heat_release: f64,
}

impl ChemicalReaction {
    /// Create a new one-step reaction.
    pub fn new(
        fuel_stoich: f64,
        oxidizer_stoich: f64,
        product_stoich: f64,
        activation_energy: f64,
        pre_exponential: f64,
        heat_release: f64,
    ) -> Self {
        Self {
            fuel_stoich,
            oxidizer_stoich,
            product_stoich,
            activation_energy,
            pre_exponential,
            heat_release,
        }
    }

    /// Compute the reaction rate (kg/m³/s) via Arrhenius:
    ///
    /// `w = A · Y_fuel · Y_ox · exp(−Eₐ / (Rᵤ · T))`
    ///
    /// Returns 0 when temperature ≤ 0 to avoid numerical issues.
    pub fn reaction_rate(&self, y_fuel: f64, y_ox: f64, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let exponent = -self.activation_energy / (R_U * temperature);
        self.pre_exponential * y_fuel * y_ox * exponent.exp()
    }
}

// ─────────────────────────────────────────────────────────────
// CombustionCell
// ─────────────────────────────────────────────────────────────

/// State of a single computational cell in the reactive flow.
#[derive(Clone)]
pub struct CombustionCell {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Fuel mass fraction.
    pub y_fuel: f64,
    /// Oxidizer mass fraction.
    pub y_oxidizer: f64,
    /// Product mass fraction.
    pub y_products: f64,
    /// 2-D velocity \[ux, uy\] (m/s).
    pub velocity: [f64; 2],
}

impl CombustionCell {
    /// Create a cell at rest with uniform initial conditions.
    pub fn new(temperature: f64, pressure: f64, density: f64) -> Self {
        Self {
            temperature,
            pressure,
            density,
            y_fuel: 0.0,
            y_oxidizer: 0.0,
            y_products: 0.0,
            velocity: [0.0; 2],
        }
    }
}

// ─────────────────────────────────────────────────────────────
// FlameProperties
// ─────────────────────────────────────────────────────────────

/// Utility functions for laminar flame properties.
pub struct FlameProperties;

impl FlameProperties {
    /// Simplified laminar flame speed (m/s):
    ///
    /// `S_L = S_L0 · (T/T_ref)^α · (P/P_ref)^β`
    ///
    /// with S_L0 = 0.4 m/s (methane), α = 2.2, β = −0.5,
    /// T_ref = 298 K, P_ref = 101325 Pa.
    pub fn laminar_flame_speed(y_fuel: f64, temperature: f64, pressure: f64) -> f64 {
        let _ = y_fuel; // fuel fraction can be used for richer correlations
        const S_L0: f64 = 0.4;
        const ALPHA: f64 = 2.2;
        const BETA: f64 = -0.5;
        const T_REF: f64 = 298.0;
        const P_REF: f64 = 101_325.0;
        S_L0 * (temperature / T_REF).powf(ALPHA) * (pressure / P_REF).powf(BETA)
    }

    /// Adiabatic flame temperature (K):
    ///
    /// `T_ad = T_initial + Y_fuel · Q / cₚ`
    pub fn adiabatic_flame_temperature(
        y_fuel: f64,
        cp: f64,
        heat_release: f64,
        t_initial: f64,
    ) -> f64 {
        t_initial + y_fuel * heat_release / cp
    }

    /// Equivalence ratio:
    ///
    /// `φ = (Y_fuel / Y_oxidizer) / stoich_ratio`
    ///
    /// Returns 0 when the oxidizer fraction is zero.
    pub fn equivalence_ratio(y_fuel: f64, y_oxidizer: f64, stoich_ratio: f64) -> f64 {
        if y_oxidizer == 0.0 {
            return 0.0;
        }
        (y_fuel / y_oxidizer) / stoich_ratio
    }
}

// ─────────────────────────────────────────────────────────────
// ReactiveFlowGrid (2D)
// ─────────────────────────────────────────────────────────────

/// 2-D reactive flow grid coupling chemistry with diffusive transport.
pub struct ReactiveFlowGrid {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Flat cell storage (row-major: idx = x + y*nx).
    pub cells: Vec<CombustionCell>,
    /// Chemical reaction model.
    pub reaction: ChemicalReaction,
    /// Species mass diffusivity (m²/s).
    pub diffusivity: f64,
    /// Thermal diffusivity (m²/s).
    pub thermal_diffusivity: f64,
}

impl ReactiveFlowGrid {
    /// Create a new grid with default air-like initial conditions.
    pub fn new(nx: usize, ny: usize, reaction: ChemicalReaction) -> Self {
        let n = nx * ny;
        let default_cell = CombustionCell::new(300.0, 101_325.0, 1.2);
        Self {
            nx,
            ny,
            cells: vec![default_cell; n],
            reaction,
            diffusivity: 2e-5,
            thermal_diffusivity: 2e-5,
        }
    }

    /// Linear index for (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        x + y * self.nx
    }

    /// Initialize a premixed flame configuration:
    /// - `x < flame_x` → fresh mixture at `t_fresh`
    /// - `x >= flame_x` → hot products at `t_burned`
    pub fn initialize_premixed(
        &mut self,
        y_fuel: f64,
        y_ox: f64,
        t_fresh: f64,
        t_burned: f64,
        flame_x: usize,
    ) {
        for cy in 0..self.ny {
            for cx in 0..self.nx {
                let idx = self.idx(cx, cy);
                if cx < flame_x {
                    self.cells[idx].y_fuel = y_fuel;
                    self.cells[idx].y_oxidizer = y_ox;
                    self.cells[idx].y_products = 0.0;
                    self.cells[idx].temperature = t_fresh;
                } else {
                    self.cells[idx].y_fuel = 0.0;
                    self.cells[idx].y_oxidizer = 0.0;
                    self.cells[idx].y_products = 1.0;
                    self.cells[idx].temperature = t_burned;
                }
            }
        }
    }

    /// Advance chemistry for all cells by `dt` (s).
    ///
    /// Uses explicit forward Euler on the reaction source terms.
    pub fn step_chemistry(&mut self, dt: f64) {
        for cell in self.cells.iter_mut() {
            let w = cell.reaction_rate_for(&self.reaction);
            let d_fuel = -w * dt;
            let d_ox = -w * self.reaction.oxidizer_stoich * dt;
            let d_prod = w * self.reaction.product_stoich * dt;
            let d_temp = w * self.reaction.heat_release * dt;

            cell.y_fuel = (cell.y_fuel + d_fuel).max(0.0);
            cell.y_oxidizer = (cell.y_oxidizer + d_ox).max(0.0);
            cell.y_products = (cell.y_products + d_prod).min(1.0);
            cell.temperature += d_temp;
        }
    }

    /// Explicit finite-difference Laplacian diffusion of species (y_fuel, y_oxidizer).
    pub fn diffuse_species(&mut self, dt: f64, dx: f64) {
        let coeff = self.diffusivity * dt / (dx * dx);
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;

        let y_fuel_old: Vec<f64> = self.cells.iter().map(|c| c.y_fuel).collect();
        let y_ox_old: Vec<f64> = self.cells.iter().map(|c| c.y_oxidizer).collect();

        for idx in 0..n {
            let x = idx % nx;
            let y = idx / nx;

            let lap_fuel = laplacian_2d(&y_fuel_old, x, y, nx, ny);
            let lap_ox = laplacian_2d(&y_ox_old, x, y, nx, ny);

            self.cells[idx].y_fuel = (self.cells[idx].y_fuel + coeff * lap_fuel).max(0.0);
            self.cells[idx].y_oxidizer = (self.cells[idx].y_oxidizer + coeff * lap_ox).max(0.0);
        }
    }

    /// Explicit finite-difference Laplacian diffusion of temperature.
    pub fn diffuse_temperature(&mut self, dt: f64, dx: f64) {
        let coeff = self.thermal_diffusivity * dt / (dx * dx);
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;

        let temp_old: Vec<f64> = self.cells.iter().map(|c| c.temperature).collect();

        for idx in 0..n {
            let x = idx % nx;
            let y = idx / nx;
            let lap_t = laplacian_2d(&temp_old, x, y, nx, ny);
            self.cells[idx].temperature += coeff * lap_t;
        }
    }

    /// One full time step: chemistry then diffusion.
    pub fn step(&mut self, dt: f64, dx: f64) {
        self.step_chemistry(dt);
        self.diffuse_species(dt, dx);
        self.diffuse_temperature(dt, dx);
    }

    /// x-coordinate where y_fuel first drops to 50 % of the maximum y_fuel.
    ///
    /// Returns `None` if no fuel is present.
    pub fn flame_position(&self) -> Option<f64> {
        let max_fuel = self.cells.iter().map(|c| c.y_fuel).fold(0.0_f64, f64::max);
        if max_fuel <= 0.0 {
            return None;
        }
        let threshold = 0.5 * max_fuel;
        // Search column-averaged y_fuel along x.
        for cx in 0..self.nx {
            let avg: f64 = (0..self.ny)
                .map(|cy| self.cells[self.idx(cx, cy)].y_fuel)
                .sum::<f64>()
                / self.ny as f64;
            if avg < threshold {
                return Some(cx as f64);
            }
        }
        None
    }

    /// Total heat release rate integrated over all cells (J/m³ equivalent).
    ///
    /// Approximated as Σ (1 − Y_fuel − Y_ox) · Q · ρ.
    pub fn total_heat_release(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| {
                let burned = (1.0 - c.y_fuel - c.y_oxidizer).max(0.0);
                burned * self.reaction.heat_release * c.density
            })
            .sum()
    }
}

// Helper: compute reaction rate for a cell without borrowing `self.reaction`.
impl CombustionCell {
    fn reaction_rate_for(&self, rxn: &ChemicalReaction) -> f64 {
        rxn.reaction_rate(self.y_fuel, self.y_oxidizer, self.temperature)
    }
}

/// Compute the discrete Laplacian with zero-flux (Neumann) boundary conditions.
fn laplacian_2d(field: &[f64], x: usize, y: usize, nx: usize, ny: usize) -> f64 {
    let idx = |cx: usize, cy: usize| cx + cy * nx;

    let xp = if x + 1 < nx { x + 1 } else { x };
    let xm = if x > 0 { x - 1 } else { x };
    let yp = if y + 1 < ny { y + 1 } else { y };
    let ym = if y > 0 { y - 1 } else { y };

    let center = field[idx(x, y)];
    field[idx(xp, y)] + field[idx(xm, y)] + field[idx(x, yp)] + field[idx(x, ym)] - 4.0 * center
}

// ─────────────────────────────────────────────────────────────
// IgnitionModel
// ─────────────────────────────────────────────────────────────

/// Ignition source models.
pub struct IgnitionModel;

impl IgnitionModel {
    /// Apply a spark ignition by raising all cells within `radius` of
    /// `(ix, iy)` to `t_ignition` (K).
    pub fn spark_ignition(
        grid: &mut ReactiveFlowGrid,
        ix: usize,
        iy: usize,
        radius: f64,
        t_ignition: f64,
    ) {
        let nx = grid.nx;
        let ny = grid.ny;
        for cy in 0..ny {
            for cx in 0..nx {
                let dx = cx as f64 - ix as f64;
                let dy = cy as f64 - iy as f64;
                if (dx * dx + dy * dy).sqrt() <= radius {
                    let idx = grid.idx(cx, cy);
                    grid.cells[idx].temperature = t_ignition;
                }
            }
        }
    }

    /// Minimum ignition energy (J) via a simplified quadratic correlation:
    ///
    /// `E_min = 0.3e-3 · (φ − 1)² + 0.1e-3`
    pub fn minimum_ignition_energy(pressure: f64, equivalence_ratio: f64) -> f64 {
        let _ = pressure; // pressure dependence neglected in this correlation
        let dphi = equivalence_ratio - 1.0;
        0.3e-3 * dphi * dphi + 0.1e-3
    }
}

// ─────────────────────────────────────────────────────────────
// Multi-species reactive transport
// ─────────────────────────────────────────────────────────────

/// A multi-species mixture with individual transport properties.
pub struct MultiSpeciesMixture {
    /// Species in the mixture.
    pub species: Vec<Species>,
    /// Mass fractions for each species (same length as `species`).
    pub mass_fractions: Vec<f64>,
}

impl MultiSpeciesMixture {
    /// Create a new mixture from species and initial mass fractions.
    pub fn new(species: Vec<Species>, mass_fractions: Vec<f64>) -> Self {
        assert_eq!(species.len(), mass_fractions.len());
        Self {
            species,
            mass_fractions,
        }
    }

    /// Mean molecular weight of the mixture (g/mol):
    ///
    /// `W_mix = 1 / sum(Y_k / W_k)`
    pub fn mean_molecular_weight(&self) -> f64 {
        let inv_sum: f64 = self
            .species
            .iter()
            .zip(self.mass_fractions.iter())
            .map(|(sp, &y)| y / sp.molecular_weight)
            .sum();
        if inv_sum <= 0.0 {
            return 0.0;
        }
        1.0 / inv_sum
    }

    /// Mixture-averaged specific heat (J/kg/K):
    ///
    /// `cp_mix = sum(Y_k * cp_k)`
    pub fn mixture_cp(&self) -> f64 {
        self.species
            .iter()
            .zip(self.mass_fractions.iter())
            .map(|(sp, &y)| y * sp.specific_heat_cp)
            .sum()
    }

    /// Mole fraction from mass fraction for species index `k`:
    ///
    /// `X_k = (Y_k / W_k) * W_mix`
    pub fn mole_fraction(&self, k: usize) -> f64 {
        let w_mix = self.mean_molecular_weight();
        if w_mix <= 0.0 {
            return 0.0;
        }
        (self.mass_fractions[k] / self.species[k].molecular_weight) * w_mix
    }

    /// Mass-weighted diffusivity for species `k` using its Lewis number:
    ///
    /// `D_k = alpha / Le_k`
    pub fn species_diffusivity(&self, k: usize, thermal_diffusivity: f64) -> f64 {
        let le = self.species[k].lewis_number;
        if le <= 0.0 {
            return thermal_diffusivity;
        }
        thermal_diffusivity / le
    }
}

// ─────────────────────────────────────────────────────────────
// Combustion modeling basics
// ─────────────────────────────────────────────────────────────

/// Multi-step reaction mechanism with multiple elementary reactions.
pub struct ReactionMechanism {
    /// List of elementary reactions.
    pub reactions: Vec<ElementaryReaction>,
}

/// An elementary reaction step with species-specific stoichiometry.
pub struct ElementaryReaction {
    /// Reactant stoichiometric coefficients (indexed by species).
    pub reactant_stoich: Vec<f64>,
    /// Product stoichiometric coefficients (indexed by species).
    pub product_stoich: Vec<f64>,
    /// Activation energy Ea (J/mol).
    pub activation_energy: f64,
    /// Pre-exponential factor A.
    pub pre_exponential: f64,
    /// Temperature exponent n in A * T^n * exp(-Ea/(RuT)).
    pub temp_exponent: f64,
    /// Heat release (J/kg).
    pub heat_release: f64,
}

impl ElementaryReaction {
    /// Create a new elementary reaction.
    pub fn new(
        reactant_stoich: Vec<f64>,
        product_stoich: Vec<f64>,
        activation_energy: f64,
        pre_exponential: f64,
        temp_exponent: f64,
        heat_release: f64,
    ) -> Self {
        Self {
            reactant_stoich,
            product_stoich,
            activation_energy,
            pre_exponential,
            temp_exponent,
            heat_release,
        }
    }

    /// Compute the modified Arrhenius reaction rate:
    ///
    /// `k(T) = A * T^n * exp(-Ea / (Ru * T))`
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        self.pre_exponential
            * temperature.powf(self.temp_exponent)
            * (-self.activation_energy / (R_U * temperature)).exp()
    }

    /// Compute forward reaction rate given species concentrations (mol/m^3):
    ///
    /// `w = k(T) * prod(C_k^nu_k)`
    pub fn forward_rate(&self, concentrations: &[f64], temperature: f64) -> f64 {
        let k = self.rate_constant(temperature);
        let mut prod = 1.0;
        for (i, &nu) in self.reactant_stoich.iter().enumerate() {
            if nu > 0.0 && i < concentrations.len() {
                prod *= concentrations[i].max(0.0).powf(nu);
            }
        }
        k * prod
    }
}

impl ReactionMechanism {
    /// Create a new reaction mechanism.
    pub fn new(reactions: Vec<ElementaryReaction>) -> Self {
        Self { reactions }
    }

    /// Compute net species production rates (mol/m^3/s).
    pub fn species_production_rates(
        &self,
        concentrations: &[f64],
        temperature: f64,
        n_species: usize,
    ) -> Vec<f64> {
        let mut rates = vec![0.0; n_species];
        for rxn in &self.reactions {
            let w = rxn.forward_rate(concentrations, temperature);
            for (k, rate_k) in rates.iter_mut().enumerate() {
                let nu_prod = if k < rxn.product_stoich.len() {
                    rxn.product_stoich[k]
                } else {
                    0.0
                };
                let nu_react = if k < rxn.reactant_stoich.len() {
                    rxn.reactant_stoich[k]
                } else {
                    0.0
                };
                *rate_k += (nu_prod - nu_react) * w;
            }
        }
        rates
    }

    /// Total heat release rate (J/m^3/s).
    pub fn total_heat_release_rate(&self, concentrations: &[f64], temperature: f64) -> f64 {
        self.reactions
            .iter()
            .map(|rxn| rxn.forward_rate(concentrations, temperature) * rxn.heat_release)
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────
// Flame speed correlations
// ─────────────────────────────────────────────────────────────

impl FlameProperties {
    /// Turbulent flame speed (m/s) using the Damkoehler correlation:
    ///
    /// `S_T = S_L + u'`
    ///
    /// where `u_prime` is the turbulence intensity (rms velocity fluctuation).
    pub fn turbulent_flame_speed(s_l: f64, u_prime: f64) -> f64 {
        s_l + u_prime
    }

    /// Flame thickness estimate (m):
    ///
    /// `delta_f = alpha / S_L`
    pub fn flame_thickness(thermal_diffusivity: f64, s_l: f64) -> f64 {
        if s_l <= 0.0 {
            return f64::MAX;
        }
        thermal_diffusivity / s_l
    }

    /// Damkoehler number:
    ///
    /// `Da = (l_t / u') * (S_L / delta_f)`
    pub fn damkoehler_number(
        integral_length: f64,
        u_prime: f64,
        s_l: f64,
        flame_thickness: f64,
    ) -> f64 {
        if u_prime <= 0.0 || flame_thickness <= 0.0 {
            return f64::MAX;
        }
        (integral_length / u_prime) * (s_l / flame_thickness)
    }

    /// Karlovitz number:
    ///
    /// `Ka = (delta_f / eta)^2`
    ///
    /// where eta is the Kolmogorov length scale.
    pub fn karlovitz_number(flame_thickness: f64, kolmogorov_scale: f64) -> f64 {
        if kolmogorov_scale <= 0.0 {
            return f64::MAX;
        }
        let ratio = flame_thickness / kolmogorov_scale;
        ratio * ratio
    }
}

// ─────────────────────────────────────────────────────────────
// Reaction rate limiting
// ─────────────────────────────────────────────────────────────

/// Rate limiters for numerical stability.
pub struct ReactionRateLimiter;

impl ReactionRateLimiter {
    /// Limit the reaction rate so that no species is consumed beyond
    /// its available mass fraction in a single timestep.
    ///
    /// Returns the limited rate.
    pub fn limit_by_available_reactant(rate: f64, y_available: f64, dt: f64) -> f64 {
        if dt <= 0.0 {
            return rate;
        }
        let max_rate = y_available / dt;
        rate.min(max_rate)
    }

    /// Temperature-based rate cutoff: returns 0 if T < T_cutoff.
    pub fn temperature_cutoff(rate: f64, temperature: f64, t_cutoff: f64) -> f64 {
        if temperature < t_cutoff { 0.0 } else { rate }
    }

    /// Clip a reaction rate to `[0, max_rate]`.
    pub fn clip_rate(rate: f64, max_rate: f64) -> f64 {
        rate.clamp(0.0, max_rate)
    }

    /// Apply all standard limiters at once.
    pub fn apply_all(
        rate: f64,
        y_fuel: f64,
        y_ox: f64,
        temperature: f64,
        dt: f64,
        t_cutoff: f64,
        max_rate: f64,
    ) -> f64 {
        let r = Self::temperature_cutoff(rate, temperature, t_cutoff);
        let r = Self::limit_by_available_reactant(r, y_fuel, dt);
        let r = Self::limit_by_available_reactant(r, y_ox, dt);
        Self::clip_rate(r, max_rate)
    }
}

// ─────────────────────────────────────────────────────────────
// Species diffusion (multi-component)
// ─────────────────────────────────────────────────────────────

/// Multi-component species diffusion on a 2D grid.
pub struct SpeciesDiffusion;

impl SpeciesDiffusion {
    /// Diffuse a single scalar field using explicit Laplacian with
    /// species-specific diffusivity.
    pub fn diffuse_field(
        field: &mut [f64],
        field_old: &[f64],
        nx: usize,
        ny: usize,
        diffusivity: f64,
        dt: f64,
        dx: f64,
    ) {
        let coeff = diffusivity * dt / (dx * dx);
        let n = nx * ny;
        for idx in 0..n {
            let x = idx % nx;
            let y = idx / nx;
            let lap = laplacian_2d(field_old, x, y, nx, ny);
            field[idx] = (field_old[idx] + coeff * lap).max(0.0);
        }
    }

    /// Compute the mass-weighted diffusion velocity for species k:
    ///
    /// `V_k = -D_k / Y_k * grad(Y_k)`
    ///
    /// Returns `[Vx, Vy]` at a single cell using central differences.
    pub fn diffusion_velocity(
        y_field: &[f64],
        x: usize,
        y: usize,
        nx: usize,
        ny: usize,
        dx: f64,
        diffusivity: f64,
    ) -> [f64; 2] {
        let idx = |cx: usize, cy: usize| cx + cy * nx;
        let y_c = y_field[idx(x, y)];
        if y_c <= 1e-30 {
            return [0.0, 0.0];
        }

        let xp = if x + 1 < nx { x + 1 } else { x };
        let xm = if x > 0 { x - 1 } else { x };
        let yp = if y + 1 < ny { y + 1 } else { y };
        let ym = if y > 0 { y - 1 } else { y };

        let dx_inv = if xp != xm {
            1.0 / ((xp - xm) as f64 * dx)
        } else {
            0.0
        };
        let dy_inv = if yp != ym {
            1.0 / ((yp - ym) as f64 * dx)
        } else {
            0.0
        };

        let grad_x = (y_field[idx(xp, y)] - y_field[idx(xm, y)]) * dx_inv;
        let grad_y = (y_field[idx(x, yp)] - y_field[idx(x, ym)]) * dy_inv;

        [-diffusivity / y_c * grad_x, -diffusivity / y_c * grad_y]
    }
}

// ─────────────────────────────────────────────────────────────
// Extended ReactiveFlowGrid methods
// ─────────────────────────────────────────────────────────────

impl ReactiveFlowGrid {
    /// Mean temperature across all cells.
    pub fn mean_temperature(&self) -> f64 {
        let n = self.cells.len();
        if n == 0 {
            return 0.0;
        }
        self.cells.iter().map(|c| c.temperature).sum::<f64>() / n as f64
    }

    /// Maximum temperature across all cells.
    pub fn max_temperature(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| c.temperature)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Total fuel mass remaining: sum of y_fuel * density over all cells.
    pub fn total_fuel_mass(&self) -> f64 {
        self.cells.iter().map(|c| c.y_fuel * c.density).sum()
    }

    /// Initialize a diffusion flame: fuel on left, oxidizer on right,
    /// with a mixing layer of width `mix_width` centered at `interface_x`.
    pub fn initialize_diffusion_flame(
        &mut self,
        y_fuel: f64,
        y_ox: f64,
        t_initial: f64,
        interface_x: usize,
        mix_width: usize,
    ) {
        for cy in 0..self.ny {
            for cx in 0..self.nx {
                let idx = self.idx(cx, cy);
                self.cells[idx].temperature = t_initial;
                if cx < interface_x.saturating_sub(mix_width / 2) {
                    self.cells[idx].y_fuel = y_fuel;
                    self.cells[idx].y_oxidizer = 0.0;
                } else if cx > interface_x + mix_width / 2 {
                    self.cells[idx].y_fuel = 0.0;
                    self.cells[idx].y_oxidizer = y_ox;
                } else {
                    // Linear blend in mixing layer
                    let span = mix_width.max(1) as f64;
                    let start = interface_x.saturating_sub(mix_width / 2) as f64;
                    let frac = (cx as f64 - start) / span;
                    self.cells[idx].y_fuel = y_fuel * (1.0 - frac);
                    self.cells[idx].y_oxidizer = y_ox * frac;
                }
                self.cells[idx].y_products = 0.0;
            }
        }
    }

    /// Compute the global heat release rate (sum of w * Q over all cells).
    pub fn global_heat_release_rate(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| {
                let w = c.reaction_rate_for(&self.reaction);
                w * self.reaction.heat_release
            })
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn methane_reaction() -> ChemicalReaction {
        ChemicalReaction::new(
            1.0,    // fuel_stoich
            2.0,    // oxidizer_stoich
            1.0,    // product_stoich
            1.26e5, // Ea (J/mol)
            1e10,   // A
            5e7,    // Q (J/kg_fuel)
        )
    }

    // 1. Reaction rate is positive for nonzero concentrations at high T.
    #[test]
    fn test_reaction_rate_positive_high_temperature() {
        let rxn = methane_reaction();
        let w = rxn.reaction_rate(0.05, 0.2, 2000.0);
        assert!(w > 0.0, "Expected positive reaction rate, got {w}");
    }

    // 2. Reaction rate approaches 0 at very low temperature.
    #[test]
    fn test_reaction_rate_low_temperature() {
        let rxn = methane_reaction();
        let w_low = rxn.reaction_rate(0.05, 0.2, 100.0);
        let w_high = rxn.reaction_rate(0.05, 0.2, 2000.0);
        assert!(
            w_low < w_high * 1e-30,
            "Rate at low T ({w_low}) should be negligible compared to high T ({w_high})"
        );
    }

    // 3. Adiabatic flame temperature > T_initial for positive heat release.
    #[test]
    fn test_adiabatic_flame_temperature() {
        let t0 = 300.0;
        let t_ad = FlameProperties::adiabatic_flame_temperature(0.05, 1200.0, 5e7, t0);
        assert!(
            t_ad > t0,
            "Adiabatic flame temp should exceed initial: {t_ad} <= {t0}"
        );
    }

    // 4. Equivalence ratio = 1 at stoichiometric conditions.
    #[test]
    fn test_equivalence_ratio_stoichiometric() {
        let stoich = 0.058;
        let phi = FlameProperties::equivalence_ratio(stoich * 0.2, 0.2, stoich);
        assert!((phi - 1.0).abs() < 1e-12, "Expected phi = 1, got {phi}");
    }

    // 5. Laminar flame speed is positive.
    #[test]
    fn test_laminar_flame_speed_positive() {
        let sl = FlameProperties::laminar_flame_speed(0.05, 1000.0, 101_325.0);
        assert!(sl > 0.0, "Laminar flame speed should be positive: {sl}");
    }

    // 6. initialize_premixed: left side has fuel, right side is hot.
    #[test]
    fn test_initialize_premixed() {
        let mut grid = ReactiveFlowGrid::new(20, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 300.0, 2200.0, 10);

        let left_idx = grid.idx(2, 0);
        assert!(
            grid.cells[left_idx].y_fuel > 0.0,
            "Left side should have fuel"
        );
        assert!(
            (grid.cells[left_idx].temperature - 300.0).abs() < 1e-10,
            "Left side should be at fresh temperature"
        );

        let right_idx = grid.idx(15, 0);
        assert!(
            grid.cells[right_idx].y_fuel == 0.0,
            "Right side should have no fuel"
        );
        assert!(
            grid.cells[right_idx].temperature > 2000.0,
            "Right side should be hot"
        );
    }

    // 7. step_chemistry reduces y_fuel when conditions are reactive.
    #[test]
    fn test_step_chemistry_reduces_fuel() {
        let mut grid = ReactiveFlowGrid::new(4, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 2000.0, 2200.0, 2);

        let fuel_before: f64 = grid.cells.iter().map(|c| c.y_fuel).sum();
        grid.step_chemistry(1e-6);
        let fuel_after: f64 = grid.cells.iter().map(|c| c.y_fuel).sum();

        assert!(
            fuel_after < fuel_before,
            "Fuel should decrease after chemistry step: before={fuel_before}, after={fuel_after}"
        );
    }

    // 8. total_heat_release >= 0.
    #[test]
    fn test_total_heat_release_nonneg() {
        let mut grid = ReactiveFlowGrid::new(10, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 300.0, 2200.0, 5);
        let q = grid.total_heat_release();
        assert!(q >= 0.0, "Total heat release should be non-negative: {q}");
    }

    // ─── Multi-species mixture tests ───

    #[test]
    fn test_multi_species_mean_molecular_weight() {
        let ch4 = Species::new("CH4", 16.0, 2200.0, -74.8e3, 1.0);
        let o2 = Species::new("O2", 32.0, 920.0, 0.0, 1.0);
        let mix = MultiSpeciesMixture::new(vec![ch4, o2], vec![0.5, 0.5]);
        let w = mix.mean_molecular_weight();
        // W = 1 / (0.5/16 + 0.5/32) = 1 / (0.03125 + 0.015625) = 1/0.046875 = 21.333...
        assert!((w - 1.0 / 0.046875).abs() < 1e-8, "Mean MW = {w}");
    }

    #[test]
    fn test_multi_species_mixture_cp() {
        let ch4 = Species::new("CH4", 16.0, 2200.0, -74.8e3, 1.0);
        let o2 = Species::new("O2", 32.0, 920.0, 0.0, 1.0);
        let mix = MultiSpeciesMixture::new(vec![ch4, o2], vec![0.3, 0.7]);
        let cp = mix.mixture_cp();
        let expected = 0.3 * 2200.0 + 0.7 * 920.0;
        assert!(
            (cp - expected).abs() < 1e-10,
            "cp_mix = {cp}, expected {expected}"
        );
    }

    #[test]
    fn test_multi_species_mole_fraction_sum() {
        let ch4 = Species::new("CH4", 16.0, 2200.0, 0.0, 1.0);
        let o2 = Species::new("O2", 32.0, 920.0, 0.0, 1.0);
        let n2 = Species::new("N2", 28.0, 1040.0, 0.0, 1.0);
        let yf = vec![0.2, 0.3, 0.5];
        let mix = MultiSpeciesMixture::new(vec![ch4, o2, n2], yf);
        let x_sum: f64 = (0..3).map(|k| mix.mole_fraction(k)).sum();
        assert!(
            (x_sum - 1.0).abs() < 1e-10,
            "Mole fractions should sum to 1: {x_sum}"
        );
    }

    #[test]
    fn test_species_diffusivity_unity_lewis() {
        let sp = Species::new("X", 28.0, 1000.0, 0.0, 1.0);
        let mix = MultiSpeciesMixture::new(vec![sp], vec![1.0]);
        let alpha = 2.0e-5;
        let d = mix.species_diffusivity(0, alpha);
        assert!((d - alpha).abs() < 1e-20, "Le=1 should give D=alpha");
    }

    // ─── Reaction mechanism tests ───

    #[test]
    fn test_elementary_reaction_rate_constant() {
        let rxn = ElementaryReaction::new(vec![1.0, 0.5], vec![0.0, 0.0], 1.0e5, 1e12, 0.0, 5e6);
        let k1 = rxn.rate_constant(1000.0);
        let k2 = rxn.rate_constant(2000.0);
        assert!(k2 > k1, "Rate constant should increase with T");
        assert!(k1 > 0.0);
    }

    #[test]
    fn test_reaction_mechanism_production_rates() {
        let rxn =
            ElementaryReaction::new(vec![1.0, 1.0], vec![0.0, 0.0, 1.0], 5.0e4, 1e8, 0.0, 1e6);
        let mech = ReactionMechanism::new(vec![rxn]);
        let conc = vec![1.0, 2.0, 0.0];
        let rates = mech.species_production_rates(&conc, 2000.0, 3);
        // Reactants consumed (negative rates), products formed (positive)
        assert!(rates[0] < 0.0, "Reactant 0 should be consumed");
        assert!(rates[1] < 0.0, "Reactant 1 should be consumed");
        assert!(rates[2] > 0.0, "Product should be formed");
    }

    #[test]
    fn test_total_heat_release_rate_positive() {
        let rxn = ElementaryReaction::new(vec![1.0], vec![0.0, 1.0], 5.0e4, 1e8, 0.0, 5e6);
        let mech = ReactionMechanism::new(vec![rxn]);
        let conc = vec![1.0, 0.0];
        let q = mech.total_heat_release_rate(&conc, 2000.0);
        assert!(q > 0.0, "Heat release rate should be positive: {q}");
    }

    // ─── Flame property extensions ───

    #[test]
    fn test_turbulent_flame_speed() {
        let s_l = 0.4;
        let u_prime = 1.0;
        let s_t = FlameProperties::turbulent_flame_speed(s_l, u_prime);
        assert!((s_t - 1.4).abs() < 1e-12, "S_T = {s_t}");
    }

    #[test]
    fn test_flame_thickness_positive() {
        let delta = FlameProperties::flame_thickness(2e-5, 0.4);
        assert!(delta > 0.0, "Flame thickness should be positive: {delta}");
        assert!((delta - 5e-5).abs() < 1e-10);
    }

    #[test]
    fn test_damkoehler_number() {
        let da = FlameProperties::damkoehler_number(0.01, 1.0, 0.4, 5e-5);
        assert!(da > 0.0, "Da should be positive: {da}");
    }

    #[test]
    fn test_karlovitz_number() {
        let ka = FlameProperties::karlovitz_number(5e-5, 1e-5);
        let expected = (5e-5 / 1e-5) * (5e-5 / 1e-5);
        assert!(
            (ka - expected).abs() < 1e-6,
            "Ka = {ka}, expected {expected}"
        );
    }

    // ─── Rate limiter tests ───

    #[test]
    fn test_rate_limiter_available_reactant() {
        let limited = ReactionRateLimiter::limit_by_available_reactant(100.0, 0.05, 1e-3);
        assert!(
            limited <= 0.05 / 1e-3 + 1e-10,
            "Rate should be limited: {limited}"
        );
    }

    #[test]
    fn test_rate_limiter_temperature_cutoff() {
        let r = ReactionRateLimiter::temperature_cutoff(100.0, 200.0, 500.0);
        assert_eq!(r, 0.0, "Rate below cutoff should be 0");
        let r2 = ReactionRateLimiter::temperature_cutoff(100.0, 600.0, 500.0);
        assert_eq!(r2, 100.0, "Rate above cutoff should pass through");
    }

    #[test]
    fn test_rate_limiter_clip() {
        let r = ReactionRateLimiter::clip_rate(1000.0, 500.0);
        assert_eq!(r, 500.0);
        let r2 = ReactionRateLimiter::clip_rate(-5.0, 500.0);
        assert_eq!(r2, 0.0);
    }

    #[test]
    fn test_rate_limiter_apply_all() {
        let r = ReactionRateLimiter::apply_all(1e10, 0.05, 0.2, 2000.0, 1e-3, 500.0, 1e8);
        assert!(
            (0.0..=1e8).contains(&r),
            "Rate should be in valid range: {r}"
        );
    }

    // ─── Species diffusion tests ───

    #[test]
    fn test_diffuse_field_uniform_no_change() {
        let nx = 5;
        let ny = 5;
        let field_old = vec![1.0; nx * ny];
        let mut field = field_old.clone();
        SpeciesDiffusion::diffuse_field(&mut field, &field_old, nx, ny, 1e-5, 1e-4, 1.0);
        for (i, &v) in field.iter().enumerate() {
            assert!(
                (v - 1.0).abs() < 1e-10,
                "Uniform field should not change: cell {i} = {v}"
            );
        }
    }

    #[test]
    fn test_diffuse_field_smooths_spike() {
        let nx = 10;
        let ny = 10;
        let mut field_old = vec![0.0; nx * ny];
        field_old[5 * nx + 5] = 1.0; // spike at center
        let mut field = field_old.clone();
        SpeciesDiffusion::diffuse_field(&mut field, &field_old, nx, ny, 0.1, 0.1, 1.0);
        // Spike should decrease
        assert!(field[5 * nx + 5] < 1.0, "Spike should be smoothed");
    }

    #[test]
    fn test_diffusion_velocity_uniform_zero() {
        let nx = 5;
        let ny = 5;
        let field = vec![1.0; nx * ny];
        let v = SpeciesDiffusion::diffusion_velocity(&field, 2, 2, nx, ny, 1.0, 1e-5);
        assert!(
            v[0].abs() < 1e-15 && v[1].abs() < 1e-15,
            "Uniform field: V = [{}, {}]",
            v[0],
            v[1]
        );
    }

    // ─── Extended grid tests ───

    #[test]
    fn test_mean_temperature() {
        let grid = ReactiveFlowGrid::new(4, 4, methane_reaction());
        let t_mean = grid.mean_temperature();
        assert!((t_mean - 300.0).abs() < 1e-10, "Mean T = {t_mean}");
    }

    #[test]
    fn test_max_temperature() {
        let mut grid = ReactiveFlowGrid::new(10, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 300.0, 2200.0, 5);
        let t_max = grid.max_temperature();
        assert!(t_max >= 2200.0, "Max T should be >= 2200: {t_max}");
    }

    #[test]
    fn test_initialize_diffusion_flame() {
        let mut grid = ReactiveFlowGrid::new(20, 4, methane_reaction());
        grid.initialize_diffusion_flame(0.1, 0.23, 300.0, 10, 4);

        // Left: fuel rich
        let left = grid.idx(2, 0);
        assert!(grid.cells[left].y_fuel > 0.05, "Left should have fuel");
        assert!(
            grid.cells[left].y_oxidizer < 1e-10,
            "Left should have no oxidizer"
        );

        // Right: oxidizer
        let right = grid.idx(18, 0);
        assert!(
            grid.cells[right].y_oxidizer > 0.1,
            "Right should have oxidizer"
        );
        assert!(
            grid.cells[right].y_fuel < 1e-10,
            "Right should have no fuel"
        );
    }

    #[test]
    fn test_global_heat_release_rate_nonneg() {
        let mut grid = ReactiveFlowGrid::new(10, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 2000.0, 2200.0, 5);
        let q = grid.global_heat_release_rate();
        assert!(
            q >= 0.0,
            "Global heat release rate should be non-negative: {q}"
        );
    }

    #[test]
    fn test_total_fuel_mass_decreases_after_chemistry() {
        let mut grid = ReactiveFlowGrid::new(6, 4, methane_reaction());
        grid.initialize_premixed(0.05, 0.2, 2000.0, 2200.0, 3);
        let fuel_before = grid.total_fuel_mass();
        grid.step_chemistry(1e-6);
        let fuel_after = grid.total_fuel_mass();
        assert!(
            fuel_after < fuel_before,
            "Fuel mass should decrease: {fuel_before} -> {fuel_after}"
        );
    }
}

// ---------------------------------------------------------------------------
// Reduced-order combustion chemistry models
// ---------------------------------------------------------------------------

/// Damköhler number classifier for reactive flow regimes.
///
/// Da = τ_flow / τ_chem.
/// - Da >> 1 → fast chemistry (well-stirred reactor / equilibrium)
/// - Da << 1 → slow chemistry (frozen flow)
pub fn damkohler_number(tau_flow: f64, tau_chem: f64) -> f64 {
    tau_flow / tau_chem
}

/// Zeldovich number Ze = Eₐ(T_ad - T_0) / (Rᵤ T_ad²).
///
/// Measures the activation energy non-dimensionally; large Ze means thin flames.
pub fn zeldovich_number(activation_energy: f64, t_ad: f64, t_0: f64) -> f64 {
    activation_energy * (t_ad - t_0) / (R_U * t_ad * t_ad)
}

/// Laminar flame speed estimate (Mallard–Le Chatelier correlation):
///
/// S_L ≈ √(2 α_th ω̄ / (ρ Y_fuel_0))
///
/// where α_th = λ/(ρ cp) is thermal diffusivity and ω̄ is the mean reaction rate.
pub fn laminar_flame_speed(
    thermal_diffusivity: f64,
    mean_reaction_rate: f64,
    rho: f64,
    y_fuel_0: f64,
) -> f64 {
    if rho <= 0.0 || y_fuel_0 <= 0.0 {
        return 0.0;
    }
    (2.0 * thermal_diffusivity * mean_reaction_rate / (rho * y_fuel_0)).sqrt()
}

// ---------------------------------------------------------------------------
// Chapman–Jouguet detonation properties
// ---------------------------------------------------------------------------

/// Compute the Chapman–Jouguet (CJ) detonation velocity.
///
/// D_CJ = c_0 + √(Q(γ²-1)/(2γ))   (simplified 1-step chemistry formula)
///
/// where c_0 = √(γ·p₀/ρ₀) is the ambient sound speed,
/// Q is the heat release per unit mass, and γ is the specific heat ratio.
pub fn cj_detonation_velocity(gamma: f64, pressure: f64, density: f64, heat_release: f64) -> f64 {
    let c0 = (gamma * pressure / density).sqrt();
    c0 + (heat_release * (gamma * gamma - 1.0) / (2.0 * gamma)).sqrt()
}

/// CJ Mach number: M_CJ = D_CJ / c_0.
pub fn cj_mach_number(gamma: f64, pressure: f64, density: f64, heat_release: f64) -> f64 {
    let c0 = (gamma * pressure / density).sqrt();
    cj_detonation_velocity(gamma, pressure, density, heat_release) / c0
}

// ---------------------------------------------------------------------------
// Mixture fraction and scalar dissipation
// ---------------------------------------------------------------------------

/// Compute the Bilger mixture fraction Z from fuel and oxidizer mass fractions.
///
/// Z = (Y_fuel - Y_ox / s + Y_ox_inf / s) / (Y_fuel_inf + Y_ox_inf / s)
///
/// where s = stoichiometric oxidizer-to-fuel ratio.
pub fn bilger_mixture_fraction(
    y_fuel: f64,
    y_ox: f64,
    y_fuel_inf: f64,
    y_ox_inf: f64,
    stoich_ratio: f64,
) -> f64 {
    let num = y_fuel - y_ox / stoich_ratio + y_ox_inf / stoich_ratio;
    let den = y_fuel_inf + y_ox_inf / stoich_ratio;
    if den.abs() < 1e-15 {
        return 0.0;
    }
    (num / den).clamp(0.0, 1.0)
}

/// Compute the scalar dissipation rate χ from mixture fraction gradient.
///
/// χ = 2 D |∇Z|²
///
/// where D is the mass diffusivity.
pub fn scalar_dissipation_rate(diffusivity: f64, grad_z: [f64; 2]) -> f64 {
    2.0 * diffusivity * (grad_z[0] * grad_z[0] + grad_z[1] * grad_z[1])
}

// ---------------------------------------------------------------------------
// Flame front tracking
// ---------------------------------------------------------------------------

/// Detect whether a cell contains an active flame front.
///
/// A flame is active when reaction rate exceeds a threshold fraction
/// of the theoretical maximum rate.
pub fn is_flame_front(reaction_rate: f64, max_rate: f64, threshold_fraction: f64) -> bool {
    if max_rate <= 0.0 {
        return false;
    }
    reaction_rate / max_rate >= threshold_fraction
}

/// Estimate the local flame thickness δ_L ≈ α_th / S_L.
pub fn flame_thickness(thermal_diffusivity: f64, flame_speed: f64) -> f64 {
    if flame_speed <= 0.0 {
        return f64::INFINITY;
    }
    thermal_diffusivity / flame_speed
}

// ---------------------------------------------------------------------------
// Extended reactive_flow tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod reactive_extra_tests {
    use super::*;

    fn methane_reaction_local() -> ChemicalReaction {
        ChemicalReaction::new(1.0, 2.0, 1.0, 125_000.0, 1e12, 50_000_000.0)
    }

    /// Damköhler number is positive and dimensionless.
    #[test]
    fn test_damkohler_number() {
        let da = damkohler_number(1e-3, 1e-5);
        assert!((da - 100.0).abs() < 1e-10, "Da={da}");
    }

    /// Zeldovich number is finite and positive for valid inputs.
    #[test]
    fn test_zeldovich_number() {
        let ze = zeldovich_number(125_000.0, 2200.0, 300.0);
        assert!(ze.is_finite() && ze > 0.0, "Ze={ze}");
    }

    /// Laminar flame speed is positive for reactive mixture.
    #[test]
    fn test_laminar_flame_speed_positive() {
        let s_l = laminar_flame_speed(1e-5, 1e6, 1.0, 0.1);
        assert!(s_l > 0.0, "S_L={s_l}");
    }

    /// CJ detonation velocity exceeds ambient sound speed.
    #[test]
    fn test_cj_velocity_exceeds_c0() {
        let gamma = 1.4_f64;
        let p0 = 101325.0_f64;
        let rho0 = 1.2_f64;
        let q = 2e6_f64;
        let c0 = (gamma * p0 / rho0).sqrt();
        let d_cj = cj_detonation_velocity(gamma, p0, rho0, q);
        assert!(d_cj > c0, "D_CJ={d_cj} should exceed c0={c0}");
    }

    /// CJ Mach number > 1 (detonation is supersonic).
    #[test]
    fn test_cj_mach_number_supersonic() {
        let m_cj = cj_mach_number(1.4, 101325.0, 1.2, 2e6);
        assert!(m_cj > 1.0, "M_CJ={m_cj} should be > 1");
    }

    /// Bilger mixture fraction: pure fuel side → Z = 1.
    #[test]
    fn test_bilger_pure_fuel() {
        // y_fuel = y_fuel_inf, y_ox = 0
        let z = bilger_mixture_fraction(0.1, 0.0, 0.1, 0.23, 17.16);
        assert!((z - 1.0).abs() < 1e-6, "Z_fuel={z}");
    }

    /// Bilger mixture fraction: pure oxidizer side → Z ≈ 0.
    #[test]
    fn test_bilger_pure_oxidizer() {
        let z = bilger_mixture_fraction(0.0, 0.23, 0.1, 0.23, 17.16);
        assert!(z.abs() < 1e-8, "Z_ox={z}");
    }

    /// Scalar dissipation rate: zero gradient → zero χ.
    #[test]
    fn test_scalar_dissipation_zero_gradient() {
        let chi = scalar_dissipation_rate(1e-5, [0.0, 0.0]);
        assert!(chi.abs() < 1e-20, "chi={chi}");
    }

    /// Scalar dissipation rate is positive for non-zero gradient.
    #[test]
    fn test_scalar_dissipation_positive() {
        let chi = scalar_dissipation_rate(1e-5, [0.1, 0.05]);
        assert!(chi > 0.0, "chi={chi}");
    }

    /// is_flame_front: high reaction rate is detected as flame.
    #[test]
    fn test_flame_front_detection() {
        assert!(is_flame_front(9e5, 1e6, 0.5), "Should be flame front");
        assert!(!is_flame_front(1e4, 1e6, 0.5), "Should NOT be flame front");
    }

    /// flame_thickness decreases with increasing flame speed.
    #[test]
    fn test_flame_thickness_decreases_with_speed() {
        let d1 = flame_thickness(1e-5, 0.1);
        let d2 = flame_thickness(1e-5, 0.5);
        assert!(d1 > d2, "δ1={d1}, δ2={d2}");
    }

    /// Arrhenius rate: reaction rate zero at T=0.
    #[test]
    fn test_arrhenius_zero_temperature() {
        let rxn = methane_reaction_local();
        let w = rxn.reaction_rate(0.1, 0.23, 0.0);
        assert!(w.abs() < 1e-15, "Rate at T=0: {w}");
    }
}
