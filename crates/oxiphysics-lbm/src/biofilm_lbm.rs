#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biofilm growth and nutrient transport simulation using the Lattice Boltzmann Method.
//!
//! This module provides:
//! - [`BiofilmCell`]: Per-grid-cell state (biomass, nutrient, EPS, fluid/solid).
//! - [`NutrientTransport`]: Diffusion and Monod-kinetic consumption model.
//! - [`BiofilmGrowth`]: Monod growth and endogenous decay.
//! - [`BiofilmSimulation`]: Full 2-D LBM-based biofilm simulator.
//! - [`BiofilmProperties`]: Derived macroscopic properties (porosity, permeability, tortuosity).
//! - [`monod_kinetics`]: Standalone Monod rate function.
//! - [`kozeny_carman_permeability`]: Kozeny–Carman equation for porous-media permeability.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Standalone free functions
// ---------------------------------------------------------------------------

/// Compute the Monod-kinetics growth rate (s⁻¹).
///
/// ```no_run
/// use oxiphysics_lbm::biofilm_lbm::monod_kinetics;
/// let rate = monod_kinetics(0.5, 0.5, 1.0); // conc == Ks → rate = 0.5 * mu_max
/// assert!((rate - 0.5).abs() < 1e-12);
/// ```
pub fn monod_kinetics(conc: f64, ks: f64, mu_max: f64) -> f64 {
    if conc <= 0.0 || ks <= 0.0 {
        return 0.0;
    }
    mu_max * conc / (ks + conc)
}

/// Kozeny–Carman permeability (m²).
///
/// `K = (ε³ d²) / (180 (1-ε)²)`
///
/// where ε = `porosity` (dimensionless, 0..1) and d = `d_particle` (m).
///
/// Returns 0 for unphysical inputs (porosity outside (0, 1)).
///
/// ```no_run
/// use oxiphysics_lbm::biofilm_lbm::kozeny_carman_permeability;
/// let k = kozeny_carman_permeability(0.4, 1e-4);
/// assert!(k > 0.0);
/// ```
pub fn kozeny_carman_permeability(porosity: f64, d_particle: f64) -> f64 {
    if porosity <= 0.0 || porosity >= 1.0 || d_particle <= 0.0 {
        return 0.0;
    }
    let eps = porosity;
    let one_m_eps = 1.0 - eps;
    (eps * eps * eps * d_particle * d_particle) / (180.0 * one_m_eps * one_m_eps)
}

// ---------------------------------------------------------------------------
// BiofilmCell
// ---------------------------------------------------------------------------

/// State of a single grid cell in the biofilm simulation.
///
/// Each cell stores the local biomass density (g/m³), dissolved-nutrient
/// concentration (g/m³), EPS fraction (dimensionless 0–1), and whether the
/// cell is occupied by fluid (true) or solid biofilm matrix (false).
#[derive(Debug, Clone, PartialEq)]
pub struct BiofilmCell {
    /// Biomass density (g m⁻³).
    pub biomass_density: f64,
    /// Dissolved-nutrient (substrate) concentration (g m⁻³).
    pub nutrient_concentration: f64,
    /// Volume fraction of extracellular polymeric substances (EPS), 0–1.
    pub eps_fraction: f64,
    /// `true` = fluid cell; `false` = biofilm/solid cell.
    pub is_fluid: bool,
}

impl BiofilmCell {
    /// Create a new fluid cell with specified initial conditions.
    pub fn fluid(biomass: f64, nutrient: f64, eps: f64) -> Self {
        Self {
            biomass_density: biomass,
            nutrient_concentration: nutrient,
            eps_fraction: eps.clamp(0.0, 1.0),
            is_fluid: true,
        }
    }

    /// Create a new solid (biofilm) cell.
    pub fn solid(biomass: f64, nutrient: f64, eps: f64) -> Self {
        Self {
            biomass_density: biomass,
            nutrient_concentration: nutrient,
            eps_fraction: eps.clamp(0.0, 1.0),
            is_fluid: false,
        }
    }

    /// Create an empty fluid cell (zero biomass, zero nutrient, zero EPS).
    pub fn empty_fluid() -> Self {
        Self::fluid(0.0, 0.0, 0.0)
    }
}

impl Default for BiofilmCell {
    fn default() -> Self {
        Self::empty_fluid()
    }
}

// ---------------------------------------------------------------------------
// NutrientTransport
// ---------------------------------------------------------------------------

/// Nutrient diffusion-reaction parameters.
///
/// Models nutrient (substrate) transport as reaction–diffusion:
///   ∂C/∂t = D ∇²C − r_consumption(C, X)
/// where `r_consumption` follows Monod kinetics.
#[derive(Debug, Clone)]
pub struct NutrientTransport {
    /// Effective diffusivity of the nutrient (m² s⁻¹).
    pub diffusivity: f64,
    /// Maximum volumetric consumption rate (g m⁻³ s⁻¹).
    pub consumption_rate: f64,
    /// Monod half-saturation constant Ks (g m⁻³).
    pub monod_constant: f64,
}

impl NutrientTransport {
    /// Create a new [`NutrientTransport`] model.
    pub fn new(diffusivity: f64, consumption_rate: f64, monod_constant: f64) -> Self {
        Self {
            diffusivity,
            consumption_rate,
            monod_constant,
        }
    }

    /// Compute the reaction sink term for nutrient:
    ///   r(C, X) = consumption_rate * C / (Ks + C) * (X / X_ref)
    ///
    /// where X_ref = 1 g m⁻³ (normalising reference).
    pub fn compute_reaction(&self, conc: f64, biomass: f64) -> f64 {
        if conc <= 0.0 || biomass <= 0.0 {
            return 0.0;
        }
        let monod = monod_kinetics(conc, self.monod_constant, 1.0);
        self.consumption_rate * monod * biomass
    }
}

// ---------------------------------------------------------------------------
// BiofilmGrowth
// ---------------------------------------------------------------------------

/// Biofilm growth kinetics (Monod model with endogenous decay).
///
/// Biomass evolution:
///   dX/dt = (μ_max · C / (Ks + C)) · X − b_d · X
///         = (net_growth_rate − decay_rate) · X
#[derive(Debug, Clone)]
pub struct BiofilmGrowth {
    /// Maximum specific growth rate μ_max (s⁻¹).
    pub max_growth_rate: f64,
    /// Yield coefficient Y (g biomass / g substrate).
    pub yield_coeff: f64,
    /// Endogenous decay rate b_d (s⁻¹).
    pub decay_rate: f64,
    /// Half-saturation constant for growth Ks (g m⁻³).
    pub monod_constant: f64,
}

impl BiofilmGrowth {
    /// Create a new [`BiofilmGrowth`] model.
    pub fn new(
        max_growth_rate: f64,
        yield_coeff: f64,
        decay_rate: f64,
        monod_constant: f64,
    ) -> Self {
        Self {
            max_growth_rate,
            yield_coeff,
            decay_rate,
            monod_constant,
        }
    }

    /// Compute net biomass change rate (g m⁻³ s⁻¹):
    ///   dX/dt = Y · q_s − b_d · X
    /// where q_s = consumption_rate · C/(Ks + C) · X.
    pub fn compute_growth(&self, conc: f64, biomass: f64) -> f64 {
        let mu = monod_kinetics(conc, self.monod_constant, self.max_growth_rate);
        let growth = mu * biomass;
        let decay = self.decay_rate * biomass;
        self.yield_coeff * growth - decay
    }
}

// ---------------------------------------------------------------------------
// BiofilmProperties
// ---------------------------------------------------------------------------

/// Derived macroscopic properties of a biofilm layer.
///
/// Computed from the local biomass density and EPS fraction.
#[derive(Debug, Clone)]
pub struct BiofilmProperties {
    /// Biofilm porosity ε (dimensionless, 0–1).
    pub porosity: f64,
    /// Darcy permeability K (m²) via Kozeny–Carman.
    pub permeability: f64,
    /// Geometric tortuosity τ (dimensionless, ≥ 1).
    pub tortuosity: f64,
}

impl BiofilmProperties {
    /// Estimate properties from `eps_fraction` and `biomass_density`.
    ///
    /// * Porosity: ε = 1 − (eps_fraction + biomass_fraction)  (floored at 0.01)
    /// * Permeability: Kozeny–Carman with `d_particle` = 1 µm
    /// * Tortuosity: τ = 1 / ε^0.5 (geometric estimate)
    pub fn from_cell(cell: &BiofilmCell, d_particle: f64) -> Self {
        let solid_fraction = (cell.eps_fraction + cell.biomass_density / 1000.0).clamp(0.0, 0.99);
        let porosity = (1.0 - solid_fraction).clamp(0.01, 1.0);
        let permeability = kozeny_carman_permeability(porosity, d_particle);
        let tortuosity = 1.0 / porosity.sqrt();
        Self {
            porosity,
            permeability,
            tortuosity,
        }
    }
}

// ---------------------------------------------------------------------------
// BiofilmSimulation (LBM-based 2D)
// ---------------------------------------------------------------------------

/// D2Q9 lattice velocities.
const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 lattice velocities (y component).
const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// D2Q9 lattice weights.
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// Opposite D2Q9 directions for bounce-back.
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// 2-D biofilm growth and nutrient transport simulation.
///
/// Couples an LBM nutrient-transport solver (D2Q9 BGK) with explicit Euler
/// time integration for biomass growth and EPS production.
///
/// # Grid layout
/// Cells are stored in row-major order: index = `iy * nx + ix`.
///
/// # Nutrient LBM
/// The LBM solves a diffusion-reaction equation for the dissolved substrate C.
/// The relaxation time τ is determined from the diffusivity:
///   τ = 3 D + 0.5
/// (in lattice units, with Δx = Δt = 1).
pub struct BiofilmSimulation {
    /// Grid cells (length = `nx * ny`).
    pub grid: Vec<BiofilmCell>,
    /// Number of grid cells in x-direction.
    pub nx: usize,
    /// Number of grid cells in y-direction.
    pub ny: usize,
    /// Nutrient transport parameters.
    pub transport: NutrientTransport,
    /// Biofilm growth parameters.
    pub growth: BiofilmGrowth,
    /// Time step Δt (s).
    pub dt: f64,
    /// LBM distribution functions for nutrient (length = 9 * nx * ny).
    f: Vec<f64>,
    /// Relaxation time τ (lattice units).
    tau: f64,
}

impl BiofilmSimulation {
    /// Create a new [`BiofilmSimulation`] on a `nx × ny` grid.
    pub fn new(
        nx: usize,
        ny: usize,
        transport: NutrientTransport,
        growth: BiofilmGrowth,
        dt: f64,
    ) -> Self {
        let n = nx * ny;
        let tau = 3.0 * transport.diffusivity + 0.5;
        // Initialise distribution functions to equilibrium at C=0
        let f = vec![0.0_f64; 9 * n];
        let grid = vec![BiofilmCell::default(); n];
        Self {
            grid,
            nx,
            ny,
            transport,
            growth,
            dt,
            f,
            tau,
        }
    }

    /// Index for cell (ix, iy).
    #[inline]
    fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Synchronise LBM distribution functions from current nutrient concentrations
    /// (initialise to equilibrium).
    pub fn init_from_grid(&mut self) {
        let n = self.nx * self.ny;
        for i in 0..n {
            let c = self.grid[i].nutrient_concentration;
            for q in 0..9 {
                self.f[q * n + i] = W[q] * c;
            }
        }
    }

    /// Advance the simulation by one time step.
    ///
    /// 1. LBM collision + streaming for nutrient transport.
    /// 2. Explicit-Euler update for biomass growth and nutrient consumption.
    pub fn step(&mut self) {
        self.lbm_step();
        self.biology_step();
    }

    /// LBM collision + streaming for nutrient transport.
    fn lbm_step(&mut self) {
        let n = self.nx * self.ny;
        let tau = self.tau;
        let inv_tau = 1.0 / tau;

        // --- collision ---
        let mut f_out = self.f.clone();
        for i in 0..n {
            let c: f64 = (0..9).map(|q| self.f[q * n + i]).sum();
            for q in 0..9 {
                let feq = W[q] * c;
                f_out[q * n + i] = self.f[q * n + i] - inv_tau * (self.f[q * n + i] - feq);
            }
        }

        // --- reaction sink (added as source term) ---
        for i in 0..n {
            let c = (0..9_usize).map(|q| self.f[q * n + i]).sum::<f64>();
            let biomass = self.grid[i].biomass_density;
            let sink = self.transport.compute_reaction(c, biomass) * self.dt;
            for q in 0..9 {
                f_out[q * n + i] -= W[q] * sink;
            }
        }

        // --- streaming ---
        let mut f_new = vec![0.0_f64; 9 * n];
        let nx = self.nx;
        let ny = self.ny;
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                for q in 0..9 {
                    // Periodic boundaries
                    let nx_i = nx as i32;
                    let ny_i = ny as i32;
                    let jx = (ix as i32 + CX[q]).rem_euclid(nx_i) as usize;
                    let jy = (iy as i32 + CY[q]).rem_euclid(ny_i) as usize;
                    let j = jy * nx + jx;
                    if self.grid[j].is_fluid {
                        f_new[q * n + j] = f_out[q * n + i];
                    } else {
                        // Bounce-back
                        f_new[OPP[q] * n + i] = f_out[q * n + i];
                    }
                }
            }
        }
        self.f = f_new;

        // Update concentrations from distribution sums
        for i in 0..n {
            let c: f64 = (0..9).map(|q| self.f[q * n + i]).sum::<f64>().max(0.0);
            self.grid[i].nutrient_concentration = c;
        }
    }

    /// Explicit-Euler update of biomass and EPS.
    fn biology_step(&mut self) {
        let dt = self.dt;
        let growth = &self.growth;
        for cell in &mut self.grid {
            let conc = cell.nutrient_concentration;
            let biomass = cell.biomass_density;
            let dx = growth.compute_growth(conc, biomass);
            cell.biomass_density = (biomass + dx * dt).max(0.0);
            // EPS production: 10 % of new biomass becomes EPS
            if dx > 0.0 {
                let eps_prod = 0.10 * dx * dt / (biomass.max(1.0));
                cell.eps_fraction = (cell.eps_fraction + eps_prod).clamp(0.0, 1.0);
            }
            // Update fluid/solid classification
            cell.is_fluid = cell.biomass_density < 1.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── monod_kinetics ───────────────────────────────────────────────────────

    #[test]
    fn test_monod_zero_conc() {
        assert_eq!(monod_kinetics(0.0, 0.5, 1.0), 0.0);
    }

    #[test]
    fn test_monod_half_saturation() {
        let rate = monod_kinetics(0.5, 0.5, 1.0);
        assert!((rate - 0.5).abs() < 1e-12, "At Ks, rate = mu_max/2");
    }

    #[test]
    fn test_monod_approaches_mu_max() {
        let rate = monod_kinetics(1e6, 0.5, 1.0);
        assert!((rate - 1.0).abs() < 1e-5, "At high conc, rate → mu_max");
    }

    #[test]
    fn test_monod_zero_ks() {
        assert_eq!(monod_kinetics(1.0, 0.0, 1.0), 0.0);
    }

    #[test]
    fn test_monod_negative_conc() {
        assert_eq!(monod_kinetics(-1.0, 0.5, 1.0), 0.0);
    }

    #[test]
    fn test_monod_proportional_to_mu_max() {
        let r1 = monod_kinetics(1.0, 1.0, 2.0);
        let r2 = monod_kinetics(1.0, 1.0, 4.0);
        assert!((r2 / r1 - 2.0).abs() < 1e-10);
    }

    // ── kozeny_carman_permeability ────────────────────────────────────────────

    #[test]
    fn test_kc_positive_permeability() {
        let k = kozeny_carman_permeability(0.4, 1e-4);
        assert!(k > 0.0);
    }

    #[test]
    fn test_kc_zero_porosity() {
        assert_eq!(kozeny_carman_permeability(0.0, 1e-4), 0.0);
    }

    #[test]
    fn test_kc_porosity_one() {
        assert_eq!(kozeny_carman_permeability(1.0, 1e-4), 0.0);
    }

    #[test]
    fn test_kc_increases_with_porosity() {
        let k1 = kozeny_carman_permeability(0.3, 1e-4);
        let k2 = kozeny_carman_permeability(0.5, 1e-4);
        assert!(k2 > k1, "Higher porosity → higher permeability");
    }

    #[test]
    fn test_kc_scales_with_d_squared() {
        let k1 = kozeny_carman_permeability(0.4, 1e-4);
        let k2 = kozeny_carman_permeability(0.4, 2e-4);
        assert!((k2 / k1 - 4.0).abs() < 1e-8, "K ∝ d²");
    }

    #[test]
    fn test_kc_negative_d() {
        assert_eq!(kozeny_carman_permeability(0.4, -1e-4), 0.0);
    }

    // ── BiofilmCell ──────────────────────────────────────────────────────────

    #[test]
    fn test_cell_fluid_defaults() {
        let c = BiofilmCell::empty_fluid();
        assert!(c.is_fluid);
        assert_eq!(c.biomass_density, 0.0);
        assert_eq!(c.nutrient_concentration, 0.0);
        assert_eq!(c.eps_fraction, 0.0);
    }

    #[test]
    fn test_cell_solid_flag() {
        let c = BiofilmCell::solid(10.0, 1.0, 0.2);
        assert!(!c.is_fluid);
    }

    #[test]
    fn test_cell_eps_clamped() {
        let c = BiofilmCell::fluid(0.0, 0.0, 2.0);
        assert_eq!(c.eps_fraction, 1.0);
        let c2 = BiofilmCell::fluid(0.0, 0.0, -0.5);
        assert_eq!(c2.eps_fraction, 0.0);
    }

    // ── NutrientTransport ─────────────────────────────────────────────────────

    #[test]
    fn test_nutrient_reaction_zero_biomass() {
        let nt = NutrientTransport::new(1e-9, 0.1, 0.5);
        assert_eq!(nt.compute_reaction(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_nutrient_reaction_zero_conc() {
        let nt = NutrientTransport::new(1e-9, 0.1, 0.5);
        assert_eq!(nt.compute_reaction(0.0, 1.0), 0.0);
    }

    #[test]
    fn test_nutrient_reaction_positive() {
        let nt = NutrientTransport::new(1e-9, 1.0, 0.5);
        let r = nt.compute_reaction(0.5, 1.0);
        assert!(r > 0.0, "Reaction rate should be positive");
    }

    #[test]
    fn test_nutrient_reaction_increases_with_biomass() {
        let nt = NutrientTransport::new(1e-9, 1.0, 0.5);
        let r1 = nt.compute_reaction(0.5, 1.0);
        let r2 = nt.compute_reaction(0.5, 2.0);
        assert!(r2 > r1, "Reaction increases with biomass");
    }

    // ── BiofilmGrowth ─────────────────────────────────────────────────────────

    #[test]
    fn test_growth_zero_biomass() {
        let g = BiofilmGrowth::new(0.5, 0.4, 0.05, 0.5);
        assert_eq!(g.compute_growth(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_growth_positive_at_high_conc() {
        let g = BiofilmGrowth::new(0.5, 0.4, 0.05, 0.1);
        let dxdt = g.compute_growth(100.0, 1.0);
        assert!(dxdt > 0.0, "Growth should be positive at high nutrient");
    }

    #[test]
    fn test_growth_decay_dominates_at_low_conc() {
        // Nearly zero nutrient → decay > growth
        let g = BiofilmGrowth::new(0.1, 0.4, 0.5, 1.0);
        let dxdt = g.compute_growth(0.001, 1.0);
        assert!(dxdt < 0.0, "Decay should dominate at near-zero nutrient");
    }

    #[test]
    fn test_growth_increases_with_conc() {
        let g = BiofilmGrowth::new(0.5, 0.5, 0.01, 0.5);
        let g1 = g.compute_growth(0.1, 1.0);
        let g2 = g.compute_growth(1.0, 1.0);
        assert!(g2 > g1, "Growth rate increases with nutrient concentration");
    }

    // ── BiofilmProperties ─────────────────────────────────────────────────────

    #[test]
    fn test_properties_fluid_cell_high_porosity() {
        let cell = BiofilmCell::fluid(0.0, 0.0, 0.0);
        let props = BiofilmProperties::from_cell(&cell, 1e-6);
        assert!(
            props.porosity > 0.9,
            "Empty fluid cell should have high porosity"
        );
    }

    #[test]
    fn test_properties_tortuosity_ge_one() {
        let cell = BiofilmCell::fluid(5.0, 1.0, 0.2);
        let props = BiofilmProperties::from_cell(&cell, 1e-6);
        assert!(props.tortuosity >= 1.0, "Tortuosity should be ≥ 1");
    }

    #[test]
    fn test_properties_permeability_positive() {
        let cell = BiofilmCell::fluid(1.0, 1.0, 0.1);
        let props = BiofilmProperties::from_cell(&cell, 1e-6);
        assert!(props.permeability > 0.0);
    }

    // ── BiofilmSimulation ─────────────────────────────────────────────────────

    #[test]
    fn test_simulation_init() {
        let transport = NutrientTransport::new(0.1, 0.01, 0.5);
        let growth = BiofilmGrowth::new(0.1, 0.5, 0.01, 0.5);
        let sim = BiofilmSimulation::new(4, 4, transport, growth, 0.1);
        assert_eq!(sim.grid.len(), 16);
        assert_eq!(sim.nx, 4);
        assert_eq!(sim.ny, 4);
    }

    #[test]
    fn test_simulation_step_runs() {
        let transport = NutrientTransport::new(0.1, 0.01, 0.5);
        let growth = BiofilmGrowth::new(0.1, 0.5, 0.01, 0.5);
        let mut sim = BiofilmSimulation::new(4, 4, transport, growth, 0.1);
        sim.grid[5].nutrient_concentration = 1.0;
        sim.grid[5].biomass_density = 0.5;
        sim.init_from_grid();
        sim.step(); // Should not panic
    }

    #[test]
    fn test_simulation_nutrient_diffuses() {
        let transport = NutrientTransport::new(0.1, 0.0, 0.5); // no consumption
        let growth = BiofilmGrowth::new(0.0, 0.5, 0.0, 0.5); // no growth
        let mut sim = BiofilmSimulation::new(5, 5, transport, growth, 1.0);
        // Set high nutrient in center cell
        sim.grid[12].nutrient_concentration = 10.0;
        sim.init_from_grid();
        sim.step();
        // After one step, center cell should have less nutrient
        let c_center = sim.grid[12].nutrient_concentration;
        assert!(c_center < 10.0, "Nutrient should diffuse away from center");
    }

    #[test]
    fn test_simulation_biomass_grows_with_nutrient() {
        let transport = NutrientTransport::new(0.1, 0.1, 0.5);
        let growth = BiofilmGrowth::new(1.0, 0.8, 0.01, 0.1);
        let mut sim = BiofilmSimulation::new(4, 4, transport, growth, 0.01);
        sim.grid[0].nutrient_concentration = 10.0;
        sim.grid[0].biomass_density = 0.1;
        sim.init_from_grid();
        let biomass_before = sim.grid[0].biomass_density;
        sim.step();
        let biomass_after = sim.grid[0].biomass_density;
        assert!(
            biomass_after >= biomass_before,
            "Biomass should not decrease at high nutrient"
        );
    }

    #[test]
    fn test_simulation_eps_fraction_bounded() {
        let transport = NutrientTransport::new(0.1, 0.01, 0.5);
        let growth = BiofilmGrowth::new(1.0, 0.8, 0.01, 0.1);
        let mut sim = BiofilmSimulation::new(4, 4, transport, growth, 0.1);
        sim.grid[0].nutrient_concentration = 5.0;
        sim.grid[0].biomass_density = 0.5;
        sim.init_from_grid();
        for _ in 0..10 {
            sim.step();
        }
        for cell in &sim.grid {
            assert!(
                (0.0..=1.0).contains(&cell.eps_fraction),
                "EPS fraction must stay in [0, 1], got {}",
                cell.eps_fraction
            );
        }
    }
}
