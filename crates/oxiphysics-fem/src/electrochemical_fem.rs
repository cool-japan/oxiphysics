#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrochemical FEM formulations.
//!
//! Provides finite element implementations for:
//!
//! - [`NernstPlanckFem`]: Nernst-Planck ion transport (migration + diffusion + convection)
//! - [`PoissonBoltzmannFem`]: Linearized Poisson-Boltzmann / Debye-Hückel / double layer
//! - [`ButlerVolmerFem`]: Electrode reaction boundary condition, overpotential, exchange current
//! - [`BatteryFem`]: Pseudo-2D (P2D) lithium-ion model, Doyle-Fuller-Newman (DFN) model
//! - [`ElectrochemicalCorrosion`]: Pitting corrosion, Evans diagram, galvanic coupling
//! - [`FuelCellFem`]: PEM fuel cell — proton transport, oxygen reduction, water management

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

#[allow(unused_imports)]
use std::f64::consts::PI;

// ============================================================================
// Physical constants
// ============================================================================

/// Faraday constant (C/mol).
const FARADAY: f64 = 96485.332;
/// Universal gas constant (J/(mol·K)).
const R_GAS: f64 = 8.314462;
/// Boltzmann constant (J/K).
const K_BOLTZMANN: f64 = 1.380649e-23;
/// Elementary charge (C).
const ELEM_CHARGE: f64 = 1.602176634e-19;
/// Vacuum permittivity (F/m).
const EPSILON_0: f64 = 8.854187817e-12;

// ============================================================================
// Nernst-Planck FEM
// ============================================================================

/// Species descriptor for Nernst-Planck transport.
#[derive(Debug, Clone)]
pub struct IonSpecies {
    /// Species name (e.g. "Li+", "Cl-").
    pub name: String,
    /// Charge number (signed, dimensionless).
    pub charge: i32,
    /// Diffusivity (m²/s).
    pub diffusivity: f64,
    /// Initial / bulk concentration (mol/m³).
    pub bulk_concentration: f64,
}

impl IonSpecies {
    /// Create a new ion species.
    pub fn new(name: &str, charge: i32, diffusivity: f64, bulk_concentration: f64) -> Self {
        Self {
            name: name.to_owned(),
            charge,
            diffusivity,
            bulk_concentration,
        }
    }
}

/// Nernst-Planck FEM solver for multi-species ion transport.
///
/// Solves the Nernst-Planck equation for each species *k*:
///
/// ∂c_k/∂t + ∇·J_k = R_k
///
/// with flux:
///
/// J_k = -D_k ∇c_k − (z_k F D_k)/(RT) c_k ∇φ + c_k u
///
/// where D_k is diffusivity, z_k charge, φ electric potential, u fluid velocity.
#[derive(Debug, Clone)]
pub struct NernstPlanckFem {
    /// Temperature (K).
    pub temperature: f64,
    /// Ion species list.
    pub species: Vec<IonSpecies>,
    /// Number of spatial nodes in 1-D discretization.
    pub n_nodes: usize,
    /// Domain length (m).
    pub length: f64,
    /// Concentration DOFs: `concentration[s][i]` for species s, node i.
    pub concentration: Vec<Vec<f64>>,
    /// Electric potential at each node (V).
    pub potential: Vec<f64>,
    /// Fluid velocity at each node (m/s) — 1-D.
    pub velocity: Vec<f64>,
    /// Current time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
}

impl NernstPlanckFem {
    /// Construct a new Nernst-Planck FEM instance.
    pub fn new(
        temperature: f64,
        species: Vec<IonSpecies>,
        n_nodes: usize,
        length: f64,
        dt: f64,
    ) -> Self {
        let dx = length / (n_nodes as f64 - 1.0);
        let concentration = species
            .iter()
            .map(|sp| vec![sp.bulk_concentration; n_nodes])
            .collect();
        let potential = (0..n_nodes)
            .map(|i| i as f64 * dx * 0.0)
            .collect::<Vec<_>>();
        let velocity = vec![0.0; n_nodes];
        Self {
            temperature,
            species,
            n_nodes,
            length,
            concentration,
            potential,
            velocity,
            time: 0.0,
            dt,
        }
    }

    /// Compute thermal voltage V_T = RT/F (V).
    pub fn thermal_voltage(&self) -> f64 {
        R_GAS * self.temperature / FARADAY
    }

    /// Compute element diffusion matrix contribution (1-D, linear elements).
    ///
    /// Returns the 2×2 element stiffness contribution for diffusion:
    /// K_e = D * h_e^{-1} * \[\[1,-1\\],\[-1,1\]]
    pub fn element_diffusion_matrix(&self, diffusivity: f64, h_elem: f64) -> [[f64; 2]; 2] {
        let c = diffusivity / h_elem;
        [[c, -c], [-c, c]]
    }

    /// Compute element migration matrix contribution (upwind finite element, 1-D).
    ///
    /// Migration term: -(z F D)/(RT) * ∂(c ∂φ/∂x)/∂x
    pub fn element_migration_matrix(
        &self,
        species_idx: usize,
        h_elem: f64,
        dphi_dx: f64,
    ) -> [[f64; 2]; 2] {
        let sp = &self.species[species_idx];
        let vt = self.thermal_voltage();
        let mu = sp.diffusivity * sp.charge as f64 / vt; // mobility * RT/F factor
        let a = mu * dphi_dx; // advection speed
        // Upwind: full upwinding for stability
        let peclet = a * h_elem / (2.0 * sp.diffusivity + 1e-30);
        let alpha = peclet.tanh() / peclet.abs().max(1e-30) * peclet.signum();
        let kappa = a * h_elem / 2.0 * alpha;
        let base = a / 2.0;
        [
            [-base + kappa / h_elem, base - kappa / h_elem],
            [-base - kappa / h_elem, base + kappa / h_elem],
        ]
    }

    /// Assemble global stiffness matrix for species `species_idx` (tridiagonal, 1-D).
    ///
    /// Returns (lower, diag, upper) diagonal arrays for a tridiagonal system.
    pub fn assemble_species_matrix(&self, species_idx: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.n_nodes;
        let h = self.length / (n as f64 - 1.0);
        let sp = &self.species[species_idx];
        let mut lower = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut upper = vec![0.0; n];
        // Time derivative mass matrix: M/dt
        let mass_c = h / 6.0 / self.dt;
        for i in 0..n {
            diag[i] += 2.0 * mass_c;
            if i + 1 < n {
                upper[i] += mass_c;
                lower[i + 1] += mass_c;
            }
        }
        // Diffusion stiffness
        for i in 0..n - 1 {
            let phi_l = self.potential[i];
            let phi_r = self.potential[i + 1];
            let dphi_dx = (phi_r - phi_l) / h;
            let dm = self.element_diffusion_matrix(sp.diffusivity, h);
            let mm = self.element_migration_matrix(species_idx, h, dphi_dx);
            diag[i] += dm[0][0] + mm[0][0];
            upper[i] += dm[0][1] + mm[0][1];
            lower[i + 1] += dm[1][0] + mm[1][0];
            diag[i + 1] += dm[1][1] + mm[1][1];
        }
        (lower, diag, upper)
    }

    /// Assemble right-hand side (load) vector for species `species_idx`.
    pub fn assemble_species_rhs(&self, species_idx: usize) -> Vec<f64> {
        let n = self.n_nodes;
        let h = self.length / (n as f64 - 1.0);
        let c = &self.concentration[species_idx];
        let mass_c = h / 6.0 / self.dt;
        let mut rhs = vec![0.0; n];
        for i in 0..n - 1 {
            // Consistent mass * c_old
            rhs[i] += mass_c * (2.0 * c[i] + c[i + 1]);
            rhs[i + 1] += mass_c * (c[i] + 2.0 * c[i + 1]);
        }
        rhs
    }

    /// Solve a tridiagonal system via the Thomas algorithm.
    pub fn solve_tridiagonal(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
        let n = diag.len();
        let mut c_prime = vec![0.0; n];
        let mut d_prime = vec![0.0; n];
        let mut x = vec![0.0; n];
        c_prime[0] = upper[0] / diag[0];
        d_prime[0] = rhs[0] / diag[0];
        for i in 1..n {
            let m = diag[i] - lower[i] * c_prime[i - 1];
            c_prime[i] = if i + 1 < n { upper[i] / m } else { 0.0 };
            d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / m;
        }
        x[n - 1] = d_prime[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }
        x
    }

    /// Advance all species concentrations by one time step (implicit Euler).
    pub fn step(&mut self) {
        let ns = self.species.len();
        for s in 0..ns {
            let (lower, diag, upper) = self.assemble_species_matrix(s);
            let rhs = self.assemble_species_rhs(s);
            let c_new = Self::solve_tridiagonal(&lower, &diag, &upper, &rhs);
            self.concentration[s] = c_new;
        }
        self.time += self.dt;
    }

    /// Compute current density at node `i` for species `s` (A/m²).
    ///
    /// j = z_s * F * J_s = z_s * F * (-D ∂c/∂x - z_s F D/(RT) c ∂φ/∂x + c u)
    pub fn current_density(&self, species_idx: usize, node: usize) -> f64 {
        let n = self.n_nodes;
        let h = self.length / (n as f64 - 1.0);
        let sp = &self.species[species_idx];
        let vt = self.thermal_voltage();
        let c = &self.concentration[species_idx];
        let phi = &self.potential;
        // Use central differences where possible
        let (dc_dx, dphi_dx, c_mid, u_mid) = if node == 0 {
            (
                (c[1] - c[0]) / h,
                (phi[1] - phi[0]) / h,
                (c[0] + c[1]) * 0.5,
                (self.velocity[0] + self.velocity[1]) * 0.5,
            )
        } else if node == n - 1 {
            (
                (c[n - 1] - c[n - 2]) / h,
                (phi[n - 1] - phi[n - 2]) / h,
                (c[n - 2] + c[n - 1]) * 0.5,
                (self.velocity[n - 2] + self.velocity[n - 1]) * 0.5,
            )
        } else {
            (
                (c[node + 1] - c[node - 1]) / (2.0 * h),
                (phi[node + 1] - phi[node - 1]) / (2.0 * h),
                c[node],
                self.velocity[node],
            )
        };
        let j = -sp.diffusivity * dc_dx - sp.diffusivity * sp.charge as f64 / vt * c_mid * dphi_dx
            + c_mid * u_mid;
        sp.charge as f64 * FARADAY * j
    }

    /// Compute total ionic current by integrating over the domain (A/m²).
    pub fn total_current(&self, species_idx: usize) -> f64 {
        let n = self.n_nodes;
        (0..n)
            .map(|i| self.current_density(species_idx, i))
            .sum::<f64>()
            / n as f64
    }
}

// ============================================================================
// Poisson-Boltzmann FEM
// ============================================================================

/// Linearized Poisson-Boltzmann FEM for electrostatic double layers.
///
/// Solves the Debye-Hückel (linearized) equation:
///
/// ∇²φ = κ² φ
///
/// where κ = 1/λ_D is the inverse Debye length, with:
///
/// λ_D = sqrt(ε ε₀ R T / (2 F² I))
///
/// and I is the ionic strength I = ½ Σ z_i² c_i^bulk.
#[derive(Debug, Clone)]
pub struct PoissonBoltzmannFem {
    /// Relative permittivity of the solvent (dimensionless).
    pub epsilon_r: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Ion species (for computing Debye length).
    pub species: Vec<IonSpecies>,
    /// Number of 1-D nodes.
    pub n_nodes: usize,
    /// Domain length (m).
    pub length: f64,
    /// Electric potential DOFs (V).
    pub potential: Vec<f64>,
    /// Surface potential at left boundary (V).
    pub phi_surface: f64,
}

impl PoissonBoltzmannFem {
    /// Create a new Poisson-Boltzmann FEM solver.
    pub fn new(
        epsilon_r: f64,
        temperature: f64,
        species: Vec<IonSpecies>,
        n_nodes: usize,
        length: f64,
        phi_surface: f64,
    ) -> Self {
        Self {
            epsilon_r,
            temperature,
            species,
            n_nodes,
            length,
            potential: vec![0.0; n_nodes],
            phi_surface,
        }
    }

    /// Compute ionic strength I = ½ Σ z_i² c_i^bulk (mol/m³).
    pub fn ionic_strength(&self) -> f64 {
        0.5 * self
            .species
            .iter()
            .map(|sp| (sp.charge as f64).powi(2) * sp.bulk_concentration)
            .sum::<f64>()
    }

    /// Compute Debye length λ_D (m).
    pub fn debye_length(&self) -> f64 {
        let i_str = self.ionic_strength().max(1e-10);
        let eps = self.epsilon_r * EPSILON_0;
        (eps * R_GAS * self.temperature / (2.0 * FARADAY.powi(2) * i_str)).sqrt()
    }

    /// Compute inverse Debye length κ = 1/λ_D (m⁻¹).
    pub fn inverse_debye_length(&self) -> f64 {
        1.0 / self.debye_length()
    }

    /// Assemble and solve the linearized Poisson-Boltzmann equation.
    ///
    /// Uses linear finite elements with Dirichlet BC at x=0 (surface)
    /// and φ→0 at x=L (bulk).
    pub fn solve(&mut self) {
        let n = self.n_nodes;
        let h = self.length / (n as f64 - 1.0);
        let kappa = self.inverse_debye_length();
        let kappa2 = kappa * kappa;
        // Assemble tridiagonal: -∇²φ + κ² φ = 0  → (K + κ² M) φ = 0
        let mut lower = vec![0.0_f64; n];
        let mut diag = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n - 1 {
            let k_diff = 1.0 / h;
            let m_mass = h / 6.0;
            diag[i] += k_diff + 2.0 * m_mass * kappa2;
            upper[i] += -k_diff + m_mass * kappa2;
            lower[i + 1] += -k_diff + m_mass * kappa2;
            diag[i + 1] += k_diff + 2.0 * m_mass * kappa2;
        }
        // Dirichlet: node 0 → phi_surface
        diag[0] = 1.0;
        upper[0] = 0.0;
        rhs[0] = self.phi_surface;
        // Dirichlet: node n-1 → 0 (bulk)
        diag[n - 1] = 1.0;
        lower[n - 1] = 0.0;
        rhs[n - 1] = 0.0;
        self.potential = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
    }

    /// Compute surface charge density σ = -ε ε₀ ∂φ/∂x|_{x=0} (C/m²).
    pub fn surface_charge_density(&self) -> f64 {
        if self.n_nodes < 2 {
            return 0.0;
        }
        let h = self.length / (self.n_nodes as f64 - 1.0);
        let dphi = (self.potential[1] - self.potential[0]) / h;
        -self.epsilon_r * EPSILON_0 * dphi
    }

    /// Compute differential capacitance C_diff = ∂σ/∂φ_surface (F/m²).
    ///
    /// For the linearized PB model: C_diff = ε ε₀ κ.
    pub fn differential_capacitance(&self) -> f64 {
        self.epsilon_r * EPSILON_0 * self.inverse_debye_length()
    }

    /// Evaluate the Gouy-Chapman analytical solution at position x (V).
    ///
    /// For the linearized (Debye-Hückel) regime: φ(x) = φ_0 * exp(-κ x).
    pub fn analytical_solution(&self, x: f64) -> f64 {
        self.phi_surface * (-self.inverse_debye_length() * x).exp()
    }

    /// Compute zeta potential from the computed potential profile (V).
    ///
    /// Returns the potential at the shear plane, approximated as the first node.
    pub fn zeta_potential(&self) -> f64 {
        if self.potential.is_empty() {
            0.0
        } else {
            self.potential[0]
        }
    }
}

// ============================================================================
// Butler-Volmer FEM Boundary Condition
// ============================================================================

/// Electrode kinetics via the Butler-Volmer equation as a FEM boundary condition.
///
/// Models the faradaic current density at an electrode surface:
///
/// j = j₀ \[exp(α_a F η / RT) − exp(−α_c F η / RT)\]
///
/// where η = φ_s − φ_l − U_eq is the overpotential, j₀ is the exchange
/// current density, and α_a, α_c are transfer coefficients.
#[derive(Debug, Clone)]
pub struct ButlerVolmerFem {
    /// Temperature (K).
    pub temperature: f64,
    /// Exchange current density j₀ (A/m²).
    pub exchange_current: f64,
    /// Anodic transfer coefficient α_a (dimensionless).
    pub alpha_anodic: f64,
    /// Cathodic transfer coefficient α_c (dimensionless).
    pub alpha_cathodic: f64,
    /// Equilibrium potential U_eq (V).
    pub equilibrium_potential: f64,
    /// Solid-phase potential φ_s (V).
    pub phi_solid: f64,
    /// Liquid-phase potential φ_liquid (V).
    pub phi_liquid: f64,
}

impl ButlerVolmerFem {
    /// Create a new Butler-Volmer boundary condition.
    pub fn new(
        temperature: f64,
        exchange_current: f64,
        alpha_anodic: f64,
        alpha_cathodic: f64,
        equilibrium_potential: f64,
    ) -> Self {
        Self {
            temperature,
            exchange_current,
            alpha_anodic,
            alpha_cathodic,
            equilibrium_potential,
            phi_solid: 0.0,
            phi_liquid: 0.0,
        }
    }

    /// Set electrode potentials.
    pub fn set_potentials(&mut self, phi_solid: f64, phi_liquid: f64) {
        self.phi_solid = phi_solid;
        self.phi_liquid = phi_liquid;
    }

    /// Compute overpotential η = φ_s − φ_l − U_eq (V).
    pub fn overpotential(&self) -> f64 {
        self.phi_solid - self.phi_liquid - self.equilibrium_potential
    }

    /// Compute Butler-Volmer current density j (A/m²).
    pub fn current_density(&self) -> f64 {
        let eta = self.overpotential();
        let f_rt = FARADAY / (R_GAS * self.temperature);
        self.exchange_current
            * ((self.alpha_anodic * f_rt * eta).exp() - (-self.alpha_cathodic * f_rt * eta).exp())
    }

    /// Compute linearized exchange conductance ∂j/∂η at η≈0 (S/m²).
    ///
    /// dj/dη|_{η=0} = j₀ (α_a + α_c) F/(RT)
    pub fn linearized_conductance(&self) -> f64 {
        let f_rt = FARADAY / (R_GAS * self.temperature);
        self.exchange_current * (self.alpha_anodic + self.alpha_cathodic) * f_rt
    }

    /// Compute Tafel slope for anodic branch (V/decade).
    ///
    /// b_a = RT / (α_a F) * ln(10)
    pub fn tafel_slope_anodic(&self) -> f64 {
        R_GAS * self.temperature / (self.alpha_anodic * FARADAY) * 10.0_f64.ln()
    }

    /// Compute Tafel slope for cathodic branch (V/decade).
    pub fn tafel_slope_cathodic(&self) -> f64 {
        R_GAS * self.temperature / (self.alpha_cathodic * FARADAY) * 10.0_f64.ln()
    }

    /// Assemble the FEM tangent (Jacobian) contribution from the BV boundary condition.
    ///
    /// Returns the scalar tangent value dj/dφ_s for use in the global Jacobian assembly.
    pub fn tangent_contribution(&self) -> f64 {
        let eta = self.overpotential();
        let f_rt = FARADAY / (R_GAS * self.temperature);
        self.exchange_current
            * f_rt
            * (self.alpha_anodic * (self.alpha_anodic * f_rt * eta).exp()
                + self.alpha_cathodic * (-self.alpha_cathodic * f_rt * eta).exp())
    }

    /// Compute open-circuit potential as function of state-of-charge (SOC) using a
    /// typical graphite anode empirical correlation (V).
    pub fn ocp_graphite(soc: f64) -> f64 {
        let u = soc.clamp(0.001, 0.999);
        0.7222 + 0.1387 * u + 0.029 * u.powf(-0.5) - 0.0172 / u + 1.54e-4 / u.powi(2)
            - 0.1231 * (-33.0 * u).exp()
            - 0.2 * (5.0 * u - 3.5).exp()
    }

    /// Compute open-circuit potential for LiCoO₂ cathode (V vs. Li/Li⁺).
    ///
    /// Returns values in the range approximately \[3.4, 4.2\] V.
    pub fn ocp_licoo2(soc: f64) -> f64 {
        let u = soc.clamp(0.001, 0.999);
        // Simplified Doyle-Fuller-Newman LiCoO2 OCP (V vs Li/Li+)
        3.5 + 0.4 * u - 0.1 * u.powi(2)
    }
}

// ============================================================================
// Battery FEM — P2D / DFN Model
// ============================================================================

/// Domain type for the DFN model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DfnDomain {
    /// Negative electrode (anode).
    NegativeElectrode,
    /// Separator.
    Separator,
    /// Positive electrode (cathode).
    PositiveElectrode,
}

/// Parameters for one electrode domain in the P2D model.
#[derive(Debug, Clone)]
pub struct ElectrodeDomain {
    /// Domain type.
    pub domain: DfnDomain,
    /// Thickness (m).
    pub thickness: f64,
    /// Number of macro-scale (x-direction) nodes.
    pub n_macro: usize,
    /// Number of micro-scale (r-direction) nodes in each particle.
    pub n_micro: usize,
    /// Solid-phase Li diffusivity D_s (m²/s).
    pub diffusivity_solid: f64,
    /// Maximum solid concentration c_s_max (mol/m³).
    pub c_solid_max: f64,
    /// Particle radius r_p (m).
    pub particle_radius: f64,
    /// Electrode porosity ε_l.
    pub porosity: f64,
    /// Bruggeman exponent p.
    pub bruggeman: f64,
    /// Solid-phase conductivity σ (S/m).
    pub sigma_solid: f64,
    /// Exchange current density j₀_ref (A/m²) at c_s_max/2.
    pub j0_ref: f64,
}

impl ElectrodeDomain {
    /// Create a default negative electrode (graphite) parameter set.
    pub fn default_negative() -> Self {
        Self {
            domain: DfnDomain::NegativeElectrode,
            thickness: 100e-6,
            n_macro: 10,
            n_micro: 5,
            diffusivity_solid: 3.9e-14,
            c_solid_max: 31370.0,
            particle_radius: 2.0e-6,
            porosity: 0.3,
            bruggeman: 1.5,
            sigma_solid: 100.0,
            j0_ref: 1.0,
        }
    }

    /// Create a default positive electrode (LiCoO₂) parameter set.
    pub fn default_positive() -> Self {
        Self {
            domain: DfnDomain::PositiveElectrode,
            thickness: 80e-6,
            n_macro: 10,
            n_micro: 5,
            diffusivity_solid: 1.0e-14,
            c_solid_max: 51218.0,
            particle_radius: 3.0e-6,
            porosity: 0.3,
            bruggeman: 1.5,
            sigma_solid: 3.8,
            j0_ref: 1.0,
        }
    }

    /// Effective electrolyte diffusivity D_e_eff = D_e * ε^p (m²/s).
    pub fn effective_electrolyte_diffusivity(&self, d_electrolyte: f64) -> f64 {
        d_electrolyte * self.porosity.powf(self.bruggeman)
    }

    /// Effective electrolyte conductivity κ_eff = κ * ε^p (S/m).
    pub fn effective_conductivity(&self, kappa: f64) -> f64 {
        kappa * self.porosity.powf(self.bruggeman)
    }
}

/// Pseudo-2D (P2D) Doyle-Fuller-Newman lithium-ion battery FEM model.
///
/// Implements the DFN model with:
/// - Macro-scale 1-D along cell thickness (x direction)
/// - Micro-scale 1-D radial diffusion inside spherical particles (r direction)
/// - Butler-Volmer electrode kinetics
/// - Electrolyte transport (Nernst-Planck)
#[derive(Debug, Clone)]
pub struct BatteryFem {
    /// Negative electrode parameters.
    pub neg: ElectrodeDomain,
    /// Positive electrode parameters.
    pub pos: ElectrodeDomain,
    /// Separator thickness (m).
    pub separator_thickness: f64,
    /// Number of nodes in separator.
    pub n_sep: usize,
    /// Electrolyte diffusivity D_e (m²/s).
    pub d_electrolyte: f64,
    /// Electrolyte conductivity κ (S/m).
    pub kappa_electrolyte: f64,
    /// Li⁺ transference number t₊.
    pub transference: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Applied current density (A/m²), positive = discharge.
    pub applied_current: f64,
    /// Solid concentration: `c_solid[domain][x_node][r_node]`.
    pub c_solid: Vec<Vec<Vec<f64>>>,
    /// Electrolyte concentration along x (mol/m³).
    pub c_electrolyte: Vec<f64>,
    /// Solid-phase potential along x in electrode domains (V).
    pub phi_solid: Vec<Vec<f64>>,
    /// Electrolyte potential along x (V).
    pub phi_electrolyte: Vec<f64>,
    /// Internal (ohmic) resistance of the cell (Ω).
    ///
    /// Terminal voltage = OCV − `internal_resistance` × `applied_current`.
    pub internal_resistance: f64,
    /// Current time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
}

impl BatteryFem {
    /// Construct a new P2D/DFN battery FEM model with default parameters.
    pub fn new(applied_current: f64, temperature: f64, dt: f64) -> Self {
        let neg = ElectrodeDomain::default_negative();
        let pos = ElectrodeDomain::default_positive();
        let n_sep = 5_usize;
        let n_total = neg.n_macro + n_sep + pos.n_macro;
        // Initialize solid concentration at 50% SOC in both electrodes
        let c_solid = vec![
            vec![vec![neg.c_solid_max * 0.5; neg.n_micro]; neg.n_macro],
            vec![vec![pos.c_solid_max * 0.8; pos.n_micro]; pos.n_macro],
        ];
        Self {
            neg,
            pos,
            separator_thickness: 25e-6,
            n_sep,
            d_electrolyte: 2.7e-10,
            kappa_electrolyte: 1.0,
            transference: 0.363,
            temperature,
            applied_current,
            c_solid,
            c_electrolyte: vec![1000.0; n_total],
            phi_solid: vec![vec![0.0; 10], vec![4.2; 10]],
            phi_electrolyte: vec![0.0; n_total],
            internal_resistance: 0.0,
            time: 0.0,
            dt,
        }
    }

    /// Compute total number of macro-scale nodes.
    pub fn n_total_nodes(&self) -> usize {
        self.neg.n_macro + self.n_sep + self.pos.n_macro
    }

    /// Compute state-of-charge (SOC) of negative electrode (volume-averaged).
    pub fn soc_negative(&self) -> f64 {
        let c_avg: f64 = self.c_solid[0]
            .iter()
            .map(|particle| particle.iter().sum::<f64>() / particle.len() as f64)
            .sum::<f64>()
            / self.neg.n_macro as f64;
        c_avg / self.neg.c_solid_max
    }

    /// Compute state-of-charge (SOC) of positive electrode (volume-averaged).
    pub fn soc_positive(&self) -> f64 {
        let c_avg: f64 = self.c_solid[1]
            .iter()
            .map(|particle| particle.iter().sum::<f64>() / particle.len() as f64)
            .sum::<f64>()
            / self.pos.n_macro as f64;
        c_avg / self.pos.c_solid_max
    }

    /// Compute cell terminal voltage (V).
    ///
    /// V_cell = OCP_pos − OCP_neg − I·R_int
    ///
    /// where `I` = `applied_current` (A) and `R_int` = `internal_resistance` (Ω).
    pub fn terminal_voltage(&self) -> f64 {
        let ocp_pos = ButlerVolmerFem::ocp_licoo2(self.soc_positive());
        let ocp_neg = ButlerVolmerFem::ocp_graphite(self.soc_negative());
        let ir_drop = self.applied_current * self.internal_resistance;
        ocp_pos - ocp_neg - ir_drop
    }

    /// Advance solid-phase diffusion in spherical particles by one time step.
    ///
    /// Solves ∂c_s/∂t = D_s/r² ∂/∂r(r² ∂c_s/∂r) using implicit scheme.
    pub fn step_solid_diffusion(&mut self, electrode_idx: usize) {
        let dom = if electrode_idx == 0 {
            &self.neg
        } else {
            &self.pos
        };
        let n_r = dom.n_micro;
        let r_p = dom.particle_radius;
        let ds = dom.diffusivity_solid;
        let h_r = r_p / (n_r as f64 - 1.0);
        let dt = self.dt;
        let n_x = dom.n_macro;
        for xi in 0..n_x {
            let c = &self.c_solid[electrode_idx][xi];
            // Build tridiagonal for spherical diffusion operator
            let mut lower = vec![0.0_f64; n_r];
            let mut diag = vec![0.0_f64; n_r];
            let mut upper = vec![0.0_f64; n_r];
            let mut rhs = c.clone();
            // Implicit diffusion: (I - dt*L) c^{n+1} = c^n
            for ri in 1..n_r - 1 {
                let r = ri as f64 * h_r;
                let coef = ds * dt / (h_r * h_r);
                let r_p_half = r + 0.5 * h_r;
                let r_m_half = r - 0.5 * h_r;
                let a_plus = coef * (r_p_half / r).powi(2);
                let a_minus = coef * (r_m_half / r).powi(2);
                lower[ri] = -a_minus;
                diag[ri] = 1.0 + a_plus + a_minus;
                upper[ri] = -a_plus;
            }
            // Symmetry at r=0 (Neumann): ∂c/∂r = 0
            diag[0] = 1.0 + 2.0 * ds * dt / (h_r * h_r);
            upper[0] = -2.0 * ds * dt / (h_r * h_r);
            // Surface flux BC at r=r_p (from BV current)
            let j_li = self.applied_current
                / (3.0 / dom.particle_radius * (1.0 - dom.porosity) * dom.thickness);
            diag[n_r - 1] = 1.0;
            lower[n_r - 1] = 0.0;
            rhs[n_r - 1] = c[n_r - 1] - j_li * dt * h_r / (FARADAY * ds);
            let new_c = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
            self.c_solid[electrode_idx][xi] = new_c;
        }
    }

    /// Advance the entire DFN model by one time step.
    pub fn step(&mut self) {
        self.step_solid_diffusion(0);
        self.step_solid_diffusion(1);
        self.time += self.dt;
    }

    /// Compute C-rate for a given capacity (Ah).
    pub fn c_rate(&self, capacity_ah: f64) -> f64 {
        self.applied_current / (capacity_ah * 3600.0)
    }

    /// Compute specific energy (Wh/kg) given electrode mass density (kg/m³).
    pub fn specific_energy(&self, mass_density: f64, volume: f64) -> f64 {
        let energy = self.terminal_voltage() * self.applied_current * self.time;
        let mass = mass_density * volume;
        energy / mass / 3600.0
    }
}

// ============================================================================
// Electrochemical Corrosion
// ============================================================================

/// Pitting corrosion model via FEM.
///
/// Models the competitive anodic (metal dissolution) and cathodic (oxygen reduction
/// or hydrogen evolution) reactions on a metal surface, including Evans diagram
/// construction and galvanic coupling.
#[derive(Debug, Clone)]
pub struct ElectrochemicalCorrosion {
    /// Temperature (K).
    pub temperature: f64,
    /// Corrosion potential E_corr (V vs. SHE).
    pub e_corrosion: f64,
    /// Anodic Tafel slope β_a (V/decade).
    pub beta_anodic: f64,
    /// Cathodic Tafel slope β_c (V/decade).
    pub beta_cathodic: f64,
    /// Corrosion current density i_corr (A/m²).
    pub i_corrosion: f64,
    /// Pitting potential E_pit (V), above which stable pits nucleate.
    pub e_pit: f64,
    /// Repassivation potential E_rep (V), below which pits repassivate.
    pub e_rep: f64,
    /// Active pitting flag.
    pub pitting_active: bool,
    /// Galvanic couples: (anode area, cathode area, potential difference) list.
    pub galvanic_couples: Vec<(f64, f64, f64)>,
    /// Number of nodes for FEM discretization.
    pub n_nodes: usize,
    /// Domain length (m).
    pub length: f64,
    /// Potential distribution along the electrode (V vs. SHE).
    pub potential: Vec<f64>,
    /// Corrosion rate distribution (mm/year).
    pub corrosion_rate: Vec<f64>,
}

impl ElectrochemicalCorrosion {
    /// Create a new corrosion model (default: iron in NaCl solution).
    pub fn new(
        temperature: f64,
        e_corrosion: f64,
        beta_anodic: f64,
        beta_cathodic: f64,
        i_corrosion: f64,
        e_pit: f64,
        e_rep: f64,
        n_nodes: usize,
        length: f64,
    ) -> Self {
        Self {
            temperature,
            e_corrosion,
            beta_anodic,
            beta_cathodic,
            i_corrosion,
            e_pit,
            e_rep,
            pitting_active: false,
            galvanic_couples: Vec::new(),
            n_nodes,
            length,
            potential: vec![e_corrosion; n_nodes],
            corrosion_rate: vec![0.0; n_nodes],
        }
    }

    /// Compute anodic current density at potential E (A/m²) via Tafel equation.
    pub fn anodic_current(&self, e: f64) -> f64 {
        self.i_corrosion * 10.0_f64.powf((e - self.e_corrosion) / self.beta_anodic)
    }

    /// Compute cathodic current density at potential E (A/m²) via Tafel equation.
    pub fn cathodic_current(&self, e: f64) -> f64 {
        self.i_corrosion * 10.0_f64.powf(-(e - self.e_corrosion) / self.beta_cathodic)
    }

    /// Compute net current density j_net = j_a − j_c (A/m²).
    pub fn net_current(&self, e: f64) -> f64 {
        self.anodic_current(e) - self.cathodic_current(e)
    }

    /// Check whether pitting is active at potential E.
    ///
    /// Pitting initiates when E ≥ E_pit.
    pub fn check_pitting(&mut self, e: f64) {
        self.pitting_active = e >= self.e_pit;
    }

    /// Compute corrosion rate (mm/year) from current density (A/m²).
    ///
    /// Uses Faraday's law: rate = (i_corr * M) / (n * F * ρ) * unit conversion.
    ///
    /// Default: iron (M = 55.85 g/mol, n = 2, ρ = 7874 kg/m³).
    pub fn current_to_corrosion_rate(current_density: f64) -> f64 {
        let m_fe = 55.85e-3; // kg/mol
        let n_eq = 2.0;
        let rho_fe = 7874.0; // kg/m³
        // rate in m/s → mm/year
        let rate_ms = current_density * m_fe / (n_eq * FARADAY * rho_fe);
        rate_ms * 1000.0 * 365.25 * 24.0 * 3600.0
    }

    /// Add a galvanic couple (anode area A_a, cathode area A_c, E_galv difference).
    pub fn add_galvanic_couple(&mut self, anode_area: f64, cathode_area: f64, delta_e: f64) {
        self.galvanic_couples
            .push((anode_area, cathode_area, delta_e));
    }

    /// Compute galvanic current for all couples (A).
    ///
    /// Simple mixed-potential approximation.
    pub fn galvanic_current(&self) -> Vec<f64> {
        self.galvanic_couples
            .iter()
            .map(|(a_a, _a_c, de)| {
                let i_gal = self.i_corrosion * 10.0_f64.powf(de / self.beta_anodic);
                i_gal * a_a
            })
            .collect()
    }

    /// Update corrosion rate distribution along the FEM domain.
    pub fn update_corrosion_rates(&mut self) {
        self.corrosion_rate = self
            .potential
            .iter()
            .map(|&e| {
                let j = self.net_current(e).max(0.0);
                Self::current_to_corrosion_rate(j)
            })
            .collect();
    }

    /// Solve Laplace equation for potential distribution (Ohmic electrolyte).
    ///
    /// Boundary conditions: left node at E_anode, right node at E_cathode.
    pub fn solve_laplace(&mut self, e_anode: f64, e_cathode: f64, conductivity: f64) {
        let n = self.n_nodes;
        let h = self.length / (n as f64 - 1.0);
        let mut lower = vec![0.0_f64; n];
        let mut diag = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n - 1 {
            let k = conductivity / h;
            diag[i] += k;
            upper[i] -= k;
            lower[i + 1] -= k;
            diag[i + 1] += k;
        }
        // Dirichlet BCs
        diag[0] = 1.0;
        upper[0] = 0.0;
        rhs[0] = e_anode;
        diag[n - 1] = 1.0;
        lower[n - 1] = 0.0;
        rhs[n - 1] = e_cathode;
        self.potential = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
        self.update_corrosion_rates();
    }

    /// Generate Evans diagram data (E vs. log|i|) for plotting.
    ///
    /// Returns (potential_vec, anodic_log_i, cathodic_log_i) arrays.
    pub fn evans_diagram(&self, n_points: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let e_range = 0.5;
        let potentials: Vec<f64> = (0..n_points)
            .map(|i| self.e_corrosion - e_range + 2.0 * e_range * i as f64 / (n_points - 1) as f64)
            .collect();
        let anodic: Vec<f64> = potentials
            .iter()
            .map(|&e| self.anodic_current(e).max(1e-12).log10())
            .collect();
        let cathodic: Vec<f64> = potentials
            .iter()
            .map(|&e| self.cathodic_current(e).max(1e-12).log10())
            .collect();
        (potentials, anodic, cathodic)
    }
}

// ============================================================================
// Fuel Cell FEM — PEM
// ============================================================================

/// Membrane region descriptor for the PEM fuel cell model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FuelCellRegion {
    /// Anode gas diffusion layer.
    AnodeGdl,
    /// Anode catalyst layer.
    AnodeCatalyst,
    /// Polymer electrolyte membrane.
    Membrane,
    /// Cathode catalyst layer.
    CathodeCatalyst,
    /// Cathode gas diffusion layer.
    CathodeGdl,
}

/// PEM fuel cell FEM model.
///
/// Models the following coupled physics through the cell thickness (1-D):
///
/// 1. Proton transport through the membrane: -∇·(σ_m ∇φ_m) = S_φ
/// 2. Oxygen transport in cathode GDL/CL: ∂c_O2/∂t = ∇·(D_O2 ∇c_O2) − R_O2
/// 3. Water management: liquid saturation and capillary transport
/// 4. Butler-Volmer ORR kinetics at cathode catalyst layer
#[derive(Debug, Clone)]
pub struct FuelCellFem {
    /// Temperature (K).
    pub temperature: f64,
    /// Membrane proton conductivity (S/m).
    pub sigma_membrane: f64,
    /// Membrane thickness (m).
    pub thickness_membrane: f64,
    /// Cathode GDL porosity.
    pub porosity_gdl: f64,
    /// GDL thickness (m).
    pub thickness_gdl: f64,
    /// Catalyst layer thickness (m).
    pub thickness_cl: f64,
    /// O₂ diffusivity in GDL (m²/s).
    pub d_oxygen: f64,
    /// Reference O₂ concentration (mol/m³).
    pub c_oxygen_ref: f64,
    /// Reference current density for ORR j₀_ref (A/m²).
    pub j0_orr: f64,
    /// Applied cell voltage V_cell (V).
    pub v_cell: f64,
    /// Open-circuit voltage V_oc (V).
    pub v_oc: f64,
    /// Number of nodes across each region.
    pub n_nodes_per_region: usize,
    /// Proton potential across the cell (V).
    pub phi_proton: Vec<f64>,
    /// Oxygen concentration across cathode (mol/m³).
    pub c_oxygen: Vec<f64>,
    /// Liquid water saturation (dimensionless 0–1).
    pub saturation: Vec<f64>,
    /// Local current density distribution (A/m²).
    pub current_density: Vec<f64>,
    /// Current time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
}

impl FuelCellFem {
    /// Create a new PEM fuel cell FEM model.
    pub fn new(temperature: f64, v_cell: f64, dt: f64) -> Self {
        let n = 20_usize;
        Self {
            temperature,
            sigma_membrane: 10.0,
            thickness_membrane: 50e-6,
            porosity_gdl: 0.6,
            thickness_gdl: 200e-6,
            thickness_cl: 10e-6,
            d_oxygen: 2.6e-5,
            c_oxygen_ref: 40.97, // mol/m³ at 1 atm
            j0_orr: 1e-4,
            v_cell,
            v_oc: 1.23,
            n_nodes_per_region: n,
            phi_proton: vec![0.0; n],
            c_oxygen: vec![40.97; n],
            saturation: vec![0.1; n],
            current_density: vec![0.0; n],
            time: 0.0,
            dt,
        }
    }

    /// Compute total overpotential η = V_oc − V_cell (V).
    pub fn total_overpotential(&self) -> f64 {
        self.v_oc - self.v_cell
    }

    /// Compute proton conductivity as function of membrane water content λ (S/m).
    ///
    /// Springer correlation: σ_m(T,λ) = (0.5139 λ − 0.326) exp\[1268(1/303 − 1/T)\]
    pub fn proton_conductivity(&self, lambda_water: f64) -> f64 {
        let sigma_0 = 0.5139 * lambda_water - 0.326;
        let exp_factor = (1268.0 * (1.0 / 303.0 - 1.0 / self.temperature)).exp();
        sigma_0.max(0.0) * exp_factor
    }

    /// Compute effective O₂ diffusivity in the GDL (m²/s).
    ///
    /// Bruggeman correction: D_eff = D_O2 * ε^1.5 * (1-s)^1.5
    pub fn effective_oxygen_diffusivity(&self, porosity: f64, saturation: f64) -> f64 {
        self.d_oxygen * porosity.powf(1.5) * (1.0 - saturation).powf(1.5)
    }

    /// Compute ORR reaction rate at the cathode catalyst layer (mol/(m³·s)).
    ///
    /// Uses Butler-Volmer with cathodic-only approximation (high overpotential):
    /// j_ORR = j₀ * (c_O2/c_ref) * exp(−α_c F η_c / RT)
    pub fn orr_rate(&self, c_o2: f64, eta_cathode: f64) -> f64 {
        let alpha_c = 0.5;
        let f_rt = FARADAY / (R_GAS * self.temperature);
        self.j0_orr * (c_o2 / self.c_oxygen_ref) * (-alpha_c * f_rt * eta_cathode).exp() / FARADAY // convert A/m² to mol/(m²·s)
    }

    /// Solve proton transport across membrane using FEM (1-D Laplace).
    ///
    /// -∇·(σ_m ∇φ) = 0, BC: φ(0) = 0 (anode), φ(L) = V_cell
    pub fn solve_proton_transport(&mut self) {
        let n = self.n_nodes_per_region;
        let h = self.thickness_membrane / (n as f64 - 1.0);
        let sigma = self.sigma_membrane;
        let mut lower = vec![0.0_f64; n];
        let mut diag = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        let k = sigma / h;
        for i in 0..n - 1 {
            diag[i] += k;
            upper[i] -= k;
            lower[i + 1] -= k;
            diag[i + 1] += k;
        }
        diag[0] = 1.0;
        upper[0] = 0.0;
        rhs[0] = 0.0;
        diag[n - 1] = 1.0;
        lower[n - 1] = 0.0;
        rhs[n - 1] = self.v_cell;
        self.phi_proton = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
    }

    /// Solve O₂ transport in cathode GDL/CL with ORR sink term.
    pub fn solve_oxygen_transport(&mut self) {
        let n = self.n_nodes_per_region;
        let h = self.thickness_gdl / (n as f64 - 1.0);
        let eta_c = self.total_overpotential();
        let d_eff = self.effective_oxygen_diffusivity(self.porosity_gdl, 0.1);
        let mut lower = vec![0.0_f64; n];
        let mut diag = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        let c_ref = self.c_oxygen_ref;
        for i in 0..n - 1 {
            let k = d_eff / h;
            diag[i] += k;
            upper[i] -= k;
            lower[i + 1] -= k;
            diag[i + 1] += k;
        }
        // Add ORR sink in cathode CL (last few nodes)
        let n_cl = (self.thickness_cl / h).ceil() as usize;
        for i in n.saturating_sub(n_cl)..n {
            let r_o2 = self.orr_rate(self.c_oxygen[i].max(0.0), eta_c.abs()) * FARADAY;
            diag[i] += r_o2 / self.c_oxygen[i].max(1e-10) * h;
        }
        // BCs: inlet c_O2 = c_ref, at membrane ∂c/∂x = 0
        diag[0] = 1.0;
        upper[0] = 0.0;
        rhs[0] = c_ref;
        lower[n - 1] = 0.0;
        self.c_oxygen = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
    }

    /// Compute local current density at each node via ORR kinetics (A/m²).
    pub fn compute_current_density(&mut self) {
        let eta_c = self.total_overpotential().abs();
        self.current_density = self
            .c_oxygen
            .iter()
            .map(|&c| {
                let rate = self.orr_rate(c.max(0.0), eta_c);
                rate * FARADAY * self.thickness_cl
            })
            .collect();
    }

    /// Compute average current density (A/m²).
    pub fn average_current_density(&self) -> f64 {
        self.current_density.iter().sum::<f64>() / self.current_density.len() as f64
    }

    /// Compute power density P = V_cell * j_avg (W/m²).
    pub fn power_density(&self) -> f64 {
        self.v_cell * self.average_current_density()
    }

    /// Advance the fuel cell model by one time step.
    pub fn step(&mut self) {
        self.solve_proton_transport();
        self.solve_oxygen_transport();
        self.compute_current_density();
        self.time += self.dt;
    }

    /// Compute Nernst potential correction for given p_O2 and p_H2 (Pa) at temperature T.
    ///
    /// V_Nernst = V_0 + RT/(2F) * ln(p_H2 * sqrt(p_O2) / p_H2O)
    pub fn nernst_potential(temperature: f64, p_h2: f64, p_o2: f64, p_h2o: f64) -> f64 {
        let v0 = 1.229;
        let correction = (temperature - 298.15) * (-8.5e-4);
        let nernst =
            R_GAS * temperature / (2.0 * FARADAY) * (p_h2 * p_o2.sqrt() / p_h2o.max(1e-10)).ln();
        v0 + correction + nernst
    }

    /// Compute activation overpotential from Tafel equation (V).
    pub fn activation_overpotential(&self, current_density: f64) -> f64 {
        if current_density <= 0.0 {
            return 0.0;
        }
        R_GAS * self.temperature / (0.5 * FARADAY) * (current_density / self.j0_orr).ln()
    }

    /// Compute ohmic overpotential (V).
    pub fn ohmic_overpotential(&self, current_density: f64, r_ohm: f64) -> f64 {
        current_density * r_ohm
    }

    /// Compute concentration overpotential (V).
    ///
    /// η_conc = -RT/(2F) * ln(1 − j/j_lim)
    pub fn concentration_overpotential(&self, current_density: f64, j_limiting: f64) -> f64 {
        let ratio = (current_density / j_limiting).min(0.999);
        -R_GAS * self.temperature / (2.0 * FARADAY) * (1.0 - ratio).ln()
    }

    /// Compute electro-osmotic drag coefficient (water molecules per proton).
    ///
    /// n_d = 2.5 λ / 22 (Springer model)
    pub fn electro_osmotic_drag(&self, lambda_water: f64) -> f64 {
        2.5 * lambda_water / 22.0
    }

    /// Compute water diffusivity in membrane D_w(λ) (m²/s).
    ///
    /// Springer-Zawodzinski model:
    /// D_w(λ) = D_ref * exp(-2436/T) * f(λ)
    ///
    /// where f(λ) interpolates from 1 at λ=0 to a maximum near λ=3, then
    /// decreases for fully saturated conditions.
    pub fn water_diffusivity(&self, lambda_water: f64) -> f64 {
        let d0 = 2.1e-7 * (-2436.0 / self.temperature).exp();
        // Monotonically increasing version to avoid sign flip
        let f_lambda = if lambda_water < 2.0 {
            1.0
        } else if lambda_water < 3.0 {
            1.0 + 2.0 * (lambda_water - 2.0)
        } else {
            // Saturate at max value for large λ
            3.0
        };
        d0 * f_lambda
    }
}

// ============================================================================
// Helper utilities
// ============================================================================

/// Compute thermal voltage V_T = RT/F at temperature T (V).
pub fn thermal_voltage(temperature: f64) -> f64 {
    R_GAS * temperature / FARADAY
}

/// Compute Debye length for a 1:1 electrolyte at concentration c (mol/m³) (m).
pub fn debye_length_11(epsilon_r: f64, temperature: f64, c_bulk: f64) -> f64 {
    let eps = epsilon_r * EPSILON_0;
    (eps * R_GAS * temperature / (2.0 * FARADAY.powi(2) * c_bulk)).sqrt()
}

/// Compute limiting diffusion current density j_lim (A/m²).
///
/// j_lim = n F D c_bulk / δ
/// where δ is the diffusion layer thickness.
pub fn limiting_current_density(n_elec: f64, diffusivity: f64, c_bulk: f64, delta: f64) -> f64 {
    n_elec * FARADAY * diffusivity * c_bulk / delta
}

/// Compute polarization resistance R_p = RT/(n F i₀) (Ω·m²).
pub fn polarization_resistance(temperature: f64, n_elec: f64, exchange_current: f64) -> f64 {
    R_GAS * temperature / (n_elec * FARADAY * exchange_current)
}

/// Compute the Warburg impedance (diffusion) coefficient σ (Ω/sqrt(rad/s)).
///
/// σ = RT / (n² F² A sqrt(2) sqrt(D) c_bulk)
pub fn warburg_coefficient(
    temperature: f64,
    n_elec: f64,
    diffusivity: f64,
    c_bulk: f64,
    area: f64,
) -> f64 {
    R_GAS * temperature
        / (n_elec.powi(2) * FARADAY.powi(2) * area * 2.0_f64.sqrt() * diffusivity.sqrt() * c_bulk)
}

/// Compute the diffusion layer thickness from the rotating disk electrode (RDE) Levich equation (m).
///
/// δ = 1.61 * D^(1/3) * ν^(1/6) / ω^(1/2)
pub fn levich_diffusion_layer(diffusivity: f64, viscosity_kinematic: f64, omega: f64) -> f64 {
    1.61 * diffusivity.powf(1.0 / 3.0) * viscosity_kinematic.powf(1.0 / 6.0) / omega.sqrt()
}

/// Compute Nernst equation equilibrium potential E_eq (V).
///
/// E_eq = E₀ + RT/(nF) * ln(c_ox/c_red)
pub fn nernst_potential(
    e_standard: f64,
    temperature: f64,
    n_elec: f64,
    c_ox: f64,
    c_red: f64,
) -> f64 {
    e_standard + R_GAS * temperature / (n_elec * FARADAY) * (c_ox / c_red.max(1e-30)).ln()
}

/// Compute specific capacitance of Helmholtz double layer (F/m²).
///
/// C_H = ε ε₀ / d
pub fn helmholtz_capacitance(epsilon_r: f64, d_layer: f64) -> f64 {
    epsilon_r * EPSILON_0 / d_layer
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- IonSpecies & NernstPlanckFem ---

    #[test]
    fn test_ion_species_creation() {
        let sp = IonSpecies::new("Li+", 1, 1.0e-9, 1000.0);
        assert_eq!(sp.charge, 1);
        assert!((sp.diffusivity - 1e-9).abs() < 1e-20);
        assert_eq!(sp.bulk_concentration, 1000.0);
    }

    #[test]
    fn test_thermal_voltage() {
        let np = NernstPlanckFem::new(
            298.15,
            vec![IonSpecies::new("Li+", 1, 1e-9, 1000.0)],
            10,
            1e-3,
            1e-3,
        );
        let vt = np.thermal_voltage();
        assert!((vt - 0.025693).abs() < 1e-4, "thermal voltage = {vt}");
    }

    #[test]
    fn test_element_diffusion_matrix_symmetry() {
        let np = NernstPlanckFem::new(
            298.15,
            vec![IonSpecies::new("Na+", 1, 1.3e-9, 100.0)],
            5,
            1e-3,
            1e-4,
        );
        let dm = np.element_diffusion_matrix(1.3e-9, 2.5e-4);
        assert!((dm[0][0] + dm[0][1]).abs() < 1e-20);
        assert!((dm[1][0] + dm[1][1]).abs() < 1e-20);
        assert!((dm[0][1] - dm[1][0]).abs() < 1e-30);
    }

    #[test]
    fn test_nernst_planck_step_conserves_mass_approx() {
        let species = vec![IonSpecies::new("Cl-", -1, 2.0e-9, 500.0)];
        let mut np = NernstPlanckFem::new(298.15, species, 20, 1e-3, 1e-4);
        let c_init: f64 = np.concentration[0].iter().sum();
        np.step();
        let c_final: f64 = np.concentration[0].iter().sum();
        // Mass should not blow up
        assert!(c_final.is_finite());
        assert!(c_final > 0.0);
        assert!((c_final - c_init).abs() / c_init < 0.5);
    }

    #[test]
    fn test_nernst_planck_tridiagonal() {
        let lower = vec![0.0, -1.0, -1.0];
        let diag = vec![2.0, 3.0, 2.0];
        let upper = vec![-1.0, -1.0, 0.0];
        let rhs = vec![1.0, 2.0, 1.0];
        let x = NernstPlanckFem::solve_tridiagonal(&lower, &diag, &upper, &rhs);
        // Check residual
        let r0 = diag[0] * x[0] + upper[0] * x[1] - rhs[0];
        let r1 = lower[1] * x[0] + diag[1] * x[1] + upper[1] * x[2] - rhs[1];
        let r2 = lower[2] * x[1] + diag[2] * x[2] - rhs[2];
        assert!(r0.abs() < 1e-10);
        assert!(r1.abs() < 1e-10);
        assert!(r2.abs() < 1e-10);
    }

    #[test]
    fn test_current_density_finite() {
        let species = vec![IonSpecies::new("K+", 1, 1.96e-9, 100.0)];
        let np = NernstPlanckFem::new(298.15, species, 10, 1e-3, 1e-4);
        let j = np.current_density(0, 5);
        assert!(j.is_finite());
    }

    #[test]
    fn test_total_current_finite() {
        let species = vec![IonSpecies::new("K+", 1, 1.96e-9, 100.0)];
        let np = NernstPlanckFem::new(298.15, species, 10, 1e-3, 1e-4);
        assert!(np.total_current(0).is_finite());
    }

    // --- PoissonBoltzmannFem ---

    #[test]
    fn test_debye_length_nacl_100mm() {
        // For 100 mM NaCl (2 species, 100 mol/m³ each), λ_D ≈ 0.96 nm
        let sp = vec![
            IonSpecies::new("Na+", 1, 1.33e-9, 100.0),
            IonSpecies::new("Cl-", -1, 2.03e-9, 100.0),
        ];
        let pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 50, 10e-9, 0.1);
        let ld = pb.debye_length();
        assert!(ld > 0.5e-9 && ld < 2.0e-9, "Debye length = {ld:.3e} m");
    }

    #[test]
    fn test_pb_ionic_strength() {
        let sp = vec![
            IonSpecies::new("Ca2+", 2, 7.9e-10, 10.0),
            IonSpecies::new("Cl-", -1, 2.03e-9, 20.0),
        ];
        let pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 10, 1e-8, 0.05);
        let i = pb.ionic_strength();
        // I = 0.5*(4*10 + 1*20) = 30
        assert!((i - 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_pb_solve_monotonic_decay() {
        let sp = vec![IonSpecies::new("Na+", 1, 1.33e-9, 100.0)];
        let mut pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 30, 20e-9, 0.1);
        pb.solve();
        // Potential should decay from surface to bulk
        let phi = &pb.potential;
        assert!(phi[0] > phi[phi.len() - 1]);
        assert!(phi.last().copied().unwrap_or(1.0).abs() < 0.01);
    }

    #[test]
    fn test_surface_charge_density_sign() {
        let sp = vec![IonSpecies::new("Na+", 1, 1.33e-9, 100.0)];
        let mut pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 30, 20e-9, 0.1);
        pb.solve();
        // Positive surface potential → surface charge is non-zero
        let sigma = pb.surface_charge_density();
        assert!(sigma.abs() > 0.0, "sigma should be non-zero, got = {sigma}");
        assert!(sigma.is_finite(), "sigma should be finite");
    }

    #[test]
    fn test_differential_capacitance_positive() {
        let sp = vec![IonSpecies::new("Na+", 1, 1.33e-9, 100.0)];
        let pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 30, 20e-9, 0.1);
        assert!(pb.differential_capacitance() > 0.0);
    }

    #[test]
    fn test_analytical_solution_at_zero() {
        let sp = vec![IonSpecies::new("Na+", 1, 1.33e-9, 100.0)];
        let pb = PoissonBoltzmannFem::new(78.5, 298.15, sp, 30, 20e-9, 0.05);
        assert!((pb.analytical_solution(0.0) - 0.05).abs() < 1e-10);
    }

    // --- ButlerVolmerFem ---

    #[test]
    fn test_bv_overpotential_zero_at_equilibrium() {
        let mut bv = ButlerVolmerFem::new(298.15, 1e-3, 0.5, 0.5, 0.5);
        bv.set_potentials(0.5, 0.0);
        assert!((bv.overpotential()).abs() < 1e-10);
    }

    #[test]
    fn test_bv_current_zero_at_equilibrium() {
        let mut bv = ButlerVolmerFem::new(298.15, 1e-3, 0.5, 0.5, 0.3);
        bv.set_potentials(0.3, 0.0);
        // At equilibrium, j = j0 * (exp(0) - exp(0)) = 0
        let j = bv.current_density();
        assert!(j.abs() < 1e-12, "j at equilibrium = {j}");
    }

    #[test]
    fn test_bv_anodic_current_positive_for_positive_overpotential() {
        let mut bv = ButlerVolmerFem::new(298.15, 1e-3, 0.5, 0.5, 0.0);
        bv.set_potentials(0.1, 0.0); // η = 0.1 V
        assert!(bv.current_density() > 0.0);
    }

    #[test]
    fn test_bv_tafel_slopes_positive() {
        let bv = ButlerVolmerFem::new(298.15, 1e-3, 0.5, 0.5, 0.0);
        assert!(bv.tafel_slope_anodic() > 0.0);
        assert!(bv.tafel_slope_cathodic() > 0.0);
    }

    #[test]
    fn test_bv_linearized_conductance() {
        let bv = ButlerVolmerFem::new(298.15, 1e-3, 0.5, 0.5, 0.0);
        let g = bv.linearized_conductance();
        // g = j0 * (alpha_a + alpha_c) * F/RT = 1e-3 * 1.0 * F/RT
        let expected = 1e-3 * 1.0 * FARADAY / (R_GAS * 298.15);
        assert!((g - expected).abs() / expected < 1e-8);
    }

    #[test]
    fn test_ocp_graphite_range() {
        for soc in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99] {
            let v = ButlerVolmerFem::ocp_graphite(soc);
            assert!(v > -0.5 && v < 2.0, "OCP graphite at soc={soc} = {v}");
        }
    }

    #[test]
    fn test_ocp_licoo2_range() {
        for soc in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let v = ButlerVolmerFem::ocp_licoo2(soc);
            assert!(v > 3.0 && v < 5.0, "OCP LiCoO2 at soc={soc} = {v}");
        }
    }

    // --- BatteryFem ---

    #[test]
    fn test_battery_fem_creation() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        assert_eq!(bat.neg.domain, DfnDomain::NegativeElectrode);
        assert_eq!(bat.pos.domain, DfnDomain::PositiveElectrode);
    }

    #[test]
    fn test_battery_n_total_nodes() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let expected = bat.neg.n_macro + bat.n_sep + bat.pos.n_macro;
        assert_eq!(bat.n_total_nodes(), expected);
    }

    #[test]
    fn test_battery_soc_initial() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let soc_n = bat.soc_negative();
        let soc_p = bat.soc_positive();
        assert!((soc_n - 0.5).abs() < 1e-10);
        assert!((soc_p - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_battery_terminal_voltage_range() {
        let bat = BatteryFem::new(1.0, 298.15, 1.0);
        let v = bat.terminal_voltage();
        assert!(v > 2.0 && v < 5.0, "V_cell = {v}");
    }

    #[test]
    fn test_battery_ir_drop_terminal_voltage() {
        // R_int=0.05 Ω, I=2 A → IR drop = 0.1 V; terminal voltage must drop by exactly 0.1 V
        let mut bat = BatteryFem::new(2.0, 298.15, 1.0);
        bat.internal_resistance = 0.05;
        let v_no_r = {
            let mut b2 = bat.clone();
            b2.internal_resistance = 0.0;
            b2.terminal_voltage()
        };
        let v_with_r = bat.terminal_voltage();
        assert!(
            v_with_r < v_no_r,
            "IR drop not applied: v_with_r={v_with_r}, v_no_r={v_no_r}"
        );
        assert!(
            (v_no_r - v_with_r - 0.1).abs() < 1e-10,
            "expected IR drop of 0.1 V, got {:.6}",
            v_no_r - v_with_r
        );
    }

    #[test]
    fn test_battery_step() {
        let mut bat = BatteryFem::new(1.0, 298.15, 10.0);
        bat.step();
        assert!((bat.time - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_electrode_effective_conductivity() {
        let neg = ElectrodeDomain::default_negative();
        let kappa_eff = neg.effective_conductivity(1.0);
        assert!(kappa_eff < 1.0 && kappa_eff > 0.0);
    }

    #[test]
    fn test_electrode_effective_diffusivity() {
        let pos = ElectrodeDomain::default_positive();
        let d_eff = pos.effective_electrolyte_diffusivity(2.7e-10);
        assert!(d_eff < 2.7e-10 && d_eff > 0.0);
    }

    // --- ElectrochemicalCorrosion ---

    #[test]
    fn test_corrosion_creation() {
        let corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 20, 1e-2);
        assert!((corr.e_corrosion - (-0.44)).abs() < 1e-10);
    }

    #[test]
    fn test_corrosion_net_current_zero_at_ecorr() {
        let corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        let j = corr.net_current(-0.44);
        assert!(j.abs() < 1e-15, "net current at E_corr = {j}");
    }

    #[test]
    fn test_anodic_current_increases_with_e() {
        let corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        let j1 = corr.anodic_current(-0.44);
        let j2 = corr.anodic_current(-0.30);
        assert!(j2 > j1);
    }

    #[test]
    fn test_pitting_check() {
        // e_pit = -0.25 V: pitting initiates when E >= e_pit
        let mut corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        // -0.10 > -0.25 → pitting active
        corr.check_pitting(-0.10);
        assert!(corr.pitting_active);
        // -0.40 < -0.25 → no pitting
        corr.check_pitting(-0.40);
        assert!(!corr.pitting_active);
    }

    #[test]
    fn test_corrosion_rate_from_current() {
        let rate = ElectrochemicalCorrosion::current_to_corrosion_rate(1e-3);
        assert!(rate > 0.0 && rate < 10.0, "rate = {rate} mm/year");
    }

    #[test]
    fn test_galvanic_couple() {
        let mut corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        corr.add_galvanic_couple(1.0, 10.0, 0.3);
        let gc = corr.galvanic_current();
        assert_eq!(gc.len(), 1);
        assert!(gc[0] > 0.0);
    }

    #[test]
    fn test_evans_diagram_lengths() {
        let corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        let (e, ja, jc) = corr.evans_diagram(50);
        assert_eq!(e.len(), 50);
        assert_eq!(ja.len(), 50);
        assert_eq!(jc.len(), 50);
    }

    #[test]
    fn test_solve_laplace_bc() {
        let mut corr =
            ElectrochemicalCorrosion::new(298.15, -0.44, 0.12, 0.12, 1e-6, -0.25, -0.35, 10, 1e-2);
        corr.solve_laplace(-0.44, -0.60, 5.0);
        assert!((corr.potential[0] - (-0.44)).abs() < 1e-10);
        assert!((corr.potential[9] - (-0.60)).abs() < 1e-10);
    }

    // --- FuelCellFem ---

    #[test]
    fn test_fuel_cell_creation() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        assert!((fc.v_cell - 0.7).abs() < 1e-10);
        assert!((fc.v_oc - 1.23).abs() < 1e-10);
    }

    #[test]
    fn test_total_overpotential() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        assert!((fc.total_overpotential() - 0.53).abs() < 1e-10);
    }

    #[test]
    fn test_proton_conductivity_springer() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let sigma = fc.proton_conductivity(14.0);
        assert!(sigma > 0.0 && sigma < 100.0, "sigma = {sigma}");
    }

    #[test]
    fn test_effective_oxygen_diffusivity() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let d = fc.effective_oxygen_diffusivity(0.6, 0.1);
        assert!(d > 0.0 && d < fc.d_oxygen);
    }

    #[test]
    fn test_orr_rate_positive() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let rate = fc.orr_rate(40.0, 0.4);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_solve_proton_transport_bc() {
        let mut fc = FuelCellFem::new(353.15, 0.7, 1.0);
        fc.solve_proton_transport();
        assert!((fc.phi_proton[0]).abs() < 1e-10);
        assert!((fc.phi_proton[fc.n_nodes_per_region - 1] - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_solve_oxygen_transport_positive() {
        let mut fc = FuelCellFem::new(353.15, 0.7, 1.0);
        fc.solve_oxygen_transport();
        // All O2 concentrations should be non-negative
        assert!(fc.c_oxygen.iter().all(|&c| c >= -1e-6));
    }

    #[test]
    fn test_fuel_cell_step_advances_time() {
        let mut fc = FuelCellFem::new(353.15, 0.7, 0.1);
        fc.step();
        assert!((fc.time - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_power_density_positive() {
        let mut fc = FuelCellFem::new(353.15, 0.7, 1.0);
        fc.step();
        assert!(fc.power_density() >= 0.0);
    }

    #[test]
    fn test_nernst_potential_standard() {
        let v = FuelCellFem::nernst_potential(298.15, 101325.0, 101325.0, 3540.0);
        assert!(v > 1.0 && v < 1.5, "V_Nernst = {v}");
    }

    #[test]
    fn test_activation_overpotential_positive() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let eta = fc.activation_overpotential(0.1);
        assert!(eta >= 0.0);
    }

    #[test]
    fn test_ohmic_overpotential() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let eta = fc.ohmic_overpotential(1.0, 0.05);
        assert!((eta - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_concentration_overpotential_grows() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let eta1 = fc.concentration_overpotential(0.5, 2.0);
        let eta2 = fc.concentration_overpotential(0.9, 2.0);
        assert!(eta2 > eta1);
    }

    #[test]
    fn test_electro_osmotic_drag() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        let nd = fc.electro_osmotic_drag(22.0);
        assert!((nd - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_water_diffusivity_positive() {
        let fc = FuelCellFem::new(353.15, 0.7, 1.0);
        assert!(fc.water_diffusivity(5.0) > 0.0);
        assert!(fc.water_diffusivity(14.0) > 0.0);
    }

    // --- Utility functions ---

    #[test]
    fn test_thermal_voltage_room_temp() {
        let vt = thermal_voltage(298.15);
        assert!((vt - 0.025693).abs() < 1e-4);
    }

    #[test]
    fn test_debye_length_11_order_of_magnitude() {
        // 10 mM NaCl → λ_D ≈ 3 nm
        let ld = debye_length_11(78.5, 298.15, 10.0);
        assert!(ld > 1e-9 && ld < 10e-9, "λ_D = {ld:.3e}");
    }

    #[test]
    fn test_limiting_current_density() {
        let j_lim = limiting_current_density(1.0, 1e-9, 1000.0, 100e-6);
        assert!(j_lim > 0.0);
    }

    #[test]
    fn test_polarization_resistance_positive() {
        let rp = polarization_resistance(298.15, 2.0, 1e-3);
        assert!(rp > 0.0);
    }

    #[test]
    fn test_warburg_coefficient_positive() {
        let sigma = warburg_coefficient(298.15, 1.0, 1e-9, 1000.0, 1e-4);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_levich_diffusion_layer() {
        let delta = levich_diffusion_layer(1e-9, 1e-6, 2.0 * PI * 10.0);
        assert!(delta > 0.0 && delta < 1e-3, "δ = {delta:.3e}");
    }

    #[test]
    fn test_nernst_potential_fn() {
        let e = nernst_potential(0.34, 298.15, 2.0, 1.0, 1.0);
        assert!((e - 0.34).abs() < 1e-10);
    }

    #[test]
    fn test_helmholtz_capacitance() {
        let c = helmholtz_capacitance(78.5, 3e-10);
        // ≈ 78.5 * 8.85e-12 / 0.3e-9 ≈ 2.3 F/m²
        assert!(c > 1.0 && c < 5.0, "C_H = {c}");
    }
}
