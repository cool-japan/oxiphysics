// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Foam and cellular material models.
//!
//! Provides Gibson-Ashby scaling laws for open/closed-cell foams,
//! Kelvin cell geometry, energy absorption (plateau stress, densification),
//! viscoelastic foam (Prony series), crushable foam plasticity,
//! foaming kinetics, thermal conductivity models, and impact attenuation.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ===========================================================================
// Cell type enum
// ===========================================================================

/// Type of foam cell structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellType {
    /// Open-cell foam (interconnected pores).
    Open,
    /// Closed-cell foam (sealed pores with gas).
    Closed,
    /// Mixed cell structure with a given fraction closed (0..1).
    Mixed,
}

// ===========================================================================
// Gibson-Ashby scaling laws
// ===========================================================================

/// Gibson-Ashby scaling law parameters for cellular solids.
///
/// Predicts effective modulus, strength, and thermal conductivity from
/// relative density using power-law scaling.
#[derive(Debug, Clone, PartialEq)]
pub struct GibsonAshby {
    /// Density of the solid cell-wall material \[kg/m³\].
    pub rho_solid: f64,
    /// Young's modulus of the solid cell-wall material \[Pa\].
    pub e_solid: f64,
    /// Compressive strength of the solid cell-wall material \[Pa\].
    pub sigma_solid: f64,
    /// Thermal conductivity of the solid material \[W/(m·K)\].
    pub k_solid: f64,
    /// Cell type.
    pub cell_type: CellType,
    /// Fraction of closed cells (only used when `cell_type == Mixed`).
    pub closed_fraction: f64,
}

impl GibsonAshby {
    /// Create a new Gibson-Ashby model.
    pub fn new(
        rho_solid: f64,
        e_solid: f64,
        sigma_solid: f64,
        k_solid: f64,
        cell_type: CellType,
        closed_fraction: f64,
    ) -> Self {
        Self {
            rho_solid,
            e_solid,
            sigma_solid,
            k_solid,
            cell_type,
            closed_fraction,
        }
    }

    /// Relative density (foam density / solid density).
    pub fn relative_density(&self, rho_foam: f64) -> f64 {
        rho_foam / self.rho_solid
    }

    /// Effective Young's modulus for open-cell foam.
    ///
    /// E* / E_s = C1 * (rho*/rho_s)^2   (C1 ≈ 1)
    pub fn modulus_open(&self, rho_foam: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        self.e_solid * rd * rd
    }

    /// Effective Young's modulus for closed-cell foam.
    ///
    /// E* / E_s ≈ phi^2 * (rho*/rho_s)^2 + (1-phi) * (rho*/rho_s)
    /// where phi is fraction of solid in cell edges (~0.6 typical).
    pub fn modulus_closed(&self, rho_foam: f64, phi_edge: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        self.e_solid * (phi_edge * phi_edge * rd * rd + (1.0 - phi_edge) * rd)
    }

    /// Effective modulus based on cell type.
    pub fn effective_modulus(&self, rho_foam: f64) -> f64 {
        match self.cell_type {
            CellType::Open => self.modulus_open(rho_foam),
            CellType::Closed => self.modulus_closed(rho_foam, 0.6),
            CellType::Mixed => {
                let e_open = self.modulus_open(rho_foam);
                let e_closed = self.modulus_closed(rho_foam, 0.6);
                let f = self.closed_fraction;
                f * e_closed + (1.0 - f) * e_open
            }
        }
    }

    /// Compressive (collapse) strength for open-cell foam.
    ///
    /// sigma* / sigma_s = C2 * (rho*/rho_s)^(3/2)   (C2 ≈ 0.3)
    pub fn strength_open(&self, rho_foam: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        0.3 * self.sigma_solid * rd.powf(1.5)
    }

    /// Compressive strength for closed-cell foam.
    ///
    /// sigma* / sigma_s = C3 * (phi * rho*/rho_s)^(3/2) + (1-phi) * (rho*/rho_s)
    pub fn strength_closed(&self, rho_foam: f64, phi_edge: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        self.sigma_solid * (0.3 * (phi_edge * rd).powf(1.5) + (1.0 - phi_edge) * rd)
    }

    /// Effective compressive strength based on cell type.
    pub fn effective_strength(&self, rho_foam: f64) -> f64 {
        match self.cell_type {
            CellType::Open => self.strength_open(rho_foam),
            CellType::Closed => self.strength_closed(rho_foam, 0.6),
            CellType::Mixed => {
                let s_open = self.strength_open(rho_foam);
                let s_closed = self.strength_closed(rho_foam, 0.6);
                let f = self.closed_fraction;
                f * s_closed + (1.0 - f) * s_open
            }
        }
    }

    /// Gibson-Ashby thermal conductivity scaling.
    ///
    /// k* / k_s ≈ (1/3) * (rho*/rho_s)  for open-cell
    pub fn thermal_conductivity_open(&self, rho_foam: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        self.k_solid * rd / 3.0
    }

    /// Gibson-Ashby Poisson's ratio estimate for open-cell foam.
    ///
    /// For open-cell foams, Poisson's ratio is approximately 1/3.
    pub fn poisson_ratio_open(&self) -> f64 {
        1.0 / 3.0
    }

    /// Fracture toughness scaling for open-cell foam.
    ///
    /// K_Ic* / (sigma_s * sqrt(l)) ≈ C * (rho*/rho_s)^(3/2)
    /// where l is cell size and C ≈ 0.65.
    pub fn fracture_toughness_open(&self, rho_foam: f64, cell_size: f64) -> f64 {
        let rd = self.relative_density(rho_foam);
        0.65 * self.sigma_solid * cell_size.sqrt() * rd.powf(1.5)
    }
}

// ===========================================================================
// Kelvin cell model
// ===========================================================================

/// Kelvin cell (truncated octahedron) foam geometry model.
///
/// The Kelvin cell is a space-filling polyhedron (truncated octahedron) that
/// provides a good approximation to real foam cells.
#[derive(Debug, Clone, PartialEq)]
pub struct KelvinCell {
    /// Cell edge length \[m\].
    pub edge_length: f64,
    /// Cell wall thickness \[m\].
    pub wall_thickness: f64,
    /// Solid material density \[kg/m³\].
    pub rho_solid: f64,
    /// Solid material modulus \[Pa\].
    pub e_solid: f64,
}

impl KelvinCell {
    /// Create a new Kelvin cell model.
    pub fn new(edge_length: f64, wall_thickness: f64, rho_solid: f64, e_solid: f64) -> Self {
        Self {
            edge_length,
            wall_thickness,
            rho_solid,
            e_solid,
        }
    }

    /// Volume of a single Kelvin cell (truncated octahedron).
    ///
    /// V = 8 * sqrt(2) * a^3 where a is the edge length.
    pub fn cell_volume(&self) -> f64 {
        8.0 * 2.0_f64.sqrt() * self.edge_length.powi(3)
    }

    /// Surface area of a single Kelvin cell.
    ///
    /// A = (6 + 12*sqrt(3)) * a^2
    pub fn cell_surface_area(&self) -> f64 {
        (6.0 + 12.0 * 3.0_f64.sqrt()) * self.edge_length.powi(2)
    }

    /// Number of faces: 8 hexagonal + 6 square = 14.
    pub fn num_faces(&self) -> usize {
        14
    }

    /// Number of edges of a Kelvin cell = 36.
    pub fn num_edges(&self) -> usize {
        36
    }

    /// Number of vertices of a Kelvin cell = 24.
    pub fn num_vertices(&self) -> usize {
        24
    }

    /// Relative density for closed-cell Kelvin foam.
    ///
    /// Approximation: rho*/rho_s ≈ surface_area * t / cell_volume
    pub fn relative_density(&self) -> f64 {
        let a = self.cell_surface_area();
        let v = self.cell_volume();
        (a * self.wall_thickness / v).min(1.0)
    }

    /// Effective foam density \[kg/m³\].
    pub fn foam_density(&self) -> f64 {
        self.rho_solid * self.relative_density()
    }

    /// Effective modulus via Gibson-Ashby with Kelvin cell relative density.
    pub fn effective_modulus(&self) -> f64 {
        let rd = self.relative_density();
        self.e_solid * rd * rd
    }

    /// Characteristic cell diameter (equivalent sphere diameter).
    pub fn equivalent_diameter(&self) -> f64 {
        let v = self.cell_volume();
        (6.0 * v / PI).powf(1.0 / 3.0)
    }

    /// Strut cross-section area for open-cell Kelvin model.
    ///
    /// For open-cell, material is concentrated in the 36 edges.
    /// A_strut ≈ (rho*/rho_s) * V_cell / (N_edges * L_edge)
    pub fn strut_area(&self, relative_density: f64) -> f64 {
        let v = self.cell_volume();
        relative_density * v / (self.num_edges() as f64 * self.edge_length)
    }
}

// ===========================================================================
// Energy absorption model
// ===========================================================================

/// Energy absorption model for foam impact.
///
/// Describes the three characteristic regions of foam compression:
/// linear elastic, plateau (constant stress), and densification.
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyAbsorption {
    /// Plateau stress \[Pa\].
    pub plateau_stress: f64,
    /// Densification strain (typically 0.5–0.9).
    pub densification_strain: f64,
    /// Young's modulus in the linear elastic region \[Pa\].
    pub elastic_modulus: f64,
    /// Yield strain (transition from elastic to plateau).
    pub yield_strain: f64,
    /// Foam density \[kg/m³\].
    pub foam_density: f64,
    /// Foam thickness \[m\].
    pub thickness: f64,
}

impl EnergyAbsorption {
    /// Create a new energy absorption model.
    pub fn new(
        plateau_stress: f64,
        densification_strain: f64,
        elastic_modulus: f64,
        yield_strain: f64,
        foam_density: f64,
        thickness: f64,
    ) -> Self {
        Self {
            plateau_stress,
            densification_strain,
            elastic_modulus,
            yield_strain,
            foam_density,
            thickness,
        }
    }

    /// Stress at a given engineering strain (simplified tri-linear model).
    pub fn stress_at_strain(&self, strain: f64) -> f64 {
        if strain < 0.0 {
            return 0.0;
        }
        if strain <= self.yield_strain {
            // Linear elastic region
            self.elastic_modulus * strain
        } else if strain <= self.densification_strain {
            // Plateau region
            self.plateau_stress
        } else {
            // Densification: exponential rise
            let excess = strain - self.densification_strain;
            let stiffening_factor = 10.0; // empirical
            self.plateau_stress * (1.0 + stiffening_factor * excess * excess)
        }
    }

    /// Energy absorbed per unit volume up to a given strain \[J/m³\].
    ///
    /// Integrates the stress-strain curve from 0 to `strain`.
    pub fn energy_density(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            return 0.0;
        }
        if strain <= self.yield_strain {
            // Triangular area under linear elastic
            0.5 * self.elastic_modulus * strain * strain
        } else if strain <= self.densification_strain {
            // Elastic area + plateau rectangle
            let elastic_energy = 0.5 * self.elastic_modulus * self.yield_strain * self.yield_strain;
            let plateau_energy = self.plateau_stress * (strain - self.yield_strain);
            elastic_energy + plateau_energy
        } else {
            // Elastic + plateau + densification
            let elastic_energy = 0.5 * self.elastic_modulus * self.yield_strain * self.yield_strain;
            let plateau_energy =
                self.plateau_stress * (self.densification_strain - self.yield_strain);
            let excess = strain - self.densification_strain;
            let k = 10.0;
            // Integral of sigma_p * (1 + k * x^2) dx from 0 to excess
            let densification_energy = self.plateau_stress * (excess + k * excess.powi(3) / 3.0);
            elastic_energy + plateau_energy + densification_energy
        }
    }

    /// Total energy absorbed through full plateau region \[J/m³\].
    pub fn plateau_energy_density(&self) -> f64 {
        self.energy_density(self.densification_strain)
    }

    /// Efficiency at a given strain.
    ///
    /// eta = W(strain) / (sigma(strain) * strain)
    pub fn efficiency(&self, strain: f64) -> f64 {
        if strain <= 1.0e-12 {
            return 0.0;
        }
        let w = self.energy_density(strain);
        let sigma = self.stress_at_strain(strain);
        if sigma.abs() < 1.0e-30 {
            return 0.0;
        }
        w / (sigma * strain)
    }

    /// Ideal energy absorption efficiency (rectangular stress-strain).
    ///
    /// eta_ideal = 1 for perfect energy absorber.
    pub fn ideal_efficiency(&self) -> f64 {
        // For the plateau region only
        let w = self.plateau_energy_density();
        let sigma_max = self.stress_at_strain(self.densification_strain);
        if sigma_max.abs() < 1.0e-30 {
            return 0.0;
        }
        w / (sigma_max * self.densification_strain)
    }

    /// Energy absorbed per unit mass up to a given strain \[J/kg\].
    pub fn specific_energy(&self, strain: f64) -> f64 {
        if self.foam_density.abs() < 1.0e-30 {
            return 0.0;
        }
        self.energy_density(strain) / self.foam_density
    }

    /// Maximum deceleration during impact for given drop parameters.
    ///
    /// G = sigma_plateau / (rho_foam * h)
    /// where h is the drop height.
    pub fn peak_g_force(&self, drop_height: f64) -> f64 {
        if self.foam_density.abs() < 1.0e-30 || self.thickness.abs() < 1.0e-30 {
            return 0.0;
        }
        self.plateau_stress / (self.foam_density * drop_height)
    }

    /// Cushion factor = peak_stress / energy_density at that strain.
    pub fn cushion_factor(&self, strain: f64) -> f64 {
        let w = self.energy_density(strain);
        if w.abs() < 1.0e-30 {
            return 0.0;
        }
        self.stress_at_strain(strain) / w
    }
}

// ===========================================================================
// Viscoelastic foam (Prony series)
// ===========================================================================

/// Single Prony series term for viscoelastic relaxation.
#[derive(Debug, Clone, PartialEq)]
pub struct PronyTerm {
    /// Relative modulus weight g_i (dimensionless).
    pub g_i: f64,
    /// Relaxation time tau_i \[s\].
    pub tau_i: f64,
}

/// Viscoelastic foam model using a Prony series representation.
///
/// The relaxation modulus is:
///   E(t) = E_inf + sum_i (E_0 * g_i * exp(-t / tau_i))
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoelasticFoam {
    /// Instantaneous modulus E_0 \[Pa\].
    pub e_instantaneous: f64,
    /// Long-term (equilibrium) modulus E_inf \[Pa\].
    pub e_equilibrium: f64,
    /// Prony series terms.
    pub prony_terms: Vec<PronyTerm>,
    /// Foam density \[kg/m³\].
    pub density: f64,
    /// Poisson's ratio (typically ~0.0 for low-density foams).
    pub poisson_ratio: f64,
}

impl ViscoelasticFoam {
    /// Create a new viscoelastic foam model.
    pub fn new(
        e_instantaneous: f64,
        e_equilibrium: f64,
        prony_terms: Vec<PronyTerm>,
        density: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            e_instantaneous,
            e_equilibrium,
            prony_terms,
            density,
            poisson_ratio,
        }
    }

    /// Create a single-term Maxwell model.
    pub fn single_term(e_inst: f64, e_eq: f64, tau: f64, density: f64) -> Self {
        let g = (e_inst - e_eq) / e_inst;
        Self::new(
            e_inst,
            e_eq,
            vec![PronyTerm { g_i: g, tau_i: tau }],
            density,
            0.0,
        )
    }

    /// Relaxation modulus at time t.
    ///
    /// E(t) = E_inf + E_0 * sum_i g_i * exp(-t / tau_i)
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let mut e = self.e_equilibrium;
        for term in &self.prony_terms {
            e += self.e_instantaneous * term.g_i * (-t / term.tau_i).exp();
        }
        e
    }

    /// Stress response to a constant-rate strain loading.
    ///
    /// For strain(t) = rate * t:
    ///   sigma(t) = rate * integral_0^t E(t - s) ds
    pub fn stress_constant_rate(&self, strain_rate: f64, t: f64) -> f64 {
        // Analytical integration for Prony series under constant rate
        let mut sigma = self.e_equilibrium * strain_rate * t;
        for term in &self.prony_terms {
            let tau = term.tau_i;
            let g = term.g_i;
            // integral_0^t g_i * E_0 * exp(-(t-s)/tau) * rate ds
            // = rate * E_0 * g_i * tau * (1 - exp(-t/tau))
            sigma += strain_rate * self.e_instantaneous * g * tau * (1.0 - (-t / tau).exp());
        }
        sigma
    }

    /// Creep compliance at time t (reciprocal relationship, approximate).
    ///
    /// J(t) ≈ 1 / E(t)   (approximate for linear viscoelastic).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let e = self.relaxation_modulus(t);
        if e.abs() < 1.0e-30 {
            return f64::INFINITY;
        }
        1.0 / e
    }

    /// Loss tangent (tan delta) at angular frequency omega.
    ///
    /// tan(delta) = E''(omega) / E'(omega)
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let (e_storage, e_loss) = self.dynamic_moduli(omega);
        if e_storage.abs() < 1.0e-30 {
            return 0.0;
        }
        e_loss / e_storage
    }

    /// Storage and loss moduli at angular frequency omega.
    ///
    /// E'(omega) = E_inf + E_0 * sum g_i * (omega*tau_i)^2 / (1 + (omega*tau_i)^2)
    /// E''(omega) = E_0 * sum g_i * omega*tau_i / (1 + (omega*tau_i)^2)
    pub fn dynamic_moduli(&self, omega: f64) -> (f64, f64) {
        let mut e_storage = self.e_equilibrium;
        let mut e_loss = 0.0;
        for term in &self.prony_terms {
            let wt = omega * term.tau_i;
            let wt2 = wt * wt;
            let denom = 1.0 + wt2;
            e_storage += self.e_instantaneous * term.g_i * wt2 / denom;
            e_loss += self.e_instantaneous * term.g_i * wt / denom;
        }
        (e_storage, e_loss)
    }

    /// Half-time for stress relaxation.
    ///
    /// Approximate: uses the longest relaxation time.
    pub fn relaxation_half_time(&self) -> f64 {
        let max_tau = self
            .prony_terms
            .iter()
            .map(|t| t.tau_i)
            .fold(0.0_f64, f64::max);
        max_tau * 2.0_f64.ln()
    }

    /// Bulk modulus from E and nu.
    pub fn bulk_modulus(&self, t: f64) -> f64 {
        let e = self.relaxation_modulus(t);
        let nu = self.poisson_ratio;
        e / (3.0 * (1.0 - 2.0 * nu))
    }

    /// Shear modulus from E and nu.
    pub fn shear_modulus(&self, t: f64) -> f64 {
        let e = self.relaxation_modulus(t);
        let nu = self.poisson_ratio;
        e / (2.0 * (1.0 + nu))
    }
}

// ===========================================================================
// Crushable foam plasticity
// ===========================================================================

/// Crushable foam plasticity model.
///
/// Uses a volumetric yield surface (elliptical in p-q space).
/// Suitable for metallic foams, polymeric foams under large compressive loads.
#[derive(Debug, Clone, PartialEq)]
pub struct CrushableFoam {
    /// Initial yield stress in uniaxial compression \[Pa\].
    pub sigma_c0: f64,
    /// Initial yield stress in hydrostatic compression \[Pa\].
    pub p_c0: f64,
    /// Ratio of hydrostatic tensile strength to compressive (k_t).
    pub tension_cutoff_ratio: f64,
    /// Hardening curve: pairs of (plastic volumetric strain, yield stress).
    pub hardening: Vec<(f64, f64)>,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}

impl CrushableFoam {
    /// Create a new crushable foam model.
    pub fn new(
        sigma_c0: f64,
        p_c0: f64,
        tension_cutoff_ratio: f64,
        elastic_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            sigma_c0,
            p_c0,
            tension_cutoff_ratio,
            hardening: Vec::new(),
            elastic_modulus,
            poisson_ratio,
        }
    }

    /// Add a point to the hardening curve.
    pub fn add_hardening_point(&mut self, plastic_vol_strain: f64, yield_stress: f64) {
        self.hardening.push((plastic_vol_strain, yield_stress));
    }

    /// Yield surface value f(p, q) for the crushable foam model.
    ///
    /// f = q^2 + alpha^2 * (p - p0)^2 - alpha^2 * (p_c - p0)^2
    /// where p0 is the center and alpha is the shape factor.
    pub fn yield_function(&self, pressure: f64, von_mises: f64) -> f64 {
        let alpha = self.shape_factor();
        let p0 = self.yield_surface_center();
        let p_c = self.p_c0;
        von_mises * von_mises + alpha * alpha * (pressure - p0).powi(2)
            - alpha * alpha * (p_c - p0).powi(2)
    }

    /// Shape factor alpha = q_c / p_c.
    pub fn shape_factor(&self) -> f64 {
        if self.p_c0.abs() < 1.0e-30 {
            return 1.0;
        }
        self.sigma_c0 / self.p_c0
    }

    /// Center of the yield ellipse in pressure space.
    pub fn yield_surface_center(&self) -> f64 {
        let k_t = self.tension_cutoff_ratio;
        // p_0 = (p_c - k_t * p_c) / 2 = p_c * (1 - k_t) / 2
        self.p_c0 * (1.0 - k_t) / 2.0
    }

    /// Check if a stress state (p, q) is inside the yield surface.
    pub fn is_elastic(&self, pressure: f64, von_mises: f64) -> bool {
        self.yield_function(pressure, von_mises) < 0.0
    }

    /// Current yield stress from hardening table by interpolation.
    pub fn current_yield_stress(&self, plastic_vol_strain: f64) -> f64 {
        if self.hardening.is_empty() {
            return self.sigma_c0;
        }
        if plastic_vol_strain <= self.hardening[0].0 {
            return self.hardening[0].1;
        }
        let last = self.hardening.len() - 1;
        if plastic_vol_strain >= self.hardening[last].0 {
            return self.hardening[last].1;
        }
        // Linear interpolation
        for i in 0..last {
            let (e0, s0) = self.hardening[i];
            let (e1, s1) = self.hardening[i + 1];
            if plastic_vol_strain >= e0 && plastic_vol_strain <= e1 {
                let t = (plastic_vol_strain - e0) / (e1 - e0);
                return s0 + t * (s1 - s0);
            }
        }
        self.sigma_c0
    }

    /// Bulk modulus.
    pub fn bulk_modulus(&self) -> f64 {
        self.elastic_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus.
    pub fn shear_modulus(&self) -> f64 {
        self.elastic_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Plastic Poisson's ratio (typically ~0 for foams).
    pub fn plastic_poisson_ratio(&self) -> f64 {
        0.0
    }
}

// ===========================================================================
// Foaming kinetics
// ===========================================================================

/// Foaming kinetics model.
///
/// Describes the nucleation, growth, and stabilization of bubbles
/// during foam processing.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamingKinetics {
    /// Initial gas concentration \[mol/m³\].
    pub c0_gas: f64,
    /// Surface tension of the foaming liquid \[N/m\].
    pub surface_tension: f64,
    /// Viscosity of the foaming liquid \[Pa·s\].
    pub viscosity: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Gas constant R = 8.314 J/(mol·K).
    pub gas_constant: f64,
    /// Activation energy for nucleation \[J/mol\].
    pub activation_energy: f64,
}

impl FoamingKinetics {
    /// Create a new foaming kinetics model.
    pub fn new(
        c0_gas: f64,
        surface_tension: f64,
        viscosity: f64,
        temperature: f64,
        activation_energy: f64,
    ) -> Self {
        Self {
            c0_gas,
            surface_tension,
            viscosity,
            temperature,
            gas_constant: 8.314,
            activation_energy,
        }
    }

    /// Classical nucleation theory: critical bubble radius.
    ///
    /// R_crit = 2 * gamma / delta_P
    pub fn critical_radius(&self, pressure_drop: f64) -> f64 {
        if pressure_drop.abs() < 1.0e-30 {
            return f64::INFINITY;
        }
        2.0 * self.surface_tension / pressure_drop
    }

    /// Nucleation energy barrier.
    ///
    /// W_crit = (16 * pi * gamma^3) / (3 * delta_P^2)
    pub fn nucleation_barrier(&self, pressure_drop: f64) -> f64 {
        if pressure_drop.abs() < 1.0e-30 {
            return f64::INFINITY;
        }
        16.0 * PI * self.surface_tension.powi(3) / (3.0 * pressure_drop.powi(2))
    }

    /// Nucleation rate (classical nucleation theory).
    ///
    /// J = J0 * exp(-W_crit / (k_B * T))
    /// Using k_B = R/N_A, simplified here with prefactor.
    pub fn nucleation_rate(&self, pressure_drop: f64, prefactor: f64) -> f64 {
        let w_crit = self.nucleation_barrier(pressure_drop);
        let k_b_t = 1.380649e-23 * self.temperature; // Boltzmann * T
        prefactor * (-w_crit / k_b_t).exp()
    }

    /// Bubble growth by diffusion (Epstein-Plesset model, simplified).
    ///
    /// R(t) = sqrt(R0^2 + 2 * D * c * t / rho_g)
    /// D is diffusivity, c is supersaturation, rho_g is gas density.
    pub fn bubble_radius_diffusion(
        &self,
        r0: f64,
        diffusivity: f64,
        supersaturation: f64,
        gas_density: f64,
        time: f64,
    ) -> f64 {
        if gas_density.abs() < 1.0e-30 {
            return r0;
        }
        let r_sq = r0 * r0 + 2.0 * diffusivity * supersaturation * time / gas_density;
        if r_sq < 0.0 {
            return r0;
        }
        r_sq.sqrt()
    }

    /// Viscous limitation on bubble growth rate.
    ///
    /// dR/dt = R * delta_P / (4 * mu)
    /// where mu is viscosity and delta_P is the internal pressure excess.
    pub fn growth_rate_viscous(&self, radius: f64, pressure_excess: f64) -> f64 {
        if self.viscosity.abs() < 1.0e-30 {
            return 0.0;
        }
        radius * pressure_excess / (4.0 * self.viscosity)
    }

    /// Avrami equation for volume fraction of foam as function of time.
    ///
    /// X(t) = 1 - exp(-k * t^n)
    pub fn avrami_fraction(&self, k: f64, n: f64, time: f64) -> f64 {
        1.0 - (-k * time.powf(n)).exp()
    }

    /// Foam expansion ratio at equilibrium.
    ///
    /// phi = V_foam / V_liquid ≈ 1 + (c0 * R * T) / P_atm
    pub fn expansion_ratio(&self, atmospheric_pressure: f64) -> f64 {
        if atmospheric_pressure.abs() < 1.0e-30 {
            return 1.0;
        }
        1.0 + self.c0_gas * self.gas_constant * self.temperature / atmospheric_pressure
    }
}

// ===========================================================================
// Thermal conductivity of foams
// ===========================================================================

/// Thermal conductivity model for foam materials.
///
/// Accounts for contributions from solid conduction, gas conduction,
/// radiation, and convection.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamThermalConductivity {
    /// Thermal conductivity of solid cell wall material \[W/(m·K)\].
    pub k_solid: f64,
    /// Thermal conductivity of the gas in cells \[W/(m·K)\].
    pub k_gas: f64,
    /// Relative density of the foam.
    pub relative_density: f64,
    /// Mean cell diameter \[m\].
    pub cell_diameter: f64,
    /// Foam temperature \[K\].
    pub temperature: f64,
    /// Cell type (open/closed).
    pub cell_type: CellType,
    /// Emissivity of cell walls (for radiation).
    pub emissivity: f64,
}

impl FoamThermalConductivity {
    /// Create a new foam thermal conductivity model.
    pub fn new(
        k_solid: f64,
        k_gas: f64,
        relative_density: f64,
        cell_diameter: f64,
        temperature: f64,
        cell_type: CellType,
        emissivity: f64,
    ) -> Self {
        Self {
            k_solid,
            k_gas,
            relative_density,
            cell_diameter,
            temperature,
            cell_type,
            emissivity,
        }
    }

    /// Solid conduction contribution.
    ///
    /// k_s_eff = (2/3) * rho_rel * k_solid   (open-cell)
    /// k_s_eff = rho_rel * k_solid            (closed-cell, simplified)
    pub fn solid_conduction(&self) -> f64 {
        match self.cell_type {
            CellType::Open => (2.0 / 3.0) * self.relative_density * self.k_solid,
            CellType::Closed | CellType::Mixed => self.relative_density * self.k_solid,
        }
    }

    /// Gas conduction contribution.
    ///
    /// k_g_eff = (1 - rho_rel) * k_gas
    pub fn gas_conduction(&self) -> f64 {
        (1.0 - self.relative_density) * self.k_gas
    }

    /// Radiation contribution (Rosseland approximation).
    ///
    /// k_rad = 16 * sigma_SB * T^3 * d / (3 * K_ext)
    /// where K_ext ≈ (1 - porosity) / d for simplicity.
    pub fn radiation_contribution(&self) -> f64 {
        let sigma_sb = 5.670374419e-8; // Stefan-Boltzmann
        let porosity = 1.0 - self.relative_density;
        if self.cell_diameter.abs() < 1.0e-30 || porosity.abs() < 1.0e-30 {
            return 0.0;
        }
        // Extinction coefficient approximation
        let k_ext = self.relative_density / self.cell_diameter;
        if k_ext.abs() < 1.0e-30 {
            return 0.0;
        }
        16.0 * sigma_sb * self.emissivity * self.temperature.powi(3) * self.cell_diameter
            / (3.0 * k_ext * self.cell_diameter)
    }

    /// Total effective thermal conductivity.
    pub fn total_conductivity(&self) -> f64 {
        self.solid_conduction() + self.gas_conduction() + self.radiation_contribution()
    }

    /// Knudsen effect correction for gas conductivity in small cells.
    ///
    /// k_gas_eff = k_gas / (1 + 2 * beta * Kn)
    /// where Kn = lambda / d is the Knudsen number,
    /// lambda is the mean free path, beta ≈ 1.64.
    pub fn knudsen_corrected_gas(&self, mean_free_path: f64) -> f64 {
        if self.cell_diameter.abs() < 1.0e-30 {
            return 0.0;
        }
        let kn = mean_free_path / self.cell_diameter;
        let beta = 1.64;
        self.k_gas / (1.0 + 2.0 * beta * kn)
    }

    /// R-value (thermal resistance) for a given foam thickness.
    ///
    /// R = thickness / k_total \[m²·K/W\]
    pub fn r_value(&self, thickness: f64) -> f64 {
        let k = self.total_conductivity();
        if k.abs() < 1.0e-30 {
            return f64::INFINITY;
        }
        thickness / k
    }
}

// ===========================================================================
// Impact attenuation
// ===========================================================================

/// Impact attenuation model for foam packaging/helmets.
///
/// Predicts peak acceleration, HIC, and energy absorption for
/// drop impact scenarios.
#[derive(Debug, Clone, PartialEq)]
pub struct ImpactAttenuation {
    /// Foam plateau stress \[Pa\].
    pub plateau_stress: f64,
    /// Foam densification strain.
    pub densification_strain: f64,
    /// Foam thickness \[m\].
    pub thickness: f64,
    /// Foam contact area \[m²\].
    pub contact_area: f64,
    /// Object mass \[kg\].
    pub impactor_mass: f64,
}

impl ImpactAttenuation {
    /// Create a new impact attenuation model.
    pub fn new(
        plateau_stress: f64,
        densification_strain: f64,
        thickness: f64,
        contact_area: f64,
        impactor_mass: f64,
    ) -> Self {
        Self {
            plateau_stress,
            densification_strain,
            thickness,
            contact_area,
            impactor_mass,
        }
    }

    /// Maximum energy the foam can absorb before densification \[J\].
    pub fn max_energy_capacity(&self) -> f64 {
        self.plateau_stress * self.contact_area * self.thickness * self.densification_strain
    }

    /// Impact velocity that saturates the foam capacity \[m/s\].
    ///
    /// 0.5 * m * v^2 = E_max  =>  v = sqrt(2 * E_max / m)
    pub fn critical_velocity(&self) -> f64 {
        if self.impactor_mass.abs() < 1.0e-30 {
            return 0.0;
        }
        let e_max = self.max_energy_capacity();
        (2.0 * e_max / self.impactor_mass).sqrt()
    }

    /// Peak deceleration in g's for a given impact velocity.
    ///
    /// a_peak = sigma_p * A / (m * g)   when v < v_crit
    pub fn peak_g(&self, _impact_velocity: f64) -> f64 {
        let g = 9.81;
        if self.impactor_mass.abs() < 1.0e-30 {
            return 0.0;
        }
        self.plateau_stress * self.contact_area / (self.impactor_mass * g)
    }

    /// Head Injury Criterion (HIC) estimate.
    ///
    /// HIC = ((1/(t2-t1)) * integral_t1^t2 a(t) dt)^2.5 * (t2 - t1)
    /// For constant deceleration a: HIC = a^2.5 * dt
    pub fn hic_estimate(&self, impact_velocity: f64) -> f64 {
        let a_g = self.peak_g(impact_velocity);
        if a_g.abs() < 1.0e-30 || self.plateau_stress.abs() < 1.0e-30 {
            return 0.0;
        }
        // Impact duration: dt = m * v / (sigma_p * A)
        let dt = self.impactor_mass * impact_velocity / (self.plateau_stress * self.contact_area);
        // HIC = (a_g * g)^2.5 * dt  (with acceleration in m/s^2)
        let a_ms2 = a_g * 9.81;
        a_ms2.powf(2.5) * dt
    }

    /// Compression depth for a given impact velocity \[m\].
    ///
    /// delta = 0.5 * m * v^2 / (sigma_p * A)
    pub fn compression_depth(&self, impact_velocity: f64) -> f64 {
        if self.plateau_stress.abs() < 1.0e-30 || self.contact_area.abs() < 1.0e-30 {
            return 0.0;
        }
        0.5 * self.impactor_mass * impact_velocity * impact_velocity
            / (self.plateau_stress * self.contact_area)
    }

    /// Whether the foam is fully crushed at a given velocity.
    pub fn is_bottomed_out(&self, impact_velocity: f64) -> bool {
        let depth = self.compression_depth(impact_velocity);
        depth >= self.thickness * self.densification_strain
    }

    /// Rebound velocity (coefficient of restitution ≈ 0 for ideal foam).
    ///
    /// For a real foam with given COR:
    ///   v_rebound = cor * v_impact
    pub fn rebound_velocity(&self, impact_velocity: f64, cor: f64) -> f64 {
        cor * impact_velocity
    }

    /// Drop height corresponding to a given impact velocity.
    ///
    /// h = v^2 / (2 * g)
    pub fn drop_height_from_velocity(impact_velocity: f64) -> f64 {
        impact_velocity * impact_velocity / (2.0 * 9.81)
    }

    /// Impact velocity from a given drop height.
    ///
    /// v = sqrt(2 * g * h)
    pub fn velocity_from_drop_height(height: f64) -> f64 {
        (2.0 * 9.81 * height).sqrt()
    }
}

// ===========================================================================
// Foam aging / degradation
// ===========================================================================

/// Simple foam aging/degradation model.
///
/// Models the reduction in foam properties over time due to gas diffusion,
/// cell wall degradation, and UV exposure.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamAging {
    /// Initial modulus \[Pa\].
    pub e_initial: f64,
    /// Initial strength \[Pa\].
    pub sigma_initial: f64,
    /// Gas loss rate constant \[1/s\].
    pub gas_loss_rate: f64,
    /// Modulus degradation rate constant \[1/s\].
    pub modulus_degradation_rate: f64,
    /// Strength degradation rate constant \[1/s\].
    pub strength_degradation_rate: f64,
}

impl FoamAging {
    /// Create a new foam aging model.
    pub fn new(
        e_initial: f64,
        sigma_initial: f64,
        gas_loss_rate: f64,
        modulus_degradation_rate: f64,
        strength_degradation_rate: f64,
    ) -> Self {
        Self {
            e_initial,
            sigma_initial,
            gas_loss_rate,
            modulus_degradation_rate,
            strength_degradation_rate,
        }
    }

    /// Gas retention fraction at time t.
    pub fn gas_retention(&self, time: f64) -> f64 {
        (-self.gas_loss_rate * time).exp()
    }

    /// Modulus at time t with aging.
    pub fn modulus_at_time(&self, time: f64) -> f64 {
        self.e_initial * (-self.modulus_degradation_rate * time).exp()
    }

    /// Strength at time t with aging.
    pub fn strength_at_time(&self, time: f64) -> f64 {
        self.sigma_initial * (-self.strength_degradation_rate * time).exp()
    }

    /// Half-life of the modulus.
    pub fn modulus_half_life(&self) -> f64 {
        if self.modulus_degradation_rate.abs() < 1.0e-30 {
            return f64::INFINITY;
        }
        2.0_f64.ln() / self.modulus_degradation_rate
    }
}

// ===========================================================================
// Foam cell size distribution
// ===========================================================================

/// Foam cell size distribution statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct CellSizeDistribution {
    /// Mean cell diameter \[m\].
    pub mean_diameter: f64,
    /// Standard deviation of cell diameter \[m\].
    pub std_diameter: f64,
    /// Number of cells per unit volume \[1/m³\].
    pub cell_density: f64,
}

impl CellSizeDistribution {
    /// Create a new cell size distribution.
    pub fn new(mean_diameter: f64, std_diameter: f64, cell_density: f64) -> Self {
        Self {
            mean_diameter,
            std_diameter,
            cell_density,
        }
    }

    /// Coefficient of variation.
    pub fn coefficient_of_variation(&self) -> f64 {
        if self.mean_diameter.abs() < 1.0e-30 {
            return 0.0;
        }
        self.std_diameter / self.mean_diameter
    }

    /// Approximate porosity from cell density and mean diameter.
    ///
    /// phi ≈ N * (pi/6) * d^3, capped at 1.0.
    pub fn estimated_porosity(&self) -> f64 {
        let vol_per_cell = PI / 6.0 * self.mean_diameter.powi(3);
        (self.cell_density * vol_per_cell).min(1.0)
    }

    /// Specific surface area (surface per unit volume) for spherical cells.
    ///
    /// S_v = N * pi * d^2
    pub fn specific_surface_area(&self) -> f64 {
        self.cell_density * PI * self.mean_diameter.powi(2)
    }

    /// Characteristic cell wall thickness for closed-cell foam.
    ///
    /// t ≈ d * (1 - (1 - rho_rel)^(1/3))
    pub fn wall_thickness(&self, relative_density: f64) -> f64 {
        let x = 1.0 - relative_density;
        if x <= 0.0 {
            return self.mean_diameter;
        }
        self.mean_diameter * (1.0 - x.powf(1.0 / 3.0))
    }
}

// ===========================================================================
// Multi-layer foam stack
// ===========================================================================

/// Multi-layer foam stack for graded energy absorption.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamStack {
    /// Layers: each described by (thickness, plateau_stress, densification_strain, modulus).
    pub layers: Vec<FoamLayer>,
}

/// A single layer in a foam stack.
#[derive(Debug, Clone, PartialEq)]
pub struct FoamLayer {
    /// Thickness \[m\].
    pub thickness: f64,
    /// Plateau stress \[Pa\].
    pub plateau_stress: f64,
    /// Densification strain.
    pub densification_strain: f64,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
}

impl FoamLayer {
    /// Create a new foam layer.
    pub fn new(
        thickness: f64,
        plateau_stress: f64,
        densification_strain: f64,
        elastic_modulus: f64,
    ) -> Self {
        Self {
            thickness,
            plateau_stress,
            densification_strain,
            elastic_modulus,
        }
    }

    /// Maximum energy this layer can absorb \[J/m²\].
    pub fn max_energy_per_area(&self) -> f64 {
        self.plateau_stress * self.thickness * self.densification_strain
    }
}

impl FoamStack {
    /// Create a new empty foam stack.
    pub fn new() -> Self {
        Self { layers: Vec::new() }
    }

    /// Add a layer to the stack.
    pub fn add_layer(&mut self, layer: FoamLayer) {
        self.layers.push(layer);
    }

    /// Total thickness of the stack.
    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// Total maximum energy absorption per unit area \[J/m²\].
    pub fn total_energy_capacity(&self) -> f64 {
        self.layers.iter().map(|l| l.max_energy_per_area()).sum()
    }

    /// Effective (series) modulus of the stack.
    ///
    /// 1/E_eff = sum (t_i / E_i) / sum(t_i)
    pub fn effective_modulus(&self) -> f64 {
        let total_t = self.total_thickness();
        if total_t.abs() < 1.0e-30 {
            return 0.0;
        }
        let compliance: f64 = self
            .layers
            .iter()
            .map(|l| {
                if l.elastic_modulus.abs() < 1.0e-30 {
                    0.0
                } else {
                    l.thickness / l.elastic_modulus
                }
            })
            .sum();
        if compliance.abs() < 1.0e-30 {
            return 0.0;
        }
        total_t / compliance
    }

    /// Number of layers.
    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }
}

impl Default for FoamStack {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-6;

    // 1. Gibson-Ashby open-cell modulus scaling
    #[test]
    fn test_gibson_ashby_modulus_open() {
        let ga = GibsonAshby::new(1000.0, 1.0e9, 1.0e6, 0.2, CellType::Open, 0.0);
        let rho_foam = 100.0; // relative density = 0.1
        let e = ga.modulus_open(rho_foam);
        let expected = 1.0e9 * 0.01; // (0.1)^2
        assert!(
            (e - expected).abs() < EPS * expected,
            "Open-cell modulus: got {e}, expected {expected}"
        );
    }

    // 2. Gibson-Ashby closed-cell modulus
    #[test]
    fn test_gibson_ashby_modulus_closed() {
        let ga = GibsonAshby::new(1000.0, 1.0e9, 1.0e6, 0.2, CellType::Closed, 0.0);
        let rho_foam = 100.0;
        let e = ga.modulus_closed(rho_foam, 0.6);
        let rd = 0.1;
        let expected = 1.0e9 * (0.36 * rd * rd + 0.4 * rd);
        assert!(
            (e - expected).abs() < EPS * expected,
            "Closed-cell modulus: got {e}, expected {expected}"
        );
    }

    // 3. Gibson-Ashby open-cell strength
    #[test]
    fn test_gibson_ashby_strength_open() {
        let ga = GibsonAshby::new(1000.0, 1.0e9, 1.0e6, 0.2, CellType::Open, 0.0);
        let rho_foam = 100.0;
        let s = ga.strength_open(rho_foam);
        let expected = 0.3 * 1.0e6 * 0.1_f64.powf(1.5);
        assert!(
            (s - expected).abs() < EPS * expected,
            "Open-cell strength: got {s}, expected {expected}"
        );
    }

    // 4. Kelvin cell volume
    #[test]
    fn test_kelvin_cell_volume() {
        let kc = KelvinCell::new(0.001, 0.0001, 1000.0, 1.0e9);
        let v = kc.cell_volume();
        let expected = 8.0 * 2.0_f64.sqrt() * 0.001_f64.powi(3);
        assert!(
            (v - expected).abs() < EPS * expected,
            "Kelvin cell volume: got {v}, expected {expected}"
        );
    }

    // 5. Kelvin cell face count
    #[test]
    fn test_kelvin_cell_faces() {
        let kc = KelvinCell::new(0.001, 0.0001, 1000.0, 1.0e9);
        assert_eq!(kc.num_faces(), 14);
        assert_eq!(kc.num_edges(), 36);
        assert_eq!(kc.num_vertices(), 24);
    }

    // 6. Energy absorption - elastic region
    #[test]
    fn test_energy_absorption_elastic() {
        let ea = EnergyAbsorption::new(1.0e5, 0.7, 1.0e6, 0.05, 50.0, 0.05);
        let strain = 0.03;
        let stress = ea.stress_at_strain(strain);
        let expected_stress = 1.0e6 * 0.03;
        assert!(
            (stress - expected_stress).abs() < EPS,
            "Elastic stress: got {stress}"
        );
    }

    // 7. Energy absorption - plateau region
    #[test]
    fn test_energy_absorption_plateau() {
        let ea = EnergyAbsorption::new(1.0e5, 0.7, 1.0e6, 0.05, 50.0, 0.05);
        let strain = 0.3; // in plateau region
        let stress = ea.stress_at_strain(strain);
        assert!(
            (stress - 1.0e5).abs() < EPS,
            "Plateau stress should be {}, got {stress}",
            1.0e5
        );
    }

    // 8. Energy absorption - densification rises
    #[test]
    fn test_energy_absorption_densification() {
        let ea = EnergyAbsorption::new(1.0e5, 0.7, 1.0e6, 0.05, 50.0, 0.05);
        let stress_plateau = ea.stress_at_strain(0.7);
        let stress_dense = ea.stress_at_strain(0.8);
        assert!(
            stress_dense > stress_plateau,
            "Densification stress should exceed plateau"
        );
    }

    // 9. Energy density monotonically increases
    #[test]
    fn test_energy_density_monotonic() {
        let ea = EnergyAbsorption::new(1.0e5, 0.7, 1.0e6, 0.05, 50.0, 0.05);
        let mut prev = 0.0;
        for i in 1..=20 {
            let strain = i as f64 * 0.05;
            let w = ea.energy_density(strain);
            assert!(
                w >= prev,
                "Energy density should be monotonic at strain={strain}"
            );
            prev = w;
        }
    }

    // 10. Viscoelastic foam relaxation modulus at t=0
    #[test]
    fn test_viscoelastic_modulus_t0() {
        let foam = ViscoelasticFoam::single_term(1.0e6, 0.5e6, 1.0, 50.0);
        let e0 = foam.relaxation_modulus(0.0);
        // At t=0: E = E_inf + E_0 * g = E_inf + (E_0 - E_inf) = E_0
        assert!(
            (e0 - 1.0e6).abs() < EPS,
            "At t=0, modulus should be E_instantaneous, got {e0}"
        );
    }

    // 11. Viscoelastic foam relaxation modulus at t→∞
    #[test]
    fn test_viscoelastic_modulus_inf() {
        let foam = ViscoelasticFoam::single_term(1.0e6, 0.5e6, 1.0, 50.0);
        let e_inf = foam.relaxation_modulus(1.0e6);
        assert!(
            (e_inf - 0.5e6).abs() < 1.0,
            "At t→∞, modulus should be E_equilibrium, got {e_inf}"
        );
    }

    // 12. Viscoelastic loss tangent is positive
    #[test]
    fn test_viscoelastic_loss_tangent() {
        let foam = ViscoelasticFoam::single_term(1.0e6, 0.5e6, 1.0, 50.0);
        let tan_delta = foam.loss_tangent(1.0);
        assert!(
            tan_delta > 0.0,
            "Loss tangent should be positive, got {tan_delta}"
        );
    }

    // 13. Crushable foam yield surface
    #[test]
    fn test_crushable_foam_yield() {
        let cf = CrushableFoam::new(1.0e5, 2.0e5, 0.1, 1.0e7, 0.3);
        // At the origin, should be elastic
        let f = cf.yield_function(0.0, 0.0);
        assert!(f < 0.0, "Origin should be inside yield surface, f={f}");
    }

    // 14. Crushable foam - large stress is plastic
    #[test]
    fn test_crushable_foam_plastic() {
        let cf = CrushableFoam::new(1.0e5, 2.0e5, 0.1, 1.0e7, 0.3);
        let f = cf.yield_function(1.0e6, 1.0e6);
        assert!(
            f > 0.0,
            "Large stress should be outside yield surface, f={f}"
        );
    }

    // 15. Foaming kinetics critical radius
    #[test]
    fn test_nucleation_critical_radius() {
        let fk = FoamingKinetics::new(100.0, 0.03, 1.0, 400.0, 50000.0);
        let r_crit = fk.critical_radius(1.0e5);
        let expected = 2.0 * 0.03 / 1.0e5;
        assert!(
            (r_crit - expected).abs() < EPS * expected,
            "Critical radius: got {r_crit}, expected {expected}"
        );
    }

    // 16. Foaming kinetics Avrami
    #[test]
    fn test_avrami_limits() {
        let fk = FoamingKinetics::new(100.0, 0.03, 1.0, 400.0, 50000.0);
        let x0 = fk.avrami_fraction(0.1, 2.0, 0.0);
        assert!((x0).abs() < EPS, "Avrami at t=0 should be 0, got {x0}");
        let x_inf = fk.avrami_fraction(0.1, 2.0, 100.0);
        assert!(
            (x_inf - 1.0).abs() < 0.01,
            "Avrami at t→∞ should be ~1, got {x_inf}"
        );
    }

    // 17. Thermal conductivity gas + solid
    #[test]
    fn test_thermal_conductivity_components() {
        let tc = FoamThermalConductivity::new(0.2, 0.026, 0.05, 0.001, 300.0, CellType::Open, 0.9);
        let k_s = tc.solid_conduction();
        let k_g = tc.gas_conduction();
        assert!(k_s > 0.0, "Solid conduction should be positive");
        assert!(k_g > 0.0, "Gas conduction should be positive");
        assert!(k_g > k_s, "For low-density foam, gas conduction dominates");
    }

    // 18. R-value increases with thickness
    #[test]
    fn test_r_value_thickness() {
        let tc =
            FoamThermalConductivity::new(0.2, 0.026, 0.05, 0.001, 300.0, CellType::Closed, 0.9);
        let r1 = tc.r_value(0.05);
        let r2 = tc.r_value(0.10);
        assert!(
            r2 > r1,
            "R-value should increase with thickness: r1={r1}, r2={r2}"
        );
        assert!(
            (r2 / r1 - 2.0).abs() < 0.01,
            "R-value should double when thickness doubles"
        );
    }

    // 19. Impact attenuation critical velocity
    #[test]
    fn test_impact_critical_velocity() {
        let ia = ImpactAttenuation::new(1.0e5, 0.7, 0.05, 0.01, 5.0);
        let v_crit = ia.critical_velocity();
        let e_max = ia.max_energy_capacity();
        let expected = (2.0 * e_max / 5.0).sqrt();
        assert!(
            (v_crit - expected).abs() < EPS,
            "Critical velocity: got {v_crit}, expected {expected}"
        );
    }

    // 20. Impact bottoming out detection
    #[test]
    fn test_impact_bottoming_out() {
        let ia = ImpactAttenuation::new(1.0e5, 0.7, 0.05, 0.01, 5.0);
        let v_low = 0.1;
        let v_high = 100.0;
        assert!(
            !ia.is_bottomed_out(v_low),
            "Low velocity should not bottom out"
        );
        assert!(
            ia.is_bottomed_out(v_high),
            "High velocity should bottom out"
        );
    }

    // 21. Foam aging modulus decreases over time
    #[test]
    fn test_foam_aging() {
        let aging = FoamAging::new(1.0e6, 1.0e5, 1.0e-8, 1.0e-7, 1.0e-7);
        let e0 = aging.modulus_at_time(0.0);
        let e1 = aging.modulus_at_time(1.0e6);
        assert!(
            (e0 - 1.0e6).abs() < EPS,
            "Initial modulus should be E_initial"
        );
        assert!(e1 < e0, "Modulus should decrease with aging");
    }

    // 22. Cell size distribution porosity
    #[test]
    fn test_cell_size_porosity() {
        let csd = CellSizeDistribution::new(0.001, 0.0002, 1.0e8);
        let phi = csd.estimated_porosity();
        assert!(
            (0.0..=1.0).contains(&phi),
            "Porosity should be in [0,1], got {phi}"
        );
    }

    // 23. Foam stack energy capacity
    #[test]
    fn test_foam_stack_energy() {
        let mut stack = FoamStack::new();
        stack.add_layer(FoamLayer::new(0.02, 5.0e4, 0.7, 5.0e5));
        stack.add_layer(FoamLayer::new(0.03, 1.0e5, 0.6, 1.0e6));
        let e_total = stack.total_energy_capacity();
        let e1 = 5.0e4 * 0.02 * 0.7;
        let e2 = 1.0e5 * 0.03 * 0.6;
        assert!(
            (e_total - (e1 + e2)).abs() < EPS,
            "Stack energy should be sum of layers"
        );
    }

    // 24. Foam stack effective modulus (series)
    #[test]
    fn test_foam_stack_modulus() {
        let mut stack = FoamStack::new();
        stack.add_layer(FoamLayer::new(0.01, 1.0e5, 0.7, 1.0e6));
        stack.add_layer(FoamLayer::new(0.01, 1.0e5, 0.7, 1.0e6));
        // Two identical layers → effective modulus = single layer modulus
        let e_eff = stack.effective_modulus();
        assert!(
            (e_eff - 1.0e6).abs() < EPS,
            "Two identical layers: E_eff should equal E_layer, got {e_eff}"
        );
    }

    // 25. Gibson-Ashby relative density
    #[test]
    fn test_relative_density() {
        let ga = GibsonAshby::new(1000.0, 1.0e9, 1.0e6, 0.2, CellType::Open, 0.0);
        let rd = ga.relative_density(100.0);
        assert!(
            (rd - 0.1).abs() < EPS,
            "Relative density should be 0.1, got {rd}"
        );
    }

    // 26. Viscoelastic dynamic moduli at zero frequency
    #[test]
    fn test_dynamic_moduli_zero_freq() {
        let foam = ViscoelasticFoam::single_term(1.0e6, 0.5e6, 1.0, 50.0);
        let (e_s, e_l) = foam.dynamic_moduli(0.0);
        assert!(
            (e_s - 0.5e6).abs() < 1.0,
            "Storage modulus at omega=0 should equal E_equilibrium, got {e_s}"
        );
        assert!(
            e_l.abs() < EPS,
            "Loss modulus at omega=0 should be zero, got {e_l}"
        );
    }

    // 27. Kelvin cell relative density bounded
    #[test]
    fn test_kelvin_relative_density_bounded() {
        let kc = KelvinCell::new(0.001, 0.5, 1000.0, 1.0e9); // very thick walls
        let rd = kc.relative_density();
        assert!(rd <= 1.0, "Relative density capped at 1.0, got {rd}");
    }

    // 28. Energy absorption efficiency in plateau
    #[test]
    fn test_efficiency_in_plateau() {
        let ea = EnergyAbsorption::new(1.0e5, 0.7, 1.0e6, 0.05, 50.0, 0.05);
        let eff = ea.efficiency(0.4); // well into plateau
        assert!(eff > 0.5, "Efficiency in plateau should be high, got {eff}");
        assert!(eff <= 1.0, "Efficiency should not exceed 1.0, got {eff}");
    }
}
