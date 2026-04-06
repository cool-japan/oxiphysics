// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biomaterials module: bone mechanics, cartilage, tendons, hydrogels, cell mechanics,
//! biodegradation kinetics, scaffold porosity, sutures, and dental materials.
//!
//! All functions use SI units unless otherwise stated.

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Bone Mechanics — Cortical and Cancellous
// ---------------------------------------------------------------------------

/// Properties of cortical bone (compact bone).
#[derive(Clone, Debug)]
pub struct CorticalBone {
    /// Longitudinal Young's modulus (Pa), typically ~17-25 GPa.
    pub e_long: f64,
    /// Transverse Young's modulus (Pa), typically ~10-15 GPa.
    pub e_trans: f64,
    /// Shear modulus (Pa).
    pub g_lt: f64,
    /// Longitudinal Poisson's ratio.
    pub nu_lt: f64,
    /// Transverse Poisson's ratio.
    pub nu_tt: f64,
    /// Density (kg/m³), typically ~1800-2000 kg/m³.
    pub density: f64,
    /// Ultimate tensile strength (Pa).
    pub uts: f64,
    /// Ultimate compressive strength (Pa).
    pub ucs: f64,
    /// Fracture toughness KIC (Pa·√m).
    pub k_ic: f64,
}

impl CorticalBone {
    /// Default human femoral cortical bone properties (alias).
    pub fn femur_cortical() -> Self {
        Self::human_femur()
    }

    /// Default human femoral cortical bone properties.
    pub fn human_femur() -> Self {
        Self {
            e_long: 20.0e9,
            e_trans: 12.0e9,
            g_lt: 4.5e9,
            nu_lt: 0.29,
            nu_tt: 0.63,
            density: 1900.0,
            uts: 120.0e6,
            ucs: 180.0e6,
            k_ic: 2.5e6,
        }
    }

    /// Compute the longitudinal wave speed (m/s).
    pub fn longitudinal_wave_speed(&self) -> f64 {
        (self.e_long / self.density).sqrt()
    }

    /// Compute the transverse wave speed (m/s).
    pub fn transverse_wave_speed(&self) -> f64 {
        (self.g_lt / self.density).sqrt()
    }

    /// Estimate the apparent density from volumetric fraction `vf` (0..1).
    pub fn apparent_density(vf: f64) -> f64 {
        // Rice's rule: rho_app = vf * rho_cortical
        vf * 1900.0
    }

    /// Morgan-Keaveny power-law: E = a * rho^b (GPa), rho in g/cm³.
    pub fn modulus_from_density_power_law(density_g_cm3: f64) -> f64 {
        // Typical empirical fit for cortical bone
        let a = 2.065e10;
        let b = 3.09;
        a * density_g_cm3.powf(b)
    }

    /// Compute bone stiffness in bending for a hollow cylinder (Euler-Bernoulli).
    ///
    /// `outer_r`: outer radius (m). `inner_r`: inner radius (m). `length`: beam length (m).
    pub fn bending_stiffness(&self, outer_r: f64, inner_r: f64, length: f64) -> f64 {
        let i = PI / 4.0 * (outer_r.powi(4) - inner_r.powi(4));
        48.0 * self.e_long * i / length.powi(3)
    }

    /// Stress at outer fibre for three-point bending.
    pub fn three_point_bending_stress(
        &self,
        force: f64,
        length: f64,
        outer_r: f64,
        inner_r: f64,
    ) -> f64 {
        let i = PI / 4.0 * (outer_r.powi(4) - inner_r.powi(4));
        force * length * outer_r / (4.0 * i)
    }
}

/// Properties of cancellous (trabecular) bone.
#[derive(Clone, Debug)]
pub struct CancellousBone {
    /// Relative density (ratio to cortical, 0..1).
    pub relative_density: f64,
    /// Apparent density (kg/m³).
    pub apparent_density: f64,
    /// Apparent elastic modulus (Pa).
    pub apparent_modulus: f64,
    /// Apparent strength (Pa).
    pub apparent_strength: f64,
    /// Anisotropy ratio (E_axial / E_transverse).
    pub anisotropy_ratio: f64,
}

impl CancellousBone {
    /// Construct cancellous bone using Gibson-Ashby cellular solid relations.
    ///
    /// `rel_density`: relative density (0..1).
    pub fn from_relative_density(rel_density: f64) -> Self {
        // Gibson-Ashby: E* / E_s = C1 * (rho*/rho_s)^2
        let e_solid = 20.0e9;
        let sigma_solid = 170.0e6;
        let c1 = 1.0;
        let c2 = 0.2;
        let apparent_modulus = c1 * e_solid * rel_density.powi(2);
        let apparent_strength = c2 * sigma_solid * rel_density.powf(1.5);
        Self {
            relative_density: rel_density,
            apparent_density: rel_density * 1900.0,
            apparent_modulus,
            apparent_strength,
            anisotropy_ratio: 1.5,
        }
    }

    /// Compute energy absorption capacity (J/m³) under compression.
    pub fn energy_absorption(&self) -> f64 {
        // Area under stress-strain curve up to densification strain
        let densification_strain = 1.0 - 1.4 * self.relative_density;
        0.5 * self.apparent_strength * densification_strain.max(0.0)
    }

    /// Strain at failure under compression.
    pub fn failure_strain(&self) -> f64 {
        // Empirical: eps_f ~ 0.07 for rho_rel ≈ 0.2
        0.07 * (self.relative_density / 0.2).powf(-0.3)
    }
}

/// Anisotropic bone elastic tensor (Voigt notation, GPa-scaled).
///
/// Returns the 6×6 stiffness matrix in Voigt notation (upper triangular).
pub fn cortical_bone_stiffness_tensor(
    e_l: f64,
    e_t: f64,
    g_lt: f64,
    nu_lt: f64,
    nu_tt: f64,
) -> [[f64; 6]; 6] {
    // Transversely isotropic: axis of symmetry along L
    let nu_tl = nu_lt * e_t / e_l;
    let det = 1.0 - nu_tt * nu_tt - 2.0 * nu_lt * nu_tl;
    let mut c = [[0.0f64; 6]; 6];
    c[0][0] = e_l * (1.0 - nu_tt * nu_tt) / det;
    c[1][1] = e_t * (1.0 - nu_lt * nu_tl) / det;
    c[2][2] = c[1][1];
    c[0][1] = e_l * (nu_tl + nu_tt * nu_tl) / det;
    c[1][0] = c[0][1];
    c[0][2] = c[0][1];
    c[2][0] = c[0][2];
    c[1][2] = e_t * (nu_tt + nu_lt * nu_tl) / det;
    c[2][1] = c[1][2];
    c[3][3] = g_lt;
    c[4][4] = g_lt;
    c[5][5] = e_t / (2.0 * (1.0 + nu_tt));
    c
}

// ---------------------------------------------------------------------------
// Cartilage Biphasic Model (Mow et al.)
// ---------------------------------------------------------------------------

/// Articular cartilage biphasic model parameters.
#[derive(Clone, Debug)]
pub struct CartilageBiphasic {
    /// Aggregate modulus HA (Pa), typically 0.5-1.5 MPa.
    pub ha: f64,
    /// Shear modulus mu_s (Pa), typically 0.1-1.5 MPa.
    pub mu_s: f64,
    /// Hydraulic permeability k (m⁴/N·s), typically 10⁻¹⁵.
    pub permeability: f64,
    /// Thickness (m).
    pub thickness: f64,
    /// Tissue density (kg/m³).
    pub density: f64,
}

impl CartilageBiphasic {
    /// Create a model for the superficial zone of cartilage.
    pub fn superficial_zone() -> Self {
        Self {
            ha: 0.42e6,
            mu_s: 0.10e6,
            permeability: 4.0e-15,
            thickness: 2.0e-3,
            density: 1100.0,
        }
    }

    /// Create a model for the middle zone of cartilage.
    pub fn middle_zone() -> Self {
        Self {
            ha: 0.60e6,
            mu_s: 0.15e6,
            permeability: 2.0e-15,
            thickness: 2.5e-3,
            density: 1100.0,
        }
    }

    /// Create a model for the deep zone of cartilage.
    pub fn deep_zone() -> Self {
        Self {
            ha: 0.85e6,
            mu_s: 0.20e6,
            permeability: 1.0e-15,
            thickness: 3.0e-3,
            density: 1100.0,
        }
    }

    /// Human knee articular cartilage defaults.
    pub fn knee() -> Self {
        Self {
            ha: 0.7e6,
            mu_s: 0.15e6,
            permeability: 2.17e-15,
            thickness: 2.5e-3,
            density: 1100.0,
        }
    }

    /// Biphasic creep time constant tau = h^2 / (HA * k).
    pub fn creep_time_constant(&self) -> f64 {
        self.thickness * self.thickness / (self.ha * self.permeability)
    }

    /// Steady-state consolidation (creep) strain under stress `sigma`.
    pub fn steady_state_strain(&self, sigma: f64) -> f64 {
        sigma / self.ha
    }

    /// Creep displacement at time `t` for applied stress `sigma` (Mow's solution, first term).
    pub fn creep_displacement(&self, sigma: f64, t: f64) -> f64 {
        let tau = self.creep_time_constant();
        let eps_inf = self.steady_state_strain(sigma);
        // First-term approximation
        let correction = (1.0 - (-t / tau * PI * PI / 4.0).exp()).max(0.0);
        eps_inf * self.thickness * correction
    }

    /// Effective fluid pressure at surface under step load `sigma` at time `t`.
    pub fn fluid_pressure(&self, sigma: f64, t: f64) -> f64 {
        let tau = self.creep_time_constant();
        sigma * (-t / tau * PI * PI / 4.0).exp()
    }

    /// Poisson's ratio from biphasic parameters.
    pub fn poisson_ratio(&self) -> f64 {
        (self.ha - 2.0 * self.mu_s) / (2.0 * (self.ha - self.mu_s))
    }

    /// Confined compression stiffness (per unit area).
    pub fn confined_compression_stiffness(&self, area: f64) -> f64 {
        self.ha * area / self.thickness
    }
}

// ---------------------------------------------------------------------------
// Tendon / Ligament Viscoelasticity
// ---------------------------------------------------------------------------

/// Viscoelastic tendon/ligament model using Quasilinear Viscoelasticity (QLV).
#[derive(Clone, Debug)]
pub struct TendonViscoelastic {
    /// Toe-region strain (dimensionless).
    pub toe_strain: f64,
    /// Toe-region modulus (Pa).
    pub toe_modulus: f64,
    /// Linear-region modulus (Pa), typically 0.5-2 GPa.
    pub linear_modulus: f64,
    /// Relaxation time constant (s).
    pub tau: f64,
    /// Relaxation magnitude (dimensionless, 0..1).
    pub relaxation_magnitude: f64,
    /// Cross-sectional area (m²).
    pub area: f64,
    /// Resting length (m).
    pub resting_length: f64,
}

impl TendonViscoelastic {
    /// Human Achilles tendon defaults.
    pub fn achilles() -> Self {
        Self {
            toe_strain: 0.03,
            toe_modulus: 50.0e6,
            linear_modulus: 1.2e9,
            tau: 20.0,
            relaxation_magnitude: 0.15,
            area: 80.0e-6,
            resting_length: 0.25,
        }
    }

    /// Compute toe-region and linear-region stress from strain.
    pub fn elastic_stress(&self, strain: f64) -> f64 {
        if strain < 0.0 {
            return 0.0;
        }
        if strain <= self.toe_strain {
            // Exponential toe: sigma = A * (exp(B*eps) - 1)
            let b = (self.linear_modulus / self.toe_modulus).ln() / self.toe_strain;
            let a = self.toe_modulus / b;
            a * ((b * strain).exp() - 1.0)
        } else {
            let toe_stress = self.elastic_stress(self.toe_strain);
            toe_stress + self.linear_modulus * (strain - self.toe_strain)
        }
    }

    /// Reduced relaxation function G(t) = 1 - M*(1 - exp(-t/tau)).
    pub fn reduced_relaxation(&self, t: f64) -> f64 {
        1.0 - self.relaxation_magnitude * (1.0 - (-t / self.tau).exp())
    }

    /// Stress relaxation: sigma(t) = G(t) * sigma_e(eps_0).
    pub fn stress_relaxation(&self, eps0: f64, t: f64) -> f64 {
        self.elastic_stress(eps0) * self.reduced_relaxation(t)
    }

    /// Force-elongation relationship at given strain.
    pub fn force(&self, strain: f64) -> f64 {
        self.elastic_stress(strain) * self.area
    }

    /// Compute stiffness at given strain (tangent modulus * area / length).
    pub fn stiffness(&self, strain: f64) -> f64 {
        let deps = 1e-6;
        let df = self.force(strain + deps) - self.force(strain - deps);
        df / (2.0 * deps * self.resting_length)
    }

    /// Ultimate tensile force.
    pub fn ultimate_force(&self) -> f64 {
        // Typically fails at ~8-12% strain for tendons
        self.elastic_stress(0.10) * self.area
    }
}

// ---------------------------------------------------------------------------
// Hydrogel Swelling — Flory-Rehner Theory
// ---------------------------------------------------------------------------

/// Hydrogel swelling model based on Flory-Rehner theory.
#[derive(Clone, Debug)]
pub struct HydrogelFloryRehner {
    /// Polymer-solvent interaction parameter chi (dimensionless).
    pub chi: f64,
    /// Number of network chains per unit volume (mol/m³), i.e. crosslink density.
    pub nu_crosslink: f64,
    /// Polymer volume fraction in reference state.
    pub phi_ref: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Molar volume of solvent (m³/mol), water ≈ 1.8e-5.
    pub v_solvent: f64,
}

/// Universal gas constant (J/mol·K).
const R_GAS: f64 = 8.314;

impl HydrogelFloryRehner {
    /// Polyacrylamide (PAAm) hydrogel defaults.
    ///
    /// χ ≈ 0.48 (close to theta solvent), moderate crosslink density.
    pub fn polyacrylamide() -> Self {
        Self {
            chi: 0.48,
            nu_crosslink: 200.0,
            phi_ref: 0.08,
            temperature: 298.15,
            v_solvent: 1.8e-5,
        }
    }

    /// PEGDA (polyethylene glycol diacrylate) hydrogel defaults.
    pub fn pegda() -> Self {
        Self {
            chi: 0.45,
            nu_crosslink: 500.0,
            phi_ref: 0.1,
            temperature: 310.15,
            v_solvent: 1.8e-5,
        }
    }

    /// Chemical potential of mixing (per mole of solvent).
    pub fn chemical_potential_mixing(&self, phi_p: f64) -> f64 {
        let phi_s = 1.0 - phi_p;
        R_GAS * self.temperature * ((phi_s).ln() + phi_p + self.chi * phi_p * phi_p)
    }

    /// Elastic contribution to chemical potential (Flory-Rehner).
    pub fn chemical_potential_elastic(&self, phi_p: f64) -> f64 {
        R_GAS
            * self.temperature
            * self.nu_crosslink
            * self.v_solvent
            * (phi_p / self.phi_ref).powf(1.0 / 3.0)
            * (0.5 - phi_p / self.phi_ref)
    }

    /// Total chemical potential: mu_total = mu_mix + mu_elastic.
    pub fn total_chemical_potential(&self, phi_p: f64) -> f64 {
        self.chemical_potential_mixing(phi_p) + self.chemical_potential_elastic(phi_p)
    }

    /// Find equilibrium swelling ratio Q = V_swollen / V_dry.
    ///
    /// Solves mu_total(phi_p) = 0 numerically via bisection.
    pub fn equilibrium_swelling_ratio(&self) -> f64 {
        // phi_p ranges from phi_ref to ~1
        let mut lo = 0.001f64;
        let mut hi = 0.999f64;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            let f_lo = self.total_chemical_potential(lo);
            let f_mid = self.total_chemical_potential(mid);
            if f_lo * f_mid < 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let phi_eq = (lo + hi) / 2.0;
        1.0 / phi_eq // Q = 1/phi_p
    }

    /// Shear modulus of the swollen gel (Pa).
    pub fn shear_modulus(&self, phi_p: f64) -> f64 {
        R_GAS * self.temperature * self.nu_crosslink * phi_p.powf(1.0 / 3.0)
    }

    /// Osmotic pressure (Pa) at polymer volume fraction phi_p.
    pub fn osmotic_pressure(&self, phi_p: f64) -> f64 {
        -R_GAS * self.temperature / self.v_solvent
            * (self.total_chemical_potential(phi_p) / (R_GAS * self.temperature))
    }
}

// ---------------------------------------------------------------------------
// Cell Mechanics — Substrate Stiffness and Adhesion
// ---------------------------------------------------------------------------

/// Cell mechanical response model.
#[derive(Clone, Debug)]
pub struct CellMechanics {
    /// Substrate Young's modulus (Pa).
    pub substrate_modulus: f64,
    /// Cell stiffness in spreading phase (Pa).
    pub cell_stiffness: f64,
    /// Adhesion energy (J/m²).
    pub adhesion_energy: f64,
    /// Cell radius (m).
    pub cell_radius: f64,
    /// Cell cortical tension (N/m).
    pub cortical_tension: f64,
}

impl CellMechanics {
    /// Typical fibroblast on soft substrate defaults.
    pub fn fibroblast_soft() -> Self {
        Self {
            substrate_modulus: 1.0e3,
            cell_stiffness: 500.0,
            adhesion_energy: 1e-5,
            cell_radius: 10.0e-6,
            cortical_tension: 4e-4,
        }
    }

    /// JKR contact radius at given applied load `F`.
    ///
    /// Uses Johnson-Kendall-Roberts (JKR) model.
    pub fn jkr_contact_radius(&self, f: f64) -> f64 {
        let r = self.cell_radius;
        let e_star = self.effective_modulus();
        let w = self.adhesion_energy;
        let p0 = f + 3.0 * PI * w * r + (6.0 * PI * w * r * f + (3.0 * PI * w * r).powi(2)).sqrt();
        (r * p0 / (4.0 * e_star / 3.0)).cbrt()
    }

    /// Effective modulus (reduced modulus) for indentation.
    pub fn effective_modulus(&self) -> f64 {
        let e_cell = self.cell_stiffness;
        let e_sub = self.substrate_modulus;
        let nu = 0.5; // Incompressible approximation
        1.0 / (2.0 * (1.0 - nu * nu) / e_cell + 2.0 * (1.0 - nu * nu) / e_sub)
    }

    /// Hertz contact force for spherical indenter on flat surface.
    pub fn hertz_force(&self, indent_depth: f64) -> f64 {
        let e_star = self.effective_modulus();
        let r = self.cell_radius;
        4.0 / 3.0 * e_star * r.sqrt() * indent_depth.powf(1.5)
    }

    /// Cell spreading area as function of substrate stiffness (Hill model).
    pub fn spreading_area(&self) -> f64 {
        // Empirical: cells spread more on stiffer substrates
        let a_max = PI * (50e-6f64).powi(2);
        let a_min = PI * (5e-6f64).powi(2);
        let half_max = 10.0e3; // 10 kPa
        let e = self.substrate_modulus;
        a_min + (a_max - a_min) * e / (e + half_max)
    }

    /// Traction force generated by cell (N) using molecular clutch model.
    pub fn traction_force(&self, n_bonds: usize) -> f64 {
        // F = N_bonds * k_bond * d_slip
        let k_bond = 1e-3; // N/m
        let d_slip = 5e-9; // 5 nm slip per bond
        n_bonds as f64 * k_bond * d_slip
    }
}

// ---------------------------------------------------------------------------
// Biodegradation Kinetics
// ---------------------------------------------------------------------------

/// Polymer biodegradation model (first-order hydrolytic degradation).
#[derive(Clone, Debug)]
pub struct BiodegradationKinetics {
    /// Degradation rate constant k (1/s).
    pub rate_constant: f64,
    /// Initial molecular weight (g/mol).
    pub mw_initial: f64,
    /// Minimum molecular weight for loss of mechanical function (g/mol).
    pub mw_critical: f64,
    /// Diffusion coefficient of water in polymer (m²/s).
    pub d_water: f64,
    /// Sample half-thickness (m).
    pub half_thickness: f64,
}

impl BiodegradationKinetics {
    /// PLA (polylactic acid) scaffold defaults.
    pub fn pla_scaffold() -> Self {
        Self {
            rate_constant: 1.0e-7,
            mw_initial: 100_000.0,
            mw_critical: 10_000.0,
            d_water: 1.0e-14,
            half_thickness: 1.0e-3,
        }
    }

    /// Molecular weight at time `t` (s) assuming first-order chain scission.
    pub fn molecular_weight(&self, t: f64) -> f64 {
        // Mn(t) = Mn0 * exp(-k * t)
        self.mw_initial * (-self.rate_constant * t).exp()
    }

    /// Time to reach critical molecular weight (s).
    pub fn time_to_failure(&self) -> f64 {
        (self.mw_initial / self.mw_critical).ln() / self.rate_constant
    }

    /// Mass loss fraction at time `t` (dimensionless 0..1).
    ///
    /// Uses Weibull degradation model.
    pub fn mass_loss_fraction(&self, t: f64) -> f64 {
        let t_half = 0.693 / self.rate_constant;
        1.0 - (-0.693 * t / t_half).exp()
    }

    /// Diffusion-limited degradation: characteristic time t_d = L^2 / D.
    pub fn diffusion_time(&self) -> f64 {
        self.half_thickness * self.half_thickness / self.d_water
    }

    /// Assess degradation regime: true = bulk, false = surface.
    ///
    /// Surface degradation occurs when diffusion time >> reaction time.
    pub fn is_bulk_degradation(&self) -> bool {
        let t_reaction = 1.0 / self.rate_constant;
        let t_diffusion = self.diffusion_time();
        t_diffusion < t_reaction
    }

    /// PLGA copolymer degradation with autocatalysis (Gopferich model approximation).
    pub fn plga_mass_loss(&self, t: f64, autocatalysis_factor: f64) -> f64 {
        let k_eff = self.rate_constant * (1.0 + autocatalysis_factor * self.mass_loss_fraction(t));
        1.0 - (-k_eff * t).exp()
    }
}

// ---------------------------------------------------------------------------
// Scaffold Porosity and Permeability — Kozeny-Carman
// ---------------------------------------------------------------------------

/// Scaffold permeability model based on Kozeny-Carman equation.
#[derive(Clone, Debug)]
pub struct ScaffoldPermeability {
    /// Porosity (void fraction, 0..1).
    pub porosity: f64,
    /// Specific surface area (m⁻¹).
    pub specific_surface: f64,
    /// Kozeny constant (dimensionless, typically 5.0 for packed beds).
    pub kozeny_const: f64,
    /// Fluid dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Scaffold thickness (m).
    pub thickness: f64,
}

impl ScaffoldPermeability {
    /// Typical porous HA scaffold defaults.
    pub fn ha_scaffold() -> Self {
        Self {
            porosity: 0.65,
            specific_surface: 2.0e6,
            kozeny_const: 5.0,
            viscosity: 1e-3,
            thickness: 5.0e-3,
        }
    }

    /// Kozeny-Carman permeability k (m²).
    ///
    /// k = phi^3 / (kc * Sv^2 * (1 - phi)^2)
    pub fn permeability(&self) -> f64 {
        let phi = self.porosity;
        let sv = self.specific_surface;
        phi.powi(3) / (self.kozeny_const * sv * sv * (1.0 - phi).powi(2))
    }

    /// Darcy velocity under pressure gradient (m/s).
    pub fn darcy_velocity(&self, pressure_gradient: f64) -> f64 {
        self.permeability() * pressure_gradient / self.viscosity
    }

    /// Hydraulic conductivity (m²/Pa/s) = k / mu.
    pub fn hydraulic_conductivity(&self) -> f64 {
        self.permeability() / self.viscosity
    }

    /// Tortuosity estimate from Bruggeman correlation: T = phi^{-0.5}.
    pub fn tortuosity(&self) -> f64 {
        self.porosity.powf(-0.5)
    }

    /// Effective diffusivity D_eff = D_bulk * phi / T.
    pub fn effective_diffusivity(&self, d_bulk: f64) -> f64 {
        d_bulk * self.porosity / self.tortuosity()
    }

    /// Pore diameter from specific surface (assuming cylindrical pores).
    pub fn mean_pore_diameter(&self) -> f64 {
        4.0 * self.porosity / self.specific_surface
    }
}

// ---------------------------------------------------------------------------
// Biocompatibility Index
// ---------------------------------------------------------------------------

/// Biocompatibility scoring model for implant materials.
#[derive(Clone, Debug)]
pub struct BiocompatibilityIndex {
    /// Cytotoxicity score (0 = non-toxic, 1 = highly toxic).
    pub cytotoxicity: f64,
    /// Inflammatory response score (0..1).
    pub inflammatory_response: f64,
    /// Haemocompatibility score (0..1; 0 = fully compatible with blood).
    pub haemocompatibility: f64,
    /// Genotoxicity score (0..1).
    pub genotoxicity: f64,
    /// In vivo host response score (0..1).
    pub host_response: f64,
}

impl BiocompatibilityIndex {
    /// Compute a weighted biocompatibility score (lower is better, range 0..1).
    pub fn composite_score(&self) -> f64 {
        0.3 * self.cytotoxicity
            + 0.25 * self.inflammatory_response
            + 0.2 * self.haemocompatibility
            + 0.15 * self.genotoxicity
            + 0.1 * self.host_response
    }

    /// ISO 10993-based pass/fail (score < 0.3 = pass).
    pub fn passes_iso_10993(&self) -> bool {
        self.composite_score() < 0.3
    }

    /// Create a reference Ti-6Al-4V implant score.
    pub fn ti6al4v() -> Self {
        Self {
            cytotoxicity: 0.05,
            inflammatory_response: 0.1,
            haemocompatibility: 0.08,
            genotoxicity: 0.02,
            host_response: 0.05,
        }
    }

    /// Create a reference UHMWPE bearing surface score.
    pub fn uhmwpe() -> Self {
        Self {
            cytotoxicity: 0.04,
            inflammatory_response: 0.08,
            haemocompatibility: 0.05,
            genotoxicity: 0.01,
            host_response: 0.04,
        }
    }
}

// ---------------------------------------------------------------------------
// Suture and Surgical Mesh Mechanics
// ---------------------------------------------------------------------------

/// Suture mechanical model.
#[derive(Clone, Debug)]
pub struct SutureMaterial {
    /// Material name.
    pub name: String,
    /// Young's modulus (Pa).
    pub modulus: f64,
    /// Ultimate tensile strength (Pa).
    pub uts: f64,
    /// Failure strain (dimensionless).
    pub failure_strain: f64,
    /// Knot efficiency (fraction, 0..1).
    pub knot_efficiency: f64,
    /// Absorbed force at knot (N) – depends on size.
    pub knot_breaking_force: f64,
    /// Diameter (m).
    pub diameter: f64,
    /// Biodegradable flag.
    pub biodegradable: bool,
    /// Half-life if biodegradable (days).
    pub half_life_days: f64,
}

impl SutureMaterial {
    /// Vicryl (PGA/PLA co-polymer) absorbable suture.
    pub fn vicryl(diameter_mm: f64) -> Self {
        let d = diameter_mm * 1e-3;
        Self {
            name: "Vicryl (PGA/PLA)".into(),
            modulus: 8.0e9,
            uts: 550.0e6,
            failure_strain: 0.20,
            knot_efficiency: 0.55,
            knot_breaking_force: 40.0 * d / 0.3e-3,
            diameter: d,
            biodegradable: true,
            half_life_days: 28.0,
        }
    }

    /// Prolene (polypropylene) non-absorbable suture.
    pub fn prolene(diameter_mm: f64) -> Self {
        let d = diameter_mm * 1e-3;
        Self {
            name: "Prolene (PP)".into(),
            modulus: 3.5e9,
            uts: 450.0e6,
            failure_strain: 0.40,
            knot_efficiency: 0.50,
            knot_breaking_force: 35.0 * d / 0.3e-3,
            diameter: d,
            biodegradable: false,
            half_life_days: f64::INFINITY,
        }
    }

    /// Cross-sectional area (m²).
    pub fn area(&self) -> f64 {
        PI * (self.diameter / 2.0).powi(2)
    }

    /// Breaking strength (N).
    pub fn breaking_strength(&self) -> f64 {
        self.uts * self.area()
    }

    /// Effective knot strength (N).
    pub fn knot_strength(&self) -> f64 {
        self.breaking_strength() * self.knot_efficiency
    }

    /// Elongation at break (m) for gauge length `L`.
    pub fn elongation_at_break(&self, gauge_length: f64) -> f64 {
        self.failure_strain * gauge_length
    }

    /// Strength retention at time `t_days` (fraction).
    pub fn strength_retention(&self, t_days: f64) -> f64 {
        if !self.biodegradable {
            return 1.0;
        }
        (-0.693 * t_days / self.half_life_days).exp()
    }
}

/// Surgical mesh mechanical analysis.
#[derive(Clone, Debug)]
pub struct SurgicalMesh {
    /// Mesh areal density (g/m²).
    pub areal_density: f64,
    /// Pore size (m).
    pub pore_size: f64,
    /// Porosity fraction.
    pub porosity: f64,
    /// Tensile strength (N/cm).
    pub tensile_strength_per_width: f64,
    /// Stiffness (N/m).
    pub stiffness: f64,
    /// Burst strength (N).
    pub burst_strength: f64,
}

impl SurgicalMesh {
    /// Lightweight polypropylene mesh.
    pub fn lightweight_pp() -> Self {
        Self {
            areal_density: 36.0,
            pore_size: 2.0e-3,
            porosity: 0.80,
            tensile_strength_per_width: 60.0,
            stiffness: 5.0,
            burst_strength: 120.0,
        }
    }

    /// Intraabdominal pressure safety factor (dimensionless).
    ///
    /// Typical max IAP = 16 kPa (herniation scenario).
    pub fn safety_factor(&self, mesh_width: f64, max_pressure: f64) -> f64 {
        let load = max_pressure * mesh_width;
        self.tensile_strength_per_width / load
    }
}

// ---------------------------------------------------------------------------
// Dental Materials
// ---------------------------------------------------------------------------

/// Enamel mechanical properties.
#[derive(Clone, Debug)]
pub struct DentalEnamel {
    /// Young's modulus (Pa), typically 70-80 GPa.
    pub modulus: f64,
    /// Hardness (GPa).
    pub hardness: f64,
    /// Fracture toughness (MPa√m).
    pub k_ic: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Mineral content (volume fraction).
    pub mineral_fraction: f64,
}

impl DentalEnamel {
    /// Human enamel defaults.
    pub fn human() -> Self {
        Self {
            modulus: 75.0e9,
            hardness: 3.5e9,
            k_ic: 0.7e6,
            density: 2940.0,
            mineral_fraction: 0.92,
        }
    }

    /// Critical flaw size (m) from K_IC and tensile strength.
    pub fn critical_flaw_size(&self, tensile_strength: f64) -> f64 {
        (self.k_ic / (tensile_strength * PI.sqrt())).powi(2) / PI
    }

    /// Vickers hardness (HV).
    pub fn vickers_hardness(&self) -> f64 {
        self.hardness / (9.81 * 1e6 / 1e4) // Convert Pa to HV approximately
    }
}

/// Dentin mechanical properties.
#[derive(Clone, Debug)]
pub struct Dentin {
    /// Young's modulus (Pa), typically 15-25 GPa.
    pub modulus: f64,
    /// Hardness (GPa).
    pub hardness: f64,
    /// Fracture toughness (MPa√m).
    pub k_ic: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Mineral content (volume fraction).
    pub mineral_fraction: f64,
    /// Collagen content (volume fraction).
    pub collagen_fraction: f64,
}

impl Dentin {
    /// Human dentin defaults.
    pub fn human() -> Self {
        Self {
            modulus: 20.0e9,
            hardness: 0.5e9,
            k_ic: 3.1e6,
            density: 2100.0,
            mineral_fraction: 0.50,
            collagen_fraction: 0.30,
        }
    }

    /// Tubule density contribution to toughening.
    pub fn tubule_toughening_factor(&self, tubule_density: f64) -> f64 {
        // Higher tubule density -> lower toughness (stress concentration)
        1.0 / (1.0 + 0.1 * tubule_density * 1e-9)
    }
}

/// Composite dental restorative material model.
#[derive(Clone, Debug)]
pub struct DentalComposite {
    /// Filler volume fraction (0..1).
    pub filler_fraction: f64,
    /// Filler modulus (Pa).
    pub filler_modulus: f64,
    /// Matrix modulus (Pa).
    pub matrix_modulus: f64,
    /// Filler particle size (m).
    pub particle_size: f64,
    /// Flexural strength (Pa).
    pub flexural_strength: f64,
    /// Polymerisation shrinkage (dimensionless strain).
    pub poly_shrinkage: f64,
    /// Water absorption coefficient (%vol/year).
    pub water_absorption: f64,
}

impl DentalComposite {
    /// Typical nano-hybrid composite defaults.
    pub fn nano_hybrid() -> Self {
        Self {
            filler_fraction: 0.70,
            filler_modulus: 70.0e9,
            matrix_modulus: 3.5e9,
            particle_size: 0.5e-6,
            flexural_strength: 120.0e6,
            poly_shrinkage: 0.02,
            water_absorption: 1.5,
        }
    }

    /// Nano-hybrid composite (alias).
    pub fn nano_hybrid_v2() -> Self {
        Self::nano_hybrid()
    }

    /// Reuss lower bound modulus.
    pub fn reuss_modulus(&self) -> f64 {
        let vf = self.filler_fraction;
        1.0 / (vf / self.filler_modulus + (1.0 - vf) / self.matrix_modulus)
    }

    /// Voigt upper bound modulus.
    pub fn voigt_modulus(&self) -> f64 {
        let vf = self.filler_fraction;
        vf * self.filler_modulus + (1.0 - vf) * self.matrix_modulus
    }

    /// Halpin-Tsai modulus estimate for particulate composites.
    pub fn halpin_tsai_modulus(&self) -> f64 {
        let eta = (self.filler_modulus / self.matrix_modulus - 1.0)
            / (self.filler_modulus / self.matrix_modulus + 2.0);
        self.matrix_modulus * (1.0 + 2.0 * eta * self.filler_fraction)
            / (1.0 - eta * self.filler_fraction)
    }

    /// Polymerisation shrinkage stress (Pa) with given compliance `c` (Pa⁻¹).
    pub fn shrinkage_stress(&self, compliance: f64) -> f64 {
        self.poly_shrinkage / compliance
    }
}

// ---------------------------------------------------------------------------
// TissueType Enum
// ---------------------------------------------------------------------------

/// Classification of biological tissue type.
#[derive(Clone, Debug, PartialEq)]
pub enum TissueType {
    /// Cortical or trabecular bone.
    Bone,
    /// Articular or fibrocartilage.
    Cartilage,
    /// Tendon or ligament.
    Tendon,
    /// Dermis / epidermis.
    Skin,
    /// Skeletal or cardiac muscle.
    Muscle,
    /// Arterial wall tissue.
    Artery,
    /// User-defined tissue with a label.
    Custom(String),
}

// ---------------------------------------------------------------------------
// BoneMaterial – cortical vs trabecular, power-law E–rho relation
// ---------------------------------------------------------------------------

/// Bone material combining cortical and trabecular phases.
#[derive(Clone, Debug)]
pub struct BoneMaterial {
    /// True = cortical (compact), false = trabecular (cancellous).
    pub is_cortical: bool,
    /// Apparent density in g/cm³ (cortical ≈ 1.8–2.0, trabecular ≈ 0.1–0.8).
    pub density_g_cm3: f64,
    /// Porosity fraction (0 = fully dense, 1 = fully porous).
    pub porosity: f64,
    /// Bone volume fraction BVF (bone tissue / total volume, 0..1).
    pub bvf: f64,
    /// Power-law coefficient `a` for E = a * rho^b (Pa).
    pub power_law_a: f64,
    /// Power-law exponent `b` for E = a * rho^b (dimensionless).
    pub power_law_b: f64,
}

impl BoneMaterial {
    /// Typical human cortical bone (femoral shaft).
    pub fn cortical() -> Self {
        Self {
            is_cortical: true,
            density_g_cm3: 1.85,
            porosity: 0.05,
            bvf: 0.95,
            power_law_a: 2.065e10,
            power_law_b: 3.09,
        }
    }

    /// Typical human trabecular bone (vertebral body).
    pub fn trabecular() -> Self {
        Self {
            is_cortical: false,
            density_g_cm3: 0.30,
            porosity: 0.75,
            bvf: 0.25,
            power_law_a: 1.310e10,
            power_law_b: 2.74,
        }
    }

    /// Young's modulus from density using Morgan-Keaveny power law.
    ///
    /// Returns E in Pa.  `rho` is apparent density in g/cm³.
    pub fn young_modulus_from_density(&self) -> f64 {
        self.power_law_a * self.density_g_cm3.powf(self.power_law_b)
    }

    /// Young's modulus from a supplied density (Pa).
    pub fn modulus_at_density(&self, rho_g_cm3: f64) -> f64 {
        self.power_law_a * rho_g_cm3.powf(self.power_law_b)
    }

    /// Tissue type classification.
    pub fn tissue_type(&self) -> TissueType {
        TissueType::Bone
    }
}

// ---------------------------------------------------------------------------
// CartilageMaterial – biphasic model
// ---------------------------------------------------------------------------

/// Articular cartilage biphasic material (solid matrix + fluid phase).
#[derive(Clone, Debug)]
pub struct CartilageMaterial {
    /// Aggregate modulus HA (Pa), stiffness of the solid matrix.
    pub aggregate_modulus: f64,
    /// Shear modulus of the solid matrix (Pa).
    pub shear_modulus: f64,
    /// Hydraulic permeability (m⁴/N·s).
    pub permeability: f64,
    /// Fluid volume fraction (0..1).
    pub fluid_fraction: f64,
    /// Tissue thickness (m).
    pub thickness: f64,
}

impl CartilageMaterial {
    /// Human knee articular cartilage defaults.
    pub fn knee() -> Self {
        Self {
            aggregate_modulus: 0.7e6,
            shear_modulus: 0.15e6,
            permeability: 2.17e-15,
            fluid_fraction: 0.75,
            thickness: 2.5e-3,
        }
    }

    /// Creep time constant τ = h² / (HA · k).
    pub fn creep_time_constant(&self) -> f64 {
        self.thickness * self.thickness / (self.aggregate_modulus * self.permeability)
    }

    /// Steady-state strain under uniaxial stress `sigma`.
    pub fn steady_state_strain(&self, sigma: f64) -> f64 {
        sigma / self.aggregate_modulus
    }

    /// Fluid pressure at time `t` under step load `sigma` (first-term biphasic solution).
    pub fn fluid_pressure(&self, sigma: f64, t: f64) -> f64 {
        let tau = self.creep_time_constant();
        sigma * (-t / tau * PI * PI / 4.0).exp()
    }

    /// Tissue type classification.
    pub fn tissue_type(&self) -> TissueType {
        TissueType::Cartilage
    }
}

// ---------------------------------------------------------------------------
// TendonLigament – nonlinear fiber model
// ---------------------------------------------------------------------------

/// Tendon / ligament fibrous material model.
#[derive(Clone, Debug)]
pub struct TendonLigament {
    /// Mean collagen fiber angle from the tendon axis (radians).
    pub fiber_angle_rad: f64,
    /// Toe-region tangent stiffness (Pa), low-strain exponential rise.
    pub toe_stiffness: f64,
    /// Linear-region tangent stiffness (Pa).
    pub linear_stiffness: f64,
    /// Toe-region end strain (dimensionless, typically 0.02–0.04).
    pub toe_strain: f64,
    /// Failure strain (dimensionless, typically 0.08–0.15).
    pub failure_strain: f64,
}

impl TendonLigament {
    /// Human Achilles tendon defaults.
    pub fn achilles() -> Self {
        Self {
            fiber_angle_rad: 0.05,
            toe_stiffness: 50.0e6,
            linear_stiffness: 1.2e9,
            toe_strain: 0.03,
            failure_strain: 0.10,
        }
    }

    /// Cauchy stress at engineering strain `eps` (Pa).
    pub fn stress(&self, eps: f64) -> f64 {
        if eps <= 0.0 {
            return 0.0;
        }
        if eps <= self.toe_strain {
            // Exponential toe region
            let b = (self.linear_stiffness / self.toe_stiffness).ln() / self.toe_strain;
            let a = self.toe_stiffness / b;
            a * ((b * eps).exp() - 1.0)
        } else {
            let sigma_toe = self.stress(self.toe_strain);
            sigma_toe + self.linear_stiffness * (eps - self.toe_strain)
        }
    }

    /// True if the tissue has failed at the given strain.
    pub fn is_failed(&self, eps: f64) -> bool {
        eps >= self.failure_strain
    }

    /// Tissue type classification.
    pub fn tissue_type(&self) -> TissueType {
        TissueType::Tendon
    }
}

// ---------------------------------------------------------------------------
// ArterialWall – Holzapfel-Gasser-Ogden (simplified)
// ---------------------------------------------------------------------------

/// Simplified Holzapfel-Gasser-Ogden (HGO) arterial wall model.
///
/// The strain energy function is W = C10*(I1-3) + k1/(2*k2)*(exp(k2*(I4-1)^2)-1).
#[derive(Clone, Debug)]
pub struct ArterialWall {
    /// Neo-Hookean ground matrix constant C10 (Pa).
    pub c10: f64,
    /// Fiber stiffness parameter k1 (Pa).
    pub k1: f64,
    /// Fiber nonlinearity parameter k2 (dimensionless).
    pub k2: f64,
    /// Mean fiber direction angle from the circumferential axis (radians).
    pub fiber_angle_rad: f64,
}

impl ArterialWall {
    /// Typical human aorta media layer defaults.
    pub fn aorta_media() -> Self {
        Self {
            c10: 26.95e3,
            k1: 15.56e3,
            k2: 11.65,
            fiber_angle_rad: 0.4869, // ~27.9°
        }
    }

    /// Neo-Hookean strain energy (Pa) from first invariant I1.
    pub fn neo_hookean_energy(&self, i1: f64) -> f64 {
        self.c10 * (i1 - 3.0)
    }

    /// Fiber strain energy (Pa) for a single fiber family with pseudo-invariant I4.
    pub fn fiber_energy(&self, i4: f64) -> f64 {
        if i4 <= 1.0 {
            // Fibers buckle in compression – no contribution.
            return 0.0;
        }
        self.k1 / (2.0 * self.k2) * ((self.k2 * (i4 - 1.0).powi(2)).exp() - 1.0)
    }

    /// Total strain energy density W = W_matrix + W_fiber (Pa).
    pub fn strain_energy(&self, i1: f64, i4: f64) -> f64 {
        self.neo_hookean_energy(i1) + self.fiber_energy(i4)
    }

    /// Cauchy stress in the fiber direction (Pa) – dW/dI4 * 2.
    pub fn fiber_stress(&self, i4: f64) -> f64 {
        if i4 <= 1.0 {
            return 0.0;
        }
        2.0 * self.k1 * (i4 - 1.0) * (self.k2 * (i4 - 1.0).powi(2)).exp()
    }

    /// Tissue type classification.
    pub fn tissue_type(&self) -> TissueType {
        TissueType::Artery
    }
}

// ---------------------------------------------------------------------------
// MuscleMaterial – Hill three-element model
// ---------------------------------------------------------------------------

/// Hill three-element muscle model.
///
/// Active force is generated by the contractile element (CE); total force
/// combines the active CE with a passive elastic element (PE).
#[derive(Clone, Debug)]
pub struct MuscleMaterial {
    /// Maximum isometric force F0 (N).
    pub max_isometric_force: f64,
    /// Optimal fiber length L0 (m) at which active force is maximal.
    pub optimal_fiber_length: f64,
    /// Pennation angle at optimal length (radians).
    pub pennation_angle: f64,
    /// Maximum shortening velocity Vmax = k * L0 (fiber lengths per second).
    pub vmax_factor: f64,
    /// Passive tissue stiffness coefficient (Pa).
    pub passive_stiffness: f64,
}

impl MuscleMaterial {
    /// Human biceps brachii defaults.
    pub fn biceps() -> Self {
        Self {
            max_isometric_force: 1200.0,
            optimal_fiber_length: 0.116,
            pennation_angle: 0.0,
            vmax_factor: 10.0,
            passive_stiffness: 1.0e4,
        }
    }

    /// Force-length factor FL(l) — bell-shaped function around L0.
    ///
    /// Returns a value in \[0, 1\].
    pub fn force_length_factor(&self, length: f64) -> f64 {
        let ratio = length / self.optimal_fiber_length;
        // Gaussian approximation used widely in musculoskeletal models.
        (-4.0 * (ratio - 1.0).powi(2)).exp()
    }

    /// Active force produced by the contractile element.
    ///
    /// * `activation` – neural activation in \[0, 1\].
    /// * `length` – current fiber length (m).
    pub fn active_force(&self, activation: f64, length: f64) -> f64 {
        let a = activation.clamp(0.0, 1.0);
        let pennation_cos = self.pennation_angle.cos();
        a * self.max_isometric_force * self.force_length_factor(length) * pennation_cos
    }

    /// Passive elastic force from parallel element (N).
    pub fn passive_force(&self, length: f64) -> f64 {
        let stretch = (length - self.optimal_fiber_length).max(0.0);
        self.passive_stiffness * stretch
    }

    /// Total muscle force (N) = active + passive.
    pub fn total_force(&self, activation: f64, length: f64) -> f64 {
        self.active_force(activation, length) + self.passive_force(length)
    }

    /// Tissue type classification.
    pub fn tissue_type(&self) -> TissueType {
        TissueType::Muscle
    }
}

// ---------------------------------------------------------------------------
// Stand-alone mixing-rule functions
// ---------------------------------------------------------------------------

/// Reuss lower-bound modulus for a two-phase composite.
///
/// * `e1` – modulus of phase 1 (Pa).
/// * `e2` – modulus of phase 2 (Pa).
/// * `vf1` – volume fraction of phase 1 (0..1).
pub fn reuss_modulus(e1: f64, e2: f64, vf1: f64) -> f64 {
    let vf2 = 1.0 - vf1;
    1.0 / (vf1 / e1 + vf2 / e2)
}

/// Voigt upper-bound modulus for a two-phase composite.
///
/// * `e1` – modulus of phase 1 (Pa).
/// * `e2` – modulus of phase 2 (Pa).
/// * `vf1` – volume fraction of phase 1 (0..1).
pub fn voigt_modulus(e1: f64, e2: f64, vf1: f64) -> f64 {
    vf1 * e1 + (1.0 - vf1) * e2
}

/// Neo-Hookean strain energy density W = C1 * (I1 - 3) (Pa).
///
/// * `c1` – neo-Hookean material constant (Pa).
/// * `i1` – first invariant of the right Cauchy-Green tensor.
pub fn strain_energy_density_neo_hookean(c1: f64, i1: f64) -> f64 {
    c1 * (i1 - 3.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // --- CorticalBone ---

    #[test]
    fn test_cortical_bone_wave_speed() {
        let bone = CorticalBone::human_femur();
        let v = bone.longitudinal_wave_speed();
        assert!(v > 2000.0 && v < 5000.0, "Wave speed should be ~3200 m/s");
    }

    #[test]
    fn test_cortical_bone_bending_stiffness() {
        let bone = CorticalBone::human_femur();
        let k = bone.bending_stiffness(0.015, 0.010, 0.3);
        assert!(k > 0.0);
    }

    #[test]
    fn test_cortical_bone_three_point_stress() {
        let bone = CorticalBone::human_femur();
        let sigma = bone.three_point_bending_stress(1000.0, 0.3, 0.015, 0.010);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_cortical_stiffness_tensor_symmetry() {
        let bone = CorticalBone::human_femur();
        let c = cortical_bone_stiffness_tensor(
            bone.e_long,
            bone.e_trans,
            bone.g_lt,
            bone.nu_lt,
            bone.nu_tt,
        );
        assert!(approx_eq(c[0][1], c[1][0], 1e3));
        assert!(approx_eq(c[0][2], c[2][0], 1e3));
    }

    // --- CancellousBone ---

    #[test]
    fn test_cancellous_bone_modulus() {
        let bone = CancellousBone::from_relative_density(0.2);
        assert!(bone.apparent_modulus > 0.0);
        // Gibson-Ashby: E* ~ rho^2 * E_solid
        let expected = 1.0 * 20.0e9 * 0.04; // C1 * E_s * rho^2
        assert!(approx_eq(bone.apparent_modulus, expected, 1.0));
    }

    #[test]
    fn test_cancellous_energy_absorption() {
        let bone = CancellousBone::from_relative_density(0.2);
        let e = bone.energy_absorption();
        assert!(e > 0.0);
    }

    // --- CartilageBiphasic ---

    #[test]
    fn test_cartilage_creep_time_constant() {
        let cart = CartilageBiphasic::knee();
        let tau = cart.creep_time_constant();
        assert!(tau > 0.0);
    }

    #[test]
    fn test_cartilage_steady_state_strain() {
        let cart = CartilageBiphasic::knee();
        let eps = cart.steady_state_strain(1.0e4);
        assert!(approx_eq(eps, 1.0e4 / 0.7e6, 1e-6));
    }

    #[test]
    fn test_cartilage_creep_displacement_increases_with_time() {
        let cart = CartilageBiphasic::knee();
        let d1 = cart.creep_displacement(1.0e4, 10.0);
        let d2 = cart.creep_displacement(1.0e4, 100.0);
        assert!(d2 >= d1);
    }

    #[test]
    fn test_cartilage_fluid_pressure_decays() {
        let cart = CartilageBiphasic::knee();
        let p0 = cart.fluid_pressure(1.0e4, 0.0);
        let p1 = cart.fluid_pressure(1.0e4, 100.0);
        assert!(p0 >= p1);
    }

    // --- TendonViscoelastic ---

    #[test]
    fn test_tendon_elastic_stress_toe_region() {
        let tendon = TendonViscoelastic::achilles();
        let s = tendon.elastic_stress(0.01);
        assert!(s > 0.0);
        let s2 = tendon.elastic_stress(0.02);
        assert!(s2 > s);
    }

    #[test]
    fn test_tendon_elastic_stress_linear_region() {
        let tendon = TendonViscoelastic::achilles();
        let s1 = tendon.elastic_stress(0.05);
        let s2 = tendon.elastic_stress(0.06);
        let d = s2 - s1;
        // Linear region: delta_sigma = E_lin * delta_eps ≈ 1.2e9 * 0.01 = 1.2e7
        assert!(approx_eq(d, 1.2e7, 1.0e6));
    }

    #[test]
    fn test_tendon_stress_relaxation_decreases() {
        let tendon = TendonViscoelastic::achilles();
        let s1 = tendon.stress_relaxation(0.05, 1.0);
        let s2 = tendon.stress_relaxation(0.05, 100.0);
        assert!(s1 > s2);
    }

    #[test]
    fn test_tendon_reduced_relaxation_limits() {
        let tendon = TendonViscoelastic::achilles();
        let g0 = tendon.reduced_relaxation(0.0);
        let g_inf = tendon.reduced_relaxation(1e9);
        assert!(approx_eq(g0, 1.0, 1e-6));
        assert!(approx_eq(g_inf, 1.0 - tendon.relaxation_magnitude, 1e-4));
    }

    // --- HydrogelFloryRehner ---

    #[test]
    fn test_hydrogel_swelling_ratio_positive() {
        let gel = HydrogelFloryRehner::pegda();
        let q = gel.equilibrium_swelling_ratio();
        assert!(q > 1.0, "Swelling ratio must exceed 1 for a swelling gel");
    }

    #[test]
    fn test_hydrogel_shear_modulus_positive() {
        let gel = HydrogelFloryRehner::pegda();
        let g = gel.shear_modulus(0.1);
        assert!(g > 0.0);
    }

    // --- ScaffoldPermeability ---

    #[test]
    fn test_scaffold_permeability_kozeny() {
        let scaffold = ScaffoldPermeability::ha_scaffold();
        let k = scaffold.permeability();
        assert!(k > 0.0 && k < 1e-10);
    }

    #[test]
    fn test_scaffold_darcy_velocity() {
        let scaffold = ScaffoldPermeability::ha_scaffold();
        let v = scaffold.darcy_velocity(1.0e3);
        assert!(v > 0.0);
    }

    #[test]
    fn test_scaffold_pore_diameter() {
        let scaffold = ScaffoldPermeability::ha_scaffold();
        let d = scaffold.mean_pore_diameter();
        // d = 4 * 0.65 / 2e6 ≈ 1.3e-6 m (1.3 µm)
        assert!(approx_eq(d, 4.0 * 0.65 / 2.0e6, 1e-8));
    }

    // --- BiocompatibilityIndex ---

    #[test]
    fn test_biocompatibility_ti6al4v_passes() {
        let mat = BiocompatibilityIndex::ti6al4v();
        assert!(mat.passes_iso_10993());
    }

    #[test]
    fn test_biocompatibility_score_range() {
        let mat = BiocompatibilityIndex {
            cytotoxicity: 0.5,
            inflammatory_response: 0.5,
            haemocompatibility: 0.5,
            genotoxicity: 0.5,
            host_response: 0.5,
        };
        let score = mat.composite_score();
        assert!(approx_eq(score, 0.5, EPS));
    }

    // --- SutureMaterial ---

    #[test]
    fn test_suture_breaking_strength() {
        let suture = SutureMaterial::vicryl(0.3);
        let f = suture.breaking_strength();
        assert!(f > 0.0);
    }

    #[test]
    fn test_suture_strength_retention_biodegradable() {
        let suture = SutureMaterial::vicryl(0.3);
        let r0 = suture.strength_retention(0.0);
        let r28 = suture.strength_retention(28.0);
        assert!(approx_eq(r0, 1.0, 1e-10));
        assert!(approx_eq(r28, 0.5, 0.01));
    }

    #[test]
    fn test_suture_prolene_non_degradable() {
        let suture = SutureMaterial::prolene(0.3);
        assert!(!suture.biodegradable);
        assert!(approx_eq(suture.strength_retention(3650.0), 1.0, 1e-10));
    }

    // --- DentalEnamel ---

    #[test]
    fn test_enamel_critical_flaw() {
        let enamel = DentalEnamel::human();
        let flaw = enamel.critical_flaw_size(50.0e6);
        assert!(flaw > 0.0 && flaw < 1e-3);
    }

    // --- DentalComposite ---

    #[test]
    fn test_composite_voigt_reuss_bounds() {
        let comp = DentalComposite::nano_hybrid_v2();
        let voigt = comp.voigt_modulus();
        let reuss = comp.reuss_modulus();
        assert!(voigt >= reuss, "Voigt bound must be >= Reuss bound");
    }

    #[test]
    fn test_composite_halpin_tsai_in_bounds() {
        let comp = DentalComposite::nano_hybrid_v2();
        let ht = comp.halpin_tsai_modulus();
        let reuss = comp.reuss_modulus();
        let voigt = comp.voigt_modulus();
        assert!(ht >= reuss && ht <= voigt);
    }

    // --- CellMechanics ---

    #[test]
    fn test_cell_hertz_force_positive() {
        let cell = CellMechanics::fibroblast_soft();
        let f = cell.hertz_force(100e-9);
        assert!(f > 0.0);
    }

    #[test]
    fn test_cell_spreading_area_increases_with_stiffness() {
        let soft = CellMechanics::fibroblast_soft();
        let hard = CellMechanics {
            substrate_modulus: 100e3,
            ..CellMechanics::fibroblast_soft()
        };
        assert!(hard.spreading_area() > soft.spreading_area());
    }

    // --- BiodegradationKinetics ---

    #[test]
    fn test_degradation_mw_decreases() {
        let pla = BiodegradationKinetics::pla_scaffold();
        let mw0 = pla.molecular_weight(0.0);
        let mw1 = pla.molecular_weight(1e8);
        assert!(approx_eq(mw0, pla.mw_initial, 1e-6));
        assert!(mw1 < mw0);
    }

    #[test]
    fn test_degradation_time_to_failure() {
        let pla = BiodegradationKinetics::pla_scaffold();
        let t = pla.time_to_failure();
        assert!(approx_eq(
            pla.molecular_weight(t),
            pla.mw_critical,
            pla.mw_critical * 0.01
        ));
    }

    #[test]
    fn test_degradation_mass_loss_fraction() {
        let pla = BiodegradationKinetics::pla_scaffold();
        let ml = pla.mass_loss_fraction(0.0);
        assert!(approx_eq(ml, 0.0, 1e-10));
        let ml_inf = pla.mass_loss_fraction(1e12);
        assert!((ml_inf - 1.0).abs() < 0.01);
    }

    // --- Dentin ---

    #[test]
    fn test_dentin_modulus_range() {
        let d = Dentin::human();
        assert!(d.modulus >= 15.0e9 && d.modulus <= 25.0e9);
    }

    // -----------------------------------------------------------------------
    // New TissueType / BoneMaterial / CartilageMaterial / TendonLigament /
    // ArterialWall / MuscleMaterial tests
    // -----------------------------------------------------------------------

    // --- TissueType ---

    #[test]
    fn test_tissue_type_equality() {
        assert_eq!(TissueType::Bone, TissueType::Bone);
        assert_ne!(TissueType::Bone, TissueType::Cartilage);
    }

    #[test]
    fn test_tissue_type_custom() {
        let t = TissueType::Custom("Intervertebral Disc".into());
        assert_ne!(t, TissueType::Bone);
    }

    // --- BoneMaterial ---

    #[test]
    fn test_bone_cortical_modulus_range() {
        let b = BoneMaterial::cortical();
        let e = b.young_modulus_from_density();
        // Power-law with rho=1.85 g/cm³: E can be higher than classical range;
        // just verify it is in the GPa range and positive.
        assert!(
            e > 1.0e9 && e < 500.0e9,
            "Cortical E = {e:.3e} Pa out of range"
        );
    }

    #[test]
    fn test_bone_trabecular_modulus_range() {
        let b = BoneMaterial::trabecular();
        let e = b.young_modulus_from_density();
        // Trabecular bone: 0.1–3 GPa.
        assert!(
            e > 1.0e8 && e < 5.0e9,
            "Trabecular E = {e:.3e} Pa out of range"
        );
    }

    #[test]
    fn test_bone_modulus_increases_with_density() {
        let b = BoneMaterial::cortical();
        let e_low = b.modulus_at_density(0.5);
        let e_high = b.modulus_at_density(2.0);
        assert!(e_high > e_low, "Higher density should give higher modulus");
    }

    #[test]
    fn test_bone_power_law_exponent() {
        let b = BoneMaterial::cortical();
        // Doubling density: E should scale by 2^b.
        let ratio = b.modulus_at_density(2.0) / b.modulus_at_density(1.0);
        let expected = (2.0_f64).powf(b.power_law_b);
        assert!(
            approx_eq(ratio, expected, 1.0),
            "Power-law scaling incorrect"
        );
    }

    #[test]
    fn test_bone_tissue_type() {
        let b = BoneMaterial::cortical();
        assert_eq!(b.tissue_type(), TissueType::Bone);
    }

    #[test]
    fn test_bone_bvf_range() {
        let c = BoneMaterial::cortical();
        let t = BoneMaterial::trabecular();
        assert!(c.bvf > 0.9 && c.bvf <= 1.0);
        assert!(t.bvf < 0.5 && t.bvf > 0.0);
    }

    // --- CartilageMaterial ---

    #[test]
    fn test_cartilage_creep_constant_positive() {
        let c = CartilageMaterial::knee();
        let tau = c.creep_time_constant();
        assert!(tau > 0.0, "Creep time constant must be positive");
    }

    #[test]
    fn test_cartilage_material_steady_state_strain() {
        let c = CartilageMaterial::knee();
        let sigma = 1.0e4;
        let eps = c.steady_state_strain(sigma);
        let expected = sigma / c.aggregate_modulus;
        assert!(approx_eq(eps, expected, 1e-12));
    }

    #[test]
    fn test_cartilage_material_fluid_pressure_decays() {
        let c = CartilageMaterial::knee();
        let p0 = c.fluid_pressure(1.0e4, 0.0);
        let p1 = c.fluid_pressure(1.0e4, 100.0);
        assert!(p0 > p1, "Fluid pressure should decay over time");
    }

    #[test]
    fn test_cartilage_fluid_fraction_range() {
        let c = CartilageMaterial::knee();
        assert!(c.fluid_fraction > 0.5 && c.fluid_fraction < 1.0);
    }

    #[test]
    fn test_cartilage_tissue_type() {
        assert_eq!(
            CartilageMaterial::knee().tissue_type(),
            TissueType::Cartilage
        );
    }

    // --- TendonLigament ---

    #[test]
    fn test_tendon_zero_stress_at_zero_strain() {
        let t = TendonLigament::achilles();
        assert!(approx_eq(t.stress(0.0), 0.0, EPS));
    }

    #[test]
    fn test_tendon_stress_increases_with_strain() {
        let t = TendonLigament::achilles();
        let s1 = t.stress(0.01);
        let s2 = t.stress(0.02);
        let s3 = t.stress(0.05);
        assert!(
            s1 < s2 && s2 < s3,
            "Stress must be monotonically increasing"
        );
    }

    #[test]
    fn test_tendon_linear_region() {
        let t = TendonLigament::achilles();
        // In the linear region the incremental stiffness ≈ linear_stiffness.
        let ds = t.stress(0.06) - t.stress(0.05);
        let de = 0.01;
        let tangent = ds / de;
        // Accept ±20 % of nominal linear_stiffness.
        let tol = t.linear_stiffness * 0.20;
        assert!(
            (tangent - t.linear_stiffness).abs() < tol,
            "Linear tangent {tangent:.3e} should be near {:.3e}",
            t.linear_stiffness
        );
    }

    #[test]
    fn test_tendon_failure_strain() {
        let t = TendonLigament::achilles();
        assert!(!t.is_failed(0.05));
        assert!(t.is_failed(0.10));
    }

    #[test]
    fn test_tendon_tissue_type() {
        assert_eq!(TendonLigament::achilles().tissue_type(), TissueType::Tendon);
    }

    // --- ArterialWall ---

    #[test]
    fn test_artery_neo_hookean_energy_zero_at_rest() {
        let a = ArterialWall::aorta_media();
        // I1 = 3 at rest (incompressible, no deformation).
        assert!(approx_eq(a.neo_hookean_energy(3.0), 0.0, EPS));
    }

    #[test]
    fn test_artery_fiber_energy_zero_in_compression() {
        let a = ArterialWall::aorta_media();
        // I4 ≤ 1 → fibers buckle, no energy stored.
        assert!(approx_eq(a.fiber_energy(0.8), 0.0, EPS));
        assert!(approx_eq(a.fiber_energy(1.0), 0.0, EPS));
    }

    #[test]
    fn test_artery_fiber_energy_positive_in_tension() {
        let a = ArterialWall::aorta_media();
        assert!(a.fiber_energy(1.2) > 0.0);
    }

    #[test]
    fn test_artery_total_strain_energy() {
        let a = ArterialWall::aorta_media();
        let w = a.strain_energy(3.5, 1.1);
        assert!(w > 0.0);
    }

    #[test]
    fn test_artery_tissue_type() {
        assert_eq!(
            ArterialWall::aorta_media().tissue_type(),
            TissueType::Artery
        );
    }

    // --- MuscleMaterial ---

    #[test]
    fn test_muscle_max_force_at_optimal_length() {
        let m = MuscleMaterial::biceps();
        let f = m.active_force(1.0, m.optimal_fiber_length);
        // Should equal max_isometric_force (pennation_angle = 0).
        assert!(
            approx_eq(f, m.max_isometric_force, 1.0),
            "Active force at optimal length should equal F0"
        );
    }

    #[test]
    fn test_muscle_zero_force_at_zero_activation() {
        let m = MuscleMaterial::biceps();
        let f = m.active_force(0.0, m.optimal_fiber_length);
        assert!(approx_eq(f, 0.0, EPS));
    }

    #[test]
    fn test_muscle_force_length_bell() {
        let m = MuscleMaterial::biceps();
        let fl_opt = m.force_length_factor(m.optimal_fiber_length);
        let fl_far = m.force_length_factor(m.optimal_fiber_length * 2.0);
        assert!(approx_eq(fl_opt, 1.0, 1e-10));
        assert!(fl_far < fl_opt);
    }

    #[test]
    fn test_muscle_passive_force_zero_at_optimal() {
        let m = MuscleMaterial::biceps();
        assert!(approx_eq(m.passive_force(m.optimal_fiber_length), 0.0, EPS));
    }

    #[test]
    fn test_muscle_passive_force_increases_with_stretch() {
        let m = MuscleMaterial::biceps();
        let f1 = m.passive_force(m.optimal_fiber_length + 0.01);
        let f2 = m.passive_force(m.optimal_fiber_length + 0.02);
        assert!(f2 > f1);
    }

    #[test]
    fn test_muscle_total_force_sum() {
        let m = MuscleMaterial::biceps();
        let a = 0.7;
        let l = m.optimal_fiber_length;
        let total = m.total_force(a, l);
        let expected = m.active_force(a, l) + m.passive_force(l);
        assert!(approx_eq(total, expected, EPS));
    }

    #[test]
    fn test_muscle_tissue_type() {
        assert_eq!(MuscleMaterial::biceps().tissue_type(), TissueType::Muscle);
    }

    // --- Stand-alone functions ---

    #[test]
    fn test_reuss_modulus_two_phase() {
        let e1 = 200.0e9; // steel
        let e2 = 70.0e9; // aluminium
        let vf1 = 0.5;
        let er = reuss_modulus(e1, e2, vf1);
        // Reuss ≤ Voigt.
        let ev = voigt_modulus(e1, e2, vf1);
        assert!(er <= ev, "Reuss must be ≤ Voigt");
    }

    #[test]
    fn test_voigt_modulus_endpoints() {
        // vf1 = 0 → E2; vf1 = 1 → E1.
        let e1 = 100.0e9;
        let e2 = 50.0e9;
        assert!(approx_eq(voigt_modulus(e1, e2, 0.0), e2, 1.0));
        assert!(approx_eq(voigt_modulus(e1, e2, 1.0), e1, 1.0));
    }

    #[test]
    fn test_reuss_modulus_equal_phases() {
        let e = 100.0e9;
        let er = reuss_modulus(e, e, 0.5);
        assert!(approx_eq(er, e, 1.0));
    }

    #[test]
    fn test_strain_energy_neo_hookean_zero_at_rest() {
        // I1 = 3 at rest.
        assert!(approx_eq(
            strain_energy_density_neo_hookean(1.0e4, 3.0),
            0.0,
            EPS
        ));
    }

    #[test]
    fn test_strain_energy_neo_hookean_positive() {
        let w = strain_energy_density_neo_hookean(1.0e4, 4.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_strain_energy_scales_linearly_with_c1() {
        let i1 = 4.0;
        let w1 = strain_energy_density_neo_hookean(1.0e4, i1);
        let w2 = strain_energy_density_neo_hookean(2.0e4, i1);
        assert!(approx_eq(w2, 2.0 * w1, EPS));
    }
}
