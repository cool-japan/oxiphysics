// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Nanoparticle MD: surface effects, quantum size effects, colloidal interactions.
//!
//! This module provides models for:
//! - Size-dependent melting point depression (Gibbs-Thomson)
//! - DLVO colloidal interaction theory
//! - Ligand shell steric repulsion (Alexander-de Gennes)
//! - Self-assembly (DNA-mediated, depletion attraction)
//! - Quantum dot optical properties (Brus equation)
//! - Nanofluid thermal conductivity (Maxwell model)
//! - Brownian dynamics (Langevin, Stokes-Einstein)
//! - Nanoparticle aggregation kinetics (DLCA/RLCA)

#![allow(dead_code)]

use std::f64::consts::PI;

// Physical constants
/// Boltzmann constant \[J/K\]
pub const K_B: f64 = 1.380_649e-23;
/// Elementary charge \[C\]
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Electron mass \[kg\]
pub const M_ELECTRON: f64 = 9.109_383_701_5e-31;
/// Planck constant \[J·s\]
pub const H_PLANCK: f64 = 6.626_070_15e-34;
/// Vacuum permittivity \[F/m\]
pub const EPS_0: f64 = 8.854_187_812_8e-12;
/// Avogadro number \[1/mol\]
pub const N_A: f64 = 6.022_140_76e23;

// ─── NanoparticleProperties ───────────────────────────────────────────────────

/// Material properties for nanoparticle calculations.
#[derive(Debug, Clone)]
pub struct NanoMaterial {
    /// Bulk melting point \[K\]
    pub t_bulk_melt: f64,
    /// Solid-liquid surface energy \[J/m²\]
    pub surface_energy: f64,
    /// Density \[kg/m³\]
    pub density: f64,
    /// Latent heat of fusion \[J/kg\]
    pub latent_heat: f64,
    /// Name of material
    pub name: String,
}

impl NanoMaterial {
    /// Create a new NanoMaterial.
    pub fn new(
        name: &str,
        t_bulk_melt: f64,
        surface_energy: f64,
        density: f64,
        latent_heat: f64,
    ) -> Self {
        Self {
            t_bulk_melt,
            surface_energy,
            density,
            latent_heat,
            name: name.to_string(),
        }
    }

    /// Gold nanoparticle preset (T_m = 1337 K).
    pub fn gold() -> Self {
        Self::new("Au", 1337.0, 1.37, 19_320.0, 63_700.0)
    }

    /// Silver nanoparticle preset (T_m = 1235 K).
    pub fn silver() -> Self {
        Self::new("Ag", 1235.0, 1.17, 10_490.0, 104_600.0)
    }
}

/// Size-dependent properties of a spherical nanoparticle.
#[derive(Debug, Clone)]
pub struct NanoparticleProperties {
    /// Particle diameter \[m\]
    pub diameter: f64,
    /// Material properties
    pub material: NanoMaterial,
}

impl NanoparticleProperties {
    /// Create a new NanoparticleProperties.
    pub fn new(diameter: f64, material: NanoMaterial) -> Self {
        Self { diameter, material }
    }

    /// Radius \[m\].
    pub fn radius(&self) -> f64 {
        self.diameter / 2.0
    }

    /// Surface-to-volume ratio for a sphere: 6/d \[1/m\].
    pub fn surface_to_volume_ratio(&self) -> f64 {
        6.0 / self.diameter
    }

    /// Surface area of the sphere \[m²\].
    pub fn surface_area(&self) -> f64 {
        PI * self.diameter * self.diameter
    }

    /// Volume of the sphere \[m³\].
    pub fn volume(&self) -> f64 {
        PI * self.diameter.powi(3) / 6.0
    }

    /// Melting point depression via Gibbs-Thomson (Pawlow) equation:
    /// T_m(d) = T_bulk * (1 - 4σ / (ρ·L·d))
    pub fn melting_point(&self) -> f64 {
        let m = &self.material;
        let correction = 4.0 * m.surface_energy / (m.density * m.latent_heat * self.diameter);
        m.t_bulk_melt * (1.0 - correction)
    }

    /// Depression relative to bulk melting point \[K\].
    pub fn melting_point_depression(&self) -> f64 {
        self.material.t_bulk_melt - self.melting_point()
    }
}

// ─── DlvoTheory ───────────────────────────────────────────────────────────────

/// DLVO interaction between two identical colloidal spheres.
///
/// Combines van der Waals attraction with electrostatic double-layer repulsion.
#[derive(Debug, Clone)]
pub struct DlvoTheory {
    /// Particle radius \[m\]
    pub radius: f64,
    /// Hamaker constant \[J\]
    pub hamaker: f64,
    /// Dielectric constant of medium (relative)
    pub epsilon_r: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Valence of electrolyte (z:z salt)
    pub valence: f64,
    /// Zeta potential \[V\]
    pub zeta_potential: f64,
    /// Inverse Debye length κ \[1/m\]
    pub kappa: f64,
}

impl DlvoTheory {
    /// Create a new DlvoTheory instance.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius: f64,
        hamaker: f64,
        epsilon_r: f64,
        temperature: f64,
        valence: f64,
        zeta_potential: f64,
        kappa: f64,
    ) -> Self {
        Self {
            radius,
            hamaker,
            epsilon_r,
            temperature,
            valence,
            zeta_potential,
            kappa,
        }
    }

    /// Van der Waals attraction energy at surface-surface separation D \[J\]:
    /// W_vdW = -A * R / (12 * D)
    pub fn vdw_energy(&self, separation: f64) -> f64 {
        -self.hamaker * self.radius / (12.0 * separation)
    }

    /// Electrostatic double-layer repulsion energy \[J\]:
    /// W_EDL = 64π ε₀ εᵣ R (kT/ze)² tanh²(zeζ/4kT) exp(-κD)
    pub fn edl_energy(&self, separation: f64) -> f64 {
        let kt = K_B * self.temperature;
        let ze = self.valence * E_CHARGE;
        let eps = EPS_0 * self.epsilon_r;
        let thermal_voltage = kt / ze;
        let tanh_arg = (ze * self.zeta_potential) / (4.0 * kt);
        let tanh_sq = tanh_arg.tanh().powi(2);
        64.0 * PI
            * eps
            * self.radius
            * thermal_voltage.powi(2)
            * tanh_sq
            * (-self.kappa * separation).exp()
    }

    /// Total DLVO interaction energy \[J\].
    pub fn total_energy(&self, separation: f64) -> f64 {
        self.vdw_energy(separation) + self.edl_energy(separation)
    }

    /// Energy barrier: maximum of total_energy in range \[d_min, d_max\], sampled at n_pts.
    pub fn energy_barrier(&self, d_min: f64, d_max: f64, n_pts: usize) -> f64 {
        let mut max_e = f64::NEG_INFINITY;
        for i in 0..n_pts {
            let d = d_min + (d_max - d_min) * (i as f64) / ((n_pts - 1) as f64);
            let e = self.total_energy(d);
            if e > max_e {
                max_e = e;
            }
        }
        max_e
    }

    /// Colloidal stability ratio W: W > 1 means stable suspension.
    /// Approximated as energy_barrier / (k_B T).
    pub fn stability_ratio(&self) -> f64 {
        let barrier = self.energy_barrier(1e-10, 50e-9, 500);
        let kt = K_B * self.temperature;
        (barrier / kt).exp()
    }
}

// ─── LigandShell ──────────────────────────────────────────────────────────────

/// Surface ligand coating on a nanoparticle (Alexander-de Gennes brush model).
#[derive(Debug, Clone)]
pub struct LigandShell {
    /// Nanoparticle radius \[m\]
    pub nanoparticle_radius: f64,
    /// Ligand chain length (number of monomers)
    pub chain_length: usize,
    /// Monomer size (Kuhn length) \[m\]
    pub monomer_size: f64,
    /// Grafting density \[chains/m²\]
    pub grafting_density: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Solvent quality parameter χ (0=good, 0.5=theta)
    pub chi: f64,
}

impl LigandShell {
    /// Create a new LigandShell.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nanoparticle_radius: f64,
        chain_length: usize,
        monomer_size: f64,
        grafting_density: f64,
        temperature: f64,
        chi: f64,
    ) -> Self {
        Self {
            nanoparticle_radius,
            chain_length,
            monomer_size,
            grafting_density,
            temperature,
            chi,
        }
    }

    /// Brush height (Alexander-de Gennes scaling): L ≈ N·b·(σ·b²)^(1/3)
    pub fn brush_height(&self) -> f64 {
        let b = self.monomer_size;
        let n = self.chain_length as f64;
        let sigma = self.grafting_density;
        n * b * (sigma * b * b).powf(1.0 / 3.0)
    }

    /// Ligand density in chains/nm²  (convenience conversion).
    pub fn grafting_density_per_nm2(&self) -> f64 {
        self.grafting_density * 1e-18
    }

    /// Steric repulsion energy per unit area (Alexander-de Gennes) \[J/m²\].
    /// G(D) = k_B T / b³ * L^5/2 * D^-5/4  (for D < L)
    pub fn steric_repulsion_energy_per_area(&self, gap: f64) -> f64 {
        let l = self.brush_height();
        if gap >= l {
            return 0.0;
        }
        let kt = K_B * self.temperature;
        let b = self.monomer_size;
        kt / b.powi(3) * l.powf(2.5) * gap.powf(-1.25)
    }

    /// Total steric repulsion between two identical brushed particles at surface gap D \[J\].
    /// Approximated via Derjaguin: W ≈ 2π R · G(D)
    pub fn steric_repulsion_total(&self, gap: f64) -> f64 {
        let g = self.steric_repulsion_energy_per_area(gap);
        2.0 * PI * self.nanoparticle_radius * g
    }
}

// ─── SelfAssembly ─────────────────────────────────────────────────────────────

/// Nanoparticle self-assembly models: DNA-mediated and depletion attraction.
#[derive(Debug, Clone)]
pub struct SelfAssembly {
    /// Particle radius \[m\]
    pub particle_radius: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Depletant (polymer/small colloid) radius \[m\]
    pub depletant_radius: f64,
    /// Depletant volume fraction
    pub depletant_phi: f64,
}

impl SelfAssembly {
    /// Create a SelfAssembly instance.
    pub fn new(
        particle_radius: f64,
        temperature: f64,
        depletant_radius: f64,
        depletant_phi: f64,
    ) -> Self {
        Self {
            particle_radius,
            temperature,
            depletant_radius,
            depletant_phi,
        }
    }

    /// Depletion (Asakura-Oosawa) attraction energy at surface gap D \[J\]:
    /// W_dep = -ρ_dep · k_B T · V_overlap(D)
    ///
    /// V_overlap for two spheres of radius R with depletion zone radius r_d:
    /// V_overlap(D) = π/6 * (2R + 2r_d - D)² * (D + 4R + 2r_d) / 4   (for D < 2*r_d)
    pub fn depletion_energy(&self, separation: f64) -> f64 {
        let rd = self.depletant_radius;
        let r = self.particle_radius;
        let kt = K_B * self.temperature;
        // Number density of depletant
        let rho_dep = self.depletant_phi / (4.0 / 3.0 * PI * rd.powi(3));
        if separation >= 2.0 * rd {
            return 0.0;
        }
        // Overlap volume: exact formula for sphere-sphere depletion zones
        let h = 2.0 * rd - separation;
        let r_eff = r + rd;
        // V_overlap = pi * h^2 / 3 * (3 * r_eff - h) ... corrected Derjaguin
        let v_overlap = PI * h.powi(2) / 3.0 * (3.0 * r_eff - h);
        -rho_dep * kt * v_overlap
    }

    /// DNA hybridization free energy (simplified): ΔG = ΔH - T·ΔS \[J\].
    pub fn dna_hybridization_energy(delta_h: f64, delta_s: f64, temperature: f64) -> f64 {
        delta_h - temperature * delta_s
    }

    /// Number of DNA sticky ends needed for stable assembly (rough estimate):
    /// N_strands ≈ |ΔG_hyb| / (k_B T) * factor
    pub fn required_strands_for_stability(hybridization_energy: f64, temperature: f64) -> f64 {
        let kt = K_B * temperature;
        hybridization_energy.abs() / kt
    }
}

// ─── QuantumDotModel ──────────────────────────────────────────────────────────

/// Quantum dot optical properties via Brus equation.
#[derive(Debug, Clone)]
pub struct QuantumDotModel {
    /// Radius of the quantum dot \[m\]
    pub radius: f64,
    /// Bulk band gap \[J\]
    pub e_gap_bulk: f64,
    /// Effective electron mass (relative to m_e)
    pub m_e_eff: f64,
    /// Effective hole mass (relative to m_e)
    pub m_h_eff: f64,
    /// Dielectric constant of material (relative)
    pub epsilon_r: f64,
}

impl QuantumDotModel {
    /// Create a new QuantumDotModel.
    pub fn new(radius: f64, e_gap_bulk: f64, m_e_eff: f64, m_h_eff: f64, epsilon_r: f64) -> Self {
        Self {
            radius,
            e_gap_bulk,
            m_e_eff,
            m_h_eff,
            epsilon_r,
        }
    }

    /// CdSe quantum dot preset.
    pub fn cdse(radius: f64) -> Self {
        // CdSe: E_gap=1.74 eV, m_e*=0.13 m_e, m_h*=0.45 m_e, ε_r=10.6
        Self::new(radius, 1.74 * E_CHARGE, 0.13, 0.45, 10.6)
    }

    /// Brus equation: ΔE = ħ²π²/(2μ r²) - 1.8 e²/(4πεε₀ r)
    /// Returns size-dependent band gap \[J\].
    pub fn band_gap(&self) -> f64 {
        let h_bar = H_PLANCK / (2.0 * PI);
        let mu_inv = 1.0 / (self.m_e_eff * M_ELECTRON) + 1.0 / (self.m_h_eff * M_ELECTRON);
        let confinement = h_bar.powi(2) * PI.powi(2) / 2.0 * mu_inv / self.radius.powi(2);
        let coulomb = 1.8 * E_CHARGE.powi(2) / (4.0 * PI * EPS_0 * self.epsilon_r * self.radius);
        self.e_gap_bulk + confinement - coulomb
    }

    /// Band gap in electron-volts \[eV\].
    pub fn band_gap_ev(&self) -> f64 {
        self.band_gap() / E_CHARGE
    }

    /// Emission wavelength \[m\] corresponding to band gap.
    pub fn emission_wavelength(&self) -> f64 {
        let c = 2.997_924_58e8_f64;
        H_PLANCK * c / self.band_gap()
    }

    /// Band gap increases as radius decreases (quantum confinement).
    pub fn confinement_energy(&self) -> f64 {
        let h_bar = H_PLANCK / (2.0 * PI);
        let mu_inv = 1.0 / (self.m_e_eff * M_ELECTRON) + 1.0 / (self.m_h_eff * M_ELECTRON);
        h_bar.powi(2) * PI.powi(2) / 2.0 * mu_inv / self.radius.powi(2)
    }
}

// ─── NanofluidConductivity ────────────────────────────────────────────────────

/// Maxwell model for effective thermal conductivity of a nanofluid.
#[derive(Debug, Clone)]
pub struct NanofluidConductivity {
    /// Fluid thermal conductivity k_f \[W/(m·K)\]
    pub k_fluid: f64,
    /// Nanoparticle thermal conductivity k_p \[W/(m·K)\]
    pub k_particle: f64,
    /// Nanoparticle volume fraction φ
    pub phi: f64,
}

impl NanofluidConductivity {
    /// Create a new NanofluidConductivity.
    pub fn new(k_fluid: f64, k_particle: f64, phi: f64) -> Self {
        Self {
            k_fluid,
            k_particle,
            phi,
        }
    }

    /// Effective thermal conductivity via Maxwell model:
    /// k_eff/k_f = (k_p + 2k_f + 2φ(k_p - k_f)) / (k_p + 2k_f - φ(k_p - k_f))
    pub fn effective_conductivity(&self) -> f64 {
        let kp = self.k_particle;
        let kf = self.k_fluid;
        let phi = self.phi;
        let num = kp + 2.0 * kf + 2.0 * phi * (kp - kf);
        let den = kp + 2.0 * kf - phi * (kp - kf);
        kf * num / den
    }

    /// Enhancement relative to pure fluid.
    pub fn conductivity_enhancement(&self) -> f64 {
        self.effective_conductivity() / self.k_fluid - 1.0
    }

    /// Hamilton-Crosser model for non-spherical particles (shape factor n=3 for spheres).
    pub fn hamilton_crosser(&self, shape_factor: f64) -> f64 {
        let kp = self.k_particle;
        let kf = self.k_fluid;
        let phi = self.phi;
        let num = kp + (shape_factor - 1.0) * kf + (shape_factor - 1.0) * phi * (kp - kf);
        let den = kp + (shape_factor - 1.0) * kf - phi * (kp - kf);
        kf * num / den
    }
}

// ─── BrownianDynamics ─────────────────────────────────────────────────────────

/// Langevin / Brownian dynamics for nanoparticles in a viscous fluid.
#[derive(Debug, Clone)]
pub struct BrownianDynamics {
    /// Particle radius \[m\]
    pub radius: f64,
    /// Fluid dynamic viscosity η \[Pa·s\]
    pub viscosity: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Particle mass \[kg\]
    pub mass: f64,
}

impl BrownianDynamics {
    /// Create a new BrownianDynamics instance.
    pub fn new(radius: f64, viscosity: f64, temperature: f64, mass: f64) -> Self {
        Self {
            radius,
            viscosity,
            temperature,
            mass,
        }
    }

    /// Stokes drag coefficient γ = 6πηr \[kg/s\].
    pub fn drag_coefficient(&self) -> f64 {
        6.0 * PI * self.viscosity * self.radius
    }

    /// Translational diffusion coefficient (Stokes-Einstein): D = k_B T / (6πηr) \[m²/s\].
    pub fn diffusion_coefficient(&self) -> f64 {
        K_B * self.temperature / self.drag_coefficient()
    }

    /// Rotational diffusion coefficient: D_r = k_B T / (8πηr³) \[1/s\].
    pub fn rotational_diffusion(&self) -> f64 {
        K_B * self.temperature / (8.0 * PI * self.viscosity * self.radius.powi(3))
    }

    /// Mean-squared displacement for time t \[m²\]: `r²` = 6 D t (3D).
    pub fn mean_squared_displacement(&self, time: f64) -> f64 {
        6.0 * self.diffusion_coefficient() * time
    }

    /// Momentum relaxation time τ = m / γ \[s\].
    pub fn relaxation_time(&self) -> f64 {
        self.mass / self.drag_coefficient()
    }

    /// Perform one Langevin step: returns new position increment Δx \[m\] (1D, for test).
    /// Δx = F_ext/γ · dt + sqrt(2 D dt) · ξ    where ξ is a unit Gaussian
    pub fn step_1d(&self, dt: f64, external_force: f64, noise: f64) -> f64 {
        let d = self.diffusion_coefficient();
        let gamma = self.drag_coefficient();
        (external_force / gamma) * dt + (2.0 * d * dt).sqrt() * noise
    }

    /// Sedimentation length: l_s = k_B T / (m_eff · g) \[m\].
    pub fn sedimentation_length(&self, effective_mass: f64, g: f64) -> f64 {
        K_B * self.temperature / (effective_mass * g)
    }
}

// ─── NanoparticleAggregation ──────────────────────────────────────────────────

/// Aggregation kinetics and fractal dimension of nanoparticle clusters.
#[derive(Debug, Clone)]
pub struct NanoparticleAggregation {
    /// Primary particle radius \[m\]
    pub primary_radius: f64,
    /// Diffusion coefficient of primary particles \[m²/s\]
    pub diffusion: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Number concentration \[1/m³\]
    pub number_concentration: f64,
    /// Regime: true = DLCA (fast), false = RLCA (slow)
    pub is_dlca: bool,
}

impl NanoparticleAggregation {
    /// Create a new NanoparticleAggregation.
    pub fn new(
        primary_radius: f64,
        diffusion: f64,
        temperature: f64,
        number_concentration: f64,
        is_dlca: bool,
    ) -> Self {
        Self {
            primary_radius,
            diffusion,
            temperature,
            number_concentration,
            is_dlca,
        }
    }

    /// Fractal dimension: DLCA ≈ 1.75, RLCA ≈ 2.1
    pub fn fractal_dimension(&self) -> f64 {
        if self.is_dlca { 1.75 } else { 2.1 }
    }

    /// Smoluchowski rate constant for fast (DLCA) aggregation k_fast \[m³/s\]:
    /// k = 8 k_B T / (3 η) ... but here we use diffusion directly: k = 8π D r
    pub fn smoluchowski_rate_fast(&self) -> f64 {
        8.0 * PI * self.diffusion * self.primary_radius
    }

    /// Smoluchowski rate constant for slow (RLCA) aggregation, with stability ratio W.
    pub fn smoluchowski_rate_slow(&self, stability_ratio: f64) -> f64 {
        self.smoluchowski_rate_fast() / stability_ratio
    }

    /// Aggregation half-time \[s\]: t_1/2 = 1 / (k * n_0)
    pub fn aggregation_half_time(&self) -> f64 {
        1.0 / (self.smoluchowski_rate_fast() * self.number_concentration)
    }

    /// Cluster radius for aggregate of n primary particles \[m\]:
    /// R_g = r_0 * n^(1/d_f)
    pub fn cluster_radius(&self, n_particles: f64) -> f64 {
        let df = self.fractal_dimension();
        self.primary_radius * n_particles.powf(1.0 / df)
    }

    /// Debye-Hückel potential at distance r from a sphere of potential ζ \[V\]:
    /// φ(r) = ζ * (a/r) * exp(-κ(r-a))
    pub fn debye_huckel_potential(&self, r: f64, zeta: f64, kappa: f64) -> f64 {
        let a = self.primary_radius;
        if r < a {
            return zeta;
        }
        zeta * (a / r) * (-(kappa * (r - a))).exp()
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ── NanoparticleProperties ──

    #[test]
    fn test_surface_to_volume_ratio() {
        let mat = NanoMaterial::gold();
        let np = NanoparticleProperties::new(10e-9, mat);
        let sv = np.surface_to_volume_ratio();
        let expected = 6.0 / 10e-9;
        assert!(
            approx_eq(sv, expected, 1.0),
            "S/V ratio mismatch: {:.6} vs {:.6}",
            sv,
            expected
        );
    }

    #[test]
    fn test_melting_point_smaller_lower() {
        let mat1 = NanoMaterial::gold();
        let mat2 = NanoMaterial::gold();
        let np_small = NanoparticleProperties::new(2e-9, mat1);
        let np_large = NanoparticleProperties::new(20e-9, mat2);
        assert!(
            np_small.melting_point() < np_large.melting_point(),
            "smaller np should have lower melting point"
        );
    }

    #[test]
    fn test_melting_point_below_bulk() {
        let mat = NanoMaterial::gold();
        let bulk_tm = mat.t_bulk_melt;
        let np = NanoparticleProperties::new(5e-9, mat);
        assert!(
            np.melting_point() < bulk_tm,
            "melting point should be below bulk"
        );
    }

    #[test]
    fn test_melting_point_depression_positive() {
        let np = NanoparticleProperties::new(5e-9, NanoMaterial::gold());
        assert!(
            np.melting_point_depression() > 0.0,
            "depression should be positive"
        );
    }

    #[test]
    fn test_surface_area_formula() {
        let d = 10e-9_f64;
        let np = NanoparticleProperties::new(d, NanoMaterial::silver());
        let sa = np.surface_area();
        let expected = PI * d * d;
        assert!(
            approx_eq(sa, expected, 1e-30),
            "surface area: {:.6e} vs {:.6e}",
            sa,
            expected
        );
    }

    #[test]
    fn test_volume_formula() {
        let d = 10e-9_f64;
        let np = NanoparticleProperties::new(d, NanoMaterial::silver());
        let vol = np.volume();
        let expected = PI * d.powi(3) / 6.0;
        assert!(
            approx_eq(vol, expected, 1e-40),
            "volume: {:.6e} vs {:.6e}",
            vol,
            expected
        );
    }

    // ── DlvoTheory ──

    fn default_dlvo() -> DlvoTheory {
        DlvoTheory::new(
            50e-9,   // radius 50 nm
            1.0e-20, // Hamaker A = 10 zJ
            80.0,    // water ε_r
            298.0,   // T = 298 K
            1.0,     // monovalent salt
            -30e-3,  // zeta = -30 mV
            1e8,     // κ = 10 nm^{-1}
        )
    }

    #[test]
    fn test_vdw_energy_negative() {
        let dlvo = default_dlvo();
        let w = dlvo.vdw_energy(1e-9);
        assert!(
            w < 0.0,
            "vdW energy should be negative (attractive): {:.6e}",
            w
        );
    }

    #[test]
    fn test_edl_energy_positive() {
        let dlvo = default_dlvo();
        let w = dlvo.edl_energy(1e-9);
        assert!(
            w > 0.0,
            "EDL energy should be positive (repulsive): {:.6e}",
            w
        );
    }

    #[test]
    fn test_dlvo_energy_barrier_finite() {
        let dlvo = default_dlvo();
        let barrier = dlvo.energy_barrier(1e-10, 50e-9, 200);
        assert!(barrier.is_finite(), "energy barrier should be finite");
    }

    #[test]
    fn test_vdw_increases_at_contact() {
        let dlvo = default_dlvo();
        let w1 = dlvo.vdw_energy(1e-9).abs();
        let w2 = dlvo.vdw_energy(5e-9).abs();
        assert!(w1 > w2, "vdW magnitude should increase at shorter distance");
    }

    #[test]
    fn test_edl_decays_exponentially() {
        let dlvo = default_dlvo();
        let w1 = dlvo.edl_energy(1e-9);
        let w2 = dlvo.edl_energy(2e-9);
        assert!(w1 > w2, "EDL energy should decay with separation");
    }

    // ── LigandShell ──

    fn default_ligand() -> LigandShell {
        LigandShell::new(
            10e-9,  // 10 nm radius
            50,     // 50-mer
            0.4e-9, // 0.4 nm monomer
            1e18,   // grafting density 1 chain/nm²
            298.0, 0.0, // good solvent
        )
    }

    #[test]
    fn test_brush_height_positive() {
        let lig = default_ligand();
        let h = lig.brush_height();
        assert!(h > 0.0, "brush height should be positive: {:.6e}", h);
    }

    #[test]
    fn test_steric_repulsion_zero_beyond_brush() {
        let lig = default_ligand();
        let h = lig.brush_height();
        let w = lig.steric_repulsion_total(h * 1.1);
        assert!(
            approx_eq(w, 0.0, 1e-50),
            "steric repulsion should be 0 beyond brush: {:.6e}",
            w
        );
    }

    #[test]
    fn test_steric_repulsion_positive_inside_brush() {
        let lig = default_ligand();
        let h = lig.brush_height();
        let w = lig.steric_repulsion_total(h * 0.5);
        assert!(w > 0.0, "steric repulsion should be positive inside brush");
    }

    #[test]
    fn test_brush_height_scales_with_chain_length() {
        let lig_short = LigandShell::new(10e-9, 20, 0.4e-9, 1e18, 298.0, 0.0);
        let lig_long = LigandShell::new(10e-9, 100, 0.4e-9, 1e18, 298.0, 0.0);
        assert!(
            lig_long.brush_height() > lig_short.brush_height(),
            "longer chain should give larger brush"
        );
    }

    // ── SelfAssembly ──

    #[test]
    fn test_depletion_energy_zero_large_separation() {
        let sa = SelfAssembly::new(50e-9, 298.0, 5e-9, 0.1);
        let e = sa.depletion_energy(20e-9); // larger than 2 * depletant_radius
        assert!(
            approx_eq(e, 0.0, 1e-40),
            "depletion energy should be 0 for large separation"
        );
    }

    #[test]
    fn test_depletion_energy_negative_at_contact() {
        let sa = SelfAssembly::new(50e-9, 298.0, 5e-9, 0.1);
        let e = sa.depletion_energy(1e-9);
        assert!(
            e < 0.0,
            "depletion energy should be negative (attractive): {:.6e}",
            e
        );
    }

    #[test]
    fn test_dna_hybridization_energy() {
        // Typical DNA: ΔH ~ -40 kJ/mol, ΔS ~ -100 J/(mol·K), T=310K
        let dh = -40_000.0 / N_A;
        let ds = -100.0 / N_A;
        let dg = SelfAssembly::dna_hybridization_energy(dh, ds, 310.0);
        // ΔG should be negative for stable hybridization
        assert!(
            dg < 0.0,
            "DNA hybridization ΔG should be negative: {:.6e}",
            dg
        );
    }

    #[test]
    fn test_required_strands_positive() {
        let n = SelfAssembly::required_strands_for_stability(-1e-20, 298.0);
        assert!(n > 0.0, "required strands should be positive");
    }

    // ── QuantumDotModel ──

    #[test]
    fn test_band_gap_larger_than_bulk() {
        let qd = QuantumDotModel::cdse(2e-9);
        assert!(
            qd.band_gap() > qd.e_gap_bulk,
            "nanoparticle band gap should exceed bulk: {:.6} > {:.6}",
            qd.band_gap_ev(),
            qd.e_gap_bulk / E_CHARGE
        );
    }

    #[test]
    fn test_smaller_dot_larger_gap() {
        let qd_small = QuantumDotModel::cdse(1e-9);
        let qd_large = QuantumDotModel::cdse(5e-9);
        assert!(
            qd_small.band_gap() > qd_large.band_gap(),
            "smaller dot should have larger band gap (quantum confinement)"
        );
    }

    #[test]
    fn test_confinement_energy_positive() {
        let qd = QuantumDotModel::cdse(3e-9);
        assert!(
            qd.confinement_energy() > 0.0,
            "confinement energy should be positive"
        );
    }

    #[test]
    fn test_emission_wavelength_blue_shifts_smaller() {
        let qd_small = QuantumDotModel::cdse(1e-9);
        let qd_large = QuantumDotModel::cdse(5e-9);
        assert!(
            qd_small.emission_wavelength() < qd_large.emission_wavelength(),
            "smaller dot should emit at shorter wavelength (blue shift)"
        );
    }

    #[test]
    fn test_band_gap_ev_reasonable() {
        // CdSe 3 nm dot: expect between 1.7 and 3.5 eV
        let qd = QuantumDotModel::cdse(3e-9);
        let eg = qd.band_gap_ev();
        assert!(eg > 1.7 && eg < 4.0, "band gap out of range: {:.6} eV", eg);
    }

    // ── NanofluidConductivity ──

    #[test]
    fn test_maxwell_model_between_limits() {
        // k_eff should be between k_fluid and k_particle
        let nc = NanofluidConductivity::new(0.6, 400.0, 0.01); // water + Cu
        let k_eff = nc.effective_conductivity();
        assert!(
            k_eff > nc.k_fluid && k_eff < nc.k_particle,
            "k_eff={:.6} should be between k_f={:.6} and k_p={:.6}",
            k_eff,
            nc.k_fluid,
            nc.k_particle
        );
    }

    #[test]
    fn test_maxwell_model_zero_phi() {
        // φ=0 → k_eff = k_fluid
        let nc = NanofluidConductivity::new(0.6, 400.0, 0.0);
        let k_eff = nc.effective_conductivity();
        assert!(
            approx_eq(k_eff, 0.6, 1e-10),
            "at φ=0, k_eff should equal k_fluid: {:.6}",
            k_eff
        );
    }

    #[test]
    fn test_maxwell_enhancement_positive() {
        let nc = NanofluidConductivity::new(0.6, 400.0, 0.05);
        assert!(
            nc.conductivity_enhancement() > 0.0,
            "enhancement should be positive"
        );
    }

    #[test]
    fn test_hamilton_crosser_sphere_same_as_maxwell() {
        // For n=3 (sphere), HC reduces to Maxwell
        let nc = NanofluidConductivity::new(0.6, 400.0, 0.02);
        let maxwell = nc.effective_conductivity();
        let hc = nc.hamilton_crosser(3.0);
        assert!(
            approx_eq(maxwell, hc, 1e-12),
            "HC with n=3 should equal Maxwell: {:.6} vs {:.6}",
            maxwell,
            hc
        );
    }

    // ── BrownianDynamics ──

    #[test]
    fn test_stokes_einstein_d_proportional_inverse_r() {
        let bd1 = BrownianDynamics::new(5e-9, 1e-3, 298.0, 1e-21);
        let bd2 = BrownianDynamics::new(10e-9, 1e-3, 298.0, 1e-21);
        let ratio = bd1.diffusion_coefficient() / bd2.diffusion_coefficient();
        assert!(
            approx_eq(ratio, 2.0, 1e-10),
            "D ∝ 1/r: ratio should be 2.0, got {:.6}",
            ratio
        );
    }

    #[test]
    fn test_diffusion_positive() {
        let bd = BrownianDynamics::new(5e-9, 1e-3, 298.0, 1e-21);
        assert!(bd.diffusion_coefficient() > 0.0);
    }

    #[test]
    fn test_msd_linear_in_time() {
        let bd = BrownianDynamics::new(5e-9, 1e-3, 298.0, 1e-21);
        let msd1 = bd.mean_squared_displacement(1.0);
        let msd2 = bd.mean_squared_displacement(2.0);
        assert!(
            approx_eq(msd2 / msd1, 2.0, 1e-10),
            "MSD should be linear in time: ratio={:.6}",
            msd2 / msd1
        );
    }

    #[test]
    fn test_drag_coefficient_formula() {
        let r = 5e-9_f64;
        let eta = 1e-3_f64;
        let bd = BrownianDynamics::new(r, eta, 298.0, 1e-21);
        let gamma = bd.drag_coefficient();
        let expected = 6.0 * PI * eta * r;
        assert!(
            approx_eq(gamma, expected, 1e-30),
            "drag: {:.6e} vs {:.6e}",
            gamma,
            expected
        );
    }

    #[test]
    fn test_rotational_diffusion_positive() {
        let bd = BrownianDynamics::new(5e-9, 1e-3, 298.0, 1e-21);
        assert!(bd.rotational_diffusion() > 0.0);
    }

    // ── NanoparticleAggregation ──

    #[test]
    fn test_fractal_dimension_dlca() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let df = agg.fractal_dimension();
        assert!(
            approx_eq(df, 1.75, 1e-10),
            "DLCA fractal dimension should be 1.75, got {:.6}",
            df
        );
    }

    #[test]
    fn test_fractal_dimension_rlca() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, false);
        let df = agg.fractal_dimension();
        assert!(
            approx_eq(df, 2.1, 1e-10),
            "RLCA fractal dimension should be 2.1, got {:.6}",
            df
        );
    }

    #[test]
    fn test_cluster_radius_grows_with_n() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let r1 = agg.cluster_radius(10.0);
        let r2 = agg.cluster_radius(100.0);
        assert!(r2 > r1, "cluster radius should grow with n_particles");
    }

    #[test]
    fn test_debye_huckel_potential_decays() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let kappa = 1e8_f64;
        let zeta = -30e-3_f64;
        let phi_near = agg.debye_huckel_potential(6e-9, zeta, kappa);
        let phi_far = agg.debye_huckel_potential(20e-9, zeta, kappa);
        assert!(
            phi_near.abs() > phi_far.abs(),
            "Debye-Hückel potential should decay: near={:.6e}, far={:.6e}",
            phi_near,
            phi_far
        );
    }

    #[test]
    fn test_aggregation_half_time_positive() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let t = agg.aggregation_half_time();
        assert!(
            t > 0.0,
            "aggregation half time should be positive: {:.6e}",
            t
        );
    }

    #[test]
    fn test_smoluchowski_rate_fast_positive() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        assert!(agg.smoluchowski_rate_fast() > 0.0);
    }

    #[test]
    fn test_smoluchowski_rate_slow_less_than_fast() {
        let agg = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let k_fast = agg.smoluchowski_rate_fast();
        let k_slow = agg.smoluchowski_rate_slow(10.0);
        assert!(k_slow < k_fast, "slow rate should be less than fast rate");
    }

    #[test]
    fn test_dlca_fractal_less_than_rlca() {
        let dlca = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, true);
        let rlca = NanoparticleAggregation::new(5e-9, 4.4e-11, 298.0, 1e17, false);
        assert!(
            dlca.fractal_dimension() < rlca.fractal_dimension(),
            "DLCA df < RLCA df"
        );
    }

    #[test]
    fn test_sedimentation_length_positive() {
        let bd = BrownianDynamics::new(50e-9, 1e-3, 298.0, 5e-19);
        let g = 9.81_f64;
        let ls = bd.sedimentation_length(5e-19, g);
        assert!(
            ls > 0.0,
            "sedimentation length should be positive: {:.6e}",
            ls
        );
    }

    #[test]
    fn test_band_gap_increases_as_r_decreases() {
        // Sweep radii: band gap should increase as r decreases
        let radii = [5e-9_f64, 4e-9, 3e-9, 2e-9, 1e-9];
        let gaps: Vec<f64> = radii
            .iter()
            .map(|&r| QuantumDotModel::cdse(r).band_gap())
            .collect();
        for i in 1..gaps.len() {
            assert!(
                gaps[i] > gaps[i - 1],
                "band gap should increase as r decreases: gaps[{}]={:.6e} <= gaps[{}]={:.6e}",
                i,
                gaps[i],
                i - 1,
                gaps[i - 1]
            );
        }
    }
}
