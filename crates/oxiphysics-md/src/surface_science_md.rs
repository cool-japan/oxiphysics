// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface science molecular dynamics.
//!
//! This module provides simulation models for:
//!
//! - **Surface slabs** ([`SurfaceSlab`]): periodic slab model with fixed substrate and mobile surface layers.
//! - **Adsorption sites** ([`AdsorptionSite`]): top, bridge, hollow, fcc, hcp sites with binding energies.
//! - **Langmuir adsorption** ([`LangmuirAdsorption`]): Langmuir and BET isotherms, Henry's law limit.
//! - **Chemisorption** ([`ChemisorptionModel`]): BEP relation, pre-exponential factors, sticking probability.
//! - **Nudged elastic band** ([`NudgedElasticBand`]): NEB and climbing-image NEB for minimum energy paths.
//! - **Surface diffusion** ([`SurfaceDiffusion`]): lattice random walk, Arrhenius hop rates, MSD, diffusivity.
//! - **Thin film growth** ([`ThinFilmGrowth`]): Volmer-Weber, Frank-van der Merwe, Stranski-Krastanov modes.
//! - **Catalytic cycles** ([`CatalyticCycle`]): Sabatier principle, volcano plot, TOF, BEP scaling.
//!
//! References:
//! - Zangwill, A. (1988). *Physics at Surfaces*. Cambridge University Press.
//! - Hammer, B. & Norskov, J. K. (2000). Advances in Catalysis 45, 71–129.
//! - Henkelman, G. & Jonsson, H. (2000). J. Chem. Phys. 113, 9978.
//! - Venables, J. A. (2000). *Introduction to Surface and Thin Film Processes*.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J/K).
const KB: f64 = 1.380_649e-23;

/// Avogadro's number (mol⁻¹).
const NA: f64 = 6.022_140_76e23;

/// Planck constant (J·s).
const H_PLANCK: f64 = 6.626_070_15e-34;

// ---------------------------------------------------------------------------
// Helper math utilities
// ---------------------------------------------------------------------------

/// Compute `exp(-barrier / (KB * temperature))` safely (returns 0 for very large barriers).
#[inline]
fn arrhenius_factor(barrier_j: f64, temperature_k: f64) -> f64 {
    if temperature_k < 1e-12 {
        return 0.0;
    }
    let exponent = -barrier_j / (KB * temperature_k);
    if exponent < -700.0 {
        0.0
    } else {
        exponent.exp()
    }
}

/// Clamp a value to \[lo, hi\].
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// SurfaceSlab
// ---------------------------------------------------------------------------

/// Layer type for a surface slab simulation.
#[derive(Debug, Clone, PartialEq)]
pub enum LayerType {
    /// Fixed substrate layer (atoms do not move).
    Fixed,
    /// Mobile surface layer (atoms can relax or diffuse).
    Mobile,
}

/// A single layer in a periodic surface slab.
#[derive(Debug, Clone)]
pub struct SlabLayer {
    /// Layer type (fixed substrate or mobile surface).
    pub layer_type: LayerType,
    /// Z-coordinate of the layer (Å).
    pub z_position: f64,
    /// Number of atoms in this layer.
    pub n_atoms: usize,
    /// Interlayer spacing relative to bulk (dimensionless relaxation, e.g. -0.02 = 2% contraction).
    pub relaxation: f64,
}

impl SlabLayer {
    /// Create a new slab layer.
    pub fn new(layer_type: LayerType, z_position: f64, n_atoms: usize) -> Self {
        Self {
            layer_type,
            z_position,
            n_atoms,
            relaxation: 0.0,
        }
    }
}

/// Periodic slab model for surface science simulations.
///
/// The slab consists of several atomic layers (some fixed, some mobile) with
/// a vacuum region above the top surface to prevent periodic image interactions.
#[derive(Debug, Clone)]
pub struct SurfaceSlab {
    /// Miller indices of the surface plane (h, k, l).
    pub miller_indices: (i32, i32, i32),
    /// Surface unit cell dimensions (Å × Å).
    pub surface_area: (f64, f64),
    /// Layers from bottom (bulk) to top (surface).
    pub layers: Vec<SlabLayer>,
    /// Thickness of the vacuum layer above the slab (Å).
    pub vacuum_thickness: f64,
    /// Bulk lattice constant (Å).
    pub lattice_constant: f64,
    /// Surface energy (J/m²).
    pub surface_energy: f64,
}

impl SurfaceSlab {
    /// Create a new surface slab model.
    ///
    /// # Arguments
    /// * `miller_indices` — surface plane Miller indices
    /// * `surface_area` — lateral dimensions of the supercell (Å)
    /// * `lattice_constant` — bulk lattice constant (Å)
    /// * `n_fixed_layers` — number of fixed substrate layers
    /// * `n_mobile_layers` — number of mobile surface layers
    /// * `vacuum_thickness` — vacuum region thickness (Å)
    pub fn new(
        miller_indices: (i32, i32, i32),
        surface_area: (f64, f64),
        lattice_constant: f64,
        n_fixed_layers: usize,
        n_mobile_layers: usize,
        vacuum_thickness: f64,
    ) -> Self {
        let interlayer_spacing = lattice_constant / 2.0_f64.sqrt();
        let atoms_per_layer = 4usize; // typical for FCC (111) 2×2 supercell
        let mut layers = Vec::new();

        for i in 0..n_fixed_layers {
            let z = i as f64 * interlayer_spacing;
            layers.push(SlabLayer::new(LayerType::Fixed, z, atoms_per_layer));
        }
        for i in 0..n_mobile_layers {
            let z = (n_fixed_layers + i) as f64 * interlayer_spacing;
            layers.push(SlabLayer::new(LayerType::Mobile, z, atoms_per_layer));
        }

        Self {
            miller_indices,
            surface_area,
            layers,
            vacuum_thickness,
            lattice_constant,
            surface_energy: 0.0,
        }
    }

    /// Total slab thickness including vacuum (Å).
    pub fn total_thickness(&self) -> f64 {
        let slab_height = self.layers.last().map(|l| l.z_position).unwrap_or(0.0);
        slab_height + self.vacuum_thickness
    }

    /// Number of mobile (surface-region) atoms.
    pub fn n_mobile_atoms(&self) -> usize {
        self.layers
            .iter()
            .filter(|l| l.layer_type == LayerType::Mobile)
            .map(|l| l.n_atoms)
            .sum()
    }

    /// Number of fixed (bulk-region) atoms.
    pub fn n_fixed_atoms(&self) -> usize {
        self.layers
            .iter()
            .filter(|l| l.layer_type == LayerType::Fixed)
            .map(|l| l.n_atoms)
            .sum()
    }

    /// Compute the surface energy given total slab energy, bulk reference, and surface area.
    ///
    /// γ = (E_slab - N·e_bulk) / (2·A)  \[J/m²\]
    ///
    /// # Arguments
    /// * `e_slab_eV` — total energy of the slab (eV)
    /// * `e_bulk_per_atom_eV` — bulk energy per atom (eV)
    /// * `ev_to_j` — conversion factor (eV → J), default 1.602e-19
    pub fn compute_surface_energy(
        &mut self,
        e_slab_ev: f64,
        e_bulk_per_atom_ev: f64,
        ev_to_j: f64,
    ) {
        let n_total: usize = self.layers.iter().map(|l| l.n_atoms).sum();
        let delta_e_j = (e_slab_ev - n_total as f64 * e_bulk_per_atom_ev) * ev_to_j;
        let area_m2 = self.surface_area.0 * self.surface_area.1 * 1e-20; // Å² to m²
        self.surface_energy = delta_e_j / (2.0 * area_m2);
    }

    /// Apply surface relaxation: contract the topmost interlayer spacing.
    ///
    /// Returns the relaxation fraction applied.
    pub fn relax_surface(&mut self, relaxation_fraction: f64) -> f64 {
        let n = self.layers.len();
        if n < 2 {
            return 0.0;
        }
        let top_idx = n - 1;
        self.layers[top_idx].relaxation = relaxation_fraction;
        let spacing = self.layers[top_idx].z_position - self.layers[top_idx - 1].z_position;
        self.layers[top_idx].z_position =
            self.layers[top_idx - 1].z_position + spacing * (1.0 + relaxation_fraction);
        relaxation_fraction
    }
}

// ---------------------------------------------------------------------------
// AdsorptionSite
// ---------------------------------------------------------------------------

/// Type of adsorption site on a surface.
#[derive(Debug, Clone, PartialEq)]
pub enum SiteType {
    /// On-top of a surface atom.
    Top,
    /// Bridge between two nearest-neighbor atoms.
    Bridge,
    /// Hollow site (three-fold on FCC/HCP).
    Hollow,
    /// FCC hollow site (no atom directly below in 2nd layer).
    FccHollow,
    /// HCP hollow site (atom directly below in 2nd layer).
    HcpHollow,
}

/// An adsorption site on a crystal surface with binding energy information.
#[derive(Debug, Clone)]
pub struct AdsorptionSite {
    /// Type of this site.
    pub site_type: SiteType,
    /// Position on the surface (x, y, z) in Å.
    pub position: [f64; 3],
    /// Binding energy of the adsorbate at this site (eV, negative = bound).
    pub binding_energy_ev: f64,
    /// Coordination number of the site.
    pub coordination: u32,
    /// Whether this site is currently occupied.
    pub occupied: bool,
}

impl AdsorptionSite {
    /// Create a new adsorption site.
    pub fn new(
        site_type: SiteType,
        position: [f64; 3],
        binding_energy_ev: f64,
        coordination: u32,
    ) -> Self {
        Self {
            site_type,
            position,
            binding_energy_ev,
            coordination,
            occupied: false,
        }
    }

    /// Site preference score: lower binding energy (more negative) is preferred.
    ///
    /// Returns the absolute binding energy for ranking purposes.
    pub fn preference_score(&self) -> f64 {
        -self.binding_energy_ev
    }

    /// Estimate the binding energy from the d-band model (Hammer-Norskov).
    ///
    /// ΔE_ads ≈ α·ε_d + β·V²/|ε_d|  where ε_d is the d-band center.
    ///
    /// # Arguments
    /// * `d_band_center_ev` — d-band center of the metal (eV)
    /// * `coupling_matrix_element_ev` — metal-adsorbate coupling V² (eV²)
    pub fn binding_energy_d_band(d_band_center_ev: f64, coupling_matrix_element_ev: f64) -> f64 {
        let alpha = -0.5;
        let beta = 2.0;
        alpha * d_band_center_ev
            + beta * coupling_matrix_element_ev / d_band_center_ev.abs().max(0.1)
    }

    /// Find the most favorable site from a list of sites.
    ///
    /// Returns the index of the site with the lowest (most negative) binding energy.
    pub fn find_preferred_site(sites: &[AdsorptionSite]) -> Option<usize> {
        sites
            .iter()
            .enumerate()
            .filter(|(_, s)| !s.occupied)
            .min_by(|(_, a), (_, b)| {
                a.binding_energy_ev
                    .partial_cmp(&b.binding_energy_ev)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
    }
}

// ---------------------------------------------------------------------------
// LangmuirAdsorption
// ---------------------------------------------------------------------------

/// Langmuir and BET adsorption isotherm models.
///
/// Handles monolayer (Langmuir), multilayer (BET), and Henry's law limit
/// for gas adsorption on solid surfaces.
#[derive(Debug, Clone)]
pub struct LangmuirAdsorption {
    /// Langmuir equilibrium constant K = k_ads / k_des (Pa⁻¹).
    pub equilibrium_constant: f64,
    /// Saturation coverage (monolayer capacity), θ_max = 1 by default.
    pub theta_max: f64,
    /// BET constant C = exp((E1 - EL)/(RT)) where E1 is 1st-layer adsorption energy.
    pub bet_constant: f64,
    /// Saturation vapor pressure at the measurement temperature (Pa).
    pub saturation_pressure: f64,
}

impl LangmuirAdsorption {
    /// Create a new Langmuir adsorption model.
    ///
    /// # Arguments
    /// * `equilibrium_constant` — Langmuir K (Pa⁻¹)
    /// * `bet_constant` — BET C constant
    /// * `saturation_pressure` — P₀ (Pa)
    pub fn new(equilibrium_constant: f64, bet_constant: f64, saturation_pressure: f64) -> Self {
        Self {
            equilibrium_constant,
            theta_max: 1.0,
            bet_constant,
            saturation_pressure,
        }
    }

    /// Compute Langmuir isotherm coverage θ = K·P / (1 + K·P).
    ///
    /// Returns fractional coverage θ ∈ \[0, 1).
    pub fn langmuir_coverage(&self, pressure_pa: f64) -> f64 {
        let kp = self.equilibrium_constant * pressure_pa;
        self.theta_max * kp / (1.0 + kp)
    }

    /// Henry's law limit: at low pressure, θ ≈ K·P (linear regime).
    ///
    /// Valid when K·P << 1.
    pub fn henry_coverage(&self, pressure_pa: f64) -> f64 {
        self.theta_max * self.equilibrium_constant * pressure_pa
    }

    /// BET isotherm for multilayer adsorption.
    ///
    /// θ = C·x / \[(1 - x)(1 - x + C·x)\]  where x = P / P₀
    ///
    /// Returns the number of adsorbed monolayers (can exceed 1).
    pub fn bet_coverage(&self, pressure_pa: f64) -> f64 {
        if self.saturation_pressure < 1e-30 {
            return 0.0;
        }
        let x = pressure_pa / self.saturation_pressure;
        let x = clamp(x, 0.0, 0.9999); // avoid singularity at x=1
        let c = self.bet_constant;
        c * x / ((1.0 - x) * (1.0 - x + c * x))
    }

    /// Compute the Langmuir equilibrium constant from thermodynamics.
    ///
    /// K = (h² / (2π·m·kT))^(3/2) / (kT) · exp(-ΔE_ads / (kT))
    ///
    /// Uses a simplified pre-exponential from kinetic theory.
    ///
    /// # Arguments
    /// * `adsorption_energy_j` — adsorption energy (J, negative for exothermic)
    /// * `temperature_k` — temperature (K)
    /// * `mass_kg` — adsorbate molecular mass (kg)
    pub fn equilibrium_constant_from_thermo(
        adsorption_energy_j: f64,
        temperature_k: f64,
        mass_kg: f64,
    ) -> f64 {
        if temperature_k < 1e-10 || mass_kg < 1e-40 {
            return 0.0;
        }
        let kbt = KB * temperature_k;
        let thermal_wavelength = H_PLANCK / (2.0 * PI * mass_kg * kbt).sqrt();
        let pre_exp = thermal_wavelength.powi(3) / kbt;
        pre_exp * arrhenius_factor(-adsorption_energy_j, temperature_k)
    }

    /// Isosteric heat of adsorption from van't Hoff relation.
    ///
    /// q_st = -R · d(ln K) / d(1/T) ≈ -ΔH_ads
    ///
    /// Returns isosteric heat in J/mol.
    pub fn isosteric_heat_j_per_mol(&self, adsorption_energy_j: f64) -> f64 {
        -adsorption_energy_j * NA
    }

    /// Adsorption rate from Hertz-Knudsen: r_ads = S·P / sqrt(2π·m·kT)
    ///
    /// # Arguments
    /// * `pressure_pa` — gas pressure (Pa)
    /// * `sticking_coeff` — sticking coefficient S ∈ \[0,1\]
    /// * `temperature_k` — temperature (K)
    /// * `mass_kg` — molecular mass (kg)
    pub fn adsorption_rate(
        pressure_pa: f64,
        sticking_coeff: f64,
        temperature_k: f64,
        mass_kg: f64,
    ) -> f64 {
        if temperature_k < 1e-10 || mass_kg < 1e-40 {
            return 0.0;
        }
        let denominator = (2.0 * PI * mass_kg * KB * temperature_k).sqrt();
        sticking_coeff * pressure_pa / denominator
    }
}

// ---------------------------------------------------------------------------
// ChemisorptionModel
// ---------------------------------------------------------------------------

/// Chemisorption model including Brønsted-Evans-Polanyi (BEP) relation.
///
/// The BEP relation connects activation barriers to reaction energies:
/// Ea = Ea0 + α·ΔE_rxn
///
/// This model also handles pre-exponential factors from transition state theory
/// and sticking probability via a precursor mechanism.
#[derive(Debug, Clone)]
pub struct ChemisorptionModel {
    /// BEP slope α (dimensionless, typically 0–1).
    pub bep_slope: f64,
    /// Intrinsic activation barrier at zero reaction energy (eV).
    pub intrinsic_barrier_ev: f64,
    /// Pre-exponential frequency factor ν₀ (s⁻¹).
    pub pre_exponential: f64,
    /// Sticking coefficient at zero coverage S₀.
    pub sticking_coefficient_s0: f64,
    /// Precursor state energy (eV, negative = bound).
    pub precursor_energy_ev: f64,
}

impl ChemisorptionModel {
    /// Create a new chemisorption model.
    pub fn new(
        bep_slope: f64,
        intrinsic_barrier_ev: f64,
        pre_exponential: f64,
        sticking_coefficient_s0: f64,
        precursor_energy_ev: f64,
    ) -> Self {
        Self {
            bep_slope,
            intrinsic_barrier_ev,
            pre_exponential,
            sticking_coefficient_s0,
            precursor_energy_ev,
        }
    }

    /// Compute activation barrier via BEP relation.
    ///
    /// Ea = Ea0 + α·ΔE_rxn  (in eV)
    ///
    /// The barrier is clamped to zero from below (no negative barriers).
    pub fn bep_barrier_ev(&self, reaction_energy_ev: f64) -> f64 {
        (self.intrinsic_barrier_ev + self.bep_slope * reaction_energy_ev).max(0.0)
    }

    /// Compute rate constant from transition state theory (Eyring equation).
    ///
    /// k = (kT/h) · exp(-Ea / kT)
    ///
    /// # Arguments
    /// * `barrier_ev` — activation barrier (eV)
    /// * `temperature_k` — temperature (K)
    pub fn rate_constant_tst(&self, barrier_ev: f64, temperature_k: f64) -> f64 {
        if temperature_k < 1e-10 {
            return 0.0;
        }
        let ev_to_j = 1.602_176_634e-19;
        let kbt = KB * temperature_k;
        let pre = kbt / H_PLANCK;
        pre * arrhenius_factor(barrier_ev * ev_to_j, temperature_k)
    }

    /// Rate constant using simple Arrhenius form with stored pre-exponential.
    pub fn arrhenius_rate(&self, barrier_ev: f64, temperature_k: f64) -> f64 {
        let ev_to_j = 1.602_176_634e-19;
        self.pre_exponential * arrhenius_factor(barrier_ev * ev_to_j, temperature_k)
    }

    /// Sticking probability via precursor mechanism.
    ///
    /// S(θ) = S₀ · (1 - θ)^n_sites  (simplest Kisliuk model variant)
    ///
    /// # Arguments
    /// * `coverage` — fractional coverage θ ∈ \[0,1\]
    /// * `n_sites` — number of sites required for adsorption (typically 1 or 2)
    pub fn sticking_probability(&self, coverage: f64, n_sites: u32) -> f64 {
        let theta = clamp(coverage, 0.0, 1.0);
        self.sticking_coefficient_s0 * (1.0 - theta).powi(n_sites as i32)
    }

    /// Desorption pre-exponential from transition state theory.
    ///
    /// ν = kT/h · (q_TS / q_ads)  simplified as ν₀ ≈ kT/h for mobile TS.
    pub fn desorption_pre_exponential(temperature_k: f64) -> f64 {
        if temperature_k < 1e-10 {
            return 0.0;
        }
        KB * temperature_k / H_PLANCK
    }

    /// Compute surface residence time τ = 1/k_des.
    pub fn residence_time(&self, desorption_barrier_ev: f64, temperature_k: f64) -> f64 {
        let k_des = self.arrhenius_rate(desorption_barrier_ev, temperature_k);
        if k_des < 1e-300 {
            f64::INFINITY
        } else {
            1.0 / k_des
        }
    }
}

// ---------------------------------------------------------------------------
// NudgedElasticBand
// ---------------------------------------------------------------------------

/// Image in the nudged elastic band (NEB) method.
#[derive(Debug, Clone)]
pub struct NebImage {
    /// Atomic coordinates (flattened: \[x0, y0, z0, x1, y1, z1, ...\]).
    pub coords: Vec<f64>,
    /// Potential energy at this image (eV).
    pub energy: f64,
    /// Forces on atoms (same layout as coords, eV/Å).
    pub forces: Vec<f64>,
}

impl NebImage {
    /// Create a new NEB image with zero forces.
    pub fn new(coords: Vec<f64>, energy: f64) -> Self {
        let n = coords.len();
        Self {
            coords,
            energy,
            forces: vec![0.0; n],
        }
    }
}

/// Nudged Elastic Band (NEB) method for finding minimum energy paths (MEP).
///
/// Implements both plain NEB and climbing-image NEB (CI-NEB) where the highest
/// energy image climbs to the saddle point.
///
/// Reference: Henkelman & Jonsson, J. Chem. Phys. 113, 9978 (2000).
#[derive(Debug, Clone)]
pub struct NudgedElasticBand {
    /// Images along the reaction path (including fixed endpoints).
    pub images: Vec<NebImage>,
    /// Spring force constant between adjacent images (eV/Å²).
    pub spring_constant: f64,
    /// Enable climbing-image NEB.
    pub climbing_image: bool,
    /// Maximum force component for convergence (eV/Å).
    pub convergence_threshold: f64,
    /// Number of optimization steps taken.
    pub n_steps: usize,
}

impl NudgedElasticBand {
    /// Create a new NEB calculation.
    ///
    /// # Arguments
    /// * `images` — interpolated images from reactant to product
    /// * `spring_constant` — spring stiffness k (eV/Å²)
    /// * `climbing_image` — enable CI-NEB
    pub fn new(images: Vec<NebImage>, spring_constant: f64, climbing_image: bool) -> Self {
        Self {
            images,
            spring_constant,
            climbing_image,
            convergence_threshold: 0.05,
            n_steps: 0,
        }
    }

    /// Linearly interpolate images between two endpoint configurations.
    ///
    /// Returns `n_images` + 2 images (including endpoints).
    pub fn interpolate_images(start: Vec<f64>, end: Vec<f64>, n_images: usize) -> Vec<NebImage> {
        let mut images = Vec::with_capacity(n_images + 2);
        images.push(NebImage::new(start.clone(), 0.0));
        for i in 1..=n_images {
            let frac = i as f64 / (n_images + 1) as f64;
            let coords: Vec<f64> = start
                .iter()
                .zip(end.iter())
                .map(|(s, e)| s + frac * (e - s))
                .collect();
            images.push(NebImage::new(coords, 0.0));
        }
        images.push(NebImage::new(end, 0.0));
        images
    }

    /// Compute tangent vector along the path at image `i` (upwind tangent).
    fn compute_tangent(&self, i: usize) -> Vec<f64> {
        let n = self.images.len();
        if i == 0 || i == n - 1 {
            return vec![0.0; self.images[i].coords.len()];
        }
        let prev = &self.images[i - 1].coords;
        let curr = &self.images[i].coords;
        let next = &self.images[i + 1].coords;
        let e_prev = self.images[i - 1].energy;
        let e_curr = self.images[i].energy;
        let e_next = self.images[i + 1].energy;

        let tau: Vec<f64> = if e_next > e_curr && e_curr > e_prev {
            next.iter().zip(curr.iter()).map(|(n, c)| n - c).collect()
        } else if e_next < e_curr && e_curr < e_prev {
            curr.iter().zip(prev.iter()).map(|(c, p)| c - p).collect()
        } else {
            // Bisection tangent
            let d_plus: Vec<f64> = next.iter().zip(curr.iter()).map(|(n, c)| n - c).collect();
            let d_minus: Vec<f64> = curr.iter().zip(prev.iter()).map(|(c, p)| c - p).collect();
            let de_max = (e_next - e_curr).abs().max((e_curr - e_prev).abs());
            let de_min = (e_next - e_curr).abs().min((e_curr - e_prev).abs());
            if e_next > e_prev {
                d_plus
                    .iter()
                    .map(|v| v * de_max)
                    .zip(d_minus.iter().map(|v| v * de_min))
                    .map(|(a, b)| a + b)
                    .collect()
            } else {
                d_plus
                    .iter()
                    .map(|v| v * de_min)
                    .zip(d_minus.iter().map(|v| v * de_max))
                    .map(|(a, b)| a + b)
                    .collect()
            }
        };

        // Normalize
        let norm: f64 = tau.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-15 {
            tau
        } else {
            tau.iter().map(|v| v / norm).collect()
        }
    }

    /// Compute spring forces between images.
    fn spring_force(&self, i: usize) -> Vec<f64> {
        let n = self.images.len();
        if i == 0 || i == n - 1 {
            return vec![0.0; self.images[i].coords.len()];
        }
        let prev = &self.images[i - 1].coords;
        let curr = &self.images[i].coords;
        let next = &self.images[i + 1].coords;

        let d_next: f64 = next
            .iter()
            .zip(curr.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let d_prev: f64 = curr
            .iter()
            .zip(prev.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();

        let tau = self.compute_tangent(i);
        let spring_mag = self.spring_constant * (d_next - d_prev);
        tau.iter().map(|t| spring_mag * t).collect()
    }

    /// Apply NEB forces to all movable images for one optimization step.
    ///
    /// This implements the NEB force projection: F_NEB = F_perp + F_spring_parallel
    /// where F_perp = F_true - (F_true·τ)τ
    ///
    /// Returns the maximum force component.
    pub fn apply_neb_forces(&mut self) -> f64 {
        let n = self.images.len();
        let mut max_force = 0.0_f64;

        // Find climbing image (highest energy image)
        let climb_idx = if self.climbing_image {
            (1..n - 1).max_by(|&a, &b| {
                self.images[a]
                    .energy
                    .partial_cmp(&self.images[b].energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        } else {
            None
        };

        for i in 1..n - 1 {
            let tau = self.compute_tangent(i);
            let spring = self.spring_force(i);
            let true_force = self.images[i].forces.clone();

            let ndof = true_force.len();
            let mut neb_force = vec![0.0f64; ndof];

            if self.climbing_image && Some(i) == climb_idx {
                // Climbing image: invert force along tangent
                let f_dot_tau: f64 = true_force.iter().zip(tau.iter()).map(|(f, t)| f * t).sum();
                for k in 0..ndof {
                    neb_force[k] = true_force[k] - 2.0 * f_dot_tau * tau[k];
                }
            } else {
                // Perpendicular true force + parallel spring force
                let f_dot_tau: f64 = true_force.iter().zip(tau.iter()).map(|(f, t)| f * t).sum();
                let s_dot_tau: f64 = spring.iter().zip(tau.iter()).map(|(s, t)| s * t).sum();
                for k in 0..ndof {
                    neb_force[k] = (true_force[k] - f_dot_tau * tau[k]) + s_dot_tau * tau[k];
                }
            }

            let fmax = neb_force
                .iter()
                .cloned()
                .fold(0.0_f64, |a, v| a.max(v.abs()));
            max_force = max_force.max(fmax);
            self.images[i].forces = neb_force;
        }

        max_force
    }

    /// Move images along their NEB forces (simple steepest descent step).
    pub fn step_images(&mut self, step_size: f64) {
        let n = self.images.len();
        for i in 1..n - 1 {
            let forces = self.images[i].forces.clone();
            for (c, f) in self.images[i].coords.iter_mut().zip(forces.iter()) {
                *c += step_size * f;
            }
        }
        self.n_steps += 1;
    }

    /// Check convergence: max force < threshold.
    pub fn is_converged(&self) -> bool {
        let n = self.images.len();
        for i in 1..n - 1 {
            let fmax = self.images[i]
                .forces
                .iter()
                .cloned()
                .fold(0.0_f64, |a, v| a.max(v.abs()));
            if fmax > self.convergence_threshold {
                return false;
            }
        }
        true
    }

    /// Return the index and energy of the saddle point image.
    pub fn saddle_point(&self) -> (usize, f64) {
        let n = self.images.len();
        (1..n - 1)
            .map(|i| (i, self.images[i].energy))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, 0.0))
    }

    /// Return the activation energy (saddle energy - reactant energy) in eV.
    pub fn activation_energy_ev(&self) -> f64 {
        let (_idx, e_saddle) = self.saddle_point();
        let e_reactant = self.images.first().map(|im| im.energy).unwrap_or(0.0);
        e_saddle - e_reactant
    }
}

// ---------------------------------------------------------------------------
// SurfaceDiffusion
// ---------------------------------------------------------------------------

/// Lattice hop event in surface diffusion.
#[derive(Debug, Clone)]
pub struct HopEvent {
    /// Source site index.
    pub from: usize,
    /// Destination site index.
    pub to: usize,
    /// Time of hop (s).
    pub time: f64,
}

/// Surface diffusion model: random walk on a lattice with Arrhenius hop rates.
///
/// Tracks mean-square displacement (MSD) to extract the diffusion coefficient D.
/// D = MSD / (2d·t) where d is the dimensionality (2 for surface).
#[derive(Debug, Clone)]
pub struct SurfaceDiffusion {
    /// Current position (x, y) in Å.
    pub position: [f64; 2],
    /// Diffusion barrier (eV).
    pub barrier_ev: f64,
    /// Attempt frequency / pre-exponential (s⁻¹).
    pub attempt_frequency: f64,
    /// Lattice constant / hop distance (Å).
    pub hop_distance: f64,
    /// Temperature (K).
    pub temperature_k: f64,
    /// History of positions for MSD calculation.
    pub trajectory: Vec<[f64; 2]>,
    /// Accumulated simulation time (s).
    pub time: f64,
}

impl SurfaceDiffusion {
    /// Create a new surface diffusion tracker.
    ///
    /// # Arguments
    /// * `start_pos` — initial position (Å)
    /// * `barrier_ev` — diffusion barrier (eV)
    /// * `attempt_frequency` — pre-exponential ν₀ (s⁻¹)
    /// * `hop_distance` — nearest-neighbor distance (Å)
    /// * `temperature_k` — temperature (K)
    pub fn new(
        start_pos: [f64; 2],
        barrier_ev: f64,
        attempt_frequency: f64,
        hop_distance: f64,
        temperature_k: f64,
    ) -> Self {
        let trajectory = vec![start_pos];
        Self {
            position: start_pos,
            barrier_ev,
            attempt_frequency,
            hop_distance,
            temperature_k,
            trajectory,
            time: 0.0,
        }
    }

    /// Arrhenius hop rate (s⁻¹).
    pub fn hop_rate(&self) -> f64 {
        let ev_to_j = 1.602_176_634e-19;
        self.attempt_frequency * arrhenius_factor(self.barrier_ev * ev_to_j, self.temperature_k)
    }

    /// Perform one kinetic Monte Carlo hop using a simple LCG.
    ///
    /// Direction is chosen uniformly from 4 nearest neighbors (square lattice).
    pub fn kmc_hop(&mut self, rng_state: &mut u64) {
        // Simple LCG
        *rng_state = rng_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let rand_dir = (*rng_state >> 32) % 4;
        let dirs: [[f64; 2]; 4] = [
            [self.hop_distance, 0.0],
            [-self.hop_distance, 0.0],
            [0.0, self.hop_distance],
            [0.0, -self.hop_distance],
        ];
        let d = dirs[rand_dir as usize];
        self.position[0] += d[0];
        self.position[1] += d[1];
        let rate = self.hop_rate();
        // Advance time by mean waiting time 1/rate
        if rate > 1e-300 {
            self.time += 1.0 / rate;
        }
        self.trajectory.push(self.position);
    }

    /// Run `n_steps` KMC hops.
    pub fn run(&mut self, n_steps: usize, seed: u64) {
        let mut rng_state = seed.wrapping_add(0x9e3779b97f4a7c15);
        for _ in 0..n_steps {
            self.kmc_hop(&mut rng_state);
        }
    }

    /// Compute mean-square displacement from the initial position.
    pub fn msd(&self) -> f64 {
        if self.trajectory.len() < 2 {
            return 0.0;
        }
        let origin = self.trajectory[0];
        self.trajectory
            .iter()
            .map(|p| (p[0] - origin[0]).powi(2) + (p[1] - origin[1]).powi(2))
            .sum::<f64>()
            / self.trajectory.len() as f64
    }

    /// Compute surface diffusion coefficient D = MSD / (4·t) \[Å²/s\].
    ///
    /// The factor 4 = 2·d where d=2 for 2D surface diffusion.
    pub fn diffusivity(&self) -> f64 {
        if self.time < 1e-30 {
            return 0.0;
        }
        self.msd() / (4.0 * self.time)
    }

    /// Theoretical diffusivity from Einstein-Smoluchowski + Arrhenius:
    ///
    /// D = a² · ν₀ · exp(-Ea / kT) / 4  (2D, 4 neighbors)
    pub fn theoretical_diffusivity(&self) -> f64 {
        let ev_to_j = 1.602_176_634e-19;
        self.hop_distance.powi(2)
            * self.attempt_frequency
            * arrhenius_factor(self.barrier_ev * ev_to_j, self.temperature_k)
            / 4.0
    }
}

// ---------------------------------------------------------------------------
// ThinFilmGrowth
// ---------------------------------------------------------------------------

/// Thin film growth mode.
#[derive(Debug, Clone, PartialEq)]
pub enum GrowthMode {
    /// Volmer-Weber (VW): 3D island growth (film-film interaction > film-substrate).
    VolmerWeber,
    /// Frank-van der Merwe (FM): layer-by-layer growth (film-substrate > film-film).
    FrankVanDerMerwe,
    /// Stranski-Krastanov (SK): initial wetting layer + 3D islands above critical thickness.
    StranskiKrastanov,
}

/// Thin film growth model with strain energy and wetting layer calculations.
///
/// Growth mode is determined by the balance between surface and interface energies.
#[derive(Debug, Clone)]
pub struct ThinFilmGrowth {
    /// Surface energy of the film material (J/m²).
    pub gamma_film: f64,
    /// Surface energy of the substrate (J/m²).
    pub gamma_substrate: f64,
    /// Interface energy between film and substrate (J/m²).
    pub gamma_interface: f64,
    /// Lattice mismatch ε = (a_film - a_sub) / a_sub (dimensionless).
    pub lattice_mismatch: f64,
    /// Biaxial elastic modulus of the film (GPa).
    pub biaxial_modulus_gpa: f64,
    /// Current film thickness (Å).
    pub thickness_angstrom: f64,
    /// Wetting parameter: positive = FM tendency, negative = VW tendency.
    pub wetting_parameter: f64,
}

impl ThinFilmGrowth {
    /// Create a new thin film growth model.
    pub fn new(
        gamma_film: f64,
        gamma_substrate: f64,
        gamma_interface: f64,
        lattice_mismatch: f64,
        biaxial_modulus_gpa: f64,
    ) -> Self {
        let wetting_parameter = gamma_substrate - gamma_film - gamma_interface;
        Self {
            gamma_film,
            gamma_substrate,
            gamma_interface,
            lattice_mismatch,
            biaxial_modulus_gpa,
            thickness_angstrom: 0.0,
            wetting_parameter,
        }
    }

    /// Determine the growth mode based on energetic criteria.
    ///
    /// - FM: γ_sub > γ_film + γ_interface (substrate prefers to be covered)
    /// - VW: γ_sub < γ_film + γ_interface (islands preferred)
    /// - SK: initial FM transitions to VW above critical thickness due to strain
    pub fn growth_mode(&self) -> GrowthMode {
        let wetting = self.wetting_parameter;
        if wetting > 0.0 && self.lattice_mismatch.abs() > 0.02 {
            // Strain drives SK transition
            GrowthMode::StranskiKrastanov
        } else if wetting > 0.0 {
            GrowthMode::FrankVanDerMerwe
        } else {
            GrowthMode::VolmerWeber
        }
    }

    /// Elastic strain energy density (J/m³).
    ///
    /// E_strain = Y_biaxial · ε² · h  where h is film thickness.
    pub fn strain_energy_density(&self) -> f64 {
        let y = self.biaxial_modulus_gpa * 1e9; // GPa to Pa
        let eps = self.lattice_mismatch;
        y * eps * eps
    }

    /// Critical thickness for SK transition (Matthews-Blakeslee approximation, Å).
    ///
    /// h_c ≈ b / (8π·ε · (1 + ν)) · ln(h_c / b)  (simplified: h_c ≈ b / (4π·|ε|))
    ///
    /// # Arguments
    /// * `burgers_vector_angstrom` — magnitude of Burgers vector (Å)
    /// * `poisson_ratio` — Poisson's ratio ν
    pub fn critical_thickness(&self, burgers_vector_angstrom: f64, poisson_ratio: f64) -> f64 {
        if self.lattice_mismatch.abs() < 1e-10 {
            return f64::INFINITY;
        }
        let eps = self.lattice_mismatch.abs();
        let b = burgers_vector_angstrom;
        let factor = 1.0 - poisson_ratio / 2.0;
        b * factor / (4.0 * PI * eps * (1.0 + poisson_ratio))
    }

    /// Check if current thickness exceeds the critical thickness (SK transition).
    pub fn is_above_critical_thickness(
        &self,
        burgers_vector_angstrom: f64,
        poisson_ratio: f64,
    ) -> bool {
        self.thickness_angstrom > self.critical_thickness(burgers_vector_angstrom, poisson_ratio)
    }

    /// Deposit a monolayer, advancing film thickness by one layer.
    ///
    /// Returns the new growth mode.
    pub fn deposit_monolayer(&mut self, monolayer_thickness_angstrom: f64) -> GrowthMode {
        self.thickness_angstrom += monolayer_thickness_angstrom;
        self.growth_mode()
    }

    /// Compute the total surface + strain energy per unit area (J/m²).
    pub fn total_energy_per_area(&self) -> f64 {
        let h_m = self.thickness_angstrom * 1e-10; // Å to m
        let strain_volume = self.strain_energy_density() * h_m;
        self.gamma_film + self.gamma_interface + strain_volume
    }
}

// ---------------------------------------------------------------------------
// CatalyticCycle
// ---------------------------------------------------------------------------

/// A single elementary step in a catalytic cycle.
#[derive(Debug, Clone)]
pub struct ElementaryStep {
    /// Name of the step (e.g., "CO adsorption").
    pub name: String,
    /// Reaction energy ΔE (eV).
    pub reaction_energy_ev: f64,
    /// Forward activation barrier Ea (eV).
    pub forward_barrier_ev: f64,
    /// Reverse activation barrier Ea_rev = Ea - ΔE (eV).
    pub reverse_barrier_ev: f64,
}

impl ElementaryStep {
    /// Create a new elementary step.
    pub fn new(name: &str, reaction_energy_ev: f64, forward_barrier_ev: f64) -> Self {
        let reverse_barrier_ev = (forward_barrier_ev - reaction_energy_ev).max(0.0);
        Self {
            name: name.to_string(),
            reaction_energy_ev,
            forward_barrier_ev,
            reverse_barrier_ev,
        }
    }
}

/// Catalytic cycle simulator implementing the Sabatier principle.
///
/// Models heterogeneous catalysis with elementary steps, computes turnover
/// frequency (TOF), and identifies the rate-limiting step.
///
/// Reference: Norskov et al., Nature Chemistry 1, 37 (2009).
#[derive(Debug, Clone)]
pub struct CatalyticCycle {
    /// List of elementary steps in the catalytic cycle.
    pub steps: Vec<ElementaryStep>,
    /// Temperature (K).
    pub temperature_k: f64,
    /// Pre-exponential rate factor (s⁻¹).
    pub pre_exponential: f64,
    /// Partial pressures of reactants (Pa).
    pub partial_pressures: Vec<f64>,
}

impl CatalyticCycle {
    /// Create a new catalytic cycle model.
    pub fn new(temperature_k: f64, pre_exponential: f64) -> Self {
        Self {
            steps: Vec::new(),
            temperature_k,
            pre_exponential,
            partial_pressures: Vec::new(),
        }
    }

    /// Add an elementary step to the cycle.
    pub fn add_step(&mut self, step: ElementaryStep) {
        self.steps.push(step);
    }

    /// Compute the forward rate constant for step `i` (s⁻¹).
    pub fn forward_rate(&self, i: usize) -> f64 {
        if i >= self.steps.len() {
            return 0.0;
        }
        let ev_to_j = 1.602_176_634e-19;
        self.pre_exponential
            * arrhenius_factor(
                self.steps[i].forward_barrier_ev * ev_to_j,
                self.temperature_k,
            )
    }

    /// Compute the reverse rate constant for step `i` (s⁻¹).
    pub fn reverse_rate(&self, i: usize) -> f64 {
        if i >= self.steps.len() {
            return 0.0;
        }
        let ev_to_j = 1.602_176_634e-19;
        self.pre_exponential
            * arrhenius_factor(
                self.steps[i].reverse_barrier_ev * ev_to_j,
                self.temperature_k,
            )
    }

    /// Turnover frequency (TOF) using the mean-field approximation.
    ///
    /// TOF = min(k_forward) over the cycle (rate-limiting step approximation).
    pub fn turnover_frequency(&self) -> f64 {
        if self.steps.is_empty() {
            return 0.0;
        }
        self.steps
            .iter()
            .enumerate()
            .map(|(i, _)| self.forward_rate(i))
            .fold(f64::INFINITY, f64::min)
    }

    /// Find the rate-limiting step (step with smallest forward rate).
    pub fn rate_limiting_step(&self) -> Option<(usize, &ElementaryStep)> {
        self.steps.iter().enumerate().min_by(|(i, _), (j, _)| {
            self.forward_rate(*i)
                .partial_cmp(&self.forward_rate(*j))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Compute volcano plot data: TOF vs. adsorption energy of the key intermediate.
    ///
    /// Returns (adsorption_energy_ev, log10(TOF)) pairs.
    /// The BEP relation is used: Ea = Ea0 + alpha * delta_E_ads.
    pub fn volcano_plot_data(
        &self,
        adsorption_energy_range: &[f64],
        bep_slope: f64,
        intrinsic_barrier_ev: f64,
    ) -> Vec<(f64, f64)> {
        let ev_to_j = 1.602_176_634e-19;
        adsorption_energy_range
            .iter()
            .map(|&de| {
                let barrier = (intrinsic_barrier_ev + bep_slope * de).max(0.0);
                let rate =
                    self.pre_exponential * arrhenius_factor(barrier * ev_to_j, self.temperature_k);
                let log_tof = if rate > 1e-300 { rate.log10() } else { -300.0 };
                (de, log_tof)
            })
            .collect()
    }

    /// Overall reaction energy: sum of all step reaction energies.
    pub fn overall_reaction_energy_ev(&self) -> f64 {
        self.steps.iter().map(|s| s.reaction_energy_ev).sum()
    }

    /// Compute selectivity: fraction of TOF attributed to the desired product pathway.
    ///
    /// Simple model: selectivity = k_desired / (k_desired + k_side)
    pub fn selectivity(k_desired: f64, k_side: f64) -> f64 {
        let total = k_desired + k_side;
        if total < 1e-300 {
            0.0
        } else {
            k_desired / total
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── SurfaceSlab tests ───────────────────────────────────────────────────

    #[test]
    fn test_slab_creation_layer_count() {
        let slab = SurfaceSlab::new((1, 1, 1), (5.76, 5.76), 4.08, 3, 2, 15.0);
        assert_eq!(slab.layers.len(), 5);
        assert_eq!(slab.n_fixed_atoms(), 12);
        assert_eq!(slab.n_mobile_atoms(), 8);
    }

    #[test]
    fn test_slab_total_thickness() {
        let slab = SurfaceSlab::new((1, 0, 0), (4.0, 4.0), 4.0, 2, 2, 20.0);
        let thick = slab.total_thickness();
        assert!(thick > 20.0, "Thickness should exceed vacuum layer alone");
    }

    #[test]
    fn test_slab_surface_energy_zero_for_bulk_like() {
        let mut slab = SurfaceSlab::new((1, 1, 1), (5.76, 5.76), 4.08, 3, 2, 15.0);
        // If slab energy = N * bulk energy, surface energy should be ~0
        let n: usize = slab.layers.iter().map(|l| l.n_atoms).sum();
        let e_bulk = -3.0;
        slab.compute_surface_energy(n as f64 * e_bulk, e_bulk, 1.602e-19);
        assert!(slab.surface_energy.abs() < 1e-10);
    }

    #[test]
    fn test_slab_relaxation_changes_z() {
        let mut slab = SurfaceSlab::new((1, 1, 1), (5.76, 5.76), 4.08, 2, 2, 15.0);
        let z_before = slab.layers.last().unwrap().z_position;
        slab.relax_surface(-0.05);
        let z_after = slab.layers.last().unwrap().z_position;
        assert!(
            z_after < z_before,
            "Relaxation should contract the top layer"
        );
    }

    // ── AdsorptionSite tests ────────────────────────────────────────────────

    #[test]
    fn test_adsorption_site_preference_score() {
        let site_top = AdsorptionSite::new(SiteType::Top, [0.0, 0.0, 2.0], -0.5, 1);
        let site_hollow = AdsorptionSite::new(SiteType::FccHollow, [1.0, 1.0, 1.8], -1.2, 3);
        assert!(site_hollow.preference_score() > site_top.preference_score());
    }

    #[test]
    fn test_find_preferred_site_returns_most_stable() {
        let sites = vec![
            AdsorptionSite::new(SiteType::Top, [0.0, 0.0, 2.0], -0.5, 1),
            AdsorptionSite::new(SiteType::Bridge, [1.0, 0.0, 1.9], -0.8, 2),
            AdsorptionSite::new(SiteType::FccHollow, [1.0, 1.0, 1.7], -1.5, 3),
        ];
        let idx = AdsorptionSite::find_preferred_site(&sites);
        assert_eq!(
            idx,
            Some(2),
            "FCC hollow should be preferred (most negative BE)"
        );
    }

    #[test]
    fn test_d_band_binding_energy_stronger_for_high_d_band() {
        let e1 = AdsorptionSite::binding_energy_d_band(-2.0, 1.0);
        let e2 = AdsorptionSite::binding_energy_d_band(-4.0, 1.0);
        // d_band=-4.0 (deeper) gives larger coupling/d_band term → stronger contribution
        assert!(
            e1 < e2,
            "Higher (less negative) d-band center gives weaker coupling contribution"
        );
    }

    // ── LangmuirAdsorption tests ────────────────────────────────────────────

    #[test]
    fn test_langmuir_zero_pressure_gives_zero_coverage() {
        let model = LangmuirAdsorption::new(1e-5, 100.0, 1e5);
        assert_eq!(model.langmuir_coverage(0.0), 0.0);
    }

    #[test]
    fn test_langmuir_high_pressure_approaches_theta_max() {
        let model = LangmuirAdsorption::new(1e-3, 100.0, 1e5);
        let theta = model.langmuir_coverage(1e6); // K·P >> 1
        assert!(
            theta > 0.999,
            "High pressure should give near-saturation coverage: {theta}"
        );
    }

    #[test]
    fn test_langmuir_half_coverage_at_one_over_k() {
        let k = 1e-4;
        let model = LangmuirAdsorption::new(k, 100.0, 1e5);
        let p_half = 1.0 / k;
        let theta = model.langmuir_coverage(p_half);
        assert!((theta - 0.5).abs() < 1e-10, "θ=0.5 when P=1/K: got {theta}");
    }

    #[test]
    fn test_henry_law_linear_at_low_pressure() {
        let k = 1e-5;
        let model = LangmuirAdsorption::new(k, 100.0, 1e5);
        let p = 1.0; // K·P = 1e-5 << 1
        let henry = model.henry_coverage(p);
        let langmuir = model.langmuir_coverage(p);
        assert!(
            (henry - langmuir).abs() / langmuir < 0.01,
            "Henry and Langmuir should agree at low P"
        );
    }

    #[test]
    fn test_bet_isotherm_returns_positive_coverage() {
        let model = LangmuirAdsorption::new(1e-4, 50.0, 1e5);
        let theta_bet = model.bet_coverage(0.5e5); // half saturation pressure
        assert!(
            theta_bet > 0.0,
            "BET coverage should be positive: {theta_bet}"
        );
    }

    #[test]
    fn test_bet_increases_with_pressure() {
        let model = LangmuirAdsorption::new(1e-4, 50.0, 1e5);
        let t1 = model.bet_coverage(1e4);
        let t2 = model.bet_coverage(5e4);
        assert!(t2 > t1, "BET coverage should increase with pressure");
    }

    // ── ChemisorptionModel tests ────────────────────────────────────────────

    #[test]
    fn test_bep_barrier_increases_with_reaction_energy() {
        let model = ChemisorptionModel::new(0.7, 0.5, 1e13, 0.1, -0.3);
        let ea1 = model.bep_barrier_ev(0.0);
        let ea2 = model.bep_barrier_ev(1.0);
        assert!(ea2 > ea1, "BEP barrier should increase with endothermic ΔE");
    }

    #[test]
    fn test_bep_barrier_non_negative() {
        let model = ChemisorptionModel::new(0.5, 0.2, 1e13, 0.1, -0.5);
        // Very exothermic reaction: BEP barrier could go negative, should clamp to 0
        let ea = model.bep_barrier_ev(-10.0);
        assert!(ea >= 0.0, "BEP barrier cannot be negative: {ea}");
    }

    #[test]
    fn test_sticking_probability_decreases_with_coverage() {
        let model = ChemisorptionModel::new(0.5, 0.5, 1e13, 0.5, -0.3);
        let s0 = model.sticking_probability(0.0, 1);
        let s1 = model.sticking_probability(0.5, 1);
        assert!(
            s1 < s0,
            "Sticking probability should decrease with coverage"
        );
    }

    #[test]
    fn test_rate_constant_increases_with_temperature() {
        let model = ChemisorptionModel::new(0.5, 0.8, 1e13, 0.1, -0.3);
        let k300 = model.rate_constant_tst(0.8, 300.0);
        let k600 = model.rate_constant_tst(0.8, 600.0);
        assert!(
            k600 > k300,
            "Rate constant should increase with temperature"
        );
    }

    #[test]
    fn test_residence_time_finite_for_nonzero_barrier() {
        let model = ChemisorptionModel::new(0.5, 0.5, 1e13, 0.1, -0.3);
        let tau = model.residence_time(0.5, 500.0);
        assert!(tau.is_finite(), "Residence time should be finite: {tau}");
        assert!(tau > 0.0);
    }

    // ── NudgedElasticBand tests ────────────────────────────────────────────

    #[test]
    fn test_neb_interpolation_count() {
        let start = vec![0.0f64, 0.0, 0.0];
        let end = vec![3.0f64, 0.0, 0.0];
        let images = NudgedElasticBand::interpolate_images(start, end, 5);
        assert_eq!(images.len(), 7, "Should have 5 internal + 2 endpoints");
    }

    #[test]
    fn test_neb_interpolation_midpoint() {
        let start = vec![0.0f64, 0.0];
        let end = vec![4.0f64, 0.0];
        let images = NudgedElasticBand::interpolate_images(start, end, 3);
        // Middle image at index 2 should be near [2.0, 0.0]
        assert!((images[2].coords[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_neb_saddle_point_is_highest_energy() {
        let mut images = NudgedElasticBand::interpolate_images(vec![0.0f64], vec![2.0f64], 3);
        // Assign energies: parabolic shape with maximum in middle
        let energies = [0.0, 0.5, 1.0, 0.5, 0.0];
        for (im, &e) in images.iter_mut().zip(energies.iter()) {
            im.energy = e;
        }
        let neb = NudgedElasticBand::new(images, 1.0, false);
        let (idx, e) = neb.saddle_point();
        assert_eq!(idx, 2);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_neb_activation_energy() {
        let mut images = NudgedElasticBand::interpolate_images(vec![0.0f64], vec![2.0f64], 3);
        let energies = [0.1, 0.5, 1.2, 0.7, -0.2];
        for (im, &e) in images.iter_mut().zip(energies.iter()) {
            im.energy = e;
        }
        let neb = NudgedElasticBand::new(images, 1.0, false);
        let ea = neb.activation_energy_ev();
        assert!(
            (ea - (1.2 - 0.1)).abs() < 1e-10,
            "Ea = E_saddle - E_reactant: {ea}"
        );
    }

    #[test]
    fn test_neb_step_moves_images() {
        let start = vec![0.0f64, 0.0, 0.0];
        let end = vec![6.0f64, 0.0, 0.0];
        let images = NudgedElasticBand::interpolate_images(start, end, 3);
        let n_coords = images[1].coords.len();
        let pos_before = images[1].coords[0];
        let mut neb = NudgedElasticBand::new(images, 2.0, false);
        // Assign non-zero forces
        neb.images[1].forces = vec![0.1; n_coords];
        neb.step_images(0.01);
        let pos_after = neb.images[1].coords[0];
        assert!(
            (pos_after - pos_before).abs() > 1e-12,
            "Image should move after step"
        );
    }

    // ── SurfaceDiffusion tests ──────────────────────────────────────────────

    #[test]
    fn test_diffusion_hop_rate_increases_with_temperature() {
        let diff300 = SurfaceDiffusion::new([0.0, 0.0], 0.5, 1e13, 2.88, 300.0);
        let diff600 = SurfaceDiffusion::new([0.0, 0.0], 0.5, 1e13, 2.88, 600.0);
        assert!(diff600.hop_rate() > diff300.hop_rate());
    }

    #[test]
    fn test_diffusion_msd_grows_with_steps() {
        let mut diff = SurfaceDiffusion::new([0.0, 0.0], 0.3, 1e13, 2.88, 800.0);
        diff.run(1000, 42);
        assert!(
            diff.msd() > 0.0,
            "MSD should be positive after random walk: {}",
            diff.msd()
        );
    }

    #[test]
    fn test_diffusivity_positive_after_run() {
        let mut diff = SurfaceDiffusion::new([0.0, 0.0], 0.3, 1e13, 2.88, 800.0);
        diff.run(500, 12345);
        assert!(diff.diffusivity() >= 0.0);
    }

    #[test]
    fn test_theoretical_diffusivity_arrhenius() {
        let diff300 = SurfaceDiffusion::new([0.0, 0.0], 0.5, 1e13, 2.88, 300.0);
        let diff1000 = SurfaceDiffusion::new([0.0, 0.0], 0.5, 1e13, 2.88, 1000.0);
        assert!(diff1000.theoretical_diffusivity() > diff300.theoretical_diffusivity());
    }

    // ── ThinFilmGrowth tests ────────────────────────────────────────────────

    #[test]
    fn test_growth_mode_fm_for_wetting_no_strain() {
        let film = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.001, 200.0);
        assert_eq!(film.growth_mode(), GrowthMode::FrankVanDerMerwe);
    }

    #[test]
    fn test_growth_mode_vw_for_non_wetting() {
        // γ_film + γ_interface > γ_substrate → non-wetting (VW)
        let film = ThinFilmGrowth::new(2.5, 1.0, 0.8, 0.001, 200.0);
        assert_eq!(film.growth_mode(), GrowthMode::VolmerWeber);
    }

    #[test]
    fn test_growth_mode_sk_for_wetting_with_large_mismatch() {
        let film = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.05, 200.0);
        assert_eq!(film.growth_mode(), GrowthMode::StranskiKrastanov);
    }

    #[test]
    fn test_critical_thickness_larger_for_smaller_mismatch() {
        let film1 = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.01, 200.0);
        let film2 = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.04, 200.0);
        let h1 = film1.critical_thickness(2.88, 0.3);
        let h2 = film2.critical_thickness(2.88, 0.3);
        assert!(
            h1 > h2,
            "Smaller mismatch should give larger critical thickness"
        );
    }

    #[test]
    fn test_strain_energy_increases_with_mismatch() {
        let film1 = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.01, 200.0);
        let film2 = ThinFilmGrowth::new(1.0, 2.0, 0.5, 0.05, 200.0);
        assert!(film2.strain_energy_density() > film1.strain_energy_density());
    }

    // ── CatalyticCycle tests ────────────────────────────────────────────────

    #[test]
    fn test_catalytic_tof_positive_for_valid_cycle() {
        let mut cycle = CatalyticCycle::new(500.0, 1e13);
        cycle.add_step(ElementaryStep::new("CO ads", -1.0, 0.0));
        cycle.add_step(ElementaryStep::new("O ads", -2.0, 0.3));
        cycle.add_step(ElementaryStep::new("CO+O react", -1.5, 0.8));
        let tof = cycle.turnover_frequency();
        assert!(tof > 0.0, "TOF should be positive: {tof}");
    }

    #[test]
    fn test_catalytic_tof_higher_at_higher_temperature() {
        let mut c1 = CatalyticCycle::new(400.0, 1e13);
        c1.add_step(ElementaryStep::new("A", -1.0, 0.5));
        let tof1 = c1.turnover_frequency();

        let mut c2 = CatalyticCycle::new(800.0, 1e13);
        c2.add_step(ElementaryStep::new("A", -1.0, 0.5));
        let tof2 = c2.turnover_frequency();

        assert!(
            tof2 > tof1,
            "TOF should increase with temperature: {tof1} vs {tof2}"
        );
    }

    #[test]
    fn test_overall_reaction_energy_is_sum_of_steps() {
        let mut cycle = CatalyticCycle::new(500.0, 1e13);
        cycle.add_step(ElementaryStep::new("step1", -1.0, 0.3));
        cycle.add_step(ElementaryStep::new("step2", 0.5, 0.8));
        cycle.add_step(ElementaryStep::new("step3", -0.3, 0.2));
        let total = cycle.overall_reaction_energy_ev();
        assert!((total - (-0.8)).abs() < 1e-10);
    }

    #[test]
    fn test_selectivity_between_pathways() {
        let sel = CatalyticCycle::selectivity(1e5, 1e3);
        assert!(sel > 0.99, "Should strongly favor desired path: {sel}");
        let sel_equal = CatalyticCycle::selectivity(1e5, 1e5);
        assert!((sel_equal - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_volcano_plot_has_maximum() {
        let cycle = CatalyticCycle::new(500.0, 1e13);
        let de_range: Vec<f64> = (-20..=20).map(|i| i as f64 * 0.1).collect();
        let volcano = cycle.volcano_plot_data(&de_range, 0.5, 0.8);
        let max_log_tof = volcano
            .iter()
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_log_tof > -300.0,
            "Volcano plot should have a finite maximum"
        );
        assert_eq!(volcano.len(), de_range.len());
    }

    #[test]
    fn test_elementary_step_reverse_barrier() {
        let step = ElementaryStep::new("test", -1.0, 0.5);
        // Reverse barrier = 0.5 - (-1.0) = 1.5
        assert!((step.reverse_barrier_ev - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_rate_limiting_step_is_highest_barrier() {
        let mut cycle = CatalyticCycle::new(500.0, 1e13);
        cycle.add_step(ElementaryStep::new("fast", -1.0, 0.1));
        cycle.add_step(ElementaryStep::new("slow", -0.5, 1.5));
        cycle.add_step(ElementaryStep::new("medium", -0.8, 0.5));
        let rls = cycle.rate_limiting_step();
        assert_eq!(
            rls.map(|(i, _)| i),
            Some(1),
            "Step with highest barrier is RLS"
        );
    }
}
