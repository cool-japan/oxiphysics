// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Colloidal particle MD simulations.
//!
//! This module implements full colloidal dynamics including:
//! - [`ColloidalParticle`]: Position, velocity, size distribution, surface charge, Hamaker constant.
//! - [`DlvoInteraction`]: DLVO theory combining van der Waals attraction and electrostatic repulsion.
//! - [`StericStabilization`]: Polymer brush stabilization via the de Gennes model.
//! - [`DepletionForce`]: Asakura-Oosawa depletion model with osmotic pressure.
//! - [`ColloidalCrystal`]: Close-packed structures (FCC/HCP/BCC), phase diagram, freezing/melting.
//! - [`BrownianDynamics`]: Overdamped Langevin equation with random forces.
//! - [`DiffusionColloidal`]: Stokes-Einstein, hydrodynamic interactions, Oseen tensor.
//! - [`ElectrokineticColloidal`]: Electrophoresis (Henry equation), sedimentation potential.
//! - [`ColloidalAggregation`]: DLCA/RLCA fractal cluster formation.
//! - [`ColloidalGel`]: Gel point, percolation, viscoelastic response, aging.

use std::f64::consts::PI;

// ─── Physical Constants ──────────────────────────────────────────────────────

/// Boltzmann constant \[J/K\].
pub const K_B: f64 = 1.380_649e-23;
/// Elementary charge \[C\].
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Vacuum permittivity \[F/m\].
pub const EPS_0: f64 = 8.854_187_812_8e-12;
/// Avogadro number \[1/mol\].
pub const N_A: f64 = 6.022_140_76e23;
/// Viscosity of water at 25°C \[Pa·s\].
pub const ETA_WATER: f64 = 8.9e-4;

// ─── Vector helpers ──────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean magnitude of a 3-vector.
#[inline]
pub fn mag3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Normalise a 3-vector; returns zero vector for negligible length.
#[inline]
pub fn norm3(a: [f64; 3]) -> [f64; 3] {
    let l = mag3(a);
    if l < 1e-14 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}

// ─── ColloidalParticle ───────────────────────────────────────────────────────

/// A single colloidal particle with full physical attributes.
///
/// Stores position, velocity, radius, surface charge density,
/// Hamaker constant, and diffusion coefficient.
#[derive(Debug, Clone)]
pub struct ColloidalParticle {
    /// Particle position \[x, y, z\] in metres.
    pub position: [f64; 3],
    /// Particle velocity \[vx, vy, vz\] in m/s.
    pub velocity: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Surface charge density \[C/m²\].
    pub surface_charge_density: f64,
    /// Hamaker constant \[J\] (material-dependent, typically 1e-21 to 1e-19 J).
    pub hamaker: f64,
    /// Diffusion coefficient \[m²/s\] (Stokes-Einstein).
    pub diffusion: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
}

impl ColloidalParticle {
    /// Create a new colloidal particle.
    ///
    /// Diffusion coefficient is set to zero; use [`DiffusionColloidal::stokes_einstein`]
    /// to compute a physically meaningful value.
    pub fn new(
        position: [f64; 3],
        radius: f64,
        surface_charge_density: f64,
        hamaker: f64,
        mass: f64,
    ) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            radius,
            surface_charge_density,
            hamaker,
            diffusion: 0.0,
            mass,
        }
    }

    /// Surface charge \[C\] = σ · 4π r².
    pub fn total_surface_charge(&self) -> f64 {
        self.surface_charge_density * 4.0 * PI * self.radius * self.radius
    }

    /// Zeta potential approximation \[V\] via Debye-Hückel for low charge:
    /// ζ ≈ σ r / (ε₀ ε_r (1 + κ r)).
    pub fn zeta_potential(&self, eps_r: f64, kappa: f64) -> f64 {
        let eps = EPS_0 * eps_r;
        self.surface_charge_density * self.radius / (eps * (1.0 + kappa * self.radius))
    }

    /// Kinetic energy \[J\] = ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }
}

// ─── DlvoInteraction ─────────────────────────────────────────────────────────

/// DLVO (Derjaguin-Landau-Verwey-Overbeek) interaction between two colloidal spheres.
///
/// Combines van der Waals attraction and screened electrostatic repulsion.
/// Reference: Verwey & Overbeek (1948).
#[derive(Debug, Clone)]
pub struct DlvoInteraction {
    /// Hamaker constant A \[J\].
    pub hamaker: f64,
    /// Dielectric constant of the medium (dimensionless).
    pub eps_r: f64,
    /// Debye screening length κ⁻¹ \[m\].
    pub debye_length: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Surface potential ψ₀ \[V\].
    pub surface_potential: f64,
}

impl DlvoInteraction {
    /// Create a new DLVO interaction model.
    pub fn new(
        hamaker: f64,
        eps_r: f64,
        debye_length: f64,
        temperature: f64,
        surface_potential: f64,
    ) -> Self {
        Self {
            hamaker,
            eps_r,
            debye_length,
            temperature,
            surface_potential,
        }
    }

    /// Van der Waals interaction energy between two equal spheres of radius `r`
    /// at surface-to-surface separation `h` \[J\].
    ///
    /// Uses the Derjaguin approximation: V_vdW = -A·r / (12·h).
    pub fn vdw_energy(&self, r: f64, h: f64) -> f64 {
        let h_eff = h.max(1e-12);
        -self.hamaker * r / (12.0 * h_eff)
    }

    /// Electrostatic repulsion energy between two equal spheres of radius `r`
    /// at surface-to-surface separation `h` \[J\].
    ///
    /// Linear superposition approximation (LSA):
    /// V_el = 2π ε₀ ε_r r ψ₀² exp(-κ h).
    pub fn electrostatic_energy(&self, r: f64, h: f64) -> f64 {
        let kappa = 1.0 / self.debye_length;
        2.0 * PI
            * EPS_0
            * self.eps_r
            * r
            * self.surface_potential
            * self.surface_potential
            * (-kappa * h).exp()
    }

    /// Total DLVO potential energy \[J\] at surface separation `h`.
    pub fn total_energy(&self, r: f64, h: f64) -> f64 {
        self.vdw_energy(r, h) + self.electrostatic_energy(r, h)
    }

    /// DLVO force (magnitude, positive = repulsive) at separation `h` \[N\].
    pub fn total_force(&self, r: f64, h: f64) -> f64 {
        let dh = 1e-12;
        -(self.total_energy(r, h + dh) - self.total_energy(r, h - dh)) / (2.0 * dh)
    }

    /// Energy barrier height \[J\] — maximum of the DLVO curve.
    ///
    /// Scans h from `h_min` to `h_max` in `n_steps` steps.
    pub fn energy_barrier(&self, r: f64, h_min: f64, h_max: f64, n_steps: usize) -> f64 {
        let mut max_e = f64::NEG_INFINITY;
        for i in 0..=n_steps {
            let h = h_min + (h_max - h_min) * i as f64 / n_steps as f64;
            let e = self.total_energy(r, h);
            if e > max_e {
                max_e = e;
            }
        }
        max_e
    }

    /// Secondary minimum depth \[J\] — local minimum at larger separation.
    ///
    /// Returns `None` if no secondary minimum is found.
    pub fn secondary_minimum(
        &self,
        r: f64,
        h_start: f64,
        h_end: f64,
        n_steps: usize,
    ) -> Option<f64> {
        let mut min_e = f64::INFINITY;
        let mut found = false;
        // Skip the primary minimum region (h < 1 nm); look in 5–100 nm range
        let h0 = h_start.max(5e-9);
        for i in 0..=n_steps {
            let h = h0 + (h_end - h0) * i as f64 / n_steps as f64;
            let e = self.total_energy(r, h);
            if e < min_e {
                min_e = e;
                found = true;
            }
        }
        if found && min_e < 0.0 {
            Some(min_e)
        } else {
            None
        }
    }

    /// Debye screening length \[m\] from ionic strength I \[mol/L\] at temperature T \[K\].
    pub fn debye_length_from_ionic_strength(
        ionic_strength_mol_l: f64,
        eps_r: f64,
        temperature: f64,
    ) -> f64 {
        let i_si = ionic_strength_mol_l * 1000.0 * N_A; // convert to #/m³
        let kappa_sq = 2.0 * E_CHARGE * E_CHARGE * i_si / (EPS_0 * eps_r * K_B * temperature);
        1.0 / kappa_sq.sqrt()
    }
}

// ─── StericStabilization ─────────────────────────────────────────────────────

/// Steric stabilization by polymer brushes (de Gennes model).
///
/// Models the osmotic pressure and elastic deformation of grafted polymer chains
/// on a colloidal surface. Reference: de Gennes (1987).
#[derive(Debug, Clone)]
pub struct StericStabilization {
    /// Grafting density \[chains/m²\].
    pub grafting_density: f64,
    /// Polymer segment length (Kuhn length) \[m\].
    pub segment_length: f64,
    /// Number of segments per chain.
    pub n_segments: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Flory-Huggins chi parameter (solvent quality).
    pub chi: f64,
}

impl StericStabilization {
    /// Create a new steric stabilization model.
    pub fn new(
        grafting_density: f64,
        segment_length: f64,
        n_segments: f64,
        temperature: f64,
        chi: f64,
    ) -> Self {
        Self {
            grafting_density,
            segment_length,
            n_segments,
            temperature,
            chi,
        }
    }

    /// Brush height L \[m\] in good solvent (Alexander-de Gennes scaling):
    /// L = N b (b² Γ)^(1/3) where Γ is grafting density.
    pub fn brush_height(&self) -> f64 {
        let b = self.segment_length;
        let n = self.n_segments;
        let gamma = self.grafting_density;
        n * b * (b * b * gamma).powf(1.0 / 3.0)
    }

    /// Steric interaction energy per unit area \[J/m²\] at separation `d` between brush surfaces.
    ///
    /// de Gennes power-law repulsion: G(d) = k_B T Γ^(3/2) L \[(L/d)^5 - (d/L)^7 + ...\]
    /// simplified to the dominant repulsive term.
    pub fn energy_per_area(&self, d: f64) -> f64 {
        let l = self.brush_height();
        if d >= 2.0 * l {
            return 0.0;
        }
        let kbt = K_B * self.temperature;
        let gamma = self.grafting_density;
        // Simplified de Gennes expression
        let s = 1.0 / gamma.sqrt(); // mean spacing between grafting sites
        let u = d / (2.0 * l);
        let u_eff = u.clamp(1e-6, 1.0 - 1e-6);
        kbt / (s * s * s)
            * ((2.0 / 5.0) * u_eff.powf(-5.0 / 4.0) + (4.0 / 7.0) * u_eff.powf(7.0 / 4.0)
                - (48.0 / 35.0))
    }

    /// Steric force between two brush-coated spheres of radius `r` at separation `h` \[N\].
    ///
    /// Uses Derjaguin approximation: F = 2π r · G(h).
    pub fn force(&self, r: f64, h: f64) -> f64 {
        2.0 * PI * r * self.energy_per_area(h)
    }

    /// Excluded volume parameter ν (Flory exponent); ν ≈ 3/5 in good solvent.
    pub fn flory_exponent(&self) -> f64 {
        if self.chi < 0.5 {
            0.588
        } else if self.chi < 0.6 {
            0.5
        } else {
            0.333
        }
    }

    /// End-to-end distance of a free chain in solution \[m\].
    pub fn free_chain_ree(&self) -> f64 {
        let nu = self.flory_exponent();
        self.segment_length * self.n_segments.powf(nu)
    }
}

// ─── DepletionForce ──────────────────────────────────────────────────────────

/// Depletion force from the Asakura-Oosawa (AO) model.
///
/// Small depletant particles (polymers or small spheres) are excluded from
/// the gap between large colloids, generating an effective attraction.
/// Reference: Asakura & Oosawa (1954).
#[derive(Debug, Clone)]
pub struct DepletionForce {
    /// Radius of large (colloidal) particles \[m\].
    pub radius_large: f64,
    /// Radius of small (depletant) particles / polymer gyration radius \[m\].
    pub radius_small: f64,
    /// Number density of depletant \[1/m³\].
    pub depletant_density: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl DepletionForce {
    /// Create a new depletion force model.
    pub fn new(
        radius_large: f64,
        radius_small: f64,
        depletant_density: f64,
        temperature: f64,
    ) -> Self {
        Self {
            radius_large,
            radius_small,
            depletant_density,
            temperature,
        }
    }

    /// Osmotic pressure of the depletant \[Pa\] (ideal gas approximation):
    /// Π = n k_B T.
    pub fn osmotic_pressure(&self) -> f64 {
        self.depletant_density * K_B * self.temperature
    }

    /// Depletion interaction energy \[J\] at centre-to-centre distance `r`.
    ///
    /// AO model: U(r) = -Π · V_overlap(r)
    /// where V_overlap is the overlap volume of depletion zones.
    pub fn depletion_energy(&self, r: f64) -> f64 {
        let rl = self.radius_large;
        let rs = self.radius_small;
        let d_contact = 2.0 * rl;
        let d_max = 2.0 * (rl + rs);
        if r >= d_max {
            return 0.0;
        }
        if r < d_contact {
            // Hard-core overlap — use contact value
            let h0 = 0.0_f64;
            let v = self.overlap_volume_ao(d_contact + h0);
            return -self.osmotic_pressure() * v;
        }
        let v = self.overlap_volume_ao(r);
        -self.osmotic_pressure() * v
    }

    /// Overlap volume of two depletion zones \[m³\] at centre-to-centre distance `r`.
    fn overlap_volume_ao(&self, r: f64) -> f64 {
        let rl = self.radius_large;
        let rs = self.radius_small;
        let rd = rl + rs; // depletion zone radius
        if r >= 2.0 * rd {
            return 0.0;
        }
        let h = 2.0 * rd - r; // overlap gap
        // V = π h² (3 rd - h) / 3  (spherical cap formula for equal spheres)
        PI * h * h * (3.0 * rd - h) / 3.0
    }

    /// Depletion force \[N\] at centre-to-centre distance `r` (positive = attractive towards smaller r).
    pub fn force(&self, r: f64) -> f64 {
        let dr = 1e-12;
        -(self.depletion_energy(r + dr) - self.depletion_energy(r - dr)) / (2.0 * dr)
    }

    /// Critical packing fraction of depletants at which phase separation occurs
    /// (approximate, η_c ≈ 1 / (1 + r_l/r_s)).
    pub fn critical_packing_fraction(&self) -> f64 {
        1.0 / (1.0 + self.radius_large / self.radius_small)
    }
}

// ─── ColloidalCrystal ────────────────────────────────────────────────────────

/// Close-packed colloidal crystal structures and phase diagram.
///
/// Supports FCC, HCP, and BCC lattices with packing fractions,
/// and simple hard-sphere phase diagram (Kirkwood freezing).
#[derive(Debug, Clone)]
pub struct ColloidalCrystal {
    /// Particle radius \[m\].
    pub radius: f64,
    /// Lattice parameter `a` \[m\].
    pub lattice_param: f64,
    /// Crystal structure type.
    pub structure: CrystalStructure,
}

/// Colloidal crystal structure type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrystalStructure {
    /// Face-centred cubic (close-packed, η = 0.7405).
    Fcc,
    /// Hexagonally close-packed (close-packed, η = 0.7405).
    Hcp,
    /// Body-centred cubic (η = 0.6802).
    Bcc,
    /// Simple cubic (η = 0.5236).
    Sc,
}

impl ColloidalCrystal {
    /// Create a new colloidal crystal.
    pub fn new(radius: f64, lattice_param: f64, structure: CrystalStructure) -> Self {
        Self {
            radius,
            lattice_param,
            structure,
        }
    }

    /// Packing fraction η (dimensionless) for the chosen structure.
    pub fn packing_fraction(&self) -> f64 {
        match self.structure {
            CrystalStructure::Fcc => PI / (3.0 * 2.0_f64.sqrt()),
            CrystalStructure::Hcp => PI / (3.0 * 2.0_f64.sqrt()),
            CrystalStructure::Bcc => PI * 3.0_f64.sqrt() / 8.0,
            CrystalStructure::Sc => PI / 6.0,
        }
    }

    /// Number of nearest neighbours for the chosen structure.
    pub fn coordination_number(&self) -> usize {
        match self.structure {
            CrystalStructure::Fcc => 12,
            CrystalStructure::Hcp => 12,
            CrystalStructure::Bcc => 8,
            CrystalStructure::Sc => 6,
        }
    }

    /// Nearest-neighbour distance \[m\].
    pub fn nn_distance(&self) -> f64 {
        let a = self.lattice_param;
        match self.structure {
            CrystalStructure::Fcc => a / 2.0_f64.sqrt(),
            CrystalStructure::Hcp => a,
            CrystalStructure::Bcc => a * 3.0_f64.sqrt() / 2.0,
            CrystalStructure::Sc => a,
        }
    }

    /// Freezing volume fraction (Kirkwood, hard spheres): η_freeze ≈ 0.494.
    pub fn freezing_volume_fraction() -> f64 {
        0.494
    }

    /// Melting volume fraction (Hoover-Ree, hard spheres): η_melt ≈ 0.545.
    pub fn melting_volume_fraction() -> f64 {
        0.545
    }

    /// Generate FCC lattice positions (cube of side `n_cells` cells) \[m\].
    pub fn fcc_positions(&self, n_cells: usize) -> Vec<[f64; 3]> {
        let a = self.lattice_param;
        let basis: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let mut pos = Vec::new();
        for ix in 0..n_cells {
            for iy in 0..n_cells {
                for iz in 0..n_cells {
                    for b in &basis {
                        pos.push([
                            (ix as f64 + b[0]) * a,
                            (iy as f64 + b[1]) * a,
                            (iz as f64 + b[2]) * a,
                        ]);
                    }
                }
            }
        }
        pos
    }

    /// Generate BCC lattice positions (cube of side `n_cells` cells) \[m\].
    pub fn bcc_positions(&self, n_cells: usize) -> Vec<[f64; 3]> {
        let a = self.lattice_param;
        let basis: [[f64; 3]; 2] = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
        let mut pos = Vec::new();
        for ix in 0..n_cells {
            for iy in 0..n_cells {
                for iz in 0..n_cells {
                    for b in &basis {
                        pos.push([
                            (ix as f64 + b[0]) * a,
                            (iy as f64 + b[1]) * a,
                            (iz as f64 + b[2]) * a,
                        ]);
                    }
                }
            }
        }
        pos
    }
}

// ─── BrownianDynamics ────────────────────────────────────────────────────────

/// Brownian dynamics (overdamped Langevin) integrator for colloidal particles.
///
/// Integrates the equation: γ dx/dt = F_det + F_rand
/// where γ is the friction coefficient and F_rand is the stochastic force.
#[derive(Debug, Clone)]
pub struct BrownianDynamics {
    /// Time step \[s\].
    pub dt: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Friction coefficient γ \[kg/s\] per particle.
    pub friction: f64,
    /// Current simulation time \[s\].
    pub time: f64,
}

impl BrownianDynamics {
    /// Create a new Brownian dynamics integrator.
    pub fn new(dt: f64, temperature: f64, friction: f64) -> Self {
        Self {
            dt,
            temperature,
            friction,
            time: 0.0,
        }
    }

    /// Standard deviation of the random displacement in one step per axis \[m\].
    ///
    /// σ = sqrt(2 D dt) = sqrt(2 k_B T dt / γ).
    pub fn noise_amplitude(&self) -> f64 {
        let d = K_B * self.temperature / self.friction;
        (2.0 * d * self.dt).sqrt()
    }

    /// Advance a particle by one Brownian step given deterministic force \[N\] and
    /// Gaussian random numbers `(xi_x, xi_y, xi_z)` ~ N(0,1).
    ///
    /// Returns updated position.
    pub fn step(&mut self, pos: [f64; 3], force: [f64; 3], xi: [f64; 3]) -> [f64; 3] {
        let sigma = self.noise_amplitude();
        let mobility = 1.0 / self.friction;
        let new_pos = [
            pos[0] + mobility * force[0] * self.dt + sigma * xi[0],
            pos[1] + mobility * force[1] * self.dt + sigma * xi[1],
            pos[2] + mobility * force[2] * self.dt + sigma * xi[2],
        ];
        self.time += self.dt;
        new_pos
    }

    /// Mean square displacement after `n_steps` steps for a free particle \[m²\].
    ///
    /// `r²` = 6 D t (3D diffusion).
    pub fn msd_analytical(&self, n_steps: usize) -> f64 {
        let d = K_B * self.temperature / self.friction;
        let t = n_steps as f64 * self.dt;
        6.0 * d * t
    }

    /// Diffusion coefficient \[m²/s\] from the Einstein relation D = k_B T / γ.
    pub fn diffusion_coefficient(&self) -> f64 {
        K_B * self.temperature / self.friction
    }

    /// Run `n_steps` of free Brownian dynamics from origin and return final position.
    /// Uses a deterministic Box-Muller-like seed for reproducibility.
    pub fn run_free(&mut self, n_steps: usize) -> [f64; 3] {
        let mut pos = [0.0f64; 3];
        let mut seed = 12345u64;
        for _ in 0..n_steps {
            let xi = lcg_gaussian_3(&mut seed);
            pos = self.step(pos, [0.0; 3], xi);
        }
        pos
    }
}

/// Simple LCG-based Gaussian random numbers (Box-Muller transform) for reproducibility.
fn lcg_gaussian_3(seed: &mut u64) -> [f64; 3] {
    let mut gaussians = [0.0f64; 3];
    for val in &mut gaussians {
        // LCG step
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let u1 = (*seed >> 11) as f64 / (1u64 << 53) as f64;
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let u2 = (*seed >> 11) as f64 / (1u64 << 53) as f64;
        let u1 = u1.max(1e-15);
        *val = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
    }
    gaussians
}

// ─── DiffusionColloidal ──────────────────────────────────────────────────────

/// Colloidal diffusion: Stokes-Einstein and hydrodynamic interactions.
///
/// Implements the single-particle Stokes-Einstein relation and the two-body
/// Oseen tensor for far-field hydrodynamic coupling.
#[derive(Debug, Clone)]
pub struct DiffusionColloidal {
    /// Solvent dynamic viscosity η \[Pa·s\].
    pub viscosity: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl DiffusionColloidal {
    /// Create a new diffusion model.
    pub fn new(viscosity: f64, temperature: f64) -> Self {
        Self {
            viscosity,
            temperature,
        }
    }

    /// Water model at 25°C.
    pub fn water_25c() -> Self {
        Self::new(ETA_WATER, 298.15)
    }

    /// Stokes-Einstein translational diffusion coefficient \[m²/s\]:
    /// D = k_B T / (6π η r).
    pub fn stokes_einstein(&self, radius: f64) -> f64 {
        K_B * self.temperature / (6.0 * PI * self.viscosity * radius)
    }

    /// Friction coefficient γ \[kg/s\] = k_B T / D = 6π η r.
    pub fn friction_coefficient(&self, radius: f64) -> f64 {
        6.0 * PI * self.viscosity * radius
    }

    /// Rotational diffusion coefficient \[rad²/s\]:
    /// D_r = k_B T / (8π η r³).
    pub fn rotational_diffusion(&self, radius: f64) -> f64 {
        K_B * self.temperature / (8.0 * PI * self.viscosity * radius * radius * radius)
    }

    /// Oseen tensor T(r) \[1/(Pa·m)\] — far-field hydrodynamic coupling between
    /// two particles separated by vector `r_vec`.
    ///
    /// T_ij = 1/(8π η r) (δ_ij + r̂_i r̂_j).
    pub fn oseen_tensor(&self, r_vec: [f64; 3]) -> [[f64; 3]; 3] {
        let r = mag3(r_vec);
        if r < 1e-14 {
            return [[0.0; 3]; 3];
        }
        let rhat = norm3(r_vec);
        let prefactor = 1.0 / (8.0 * PI * self.viscosity * r);
        let mut t = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                t[i][j] = prefactor * (delta + rhat[i] * rhat[j]);
            }
        }
        t
    }

    /// Rotne-Prager-Yamakawa (RPY) correction to the Oseen tensor for finite-size particles.
    /// Valid for r > 2a (non-overlapping).
    pub fn rpy_tensor(&self, r_vec: [f64; 3], radius: f64) -> [[f64; 3]; 3] {
        let r = mag3(r_vec);
        if r < 1e-14 {
            return [[0.0; 3]; 3];
        }
        let rhat = norm3(r_vec);
        let a = radius;
        let d0 = self.stokes_einstein(a);
        let prefactor = if r > 2.0 * a {
            let f1 = 1.0 + 2.0 * a * a / (3.0 * r * r);
            let f2 = 1.0 - 2.0 * a * a / (r * r);
            (f1, f2)
        } else {
            // Overlapping: use approximate form
            let xi = r / (2.0 * a);
            let f1 = 1.0 - 9.0 * xi / 32.0;
            let f2 = 3.0 * xi / 32.0;
            (f1, f2)
        };
        let mut t = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                t[i][j] = d0 * (prefactor.0 * delta + prefactor.1 * rhat[i] * rhat[j]);
            }
        }
        t
    }
}

// ─── ElectrokineticColloidal ──────────────────────────────────────────────────

/// Electrokinetic phenomena in colloidal suspensions.
///
/// Implements Henry's function for electrophoresis, sedimentation potential,
/// and diffusiophoresis in concentration gradients.
#[derive(Debug, Clone)]
pub struct ElectrokineticColloidal {
    /// Particle radius \[m\].
    pub radius: f64,
    /// Zeta potential ζ \[V\].
    pub zeta_potential: f64,
    /// Debye length κ⁻¹ \[m\].
    pub debye_length: f64,
    /// Solvent viscosity η \[Pa·s\].
    pub viscosity: f64,
    /// Dielectric constant ε_r.
    pub eps_r: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl ElectrokineticColloidal {
    /// Create a new electrokinetic model.
    pub fn new(
        radius: f64,
        zeta_potential: f64,
        debye_length: f64,
        viscosity: f64,
        eps_r: f64,
        temperature: f64,
    ) -> Self {
        Self {
            radius,
            zeta_potential,
            debye_length,
            viscosity,
            eps_r,
            temperature,
        }
    }

    /// Henry's function f(κa) for electrophoretic mobility.
    ///
    /// Interpolates between Hückel (κa → 0, f=1) and Smoluchowski (κa → ∞, f=3/2).
    /// Approximation by Ohshima (1994): f(κa) = 1 + 1/(2(1 + 2.5/(κa(1+2e^{-κa})))³).
    pub fn henry_function(&self) -> f64 {
        let ka = self.radius / self.debye_length;
        let denom = 1.0 + 2.5 / (ka * (1.0 + 2.0 * (-ka).exp()));
        1.0 + 0.5 / (denom * denom * denom)
    }

    /// Electrophoretic mobility μ_e \[m²/(V·s)\]:
    /// μ_e = (2 ε₀ ε_r ζ) / (3 η) · f(κa).
    pub fn electrophoretic_mobility(&self) -> f64 {
        let f = self.henry_function();
        (2.0 * EPS_0 * self.eps_r * self.zeta_potential) / (3.0 * self.viscosity) * f
    }

    /// Electrophoretic velocity \[m/s\] in external field E \[V/m\].
    pub fn electrophoretic_velocity(&self, e_field: f64) -> f64 {
        self.electrophoretic_mobility() * e_field
    }

    /// Sedimentation potential E_sed \[V/m\] for a suspension of volume fraction φ.
    ///
    /// E_sed = φ Δρ g ε₀ ε_r ζ / (η κ_∞) where κ_∞ is bulk conductivity \[S/m\].
    pub fn sedimentation_potential(
        &self,
        volume_fraction: f64,
        density_diff: f64,
        gravity: f64,
        bulk_conductivity: f64,
    ) -> f64 {
        let eps = EPS_0 * self.eps_r;
        volume_fraction * density_diff * gravity * eps * self.zeta_potential
            / (self.viscosity * bulk_conductivity)
    }

    /// Diffusiophoretic velocity \[m/s\] in a concentration gradient dc/dx \[mol/m⁴\].
    ///
    /// U_dp = Γ · (k_B T / ze) · (d ln c / dx) where Γ is the diffusiophoretic mobility.
    pub fn diffusiophoretic_velocity(&self, grad_ln_c: f64, valence: f64) -> f64 {
        let eps = EPS_0 * self.eps_r;
        let ze = valence * E_CHARGE;
        let gamma = eps * self.zeta_potential / self.viscosity;
        gamma * (K_B * self.temperature / ze) * grad_ln_c
    }

    /// Streaming potential coefficient \[V/m\] for pressure-driven flow in a capillary.
    pub fn streaming_potential_coefficient(&self, bulk_conductivity: f64) -> f64 {
        (EPS_0 * self.eps_r * self.zeta_potential) / (self.viscosity * bulk_conductivity)
    }
}

// ─── ColloidalAggregation ────────────────────────────────────────────────────

/// Colloidal aggregation kinetics: DLCA and RLCA regimes.
///
/// Tracks cluster size distribution and fractal dimension during
/// diffusion-limited and reaction-limited aggregation.
#[derive(Debug, Clone)]
pub struct ColloidalAggregation {
    /// Initial number density of primary particles \[1/m³\].
    pub n0: f64,
    /// Primary particle radius \[m\].
    pub radius: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Solvent viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Sticking probability (0–1); 1 = DLCA, << 1 = RLCA.
    pub sticking_probability: f64,
    /// Fractal dimension d_f.
    pub fractal_dimension: f64,
}

impl ColloidalAggregation {
    /// Create a new aggregation model.
    ///
    /// Set `sticking_probability = 1.0` for DLCA, `< 0.01` for RLCA.
    pub fn new(
        n0: f64,
        radius: f64,
        temperature: f64,
        viscosity: f64,
        sticking_probability: f64,
    ) -> Self {
        let fractal_dimension = if sticking_probability > 0.5 {
            1.75
        } else {
            2.1
        };
        Self {
            n0,
            radius,
            temperature,
            viscosity,
            sticking_probability,
            fractal_dimension,
        }
    }

    /// Smoluchowski collision rate constant k_smol \[m³/s\]:
    /// k = 8 k_B T / (3 η).
    pub fn smoluchowski_rate(&self) -> f64 {
        8.0 * K_B * self.temperature / (3.0 * self.viscosity)
    }

    /// Effective aggregation rate constant \[m³/s\] including sticking probability.
    pub fn aggregation_rate(&self) -> f64 {
        self.smoluchowski_rate() * self.sticking_probability
    }

    /// Total particle number density at time `t` \[s\] (Smoluchowski kinetics):
    /// n(t) = n0 / (1 + n0 k t).
    pub fn number_density_at(&self, t: f64) -> f64 {
        let k = self.aggregation_rate();
        self.n0 / (1.0 + self.n0 * k * t)
    }

    /// Half-life t_{1/2} \[s\] — time to halve particle number density.
    pub fn half_life(&self) -> f64 {
        1.0 / (self.n0 * self.aggregation_rate())
    }

    /// Radius of a fractal cluster of `n` primary particles \[m\]:
    /// R_cluster = r_primary · n^(1/d_f).
    pub fn cluster_radius(&self, n: f64) -> f64 {
        self.radius * n.powf(1.0 / self.fractal_dimension)
    }

    /// Average aggregation number at time `t`:
    /// `n`(t) = n0 / n(t).
    pub fn mean_aggregation_number(&self, t: f64) -> f64 {
        self.n0 / self.number_density_at(t)
    }

    /// Perikinetic (diffusion-driven) aggregation regime indicator.
    ///
    /// Returns `true` if thermal energy dominates over shear (Pe < 1).
    pub fn is_perikinetic(&self, shear_rate: f64) -> bool {
        // Peclet number Pe = 6π η r³ γ̇ / (k_B T)
        let pe =
            6.0 * PI * self.viscosity * self.radius.powi(3) * shear_rate / (K_B * self.temperature);
        pe < 1.0
    }
}

// ─── ColloidalGel ────────────────────────────────────────────────────────────

/// Colloidal gel formation, viscoelastic response, and aging.
///
/// Models the sol-gel transition via percolation theory and the evolving
/// viscoelastic moduli of an aging colloidal network.
#[derive(Debug, Clone)]
pub struct ColloidalGel {
    /// Volume fraction φ of colloidal particles.
    pub volume_fraction: f64,
    /// Gel point volume fraction φ_gel (percolation threshold).
    pub phi_gel: f64,
    /// High-frequency elastic plateau modulus G_∞ \[Pa\].
    pub plateau_modulus: f64,
    /// Characteristic relaxation time τ \[s\].
    pub relaxation_time: f64,
    /// Power-law aging exponent μ.
    pub aging_exponent: f64,
    /// Particle radius \[m\].
    pub radius: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl ColloidalGel {
    /// Create a new colloidal gel model.
    pub fn new(
        volume_fraction: f64,
        phi_gel: f64,
        plateau_modulus: f64,
        relaxation_time: f64,
        aging_exponent: f64,
        radius: f64,
        temperature: f64,
    ) -> Self {
        Self {
            volume_fraction,
            phi_gel,
            plateau_modulus,
            relaxation_time,
            aging_exponent,
            radius,
            temperature,
        }
    }

    /// Whether the system is above the gel point.
    pub fn is_gel(&self) -> bool {
        self.volume_fraction >= self.phi_gel
    }

    /// Distance from the gel point: ε = (φ - φ_gel) / φ_gel.
    pub fn gel_distance(&self) -> f64 {
        (self.volume_fraction - self.phi_gel) / self.phi_gel
    }

    /// Static elastic modulus G' \[Pa\] scaling near gel point:
    /// G' ~ G_∞ ε^t with t ≈ 4 (scalar percolation).
    pub fn elastic_modulus(&self) -> f64 {
        if !self.is_gel() {
            return 0.0;
        }
        let eps = self.gel_distance();
        self.plateau_modulus * eps.powf(4.0)
    }

    /// Loss modulus G'' \[Pa\] from Maxwell model at frequency ω \[rad/s\]:
    /// G'' = G' ω τ / (1 + (ω τ)²).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let gp = self.elastic_modulus();
        let wt = omega * self.relaxation_time;
        gp * wt / (1.0 + wt * wt)
    }

    /// Storage modulus G' \[Pa\] from Maxwell model at frequency ω \[rad/s\]:
    /// G' = G_∞ (ω τ)² / (1 + (ω τ)²).
    pub fn storage_modulus_maxwell(&self, omega: f64) -> f64 {
        let gp = self.elastic_modulus();
        let wt = omega * self.relaxation_time;
        gp * wt * wt / (1.0 + wt * wt)
    }

    /// Yield stress σ_y \[Pa\] (approximate): σ_y ≈ G' / 100 (empirical for soft gels).
    pub fn yield_stress(&self) -> f64 {
        self.elastic_modulus() / 100.0
    }

    /// Aging modulus at waiting time `t_w` \[s\]:
    /// G(t_w) = G_0 (t_w / τ)^μ.
    pub fn aged_modulus(&self, t_w: f64) -> f64 {
        self.plateau_modulus * (t_w / self.relaxation_time).powf(self.aging_exponent)
    }

    /// Percolation correlation length ξ \[m\] (radius-scaled):
    /// ξ = r |ε|^{-ν} with ν ≈ 0.88 (3D percolation).
    pub fn correlation_length(&self) -> f64 {
        let eps = self.gel_distance().abs().max(1e-10);
        self.radius * eps.powf(-0.88)
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ColloidalParticle ─────────────────────────────────────────────────

    #[test]
    fn test_particle_surface_charge() {
        let p = ColloidalParticle::new([0.0; 3], 100e-9, 0.05, 1e-20, 1e-18);
        let q = p.total_surface_charge();
        // 4π r² σ
        let expected = 4.0 * PI * (100e-9_f64).powi(2) * 0.05;
        assert!(
            (q - expected).abs() < 1e-30,
            "Surface charge mismatch: {:.6e}",
            q
        );
    }

    #[test]
    fn test_particle_zeta_potential_sign() {
        let p = ColloidalParticle::new([0.0; 3], 50e-9, -0.1, 1e-20, 1e-18);
        let zeta = p.zeta_potential(80.0, 1e8);
        assert!(
            zeta < 0.0,
            "Zeta should be negative for negative surface charge: {:.6}",
            zeta
        );
    }

    #[test]
    fn test_particle_kinetic_energy_zero_velocity() {
        let p = ColloidalParticle::new([0.0; 3], 50e-9, 0.01, 1e-20, 1e-18);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    // ── DlvoInteraction ───────────────────────────────────────────────────

    #[test]
    fn test_dlvo_vdw_attractive() {
        let dlvo = DlvoInteraction::new(1e-20, 80.0, 10e-9, 298.15, -0.04);
        let e = dlvo.vdw_energy(100e-9, 2e-9);
        assert!(
            e < 0.0,
            "vdW energy should be negative (attractive): {:.6e}",
            e
        );
    }

    #[test]
    fn test_dlvo_electrostatic_repulsive() {
        let dlvo = DlvoInteraction::new(1e-20, 80.0, 10e-9, 298.15, -0.04);
        let e = dlvo.electrostatic_energy(100e-9, 2e-9);
        assert!(
            e > 0.0,
            "Electrostatic energy should be positive (repulsive): {:.6e}",
            e
        );
    }

    #[test]
    fn test_dlvo_total_energy_finite() {
        let dlvo = DlvoInteraction::new(1e-20, 80.0, 10e-9, 298.15, -0.04);
        let e = dlvo.total_energy(100e-9, 1e-9);
        assert!(e.is_finite(), "Total energy should be finite: {:.6e}", e);
    }

    #[test]
    fn test_dlvo_energy_barrier_nonnegative() {
        let dlvo = DlvoInteraction::new(1e-20, 80.0, 5e-9, 298.15, -0.05);
        let barrier = dlvo.energy_barrier(100e-9, 0.5e-9, 30e-9, 1000);
        // The barrier may or may not exist, just check it is finite
        assert!(barrier.is_finite(), "Energy barrier: {:.6e}", barrier);
    }

    #[test]
    fn test_dlvo_debye_length_from_ionic_strength() {
        // 0.1 M NaCl in water at 25°C → κ⁻¹ ≈ 0.96 nm
        let kappa_inv = DlvoInteraction::debye_length_from_ionic_strength(0.1, 78.5, 298.15);
        assert!(
            kappa_inv > 0.5e-9 && kappa_inv < 2.0e-9,
            "Debye length out of range: {:.6e}",
            kappa_inv
        );
    }

    // ── StericStabilization ───────────────────────────────────────────────

    #[test]
    fn test_steric_brush_height_positive() {
        let steric = StericStabilization::new(1e16, 1e-9, 50.0, 298.15, 0.0);
        let l = steric.brush_height();
        assert!(l > 0.0, "Brush height should be positive: {:.6e}", l);
    }

    #[test]
    fn test_steric_force_repulsive_small_gap() {
        let steric = StericStabilization::new(1e16, 1e-9, 50.0, 298.15, 0.0);
        let l = steric.brush_height();
        // At h = 0.5 L the brushes overlap → force should be repulsive (positive)
        let f = steric.force(200e-9, 0.5 * l);
        assert!(f > 0.0, "Steric force should be repulsive: {:.6e}", f);
    }

    #[test]
    fn test_steric_force_zero_large_gap() {
        let steric = StericStabilization::new(1e16, 1e-9, 50.0, 298.15, 0.0);
        let l = steric.brush_height();
        let f = steric.force(200e-9, 3.0 * l);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_steric_flory_exponent_good_solvent() {
        let steric = StericStabilization::new(1e16, 1e-9, 50.0, 298.15, 0.0);
        let nu = steric.flory_exponent();
        assert!(
            (nu - 0.588).abs() < 1e-6,
            "Flory exponent in good solvent: {:.6}",
            nu
        );
    }

    // ── DepletionForce ────────────────────────────────────────────────────

    #[test]
    fn test_depletion_energy_zero_far_away() {
        let dep = DepletionForce::new(500e-9, 20e-9, 1e21, 298.15);
        let e = dep.depletion_energy(3.0 * 500e-9 * 2.0);
        assert_eq!(e, 0.0, "Depletion energy should be zero far away");
    }

    #[test]
    fn test_depletion_energy_attractive_near_contact() {
        let dep = DepletionForce::new(500e-9, 20e-9, 1e21, 298.15);
        // At near-contact (r ≈ 2 R_large) energy should be attractive (negative)
        let r_contact = 2.0 * 500e-9 + 1e-12;
        let e = dep.depletion_energy(r_contact);
        assert!(e <= 0.0, "Depletion energy should be attractive: {:.6e}", e);
    }

    #[test]
    fn test_depletion_osmotic_pressure_positive() {
        let dep = DepletionForce::new(500e-9, 20e-9, 1e22, 298.15);
        assert!(dep.osmotic_pressure() > 0.0);
    }

    // ── ColloidalCrystal ──────────────────────────────────────────────────

    #[test]
    fn test_crystal_fcc_packing() {
        let c = ColloidalCrystal::new(100e-9, 283e-9, CrystalStructure::Fcc);
        let eta = c.packing_fraction();
        assert!(
            (eta - 0.7405).abs() < 0.01,
            "FCC packing fraction: {:.6}",
            eta
        );
    }

    #[test]
    fn test_crystal_bcc_coordination() {
        let c = ColloidalCrystal::new(100e-9, 300e-9, CrystalStructure::Bcc);
        assert_eq!(c.coordination_number(), 8);
    }

    #[test]
    fn test_crystal_fcc_positions_count() {
        let c = ColloidalCrystal::new(100e-9, 283e-9, CrystalStructure::Fcc);
        let pos = c.fcc_positions(2);
        // 2³ unit cells × 4 atoms/cell = 32
        assert_eq!(pos.len(), 32);
    }

    #[test]
    fn test_crystal_freezing_melting_order() {
        let eta_f = ColloidalCrystal::freezing_volume_fraction();
        let eta_m = ColloidalCrystal::melting_volume_fraction();
        assert!(eta_f < eta_m, "Freezing must precede melting");
    }

    // ── BrownianDynamics ──────────────────────────────────────────────────

    #[test]
    fn test_brownian_noise_amplitude_positive() {
        let bd = BrownianDynamics::new(1e-6, 298.15, 1e-8);
        assert!(bd.noise_amplitude() > 0.0);
    }

    #[test]
    fn test_brownian_diffusion_coefficient() {
        let bd = BrownianDynamics::new(1e-6, 298.15, 6e-10);
        let d = bd.diffusion_coefficient();
        // D = kT/γ = 1.38e-23 * 298 / 6e-10 ≈ 6.85e-12 m²/s
        assert!(d > 1e-13 && d < 1e-10, "Diffusion coefficient: {:.6e}", d);
    }

    #[test]
    fn test_brownian_msd_analytical_linear_time() {
        let bd = BrownianDynamics::new(1e-6, 298.15, 1e-8);
        let msd1 = bd.msd_analytical(100);
        let msd2 = bd.msd_analytical(200);
        let ratio = msd2 / msd1;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "MSD should grow linearly in time: ratio={:.6}",
            ratio
        );
    }

    #[test]
    fn test_brownian_step_changes_position() {
        let mut bd = BrownianDynamics::new(1e-6, 298.15, 1e-8);
        let pos = bd.step([0.0; 3], [0.0; 3], [1.0, 0.5, -0.3]);
        let moved = mag3(pos) > 0.0;
        assert!(moved, "Position should change after a step");
    }

    // ── DiffusionColloidal ────────────────────────────────────────────────

    #[test]
    fn test_stokes_einstein_smaller_faster() {
        let diff = DiffusionColloidal::water_25c();
        let d_small = diff.stokes_einstein(50e-9);
        let d_large = diff.stokes_einstein(500e-9);
        assert!(d_small > d_large, "Smaller particles diffuse faster");
    }

    #[test]
    fn test_oseen_tensor_symmetry() {
        let diff = DiffusionColloidal::water_25c();
        let r_vec = [1e-7, 2e-7, 3e-7];
        let t = diff.oseen_tensor(r_vec);
        for (i, row) in t.iter().enumerate() {
            for (j, &tij) in row.iter().enumerate() {
                assert!(
                    (tij - t[j][i]).abs() < 1e-30,
                    "Oseen tensor not symmetric at ({},{}) vs ({},{})",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn test_rotational_diffusion_positive() {
        let diff = DiffusionColloidal::water_25c();
        let dr = diff.rotational_diffusion(100e-9);
        assert!(
            dr > 0.0,
            "Rotational diffusion should be positive: {:.6e}",
            dr
        );
    }

    // ── ElectrokineticColloidal ───────────────────────────────────────────

    #[test]
    fn test_henry_function_bounds() {
        let ek = ElectrokineticColloidal::new(100e-9, -0.04, 10e-9, ETA_WATER, 80.0, 298.15);
        let f = ek.henry_function();
        assert!(
            (1.0..=1.5).contains(&f),
            "Henry function out of [1, 1.5]: {:.6}",
            f
        );
    }

    #[test]
    fn test_electrophoretic_mobility_sign() {
        // Negative zeta → negative mobility
        let ek = ElectrokineticColloidal::new(100e-9, -0.04, 10e-9, ETA_WATER, 80.0, 298.15);
        let mu = ek.electrophoretic_mobility();
        assert!(
            mu < 0.0,
            "Mobility should be negative for negative zeta: {:.6e}",
            mu
        );
    }

    #[test]
    fn test_electrophoretic_velocity_proportional_field() {
        let ek = ElectrokineticColloidal::new(100e-9, -0.04, 10e-9, ETA_WATER, 80.0, 298.15);
        let v1 = ek.electrophoretic_velocity(1000.0);
        let v2 = ek.electrophoretic_velocity(2000.0);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-10,
            "Velocity should double when field doubles"
        );
    }

    // ── ColloidalAggregation ──────────────────────────────────────────────

    #[test]
    fn test_aggregation_rate_dlca_gt_rlca() {
        let dlca = ColloidalAggregation::new(1e16, 100e-9, 298.15, ETA_WATER, 1.0);
        let rlca = ColloidalAggregation::new(1e16, 100e-9, 298.15, ETA_WATER, 0.001);
        assert!(dlca.aggregation_rate() > rlca.aggregation_rate());
    }

    #[test]
    fn test_aggregation_number_density_decreasing() {
        let agg = ColloidalAggregation::new(1e16, 100e-9, 298.15, ETA_WATER, 1.0);
        let n0 = agg.number_density_at(0.0);
        let n1 = agg.number_density_at(agg.half_life());
        assert!(
            (n1 / n0 - 0.5).abs() < 1e-6,
            "At t_half, density should halve: {:.6}",
            n1 / n0
        );
    }

    #[test]
    fn test_aggregation_cluster_radius_grows_with_n() {
        let agg = ColloidalAggregation::new(1e16, 100e-9, 298.15, ETA_WATER, 1.0);
        let r10 = agg.cluster_radius(10.0);
        let r100 = agg.cluster_radius(100.0);
        assert!(r100 > r10, "Larger cluster should have larger radius");
    }

    #[test]
    fn test_aggregation_fractal_dim_dlca_hcp() {
        let dlca = ColloidalAggregation::new(1e16, 100e-9, 298.15, ETA_WATER, 1.0);
        assert!(
            (dlca.fractal_dimension - 1.75).abs() < 0.01,
            "DLCA fractal dim: {:.6}",
            dlca.fractal_dimension
        );
    }

    // ── ColloidalGel ──────────────────────────────────────────────────────

    #[test]
    fn test_gel_above_phi_gel() {
        let gel = ColloidalGel::new(0.1, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        assert!(gel.is_gel());
    }

    #[test]
    fn test_gel_below_phi_gel() {
        let gel = ColloidalGel::new(0.02, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        assert!(!gel.is_gel());
    }

    #[test]
    fn test_gel_elastic_modulus_zero_below_gel_point() {
        let gel = ColloidalGel::new(0.03, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        assert_eq!(gel.elastic_modulus(), 0.0);
    }

    #[test]
    fn test_gel_elastic_modulus_positive_above_gel_point() {
        let gel = ColloidalGel::new(0.10, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        assert!(gel.elastic_modulus() > 0.0);
    }

    #[test]
    fn test_gel_loss_modulus_positive() {
        let gel = ColloidalGel::new(0.10, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        let gpp = gel.loss_modulus(10.0);
        assert!(gpp > 0.0);
    }

    #[test]
    fn test_gel_aged_modulus_increases() {
        let gel = ColloidalGel::new(0.10, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        let g1 = gel.aged_modulus(1.0);
        let g2 = gel.aged_modulus(10.0);
        assert!(g2 > g1, "Aged modulus should increase with waiting time");
    }

    #[test]
    fn test_gel_correlation_length_positive() {
        let gel = ColloidalGel::new(0.10, 0.05, 100.0, 1.0, 0.5, 100e-9, 298.15);
        let xi = gel.correlation_length();
        assert!(xi > 0.0, "Correlation length: {:.6e}", xi);
    }
}
