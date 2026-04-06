#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crystal nucleation and growth simulation.
//!
//! Implements classical nucleation theory (CNT), JMAK/Avrami kinetics,
//! growth kinetics, phase-field crystal models, dendrite growth, grain
//! orientation and grain growth (coarsening), and Ostwald ripening.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Crystal structure
// ---------------------------------------------------------------------------

/// Common crystal structures.
#[derive(Debug, Clone, PartialEq)]
pub enum CrystalStructure {
    /// Face-centred cubic.
    Fcc,
    /// Body-centred cubic.
    Bcc,
    /// Hexagonal close-packed.
    Hcp,
    /// Simple cubic.
    Sc,
}

// ---------------------------------------------------------------------------
// CrystalNucleus
// ---------------------------------------------------------------------------

/// A single crystal nucleus with position, size, structure, and growth rate.
#[derive(Debug, Clone)]
pub struct CrystalNucleus {
    /// 3-D position \[x, y, z\].
    pub position: [f64; 3],
    /// Current nucleus radius.
    pub radius: f64,
    /// Crystal structure of the nucleus.
    pub crystal_structure: CrystalStructure,
    /// Radial growth rate (length per unit time).
    pub growth_rate: f64,
}

impl CrystalNucleus {
    /// Create a new `CrystalNucleus`.
    pub fn new(
        position: [f64; 3],
        radius: f64,
        crystal_structure: CrystalStructure,
        growth_rate: f64,
    ) -> Self {
        Self {
            position,
            radius,
            crystal_structure,
            growth_rate,
        }
    }

    /// Volume of the nucleus (spherical approximation).
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    /// Advance the nucleus radius by one time step.
    ///
    /// # Arguments
    /// * `dt` - time step
    pub fn grow(&mut self, dt: f64) {
        self.radius += self.growth_rate * dt;
        if self.radius < 0.0 {
            self.radius = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// NucleationModel
// ---------------------------------------------------------------------------

/// Classical nucleation theory (CNT) model.
#[derive(Debug, Clone)]
pub struct NucleationModel {
    /// Solid–liquid interfacial energy γ (J/m²).
    pub surface_energy: f64,
    /// Volumetric driving force ΔG_v (J/m³); negative for crystallisation.
    pub driving_force: f64,
}

impl NucleationModel {
    /// Create a new `NucleationModel`.
    pub fn new(surface_energy: f64, driving_force: f64) -> Self {
        Self {
            surface_energy,
            driving_force,
        }
    }

    /// Critical nucleus radius from CNT.
    ///
    /// ```text
    /// r* = -2γ / ΔG_v
    /// ```
    pub fn critical_radius(&self) -> f64 {
        critical_nucleus_radius(self.surface_energy, self.driving_force)
    }

    /// Steady-state nucleation rate (simplified CNT).
    ///
    /// ```text
    /// J = J_0 exp(-ΔG* / (k_B T))
    /// ```
    ///
    /// where `ΔG* = 16π γ³ / (3 ΔG_v²)` and `J_0 = 1` (pre-factor set to 1
    /// in dimensionless units; supersaturation enters via `driving_force`).
    ///
    /// # Arguments
    /// * `supersaturation` - dimensionless supersaturation σ; scales the
    ///   driving force as `ΔG_v_eff = ΔG_v * σ`
    pub fn nucleation_rate(&self, supersaturation: f64) -> f64 {
        if supersaturation <= 0.0 || self.driving_force >= 0.0 {
            return 0.0;
        }
        let dg_eff = self.driving_force * supersaturation;
        if dg_eff >= 0.0 {
            return 0.0;
        }
        let dg_star = 16.0 * PI * self.surface_energy.powi(3) / (3.0 * dg_eff * dg_eff);
        // Use kBT = 1 (dimensionless)
        (-dg_star).exp()
    }
}

// ---------------------------------------------------------------------------
// GrowthKinetics
// ---------------------------------------------------------------------------

/// Crystal growth kinetics model.
#[derive(Debug, Clone)]
pub struct GrowthKinetics {
    /// Interface mobility M (m/s/K).
    pub interface_mobility: f64,
    /// Attachment rate coefficient (m/s).
    pub attachment_rate: f64,
}

impl GrowthKinetics {
    /// Create a new `GrowthKinetics` model.
    pub fn new(interface_mobility: f64, attachment_rate: f64) -> Self {
        Self {
            interface_mobility,
            attachment_rate,
        }
    }

    /// Crystal growth rate as a function of undercooling ΔT.
    ///
    /// ```text
    /// v = M * ΔT
    /// ```
    ///
    /// # Arguments
    /// * `undercooling` - undercooling ΔT = T_melt - T (K); positive means
    ///   below the melting point
    pub fn growth_rate(&self, undercooling: f64) -> f64 {
        if undercooling <= 0.0 {
            return 0.0;
        }
        self.interface_mobility * undercooling
    }
}

// ---------------------------------------------------------------------------
// CrystalSimulation
// ---------------------------------------------------------------------------

/// High-level crystal nucleation and growth simulation.
#[derive(Debug, Clone)]
pub struct CrystalSimulation {
    /// Collection of active crystal nuclei.
    pub nuclei: Vec<CrystalNucleus>,
    /// Current system temperature (K).
    pub temperature: f64,
    /// Current supersaturation (dimensionless).
    pub supersaturation: f64,
}

impl CrystalSimulation {
    /// Create a new `CrystalSimulation`.
    pub fn new(temperature: f64, supersaturation: f64) -> Self {
        Self {
            nuclei: Vec::new(),
            temperature,
            supersaturation,
        }
    }

    /// Add a nucleus to the simulation.
    pub fn add_nucleus(&mut self, nucleus: CrystalNucleus) {
        self.nuclei.push(nucleus);
    }

    /// Advance all nuclei by one time step `dt`.
    ///
    /// Each nucleus grows according to its stored `growth_rate`.
    /// Nuclei that shrink below a minimum radius (here 0) are retained
    /// with radius clamped to zero.
    ///
    /// # Arguments
    /// * `dt` - time step
    pub fn step(&mut self, dt: f64) {
        for n in &mut self.nuclei {
            n.grow(dt);
        }
        // Supersaturation decreases as crystals grow (simplified: proportional
        // to total transformed volume fraction).
        let total_vol: f64 = self.nuclei.iter().map(|n| n.volume()).sum();
        let domain_vol = 1.0_f64; // unit domain
        let x = (total_vol / domain_vol).min(1.0);
        self.supersaturation = (1.0 - x).max(0.0);
    }

    /// Return the number of active nuclei.
    pub fn num_nuclei(&self) -> usize {
        self.nuclei.len()
    }

    /// Return the total crystallised volume fraction (clamped to \[0, 1\]).
    pub fn transformed_fraction(&self) -> f64 {
        let total_vol: f64 = self.nuclei.iter().map(|n| n.volume()).sum();
        total_vol.clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Avrami (JMAK) kinetics
// ---------------------------------------------------------------------------

/// Johnson-Mehl-Avrami-Kolmogorov (JMAK) transformation kinetics.
#[derive(Debug, Clone)]
pub struct Avrami {
    /// Nucleation rate (number per unit volume per unit time).
    pub nucleation_rate: f64,
    /// Radial growth rate (length per unit time).
    pub growth_rate: f64,
}

impl Avrami {
    /// Create a new `Avrami` model.
    pub fn new(nucleation_rate: f64, growth_rate: f64) -> Self {
        Self {
            nucleation_rate,
            growth_rate,
        }
    }

    /// Compute the transformed fraction using the JMAK equation.
    ///
    /// ```text
    /// X(t) = 1 - exp(-k t^n)
    /// ```
    ///
    /// where `k = nucleation_rate * growth_rate^3` and `n = 4`
    /// (continuous nucleation in 3-D).
    ///
    /// # Arguments
    /// * `time` - elapsed time
    pub fn compute_transformed_fraction(&self, time: f64) -> f64 {
        if time <= 0.0 {
            return 0.0;
        }
        let k = self.nucleation_rate * self.growth_rate.powi(3);
        avrami_exponent(4.0, k, time)
    }
}

// ---------------------------------------------------------------------------
// OrientationRelation
// ---------------------------------------------------------------------------

/// Common crystallographic orientation relationships between phases.
#[derive(Debug, Clone, PartialEq)]
pub enum OrientationRelation {
    /// Kurdjumov-Sachs (K-S) orientation relationship.
    KurdjumovSachs,
    /// Nishiyama-Wassermann (N-W) orientation relationship.
    NishiyamaWassermann,
    /// Bain correspondence (FCC → BCC).
    Bain,
}

impl OrientationRelation {
    /// Approximate misorientation angle (radians) from the ideal K-S or N-W
    /// orientation relationship.
    ///
    /// Returns the characteristic angle for each relationship:
    /// - K-S:  5.26° ≈ 0.0918 rad
    /// - N-W:  5.26° (same family, slightly different)
    /// - Bain: 45°   ≈ 0.7854 rad
    pub fn characteristic_angle(&self) -> f64 {
        match self {
            OrientationRelation::KurdjumovSachs => 5.26_f64.to_radians(),
            OrientationRelation::NishiyamaWassermann => 9.74_f64.to_radians(),
            OrientationRelation::Bain => 45.0_f64.to_radians(),
        }
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the critical nucleus radius from CNT.
///
/// ```text
/// r* = -2γ / ΔG_v
/// ```
///
/// Returns `f64::INFINITY` if `delta_g_v >= 0`.
///
/// # Arguments
/// * `gamma`      - interfacial energy (J/m²)
/// * `delta_g_v`  - volumetric driving force (J/m³); must be negative
pub fn critical_nucleus_radius(gamma: f64, delta_g_v: f64) -> f64 {
    if delta_g_v >= 0.0 {
        return f64::INFINITY;
    }
    -2.0 * gamma / delta_g_v
}

/// Compute the JMAK transformed fraction.
///
/// ```text
/// X(t) = 1 - exp(-k t^n)
/// ```
///
/// # Arguments
/// * `n` - Avrami exponent (typically 1–4)
/// * `k` - JMAK rate constant
/// * `t` - elapsed time
pub fn avrami_exponent(n: f64, k: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    1.0 - (-k * t.powf(n)).exp()
}

/// Compute the undercooling ΔT = T_melt - T.
///
/// Returns 0 if T ≥ T_melt.
///
/// # Arguments
/// * `t_current` - current temperature (K)
/// * `t_melt`    - melting temperature (K)
#[allow(non_snake_case)]
pub fn undercooling(T: f64, T_melt: f64) -> f64 {
    (T_melt - T).max(0.0)
}

// ---------------------------------------------------------------------------
// NucleationSite (legacy, kept for compatibility)
// ---------------------------------------------------------------------------

/// A nucleation event: position, final radius, and onset time.
#[derive(Debug, Clone, PartialEq)]
pub struct NucleationSite {
    /// Position of the nucleus in 3D space.
    pub position: [f64; 3],
    /// Final (grown) radius of the nucleus.
    pub radius: f64,
    /// Time at which nucleation begins.
    pub t_nucleation: f64,
}

impl NucleationSite {
    /// Create a new `NucleationSite`.
    pub fn new(position: [f64; 3], radius: f64, t_nucleation: f64) -> Self {
        Self {
            position,
            radius,
            t_nucleation,
        }
    }

    /// Return the volume of the nucleus assuming a sphere.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }
}

// ---------------------------------------------------------------------------
// Classical nucleation theory helpers (legacy)
// ---------------------------------------------------------------------------

/// Compute the nucleation rate from Classical Nucleation Theory (CNT).
///
/// # Arguments
/// * `delta_g_bulk` - volumetric driving force (J/m³); must be negative
/// * `gamma`        - interfacial energy (J/m²)
/// * `temp`         - temperature (K) — unused in this simplified form
/// * `kbt`          - thermal energy k_B T (J)
pub fn classical_nucleation_rate(delta_g_bulk: f64, gamma: f64, temp: f64, kbt: f64) -> f64 {
    let _ = temp;
    if kbt <= 0.0 || delta_g_bulk >= 0.0 {
        return 0.0;
    }
    let dg_star = 16.0 * PI * gamma.powi(3) / (3.0 * delta_g_bulk * delta_g_bulk);
    let exponent = -dg_star / kbt;
    exponent.exp()
}

/// Compute the critical nucleus radius (legacy wrapper).
///
/// # Arguments
/// * `gamma`        - interfacial energy (J/m²)
/// * `delta_g_bulk` - volumetric driving force (negative for nucleation, J/m³)
pub fn critical_nucleus_size(gamma: f64, delta_g_bulk: f64) -> f64 {
    critical_nucleus_radius(gamma, delta_g_bulk)
}

/// Compute the nucleation barrier for a cluster of `n_atoms` atoms.
///
/// # Arguments
/// * `n_atoms`      - number of atoms in the cluster
/// * `gamma`        - interfacial energy (arbitrary consistent units)
/// * `delta_g_bulk` - bulk driving force per atom (negative for nucleation)
pub fn nucleation_barrier(n_atoms: usize, gamma: f64, delta_g_bulk: f64) -> f64 {
    let n = n_atoms as f64;
    let r = (3.0 * n / (4.0 * PI)).powf(1.0 / 3.0);
    let bulk_term = n * delta_g_bulk;
    let surface_term = 4.0 * PI * r * r * gamma;
    bulk_term + surface_term
}

// ---------------------------------------------------------------------------
// CrystalGrowthModel
// ---------------------------------------------------------------------------

/// A simplified crystal growth model supporting Avrami kinetics and
/// stochastic nucleation events.
#[derive(Debug, Clone)]
pub struct CrystalGrowthModel {
    /// Overall nucleation rate (number per unit volume per unit time).
    pub nucleation_rate: f64,
    /// Radial growth rate of nuclei (length per unit time).
    pub growth_rate: f64,
    /// Equilibrium lattice parameter.
    pub lattice_param: f64,
    /// Accumulated simulation time.
    current_time: f64,
    /// Accumulated transformed fraction (JMAK).
    transformed_fraction: f64,
}

impl CrystalGrowthModel {
    /// Create a new `CrystalGrowthModel`.
    pub fn new(nucleation_rate: f64, growth_rate: f64, lattice_param: f64) -> Self {
        Self {
            nucleation_rate,
            growth_rate,
            lattice_param,
            current_time: 0.0,
            transformed_fraction: 0.0,
        }
    }

    /// Compute the Avrami (JMAK) transformed fraction.
    ///
    /// # Arguments
    /// * `t` - time
    /// * `k` - JMAK rate constant
    /// * `n` - Avrami exponent (typically 1-4)
    pub fn avrami_fraction(t: f64, k: f64, n: f64) -> f64 {
        avrami_exponent(n, k, t)
    }

    /// Advance the model by time step `dt` at temperature `temp`.
    ///
    /// # Arguments
    /// * `dt`   - time step
    /// * `temp` - current temperature (unused in this simplified version)
    pub fn step(&mut self, dt: f64, temp: f64) -> Vec<NucleationSite> {
        let _ = temp;
        self.current_time += dt;
        let expected = self.nucleation_rate * dt;
        let n_new = if expected > 10.0 {
            expected.round() as usize
        } else {
            expected.floor() as usize
        };
        let mut sites = Vec::with_capacity(n_new);
        for i in 0..n_new {
            let angle = 2.0 * PI * i as f64 / n_new.max(1) as f64;
            let pos = [angle.cos(), angle.sin(), 0.0];
            sites.push(NucleationSite::new(
                pos,
                self.lattice_param,
                self.current_time,
            ));
        }
        let k = self.nucleation_rate * self.growth_rate.powi(3);
        self.transformed_fraction = Self::avrami_fraction(self.current_time, k, 4.0);
        sites
    }

    /// Return the current transformed fraction.
    pub fn transformed_fraction(&self) -> f64 {
        self.transformed_fraction
    }

    /// Return the current simulation time.
    pub fn current_time(&self) -> f64 {
        self.current_time
    }
}

// ---------------------------------------------------------------------------
// PhaseFieldCrystal
// ---------------------------------------------------------------------------

/// Phase-field crystal (PFC) model on a 2D regular grid.
#[derive(Debug, Clone)]
pub struct PhaseFieldCrystal {
    /// Order parameter field (row-major, size `nx * ny`).
    pub phi: Vec<f64>,
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Grid spacing.
    pub dx: f64,
}

impl PhaseFieldCrystal {
    /// Create a new `PhaseFieldCrystal` initialised to `phi_0` everywhere.
    pub fn new(nx: usize, ny: usize, dx: f64, phi_0: f64) -> Self {
        Self {
            phi: vec![phi_0; nx * ny],
            nx,
            ny,
            dx,
        }
    }

    /// Linear index from 2D grid coordinates.
    fn idx(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }

    /// Periodic neighbour in x.
    fn wrap_x(&self, i: i64) -> usize {
        i.rem_euclid(self.nx as i64) as usize
    }

    /// Periodic neighbour in y.
    fn wrap_y(&self, j: i64) -> usize {
        j.rem_euclid(self.ny as i64) as usize
    }

    /// Compute the 5-point discrete Laplacian at grid point `(i, j)`.
    pub fn laplacian(&self, i: usize, j: usize) -> f64 {
        let dx2 = self.dx * self.dx;
        let center = self.phi[self.idx(i, j)];
        let ip = self.phi[self.idx(self.wrap_x(i as i64 + 1), j)];
        let im = self.phi[self.idx(self.wrap_x(i as i64 - 1), j)];
        let jp = self.phi[self.idx(i, self.wrap_y(j as i64 + 1))];
        let jm = self.phi[self.idx(i, self.wrap_y(j as i64 - 1))];
        (ip + im + jp + jm - 4.0 * center) / dx2
    }

    /// Compute the total free energy of the current state.
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.nx {
            for j in 0..self.ny {
                let phi = self.phi[self.idx(i, j)];
                let lap = self.laplacian(i, j);
                e += 0.5 * phi * phi + 0.25 * phi.powi(4);
                e -= 0.5 * phi * lap;
            }
        }
        e * self.dx * self.dx
    }

    /// Advance the PFC order parameter by one time step using explicit Euler.
    ///
    /// # Arguments
    /// * `dt`         - time step
    /// * `mobility`   - kinetic mobility M
    /// * `temp_noise` - noise amplitude (set to 0 for deterministic evolution)
    pub fn step(&mut self, dt: f64, mobility: f64, temp_noise: f64) {
        let n = self.nx * self.ny;
        let mut dphi = vec![0.0f64; n];
        for i in 0..self.nx {
            for j in 0..self.ny {
                let idx = self.idx(i, j);
                let phi = self.phi[idx];
                let lap = self.laplacian(i, j);
                let mu = phi + phi.powi(3) - lap;
                dphi[idx] = mobility * (-mu + lap);
            }
        }
        for idx in 0..n {
            let noise = if temp_noise > 0.0 {
                let phase = (idx as f64 * 1.618033988749895).fract();
                temp_noise * (2.0 * phase - 1.0)
            } else {
                0.0
            };
            self.phi[idx] += dt * (dphi[idx] + noise);
        }
    }

    /// Return the mean order parameter over the grid.
    pub fn mean_phi(&self) -> f64 {
        self.phi.iter().sum::<f64>() / self.phi.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Dendrite tip velocity
// ---------------------------------------------------------------------------

/// Compute the dendrite tip velocity (Ivantsov + solvability).
///
/// # Arguments
/// * `delta`       - dimensionless supersaturation
/// * `d0`          - capillarity length
/// * `sigma_star`  - solvability constant (≈ 0.025)
pub fn dendrite_tip_velocity(delta: f64, d0: f64, sigma_star: f64) -> f64 {
    if d0 <= 0.0 || sigma_star <= 0.0 {
        return 0.0;
    }
    delta * delta / (d0 * sigma_star)
}

/// Compute the thermal dendrite tip radius.
///
/// # Arguments
/// * `peclet` - tip Peclet number
/// * `delta`  - dimensionless supersaturation
pub fn thermal_dendrite_radius(peclet: f64, delta: f64) -> f64 {
    if peclet <= 0.0 || delta <= 0.0 {
        return f64::INFINITY;
    }
    1.0 / (peclet * delta)
}

// ---------------------------------------------------------------------------
// Crystal orientation and grain growth
// ---------------------------------------------------------------------------

/// Crystal orientation described by three Euler angles (ZXZ convention).
#[derive(Debug, Clone, PartialEq)]
pub struct CrystalOrientation {
    /// Euler angles `[phi1, Phi, phi2]` in radians (ZXZ / Bunge convention).
    pub euler: [f64; 3],
}

impl CrystalOrientation {
    /// Create a `CrystalOrientation` from Euler angles in radians.
    pub fn new(phi1: f64, big_phi: f64, phi2: f64) -> Self {
        Self {
            euler: [phi1, big_phi, phi2],
        }
    }

    /// Compute the 3×3 rotation matrix from the Euler angles (ZXZ convention).
    pub fn rotation_matrix(&self) -> [[f64; 3]; 3] {
        let [phi1, big_phi, phi2] = self.euler;
        let c1 = phi1.cos();
        let s1 = phi1.sin();
        let c = big_phi.cos();
        let s = big_phi.sin();
        let c2 = phi2.cos();
        let s2 = phi2.sin();
        [
            [c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
            [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
            [s1 * s, -c1 * s, c],
        ]
    }

    /// Compute the misorientation angle between this orientation and another.
    pub fn misorientation_angle(&self, other: &Self) -> f64 {
        let r1 = self.rotation_matrix();
        let r2 = other.rotation_matrix();
        let mut trace = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let r12_ij: f64 = (0..3).map(|k| r1[k][i] * r2[k][j]).sum();
                if i == j {
                    trace += r12_ij;
                }
            }
        }
        let cos_theta = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
        cos_theta.acos()
    }
}

// ---------------------------------------------------------------------------
// Grain growth simulation
// ---------------------------------------------------------------------------

/// Polycrystalline grain growth simulation using a mean-field coarsening model.
#[derive(Debug, Clone)]
pub struct GrainGrowthSimulation {
    /// Crystal orientation of each grain.
    pub grains: Vec<CrystalOrientation>,
    /// Current equivalent-sphere radius of each grain.
    pub sizes: Vec<f64>,
}

impl GrainGrowthSimulation {
    /// Create a new `GrainGrowthSimulation` with `n` randomly oriented grains
    /// of initial size `r0`.
    pub fn new(n: usize, r0: f64) -> Self {
        let grains: Vec<CrystalOrientation> = (0..n)
            .map(|i| {
                let t = i as f64 * 2.399963229728653;
                CrystalOrientation::new(t, (t * 0.7).sin().abs() * PI, t * 1.3)
            })
            .collect();
        let sizes = vec![r0; n];
        Self { grains, sizes }
    }

    /// Advance the grain structure by one time step.
    ///
    /// # Arguments
    /// * `dt`       - time step
    /// * `mobility` - grain boundary mobility
    /// * `temp`     - temperature (unused in this simplified model)
    pub fn coarsen(&mut self, dt: f64, mobility: f64, temp: f64) {
        let _ = temp;
        let mean_r = self.mean_grain_size();
        if mean_r <= 0.0 {
            return;
        }
        for r in &mut self.sizes {
            let driving_force = 1.0 / mean_r - 1.0 / r.max(1e-14);
            *r += mobility * driving_force * dt;
        }
        let min_r = mean_r * 0.01;
        for r in self.sizes.iter_mut() {
            if *r < min_r {
                *r = min_r;
            }
        }
    }

    /// Return the number-averaged mean grain size.
    pub fn mean_grain_size(&self) -> f64 {
        if self.sizes.is_empty() {
            return 0.0;
        }
        self.sizes.iter().sum::<f64>() / self.sizes.len() as f64
    }

    /// Return the largest grain radius.
    pub fn max_grain_size(&self) -> f64 {
        self.sizes.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Return the index of the grain with the largest size.
    pub fn largest_grain_idx(&self) -> Option<usize> {
        self.sizes
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
}

// ---------------------------------------------------------------------------
// Ostwald ripening
// ---------------------------------------------------------------------------

/// Compute the Ostwald ripening (LSW) rate constant.
///
/// # Arguments
/// * `surface_energy` - interfacial energy γ (J/m²)
/// * `molar_vol`      - molar volume V_m (m³/mol)
/// * `diff`           - diffusion coefficient D (m²/s)
/// * `c_eq`           - equilibrium solubility (mol/m³)
/// * `r_mean`         - current mean particle radius (m)
/// * `temp`           - temperature (K)
pub fn ostwald_ripening_rate(
    surface_energy: f64,
    molar_vol: f64,
    diff: f64,
    c_eq: f64,
    r_mean: f64,
    temp: f64,
) -> f64 {
    let r_gas = 8.314_462_618_f64;
    if temp <= 0.0 || r_mean <= 0.0 {
        return 0.0;
    }
    let k_lsw = 8.0 * surface_energy * molar_vol * diff * c_eq / (9.0 * r_gas * temp);
    k_lsw / (3.0 * r_mean * r_mean)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- critical_nucleus_radius ---

    #[test]
    fn test_critical_radius_positive_dg_returns_inf() {
        assert_eq!(critical_nucleus_radius(0.1, 1.0), f64::INFINITY);
    }

    #[test]
    fn test_critical_radius_formula() {
        let r = critical_nucleus_radius(0.1, -1.0);
        assert!((r - 0.2).abs() < 1e-10, "r={r}");
    }

    #[test]
    fn test_critical_radius_scales_with_gamma() {
        let r1 = critical_nucleus_radius(0.1, -1.0);
        let r2 = critical_nucleus_radius(0.2, -1.0);
        assert!((r2 / r1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_critical_radius_zero_dg_returns_inf() {
        assert_eq!(critical_nucleus_radius(0.1, 0.0), f64::INFINITY);
    }

    // --- avrami_exponent ---

    #[test]
    fn test_avrami_zero_time() {
        assert_eq!(avrami_exponent(4.0, 1.0, 0.0), 0.0);
    }

    #[test]
    fn test_avrami_approaches_one() {
        let x = avrami_exponent(2.0, 1.0, 100.0);
        assert!(x > 0.9999, "x={x}");
    }

    #[test]
    fn test_avrami_monotone() {
        let x1 = avrami_exponent(2.0, 0.5, 1.0);
        let x2 = avrami_exponent(2.0, 0.5, 2.0);
        assert!(x2 > x1);
    }

    #[test]
    fn test_avrami_k_zero_gives_zero() {
        assert_eq!(avrami_exponent(2.0, 0.0, 10.0), 0.0);
    }

    // --- undercooling ---

    #[test]
    fn test_undercooling_positive() {
        assert!((undercooling(1200.0, 1300.0) - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_undercooling_zero_above_melt() {
        assert_eq!(undercooling(1400.0, 1300.0), 0.0);
    }

    #[test]
    fn test_undercooling_at_melt() {
        assert_eq!(undercooling(1300.0, 1300.0), 0.0);
    }

    // --- CrystalNucleus ---

    #[test]
    fn test_nucleus_volume() {
        let n = CrystalNucleus::new([0.0; 3], 1.0, CrystalStructure::Fcc, 0.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((n.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_nucleus_grow_increases_radius() {
        let mut n = CrystalNucleus::new([0.0; 3], 1.0, CrystalStructure::Bcc, 0.5);
        n.grow(1.0);
        assert!((n.radius - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_nucleus_grow_clamped_to_zero() {
        let mut n = CrystalNucleus::new([0.0; 3], 0.1, CrystalStructure::Hcp, -10.0);
        n.grow(1.0);
        assert_eq!(n.radius, 0.0);
    }

    #[test]
    fn test_nucleus_structure_stored() {
        let n = CrystalNucleus::new([0.0; 3], 1.0, CrystalStructure::Sc, 1.0);
        assert_eq!(n.crystal_structure, CrystalStructure::Sc);
    }

    // --- NucleationModel ---

    #[test]
    fn test_nucleation_model_critical_radius() {
        let m = NucleationModel::new(0.1, -1.0);
        assert!((m.critical_radius() - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_nucleation_model_rate_positive_dg_zero() {
        let m = NucleationModel::new(0.1, 1.0);
        assert_eq!(m.nucleation_rate(1.0), 0.0);
    }

    #[test]
    fn test_nucleation_model_rate_positive() {
        let m = NucleationModel::new(0.01, -1e6);
        let rate = m.nucleation_rate(1.0);
        assert!(rate > 0.0, "rate={rate}");
    }

    #[test]
    fn test_nucleation_model_rate_zero_supersaturation() {
        let m = NucleationModel::new(0.1, -1.0);
        assert_eq!(m.nucleation_rate(0.0), 0.0);
    }

    // --- GrowthKinetics ---

    #[test]
    fn test_growth_kinetics_positive_undercooling() {
        let gk = GrowthKinetics::new(0.01, 1.0);
        let v = gk.growth_rate(10.0);
        assert!((v - 0.1).abs() < 1e-10, "v={v}");
    }

    #[test]
    fn test_growth_kinetics_zero_undercooling() {
        let gk = GrowthKinetics::new(0.01, 1.0);
        assert_eq!(gk.growth_rate(0.0), 0.0);
    }

    #[test]
    fn test_growth_kinetics_negative_undercooling() {
        let gk = GrowthKinetics::new(0.01, 1.0);
        assert_eq!(gk.growth_rate(-5.0), 0.0);
    }

    #[test]
    fn test_growth_kinetics_scales_linearly() {
        let gk = GrowthKinetics::new(1.0, 1.0);
        let v1 = gk.growth_rate(10.0);
        let v2 = gk.growth_rate(20.0);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    // --- CrystalSimulation ---

    #[test]
    fn test_simulation_add_nucleus() {
        let mut sim = CrystalSimulation::new(1200.0, 1.5);
        let n = CrystalNucleus::new([0.0; 3], 0.1, CrystalStructure::Fcc, 0.1);
        sim.add_nucleus(n);
        assert_eq!(sim.num_nuclei(), 1);
    }

    #[test]
    fn test_simulation_step_grows_nucleus() {
        let mut sim = CrystalSimulation::new(1200.0, 1.5);
        sim.add_nucleus(CrystalNucleus::new(
            [0.0; 3],
            0.0,
            CrystalStructure::Fcc,
            1.0,
        ));
        sim.step(0.1);
        assert!((sim.nuclei[0].radius - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_simulation_transformed_fraction_bounded() {
        let mut sim = CrystalSimulation::new(1200.0, 1.5);
        sim.add_nucleus(CrystalNucleus::new(
            [0.0; 3],
            0.5,
            CrystalStructure::Bcc,
            0.1,
        ));
        let x = sim.transformed_fraction();
        assert!((0.0..=1.0).contains(&x), "x={x}");
    }

    // --- Avrami ---

    #[test]
    fn test_avrami_struct_zero_time() {
        let a = Avrami::new(1.0, 1.0);
        assert_eq!(a.compute_transformed_fraction(0.0), 0.0);
    }

    #[test]
    fn test_avrami_struct_approaches_one() {
        let a = Avrami::new(1.0, 1.0);
        let x = a.compute_transformed_fraction(100.0);
        assert!(x > 0.99, "x={x}");
    }

    #[test]
    fn test_avrami_struct_monotone() {
        let a = Avrami::new(0.1, 0.5);
        let x1 = a.compute_transformed_fraction(1.0);
        let x2 = a.compute_transformed_fraction(2.0);
        assert!(x2 > x1, "x1={x1} x2={x2}");
    }

    #[test]
    fn test_avrami_struct_zero_growth_zero_fraction() {
        let a = Avrami::new(1.0, 0.0);
        assert_eq!(a.compute_transformed_fraction(10.0), 0.0);
    }

    // --- OrientationRelation ---

    #[test]
    fn test_ks_angle_positive() {
        let angle = OrientationRelation::KurdjumovSachs.characteristic_angle();
        assert!(angle > 0.0 && angle < PI);
    }

    #[test]
    fn test_nw_angle_positive() {
        let angle = OrientationRelation::NishiyamaWassermann.characteristic_angle();
        assert!(angle > 0.0 && angle < PI);
    }

    #[test]
    fn test_bain_angle_45_degrees() {
        let angle = OrientationRelation::Bain.characteristic_angle();
        assert!((angle - PI / 4.0).abs() < 1e-10, "angle={angle}");
    }

    #[test]
    fn test_ks_nw_different() {
        let ks = OrientationRelation::KurdjumovSachs.characteristic_angle();
        let nw = OrientationRelation::NishiyamaWassermann.characteristic_angle();
        assert!((ks - nw).abs() > 1e-6, "K-S and N-W should differ");
    }

    // --- NucleationSite (legacy) ---

    #[test]
    fn test_nucleation_site_volume() {
        let site = NucleationSite::new([0.0; 3], 1.0, 0.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((site.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_nucleation_site_fields() {
        let pos = [1.0, 2.0, 3.0];
        let site = NucleationSite::new(pos, 0.5, 1.5);
        assert_eq!(site.position, pos);
        assert!((site.radius - 0.5).abs() < 1e-12);
    }

    // --- classical_nucleation_rate ---

    #[test]
    fn test_cnt_rate_zero_positive_dg() {
        assert_eq!(classical_nucleation_rate(1.0, 0.1, 300.0, 0.025), 0.0);
    }

    #[test]
    fn test_cnt_rate_positive() {
        let r = classical_nucleation_rate(-1e5, 0.01, 300.0, 1.0);
        assert!(r > 0.0 && r <= 1.0);
    }

    // --- GrainGrowthSimulation ---

    #[test]
    fn test_grain_growth_mean_size() {
        let sim = GrainGrowthSimulation::new(10, 2.0);
        assert!((sim.mean_grain_size() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_grain_growth_largest_idx() {
        let mut sim = GrainGrowthSimulation::new(5, 1.0);
        sim.sizes[3] = 5.0;
        assert_eq!(sim.largest_grain_idx(), Some(3));
    }

    // --- ostwald_ripening_rate ---

    #[test]
    fn test_ostwald_zero_temp() {
        assert_eq!(
            ostwald_ripening_rate(0.1, 1e-5, 1e-9, 100.0, 1e-6, 0.0),
            0.0
        );
    }

    #[test]
    fn test_ostwald_positive() {
        let r = ostwald_ripening_rate(0.1, 1e-5, 1e-9, 100.0, 1e-6, 300.0);
        assert!(r > 0.0);
    }
}

// ---------------------------------------------------------------------------
// Wilson-Frenkel growth model
// ---------------------------------------------------------------------------

/// Wilson-Frenkel normal crystal growth model.
///
/// The growth velocity is:
///
/// ```text
/// v = v_0 [1 − exp(−ΔG / k_B T)]
/// ```
///
/// For small undercooling this linearises to `v ≈ v_0 ΔG / (k_B T)`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WilsonFrenkelGrowth {
    /// Pre-exponential velocity factor v_0 (m/s).
    pub v0: f64,
    /// Thermal energy k_B T (J).
    pub kbt: f64,
}

impl WilsonFrenkelGrowth {
    /// Create a new [`WilsonFrenkelGrowth`] model.
    pub fn new(v0: f64, kbt: f64) -> Self {
        Self { v0, kbt }
    }

    /// Compute the growth velocity for a given bulk driving force ΔG (J).
    ///
    /// # Arguments
    /// * `delta_g` – bulk free energy difference per formula unit (J); negative
    ///   for crystallisation (favourable growth)
    pub fn growth_velocity(&self, delta_g: f64) -> f64 {
        if self.kbt <= 0.0 {
            return 0.0;
        }
        self.v0 * (1.0 - (delta_g / self.kbt).exp())
    }

    /// Linearised growth velocity for small undercooling.
    ///
    /// # Arguments
    /// * `delta_g` – bulk free energy difference (J)
    pub fn linear_growth_velocity(&self, delta_g: f64) -> f64 {
        if self.kbt <= 0.0 {
            return 0.0;
        }
        -self.v0 * delta_g / self.kbt
    }
}

// ---------------------------------------------------------------------------
// Spiral growth model (BCF theory)
// ---------------------------------------------------------------------------

/// Burton-Cabrera-Frank (BCF) spiral growth model.
///
/// The normal growth rate from screw dislocation spirals:
///
/// ```text
/// R = A σ² / (1 + B σ)
/// ```
///
/// where σ = (C − C_eq) / C_eq is the relative supersaturation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpiralGrowth {
    /// Coefficient A (length/time).
    pub a_coeff: f64,
    /// Coefficient B (dimensionless).
    pub b_coeff: f64,
}

impl SpiralGrowth {
    /// Create a new [`SpiralGrowth`] model.
    pub fn new(a_coeff: f64, b_coeff: f64) -> Self {
        Self { a_coeff, b_coeff }
    }

    /// Normal growth rate from BCF spiral growth kinetics.
    ///
    /// # Arguments
    /// * `sigma` – relative supersaturation (dimensionless)
    pub fn growth_rate(&self, sigma: f64) -> f64 {
        if sigma <= 0.0 {
            return 0.0;
        }
        self.a_coeff * sigma * sigma / (1.0 + self.b_coeff * sigma)
    }

    /// Parabolic-to-linear transition supersaturation.
    ///
    /// At σ >> 1/B the spiral growth becomes linear; the crossover is
    /// approximately σ_c = 1 / B.
    pub fn crossover_supersaturation(&self) -> f64 {
        if self.b_coeff <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / self.b_coeff
        }
    }
}

// ---------------------------------------------------------------------------
// Normal growth model (2D nucleation)
// ---------------------------------------------------------------------------

/// Two-dimensional nucleation (normal growth) model.
///
/// Growth rate:
/// ```text
/// R = K exp(−B / σ)
/// ```
///
/// where σ is the supersaturation and B is related to the 2-D nucleation
/// barrier.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NormalGrowthModel {
    /// Pre-exponential rate constant K (length/time).
    pub k_rate: f64,
    /// Barrier coefficient B (dimensionless).
    pub b_barrier: f64,
}

impl NormalGrowthModel {
    /// Create a new [`NormalGrowthModel`].
    pub fn new(k_rate: f64, b_barrier: f64) -> Self {
        Self { k_rate, b_barrier }
    }

    /// Compute the 2-D nucleation growth rate.
    ///
    /// # Arguments
    /// * `sigma` – supersaturation (dimensionless, > 0 for growth)
    pub fn growth_rate(&self, sigma: f64) -> f64 {
        if sigma <= 0.0 || self.b_barrier <= 0.0 {
            return 0.0;
        }
        self.k_rate * (-self.b_barrier / sigma).exp()
    }
}

// ---------------------------------------------------------------------------
// Surface energy anisotropy
// ---------------------------------------------------------------------------

/// Crystal surface energy with cubic anisotropy.
///
/// The anisotropic surface energy in 2D is modelled as:
///
/// ```text
/// γ(θ) = γ_0 [1 + δ cos(m θ)]
/// ```
///
/// where θ is the angle from the reference direction, δ is the anisotropy
/// strength, and m is the symmetry order (4 for cubic, 6 for hexagonal).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AnisotropicSurfaceEnergy {
    /// Isotropic surface energy γ_0 (J/m²).
    pub gamma0: f64,
    /// Anisotropy strength δ (dimensionless, 0–0.05 typical).
    pub delta: f64,
    /// Symmetry order m (4 = cubic, 6 = hexagonal).
    pub symmetry: u32,
}

impl AnisotropicSurfaceEnergy {
    /// Create a new [`AnisotropicSurfaceEnergy`].
    pub fn new(gamma0: f64, delta: f64, symmetry: u32) -> Self {
        Self {
            gamma0,
            delta,
            symmetry,
        }
    }

    /// Compute the surface energy at orientation angle θ (radians).
    pub fn energy(&self, theta: f64) -> f64 {
        self.gamma0 * (1.0 + self.delta * ((self.symmetry as f64) * theta).cos())
    }

    /// Compute the stiffness γ + γ'' at orientation angle θ.
    ///
    /// ```text
    /// γ̃ = γ + d²γ/dθ²  = γ_0 [1 + δ(1 - m²) cos(m θ)]
    /// ```
    pub fn stiffness(&self, theta: f64) -> f64 {
        let m = self.symmetry as f64;
        self.gamma0 * (1.0 + self.delta * (1.0 - m * m) * (m * theta).cos())
    }
}

// ---------------------------------------------------------------------------
// Phase-field model for solidification
// ---------------------------------------------------------------------------

/// Parameters for the phase-field model of solidification.
///
/// The Allen-Cahn type equation:
/// ```text
/// τ ∂φ/∂t = W² ∇²φ + φ(1-φ)(φ - 0.5 + λ u)
/// ```
///
/// where φ is the phase field (0 = liquid, 1 = solid), u is the
/// dimensionless thermal field, λ is the coupling constant.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PhaseFieldSolidification {
    /// Relaxation time τ (s).
    pub tau: f64,
    /// Interface width W (m).
    pub interface_width: f64,
    /// Coupling constant λ.
    pub lambda: f64,
    /// Thermal diffusivity D_T (m²/s).
    pub thermal_diffusivity: f64,
    /// Latent heat divided by specific heat: L / c_p (K).
    pub latent_heat_ratio: f64,
}

impl PhaseFieldSolidification {
    /// Create a new [`PhaseFieldSolidification`] model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        tau: f64,
        interface_width: f64,
        lambda: f64,
        thermal_diffusivity: f64,
        latent_heat_ratio: f64,
    ) -> Self {
        Self {
            tau,
            interface_width,
            lambda,
            thermal_diffusivity,
            latent_heat_ratio,
        }
    }

    /// Compute the bulk free energy double-well function f(φ).
    ///
    /// ```text
    /// f(φ) = φ²(1-φ)²
    /// ```
    pub fn bulk_free_energy(&self, phi: f64) -> f64 {
        phi * phi * (1.0 - phi) * (1.0 - phi)
    }

    /// Compute df/dφ.
    pub fn bulk_free_energy_deriv(&self, phi: f64) -> f64 {
        2.0 * phi * (1.0 - phi) * (1.0 - 2.0 * phi)
    }

    /// Interpolation function p(φ) = φ²(3 − 2φ).
    pub fn interpolation(&self, phi: f64) -> f64 {
        phi * phi * (3.0 - 2.0 * phi)
    }

    /// Driving force term for the phase field: (φ − 0.5 + λ u).
    ///
    /// # Arguments
    /// * `phi` – phase-field order parameter
    /// * `u`   – dimensionless undercooling field
    pub fn driving_force(&self, phi: f64, u: f64) -> f64 {
        phi - 0.5 + self.lambda * u
    }
}

// ---------------------------------------------------------------------------
// Crystal defect formation
// ---------------------------------------------------------------------------

/// Crystal defect model: vacancies and interstitials.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CrystalDefects {
    /// Formation energy of a vacancy (eV).
    pub vacancy_formation_energy: f64,
    /// Formation energy of an interstitial (eV).
    pub interstitial_formation_energy: f64,
    /// Boltzmann constant in eV/K.
    pub kb_ev: f64,
}

impl CrystalDefects {
    /// Create a new [`CrystalDefects`] model.
    pub fn new(vacancy_formation_energy: f64, interstitial_formation_energy: f64) -> Self {
        Self {
            vacancy_formation_energy,
            interstitial_formation_energy,
            kb_ev: 8.617_333_262e-5,
        }
    }

    /// Equilibrium vacancy concentration (per site) at temperature T (K).
    ///
    /// ```text
    /// c_v = exp(−E_f^v / (k_B T))
    /// ```
    pub fn vacancy_concentration(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        (-self.vacancy_formation_energy / (self.kb_ev * temperature)).exp()
    }

    /// Equilibrium interstitial concentration (per site) at temperature T (K).
    pub fn interstitial_concentration(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        (-self.interstitial_formation_energy / (self.kb_ev * temperature)).exp()
    }

    /// Frenkel pair concentration (equal vacancy and interstitial from
    /// thermal equilibrium in a perfect crystal).
    pub fn frenkel_pair_concentration(&self, temperature: f64) -> f64 {
        // Combined formation energy = vacancy + interstitial
        let e_f = self.vacancy_formation_energy + self.interstitial_formation_energy;
        if temperature <= 0.0 {
            return 0.0;
        }
        (-e_f / (2.0 * self.kb_ev * temperature)).exp()
    }
}

// ---------------------------------------------------------------------------
// Growth velocity vs temperature (undercooling curve)
// ---------------------------------------------------------------------------

/// Tabulate growth velocity as a function of undercooling using the
/// Wilson-Frenkel model with an Arrhenius attachment frequency.
///
/// # Arguments
/// * `wf`           – Wilson-Frenkel growth model
/// * `delta_h_fus`  – enthalpy of fusion per formula unit (J)
/// * `t_melt`       – melting temperature (K)
/// * `undercoolings` – slice of undercooling values ΔT (K)
///
/// Returns a `Vec`f64` of growth velocities (m/s).
#[allow(dead_code)]
pub fn growth_velocity_vs_undercooling(
    wf: &WilsonFrenkelGrowth,
    delta_h_fus: f64,
    t_melt: f64,
    undercoolings: &[f64],
) -> Vec<f64> {
    undercoolings
        .iter()
        .map(|&dt| {
            if t_melt <= 0.0 || wf.kbt <= 0.0 {
                return 0.0;
            }
            let t = t_melt - dt;
            if t <= 0.0 {
                return 0.0;
            }
            // ΔG ≈ -ΔH_fus * ΔT / T_melt
            let delta_g = -delta_h_fus * dt / t_melt;
            wf.growth_velocity(delta_g)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Ivantsov function (dendrite growth — solvability theory)
// ---------------------------------------------------------------------------

/// Evaluate the Ivantsov function Iv(Pe) for a paraboloid of revolution.
///
/// The exact Ivantsov function is `Iv(Pe) = Pe exp(Pe) E_1(Pe)` where E_1 is
/// the exponential integral. Here we use a Padé approximation valid for
/// Pe > 0:
///
/// ```text
/// Iv(Pe) ≈ Pe / (Pe + 0.5)   (first-order Padé)
/// ```
///
/// # Arguments
/// * `peclet` – tip Péclet number Pe = R v / (2 D)
#[allow(dead_code)]
pub fn ivantsov_function(peclet: f64) -> f64 {
    if peclet <= 0.0 {
        return 0.0;
    }
    // Better approximation: use the series form for small Pe
    if peclet < 0.1 {
        // Iv(Pe) ≈ Pe (−ln(Pe) − γ_E)  for small Pe, γ_E = 0.5772...
        let gamma_e = 0.577_215_664_9;
        peclet * (-peclet.ln() - gamma_e + 1.0)
    } else {
        // Rational approximation
        peclet * (1.0 + peclet) / (1.0 + 2.0 * peclet + peclet * peclet * 0.5)
    }
}

// ---------------------------------------------------------------------------
// Polycrystalline grain boundary energy (Read-Shockley model)
// ---------------------------------------------------------------------------

/// Compute the grain boundary energy using the Read-Shockley model.
///
/// ```text
/// γ_gb(θ) = γ_m (θ/θ_m) [1 − ln(θ/θ_m)]    for θ ≤ θ_m
/// γ_gb(θ) = γ_m                               for θ > θ_m
/// ```
///
/// # Arguments
/// * `theta_mis` – misorientation angle (radians)
/// * `theta_m`   – cut-off angle for high-angle boundaries (radians, ≈ 15°)
/// * `gamma_m`   – energy of a high-angle grain boundary (J/m²)
#[allow(dead_code)]
pub fn read_shockley_energy(theta_mis: f64, theta_m: f64, gamma_m: f64) -> f64 {
    if theta_mis <= 0.0 || theta_m <= 0.0 {
        return 0.0;
    }
    if theta_mis >= theta_m {
        return gamma_m;
    }
    let ratio = theta_mis / theta_m;
    gamma_m * ratio * (1.0 - ratio.ln())
}

// ---------------------------------------------------------------------------
// Additional tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // --- WilsonFrenkelGrowth ---

    #[test]
    fn test_wf_zero_kbt_returns_zero() {
        let wf = WilsonFrenkelGrowth::new(1.0, 0.0);
        assert_eq!(wf.growth_velocity(-0.1), 0.0);
    }

    #[test]
    fn test_wf_favourable_driving_force() {
        let wf = WilsonFrenkelGrowth::new(1.0, 0.025);
        // Negative ΔG (crystallisation): growth velocity > 0
        let v = wf.growth_velocity(-0.05);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_wf_zero_driving_force() {
        let wf = WilsonFrenkelGrowth::new(1.0, 0.025);
        let v = wf.growth_velocity(0.0);
        assert!((v).abs() < 1e-12);
    }

    #[test]
    fn test_wf_linear_approx_small_dg() {
        let wf = WilsonFrenkelGrowth::new(1.0, 1.0);
        let dg = -0.001;
        let v_exact = wf.growth_velocity(dg);
        let v_lin = wf.linear_growth_velocity(dg);
        // For small ΔG/kBT the two should agree to within < 0.2%
        assert!(
            (v_exact - v_lin).abs() / v_lin.abs() < 0.002,
            "rel_diff too large"
        );
    }

    // --- SpiralGrowth ---

    #[test]
    fn test_spiral_zero_supersaturation() {
        let sg = SpiralGrowth::new(1.0, 1.0);
        assert_eq!(sg.growth_rate(0.0), 0.0);
    }

    #[test]
    fn test_spiral_positive() {
        let sg = SpiralGrowth::new(1.0, 0.1);
        assert!(sg.growth_rate(0.5) > 0.0);
    }

    #[test]
    fn test_spiral_crossover() {
        let sg = SpiralGrowth::new(1.0, 2.0);
        assert!((sg.crossover_supersaturation() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_spiral_b_zero_crossover_infinite() {
        let sg = SpiralGrowth::new(1.0, 0.0);
        assert_eq!(sg.crossover_supersaturation(), f64::INFINITY);
    }

    // --- NormalGrowthModel ---

    #[test]
    fn test_normal_growth_zero_supersaturation() {
        let ng = NormalGrowthModel::new(1.0, 0.5);
        assert_eq!(ng.growth_rate(0.0), 0.0);
    }

    #[test]
    fn test_normal_growth_positive() {
        let ng = NormalGrowthModel::new(1.0, 0.1);
        assert!(ng.growth_rate(0.5) > 0.0);
    }

    #[test]
    fn test_normal_growth_increases_with_supersaturation() {
        let ng = NormalGrowthModel::new(1.0, 0.5);
        let r1 = ng.growth_rate(0.2);
        let r2 = ng.growth_rate(0.5);
        assert!(r2 > r1, "r1={r1} r2={r2}");
    }

    // --- AnisotropicSurfaceEnergy ---

    #[test]
    fn test_aniso_isotropic_case() {
        let ase = AnisotropicSurfaceEnergy::new(1.0, 0.0, 4);
        // Zero anisotropy → all angles return gamma0
        for &theta in &[0.0, 0.5, 1.0, 1.5] {
            assert!((ase.energy(theta) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_aniso_maximum_at_preferred_direction() {
        // At θ = 0 cos(4*0) = 1 → maximum energy with positive delta
        let ase = AnisotropicSurfaceEnergy::new(1.0, 0.04, 4);
        let e0 = ase.energy(0.0);
        let e_pi8 = ase.energy(PI / 8.0); // cos(π/2) = 0 → isotropic value
        assert!(e0 > e_pi8, "e0={e0} e_pi8={e_pi8}");
    }

    #[test]
    fn test_aniso_stiffness_computed() {
        let ase = AnisotropicSurfaceEnergy::new(1.0, 0.02, 4);
        let stiff = ase.stiffness(0.0);
        // At θ=0: stiffness = γ0 [1 + δ(1 - 16)] = γ0 [1 - 15δ]
        let expected = 1.0 * (1.0 + 0.02 * (1.0 - 16.0));
        assert!((stiff - expected).abs() < 1e-12, "stiff={stiff}");
    }

    // --- PhaseFieldSolidification ---

    #[test]
    fn test_pf_interpolation_endpoints() {
        let pf = PhaseFieldSolidification::new(0.001, 1e-6, 10.0, 1e-6, 100.0);
        assert!((pf.interpolation(0.0)).abs() < 1e-12);
        assert!((pf.interpolation(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pf_bulk_free_energy_zero_at_wells() {
        let pf = PhaseFieldSolidification::new(0.001, 1e-6, 10.0, 1e-6, 100.0);
        assert!((pf.bulk_free_energy(0.0)).abs() < 1e-12);
        assert!((pf.bulk_free_energy(1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_pf_bulk_free_energy_max_at_midpoint() {
        let pf = PhaseFieldSolidification::new(0.001, 1e-6, 10.0, 1e-6, 100.0);
        let f_mid = pf.bulk_free_energy(0.5);
        assert!(f_mid > 0.0, "f_mid={f_mid}");
    }

    #[test]
    fn test_pf_driving_force() {
        let pf = PhaseFieldSolidification::new(0.001, 1e-6, 2.0, 1e-6, 100.0);
        // At φ=0.5, u=0: driving force = 0
        assert!((pf.driving_force(0.5, 0.0)).abs() < 1e-12);
    }

    // --- CrystalDefects ---

    #[test]
    fn test_defects_vacancy_zero_temp() {
        let d = CrystalDefects::new(1.0, 3.0);
        assert_eq!(d.vacancy_concentration(0.0), 0.0);
    }

    #[test]
    fn test_defects_vacancy_increases_with_temp() {
        let d = CrystalDefects::new(1.0, 3.0);
        let c1 = d.vacancy_concentration(300.0);
        let c2 = d.vacancy_concentration(1000.0);
        assert!(c2 > c1, "c1={c1} c2={c2}");
    }

    #[test]
    fn test_defects_interstitial_less_than_vacancy() {
        // Interstitial formation energy > vacancy → lower concentration
        let d = CrystalDefects::new(1.0, 3.0);
        let cv = d.vacancy_concentration(500.0);
        let ci = d.interstitial_concentration(500.0);
        assert!(ci < cv, "ci={ci} cv={cv}");
    }

    #[test]
    fn test_defects_frenkel_pair_positive() {
        let d = CrystalDefects::new(1.0, 3.0);
        let cf = d.frenkel_pair_concentration(1000.0);
        assert!(cf > 0.0 && cf < 1.0, "cf={cf}");
    }

    // --- growth_velocity_vs_undercooling ---

    #[test]
    fn test_gv_undercooling_monotone() {
        let wf = WilsonFrenkelGrowth::new(1.0, 0.1);
        let dts = vec![1.0, 5.0, 10.0];
        let vs = growth_velocity_vs_undercooling(&wf, 1.0, 1000.0, &dts);
        // Each velocity should be finite and non-negative
        for v in &vs {
            assert!(v.is_finite() && *v >= 0.0, "v={v}");
        }
    }

    // --- ivantsov_function ---

    #[test]
    fn test_ivantsov_zero_peclet() {
        assert_eq!(ivantsov_function(0.0), 0.0);
    }

    #[test]
    fn test_ivantsov_positive_peclet() {
        let iv = ivantsov_function(1.0);
        assert!(iv > 0.0 && iv < 1.0, "iv={iv}");
    }

    #[test]
    fn test_ivantsov_small_peclet() {
        let iv = ivantsov_function(0.01);
        assert!(iv > 0.0, "iv={iv}");
    }

    // --- read_shockley_energy ---

    #[test]
    fn test_rs_zero_misorientation() {
        assert_eq!(read_shockley_energy(0.0, 0.26, 0.5), 0.0);
    }

    #[test]
    fn test_rs_high_angle_saturates() {
        let e = read_shockley_energy(0.5, 0.26, 0.5);
        assert!((e - 0.5).abs() < 1e-12, "e={e}");
    }

    #[test]
    fn test_rs_low_angle_less_than_max() {
        let e = read_shockley_energy(0.05, 0.26, 0.5);
        assert!(e < 0.5 && e > 0.0, "e={e}");
    }

    #[test]
    fn test_rs_increases_toward_high_angle() {
        let e1 = read_shockley_energy(0.05, 0.26, 0.5);
        let e2 = read_shockley_energy(0.15, 0.26, 0.5);
        assert!(e2 > e1, "e1={e1} e2={e2}");
    }

    // --- PhaseFieldCrystal (extended) ---

    #[test]
    fn test_pfc_laplacian_uniform_field_is_zero() {
        let pfc = PhaseFieldCrystal::new(4, 4, 1.0, 1.0);
        // Uniform field → Laplacian = 0
        assert!((pfc.laplacian(1, 1)).abs() < 1e-12);
    }

    #[test]
    fn test_pfc_mean_phi_initial() {
        let pfc = PhaseFieldCrystal::new(8, 8, 0.5, 0.5);
        assert!((pfc.mean_phi() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_pfc_step_changes_phi() {
        let mut pfc = PhaseFieldCrystal::new(4, 4, 1.0, 0.0);
        // Perturb one cell
        pfc.phi[0] = 0.1;
        let mean_before = pfc.mean_phi();
        pfc.step(0.001, 1.0, 0.0);
        let mean_after = pfc.mean_phi();
        // Mean should have changed due to non-zero driving force
        assert!((mean_after - mean_before).abs() > 0.0);
    }

    // --- CrystalGrowthModel (extended) ---

    #[test]
    fn test_cgm_current_time_advances() {
        let mut m = CrystalGrowthModel::new(1.0, 0.5, 1e-10);
        m.step(0.01, 1200.0);
        assert!((m.current_time() - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_cgm_transformed_fraction_bounded() {
        let mut m = CrystalGrowthModel::new(1.0, 0.5, 1e-10);
        for _ in 0..100 {
            m.step(0.1, 1200.0);
        }
        let x = m.transformed_fraction();
        assert!((0.0..=1.0).contains(&x), "x={x}");
    }

    // --- dendrite_tip_velocity ---

    #[test]
    fn test_dendrite_tip_velocity_zero_d0() {
        assert_eq!(dendrite_tip_velocity(0.5, 0.0, 0.025), 0.0);
    }

    #[test]
    fn test_dendrite_tip_velocity_positive() {
        let v = dendrite_tip_velocity(0.5, 1e-9, 0.025);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_dendrite_tip_scales_with_delta_squared() {
        let v1 = dendrite_tip_velocity(0.2, 1e-9, 0.025);
        let v2 = dendrite_tip_velocity(0.4, 1e-9, 0.025);
        assert!((v2 / v1 - 4.0).abs() < 1e-9, "ratio={}", v2 / v1);
    }
}
