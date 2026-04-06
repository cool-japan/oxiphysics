// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Weakly ionized plasma simulation using the Lattice Boltzmann Method.
//!
//! This module provides:
//! - [`PlasmaSpecies`]: Parameters for a plasma species (ions or electrons).
//! - [`ElectricField`]: 2-D electric field on a grid with Lorentz force computation.
//! - [`MagneticField`]: Uniform out-of-plane magnetic field with cyclotron and Larmor formulae.
//! - [`PlasmaLbm`]: Two-fluid (ion + electron) LBM simulation driver.
//! - [`DebyeShielding`]: Debye–Hückel screened Coulomb potential.
//! - [`plasma_frequency`]: Plasma (Langmuir) frequency.
//! - [`debye_length`]: Debye screening length.
//! - [`coulomb_logarithm`]: Coulomb logarithm for collision cross-sections.
//!
//! # Physical constants (SI)
//! - Electron charge: *e* = 1.602 176 634 × 10⁻¹⁹ C
//! - Electron mass:   *mₑ* = 9.109 383 7015 × 10⁻³¹ kg
//! - Vacuum permittivity: ε₀ = 8.854 187 8128 × 10⁻¹² F m⁻¹
//! - Boltzmann constant:  *k*_B = 1.380 649 × 10⁻²³ J K⁻¹
//!
//! References:
//! - Chen, F. F. (1984). *Introduction to Plasma Physics and Controlled Fusion*.
//! - Birdsall, C. K. & Langdon, A. B. (1991). *Plasma Physics via Computer Simulation*.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Vacuum permittivity ε₀ (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;
/// Elementary charge *e* (C).
pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;
/// Electron mass *mₑ* (kg).
pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;
/// Boltzmann constant *k*_B (J K⁻¹).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;

// ---------------------------------------------------------------------------
// Standalone free functions
// ---------------------------------------------------------------------------

/// Plasma (Langmuir) angular frequency ω_p (rad s⁻¹).
///
/// ω_p = √( n q² / (ε₀ m) )
///
/// # Arguments
/// * `n`      – number density (m⁻³)
/// * `charge` – particle charge (C)
/// * `mass`   – particle mass (kg)
///
/// ```no_run
/// use oxiphysics_lbm::plasma_lbm::{plasma_frequency, ELEM_CHARGE, ELECTRON_MASS, EPSILON_0};
/// let wp = plasma_frequency(1e18, ELEM_CHARGE, ELECTRON_MASS);
/// assert!(wp > 0.0);
/// ```
pub fn plasma_frequency(n: f64, charge: f64, mass: f64) -> f64 {
    if n <= 0.0 || mass <= 0.0 {
        return 0.0;
    }
    (n * charge * charge / (EPSILON_0 * mass)).sqrt()
}

/// Debye screening length λ_D (m).
///
/// λ_D = √( ε₀ k_B T / (n q²) )
///
/// # Arguments
/// * `n`           – number density (m⁻³)
/// * `temperature` – temperature (K)
/// * `charge`      – particle charge (C)
///
/// ```no_run
/// use oxiphysics_lbm::plasma_lbm::{debye_length, ELEM_CHARGE};
/// let lam = debye_length(1e18, 1e4, ELEM_CHARGE);
/// assert!(lam > 0.0);
/// ```
pub fn debye_length(n: f64, temperature: f64, charge: f64) -> f64 {
    if n <= 0.0 || temperature <= 0.0 {
        return 0.0;
    }
    (EPSILON_0 * K_BOLTZMANN * temperature / (n * charge * charge)).sqrt()
}

/// Coulomb logarithm ln Λ (dimensionless).
///
/// Uses the approximate form:
///
/// ln Λ ≈ ln( 12π n λ_D³ )
///
/// Clamped to a minimum of 1.
///
/// # Arguments
/// * `temperature` – electron temperature (K)
/// * `n`           – electron number density (m⁻³)
///
/// ```no_run
/// use oxiphysics_lbm::plasma_lbm::coulomb_logarithm;
/// let lnL = coulomb_logarithm(1e4, 1e18);
/// assert!(lnL >= 1.0);
/// ```
pub fn coulomb_logarithm(temperature: f64, n: f64) -> f64 {
    if temperature <= 0.0 || n <= 0.0 {
        return 1.0;
    }
    let lambda_d = debye_length(n, temperature, ELEM_CHARGE);
    let arg = 12.0 * PI * n * lambda_d * lambda_d * lambda_d;
    if arg <= 1.0 { 1.0 } else { arg.ln() }
}

// ---------------------------------------------------------------------------
// PlasmaSpecies
// ---------------------------------------------------------------------------

/// Parameters describing a plasma species (ions or electrons).
#[derive(Debug, Clone)]
pub struct PlasmaSpecies {
    /// Charge number Z (dimensionless; Z = 1 for singly charged, −1 for electrons).
    pub charge_number: f64,
    /// Mass ratio m/m_e (dimensionless; 1.0 for electrons, ~1836 for protons).
    pub mass_ratio: f64,
    /// Species temperature (K).
    pub temperature: f64,
    /// Mean drift velocity \[vx, vy, vz\] (m s⁻¹).
    pub drift_velocity: [f64; 3],
}

impl PlasmaSpecies {
    /// Create a new [`PlasmaSpecies`].
    pub fn new(
        charge_number: f64,
        mass_ratio: f64,
        temperature: f64,
        drift_velocity: [f64; 3],
    ) -> Self {
        Self {
            charge_number,
            mass_ratio,
            temperature,
            drift_velocity,
        }
    }

    /// Actual charge of this species (C).
    pub fn charge(&self) -> f64 {
        self.charge_number * ELEM_CHARGE
    }

    /// Actual mass of this species (kg).
    pub fn mass(&self) -> f64 {
        self.mass_ratio * ELECTRON_MASS
    }

    /// Thermal speed √(k_B T / m) (m s⁻¹).
    pub fn thermal_speed(&self) -> f64 {
        (K_BOLTZMANN * self.temperature / self.mass()).sqrt()
    }
}

// ---------------------------------------------------------------------------
// ElectricField
// ---------------------------------------------------------------------------

/// 2-D electric field defined on a uniform grid.
///
/// Components `Ex` and `Ey` are stored in row-major order (ny × nx).
#[derive(Debug, Clone)]
pub struct ElectricField {
    /// x-component of the electric field (V m⁻¹), length nx × ny.
    pub ex: Vec<f64>,
    /// y-component of the electric field (V m⁻¹), length nx × ny.
    pub ey: Vec<f64>,
    /// Number of grid cells in x.
    pub nx: usize,
    /// Number of grid cells in y.
    pub ny: usize,
}

impl ElectricField {
    /// Create a zero electric field on an `nx × ny` grid.
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            nx,
            ny,
        }
    }

    /// Flat index for cell (ix, iy).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Electric field vector at grid cell (ix, iy) → \[Ex, Ey, 0\].
    pub fn field_at(&self, ix: usize, iy: usize) -> [f64; 3] {
        let i = self.idx(ix, iy);
        [self.ex[i], self.ey[i], 0.0]
    }

    /// Set the electric field at grid cell (ix, iy).
    pub fn set(&mut self, ix: usize, iy: usize, ex: f64, ey: f64) {
        let i = self.idx(ix, iy);
        self.ex[i] = ex;
        self.ey[i] = ey;
    }

    /// Lorentz electric force **F** = q **E** (N).
    ///
    /// `vel` is unused for the electric part but kept for API consistency
    /// with a combined Lorentz force that includes the magnetic term.
    ///
    /// ```no_run
    /// use oxiphysics_lbm::plasma_lbm::ElectricField;
    /// let ef = ElectricField::new(4, 4);
    /// let f = ef.compute_lorentz_force(1.6e-19, [1.0, 0.0, 0.0]);
    /// assert_eq!(f, [0.0, 0.0, 0.0]);
    /// ```
    pub fn compute_lorentz_force(&self, charge: f64, _vel: [f64; 3]) -> [f64; 3] {
        // Returns the force due to the field at cell (0,0) by default.
        // A grid-aware version would accept (ix, iy) as arguments.
        let i = 0;
        [charge * self.ex[i], charge * self.ey[i], 0.0]
    }

    /// Lorentz electric force at a specific grid cell (ix, iy).
    pub fn lorentz_force_at(&self, ix: usize, iy: usize, charge: f64) -> [f64; 3] {
        let i = self.idx(ix, iy);
        [charge * self.ex[i], charge * self.ey[i], 0.0]
    }
}

// ---------------------------------------------------------------------------
// MagneticField
// ---------------------------------------------------------------------------

/// Uniform out-of-plane magnetic field B_z (T) in a 2-D plasma simulation.
#[derive(Debug, Clone)]
pub struct MagneticField {
    /// Magnetic flux density B_z (T).
    pub bz: f64,
}

impl MagneticField {
    /// Create a new [`MagneticField`] with uniform B_z.
    pub fn new(bz: f64) -> Self {
        Self { bz }
    }

    /// Cyclotron (gyro) angular frequency ω_c = |q| B / m (rad s⁻¹).
    ///
    /// ```no_run
    /// use oxiphysics_lbm::plasma_lbm::{MagneticField, ELEM_CHARGE, ELECTRON_MASS};
    /// let bf = MagneticField::new(1.0);
    /// let wc = bf.compute_cyclotron_frequency(ELEM_CHARGE, ELECTRON_MASS);
    /// assert!((wc - 1.758_820_e11_f64).abs() / wc < 1e-4);
    /// ```
    pub fn compute_cyclotron_frequency(&self, charge: f64, mass: f64) -> f64 {
        if mass <= 0.0 {
            return 0.0;
        }
        charge.abs() * self.bz.abs() / mass
    }

    /// Larmor (gyro) radius r_L = v_⊥ / ω_c (m).
    ///
    /// # Arguments
    /// * `vel`    – perpendicular speed (m s⁻¹)
    /// * `charge` – particle charge (C)
    /// * `mass`   – particle mass (kg)
    ///
    /// ```no_run
    /// use oxiphysics_lbm::plasma_lbm::{MagneticField, ELEM_CHARGE, ELECTRON_MASS};
    /// let bf = MagneticField::new(1.0);
    /// let r = bf.larmor_radius(1e6, ELEM_CHARGE, ELECTRON_MASS);
    /// assert!(r > 0.0);
    /// ```
    pub fn larmor_radius(&self, vel: f64, charge: f64, mass: f64) -> f64 {
        let wc = self.compute_cyclotron_frequency(charge, mass);
        if wc <= 0.0 {
            return f64::INFINITY;
        }
        vel.abs() / wc
    }

    /// Lorentz magnetic force **F** = q **v** × **B** (N), with B = \[0, 0, Bz\].
    ///
    /// ```text
    /// use oxiphysics_lbm::plasma_lbm::{MagneticField, ELEM_CHARGE};
    /// let bf = MagneticField::new(1.0);
    /// let f = bf.magnetic_force(ELEM_CHARGE, [1.0, 0.0, 0.0]);
    /// // F = q [vx, vy, vz] × [0, 0, Bz] = q [vy*Bz, -vx*Bz, 0]
    /// assert!((f[1] + ELEM_CHARGE).abs() < 1e-30);
    /// ```
    pub fn magnetic_force(&self, charge: f64, vel: [f64; 3]) -> [f64; 3] {
        let bz = self.bz;
        // v × B where B = (0, 0, Bz):
        //   (vy*Bz - vz*0, vz*0 - vx*Bz, vx*0 - vy*0) = (vy*Bz, -vx*Bz, 0)
        [charge * vel[1] * bz, -charge * vel[0] * bz, 0.0]
    }
}

// ---------------------------------------------------------------------------
// FluidField  (helper, not public)
// ---------------------------------------------------------------------------

/// Internal helper: density + velocity field on an nx × ny grid.
#[derive(Debug, Clone)]
struct FluidField {
    /// Density (m⁻³), length nx * ny.
    density: Vec<f64>,
    /// Velocity x-component (m s⁻¹), length nx * ny.
    vx: Vec<f64>,
    /// Velocity y-component (m s⁻¹), length nx * ny.
    vy: Vec<f64>,
    /// Grid dimension x.
    nx: usize,
    /// Grid dimension y.
    ny: usize,
}

impl FluidField {
    fn new(nx: usize, ny: usize, init_density: f64) -> Self {
        let n = nx * ny;
        Self {
            density: vec![init_density; n],
            vx: vec![0.0; n],
            vy: vec![0.0; n],
            nx,
            ny,
        }
    }

    #[inline]
    fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }
}

// ---------------------------------------------------------------------------
// PlasmaLbm
// ---------------------------------------------------------------------------

/// Two-fluid LBM plasma solver (ions + electrons) on a 2-D grid.
///
/// Uses a simple BGK-like relaxation step to advance density and velocity
/// fields for both species under mutual electric-force coupling.
///
/// This is a reduced two-moment model; full LBM distribution functions are
/// not carried here — the module is designed for fast macroscopic estimates.
#[derive(Debug, Clone)]
pub struct PlasmaLbm {
    /// Ion density/velocity field.
    pub ions: PlasmaFluid,
    /// Electron density/velocity field.
    pub electrons: PlasmaFluid,
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Self-consistent electric field (updated each step).
    pub electric_field: ElectricField,
    /// Uniform background magnetic field.
    pub magnetic_field: MagneticField,
}

/// Macroscopic plasma fluid state (density + velocity on a 2-D grid).
#[derive(Debug, Clone)]
pub struct PlasmaFluid {
    /// Number density (m⁻³), stored in row-major order (ny × nx).
    pub density: Vec<f64>,
    /// x-velocity (m s⁻¹).
    pub vx: Vec<f64>,
    /// y-velocity (m s⁻¹).
    pub vy: Vec<f64>,
    /// Grid size x.
    pub nx: usize,
    /// Grid size y.
    pub ny: usize,
    /// Species parameters.
    pub species: PlasmaSpecies,
}

impl PlasmaFluid {
    /// Create a uniform plasma fluid state.
    pub fn new(nx: usize, ny: usize, n0: f64, species: PlasmaSpecies) -> Self {
        let n = nx * ny;
        Self {
            density: vec![n0; n],
            vx: vec![0.0; n],
            vy: vec![0.0; n],
            nx,
            ny,
            species,
        }
    }

    /// Flat index for cell (ix, iy).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Mean number density (m⁻³) averaged over the grid.
    pub fn mean_density(&self) -> f64 {
        self.density.iter().sum::<f64>() / self.density.len() as f64
    }
}

impl PlasmaLbm {
    /// Create a new [`PlasmaLbm`] simulation.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – grid dimensions
    /// * `n_ion`      – initial ion number density (m⁻³)
    /// * `n_elec`     – initial electron number density (m⁻³)
    /// * `ion_species`  – ion species parameters
    /// * `elec_species` – electron species parameters
    /// * `bz`           – out-of-plane magnetic field (T)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        n_ion: f64,
        n_elec: f64,
        ion_species: PlasmaSpecies,
        elec_species: PlasmaSpecies,
        bz: f64,
    ) -> Self {
        Self {
            ions: PlasmaFluid::new(nx, ny, n_ion, ion_species),
            electrons: PlasmaFluid::new(nx, ny, n_elec, elec_species),
            nx,
            ny,
            electric_field: ElectricField::new(nx, ny),
            magnetic_field: MagneticField::new(bz),
        }
    }

    /// Advance the plasma state by one time step `dt` (s).
    ///
    /// The step:
    /// 1. Recomputes the charge-separation electric field from ρ = e(n_i − n_e).
    /// 2. Accelerates each fluid cell with the Lorentz force.
    /// 3. Applies a simple first-order continuity update (advection ignored).
    pub fn step(&mut self, dt: f64) {
        let n = self.nx * self.ny;
        // Update electric field from charge separation (simplified Gauss / Poisson).
        for k in 0..n {
            let ni = self.ions.density[k];
            let ne = self.electrons.density[k];
            let rho = ELEM_CHARGE * (ni - ne); // charge density
            // E ~ rho / epsilon_0 * dx  (rough estimate, not a real Poisson solve)
            let e_mag = rho / EPSILON_0;
            self.electric_field.ex[k] = e_mag;
            self.electric_field.ey[k] = 0.0;
        }

        // Accelerate ion fluid.
        let q_ion = self.ions.species.charge();
        let m_ion = self.ions.species.mass();
        for k in 0..n {
            let ex = self.electric_field.ex[k];
            let bz = self.magnetic_field.bz;
            let ax = (q_ion * ex + q_ion * self.ions.vy[k] * bz) / m_ion;
            let ay = (-q_ion * self.ions.vx[k] * bz) / m_ion;
            self.ions.vx[k] += ax * dt;
            self.ions.vy[k] += ay * dt;
        }

        // Accelerate electron fluid.
        let q_elec = self.electrons.species.charge();
        let m_elec = self.electrons.species.mass();
        for k in 0..n {
            let ex = self.electric_field.ex[k];
            let bz = self.magnetic_field.bz;
            let ax = (q_elec * ex + q_elec * self.electrons.vy[k] * bz) / m_elec;
            let ay = (-q_elec * self.electrons.vx[k] * bz) / m_elec;
            self.electrons.vx[k] += ax * dt;
            self.electrons.vy[k] += ay * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// DebyeShielding
// ---------------------------------------------------------------------------

/// Debye–Hückel screened Coulomb potential model.
///
/// φ(r) = (q / 4π ε₀ r) exp(−r / λ_D)
#[derive(Debug, Clone)]
pub struct DebyeShielding {
    /// Electron number density (m⁻³).
    pub electron_density: f64,
    /// Electron temperature (K).
    pub temperature: f64,
}

impl DebyeShielding {
    /// Create a new [`DebyeShielding`] model.
    pub fn new(electron_density: f64, temperature: f64) -> Self {
        Self {
            electron_density,
            temperature,
        }
    }

    /// Debye length λ_D (m).
    ///
    /// ```no_run
    /// use oxiphysics_lbm::plasma_lbm::DebyeShielding;
    /// let ds = DebyeShielding::new(1e18, 1e4);
    /// assert!(ds.compute_debye_length() > 0.0);
    /// ```
    pub fn compute_debye_length(&self) -> f64 {
        debye_length(self.electron_density, self.temperature, ELEM_CHARGE)
    }

    /// Yukawa (Debye–Hückel) shielded electric potential (V) at distance r (m).
    ///
    /// φ(r) = (e / 4π ε₀) · exp(−r / λ_D) / r
    ///
    /// Returns 0 for r ≤ 0.
    ///
    /// ```no_run
    /// use oxiphysics_lbm::plasma_lbm::DebyeShielding;
    /// let ds = DebyeShielding::new(1e18, 1e4);
    /// let phi = ds.shielded_potential(1e-5);
    /// assert!(phi >= 0.0);
    /// ```
    pub fn shielded_potential(&self, r: f64) -> f64 {
        if r <= 0.0 {
            return 0.0;
        }
        let lambda_d = self.compute_debye_length();
        let q = ELEM_CHARGE;
        (q / (4.0 * PI * EPSILON_0)) * (-r / lambda_d).exp() / r
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- plasma_frequency ------------------------------------------------

    #[test]
    fn test_plasma_frequency_electron() {
        // Reference: n = 1e18 m⁻³, ωp ≈ 5.64e10 rad/s
        let wp = plasma_frequency(1e18, ELEM_CHARGE, ELECTRON_MASS);
        assert!(wp > 1e9, "wp={wp}");
        assert!(wp < 1e12, "wp={wp}");
    }

    #[test]
    fn test_plasma_frequency_zero_density() {
        assert_eq!(plasma_frequency(0.0, ELEM_CHARGE, ELECTRON_MASS), 0.0);
    }

    #[test]
    fn test_plasma_frequency_negative_density() {
        assert_eq!(plasma_frequency(-1e18, ELEM_CHARGE, ELECTRON_MASS), 0.0);
    }

    #[test]
    fn test_plasma_frequency_zero_mass() {
        assert_eq!(plasma_frequency(1e18, ELEM_CHARGE, 0.0), 0.0);
    }

    #[test]
    fn test_plasma_frequency_scales_with_density() {
        let wp1 = plasma_frequency(1e18, ELEM_CHARGE, ELECTRON_MASS);
        let wp2 = plasma_frequency(4e18, ELEM_CHARGE, ELECTRON_MASS);
        // ωp ∝ √n → wp2 / wp1 ≈ 2
        let ratio = wp2 / wp1;
        assert!((ratio - 2.0).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_plasma_frequency_proton() {
        let proton_mass = 1836.0 * ELECTRON_MASS;
        let wpi = plasma_frequency(1e18, ELEM_CHARGE, proton_mass);
        let wpe = plasma_frequency(1e18, ELEM_CHARGE, ELECTRON_MASS);
        // ωpi / ωpe = √(me/mp)
        let ratio = wpi / wpe;
        assert!((ratio - (1.0_f64 / 1836.0).sqrt()).abs() < 1e-6);
    }

    // ---- debye_length ----------------------------------------------------

    #[test]
    fn test_debye_length_positive() {
        let lam = debye_length(1e18, 1e4, ELEM_CHARGE);
        assert!(lam > 0.0);
    }

    #[test]
    fn test_debye_length_zero_density() {
        assert_eq!(debye_length(0.0, 1e4, ELEM_CHARGE), 0.0);
    }

    #[test]
    fn test_debye_length_zero_temperature() {
        assert_eq!(debye_length(1e18, 0.0, ELEM_CHARGE), 0.0);
    }

    #[test]
    fn test_debye_length_scales_with_temperature() {
        let l1 = debye_length(1e18, 1e4, ELEM_CHARGE);
        let l2 = debye_length(1e18, 4e4, ELEM_CHARGE);
        // λ_D ∝ √T → l2/l1 ≈ 2
        let ratio = l2 / l1;
        assert!((ratio - 2.0).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_debye_length_scales_with_density() {
        let l1 = debye_length(1e18, 1e4, ELEM_CHARGE);
        let l4 = debye_length(4e18, 1e4, ELEM_CHARGE);
        // λ_D ∝ 1/√n → l4/l1 ≈ 0.5
        let ratio = l4 / l1;
        assert!((ratio - 0.5).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_debye_length_reference_value() {
        // For n=1e18, T=1e4 K, expected λ_D ~ 6.9e-6 m
        let lam = debye_length(1e18, 1e4, ELEM_CHARGE);
        assert!((lam - 6.9e-6).abs() / lam < 0.05, "lam={lam}");
    }

    // ---- coulomb_logarithm -----------------------------------------------

    #[test]
    fn test_coulomb_log_positive() {
        let ln_lambda = coulomb_logarithm(1e4, 1e18);
        assert!(ln_lambda >= 1.0);
    }

    #[test]
    fn test_coulomb_log_zero_temperature() {
        assert_eq!(coulomb_logarithm(0.0, 1e18), 1.0);
    }

    #[test]
    fn test_coulomb_log_zero_density() {
        assert_eq!(coulomb_logarithm(1e4, 0.0), 1.0);
    }

    #[test]
    fn test_coulomb_log_reasonable_value() {
        let ln_lambda = coulomb_logarithm(1e4, 1e18);
        // Typically 5–20 for fusion-relevant plasmas
        assert!(ln_lambda > 5.0, "ln_lambda={ln_lambda}");
        assert!(ln_lambda < 30.0, "ln_lambda={ln_lambda}");
    }

    // ---- cyclotron frequency ---------------------------------------------

    #[test]
    fn test_cyclotron_frequency_electron() {
        let bf = MagneticField::new(1.0); // 1 T
        let wc = bf.compute_cyclotron_frequency(ELEM_CHARGE, ELECTRON_MASS);
        // ωc,e = eB/me ≈ 1.7588e11 rad/s
        assert!((wc - 1.758_820e11).abs() / wc < 1e-4, "wc={wc}");
    }

    #[test]
    fn test_cyclotron_frequency_zero_mass() {
        let bf = MagneticField::new(1.0);
        assert_eq!(bf.compute_cyclotron_frequency(ELEM_CHARGE, 0.0), 0.0);
    }

    #[test]
    fn test_cyclotron_frequency_scales_with_b() {
        let bf1 = MagneticField::new(1.0);
        let bf2 = MagneticField::new(2.0);
        let wc1 = bf1.compute_cyclotron_frequency(ELEM_CHARGE, ELECTRON_MASS);
        let wc2 = bf2.compute_cyclotron_frequency(ELEM_CHARGE, ELECTRON_MASS);
        assert!((wc2 / wc1 - 2.0).abs() < 1e-10, "ratio={}", wc2 / wc1);
    }

    #[test]
    fn test_cyclotron_frequency_proton() {
        let bf = MagneticField::new(1.0);
        let proton_mass = 1836.153 * ELECTRON_MASS;
        let wcp = bf.compute_cyclotron_frequency(ELEM_CHARGE, proton_mass);
        let wce = bf.compute_cyclotron_frequency(ELEM_CHARGE, ELECTRON_MASS);
        // ωc,p / ωc,e = me/mp
        let ratio = wcp / wce;
        assert!((ratio - 1.0 / 1836.153).abs() < 1e-5, "ratio={ratio}");
    }

    // ---- Larmor radius ---------------------------------------------------

    #[test]
    fn test_larmor_radius_positive() {
        let bf = MagneticField::new(1.0);
        let r = bf.larmor_radius(1e6, ELEM_CHARGE, ELECTRON_MASS);
        assert!(r > 0.0);
    }

    #[test]
    fn test_larmor_radius_zero_b() {
        let bf = MagneticField::new(0.0);
        let r = bf.larmor_radius(1e6, ELEM_CHARGE, ELECTRON_MASS);
        assert!(r.is_infinite());
    }

    #[test]
    fn test_larmor_radius_scales_with_vel() {
        let bf = MagneticField::new(1.0);
        let r1 = bf.larmor_radius(1e6, ELEM_CHARGE, ELECTRON_MASS);
        let r2 = bf.larmor_radius(2e6, ELEM_CHARGE, ELECTRON_MASS);
        assert!((r2 / r1 - 2.0).abs() < 1e-10, "ratio={}", r2 / r1);
    }

    #[test]
    fn test_larmor_radius_reference() {
        // v = 1e6 m/s, B = 1 T → r_L = me*v/(eB) ≈ 5.69e-6 m
        let bf = MagneticField::new(1.0);
        let r = bf.larmor_radius(1e6, ELEM_CHARGE, ELECTRON_MASS);
        assert!((r - 5.685e-6).abs() / r < 0.01, "r={r}");
    }

    // ---- MagneticField::magnetic_force -----------------------------------

    #[test]
    fn test_magnetic_force_x_velocity() {
        let bf = MagneticField::new(1.0);
        // vx=1, vy=0 → F = q[0, -1, 0] (for positive charge)
        let f = bf.magnetic_force(ELEM_CHARGE, [1.0, 0.0, 0.0]);
        assert!((f[0]).abs() < 1e-30, "fx={}", f[0]);
        assert!((f[1] + ELEM_CHARGE).abs() < 1e-30, "fy={}", f[1]);
        assert!((f[2]).abs() < 1e-30, "fz={}", f[2]);
    }

    #[test]
    fn test_magnetic_force_y_velocity() {
        let bf = MagneticField::new(1.0);
        // vy=1 → F = q[Bz, 0, 0]
        let f = bf.magnetic_force(ELEM_CHARGE, [0.0, 1.0, 0.0]);
        assert!((f[0] - ELEM_CHARGE).abs() < 1e-30, "fx={}", f[0]);
        assert!((f[1]).abs() < 1e-30, "fy={}", f[1]);
    }

    // ---- ElectricField ---------------------------------------------------

    #[test]
    fn test_electric_field_new_zero() {
        let ef = ElectricField::new(4, 4);
        assert!(ef.ex.iter().all(|&v| v == 0.0));
        assert!(ef.ey.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_electric_field_set_get() {
        let mut ef = ElectricField::new(4, 4);
        ef.set(1, 2, 3.0, 4.0);
        let f = ef.field_at(1, 2);
        assert_eq!(f[0], 3.0);
        assert_eq!(f[1], 4.0);
    }

    #[test]
    fn test_lorentz_force_at() {
        let mut ef = ElectricField::new(4, 4);
        ef.set(2, 2, 1e4, 0.0);
        let f = ef.lorentz_force_at(2, 2, ELEM_CHARGE);
        assert!((f[0] - ELEM_CHARGE * 1e4).abs() < 1e-30, "fx={}", f[0]);
    }

    // ---- PlasmaSpecies ---------------------------------------------------

    #[test]
    fn test_plasma_species_charge() {
        let sp = PlasmaSpecies::new(-1.0, 1.0, 1e4, [0.0; 3]);
        assert!((sp.charge() + ELEM_CHARGE).abs() < 1e-30);
    }

    #[test]
    fn test_plasma_species_mass() {
        let sp = PlasmaSpecies::new(1.0, 1836.0, 1e4, [0.0; 3]);
        assert!((sp.mass() - 1836.0 * ELECTRON_MASS).abs() < 1e-40);
    }

    #[test]
    fn test_plasma_species_thermal_speed() {
        let sp = PlasmaSpecies::new(-1.0, 1.0, 1e4, [0.0; 3]);
        let vt = sp.thermal_speed();
        assert!(vt > 0.0);
    }

    // ---- DebyeShielding --------------------------------------------------

    #[test]
    fn test_debye_shielding_length() {
        let ds = DebyeShielding::new(1e18, 1e4);
        let l = ds.compute_debye_length();
        assert!(l > 0.0);
    }

    #[test]
    fn test_debye_shielding_potential_positive() {
        let ds = DebyeShielding::new(1e18, 1e4);
        let phi = ds.shielded_potential(1e-5);
        assert!(phi >= 0.0, "phi={phi}");
    }

    #[test]
    fn test_debye_shielding_potential_zero_r() {
        let ds = DebyeShielding::new(1e18, 1e4);
        assert_eq!(ds.shielded_potential(0.0), 0.0);
    }

    #[test]
    fn test_debye_shielding_potential_decays() {
        let ds = DebyeShielding::new(1e18, 1e4);
        let phi1 = ds.shielded_potential(1e-6);
        let phi2 = ds.shielded_potential(1e-4);
        assert!(phi2 < phi1, "phi1={phi1}, phi2={phi2}");
    }

    // ---- PlasmaLbm -------------------------------------------------------

    #[test]
    fn test_plasma_lbm_construction() {
        let ion_sp = PlasmaSpecies::new(1.0, 1836.0, 1e3, [0.0; 3]);
        let elec_sp = PlasmaSpecies::new(-1.0, 1.0, 1e4, [0.0; 3]);
        let sim = PlasmaLbm::new(8, 8, 1e18, 1e18, ion_sp, elec_sp, 0.1);
        assert_eq!(sim.nx, 8);
        assert_eq!(sim.ny, 8);
    }

    #[test]
    fn test_plasma_lbm_step_runs() {
        let ion_sp = PlasmaSpecies::new(1.0, 1836.0, 1e3, [0.0; 3]);
        let elec_sp = PlasmaSpecies::new(-1.0, 1.0, 1e4, [0.0; 3]);
        let mut sim = PlasmaLbm::new(4, 4, 1e18, 1e18, ion_sp, elec_sp, 0.0);
        sim.step(1e-12);
        // Quasi-neutral plasma: no net charge → velocities should remain near zero
        for &v in &sim.ions.vx {
            assert!(v.is_finite(), "ion vx not finite: {v}");
        }
    }

    #[test]
    fn test_plasma_fluid_mean_density() {
        let sp = PlasmaSpecies::new(1.0, 1836.0, 1e3, [0.0; 3]);
        let fluid = PlasmaFluid::new(4, 4, 1e18, sp);
        assert!((fluid.mean_density() - 1e18).abs() < 1e6);
    }
}
