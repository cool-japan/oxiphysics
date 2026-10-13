// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Weakly ionized plasma simulation using the Lattice Boltzmann Method.
//!
//! This module provides:
//! - [`PlasmaSpecies`]: Parameters for a plasma species (ions or electrons).
//! - [`ElectricField`]: 2-D electric field on a grid with Lorentz force computation.
//! - [`MagneticField`]: Uniform out-of-plane magnetic field with cyclotron and Larmor formulae.
//! - [`PlasmaLbm`]: Two-fluid (ion + electron) driver with a spectral Poisson solve,
//!   a Boris push and an upwind finite-volume continuity update.
//! - [`poisson_electric_field_periodic`]: Electrostatic field of a periodic charge density.
//! - [`boris_push_2d`]: Boris integration of the Lorentz force.
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

use crate::{Error, Result};
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
// PlasmaLbm
// ---------------------------------------------------------------------------

/// Largest Courant number (sum over both directions, counting the in- and
/// outflow face of a cell) for which one upwind sub-step of the continuity
/// equation is taken; larger time steps are sub-cycled.
const CONTINUITY_CFL: f64 = 0.9;

/// Largest number of continuity sub-steps taken for one plasma step.
const MAX_CONTINUITY_SUBSTEPS: usize = 10_000;

/// Fallback grid spacing (m) when no electron Debye length is available.
const FALLBACK_GRID_SPACING: f64 = 1.0e-6;

/// Wavenumbers of FFT bin `index` on an `n`-point periodic axis with grid
/// spacing `dx`: `(k, k_derivative)`.
///
/// `k` is the plain wavenumber `2 pi m / (n dx)` with `m` the signed mode
/// number.  `k_derivative` is the wavenumber to use when differentiating a
/// real field: it is zero for the Nyquist mode of an even-length axis, whose
/// derivative is not representable by a real field.
fn fft_wavenumbers(index: usize, n: usize, dx: f64) -> (f64, f64) {
    let signed = if index <= n / 2 {
        index as f64
    } else {
        index as f64 - n as f64
    };
    let k = 2.0 * PI * signed / (n as f64 * dx);
    let k_derivative = if n % 2 == 0 && index == n / 2 { 0.0 } else { k };
    (k, k_derivative)
}

/// Electrostatic field `E = -∇φ` of a charge density on a periodic grid.
///
/// Solves Gauss's law `∇²φ = -ρ_c / ε₀` spectrally: the charge density is
/// transformed with a 2-D FFT, divided by `ε₀ k²` (the `k = 0` mode is dropped,
/// i.e. the net charge is assumed to vanish or to be compensated by a uniform
/// background), differentiated in Fourier space and transformed back.  The
/// result is exact for every Fourier mode that fits on the grid.
///
/// `charge_density` is stored in row-major order (index `iy * nx + ix`) in
/// C m⁻³, `dx` is the grid spacing (m) in both directions.  Returns
/// `(E_x, E_y)` in V m⁻¹.
///
/// # Errors
/// Empty grids, a length mismatch, a non-positive or non-finite `dx`, or
/// non-finite charge values.
pub fn poisson_electric_field_periodic(
    charge_density: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Result<(Vec<f64>, Vec<f64>)> {
    Error::check_grid(nx, ny)?;
    let n = nx * ny;
    if charge_density.len() != n {
        return Err(Error::ComponentCountMismatch {
            expected: n,
            found: charge_density.len(),
            field: "charge_density".to_string(),
        });
    }
    if !(dx.is_finite() && dx > 0.0) {
        return Err(Error::ParameterOutOfBounds {
            name: "dx".to_string(),
            value: dx,
            min: f64::MIN_POSITIVE,
            max: f64::MAX,
        });
    }
    for (index, &q) in charge_density.iter().enumerate() {
        Error::check_finite(q, "charge_density", index)?;
    }

    let zeros = vec![0.0_f64; n];
    // Row-major (iy, ix): `ny` rows of `nx` columns.
    let (rho_re, rho_im) = oxifft::fft2d_split::<f64>(charge_density, &zeros, ny, nx);

    let mut ex_re = vec![0.0_f64; n];
    let mut ex_im = vec![0.0_f64; n];
    let mut ey_re = vec![0.0_f64; n];
    let mut ey_im = vec![0.0_f64; n];
    for iy in 0..ny {
        let (ky, ky_d) = fft_wavenumbers(iy, ny, dx);
        for ix in 0..nx {
            let (kx, kx_d) = fft_wavenumbers(ix, nx, dx);
            let k2 = kx * kx + ky * ky;
            if k2 == 0.0 {
                continue;
            }
            let idx = iy * nx + ix;
            // φ̂ = ρ̂ / (ε₀ k²)
            let scale = 1.0 / (EPSILON_0 * k2);
            let phi_re = rho_re[idx] * scale;
            let phi_im = rho_im[idx] * scale;
            // Ê = -i k φ̂
            ex_re[idx] = kx_d * phi_im;
            ex_im[idx] = -kx_d * phi_re;
            ey_re[idx] = ky_d * phi_im;
            ey_im[idx] = -ky_d * phi_re;
        }
    }
    let (ex, _) = oxifft::ifft2d_split::<f64>(&ex_re, &ex_im, ny, nx);
    let (ey, _) = oxifft::ifft2d_split::<f64>(&ey_re, &ey_im, ny, nx);
    Ok((ex, ey))
}

/// Boris push of an in-plane velocity in an electric field `e` and a uniform
/// out-of-plane magnetic field `B = (0, 0, bz)`.
///
/// Integrates `m dv/dt = q (E + v × B)` over `dt` with the Boris scheme: half an
/// electric kick, an exact rotation about `B` by the angle `2 atan(q B dt / 2m)`
/// and a second half kick.  Without an electric field the speed `|v|` is
/// conserved to rounding error for any `dt` (an explicit Euler step multiplies
/// it by `sqrt(1 + (ω_c dt)²)` every step), and a particle in crossed fields
/// keeps its `E × B / B²` drift exactly.
///
/// # Arguments
/// * `v`      – velocity `[vx, vy]` (m s⁻¹)
/// * `e`      – electric field `[Ex, Ey]` (V m⁻¹)
/// * `bz`     – magnetic flux density B_z (T)
/// * `charge` – particle charge (C)
/// * `mass`   – particle mass (kg)
/// * `dt`     – time step (s)
pub fn boris_push_2d(
    v: [f64; 2],
    e: [f64; 2],
    bz: f64,
    charge: f64,
    mass: f64,
    dt: f64,
) -> [f64; 2] {
    let q_over_m = charge / mass;
    let kick = [q_over_m * e[0] * 0.5 * dt, q_over_m * e[1] * 0.5 * dt];
    // First half electric kick.
    let v_minus = [v[0] + kick[0], v[1] + kick[1]];
    // Magnetic rotation about z: v' = v⁻ + v⁻ × t,  v⁺ = v⁻ + v' × s.
    let t = q_over_m * bz * 0.5 * dt;
    let s = 2.0 * t / (1.0 + t * t);
    let v_prime = [v_minus[0] + v_minus[1] * t, v_minus[1] - v_minus[0] * t];
    let v_plus = [v_minus[0] + v_prime[1] * s, v_minus[1] - v_prime[0] * s];
    // Second half electric kick.
    [v_plus[0] + kick[0], v_plus[1] + kick[1]]
}

/// Advance a density field by `dt` with the conservative first-order upwind
/// finite-volume discretisation of `∂n/∂t + ∇·(n v) = 0` on a periodic grid.
///
/// The face velocity is the average of the two adjacent cell velocities and the
/// density is taken from the upwind cell.  The scheme conserves the total
/// particle number exactly and keeps the density non-negative; time steps
/// that violate the CFL condition are split into sub-steps.
fn advance_density_upwind(
    density: &mut [f64],
    vx: &[f64],
    vy: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dt: f64,
) -> Result<()> {
    let n = nx * ny;
    let east = |ix: usize| (ix + 1) % nx;
    let north = |iy: usize| (iy + 1) % ny;

    let mut max_face_x = 0.0_f64;
    let mut max_face_y = 0.0_f64;
    for iy in 0..ny {
        for ix in 0..nx {
            let k = iy * nx + ix;
            max_face_x = max_face_x.max((0.5 * (vx[k] + vx[iy * nx + east(ix)])).abs());
            max_face_y = max_face_y.max((0.5 * (vy[k] + vy[north(iy) * nx + ix])).abs());
        }
    }
    // Positivity needs dt * (sum of the outgoing face speeds) / dx <= 1 and a
    // cell can lose through both faces of each direction.
    let courant = 2.0 * dt * (max_face_x + max_face_y) / dx;
    let n_sub_f = (courant / CONTINUITY_CFL).ceil().max(1.0);
    if !n_sub_f.is_finite() || n_sub_f > MAX_CONTINUITY_SUBSTEPS as f64 {
        return Err(Error::TimeStepTooLarge {
            dt,
            dt_max: dt * CONTINUITY_CFL * MAX_CONTINUITY_SUBSTEPS as f64 / courant,
        });
    }
    let n_sub = n_sub_f as usize;
    let h = dt / n_sub_f;

    let mut flux_x = vec![0.0_f64; n];
    let mut flux_y = vec![0.0_f64; n];
    for _ in 0..n_sub {
        // Flux through the east / north face of every cell.
        for iy in 0..ny {
            for ix in 0..nx {
                let k = iy * nx + ix;
                let ke = iy * nx + east(ix);
                let kn = north(iy) * nx + ix;
                let ux = 0.5 * (vx[k] + vx[ke]);
                flux_x[k] = ux * if ux >= 0.0 { density[k] } else { density[ke] };
                let uy = 0.5 * (vy[k] + vy[kn]);
                flux_y[k] = uy * if uy >= 0.0 { density[k] } else { density[kn] };
            }
        }
        for iy in 0..ny {
            for ix in 0..nx {
                let k = iy * nx + ix;
                let kw = iy * nx + (ix + nx - 1) % nx;
                let ks = ((iy + ny - 1) % ny) * nx + ix;
                density[k] -= h / dx * (flux_x[k] - flux_x[kw] + flux_y[k] - flux_y[ks]);
            }
        }
    }
    Ok(())
}

/// Two-fluid plasma solver (ions + electrons) on a periodic 2-D grid.
///
/// Each species carries a number density and an in-plane velocity field; the
/// grid has the uniform spacing [`PlasmaLbm::dx`].  One step of
/// [`PlasmaLbm::step`] / [`PlasmaLbm::try_step`]
///
/// 1. solves Gauss's law `∇²φ = -ρ_c/ε₀` spectrally for the self-consistent
///    electrostatic field `E = -∇φ` of the charge density
///    `ρ_c = q_i n_i + q_e n_e` ([`poisson_electric_field_periodic`]),
/// 2. accelerates both species with the Lorentz force `q (E + v × B)` using the
///    Boris scheme ([`boris_push_2d`]), which rotates the velocity about the
///    uniform out-of-plane field exactly instead of letting `|v|` grow,
/// 3. advances the number densities with the continuity equation
///    `∂n/∂t + ∇·(n v) = 0` (conservative first-order upwind finite volumes,
///    sub-cycled to respect the CFL condition).
///
/// Updating the velocity before the densities makes the coupled
/// density-velocity oscillation (the Langmuir wave) stable and, for
/// `ω_p dt ≲ 0.2`, accurate to second order in `dt`.
///
/// This is a cold, collisionless, reduced two-moment model: pressure,
/// collision and viscous terms are not included, and no LBM distribution
/// functions are carried — the module is designed for fast macroscopic
/// estimates.  The net charge of the periodic domain should vanish (the
/// uniform part of the charge density does not produce a field).
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
    /// Grid spacing in x and y (m).
    pub dx: f64,
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

    /// Check that the field arrays match the grid and the species can be
    /// integrated (positive finite mass, finite charge and fields).
    fn validate(&self, nx: usize, ny: usize, name: &str) -> Result<()> {
        let n = nx * ny;
        for (field, values) in [
            ("density", &self.density),
            ("vx", &self.vx),
            ("vy", &self.vy),
        ] {
            if values.len() != n {
                return Err(Error::ComponentCountMismatch {
                    expected: n,
                    found: values.len(),
                    field: format!("{name}.{field}"),
                });
            }
            for (index, &value) in values.iter().enumerate() {
                Error::check_finite(value, &format!("{name}.{field}"), index)?;
            }
        }
        let mass = self.species.mass();
        if !(mass.is_finite() && mass > 0.0) {
            return Err(Error::ParameterOutOfBounds {
                name: format!("{name}.mass"),
                value: mass,
                min: f64::MIN_POSITIVE,
                max: f64::MAX,
            });
        }
        Error::check_finite(self.species.charge(), &format!("{name}.charge"), 0)
    }
}

impl PlasmaLbm {
    /// Create a new [`PlasmaLbm`] simulation.
    ///
    /// The grid spacing [`PlasmaLbm::dx`] defaults to the electron Debye length
    /// of the initial electron density and temperature (1 µm if that is
    /// undefined), so the default grid resolves charge separation; use
    /// [`PlasmaLbm::with_grid_spacing`] to choose another one.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – grid dimensions
    /// * `n_ion`      – initial ion number density (m⁻³)
    /// * `n_elec`     – initial electron number density (m⁻³)
    /// * `ion_species`  – ion species parameters
    /// * `elec_species` – electron species parameters
    /// * `bz`           – out-of-plane magnetic field (T)
    pub fn new(
        nx: usize,
        ny: usize,
        n_ion: f64,
        n_elec: f64,
        ion_species: PlasmaSpecies,
        elec_species: PlasmaSpecies,
        bz: f64,
    ) -> Self {
        let debye = debye_length(n_elec, elec_species.temperature, ELEM_CHARGE);
        let dx = if debye.is_finite() && debye > 0.0 {
            debye
        } else {
            FALLBACK_GRID_SPACING
        };
        Self {
            ions: PlasmaFluid::new(nx, ny, n_ion, ion_species),
            electrons: PlasmaFluid::new(nx, ny, n_elec, elec_species),
            nx,
            ny,
            dx,
            electric_field: ElectricField::new(nx, ny),
            magnetic_field: MagneticField::new(bz),
        }
    }

    /// Set the grid spacing (m) in x and y.
    pub fn with_grid_spacing(mut self, dx: f64) -> Self {
        self.dx = dx;
        self
    }

    /// Charge density `ρ_c = q_i n_i + q_e n_e` (C m⁻³) of the current state.
    pub fn charge_density(&self) -> Vec<f64> {
        let q_ion = self.ions.species.charge();
        let q_elec = self.electrons.species.charge();
        self.ions
            .density
            .iter()
            .zip(self.electrons.density.iter())
            .map(|(&ni, &ne)| q_ion * ni + q_elec * ne)
            .collect()
    }

    /// Advance the plasma state by one time step `dt` (s).
    ///
    /// See [`PlasmaLbm::try_step`] for the algorithm.  If the inputs are
    /// invalid (non-positive or non-finite `dt` or `dx`, inconsistent field
    /// sizes, non-finite fields, a time step too large for the sub-cycling
    /// limit) the state is left unchanged; use [`PlasmaLbm::try_step`] to learn
    /// why.
    pub fn step(&mut self, dt: f64) {
        let _ = self.try_step(dt);
    }

    /// Advance the plasma state by one time step `dt` (s), reporting invalid
    /// input.
    ///
    /// The step:
    /// 1. Solves the Poisson equation for the electrostatic field of the charge
    ///    density `ρ_c = q_i n_i + q_e n_e` on the periodic grid and stores it in
    ///    [`PlasmaLbm::electric_field`].
    /// 2. Accelerates each species with the Lorentz force `q (E + v × B)` (Boris
    ///    scheme).
    /// 3. Advances both densities with the continuity equation using the
    ///    updated velocities (upwind finite volumes, total particle number
    ///    conserved).
    ///
    /// # Errors
    /// Returns an error — and leaves the state untouched — if `dt` or `dx` is
    /// not positive and finite, the grid is empty, the field arrays do not match
    /// the grid or contain non-finite values, a species has a non-positive mass,
    /// or the continuity step would need more than 10 000 CFL sub-steps.
    pub fn try_step(&mut self, dt: f64) -> Result<()> {
        Error::check_grid(self.nx, self.ny)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(Error::ParameterOutOfBounds {
                name: "dt".to_string(),
                value: dt,
                min: f64::MIN_POSITIVE,
                max: f64::MAX,
            });
        }
        self.ions.validate(self.nx, self.ny, "ions")?;
        self.electrons.validate(self.nx, self.ny, "electrons")?;

        // 1. Self-consistent electrostatic field.
        let charge = self.charge_density();
        let (ex, ey) = poisson_electric_field_periodic(&charge, self.nx, self.ny, self.dx)?;

        // Work on copies so that a failing continuity step leaves no trace.
        let mut ions = self.ions.clone();
        let mut electrons = self.electrons.clone();
        let bz = self.magnetic_field.bz;
        for fluid in [&mut ions, &mut electrons] {
            let q = fluid.species.charge();
            let m = fluid.species.mass();
            // 2. Boris push of every cell.
            for k in 0..ex.len() {
                let v = boris_push_2d([fluid.vx[k], fluid.vy[k]], [ex[k], ey[k]], bz, q, m, dt);
                fluid.vx[k] = v[0];
                fluid.vy[k] = v[1];
            }
            // 3. Continuity with the updated velocities.
            advance_density_upwind(
                &mut fluid.density,
                &fluid.vx,
                &fluid.vy,
                self.nx,
                self.ny,
                self.dx,
                dt,
            )?;
        }

        self.ions = ions;
        self.electrons = electrons;
        self.electric_field.ex = ex;
        self.electric_field.ey = ey;
        Ok(())
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
mod field_tests;

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
