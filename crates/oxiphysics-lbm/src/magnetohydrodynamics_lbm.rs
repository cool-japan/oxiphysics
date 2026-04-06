// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Magnetohydrodynamics Lattice Boltzmann Method (MHD-LBM).
//!
//! This module provides:
//! - [`MhdState`]: Full MHD state (density, velocity, magnetic field, pressure).
//! - [`MhdLbm`]: Dual-distribution MHD-LBM solver on a D2Q9 lattice.
//! - [`LorentzForce`]: Computes the J×B Lorentz body force.
//! - [`InductionEquation`]: Magnetic field advection-diffusion solver.
//! - [`AlfvenWave`]: Alfvén wave speed and dispersion.
//! - [`MagneticReconnection`]: Sweet-Parker reconnection diagnostics.
//! - [`MhdEquilibrium`]: Force-free equilibrium and Grad-Shafranov solver.
//! - [`MhdAnalysis`]: Dimensionless MHD numbers and plasma diagnostics.
//!
//! # Physics Background
//!
//! The MHD equations couple the Navier-Stokes equations to Maxwell's
//! equations in the low-frequency, non-relativistic limit:
//!
//! ∂ρ/∂t + ∇·(ρu) = 0
//! ρ(∂u/∂t + u·∇u) = −∇p + J×B + ν∇²u
//! ∂B/∂t = ∇×(u×B) + η∇²B
//! ∇·B = 0
//!
//! where J = ∇×B/μ₀ is the current density, η the magnetic diffusivity,
//! and ν the kinematic viscosity.
//!
//! The dual-distribution LBM uses two sets of distribution functions:
//! - f_i: standard fluid distributions recovering mass and momentum
//! - g_i: magnetic distributions recovering the magnetic field evolution
//!
//! # References
//! - Dellar, P.J. (2002). Lattice kinetic schemes for magnetohydrodynamics.
//!   *J. Comput. Phys.*, 179, 95–126.
//! - Succi, S., Vergassola, M., & Benzi, R. (1991). Lattice Boltzmann scheme
//!   for two-dimensional magnetohydrodynamics. *Phys. Rev. A*, 43, 4521.
//! - Sweet, P.A. (1958). The neutral point theory of solar flares.
//!   *IAU Symp.*, 6, 123–134.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Magnetic permeability of free space μ₀ (H m⁻¹).
pub const MU_0: f64 = 1.256_637_061_4e-6;

/// Speed of light in vacuum c (m s⁻¹).
pub const SPEED_OF_LIGHT: f64 = 2.997_924_58e8;

/// Number of discrete velocities in D2Q9 lattice.
pub const D2Q9_Q: usize = 9;

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 discrete velocity components (x direction).
pub const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 discrete velocity components (y direction).
pub const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 lattice weights.
pub const D2Q9_W: [f64; 9] = [
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

// ---------------------------------------------------------------------------
// MhdState
// ---------------------------------------------------------------------------

/// Complete MHD state at a single spatial point.
///
/// Stores density ρ, velocity **u**, magnetic field **B**, and pressure p.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MhdState {
    /// Mass density ρ (kg m⁻³).
    pub rho: f64,
    /// Velocity vector \[ux, uy, uz\] (m s⁻¹).
    pub velocity: [f64; 3],
    /// Magnetic field vector \[Bx, By, Bz\] (T).
    pub b_field: [f64; 3],
    /// Thermal pressure p (Pa).
    pub pressure: f64,
}

impl MhdState {
    /// Create a new MHD state with specified values.
    ///
    /// # Arguments
    /// * `rho`      – mass density (kg m⁻³)
    /// * `velocity` – velocity vector \[ux, uy, uz\] (m s⁻¹)
    /// * `b_field`  – magnetic field \[Bx, By, Bz\] (T)
    /// * `pressure` – thermal pressure (Pa)
    pub fn new(rho: f64, velocity: [f64; 3], b_field: [f64; 3], pressure: f64) -> Self {
        Self {
            rho,
            velocity,
            b_field,
            pressure,
        }
    }

    /// Compute total energy density: kinetic + magnetic + thermal.
    ///
    /// E_total = ½ρu² + B²/(2μ₀) + p/(γ-1)
    /// For convenience γ is taken as 5/3 (monatomic ideal gas).
    pub fn total_energy_density(&self) -> f64 {
        let ke = 0.5 * self.rho * self.velocity.iter().map(|&v| v * v).sum::<f64>();
        let b2 = self.b_field.iter().map(|&b| b * b).sum::<f64>();
        let me = b2 / (2.0 * MU_0);
        let te = self.pressure / (5.0 / 3.0 - 1.0);
        ke + me + te
    }

    /// Compute magnetic pressure p_B = B²/(2μ₀).
    pub fn magnetic_pressure(&self) -> f64 {
        let b2 = self.b_field.iter().map(|&b| b * b).sum::<f64>();
        b2 / (2.0 * MU_0)
    }

    /// Compute total pressure: thermal + magnetic.
    pub fn total_pressure(&self) -> f64 {
        self.pressure + self.magnetic_pressure()
    }
}

impl Default for MhdState {
    fn default() -> Self {
        Self {
            rho: 1.0,
            velocity: [0.0; 3],
            b_field: [0.0; 3],
            pressure: 1.0,
        }
    }
}

// ---------------------------------------------------------------------------
// LorentzForce
// ---------------------------------------------------------------------------

/// Computes the Lorentz body force **F** = J×B on a conducting fluid.
///
/// The current density is estimated from a finite-difference curl of B:
/// J = ∇×B / μ₀
///
/// # References
/// - Priest, E. & Forbes, T. (2000). *Magnetic Reconnection*. Cambridge UP.
#[derive(Debug, Clone)]
pub struct LorentzForce {
    /// Magnetic permeability (H m⁻¹), defaults to μ₀.
    pub mu: f64,
    /// Grid spacing used for finite-difference curl estimate (m).
    pub dx: f64,
}

impl LorentzForce {
    /// Construct a new `LorentzForce` helper.
    ///
    /// # Arguments
    /// * `mu` – permeability (H m⁻¹)
    /// * `dx` – spatial resolution (m)
    pub fn new(mu: f64, dx: f64) -> Self {
        Self { mu, dx }
    }

    /// Compute J×B given the local magnetic field **B** and its neighbours.
    ///
    /// Uses a central-difference 2D curl (in x-y plane):
    /// Jz = (∂By/∂x − ∂Bx/∂y) / μ
    ///
    /// Returns the Lorentz force vector \[Fx, Fy, Fz\] (N m⁻³).
    ///
    /// # Arguments
    /// * `b_center`   – B at (i,j)
    /// * `b_xp`       – B at (i+1,j)
    /// * `b_xm`       – B at (i-1,j)
    /// * `b_yp`       – B at (i,j+1)
    /// * `b_ym`       – B at (i,j-1)
    #[allow(clippy::too_many_arguments)]
    pub fn compute(
        &self,
        b_center: [f64; 3],
        b_xp: [f64; 3],
        b_xm: [f64; 3],
        b_yp: [f64; 3],
        b_ym: [f64; 3],
    ) -> [f64; 3] {
        let two_dx = 2.0 * self.dx;
        // Central differences for curl in 2-D (Jz component)
        let jz = (b_yp[0] - b_ym[0]) / two_dx   // ∂Bx/∂y  (note: curl sign)
            - (b_xp[1] - b_xm[1]) / two_dx; // ∂By/∂x
        let jz = jz / self.mu;

        // Jx and Jy from out-of-plane gradients (simplified, set zero in 2D)
        let jx = 0.0_f64;
        let jy = 0.0_f64;

        let j = [jx, jy, jz];
        let b = b_center;
        // J × B
        [
            j[1] * b[2] - j[2] * b[1],
            j[2] * b[0] - j[0] * b[2],
            j[0] * b[1] - j[1] * b[0],
        ]
    }

    /// Compute the current density J = ∇×B/μ from finite differences.
    ///
    /// Returns J vector (A m⁻²).
    #[allow(clippy::too_many_arguments)]
    pub fn current_density(
        &self,
        b_xp: [f64; 3],
        b_xm: [f64; 3],
        b_yp: [f64; 3],
        b_ym: [f64; 3],
    ) -> [f64; 3] {
        let two_dx = 2.0 * self.dx;
        let jz = ((b_xp[1] - b_xm[1]) - (b_yp[0] - b_ym[0])) / (two_dx * self.mu);
        [0.0, 0.0, jz]
    }
}

// ---------------------------------------------------------------------------
// InductionEquation
// ---------------------------------------------------------------------------

/// Solves the magnetic induction equation on a 2D grid.
///
/// ∂B/∂t = ∇×(u×B) + η∇²B
///
/// Implements a simple explicit finite-difference scheme (Euler in time,
/// central differences in space).
#[derive(Debug, Clone)]
pub struct InductionEquation {
    /// Grid width (number of cells).
    pub nx: usize,
    /// Grid height (number of cells).
    pub ny: usize,
    /// Cell size (m).
    pub dx: f64,
    /// Magnetic diffusivity η = 1/(μ₀σ) (m² s⁻¹).
    pub eta: f64,
    /// Magnetic field x-component on the grid (T).
    pub bx: Vec<f64>,
    /// Magnetic field y-component on the grid (T).
    pub by: Vec<f64>,
    /// Magnetic field z-component on the grid (T).
    pub bz: Vec<f64>,
}

impl InductionEquation {
    /// Create a new induction equation solver with zero initial B.
    ///
    /// # Arguments
    /// * `nx`  – grid width
    /// * `ny`  – grid height
    /// * `dx`  – cell size (m)
    /// * `eta` – magnetic diffusivity (m² s⁻¹)
    pub fn new(nx: usize, ny: usize, dx: f64, eta: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            dx,
            eta,
            bx: vec![0.0; n],
            by: vec![0.0; n],
            bz: vec![0.0; n],
        }
    }

    /// Set the magnetic field at cell (i, j).
    pub fn set_b(&mut self, i: usize, j: usize, bx: f64, by: f64, bz: f64) {
        let idx = j * self.nx + i;
        self.bx[idx] = bx;
        self.by[idx] = by;
        self.bz[idx] = bz;
    }

    /// Get the magnetic field at cell (i, j).
    pub fn get_b(&self, i: usize, j: usize) -> [f64; 3] {
        let idx = j * self.nx + i;
        [self.bx[idx], self.by[idx], self.bz[idx]]
    }

    /// Compute diffusion of B (Laplacian term η∇²B).
    ///
    /// Returns dB/dt from the diffusive term only.
    pub fn diffusion_rhs(&self, i: usize, j: usize) -> [f64; 3] {
        let nx = self.nx;
        let ip = if i + 1 < nx { i + 1 } else { 0 };
        let im = if i > 0 { i - 1 } else { nx - 1 };
        let jp = if j + 1 < self.ny { j + 1 } else { 0 };
        let jm = if j > 0 { j - 1 } else { self.ny - 1 };

        let center = self.get_b(i, j);
        let xp = self.get_b(ip, j);
        let xm = self.get_b(im, j);
        let yp = self.get_b(i, jp);
        let ym = self.get_b(i, jm);

        let dx2 = self.dx * self.dx;
        let mut rhs = [0.0f64; 3];
        for k in 0..3 {
            let arr = [center[k], xp[k], xm[k], yp[k], ym[k]];
            rhs[k] = self.eta * (arr[1] + arr[2] + arr[3] + arr[4] - 4.0 * arr[0]) / dx2;
        }
        rhs
    }

    /// Advance the induction equation by one time step (pure diffusion limit).
    ///
    /// Uses explicit Euler time integration. For production use a more
    /// stable implicit or Runge-Kutta scheme.
    ///
    /// # Arguments
    /// * `dt` – time step (s)
    pub fn step_diffusion(&mut self, dt: f64) {
        let mut new_bx = self.bx.clone();
        let mut new_by = self.by.clone();
        let mut new_bz = self.bz.clone();
        for j in 0..self.ny {
            for i in 0..self.nx {
                let rhs = self.diffusion_rhs(i, j);
                let idx = j * self.nx + i;
                new_bx[idx] += dt * rhs[0];
                new_by[idx] += dt * rhs[1];
                new_bz[idx] += dt * rhs[2];
            }
        }
        self.bx = new_bx;
        self.by = new_by;
        self.bz = new_bz;
    }

    /// Total magnetic energy ∫B²/(2μ₀) dV (J m⁻¹ for 2D).
    pub fn magnetic_energy(&self) -> f64 {
        let mut e = 0.0;
        for k in 0..self.bx.len() {
            let b2 = self.bx[k] * self.bx[k] + self.by[k] * self.by[k] + self.bz[k] * self.bz[k];
            e += b2;
        }
        e * self.dx * self.dx / (2.0 * MU_0)
    }
}

// ---------------------------------------------------------------------------
// AlfvenWave
// ---------------------------------------------------------------------------

/// Alfvén wave propagation properties in a magnetised plasma.
///
/// Alfvén waves are transverse MHD waves propagating along magnetic field
/// lines with phase velocity v_A = B/√(μ₀ρ).
#[derive(Debug, Clone, Copy)]
pub struct AlfvenWave {
    /// Magnetic field magnitude |B| (T).
    pub b_mag: f64,
    /// Mass density ρ (kg m⁻³).
    pub rho: f64,
    /// Magnetic permeability (H m⁻¹).
    pub mu: f64,
    /// Wave number k (rad m⁻¹).
    pub k: f64,
}

impl AlfvenWave {
    /// Construct a new `AlfvenWave` instance.
    ///
    /// # Arguments
    /// * `b_mag` – magnetic field strength (T)
    /// * `rho`   – mass density (kg m⁻³)
    /// * `mu`    – permeability (H m⁻¹); pass `MU_0` for vacuum
    /// * `k`     – wave number (rad m⁻¹)
    pub fn new(b_mag: f64, rho: f64, mu: f64, k: f64) -> Self {
        Self { b_mag, rho, mu, k }
    }

    /// Alfvén speed v_A = B / √(μ ρ) (m s⁻¹).
    ///
    /// ```no_run
    /// use oxiphysics_lbm::magnetohydrodynamics_lbm::{AlfvenWave, MU_0};
    /// let wave = AlfvenWave::new(1e-3, 1e-6, MU_0, 1.0);
    /// let va = wave.alfven_speed();
    /// assert!(va > 0.0);
    /// ```
    pub fn alfven_speed(&self) -> f64 {
        if self.rho <= 0.0 || self.mu <= 0.0 {
            return 0.0;
        }
        self.b_mag / (self.mu * self.rho).sqrt()
    }

    /// Angular frequency ω = v_A · k (rad s⁻¹).
    pub fn angular_frequency(&self) -> f64 {
        self.alfven_speed() * self.k
    }

    /// Wave period T = 2π/ω (s).
    pub fn period(&self) -> f64 {
        let omega = self.angular_frequency();
        if omega <= 0.0 {
            return f64::INFINITY;
        }
        2.0 * PI / omega
    }

    /// Phase velocity (same as Alfvén speed for pure Alfvén waves).
    pub fn phase_velocity(&self) -> f64 {
        self.alfven_speed()
    }

    /// Group velocity (same as phase velocity for non-dispersive Alfvén waves).
    pub fn group_velocity(&self) -> f64 {
        self.alfven_speed()
    }

    /// Magnetic tension force per unit volume τ = B²/(μρ) (N m⁻³).
    pub fn magnetic_tension(&self) -> f64 {
        if self.mu <= 0.0 || self.rho <= 0.0 {
            return 0.0;
        }
        self.b_mag * self.b_mag / (self.mu * self.rho)
    }

    /// Fast magnetosonic speed v_f = √(v_A² + c_s²) in the parallel direction.
    ///
    /// # Arguments
    /// * `c_s` – sound speed (m s⁻¹)
    pub fn fast_magnetosonic_speed(&self, c_s: f64) -> f64 {
        let va = self.alfven_speed();
        (va * va + c_s * c_s).sqrt()
    }

    /// Slow magnetosonic speed (perpendicular component).
    ///
    /// # Arguments
    /// * `c_s` – sound speed (m s⁻¹)
    pub fn slow_magnetosonic_speed(&self, c_s: f64) -> f64 {
        let va = self.alfven_speed();
        // v_slow = |v_A - c_s| / √2 (simplified)
        ((va - c_s).abs()) / 2.0_f64.sqrt()
    }
}

// ---------------------------------------------------------------------------
// MagneticReconnection
// ---------------------------------------------------------------------------

/// Sweet-Parker magnetic reconnection diagnostics.
///
/// In the Sweet-Parker model, antiparallel magnetic fields annihilate in a
/// thin current sheet of thickness δ = L/√Rm, releasing stored magnetic
/// energy as kinetic energy and heat.
///
/// # References
/// - Parker, E.N. (1957). Sweet's mechanism for merging magnetic fields.
///   *J. Geophys. Res.*, 62, 509.
#[derive(Debug, Clone, Copy)]
pub struct MagneticReconnection {
    /// Upstream magnetic field B₀ (T).
    pub b0: f64,
    /// Upstream density ρ₀ (kg m⁻³).
    pub rho0: f64,
    /// Current sheet half-length L (m).
    pub length: f64,
    /// Magnetic diffusivity η (m² s⁻¹).
    pub eta: f64,
    /// Magnetic permeability (H m⁻¹).
    pub mu: f64,
}

impl MagneticReconnection {
    /// Construct a new reconnection model.
    ///
    /// # Arguments
    /// * `b0`     – upstream B field (T)
    /// * `rho0`   – upstream density (kg m⁻³)
    /// * `length` – current sheet half-length (m)
    /// * `eta`    – magnetic diffusivity (m² s⁻¹)
    /// * `mu`     – permeability (H m⁻¹)
    pub fn new(b0: f64, rho0: f64, length: f64, eta: f64, mu: f64) -> Self {
        Self {
            b0,
            rho0,
            length,
            eta,
            mu,
        }
    }

    /// Alfvén speed v_A = B₀/√(μρ₀) (m s⁻¹).
    pub fn alfven_speed(&self) -> f64 {
        if self.rho0 <= 0.0 || self.mu <= 0.0 {
            return 0.0;
        }
        self.b0 / (self.mu * self.rho0).sqrt()
    }

    /// Magnetic Reynolds number Rm = v_A L / η.
    pub fn magnetic_reynolds(&self) -> f64 {
        if self.eta <= 0.0 {
            return f64::INFINITY;
        }
        self.alfven_speed() * self.length / self.eta
    }

    /// Sweet-Parker reconnection rate E_SP = v_A / √Rm (m s⁻¹).
    ///
    /// The inflow speed (reconnection rate) is:
    /// v_in = v_A / √Rm
    pub fn reconnection_rate(&self) -> f64 {
        let rm = self.magnetic_reynolds();
        if rm <= 0.0 || rm.is_infinite() {
            return 0.0;
        }
        self.alfven_speed() / rm.sqrt()
    }

    /// Sweet-Parker current sheet thickness δ = L / √Rm (m).
    pub fn sheet_thickness(&self) -> f64 {
        let rm = self.magnetic_reynolds();
        if rm <= 0.0 || rm.is_infinite() {
            return 0.0;
        }
        self.length / rm.sqrt()
    }

    /// Magnetic energy density upstream B₀²/(2μ) (J m⁻³).
    pub fn upstream_magnetic_energy_density(&self) -> f64 {
        self.b0 * self.b0 / (2.0 * self.mu)
    }

    /// Energy release rate per unit area (W m⁻²).
    ///
    /// Power = v_in × E_mag_density × area, per unit area:
    /// P = v_in × B₀²/(2μ)
    pub fn energy_release_rate(&self) -> f64 {
        self.reconnection_rate() * self.upstream_magnetic_energy_density()
    }

    /// Outflow velocity (approximately Alfvén speed) (m s⁻¹).
    pub fn outflow_velocity(&self) -> f64 {
        self.alfven_speed()
    }
}

// ---------------------------------------------------------------------------
// MhdEquilibrium
// ---------------------------------------------------------------------------

/// Force-free and Grad-Shafranov MHD equilibrium solver.
///
/// A force-free equilibrium satisfies J×B = 0, i.e., ∇×B = αB.
/// The Grad-Shafranov equation is a 2D MHD equilibrium for axisymmetric
/// configurations (tokamaks, flux tubes):
///
/// R ∂/∂R (1/R ∂ψ/∂R) + ∂²ψ/∂Z² = −μ₀ R² dp/dψ − F dF/dψ
///
/// where ψ is the poloidal flux function.
#[derive(Debug, Clone)]
pub struct MhdEquilibrium {
    /// Grid size in R direction.
    pub nr: usize,
    /// Grid size in Z direction.
    pub nz: usize,
    /// Radial grid spacing ΔR (m).
    pub dr: f64,
    /// Vertical grid spacing ΔZ (m).
    pub dz: f64,
    /// Poloidal flux function ψ(R,Z) (Wb).
    pub psi: Vec<f64>,
    /// Plasma pressure profile p(ψ).
    pub pressure_coeff: f64,
    /// Toroidal field function F(ψ).
    pub f_coeff: f64,
}

impl MhdEquilibrium {
    /// Construct a new `MhdEquilibrium` grid.
    ///
    /// # Arguments
    /// * `nr` – radial resolution
    /// * `nz` – vertical resolution
    /// * `dr` – radial cell size (m)
    /// * `dz` – vertical cell size (m)
    pub fn new(nr: usize, nz: usize, dr: f64, dz: f64) -> Self {
        Self {
            nr,
            nz,
            dr,
            dz,
            psi: vec![0.0; nr * nz],
            pressure_coeff: 1.0,
            f_coeff: 1.0,
        }
    }

    /// Set ψ at grid point (ir, iz).
    pub fn set_psi(&mut self, ir: usize, iz: usize, val: f64) {
        self.psi[iz * self.nr + ir] = val;
    }

    /// Get ψ at grid point (ir, iz).
    pub fn get_psi(&self, ir: usize, iz: usize) -> f64 {
        self.psi[iz * self.nr + ir]
    }

    /// Apply a simple circular flux-surface initialisation:
    /// ψ(R,Z) = ψ₀ × exp(−(R² + Z²) / a²)
    ///
    /// # Arguments
    /// * `psi0` – peak flux (Wb)
    /// * `a`    – width parameter (m)
    /// * `r0`   – major radius offset (m)
    pub fn init_circular(&mut self, psi0: f64, a: f64, r0: f64) {
        for iz in 0..self.nz {
            for ir in 0..self.nr {
                let r = (ir as f64 + 0.5) * self.dr - r0;
                let z = (iz as f64 + 0.5) * self.dz - (self.nz as f64 * self.dz) / 2.0;
                self.psi[iz * self.nr + ir] = psi0 * (-(r * r + z * z) / (a * a)).exp();
            }
        }
    }

    /// Evaluate Grad-Shafranov RHS at a single grid point.
    ///
    /// Returns the residual of the GS equation.
    pub fn gs_residual(&self, ir: usize, iz: usize) -> f64 {
        if ir == 0 || ir + 1 >= self.nr || iz == 0 || iz + 1 >= self.nz {
            return 0.0;
        }
        let r = (ir as f64 + 0.5) * self.dr;
        let psi_c = self.get_psi(ir, iz);
        let psi_rp = self.get_psi(ir + 1, iz);
        let psi_rm = self.get_psi(ir - 1, iz);
        let psi_zp = self.get_psi(ir, iz + 1);
        let psi_zm = self.get_psi(ir, iz - 1);

        // Elliptic operator: R ∂/∂R(1/R ∂ψ/∂R) + ∂²ψ/∂Z²
        let d2psi_rr = (psi_rp - 2.0 * psi_c + psi_rm) / (self.dr * self.dr);
        let dpsi_r = (psi_rp - psi_rm) / (2.0 * self.dr);
        let d2psi_zz = (psi_zp - 2.0 * psi_c + psi_zm) / (self.dz * self.dz);

        let lhs = d2psi_rr - dpsi_r / r + d2psi_zz;
        let rhs = -MU_0 * r * r * self.pressure_coeff * psi_c - self.f_coeff * self.f_coeff * psi_c;
        lhs - rhs
    }

    /// Compute the force-free parameter α = (∇×B · B) / B² at a grid point.
    ///
    /// For a truly force-free field, α is constant along field lines.
    pub fn force_free_alpha(&self, ir: usize, iz: usize) -> f64 {
        let psi = self.get_psi(ir, iz);
        if psi.abs() < 1e-30 {
            return 0.0;
        }
        // Simplified: α ≈ F(ψ) / R²
        let r = (ir as f64 + 0.5) * self.dr;
        self.f_coeff / (r * r) * psi
    }
}

// ---------------------------------------------------------------------------
// MhdAnalysis
// ---------------------------------------------------------------------------

/// Dimensionless MHD numbers and plasma diagnostics.
///
/// Provides methods to compute the key dimensionless parameters that
/// characterise MHD flows:
/// - Magnetic Reynolds number Rm = UL/η
/// - Plasma beta β = p/(B²/2μ)
/// - Alfvén Mach number M_A = U/v_A
/// - Hartmann number Ha = BL√(σ/ρν)
#[derive(Debug, Clone, Copy)]
pub struct MhdAnalysis {
    /// Characteristic velocity U (m s⁻¹).
    pub velocity: f64,
    /// Characteristic length L (m).
    pub length: f64,
    /// Magnetic diffusivity η (m² s⁻¹).
    pub eta: f64,
    /// Kinematic viscosity ν (m² s⁻¹).
    pub nu: f64,
    /// Magnetic field magnitude B (T).
    pub b_mag: f64,
    /// Mass density ρ (kg m⁻³).
    pub rho: f64,
    /// Thermal pressure p (Pa).
    pub pressure: f64,
    /// Magnetic permeability μ (H m⁻¹).
    pub mu: f64,
}

impl MhdAnalysis {
    /// Construct a new `MhdAnalysis` instance.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        velocity: f64,
        length: f64,
        eta: f64,
        nu: f64,
        b_mag: f64,
        rho: f64,
        pressure: f64,
        mu: f64,
    ) -> Self {
        Self {
            velocity,
            length,
            eta,
            nu,
            b_mag,
            rho,
            pressure,
            mu,
        }
    }

    /// Magnetic Reynolds number Rm = U L / η.
    ///
    /// Rm ≫ 1: advection dominates (ideal MHD limit).
    /// Rm ≪ 1: diffusion dominates (resistive limit).
    pub fn magnetic_reynolds(&self) -> f64 {
        if self.eta <= 0.0 {
            return f64::INFINITY;
        }
        self.velocity * self.length / self.eta
    }

    /// Plasma beta β = p / (B²/2μ).
    ///
    /// β ≪ 1: magnetically dominated plasma.
    /// β ≫ 1: thermally dominated plasma.
    ///
    /// ```no_run
    /// use oxiphysics_lbm::magnetohydrodynamics_lbm::{MhdAnalysis, MU_0};
    /// let m = MhdAnalysis::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, MU_0);
    /// let beta = m.plasma_beta();
    /// assert!((beta - 2.0 * MU_0).abs() < 1e-20);
    /// ```
    pub fn plasma_beta(&self) -> f64 {
        let b2_2mu = self.b_mag * self.b_mag / (2.0 * self.mu);
        if b2_2mu <= 0.0 {
            return f64::INFINITY;
        }
        self.pressure / b2_2mu
    }

    /// Alfvén Mach number M_A = U / v_A.
    ///
    /// M_A < 1: sub-Alfvénic flow.
    /// M_A > 1: super-Alfvénic flow.
    pub fn alfven_mach(&self) -> f64 {
        let va = self.alfven_speed();
        if va <= 0.0 {
            return f64::INFINITY;
        }
        self.velocity / va
    }

    /// Alfvén speed v_A = B / √(μρ) (m s⁻¹).
    pub fn alfven_speed(&self) -> f64 {
        if self.rho <= 0.0 || self.mu <= 0.0 {
            return 0.0;
        }
        self.b_mag / (self.mu * self.rho).sqrt()
    }

    /// Hartmann number Ha = B L √(1/(ρ ν η)).
    ///
    /// Ha characterises the ratio of magnetic to viscous forces.
    pub fn hartmann_number(&self) -> f64 {
        if self.rho <= 0.0 || self.nu <= 0.0 || self.eta <= 0.0 {
            return 0.0;
        }
        self.b_mag * self.length / (self.rho * self.nu * self.eta).sqrt()
    }

    /// Lundquist number S = v_A L / η (magnetic Rm based on Alfvén speed).
    pub fn lundquist_number(&self) -> f64 {
        if self.eta <= 0.0 {
            return f64::INFINITY;
        }
        self.alfven_speed() * self.length / self.eta
    }

    /// Magnetic Prandtl number Pm = ν / η.
    pub fn magnetic_prandtl(&self) -> f64 {
        if self.eta <= 0.0 {
            return f64::INFINITY;
        }
        self.nu / self.eta
    }

    /// Hydrodynamic Reynolds number Re = U L / ν.
    pub fn reynolds_number(&self) -> f64 {
        if self.nu <= 0.0 {
            return f64::INFINITY;
        }
        self.velocity * self.length / self.nu
    }
}

// ---------------------------------------------------------------------------
// MhdLbm
// ---------------------------------------------------------------------------

/// Dual-distribution MHD-LBM solver on a D2Q9 lattice.
///
/// Uses two sets of distribution functions:
/// - `f[i][q]`: fluid distributions → density ρ, momentum ρu
/// - `g[i][q]`: magnetic distributions → magnetic field B
///
/// The collision step uses BGK relaxation for both sets:
/// f_i* = f_i − (f_i − f_i^eq) / τ_f
/// g_i* = g_i − (g_i − g_i^eq) / τ_g
///
/// # References
/// - Dellar, P.J. (2002). *J. Comput. Phys.*, 179, 95–126.
#[derive(Debug, Clone)]
pub struct MhdLbm {
    /// Grid width (cells).
    pub nx: usize,
    /// Grid height (cells).
    pub ny: usize,
    /// Fluid relaxation time τ_f (lattice units, τ > 0.5).
    pub tau_f: f64,
    /// Magnetic relaxation time τ_g (lattice units).
    pub tau_g: f64,
    /// Fluid distribution functions f\[cell\]\[q\].
    pub f: Vec<[f64; D2Q9_Q]>,
    /// Magnetic distribution functions g\[cell\]\[q\].
    pub g: Vec<[f64; D2Q9_Q]>,
    /// Density field ρ (lattice units).
    pub rho: Vec<f64>,
    /// Velocity field ux (lattice units).
    pub ux: Vec<f64>,
    /// Velocity field uy (lattice units).
    pub uy: Vec<f64>,
    /// Magnetic field Bx (lattice units).
    pub bx: Vec<f64>,
    /// Magnetic field By (lattice units).
    pub by: Vec<f64>,
    /// Simulation step counter.
    pub step: u64,
}

impl MhdLbm {
    /// Create a new MHD-LBM grid initialised with uniform state.
    ///
    /// # Arguments
    /// * `nx`    – grid width
    /// * `ny`    – grid height
    /// * `tau_f` – fluid BGK relaxation time (must be > 0.5)
    /// * `tau_g` – magnetic BGK relaxation time (must be > 0.5)
    pub fn new(nx: usize, ny: usize, tau_f: f64, tau_g: f64) -> Self {
        let n = nx * ny;
        let rho = vec![1.0; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let bx = vec![0.1; n];
        let by = vec![0.0; n];

        let mut f = vec![[0.0f64; D2Q9_Q]; n];
        let mut g = vec![[0.0f64; D2Q9_Q]; n];
        for i in 0..n {
            f[i] = fluid_equilibrium(rho[i], ux[i], uy[i]);
            g[i] = magnetic_equilibrium(bx[i], by[i], ux[i], uy[i]);
        }
        Self {
            nx,
            ny,
            tau_f,
            tau_g,
            f,
            g,
            rho,
            ux,
            uy,
            bx,
            by,
            step: 0,
        }
    }

    /// Perform one LBM time step: collision + streaming + macroscopic update.
    pub fn advance(&mut self) {
        self.collide();
        self.stream_f();
        self.stream_g();
        self.update_macroscopic();
        self.step += 1;
    }

    /// BGK collision for both f and g distribution functions.
    fn collide(&mut self) {
        let inv_tauf = 1.0 / self.tau_f;
        let inv_taug = 1.0 / self.tau_g;
        for i in 0..self.nx * self.ny {
            let feq = fluid_equilibrium(self.rho[i], self.ux[i], self.uy[i]);
            let geq = magnetic_equilibrium(self.bx[i], self.by[i], self.ux[i], self.uy[i]);
            for q in 0..D2Q9_Q {
                self.f[i][q] += (feq[q] - self.f[i][q]) * inv_tauf;
                self.g[i][q] += (geq[q] - self.g[i][q]) * inv_taug;
            }
        }
    }

    /// Stream fluid distributions (periodic boundary).
    fn stream_f(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![[0.0f64; D2Q9_Q]; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let src = j * nx + i;
                for q in 0..D2Q9_Q {
                    let ni = (i as isize + D2Q9_EX[q] as isize).rem_euclid(nx as isize) as usize;
                    let nj = (j as isize + D2Q9_EY[q] as isize).rem_euclid(ny as isize) as usize;
                    let dst = nj * nx + ni;
                    f_new[dst][q] = self.f[src][q];
                }
            }
        }
        self.f = f_new;
    }

    /// Stream magnetic distributions (periodic boundary).
    fn stream_g(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut g_new = vec![[0.0f64; D2Q9_Q]; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let src = j * nx + i;
                for q in 0..D2Q9_Q {
                    let ni = (i as isize + D2Q9_EX[q] as isize).rem_euclid(nx as isize) as usize;
                    let nj = (j as isize + D2Q9_EY[q] as isize).rem_euclid(ny as isize) as usize;
                    let dst = nj * nx + ni;
                    g_new[dst][q] = self.g[src][q];
                }
            }
        }
        self.g = g_new;
    }

    /// Update macroscopic variables ρ, u, B from the distributions.
    fn update_macroscopic(&mut self) {
        for i in 0..self.nx * self.ny {
            let mut rho = 0.0;
            let mut jx = 0.0;
            let mut jy = 0.0;
            let mut bx = 0.0;
            let mut by = 0.0;
            for q in 0..D2Q9_Q {
                rho += self.f[i][q];
                jx += self.f[i][q] * D2Q9_EX[q];
                jy += self.f[i][q] * D2Q9_EY[q];
                bx += self.g[i][q] * D2Q9_EX[q];
                by += self.g[i][q] * D2Q9_EY[q];
            }
            self.rho[i] = rho.max(1e-12);
            self.ux[i] = jx / self.rho[i];
            self.uy[i] = jy / self.rho[i];
            self.bx[i] = bx;
            self.by[i] = by;
        }
    }

    /// Get the MHD state at grid cell (i, j).
    pub fn state_at(&self, i: usize, j: usize) -> MhdState {
        let idx = j * self.nx + i;
        let cs2 = 1.0 / 3.0;
        MhdState {
            rho: self.rho[idx],
            velocity: [self.ux[idx], self.uy[idx], 0.0],
            b_field: [self.bx[idx], self.by[idx], 0.0],
            pressure: self.rho[idx] * cs2,
        }
    }

    /// Total fluid kinetic energy (lattice units).
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.nx * self.ny {
            ke += 0.5 * self.rho[i] * (self.ux[i] * self.ux[i] + self.uy[i] * self.uy[i]);
        }
        ke
    }

    /// Total magnetic energy (lattice units).
    pub fn magnetic_energy(&self) -> f64 {
        let mut me = 0.0;
        for i in 0..self.nx * self.ny {
            me += 0.5 * (self.bx[i] * self.bx[i] + self.by[i] * self.by[i]);
        }
        me
    }

    /// Inject velocity perturbation at cell (i, j).
    ///
    /// Useful for Alfvén wave tests.
    pub fn perturb_velocity(&mut self, i: usize, j: usize, dux: f64, duy: f64) {
        let idx = j * self.nx + i;
        let new_ux = self.ux[idx] + dux;
        let new_uy = self.uy[idx] + duy;
        self.f[idx] = fluid_equilibrium(self.rho[idx], new_ux, new_uy);
        self.ux[idx] = new_ux;
        self.uy[idx] = new_uy;
    }
}

// ---------------------------------------------------------------------------
// Free functions — equilibria
// ---------------------------------------------------------------------------

/// Fluid D2Q9 equilibrium distribution function.
///
/// f_i^eq = w_i ρ \[1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) − u²/(2cs²)\]
///
/// Returns the 9 equilibrium populations.
pub fn fluid_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; D2Q9_Q] {
    let cs2 = 1.0 / 3.0;
    let cs4 = cs2 * cs2;
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; D2Q9_Q];
    for q in 0..D2Q9_Q {
        let eu = D2Q9_EX[q] * ux + D2Q9_EY[q] * uy;
        feq[q] = D2Q9_W[q] * rho * (1.0 + eu / cs2 + eu * eu / (2.0 * cs4) - u2 / (2.0 * cs2));
    }
    feq
}

/// Magnetic D2Q9 equilibrium distribution function (Dellar 2002).
///
/// The magnetic equilibrium is constructed to reproduce the induction
/// equation in the hydrodynamic limit:
///
/// g_i^eq = w_i \[B + (e_i · B)(e_i · u)/cs² − (B · u) e_i / cs²\]
///
/// Returns the 9 equilibrium populations for Bx and By as a single
/// flat array of length 9 (packing Bx component only; By stored separately).
pub fn magnetic_equilibrium(bx: f64, by: f64, ux: f64, uy: f64) -> [f64; D2Q9_Q] {
    let cs2 = 1.0 / 3.0;
    let mut geq = [0.0f64; D2Q9_Q];
    for q in 0..D2Q9_Q {
        let eu = D2Q9_EX[q] * ux + D2Q9_EY[q] * uy;
        let eb = D2Q9_EX[q] * bx + D2Q9_EY[q] * by;
        let bu = bx * ux + by * uy;
        geq[q] = D2Q9_W[q] * (eb + (eb * eu - bu * D2Q9_EX[q]) / cs2);
    }
    geq
}

/// Compute the D2Q9 sound speed in lattice units.
///
/// In the D2Q9 model, cs = 1/√3 in lattice units.
pub fn lattice_sound_speed() -> f64 {
    (1.0_f64 / 3.0).sqrt()
}

/// Estimate the MHD CFL time step limit.
///
/// Δt ≤ Δx / (|u| + v_A + c_s)
///
/// # Arguments
/// * `dx` – cell size (m)
/// * `u`  – flow speed (m s⁻¹)
/// * `va` – Alfvén speed (m s⁻¹)
/// * `cs` – sound speed (m s⁻¹)
pub fn mhd_cfl_timestep(dx: f64, u: f64, va: f64, cs: f64) -> f64 {
    let denom = u.abs() + va + cs;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    dx / denom
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- MhdState tests ---

    #[test]
    fn test_mhd_state_default() {
        let s = MhdState::default();
        assert_eq!(s.rho, 1.0);
        assert_eq!(s.pressure, 1.0);
        assert_eq!(s.b_field, [0.0; 3]);
    }

    #[test]
    fn test_mhd_state_magnetic_pressure() {
        let b = [1.0, 0.0, 0.0];
        let s = MhdState::new(1.0, [0.0; 3], b, 0.0);
        let p_b = s.magnetic_pressure();
        let expected = 1.0 / (2.0 * MU_0);
        assert!((p_b - expected).abs() < 1e-6);
    }

    #[test]
    fn test_mhd_state_total_pressure() {
        let b = [1.0, 0.0, 0.0];
        let s = MhdState::new(1.0, [0.0; 3], b, 2.0);
        assert!(s.total_pressure() > s.pressure);
    }

    #[test]
    fn test_mhd_state_total_energy_positive() {
        let s = MhdState::new(1.0, [1.0, 0.0, 0.0], [0.1, 0.0, 0.0], 1.0);
        assert!(s.total_energy_density() > 0.0);
    }

    // --- LorentzForce tests ---

    #[test]
    fn test_lorentz_force_zero_uniform_field() {
        let lf = LorentzForce::new(MU_0, 1.0);
        let b = [1.0, 0.0, 0.0];
        let force = lf.compute(b, b, b, b, b);
        // Uniform B → zero current → zero force
        for c in force {
            assert!(c.abs() < EPS);
        }
    }

    #[test]
    fn test_lorentz_force_nonzero_gradient() {
        let lf = LorentzForce::new(MU_0, 0.1);
        let b_center = [1.0, 0.0, 0.0];
        let b_xp = [1.0, 0.2, 0.0];
        let b_xm = [1.0, -0.2, 0.0];
        let b_yp = [0.5, 0.0, 0.0];
        let b_ym = [0.5, 0.0, 0.0];
        let force = lf.compute(b_center, b_xp, b_xm, b_yp, b_ym);
        // With gradient, force should be non-zero
        let mag = force.iter().map(|&f| f * f).sum::<f64>().sqrt();
        assert!(mag > 0.0);
    }

    #[test]
    fn test_current_density_uniform_zero() {
        let lf = LorentzForce::new(MU_0, 1.0);
        let b = [1.0, 1.0, 0.0];
        let j = lf.current_density(b, b, b, b);
        for c in j {
            assert!(c.abs() < EPS);
        }
    }

    // --- AlfvenWave tests ---

    #[test]
    fn test_alfven_speed_formula() {
        // v_A = B / sqrt(mu * rho)
        let b = 1.0;
        let rho = 1.0;
        let mu = 1.0; // normalised
        let wave = AlfvenWave::new(b, rho, mu, 1.0);
        let va = wave.alfven_speed();
        assert!((va - 1.0).abs() < EPS);
    }

    #[test]
    fn test_alfven_speed_proportional_to_b() {
        let mu = 1.0;
        let rho = 1.0;
        let w1 = AlfvenWave::new(1.0, rho, mu, 1.0);
        let w2 = AlfvenWave::new(2.0, rho, mu, 1.0);
        assert!((w2.alfven_speed() / w1.alfven_speed() - 2.0).abs() < EPS);
    }

    #[test]
    fn test_alfven_speed_inversely_sqrt_rho() {
        let mu = 1.0;
        let b = 1.0;
        let w1 = AlfvenWave::new(b, 1.0, mu, 1.0);
        let w2 = AlfvenWave::new(b, 4.0, mu, 1.0);
        assert!((w1.alfven_speed() / w2.alfven_speed() - 2.0).abs() < EPS);
    }

    #[test]
    fn test_alfven_wave_frequency() {
        let wave = AlfvenWave::new(1.0, 1.0, 1.0, 2.0 * PI);
        let omega = wave.angular_frequency();
        assert!((omega - 2.0 * PI).abs() < EPS);
    }

    #[test]
    fn test_alfven_period() {
        let wave = AlfvenWave::new(1.0, 1.0, 1.0, 2.0 * PI);
        let t = wave.period();
        assert!((t - 1.0).abs() < EPS);
    }

    #[test]
    fn test_fast_magnetosonic_exceeds_alfven() {
        let wave = AlfvenWave::new(1.0, 1.0, 1.0, 1.0);
        let vf = wave.fast_magnetosonic_speed(0.5);
        assert!(vf > wave.alfven_speed());
    }

    // --- MagneticReconnection tests ---

    #[test]
    fn test_reconnection_rate_positive() {
        let rec = MagneticReconnection::new(1.0, 1.0, 1.0, 0.01, 1.0);
        assert!(rec.reconnection_rate() > 0.0);
    }

    #[test]
    fn test_reconnection_rate_decreases_with_eta() {
        // Lower η → higher Rm → lower reconnection rate
        let rec_low = MagneticReconnection::new(1.0, 1.0, 1.0, 0.001, 1.0);
        let rec_hi = MagneticReconnection::new(1.0, 1.0, 1.0, 0.1, 1.0);
        assert!(rec_hi.reconnection_rate() > rec_low.reconnection_rate());
    }

    #[test]
    fn test_sheet_thickness_positive() {
        let rec = MagneticReconnection::new(1.0, 1.0, 10.0, 0.01, 1.0);
        assert!(rec.sheet_thickness() > 0.0);
    }

    #[test]
    fn test_energy_release_rate_positive() {
        let rec = MagneticReconnection::new(1.0, 1.0, 1.0, 0.01, 1.0);
        assert!(rec.energy_release_rate() > 0.0);
    }

    // --- MhdAnalysis tests ---

    #[test]
    fn test_plasma_beta() {
        // β = p / (B²/2μ)
        // With p = 1, B = 1, μ = 1: β = 2
        let m = MhdAnalysis::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5);
        let beta = m.plasma_beta();
        assert!((beta - 1.0).abs() < EPS);
    }

    #[test]
    fn test_plasma_beta_mu0() {
        let m = MhdAnalysis::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, MU_0);
        let beta = m.plasma_beta();
        let expected = 1.0 / (1.0 / (2.0 * MU_0));
        assert!((beta - expected).abs() < 1e-20);
    }

    #[test]
    fn test_magnetic_reynolds() {
        let m = MhdAnalysis::new(2.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!((m.magnetic_reynolds() - 6.0).abs() < EPS);
    }

    #[test]
    fn test_alfven_mach_sub() {
        // v_A = B/sqrt(mu*rho) = 1/1 = 1.0; U = 0.5 → sub-Alfvénic
        let m = MhdAnalysis::new(0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(m.alfven_mach() < 1.0);
    }

    #[test]
    fn test_alfven_mach_super() {
        let m = MhdAnalysis::new(2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(m.alfven_mach() > 1.0);
    }

    #[test]
    fn test_hartmann_number_positive() {
        let m = MhdAnalysis::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(m.hartmann_number() > 0.0);
    }

    #[test]
    fn test_magnetic_prandtl() {
        let m = MhdAnalysis::new(1.0, 1.0, 2.0, 4.0, 1.0, 1.0, 1.0, 1.0);
        assert!((m.magnetic_prandtl() - 2.0).abs() < EPS);
    }

    #[test]
    fn test_reynolds_number() {
        let m = MhdAnalysis::new(3.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!((m.reynolds_number() - 15.0).abs() < EPS);
    }

    // --- InductionEquation tests ---

    #[test]
    fn test_induction_diffusion_reduces_field() {
        let mut ind = InductionEquation::new(5, 5, 1.0, 1.0);
        // Set high B in centre
        ind.set_b(2, 2, 10.0, 0.0, 0.0);
        let e0 = ind.magnetic_energy();
        // Take several diffusion steps
        for _ in 0..10 {
            ind.step_diffusion(0.01);
        }
        let e1 = ind.magnetic_energy();
        // Energy should decrease due to diffusion
        assert!(e1 < e0);
    }

    #[test]
    fn test_induction_set_get_b() {
        let mut ind = InductionEquation::new(4, 4, 1.0, 0.1);
        ind.set_b(1, 2, 3.0, 4.0, 5.0);
        let b = ind.get_b(1, 2);
        assert!((b[0] - 3.0).abs() < EPS);
        assert!((b[1] - 4.0).abs() < EPS);
        assert!((b[2] - 5.0).abs() < EPS);
    }

    #[test]
    fn test_induction_diffusion_uniform_stable() {
        let mut ind = InductionEquation::new(4, 4, 1.0, 0.01);
        // Set uniform B
        for j in 0..4 {
            for i in 0..4 {
                ind.set_b(i, j, 1.0, 0.0, 0.0);
            }
        }
        let e0 = ind.magnetic_energy();
        ind.step_diffusion(0.1);
        let e1 = ind.magnetic_energy();
        // Uniform field: diffusion should not change it
        assert!((e0 - e1).abs() / e0 < 1e-12);
    }

    // --- MhdLbm tests ---

    #[test]
    fn test_mhd_lbm_initialises() {
        let lbm = MhdLbm::new(8, 8, 1.0, 1.0);
        assert_eq!(lbm.nx, 8);
        assert_eq!(lbm.ny, 8);
        assert_eq!(lbm.step, 0);
    }

    #[test]
    fn test_mhd_lbm_density_conservation() {
        let mut lbm = MhdLbm::new(8, 8, 1.0, 1.0);
        let total0: f64 = lbm.rho.iter().sum();
        lbm.advance();
        let total1: f64 = lbm.rho.iter().sum();
        assert!((total0 - total1).abs() < 1e-6);
    }

    #[test]
    fn test_mhd_lbm_step_increments() {
        let mut lbm = MhdLbm::new(4, 4, 1.0, 1.0);
        lbm.advance();
        lbm.advance();
        assert_eq!(lbm.step, 2);
    }

    #[test]
    fn test_mhd_lbm_kinetic_energy_non_negative() {
        let lbm = MhdLbm::new(4, 4, 1.0, 1.0);
        assert!(lbm.kinetic_energy() >= 0.0);
    }

    #[test]
    fn test_mhd_lbm_magnetic_energy_non_negative() {
        let lbm = MhdLbm::new(4, 4, 1.0, 1.0);
        assert!(lbm.magnetic_energy() >= 0.0);
    }

    #[test]
    fn test_fluid_equilibrium_sum_to_rho() {
        let rho = 1.5;
        let feq = fluid_equilibrium(rho, 0.1, -0.05);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-12);
    }

    #[test]
    fn test_mhd_cfl_timestep() {
        let dt = mhd_cfl_timestep(1.0, 0.5, 0.3, 0.2);
        assert!(dt > 0.0);
        assert!(dt <= 1.0);
    }

    // --- MhdEquilibrium tests ---

    #[test]
    fn test_gs_init_circular() {
        let mut eq = MhdEquilibrium::new(10, 10, 0.1, 0.1);
        eq.init_circular(1.0, 0.5, 0.5);
        // Peak should be close to the centre
        let centre = eq.get_psi(5, 5);
        let edge = eq.get_psi(0, 0);
        assert!(centre > edge);
    }

    #[test]
    fn test_mhd_lbm_perturb_velocity() {
        let mut lbm = MhdLbm::new(4, 4, 1.0, 1.0);
        lbm.perturb_velocity(2, 2, 0.1, 0.0);
        assert!(lbm.ux[2 * 4 + 2].abs() > 0.0);
    }

    #[test]
    fn test_state_at_returns_correct_rho() {
        let lbm = MhdLbm::new(4, 4, 1.0, 1.0);
        let s = lbm.state_at(0, 0);
        assert!((s.rho - 1.0).abs() < 0.01);
    }
}
