// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrokinetic Lattice Boltzmann Method: electroosmotic flow in a slit
//! channel.
//!
//! The model couples
//! - the **linearised Poisson–Boltzmann (Debye–Hückel) equation**
//!   `∇²φ = κ² φ` (`κ = 1/λ_D`) for the electric double-layer potential with
//!   `φ = ζ` on both walls, solved by successive over-relaxation (SOR); the
//!   net charge density follows from Poisson's equation, `ρ_e = −ε κ² φ`;
//! - a **D2Q9 BGK** fluid with **Guo forcing** driven by the Coulomb body
//!   force `F = ρ_e (E_ext − ∇φ)`.
//!
//! Geometry: the domain is periodic in x and bounded in y by two no-slip,
//! charged walls located half a lattice spacing below row `0` and above row
//! `ny − 1` (half-way bounce-back).  The zeta potential is imposed at the
//! same wall positions through ghost nodes, so the electrical and
//! hydrodynamic walls coincide.
//!
//! Units: [`ElectrokineticParams`] are physical (SI).  The lattice spacing
//! `dx`, time step `dt` and reference density `ρ₀` convert between the two:
//! `τ = 1/2 + 3 ν dt/dx²` with `ν = μ/ρ₀`, lattice force
//! `F_lu = F dt²/(ρ₀ dx)`, physical velocity `u = u_lu dx/dt`.
//!
//! Key items:
//! - [`ElectrokineticLBM`] / [`init_electrokinetic`] / [`step_electrokinetic`].
//! - [`solve_potential`]: SOR solution of the linearised Poisson–Boltzmann problem.
//! - [`electroosmotic_flow_velocity`]: Helmholtz–Smoluchowski slip velocity.
//! - [`poisson_boltzmann_eq`], [`debye_huckel_potential`], [`electrophoretic_mobility`].
//!
//! References:
//! - Probstein, R. F. (1994). *Physicochemical Hydrodynamics*.
//! - Hunter, R. J. (2001). *Foundations of Colloid Science*.
//! - Guo, Z., Zheng, C., Shi, B. (2002). *Phys. Rev. E* 65, 046308.

use crate::error::{Error, Result};

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 discrete velocity vectors `[cx, cy]`.
const D2Q9_C: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];

/// D2Q9 equilibrium weights.
const D2Q9_W: [f64; 9] = [
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

/// Opposite directions of [`D2Q9_C`].
const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// D2Q9 speed of sound squared (lattice units): cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

/// SOR relaxation factor for the Poisson–Boltzmann solve.
const SOR_OMEGA: f64 = 1.7;

/// Relative convergence tolerance of the Poisson–Boltzmann solve.
const POISSON_TOL: f64 = 1e-11;

/// Maximum number of SOR sweeps per solve.
const POISSON_MAX_SWEEPS: usize = 200_000;

// ---------------------------------------------------------------------------
// ElectrokineticParams
// ---------------------------------------------------------------------------

/// Physical and driving parameters for the electrokinetic LBM simulation.
#[derive(Debug, Clone)]
pub struct ElectrokineticParams {
    /// Fluid permittivity ε (F m⁻¹).
    pub permittivity: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Zeta potential ζ (V) at the solid–fluid interface.
    pub zeta_potential: f64,
    /// Debye screening length λ_D (m).
    pub debye_length: f64,
    /// External electric field `[Ex, Ey]` (V m⁻¹).
    pub external_field: [f64; 2],
    /// Lattice spacing Δx (m).
    pub dx: f64,
    /// Time step Δt (s).
    pub dt: f64,
    /// Reference fluid density ρ₀ (kg m⁻³).
    pub rho0: f64,
}

impl ElectrokineticParams {
    /// Default reference density (water, kg m⁻³) used by [`Self::new`].
    pub const DEFAULT_RHO0: f64 = 1000.0;

    /// Default number of lattice spacings per Debye length used by [`Self::new`].
    pub const DEFAULT_NODES_PER_DEBYE_LENGTH: f64 = 4.0;

    /// Construct parameters with an automatic lattice scaling.
    ///
    /// The scaling resolves the double layer with
    /// [`Self::DEFAULT_NODES_PER_DEBYE_LENGTH`] nodes per Debye length
    /// (`dx = λ_D / 4`), uses `ρ₀ =` [`Self::DEFAULT_RHO0`], and chooses `dt` so
    /// that `τ = 1` (`dt = dx² / (6 ν)`).  Use [`Self::with_lattice_scaling`]
    /// to override it.  Non-positive `debye_length` or `viscosity` yield a
    /// zero scaling that [`init_electrokinetic`] rejects.
    ///
    /// # Arguments
    /// * `permittivity` — fluid permittivity (F m⁻¹)
    /// * `viscosity` — dynamic viscosity (Pa·s)
    /// * `zeta_potential` — surface zeta potential (V)
    /// * `debye_length` — Debye screening length (m)
    /// * `external_field` — applied field `[Ex, Ey]` (V m⁻¹)
    pub fn new(
        permittivity: f64,
        viscosity: f64,
        zeta_potential: f64,
        debye_length: f64,
        external_field: [f64; 2],
    ) -> Self {
        let rho0 = Self::DEFAULT_RHO0;
        let dx = if debye_length > 0.0 {
            debye_length / Self::DEFAULT_NODES_PER_DEBYE_LENGTH
        } else {
            0.0
        };
        let nu = viscosity / rho0;
        let dt = if nu > 0.0 { dx * dx / (6.0 * nu) } else { 0.0 };
        Self {
            permittivity,
            viscosity,
            zeta_potential,
            debye_length,
            external_field,
            dx,
            dt,
            rho0,
        }
    }

    /// Override the lattice scaling (`dx` in m, `dt` in s, `rho0` in kg m⁻³).
    pub fn with_lattice_scaling(mut self, dx: f64, dt: f64, rho0: f64) -> Self {
        self.dx = dx;
        self.dt = dt;
        self.rho0 = rho0;
        self
    }

    /// Kinematic viscosity ν = μ/ρ₀ (m² s⁻¹).
    pub fn kinematic_viscosity(&self) -> f64 {
        self.viscosity / self.rho0
    }

    /// BGK relaxation time `τ = 1/2 + 3 ν dt/dx²` (lattice units).
    pub fn relaxation_time(&self) -> f64 {
        0.5 + 3.0 * self.kinematic_viscosity() * self.dt / (self.dx * self.dx)
    }

    /// Inverse Debye length κ = 1/λ_D (m⁻¹).
    pub fn kappa(&self) -> f64 {
        1.0 / self.debye_length
    }

    /// Validate the parameters and return the relaxation time.
    ///
    /// # Errors
    /// Non-finite values, non-positive `permittivity`, `viscosity`,
    /// `debye_length`, `dx`, `dt` or `rho0`, or `τ ≤ 1/2`.
    pub fn validate(&self) -> Result<f64> {
        let positive = [
            ("permittivity", self.permittivity),
            ("viscosity", self.viscosity),
            ("debye_length", self.debye_length),
            ("dx", self.dx),
            ("dt", self.dt),
            ("rho0", self.rho0),
        ];
        for (name, value) in positive {
            if !value.is_finite() || value <= 0.0 {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min: f64::MIN_POSITIVE,
                    max: f64::MAX,
                });
            }
        }
        let finite = [
            ("zeta_potential", self.zeta_potential),
            ("external_field[0]", self.external_field[0]),
            ("external_field[1]", self.external_field[1]),
        ];
        for (name, value) in finite {
            if !value.is_finite() {
                return Err(Error::NonFiniteField {
                    field: name.to_string(),
                    index: 0,
                    value,
                });
            }
        }
        let tau = self.relaxation_time();
        Error::check_tau(tau)?;
        Ok(tau)
    }
}

// ---------------------------------------------------------------------------
// ElectrokineticLBM
// ---------------------------------------------------------------------------

/// Electrokinetic LBM simulation state on a 2D D2Q9 slit-channel grid.
#[derive(Debug, Clone)]
pub struct ElectrokineticLBM {
    /// Grid dimension in x (periodic).
    pub nx: usize,
    /// Grid dimension in y (walls below row 0 and above row `ny − 1`).
    pub ny: usize,
    /// D2Q9 velocity distribution functions, length `nx × ny × 9`.
    pub f: Vec<f64>,
    /// Net charge density ρ_e (C m⁻³), length `nx × ny`.
    pub charge_density: Vec<f64>,
    /// Electric double-layer potential φ (V), length `nx × ny`.
    pub potential: Vec<f64>,
    /// Body force of the last step in lattice units `[Fx, Fy]`, length `nx × ny`.
    pub body_force: Vec<[f64; 2]>,
    /// LBM relaxation time τ (lattice units).
    pub tau: f64,
    /// Lattice spacing Δx (m).
    pub dx: f64,
    /// Time step Δt (s).
    pub dt: f64,
}

impl ElectrokineticLBM {
    /// Flat node index `j * nx + i`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Flat distribution-function index `(j * nx + i) * 9 + q`.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }

    /// D2Q9 equilibrium distribution at (rho, ux, uy) for direction `q`.
    #[inline]
    pub fn feq(&self, rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
        let cx = D2Q9_C[q][0];
        let cy = D2Q9_C[q][1];
        let cu = cx * ux + cy * uy;
        let u2 = ux * ux + uy * uy;
        D2Q9_W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }

    /// Lattice density and velocity `(ρ, u_x, u_y)` of a node, with the Guo
    /// half-force correction `u = (Σ c f + F/2)/ρ`.
    pub fn macroscopic(&self, i: usize, j: usize) -> (f64, f64, f64) {
        let cell = self.idx(i, j);
        let force = self.body_force.get(cell).copied().unwrap_or([0.0; 2]);
        self.f
            .get(cell * 9..cell * 9 + 9)
            .map_or((0.0, 0.0, 0.0), |node| node_macroscopic(node, force))
    }

    /// Physical fluid velocity `[u_x, u_y]` (m s⁻¹) of a node.
    pub fn velocity(&self, i: usize, j: usize) -> [f64; 2] {
        let (_, ux, uy) = self.macroscopic(i, j);
        let scale = self.dx / self.dt;
        [ux * scale, uy * scale]
    }
}

/// `(ρ, u_x, u_y)` from one node's populations with half-force correction.
fn node_macroscopic(f: &[f64], force: [f64; 2]) -> (f64, f64, f64) {
    let mut rho = 0.0_f64;
    let mut jx = 0.0_f64;
    let mut jy = 0.0_f64;
    for (fq, c) in f.iter().zip(D2Q9_C.iter()) {
        rho += fq;
        jx += c[0] * fq;
        jy += c[1] * fq;
    }
    if rho > 1e-20 {
        (
            rho,
            (jx + 0.5 * force[0]) / rho,
            (jy + 0.5 * force[1]) / rho,
        )
    } else {
        (rho, 0.0, 0.0)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Linearized Poisson–Boltzmann potential profile.
///
/// In the Debye–Hückel limit: φ(x) = ζ exp(−x / λ_D).
///
/// # Arguments
/// * `potential` — distance from the surface (m) \[used as x\]
/// * `zeta`      — surface zeta potential (V)
/// * `debye`     — Debye screening length (m)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::poisson_boltzmann_eq;
/// let phi = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
/// assert!((phi - (-0.05)).abs() < 1e-12);
/// ```
pub fn poisson_boltzmann_eq(potential: f64, zeta: f64, debye: f64) -> f64 {
    if debye <= 0.0 {
        return 0.0;
    }
    zeta * (-potential / debye).exp()
}

/// Helmholtz–Smoluchowski electroosmotic flow velocity.
///
/// v_eo = −ε ζ E / μ
///
/// # Arguments
/// * `zeta`        — zeta potential (V)
/// * `field`       — applied electric field (V m⁻¹)
/// * `viscosity`   — dynamic viscosity (Pa·s)
/// * `permittivity`— fluid permittivity (F m⁻¹)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::electroosmotic_flow_velocity;
/// let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
/// assert!(v > 0.0);
/// ```
pub fn electroosmotic_flow_velocity(
    zeta: f64,
    field: f64,
    viscosity: f64,
    permittivity: f64,
) -> f64 {
    if viscosity == 0.0 {
        return 0.0;
    }
    -permittivity * zeta * field / viscosity
}

/// Initialize an [`ElectrokineticLBM`] slit channel at rest with unit density.
///
/// The double-layer potential and charge density are solved immediately
/// (see [`solve_potential`]), so the returned state is electrostatically
/// consistent.
///
/// # Arguments
/// * `nx`, `ny` — grid dimensions (x periodic, walls bounding y)
/// * `params`   — electrokinetic parameters
///
/// # Errors
/// Invalid grid dimensions, invalid parameters (see
/// [`ElectrokineticParams::validate`]) or a non-converging potential solve.
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::{ElectrokineticParams, init_electrokinetic};
/// let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
/// let state = init_electrokinetic(8, 6, &p).expect("valid parameters");
/// assert_eq!(state.nx, 8);
/// assert_eq!(state.ny, 6);
/// assert!((state.tau - 1.0).abs() < 1e-12);
/// ```
pub fn init_electrokinetic(
    nx: usize,
    ny: usize,
    params: &ElectrokineticParams,
) -> Result<ElectrokineticLBM> {
    Error::check_grid(nx, ny)?;
    let tau = params.validate()?;
    let ncells = nx * ny;
    let mut f = vec![0.0_f64; ncells * 9];
    for node in f.chunks_exact_mut(9) {
        node.copy_from_slice(&D2Q9_W);
    }
    let mut state = ElectrokineticLBM {
        nx,
        ny,
        f,
        charge_density: vec![0.0; ncells],
        potential: vec![0.0; ncells],
        body_force: vec![[0.0; 2]; ncells],
        tau,
        dx: params.dx,
        dt: params.dt,
    };
    solve_potential(&mut state, params)?;
    Ok(state)
}

/// Check that every field of `state` has the size implied by `nx × ny`.
fn check_state_sizes(state: &ElectrokineticLBM) -> Result<()> {
    Error::check_grid(state.nx, state.ny)?;
    let ncells = state.nx * state.ny;
    let checks = [
        ("f", ncells * 9, state.f.len()),
        ("potential", ncells, state.potential.len()),
        ("charge_density", ncells, state.charge_density.len()),
        ("body_force", ncells, state.body_force.len()),
    ];
    for (field, expected, found) in checks {
        if expected != found {
            return Err(Error::ComponentCountMismatch {
                expected,
                found,
                field: field.to_string(),
            });
        }
    }
    Ok(())
}

/// Potential of the neighbours `(south, north)` of row `j`, using the
/// ghost value `2ζ − φ_P` beyond a wall so that `φ = ζ` half-way between.
#[inline]
fn vertical_neighbours(state: &ElectrokineticLBM, i: usize, j: usize, zeta: f64) -> (f64, f64) {
    let phi_p = state.potential[state.idx(i, j)];
    let south = if j == 0 {
        2.0 * zeta - phi_p
    } else {
        state.potential[state.idx(i, j - 1)]
    };
    let north = if j + 1 == state.ny {
        2.0 * zeta - phi_p
    } else {
        state.potential[state.idx(i, j + 1)]
    };
    (south, north)
}

/// Solve the linearised Poisson–Boltzmann (Debye–Hückel) problem
/// `∇²φ = κ² φ` with `φ = ζ` on both walls by SOR, warm-started from the
/// current `state.potential`, and update `state.charge_density = −ε κ² φ`.
///
/// Returns the number of SOR sweeps performed (1 when already converged).
///
/// # Errors
/// Invalid parameters, inconsistent state sizes, or no convergence within
/// the sweep limit ([`Error::ConvergenceFailure`]).
pub fn solve_potential(
    state: &mut ElectrokineticLBM,
    params: &ElectrokineticParams,
) -> Result<usize> {
    params.validate()?;
    check_state_sizes(state)?;
    let (nx, ny) = (state.nx, state.ny);
    let zeta = params.zeta_potential;
    let k2 = (params.kappa() * params.dx).powi(2);
    // Relative tolerance on the larger of the boundary value and the
    // warm-start field, so that switching ζ to zero on a charged state still
    // converges in a bounded number of sweeps.  A zero scale means the exact
    // solution φ ≡ 0 is already reached after the first sweep.
    let phi_scale = state
        .potential
        .iter()
        .fold(zeta.abs(), |acc, phi| acc.max(phi.abs()));
    let tol = POISSON_TOL * phi_scale;
    let mut max_delta = f64::INFINITY;
    for sweep in 1..=POISSON_MAX_SWEEPS {
        max_delta = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let mut sum = state.potential[state.idx(ip, j)] + state.potential[state.idx(im, j)];
                let mut diag = 4.0 + k2;
                if j == 0 {
                    sum += 2.0 * zeta;
                    diag += 1.0;
                } else {
                    sum += state.potential[state.idx(i, j - 1)];
                }
                if j + 1 == ny {
                    sum += 2.0 * zeta;
                    diag += 1.0;
                } else {
                    sum += state.potential[state.idx(i, j + 1)];
                }
                let cell = state.idx(i, j);
                let old = state.potential[cell];
                let new = old + SOR_OMEGA * (sum / diag - old);
                state.potential[cell] = new;
                max_delta = max_delta.max((new - old).abs());
            }
        }
        if max_delta <= tol {
            let eps_k2 = params.permittivity * params.kappa().powi(2);
            for (rho_e, phi) in state.charge_density.iter_mut().zip(state.potential.iter()) {
                *rho_e = -eps_k2 * phi;
            }
            return Ok(sweep);
        }
    }
    Err(Error::ConvergenceFailure {
        residual: max_delta,
        iterations: POISSON_MAX_SWEEPS,
        tolerance: tol,
    })
}

/// Advance the electrokinetic simulation by one time step.
///
/// 1. Re-solve the double-layer potential (warm start; one sweep once
///    converged) and the charge density `ρ_e = −ε κ² φ`.
/// 2. Coulomb body force `F = ρ_e (E_ext − ∇φ)`, converted to lattice units
///    `F_lu = F dt²/(ρ₀ dx)`.
/// 3. BGK collision with Guo forcing: `u = (Σ c f + F/2)/ρ` and source
///    `(1 − 1/(2τ)) w_i [(c_i − u)/c_s² + (c_i·u) c_i/c_s⁴]·F`.
/// 4. Streaming, periodic in x, half-way bounce-back at the y walls.
///
/// # Errors
/// Invalid parameters or state sizes, or a non-converging potential solve.
pub fn step_electrokinetic(
    state: &mut ElectrokineticLBM,
    params: &ElectrokineticParams,
) -> Result<()> {
    let tau = params.validate()?;
    check_state_sizes(state)?;
    state.tau = tau;
    state.dx = params.dx;
    state.dt = params.dt;
    solve_potential(state, params)?;

    let (nx, ny) = (state.nx, state.ny);
    let zeta = params.zeta_potential;
    let force_scale = params.dt * params.dt / (params.rho0 * params.dx);
    let inv_2dx = 0.5 / params.dx;
    for j in 0..ny {
        for i in 0..nx {
            let cell = state.idx(i, j);
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let dphi_dx =
                (state.potential[state.idx(ip, j)] - state.potential[state.idx(im, j)]) * inv_2dx;
            let (south, north) = vertical_neighbours(state, i, j, zeta);
            let dphi_dy = (north - south) * inv_2dx;
            let ex = params.external_field[0] - dphi_dx;
            let ey = params.external_field[1] - dphi_dy;
            let rho_e = state.charge_density[cell];
            state.body_force[cell] = [rho_e * ex * force_scale, rho_e * ey * force_scale];
        }
    }

    // Collision (BGK + Guo forcing)
    let omega = 1.0 / tau;
    let source_prefactor = 1.0 - 0.5 * omega;
    let mut f_post = vec![0.0_f64; nx * ny * 9];
    for cell in 0..nx * ny {
        let base = cell * 9;
        let force = state.body_force[cell];
        let (rho, ux, uy) = node_macroscopic(&state.f[base..base + 9], force);
        for q in 0..9 {
            let [cx, cy] = D2Q9_C[q];
            let cu = cx * ux + cy * uy;
            let source = source_prefactor
                * D2Q9_W[q]
                * (((cx - ux) / CS2 + cu * cx / (CS2 * CS2)) * force[0]
                    + ((cy - uy) / CS2 + cu * cy / (CS2 * CS2)) * force[1]);
            let feq = state.feq(rho, ux, uy, q);
            f_post[base + q] = state.f[base + q] - omega * (state.f[base + q] - feq) + source;
        }
    }

    // Streaming: periodic in x, half-way bounce-back on the y walls.
    for j in 0..ny {
        for i in 0..nx {
            let src = (j * nx + i) * 9;
            for q in 0..9 {
                let val = f_post[src + q];
                let nj = j as isize + D2Q9_C[q][1] as isize;
                if nj < 0 || nj >= ny as isize {
                    state.f[src + D2Q9_OPP[q]] = val;
                } else {
                    let ni = (i as isize + D2Q9_C[q][0] as isize).rem_euclid(nx as isize) as usize;
                    state.f[(nj as usize * nx + ni) * 9 + q] = val;
                }
            }
        }
    }
    Ok(())
}

/// Electric current density magnitude carried by convection, |ρ_e **u**|
/// (C m⁻² s⁻¹), at every cell.
///
/// Uses the physical fluid velocity of [`ElectrokineticLBM::velocity`]
/// (including the Guo half-force correction).  Returns a vector of length
/// `nx × ny`.
pub fn compute_charge_flux(state: &ElectrokineticLBM) -> Vec<f64> {
    let mut flux = vec![0.0_f64; state.nx * state.ny];
    for j in 0..state.ny {
        for i in 0..state.nx {
            let cell = state.idx(i, j);
            let u = state.velocity(i, j);
            let rho_e = state.charge_density.get(cell).copied().unwrap_or(0.0);
            flux[cell] = rho_e.abs() * u[0].hypot(u[1]);
        }
    }
    flux
}

/// Debye–Hückel screened electrostatic potential.
///
/// φ(r) = (charge / (4π ε₀ r)) exp(−r / λ_D)
///
/// Using a dimensionless convention (divided by 4π ε₀):
/// φ(r) = charge · exp(−r / λ_D) / r
///
/// # Arguments
/// * `r`           — radial distance (m)
/// * `charge`      — source charge (C or dimensionless)
/// * `debye_length`— Debye screening length (m)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::debye_huckel_potential;
/// let phi = debye_huckel_potential(1e-9, 1.0, 1e-8);
/// assert!(phi > 0.0);
/// ```
pub fn debye_huckel_potential(r: f64, charge: f64, debye_length: f64) -> f64 {
    if r <= 0.0 || debye_length <= 0.0 {
        return 0.0;
    }
    charge * (-r / debye_length).exp() / r
}

/// Electrophoretic mobility (Helmholtz–Smoluchowski).
///
/// μ_e = ε ζ / μ
///
/// # Arguments
/// * `zeta`        — zeta potential (V)
/// * `viscosity`   — dynamic viscosity (Pa·s)
/// * `permittivity`— fluid permittivity (F m⁻¹)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::electrophoretic_mobility;
/// let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
/// assert!(mob < 0.0);
/// ```
pub fn electrophoretic_mobility(zeta: f64, viscosity: f64, permittivity: f64) -> f64 {
    if viscosity == 0.0 {
        return 0.0;
    }
    permittivity * zeta / viscosity
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── poisson_boltzmann_eq ─────────────────────────────────────────────

    #[test]
    fn test_pb_eq_at_zero_equals_zeta() {
        let phi = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
        assert!((phi - (-0.05)).abs() < 1e-12);
    }

    #[test]
    fn test_pb_eq_decays_with_distance() {
        let phi0 = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
        let phi1 = poisson_boltzmann_eq(1e-8, -0.05, 1e-8);
        assert!(phi1.abs() < phi0.abs());
    }

    #[test]
    fn test_pb_eq_zero_debye_returns_zero() {
        let phi = poisson_boltzmann_eq(1.0, -0.05, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_pb_eq_e_fold_at_debye_length() {
        let zeta = -0.05;
        let debye = 1e-8;
        let phi = poisson_boltzmann_eq(debye, zeta, debye);
        let expected = zeta / std::f64::consts::E;
        assert!((phi - expected).abs() < 1e-14);
    }

    #[test]
    fn test_pb_eq_sign_preserved() {
        let phi = poisson_boltzmann_eq(5e-9, -0.05, 1e-8);
        assert!(phi < 0.0);
    }

    #[test]
    fn test_pb_eq_positive_zeta() {
        let phi = poisson_boltzmann_eq(0.0, 0.05, 1e-8);
        assert!((phi - 0.05).abs() < 1e-12);
    }

    // ── electroosmotic_flow_velocity ─────────────────────────────────────

    #[test]
    fn test_eoflow_negative_zeta_positive_e() {
        // Negative zeta, positive field → positive velocity
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_eoflow_positive_zeta_positive_e() {
        let v = electroosmotic_flow_velocity(0.05, 1e4, 1e-3, 7.1e-10);
        assert!(v < 0.0, "v={v}");
    }

    #[test]
    fn test_eoflow_zero_viscosity() {
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 0.0, 7.1e-10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eoflow_zero_field() {
        let v = electroosmotic_flow_velocity(-0.05, 0.0, 1e-3, 7.1e-10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eoflow_helmholtz_smoluchowski_value() {
        // v = -ε ζ E / μ = -(7.1e-10)*(-0.05)*(1e4)/(1e-3)
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let expected = 7.1e-10 * 0.05 * 1e4 / 1e-3;
        assert!((v - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_eoflow_linear_in_field() {
        let v1 = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let v2 = electroosmotic_flow_velocity(-0.05, 2e4, 1e-3, 7.1e-10);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_eoflow_linear_in_zeta() {
        let v1 = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let v2 = electroosmotic_flow_velocity(-0.10, 1e4, 1e-3, 7.1e-10);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    // ── init_electrokinetic ──────────────────────────────────────────────

    #[test]
    fn test_init_dimensions() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
        let state = init_electrokinetic(8, 6, &p).expect("init");
        assert_eq!(state.nx, 8);
        assert_eq!(state.ny, 6);
        assert_eq!(state.charge_density.len(), 48);
        assert_eq!(state.potential.len(), 48);
        assert_eq!(state.f.len(), 48 * 9);
    }

    #[test]
    fn test_init_unit_density() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p).expect("init");
        for j in 0..4 {
            for i in 0..4 {
                let rho: f64 = (0..9).map(|q| state.f[state.fi(i, j, q)]).sum();
                assert!((rho - 1.0).abs() < 1e-12, "rho at ({i},{j})={rho}");
            }
        }
    }

    #[test]
    fn test_init_solves_double_layer() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 40, &p).expect("init");
        // φ is largest in magnitude next to the walls and decays to the centre.
        let near_wall = state.potential[state.idx(0, 0)];
        let centre = state.potential[state.idx(0, 20)];
        assert!(near_wall < centre && centre < 0.0, "{near_wall} {centre}");
        // Debye–Hückel charge density ρ_e = −ε κ² φ (counter-ions: ρ_e > 0 for ζ < 0).
        let eps_k2 = 7.1e-10 / (1e-8 * 1e-8);
        for (rho_e, phi) in state.charge_density.iter().zip(state.potential.iter()) {
            assert!((rho_e + eps_k2 * phi).abs() <= 1e-12 * eps_k2 * 0.05);
            assert!(*rho_e > 0.0);
        }
    }

    // ── step_electrokinetic ──────────────────────────────────────────────

    #[test]
    fn test_step_conserves_mass() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e3, 0.0]);
        let mut state = init_electrokinetic(8, 8, &p).expect("init");
        let m_before: f64 = state.f.iter().sum();
        step_electrokinetic(&mut state, &p).expect("step");
        let m_after: f64 = state.f.iter().sum();
        assert!(
            (m_after - m_before).abs() / m_before < 1e-10,
            "m_before={m_before} m_after={m_after}"
        );
    }

    #[test]
    fn test_step_finite_values() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let mut state = init_electrokinetic(4, 4, &p).expect("init");
        step_electrokinetic(&mut state, &p).expect("step");
        assert!(state.f.iter().all(|&x| x.is_finite()));
    }

    #[test]
    fn test_step_multiple_iterations_finite() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e2, 0.0]);
        let mut state = init_electrokinetic(8, 8, &p).expect("init");
        for _ in 0..20 {
            step_electrokinetic(&mut state, &p).expect("step");
        }
        assert!(state.f.iter().all(|&x| x.is_finite()));
    }

    // ── compute_charge_flux ──────────────────────────────────────────────

    #[test]
    fn test_charge_flux_zero_at_rest() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p).expect("init");
        let flux = compute_charge_flux(&state);
        assert!(flux.iter().all(|&f| f == 0.0));
    }

    #[test]
    fn test_charge_flux_length() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(6, 5, &p).expect("init");
        let flux = compute_charge_flux(&state);
        assert_eq!(flux.len(), 30);
    }

    #[test]
    fn test_charge_flux_nonnegative() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e3, 0.0]);
        let mut state = init_electrokinetic(4, 4, &p).expect("init");
        for v in state.charge_density.iter_mut() {
            *v = 0.1;
        }
        step_electrokinetic(&mut state, &p).expect("step");
        let flux = compute_charge_flux(&state);
        assert!(flux.iter().all(|&f| f >= 0.0));
    }

    // ── debye_huckel_potential ───────────────────────────────────────────

    #[test]
    fn test_dh_positive_charge() {
        let phi = debye_huckel_potential(1e-9, 1.0, 1e-8);
        assert!(phi > 0.0);
    }

    #[test]
    fn test_dh_decays_with_r() {
        let phi1 = debye_huckel_potential(1e-9, 1.0, 1e-8);
        let phi2 = debye_huckel_potential(2e-9, 1.0, 1e-8);
        assert!(phi2 < phi1);
    }

    #[test]
    fn test_dh_zero_r_returns_zero() {
        let phi = debye_huckel_potential(0.0, 1.0, 1e-8);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_dh_zero_debye_returns_zero() {
        let phi = debye_huckel_potential(1e-9, 1.0, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_dh_negative_charge() {
        let phi = debye_huckel_potential(1e-9, -1.0, 1e-8);
        assert!(phi < 0.0);
    }

    #[test]
    fn test_dh_antisymmetric_charge() {
        let phi_pos = debye_huckel_potential(1e-9, 1.0, 1e-8);
        let phi_neg = debye_huckel_potential(1e-9, -1.0, 1e-8);
        assert!((phi_pos + phi_neg).abs() < 1e-30);
    }

    #[test]
    fn test_dh_finite_value() {
        let phi = debye_huckel_potential(5e-10, 1.602e-19, 3e-9);
        assert!(phi.is_finite());
    }

    // ── electrophoretic_mobility ─────────────────────────────────────────

    #[test]
    fn test_mob_negative_zeta() {
        let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        assert!(mob < 0.0);
    }

    #[test]
    fn test_mob_zero_viscosity() {
        let mob = electrophoretic_mobility(-0.05, 0.0, 7.1e-10);
        assert_eq!(mob, 0.0);
    }

    #[test]
    fn test_mob_formula() {
        let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let expected = 7.1e-10 * (-0.05) / 1e-3;
        assert!((mob - expected).abs() < 1e-30);
    }

    #[test]
    fn test_mob_proportional_to_permittivity() {
        let mob1 = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let mob2 = electrophoretic_mobility(-0.05, 1e-3, 2.0 * 7.1e-10);
        assert!((mob2 / mob1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mob_inversely_proportional_to_viscosity() {
        let mob1 = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let mob2 = electrophoretic_mobility(-0.05, 2e-3, 7.1e-10);
        assert!((mob1 / mob2 - 2.0).abs() < 1e-10);
    }

    // ── ElectrokineticParams ─────────────────────────────────────────────

    #[test]
    fn test_params_construction() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
        assert!((p.permittivity - 7.1e-10).abs() < 1e-20);
        assert!((p.viscosity - 1e-3).abs() < 1e-10);
        assert!((p.zeta_potential - (-0.05)).abs() < 1e-12);
        assert!((p.debye_length - 1e-8).abs() < 1e-18);
        assert!((p.external_field[0] - 1e4).abs() < 1e-6);
        assert_eq!(p.external_field[1], 0.0);
    }

    // ── ElectrokineticLBM helpers ─────────────────────────────────────────

    #[test]
    fn test_ek_idx() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(5, 4, &p).expect("init");
        assert_eq!(state.idx(2, 3), 3 * 5 + 2);
    }

    #[test]
    fn test_ek_fi() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(5, 4, &p).expect("init");
        assert_eq!(state.fi(2, 3, 5), (3 * 5 + 2) * 9 + 5);
    }

    #[test]
    fn test_ek_feq_sum_to_rho() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p).expect("init");
        let rho = 1.5_f64;
        let sum: f64 = (0..9).map(|q| state.feq(rho, 0.1, -0.05, q)).sum();
        assert!((sum - rho).abs() < 1e-12);
    }

    #[test]
    fn test_ek_feq_momentum_x() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p).expect("init");
        let rho = 1.0;
        let ux = 0.1;
        let jx: f64 = (0..9)
            .map(|q| D2Q9_C[q][0] * state.feq(rho, ux, 0.0, q))
            .sum();
        assert!((jx - rho * ux).abs() < 1e-12);
    }

    // ── lattice scaling and relaxation time ──────────────────────────────

    #[test]
    fn test_tau_from_physical_scaling() {
        // ν = μ/ρ₀ = 1e-6 m²/s, dt/dx² = 1e-7/1e-12 → τ = 0.5 + 3·0.1 = 0.8
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0])
            .with_lattice_scaling(1e-6, 1e-7, 1000.0);
        assert!((p.relaxation_time() - 0.8).abs() < 1e-12);
        let state = init_electrokinetic(3, 4, &p).expect("init");
        assert!((state.tau - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_default_scaling_is_stable() {
        // The Pa·s viscosity must not be used directly as a lattice viscosity
        // (that gave τ ≈ 0.503); the default scaling yields τ = 1.
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
        assert!((p.dx - 2.5e-9).abs() < 1e-24);
        assert!((p.relaxation_time() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_invalid_parameters_rejected() {
        let base = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        // 3 ν dt/dx² below the resolution of 0.5 → τ rounds to exactly 1/2.
        let unstable = base.clone().with_lattice_scaling(1e-6, 1e-30, 1000.0);
        assert!(matches!(
            init_electrokinetic(4, 4, &unstable),
            Err(Error::UnstableRelaxationTime { .. })
        ));
        let bad_debye = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 0.0, [0.0, 0.0]);
        assert!(init_electrokinetic(4, 4, &bad_debye).is_err());
        assert!(init_electrokinetic(0, 4, &base).is_err());
        let mut state = init_electrokinetic(4, 4, &base).expect("init");
        state.potential.pop();
        assert!(step_electrokinetic(&mut state, &base).is_err());
    }

    // ── Poisson–Boltzmann solve ──────────────────────────────────────────

    /// Analytical slit-channel double-layer potential with walls half a
    /// spacing outside rows 0 and ny − 1.
    fn analytic_phi(j: usize, ny: usize, p: &ElectrokineticParams) -> f64 {
        let y = j as f64 * p.dx;
        let y_c = 0.5 * (ny as f64 - 1.0) * p.dx;
        let h = 0.5 * ny as f64 * p.dx;
        let k = p.kappa();
        p.zeta_potential * (k * (y - y_c)).cosh() / (k * h).cosh()
    }

    #[test]
    fn test_potential_matches_debye_huckel_slit_solution() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let ny = 40;
        let mut state = init_electrokinetic(3, ny, &p).expect("init");
        for j in 0..ny {
            let exact = analytic_phi(j, ny, &p);
            for i in 0..3 {
                let phi = state.potential[state.idx(i, j)];
                assert!(
                    (phi - exact).abs() < 1e-2 * 0.05,
                    "j={j}: phi={phi} exact={exact}"
                );
            }
        }
        // Warm start: an already converged field needs a single sweep.
        let sweeps = solve_potential(&mut state, &p).expect("solve");
        assert_eq!(sweeps, 1);
    }

    // ── electroosmotic flow ──────────────────────────────────────────────

    #[test]
    fn test_electroosmotic_slit_flow_matches_helmholtz_smoluchowski() {
        let field = 1e5;
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [field, 0.0]);
        let (nx, ny) = (3, 40);
        let mut state = init_electrokinetic(nx, ny, &p).expect("init");
        let u_hs = electroosmotic_flow_velocity(-0.05, field, 1e-3, 7.1e-10);
        let centre = ny / 2;
        let mut prev = 0.0;
        for step in 0..60_000 {
            step_electrokinetic(&mut state, &p).expect("step");
            if step % 500 == 499 {
                let u = state.velocity(0, centre)[0];
                if (u - prev).abs() < 1e-7 * u_hs.abs() {
                    break;
                }
                prev = u;
            }
        }
        // Plug velocity at the channel centre within 5 % of the
        // Helmholtz–Smoluchowski value (the exact centre value differs from
        // it by 1/cosh(κh) ≈ 1.3 %).
        let u_centre = state.velocity(0, centre)[0];
        assert!(u_hs > 0.0 && u_centre > 0.0);
        assert!(
            (u_centre - u_hs).abs() < 0.05 * u_hs,
            "u_centre={u_centre} u_hs={u_hs}"
        );
        // Whole profile against u_HS (1 − φ/ζ).
        for j in 0..ny {
            let exact = u_hs * (1.0 - analytic_phi(j, ny, &p) / p.zeta_potential);
            for i in 0..nx {
                let u = state.velocity(i, j);
                assert!(
                    (u[0] - exact).abs() < 0.01 * u_hs,
                    "j={j}: u={} exact={exact}",
                    u[0]
                );
                // The wall-normal Coulomb force is balanced by the osmotic
                // pressure; half-way bounce-back leaves an O(1 %) residual.
                assert!(u[1].abs() < 1e-2 * u_hs, "j={j}: uy={}", u[1]);
            }
        }
        // Convective current is carried by the counter-ion layer.
        let flux = compute_charge_flux(&state);
        assert!(flux.iter().all(|v| v.is_finite() && *v >= 0.0));
        assert!(flux[state.idx(0, 2)] > flux[state.idx(0, centre)]);
    }

    #[test]
    fn test_reversing_field_reverses_flow() {
        let p_pos = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e5, 0.0]);
        let p_neg = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [-1e5, 0.0]);
        let mut a = init_electrokinetic(3, 20, &p_pos).expect("init");
        let mut b = init_electrokinetic(3, 20, &p_neg).expect("init");
        for _ in 0..200 {
            step_electrokinetic(&mut a, &p_pos).expect("step");
            step_electrokinetic(&mut b, &p_neg).expect("step");
        }
        for j in 0..20 {
            let ua = a.velocity(0, j)[0];
            let ub = b.velocity(0, j)[0];
            assert!(ua > 0.0, "j={j}");
            assert!(
                (ua + ub).abs() < 1e-9 * ua.abs().max(1e-30),
                "j={j}: {ua} {ub}"
            );
        }
    }

    #[test]
    fn test_discharging_walls_relaxes_potential_to_zero() {
        let charged = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let neutral = ElectrokineticParams::new(7.1e-10, 1e-3, 0.0, 1e-8, [0.0, 0.0]);
        let mut state = init_electrokinetic(3, 40, &charged).expect("init");
        assert!(state.potential.iter().any(|phi| phi.abs() > 1e-3));
        let sweeps = solve_potential(&mut state, &neutral).expect("solve");
        assert!(sweeps < 20_000, "took {sweeps} sweeps");
        for (phi, rho_e) in state.potential.iter().zip(state.charge_density.iter()) {
            assert!(phi.abs() < 1e-9, "phi={phi}");
            assert!(rho_e.abs() < 1e-9 * 7.1e-10 * 1e16, "rho_e={rho_e}");
        }
    }

    #[test]
    fn test_no_field_no_zeta_stays_at_rest() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, 0.0, 1e-8, [1e5, 0.0]);
        let mut state = init_electrokinetic(4, 6, &p).expect("init");
        for _ in 0..50 {
            step_electrokinetic(&mut state, &p).expect("step");
        }
        for j in 0..6 {
            for i in 0..4 {
                let u = state.velocity(i, j);
                assert_eq!(u, [0.0, 0.0]);
            }
        }
    }
}
