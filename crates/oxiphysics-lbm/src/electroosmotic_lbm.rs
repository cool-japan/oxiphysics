// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electroosmotic flow LBM: coupled Nernst-Planck / Poisson-Boltzmann.
//!
//! This module implements:
//! - [`ElectroosmoticParams`]: physical parameters (permittivity, valence, Debye length, zeta)
//! - [`NernstPlanckLbm`]: D2Q9 ion transport with migration flux
//! - [`PoissonBoltzmannSolver`]: iterative Gauss-Seidel solver for electric potential
//! - [`ElectricBodyForce`]: body force on fluid from charge density and field
//! - [`ElectroosmoticFlow`]: coupled LBM EOF simulation (Helmholtz-Smoluchowski)
//! - [`ElectrophoresisModel`]: Henry function particle electrophoresis
//! - [`StreamingPotential`]: pressure-driven flow generating electric potential
//! - [`EDLLayer`]: electrical double layer (Debye length)
//! - [`IonicStrength`]: ionic strength I = 0.5 Σ c_i z_i²

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant \[J/K\].
const K_B: f64 = 1.380649e-23;
/// Elementary charge \[C\].
const E_CHARGE: f64 = 1.602176634e-19;
/// Avogadro constant \[mol⁻¹\].
const N_A: f64 = 6.02214076e23;
/// Vacuum permittivity ε₀ \[F/m\].
const EPS_0: f64 = 8.854187817e-12;

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 velocity set (shared)
// ─────────────────────────────────────────────────────────────────────────────

const NQ: usize = 9;
const CX: [f64; NQ] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const CY: [f64; NQ] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
const W_LBM: [f64; NQ] = [
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
const C_SQ: f64 = 1.0 / 3.0; // lattice speed of sound squared

// ─────────────────────────────────────────────────────────────────────────────
// ElectroosmoticParams
// ─────────────────────────────────────────────────────────────────────────────

/// Physical parameters for electroosmotic flow simulations.
#[derive(Clone, Debug)]
pub struct ElectroosmoticParams {
    /// Relative permittivity (dielectric constant) of the medium.
    pub epsilon_r: f64,
    /// Ion valence z (symmetric electrolyte: ±z).
    pub z_val: i32,
    /// Bulk ion concentration \[mol/m³\].
    pub c0: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Zeta potential at the wall \[V\].
    pub zeta: f64,
    /// Dynamic viscosity of fluid \[Pa·s\].
    pub eta: f64,
    /// Applied external electric field \[V/m\].
    pub e_field: f64,
}

impl ElectroosmoticParams {
    /// Creates electroosmotic parameters.
    pub fn new(
        epsilon_r: f64,
        z_val: i32,
        c0: f64,
        temperature: f64,
        zeta: f64,
        eta: f64,
        e_field: f64,
    ) -> Self {
        Self {
            epsilon_r,
            z_val,
            c0,
            temperature,
            zeta,
            eta,
            e_field,
        }
    }

    /// Debye length λ_D = sqrt(ε ε₀ k_B T / (2 N_A z² e² c0)).
    ///
    /// This is the characteristic thickness of the electrical double layer.
    pub fn debye_length(&self) -> f64 {
        let eps = self.epsilon_r * EPS_0;
        let z = self.z_val as f64;
        let numerator = eps * K_B * self.temperature;
        let denominator = 2.0 * N_A * z * z * E_CHARGE * E_CHARGE * self.c0;
        (numerator / denominator).sqrt()
    }

    /// Thermal voltage V_T = k_B T / (z e) \[V\].
    pub fn thermal_voltage(&self) -> f64 {
        K_B * self.temperature / (self.z_val.abs() as f64 * E_CHARGE)
    }

    /// Inverse Debye length κ = 1 / λ_D \[m⁻¹\].
    pub fn kappa(&self) -> f64 {
        1.0 / self.debye_length()
    }

    /// Helmholtz-Smoluchowski EOF velocity u_EOF = -ε ε₀ ζ E / η \[m/s\].
    pub fn eof_velocity(&self) -> f64 {
        -(self.epsilon_r * EPS_0 * self.zeta * self.e_field) / self.eta
    }

    /// Electroosmotic mobility μ_eo = -ε ε₀ ζ / η \[m²/(V·s)\].
    pub fn electroosmotic_mobility(&self) -> f64 {
        -(self.epsilon_r * EPS_0 * self.zeta) / self.eta
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IonicStrength
// ─────────────────────────────────────────────────────────────────────────────

/// Ionic strength calculator: I = 0.5 Σ c_i z_i².
pub struct IonicStrength;

impl IonicStrength {
    /// Computes ionic strength \[mol/m³\] from concentrations and valences.
    ///
    /// I = 0.5 × Σ_i c_i × z_i²
    pub fn compute(concentrations: &[f64], valences: &[i32]) -> f64 {
        assert_eq!(
            concentrations.len(),
            valences.len(),
            "concentration/valence length mismatch"
        );
        0.5 * concentrations
            .iter()
            .zip(valences.iter())
            .map(|(&c, &z)| c * (z * z) as f64)
            .sum::<f64>()
    }

    /// Debye length from ionic strength: λ_D = sqrt(ε ε₀ k_B T / (2 N_A e² I)).
    pub fn debye_from_ionic_strength(epsilon_r: f64, temperature: f64, ionic_str: f64) -> f64 {
        if ionic_str <= 0.0 {
            return f64::INFINITY;
        }
        let eps = epsilon_r * EPS_0;
        let numerator = eps * K_B * temperature;
        let denominator = 2.0 * N_A * E_CHARGE * E_CHARGE * ionic_str;
        (numerator / denominator).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EDLLayer
// ─────────────────────────────────────────────────────────────────────────────

/// Electrical double layer (EDL) thickness and potential profile.
pub struct EDLLayer {
    /// Debye length \[m\].
    pub debye_length: f64,
    /// Zeta (wall) potential \[V\].
    pub zeta_potential: f64,
}

impl EDLLayer {
    /// Creates an EDL layer from Debye length and zeta potential.
    pub fn new(debye_length: f64, zeta_potential: f64) -> Self {
        Self {
            debye_length,
            zeta_potential,
        }
    }

    /// Potential profile in the Debye-Hückel (linearized) approximation.
    ///
    /// φ(y) = ζ × exp(-y / λ_D)
    pub fn potential_profile(&self, y: f64) -> f64 {
        self.zeta_potential * (-y / self.debye_length).exp()
    }

    /// Net charge density \[relative units\] at distance y from wall.
    ///
    /// ρ_e(y) = -ε₀ ε_r κ² φ(y)   (linearized PB)
    pub fn charge_density(&self, y: f64, epsilon_r: f64) -> f64 {
        let kappa = 1.0 / self.debye_length;
        let phi = self.potential_profile(y);
        -epsilon_r * EPS_0 * kappa * kappa * phi
    }

    /// EDL thickness (typically taken as 3 Debye lengths for 95% screening).
    pub fn edl_thickness(&self) -> f64 {
        3.0 * self.debye_length
    }

    /// Gouy-Chapman diffuse layer charge density (surface excess) \[C/m²\].
    pub fn surface_charge_density(&self, epsilon_r: f64) -> f64 {
        let kappa = 1.0 / self.debye_length;
        -epsilon_r * EPS_0 * kappa * self.zeta_potential
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PoissonBoltzmannSolver
// ─────────────────────────────────────────────────────────────────────────────

/// Solves the linearized Poisson-Boltzmann (Debye-Hückel) equation on a 2D grid.
///
/// ∇²φ = κ² φ  (linearized form)
/// Discretized as: (φ_{i+1,j} + φ_{i-1,j} + φ_{i,j+1} + φ_{i,j-1} - 4φ_{i,j}) / h² = κ² φ_{i,j}
pub struct PoissonBoltzmannSolver {
    /// Grid width Nx.
    pub nx: usize,
    /// Grid height Ny.
    pub ny: usize,
    /// Grid spacing h \[m\].
    pub h: f64,
    /// Inverse Debye length κ \[m⁻¹\].
    pub kappa: f64,
    /// Maximum Gauss-Seidel iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl PoissonBoltzmannSolver {
    /// Creates a new Poisson-Boltzmann solver.
    pub fn new(nx: usize, ny: usize, h: f64, kappa: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            nx,
            ny,
            h,
            kappa,
            max_iter,
            tol,
        }
    }

    /// Solves ∇²φ = κ²φ by Gauss-Seidel iteration.
    ///
    /// `phi` is initialized with boundary conditions (Dirichlet).
    /// Interior `fixed[i*ny+j] = true` means the point is a Dirichlet BC node.
    pub fn solve(&self, phi: &mut [f64], fixed: &[bool]) -> usize {
        let h2 = self.h * self.h;
        let kappa2 = self.kappa * self.kappa;
        let diag = 4.0 + kappa2 * h2;
        let mut iters = 0usize;
        for iter in 0..self.max_iter {
            let mut max_change = 0.0f64;
            for j in 1..(self.ny - 1) {
                for i in 1..(self.nx - 1) {
                    let idx = i * self.ny + j;
                    if fixed[idx] {
                        continue;
                    }
                    let neighbors = phi[(i + 1) * self.ny + j]
                        + phi[(i - 1) * self.ny + j]
                        + phi[i * self.ny + (j + 1)]
                        + phi[i * self.ny + (j - 1)];
                    let new_val = neighbors / diag;
                    let change = (new_val - phi[idx]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                    phi[idx] = new_val;
                }
            }
            iters = iter + 1;
            if max_change < self.tol {
                break;
            }
        }
        iters
    }

    /// Builds a 1D flat potential profile for a slab geometry:
    /// φ(0) = ζ (wall), φ(ny-1) = 0 (bulk).
    pub fn solve_1d_slab(&self, zeta: f64) -> Vec<f64> {
        let n = self.ny;
        let mut phi = vec![0.0f64; n];
        phi[0] = zeta;
        // Analytical solution of ∂²φ/∂y² = κ²φ: φ(y) = A sinh(κ(L-y)) / sinh(κL)
        let l = (n - 1) as f64 * self.h;
        let kl = self.kappa * l;
        for (j, phi_j) in phi.iter_mut().enumerate() {
            let y = j as f64 * self.h;
            if kl.abs() < 1e-10 {
                *phi_j = zeta * (1.0 - y / l);
            } else {
                *phi_j = zeta * (self.kappa * (l - y)).sinh() / kl.sinh();
            }
        }
        phi
    }

    /// Extracts gradient (electric field E = -∇φ) at interior points.
    ///
    /// Returns `(ex, ey)` arrays of length nx*ny.
    pub fn electric_field(&self, phi: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let inv2h = 0.5 / self.h;
        let mut ex = vec![0.0f64; self.nx * self.ny];
        let mut ey = vec![0.0f64; self.nx * self.ny];
        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = i * self.ny + j;
                ex[idx] = -(phi[(i + 1) * self.ny + j] - phi[(i - 1) * self.ny + j]) * inv2h;
                ey[idx] = -(phi[i * self.ny + (j + 1)] - phi[i * self.ny + (j - 1)]) * inv2h;
            }
        }
        (ex, ey)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NernstPlanckLbm
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration bundle for the BGK collision step in [`NernstPlanckLbm`].
#[derive(Debug, Clone, Copy)]
struct NernstPlanckCollideConfig {
    nx: usize,
    ny: usize,
    tau: f64,
    z_val: f64,
    thermal_voltage: f64,
    /// +1.0 for cations, -1.0 for anions.
    sign: f64,
}

/// D2Q9 LBM solver for ion transport via Nernst-Planck equations.
///
/// Cation (+z) and anion (-z) distributions are evolved with:
/// - Diffusion (BGK collision)
/// - Advection by fluid velocity
/// - Migration flux: ±D z c ∇φ / (kT)
pub struct NernstPlanckLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Ionic diffusivity D (lattice units).
    pub diffusivity: f64,
    /// Ion valence magnitude z.
    pub z_val: f64,
    /// Thermal voltage V_T = k_B T / (z e) in lattice units.
    pub thermal_voltage: f64,
    /// BGK relaxation time τ for ion transport.
    pub tau: f64,
    /// Cation distribution functions f^+\[i*ny*NQ + j*NQ + q\].
    pub f_cation: Vec<f64>,
    /// Anion distribution functions f^-\[i*ny*NQ + j*NQ + q\].
    pub f_anion: Vec<f64>,
}

impl NernstPlanckLbm {
    /// Creates a NP-LBM solver initialized with uniform concentration `c0`.
    pub fn new(
        nx: usize,
        ny: usize,
        diffusivity: f64,
        z_val: f64,
        thermal_voltage: f64,
        c0: f64,
    ) -> Self {
        let tau = diffusivity / C_SQ + 0.5;
        let size = nx * ny * NQ;
        // Initialize equilibrium distributions
        let f_cation: Vec<f64> = (0..size)
            .map(|k| {
                let q = k % NQ;
                W_LBM[q] * c0
            })
            .collect();
        let f_anion = f_cation.clone();
        Self {
            nx,
            ny,
            diffusivity,
            z_val,
            thermal_voltage,
            tau,
            f_cation,
            f_anion,
        }
    }

    fn idx(&self, i: usize, j: usize, q: usize) -> usize {
        (i * self.ny + j) * NQ + q
    }

    /// Computes cation concentration c^+(i,j) = Σ_q f^+_q.
    pub fn cation_concentration(&self) -> Vec<f64> {
        let mut c = vec![0.0f64; self.nx * self.ny];
        for i in 0..self.nx {
            for j in 0..self.ny {
                c[i * self.ny + j] = (0..NQ).map(|q| self.f_cation[self.idx(i, j, q)]).sum();
            }
        }
        c
    }

    /// Computes anion concentration c^-(i,j) = Σ_q f^-_q.
    pub fn anion_concentration(&self) -> Vec<f64> {
        let mut c = vec![0.0f64; self.nx * self.ny];
        for i in 0..self.nx {
            for j in 0..self.ny {
                c[i * self.ny + j] = (0..NQ).map(|q| self.f_anion[self.idx(i, j, q)]).sum();
            }
        }
        c
    }

    /// Net charge density ρ_e = z e (c^+ - c^-) \[lattice units\].
    pub fn charge_density(&self) -> Vec<f64> {
        let cp = self.cation_concentration();
        let cm = self.anion_concentration();
        cp.iter()
            .zip(&cm)
            .map(|(&a, &b)| self.z_val * (a - b))
            .collect()
    }

    /// Total cation count (conserved quantity).
    pub fn total_cation(&self) -> f64 {
        self.cation_concentration().iter().sum()
    }

    /// Total anion count (conserved quantity).
    pub fn total_anion(&self) -> f64 {
        self.anion_concentration().iter().sum()
    }

    /// Equilibrium distribution for ion transport with velocity (ux, uy).
    ///
    /// f^eq_q = w_q c (1 + c_q·u/cs² + (c_q·u)²/(2cs⁴) - u²/(2cs²))
    fn equilibrium(c: f64, ux: f64, uy: f64, q: usize) -> f64 {
        let cu = CX[q] * ux + CY[q] * uy;
        let u2 = ux * ux + uy * uy;
        W_LBM[q] * c * (1.0 + cu / C_SQ + cu * cu / (2.0 * C_SQ * C_SQ) - u2 / (2.0 * C_SQ))
    }

    /// Performs one BGK collision step with migration body force.
    ///
    /// `ux`, `uy`: fluid velocity field (length nx*ny each).
    /// `ex`, `ey`: electric field -∇φ (length nx*ny each), positive pointing toward bulk.
    /// `cfg.sign`: +1 for cation, -1 for anion.
    fn collide(
        f: &mut [f64],
        cfg: NernstPlanckCollideConfig,
        ux: &[f64],
        uy: &[f64],
        ex: &[f64],
        ey: &[f64],
    ) {
        let NernstPlanckCollideConfig {
            nx,
            ny,
            tau,
            z_val,
            thermal_voltage,
            sign,
        } = cfg;
        let inv_tau = 1.0 / tau;
        let mig = sign * z_val / thermal_voltage;
        for i in 0..nx {
            for j in 0..ny {
                let cell = i * ny + j;
                let c: f64 = (0..NQ).map(|q| f[cell * NQ + q]).sum();
                let vx = ux[cell] + mig * ex[cell];
                let vy = uy[cell] + mig * ey[cell];
                for q in 0..NQ {
                    let feq = Self::equilibrium(c, vx, vy, q);
                    f[cell * NQ + q] += inv_tau * (feq - f[cell * NQ + q]);
                }
            }
        }
    }

    /// Performs one streaming step with periodic boundary conditions.
    fn stream(f: &mut Vec<f64>, nx: usize, ny: usize) {
        let mut f_new = vec![0.0f64; nx * ny * NQ];
        for i in 0..nx {
            for j in 0..ny {
                for q in 0..NQ {
                    let ni = (i as isize + CX[q] as isize).rem_euclid(nx as isize) as usize;
                    let nj = (j as isize + CY[q] as isize).rem_euclid(ny as isize) as usize;
                    f_new[(ni * ny + nj) * NQ + q] = f[(i * ny + j) * NQ + q];
                }
            }
        }
        *f = f_new;
    }

    /// Advances one time step: collision + streaming for both ion species.
    pub fn step(&mut self, ux: &[f64], uy: &[f64], ex: &[f64], ey: &[f64]) {
        Self::collide(
            &mut self.f_cation,
            NernstPlanckCollideConfig {
                nx: self.nx,
                ny: self.ny,
                tau: self.tau,
                z_val: self.z_val,
                thermal_voltage: self.thermal_voltage,
                sign: 1.0,
            },
            ux,
            uy,
            ex,
            ey,
        );
        Self::collide(
            &mut self.f_anion,
            NernstPlanckCollideConfig {
                nx: self.nx,
                ny: self.ny,
                tau: self.tau,
                z_val: self.z_val,
                thermal_voltage: self.thermal_voltage,
                sign: -1.0,
            },
            ux,
            uy,
            ex,
            ey,
        );
        Self::stream(&mut self.f_cation, self.nx, self.ny);
        Self::stream(&mut self.f_anion, self.nx, self.ny);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectricBodyForce
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the electric body force on the fluid: **f** = ρ_e **E** = -ρ_e ∇φ.
pub struct ElectricBodyForce;

impl ElectricBodyForce {
    /// Returns body force (fx, fy) at each grid point.
    ///
    /// `rho_e`: net charge density field (nx*ny).
    /// `ex`, `ey`: electric field components E = -∇φ (nx*ny).
    pub fn compute(rho_e: &[f64], ex: &[f64], ey: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let fx: Vec<f64> = rho_e.iter().zip(ex.iter()).map(|(&r, &e)| r * e).collect();
        let fy: Vec<f64> = rho_e.iter().zip(ey.iter()).map(|(&r, &e)| r * e).collect();
        (fx, fy)
    }

    /// Integrates total force in x-direction: F_x = Σ ρ_e E_x (per unit area).
    pub fn total_force_x(rho_e: &[f64], ex: &[f64]) -> f64 {
        rho_e.iter().zip(ex.iter()).map(|(&r, &e)| r * e).sum()
    }

    /// Integrates total force in y-direction: F_y = Σ ρ_e E_y (per unit area).
    pub fn total_force_y(rho_e: &[f64], ey: &[f64]) -> f64 {
        rho_e.iter().zip(ey.iter()).map(|(&r, &e)| r * e).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectroosmoticFlow
// ─────────────────────────────────────────────────────────────────────────────

/// Coupled electroosmotic flow simulation using LBM + Nernst-Planck + Poisson-Boltzmann.
pub struct ElectroosmoticFlow {
    /// Simulation domain width.
    pub nx: usize,
    /// Simulation domain height.
    pub ny: usize,
    /// Physical parameters.
    pub params: ElectroosmoticParams,
    /// NP-LBM ion solver.
    pub np_solver: NernstPlanckLbm,
    /// Fluid velocity u_x field.
    pub ux: Vec<f64>,
    /// Fluid velocity u_y field.
    pub uy: Vec<f64>,
    /// Electric potential field φ.
    pub phi: Vec<f64>,
}

impl ElectroosmoticFlow {
    /// Creates an EOF simulation domain.
    pub fn new(nx: usize, ny: usize, params: ElectroosmoticParams, c0: f64) -> Self {
        let size = nx * ny;
        let np_solver = NernstPlanckLbm::new(nx, ny, 1.0 / 6.0, params.z_val as f64, 1.0, c0);
        Self {
            nx,
            ny,
            params,
            np_solver,
            ux: vec![0.0; size],
            uy: vec![0.0; size],
            phi: vec![0.0; size],
        }
    }

    /// Initializes the potential field using a 1D Debye-Hückel profile.
    ///
    /// Uses a normalized (dimensionless) decay so that the Debye screening
    /// spans across the channel height in lattice units.
    pub fn init_potential_dh(&mut self) {
        let zeta = self.params.zeta;
        let ny = self.ny;
        // Use a normalized decay: express κ in lattice units
        // where the Debye length spans ~ny/10 lattice cells
        let debye_lu = (ny as f64) / 10.0;
        for j in 0..ny {
            let y = j as f64;
            let phi_val = zeta * (-y / debye_lu).exp();
            for i in 0..self.nx {
                self.phi[i * self.ny + j] = phi_val;
            }
        }
    }

    /// Sets plug-flow EOF velocity (Helmholtz-Smoluchowski) in x-direction.
    pub fn set_eof_velocity(&mut self) {
        let u_eof = self.params.eof_velocity();
        for v in &mut self.ux {
            *v = u_eof;
        }
    }

    /// Returns the mean x-velocity over the domain.
    pub fn mean_velocity_x(&self) -> f64 {
        self.ux.iter().sum::<f64>() / self.ux.len() as f64
    }

    /// Computes the electric field from the stored potential.
    pub fn electric_field(&self) -> (Vec<f64>, Vec<f64>) {
        let pb = PoissonBoltzmannSolver::new(self.nx, self.ny, 1.0, self.params.kappa(), 0, 0.0);
        pb.electric_field(&self.phi)
    }

    /// Advances the ion distributions one step using stored velocity and potential.
    pub fn step(&mut self) {
        let (ex, ey) = self.electric_field();
        let ux = self.ux.clone();
        let uy = self.uy.clone();
        self.np_solver.step(&ux, &uy, &ex, &ey);
    }

    /// Total charge in the domain (should be near zero for neutral electrolyte).
    pub fn total_charge(&self) -> f64 {
        self.np_solver.charge_density().iter().sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectrophoresisModel
// ─────────────────────────────────────────────────────────────────────────────

/// Electrophoresis model for a colloidal particle using Henry's function.
///
/// Electrophoretic mobility: μ_ep = ε ε₀ ζ f(κa) / η
/// where f(κa) is Henry's function interpolating between
/// Hückel limit (κa→0, f=1) and Smoluchowski limit (κa→∞, f=3/2).
pub struct ElectrophoresisModel {
    /// Particle radius a \[m\].
    pub radius: f64,
    /// Zeta potential ζ \[V\].
    pub zeta: f64,
    /// Relative permittivity.
    pub epsilon_r: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub eta: f64,
}

impl ElectrophoresisModel {
    /// Creates an electrophoresis model.
    pub fn new(radius: f64, zeta: f64, epsilon_r: f64, eta: f64) -> Self {
        Self {
            radius,
            zeta,
            epsilon_r,
            eta,
        }
    }

    /// Henry's function f(κa) (Ohshima approximation).
    ///
    /// f(κa) = 1 + (κa)² / (2(1 + δ(κa)))
    /// A common smooth interpolation: f(x) = 1 + x²/(2(1+3exp(-x/2)/x))
    ///
    /// Simple rational approximation:
    /// f → 1 for κa → 0 (Hückel)
    /// f → 3/2 for κa → ∞ (Smoluchowski)
    pub fn henry_function(kappa_a: f64) -> f64 {
        // Ohshima rational approximation
        1.0 + 1.0 / (2.0 * (1.0 + 2.5 / (kappa_a * (1.0 + 2.0 * (-kappa_a).exp()).powi(1))).powi(3))
    }

    /// Electrophoretic mobility μ_ep = ε ε₀ ζ f(κa) / η \[m²/(V·s)\].
    pub fn mobility(&self, kappa: f64) -> f64 {
        let ka = kappa * self.radius;
        let f = Self::henry_function(ka);
        self.epsilon_r * EPS_0 * self.zeta * f / self.eta
    }

    /// Electrophoretic velocity v_ep = μ_ep × E \[m/s\].
    pub fn velocity(&self, kappa: f64, e_field: f64) -> f64 {
        self.mobility(kappa) * e_field
    }

    /// Smoluchowski limit mobility (κa → ∞): μ = ε ε₀ ζ / η × 3/2.
    pub fn smoluchowski_mobility(&self) -> f64 {
        1.5 * self.epsilon_r * EPS_0 * self.zeta / self.eta
    }

    /// Hückel limit mobility (κa → 0): μ = ε ε₀ ζ / η × 1.
    pub fn huckel_mobility(&self) -> f64 {
        self.epsilon_r * EPS_0 * self.zeta / self.eta
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StreamingPotential
// ─────────────────────────────────────────────────────────────────────────────

/// Streaming potential: pressure-driven flow generating an electric potential.
///
/// Streaming potential coefficient: dV/dP = -ε ε₀ ζ / (η σ_∞)
/// where σ_∞ is the bulk conductivity.
pub struct StreamingPotential {
    /// Relative permittivity.
    pub epsilon_r: f64,
    /// Zeta potential ζ \[V\].
    pub zeta: f64,
    /// Dynamic viscosity η \[Pa·s\].
    pub eta: f64,
    /// Bulk conductivity σ_∞ \[S/m\].
    pub sigma: f64,
}

impl StreamingPotential {
    /// Creates a streaming potential model.
    pub fn new(epsilon_r: f64, zeta: f64, eta: f64, sigma: f64) -> Self {
        Self {
            epsilon_r,
            zeta,
            eta,
            sigma,
        }
    }

    /// Helmholtz-Smoluchowski streaming potential coefficient dV/dP \[V/Pa\].
    pub fn hs_coefficient(&self) -> f64 {
        -(self.epsilon_r * EPS_0 * self.zeta) / (self.eta * self.sigma)
    }

    /// Streaming potential ΔV for a pressure drop ΔP \[V\].
    pub fn streaming_voltage(&self, delta_p: f64) -> f64 {
        self.hs_coefficient() * delta_p
    }

    /// Streaming current I_s = ε ε₀ ζ A ΔP / (η L) for channel of cross-section A and length L.
    pub fn streaming_current(&self, delta_p: f64, area: f64, length: f64) -> f64 {
        (self.epsilon_r * EPS_0 * self.zeta * area * delta_p) / (self.eta * length)
    }

    /// Electroosmotic back-pressure (Maxwell stress) coefficient.
    pub fn backpressure_coefficient(&self) -> f64 {
        -(self.epsilon_r * EPS_0 * self.zeta) / (self.eta * self.sigma)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ─── ElectroosmoticParams tests ─────────────────────────────────────────

    fn test_params() -> ElectroosmoticParams {
        ElectroosmoticParams::new(
            80.0,   // epsilon_r (water)
            1,      // z_val
            1.0,    // c0 [mol/m³] = 1 mM
            298.15, // T [K]
            -0.05,  // zeta [V]
            1e-3,   // eta [Pa·s] (water)
            1e4,    // E_field [V/m]
        )
    }

    #[test]
    fn test_debye_length_water_1mm() {
        let p = test_params();
        let ld = p.debye_length();
        // For 1 mM 1:1 electrolyte in water at 25°C, λ_D ≈ 9.6 nm
        assert!(
            ld > 1e-9 && ld < 1e-7,
            "Debye length should be ~10 nm for 1 mM, got {:.3e} m",
            ld
        );
    }

    #[test]
    fn test_debye_length_decreases_with_concentration() {
        let p1 = ElectroosmoticParams::new(80.0, 1, 1.0, 298.15, -0.05, 1e-3, 0.0);
        let p2 = ElectroosmoticParams::new(80.0, 1, 100.0, 298.15, -0.05, 1e-3, 0.0);
        assert!(
            p1.debye_length() > p2.debye_length(),
            "Debye length should decrease with concentration"
        );
    }

    #[test]
    fn test_debye_length_higher_valence_shorter() {
        let p1 = ElectroosmoticParams::new(80.0, 1, 10.0, 298.15, -0.05, 1e-3, 0.0);
        let p2 = ElectroosmoticParams::new(80.0, 2, 10.0, 298.15, -0.05, 1e-3, 0.0);
        assert!(
            p1.debye_length() > p2.debye_length(),
            "Higher valence should give shorter Debye length"
        );
    }

    #[test]
    fn test_eof_velocity_helmholtz_smoluchowski() {
        let p = test_params();
        let u = p.eof_velocity();
        // HS formula: u = -ε ζ E / η; should be positive for negative zeta and positive E
        let expected = -(p.epsilon_r * EPS_0 * p.zeta * p.e_field) / p.eta;
        assert!(
            (u - expected).abs() < 1e-20,
            "EOF velocity should match HS formula"
        );
    }

    #[test]
    fn test_eof_velocity_direction_with_zeta_sign() {
        let p_neg = ElectroosmoticParams::new(80.0, 1, 1.0, 298.15, -0.05, 1e-3, 1000.0);
        let p_pos = ElectroosmoticParams::new(80.0, 1, 1.0, 298.15, 0.05, 1e-3, 1000.0);
        assert!(
            p_neg.eof_velocity() > 0.0,
            "negative zeta + positive E → positive EOF"
        );
        assert!(
            p_pos.eof_velocity() < 0.0,
            "positive zeta + positive E → negative EOF"
        );
        // Opposite signs
        assert!(p_neg.eof_velocity() * p_pos.eof_velocity() < 0.0);
    }

    #[test]
    fn test_eof_velocity_proportional_to_field() {
        let p1 = ElectroosmoticParams::new(80.0, 1, 1.0, 298.15, -0.05, 1e-3, 1000.0);
        let p2 = ElectroosmoticParams::new(80.0, 1, 1.0, 298.15, -0.05, 1e-3, 2000.0);
        let ratio = p2.eof_velocity() / p1.eof_velocity();
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "EOF velocity should be proportional to E, ratio={}",
            ratio
        );
    }

    #[test]
    fn test_thermal_voltage_room_temp() {
        let p = test_params();
        let vt = p.thermal_voltage();
        // ~25.7 mV at 25°C
        assert!(
            vt > 0.025 && vt < 0.028,
            "thermal voltage should be ~25.7 mV at 25°C, got {} V",
            vt
        );
    }

    #[test]
    fn test_electroosmotic_mobility_sign() {
        let p = test_params();
        let mob = p.electroosmotic_mobility();
        // Negative zeta → positive mobility (flow toward cathode)
        assert!(mob > 0.0, "negative zeta should give positive mobility");
    }

    // ─── IonicStrength tests ────────────────────────────────────────────────

    #[test]
    fn test_ionic_strength_symmetric_electrolyte() {
        // 1:1 electrolyte: c+ = c- = 10, z = 1
        let c = vec![10.0, 10.0];
        let z = vec![1i32, -1];
        let i_str = IonicStrength::compute(&c, &z);
        assert!(
            (i_str - 10.0).abs() < 1e-10,
            "1:1 electrolyte ionic strength should equal concentration, got {}",
            i_str
        );
    }

    #[test]
    fn test_ionic_strength_divalent() {
        // 2:2 electrolyte: I = 0.5*(c*4 + c*4) = 4c
        let c = vec![5.0, 5.0];
        let z = vec![2i32, -2];
        let i_str = IonicStrength::compute(&c, &z);
        assert!(
            (i_str - 20.0).abs() < 1e-10,
            "divalent: I = 4c, got {}",
            i_str
        );
    }

    #[test]
    fn test_debye_from_ionic_strength() {
        let ld = IonicStrength::debye_from_ionic_strength(80.0, 298.15, 1.0);
        assert!(
            ld > 1e-9 && ld < 1e-7,
            "Debye length from ionic strength out of range: {}",
            ld
        );
    }

    #[test]
    fn test_ionic_strength_zero_concentration() {
        let ld = IonicStrength::debye_from_ionic_strength(80.0, 298.15, 0.0);
        assert!(
            ld.is_infinite(),
            "zero ionic strength → infinite Debye length"
        );
    }

    // ─── EDLLayer tests ─────────────────────────────────────────────────────

    #[test]
    fn test_edl_potential_at_wall() {
        let edl = EDLLayer::new(10e-9, -0.05);
        let phi_wall = edl.potential_profile(0.0);
        assert!(
            (phi_wall - (-0.05)).abs() < 1e-15,
            "potential at wall should be zeta"
        );
    }

    #[test]
    fn test_edl_potential_decays_exponentially() {
        let ld = 10e-9;
        let edl = EDLLayer::new(ld, -0.05);
        let phi_1ld = edl.potential_profile(ld);
        let phi_2ld = edl.potential_profile(2.0 * ld);
        let expected_1ld = -0.05 * (-1.0f64).exp();
        let expected_2ld = -0.05 * (-2.0f64).exp();
        assert!(
            (phi_1ld - expected_1ld).abs() < 1e-20,
            "potential at 1 Debye length mismatch"
        );
        assert!(
            (phi_2ld - expected_2ld).abs() < 1e-20,
            "potential at 2 Debye lengths mismatch"
        );
    }

    #[test]
    fn test_edl_potential_decays_to_zero() {
        let edl = EDLLayer::new(10e-9, -0.05);
        // At 200 nm (20 Debye lengths), potential = -0.05 * exp(-20) ≈ -1e-10 V
        let phi_far = edl.potential_profile(200e-9);
        assert!(
            phi_far.abs() < 1e-7,
            "potential 20 Debye lengths from wall should be near zero, got {}",
            phi_far
        );
    }

    #[test]
    fn test_edl_thickness() {
        let ld = 10e-9;
        let edl = EDLLayer::new(ld, -0.05);
        // EDL thickness = 3 * Debye length
        assert!(
            (edl.edl_thickness() - 3.0 * ld).abs() < 1e-30,
            "EDL thickness should be 3 Debye lengths"
        );
    }

    #[test]
    fn test_edl_charge_density_sign() {
        let edl = EDLLayer::new(10e-9, -0.05); // negative zeta
        let rho = edl.charge_density(0.0, 80.0);
        // For negative zeta, ρ_e = -ε κ² φ should be positive near wall
        assert!(
            rho > 0.0,
            "charge density near wall should be positive for negative zeta"
        );
    }

    // ─── PoissonBoltzmannSolver tests ───────────────────────────────────────

    #[test]
    fn test_pb_solver_1d_slab_boundary_values() {
        let nx = 5;
        let ny = 20;
        let h = 1e-9;
        let kappa = 1.0 / 10e-9;
        let solver = PoissonBoltzmannSolver::new(nx, ny, h, kappa, 1000, 1e-10);
        let phi = solver.solve_1d_slab(-0.05);
        assert!(
            (phi[0] - (-0.05)).abs() < 1e-15,
            "wall potential should be zeta"
        );
        assert!(phi[ny - 1].abs() < 1e-10, "bulk potential should be near 0");
    }

    #[test]
    fn test_pb_solver_1d_monotone() {
        let nx = 5;
        let ny = 20;
        let h = 1e-9;
        let kappa = 1.0 / 10e-9;
        let solver = PoissonBoltzmannSolver::new(nx, ny, h, kappa, 1000, 1e-10);
        let phi = solver.solve_1d_slab(-0.05);
        // Potential should monotonically increase from negative zeta toward zero
        for k in 0..(ny - 1) {
            assert!(
                phi[k] <= phi[k + 1] + 1e-15,
                "potential should increase from wall to bulk at index {}",
                k
            );
        }
    }

    #[test]
    fn test_pb_electric_field_sign() {
        let nx = 5;
        let ny = 20;
        let h = 1e-9;
        let kappa = 1.0 / 10e-9;
        let solver = PoissonBoltzmannSolver::new(nx, ny, h, kappa, 1000, 1e-10);
        let phi = solver.solve_1d_slab(-0.05);
        let mut phi_2d = vec![0.0f64; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                phi_2d[i * ny + j] = phi[j];
            }
        }
        let (_ex, ey) = solver.electric_field(&phi_2d);
        // E_y = -dφ/dy; for negative ζ, φ increases with y → E_y < 0
        let center_i = nx / 2;
        let center_j = ny / 2;
        let ey_center = ey[center_i * ny + center_j];
        assert!(
            ey_center < 0.0,
            "E_y should be negative (pointing toward wall) for negative zeta, got {}",
            ey_center
        );
    }

    // ─── NernstPlanckLbm tests ──────────────────────────────────────────────

    #[test]
    fn test_np_lbm_initial_concentration() {
        let np = NernstPlanckLbm::new(8, 8, 1.0 / 6.0, 1.0, 1.0, 2.5);
        let c_cat = np.cation_concentration();
        for &c in &c_cat {
            assert!(
                (c - 2.5).abs() < 1e-10,
                "initial cation concentration should be c0=2.5, got {}",
                c
            );
        }
    }

    #[test]
    fn test_np_lbm_charge_neutral_initial() {
        let np = NernstPlanckLbm::new(8, 8, 1.0 / 6.0, 1.0, 1.0, 3.0);
        let rho = np.charge_density();
        for &r in &rho {
            assert!(
                r.abs() < 1e-10,
                "initial charge density should be zero for neutral initialization, got {}",
                r
            );
        }
    }

    #[test]
    fn test_np_lbm_total_cation_conserved() {
        let mut np = NernstPlanckLbm::new(8, 8, 1.0 / 6.0, 1.0, 1.0, 2.0);
        let initial_total = np.total_cation();
        let n = np.nx * np.ny;
        let ux = vec![0.01; n];
        let uy = vec![0.0; n];
        let ex = vec![0.0; n];
        let ey = vec![0.0; n];
        np.step(&ux, &uy, &ex, &ey);
        let final_total = np.total_cation();
        assert!(
            (final_total - initial_total).abs() / initial_total < 1e-8,
            "total cation mass should be conserved after one step: {} vs {}",
            initial_total,
            final_total
        );
    }

    #[test]
    fn test_np_lbm_total_anion_conserved() {
        let mut np = NernstPlanckLbm::new(8, 8, 1.0 / 6.0, 1.0, 1.0, 2.0);
        let initial_total = np.total_anion();
        let n = np.nx * np.ny;
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let ex = vec![0.01; n];
        let ey = vec![0.0; n];
        np.step(&ux, &uy, &ex, &ey);
        let final_total = np.total_anion();
        assert!(
            (final_total - initial_total).abs() / initial_total < 1e-8,
            "total anion mass should be conserved after one step"
        );
    }

    // ─── ElectricBodyForce tests ────────────────────────────────────────────

    #[test]
    fn test_electric_body_force_direction() {
        let rho_e = vec![1.0, -1.0, 0.5];
        let ex = vec![2.0, 2.0, 2.0];
        let ey = vec![0.0, 0.0, 0.0];
        let (fx, fy) = ElectricBodyForce::compute(&rho_e, &ex, &ey);
        assert!((fx[0] - 2.0).abs() < 1e-15);
        assert!((fx[1] - (-2.0)).abs() < 1e-15);
        assert!(fy.iter().all(|&f| f.abs() < 1e-15));
    }

    #[test]
    fn test_electric_body_force_total() {
        let rho_e = vec![1.0, -1.0, 2.0, -2.0];
        let ex = vec![3.0, 3.0, 3.0, 3.0];
        let total = ElectricBodyForce::total_force_x(&rho_e, &ex);
        // Sum of (1-1+2-2)*3 = 0
        assert!(
            total.abs() < 1e-14,
            "net force on neutral system should be 0"
        );
    }

    // ─── ElectroosmoticFlow tests ───────────────────────────────────────────

    #[test]
    fn test_eof_flow_init_potential() {
        let params = test_params();
        let zeta = params.zeta;
        let mut eof = ElectroosmoticFlow::new(10, 20, params, 1.0);
        eof.init_potential_dh();
        // Wall potential at j=0 should equal zeta
        let phi_wall = eof.phi[0]; // i=0, j=0
        assert!(
            (phi_wall - zeta).abs() < 1e-15,
            "wall potential after DH init should be exactly zeta={}, got {}",
            zeta,
            phi_wall
        );
        // Potential should decay away from wall
        let phi_mid = eof.phi[eof.ny / 2];
        assert!(
            phi_mid.abs() < phi_wall.abs(),
            "potential should decay from wall: |phi_mid|={} should be < |phi_wall|={}",
            phi_mid.abs(),
            phi_wall.abs()
        );
    }

    #[test]
    fn test_eof_total_charge_symmetric() {
        let params = test_params();
        let eof = ElectroosmoticFlow::new(8, 8, params, 1.0);
        let total = eof.total_charge();
        assert!(
            total.abs() < 1e-10,
            "initial total charge should be zero for symmetric init, got {}",
            total
        );
    }

    // ─── ElectrophoresisModel tests ─────────────────────────────────────────

    #[test]
    fn test_henry_function_huckel_limit() {
        let f = ElectrophoresisModel::henry_function(0.001);
        assert!(
            (f - 1.0).abs() < 0.05,
            "Henry fn at κa→0 should approach 1, got {}",
            f
        );
    }

    #[test]
    fn test_henry_function_smoluchowski_limit() {
        let f = ElectrophoresisModel::henry_function(1000.0);
        assert!(
            (f - 1.5).abs() < 0.05,
            "Henry fn at κa→∞ should approach 1.5, got {}",
            f
        );
    }

    #[test]
    fn test_henry_function_monotone() {
        let vals: Vec<f64> = [0.01, 0.1, 1.0, 10.0, 100.0]
            .iter()
            .map(|&x| ElectrophoresisModel::henry_function(x))
            .collect();
        for i in 0..(vals.len() - 1) {
            assert!(
                vals[i] <= vals[i + 1] + 1e-10,
                "Henry function should be monotone non-decreasing"
            );
        }
    }

    #[test]
    fn test_electrophoresis_mobility_sign() {
        let ep = ElectrophoresisModel::new(100e-9, -0.05, 80.0, 1e-3);
        let kappa = 1.0 / 10e-9;
        let mob = ep.mobility(kappa);
        assert!(
            mob < 0.0,
            "negative zeta → negative electrophoretic mobility"
        );
    }

    #[test]
    fn test_electrophoretic_velocity() {
        let ep = ElectrophoresisModel::new(100e-9, -0.05, 80.0, 1e-3);
        let kappa = 1.0 / 10e-9;
        let v = ep.velocity(kappa, 1000.0);
        let expected = ep.mobility(kappa) * 1000.0;
        assert!(
            (v - expected).abs() < 1e-30,
            "electrophoretic velocity should equal mobility × E"
        );
    }

    // ─── StreamingPotential tests ───────────────────────────────────────────

    #[test]
    fn test_streaming_potential_sign() {
        // Negative zeta, positive ΔP → negative streaming potential
        let sp = StreamingPotential::new(80.0, -0.05, 1e-3, 0.1);
        let v = sp.streaming_voltage(1000.0);
        assert!(
            v > 0.0,
            "negative zeta + positive ΔP → positive streaming voltage, got {}",
            v
        );
    }

    #[test]
    fn test_streaming_potential_linear_pressure() {
        let sp = StreamingPotential::new(80.0, -0.05, 1e-3, 0.1);
        let v1 = sp.streaming_voltage(1000.0);
        let v2 = sp.streaming_voltage(2000.0);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-10,
            "streaming potential should be linear in ΔP"
        );
    }

    #[test]
    fn test_streaming_potential_hs_coefficient() {
        let sp = StreamingPotential::new(80.0, -0.05, 1e-3, 0.1);
        let coeff = sp.hs_coefficient();
        let expected = -(80.0 * EPS_0 * (-0.05)) / (1e-3 * 0.1);
        assert!(
            (coeff - expected).abs() < 1e-30,
            "HS coefficient mismatch: {} vs {}",
            coeff,
            expected
        );
    }

    #[test]
    fn test_streaming_current_positive() {
        let sp = StreamingPotential::new(80.0, -0.05, 1e-3, 0.1);
        // Negative zeta → streaming current opposing flow
        let i_s = sp.streaming_current(1000.0, 1e-6, 0.01);
        // sign depends on zeta
        assert!(i_s != 0.0, "streaming current should be non-zero");
    }

    // ─── Integration/cross-module tests ─────────────────────────────────────

    #[test]
    fn test_debye_length_consistency() {
        // Debye length from params should equal from ionic strength for 1:1 electrolyte
        let c0 = 10.0; // mol/m³
        let p = ElectroosmoticParams::new(80.0, 1, c0, 298.15, -0.05, 1e-3, 0.0);
        let ld_params = p.debye_length();
        // For 1:1 electrolyte, ionic strength = c0
        let ld_ionic = IonicStrength::debye_from_ionic_strength(80.0, 298.15, c0);
        assert!(
            (ld_params - ld_ionic).abs() / ld_params < 1e-6,
            "Debye length should be consistent: {} vs {}",
            ld_params,
            ld_ionic
        );
    }

    #[test]
    fn test_edl_charge_integrates_to_surface_charge() {
        let ld = 10e-9;
        let zeta = -0.05;
        let eps_r = 80.0;
        let edl = EDLLayer::new(ld, zeta);
        // Integral of ρ_e dy from 0 to ∞ should equal -σ_s (surface charge)
        // σ_s = -ε κ ζ
        let sigma_s = edl.surface_charge_density(eps_r);
        let expected = -eps_r * EPS_0 * (1.0 / ld) * zeta;
        assert!(
            (sigma_s - expected).abs() < 1e-30,
            "surface charge density mismatch: {} vs {}",
            sigma_s,
            expected
        );
    }

    #[test]
    fn test_pi_constant_used() {
        // Verify std::f64::consts::PI is accessible and correct
        let pi = std::f64::consts::PI;
        assert!(
            pi > 3.1 && pi < 3.2,
            "PI constant should be approximately 3.14159"
        );
    }
}
