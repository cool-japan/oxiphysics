#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced porous media LBM: Brinkman–Darcy–Forchheimer with capillary effects.
//!
//! This module extends the basic [`crate::porous_media`] with:
//! - **[`PorousGeometry`]**: porosity field, permeability tensor, tortuosity factor.
//! - **[`BrinkmanModel`]**: effective viscosity and Darcy–Brinkman momentum equation.
//! - **[`DarcyFlow`]**: pressure-driven Darcy velocity computation.
//! - **[`ForcheimerCorrection`]**: nonlinear (inertial) drag at elevated Reynolds numbers.
//! - **[`MultiscalePorous`]**: Representative Elementary Volume (REV) averaging and
//!   upscaling of effective transport properties.
//! - **[`PorousLBM`]**: D2Q9 LBM solver with porosity-modified BGK collision and
//!   Darcy–Brinkman forcing.
//! - **[`DrainageImbibition`]**: capillary pressure curves, wetting/drainage relative
//!   permeability, and Brooks–Corey parameterization.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// D2Q9 lattice constants (private)
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights.
const W9: [f64; 9] = [
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

/// D2Q9 x-velocity components.
const CX9: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 y-velocity components.
const CY9: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 speed-of-sound squared: cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// PorousGeometry
// ---------------------------------------------------------------------------

/// Geometry descriptor for a porous medium.
///
/// Stores the spatially varying porosity ε, permeability tensor K (stored as a
/// symmetric 2×2 matrix in row-major order), and tortuosity τ at each node of a
/// 2-D lattice.
pub struct PorousGeometry {
    /// Number of nodes in x-direction.
    pub nx: usize,
    /// Number of nodes in y-direction.
    pub ny: usize,
    /// Porosity field ε ∈ \[0, 1\] per node.
    pub porosity: Vec<f64>,
    /// Permeability tensor (2×2 row-major) per node, in lattice units².
    pub permeability: Vec<[f64; 4]>,
    /// Tortuosity factor τ ≥ 1 per node (1 = straight channel).
    pub tortuosity: Vec<f64>,
}

impl PorousGeometry {
    /// Create a uniform porous geometry.
    ///
    /// # Arguments
    /// * `nx`, `ny` – grid dimensions
    /// * `eps`      – uniform porosity ε ∈ (0, 1]
    /// * `kxx`      – diagonal permeability K_xx (lattice units²)
    /// * `kyy`      – diagonal permeability K_yy (lattice units²)
    /// * `tau_t`    – uniform tortuosity τ ≥ 1
    pub fn uniform(nx: usize, ny: usize, eps: f64, kxx: f64, kyy: f64, tau_t: f64) -> Self {
        let n = nx * ny;
        let k_tensor = [kxx, 0.0, 0.0, kyy];
        Self {
            nx,
            ny,
            porosity: vec![eps; n],
            permeability: vec![k_tensor; n],
            tortuosity: vec![tau_t; n],
        }
    }

    /// Set porosity at node (ix, iy).
    pub fn set_porosity(&mut self, ix: usize, iy: usize, eps: f64) {
        self.porosity[iy * self.nx + ix] = eps;
    }

    /// Get porosity at node (ix, iy).
    pub fn get_porosity(&self, ix: usize, iy: usize) -> f64 {
        self.porosity[iy * self.nx + ix]
    }

    /// Mean porosity over all nodes.
    pub fn mean_porosity(&self) -> f64 {
        self.porosity.iter().sum::<f64>() / self.porosity.len() as f64
    }

    /// Effective diffusivity coefficient D_eff = D_bulk * eps / tau.
    ///
    /// # Arguments
    /// * `d_bulk` – bulk (free-fluid) diffusivity
    /// * `ix`, `iy` – node indices
    pub fn effective_diffusivity(&self, d_bulk: f64, ix: usize, iy: usize) -> f64 {
        let eps = self.get_porosity(ix, iy);
        let tau = self.tortuosity[iy * self.nx + ix];
        if tau <= 0.0 {
            return 0.0;
        }
        d_bulk * eps / tau
    }

    /// Scalar (isotropic) permeability K_xx at a node.
    pub fn permeability_xx(&self, ix: usize, iy: usize) -> f64 {
        self.permeability[iy * self.nx + ix][0]
    }

    /// Permeability ratio K_xx / K_yy at a node (anisotropy ratio).
    pub fn anisotropy_ratio(&self, ix: usize, iy: usize) -> f64 {
        let k = &self.permeability[iy * self.nx + ix];
        if k[3] == 0.0 {
            f64::INFINITY
        } else {
            k[0] / k[3]
        }
    }
}

// ---------------------------------------------------------------------------
// BrinkmanModel
// ---------------------------------------------------------------------------

/// Brinkman extension of Darcy's law.
///
/// Adds a viscous Laplacian term to the Darcy drag to resolve boundary layers:
///
/// ```text
/// −∇p + μ_eff ∇²u − (μ/K) u = 0
/// ```
///
/// where μ_eff = μ / ε is the effective (Brinkman) viscosity.
pub struct BrinkmanModel {
    /// Fluid dynamic viscosity μ (Pa·s or lattice units).
    pub viscosity: f64,
    /// Local porosity ε ∈ (0, 1].
    pub porosity: f64,
    /// Scalar permeability K (m² or lattice units²).
    pub permeability: f64,
}

impl BrinkmanModel {
    /// Create a [`BrinkmanModel`].
    pub fn new(viscosity: f64, porosity: f64, permeability: f64) -> Self {
        Self {
            viscosity,
            porosity,
            permeability,
        }
    }

    /// Effective (Brinkman) viscosity μ_eff = μ / ε.
    pub fn effective_viscosity(&self) -> f64 {
        if self.porosity <= 0.0 {
            return f64::INFINITY;
        }
        self.viscosity / self.porosity
    }

    /// Darcy drag coefficient α = μ / K (Pa s m⁻²).
    pub fn darcy_drag_coefficient(&self) -> f64 {
        if self.permeability <= 0.0 {
            return f64::INFINITY;
        }
        self.viscosity / self.permeability
    }

    /// BGK relaxation parameter ω for Brinkman–LBM.
    ///
    /// ω = 1 / (3 ν_eff + 0.5),  ν_eff = μ_eff / ρ (ρ = 1 in lattice units).
    pub fn brinkman_omega(&self, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let nu_eff = self.effective_viscosity() / density;
        1.0 / (3.0 * nu_eff + 0.5)
    }

    /// Brinkman-modified Darcy velocity for a 1D pressure gradient dp/dx.
    ///
    /// In the Darcy limit (no Laplacian): u = −(K/μ) dp/dx.
    pub fn darcy_brinkman_velocity_1d(&self, pressure_gradient: f64) -> f64 {
        if self.viscosity <= 0.0 || self.permeability <= 0.0 {
            return 0.0;
        }
        -(self.permeability / self.viscosity) * pressure_gradient
    }
}

// ---------------------------------------------------------------------------
// DarcyFlow
// ---------------------------------------------------------------------------

/// Darcy flow through a porous domain.
///
/// Computes the seepage (Darcy) velocity u = −(K/μ) ∇p and related quantities.
pub struct DarcyFlow {
    /// Dynamic viscosity μ of the fluid.
    pub viscosity: f64,
    /// Scalar permeability K of the medium.
    pub permeability: f64,
    /// Porosity ε for pore velocity correction.
    pub porosity: f64,
}

impl DarcyFlow {
    /// Create a [`DarcyFlow`] model.
    pub fn new(viscosity: f64, permeability: f64, porosity: f64) -> Self {
        Self {
            viscosity,
            permeability,
            porosity,
        }
    }

    /// Darcy (seepage) velocity components (ux, uy) from a pressure gradient.
    ///
    /// u = −(K/μ) ∇p
    ///
    /// # Arguments
    /// * `grad_px` – pressure gradient in x-direction (Pa m⁻¹)
    /// * `grad_py` – pressure gradient in y-direction (Pa m⁻¹)
    pub fn darcy_velocity(&self, grad_px: f64, grad_py: f64) -> (f64, f64) {
        if self.viscosity <= 0.0 || self.permeability < 0.0 {
            return (0.0, 0.0);
        }
        let factor = self.permeability / self.viscosity;
        (-factor * grad_px, -factor * grad_py)
    }

    /// Pore velocity (interstitial velocity): v = u / ε.
    pub fn pore_velocity(&self, grad_px: f64, grad_py: f64) -> (f64, f64) {
        let (ux, uy) = self.darcy_velocity(grad_px, grad_py);
        if self.porosity <= 0.0 {
            return (0.0, 0.0);
        }
        (ux / self.porosity, uy / self.porosity)
    }

    /// Hydraulic conductivity K_h = K ρ g / μ (m/s) with gravitational acceleration g.
    pub fn hydraulic_conductivity(&self, density: f64, gravity: f64) -> f64 {
        if self.viscosity <= 0.0 {
            return 0.0;
        }
        self.permeability * density * gravity / self.viscosity
    }

    /// Darcy number Da = K / L² for characteristic length L.
    pub fn darcy_number(&self, char_length: f64) -> f64 {
        if char_length <= 0.0 {
            return 0.0;
        }
        self.permeability / (char_length * char_length)
    }

    /// Reynolds number for porous flow: Re_p = ρ |u| L / μ.
    pub fn pore_reynolds(&self, density: f64, speed: f64, char_length: f64) -> f64 {
        if self.viscosity <= 0.0 {
            return 0.0;
        }
        density * speed * char_length / self.viscosity
    }
}

// ---------------------------------------------------------------------------
// ForcheimerCorrection
// ---------------------------------------------------------------------------

/// Forchheimer (nonlinear inertial) correction to Darcy's law.
///
/// The extended Darcy–Forchheimer equation in 1D:
///
/// ```text
/// −dp/dx = (μ/K) u + β_F ρ |u| u
/// ```
///
/// where β_F is the Forchheimer coefficient (m⁻¹), also written as
/// C_F / √K with C_F ≈ 0.55 (Ergun constant).
pub struct ForcheimerCorrection {
    /// Scalar permeability K (m²).
    pub permeability: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Forchheimer (inertia) coefficient β_F (m⁻¹).
    pub beta_f: f64,
    /// Fluid density ρ (kg m⁻³).
    pub density: f64,
}

impl ForcheimerCorrection {
    /// Create a [`ForcheimerCorrection`] model.
    pub fn new(permeability: f64, viscosity: f64, beta_f: f64, density: f64) -> Self {
        Self {
            permeability,
            viscosity,
            beta_f,
            density,
        }
    }

    /// Estimate β_F from permeability using the Ergun relation:
    /// β_F = C_F / √K,  C_F = 0.55.
    pub fn ergun_beta(permeability: f64) -> f64 {
        if permeability <= 0.0 {
            return 0.0;
        }
        0.55 / permeability.sqrt()
    }

    /// Pressure gradient for given Darcy velocity u in 1D (Darcy–Forchheimer).
    pub fn pressure_gradient_1d(&self, velocity: f64) -> f64 {
        let darcy_term = (self.viscosity / self.permeability) * velocity;
        let forchheimer_term = self.beta_f * self.density * velocity.abs() * velocity;
        darcy_term + forchheimer_term
    }

    /// Additional drag force per unit volume from inertial effects: F_F = β_F ρ |u| u.
    pub fn inertial_drag(&self, velocity: f64) -> f64 {
        self.beta_f * self.density * velocity.abs() * velocity
    }

    /// Effective (velocity-dependent) permeability K_eff(Re):
    ///
    /// 1/K_eff = 1/K + (β_F ρ |u|) / μ
    pub fn effective_permeability(&self, velocity: f64) -> f64 {
        if self.permeability <= 0.0 || self.viscosity <= 0.0 {
            return 0.0;
        }
        let inertia = self.beta_f * self.density * velocity.abs() / self.viscosity;
        1.0 / (1.0 / self.permeability + inertia)
    }

    /// Forchheimer number Fo = β_F ρ K |u| / μ (ratio of inertial to viscous drag).
    pub fn forchheimer_number(&self, velocity: f64) -> f64 {
        if self.viscosity <= 0.0 {
            return 0.0;
        }
        self.beta_f * self.density * self.permeability * velocity.abs() / self.viscosity
    }
}

// ---------------------------------------------------------------------------
// MultiscalePorous
// ---------------------------------------------------------------------------

/// Multiscale porous media: REV averaging and upscaling of effective properties.
///
/// Performs volume (REV) averaging over a set of micro-scale cells to obtain
/// macro-scale effective transport properties.
pub struct MultiscalePorous {
    /// Microscale porosity values (one per micro-cell).
    pub micro_porosity: Vec<f64>,
    /// Microscale permeability values K_xx (one per micro-cell).
    pub micro_permeability: Vec<f64>,
    /// Microscale tortuosity values (one per micro-cell).
    pub micro_tortuosity: Vec<f64>,
}

impl MultiscalePorous {
    /// Create a [`MultiscalePorous`] model from microscale data.
    pub fn new(
        micro_porosity: Vec<f64>,
        micro_permeability: Vec<f64>,
        micro_tortuosity: Vec<f64>,
    ) -> Self {
        Self {
            micro_porosity,
            micro_permeability,
            micro_tortuosity,
        }
    }

    /// REV-averaged (volume-averaged) porosity.
    pub fn rev_porosity(&self) -> f64 {
        if self.micro_porosity.is_empty() {
            return 0.0;
        }
        self.micro_porosity.iter().sum::<f64>() / self.micro_porosity.len() as f64
    }

    /// Upscaled permeability by arithmetic averaging (parallel flow).
    pub fn upscaled_permeability_parallel(&self) -> f64 {
        if self.micro_permeability.is_empty() {
            return 0.0;
        }
        self.micro_permeability.iter().sum::<f64>() / self.micro_permeability.len() as f64
    }

    /// Upscaled permeability by harmonic averaging (series/perpendicular flow).
    pub fn upscaled_permeability_series(&self) -> f64 {
        if self.micro_permeability.is_empty() {
            return 0.0;
        }
        let n = self.micro_permeability.len() as f64;
        let harm_sum: f64 = self
            .micro_permeability
            .iter()
            .map(|&k| if k > 0.0 { 1.0 / k } else { f64::INFINITY })
            .sum();
        if harm_sum <= 0.0 || harm_sum.is_infinite() {
            return 0.0;
        }
        n / harm_sum
    }

    /// Effective tortuosity factor (geometric mean over REV).
    pub fn effective_tortuosity(&self) -> f64 {
        if self.micro_tortuosity.is_empty() {
            return 1.0;
        }
        let log_sum: f64 = self.micro_tortuosity.iter().map(|&t| t.max(1.0).ln()).sum();
        (log_sum / self.micro_tortuosity.len() as f64).exp()
    }

    /// Effective diffusivity D_eff = D_bulk * ε_REV / τ_eff.
    pub fn effective_diffusivity(&self, d_bulk: f64) -> f64 {
        let eps = self.rev_porosity();
        let tau = self.effective_tortuosity();
        if tau <= 0.0 {
            return 0.0;
        }
        d_bulk * eps / tau
    }

    /// Standard deviation of porosity over the REV (heterogeneity measure).
    pub fn porosity_std_dev(&self) -> f64 {
        if self.micro_porosity.is_empty() {
            return 0.0;
        }
        let mean = self.rev_porosity();
        let var = self
            .micro_porosity
            .iter()
            .map(|&e| (e - mean).powi(2))
            .sum::<f64>()
            / self.micro_porosity.len() as f64;
        var.sqrt()
    }
}

// ---------------------------------------------------------------------------
// PorousLBM
// ---------------------------------------------------------------------------

/// LBM solver for porous media flow (D2Q9, BGK, Darcy–Brinkman forcing).
///
/// The collision term includes an additional Darcy drag body force per unit
/// volume **F** = −(μ/K) ε **u**, applied via the Guo forcing scheme.
pub struct PorousLBM {
    /// Grid x-size.
    pub nx: usize,
    /// Grid y-size.
    pub ny: usize,
    /// Distribution functions f_i per node (9-velocity D2Q9).
    pub f: Vec<[f64; 9]>,
    /// Post-streaming buffer.
    f_buf: Vec<[f64; 9]>,
    /// BGK relaxation frequency ω (1/τ_BGK).
    pub omega: f64,
    /// Porous geometry (porosity + permeability + tortuosity).
    pub geometry: PorousGeometry,
    /// External body-force per unit mass in x (e.g. gravity or pressure gradient).
    pub force_x: f64,
    /// External body-force per unit mass in y.
    pub force_y: f64,
}

impl PorousLBM {
    /// Create a porous LBM domain at rest.
    ///
    /// # Arguments
    /// * `nx`, `ny` – domain size
    /// * `omega`    – BGK relaxation frequency
    /// * `geometry` – porosity / permeability / tortuosity field
    /// * `force_x`, `force_y` – external body-force per unit mass
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        geometry: PorousGeometry,
        force_x: f64,
        force_y: f64,
    ) -> Self {
        let n = nx * ny;
        // Equilibrium at rest with ρ = 1: f_i = w_i
        let f0: [f64; 9] = W9;
        Self {
            nx,
            ny,
            f: vec![f0; n],
            f_buf: vec![[0.0; 9]; n],
            omega,
            geometry,
            force_x,
            force_y,
        }
    }

    /// Compute macroscopic density at node index `n`.
    fn node_density(&self, n: usize) -> f64 {
        self.f[n].iter().sum()
    }

    /// Compute macroscopic velocity at node `n` with porosity-modified Darcy drag.
    fn node_velocity(&self, n: usize) -> (f64, f64) {
        let rho = self.node_density(n);
        if rho == 0.0 {
            return (0.0, 0.0);
        }
        let ix = n % self.nx;
        let iy = n / self.nx;
        let eps = self.geometry.get_porosity(ix, iy);
        let k = self.geometry.permeability[n][0]; // K_xx scalar approximation

        // Darcy–Brinkman: effective velocity accounts for porosity
        let nu = (1.0 / self.omega - 0.5) * CS2;
        let mu = nu * rho;

        // Solve: (mu/K + 0.5 F_drag) * u = momentum / ε
        let drag_coeff = if k > 0.0 { mu / k } else { 0.0 };
        let alpha = eps / (eps + 0.5 * drag_coeff);

        let mom_x: f64 = self.f[n]
            .iter()
            .zip(CX9.iter())
            .map(|(fi, ci)| fi * ci)
            .sum();
        let mom_y: f64 = self.f[n]
            .iter()
            .zip(CY9.iter())
            .map(|(fi, ci)| fi * ci)
            .sum();

        // Add half-step forcing
        let ux = alpha * (mom_x / rho + 0.5 * self.force_x);
        let uy = alpha * (mom_y / rho + 0.5 * self.force_y);
        (ux, uy)
    }

    /// D2Q9 equilibrium distribution at density ρ, velocity (ux, uy).
    fn feq(i: usize, rho: f64, ux: f64, uy: f64) -> f64 {
        let cu = CX9[i] * ux + CY9[i] * uy;
        let u2 = ux * ux + uy * uy;
        W9[i] * rho * (1.0 + cu / CS2 + 0.5 * cu * cu / (CS2 * CS2) - 0.5 * u2 / CS2)
    }

    /// Guo forcing term (accounts for Darcy drag + external force).
    fn guo_force(i: usize, rho: f64, ux: f64, uy: f64, fx: f64, fy: f64, omega: f64) -> f64 {
        let cu = CX9[i] * ux + CY9[i] * uy;
        let u2 = ux * ux + uy * uy;
        let e_dot_f = CX9[i] * fx + CY9[i] * fy;
        let u_dot_f = ux * fx + uy * fy;
        let coeff = (1.0 - 0.5 * omega) * W9[i] / CS2;
        coeff
            * rho
            * (e_dot_f - u_dot_f / CS2 + cu * e_dot_f / CS2 - cu * u_dot_f / (CS2 * CS2)
                + 0.5 * (CX9[i] * CX9[i] - u2 / 3.0) * fx
                + 0.5 * (CY9[i] * CY9[i] - u2 / 3.0) * fy)
    }

    /// Perform a single LBM time step (collision + streaming).
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n_nodes = nx * ny;
        let omega = self.omega;

        // ── Collision ───────────────────────────────────────────────────────
        for node in 0..n_nodes {
            let rho = self.node_density(node);
            if rho == 0.0 {
                continue;
            }
            let (ux, uy) = self.node_velocity(node);
            let ix = node % nx;
            let iy = node / nx;
            let k = self.geometry.permeability[node][0];
            let nu = (1.0 / omega - 0.5) * CS2;
            let mu = nu * rho;
            let drag = if k > 0.0 { mu / k } else { 0.0 };
            // Net body force: external − Darcy drag
            let fx_net = self.force_x - drag * ux / rho;
            let fy_net = self.force_y - drag * uy / rho;
            let _eps = self.geometry.get_porosity(ix, iy);

            for i in 0..9 {
                let feq_i = Self::feq(i, rho, ux, uy);
                let g_i = Self::guo_force(i, rho, ux, uy, fx_net, fy_net, omega);
                self.f_buf[node][i] = self.f[node][i] - omega * (self.f[node][i] - feq_i) + g_i;
            }
        }

        // ── Streaming (periodic BC) ─────────────────────────────────────────
        let mut f_stream = vec![[0.0f64; 9]; n_nodes];
        for y in 0..ny {
            for x in 0..nx {
                let node = y * nx + x;
                for i in 0..9 {
                    let xd = ((x as isize + CX9[i] as isize).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + CY9[i] as isize).rem_euclid(ny as isize)) as usize;
                    f_stream[yd * nx + xd][i] = self.f_buf[node][i];
                }
            }
        }
        self.f = f_stream;
    }

    /// Macroscopic density at (ix, iy).
    pub fn density(&self, ix: usize, iy: usize) -> f64 {
        self.node_density(iy * self.nx + ix)
    }

    /// Macroscopic velocity at (ix, iy).
    pub fn velocity(&self, ix: usize, iy: usize) -> (f64, f64) {
        self.node_velocity(iy * self.nx + ix)
    }

    /// Total mass in the domain.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
}

// ---------------------------------------------------------------------------
// DrainageImbibition
// ---------------------------------------------------------------------------

/// Drainage and imbibition capillary pressure and relative permeability curves.
///
/// Uses the Brooks–Corey model:
///
/// ```text
/// P_c(S_w) = P_e · S_e^(−1/λ_bc)
/// k_rw(S_w) = S_e^((2 + 3λ_bc)/λ_bc)
/// k_rnw(S_w) = (1 − S_e)² (1 − S_e^((2+λ_bc)/λ_bc))
/// ```
///
/// where S_e = (S_w − S_wr) / (1 − S_wr − S_nr) is the effective water saturation.
pub struct DrainageImbibition {
    /// Entry (bubbling) capillary pressure P_e (Pa).
    pub entry_pressure: f64,
    /// Brooks–Corey pore-size distribution index λ.
    pub lambda_bc: f64,
    /// Residual wetting-phase saturation S_wr ∈ \[0, 1).
    pub s_wr: f64,
    /// Residual non-wetting phase saturation S_nr ∈ \[0, 1).
    pub s_nr: f64,
}

impl DrainageImbibition {
    /// Create a [`DrainageImbibition`] model.
    pub fn new(entry_pressure: f64, lambda_bc: f64, s_wr: f64, s_nr: f64) -> Self {
        Self {
            entry_pressure,
            lambda_bc,
            s_wr,
            s_nr,
        }
    }

    /// Effective wetting saturation S_e ∈ \[0, 1\].
    pub fn effective_saturation(&self, s_w: f64) -> f64 {
        let denom = 1.0 - self.s_wr - self.s_nr;
        if denom <= 0.0 {
            return 0.0;
        }
        ((s_w - self.s_wr) / denom).clamp(0.0, 1.0)
    }

    /// Brooks–Corey capillary pressure P_c(S_w) (Pa).
    ///
    /// Returns 0 if S_w ≥ 1 − S_nr (fully saturated).
    pub fn capillary_pressure(&self, s_w: f64) -> f64 {
        let se = self.effective_saturation(s_w);
        if se <= 0.0 {
            return f64::INFINITY;
        }
        if se >= 1.0 {
            return 0.0;
        }
        if self.lambda_bc <= 0.0 {
            return self.entry_pressure;
        }
        self.entry_pressure * se.powf(-1.0 / self.lambda_bc)
    }

    /// Relative permeability of the wetting phase k_rw(S_w).
    pub fn relative_permeability_wetting(&self, s_w: f64) -> f64 {
        let se = self.effective_saturation(s_w);
        if self.lambda_bc <= 0.0 {
            return se;
        }
        se.powf((2.0 + 3.0 * self.lambda_bc) / self.lambda_bc)
    }

    /// Relative permeability of the non-wetting phase k_rnw(S_w).
    pub fn relative_permeability_nonwetting(&self, s_w: f64) -> f64 {
        let se = self.effective_saturation(s_w);
        if self.lambda_bc <= 0.0 {
            return 1.0 - se;
        }
        (1.0 - se).powi(2) * (1.0 - se.powf((2.0 + self.lambda_bc) / self.lambda_bc))
    }

    /// Fractional flow of the wetting phase f_w (Buckley–Leverett).
    ///
    /// f_w = k_rw / (k_rw + μ_w / μ_nw · k_rnw)
    pub fn fractional_flow(&self, s_w: f64, viscosity_ratio: f64) -> f64 {
        let krw = self.relative_permeability_wetting(s_w);
        let krnw = self.relative_permeability_nonwetting(s_w);
        let denom = krw + viscosity_ratio * krnw;
        if denom <= 0.0 {
            return 0.0;
        }
        krw / denom
    }

    /// Drainage curve: capillary pressure at uniformly sampled saturations.
    ///
    /// Returns `n_points` (S_w, P_c) pairs from s_w = s_wr to s_w = 1 − s_nr.
    pub fn drainage_curve(&self, n_points: usize) -> Vec<(f64, f64)> {
        if n_points < 2 {
            return Vec::new();
        }
        let s_min = self.s_wr + 1e-6;
        let s_max = 1.0 - self.s_nr;
        (0..n_points)
            .map(|i| {
                let s_w = s_min + (s_max - s_min) * i as f64 / (n_points - 1) as f64;
                (s_w, self.capillary_pressure(s_w))
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── PorousGeometry ────────────────────────────────────────────────────

    #[test]
    fn test_geometry_uniform_porosity() {
        let g = PorousGeometry::uniform(4, 4, 0.4, 1e-3, 1e-3, 1.5);
        assert!((g.mean_porosity() - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_geometry_get_set_porosity() {
        let mut g = PorousGeometry::uniform(4, 4, 0.4, 1e-3, 1e-3, 1.5);
        g.set_porosity(1, 2, 0.7);
        assert!((g.get_porosity(1, 2) - 0.7).abs() < 1e-12);
    }

    #[test]
    fn test_geometry_effective_diffusivity() {
        let g = PorousGeometry::uniform(4, 4, 0.5, 1e-3, 1e-3, 2.0);
        let d_eff = g.effective_diffusivity(1.0, 0, 0);
        // D_eff = D * eps / tau = 1 * 0.5 / 2 = 0.25
        assert!((d_eff - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_geometry_permeability_xx() {
        let g = PorousGeometry::uniform(4, 4, 0.4, 2e-3, 1e-3, 1.0);
        assert!((g.permeability_xx(0, 0) - 2e-3).abs() < 1e-15);
    }

    #[test]
    fn test_geometry_anisotropy_ratio_isotropic() {
        let g = PorousGeometry::uniform(4, 4, 0.4, 1e-3, 1e-3, 1.0);
        assert!((g.anisotropy_ratio(0, 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_geometry_anisotropy_ratio_2x() {
        let g = PorousGeometry::uniform(4, 4, 0.4, 2e-3, 1e-3, 1.0);
        assert!((g.anisotropy_ratio(0, 0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_geometry_node_count() {
        let g = PorousGeometry::uniform(5, 3, 0.5, 1e-3, 1e-3, 1.0);
        assert_eq!(g.porosity.len(), 15);
    }

    // ── BrinkmanModel ─────────────────────────────────────────────────────

    #[test]
    fn test_brinkman_effective_viscosity() {
        let b = BrinkmanModel::new(0.01, 0.5, 1e-3);
        // mu_eff = 0.01 / 0.5 = 0.02
        assert!((b.effective_viscosity() - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_brinkman_darcy_drag_coeff() {
        let b = BrinkmanModel::new(0.01, 0.5, 2e-3);
        // alpha = mu / K = 0.01 / 2e-3 = 5.0
        assert!((b.darcy_drag_coefficient() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_brinkman_velocity_1d() {
        // u = -(K/mu) * dp/dx = -(1e-3/0.01)*(-10) = 1.0
        let b = BrinkmanModel::new(0.01, 0.5, 1e-3);
        let u = b.darcy_brinkman_velocity_1d(-10.0);
        assert!((u - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_brinkman_zero_permeability() {
        let b = BrinkmanModel::new(0.01, 0.5, 0.0);
        let u = b.darcy_brinkman_velocity_1d(-10.0);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_brinkman_omega_positive() {
        let b = BrinkmanModel::new(0.001, 0.5, 1e-3);
        let om = b.brinkman_omega(1.0);
        assert!(om > 0.0 && om <= 2.0);
    }

    // ── DarcyFlow ─────────────────────────────────────────────────────────

    #[test]
    fn test_darcy_velocity_x() {
        // u_x = -(K/mu)*dp/dx = -(1e-3/0.01)*(-1.0) = 0.1
        let d = DarcyFlow::new(0.01, 1e-3, 0.4);
        let (ux, uy) = d.darcy_velocity(-1.0, 0.0);
        assert!((ux - 0.1).abs() < 1e-12);
        assert_eq!(uy, 0.0);
    }

    #[test]
    fn test_darcy_velocity_y() {
        let d = DarcyFlow::new(0.01, 1e-3, 0.4);
        let (ux, uy) = d.darcy_velocity(0.0, -1.0);
        assert_eq!(ux, 0.0);
        assert!((uy - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_darcy_pore_velocity() {
        // v = u / eps = 0.1 / 0.4 = 0.25
        let d = DarcyFlow::new(0.01, 1e-3, 0.4);
        let (vx, _) = d.pore_velocity(-1.0, 0.0);
        assert!((vx - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_darcy_hydraulic_conductivity() {
        // K_h = K * rho * g / mu = 1e-3 * 1000 * 9.81 / 0.001 = 9810
        let d = DarcyFlow::new(0.001, 1e-3, 0.4);
        let kh = d.hydraulic_conductivity(1000.0, 9.81);
        assert!((kh - 9810.0).abs() < 1e-6);
    }

    #[test]
    fn test_darcy_number() {
        let d = DarcyFlow::new(0.01, 1e-6, 0.4);
        // Da = K/L^2 = 1e-6 / 1^2 = 1e-6
        let da = d.darcy_number(1.0);
        assert!((da - 1e-6).abs() < 1e-20);
    }

    #[test]
    fn test_darcy_zero_viscosity() {
        let d = DarcyFlow::new(0.0, 1e-3, 0.4);
        let (ux, uy) = d.darcy_velocity(-1.0, 0.0);
        assert_eq!(ux, 0.0);
        assert_eq!(uy, 0.0);
    }

    // ── ForcheimerCorrection ───────────────────────────────────────────────

    #[test]
    fn test_forchheimer_ergun_beta() {
        // beta = 0.55 / sqrt(K)
        let beta = ForcheimerCorrection::ergun_beta(1e-4);
        assert!((beta - 0.55 / (1e-4_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_forchheimer_pressure_gradient_low_velocity() {
        // At low Re, pressure gradient ≈ Darcy term
        let f = ForcheimerCorrection::new(1e-3, 0.01, 0.01, 1.0);
        let dp = f.pressure_gradient_1d(0.001);
        let darcy = (0.01 / 1e-3) * 0.001;
        // Inertial term = 0.01 * 1.0 * 0.001 * 0.001 = 1e-8 (negligible)
        assert!((dp - darcy).abs() / darcy < 0.01);
    }

    #[test]
    fn test_forchheimer_inertial_drag_positive() {
        let f = ForcheimerCorrection::new(1e-3, 0.01, 1.0, 1.0);
        let drag = f.inertial_drag(1.0);
        assert!(drag > 0.0);
    }

    #[test]
    fn test_forchheimer_inertial_drag_zero_velocity() {
        let f = ForcheimerCorrection::new(1e-3, 0.01, 1.0, 1.0);
        assert_eq!(f.inertial_drag(0.0), 0.0);
    }

    #[test]
    fn test_forchheimer_effective_permeability_reduces() {
        // Effective K < intrinsic K at non-zero velocity
        let f = ForcheimerCorrection::new(1e-3, 0.01, 1.0, 1.0);
        let k_eff = f.effective_permeability(1.0);
        assert!(k_eff < 1e-3);
    }

    #[test]
    fn test_forchheimer_number_zero_velocity() {
        let f = ForcheimerCorrection::new(1e-3, 0.01, 1.0, 1.0);
        assert_eq!(f.forchheimer_number(0.0), 0.0);
    }

    #[test]
    fn test_forchheimer_ergun_beta_zero_permeability() {
        assert_eq!(ForcheimerCorrection::ergun_beta(0.0), 0.0);
    }

    // ── MultiscalePorous ─────────────────────────────────────────────────

    #[test]
    fn test_multiscale_rev_porosity() {
        let mp = MultiscalePorous::new(vec![0.3, 0.5, 0.7], vec![1e-3; 3], vec![1.5; 3]);
        assert!((mp.rev_porosity() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_multiscale_parallel_perm() {
        let mp = MultiscalePorous::new(vec![0.4; 3], vec![1e-3, 2e-3, 3e-3], vec![1.0; 3]);
        assert!((mp.upscaled_permeability_parallel() - 2e-3).abs() < 1e-18);
    }

    #[test]
    fn test_multiscale_series_perm_less_than_parallel() {
        let mp = MultiscalePorous::new(vec![0.4; 3], vec![1e-3, 2e-3, 3e-3], vec![1.0; 3]);
        let par = mp.upscaled_permeability_parallel();
        let ser = mp.upscaled_permeability_series();
        assert!(ser < par, "ser={ser}, par={par}");
    }

    #[test]
    fn test_multiscale_uniform_tortuosity() {
        let mp = MultiscalePorous::new(vec![0.4; 3], vec![1e-3; 3], vec![2.0; 3]);
        assert!((mp.effective_tortuosity() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_multiscale_effective_diffusivity() {
        // D_eff = D * eps / tau = 1.0 * 0.5 / 2.0 = 0.25
        let mp = MultiscalePorous::new(vec![0.5; 4], vec![1e-3; 4], vec![2.0; 4]);
        let d = mp.effective_diffusivity(1.0);
        assert!((d - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_multiscale_std_dev_uniform() {
        let mp = MultiscalePorous::new(vec![0.4; 5], vec![1e-3; 5], vec![1.0; 5]);
        assert!(mp.porosity_std_dev() < 1e-12);
    }

    #[test]
    fn test_multiscale_std_dev_nonzero() {
        let mp = MultiscalePorous::new(vec![0.2, 0.8], vec![1e-3; 2], vec![1.0; 2]);
        assert!(mp.porosity_std_dev() > 0.0);
    }

    // ── PorousLBM ─────────────────────────────────────────────────────────

    #[test]
    fn test_porous_lbm_initial_density() {
        let geo = PorousGeometry::uniform(4, 4, 0.5, 1e-2, 1e-2, 1.5);
        let lbm = PorousLBM::new(4, 4, 1.0, geo, 0.0, 0.0);
        let rho = lbm.density(0, 0);
        assert!((rho - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_porous_lbm_initial_velocity_zero() {
        let geo = PorousGeometry::uniform(4, 4, 0.5, 1e-2, 1e-2, 1.5);
        let lbm = PorousLBM::new(4, 4, 1.0, geo, 0.0, 0.0);
        let (ux, uy) = lbm.velocity(0, 0);
        assert!(ux.abs() < 1e-10);
        assert!(uy.abs() < 1e-10);
    }

    #[test]
    fn test_porous_lbm_step_runs() {
        let geo = PorousGeometry::uniform(4, 4, 0.5, 1e-2, 1e-2, 1.5);
        let mut lbm = PorousLBM::new(4, 4, 1.0, geo, 1e-4, 0.0);
        lbm.step();
    }

    #[test]
    fn test_porous_lbm_mass_conservation() {
        let geo = PorousGeometry::uniform(6, 6, 0.5, 1e-2, 1e-2, 1.5);
        let mut lbm = PorousLBM::new(6, 6, 1.0, geo, 0.0, 0.0);
        let m0 = lbm.total_mass();
        lbm.step();
        let m1 = lbm.total_mass();
        assert!((m0 - m1).abs() < 1e-8, "mass change: {}", (m0 - m1).abs());
    }

    #[test]
    fn test_porous_lbm_force_generates_velocity() {
        let geo = PorousGeometry::uniform(8, 4, 0.5, 1e-1, 1e-1, 1.0);
        let mut lbm = PorousLBM::new(8, 4, 1.0, geo, 1e-3, 0.0);
        for _ in 0..50 {
            lbm.step();
        }
        // With a body force in x, ux should be positive somewhere
        let (ux, _) = lbm.velocity(4, 2);
        assert!(ux > 0.0, "ux={ux}");
    }

    #[test]
    fn test_porous_lbm_node_count() {
        let geo = PorousGeometry::uniform(5, 3, 0.5, 1e-2, 1e-2, 1.5);
        let lbm = PorousLBM::new(5, 3, 1.0, geo, 0.0, 0.0);
        assert_eq!(lbm.f.len(), 15);
    }

    // ── DrainageImbibition ────────────────────────────────────────────────

    #[test]
    fn test_drainage_capillary_pressure_at_residual() {
        // At S_w just above S_wr, Se → 0+, Pc should be very large
        // S_e = (S_w - S_wr) / (1 - S_wr - S_nr) = 1e-9 / 0.85
        // Pc = Pe * Se^(-1/lambda) = 1000 * (1e-9/0.85)^(-0.5) ≈ very large
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let pc = di.capillary_pressure(0.1 + 1e-9);
        assert!(pc > 1e6, "pc={pc}");
    }

    #[test]
    fn test_drainage_capillary_pressure_fully_saturated() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let pc = di.capillary_pressure(1.0 - 0.05);
        assert_eq!(pc, 0.0);
    }

    #[test]
    fn test_drainage_relative_permeability_wetting_range() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let krw = di.relative_permeability_wetting(0.5);
        assert!(krw >= 0.0 && krw <= 1.0, "krw={krw}");
    }

    #[test]
    fn test_drainage_relative_permeability_nonwetting_range() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let krnw = di.relative_permeability_nonwetting(0.5);
        assert!(krnw >= 0.0 && krnw <= 1.0, "krnw={krnw}");
    }

    #[test]
    fn test_drainage_wetting_krw_increases_with_sw() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let krw1 = di.relative_permeability_wetting(0.4);
        let krw2 = di.relative_permeability_wetting(0.7);
        assert!(krw2 > krw1, "krw1={krw1}, krw2={krw2}");
    }

    #[test]
    fn test_drainage_nonwetting_krnw_decreases_with_sw() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let krnw1 = di.relative_permeability_nonwetting(0.3);
        let krnw2 = di.relative_permeability_nonwetting(0.7);
        assert!(krnw2 < krnw1, "krnw1={krnw1}, krnw2={krnw2}");
    }

    #[test]
    fn test_drainage_effective_saturation_clamp() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let se = di.effective_saturation(2.0); // beyond physical range
        assert!(se <= 1.0, "se={se}");
        let se2 = di.effective_saturation(-1.0);
        assert!(se2 >= 0.0, "se2={se2}");
    }

    #[test]
    fn test_drainage_fractional_flow_unit_viscosity_ratio() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let fw = di.fractional_flow(0.5, 1.0);
        assert!(fw >= 0.0 && fw <= 1.0, "fw={fw}");
    }

    #[test]
    fn test_drainage_curve_length() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let curve = di.drainage_curve(10);
        assert_eq!(curve.len(), 10);
    }

    #[test]
    fn test_drainage_curve_pressure_decreasing() {
        let di = DrainageImbibition::new(1000.0, 2.0, 0.1, 0.05);
        let curve = di.drainage_curve(5);
        // Capillary pressure should decrease as S_w increases
        for w in curve.windows(2) {
            assert!(w[1].1 <= w[0].1, "Pc not monotone: {:?} {:?}", w[0], w[1]);
        }
    }

    // ── Cross-model integration ───────────────────────────────────────────

    #[test]
    fn test_brinkman_darcy_consistency() {
        // Brinkman velocity should match Darcy in the eps → 1 limit
        let b = BrinkmanModel::new(0.01, 1.0, 1e-3);
        let d = DarcyFlow::new(0.01, 1e-3, 1.0);
        let u_b = b.darcy_brinkman_velocity_1d(-10.0);
        let (u_d, _) = d.darcy_velocity(-10.0, 0.0);
        assert!((u_b - u_d).abs() < 1e-12);
    }

    #[test]
    fn test_forchheimer_reduces_to_darcy_at_zero_beta() {
        let f = ForcheimerCorrection::new(1e-3, 0.01, 0.0, 1.0);
        let dp = f.pressure_gradient_1d(0.1);
        let darcy = (0.01 / 1e-3) * 0.1;
        assert!((dp - darcy).abs() < 1e-12);
    }
}
