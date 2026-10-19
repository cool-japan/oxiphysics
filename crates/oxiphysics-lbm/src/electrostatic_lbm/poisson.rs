// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! D2Q5 lattice-Boltzmann Poisson solver with Dirichlet cells.

use super::{CX5, CY5, W5, constants, lattice_diffusivity, stream_periodic};
use crate::error::{Error, Result};

/// Relative tolerance below which the net charge of a periodic problem is
/// treated as zero (compatibility condition of the periodic Poisson problem).
const NET_CHARGE_REL_TOL: f64 = 1.0e-9;

/// LBM solver for the 2-D Poisson equation ∇²φ = −ρ/ε.
///
/// The solver evolves the pseudo-transient diffusion equation
/// `∂φ/∂t = D_lu·∇²_lu φ + S` on a D2Q5 BGK lattice whose diffusivity is
/// `D_lu = (1/ω − ½)/3`.  The per-step source `S = D_lu·Δx²·ρ/ε` is chosen so
/// that the steady state satisfies the physical Poisson equation on a grid of
/// spacing [`PoissonLbm::dx`].
///
/// Streaming is periodic.  Cells pinned with [`PoissonLbm::set_phi`] are
/// Dirichlet nodes: their potential is held at the prescribed value on every
/// step.  A purely periodic problem (no Dirichlet cell) only has a solution
/// when the net charge vanishes; [`PoissonLbm::solve`] reports an error
/// otherwise.
pub struct PoissonLbm {
    /// Grid width (number of cells in x-direction).
    pub nx: usize,
    /// Grid height (number of cells in y-direction).
    pub ny: usize,
    /// Distribution functions f\[y * nx * 5 + x * 5 + q\].
    pub f: Vec<f64>,
    /// Charge density ρ (C/m³) at each cell.
    pub rho_charge: Vec<f64>,
    /// Electric potential φ (V) at each cell.
    pub phi: Vec<f64>,
    /// Relative permittivity (dimensionless).
    pub epsilon_r: f64,
    /// BGK relaxation parameter ω ∈ (0, 2).
    pub omega: f64,
    /// Physical cell spacing Δx (m).
    pub dx: f64,
    /// Dirichlet constraints: `Some(φ)` pins the potential of that cell.
    pub fixed: Vec<Option<f64>>,
    /// Streaming buffer reused between steps.
    scratch: Vec<f64>,
}

impl PoissonLbm {
    /// Create a new Poisson LBM solver on an nx × ny grid with unit spacing
    /// (Δx = 1 m; use [`PoissonLbm::with_spacing`] for another spacing).
    ///
    /// # Arguments
    /// * `nx` - grid width
    /// * `ny` - grid height
    /// * `epsilon_r` - relative permittivity of the medium
    /// * `omega` - BGK relaxation rate (1.0 gives the fastest convergence)
    pub fn new(nx: usize, ny: usize, epsilon_r: f64, omega: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            f: vec![0.0; n * 5],
            rho_charge: vec![0.0; n],
            phi: vec![0.0; n],
            epsilon_r,
            omega,
            dx: 1.0,
            fixed: vec![None; n],
            scratch: vec![0.0; n * 5],
        }
    }

    /// Set the physical cell spacing Δx (m) and return the solver.
    pub fn with_spacing(mut self, dx: f64) -> Self {
        self.dx = dx;
        self
    }

    /// Linear index for cell (x, y).
    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Lattice diffusivity `D_lu = (1/ω − ½)/3` of the pseudo-transient
    /// relaxation.
    pub fn lattice_diffusivity(&self) -> f64 {
        lattice_diffusivity(self.omega)
    }

    /// Whether at least one Dirichlet cell has been set.
    pub fn has_dirichlet(&self) -> bool {
        self.fixed.iter().any(Option::is_some)
    }

    /// Initialize distributions from the current φ field (Dirichlet cells take
    /// their prescribed value).
    pub fn init_from_phi(&mut self) {
        for cell in 0..self.nx * self.ny {
            if let Some(val) = self.fixed[cell] {
                self.phi[cell] = val;
            }
            let p = self.phi[cell];
            for (q, w5_q) in W5.iter().enumerate() {
                self.f[cell * 5 + q] = w5_q * p;
            }
        }
    }

    /// Set charge density (C/m³) at cell (x, y).
    pub fn set_charge(&mut self, x: usize, y: usize, rho: f64) {
        let i = self.idx(x, y);
        self.rho_charge[i] = rho;
    }

    /// Pin the potential of cell (x, y) to `val` (Dirichlet condition).
    ///
    /// The value is held on every subsequent [`PoissonLbm::step`] until the
    /// cell is released with [`PoissonLbm::release_phi`].
    pub fn set_phi(&mut self, x: usize, y: usize, val: f64) {
        let cell = self.idx(x, y);
        self.phi[cell] = val;
        self.fixed[cell] = Some(val);
        for (q, w5_q) in W5.iter().enumerate() {
            self.f[cell * 5 + q] = w5_q * val;
        }
    }

    /// Remove the Dirichlet condition of cell (x, y); its potential evolves
    /// freely afterwards.
    pub fn release_phi(&mut self, x: usize, y: usize) {
        let cell = self.idx(x, y);
        self.fixed[cell] = None;
    }

    /// Lattice source `S = D_lu·Δx²/ε` per unit charge density.
    fn source_scale(&self) -> f64 {
        let eps = self.epsilon_r * constants::EPSILON_0;
        self.lattice_diffusivity() * self.dx * self.dx / eps
    }

    /// Periodic neighbour of cell (x, y) along lattice direction `q`.
    #[inline]
    fn neighbour(&self, x: usize, y: usize, q: usize) -> (usize, usize) {
        let xn = (x as i64 + i64::from(CX5[q])).rem_euclid(self.nx as i64) as usize;
        let yn = (y as i64 + i64::from(CY5[q])).rem_euclid(self.ny as i64) as usize;
        (xn, yn)
    }

    /// Lattice derivative of φ from a Dirichlet cell (x, y) holding `val`
    /// towards its free neighbour along direction `q`: second-order one-sided
    /// when the next cell is free as well, first-order otherwise.  `None`
    /// when the neighbour is itself a Dirichlet cell.
    fn directional_derivative(&self, x: usize, y: usize, q: usize, val: f64) -> Option<f64> {
        let (x1, y1) = self.neighbour(x, y, q);
        let i1 = self.idx(x1, y1);
        if self.fixed[i1].is_some() {
            return None;
        }
        let (x2, y2) = self.neighbour(x1, y1, q);
        let i2 = self.idx(x2, y2);
        let own = self.idx(x, y);
        if self.fixed[i2].is_none() && i2 != own && i2 != i1 {
            Some(0.5 * (-3.0 * val + 4.0 * self.phi[i1] - self.phi[i2]))
        } else {
            Some(self.phi[i1] - val)
        }
    }

    /// Perform one BGK collision + streaming step.
    ///
    /// Free cells relax towards `w_q·φ` and receive the source `w_q·S`.  A
    /// Dirichlet cell emits the post-collision state of a cell whose potential
    /// is exactly the prescribed value: `w_q·(φ_D + S)` plus the first-order
    /// non-equilibrium part `−(1/ω − 1)·w_q·c_q·∇φ` towards each free
    /// neighbour (non-equilibrium extrapolation), which keeps the boundary
    /// second-order accurate for every ω.
    ///
    /// After the step [`PoissonLbm::phi`] holds the zeroth moment of the
    /// streamed distributions (Dirichlet cells hold their prescribed value).
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let omega = self.omega;
        let scale = self.source_scale();
        let neq_factor = 1.0 / omega - 1.0;

        // --- macroscopic potential of the current state ---
        for i in 0..n {
            self.phi[i] = match self.fixed[i] {
                Some(val) => val,
                None => self.f[i * 5..i * 5 + 5].iter().sum(),
            };
        }

        // --- collision (with source) ---
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                let src = scale * self.rho_charge[i];
                let base = i * 5;
                if let Some(val) = self.fixed[i] {
                    let mut post = [0.0_f64; 5];
                    for (q, w5_q) in W5.iter().enumerate() {
                        post[q] = w5_q * (val + src);
                        if q > 0
                            && neq_factor != 0.0
                            && let Some(grad) = self.directional_derivative(x, y, q, val)
                        {
                            post[q] -= neq_factor * w5_q * grad;
                        }
                    }
                    self.f[base..base + 5].copy_from_slice(&post);
                } else {
                    let phi_loc = self.phi[i];
                    for (q, w5_q) in W5.iter().enumerate() {
                        let feq = w5_q * phi_loc;
                        self.f[base + q] += omega * (feq - self.f[base + q]) + w5_q * src;
                    }
                }
            }
        }

        // --- streaming (periodic) ---
        if self.scratch.len() != self.f.len() {
            self.scratch = vec![0.0; self.f.len()];
        }
        stream_periodic(&self.f, nx, ny, &mut self.scratch);
        std::mem::swap(&mut self.f, &mut self.scratch);

        // --- macroscopic potential of the streamed state ---
        for i in 0..n {
            self.phi[i] = match self.fixed[i] {
                Some(val) => val,
                None => self.f[i * 5..i * 5 + 5].iter().sum(),
            };
        }
    }

    /// Run until convergence or `max_iter` steps.
    ///
    /// Returns the number of iterations performed and the final L∞ change of
    /// φ between the last two iterations.  Reaching `max_iter` without meeting
    /// `tol` is not an error: the caller inspects the returned residual.
    ///
    /// # Errors
    /// * [`Error::ParameterOutOfBounds`] if ω ∉ (0, 2) or Δx/ε_r is not
    ///   positive and finite.
    /// * [`Error::ConvergenceFailure`] if there is no Dirichlet cell and the
    ///   net charge is non-zero: the periodic Poisson problem then has no
    ///   solution and the mean potential grows without bound.  The reported
    ///   residual is that per-iteration growth of the mean potential.
    pub fn solve(&mut self, max_iter: usize, tol: f64) -> Result<(usize, f64)> {
        if !(self.omega > 0.0 && self.omega < 2.0) {
            return Err(Error::ParameterOutOfBounds {
                name: "omega".to_string(),
                value: self.omega,
                min: 0.0,
                max: 2.0,
            });
        }
        if !(self.dx > 0.0 && self.dx.is_finite()) {
            return Err(Error::ParameterOutOfBounds {
                name: "dx".to_string(),
                value: self.dx,
                min: 0.0,
                max: f64::MAX,
            });
        }
        if !(self.epsilon_r > 0.0 && self.epsilon_r.is_finite()) {
            return Err(Error::ParameterOutOfBounds {
                name: "epsilon_r".to_string(),
                value: self.epsilon_r,
                min: 0.0,
                max: f64::MAX,
            });
        }
        if !self.has_dirichlet() {
            let net: f64 = self.rho_charge.iter().sum();
            let total: f64 = self.rho_charge.iter().map(|r| r.abs()).sum();
            if net.abs() > NET_CHARGE_REL_TOL * total {
                let cells = (self.nx * self.ny).max(1) as f64;
                return Err(Error::ConvergenceFailure {
                    residual: self.source_scale() * net.abs() / cells,
                    iterations: 0,
                    tolerance: tol,
                });
            }
        }
        self.init_from_phi();
        let mut residual = f64::MAX;
        let mut iter = 0;
        let mut phi_old = self.phi.clone();
        while iter < max_iter && residual > tol {
            phi_old.copy_from_slice(&self.phi);
            self.step();
            residual = self
                .phi
                .iter()
                .zip(phi_old.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            iter += 1;
        }
        Ok((iter, residual))
    }

    /// Return the electric potential at cell (x, y).
    pub fn get_phi(&self, x: usize, y: usize) -> f64 {
        self.phi[self.idx(x, y)]
    }
}
