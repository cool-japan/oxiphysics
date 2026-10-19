// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Passive-scalar advection-diffusion LBM on the D2Q9 lattice.
//!
//! [`PassiveScalarLbm`] evolves one scalar `φ` (a species mass fraction or the temperature)
//! with a set of populations `g_q` that follow the flow velocity `u` of the host LBM solver:
//!
//! ```text
//! geq_q = w_q φ (1 + c_q·u / cs²)                      (first-order equilibrium)
//! g_q  ← g_q − ω (g_q − geq_q)                         (BGK collision)
//! g_q(x + c_q) ← g_q(x)                                (periodic streaming)
//! φ = Σ_q g_q
//! ```
//!
//! The Chapman-Enskog expansion of this scheme recovers
//! `∂φ/∂t + ∇·(φ u) = D ∇²φ` with the diffusivity `D = cs² (1/ω − 1/2)` in lattice units
//! (cells² per step), so `ω = 1 / (D/cs² + 1/2)`. Because the first moment of the
//! populations relaxes towards `φ u`, a uniform flow moves a blob of scalar with exactly the flow
//! speed, and the total amount `Σ φ` is conserved to round-off.
//!
//! Source terms (chemical reaction, heat release) and external edits of the macroscopic field
//! are fed into the populations by [`PassiveScalarLbm::sync_to_field`], which distributes the
//! change `Δφ = φ_field − Σ g_q` over the lattice directions with the weights `w_q`, i.e. the
//! source enters as `w_q S Δt`.

use super::{C9, CS2, W9, scalar_equilibrium_d2q9};

/// Number of D2Q9 lattice directions.
const NQ: usize = 9;

/// D2Q9 advection-diffusion LBM for one passive scalar on a periodic `nx × ny` grid.
///
/// The cell index is `iy * nx + ix`, like the macroscopic fields of the combustion solver.
#[derive(Clone, Debug)]
pub struct PassiveScalarLbm {
    /// Number of cells in x.
    nx: usize,
    /// Number of cells in y.
    ny: usize,
    /// BGK relaxation rate `ω = 1 / (D/cs² + 1/2)`.
    omega: f64,
    /// Scalar populations `g_q` per cell.
    g: Vec<[f64; NQ]>,
    /// Streaming buffer (swapped with `g` every step).
    scratch: Vec<[f64; NQ]>,
}

/// BGK relaxation rate that realises the lattice diffusivity `diffusivity` (cells² per step).
///
/// A negative diffusivity is treated as zero, which gives the marginal rate `ω = 2`.
pub fn scalar_relaxation_rate(diffusivity: f64) -> f64 {
    1.0 / (diffusivity.max(0.0) / CS2 + 0.5)
}

impl PassiveScalarLbm {
    /// Create a solver on an `nx × ny` grid holding the uniform scalar value `value`, at rest
    /// equilibrium, with the lattice diffusivity `diffusivity` (cells² per step).
    pub fn uniform(nx: usize, ny: usize, diffusivity: f64, value: f64) -> Self {
        let rest: [f64; NQ] = std::array::from_fn(|q| W9[q] * value);
        let n = nx * ny;
        Self {
            nx,
            ny,
            omega: scalar_relaxation_rate(diffusivity),
            g: vec![rest; n],
            scratch: vec![rest; n],
        }
    }

    /// Create a solver whose scalar field is `phi` (length `nx * ny`), at rest equilibrium.
    ///
    /// # Errors
    ///
    /// [`Error::ComponentCountMismatch`](crate::error::Error::ComponentCountMismatch) if `phi`
    /// does not hold `nx * ny` values.
    pub fn from_field(
        nx: usize,
        ny: usize,
        diffusivity: f64,
        phi: &[f64],
    ) -> crate::error::Result<Self> {
        let n = nx * ny;
        if phi.len() != n {
            return Err(crate::error::Error::ComponentCountMismatch {
                expected: n,
                found: phi.len(),
                field: "scalar field".to_string(),
            });
        }
        let mut solver = Self::uniform(nx, ny, diffusivity, 0.0);
        solver.sync_to_field(phi);
        Ok(solver)
    }

    /// Number of cells in x.
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Number of cells in y.
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// BGK relaxation rate `ω` of the scalar populations.
    pub fn omega(&self) -> f64 {
        self.omega
    }

    /// Lattice diffusivity `D = cs² (1/ω − 1/2)` (cells² per step).
    pub fn diffusivity(&self) -> f64 {
        CS2 * (1.0 / self.omega - 0.5)
    }

    /// Change the lattice diffusivity (cells² per step); the populations are kept.
    pub fn set_diffusivity(&mut self, diffusivity: f64) {
        self.omega = scalar_relaxation_rate(diffusivity);
    }

    /// Scalar populations `g_q` of every cell (`iy * nx + ix` layout).
    pub fn populations(&self) -> &[[f64; NQ]] {
        &self.g
    }

    /// Scalar value `Σ_q g_q` of cell `cell`, or `None` outside the grid.
    pub fn value(&self, cell: usize) -> Option<f64> {
        self.g.get(cell).map(|g| g.iter().sum())
    }

    /// Total amount of scalar `Σ φ` over the grid.
    pub fn total(&self) -> f64 {
        self.g.iter().map(|g| g.iter().sum::<f64>()).sum()
    }

    /// Make the populations agree with the macroscopic field `phi`.
    ///
    /// For every cell the difference `Δφ = phi − Σ_q g_q` is added to the populations as
    /// `g_q += w_q Δφ`. This keeps the stored non-equilibrium part (and with it the diffusive
    /// flux), imports external edits of the field and realises a source term `S Δt` as
    /// `w_q S Δt`. Only the common prefix of the grid and `phi` is synchronised.
    pub fn sync_to_field(&mut self, phi: &[f64]) {
        for (g, &value) in self.g.iter_mut().zip(phi.iter()) {
            let delta = value - g.iter().sum::<f64>();
            for (g_q, &w_q) in g.iter_mut().zip(W9.iter()) {
                *g_q += w_q * delta;
            }
        }
    }

    /// Write the macroscopic scalar `φ = Σ_q g_q` into `phi` (common prefix of the grid and `phi`).
    pub fn read_field(&self, phi: &mut [f64]) {
        for (value, g) in phi.iter_mut().zip(self.g.iter()) {
            *value = g.iter().sum();
        }
    }

    /// Advance the scalar by one lattice step in the velocity field `(ux, uy)`: BGK collision
    /// towards the advected equilibrium followed by periodic pull streaming.
    ///
    /// A velocity slice that is shorter than the grid is padded with zero velocity.
    pub fn step(&mut self, ux: &[f64], uy: &[f64]) {
        let omega = self.omega;
        for (cell, g) in self.g.iter_mut().enumerate() {
            let phi: f64 = g.iter().sum();
            let u_x = ux.get(cell).copied().unwrap_or(0.0);
            let u_y = uy.get(cell).copied().unwrap_or(0.0);
            let geq = scalar_equilibrium_d2q9(phi, u_x, u_y);
            for (g_q, &geq_q) in g.iter_mut().zip(geq.iter()) {
                *g_q -= omega * (*g_q - geq_q);
            }
        }

        let nx = self.nx as i64;
        let ny = self.ny as i64;
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let dst = iy * self.nx + ix;
                for (q, &(cx, cy)) in C9.iter().enumerate() {
                    let src_x = (ix as i64 - i64::from(cx)).rem_euclid(nx) as usize;
                    let src_y = (iy as i64 - i64::from(cy)).rem_euclid(ny) as usize;
                    self.scratch[dst][q] = self.g[src_y * self.nx + src_x][q];
                }
            }
        }
        std::mem::swap(&mut self.g, &mut self.scratch);
    }
}
