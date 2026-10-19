// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Three-scale (nano → meso → macro) D1Q3 lattice Boltzmann solver
//! ([`MesoscaleLbm`]).

use super::d1q3::{
    CS2, Q, bgk_relax, hydrodynamic_equilibrium, nonequilibrium, stream_periodic,
    validate_relaxation_rate, velocity, zeroth_moment,
};
use crate::error::Result;

/// Smallest density kept in a meso cell when the nano scale removes mass.
const MIN_DENSITY: f64 = 1e-6;

/// Default BGK relaxation rate of the meso lattice (`τ = 2/3`).
const DEFAULT_OMEGA: f64 = 1.5;

/// Three-scale (nano → meso → macro) lattice Boltzmann coupling.
///
/// The meso scale is a periodic D1Q3 BGK lattice.  Every meso cell carries the
/// three populations `f_i` of the lattice velocities `c = {0, +1, -1}` and
/// relaxes towards the isothermal equilibrium
/// `f_i^eq = w_i ρ (1 + 3 c_i u + 4.5 c_i² u² - 1.5 u²)` with `w = {2/3, 1/6, 1/6}`
/// and `cs² = 1/3`.  Collision and streaming conserve mass `Σ f_i` and
/// momentum `Σ c_i f_i` exactly.  The kinematic viscosity of the lattice is
/// `ν = cs² (τ - 1/2)` with `τ = 1/ω`; small density waves travel at the
/// speed `cs` and every pulse broadens with the variance `2 ν t`.
///
/// - **nano → meso** ([`MesoscaleLbm::nano_to_meso`]): the nano values
///   (particle-number-density excess over the reference lattice density `1`)
///   are averaged over the `nano_per_meso` nano cells of each meso cell and
///   stored in `meso_cells`.  The *change* of that average since the previous
///   transfer is injected into the meso populations as mass that moves with
///   the local fluid velocity (the velocity and the non-equilibrium part of the
///   cell are unchanged).  A removal never drives the density below `1e-6`.
/// - **meso scale** ([`MesoscaleLbm::lbm_step`]): BGK collision and periodic
///   streaming.
/// - **meso → macro** ([`MesoscaleLbm::meso_to_macro`]): the macroscopic
///   fields are the moments of the populations, `ρ = Σ f_i` and
///   `u = Σ c_i f_i / ρ`.
///
/// `macro_density` and `macro_velocity` can also be written directly (e.g. to
/// seed an initial condition); the next method call imposes the written values
/// on the populations by replacing the equilibrium part and keeping the
/// non-equilibrium part.
///
/// A freshly created solver is a fluid at rest with the reference density `1`.
#[derive(Debug, Clone)]
pub struct MesoscaleLbm {
    /// Coarse-grained nano field: the mean of the `nano_per_meso` nano values
    /// of each meso cell, set by [`MesoscaleLbm::nano_to_meso`].
    pub meso_cells: Vec<f64>,
    /// Nano-scale input values (particle-number-density excess over the
    /// reference lattice density `1`).
    pub nano_cells: Vec<f64>,
    /// Number of nano cells per meso cell (coarse-graining ratio).
    pub nano_per_meso: usize,
    /// Macro velocity field `u = Σ c_i f_i / ρ` (one value per meso cell).
    pub macro_velocity: Vec<f64>,
    /// Macro density field `ρ = Σ f_i` (one value per meso cell).
    pub macro_density: Vec<f64>,
    /// BGK relaxation rate `ω = 1/τ` of the meso lattice (`0 < ω < 2`).
    pub omega: f64,
    /// D1Q3 populations of the meso cells.
    pub(super) meso_pops: Vec<[f64; Q]>,
    /// Coarse-grained nano field that has already been transferred to the
    /// meso populations.
    nano_transferred: Vec<f64>,
}

impl MesoscaleLbm {
    /// Create a new three-scale LBM.
    ///
    /// `n_meso` is the number of meso cells; `nano_per_meso` is how many nano
    /// cells map to one meso cell.  The fluid starts at rest with `ρ = 1`.
    pub fn new(n_meso: usize, nano_per_meso: usize) -> Self {
        let rest = hydrodynamic_equilibrium(1.0, 0.0);
        Self {
            meso_cells: vec![0.0; n_meso],
            nano_cells: vec![0.0; n_meso * nano_per_meso],
            nano_per_meso,
            macro_velocity: vec![velocity(&rest); n_meso],
            macro_density: vec![zeroth_moment(&rest); n_meso],
            omega: DEFAULT_OMEGA,
            meso_pops: vec![rest; n_meso],
            nano_transferred: vec![0.0; n_meso],
        }
    }

    /// Number of meso cells.
    pub fn n_meso(&self) -> usize {
        self.meso_cells.len()
    }

    /// Set the BGK relaxation rate `ω`.
    ///
    /// # Errors
    /// Fails unless `0 < ω < 2` (`τ > 1/2`).
    pub fn set_relaxation_rate(&mut self, omega: f64) -> Result<()> {
        validate_relaxation_rate(omega, "omega")?;
        self.omega = omega;
        Ok(())
    }

    /// Relaxation time `τ = 1/ω` of the meso lattice.
    pub fn relaxation_time(&self) -> f64 {
        1.0 / self.omega
    }

    /// Kinematic viscosity `ν = cs² (τ - 1/2)` of the meso lattice (lattice units).
    pub fn kinematic_viscosity(&self) -> f64 {
        CS2 * (self.relaxation_time() - 0.5)
    }

    /// Total mass `Σ ρ` of the macro density field.
    pub fn total_mass(&self) -> f64 {
        self.macro_density.iter().sum()
    }

    /// Total momentum `Σ ρ u` of the macro fields.
    pub fn total_momentum(&self) -> f64 {
        self.macro_density
            .iter()
            .zip(self.macro_velocity.iter())
            .map(|(&rho, &u)| rho * u)
            .sum()
    }

    /// Impose externally written macro fields on the populations: wherever the
    /// written `(ρ, u)` differs from the moments of a cell, replace the
    /// equilibrium part of its populations and keep the non-equilibrium part.
    fn absorb_macro_writes(&mut self) {
        for ((g, &rho_new), &u_new) in self
            .meso_pops
            .iter_mut()
            .zip(self.macro_density.iter())
            .zip(self.macro_velocity.iter())
        {
            let rho = zeroth_moment(g);
            let u = velocity(g);
            if rho_new != rho || u_new != u {
                let neq = nonequilibrium(g, &hydrodynamic_equilibrium(rho, u));
                let feq = hydrodynamic_equilibrium(rho_new, u_new);
                for q in 0..Q {
                    g[q] = feq[q] + neq[q];
                }
            }
        }
    }

    /// Refresh the macro fields from the populations.
    fn refresh_macro_fields(&mut self) {
        for ((rho, u), g) in self
            .macro_density
            .iter_mut()
            .zip(self.macro_velocity.iter_mut())
            .zip(self.meso_pops.iter())
        {
            *rho = zeroth_moment(g);
            *u = velocity(g);
        }
    }

    /// Nano-to-meso coupling: coarse-grain the nano values and transfer them to
    /// the meso populations.
    ///
    /// The `nano_per_meso` nano values of every meso cell are averaged into
    /// `meso_cells`.  The change of this average since the previous call is
    /// added to the cell as mass moving with the local fluid velocity, so the
    /// meso mass follows the nano mass and calling the method twice with
    /// unchanged nano data transfers nothing the second time.  Removing more
    /// mass than is available leaves the density at `1e-6`; the missing amount
    /// stays pending.
    pub fn nano_to_meso(&mut self) {
        self.absorb_macro_writes();
        let r = self.nano_per_meso;
        for mi in 0..self.n_meso() {
            let mut sum = 0.0;
            for k in 0..r {
                if let Some(&value) = self.nano_cells.get(mi * r + k) {
                    sum += value;
                }
            }
            let mean = if r == 0 { 0.0 } else { sum / r as f64 };
            self.meso_cells[mi] = mean;
            let pending = mean - self.nano_transferred[mi];
            if pending == 0.0 {
                continue;
            }
            let g = &mut self.meso_pops[mi];
            let rho = zeroth_moment(g);
            let floor = (MIN_DENSITY - rho).min(0.0);
            let applied = pending.max(floor);
            let shape = hydrodynamic_equilibrium(applied, velocity(g));
            for q in 0..Q {
                g[q] += shape[q];
            }
            self.nano_transferred[mi] += applied;
        }
        self.refresh_macro_fields();
    }

    /// One meso-scale time step: BGK collision towards the isothermal
    /// equilibrium of the local `(ρ, u)` followed by periodic streaming.
    ///
    /// Mass and momentum are conserved exactly.
    pub fn lbm_step(&mut self) {
        self.absorb_macro_writes();
        let omega = self.omega;
        for g in self.meso_pops.iter_mut() {
            let feq = hydrodynamic_equilibrium(zeroth_moment(g), velocity(g));
            bgk_relax(g, &feq, omega);
        }
        stream_periodic(&mut self.meso_pops);
        self.refresh_macro_fields();
    }

    /// Meso-to-macro projection: the macro density and velocity are the
    /// moments of the meso populations, `ρ = Σ f_i` and `u = Σ c_i f_i / ρ`
    /// (with `c_i ∈ {0, +1, -1}`).
    pub fn meso_to_macro(&mut self) {
        self.absorb_macro_writes();
        self.refresh_macro_fields();
    }

    /// Perform one heterogeneous time step:
    /// 1. [`MesoscaleLbm::nano_to_meso`] coarse-graining and mass transfer,
    /// 2. [`MesoscaleLbm::lbm_step`] BGK collision and streaming on the meso scale,
    /// 3. [`MesoscaleLbm::meso_to_macro`] moment projection.
    pub fn heterogeneous_step(&mut self) {
        self.nano_to_meso();
        self.lbm_step();
        self.meso_to_macro();
    }
}
