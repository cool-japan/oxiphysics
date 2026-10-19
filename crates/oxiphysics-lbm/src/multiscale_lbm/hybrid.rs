// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hybrid LBM / Navier-Stokes domain decomposition ([`HybridLbmNavier`]).

use super::d1q3::{
    CS2, Q, advection_equilibrium, bgk_relax, nonequilibrium, stream_open, zeroth_moment,
};

/// Hybrid LBM / Navier-Stokes domain decomposition solver.
///
/// The solver treats the unidirectional shear flow `u(x)` of a fluid of
/// uniform density: the streamwise velocity diffuses across the cells,
/// `∂u/∂t = ν ∂²u/∂x²`.  The domain is split into an LBM sub-domain and a
/// Navier-Stokes (NS) sub-domain with an overlap of `overlap_width` cells on
/// each side of the interface.  Both sub-domains solve the same equation:
///
/// - the **NS side** uses an explicit finite-difference Laplacian
///   (`ν_ns dt/dx²` per step);
/// - the **LBM side** is a D1Q3 BGK lattice (`w = {2/3, 1/6, 1/6}`,
///   `c = {0, +1, -1}`, `cs² = 1/3`) whose populations carry the momentum
///   density `q = ρ u` of the shear flow and relax towards `w_i q`.  Its
///   relaxation time follows from the time step, `τ = 1/2 + ν_lbm dt/(cs² dx²)`,
///   so the lattice diffusivity `cs² (τ - 1/2) dx²/dt` equals `ν_lbm`.  A
///   spike therefore splits into `(1/6, 2/3, 1/6)` of its value in one step
///   whatever the viscosity, and the total momentum is conserved exactly.
///
/// The sub-domains exchange data in a Schwarz fashion: before a step the cells
/// outside the LBM region receive the NS velocity (Dirichlet data of the LBM
/// side, imposed with non-equilibrium extrapolation) and the cells outside the
/// NS region receive the LBM velocity; after the step the velocities of the
/// overlap are averaged ([`HybridLbmNavier::schwarz_coupling`]).  The first and
/// the last cell of the domain keep their values.
///
/// `lbm_velocity` is the macroscopic view of the LBM populations and is
/// refreshed by [`HybridLbmNavier::step`].  It may be written directly: the
/// next step imposes the written values on the populations and keeps their
/// non-equilibrium part.
#[derive(Debug, Clone)]
pub struct HybridLbmNavier {
    /// Boolean mask: `true` where LBM is active.
    pub lbm_region: Vec<bool>,
    /// Boolean mask: `true` where Navier-Stokes is active.
    pub ns_region: Vec<bool>,
    /// Width (number of cells) of the overlap / coupling region.
    pub overlap_width: usize,
    /// LBM velocity field (one value per cell, x-component representative).
    pub lbm_velocity: Vec<f64>,
    /// NS velocity field (one value per cell, x-component representative).
    pub ns_velocity: Vec<f64>,
    /// LBM density field (uniform reference density of the shear flow).
    pub lbm_density: Vec<f64>,
    /// NS density field.
    pub ns_density: Vec<f64>,
    /// LBM kinematic viscosity ν (m²/s).
    pub lbm_viscosity: f64,
    /// NS kinematic viscosity ν (m²/s).
    pub ns_viscosity: f64,
    /// Uniform cell size `dx` (m) of both sub-domains.
    pub cell_size: f64,
    /// D1Q3 populations of the momentum density `ρ u` on the LBM side.
    pub(super) lbm_pops: Vec<[f64; Q]>,
}

/// Velocity `Σ g_i / ρ` carried by the populations of a node of density `rho`.
fn population_velocity(g: &[f64; Q], rho: f64) -> f64 {
    if rho > 0.0 {
        zeroth_moment(g) / rho
    } else {
        0.0
    }
}

impl HybridLbmNavier {
    /// Create a new hybrid solver with `n` cells.
    ///
    /// The first `n/2` cells are assigned to LBM, the second half to NS.
    /// The overlap spans `overlap_width` cells on each side of the interface.
    /// The cell size is `1`; `viscosity` is the kinematic viscosity of both
    /// sub-domains.
    pub fn new(n: usize, overlap_width: usize, viscosity: f64) -> Self {
        let half = n / 2;
        let mut lbm_region = vec![false; n];
        let mut ns_region = vec![false; n];
        for i in 0..n {
            if i < half + overlap_width {
                lbm_region[i] = true;
            }
            if i + overlap_width >= half {
                ns_region[i] = true;
            }
        }
        Self {
            lbm_region,
            ns_region,
            overlap_width,
            lbm_velocity: vec![0.0; n],
            ns_velocity: vec![0.0; n],
            lbm_density: vec![1.0; n],
            ns_density: vec![1.0; n],
            lbm_viscosity: viscosity,
            ns_viscosity: viscosity,
            cell_size: 1.0,
            lbm_pops: vec![[0.0; Q]; n],
        }
    }

    /// Number of cells in the domain.
    pub fn n_cells(&self) -> usize {
        self.lbm_region.len()
    }

    /// LBM relaxation time `τ = 1/2 + ν_lbm dt/(cs² dx²)` for the time step `dt`.
    pub fn lbm_relaxation_time(&self, dt: f64) -> f64 {
        0.5 + self.lbm_viscosity * dt / (CS2 * self.cell_size * self.cell_size)
    }

    /// Impose externally written LBM velocities on the populations: wherever
    /// `lbm_velocity[i]` differs from the velocity carried by the populations
    /// shift them by the equilibrium shape `w_i Δq`, which keeps their
    /// non-equilibrium part.
    fn absorb_lbm_writes(&mut self) {
        for ((g, &u_written), &rho) in self
            .lbm_pops
            .iter_mut()
            .zip(self.lbm_velocity.iter())
            .zip(self.lbm_density.iter())
        {
            if u_written != population_velocity(g, rho) {
                let shift = advection_equilibrium(rho * u_written - zeroth_moment(g), 0.0);
                for (gi, s) in g.iter_mut().zip(shift.iter()) {
                    *gi += s;
                }
            }
        }
    }

    /// Refresh the LBM velocity field from the populations.
    fn refresh_lbm_velocity(&mut self) {
        for ((u, g), &rho) in self
            .lbm_velocity
            .iter_mut()
            .zip(self.lbm_pops.iter())
            .zip(self.lbm_density.iter())
        {
            *u = population_velocity(g, rho);
        }
    }

    /// Schwarz data exchange at the artificial interfaces: cells outside the
    /// LBM region take the NS velocity, cells outside the NS region take the
    /// LBM velocity.
    fn exchange_interface_data(&mut self) {
        for i in 0..self.n_cells() {
            if !self.lbm_region[i] {
                self.lbm_velocity[i] = self.ns_velocity[i];
            }
            if !self.ns_region[i] {
                self.ns_velocity[i] = self.lbm_velocity[i];
            }
        }
    }

    /// Advance the LBM sub-domain by one BGK collision / open streaming step.
    ///
    /// Cells outside the LBM region and the two end cells are Dirichlet cells:
    /// their momentum density is held and their populations are rebuilt after
    /// streaming as `w_i q + f_i^neq` with the non-equilibrium part
    /// extrapolated from the nearest evolving neighbour.
    fn advance_lbm(&mut self, dt: f64) {
        let n = self.n_cells();
        if n < 3 {
            return;
        }
        let omega = 1.0 / self.lbm_relaxation_time(dt);
        let is_fixed: Vec<bool> = (0..n)
            .map(|i| i == 0 || i == n - 1 || !self.lbm_region[i])
            .collect();
        let held: Vec<f64> = self.lbm_pops.iter().map(zeroth_moment).collect();
        for g in self.lbm_pops.iter_mut() {
            let geq = advection_equilibrium(zeroth_moment(g), 0.0);
            bgk_relax(g, &geq, omega);
        }
        stream_open(&mut self.lbm_pops);
        for (i, &fixed) in is_fixed.iter().enumerate() {
            if !fixed {
                continue;
            }
            let neighbour = [i + 1, i.wrapping_sub(1)]
                .into_iter()
                .find(|&j| j < n && !is_fixed[j]);
            let neq = match neighbour {
                Some(j) => {
                    let g = &self.lbm_pops[j];
                    nonequilibrium(g, &advection_equilibrium(zeroth_moment(g), 0.0))
                }
                None => [0.0; Q],
            };
            let eq = advection_equilibrium(held[i], 0.0);
            for q in 0..Q {
                self.lbm_pops[i][q] = eq[q] + neq[q];
            }
        }
        self.refresh_lbm_velocity();
    }

    /// Advance the NS sub-domain by one explicit finite-difference step.
    fn advance_ns(&mut self, dt: f64) {
        let n = self.n_cells();
        if n < 3 {
            return;
        }
        let coefficient = self.ns_viscosity * dt / (self.cell_size * self.cell_size);
        let mut ns_new = self.ns_velocity.clone();
        for (ns_out, i) in ns_new[1..n - 1].iter_mut().zip(1..n - 1) {
            if self.ns_region[i] {
                let lap =
                    self.ns_velocity[i + 1] - 2.0 * self.ns_velocity[i] + self.ns_velocity[i - 1];
                *ns_out += coefficient * lap;
            }
        }
        self.ns_velocity = ns_new;
    }

    /// Perform one time step of size `dt`.
    ///
    /// The sub-domains first exchange their interface data, then the LBM
    /// sub-domain performs one D1Q3 BGK collision and streaming step with the
    /// relaxation time [`HybridLbmNavier::lbm_relaxation_time`], the NS
    /// sub-domain performs one explicit finite-difference step, and finally
    /// the overlap is coupled via [`HybridLbmNavier::schwarz_coupling`].
    pub fn step(&mut self, dt: f64) {
        self.exchange_interface_data();
        self.absorb_lbm_writes();
        self.advance_lbm(dt);
        self.advance_ns(dt);
        self.schwarz_coupling();
    }

    /// Schwarz alternating coupling: enforce continuity of velocity across the
    /// overlap region by averaging LBM and NS solutions.
    pub fn schwarz_coupling(&mut self) {
        let n = self.n_cells();
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let avg = 0.5 * (self.lbm_velocity[i] + self.ns_velocity[i]);
                self.lbm_velocity[i] = avg;
                self.ns_velocity[i] = avg;
            }
        }
    }

    /// Flux matching: enforce continuity of momentum flux at the interface.
    ///
    /// Returns the interface flux (ρ·u) at the first overlap cell.
    pub fn flux_matching(&self) -> f64 {
        let n = self.n_cells();
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let flux_lbm = self.lbm_density[i] * self.lbm_velocity[i];
                let flux_ns = self.ns_density[i] * self.ns_velocity[i];
                return 0.5 * (flux_lbm + flux_ns);
            }
        }
        0.0
    }

    /// Velocity mismatch L2-norm across the overlap region.
    ///
    /// Returns `√( Σ (u_lbm - u_ns)² / N_overlap )`.
    pub fn velocity_mismatch(&self) -> f64 {
        let n = self.n_cells();
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let d = self.lbm_velocity[i] - self.ns_velocity[i];
                sum += d * d;
                count += 1;
            }
        }
        if count == 0 {
            return 0.0;
        }
        (sum / count as f64).sqrt()
    }
}
