// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multi-level D1Q3 advection-diffusion lattice Boltzmann solver with
//! subcycled refinement levels ([`AdaptiveLbm`]).

use super::clamp;
use super::d1q3::{
    CS2, Q, advection_equilibrium, bgk_relax, nonequilibrium, stream_periodic,
    validate_relaxation_rate, zeroth_moment,
};
use crate::error::{Error, Result};
use crate::multiscale::{coarse_to_fine_scale, fine_relaxation_time, fine_to_coarse_scale};

/// Largest advection velocity (lattice units) for which the equilibrium
/// `w_i φ (1 + 3 c_i u)` stays non-negative.
const MAX_ADVECTION_VELOCITY: f64 = 1.0 / 3.0;

/// Multi-level D1Q3 advection-diffusion LBM with spatial and temporal
/// refinement.
///
/// The transported scalar `φ` obeys `∂φ/∂t + u ∂φ/∂x = D ∂²φ/∂x²`.  Every
/// level is a complete periodic D1Q3 BGK lattice of the *whole* domain:
///
/// - the coarse level has `n_coarse` cells of size `dx = 1`, time step
///   `dt = 1` and relaxation time `τ_c = 1/ω_c`;
/// - refinement level `l` has `r_l = refinement_levels[l]` times finer cells
///   (`n_coarse * r_l` of them) and, following acoustic scaling, a time step
///   `dt / r_l`, so the lattice speed (and hence the advection velocity `u`)
///   is the same on every level.  Its relaxation time is
///   `τ_f = r_l (τ_c - 1/2) + 1/2`, which gives all levels the same diffusivity
///   `D = cs² (τ - 1/2) dx² / dt = cs² (τ_c - 1/2)`.
///
/// [`AdaptiveLbm::step`] advances the coarse level once and every refinement
/// level `r_l` times.  The levels are coupled only through the two transfer
/// operators, which also rescale the non-equilibrium populations with the
/// Dupuis-Chopard factors `τ_f / (r τ_c)` (coarse to fine) and `r τ_c / τ_f`
/// (fine to coarse) and conserve the scalar mass exactly:
///
/// - [`AdaptiveLbm::refine_region`] (re)initialises the refinement level in a
///   window of coarse cells from the coarse state (conservative, slope-limited
///   piecewise-linear prolongation);
/// - [`AdaptiveLbm::coarsen_region`] writes the refinement level back into the
///   coarse level (cell averages).
///
/// Cells of a refinement level that were never refined hold `φ = 0`; seed a
/// complete replica with `refine_region(level, 0, n_coarse)`.
///
/// # Macroscopic fields
///
/// `coarse_grid` and `fine_grids` expose the zeroth moment `φ = Σ_i g_i` of
/// the (private) populations and are refreshed by every method.  They may also
/// be written directly, e.g. to seed an initial condition: the next method call
/// imposes the written values on the populations by adding an equilibrium-shaped
/// correction, which keeps the non-equilibrium part of the populations.
///
/// # References
/// - Dupuis, A. & Chopard, B. (2003). *Phys. Rev. E* 67, 066707.
/// - Krüger, T. et al. (2017). *The Lattice Boltzmann Method*. Springer.
#[derive(Debug, Clone)]
pub struct AdaptiveLbm {
    /// Refinement factor at each level (e.g. `[2, 4]` means 2× and 4×).
    pub refinement_levels: Vec<usize>,
    /// Coarse-level scalar field `φ = Σ g_i` (one value per coarse cell).
    pub coarse_grid: Vec<f64>,
    /// Scalar field `φ = Σ g_i` of every refinement level, covering the whole
    /// domain with `n_coarse * refinement_levels[level]` cells.
    pub fine_grids: Vec<Vec<f64>>,
    /// Number of coarse cells.
    pub n_coarse: usize,
    /// BGK relaxation rate `ω = 1/τ` on the coarse level (`0 < ω < 2`).
    pub omega_coarse: f64,
    /// Advection velocity in lattice units, identical on every level
    /// (`|u| <= 1/3` keeps the equilibrium populations non-negative).
    pub advection_velocity: f64,
    /// D1Q3 populations of the coarse level.
    pub(super) coarse_pops: Vec<[f64; Q]>,
    /// D1Q3 populations of every refinement level.
    pub(super) fine_pops: Vec<Vec<[f64; Q]>>,
}

/// Minmod slope limiter: the smaller one-sided difference, or zero at a
/// local extremum.
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// Impose the scalar `field` on `pops` wherever it differs from the zeroth
/// moment of the populations (an external write), keeping the
/// non-equilibrium part.
fn absorb_level_writes(pops: &mut [[f64; Q]], field: &[f64], u: f64) {
    for (g, &phi) in pops.iter_mut().zip(field.iter()) {
        let current = zeroth_moment(g);
        if phi != current {
            let shift = advection_equilibrium(phi - current, u);
            for (gi, s) in g.iter_mut().zip(shift.iter()) {
                *gi += s;
            }
        }
    }
}

/// One BGK collision followed by periodic streaming of a whole level.
fn advance_level(pops: &mut [[f64; Q]], omega: f64, u: f64) {
    for g in pops.iter_mut() {
        let geq = advection_equilibrium(zeroth_moment(g), u);
        bgk_relax(g, &geq, omega);
    }
    stream_periodic(pops);
}

impl AdaptiveLbm {
    /// Create a new adaptive LBM with `n_coarse` cells and given refinement
    /// levels.  All levels start with `φ = 0` and no advection.
    ///
    /// # Errors
    /// Fails if `n_coarse` is zero, a refinement factor is zero, `omega_coarse`
    /// is not finite and positive, or `ω >= 2` (`τ <= 1/2`, unstable BGK).
    pub fn new(n_coarse: usize, refinement_levels: Vec<usize>, omega_coarse: f64) -> Result<Self> {
        if n_coarse == 0 {
            return Err(Error::invalid_grid(n_coarse, 1));
        }
        if let Some(&factor) = refinement_levels.iter().find(|&&r| r == 0) {
            return Err(Error::IncompatibleParameters {
                message: format!("refinement factor {factor} must be at least 1"),
            });
        }
        validate_relaxation_rate(omega_coarse, "omega_coarse")?;
        let fine_grids = refinement_levels
            .iter()
            .map(|&r| vec![0.0; n_coarse * r])
            .collect();
        let fine_pops = refinement_levels
            .iter()
            .map(|&r| vec![[0.0; Q]; n_coarse * r])
            .collect();
        Ok(Self {
            refinement_levels,
            coarse_grid: vec![0.0; n_coarse],
            fine_grids,
            n_coarse,
            omega_coarse,
            advection_velocity: 0.0,
            coarse_pops: vec![[0.0; Q]; n_coarse],
            fine_pops,
        })
    }

    /// Set the advection velocity (lattice units, shared by all levels).
    ///
    /// # Errors
    /// Fails unless `|u| <= 1/3`.
    pub fn set_advection_velocity(&mut self, u: f64) -> Result<()> {
        Error::check_range(
            u,
            -MAX_ADVECTION_VELOCITY,
            MAX_ADVECTION_VELOCITY,
            "advection_velocity",
        )?;
        self.advection_velocity = u;
        Ok(())
    }

    /// Relaxation time `τ_c = 1/ω_c` of the coarse level.
    pub fn relaxation_time(&self) -> f64 {
        1.0 / self.omega_coarse
    }

    /// Relaxation time `τ_f = r (τ_c - 1/2) + 1/2` of refinement level `level`
    /// (acoustic scaling), or `None` for an unknown level.
    pub fn fine_relaxation_time(&self, level: usize) -> Option<f64> {
        let &r = self.refinement_levels.get(level)?;
        Some(fine_relaxation_time(self.relaxation_time(), r))
    }

    /// Diffusivity `D = cs² (τ_c - 1/2)` in coarse lattice units; it is the
    /// same on every refinement level.
    pub fn diffusivity(&self) -> f64 {
        CS2 * (self.relaxation_time() - 0.5)
    }

    /// Total scalar mass `Σ φ dx` of the coarse level.
    pub fn total_mass(&self) -> f64 {
        self.coarse_grid.iter().sum()
    }

    /// Total scalar mass of refinement level `level`: `Σ φ dx` with
    /// `dx = 1/r` in coarse cell units (`None` for an unknown level).
    pub fn level_mass(&self, level: usize) -> Option<f64> {
        let &r = self.refinement_levels.get(level)?;
        let field = self.fine_grids.get(level)?;
        Some(field.iter().sum::<f64>() / r as f64)
    }

    /// Impose externally written values of the public scalar fields on the
    /// populations of every level.
    fn absorb_field_writes(&mut self) {
        let u = self.advection_velocity;
        absorb_level_writes(&mut self.coarse_pops, &self.coarse_grid, u);
        for (pops, field) in self.fine_pops.iter_mut().zip(self.fine_grids.iter()) {
            absorb_level_writes(pops, field, u);
        }
    }

    /// Refresh the public scalar fields from the populations.
    fn refresh_fields(&mut self) {
        for (phi, g) in self.coarse_grid.iter_mut().zip(self.coarse_pops.iter()) {
            *phi = zeroth_moment(g);
        }
        for (field, pops) in self.fine_grids.iter_mut().zip(self.fine_pops.iter()) {
            for (phi, g) in field.iter_mut().zip(pops.iter()) {
                *phi = zeroth_moment(g);
            }
        }
    }

    /// Refine coarse cells `start..end` on level `level`: initialise the
    /// refinement level in this window from the coarse level.
    ///
    /// Each coarse cell becomes `r` fine cells.  The scalar is prolonged with a
    /// conservative, minmod-limited piecewise-linear reconstruction (the `r`
    /// fine values average to the coarse value and create no new extrema), and
    /// the non-equilibrium populations are multiplied by the Dupuis-Chopard
    /// factor `τ_f / (r τ_c)`.  An unknown `level` is ignored.
    pub fn refine_region(&mut self, level: usize, start: usize, end: usize) {
        self.absorb_field_writes();
        let tau_c = self.relaxation_time();
        let u = self.advection_velocity;
        let n = self.n_coarse.min(self.coarse_pops.len());
        if let (Some(&r), Some(fine)) = (
            self.refinement_levels.get(level),
            self.fine_pops.get_mut(level),
        ) {
            let tau_f = fine_relaxation_time(tau_c, r);
            let scale = coarse_to_fine_scale(tau_c, tau_f, r);
            for ci in start..end.min(n) {
                let g = self.coarse_pops[ci];
                let phi = zeroth_moment(&g);
                let neq = nonequilibrium(&g, &advection_equilibrium(phi, u));
                let left = zeroth_moment(&self.coarse_pops[(ci + n - 1) % n]);
                let right = zeroth_moment(&self.coarse_pops[(ci + 1) % n]);
                let slope = minmod(phi - left, right - phi);
                for k in 0..r {
                    // Offsets of the fine cell centres from the coarse centre,
                    // in coarse cells; they are antisymmetric and sum to zero.
                    let offset = (2 * k + 1) as f64 / (2 * r) as f64 - 0.5;
                    let eq = advection_equilibrium(phi + slope * offset, u);
                    if let Some(node) = fine.get_mut(ci * r + k) {
                        for q in 0..Q {
                            node[q] = eq[q] + scale * neq[q];
                        }
                    }
                }
            }
        }
        self.refresh_fields();
    }

    /// Coarsen cells `start..end` on level `level`: write the refinement level
    /// back into the coarse level.
    ///
    /// The scalar of a coarse cell is the average of its `r` fine cells (mass
    /// conserving) and the average non-equilibrium population is multiplied by
    /// the Dupuis-Chopard factor `r τ_c / τ_f`.  An unknown `level` is ignored.
    pub fn coarsen_region(&mut self, level: usize, start: usize, end: usize) {
        self.absorb_field_writes();
        let tau_c = self.relaxation_time();
        let u = self.advection_velocity;
        let n = self.n_coarse.min(self.coarse_pops.len());
        if let (Some(&r), Some(fine)) =
            (self.refinement_levels.get(level), self.fine_pops.get(level))
        {
            let tau_f = fine_relaxation_time(tau_c, r);
            let scale = fine_to_coarse_scale(tau_c, tau_f, r);
            for ci in start..end.min(n) {
                let mut phi_sum = 0.0;
                let mut neq_sum = [0.0_f64; Q];
                for k in 0..r {
                    if let Some(g) = fine.get(ci * r + k) {
                        let phi = zeroth_moment(g);
                        let neq = nonequilibrium(g, &advection_equilibrium(phi, u));
                        phi_sum += phi;
                        for q in 0..Q {
                            neq_sum[q] += neq[q];
                        }
                    }
                }
                let inv_r = 1.0 / r as f64;
                let eq = advection_equilibrium(phi_sum * inv_r, u);
                for q in 0..Q {
                    self.coarse_pops[ci][q] = eq[q] + scale * neq_sum[q] * inv_r;
                }
            }
        }
        self.refresh_fields();
    }

    /// Linear time interpolation between two coarse-grid states.
    ///
    /// `alpha` ∈ \[0, 1\] gives the interpolation weight (0 = old, 1 = new).
    pub fn time_interpolation(&self, old_state: &[f64], new_state: &[f64], alpha: f64) -> Vec<f64> {
        let a = clamp(alpha, 0.0, 1.0);
        old_state
            .iter()
            .zip(new_state.iter())
            .map(|(&o, &n)| (1.0 - a) * o + a * n)
            .collect()
    }

    /// Spatial interpolation from coarse to fine grid using linear interpolation.
    pub fn spatial_interpolation(&self, coarse: &[f64], refinement: usize) -> Vec<f64> {
        let n = coarse.len();
        if refinement == 0 {
            return coarse.to_vec();
        }
        let fine_n = n * refinement;
        let mut fine = vec![0.0; fine_n];
        for ci in 0..n {
            let v_left = coarse[ci];
            let v_right = coarse[(ci + 1).min(n - 1)];
            for k in 0..refinement {
                let t = k as f64 / refinement as f64;
                fine[ci * refinement + k] = (1.0 - t) * v_left + t * v_right;
            }
        }
        fine
    }

    /// Advance the solver by one coarse time step.
    ///
    /// The coarse level performs one BGK collision and periodic streaming
    /// step, and every refinement level performs `r` of them with its own
    /// relaxation time `τ_f = r (τ_c - 1/2) + 1/2`, so that all levels reach
    /// the same physical time.  Each level conserves its scalar mass exactly.
    pub fn step(&mut self) {
        self.absorb_field_writes();
        let u = self.advection_velocity;
        advance_level(&mut self.coarse_pops, self.omega_coarse, u);
        let tau_c = self.relaxation_time();
        for (pops, &r) in self.fine_pops.iter_mut().zip(self.refinement_levels.iter()) {
            let omega_fine = 1.0 / fine_relaxation_time(tau_c, r);
            for _ in 0..r {
                advance_level(pops, omega_fine, u);
            }
        }
        self.refresh_fields();
    }
}
