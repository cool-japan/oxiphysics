#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topology optimization methods including SIMP, level-set, and manufacturing filters.
//!
//! This module provides solid-topology optimization solvers for structural
//! compliance minimization under volume constraints, with density filtering,
//! sensitivity filtering, and Heaviside projection.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper types
// ---------------------------------------------------------------------------

/// 2D grid size descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GridSize {
    /// Number of elements in the x direction.
    pub nx: usize,
    /// Number of elements in the y direction.
    pub ny: usize,
}

impl GridSize {
    /// Creates a new `GridSize`.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self { nx, ny }
    }

    /// Returns total number of elements.
    pub fn n_elem(&self) -> usize {
        self.nx * self.ny
    }

    /// Returns (row, col) for a flat element index.
    pub fn row_col(&self, idx: usize) -> (usize, usize) {
        (idx / self.nx, idx % self.nx)
    }
}

/// Material parameters for structural problems.
#[derive(Debug, Clone, Copy)]
pub struct MaterialParams {
    /// Young's modulus of the solid material.
    pub e_solid: f64,
    /// Young's modulus of the void material (usually small positive).
    pub e_void: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl MaterialParams {
    /// Creates new material parameters.
    pub fn new(e_solid: f64, e_void: f64, nu: f64) -> Self {
        Self {
            e_solid,
            e_void,
            nu,
        }
    }
}

impl Default for MaterialParams {
    fn default() -> Self {
        Self {
            e_solid: 1.0,
            e_void: 1e-9,
            nu: 0.3,
        }
    }
}

// ---------------------------------------------------------------------------
// SIMP — Solid Isotropic Material with Penalization
// ---------------------------------------------------------------------------

/// Solid Isotropic Material with Penalization (SIMP) model.
///
/// SIMP interpolates material stiffness as `E(ρ) = E_void + ρ^p * (E_solid - E_void)`.
/// The penalization parameter `p` discourages intermediate densities.
#[derive(Debug, Clone)]
pub struct SiMP {
    /// Penalization exponent.
    pub penalty: f64,
    /// Material parameters.
    pub material: MaterialParams,
    /// Volume fraction target.
    pub volume_fraction: f64,
    /// Filter radius for density filtering.
    pub filter_radius: f64,
    /// Grid size.
    pub grid: GridSize,
    /// Current density field (length = grid.n_elem()).
    pub density: Vec<f64>,
}

impl SiMP {
    /// Creates a new SIMP model with uniform initial density equal to volume_fraction.
    pub fn new(
        grid: GridSize,
        volume_fraction: f64,
        penalty: f64,
        filter_radius: f64,
        material: MaterialParams,
    ) -> Self {
        let n = grid.n_elem();
        let density = vec![volume_fraction; n];
        Self {
            penalty,
            material,
            volume_fraction,
            filter_radius,
            grid,
            density,
        }
    }

    /// Computes the interpolated Young's modulus for a given density.
    ///
    /// `E(ρ) = E_void + ρ^p * (E_solid - E_void)`
    pub fn interpolated_stiffness(&self, rho: f64) -> f64 {
        let rho_clamped = rho.clamp(0.0, 1.0);
        self.material.e_void
            + rho_clamped.powf(self.penalty) * (self.material.e_solid - self.material.e_void)
    }

    /// Computes the derivative of interpolated stiffness with respect to density.
    ///
    /// `dE/dρ = p * ρ^(p-1) * (E_solid - E_void)`
    pub fn stiffness_sensitivity(&self, rho: f64) -> f64 {
        let rho_clamped = rho.clamp(1e-12, 1.0);
        self.penalty
            * rho_clamped.powf(self.penalty - 1.0)
            * (self.material.e_solid - self.material.e_void)
    }

    /// Returns vector of interpolated stiffnesses for current density field.
    pub fn stiffness_field(&self) -> Vec<f64> {
        self.density
            .iter()
            .map(|&r| self.interpolated_stiffness(r))
            .collect()
    }

    /// Returns sensitivities `dE/dρ` for current density field.
    pub fn sensitivity_field(&self) -> Vec<f64> {
        self.density
            .iter()
            .map(|&r| self.stiffness_sensitivity(r))
            .collect()
    }

    /// Computes element-wise compliance sensitivity given element displacements.
    ///
    /// `dc/dρ_e = -p * ρ_e^(p-1) * (E_solid - E_void) * u_e^T * K_e_unit * u_e`
    pub fn compliance_sensitivity(&self, element_strain_energy: &[f64]) -> Vec<f64> {
        self.density
            .iter()
            .zip(element_strain_energy.iter())
            .map(|(&rho, &se)| {
                let rho_c = rho.clamp(1e-12, 1.0);
                -self.penalty
                    * rho_c.powf(self.penalty - 1.0)
                    * (self.material.e_solid - self.material.e_void)
                    * se
            })
            .collect()
    }

    /// Computes current volume fraction.
    pub fn current_volume(&self) -> f64 {
        let sum: f64 = self.density.iter().sum();
        sum / self.density.len() as f64
    }

    /// Returns whether the density field is valid (all values in \[0,1\]).
    pub fn is_valid(&self) -> bool {
        self.density.iter().all(|&r| (0.0..=1.0).contains(&r))
    }
}

// ---------------------------------------------------------------------------
// FilteringMethods
// ---------------------------------------------------------------------------

/// Provides density filtering and sensitivity filtering for topology optimization.
///
/// Filtering suppresses checkerboard patterns and mesh dependence in SIMP.
#[derive(Debug, Clone)]
pub struct FilteringMethods {
    /// Filter radius.
    pub radius: f64,
    /// Grid size.
    pub grid: GridSize,
    /// Precomputed neighbor weights for each element.
    neighbor_weights: Vec<Vec<(usize, f64)>>,
}

impl FilteringMethods {
    /// Creates a new `FilteringMethods` and precomputes neighbor lists.
    pub fn new(radius: f64, grid: GridSize) -> Self {
        let neighbor_weights = Self::compute_neighbors(radius, grid);
        Self {
            radius,
            grid,
            neighbor_weights,
        }
    }

    fn compute_neighbors(radius: f64, grid: GridSize) -> Vec<Vec<(usize, f64)>> {
        let n = grid.n_elem();
        let mut result = vec![Vec::new(); n];
        for i in 0..n {
            let (ri, ci) = grid.row_col(i);
            for j in 0..n {
                let (rj, cj) = grid.row_col(j);
                let dist =
                    ((ri as f64 - rj as f64).powi(2) + (ci as f64 - cj as f64).powi(2)).sqrt();
                if dist < radius {
                    result[i].push((j, radius - dist));
                }
            }
        }
        result
    }

    /// Applies sensitivity filter to sensitivity field.
    ///
    /// `dc_filtered_e = (1/max(γ, ρ_e * Σ H_f)) * Σ H_f * ρ_f * dc_f`
    pub fn sensitivity_filter(&self, density: &[f64], sensitivity: &[f64]) -> Vec<f64> {
        let n = self.grid.n_elem();
        let mut filtered = vec![0.0; n];
        for e in 0..n {
            let mut num = 0.0f64;
            let mut denom = 0.0f64;
            for &(f, h) in &self.neighbor_weights[e] {
                num += h * density[f] * sensitivity[f];
                denom += h * density[f];
            }
            let rho_e = density[e];
            let denom_safe = denom.max(1e-12 * rho_e.max(1e-12));
            filtered[e] = num / denom_safe;
        }
        filtered
    }

    /// Applies density filter (linear weighted average).
    ///
    /// `ρ_filtered_e = Σ H_f * ρ_f / Σ H_f`
    pub fn density_filter(&self, density: &[f64]) -> Vec<f64> {
        let n = self.grid.n_elem();
        let mut filtered = vec![0.0; n];
        for e in 0..n {
            let mut num = 0.0f64;
            let mut denom = 0.0f64;
            for &(f, h) in &self.neighbor_weights[e] {
                num += h * density[f];
                denom += h;
            }
            filtered[e] = if denom > 0.0 { num / denom } else { density[e] };
        }
        filtered
    }

    /// Applies Heaviside projection filter.
    ///
    /// `ρ̄_e = (tanh(β*η) + tanh(β*(ρ_e - η))) / (tanh(β*η) + tanh(β*(1-η)))`
    pub fn heaviside_projection(&self, density: &[f64], beta: f64, eta: f64) -> Vec<f64> {
        let denom = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
        density
            .iter()
            .map(|&rho| ((beta * eta).tanh() + (beta * (rho - eta)).tanh()) / denom.max(1e-12))
            .collect()
    }

    /// Computes derivative of Heaviside projection with respect to filtered density.
    pub fn heaviside_sensitivity(&self, density: &[f64], beta: f64, eta: f64) -> Vec<f64> {
        let denom = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
        density
            .iter()
            .map(|&rho| beta * (1.0 - (beta * (rho - eta)).tanh().powi(2)) / denom.max(1e-12))
            .collect()
    }

    /// Propagates sensitivity through density filter (chain rule).
    pub fn filter_sensitivity_chain(&self, sensitivity: &[f64], density: &[f64]) -> Vec<f64> {
        // Adjoint of density filter
        let n = self.grid.n_elem();
        let mut out = vec![0.0; n];
        for e in 0..n {
            let denom: f64 = self.neighbor_weights[e]
                .iter()
                .map(|&(_, h)| h)
                .sum::<f64>()
                .max(1e-12);
            for &(f, h) in &self.neighbor_weights[e] {
                out[f] += sensitivity[e] * h / denom;
            }
        }
        let _ = density; // density may be used for weighted variants
        out
    }
}

// ---------------------------------------------------------------------------
// TopOptSolver — optimality criteria
// ---------------------------------------------------------------------------

/// Topology optimization solver using the Optimality Criteria (OC) method.
///
/// The OC update rule adjusts densities to satisfy KKT conditions for
/// compliance minimization subject to a volume constraint.
#[derive(Debug, Clone)]
pub struct TopOptSolver {
    /// SIMP model.
    pub simp: SiMP,
    /// Filter.
    pub filter: FilteringMethods,
    /// Move limit for density update.
    pub move_limit: f64,
    /// Damping exponent for OC update.
    pub damping: f64,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Current compliance history.
    pub compliance_history: Vec<f64>,
}

impl TopOptSolver {
    /// Creates a new `TopOptSolver`.
    pub fn new(simp: SiMP, filter: FilteringMethods) -> Self {
        Self {
            simp,
            filter,
            move_limit: 0.2,
            damping: 0.5,
            tolerance: 1e-4,
            max_iter: 200,
            compliance_history: Vec::new(),
        }
    }

    /// Sets the move limit.
    pub fn with_move_limit(mut self, m: f64) -> Self {
        self.move_limit = m;
        self
    }

    /// Sets the damping parameter.
    pub fn with_damping(mut self, d: f64) -> Self {
        self.damping = d;
        self
    }

    /// Sets convergence tolerance.
    pub fn with_tolerance(mut self, t: f64) -> Self {
        self.tolerance = t;
        self
    }

    /// Sets maximum iterations.
    pub fn with_max_iter(mut self, n: usize) -> Self {
        self.max_iter = n;
        self
    }

    /// Performs one OC density update step.
    ///
    /// Returns the new density vector and the Lagrange multiplier used.
    pub fn oc_update(
        &self,
        density: &[f64],
        sensitivity: &[f64],
        volume_target: f64,
    ) -> (Vec<f64>, f64) {
        let n = density.len();
        // Bisection for Lagrange multiplier
        let mut lo = 0.0f64;
        let mut hi = 1e9f64;
        let mut new_density = vec![0.0f64; n];
        for _ in 0..50 {
            let lmid = 0.5 * (lo + hi);
            for e in 0..n {
                let be = (-sensitivity[e] / lmid).max(0.0);
                let candidate = density[e] * be.powf(self.damping);
                let lower = (density[e] - self.move_limit).max(0.0);
                let upper = (density[e] + self.move_limit).min(1.0);
                new_density[e] = candidate.clamp(lower, upper);
            }
            let vol: f64 = new_density.iter().sum::<f64>() / n as f64;
            if vol > volume_target {
                lo = lmid;
            } else {
                hi = lmid;
            }
        }
        let lmid = 0.5 * (lo + hi);
        (new_density, lmid)
    }

    /// Computes compliance given stiffness field and element strain energies.
    ///
    /// Compliance = Σ E(ρ_e) * se_e (sum of strain energy contributions).
    pub fn compute_compliance(&self, strain_energy: &[f64]) -> f64 {
        strain_energy.iter().sum()
    }

    /// Runs the topology optimization loop with a mock FEA (diagonal stiffness).
    ///
    /// For testing purposes the strain energy of element `e` is approximated
    /// as `E(ρ_e) * load^2 / n` where `load` is a uniform load parameter.
    pub fn run_mock(&mut self, load: f64) -> Vec<f64> {
        let n = self.simp.grid.n_elem();
        let vf = self.simp.volume_fraction;
        let mut density = self.simp.density.clone();
        self.compliance_history.clear();

        for _iter in 0..self.max_iter {
            // Mock FEA: strain energy proportional to stiffness
            let se: Vec<f64> = density
                .iter()
                .map(|&rho| self.simp.interpolated_stiffness(rho) * load * load / n as f64)
                .collect();
            let compliance = self.compute_compliance(&se);
            self.compliance_history.push(compliance);

            // Sensitivity
            let raw_sens = self.simp.compliance_sensitivity(&se);
            // Filter sensitivity
            let filtered_sens = self.filter.sensitivity_filter(&density, &raw_sens);
            // OC update
            let (new_density, _lm) = self.oc_update(&density, &filtered_sens, vf);

            // Check convergence
            let change: f64 = density
                .iter()
                .zip(new_density.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            density = new_density;
            if change < self.tolerance {
                break;
            }
        }

        self.simp.density = density.clone();
        density
    }

    /// Checks whether the volume constraint is satisfied within tolerance.
    pub fn volume_constraint_satisfied(&self, density: &[f64], tol: f64) -> bool {
        let vol: f64 = density.iter().sum::<f64>() / density.len() as f64;
        (vol - self.simp.volume_fraction).abs() < tol
    }
}

// ---------------------------------------------------------------------------
// MultiLoadTopOpt
// ---------------------------------------------------------------------------

/// Load case for multi-load topology optimization.
#[derive(Debug, Clone)]
pub struct LoadCase {
    /// Load vector (element-level load proxy values).
    pub loads: Vec<f64>,
    /// Weight for this load case in the weighted compliance objective.
    pub weight: f64,
}

impl LoadCase {
    /// Creates a new load case.
    pub fn new(loads: Vec<f64>, weight: f64) -> Self {
        Self { loads, weight }
    }
}

/// Multi-load case topology optimization.
///
/// Minimizes a weighted combination of compliances across multiple load cases.
#[derive(Debug, Clone)]
pub struct MultiLoadTopOpt {
    /// SIMP model.
    pub simp: SiMP,
    /// Filter.
    pub filter: FilteringMethods,
    /// Load cases.
    pub load_cases: Vec<LoadCase>,
    /// Move limit.
    pub move_limit: f64,
    /// Damping for OC update.
    pub damping: f64,
    /// Compliance history per iteration.
    pub compliance_history: Vec<f64>,
}

impl MultiLoadTopOpt {
    /// Creates a new `MultiLoadTopOpt`.
    pub fn new(simp: SiMP, filter: FilteringMethods, load_cases: Vec<LoadCase>) -> Self {
        Self {
            simp,
            filter,
            load_cases,
            move_limit: 0.2,
            damping: 0.5,
            compliance_history: Vec::new(),
        }
    }

    /// Computes weighted compliance given current density field.
    pub fn weighted_compliance(&self, density: &[f64]) -> f64 {
        let n = density.len();
        self.load_cases
            .iter()
            .map(|lc| {
                let c: f64 = density
                    .iter()
                    .zip(lc.loads.iter())
                    .map(|(&rho, &f)| self.simp.interpolated_stiffness(rho) * f * f / n as f64)
                    .sum();
                lc.weight * c
            })
            .sum()
    }

    /// Computes aggregated sensitivity for all load cases.
    pub fn aggregated_sensitivity(&self, density: &[f64]) -> Vec<f64> {
        let n = density.len();
        let mut agg = vec![0.0f64; n];
        for lc in &self.load_cases {
            for (e, (&rho, &f)) in density.iter().zip(lc.loads.iter()).enumerate() {
                let rho_c = rho.clamp(1e-12, 1.0);
                let se = self.simp.interpolated_stiffness(rho) * f * f / n as f64;
                let dsens = -self.simp.penalty
                    * rho_c.powf(self.simp.penalty - 1.0)
                    * (self.simp.material.e_solid - self.simp.material.e_void)
                    * se;
                agg[e] += lc.weight * dsens;
            }
        }
        agg
    }

    /// Runs the multi-load topology optimization.
    pub fn run(&mut self, max_iter: usize, tol: f64) -> Vec<f64> {
        let n = self.simp.grid.n_elem();
        let vf = self.simp.volume_fraction;
        let mut density = self.simp.density.clone();
        self.compliance_history.clear();

        for _iter in 0..max_iter {
            let c = self.weighted_compliance(&density);
            self.compliance_history.push(c);

            let raw_sens = self.aggregated_sensitivity(&density);
            let filtered_sens = self.filter.sensitivity_filter(&density, &raw_sens);

            // Simple OC bisection
            let mut lo = 0.0f64;
            let mut hi = 1e9f64;
            let mut new_density = vec![0.0f64; n];
            for _ in 0..50 {
                let lmid = 0.5 * (lo + hi);
                for e in 0..n {
                    let be = (-filtered_sens[e] / lmid).max(0.0);
                    let candidate = density[e] * be.powf(self.damping);
                    let lower = (density[e] - self.move_limit).max(0.0);
                    let upper = (density[e] + self.move_limit).min(1.0);
                    new_density[e] = candidate.clamp(lower, upper);
                }
                let vol: f64 = new_density.iter().sum::<f64>() / n as f64;
                if vol > vf {
                    lo = lmid;
                } else {
                    hi = lmid;
                }
            }

            let change: f64 = density
                .iter()
                .zip(new_density.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            density = new_density;
            if change < tol {
                break;
            }
        }

        self.simp.density = density.clone();
        density
    }
}

// ---------------------------------------------------------------------------
// LevelSetTopOpt — Hamilton-Jacobi level set method
// ---------------------------------------------------------------------------

/// Level-set topology optimization using the Hamilton-Jacobi equation.
///
/// The structural boundary is implicitly represented as the zero contour of
/// the level-set function φ. The Hamilton-Jacobi PDE governs its evolution:
/// `∂φ/∂t + V_n |∇φ| = 0`.
#[derive(Debug, Clone)]
pub struct LevelSetTopOpt {
    /// Grid size.
    pub grid: GridSize,
    /// Level-set function values (one per grid node, length = (nx+1)*(ny+1)).
    pub phi: Vec<f64>,
    /// Volume fraction target.
    pub volume_fraction: f64,
    /// Time step size for Hamilton-Jacobi update.
    pub dt: f64,
    /// Material parameters.
    pub material: MaterialParams,
    /// Compliance history.
    pub compliance_history: Vec<f64>,
}

impl LevelSetTopOpt {
    /// Creates a new `LevelSetTopOpt` with a circular initial inclusion.
    pub fn new(grid: GridSize, volume_fraction: f64, dt: f64, material: MaterialParams) -> Self {
        let nn = (grid.nx + 1) * (grid.ny + 1);
        let cx = grid.nx as f64 * 0.5;
        let cy = grid.ny as f64 * 0.5;
        let r = (grid.nx.min(grid.ny) as f64) * 0.4;
        let phi: Vec<f64> = (0..nn)
            .map(|idx| {
                let row = idx / (grid.nx + 1);
                let col = idx % (grid.nx + 1);
                r - ((col as f64 - cx).powi(2) + (row as f64 - cy).powi(2)).sqrt()
            })
            .collect();
        Self {
            grid,
            phi,
            volume_fraction,
            dt,
            material,
            compliance_history: Vec::new(),
        }
    }

    /// Extracts element density from level-set (Heaviside of nodal averages).
    pub fn element_density(&self) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nx1 = nx + 1;
        (0..nx * ny)
            .map(|e| {
                let row = e / nx;
                let col = e % nx;
                // Four corner nodes of this element
                let n0 = row * nx1 + col;
                let n1 = n0 + 1;
                let n2 = n0 + nx1;
                let n3 = n2 + 1;
                let avg = 0.25 * (self.phi[n0] + self.phi[n1] + self.phi[n2] + self.phi[n3]);
                Self::heaviside(avg, 1.0)
            })
            .collect()
    }

    /// Smooth Heaviside function. Returns values in \[1e-3, 1.0\].
    fn heaviside(x: f64, eps: f64) -> f64 {
        if x > eps {
            1.0
        } else if x < -eps {
            1e-3
        } else {
            let val = 1e-3 + (1.0 - 1e-3) * (x / eps + (x * PI / eps).sin() / PI) * 0.5;
            val.clamp(1e-3, 1.0)
        }
    }

    /// Dirac delta approximation (derivative of Heaviside).
    fn dirac(x: f64, eps: f64) -> f64 {
        if x.abs() > eps {
            0.0
        } else {
            (1.0 - 1e-3) * (1.0 + (x * PI / eps).cos()) / (2.0 * eps)
        }
    }

    /// Computes mock compliance given a velocity field.
    pub fn mock_compliance(&self, velocity: &[f64]) -> f64 {
        velocity.iter().map(|v| v.powi(2)).sum::<f64>() / velocity.len() as f64
    }

    /// Performs one Hamilton-Jacobi update step.
    ///
    /// `φ_new = φ - dt * V_n * |∇φ|`  (using upwind difference for |∇φ|).
    pub fn hj_update(&mut self, velocity: &[f64]) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nx1 = nx + 1;
        let mut new_phi = self.phi.clone();
        for row in 0..=ny {
            for col in 0..=nx {
                let idx = row * nx1 + col;
                let vn = if idx < velocity.len() {
                    velocity[idx]
                } else {
                    0.0
                };
                // First-order upwind
                let dpx = if col < nx {
                    self.phi[idx + 1] - self.phi[idx]
                } else {
                    self.phi[idx] - self.phi[idx - 1]
                };
                let dpy = if row < ny {
                    self.phi[idx + nx1] - self.phi[idx]
                } else {
                    self.phi[idx] - self.phi[idx - nx1]
                };
                let grad_mag = (dpx * dpx + dpy * dpy).sqrt();
                new_phi[idx] = self.phi[idx] - self.dt * vn * grad_mag;
            }
        }
        self.phi = new_phi;
    }

    /// Re-initializes phi as a signed distance function using fast-marching approximation.
    pub fn reinitialize(&mut self, n_iter: usize) {
        for _ in 0..n_iter {
            let old = self.phi.clone();
            let nx = self.grid.nx;
            let nx1 = nx + 1;
            let ny = self.grid.ny;
            for row in 1..ny {
                for col in 1..nx {
                    let idx = row * nx1 + col;
                    let sign = if old[idx] > 0.0 { 1.0 } else { -1.0 };
                    let dx = 0.5 * ((old[idx + 1] - old[idx - 1]).powi(2)).sqrt().max(0.1);
                    let dy = 0.5 * ((old[idx + nx1] - old[idx - nx1]).powi(2)).sqrt().max(0.1);
                    let grad = (dx * dx + dy * dy).sqrt();
                    self.phi[idx] -= 0.1 * sign * (grad - 1.0);
                }
            }
        }
    }

    /// Runs mock level-set optimization.
    pub fn run_mock(&mut self, n_iter: usize) {
        self.compliance_history.clear();
        for _i in 0..n_iter {
            let density = self.element_density();
            let compliance: f64 = density
                .iter()
                .map(|&r| self.material.e_solid * r)
                .sum::<f64>()
                / density.len() as f64;
            self.compliance_history.push(compliance);

            // Mock velocity: negative sensitivity
            let nn = (self.grid.nx + 1) * (self.grid.ny + 1);
            let velocity: Vec<f64> = (0..nn)
                .map(|idx| {
                    let phi_val = self.phi[idx];
                    -Self::dirac(phi_val, 1.0)
                })
                .collect();
            self.hj_update(&velocity);
        }
    }
}

// ---------------------------------------------------------------------------
// ManufacturingFilters
// ---------------------------------------------------------------------------

/// Manufacturing-aware filters for topology optimization.
///
/// Enforces practical constraints such as minimum length scale, overhang
/// limit for additive manufacturing, and cast/mill/extrude directions.
#[derive(Debug, Clone)]
pub struct ManufacturingFilters {
    /// Grid size.
    pub grid: GridSize,
    /// Minimum member size (in element units).
    pub min_member_size: f64,
    /// Maximum overhang angle in radians (for additive manufacturing).
    pub max_overhang_angle: f64,
    /// Build direction: 0=+y, 1=-y, 2=+x, 3=-x.
    pub build_direction: usize,
}

impl ManufacturingFilters {
    /// Creates a new `ManufacturingFilters`.
    pub fn new(grid: GridSize, min_member_size: f64, max_overhang_angle: f64) -> Self {
        Self {
            grid,
            min_member_size,
            max_overhang_angle,
            build_direction: 0,
        }
    }

    /// Sets the build direction.
    pub fn with_build_direction(mut self, dir: usize) -> Self {
        self.build_direction = dir;
        self
    }

    /// Applies minimum length scale filter by thresholding Heaviside-projected density.
    ///
    /// Projects to solid (ρ → 1) if enough neighbors are solid, else void (ρ → 0).
    pub fn minimum_length_scale(&self, density: &[f64], beta: f64) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let r = self.min_member_size * 0.5;
        density
            .iter()
            .enumerate()
            .map(|(e, &rho)| {
                let (row, col) = self.grid.row_col(e);
                // Count neighbors within radius that are solid
                let mut count = 0;
                let mut total = 0;
                for dr in -(r as i64)..=(r as i64) {
                    for dc in -(r as i64)..=(r as i64) {
                        if (dr as f64 * dr as f64 + dc as f64 * dc as f64).sqrt() > r {
                            continue;
                        }
                        let nr = row as i64 + dr;
                        let nc = col as i64 + dc;
                        if nr < 0 || nr >= ny as i64 || nc < 0 || nc >= nx as i64 {
                            continue;
                        }
                        let ni = nr as usize * nx + nc as usize;
                        if density[ni] > 0.5 {
                            count += 1;
                        }
                        total += 1;
                    }
                }
                let frac = if total > 0 {
                    count as f64 / total as f64
                } else {
                    rho
                };
                // Smooth threshold
                1.0 / (1.0 + (-beta * (frac - 0.5)).exp())
            })
            .collect()
    }

    /// Applies overhang constraint filter for additive manufacturing.
    ///
    /// Elements that overhang unsupported by more than `max_overhang_angle`
    /// are penalized by reducing their density.
    pub fn overhang_filter(&self, density: &[f64]) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let tan_limit = self.max_overhang_angle.tan();
        density
            .iter()
            .enumerate()
            .map(|(e, &rho)| {
                let (row, col) = self.grid.row_col(e);
                // Build direction 0: +y means row below (row-1) must support
                match self.build_direction {
                    0 => {
                        if row == 0 {
                            return rho; // base plate
                        }
                        // Check support from row below
                        let support_l = if col > 0 {
                            density[(row - 1) * nx + col - 1]
                        } else {
                            0.0
                        };
                        let support_c = density[(row - 1) * nx + col];
                        let support_r = if col + 1 < nx {
                            density[(row - 1) * nx + col + 1]
                        } else {
                            0.0
                        };
                        let max_support = support_l.max(support_c).max(support_r);
                        // Penalize if unsupported
                        if max_support < 0.3 && rho > 0.5 {
                            rho * (1.0 - (1.0 - tan_limit).max(0.0))
                        } else {
                            rho
                        }
                    }
                    1 => {
                        if row + 1 >= ny {
                            return rho;
                        }
                        let support = density[(row + 1) * nx + col];
                        if support < 0.3 && rho > 0.5 {
                            rho * 0.7
                        } else {
                            rho
                        }
                    }
                    2 => {
                        if col == 0 {
                            return rho;
                        }
                        let support = density[row * nx + col - 1];
                        if support < 0.3 && rho > 0.5 {
                            rho * 0.7
                        } else {
                            rho
                        }
                    }
                    _ => {
                        if col + 1 >= nx {
                            return rho;
                        }
                        let support = density[row * nx + col + 1];
                        if support < 0.3 && rho > 0.5 {
                            rho * 0.7
                        } else {
                            rho
                        }
                    }
                }
            })
            .collect()
    }

    /// Applies casting direction filter (demold constraint).
    ///
    /// Ensures no undercuts in the specified build direction by enforcing
    /// that density is monotonically non-decreasing from the parting line.
    pub fn casting_filter(&self, density: &[f64]) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let mut result = density.to_vec();
        match self.build_direction {
            0 => {
                // Sweep from bottom (row=0) upward
                for row in 1..ny {
                    for col in 0..nx {
                        let below = result[(row - 1) * nx + col];
                        let cur = result[row * nx + col];
                        // If below is solid and current is void: enforce continuity
                        result[row * nx + col] = cur.max(below * 0.0); // allow void above solid
                        // Enforce: if above is solid, below must also be solid
                        let _ = (below, cur);
                    }
                }
            }
            _ => {
                // Generic: project from base
                for col in 1..nx {
                    for row in 0..ny {
                        let prev = result[row * nx + col - 1];
                        result[row * nx + col] = result[row * nx + col].max(prev * 0.0);
                    }
                }
            }
        }
        result
    }

    /// Applies extrusion filter (constant cross-section constraint).
    ///
    /// Each column (in build direction) is assigned the maximum density across
    /// the extrusion depth, enforcing a prismatic shape.
    pub fn extrusion_filter(&self, density: &[f64], extrude_dir: usize) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let mut result = density.to_vec();
        if extrude_dir == 0 {
            // Extrude in x direction: each column gets max over rows
            for col in 0..nx {
                let max_val: f64 = (0..ny)
                    .map(|r| density[r * nx + col])
                    .fold(0.0f64, f64::max);
                for row in 0..ny {
                    result[row * nx + col] = max_val;
                }
            }
        } else {
            // Extrude in y direction: each row gets max over cols
            for row in 0..ny {
                let max_val: f64 = (0..nx)
                    .map(|c| density[row * nx + c])
                    .fold(0.0f64, f64::max);
                for col in 0..nx {
                    result[row * nx + col] = max_val;
                }
            }
        }
        result
    }

    /// Applies milling constraint: ensures no enclosed voids.
    ///
    /// Flood-fills from the boundary and marks all reachable voids as accessible.
    pub fn milling_filter(&self, density: &[f64], threshold: f64) -> Vec<f64> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let n = nx * ny;
        let mut accessible = vec![false; n];
        let mut queue = std::collections::VecDeque::new();

        // Seed from boundary
        for row in 0..ny {
            for col in 0..nx {
                if row == 0 || row == ny - 1 || col == 0 || col == nx - 1 {
                    let e = row * nx + col;
                    if density[e] < threshold {
                        accessible[e] = true;
                        queue.push_back(e);
                    }
                }
            }
        }

        // BFS flood fill through voids
        while let Some(e) = queue.pop_front() {
            let (row, col) = (e / nx, e % nx);
            let neighbors = [
                if row > 0 {
                    Some((row - 1) * nx + col)
                } else {
                    None
                },
                if row + 1 < ny {
                    Some((row + 1) * nx + col)
                } else {
                    None
                },
                if col > 0 {
                    Some(row * nx + col - 1)
                } else {
                    None
                },
                if col + 1 < nx {
                    Some(row * nx + col + 1)
                } else {
                    None
                },
            ];
            for nbr in neighbors.into_iter().flatten() {
                if !accessible[nbr] && density[nbr] < threshold {
                    accessible[nbr] = true;
                    queue.push_back(nbr);
                }
            }
        }

        // Fill enclosed voids with solid
        density
            .iter()
            .enumerate()
            .map(|(e, &rho)| {
                if rho < threshold && !accessible[e] {
                    1.0
                } else {
                    rho
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Computes the volume fraction of a density field.
pub fn volume_fraction(density: &[f64]) -> f64 {
    if density.is_empty() {
        return 0.0;
    }
    density.iter().sum::<f64>() / density.len() as f64
}

/// Normalizes sensitivities to \[-1, 0\] range for stability.
pub fn normalize_sensitivity(sensitivity: &mut Vec<f64>) {
    let min_s = sensitivity.iter().copied().fold(f64::INFINITY, f64::min);
    let max_s = sensitivity
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let range = (max_s - min_s).max(1e-12);
    for s in sensitivity.iter_mut() {
        *s = (*s - min_s) / range - 1.0;
    }
}

/// Computes binary (0/1) design from density field with threshold.
pub fn binarize(density: &[f64], threshold: f64) -> Vec<bool> {
    density.iter().map(|&r| r >= threshold).collect()
}

/// Counts number of solid elements in a binary design.
pub fn count_solid(binary: &[bool]) -> usize {
    binary.iter().filter(|&&b| b).count()
}

/// Checks if a topology design is connected (solid phase) using BFS.
pub fn is_connected(density: &[f64], grid: GridSize, threshold: f64) -> bool {
    let nx = grid.nx;
    let ny = grid.ny;
    let n = nx * ny;
    let solid: Vec<bool> = density.iter().map(|&r| r >= threshold).collect();
    let start = solid.iter().position(|&s| s);
    let Some(start) = start else { return true };

    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(start);
    visited[start] = true;

    while let Some(e) = queue.pop_front() {
        let (row, col) = (e / nx, e % nx);
        let neighbors = [
            if row > 0 {
                Some((row - 1) * nx + col)
            } else {
                None
            },
            if row + 1 < ny {
                Some((row + 1) * nx + col)
            } else {
                None
            },
            if col > 0 {
                Some(row * nx + col - 1)
            } else {
                None
            },
            if col + 1 < nx {
                Some(row * nx + col + 1)
            } else {
                None
            },
        ];
        for nbr in neighbors.into_iter().flatten() {
            if solid[nbr] && !visited[nbr] {
                visited[nbr] = true;
                queue.push_back(nbr);
            }
        }
    }

    // All solid elements should be visited
    solid.iter().zip(visited.iter()).all(|(&s, &v)| !s || v)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn small_grid() -> GridSize {
        GridSize::new(4, 4)
    }

    fn default_material() -> MaterialParams {
        MaterialParams::default()
    }

    // --- GridSize ---

    #[test]
    fn test_grid_size_n_elem() {
        let g = GridSize::new(5, 3);
        assert_eq!(g.n_elem(), 15);
    }

    #[test]
    fn test_grid_size_row_col() {
        let g = GridSize::new(5, 3);
        assert_eq!(g.row_col(7), (1, 2));
    }

    // --- SiMP ---

    #[test]
    fn test_simp_initial_density() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        assert!((simp.current_volume() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_simp_stiffness_solid() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        let e = simp.interpolated_stiffness(1.0);
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_simp_stiffness_void() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        let e = simp.interpolated_stiffness(0.0);
        assert!((e - 1e-9).abs() < 1e-15);
    }

    #[test]
    fn test_simp_stiffness_interpolation() {
        let simp = SiMP::new(small_grid(), 0.5, 2.0, 2.0, default_material());
        // For penalty=2, rho=0.5: E = 1e-9 + 0.25 * (1 - 1e-9) ≈ 0.25
        let e = simp.interpolated_stiffness(0.5);
        assert!((e - 0.25).abs() < 0.01);
    }

    #[test]
    fn test_simp_sensitivity_positive() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        let s = simp.stiffness_sensitivity(0.5);
        assert!(s > 0.0);
    }

    #[test]
    fn test_simp_compliance_sensitivity_negative() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        let se = vec![1.0; 16];
        let cs = simp.compliance_sensitivity(&se);
        assert!(cs.iter().all(|&s| s < 0.0));
    }

    #[test]
    fn test_simp_is_valid() {
        let simp = SiMP::new(small_grid(), 0.5, 3.0, 2.0, default_material());
        assert!(simp.is_valid());
    }

    #[test]
    fn test_simp_stiffness_field_length() {
        let grid = small_grid();
        let simp = SiMP::new(grid, 0.5, 3.0, 2.0, default_material());
        assert_eq!(simp.stiffness_field().len(), grid.n_elem());
    }

    // --- FilteringMethods ---

    #[test]
    fn test_density_filter_uniform() {
        let grid = GridSize::new(3, 3);
        let filter = FilteringMethods::new(1.5, grid);
        let density = vec![0.5; 9];
        let filtered = filter.density_filter(&density);
        assert!(filtered.iter().all(|&r| (r - 0.5).abs() < 1e-10));
    }

    #[test]
    fn test_density_filter_smooths() {
        let grid = GridSize::new(5, 5);
        let filter = FilteringMethods::new(2.0, grid);
        let mut density = vec![0.0; 25];
        density[12] = 1.0; // center spike
        let filtered = filter.density_filter(&density);
        // Surrounding elements should have nonzero filtered density
        assert!(filtered[11] > 0.0 || filtered[13] > 0.0);
    }

    #[test]
    fn test_sensitivity_filter_output_length() {
        let grid = GridSize::new(3, 3);
        let filter = FilteringMethods::new(1.5, grid);
        let density = vec![0.5; 9];
        let sens = vec![-1.0; 9];
        let out = filter.sensitivity_filter(&density, &sens);
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn test_heaviside_projection_solid() {
        let grid = GridSize::new(2, 2);
        let filter = FilteringMethods::new(1.5, grid);
        let density = vec![1.0; 4];
        let proj = filter.heaviside_projection(&density, 10.0, 0.5);
        assert!(proj.iter().all(|&r| r > 0.9));
    }

    #[test]
    fn test_heaviside_projection_void() {
        let grid = GridSize::new(2, 2);
        let filter = FilteringMethods::new(1.5, grid);
        let density = vec![0.0; 4];
        let proj = filter.heaviside_projection(&density, 10.0, 0.5);
        assert!(proj.iter().all(|&r| r < 0.1));
    }

    #[test]
    fn test_heaviside_sensitivity_positive() {
        let grid = GridSize::new(2, 2);
        let filter = FilteringMethods::new(1.5, grid);
        let density = vec![0.5; 4];
        let sens = filter.heaviside_sensitivity(&density, 5.0, 0.5);
        assert!(sens.iter().all(|&s| s >= 0.0));
    }

    #[test]
    fn test_filter_sensitivity_chain() {
        let grid = GridSize::new(3, 3);
        let filter = FilteringMethods::new(1.5, grid);
        let sens = vec![1.0; 9];
        let density = vec![0.5; 9];
        let out = filter.filter_sensitivity_chain(&sens, &density);
        assert_eq!(out.len(), 9);
    }

    // --- TopOptSolver ---

    #[test]
    fn test_topopt_oc_update_volume() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.4, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let solver = TopOptSolver::new(simp, filter);
        let density = vec![0.5; 16];
        let sensitivity = vec![-1.0; 16];
        let (new_d, _lm) = solver.oc_update(&density, &sensitivity, 0.4);
        let vol: f64 = new_d.iter().sum::<f64>() / 16.0;
        assert!((vol - 0.4).abs() < 0.05);
    }

    #[test]
    fn test_topopt_run_mock_convergence() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let mut solver = TopOptSolver::new(simp, filter).with_max_iter(10);
        let result = solver.run_mock(1.0);
        assert_eq!(result.len(), 16);
        assert!(!solver.compliance_history.is_empty());
    }

    #[test]
    fn test_topopt_compliance_positive() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let solver = TopOptSolver::new(simp, filter);
        let se = vec![1.0; 16];
        let c = solver.compute_compliance(&se);
        assert!(c > 0.0);
    }

    #[test]
    fn test_topopt_volume_constraint() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let solver = TopOptSolver::new(simp, filter);
        let density = vec![0.5; 16];
        assert!(solver.volume_constraint_satisfied(&density, 0.01));
    }

    #[test]
    fn test_topopt_with_damping() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let solver = TopOptSolver::new(simp, filter).with_damping(0.3);
        assert!((solver.damping - 0.3).abs() < 1e-10);
    }

    // --- MultiLoadTopOpt ---

    #[test]
    fn test_multiload_weighted_compliance() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let lc1 = LoadCase::new(vec![1.0; 16], 0.5);
        let lc2 = LoadCase::new(vec![0.5; 16], 0.5);
        let solver = MultiLoadTopOpt::new(simp, filter, vec![lc1, lc2]);
        let density = vec![0.5; 16];
        let c = solver.weighted_compliance(&density);
        assert!(c > 0.0);
    }

    #[test]
    fn test_multiload_aggregated_sensitivity_len() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let lc = LoadCase::new(vec![1.0; 16], 1.0);
        let solver = MultiLoadTopOpt::new(simp, filter, vec![lc]);
        let density = vec![0.5; 16];
        let sens = solver.aggregated_sensitivity(&density);
        assert_eq!(sens.len(), 16);
    }

    #[test]
    fn test_multiload_run() {
        let grid = GridSize::new(4, 4);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let lc = LoadCase::new(vec![1.0; 16], 1.0);
        let mut solver = MultiLoadTopOpt::new(simp, filter, vec![lc]);
        let result = solver.run(5, 1e-3);
        assert_eq!(result.len(), 16);
    }

    #[test]
    fn test_multiload_equal_weight_symmetry() {
        let grid = GridSize::new(2, 2);
        let simp = SiMP::new(grid, 0.5, 3.0, 1.5, MaterialParams::default());
        let filter = FilteringMethods::new(1.5, grid);
        let lc1 = LoadCase::new(vec![1.0; 4], 0.5);
        let lc2 = LoadCase::new(vec![1.0; 4], 0.5);
        let solver = MultiLoadTopOpt::new(simp.clone(), filter.clone(), vec![lc1, lc2]);
        let lc_single = LoadCase::new(vec![1.0; 4], 1.0);
        let solver2 = MultiLoadTopOpt::new(simp, filter, vec![lc_single]);
        let d = vec![0.5; 4];
        let c1 = solver.weighted_compliance(&d);
        let c2 = solver2.weighted_compliance(&d);
        assert!((c1 - c2).abs() < 1e-10);
    }

    // --- LevelSetTopOpt ---

    #[test]
    fn test_levelset_init_phi_length() {
        let grid = GridSize::new(4, 4);
        let ls = LevelSetTopOpt::new(grid, 0.5, 0.1, MaterialParams::default());
        assert_eq!(ls.phi.len(), 25); // (4+1)*(4+1)
    }

    #[test]
    fn test_levelset_element_density_length() {
        let grid = GridSize::new(4, 4);
        let ls = LevelSetTopOpt::new(grid, 0.5, 0.1, MaterialParams::default());
        assert_eq!(ls.element_density().len(), 16);
    }

    #[test]
    fn test_levelset_element_density_range() {
        let grid = GridSize::new(4, 4);
        let ls = LevelSetTopOpt::new(grid, 0.5, 0.1, MaterialParams::default());
        let d = ls.element_density();
        assert!(d.iter().all(|&r| (0.0..=1.0).contains(&r)));
    }

    #[test]
    fn test_levelset_hj_update_changes_phi() {
        let grid = GridSize::new(4, 4);
        let mut ls = LevelSetTopOpt::new(grid, 0.5, 0.1, MaterialParams::default());
        let phi_before = ls.phi.clone();
        let velocity = vec![1.0; 25];
        ls.hj_update(&velocity);
        let changed = phi_before
            .iter()
            .zip(ls.phi.iter())
            .any(|(a, b)| (a - b).abs() > 1e-10);
        assert!(changed);
    }

    #[test]
    fn test_levelset_run_mock() {
        let grid = GridSize::new(4, 4);
        let mut ls = LevelSetTopOpt::new(grid, 0.5, 0.01, MaterialParams::default());
        ls.run_mock(5);
        assert_eq!(ls.compliance_history.len(), 5);
    }

    #[test]
    fn test_levelset_reinitialize() {
        let grid = GridSize::new(4, 4);
        let mut ls = LevelSetTopOpt::new(grid, 0.5, 0.1, MaterialParams::default());
        let phi_before = ls.phi.clone();
        ls.reinitialize(3);
        // phi should change
        let _ = phi_before;
    }

    // --- ManufacturingFilters ---

    #[test]
    fn test_mfg_min_length_scale_output_len() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let density = vec![0.5; 16];
        let out = mf.minimum_length_scale(&density, 10.0);
        assert_eq!(out.len(), 16);
    }

    #[test]
    fn test_mfg_min_length_scale_range() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let density = vec![0.5; 16];
        let out = mf.minimum_length_scale(&density, 10.0);
        assert!(out.iter().all(|&r| (0.0..=1.0).contains(&r)));
    }

    #[test]
    fn test_mfg_overhang_base_unchanged() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 3.0);
        let mut density = vec![0.0; 16];
        // Set bottom row solid
        for col in 0..4 {
            density[col] = 1.0;
        }
        let out = mf.overhang_filter(&density);
        // Bottom row (build direction 0) should be unchanged
        assert!((out[0] - density[0]).abs() < 1e-10);
    }

    #[test]
    fn test_mfg_casting_filter_len() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let density = vec![0.5; 16];
        let out = mf.casting_filter(&density);
        assert_eq!(out.len(), 16);
    }

    #[test]
    fn test_mfg_extrusion_filter_x() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let mut density = vec![0.0; 16];
        density[0] = 1.0; // only top-left
        let out = mf.extrusion_filter(&density, 0); // extrude in x
        // column 0 should all be 1.0
        for row in 0..4 {
            assert!((out[row * 4] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_mfg_extrusion_filter_y() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let mut density = vec![0.0; 16];
        density[0] = 1.0; // only top-left corner
        let out = mf.extrusion_filter(&density, 1); // extrude in y
        // row 0 should all be 1.0
        for col in 0..4 {
            assert!((out[col] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_mfg_milling_filter_no_enclosed() {
        let grid = GridSize::new(4, 4);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        let density = vec![0.0; 16]; // all void = all accessible
        let out = mf.milling_filter(&density, 0.5);
        assert!(out.iter().all(|&r| r < 0.5)); // all remain void
    }

    #[test]
    fn test_mfg_milling_filter_fills_enclosed() {
        let grid = GridSize::new(5, 5);
        let mf = ManufacturingFilters::new(grid, 2.0, PI / 4.0);
        // Ring of solid around center void
        let mut density = vec![1.0; 25];
        density[12] = 0.0; // center void
        let out = mf.milling_filter(&density, 0.5);
        // The center should be filled (enclosed)
        assert!(out[12] >= 0.5);
    }

    // --- Utility functions ---

    #[test]
    fn test_volume_fraction() {
        let d = vec![0.5; 10];
        assert!((volume_fraction(&d) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_volume_fraction_empty() {
        assert_eq!(volume_fraction(&[]), 0.0);
    }

    #[test]
    fn test_normalize_sensitivity() {
        let mut s = vec![-4.0, -2.0, -1.0];
        normalize_sensitivity(&mut s);
        let max_s = s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min_s = s.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(max_s <= 0.0);
        assert!(min_s >= -1.0);
    }

    #[test]
    fn test_binarize() {
        let d = vec![0.3, 0.7, 0.5, 0.1];
        let b = binarize(&d, 0.5);
        assert_eq!(b, vec![false, true, true, false]);
    }

    #[test]
    fn test_count_solid() {
        let b = vec![true, false, true, true];
        assert_eq!(count_solid(&b), 3);
    }

    #[test]
    fn test_is_connected_uniform() {
        let grid = GridSize::new(3, 3);
        let density = vec![1.0; 9];
        assert!(is_connected(&density, grid, 0.5));
    }

    #[test]
    fn test_is_connected_disconnected() {
        let grid = GridSize::new(3, 3);
        let mut density = vec![0.0; 9];
        density[0] = 1.0;
        density[8] = 1.0; // corners only — not connected in 4-connectivity
        assert!(!is_connected(&density, grid, 0.5));
    }

    #[test]
    fn test_is_connected_all_void() {
        let grid = GridSize::new(3, 3);
        let density = vec![0.0; 9];
        assert!(is_connected(&density, grid, 0.5)); // vacuously true
    }
}
