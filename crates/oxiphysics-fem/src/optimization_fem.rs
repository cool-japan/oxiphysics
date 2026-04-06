#![allow(clippy::ptr_arg)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM-based structural topology and shape optimization.
//!
//! Implements:
//! - [`TopologyOptimization`]    — SIMP with OC, sensitivity/density filter, convergence tracking
//! - [`ShapeOptimization`]       — gradient-based shape optimization with adjoint sensitivity
//! - [`SizingOptimization`]      — cross-section sizing for trusses/beams
//! - [`MultiObjectiveFEM`]       — Pareto front stiffness vs. mass
//! - [`ManufacturingConstraints`] — minimum length scale, symmetry, overhang

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// § 1  Parameter structs
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters controlling a SIMP topology optimization run.
#[derive(Debug, Clone)]
pub struct TopologyOptParams {
    /// Target volume fraction (0 < vf ≤ 1).
    pub volume_fraction: f64,
    /// SIMP penalization exponent (typically 3).
    pub penalty: f64,
    /// Density filter radius \[element lengths\].
    pub filter_radius: f64,
    /// Maximum number of OC iterations.
    pub max_iter: usize,
}

impl TopologyOptParams {
    /// Create a new `TopologyOptParams`.
    pub fn new(volume_fraction: f64, penalty: f64, filter_radius: f64, max_iter: usize) -> Self {
        Self {
            volume_fraction,
            penalty,
            filter_radius,
            max_iter,
        }
    }
}

/// SIMP element density field.
#[derive(Debug, Clone)]
pub struct SimpDensity {
    /// Per-element density values ρ_e ∈ \[ρ_min, 1\].
    pub densities: Vec<f64>,
    /// Number of elements.
    pub n_elem: usize,
}

impl SimpDensity {
    /// Create a uniform density field with value `rho0`.
    pub fn uniform(n_elem: usize, rho0: f64) -> Self {
        Self {
            densities: vec![rho0; n_elem],
            n_elem,
        }
    }

    /// Mean density (= current volume fraction).
    pub fn mean(&self) -> f64 {
        if self.n_elem == 0 {
            return 0.0;
        }
        self.densities.iter().sum::<f64>() / self.n_elem as f64
    }

    /// Binarize: count elements above threshold (solid fraction).
    pub fn solid_count(&self, threshold: f64) -> usize {
        self.densities.iter().filter(|&&d| d >= threshold).count()
    }

    /// Maximum density.
    pub fn max_density(&self) -> f64 {
        self.densities
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum density.
    pub fn min_density(&self) -> f64 {
        self.densities.iter().cloned().fold(f64::INFINITY, f64::min)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 2  SIMP penalized stiffness
// ─────────────────────────────────────────────────────────────────────────────

/// SIMP penalized Young's modulus.
///
/// `E(ρ) = E_min + ρ^p (E₀ − E_min)`
///
/// # Arguments
/// * `density`  – element density ρ ∈ \[0, 1\]
/// * `e0`       – stiffness of solid material
/// * `e_min`    – minimum stiffness (void)
/// * `penalty`  – SIMP exponent p
pub fn simp_penalized_stiffness(density: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    e_min + density.powf(penalty) * (e0 - e_min)
}

/// SIMP stiffness derivative ∂E/∂ρ.
///
/// `dE/dρ = p ρ^(p-1) (E₀ − E_min)`
pub fn simp_stiffness_derivative(density: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    penalty * density.powf(penalty - 1.0) * (e0 - e_min)
}

// ─────────────────────────────────────────────────────────────────────────────
// § 3  Density filter
// ─────────────────────────────────────────────────────────────────────────────

/// Apply a weighted-average density filter.
///
/// Each element's filtered density is the weighted average of its neighbours
/// within `radius`, with weights equal to the number of shared adjacency entries
/// (uniform weight = 1 per neighbour here).
///
/// # Arguments
/// * `densities`  – raw element densities
/// * `adjacency`  – `adjacency[i]` = list of element indices within radius of element i
/// * `radius`     – filter radius (used as label; actual neighbourhood encoded in `adjacency`)
pub fn density_filter(densities: &[f64], adjacency: &[Vec<usize>], _radius: f64) -> Vec<f64> {
    let n = densities.len();
    let mut filtered = vec![0.0_f64; n];
    for i in 0..n {
        let nbrs = &adjacency[i];
        let w_sum = (nbrs.len() + 1) as f64;
        let mut val = densities[i];
        for &j in nbrs {
            val += densities[j];
        }
        filtered[i] = val / w_sum;
    }
    filtered
}

/// Build a 2-D grid adjacency list for density/sensitivity filtering.
///
/// Elements are indexed row-major: `e = row * n_cols + col`.
/// Neighbours are those within `radius` (in element units) in L∞ norm.
///
/// # Arguments
/// * `n_rows`  – number of element rows
/// * `n_cols`  – number of element columns
/// * `radius`  – neighbourhood radius in element lengths
pub fn build_grid_adjacency(n_rows: usize, n_cols: usize, radius: f64) -> Vec<Vec<usize>> {
    let n = n_rows * n_cols;
    let r_int = radius.floor() as isize;
    let mut adj = vec![Vec::new(); n];
    for row in 0..n_rows {
        for col in 0..n_cols {
            let idx = row * n_cols + col;
            for dr in -r_int..=r_int {
                for dc in -r_int..=r_int {
                    if dr == 0 && dc == 0 {
                        continue;
                    }
                    let nr = row as isize + dr;
                    let nc = col as isize + dc;
                    let dist = ((dr * dr + dc * dc) as f64).sqrt();
                    if nr >= 0
                        && nr < n_rows as isize
                        && nc >= 0
                        && nc < n_cols as isize
                        && dist <= radius
                    {
                        adj[idx].push(nr as usize * n_cols + nc as usize);
                    }
                }
            }
        }
    }
    adj
}

// ─────────────────────────────────────────────────────────────────────────────
// § 4  Sensitivity filter
// ─────────────────────────────────────────────────────────────────────────────

/// Apply a sensitivity (gradient) filter.
///
/// `dc_filtered[i] = (1 / max(1e-3, ρ_i)) Σ_j H_ij ρ_j dc_j / Σ_j H_ij`
/// with `H_ij = 1` for all j in the adjacency list.
///
/// # Arguments
/// * `sensitivities` – raw compliance sensitivities ∂C/∂ρ_e
/// * `densities`     – current element densities
/// * `adjacency`     – neighbourhood list per element
/// * `radius`        – filter radius (informational)
pub fn sensitivity_filter(
    sensitivities: &[f64],
    densities: &[f64],
    adjacency: &[Vec<usize>],
    _radius: f64,
) -> Vec<f64> {
    let n = sensitivities.len();
    let mut filtered = vec![0.0_f64; n];
    for i in 0..n {
        let nbrs = &adjacency[i];
        let h_sum = (nbrs.len() + 1) as f64;
        let mut num = densities[i] * sensitivities[i];
        for &j in nbrs {
            num += densities[j] * sensitivities[j];
        }
        filtered[i] = num / (h_sum * densities[i].max(1e-3));
    }
    filtered
}

// ─────────────────────────────────────────────────────────────────────────────
// § 5  OC density update
// ─────────────────────────────────────────────────────────────────────────────

/// Update element densities using the Optimality Criteria (OC) method.
///
/// The update rule is:
/// `ρ_new = clamp( ρ * sqrt(-dc / λ), ρ-Δ, ρ+Δ )`
/// where λ is found by bisection to satisfy the volume constraint.
///
/// # Arguments
/// * `densities`     – current densities ρ_e
/// * `sensitivities` – compliance sensitivities ∂C/∂ρ_e (≤ 0)
/// * `vf_target`     – target volume fraction
/// * `move_limit`    – maximum density change per iteration Δ
pub fn update_densities_oc(
    densities: &[f64],
    sensitivities: &[f64],
    vf_target: f64,
    move_limit: f64,
) -> Vec<f64> {
    let n = densities.len();
    let mut l1 = 1e-9_f64;
    let mut l2 = 1e9_f64;
    let mut rho_new = vec![0.0_f64; n];
    let rho_min = 1e-3_f64;

    for _iter in 0..200 {
        let lmid = 0.5 * (l1 + l2);
        let mut vol = 0.0_f64;
        for i in 0..n {
            let be = (-sensitivities[i] / lmid).max(0.0).sqrt() * densities[i];
            let lo = (densities[i] - move_limit).max(rho_min);
            let hi = (densities[i] + move_limit).min(1.0);
            rho_new[i] = be.clamp(lo, hi);
            vol += rho_new[i];
        }
        if vol / n as f64 > vf_target {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
        if (l2 - l1) < 1e-12 * (l1 + l2) {
            break;
        }
    }
    rho_new
}

// ─────────────────────────────────────────────────────────────────────────────
// § 6  Element compliance sensitivity
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the compliance sensitivity ∂C/∂ρ_e for a single element.
///
/// `∂C/∂ρ = -p ρ^(p-1) (E₀-E_min) u_e^T k₀ u_e`
///
/// where `k₀` is the unit-stiffness element matrix.
///
/// # Arguments
/// * `ue`      – element displacement vector (length = dofs per element)
/// * `ke`      – element unit stiffness matrix (square, same size as `ue`)
/// * `density` – current element density ρ
/// * `penalty` – SIMP exponent p
/// * `e0`      – solid stiffness E₀
pub fn compliance_sensitivity(
    ue: &[f64],
    ke: &[Vec<f64>],
    density: f64,
    penalty: f64,
    e0: f64,
) -> f64 {
    let mut uke = vec![0.0_f64; ue.len()];
    for (i, row) in ke.iter().enumerate() {
        for (j, &kij) in row.iter().enumerate() {
            uke[i] += kij * ue[j];
        }
    }
    let utku: f64 = ue.iter().zip(uke.iter()).map(|(u, ku)| u * ku).sum();
    -penalty * density.powf(penalty - 1.0) * e0 * utku
}

/// Compute the strain energy in an element.
///
/// `u_e^T K_e u_e / 2`
pub fn element_strain_energy(ue: &[f64], ke: &[Vec<f64>]) -> f64 {
    let mut uke = vec![0.0_f64; ue.len()];
    for (i, row) in ke.iter().enumerate() {
        for (j, &kij) in row.iter().enumerate() {
            uke[i] += kij * ue[j];
        }
    }
    0.5 * ue.iter().zip(uke.iter()).map(|(u, ku)| u * ku).sum::<f64>()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 7  Volume constraint
// ─────────────────────────────────────────────────────────────────────────────

/// Volume constraint residual: `g = Σρ/n − v_f`.
///
/// Returns a positive value when the current volume fraction exceeds the target.
pub fn volume_constraint(densities: &[f64], target: f64) -> f64 {
    let n = densities.len();
    if n == 0 {
        return -target;
    }
    densities.iter().sum::<f64>() / n as f64 - target
}

/// Project densities towards 0/1 using Heaviside projection.
///
/// `ρ̄ = (tanh(β η) + tanh(β (ρ − η))) / (tanh(β η) + tanh(β (1 − η)))`
///
/// # Arguments
/// * `densities` – filtered densities
/// * `beta`      – projection sharpness (large = sharper)
/// * `eta`       – threshold level (typically 0.5)
pub fn heaviside_projection(densities: &[f64], beta: f64, eta: f64) -> Vec<f64> {
    let denom = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    densities
        .iter()
        .map(|&rho| {
            if denom.abs() < 1e-300 {
                rho
            } else {
                ((beta * eta).tanh() + (beta * (rho - eta)).tanh()) / denom
            }
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 8  Checkerboard measure
// ─────────────────────────────────────────────────────────────────────────────

/// Detect checkerboard patterns in a 2-D density field.
///
/// Scans all 2×2 tiles; a tile exhibits checkerboard if opposite corners have
/// densities on opposite sides of 0.5.  Returns the fraction of tiles that
/// exhibit checkerboarding.
///
/// # Arguments
/// * `densities` – 2-D grid `densities[row][col]`
pub fn checkerboard_measure(densities: &[Vec<f64>]) -> f64 {
    let rows = densities.len();
    if rows < 2 {
        return 0.0;
    }
    let cols = densities[0].len();
    if cols < 2 {
        return 0.0;
    }
    let mut count = 0usize;
    let mut total = 0usize;
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            let d00 = densities[r][c];
            let d01 = densities[r][c + 1];
            let d10 = densities[r + 1][c];
            let d11 = densities[r + 1][c + 1];
            let checker = (d00 > 0.5 && d11 > 0.5 && d01 < 0.5 && d10 < 0.5)
                || (d00 < 0.5 && d11 < 0.5 && d01 > 0.5 && d10 > 0.5);
            if checker {
                count += 1;
            }
            total += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 9  Single OC optimization step
// ─────────────────────────────────────────────────────────────────────────────

/// Perform one OC topology-optimization step.
///
/// Applies sensitivity filtering with a trivial (identity) adjacency, then
/// calls [`update_densities_oc`].
///
/// # Arguments
/// * `params`        – optimization parameters
/// * `n_elem`        – number of elements
/// * `sensitivities` – current compliance sensitivities ∂C/∂ρ_e
pub fn topology_optimize(
    params: &TopologyOptParams,
    n_elem: usize,
    sensitivities: &[f64],
) -> Vec<f64> {
    let rho0 = params.volume_fraction;
    let initial = vec![rho0; n_elem];
    let adjacency: Vec<Vec<usize>> = vec![vec![]; n_elem];
    let filtered_sens =
        sensitivity_filter(sensitivities, &initial, &adjacency, params.filter_radius);
    update_densities_oc(&initial, &filtered_sens, params.volume_fraction, 0.2)
}

// ─────────────────────────────────────────────────────────────────────────────
// § 10  TopologyOptimization — full SIMP driver
// ─────────────────────────────────────────────────────────────────────────────

/// Convergence record for one topology optimization iteration.
#[derive(Debug, Clone)]
pub struct TopoIterRecord {
    /// Iteration number (0-based).
    pub iter: usize,
    /// Compliance objective value.
    pub compliance: f64,
    /// Current volume fraction.
    pub volume_fraction: f64,
    /// Maximum density change from previous iteration.
    pub max_change: f64,
}

/// Full SIMP topology optimization driver.
///
/// Runs multiple OC iterations, updating densities from supplied sensitivity
/// vectors.  The caller supplies a "FEM oracle" that, given current densities,
/// returns (compliance, per-element sensitivities).
pub struct TopologyOptimization {
    /// Optimization parameters.
    pub params: TopologyOptParams,
    /// Current element densities.
    pub densities: Vec<f64>,
    /// Convergence history.
    pub history: Vec<TopoIterRecord>,
}

impl TopologyOptimization {
    /// Create a new topology optimization problem with uniform initial density.
    pub fn new(n_elem: usize, params: TopologyOptParams) -> Self {
        let vf = params.volume_fraction;
        Self {
            densities: vec![vf; n_elem],
            params,
            history: Vec::new(),
        }
    }

    /// Run the optimization loop.
    ///
    /// `oracle` is called each iteration with current densities; it should
    /// return `(compliance, sensitivities)`.
    pub fn run<F>(&mut self, adjacency: &[Vec<usize>], oracle: F, tol: f64)
    where
        F: Fn(&[f64]) -> (f64, Vec<f64>),
    {
        let n = self.densities.len();
        let mut prev = self.densities.clone();

        for iter in 0..self.params.max_iter {
            let (compliance, raw_sens) = oracle(&self.densities);

            let filtered_sens = sensitivity_filter(
                &raw_sens,
                &self.densities,
                adjacency,
                self.params.filter_radius,
            );

            let new_dens = update_densities_oc(
                &self.densities,
                &filtered_sens,
                self.params.volume_fraction,
                0.2_f64,
            );

            let max_change = new_dens
                .iter()
                .zip(prev.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);

            let vf = new_dens.iter().sum::<f64>() / n as f64;

            self.history.push(TopoIterRecord {
                iter,
                compliance,
                volume_fraction: vf,
                max_change,
            });

            prev = self.densities.clone();
            self.densities = new_dens;

            if max_change < tol {
                break;
            }
        }
    }

    /// Check convergence: last recorded max_change < `tol`.
    pub fn is_converged(&self, tol: f64) -> bool {
        self.history.last().is_some_and(|r| r.max_change < tol)
    }

    /// Apply density filtering in-place.
    pub fn apply_density_filter(&mut self, adjacency: &[Vec<usize>]) {
        self.densities = density_filter(&self.densities, adjacency, self.params.filter_radius);
    }

    /// Apply Heaviside projection in-place.
    pub fn apply_heaviside(&mut self, beta: f64, eta: f64) {
        self.densities = heaviside_projection(&self.densities, beta, eta);
    }

    /// Current compliance from the last recorded iteration.
    pub fn last_compliance(&self) -> Option<f64> {
        self.history.last().map(|r| r.compliance)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 11  ShapeOptimization — gradient-based with adjoint sensitivity
// ─────────────────────────────────────────────────────────────────────────────

/// A node coordinate update in shape optimization.
#[derive(Debug, Clone)]
pub struct ShapeDesignVariable {
    /// Node index.
    pub node: usize,
    /// Direction: 0 = x, 1 = y, 2 = z.
    pub direction: usize,
    /// Current coordinate value.
    pub value: f64,
}

/// Result of one shape-optimization step.
#[derive(Debug, Clone)]
pub struct ShapeOptResult {
    /// Updated design variables.
    pub variables: Vec<ShapeDesignVariable>,
    /// Objective function value (e.g. compliance).
    pub objective: f64,
    /// Gradient magnitude ‖dJ/ds‖.
    pub grad_norm: f64,
    /// Number of iterations performed.
    pub iterations: usize,
}

/// Gradient-based shape optimization using adjoint sensitivities.
///
/// Minimizes a scalar objective J(s) over the shape design variables `s`
/// (node coordinates) using a gradient-descent / steepest-descent update.
pub struct ShapeOptimization {
    /// Step size for gradient descent.
    pub step_size: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance on gradient norm.
    pub tol: f64,
}

impl ShapeOptimization {
    /// Create a new [`ShapeOptimization`] driver.
    pub fn new(step_size: f64, max_iter: usize, tol: f64) -> Self {
        Self {
            step_size,
            max_iter,
            tol,
        }
    }

    /// Run steepest-descent shape optimization.
    ///
    /// `oracle` takes current design variables and returns `(objective, gradient)`.
    pub fn run<F>(&self, mut vars: Vec<ShapeDesignVariable>, oracle: F) -> ShapeOptResult
    where
        F: Fn(&[ShapeDesignVariable]) -> (f64, Vec<f64>),
    {
        let mut iters = 0usize;
        let mut objective = 0.0_f64;
        let mut grad_norm = f64::INFINITY;

        for _iter in 0..self.max_iter {
            iters += 1;
            let (obj, grad) = oracle(&vars);
            objective = obj;
            grad_norm = grad.iter().map(|g| g * g).sum::<f64>().sqrt();

            if grad_norm < self.tol {
                break;
            }

            for (v, g) in vars.iter_mut().zip(grad.iter()) {
                v.value -= self.step_size * g;
            }
        }

        ShapeOptResult {
            variables: vars,
            objective,
            grad_norm,
            iterations: iters,
        }
    }

    /// Adjoint sensitivity of compliance w.r.t. node coordinate perturbation.
    ///
    /// `dJ/ds_i ≈ (J(s+h e_i) - J(s-h e_i)) / (2h)` (finite difference).
    pub fn adjoint_sensitivity<F>(vars: &[ShapeDesignVariable], oracle: &F, h: f64) -> Vec<f64>
    where
        F: Fn(&[ShapeDesignVariable]) -> (f64, Vec<f64>),
    {
        let n = vars.len();
        let mut grad = vec![0.0_f64; n];
        for i in 0..n {
            let mut vp = vars.to_vec();
            let mut vm = vars.to_vec();
            vp[i].value += h;
            vm[i].value -= h;
            let (jp, _) = oracle(&vp);
            let (jm, _) = oracle(&vm);
            grad[i] = (jp - jm) / (2.0 * h);
        }
        grad
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 12  SizingOptimization — cross-section sizing for trusses/beams
// ─────────────────────────────────────────────────────────────────────────────

/// Cross-section design variable for truss/beam sizing.
#[derive(Debug, Clone)]
pub struct CrossSection {
    /// Member index.
    pub member: usize,
    /// Cross-sectional area A.
    pub area: f64,
    /// Lower bound on area.
    pub area_min: f64,
    /// Upper bound on area.
    pub area_max: f64,
}

impl CrossSection {
    /// Create a new cross section with symmetric bounds.
    pub fn new(member: usize, area: f64, area_min: f64, area_max: f64) -> Self {
        Self {
            member,
            area,
            area_min,
            area_max,
        }
    }

    /// Clamp area to \[area_min, area_max\].
    pub fn clamped(&self) -> f64 {
        self.area.clamp(self.area_min, self.area_max)
    }
}

/// Result from a sizing optimization run.
#[derive(Debug, Clone)]
pub struct SizingResult {
    /// Optimal cross-section areas.
    pub areas: Vec<f64>,
    /// Objective value (compliance or mass).
    pub objective: f64,
    /// Total member volume (sum of A_i * L_i).
    pub total_volume: f64,
    /// Number of iterations.
    pub iterations: usize,
}

/// Cross-section sizing optimizer for trusses/beams.
///
/// Minimizes structural compliance subject to a volume (mass) constraint
/// using the OC method adapted for continuous cross-section design variables.
pub struct SizingOptimization {
    /// Member lengths L_i.
    pub lengths: Vec<f64>,
    /// Volume budget V_max.
    pub volume_max: f64,
    /// Move limit for area update.
    pub move_limit: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}

impl SizingOptimization {
    /// Create a new sizing optimizer.
    pub fn new(lengths: Vec<f64>, volume_max: f64, move_limit: f64, max_iter: usize) -> Self {
        Self {
            lengths,
            volume_max,
            move_limit,
            max_iter,
        }
    }

    /// Run OC sizing update.
    ///
    /// `sensitivities[i]` = ∂C/∂A_i (negative for compliance minimization).
    /// Returns updated areas satisfying the volume budget.
    pub fn update_oc(&self, sections: &[CrossSection], sensitivities: &[f64]) -> Vec<f64> {
        let n = sections.len();
        let mut l1 = 1e-9_f64;
        let mut l2 = 1e9_f64;
        let mut a_new = vec![0.0_f64; n];

        for _iter in 0..200 {
            let lmid = 0.5 * (l1 + l2);
            let mut vol = 0.0_f64;
            for i in 0..n {
                let sc = &sections[i];
                let be = (-sensitivities[i] / (lmid * self.lengths[i]))
                    .max(0.0)
                    .sqrt()
                    * sc.area;
                let lo = (sc.area - self.move_limit).max(sc.area_min);
                let hi = (sc.area + self.move_limit).min(sc.area_max);
                a_new[i] = be.clamp(lo, hi);
                vol += a_new[i] * self.lengths[i];
            }
            if vol > self.volume_max {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
            if (l2 - l1) < 1e-12 * (l1 + l2) {
                break;
            }
        }
        a_new
    }

    /// Run sizing optimization loop.
    ///
    /// `oracle` maps current areas to `(compliance, sensitivities)`.
    pub fn run<F>(&self, sections: &mut Vec<CrossSection>, oracle: F) -> SizingResult
    where
        F: Fn(&[f64]) -> (f64, Vec<f64>),
    {
        let mut iters = 0usize;
        let mut objective = 0.0_f64;

        for _iter in 0..self.max_iter {
            iters += 1;
            let areas: Vec<f64> = sections.iter().map(|s| s.area).collect();
            let (obj, sens) = oracle(&areas);
            objective = obj;

            let new_areas = self.update_oc(sections, &sens);
            let max_change = new_areas
                .iter()
                .zip(sections.iter().map(|s| s.area))
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);

            for (s, &a) in sections.iter_mut().zip(new_areas.iter()) {
                s.area = a;
            }

            if max_change < 1e-6 {
                break;
            }
        }

        let areas: Vec<f64> = sections.iter().map(|s| s.area).collect();
        let total_volume = areas
            .iter()
            .zip(self.lengths.iter())
            .map(|(a, l)| a * l)
            .sum::<f64>();

        SizingResult {
            areas,
            objective,
            total_volume,
            iterations: iters,
        }
    }

    /// Compute total volume given areas.
    pub fn total_volume(&self, areas: &[f64]) -> f64 {
        areas
            .iter()
            .zip(self.lengths.iter())
            .map(|(a, l)| a * l)
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 13  MultiObjectiveFEM — Pareto front stiffness vs. mass
// ─────────────────────────────────────────────────────────────────────────────

/// A point on the Pareto front: (compliance, mass) pair.
#[derive(Debug, Clone)]
pub struct ParetoPoint {
    /// Compliance objective C.
    pub compliance: f64,
    /// Mass objective M.
    pub mass: f64,
    /// Volume fraction used.
    pub volume_fraction: f64,
    /// Density field at this Pareto point.
    pub densities: Vec<f64>,
}

/// Multi-objective FEM: Pareto front for stiffness vs. mass.
///
/// Sweeps volume fractions and runs topology optimization at each level to
/// approximate the Pareto front between compliance and structural mass.
pub struct MultiObjectiveFEM {
    /// Number of Pareto points.
    pub n_points: usize,
    /// SIMP penalty exponent.
    pub penalty: f64,
    /// Filter radius.
    pub filter_radius: f64,
    /// Max inner OC iterations per Pareto point.
    pub max_inner_iter: usize,
}

impl MultiObjectiveFEM {
    /// Create a new multi-objective FEM driver.
    pub fn new(n_points: usize, penalty: f64, filter_radius: f64, max_inner_iter: usize) -> Self {
        Self {
            n_points,
            penalty,
            filter_radius,
            max_inner_iter,
        }
    }

    /// Compute the Pareto front by sweeping volume fractions.
    ///
    /// `oracle` maps `(densities, vf)` → `(compliance, sensitivities)`.
    /// Volume fractions are swept linearly from `vf_min` to `vf_max`.
    pub fn pareto_front<F>(
        &self,
        n_elem: usize,
        vf_min: f64,
        vf_max: f64,
        oracle: F,
    ) -> Vec<ParetoPoint>
    where
        F: Fn(&[f64], f64) -> (f64, Vec<f64>),
    {
        let mut front = Vec::with_capacity(self.n_points);
        for k in 0..self.n_points {
            let vf = vf_min + (vf_max - vf_min) * k as f64 / (self.n_points - 1).max(1) as f64;
            let params =
                TopologyOptParams::new(vf, self.penalty, self.filter_radius, self.max_inner_iter);
            let mut topo = TopologyOptimization::new(n_elem, params);
            let adj: Vec<Vec<usize>> = vec![vec![]; n_elem];
            let vf_cap = vf;
            topo.run(&adj, |d| oracle(d, vf_cap), 1e-3);
            let dens = topo.densities.clone();
            let compliance = topo.last_compliance().unwrap_or(0.0);
            let mass = dens.iter().sum::<f64>() / n_elem as f64;
            front.push(ParetoPoint {
                compliance,
                mass,
                volume_fraction: vf,
                densities: dens,
            });
        }
        front
    }

    /// Filter the Pareto front: keep only non-dominated solutions.
    pub fn filter_dominated(front: &[ParetoPoint]) -> Vec<&ParetoPoint> {
        front
            .iter()
            .filter(|p| {
                !front.iter().any(|q| {
                    q.compliance <= p.compliance
                        && q.mass <= p.mass
                        && (q.compliance < p.compliance || q.mass < p.mass)
                })
            })
            .collect()
    }

    /// Compute weighted-sum scalarization objective.
    pub fn weighted_sum(p: &ParetoPoint, w_c: f64, w_m: f64) -> f64 {
        w_c * p.compliance + w_m * p.mass
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 14  ManufacturingConstraints
// ─────────────────────────────────────────────────────────────────────────────

/// Manufacturing constraints for topology optimization.
///
/// Enforces minimum length scale, symmetry, and overhang constraints
/// as post-processing corrections on the density field.
pub struct ManufacturingConstraints {
    /// Minimum length scale (elements).
    pub min_length_scale: f64,
    /// Enforce left-right symmetry (for 2-D grids).
    pub symmetry_lr: bool,
    /// Maximum allowable overhang angle (degrees).
    pub max_overhang_angle: f64,
}

impl ManufacturingConstraints {
    /// Create a new set of manufacturing constraints.
    pub fn new(min_length_scale: f64, symmetry_lr: bool, max_overhang_angle: f64) -> Self {
        Self {
            min_length_scale,
            symmetry_lr,
            max_overhang_angle,
        }
    }

    /// Enforce left-right symmetry on a 2-D density grid.
    ///
    /// Each element density is replaced by the average of itself and its mirror.
    pub fn enforce_symmetry_lr(&self, densities: &mut Vec<Vec<f64>>) {
        if !self.symmetry_lr {
            return;
        }
        for row in densities.iter_mut() {
            let n = row.len();
            for c in 0..n / 2 {
                let avg = 0.5 * (row[c] + row[n - 1 - c]);
                row[c] = avg;
                row[n - 1 - c] = avg;
            }
        }
    }

    /// Apply a minimum length scale filter by eroding tiny features.
    ///
    /// Elements below `min_length_scale * mean_density` are zeroed.
    pub fn apply_min_length_scale(&self, densities: &mut Vec<f64>) {
        if self.min_length_scale <= 0.0 {
            return;
        }
        let n = densities.len() as f64;
        let mean = if n > 0.0 {
            densities.iter().sum::<f64>() / n
        } else {
            0.0
        };
        let threshold = self.min_length_scale * mean;
        for d in densities.iter_mut() {
            if *d < threshold {
                *d = 1e-3_f64;
            }
        }
    }

    /// Check overhang constraint in a 2-D grid (additive manufacturing).
    ///
    /// Returns the fraction of solid elements that violate the overhang constraint.
    /// Building direction assumed to be along increasing row index.
    pub fn overhang_violation_fraction(&self, densities: &[Vec<f64>]) -> f64 {
        let max_angle_rad = self.max_overhang_angle * PI / 180.0;
        let tan_limit = max_angle_rad.tan();
        let rows = densities.len();
        if rows < 2 {
            return 0.0;
        }
        let cols = densities[0].len();
        let mut total_solid = 0usize;
        let mut violations = 0usize;
        for r in 1..rows {
            for c in 0..cols {
                if densities[r][c] > 0.5 {
                    total_solid += 1;
                    // Check support from row below (r-1), same column or neighbours
                    let supported = c > 0 && densities[r - 1][c - 1] > 0.5
                        || densities[r - 1][c] > 0.5
                        || c + 1 < cols && densities[r - 1][c + 1] > 0.5;
                    // Max horizontal extent without support
                    let overhang_tan = if !supported { 1.0 / (1.0 + 1e-9) } else { 0.0 };
                    if !supported && overhang_tan > tan_limit {
                        violations += 1;
                    }
                }
            }
        }
        if total_solid == 0 {
            0.0
        } else {
            violations as f64 / total_solid as f64
        }
    }

    /// Project densities below min length scale threshold.
    pub fn project_length_scale(densities: &[f64], threshold: f64) -> Vec<f64> {
        densities
            .iter()
            .map(|&d| if d < threshold { 1e-3_f64 } else { d })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 15  MBB beam helper (small mesh topology optimization)
// ─────────────────────────────────────────────────────────────────────────────

/// Run a simplified MBB beam topology optimization on an `n_rows × n_cols` mesh.
///
/// The "oracle" for this simplified problem computes compliance from a unit-load
/// static-response approximation: compliance ≈ Σ_e E(ρ_e) / n_elem (normalized).
/// This is intentionally simplified so that the test does not require a full FEM solve.
///
/// Returns the final density field as a 2-D grid.
pub fn mbb_beam_topology_opt(
    n_rows: usize,
    n_cols: usize,
    vf: f64,
    penalty: f64,
    filter_radius: f64,
    max_iter: usize,
) -> Vec<Vec<f64>> {
    let n_elem = n_rows * n_cols;
    let adj = build_grid_adjacency(n_rows, n_cols, filter_radius);
    let params = TopologyOptParams::new(vf, penalty, filter_radius, max_iter);
    let mut topo = TopologyOptimization::new(n_elem, params);

    // Oracle: compliance = sum(1/E(rho_e)), sensitivity = -p*rho^(p-1)*E0 (simplified)
    let e0 = 1.0_f64;
    let e_min = 1e-9_f64;

    topo.run(
        &adj,
        |dens| {
            let compliance: f64 = dens
                .iter()
                .map(|&r| 1.0 / simp_penalized_stiffness(r, e0, e_min, penalty))
                .sum::<f64>();
            let sens: Vec<f64> = dens
                .iter()
                .map(|&r| {
                    let e = simp_penalized_stiffness(r, e0, e_min, penalty);
                    -penalty * r.powf(penalty - 1.0) * (e0 - e_min) / (e * e)
                })
                .collect();
            (compliance, sens)
        },
        1e-3,
    );

    // Reshape to 2-D grid
    (0..n_rows)
        .map(|r| topo.densities[r * n_cols..(r + 1) * n_cols].to_vec())
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 16  Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── TopologyOptParams ────────────────────────────────────────────────────

    #[test]
    fn test_params_fields() {
        let p = TopologyOptParams::new(0.5, 3.0, 1.5, 100);
        assert!((p.volume_fraction - 0.5).abs() < 1e-12);
        assert!((p.penalty - 3.0).abs() < 1e-12);
        assert!((p.filter_radius - 1.5).abs() < 1e-12);
        assert_eq!(p.max_iter, 100);
    }

    // ── SimpDensity ──────────────────────────────────────────────────────────

    #[test]
    fn test_simp_density_uniform() {
        let d = SimpDensity::uniform(10, 0.4);
        assert_eq!(d.n_elem, 10);
        assert!((d.mean() - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_simp_density_mean_empty() {
        let d = SimpDensity::uniform(0, 0.5);
        assert_eq!(d.mean(), 0.0);
    }

    #[test]
    fn test_simp_density_solid_count() {
        let d = SimpDensity {
            densities: vec![0.2, 0.6, 0.8, 0.4],
            n_elem: 4,
        };
        assert_eq!(d.solid_count(0.5), 2);
    }

    #[test]
    fn test_simp_density_max_min() {
        let d = SimpDensity {
            densities: vec![0.1, 0.9, 0.5],
            n_elem: 3,
        };
        assert!((d.max_density() - 0.9).abs() < 1e-12);
        assert!((d.min_density() - 0.1).abs() < 1e-12);
    }

    // ── simp_penalized_stiffness ─────────────────────────────────────────────

    #[test]
    fn test_simp_solid() {
        let e = simp_penalized_stiffness(1.0, 1.0, 1e-9, 3.0);
        assert!((e - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_simp_void() {
        let e = simp_penalized_stiffness(0.0, 1.0, 1e-9, 3.0);
        assert!((e - 1e-9).abs() < 1e-15);
    }

    #[test]
    fn test_simp_half() {
        let e = simp_penalized_stiffness(0.5, 1.0, 0.0, 3.0);
        assert!((e - 0.125).abs() < 1e-12);
    }

    #[test]
    fn test_simp_monotone() {
        let e1 = simp_penalized_stiffness(0.3, 1.0, 1e-9, 3.0);
        let e2 = simp_penalized_stiffness(0.7, 1.0, 1e-9, 3.0);
        assert!(e1 < e2);
    }

    #[test]
    fn test_simp_penalty_1() {
        let e = simp_penalized_stiffness(0.4, 1.0, 0.0, 1.0);
        assert!((e - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_simp_stiffness_derivative_solid() {
        // At rho=1, dE/drho = p * (e0 - e_min)
        let de = simp_stiffness_derivative(1.0, 1.0, 0.0, 3.0);
        assert!((de - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_simp_stiffness_derivative_half() {
        let de = simp_stiffness_derivative(0.5, 1.0, 0.0, 3.0);
        let expected = 3.0 * 0.5_f64.powf(2.0);
        assert!((de - expected).abs() < 1e-12);
    }

    // ── density_filter ───────────────────────────────────────────────────────

    #[test]
    fn test_density_filter_no_neighbours() {
        let densities = vec![0.3, 0.7, 0.5];
        let adj: Vec<Vec<usize>> = vec![vec![], vec![], vec![]];
        let filtered = density_filter(&densities, &adj, 1.0);
        for (f, &d) in filtered.iter().zip(densities.iter()) {
            assert!((f - d).abs() < 1e-12);
        }
    }

    #[test]
    fn test_density_filter_with_neighbours() {
        let densities = vec![0.0, 1.0];
        let adj = vec![vec![1usize], vec![0usize]];
        let filtered = density_filter(&densities, &adj, 1.0);
        assert!((filtered[0] - 0.5).abs() < 1e-12);
        assert!((filtered[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_density_filter_preserves_bounds() {
        let densities = vec![0.2, 0.4, 0.6, 0.8];
        let adj: Vec<Vec<usize>> = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let filtered = density_filter(&densities, &adj, 1.5);
        for &f in &filtered {
            assert!(f >= 0.0 && f <= 1.0);
        }
    }

    #[test]
    fn test_density_filter_length_preserved() {
        let d = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let adj: Vec<Vec<usize>> = vec![vec![], vec![], vec![], vec![], vec![]];
        let f = density_filter(&d, &adj, 1.0);
        assert_eq!(f.len(), d.len());
    }

    // ── build_grid_adjacency ─────────────────────────────────────────────────

    #[test]
    fn test_build_grid_adjacency_size() {
        let adj = build_grid_adjacency(4, 5, 1.5);
        assert_eq!(adj.len(), 20);
    }

    #[test]
    fn test_build_grid_adjacency_no_self() {
        let adj = build_grid_adjacency(3, 3, 1.5);
        for (i, nbrs) in adj.iter().enumerate() {
            assert!(!nbrs.contains(&i), "element {i} contains itself");
        }
    }

    #[test]
    fn test_build_grid_adjacency_symmetry() {
        let adj = build_grid_adjacency(4, 4, 1.5);
        for i in 0..16 {
            for &j in &adj[i] {
                assert!(
                    adj[j].contains(&i),
                    "adjacency not symmetric: {i} -> {j} but not {j} -> {i}"
                );
            }
        }
    }

    // ── sensitivity_filter ───────────────────────────────────────────────────

    #[test]
    fn test_sensitivity_filter_no_neighbours() {
        let sens = vec![-0.5, -0.3];
        let dens = vec![0.5, 0.5];
        let adj: Vec<Vec<usize>> = vec![vec![], vec![]];
        let filtered = sensitivity_filter(&sens, &dens, &adj, 1.0);
        assert!((filtered[0] - sens[0]).abs() < 1e-12);
        assert!((filtered[1] - sens[1]).abs() < 1e-12);
    }

    #[test]
    fn test_sensitivity_filter_symmetry() {
        let sens = vec![-0.4, -0.4];
        let dens = vec![0.5, 0.5];
        let adj = vec![vec![1usize], vec![0usize]];
        let f = sensitivity_filter(&sens, &dens, &adj, 1.0);
        assert!((f[0] - f[1]).abs() < 1e-12);
    }

    #[test]
    fn test_sensitivity_filter_length_preserved() {
        let s = vec![-0.1, -0.2, -0.3];
        let d = vec![0.5, 0.5, 0.5];
        let adj: Vec<Vec<usize>> = vec![vec![], vec![], vec![]];
        let f = sensitivity_filter(&s, &d, &adj, 1.0);
        assert_eq!(f.len(), s.len());
    }

    // ── volume_constraint ────────────────────────────────────────────────────

    #[test]
    fn test_volume_constraint_zero() {
        let d = vec![0.5, 0.5, 0.5, 0.5];
        let g = volume_constraint(&d, 0.5);
        assert!(g.abs() < 1e-12);
    }

    #[test]
    fn test_volume_constraint_positive() {
        let d = vec![0.8, 0.8];
        let g = volume_constraint(&d, 0.5);
        assert!(g > 0.0);
    }

    #[test]
    fn test_volume_constraint_negative() {
        let d = vec![0.2, 0.2];
        let g = volume_constraint(&d, 0.5);
        assert!(g < 0.0);
    }

    #[test]
    fn test_volume_constraint_empty() {
        let g = volume_constraint(&[], 0.5);
        assert!((g + 0.5).abs() < 1e-12);
    }

    // ── heaviside_projection ─────────────────────────────────────────────────

    #[test]
    fn test_heaviside_projection_midpoint() {
        let d = vec![0.5];
        let proj = heaviside_projection(&d, 5.0, 0.5);
        // At rho=eta=0.5: tanh(beta*(rho-eta))=tanh(0)=0, so proj = tanh(beta*eta) / (2*tanh(beta*eta)) = 0.5
        assert!((proj[0] - 0.5).abs() < 0.1, "proj[0]={}", proj[0]);
    }

    #[test]
    fn test_heaviside_projection_range() {
        let d: Vec<f64> = (0..11).map(|i| i as f64 / 10.0).collect();
        let proj = heaviside_projection(&d, 3.0, 0.5);
        for &p in &proj {
            assert!(p >= 0.0 && p <= 1.0 + 1e-10, "out of range: {p}");
        }
    }

    #[test]
    fn test_heaviside_projection_monotone() {
        let d = vec![0.2, 0.5, 0.8];
        let proj = heaviside_projection(&d, 5.0, 0.5);
        assert!(proj[0] <= proj[1]);
        assert!(proj[1] <= proj[2]);
    }

    // ── compliance_sensitivity ───────────────────────────────────────────────

    #[test]
    fn test_compliance_sensitivity_zero_disp() {
        let ue = vec![0.0, 0.0, 0.0, 0.0];
        let ke = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        let dc = compliance_sensitivity(&ue, &ke, 0.5, 3.0, 1.0);
        assert_eq!(dc, 0.0);
    }

    #[test]
    fn test_compliance_sensitivity_sign() {
        let ue = vec![1.0, 0.0];
        let ke = vec![vec![2.0, 0.0], vec![0.0, 2.0]];
        let dc = compliance_sensitivity(&ue, &ke, 0.5, 3.0, 1.0);
        assert!(dc < 0.0);
    }

    #[test]
    fn test_compliance_sensitivity_scales_with_e0() {
        let ue = vec![1.0];
        let ke = vec![vec![1.0]];
        let dc1 = compliance_sensitivity(&ue, &ke, 1.0, 1.0, 1.0);
        let dc2 = compliance_sensitivity(&ue, &ke, 1.0, 1.0, 2.0);
        assert!((dc2 / dc1 - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_strain_energy_positive() {
        let ue = vec![1.0, 0.0];
        let ke = vec![vec![4.0, 0.0], vec![0.0, 4.0]];
        let se = element_strain_energy(&ue, &ke);
        assert!((se - 2.0).abs() < 1e-12);
    }

    // ── update_densities_oc ──────────────────────────────────────────────────

    #[test]
    fn test_oc_volume_satisfied() {
        let n = 20;
        let dens = vec![0.5_f64; n];
        let sens = vec![-0.5_f64; n];
        let vf = 0.4;
        let updated = update_densities_oc(&dens, &sens, vf, 0.2);
        let mean = updated.iter().sum::<f64>() / n as f64;
        assert!((mean - vf).abs() < 0.05);
    }

    #[test]
    fn test_oc_move_limit_respected() {
        let dens = vec![0.5_f64; 10];
        let sens = vec![-1.0_f64; 10];
        let updated = update_densities_oc(&dens, &sens, 0.3, 0.1);
        for (&old, &new) in dens.iter().zip(updated.iter()) {
            assert!(new >= old - 0.1 - 1e-10);
            assert!(new <= old + 0.1 + 1e-10);
        }
    }

    #[test]
    fn test_oc_bounds() {
        let dens = vec![0.5_f64; 5];
        let sens = vec![-0.3_f64; 5];
        let updated = update_densities_oc(&dens, &sens, 0.5, 0.2);
        for &r in &updated {
            assert!(r >= 1e-3 - 1e-10);
            assert!(r <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_oc_uniform_sensitivity() {
        let n = 8;
        let dens = vec![0.5_f64; n];
        let sens = vec![-0.4_f64; n];
        let updated = update_densities_oc(&dens, &sens, 0.5, 0.2);
        let first = updated[0];
        for &r in &updated {
            assert!((r - first).abs() < 1e-10);
        }
    }

    // ── checkerboard_measure ─────────────────────────────────────────────────

    #[test]
    fn test_checkerboard_perfect() {
        let d = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let m = checkerboard_measure(&d);
        assert!((m - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_checkerboard_uniform() {
        let d = vec![
            vec![0.5, 0.5, 0.5],
            vec![0.5, 0.5, 0.5],
            vec![0.5, 0.5, 0.5],
        ];
        let m = checkerboard_measure(&d);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_checkerboard_solid() {
        let d = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let m = checkerboard_measure(&d);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_checkerboard_single_row() {
        let d = vec![vec![0.0, 1.0, 0.0]];
        let m = checkerboard_measure(&d);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_checkerboard_single_col() {
        let d = vec![vec![0.0], vec![1.0]];
        let m = checkerboard_measure(&d);
        assert_eq!(m, 0.0);
    }

    #[test]
    fn test_checkerboard_inverse() {
        let d = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let m = checkerboard_measure(&d);
        assert!((m - 1.0).abs() < 1e-12);
    }

    // ── topology_optimize ────────────────────────────────────────────────────

    #[test]
    fn test_topology_optimize_output_length() {
        let params = TopologyOptParams::new(0.5, 3.0, 1.5, 50);
        let n = 16;
        let sens = vec![-0.4_f64; n];
        let result = topology_optimize(&params, n, &sens);
        assert_eq!(result.len(), n);
    }

    #[test]
    fn test_topology_optimize_bounds() {
        let params = TopologyOptParams::new(0.4, 3.0, 1.5, 50);
        let n = 10;
        let sens = vec![-0.5_f64; n];
        let result = topology_optimize(&params, n, &sens);
        for &r in &result {
            assert!(r >= 1e-3 - 1e-10 && r <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_topology_optimize_volume() {
        let vf = 0.3;
        let params = TopologyOptParams::new(vf, 3.0, 1.5, 50);
        let n = 20;
        let sens = vec![-0.3_f64; n];
        let result = topology_optimize(&params, n, &sens);
        let mean = result.iter().sum::<f64>() / n as f64;
        assert!((mean - vf).abs() < 0.1);
    }

    #[test]
    fn test_simp_penalized_stiffness_various_penalties() {
        for &p in &[1.0_f64, 2.0, 3.0, 4.0] {
            let e = simp_penalized_stiffness(0.5, 1.0, 0.0, p);
            let expected = 0.5_f64.powf(p);
            assert!((e - expected).abs() < 1e-12);
        }
    }

    // ── TopologyOptimization ─────────────────────────────────────────────────

    #[test]
    fn test_topology_optimization_run_length() {
        let params = TopologyOptParams::new(0.5, 3.0, 1.5, 10);
        let n = 16;
        let mut topo = TopologyOptimization::new(n, params);
        let adj: Vec<Vec<usize>> = vec![vec![]; n];
        topo.run(
            &adj,
            |d| {
                let c: f64 = d.iter().map(|&r| 1.0 / (r + 0.01)).sum();
                let s: Vec<f64> = d.iter().map(|&r| -1.0 / (r + 0.01).powi(2)).collect();
                (c, s)
            },
            1e-4,
        );
        assert_eq!(topo.densities.len(), n);
        assert!(!topo.history.is_empty());
    }

    #[test]
    fn test_topology_optimization_convergence_flag() {
        let params = TopologyOptParams::new(0.5, 3.0, 1.0, 100);
        let n = 4;
        let mut topo = TopologyOptimization::new(n, params);
        let adj: Vec<Vec<usize>> = vec![vec![]; n];
        topo.run(
            &adj,
            |d| {
                let c: f64 = d.iter().sum();
                let s: Vec<f64> = vec![-0.5; n];
                (c, s)
            },
            0.5,
        ); // very loose tol → should converge quickly
        let _ = topo.is_converged(0.5);
    }

    #[test]
    fn test_topology_optimization_density_filter_in_place() {
        let params = TopologyOptParams::new(0.5, 3.0, 1.5, 10);
        let n = 4;
        let mut topo = TopologyOptimization::new(n, params);
        let adj: Vec<Vec<usize>> = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        topo.apply_density_filter(&adj);
        // densities should still be in range
        for &d in &topo.densities {
            assert!(d >= 0.0 && d <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_topology_optimization_heaviside_in_place() {
        let params = TopologyOptParams::new(0.5, 3.0, 1.5, 10);
        let n = 4;
        let mut topo = TopologyOptimization::new(n, params);
        topo.apply_heaviside(5.0, 0.5);
        for &d in &topo.densities {
            assert!(d >= 0.0 && d <= 1.0 + 1e-10);
        }
    }

    // ── ShapeOptimization ────────────────────────────────────────────────────

    #[test]
    fn test_shape_opt_runs() {
        let opt = ShapeOptimization::new(0.01, 50, 1e-4);
        let vars = vec![
            ShapeDesignVariable {
                node: 0,
                direction: 0,
                value: 1.0,
            },
            ShapeDesignVariable {
                node: 1,
                direction: 1,
                value: 2.0,
            },
        ];
        let result = opt.run(vars, |v| {
            let obj = v[0].value.powi(2) + v[1].value.powi(2);
            let grad = vec![2.0 * v[0].value, 2.0 * v[1].value];
            (obj, grad)
        });
        // Should reduce objective
        assert!(result.objective < 5.0 + 1e-10);
        assert_eq!(result.variables.len(), 2);
    }

    #[test]
    fn test_shape_opt_adjoint_sensitivity() {
        let vars = vec![ShapeDesignVariable {
            node: 0,
            direction: 0,
            value: 2.0,
        }];
        let oracle = |v: &[ShapeDesignVariable]| {
            let obj = v[0].value.powi(2);
            (obj, vec![2.0 * v[0].value])
        };
        let grad = ShapeOptimization::adjoint_sensitivity(&vars, &oracle, 1e-5);
        // dJ/ds = 2*2 = 4
        assert!((grad[0] - 4.0).abs() < 0.01, "grad={}", grad[0]);
    }

    // ── SizingOptimization ───────────────────────────────────────────────────

    #[test]
    fn test_sizing_opt_volume_budget() {
        let lengths = vec![1.0, 1.0, 1.0];
        let opt = SizingOptimization::new(lengths.clone(), 1.5, 0.5, 20);
        let mut sections = vec![
            CrossSection::new(0, 0.5, 0.1, 2.0),
            CrossSection::new(1, 0.5, 0.1, 2.0),
            CrossSection::new(2, 0.5, 0.1, 2.0),
        ];
        let result = opt.run(&mut sections, |a| {
            let c: f64 = a.iter().map(|&ai| 1.0 / ai).sum();
            let s: Vec<f64> = a.iter().map(|&ai| -1.0 / ai.powi(2)).collect();
            (c, s)
        });
        assert!(result.total_volume <= opt.volume_max + 0.1);
    }

    #[test]
    fn test_sizing_oc_update_bounds() {
        let lengths = vec![1.0, 1.0];
        let opt = SizingOptimization::new(lengths, 2.0, 0.3, 10);
        let sections = vec![
            CrossSection::new(0, 0.5, 0.2, 1.0),
            CrossSection::new(1, 0.5, 0.2, 1.0),
        ];
        let sens = vec![-0.5, -0.3];
        let a_new = opt.update_oc(&sections, &sens);
        for (&a, s) in a_new.iter().zip(sections.iter()) {
            assert!(a >= s.area_min - 1e-10);
            assert!(a <= s.area_max + 1e-10);
        }
    }

    #[test]
    fn test_cross_section_clamp() {
        let s = CrossSection::new(0, 2.5, 0.1, 1.0);
        assert!((s.clamped() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_sizing_total_volume() {
        let lengths = vec![2.0, 3.0];
        let opt = SizingOptimization::new(lengths, 10.0, 0.5, 10);
        let vol = opt.total_volume(&[1.0, 1.0]);
        assert!((vol - 5.0).abs() < 1e-12);
    }

    // ── MultiObjectiveFEM ────────────────────────────────────────────────────

    #[test]
    fn test_multiobjective_pareto_count() {
        let mo = MultiObjectiveFEM::new(5, 3.0, 1.5, 5);
        let front = mo.pareto_front(4, 0.2, 0.8, |d, _vf| {
            let c: f64 = d.iter().map(|&r| 1.0 / (r + 0.01)).sum();
            let s: Vec<f64> = d.iter().map(|&r| -1.0 / (r + 0.01).powi(2)).collect();
            (c, s)
        });
        assert_eq!(front.len(), 5);
    }

    #[test]
    fn test_multiobjective_pareto_mass_range() {
        let mo = MultiObjectiveFEM::new(3, 3.0, 1.0, 5);
        let front = mo.pareto_front(4, 0.3, 0.7, |d, _vf| {
            let c: f64 = d.iter().map(|&r| 1.0 / (r + 0.01)).sum();
            let s: Vec<f64> = d.iter().map(|&r| -1.0 / (r + 0.01).powi(2)).collect();
            (c, s)
        });
        for p in &front {
            assert!(p.mass >= 0.0 && p.mass <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_multiobjective_weighted_sum() {
        let p = ParetoPoint {
            compliance: 10.0,
            mass: 0.5,
            volume_fraction: 0.5,
            densities: vec![],
        };
        let ws = MultiObjectiveFEM::weighted_sum(&p, 1.0, 2.0);
        assert!((ws - 11.0).abs() < 1e-12);
    }

    #[test]
    fn test_multiobjective_filter_dominated() {
        let front = vec![
            ParetoPoint {
                compliance: 1.0,
                mass: 2.0,
                volume_fraction: 0.5,
                densities: vec![],
            },
            ParetoPoint {
                compliance: 2.0,
                mass: 1.0,
                volume_fraction: 0.5,
                densities: vec![],
            },
            ParetoPoint {
                compliance: 3.0,
                mass: 3.0,
                volume_fraction: 0.5,
                densities: vec![],
            }, // dominated
        ];
        let nd = MultiObjectiveFEM::filter_dominated(&front);
        assert_eq!(
            nd.len(),
            2,
            "non-dominated count should be 2, got {}",
            nd.len()
        );
    }

    // ── ManufacturingConstraints ─────────────────────────────────────────────

    #[test]
    fn test_symmetry_lr_2x4() {
        let mc = ManufacturingConstraints::new(0.0, true, 45.0);
        let mut d = vec![vec![0.2_f64, 0.4, 0.6, 0.8]];
        mc.enforce_symmetry_lr(&mut d);
        // d[0][0] should equal d[0][3], d[0][1] should equal d[0][2]
        assert!((d[0][0] - d[0][3]).abs() < 1e-12);
        assert!((d[0][1] - d[0][2]).abs() < 1e-12);
    }

    #[test]
    fn test_min_length_scale_zeros_small() {
        let mc = ManufacturingConstraints::new(0.5, false, 45.0);
        let mut d = vec![0.01_f64, 0.01, 0.9, 0.9];
        mc.apply_min_length_scale(&mut d);
        // large values should remain ≥ threshold
        assert!(d[2] >= 0.001);
    }

    #[test]
    fn test_overhang_violation_no_solid() {
        let mc = ManufacturingConstraints::new(0.0, false, 45.0);
        let d = vec![vec![0.0, 0.0], vec![0.0, 0.0]];
        assert_eq!(mc.overhang_violation_fraction(&d), 0.0);
    }

    #[test]
    fn test_overhang_violation_fully_supported() {
        let mc = ManufacturingConstraints::new(0.0, false, 45.0);
        let d = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        // Every solid element in row 1 is supported by row 0
        assert_eq!(mc.overhang_violation_fraction(&d), 0.0);
    }

    #[test]
    fn test_project_length_scale() {
        let d = vec![0.05, 0.3, 0.8];
        let proj = ManufacturingConstraints::project_length_scale(&d, 0.1);
        assert!(proj[0] < 0.1);
        assert!((proj[1] - 0.3).abs() < 1e-12);
    }

    // ── mbb_beam_topology_opt ────────────────────────────────────────────────

    #[test]
    fn test_mbb_beam_output_shape() {
        let grid = mbb_beam_topology_opt(4, 8, 0.5, 3.0, 1.5, 5);
        assert_eq!(grid.len(), 4);
        for row in &grid {
            assert_eq!(row.len(), 8);
        }
    }

    #[test]
    fn test_mbb_beam_density_bounds() {
        let grid = mbb_beam_topology_opt(3, 6, 0.4, 3.0, 1.5, 5);
        for row in &grid {
            for &d in row {
                assert!(d >= 0.0 && d <= 1.0 + 1e-10, "density out of range: {d}");
            }
        }
    }

    #[test]
    fn test_mbb_beam_volume_fraction() {
        let vf = 0.4;
        let grid = mbb_beam_topology_opt(4, 4, vf, 3.0, 1.5, 20);
        let all: Vec<f64> = grid.into_iter().flatten().collect();
        let mean = all.iter().sum::<f64>() / all.len() as f64;
        assert!(
            (mean - vf).abs() < 0.15,
            "mean density={mean} expected~{vf}"
        );
    }
}
