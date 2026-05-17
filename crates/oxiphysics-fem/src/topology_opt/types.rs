//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    augmented_lagrangian_volume_penalty, heaviside_derivative, heaviside_projection,
};

/// Configuration for a topology optimization run.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TopOptConfig {
    /// Target volume fraction.
    pub vf_target: f64,
    /// SIMP penalty exponent.
    pub penalty: f64,
    /// Filter radius.
    pub r_min: f64,
    /// Maximum number of iterations.
    pub n_iter: usize,
    /// Convergence tolerance on density change.
    pub tol: f64,
}
/// SIMP interpolation model: E(ρ) = E_min + ρ^p * (E_0 - E_min).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SimpModel {
    /// Penalization exponent.
    pub penalty: f64,
    /// Void Young's modulus (minimum).
    pub e_min: f64,
    /// Solid Young's modulus.
    pub e_0: f64,
}
#[allow(dead_code)]
impl SimpModel {
    /// Effective Young's modulus at density `rho`.
    pub fn young_modulus(&self, rho: f64) -> f64 {
        self.e_min + rho.powf(self.penalty) * (self.e_0 - self.e_min)
    }
    /// Sensitivity of compliance w.r.t. density: dC/dρ = -p * ρ^(p-1) * (E_0 - E_min) * u_e.
    ///
    /// `strain_energy` is the element strain energy u_e^T K_e u_e.
    pub fn sensitivity(&self, rho: f64, strain_energy: f64) -> f64 {
        -self.penalty * rho.powf(self.penalty - 1.0) * (self.e_0 - self.e_min) * strain_energy
    }
}
/// Projection filter combining density filtering with Heaviside thresholding.
///
/// Applies the chain: ρ_tilde → (density filter) → ρ_bar → (Heaviside) → ρ_phys.
///
/// The chain rule for the sensitivity is:
/// `∂c/∂ρ = ∂c/∂ρ_phys * ∂H/∂ρ_bar * ∂ρ_bar/∂ρ`
pub struct ProjectionFilter {
    /// Filter radius.
    pub radius: f64,
    /// Heaviside projection sharpness.
    pub beta: f64,
    /// Heaviside threshold.
    pub eta: f64,
}
impl ProjectionFilter {
    /// Create a new projection filter.
    pub fn new(radius: f64, beta: f64, eta: f64) -> Self {
        Self { radius, beta, eta }
    }
    /// Apply the full projection filter chain to the grid.
    ///
    /// Returns `(rho_bar, rho_phys)` where:
    /// - `rho_bar`  is the density-filtered field
    /// - `rho_phys` is the Heaviside-projected physical density
    pub fn apply(&self, grid: &TopologyGrid) -> (Vec<f64>, Vec<f64>) {
        let rho_bar = grid.filter_densities(self.radius);
        let rho_phys: Vec<f64> = rho_bar
            .iter()
            .map(|&r| heaviside_projection(r, self.beta, self.eta))
            .collect();
        (rho_bar, rho_phys)
    }
    /// Compute the sensitivity chain rule for the projection filter.
    ///
    /// `dc_drho[i] = dc_drho_phys[i] * dH_drho_bar[i] * w_sum[i]`
    ///
    /// where `w_sum[i] = sum_j w_ij / (sum_j w_ij * rho_j)` is the density
    /// filter weight sum.
    #[allow(clippy::too_many_arguments)]
    pub fn chain_sensitivity(
        &self,
        grid: &TopologyGrid,
        dc_drho_phys: &[f64],
        rho_bar: &[f64],
    ) -> Vec<f64> {
        let n = grid.elements.len();
        let dh: Vec<f64> = rho_bar
            .iter()
            .map(|&r| heaviside_derivative(r, self.beta, self.eta))
            .collect();
        let dc_times_dh: Vec<f64> = (0..n).map(|i| dc_drho_phys[i] * dh[i]).collect();
        let mut filtered = vec![0.0f64; n];
        for i in 0..n {
            let ei = &grid.elements[i];
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            for j in 0..n {
                let ej = &grid.elements[j];
                let dist = ((ei.x0 - ej.x0).powi(2) + (ei.y0 - ej.y0).powi(2)).sqrt();
                if dist < self.radius {
                    let w = self.radius - dist;
                    num += w * dc_times_dh[j];
                    den += w;
                }
            }
            filtered[i] = if den > 1e-30 {
                num / den
            } else {
                dc_times_dh[i]
            };
        }
        filtered
    }
}
/// Stores the last `n_history` sensitivity vectors and provides their average.
///
/// Averaging sensitivities over multiple iterations reduces numerical
/// oscillations and improves convergence of the OC method.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SensitivityHistory {
    /// Number of elements.
    pub n_elements: usize,
    /// Maximum number of stored iterations.
    pub n_history: usize,
    /// Circular buffer of sensitivity vectors.
    pub(super) buffer: Vec<Vec<f64>>,
    /// Write cursor into the circular buffer.
    pub(super) cursor: usize,
    /// Number of entries actually filled (≤ n_history).
    pub(super) filled: usize,
}
impl SensitivityHistory {
    /// Create a new sensitivity history buffer.
    pub fn new(n_elements: usize, n_history: usize) -> Self {
        let n_hist = n_history.max(1);
        Self {
            n_elements,
            n_history: n_hist,
            buffer: vec![vec![0.0; n_elements]; n_hist],
            cursor: 0,
            filled: 0,
        }
    }
    /// Insert the sensitivities for the current iteration.
    pub fn add_iteration(&mut self, sensitivities: &[f64]) {
        assert_eq!(sensitivities.len(), self.n_elements);
        self.buffer[self.cursor].copy_from_slice(sensitivities);
        self.cursor = (self.cursor + 1) % self.n_history;
        if self.filled < self.n_history {
            self.filled += 1;
        }
    }
    /// Compute element-wise average over stored iterations.
    pub fn average_sensitivity(&self) -> Vec<f64> {
        if self.filled == 0 {
            return vec![0.0; self.n_elements];
        }
        let mut avg = vec![0.0; self.n_elements];
        for iter_buf in self.buffer.iter().take(self.filled) {
            for (a, &s) in avg.iter_mut().zip(iter_buf.iter()) {
                *a += s;
            }
        }
        let inv = 1.0 / self.filled as f64;
        avg.iter_mut().for_each(|v| *v *= inv);
        avg
    }
    /// Clear all stored history.
    pub fn reset(&mut self) {
        self.cursor = 0;
        self.filled = 0;
        for buf in self.buffer.iter_mut() {
            buf.iter_mut().for_each(|v| *v = 0.0);
        }
    }
}
/// Parameters for the BESO method.
#[derive(Debug, Clone)]
pub struct BesoParams {
    /// Target volume fraction.
    pub target_vf: f64,
    /// Evolutionary rate (fraction of elements switched per iteration).
    pub er: f64,
    /// Penalization exponent (same meaning as SIMP penalty).
    pub penalty: f64,
    /// Minimum density for "solid" elements.
    pub x_max: f64,
    /// Density for "void" elements.
    pub x_min: f64,
}
impl BesoParams {
    /// Default BESO parameters.
    pub fn default_structural() -> Self {
        Self {
            target_vf: 0.5,
            er: 0.02,
            penalty: 3.0,
            x_max: 1.0,
            x_min: 1e-3,
        }
    }
}
/// Parameters for the SIMP topology-optimization method.
#[derive(Debug, Clone)]
pub struct SimpParams {
    /// Full (solid) Young's modulus.
    pub e0: f64,
    /// Void modulus (numerically small to avoid singularity).
    pub emin: f64,
    /// Penalization exponent (typically 3).
    pub penalty: f64,
    /// Target volume fraction (0 < vf ≤ 1).
    pub volume_fraction: f64,
    /// Filter radius for sensitivity smoothing.
    pub filter_radius: f64,
    /// Maximum allowable density change per OC step.
    pub move_limit: f64,
}
impl SimpParams {
    /// Default parameters representative of structural steel.
    pub fn default_steel() -> Self {
        let e0 = 210e9_f64;
        Self {
            e0,
            emin: 1e-9 * e0,
            penalty: 3.0,
            volume_fraction: 0.5,
            filter_radius: 1.5,
            move_limit: 0.2,
        }
    }
    /// Parameters for aluminum.
    pub fn default_aluminum() -> Self {
        let e0 = 70e9_f64;
        Self {
            e0,
            emin: 1e-9 * e0,
            penalty: 3.0,
            volume_fraction: 0.4,
            filter_radius: 1.5,
            move_limit: 0.2,
        }
    }
}
/// A structured 2-D grid of [`SimpElement`]s.
pub struct TopologyGrid {
    /// Row-major flat array of elements (row = y, col = x).
    pub elements: Vec<SimpElement>,
    /// Number of columns.
    pub nx: usize,
    /// Number of rows.
    pub ny: usize,
}
impl TopologyGrid {
    /// Create a uniform grid with all densities set to `rho0`.
    pub fn new_uniform(nx: usize, ny: usize, rho0: f64, e0: f64, dx: f64, dy: f64) -> Self {
        let emin = 1e-9 * e0;
        let penalty = 3.0;
        let mut elements = Vec::with_capacity(nx * ny);
        for row in 0..ny {
            for col in 0..nx {
                let x0 = (col as f64 + 0.5) * dx;
                let y0 = (row as f64 + 0.5) * dy;
                let stiffness = emin + rho0.powf(penalty) * (e0 - emin);
                elements.push(SimpElement {
                    rho: rho0,
                    stiffness,
                    sensitivity: 0.0,
                    x0,
                    y0,
                    dx,
                    dy,
                });
            }
        }
        Self { elements, nx, ny }
    }
    /// Average density across all elements.
    pub fn volume_fraction(&self) -> f64 {
        let sum: f64 = self.elements.iter().map(|e| e.rho).sum();
        sum / self.elements.len() as f64
    }
    /// Volume-weighted average density (accounts for different element sizes).
    pub fn volume_fraction_weighted(&self) -> f64 {
        let total_area: f64 = self.elements.iter().map(|e| e.area()).sum();
        if total_area < 1e-30 {
            return 0.0;
        }
        let weighted: f64 = self.elements.iter().map(|e| e.rho * e.area()).sum();
        weighted / total_area
    }
    /// Proxy compliance: Σ sensitivity_i · ρ_i (used when full FEA is skipped).
    pub fn total_compliance(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| e.sensitivity * e.rho)
            .sum::<f64>()
            .abs()
    }
    /// Density-weighted sensitivity filter of radius `radius` (hat function).
    pub fn filter_sensitivities(&self, radius: f64) -> Vec<f64> {
        let n = self.elements.len();
        let mut filtered = vec![0.0_f64; n];
        for i in 0..n {
            let ei = &self.elements[i];
            let mut num = 0.0_f64;
            let mut den = 0.0_f64;
            for j in 0..n {
                let ej = &self.elements[j];
                let dist = ((ei.x0 - ej.x0).powi(2) + (ei.y0 - ej.y0).powi(2)).sqrt();
                if dist < radius {
                    let w = radius - dist;
                    num += w * ej.rho * ej.sensitivity;
                    den += w * ej.rho;
                }
            }
            filtered[i] = if den.abs() > 1e-30 {
                num / den
            } else {
                ei.sensitivity
            };
        }
        filtered
    }
    /// Density filter (Bruns & Tortorelli): ρ_filtered = Σ w_j ρ_j / Σ w_j.
    ///
    /// Unlike the sensitivity filter, this directly filters the density field.
    pub fn filter_densities(&self, radius: f64) -> Vec<f64> {
        let n = self.elements.len();
        let mut filtered = vec![0.0_f64; n];
        for i in 0..n {
            let ei = &self.elements[i];
            let mut num = 0.0_f64;
            let mut den = 0.0_f64;
            for j in 0..n {
                let ej = &self.elements[j];
                let dist = ((ei.x0 - ej.x0).powi(2) + (ei.y0 - ej.y0).powi(2)).sqrt();
                if dist < radius {
                    let w = radius - dist;
                    num += w * ej.rho;
                    den += w;
                }
            }
            filtered[i] = if den.abs() > 1e-30 { num / den } else { ei.rho };
        }
        filtered
    }
    /// OC density update with bisection to enforce the volume constraint.
    pub fn update_densities_oc(&mut self, sensitivities: &[f64], params: &SimpParams) {
        let mut l1 = 0.0_f64;
        let mut l2 = 1e9_f64;
        let target = params.volume_fraction * self.elements.len() as f64;
        for _ in 0..50 {
            let lmid = 0.5 * (l1 + l2);
            let vol: f64 = self
                .elements
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    let be = (-sensitivities[i] / lmid).max(0.0).sqrt();
                    (e.rho * be)
                        .max(e.rho - params.move_limit)
                        .min(e.rho + params.move_limit)
                        .clamp(1e-3, 1.0)
                })
                .sum();
            if vol > target {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
        }
        let lmid = 0.5 * (l1 + l2);
        for (i, e) in self.elements.iter_mut().enumerate() {
            let be = (-sensitivities[i] / lmid).max(0.0).sqrt();
            e.rho = (e.rho * be)
                .max(e.rho - params.move_limit)
                .min(e.rho + params.move_limit)
                .clamp(1e-3, 1.0);
            e.stiffness = params.emin + e.rho.powf(params.penalty) * (params.e0 - params.emin);
        }
    }
    /// Method of Moving Asymptotes (MMA) simplified density update.
    ///
    /// Uses a conservative convex separable approximation.
    pub fn update_densities_mma(
        &mut self,
        sensitivities: &[f64],
        params: &SimpParams,
        lower_asymptotes: &[f64],
        upper_asymptotes: &[f64],
    ) {
        let n = self.elements.len();
        let target_vol = params.volume_fraction * n as f64;
        let mut l1 = 0.0_f64;
        let mut l2 = 1e9_f64;
        for _ in 0..50 {
            let lmid = 0.5 * (l1 + l2);
            let vol: f64 = (0..n)
                .map(|i| {
                    let e = &self.elements[i];
                    let lo = lower_asymptotes.get(i).copied().unwrap_or(0.0);
                    let hi = upper_asymptotes.get(i).copied().unwrap_or(1.0);
                    let be = (-sensitivities[i] / lmid).max(0.0).sqrt();
                    (e.rho * be)
                        .max(e.rho - params.move_limit)
                        .min(e.rho + params.move_limit)
                        .clamp(lo.max(1e-3), hi.min(1.0))
                })
                .sum();
            if vol > target_vol {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
        }
        let lmid = 0.5 * (l1 + l2);
        for i in 0..n {
            let lo = lower_asymptotes.get(i).copied().unwrap_or(0.0);
            let hi = upper_asymptotes.get(i).copied().unwrap_or(1.0);
            let e = &mut self.elements[i];
            let be = (-sensitivities[i] / lmid).max(0.0).sqrt();
            e.rho = (e.rho * be)
                .max(e.rho - params.move_limit)
                .min(e.rho + params.move_limit)
                .clamp(lo.max(1e-3), hi.min(1.0));
            e.stiffness = params.emin + e.rho.powf(params.penalty) * (params.e0 - params.emin);
        }
    }
    /// Get element index from row, col.
    pub fn element_index(&self, row: usize, col: usize) -> usize {
        row * self.nx + col
    }
    /// Get the (row, col) from a flat index.
    pub fn row_col(&self, index: usize) -> (usize, usize) {
        (index / self.nx, index % self.nx)
    }
    /// Compute gradient magnitude of the density field (central differences).
    pub fn density_gradient_magnitude(&self) -> Vec<f64> {
        let n = self.elements.len();
        let mut grad = vec![0.0f64; n];
        for row in 0..self.ny {
            for col in 0..self.nx {
                let idx = self.element_index(row, col);
                let rho = self.elements[idx].rho;
                let dx_grad = if col > 0 && col + 1 < self.nx {
                    let left = self.elements[self.element_index(row, col - 1)].rho;
                    let right = self.elements[self.element_index(row, col + 1)].rho;
                    (right - left) / (2.0 * self.elements[idx].dx)
                } else if col + 1 < self.nx {
                    let right = self.elements[self.element_index(row, col + 1)].rho;
                    (right - rho) / self.elements[idx].dx
                } else if col > 0 {
                    let left = self.elements[self.element_index(row, col - 1)].rho;
                    (rho - left) / self.elements[idx].dx
                } else {
                    0.0
                };
                let dy_grad = if row > 0 && row + 1 < self.ny {
                    let below = self.elements[self.element_index(row - 1, col)].rho;
                    let above = self.elements[self.element_index(row + 1, col)].rho;
                    (above - below) / (2.0 * self.elements[idx].dy)
                } else if row + 1 < self.ny {
                    let above = self.elements[self.element_index(row + 1, col)].rho;
                    (above - rho) / self.elements[idx].dy
                } else if row > 0 {
                    let below = self.elements[self.element_index(row - 1, col)].rho;
                    (rho - below) / self.elements[idx].dy
                } else {
                    0.0
                };
                grad[idx] = (dx_grad * dx_grad + dy_grad * dy_grad).sqrt();
            }
        }
        grad
    }
}
/// Manages the augmented Lagrangian multiplier for volume constraints.
///
/// Dual update: λ ← λ + ρ_aug * (V(ρ) - V_target)
/// Penalty update: ρ_aug ← min(γ * ρ_aug, ρ_max)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AugmentedLagrangianVolumeConstraint {
    /// Target volume fraction.
    pub target_vf: f64,
    /// Augmented Lagrangian penalty parameter.
    pub rho_aug: f64,
    /// Lagrange multiplier estimate.
    pub lambda: f64,
    /// Maximum allowed penalty parameter.
    pub rho_max: f64,
    /// Penalty growth factor.
    pub gamma: f64,
}
impl AugmentedLagrangianVolumeConstraint {
    /// Create a new augmented Lagrangian constraint handler.
    pub fn new(target_vf: f64, rho_aug: f64) -> Self {
        Self {
            target_vf,
            rho_aug,
            lambda: 0.0,
            rho_max: 1e8,
            gamma: 2.0,
        }
    }
    /// Update the Lagrange multiplier after computing the current volume fraction.
    pub fn update_multiplier(&mut self, current_vf: f64) {
        let g = current_vf - self.target_vf;
        self.lambda += self.rho_aug * g;
    }
    /// Grow the penalty parameter (called when constraint violation does not improve).
    pub fn grow_penalty(&mut self) {
        self.rho_aug = (self.rho_aug * self.gamma).min(self.rho_max);
    }
    /// Compute effective gradient contribution for this constraint.
    ///
    /// Returns the multiplier-weighted gradient entry for element `e`.
    pub fn gradient_contribution(
        &self,
        element_volume: f64,
        total_volume: f64,
        current_vf: f64,
    ) -> f64 {
        let g = current_vf - self.target_vf;
        let dg_drho_e = element_volume / total_volume;
        (self.lambda + self.rho_aug * g) * dg_drho_e
    }
    /// Compute the augmented Lagrangian penalty value.
    pub fn penalty_value(&self, current_vf: f64) -> f64 {
        augmented_lagrangian_volume_penalty(current_vf, self.target_vf, self.lambda, self.rho_aug)
    }
}
/// Level-set function field for topology optimization.
///
/// The material domain is defined by `phi > 0` (solid) and `phi < 0` (void).
/// The interface is the zero-level set `phi = 0`.
#[derive(Debug, Clone)]
pub struct LevelSetField {
    /// Flat (row-major) level-set function values.
    pub phi: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Cell size.
    pub dx: f64,
    /// Cell size.
    pub dy: f64,
}
impl LevelSetField {
    /// Create a level-set field initialized to a signed-distance function
    /// from the boundary of the domain (`phi = vf - 0.5` uniform).
    pub fn new_uniform(nx: usize, ny: usize, dx: f64, dy: f64, vf: f64) -> Self {
        let phi_init = vf - 0.5;
        Self {
            phi: vec![phi_init; nx * ny],
            nx,
            ny,
            dx,
            dy,
        }
    }
    /// Create a level-set field from a pre-computed signed-distance function.
    pub fn from_sdf(nx: usize, ny: usize, dx: f64, dy: f64, phi: Vec<f64>) -> Self {
        assert_eq!(phi.len(), nx * ny);
        Self {
            phi,
            nx,
            ny,
            dx,
            dy,
        }
    }
    /// Compute the approximate volume fraction (fraction of cells with phi ≥ 0).
    pub fn volume_fraction(&self) -> f64 {
        let solid = self.phi.iter().filter(|&&p| p >= 0.0).count();
        solid as f64 / self.phi.len() as f64
    }
    /// Evaluate the Heaviside function H(phi) for material interpolation.
    ///
    /// Uses a smooth regularized Heaviside:
    /// `H(phi) = 0.5 * (1 + phi/eps + sin(pi*phi/eps)/pi)` for |phi| <= eps
    pub fn heaviside_smooth(phi: f64, eps: f64) -> f64 {
        if phi < -eps {
            0.0
        } else if phi > eps {
            1.0
        } else {
            0.5 * (1.0
                + phi / eps
                + (std::f64::consts::PI * phi / eps).sin() / std::f64::consts::PI)
        }
    }
    /// Evaluate the smooth Dirac delta δ(phi).
    ///
    /// `delta(phi) = 0.5/eps * (1 + cos(pi*phi/eps))` for |phi| <= eps
    pub fn dirac_smooth(phi: f64, eps: f64) -> f64 {
        if phi.abs() > eps {
            0.0
        } else {
            0.5 / eps * (1.0 + (std::f64::consts::PI * phi / eps).cos())
        }
    }
    /// Compute the element-wise material indicator from the level-set field.
    ///
    /// Returns H(phi_i) for each element `i`.
    pub fn material_indicator(&self, eps: f64) -> Vec<f64> {
        self.phi
            .iter()
            .map(|&p| Self::heaviside_smooth(p, eps))
            .collect()
    }
    /// Advance the level-set function by one pseudo-time step using
    /// a simple upwind scheme for the Hamilton-Jacobi equation:
    ///
    /// `∂φ/∂t + V_n |∇φ| = 0`
    ///
    /// where `V_n` is the normal velocity (proportional to the sensitivity).
    pub fn advect(&mut self, velocity: &[f64], dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let dy = self.dy;
        let mut phi_new = self.phi.clone();
        for row in 0..ny {
            for col in 0..nx {
                let idx = row * nx + col;
                let vn = velocity[idx];
                let dphidx = if vn >= 0.0 {
                    if col > 0 {
                        (self.phi[idx] - self.phi[row * nx + col - 1]) / dx
                    } else {
                        0.0
                    }
                } else if col + 1 < nx {
                    (self.phi[row * nx + col + 1] - self.phi[idx]) / dx
                } else {
                    0.0
                };
                let dphidy = if vn >= 0.0 {
                    if row > 0 {
                        (self.phi[idx] - self.phi[(row - 1) * nx + col]) / dy
                    } else {
                        0.0
                    }
                } else if row + 1 < ny {
                    (self.phi[(row + 1) * nx + col] - self.phi[idx]) / dy
                } else {
                    0.0
                };
                let grad_phi = (dphidx * dphidx + dphidy * dphidy).sqrt();
                phi_new[idx] = self.phi[idx] - dt * vn * grad_phi;
            }
        }
        self.phi = phi_new;
    }
    /// Re-initialize the level-set field to a signed-distance function using
    /// a simple fast sweeping method (one pass, approximate).
    pub fn reinitialize(&mut self, max_iter: usize) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        for _ in 0..max_iter {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = j * nx + i;
                    let phi_ij = self.phi[idx];
                    let mut candidates = vec![phi_ij];
                    if i > 0 {
                        candidates.push(self.phi[j * nx + i - 1] + dx * phi_ij.signum());
                    }
                    if j > 0 {
                        candidates.push(self.phi[(j - 1) * nx + i] + dx * phi_ij.signum());
                    }
                    if phi_ij >= 0.0 {
                        self.phi[idx] = candidates.iter().cloned().fold(f64::INFINITY, f64::min);
                    } else {
                        self.phi[idx] =
                            candidates.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    }
                }
            }
            for j in (0..ny).rev() {
                for i in (0..nx).rev() {
                    let idx = j * nx + i;
                    let phi_ij = self.phi[idx];
                    let mut candidates = vec![phi_ij];
                    if i + 1 < nx {
                        candidates.push(self.phi[j * nx + i + 1] + dx * phi_ij.signum());
                    }
                    if j + 1 < ny {
                        candidates.push(self.phi[(j + 1) * nx + i] + dx * phi_ij.signum());
                    }
                    if phi_ij >= 0.0 {
                        self.phi[idx] = candidates.iter().cloned().fold(f64::INFINITY, f64::min);
                    } else {
                        self.phi[idx] =
                            candidates.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    }
                }
            }
        }
    }
    /// Compute the level-set sensitivity velocity from element sensitivities.
    ///
    /// For compliance minimization: `V_n = -dC/dφ`.
    /// Uses the smooth Dirac delta to localize sensitivity to the interface.
    pub fn sensitivity_velocity(&self, sensitivities: &[f64], eps: f64) -> Vec<f64> {
        let n = self.phi.len();
        let mut vel = vec![0.0f64; n];
        for i in 0..n {
            let delta = Self::dirac_smooth(self.phi[i], eps);
            vel[i] = -sensitivities[i] * delta;
        }
        vel
    }
}
/// RAMP interpolation: E(ρ) = E_min + ρ / (1 + q*(1-ρ)) * (E_0 - E_min).
///
/// The RAMP scheme avoids the vanishing derivative at ρ=0 that SIMP exhibits,
/// making it better suited for problems where intermediate densities occur.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RampModel {
    /// Penalization parameter q (typically 3–8).
    pub q: f64,
    /// Void Young's modulus (numerical floor).
    pub e_min: f64,
    /// Solid Young's modulus.
    pub e_0: f64,
}
impl RampModel {
    /// Construct a new RAMP model.
    pub fn new(q: f64, e_min: f64, e_0: f64) -> Self {
        assert!(q > 0.0, "q must be positive");
        assert!(e_min > 0.0, "E_min must be positive");
        assert!(e_0 >= e_min, "E_0 must be >= E_min");
        Self { q, e_min, e_0 }
    }
    /// Interpolated Young's modulus at density `rho ∈ [0, 1]`.
    pub fn young_modulus(&self, rho: f64) -> f64 {
        let rho = rho.clamp(0.0, 1.0);
        self.e_min + rho / (1.0 + self.q * (1.0 - rho)) * (self.e_0 - self.e_min)
    }
    /// Derivative ∂E/∂ρ at density `rho`.
    pub fn derivative(&self, rho: f64) -> f64 {
        let rho = rho.clamp(0.0, 1.0);
        let denom = 1.0 + self.q * (1.0 - rho);
        (1.0 + self.q) / (denom * denom) * (self.e_0 - self.e_min)
    }
    /// Sensitivity of compliance w.r.t. density for a single element.
    ///
    /// dC/dρ = -dE/dρ * u^T * k_0 * u  (element compliance)
    pub fn compliance_sensitivity(&self, rho: f64, element_compliance: f64) -> f64 {
        let denom = 1.0 + self.q * (1.0 - rho);
        let de_drho = (1.0 + self.q) / (denom * denom) * (self.e_0 - self.e_min);
        let e_rho = self.young_modulus(rho);
        if e_rho.abs() < 1e-30 {
            return 0.0;
        }
        -de_drho / e_rho * element_compliance
    }
}
/// Heaviside threshold projection with continuation parameter β.
///
/// Maps filtered densities \tilde{ρ} to physical densities \bar{ρ}.
/// As β → ∞ the projection converges to a 0/1 Heaviside.
#[allow(dead_code)]
#[allow(non_snake_case)]
#[derive(Debug, Clone, Copy)]
pub struct DensityProjection {
    /// Projection sharpness (β); must be positive.
    pub beta: f64,
    /// Threshold parameter (η); typically 0.5.
    pub eta: f64,
}
impl DensityProjection {
    /// Create a new density projection.
    pub fn new(beta: f64, eta: f64) -> Self {
        assert!(beta > 0.0, "beta must be positive");
        assert!(eta > 0.0 && eta < 1.0, "eta must be in (0,1)");
        Self { beta, eta }
    }
    /// Project filtered density `rho_tilde` to physical density.
    #[allow(non_snake_case)]
    pub fn project(&self, rho_tilde: f64) -> f64 {
        let β = self.beta;
        let η = self.eta;
        let tanh_b_eta = (β * η).tanh();
        let tanh_b_1_eta = (β * (1.0 - η)).tanh();
        let num = (tanh_b_eta + (β * (rho_tilde - η)).tanh()).max(0.0);
        let denom = tanh_b_eta + tanh_b_1_eta;
        if denom < 1e-30 {
            return rho_tilde;
        }
        num / denom
    }
    /// Derivative of projection w.r.t. `rho_tilde`.
    #[allow(non_snake_case)]
    pub fn derivative(&self, rho_tilde: f64) -> f64 {
        let β = self.beta;
        let η = self.eta;
        let tanh_b_eta = (β * η).tanh();
        let tanh_b_1_eta = (β * (1.0 - η)).tanh();
        let denom = tanh_b_eta + tanh_b_1_eta;
        if denom < 1e-30 {
            return 1.0;
        }
        let sech2 = 1.0 - (β * (rho_tilde - η)).tanh().powi(2);
        β * sech2 / denom
    }
    /// Increase β for continuation (return new projection).
    pub fn continuation_step(&self, factor: f64, beta_max: f64) -> Self {
        Self {
            beta: (self.beta * factor).min(beta_max),
            eta: self.eta,
        }
    }
}
/// Drives the OC topology optimization loop.
pub struct OcOptimizer {
    /// SIMP parameters.
    pub params: SimpParams,
    /// Current design grid.
    pub grid: TopologyGrid,
    /// Completed iteration count.
    pub iteration: u32,
    /// Compliance history (one value per iteration).
    pub compliance_history: Vec<f64>,
    /// Volume fraction history.
    pub volume_history: Vec<f64>,
}
impl OcOptimizer {
    /// Construct a new optimizer with a uniform initial density of `params.volume_fraction`.
    pub fn new(nx: usize, ny: usize, params: SimpParams) -> Self {
        let rho0 = params.volume_fraction;
        let e0 = params.e0;
        let grid = TopologyGrid::new_uniform(nx, ny, rho0, e0, 1.0, 1.0);
        Self {
            params,
            grid,
            iteration: 0,
            compliance_history: Vec::new(),
            volume_history: Vec::new(),
        }
    }
    /// Perform one OC iteration (filter → update densities).
    ///
    /// In a full FEA-coupled loop the caller would compute element strain energies
    /// and call `SimpElement::compliance_sensitivity` first; here we use the stored
    /// `sensitivity` field directly so the struct is FEA-agnostic.
    pub fn step(&mut self) {
        let radius = self.params.filter_radius;
        let filtered = self.grid.filter_sensitivities(radius);
        for (e, s) in self.grid.elements.iter_mut().zip(filtered.iter()) {
            e.sensitivity = *s;
        }
        let sens = filtered;
        self.grid.update_densities_oc(&sens, &self.params);
        self.compliance_history.push(self.grid.total_compliance());
        self.volume_history.push(self.grid.volume_fraction());
        self.iteration += 1;
    }
    /// Returns per-element density changes since last step.
    /// Returns an empty `Vec` when all changes are below `tol` (converged).
    pub fn has_converged(&self, tol: f64) -> Vec<f64> {
        let changes: Vec<f64> = self
            .grid
            .elements
            .iter()
            .map(|e| (e.rho - self.params.volume_fraction).abs())
            .collect();
        if changes.iter().all(|&c| c < tol) {
            Vec::new()
        } else {
            changes
        }
    }
    /// Current flat density field (row-major).
    pub fn density_field(&self) -> Vec<f64> {
        self.grid.elements.iter().map(|e| e.rho).collect()
    }
    /// Compliance change between the last two iterations.
    /// Returns `None` if fewer than 2 iterations have been performed.
    pub fn compliance_change(&self) -> Option<f64> {
        if self.compliance_history.len() < 2 {
            return None;
        }
        let n = self.compliance_history.len();
        Some((self.compliance_history[n - 1] - self.compliance_history[n - 2]).abs())
    }
}
/// Monitors density change convergence across iterations.
///
/// Declares convergence when the last `window` recorded changes are all
/// below `tolerance`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TopOptConvergenceMonitor {
    /// Number of consecutive iterations required below tolerance.
    pub window: usize,
    /// Convergence tolerance on the maximum density change.
    pub tolerance: f64,
    /// History of recorded changes (most recent at the end).
    pub(super) history: Vec<f64>,
}
impl TopOptConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(window: usize, tolerance: f64) -> Self {
        Self {
            window,
            tolerance,
            history: Vec::new(),
        }
    }
    /// Record the maximum density change for this iteration.
    pub fn record_change(&mut self, max_change: f64) {
        self.history.push(max_change);
    }
    /// Returns `true` when the last `window` changes are all < `tolerance`.
    pub fn is_converged(&self) -> bool {
        if self.history.len() < self.window {
            return false;
        }
        self.history
            .iter()
            .rev()
            .take(self.window)
            .all(|&c| c < self.tolerance)
    }
    /// Maximum density change recorded over all iterations.
    pub fn max_change_overall(&self) -> f64 {
        self.history
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Number of iterations recorded.
    pub fn iteration_count(&self) -> usize {
        self.history.len()
    }
    /// Reset the monitor.
    pub fn reset(&mut self) {
        self.history.clear();
    }
}
/// A full level-set topology optimizer coupling the level-set field with the
/// Hamilton-Jacobi advection.
pub struct LevelSetOptimizer {
    /// The level-set field.
    pub field: LevelSetField,
    /// Pseudo-time step size.
    pub dt: f64,
    /// Interface regularization width.
    pub eps: f64,
    /// Volume fraction target.
    pub target_vf: f64,
    /// Iteration count.
    pub iteration: u32,
}
impl LevelSetOptimizer {
    /// Create a new level-set optimizer.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, target_vf: f64, dt: f64, eps: f64) -> Self {
        let field = LevelSetField::new_uniform(nx, ny, dx, dy, target_vf);
        Self {
            field,
            dt,
            eps,
            target_vf,
            iteration: 0,
        }
    }
    /// Perform one optimization step.
    ///
    /// `sensitivities` should be compliance sensitivities ∂c/∂ρ.
    pub fn step(&mut self, sensitivities: &[f64]) {
        let vel = self.field.sensitivity_velocity(sensitivities, self.eps);
        self.field.advect(&vel, self.dt);
        self.enforce_volume_constraint();
        self.iteration += 1;
    }
    /// Enforce the volume constraint by shifting phi uniformly.
    fn enforce_volume_constraint(&mut self) {
        let current_vf = self.field.volume_fraction();
        let error = current_vf - self.target_vf;
        let shift = -error * 0.1;
        for phi in self.field.phi.iter_mut() {
            *phi += shift;
        }
    }
    /// Return the current material indicator vector.
    pub fn material_indicator(&self) -> Vec<f64> {
        self.field.material_indicator(self.eps)
    }
}
/// Optimality Criteria (OC) density update solver.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OcSolver {
    /// Current element densities.
    pub densities: Vec<f64>,
    /// Lagrange multiplier for volume constraint.
    pub lagrange: f64,
    /// Maximum density move per iteration.
    pub move_limit: f64,
}
#[allow(dead_code)]
impl OcSolver {
    /// Create a new OC solver with uniform initial density.
    pub fn new(n_elements: usize, initial_density: f64, move_limit: f64) -> Self {
        Self {
            densities: vec![initial_density; n_elements],
            lagrange: 1.0,
            move_limit,
        }
    }
    /// Perform one OC density update step.
    ///
    /// `sensitivities[i]` = |dC/dρ_i| (positive sensitivity magnitude).
    /// `target_volume` = desired volume fraction in \[0,1\].
    ///
    /// Uses bisection to find the Lagrange multiplier that satisfies the
    /// volume constraint, then updates densities.
    pub fn update(&mut self, sensitivities: &[f64], target_volume: f64) {
        assert_eq!(sensitivities.len(), self.densities.len());
        let n = self.densities.len();
        let move_limit = self.move_limit;
        let mut lo = 1e-20_f64;
        let mut hi = 1e20_f64;
        let target_vol_total = target_volume * n as f64;
        for _ in 0..50 {
            let mid = (lo * hi).sqrt();
            let vol: f64 = self
                .densities
                .iter()
                .zip(sensitivities.iter())
                .map(|(&rho, &s)| {
                    let rho_new = rho * (s / mid).sqrt();
                    rho_new.clamp((rho - move_limit).max(1e-3), (rho + move_limit).min(1.0))
                })
                .sum();
            if vol > target_vol_total {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        self.lagrange = (lo * hi).sqrt();
        let lam = self.lagrange;
        let old = self.densities.clone();
        for (rho, &s) in self.densities.iter_mut().zip(sensitivities.iter()) {
            let rho_old = *rho;
            let rho_new = rho_old * (s / lam).sqrt();
            *rho = rho_new.clamp(
                (old[old.iter().position(|&x| x == rho_old).unwrap_or(0)] - move_limit).max(1e-3),
                (old[old.iter().position(|&x| x == rho_old).unwrap_or(0)] + move_limit).min(1.0),
            );
        }
    }
}
/// Multi-load case topology optimization using weighted compliance.
///
/// Given `n_loads` load cases, minimizes:
/// `C_total = sum_l weight_l * compliance_l`
///
/// The sensitivity is accumulated across all load cases.
#[derive(Debug, Clone)]
pub struct MultiLoadOptimizer {
    /// Weights for each load case.
    pub weights: Vec<f64>,
    /// OC optimizer per load case (shares the same grid via indices).
    pub params: SimpParams,
    /// Number of elements.
    pub n_elements: usize,
    /// Accumulated sensitivity across all load cases.
    pub accumulated_sensitivity: Vec<f64>,
}
impl MultiLoadOptimizer {
    /// Create a multi-load optimizer.
    pub fn new(n_loads: usize, weights: Vec<f64>, params: SimpParams, n_elements: usize) -> Self {
        assert_eq!(weights.len(), n_loads);
        let norm: f64 = weights.iter().sum();
        assert!(norm > 0.0, "Weights must sum to a positive value");
        Self {
            weights,
            params,
            n_elements,
            accumulated_sensitivity: vec![0.0; n_elements],
        }
    }
    /// Accumulate sensitivities from a single load case.
    ///
    /// `load_case_idx` identifies which load case.
    /// `sensitivities`  are the per-element sensitivities for this load case.
    pub fn accumulate_sensitivity(&mut self, load_case_idx: usize, sensitivities: &[f64]) {
        assert_eq!(sensitivities.len(), self.n_elements);
        let w = self.weights[load_case_idx];
        for (s_acc, &s) in self.accumulated_sensitivity.iter_mut().zip(sensitivities) {
            *s_acc += w * s;
        }
    }
    /// Reset accumulated sensitivities for a new iteration.
    pub fn reset_sensitivity(&mut self) {
        for s in self.accumulated_sensitivity.iter_mut() {
            *s = 0.0;
        }
    }
    /// Apply OC update using the accumulated (weighted) sensitivities.
    pub fn update_densities_oc(&self, grid: &mut TopologyGrid) {
        let sens = self.accumulated_sensitivity.clone();
        grid.update_densities_oc(&sens, &self.params);
    }
    /// Total weighted compliance given per-load-case compliances.
    pub fn total_compliance(&self, compliances: &[f64]) -> f64 {
        self.weights
            .iter()
            .zip(compliances.iter())
            .map(|(w, c)| w * c)
            .sum()
    }
}
/// A single 2-D rectangular element in the SIMP discretization.
#[derive(Debug, Clone)]
pub struct SimpElement {
    /// Current pseudo-density (0 = void, 1 = solid).
    pub rho: f64,
    /// Effective stiffness computed from `rho`.
    pub stiffness: f64,
    /// Raw compliance sensitivity ∂c/∂ρ.
    pub sensitivity: f64,
    /// Centroid x-coordinate.
    pub x0: f64,
    /// Centroid y-coordinate.
    pub y0: f64,
    /// Element width.
    pub dx: f64,
    /// Element height.
    pub dy: f64,
}
impl SimpElement {
    /// SIMP interpolation: E(ρ) = Emin + ρ^p · (E0 − Emin).
    pub fn effective_modulus(&self, params: &SimpParams) -> f64 {
        params.emin + self.rho.powf(params.penalty) * (params.e0 - params.emin)
    }
    /// Compliance sensitivity: ∂c/∂ρ = −p · ρ^(p−1) · (E0 − Emin) · u^T k_e u.
    ///
    /// `strain_energy` is the element strain energy u^T k_e u.
    pub fn compliance_sensitivity(&self, strain_energy: f64, params: &SimpParams) -> f64 {
        -params.penalty
            * self.rho.powf(params.penalty - 1.0)
            * (params.e0 - params.emin)
            * strain_energy
    }
    /// RAMP interpolation: E(ρ) = Emin + ρ/(1 + p*(1-ρ)) · (E0 − Emin).
    pub fn ramp_modulus(&self, params: &SimpParams) -> f64 {
        let denom = 1.0 + params.penalty * (1.0 - self.rho);
        params.emin + (self.rho / denom) * (params.e0 - params.emin)
    }
    /// RAMP sensitivity: ∂E/∂ρ = (1+p)/((1+p*(1-ρ))^2) · (E0 − Emin).
    pub fn ramp_sensitivity(&self, strain_energy: f64, params: &SimpParams) -> f64 {
        let p = params.penalty;
        let denom = 1.0 + p * (1.0 - self.rho);
        let de_drho = (1.0 + p) / (denom * denom) * (params.e0 - params.emin);
        -de_drho * strain_energy
    }
    /// Element area.
    pub fn area(&self) -> f64 {
        self.dx * self.dy
    }
}
/// Summary of a completed topology optimization run.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TopOptResult {
    /// Final density field (n_elements).
    pub densities: Vec<f64>,
    /// Final compliance objective value.
    pub compliance: f64,
    /// Final volume fraction.
    pub volume_fraction: f64,
    /// Iteration at which convergence was declared (0 = not converged).
    pub converged_iter: usize,
    /// Compliance history (one value per iteration).
    pub compliance_history: Vec<f64>,
    /// Density change history (max change per iteration).
    pub change_history: Vec<f64>,
}
impl TopOptResult {
    /// Create a result from optimizer state.
    pub fn new(
        densities: Vec<f64>,
        compliance: f64,
        volume_fraction: f64,
        converged_iter: usize,
        compliance_history: Vec<f64>,
        change_history: Vec<f64>,
    ) -> Self {
        Self {
            densities,
            compliance,
            volume_fraction,
            converged_iter,
            compliance_history,
            change_history,
        }
    }
    /// Discreteness measure: fraction of elements with density > 0.9 or < 0.1.
    pub fn binary_fraction(&self) -> f64 {
        let n = self.densities.len();
        if n == 0 {
            return 0.0;
        }
        let binary = self
            .densities
            .iter()
            .filter(|&&d| !(0.1..=0.9).contains(&d))
            .count();
        binary as f64 / n as f64
    }
    /// Returns `true` if the optimizer converged.
    pub fn converged(&self) -> bool {
        self.converged_iter > 0
    }
    /// Mean density over all elements.
    pub fn mean_density(&self) -> f64 {
        if self.densities.is_empty() {
            return 0.0;
        }
        self.densities.iter().sum::<f64>() / self.densities.len() as f64
    }
}
