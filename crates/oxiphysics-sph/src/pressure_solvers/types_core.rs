//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Identifies which pressure solver variant to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverType {
    /// Weakly Compressible SPH (Tait EOS).
    Wcsph,
    /// Predictive-Corrective Incompressible SPH.
    Pcisph,
    /// Implicit Incompressible SPH (Jacobi PPE).
    Iisph,
}
/// Adaptive solver that switches between WCSPH and PCISPH based on
/// the current density error.
pub struct AdaptiveSolverSwitch {
    /// Threshold for switching from WCSPH to PCISPH.
    pub error_threshold: f64,
    /// Rest density.
    pub rho0: f64,
    /// Speed of sound.
    pub c0: f64,
    /// Tait exponent.
    pub gamma: f64,
    /// PCISPH delta.
    pub delta: f64,
    /// PCISPH max iterations.
    pub max_iter: usize,
    /// PCISPH tolerance.
    pub tolerance: f64,
    /// Currently selected solver.
    pub current_solver: SolverType,
    /// Number of switches performed.
    pub switch_count: usize,
}
impl AdaptiveSolverSwitch {
    /// Create a new adaptive solver.
    pub fn new(rho0: f64, c0: f64, error_threshold: f64) -> Self {
        Self {
            error_threshold,
            rho0,
            c0,
            gamma: 7.0,
            delta: 0.5,
            max_iter: 50,
            tolerance: 0.001,
            current_solver: SolverType::Wcsph,
            switch_count: 0,
        }
    }
    /// Compute pressures with automatic solver selection.
    ///
    /// If the max density error exceeds `error_threshold`, switches to
    /// PCISPH; otherwise uses WCSPH.
    pub fn compute_pressures(&mut self, densities: &[f64]) -> Vec<f64> {
        let max_err = densities
            .iter()
            .map(|&rho| ((rho - self.rho0) / self.rho0).abs())
            .fold(0.0_f64, f64::max);
        let new_solver = if max_err > self.error_threshold {
            SolverType::Pcisph
        } else {
            SolverType::Wcsph
        };
        if new_solver != self.current_solver {
            self.current_solver = new_solver;
            self.switch_count += 1;
        }
        match self.current_solver {
            SolverType::Wcsph => {
                let solver = WcSphSolver::new(self.rho0, self.c0, self.gamma);
                let mut pressures = vec![0.0; densities.len()];
                solver.compute_all_pressures(densities, &mut pressures);
                pressures
            }
            SolverType::Pcisph => {
                let solver = PcisphSolverSimple {
                    rho0: self.rho0,
                    delta: self.delta,
                    max_iter: self.max_iter,
                    tolerance: self.tolerance,
                };
                let (pressures, _) = solver.solve(densities);
                pressures
            }
            SolverType::Iisph => {
                let solver = WcSphSolver::new(self.rho0, self.c0, self.gamma);
                let mut pressures = vec![0.0; densities.len()];
                solver.compute_all_pressures(densities, &mut pressures);
                pressures
            }
        }
    }
}
/// Helper for building the IISPH pressure Poisson system.
///
/// Each particle contributes off-diagonal coefficients from the SPH kernel
/// gradient sum.
#[derive(Debug, Clone)]
pub struct IisphMatrixBuilder {
    /// Smoothing length.
    pub h: f64,
    /// Rest density.
    pub rho0: f64,
    /// Time step.
    pub dt: f64,
}
impl IisphMatrixBuilder {
    /// Create a matrix builder.
    pub fn new(h: f64, rho0: f64, dt: f64) -> Self {
        Self { h, rho0, dt }
    }
    /// Build the coefficient matrix (off-diagonal) for the IISPH PPE.
    ///
    /// `aij[i]` lists `(j, coeff)` off-diagonal pairs.
    /// The SPH kernel gradient is approximated as `r_ij / (|r_ij| * h²)`.
    pub fn build_matrix(
        &self,
        positions: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
    ) -> Vec<Vec<(usize, f64)>> {
        let n = positions.len();
        let h = self.h;
        let dt = self.dt;
        let mut aij: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        for i in 0..n {
            let rho_i = densities[i].max(1e-20);
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let grad_w_mag = 1.0 / (r * h * h);
                let rho_j = densities[j].max(1e-20);
                let coeff = -dt
                    * dt
                    * masses[j]
                    * (masses[i] / (rho_i * rho_i) + masses[j] / (rho_j * rho_j))
                    * grad_w_mag
                    * grad_w_mag;
                aij[i].push((j, coeff));
            }
        }
        aij
    }
    /// Compute the RHS `bᵢ = 1 - ρᵢ/ρ₀`.
    pub fn build_rhs(&self, densities: &[f64]) -> Vec<f64> {
        densities.iter().map(|&d| 1.0 - d / self.rho0).collect()
    }
    /// Diagonal coefficients `aᵢᵢ = -Σⱼ aᵢⱼ`.
    pub fn diagonal(&self, aij: &[Vec<(usize, f64)>]) -> Vec<f64> {
        aij.iter()
            .map(|row| -row.iter().map(|&(_, a)| a).sum::<f64>())
            .collect()
    }
}
/// Simplified PCISPH solver for the pressure_solvers comparison module.
///
/// Unlike the full [`crate::pcisph::PcisphSolver`], this version operates
/// on plain arrays and stores per-particle pressure data internally.
#[derive(Debug, Clone)]
pub struct PcisphSolverSimple {
    /// Rest density (kg/m³).
    pub rho0: f64,
    /// Pressure correction factor δ.
    pub delta: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance (fraction of rho0).
    pub tolerance: f64,
}
impl PcisphSolverSimple {
    /// Create a new simplified PCISPH solver.
    pub fn new(rho0: f64, delta: f64) -> Self {
        Self {
            rho0,
            delta,
            max_iter: 50,
            tolerance: 0.001,
        }
    }
    /// One prediction-correction iteration.
    ///
    /// Given predicted densities, updates pressures and returns the maximum
    /// absolute density error.
    pub fn iteration(&self, predicted_densities: &[f64], pressures: &mut [f64]) -> f64 {
        let mut max_err = 0.0_f64;
        for (i, &rho_pred) in predicted_densities.iter().enumerate() {
            let err = rho_pred - self.rho0;
            pressures[i] = (pressures[i] + self.delta * err).max(0.0);
            let abs_err = err.abs();
            if abs_err > max_err {
                max_err = abs_err;
            }
        }
        max_err
    }
    /// Run the convergence loop until the density error drops below tolerance
    /// or max_iter is reached.  Returns `(final_pressures, iterations_used)`.
    pub fn solve(&self, predicted_densities: &[f64]) -> (Vec<f64>, usize) {
        let n = predicted_densities.len();
        let mut pressures = vec![0.0_f64; n];
        let tol = self.tolerance * self.rho0;
        for iter in 0..self.max_iter {
            let max_err = self.iteration(predicted_densities, &mut pressures);
            if max_err < tol {
                return (pressures, iter + 1);
            }
        }
        (pressures, self.max_iter)
    }
}
/// Encapsulates the full PCISPH prediction-correction loop parameters and
/// provides a density-error convergence monitor with early exit.
#[derive(Debug, Clone)]
pub struct PcisphIterativeLoop {
    /// Rest density ρ₀ \[kg/m³\].
    pub rho0: f64,
    /// Pressure correction factor δ.
    pub delta: f64,
    /// Minimum iterations (avoids premature exit on first step).
    pub min_iter: usize,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence criterion: max |ρ_err| / ρ₀ < eta.
    pub eta: f64,
}
impl PcisphIterativeLoop {
    /// Create a loop controller with default min/max iterations.
    pub fn new(rho0: f64, delta: f64, eta: f64) -> Self {
        Self {
            rho0,
            delta,
            min_iter: 2,
            max_iter: 100,
            eta,
        }
    }
    /// Run the prediction-correction loop on an already-predicted density array.
    ///
    /// Returns `(pressures, iterations_used, diagnostics)`.
    pub fn solve_from_predicted_densities(
        &self,
        predicted_densities: &[f64],
    ) -> (Vec<f64>, usize, ConvergenceDiagnostics) {
        let n = predicted_densities.len();
        let mut pressures = vec![0.0_f64; n];
        let mut diag = ConvergenceDiagnostics::new();
        for iter in 0..self.max_iter {
            let mut max_err = 0.0_f64;
            let mut sum_err = 0.0_f64;
            for i in 0..n {
                let err = predicted_densities[i] - self.rho0;
                pressures[i] = (pressures[i] + self.delta * err).max(0.0);
                let ae = err.abs();
                if ae > max_err {
                    max_err = ae;
                }
                sum_err += ae;
            }
            let mean_err = if n > 0 { sum_err / n as f64 } else { 0.0 };
            let max_p = pressures.iter().cloned().fold(0.0_f64, f64::max);
            let converged = max_err < self.eta * self.rho0;
            diag.history.push(IterationDiagnostic {
                iteration: iter + 1,
                max_density_error: max_err,
                mean_density_error: mean_err,
                max_pressure: max_p,
                converged,
            });
            if converged && iter + 1 >= self.min_iter {
                return (pressures, iter + 1, diag);
            }
        }
        (pressures, self.max_iter, diag)
    }
    /// Compute max and mean density errors for a given predicted density array.
    pub fn density_error_stats(&self, predicted_densities: &[f64]) -> (f64, f64) {
        let n = predicted_densities.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        let errors: Vec<f64> = predicted_densities
            .iter()
            .map(|&r| (r - self.rho0).abs())
            .collect();
        let max_err = errors.iter().cloned().fold(0.0_f64, f64::max);
        let mean_err = errors.iter().sum::<f64>() / n as f64;
        (max_err, mean_err)
    }
}
/// Convenience wrapper for WCSPH pressure computation with a
/// pre-configured equation of state.
#[derive(Debug, Clone)]
pub struct WcSphSolver {
    /// Rest density (kg/m³).
    pub rho0: f64,
    /// Tait bulk modulus B.
    pub big_b: f64,
    /// Tait exponent.
    pub gamma: f64,
}
impl WcSphSolver {
    /// Create a solver from rest density and sound speed.
    pub fn new(rho0: f64, c0: f64, gamma: f64) -> Self {
        Self {
            rho0,
            big_b: tait_eos_b(rho0, c0, gamma),
            gamma,
        }
    }
    /// Compute Tait pressure from a given density.
    pub fn compute_tait_pressure(density: f64, rest_density: f64, b: f64, gamma: f64) -> f64 {
        b * ((density / rest_density).powf(gamma) - 1.0)
    }
    /// Compute pressure for all particles in-place.
    pub fn compute_all_pressures(&self, densities: &[f64], pressures: &mut [f64]) {
        for (i, &rho) in densities.iter().enumerate() {
            pressures[i] = Self::compute_tait_pressure(rho, self.rho0, self.big_b, self.gamma);
        }
    }
    /// Pressure for a single density value using this solver's parameters.
    pub fn pressure(&self, rho: f64) -> f64 {
        Self::compute_tait_pressure(rho, self.rho0, self.big_b, self.gamma)
    }
    /// Speed of sound at a given density: `c = c0 * (rho/rho0)^((gamma-1)/2)`.
    pub fn sound_speed(&self, rho: f64) -> f64 {
        let c0 = (self.big_b * self.gamma / self.rho0).sqrt();
        c0 * (rho / self.rho0).powf((self.gamma - 1.0) / 2.0)
    }
}
/// Results from a solver benchmark run.
#[derive(Debug, Clone)]
pub struct SolverBenchmarkResult {
    /// Which solver was used.
    pub solver_type: SolverType,
    /// Number of iterations to converge (0 for explicit solvers).
    pub iterations: usize,
    /// Final maximum density error (absolute).
    pub max_density_error: f64,
    /// Final average density error.
    pub avg_density_error: f64,
    /// Resulting pressures.
    pub pressures: Vec<f64>,
}
/// Diagnostic record for one pressure-solver iteration.
#[derive(Debug, Clone)]
pub struct IterationDiagnostic {
    /// Iteration number (1-based).
    pub iteration: usize,
    /// Maximum absolute density error across all particles \[kg/m³\].
    pub max_density_error: f64,
    /// Mean absolute density error \[kg/m³\].
    pub mean_density_error: f64,
    /// Maximum absolute pressure \[Pa\].
    pub max_pressure: f64,
    /// Whether convergence was reached at this iteration.
    pub converged: bool,
}
/// Accumulates per-iteration convergence data for pressure solvers.
#[derive(Debug, Clone, Default)]
pub struct ConvergenceDiagnostics {
    /// History of all iterations recorded.
    pub history: Vec<IterationDiagnostic>,
}
impl ConvergenceDiagnostics {
    /// Create a new empty diagnostics collector.
    pub fn new() -> Self {
        Self {
            history: Vec::new(),
        }
    }
    /// Record one iteration.
    pub fn record(
        &mut self,
        iteration: usize,
        densities: &[f64],
        pressures: &[f64],
        rho0: f64,
        tol: f64,
    ) {
        let n = densities.len();
        if n == 0 {
            return;
        }
        let errors: Vec<f64> = densities.iter().map(|&rho| (rho - rho0).abs()).collect();
        let max_err = errors.iter().cloned().fold(0.0_f64, f64::max);
        let mean_err = errors.iter().sum::<f64>() / n as f64;
        let max_p = pressures
            .iter()
            .cloned()
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        let converged = max_err < tol * rho0;
        self.history.push(IterationDiagnostic {
            iteration,
            max_density_error: max_err,
            mean_density_error: mean_err,
            max_pressure: max_p,
            converged,
        });
    }
    /// Return the iteration at which convergence was first achieved, if any.
    pub fn convergence_iteration(&self) -> Option<usize> {
        self.history
            .iter()
            .find(|d| d.converged)
            .map(|d| d.iteration)
    }
    /// Maximum density error ever seen across all iterations.
    pub fn worst_density_error(&self) -> f64 {
        self.history
            .iter()
            .map(|d| d.max_density_error)
            .fold(0.0_f64, f64::max)
    }
    /// Number of iterations recorded.
    pub fn len(&self) -> usize {
        self.history.len()
    }
    /// Whether no iterations have been recorded.
    pub fn is_empty(&self) -> bool {
        self.history.is_empty()
    }
    /// Final (last recorded) max density error.
    pub fn final_max_error(&self) -> f64 {
        self.history
            .last()
            .map(|d| d.max_density_error)
            .unwrap_or(0.0)
    }
}
/// Utility struct wrapping multi-step Jacobi pressure iteration.
///
/// Executes up to `max_iter` relaxed-Jacobi steps on the pressure Poisson
/// system `A p = rhs` and records convergence diagnostics.
#[derive(Debug, Clone)]
pub struct JacobiPressureIterator {
    /// Relaxation factor ω (default 0.7).
    pub omega: f64,
    /// Maximum number of steps.
    pub max_iter: usize,
    /// Convergence tolerance (absolute residual).
    pub tolerance: f64,
    /// Rest density (for recording diagnostics).
    pub rho0: f64,
}
impl JacobiPressureIterator {
    /// Create a new iterator with given parameters.
    pub fn new(rho0: f64, omega: f64, max_iter: usize, tolerance: f64) -> Self {
        Self {
            omega,
            max_iter,
            tolerance,
            rho0,
        }
    }
    /// Run Jacobi iterations and return `(final_pressures, residual, iters)`.
    ///
    /// The system matrix is given in sparse CSR-like form: `aij[i]` is a list
    /// of `(j, a_ij)` off-diagonal pairs.  Diagonal is inferred as `-Σ a_ij`.
    pub fn run(
        &self,
        rhs: &[f64],
        aij: &[Vec<(usize, f64)>],
        mut diagnostics: Option<&mut ConvergenceDiagnostics>,
    ) -> (Vec<f64>, f64, usize) {
        let n = rhs.len();
        let mut p = vec![0.0_f64; n];
        let mut residual = 0.0_f64;
        for iter in 0..self.max_iter {
            let mut p_new = vec![0.0_f64; n];
            residual = 0.0;
            for i in 0..n {
                let aii: f64 = -aij[i].iter().map(|&(_, a)| a).sum::<f64>();
                if aii.abs() < 1e-20 {
                    p_new[i] = p[i];
                    continue;
                }
                let off: f64 = aij[i].iter().map(|&(j, a)| a * p[j]).sum();
                let target = (rhs[i] - off) / aii;
                p_new[i] = ((1.0 - self.omega) * p[i] + self.omega * target).max(0.0);
                let res_i = (aii * p_new[i] + off - rhs[i]).abs();
                if res_i > residual {
                    residual = res_i;
                }
            }
            p = p_new;
            if let Some(ref mut diag) = diagnostics.as_deref_mut() {
                let densities: Vec<f64> = rhs.iter().map(|&r| self.rho0 + r).collect();
                diag.record(iter + 1, &densities, &p, self.rho0, self.tolerance);
            }
            if residual < self.tolerance {
                return (p, residual, iter + 1);
            }
        }
        (p, residual, self.max_iter)
    }
}
/// Integrate the SPH continuity equation forward in time.
///
/// The density evolution follows `Dρᵢ/Dt = Σⱼ mⱼ (vᵢ - vⱼ) · ∇Wᵢⱼ`.
#[derive(Debug, Clone)]
pub struct SphDensityIntegrator {
    /// Rest density ρ₀ \[kg/m³\].
    pub rho0: f64,
    /// Smoothing length h \[m\].
    pub h: f64,
}
impl SphDensityIntegrator {
    /// Create a new integrator.
    pub fn new(rho0: f64, h: f64) -> Self {
        Self { rho0, h }
    }
    /// Compute the density rate-of-change `Dρᵢ/Dt` for particle `i`.
    ///
    /// Uses an approximate cubic-spline gradient `-r_ij / (h * |r_ij|²)`.
    pub fn drho_dt(
        &self,
        i: usize,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
    ) -> f64 {
        let n = positions.len();
        let pi = positions[i];
        let vi = velocities[i];
        let mut drho = 0.0_f64;
        for j in 0..n {
            if j == i {
                continue;
            }
            let r = [
                pi[0] - positions[j][0],
                pi[1] - positions[j][1],
                pi[2] - positions[j][2],
            ];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            let r_len = r2.sqrt();
            if r_len < 1e-14 || r_len > 2.0 * self.h {
                continue;
            }
            let dv = [
                vi[0] - velocities[j][0],
                vi[1] - velocities[j][1],
                vi[2] - velocities[j][2],
            ];
            let inv_hr = 1.0 / (self.h * r2);
            let dv_dot_r = dv[0] * r[0] + dv[1] * r[1] + dv[2] * r[2];
            drho += masses[j] * dv_dot_r * inv_hr;
        }
        drho
    }
    /// Forward-Euler density update: `ρᵢ(t+dt) = ρᵢ(t) + (Dρᵢ/Dt) * dt`.
    pub fn integrate(
        &self,
        densities: &mut [f64],
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        dt: f64,
    ) {
        let n = densities.len();
        let rates: Vec<f64> = (0..n)
            .map(|i| self.drho_dt(i, positions, velocities, masses))
            .collect();
        for i in 0..n {
            densities[i] = (densities[i] + rates[i] * dt).max(0.0);
        }
    }
    /// Compute density errors relative to `rho0`.
    pub fn density_errors(&self, densities: &[f64]) -> Vec<f64> {
        densities.iter().map(|&d| d - self.rho0).collect()
    }
    /// Maximum absolute density error.
    pub fn max_density_error(&self, densities: &[f64]) -> f64 {
        densities
            .iter()
            .map(|&d| (d - self.rho0).abs())
            .fold(0.0_f64, f64::max)
    }
}
/// Applies pressure-gradient corrections to predicted velocities and positions
/// at the end of the PCISPH loop.
#[derive(Debug, Clone)]
pub struct PcisphCorrector {
    /// Time step.
    pub dt: f64,
    /// Rest density.
    pub rho0: f64,
}
impl PcisphCorrector {
    /// Create a corrector.
    pub fn new(dt: f64, rho0: f64) -> Self {
        Self { dt, rho0 }
    }
    /// Apply pressure-gradient corrections to predicted velocities.
    ///
    /// `v_corr_i = v_pred_i - dt/m_i * Σⱼ m_j (p_i/ρᵢ² + p_j/ρⱼ²) ∇Wᵢⱼ`
    ///
    /// The kernel gradient is approximated as `r_ij / (h * |r_ij|²)`.
    pub fn correct_velocities(
        &self,
        positions: &[[f64; 3]],
        pred_velocities: &[[f64; 3]],
        pressures: &[f64],
        masses: &[f64],
        densities: &[f64],
        h: f64,
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut v_corr = pred_velocities.to_vec();
        for i in 0..n {
            let rho_i = densities[i].max(1e-20);
            let p_i = pressures[i];
            let m_i = masses[i];
            let mut dv = [0.0_f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = [
                    positions[i][0] - positions[j][0],
                    positions[i][1] - positions[j][1],
                    positions[i][2] - positions[j][2],
                ];
                let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                let r_len = r2.sqrt();
                if r_len < 1e-14 || r_len > 2.0 * h {
                    continue;
                }
                let rho_j = densities[j].max(1e-20);
                let p_j = pressures[j];
                let grad_w = [-r[0] / (h * r2), -r[1] / (h * r2), -r[2] / (h * r2)];
                let coeff = masses[j] * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j));
                dv[0] += coeff * grad_w[0];
                dv[1] += coeff * grad_w[1];
                dv[2] += coeff * grad_w[2];
            }
            v_corr[i][0] -= self.dt / m_i * dv[0];
            v_corr[i][1] -= self.dt / m_i * dv[1];
            v_corr[i][2] -= self.dt / m_i * dv[2];
        }
        v_corr
    }
    /// Apply position corrections: `x_corr_i = x_pred_i + v_corr_i * dt`.
    pub fn correct_positions(
        &self,
        pred_positions: &[[f64; 3]],
        v_corr: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        pred_positions
            .iter()
            .zip(v_corr.iter())
            .map(|(x, v)| {
                [
                    x[0] + v[0] * self.dt,
                    x[1] + v[1] * self.dt,
                    x[2] + v[2] * self.dt,
                ]
            })
            .collect()
    }
    /// Verify position-correction energy: corrected KE should generally be
    /// less than predicted KE for positive pressures (pressure does negative
    /// work on approaching particles).
    pub fn ke_change(&self, masses: &[f64], v_pred: &[[f64; 3]], v_corr: &[[f64; 3]]) -> f64 {
        let ke_pred: f64 = masses
            .iter()
            .zip(v_pred.iter())
            .map(|(&m, v)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum();
        let ke_corr: f64 = masses
            .iter()
            .zip(v_corr.iter())
            .map(|(&m, v)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum();
        ke_corr - ke_pred
    }
}
/// IISPH solver (Ihmsen et al. 2014) with diagonal preconditioning.
///
/// Iteratively solves the pressure Poisson equation (PPE) using a relaxed
/// Jacobi scheme with diagonal (Jacobi) preconditioning.
///
/// This is distinct from [`IisphSolverSimple`] in that it exposes the diagonal
/// element computation and the main iterate method as separate `pub fn`s.
#[derive(Debug, Clone)]
pub struct IisphSolver {
    /// Rest density ρ₀ \[kg/m³\].
    pub rest_density: f64,
    /// Relaxation factor ω ∈ (0, 2).
    pub omega: f64,
    /// Maximum number of Jacobi iterations.
    pub max_iter: usize,
    /// Convergence tolerance (fraction of ρ₀).
    pub tolerance: f64,
}
impl IisphSolver {
    /// Create a new `IisphSolver` with given parameters.
    pub fn new(rest_density: f64, omega: f64, max_iter: usize, tolerance: f64) -> Self {
        Self {
            rest_density,
            omega,
            max_iter,
            tolerance,
        }
    }
    /// Compute the diagonal coefficient `a_ii` from the off-diagonal row entries.
    ///
    /// Standard Laplacian assembly: `a_ii = -Σ_j a_ij`.
    pub fn compute_aii(off_diag: &[(usize, f64)]) -> f64 {
        -off_diag.iter().map(|&(_, a)| a).sum::<f64>()
    }
    /// One pass of the relaxed-Jacobi pressure iteration.
    ///
    /// # Arguments
    /// * `rhs`  – Right-hand side of the PPE (density error / dt²).
    /// * `aij`  – Sparse off-diagonal coefficients for each row.
    ///
    /// # Returns
    /// `(pressures, iterations_used, final_residual)`
    pub fn iterate_pressure(
        &self,
        rhs: &[f64],
        aij: &[Vec<(usize, f64)>],
    ) -> (Vec<f64>, usize, f64) {
        let n = rhs.len();
        let mut p = vec![0.0_f64; n];
        let mut residual = 0.0_f64;
        for iter in 0..self.max_iter {
            let mut p_new = vec![0.0_f64; n];
            residual = 0.0;
            for i in 0..n {
                let aii = Self::compute_aii(&aij[i]);
                if aii.abs() < 1e-20 {
                    p_new[i] = p[i];
                    continue;
                }
                let off: f64 = aij[i].iter().map(|&(j, a)| a * p[j]).sum();
                let target = (rhs[i] - off) / aii;
                p_new[i] = ((1.0 - self.omega) * p[i] + self.omega * target).max(0.0);
                let res_i = (aii * p_new[i] + off - rhs[i]).abs();
                if res_i > residual {
                    residual = res_i;
                }
            }
            p = p_new;
            if residual < self.tolerance * self.rest_density {
                return (p, iter + 1, residual);
            }
        }
        (p, self.max_iter, residual)
    }
}
impl IisphSolver {
    /// Check whether the velocity field satisfies the divergence-free constraint.
    ///
    /// The IISPH divergence-free condition requires that for all particles:
    ///   `|D_i| = |sum_j m_j (v_i - v_j) · ∇W_ij / ρ_i| < tol`
    ///
    /// This function computes the maximum absolute divergence across all particles
    /// and returns `true` if it is below `tolerance`.
    ///
    /// # Arguments
    /// * `velocities` – Current velocities.
    /// * `masses`     – Particle masses.
    /// * `densities`  – Current densities.
    /// * `neighbors`  – For each particle i, list of `(j, r_ij, grad_W_ij)`.
    /// * `tolerance`  – Divergence threshold \[1/s\].
    pub fn compute_divergence_free_constraint(
        &self,
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<(usize, f64, [f64; 3])>],
        tolerance: f64,
    ) -> bool {
        let n = velocities.len();
        if n == 0 {
            return true;
        }
        for i in 0..n {
            let rho_i = densities[i].max(1e-14);
            let mut div = 0.0_f64;
            for &(j, _r_ij, grad_w) in &neighbors[i] {
                let mj = masses[j];
                let dv = [
                    velocities[i][0] - velocities[j][0],
                    velocities[i][1] - velocities[j][1],
                    velocities[i][2] - velocities[j][2],
                ];
                div += mj / rho_i * (dv[0] * grad_w[0] + dv[1] * grad_w[1] + dv[2] * grad_w[2]);
            }
            if div.abs() >= tolerance {
                return false;
            }
        }
        true
    }
}
/// Projects the velocity field onto the divergence-free subspace using the
/// DFSPH approach: correct velocity until `|∇·v| < tol`.
#[derive(Debug, Clone)]
pub struct DivergenceFreeProjector {
    /// Rest density.
    pub rho0: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Divergence tolerance.
    pub tolerance: f64,
    /// Relaxation factor.
    pub omega: f64,
}
impl DivergenceFreeProjector {
    /// Create a new projector.
    pub fn new(rho0: f64, max_iter: usize, tolerance: f64, omega: f64) -> Self {
        Self {
            rho0,
            max_iter,
            tolerance,
            omega,
        }
    }
    /// Compute the velocity divergence at particle `i` using SPH sum.
    ///
    /// `div_i = Σⱼ (m_j / ρ_j) (v_i - v_j) · ∇W_ij`
    pub fn divergence_at(
        &self,
        i: usize,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        h: f64,
    ) -> f64 {
        let n = positions.len();
        let pi = positions[i];
        let vi = velocities[i];
        let mut div = 0.0_f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let r = [
                pi[0] - positions[j][0],
                pi[1] - positions[j][1],
                pi[2] - positions[j][2],
            ];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            let r_len = r2.sqrt();
            if r_len < 1e-14 || r_len > 2.0 * h {
                continue;
            }
            let dv = [
                vi[0] - velocities[j][0],
                vi[1] - velocities[j][1],
                vi[2] - velocities[j][2],
            ];
            let grad_w = [-r[0] / (h * r2), -r[1] / (h * r2), -r[2] / (h * r2)];
            let m_rho = masses[j] / densities[j].max(1e-20);
            div += m_rho * (dv[0] * grad_w[0] + dv[1] * grad_w[1] + dv[2] * grad_w[2]);
        }
        div
    }
    /// Check whether the velocity field is divergence-free within `tolerance`.
    pub fn is_divergence_free(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        h: f64,
    ) -> bool {
        let n = positions.len();
        for i in 0..n {
            let d = self
                .divergence_at(i, positions, velocities, masses, densities, h)
                .abs();
            if d > self.tolerance {
                return false;
            }
        }
        true
    }
    /// Max absolute divergence across all particles.
    pub fn max_divergence(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        h: f64,
    ) -> f64 {
        let n = positions.len();
        (0..n)
            .map(|i| {
                self.divergence_at(i, positions, velocities, masses, densities, h)
                    .abs()
            })
            .fold(0.0_f64, f64::max)
    }
}
/// Per-particle state used during the PCISPH prediction-correction loop.
///
/// Holds predicted (intermediate) quantities updated in each iteration, as well
/// as the accumulated pressure corrections.
#[derive(Debug, Clone)]
pub struct PcisphState {
    /// Predicted positions after the current sub-step \[m\].
    pub predicted_positions: Vec<[f64; 3]>,
    /// Predicted velocities \[m/s\].
    pub predicted_velocities: Vec<[f64; 3]>,
    /// Predicted densities \[kg/m³\] computed from predicted positions.
    pub predicted_densities: Vec<f64>,
    /// Per-particle pressure correction accumulated this iteration \[Pa\].
    pub pressure_corrections: Vec<f64>,
}
impl PcisphState {
    /// Create a zero-initialised state for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            predicted_positions: vec![[0.0; 3]; n],
            predicted_velocities: vec![[0.0; 3]; n],
            predicted_densities: vec![0.0; n],
            pressure_corrections: vec![0.0; n],
        }
    }
    /// Reset all pressure corrections to zero for the next iteration.
    pub fn reset_corrections(&mut self) {
        for c in &mut self.pressure_corrections {
            *c = 0.0;
        }
    }
    /// Number of particles tracked.
    pub fn len(&self) -> usize {
        self.predicted_positions.len()
    }
    /// Returns `true` if no particles are tracked.
    pub fn is_empty(&self) -> bool {
        self.predicted_positions.is_empty()
    }
}
/// Simplified IISPH solver for the pressure_solvers comparison module.
///
/// Operates on a pre-assembled PPE system (diagonal + off-diagonal).
#[derive(Debug, Clone)]
pub struct IisphSolverSimple {
    /// Rest density (kg/m³).
    pub rho0: f64,
    /// Relaxation factor (default 0.5).
    pub omega: f64,
    /// Maximum iterations (default 100).
    pub max_iter: usize,
    /// Convergence tolerance as fraction of rho0 (default 0.001).
    pub tolerance: f64,
}
impl IisphSolverSimple {
    /// Create a new simplified IISPH solver.
    pub fn new(rho0: f64) -> Self {
        Self {
            rho0,
            omega: 0.5,
            max_iter: 100,
            tolerance: 0.001,
        }
    }
    /// Compute the diagonal element `a_ii` from off-diagonal coefficients.
    ///
    /// In the standard Laplacian discretisation: `a_ii = -Σ_j a_ij`.
    pub fn compute_aii(off_diag: &[(usize, f64)]) -> f64 {
        -off_diag.iter().map(|&(_, a)| a).sum::<f64>()
    }
    /// Solve the pressure Poisson equation using relaxed Jacobi.
    ///
    /// # Arguments
    /// * `rhs` – right-hand side vector (density error / dt²).
    /// * `aij` – sparse off-diagonal coefficients for each row.
    ///
    /// # Returns
    /// `(pressures, iterations_used, final_residual)`
    pub fn solve_pressure_poisson(
        &self,
        rhs: &[f64],
        aij: &[Vec<(usize, f64)>],
    ) -> (Vec<f64>, usize, f64) {
        let n = rhs.len();
        let mut p = vec![0.0_f64; n];
        let mut residual = 0.0_f64;
        for iter in 0..self.max_iter {
            let mut p_new = vec![0.0_f64; n];
            residual = 0.0;
            for i in 0..n {
                let aii = Self::compute_aii(&aij[i]);
                if aii.abs() < 1e-20 {
                    p_new[i] = p[i];
                    continue;
                }
                let off: f64 = aij[i].iter().map(|&(j, a)| a * p[j]).sum();
                let target = (rhs[i] - off) / aii;
                p_new[i] = ((1.0 - self.omega) * p[i] + self.omega * target).max(0.0);
                let res_i = (aii * p_new[i] + off - rhs[i]).abs();
                if res_i > residual {
                    residual = res_i;
                }
            }
            p = p_new;
            if residual < self.tolerance * self.rho0 {
                return (p, iter + 1, residual);
            }
        }
        (p, self.max_iter, residual)
    }
}
/// Successive-over-relaxation (SOR) solver for the SPH pressure Poisson
/// equation.
///
/// Solves `A p = b` where `A` is a Laplacian-like matrix encoded as sparse
/// off-diagonal coefficients.
#[derive(Debug, Clone)]
pub struct PressurePoissonSolver {
    /// Over-relaxation factor ω ∈ (0, 2).
    pub omega: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Absolute residual tolerance.
    pub tolerance: f64,
    /// Whether to clamp pressures ≥ 0.
    pub clamp_non_negative: bool,
}
impl PressurePoissonSolver {
    /// Create a solver with given parameters.
    pub fn new(omega: f64, max_iter: usize, tolerance: f64) -> Self {
        Self {
            omega,
            max_iter,
            tolerance,
            clamp_non_negative: true,
        }
    }
    /// Run SOR iterations.
    ///
    /// Returns `(pressures, residual, iterations)`.
    pub fn solve(&self, rhs: &[f64], aij: &[Vec<(usize, f64)>]) -> (Vec<f64>, f64, usize) {
        let n = rhs.len();
        if n == 0 {
            return (vec![], 0.0, 0);
        }
        let mut p = vec![0.0_f64; n];
        let mut residual = 0.0_f64;
        for iter in 0..self.max_iter {
            residual = 0.0;
            for i in 0..n {
                let aii: f64 = -aij[i].iter().map(|&(_, a)| a).sum::<f64>();
                if aii.abs() < 1e-20 {
                    continue;
                }
                let off: f64 = aij[i].iter().map(|&(j, a)| a * p[j]).sum();
                let p_gs = (rhs[i] - off) / aii;
                let p_new = (1.0 - self.omega) * p[i] + self.omega * p_gs;
                let p_new = if self.clamp_non_negative {
                    p_new.max(0.0)
                } else {
                    p_new
                };
                let res_i = (aii * p_new + off - rhs[i]).abs();
                if res_i > residual {
                    residual = res_i;
                }
                p[i] = p_new;
            }
            if residual < self.tolerance {
                return (p, residual, iter + 1);
            }
        }
        (p, residual, self.max_iter)
    }
    /// Compute the residual `|A p - b|_max` for a given pressure vector.
    pub fn residual(&self, p: &[f64], rhs: &[f64], aij: &[Vec<(usize, f64)>]) -> f64 {
        let n = p.len();
        let mut max_res = 0.0_f64;
        for i in 0..n {
            let aii: f64 = -aij[i].iter().map(|&(_, a)| a).sum::<f64>();
            let off: f64 = aij[i].iter().map(|&(j, a)| a * p[j]).sum();
            let res_i = (aii * p[i] + off - rhs[i]).abs();
            if res_i > max_res {
                max_res = res_i;
            }
        }
        max_res
    }
}
