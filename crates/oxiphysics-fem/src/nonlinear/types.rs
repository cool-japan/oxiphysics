//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Convergence history record for a Newton-Raphson solve.
#[derive(Debug, Clone)]
pub struct ConvergenceHistory {
    /// Residual norm at each iteration.
    pub residuals: Vec<f64>,
    /// Displacement increment norm at each iteration.
    pub du_norms: Vec<f64>,
}
impl ConvergenceHistory {
    /// Create an empty convergence history.
    pub fn new() -> Self {
        Self {
            residuals: Vec::new(),
            du_norms: Vec::new(),
        }
    }
    /// Record an iteration.
    pub fn record(&mut self, residual: f64, du_norm: f64) {
        self.residuals.push(residual);
        self.du_norms.push(du_norm);
    }
    /// Check whether the residuals are monotonically decreasing.
    pub fn is_monotone_decreasing(&self) -> bool {
        self.residuals.windows(2).all(|w| w[1] <= w[0])
    }
    /// Convergence rate estimate (ratio of last two residuals).
    pub fn convergence_rate(&self) -> Option<f64> {
        let n = self.residuals.len();
        if n < 2 || self.residuals[n - 2].abs() < 1e-60 {
            return None;
        }
        Some(self.residuals[n - 1] / self.residuals[n - 2])
    }
}
/// State for the Riks arc-length method.
#[allow(dead_code)]
pub struct RiksState {
    /// Current displacement vector.
    pub u: Vec<f64>,
    /// Current load factor.
    pub lambda: f64,
    /// Displacement increment from last step.
    pub du: Vec<f64>,
    /// Load factor increment from last step.
    pub dlambda: f64,
    /// Arc-length step size.
    pub arc_length: f64,
}
impl RiksState {
    /// Create the initial Riks state.
    #[allow(dead_code)]
    pub fn new(n_dof: usize, arc_length: f64) -> Self {
        Self {
            u: vec![0.0; n_dof],
            lambda: 0.0,
            du: vec![0.0; n_dof],
            dlambda: 0.01,
            arc_length,
        }
    }
    /// Compute the current arc length from increments.
    pub fn current_arc(&self) -> f64 {
        let du_sq: f64 = self.du.iter().map(|d| d * d).sum::<f64>();
        (du_sq + self.dlambda * self.dlambda).sqrt()
    }
}
/// Snap-through detector based on sign change of tangent stiffness determinant.
///
/// Monitors the sign of the determinant of the tangent stiffness matrix across
/// load steps. A sign change indicates a critical (buckling/snap-through) point.
#[allow(dead_code)]
pub struct SnapThroughDetector {
    /// Determinant signs from previous steps.
    pub det_signs: Vec<f64>,
    /// Load levels at which sign changes occurred.
    pub critical_loads: Vec<f64>,
}
impl SnapThroughDetector {
    /// Create a new detector.
    pub fn new() -> Self {
        Self {
            det_signs: Vec::new(),
            critical_loads: Vec::new(),
        }
    }
    /// Record a tangent stiffness matrix (2×2 dense) and load level.
    ///
    /// For a 2×2 matrix we can compute the determinant directly.
    /// Returns `true` if a snap-through was detected.
    pub fn check_2x2(&mut self, k: [[f64; 2]; 2], lambda: f64) -> bool {
        let det = k[0][0] * k[1][1] - k[0][1] * k[1][0];
        let sign = if det >= 0.0 { 1.0 } else { -1.0 };
        let snap = if let Some(&prev_sign) = self.det_signs.last() {
            prev_sign * sign < 0.0
        } else {
            false
        };
        if snap {
            self.critical_loads.push(lambda);
        }
        self.det_signs.push(sign);
        snap
    }
    /// Whether any snap-through has been detected.
    pub fn has_snap_through(&self) -> bool {
        !self.critical_loads.is_empty()
    }
}
/// Energy-based convergence criterion.
///
/// Convergence is declared when the incremental work `Δu · R` is less than
/// `energy_tol * (u · R_0)`, where `R_0` is the initial residual.
#[allow(dead_code)]
pub struct EnergyConvergenceCriteria {
    /// Maximum iterations.
    pub max_iter: usize,
    /// Energy ratio tolerance.
    pub energy_tol: f64,
}
impl EnergyConvergenceCriteria {
    /// Create new energy convergence criteria.
    pub fn new(max_iter: usize, energy_tol: f64) -> Self {
        Self {
            max_iter,
            energy_tol,
        }
    }
}
/// Load-step controller for incremental nonlinear analysis.
///
/// Applies load in a series of steps, adaptively adjusting step size based on
/// convergence history.
pub struct LoadStepper {
    /// Total load (scalar multiplier).
    pub total_load: f64,
    /// Current load level (λ in \[0, 1\]).
    pub current_lambda: f64,
    /// Step size for next load increment.
    pub step_size: f64,
    /// Minimum allowed step size.
    pub min_step: f64,
    /// Maximum allowed step size.
    pub max_step: f64,
    /// Target number of Newton iterations per step.
    pub target_iters: usize,
}
impl LoadStepper {
    /// Create a new load stepper with initial uniform step.
    pub fn new(total_load: f64, n_steps: usize) -> Self {
        let step = 1.0 / (n_steps.max(1) as f64);
        Self {
            total_load,
            current_lambda: 0.0,
            step_size: step,
            min_step: step * 0.1,
            max_step: step * 4.0,
            target_iters: 5,
        }
    }
    /// Advance to the next load level.
    ///
    /// Returns `None` when the full load has been applied.
    pub fn advance(&mut self) -> Option<f64> {
        if self.current_lambda >= 1.0 - 1e-14 {
            return None;
        }
        self.current_lambda = (self.current_lambda + self.step_size).min(1.0);
        Some(self.current_lambda * self.total_load)
    }
    /// Adapt the step size based on the number of Newton iterations required.
    ///
    /// If `actual_iters < target_iters`, increase the step; otherwise reduce.
    pub fn adapt_step(&mut self, actual_iters: usize) {
        let ratio = if actual_iters == 0 {
            2.0
        } else {
            (self.target_iters as f64 / actual_iters as f64).sqrt()
        };
        self.step_size = (self.step_size * ratio).clamp(self.min_step, self.max_step);
    }
    /// Check whether the full load has been applied.
    pub fn is_complete(&self) -> bool {
        self.current_lambda >= 1.0 - 1e-14
    }
}
/// Arc-length control for nonlinear buckling.
///
/// Simplified Riks method: controls arc length `Δs = sqrt(Δu·Δu + Δλ²)`.
pub struct ArcLengthControl {
    /// Arc-length step size `Δs`.
    pub arc_length: f64,
    /// Current load factor λ.
    pub lambda: f64,
    /// Load factor increment Δλ.
    pub d_lambda: f64,
}
impl ArcLengthControl {
    /// Create a new arc-length controller with the given step size.
    pub fn new(arc_length: f64) -> Self {
        Self {
            arc_length,
            lambda: 0.0,
            d_lambda: 0.0,
        }
    }
    /// Return the current arc-length step size.
    pub fn step_size(&self) -> f64 {
        self.arc_length
    }
}
/// Convergence criterion for Newton-Raphson.
#[derive(Debug, Clone)]
pub struct ConvergenceCriteria {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Residual tolerance: iteration stops when `||R|| < residual_tol`.
    pub residual_tol: f64,
    /// Displacement increment tolerance: iteration stops when `||Δu|| < displacement_tol`.
    pub displacement_tol: f64,
}
/// Result of a single continuation step.
#[allow(dead_code)]
pub struct ContinuationStepResult {
    /// Displacement at this load level.
    pub displacement: Vec<f64>,
    /// Load factor lambda at this step.
    pub lambda: f64,
    /// Whether Newton-Raphson converged at this step.
    pub converged: bool,
    /// Newton iterations used.
    pub iterations: usize,
}
/// Result of a Newton-Raphson solve.
pub struct NrResult {
    /// Computed displacement vector.
    pub displacement: Vec<f64>,
    /// Whether the solver converged within the allowed iterations.
    pub converged: bool,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Residual norm at final iteration.
    pub final_residual: f64,
}
