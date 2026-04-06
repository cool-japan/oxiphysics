//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Result of an LCP solve.
#[derive(Clone, Debug)]
pub struct LcpResult {
    /// Solution vector z.
    pub z: Vec<f64>,
    /// Slack vector w = M*z + q.
    pub w: Vec<f64>,
    /// Number of pivot steps.
    pub pivots: usize,
    /// Whether the solver found a solution.
    pub solved: bool,
}
/// Result of a barrier method solve.
#[derive(Clone, Debug)]
pub struct BarrierResult {
    /// Optimal primal variable.
    pub x: Vec<f64>,
    /// Objective value.
    pub objective: f64,
    /// Final barrier parameter.
    pub mu: f64,
    /// Number of outer iterations.
    pub outer_iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Result of a penalty method solve.
#[derive(Clone, Debug)]
pub struct PenaltyResult {
    /// Solution.
    pub x: Vec<f64>,
    /// Objective value.
    pub objective: f64,
    /// Final penalty parameter.
    pub penalty: f64,
    /// Constraint violation.
    pub violation: f64,
    /// Whether converged.
    pub converged: bool,
}
/// Result of a proximal gradient solve.
#[derive(Clone, Debug)]
pub struct ProxGradResult {
    /// Optimal variable.
    pub x: Vec<f64>,
    /// Objective value (smooth part).
    pub objective: f64,
    /// Number of iterations.
    pub iterations: usize,
    /// Whether converged.
    pub converged: bool,
}
/// Result of an ADMM solve.
#[derive(Clone, Debug)]
pub struct AdmmResult {
    /// Primal variable x.
    pub x: Vec<f64>,
    /// Consensus variable z.
    pub z: Vec<f64>,
    /// Dual variable y.
    pub y: Vec<f64>,
    /// Number of iterations.
    pub iterations: usize,
    /// Primal residual norm.
    pub primal_residual: f64,
    /// Dual residual norm.
    pub dual_residual: f64,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Result of a quadratic programming solve.
#[derive(Clone, Debug)]
pub struct QpResult {
    /// Optimal primal variable.
    pub x: Vec<f64>,
    /// Lagrange multipliers for inequality constraints.
    pub multipliers: Vec<f64>,
    /// Number of iterations used.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Result of dual decomposition.
#[derive(Clone, Debug)]
pub struct DualDecompResult {
    /// Local solutions for each subproblem.
    pub x_locals: Vec<Vec<f64>>,
    /// Consensus variable.
    pub z: Vec<f64>,
    /// Dual variables.
    pub dual: Vec<f64>,
    /// Number of iterations.
    pub iterations: usize,
    /// Whether converged.
    pub converged: bool,
}
/// A second-order cone constraint: ||A_i * x + b_i|| <= c_i^T * x + d_i.
#[derive(Clone, Debug)]
pub struct SocConstraint {
    /// Matrix A_i (k x n, row-major).
    pub a_mat: Vec<f64>,
    /// Offset vector b_i (length k).
    pub b_vec: Vec<f64>,
    /// Linear coefficient c_i (length n).
    pub c_vec: Vec<f64>,
    /// Scalar offset d_i.
    pub d_val: f64,
    /// Dimension of the cone (k).
    pub cone_dim: usize,
}
/// Result of an SOCP solve.
#[derive(Clone, Debug)]
pub struct SocpResult {
    /// Optimal primal variable.
    pub x: Vec<f64>,
    /// Objective value.
    pub objective: f64,
    /// Number of iterations.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Result of an augmented Lagrangian solve.
#[derive(Clone, Debug)]
pub struct AugLagResult {
    /// Optimal primal variable.
    pub x: Vec<f64>,
    /// Lagrange multiplier estimates.
    pub lambda: Vec<f64>,
    /// Final penalty parameter.
    pub rho: f64,
    /// Number of outer iterations.
    pub outer_iterations: usize,
    /// Constraint violation norm.
    pub violation: f64,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Result of a KKT conditions check.
#[derive(Clone, Debug)]
pub struct KktCheck {
    /// Stationarity violation (norm of gradient of Lagrangian).
    pub stationarity: f64,
    /// Primal feasibility violation.
    pub primal_feasibility: f64,
    /// Dual feasibility violation (min of multipliers for inequalities).
    pub dual_feasibility: f64,
    /// Complementary slackness violation.
    pub complementarity: f64,
    /// Whether KKT conditions are satisfied within tolerance.
    pub satisfied: bool,
}
