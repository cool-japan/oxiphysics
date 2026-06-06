//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Detected bifurcation point.
#[derive(Debug, Clone)]
pub struct BifurcationPoint {
    /// State at bifurcation
    pub x: f64,
    /// Parameter at bifurcation
    pub lambda: f64,
    /// Bifurcation type
    pub bif_type: BifurcationType,
}
/// Bifurcation type.
#[derive(Debug, Clone, PartialEq)]
pub enum BifurcationType {
    /// Saddle-node (fold) bifurcation
    SaddleNode,
    /// Pitchfork bifurcation
    Pitchfork,
    /// Transcritical bifurcation
    Transcritical,
    /// Hopf bifurcation
    Hopf,
}
/// Result of a root-finding computation.
#[derive(Debug, Clone)]
pub struct RootResult {
    /// Estimated root
    pub root: f64,
    /// Function value at root
    pub f_val: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Converged?
    pub converged: bool,
}
/// Dormand-Prince RK45 adaptive step integrator.
///
/// Returns (t, y) pairs with local error control.
#[derive(Debug, Clone)]
pub struct Rk45State {
    /// Current time
    pub t: f64,
    /// Current state
    pub y: Vec<f64>,
    /// Step size
    pub h: f64,
    /// Tolerance
    pub rtol: f64,
    /// Absolute tolerance
    pub atol: f64,
}
