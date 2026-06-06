//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Result of a modal analysis.
pub struct ModalAnalysis {
    /// Natural frequencies in Hz (cycles per second), sorted ascending.
    pub frequencies_hz: Vec<f64>,
    /// Mode shapes (each is a unit vector of length n_dof), sorted by ascending frequency.
    pub mode_shapes: Vec<Vec<f64>>,
    /// Natural angular frequencies squared (omega^2), in (rad/s)^2.
    pub omega_sq: Vec<f64>,
}
/// Convergence diagnostics for subspace iteration.
pub struct SubspaceConvergence {
    /// Residuals per mode: ||A*v - lambda*v||.
    pub residuals: Vec<f64>,
}
