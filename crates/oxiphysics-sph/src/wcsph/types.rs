//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::{Real, Vec3};

/// Parameters for the WCSPH solver.
#[derive(Debug, Clone)]
pub struct WcsphParams {
    /// Rest density of the fluid (kg/m³).
    pub rest_density: f64,
    /// Stiffness parameter for Tait EOS.
    pub stiffness: f64,
    /// Exponent for Tait EOS (typically 7 for water).
    pub gamma: f64,
    /// Kinematic viscosity coefficient.
    pub viscosity: f64,
    /// Smoothing length.
    pub smoothing_length: f64,
}
/// WCSPH solver state, providing per-step utilities such as adaptive timestepping.
#[derive(Debug, Clone)]
pub struct WcsphSolver {
    /// Smoothing length used for this solver instance.
    pub smoothing_length: f64,
}
impl WcsphSolver {
    /// Create a new WCSPH solver with the given smoothing length.
    pub fn new(smoothing_length: f64) -> Self {
        Self { smoothing_length }
    }
    /// Compute a CFL-limited adaptive timestep.
    ///
    /// Returns `min(dt_max, dt_cfl, dt_acc)` where:
    /// - `dt_cfl = cfl_factor * h / v_max`  (velocity constraint)
    /// - `dt_acc = sqrt(h / a_max)`          (acceleration constraint)
    ///
    /// If all velocities and forces are zero the method returns `dt_max`.
    pub fn compute_adaptive_dt(
        &self,
        velocities: &[Vec3],
        forces: &[Vec3],
        dt_max: f64,
        cfl_factor: f64,
    ) -> f64 {
        let h = self.smoothing_length;
        let v_max: Real = velocities.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
        let a_max: Real = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        let dt_cfl = if v_max > 1e-14 {
            cfl_factor * h / v_max
        } else {
            dt_max
        };
        let dt_acc = if a_max > 1e-14 {
            (h / a_max).sqrt()
        } else {
            dt_max
        };
        dt_max.min(dt_cfl).min(dt_acc)
    }
}
