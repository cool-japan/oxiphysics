//! # SphConfig - Trait Implementations
//!
//! This module contains trait implementations for `SphConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::SphConfig;

/// The default `surface_tension` of `0.0728` is the Akinci coefficient `γ` of
/// [`sph_compute_forces`](super::functions::sph_compute_forces): the number
/// coincides with the surface tension of water in N/m, but it is only a tuning
/// coefficient of the cohesion and curvature model, not a physical calibration.
impl Default for SphConfig {
    fn default() -> Self {
        Self {
            h: 0.1,
            rest_density: 1000.0,
            pressure_k: 100.0,
            viscosity: 0.1,
            surface_tension: 0.0728,
            gravity: [0.0, -9.81, 0.0],
            dt: 0.001,
        }
    }
}
