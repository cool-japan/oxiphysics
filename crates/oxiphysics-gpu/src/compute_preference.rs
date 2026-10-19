// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compute-device preference shared by the simulation configurations.
//!
//! [`SphConfig`](crate::sph_gpu::SphConfig) and
//! [`LbmConfig`](crate::lbm_gpu::LbmConfig) carry a [`ComputePreference`] that
//! tells the simulation whether it may, must not, or must use the GPU path.
//! Benchmarks and parity tests use it to be sure which device actually ran.

/// Which compute device a simulation is allowed or required to use.
///
/// ```
/// use oxiphysics_gpu::compute_preference::ComputePreference;
/// use oxiphysics_gpu::sph_gpu::{SphConfig, SphSimulation};
///
/// // A CPU reference run: no GPU adapter is requested, no device is created.
/// let cfg = SphConfig {
///     n_particles: 8,
///     compute: ComputePreference::CpuOnly,
///     ..SphConfig::default()
/// };
/// let sim = SphSimulation::new(cfg);
/// assert!(!sim.has_gpu());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ComputePreference {
    /// Use the GPU when an adapter is available and fall back to the CPU
    /// reference otherwise (the default).
    #[default]
    Auto,
    /// Never touch the GPU: no adapter is requested and no device resources
    /// are created; every step runs on the CPU reference.
    CpuOnly,
    /// The GPU path must be used.  A missing adapter or a failed dispatch is
    /// reported as a [`GpuError`](crate::GpuError) by the fallible entry points
    /// (`try_new`, `try_step`) and is never replaced by a CPU step; the
    /// infallible `step` leaves the simulation unchanged in that case.
    GpuRequired,
}

impl ComputePreference {
    /// `true` when a CPU step may stand in for a GPU step that cannot run
    /// (`Auto` and `CpuOnly`).
    pub fn allows_cpu(self) -> bool {
        !matches!(self, Self::GpuRequired)
    }

    /// `true` when the GPU path may be used at all (`Auto` and `GpuRequired`).
    pub fn allows_gpu(self) -> bool {
        !matches!(self, Self::CpuOnly)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_preference_is_auto() {
        assert_eq!(ComputePreference::default(), ComputePreference::Auto);
    }

    #[test]
    fn preference_permissions() {
        assert!(ComputePreference::Auto.allows_cpu());
        assert!(ComputePreference::Auto.allows_gpu());
        assert!(ComputePreference::CpuOnly.allows_cpu());
        assert!(!ComputePreference::CpuOnly.allows_gpu());
        assert!(!ComputePreference::GpuRequired.allows_cpu());
        assert!(ComputePreference::GpuRequired.allows_gpu());
    }
}
