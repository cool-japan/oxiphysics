// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU/CPU compute kernels for physics simulation.
//!
//! This module groups all low-level compute kernels.  Each sub-module exposes
//! a CPU-mock implementation that mirrors a GPU kernel in its data layout and
//! dispatch model, but executes on the CPU using Rayon for parallelism.

#![allow(dead_code)]

pub mod broadphase;
pub mod md_force;
pub mod rigid;
pub mod sph;

// ── Kernel registry helpers ──────────────────────────────────────────────────

/// Identifier for a built-in kernel family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KernelFamily {
    /// Smoothed Particle Hydrodynamics kernels.
    Sph,
    /// Rigid-body integration and collision kernels.
    Rigid,
    /// Broad-phase AABB/BVH traversal kernels.
    Broadphase,
    /// Molecular dynamics force kernels.
    MdForce,
    /// Signed distance field evaluation kernels.
    SdfCompute,
    /// Neural-network inference kernels.
    NeuralCompute,
    /// Grid-reduce / scan kernels.
    GridReduce,
}

impl std::fmt::Display for KernelFamily {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            KernelFamily::Sph => "sph",
            KernelFamily::Rigid => "rigid",
            KernelFamily::Broadphase => "broadphase",
            KernelFamily::MdForce => "md_force",
            KernelFamily::SdfCompute => "sdf_compute",
            KernelFamily::NeuralCompute => "neural_compute",
            KernelFamily::GridReduce => "grid_reduce",
        };
        write!(f, "{name}")
    }
}

// ── Dispatch descriptor ──────────────────────────────────────────────────────

/// Describes the 3-D work-group dispatch dimensions for a kernel launch.
///
/// Mirrors the `(group_count_x, group_count_y, group_count_z)` triple passed
/// to `vkCmdDispatch` / `wgpuComputePassEncoderDispatchWorkgroups`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DispatchDims {
    /// Number of work-groups in X.
    pub x: u32,
    /// Number of work-groups in Y.
    pub y: u32,
    /// Number of work-groups in Z.
    pub z: u32,
}

impl DispatchDims {
    /// Create a 1-D dispatch of `n` work-groups.
    pub fn linear(n: u32) -> Self {
        Self { x: n, y: 1, z: 1 }
    }

    /// Create a 2-D dispatch.
    pub fn grid2d(x: u32, y: u32) -> Self {
        Self { x, y, z: 1 }
    }

    /// Create a 3-D dispatch.
    pub fn grid3d(x: u32, y: u32, z: u32) -> Self {
        Self { x, y, z }
    }

    /// Total number of work-groups.
    pub fn total_groups(&self) -> u64 {
        self.x as u64 * self.y as u64 * self.z as u64
    }

    /// Total threads given `threads_per_group`.
    pub fn total_threads(&self, threads_per_group: u32) -> u64 {
        self.total_groups() * threads_per_group as u64
    }
}

/// Compute the 1-D dispatch size needed to cover `n` items with `group_size`
/// threads per work-group.
pub fn dispatch_size_1d(n: u32, group_size: u32) -> u32 {
    if group_size == 0 {
        return 0;
    }
    n.div_ceil(group_size)
}

// ── Kernel performance counters ──────────────────────────────────────────────

/// Lightweight performance counters attached to a single kernel invocation.
#[derive(Debug, Clone, Default)]
pub struct KernelPerfCounters {
    /// Number of times the kernel was dispatched.
    pub dispatch_count: u64,
    /// Total elements processed across all dispatches.
    pub elements_processed: u64,
    /// Estimated floating-point operations (MACs counted as 2 FLOPs).
    pub flop_count: u64,
    /// Total bytes read from global memory (mock).
    pub bytes_read: u64,
    /// Total bytes written to global memory (mock).
    pub bytes_written: u64,
}

impl KernelPerfCounters {
    /// Record one dispatch that processed `n` elements.
    pub fn record_dispatch(&mut self, elements: u64, flops: u64, bytes_r: u64, bytes_w: u64) {
        self.dispatch_count += 1;
        self.elements_processed += elements;
        self.flop_count += flops;
        self.bytes_read += bytes_r;
        self.bytes_written += bytes_w;
    }

    /// Arithmetic intensity (FLOPs per byte).
    pub fn arithmetic_intensity(&self) -> f64 {
        let bytes = self.bytes_read + self.bytes_written;
        if bytes == 0 {
            return 0.0;
        }
        self.flop_count as f64 / bytes as f64
    }

    /// Reset all counters.
    pub fn reset(&mut self) {
        *self = KernelPerfCounters::default();
    }
}

// ── Shared-memory size helper ────────────────────────────────────────────────

/// Calculate the shared-memory footprint (bytes) for a tiled matrix-multiply
/// kernel with tiles of size `tile` × `tile` of `T`-sized elements.
pub fn smem_bytes_matmul<T>(tile: usize) -> usize {
    2 * tile * tile * std::mem::size_of::<T>()
}

// ── Barrier simulation ───────────────────────────────────────────────────────

/// Simulated GPU barrier: in CPU mock this is a no-op but documents
/// synchronisation points for future GPU backend porting.
#[inline(always)]
pub fn workgroup_barrier() {
    // CPU: no-op — Rayon fork-join already provides synchronisation.
    std::sync::atomic::fence(std::sync::atomic::Ordering::SeqCst);
}

// ── Predefined group sizes ───────────────────────────────────────────────────

/// Typical work-group sizes used by NVIDIA/AMD GPUs.
pub mod group_sizes {
    /// 64 threads — common on AMD RDNA and for register-heavy kernels.
    pub const WG_64: u32 = 64;
    /// 128 threads — common general-purpose choice.
    pub const WG_128: u32 = 128;
    /// 256 threads — default for many CUDA/Vulkan kernels.
    pub const WG_256: u32 = 256;
    /// 512 threads — useful for reduction passes.
    pub const WG_512: u32 = 512;
    /// 1024 threads — maximum work-group size on most hardware.
    pub const WG_1024: u32 = 1024;
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod kernel_mod_tests {
    use super::*;

    #[test]
    fn test_kernel_family_display() {
        assert_eq!(KernelFamily::Sph.to_string(), "sph");
        assert_eq!(KernelFamily::NeuralCompute.to_string(), "neural_compute");
        assert_eq!(KernelFamily::GridReduce.to_string(), "grid_reduce");
    }

    #[test]
    fn test_dispatch_dims_linear() {
        let d = DispatchDims::linear(128);
        assert_eq!(d.total_groups(), 128);
        assert_eq!(d.total_threads(256), 128 * 256);
    }

    #[test]
    fn test_dispatch_dims_grid3d() {
        let d = DispatchDims::grid3d(4, 4, 4);
        assert_eq!(d.total_groups(), 64);
    }

    #[test]
    fn test_dispatch_size_1d_exact() {
        assert_eq!(dispatch_size_1d(256, 64), 4);
    }

    #[test]
    fn test_dispatch_size_1d_remainder() {
        assert_eq!(dispatch_size_1d(257, 64), 5);
    }

    #[test]
    fn test_dispatch_size_1d_zero_group() {
        assert_eq!(dispatch_size_1d(100, 0), 0);
    }

    #[test]
    fn test_perf_counters_arithmetic_intensity() {
        let mut c = KernelPerfCounters::default();
        c.record_dispatch(1024, 8192, 4096, 4096);
        // intensity = 8192 / 8192 = 1.0
        assert!((c.arithmetic_intensity() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_perf_counters_reset() {
        let mut c = KernelPerfCounters::default();
        c.record_dispatch(512, 1024, 512, 512);
        c.reset();
        assert_eq!(c.dispatch_count, 0);
        assert_eq!(c.flop_count, 0);
    }

    #[test]
    fn test_smem_bytes_matmul_f32() {
        // 2 * 16 * 16 * 4 bytes = 2048
        let bytes = smem_bytes_matmul::<f32>(16);
        assert_eq!(bytes, 2048);
    }

    #[test]
    fn test_smem_bytes_matmul_f64() {
        // 2 * 16 * 16 * 8 bytes = 4096
        let bytes = smem_bytes_matmul::<f64>(16);
        assert_eq!(bytes, 4096);
    }

    #[test]
    fn test_workgroup_barrier_no_panic() {
        workgroup_barrier(); // must not panic
    }

    #[allow(clippy::assertions_on_constants)]
    #[test]
    fn test_group_sizes_constants() {
        use group_sizes::*;
        assert!(WG_64 < WG_128);
        assert!(WG_128 < WG_256);
        assert!(WG_256 < WG_512);
        assert!(WG_512 < WG_1024);
    }
}
