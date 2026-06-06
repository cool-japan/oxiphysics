// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU acceleration backends for the OxiPhysics engine.
//!
//! This crate provides a GPU compute abstraction layer that can work with any
//! backend, with a CPU fallback as the default implementation. No heavy GPU
//! dependencies (such as wgpu) are required.

mod error;
pub use error::*;

pub mod bvh;
pub mod cell_list;
pub mod compute;
pub mod compute_pipeline;
pub mod flux_compute;
pub mod gpu_bench;
pub mod grid_reduce;
pub mod kernels;
pub mod lbm_gpu;
pub mod neural_compute;
pub mod parallel;
pub mod parallel_sort;
pub mod particle_system;
pub mod pipeline;
pub mod sdf_compute;
pub mod shader_registry;
pub mod shaders;
pub mod sparse_gpu;
pub mod sph_gpu;

pub use compute::{BufferHandle, ComputeBackend, ComputeKernel, CpuBackend};
pub use neural_compute::*;
pub use particle_system::functions::*;
pub use particle_system::types::*;
pub use sparse_gpu::*;

// ── GPU compute utility functions ───────────────────────────────────────────

/// Compute the optimal work group size for a given total work item count.
///
/// Rounds up `total` to the next multiple of `group_size`.
pub fn dispatch_count(total: usize, group_size: usize) -> usize {
    if group_size == 0 {
        return 0;
    }
    total.div_ceil(group_size)
}

/// Compute the padded buffer size to meet alignment requirements.
///
/// Returns the smallest multiple of `alignment` that is >= `size`.
pub fn aligned_size(size: usize, alignment: usize) -> usize {
    if alignment == 0 {
        return size;
    }
    size.div_ceil(alignment) * alignment
}

/// Flatten a 3D dispatch (x, y, z) into a linear index, given grid dimensions.
pub fn linear_index_3d(x: usize, y: usize, z: usize, dim_x: usize, dim_y: usize) -> usize {
    z * dim_x * dim_y + y * dim_x + x
}

/// Convert a linear index back to 3D coordinates.
pub fn index_3d_from_linear(index: usize, dim_x: usize, dim_y: usize) -> (usize, usize, usize) {
    let z = index / (dim_x * dim_y);
    let rem = index % (dim_x * dim_y);
    let y = rem / dim_x;
    let x = rem % dim_x;
    (x, y, z)
}

/// A simple timer utility for profiling GPU-like dispatches.
#[derive(Debug, Clone)]
pub struct DispatchTimer {
    /// Label for this dispatch.
    pub label: String,
    /// Elapsed time in seconds (set after timing).
    pub elapsed_secs: f64,
}

impl DispatchTimer {
    /// Create a new timer with the given label.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            elapsed_secs: 0.0,
        }
    }

    /// Record elapsed time.
    pub fn record(&mut self, elapsed: f64) {
        self.elapsed_secs = elapsed;
    }
}

/// Estimate memory bandwidth in GB/s.
///
/// * `bytes_transferred` - Total bytes read + written.
/// * `elapsed_secs` - Elapsed time in seconds.
pub fn bandwidth_gb_s(bytes_transferred: usize, elapsed_secs: f64) -> f64 {
    if elapsed_secs <= 0.0 {
        return 0.0;
    }
    (bytes_transferred as f64) / elapsed_secs / 1e9
}

/// Compute the number of elements that fit in a given memory budget.
///
/// * `budget_bytes` - Available memory in bytes.
/// * `element_size` - Size of one element in bytes.
pub fn elements_in_budget(budget_bytes: usize, element_size: usize) -> usize {
    if element_size == 0 {
        return 0;
    }
    budget_bytes / element_size
}

// ── GPU buffer utilities ─────────────────────────────────────────────────────

/// Stride (in bytes) of a row in a 2-D buffer, given the element count per row
/// and the required alignment.
///
/// This mirrors `wgpuDeviceGetSupportedSurfaceFormats` style pitch calculation.
pub fn row_pitch(elements_per_row: usize, element_size: usize, alignment: usize) -> usize {
    let raw = elements_per_row * element_size;
    aligned_size(raw, alignment)
}

/// Compute the 2-D buffer size (rows × pitch) for a texture-like allocation.
pub fn buffer_size_2d(
    width: usize,
    height: usize,
    element_size: usize,
    row_alignment: usize,
) -> usize {
    row_pitch(width, element_size, row_alignment) * height
}

/// Round `value` up to the next power of two.
///
/// Returns `value` unchanged when it is already a power of two.
/// Returns 1 when `value` is 0.
pub fn next_power_of_two(value: usize) -> usize {
    if value == 0 {
        return 1;
    }
    let mut p = 1usize;
    while p < value {
        p <<= 1;
    }
    p
}

/// True when `value` is a power of two (including 1).
pub fn is_power_of_two(value: usize) -> bool {
    value != 0 && (value & (value - 1)) == 0
}

/// Log2 of a power-of-two value.  Panics in debug mode if `v` is not a power
/// of two.
pub fn log2_pow2(v: usize) -> u32 {
    debug_assert!(is_power_of_two(v), "{v} is not a power of two");
    v.trailing_zeros()
}

// ── Work-group tiling helpers ─────────────────────────────────────────────────

/// Divides a 2-D problem of `(width × height)` into tiles of `(tw × th)` and
/// returns `(tiles_x, tiles_y)`.
///
/// Each dimension is rounded up so the full problem is covered.
pub fn tile_count_2d(width: usize, height: usize, tw: usize, th: usize) -> (usize, usize) {
    let tx = width.div_ceil(tw);
    let ty = height.div_ceil(th);
    (tx, ty)
}

/// Total number of tiles for a 2-D problem.
pub fn total_tiles_2d(width: usize, height: usize, tw: usize, th: usize) -> usize {
    let (tx, ty) = tile_count_2d(width, height, tw, th);
    tx * ty
}

/// Convert a flat tile index back to `(tile_x, tile_y)` for a grid with
/// `tiles_x` columns.
pub fn tile_index_to_2d(flat: usize, tiles_x: usize) -> (usize, usize) {
    (flat % tiles_x, flat / tiles_x)
}

// ── Numeric helpers used across GPU kernels ──────────────────────────────────

/// Clamp `v` to `[lo, hi]`.
pub fn clamp_f64(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// Smooth-step function: `3t² - 2t³` with `t = (v - lo) / (hi - lo)`.
pub fn smoothstep(lo: f64, hi: f64, v: f64) -> f64 {
    let t = clamp_f64((v - lo) / (hi - lo), 0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Smoother-step (Ken Perlin's quintic): `6t⁵ − 15t⁴ + 10t³`.
pub fn smootherstep(lo: f64, hi: f64, v: f64) -> f64 {
    let t = clamp_f64((v - lo) / (hi - lo), 0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

/// Linear interpolation: `a + t*(b-a)`.
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Inverse lerp: returns `t` such that `lerp(a, b, t) == v`, or `0` if `a==b`.
pub fn inv_lerp(a: f64, b: f64, v: f64) -> f64 {
    if (b - a).abs() < f64::EPSILON {
        return 0.0;
    }
    (v - a) / (b - a)
}

// ── FP utilities ─────────────────────────────────────────────────────────────

/// Safe reciprocal: returns `1/x` when `|x| > eps`, else `0`.
pub fn safe_recip(x: f64, eps: f64) -> f64 {
    if x.abs() > eps { 1.0 / x } else { 0.0 }
}

/// Safe square root: clamps negative values to 0 before taking sqrt.
pub fn safe_sqrt(x: f64) -> f64 {
    x.max(0.0).sqrt()
}

/// Wrap an angle in radians to `(-π, π]`.
pub fn wrap_angle(theta: f64) -> f64 {
    use std::f64::consts::PI;
    let mut t = theta % (2.0 * PI);
    if t > PI {
        t -= 2.0 * PI;
    }
    if t <= -PI {
        t += 2.0 * PI;
    }
    t
}

// ── Vector math (3-D, f64) ────────────────────────────────────────────────────

/// Compute the dot product of two 3-element arrays.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the cross product of two 3-element arrays.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Length of a 3-D vector.
pub fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise a 3-D vector.  Returns the zero vector if the length is < eps.
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = length3(v);
    if len < 1e-15 {
        return [0.0; 3];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}

/// Reflect vector `d` about normal `n` (both assumed normalised).
pub fn reflect3(d: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    let dn2 = 2.0 * dot3(d, n);
    [d[0] - dn2 * n[0], d[1] - dn2 * n[1], d[2] - dn2 * n[2]]
}

// ── Parallel prefix sum (scan) ───────────────────────────────────────────────

/// Parallel prefix sum (scan) on a slice of f64 values.
///
/// Returns a new vector where `result[i] = sum(data[0..i])`.
/// This is the exclusive scan variant.
pub fn exclusive_scan(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0;
    for &v in data {
        result.push(acc);
        acc += v;
    }
    result
}

/// Inclusive scan: `result[i] = sum(data[0..=i])`.
pub fn inclusive_scan(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0;
    for &v in data {
        acc += v;
        result.push(acc);
    }
    result
}

/// Parallel reduce: compute the sum of all elements.
pub fn reduce_sum(data: &[f64]) -> f64 {
    data.iter().copied().sum()
}

/// Parallel reduce: compute the maximum of all elements.
pub fn reduce_max(data: &[f64]) -> f64 {
    data.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Parallel reduce: compute the minimum of all elements.
pub fn reduce_min(data: &[f64]) -> f64 {
    data.iter().copied().fold(f64::INFINITY, f64::min)
}

#[cfg(test)]
mod gpu_util_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_dispatch_count_exact() {
        assert_eq!(dispatch_count(256, 64), 4);
    }

    #[test]
    fn test_dispatch_count_remainder() {
        assert_eq!(dispatch_count(257, 64), 5);
    }

    #[test]
    fn test_dispatch_count_zero_group() {
        assert_eq!(dispatch_count(100, 0), 0);
    }

    #[test]
    fn test_aligned_size_exact() {
        assert_eq!(aligned_size(256, 64), 256);
    }

    #[test]
    fn test_aligned_size_pad() {
        assert_eq!(aligned_size(257, 64), 320);
    }

    #[test]
    fn test_aligned_size_zero_alignment() {
        assert_eq!(aligned_size(100, 0), 100);
    }

    #[test]
    fn test_linear_index_3d() {
        // Grid 4x3x2
        assert_eq!(linear_index_3d(0, 0, 0, 4, 3), 0);
        assert_eq!(linear_index_3d(3, 2, 1, 4, 3), 12 + 2 * 4 + 3);
    }

    #[test]
    fn test_index_3d_roundtrip() {
        let (dx, dy) = (4, 3);
        for z in 0..2 {
            for y in 0..dy {
                for x in 0..dx {
                    let idx = linear_index_3d(x, y, z, dx, dy);
                    let (rx, ry, rz) = index_3d_from_linear(idx, dx, dy);
                    assert_eq!((rx, ry, rz), (x, y, z));
                }
            }
        }
    }

    #[test]
    fn test_dispatch_timer() {
        let mut timer = DispatchTimer::new("test");
        assert_eq!(timer.label, "test");
        timer.record(0.5);
        assert!((timer.elapsed_secs - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_bandwidth_gb_s() {
        // 1 GB in 1 second = 1 GB/s
        let bw = bandwidth_gb_s(1_000_000_000, 1.0);
        assert!((bw - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bandwidth_zero_time() {
        assert!((bandwidth_gb_s(1000, 0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_elements_in_budget() {
        assert_eq!(elements_in_budget(1024, 4), 256);
        assert_eq!(elements_in_budget(1024, 0), 0);
    }

    #[test]
    fn test_exclusive_scan() {
        let data = [1.0, 2.0, 3.0, 4.0];
        let result = exclusive_scan(&data);
        assert_eq!(result, vec![0.0, 1.0, 3.0, 6.0]);
    }

    #[test]
    fn test_inclusive_scan() {
        let data = [1.0, 2.0, 3.0, 4.0];
        let result = inclusive_scan(&data);
        assert_eq!(result, vec![1.0, 3.0, 6.0, 10.0]);
    }

    #[test]
    fn test_reduce_sum() {
        assert!((reduce_sum(&[1.0, 2.0, 3.0]) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_reduce_max() {
        assert!((reduce_max(&[1.0, 5.0, 3.0]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_reduce_min() {
        assert!((reduce_min(&[1.0, 5.0, 3.0]) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_exclusive_scan_empty() {
        let result = exclusive_scan(&[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_inclusive_scan_single() {
        let result = inclusive_scan(&[42.0]);
        assert_eq!(result, vec![42.0]);
    }

    // ── Buffer utility tests ─────────────────────────────────────────────

    #[test]
    fn test_row_pitch_aligned() {
        // 128 elements × 4 bytes = 512 bytes, already aligned to 256
        assert_eq!(row_pitch(128, 4, 256), 512);
    }

    #[test]
    fn test_row_pitch_needs_padding() {
        // 100 × 4 = 400; aligned to 256 => 512
        assert_eq!(row_pitch(100, 4, 256), 512);
    }

    #[test]
    fn test_buffer_size_2d() {
        // 4 rows × pitch(64 elems × 4 bytes, align 256) = 4 × 256 = 1024
        assert_eq!(buffer_size_2d(64, 4, 4, 256), 1024);
    }

    #[test]
    fn test_next_power_of_two() {
        assert_eq!(next_power_of_two(0), 1);
        assert_eq!(next_power_of_two(1), 1);
        assert_eq!(next_power_of_two(5), 8);
        assert_eq!(next_power_of_two(8), 8);
        assert_eq!(next_power_of_two(9), 16);
    }

    #[test]
    fn test_is_power_of_two() {
        assert!(is_power_of_two(1));
        assert!(is_power_of_two(16));
        assert!(!is_power_of_two(0));
        assert!(!is_power_of_two(7));
    }

    #[test]
    fn test_log2_pow2() {
        assert_eq!(log2_pow2(1), 0);
        assert_eq!(log2_pow2(2), 1);
        assert_eq!(log2_pow2(256), 8);
    }

    #[test]
    fn test_tile_count_2d_exact() {
        let (tx, ty) = tile_count_2d(64, 64, 16, 16);
        assert_eq!(tx, 4);
        assert_eq!(ty, 4);
    }

    #[test]
    fn test_tile_count_2d_remainder() {
        let (tx, ty) = tile_count_2d(65, 65, 16, 16);
        assert_eq!(tx, 5);
        assert_eq!(ty, 5);
    }

    #[test]
    fn test_total_tiles_2d() {
        assert_eq!(total_tiles_2d(64, 64, 16, 16), 16);
    }

    #[test]
    fn test_tile_index_to_2d() {
        // tiles_x = 4; flat=5 => (1, 1)
        assert_eq!(tile_index_to_2d(5, 4), (1, 1));
        assert_eq!(tile_index_to_2d(0, 4), (0, 0));
    }

    // ── Numeric helpers tests ────────────────────────────────────────────

    #[test]
    fn test_smoothstep_edges() {
        assert!((smoothstep(0.0, 1.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((smoothstep(0.0, 1.0, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_smoothstep_midpoint() {
        // at t=0.5: 3*(0.25) - 2*(0.125) = 0.75 - 0.25 = 0.5
        assert!((smoothstep(0.0, 1.0, 0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_smootherstep_edges() {
        assert!((smootherstep(0.0, 1.0, 0.0)).abs() < 1e-12);
        assert!((smootherstep(0.0, 1.0, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_inv_lerp_roundtrip() {
        let a = 10.0;
        let b = 20.0;
        let t = 0.3;
        let v = lerp(a, b, t);
        assert!((inv_lerp(a, b, v) - t).abs() < 1e-12);
    }

    #[test]
    fn test_safe_recip_normal() {
        assert!((safe_recip(2.0, 1e-9) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_safe_recip_near_zero() {
        assert!((safe_recip(1e-15, 1e-9)).abs() < 1e-12);
    }

    #[test]
    fn test_safe_sqrt_positive() {
        assert!((safe_sqrt(9.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_safe_sqrt_negative() {
        assert!((safe_sqrt(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_wrap_angle_in_range() {
        let wrapped = wrap_angle(3.0 * PI);
        assert!(wrapped.abs() <= PI + 1e-12, "wrapped = {wrapped}");
    }

    // ── Vector math tests ────────────────────────────────────────────────

    #[test]
    fn test_dot3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot3(a, b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_cross3() {
        let i = [1.0, 0.0, 0.0];
        let j = [0.0, 1.0, 0.0];
        let k = cross3(i, j);
        assert!((k[0]).abs() < 1e-12);
        assert!((k[1]).abs() < 1e-12);
        assert!((k[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_length3() {
        let v = [3.0, 4.0, 0.0];
        assert!((length3(v) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3() {
        let v = [0.0, 0.0, 5.0];
        let n = normalize3(v);
        assert!((length3(n) - 1.0).abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize3_zero_vec() {
        let n = normalize3([0.0; 3]);
        assert_eq!(n, [0.0; 3]);
    }

    #[test]
    fn test_reflect3() {
        // Reflect (1,0,0) about (0,1,0) => still (1,0,0) but inverted y
        let d = [0.0, -1.0, 0.0]; // pointing down
        let n = [0.0, 1.0, 0.0]; // surface normal up
        let r = reflect3(d, n);
        // r = d - 2*(d·n)*n = [0,-1,0] - 2*(-1)*[0,1,0] = [0,1,0]
        assert!((r[1] - 1.0).abs() < 1e-12);
    }
}
pub mod collision_gpu;
pub mod deformable_gpu;
pub mod fluid_gpu;
pub mod fluid_sim_gpu;
pub mod gpu_cloth;
pub mod gpu_collision_detection;
pub mod gpu_collision_ext;
pub mod gpu_fem_assembly;
pub mod gpu_fluid;
pub mod gpu_fluid_euler;
pub mod gpu_lbm;
pub mod gpu_md_solver;
pub mod gpu_mesh_processing;
pub mod gpu_neural_solver;
pub mod gpu_nn;
pub mod gpu_particle_system;
pub mod gpu_particles;
pub mod gpu_ray_tracing;
pub mod gpu_reduction;
pub mod gpu_rigid;
pub mod gpu_sdf;
pub mod gpu_sort;
pub mod gpu_sparse_solver;
pub mod gpu_sph_density;
pub mod gpu_sph_pressure;
pub mod gpu_sph_solver;
pub mod gpu_thermal;
pub mod gpu_voxel;
pub mod memory;
pub mod neural_physics;
pub mod path_tracer;
pub mod ray_marching;
pub mod ray_tracing_gpu;
pub mod raytracing;
pub mod scheduler;
