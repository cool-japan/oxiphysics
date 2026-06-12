// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WGSL source for the GPU reduction-primitive kernels.
//!
//! Each constant embeds a `.wgsl` shader via [`include_str!`]; the shaders are
//! dispatched by [`crate::gpu_primitives`] against the real wgpu backend.

/// Exclusive prefix-scan kernel (`scan_block` entry point) over `u32`.
pub const SCAN_WGSL: &str = include_str!("scan.wgsl");

/// Uniform block-offset add kernel (`add_block_offsets` entry point).
pub const SCAN_ADD_WGSL: &str = include_str!("scan_add.wgsl");

/// Tree-reduction kernels (`reduce_sum_block`, `reduce_max_block`) over `f32`.
pub const REDUCE_WGSL: &str = include_str!("reduce.wgsl");

/// Stream-compaction scatter kernel (`scatter_kept` entry point).
pub const COMPACT_WGSL: &str = include_str!("compact.wgsl");

/// Histogram kernel (`histogram_main` entry point) over `u32`.
pub const HISTOGRAM_WGSL: &str = include_str!("histogram.wgsl");

/// Radix-sort digit histogram kernel (`radix_histogram` entry point) over `u32`.
pub const RADIX_HISTOGRAM_WGSL: &str = include_str!("radix_histogram.wgsl");

/// Radix-sort stable scatter kernel, keys only (`radix_scatter` entry point).
pub const RADIX_SCATTER_WGSL: &str = include_str!("radix_scatter.wgsl");

/// Radix-sort stable scatter kernel for key+payload pairs (`radix_scatter_pairs` entry point).
pub const RADIX_SCATTER_PAIRS_WGSL: &str = include_str!("radix_scatter_pairs.wgsl");

/// SPH spatial-hash cell-key generation (`sph_cell_list` entry point).
pub const SPH_CELL_LIST_WGSL: &str = include_str!("sph_cell_list.wgsl");

/// SPH density summation over the cell-list (`sph_density` entry point).
pub const SPH_DENSITY_GRID_WGSL: &str = include_str!("sph_density.wgsl");

/// SPH pressure + viscosity + gravity acceleration (`sph_force` entry point).
pub const SPH_FORCE_GRID_WGSL: &str = include_str!("sph_force.wgsl");

/// SPH symplectic-Euler integration (`sph_integrate` entry point).
pub const SPH_INTEGRATE_WGSL: &str = include_str!("sph_integrate.wgsl");

/// SPH axis-aligned box boundary enforcement (`sph_boundary` entry point).
pub const SPH_BOUNDARY_WGSL: &str = include_str!("sph_boundary.wgsl");
