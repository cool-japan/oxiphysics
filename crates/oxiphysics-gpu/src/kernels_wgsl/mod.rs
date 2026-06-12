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
