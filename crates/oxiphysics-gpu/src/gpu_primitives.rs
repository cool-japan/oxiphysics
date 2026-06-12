// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Safe GPU reduction primitives with transparent CPU fallback.
//!
//! Public entry points ([`exclusive_scan_u32`], [`reduce_sum_f32`],
//! [`reduce_max_f32`], [`compact_u32`], [`histogram_u32`]) attempt the real
//! wgpu backend (when the `wgpu-backend` feature is enabled and a GPU adapter
//! is present) and fall back to the CPU reference implementation otherwise.
//! The public signatures never change between feature configurations.
//!
//! The CPU implementations in [`crate::gpu_reduction`] and
//! [`crate::grid_reduce`] act as the parity oracles for the GPU kernels.

/// Exclusive prefix sum over `u32` (GPU when available, else CPU).
///
/// Overflow wraps (matching the GPU's modular `u32` arithmetic).
pub fn exclusive_scan_u32(data: &[u32]) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) =
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::exclusive_scan_u32(&mut backend, data)
        {
            return result;
        }
    }
    cpu::exclusive_scan_u32(data)
}

/// Sum-reduce an `f32` slice (GPU when available, else CPU). Empty → `0.0`.
pub fn reduce_sum_f32(data: &[f32]) -> f32 {
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) =
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::reduce_sum_f32(&mut backend, data)
        {
            return result;
        }
    }
    cpu::reduce_sum_f32(data)
}

/// Max-reduce an `f32` slice (GPU when available, else CPU). Empty → `NEG_INFINITY`.
pub fn reduce_max_f32(data: &[f32]) -> f32 {
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) =
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::reduce_max_f32(&mut backend, data)
        {
            return result;
        }
    }
    cpu::reduce_max_f32(data)
}

/// Stable stream compaction: keep `values[i]` where `keep[i] != 0`.
///
/// `keep.len()` must equal `values.len()`; otherwise an empty vec is returned.
pub fn compact_u32(values: &[u32], keep: &[u32]) -> Vec<u32> {
    if values.len() != keep.len() {
        return Vec::new();
    }
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) =
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::compact_u32(&mut backend, values, keep)
        {
            return result;
        }
    }
    cpu::compact_u32(values, keep)
}

/// Histogram over `u32`: value `v` increments bin `min(v, num_bins - 1)`.
pub fn histogram_u32(data: &[u32], num_bins: usize) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        if let Ok(mut backend) =
            crate::compute::wgpu_backend::real::WgpuBackendReal::try_new()
            && let Ok(result) = gpu::histogram_u32(&mut backend, data, num_bins)
        {
            return result;
        }
    }
    cpu::histogram_u32(data, num_bins)
}

// ── CPU reference (parity oracle) ─────────────────────────────────────────────

mod cpu {
    pub(super) fn exclusive_scan_u32(data: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(data.len());
        let mut acc: u32 = 0;
        for &v in data {
            out.push(acc);
            acc = acc.wrapping_add(v);
        }
        out
    }

    pub(super) fn reduce_sum_f32(data: &[f32]) -> f32 {
        let mut s = 0f32;
        for &v in data {
            s += v;
        }
        s
    }

    pub(super) fn reduce_max_f32(data: &[f32]) -> f32 {
        data.iter().copied().fold(f32::NEG_INFINITY, f32::max)
    }

    pub(super) fn compact_u32(values: &[u32], keep: &[u32]) -> Vec<u32> {
        values
            .iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k != 0)
            .map(|(&v, _)| v)
            .collect()
    }

    pub(super) fn histogram_u32(data: &[u32], num_bins: usize) -> Vec<u32> {
        if num_bins == 0 {
            return Vec::new();
        }
        let mut bins = vec![0u32; num_bins];
        let last = num_bins - 1;
        for &v in data {
            let idx = (v as usize).min(last);
            bins[idx] = bins[idx].wrapping_add(1);
        }
        bins
    }
}

// ── GPU implementation (real wgpu) ────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
mod gpu {
    use crate::compute::wgpu_backend::real::WgpuBackendReal;
    use crate::compute::wgpu_backend::WgpuInitError;
    use crate::kernels_wgsl::{
        COMPACT_WGSL, HISTOGRAM_WGSL, REDUCE_WGSL, SCAN_ADD_WGSL, SCAN_WGSL,
    };

    const BLOCK: usize = 256;

    fn params_bytes(n: u32, num_bins: u32) -> Vec<u8> {
        let p: [u32; 4] = [n, num_bins, 0, 0];
        bytemuck::cast_slice(&p).to_vec()
    }

    fn ro() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: false }
    }

    pub(super) fn exclusive_scan_u32(
        backend: &mut WgpuBackendReal,
        data: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = data.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if n == 1 {
            return Ok(vec![0u32]);
        }

        let num_blocks = n.div_ceil(BLOCK);
        let padded = num_blocks * BLOCK;

        let mut padded_in = vec![0u32; padded];
        padded_in[..n].copy_from_slice(data);

        let input_buf = backend.create_buffer_storage((padded * 4) as u64);
        let output_buf = backend.create_buffer_storage((padded * 4) as u64);
        let block_sums_buf = backend.create_buffer_storage((num_blocks.max(1) * 4) as u64);
        let params_buf = backend.create_buffer_storage(16);

        backend.queue_write_buffer_raw(&input_buf, bytemuck::cast_slice(&padded_in));
        backend.queue_write_buffer_raw(&params_buf, &params_bytes(n as u32, 0));

        backend.dispatch_wgsl(
            SCAN_WGSL,
            "scan_block",
            &[
                (input_buf, ro()),
                (output_buf, rw()),
                (block_sums_buf, rw()),
                (params_buf, ro()),
            ],
            [num_blocks as u32, 1, 1],
        )?;

        if num_blocks > 1 {
            let block_sums = backend.read_buffer_u32(block_sums_buf);
            let block_sums = &block_sums[..num_blocks.min(block_sums.len())];
            let block_offsets = exclusive_scan_u32(backend, block_sums)?;

            let offsets_buf = backend.create_buffer_storage((num_blocks * 4) as u64);
            backend.queue_write_buffer_raw(&offsets_buf, bytemuck::cast_slice(&block_offsets));

            backend.dispatch_wgsl(
                SCAN_ADD_WGSL,
                "add_block_offsets",
                &[(output_buf, rw()), (offsets_buf, ro()), (params_buf, ro())],
                [num_blocks as u32, 1, 1],
            )?;
        }

        let mut out = backend.read_buffer_u32(output_buf);
        out.truncate(n);
        Ok(out)
    }

    pub(super) fn reduce_sum_f32(
        backend: &mut WgpuBackendReal,
        data: &[f32],
    ) -> Result<f32, WgpuInitError> {
        if data.is_empty() {
            return Ok(0.0);
        }
        let mut current: Vec<f32> = data.to_vec();
        loop {
            let n = current.len();
            let num_blocks = n.div_ceil(BLOCK);
            let padded = num_blocks * BLOCK;

            let mut padded_in = vec![0f32; padded];
            padded_in[..n].copy_from_slice(&current);

            let in_buf = backend.create_buffer_storage((padded * 4) as u64);
            let partials_buf = backend.create_buffer_storage((num_blocks * 4) as u64);
            let params_buf = backend.create_buffer_storage(16);

            backend.queue_write_buffer_raw(&in_buf, bytemuck::cast_slice(&padded_in));
            backend.queue_write_buffer_raw(&params_buf, &params_bytes(n as u32, 0));

            backend.dispatch_wgsl(
                REDUCE_WGSL,
                "reduce_sum_block",
                &[(in_buf, ro()), (partials_buf, rw()), (params_buf, ro())],
                [num_blocks as u32, 1, 1],
            )?;

            let partials = backend.read_buffer_f32(partials_buf);
            if num_blocks <= 1 {
                return Ok(partials.first().copied().unwrap_or(0.0));
            }
            current = partials;
            current.truncate(num_blocks);
        }
    }

    pub(super) fn reduce_max_f32(
        backend: &mut WgpuBackendReal,
        data: &[f32],
    ) -> Result<f32, WgpuInitError> {
        if data.is_empty() {
            return Ok(f32::NEG_INFINITY);
        }
        const NEG_FLT_MAX: f32 = -3.4028235e38;
        let mut current: Vec<f32> = data.to_vec();
        loop {
            let n = current.len();
            let num_blocks = n.div_ceil(BLOCK);
            let padded = num_blocks * BLOCK;

            let mut padded_in = vec![NEG_FLT_MAX; padded];
            padded_in[..n].copy_from_slice(&current);

            let in_buf = backend.create_buffer_storage((padded * 4) as u64);
            let partials_buf = backend.create_buffer_storage((num_blocks * 4) as u64);
            let params_buf = backend.create_buffer_storage(16);

            backend.queue_write_buffer_raw(&in_buf, bytemuck::cast_slice(&padded_in));
            backend.queue_write_buffer_raw(&params_buf, &params_bytes(n as u32, 0));

            backend.dispatch_wgsl(
                REDUCE_WGSL,
                "reduce_max_block",
                &[(in_buf, ro()), (partials_buf, rw()), (params_buf, ro())],
                [num_blocks as u32, 1, 1],
            )?;

            let partials = backend.read_buffer_f32(partials_buf);
            if num_blocks <= 1 {
                return Ok(partials.first().copied().unwrap_or(f32::NEG_INFINITY));
            }
            current = partials;
            current.truncate(num_blocks);
        }
    }

    pub(super) fn compact_u32(
        backend: &mut WgpuBackendReal,
        values: &[u32],
        keep: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = values.len();
        if n == 0 {
            return Ok(Vec::new());
        }

        let positions = exclusive_scan_u32(backend, keep)?;
        let last_pos = positions.get(n - 1).copied().unwrap_or(0) as usize;
        let last_keep = keep.get(n - 1).copied().unwrap_or(0) as usize;
        let total_kept = last_pos + if last_keep != 0 { 1 } else { 0 };
        if total_kept == 0 {
            return Ok(Vec::new());
        }

        let values_buf = backend.create_buffer_storage((n * 4) as u64);
        let keep_buf = backend.create_buffer_storage((n * 4) as u64);
        let positions_buf = backend.create_buffer_storage((n * 4) as u64);
        let output_buf = backend.create_buffer_storage((total_kept * 4) as u64);
        let params_buf = backend.create_buffer_storage(16);

        backend.queue_write_buffer_raw(&values_buf, bytemuck::cast_slice(values));
        backend.queue_write_buffer_raw(&keep_buf, bytemuck::cast_slice(keep));
        backend.queue_write_buffer_raw(&positions_buf, bytemuck::cast_slice(&positions));
        backend.queue_write_buffer_raw(&params_buf, &params_bytes(n as u32, 0));

        let workgroups = WgpuBackendReal::dispatch_count_for(n, 256);
        backend.dispatch_wgsl(
            COMPACT_WGSL,
            "scatter_kept",
            &[
                (values_buf, ro()),
                (keep_buf, ro()),
                (positions_buf, ro()),
                (output_buf, rw()),
                (params_buf, ro()),
            ],
            workgroups,
        )?;

        let mut out = backend.read_buffer_u32(output_buf);
        out.truncate(total_kept);
        Ok(out)
    }

    pub(super) fn histogram_u32(
        backend: &mut WgpuBackendReal,
        data: &[u32],
        num_bins: usize,
    ) -> Result<Vec<u32>, WgpuInitError> {
        if num_bins == 0 {
            return Ok(Vec::new());
        }
        if num_bins > 256 {
            return Err(WgpuInitError::FeatureDisabled);
        }
        let n = data.len();
        if n == 0 {
            return Ok(vec![0u32; num_bins]);
        }

        let data_buf = backend.create_buffer_storage((n.max(1) * 4) as u64);
        let global_bins_buf = backend.create_buffer_storage((num_bins * 4) as u64);
        let params_buf = backend.create_buffer_storage(16);

        backend.queue_write_buffer_raw(&data_buf, bytemuck::cast_slice(data));
        let zero_bins = vec![0u32; num_bins];
        backend.queue_write_buffer_raw(&global_bins_buf, bytemuck::cast_slice(&zero_bins));
        backend.queue_write_buffer_raw(&params_buf, &params_bytes(n as u32, num_bins as u32));

        let workgroups = WgpuBackendReal::dispatch_count_for(n, 256);
        backend.dispatch_wgsl(
            HISTOGRAM_WGSL,
            "histogram_main",
            &[(data_buf, ro()), (global_bins_buf, rw()), (params_buf, ro())],
            workgroups,
        )?;

        let mut out = backend.read_buffer_u32(global_bins_buf);
        out.truncate(num_bins);
        Ok(out)
    }
}

// ── Tests (CPU oracle unit tests; GPU parity lives in tests/gpu_primitives.rs) ─

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_exclusive_scan_basic() {
        assert_eq!(cpu::exclusive_scan_u32(&[1, 2, 3, 4]), vec![0, 1, 3, 6]);
        assert_eq!(cpu::exclusive_scan_u32(&[]), Vec::<u32>::new());
        assert_eq!(cpu::exclusive_scan_u32(&[7]), vec![0]);
    }

    #[test]
    fn cpu_reduce_sum_and_max() {
        assert_eq!(cpu::reduce_sum_f32(&[1.0, 2.0, 3.0]), 6.0);
        assert_eq!(cpu::reduce_sum_f32(&[]), 0.0);
        assert_eq!(cpu::reduce_max_f32(&[1.0, 5.0, 2.0]), 5.0);
        assert_eq!(cpu::reduce_max_f32(&[]), f32::NEG_INFINITY);
    }

    #[test]
    fn cpu_compact_basic() {
        let v = [10u32, 20, 30, 40];
        let k = [1u32, 0, 1, 0];
        assert_eq!(cpu::compact_u32(&v, &k), vec![10, 30]);
    }

    #[test]
    fn cpu_histogram_basic() {
        assert_eq!(cpu::histogram_u32(&[0, 1, 1, 2, 2, 2, 5], 4), vec![1, 2, 3, 1]);
        assert_eq!(cpu::histogram_u32(&[0, 1], 0), Vec::<u32>::new());
    }

    #[test]
    fn public_api_matches_cpu_smoke() {
        let data: Vec<u32> = (0..1000u32).map(|i| i % 7).collect();
        assert_eq!(exclusive_scan_u32(&data), cpu::exclusive_scan_u32(&data));
        let f: Vec<f32> = (0..1000u32).map(|i| (i % 5) as f32).collect();
        assert_eq!(reduce_max_f32(&f), cpu::reduce_max_f32(&f));
    }
}
