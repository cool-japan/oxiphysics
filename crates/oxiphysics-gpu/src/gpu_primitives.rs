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
//! One dispatch dimension holds at most 65 535 workgroups of 256 elements
//! (16 776 960 elements).  Longer inputs are not refused: their dispatch is
//! folded into a 2-D grid ([`crate::compute::timestamp::dispatch_count_for`])
//! and every kernel linearises its invocation and workgroup indices with
//! `num_workgroups`, so only the device's memory limits and the 32-bit index
//! space bound the input.  Histograms with more than 256 bins are computed by
//! the CPU reference implementation, and so is any call the GPU rejects (for
//! example an allocation beyond the device limits).
//!
//! All entry points of this module and of [`crate::gpu_radix`] share one
//! process-wide backend that is created on first use: the wgpu instance,
//! adapter and device and every compiled pipeline are reused by later calls
//! instead of being rebuilt each time, and a host without a GPU is probed only
//! once.  Every call releases the scratch buffers it allocated before it
//! returns, so the long-lived backend does not accumulate buffers.  Calls from
//! several threads are serialised on the backend.
//!
//! The CPU implementations in [`crate::gpu_reduction`] and
//! [`crate::grid_reduce`] act as the parity oracles for the GPU kernels.

/// Exclusive prefix sum over `u32` (GPU when available, else CPU).
///
/// Overflow wraps (matching the GPU's modular `u32` arithmetic).
pub fn exclusive_scan_u32(data: &[u32]) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        if data.len() >= 2
            && let Some(result) = shared::run(|gpu| gpu::exclusive_scan_u32(gpu, data))
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
        if !data.is_empty()
            && let Some(result) = shared::run(|gpu| gpu::reduce_sum_f32(gpu, data))
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
        if !data.is_empty()
            && let Some(result) = shared::run(|gpu| gpu::reduce_max_f32(gpu, data))
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
        if !values.is_empty()
            && let Some(result) = shared::run(|gpu| gpu::compact_u32(gpu, values, keep))
        {
            return result;
        }
    }
    cpu::compact_u32(values, keep)
}

/// Histogram over `u32`: value `v` increments bin `min(v, num_bins - 1)`.
///
/// The GPU kernel supports up to 256 bins; a histogram with more bins is
/// computed on the CPU.
pub fn histogram_u32(data: &[u32], num_bins: usize) -> Vec<u32> {
    #[cfg(feature = "wgpu-backend")]
    {
        // The histogram kernel keeps its bins in 256 workgroup-local counters.
        if !data.is_empty()
            && (1..=256).contains(&num_bins)
            && let Some(result) = shared::run(|gpu| gpu::histogram_u32(gpu, data, num_bins))
        {
            return result;
        }
    }
    cpu::histogram_u32(data, num_bins)
}

// ── CPU reference (parity oracle) ─────────────────────────────────────────────

pub(crate) mod cpu {
    /// Wrapping exclusive prefix sum; also the host-side scan of the radix sort.
    pub(crate) fn exclusive_scan_u32(data: &[u32]) -> Vec<u32> {
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

// ── Shared backend (real wgpu) ────────────────────────────────────────────────

/// The process-wide wgpu backend and the per-call scratch buffers used by the public
/// primitives of this module and of [`crate::gpu_radix`].
///
/// Creating a wgpu instance, adapter and device and compiling the WGSL pipelines is
/// orders of magnitude more expensive than the work of a primitive on a few thousand
/// elements, so the backend is created once, on first use, and kept for the life of
/// the process.  The pipeline cache of the backend then makes every later call a plain
/// upload / dispatch / readback.  A missing adapter is remembered as well: a host
/// without a GPU probes for one exactly once.
///
/// Calls serialise on the backend's mutex.  A closure given to [`shared::run`] holds
/// that mutex, so it must never call a public primitive itself (that would deadlock);
/// nested GPU work has to go through the [`shared::Scratch`] it receives.
#[cfg(feature = "wgpu-backend")]
pub(crate) mod shared {
    use crate::compute::timestamp::{MAX_WORKGROUPS_PER_DIMENSION, dispatch_count_with_limit};
    use crate::compute::wgpu_backend::real::WgpuBackendReal;
    use crate::compute::wgpu_backend::{WgpuBufferHandle, WgpuInitError};
    use std::ops::Deref;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Mutex, OnceLock, PoisonError};

    /// The one backend of the process; `None` records that no usable adapter exists.
    static BACKEND: OnceLock<Option<Mutex<WgpuBackendReal>>> = OnceLock::new();

    /// How many times the backend creation ran (at most once).
    static INIT_ATTEMPTS: AtomicUsize = AtomicUsize::new(0);

    /// How many primitive calls completed on the GPU.
    static GPU_RUNS: AtomicUsize = AtomicUsize::new(0);

    /// The backend shared by all GPU primitives, created on first use.
    ///
    /// Returns `None` when no GPU adapter or device is available; that outcome is
    /// cached, so callers can ask on every call without probing the host again.
    pub(crate) fn shared_backend() -> Option<&'static Mutex<WgpuBackendReal>> {
        BACKEND
            .get_or_init(|| {
                INIT_ATTEMPTS.fetch_add(1, Ordering::Relaxed);
                WgpuBackendReal::try_new().ok().map(Mutex::new)
            })
            .as_ref()
    }

    /// Run `body` on the shared backend inside a [`Scratch`] scope.
    ///
    /// Every scratch buffer `body` allocates is destroyed when it returns — also when
    /// it fails or panics — so the long-lived backend does not accumulate buffers.
    ///
    /// Returns `None` when there is no GPU or when `body` fails; the caller then
    /// takes its CPU path.  A mutex poisoned by a panic in an earlier call is
    /// recovered: the scope of that call has already released its buffers.
    pub(crate) fn run<T>(
        body: impl FnOnce(&mut Scratch<'_>) -> Result<T, WgpuInitError>,
    ) -> Option<T> {
        let mutex = shared_backend()?;
        let mut guard = mutex.lock().unwrap_or_else(PoisonError::into_inner);
        let mut scratch = Scratch::new(&mut guard);
        let result = body(&mut scratch).ok();
        if result.is_some() {
            GPU_RUNS.fetch_add(1, Ordering::Relaxed);
        }
        result
    }

    /// Number of times the shared backend was created: `1` once it has been asked
    /// for, whether or not an adapter was found.
    #[cfg(test)]
    pub(crate) fn init_attempts() -> usize {
        INIT_ATTEMPTS.load(Ordering::Relaxed)
    }

    /// Number of primitive calls that ran to completion on the GPU.
    #[cfg(test)]
    pub(crate) fn gpu_runs() -> usize {
        GPU_RUNS.load(Ordering::Relaxed)
    }

    /// Scope that owns the scratch buffers of one primitive call.
    ///
    /// Buffers come from [`alloc`](Self::alloc) and are destroyed when the scope is
    /// dropped.  Every other backend operation (uploads, dispatches) is available
    /// through `Deref`.
    pub(crate) struct Scratch<'a> {
        backend: &'a mut WgpuBackendReal,
        owned: Vec<WgpuBufferHandle>,
        /// Most workgroups one dispatch of this scope may span along a single
        /// dimension (see [`dispatch_shape`](Self::dispatch_shape)).
        workgroup_limit: u32,
    }

    impl<'a> Scratch<'a> {
        /// Open a scope over `backend`.
        pub(crate) fn new(backend: &'a mut WgpuBackendReal) -> Self {
            Self {
                backend,
                owned: Vec::new(),
                workgroup_limit: MAX_WORKGROUPS_PER_DIMENSION,
            }
        }

        /// Fold dispatches under a per-dimension limit of `limit` workgroups.
        ///
        /// Production scopes use the portable WebGPU limit; the tests lower it so
        /// that small inputs already need a multi-row (folded) grid.
        #[cfg(test)]
        pub(crate) fn with_workgroup_limit(mut self, limit: u32) -> Self {
            self.workgroup_limit = limit;
            self
        }

        /// Dispatch shape that covers `items` threads with workgroups of
        /// `group_size` threads.
        ///
        /// A shape wider than one dimension allows is folded into y; every kernel
        /// dispatched with it linearises its indices with `num_workgroups`
        /// (see [`dispatch_count_with_limit`]).  Fails only when even the folded
        /// grid cannot cover `items` or its thread indices would not fit 32 bits.
        pub(crate) fn dispatch_shape(
            &self,
            items: usize,
            group_size: u32,
        ) -> Result<[u32; 3], WgpuInitError> {
            Ok(dispatch_count_with_limit(
                items,
                group_size,
                self.workgroup_limit,
            )?)
        }

        /// Allocate a storage buffer of `size_bytes` bytes that lives until the scope ends.
        pub(crate) fn alloc(&mut self, size_bytes: u64) -> WgpuBufferHandle {
            let handle = self.backend.create_buffer_storage(size_bytes);
            self.owned.push(handle);
            handle
        }

        /// Read the first `words` `u32` values of `handle`.
        ///
        /// A readback that fails, or that comes back shorter than `words`, is an
        /// error - never silently padded or truncated data.  Reading zero words
        /// needs no readback and always succeeds.
        pub(crate) fn read_u32(
            &self,
            handle: WgpuBufferHandle,
            words: usize,
        ) -> Result<Vec<u32>, WgpuInitError> {
            if words == 0 {
                return Ok(Vec::new());
            }
            let mut data = self.backend.read_buffer_u32(handle)?;
            if data.len() < words {
                return Err(short_readback(handle, data.len(), words));
            }
            data.truncate(words);
            Ok(data)
        }

        /// Read the first `words` `f32` values of `handle`; see [`read_u32`](Self::read_u32).
        pub(crate) fn read_f32(
            &self,
            handle: WgpuBufferHandle,
            words: usize,
        ) -> Result<Vec<f32>, WgpuInitError> {
            if words == 0 {
                return Ok(Vec::new());
            }
            let mut data = self.backend.read_buffer_f32(handle)?;
            if data.len() < words {
                return Err(short_readback(handle, data.len(), words));
            }
            data.truncate(words);
            Ok(data)
        }
    }

    fn short_readback(handle: WgpuBufferHandle, got: usize, expected: usize) -> WgpuInitError {
        WgpuInitError::DeviceRequestFailed(format!(
            "readback of GPU buffer {} returned {got} words, expected {expected}",
            handle.0
        ))
    }

    impl Deref for Scratch<'_> {
        type Target = WgpuBackendReal;

        fn deref(&self) -> &WgpuBackendReal {
            self.backend
        }
    }

    impl Drop for Scratch<'_> {
        fn drop(&mut self) {
            for handle in self.owned.drain(..) {
                self.backend.destroy_buffer(handle);
            }
        }
    }
}

// ── GPU implementation (real wgpu) ────────────────────────────────────────────

#[cfg(feature = "wgpu-backend")]
mod gpu {
    use super::shared::Scratch;
    use crate::compute::wgpu_backend::WgpuInitError;
    use crate::kernels_wgsl::{
        COMPACT_WGSL, HISTOGRAM_WGSL, REDUCE_WGSL, SCAN_ADD_WGSL, SCAN_WGSL,
    };

    const BLOCK: usize = 256;

    /// The parameter words of a kernel: element count, bin count, two padding words.
    ///
    /// Counts reach the kernels as 32-bit words; a count that does not fit is
    /// refused instead of wrapping around.
    fn params_bytes(n: usize, num_bins: usize) -> Result<Vec<u8>, WgpuInitError> {
        let word = |value: usize| {
            u32::try_from(value).map_err(|_| {
                WgpuInitError::Validation(format!("{value} does not fit a 32-bit kernel parameter"))
            })
        };
        let p: [u32; 4] = [word(n)?, word(num_bins)?, 0, 0];
        Ok(bytemuck::cast_slice(&p).to_vec())
    }

    fn ro() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: true }
    }

    fn rw() -> wgpu::BufferBindingType {
        wgpu::BufferBindingType::Storage { read_only: false }
    }

    pub(super) fn exclusive_scan_u32(
        gpu: &mut Scratch<'_>,
        data: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = data.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        if n == 1 {
            return Ok(vec![0u32]);
        }

        // One workgroup per 256-element block; more blocks than one dispatch
        // dimension holds are folded into y (the kernels linearise the block index).
        let workgroups = gpu.dispatch_shape(n, BLOCK as u32)?;
        let num_blocks = n.div_ceil(BLOCK);
        let padded = num_blocks * BLOCK;

        let mut padded_in = vec![0u32; padded];
        padded_in[..n].copy_from_slice(data);

        let input_buf = gpu.alloc((padded * 4) as u64);
        let output_buf = gpu.alloc((padded * 4) as u64);
        let block_sums_buf = gpu.alloc((num_blocks.max(1) * 4) as u64);
        let params_buf = gpu.alloc(16);

        gpu.queue_write_buffer_raw(&input_buf, bytemuck::cast_slice(&padded_in));
        gpu.queue_write_buffer_raw(&params_buf, &params_bytes(n, 0)?);

        gpu.dispatch_wgsl(
            SCAN_WGSL,
            "scan_block",
            &[
                (input_buf, ro()),
                (output_buf, rw()),
                (block_sums_buf, rw()),
                (params_buf, ro()),
            ],
            workgroups,
        )?;

        if num_blocks > 1 {
            let block_sums = gpu.read_u32(block_sums_buf, num_blocks)?;
            let block_offsets = exclusive_scan_u32(gpu, &block_sums)?;

            let offsets_buf = gpu.alloc((num_blocks * 4) as u64);
            gpu.queue_write_buffer_raw(&offsets_buf, bytemuck::cast_slice(&block_offsets));

            gpu.dispatch_wgsl(
                SCAN_ADD_WGSL,
                "add_block_offsets",
                &[(output_buf, rw()), (offsets_buf, ro()), (params_buf, ro())],
                workgroups,
            )?;
        }

        gpu.read_u32(output_buf, n)
    }

    pub(super) fn reduce_sum_f32(
        gpu: &mut Scratch<'_>,
        data: &[f32],
    ) -> Result<f32, WgpuInitError> {
        if data.is_empty() {
            return Ok(0.0);
        }
        let mut current: Vec<f32> = data.to_vec();
        loop {
            let n = current.len();
            // One workgroup per block; a grid wider than one dimension is folded.
            let workgroups = gpu.dispatch_shape(n, BLOCK as u32)?;
            let num_blocks = n.div_ceil(BLOCK);
            let padded = num_blocks * BLOCK;

            let mut padded_in = vec![0f32; padded];
            padded_in[..n].copy_from_slice(&current);

            let in_buf = gpu.alloc((padded * 4) as u64);
            let partials_buf = gpu.alloc((num_blocks * 4) as u64);
            let params_buf = gpu.alloc(16);

            gpu.queue_write_buffer_raw(&in_buf, bytemuck::cast_slice(&padded_in));
            gpu.queue_write_buffer_raw(&params_buf, &params_bytes(n, 0)?);

            gpu.dispatch_wgsl(
                REDUCE_WGSL,
                "reduce_sum_block",
                &[(in_buf, ro()), (partials_buf, rw()), (params_buf, ro())],
                workgroups,
            )?;

            let partials = gpu.read_f32(partials_buf, num_blocks)?;
            if num_blocks <= 1 {
                return Ok(partials.first().copied().unwrap_or(0.0));
            }
            current = partials;
        }
    }

    pub(super) fn reduce_max_f32(
        gpu: &mut Scratch<'_>,
        data: &[f32],
    ) -> Result<f32, WgpuInitError> {
        if data.is_empty() {
            return Ok(f32::NEG_INFINITY);
        }
        const NEG_FLT_MAX: f32 = -3.4028235e38;
        let mut current: Vec<f32> = data.to_vec();
        loop {
            let n = current.len();
            // One workgroup per block; a grid wider than one dimension is folded.
            let workgroups = gpu.dispatch_shape(n, BLOCK as u32)?;
            let num_blocks = n.div_ceil(BLOCK);
            let padded = num_blocks * BLOCK;

            let mut padded_in = vec![NEG_FLT_MAX; padded];
            padded_in[..n].copy_from_slice(&current);

            let in_buf = gpu.alloc((padded * 4) as u64);
            let partials_buf = gpu.alloc((num_blocks * 4) as u64);
            let params_buf = gpu.alloc(16);

            gpu.queue_write_buffer_raw(&in_buf, bytemuck::cast_slice(&padded_in));
            gpu.queue_write_buffer_raw(&params_buf, &params_bytes(n, 0)?);

            gpu.dispatch_wgsl(
                REDUCE_WGSL,
                "reduce_max_block",
                &[(in_buf, ro()), (partials_buf, rw()), (params_buf, ro())],
                workgroups,
            )?;

            let partials = gpu.read_f32(partials_buf, num_blocks)?;
            if num_blocks <= 1 {
                return Ok(partials.first().copied().unwrap_or(f32::NEG_INFINITY));
            }
            current = partials;
        }
    }

    pub(super) fn compact_u32(
        gpu: &mut Scratch<'_>,
        values: &[u32],
        keep: &[u32],
    ) -> Result<Vec<u32>, WgpuInitError> {
        let n = values.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        // Fail before any work if even the folded grid cannot cover the input.
        let workgroups = gpu.dispatch_shape(n, BLOCK as u32)?;

        // The scan sums raw values and the kernel tests `keep == 1`, while the CPU
        // reference keeps every non-zero flag: normalise the flags to 0/1 so both
        // paths agree for any flag value.
        let flags: Vec<u32> = keep.iter().map(|&k| u32::from(k != 0)).collect();
        let positions = exclusive_scan_u32(gpu, &flags)?;
        let last_pos = positions.get(n - 1).copied().unwrap_or(0) as usize;
        let last_keep = flags.get(n - 1).copied().unwrap_or(0) as usize;
        let total_kept = last_pos + last_keep;
        if total_kept == 0 {
            return Ok(Vec::new());
        }

        let values_buf = gpu.alloc((n * 4) as u64);
        let keep_buf = gpu.alloc((n * 4) as u64);
        let positions_buf = gpu.alloc((n * 4) as u64);
        let output_buf = gpu.alloc((total_kept * 4) as u64);
        let params_buf = gpu.alloc(16);

        gpu.queue_write_buffer_raw(&values_buf, bytemuck::cast_slice(values));
        gpu.queue_write_buffer_raw(&keep_buf, bytemuck::cast_slice(&flags));
        gpu.queue_write_buffer_raw(&positions_buf, bytemuck::cast_slice(&positions));
        gpu.queue_write_buffer_raw(&params_buf, &params_bytes(n, 0)?);

        gpu.dispatch_wgsl(
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

        gpu.read_u32(output_buf, total_kept)
    }

    pub(super) fn histogram_u32(
        gpu: &mut Scratch<'_>,
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
        // Fail before any work if even the folded grid cannot cover the input.
        let workgroups = gpu.dispatch_shape(n, BLOCK as u32)?;

        let data_buf = gpu.alloc((n.max(1) * 4) as u64);
        let global_bins_buf = gpu.alloc((num_bins * 4) as u64);
        let params_buf = gpu.alloc(16);

        gpu.queue_write_buffer_raw(&data_buf, bytemuck::cast_slice(data));
        let zero_bins = vec![0u32; num_bins];
        gpu.queue_write_buffer_raw(&global_bins_buf, bytemuck::cast_slice(&zero_bins));
        gpu.queue_write_buffer_raw(&params_buf, &params_bytes(n, num_bins)?);

        gpu.dispatch_wgsl(
            HISTOGRAM_WGSL,
            "histogram_main",
            &[
                (data_buf, ro()),
                (global_bins_buf, rw()),
                (params_buf, ro()),
            ],
            workgroups,
        )?;

        gpu.read_u32(global_bins_buf, num_bins)
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
        assert_eq!(
            cpu::histogram_u32(&[0, 1, 1, 2, 2, 2, 5], 4),
            vec![1, 2, 3, 1]
        );
        assert_eq!(cpu::histogram_u32(&[0, 1], 0), Vec::<u32>::new());
    }

    #[test]
    fn public_api_matches_cpu_smoke() {
        let data: Vec<u32> = (0..1000u32).map(|i| i % 7).collect();
        assert_eq!(exclusive_scan_u32(&data), cpu::exclusive_scan_u32(&data));
        let f: Vec<f32> = (0..1000u32).map(|i| (i % 5) as f32).collect();
        assert_eq!(reduce_max_f32(&f), cpu::reduce_max_f32(&f));
    }

    /// Deterministic pseudo-random data (numerical-recipes LCG) in `0..modulus`.
    fn lcg_data(n: usize, modulus: u32) -> Vec<u32> {
        let mut state = 0x9E37_79B9u32;
        (0..n)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                state % modulus
            })
            .collect()
    }

    #[test]
    fn cpu_exclusive_scan_wraps_on_overflow() {
        assert_eq!(
            cpu::exclusive_scan_u32(&[u32::MAX, 2, 3]),
            vec![0, u32::MAX, 1]
        );
    }

    #[test]
    fn compact_treats_every_nonzero_flag_as_keep() {
        let values = [10u32, 20, 30, 40, 50];
        // The CPU reference keeps `keep[i] != 0`; the GPU path must agree even
        // for flags other than 1.
        assert_eq!(compact_u32(&values, &[2, 0, 7, 0, 1]), vec![10, 30, 50]);
        assert_eq!(compact_u32(&values, &[0, 0, 0, 0, 0]), Vec::<u32>::new());
        assert_eq!(compact_u32(&values, &[u32::MAX; 5]), values.to_vec());
        assert_eq!(
            compact_u32(&values, &[1, 2]),
            Vec::<u32>::new(),
            "length mismatch yields an empty vec"
        );
    }

    #[test]
    fn public_primitives_match_the_cpu_reference_across_block_boundaries() {
        for &n in &[0usize, 1, 2, 255, 256, 257, 511, 512, 513, 1000, 5000] {
            let data = lcg_data(n, 1000);
            assert_eq!(
                exclusive_scan_u32(&data),
                cpu::exclusive_scan_u32(&data),
                "scan n={n}"
            );
            let floats: Vec<f32> = data.iter().map(|&v| (v % 11) as f32).collect();
            assert_eq!(
                reduce_sum_f32(&floats),
                cpu::reduce_sum_f32(&floats),
                "sum n={n}"
            );
            assert_eq!(
                reduce_max_f32(&floats),
                cpu::reduce_max_f32(&floats),
                "max n={n}"
            );
            let keep: Vec<u32> = data.iter().map(|&v| u32::from(v % 3 == 0)).collect();
            assert_eq!(
                compact_u32(&data, &keep),
                cpu::compact_u32(&data, &keep),
                "compact n={n}"
            );
            for &bins in &[1usize, 7, 200, 256, 300] {
                assert_eq!(
                    histogram_u32(&data, bins),
                    cpu::histogram_u32(&data, bins),
                    "histogram n={n} bins={bins}"
                );
            }
        }
    }

    #[cfg(feature = "wgpu-backend")]
    mod shared_backend {
        use super::*;
        use crate::compute::wgpu_backend::{WgpuBufferHandle, WgpuInitError};
        use std::sync::PoisonError;

        #[test]
        fn backend_is_created_once_for_all_threads_and_calls() {
            let workers: Vec<_> = (0..4usize)
                .map(|t| {
                    std::thread::spawn(move || {
                        for round in 0..3usize {
                            let data = lcg_data(3000 + 17 * t + round, 97);
                            assert_eq!(exclusive_scan_u32(&data), cpu::exclusive_scan_u32(&data));
                            assert_eq!(histogram_u32(&data, 97), cpu::histogram_u32(&data, 97));
                        }
                    })
                })
                .collect();
            for worker in workers {
                worker.join().expect("worker thread panicked");
            }
            let first = shared::shared_backend().map(|b| b as *const _);
            let second = shared::shared_backend().map(|b| b as *const _);
            assert_eq!(first, second, "every call must see the same backend");
            assert_eq!(
                shared::init_attempts(),
                1,
                "the backend (or the failure to find one) must be created exactly once"
            );
        }

        #[test]
        fn pipelines_are_compiled_once_and_reused_by_later_calls() {
            // The whole measurement runs inside one scope.  It holds the lock of the
            // shared backend, so no other test thread can compile a pipeline in between.
            let measured = shared::run(|gpu| {
                let data = lcg_data(3000, 100);
                let floats: Vec<f32> = data.iter().map(|&v| v as f32).collect();
                let keep: Vec<u32> = data.iter().map(|&v| u32::from(v % 2 == 0)).collect();
                // Warm-up: one call of every kernel compiles all of its pipelines.
                gpu::exclusive_scan_u32(gpu, &data)?;
                gpu::reduce_sum_f32(gpu, &floats)?;
                gpu::reduce_max_f32(gpu, &floats)?;
                gpu::compact_u32(gpu, &data, &keep)?;
                gpu::histogram_u32(gpu, &data, 100)?;
                let warm = gpu.compiled_pipeline_count();
                // Later calls with other sizes (including a multi-level scan that needs
                // the block-offset pass) must find every pipeline in the cache.
                for n in [1000usize, 5000, 70_000] {
                    let data = lcg_data(n, 100);
                    let floats: Vec<f32> = data.iter().map(|&v| v as f32).collect();
                    let keep: Vec<u32> = data.iter().map(|&v| u32::from(v % 3 == 0)).collect();
                    gpu::exclusive_scan_u32(gpu, &data)?;
                    gpu::reduce_sum_f32(gpu, &floats)?;
                    gpu::reduce_max_f32(gpu, &floats)?;
                    gpu::compact_u32(gpu, &data, &keep)?;
                    gpu::histogram_u32(gpu, &data, 100)?;
                }
                Ok::<_, WgpuInitError>((warm, gpu.compiled_pipeline_count()))
            });
            let Some((warm, after_reuse)) = measured else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            // scan_block, add_block_offsets, reduce_sum_block, reduce_max_block,
            // scatter_kept and histogram_main.
            assert!(warm >= 6, "expected at least six pipelines, got {warm}");
            assert_eq!(
                warm, after_reuse,
                "repeated calls must reuse the cached pipelines instead of recompiling"
            );
        }

        #[test]
        fn scratch_buffers_are_released_after_every_call() {
            let Some(backend) = shared::shared_backend() else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            let runs_before = shared::gpu_runs();
            let data = lcg_data(5000, 255);
            let floats: Vec<f32> = data.iter().map(|&v| v as f32).collect();
            let keep: Vec<u32> = data.iter().map(|&v| u32::from(v % 2 == 0)).collect();
            for _ in 0..3 {
                assert_eq!(exclusive_scan_u32(&data), cpu::exclusive_scan_u32(&data));
                assert_eq!(reduce_sum_f32(&floats), cpu::reduce_sum_f32(&floats));
                assert_eq!(reduce_max_f32(&floats), cpu::reduce_max_f32(&floats));
                assert_eq!(compact_u32(&data, &keep), cpu::compact_u32(&data, &keep));
                assert_eq!(histogram_u32(&data, 255), cpu::histogram_u32(&data, 255));
            }
            assert!(
                shared::gpu_runs() >= runs_before + 15,
                "all fifteen calls must have run on the GPU"
            );
            let guard = backend.lock().unwrap_or_else(PoisonError::into_inner);
            assert_eq!(guard.live_buffer_count(), 0, "scratch buffers leaked");
            assert_eq!(guard.live_buffer_bytes(), 0, "scratch buffers leaked");
        }

        #[test]
        fn scope_releases_buffers_when_the_body_panics_and_the_lock_recovers() {
            let Some(backend) = shared::shared_backend() else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            let outcome = std::thread::spawn(|| {
                shared::run(|gpu| -> Result<(), WgpuInitError> {
                    let _scratch_buffer = gpu.alloc(4096);
                    panic!("simulated failure inside a GPU scope");
                })
            })
            .join();
            assert!(outcome.is_err(), "the panic must propagate to the joiner");
            // The panic poisoned the mutex; later calls recover it and still work.
            let data = lcg_data(1000, 50);
            let runs_before = shared::gpu_runs();
            assert_eq!(exclusive_scan_u32(&data), cpu::exclusive_scan_u32(&data));
            assert!(shared::gpu_runs() > runs_before);
            let guard = backend.lock().unwrap_or_else(PoisonError::into_inner);
            assert_eq!(guard.live_buffer_count(), 0, "the unwound scope leaked");
        }

        /// Run `body` in a scratch scope of the shared backend whose dispatches are
        /// folded under `limit` workgroups per dimension.
        ///
        /// `None` when the host has no GPU adapter.  Unlike the public primitives this
        /// surfaces every GPU error, so a result that comes back really was computed
        /// by the GPU and not by the CPU fallback.
        fn on_gpu<T>(
            limit: u32,
            body: impl FnOnce(&mut shared::Scratch<'_>) -> Result<T, WgpuInitError>,
        ) -> Option<Result<T, WgpuInitError>> {
            let mutex = shared::shared_backend()?;
            let mut guard = mutex.lock().unwrap_or_else(PoisonError::into_inner);
            let mut scratch = shared::Scratch::new(&mut guard).with_workgroup_limit(limit);
            Some(body(&mut scratch))
        }

        #[test]
        fn dispatch_shape_folds_exactly_beyond_the_workgroup_limit() {
            let shapes = on_gpu(3, |gpu| {
                Ok([0usize, 1, 700, 768, 769, 1000, 1100, 2000, 2304]
                    .map(|items| gpu.dispatch_shape(items, 256).expect("shape")))
            });
            let Some(Ok(shapes)) = shapes else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            assert_eq!(
                shapes,
                [
                    [0, 1, 1], // nothing to do
                    [1, 1, 1],
                    [3, 1, 1], // 3 workgroups still fit one dimension
                    [3, 1, 1],
                    [2, 2, 1], // 4 workgroups: folded
                    [2, 2, 1],
                    [3, 2, 1], // 5 workgroups in a 3 x 2 grid with one spare
                    [3, 3, 1], // 8 workgroups in a 3 x 3 grid with one spare
                    [3, 3, 1], // 9 workgroups fill the grid
                ]
            );
            // One workgroup more than a 3 x 3 grid holds is an error, not a panic.
            let beyond = on_gpu(3, |gpu| Ok(gpu.dispatch_shape(2305, 256)));
            assert!(matches!(beyond, Some(Ok(Err(_)))), "{beyond:?}");
            // At the portable limit the same inputs are not folded at all.
            let portable = on_gpu(shared_limit(), |gpu| gpu.dispatch_shape(5000, 256));
            assert!(matches!(portable, Some(Ok([20, 1, 1]))), "{portable:?}");
        }

        /// The limit production scopes use.
        fn shared_limit() -> u32 {
            crate::compute::timestamp::MAX_WORKGROUPS_PER_DIMENSION
        }

        #[test]
        fn folded_primitives_match_the_cpu_reference_when_the_grid_is_wider_than_one_row() {
            // A limit of 3 workgroups per dimension folds every input of 4 to 9
            // blocks into a 2-D grid, most of them with spare workgroups that lie
            // behind the data.
            for &n in &[257usize, 700, 1000, 1100, 1280, 2000, 2304] {
                let data = lcg_data(n, 1000);
                let floats: Vec<f32> = data.iter().map(|&v| (v % 13) as f32).collect();
                let keep: Vec<u32> = data.iter().map(|&v| u32::from(v % 3 == 0)).collect();
                let outcome = on_gpu(3, |gpu| {
                    let scan = gpu::exclusive_scan_u32(gpu, &data)?;
                    let sum = gpu::reduce_sum_f32(gpu, &floats)?;
                    let max = gpu::reduce_max_f32(gpu, &floats)?;
                    let compacted = gpu::compact_u32(gpu, &data, &keep)?;
                    let bins = gpu::histogram_u32(gpu, &data, 200)?;
                    Ok((scan, sum, max, compacted, bins))
                });
                let Some(outcome) = outcome else {
                    eprintln!("SKIPPED: no GPU adapter");
                    return;
                };
                let (scan, sum, max, compacted, bins) =
                    outcome.unwrap_or_else(|e| panic!("n={n}: the folded dispatch failed: {e}"));
                assert_eq!(scan, cpu::exclusive_scan_u32(&data), "scan n={n}");
                assert_eq!(sum, cpu::reduce_sum_f32(&floats), "sum n={n}");
                assert_eq!(max, cpu::reduce_max_f32(&floats), "max n={n}");
                assert_eq!(compacted, cpu::compact_u32(&data, &keep), "compact n={n}");
                assert_eq!(bins, cpu::histogram_u32(&data, 200), "histogram n={n}");
            }
        }

        #[test]
        fn folded_primitives_reject_inputs_beyond_the_folded_grid_without_allocating() {
            // 10 blocks do not fit a 3 x 3 grid.
            let data = vec![1u32; 256 * 10];
            let floats = vec![1.0f32; 256 * 10];
            let refused = on_gpu(3, |gpu| {
                let scan = gpu::exclusive_scan_u32(gpu, &data).is_err();
                let hist = gpu::histogram_u32(gpu, &data, 4).is_err();
                let keep = gpu::compact_u32(gpu, &data, &data).is_err();
                let sum = gpu::reduce_sum_f32(gpu, &floats).is_err();
                let max = gpu::reduce_max_f32(gpu, &floats).is_err();
                Ok((scan, hist, keep, sum, max, gpu.live_buffer_count()))
            });
            let Some(refused) = refused else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            assert_eq!(
                refused.expect("the closure never fails"),
                (true, true, true, true, true, 0),
                "an input the grid cannot cover must be refused before allocating anything"
            );
        }

        #[test]
        fn inputs_wider_than_one_dispatch_dimension_run_on_the_gpu() {
            // 65 537 blocks of 256 elements: one more than 65 536, which already
            // needs a second row of the portable 65 535-wide dispatch grid, and the
            // last block is partial.
            let n = 256 * 65_536 + 1;
            let data: Vec<u32> = (0..n as u32).map(|i| i % 4).collect();
            let limit = shared_limit();
            let scan = on_gpu(limit, |gpu| gpu::exclusive_scan_u32(gpu, &data));
            let Some(scan) = scan else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            assert_eq!(
                scan.expect("the scan must run on the GPU"),
                cpu::exclusive_scan_u32(&data),
                "scan of {n} elements"
            );
            let hist = on_gpu(limit, |gpu| gpu::histogram_u32(gpu, &data, 4));
            assert_eq!(
                hist.expect("adapter")
                    .expect("the histogram must run on the GPU"),
                cpu::histogram_u32(&data, 4),
                "histogram of {n} elements"
            );
            // 65 537 full blocks: every partial sum is exactly representable in f32
            // (16 777 472 = 2^24 + 256), so the order of summation cannot matter.
            let full_blocks = 256 * 65_537;
            let floats = vec![1.0f32; full_blocks];
            let sum = on_gpu(limit, |gpu| gpu::reduce_sum_f32(gpu, &floats));
            assert_eq!(
                sum.expect("adapter")
                    .expect("the reduction must run on the GPU"),
                16_777_472.0,
                "sum of {full_blocks} ones"
            );
        }

        #[test]
        fn short_readback_is_an_error_not_silent_data() {
            let seen = shared::run(|gpu| {
                Ok::<_, WgpuInitError>((
                    gpu.read_u32(WgpuBufferHandle(usize::MAX), 4).is_err(),
                    gpu.read_f32(WgpuBufferHandle(usize::MAX), 1).is_err(),
                    gpu.read_u32(WgpuBufferHandle(usize::MAX), 0).is_ok(),
                ))
            });
            match seen {
                Some(flags) => assert_eq!(flags, (true, true, true)),
                None => eprintln!("SKIPPED: no GPU adapter"),
            }
            // `run` turns a failing body into `None` so that callers take the CPU path.
            if shared::shared_backend().is_some() {
                let failed = shared::run(|gpu| gpu.read_u32(WgpuBufferHandle(usize::MAX), 4));
                assert!(failed.is_none());
            }
        }

        #[test]
        fn readback_shorter_than_requested_is_an_error_naming_both_lengths() {
            let outcome = shared::run(|gpu| {
                let buffer = gpu.alloc(16);
                gpu.queue_write_buffer_raw(&buffer, bytemuck::cast_slice(&[1u32, 2, 3, 4]));
                let prefix = gpu.read_u32(buffer, 2)?;
                let whole = gpu.read_u32(buffer, 4)?;
                let float_count = gpu.read_f32(buffer, 4)?.len();
                let too_many_words = gpu.read_u32(buffer, 5).map(|words| words.len());
                let too_many_floats = gpu.read_f32(buffer, 8).map(|floats| floats.len());
                Ok::<_, WgpuInitError>((
                    prefix,
                    whole,
                    float_count,
                    too_many_words,
                    too_many_floats,
                ))
            });
            let Some((prefix, whole, float_count, too_many_words, too_many_floats)) = outcome
            else {
                eprintln!("SKIPPED: no GPU adapter");
                return;
            };
            assert_eq!(prefix, vec![1, 2], "a shorter request reads a prefix");
            assert_eq!(whole, vec![1, 2, 3, 4]);
            assert_eq!(float_count, 4);
            for result in [too_many_words, too_many_floats] {
                match result {
                    Err(WgpuInitError::DeviceRequestFailed(message)) => {
                        assert!(
                            message.contains("returned 4 words, expected"),
                            "the error must name both lengths: {message}"
                        );
                    }
                    other => panic!("expected a short-readback error, got {other:?}"),
                }
            }
        }
    }
}
