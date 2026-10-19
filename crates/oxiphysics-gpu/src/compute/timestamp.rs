// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-side dispatch timing utilities for the GPU compute layer.
//!
//! [`ComputeDispatchTimer`] measures wall-clock time around a compute dispatch
//! using [`std::time::Instant`].  On platforms where GPU timestamp queries are
//! available (via the `wgpu-backend` feature), the companion [`GpuTimestamp`]
//! type can record hardware-level start/end values from a timestamp query set.

use std::time::Instant;

/// Wall-clock timer for a single compute dispatch.
///
/// Usage:
/// ```
/// use oxiphysics_gpu::compute::timestamp::ComputeDispatchTimer;
///
/// let start = ComputeDispatchTimer::start_cpu();
/// // ... do work ...
/// let timer = ComputeDispatchTimer::stop_cpu(start);
/// if let Some(ms) = timer.elapsed_ms() {
///     println!("dispatch took {ms:.3} ms");
/// }
/// ```
#[derive(Debug, Clone, Default)]
pub struct ComputeDispatchTimer {
    /// Elapsed time in nanoseconds, set after [`stop_cpu`](Self::stop_cpu).
    pub elapsed_ns: Option<u64>,
}

impl ComputeDispatchTimer {
    /// Create a new, unstopped timer.
    pub fn new() -> Self {
        Self { elapsed_ns: None }
    }

    /// Record the current instant as the start of a timed region.
    ///
    /// Pass the returned `Instant` to [`stop_cpu`](Self::stop_cpu) when the
    /// region ends.
    pub fn start_cpu() -> Instant {
        Instant::now()
    }

    /// Stop the timer and return a `ComputeDispatchTimer` with elapsed time.
    pub fn stop_cpu(start: Instant) -> Self {
        Self {
            elapsed_ns: Some(start.elapsed().as_nanos() as u64),
        }
    }

    /// Return elapsed time in milliseconds, or `None` if the timer was never stopped.
    pub fn elapsed_ms(&self) -> Option<f64> {
        self.elapsed_ns.map(|ns| ns as f64 / 1_000_000.0)
    }

    /// Return elapsed time in microseconds, or `None` if the timer was never stopped.
    pub fn elapsed_us(&self) -> Option<f64> {
        self.elapsed_ns.map(|ns| ns as f64 / 1_000.0)
    }

    /// Return `true` if the timer has been stopped and holds a measurement.
    pub fn has_measurement(&self) -> bool {
        self.elapsed_ns.is_some()
    }
}

/// GPU hardware timestamp pair (start / end), feature-gated to `wgpu-backend`.
///
/// In production use, these values are read back from a `wgpu::QuerySet` of
/// type `Timestamp` after the GPU has signalled completion.  The raw values
/// are in nanoseconds (after scaling by the adapter's timestamp period).
#[cfg(feature = "wgpu-backend")]
#[derive(Debug, Clone, Copy, Default)]
pub struct GpuTimestamp {
    /// Raw GPU timestamp at the start of the pass (nanoseconds).
    pub start_ns: u64,
    /// Raw GPU timestamp at the end of the pass (nanoseconds).
    pub end_ns: u64,
}

#[cfg(feature = "wgpu-backend")]
impl GpuTimestamp {
    /// Create a new GPU timestamp pair.
    pub fn new(start_ns: u64, end_ns: u64) -> Self {
        Self { start_ns, end_ns }
    }

    /// Elapsed GPU time in nanoseconds (saturating subtraction).
    pub fn elapsed_ns(&self) -> u64 {
        self.end_ns.saturating_sub(self.start_ns)
    }

    /// Elapsed GPU time in milliseconds.
    pub fn elapsed_ms(&self) -> f64 {
        self.elapsed_ns() as f64 / 1_000_000.0
    }
}

/// Largest workgroup count one dispatch dimension may have on every WebGPU
/// implementation (`maxComputeWorkgroupsPerDimension` of the default limits).
pub const MAX_WORKGROUPS_PER_DIMENSION: u32 = 65_535;

/// Compute the dispatch shape that covers `n_items` with workgroups of
/// `workgroup_size` threads.
///
/// Returns `[0, 1, 1]` when `n_items` is zero to produce a no-op dispatch
/// without panicking.  Up to [`MAX_WORKGROUPS_PER_DIMENSION`] workgroups
/// (`65535 × workgroup_size` items) the shape is the plain `[x, 1, 1]`.  Larger
/// inputs are *folded* into the y dimension — `y = ceil(groups / 65535)`,
/// `x = ceil(groups / y)` — so that no dimension exceeds the limit; the shape
/// then holds up to `65535² × workgroup_size` threads, a few of which lie past
/// the end of the data.
///
/// A kernel launched with a folded shape must linearise its thread and workgroup
/// indices with `@builtin(num_workgroups)`:
///
/// ```wgsl
/// let idx = gid.x + gid.y * nwg.x * WORKGROUP_SIZE;   // global invocation
/// let block = wid.x + wid.y * nwg.x;                  // workgroup
/// ```
///
/// and guard `idx < n` (or `block < num_blocks`).  For `y == 1` both formulas
/// reduce to the unfolded indices.  Passing `workgroup_size = 1` folds a plain
/// workgroup count.
///
/// This is a pure helper used by both the feature-gated real backend and any
/// CPU-side utilities that need to replicate the same dispatch sizing logic.
///
/// # Errors
///
/// [`Error::DispatchLimitExceeded`] when the folded shape would need more than
/// 65535 workgroups in y, when the threads of the shape cannot be addressed by
/// a 32-bit index (more than `u32::MAX + 1` threads), or when `n_items` itself
/// exceeds `u32::MAX` (the kernels receive the item count as a 32-bit word).
///
/// # Examples
///
/// ```
/// use oxiphysics_gpu::compute::timestamp::dispatch_count_for;
///
/// assert_eq!(dispatch_count_for(0, 64).unwrap(), [0, 1, 1]);
/// assert_eq!(dispatch_count_for(64, 64).unwrap(), [1, 1, 1]);
/// assert_eq!(dispatch_count_for(65, 64).unwrap(), [2, 1, 1]);
/// // 5 million items need 78 125 workgroups of 64: folded into a 2-D grid.
/// let [x, y, z] = dispatch_count_for(5_000_000, 64).unwrap();
/// assert!(x <= 65_535 && y <= 65_535 && z == 1);
/// assert!(u64::from(x) * u64::from(y) * 64 >= 5_000_000);
/// ```
pub fn dispatch_count_for(n_items: usize, workgroup_size: u32) -> crate::Result<[u32; 3]> {
    dispatch_count_with_limit(n_items, workgroup_size, MAX_WORKGROUPS_PER_DIMENSION)
}

/// Like [`dispatch_count_for`] with an explicit per-dimension workgroup limit.
///
/// Use this with `max_compute_workgroups_per_dimension` of the actual device
/// when it is known to exceed the portable [`MAX_WORKGROUPS_PER_DIMENSION`].  A
/// `max_per_dimension` of zero is treated as one.
///
/// # Errors
///
/// See [`dispatch_count_for`].
pub fn dispatch_count_with_limit(
    n_items: usize,
    workgroup_size: u32,
    max_per_dimension: u32,
) -> crate::Result<[u32; 3]> {
    if n_items == 0 {
        return Ok([0, 1, 1]);
    }
    // The kernels take the item count as a 32-bit word: a longer input cannot be
    // processed whatever the shape of the dispatch.
    if u32::try_from(n_items).is_err() {
        return Err(crate::Error::DispatchLimitExceeded {
            dispatch_size: n_items,
            limit: u32::MAX as usize,
        });
    }
    let threads_per_group = u64::from(workgroup_size.max(1));
    let limit = u64::from(max_per_dimension.max(1));
    let groups = (n_items as u64).div_ceil(threads_per_group);
    let (x, y) = if groups <= limit {
        (groups, 1)
    } else {
        let y = groups.div_ceil(limit);
        (groups.div_ceil(y), y)
    };
    let addressable_threads = u64::from(u32::MAX) + 1;
    if y > limit || x.saturating_mul(y).saturating_mul(threads_per_group) > addressable_threads {
        return Err(crate::Error::DispatchLimitExceeded {
            dispatch_size: usize::try_from(groups).unwrap_or(usize::MAX),
            limit: usize::try_from(limit.saturating_mul(limit)).unwrap_or(usize::MAX),
        });
    }
    // Both factors are at most `limit <= u32::MAX`, so the casts are lossless.
    Ok([x as u32, y as u32, 1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_cpu_timer_new_has_no_measurement() {
        let t = ComputeDispatchTimer::new();
        assert!(!t.has_measurement());
        assert!(t.elapsed_ms().is_none());
    }

    #[test]
    fn test_cpu_timer_stop_records_elapsed() {
        // Use a tiny sleep to ensure elapsed > 0.
        let start = ComputeDispatchTimer::start_cpu();
        std::thread::sleep(Duration::from_millis(1));
        let timer = ComputeDispatchTimer::stop_cpu(start);
        assert!(timer.has_measurement());
        let ns = timer.elapsed_ns.unwrap();
        assert!(ns > 0, "elapsed_ns should be > 0, got {ns}");
    }

    #[test]
    fn test_cpu_timer_elapsed_ms_positive() {
        let start = ComputeDispatchTimer::start_cpu();
        std::thread::sleep(Duration::from_millis(1));
        let timer = ComputeDispatchTimer::stop_cpu(start);
        let ms = timer.elapsed_ms().unwrap();
        assert!(ms > 0.0, "elapsed_ms should be positive, got {ms}");
    }

    #[test]
    fn test_dispatch_count_for_zero() {
        assert_eq!(dispatch_count_for(0, 64).unwrap(), [0, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_exact() {
        assert_eq!(dispatch_count_for(64, 64).unwrap(), [1, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_overflow() {
        assert_eq!(dispatch_count_for(65, 64).unwrap(), [2, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_one() {
        assert_eq!(dispatch_count_for(1, 64).unwrap(), [1, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_zero_workgroup_size_is_clamped() {
        assert_eq!(dispatch_count_for(5, 0).unwrap(), [5, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_stays_one_dimensional_up_to_the_limit() {
        // 65535 workgroups of 64 threads = 4_194_240 items fit in x.
        assert_eq!(dispatch_count_for(4_194_240, 64).unwrap(), [65_535, 1, 1]);
        // One item more needs a 65 536th workgroup -> folded.
        let [x, y, z] = dispatch_count_for(4_194_241, 64).unwrap();
        assert!(y == 2 && z == 1 && x <= 65_535, "{x} x {y}");
    }

    #[test]
    fn test_dispatch_count_for_folds_five_million_items() {
        let [x, y, z] = dispatch_count_for(5_000_000, 64).unwrap();
        assert!(x <= MAX_WORKGROUPS_PER_DIMENSION);
        assert!(y <= MAX_WORKGROUPS_PER_DIMENSION);
        assert_eq!(z, 1);
        assert!(y > 1, "78 125 workgroups cannot fit in x");
        // Every item is covered, and the grid is not wastefully large.
        let threads = u64::from(x) * u64::from(y) * 64;
        assert!(threads >= 5_000_000);
        assert!(threads < 5_000_000 + 64 * u64::from(y));
    }

    #[test]
    fn test_dispatch_count_for_folded_shape_covers_every_size() {
        for &(n, ws) in &[
            (4_194_241usize, 64u32),
            (16_776_961, 256),
            (100_000_000, 256),
            (1_000_000_007, 64),
            (4_000_000_000, 256),
        ] {
            let [x, y, z] = dispatch_count_for(n, ws).unwrap();
            assert!((1..=MAX_WORKGROUPS_PER_DIMENSION).contains(&x), "x for {n}");
            assert!((1..=MAX_WORKGROUPS_PER_DIMENSION).contains(&y), "y for {n}");
            assert_eq!(z, 1);
            let needed = (n as u64).div_ceil(u64::from(ws));
            let groups = u64::from(x) * u64::from(y);
            assert!(groups >= needed, "coverage for {n}");
            // The fold wastes fewer than one workgroup per row.
            assert!(groups - needed < u64::from(y), "waste for {n}");
            // Every thread index fits a 32-bit word.
            assert!(
                groups * u64::from(ws) <= u64::from(u32::MAX) + 1,
                "index space for {n}"
            );
        }
    }

    #[test]
    fn test_dispatch_count_for_rejects_sizes_whose_thread_index_would_wrap() {
        // 2^32 - 1 items in 256-thread groups need 16_777_217 workgroups in the
        // 65535-wide fold: the last threads would index past 2^32 and wrap.
        assert!(dispatch_count_for(4_294_967_295, 256).is_err());
    }

    #[test]
    fn test_dispatch_count_with_limit_honours_a_custom_limit() {
        assert_eq!(dispatch_count_with_limit(10, 1, 10).unwrap(), [10, 1, 1]);
        // 11 groups with a limit of 4 -> y = 3, x = 4.
        assert_eq!(dispatch_count_with_limit(11, 1, 4).unwrap(), [4, 3, 1]);
        assert_eq!(dispatch_count_with_limit(16, 1, 4).unwrap(), [4, 4, 1]);
        assert!(dispatch_count_with_limit(17, 1, 4).is_err());
        // A zero limit behaves like a limit of one.
        assert_eq!(dispatch_count_with_limit(1, 1, 0).unwrap(), [1, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_rejects_oversized_inputs() {
        // 2^32 items of one thread each cannot be addressed by a 32-bit index.
        let too_many = usize::try_from(u64::from(u32::MAX) + 2).unwrap_or(usize::MAX);
        match dispatch_count_for(too_many, 1) {
            Err(crate::Error::DispatchLimitExceeded {
                dispatch_size,
                limit,
            }) => {
                assert!(dispatch_size > 0 && limit > 0);
            }
            other => panic!("expected DispatchLimitExceeded, got {other:?}"),
        }
        // More workgroups than a 65535 x 65535 grid holds.
        let beyond_grid = 65_535usize * 65_535 + 1;
        assert!(dispatch_count_for(beyond_grid, 1).is_err());
        // The index-space guard also catches the rounding of the last row.
        assert!(dispatch_count_for(usize::MAX, 256).is_err());
    }

    #[test]
    fn test_dispatch_count_for_error_names_the_limit() {
        let err = dispatch_count_for(usize::MAX, 1).unwrap_err();
        let message = err.to_string();
        assert!(message.contains("exceeds hardware limit"), "{message}");
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn test_dispatch_count_rejects_item_counts_beyond_a_32_bit_word() {
        // With a huge custom limit the grid shape alone would be addressable
        // (2 x 2^31 threads), yet the kernels cannot take 2^32 items as their
        // 32-bit count parameter.
        let beyond = usize::try_from(u64::from(u32::MAX) + 1).unwrap_or(usize::MAX);
        assert!(dispatch_count_with_limit(beyond, 1, u32::MAX).is_err());
        // The largest countable input is still dispatched.
        let largest = u32::MAX as usize;
        assert_eq!(
            dispatch_count_with_limit(largest, 1, u32::MAX).unwrap(),
            [u32::MAX, 1, 1]
        );
    }

    #[cfg(feature = "wgpu-backend")]
    #[test]
    fn test_gpu_timestamp_elapsed() {
        let ts = GpuTimestamp::new(1000, 5000);
        assert_eq!(ts.elapsed_ns(), 4000);
        assert!((ts.elapsed_ms() - 0.004).abs() < 1e-9);
    }

    #[cfg(feature = "wgpu-backend")]
    #[test]
    fn test_gpu_timestamp_saturating_sub() {
        // start > end should not panic
        let ts = GpuTimestamp::new(5000, 1000);
        assert_eq!(ts.elapsed_ns(), 0);
    }
}
