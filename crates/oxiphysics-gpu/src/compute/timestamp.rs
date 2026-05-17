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

/// Compute the number of workgroups needed to cover `n_items` in the X dimension.
///
/// Returns `[0, 1, 1]` when `n_items` is zero to produce a no-op dispatch
/// without panicking.
///
/// This is a pure helper used by both the feature-gated real backend and any
/// CPU-side utilities that need to replicate the same dispatch sizing logic.
///
/// # Examples
///
/// ```
/// use oxiphysics_gpu::compute::timestamp::dispatch_count_for;
///
/// assert_eq!(dispatch_count_for(0, 64), [0, 1, 1]);
/// assert_eq!(dispatch_count_for(64, 64), [1, 1, 1]);
/// assert_eq!(dispatch_count_for(65, 64), [2, 1, 1]);
/// ```
pub fn dispatch_count_for(n_items: usize, workgroup_size: u32) -> [u32; 3] {
    if n_items == 0 {
        return [0, 1, 1];
    }
    let ws = workgroup_size.max(1);
    let x = (n_items as u32).div_ceil(ws);
    [x, 1, 1]
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
        assert_eq!(dispatch_count_for(0, 64), [0, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_exact() {
        assert_eq!(dispatch_count_for(64, 64), [1, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_overflow() {
        assert_eq!(dispatch_count_for(65, 64), [2, 1, 1]);
    }

    #[test]
    fn test_dispatch_count_for_one() {
        assert_eq!(dispatch_count_for(1, 64), [1, 1, 1]);
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
