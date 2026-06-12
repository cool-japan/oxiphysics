// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compensated and deterministic floating-point summation algorithms.
//!
//! This module provides numerically stable summation routines based on the
//! Neumaier (improved Kahan-Babuška) algorithm and pairwise (tree) summation,
//! as well as a [`DeterministicReducer`] that produces bit-identical results
//! regardless of the number of parallel threads used.

use std::sync::{Arc, Mutex};

use crate::parallel::WorkStealingPool;

/// Computes the Neumaier (improved Kahan-Babuška) compensated sum of `xs`.
///
/// This algorithm corrects for catastrophic cancellation in both directions:
/// when the running sum dominates the incoming value *and* when the incoming
/// value dominates the running sum.  It is strictly more accurate than the
/// classic Kahan summation algorithm.
///
/// # Algorithm
///
/// For each element `x` in `xs`:
/// - If `|sum| >= |x|`: `x` is small relative to `sum`; the low bits of `x`
///   are lost in `sum + x`, so we track the error in `compensation`.
/// - Otherwise: `sum` is small relative to `x`; the low bits of `sum` are
///   lost in `x + sum`, so we track the error in `compensation`.
///
/// The final result is `sum + compensation`.
///
/// # Examples
///
/// ```
/// # use oxiphysics_core::compensated::neumaier_sum;
/// // Naive IEEE 754: 1e16 + 1.0 - 1e16 = 0.0 due to float cancellation
/// assert_eq!(neumaier_sum(&[1e16_f64, 1.0, -1e16]), 1.0);
/// ```
pub fn neumaier_sum(xs: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    let mut compensation = 0.0_f64;

    for &x in xs {
        let t = sum + x;
        if sum.abs() >= x.abs() {
            // `x` is small; low bits of `x` lost when added to `sum`
            compensation += (sum - t) + x;
        } else {
            // `sum` is small; low bits of `sum` lost when added to `x`
            compensation += (x - t) + sum;
        }
        sum = t;
    }

    sum + compensation
}

/// Computes the pairwise (recursive tree) sum of `xs`.
///
/// Pairwise summation divides the input recursively into halves, summing
/// each half independently and then combining the results.  For inputs longer
/// than 128 elements the recursion bottoms out at [`neumaier_sum`], combining
/// both accuracy and efficiency.
///
/// The error bound of pairwise summation grows as O(log n) rather than the
/// O(n) bound of naive summation, making it well-suited for large arrays.
///
/// # Examples
///
/// ```
/// # use oxiphysics_core::compensated::pairwise_sum;
/// let xs: Vec<f64> = (1..=100).map(|i| i as f64).collect();
/// assert_eq!(pairwise_sum(&xs), 5050.0);
/// ```
pub fn pairwise_sum(xs: &[f64]) -> f64 {
    const BASE: usize = 128;

    if xs.len() <= BASE {
        return neumaier_sum(xs);
    }

    let mid = xs.len() / 2;
    let left = pairwise_sum(&xs[..mid]);
    let right = pairwise_sum(&xs[mid..]);
    left + right
}

/// Deterministic parallel floating-point reducer.
///
/// `DeterministicReducer` partitions an input slice into fixed-size chunks,
/// computes [`neumaier_sum`] on each chunk, and then combines the per-chunk
/// partial sums via [`pairwise_sum`].  Because the chunk boundaries are fixed
/// and the partial sums are always reduced in index order, results are
/// **bit-identical** regardless of whether one thread or many are used.
///
/// # Design
///
/// - Each chunk is assigned a stable slot in a pre-allocated `partials` array.
/// - Parallel threads write exactly one slot each via a `Mutex<f64>`.
/// - The caller collects slots in deterministic index order before the final
///   `pairwise_sum`, ensuring reproducibility.
///
/// # Examples
///
/// ```
/// # use oxiphysics_core::compensated::DeterministicReducer;
/// let xs: Vec<f64> = (0..10_000).map(|i| i as f64 * 0.001).collect();
/// let reducer = DeterministicReducer::default();
/// let serial   = reducer.sum(&xs);
/// let parallel = reducer.par_sum(&xs, 4);
/// assert_eq!(serial.to_bits(), parallel.to_bits());
/// ```
#[derive(Debug, Clone)]
pub struct DeterministicReducer {
    /// Number of input elements per processing chunk.
    pub chunk: usize,
}

impl Default for DeterministicReducer {
    fn default() -> Self {
        Self { chunk: 4096 }
    }
}

impl DeterministicReducer {
    /// Creates a `DeterministicReducer` with the given chunk size.
    ///
    /// A larger `chunk` reduces overhead but may slightly increase numerical
    /// error in the partial sums (still compensated).  The default of 4096
    /// is a reasonable starting point for most workloads.
    pub fn new(chunk: usize) -> Self {
        Self { chunk }
    }

    /// Computes the deterministic compensated sum of `xs` serially.
    ///
    /// The slice is partitioned into chunks of `self.chunk` elements (the last
    /// chunk may be smaller).  Each chunk is reduced by [`neumaier_sum`], and
    /// the resulting partial sums are combined via [`pairwise_sum`].
    pub fn sum(&self, xs: &[f64]) -> f64 {
        if xs.is_empty() {
            return 0.0;
        }
        let partials: Vec<f64> = xs
            .chunks(self.chunk)
            .map(neumaier_sum)
            .collect();
        pairwise_sum(&partials)
    }

    /// Computes the deterministic compensated sum of `xs` using `n_threads` worker threads.
    ///
    /// The result is **bit-identical** to [`sum`](DeterministicReducer::sum) for the same
    /// input regardless of `n_threads`, because:
    /// 1. Chunk boundaries are fixed by `self.chunk` (not by thread count).
    /// 2. Each chunk writes to a pre-assigned slot in `partials`.
    /// 3. The final reduction reads slots in index order.
    ///
    /// Poisoned mutex guards are recovered via `unwrap_or_else(|e| e.into_inner())`
    /// so that a panicking worker thread does not propagate into the caller.
    pub fn par_sum(&self, xs: &[f64], n_threads: usize) -> f64 {
        if xs.is_empty() {
            return 0.0;
        }

        let n_chunks = xs.len().div_ceil(self.chunk);

        // Share the data read-only across threads.
        let data: Arc<Vec<f64>> = Arc::new(xs.to_vec());

        // Pre-allocate one mutex-protected slot per chunk.
        let partials: Arc<Vec<Mutex<f64>>> = Arc::new(
            (0..n_chunks).map(|_| Mutex::new(0.0_f64)).collect(),
        );

        let pool = WorkStealingPool::new(n_threads);

        for i in 0..n_chunks {
            let start = i * self.chunk;
            let end = ((i + 1) * self.chunk).min(xs.len());

            let data_clone = Arc::clone(&data);
            let partials_clone = Arc::clone(&partials);

            pool.submit(move || {
                let partial = neumaier_sum(&data_clone[start..end]);
                let mut slot = partials_clone[i]
                    .lock()
                    .unwrap_or_else(|e| e.into_inner());
                *slot = partial;
            });
        }

        // Consume the pool: drains the queue on the calling thread, then joins workers.
        pool.join();

        // Collect partials in deterministic index order.
        let collected: Vec<f64> = partials
            .iter()
            .map(|m| *m.lock().unwrap_or_else(|e| e.into_inner()))
            .collect();

        pairwise_sum(&collected)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neumaier_beats_naive() {
        // naive sum: 1e16 + 1.0 - 1e16 = 0.0 due to float cancellation
        let xs = [1e16_f64, 1.0, -1e16];
        assert_eq!(neumaier_sum(&xs), 1.0);
        let naive: f64 = xs.iter().sum();
        assert_eq!(naive, 0.0); // confirm naive fails
    }

    #[test]
    fn pairwise_close_to_neumaier() {
        // Generate values using simple LCG pattern
        let xs: Vec<f64> = (0_u64..1024)
            .scan(6364136223846793005_u64, |state, _| {
                *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                Some((*state as f64) / u64::MAX as f64 - 0.5)
            })
            .collect();
        let n = neumaier_sum(&xs);
        let p = pairwise_sum(&xs);
        assert!((n - p).abs() < 1e-10, "neumaier={n} pairwise={p}");
    }

    #[test]
    fn deterministic_reducer_bit_identical() {
        // Build 1_000_000 values with a simple pattern
        let xs: Vec<f64> = (0_u64..1_000_000)
            .scan(1234567890_u64, |state, _| {
                *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                Some((*state as f64) / u64::MAX as f64 * 200.0 - 100.0)
            })
            .collect();
        let reducer = DeterministicReducer::default();
        let serial = reducer.sum(&xs);
        for &n in &[1_usize, 2, 4, 8] {
            let parallel = reducer.par_sum(&xs, n);
            assert_eq!(
                serial.to_bits(), parallel.to_bits(),
                "par_sum(n_threads={n}) not bit-identical to sum(): serial={serial} parallel={parallel}"
            );
        }
    }

    #[test]
    fn ill_conditioned_sum() {
        // Alternating large cancellations
        let xs: Vec<f64> = (0..1000)
            .flat_map(|i| [1e15_f64 * (i as f64 + 1.0), -1e15 * (i as f64 + 1.0)])
            .collect();
        let result = neumaier_sum(&xs);
        assert!(result.abs() < 1.0, "Expected near-zero, got {result}");
    }
}
