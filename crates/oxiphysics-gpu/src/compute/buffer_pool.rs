// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! A power-of-two bucket pool for reusing CPU-side `Vec<f64>` buffers.
//!
//! Allocations are rounded up to the next power of two so that freed buffers
//! can be reused for similarly-sized requests without reallocating.

use std::collections::HashMap;

/// Pool of reusable `Vec<f64>` allocations, bucketed by capacity (power of two).
///
/// This is used by the CPU-shadow (stub) path to avoid repeated heap allocation
/// when the same buffer sizes are recycled across simulation steps.
pub struct BufferPool {
    /// Map from power-of-two bucket size → list of free buffers.
    pub pools: HashMap<usize, Vec<Vec<f64>>>,
    /// Maximum number of free buffers to retain per bucket before dropping.
    pub cap_per_bucket: usize,
}

impl BufferPool {
    /// Create a new pool with a default per-bucket cap of 64.
    pub fn new() -> Self {
        Self {
            pools: HashMap::new(),
            cap_per_bucket: 64,
        }
    }

    /// Allocate a buffer of at least `len` elements.
    ///
    /// If a free buffer of the correct power-of-two bucket size exists it is
    /// returned directly (contents may be stale); otherwise a fresh
    /// zero-filled allocation is made.
    pub fn alloc(&mut self, len: usize) -> Vec<f64> {
        let bucket = len.next_power_of_two();
        self.pools
            .entry(bucket)
            .or_default()
            .pop()
            .unwrap_or_else(|| vec![0.0; bucket])
    }

    /// Return a buffer to the pool for future reuse.
    ///
    /// If the pool for this bucket is already at capacity, the buffer is
    /// silently dropped (deallocated).
    pub fn free(&mut self, len: usize, buf: Vec<f64>) {
        let bucket = len.next_power_of_two();
        let pool = self.pools.entry(bucket).or_default();
        if pool.len() < self.cap_per_bucket {
            pool.push(buf);
        }
        // Otherwise drop `buf` — no-op, it is deallocated here.
    }

    /// Return the number of free buffers currently held by the pool.
    pub fn total_free(&self) -> usize {
        self.pools.values().map(|v| v.len()).sum()
    }

    /// Clear all pooled buffers, releasing memory back to the allocator.
    pub fn clear(&mut self) {
        self.pools.clear();
    }
}

impl Default for BufferPool {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pool_alloc_returns_correctly_sized_vec() {
        let mut pool = BufferPool::new();
        let buf = pool.alloc(10);
        // bucket = 16 (next power of two ≥ 10)
        assert!(buf.len() >= 10);
        assert_eq!(buf.capacity(), 16);
    }

    #[test]
    fn test_pool_alloc_free_reuses_buffer() {
        let mut pool = BufferPool::new();
        let buf = pool.alloc(8);
        let ptr = buf.as_ptr();
        pool.free(8, buf);
        let buf2 = pool.alloc(8);
        // The reused buffer should be the same allocation.
        assert_eq!(buf2.as_ptr(), ptr);
    }

    #[test]
    fn test_pool_bucket_cap_does_not_grow_unbounded() {
        let mut pool = BufferPool::new();
        pool.cap_per_bucket = 4;
        // Free 8 buffers into the same bucket — only 4 should be retained.
        for _ in 0..8 {
            let buf = vec![0.0f64; 8];
            pool.free(8, buf);
        }
        assert_eq!(pool.total_free(), 4);
    }

    #[test]
    fn test_pool_clear_releases_all() {
        let mut pool = BufferPool::new();
        pool.free(4, vec![0.0; 4]);
        pool.free(8, vec![0.0; 8]);
        assert_eq!(pool.total_free(), 2);
        pool.clear();
        assert_eq!(pool.total_free(), 0);
    }

    #[test]
    fn test_pool_alloc_zero_len_returns_power_of_two() {
        let mut pool = BufferPool::new();
        // 0.next_power_of_two() == 1 in Rust
        let buf = pool.alloc(0);
        assert!(!buf.is_empty());
    }
}
