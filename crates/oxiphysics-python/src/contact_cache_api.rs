// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics contact_cache module.
//!
//! Exposes persistent contact pair cache with warm-start impulse data to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::contact_cache::{ContactCache, ContactPoint};
use pyo3::prelude::*;

/// `(pos_a, pos_b, normal, depth)` contact point input tuple.
type ContactPointInput = ([f64; 3], [f64; 3], [f64; 3], f64);

// ─────────────────────────────────────────────────────────────────────────────
// PyContactCache
// ─────────────────────────────────────────────────────────────────────────────

/// Persistent contact pair cache for warm-starting the constraint solver.
#[pyclass(name = "ContactCache")]
pub struct PyContactCache {
    inner: ContactCache,
}

#[pymethods]
impl PyContactCache {
    /// Create a new contact cache with the given max TTL in steps.
    #[new]
    pub fn new(max_ttl: u32) -> Self {
        Self {
            inner: ContactCache::new(max_ttl),
        }
    }

    /// Begin a new simulation step (increments lifetimes).
    pub fn begin_step(&mut self) {
        self.inner.begin_step();
    }

    /// Register or refresh the contact manifold for body pair (a, b).
    ///
    /// `points` is a list of `(pos_a [x,y,z], pos_b [x,y,z], normal [x,y,z], depth)` tuples.
    pub fn update_pair(&mut self, a: u32, b: u32, points: Vec<ContactPointInput>) {
        let contact_points: Vec<ContactPoint> = points
            .into_iter()
            .map(|(pos_a, pos_b, normal, depth)| ContactPoint {
                pos_a,
                pos_b,
                normal,
                depth,
            })
            .collect();
        self.inner.update_pair(a, b, contact_points);
    }

    /// Evict stale entries. Returns number evicted.
    pub fn evict_stale(&mut self) -> usize {
        self.inner.evict_stale()
    }

    /// Remove the contact entry for a specific pair.
    pub fn remove(&mut self, a: u32, b: u32) {
        self.inner.remove(a, b);
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Number of contact entries.
    pub fn entry_count(&self) -> usize {
        self.inner.entry_count()
    }

    /// `true` if the cache has no entries.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Find the pair with the highest impact (max impulse sum). Returns (a, b, impulse) or None.
    pub fn highest_impact_pair(&self) -> Option<(u32, u32, f64)> {
        self.inner.highest_impact_pair()
    }

    /// Get all pairs with summed normal impulse above `threshold` as a JSON array.
    ///
    /// Each element: `[body_a, body_b, total_impulse]`
    pub fn pairs_above_impulse_threshold_json(&self, threshold: f64) -> String {
        let pairs = self.inner.pairs_above_impulse_threshold(threshold);
        let parts: Vec<String> = pairs
            .iter()
            .map(|(a, b, imp)| format!("[{},{},{}]", a, b, imp))
            .collect();
        format!("[{}]", parts.join(","))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyContactCache>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contact_cache_instantiation() {
        let mut cache = PyContactCache::new(4);
        assert!(cache.is_empty());
        cache.begin_step();
        cache.update_pair(
            0,
            1,
            vec![([0.0, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 1.0, 0.0], 0.01)],
        );
        assert_eq!(cache.entry_count(), 1);
        cache.clear();
        assert!(cache.is_empty());
    }

    #[test]
    fn test_contact_cache_evict_stale() {
        let mut cache = PyContactCache::new(2);
        cache.begin_step();
        cache.update_pair(
            0,
            1,
            vec![([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01)],
        );
        // Age the entry beyond max_lifetime
        cache.begin_step();
        cache.begin_step();
        cache.begin_step();
        let evicted = cache.evict_stale();
        assert!(evicted >= 1);
    }
}
