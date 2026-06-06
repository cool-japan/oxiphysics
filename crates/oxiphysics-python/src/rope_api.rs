// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics rope module.
//!
//! Exposes rope/chain distance constraints with Verlet integration to Python.

use oxiphysics::rope::Rope;
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyRope
// ─────────────────────────────────────────────────────────────────────────────

/// A rope defined by a chain of links with Verlet integration.
#[pyclass(name = "Rope")]
pub struct PyRope {
    inner: Rope,
}

#[pymethods]
impl PyRope {
    /// Create a straight hanging rope from `root` with the given number of links.
    ///
    /// The rope hangs downward (−Y) from `root`, with `segment_length` per segment.
    #[new]
    pub fn new(
        root: [f64; 3],
        num_links: usize,
        segment_length: f64,
        mass_per_link: f64,
    ) -> PyResult<Self> {
        Rope::new(root, num_links, segment_length, mass_per_link)
            .map(|inner| Self { inner })
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))
    }

    /// Total rope length (sum of segment lengths).
    pub fn total_length(&self) -> f64 {
        self.inner.total_length()
    }

    /// Get the segment endpoints as a list of (start, end) pairs.
    pub fn segments(&self) -> Vec<([f64; 3], [f64; 3])> {
        self.inner.iter_segments().collect()
    }

    /// Number of links in the rope.
    pub fn link_count(&self) -> usize {
        self.inner.links.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRope>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rope_instantiation() {
        let rope = PyRope::new([0.0, 5.0, 0.0], 5, 0.5, 1.0).expect("rope creation failed");
        assert_eq!(rope.link_count(), 5);
        assert!(rope.total_length() > 0.0);
    }

    #[test]
    fn test_rope_segments() {
        let rope = PyRope::new([0.0, 0.0, 0.0], 4, 0.25, 1.0).expect("rope creation failed");
        let segs = rope.segments();
        assert_eq!(segs.len(), 3); // 4 links => 3 segments
    }
}
