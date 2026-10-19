// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics spatial_grid module.
//!
//! Exposes the uniform spatial hash grid for fast neighbourhood queries.

use oxiphysics::spatial_grid::SpatialGrid;
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PySpatialGrid
// ─────────────────────────────────────────────────────────────────────────────

/// Uniform spatial hash grid for fast neighbourhood and overlap queries.
#[pyclass(name = "SpatialGrid")]
pub struct PySpatialGrid {
    inner: SpatialGrid,
}

#[pymethods]
impl PySpatialGrid {
    /// Create a new spatial grid with the given cell size.
    ///
    /// Raises `ValueError` if `cell_size` is not a positive, finite number.
    #[new]
    pub fn new(cell_size: f64) -> PyResult<Self> {
        let inner = SpatialGrid::try_new(cell_size)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        Ok(Self { inner })
    }

    /// Insert or update a body at the given position.
    pub fn insert(&mut self, id: u32, position: [f64; 3]) {
        self.inner.insert(id, position);
    }

    /// Update an existing body's position.
    pub fn update(&mut self, id: u32, position: [f64; 3]) {
        self.inner.update(id, position);
    }

    /// Remove a body.
    pub fn remove(&mut self, id: u32) {
        self.inner.remove(id);
    }

    /// `true` if the body is registered.
    pub fn contains(&self, id: u32) -> bool {
        self.inner.contains(id)
    }

    /// Get the stored position of a body, or None.
    pub fn position(&self, id: u32) -> Option<[f64; 3]> {
        self.inner.position(id)
    }

    /// Clear all entries.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Number of bodies.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if grid is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Find all body IDs (sorted) within a sphere of the given radius centered at `center`.
    pub fn query_radius(&self, center: [f64; 3], radius: f64) -> Vec<u32> {
        self.inner.query_radius(center, radius)
    }

    /// Find all body IDs (sorted) within an AABB.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<u32> {
        self.inner.query_aabb(min, max)
    }

    /// Find the nearest body ID to a point, or None.
    pub fn nearest(&self, point: [f64; 3]) -> Option<(u32, f64)> {
        self.inner.nearest(point)
    }

    /// Find the k nearest bodies. Returns list of (id, distance).
    pub fn k_nearest(&self, point: [f64; 3], k: usize) -> Vec<(u32, f64)> {
        self.inner.k_nearest(point, k)
    }

    /// Find all pairs `(a, b)` where `a < b` whose positions are within `radius` of each other.
    pub fn pairs_within_radius(&self, radius: f64) -> Vec<(u32, u32)> {
        self.inner.pairs_within_radius(radius)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySpatialGrid>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grid(cell_size: f64) -> PySpatialGrid {
        PySpatialGrid::new(cell_size).expect("valid cell size")
    }

    #[test]
    fn test_spatial_grid_rejects_invalid_cell_sizes() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                PySpatialGrid::new(bad).is_err(),
                "cell size {bad} must be rejected with a ValueError"
            );
        }
        assert!(PySpatialGrid::new(0.25).is_ok());
    }

    #[test]
    fn test_spatial_grid_instantiation() {
        let mut grid = make_grid(1.0);
        assert!(grid.is_empty());
        grid.insert(0, [0.0, 0.0, 0.0]);
        grid.insert(1, [0.5, 0.0, 0.0]);
        assert_eq!(grid.len(), 2);
        grid.remove(0);
        assert_eq!(grid.len(), 1);
    }

    #[test]
    fn test_spatial_grid_radius_query() {
        let mut grid = make_grid(1.0);
        grid.insert(0, [0.0, 0.0, 0.0]);
        grid.insert(1, [100.0, 0.0, 0.0]);
        let results = grid.query_radius([0.0, 0.0, 0.0], 2.0);
        assert!(results.contains(&0));
        assert!(!results.contains(&1));
    }

    #[test]
    fn test_spatial_grid_results_are_sorted_and_deterministic() {
        let mut grid = make_grid(2.0);
        for (id, x) in [(9_u32, 4.0), (3, 0.0), (7, 2.0), (1, 6.0)] {
            grid.insert(id, [x, 0.0, 0.0]);
        }
        // IDs come back sorted, whatever the insertion order.
        assert_eq!(grid.query_radius([2.0, 0.0, 0.0], 2.0), vec![3, 7, 9]);
        assert_eq!(grid.query_aabb([0.0; 3], [4.0, 0.0, 0.0]), vec![3, 7, 9]);
        // Pairs are `(a, b)` with `a < b`, sorted.
        assert_eq!(grid.pairs_within_radius(2.0), vec![(1, 9), (3, 7), (7, 9)]);
        // Nearest neighbours come back nearest first.
        assert_eq!(
            grid.k_nearest([0.0; 3], 3),
            vec![(3, 0.0), (7, 2.0), (9, 4.0)]
        );
    }

    #[test]
    fn test_spatial_grid_nearest() {
        let mut grid = make_grid(1.0);
        grid.insert(42, [1.0, 0.0, 0.0]);
        let nearest = grid.nearest([0.0, 0.0, 0.0]);
        assert!(nearest.is_some());
        let (id, _dist) = nearest.expect("nearest should not be None");
        assert_eq!(id, 42);
    }
}
