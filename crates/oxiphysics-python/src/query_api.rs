// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics query module.
//!
//! Exposes spatial query API (raycasting, sphere/AABB overlap, k-nearest) to Python.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::query::{QueryFilter, QueryWorld, Ray, RayHit};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn rayhit_to_json(hit: &RayHit) -> String {
    format!(
        "{{\"t\":{},\"body_index\":{},\"normal\":[{},{},{}],\"point\":[{},{},{}]}}",
        hit.t,
        hit.body_index,
        hit.normal[0],
        hit.normal[1],
        hit.normal[2],
        hit.point[0],
        hit.point[1],
        hit.point[2],
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// PyQueryWorld
// ─────────────────────────────────────────────────────────────────────────────

/// Spatial query world — accelerated raycasting and overlap tests.
#[pyclass(name = "QueryWorld")]
pub struct PyQueryWorld {
    inner: QueryWorld,
}

impl Default for PyQueryWorld {
    fn default() -> Self {
        Self {
            inner: QueryWorld::new(),
        }
    }
}

#[pymethods]
impl PyQueryWorld {
    /// Create a new empty query world.
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a sphere body at the given position.
    pub fn add_sphere(&mut self, id: usize, center: [f64; 3], radius: f64, sleeping: bool) {
        self.inner.add_sphere(id, center, radius, sleeping);
    }

    /// Register an AABB body.
    pub fn add_aabb(&mut self, id: usize, min: [f64; 3], max: [f64; 3], sleeping: bool) {
        self.inner.add_aabb(id, min, max, sleeping);
    }

    /// Remove all entries for a body index. Returns how many were removed.
    pub fn remove(&mut self, id: usize) -> usize {
        self.inner.remove(id)
    }

    /// Clear all registered bodies.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Number of registered bodies.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// `true` if no bodies are registered.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Cast a ray and return the closest hit as JSON, or None.
    ///
    /// JSON shape: `{"t": f64, "body_index": usize, "normal": [x,y,z], "point": [x,y,z]}`
    pub fn raycast_closest_json(&self, origin: [f64; 3], direction: [f64; 3]) -> Option<String> {
        let ray = Ray::new(origin, direction);
        let filter = QueryFilter::default();
        self.inner
            .raycast(&ray, &filter)
            .map(|hit| rayhit_to_json(&hit))
    }

    /// Cast a ray and return all hits as a JSON array, sorted nearest-first.
    pub fn raycast_all_json(&self, origin: [f64; 3], direction: [f64; 3]) -> String {
        let ray = Ray::new(origin, direction);
        let filter = QueryFilter::default();
        let hits = self.inner.raycast_all(&ray, &filter);
        let parts: Vec<String> = hits.iter().map(rayhit_to_json).collect();
        format!("[{}]", parts.join(","))
    }

    /// Return all body IDs whose shapes overlap an AABB.
    pub fn overlap_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<usize> {
        let filter = QueryFilter::default();
        self.inner.overlap_aabb(min, max, &filter)
    }

    /// Return all body IDs whose shapes overlap a sphere.
    pub fn overlap_sphere(&self, center: [f64; 3], radius: f64) -> Vec<usize> {
        let filter = QueryFilter::default();
        self.inner.overlap_sphere(center, radius, &filter)
    }

    /// Return the k nearest bodies to `point` as a JSON array of `[id, distance]` pairs.
    pub fn k_nearest_json(&self, point: [f64; 3], k: usize) -> String {
        let filter = QueryFilter::default();
        let results = self.inner.k_nearest(point, k, &filter);
        let parts: Vec<String> = results
            .iter()
            .map(|(id, dist)| format!("[{},{}]", id, dist))
            .collect();
        format!("[{}]", parts.join(","))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyQueryWorld>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_world_instantiation() {
        let mut qw = PyQueryWorld::new();
        assert!(qw.is_empty());
        qw.add_sphere(0, [0.0, 0.0, 0.0], 1.0, false);
        assert_eq!(qw.len(), 1);
        qw.clear();
        assert!(qw.is_empty());
    }

    #[test]
    fn test_raycast_closest() {
        let mut qw = PyQueryWorld::new();
        qw.add_sphere(0, [0.0, 0.0, 5.0], 1.0, false);
        let hit = qw.raycast_closest_json([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert!(hit.is_some());
        assert!(hit.expect("should exist").contains("body_index"));
    }

    #[test]
    fn test_raycast_miss() {
        let mut qw = PyQueryWorld::new();
        qw.add_sphere(0, [0.0, 10.0, 0.0], 1.0, false);
        let hit = qw.raycast_closest_json([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_none());
    }

    #[test]
    fn test_overlap_sphere() {
        let mut qw = PyQueryWorld::new();
        qw.add_sphere(42, [0.0, 0.0, 0.0], 1.0, false);
        let hits = qw.overlap_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(hits.contains(&42));
    }
}
