// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the oxiphysics navmesh module.
//!
//! Exposes A*-based pathfinding over a triangle navigation mesh.

#![allow(missing_docs)]
#![allow(dead_code)]

use oxiphysics::navmesh::{Circle2D, NavMesh};
use pyo3::prelude::*;

// ─────────────────────────────────────────────────────────────────────────────
// PyNavMesh
// ─────────────────────────────────────────────────────────────────────────────

/// A triangle navigation mesh for pathfinding.
///
/// Build from raw vertex and triangle data via the constructor.
/// Use `find_path_json` to compute A*-smoothed paths.
#[pyclass(name = "NavMesh")]
pub struct PyNavMesh {
    inner: NavMesh,
}

#[pymethods]
impl PyNavMesh {
    /// Build a `NavMesh` from vertices and triangles.
    ///
    /// `vertices` — list of `[x, y, z]` positions.
    /// `tris` — list of `[v0, v1, v2]` triangle vertex indices (u32).
    /// Adjacency is computed automatically.
    #[new]
    pub fn new(vertices: Vec<[f64; 3]>, tris: Vec<[u32; 3]>) -> Self {
        Self {
            inner: NavMesh::from_triangles(vertices, tris),
        }
    }

    /// Number of vertices in the mesh.
    pub fn vertex_count(&self) -> usize {
        self.inner.vertices.len()
    }

    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.inner.tris.len()
    }

    /// Find which triangle index contains `point` (tested in XZ plane), or `None`.
    pub fn find_triangle(&self, point: [f64; 3]) -> Option<usize> {
        self.inner.find_triangle(point)
    }

    /// Find a path from `start` to `goal` and return it as JSON.
    ///
    /// Returns `{"waypoints": [[x,y,z],...], "length": f64}` on success,
    /// or raises `ValueError` if no path is found or start/goal are outside
    /// the mesh.
    ///
    /// `agent_radius` — minimum clearance from obstacle edges.
    /// `obstacles` — list of `[cx, cz, radius]` circular obstacles (in XZ).
    pub fn find_path_json(
        &self,
        start: [f64; 3],
        goal: [f64; 3],
        agent_radius: f64,
        obstacles: Vec<[f64; 3]>,
    ) -> PyResult<String> {
        let obs: Vec<Circle2D> = obstacles
            .into_iter()
            .map(|o| Circle2D {
                center: [o[0], o[1]],
                radius: o[2],
            })
            .collect();
        self.inner
            .find_path(start, goal, agent_radius, &obs)
            .map(|path| {
                serde_json::to_string(&path).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
            })
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module registration
// ─────────────────────────────────────────────────────────────────────────────

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNavMesh>()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a simple flat 2-triangle mesh in the XZ plane (Y=0).
    ///
    /// ```
    ///  (0,0,0)─────(2,0,0)
    ///     │  ╲   tri1 │
    ///     │  tri0  ╲  │
    ///  (0,0,2)─────(2,0,2)
    /// ```
    fn two_tri_mesh() -> PyNavMesh {
        let verts = vec![
            [0.0f64, 0.0, 0.0], // 0
            [2.0, 0.0, 0.0],    // 1
            [0.0, 0.0, 2.0],    // 2
            [2.0, 0.0, 2.0],    // 3
        ];
        // tri0: 0,1,2   tri1: 1,3,2
        let tris = vec![[0u32, 1, 2], [1, 3, 2]];
        PyNavMesh::new(verts, tris)
    }

    #[test]
    fn test_navmesh_instantiation() {
        let nm = two_tri_mesh();
        assert_eq!(nm.vertex_count(), 4);
        assert_eq!(nm.triangle_count(), 2);
    }

    #[test]
    fn test_navmesh_find_triangle() {
        let nm = two_tri_mesh();
        // Point inside tri0 (centroid ≈ (0.67, 0, 0.67))
        let tri = nm.find_triangle([0.5, 0.0, 0.5]);
        assert!(tri.is_some());
    }

    #[test]
    fn test_navmesh_find_path() {
        let nm = two_tri_mesh();
        // Path from inside tri0 to inside tri1
        let result = nm.find_path_json([0.3, 0.0, 0.3], [1.7, 0.0, 1.7], 0.0, vec![]);
        assert!(result.is_ok(), "path should succeed: {:?}", result);
        let json = result.expect("should be ok");
        assert!(json.contains("waypoints"));
    }

    #[test]
    fn test_navmesh_no_path_outside() {
        let nm = two_tri_mesh();
        let result = nm.find_path_json([10.0, 0.0, 10.0], [0.5, 0.0, 0.5], 0.0, vec![]);
        assert!(result.is_err());
    }
}
