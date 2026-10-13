// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `#[wasm_bindgen]` wrapper for [`WasmNavMesh`] — exposes the A* navigation
//! mesh to JavaScript with flat-array constructor and query API.

use wasm_bindgen::prelude::*;

use crate::navmesh_bridge::{WasmCircle2D, WasmNavError, WasmNavMesh, make_grid_navmesh};
use crate::wasm_helpers::err_to_jsvalue;

/// JavaScript-accessible wrapper for the A* navigation mesh.
///
/// ## JavaScript example
///
/// ```js
/// // Build a 10×10 grid navmesh, cell size 1 m
/// const nav = WasmNavMeshJs.grid(10, 10, 1.0);
/// // Find path from (0.5, 0, 0.5) to (9.5, 0, 9.5)
/// const path = nav.find_path_flat(0.5, 0.0, 0.5,  9.5, 0.0, 9.5,  0.0);
/// // path is Float64Array of [x0,y0,z0, x1,y1,z1, ...] or empty on failure
/// ```
#[wasm_bindgen]
pub struct WasmNavMeshJs {
    inner: WasmNavMesh,
}

impl WasmNavMeshJs {
    /// Validated, Rust-side form of [`WasmNavMeshJs::from_flat`].
    ///
    /// # Errors
    ///
    /// * [`WasmNavError::MalformedInput`] - `vertices_flat` or `indices_flat`
    ///   is not a whole number of triples;
    /// * [`WasmNavError::InvalidIndex`] - a triangle references a vertex that
    ///   does not exist;
    /// * [`WasmNavError::Degenerate`] - a triangle repeats a vertex index.
    pub fn try_from_flat(
        vertices_flat: &[f64],
        indices_flat: &[u32],
    ) -> Result<WasmNavMeshJs, WasmNavError> {
        if !vertices_flat.len().is_multiple_of(3) {
            return Err(WasmNavError::MalformedInput {
                array: "vertices_flat".to_string(),
                len: vertices_flat.len(),
            });
        }
        if !indices_flat.len().is_multiple_of(3) {
            return Err(WasmNavError::MalformedInput {
                array: "indices_flat".to_string(),
                len: indices_flat.len(),
            });
        }
        let vertices: Vec<[f64; 3]> = vertices_flat
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        let tris: Vec<[u32; 3]> = indices_flat
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        Ok(WasmNavMeshJs {
            inner: WasmNavMesh::from_triangles(vertices, tris)?,
        })
    }
}

#[wasm_bindgen]
impl WasmNavMeshJs {
    /// Build a uniform grid navigation mesh of `rows × cols` quad cells.
    ///
    /// Each cell is subdivided into two triangles. `cell_size` is the edge
    /// length in world units.
    pub fn grid(rows: u32, cols: u32, cell_size: f64) -> WasmNavMeshJs {
        WasmNavMeshJs {
            inner: make_grid_navmesh(rows as usize, cols as usize, cell_size),
        }
    }

    /// Build a navmesh from raw triangle data.
    ///
    /// `vertices_flat` — interleaved `[x0,y0,z0, x1,y1,z1, ...]`.
    /// `indices_flat`  — triangle index triples `[a0,b0,c0, a1,b1,c1, ...]`.
    ///
    /// Throws if the arrays are malformed (`vertices_flat.length % 3 != 0` or
    /// `indices_flat.length % 3 != 0`), if a triangle references a vertex
    /// that does not exist, or if a triangle repeats a vertex index.  Queries
    /// on the returned mesh can therefore never run off the vertex list.
    pub fn from_flat(
        vertices_flat: Vec<f64>,
        indices_flat: Vec<u32>,
    ) -> Result<WasmNavMeshJs, JsValue> {
        WasmNavMeshJs::try_from_flat(&vertices_flat, &indices_flat).map_err(err_to_jsvalue)
    }

    /// Number of triangles in this navmesh.
    pub fn triangle_count(&self) -> u32 {
        self.inner.tris.len() as u32
    }

    /// Number of vertices in this navmesh.
    pub fn vertex_count(&self) -> u32 {
        self.inner.vertices.len() as u32
    }

    /// Find a path from `(sx,sy,sz)` to `(gx,gy,gz)` for an agent of `radius`.
    ///
    /// Returns a flat `Float64Array` of waypoints `[x0,y0,z0, ...]`, or an
    /// empty array if there is no path or start/goal are outside the mesh.
    pub fn find_path_flat(
        &self,
        sx: f64,
        sy: f64,
        sz: f64,
        gx: f64,
        gy: f64,
        gz: f64,
        radius: f64,
    ) -> Vec<f64> {
        match self
            .inner
            .find_path([sx, sy, sz], [gx, gy, gz], radius, &[])
        {
            Ok(path) => path
                .waypoints
                .iter()
                .flat_map(|p| [p[0], p[1], p[2]])
                .collect(),
            Err(_) => vec![],
        }
    }

    /// Find a path avoiding circular obstacles.
    ///
    /// `obstacles_flat` — interleaved `[cx0,cz0,r0, cx1,cz1,r1, ...]` (x and
    /// z are the obstacle horizontal coordinates; y is ignored for 2D obstacle
    /// avoidance).
    pub fn find_path_with_obstacles_flat(
        &self,
        sx: f64,
        sy: f64,
        sz: f64,
        gx: f64,
        gy: f64,
        gz: f64,
        radius: f64,
        obstacles_flat: Vec<f64>,
    ) -> Vec<f64> {
        let obstacles: Vec<WasmCircle2D> = obstacles_flat
            .chunks_exact(3)
            .map(|c| WasmCircle2D {
                center: [c[0], c[1]],
                radius: c[2],
            })
            .collect();
        match self
            .inner
            .find_path([sx, sy, sz], [gx, gy, gz], radius, &obstacles)
        {
            Ok(path) => path
                .waypoints
                .iter()
                .flat_map(|p| [p[0], p[1], p[2]])
                .collect(),
            Err(_) => vec![],
        }
    }

    /// All vertex positions as a flat `Float64Array` `[x0,y0,z0, ...]`.
    pub fn all_vertices_flat(&self) -> Vec<f64> {
        self.inner
            .vertices
            .iter()
            .flat_map(|v| [v[0], v[1], v[2]])
            .collect()
    }

    /// All triangle indices as a flat `Uint32Array` `[a0,b0,c0, a1,b1,c1, ...]`.
    pub fn all_indices_flat(&self) -> Vec<u32> {
        self.inner
            .tris
            .iter()
            .flat_map(|t| [t[0], t[1], t[2]])
            .collect()
    }

    /// Serialise the navmesh as a JSON string.
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Corners of the unit square in the x-z plane.
    const SQUARE: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0];

    #[test]
    fn test_try_from_flat_builds_a_queryable_mesh() {
        let nav = WasmNavMeshJs::try_from_flat(&SQUARE, &[0, 1, 2, 1, 3, 2]).expect("valid mesh");
        assert_eq!(nav.triangle_count(), 2);
        assert_eq!(nav.vertex_count(), 4);
        let path = nav.find_path_flat(0.1, 0.0, 0.1, 0.9, 0.0, 0.9, 0.0);
        assert!(path.len() >= 6, "expected at least two waypoints: {path:?}");
        assert_eq!(path.len() % 3, 0);
    }

    #[test]
    fn test_try_from_flat_rejects_malformed_lengths() {
        let err = WasmNavMeshJs::try_from_flat(&SQUARE[..11], &[0, 1, 2])
            .err()
            .expect("11 vertex floats are not whole points");
        assert_eq!(
            err,
            WasmNavError::MalformedInput {
                array: "vertices_flat".to_string(),
                len: 11
            }
        );
        let err = WasmNavMeshJs::try_from_flat(&SQUARE, &[0, 1])
            .err()
            .expect("2 indices are not a whole triangle");
        assert_eq!(
            err,
            WasmNavError::MalformedInput {
                array: "indices_flat".to_string(),
                len: 2
            }
        );
        assert!(err.to_string().contains("indices_flat"));
    }

    #[test]
    fn test_try_from_flat_rejects_out_of_range_indices() {
        // Index 4 of 4 vertices passed the old length-only check and then
        // aborted the instance during the first path query.
        let err = WasmNavMeshJs::try_from_flat(&SQUARE, &[0, 1, 4])
            .err()
            .expect("vertex 4 does not exist");
        assert_eq!(
            err,
            WasmNavError::InvalidIndex {
                tri: 0,
                index: 4,
                vertex_count: 4
            }
        );
        let err = WasmNavMeshJs::try_from_flat(&SQUARE, &[0, 1, 2, 1, 3, u32::MAX])
            .err()
            .expect("u32::MAX does not exist");
        assert!(matches!(err, WasmNavError::InvalidIndex { tri: 1, .. }));
    }

    #[test]
    fn test_try_from_flat_rejects_degenerate_triangles() {
        let err = WasmNavMeshJs::try_from_flat(&SQUARE, &[0, 1, 1])
            .err()
            .expect("repeated corner");
        assert_eq!(err, WasmNavError::Degenerate { tri: 0 });
    }

    #[test]
    fn test_from_flat_js_ok_path() {
        let nav = WasmNavMeshJs::from_flat(SQUARE.to_vec(), vec![0, 1, 2, 1, 3, 2]);
        assert!(nav.is_ok());
    }
}
