// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the navigation mesh and A* pathfinding system.
//!
//! Exposes a JSON-oriented surface around `NavMesh`, `Path`, and obstacle
//! types suitable for use across the WASM boundary.

use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, flatten_vec3s, to_js_value};

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Centroid of a triangle.
fn centroid(verts: &[[f64; 3]], tri: [u32; 3]) -> [f64; 3] {
    let v0 = verts[tri[0] as usize];
    let v1 = verts[tri[1] as usize];
    let v2 = verts[tri[2] as usize];
    [
        (v0[0] + v1[0] + v2[0]) / 3.0,
        (v0[1] + v1[1] + v2[1]) / 3.0,
        (v0[2] + v1[2] + v2[2]) / 3.0,
    ]
}

/// Barycentric point-in-triangle test (x-z plane).
fn point_in_tri_xz(p: [f64; 3], verts: &[[f64; 3]], tri: [u32; 3]) -> bool {
    let v0 = verts[tri[0] as usize];
    let v1 = verts[tri[1] as usize];
    let v2 = verts[tri[2] as usize];
    let pxz = [p[0], p[2]];
    let a = [v0[0], v0[2]];
    let b = [v1[0], v1[2]];
    let c = [v2[0], v2[2]];
    let ab = sub2(b, a);
    let bc = sub2(c, b);
    let ca = sub2(a, c);
    let ap = sub2(pxz, a);
    let bp = sub2(pxz, b);
    let cp = sub2(pxz, c);
    let d1 = cross2(ab, ap);
    let d2 = cross2(bc, bp);
    let d3 = cross2(ca, cp);
    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    !(has_neg && has_pos)
}

// ---------------------------------------------------------------------------
// WasmNavError
// ---------------------------------------------------------------------------

/// Error variants returned by pathfinding.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmNavError {
    EmptyMesh,
    StartOutsideMesh,
    GoalOutsideMesh,
    NoPath,
}

impl std::fmt::Display for WasmNavError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WasmNavError::EmptyMesh => write!(f, "navigation mesh is empty"),
            WasmNavError::StartOutsideMesh => write!(f, "start point is outside the mesh"),
            WasmNavError::GoalOutsideMesh => write!(f, "goal point is outside the mesh"),
            WasmNavError::NoPath => write!(f, "no path exists between start and goal"),
        }
    }
}

// ---------------------------------------------------------------------------
// WasmCircle2D — 2D circular obstacle
// ---------------------------------------------------------------------------

/// 2D circular obstacle used to narrow portal widths during path queries.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCircle2D {
    #[wasm_bindgen(skip)]
    pub center: [f64; 2],
    pub radius: f64,
}

#[wasm_bindgen]
impl WasmCircle2D {
    /// Construct a circle at `(cx, cz)` with the given radius (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(cx: f64, cz: f64, radius: f64) -> WasmCircle2D {
        WasmCircle2D {
            center: [cx, cz],
            radius,
        }
    }

    /// Centre of the circle as a flat `Vec<f64>` of length 2.
    #[wasm_bindgen(js_name = "get_center")]
    pub fn get_center_js(&self) -> Vec<f64> {
        self.center.to_vec()
    }

    /// Replace the centre with `(cx, cz)`.
    #[wasm_bindgen(js_name = "set_center")]
    pub fn set_center_js(&mut self, cx: f64, cz: f64) {
        self.center = [cx, cz];
    }

    /// Serialise as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }
}

// ---------------------------------------------------------------------------
// WasmPath
// ---------------------------------------------------------------------------

/// Computed path: list of waypoints and total arc length.
#[wasm_bindgen]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WasmPath {
    #[wasm_bindgen(skip)]
    pub waypoints: Vec<[f64; 3]>,
    pub length: f64,
}

impl WasmPath {
    /// Serialise the path as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

#[wasm_bindgen]
impl WasmPath {
    /// Number of waypoints in the path.
    #[wasm_bindgen(js_name = "waypoint_count")]
    pub fn waypoint_count_js(&self) -> usize {
        self.waypoints.len()
    }

    /// All waypoints flattened to `[x0, y0, z0, x1, y1, z1, ...]`.
    #[wasm_bindgen(js_name = "waypoints_flat")]
    pub fn waypoints_flat_js(&self) -> Vec<f64> {
        flatten_vec3s(&self.waypoints)
    }

    /// Single waypoint at `index` as a flat `Vec<f64>` of length 3, or empty
    /// if `index` is out of range.
    #[wasm_bindgen(js_name = "waypoint_at")]
    pub fn waypoint_at_js(&self, index: usize) -> Vec<f64> {
        match self.waypoints.get(index) {
            Some(w) => w.to_vec(),
            None => Vec::new(),
        }
    }

    /// Serialise the path as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Serialise the path as a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmNavMesh
// ---------------------------------------------------------------------------

/// WASM wrapper for a static triangle navigation mesh with A* pathfinding.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmNavMesh {
    #[wasm_bindgen(skip)]
    pub vertices: Vec<[f64; 3]>,
    #[wasm_bindgen(skip)]
    pub tris: Vec<[u32; 3]>,
    /// Adjacency list: for triangle `i`, `adjacency[i*3 + j]` = adjacent triangle
    /// index across edge j, or -1 if none.
    #[wasm_bindgen(skip)]
    pub adjacency: Vec<i32>,
}

impl WasmNavMesh {
    /// Construct a nav mesh from vertex positions and triangle index triples.
    pub fn from_triangles(vertices: Vec<[f64; 3]>, tris: Vec<[u32; 3]>) -> Self {
        let adjacency = Self::build_adjacency(&vertices, &tris);
        WasmNavMesh {
            vertices,
            tris,
            adjacency,
        }
    }

    fn build_adjacency(vertices: &[[f64; 3]], tris: &[[u32; 3]]) -> Vec<i32> {
        let n = tris.len();
        let mut adj = vec![-1_i32; n * 3];
        // Build edge → triangle map
        let mut edge_map: HashMap<(u32, u32), usize> = HashMap::new();
        for (ti, tri) in tris.iter().enumerate() {
            for e in 0..3 {
                let a = tri[e];
                let b = tri[(e + 1) % 3];
                let key = (a.min(b), a.max(b));
                if let Some(&other_ti) = edge_map.get(&key) {
                    // Find edge indices on both triangles
                    for oe in 0..3 {
                        let oa = tris[other_ti][oe];
                        let ob = tris[other_ti][(oe + 1) % 3];
                        let other_key = (oa.min(ob), oa.max(ob));
                        if other_key == key {
                            adj[ti * 3 + e] = other_ti as i32;
                            adj[other_ti * 3 + oe] = ti as i32;
                        }
                    }
                } else {
                    let _ = vertices; // silence unused warning
                    edge_map.insert(key, ti);
                }
            }
        }
        adj
    }

    /// Find the triangle index containing `p` (x-z plane), or `None`.
    pub fn find_triangle(&self, p: [f64; 3]) -> Option<usize> {
        self.tris
            .iter()
            .enumerate()
            .find(|(_, tri)| point_in_tri_xz(p, &self.vertices, **tri))
            .map(|(i, _)| i)
    }

    /// Find a path from `start` to `goal` avoiding `obstacles`.
    pub fn find_path(
        &self,
        start: [f64; 3],
        goal: [f64; 3],
        agent_radius: f64,
        obstacles: &[WasmCircle2D],
    ) -> Result<WasmPath, WasmNavError> {
        if self.tris.is_empty() {
            return Err(WasmNavError::EmptyMesh);
        }
        let start_tri = self
            .find_triangle(start)
            .ok_or(WasmNavError::StartOutsideMesh)?;
        let goal_tri = self
            .find_triangle(goal)
            .ok_or(WasmNavError::GoalOutsideMesh)?;

        if start_tri == goal_tri {
            let length = len3(sub3(goal, start));
            return Ok(WasmPath {
                waypoints: vec![start, goal],
                length,
            });
        }

        // A* over triangle dual graph
        let n = self.tris.len();
        let mut dist = vec![f64::INFINITY; n];
        let mut prev = vec![usize::MAX; n];
        dist[start_tri] = 0.0;

        // (neg_dist, tri_index)
        let mut heap: BinaryHeap<(Reverse<u64>, usize)> = BinaryHeap::new();
        heap.push((Reverse(0), start_tri));

        let goal_c = centroid(&self.vertices, self.tris[goal_tri]);

        while let Some((_, current)) = heap.pop() {
            if current == goal_tri {
                break;
            }
            let current_c = centroid(&self.vertices, self.tris[current]);
            for e in 0..3 {
                let adj_idx = self.adjacency[current * 3 + e];
                if adj_idx < 0 {
                    continue;
                }
                let neighbor = adj_idx as usize;

                // Check if the shared portal is blocked by any obstacle
                if self.portal_blocked(current, e, agent_radius, obstacles) {
                    continue;
                }

                let neighbor_c = centroid(&self.vertices, self.tris[neighbor]);
                let edge_cost = len3(sub3(neighbor_c, current_c));
                let new_dist = dist[current] + edge_cost;
                if new_dist < dist[neighbor] {
                    dist[neighbor] = new_dist;
                    prev[neighbor] = current;
                    let h = len3(sub3(goal_c, neighbor_c));
                    let priority = ((new_dist + h) * 1_000_000.0) as u64;
                    heap.push((Reverse(priority), neighbor));
                }
            }
        }

        if dist[goal_tri].is_infinite() {
            return Err(WasmNavError::NoPath);
        }

        // Reconstruct triangle corridor
        let mut corridor = Vec::new();
        let mut cur = goal_tri;
        while cur != usize::MAX {
            corridor.push(cur);
            cur = prev[cur];
        }
        corridor.reverse();

        // Funnel smoothing (simplified string-pull)
        let waypoints = self.funnel_smooth(start, goal, &corridor);
        let length: f64 = waypoints.windows(2).map(|w| len3(sub3(w[1], w[0]))).sum();

        Ok(WasmPath { waypoints, length })
    }

    /// Check if the portal at edge `e` of triangle `tri` is too narrow for `agent_radius`.
    fn portal_blocked(
        &self,
        tri: usize,
        edge: usize,
        agent_radius: f64,
        obstacles: &[WasmCircle2D],
    ) -> bool {
        let tri_idx = &self.tris[tri];
        let va = self.vertices[tri_idx[edge] as usize];
        let vb = self.vertices[tri_idx[(edge + 1) % 3] as usize];
        let portal_mid = [
            (va[0] + vb[0]) / 2.0,
            (va[1] + vb[1]) / 2.0,
            (va[2] + vb[2]) / 2.0,
        ];
        let portal_half_len = len3(sub3(vb, va)) / 2.0;
        if portal_half_len < agent_radius {
            return true;
        }
        for obs in obstacles {
            let dx = portal_mid[0] - obs.center[0];
            let dz = portal_mid[2] - obs.center[1];
            let dist = (dx * dx + dz * dz).sqrt();
            if dist < obs.radius + agent_radius {
                // Check if obstacle blocks the full portal width
                let portal_dir = sub3(vb, va);
                let portal_len = len3(portal_dir);
                if portal_len < obs.radius * 2.0 + agent_radius * 2.0 {
                    return true;
                }
            }
        }
        false
    }

    /// Simplified string-pull: take centroids of each triangle in the corridor.
    fn funnel_smooth(&self, start: [f64; 3], goal: [f64; 3], corridor: &[usize]) -> Vec<[f64; 3]> {
        let mut waypoints = vec![start];
        let last_idx = corridor.last().copied().unwrap_or(0);
        // Add centroids of intermediate triangles (not first/last)
        for &tri_idx in corridor.iter().skip(1) {
            if tri_idx != last_idx {
                let c = centroid(&self.vertices, self.tris[tri_idx]);
                // Only add if it's meaningfully different from last waypoint
                if let Some(&last) = waypoints.last()
                    && len3(sub3(c, last)) > 0.01
                {
                    waypoints.push(c);
                }
            }
        }
        waypoints.push(goal);
        waypoints
    }

    /// Serialise the nav mesh as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// WasmNavMesh — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmNavMesh {
    /// Construct a nav mesh from JSON-encoded vertex and triangle arrays.
    ///
    /// `vertices_json` must decode to `Vec<[f64; 3]>` and `tris_json` to
    /// `Vec<[u32; 3]>`.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(vertices_json: &str, tris_json: &str) -> Result<WasmNavMesh, JsValue> {
        let vertices: Vec<[f64; 3]> =
            serde_json::from_str(vertices_json).map_err(err_to_jsvalue)?;
        let tris: Vec<[u32; 3]> = serde_json::from_str(tris_json).map_err(err_to_jsvalue)?;
        Ok(WasmNavMesh::from_triangles(vertices, tris))
    }

    /// Number of vertices in the mesh.
    #[wasm_bindgen(js_name = "vertex_count")]
    pub fn vertex_count_js(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles in the mesh.
    #[wasm_bindgen(js_name = "triangle_count")]
    pub fn triangle_count_js(&self) -> usize {
        self.tris.len()
    }

    /// All vertices flattened to `[x0, y0, z0, x1, y1, z1, ...]`.
    #[wasm_bindgen(js_name = "vertices_flat")]
    pub fn vertices_flat_js(&self) -> Vec<f64> {
        flatten_vec3s(&self.vertices)
    }

    /// All triangles flattened to `[a0, b0, c0, a1, b1, c1, ...]`.
    #[wasm_bindgen(js_name = "tris_flat")]
    pub fn tris_flat_js(&self) -> Vec<u32> {
        let mut out = Vec::with_capacity(self.tris.len() * 3);
        for tri in &self.tris {
            out.extend_from_slice(tri);
        }
        out
    }

    /// Adjacency list (one `i32` per triangle edge, `-1` if none).
    #[wasm_bindgen(js_name = "adjacency_flat")]
    pub fn adjacency_flat_js(&self) -> Vec<i32> {
        self.adjacency.clone()
    }

    /// Find the triangle containing `(x, _, z)`, or `-1` if outside the mesh.
    #[wasm_bindgen(js_name = "find_triangle")]
    pub fn find_triangle_js(&self, x: f64, y: f64, z: f64) -> i32 {
        match self.find_triangle([x, y, z]) {
            Some(i) => i as i32,
            None => -1,
        }
    }

    /// Find a path between `(sx, sy, sz)` and `(gx, gy, gz)` avoiding the
    /// JSON-encoded list of [`WasmCircle2D`] obstacles. Returns the path as
    /// a JSON string, or an error string.
    #[wasm_bindgen(js_name = "find_path_json")]
    pub fn find_path_json_js(
        &self,
        sx: f64,
        sy: f64,
        sz: f64,
        gx: f64,
        gy: f64,
        gz: f64,
        agent_radius: f64,
        obstacles_json: &str,
    ) -> Result<String, JsValue> {
        let obstacles: Vec<WasmCircle2D> = if obstacles_json.is_empty() {
            Vec::new()
        } else {
            serde_json::from_str(obstacles_json).map_err(err_to_jsvalue)?
        };
        let path = self
            .find_path([sx, sy, sz], [gx, gy, gz], agent_radius, &obstacles)
            .map_err(err_to_jsvalue)?;
        serde_json::to_string(&path).map_err(err_to_jsvalue)
    }

    /// Find a path between two world-space points and return it as a `WasmPath` value.
    #[wasm_bindgen(js_name = "find_path")]
    pub fn find_path_native_js(
        &self,
        sx: f64,
        sy: f64,
        sz: f64,
        gx: f64,
        gy: f64,
        gz: f64,
        agent_radius: f64,
        obstacles_json: &str,
    ) -> Result<WasmPath, JsValue> {
        let obstacles: Vec<WasmCircle2D> = if obstacles_json.is_empty() {
            Vec::new()
        } else {
            serde_json::from_str(obstacles_json).map_err(err_to_jsvalue)?
        };
        self.find_path([sx, sy, sz], [gx, gy, gz], agent_radius, &obstacles)
            .map_err(err_to_jsvalue)
    }

    /// Serialise the nav mesh as a JSON string.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> String {
        self.to_json()
    }

    /// Serialise the nav mesh as a `JsValue`.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// Helper: generate a simple grid mesh for testing
// ---------------------------------------------------------------------------

/// Build a uniform grid navmesh of `rows × cols` quads, each with size `cell_size`.
pub fn make_grid_navmesh(rows: usize, cols: usize, cell_size: f64) -> WasmNavMesh {
    let mut vertices = Vec::new();
    for r in 0..=rows {
        for c in 0..=cols {
            vertices.push([c as f64 * cell_size, 0.0, r as f64 * cell_size]);
        }
    }
    let stride = cols + 1;
    let mut tris = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            let tl = (r * stride + c) as u32;
            let tr = tl + 1;
            let bl = tl + stride as u32;
            let br = bl + 1;
            tris.push([tl, tr, bl]);
            tris.push([tr, br, bl]);
        }
    }
    WasmNavMesh::from_triangles(vertices, tris)
}

/// Build a uniform grid navmesh and return it (free function exposed to JS).
#[wasm_bindgen(js_name = "make_grid_navmesh")]
pub fn make_grid_navmesh_js(rows: u32, cols: u32, cell_size: f64) -> WasmNavMesh {
    make_grid_navmesh(rows as usize, cols as usize, cell_size)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_navmesh_bridge_instantiation() {
        let mesh = make_grid_navmesh(2, 2, 1.0);
        let json = mesh.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_navmesh_bridge_empty_mesh_error() {
        let mesh = WasmNavMesh {
            vertices: vec![],
            tris: vec![],
            adjacency: vec![],
        };
        let result = mesh.find_path([0.0; 3], [1.0, 0.0, 1.0], 0.0, &[]);
        assert_eq!(result, Err(WasmNavError::EmptyMesh));
    }

    #[test]
    fn test_navmesh_bridge_start_outside() {
        let mesh = make_grid_navmesh(10, 10, 1.0);
        let result = mesh.find_path([100.0, 0.0, 100.0], [5.0, 0.0, 5.0], 0.0, &[]);
        assert_eq!(result, Err(WasmNavError::StartOutsideMesh));
    }

    #[test]
    fn test_navmesh_bridge_simple_path() {
        let mesh = make_grid_navmesh(10, 10, 1.0);
        let path = mesh
            .find_path([0.5, 0.0, 0.5], [9.5, 0.0, 9.5], 0.0, &[])
            .expect("path should exist");
        assert!(path.waypoints.len() >= 2);
        assert!(path.length > 0.0);
    }

    #[test]
    fn test_navmesh_bridge_same_triangle() {
        let mesh = make_grid_navmesh(1, 1, 1.0);
        let path = mesh
            .find_path([0.1, 0.0, 0.1], [0.2, 0.0, 0.2], 0.0, &[])
            .expect("trivial path");
        assert_eq!(path.waypoints.len(), 2);
    }

    #[test]
    fn test_navmesh_bridge_serde_roundtrip() {
        let mesh = make_grid_navmesh(3, 3, 1.0);
        let json = mesh.to_json();
        let mesh2: WasmNavMesh = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(mesh.vertices.len(), mesh2.vertices.len());
        assert_eq!(mesh.tris.len(), mesh2.tris.len());
    }
}
