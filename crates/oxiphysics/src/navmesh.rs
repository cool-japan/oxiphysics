// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Navigation mesh with A\* pathfinding and funnel-algorithm path smoothing.
//!
//! This module implements a static triangle navigation mesh that supports
//! efficient A\* graph search over the dual graph of triangles, followed by
//! the Simple Stupid Funnel (SSF) algorithm to produce smooth, string-pulled
//! waypoint sequences through the corridor of triangles returned by A\*.
//!
//! ## Types
//!
//! - `NavMesh` — static triangle mesh with precomputed triangle adjacency.
//! - `Circle2D` — 2D circular obstacle (x-z plane) used to narrow portals.
//! - `Path` — found path with waypoints and total length.
//! - `NavError` — errors from pathfinding (no path, point outside mesh, …).
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::navmesh::{NavMesh, Circle2D};
//!
//! let verts = vec![
//!     [0.0, 0.0, 0.0f64],
//!     [1.0, 0.0, 0.0],
//!     [0.0, 0.0, 1.0],
//!     [1.0, 0.0, 1.0],
//! ];
//! let tris = vec![[0u32, 1, 2], [1, 3, 2]];
//! let mesh = NavMesh::from_triangles(verts, tris);
//! let path = mesh.find_path([0.1, 0.0, 0.1], [0.9, 0.0, 0.9], 0.0, &[]);
//! assert!(path.is_ok());
//! ```

use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

// ---------------------------------------------------------------------------
// Helper math (private)
// ---------------------------------------------------------------------------

#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

/// Signed 2D cross product (positive = b is CCW from a).
#[inline]
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Centroid of a triangle.
#[inline]
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

/// Point-in-triangle test using barycentric coordinates (XZ plane).
fn point_in_tri_xz(p: [f64; 3], verts: &[[f64; 3]], tri: [u32; 3]) -> bool {
    let v0 = verts[tri[0] as usize];
    let v1 = verts[tri[1] as usize];
    let v2 = verts[tri[2] as usize];

    let e0 = sub2([v1[0] - v0[0], v1[2] - v0[2]], [0.0, 0.0]);
    let _ = e0; // suppress; computed differently below
    let d0 = [v1[0] - v0[0], v1[2] - v0[2]];
    let d1 = [v2[0] - v0[0], v2[2] - v0[2]];
    let dp = [p[0] - v0[0], p[2] - v0[2]];

    let d00 = dot2(d0, d0);
    let d01 = dot2(d0, d1);
    let d11 = dot2(d1, d1);
    let d20 = dot2(dp, d0);
    let d21 = dot2(dp, d1);

    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-14 {
        return false; // degenerate triangle
    }
    let inv = 1.0 / denom;
    let v = (d11 * d20 - d01 * d21) * inv;
    let w = (d00 * d21 - d01 * d20) * inv;
    let u = 1.0 - v - w;

    // Inside if all barycentric coords are in [0, 1].
    const EPS: f64 = -1e-8;
    u >= EPS && v >= EPS && w >= EPS
}

/// Midpoint between two 3D points.
#[inline]
fn midpoint3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}

// ---------------------------------------------------------------------------
// Obstacle check
// ---------------------------------------------------------------------------

/// Returns true if the XZ midpoint of edge (a, b) is inside any of the obstacles
/// extended by `agent_radius`.
fn edge_blocked(a: [f64; 3], b: [f64; 3], agent_radius: f64, obstacles: &[Circle2D]) -> bool {
    let mid = midpoint3(a, b);
    let mx = mid[0];
    let mz = mid[2];
    for obs in obstacles {
        let dx = mx - obs.center[0];
        let dz = mz - obs.center[1];
        let dist = (dx * dx + dz * dz).sqrt();
        if dist < obs.radius + agent_radius {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// A* open set entry
// ---------------------------------------------------------------------------

/// Newtype for `f64` cost that implements `Ord` via `total_cmp` and wraps in
/// `Reverse` so the BinaryHeap is a min-heap.
#[derive(Clone, Copy, PartialEq)]
struct Cost(f64);

impl Eq for Cost {}

impl PartialOrd for Cost {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cost {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // We wrap in Reverse in the heap entry, so forward total_cmp here.
        self.0.total_cmp(&other.0)
    }
}

/// Entry stored in the A* open set.  `Reverse` makes BinaryHeap a min-heap.
#[derive(Clone, Copy, Eq, PartialEq)]
struct OpenEntry {
    priority: Reverse<Cost>,
    tri: usize,
}

impl PartialOrd for OpenEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OpenEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.priority.cmp(&other.priority)
    }
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A static triangle navigation mesh.
///
/// Build with [`NavMesh::from_triangles`]; then call [`NavMesh::find_path`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NavMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Each triangle: `[v0, v1, v2]` — indices into `vertices`.
    pub tris: Vec<[u32; 3]>,
    /// For each triangle, `[adj0, adj1, adj2]` — adjacent triangle index
    /// across the edge **opposite** to vertex 0, 1, 2 respectively.
    /// `-1` means a border edge.
    pub adjacency: Vec<[i32; 3]>,
}

/// A 2D circular obstacle in the x-z plane (y is ignored).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Circle2D {
    /// Center `[x, z]`.
    pub center: [f64; 2],
    /// Obstacle radius.
    pub radius: f64,
}

/// A found path.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Path {
    /// Sequence of world-space waypoints, including start and goal.
    pub waypoints: Vec<[f64; 3]>,
    /// Sum of distances between consecutive waypoints.
    pub length: f64,
}

/// Errors that can occur during pathfinding.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum NavError {
    NoPath,
    StartOutsideMesh,
    GoalOutsideMesh,
    EmptyMesh,
}

impl std::fmt::Display for NavError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NavError::NoPath => write!(f, "no path found"),
            NavError::StartOutsideMesh => write!(f, "start point is outside the navigation mesh"),
            NavError::GoalOutsideMesh => write!(f, "goal point is outside the navigation mesh"),
            NavError::EmptyMesh => write!(f, "navigation mesh has no triangles"),
        }
    }
}

impl std::error::Error for NavError {}

// ---------------------------------------------------------------------------
// NavMesh implementation
// ---------------------------------------------------------------------------

impl NavMesh {
    /// Build a `NavMesh` from raw vertex and triangle data.
    ///
    /// Computes the adjacency table automatically.
    pub fn from_triangles(vertices: Vec<[f64; 3]>, tris: Vec<[u32; 3]>) -> Self {
        let n = tris.len();
        let mut adjacency = vec![[-1i32; 3]; n];

        // Map from canonical edge (min_v, max_v) → list of (tri_idx, edge_slot).
        // edge_slot k means the edge opposite vertex k:
        //   slot 0: edge (v1, v2)
        //   slot 1: edge (v0, v2)
        //   slot 2: edge (v0, v1)
        let mut edge_map: HashMap<(u32, u32), Vec<(usize, usize)>> = HashMap::new();

        for (ti, tri) in tris.iter().enumerate() {
            // Edges opposite each vertex:
            let edges = [
                (tri[1], tri[2]), // opposite v0
                (tri[0], tri[2]), // opposite v1
                (tri[0], tri[1]), // opposite v2
            ];
            for (slot, &(va, vb)) in edges.iter().enumerate() {
                let key = if va < vb { (va, vb) } else { (vb, va) };
                edge_map.entry(key).or_default().push((ti, slot));
            }
        }

        // For each edge that has two triangles sharing it, set mutual adjacency.
        for entries in edge_map.values() {
            if entries.len() == 2 {
                let (ti, si) = entries[0];
                let (tj, sj) = entries[1];
                adjacency[ti][si] = tj as i32;
                adjacency[tj][sj] = ti as i32;
            }
            // len == 1 → border edge, remains -1
        }

        Self {
            vertices,
            tris,
            adjacency,
        }
    }

    /// Find which triangle (by index) contains `point` (tested in XZ plane).
    ///
    /// Returns `None` if `point` is outside all triangles.
    pub fn find_triangle(&self, point: [f64; 3]) -> Option<usize> {
        for (i, &tri) in self.tris.iter().enumerate() {
            if point_in_tri_xz(point, &self.vertices, tri) {
                return Some(i);
            }
        }
        None
    }

    /// Find a path from `start` to `goal`.
    ///
    /// Uses A\* over the triangle graph, then the Simple Stupid Funnel
    /// algorithm to smooth the corridor into waypoints.
    ///
    /// `agent_radius` and `obstacles` are used to mark portals as impassable.
    pub fn find_path(
        &self,
        start: [f64; 3],
        goal: [f64; 3],
        agent_radius: f64,
        obstacles: &[Circle2D],
    ) -> Result<Path, NavError> {
        if self.tris.is_empty() {
            return Err(NavError::EmptyMesh);
        }

        let start_tri = self
            .find_triangle(start)
            .ok_or(NavError::StartOutsideMesh)?;
        let goal_tri = self.find_triangle(goal).ok_or(NavError::GoalOutsideMesh)?;

        // Trivial case: same triangle.
        if start_tri == goal_tri {
            let length = len3(sub3(goal, start));
            return Ok(Path {
                waypoints: vec![start, goal],
                length,
            });
        }

        // ---------------------------------------------------------------
        // A* over triangle graph
        // ---------------------------------------------------------------

        // came_from[tri] = (parent_tri, edge_slot_in_parent)
        // edge_slot_in_parent is the slot in `parent_tri` that led to `tri`.
        let mut came_from: HashMap<usize, (usize, usize)> = HashMap::new();
        let mut g_score: HashMap<usize, f64> = HashMap::new();
        let mut open: BinaryHeap<OpenEntry> = BinaryHeap::new();

        g_score.insert(start_tri, 0.0);
        open.push(OpenEntry {
            priority: Reverse(Cost(0.0)),
            tri: start_tri,
        });

        let goal_centroid = centroid(&self.vertices, self.tris[goal_tri]);

        let mut found = false;

        while let Some(entry) = open.pop() {
            let current = entry.tri;

            if current == goal_tri {
                found = true;
                break;
            }

            let current_g = *g_score.get(&current).unwrap_or(&f64::INFINITY);

            let adj = self.adjacency[current];
            for (slot, &neighbor_idx) in adj.iter().enumerate() {
                if neighbor_idx < 0 {
                    continue; // border edge
                }
                let neighbor = neighbor_idx as usize;

                // The edge in `current` at `slot` is opposite vertex `slot`,
                // so it connects vertices (slot+1)%3 and (slot+2)%3.
                let va = self.vertices[self.tris[current][(slot + 1) % 3] as usize];
                let vb = self.vertices[self.tris[current][(slot + 2) % 3] as usize];

                if edge_blocked(va, vb, agent_radius, obstacles) {
                    continue;
                }

                // Cost: centroid-to-centroid distance.
                let c_cur = centroid(&self.vertices, self.tris[current]);
                let c_nbr = centroid(&self.vertices, self.tris[neighbor]);
                let step_cost = len3(sub3(c_nbr, c_cur));
                let tent_g = current_g + step_cost;

                let prev_g = *g_score.get(&neighbor).unwrap_or(&f64::INFINITY);
                if tent_g < prev_g {
                    g_score.insert(neighbor, tent_g);
                    came_from.insert(neighbor, (current, slot));

                    let h = len3(sub3(goal_centroid, c_nbr));
                    let f = tent_g + h;
                    open.push(OpenEntry {
                        priority: Reverse(Cost(f)),
                        tri: neighbor,
                    });

                    if neighbor == goal_tri {
                        // We could break here; let the pop handle it properly.
                        let _ = found; // will be set on pop
                    }
                }
            }
        }

        if !found && !came_from.contains_key(&goal_tri) && start_tri != goal_tri {
            return Err(NavError::NoPath);
        }

        // ---------------------------------------------------------------
        // Reconstruct corridor (list of triangle indices, start → goal)
        // ---------------------------------------------------------------

        let mut corridor: Vec<usize> = Vec::new();
        {
            let mut cur = goal_tri;
            loop {
                corridor.push(cur);
                if cur == start_tri {
                    break;
                }
                match came_from.get(&cur) {
                    Some(&(parent, _)) => cur = parent,
                    None => return Err(NavError::NoPath),
                }
            }
        }
        corridor.reverse(); // now start_tri … goal_tri

        // ---------------------------------------------------------------
        // Build portal list from corridor
        // ---------------------------------------------------------------
        // A portal is the shared edge between corridor[i] and corridor[i+1].
        // For triangle corridor[i] with adjacency slot `edge_slot_in_parent`
        // leading to corridor[i+1], the portal vertices are:
        //   left  = tris[corridor[i]][(edge_slot_in_parent + 1) % 3]
        //   right = tris[corridor[i]][(edge_slot_in_parent + 2) % 3]
        //
        // "Left" and "right" here correspond to a CCW-wound triangle viewed
        // from above (+Y).  The funnel code uses signed cross products (XZ)
        // to determine which side of the funnel a vertex falls on.

        // portals: (left_vertex_3d, right_vertex_3d)
        let mut portals: Vec<([f64; 3], [f64; 3])> = Vec::new();

        for i in 0..corridor.len().saturating_sub(1) {
            let from_tri = corridor[i];
            let to_tri = corridor[i + 1];

            // Find which slot in from_tri points to to_tri.
            let adj = self.adjacency[from_tri];
            let found_slot = adj.iter().position(|&a| a == to_tri as i32);
            let slot = match found_slot {
                Some(s) => s,
                None => return Err(NavError::NoPath),
            };

            let lv = self.tris[from_tri][(slot + 1) % 3] as usize;
            let rv = self.tris[from_tri][(slot + 2) % 3] as usize;
            portals.push((self.vertices[lv], self.vertices[rv]));
        }
        // Append goal as a zero-width portal.
        portals.push((goal, goal));

        // ---------------------------------------------------------------
        // Simple Stupid Funnel (SSF) algorithm
        // ---------------------------------------------------------------
        let waypoints = funnel_smooth(start, &portals);

        let length = compute_path_length(&waypoints);
        Ok(Path { waypoints, length })
    }
}

// ---------------------------------------------------------------------------
// Funnel algorithm
// ---------------------------------------------------------------------------

/// 2D XZ representation of a point.
#[inline]
fn xz(p: [f64; 3]) -> [f64; 2] {
    [p[0], p[2]]
}

/// Reconstruct 3D point, preserving Y from original 3D point.
#[inline]
fn xz_to_3d(xz_pt: [f64; 2], y: f64) -> [f64; 3] {
    [xz_pt[0], y, xz_pt[1]]
}

/// Run the Simple Stupid Funnel algorithm.
///
/// `start` is the apex start point.
/// `portals` is a list of `(left, right)` 3D portal vertices ending with `(goal, goal)`.
///
/// Returns a smoothed list of waypoints `[start, apex1, ..., goal]`.
fn funnel_smooth(start: [f64; 3], portals: &[([f64; 3], [f64; 3])]) -> Vec<[f64; 3]> {
    if portals.is_empty() {
        return vec![start];
    }

    let mut waypoints: Vec<[f64; 3]> = vec![start];

    let apex = xz(start);
    let apex_3d = start;
    let mut left_leg = apex;
    let mut right_leg = apex;

    let mut portal_idx = 0usize;

    while portal_idx < portals.len() {
        let (pl, pr) = portals[portal_idx];
        let new_left = xz(pl);
        let new_right = xz(pr);

        let apex_to_left = sub2(left_leg, apex);
        let apex_to_right = sub2(right_leg, apex);
        let apex_to_new_left = sub2(new_left, apex);
        let apex_to_new_right = sub2(new_right, apex);

        // --- Update right leg ---
        // New right narrows funnel (tightens right) if it's to the LEFT of current right.
        // In XZ: cross2(apex_to_right, apex_to_new_right) > 0 means new_right is CCW from right.
        if cross2(apex_to_right, apex_to_new_right) >= 0.0 {
            // Cross right over left: apex has moved to old right leg.
            if cross2(apex_to_left, apex_to_new_right) > 0.0 {
                // New right is past the left leg → emit left leg as waypoint and restart.
                waypoints.push(apex_3d);
                // Advance apex to left leg position, find its portal.
                // Find the portal where left_leg is a vertex.
                let new_apex_3d = find_portal_vertex_3d(portals, left_leg, portal_idx);
                // Reset funnel to new apex, restart from portal_idx = the one after apex.
                let restart_idx = find_apex_portal_idx(portals, left_leg, portal_idx);
                let funnel_tail = funnel_smooth(new_apex_3d, &portals[restart_idx..]);
                // Append tail (skip first element which would be new_apex_3d duplicated).
                for w in funnel_tail.into_iter().skip(1) {
                    waypoints.push(w);
                }
                return waypoints;
            }
            // Tighten right.
            right_leg = new_right;
        }

        // --- Update left leg ---
        // New left narrows funnel if it's to the RIGHT of current left.
        // cross2(apex_to_left, apex_to_new_left) < 0 means new_left is CW from left.
        if cross2(apex_to_left, apex_to_new_left) <= 0.0 {
            // New left crossed over right leg → emit right leg as waypoint and restart.
            if cross2(apex_to_right, apex_to_new_left) < 0.0 {
                waypoints.push(apex_3d);
                let new_apex_3d = find_portal_vertex_3d(portals, right_leg, portal_idx);
                let restart_idx = find_apex_portal_idx(portals, right_leg, portal_idx);
                let funnel_tail = funnel_smooth(new_apex_3d, &portals[restart_idx..]);
                for w in funnel_tail.into_iter().skip(1) {
                    waypoints.push(w);
                }
                return waypoints;
            }
            // Tighten left.
            left_leg = new_left;
        }

        portal_idx += 1;
    }

    // Append the last waypoint (goal), which is portals.last().0 == portals.last().1.
    if let Some(&(goal, _)) = portals.last() {
        waypoints.push(goal);
    }

    waypoints
}

/// Find the 3D position of a portal vertex that matches `xz_pt` within `portals[0..limit]`.
fn find_portal_vertex_3d(
    portals: &[([f64; 3], [f64; 3])],
    xz_pt: [f64; 2],
    limit: usize,
) -> [f64; 3] {
    for &(pl, pr) in portals.iter().take(limit + 1) {
        if approx_eq2(xz(pl), xz_pt) {
            return pl;
        }
        if approx_eq2(xz(pr), xz_pt) {
            return pr;
        }
    }
    // Fallback: reconstruct from xz with Y=0.
    xz_to_3d(xz_pt, 0.0)
}

/// Find the portal index just after the portal that contains `xz_pt` as a vertex.
fn find_apex_portal_idx(portals: &[([f64; 3], [f64; 3])], xz_pt: [f64; 2], limit: usize) -> usize {
    for (i, &(pl, pr)) in portals.iter().take(limit + 1).enumerate() {
        if approx_eq2(xz(pl), xz_pt) || approx_eq2(xz(pr), xz_pt) {
            return i + 1;
        }
    }
    limit + 1
}

#[inline]
fn approx_eq2(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9
}

// ---------------------------------------------------------------------------
// Path length
// ---------------------------------------------------------------------------

fn compute_path_length(waypoints: &[[f64; 3]]) -> f64 {
    if waypoints.len() < 2 {
        return 0.0;
    }
    waypoints.windows(2).map(|w| len3(sub3(w[1], w[0]))).sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a flat grid nav mesh.
    ///
    /// Produces `rows * cols` quads → `2 * rows * cols` triangles.
    /// Vertex `(i, j)` is at `(i*cell_size, 0, j*cell_size)` for i in
    /// `0..=cols`, j in `0..=rows`.
    /// Each quad is split into two CCW triangles (viewed from +Y):
    ///   tri A: (i,j), (i+1,j), (i,j+1)
    ///   tri B: (i+1,j), (i+1,j+1), (i,j+1)
    fn make_grid_mesh(rows: usize, cols: usize, cell_size: f64) -> NavMesh {
        let mut verts: Vec<[f64; 3]> = Vec::new();
        // vertex index: i*(rows+1) + j
        for i in 0..=cols {
            for j in 0..=rows {
                verts.push([i as f64 * cell_size, 0.0, j as f64 * cell_size]);
            }
        }

        let vi = |i: usize, j: usize| -> u32 { (i * (rows + 1) + j) as u32 };

        let mut tris: Vec<[u32; 3]> = Vec::new();
        for i in 0..cols {
            for j in 0..rows {
                // Lower-left triangle (CCW from above): (i,j), (i+1,j), (i,j+1)
                tris.push([vi(i, j), vi(i + 1, j), vi(i, j + 1)]);
                // Upper-right triangle (CCW from above): (i+1,j), (i+1,j+1), (i,j+1)
                tris.push([vi(i + 1, j), vi(i + 1, j + 1), vi(i, j + 1)]);
            }
        }

        NavMesh::from_triangles(verts, tris)
    }

    // -----------------------------------------------------------------------
    // 1. Adjacency counts for a 2×2 grid (8 triangles)
    // -----------------------------------------------------------------------
    #[test]
    fn test_adjacency_counts() {
        let mesh = make_grid_mesh(2, 2, 1.0);
        assert_eq!(mesh.tris.len(), 8, "2x2 grid should have 8 triangles");

        let border_count: usize = mesh
            .adjacency
            .iter()
            .flat_map(|adj| adj.iter())
            .filter(|&&a| a == -1)
            .count();
        let interior_count: usize = mesh
            .adjacency
            .iter()
            .flat_map(|adj| adj.iter())
            .filter(|&&a| a >= 0)
            .count();

        // Total edge slots = 8 * 3 = 24.
        // A 2x2 grid has 4 quads = 8 triangles.
        // Perimeter edges = 2*(cols + rows) = 2*(2+2) = 8 → 8 border slots.
        // But each internal edge is shared → each contributes 2 interior slots.
        // Internal edges: total edges - perimeter.
        // Total triangle edges (unshared count): 8*3/2 with adjustment for border.
        // Easier: border_count + interior_count = 24.
        assert_eq!(border_count + interior_count, 24, "total edge slots = 24");
        assert!(border_count > 0, "grid must have border edges");
        assert!(interior_count > 0, "grid must have interior edges");

        // In a 2x2 grid: 8 border slots (perimeter has 8 half-edges from triangles).
        // Let's just verify the constraint: at least some interior adjacency exists.
        assert!(
            interior_count >= 2,
            "should have at least some interior adjacency"
        );
    }

    // -----------------------------------------------------------------------
    // 2. find_triangle: point inside a 1×1 mesh
    // -----------------------------------------------------------------------
    #[test]
    fn test_find_triangle_inside() {
        let mesh = make_grid_mesh(1, 1, 1.0);
        let result = mesh.find_triangle([0.5, 0.0, 0.5]);
        assert!(result.is_some(), "point [0.5,0,0.5] should be inside mesh");
    }

    #[test]
    fn test_find_triangle_outside() {
        let mesh = make_grid_mesh(1, 1, 1.0);
        let result = mesh.find_triangle([2.0, 0.0, 2.0]);
        assert!(result.is_none(), "point [2,0,2] should be outside mesh");
    }

    // -----------------------------------------------------------------------
    // 3. Simple path: 10×10 grid
    // -----------------------------------------------------------------------
    #[test]
    fn test_simple_path() {
        let mesh = make_grid_mesh(10, 10, 1.0);
        let path = mesh
            .find_path([0.5, 0.0, 0.5], [9.5, 0.0, 9.5], 0.0, &[])
            .expect("path should exist");
        assert!(
            path.waypoints.len() >= 2,
            "path should have at least start and goal"
        );
        assert!(path.length > 0.0, "path length should be positive");
    }

    // -----------------------------------------------------------------------
    // 4. Obstacle forcing detour (or NoPath if truly blocked)
    // -----------------------------------------------------------------------
    #[test]
    fn test_obstacle_forces_detour_or_nopath() {
        let mesh = make_grid_mesh(10, 10, 1.0);
        // Block the direct diagonal corridor with obstacles.
        let obstacles: Vec<Circle2D> = (1..9)
            .map(|k| Circle2D {
                center: [k as f64 + 0.5, k as f64 + 0.5],
                radius: 1.5,
            })
            .collect();
        let result = mesh.find_path([0.5, 0.0, 0.5], [9.5, 0.0, 9.5], 0.0, &obstacles);
        // Either a detour is found or no path (if fully blocked); both are valid.
        match result {
            Ok(path) => {
                assert!(path.waypoints.len() >= 2);
                assert!(path.length > 0.0);
            }
            Err(NavError::NoPath) => {
                // Acceptable: obstacles fully blocked all routes.
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    // -----------------------------------------------------------------------
    // 5. Obstacle blocks narrow corridor → NoPath
    // -----------------------------------------------------------------------
    #[test]
    fn test_obstacle_blocks_corridor() {
        // 1×5 corridor mesh.
        let mesh = make_grid_mesh(1, 5, 1.0);
        // Place a large obstacle in the middle of the corridor.
        let obstacles = vec![Circle2D {
            center: [2.5, 0.5],
            radius: 2.0, // larger than cell_size
        }];
        let result = mesh.find_path([0.5, 0.0, 0.5], [4.5, 0.0, 0.5], 0.0, &obstacles);
        assert!(
            matches!(result, Err(NavError::NoPath)),
            "large obstacle should block all portals → NoPath"
        );
    }

    // -----------------------------------------------------------------------
    // 6. Funnel smoothing: straight corridor has fewer waypoints than triangles
    // -----------------------------------------------------------------------
    #[test]
    fn test_funnel_smoothing_straight_corridor() {
        // A 1-row, 10-column mesh is a straight horizontal corridor.
        let mesh = make_grid_mesh(1, 10, 1.0);
        let path = mesh
            .find_path([0.5, 0.0, 0.5], [9.5, 0.0, 0.5], 0.0, &[])
            .expect("straight corridor path should exist");
        let corridor_tri_count = 20; // 1x10 → 20 triangles
        assert!(
            path.waypoints.len() < corridor_tri_count,
            "funnel should reduce waypoints; got {} for {} triangles",
            path.waypoints.len(),
            corridor_tri_count
        );
    }

    // -----------------------------------------------------------------------
    // 7. Serde round-trip
    // -----------------------------------------------------------------------
    #[test]
    fn test_serde_round_trip() {
        let mesh = make_grid_mesh(3, 3, 1.0);
        let json = serde_json::to_string(&mesh).expect("serialize mesh");
        let mesh2: NavMesh = serde_json::from_str(&json).expect("deserialize mesh");
        assert_eq!(
            mesh.vertices.len(),
            mesh2.vertices.len(),
            "vertex count should survive round-trip"
        );
        assert_eq!(
            mesh.tris.len(),
            mesh2.tris.len(),
            "tri count should survive round-trip"
        );

        let path = mesh
            .find_path([0.5, 0.0, 0.5], [2.5, 0.0, 2.5], 0.0, &[])
            .expect("path in 3x3 mesh");
        let pjson = serde_json::to_string(&path).expect("serialize path");
        let path2: Path = serde_json::from_str(&pjson).expect("deserialize path");
        assert_eq!(
            path.waypoints.len(),
            path2.waypoints.len(),
            "waypoint count should survive round-trip"
        );
    }

    // -----------------------------------------------------------------------
    // 8. Start outside mesh → StartOutsideMesh error
    // -----------------------------------------------------------------------
    #[test]
    fn test_start_outside_mesh() {
        let mesh = make_grid_mesh(10, 10, 1.0);
        let result = mesh.find_path([100.0, 0.0, 100.0], [5.0, 0.0, 5.0], 0.0, &[]);
        assert!(
            matches!(result, Err(NavError::StartOutsideMesh)),
            "point far outside should give StartOutsideMesh"
        );
    }

    // -----------------------------------------------------------------------
    // Extra: empty mesh → EmptyMesh error
    // -----------------------------------------------------------------------
    #[test]
    fn test_empty_mesh_error() {
        let mesh = NavMesh {
            vertices: vec![],
            tris: vec![],
            adjacency: vec![],
        };
        let result = mesh.find_path([0.0, 0.0, 0.0], [1.0, 0.0, 1.0], 0.0, &[]);
        assert!(matches!(result, Err(NavError::EmptyMesh)));
    }

    // -----------------------------------------------------------------------
    // Extra: same-triangle path is trivial
    // -----------------------------------------------------------------------
    #[test]
    fn test_same_triangle_path() {
        let mesh = make_grid_mesh(1, 1, 1.0);
        // Both points in the same bottom-left triangle.
        let path = mesh
            .find_path([0.1, 0.0, 0.1], [0.2, 0.0, 0.2], 0.0, &[])
            .expect("trivial path");
        assert_eq!(path.waypoints.len(), 2, "trivial path = just start + goal");
    }

    // -----------------------------------------------------------------------
    // Extra: NavError Display
    // -----------------------------------------------------------------------
    #[test]
    fn test_nav_error_display() {
        let e = NavError::NoPath;
        assert!(!format!("{e}").is_empty());
    }
}
