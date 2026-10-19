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
    ///
    /// The returned waypoints begin with `start` and end with `goal`.  In
    /// between, the path bends only at mesh vertices -- the corners of the
    /// corridor found by A\* -- and is the shortest path that stays inside that
    /// corridor.  A corner never repeats the waypoint before it.
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

        let corridor = self.find_corridor(start_tri, goal_tri, agent_radius, obstacles)?;
        let mut portals = self.corridor_portals(&corridor)?;

        // A start point that lies on a shared edge (for example on the diagonal
        // of a grid cell) is already past that edge and past every edge before
        // it.  Such portals do not constrain the path, and keeping them would
        // put the funnel apex on the funnel boundary.
        let first_constraining = portals
            .iter()
            .rposition(|&(left, right)| point_on_segment_xz(start, left, right))
            .map_or(0, |index| index + 1);
        let mut funnel_portals = portals.split_off(first_constraining);

        // The goal closes the funnel as a zero-width portal.
        funnel_portals.push((goal, goal));

        let waypoints = funnel_smooth(start, &funnel_portals);
        let length = compute_path_length(&waypoints);
        Ok(Path { waypoints, length })
    }

    /// A\* over the triangle dual graph: the corridor of triangle indices that
    /// leads from `start_tri` to `goal_tri` (both inclusive, in walking order).
    fn find_corridor(
        &self,
        start_tri: usize,
        goal_tri: usize,
        agent_radius: f64,
        obstacles: &[Circle2D],
    ) -> Result<Vec<usize>, NavError> {
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

        Ok(corridor)
    }

    /// The portals crossed by `corridor` (triangle indices of consecutive
    /// adjacent triangles, in walking order).
    ///
    /// Portal `i` is the edge shared by `corridor[i]` and `corridor[i + 1]`.
    /// Its endpoints are returned as `(left, right)` from the viewpoint of a
    /// walker leaving `corridor[i]` through the edge: `left` is the endpoint on
    /// the counter-clockwise side of the walking direction in the x-z plane, which is
    /// the convention the funnel algorithm relies on.  The side is decided
    /// geometrically (from the centroid of `corridor[i]`), so it does not
    /// depend on the winding order of the triangles.
    fn corridor_portals(&self, corridor: &[usize]) -> Result<Vec<Portal>, NavError> {
        let mut portals = Vec::with_capacity(corridor.len().saturating_sub(1));
        for (&from_tri, &to_tri) in corridor.iter().zip(corridor.iter().skip(1)) {
            // Find which slot in from_tri points to to_tri.
            let slot = self.adjacency[from_tri]
                .iter()
                .position(|&adjacent| adjacent == to_tri as i32)
                .ok_or(NavError::NoPath)?;

            // The portal is the edge opposite vertex `slot`.
            let tri = self.tris[from_tri];
            let a = self.vertices[tri[(slot + 1) % 3] as usize];
            let b = self.vertices[tri[(slot + 2) % 3] as usize];

            // `b` is counter-clockwise of `a` as seen from inside the triangle,
            // i.e. on the left of a walker heading for the edge.
            let inside = xz(centroid(&self.vertices, tri));
            let portal = if cross2(sub2(xz(a), inside), sub2(xz(b), inside)) > 0.0 {
                (b, a)
            } else {
                (a, b)
            };
            portals.push(portal);
        }
        Ok(portals)
    }
}

// ---------------------------------------------------------------------------
// Funnel algorithm
// ---------------------------------------------------------------------------

/// A portal between two consecutive corridor triangles: its `(left, right)`
/// endpoints as seen by a walker crossing it (`left` is the endpoint on the
/// counter-clockwise side in the x-z plane).
type Portal = ([f64; 3], [f64; 3]);

/// 2D XZ representation of a point.
#[inline]
fn xz(p: [f64; 3]) -> [f64; 2] {
    [p[0], p[2]]
}

#[inline]
fn approx_eq2(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9
}

#[inline]
fn approx_eq3(a: [f64; 3], b: [f64; 3]) -> bool {
    (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9 && (a[2] - b[2]).abs() < 1e-9
}

/// Append `point` to `waypoints` unless it coincides with the last waypoint, so
/// that a path never contains two consecutive identical points.
fn push_waypoint(waypoints: &mut Vec<[f64; 3]>, point: [f64; 3]) {
    if waypoints
        .last()
        .is_some_and(|&last| approx_eq3(last, point))
    {
        return;
    }
    waypoints.push(point);
}

/// `true` if the XZ projection of `p` lies on the segment `a`-`b` (within a
/// small tolerance).
fn point_on_segment_xz(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> bool {
    let ab = sub2(xz(b), xz(a));
    let ap = sub2(xz(p), xz(a));
    let len_sq = dot2(ab, ab);
    if len_sq < 1e-24 {
        // Degenerate segment: it contains `p` only if the points coincide.
        return dot2(ap, ap) < 1e-18;
    }
    // Distance from the supporting line: |ab x ap| / |ab|.
    if cross2(ab, ap).abs() > 1e-9 * len_sq.sqrt() {
        return false;
    }
    let along = dot2(ab, ap) / len_sq;
    (-1e-9..=1.0 + 1e-9).contains(&along)
}

/// Run the Simple Stupid Funnel (string-pulling) algorithm.
///
/// `start` is the initial apex.  `portals` lists the `(left, right)` portal
/// endpoints in walking order and ends with the zero-width portal
/// `(goal, goal)`; `left` is the endpoint on the counter-clockwise side of the
/// walking direction in the x-z plane (see [`NavMesh::corridor_portals`]).
///
/// Returns the taut path `[start, corner_1, ..., corner_k, goal]`, where the
/// corners are the portal endpoints at which the shortest path through all
/// portals bends.  The result never contains two consecutive coincident points.
///
/// The funnel is described by an apex and a left and a right leg, each with the
/// index of the portal that supplied it.  Every portal narrows the funnel from
/// one or both sides.  When a leg would cross over the opposite one, the
/// opposite leg's endpoint is a corner of the path: it is emitted, becomes the
/// new apex, and scanning restarts with the portal *after* the one that supplied
/// it.  The restart index strictly increases, so the scan terminates.
fn funnel_smooth(start: [f64; 3], portals: &[Portal]) -> Vec<[f64; 3]> {
    let mut waypoints: Vec<[f64; 3]> = vec![start];
    let Some(&(goal, _)) = portals.last() else {
        return waypoints;
    };

    // Portal 0 is the zero-width portal at `start`, so that the apex and both
    // legs can always be identified by a portal index.
    let mut all: Vec<Portal> = Vec::with_capacity(portals.len() + 1);
    all.push((start, start));
    all.extend_from_slice(portals);

    let mut apex = start;
    let mut left = start;
    let mut right = start;
    let mut left_idx = 0_usize;
    let mut right_idx = 0_usize;

    let mut i = 1_usize;
    while let Some(&(new_left, new_right)) = all.get(i) {
        let apex_xz = xz(apex);

        // --- Right leg ---
        // It tightens (moves towards the left leg) when the new right endpoint
        // lies counter-clockwise of it.
        if cross2(sub2(xz(right), apex_xz), sub2(xz(new_right), apex_xz)) >= 0.0 {
            let right_at_apex = approx_eq2(apex_xz, xz(right));
            if right_at_apex || cross2(sub2(xz(left), apex_xz), sub2(xz(new_right), apex_xz)) < 0.0
            {
                right = new_right;
                right_idx = i;
            } else {
                // The right leg crossed over the left leg: the left endpoint is
                // a corner.  It becomes the apex; restart after its portal.
                push_waypoint(&mut waypoints, left);
                apex = left;
                right = left;
                right_idx = left_idx;
                i = left_idx + 1;
                continue;
            }
        }

        // --- Left leg ---
        // It tightens (moves towards the right leg) when the new left endpoint
        // lies clockwise of it.
        if cross2(sub2(xz(left), apex_xz), sub2(xz(new_left), apex_xz)) <= 0.0 {
            let left_at_apex = approx_eq2(apex_xz, xz(left));
            if left_at_apex || cross2(sub2(xz(right), apex_xz), sub2(xz(new_left), apex_xz)) > 0.0 {
                left = new_left;
                left_idx = i;
            } else {
                // The left leg crossed over the right leg: the right endpoint
                // is a corner.  It becomes the apex; restart after its portal.
                push_waypoint(&mut waypoints, right);
                apex = right;
                left = right;
                left_idx = right_idx;
                i = right_idx + 1;
                continue;
            }
        }

        i += 1;
    }

    push_waypoint(&mut waypoints, goal);
    waypoints
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

    // -----------------------------------------------------------------------
    // Funnel regression tests: the corner waypoint used to be dropped and the
    // apex duplicated ([start, start, goal]).
    // -----------------------------------------------------------------------

    /// Mesh of the unit cells `(i, j)` (x = i, z = j) for which `is_open(i, j)`
    /// holds; every cell is split into two triangles.
    ///
    /// The triangles are wound counter-clockwise in the x-z plane (positive 2D
    /// cross product) unless `reversed` is set, which flips every triangle so
    /// that both windings can be exercised.
    fn make_cells_mesh(
        cols: usize,
        rows: usize,
        is_open: impl Fn(usize, usize) -> bool,
        reversed: bool,
    ) -> NavMesh {
        let mut verts: Vec<[f64; 3]> = Vec::new();
        for i in 0..=cols {
            for j in 0..=rows {
                verts.push([i as f64, 0.0, j as f64]);
            }
        }
        let vi = |i: usize, j: usize| -> u32 { (i * (rows + 1) + j) as u32 };
        let mut tris: Vec<[u32; 3]> = Vec::new();
        for i in 0..cols {
            for j in 0..rows {
                if !is_open(i, j) {
                    continue;
                }
                let mut lower = [vi(i, j), vi(i + 1, j), vi(i, j + 1)];
                let mut upper = [vi(i + 1, j), vi(i + 1, j + 1), vi(i, j + 1)];
                if reversed {
                    lower.swap(1, 2);
                    upper.swap(1, 2);
                }
                tris.push(lower);
                tris.push(upper);
            }
        }
        NavMesh::from_triangles(verts, tris)
    }

    /// 5x5 grid of which only the bottom row (`z` in `[0, 1]`) and the right
    /// column (`x` in `[4, 5]`) are walkable: an L-shaped corridor.
    fn make_l_mesh(reversed: bool) -> NavMesh {
        make_cells_mesh(5, 5, |i, j| j == 0 || i == 4, reversed)
    }

    fn assert_no_consecutive_duplicates(waypoints: &[[f64; 3]]) {
        for pair in waypoints.windows(2) {
            assert!(
                len3(sub3(pair[1], pair[0])) > 1e-9,
                "consecutive coincident waypoints in {waypoints:?}"
            );
        }
    }

    #[test]
    fn test_l_corridor_path_bends_at_the_inner_corner() {
        let start = [0.3, 0.0, 0.2];
        let goal = [4.7, 0.0, 4.6];
        let corner = [4.0, 0.0, 1.0];
        for reversed in [false, true] {
            let mesh = make_l_mesh(reversed);
            let path = mesh
                .find_path(start, goal, 0.0, &[])
                .expect("L corridor must be passable");
            assert_eq!(
                path.waypoints,
                vec![start, corner, goal],
                "reversed winding = {reversed}"
            );
            assert_no_consecutive_duplicates(&path.waypoints);
            let expected = len3(sub3(corner, start)) + len3(sub3(goal, corner));
            assert!((path.length - expected).abs() < 1e-12);
            // The straight line would cut through the closed part of the grid.
            assert!(path.length > len3(sub3(goal, start)));
        }
    }

    #[test]
    fn test_mirrored_l_corridor_path_bends_at_the_inner_corner() {
        // Bottom row plus the *left* column: the corridor turns the other way.
        let start = [4.7, 0.0, 0.2];
        let goal = [0.3, 0.0, 4.6];
        let corner = [1.0, 0.0, 1.0];
        for reversed in [false, true] {
            let mesh = make_cells_mesh(5, 5, |i, j| j == 0 || i == 0, reversed);
            let path = mesh
                .find_path(start, goal, 0.0, &[])
                .expect("mirrored L corridor must be passable");
            assert_eq!(
                path.waypoints,
                vec![start, corner, goal],
                "reversed winding = {reversed}"
            );
            assert_no_consecutive_duplicates(&path.waypoints);
        }
    }

    #[test]
    fn test_u_corridor_path_has_two_corners() {
        // Bottom row, top row and right column of a 5x4 grid: a U opening left.
        let start = [0.3, 0.0, 0.2];
        let goal = [0.3, 0.0, 3.7];
        for reversed in [false, true] {
            let mesh = make_cells_mesh(5, 4, |i, j| j == 0 || j == 3 || i == 4, reversed);
            let path = mesh
                .find_path(start, goal, 0.0, &[])
                .expect("U corridor must be passable");
            assert_eq!(
                path.waypoints,
                vec![start, [4.0, 0.0, 1.0], [4.0, 0.0, 3.0], goal],
                "reversed winding = {reversed}"
            );
            assert_no_consecutive_duplicates(&path.waypoints);
        }
    }

    #[test]
    fn test_start_on_a_shared_edge_does_not_duplicate_the_apex() {
        // (0.5, 0.5) and (4.5, 4.5) lie on the diagonal shared by the two
        // triangles of their cells, i.e. exactly on a portal.
        let start = [0.5, 0.0, 0.5];
        let goal = [4.5, 0.0, 4.5];
        for reversed in [false, true] {
            let mesh = make_l_mesh(reversed);
            let path = mesh.find_path(start, goal, 0.0, &[]).expect("path");
            assert_eq!(
                path.waypoints,
                vec![start, [4.0, 0.0, 1.0], goal],
                "reversed winding = {reversed}"
            );
            assert_no_consecutive_duplicates(&path.waypoints);
        }
    }

    #[test]
    fn test_start_on_a_mesh_vertex() {
        let start = [1.0, 0.0, 0.0];
        let goal = [4.6, 0.0, 4.4];
        for reversed in [false, true] {
            let mesh = make_l_mesh(reversed);
            let path = mesh.find_path(start, goal, 0.0, &[]).expect("path");
            assert_eq!(path.waypoints.first(), Some(&start));
            assert_eq!(path.waypoints.last(), Some(&goal));
            assert_no_consecutive_duplicates(&path.waypoints);
            assert_eq!(path.waypoints.len(), 3, "one corner: {:?}", path.waypoints);
            assert_eq!(path.waypoints[1], [4.0, 0.0, 1.0]);
        }
    }

    #[test]
    fn test_straight_corridor_collapses_to_start_and_goal() {
        let start = [0.3, 0.0, 0.2];
        let goal = [9.4, 0.0, 0.8];
        for reversed in [false, true] {
            let mesh = make_cells_mesh(10, 1, |_, _| true, reversed);
            let path = mesh.find_path(start, goal, 0.0, &[]).expect("path");
            assert_eq!(path.waypoints, vec![start, goal]);
            assert!((path.length - len3(sub3(goal, start))).abs() < 1e-12);
        }
    }

    #[test]
    fn test_funnel_without_portals_returns_only_the_start() {
        assert_eq!(funnel_smooth([1.0, 2.0, 3.0], &[]), vec![[1.0, 2.0, 3.0]]);
    }

    #[test]
    fn test_funnel_keeps_the_height_of_the_corner() {
        // Same L corridor, lifted: the corner keeps its own y coordinate.
        let mut mesh = make_l_mesh(false);
        for vertex in &mut mesh.vertices {
            vertex[1] = 0.5 * vertex[0] + 2.0 * vertex[2];
        }
        let start = [0.3, 0.0, 0.2];
        let goal = [4.7, 0.0, 4.6];
        let path = mesh.find_path(start, goal, 0.0, &[]).expect("path");
        assert_eq!(path.waypoints.len(), 3);
        assert_eq!(path.waypoints[1], [4.0, 0.5 * 4.0 + 2.0 * 1.0, 1.0]);
    }

    #[test]
    fn test_portals_are_oriented_independently_of_the_winding() {
        let start = [0.3, 0.0, 0.2];
        let goal = [4.7, 0.0, 4.6];
        for reversed in [false, true] {
            let mesh = make_l_mesh(reversed);
            let start_tri = mesh.find_triangle(start).expect("start in mesh");
            let goal_tri = mesh.find_triangle(goal).expect("goal in mesh");
            let corridor = mesh
                .find_corridor(start_tri, goal_tri, 0.0, &[])
                .expect("corridor");
            let portals = mesh.corridor_portals(&corridor).expect("portals");
            assert_eq!(portals.len(), corridor.len() - 1);
            for (i, &(left, right)) in portals.iter().enumerate() {
                // Walking from corridor[i] through the portal, `right` must be
                // clockwise of `left` as seen from the triangle's centroid.
                let inside = xz(centroid(&mesh.vertices, mesh.tris[corridor[i]]));
                let cross = cross2(sub2(xz(right), inside), sub2(xz(left), inside));
                assert!(
                    cross > 0.0,
                    "portal {i} (reversed = {reversed}): left {left:?} / right {right:?}"
                );
            }
        }
    }

    /// Minimal xorshift64* generator so the property test is deterministic.
    struct XorShift(u64);

    impl XorShift {
        fn next_u64(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }

        /// Uniform in `[0, 1)`.
        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
        }

        /// Uniform in `0..n`.
        fn below(&mut self, n: usize) -> usize {
            (self.unit() * n as f64) as usize % n
        }
    }

    /// Length of the shortest path from `start` through every portal in order,
    /// found by dynamic programming over `samples + 1` points per portal
    /// (portal endpoints included).  It converges to the true optimum from
    /// above as `samples` grows.
    fn dp_shortest_through_portals(start: [f64; 3], portals: &[Portal], samples: usize) -> f64 {
        let mut previous: Vec<([f64; 3], f64)> = vec![(start, 0.0)];
        for &(left, right) in portals {
            let points: Vec<[f64; 3]> = (0..=samples)
                .map(|k| {
                    let t = k as f64 / samples as f64;
                    [
                        left[0] + (right[0] - left[0]) * t,
                        left[1] + (right[1] - left[1]) * t,
                        left[2] + (right[2] - left[2]) * t,
                    ]
                })
                .collect();
            previous = points
                .into_iter()
                .map(|point| {
                    let best = previous
                        .iter()
                        .map(|&(from, cost)| cost + len3(sub3(point, from)))
                        .fold(f64::INFINITY, f64::min);
                    (point, best)
                })
                .collect();
        }
        previous
            .iter()
            .map(|&(_, cost)| cost)
            .fold(f64::INFINITY, f64::min)
    }

    /// Random mazes: whatever A\* returns, the funnel path must start and end
    /// at the query points, stay inside the mesh, contain no coincident
    /// consecutive points, and be the shortest path through the corridor.
    #[test]
    fn test_random_corridors_give_taut_paths_inside_the_mesh() {
        let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
        let mut checked = 0_usize;
        for _ in 0..150 {
            let cols = 3 + rng.below(5);
            let rows = 3 + rng.below(5);
            let mut open = vec![vec![false; rows]; cols];
            for column in &mut open {
                for cell in column.iter_mut() {
                    *cell = rng.unit() < 0.7;
                }
            }
            let reversed = rng.below(2) == 1;
            let mesh = make_cells_mesh(cols, rows, |i, j| open[i][j], reversed);
            if mesh.tris.is_empty() {
                continue;
            }

            let mut random_point = || {
                let tri = mesh.tris[rng.below(mesh.tris.len())];
                let a = mesh.vertices[tri[0] as usize];
                let b = mesh.vertices[tri[1] as usize];
                let c = mesh.vertices[tri[2] as usize];
                let (mut u, mut v) = (rng.unit(), rng.unit());
                if u + v > 1.0 {
                    u = 1.0 - u;
                    v = 1.0 - v;
                }
                [
                    a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]),
                    0.0,
                    a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]),
                ]
            };
            let start = random_point();
            let goal = random_point();

            let Ok(path) = mesh.find_path(start, goal, 0.0, &[]) else {
                continue; // disconnected maze
            };
            checked += 1;

            assert_eq!(path.waypoints.first(), Some(&start));
            assert_eq!(path.waypoints.last(), Some(&goal));
            assert_no_consecutive_duplicates(&path.waypoints);
            assert!((path.length - compute_path_length(&path.waypoints)).abs() < 1e-12);
            assert!(path.length >= len3(sub3(goal, start)) - 1e-9);

            // The polyline never leaves the mesh.
            for pair in path.waypoints.windows(2) {
                for k in 0..=32 {
                    let t = k as f64 / 32.0;
                    let point = add3(pair[0], scale3(sub3(pair[1], pair[0]), t));
                    assert!(
                        mesh.find_triangle(point).is_some(),
                        "path leaves the mesh at {point:?}: {:?}",
                        path.waypoints
                    );
                }
            }

            // ... and it is the shortest path through the A* corridor.
            let (Some(start_tri), Some(goal_tri)) =
                (mesh.find_triangle(start), mesh.find_triangle(goal))
            else {
                panic!("query points are inside the mesh");
            };
            if start_tri == goal_tri {
                assert_eq!(path.waypoints, vec![start, goal]);
                continue;
            }
            let corridor = mesh
                .find_corridor(start_tri, goal_tri, 0.0, &[])
                .expect("find_path succeeded, so a corridor exists");
            let mut portals = mesh.corridor_portals(&corridor).expect("portals");
            portals.push((goal, goal));
            let dp = dp_shortest_through_portals(start, &portals, 128);
            assert!(
                path.length <= dp + 1e-9,
                "funnel path ({}) longer than the DP path ({dp})",
                path.length
            );
            assert!(
                path.length >= dp - 5e-3,
                "funnel path ({}) much shorter than the DP path ({dp}): not through the corridor",
                path.length
            );
        }
        assert!(
            checked > 60,
            "only {checked} connected random mazes checked"
        );
    }
}
