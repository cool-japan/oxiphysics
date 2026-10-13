// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Incremental quadric-error-metric edge-collapse engine (Garland & Heckbert 1997).
//!
//! One engine run shares its state between all collapses of a simplification
//! pass instead of rebuilding it for every collapse:
//!
//! * the vertex quadrics are computed once and merged on a collapse,
//!   `Q[keep] += Q[removed]`;
//! * a vertex -> faces adjacency lets a collapse rewrite only the faces of the two
//!   end points;
//! * a binary heap keyed by `(cost, v0, v1)` holds the candidate edges.  Entries
//!   carry the version stamps of both end points and are discarded lazily when a
//!   collapse has invalidated them, so after a collapse only the edges incident to
//!   the surviving vertex are re-costed and pushed again.
//!
//! The cost of a pass over a mesh with `F` triangles is therefore
//! `O(F log F)` instead of the `O(F^2)` of recomputing every quadric and edge for
//! each collapse.  Ties are broken on the vertex pair, so the result does not
//! depend on hash-map iteration order.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::functions::{cross, edge_collapse_cost, length};
use super::types::{CollapseRecord, QuadricMatrix, SimpleMesh};

/// Which edges may collapse, and which end point survives a collapse.
pub(super) enum CollapsePolicy {
    /// Every edge may collapse; the end point with the smaller index survives.
    Unrestricted,
    /// Flagged vertices are never removed.  An edge between two flagged vertices
    /// never collapses; when exactly one end point is flagged it survives,
    /// otherwise the end point with the smaller index does.
    KeepFlagged(Vec<bool>),
    /// An edge that belongs to exactly one live triangle (a boundary edge) never
    /// collapses.
    SkipBoundaryEdges,
}

/// Result of one applied collapse.
#[derive(Debug, Clone, Copy)]
pub(super) struct CollapseOutcome {
    /// The surviving vertex.
    pub(super) kept: usize,
    /// The vertex that was merged into `kept`.
    pub(super) removed: usize,
    /// Position of the surviving vertex before the collapse.
    pub(super) old_pos: [f64; 3],
    /// Position of the surviving vertex after the collapse.
    pub(super) new_pos: [f64; 3],
}

impl CollapseOutcome {
    /// The progressive-mesh record of this collapse.
    pub(super) fn record(&self) -> CollapseRecord {
        CollapseRecord {
            v_kept: self.kept,
            v_removed: self.removed,
            old_pos: self.old_pos,
            new_pos: self.new_pos,
        }
    }
}

/// A candidate collapse in the priority queue.
struct HeapEntry {
    /// Quadric error of the collapse.
    cost: f64,
    /// Smaller vertex index of the edge.
    lo: usize,
    /// Larger vertex index of the edge.
    hi: usize,
    /// `true` when `hi` survives and `lo` is removed.
    keep_hi: bool,
    /// Version stamp of `lo` when the entry was created.
    lo_version: u32,
    /// Version stamp of `hi` when the entry was created.
    hi_version: u32,
    /// Optimal position of the surviving vertex.
    pos: [f64; 3],
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for HeapEntry {}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEntry {
    /// Reversed so that `BinaryHeap` pops the smallest `(cost, lo, hi)` first.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.lo.cmp(&self.lo))
            .then_with(|| other.hi.cmp(&self.hi))
    }
}

/// The incremental edge-collapse state of one simplification pass.
pub(super) struct QemEngine {
    /// Working vertex positions.
    positions: Vec<[f64; 3]>,
    /// Accumulated quadric of every vertex.
    quadrics: Vec<QuadricMatrix>,
    /// Working triangle list (indices never move; dead triangles stay in place).
    faces: Vec<[usize; 3]>,
    /// Whether each triangle is still part of the mesh.
    face_alive: Vec<bool>,
    /// Number of triangles that are still alive.
    live_faces: usize,
    /// Live triangles incident to every vertex.
    vertex_faces: Vec<Vec<usize>>,
    /// Whether each vertex has not been merged into another one.
    vertex_alive: Vec<bool>,
    /// Version stamp of each vertex, bumped whenever its quadric or position changes.
    version: Vec<u32>,
    /// Candidate collapses, cheapest first.
    heap: BinaryHeap<HeapEntry>,
    /// Edge admissibility and end point selection.
    policy: CollapsePolicy,
    /// Number of edge-cost evaluations (complexity regression tests).
    #[cfg(test)]
    cost_evaluations: usize,
}

impl QemEngine {
    /// Build the engine for `mesh`.
    ///
    /// Returns `None` when a triangle references a vertex outside the mesh.
    /// Triangles with a repeated vertex index are dropped (they carry no area).
    pub(super) fn new(mesh: &SimpleMesh, policy: CollapsePolicy) -> Option<Self> {
        let n = mesh.vertices.len();
        if mesh.triangles.iter().flatten().any(|&v| v >= n) {
            return None;
        }
        let faces = mesh.triangles.clone();
        let face_alive: Vec<bool> = faces
            .iter()
            .map(|&[a, b, c]| a != b && b != c && a != c)
            .collect();
        let live_faces = face_alive.iter().filter(|&&alive| alive).count();

        let mut quadrics: Vec<QuadricMatrix> = (0..n).map(|_| QuadricMatrix::zero()).collect();
        let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (f, &[a, b, c]) in faces.iter().enumerate() {
            if !face_alive[f] {
                continue;
            }
            for v in [a, b, c] {
                vertex_faces[v].push(f);
            }
            let (va, vb, vc) = (mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]);
            let normal = cross(
                [vb[0] - va[0], vb[1] - va[1], vb[2] - va[2]],
                [vc[0] - va[0], vc[1] - va[1], vc[2] - va[2]],
            );
            let len = length(normal);
            if len < 1e-300 {
                // A sliver without a plane contributes no error term.
                continue;
            }
            let unit = [normal[0] / len, normal[1] / len, normal[2] / len];
            let d = -(unit[0] * va[0] + unit[1] * va[1] + unit[2] * va[2]);
            let plane = QuadricMatrix::from_plane(unit[0], unit[1], unit[2], d);
            for v in [a, b, c] {
                add_assign(&mut quadrics[v], &plane);
            }
        }

        let mut engine = Self {
            positions: mesh.vertices.clone(),
            quadrics,
            faces,
            face_alive,
            live_faces,
            vertex_faces,
            vertex_alive: vec![true; n],
            version: vec![0; n],
            heap: BinaryHeap::new(),
            policy,
            #[cfg(test)]
            cost_evaluations: 0,
        };

        // Seed the heap in a deterministic order.
        let mut edges: Vec<(usize, usize)> = Vec::with_capacity(engine.faces.len() * 3);
        for (f, &[a, b, c]) in engine.faces.iter().enumerate() {
            if engine.face_alive[f] {
                for (p, q) in [(a, b), (b, c), (c, a)] {
                    edges.push((p.min(q), p.max(q)));
                }
            }
        }
        edges.sort_unstable();
        edges.dedup();
        for (lo, hi) in edges {
            engine.push_edge(lo, hi);
        }
        Some(engine)
    }

    /// Number of triangles that are still part of the mesh.
    pub(super) fn live_triangles(&self) -> usize {
        self.live_faces
    }

    /// Number of edge-cost evaluations performed so far.
    #[cfg(test)]
    pub(super) fn cost_evaluations(&self) -> usize {
        self.cost_evaluations
    }

    /// The accumulated quadric of vertex `v`.
    #[cfg(test)]
    pub(super) fn quadric(&self, v: usize) -> &QuadricMatrix {
        &self.quadrics[v]
    }

    /// Whether the edge `(a, b)` belongs to exactly one live triangle.
    fn is_boundary_edge(&self, a: usize, b: usize) -> bool {
        self.vertex_faces[a]
            .iter()
            .filter(|&&f| self.faces[f].contains(&b))
            .count()
            == 1
    }

    /// `(kept, removed)` for the edge `(lo, hi)` under the active policy, or
    /// `None` when the edge must not collapse.
    fn choose_endpoints(&self, lo: usize, hi: usize) -> Option<(usize, usize)> {
        match &self.policy {
            CollapsePolicy::Unrestricted => Some((lo, hi)),
            CollapsePolicy::KeepFlagged(flags) => {
                let flagged = |v: usize| flags.get(v).copied().unwrap_or(false);
                match (flagged(lo), flagged(hi)) {
                    (true, true) => None,
                    (_, false) => Some((lo, hi)),
                    (false, true) => Some((hi, lo)),
                }
            }
            CollapsePolicy::SkipBoundaryEdges => {
                if self.is_boundary_edge(lo, hi) {
                    None
                } else {
                    Some((lo, hi))
                }
            }
        }
    }

    /// Cost the edge `(lo, hi)` and queue it if it may collapse.
    fn push_edge(&mut self, lo: usize, hi: usize) {
        let Some((kept, removed)) = self.choose_endpoints(lo, hi) else {
            return;
        };
        #[cfg(test)]
        {
            self.cost_evaluations += 1;
        }
        let (cost, pos) = edge_collapse_cost(
            self.positions[kept],
            self.positions[removed],
            &self.quadrics[kept],
            &self.quadrics[removed],
        );
        if !cost.is_finite() {
            return;
        }
        self.heap.push(HeapEntry {
            cost,
            lo,
            hi,
            keep_hi: kept == hi,
            lo_version: self.version[lo],
            hi_version: self.version[hi],
            pos,
        });
    }

    /// Whether a queued entry still describes the current mesh.
    fn is_current(&self, entry: &HeapEntry) -> bool {
        self.vertex_alive[entry.lo]
            && self.vertex_alive[entry.hi]
            && self.version[entry.lo] == entry.lo_version
            && self.version[entry.hi] == entry.hi_version
    }

    /// Collapse the cheapest admissible edge whose cost does not exceed `max_cost`.
    ///
    /// Returns `None` when no such edge is left.
    pub(super) fn collapse_cheapest(&mut self, max_cost: f64) -> Option<CollapseOutcome> {
        loop {
            let cost = self.heap.peek()?.cost;
            // `cost <= NaN` is false, so a NaN limit admits nothing.
            let within_limit = cost <= max_cost;
            if !within_limit {
                return None;
            }
            let entry = self.heap.pop()?;
            if !self.is_current(&entry) {
                continue;
            }
            let (kept, removed) = if entry.keep_hi {
                (entry.hi, entry.lo)
            } else {
                (entry.lo, entry.hi)
            };
            let old_pos = self.positions[kept];
            self.apply_collapse(kept, removed, entry.pos);
            return Some(CollapseOutcome {
                kept,
                removed,
                old_pos,
                new_pos: entry.pos,
            });
        }
    }

    /// Merge `removed` into `kept` at `new_pos` and refresh the affected edges.
    fn apply_collapse(&mut self, kept: usize, removed: usize, new_pos: [f64; 3]) {
        let incident = std::mem::take(&mut self.vertex_faces[removed]);
        for f in incident {
            let tri = self.faces[f];
            if tri.contains(&kept) {
                // The triangle contained the collapsed edge: it vanishes.
                self.face_alive[f] = false;
                self.live_faces -= 1;
                for v in tri {
                    if v != removed {
                        self.vertex_faces[v].retain(|&g| g != f);
                    }
                }
            } else {
                self.faces[f] = tri.map(|v| if v == removed { kept } else { v });
                self.vertex_faces[kept].push(f);
            }
        }
        self.vertex_alive[removed] = false;
        self.positions[kept] = new_pos;
        let merged = std::mem::replace(&mut self.quadrics[removed], QuadricMatrix::zero());
        add_assign(&mut self.quadrics[kept], &merged);
        self.version[kept] = self.version[kept].wrapping_add(1);

        // Only the edges incident to the surviving vertex changed.
        let mut neighbours: Vec<usize> = Vec::new();
        for &f in &self.vertex_faces[kept] {
            for v in self.faces[f] {
                if v != kept {
                    neighbours.push(v);
                }
            }
        }
        neighbours.sort_unstable();
        neighbours.dedup();
        for v in neighbours {
            self.push_edge(kept.min(v), kept.max(v));
        }
    }

    /// Write the simplified geometry back into `mesh`: surviving vertices get their
    /// new positions and the triangle list is replaced by the live triangles (in
    /// their original order).  Removed vertices stay in the vertex array unreferenced
    /// so that vertex indices remain stable.
    pub(super) fn write_back(&self, mesh: &mut SimpleMesh) {
        for (v, &p) in self.positions.iter().enumerate() {
            if self.vertex_alive[v] {
                mesh.vertices[v] = p;
            }
        }
        mesh.triangles = self
            .faces
            .iter()
            .zip(self.face_alive.iter())
            .filter_map(|(&tri, &alive)| alive.then_some(tri))
            .collect();
    }
}

/// `target += other`, element-wise.
fn add_assign(target: &mut QuadricMatrix, other: &QuadricMatrix) {
    for (row, other_row) in target.q.iter_mut().zip(other.q.iter()) {
        for (value, &add) in row.iter_mut().zip(other_row.iter()) {
            *value += add;
        }
    }
}
