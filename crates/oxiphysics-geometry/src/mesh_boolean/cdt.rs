// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constrained Delaunay triangulation of one cut triangle.
//!
//! The domain is always a single parent triangle `T` of a boolean operand,
//! already mapped to 2-D by its dominant-axis projection. The input is a
//! fully *noded* planar straight-line graph:
//!
//! * a counter-clockwise `boundary` cycle made of `T`'s three corners and every
//!   point that lies on one of `T`'s edges, ordered along the edge;
//! * `interior` points strictly inside `T`;
//! * constraint `segments` (intersection-curve pieces) that meet only at shared
//!   endpoints and never run along `T`'s boundary.
//!
//! Construction:
//!
//! 1. the boundary cycle is fanned from an interior `apex` point. This needs no
//!    geometric predicate on the (nearly collinear) boundary points, so points
//!    that floating-point noise puts a hair off `T`'s edges never invert a
//!    triangle;
//! 2. the remaining interior points are inserted incrementally (exact
//!    `orient2d` point location, 1→3 or 2→4 splits) with Lawson legalisation;
//! 3. every constraint is recovered by Sloan's edge-flip algorithm, which
//!    terminates for any valid triangulation when decided with exact
//!    predicates;
//! 4. interior edges that join two points of the same parent edge (they would
//!    bound zero-area slivers) are flipped away, then the unconstrained edges
//!    are made locally Delaunay with the exact `incircle` predicate.
//!
//! Every output triangle is counter-clockwise in the projected frame.

use std::collections::{HashMap, HashSet, VecDeque};

use oxiphysics_core::exact_predicates::{Orientation, incircle, orient2d};

/// Input of [`triangulate`]. All indices refer to `pts`.
pub(super) struct CdtInput<'a> {
    /// Projected point positions.
    pub pts: &'a [[f64; 2]],
    /// Bit `k` is set when the point lies on edge `k` of the parent triangle
    /// (corners carry two bits, interior points none).
    pub masks: &'a [u8],
    /// Counter-clockwise boundary cycle (corners plus on-edge points).
    pub boundary: &'a [usize],
    /// Points strictly inside the parent triangle, including `apex`.
    pub interior: &'a [usize],
    /// Interior point the boundary cycle is fanned from.
    pub apex: usize,
    /// Constraint segments (non-crossing, noded).
    pub segments: &'a [[usize; 2]],
}

/// Triangulate `input`, returning counter-clockwise index triangles, or `None`
/// when the fan from `apex` is not a valid triangulation (the apex is not
/// strictly inside the boundary cycle).
pub(super) fn triangulate(input: &CdtInput<'_>) -> Option<Vec<[usize; 3]>> {
    let mut cdt = Cdt::new(input.pts, input.masks);
    if !cdt.fan(input.apex, input.boundary) {
        return None;
    }
    for &p in input.interior {
        if p != input.apex {
            cdt.insert_point(p);
        }
    }
    for s in input.segments {
        cdt.insert_segment(s[0], s[1], 0);
    }
    cdt.repair_boundary_chords();
    cdt.make_delaunay();
    Some(cdt.triangles())
}

/// Recursion bound for splitting a constraint at collinear vertices.
const MAX_SPLIT_DEPTH: usize = 64;

/// Working triangulation with directed-edge adjacency.
struct Cdt<'a> {
    pts: &'a [[f64; 2]],
    masks: &'a [u8],
    tris: Vec<[usize; 3]>,
    alive: Vec<bool>,
    /// Directed edge `(u, v)` → triangle that traverses it counter-clockwise.
    edges: HashMap<(usize, usize), usize>,
    /// Undirected constrained edges, stored as `(min, max)`.
    constrained: HashSet<(usize, usize)>,
    inserted: Vec<bool>,
}

/// Result of locating a point in the triangulation.
enum Located {
    /// Strictly inside triangle `t`.
    Inside(usize),
    /// On the interior of the directed edge `(u, v)` of triangle `t`.
    OnEdge(usize, usize, usize),
    /// Coincides with an existing vertex.
    OnVertex,
    /// Not inside the triangulated domain.
    Outside,
}

fn key(u: usize, v: usize) -> (usize, usize) {
    if u < v { (u, v) } else { (v, u) }
}

impl<'a> Cdt<'a> {
    fn new(pts: &'a [[f64; 2]], masks: &'a [u8]) -> Self {
        Self {
            pts,
            masks,
            tris: Vec::new(),
            alive: Vec::new(),
            edges: HashMap::new(),
            constrained: HashSet::new(),
            inserted: vec![false; pts.len()],
        }
    }

    fn orient(&self, a: usize, b: usize, c: usize) -> Orientation {
        orient2d(self.pts[a], self.pts[b], self.pts[c])
    }

    fn is_ccw(&self, a: usize, b: usize, c: usize) -> bool {
        self.orient(a, b, c) == Orientation::Positive
    }

    /// Whether `u` and `v` lie on a common edge of the parent triangle.
    fn share_parent_edge(&self, u: usize, v: usize) -> bool {
        let mu = self.masks.get(u).copied().unwrap_or(0);
        let mv = self.masks.get(v).copied().unwrap_or(0);
        mu & mv != 0
    }

    fn add_tri(&mut self, v: [usize; 3]) -> usize {
        let t = self.tris.len();
        self.tris.push(v);
        self.alive.push(true);
        for k in 0..3 {
            self.edges.insert((v[k], v[(k + 1) % 3]), t);
        }
        t
    }

    fn kill_tri(&mut self, t: usize) {
        if !self.alive[t] {
            return;
        }
        self.alive[t] = false;
        let v = self.tris[t];
        for k in 0..3 {
            let e = (v[k], v[(k + 1) % 3]);
            if self.edges.get(&e) == Some(&t) {
                self.edges.remove(&e);
            }
        }
    }

    /// The vertex of triangle `t` that is not `u` or `v`.
    fn third(&self, t: usize, u: usize, v: usize) -> usize {
        let tri = self.tris[t];
        tri.iter()
            .copied()
            .find(|&w| w != u && w != v)
            .unwrap_or(tri[0])
    }

    fn has_edge(&self, u: usize, v: usize) -> bool {
        self.edges.contains_key(&(u, v)) || self.edges.contains_key(&(v, u))
    }

    /// Fan the boundary cycle from `apex`. Returns `false` if a fan triangle is
    /// not strictly counter-clockwise.
    fn fan(&mut self, apex: usize, boundary: &[usize]) -> bool {
        let m = boundary.len();
        if m < 3 {
            return false;
        }
        for i in 0..m {
            let (b0, b1) = (boundary[i], boundary[(i + 1) % m]);
            if !self.is_ccw(apex, b0, b1) {
                return false;
            }
        }
        for i in 0..m {
            let (b0, b1) = (boundary[i], boundary[(i + 1) % m]);
            self.add_tri([apex, b0, b1]);
            self.inserted[b0] = true;
        }
        self.inserted[apex] = true;
        true
    }

    /// Flip the interior edge `(u, v)` if its quadrilateral is strictly convex.
    /// Returns the new diagonal `(x, y)` on success.
    fn flip(&mut self, u: usize, v: usize) -> Option<(usize, usize)> {
        let t1 = *self.edges.get(&(u, v))?;
        let t2 = *self.edges.get(&(v, u))?;
        let x = self.third(t1, u, v);
        let y = self.third(t2, v, u);
        if !(self.is_ccw(u, y, x) && self.is_ccw(y, v, x)) {
            return None;
        }
        let was_constrained = self.constrained.remove(&key(u, v));
        self.kill_tri(t1);
        self.kill_tri(t2);
        self.add_tri([u, y, x]);
        self.add_tri([y, v, x]);
        if was_constrained {
            // Never reached for constrained edges (callers skip them), but keep
            // the constraint set consistent with the edge set regardless.
            self.constrained.insert(key(x, y));
        }
        Some((x, y))
    }

    /// Lawson legalisation after inserting a point: each queued directed edge
    /// has the new point on its left.
    fn legalize(&mut self, mut stack: Vec<(usize, usize)>) {
        let mut budget = 64 * (self.tris.len() + 16);
        while let Some((u, v)) = stack.pop() {
            if budget == 0 {
                break;
            }
            budget -= 1;
            if self.constrained.contains(&key(u, v)) {
                continue;
            }
            let (Some(&t1), Some(&t2)) = (self.edges.get(&(u, v)), self.edges.get(&(v, u))) else {
                continue;
            };
            let x = self.third(t1, u, v);
            let y = self.third(t2, v, u);
            if self.share_parent_edge(x, y) {
                continue;
            }
            if incircle(self.pts[u], self.pts[v], self.pts[x], self.pts[y]) == Orientation::Positive
                && self.flip(u, v).is_some()
            {
                stack.push((u, y));
                stack.push((y, v));
            }
        }
    }

    fn locate(&self, p: usize) -> Located {
        for (t, tri) in self.tris.iter().enumerate() {
            if !self.alive[t] {
                continue;
            }
            let [a, b, c] = *tri;
            let o = [
                self.orient(a, b, p),
                self.orient(b, c, p),
                self.orient(c, a, p),
            ];
            if o.contains(&Orientation::Negative) {
                continue;
            }
            let zeros = o.iter().filter(|&&s| s == Orientation::Degenerate).count();
            match zeros {
                0 => return Located::Inside(t),
                1 => {
                    let (u, v) = if o[0] == Orientation::Degenerate {
                        (a, b)
                    } else if o[1] == Orientation::Degenerate {
                        (b, c)
                    } else {
                        (c, a)
                    };
                    return Located::OnEdge(t, u, v);
                }
                // On two edge lines: p coincides with their shared vertex.
                _ => return Located::OnVertex,
            }
        }
        Located::Outside
    }

    fn insert_point(&mut self, p: usize) {
        if self.inserted[p] {
            return;
        }
        match self.locate(p) {
            Located::Inside(t) => {
                let [a, b, c] = self.tris[t];
                self.kill_tri(t);
                self.add_tri([a, b, p]);
                self.add_tri([b, c, p]);
                self.add_tri([c, a, p]);
                self.inserted[p] = true;
                self.legalize(vec![(a, b), (b, c), (c, a)]);
            }
            Located::OnEdge(t1, a, b) => {
                let c = self.third(t1, a, b);
                let t2 = self.edges.get(&(b, a)).copied();
                let was_constrained = self.constrained.remove(&key(a, b));
                self.kill_tri(t1);
                self.add_tri([a, p, c]);
                self.add_tri([p, b, c]);
                let mut stack = vec![(b, c), (c, a)];
                if let Some(t2) = t2 {
                    let d = self.third(t2, b, a);
                    self.kill_tri(t2);
                    self.add_tri([b, p, d]);
                    self.add_tri([p, a, d]);
                    stack.push((a, d));
                    stack.push((d, b));
                }
                if was_constrained {
                    self.constrained.insert(key(a, p));
                    self.constrained.insert(key(p, b));
                }
                self.inserted[p] = true;
                self.legalize(stack);
            }
            // A coincident vertex cannot occur for noded input; an outside point
            // cannot be represented. Either way the point stays un-inserted and
            // constraints through it are skipped.
            Located::OnVertex | Located::Outside => {}
        }
    }

    /// Whether segment `(a, b)` and edge `(u, v)` cross at a single point
    /// interior to both.
    fn crosses(&self, a: usize, b: usize, u: usize, v: usize) -> bool {
        if u == a || u == b || v == a || v == b {
            return false;
        }
        let o1 = self.orient(a, b, u);
        let o2 = self.orient(a, b, v);
        if o1 == Orientation::Degenerate || o2 == Orientation::Degenerate || o1 == o2 {
            return false;
        }
        let o3 = self.orient(u, v, a);
        let o4 = self.orient(u, v, b);
        o3 != Orientation::Degenerate && o4 != Orientation::Degenerate && o3 != o4
    }

    /// A vertex strictly inside segment `(a, b)`, if any.
    fn vertex_on_segment(&self, a: usize, b: usize) -> Option<usize> {
        let (pa, pb) = (self.pts[a], self.pts[b]);
        let d = [pb[0] - pa[0], pb[1] - pa[1]];
        (0..self.pts.len()).find(|&w| {
            if w == a || w == b || !self.inserted[w] {
                return false;
            }
            if self.orient(a, b, w) != Orientation::Degenerate {
                return false;
            }
            let pw = self.pts[w];
            let s = (pw[0] - pa[0]) * d[0] + (pw[1] - pa[1]) * d[1];
            let e = (pw[0] - pb[0]) * d[0] + (pw[1] - pb[1]) * d[1];
            s > 0.0 && e < 0.0
        })
    }

    fn crossing_edges(&self, a: usize, b: usize) -> Vec<(usize, usize)> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for (t, tri) in self.tris.iter().enumerate() {
            if !self.alive[t] {
                continue;
            }
            for k in 0..3 {
                let (u, v) = (tri[k], tri[(k + 1) % 3]);
                if seen.insert(key(u, v)) && self.crosses(a, b, u, v) {
                    out.push((u, v));
                }
            }
        }
        out
    }

    /// Recover constraint `(a, b)` with Sloan's flip algorithm.
    fn insert_segment(&mut self, a: usize, b: usize, depth: usize) {
        if a == b || a >= self.pts.len() || b >= self.pts.len() {
            return;
        }
        if !self.inserted[a] || !self.inserted[b] {
            return;
        }
        if self.has_edge(a, b) {
            self.constrained.insert(key(a, b));
            return;
        }
        if depth < MAX_SPLIT_DEPTH
            && let Some(w) = self.vertex_on_segment(a, b)
        {
            self.insert_segment(a, w, depth + 1);
            self.insert_segment(w, b, depth + 1);
            return;
        }
        let crossing = self.crossing_edges(a, b);
        // Crossing an existing constraint means two constraints intersect
        // without a shared vertex, which noding rules out; such a segment is
        // left unrecovered rather than breaking the earlier constraint.
        if crossing
            .iter()
            .any(|&(u, v)| self.constrained.contains(&key(u, v)))
        {
            return;
        }
        let mut queue: VecDeque<(usize, usize)> = crossing.into();
        let n = queue.len() + 1;
        let mut budget = 64 * n * n + 1024;
        while let Some((u, v)) = queue.pop_front() {
            if budget == 0 {
                break;
            }
            budget -= 1;
            if !self.edges.contains_key(&(u, v)) || !self.edges.contains_key(&(v, u)) {
                continue;
            }
            match self.flip(u, v) {
                Some((x, y)) => {
                    if self.crosses(a, b, x, y) {
                        queue.push_back((x, y));
                    }
                }
                None => queue.push_back((u, v)),
            }
        }
        if self.has_edge(a, b) {
            self.constrained.insert(key(a, b));
        }
    }

    fn interior_edges(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (&(u, v), _) in self.edges.iter() {
            if u < v && self.edges.contains_key(&(v, u)) {
                out.push((u, v));
            }
        }
        out.sort_unstable();
        out
    }

    /// Flip away unconstrained interior edges whose endpoints lie on a common
    /// parent edge: such an edge runs along the boundary and bounds a sliver.
    fn repair_boundary_chords(&mut self) {
        let cap = 4 * (self.tris.len() + 16);
        for _ in 0..cap {
            let mut changed = false;
            for (u, v) in self.interior_edges() {
                if self.constrained.contains(&key(u, v)) || !self.share_parent_edge(u, v) {
                    continue;
                }
                if self.flip(u, v).is_some() {
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// Lawson flips until every unconstrained interior edge is locally
    /// Delaunay (never creating an edge along the parent boundary).
    fn make_delaunay(&mut self) {
        let mut stack = self.interior_edges();
        let mut budget = 64 * (self.tris.len() + 16) * 4;
        while let Some((u, v)) = stack.pop() {
            if budget == 0 {
                break;
            }
            budget -= 1;
            if self.constrained.contains(&key(u, v)) {
                continue;
            }
            let (Some(&t1), Some(&t2)) = (self.edges.get(&(u, v)), self.edges.get(&(v, u))) else {
                continue;
            };
            let x = self.third(t1, u, v);
            let y = self.third(t2, v, u);
            if self.share_parent_edge(x, y) {
                continue;
            }
            if incircle(self.pts[u], self.pts[v], self.pts[x], self.pts[y]) == Orientation::Positive
                && self.flip(u, v).is_some()
            {
                stack.push((u, y));
                stack.push((y, v));
                stack.push((v, x));
                stack.push((x, u));
            }
        }
    }

    fn triangles(&self) -> Vec<[usize; 3]> {
        self.tris
            .iter()
            .zip(&self.alive)
            .filter(|(_, alive)| **alive)
            .map(|(t, _)| *t)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area2(pts: &[[f64; 2]], t: [usize; 3]) -> f64 {
        let [a, b, c] = t.map(|i| pts[i]);
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    }

    fn has_edge(tris: &[[usize; 3]], u: usize, v: usize) -> bool {
        tris.iter().any(|t| {
            (0..3).any(|k| {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                (a == u && b == v) || (a == v && b == u)
            })
        })
    }

    #[test]
    fn cdt_recovers_constraint_and_tiles_triangle() {
        // Parent triangle (0,0),(4,0),(0,4); one point on each of edges 0 and 2,
        // two interior points, and a chord plus an interior constraint.
        let pts = vec![
            [0.0, 0.0], // 0 corner
            [4.0, 0.0], // 1 corner
            [0.0, 4.0], // 2 corner
            [2.0, 0.0], // 3 on edge 0
            [0.0, 1.0], // 4 on edge 2
            [1.0, 1.0], // 5 interior (apex)
            [0.6, 2.2], // 6 interior
        ];
        let masks = [0b101, 0b011, 0b110, 0b001, 0b100, 0, 0];
        let boundary = [0, 3, 1, 2, 4];
        let segments = [[3, 4], [5, 6]];
        let out = triangulate(&CdtInput {
            pts: &pts,
            masks: &masks,
            boundary: &boundary,
            interior: &[5, 6],
            apex: 5,
            segments: &segments,
        })
        .expect("valid fan");
        let total: f64 = out.iter().map(|&t| area2(&pts, t)).sum();
        assert!((total - 16.0).abs() < 1e-12, "area {total}");
        assert!(out.iter().all(|&t| area2(&pts, t) > 0.0));
        assert!(has_edge(&out, 3, 4));
        assert!(has_edge(&out, 5, 6));
        // Every point is used.
        for p in 0..pts.len() {
            assert!(out.iter().any(|t| t.contains(&p)), "point {p} unused");
        }
    }

    #[test]
    fn cdt_rejects_apex_outside_boundary() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [2.0, 2.0]];
        let masks = [0b101, 0b011, 0b110, 0];
        let out = triangulate(&CdtInput {
            pts: &pts,
            masks: &masks,
            boundary: &[0, 1, 2],
            interior: &[3],
            apex: 3,
            segments: &[],
        });
        assert!(out.is_none());
    }
}
