// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface arrangement behind [`super::mesh_boolean`].
//!
//! 1. **Intersection.** Triangle pairs of the two operands are culled by their
//!    bounding boxes; each surviving pair yields its intersection segment
//!    (transversal case) or, for coplanar pairs, each triangle's edges clipped
//!    to the other triangle. All coordinates of a crossing are computed
//!    canonically (endpoints ordered lexicographically), so both triangles of a
//!    pair, and both triangles sharing a cut edge, receive bit-identical points.
//! 2. **Noding.** Inside every cut triangle the constraint segments are split at
//!    their mutual crossings and at points lying on them; points are snapped to
//!    the triangle's corners and edges.
//! 3. **Edge conformity.** Every point on a mesh edge is propagated to *both*
//!    triangles sharing that edge (keyed by the edge's exact endpoint
//!    coordinates), so a triangle that is not cut itself but borders a cut edge
//!    is split too and no T-junction is left behind.
//! 4. **Re-triangulation.** Each cut triangle is replaced by a constrained
//!    Delaunay triangulation of its noded graph ([`super::cdt`]), keeping the
//!    parent orientation.
//! 5. **Classification.** Every piece is classified by its centroid: coplanar
//!    pieces lying in a face of the other operand become
//!    [`TriClass::CoplanarSame`] / [`TriClass::CoplanarOpposite`]; all others go
//!    through the exact `orient3d` ray-parity test. Pieces are then selected per
//!    operation ([`keep_triangle`]) and the seam is welded.

use std::cmp::Ordering;
use std::collections::HashMap;

use oxiphysics_core::exact_predicates::{Orientation, orient2d};

use super::cdt::{CdtInput, triangulate};
use super::{
    MeshBooleanOp, SimpleMesh, TriClass, V3, add, cross, dot, exact_classify_point_vs_mesh,
    keep_triangle, length, remove_unused_vertices, scale, sub, weld_vertices,
};

/// Snapping tolerance relative to the coordinate scale of the operands.
const SNAP_REL: f64 = 1e-13;
/// Seam welding tolerance as a multiple of the snapping tolerance.
const WELD_FACTOR: f64 = 4.0;
/// Upper bound on noding sweeps per triangle (convergence takes 1–3).
const MAX_NODING_PASSES: usize = 64;
/// Lattice resolution for choosing a Steiner fan apex.
const APEX_LATTICE: usize = 9;
/// Minimum boundary clearance of an existing point used as fan apex, relative
/// to the longest parent edge.
const APEX_MIN_CLEARANCE: f64 = 1e-4;

type P2 = [f64; 2];

// ── small vector helpers ────────────────────────────────────────────────────

fn sub2(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn dot2(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn cross2(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn norm2(a: P2) -> f64 {
    dot2(a, a).sqrt()
}

fn dist3(a: V3, b: V3) -> f64 {
    length(sub(a, b))
}

/// Lexicographic order on coordinates (total, NaN-free input assumed).
fn lex_cmp(a: V3, b: V3) -> Ordering {
    a[0].total_cmp(&b[0])
        .then(a[1].total_cmp(&b[1]))
        .then(a[2].total_cmp(&b[2]))
}

/// `a + t (b - a)` that returns the endpoints exactly at `t = 0` and `t = 1`.
fn lerp_exact(a: V3, b: V3, t: f64) -> V3 {
    if t <= 0.0 {
        a
    } else if t >= 1.0 {
        b
    } else {
        add(a, scale(sub(b, a), t))
    }
}

/// Orthogonal projection of `p` onto the line through edge `a b`, evaluated
/// in the edge's canonical (lexicographic) direction so both triangles sharing
/// the edge obtain bit-identical results.
fn project_on_edge(a: V3, b: V3, p: V3) -> V3 {
    let (a, b) = if lex_cmp(a, b) == Ordering::Greater {
        (b, a)
    } else {
        (a, b)
    };
    let e = sub(b, a);
    let ee = dot(e, e);
    if ee > 0.0 {
        lerp_exact(a, b, dot(sub(p, a), e) / ee)
    } else {
        a
    }
}

fn centroid3(p: [V3; 3]) -> V3 {
    scale(add(add(p[0], p[1]), p[2]), 1.0 / 3.0)
}

/// Signed distance of `p` from the directed line `a → b` (positive on the left).
fn signed_dist2(a: P2, b: P2, p: P2) -> f64 {
    let e = sub2(b, a);
    let len = norm2(e);
    if len > 0.0 {
        cross2(e, sub2(p, a)) / len
    } else {
        norm2(sub2(p, a))
    }
}

/// Distance from `p` to the closed segment `[a, b]` in 2-D.
fn seg_dist2(a: P2, b: P2, p: P2) -> f64 {
    let e = sub2(b, a);
    let ee = dot2(e, e);
    let t = if ee > 0.0 {
        (dot2(sub2(p, a), e) / ee).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm2(sub2(p, [a[0] + t * e[0], a[1] + t * e[1]]))
}

// ── operands ────────────────────────────────────────────────────────────────

/// A boolean operand restricted to its valid, non-degenerate triangles.
struct Operand {
    mesh: SimpleMesh,
    normals: Vec<V3>,
    boxes: Vec<(V3, V3)>,
}

impl Operand {
    /// Drop triangles with out-of-range or repeated indices, non-finite
    /// corners or zero area, then compact the vertex list.
    fn new(mesh: &SimpleMesh) -> Self {
        let n = mesh.vertices.len();
        let finite = |v: &V3| v.iter().all(|c| c.is_finite());
        let mut tris = Vec::with_capacity(mesh.triangles.len());
        for t in &mesh.triangles {
            if t.iter().any(|&i| i >= n) || t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
                continue;
            }
            let [a, b, c] = t.map(|i| mesh.vertices[i]);
            if !(finite(&a) && finite(&b) && finite(&c)) {
                continue;
            }
            let len = length(cross(sub(b, a), sub(c, a)));
            if len.is_finite() && len > 0.0 {
                tris.push(*t);
            }
        }
        let mut clean = SimpleMesh::from_data(mesh.vertices.clone(), tris);
        remove_unused_vertices(&mut clean);
        let mut normals = Vec::with_capacity(clean.triangles.len());
        let mut boxes = Vec::with_capacity(clean.triangles.len());
        for t in &clean.triangles {
            let [a, b, c] = t.map(|i| clean.vertices[i]);
            let nrm = cross(sub(b, a), sub(c, a));
            normals.push(scale(nrm, 1.0 / length(nrm)));
            let mut lo = a;
            let mut hi = a;
            for p in [b, c] {
                for d in 0..3 {
                    lo[d] = lo[d].min(p[d]);
                    hi[d] = hi[d].max(p[d]);
                }
            }
            boxes.push((lo, hi));
        }
        Self {
            mesh: clean,
            normals,
            boxes,
        }
    }

    fn len(&self) -> usize {
        self.mesh.triangles.len()
    }

    fn corners(&self, t: usize) -> [V3; 3] {
        self.mesh.triangles[t].map(|i| self.mesh.vertices[i])
    }
}

/// Coordinate scale of both operands: the larger of their joint extent and
/// their largest absolute coordinate (rounding error grows with the latter).
fn coordinate_scale(a: &Operand, b: &Operand) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut mag = 0.0_f64;
    for v in a.mesh.vertices.iter().chain(&b.mesh.vertices) {
        for d in 0..3 {
            lo[d] = lo[d].min(v[d]);
            hi[d] = hi[d].max(v[d]);
            mag = mag.max(v[d].abs());
        }
    }
    let extent = (0..3).map(|d| hi[d] - lo[d]).fold(0.0_f64, f64::max);
    extent.max(mag).max(f64::MIN_POSITIVE)
}

/// Pairs `(i, j)` whose bounding boxes overlap (padded by `pad`), found with a
/// sweep along `x`.
fn candidate_pairs(a: &Operand, b: &Operand, pad: f64) -> Vec<(usize, usize)> {
    let mut order: Vec<usize> = (0..b.len()).collect();
    order.sort_by(|&x, &y| b.boxes[x].0[0].total_cmp(&b.boxes[y].0[0]));
    let mins: Vec<f64> = order.iter().map(|&j| b.boxes[j].0[0]).collect();
    let mut pairs = Vec::new();
    for (i, &(alo, ahi)) in a.boxes.iter().enumerate() {
        let end = mins.partition_point(|&m| m <= ahi[0] + pad);
        for &j in &order[..end] {
            let (blo, bhi) = b.boxes[j];
            if (0..3).all(|d| alo[d] <= bhi[d] + pad && blo[d] <= ahi[d] + pad) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

// ── triangle frames ─────────────────────────────────────────────────────────

/// Dominant-axis projection of a triangle that keeps its orientation
/// counter-clockwise in 2-D.
#[derive(Debug, Clone, Copy)]
struct Frame {
    axes: [usize; 2],
    c3: [V3; 3],
    c2: [P2; 3],
}

impl Frame {
    fn new(c3: [V3; 3], normal: V3) -> Self {
        let an = normal.map(f64::abs);
        let k = if an[0] >= an[1] && an[0] >= an[2] {
            0
        } else if an[1] >= an[2] {
            1
        } else {
            2
        };
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        let axes = if normal[k] >= 0.0 { [i, j] } else { [j, i] };
        let c2 = c3.map(|p| [p[axes[0]], p[axes[1]]]);
        Self { axes, c3, c2 }
    }

    fn proj(&self, p: V3) -> P2 {
        [p[self.axes[0]], p[self.axes[1]]]
    }

    fn is_valid(&self) -> bool {
        orient2d(self.c2[0], self.c2[1], self.c2[2]) == Orientation::Positive
    }
}

// ── pairwise intersection ───────────────────────────────────────────────────

/// Constraints collected for one triangle from the other operand.
#[derive(Debug, Clone, Default)]
struct TriWork {
    segs: Vec<(V3, V3)>,
    points: Vec<V3>,
    /// Triangles of the other operand that are coplanar with this one.
    coplanar: Vec<usize>,
}

impl TriWork {
    fn add(&mut self, p: V3, q: V3, tol: f64) {
        if dist3(p, q) > tol {
            self.segs.push((p, q));
        } else {
            self.points.push(p);
        }
    }
}

/// Outcome of intersecting two triangles.
enum PairHit {
    Disjoint,
    Coplanar,
    Segment(V3, V3),
    Point(V3),
}

/// Signed distances of `pts` from the plane through `origin` with unit normal
/// `n`, snapped to zero within `tol`.
fn plane_dists(n: V3, origin: V3, pts: [V3; 3], tol: f64) -> [f64; 3] {
    pts.map(|p| {
        let d = dot(n, sub(p, origin));
        if d.abs() <= tol { 0.0 } else { d }
    })
}

fn strictly_one_side(d: [f64; 3]) -> bool {
    d.iter().all(|&x| x > 0.0) || d.iter().all(|&x| x < 0.0)
}

/// Points where triangle `v` meets the other plane (`d` = snapped distances).
fn plane_crossings(v: [V3; 3], d: [f64; 3], out: &mut Vec<V3>) {
    for k in 0..3 {
        if d[k] == 0.0 {
            out.push(v[k]);
        }
    }
    for (k, l) in [(0, 1), (1, 2), (2, 0)] {
        if d[k] * d[l] < 0.0 {
            // Canonical evaluation: the same edge always yields the same point.
            let (p, dp, q, dq) = if lex_cmp(v[k], v[l]) == Ordering::Greater {
                (v[l], d[l], v[k], d[k])
            } else {
                (v[k], d[k], v[l], d[l])
            };
            out.push(lerp_exact(p, q, dp / (dp - dq)));
        }
    }
}

/// Extreme points of `pts` along `dir` as `((t_min, p_min), (t_max, p_max))`.
fn extremes(pts: &[V3], dir: V3) -> ((f64, V3), (f64, V3)) {
    let first = (dot(pts[0], dir), pts[0]);
    let mut lo = first;
    let mut hi = first;
    for &p in &pts[1..] {
        let t = dot(p, dir);
        if t < lo.0 {
            lo = (t, p);
        }
        if t > hi.0 {
            hi = (t, p);
        }
    }
    (lo, hi)
}

fn intersect_pair(ta: [V3; 3], na: V3, tb: [V3; 3], nb: V3, tol: f64) -> PairHit {
    let db = plane_dists(na, ta[0], tb, tol);
    if strictly_one_side(db) {
        return PairHit::Disjoint;
    }
    let da = plane_dists(nb, tb[0], ta, tol);
    if strictly_one_side(da) {
        return PairHit::Disjoint;
    }
    if db.iter().all(|&d| d == 0.0) || da.iter().all(|&d| d == 0.0) {
        return PairHit::Coplanar;
    }
    let dir = cross(na, nb);
    let len = length(dir);
    if len.is_nan() || len <= 0.0 {
        return PairHit::Coplanar;
    }
    let dir = scale(dir, 1.0 / len);
    let mut pa = Vec::with_capacity(4);
    plane_crossings(ta, da, &mut pa);
    let mut pb = Vec::with_capacity(4);
    plane_crossings(tb, db, &mut pb);
    if pa.is_empty() || pb.is_empty() {
        return PairHit::Disjoint;
    }
    let (a_lo, a_hi) = extremes(&pa, dir);
    let (b_lo, b_hi) = extremes(&pb, dir);
    let lo = if b_lo.0 > a_lo.0 { b_lo } else { a_lo };
    let hi = if b_hi.0 < a_hi.0 { b_hi } else { a_hi };
    if lo.0 > hi.0 + tol {
        return PairHit::Disjoint;
    }
    if lo.0 >= hi.0 || dist3(lo.1, hi.1) <= tol {
        return PairHit::Point(lo.1);
    }
    PairHit::Segment(lo.1, hi.1)
}

/// Clip segment `p q` to the triangle of `frame`, keeping ends that lie within
/// `tol` outside it.
fn clip_to_triangle(frame: &Frame, p: V3, q: V3, tol: f64) -> Option<(V3, V3)> {
    let (p, q) = if lex_cmp(p, q) == Ordering::Greater {
        (q, p)
    } else {
        (p, q)
    };
    let (p2, q2) = (frame.proj(p), frame.proj(q));
    let mut t0 = 0.0_f64;
    let mut t1 = 1.0_f64;
    for k in 0..3 {
        let a = frame.c2[k];
        let b = frame.c2[(k + 1) % 3];
        let s0 = signed_dist2(a, b, p2);
        let s1 = signed_dist2(a, b, q2);
        if s0 < -tol && s1 < -tol {
            return None;
        }
        // An end within `tol` outside is kept (it is snapped onto the edge
        // later); a farther end is cut where the segment crosses the *exact*
        // edge line, or collapses onto the near end when it never does.
        if s0 < -tol {
            t0 = t0.max(if s1 > 0.0 { s0 / (s0 - s1) } else { 1.0 });
        } else if s1 < -tol {
            t1 = t1.min(if s0 > 0.0 { s0 / (s0 - s1) } else { 0.0 });
        }
    }
    if t0 > t1 {
        return None;
    }
    Some((lerp_exact(p, q, t0), lerp_exact(p, q, t1)))
}

/// Intersect every candidate pair and collect per-triangle constraints.
fn collect_constraints(a: &Operand, b: &Operand, tol: f64) -> (Vec<TriWork>, Vec<TriWork>) {
    let mut work_a = vec![TriWork::default(); a.len()];
    let mut work_b = vec![TriWork::default(); b.len()];
    for (i, j) in candidate_pairs(a, b, tol) {
        let ta = a.corners(i);
        let tb = b.corners(j);
        match intersect_pair(ta, a.normals[i], tb, b.normals[j], tol) {
            PairHit::Disjoint => {}
            PairHit::Point(p) => {
                work_a[i].points.push(p);
                work_b[j].points.push(p);
            }
            PairHit::Segment(p, q) => {
                work_a[i].segs.push((p, q));
                work_b[j].segs.push((p, q));
            }
            PairHit::Coplanar => {
                let fa = Frame::new(ta, a.normals[i]);
                let fb = Frame::new(tb, b.normals[j]);
                for (k, l) in [(0, 1), (1, 2), (2, 0)] {
                    if let Some((p, q)) = clip_to_triangle(&fa, tb[k], tb[l], tol) {
                        work_a[i].add(p, q, tol);
                    }
                    if let Some((p, q)) = clip_to_triangle(&fb, ta[k], ta[l], tol) {
                        work_b[j].add(p, q, tol);
                    }
                }
                work_a[i].coplanar.push(j);
                work_b[j].coplanar.push(i);
            }
        }
    }
    (work_a, work_b)
}

// ── per-triangle planar graph ───────────────────────────────────────────────

/// Where a point lies relative to its parent triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Loc {
    /// Corner `k`.
    Corner(usize),
    /// Interior of edge `k` (corner `k` → corner `k + 1`) at parameter `t`.
    Edge(usize, f64),
    /// Strictly inside.
    Interior,
}

impl Loc {
    /// Bit `k` set iff the point lies on edge `k`.
    fn mask(self) -> u8 {
        match self {
            Loc::Corner(k) => (1 << k) | (1 << ((k + 2) % 3)),
            Loc::Edge(k, _) => 1 << k,
            Loc::Interior => 0,
        }
    }
}

/// Noded planar straight-line graph inside one parent triangle.
struct Local {
    frame: Frame,
    tol: f64,
    pts3: Vec<V3>,
    pts2: Vec<P2>,
    locs: Vec<Loc>,
    segs: Vec<[usize; 2]>,
}

impl Local {
    fn new(frame: Frame, tol: f64) -> Self {
        Self {
            frame,
            tol,
            pts3: frame.c3.to_vec(),
            pts2: frame.c2.to_vec(),
            locs: vec![Loc::Corner(0), Loc::Corner(1), Loc::Corner(2)],
            segs: Vec::new(),
        }
    }

    /// Position of `p` relative to the parent triangle. Points within `tol` of
    /// the boundary (or outside it) are snapped onto the boundary, so on-edge
    /// points lie on the parent edge itself.
    fn locate(&self, p: V3) -> (V3, Loc) {
        let p2 = self.frame.proj(p);
        let mut best = (f64::INFINITY, 0usize);
        for k in 0..3 {
            let sd = signed_dist2(self.frame.c2[k], self.frame.c2[(k + 1) % 3], p2);
            if sd < best.0 {
                best = (sd, k);
            }
        }
        if best.0 > self.tol {
            return (p, Loc::Interior);
        }
        let k = best.1;
        let (a, b) = (self.frame.c3[k], self.frame.c3[(k + 1) % 3]);
        let e = sub(b, a);
        let ee = dot(e, e);
        let len = ee.sqrt();
        let t = if ee > 0.0 {
            (dot(sub(p, a), e) / ee).clamp(0.0, 1.0)
        } else {
            0.0
        };
        if t * len <= self.tol {
            return (a, Loc::Corner(k));
        }
        if (1.0 - t) * len <= self.tol {
            return (b, Loc::Corner((k + 1) % 3));
        }
        (project_on_edge(a, b, p), Loc::Edge(k, t))
    }

    fn push(&mut self, p: V3, loc: Loc) -> usize {
        self.pts3.push(p);
        self.pts2.push(self.frame.proj(p));
        self.locs.push(loc);
        self.pts3.len() - 1
    }

    /// Insert `p`, merging it with an existing point within `tol`.
    fn add_point(&mut self, p: V3) -> usize {
        for k in 0..3 {
            if dist3(p, self.frame.c3[k]) <= self.tol {
                return k;
            }
        }
        let (q, loc) = self.locate(p);
        if let Loc::Corner(k) = loc {
            return k;
        }
        if let Some(i) = (3..self.pts3.len()).find(|&i| dist3(self.pts3[i], q) <= self.tol) {
            return i;
        }
        self.push(q, loc)
    }

    fn add_segment(&mut self, p: V3, q: V3) {
        let a = self.add_point(p);
        let b = self.add_point(q);
        if a != b {
            self.segs.push([a, b]);
        }
    }

    /// Insert the edge-`k` representative `r` shared with the neighbouring
    /// triangle, adopting its exact coordinates for a coincident point.
    fn add_edge_rep(&mut self, k: usize, r: V3) {
        for c in 0..3 {
            if dist3(r, self.frame.c3[c]) <= self.tol {
                return;
            }
        }
        let (a, b) = (self.frame.c3[k], self.frame.c3[(k + 1) % 3]);
        let e = sub(b, a);
        let ee = dot(e, e);
        let t = if ee > 0.0 {
            (dot(sub(r, a), e) / ee).clamp(0.0, 1.0)
        } else {
            0.0
        };
        if let Some(i) = (3..self.pts3.len()).find(|&i| dist3(self.pts3[i], r) <= self.tol) {
            if matches!(self.locs[i], Loc::Edge(..)) {
                self.pts3[i] = r;
                self.pts2[i] = self.frame.proj(r);
                self.locs[i] = Loc::Edge(k, t);
            }
            return;
        }
        self.push(r, Loc::Edge(k, t));
    }

    /// Sort/dedupe segments and drop the ones running along a parent edge
    /// (the boundary is constrained implicitly).
    fn canonicalize(&mut self) {
        let locs = &self.locs;
        let mut segs: Vec<[usize; 2]> = self
            .segs
            .iter()
            .filter(|s| s[0] != s[1] && locs[s[0]].mask() & locs[s[1]].mask() == 0)
            .map(|s| if s[0] < s[1] { *s } else { [s[1], s[0]] })
            .collect();
        segs.sort_unstable();
        segs.dedup();
        self.segs = segs;
    }

    /// Replace segment `[a, b]` by the chain through `cuts` (point ids), ordered
    /// along the segment.
    fn chain(&self, a: usize, b: usize, cuts: &mut Vec<usize>, out: &mut Vec<[usize; 2]>) {
        let pa = self.pts2[a];
        let d = sub2(self.pts2[b], pa);
        cuts.sort_by(|&x, &y| {
            dot2(sub2(self.pts2[x], pa), d).total_cmp(&dot2(sub2(self.pts2[y], pa), d))
        });
        cuts.dedup();
        let mut prev = a;
        for &w in cuts.iter().chain(std::iter::once(&b)) {
            if w != prev {
                out.push([prev, w]);
                prev = w;
            }
        }
    }

    /// Split segments at points lying on their interior. Returns whether any
    /// segment was split.
    fn split_at_points(&mut self) -> bool {
        let mut out = Vec::with_capacity(self.segs.len());
        let mut changed = false;
        let segs = std::mem::take(&mut self.segs);
        for [a, b] in segs {
            let pa = self.pts2[a];
            let d = sub2(self.pts2[b], pa);
            let len = norm2(d);
            let mut cuts = Vec::new();
            if len > 0.0 {
                for w in 0..self.pts2.len() {
                    if w == a || w == b {
                        continue;
                    }
                    let r = sub2(self.pts2[w], pa);
                    let along = dot2(r, d) / len;
                    if along <= self.tol || along >= len - self.tol {
                        continue;
                    }
                    if (cross2(d, r) / len).abs() <= self.tol {
                        cuts.push(w);
                    }
                }
            }
            if cuts.is_empty() {
                out.push([a, b]);
            } else {
                changed = true;
                self.chain(a, b, &mut cuts, &mut out);
            }
        }
        self.segs = out;
        changed
    }

    /// Insert the crossing points of properly crossing segment pairs and split
    /// both segments there. Returns whether anything changed.
    fn split_crossings(&mut self) -> bool {
        let n = self.segs.len();
        let mut cuts: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut changed = false;
        for i in 0..n {
            for j in (i + 1)..n {
                let [a, b] = self.segs[i];
                let [c, d] = self.segs[j];
                if a == c || a == d || b == c || b == d {
                    continue;
                }
                let (pa, pb, pc, pd) = (self.pts2[a], self.pts2[b], self.pts2[c], self.pts2[d]);
                let sc = signed_dist2(pa, pb, pc);
                let sd = signed_dist2(pa, pb, pd);
                if !((sc > self.tol && sd < -self.tol) || (sc < -self.tol && sd > self.tol)) {
                    continue;
                }
                let sa = signed_dist2(pc, pd, pa);
                let sb = signed_dist2(pc, pd, pb);
                if !((sa > self.tol && sb < -self.tol) || (sa < -self.tol && sb > self.tol)) {
                    continue;
                }
                let x = lerp_exact(self.pts3[a], self.pts3[b], sa / (sa - sb));
                let w = self.add_point(x);
                if w != a && w != b {
                    cuts[i].push(w);
                }
                if w != c && w != d {
                    cuts[j].push(w);
                }
                changed = true;
            }
        }
        if changed {
            let segs = std::mem::take(&mut self.segs);
            let mut out = Vec::with_capacity(segs.len() * 2);
            for (s, mut cut) in segs.into_iter().zip(cuts) {
                if cut.is_empty() {
                    out.push(s);
                } else {
                    self.chain(s[0], s[1], &mut cut, &mut out);
                }
            }
            self.segs = out;
        }
        changed
    }

    /// Node the graph: split until segments meet only at shared endpoints.
    fn node(&mut self) {
        for _ in 0..MAX_NODING_PASSES {
            self.canonicalize();
            let split = self.split_at_points();
            let crossed = self.split_crossings();
            if !split && !crossed {
                break;
            }
        }
        self.canonicalize();
    }

    /// Clearance of the 2-D point `c` from the boundary, segments and points.
    fn clearance(&self, c: P2) -> f64 {
        let mut m = f64::INFINITY;
        for k in 0..3 {
            m = m.min(signed_dist2(
                self.frame.c2[k],
                self.frame.c2[(k + 1) % 3],
                c,
            ));
        }
        for s in &self.segs {
            m = m.min(seg_dist2(self.pts2[s[0]], self.pts2[s[1]], c));
        }
        for p in &self.pts2 {
            m = m.min(norm2(sub2(*p, c)));
        }
        m
    }

    /// Smallest signed distance of the 2-D point `c` to the parent edges.
    fn boundary_clearance(&self, c: P2) -> f64 {
        (0..3)
            .map(|k| signed_dist2(self.frame.c2[k], self.frame.c2[(k + 1) % 3], c))
            .fold(f64::INFINITY, f64::min)
    }

    /// Counter-clockwise boundary cycle: corners and on-edge points in order.
    fn boundary_cycle(&self) -> Vec<usize> {
        let mut boundary = Vec::with_capacity(self.pts3.len());
        for k in 0..3 {
            boundary.push(k);
            let mut on_edge: Vec<(f64, usize)> = (3..self.pts3.len())
                .filter_map(|i| match self.locs[i] {
                    Loc::Edge(e, t) if e == k => Some((t, i)),
                    _ => None,
                })
                .collect();
            on_edge.sort_by(|x, y| x.0.total_cmp(&y.0));
            boundary.extend(on_edge.into_iter().map(|(_, i)| i));
        }
        boundary
    }

    /// Whether every fan triangle `(c, b_i, b_{i+1})` is strictly
    /// counter-clockwise (exact predicate).
    fn fan_is_valid(&self, c: P2, boundary: &[usize]) -> bool {
        let m = boundary.len();
        (0..m).all(|i| {
            orient2d(c, self.pts2[boundary[i]], self.pts2[boundary[(i + 1) % m]])
                == Orientation::Positive
        })
    }

    /// Choose the fan apex. An existing interior point is used when it keeps
    /// a clearance of at least `APEX_MIN_CLEARANCE ×` the longest parent edge
    /// from the boundary (a fan from a point hugging an edge is
    /// ill-conditioned against closely spaced boundary points); otherwise the
    /// interior lattice point with the best clearance from the boundary,
    /// segments and points whose fan is valid is inserted as a Steiner point.
    fn choose_apex(&mut self, boundary: &[usize]) -> Option<usize> {
        let longest = (0..3)
            .map(|k| norm2(sub2(self.frame.c2[(k + 1) % 3], self.frame.c2[k])))
            .fold(0.0_f64, f64::max);
        let min_clearance = APEX_MIN_CLEARANCE * longest;
        let mut interior: Vec<(f64, usize)> = (3..self.pts3.len())
            .filter(|&i| self.locs[i] == Loc::Interior)
            .map(|i| (self.boundary_clearance(self.pts2[i]), i))
            .collect();
        interior.sort_by(|x, y| y.0.total_cmp(&x.0));
        if let Some(&(_, i)) = interior.iter().find(|&&(clear, i)| {
            clear >= min_clearance && self.fan_is_valid(self.pts2[i], boundary)
        }) {
            return Some(i);
        }
        let n = APEX_LATTICE;
        let mut candidates: Vec<(f64, V3)> = Vec::new();
        for i in 1..n {
            for j in 1..(n - i) {
                let k = n - i - j;
                let w = [i as f64, j as f64, k as f64].map(|x| x / n as f64);
                let c3 = add(
                    add(scale(self.frame.c3[0], w[0]), scale(self.frame.c3[1], w[1])),
                    scale(self.frame.c3[2], w[2]),
                );
                candidates.push((self.clearance(self.frame.proj(c3)), c3));
            }
        }
        candidates.sort_by(|x, y| y.0.total_cmp(&x.0));
        for (_, c3) in candidates {
            if !self.fan_is_valid(self.frame.proj(c3), boundary) {
                continue;
            }
            let (_, loc) = self.locate(c3);
            if loc != Loc::Interior {
                continue;
            }
            let id = self.add_point(c3);
            if self.locs[id] == Loc::Interior {
                // Split any constraint the Steiner point happens to touch.
                self.node();
                return Some(id);
            }
        }
        // Fall back to any interior point whose fan is valid.
        interior
            .iter()
            .find(|&&(_, i)| self.fan_is_valid(self.pts2[i], boundary))
            .map(|&(_, i)| i)
    }

    /// Constrained triangulation of the graph, or `None` when the parent
    /// triangle needs no split (or has no valid fan apex, which only happens
    /// for needle-like parents).
    fn triangulate(&mut self) -> Option<Vec<[usize; 3]>> {
        if self.pts3.len() <= 3 || !self.frame.is_valid() {
            return None;
        }
        let boundary = self.boundary_cycle();
        let apex = self.choose_apex(&boundary)?;
        // A Steiner apex re-nodes the graph; refresh the cycle in case that
        // touched the boundary.
        let boundary = self.boundary_cycle();
        let interior: Vec<usize> = (3..self.pts3.len())
            .filter(|&i| self.locs[i] == Loc::Interior)
            .collect();
        let masks: Vec<u8> = self.locs.iter().map(|l| l.mask()).collect();
        triangulate(&CdtInput {
            pts: &self.pts2,
            masks: &masks,
            boundary: &boundary,
            interior: &interior,
            apex,
            segments: &self.segs,
        })
    }
}

// ── subdivision of a whole operand ──────────────────────────────────────────

/// One output piece of a parent triangle.
struct Piece {
    p3: [V3; 3],
    c2: P2,
}

type EdgeKey = [u64; 6];

/// Exact-coordinate key of an undirected edge plus its canonical endpoints.
fn edge_key(p: V3, q: V3) -> (EdgeKey, V3, V3) {
    let (a, b) = if lex_cmp(p, q) == Ordering::Greater {
        (q, p)
    } else {
        (p, q)
    };
    // `+ 0.0` folds -0.0 into +0.0 so equal coordinates share a key.
    let bits = |x: f64| (x + 0.0).to_bits();
    (
        [
            bits(a[0]),
            bits(a[1]),
            bits(a[2]),
            bits(b[0]),
            bits(b[1]),
            bits(b[2]),
        ],
        a,
        b,
    )
}

/// Split every triangle of `op` along its constraints. Returns the frames and
/// the pieces of each triangle (one whole-triangle piece when uncut).
fn subdivide(op: &Operand, work: &[TriWork], tol: f64) -> (Vec<Frame>, Vec<Vec<Piece>>) {
    let n = op.len();
    let frames: Vec<Frame> = (0..n)
        .map(|t| Frame::new(op.corners(t), op.normals[t]))
        .collect();
    let mut locals: Vec<Option<Local>> = (0..n).map(|_| None).collect();
    for t in 0..n {
        let w = &work[t];
        if w.segs.is_empty() && w.points.is_empty() {
            continue;
        }
        let mut local = Local::new(frames[t], tol);
        for &(p, q) in &w.segs {
            local.add_segment(p, q);
        }
        for &p in &w.points {
            local.add_point(p);
        }
        local.node();
        locals[t] = Some(local);
    }

    // Edge conformity: gather the on-edge points of every mesh edge and hand
    // one representative per cluster to every triangle sharing the edge.
    let mut registry: HashMap<EdgeKey, (V3, V3, Vec<V3>)> = HashMap::new();
    for local in locals.iter().flatten() {
        for (i, loc) in local.locs.iter().enumerate() {
            if let Loc::Edge(k, _) = *loc {
                let (key, a, b) = edge_key(local.frame.c3[k], local.frame.c3[(k + 1) % 3]);
                registry
                    .entry(key)
                    .or_insert_with(|| (a, b, Vec::new()))
                    .2
                    .push(local.pts3[i]);
            }
        }
    }
    let mut reps: HashMap<EdgeKey, Vec<V3>> = HashMap::with_capacity(registry.len());
    for (key, (a, b, mut pts)) in registry {
        let e = sub(b, a);
        pts.sort_by(|x, y| dot(sub(*x, a), e).total_cmp(&dot(sub(*y, a), e)));
        let mut cluster: Vec<V3> = Vec::with_capacity(pts.len());
        for p in pts {
            if cluster.last().is_none_or(|&last| dist3(last, p) > tol) {
                cluster.push(p);
            }
        }
        reps.insert(key, cluster);
    }
    for t in 0..n {
        for k in 0..3 {
            let (key, _, _) = edge_key(frames[t].c3[k], frames[t].c3[(k + 1) % 3]);
            if let Some(list) = reps.get(&key) {
                let local = locals[t].get_or_insert_with(|| Local::new(frames[t], tol));
                for &r in list {
                    local.add_edge_rep(k, r);
                }
            }
        }
    }

    let mut pieces = Vec::with_capacity(n);
    for t in 0..n {
        let whole = || {
            vec![Piece {
                p3: frames[t].c3,
                c2: centroid2(frames[t].c2),
            }]
        };
        let split = match locals[t].as_mut() {
            Some(local) => {
                local.node();
                local.triangulate().map(|tris| {
                    tris.iter()
                        .map(|tri| Piece {
                            p3: tri.map(|i| local.pts3[i]),
                            c2: centroid2(tri.map(|i| local.pts2[i])),
                        })
                        .collect::<Vec<_>>()
                })
            }
            None => None,
        };
        pieces.push(split.unwrap_or_else(whole));
    }
    (frames, pieces)
}

fn centroid2(p: [P2; 3]) -> P2 {
    [
        (p[0][0] + p[1][0] + p[2][0]) / 3.0,
        (p[0][1] + p[1][1] + p[2][1]) / 3.0,
    ]
}

/// Inclusive exact point-in-triangle test (either orientation).
fn in_triangle_2d(p: P2, t: [P2; 3]) -> bool {
    let o = [
        orient2d(t[0], t[1], p),
        orient2d(t[1], t[2], p),
        orient2d(t[2], t[0], p),
    ];
    !(o.contains(&Orientation::Positive) && o.contains(&Orientation::Negative))
}

/// Classify a piece of triangle `t` of `own` against `other`.
fn classify_piece(
    piece: &Piece,
    frame: &Frame,
    normal: V3,
    coplanar: &[usize],
    other: &Operand,
) -> TriClass {
    for &u in coplanar {
        let corners = other.corners(u).map(|p| frame.proj(p));
        if in_triangle_2d(piece.c2, corners) {
            return if dot(normal, other.normals[u]) >= 0.0 {
                TriClass::CoplanarSame
            } else {
                TriClass::CoplanarOpposite
            };
        }
    }
    if exact_classify_point_vs_mesh(&other.mesh, centroid3(piece.p3)) {
        TriClass::Inside
    } else {
        TriClass::Outside
    }
}

/// One operand after subdivision.
struct Side<'a> {
    own: &'a Operand,
    frames: Vec<Frame>,
    pieces: Vec<Vec<Piece>>,
    work: Vec<TriWork>,
}

impl Side<'_> {
    /// Append the pieces that `op` keeps against `other`, flipped when `flip`.
    fn emit(
        &self,
        other: &Operand,
        op: MeshBooleanOp,
        is_a: bool,
        flip: bool,
        out: &mut SimpleMesh,
    ) {
        for t in 0..self.own.len() {
            for piece in &self.pieces[t] {
                let class = classify_piece(
                    piece,
                    &self.frames[t],
                    self.own.normals[t],
                    &self.work[t].coplanar,
                    other,
                );
                if !keep_triangle(op, is_a, class) {
                    continue;
                }
                let base = out.vertices.len();
                out.vertices.extend_from_slice(&piece.p3);
                out.triangles.push(if flip {
                    [base + 1, base, base + 2]
                } else {
                    [base, base + 1, base + 2]
                });
            }
        }
    }
}

/// Boolean `op` of `mesh_a` and `mesh_b` with conforming seams.
pub(super) fn boolean(mesh_a: &SimpleMesh, mesh_b: &SimpleMesh, op: MeshBooleanOp) -> SimpleMesh {
    let a = Operand::new(mesh_a);
    let b = Operand::new(mesh_b);
    let tol = SNAP_REL * coordinate_scale(&a, &b);
    let (work_a, work_b) = collect_constraints(&a, &b, tol);
    let (frames_a, pieces_a) = subdivide(&a, &work_a, tol);
    let (frames_b, pieces_b) = subdivide(&b, &work_b, tol);
    let side_a = Side {
        own: &a,
        frames: frames_a,
        pieces: pieces_a,
        work: work_a,
    };
    let side_b = Side {
        own: &b,
        frames: frames_b,
        pieces: pieces_b,
        work: work_b,
    };

    let mut out = SimpleMesh::new();
    // For a difference, B's pieces inside A bound the carved cavity and are
    // emitted with reversed winding.
    let flip_b = matches!(op, MeshBooleanOp::Difference);
    side_a.emit(&b, op, true, false, &mut out);
    side_b.emit(&a, op, false, flip_b, &mut out);

    let mut welded = weld_vertices(&out, WELD_FACTOR * tol);
    welded
        .triangles
        .retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
    cancel_opposite_pairs(&mut welded.triangles);
    remove_unused_vertices(&mut welded);
    welded
}

/// Remove pairs of triangles with the same three vertices and opposite
/// orientation. As 2-chains they sum to zero (`∂` is linear and
/// `[a,b,c] + [a,c,b] = 0`), so dropping both leaves the boundary of the
/// result unchanged; such pairs appear when welding collapses two
/// near-coincident slivers of a near-degenerate configuration.
fn cancel_opposite_pairs(tris: &mut Vec<[usize; 3]>) {
    /// Sorted vertex triple plus the cyclic-order parity.
    type CanonTri = ((usize, usize, usize), bool);
    // Canonical form: rotate the smallest index first; `parity` tells which
    // of the two cyclic orders the triangle has.
    let canon = |t: &[usize; 3]| -> CanonTri {
        let r = if t[0] <= t[1] && t[0] <= t[2] {
            *t
        } else if t[1] <= t[2] {
            [t[1], t[2], t[0]]
        } else {
            [t[2], t[0], t[1]]
        };
        ((r[0], r[1].min(r[2]), r[1].max(r[2])), r[1] < r[2])
    };
    let mut open: HashMap<CanonTri, Vec<usize>> = HashMap::new();
    let mut remove = vec![false; tris.len()];
    for (i, t) in tris.iter().enumerate() {
        let (key, parity) = canon(t);
        if let Some(j) = open.get_mut(&(key, !parity)).and_then(Vec::pop) {
            remove[i] = true;
            remove[j] = true;
        } else {
            open.entry((key, parity)).or_default().push(i);
        }
    }
    let mut keep = remove.iter().map(|r| !r);
    tris.retain(|_| keep.next().unwrap_or(true));
}
