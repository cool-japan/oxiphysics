// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Edge flips for triangle meshes.
//!
//! One engine serves the three flip operators of the remeshing module:
//!
//! * [`flip_for_valence`] equalises vertex valences (isotropic remeshing),
//! * [`flip_edges_for_quality`] maximises the smallest angle,
//! * [`delaunay_edge_flip`] establishes the Delaunay criterion in the X-Z plane.
//!
//! The engine works in passes. Each pass rebuilds the edge topology of the current
//! triangles, visits the interior edges in order of first appearance and flips
//! every edge its policy accepts. A flip rewrites exactly two triangles; those are
//! marked dirty, and edges that touch a dirty triangle are skipped until the next
//! pass, because the topology snapshot no longer describes them. Passes repeat
//! until nothing flips.
//!
//! An edge `a -> b` of triangle `t_ab` with opposite vertex `c`, shared with the
//! triangle `t_ba` that holds `b -> a` with opposite vertex `d`, is replaced by the
//! edge `c -> d`:
//!
//! ```text
//!        c                      c
//!       / \                    /|\
//!      /   \                  / | \
//!     a --- b      ==>       a  |  b
//!      \   /                  \ | /
//!       \ /                    \|/
//!        d                      d
//!   (a,b,c) and (b,a,d)    (c,a,d) and (d,b,c)
//! ```
//!
//! The directed boundary edges `b -> c`, `c -> a`, `a -> d` and `d -> b` of the quad
//! are kept, so the winding of the triangles is preserved.
//!
//! A flip is only performed when
//!
//! * both triangles carry the edge in opposite directions (consistent winding),
//! * the opposite vertices differ and are not already joined by an edge, so the
//!   mesh never gains duplicate faces or non-manifold edges, and
//! * the quad is convex, so neither new triangle is inverted or degenerate.

use std::collections::{HashMap, HashSet};

use super::{cross3, dot3, len3, min_triangle_angle, sub3};

/// Upper bound on the number of passes. Every pass that flips something strictly
/// improves the quantity the policy optimises, so this only guards against
/// pathological floating-point cycles.
const MAX_FLIP_PASSES: usize = 1_000;

/// A flip must raise the smallest angle of the two triangles by more than this
/// (radians) to be performed.
const MIN_ANGLE_GAIN: f64 = 1e-6;

/// Relative rounding-error bound of the in-circle determinant. Shewchuk's bound
/// for the plain floating-point evaluation is about `1.1e-15` times the
/// permanent of the determinant; a margin is added.
const INCIRCLE_ERROR_FACTOR: f64 = 4e-15;

/// Ideal vertex valence in the interior of a triangle mesh.
const INTERIOR_VALENCE: i64 = 6;

/// Ideal vertex valence on the boundary of a triangle mesh.
const BOUNDARY_VALENCE: i64 = 4;

/// One use of a directed edge by a triangle.
#[derive(Debug, Clone, Copy)]
struct HalfEdgeUse {
    /// Index of the triangle.
    tri: usize,
    /// Start vertex of the directed edge.
    from: usize,
    /// End vertex of the directed edge.
    to: usize,
    /// Third vertex of the triangle.
    opposite: usize,
}

/// An undirected edge with all triangles that use it.
#[derive(Debug)]
struct EdgeRecord {
    /// Vertex pair, smaller index first.
    key: (usize, usize),
    /// Triangles using the edge, in triangle order.
    uses: Vec<HalfEdgeUse>,
}

/// Edge topology of a triangle list, with edges in order of first appearance.
pub(super) struct Topology {
    edges: Vec<EdgeRecord>,
    vertex_count: usize,
}

impl Topology {
    /// Collects the edges of `tris`. Triangles with an index outside
    /// `0..vertex_count` or with a repeated vertex are ignored.
    fn build(tris: &[[usize; 3]], vertex_count: usize) -> Self {
        let mut index: HashMap<(usize, usize), usize> = HashMap::new();
        let mut edges: Vec<EdgeRecord> = Vec::new();
        for (tri, &[a, b, c]) in tris.iter().enumerate() {
            let valid = a < vertex_count
                && b < vertex_count
                && c < vertex_count
                && a != b
                && b != c
                && a != c;
            if !valid {
                continue;
            }
            for (from, to, opposite) in [(a, b, c), (b, c, a), (c, a, b)] {
                let key = (from.min(to), from.max(to));
                let slot = *index.entry(key).or_insert_with(|| {
                    edges.push(EdgeRecord {
                        key,
                        uses: Vec::new(),
                    });
                    edges.len() - 1
                });
                edges[slot].uses.push(HalfEdgeUse {
                    tri,
                    from,
                    to,
                    opposite,
                });
            }
        }
        Self {
            edges,
            vertex_count,
        }
    }
}

/// An interior edge that may be flipped (see the module documentation).
#[derive(Debug, Clone, Copy)]
pub(super) struct FlipCandidate {
    /// Start of the directed edge in `tri_ab`.
    pub a: usize,
    /// End of the directed edge in `tri_ab`.
    pub b: usize,
    /// Vertex of `tri_ab` opposite to the edge.
    pub c: usize,
    /// Vertex of `tri_ba` opposite to the edge.
    pub d: usize,
    /// Triangle holding `a -> b`.
    pub tri_ab: usize,
    /// Triangle holding `b -> a`.
    pub tri_ba: usize,
}

impl FlipCandidate {
    /// The candidate of an edge shared by exactly two consistently oriented
    /// triangles, `None` for boundary, non-manifold and inconsistently oriented
    /// edges.
    fn of(record: &EdgeRecord) -> Option<Self> {
        let [first, second] = record.uses[..] else {
            return None;
        };
        if first.from != second.to || first.to != second.from {
            return None;
        }
        Some(Self {
            a: first.from,
            b: first.to,
            c: first.opposite,
            d: second.opposite,
            tri_ab: first.tri,
            tri_ba: second.tri,
        })
    }
}

/// Decides which edges the engine flips.
pub(super) trait FlipPolicy {
    /// Called at the start of every pass with the topology of that pass.
    fn begin_pass(&mut self, _topology: &Topology) {}

    /// Whether flipping the edge is valid and improves the mesh.
    fn accepts(&self, flip: &FlipCandidate) -> bool;

    /// Called after the engine flipped an edge this policy accepted.
    fn flipped(&mut self, _flip: &FlipCandidate) {}
}

/// Runs flip passes until no edge flips and returns the number of flips.
fn flip_edges<P: FlipPolicy>(
    tris: &mut [[usize; 3]],
    vertex_count: usize,
    policy: &mut P,
) -> usize {
    let mut total = 0;
    for _ in 0..MAX_FLIP_PASSES {
        let flips = flip_pass(tris, vertex_count, policy);
        if flips == 0 {
            break;
        }
        total += flips;
    }
    total
}

/// One pass over the edges of the current triangles.
fn flip_pass<P: FlipPolicy>(tris: &mut [[usize; 3]], vertex_count: usize, policy: &mut P) -> usize {
    let topology = Topology::build(tris, vertex_count);
    policy.begin_pass(&topology);
    // Every undirected edge of the mesh, kept in step with the flips of this pass.
    let mut existing: HashSet<(usize, usize)> =
        topology.edges.iter().map(|record| record.key).collect();
    let mut dirty = vec![false; tris.len()];
    let mut flips = 0;
    for record in &topology.edges {
        let Some(flip) = FlipCandidate::of(record) else {
            continue;
        };
        let FlipCandidate {
            a,
            b,
            c,
            d,
            tri_ab,
            tri_ba,
        } = flip;
        if dirty[tri_ab] || dirty[tri_ba] || c == d {
            continue;
        }
        let new_edge = (c.min(d), c.max(d));
        if existing.contains(&new_edge) || !policy.accepts(&flip) {
            continue;
        }
        tris[tri_ab] = [c, a, d];
        tris[tri_ba] = [d, b, c];
        dirty[tri_ab] = true;
        dirty[tri_ba] = true;
        existing.remove(&record.key);
        existing.insert(new_edge);
        policy.flipped(&flip);
        flips += 1;
    }
    flips
}

// ---------------------------------------------------------------------------
// Geometric validity
// ---------------------------------------------------------------------------

/// Whether the flip keeps the surface valid in 3-D: both new triangles are
/// non-degenerate and face the same side as both old triangles (positive dot
/// product of their normals). For a planar quad this is exactly strict
/// convexity of the quad.
fn flip_keeps_surface_valid(verts: &[[f64; 3]], flip: &FlipCandidate) -> bool {
    let [pa, pb, pc, pd] = [flip.a, flip.b, flip.c, flip.d].map(|v| verts[v]);
    let old = [
        cross3(sub3(pb, pa), sub3(pc, pa)), // (a, b, c)
        cross3(sub3(pa, pb), sub3(pd, pb)), // (b, a, d)
    ];
    let new = [
        cross3(sub3(pa, pc), sub3(pd, pc)), // (c, a, d)
        cross3(sub3(pb, pd), sub3(pc, pd)), // (d, b, c)
    ];
    new.iter().all(|&n| {
        old.iter()
            .all(|&o| len3(n) * len3(o) > 0.0 && dot3(n, o) > 0.0)
    })
}

/// Twice the signed area of `(a, b, c)` in the X-Z projection, positive when the
/// triangle is counter-clockwise in the `(x, z)` coordinate plane.
fn orient_xz(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    (a[0] - c[0]) * (b[2] - c[2]) - (a[2] - c[2]) * (b[0] - c[0])
}

/// In-circle determinant of `d` against `(a, b, c)` in the X-Z projection
/// together with the permanent that bounds its rounding error. The determinant is
/// positive when `d` lies inside the circle through `a`, `b`, `c` and the triangle
/// `(a, b, c)` is counter-clockwise in the `(x, z)` plane, negative when it is
/// clockwise.
fn incircle_xz(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> (f64, f64) {
    let (adx, adz) = (a[0] - d[0], a[2] - d[2]);
    let (bdx, bdz) = (b[0] - d[0], b[2] - d[2]);
    let (cdx, cdz) = (c[0] - d[0], c[2] - d[2]);
    let alift = adx * adx + adz * adz;
    let blift = bdx * bdx + bdz * bdz;
    let clift = cdx * cdx + cdz * cdz;
    let determinant = alift * (bdx * cdz - cdx * bdz)
        + blift * (cdx * adz - adx * cdz)
        + clift * (adx * bdz - bdx * adz);
    let permanent = ((bdx * cdz).abs() + (cdx * bdz).abs()) * alift
        + ((cdx * adz).abs() + (adx * cdz).abs()) * blift
        + ((adx * bdz).abs() + (bdx * adz).abs()) * clift;
    (determinant, permanent)
}

/// Returns `true` if `d` lies strictly inside the circumcircle of triangle
/// `(a, b, c)` in the X-Z projection.
///
/// The result does not depend on the winding of `(a, b, c)`. It is `false` for a
/// triangle that is degenerate in the projection and when `d` is so close to the
/// circle that the rounding error of the determinant could decide the answer.
pub(super) fn point_in_circumcircle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    let orientation = orient_xz(a, b, c);
    if orientation == 0.0 || !orientation.is_finite() {
        return false;
    }
    let (determinant, permanent) = incircle_xz(a, b, c, d);
    determinant * orientation.signum() > INCIRCLE_ERROR_FACTOR * permanent
}

// ---------------------------------------------------------------------------
// Valence flips
// ---------------------------------------------------------------------------

/// Flips edges so that vertex valences approach the ideal values (6 in the
/// interior, 4 on the boundary).
struct ValencePolicy<'a, L> {
    verts: &'a [[f64; 3]],
    locked: L,
    /// Number of edges at every vertex.
    degree: Vec<usize>,
    /// Whether a vertex lies on a boundary or non-manifold edge.
    on_boundary: Vec<bool>,
}

impl<L: Fn(usize, usize) -> bool> ValencePolicy<'_, L> {
    /// Squared deviation of `degree` from the ideal valence of vertex `v`.
    fn deviation(&self, v: usize, degree: usize) -> i64 {
        let ideal = if self.on_boundary[v] {
            BOUNDARY_VALENCE
        } else {
            INTERIOR_VALENCE
        };
        let delta = degree as i64 - ideal;
        delta * delta
    }
}

impl<L: Fn(usize, usize) -> bool> FlipPolicy for ValencePolicy<'_, L> {
    fn begin_pass(&mut self, topology: &Topology) {
        self.degree = vec![0; topology.vertex_count];
        self.on_boundary = vec![false; topology.vertex_count];
        for record in &topology.edges {
            let (u, v) = record.key;
            self.degree[u] += 1;
            self.degree[v] += 1;
            if record.uses.len() != 2 {
                self.on_boundary[u] = true;
                self.on_boundary[v] = true;
            }
        }
    }

    fn accepts(&self, flip: &FlipCandidate) -> bool {
        let FlipCandidate { a, b, c, d, .. } = *flip;
        // A vertex must keep at least three edges.
        if (self.locked)(a, b) || self.degree[a] <= 3 || self.degree[b] <= 3 {
            return false;
        }
        let before = self.deviation(a, self.degree[a])
            + self.deviation(b, self.degree[b])
            + self.deviation(c, self.degree[c])
            + self.deviation(d, self.degree[d]);
        let after = self.deviation(a, self.degree[a] - 1)
            + self.deviation(b, self.degree[b] - 1)
            + self.deviation(c, self.degree[c] + 1)
            + self.deviation(d, self.degree[d] + 1);
        after < before && flip_keeps_surface_valid(self.verts, flip)
    }

    fn flipped(&mut self, flip: &FlipCandidate) {
        self.degree[flip.a] -= 1;
        self.degree[flip.b] -= 1;
        self.degree[flip.c] += 1;
        self.degree[flip.d] += 1;
    }
}

/// Flips edges toward valence 6 (4 on the boundary) and returns the number of
/// flips.
///
/// Only edges shared by two consistently oriented triangles are flipped, and only
/// when the flip lowers the summed squared valence deviation of the four vertices
/// involved and keeps the quad convex (see the module documentation). `locked`
/// marks edges, given as `(a, b)` in either order, that must stay.
pub(super) fn flip_for_valence<L: Fn(usize, usize) -> bool>(
    tris: &mut [[usize; 3]],
    verts: &[[f64; 3]],
    locked: L,
) -> usize {
    let mut policy = ValencePolicy {
        verts,
        locked,
        degree: Vec::new(),
        on_boundary: Vec::new(),
    };
    flip_edges(tris, verts.len(), &mut policy)
}

// ---------------------------------------------------------------------------
// Quality flips
// ---------------------------------------------------------------------------

/// Flips edges that raise the smallest angle of the two adjacent triangles.
struct QualityPolicy<'a> {
    verts: &'a [[f64; 3]],
}

impl FlipPolicy for QualityPolicy<'_> {
    fn accepts(&self, flip: &FlipCandidate) -> bool {
        let [pa, pb, pc, pd] = [flip.a, flip.b, flip.c, flip.d].map(|v| self.verts[v]);
        let before = min_triangle_angle(pa, pb, pc).min(min_triangle_angle(pa, pb, pd));
        let after = min_triangle_angle(pc, pd, pa).min(min_triangle_angle(pc, pd, pb));
        after > before + MIN_ANGLE_GAIN && flip_keeps_surface_valid(self.verts, flip)
    }
}

/// Flip edges to improve triangle quality (maximize minimum angle).
///
/// An interior edge shared by two triangles is flipped when the flipped
/// configuration has a larger minimum angle (by more than `1e-6` rad). Flips keep
/// the winding of the triangles, never create an edge that already exists and are
/// refused unless the quad is convex (both new triangles face the same side as
/// both old ones). Passes repeat until no edge flips; the number of flips is
/// returned.
pub fn flip_edges_for_quality(tris: &mut [[usize; 3]], verts: &[[f64; 3]]) -> usize {
    flip_edges(tris, verts.len(), &mut QualityPolicy { verts })
}

// ---------------------------------------------------------------------------
// Delaunay flips
// ---------------------------------------------------------------------------

/// Flips edges whose opposite vertex lies inside the circumcircle (X-Z plane).
struct DelaunayPolicy<'a> {
    verts: &'a [[f64; 3]],
}

impl FlipPolicy for DelaunayPolicy<'_> {
    fn accepts(&self, flip: &FlipCandidate) -> bool {
        let [pa, pb, pc, pd] = [flip.a, flip.b, flip.c, flip.d].map(|v| self.verts[v]);
        // Both triangles of the current and of the flipped triangulation must have
        // the same orientation in the projection: the quad is strictly convex.
        let side = orient_xz(pa, pb, pc);
        if side == 0.0 || !side.is_finite() {
            return false;
        }
        let same_side = |orientation: f64| orientation * side > 0.0;
        same_side(orient_xz(pb, pa, pd))
            && same_side(orient_xz(pc, pa, pd))
            && same_side(orient_xz(pd, pb, pc))
            && point_in_circumcircle(pa, pb, pc, pd)
    }
}

/// Flip edges in a 2D-projected triangle mesh to satisfy the Delaunay criterion.
///
/// The mesh is projected onto the X-Z plane (the Y coordinate is ignored). An
/// interior edge is flipped if the opposite vertex of one adjacent triangle lies
/// strictly inside the circumcircle of the other, whatever the winding of the
/// triangles in the projection. Cocircular configurations (within rounding error)
/// are left alone. A flip is refused unless the quad is strictly convex in the
/// projection, so folded or degenerate triangle pairs are never touched. Flips
/// keep the winding of the triangles and never create an edge that already exists.
/// Passes repeat until no edge flips; the number of flips is returned. Operates
/// directly on the `tris` array (in-place).
pub fn delaunay_edge_flip(tris: &mut [[usize; 3]], verts: &[[f64; 3]]) -> usize {
    flip_edges(tris, verts.len(), &mut DelaunayPolicy { verts })
}

#[cfg(test)]
mod tests;
