// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact 2-D set difference `a \ b` of convex polygons.
//!
//! The difference of two convex polygons is in general **not** a single
//! polygon: `b` may cut `a` in two, or lie strictly inside `a` and leave a
//! ring.  [`polygon2d_difference_pieces`] therefore returns the result as a
//! list of pairwise interior-disjoint convex polygons.
//! [`polygon2d_difference`] returns the outline of the result for the common
//! case in which it is a single simple polygon and reports an error when it
//! is not.
//!
//! Both functions accept either winding order and work in a normalised
//! coordinate frame, so they are independent of the absolute scale and
//! position of the input.

use std::collections::BTreeMap;

use super::functions::{
    clip_polygon_2d, clip_polygon_half_plane_2d, polygon2d_area, polygon2d_is_convex,
    polygon2d_signed_area,
};
use crate::{Error, Result};

/// Vertices closer than this (in the normalised frame, where the input has
/// unit extent) are treated as the same vertex when pieces are merged.
const COINCIDENT_TOLERANCE: f64 = 1e-9;
/// Pieces smaller than this area (normalised frame) are degenerate slivers,
/// for example the zero-width strip left when an edge of `b` is collinear
/// with an edge of `a`.
const SLIVER_AREA_TOLERANCE: f64 = 1e-12;
/// Vertices of a clipped piece closer than this (normalised frame) are
/// duplicates produced by clipping through an existing vertex.
const REPEATED_VERTEX_TOLERANCE: f64 = 1e-12;
/// Maximum `sin(angle)` between consecutive outline edges for the shared
/// vertex to be dropped as collinear.
const COLLINEAR_SINE_TOLERANCE: f64 = 1e-10;

/// Affine map between the input coordinates and a frame in which the
/// combined bounding box of both polygons has unit extent.
struct Frame {
    origin: [f64; 2],
    extent: f64,
}

impl Frame {
    /// Fit a frame to the combined bounding box of `a` and `b`.
    ///
    /// Returns `None` if any coordinate is not finite or if all points
    /// coincide.
    fn fit(a: &[[f64; 2]], b: &[[f64; 2]]) -> Option<Self> {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for p in a.iter().chain(b.iter()) {
            if !(p[0].is_finite() && p[1].is_finite()) {
                return None;
            }
            for ((low, high), &coordinate) in lo.iter_mut().zip(hi.iter_mut()).zip(p.iter()) {
                *low = low.min(coordinate);
                *high = high.max(coordinate);
            }
        }
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]);
        if !(extent.is_finite() && extent > 0.0) {
            return None;
        }
        Some(Self { origin: lo, extent })
    }

    fn to_unit(&self, p: [f64; 2]) -> [f64; 2] {
        [
            (p[0] - self.origin[0]) / self.extent,
            (p[1] - self.origin[1]) / self.extent,
        ]
    }

    fn to_input(&self, p: [f64; 2]) -> [f64; 2] {
        [
            p[0] * self.extent + self.origin[0],
            p[1] * self.extent + self.origin[1],
        ]
    }

    fn ring_to_unit(&self, ring: &[[f64; 2]]) -> Vec<[f64; 2]> {
        ring.iter().map(|&p| self.to_unit(p)).collect()
    }

    fn ring_to_input(&self, ring: &[[f64; 2]]) -> Vec<[f64; 2]> {
        ring.iter().map(|&p| self.to_input(p)).collect()
    }
}

/// Copy of `poly` wound counter-clockwise.
fn counter_clockwise(poly: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut ring = poly.to_vec();
    if polygon2d_signed_area(&ring) < 0.0 {
        ring.reverse();
    }
    ring
}

fn distance_2d(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

/// Remove consecutive (and wrap-around) vertices closer than `tol`.
fn drop_repeated_vertices(ring: Vec<[f64; 2]>, tol: f64) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(ring.len());
    for p in ring {
        if out.last().is_none_or(|&q| distance_2d(q, p) > tol) {
            out.push(p);
        }
    }
    while out.len() > 1 && distance_2d(out[0], out[out.len() - 1]) <= tol {
        out.pop();
    }
    out
}

/// Decompose `a \ b` in the normalised frame.  Both rings must already be
/// counter-clockwise and fit in the unit square.
fn difference_pieces_unit(subject: Vec<[f64; 2]>, cutter: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    let subject_area = polygon2d_area(&subject);
    if subject_area <= SLIVER_AREA_TOLERANCE {
        return Vec::new();
    }
    let overlap_area = polygon2d_area(&clip_polygon_2d(&subject, cutter));
    if overlap_area <= SLIVER_AREA_TOLERANCE {
        // Disjoint or only touching: nothing is removed.
        return vec![subject];
    }
    if subject_area - overlap_area <= SLIVER_AREA_TOLERANCE {
        // `a` lies inside `b`: nothing is left.
        return Vec::new();
    }
    // piece_i = subject ∩ inside(edge_0 .. edge_{i-1}) ∩ outside(edge_i).
    // A point of `a \ b` is outside at least one edge of the convex `b`; it
    // belongs to the piece of the first edge it is outside of, so the pieces
    // are disjoint and their union is `a \ b`.
    let n = cutter.len();
    let mut pieces: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut inside_so_far = subject;
    for i in 0..n {
        let edge_start = cutter[i];
        let edge_end = cutter[(i + 1) % n];
        // The closed half-plane to the right of `start → end` is the closed
        // half-plane to the left of the reversed edge.
        let outside = clip_polygon_half_plane_2d(&inside_so_far, edge_end, edge_start);
        let outside = drop_repeated_vertices(outside, REPEATED_VERTEX_TOLERANCE);
        if polygon2d_area(&outside) > SLIVER_AREA_TOLERANCE {
            pieces.push(outside);
        }
        inside_so_far = clip_polygon_half_plane_2d(&inside_so_far, edge_start, edge_end);
        if inside_so_far.len() < 3 {
            break;
        }
    }
    pieces
}

/// Set difference `a \ b` of two convex polygons as a list of pairwise
/// interior-disjoint convex polygons whose union is `a \ b`.
///
/// Both polygons must be convex; either winding order is accepted and the
/// returned pieces are counter-clockwise.  The result is:
///
/// * empty if `a` has fewer than three vertices, is degenerate, is contained
///   in `b`, or any coordinate is not finite;
/// * `[a]` (re-wound counter-clockwise) if `b` has fewer than three vertices
///   or does not overlap `a` (disjoint or merely touching);
/// * otherwise one piece per edge of `b` that cuts a non-empty region out of
///   `a`: the part of `a` outside that edge and inside all earlier edges.
///   Zero-area slivers (for example from collinear edges) are dropped.
pub fn polygon2d_difference_pieces(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<Vec<[f64; 2]>> {
    if a.len() < 3 {
        return Vec::new();
    }
    let Some(frame) = Frame::fit(a, b) else {
        return Vec::new();
    };
    let subject = counter_clockwise(&frame.ring_to_unit(a));
    let pieces = if b.len() < 3 {
        vec![subject]
    } else {
        let cutter = counter_clockwise(&frame.ring_to_unit(b));
        difference_pieces_unit(subject, &cutter)
    };
    pieces
        .iter()
        .map(|piece| frame.ring_to_input(piece))
        .collect()
}

/// Set difference `a \ b` of two convex polygons as a single outline.
///
/// This is the outline (counter-clockwise, collinear vertices removed) of the
/// union of [`polygon2d_difference_pieces`] when that union is one simple
/// polygon, which covers the usual case of `b` removing a corner or a side
/// of `a`.  The outline of a non-convex result (an L-shape, for instance) is
/// returned as is.  An empty vector means `a \ b` is empty.
///
/// # Errors
///
/// Returns [`Error::Unsupported`] when
///
/// * `a` or `b` is not convex (the decomposition into pieces is only valid
///   for convex operands), or
/// * `a \ b` cannot be written as one vertex list: `b` splits `a` into
///   several regions, or lies strictly inside `a` and leaves a hole.  Use
///   [`polygon2d_difference_pieces`] for those.
pub fn polygon2d_difference(a: &[[f64; 2]], b: &[[f64; 2]]) -> Result<Vec<[f64; 2]>> {
    if a.len() < 3 {
        return Ok(Vec::new());
    }
    let Some(frame) = Frame::fit(a, b) else {
        return Ok(Vec::new());
    };
    let subject = counter_clockwise(&frame.ring_to_unit(a));
    if !polygon2d_is_convex(&subject) {
        return Err(Error::unsupported(
            "polygon2d_difference",
            "the polygon a is not convex",
        ));
    }
    let pieces = if b.len() < 3 {
        vec![subject]
    } else {
        let cutter = counter_clockwise(&frame.ring_to_unit(b));
        if !polygon2d_is_convex(&cutter) {
            return Err(Error::unsupported(
                "polygon2d_difference",
                "the polygon b is not convex",
            ));
        }
        difference_pieces_unit(subject, &cutter)
    };
    if pieces.is_empty() {
        return Ok(Vec::new());
    }
    match merge_pieces_into_outline(&pieces) {
        Some(outline) => Ok(frame.ring_to_input(&outline)),
        None => Err(Error::unsupported(
            "polygon2d_difference",
            format!(
                "a \\ b consists of {} convex pieces that do not form one simple polygon \
                 (b splits a or leaves a hole); use polygon2d_difference_pieces",
                pieces.len()
            ),
        )),
    }
}

/// Add the directed edge `from → to`, cancelling an existing opposite edge.
///
/// Two pieces that share a boundary contribute the shared edge once in each
/// direction; cancelling those pairs leaves only the outer boundary.
fn add_directed_edge(edges: &mut BTreeMap<(usize, usize), u32>, from: usize, to: usize) {
    if from == to {
        return;
    }
    match edges.get(&(to, from)).copied() {
        Some(1) => {
            edges.remove(&(to, from));
        }
        Some(count) => {
            edges.insert((to, from), count - 1);
        }
        None => {
            *edges.entry((from, to)).or_insert(0) += 1;
        }
    }
}

/// Parameter `t ∈ (0, 1)` of `p` along `a → b` if `p` lies on the interior of
/// that segment within `tol`.
fn interior_segment_parameter(a: [f64; 2], b: [f64; 2], p: [f64; 2], tol: f64) -> Option<f64> {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];
    if len_sq <= tol * tol {
        return None;
    }
    let len = len_sq.sqrt();
    let ap = [p[0] - a[0], p[1] - a[1]];
    let t = (ap[0] * ab[0] + ap[1] * ab[1]) / len_sq;
    let off_line = (ap[0] * ab[1] - ap[1] * ab[0]).abs() / len;
    (off_line <= tol && t * len > tol && (1.0 - t) * len > tol).then_some(t)
}

/// Remove vertices whose two incident edges are collinear and point the same
/// way.  A triangle is never reduced further.
fn remove_collinear_vertices(mut ring: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    loop {
        let n = ring.len();
        if n < 4 {
            return ring;
        }
        let removable = (0..n).find(|&i| {
            let prev = ring[(i + n - 1) % n];
            let cur = ring[i];
            let next = ring[(i + 1) % n];
            let u = [cur[0] - prev[0], cur[1] - prev[1]];
            let v = [next[0] - cur[0], next[1] - cur[1]];
            let (len_u, len_v) = (u[0].hypot(u[1]), v[0].hypot(v[1]));
            let cross = u[0] * v[1] - u[1] * v[0];
            len_u * len_v > 0.0
                && cross.abs() <= COLLINEAR_SINE_TOLERANCE * len_u * len_v
                && u[0] * v[0] + u[1] * v[1] > 0.0
        });
        match removable {
            Some(i) => {
                ring.remove(i);
            }
            None => return ring,
        }
    }
}

/// Outline of the union of counter-clockwise polygons that tile a region
/// without overlapping, if that region is a single simple polygon.
///
/// Vertices are merged within [`COINCIDENT_TOLERANCE`], edges are split at
/// T-junctions, edges shared by two pieces (present in opposite directions)
/// cancel, and the remaining directed edges must form exactly one simple
/// closed loop.
fn merge_pieces_into_outline(pieces: &[Vec<[f64; 2]>]) -> Option<Vec<[f64; 2]>> {
    let tol = COINCIDENT_TOLERANCE;
    // Canonical vertex table and vertex-index rings.
    let mut table: Vec<[f64; 2]> = Vec::new();
    let mut rings: Vec<Vec<usize>> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let mut ids: Vec<usize> = Vec::with_capacity(piece.len());
        for &p in piece {
            let id = match table.iter().position(|&q| distance_2d(p, q) <= tol) {
                Some(id) => id,
                None => {
                    table.push(p);
                    table.len() - 1
                }
            };
            if ids.last() != Some(&id) {
                ids.push(id);
            }
        }
        while ids.len() > 1 && ids.first() == ids.last() {
            ids.pop();
        }
        if ids.len() >= 3 {
            rings.push(ids);
        }
    }
    // Directed boundary edges with T-junction splitting and cancellation.
    let mut directed: BTreeMap<(usize, usize), u32> = BTreeMap::new();
    for ring in &rings {
        for (k, &u) in ring.iter().enumerate() {
            let v = ring[(k + 1) % ring.len()];
            let mut on_edge: Vec<(f64, usize)> = table
                .iter()
                .enumerate()
                .filter(|&(id, _)| id != u && id != v)
                .filter_map(|(id, &p)| {
                    interior_segment_parameter(table[u], table[v], p, tol).map(|t| (t, id))
                })
                .collect();
            on_edge.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut from = u;
            for id in on_edge.iter().map(|&(_, id)| id).chain(std::iter::once(v)) {
                add_directed_edge(&mut directed, from, id);
                from = id;
            }
        }
    }
    // The boundary must be a permutation: every vertex has at most one
    // outgoing edge (a second one is a pinch point or an overlap).
    let mut next_vertex: BTreeMap<usize, usize> = BTreeMap::new();
    for (&(from, to), &count) in &directed {
        if count != 1 || next_vertex.insert(from, to).is_some() {
            return None;
        }
    }
    let start = *next_vertex.keys().next()?;
    let mut loop_ids = vec![start];
    let mut current = *next_vertex.get(&start)?;
    while current != start {
        loop_ids.push(current);
        if loop_ids.len() > next_vertex.len() {
            return None;
        }
        current = *next_vertex.get(&current)?;
    }
    if loop_ids.len() != next_vertex.len() {
        // More than one boundary loop: several regions or a hole.
        return None;
    }
    let mut outline: Vec<[f64; 2]> = loop_ids.iter().map(|&id| table[id]).collect();
    if polygon2d_signed_area(&outline) < 0.0 {
        outline.reverse();
    }
    let outline = remove_collinear_vertices(outline);
    (outline.len() >= 3).then_some(outline)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_ops::functions::polygon2d_is_convex;

    /// Convex hull by Andrew's monotone chain (counter-clockwise, no
    /// collinear points).
    fn monotone_chain_hull(points: &[[f64; 2]]) -> Vec<[f64; 2]> {
        let mut pts = points.to_vec();
        pts.sort_by(|p, q| p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1])));
        pts.dedup();
        if pts.len() < 3 {
            return pts;
        }
        let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| {
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        };
        let mut hull: Vec<[f64; 2]> = Vec::new();
        for &p in &pts {
            while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
                hull.pop();
            }
            hull.push(p);
        }
        let lower_len = hull.len() + 1;
        for &p in pts.iter().rev().skip(1) {
            while hull.len() >= lower_len
                && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
        hull
    }

    /// Nine pseudo-random points around `(cx, cy)` and their convex hull.
    fn random_convex_polygon(next: &mut impl FnMut() -> f64, cx: f64, cy: f64) -> Vec<[f64; 2]> {
        let pts: Vec<[f64; 2]> = (0..9)
            .map(|_| [cx + 2.0 * next() - 1.0, cy + 2.0 * next() - 1.0])
            .collect();
        monotone_chain_hull(&pts)
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<[f64; 2]> {
        vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    }

    fn total_area(pieces: &[Vec<[f64; 2]>]) -> f64 {
        pieces.iter().map(|p| polygon2d_area(p)).sum()
    }

    #[test]
    fn test_pieces_overlapping_squares_have_area_two() {
        // a = [0,2]², b = [1,3]×[0,2]: a \ b = [0,1]×[0,2], area 2.  b's
        // bottom and top edges are collinear with a's, which produces
        // zero-width slivers that must be dropped.
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 0.0, 3.0, 2.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert_eq!(pieces.len(), 1, "pieces: {pieces:?}");
        assert!((total_area(&pieces) - 2.0).abs() < 1e-12);
        for p in &pieces[0] {
            assert!(p[0] >= -1e-12 && p[0] <= 1.0 + 1e-12, "x = {}", p[0]);
            assert!(p[1] >= -1e-12 && p[1] <= 2.0 + 1e-12, "y = {}", p[1]);
        }
        let outline = polygon2d_difference(&a, &b).expect("single rectangle");
        assert_eq!(outline.len(), 4);
        assert!((polygon2d_signed_area(&outline) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_pieces_disjoint_returns_subject_unchanged() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(5.0, 0.0, 6.0, 1.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert_eq!(pieces.len(), 1);
        assert!((total_area(&pieces) - 4.0).abs() < 1e-12);
        let outline = polygon2d_difference(&a, &b).expect("a is returned");
        assert_eq!(outline.len(), 4);
        for (got, want) in outline.iter().zip(a.iter()) {
            assert!(distance_2d(*got, *want) < 1e-12, "{got:?} vs {want:?}");
        }
    }

    #[test]
    fn test_pieces_touching_returns_subject_unchanged() {
        let a = rect(0.0, 0.0, 1.0, 1.0);
        let b = rect(1.0, 0.0, 2.0, 1.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert_eq!(pieces.len(), 1);
        assert!((total_area(&pieces) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pieces_contained_subject_is_empty() {
        let a = rect(1.0, 1.0, 2.0, 2.0);
        let b = rect(0.0, 0.0, 3.0, 3.0);
        assert!(polygon2d_difference_pieces(&a, &b).is_empty());
        assert!(polygon2d_difference(&a, &b).expect("empty").is_empty());
        // Identical polygons.
        assert!(polygon2d_difference_pieces(&a, &a).is_empty());
    }

    #[test]
    fn test_corner_overlap_gives_l_shaped_outline() {
        // a = [0,4]², b = [2,6]²: a \ b is an L-shape of area 12.
        let a = rect(0.0, 0.0, 4.0, 4.0);
        let b = rect(2.0, 2.0, 6.0, 6.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert_eq!(pieces.len(), 2, "pieces: {pieces:?}");
        assert!((total_area(&pieces) - 12.0).abs() < 1e-12);
        for piece in &pieces {
            assert!(polygon2d_is_convex(piece));
        }
        let outline = polygon2d_difference(&a, &b).expect("L-shape is one polygon");
        assert_eq!(outline.len(), 6, "outline: {outline:?}");
        assert!((polygon2d_signed_area(&outline) - 12.0).abs() < 1e-12);
        assert!(!polygon2d_is_convex(&outline));
        let expected = [
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 2.0],
            [2.0, 2.0],
            [2.0, 4.0],
            [0.0, 4.0],
        ];
        for corner in expected {
            assert!(
                outline.iter().any(|&p| distance_2d(p, corner) < 1e-12),
                "corner {corner:?} missing from {outline:?}"
            );
        }
    }

    #[test]
    fn test_strictly_inside_leaves_a_ring_reported_as_error() {
        let a = rect(0.0, 0.0, 4.0, 4.0);
        let b = rect(1.0, 1.0, 3.0, 3.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert!((total_area(&pieces) - 12.0).abs() < 1e-12);
        assert!(
            pieces.len() >= 4,
            "a ring needs at least four convex pieces"
        );
        let err = polygon2d_difference(&a, &b).expect_err("a ring is not a polygon");
        assert!(
            err.to_string().contains("polygon2d_difference_pieces"),
            "{err}"
        );
    }

    #[test]
    fn test_strip_splits_subject_into_two_pieces() {
        let a = rect(0.0, 0.0, 3.0, 1.0);
        let b = rect(1.0, -1.0, 2.0, 2.0);
        let pieces = polygon2d_difference_pieces(&a, &b);
        assert_eq!(pieces.len(), 2, "pieces: {pieces:?}");
        assert!((total_area(&pieces) - 2.0).abs() < 1e-12);
        assert!(polygon2d_difference(&a, &b).is_err());
    }

    #[test]
    fn test_outline_function_rejects_non_convex_operands() {
        let l_shape = vec![
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 2.0],
            [2.0, 2.0],
            [2.0, 4.0],
            [0.0, 4.0],
        ];
        let square = rect(1.0, 1.0, 3.0, 3.0);
        for (a, b, which) in [(&l_shape, &square, "a"), (&square, &l_shape, "b")] {
            let err = polygon2d_difference(a, b).expect_err("operand is not convex");
            assert!(matches!(err, Error::Unsupported { .. }), "{err}");
            assert!(
                err.to_string()
                    .contains(&format!("polygon {which} is not convex")),
                "{err}"
            );
        }
        // Convex operands (either winding) are still accepted.
        let reversed: Vec<[f64; 2]> = square.iter().rev().copied().collect();
        assert!(polygon2d_difference(&reversed, &rect(5.0, 5.0, 6.0, 6.0)).is_ok());
    }

    #[test]
    fn test_clockwise_inputs_give_the_same_result() {
        let a = rect(0.0, 0.0, 4.0, 4.0);
        let b = rect(2.0, 2.0, 6.0, 6.0);
        let a_cw: Vec<[f64; 2]> = a.iter().rev().copied().collect();
        let b_cw: Vec<[f64; 2]> = b.iter().rev().copied().collect();
        let pieces = polygon2d_difference_pieces(&a_cw, &b_cw);
        assert!((total_area(&pieces) - 12.0).abs() < 1e-12);
        for piece in &pieces {
            assert!(polygon2d_signed_area(piece) > 0.0, "pieces are CCW");
        }
        let outline = polygon2d_difference(&a_cw, &b_cw).expect("L-shape");
        assert_eq!(outline.len(), 6);
        assert!(polygon2d_signed_area(&outline) > 0.0);
    }

    #[test]
    fn test_result_is_independent_of_scale_and_position() {
        for &(scale, shift) in &[(1.0e-7, [3.0e-7, -2.0e-7]), (1.0e6, [5.0e6, 4.0e6])] {
            let place = |poly: Vec<[f64; 2]>| -> Vec<[f64; 2]> {
                poly.into_iter()
                    .map(|p| [p[0] * scale + shift[0], p[1] * scale + shift[1]])
                    .collect()
            };
            let a = place(rect(0.0, 0.0, 4.0, 4.0));
            let b = place(rect(2.0, 2.0, 6.0, 6.0));
            let pieces = polygon2d_difference_pieces(&a, &b);
            let got = total_area(&pieces);
            let want = 12.0 * scale * scale;
            assert!(
                (got - want).abs() <= 1e-9 * want,
                "scale {scale}: {got} vs {want}"
            );
            let outline = polygon2d_difference(&a, &b).expect("L-shape");
            assert_eq!(outline.len(), 6, "scale {scale}");
        }
    }

    #[test]
    fn test_non_finite_and_degenerate_inputs() {
        let a = rect(0.0, 0.0, 1.0, 1.0);
        let nan = vec![[f64::NAN, 0.0], [1.0, 0.0], [0.0, 1.0]];
        assert!(polygon2d_difference_pieces(&a, &nan).is_empty());
        assert!(polygon2d_difference_pieces(&nan, &a).is_empty());
        assert!(polygon2d_difference(&nan, &a).expect("empty").is_empty());
        // Fewer than three vertices: `b` removes nothing, `a` is no polygon.
        let segment = vec![[0.0, 0.0], [1.0, 1.0]];
        let kept = polygon2d_difference_pieces(&a, &segment);
        assert_eq!(kept.len(), 1);
        assert!((total_area(&kept) - 1.0).abs() < 1e-12);
        assert!(polygon2d_difference_pieces(&segment, &a).is_empty());
        assert!(polygon2d_difference_pieces(&[], &a).is_empty());
        // Collinear cutter has zero area and removes nothing.
        let flat = vec![[0.2, 0.5], [0.5, 0.5], [0.8, 0.5]];
        let kept = polygon2d_difference_pieces(&a, &flat);
        assert!((total_area(&kept) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_random_convex_polygons_partition_the_difference() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut checked_overlaps = 0;
        for case in 0..80 {
            let a = random_convex_polygon(&mut next, 0.0, 0.0);
            let (shift_x, shift_y) = (1.5 * next() - 0.75, 1.5 * next() - 0.75);
            let b = random_convex_polygon(&mut next, shift_x, shift_y);
            if a.len() < 3 || b.len() < 3 {
                continue;
            }
            let pieces = polygon2d_difference_pieces(&a, &b);
            let area_a = polygon2d_area(&a);
            let overlap = polygon2d_area(&clip_polygon_2d(
                &counter_clockwise(&a),
                &counter_clockwise(&b),
            ));
            // The pieces tile exactly a \ b.
            assert!(
                (total_area(&pieces) - (area_a - overlap)).abs() < 1e-9,
                "case {case}: pieces {} vs a \\ b {}",
                total_area(&pieces),
                area_a - overlap
            );
            if overlap > 1e-9 {
                checked_overlaps += 1;
            }
            for (i, piece) in pieces.iter().enumerate() {
                assert!(polygon2d_is_convex(piece), "case {case}: piece {i}");
                // No piece overlaps b ...
                let with_b = polygon2d_area(&clip_polygon_2d(piece, &counter_clockwise(&b)));
                assert!(
                    with_b < 1e-9,
                    "case {case}: piece {i} overlaps b by {with_b}"
                );
                // ... and pieces do not overlap each other.
                for other in &pieces[i + 1..] {
                    let shared = polygon2d_area(&clip_polygon_2d(piece, other));
                    assert!(shared < 1e-9, "case {case}: pieces overlap by {shared}");
                }
            }
        }
        assert!(
            checked_overlaps > 20,
            "only {checked_overlaps} overlapping cases"
        );
    }
}
