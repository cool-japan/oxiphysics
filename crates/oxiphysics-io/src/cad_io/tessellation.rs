// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tessellation of BREP faces into triangles.
//!
//! Every face is triangulated in the plane of its outer boundary loop:
//!
//! 1. each edge loop is turned into a closed polygon of mesh vertices.  The
//!    stored edges are *undirected*: their direction is recovered from the
//!    connectivity of the loop, each edge being walked from the vertex where
//!    the previous edge ended;
//! 2. curved edges ([`CurveType::CircularArc`]) contribute interior sample
//!    points, which are created once per edge and shared by the two faces that
//!    meet there, so the mesh stays closed; straight edges are not subdivided;
//! 3. the polygons are projected into the face plane and triangulated by ear
//!    clipping, holes being joined to the outer polygon by bridges;
//! 4. the triangles are wound counter-clockwise about the face normal: the
//!    declared outward normal for a planar face, the normal given by the
//!    winding of the outer loop for other surface types.
//!
//! Curved surfaces are approximated by the flat facets spanned by their
//! boundary; the BREP model stores no surface geometry to refine them with.

use std::collections::HashMap;

use super::functions::{cross3, len3, normalize3, sub3};
use super::types::{BrepEdge, BrepFace, BrepSolid, CurveType, SurfaceType, TriangleMesh};

/// A point in the face plane.
type Point2 = [f64; 2];

/// Fraction of the extent of the loop concerned (its bounding box, for an
/// area) below which a polygon area or a ray-casting coordinate counts as zero.
const EPS: f64 = 1e-12;

/// Relative tolerance of the local orientation tests (corner turns and
/// point-in-triangle), measured against the lengths of the edges involved so
/// that very fine or very thin polygons keep working.
const REL_EPS: f64 = 1e-12;

/// Smallest `|cos|` between the declared normal of a planar face and the
/// normal of its boundary loop for the two to count as consistent.
const MIN_NORMAL_ALIGNMENT: f64 = 1e-6;

/// Relative tolerance of the area check on a finished triangulation.
const AREA_TOLERANCE: f64 = 1e-8;

/// Largest `segments_per_edge` accepted: every circular arc edge yields that
/// many vertices, so an absurd value would exhaust the memory instead of
/// producing a useful mesh.
const MAX_SEGMENTS_PER_EDGE: usize = 1 << 16;

/// Tessellate every face of `solid`; see the module documentation.
///
/// `segments_per_edge` is the number of straight pieces a circular arc edge is
/// divided into; 0 is treated as 1 and values above [`MAX_SEGMENTS_PER_EDGE`]
/// are rejected.
pub(super) fn tessellate(
    solid: &BrepSolid,
    segments_per_edge: usize,
) -> Result<TriangleMesh, String> {
    if segments_per_edge > MAX_SEGMENTS_PER_EDGE {
        return Err(format!(
            "segments_per_edge is {segments_per_edge}, more than the supported {MAX_SEGMENTS_PER_EDGE}"
        ));
    }
    solid.validate()?;
    let mut mesh = TriangleMesh::new();
    let vertex_ids: Vec<usize> = solid
        .vertices
        .iter()
        .map(|v| mesh.add_vertex(v.position))
        .collect();
    let mut tessellator = Tessellator {
        solid,
        segments: segments_per_edge.max(1),
        vertex_ids,
        samples: HashMap::new(),
        mesh,
    };
    for face in &solid.faces {
        tessellator.face(face)?;
    }
    let mut mesh = tessellator.mesh;
    mesh.compute_normals();
    Ok(mesh)
}

/// One edge of a loop, walked in a definite direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OrientedEdge {
    index: usize,
    /// `true` when walked from `start_vertex` to `end_vertex`.
    forward: bool,
}

/// The edge `index` of `solid`, or an error naming the loop that refers to it.
fn edge_of(
    solid: &BrepSolid,
    face_id: usize,
    loop_index: usize,
    index: usize,
) -> Result<&BrepEdge, String> {
    solid.edges.get(index).ok_or_else(|| {
        format!("face {face_id}: loop {loop_index} refers to the missing edge {index}")
    })
}

/// Decide in which direction each edge of a loop is walked.
///
/// The first edge is oriented by the edge that follows it (its end vertex must
/// touch the second edge, otherwise it is walked backwards); every other edge
/// is walked away from the vertex where the previous edge ended.
fn orient_loop(
    solid: &BrepSolid,
    face_id: usize,
    loop_index: usize,
    edge_loop: &[usize],
) -> Result<Vec<OrientedEdge>, String> {
    let (&first_index, rest) = edge_loop
        .split_first()
        .ok_or_else(|| format!("face {face_id}: loop {loop_index} has no edges"))?;
    let first = edge_of(solid, face_id, loop_index, first_index)?;

    let first_forward = match rest.first() {
        None => true,
        Some(&second_index) => {
            let second = edge_of(solid, face_id, loop_index, second_index)?;
            let touches_second =
                |vertex: usize| vertex == second.start_vertex || vertex == second.end_vertex;
            if touches_second(first.end_vertex) {
                true
            } else if touches_second(first.start_vertex) {
                false
            } else {
                return Err(format!(
                    "face {face_id}: the edges of loop {loop_index} are not connected: edge {first_index} does not touch edge {second_index}"
                ));
            }
        }
    };

    let (loop_start, mut tail) = if first_forward {
        (first.start_vertex, first.end_vertex)
    } else {
        (first.end_vertex, first.start_vertex)
    };
    let mut oriented = Vec::with_capacity(edge_loop.len());
    oriented.push(OrientedEdge {
        index: first_index,
        forward: first_forward,
    });
    for &index in rest {
        let next = edge_of(solid, face_id, loop_index, index)?;
        if next.start_vertex == tail {
            oriented.push(OrientedEdge {
                index,
                forward: true,
            });
            tail = next.end_vertex;
        } else if next.end_vertex == tail {
            oriented.push(OrientedEdge {
                index,
                forward: false,
            });
            tail = next.start_vertex;
        } else {
            return Err(format!(
                "face {face_id}: the edges of loop {loop_index} are not connected: edge {index} does not touch vertex {tail}"
            ));
        }
    }
    if tail != loop_start {
        return Err(format!(
            "face {face_id}: loop {loop_index} is not closed: it ends at vertex {tail} but starts at vertex {loop_start}"
        ));
    }
    Ok(oriented)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Vector area (twice the area times the normal) of a closed polygon, from
/// Newell's formula relative to the first vertex.
fn newell_normal(points: &[[f64; 3]]) -> [f64; 3] {
    let Some(&origin) = points.first() else {
        return [0.0; 3];
    };
    let mut sum = [0.0; 3];
    for pair in points.windows(2) {
        let c = cross3(sub3(pair[0], origin), sub3(pair[1], origin));
        sum = [sum[0] + c[0], sum[1] + c[1], sum[2] + c[2]];
    }
    sum
}

/// The two largest extents of the bounding box of `points`, longest first.
fn two_largest_extents3(points: &[[f64; 3]]) -> [f64; 2] {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in points {
        for ((lo, hi), &value) in min.iter_mut().zip(max.iter_mut()).zip(p) {
            *lo = lo.min(value);
            *hi = hi.max(value);
        }
    }
    let mut extents: Vec<f64> = min.iter().zip(&max).map(|(lo, hi)| hi - lo).collect();
    extents.sort_by(|a, b| b.total_cmp(a));
    [
        extents.first().copied().unwrap_or(0.0),
        extents.get(1).copied().unwrap_or(0.0),
    ]
}

/// Unit normal of the plane in which `face` is triangulated, the triangles
/// facing along it.
///
/// For a planar face with a usable declared normal this is the normal of the
/// boundary loop, oriented to agree with the declared one; for every other
/// face it is the normal given by the winding of the outer loop.
fn face_normal(face: &BrepFace, outer: &[[f64; 3]]) -> Result<[f64; 3], String> {
    let vector_area = newell_normal(outer);
    let magnitude = len3(vector_area);
    // The area of a planar polygon is at most the product of the two largest
    // extents of its bounding box; relative to it, a sliver is not "empty".
    let [longest, middle] = two_largest_extents3(outer);
    if !magnitude.is_finite() || !longest.is_finite() || magnitude <= EPS * longest * middle {
        return Err(format!("face {}: the outer loop encloses no area", face.id));
    }
    let loop_normal = [
        vector_area[0] / magnitude,
        vector_area[1] / magnitude,
        vector_area[2] / magnitude,
    ];
    if face.surface_type != SurfaceType::Plane {
        return Ok(loop_normal);
    }
    let declared = normalize3(face.normal);
    if len3(declared) < 0.5 {
        return Ok(loop_normal);
    }
    let alignment = dot3(loop_normal, declared);
    if alignment.abs() < MIN_NORMAL_ALIGNMENT {
        return Err(format!(
            "face {}: the declared normal is perpendicular to the plane of the outer loop",
            face.id
        ));
    }
    Ok(if alignment > 0.0 {
        loop_normal
    } else {
        [-loop_normal[0], -loop_normal[1], -loop_normal[2]]
    })
}

/// Orthonormal basis `(u, v)` of the plane with unit normal `n`, such that
/// `u x v = n`.
fn plane_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let helper = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize3(cross3(n, helper));
    let v = cross3(n, u);
    (u, v)
}

/// The mesh under construction and the state shared by its faces.
struct Tessellator<'a> {
    solid: &'a BrepSolid,
    /// Pieces a circular arc edge is divided into (at least 1).
    segments: usize,
    /// Mesh vertex of every BREP vertex.
    vertex_ids: Vec<usize>,
    /// Interior sample vertices of curved edges, by (edge, sample number).
    samples: HashMap<(usize, usize), usize>,
    mesh: TriangleMesh,
}

impl Tessellator<'_> {
    /// The closed polygon of mesh vertices along one edge loop.
    fn loop_polygon(
        &mut self,
        face_id: usize,
        loop_index: usize,
        edge_loop: &[usize],
    ) -> Result<Vec<usize>, String> {
        let mut polygon = Vec::with_capacity(edge_loop.len());
        for step in orient_loop(self.solid, face_id, loop_index, edge_loop)? {
            let edge = self
                .solid
                .edges
                .get(step.index)
                .ok_or_else(|| format!("face {face_id}: missing edge {}", step.index))?;
            let start_vertex = if step.forward {
                edge.start_vertex
            } else {
                edge.end_vertex
            };
            let start = self
                .vertex_ids
                .get(start_vertex)
                .copied()
                .ok_or_else(|| format!("face {face_id}: missing vertex {start_vertex}"))?;
            polygon.push(start);
            if self.segments > 1 && edge.curve_type == CurveType::CircularArc {
                for piece in 1..self.segments {
                    // The sample number counts along the edge itself, so both
                    // faces that use the edge find the same vertex.
                    let sample = if step.forward {
                        piece
                    } else {
                        self.segments - piece
                    };
                    let id = match self.samples.get(&(step.index, sample)) {
                        Some(&id) => id,
                        None => {
                            let t = sample as f64 / self.segments as f64;
                            let id = self.mesh.add_vertex(edge.evaluate(t, &self.solid.vertices));
                            self.samples.insert((step.index, sample), id);
                            id
                        }
                    };
                    polygon.push(id);
                }
            }
        }
        Ok(polygon)
    }

    /// Triangulate one face and add its triangles to the mesh.
    fn face(&mut self, face: &BrepFace) -> Result<(), String> {
        if face.edge_loops.is_empty() {
            return Ok(());
        }
        let mut ids: Vec<usize> = Vec::new();
        let mut loop_lengths: Vec<usize> = Vec::with_capacity(face.edge_loops.len());
        for (loop_index, edge_loop) in face.edge_loops.iter().enumerate() {
            let polygon = self.loop_polygon(face.id, loop_index, edge_loop)?;
            if polygon.len() < 3 {
                return Err(format!(
                    "face {}: loop {loop_index} has only {} vertices",
                    face.id,
                    polygon.len()
                ));
            }
            loop_lengths.push(polygon.len());
            ids.extend(polygon);
        }
        let positions: Vec<[f64; 3]> = ids
            .iter()
            .map(|&id| {
                self.mesh
                    .vertices
                    .get(id)
                    .copied()
                    .ok_or_else(|| format!("face {}: missing mesh vertex {id}", face.id))
            })
            .collect::<Result<_, _>>()?;
        let outer_len = loop_lengths.first().copied().unwrap_or(0);
        let outer = positions
            .get(..outer_len)
            .ok_or_else(|| format!("face {}: inconsistent loop lengths", face.id))?;
        let normal = face_normal(face, outer)?;
        let (u, v) = plane_basis(normal);
        let origin = positions.first().copied().unwrap_or([0.0; 3]);
        let projected: Vec<Point2> = positions
            .iter()
            .map(|&p| {
                let d = sub3(p, origin);
                [dot3(d, u), dot3(d, v)]
            })
            .collect();

        let triangles = triangulate(&projected, &loop_lengths)
            .map_err(|message| format!("face {}: {message}", face.id))?;
        for [a, b, c] in triangles {
            let (Some(&a), Some(&b), Some(&c)) = (ids.get(a), ids.get(b), ids.get(c)) else {
                return Err(format!(
                    "face {}: the triangulation refers to a vertex outside the face",
                    face.id
                ));
            };
            self.mesh.add_triangle(a, b, c);
        }
        Ok(())
    }
}

// ── 2D triangulation ───────────────────────────────────────────────────────

/// Twice the signed area of the triangle `a, b, c`: positive when it runs
/// counter-clockwise.
fn orient2d(a: Point2, b: Point2, c: Point2) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Twice the signed area of the polygon through `points[ring[i]]`.
fn signed_area2(points: &[Point2], ring: &[usize]) -> f64 {
    let mut sum = 0.0;
    for (i, &a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        sum += points[a][0] * points[b][1] - points[b][0] * points[a][1];
    }
    sum
}

/// Width and height of the bounding box of the loop through `points[ring[i]]`.
fn loop_extents(points: &[Point2], ring: &[usize]) -> (f64, f64) {
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    for &i in ring {
        for ((lo, hi), &value) in min.iter_mut().zip(max.iter_mut()).zip(&points[i]) {
            *lo = lo.min(value);
            *hi = hi.max(value);
        }
    }
    (max[0] - min[0], max[1] - min[1])
}

/// Whether the loop encloses (almost) no area compared with its bounding box:
/// all its points on one line, or a sliver thinner than [`EPS`] of its extent.
fn encloses_no_area(points: &[Point2], ring: &[usize]) -> bool {
    let (width, height) = loop_extents(points, ring);
    signed_area2(points, ring).abs() <= EPS * width * height
}

fn distance2d(a: Point2, b: Point2) -> f64 {
    (b[0] - a[0]).hypot(b[1] - a[1])
}

/// Whether the path `a -> b -> c` turns strictly left at `b`; a turn whose
/// sine is below [`REL_EPS`] counts as straight.
fn turns_left(a: Point2, b: Point2, c: Point2) -> bool {
    orient2d(a, b, c) > REL_EPS * distance2d(a, b) * distance2d(b, c)
}

/// Whether `a -> b -> c` is straight (or doubles back) at `b`, up to
/// [`REL_EPS`].
fn is_straight(a: Point2, b: Point2, c: Point2) -> bool {
    orient2d(a, b, c).abs() <= REL_EPS * distance2d(a, b) * distance2d(b, c)
}

/// Whether `q` lies in the closed triangle `a, b, c` (either orientation).
fn in_triangle(a: Point2, b: Point2, c: Point2, q: Point2) -> bool {
    let longest = distance2d(a, b).max(distance2d(b, c)).max(distance2d(c, a));
    let tolerance = REL_EPS * longest * longest;
    let d1 = orient2d(a, b, q);
    let d2 = orient2d(b, c, q);
    let d3 = orient2d(c, a, q);
    let has_negative = d1 < -tolerance || d2 < -tolerance || d3 < -tolerance;
    let has_positive = d1 > tolerance || d2 > tolerance || d3 > tolerance;
    !(has_negative && has_positive)
}

/// Whether the direction from `cur` towards `target` points into the interior
/// angle of a counter-clockwise polygon at `cur`, whose neighbours are `prev`
/// and `next`.
fn inside_wedge(prev: Point2, cur: Point2, next: Point2, target: Point2) -> bool {
    let turn = orient2d(prev, cur, next);
    let left_of_incoming = orient2d(prev, cur, target) >= 0.0;
    let left_of_outgoing = orient2d(cur, next, target) >= 0.0;
    if turn > 0.0 {
        left_of_incoming && left_of_outgoing
    } else if turn < 0.0 {
        left_of_incoming || left_of_outgoing
    } else {
        left_of_outgoing
    }
}

/// Position in `ring` of the occurrence of the vertex at `p_pos` (a vertex can
/// occur several times after earlier bridges) whose interior angle contains
/// the direction towards `target`.
fn visible_occurrence(points: &[Point2], ring: &[usize], p_pos: usize, target: Point2) -> usize {
    let n = ring.len();
    let p = points[ring[p_pos]];
    for pos in 0..n {
        if points[ring[pos]] != p {
            continue;
        }
        let prev = points[ring[(pos + n - 1) % n]];
        let next = points[ring[(pos + 1) % n]];
        if inside_wedge(prev, p, next, target) {
            return pos;
        }
    }
    p_pos
}

/// Join the clockwise `hole` to the counter-clockwise `ring` through a bridge
/// between the right-most hole vertex and a ring vertex it can see.
fn bridge_hole(points: &[Point2], ring: &mut Vec<usize>, hole: &[usize]) -> Result<(), String> {
    let (m_pos, m_idx) = hole
        .iter()
        .copied()
        .enumerate()
        .max_by(|x, y| {
            points[x.1][0]
                .total_cmp(&points[y.1][0])
                .then(points[y.1][1].total_cmp(&points[x.1][1]))
        })
        .ok_or_else(|| "a hole has no vertices".to_string())?;
    let m = points[m_idx];
    let n = ring.len();
    // Tolerance along x, relative to the extent of the polygon in x (which
    // can be far below the unit frame's scale for a tall, thin polygon).
    let x_extent = loop_extents(points, ring)
        .0
        .max(loop_extents(points, hole).0);
    let eps_x = EPS * x_extent;

    // The nearest ring edge crossed by the ray from `m` towards +x.
    let mut nearest: Option<(f64, usize)> = None;
    for e in 0..n {
        let a = points[ring[e]];
        let b = points[ring[(e + 1) % n]];
        if (a[1] - m[1]) * (b[1] - m[1]) > 0.0 || a[1] == b[1] {
            continue;
        }
        let t = (m[1] - a[1]) / (b[1] - a[1]);
        let x = a[0] + t * (b[0] - a[0]);
        if x < m[0] - eps_x {
            continue;
        }
        if nearest.is_none_or(|(best_x, _)| x < best_x) {
            nearest = Some((x, e));
        }
    }
    let (hit_x, edge) =
        nearest.ok_or_else(|| "a hole lies outside the outer boundary".to_string())?;

    // Candidate: the end of the hit edge with the larger x.
    let (a_pos, b_pos) = (edge, (edge + 1) % n);
    let mut p_pos = if points[ring[a_pos]][0] >= points[ring[b_pos]][0] {
        a_pos
    } else {
        b_pos
    };
    let hit = [hit_x, m[1]];
    let p = points[ring[p_pos]];
    if (p[0] - hit_x).abs() > eps_x {
        // Vertices inside the triangle `m, hit, p` can block the view of `p`;
        // the one making the smallest angle with the ray is visible.
        let mut best_tan = f64::INFINITY;
        let mut best_dist = f64::INFINITY;
        let mut best_pos = None;
        for (pos, &idx) in ring.iter().enumerate() {
            let q = points[idx];
            if q == p || q[0] <= m[0] || !in_triangle(m, hit, p, q) {
                continue;
            }
            let dx = q[0] - m[0];
            let dy = q[1] - m[1];
            let tan = dy.abs() / dx;
            let dist = dx * dx + dy * dy;
            if tan < best_tan || (tan == best_tan && dist < best_dist) {
                best_tan = tan;
                best_dist = dist;
                best_pos = Some(pos);
            }
        }
        if let Some(pos) = best_pos {
            p_pos = pos;
        }
    }
    let p_pos = visible_occurrence(points, ring, p_pos, m);

    let mut merged = Vec::with_capacity(ring.len() + hole.len() + 2);
    merged.extend_from_slice(&ring[..=p_pos]);
    for k in 0..hole.len() {
        merged.push(hole[(m_pos + k) % hole.len()]);
    }
    merged.push(m_idx);
    merged.push(ring[p_pos]);
    merged.extend_from_slice(&ring[p_pos + 1..]);
    *ring = merged;
    Ok(())
}

/// Whether the ring vertex at `c` is the tip of an ear: a strictly convex
/// corner whose triangle contains no other vertex of the ring.
fn is_ear(points: &[Point2], ring: &[usize], c: usize) -> bool {
    let n = ring.len();
    let prev = (c + n - 1) % n;
    let next = (c + 1) % n;
    let (a, b, d) = (points[ring[prev]], points[ring[c]], points[ring[next]]);
    if !turns_left(a, b, d) {
        return false;
    }
    ring.iter().enumerate().all(|(k, &idx)| {
        if k == prev || k == c || k == next {
            return true;
        }
        let q = points[idx];
        q == a || q == b || q == d || !in_triangle(a, b, d, q)
    })
}

/// Cut the weakly simple polygon `ring` into triangles by ear clipping.
fn clip_ears(points: &[Point2], mut ring: Vec<usize>) -> Result<Vec<[usize; 3]>, String> {
    let mut triangles = Vec::with_capacity(ring.len().saturating_sub(2));
    let mut cursor = 0usize;
    while ring.len() > 3 {
        let n = ring.len();
        let ear = (0..n)
            .map(|step| (cursor + step) % n)
            .find(|&c| is_ear(points, &ring, c));
        if let Some(c) = ear {
            let prev = (c + n - 1) % n;
            let next = (c + 1) % n;
            triangles.push([ring[prev], ring[c], ring[next]]);
            ring.remove(c);
            cursor = c.saturating_sub(1) % ring.len();
            continue;
        }
        // No ear: drop a vertex that adds no area (a collinear point or a
        // spike), which unblocks the corners next to it.
        let flat = (0..n).find(|&c| {
            let a = points[ring[(c + n - 1) % n]];
            let b = points[ring[c]];
            let d = points[ring[(c + 1) % n]];
            is_straight(a, b, d)
        });
        match flat {
            Some(c) => {
                ring.remove(c);
                cursor = c.saturating_sub(1) % ring.len();
            }
            None => {
                return Err(
                    "no ear found: the boundary probably intersects itself or has features thinner than the relative tolerance of its extent"
                        .to_string(),
                );
            }
        }
    }
    if let &[a, b, c] = ring.as_slice()
        && !is_straight(points[a], points[b], points[c])
    {
        triangles.push([a, b, c]);
    }
    Ok(triangles)
}

/// Triangulate a polygon with holes.
///
/// `points` holds the vertices of all loops one after the other;
/// `loop_lengths[0]` is the length of the outer loop and the other lengths
/// those of the holes.  The orientation of the loops does not matter.  The
/// result holds counter-clockwise triangles as triples of indices into
/// `points`.
///
/// # Errors
///
/// Fails when the loops do not describe a polygon with holes: a loop of zero
/// area, a hole outside the outer loop, a self-intersecting boundary, or a
/// triangulation whose area does not add up to the area of the polygon.
fn triangulate(points: &[Point2], loop_lengths: &[usize]) -> Result<Vec<[usize; 3]>, String> {
    // Work in a unit-sized frame so that the tolerances are relative.
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    for p in points {
        for ((lo, hi), &value) in min.iter_mut().zip(max.iter_mut()).zip(p) {
            *lo = lo.min(value);
            *hi = hi.max(value);
        }
    }
    let scale = (max[0] - min[0]).max(max[1] - min[1]);
    if !scale.is_finite() || scale <= 0.0 {
        return Err("the boundary has no extent".to_string());
    }
    let unit: Vec<Point2> = points
        .iter()
        .map(|p| [(p[0] - min[0]) / scale, (p[1] - min[1]) / scale])
        .collect();

    let mut loops: Vec<Vec<usize>> = Vec::with_capacity(loop_lengths.len());
    let mut start = 0usize;
    for &len in loop_lengths {
        let end = start + len;
        if end > unit.len() {
            return Err("loop lengths exceed the number of points".to_string());
        }
        loops.push((start..end).collect());
        start = end;
    }
    if loops.is_empty() {
        return Err("there is no outer loop".to_string());
    }
    let mut outer = loops.remove(0);
    let outer_area2 = signed_area2(&unit, &outer);
    if encloses_no_area(&unit, &outer) {
        return Err("the outer loop encloses no area".to_string());
    }
    if outer_area2 < 0.0 {
        outer.reverse();
    }
    let mut expected_area2 = outer_area2.abs();

    // Holes run clockwise; the right-most one is bridged first.
    let mut holes: Vec<(f64, Vec<usize>)> = Vec::with_capacity(loops.len());
    for mut hole in loops {
        let area2 = signed_area2(&unit, &hole);
        if encloses_no_area(&unit, &hole) {
            return Err("a hole encloses no area".to_string());
        }
        if area2 > 0.0 {
            hole.reverse();
        }
        expected_area2 -= area2.abs();
        let right_most = hole
            .iter()
            .map(|&i| unit[i][0])
            .fold(f64::NEG_INFINITY, f64::max);
        holes.push((right_most, hole));
    }
    holes.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut ring = outer;
    for (_, hole) in &holes {
        bridge_hole(&unit, &mut ring, hole)?;
    }

    let triangles = clip_ears(&unit, ring)?;
    let covered_area2: f64 = triangles
        .iter()
        .map(|&[a, b, c]| orient2d(unit[a], unit[b], unit[c]).abs())
        .sum();
    if (covered_area2 - expected_area2).abs() > AREA_TOLERANCE * outer_area2.abs() {
        return Err(format!(
            "the triangulation covers {covered_area2} instead of {expected_area2} (twice the area in the unit frame): the boundary is probably invalid"
        ));
    }
    Ok(triangles)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(points: &[Point2], triangles: &[[usize; 3]]) -> f64 {
        triangles
            .iter()
            .map(|&[a, b, c]| 0.5 * orient2d(points[a], points[b], points[c]))
            .sum()
    }

    fn all_counter_clockwise(points: &[Point2], triangles: &[[usize; 3]]) -> bool {
        triangles
            .iter()
            .all(|&[a, b, c]| orient2d(points[a], points[b], points[c]) > 0.0)
    }

    fn square(x0: f64, y0: f64, side: f64) -> Vec<Point2> {
        vec![
            [x0, y0],
            [x0 + side, y0],
            [x0 + side, y0 + side],
            [x0, y0 + side],
        ]
    }

    #[test]
    fn test_square_gives_two_counter_clockwise_triangles() {
        let points = square(0.0, 0.0, 2.0);
        let triangles = triangulate(&points, &[4]).unwrap();
        assert_eq!(triangles.len(), 2);
        assert!(all_counter_clockwise(&points, &triangles));
        assert!((area(&points, &triangles) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_clockwise_outer_loop_is_accepted() {
        let mut points = square(0.0, 0.0, 2.0);
        points.reverse();
        let triangles = triangulate(&points, &[4]).unwrap();
        assert_eq!(triangles.len(), 2);
        assert!(all_counter_clockwise(&points, &triangles));
    }

    #[test]
    fn test_l_shape_with_a_reflex_corner() {
        let points = vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let triangles = triangulate(&points, &[6]).unwrap();
        assert_eq!(triangles.len(), 4);
        assert!(all_counter_clockwise(&points, &triangles));
        assert!((area(&points, &triangles) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_star_polygon() {
        let n = 24;
        let points: Vec<Point2> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * i as f64 / n as f64;
                let radius = if i % 2 == 0 { 3.0 } else { 1.2 };
                [radius * angle.cos(), radius * angle.sin()]
            })
            .collect();
        let triangles = triangulate(&points, &[n]).unwrap();
        assert_eq!(triangles.len(), n - 2);
        assert!(all_counter_clockwise(&points, &triangles));
        let ring: Vec<usize> = (0..n).collect();
        let expected = 0.5 * signed_area2(&points, &ring);
        assert!((area(&points, &triangles) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_square_with_a_square_hole() {
        let mut points = square(0.0, 0.0, 4.0);
        points.extend(square(1.0, 1.0, 2.0));
        let triangles = triangulate(&points, &[4, 4]).unwrap();
        assert!(all_counter_clockwise(&points, &triangles));
        assert!((area(&points, &triangles) - 12.0).abs() < 1e-12);
        // 8 polygon vertices and the two ends of the bridge counted twice.
        assert_eq!(triangles.len(), 8);
        // No triangle may cover the hole: its centre must stay uncovered.
        let centre = [2.0, 2.0];
        for &[a, b, c] in &triangles {
            assert!(
                !(orient2d(points[a], points[b], centre) > 1e-9
                    && orient2d(points[b], points[c], centre) > 1e-9
                    && orient2d(points[c], points[a], centre) > 1e-9),
                "triangle {a} {b} {c} covers the hole"
            );
        }
    }

    #[test]
    fn test_hole_orientation_does_not_matter() {
        let mut points = square(0.0, 0.0, 4.0);
        let mut hole = square(1.0, 1.0, 2.0);
        hole.reverse();
        points.extend(hole);
        let triangles = triangulate(&points, &[4, 4]).unwrap();
        assert!((area(&points, &triangles) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn test_two_holes_and_a_concave_hole() {
        let mut points = square(0.0, 0.0, 10.0);
        points.extend(square(1.0, 1.0, 2.0));
        // A triangular hole and an L-shaped one.
        points.extend([[5.0, 1.0], [8.0, 1.0], [6.0, 4.0]]);
        points.extend([
            [2.0, 5.0],
            [6.0, 5.0],
            [6.0, 6.0],
            [3.0, 6.0],
            [3.0, 8.0],
            [2.0, 8.0],
        ]);
        let triangles = triangulate(&points, &[4, 4, 3, 6]).unwrap();
        // The square hole (4), the triangle (3 * 3 / 2) and the L (4 + 2).
        let expected = 100.0 - 4.0 - 4.5 - (4.0 + 2.0);
        assert!(all_counter_clockwise(&points, &triangles));
        assert!(
            (area(&points, &triangles) - expected).abs() < 1e-9,
            "area {} vs {expected}",
            area(&points, &triangles)
        );
    }

    #[test]
    fn test_holes_stacked_in_the_same_column() {
        let mut points = square(0.0, 0.0, 10.0);
        points.extend(square(4.0, 1.0, 2.0));
        points.extend(square(4.0, 4.0, 2.0));
        points.extend(square(4.0, 7.0, 2.0));
        let triangles = triangulate(&points, &[4, 4, 4, 4]).unwrap();
        assert!((area(&points, &triangles) - 88.0).abs() < 1e-9);
    }

    #[test]
    fn test_invalid_boundaries_are_errors() {
        // A hole outside the outer loop.
        let mut points = square(0.0, 0.0, 2.0);
        points.extend(square(10.0, 10.0, 1.0));
        assert!(triangulate(&points, &[4, 4]).is_err());
        // A collinear outer loop.
        let line = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        assert!(triangulate(&line, &[3]).unwrap_err().contains("no area"));
        // A bow tie.
        let bow_tie = vec![[0.0, 0.0], [2.0, 2.0], [2.0, 0.0], [0.0, 2.0]];
        assert!(triangulate(&bow_tie, &[4]).is_err());
        // All points identical.
        let same = vec![[1.0, 1.0]; 4];
        assert!(triangulate(&same, &[4]).unwrap_err().contains("no extent"));
    }

    #[test]
    fn test_collinear_boundary_vertices_leave_no_gaps() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let triangles = triangulate(&points, &[5]).unwrap();
        assert!((area(&points, &triangles) - 4.0).abs() < 1e-12);
        assert!(all_counter_clockwise(&points, &triangles));
    }

    #[test]
    fn test_very_thin_polygon_is_triangulated() {
        // Aspect ratio 1e9: in the unit frame the strip is 1e-9 thick and the
        // notch in its top edge turns by about 2e-9 radians.  That is far below
        // an absolute tolerance but above the relative one.
        let points = vec![
            [0.0, 0.0],
            [1e6, 0.0],
            [1e6, 1e-3],
            [5e5, 5e-4],
            [0.0, 1e-3],
        ];
        let triangles = triangulate(&points, &[5]).unwrap();
        assert_eq!(triangles.len(), 3);
        let expected = 1e6 * 1e-3 - 0.5 * 1e6 * 5e-4;
        assert!(
            (area(&points, &triangles) - expected).abs() < 1e-9 * expected,
            "area {} vs {expected}",
            area(&points, &triangles)
        );
        assert!(all_counter_clockwise(&points, &triangles));
    }

    #[test]
    fn test_polygon_thinner_than_the_tolerance_is_rejected() {
        // Aspect ratio 1e13: the notch turns by about 2e-13 radians, below the
        // relative tolerance of 1e-12, so it cannot be told from a straight
        // edge.  Such a polygon is reported, never silently mis-triangulated.
        let points = vec![
            [0.0, 0.0],
            [1e6, 0.0],
            [1e6, 1e-7],
            [5e5, 5e-8],
            [0.0, 1e-7],
        ];
        let message = triangulate(&points, &[5]).unwrap_err();
        assert!(message.contains("thinner"), "{message}");
    }

    #[test]
    fn test_finely_sampled_circle() {
        let n = 6000;
        let points: Vec<Point2> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * i as f64 / n as f64;
                [angle.cos(), angle.sin()]
            })
            .collect();
        let triangles = triangulate(&points, &[n]).unwrap();
        assert_eq!(triangles.len(), n - 2);
        let polygon = 0.5 * n as f64 * (std::f64::consts::TAU / n as f64).sin();
        assert!((area(&points, &triangles) - polygon).abs() < 1e-9);
    }

    #[test]
    fn test_tiny_and_huge_polygons_use_relative_tolerances() {
        for scale in [1e-6, 1.0, 1e6] {
            let mut points = square(0.0, 0.0, 4.0 * scale);
            points.extend(square(scale, scale, 2.0 * scale));
            let triangles = triangulate(&points, &[4, 4]).unwrap();
            let expected = 12.0 * scale * scale;
            assert!(
                (area(&points, &triangles) - expected).abs() < 1e-9 * expected,
                "scale {scale}"
            );
        }
    }
}

#[cfg(test)]
mod mesh_tests {
    use super::*;
    use crate::cad_io::{Assembly, AssemblyComponent, AssemblyTransform, LengthUnit, StlExporter};

    fn volume(mesh: &TriangleMesh) -> f64 {
        mesh.triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (
                    mesh.vertices[t[0]],
                    mesh.vertices[t[1]],
                    mesh.vertices[t[2]],
                );
                dot3(a, cross3(b, c)) / 6.0
            })
            .sum()
    }

    fn net_vector_area(mesh: &TriangleMesh) -> [f64; 3] {
        let mut sum = [0.0; 3];
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let n = cross3(sub3(b, a), sub3(c, a));
            sum = [
                sum[0] + 0.5 * n[0],
                sum[1] + 0.5 * n[1],
                sum[2] + 0.5 * n[2],
            ];
        }
        sum
    }

    /// Every directed edge occurs exactly once and so does its reverse: the
    /// mesh is closed and consistently oriented.
    fn assert_closed_and_consistently_oriented(mesh: &TriangleMesh) {
        let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &mesh.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *directed.entry((a, b)).or_insert(0) += 1;
            }
        }
        for (&(a, b), &count) in &directed {
            assert_eq!(count, 1, "directed edge {a}->{b} used {count} times");
            assert_eq!(
                directed.get(&(b, a)),
                Some(&1),
                "edge {a}->{b} has no opposite partner"
            );
        }
    }

    /// Every directed edge is matched by as many triangles running the other
    /// way: the surface has no boundary and is oriented, even if two faces
    /// happen to share diagonals.
    fn assert_balanced(mesh: &TriangleMesh) {
        let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &mesh.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *directed.entry((a, b)).or_insert(0) += 1;
            }
        }
        for (&(a, b), &count) in &directed {
            assert_eq!(
                directed.get(&(b, a)).copied().unwrap_or(0),
                count,
                "edge {a}->{b} is not matched by the opposite edge"
            );
        }
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() < tolerance,
            "{actual} vs {expected}"
        );
    }

    #[test]
    fn test_box_is_a_closed_outward_oriented_mesh() {
        for segments in [0, 1, 10] {
            let solid = BrepSolid::create_box("box", 1.0, 2.0, 3.0);
            let mesh = tessellate(&solid, segments).unwrap();
            assert_eq!(mesh.triangle_count(), 12, "segments {segments}");
            assert_eq!(mesh.vertex_count(), 8, "segments {segments}");
            assert_closed_and_consistently_oriented(&mesh);
            assert_close(volume(&mesh), 6.0, 1e-12);
            assert_close(mesh.surface_area(), 22.0, 1e-12);
            let net = net_vector_area(&mesh);
            assert!(len3(net) < 1e-12, "net vector area {net:?}");
        }
    }

    #[test]
    fn test_slender_boxes_are_closed_meshes_of_the_right_volume() {
        // A long thin slab and a needle: the faces have aspect ratios of
        // 1e7 and 1e8.
        for (sx, sy, sz) in [(1e4, 1e-3, 1.0), (1e-4, 1.0, 1e4)] {
            let solid = BrepSolid::create_box("slender", sx, sy, sz);
            let mesh = tessellate(&solid, 1).unwrap();
            assert_eq!(mesh.triangle_count(), 12, "box {sx} x {sy} x {sz}");
            assert_closed_and_consistently_oriented(&mesh);
            let expected = sx * sy * sz;
            assert_close(volume(&mesh), expected, 1e-9 * expected);
        }
    }

    #[test]
    fn test_box_volume_does_not_depend_on_where_the_box_is() {
        // A face lying in a plane through the origin contributes nothing to
        // the volume whichever way it faces, so move the box away.
        let mut solid = BrepSolid::create_box("box", 1.0, 2.0, 3.0);
        for vertex in &mut solid.vertices {
            vertex.position = [
                vertex.position[0] + 10.0,
                vertex.position[1] - 20.0,
                vertex.position[2] + 30.0,
            ];
        }
        let mesh = tessellate(&solid, 4).unwrap();
        assert_closed_and_consistently_oriented(&mesh);
        assert_close(volume(&mesh), 6.0, 1e-9);
    }

    #[test]
    fn test_every_box_face_triangle_faces_along_the_declared_normal() {
        let solid = BrepSolid::create_box("box", 1.0, 2.0, 3.0);
        let mesh = tessellate(&solid, 1).unwrap();
        let mut seen = vec![0usize; solid.faces.len()];
        for (triangle, normal) in mesh.triangles.iter().zip(&mesh.normals) {
            let centroid = {
                let (a, b, c) = (
                    mesh.vertices[triangle[0]],
                    mesh.vertices[triangle[1]],
                    mesh.vertices[triangle[2]],
                );
                [
                    (a[0] + b[0] + c[0]) / 3.0,
                    (a[1] + b[1] + c[1]) / 3.0,
                    (a[2] + b[2] + c[2]) / 3.0,
                ]
            };
            // The triangle belongs to the face whose plane contains it and
            // whose declared normal it must face.
            let face_index = solid
                .faces
                .iter()
                .position(|f| dot3(sub3(centroid, f.origin), f.normal).abs() < 1e-12)
                .expect("every triangle lies in a face plane");
            assert!(
                dot3(*normal, solid.faces[face_index].normal) > 0.999,
                "triangle {triangle:?} faces {normal:?}, face {face_index} wants {:?}",
                solid.faces[face_index].normal
            );
            seen[face_index] += 1;
        }
        assert!(seen.iter().all(|&n| n == 2), "{seen:?}");
    }

    #[test]
    fn test_box_face_loops_wind_counter_clockwise_about_their_normals() {
        let solid = BrepSolid::create_box("box", 1.0, 2.0, 3.0);
        for face in &solid.faces {
            let oriented = orient_loop(&solid, face.id, 0, &face.edge_loops[0]).unwrap();
            let positions: Vec<[f64; 3]> = oriented
                .iter()
                .map(|step| {
                    let edge = &solid.edges[step.index];
                    let vertex = if step.forward {
                        edge.start_vertex
                    } else {
                        edge.end_vertex
                    };
                    solid.vertices[vertex].position
                })
                .collect();
            let normal = normalize3(newell_normal(&positions));
            assert!(
                dot3(normal, face.normal) > 0.999,
                "face {}: loop normal {normal:?}, declared {:?}",
                face.id,
                face.normal
            );
        }
    }

    #[test]
    fn test_edges_are_oriented_by_connectivity() {
        let solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        // Face 2 is [e0, e9, e4, e8]: e0 and e9 run with the loop, e4 and e8
        // against it.
        let oriented = orient_loop(&solid, 2, 0, &solid.faces[2].edge_loops[0]).unwrap();
        let flags: Vec<bool> = oriented.iter().map(|s| s.forward).collect();
        assert_eq!(flags, vec![true, true, false, false]);
        // The same loop started from another edge and read backwards.
        let reversed: Vec<usize> = solid.faces[2].edge_loops[0].iter().rev().copied().collect();
        let oriented = orient_loop(&solid, 2, 0, &reversed).unwrap();
        let flags: Vec<bool> = oriented.iter().map(|s| s.forward).collect();
        assert_eq!(flags, vec![true, true, false, false]);
    }

    #[test]
    fn test_loop_order_and_direction_do_not_matter() {
        let reference = tessellate(&BrepSolid::create_box("box", 1.0, 2.0, 3.0), 1).unwrap();
        let mut solid = BrepSolid::create_box("box", 1.0, 2.0, 3.0);
        for (index, face) in solid.faces.iter_mut().enumerate() {
            let edge_loop = &mut face.edge_loops[0];
            if index % 2 == 0 {
                edge_loop.reverse();
            }
            edge_loop.rotate_left(index % 4);
        }
        let mesh = tessellate(&solid, 1).unwrap();
        assert_eq!(mesh.triangle_count(), 12);
        assert_closed_and_consistently_oriented(&mesh);
        assert_close(volume(&mesh), volume(&reference), 1e-12);
    }

    #[test]
    fn test_defective_loops_are_errors_with_the_face_named() {
        // Not connected: the edges of face 0 are replaced by unrelated ones.
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[0].edge_loops[0] = vec![0, 2, 4, 6];
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(
            message.contains("face 0") && message.contains("not connected"),
            "{message}"
        );

        // Connected but not closed.
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[0].edge_loops[0] = vec![3, 2, 1];
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(message.contains("not closed"), "{message}");

        // A reference to an edge that does not exist.
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[3].edge_loops[0][1] = 99;
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(message.contains("99"), "{message}");

        // A loop with fewer than three vertices.
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[1].edge_loops[0] = vec![4, 4];
        assert!(tessellate(&solid, 1).is_err());

        // An empty loop.
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[1].edge_loops[0].clear();
        assert!(tessellate(&solid, 1).unwrap_err().contains("no edges"));
    }

    #[test]
    fn test_absurd_segment_counts_are_rejected_not_allocated() {
        let solid = pillow();
        let message = tessellate(&solid, usize::MAX).unwrap_err();
        assert!(message.contains("segments_per_edge"), "{message}");
        assert!(tessellate(&solid, MAX_SEGMENTS_PER_EDGE + 1).is_err());
        assert!(tessellate(&BrepSolid::create_box("b", 1.0, 1.0, 1.0), 1 << 16).is_ok());
    }

    #[test]
    fn test_faces_without_loops_are_skipped() {
        let mut solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        solid.faces[5].edge_loops.clear();
        let mesh = tessellate(&solid, 1).unwrap();
        assert_eq!(mesh.triangle_count(), 10);
        let empty = tessellate(&BrepSolid::new("empty"), 8).unwrap();
        assert_eq!(empty.triangle_count(), 0);
        assert_eq!(empty.vertex_count(), 0);
    }

    /// A single planar face over the unit square with the given declared normal
    /// and the vertices in the given order.
    fn square_face(order: [usize; 4], normal: [f64; 3], surface_type: SurfaceType) -> BrepSolid {
        let mut solid = BrepSolid::new("face");
        let corners = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let ids: Vec<usize> = order
            .iter()
            .map(|&i| solid.add_vertex(corners[i]))
            .collect();
        let edges: Vec<usize> = (0..4)
            .map(|i| solid.add_line_edge(ids[i], ids[(i + 1) % 4]))
            .collect();
        solid.add_planar_face(edges, normal, [0.0; 3]);
        solid.faces[0].surface_type = surface_type;
        solid
    }

    fn mean_normal_z(mesh: &TriangleMesh) -> f64 {
        mesh.normals.iter().map(|n| n[2]).sum::<f64>() / mesh.normals.len() as f64
    }

    #[test]
    fn test_planar_face_follows_the_declared_normal_not_the_loop_direction() {
        let counter_clockwise = [0, 1, 2, 3];
        let clockwise = [0, 3, 2, 1];
        for (order, normal, expected) in [
            (counter_clockwise, [0.0, 0.0, 1.0], 1.0),
            (counter_clockwise, [0.0, 0.0, -1.0], -1.0),
            (clockwise, [0.0, 0.0, 1.0], 1.0),
            (clockwise, [0.0, 0.0, -1.0], -1.0),
        ] {
            let solid = square_face(order, normal, SurfaceType::Plane);
            let mesh = tessellate(&solid, 1).unwrap();
            assert_eq!(mesh.triangle_count(), 2);
            assert_close(mean_normal_z(&mesh), expected, 1e-12);
        }
    }

    #[test]
    fn test_planar_face_without_a_usable_normal_follows_the_loop() {
        for (order, expected) in [([0, 1, 2, 3], 1.0), ([0, 3, 2, 1], -1.0)] {
            let solid = square_face(order, [0.0; 3], SurfaceType::Plane);
            let mesh = tessellate(&solid, 1).unwrap();
            assert_close(mean_normal_z(&mesh), expected, 1e-12);
        }
    }

    #[test]
    fn test_curved_surface_types_follow_the_loop_winding() {
        for surface in [
            SurfaceType::Cylinder,
            SurfaceType::Sphere,
            SurfaceType::BSpline,
        ] {
            // Their `normal` field is an axis, not a facing direction.
            let up = tessellate(&square_face([0, 1, 2, 3], [0.0, 0.0, -1.0], surface), 1).unwrap();
            assert_close(mean_normal_z(&up), 1.0, 1e-12);
            let down = tessellate(&square_face([0, 3, 2, 1], [0.0, 0.0, 1.0], surface), 1).unwrap();
            assert_close(mean_normal_z(&down), -1.0, 1e-12);
        }
    }

    #[test]
    fn test_declared_normal_perpendicular_to_the_loop_is_an_error() {
        let solid = square_face([0, 1, 2, 3], [1.0, 0.0, 0.0], SurfaceType::Plane);
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(message.contains("perpendicular"), "{message}");
    }

    #[test]
    fn test_degenerate_outer_loop_is_an_error() {
        let mut solid = BrepSolid::new("flat");
        let a = solid.add_vertex([0.0, 0.0, 0.0]);
        let b = solid.add_vertex([1.0, 0.0, 0.0]);
        let c = solid.add_vertex([2.0, 0.0, 0.0]);
        let e0 = solid.add_line_edge(a, b);
        let e1 = solid.add_line_edge(b, c);
        let e2 = solid.add_line_edge(c, a);
        solid.add_planar_face(vec![e0, e1, e2], [0.0, 0.0, 1.0], [0.0; 3]);
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(message.contains("encloses no area"), "{message}");
    }

    /// A 4 x 4 plate with a 2 x 2 hole, placed in the plane through
    /// `(5, -3, 2)` with normal `(1, 1, 1) / sqrt(3)`.
    fn tilted_plate_with_hole() -> BrepSolid {
        let normal = [1.0, 1.0, 1.0].map(|c: f64| c / 3.0_f64.sqrt());
        let u = [1.0, -1.0, 0.0].map(|c: f64| c / 2.0_f64.sqrt());
        let v = [1.0, 1.0, -2.0].map(|c: f64| c / 6.0_f64.sqrt());
        let origin = [5.0, -3.0, 2.0];
        let place = |a: f64, b: f64| {
            [
                origin[0] + a * u[0] + b * v[0],
                origin[1] + a * u[1] + b * v[1],
                origin[2] + a * u[2] + b * v[2],
            ]
        };
        let mut solid = BrepSolid::new("plate");
        let mut loops = Vec::new();
        for (low, high) in [(0.0, 4.0), (1.0, 3.0)] {
            let ids = [
                solid.add_vertex(place(low, low)),
                solid.add_vertex(place(high, low)),
                solid.add_vertex(place(high, high)),
                solid.add_vertex(place(low, high)),
            ];
            let edges: Vec<usize> = (0..4)
                .map(|i| solid.add_line_edge(ids[i], ids[(i + 1) % 4]))
                .collect();
            loops.push(edges);
        }
        let hole = loops.pop().unwrap_or_default();
        solid.add_planar_face(loops.pop().unwrap_or_default(), normal, origin);
        solid.faces[0].edge_loops.push(hole);
        solid
    }

    #[test]
    fn test_face_with_a_hole_in_a_tilted_plane() {
        let solid = tilted_plate_with_hole();
        let mesh = tessellate(&solid, 1).unwrap();
        assert_close(mesh.surface_area(), 12.0, 1e-10);
        assert_eq!(mesh.triangle_count(), 8);
        let normal = solid.faces[0].normal;
        for n in &mesh.normals {
            assert!(dot3(*n, normal) > 0.999_999, "{n:?}");
        }
        // The hole is empty: no triangle contains the centre of the plate.
        let centre = {
            let a = solid.vertices[0].position;
            let c = solid.vertices[2].position;
            [
                (a[0] + c[0]) / 2.0,
                (a[1] + c[1]) / 2.0,
                (a[2] + c[2]) / 2.0,
            ]
        };
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]],
            );
            let n = cross3(sub3(b, a), sub3(c, a));
            let inside = [(a, b), (b, c), (c, a)]
                .iter()
                .all(|&(p, q)| dot3(cross3(sub3(q, p), sub3(centre, p)), n) > 1e-9);
            assert!(!inside, "triangle {t:?} covers the hole");
        }
    }

    #[test]
    fn test_plate_with_a_circular_hole_made_of_four_quarter_arcs() {
        let mut solid = BrepSolid::new("plate");
        let corners: Vec<usize> = [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
        ]
        .iter()
        .map(|&p| solid.add_vertex(p))
        .collect();
        let outer: Vec<usize> = (0..4)
            .map(|i| solid.add_line_edge(corners[i], corners[(i + 1) % 4]))
            .collect();
        let centre = [2.0, 2.0, 0.0];
        let rim: Vec<usize> = [
            [3.0, 2.0, 0.0],
            [2.0, 3.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
        .iter()
        .map(|&p| solid.add_vertex(p))
        .collect();
        let hole: Vec<usize> = (0..4)
            .map(|i| {
                let id = solid.edges.len();
                solid
                    .edges
                    .push(BrepEdge::arc(id, rim[i], rim[(i + 1) % 4], centre));
                id
            })
            .collect();
        solid.add_planar_face(outer, [0.0, 0.0, 1.0], [0.0; 3]);
        solid.faces[0].edge_loops.push(hole);

        let segments = 16;
        let mesh = tessellate(&solid, segments).unwrap();
        let hole_polygon = 4 * segments;
        assert_eq!(mesh.vertex_count(), 8 + 4 * (segments - 1));
        // Outer 4 + hole 64 vertices, the bridge adds two more, minus two.
        assert_eq!(mesh.triangle_count(), 4 + hole_polygon);
        let disc = 0.5 * hole_polygon as f64 * (std::f64::consts::TAU / hole_polygon as f64).sin();
        assert_close(mesh.surface_area(), 16.0 - disc, 1e-9);
        for normal in &mesh.normals {
            assert!(normal[2] > 0.999_999, "{normal:?}");
        }
        for p in &mesh.vertices[4..] {
            assert_close(len3(sub3(*p, centre)), 1.0, 1e-12);
        }
    }

    #[test]
    fn test_hole_outside_the_outer_loop_is_an_error() {
        let mut solid = tilted_plate_with_hole();
        for vertex in &mut solid.vertices[4..] {
            let p = vertex.position;
            vertex.position = [p[0] + 40.0, p[1] - 40.0, p[2]];
        }
        let message = tessellate(&solid, 1).unwrap_err();
        assert!(message.contains("face 0"), "{message}");
    }

    /// Two faces back to back over the same boundary: the quarter disc with
    /// centre `v0`, one straight edge along x, a quarter-circle arc and one
    /// straight edge along y.  The pair is a closed (flat) surface.
    fn pillow() -> BrepSolid {
        let mut solid = BrepSolid::new("pillow");
        let centre = solid.add_vertex([0.0, 0.0, 0.0]);
        let start = solid.add_vertex([1.0, 0.0, 0.0]);
        let end = solid.add_vertex([0.0, 1.0, 0.0]);
        let along_x = solid.add_line_edge(centre, start);
        let arc_id = solid.edges.len();
        solid
            .edges
            .push(BrepEdge::arc(arc_id, start, end, [0.0, 0.0, 0.0]));
        let along_y = solid.add_line_edge(end, centre);
        solid.add_planar_face(vec![along_x, arc_id, along_y], [0.0, 0.0, 1.0], [0.0; 3]);
        solid.add_planar_face(vec![along_y, arc_id, along_x], [0.0, 0.0, -1.0], [0.0; 3]);
        solid
    }

    #[test]
    fn test_arc_edges_are_sampled_and_the_samples_are_shared() {
        let segments = 64;
        let mesh = tessellate(&pillow(), segments).unwrap();
        // 3 corners plus the interior samples of the arc, created once.
        assert_eq!(mesh.vertex_count(), 3 + (segments - 1));
        assert_eq!(mesh.triangle_count(), 2 * (3 + (segments - 1) - 2));
        assert_balanced(&mesh);
        assert_close(volume(&mesh), 0.0, 1e-12);
        // Every sample lies on the unit circle.
        for p in &mesh.vertices[3..] {
            assert_close(len3(*p), 1.0, 1e-12);
        }
        // Area of the fan polygon, twice (two faces).
        let polygon_area =
            0.5 * segments as f64 * (std::f64::consts::FRAC_PI_2 / segments as f64).sin();
        assert_close(mesh.surface_area(), 2.0 * polygon_area, 1e-12);
        assert_close(mesh.surface_area(), std::f64::consts::FRAC_PI_2, 1e-3);
    }

    #[test]
    fn test_faces_sharing_an_arc_agree_for_every_segment_count() {
        for segments in [1, 2, 3, 7, 16] {
            let mesh = tessellate(&pillow(), segments).unwrap();
            assert_eq!(
                mesh.vertex_count(),
                3 + (segments - 1),
                "segments {segments}"
            );
            assert_balanced(&mesh);
        }
    }

    #[test]
    fn test_assembly_and_exporter_report_tessellation_errors() {
        let mut broken = BrepSolid::create_box("broken", 1.0, 1.0, 1.0);
        broken.faces[0].edge_loops[0] = vec![0, 2, 4, 6];
        let message = StlExporter::new("x")
            .export_brep_ascii(&broken)
            .unwrap_err();
        assert!(message.contains("not connected"), "{message}");

        let good = BrepSolid::create_box("good", 1.0, 1.0, 1.0);
        let stl = StlExporter::new("x").export_brep_ascii(&good).unwrap();
        assert_eq!(stl.matches("facet normal").count(), 12);

        let mut root = AssemblyComponent::new("root", good, AssemblyTransform::identity());
        root.add_child(AssemblyComponent::new(
            "bad_child",
            broken,
            AssemblyTransform::translation([5.0, 0.0, 0.0]),
        ));
        let assembly = Assembly::new("asm", root, LengthUnit::Millimeter);
        let message = assembly.to_stl().unwrap_err();
        assert!(message.contains("bad_child"), "{message}");
    }

    #[test]
    fn test_assembly_stl_applies_the_component_transform() {
        let root = AssemblyComponent::new(
            "root",
            BrepSolid::create_box("box", 1.0, 1.0, 1.0),
            AssemblyTransform::translation([10.0, 0.0, 0.0]),
        );
        let assembly = Assembly::new("moved", root, LengthUnit::Millimeter);
        let stl = assembly.to_stl().unwrap();
        let xs: Vec<f64> = stl
            .lines()
            .filter_map(|l| l.trim().strip_prefix("vertex "))
            .filter_map(|rest| rest.split_whitespace().next()?.parse().ok())
            .collect();
        assert_eq!(xs.len(), 36);
        assert!(xs.iter().all(|x| (10.0..=11.0).contains(x)), "{xs:?}");
    }
}
