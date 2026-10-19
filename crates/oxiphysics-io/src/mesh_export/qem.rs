// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quadric-error-metric edge-collapse simplification (Garland and Heckbert, 1997).
//!
//! Every vertex carries a *quadric*: a symmetric 4x4 matrix `Q` for which
//! `v^T Q v`, with `v = (x, y, z, 1)`, is the weighted sum of the squared distances
//! from the point `(x, y, z)` to the planes of the triangles that meet at the
//! vertex.  Collapsing the edge `(a, b)` replaces both vertices by one vertex at the
//! position `v` and gives it the quadric `Q_a + Q_b`, so the error of a vertex
//! accumulates over the whole sequence of collapses instead of being measured against
//! the original surface of one vertex only.
//!
//! The simplifier
//!
//! 1. builds the quadrics from the area-weighted face planes and adds one heavily
//!    weighted plane per boundary edge (perpendicular to the face, through the edge),
//!    so open boundaries keep their outline;
//! 2. costs every edge at its optimal placement: the solution of the 3x3 system
//!    `A v = -b`; when that system is (nearly) singular, or its solution lies more
//!    than one edge length away from the edge midpoint, the best point on the edge
//!    itself (the minimum of the quadric along the segment, or else the cheapest of
//!    the two endpoints and the midpoint);
//! 3. repeatedly collapses the cheapest edge taken from a min-heap.  Heap entries carry
//!    the version stamps of both endpoints; a collapse bumps the survivor's version, so
//!    every entry that mentions it becomes stale and is skipped, and only the edges of
//!    the survivor are costed again;
//! 4. refuses collapses that would change the topology (the link condition) or fold a
//!    surrounding triangle over (its normal would turn by 90 degrees or more).
//!
//! Vertex-to-face adjacency is kept up to date incrementally, so a collapse costs time
//! proportional to the valence of the two vertices rather than to the size of the mesh.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::{MeshData, compact_mesh};

type Vec3 = [f64; 3];

/// Weight of the plane that pins a boundary edge, relative to the squared edge length.
///
/// The area-weighted face quadrics are of the order of the squared edge length, so a
/// weight of 1000 makes leaving the boundary about three orders of magnitude more
/// expensive than leaving the surface.
const BOUNDARY_WEIGHT: f64 = 1000.0;

/// A 3x3 system is singular when its determinant is at most this fraction of the cube
/// of its largest entry.
const SINGULARITY_TOLERANCE: f64 = 1e-10;

/// The quadric is flat along an edge when the curvature along it is at most this
/// fraction of `max|A_ij| * |edge|^2`.
const FLAT_CURVATURE_TOLERANCE: f64 = 1e-12;

/// The optimal position may be at most this many edge lengths away from the midpoint
/// of its edge; further away it is the symptom of an ill-conditioned system.
const MAX_PLACEMENT_SHIFT: f64 = 1.0;

/// A triangle whose smallest angle has a sine below this is degenerate.
const DEGENERATE_SINE: f64 = 1e-9;

/// The cosine of the angle between the old and the new normal of a surrounding triangle
/// must exceed this value; zero refuses a turn of 90 degrees or more.
const MIN_NORMAL_COSINE: f64 = 0.0;

/// A position this close to an endpoint of the edge (as a fraction of the edge length) is
/// rounding noise of a minimum that lies exactly at that endpoint.
const ENDPOINT_SNAP: f64 = 1e-12;

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

fn lerp(a: Vec3, b: Vec3, t: f64) -> Vec3 {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

/// The point `from + t (to - from)` for `t` in `[0, 1]`; the endpoints themselves (bit
/// for bit) when `t` is within `ENDPOINT_SNAP` of 0 or 1.
fn point_on_segment(from: Vec3, to: Vec3, t: f64) -> Vec3 {
    if t <= ENDPOINT_SNAP {
        from
    } else if t >= 1.0 - ENDPOINT_SNAP {
        to
    } else {
        lerp(from, to, t)
    }
}

/// `point`, or the endpoint of the edge `from`-`to` (of the given `length`) that it
/// coincides with up to rounding noise.
///
/// The solution of a 3x3 system that is meant to be an existing vertex, such as the
/// corner where three planes of a box meet, is off by a few units in the last place;
/// snapping keeps the coordinates of such vertices exact.
fn snap_to_endpoint(from: Vec3, to: Vec3, point: Vec3, length: f64) -> Vec3 {
    let tolerance = ENDPOINT_SNAP * length;
    if norm(sub(point, from)) <= tolerance {
        from
    } else if norm(sub(point, to)) <= tolerance {
        to
    } else {
        point
    }
}

/// Symmetric 4x4 quadric `[[A, b], [b^T, c]]` stored as its 10 upper-triangular
/// entries `[A00, A01, A02, b0, A11, A12, b1, A22, b2, c]`.
///
/// The error of the point `v` is `v^T A v + 2 b^T v + c`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct Quadric {
    q: [f64; 10],
}

impl Quadric {
    /// The quadric of the plane `normal . x + offset = 0` (with a unit `normal`),
    /// scaled by `weight`.
    fn from_plane(weight: f64, normal: Vec3, offset: f64) -> Self {
        let [a, b, c] = normal;
        Self {
            q: [
                weight * a * a,
                weight * a * b,
                weight * a * c,
                weight * a * offset,
                weight * b * b,
                weight * b * c,
                weight * b * offset,
                weight * c * c,
                weight * c * offset,
                weight * offset * offset,
            ],
        }
    }

    fn add(&mut self, other: &Quadric) {
        for (sum, term) in self.q.iter_mut().zip(other.q.iter()) {
            *sum += *term;
        }
    }

    fn sum(&self, other: &Quadric) -> Quadric {
        let mut total = *self;
        total.add(other);
        total
    }

    /// `A v` for the upper-left 3x3 block `A`.
    fn apply_matrix(&self, v: Vec3) -> Vec3 {
        let q = &self.q;
        [
            q[0] * v[0] + q[1] * v[1] + q[2] * v[2],
            q[1] * v[0] + q[4] * v[1] + q[5] * v[2],
            q[2] * v[0] + q[5] * v[1] + q[7] * v[2],
        ]
    }

    /// The vector `b` of the linear term.
    fn linear_term(&self) -> Vec3 {
        [self.q[3], self.q[6], self.q[8]]
    }

    /// Largest absolute entry of `A`.
    fn matrix_scale(&self) -> f64 {
        let q = &self.q;
        [q[0], q[1], q[2], q[4], q[5], q[7]]
            .iter()
            .fold(0.0_f64, |scale, entry| scale.max(entry.abs()))
    }

    /// Sum of the weighted squared plane distances at `v`.
    fn error_at(&self, v: Vec3) -> f64 {
        dot(v, self.apply_matrix(v)) + 2.0 * dot(self.linear_term(), v) + self.q[9]
    }

    /// The point that minimizes the error: the solution of `A v = -b`.
    ///
    /// `None` when `A` is singular (the planes do not meet in a single point), judged
    /// relative to the size of its entries, or when the solution is not finite.
    fn minimizer(&self) -> Option<Vec3> {
        let q = &self.q;
        let (a00, a01, a02, a11, a12, a22) = (q[0], q[1], q[2], q[4], q[5], q[7]);
        let rhs = [-q[3], -q[6], -q[8]];
        // Cofactors of the symmetric matrix `A`; `adj(A) = C`.
        let c00 = a11 * a22 - a12 * a12;
        let c01 = a02 * a12 - a01 * a22;
        let c02 = a01 * a12 - a02 * a11;
        let c11 = a00 * a22 - a02 * a02;
        let c12 = a01 * a02 - a00 * a12;
        let c22 = a00 * a11 - a01 * a01;
        let det = a00 * c00 + a01 * c01 + a02 * c02;
        let scale = self.matrix_scale();
        // Written so that a NaN determinant is also rejected.
        if det.abs() <= SINGULARITY_TOLERANCE * scale * scale * scale || det.is_nan() {
            return None;
        }
        let point = [
            (c00 * rhs[0] + c01 * rhs[1] + c02 * rhs[2]) / det,
            (c01 * rhs[0] + c11 * rhs[1] + c12 * rhs[2]) / det,
            (c02 * rhs[0] + c12 * rhs[1] + c22 * rhs[2]) / det,
        ];
        point.iter().all(|c| c.is_finite()).then_some(point)
    }

    /// The parameter `t` in `[0, 1]` that minimizes the error along the segment
    /// `from + t (to - from)`.
    ///
    /// `None` when the error is constant (or numerically so) along the segment, in
    /// which case every point of it is as good as the others.
    fn minimum_on_segment(&self, from: Vec3, to: Vec3) -> Option<f64> {
        let direction = sub(to, from);
        let curvature = dot(direction, self.apply_matrix(direction));
        let flat = FLAT_CURVATURE_TOLERANCE * self.matrix_scale() * dot(direction, direction);
        if curvature <= flat || curvature.is_nan() {
            return None;
        }
        let slope = dot(add(self.apply_matrix(from), self.linear_term()), direction);
        let t = -slope / curvature;
        t.is_finite().then(|| t.clamp(0.0, 1.0))
    }
}

/// Unit normal and twice the area of the triangle `p0 p1 p2`; `None` for a triangle
/// whose three corners are collinear.
fn face_plane(p0: Vec3, p1: Vec3, p2: Vec3) -> Option<(Vec3, f64)> {
    let n = cross(sub(p1, p0), sub(p2, p0));
    let double_area = norm(n);
    (double_area > 0.0 && double_area.is_finite()).then(|| {
        (
            [n[0] / double_area, n[1] / double_area, n[2] / double_area],
            double_area,
        )
    })
}

/// An edge collapse waiting in the priority queue.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    /// Quadric error at `target`.
    cost: f64,
    /// Squared edge length; shorter edges win ties of equal cost.
    length_sq: f64,
    /// Lower vertex index of the edge.
    a: usize,
    /// Higher vertex index of the edge.
    b: usize,
    /// Version of vertex `a` when the candidate was costed.
    version_a: u32,
    /// Version of vertex `b` when the candidate was costed.
    version_b: u32,
    /// Position of the merged vertex.
    target: Vec3,
}

impl Ord for Candidate {
    /// A `BinaryHeap` pops its greatest element, so the order is reversed: the
    /// cheapest edge is the greatest.  `total_cmp` and the vertex pair make the order
    /// total, so the result never depends on hashing or insertion order.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.length_sq.total_cmp(&self.length_sq))
            .then_with(|| other.a.cmp(&self.a))
            .then_with(|| other.b.cmp(&self.b))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Candidate {}

/// Costs are never NaN or negative: NaN sorts last as `+inf`, rounding noise below zero
/// (and `-0.0`) is `0.0`.
fn finite_cost(cost: f64) -> f64 {
    if cost.is_nan() {
        f64::INFINITY
    } else if cost > 0.0 {
        cost
    } else {
        0.0
    }
}

/// Whether a triangle with the edge vectors `e1` and `e2` from one corner, whose cross
/// product has the length `cross_length`, has collapsed to a sliver or a line.
fn is_degenerate(cross_length: f64, e1: Vec3, e2: Vec3) -> bool {
    cross_length <= DEGENERATE_SINE * norm(e1) * norm(e2)
}

/// Whether moving a triangle from `before` to `after` keeps its normal within the
/// allowed turn and keeps it a proper triangle.
///
/// A triangle that already is degenerate has no normal to preserve.
fn keeps_orientation(before: [Vec3; 3], after: [Vec3; 3]) -> bool {
    let old_edges = [sub(before[1], before[0]), sub(before[2], before[0])];
    let new_edges = [sub(after[1], after[0]), sub(after[2], after[0])];
    let old_normal = cross(old_edges[0], old_edges[1]);
    let new_normal = cross(new_edges[0], new_edges[1]);
    let old_length = norm(old_normal);
    let new_length = norm(new_normal);
    if is_degenerate(old_length, old_edges[0], old_edges[1]) {
        return true;
    }
    if is_degenerate(new_length, new_edges[0], new_edges[1]) {
        return false;
    }
    dot(old_normal, new_normal) > MIN_NORMAL_COSINE * old_length * new_length
}

/// The state of one simplification run.
struct Simplifier {
    /// Vertex positions, relative to `centre`.
    positions: Vec<Vec3>,
    /// The coordinates of every vertex as the caller sees them: the input coordinates
    /// until a collapse moves the vertex, and bit for bit the coordinates of an endpoint
    /// when the collapse puts the merged vertex on that endpoint.
    absolute: Vec<Vec3>,
    /// The origin that `positions` are relative to.
    centre: Vec3,
    /// Accumulated quadric of every vertex.
    quadrics: Vec<Quadric>,
    /// Bumped whenever the position or quadric of a vertex changes.
    versions: Vec<u32>,
    /// Vertices that have not been merged into another one.
    vertex_alive: Vec<bool>,
    /// Triangles as vertex indices; the corners of live triangles are live vertices.
    faces: Vec<[usize; 3]>,
    /// Triangles that have not vanished in a collapse.
    face_alive: Vec<bool>,
    /// For every vertex the live triangles that use it.
    vertex_faces: Vec<Vec<usize>>,
    /// Number of live triangles.
    live_faces: usize,
    /// Edges waiting to be collapsed, cheapest first.
    heap: BinaryHeap<Candidate>,
    /// Edges whose collapse was refused; they are tried again once the neighbourhood
    /// has changed.
    deferred: Vec<(usize, usize)>,
}

impl Simplifier {
    /// Set up quadrics and adjacency for `faces`, whose corners must all be valid,
    /// pairwise distinct indices into `absolute`.
    ///
    /// The optimisation works on the coordinates relative to `centre`: quadrics hold
    /// squared plane offsets, which would lose their precision for a mesh far away from
    /// the origin.
    fn new(absolute: Vec<Vec3>, centre: Vec3, faces: Vec<[usize; 3]>) -> Self {
        let positions: Vec<Vec3> = absolute.iter().map(|&p| sub(p, centre)).collect();
        let vertex_count = positions.len();
        let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
        let mut quadrics = vec![Quadric::default(); vertex_count];
        // One record per triangle edge: (low vertex, high vertex, triangle).
        let mut edge_records: Vec<(usize, usize, usize)> = Vec::with_capacity(3 * faces.len());
        for (index, face) in faces.iter().enumerate() {
            let [p0, p1, p2] = face.map(|v| positions[v]);
            if let Some((normal, double_area)) = face_plane(p0, p1, p2) {
                let plane = Quadric::from_plane(0.5 * double_area, normal, -dot(normal, p0));
                for &v in face {
                    quadrics[v].add(&plane);
                }
            }
            for (slot, &v) in face.iter().enumerate() {
                vertex_faces[v].push(index);
                let w = face[(slot + 1) % 3];
                edge_records.push((v.min(w), v.max(w), index));
            }
        }
        edge_records.sort_unstable();
        let mut start = 0;
        while start < edge_records.len() {
            let (low, high, face) = edge_records[start];
            let mut end = start + 1;
            while end < edge_records.len()
                && (edge_records[end].0, edge_records[end].1) == (low, high)
            {
                end += 1;
            }
            if end - start == 1
                && let Some(constraint) = boundary_constraint(
                    positions[low],
                    positions[high],
                    faces[face].map(|v| positions[v]),
                )
            {
                quadrics[low].add(&constraint);
                quadrics[high].add(&constraint);
            }
            start = end;
        }
        let live_faces = faces.len();
        Self {
            positions,
            absolute,
            centre,
            quadrics,
            versions: vec![0; vertex_count],
            vertex_alive: vec![true; vertex_count],
            face_alive: vec![true; live_faces],
            faces,
            vertex_faces,
            live_faces,
            heap: BinaryHeap::new(),
            deferred: Vec::new(),
        }
    }

    /// Collapse edges, cheapest first, until at most `target_faces` triangles are left
    /// or no collapse is allowed any more.
    fn run(&mut self, target_faces: usize) {
        if self.live_faces <= target_faces {
            return;
        }
        self.seed_all_edges();
        loop {
            let mut progressed = false;
            while self.live_faces > target_faces {
                let Some(candidate) = self.heap.pop() else {
                    break;
                };
                if !self.is_current(&candidate) {
                    continue;
                }
                if !self.allows_collapse(candidate.a, candidate.b, candidate.target) {
                    self.deferred.push((candidate.a, candidate.b));
                    continue;
                }
                self.collapse(candidate.a, candidate.b, candidate.target);
                self.refresh_edges(candidate.a);
                progressed = true;
            }
            if self.live_faces <= target_faces || !progressed || self.deferred.is_empty() {
                break;
            }
            self.reseed_deferred();
        }
    }

    /// Push every edge of the live triangles.
    fn seed_all_edges(&mut self) {
        let mut edges: Vec<(usize, usize)> = Vec::with_capacity(3 * self.live_faces);
        for (face, &alive) in self.faces.iter().zip(&self.face_alive) {
            if !alive {
                continue;
            }
            for (slot, &v) in face.iter().enumerate() {
                let w = face[(slot + 1) % 3];
                edges.push((v.min(w), v.max(w)));
            }
        }
        edges.sort_unstable();
        edges.dedup();
        for (a, b) in edges {
            self.push_edge(a, b);
        }
    }

    /// Cost the edges whose collapse was refused again, now that the mesh around them
    /// may have changed.
    fn reseed_deferred(&mut self) {
        let mut edges = std::mem::take(&mut self.deferred);
        edges.sort_unstable();
        edges.dedup();
        for (a, b) in edges {
            if self.vertex_alive[a] && self.vertex_alive[b] && self.share_face(a, b) {
                self.push_edge(a, b);
            }
        }
    }

    /// Cost the edges of the survivor of a collapse again.
    fn refresh_edges(&mut self, survivor: usize) {
        for (neighbour, _) in self.neighbour_counts(survivor) {
            self.push_edge(survivor, neighbour);
        }
    }

    /// Whether `candidate` was costed against the current state of its endpoints.
    fn is_current(&self, candidate: &Candidate) -> bool {
        self.vertex_alive[candidate.a]
            && self.vertex_alive[candidate.b]
            && self.versions[candidate.a] == candidate.version_a
            && self.versions[candidate.b] == candidate.version_b
    }

    fn share_face(&self, a: usize, b: usize) -> bool {
        self.vertex_faces[a]
            .iter()
            .any(|&f| self.faces[f].contains(&b))
    }

    /// Whether some live triangle has the corners `a`, `p` and `q`.
    fn has_face(&self, a: usize, p: usize, q: usize) -> bool {
        self.vertex_faces[a].iter().any(|&f| {
            let face = &self.faces[f];
            face.contains(&p) && face.contains(&q)
        })
    }

    /// Cost the edge `(a, b)` and put it in the queue.
    fn push_edge(&mut self, a: usize, b: usize) {
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        let (target, cost) = self.best_target(a, b);
        let span = sub(self.positions[b], self.positions[a]);
        self.heap.push(Candidate {
            cost: finite_cost(cost),
            length_sq: dot(span, span),
            a,
            b,
            version_a: self.versions[a],
            version_b: self.versions[b],
            target,
        });
    }

    /// Where the merged vertex of the edge `(a, b)` goes, and what that costs.
    fn best_target(&self, a: usize, b: usize) -> (Vec3, f64) {
        let quadric = self.quadrics[a].sum(&self.quadrics[b]);
        let (from, to) = (self.positions[a], self.positions[b]);
        let midpoint = lerp(from, to, 0.5);
        let length = norm(sub(to, from));
        if let Some(optimum) = quadric.minimizer()
            && norm(sub(optimum, midpoint)) <= MAX_PLACEMENT_SHIFT * length
        {
            let optimum = snap_to_endpoint(from, to, optimum, length);
            return (optimum, quadric.error_at(optimum));
        }
        if let Some(t) = quadric.minimum_on_segment(from, to) {
            let point = point_on_segment(from, to, t);
            return (point, quadric.error_at(point));
        }
        // Flat along the edge: the cheapest of the endpoints and the midpoint, the
        // earlier one of them on ties.
        let mut best = (from, quadric.error_at(from));
        for point in [to, midpoint] {
            let error = quadric.error_at(point);
            if error < best.1 {
                best = (point, error);
            }
        }
        best
    }

    /// Adjacent vertices of `vertex` with the number of live triangles that contain the
    /// edge to them, sorted by vertex index.
    fn neighbour_counts(&self, vertex: usize) -> Vec<(usize, usize)> {
        let mut others: Vec<usize> = Vec::with_capacity(2 * self.vertex_faces[vertex].len());
        for &f in &self.vertex_faces[vertex] {
            others.extend(self.faces[f].iter().copied().filter(|&u| u != vertex));
        }
        others.sort_unstable();
        let mut counts: Vec<(usize, usize)> = Vec::new();
        for u in others {
            match counts.last_mut() {
                Some((last, count)) if *last == u => *count += 1,
                _ => counts.push((u, 1)),
            }
        }
        counts
    }

    /// Whether the collapse of `(a, b)` to `target` keeps the surface a manifold and
    /// does not fold over a surrounding triangle.
    fn allows_collapse(&self, a: usize, b: usize, target: Vec3) -> bool {
        self.preserves_topology(a, b) && self.preserves_orientation(a, b, target)
    }

    /// The link condition of Dey, Edelsbrunner, Guha and Nekhayev for the surface closed
    /// up by a cone over its boundary: collapsing the edge `(a, b)` keeps the topology of
    /// the surface exactly when
    ///
    /// * the vertices adjacent to both ends are exactly the apexes of the triangles on the
    ///   edge (otherwise the collapse would pinch the surface);
    /// * an interior edge does not join two boundary vertices (that would glue two
    ///   boundary pieces together);
    /// * a boundary edge is not an edge of a triangle whose other two edges are boundary
    ///   edges as well, an ear that touches the rest of the surface at vertices only or is a
    ///   whole component, which the collapse would cut off or delete;
    /// * the two triangles at an interior edge are not part of a tetrahedron, whose
    ///   remaining two triangles would coincide.
    ///
    /// The number of connected pieces and the Euler characteristic therefore never change.
    /// Edges next to non-manifold edges are left alone.
    fn preserves_topology(&self, a: usize, b: usize) -> bool {
        let around_a = self.neighbour_counts(a);
        let around_b = self.neighbour_counts(b);
        if around_a
            .iter()
            .chain(&around_b)
            .any(|&(_, count)| count > 2)
        {
            return false;
        }
        let shared = around_a
            .iter()
            .find(|&&(u, _)| u == b)
            .map_or(0, |&(_, count)| count);
        if shared == 0 {
            return false;
        }
        let mut apexes: Vec<usize> = Vec::with_capacity(shared);
        for &f in &self.vertex_faces[a] {
            let face = &self.faces[f];
            if face.contains(&b)
                && let Some(&apex) = face.iter().find(|&&v| v != a && v != b)
            {
                apexes.push(apex);
            }
        }
        apexes.sort_unstable();
        let common: Vec<usize> = around_a
            .iter()
            .map(|&(u, _)| u)
            .filter(|&u| u != b && around_b.binary_search_by_key(&u, |&(w, _)| w).is_ok())
            .collect();
        if apexes != common {
            return false;
        }
        if let [apex] = apexes[..] {
            // A boundary edge whose triangle has no other neighbour across an edge is
            // an ear that touches the rest of the surface at vertices only (or a whole
            // component): collapsing it would cut or delete a piece of the surface.
            let on_boundary_edge = |counts: &[(usize, usize)]| {
                counts.iter().any(|&(u, count)| u == apex && count == 1)
            };
            if on_boundary_edge(&around_a) && on_boundary_edge(&around_b) {
                return false;
            }
        }
        if shared == 2 {
            let on_boundary = |counts: &[(usize, usize)]| counts.iter().any(|&(_, c)| c == 1);
            if on_boundary(&around_a) && on_boundary(&around_b) {
                return false;
            }
            if let [p, q] = apexes[..]
                && self.has_face(a, p, q)
                && self.has_face(b, p, q)
            {
                return false;
            }
        }
        true
    }

    /// No triangle that survives the collapse may turn its normal by 90 degrees or
    /// more, or shrink to a sliver.
    fn preserves_orientation(&self, a: usize, b: usize, target: Vec3) -> bool {
        for (moving, other) in [(a, b), (b, a)] {
            for &f in &self.vertex_faces[moving] {
                let face = self.faces[f];
                if face.contains(&other) {
                    continue;
                }
                let before = face.map(|v| self.positions[v]);
                let after = face.map(|v| {
                    if v == moving {
                        target
                    } else {
                        self.positions[v]
                    }
                });
                if !keeps_orientation(before, after) {
                    return false;
                }
            }
        }
        true
    }

    /// Remove a triangle from the mesh.
    fn kill_face(&mut self, face_index: usize) {
        self.face_alive[face_index] = false;
        self.live_faces -= 1;
        for v in self.faces[face_index] {
            self.vertex_faces[v].retain(|&g| g != face_index);
        }
    }

    /// Merge `b` into `a`, which moves to `target` and takes over the quadric of `b`.
    ///
    /// The triangles on the edge vanish; the other triangles of `b` become triangles of
    /// `a`.
    fn collapse(&mut self, a: usize, b: usize, target: Vec3) {
        let vanishing: Vec<usize> = self.vertex_faces[a]
            .iter()
            .copied()
            .filter(|&f| self.faces[f].contains(&b))
            .collect();
        for f in vanishing {
            self.kill_face(f);
        }
        for f in std::mem::take(&mut self.vertex_faces[b]) {
            for corner in self.faces[f].iter_mut() {
                if *corner == b {
                    *corner = a;
                }
            }
            self.vertex_faces[a].push(f);
        }
        self.vertex_alive[b] = false;
        let absorbed = self.quadrics[b];
        self.quadrics[a].add(&absorbed);
        self.absolute[a] = if target == self.positions[a] {
            self.absolute[a]
        } else if target == self.positions[b] {
            self.absolute[b]
        } else {
            add(target, self.centre)
        };
        self.positions[a] = target;
        self.versions[a] = self.versions[a].wrapping_add(1);
    }

    /// The live triangles in their original order, with the index each had in the input.
    fn surviving_faces(&self, origin: &[usize]) -> (Vec<[usize; 3]>, Vec<usize>) {
        let mut faces = Vec::with_capacity(self.live_faces);
        let mut indices = Vec::with_capacity(self.live_faces);
        for ((face, &alive), &input_index) in self.faces.iter().zip(&self.face_alive).zip(origin) {
            if alive {
                faces.push(*face);
                indices.push(input_index);
            }
        }
        (faces, indices)
    }
}

/// The plane that pins the boundary edge `(from, to)` of `triangle`: it contains the
/// edge, is perpendicular to the triangle, and is weighted by `BOUNDARY_WEIGHT` times
/// the squared edge length.
fn boundary_constraint(from: Vec3, to: Vec3, triangle: [Vec3; 3]) -> Option<Quadric> {
    let edge = sub(to, from);
    let face_normal = cross(sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]));
    let across = cross(edge, face_normal);
    let length = norm(across);
    if !(length > 0.0 && length.is_finite()) {
        return None;
    }
    let unit = [across[0] / length, across[1] / length, across[2] / length];
    Some(Quadric::from_plane(
        BOUNDARY_WEIGHT * dot(edge, edge),
        unit,
        -dot(unit, from),
    ))
}

/// Centre of the bounding box of the corners of `faces`.
fn bounding_box_centre(vertices: &[Vec3], faces: &[[usize; 3]]) -> Vec3 {
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for &v in faces.iter().flatten() {
        let Some(point) = vertices.get(v) else {
            continue;
        };
        for ((lo, hi), &coordinate) in low.iter_mut().zip(high.iter_mut()).zip(point) {
            *lo = lo.min(coordinate);
            *hi = hi.max(coordinate);
        }
    }
    if low.iter().chain(&high).all(|c| c.is_finite()) {
        [
            0.5 * (low[0] + high[0]),
            0.5 * (low[1] + high[1]),
            0.5 * (low[2] + high[2]),
        ]
    } else {
        [0.0; 3]
    }
}

/// Simplify `mesh` to at most `target_triangles` triangles where the collapse rules
/// allow it; see [`MeshSimplification::simplify`](super::MeshSimplification::simplify).
pub(super) fn simplify(mesh: &MeshData, target_triangles: usize) -> MeshData {
    let usable = |face: &[usize; 3]| {
        face[0] != face[1]
            && face[1] != face[2]
            && face[0] != face[2]
            && face.iter().all(|&v| {
                mesh.vertices
                    .get(v)
                    .is_some_and(|p| p.iter().all(|c| c.is_finite()))
            })
    };
    let mut faces: Vec<[usize; 3]> = Vec::with_capacity(mesh.faces.len());
    let mut origin: Vec<usize> = Vec::with_capacity(mesh.faces.len());
    for (index, face) in mesh.faces.iter().enumerate() {
        if usable(face) {
            faces.push(*face);
            origin.push(index);
        }
    }
    let centre = bounding_box_centre(&mesh.vertices, &faces);
    let mut simplifier = Simplifier::new(mesh.vertices.clone(), centre, faces);
    simplifier.run(target_triangles);
    let (kept_faces, kept_origin) = simplifier.surviving_faces(&origin);
    let mut simplified = compact_mesh(&simplifier.absolute, &kept_faces);
    if mesh.material_ids.len() == mesh.faces.len() {
        simplified.material_ids = kept_origin
            .iter()
            .filter_map(|&f| mesh.material_ids.get(f).copied())
            .collect();
    }
    simplified
}

#[cfg(test)]
mod tests;
