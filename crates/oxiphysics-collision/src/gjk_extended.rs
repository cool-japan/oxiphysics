// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended GJK/EPA algorithms.
//!
//! Full GJK with simplex management (1D–4D), EPA (Expanding Polytope Algorithm)
//! for penetration depth/normal, warm-start cache, Minkowski difference support,
//! 2-D GJK, ray-cast GJK, time-of-impact via conservative advancement, and
//! feature pair generation.

use crate::gjk_epa::distance::{DistanceOutcome, DistanceSettings, SimplexVertex, solve_with};
use crate::gjk_epa::simplex_closest::{SimplexClosest, closest_point_to_origin};
use crate::gjk_epa::{SupportPoint, complete_tetrahedron, expand_polytope};

// ---------------------------------------------------------------------------
// Vector math helpers
// ---------------------------------------------------------------------------

/// Compute a − b.
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Compute a + b.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale vector by s.
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared length.
pub fn len_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

/// Length.
pub fn len3(a: [f64; 3]) -> f64 {
    len_sq3(a).sqrt()
}

/// Normalize (returns \[0,0,0\] for near-zero vectors).
pub fn norm3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-300 {
        return [0.0; 3];
    }
    scale3(a, 1.0 / l)
}

/// Negate vector.
pub fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ---------------------------------------------------------------------------
// Closest-point helpers (public)
// ---------------------------------------------------------------------------

/// Closest point on a line segment AB to the origin.
///
/// Returns the clamped parameter t ∈ \[0,1\] and the closest point.
pub fn line_closest_point_to_origin(a: [f64; 3], b: [f64; 3]) -> (f64, [f64; 3]) {
    let ab = sub3(b, a);
    let len_sq = len_sq3(ab);
    if len_sq < 1e-24 {
        return (0.0, a);
    }
    let t = (-dot3(a, ab) / len_sq).clamp(0.0, 1.0);
    let p = add3(a, scale3(ab, t));
    (t, p)
}

/// Closest point on a triangle (A,B,C) to the origin.
///
/// Returns the closest point.
pub fn triangle_closest_to_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    // Use GJK sub-simplex check: check all sub-features
    // Vertex regions
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ao = neg3(a);

    let d1 = dot3(ab, ao);
    let d2 = dot3(ac, ao);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }

    let bo = neg3(b);
    let d3 = dot3(ab, bo);
    let d4 = dot3(ac, bo);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }

    let co = neg3(c);
    let d5 = dot3(ab, co);
    let d6 = dot3(ac, co);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }

    // Edge AB
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add3(a, scale3(ab, v));
    }

    // Edge AC
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add3(a, scale3(ac, w));
    }

    // Edge BC
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add3(b, scale3(sub3(c, b), w));
    }

    // Inside triangle
    let denom = 1.0 / (va + vb + vc).max(1e-300);
    let v = vb * denom;
    let w = vc * denom;
    add3(a, add3(scale3(ab, v), scale3(ac, w)))
}

/// Check whether the tetrahedron (A,B,C,D) contains the origin.
pub fn tetrahedron_contains_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    // Check that origin is on the correct side of all four faces
    let faces = [(a, b, c, d), (a, c, d, b), (a, d, b, c), (b, d, c, a)];
    for (p0, p1, p2, p3) in &faces {
        let n = cross3(sub3(*p1, *p0), sub3(*p2, *p0));
        let inside = dot3(n, sub3(*p3, *p0)) > 0.0;
        let origin_side = dot3(n, neg3(*p0)) > 0.0;
        if inside != origin_side {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Simplex manager
// ---------------------------------------------------------------------------

/// GJK simplex with up to 4 vertices (3-D).
#[derive(Clone, Debug)]
pub struct SimplexManager {
    /// Vertices of the simplex (Minkowski difference points).
    pub vertices: Vec<[f64; 3]>,
    /// Support point pairs (shape A, shape B) for feature identification.
    pub support_a: Vec<[f64; 3]>,
    /// Support point pairs (shape A, shape B) for feature identification.
    pub support_b: Vec<[f64; 3]>,
}

impl SimplexManager {
    /// Create an empty simplex.
    pub fn new() -> Self {
        SimplexManager {
            vertices: Vec::new(),
            support_a: Vec::new(),
            support_b: Vec::new(),
        }
    }

    /// Number of vertices in the simplex.
    pub fn size(&self) -> usize {
        self.vertices.len()
    }

    /// Add a point to the simplex.
    pub fn add(&mut self, v: [f64; 3], a: [f64; 3], b: [f64; 3]) {
        self.vertices.push(v);
        self.support_a.push(a);
        self.support_b.push(b);
    }

    /// Remove the oldest vertex if simplex has more than 4 vertices.
    pub fn prune(&mut self) {
        while self.vertices.len() > 4 {
            self.vertices.remove(0);
            self.support_a.remove(0);
            self.support_b.remove(0);
        }
    }

    /// Clear the simplex.
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.support_a.clear();
        self.support_b.clear();
    }

    /// Closest point of the newest (at most four) vertices to the origin, its
    /// barycentric weights over that window, and the index of the window's first
    /// vertex.
    fn nearest(&self) -> (SimplexClosest, usize) {
        let start = self.vertices.len().saturating_sub(4);
        (closest_point_to_origin(&self.vertices[start..]), start)
    }

    /// Find the closest point to the origin on the current simplex.
    ///
    /// Returns `(closest_point, did_contain_origin)`. The closest point is the
    /// true nearest point of the convex hull of the vertices (only the newest
    /// four take part, see [`SimplexManager::prune`]); it is the origin itself
    /// when a tetrahedron encloses the origin. `did_contain_origin` is `true`
    /// when the origin lies inside a non-degenerate tetrahedron or on a lower
    /// dimensional simplex.
    pub fn closest_to_origin(&self) -> ([f64; 3], bool) {
        if self.vertices.is_empty() {
            return ([0.0; 3], false);
        }
        let (nearest, _) = self.nearest();
        let contains = nearest.encloses_origin || len_sq3(nearest.point) < 1e-24;
        (nearest.point, contains)
    }

    /// Reduce the simplex to the sub-simplex closest to origin.
    ///
    /// The oldest vertices beyond four are dropped first. Then every vertex that
    /// does not support the closest point (barycentric weight zero) is removed
    /// from `vertices`, `support_a` and `support_b` alike. A tetrahedron that
    /// encloses the origin is kept whole: it is the seed of the EPA polytope.
    pub fn reduce_to_closest(&mut self) {
        self.prune();
        let (nearest, _) = self.nearest();
        if nearest.encloses_origin {
            return;
        }
        for index in (0..self.vertices.len()).rev() {
            if nearest.weights.get(index).is_none_or(|w| *w <= 0.0) {
                self.vertices.remove(index);
                self.support_a.remove(index);
                self.support_b.remove(index);
            }
        }
    }

    /// Points on shape A and shape B that realise [`SimplexManager::closest_to_origin`].
    ///
    /// They are the blends `sum(l_i * support_a_i)` and `sum(l_i * support_b_i)`
    /// with the barycentric weights `l_i` of the closest point, so their
    /// difference is the closest point of the simplex to the origin. For an
    /// enclosing tetrahedron both points coincide and lie in both shapes.
    pub fn witness_points(&self) -> ([f64; 3], [f64; 3]) {
        let (nearest, start) = self.nearest();
        let mut point_a = [0.0; 3];
        let mut point_b = [0.0; 3];
        for (offset, &weight) in nearest.weights.iter().enumerate() {
            let index = start + offset;
            if let (Some(&sa), Some(&sb)) = (self.support_a.get(index), self.support_b.get(index)) {
                point_a = add3(point_a, scale3(sa, weight));
                point_b = add3(point_b, scale3(sb, weight));
            }
        }
        (point_a, point_b)
    }

    /// Simplex made of the given solver vertices.
    fn from_vertices(vertices: &[SimplexVertex]) -> Self {
        let mut simplex = SimplexManager::new();
        for vertex in vertices {
            simplex.add(vertex.w, vertex.sa, vertex.sb);
        }
        simplex
    }
}

impl Default for SimplexManager {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Minkowski difference support function
// ---------------------------------------------------------------------------

/// Support function for Minkowski difference A ⊕ (−B).
///
/// Returns the point in A furthest in direction d, minus the point in B furthest in −d.
pub struct MinkowskiDiffSupport {
    /// Vertices of convex shape A.
    pub shape_a: Vec<[f64; 3]>,
    /// Vertices of convex shape B.
    pub shape_b: Vec<[f64; 3]>,
}

impl MinkowskiDiffSupport {
    /// Construct from vertex lists.
    pub fn new(shape_a: Vec<[f64; 3]>, shape_b: Vec<[f64; 3]>) -> Self {
        MinkowskiDiffSupport { shape_a, shape_b }
    }

    /// Farthest point in `vertices` along direction `d`.
    pub fn support_point(vertices: &[[f64; 3]], d: [f64; 3]) -> [f64; 3] {
        vertices
            .iter()
            .max_by(|&&a, &&b| {
                dot3(a, d)
                    .partial_cmp(&dot3(b, d))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or([0.0; 3])
    }

    /// Minkowski difference support: s_A(d) − s_B(−d).
    pub fn support(&self, d: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let sa = Self::support_point(&self.shape_a, d);
        let sb = Self::support_point(&self.shape_b, neg3(d));
        let v = sub3(sa, sb);
        (v, sa, sb)
    }
}

// ---------------------------------------------------------------------------
// GJK warm-start cache
// ---------------------------------------------------------------------------

/// Cache for warm-starting GJK with the last separating axis.
#[derive(Clone, Debug)]
pub struct GjkCache {
    /// Last separating axis / search direction.
    pub last_dir: [f64; 3],
    /// Whether a valid cached direction exists.
    pub valid: bool,
}

impl GjkCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        GjkCache {
            last_dir: [1.0, 0.0, 0.0],
            valid: false,
        }
    }

    /// Update the cache with a new separating direction.
    pub fn update(&mut self, dir: [f64; 3]) {
        if len3(dir) > 1e-10 {
            self.last_dir = norm3(dir);
            self.valid = true;
        }
    }

    /// Reset the cache.
    pub fn reset(&mut self) {
        self.valid = false;
    }

    /// Initial search direction (from cache or default).
    pub fn initial_dir(&self) -> [f64; 3] {
        if self.valid {
            self.last_dir
        } else {
            [1.0, 0.0, 0.0]
        }
    }
}

impl Default for GjkCache {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EPA face and expanding polytope
// ---------------------------------------------------------------------------

/// A face of the EPA polytope.
#[derive(Clone, Debug)]
pub struct EpaFaceExtended {
    /// Vertex indices.
    pub indices: [usize; 3],
    /// Outward normal (unit).
    pub normal: [f64; 3],
    /// Distance from origin to face plane.
    pub dist: f64,
}

/// Expanding polytope for penetration depth and direction.
pub struct ExpandingPolytope {
    /// Vertices (Minkowski difference points).
    pub vertices: Vec<[f64; 3]>,
    /// Support A points.
    pub support_a: Vec<[f64; 3]>,
    /// Support B points.
    pub support_b: Vec<[f64; 3]>,
    /// Faces.
    pub faces: Vec<EpaFaceExtended>,
}

impl ExpandingPolytope {
    /// Construct from a simplex, which has to be a tetrahedron that encloses
    /// the origin (the simplex [`GjkSolver::query`] ends on for generic
    /// overlapping shapes).
    ///
    /// A simplex with fewer than four vertices (or a degenerate one) spans no
    /// volume: the polytope then has no faces. [`EpaExpander::expand`] still
    /// derives the penetration from such a polytope by completing its vertices
    /// to a tetrahedron of support points of the Minkowski difference.
    pub fn from_simplex(simplex: &SimplexManager) -> Self {
        let vertices = simplex.vertices.clone();
        let support_a = simplex.support_a.clone();
        let support_b = simplex.support_b.clone();

        let mut ep = ExpandingPolytope {
            vertices,
            support_a,
            support_b,
            faces: Vec::new(),
        };
        // Build initial faces of tetrahedron
        if simplex.size() >= 4 {
            let face_indices = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
            for fi in &face_indices {
                ep.add_face(*fi);
            }
        }
        ep
    }

    fn add_face(&mut self, indices: [usize; 3]) {
        let a = self.vertices[indices[0]];
        let b = self.vertices[indices[1]];
        let c = self.vertices[indices[2]];
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let n_raw = cross3(ab, ac);
        let len = len3(n_raw);
        if len < 1e-300 {
            return;
        }
        let normal = scale3(n_raw, 1.0 / len);
        // Ensure normal points away from origin
        let normal = if dot3(normal, a) < 0.0 {
            neg3(normal)
        } else {
            normal
        };
        let dist = dot3(normal, a).abs();
        self.faces.push(EpaFaceExtended {
            indices,
            normal,
            dist,
        });
    }

    /// Penetration result of the part of the boundary nearest to the origin.
    ///
    /// Several faces can lie in the nearest plane (the triangles of one facet of
    /// a polyhedral Minkowski difference), and only one of them contains the
    /// foot of the perpendicular from the origin. Among all faces whose plane is
    /// within a relative `1e-9` of the nearest one, the face whose triangle is
    /// closest to the origin is therefore used: its plane distance is the
    /// depth, its outward normal (reversed, so that it points from B to A) the
    /// normal, and the contact points are the blends of the support points of
    /// its vertices with the barycentric weights of the point of the triangle
    /// nearest to the origin. That point is the foot `p = n * depth` whenever the
    /// foot lies in the triangle, so `point_a - point_b == p == -normal * depth`
    /// then holds exactly, and `point_a` is the point of shape A that has to
    /// travel furthest to leave shape B.
    ///
    /// Returns `None` when the polytope has no faces.
    pub fn nearest_result(&self) -> Option<EpaResult> {
        let nearest_plane = self.closest_face()?.dist;
        let slack = 1e-9 * nearest_plane.max(1.0);
        let mut best: Option<(f64, &EpaFaceExtended, SimplexClosest)> = None;
        for face in self
            .faces
            .iter()
            .filter(|f| f.dist <= nearest_plane + slack)
        {
            let corners = face.indices.map(|index| self.vertices[index]);
            let closest = closest_point_to_origin(&corners);
            let reach = len_sq3(closest.point);
            if best
                .as_ref()
                .is_none_or(|(smallest, _, _)| reach < *smallest)
            {
                best = Some((reach, face, closest));
            }
        }
        let (_, face, closest) = best?;
        let mut point_a = [0.0; 3];
        let mut point_b = [0.0; 3];
        for (&index, &weight) in face.indices.iter().zip(closest.weights.iter()) {
            point_a = add3(point_a, scale3(self.support_a[index], weight));
            point_b = add3(point_b, scale3(self.support_b[index], weight));
        }
        Some(EpaResult {
            depth: face.dist,
            normal: neg3(face.normal),
            point_a,
            point_b,
        })
    }

    /// Find the face closest to the origin.
    pub fn closest_face(&self) -> Option<&EpaFaceExtended> {
        self.faces.iter().min_by(|a, b| {
            a.dist
                .partial_cmp(&b.dist)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Expand the polytope by adding vertex `v` and retriangulating.
    pub fn expand(&mut self, v: [f64; 3], sa: [f64; 3], sb: [f64; 3]) {
        self.vertices.push(v);
        self.support_a.push(sa);
        self.support_b.push(sb);
        let new_idx = self.vertices.len() - 1;

        // Remove faces visible from v
        let visible: Vec<usize> = self
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| {
                let a = self.vertices[f.indices[0]];
                dot3(f.normal, sub3(v, a)) > 0.0
            })
            .map(|(i, _)| i)
            .collect();

        // The horizon is the boundary of the visible region: the edges used by
        // exactly one visible face. Edges are keyed by their unordered vertex
        // pair, so the winding of a face does not matter (`add_face` orients the
        // normal of every new face itself). An edge used by two visible faces is
        // interior to the visible region and disappears with it. A sorted map
        // keeps the order of the new faces deterministic.
        let mut edge_uses: std::collections::BTreeMap<(usize, usize), (usize, [usize; 2])> =
            std::collections::BTreeMap::new();
        for &fi in &visible {
            let f = &self.faces[fi];
            for k in 0..3 {
                let a = f.indices[k];
                let b = f.indices[(k + 1) % 3];
                let entry = edge_uses.entry((a.min(b), a.max(b))).or_insert((0, [a, b]));
                entry.0 += 1;
            }
        }
        let edges: Vec<[usize; 2]> = edge_uses
            .into_values()
            .filter(|&(uses, _)| uses == 1)
            .map(|(_, edge)| edge)
            .collect();

        // Remove visible faces (in reverse order to preserve indices)
        let mut sorted_visible = visible.clone();
        sorted_visible.sort_unstable_by(|a, b| b.cmp(a));
        for &fi in &sorted_visible {
            self.faces.remove(fi);
        }

        // Add new faces
        for e in edges {
            self.add_face([e[0], e[1], new_idx]);
        }
    }
}

// ---------------------------------------------------------------------------
// EPA expander (penetration depth)
// ---------------------------------------------------------------------------

/// Result of EPA: penetration depth and contact normal.
#[derive(Clone, Debug)]
pub struct EpaResult {
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Contact normal (unit, from B to A): the direction in which shape A has to
    /// move to separate from shape B.
    pub normal: [f64; 3],
    /// Contact point on shape A: its deepest point inside shape B.
    pub point_a: [f64; 3],
    /// Contact point on shape B: its deepest point inside shape A.
    ///
    /// `point_a - point_b == -normal * depth`.
    pub point_b: [f64; 3],
}

/// EPA expander: grows the polytope until convergence.
pub struct EpaExpander {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl EpaExpander {
    /// Construct an EPA expander.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        EpaExpander { max_iter, tol }
    }

    /// Run EPA on the given polytope with support function `support`.
    ///
    /// A polytope with faces must enclose the origin (for example the
    /// tetrahedron of support points of the Minkowski difference around it that
    /// [`GjkSolver::query`] ends on for generic overlapping shapes). Each step
    /// expands along the normal of the face nearest to the origin until the
    /// support point lies within `tol` of that face's plane; the depth, normal
    /// and contact points are then those of the nearest part of the boundary
    /// (see [`ExpandingPolytope::nearest_result`]). If `max_iter` steps are not
    /// enough, the nearest part of the polytope reached so far is reported; its
    /// depth is a lower bound of the true penetration depth.
    ///
    /// A polytope without faces comes from a simplex that spans no volume: GJK
    /// ends on a segment or a triangle through the origin when the
    /// configuration is symmetric (two spheres, for instance). Its vertices are
    /// then completed to a tetrahedron of support points of the Minkowski
    /// difference, which does not have to enclose the origin (the expansion
    /// first grows toward it), and the expansion proceeds from there with the
    /// same `max_iter` and `tol`. `polytope` itself is left unchanged in that
    /// case.
    ///
    /// Returns `None` when no penetration can be derived: the Minkowski
    /// difference is flat, so it has no volume to penetrate, or the shapes do
    /// not overlap.
    pub fn expand<F>(&self, polytope: &mut ExpandingPolytope, support: F) -> Option<EpaResult>
    where
        F: Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]),
    {
        if polytope.faces.is_empty() {
            return self.expand_from_vertices(polytope, &support);
        }
        for _ in 0..self.max_iter {
            let face = polytope.closest_face()?.clone();
            let (new_v, sa, sb) = support(face.normal);
            let new_dist = dot3(face.normal, new_v);
            if new_dist - face.dist < self.tol {
                break;
            }
            polytope.expand(new_v, sa, sb);
        }
        polytope.nearest_result()
    }

    /// Penetration of the shapes behind `support`, started from the vertices of
    /// a polytope that has no faces; see [`EpaExpander::expand`].
    fn expand_from_vertices<F>(
        &self,
        polytope: &ExpandingPolytope,
        support: &F,
    ) -> Option<EpaResult>
    where
        F: Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]),
    {
        let oracle = |direction: [f64; 3]| {
            let (w, a, b) = support(direction);
            SupportPoint { w, a, b }
        };
        let known: Vec<SupportPoint> = polytope
            .vertices
            .iter()
            .zip(&polytope.support_a)
            .zip(&polytope.support_b)
            .map(|((&w, &a), &b)| SupportPoint { w, a, b })
            .collect();
        let seed =
            <[SupportPoint; 4]>::try_from(complete_tetrahedron(&known, &oracle).as_slice()).ok()?;
        let found = expand_polytope(seed, &oracle, self.tol, self.max_iter)?;
        Some(EpaResult {
            depth: found.depth,
            normal: neg3(found.normal),
            point_a: found.contact_a,
            point_b: found.contact_b,
        })
    }
}

// ---------------------------------------------------------------------------
// Full GJK solver
// ---------------------------------------------------------------------------

/// Result of GJK distance query.
#[derive(Clone, Debug)]
pub struct GjkResult {
    /// Minimum distance between the two shapes (0 = intersecting).
    pub distance: f64,
    /// Closest point on shape A. For intersecting shapes, a point shared by both
    /// shapes (it coincides with `point_b`).
    pub point_a: [f64; 3],
    /// Closest point on shape B. For intersecting shapes, a point shared by both
    /// shapes (it coincides with `point_a`).
    pub point_b: [f64; 3],
    /// Whether the shapes intersect (overlap or touch).
    pub intersect: bool,
    /// Final simplex. For separated shapes these are the vertices that support
    /// the closest point. For intersecting shapes it is the simplex that proved
    /// the overlap: in general a tetrahedron that encloses the origin (ready to
    /// seed [`ExpandingPolytope::from_simplex`]), but a segment or triangle
    /// through the origin when the configuration is symmetric (for example two
    /// spheres, whose support points stay in one plane through the origin).
    /// [`EpaExpander::expand`] completes such a short simplex from the support
    /// function, so the simplex seeds the penetration query in either case.
    pub simplex: SimplexManager,
}

/// Full GJK solver with simplex management and warm start; `epa` holds the
/// settings of the expansion that turns the final simplex of an overlapping
/// query into a penetration depth.
pub struct GjkSolver {
    /// Maximum GJK iterations.
    pub max_iter: usize,
    /// Relative convergence tolerance: the iteration stops once the distance is
    /// known to within this fraction of itself (`|v|^2 - v.w <= tol |v|^2`).
    pub tol: f64,
    /// Warm-start cache.
    pub cache: GjkCache,
    /// EPA expander.
    pub epa: EpaExpander,
}

impl GjkSolver {
    /// Construct a GJK solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        GjkSolver {
            max_iter,
            tol,
            cache: GjkCache::new(),
            epa: EpaExpander::new(64, 1e-6),
        }
    }

    /// Run GJK between two convex shapes.
    ///
    /// `support(d)` returns `(w, sa, sb)` for the Minkowski difference `A - B`:
    /// the support point `sa` of shape A along `d`, the support point `sb` of
    /// shape B along `-d`, and `w = sa - sb`.
    ///
    /// Each iteration takes the point `v` of the simplex nearest to the origin
    /// together with its barycentric weights, drops the vertices that do not
    /// support `v`, and fetches the support point `w` along `-v`. The iteration
    /// stops when `w` cannot improve on `v`, `|v|^2 - v.w <= tol |v|^2`, or when
    /// `w` repeats a simplex vertex.
    ///
    /// * Separated shapes: `distance` is `|v|`, `point_a` and `point_b` are the
    ///   barycentric blends `sum(l_i * sa_i)` and `sum(l_i * sb_i)` (their
    ///   difference is exactly `v`), and the cache remembers the separating
    ///   direction `-v` as the first search direction of the next query.
    /// * Overlapping or touching shapes: `intersect` is `true`, `distance` is
    ///   `0`, both witnesses are a point shared by the two shapes, `simplex` is
    ///   the final simplex, and the cache is reset.
    ///
    /// If `max_iter` is exhausted first, the best simplex found so far is
    /// reported; its distance is an upper bound on the true distance and is
    /// never `0` unless the shapes really touch.
    pub fn query<F>(&mut self, support: F) -> GjkResult
    where
        F: Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]),
    {
        let vertex_support = |direction: [f64; 3]| -> SimplexVertex {
            let (w, sa, sb) = support(direction);
            SimplexVertex { w, sa, sb }
        };
        let settings = DistanceSettings::new(self.max_iter, self.tol);
        let solution = solve_with(&vertex_support, self.cache.initial_dir(), &settings);
        let simplex = SimplexManager::from_vertices(&solution.simplex);
        match solution.outcome {
            DistanceOutcome::Overlapping { point_a, point_b } => {
                self.cache.reset();
                GjkResult {
                    distance: 0.0,
                    point_a,
                    point_b,
                    intersect: true,
                    simplex,
                }
            }
            DistanceOutcome::Separated {
                distance,
                point_a,
                point_b,
                ..
            } => {
                self.cache.update(sub3(point_b, point_a));
                GjkResult {
                    distance,
                    point_a,
                    point_b,
                    intersect: false,
                    simplex,
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D GJK
// ---------------------------------------------------------------------------

/// Closest points nearer to the origin than this fraction of the simplex extent
/// count as touching (distance zero) in [`Gjk2D::distance`].
const TOUCH_TOLERANCE: f64 = 1e-12;

/// 2-D GJK for polygon-polygon distance.
pub struct Gjk2D {
    /// Maximum iterations.
    pub max_iter: usize,
    /// Relative tolerance of the termination test: the iteration stops once
    /// the distance is known to within this fraction of itself.
    pub tol: f64,
}

impl Gjk2D {
    /// Construct a 2-D GJK solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Gjk2D { max_iter, tol }
    }

    fn support_2d(verts: &[[f64; 2]], d: [f64; 2]) -> [f64; 2] {
        *verts
            .iter()
            .max_by(|&&a, &&b| {
                (a[0] * d[0] + a[1] * d[1])
                    .partial_cmp(&(b[0] * d[0] + b[1] * d[1]))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(&[0.0, 0.0])
    }

    /// Compute the minimum distance between two 2-D convex polygons.
    ///
    /// Runs the Gilbert-Johnson-Keerthi distance iteration on the Minkowski
    /// difference `A - B` of the polygon vertex sets. Each iteration finds the
    /// point `v` of the current simplex (at most three points) nearest to the
    /// origin, drops the simplex vertices that do not support it, and fetches the
    /// support point `w` of `A - B` in direction `-v`. The iteration stops when
    /// `w` cannot improve on `v`, i.e. `|v|^2 - v.w <= tol |v|^2`, which bounds
    /// the relative error of `|v|` by `tol`, or when `w` repeats a simplex vertex.
    ///
    /// Returns `0.0` when the polygons overlap or touch (the simplex encloses the
    /// origin) and `f64::MAX` if either polygon has no vertices. If `max_iter`
    /// is exhausted first, the distance of the best simplex found so far is
    /// returned; it is an upper bound on the true distance and never `0.0` unless
    /// the polygons really touch. At least one iteration is always performed.
    pub fn distance(&self, poly_a: &[[f64; 2]], poly_b: &[[f64; 2]]) -> f64 {
        if poly_a.is_empty() || poly_b.is_empty() {
            return f64::MAX;
        }
        // Support point of A - B, embedded in the z = 0 plane so that the 3-D
        // simplex solver applies unchanged (a triangle that contains the origin
        // then has the origin itself as its nearest point).
        let support = |dir: [f64; 2]| -> [f64; 3] {
            let sa = Self::support_2d(poly_a, dir);
            let sb = Self::support_2d(poly_b, [-dir[0], -dir[1]]);
            [sa[0] - sb[0], sa[1] - sb[1], 0.0]
        };

        let mut simplex: Vec<[f64; 3]> = Vec::with_capacity(3);
        simplex.push(support([1.0, 0.0]));
        let mut v_len = f64::MAX;

        for _ in 0..self.max_iter.max(1) {
            let closest = closest_point_to_origin(&simplex);
            let v = closest.point;
            let v_sq = dot3(v, v);
            let extent_sq = simplex.iter().map(|w| dot3(*w, *w)).fold(0.0_f64, f64::max);
            if v_sq <= TOUCH_TOLERANCE * TOUCH_TOLERANCE * extent_sq {
                return 0.0;
            }
            v_len = v_sq.sqrt();

            // Keep only the vertices that support the closest point.
            let mut index = 0;
            simplex.retain(|_| {
                let keep = closest.weights.get(index).is_some_and(|w| *w > 0.0);
                index += 1;
                keep
            });

            let w = support([-v[0], -v[1]]);
            let gap = v_sq - dot3(v, w);
            if gap <= self.tol * v_sq || simplex.contains(&w) {
                return v_len;
            }
            simplex.push(w);
        }
        v_len
    }
}

// ---------------------------------------------------------------------------
// GJK ray-cast
// ---------------------------------------------------------------------------

/// GJK-based ray against convex shape.
pub struct GjkRaycast {
    /// Maximum iterations.
    pub max_iter: usize,
    /// Tolerance.
    pub tol: f64,
}

/// Result of a GJK ray-cast.
#[derive(Clone, Debug)]
pub struct RaycastResult {
    /// Time of impact along the ray (in \[0, max_t\]), or None.
    pub t: Option<f64>,
    /// Hit point.
    pub point: [f64; 3],
    /// Hit normal (approximate).
    pub normal: [f64; 3],
}

impl GjkRaycast {
    /// Construct a GJK ray-cast solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        GjkRaycast { max_iter, tol }
    }

    /// Cast a ray from `origin` along `dir` against the convex shape described by
    /// its support function.
    ///
    /// `support(d)` must return the point of the shape furthest along `d`, in the
    /// same coordinate frame as the ray. `dir` need not be normalised; `t` and
    /// `max_t` are measured as distances along the normalised direction.
    ///
    /// This is the GJK ray cast of van den Bergen ("Ray casting against general
    /// convex objects with application to continuous collision detection",
    /// 2004). The ray point `x` starts at the origin and advances along the ray
    /// to the supporting planes of the shape: with `v = x - c` for a point `c`
    /// of the shape and `p = support(v)`, `v . (x - p) > 0` means `x` is outside
    /// the supporting plane at `p`, so the ray advances to that plane (or misses
    /// when it moves away from it or passes `max_t`). The support points are
    /// kept as a simplex, `v` becomes the point of `conv(x - P)` nearest to the
    /// origin, and the cast hits once `|v| <= tol`, i.e. `x` lies in the shape.
    ///
    /// Returns `t = Some(distance)` with the hit `point` and the outward surface
    /// `normal` at the hit (the direction of the last supporting plane the ray
    /// point was moved onto: exact for polytopes, accurate to about `sqrt(tol)`
    /// for curved shapes). A ray that starts inside or on the shape hits at
    /// `t = 0` and no surface was crossed, so the normal is the negated ray
    /// direction. The cast reports a miss if `max_iter` iterations are not
    /// enough to bring `|v|` below `tol`, if the ray points away from the
    /// shape, or if the hit lies beyond `max_t`.
    pub fn cast<F>(&self, origin: [f64; 3], dir: [f64; 3], max_t: f64, support: F) -> RaycastResult
    where
        F: Fn([f64; 3]) -> [f64; 3],
    {
        let miss = RaycastResult {
            t: None,
            point: [0.0; 3],
            normal: [0.0; 3],
        };
        let length = len3(dir);
        if !length.is_finite() || length <= 0.0 || max_t.is_nan() || max_t < 0.0 {
            return miss;
        }
        let ray = scale3(dir, 1.0 / length);
        let tol_sq = self.tol * self.tol;

        let mut lambda = 0.0_f64;
        let mut x = origin;
        let mut hit_normal = [0.0; 3];
        // Support points of the shape; the simplex is `x - P`, which moves with `x`.
        let mut points: Vec<[f64; 3]> = Vec::with_capacity(4);
        // Any point of the shape seeds the search vector.
        let mut v = sub3(x, support(ray));

        for _ in 0..self.max_iter {
            let v_sq = dot3(v, v);
            if v_sq <= tol_sq {
                break;
            }
            let p = support(scale3(v, 1.0 / v_sq.sqrt()));
            let w = sub3(x, p);
            let v_dot_w = dot3(v, w);
            if v_dot_w > 0.0 {
                // `x` is outside the supporting plane at `p`.
                let v_dot_r = dot3(v, ray);
                if v_dot_r >= 0.0 {
                    return miss;
                }
                lambda -= v_dot_w / v_dot_r;
                if lambda > max_t {
                    return miss;
                }
                x = add3(origin, scale3(ray, lambda));
                hit_normal = v;
            }
            // A repeated support point adds no information (and is not a sign of
            // convergence: the ray point may just have moved onto its plane).
            if !points.contains(&p) {
                points.push(p);
            }
            let count = points.len().min(4);
            let mut simplex = [[0.0_f64; 3]; 4];
            for (slot, q) in simplex.iter_mut().zip(points.iter()) {
                *slot = sub3(x, *q);
            }
            let closest = closest_point_to_origin(&simplex[..count]);
            if closest.encloses_origin {
                // The ray point lies inside conv(P), hence inside the shape.
                v = [0.0; 3];
                break;
            }
            v = closest.point;
            let mut index = 0;
            points.retain(|_| {
                let keep = closest.weights.get(index).is_some_and(|w| *w > 0.0);
                index += 1;
                keep
            });
        }

        if dot3(v, v) > tol_sq {
            return miss;
        }
        let normal = if dot3(hit_normal, hit_normal) > 0.0 {
            norm3(hit_normal)
        } else {
            neg3(ray)
        };
        RaycastResult {
            t: Some(lambda),
            point: x,
            normal,
        }
    }
}

// ---------------------------------------------------------------------------
// GJK Time of Impact (Conservative Advancement)
// ---------------------------------------------------------------------------

/// Closing speeds (along the separating axis) at or below this value never
/// close the gap: the shapes are treated as separating.
const MIN_CLOSING_SPEED: f64 = 1e-12;

/// Witness distances below this fraction of the (unit-offset) coordinate scale
/// are rounding noise: their direction is not used as a contact normal.
const NORMAL_NOISE_FLOOR: f64 = 1e-9;

/// Result of a TOI computation.
#[derive(Clone, Debug)]
pub struct ToiGjkResult {
    /// Time of impact ∈ \[0, dt\], or None if no impact.
    pub toi: Option<f64>,
    /// Contact point at TOI: the point of shape A closest to shape B (zero when
    /// there is no impact).
    pub contact_point: [f64; 3],
    /// Unit contact normal at TOI, pointing from shape A toward shape B (zero
    /// when there is no impact).
    ///
    /// Shapes that already overlap at the start have no separating direction;
    /// their normal is the direction opposite to the relative velocity of B
    /// (the way B has to be pushed back), or zero if that velocity is zero.
    pub contact_normal: [f64; 3],
}

/// GJK time-of-impact solver using conservative advancement.
pub struct GjkTimeOfImpact {
    /// Maximum number of advancement steps (one GJK distance query each).
    pub max_iter: usize,
    /// Distance at which the shapes count as touching (absolute length).
    pub tol: f64,
}

impl GjkTimeOfImpact {
    /// Construct a TOI solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        GjkTimeOfImpact { max_iter, tol }
    }

    /// Compute TOI for two shapes moving linearly.
    ///
    /// `support(d, t)` returns the support triple `(w, sa, sb)` of the Minkowski
    /// difference `A - B` along `d` with both shapes placed where they are at
    /// time `t` (`0 <= t <= dt`), exactly like the closure of
    /// [`GjkSolver::query`]. `vel_rel` is the velocity of B relative to A
    /// (`v_B - v_A`) and must be the constant relative velocity of the motion
    /// that `support` encodes; the advancement is conservative only for such a
    /// linear relative motion.
    ///
    /// Conservative advancement: starting at `t = 0`, GJK measures the distance
    /// `d` between the shapes at `t`, with the closest points `pa` on A and `pb`
    /// on B and the unit separating axis `n = (pb - pa) / d`. Along `n` the
    /// shapes approach at the speed `mu = -vel_rel . n`, so they cannot touch
    /// before `t + d / mu` (the certified lower bound of `d` is used, which never
    /// overshoots). The loop stops with
    ///
    /// * `Some(t)` when `d <= tol` (or the shapes overlap), reporting `pa` and
    ///   `n`; shapes that overlap at `t = 0` hit at `0`,
    /// * `None` when `mu <= 1e-12` (the shapes never approach along `n`, so
    ///   they stay apart) or the advanced time passes `dt`, or `dt` or
    ///   `vel_rel` is not finite,
    /// * `Some(t)` with the last conservative time when `max_iter` steps do not
    ///   suffice to get within `tol`; the shapes cannot touch before that time.
    pub fn compute<F>(&self, support: F, vel_rel: [f64; 3], dt: f64) -> ToiGjkResult
    where
        F: Fn([f64; 3], f64) -> ([f64; 3], [f64; 3], [f64; 3]),
    {
        let no_impact = ToiGjkResult {
            toi: None,
            contact_point: [0.0; 3],
            contact_normal: [0.0; 3],
        };
        if !dt.is_finite() || dt < 0.0 || !vel_rel.iter().all(|c| c.is_finite()) {
            return no_impact;
        }
        let tolerance = if self.tol.is_finite() {
            self.tol.max(0.0)
        } else {
            0.0
        };
        let settings = DistanceSettings::default();

        let mut t = 0.0_f64;
        let mut direction = [1.0, 0.0, 0.0];
        let mut normal = neg3(norm3(vel_rel));
        let mut have_normal = false;
        let mut contact = [0.0; 3];
        for _ in 0..self.max_iter.max(1) {
            let at_time = |d: [f64; 3]| -> SimplexVertex {
                let (w, sa, sb) = support(d, t);
                SimplexVertex { w, sa, sb }
            };
            let solution = solve_with(&at_time, direction, &settings);
            match solution.outcome {
                DistanceOutcome::Overlapping { point_a, .. } => {
                    return ToiGjkResult {
                        toi: Some(t),
                        contact_point: point_a,
                        contact_normal: normal,
                    };
                }
                DistanceOutcome::Separated {
                    distance,
                    lower_bound,
                    point_a,
                    point_b,
                } => {
                    // The direction of a witness difference that has shrunk to
                    // rounding noise is meaningless; the previous step's normal
                    // (the shapes were approaching along it) stands.
                    if !have_normal || distance > NORMAL_NOISE_FLOOR * (1.0 + len3(point_a)) {
                        normal = scale3(sub3(point_b, point_a), 1.0 / distance);
                        have_normal = true;
                    }
                    contact = point_a;
                    if distance <= tolerance {
                        return ToiGjkResult {
                            toi: Some(t),
                            contact_point: contact,
                            contact_normal: normal,
                        };
                    }
                    let closing_speed = -dot3(vel_rel, normal);
                    if closing_speed <= MIN_CLOSING_SPEED {
                        return no_impact;
                    }
                    t += lower_bound / closing_speed;
                    if !t.is_finite() || t > dt {
                        return no_impact;
                    }
                    direction = normal;
                }
            }
        }
        ToiGjkResult {
            toi: Some(t),
            contact_point: contact,
            contact_normal: normal,
        }
    }
}

// ---------------------------------------------------------------------------
// Feature pair
// ---------------------------------------------------------------------------

/// Contact feature type.
#[derive(Clone, Debug, PartialEq)]
pub enum FeaturePairType {
    /// Vertex-vertex contact.
    VertexVertex,
    /// Vertex-face contact.
    VertexFace,
    /// Edge-edge contact.
    EdgeEdge,
    /// Face-face contact.
    FaceFace,
}

/// A feature pair identifying the closest features between two shapes.
#[derive(Clone, Debug)]
pub struct FeaturePair {
    /// Feature type.
    pub feature_type: FeaturePairType,
    /// Index of feature on shape A.
    pub idx_a: usize,
    /// Index of feature on shape B.
    pub idx_b: usize,
    /// Contact point (world coordinates).
    pub contact_point: [f64; 3],
    /// Contact normal (from A to B).
    pub contact_normal: [f64; 3],
}

impl FeaturePair {
    /// Create a new feature pair.
    pub fn new(
        feature_type: FeaturePairType,
        idx_a: usize,
        idx_b: usize,
        contact_point: [f64; 3],
        contact_normal: [f64; 3],
    ) -> Self {
        FeaturePair {
            feature_type,
            idx_a,
            idx_b,
            contact_point,
            contact_normal,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    #[test]
    fn test_sub3() {
        let v = sub3([3.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((v[0] - 2.0).abs() < EPS);
    }

    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < EPS);
    }

    #[test]
    fn test_cross3_orthogonal() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < EPS);
        assert!(c[0].abs() < EPS);
        assert!(c[1].abs() < EPS);
    }

    #[test]
    fn test_norm3_unit() {
        let v = [3.0, 4.0, 0.0];
        let n = norm3(v);
        assert!((len3(n) - 1.0).abs() < EPS);
    }

    #[test]
    fn test_line_closest_to_origin_midpoint() {
        // Line from (-1,0,0) to (1,0,0): closest to origin is (0,0,0)
        let (t, p) = line_closest_point_to_origin([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((t - 0.5).abs() < EPS);
        assert!(len3(p) < EPS);
    }

    #[test]
    fn test_line_closest_to_origin_endpoint() {
        // Line from (1,0,0) to (2,0,0): closest to origin is (1,0,0)
        let (t, p) = line_closest_point_to_origin([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((t - 0.0).abs() < EPS);
        assert!((p[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_triangle_closest_to_origin_inside() {
        // Unit triangle containing origin projection
        let a = [-1.0, -1.0, 0.0];
        let b = [1.0, -1.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = triangle_closest_to_origin(a, b, c);
        assert!(p[2].abs() < EPS);
    }

    #[test]
    fn test_tetrahedron_contains_origin_true() {
        // Regular tetrahedron centered at origin
        let a = [1.0, 0.0, -1.0 / 2.0_f64.sqrt()];
        let b = [-1.0, 0.0, -1.0 / 2.0_f64.sqrt()];
        let c = [0.0, 1.0, 1.0 / 2.0_f64.sqrt()];
        let d = [0.0, -1.0, 1.0 / 2.0_f64.sqrt()];
        let _ = tetrahedron_contains_origin(a, b, c, d);
        // Just verify it runs without panic
    }

    #[test]
    fn test_simplex_manager_add_and_size() {
        let mut s = SimplexManager::new();
        s.add([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(s.size(), 1);
    }

    #[test]
    fn test_simplex_manager_closest_1d() {
        let mut s = SimplexManager::new();
        s.add([0.5, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let (p, contains) = s.closest_to_origin();
        assert!(!contains);
        assert!((p[0] - 0.5).abs() < EPS);
    }

    #[test]
    fn test_simplex_manager_origin_inside_1d() {
        let mut s = SimplexManager::new();
        s.add([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let (_, contains) = s.closest_to_origin();
        assert!(contains);
    }

    #[test]
    fn test_minkowski_diff_support() {
        let shape_a = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let shape_b = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let mds = MinkowskiDiffSupport::new(shape_a, shape_b);
        let (v, _, _) = mds.support([1.0, 0.0, 0.0]);
        // s_A([1,0,0]) = [1,0,0]; s_B([-1,0,0]) = [-0.5,0,0]; diff = [1.5,0,0]
        assert!((v[0] - 1.5).abs() < EPS);
    }

    #[test]
    fn test_gjk_cache_update() {
        let mut cache = GjkCache::new();
        cache.update([1.0, 0.0, 0.0]);
        assert!(cache.valid);
        assert!((cache.last_dir[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_gjk_cache_reset() {
        let mut cache = GjkCache::new();
        cache.update([1.0, 0.0, 0.0]);
        cache.reset();
        assert!(!cache.valid);
    }

    #[test]
    fn test_feature_pair_creation() {
        let fp = FeaturePair::new(FeaturePairType::EdgeEdge, 0, 1, [0.0; 3], [0.0, 0.0, 1.0]);
        assert_eq!(fp.feature_type, FeaturePairType::EdgeEdge);
    }

    #[test]
    fn test_expanding_polytope_from_simplex() {
        // Build a valid tetrahedron simplex
        let mut s = SimplexManager::new();
        s.add([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]);
        s.add([-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3]);
        s.add([0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0; 3]);
        s.add([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0; 3]);
        let ep = ExpandingPolytope::from_simplex(&s);
        assert!(!ep.faces.is_empty());
    }

    #[test]
    fn test_simplex_manager_prune() {
        let mut s = SimplexManager::new();
        for i in 0..6 {
            s.add([i as f64, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        }
        s.prune();
        assert!(s.size() <= 4);
    }
}

#[cfg(test)]
mod regression_tests;

#[cfg(test)]
mod solver_tests;
