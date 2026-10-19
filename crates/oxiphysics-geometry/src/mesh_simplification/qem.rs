//! Quadric error metric mesh simplification (Garland & Heckbert, 1997).
//!
//! The simplifier collapses edges in order of increasing quadric error. Every
//! vertex carries the sum of the squared-distance quadrics of the planes of its
//! incident faces; collapsing the edge `(a, b)` places the merged vertex at the
//! minimiser of the summed quadric, and the value of the quadric there is the cost
//! of the collapse.
//!
//! Bookkeeping that keeps the result a valid mesh:
//!
//! - The surface is kept as vertex-to-face adjacency lists that are updated in
//!   place, so a collapse touches only the faces around the two vertices.
//! - Every vertex has a version number that grows when the vertex survives a
//!   collapse. Heap entries remember the versions of their end points and are
//!   dropped when they are popped after either end point has changed or is gone,
//!   so a position or cost computed for an old neighbourhood is never applied.
//! - A collapse is only performed when the *link condition* holds (the end points
//!   share exactly the vertices opposite to the edge, with the virtual vertex of a
//!   boundary counted for boundary vertices), when the merged vertex keeps a valid
//!   fan of faces, and when none of the surrounding faces folds over, becomes
//!   degenerate or turns by more than `max_normal_deviation`.
//! - A rejected edge is remembered and tried again as soon as a collapse changes
//!   the neighbourhood of one of its end points, so the simplifier only stops
//!   short of the target when no valid collapse is left.

use std::collections::BinaryHeap;

use super::functions::{
    compute_vertex_quadrics, corner_move_keeps_normal, cross3, dist3, dot3, mid3, normalize3, sub3,
};
use super::types::{CollapseCandidate, QemConfig, SymMat4, TextureSeamFlags};

/// Largest distance between the optimal position of a collapse and the middle of
/// the edge, in edge lengths. A minimiser of an almost singular quadric can lie
/// far away from the edge; such a position is not trusted.
const MAX_PLACEMENT_EDGES: f64 = 2.0;

/// Smallest cosine between the old and the new normal of a face that is accepted
/// whatever `max_normal_deviation` says: a face whose normal turns by (almost) a
/// right angle or more is folded over.
const FLIP_COSINE: f64 = 1.0e-9;

/// Outcome of checking a popped collapse against the current mesh.
enum Verdict {
    /// The collapse can be performed.
    Valid,
    /// The collapse is not possible now but may become possible when the
    /// neighbourhood of the edge changes.
    Retry,
    /// The collapse is not possible and never will be.
    Never,
}

/// Perform QEM-based mesh simplification.
///
/// Returns simplified vertices and triangles. The triangle orientation of the
/// input is kept, vertices that no triangle uses are dropped, and triangles that
/// reference a missing vertex, repeat a vertex or have a non-finite corner are
/// ignored.
///
/// Edges are collapsed in order of increasing quadric error until the number of
/// triangles drops to `config.target_triangles` (when it is non-zero), until the
/// cheapest remaining collapse costs more than `config.max_error` (when it is
/// positive), or until no valid collapse is left. A collapse is valid when it
/// keeps the surface a manifold (link condition), does not fold over or
/// degenerate a face, and, when `config.max_normal_deviation` is positive, turns
/// no face by more than that angle. A closed surface is never reduced below a
/// tetrahedron.
///
/// - `config.preserve_boundary`: the boundary is kept exactly. Boundary vertices
///   never move and are never merged with each other; a boundary vertex only
///   absorbs interior vertices. When it is `false`, the boundary may be
///   simplified, and `config.boundary_weight` scales the quadrics of the planes
///   through each boundary edge perpendicular to its face, which keep the
///   boundary near its original shape.
/// - `config.preserve_texture_seams` has an effect through
///   [`qem_simplify_with_seams`] only; this function knows no seams.
pub fn qem_simplify(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    config: &QemConfig,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    qem_simplify_with_seams(vertices, triangles, config, None)
}

/// Perform QEM-based mesh simplification that keeps texture seams.
///
/// Behaves like [`qem_simplify`]. When `config.preserve_texture_seams` is set and
/// `seams` is given, every vertex of a seam edge is kept in place and seam
/// vertices are never merged with each other, so the seam edges survive
/// unchanged; a seam vertex only absorbs vertices that are not on a seam. Seam
/// edges that are not edges of the mesh are ignored.
pub fn qem_simplify_with_seams(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    config: &QemConfig,
    seams: Option<&TextureSeamFlags>,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut simplifier = Simplifier::new(vertices, triangles, config, seams);
    simplifier.run();
    simplifier.into_mesh()
}

/// State of one simplification run.
struct Simplifier<'a> {
    config: &'a QemConfig,
    /// Smallest accepted cosine between the old and the new normal of a face.
    min_cosine: f64,
    position: Vec<[f64; 3]>,
    quadric: Vec<SymMat4>,
    /// Grows each time the vertex survives a collapse.
    version: Vec<u32>,
    vertex_alive: Vec<bool>,
    /// The vertex lies on the boundary of the surface.
    boundary: Vec<bool>,
    /// The vertex is an end point of a texture seam that is preserved.
    seam_vertex: Vec<bool>,
    /// The neighbourhood of the vertex is not a manifold; it is left alone.
    frozen: Vec<bool>,
    /// Live faces around each vertex.
    vertex_faces: Vec<Vec<usize>>,
    /// For each vertex the other end points of edges that were rejected and may
    /// become collapsible when the neighbourhood changes.
    deferred: Vec<Vec<usize>>,
    faces: Vec<[usize; 3]>,
    face_alive: Vec<bool>,
    live_faces: usize,
    heap: BinaryHeap<CollapseCandidate>,
}

impl<'a> Simplifier<'a> {
    fn new(
        vertices: &[[f64; 3]],
        triangles: &[[usize; 3]],
        config: &'a QemConfig,
        seams: Option<&TextureSeamFlags>,
    ) -> Self {
        let vertex_count = vertices.len();
        let faces: Vec<[usize; 3]> = triangles
            .iter()
            .copied()
            .filter(|tri| {
                tri.iter().all(|&v| v < vertex_count)
                    && tri[0] != tri[1]
                    && tri[1] != tri[2]
                    && tri[0] != tri[2]
                    && tri
                        .iter()
                        .all(|&v| vertices[v].iter().all(|x| x.is_finite()))
            })
            .collect();

        let mut vertex_faces = vec![Vec::new(); vertex_count];
        for (fi, tri) in faces.iter().enumerate() {
            for &v in tri {
                vertex_faces[v].push(fi);
            }
        }

        let deviation_cosine = if config.max_normal_deviation > 0.0 {
            config.max_normal_deviation.cos()
        } else {
            -1.0
        };
        let mut simplifier = Self {
            config,
            min_cosine: deviation_cosine.max(FLIP_COSINE),
            position: vertices.to_vec(),
            quadric: compute_vertex_quadrics(vertices, &faces),
            version: vec![0; vertex_count],
            vertex_alive: vec![true; vertex_count],
            boundary: vec![false; vertex_count],
            seam_vertex: vec![false; vertex_count],
            frozen: vec![false; vertex_count],
            vertex_faces,
            deferred: vec![Vec::new(); vertex_count],
            live_faces: faces.len(),
            face_alive: vec![true; faces.len()],
            faces,
            heap: BinaryHeap::new(),
        };
        simplifier.classify_vertices();
        if !config.preserve_boundary && config.boundary_weight > 0.0 {
            simplifier.add_boundary_quadrics(config.boundary_weight);
        }
        if config.preserve_texture_seams
            && let Some(seams) = seams
        {
            simplifier.mark_seams(seams);
        }
        simplifier
    }

    /// Find the boundary vertices and the vertices with a non-manifold
    /// neighbourhood.
    fn classify_vertices(&mut self) {
        for v in 0..self.position.len() {
            let faces = &self.vertex_faces[v];
            if faces.is_empty() {
                continue;
            }
            let mut others: Vec<usize> = faces
                .iter()
                .flat_map(|&f| self.faces[f].iter().copied().filter(move |&u| u != v))
                .collect();
            others.sort_unstable();
            let mut boundary_edges = 0;
            let mut non_manifold = false;
            let mut i = 0;
            while i < others.len() {
                let run = others[i..].iter().take_while(|&&u| u == others[i]).count();
                match run {
                    1 => boundary_edges += 1,
                    2 => {}
                    _ => non_manifold = true,
                }
                i += run;
            }
            // The faces around a manifold vertex form one closed or one open fan.
            if boundary_edges != 0 && boundary_edges != 2 {
                non_manifold = true;
            }
            if !non_manifold && !self.faces_form_one_fan(v) {
                non_manifold = true;
            }
            self.boundary[v] = boundary_edges > 0;
            self.frozen[v] = non_manifold;
        }
    }

    /// Whether the faces around `v` are connected through edges at `v`.
    fn faces_form_one_fan(&self, v: usize) -> bool {
        let faces = &self.vertex_faces[v];
        if faces.len() <= 1 {
            return true;
        }
        // Faces that share a neighbour of `v` share an edge at `v`.
        let mut by_neighbour: Vec<(usize, usize)> = Vec::with_capacity(2 * faces.len());
        for (local, &f) in faces.iter().enumerate() {
            for &u in &self.faces[f] {
                if u != v {
                    by_neighbour.push((u, local));
                }
            }
        }
        by_neighbour.sort_unstable();
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn root(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for pair in by_neighbour.windows(2) {
            if pair[0].0 == pair[1].0 {
                let (a, b) = (root(&mut parent, pair[0].1), root(&mut parent, pair[1].1));
                parent[a] = b;
            }
        }
        let first = root(&mut parent, 0);
        (1..faces.len()).all(|local| root(&mut parent, local) == first)
    }

    /// Add the quadric of the plane through every boundary edge and perpendicular
    /// to its face, scaled by `weight`, to both end points of the edge.
    fn add_boundary_quadrics(&mut self, weight: f64) {
        let mut constraints: Vec<(usize, usize, SymMat4)> = Vec::new();
        for tri in &self.faces {
            let normal = normalize3(cross3(
                sub3(self.position[tri[1]], self.position[tri[0]]),
                sub3(self.position[tri[2]], self.position[tri[0]]),
            ));
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                if self.faces_with_both(a, b).len() != 1 {
                    continue;
                }
                let plane = normalize3(cross3(normal, sub3(self.position[b], self.position[a])));
                let offset = -dot3(plane, self.position[a]);
                let constraint =
                    SymMat4::from_plane(plane[0], plane[1], plane[2], offset).scale(weight);
                constraints.push((a, b, constraint));
            }
        }
        for (a, b, constraint) in constraints {
            self.quadric[a] = self.quadric[a].add(&constraint);
            self.quadric[b] = self.quadric[b].add(&constraint);
        }
    }

    /// Mark the end points of the seam edges that exist in the mesh.
    fn mark_seams(&mut self, seams: &TextureSeamFlags) {
        for &(a, b) in &seams.seam_edges {
            if a < self.position.len()
                && b < self.position.len()
                && !self.faces_with_both(a, b).is_empty()
            {
                self.seam_vertex[a] = true;
                self.seam_vertex[b] = true;
            }
        }
    }

    /// Live faces that contain both `a` and `b`.
    fn faces_with_both(&self, a: usize, b: usize) -> Vec<usize> {
        self.vertex_faces[a]
            .iter()
            .copied()
            .filter(|&f| self.faces[f].contains(&b))
            .collect()
    }

    /// Sorted vertices that share a live face with `v`.
    fn neighbours(&self, v: usize) -> Vec<usize> {
        let mut result: Vec<usize> = self.vertex_faces[v]
            .iter()
            .flat_map(|&f| self.faces[f].iter().copied().filter(move |&u| u != v))
            .collect();
        result.sort_unstable();
        result.dedup();
        result
    }

    /// A pinned vertex must stay where it is.
    fn is_pinned(&self, v: usize) -> bool {
        (self.config.preserve_boundary && self.boundary[v]) || self.seam_vertex[v]
    }

    /// Position of the merged vertex of the collapse of the edge `(a, b)`: the
    /// minimiser of the summed quadric if it is trustworthy, else the best of the
    /// middle of the edge and its two end points.
    fn placement(&self, a: usize, b: usize, quadric: &SymMat4) -> [f64; 3] {
        let (pa, pb) = (self.position[a], self.position[b]);
        let middle = mid3(pa, pb);
        let reach = MAX_PLACEMENT_EDGES * dist3(pa, pb);
        if let Some(optimum) = quadric.optimal_point()
            && optimum.iter().all(|x| x.is_finite())
            && dist3(optimum, middle) <= reach
        {
            return optimum;
        }
        let mut best = middle;
        let mut best_error = quadric.evaluate(middle);
        for point in [pa, pb] {
            let error = quadric.evaluate(point);
            if error < best_error {
                best = point;
                best_error = error;
            }
        }
        best
    }

    /// The heap entry for collapsing the edge `(a, b)`, if the edge may ever be
    /// collapsed.
    fn candidate(&self, a: usize, b: usize) -> Option<CollapseCandidate> {
        if a == b
            || !self.vertex_alive[a]
            || !self.vertex_alive[b]
            || self.frozen[a]
            || self.frozen[b]
        {
            return None;
        }
        let (pinned_a, pinned_b) = (self.is_pinned(a), self.is_pinned(b));
        if pinned_a && pinned_b {
            return None;
        }
        let (survivor, removed) = if pinned_a {
            (a, b)
        } else if pinned_b {
            (b, a)
        } else {
            (a.min(b), a.max(b))
        };
        let quadric = self.quadric[survivor].add(&self.quadric[removed]);
        let optimal_pos = if pinned_a || pinned_b {
            self.position[survivor]
        } else {
            self.placement(survivor, removed, &quadric)
        };
        let cost = quadric.evaluate(optimal_pos);
        if !cost.is_finite() {
            return None;
        }
        Some(CollapseCandidate {
            v0: survivor,
            v1: removed,
            optimal_pos,
            cost: cost.max(0.0),
            v0_version: self.version[survivor],
            v1_version: self.version[removed],
        })
    }

    fn push_candidate(&mut self, a: usize, b: usize) {
        if let Some(candidate) = self.candidate(a, b) {
            self.heap.push(candidate);
        }
    }

    /// Whether the end points of the candidate are unchanged since it was made.
    fn is_current(&self, candidate: &CollapseCandidate) -> bool {
        self.vertex_alive[candidate.v0]
            && self.vertex_alive[candidate.v1]
            && self.version[candidate.v0] == candidate.v0_version
            && self.version[candidate.v1] == candidate.v1_version
    }

    /// Check a candidate against the current mesh.
    fn assess(&self, candidate: &CollapseCandidate) -> Verdict {
        let (survivor, removed) = (candidate.v0, candidate.v1);
        let shared = self.faces_with_both(survivor, removed);
        if shared.is_empty() || shared.len() > 2 {
            return Verdict::Never;
        }

        // Link condition: the end points share exactly the vertices opposite to
        // the edge; a boundary counts as one more shared (virtual) vertex.
        let around_survivor = self.neighbours(survivor);
        let around_removed = self.neighbours(removed);
        let common = around_survivor
            .iter()
            .filter(|u| around_removed.binary_search(u).is_ok())
            .count();
        let virtual_vertex = usize::from(self.boundary[survivor] && self.boundary[removed]);
        if common + virtual_vertex != 2 {
            return Verdict::Retry;
        }

        // The merged vertex needs a fan of faces: a closed one of at least three
        // faces, or an open one of at least one face.
        let merged_on_boundary = self.boundary[survivor] || self.boundary[removed];
        let remaining =
            self.vertex_faces[survivor].len() + self.vertex_faces[removed].len() - 2 * shared.len();
        let needed = if merged_on_boundary { 1 } else { 3 };
        if remaining < needed {
            return Verdict::Retry;
        }

        // No surrounding face may fold over, degenerate or turn too far.
        let target = candidate.optimal_pos;
        for &f in self.vertex_faces[survivor]
            .iter()
            .chain(&self.vertex_faces[removed])
        {
            if shared.contains(&f) {
                continue;
            }
            let tri = self.faces[f];
            let before = tri.map(|v| self.position[v]);
            let after = tri.map(|v| {
                if v == survivor || v == removed {
                    target
                } else {
                    self.position[v]
                }
            });
            if !corner_move_keeps_normal(before, after, self.min_cosine) {
                return Verdict::Retry;
            }
        }
        Verdict::Valid
    }

    /// Remember a rejected edge so that it is tried again when the neighbourhood
    /// of one of its end points changes.
    fn defer(&mut self, a: usize, b: usize) {
        self.deferred[a].push(b);
        self.deferred[b].push(a);
    }

    /// Merge `removed` into `survivor`, which moves to `position`.
    fn collapse(&mut self, survivor: usize, removed: usize, position: [f64; 3]) {
        let shared = self.faces_with_both(survivor, removed);
        for &f in &shared {
            self.face_alive[f] = false;
            if let Some(third) = self.faces[f]
                .iter()
                .copied()
                .find(|&v| v != survivor && v != removed)
            {
                self.vertex_faces[third].retain(|&g| g != f);
            }
        }
        self.live_faces -= shared.len();

        for f in std::mem::take(&mut self.vertex_faces[removed]) {
            if shared.contains(&f) {
                continue;
            }
            for v in self.faces[f].iter_mut() {
                if *v == removed {
                    *v = survivor;
                }
            }
            self.vertex_faces[survivor].push(f);
        }
        let face_alive = &self.face_alive;
        self.vertex_faces[survivor].retain(|&f| face_alive[f]);

        self.position[survivor] = position;
        self.quadric[survivor] = self.quadric[survivor].add(&self.quadric[removed]);
        self.boundary[survivor] = self.boundary[survivor] || self.boundary[removed];
        self.seam_vertex[survivor] = self.seam_vertex[survivor] || self.seam_vertex[removed];
        self.vertex_alive[removed] = false;
        self.version[survivor] += 1;
        self.deferred[removed].clear();
    }

    /// Queue the collapses that the collapse into `survivor` made new or possible.
    fn refresh_around(&mut self, survivor: usize) {
        // Every edge of the survivor has a new cost.
        let ring = self.neighbours(survivor);
        for &u in &ring {
            self.push_candidate(survivor, u);
        }
        // Whether an edge may be collapsed depends on the faces and neighbours of
        // its end points only, so a rejected edge can have become collapsible only
        // if one of its end points is the survivor or one of its neighbours.
        let mut retry: Vec<(usize, usize)> = Vec::new();
        for x in std::iter::once(survivor).chain(ring.iter().copied()) {
            for y in std::mem::take(&mut self.deferred[x]) {
                if self.vertex_alive[y] {
                    retry.push((x.min(y), x.max(y)));
                }
            }
        }
        retry.sort_unstable();
        retry.dedup();
        for (a, b) in retry {
            self.deferred[a].retain(|&z| z != b);
            self.deferred[b].retain(|&z| z != a);
            if a == survivor || b == survivor || self.faces_with_both(a, b).is_empty() {
                continue;
            }
            self.push_candidate(a, b);
        }
    }

    /// Collapse edges until the target is met or no valid collapse is left.
    fn run(&mut self) {
        let mut edges: Vec<(usize, usize)> = Vec::with_capacity(3 * self.faces.len());
        for tri in &self.faces {
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                edges.push((a.min(b), a.max(b)));
            }
        }
        edges.sort_unstable();
        edges.dedup();
        for (a, b) in edges {
            self.push_candidate(a, b);
        }

        let target = self.config.target_triangles;
        while let Some(candidate) = self.heap.pop() {
            if target > 0 && self.live_faces <= target {
                break;
            }
            if !self.is_current(&candidate) {
                continue;
            }
            if self.config.max_error > 0.0 && candidate.cost > self.config.max_error {
                break;
            }
            match self.assess(&candidate) {
                Verdict::Valid => {
                    self.collapse(candidate.v0, candidate.v1, candidate.optimal_pos);
                    self.refresh_around(candidate.v0);
                }
                Verdict::Retry => self.defer(candidate.v0, candidate.v1),
                Verdict::Never => {}
            }
        }
    }

    /// The simplified mesh: live faces and the vertices they use.
    fn into_mesh(self) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let mut used = vec![false; self.position.len()];
        for (fi, tri) in self.faces.iter().enumerate() {
            if self.face_alive[fi] {
                for &v in tri {
                    used[v] = true;
                }
            }
        }
        let mut remap = vec![usize::MAX; self.position.len()];
        let mut vertices = Vec::new();
        for (v, &position) in self.position.iter().enumerate() {
            if self.vertex_alive[v] && used[v] {
                remap[v] = vertices.len();
                vertices.push(position);
            }
        }
        let triangles = self
            .faces
            .iter()
            .enumerate()
            .filter(|&(fi, _)| self.face_alive[fi])
            .map(|(_, tri)| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
            .collect();
        (vertices, triangles)
    }
}

#[cfg(test)]
mod invariant_tests;
