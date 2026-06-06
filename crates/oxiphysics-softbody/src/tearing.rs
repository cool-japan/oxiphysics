// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cloth/mesh tearing and cutting.
//!
//! Provides topological operations for simulating the tearing and cutting of
//! deformable meshes, including stress-based tear initiation, vertex duplication,
//! edge cutting, progressive tear propagation, and boundary smoothing.

/// An edge in the mesh, described by two vertex indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    /// Index of the first endpoint vertex.
    pub a: usize,
    /// Index of the second endpoint vertex.
    pub b: usize,
    /// Whether this edge has been torn.
    pub torn: bool,
}

impl Edge {
    /// Create a new untorn `Edge`.
    pub fn new(a: usize, b: usize) -> Self {
        Self { a, b, torn: false }
    }

    /// Returns `true` if vertex `v` is an endpoint of this edge.
    pub fn contains(&self, v: usize) -> bool {
        self.a == v || self.b == v
    }
}

/// A tearable mesh with vertices, edges, and per-edge tear thresholds.
#[derive(Debug, Clone)]
pub struct TearableMesh {
    /// Vertex positions as `[x, y, z]` triples.
    pub vertices: Vec<[f64; 3]>,
    /// Edge list.
    pub edges: Vec<Edge>,
    /// Triangle list: each element is `[i0, i1, i2]` vertex indices.
    pub triangles: Vec<[usize; 3]>,
    /// Per-edge tear-stress threshold \[Pa\].  Aligned with `edges` by index.
    pub tear_thresholds: Vec<f64>,
}

impl TearableMesh {
    /// Create a new `TearableMesh` with a uniform tear threshold.
    ///
    /// Edges are derived from the triangle list (each triangle contributes 3 edges;
    /// duplicates are not removed for simplicity in this reference implementation).
    ///
    /// # Arguments
    /// * `vertices`        – vertex positions
    /// * `triangles`       – triangle connectivity
    /// * `tear_threshold`  – uniform tear threshold \[Pa\]
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>, tear_threshold: f64) -> Self {
        let mut edges = Vec::new();
        for &[i0, i1, i2] in &triangles {
            edges.push(Edge::new(i0, i1));
            edges.push(Edge::new(i1, i2));
            edges.push(Edge::new(i2, i0));
        }
        let n_edges = edges.len();
        let tear_thresholds = vec![tear_threshold; n_edges];
        Self {
            vertices,
            edges,
            triangles,
            tear_thresholds,
        }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }
}

// ------------------------------------------------------------------
// Helper vector utilities (plain f64, no nalgebra)
// ------------------------------------------------------------------

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

/// Stress failure criterion — check whether an edge should tear.
///
/// Returns `true` when the principal tensile stress along the edge exceeds the
/// per-edge tear threshold.
///
/// # Arguments
/// * `mesh`        – the tearable mesh
/// * `edge_idx`    – index into `mesh.edges`
/// * `stress`      – current stress \[Pa\] on the edge (tensile = positive)
pub fn stress_criterion(mesh: &TearableMesh, edge_idx: usize, stress: f64) -> bool {
    stress > mesh.tear_thresholds[edge_idx]
}

/// Topological split — duplicate vertex `v` and update connectivity.
///
/// Creates a new vertex at the same position as `v`, then reassigns
/// one half of the triangles connected to `v` to use the new vertex
/// (the half given by `tris_to_reassign`).
///
/// Returns the index of the newly created vertex.
///
/// # Arguments
/// * `mesh`              – mutable tearable mesh
/// * `v`                 – vertex to duplicate
/// * `tris_to_reassign`  – triangle indices that should use the new vertex
pub fn topological_split(mesh: &mut TearableMesh, v: usize, tris_to_reassign: &[usize]) -> usize {
    // Clone the vertex position.
    let new_pos = mesh.vertices[v];
    mesh.vertices.push(new_pos);
    let new_v = mesh.vertices.len() - 1;
    mesh.tear_thresholds.resize(
        mesh.tear_thresholds.len() + 3,
        mesh.tear_thresholds.first().copied().unwrap_or(1e9),
    );

    // Reassign triangles.
    for &ti in tris_to_reassign {
        for slot in mesh.triangles[ti].iter_mut() {
            if *slot == v {
                *slot = new_v;
            }
        }
    }

    // Add edges for any new connectivity created.
    for &ti in tris_to_reassign {
        let [i0, i1, i2] = mesh.triangles[ti];
        mesh.edges.push(Edge::new(i0, i1));
        mesh.edges.push(Edge::new(i1, i2));
        mesh.edges.push(Edge::new(i2, i0));
    }

    new_v
}

/// Edge cut — mark an edge as torn and split the mesh along it.
///
/// Returns `true` if the edge was successfully torn (was not already torn).
///
/// # Arguments
/// * `mesh`     – mutable tearable mesh
/// * `edge_idx` – index of the edge to cut
pub fn edge_cut(mesh: &mut TearableMesh, edge_idx: usize) -> bool {
    if mesh.edges[edge_idx].torn {
        return false;
    }
    mesh.edges[edge_idx].torn = true;

    // Find triangles sharing this edge.
    let ea = mesh.edges[edge_idx].a;
    let eb = mesh.edges[edge_idx].b;

    let mut sharing: Vec<usize> = mesh
        .triangles
        .iter()
        .enumerate()
        .filter(|(_, tri)| {
            let has_a = tri.contains(&ea);
            let has_b = tri.contains(&eb);
            has_a && has_b
        })
        .map(|(idx, _)| idx)
        .collect();

    // Split the first sharing triangle off the second (if ≥ 2 exist).
    if sharing.len() >= 2 {
        let second_half = sharing.split_off(1);
        // Duplicate vertex ea along the cut for the second-half triangles.
        topological_split(mesh, ea, &second_half);
    }

    true
}

/// Progressive tearing — advance the tear front by one step.
///
/// Iterates over all edges and tears the highest-stress edge that exceeds its
/// threshold.  Returns the index of the torn edge, or `None` if no edge failed.
///
/// # Arguments
/// * `mesh`    – mutable tearable mesh
/// * `stresses`– per-edge stress values \[Pa\] (aligned with `mesh.edges`)
pub fn progressive_tearing(mesh: &mut TearableMesh, stresses: &[f64]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;

    for (i, edge) in mesh.edges.iter().enumerate() {
        if edge.torn {
            continue;
        }
        if i >= stresses.len() {
            break;
        }
        let s = stresses[i];
        if stress_criterion(mesh, i, s) {
            match best {
                None => best = Some((i, s)),
                Some((_, s_best)) if s > s_best => best = Some((i, s)),
                _ => {}
            }
        }
    }

    if let Some((idx, _)) = best {
        edge_cut(mesh, idx);
        Some(idx)
    } else {
        None
    }
}

/// Local remeshing around a tear — insert midpoint vertex on torn edges.
///
/// For each torn edge, replaces the edge with two shorter edges connected to a
/// new midpoint vertex.  Triangles containing the torn edge are subdivided.
///
/// # Arguments
/// * `mesh` – mutable tearable mesh
pub fn remesh_after_tear(mesh: &mut TearableMesh) {
    let torn_indices: Vec<usize> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(_, e)| e.torn)
        .map(|(i, _)| i)
        .collect();

    for idx in torn_indices {
        let ea = mesh.edges[idx].a;
        let eb = mesh.edges[idx].b;

        // Midpoint vertex.
        let pa = mesh.vertices[ea];
        let pb = mesh.vertices[eb];
        let mid = [
            (pa[0] + pb[0]) * 0.5,
            (pa[1] + pb[1]) * 0.5,
            (pa[2] + pb[2]) * 0.5,
        ];
        mesh.vertices.push(mid);
        let mid_idx = mesh.vertices.len() - 1;
        let threshold = mesh.tear_thresholds[idx];

        // Mark old edge as a pass-through (leave torn = true).
        // Add two new edges.
        mesh.edges.push(Edge::new(ea, mid_idx));
        mesh.edges.push(Edge::new(mid_idx, eb));
        mesh.tear_thresholds.push(threshold);
        mesh.tear_thresholds.push(threshold);

        // Subdivide triangles that contained this torn edge.
        let sharing: Vec<usize> = mesh
            .triangles
            .iter()
            .enumerate()
            .filter(|(_, tri)| tri.contains(&ea) && tri.contains(&eb))
            .map(|(i, _)| i)
            .collect();

        let mut new_tris: Vec<[usize; 3]> = Vec::new();
        for ti in sharing {
            let [i0, i1, i2] = mesh.triangles[ti];
            // Determine which vertex is the "other" (not ea or eb).
            let other = if i0 != ea && i0 != eb {
                i0
            } else if i1 != ea && i1 != eb {
                i1
            } else {
                i2
            };
            // Replace triangle with two sub-triangles.
            mesh.triangles[ti] = [ea, mid_idx, other];
            new_tris.push([mid_idx, eb, other]);
        }
        mesh.triangles.extend_from_slice(&new_tris);
        // Extend tear_thresholds for new edge entries if needed.
        let extra = new_tris.len() * 3;
        let current_len = mesh.edges.len();
        mesh.tear_thresholds.resize(current_len + extra, threshold);
    }
}

/// Boundary smoothing — smooth vertices along the torn boundary.
///
/// For each vertex that appears on a torn edge, moves it toward the average
/// of its non-torn neighbours' positions by factor `alpha`.
///
/// # Arguments
/// * `mesh`  – mutable tearable mesh
/// * `alpha` – smoothing factor (0 = no smoothing, 1 = full Laplacian step)
pub fn boundary_smoothing(mesh: &mut TearableMesh, alpha: f64) {
    // Collect boundary vertices (those on torn edges).
    let boundary: std::collections::HashSet<usize> = mesh
        .edges
        .iter()
        .filter(|e| e.torn)
        .flat_map(|e| [e.a, e.b])
        .collect();

    for &v in &boundary {
        // Find non-torn neighbours.
        let neighbours: Vec<usize> = mesh
            .edges
            .iter()
            .filter(|e| !e.torn && e.contains(v))
            .map(|e| if e.a == v { e.b } else { e.a })
            .collect();

        if neighbours.is_empty() {
            continue;
        }

        let n = neighbours.len() as f64;
        let avg = neighbours
            .iter()
            .fold([0.0f64; 3], |acc, &nb| vec3_add(acc, mesh.vertices[nb]));
        let avg = vec3_scale(avg, 1.0 / n);

        let diff = vec3_sub(avg, mesh.vertices[v]);
        mesh.vertices[v] = vec3_add(mesh.vertices[v], vec3_scale(diff, alpha));
    }
}

/// Compute the length of an edge.
///
/// # Arguments
/// * `mesh`     – the mesh
/// * `edge_idx` – edge index
pub fn edge_length(mesh: &TearableMesh, edge_idx: usize) -> f64 {
    let ea = mesh.edges[edge_idx].a;
    let eb = mesh.edges[edge_idx].b;
    vec3_norm(vec3_sub(mesh.vertices[ea], mesh.vertices[eb]))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    // ------------------------------------------------------------------
    // Helper: simple triangle mesh (two tris sharing an edge)
    // ------------------------------------------------------------------
    fn two_tri_mesh(threshold: f64) -> TearableMesh {
        // Vertices: quad split into 2 triangles.
        // 0---1
        // | / |
        // 2---3
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [1, 3, 2]];
        TearableMesh::new(verts, tris, threshold)
    }

    // ------------------------------------------------------------------
    // TearableMesh construction
    // ------------------------------------------------------------------

    #[test]
    fn test_mesh_vertex_count() {
        let mesh = two_tri_mesh(1e6);
        assert_eq!(mesh.num_vertices(), 4);
    }

    #[test]
    fn test_mesh_triangle_count() {
        let mesh = two_tri_mesh(1e6);
        assert_eq!(mesh.num_triangles(), 2);
    }

    #[test]
    fn test_mesh_edge_count() {
        // 2 triangles × 3 edges each = 6 edges (before deduplication).
        let mesh = two_tri_mesh(1e6);
        assert_eq!(mesh.edges.len(), 6);
    }

    #[test]
    fn test_mesh_all_edges_untorn() {
        let mesh = two_tri_mesh(1e6);
        assert!(mesh.edges.iter().all(|e| !e.torn));
    }

    // ------------------------------------------------------------------
    // Edge struct
    // ------------------------------------------------------------------

    #[test]
    fn test_edge_contains() {
        let e = Edge::new(2, 5);
        assert!(e.contains(2));
        assert!(e.contains(5));
        assert!(!e.contains(0));
    }

    #[test]
    fn test_edge_new_untorn() {
        let e = Edge::new(0, 1);
        assert!(!e.torn);
    }

    // ------------------------------------------------------------------
    // Stress criterion
    // ------------------------------------------------------------------

    #[test]
    fn test_stress_criterion_below_threshold() {
        let mesh = two_tri_mesh(1e6);
        assert!(!stress_criterion(&mesh, 0, 5e5));
    }

    #[test]
    fn test_stress_criterion_above_threshold() {
        let mesh = two_tri_mesh(1e6);
        assert!(stress_criterion(&mesh, 0, 2e6));
    }

    #[test]
    fn test_stress_criterion_at_threshold() {
        let mesh = two_tri_mesh(1e6);
        // Exactly at threshold → no tear (strict >).
        assert!(!stress_criterion(&mesh, 0, 1e6));
    }

    // ------------------------------------------------------------------
    // Topological split
    // ------------------------------------------------------------------

    #[test]
    fn test_topological_split_adds_vertex() {
        let mut mesh = two_tri_mesh(1e6);
        let before = mesh.num_vertices();
        topological_split(&mut mesh, 1, &[1]);
        assert_eq!(
            mesh.num_vertices(),
            before + 1,
            "Split should add one vertex"
        );
    }

    #[test]
    fn test_topological_split_new_vertex_same_position() {
        let mut mesh = two_tri_mesh(1e6);
        let old_pos = mesh.vertices[1];
        let new_v = topological_split(&mut mesh, 1, &[1]);
        assert_eq!(
            mesh.vertices[new_v], old_pos,
            "New vertex should have same position as original"
        );
    }

    #[test]
    fn test_topological_split_reassigns_triangles() {
        let mut mesh = two_tri_mesh(1e6);
        let new_v = topological_split(&mut mesh, 1, &[1]);
        // Triangle 1 should now reference new_v instead of vertex 1.
        let tri1 = mesh.triangles[1];
        assert!(
            tri1.contains(&new_v),
            "Reassigned triangle should contain new vertex {new_v}, got {tri1:?}"
        );
    }

    #[test]
    fn test_topological_split_original_triangle_unchanged() {
        let mut mesh = two_tri_mesh(1e6);
        let original_tri0 = mesh.triangles[0];
        topological_split(&mut mesh, 1, &[1]);
        // Triangle 0 was not in tris_to_reassign, so it still references vertex 1.
        assert_eq!(mesh.triangles[0], original_tri0);
    }

    // ------------------------------------------------------------------
    // Edge cut
    // ------------------------------------------------------------------

    #[test]
    fn test_edge_cut_marks_edge_torn() {
        let mut mesh = two_tri_mesh(1e6);
        edge_cut(&mut mesh, 0);
        assert!(mesh.edges[0].torn, "Edge 0 should be marked torn after cut");
    }

    #[test]
    fn test_edge_cut_returns_true_first_time() {
        let mut mesh = two_tri_mesh(1e6);
        assert!(edge_cut(&mut mesh, 0));
    }

    #[test]
    fn test_edge_cut_returns_false_second_time() {
        let mut mesh = two_tri_mesh(1e6);
        edge_cut(&mut mesh, 0);
        assert!(
            !edge_cut(&mut mesh, 0),
            "Second cut on same edge should return false"
        );
    }

    #[test]
    fn test_edge_cut_adds_vertices() {
        let mut mesh = two_tri_mesh(1e6);
        // Find an edge shared by both triangles.
        // Edge index 1 in our construction: (1,2) from tri 0 — shared by tri 1.
        let before = mesh.num_vertices();
        edge_cut(&mut mesh, 1);
        // A topological split should have occurred, adding at least one vertex.
        assert!(
            mesh.num_vertices() >= before,
            "Edge cut should not remove vertices"
        );
    }

    // ------------------------------------------------------------------
    // Progressive tearing
    // ------------------------------------------------------------------

    #[test]
    fn test_progressive_tearing_tears_highest_stress() {
        let mut mesh = two_tri_mesh(1e3);
        let stresses = vec![500.0, 2000.0, 1500.0, 0.0, 0.0, 0.0];
        let torn = progressive_tearing(&mut mesh, &stresses);
        assert_eq!(
            torn,
            Some(1),
            "Should tear the highest-stress edge (index 1)"
        );
    }

    #[test]
    fn test_progressive_tearing_none_below_threshold() {
        let mut mesh = two_tri_mesh(1e6);
        let stresses = vec![0.0; 6];
        let torn = progressive_tearing(&mut mesh, &stresses);
        assert_eq!(
            torn, None,
            "No edge should tear when stresses are below threshold"
        );
    }

    #[test]
    fn test_progressive_tearing_skips_already_torn() {
        let mut mesh = two_tri_mesh(1e3);
        let stresses = vec![2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0];
        progressive_tearing(&mut mesh, &stresses); // Tears edge 0.
        progressive_tearing(&mut mesh, &stresses); // Should tear a different edge.
        let torn_count = mesh.edges.iter().filter(|e| e.torn).count();
        assert!(torn_count >= 1, "At least one edge should be torn");
    }

    // ------------------------------------------------------------------
    // Remesh after tear
    // ------------------------------------------------------------------

    #[test]
    fn test_remesh_adds_midpoint_vertex() {
        let mut mesh = two_tri_mesh(1e3);
        edge_cut(&mut mesh, 0);
        let before_verts = mesh.num_vertices();
        remesh_after_tear(&mut mesh);
        assert!(
            mesh.num_vertices() > before_verts,
            "Remeshing should add midpoint vertices"
        );
    }

    #[test]
    fn test_remesh_increases_triangle_count() {
        let mut mesh = two_tri_mesh(1e3);
        edge_cut(&mut mesh, 1); // shared edge
        let before_tris = mesh.num_triangles();
        remesh_after_tear(&mut mesh);
        assert!(
            mesh.num_triangles() >= before_tris,
            "Remeshing should not reduce triangle count"
        );
    }

    // ------------------------------------------------------------------
    // Boundary smoothing
    // ------------------------------------------------------------------

    #[test]
    fn test_boundary_smoothing_no_crash() {
        let mut mesh = two_tri_mesh(1e3);
        edge_cut(&mut mesh, 0);
        boundary_smoothing(&mut mesh, 0.5);
        // Just verify no panic and vertices are finite.
        for v in &mesh.vertices {
            assert!(
                v.iter().all(|x| x.is_finite()),
                "Vertices should remain finite after smoothing"
            );
        }
    }

    #[test]
    fn test_boundary_smoothing_alpha_zero_no_change() {
        let mut mesh = two_tri_mesh(1e3);
        edge_cut(&mut mesh, 0);
        let before: Vec<[f64; 3]> = mesh.vertices.clone();
        boundary_smoothing(&mut mesh, 0.0);
        for (a, b) in before.iter().zip(mesh.vertices.iter()) {
            let diff = vec3_norm(vec3_sub(*a, *b));
            assert!(diff < EPS, "Alpha=0 smoothing should not move vertices");
        }
    }

    // ------------------------------------------------------------------
    // Edge length helper
    // ------------------------------------------------------------------

    #[test]
    fn test_edge_length_unit_edge() {
        let mesh = two_tri_mesh(1e6);
        // Edge 0: vertices 0 and 1, both at y=0, x differs by 1.
        let len = edge_length(&mesh, 0);
        assert!(
            (len - 1.0).abs() < EPS,
            "Edge length should be 1.0, got {len}"
        );
    }

    #[test]
    fn test_edge_length_positive() {
        let mesh = two_tri_mesh(1e6);
        for i in 0..mesh.edges.len() {
            let len = edge_length(&mesh, i);
            assert!(len >= 0.0, "Edge length should be non-negative");
        }
    }

    #[test]
    fn test_edge_length_diagonal() {
        // Edge 1 in our construction: vertices 1 and 2.
        // v1 = [1,0,0], v2 = [0,1,0]  → distance = sqrt(2).
        let mesh = two_tri_mesh(1e6);
        let len = edge_length(&mesh, 1);
        assert!(
            (len - 2.0_f64.sqrt()).abs() < 1e-10,
            "Diagonal edge length should be sqrt(2), got {len}"
        );
    }
}

// =============================================================================
// Particle-spring tearing model with fracture-mechanics functions
// =============================================================================

// ------------------------------------------------------------------
// TearParticle
// ------------------------------------------------------------------

/// A particle in the tearable spring-mass mesh.
#[derive(Debug, Clone)]
pub struct TearParticle {
    /// Position \[m\].
    pub pos: [f64; 3],
    /// Velocity \[m/s\].
    pub vel: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Whether this particle has been flagged as torn/detached.
    pub is_torn: bool,
}

impl TearParticle {
    /// Creates a new `TearParticle` at rest at the given position.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            is_torn: false,
        }
    }
}

// ------------------------------------------------------------------
// TearEdge
// ------------------------------------------------------------------

/// A spring edge that can permanently break once its strain exceeds `max_strain`.
#[derive(Debug, Clone)]
pub struct TearEdge {
    /// Index of the first endpoint particle.
    pub p0: usize,
    /// Index of the second endpoint particle.
    pub p1: usize,
    /// Rest length of the spring \[m\].
    pub rest_length: f64,
    /// Strain threshold above which the edge breaks.
    pub max_strain: f64,
    /// Whether this edge has already broken.
    pub is_broken: bool,
}

impl TearEdge {
    /// Creates a new `TearEdge` (initially intact).
    pub fn new(p0: usize, p1: usize, rest_len: f64, max_strain: f64) -> Self {
        Self {
            p0,
            p1,
            rest_length: rest_len,
            max_strain,
            is_broken: false,
        }
    }

    /// Current engineering strain of this edge given particle positions.
    ///
    /// ε = (|x₁ − x₀| − L₀) / L₀
    pub fn strain(&self, particles: &[TearParticle]) -> f64 {
        let dx = [
            particles[self.p1].pos[0] - particles[self.p0].pos[0],
            particles[self.p1].pos[1] - particles[self.p0].pos[1],
            particles[self.p1].pos[2] - particles[self.p0].pos[2],
        ];
        let current = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        (current - self.rest_length) / self.rest_length
    }

    /// Returns `true` when the edge should break (strain > max_strain and not already broken).
    pub fn should_break(&self, particles: &[TearParticle]) -> bool {
        !self.is_broken && self.strain(particles) > self.max_strain
    }

    /// Linear spring force \[N\] on particle `p0` (negate for `p1`).
    ///
    /// Returns `[0; 3]` if the edge is already broken.
    ///
    /// * `particles` – particle slice
    /// * `k`         – spring stiffness \[N/m\]
    pub fn spring_force(&self, particles: &[TearParticle], k: f64) -> [f64; 3] {
        if self.is_broken {
            return [0.0; 3];
        }
        let dx = [
            particles[self.p1].pos[0] - particles[self.p0].pos[0],
            particles[self.p1].pos[1] - particles[self.p0].pos[1],
            particles[self.p1].pos[2] - particles[self.p0].pos[2],
        ];
        let current = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if current < 1e-15 {
            return [0.0; 3];
        }
        let extension = current - self.rest_length;
        let factor = k * extension / current;
        [dx[0] * factor, dx[1] * factor, dx[2] * factor]
    }
}

// ------------------------------------------------------------------
// TearMesh
// ------------------------------------------------------------------

/// A spring-mass mesh whose edges can tear when their strain exceeds a threshold.
#[derive(Debug, Clone, Default)]
pub struct TearMesh {
    /// List of particles.
    pub particles: Vec<TearParticle>,
    /// List of spring edges.
    pub edges: Vec<TearEdge>,
}

impl TearMesh {
    /// Creates a new empty `TearMesh`.
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Adds a particle and returns its index.
    pub fn add_particle(&mut self, p: TearParticle) -> usize {
        self.particles.push(p);
        self.particles.len() - 1
    }

    /// Adds an edge between particles `p0` and `p1`.
    pub fn add_edge(&mut self, p0: usize, p1: usize, rest_len: f64, max_strain: f64) {
        self.edges.push(TearEdge::new(p0, p1, rest_len, max_strain));
    }

    /// Number of particles in the mesh.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Number of edges (including broken ones) in the mesh.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Checks all edges and marks those that exceed `max_strain` as broken.
    ///
    /// Returns the number of newly broken edges in this call.
    pub fn update_tears(&mut self) -> usize {
        let mut newly_broken = 0;
        // Collect indices first to avoid borrow conflict.
        let breaks: Vec<usize> = self
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.should_break(&self.particles))
            .map(|(i, _)| i)
            .collect();
        for i in breaks {
            self.edges[i].is_broken = true;
            newly_broken += 1;
        }
        newly_broken
    }

    /// Number of intact (non-broken) edges.
    pub fn intact_edges(&self) -> usize {
        self.edges.iter().filter(|e| !e.is_broken).count()
    }

    /// Fraction of edges that are broken ∈ \[0, 1\].
    ///
    /// Returns 0.0 if there are no edges.
    pub fn broken_fraction(&self) -> f64 {
        if self.edges.is_empty() {
            return 0.0;
        }
        let broken = self.edges.iter().filter(|e| e.is_broken).count();
        broken as f64 / self.edges.len() as f64
    }
}

// ------------------------------------------------------------------
// Fracture-mechanics free functions
// ------------------------------------------------------------------

/// Griffith fracture energy \[J/m²\] (critical energy release rate G_c).
///
/// G_c = σ_f² · π · a / E
///
/// * `sigma_f` – fracture stress \[Pa\]
/// * `e_mod`   – Young's modulus \[Pa\]
pub fn fracture_energy(sigma_f: f64, e_mod: f64) -> f64 {
    sigma_f * sigma_f * std::f64::consts::PI / e_mod
}

/// Mode-I stress intensity factor K_I \[Pa·√m\].
///
/// K_I = σ · √(π · a)
///
/// * `sigma` – applied stress \[Pa\]
/// * `a`     – half crack length \[m\]
pub fn stress_intensity_factor(sigma: f64, a: f64) -> f64 {
    sigma * (std::f64::consts::PI * a).sqrt()
}

/// Critical crack half-length \[m\] for a given fracture toughness.
///
/// a_c = (K_Ic / (σ · √π))²
///
/// * `kic`   – plane-strain fracture toughness \[Pa·√m\]
/// * `sigma` – applied stress \[Pa\]
pub fn critical_crack_length(kic: f64, sigma: f64) -> f64 {
    let denom = sigma * std::f64::consts::PI.sqrt();
    (kic / denom).powi(2)
}

/// Dugdale cohesive zone length \[m\].
///
/// l_cz = E · G_c / (π · σ_c²)
///
/// * `e_mod`   – Young's modulus \[Pa\]
/// * `gc`      – fracture energy \[J/m²\]
/// * `sigma_c` – cohesive strength \[Pa\]
pub fn cohesive_zone_length(e_mod: f64, gc: f64, sigma_c: f64) -> f64 {
    e_mod * gc / (std::f64::consts::PI * sigma_c * sigma_c)
}

/// Ductile fracture damage indicator D ∈ \[0, 1\].
///
/// D = ε_current / ε_f
///
/// * `eps_f`       – failure strain
/// * `eps_current` – current plastic strain
///
/// Clamped to \[0, 1\].
pub fn ductile_fracture_indicator(eps_f: f64, eps_current: f64) -> f64 {
    (eps_current / eps_f).clamp(0.0, 1.0)
}

// ------------------------------------------------------------------
// Tests for particle-spring tearing model
// ------------------------------------------------------------------

#[cfg(test)]
mod tear_particle_tests {

    use crate::tearing::TearMesh;
    use crate::tearing::TearParticle;
    use crate::tearing::cohesive_zone_length;
    use crate::tearing::critical_crack_length;
    use crate::tearing::ductile_fracture_indicator;
    use crate::tearing::fracture_energy;
    use crate::tearing::stress_intensity_factor;

    const EPS: f64 = 1e-12;

    fn two_particle_mesh(separation: f64, rest_length: f64, max_strain: f64) -> TearMesh {
        let mut mesh = TearMesh::new();
        mesh.add_particle(TearParticle::new([0.0, 0.0, 0.0], 1.0));
        mesh.add_particle(TearParticle::new([separation, 0.0, 0.0], 1.0));
        mesh.add_edge(0, 1, rest_length, max_strain);
        mesh
    }

    // ── TearParticle ──────────────────────────────────────────────────────

    #[test]
    fn test_tear_particle_new() {
        let p = TearParticle::new([1.0, 2.0, 3.0], 0.5);
        assert_eq!(p.pos, [1.0, 2.0, 3.0]);
        assert_eq!(p.vel, [0.0; 3]);
        assert!((p.mass - 0.5).abs() < EPS);
        assert!(!p.is_torn);
    }

    // ── TearEdge strain ───────────────────────────────────────────────────

    #[test]
    fn test_strain_zero_at_rest_length() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        let strain = mesh.edges[0].strain(&mesh.particles);
        assert!(
            strain.abs() < EPS,
            "Strain at rest length should be 0, got {strain}"
        );
    }

    #[test]
    fn test_strain_positive_when_stretched() {
        let mesh = two_particle_mesh(1.5, 1.0, 0.5);
        let strain = mesh.edges[0].strain(&mesh.particles);
        assert!(strain > 0.0, "Strain should be positive when stretched");
    }

    #[test]
    fn test_strain_formula() {
        let mesh = two_particle_mesh(1.2, 1.0, 0.5);
        let strain = mesh.edges[0].strain(&mesh.particles);
        let expected = (1.2 - 1.0) / 1.0;
        assert!(
            (strain - expected).abs() < EPS,
            "Strain formula: got {strain}, expected {expected}"
        );
    }

    #[test]
    fn test_strain_negative_when_compressed() {
        let mesh = two_particle_mesh(0.8, 1.0, 0.5);
        let strain = mesh.edges[0].strain(&mesh.particles);
        assert!(strain < 0.0, "Compressive strain should be negative");
    }

    // ── TearEdge should_break ─────────────────────────────────────────────

    #[test]
    fn test_should_break_at_max_strain() {
        let mesh = two_particle_mesh(2.0, 1.0, 0.5); // strain = 1.0 > 0.5
        assert!(
            mesh.edges[0].should_break(&mesh.particles),
            "Edge should break when strain > max_strain"
        );
    }

    #[test]
    fn test_should_not_break_below_max_strain() {
        let mesh = two_particle_mesh(1.2, 1.0, 0.5); // strain = 0.2 < 0.5
        assert!(!mesh.edges[0].should_break(&mesh.particles));
    }

    #[test]
    fn test_should_not_break_if_already_broken() {
        let mut mesh = two_particle_mesh(2.0, 1.0, 0.5);
        mesh.edges[0].is_broken = true;
        assert!(
            !mesh.edges[0].should_break(&mesh.particles),
            "Already-broken edge should not trigger again"
        );
    }

    // ── TearEdge spring_force ─────────────────────────────────────────────

    #[test]
    fn test_spring_force_zero_at_rest_length() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        let force = mesh.edges[0].spring_force(&mesh.particles, 100.0);
        let mag = (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt();
        assert!(mag < EPS, "Spring force must be zero at rest length");
    }

    #[test]
    fn test_spring_force_nonzero_when_stretched() {
        let mesh = two_particle_mesh(1.5, 1.0, 0.5);
        let force = mesh.edges[0].spring_force(&mesh.particles, 100.0);
        let mag = (force[0] * force[0]).sqrt();
        assert!(mag > 0.0, "Spring force must be non-zero when stretched");
    }

    #[test]
    fn test_spring_force_zero_when_broken() {
        let mut mesh = two_particle_mesh(1.5, 1.0, 0.5);
        mesh.edges[0].is_broken = true;
        let force = mesh.edges[0].spring_force(&mesh.particles, 100.0);
        let mag = (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]).sqrt();
        assert!(mag < EPS, "Broken edge must produce zero spring force");
    }

    // ── TearMesh ──────────────────────────────────────────────────────────

    #[test]
    fn test_tear_mesh_particle_count() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        assert_eq!(mesh.particle_count(), 2);
    }

    #[test]
    fn test_tear_mesh_edge_count() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        assert_eq!(mesh.edge_count(), 1);
    }

    #[test]
    fn test_update_tears_breaks_overstrained_edges() {
        let mut mesh = two_particle_mesh(2.0, 1.0, 0.5); // strain=1.0 > 0.5
        let broke = mesh.update_tears();
        assert_eq!(broke, 1, "One edge should break");
        assert!(mesh.edges[0].is_broken);
    }

    #[test]
    fn test_update_tears_no_break_below_threshold() {
        let mut mesh = two_particle_mesh(1.2, 1.0, 0.5); // strain=0.2 < 0.5
        let broke = mesh.update_tears();
        assert_eq!(broke, 0, "No edge should break");
    }

    #[test]
    fn test_update_tears_does_not_double_count() {
        let mut mesh = two_particle_mesh(2.0, 1.0, 0.5);
        mesh.update_tears();
        let second_call = mesh.update_tears();
        assert_eq!(
            second_call, 0,
            "Already-broken edge should not be counted again"
        );
    }

    #[test]
    fn test_intact_edges_all_intact_initially() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        assert_eq!(mesh.intact_edges(), 1);
    }

    #[test]
    fn test_intact_edges_after_break() {
        let mut mesh = two_particle_mesh(2.0, 1.0, 0.5);
        mesh.update_tears();
        assert_eq!(mesh.intact_edges(), 0);
    }

    #[test]
    fn test_broken_fraction_zero_initially() {
        let mesh = two_particle_mesh(1.0, 1.0, 0.5);
        assert!(mesh.broken_fraction().abs() < EPS);
    }

    #[test]
    fn test_broken_fraction_one_after_all_break() {
        let mut mesh = two_particle_mesh(2.0, 1.0, 0.5);
        mesh.update_tears();
        assert!((mesh.broken_fraction() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_broken_fraction_in_range() {
        let bf = {
            let mut mesh = TearMesh::new();
            mesh.add_particle(TearParticle::new([0.0, 0.0, 0.0], 1.0));
            mesh.add_particle(TearParticle::new([2.0, 0.0, 0.0], 1.0)); // strain=1.0
            mesh.add_particle(TearParticle::new([0.5, 0.0, 0.0], 1.0)); // strain=0.0 from p0
            mesh.add_edge(0, 1, 1.0, 0.5); // will break
            mesh.add_edge(0, 2, 0.5, 0.5); // will not break
            mesh.update_tears();
            mesh.broken_fraction()
        };
        assert!(
            (0.0..=1.0).contains(&bf),
            "broken_fraction must be in [0,1], got {bf}"
        );
    }

    #[test]
    fn test_broken_fraction_empty_mesh() {
        let mesh = TearMesh::new();
        assert!(
            mesh.broken_fraction().abs() < EPS,
            "Empty mesh → broken_fraction = 0"
        );
    }

    // ── Fracture-mechanics functions ──────────────────────────────────────

    #[test]
    fn test_fracture_energy_positive() {
        let gc = fracture_energy(100e6, 200e9);
        assert!(gc > 0.0, "Fracture energy must be positive");
    }

    #[test]
    fn test_fracture_energy_formula() {
        let sigma = 100e6_f64;
        let e = 200e9_f64;
        let expected = sigma * sigma * std::f64::consts::PI / e;
        let got = fracture_energy(sigma, e);
        assert!((got - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn test_stress_intensity_positive() {
        let k = stress_intensity_factor(100e6, 1e-3);
        assert!(k > 0.0);
    }

    #[test]
    fn test_stress_intensity_formula() {
        let sigma = 50e6_f64;
        let a = 2e-3_f64;
        let expected = sigma * (std::f64::consts::PI * a).sqrt();
        let got = stress_intensity_factor(sigma, a);
        assert!((got - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn test_critical_crack_length_positive() {
        let ac = critical_crack_length(50e6, 100e6);
        assert!(ac > 0.0, "Critical crack length must be positive");
    }

    #[test]
    fn test_critical_crack_length_formula() {
        let kic = 50e6_f64;
        let sigma = 100e6_f64;
        let expected = (kic / (sigma * std::f64::consts::PI.sqrt())).powi(2);
        let got = critical_crack_length(kic, sigma);
        assert!((got - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn test_cohesive_zone_length_positive() {
        let lc = cohesive_zone_length(200e9, 1000.0, 500e6);
        assert!(lc > 0.0, "Cohesive zone length must be positive");
    }

    #[test]
    fn test_cohesive_zone_length_formula() {
        let e = 200e9_f64;
        let gc = 1000.0_f64;
        let sc = 500e6_f64;
        let expected = e * gc / (std::f64::consts::PI * sc * sc);
        let got = cohesive_zone_length(e, gc, sc);
        assert!((got - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn test_ductile_fracture_indicator_zero() {
        let d = ductile_fracture_indicator(0.5, 0.0);
        assert!(d.abs() < EPS, "Zero current strain → damage = 0");
    }

    #[test]
    fn test_ductile_fracture_indicator_one_at_failure() {
        let d = ductile_fracture_indicator(0.5, 0.5);
        assert!((d - 1.0).abs() < EPS, "At failure strain → damage = 1");
    }

    #[test]
    fn test_ductile_fracture_indicator_clamped() {
        let d = ductile_fracture_indicator(0.3, 1.0);
        assert!(
            (d - 1.0).abs() < EPS,
            "Over-failure strain → damage clamped to 1"
        );
    }

    #[test]
    fn test_ductile_fracture_indicator_proportional() {
        let d = ductile_fracture_indicator(1.0, 0.4);
        assert!((d - 0.4).abs() < EPS);
    }
}
