//! Tests for [`TopologicalSurgery`]: every operation is checked for the Euler
//! characteristic, closedness, genus and the half-edge structure invariants,
//! and every failure for leaving the mesh untouched.

use super::fixtures::{
    Fixture, connected_sum, disjoint_union, grid, grid_torus, icosahedron, octahedron,
    seven_vertex_torus, tetrahedron, tube, without_faces,
};
use super::functions::extract_boundary_loops;
use super::homology::betti_numbers_by_dimension;
use super::surface::Surface;
use super::surgery::TopologicalSurgery;
use super::types::{HalfEdge, HalfEdgeMesh, Simplex, SimplicialComplex};
use crate::error::Error;

fn build(fixture: &Fixture) -> HalfEdgeMesh {
    HalfEdgeMesh::from_triangles(fixture.0.clone(), &fixture.1)
}

/// Everything stored in the mesh, in a comparable form.
type Snapshot = (Vec<[f64; 3]>, Vec<HalfEdge>, Vec<usize>, Vec<(usize, bool)>);

fn snapshot(mesh: &HalfEdgeMesh) -> Snapshot {
    (
        mesh.vertices.clone(),
        mesh.half_edges.clone(),
        mesh.vertex_he.clone(),
        mesh.faces
            .iter()
            .map(|f| (f.start_he, f.is_outer))
            .collect(),
    )
}

/// `(vertices, unique edges, faces)`.
fn counts(mesh: &HalfEdgeMesh) -> (usize, usize, usize) {
    (
        mesh.num_vertices(),
        mesh.count_unique_edges(),
        mesh.num_faces(),
    )
}

/// Structural invariants of the half-edge arrays.
fn assert_half_edge_invariants(mesh: &HalfEdgeMesh) {
    let h = mesh.half_edges.len();
    assert_eq!(h, 3 * mesh.faces.len(), "three half-edges per triangle");
    assert_eq!(mesh.vertex_he.len(), mesh.vertices.len());
    for (i, he) in mesh.half_edges.iter().enumerate() {
        assert!(
            he.origin < mesh.vertices.len(),
            "half-edge {i}: origin range"
        );
        assert!(he.next < h && he.prev < h, "half-edge {i}: pointer range");
        assert_eq!(mesh.half_edges[he.next].prev, i, "next/prev inverse at {i}");
        assert_eq!(mesh.half_edges[he.prev].next, i, "prev/next inverse at {i}");
        assert_eq!(mesh.half_edges[he.next].face, he.face, "face loop at {i}");
        let after_two = mesh.half_edges[mesh.half_edges[he.next].next].next;
        assert_eq!(after_two, i, "face of half-edge {i} is not a triangle");
        if he.twin != usize::MAX {
            let twin = mesh.half_edges[he.twin];
            assert_eq!(twin.twin, i, "twin symmetry at {i}");
            assert_eq!(
                twin.origin, mesh.half_edges[he.next].origin,
                "twin origin at {i}"
            );
            assert_eq!(
                mesh.half_edges[twin.next].origin, he.origin,
                "twin destination at {i}"
            );
        }
    }
    for (f, face) in mesh.faces.iter().enumerate() {
        assert_eq!(mesh.half_edges[face.start_he].face, f, "face {f} start");
    }
    for (v, &he) in mesh.vertex_he.iter().enumerate() {
        assert_ne!(he, usize::MAX, "vertex {v} is isolated");
        assert_eq!(mesh.half_edges[he].origin, v, "vertex_he of {v}");
    }
    Surface::from_mesh(mesh).expect("surgery result must be a valid oriented 2-manifold");
}

fn half_edge_between(mesh: &HalfEdgeMesh, u: usize, v: usize) -> usize {
    (0..mesh.half_edges.len())
        .find(|&i| {
            mesh.half_edges[i].origin == u && mesh.half_edges[mesh.half_edges[i].next].origin == v
        })
        .unwrap_or_else(|| panic!("no half-edge {u} -> {v}"))
}

fn components(mesh: &HalfEdgeMesh) -> usize {
    let (surface, _) = Surface::from_mesh(mesh).expect("valid mesh");
    surface.component_count()
}

fn homology_of(mesh: &HalfEdgeMesh) -> Vec<usize> {
    let mut complex = SimplicialComplex::new(2);
    for face in 0..mesh.num_faces() {
        let vertices: Vec<usize> = mesh
            .face_half_edges(face)
            .iter()
            .map(|&he| mesh.half_edges[he].origin)
            .collect();
        complex.add_simplex(Simplex::new(vertices));
    }
    betti_numbers_by_dimension(&complex)
}

struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) as usize) % n
    }
}

// ---------------------------------------------------------------------------
// edge_split
// ---------------------------------------------------------------------------

#[test]
fn edge_split_every_half_edge_of_a_tetrahedron() {
    let fixture = tetrahedron();
    for he in 0..12 {
        let mut mesh = build(&fixture);
        let origin = mesh.half_edges[he].origin;
        let destination = mesh.half_edges[mesh.half_edges[he].next].origin;
        let new_vertex = TopologicalSurgery::new(&mut mesh)
            .edge_split(he)
            .expect("splitting a tetrahedron edge");
        assert_eq!(new_vertex, 4);
        assert_eq!(counts(&mesh), (5, 9, 6), "half-edge {he}");
        assert_eq!(mesh.euler_characteristic(), 2);
        assert!(mesh.is_closed());
        assert_eq!(mesh.genus(), 0);
        assert_eq!(&mesh.vertices[..4], &fixture.0[..]);
        for k in 0..3 {
            let midpoint = (fixture.0[origin][k] + fixture.0[destination][k]) / 2.0;
            assert!((mesh.vertices[new_vertex][k] - midpoint).abs() < 1e-15);
        }
        assert_half_edge_invariants(&mesh);
    }
}

#[test]
fn edge_split_boundary_edge_splits_one_triangle() {
    let mut mesh = build(&(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
        vec![[0, 1, 2]],
    ));
    let new_vertex = TopologicalSurgery::new(&mut mesh)
        .edge_split(0)
        .expect("boundary split");
    assert_eq!(new_vertex, 3);
    assert_eq!(counts(&mesh), (4, 5, 2));
    assert_eq!(mesh.euler_characteristic(), 1);
    assert!(!mesh.is_closed());
    assert_eq!(extract_boundary_loops(&mesh).len(), 1);
    assert_eq!(extract_boundary_loops(&mesh)[0].len(), 4);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn edge_split_interior_edge_of_an_open_mesh() {
    let fixture = grid(3, 3, false);
    let mut mesh = build(&fixture);
    let interior = half_edge_between(&mesh, 0, 4);
    let before = counts(&mesh);
    TopologicalSurgery::new(&mut mesh)
        .edge_split(interior)
        .expect("interior split");
    assert_eq!(
        counts(&mesh),
        (before.0 + 1, before.1 + 3, before.2 + 2),
        "an interior edge split adds one vertex, three edges and two faces"
    );
    assert_eq!(mesh.euler_characteristic(), 1);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn edge_split_rejects_a_bad_half_edge_and_keeps_the_mesh() {
    let mut mesh = build(&tetrahedron());
    let before = snapshot(&mesh);
    let err = TopologicalSurgery::new(&mut mesh)
        .edge_split(12)
        .expect_err("half-edge 12 does not exist");
    assert!(matches!(
        err,
        Error::IndexOutOfBounds { index: 12, len: 12 }
    ));
    assert_eq!(snapshot(&mesh), before);
}

// ---------------------------------------------------------------------------
// edge_collapse
// ---------------------------------------------------------------------------

#[test]
fn edge_collapse_every_octahedron_edge() {
    let fixture = octahedron();
    for he in 0..24 {
        let mut mesh = build(&fixture);
        let origin = mesh.half_edges[he].origin;
        let destination = mesh.half_edges[mesh.half_edges[he].next].origin;
        let survivor = TopologicalSurgery::new(&mut mesh)
            .edge_collapse(he)
            .unwrap_or_else(|e| panic!("half-edge {he}: {e}"));
        assert_eq!(counts(&mesh), (5, 9, 6), "half-edge {he}");
        assert_eq!(mesh.euler_characteristic(), 2);
        assert!(mesh.is_closed());
        assert_eq!(mesh.genus(), 0);
        let expected_survivor = if destination > origin {
            destination - 1
        } else {
            destination
        };
        assert_eq!(survivor, expected_survivor);
        for k in 0..3 {
            let midpoint = (fixture.0[origin][k] + fixture.0[destination][k]) / 2.0;
            assert!((mesh.vertices[survivor][k] - midpoint).abs() < 1e-15);
        }
        assert_half_edge_invariants(&mesh);
    }
}

#[test]
fn edge_collapse_rejects_every_tetrahedron_edge() {
    for he in 0..12 {
        let mut mesh = build(&tetrahedron());
        let before = snapshot(&mesh);
        let err = TopologicalSurgery::new(&mut mesh)
            .edge_collapse(he)
            .expect_err("a tetrahedron cannot be collapsed");
        assert!(
            matches!(err, Error::Unsupported { .. }),
            "half-edge {he}: {err}"
        );
        assert_eq!(snapshot(&mesh), before, "half-edge {he} changed the mesh");
    }
}

#[test]
fn edge_collapse_link_condition_on_a_bipyramid() {
    // Collapsing one octahedron edge leaves a triangular bipyramid: two apexes
    // of valence 3 over a triangle of valence-4 equator vertices.
    let mut mesh = build(&octahedron());
    TopologicalSurgery::new(&mut mesh)
        .edge_collapse(0)
        .expect("octahedron edge");
    let snapshot_before = snapshot(&mesh);
    let (surface, _) = Surface::from_mesh(&mesh).expect("valid");
    let valence: Vec<usize> = (0..mesh.num_vertices())
        .map(|v| surface.neighbours(v).len())
        .collect();
    let apex = (0..5).find(|&v| valence[v] == 3).expect("an apex");
    let equator = (0..5)
        .find(|&v| valence[v] == 4)
        .expect("an equator vertex");
    let equator_edge = {
        let other = (0..5)
            .find(|&v| v != equator && valence[v] == 4 && surface.neighbours(equator).contains(&v))
            .expect("another equator vertex");
        half_edge_between(&mesh, equator, other)
    };
    // The equator edge shares a third equator neighbour: link condition fails.
    let err = TopologicalSurgery::new(&mut mesh)
        .edge_collapse(equator_edge)
        .expect_err("link condition");
    assert!(
        matches!(err, Error::Unsupported { .. }),
        "unexpected: {err}"
    );
    assert_eq!(snapshot(&mesh), snapshot_before);
    // Collapsing an apex onto an equator vertex is fine and yields a tetrahedron.
    let apex_edge = half_edge_between(&mesh, apex, equator);
    TopologicalSurgery::new(&mut mesh)
        .edge_collapse(apex_edge)
        .expect("apex collapse");
    assert_eq!(counts(&mesh), (4, 6, 4));
    assert_half_edge_invariants(&mesh);
}

#[test]
fn edge_collapse_boundary_edge_of_an_open_mesh() {
    let mut mesh = build(&grid(3, 3, false));
    let boundary = half_edge_between(&mesh, 0, 1);
    assert_eq!(mesh.half_edges[boundary].twin, usize::MAX);
    let before = counts(&mesh);
    let survivor = TopologicalSurgery::new(&mut mesh)
        .edge_collapse(boundary)
        .expect("boundary edge collapse");
    assert_eq!(survivor, 0);
    assert_eq!(counts(&mesh), (before.0 - 1, before.1 - 2, before.2 - 1));
    assert_eq!(mesh.euler_characteristic(), 1);
    assert_eq!(extract_boundary_loops(&mesh).len(), 1);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn edge_collapse_rejects_an_interior_edge_between_boundary_vertices() {
    let mut mesh = build(&grid(3, 3, false));
    // Vertices 1 (top edge) and 5 (right edge) are boundary vertices joined by
    // an interior diagonal; collapsing it would pinch the surface.
    let diagonal = half_edge_between(&mesh, 1, 5);
    assert_ne!(mesh.half_edges[diagonal].twin, usize::MAX);
    let before = snapshot(&mesh);
    let err = TopologicalSurgery::new(&mut mesh)
        .edge_collapse(diagonal)
        .expect_err("link condition with the vertex at infinity");
    assert!(
        matches!(err, Error::Unsupported { .. }),
        "unexpected: {err}"
    );
    assert_eq!(snapshot(&mesh), before);
}

#[test]
fn edge_collapse_on_a_torus_preserves_genus() {
    let fixture = grid_torus(4, 4);
    let mut collapsed = 0;
    for he in 0..96 {
        let mut mesh = build(&fixture);
        let before = snapshot(&mesh);
        match TopologicalSurgery::new(&mut mesh).edge_collapse(he) {
            Ok(_) => {
                collapsed += 1;
                assert_eq!(counts(&mesh), (15, 45, 30), "half-edge {he}");
                assert_eq!(mesh.euler_characteristic(), 0);
                assert_eq!(mesh.genus(), 1);
                assert!(mesh.is_closed());
                assert_eq!(homology_of(&mesh), vec![1, 2, 1], "half-edge {he}");
                assert_half_edge_invariants(&mesh);
            }
            Err(e) => {
                assert!(
                    matches!(e, Error::Unsupported { .. }),
                    "half-edge {he}: {e}"
                );
                assert_eq!(snapshot(&mesh), before);
            }
        }
    }
    assert!(collapsed > 0, "some torus edges must be collapsible");
}

#[test]
fn split_and_collapse_random_walk_keeps_a_valid_sphere() {
    let mut mesh = build(&icosahedron());
    let mut rng = Lcg(0x1234_5678_9abc_def0);
    let (mut splits, mut collapses, mut rejected) = (0, 0, 0);
    for step in 0..400 {
        let he = rng.below(mesh.half_edges.len());
        let prefer_collapse = mesh.num_vertices() > 40 || rng.below(2) == 0;
        let before = snapshot(&mesh);
        if prefer_collapse {
            match TopologicalSurgery::new(&mut mesh).edge_collapse(he) {
                Ok(_) => collapses += 1,
                Err(_) => {
                    rejected += 1;
                    assert_eq!(
                        snapshot(&mesh),
                        before,
                        "step {step}: failed collapse changed the mesh"
                    );
                    continue;
                }
            }
        } else {
            TopologicalSurgery::new(&mut mesh)
                .edge_split(he)
                .expect("a split of a valid mesh never fails");
            splits += 1;
        }
        assert_eq!(mesh.euler_characteristic(), 2, "step {step}");
        assert!(mesh.is_closed(), "step {step}");
        assert_eq!(mesh.genus(), 0);
        assert_half_edge_invariants(&mesh);
    }
    assert!(
        splits > 50 && collapses > 50,
        "splits {splits}, collapses {collapses}"
    );
    assert!(rejected < 400);
}

#[test]
fn split_and_collapse_random_walk_keeps_a_valid_disc() {
    // A disc: chi = 1 and exactly one boundary loop must survive every step,
    // which exercises the link condition at boundary vertices and edges.
    let mut mesh = build(&grid(6, 6, true));
    let mut rng = Lcg(0x0dd5_ba11_cafe_f00d);
    let (mut splits, mut collapses) = (0, 0);
    for step in 0..500 {
        let he = rng.below(mesh.half_edges.len());
        let before = snapshot(&mesh);
        if mesh.num_vertices() > 60 || rng.below(2) == 0 {
            match TopologicalSurgery::new(&mut mesh).edge_collapse(he) {
                Ok(_) => collapses += 1,
                Err(e) => {
                    assert!(matches!(e, Error::Unsupported { .. }), "step {step}: {e}");
                    assert_eq!(snapshot(&mesh), before, "step {step}");
                    continue;
                }
            }
        } else {
            TopologicalSurgery::new(&mut mesh)
                .edge_split(he)
                .expect("split");
            splits += 1;
        }
        assert_eq!(mesh.euler_characteristic(), 1, "step {step}");
        assert_eq!(extract_boundary_loops(&mesh).len(), 1, "step {step}");
        assert_eq!(components(&mesh), 1, "step {step}");
        assert_half_edge_invariants(&mesh);
    }
    assert!(
        splits > 50 && collapses > 50,
        "splits {splits}, collapses {collapses}"
    );
}

#[test]
fn split_and_collapse_random_walk_keeps_a_torus() {
    let mut mesh = build(&grid_torus(5, 5));
    let mut rng = Lcg(0x7011_5ca1_ab1e_0001);
    let (mut splits, mut collapses) = (0, 0);
    for step in 0..300 {
        let he = rng.below(mesh.half_edges.len());
        if mesh.num_vertices() > 50 || rng.below(2) == 0 {
            if TopologicalSurgery::new(&mut mesh).edge_collapse(he).is_ok() {
                collapses += 1;
            } else {
                continue;
            }
        } else {
            TopologicalSurgery::new(&mut mesh)
                .edge_split(he)
                .expect("split");
            splits += 1;
        }
        assert_eq!(mesh.euler_characteristic(), 0, "step {step}");
        assert!(mesh.is_closed(), "step {step}");
        if step % 25 == 0 {
            assert_eq!(homology_of(&mesh), vec![1, 2, 1], "step {step}");
        }
        assert_half_edge_invariants(&mesh);
    }
    assert!(
        splits > 30 && collapses > 30,
        "splits {splits}, collapses {collapses}"
    );
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
}

#[test]
fn surgery_tolerates_isolated_vertices_and_rotated_face_starts() {
    // An octahedron plus a vertex that no face uses: collapsing an edge shifts
    // the isolated vertex's index down by one and keeps its position.
    let (mut positions, triangles) = octahedron();
    positions.push([7.0, 8.0, 9.0]);
    let mut mesh = HalfEdgeMesh::from_triangles(positions, &triangles);
    // Start every face loop at a different half-edge than the first one.
    for face in &mut mesh.faces {
        face.start_he = mesh.half_edges[face.start_he].next;
    }
    let isolated = 6;
    assert_eq!(mesh.vertex_he[isolated], usize::MAX);
    let split_vertex = TopologicalSurgery::new(&mut mesh)
        .edge_split(0)
        .expect("split with an isolated vertex present");
    assert_eq!(
        split_vertex, 7,
        "the new vertex goes after the isolated one"
    );
    assert_eq!(mesh.vertices[isolated], [7.0, 8.0, 9.0]);
    assert_eq!(mesh.vertex_he[isolated], usize::MAX);
    let survivor = TopologicalSurgery::new(&mut mesh)
        .edge_collapse(0)
        .expect("collapse with an isolated vertex present");
    assert!(survivor < mesh.num_vertices());
    let moved = mesh
        .vertices
        .iter()
        .position(|p| *p == [7.0, 8.0, 9.0])
        .expect("the isolated vertex is kept");
    assert_eq!(mesh.vertex_he[moved], usize::MAX);
    assert_eq!(
        mesh.num_faces(),
        8,
        "split then collapse restores the face count"
    );
    let (surface, _) = Surface::from_mesh(&mesh).expect("still a valid surface");
    assert_eq!(surface.euler_characteristic(), 2);
}

// ---------------------------------------------------------------------------
// boundary loops
// ---------------------------------------------------------------------------

#[test]
fn boundary_loops_are_deterministic_and_oriented() {
    let mesh = build(&tube(3, 4));
    let loops = TopologicalSurgery::new(&mut mesh.clone())
        .boundary_loops()
        .expect("valid mesh");
    assert_eq!(loops.len(), 2);
    assert!(loops.iter().all(|l| l.len() == 4));
    for _ in 0..5 {
        let again = TopologicalSurgery::new(&mut mesh.clone())
            .boundary_loops()
            .expect("valid mesh");
        assert_eq!(again, loops);
    }
    assert_eq!(extract_boundary_loops(&mesh), loops);
    // Boundary direction: the surface lies on the left, so the reverse of every
    // boundary edge is not an edge of any face.
    let (surface, _) = Surface::from_mesh(&mesh).expect("valid");
    let directed = surface.directed_edges();
    for cycle in &loops {
        for i in 0..cycle.len() {
            let (u, v) = (cycle[i], cycle[(i + 1) % cycle.len()]);
            assert!(directed.contains_key(&(u, v)));
            assert!(!directed.contains_key(&(v, u)));
        }
    }
    let closed = build(&octahedron());
    assert!(extract_boundary_loops(&closed).is_empty());
}

// ---------------------------------------------------------------------------
// attach_handle
// ---------------------------------------------------------------------------

#[test]
fn attach_handle_turns_a_tube_into_a_torus() {
    for segments in 3..=8 {
        let fixture = tube(3, segments);
        let mut mesh = build(&fixture);
        let faces_before = mesh.num_faces();
        let vertices_before = mesh.num_vertices();
        assert_eq!(mesh.euler_characteristic(), 0, "a tube has chi = 0");
        TopologicalSurgery::new(&mut mesh)
            .attach_handle(0, 1)
            .unwrap_or_else(|e| panic!("{segments} segments: {e}"));
        assert_eq!(mesh.num_vertices(), vertices_before, "no new vertices");
        assert_eq!(mesh.num_faces(), faces_before + 2 * segments);
        assert_eq!(mesh.euler_characteristic(), 0);
        assert!(mesh.is_closed(), "{segments} segments");
        assert_eq!(mesh.genus(), 1);
        assert_eq!(homology_of(&mesh), vec![1, 2, 1], "{segments} segments");
        assert_half_edge_invariants(&mesh);
    }
}

#[test]
fn attach_handle_accepts_the_loops_in_either_order() {
    let fixture = tube(3, 5);
    let mut forward = build(&fixture);
    let mut backward = build(&fixture);
    TopologicalSurgery::new(&mut forward)
        .attach_handle(0, 1)
        .expect("0, 1");
    TopologicalSurgery::new(&mut backward)
        .attach_handle(1, 0)
        .expect("1, 0");
    for mesh in [&forward, &backward] {
        assert!(mesh.is_closed());
        assert_eq!(mesh.euler_characteristic(), 0);
        assert_eq!(homology_of(mesh), vec![1, 2, 1]);
    }
}

#[test]
fn attach_handle_joins_loops_of_different_lengths() {
    // An icosahedron with one triangular hole and one quadrilateral hole around
    // the antipode of the triangular one.
    let ico = icosahedron();
    let single = 0usize;
    let far = antipodal_face(&ico, single);
    let neighbour = (0..20)
        .find(|&f| f != far && ico.1[f].iter().filter(|v| ico.1[far].contains(v)).count() == 2)
        .expect("a face next to the antipodal face");
    let fixture = without_faces(&ico, &[single, far, neighbour]);
    let mut mesh = build(&fixture);
    let mut lengths: Vec<usize> = extract_boundary_loops(&mesh).iter().map(Vec::len).collect();
    lengths.sort_unstable();
    assert_eq!(lengths, vec![3, 4]);
    assert_eq!(mesh.euler_characteristic(), 0);
    let faces_before = mesh.num_faces();
    TopologicalSurgery::new(&mut mesh)
        .attach_handle(0, 1)
        .expect("unequal loops");
    assert_eq!(mesh.num_faces(), faces_before + 7, "3 + 4 strip triangles");
    assert!(mesh.is_closed());
    assert_eq!(mesh.euler_characteristic(), 0);
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
    assert_half_edge_invariants(&mesh);
}

/// Index of the face of the icosahedron whose vertices are the antipodes of
/// the vertices of face `face`.
fn antipodal_face(ico: &Fixture, face: usize) -> usize {
    let antipode = |v: usize| {
        let p = ico.0[v];
        (0..ico.0.len())
            .find(|&w| (0..3).all(|k| (ico.0[w][k] + p[k]).abs() < 1e-9))
            .expect("the icosahedron is centrally symmetric")
    };
    let mut target: Vec<usize> = ico.1[face].iter().map(|&v| antipode(v)).collect();
    target.sort_unstable();
    (0..ico.1.len())
        .find(|&f| {
            let mut vertices = ico.1[f].to_vec();
            vertices.sort_unstable();
            vertices == target
        })
        .expect("antipodal face")
}

#[test]
fn attach_handle_on_a_sphere_with_two_holes() {
    // Icosahedron without two antipodal faces. (Two opposite faces of an
    // octahedron are joined by too many existing edges to be stitched without
    // creating a double edge.)
    let ico = icosahedron();
    let opposite = antipodal_face(&ico, 0);
    let mut mesh = build(&without_faces(&ico, &[0, opposite]));
    assert_eq!(mesh.euler_characteristic(), 0);
    assert_eq!(extract_boundary_loops(&mesh).len(), 2);
    TopologicalSurgery::new(&mut mesh)
        .attach_handle(0, 1)
        .expect("sphere with two holes");
    assert!(mesh.is_closed());
    assert_eq!(mesh.euler_characteristic(), 0);
    assert_eq!(mesh.genus(), 1);
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
    assert_half_edge_invariants(&mesh);
}

/// A disc (8 x 8 grid) with triangular holes punched at the given grid vertices.
fn disc_with_holes(holes: &[(usize, usize)]) -> Fixture {
    let disc = grid(8, 8, true);
    let removed: Vec<usize> = holes
        .iter()
        .map(|&(r, c)| {
            let vertex = r * 8 + c;
            disc.1
                .iter()
                .position(|t| t.contains(&vertex))
                .expect("a face at the vertex")
        })
        .collect();
    without_faces(&disc, &removed)
}

#[test]
fn attach_handle_between_two_holes_of_a_disc() {
    let mut mesh = build(&disc_with_holes(&[(2, 2), (5, 5)]));
    assert_eq!(mesh.euler_characteristic(), -1);
    let loops = extract_boundary_loops(&mesh);
    assert_eq!(loops.len(), 3);
    let holes: Vec<usize> = (0..3).filter(|&i| loops[i].len() == 3).collect();
    assert_eq!(holes.len(), 2, "two triangular holes");
    TopologicalSurgery::new(&mut mesh)
        .attach_handle(holes[0], holes[1])
        .expect("handle between the holes");
    // A torus with one hole: chi = -1, one boundary loop, beta = (1, 2, 0).
    assert_eq!(mesh.euler_characteristic(), -1);
    assert_eq!(extract_boundary_loops(&mesh).len(), 1);
    assert_eq!(homology_of(&mesh), vec![1, 2, 0]);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn two_handles_make_a_genus_two_surface_with_one_hole() {
    let mut mesh = build(&disc_with_holes(&[(2, 2), (5, 5), (2, 5), (5, 2)]));
    assert_eq!(mesh.euler_characteristic(), -3);
    for _ in 0..2 {
        let loops = extract_boundary_loops(&mesh);
        let holes: Vec<usize> = (0..loops.len()).filter(|&i| loops[i].len() == 3).collect();
        TopologicalSurgery::new(&mut mesh)
            .attach_handle(holes[0], holes[1])
            .expect("handle between two holes");
        assert_half_edge_invariants(&mesh);
    }
    // genus 2 with one boundary loop: chi = 2 - 4 - 1 = -3, beta_1 = 4.
    assert_eq!(mesh.euler_characteristic(), -3);
    assert_eq!(extract_boundary_loops(&mesh).len(), 1);
    assert_eq!(homology_of(&mesh), vec![1, 4, 0]);
}

#[test]
fn attach_handle_is_deterministic() {
    let fixture = tube(4, 6);
    let mut first = build(&fixture);
    let mut second = build(&fixture);
    TopologicalSurgery::new(&mut first)
        .attach_handle(0, 1)
        .expect("first");
    TopologicalSurgery::new(&mut second)
        .attach_handle(0, 1)
        .expect("second");
    assert_eq!(snapshot(&first), snapshot(&second));
}

#[test]
fn attach_handle_rejects_bad_requests_and_keeps_the_mesh() {
    // A single triangle has one boundary loop, so there is no second loop.
    let triangle = (
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
        vec![[0, 1, 2]],
    );
    let mut mesh = build(&triangle);
    let before = snapshot(&mesh);
    let err = TopologicalSurgery::new(&mut mesh)
        .attach_handle(0, 1)
        .expect_err("only one loop");
    assert!(
        matches!(err, Error::IndexOutOfBounds { index: 1, len: 1 }),
        "{err}"
    );
    let err = TopologicalSurgery::new(&mut mesh)
        .attach_handle(0, 0)
        .expect_err("same loop twice");
    assert!(matches!(err, Error::Unsupported { .. }), "{err}");
    assert_eq!(snapshot(&mesh), before);

    // A closed mesh has no boundary loops at all.
    let mut closed = build(&octahedron());
    let err = TopologicalSurgery::new(&mut closed)
        .attach_handle(0, 1)
        .expect_err("closed mesh");
    assert!(
        matches!(err, Error::IndexOutOfBounds { index: 0, len: 0 }),
        "{err}"
    );

    // Loops in different components cannot form a handle.
    let two_triangles = disjoint_union(&triangle, &triangle);
    let mut apart = build(&two_triangles);
    let before = snapshot(&apart);
    let err = TopologicalSurgery::new(&mut apart)
        .attach_handle(0, 1)
        .expect_err("different components");
    assert!(matches!(err, Error::Unsupported { .. }), "{err}");
    assert_eq!(snapshot(&apart), before);
}

// ---------------------------------------------------------------------------
// delete_handle
// ---------------------------------------------------------------------------

#[test]
fn delete_handle_turns_a_grid_torus_into_a_sphere() {
    let mut mesh = build(&grid_torus(4, 4));
    let (v0, e0, f0) = counts(&mesh);
    assert_eq!(mesh.euler_characteristic(), 0);
    TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect("torus has a handle");
    let (v1, e1, f1) = counts(&mesh);
    // k loop vertices are duplicated, two cone vertices are added.
    let k = v1 - v0 - 2;
    assert!(k >= 3, "cut loop of length {k}");
    assert_eq!(e1, e0 + k + 2 * k, "loop edges duplicated plus 2k spokes");
    assert_eq!(f1, f0 + 2 * k);
    assert_eq!(mesh.euler_characteristic(), 2);
    assert!(mesh.is_closed());
    assert_eq!(mesh.genus(), 0);
    assert_eq!(components(&mesh), 1);
    assert_eq!(homology_of(&mesh), vec![1, 0, 1]);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn delete_handle_on_the_seven_vertex_torus() {
    let mut mesh = build(&seven_vertex_torus());
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
    TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect("torus has a handle");
    assert_eq!(mesh.euler_characteristic(), 2);
    assert!(mesh.is_closed());
    assert_eq!(homology_of(&mesh), vec![1, 0, 1]);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn delete_handle_reduces_genus_two_step_by_step() {
    let double_torus = connected_sum(&grid_torus(4, 4), &grid_torus(4, 4));
    let mut mesh = build(&double_torus);
    assert!(
        Surface::from_mesh(&mesh).is_ok(),
        "connected sum must be a manifold"
    );
    assert_eq!(mesh.euler_characteristic(), -2);
    assert_eq!(homology_of(&mesh), vec![1, 4, 1]);
    for expected_chi in [0, 2] {
        TopologicalSurgery::new(&mut mesh)
            .delete_handle()
            .expect("handle to delete");
        assert_eq!(mesh.euler_characteristic(), expected_chi);
        assert!(mesh.is_closed());
        assert_eq!(components(&mesh), 1);
        assert_half_edge_invariants(&mesh);
    }
    assert_eq!(homology_of(&mesh), vec![1, 0, 1]);
    let before = snapshot(&mesh);
    let err = TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect_err("a sphere has no handle");
    assert!(matches!(err, Error::Unsupported { .. }), "{err}");
    assert_eq!(snapshot(&mesh), before);
}

#[test]
fn delete_handle_after_random_remeshing() {
    let mut mesh = build(&grid_torus(5, 5));
    let mut rng = Lcg(0xa11c_e5ea_5011_0007);
    for _ in 0..120 {
        let he = rng.below(mesh.half_edges.len());
        if mesh.num_vertices() > 45 || rng.below(2) == 0 {
            let _ = TopologicalSurgery::new(&mut mesh).edge_collapse(he);
        } else {
            TopologicalSurgery::new(&mut mesh)
                .edge_split(he)
                .expect("split");
        }
    }
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
    TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect("the remeshed torus still has a handle");
    assert_eq!(mesh.euler_characteristic(), 2);
    assert!(mesh.is_closed());
    assert_eq!(homology_of(&mesh), vec![1, 0, 1]);
    assert_half_edge_invariants(&mesh);
}

#[test]
fn delete_handle_leaves_other_components_alone() {
    let torus = grid_torus(4, 4);
    let sphere = octahedron();
    let fixture = disjoint_union(&torus, &sphere);
    let mut mesh = build(&fixture);
    assert_eq!(components(&mesh), 2);
    assert_eq!(mesh.euler_characteristic(), 2);
    TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect("the torus component has a handle");
    assert_eq!(components(&mesh), 2);
    assert_eq!(mesh.euler_characteristic(), 4, "two spheres");
    assert!(mesh.is_closed());
    // Vertex indices are preserved, so the octahedron keeps its place.
    for (i, p) in sphere.0.iter().enumerate() {
        let moved = mesh.vertices[torus.0.len() + i];
        assert_eq!([moved[0] - 20.0, moved[1], moved[2]], *p);
    }
    assert_half_edge_invariants(&mesh);
}

#[test]
fn delete_handle_rejects_surfaces_without_a_closed_handle() {
    for (name, fixture) in [
        ("tetrahedron", tetrahedron()),
        ("octahedron", octahedron()),
        ("tube", tube(3, 4)),
        ("grid", grid(4, 4, true)),
    ] {
        let mut mesh = build(&fixture);
        let before = snapshot(&mesh);
        let err = TopologicalSurgery::new(&mut mesh)
            .delete_handle()
            .expect_err(name);
        assert!(matches!(err, Error::Unsupported { .. }), "{name}: {err}");
        assert_eq!(snapshot(&mesh), before, "{name} changed");
    }
}

#[test]
fn attach_then_delete_handle_round_trips_the_topology() {
    let mut mesh = build(&tube(3, 4));
    assert_eq!(homology_of(&mesh), vec![1, 1, 0]);
    TopologicalSurgery::new(&mut mesh)
        .attach_handle(0, 1)
        .expect("attach");
    assert_eq!(homology_of(&mesh), vec![1, 2, 1]);
    TopologicalSurgery::new(&mut mesh)
        .delete_handle()
        .expect("delete");
    assert_eq!(mesh.euler_characteristic(), 2);
    assert!(mesh.is_closed());
    assert_eq!(homology_of(&mesh), vec![1, 0, 1]);
    assert_half_edge_invariants(&mesh);
}

// ---------------------------------------------------------------------------
// invalid input meshes
// ---------------------------------------------------------------------------

fn assert_rejected_untouched(mesh: &mut HalfEdgeMesh, what: &str) {
    let before = snapshot(mesh);
    let err = TopologicalSurgery::new(mesh).edge_split(0).expect_err(what);
    assert!(matches!(err, Error::InvalidMesh { .. }), "{what}: {err}");
    let err = TopologicalSurgery::new(mesh)
        .edge_collapse(0)
        .expect_err(what);
    assert!(matches!(err, Error::InvalidMesh { .. }), "{what}: {err}");
    let err = TopologicalSurgery::new(mesh)
        .delete_handle()
        .expect_err(what);
    assert!(matches!(err, Error::InvalidMesh { .. }), "{what}: {err}");
    assert_eq!(snapshot(mesh), before, "{what} changed the mesh");
}

fn points(n: usize) -> Vec<[f64; 3]> {
    (0..n).map(|i| [i as f64, (i * i) as f64, 0.0]).collect()
}

#[test]
fn surgery_rejects_invalid_meshes() {
    // Three faces on one edge (non-manifold edge).
    let mut fan = HalfEdgeMesh::from_triangles(points(5), &[[0, 1, 2], [0, 1, 3], [0, 1, 4]]);
    assert_rejected_untouched(&mut fan, "three faces on an edge");
    // Inconsistent orientation: both faces traverse 0 -> 1.
    let mut flipped = HalfEdgeMesh::from_triangles(points(4), &[[0, 1, 2], [0, 1, 3]]);
    assert_rejected_untouched(&mut flipped, "inconsistent orientation");
    // Bow-tie: two faces meeting at a single vertex.
    let mut bowtie = HalfEdgeMesh::from_triangles(points(5), &[[0, 1, 2], [0, 3, 4]]);
    assert_rejected_untouched(&mut bowtie, "bow-tie vertex");
    // The same vertex set on two opposite faces.
    let mut pillow = HalfEdgeMesh::from_triangles(points(3), &[[0, 1, 2], [0, 2, 1]]);
    assert_rejected_untouched(&mut pillow, "duplicate face");
    // Degenerate face.
    let mut degenerate = HalfEdgeMesh::from_triangles(points(3), &[[0, 0, 1]]);
    assert_rejected_untouched(&mut degenerate, "degenerate face");
    // No faces at all.
    let mut empty = HalfEdgeMesh::new();
    let err = TopologicalSurgery::new(&mut empty)
        .edge_split(0)
        .expect_err("empty mesh");
    assert!(matches!(err, Error::InvalidMesh { .. }), "{err}");
}

#[test]
fn surgery_rejects_non_triangular_faces_and_broken_pointers() {
    use super::types::Face;
    // One quadrilateral face.
    let mut quad = HalfEdgeMesh::new();
    quad.vertices = points(4);
    quad.vertex_he = vec![0, 1, 2, 3];
    for i in 0..4 {
        let mut he = HalfEdge::new(i);
        he.next = (i + 1) % 4;
        he.prev = (i + 3) % 4;
        he.face = 0;
        quad.half_edges.push(he);
    }
    quad.faces.push(Face {
        start_he: 0,
        is_outer: false,
    });
    let err = TopologicalSurgery::new(&mut quad)
        .edge_split(0)
        .expect_err("quad face");
    assert!(matches!(err, Error::InvalidMesh { .. }), "{err}");

    // A `next` pointer that leaves the half-edge array.
    let mut broken = build(&tetrahedron());
    broken.half_edges[1].next = 999;
    let err = TopologicalSurgery::new(&mut broken)
        .edge_collapse(0)
        .expect_err("dangling pointer");
    assert!(matches!(err, Error::InvalidMesh { .. }), "{err}");
}
