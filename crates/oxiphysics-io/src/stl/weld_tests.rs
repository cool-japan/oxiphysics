//! Regression tests for [`weld_vertices`]: first-match semantics, tolerance
//! corner cases, and linear-time behaviour on production-size meshes.

use super::functions_2::{
    count_unique_vertices, merge_and_weld, stl_to_triangle_mesh, weld_vertices,
};
use super::types::StlMesh;

/// Deterministic xorshift generator (no dev-dependency needed).
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next_u64() % bound
    }
}

/// `nx` x `ny` quads in the z = 0 plane, two triangles per quad, every triangle
/// listing its own copies of the corner coordinates (as an STL file does).
fn grid_mesh(nx: usize, ny: usize, spacing: f32) -> StlMesh {
    let mut mesh = StlMesh::new("grid");
    let point = |i: usize, j: usize| [i as f32 * spacing, j as f32 * spacing, 0.0];
    for j in 0..ny {
        for i in 0..nx {
            mesh.add_triangle(point(i, j), point(i + 1, j), point(i, j + 1));
            mesh.add_triangle(point(i + 1, j), point(i + 1, j + 1), point(i, j + 1));
        }
    }
    mesh
}

/// The pre-fix algorithm: linear scan over the canonical vertices, first match.
fn reference_weld(mesh: &StlMesh, tolerance: f32) -> (Vec<[f32; 3]>, Vec<usize>) {
    let mut unique: Vec<[f32; 3]> = Vec::new();
    let mut remap = Vec::new();
    for tri in &mesh.triangles {
        for v in [tri.v0, tri.v1, tri.v2] {
            let found = unique.iter().position(|u| {
                let dx = u[0] - v[0];
                let dy = u[1] - v[1];
                let dz = u[2] - v[2];
                (dx * dx + dy * dy + dz * dz).sqrt() <= tolerance
            });
            match found {
                Some(i) => remap.push(i),
                None => {
                    remap.push(unique.len());
                    unique.push(v);
                }
            }
        }
    }
    (unique, remap)
}

fn flat_remap(triangles: &[[usize; 3]]) -> Vec<usize> {
    triangles.iter().flat_map(|t| t.iter().copied()).collect()
}

#[test]
fn weld_of_a_200k_triangle_grid_is_fast_and_exact() {
    // Regression: the linear scan over the canonical vertices made a mesh of
    // this size take minutes (about 6e10 distance evaluations).
    let (nx, ny) = (320usize, 320usize);
    let mesh = grid_mesh(nx, ny, 0.5);
    assert_eq!(mesh.triangles.len(), 2 * nx * ny);
    let start = std::time::Instant::now();
    let welded = weld_vertices(&mesh, 1e-4);
    let elapsed = start.elapsed();
    assert!(
        elapsed < std::time::Duration::from_secs(15),
        "welding {} triangles took {elapsed:?}",
        mesh.triangles.len()
    );
    assert_eq!(welded.vertices.len(), (nx + 1) * (ny + 1));
    assert_eq!(welded.triangles.len(), mesh.triangles.len());
    for (tri, indices) in mesh.triangles.iter().zip(&welded.triangles) {
        for (corner, &index) in [tri.v0, tri.v1, tri.v2].iter().zip(indices) {
            assert_eq!(&welded.vertices[index], corner);
        }
        assert!(indices[0] != indices[1] && indices[1] != indices[2] && indices[0] != indices[2]);
    }
}

#[test]
fn unique_vertex_count_and_indexed_mesh_follow_the_weld() {
    let (nx, ny) = (40usize, 25usize);
    let mesh = grid_mesh(nx, ny, 0.25);
    assert_eq!(count_unique_vertices(&mesh, 1e-4), (nx + 1) * (ny + 1));
    // Zero tolerance merges the exact duplicates only, which is all of them here.
    assert_eq!(count_unique_vertices(&mesh, 0.0), (nx + 1) * (ny + 1));
    let indexed = stl_to_triangle_mesh(&mesh, 1e-4);
    assert_eq!(indexed.positions.len(), (nx + 1) * (ny + 1));
    assert_eq!(indexed.indices.len(), 2 * nx * ny);
    assert_eq!(indexed.normals.len(), indexed.positions.len());
    // The mesh is flat in z = 0 with counter-clockwise triangles: +z everywhere.
    for normal in &indexed.normals {
        assert!((normal[2] - 1.0).abs() < 1e-5, "normal {normal:?}");
    }
}

#[test]
fn weld_matches_the_linear_scan_on_a_jittered_lattice() {
    // Coordinates are multiples of 2^-10, so the squared distances are exact in
    // both f32 and f64 and the grid and the scan must agree on every corner.
    let unit = 1.0f32 / 1024.0;
    let mut rng = XorShift(0xD1B5_4A32_D192_ED03);
    let mut mesh = StlMesh::new("jitter");
    let random_point = |rng: &mut XorShift| {
        [
            rng.below(24) as f32 * unit,
            rng.below(24) as f32 * unit,
            rng.below(24) as f32 * unit,
        ]
    };
    for _ in 0..2000 {
        let (a, b, c) = (
            random_point(&mut rng),
            random_point(&mut rng),
            random_point(&mut rng),
        );
        mesh.add_triangle(a, b, c);
    }
    for tolerance in [0.0f32, 1.0 * unit, 1.5 * unit, 2.0 * unit, 4.5 * unit] {
        let welded = weld_vertices(&mesh, tolerance);
        let (want_vertices, want_remap) = reference_weld(&mesh, tolerance);
        assert_eq!(welded.vertices, want_vertices, "tolerance = {tolerance}");
        assert_eq!(
            flat_remap(&welded.triangles),
            want_remap,
            "tolerance = {tolerance}"
        );
    }
}

#[test]
fn weld_merges_into_the_first_canonical_vertex_not_into_merged_corners() {
    let mut mesh = StlMesh::new("chain");
    // a and b are 0.8 apart, b and c are 0.8 apart, a and c are 1.6 apart.
    mesh.add_triangle([0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [1.6, 0.0, 0.0]);
    let welded = weld_vertices(&mesh, 1.0);
    assert_eq!(welded.vertices, vec![[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]]);
    assert_eq!(welded.triangles, vec![[0, 0, 1]]);
}

#[test]
fn weld_picks_the_smallest_matching_canonical_index() {
    let mut mesh = StlMesh::new("two_candidates");
    // u0 and u1 are 1.5 apart (distinct); p lies within 1.0 of both.
    mesh.add_triangle([0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.75, 0.0, 0.0]);
    let welded = weld_vertices(&mesh, 1.0);
    assert_eq!(welded.vertices.len(), 2);
    assert_eq!(welded.triangles, vec![[0, 1, 0]]);
}

#[test]
fn weld_with_zero_or_negative_tolerance_merges_identical_corners_only() {
    let mut mesh = StlMesh::new("exact");
    mesh.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    mesh.add_triangle([-0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 1e-7]);
    for tolerance in [0.0f32, -1.0, f32::NAN] {
        let welded = weld_vertices(&mesh, tolerance);
        assert_eq!(welded.vertices.len(), 4, "tolerance = {tolerance}");
        assert_eq!(welded.triangles, vec![[0, 1, 2], [0, 1, 3]]);
    }
}

#[test]
fn weld_never_merges_non_finite_corners() {
    let mut mesh = StlMesh::new("nan");
    mesh.add_triangle(
        [f32::NAN, 0.0, 0.0],
        [f32::NAN, 0.0, 0.0],
        [f32::INFINITY, 0.0, 0.0],
    );
    mesh.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [f32::INFINITY, 0.0, 0.0]);
    for tolerance in [0.0f32, 1e-3, 10.0] {
        let welded = weld_vertices(&mesh, tolerance);
        // The four non-finite corners stay separate; the two finite origins merge.
        assert_eq!(welded.vertices.len(), 5, "tolerance = {tolerance}");
        assert_eq!(welded.triangles, vec![[0, 1, 2], [3, 3, 4]]);
    }
}

#[test]
fn weld_handles_huge_coordinates_and_tiny_tolerances() {
    let mut mesh = StlMesh::new("huge");
    mesh.add_triangle(
        [1.0e30, -1.0e30, 5.0e29],
        [1.0e30, -1.0e30, 5.0e29],
        [0.0, 0.0, 0.0],
    );
    mesh.add_triangle(
        [f32::MAX, f32::MAX, f32::MAX],
        [0.0, 0.0, 0.0],
        [1.0e-30, 0.0, 0.0],
    );
    for tolerance in [1e-30f32, 1e-12, 1.0] {
        let welded = weld_vertices(&mesh, tolerance);
        assert_eq!(welded.triangles.len(), 2);
        assert_eq!(welded.triangles[0][0], welded.triangles[0][1]);
        assert_ne!(welded.triangles[0][0], welded.triangles[0][2]);
        assert_eq!(welded.triangles[1][1], welded.triangles[0][2]);
    }
}

#[test]
fn weld_of_an_empty_mesh_is_empty() {
    let welded = weld_vertices(&StlMesh::new("empty"), 1e-3);
    assert!(welded.vertices.is_empty());
    assert!(welded.triangles.is_empty());
}

#[test]
fn merge_and_weld_shares_the_common_edge() {
    let mut a = StlMesh::new("a");
    a.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let mut b = StlMesh::new("b");
    b.add_triangle([1.0, 0.0, 1e-6], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
    let welded = merge_and_weld(&a, &b, 1e-4);
    assert_eq!(welded.vertices.len(), 4);
    assert_eq!(welded.triangles, vec![[0, 1, 2], [1, 3, 2]]);
}
