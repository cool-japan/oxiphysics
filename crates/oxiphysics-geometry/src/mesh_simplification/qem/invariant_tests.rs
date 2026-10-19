//! Invariants of the simplifier state.
//!
//! These tests look inside [`Simplifier`]: after a run that is not stopped by a
//! target or an error limit, (1) the adjacency lists must describe exactly the
//! live faces, (2) the boundary flags must be those of the surface that is left,
//! and (3) no valid collapse may remain. The last one only holds if every
//! rejected edge is retried when its neighbourhood changes, so it checks the
//! retry bookkeeping.

use super::*;
use crate::test_meshes::{flat_grid, icosphere};

type Mesh = (Vec<[f64; 3]>, Vec<[usize; 3]>);

fn config(preserve_boundary: bool) -> QemConfig {
    QemConfig {
        target_triangles: 0,
        max_error: 0.0,
        preserve_boundary,
        max_normal_deviation: 0.0,
        preserve_texture_seams: false,
        boundary_weight: 100.0,
    }
}

fn run_to_the_end<'a>(
    mesh: &Mesh,
    config: &'a QemConfig,
    seams: Option<&TextureSeamFlags>,
) -> Simplifier<'a> {
    let mut simplifier = Simplifier::new(&mesh.0, &mesh.1, config, seams);
    simplifier.run();
    simplifier
}

fn assert_consistent(s: &Simplifier<'_>) {
    let alive = (0..s.faces.len()).filter(|&f| s.face_alive[f]).count();
    assert_eq!(alive, s.live_faces, "live face counter");
    for (f, tri) in s.faces.iter().enumerate() {
        if !s.face_alive[f] {
            continue;
        }
        assert!(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
        for &v in tri {
            assert!(s.vertex_alive[v], "live face {f} uses the dead vertex {v}");
            assert!(s.vertex_faces[v].contains(&f), "vertex {v} misses face {f}");
        }
    }
    for (v, faces) in s.vertex_faces.iter().enumerate() {
        if !s.vertex_alive[v] {
            assert!(faces.is_empty(), "dead vertex {v} keeps faces");
            continue;
        }
        let mut sorted = faces.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), faces.len(), "duplicate face at vertex {v}");
        for &f in faces {
            assert!(s.face_alive[f] && s.faces[f].contains(&v));
        }
        let on_boundary = s
            .neighbours(v)
            .into_iter()
            .any(|u| s.faces_with_both(v, u).len() == 1);
        assert_eq!(s.boundary[v], on_boundary, "boundary flag of vertex {v}");
    }
}

fn assert_no_valid_collapse_left(s: &Simplifier<'_>) {
    for a in 0..s.position.len() {
        if !s.vertex_alive[a] {
            continue;
        }
        for b in s.neighbours(a) {
            if b < a {
                continue;
            }
            if let Some(candidate) = s.candidate(a, b) {
                assert!(
                    !matches!(s.assess(&candidate), Verdict::Valid),
                    "the valid collapse of ({a}, {b}) was left over"
                );
            }
        }
    }
}

fn check(mesh: &Mesh, config: &QemConfig, seams: Option<&TextureSeamFlags>) {
    let simplifier = run_to_the_end(mesh, config, seams);
    assert!(simplifier.heap.is_empty());
    assert_consistent(&simplifier);
    assert_no_valid_collapse_left(&simplifier);
}

/// Pseudo-random numbers in `[-1, 1)` without a dependency.
fn jitter(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
}

fn jittered_icosphere(subdivisions: usize, amount: f64) -> Mesh {
    let (vertices, triangles) = icosphere(subdivisions);
    let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
    let vertices = vertices
        .iter()
        .map(|p| {
            let scale = 1.0 + amount * jitter(&mut seed);
            [p[0] * scale, p[1] * scale, p[2] * scale]
        })
        .collect();
    (vertices, triangles)
}

/// A grid with a rectangular hole: the cells `i, j` in `hole` are left out.
fn grid_with_hole(n: usize, hole: std::ops::Range<usize>) -> Mesh {
    let (vertices, triangles) = flat_grid(n, n, 1.0);
    let kept = triangles
        .iter()
        .enumerate()
        .filter(|(index, _)| {
            let cell = index / 2;
            let (i, j) = (cell % n, cell / n);
            !(hole.contains(&i) && hole.contains(&j))
        })
        .map(|(_, tri)| *tri)
        .collect();
    (vertices, kept)
}

#[test]
fn a_closed_regular_surface_ends_with_no_valid_collapse() {
    for subdivisions in 0..=3 {
        check(&icosphere(subdivisions), &config(true), None);
    }
}

#[test]
fn a_closed_irregular_surface_ends_with_no_valid_collapse() {
    for amount in [0.02, 0.1, 0.3] {
        check(&jittered_icosphere(2, amount), &config(true), None);
    }
}

#[test]
fn an_open_surface_ends_with_no_valid_collapse() {
    let grid = flat_grid(8, 8, 1.0);
    check(&grid, &config(true), None);
    check(&grid, &config(false), None);
    let bumpy = (
        grid.0
            .iter()
            .map(|p| [p[0], p[1], 0.4 * (0.7 * p[0]).sin() * (0.5 * p[1]).cos()])
            .collect::<Vec<_>>(),
        grid.1.clone(),
    );
    check(&bumpy, &config(true), None);
    check(&bumpy, &config(false), None);
}

#[test]
fn a_surface_with_a_hole_ends_with_no_valid_collapse() {
    let holed = grid_with_hole(9, 3..6);
    check(&holed, &config(true), None);
    check(&holed, &config(false), None);
}

#[test]
fn limits_do_not_break_the_maximality() {
    let sphere = jittered_icosphere(2, 0.05);
    let limited = QemConfig {
        max_normal_deviation: 0.5,
        ..config(true)
    };
    check(&sphere, &limited, None);
    let tight = QemConfig {
        max_normal_deviation: 0.15,
        ..config(false)
    };
    check(&sphere, &tight, None);
}

#[test]
fn seams_do_not_break_the_maximality() {
    let grid = flat_grid(8, 8, 1.0);
    let mut seams = TextureSeamFlags::new();
    for j in 0..8 {
        seams.mark_seam(j * 9 + 4, (j + 1) * 9 + 4);
    }
    for i in 0..8 {
        seams.mark_seam(4 * 9 + i, 4 * 9 + i + 1);
    }
    let config = QemConfig {
        preserve_texture_seams: true,
        ..config(false)
    };
    check(&grid, &config, Some(&seams));
}

#[test]
fn the_state_is_consistent_at_every_stage_of_a_run() {
    let mesh = jittered_icosphere(2, 0.1);
    for target in [300, 250, 200, 150, 100, 60, 30] {
        let config = QemConfig {
            target_triangles: target,
            ..config(true)
        };
        let simplifier = run_to_the_end(&mesh, &config, None);
        assert_consistent(&simplifier);
        if target >= 100 {
            assert!(simplifier.live_faces <= target, "target {target}");
        }
    }
}

#[test]
fn rejected_edges_are_retried_when_the_neighbourhood_changes() {
    // A face of the octahedron is split by a vertex far above it. The edges of the
    // original face share three neighbours with the new vertex around them (the
    // face is a separating triangle), so the link condition forbids their
    // collapse until the spike is gone.
    let vertices = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.5, 1.5, 1.5],
    ];
    let triangles = vec![
        [0, 2, 6],
        [2, 4, 6],
        [4, 0, 6],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let mesh = (vertices, triangles);
    let config = config(true);
    let simplifier = run_to_the_end(&mesh, &config, None);
    assert_consistent(&simplifier);
    assert_no_valid_collapse_left(&simplifier);
    // The spike and the edges around the separating triangle can all go.
    assert!(simplifier.live_faces < 10);
}
