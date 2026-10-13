// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the hashed vertex welding behind `remove_duplicate_vertices` and
//! `suggest_improvements`.

use super::functions::{remove_duplicate_vertices, suggest_improvements, weld_vertices};
use super::types::{QualitySuggestion, QualityThresholds, TriangleMesh};

/// Small deterministic xorshift generator (the crate has no `rand` dependency).
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

    /// Uniform value in `[0, 1)`.
    fn unit(&mut self) -> f32 {
        (self.next_u64() >> 40) as f32 / (1u64 << 24) as f32
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Reference implementation: the linear first-match scan with the same arithmetic.
fn reference_weld(verts: &[[f32; 3]], tol: f32) -> (Vec<[f32; 3]>, Vec<u32>) {
    let tol_sq = f64::from(tol) * f64::from(tol);
    let mut unique: Vec<[f32; 3]> = Vec::new();
    let mut remap = Vec::new();
    for v in verts {
        let found = unique.iter().position(|u| {
            let ex = f64::from(v[0]) - f64::from(u[0]);
            let ey = f64::from(v[1]) - f64::from(u[1]);
            let ez = f64::from(v[2]) - f64::from(u[2]);
            ex * ex + ey * ey + ez * ez <= tol_sq
        });
        match found {
            Some(i) => remap.push(i as u32),
            None => {
                remap.push(unique.len() as u32);
                unique.push(*v);
            }
        }
    }
    (unique, remap)
}

fn grid_vertices(n: usize) -> Vec<[f32; 3]> {
    (0..n * n)
        .map(|i| [(i % n) as f32, (i / n) as f32, 0.0])
        .collect()
}

/// A grid whose vertices are all present twice, with triangles on the second copy.
fn duplicated_grid_mesh(n: usize) -> TriangleMesh {
    let mut verts = grid_vertices(n);
    verts.extend(grid_vertices(n));
    let offset = (n * n) as u32;
    let mut tris = Vec::new();
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            let a = offset + (j * n + i) as u32;
            let b = a + 1;
            let c = a + n as u32;
            let d = c + 1;
            tris.push([a, b, d]);
            tris.push([a, d, c]);
        }
    }
    TriangleMesh::from_raw(verts, tris)
}

#[test]
fn duplicated_200x200_grid_has_40000_unique_vertices() {
    let mesh = duplicated_grid_mesh(200);
    assert_eq!(mesh.vertices.len(), 80_000);
    let (unique, remap) = weld_vertices(&mesh.vertices, 1e-6);
    assert_eq!(unique.len(), 40_000);
    assert_eq!(remap.len(), 80_000);
    for i in 0..40_000usize {
        assert_eq!(remap[i], i as u32, "first copy keeps its index");
        assert_eq!(
            remap[40_000 + i],
            i as u32,
            "second copy folds onto the first"
        );
    }
    assert_eq!(unique, grid_vertices(200));
}

#[test]
fn remove_duplicate_vertices_keeps_triangle_geometry() {
    let mesh = duplicated_grid_mesh(60);
    let welded = remove_duplicate_vertices(&mesh, 1e-6);
    assert_eq!(welded.vertex_count(), 60 * 60);
    assert_eq!(welded.triangle_count(), mesh.triangle_count());
    for (old, new) in mesh.triangles.iter().zip(&welded.triangles) {
        for k in 0..3 {
            assert!((new[k] as usize) < welded.vertices.len());
            assert_eq!(
                welded.vertices[new[k] as usize], mesh.vertices[old[k] as usize],
                "triangle corner moved"
            );
        }
    }
}

#[test]
fn suggest_improvements_counts_duplicates_on_a_large_mesh() {
    let mesh = duplicated_grid_mesh(200);
    let suggestions = suggest_improvements(&mesh, &QualityThresholds::default());
    let reported = suggestions.iter().find_map(|s| match s {
        QualitySuggestion::RemoveDuplicateVertices {
            estimated_duplicates,
        } => Some(*estimated_duplicates),
        _ => None,
    });
    assert_eq!(reported, Some(40_000), "{suggestions:?}");
}

#[test]
fn suggest_improvements_reports_no_duplicates_for_a_clean_mesh() {
    let mesh = TriangleMesh::from_raw(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]],
        vec![[0, 1, 2]],
    );
    let suggestions = suggest_improvements(&mesh, &QualityThresholds::default());
    assert!(
        !suggestions
            .iter()
            .any(|s| matches!(s, QualitySuggestion::RemoveDuplicateVertices { .. })),
        "{suggestions:?}"
    );
}

#[test]
fn matches_the_linear_first_match_scan_on_random_clouds() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    for &tol in &[0.02_f32, 0.08, 0.3] {
        let verts: Vec<[f32; 3]> = (0..900)
            .map(|_| [rng.unit(), rng.unit(), rng.unit()])
            .collect();
        let (got_unique, got_remap) = weld_vertices(&verts, tol);
        let (want_unique, want_remap) = reference_weld(&verts, tol);
        assert_eq!(got_unique, want_unique, "tol {tol}");
        assert_eq!(got_remap, want_remap, "tol {tol}");
        assert!(
            got_unique.len() < verts.len() || tol < 0.05,
            "tol {tol} merged nothing"
        );
    }
}

#[test]
fn matches_the_linear_scan_on_overlapping_chains() {
    // Neighbours 0.8 tol apart chain together without being transitive, and the
    // shuffled order makes the lowest matching index differ from the nearest one.
    let tol = 0.05_f32;
    let mut rng = XorShift(12345);
    let mut verts: Vec<[f32; 3]> = (0..200)
        .map(|k| [k as f32 * 0.8 * tol, 2.0, -1.0])
        .collect();
    for i in (1..verts.len()).rev() {
        let j = rng.below(i + 1);
        verts.swap(i, j);
    }
    let (got_unique, got_remap) = weld_vertices(&verts, tol);
    let (want_unique, want_remap) = reference_weld(&verts, tol);
    assert_eq!(got_unique, want_unique);
    assert_eq!(got_remap, want_remap);
}

#[test]
fn lowest_matching_index_wins_across_cells() {
    // B (index 0) and A (index 1) are 1.5 apart, so both stay; P is within 0.75 of
    // both and must join the earlier one, B, even though A's cell is probed first.
    let verts = [[1.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.75, 0.0, 0.0]];
    let (unique, remap) = weld_vertices(&verts, 1.0);
    assert_eq!(unique.len(), 2);
    assert_eq!(remap, vec![0, 1, 0]);
}

#[test]
fn points_across_a_cell_boundary_still_merge() {
    let tol = 1.0;
    let straddle_positive = [[0.9999_f32, 0.0, 0.0], [1.0001, 0.0, 0.0]];
    assert_eq!(weld_vertices(&straddle_positive, tol).0.len(), 1);
    let straddle_zero = [[-0.0001_f32, 0.0, 0.0], [0.0001, 0.0, 0.0]];
    assert_eq!(weld_vertices(&straddle_zero, tol).0.len(), 1);
    let corner = [[0.9_f32, 0.9, 0.9], [1.2, 1.2, 1.2]];
    assert_eq!(weld_vertices(&corner, tol).0.len(), 1);
    let negative_corner = [[-0.9_f32, -0.9, -0.9], [-1.2, -1.2, -1.2]];
    assert_eq!(weld_vertices(&negative_corner, tol).0.len(), 1);
}

#[test]
fn tolerance_is_a_euclidean_inclusive_bound() {
    let tol = 1.0;
    assert_eq!(weld_vertices(&[[0.0; 3], [1.0, 0.0, 0.0]], tol).0.len(), 1);
    assert_eq!(
        weld_vertices(&[[0.0; 3], [1.001, 0.0, 0.0]], tol).0.len(),
        2
    );
    // Inside the bounding cube of side 2 * tol but outside the sphere.
    assert_eq!(weld_vertices(&[[0.0; 3], [0.6, 0.6, 0.6]], tol).0.len(), 2);
    assert_eq!(
        weld_vertices(&[[0.0; 3], [0.55, 0.55, 0.55]], tol).0.len(),
        1
    );
}

#[test]
fn zero_negative_and_nan_tolerance_merge_only_identical_points() {
    let verts = [
        [0.0, 0.0, 0.0],
        [-0.0, 0.0, -0.0],
        [1e-9, 0.0, 0.0],
        [f32::NAN, 0.0, 0.0],
        [f32::NAN, 0.0, 0.0],
        [f32::INFINITY, 1.0, 1.0],
        [f32::INFINITY, 1.0, 1.0],
        [0.0, 0.0, 0.0],
    ];
    for tol in [0.0_f32, -1.0, f32::NAN] {
        let (unique, remap) = weld_vertices(&verts, tol);
        assert_eq!(unique.len(), 6, "tol {tol}");
        assert_eq!(remap, vec![0, 0, 1, 2, 3, 4, 5, 0], "tol {tol}");
    }
}

#[test]
fn non_finite_vertices_never_merge() {
    let verts = [
        [0.0, 0.0, 0.0],
        [f32::INFINITY, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [f32::INFINITY, 0.0, 0.0],
        [f32::NAN, f32::NAN, f32::NAN],
        [f32::NEG_INFINITY, 1.0, 1.0],
    ];
    let (unique, remap) = weld_vertices(&verts, 0.5);
    assert_eq!(unique.len(), 5);
    assert_eq!(remap, vec![0, 1, 0, 2, 3, 4]);
}

#[test]
fn extreme_coordinates_and_tolerances_do_not_panic() {
    let verts = [
        [1e30_f32, 0.0, 0.0],
        [1e30, 0.0, 0.0],
        [-1e30, 5.0, 0.0],
        [-1e30, 5.0, 0.0],
        [0.0, 0.0, 0.0],
    ];
    let (unique, remap) = weld_vertices(&verts, 1e-12);
    assert_eq!(unique.len(), 3);
    assert_eq!(remap, vec![0, 0, 1, 1, 2]);
    // An enormous tolerance merges every finite point into the first one.
    let (unique, remap) = weld_vertices(&verts, f32::MAX);
    assert_eq!(unique.len(), 1);
    assert_eq!(remap, vec![0; 5]);
    let (unique, _) = weld_vertices(&verts, f32::INFINITY);
    assert_eq!(unique.len(), 1);
}

#[test]
fn empty_input_is_fine() {
    let (unique, remap) = weld_vertices(&[], 1e-3);
    assert!(unique.is_empty() && remap.is_empty());
    let mesh = TriangleMesh::new();
    assert_eq!(remove_duplicate_vertices(&mesh, 1e-3).vertex_count(), 0);
}

#[test]
fn remove_duplicate_vertices_remaps_triangles() {
    let verts = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let mesh = TriangleMesh::from_raw(verts, vec![[0, 1, 3], [2, 1, 3]]);
    let welded = remove_duplicate_vertices(&mesh, 1e-5);
    assert_eq!(
        welded.vertices,
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    );
    assert_eq!(welded.triangles, vec![[0, 1, 2], [0, 1, 2]]);
}
