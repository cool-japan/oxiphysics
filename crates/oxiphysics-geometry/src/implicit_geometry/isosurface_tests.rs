// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the marching-cubes extraction of [`SdfGrid`] fields
//! ([`extract_isosurface`]).

use super::functions::{extract_isosurface, interp_edge};
use super::types::{IsoMeshResult, SdfGrid, SdfSphere};
use crate::mc_tables::test_support::{
    area_and_volume, degenerate_triangle_count, edge_valence_histogram,
    is_closed_and_consistently_oriented, weld_soup,
};
use crate::mc_tables::{MC_CORNER_OFFSETS, case_triangles};
use std::f64::consts::PI;

fn soup_of(mesh: &IsoMeshResult) -> Vec<[[f64; 3]; 3]> {
    mesh.triangles.iter().map(|t| [t.v0, t.v1, t.v2]).collect()
}

/// A cubic grid of `n` nodes per axis around an origin-centred sphere.  The grid is
/// shifted off the sphere centre by a fraction of the spacing so that no node lies
/// exactly on the surface (a zero node value would put several crossing points on
/// the same node and create zero-area triangles).
fn sphere_grid(n: usize, spacing: f64, radius: f64) -> SdfGrid {
    let half = 0.5 * (n - 1) as f64 * spacing;
    let origin = [-half + 0.013, -half - 0.007, -half + 0.0031];
    SdfGrid::from_sdf(&SdfSphere::new(radius), n, n, n, origin, spacing)
}

#[test]
fn interp_edge_returns_the_zero_crossing() {
    let p0 = [0.0, 0.0, 0.0];
    let p1 = [2.0, 0.0, 0.0];
    let q = interp_edge(p0, -1.0, p1, 3.0);
    assert!((q[0] - 0.5).abs() < 1e-12, "q = {q:?}");
    let q = interp_edge(p0, 3.0, p1, -1.0);
    assert!((q[0] - 1.5).abs() < 1e-12, "q = {q:?}");
}

#[test]
fn interp_edge_returns_the_midpoint_for_a_flat_edge() {
    let p0 = [1.0, 2.0, 3.0];
    let p1 = [3.0, 2.0, 5.0];
    // Equal values and values closer than 1e-12 carry no crossing information.
    for (d0, d1) in [(0.0, 0.0), (-2.0e-13, 3.0e-13), (0.5, 0.5)] {
        let q = interp_edge(p0, d0, p1, d1);
        assert_eq!(q, [2.0, 2.0, 4.0], "d0 = {d0}, d1 = {d1}");
    }
}

#[test]
fn every_nontrivial_cube_case_emits_its_table_triangles() {
    for case in 0..256usize {
        let mut grid = SdfGrid::new(2, 2, 2, [0.0; 3], 1.0);
        for (corner, &[di, dj, dk]) in MC_CORNER_OFFSETS.iter().enumerate() {
            let value = if (case >> corner) & 1 == 1 { -1.0 } else { 1.0 };
            grid.set(di, dj, dk, value);
        }
        let mesh = extract_isosurface(&grid);
        let expected = case_triangles(case).count();
        assert_eq!(mesh.triangles.len(), expected, "case {case}");
        assert_eq!(mesh.vertex_count, 3 * expected, "case {case}");
        if (1..=254).contains(&case) {
            assert!(expected > 0, "case {case} must not be empty");
        }
        for tri in &mesh.triangles {
            for p in [tri.v0, tri.v1, tri.v2] {
                // Crossings of a +-1 field are edge midpoints: exactly one
                // coordinate is 0.5 and the other two are 0 or 1.
                let halves = p.iter().filter(|&&c| (c - 0.5).abs() < 1e-12).count();
                let ends = p
                    .iter()
                    .filter(|&&c| c.abs() < 1e-12 || (c - 1.0).abs() < 1e-12)
                    .count();
                assert_eq!((halves, ends), (1, 2), "case {case}: vertex {p:?}");
            }
        }
    }
}

#[test]
fn sphere_grid_gives_a_closed_outward_mesh_with_the_right_area() {
    let radius = 1.0;
    let grid = sphere_grid(31, 0.1, radius);
    let mesh = extract_isosurface(&grid);
    assert!(mesh.triangles.len() > 1000, "{}", mesh.triangles.len());
    assert_eq!(mesh.vertex_count, 3 * mesh.triangles.len());
    let (vertices, tris) = weld_soup(&soup_of(&mesh));
    assert_eq!(
        edge_valence_histogram(&tris),
        vec![(2, tris.len() * 3 / 2)],
        "every edge must be shared by exactly two triangles"
    );
    assert!(is_closed_and_consistently_oriented(&tris));
    assert_eq!(degenerate_triangle_count(&vertices, &tris), 0);
    let (area, volume) = area_and_volume(&vertices, &tris);
    let exact_area = 4.0 * PI * radius * radius;
    let exact_volume = 4.0 / 3.0 * PI * radius.powi(3);
    assert!(
        (area / exact_area - 1.0).abs() < 0.02,
        "area {area} vs {exact_area}"
    );
    assert!(
        volume > 0.0,
        "triangles must face outwards, volume {volume}"
    );
    assert!(
        (volume / exact_volume - 1.0).abs() < 0.02,
        "volume {volume} vs {exact_volume}"
    );
}

#[test]
fn triangle_normals_point_away_from_the_sphere_centre() {
    let grid = sphere_grid(21, 0.15, 1.0);
    let mesh = extract_isosurface(&grid);
    assert!(!mesh.triangles.is_empty());
    for tri in &mesh.triangles {
        let n = tri.normal();
        let centroid = [
            (tri.v0[0] + tri.v1[0] + tri.v2[0]) / 3.0,
            (tri.v0[1] + tri.v1[1] + tri.v2[1]) / 3.0,
            (tri.v0[2] + tri.v1[2] + tri.v2[2]) / 3.0,
        ];
        let outward = n[0] * centroid[0] + n[1] * centroid[1] + n[2] * centroid[2];
        assert!(outward > 0.0, "inward-facing triangle {tri:?}");
    }
}

#[test]
fn only_edges_with_a_sign_change_produce_vertices() {
    // The old table extrapolated along edges that do not cross the surface; every
    // emitted vertex must lie inside the grid cell it came from, i.e. inside the
    // bounding box of the whole grid.
    let grid = sphere_grid(15, 0.2, 0.9);
    let half = 0.5 * 14.0 * 0.2;
    let mesh = extract_isosurface(&grid);
    for tri in &mesh.triangles {
        for p in [tri.v0, tri.v1, tri.v2] {
            for c in p {
                assert!(c.abs() <= half + 1e-12, "vertex outside the grid: {p:?}");
            }
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 0.9).abs() < 1e-2, "vertex off the sphere: r = {r}");
        }
    }
}

#[test]
fn degenerate_grids_give_empty_meshes() {
    let grid = SdfGrid::new(1, 4, 4, [0.0; 3], 1.0);
    let mesh = extract_isosurface(&grid);
    assert!(mesh.triangles.is_empty());
    assert_eq!(mesh.vertex_count, 0);
}
