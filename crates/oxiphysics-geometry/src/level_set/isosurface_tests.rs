// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the marching-cubes isosurface extraction of level-set
//! fields ([`IsosurfaceExtraction`]).

use super::functions::interp_vertex;
use super::types::{IsosurfaceExtraction, LevelSetField, McTriangle};
use crate::mc_tables::test_support::{
    area_and_volume, degenerate_triangle_count, edge_valence_histogram,
    is_closed_and_consistently_oriented, weld_soup,
};
use crate::mc_tables::{MC_CORNER_OFFSETS, case_triangles};
use std::f64::consts::PI;

fn soup_of(triangles: &[McTriangle]) -> Vec<[[f64; 3]; 3]> {
    triangles.iter().map(|t| t.vertices).collect()
}

/// Centre of the test spheres.  It is off the grid lattice, so that no grid node
/// lies exactly on a sphere surface (a zero node value would put several crossing
/// points on the same node and create zero-area triangles).
const SPHERE_CENTRE: [f64; 3] = [0.013, -0.007, 0.0031];

/// A cubic grid of `n` nodes per axis centred on the origin holding the signed
/// distance to a sphere of the given radius.
fn sphere_field(n: usize, dx: f64, radius: f64) -> LevelSetField {
    let half = 0.5 * (n - 1) as f64 * dx;
    let mut field = LevelSetField::new(n, n, n, dx, [-half; 3]);
    field.init_sphere(SPHERE_CENTRE, radius);
    field
}

#[test]
fn interp_vertex_clamps_to_the_edge_and_handles_flat_edges() {
    let a = [0.0, 0.0, 0.0];
    let b = [2.0, 0.0, 0.0];
    assert_eq!(interp_vertex(a, b, -1.0, 1.0, 0.0), [1.0, 0.0, 0.0]);
    // An isovalue outside the range of the end point values is clamped to the edge.
    assert_eq!(interp_vertex(a, b, 1.0, 2.0, 0.0), a);
    assert_eq!(interp_vertex(a, b, 1.0, 2.0, 5.0), b);
    // A constant field carries no crossing information: take the midpoint.
    assert_eq!(interp_vertex(a, b, 3.0, 3.0, 3.0), [1.0, 0.0, 0.0]);
}

#[test]
fn every_nontrivial_cube_case_emits_its_table_triangles() {
    // One 2x2x2 cell per case; corner `c` is inside (-1) when bit `c` is set.
    for case in 0..256usize {
        let mut field = LevelSetField::new(2, 2, 2, 1.0, [0.0; 3]);
        for (corner, &[di, dj, dk]) in MC_CORNER_OFFSETS.iter().enumerate() {
            let value = if (case >> corner) & 1 == 1 { -1.0 } else { 1.0 };
            field.set(di, dj, dk, value);
        }
        let triangles = IsosurfaceExtraction::new(0.0).extract(&field);
        let expected = case_triangles(case).count();
        assert_eq!(
            triangles.len(),
            expected,
            "case {case}: wrong number of triangles"
        );
        if (1..=254).contains(&case) {
            assert!(expected > 0, "case {case} must not be empty");
        }
        for tri in &triangles {
            for p in &tri.vertices {
                // With values +-1 every crossing is the midpoint of a cube edge:
                // exactly one coordinate is 0.5 and the others are 0 or 1.
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
fn single_inside_node_gives_a_closed_octahedron() {
    let mut field = LevelSetField::new(5, 5, 5, 1.0, [0.0; 3]);
    for v in field.data.iter_mut() {
        *v = 1.0;
    }
    field.set(2, 2, 2, -1.0);
    let triangles = IsosurfaceExtraction::new(0.0).extract(&field);
    assert_eq!(triangles.len(), 8, "one triangle per cell around the node");
    let (vertices, tris) = weld_soup(&soup_of(&triangles));
    assert_eq!(vertices.len(), 6, "octahedron has six vertices");
    assert_eq!(edge_valence_histogram(&tris), vec![(2, 12)]);
    assert!(is_closed_and_consistently_oriented(&tris));
    let (area, volume) = area_and_volume(&vertices, &tris);
    // Octahedron with vertices at +-0.5 along the axes through the node (2,2,2).
    let a = 0.5_f64;
    assert!(
        (area - 4.0 * 3.0_f64.sqrt() * a * a).abs() < 1e-12,
        "area {area}"
    );
    assert!(
        (volume - 4.0 / 3.0 * a * a * a).abs() < 1e-12,
        "volume {volume}"
    );
}

#[test]
fn sphere_isosurface_is_closed_outward_and_accurate() {
    let radius = 1.0;
    let field = sphere_field(31, 0.1, radius);
    let triangles = IsosurfaceExtraction::new(0.0).extract(&field);
    assert!(triangles.len() > 1000, "{} triangles", triangles.len());
    let (vertices, tris) = weld_soup(&soup_of(&triangles));
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
    for v in &vertices {
        let d = [
            v[0] - SPHERE_CENTRE[0],
            v[1] - SPHERE_CENTRE[1],
            v[2] - SPHERE_CENTRE[2],
        ];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!((r - radius).abs() < 2e-3, "vertex off the sphere: r = {r}");
    }
}

#[test]
fn neighbouring_cells_share_bit_identical_vertices() {
    let field = sphere_field(17, 0.2, 1.0);
    let triangles = IsosurfaceExtraction::new(0.0).extract(&field);
    let (vertices, _) = weld_soup(&soup_of(&triangles));
    // One surface vertex per grid edge that crosses the isovalue.
    let inside = |ix: usize, iy: usize, iz: usize| field.get(ix, iy, iz) < 0.0;
    let mut crossed = 0usize;
    for iz in 0..field.nz {
        for iy in 0..field.ny {
            for ix in 0..field.nx {
                if ix + 1 < field.nx && inside(ix, iy, iz) != inside(ix + 1, iy, iz) {
                    crossed += 1;
                }
                if iy + 1 < field.ny && inside(ix, iy, iz) != inside(ix, iy + 1, iz) {
                    crossed += 1;
                }
                if iz + 1 < field.nz && inside(ix, iy, iz) != inside(ix, iy, iz + 1) {
                    crossed += 1;
                }
            }
        }
    }
    assert_eq!(vertices.len(), crossed);
}

#[test]
fn isovalue_shifts_the_extracted_sphere() {
    let field = sphere_field(31, 0.1, 1.0);
    let (vertices, tris) = {
        let triangles = IsosurfaceExtraction::new(0.3).extract(&field);
        weld_soup(&soup_of(&triangles))
    };
    assert!(is_closed_and_consistently_oriented(&tris));
    let (_, volume) = area_and_volume(&vertices, &tris);
    let exact = 4.0 / 3.0 * PI * 1.3_f64.powi(3);
    assert!(
        (volume / exact - 1.0).abs() < 0.02,
        "volume {volume} vs {exact}"
    );
}

#[test]
fn field_without_a_sign_change_gives_no_triangles() {
    let mut field = LevelSetField::new(6, 6, 6, 0.5, [0.0; 3]);
    for v in field.data.iter_mut() {
        *v = 0.25;
    }
    assert!(IsosurfaceExtraction::new(0.0).extract(&field).is_empty());
    assert_eq!(IsosurfaceExtraction::new(0.5).count_triangles(&field), 0);
}
