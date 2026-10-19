// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the mixed Voronoi areas and of the curvatures built on them.

use super::test_support::{discrete, jittered_grid, open_cylinder};
use super::*;
use crate::test_meshes::icosphere;

/// The Voronoi area formula with the cotangents paired to the wrong edges, as the
/// code computed it before the fix: the cotangent at `p_j` multiplies the squared
/// length of the edge `(p_i, p_j)`, which is the edge *opposite* `p_k`.
fn swapped_pairing_area(mesh: &DiscreteMesh, v: usize) -> f64 {
    let mut area = 0.0;
    for t in &mesh.triangles {
        let Some(i) = t.iter().position(|&w| w == v) else {
            continue;
        };
        let pi = mesh.vertices[t[i]];
        let pj = mesh.vertices[t[(i + 1) % 3]];
        let pk = mesh.vertices[t[(i + 2) % 3]];
        let cot_j = cotan_at(pi, pj, pk);
        let cot_k = cotan_at(pi, pk, pj);
        let ek_sq = dot3(sub3(pi, pj), sub3(pi, pj));
        let ej_sq = dot3(sub3(pi, pk), sub3(pi, pk));
        area += 0.125 * (cot_j * ek_sq + cot_k * ej_sq);
    }
    area
}

/// Area of the part of a planar acute triangle closer to corner `i` than to the
/// other two: the quadrilateral corner, edge midpoint, circumcentre, edge midpoint.
fn geometric_voronoi_region(p: [[f64; 3]; 3], i: usize) -> f64 {
    let (a, b, c) = (p[0], p[1], p[2]);
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    let (sa, sb, sc) = (
        a[0] * a[0] + a[1] * a[1],
        b[0] * b[0] + b[1] * b[1],
        c[0] * c[0] + c[1] * c[1],
    );
    let centre = [
        (sa * (b[1] - c[1]) + sb * (c[1] - a[1]) + sc * (a[1] - b[1])) / d,
        (sa * (c[0] - b[0]) + sb * (a[0] - c[0]) + sc * (b[0] - a[0])) / d,
    ];
    let corner = p[i];
    let mid_next = [
        0.5 * (corner[0] + p[(i + 1) % 3][0]),
        0.5 * (corner[1] + p[(i + 1) % 3][1]),
    ];
    let mid_prev = [
        0.5 * (corner[0] + p[(i + 2) % 3][0]),
        0.5 * (corner[1] + p[(i + 2) % 3][1]),
    ];
    let polygon = [[corner[0], corner[1]], mid_next, centre, mid_prev];
    let mut twice_area = 0.0;
    for k in 0..4 {
        let (u, w) = (polygon[k], polygon[(k + 1) % 4]);
        twice_area += u[0] * w[1] - w[0] * u[1];
    }
    0.5 * twice_area.abs()
}

#[test]
fn right_isosceles_triangle_has_the_exact_corner_areas() {
    // Legs of length one: the right-angle corner owns 1/4, each 45 degree corner 1/8.
    // The swapped pairing gave 1/4 to every corner, 3/4 in total.
    let m = DiscreteMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    );
    let areas = m.mixed_voronoi_areas();
    assert!((areas[0] - 0.25).abs() < 1e-14, "{areas:?}");
    assert!((areas[1] - 0.125).abs() < 1e-14, "{areas:?}");
    assert!((areas[2] - 0.125).abs() < 1e-14, "{areas:?}");
    assert!((areas.iter().sum::<f64>() - 0.5).abs() < 1e-14);
    for (v, &a) in areas.iter().enumerate() {
        assert!((m.voronoi_area(v) - a).abs() < 1e-15, "vertex {v}");
    }
    let old: f64 = (0..3).map(|v| swapped_pairing_area(&m, v)).sum();
    assert!(
        (old - 0.75).abs() < 1e-14,
        "the swapped pairing sums to {old}"
    );
}

#[test]
fn acute_triangle_corner_areas_are_the_geometric_voronoi_regions() {
    let p = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [1.0, 3.0, 0.0]];
    let m = DiscreteMesh::new(p.to_vec(), vec![[0, 1, 2]]);
    let areas = m.mixed_voronoi_areas();
    for (i, &area) in areas.iter().enumerate() {
        let expected = geometric_voronoi_region(p, i);
        assert!(
            (area - expected).abs() < 1e-12,
            "corner {i}: {area} vs {expected}"
        );
    }
    assert!((areas.iter().sum::<f64>() - 6.0).abs() < 1e-12);
    // This triangle is far from equilateral, so the swapped pairing is visibly wrong:
    // its corner areas add up to 7.08 instead of the triangle area 6.
    let wrong: f64 = (0..3).map(|v| swapped_pairing_area(&m, v)).sum();
    assert!((wrong - 6.0).abs() > 1.0, "swapped pairing sums to {wrong}");
}

#[test]
fn obtuse_triangle_gives_half_to_the_obtuse_corner_and_a_quarter_to_the_others() {
    // The angle at (1, 0.5, 0) is obtuse; the triangle has area 1.
    let m = DiscreteMesh::new(
        vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [1.0, 0.5, 0.0]],
        vec![[0, 1, 2]],
    );
    let areas = m.mixed_voronoi_areas();
    assert!((areas[0] - 0.25).abs() < 1e-14, "{areas:?}");
    assert!((areas[1] - 0.25).abs() < 1e-14, "{areas:?}");
    assert!((areas[2] - 0.5).abs() < 1e-14, "{areas:?}");
    // The obtuse corner also wins for every cyclic rotation of the vertex order.
    let rotated = DiscreteMesh::new(m.vertices.clone(), vec![[2, 0, 1]]);
    let areas = rotated.mixed_voronoi_areas();
    assert!((areas[2] - 0.5).abs() < 1e-14, "{areas:?}");
}

#[test]
fn vertex_areas_add_up_to_the_surface_area_of_any_mesh() {
    let meshes = [
        discrete(icosphere(3)),
        discrete(jittered_grid(12, 0.3, 7)),
        discrete(open_cylinder(16, 6, 2.0, 3.0)),
        unit_tetrahedron(),
        unit_icosahedron(),
    ];
    for (k, m) in meshes.iter().enumerate() {
        let total = total_surface_area(m);
        let sum: f64 = m.mixed_voronoi_areas().iter().sum();
        assert!(
            (sum - total).abs() < 1e-12 * total,
            "mesh {k}: vertex areas {sum} vs surface area {total}"
        );
    }
    // The jittered grid has obtuse triangles; the swapped pairing does not partition it.
    let m = discrete(jittered_grid(12, 0.3, 7));
    let old: f64 = (0..m.num_vertices())
        .map(|v| swapped_pairing_area(&m, v))
        .sum();
    assert!((old - total_surface_area(&m)).abs() > 1e-3 * total_surface_area(&m));
}

#[test]
fn single_vertex_query_matches_the_bulk_computation() {
    let m = discrete(jittered_grid(10, 0.3, 11));
    let bulk = m.mixed_voronoi_areas();
    for (v, &a) in bulk.iter().enumerate() {
        assert!(
            (m.voronoi_area(v) - a).abs() < 1e-14,
            "vertex {v}: {} vs {a}",
            m.voronoi_area(v)
        );
    }
}

#[test]
fn vertices_without_triangles_get_the_area_floor() {
    let m = DiscreteMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [5.0, 5.0, 5.0],
        ],
        vec![[0, 1, 2]],
    );
    assert_eq!(m.mixed_voronoi_areas()[3], MIN_VORONOI_AREA);
    assert_eq!(m.voronoi_area(3), MIN_VORONOI_AREA);
    assert_eq!(m.voronoi_area(99), MIN_VORONOI_AREA);
}

#[test]
fn degenerate_triangles_carry_no_area() {
    let m = DiscreteMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 1, 3]],
    );
    let areas = m.mixed_voronoi_areas();
    assert!(areas.iter().all(|a| a.is_finite() && *a > 0.0));
    assert!((areas.iter().sum::<f64>() - 0.5).abs() < 1e-12);
}

#[test]
fn angle_defects_add_up_to_four_pi_on_a_sphere() {
    // Gauss-Bonnet: sum K_v A_v = 2 pi chi = 4 pi, exactly, on every sphere mesh.
    for m in [
        unit_tetrahedron(),
        unit_icosahedron(),
        discrete(icosphere(3)),
    ] {
        let k = m.gaussian_curvature();
        let a = m.mixed_voronoi_areas();
        let total: f64 = k.iter().zip(a.iter()).map(|(k, a)| k * a).sum();
        assert!(
            (total - 4.0 * PI).abs() < 1e-11,
            "integrated curvature {total}"
        );
    }
}

#[test]
fn sphere_curvatures_converge_to_one_with_the_mixed_area() {
    let m = discrete(icosphere(3));
    let k = m.gaussian_curvature();
    let h = m.mean_curvature();
    let mean_k = k.iter().sum::<f64>() / k.len() as f64;
    assert!(
        (mean_k - 1.0).abs() < 0.02,
        "mean Gaussian curvature {mean_k}"
    );
    for (v, (&kv, &hv)) in k.iter().zip(h.iter()).enumerate() {
        assert!((kv - 1.0).abs() < 0.08, "K[{v}] = {kv}");
        assert!((hv - 1.0).abs() < 0.05, "H[{v}] = {hv}");
    }
    let principal = m.principal_curvatures();
    for (v, &(k1, k2)) in principal.iter().enumerate() {
        assert!(k1 >= k2);
        assert!(
            (k1 - 1.0).abs() < 0.35 && (k2 - 1.0).abs() < 0.35,
            "{v}: {k1} {k2}"
        );
    }
}

#[test]
fn cylinder_has_mean_curvature_half_over_radius_and_no_gaussian_curvature() {
    // A cylinder of radius two: H = 1 / (2 R) = 0.25 and K = 0 away from the rims.
    let (nu, nv) = (48, 10);
    let m = discrete(open_cylinder(nu, nv, 2.0, 5.0));
    let h = m.mean_curvature();
    let k = m.gaussian_curvature();
    for j in 1..nv {
        for i in 0..nu {
            let v = j * nu + i;
            assert!((h[v] - 0.25).abs() < 0.01, "H[{v}] = {}", h[v]);
            assert!(k[v].abs() < 1e-9, "K[{v}] = {}", k[v]);
        }
    }
}

#[test]
fn cotan_laplacian_is_reproducible_and_symmetric() {
    let m = discrete(jittered_grid(8, 0.3, 3));
    let first = m.cotan_laplacian();
    let second = m.cotan_laplacian();
    assert_eq!(first.len(), second.len());
    for (a, b) in first.iter().zip(second.iter()) {
        assert_eq!((a.0, a.1), (b.0, b.1));
        assert_eq!(a.2.to_bits(), b.2.to_bits());
    }
    let weight = |i: usize, j: usize| -> f64 {
        first
            .iter()
            .filter(|e| e.0 == i && e.1 == j)
            .map(|e| e.2)
            .sum()
    };
    for &(i, j, w) in first.iter().take(200) {
        assert!((weight(j, i) - w).abs() < 1e-14);
    }
}
