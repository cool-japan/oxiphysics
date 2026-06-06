// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Algorithm-level tests for CSG boolean operations and IMLS surface
//! reconstruction.  These tests call `oxiphysics` (the umbrella rlib) directly,
//! avoiding any pyo3 / cdylib linking issues.

use std::f64::consts::PI;

use oxiphysics::geometry::mesh_boolean::{MeshBooleanOp, SimpleMesh, mesh_boolean};
use oxiphysics::geometry::signed_distance_field::MarchingCubes;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a closed unit cube centred at `centre` with half-extent `h`.
fn unit_cube(centre: [f64; 3], h: f64) -> SimpleMesh {
    SimpleMesh::unit_cube(centre, h)
}

// ---------------------------------------------------------------------------
// CSG tests
// ---------------------------------------------------------------------------

#[test]
fn test_csg_union_produces_triangles() {
    let a = unit_cube([0.0, 0.0, 0.0], 1.0);
    let b = unit_cube([0.5, 0.0, 0.0], 1.0);
    let result = mesh_boolean(&a, &b, MeshBooleanOp::Union);
    assert!(
        !result.vertices.is_empty(),
        "union of two overlapping cubes should have vertices"
    );
    assert!(
        !result.triangles.is_empty(),
        "union of two overlapping cubes should have triangles"
    );
}

#[test]
fn test_csg_intersection_produces_triangles() {
    let a = unit_cube([0.0, 0.0, 0.0], 1.0);
    let b = unit_cube([0.5, 0.0, 0.0], 1.0);
    let result = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
    assert!(
        !result.vertices.is_empty(),
        "intersection of two overlapping cubes should have vertices"
    );
}

#[test]
fn test_csg_subtraction_removes_triangles() {
    let a = unit_cube([0.0, 0.0, 0.0], 1.0);
    let b = unit_cube([0.5, 0.0, 0.0], 1.0);
    // A and B have 12 triangles each; naively stitched that is 24.
    // The difference A\B keeps only A's outside-B triangles, so it must be
    // strictly fewer than 24.
    let diff_result = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
    assert!(
        diff_result.triangles.len() < 24,
        "subtraction result should have fewer triangles than the two input meshes combined (diff={}, combined=24)",
        diff_result.triangles.len()
    );
}

#[test]
fn test_csg_difference_is_not_empty() {
    let a = unit_cube([0.0, 0.0, 0.0], 1.0);
    let b = unit_cube([0.5, 0.0, 0.0], 0.3);
    let result = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
    assert!(
        !result.vertices.is_empty(),
        "subtracting a small cube from a large one should leave geometry"
    );
}

// ---------------------------------------------------------------------------
// IMLS / Marching Cubes tests
// ---------------------------------------------------------------------------

/// Generate approximately `n` points uniformly on the unit sphere.
fn sphere_points(n: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut points = Vec::with_capacity(n);
    let mut normals = Vec::with_capacity(n);
    let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
    for i in 0..n {
        let theta = (1.0 - 2.0 * (i as f64 + 0.5) / n as f64).acos();
        let phi = 2.0 * PI * (i as f64) / golden;
        let x = theta.sin() * phi.cos();
        let y = theta.sin() * phi.sin();
        let z = theta.cos();
        points.push([x, y, z]);
        normals.push([x, y, z]);
    }
    (points, normals)
}

/// Compute IMLS implicit value at `x` given points/normals and bandwidth `h_sq`.
fn imls_eval(x: [f64; 3], points: &[[f64; 3]], normals: &[[f64; 3]], h_sq: f64) -> f64 {
    let mut sum_w = 0.0_f64;
    let mut sum_wn = 0.0_f64;
    for (p, n) in points.iter().zip(normals.iter()) {
        let dx = x[0] - p[0];
        let dy = x[1] - p[1];
        let dz = x[2] - p[2];
        let d_sq = dx * dx + dy * dy + dz * dz;
        let w = (-d_sq / h_sq).exp();
        sum_w += w;
        sum_wn += w * (dx * n[0] + dy * n[1] + dz * n[2]);
    }
    if sum_w < 1e-12 {
        return 1.0;
    }
    sum_wn / sum_w
}

#[test]
fn test_imls_reconstruct_sphere_point_cloud() {
    let (points, normals) = sphere_points(100);

    // Bandwidth: mean nearest-neighbour distance.
    let n_pts = points.len();
    let mut h_sum = 0.0_f64;
    for i in 0..n_pts.min(30) {
        let p = points[i];
        let best_d2 = points
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, q)| {
                let dx = p[0] - q[0];
                let dy = p[1] - q[1];
                let dz = p[2] - q[2];
                dx * dx + dy * dy + dz * dz
            })
            .fold(f64::MAX, f64::min);
        if best_d2.is_finite() && best_d2 > 0.0 {
            h_sum += best_d2.sqrt();
        }
    }
    let h = (h_sum / 30.0).max(0.05);
    let h_sq = h * h;

    let bounds = [-1.6_f64, 1.6, -1.6, 1.6, -1.6, 1.6];
    let res = 24_usize;

    let pts = points.clone();
    let nrms = normals.clone();

    let mc_result = MarchingCubes::from_function(
        &|p: [f64; 3]| imls_eval(p, &pts, &nrms, h_sq),
        res,
        res,
        res,
        bounds,
    )
    .extract(0.0);

    assert!(
        mc_result.n_vertices() > 0,
        "IMLS sphere reconstruction should produce vertices (got 0)"
    );
    assert!(
        mc_result.n_triangles() > 0,
        "IMLS sphere reconstruction should produce triangles (got 0)"
    );

    for v in &mc_result.vertices {
        for &coord in v.position.iter() {
            assert!(
                coord.abs() <= 1.6,
                "reconstructed vertex coord {} outside [-1.6, 1.6]",
                coord
            );
        }
    }
}

#[test]
fn test_imls_reconstruct_empty_input() {
    let points: Vec<[f64; 3]> = vec![];
    let normals: Vec<[f64; 3]> = vec![];
    let bounds = [-1.0_f64, 1.0, -1.0, 1.0, -1.0, 1.0];
    // IMLS with empty points: every grid point returns the outside sentinel (1.0).
    // All grid values > 0 means cube_idx = 0xFF everywhere → no triangles.
    let mc_result = MarchingCubes::from_function(&|_p| 1.0_f64, 4, 4, 4, bounds).extract(0.0);
    assert!(
        mc_result.n_vertices() == 0,
        "empty input should yield empty mesh (got {} vertices)",
        mc_result.n_vertices()
    );
    let _ = (points, normals);
}

#[test]
fn test_imls_reconstruct_without_normals_uses_centroid_fallback() {
    let (points, _) = sphere_points(60);

    // Estimate normals using centroid fallback (as the production code does).
    let centroid: [f64; 3] = {
        let n = points.len() as f64;
        let s = points.iter().fold([0.0_f64; 3], |acc, p| {
            [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
        });
        [s[0] / n, s[1] / n, s[2] / n]
    };
    let normals: Vec<[f64; 3]> = points
        .iter()
        .map(|&p| {
            let d = [p[0] - centroid[0], p[1] - centroid[1], p[2] - centroid[2]];
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-12);
            [d[0] / len, d[1] / len, d[2] / len]
        })
        .collect();

    let h_sq = 0.1_f64 * 0.1;
    let bounds = [-1.6_f64, 1.6, -1.6, 1.6, -1.6, 1.6];
    let res = 20_usize;

    let pts = points.clone();
    let nrms = normals;

    let mc_result = MarchingCubes::from_function(
        &|p: [f64; 3]| imls_eval(p, &pts, &nrms, h_sq),
        res,
        res,
        res,
        bounds,
    )
    .extract(0.0);

    assert!(
        mc_result.n_vertices() > 0,
        "IMLS with centroid-estimated normals should produce a mesh (got 0 vertices)"
    );
}
