// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`heat_method_geodesic`].
//!
//! The old implementation added only one of the two Crane divergence terms of every
//! edge, so the right-hand side of the Poisson problem did not sum to zero, and
//! solved the singular cotangent Laplacian with 500 Jacobi sweeps and no
//! convergence check.  The distances were wrong by tens of percent.  The tests
//! compare against exact Euclidean and great-circle distances.

use super::functions::heat_method_geodesic;
use super::functions_2::build_grid_mesh;
use super::types::{GeodesicMesh, HeatMethodParams};
use crate::Error;
use crate::test_meshes::icosphere;
use std::f64::consts::PI;

fn euclid(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Largest and root-mean-square deviation of `dist` from the Euclidean distance to
/// the nearest of `sources`, both relative to the diameter of the mesh.
fn deviation_from_euclid(mesh: &GeodesicMesh, dist: &[f64], sources: &[usize]) -> (f64, f64) {
    let mut diameter = 0.0_f64;
    for a in &mesh.vertices {
        for b in &mesh.vertices {
            diameter = diameter.max(euclid(*a, *b));
        }
    }
    let (mut worst, mut sum_sq) = (0.0_f64, 0.0_f64);
    for (v, &d) in dist.iter().enumerate() {
        let exact = sources
            .iter()
            .map(|&s| euclid(mesh.vertices[v], mesh.vertices[s]))
            .fold(f64::INFINITY, f64::min);
        let err = (d - exact).abs();
        worst = worst.max(err);
        sum_sq += err * err;
    }
    (
        worst / diameter,
        (sum_sq / dist.len() as f64).sqrt() / diameter,
    )
}

#[test]
fn flat_grid_30x30_corner_source_is_within_five_percent_of_the_diameter() {
    let mesh = build_grid_mesh(30, 30, 1.0);
    let dist = heat_method_geodesic(&mesh, &[0], &HeatMethodParams::default()).expect("valid mesh");
    assert_eq!(dist.len(), mesh.num_vertices());
    let (worst, rms) = deviation_from_euclid(&mesh, &dist, &[0]);
    assert!(worst < 0.05, "worst error {worst} of the diameter");
    assert!(rms < 0.02, "rms error {rms} of the diameter");
    // The source is at distance exactly zero and nothing is nearer to it.
    assert_eq!(dist[0], 0.0);
    assert!(dist.iter().all(|&d| d >= 0.0));
    // Along the diagonal the distance grows monotonically.
    let mut previous = -1.0;
    for i in 0..30 {
        let d = dist[i * 30 + i];
        assert!(d > previous, "diagonal distance not increasing at {i}");
        previous = d;
    }
}

#[test]
fn flat_grid_central_source_is_within_five_percent_of_the_diameter() {
    let mesh = build_grid_mesh(30, 30, 1.0);
    let source = 15 * 30 + 15;
    let dist =
        heat_method_geodesic(&mesh, &[source], &HeatMethodParams::default()).expect("valid mesh");
    let (worst, rms) = deviation_from_euclid(&mesh, &dist, &[source]);
    assert!(worst < 0.05, "worst error {worst} of the diameter");
    assert!(rms < 0.02, "rms error {rms} of the diameter");
}

#[test]
fn distances_do_not_depend_on_the_scale_of_the_mesh() {
    let unit = build_grid_mesh(12, 12, 1.0);
    let big = build_grid_mesh(12, 12, 250.0);
    let params = HeatMethodParams::default();
    let d_unit = heat_method_geodesic(&unit, &[5], &params).expect("valid mesh");
    let d_big = heat_method_geodesic(&big, &[5], &params).expect("valid mesh");
    for (a, b) in d_unit.iter().zip(d_big.iter()) {
        assert!((250.0 * a - b).abs() <= 1e-6 * b.max(1.0), "{a} vs {b}");
    }
}

#[test]
fn several_sources_give_the_distance_to_the_nearest() {
    let mesh = build_grid_mesh(30, 30, 1.0);
    let sources = [0, 29, 29 * 30];
    let dist =
        heat_method_geodesic(&mesh, &sources, &HeatMethodParams::default()).expect("valid mesh");
    let (worst, _) = deviation_from_euclid(&mesh, &dist, &sources);
    assert!(worst < 0.05, "worst error {worst} of the diameter");
    // The sources themselves are the nearest points.
    for &s in &sources {
        assert!(dist[s] < 0.8, "source {s} at distance {}", dist[s]);
    }
}

#[test]
fn the_time_factor_changes_the_smoothing_but_not_the_geometry() {
    let mesh = build_grid_mesh(24, 24, 1.0);
    for factor in [0.5, 1.0, 4.0] {
        let params = HeatMethodParams {
            time_factor: factor,
        };
        let dist = heat_method_geodesic(&mesh, &[0], &params).expect("valid mesh");
        let (worst, _) = deviation_from_euclid(&mesh, &dist, &[0]);
        assert!(worst < 0.08, "time factor {factor}: worst error {worst}");
    }
}

#[test]
fn strip_distances_are_monotone() {
    // build_grid_mesh(5, 1, ..) has no faces; a strip needs two rows.
    let mesh = build_grid_mesh(12, 2, 1.0);
    let dist = heat_method_geodesic(&mesh, &[0], &HeatMethodParams::default()).expect("valid mesh");
    for i in 1..12 {
        assert!(
            dist[i] > dist[i - 1],
            "bottom row not increasing: d[{i}] = {} <= d[{}] = {}",
            dist[i],
            i - 1,
            dist[i - 1]
        );
        assert!(
            (dist[i] - i as f64).abs() < 0.4,
            "d[{i}] = {} should be close to {i}",
            dist[i]
        );
    }
}

#[test]
fn icosphere_distances_follow_the_great_circle_distance() {
    let (vertices, faces) = icosphere(3);
    let mesh = GeodesicMesh::new(vertices, faces);
    let source = 0;
    let dist =
        heat_method_geodesic(&mesh, &[source], &HeatMethodParams::default()).expect("valid mesh");
    let mut worst = 0.0_f64;
    for (v, &d) in dist.iter().enumerate() {
        let a = mesh.vertices[source];
        let b = mesh.vertices[v];
        let cos = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).clamp(-1.0, 1.0);
        worst = worst.max((d - cos.acos()).abs());
    }
    // The sphere is the unit sphere: the farthest point is at distance pi.
    assert!(worst < 0.06 * PI, "worst error {worst}");
    let farthest = dist.iter().copied().fold(0.0_f64, f64::max);
    assert!((farthest - PI).abs() < 0.06 * PI, "farthest {farthest}");
}

#[test]
fn a_long_strip_beyond_the_double_precision_range_of_the_heat_still_works() {
    // The heat decays by a factor ~2.4 per edge: after ~650 edges it underflows.
    // The diffusion time is enlarged until it reaches every vertex.
    let mesh = build_grid_mesh(800, 2, 1.0);
    let dist = heat_method_geodesic(&mesh, &[0], &HeatMethodParams::default()).expect("valid mesh");
    assert_eq!(dist.len(), 1600);
    let mut previous = f64::NEG_INFINITY;
    for (i, &d) in dist.iter().take(800).enumerate() {
        assert!(d.is_finite());
        assert!(d >= previous - 1e-6, "not monotone at {i}");
        previous = d;
    }
    assert!((dist[799] - 799.0).abs() < 0.05 * 799.0, "{}", dist[799]);
}

#[test]
fn vertices_in_other_components_are_infinitely_far() {
    // Two separate 4x4 grids in one mesh.
    let a = build_grid_mesh(4, 4, 1.0);
    let b = build_grid_mesh(4, 4, 1.0);
    let offset = a.num_vertices();
    let mut vertices = a.vertices.clone();
    vertices.extend(b.vertices.iter().map(|p| [p[0] + 100.0, p[1], p[2]]));
    let mut faces = a.faces.clone();
    faces.extend(
        b.faces
            .iter()
            .map(|f| [f[0] + offset, f[1] + offset, f[2] + offset]),
    );
    let mesh = GeodesicMesh::new(vertices, faces);
    let dist = heat_method_geodesic(&mesh, &[0], &HeatMethodParams::default()).expect("valid mesh");
    assert!(dist[..offset].iter().all(|d| d.is_finite()));
    assert!(dist[offset..].iter().all(|d| *d == f64::INFINITY));
    // A source in each component reaches everything.
    let both = heat_method_geodesic(&mesh, &[0, offset], &HeatMethodParams::default())
        .expect("valid mesh");
    assert!(both.iter().all(|d| d.is_finite()));
}

#[test]
fn without_sources_nothing_is_reachable() {
    let mesh = build_grid_mesh(5, 5, 1.0);
    let dist = heat_method_geodesic(&mesh, &[], &HeatMethodParams::default()).expect("valid mesh");
    assert_eq!(dist.len(), 25);
    assert!(dist.iter().all(|d| *d == f64::INFINITY));
    let empty = GeodesicMesh::new(Vec::new(), Vec::new());
    let none = heat_method_geodesic(&empty, &[], &HeatMethodParams::default()).expect("valid mesh");
    assert!(none.is_empty());
}

#[test]
fn invalid_input_is_reported_instead_of_panicking() {
    let mesh = build_grid_mesh(4, 4, 1.0);
    for factor in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let params = HeatMethodParams {
            time_factor: factor,
        };
        let result = heat_method_geodesic(&mesh, &[0], &params);
        assert!(
            matches!(
                result,
                Err(Error::OutOfRange {
                    name: "time_factor",
                    ..
                })
            ),
            "time factor {factor}: {result:?}"
        );
    }
    let params = HeatMethodParams::default();
    let result = heat_method_geodesic(&mesh, &[99], &params);
    assert!(
        matches!(result, Err(Error::IndexOutOfBounds { index: 99, len: 16 })),
        "{result:?}"
    );
    let mut broken = mesh.clone();
    broken.faces.push([0, 1, 40]);
    let result = heat_method_geodesic(&broken, &[0], &params);
    assert!(
        matches!(result, Err(Error::InvalidMesh { .. })),
        "{result:?}"
    );
    let mut nan = mesh.clone();
    nan.vertices[3][1] = f64::NAN;
    let result = heat_method_geodesic(&nan, &[0], &params);
    assert!(
        matches!(result, Err(Error::InvalidMesh { .. })),
        "{result:?}"
    );
}
