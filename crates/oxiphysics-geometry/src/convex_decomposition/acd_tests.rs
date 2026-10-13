// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the mesh-driven approximate convex decomposition.
//!
//! `ConvexPart::is_convex` only checks that a point set lies inside its own
//! hull, which always holds, so these tests check what decomposition
//! actually promises: the right number of parts, hull volumes that add up to
//! the volume of the solid, concave regions left uncovered, solid regions
//! covered, and concavity below the threshold.

use super::*;
use crate::Error;
use crate::convex_hull::ConvexHull3D;

type Mesh = (Vec<Vec3>, Vec<[usize; 3]>);

/// Axis-aligned box `[lo, hi]` with outward counter-clockwise triangles.
fn box_mesh(lo: Vec3, hi: Vec3) -> Mesh {
    let vertices = vec![
        [lo[0], lo[1], lo[2]],
        [hi[0], lo[1], lo[2]],
        [hi[0], hi[1], lo[2]],
        [lo[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]],
        [hi[0], lo[1], hi[2]],
        [hi[0], hi[1], hi[2]],
        [lo[0], hi[1], hi[2]],
    ];
    let triangles = vec![
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    (vertices, triangles)
}

/// Prism over a counter-clockwise polygon (`outline`, triangulated by
/// `cap_triangles`, indices into `outline`) between `z = 0` and `z = height`.
fn prism_mesh(outline: &[[f64; 2]], cap_triangles: &[[usize; 3]], height: f64) -> Mesh {
    let n = outline.len();
    let mut vertices: Vec<Vec3> = outline.iter().map(|p| [p[0], p[1], 0.0]).collect();
    vertices.extend(outline.iter().map(|p| [p[0], p[1], height]));
    let mut triangles = Vec::new();
    for t in cap_triangles {
        triangles.push([t[0], t[2], t[1]]);
        triangles.push([n + t[0], n + t[1], n + t[2]]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        triangles.push([i, j, n + j]);
        triangles.push([i, n + j, n + i]);
    }
    (vertices, triangles)
}

/// L-shaped prism: `[0,2]×[0,1] ∪ [0,1]×[0,2]`, thickness 1, volume 3.
fn l_prism() -> Mesh {
    let outline = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    prism_mesh(&outline, &[[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]], 1.0)
}

/// U-shaped prism: `[0,3]×[0,2]` minus the notch `[1,2]×[1,2]`, thickness 1,
/// volume 5.
fn u_prism() -> Mesh {
    let outline = [
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 2.0],
        [2.0, 2.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let caps = [
        [0, 1, 4],
        [0, 4, 5],
        [1, 2, 3],
        [1, 3, 4],
        [0, 5, 6],
        [0, 6, 7],
    ];
    prism_mesh(&outline, &caps, 1.0)
}

fn in_l(x: f64, y: f64) -> bool {
    (0.0..=2.0).contains(&x) && (0.0..=1.0).contains(&y)
        || (0.0..=1.0).contains(&x) && (0.0..=2.0).contains(&y)
}

fn in_u(x: f64, y: f64) -> bool {
    (0.0..=3.0).contains(&x) && (0.0..=2.0).contains(&y) && !(x > 1.0 && x < 2.0 && y > 1.0)
}

fn rotate(mesh: &Mesh, axis: Vec3, angle: f64) -> Mesh {
    let axis = normalize(axis);
    let (sin, cos) = angle.sin_cos();
    let rotated = mesh
        .0
        .iter()
        .map(|&v| {
            // Rodrigues' rotation formula.
            let k_cross_v = cross(axis, v);
            let k_dot_v = dot(axis, v);
            add(
                add(scale(v, cos), scale(k_cross_v, sin)),
                scale(axis, k_dot_v * (1.0 - cos)),
            )
        })
        .collect();
    (rotated, mesh.1.clone())
}

fn flip_winding(mesh: &Mesh) -> Mesh {
    (
        mesh.0.clone(),
        mesh.1.iter().map(|t| [t[0], t[2], t[1]]).collect(),
    )
}

/// Whether `p` lies in the hull described by the part's outward triangles.
fn hull_contains(part: &ConvexPart, p: Vec3, tol: f64) -> bool {
    part.indices.iter().all(|t| {
        let a = part.vertices[t[0] as usize];
        let b = part.vertices[t[1] as usize];
        let c = part.vertices[t[2] as usize];
        let n = normalize(cross(sub(b, a), sub(c, a)));
        dot(n, sub(p, a)) <= tol
    })
}

fn covered(parts: &ConvexParts, p: Vec3) -> bool {
    parts.parts.iter().any(|part| hull_contains(part, p, 1e-9))
}

fn total_volume(parts: &ConvexParts) -> f64 {
    parts.parts.iter().map(|p| p.volume).sum()
}

/// Every part's recorded hull data must agree with an independent hull of
/// its vertices, and the hull triangles must be outward oriented.
fn assert_parts_are_consistent_hulls(parts: &ConvexParts) {
    for (i, part) in parts.parts.iter().enumerate() {
        let hull = ConvexHull3D::build(&part.vertices).expect("part spans a volume");
        assert!(
            (hull.compute_volume() - part.volume).abs() <= 1e-9 * part.volume,
            "part {i}: volume {} vs hull {}",
            part.volume,
            hull.compute_volume()
        );
        let com = hull.compute_center_of_mass();
        assert!(
            norm(sub(com, part.centroid)) <= 1e-9 * (1.0 + norm(com)),
            "part {i}: centroid {:?} vs hull {:?}",
            part.centroid,
            com
        );
        assert!(!part.indices.is_empty(), "part {i} has no hull triangles");
        assert!(
            part.is_convex(1e-9),
            "part {i}: vertices are not in convex position"
        );
        for &v in &part.vertices {
            assert!(
                hull_contains(part, v, 1e-9),
                "part {i}: vertex outside its hull"
            );
        }
        assert!(
            hull_contains(part, part.centroid, 1e-9),
            "part {i}: centroid outside"
        );
    }
}

/// A part that is an axis-aligned box: hull volume equals the bounding box
/// volume and the centroid is the box centre.
fn assert_part_is_axis_aligned_box(part: &ConvexPart) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in &part.vertices {
        for k in 0..3 {
            lo[k] = lo[k].min(v[k]);
            hi[k] = hi[k].max(v[k]);
        }
    }
    let box_volume = (0..3).map(|k| hi[k] - lo[k]).product::<f64>();
    assert!(
        (part.volume - box_volume).abs() < 1e-9,
        "{} vs {box_volume}",
        part.volume
    );
    for k in 0..3 {
        assert!((part.centroid[k] - 0.5 * (lo[k] + hi[k])).abs() < 1e-9);
    }
}

#[test]
fn test_cube_is_a_single_convex_part() {
    let (vertices, triangles) = box_mesh([-1.0; 3], [1.0; 3]);
    let parts = ApproxConvexDecomp::new()
        .decompose(&vertices, &triangles)
        .expect("cube decomposes");
    assert_eq!(parts.len(), 1);
    let part = &parts.parts[0];
    assert!((part.volume - 8.0).abs() < 1e-12);
    assert!(norm(part.centroid) < 1e-12);
    assert!(part.concavity <= 1e-9, "concavity = {}", part.concavity);
    assert_eq!(part.vertices.len(), 8);
    assert_eq!(part.indices.len(), 12);
    assert_parts_are_consistent_hulls(&parts);
}

#[test]
fn test_l_prism_splits_into_two_boxes_covering_the_solid() {
    let (vertices, triangles) = l_prism();
    let parts = ApproxConvexDecomp::new()
        .decompose(&vertices, &triangles)
        .expect("L prism decomposes");
    // The hull of the L is a 3.5 unit slab, so one part is not enough.
    assert_eq!(
        parts.len(),
        2,
        "parts: {:?}",
        parts.parts.iter().map(|p| p.volume).collect::<Vec<_>>()
    );
    assert_parts_are_consistent_hulls(&parts);
    // Two boxes tile the solid exactly: hull volumes add up to 3.
    assert!(
        (total_volume(&parts) - 3.0).abs() < 1e-9,
        "total = {}",
        total_volume(&parts)
    );
    for part in &parts.parts {
        assert!(part.concavity <= 0.01, "concavity = {}", part.concavity);
        assert!(part.is_convex(1e-9), "part is not convex");
        assert_part_is_axis_aligned_box(part);
    }
    // The notch is not covered by any part ...
    for p in [[1.5, 1.5, 0.5], [1.2, 1.8, 0.5], [1.8, 1.2, 0.2]] {
        assert!(!covered(&parts, p), "notch point {p:?} is covered");
    }
    // ... and every point of the solid is.
    let mut checked = 0;
    for i in 0..=20 {
        for j in 0..=20 {
            for k in 0..=2 {
                let (x, y, z) = (0.1 * i as f64, 0.1 * j as f64, 0.5 * k as f64);
                if in_l(x, y) {
                    checked += 1;
                    assert!(
                        covered(&parts, [x, y, z]),
                        "solid point {:?} uncovered",
                        [x, y, z]
                    );
                }
            }
        }
    }
    assert!(checked > 100);
}

#[test]
fn test_u_prism_needs_three_parts() {
    let (vertices, triangles) = u_prism();
    let parts = ApproxConvexDecomp::new()
        .decompose(&vertices, &triangles)
        .expect("U prism decomposes");
    // Splitting stops at three convex parts instead of running on to
    // `max_parts`: the planar faces created by the first cut are not mistaken
    // for concave regions.
    assert_eq!(
        parts.len(),
        3,
        "volumes: {:?}",
        parts.parts.iter().map(|p| p.volume).collect::<Vec<_>>()
    );
    assert_parts_are_consistent_hulls(&parts);
    assert!(
        (total_volume(&parts) - 5.0).abs() < 1e-9,
        "total = {}",
        total_volume(&parts)
    );
    for part in &parts.parts {
        assert!(part.concavity <= 0.01);
        assert_part_is_axis_aligned_box(part);
    }
    for p in [[1.5, 1.5, 0.5], [1.1, 1.9, 0.5], [1.9, 1.1, 0.9]] {
        assert!(!covered(&parts, p), "notch point {p:?} is covered");
    }
    for i in 0..=30 {
        for j in 0..=20 {
            let (x, y) = (0.1 * i as f64, 0.1 * j as f64);
            if in_u(x, y) {
                assert!(
                    covered(&parts, [x, y, 0.5]),
                    "solid point ({x}, {y}) uncovered"
                );
            }
        }
    }
}

#[test]
fn test_rotated_l_prism_still_splits_into_two_parts() {
    let mesh = rotate(&l_prism(), [1.0, 2.0, 3.0], 0.7);
    let parts = ApproxConvexDecomp::new()
        .decompose(&mesh.0, &mesh.1)
        .expect("rotated L prism decomposes");
    assert_eq!(parts.len(), 2);
    assert_parts_are_consistent_hulls(&parts);
    assert!(
        (total_volume(&parts) - 3.0).abs() < 1e-8,
        "total = {}",
        total_volume(&parts)
    );
    for part in &parts.parts {
        assert!(part.concavity <= 0.01, "concavity = {}", part.concavity);
    }
    // The notch centre (1.5, 1.5, 0.5), rotated with the mesh, stays uncovered.
    let rotated_notch = rotate(&(vec![[1.5, 1.5, 0.5]], Vec::new()), [1.0, 2.0, 3.0], 0.7).0[0];
    assert!(!covered(&parts, rotated_notch));
}

#[test]
fn test_result_does_not_depend_on_triangle_winding() {
    let mesh = flip_winding(&l_prism());
    let parts = ApproxConvexDecomp::new()
        .decompose(&mesh.0, &mesh.1)
        .expect("inward wound L prism decomposes");
    assert_eq!(parts.len(), 2);
    assert!((total_volume(&parts) - 3.0).abs() < 1e-9);
}

#[test]
fn test_max_parts_limits_the_split_and_reports_remaining_concavity() {
    let (vertices, triangles) = u_prism();
    let config = AcdConfig {
        max_parts: 2,
        ..AcdConfig::default()
    };
    let parts = ApproxConvexDecomp::with_config(config)
        .decompose(&vertices, &triangles)
        .expect("U prism decomposes");
    assert_eq!(parts.len(), 2);
    let worst = parts
        .parts
        .iter()
        .map(|p| p.concavity)
        .fold(0.0_f64, f64::max);
    assert!(
        (worst - 0.5).abs() < 1e-9,
        "the L-shaped remainder keeps its concavity: {worst}"
    );
    let single = ApproxConvexDecomp::with_config(AcdConfig {
        max_parts: 1,
        ..AcdConfig::default()
    })
    .decompose(&vertices, &triangles)
    .expect("one part");
    assert_eq!(single.len(), 1);
    // The single part is the whole hull of the U: 3 x 2 x 1.
    assert!((single.parts[0].volume - 6.0).abs() < 1e-9);
    assert!(
        (single.parts[0].concavity - 0.5).abs() < 1e-9,
        "{}",
        single.parts[0].concavity
    );
}

#[test]
fn test_loose_concavity_threshold_keeps_the_hull() {
    let (vertices, triangles) = l_prism();
    let config = AcdConfig {
        max_concavity: 1.0,
        ..AcdConfig::default()
    };
    let parts = ApproxConvexDecomp::with_config(config)
        .decompose(&vertices, &triangles)
        .expect("decomposes");
    assert_eq!(parts.len(), 1);
    // The farthest hull point from the solid is the middle of the hull
    // diagonal, half a unit from each notch wall.
    assert!(
        (parts.parts[0].concavity - 0.5).abs() < 1e-9,
        "concavity = {}",
        parts.parts[0].concavity
    );
    assert!((parts.parts[0].volume - 3.5).abs() < 1e-9);
}

#[test]
fn test_min_volume_fraction_drops_small_parts() {
    let (vertices, triangles) = l_prism();
    // Hull of the whole L has volume 3.5; the small box has volume 1 (< 40 %).
    let config = AcdConfig {
        min_volume_fraction: 0.4,
        ..AcdConfig::default()
    };
    let parts = ApproxConvexDecomp::with_config(config)
        .decompose(&vertices, &triangles)
        .expect("decomposes");
    assert_eq!(parts.len(), 1);
    assert!((parts.parts[0].volume - 2.0).abs() < 1e-9);
}

#[test]
fn test_invalid_meshes_are_rejected() {
    let acd = ApproxConvexDecomp::new();
    let (vertices, triangles) = box_mesh([0.0; 3], [1.0; 3]);
    assert!(matches!(
        acd.decompose(&[], &triangles),
        Err(Error::InvalidMesh { .. })
    ));
    assert!(matches!(
        acd.decompose(&vertices, &[]),
        Err(Error::InvalidMesh { .. })
    ));
    let mut bad_index = triangles.clone();
    bad_index[3] = [0, 1, 8];
    assert!(matches!(
        acd.decompose(&vertices, &bad_index),
        Err(Error::IndexOutOfBounds { index: 8, len: 8 })
    ));
    let mut nan = vertices.clone();
    nan[2][1] = f64::NAN;
    assert!(matches!(
        acd.decompose(&nan, &triangles),
        Err(Error::InvalidMesh { .. })
    ));
    // Every triangle has zero area.
    let collapsed = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert!(matches!(
        acd.decompose(&collapsed, &[[0, 1, 2]]),
        Err(Error::DegenerateGeometry { .. })
    ));
    // A single triangle spans no volume.
    let flat = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(matches!(
        acd.decompose(&flat, &[[0, 1, 2]]),
        Err(Error::DegenerateGeometry { .. })
    ));
}

#[test]
fn test_invalid_configuration_is_rejected() {
    let (vertices, triangles) = box_mesh([0.0; 3], [1.0; 3]);
    let run = |config: AcdConfig| {
        ApproxConvexDecomp::with_config(config).decompose(&vertices, &triangles)
    };
    let base = AcdConfig::default();
    assert!(matches!(
        run(AcdConfig {
            max_parts: 0,
            ..base.clone()
        }),
        Err(Error::OutOfRange {
            name: "max_parts",
            ..
        })
    ));
    assert!(matches!(
        run(AcdConfig {
            max_concavity: -0.1,
            ..base.clone()
        }),
        Err(Error::OutOfRange {
            name: "max_concavity",
            ..
        })
    ));
    assert!(matches!(
        run(AcdConfig {
            max_concavity: f64::NAN,
            ..base.clone()
        }),
        Err(Error::OutOfRange {
            name: "max_concavity",
            ..
        })
    ));
    assert!(matches!(
        run(AcdConfig {
            min_volume_fraction: 1.5,
            ..base.clone()
        }),
        Err(Error::OutOfRange {
            name: "min_volume_fraction",
            ..
        })
    ));
    assert!(matches!(
        run(AcdConfig {
            min_volume_fraction: f64::NAN,
            ..base
        }),
        Err(Error::OutOfRange {
            name: "min_volume_fraction",
            ..
        })
    ));
}

#[test]
fn test_unused_vertices_and_duplicate_triangles_are_harmless() {
    let (mut vertices, mut triangles) = box_mesh([0.0; 3], [2.0; 3]);
    vertices.push([50.0, 50.0, 50.0]); // referenced by no triangle
    let duplicate = triangles[0];
    triangles.push(duplicate);
    let parts = ApproxConvexDecomp::new()
        .decompose(&vertices, &triangles)
        .expect("cube decomposes");
    assert_eq!(parts.len(), 1);
    assert!((parts.parts[0].volume - 8.0).abs() < 1e-12);
}
