// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for V-HACD voxel decomposition.

use oxiphysics_geometry::vhacd::{VHacdConfig, VHacdVoxel};

// ---------------------------------------------------------------------------
// Helper: axis-aligned box triangle soup (12 triangles, 6 faces × 2)
// ---------------------------------------------------------------------------

fn box_tris(min: [f64; 3], max: [f64; 3]) -> Vec<[f64; 9]> {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    vec![
        // -X face
        [x0, y0, z0, x0, y1, z0, x0, y1, z1],
        [x0, y0, z0, x0, y1, z1, x0, y0, z1],
        // +X face
        [x1, y0, z0, x1, y1, z1, x1, y1, z0],
        [x1, y0, z0, x1, y0, z1, x1, y1, z1],
        // -Y face
        [x0, y0, z0, x1, y0, z1, x1, y0, z0],
        [x0, y0, z0, x0, y0, z1, x1, y0, z1],
        // +Y face
        [x0, y1, z0, x1, y1, z0, x1, y1, z1],
        [x0, y1, z0, x1, y1, z1, x0, y1, z1],
        // -Z face
        [x0, y0, z0, x1, y0, z0, x1, y1, z0],
        [x0, y0, z0, x1, y1, z0, x0, y1, z0],
        // +Z face
        [x0, y0, z1, x1, y1, z1, x1, y0, z1],
        [x0, y0, z1, x0, y1, z1, x1, y1, z1],
    ]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// A single convex box should decompose into exactly one part.
/// The hull volume reported by the part should be within 10 % of the analytic
/// volume (1.0 for a unit box).
#[test]
fn test_vhacd_convex_box_one_part() {
    let tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let cfg = VHacdConfig {
        resolution: 32,
        ..Default::default()
    };
    let decomposer = VHacdVoxel::new(cfg);
    let result = decomposer
        .decompose(&tris)
        .expect("decompose should succeed");

    assert_eq!(
        result.parts.len(),
        1,
        "single convex box should produce 1 part, got {}",
        result.parts.len()
    );

    let total_vol: f64 = result.parts.iter().map(|p| p.volume).sum();
    assert!(
        (total_vol - 1.0).abs() / 1.0 <= 0.10,
        "volume should be within 10% of 1.0, got {total_vol}"
    );
}

/// An L-shaped mesh (two boxes joined at a corner) is not convex.
/// The decomposer must split it into at least two parts that are each close to
/// convex, and the total hull volume must be within 10 % of the analytic
/// volume (1.5).
#[test]
fn test_vhacd_l_shape_two_parts() {
    // Box 1: [0,0,0] to [1,1,1]  (volume = 1.0)
    // Box 2: [1,0,0] to [2,1,0.5] (volume = 0.5)
    // Together they form an L-shape that is not convex.
    let mut tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    tris.extend(box_tris([1.0, 0.0, 0.0], [2.0, 1.0, 0.5]));

    let cfg = VHacdConfig {
        resolution: 32,
        concavity_threshold: 0.03,
        min_voxels_per_cluster: 8,
        ..Default::default()
    };
    let decomposer = VHacdVoxel::new(cfg.clone());
    let result = decomposer
        .decompose(&tris)
        .expect("decompose should succeed");

    assert!(
        result.parts.len() >= 2,
        "L-shape must be split into at least 2 parts, got {}",
        result.parts.len()
    );

    // Voxel hulls cover the voxelised mesh and exceed it by up to one voxel per
    // side, so the total converges to the analytic volume from above.
    let analytic_vol = 1.5_f64;
    let total_vol: f64 = result.parts.iter().map(|p| p.volume).sum();
    assert!(
        total_vol >= analytic_vol * 0.99 && total_vol <= analytic_vol * 1.10,
        "total volume should be within 10% of {analytic_vol}, got {total_vol}"
    );

    for (i, part) in result.parts.iter().enumerate() {
        assert!(
            part.concavity <= cfg.concavity_threshold + 1e-12,
            "part {i} is more concave ({}) than the threshold",
            part.concavity
        );
        assert!(part.volume > 0.0, "part {i} has no volume");
        assert!(!part.indices.is_empty(), "part {i} has no hull triangles");
    }
}

/// A convex solid that is not aligned with the voxel grid must stay in one
/// part: a cube rotated about all three axes.  The hull follows the solid (many
/// more than the 8 vertices of a box) and covers its volume.
#[test]
fn test_vhacd_rotated_convex_solid_is_one_part() {
    let rotate = |p: [f64; 3]| {
        let angles = [0.5_f64, 0.3, 0.7];
        let (sx, cx) = angles[0].sin_cos();
        let (sy, cy) = angles[1].sin_cos();
        let (sz, cz) = angles[2].sin_cos();
        let (y1, z1) = (cx * p[1] - sx * p[2], sx * p[1] + cx * p[2]);
        let (x2, z2) = (cy * p[0] + sy * z1, -sy * p[0] + cy * z1);
        [cz * x2 - sz * y1, sz * x2 + cz * y1, z2]
    };
    let tris: Vec<[f64; 9]> = box_tris([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5])
        .iter()
        .map(|t| {
            let a = rotate([t[0], t[1], t[2]]);
            let b = rotate([t[3], t[4], t[5]]);
            let c = rotate([t[6], t[7], t[8]]);
            [a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]]
        })
        .collect();
    let cfg = VHacdConfig {
        resolution: 32,
        ..Default::default()
    };
    let result = VHacdVoxel::new(cfg)
        .decompose(&tris)
        .expect("decompose should succeed");

    assert_eq!(
        result.parts.len(),
        1,
        "a convex solid is one part, got {}",
        result.parts.len()
    );
    let part = &result.parts[0];
    assert!(
        part.vertices.len() > 8,
        "{} hull vertices",
        part.vertices.len()
    );
    assert!(
        part.volume >= 1.0 && part.volume <= 1.8,
        "hull volume should cover the unit cube, got {}",
        part.volume
    );
}

/// Empty input must return an error (not panic).
#[test]
fn test_vhacd_empty_input_error() {
    let decomposer = VHacdVoxel::new(VHacdConfig::default());
    let result = decomposer.decompose(&[]);
    assert!(result.is_err(), "empty input should return an error");
}

/// All parts must have non-negative volume.
#[test]
fn test_vhacd_parts_non_negative_volume() {
    let tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let cfg = VHacdConfig {
        resolution: 16,
        ..Default::default()
    };
    let decomposer = VHacdVoxel::new(cfg);
    let result = decomposer
        .decompose(&tris)
        .expect("decompose should succeed");

    for (i, part) in result.parts.iter().enumerate() {
        assert!(
            part.volume >= 0.0,
            "part {} has negative volume: {}",
            i,
            part.volume
        );
    }
}

/// Parts count must not exceed the configured `max_parts`.
#[test]
fn test_vhacd_respects_max_parts() {
    let mut tris = box_tris([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    // Add more boxes to stress the splitter
    tris.extend(box_tris([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]));
    tris.extend(box_tris([0.0, 2.0, 0.0], [1.0, 3.0, 1.0]));

    let max_parts = 4;
    let cfg = VHacdConfig {
        resolution: 16,
        max_parts,
        concavity_threshold: 0.001,
        min_voxels_per_cluster: 4,
        ..Default::default()
    };
    let decomposer = VHacdVoxel::new(cfg);
    let result = decomposer
        .decompose(&tris)
        .expect("decompose should succeed");

    assert!(
        result.parts.len() <= max_parts,
        "should not exceed max_parts={max_parts}, got {}",
        result.parts.len()
    );
}
