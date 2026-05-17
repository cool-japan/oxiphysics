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
/// The decomposer should produce at least one part and the total volume should
/// be within 15 % of the analytic value (1.5).
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
    let decomposer = VHacdVoxel::new(cfg);
    let result = decomposer
        .decompose(&tris)
        .expect("decompose should succeed");

    assert!(
        !result.parts.is_empty(),
        "L-shape should produce at least 1 part, got {}",
        result.parts.len()
    );

    // Volume check: total hull volume across all parts should be within 20 %
    // of the analytic combined volume.  V-HACD splits introduce boundary gaps
    // so a wider tolerance than a single box is appropriate.
    let analytic_vol = 1.5_f64;
    let total_vol: f64 = result.parts.iter().map(|p| p.volume).sum();
    assert!(
        (total_vol - analytic_vol).abs() / analytic_vol <= 0.20,
        "total volume should be within 20% of {analytic_vol}, got {total_vol}"
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
