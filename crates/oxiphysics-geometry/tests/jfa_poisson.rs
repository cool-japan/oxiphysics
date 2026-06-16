// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for JFA signed-distance field and screened-Poisson reconstruction.

use oxiphysics_geometry::signed_distance_field::JfaGrid;
use oxiphysics_geometry::{OrientedPoint, screened_poisson_reconstruct};

// ─────────────────────────────────────────────────────────────────────────────
// Test 1: JFA sphere distance accuracy
// ─────────────────────────────────────────────────────────────────────────────

/// Build a 41×41×41 JFA grid with seeds on a thin shell at radius 15 centered at (20,20,20).
/// For a sampled subset of non-seed interior voxels, check that JFA |distance|
/// is within 1.0 voxel of the brute-force nearest-seed distance.
#[test]
fn test_jfa_sphere_distance_accuracy() {
    let n = 41usize;
    let dims = [n, n, n];
    let center = (20usize, 20usize, 20usize);
    let radius = 15.0_f64;
    let shell_thickness = 1.5_f64;
    let voxel_size = 1.0_f64;

    // Collect seed voxels: thin shell at radius 15
    let mut seeds: Vec<(usize, usize, usize)> = Vec::new();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let dx = i as f64 - center.0 as f64;
                let dy = j as f64 - center.1 as f64;
                let dz = k as f64 - center.2 as f64;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if (dist - radius).abs() < shell_thickness {
                    seeds.push((i, j, k));
                }
            }
        }
    }

    let grid = JfaGrid::build(dims, voxel_size, &seeds);

    // Build a set for quick seed membership lookup
    let seed_set: std::collections::HashSet<(usize, usize, usize)> = seeds.iter().copied().collect();

    // Sample ~200 voxels stepping by 3 in each axis, check non-seed ones
    let mut checked = 0usize;
    let step = 3usize;
    'outer: for k in (0..n).step_by(step) {
        for j in (0..n).step_by(step) {
            for i in (0..n).step_by(step) {
                if seed_set.contains(&(i, j, k)) {
                    continue;
                }

                // Brute-force: minimum Euclidean distance to any seed voxel
                let mut brute_min = f64::MAX;
                for &(si, sj, sk) in &seeds {
                    let dx = i as f64 - si as f64;
                    let dy = j as f64 - sj as f64;
                    let dz = k as f64 - sk as f64;
                    let d = voxel_size * (dx * dx + dy * dy + dz * dz).sqrt();
                    if d < brute_min {
                        brute_min = d;
                    }
                }

                let jfa_abs = grid.get(i, j, k).abs();
                assert!(
                    (jfa_abs - brute_min).abs() <= voxel_size + 1e-9,
                    "voxel ({i},{j},{k}): JFA |dist|={jfa_abs:.4} vs brute={brute_min:.4}, diff > 1 voxel"
                );

                checked += 1;
                if checked >= 200 {
                    break 'outer;
                }
            }
        }
    }

    assert!(checked > 0, "No voxels were checked");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2: JFA sphere sign correctness
// ─────────────────────────────────────────────────────────────────────────────

/// Verify that voxels clearly inside the sphere (dist < 13) have negative JFA distance
/// and voxels clearly outside (dist > 17) have positive JFA distance.
#[test]
fn test_jfa_sphere_sign_correctness() {
    let n = 41usize;
    let dims = [n, n, n];
    let cx = 20usize;
    let cy = 20usize;
    let cz = 20usize;
    let radius = 15.0_f64;
    let shell_thickness = 1.5_f64;
    let voxel_size = 1.0_f64;

    // Seed the thin shell at radius 15
    let mut seeds: Vec<(usize, usize, usize)> = Vec::new();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let dx = i as f64 - cx as f64;
                let dy = j as f64 - cy as f64;
                let dz = k as f64 - cz as f64;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if (dist - radius).abs() < shell_thickness {
                    seeds.push((i, j, k));
                }
            }
        }
    }

    let grid = JfaGrid::build(dims, voxel_size, &seeds);

    // Check a 5×5×5 sample of interior voxels (dist_from_center < 13 → clearly inside)
    // and exterior voxels (dist_from_center > 17 → clearly outside).
    let mut inside_checked = 0usize;
    let mut outside_checked = 0usize;

    // Interior voxels: offset from center within ±8 voxels ensuring dist < 13
    for dz in [-8i64, -5, -2, 1, 4].iter() {
        for dy in [-8i64, -5, -2, 1, 4].iter() {
            for dx in [-8i64, -5, -2, 1, 4].iter() {
                let ix = (cx as i64 + dx) as usize;
                let iy = (cy as i64 + dy) as usize;
                let iz = (cz as i64 + dz) as usize;
                let dist_from_center = ((*dx * *dx + *dy * *dy + *dz * *dz) as f64).sqrt();
                if dist_from_center < 13.0 {
                    let d = grid.get(ix, iy, iz);
                    assert!(
                        d < 0.0,
                        "voxel ({ix},{iy},{iz}) clearly inside (dist_from_center={dist_from_center:.2}) \
                         should have negative JFA dist, got {d:.4}"
                    );
                    inside_checked += 1;
                }
            }
        }
    }

    // Exterior voxels: pick corners far from center (dist > 17)
    // Use offsets ensuring dist > 17 from center at (20,20,20)
    let exterior_offsets: [(i64, i64, i64); 8] = [
        (-19, -19, -19),
        (19, -19, -19),
        (-19, 19, -19),
        (19, 19, -19),
        (-19, -19, 19),
        (19, -19, 19),
        (-19, 19, 19),
        (19, 19, 19),
    ];
    for &(dx, dy, dz) in exterior_offsets.iter() {
        let ix = (cx as i64 + dx).clamp(0, n as i64 - 1) as usize;
        let iy = (cy as i64 + dy).clamp(0, n as i64 - 1) as usize;
        let iz = (cz as i64 + dz).clamp(0, n as i64 - 1) as usize;
        let dist_from_center = ((dx * dx + dy * dy + dz * dz) as f64).sqrt();
        if dist_from_center > 17.0 {
            let d = grid.get(ix, iy, iz);
            assert!(
                d > 0.0,
                "voxel ({ix},{iy},{iz}) clearly outside (dist_from_center={dist_from_center:.2}) \
                 should have positive JFA dist, got {d:.4}"
            );
            outside_checked += 1;
        }
    }

    assert!(inside_checked > 0, "No interior voxels were checked");
    assert!(outside_checked > 0, "No exterior voxels were checked");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3: Screened Poisson sphere vertex reconstruction
// ─────────────────────────────────────────────────────────────────────────────

/// Sample 400 oriented points on a unit sphere via Fibonacci/golden-ratio distribution,
/// reconstruct via screened Poisson, and check the average radius of resulting vertices
/// is within 20% of 1.0.
#[test]
fn test_screened_poisson_sphere_vertices() {
    let n = 400usize;
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());

    let mut pts: Vec<OrientedPoint> = Vec::with_capacity(n);
    for i in 0..n {
        let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
        let r = (1.0 - y * y).max(0.0).sqrt();
        let theta = i as f64 * golden_angle;
        let pos = [r * theta.cos(), y, r * theta.sin()];
        let len = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        let normal = if len > 1e-12 {
            [pos[0] / len, pos[1] / len, pos[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        };
        pts.push(OrientedPoint { position: pos, normal });
    }

    let result = screened_poisson_reconstruct(&pts, 24, 10.0, 0.2);
    assert!(result.is_ok(), "reconstruction should succeed: {:?}", result.as_ref().err().map(|e| e.to_string()));

    let surface = result.unwrap();
    assert!(!surface.vertices.is_empty(), "reconstructed surface must have vertices");

    let vertex_count = surface.vertices.len();
    let sum_radius: f64 = surface.vertices.iter().map(|v| {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }).sum();
    let mean_radius = sum_radius / vertex_count as f64;

    assert!(
        (mean_radius - 1.0).abs() < 0.20,
        "mean vertex radius {mean_radius:.4} should be within 20% of 1.0"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4: Laplacian FD on quadratic function
// ─────────────────────────────────────────────────────────────────────────────

/// Verify 7-point finite-difference Laplacian of f(x,y,z) = x²+y²+z² equals 6.0
/// at interior point (2,2,2) on a 5×5×5 grid with h=1.0.
#[test]
fn test_laplacian_fd_quadratic() {
    let dims = [5usize, 5, 5];
    let h = 1.0_f64;

    // f(x,y,z) = x² + y² + z² at integer coordinates
    let f = |ix: usize, iy: usize, iz: usize| -> f64 {
        let x = ix as f64;
        let y = iy as f64;
        let z = iz as f64;
        x * x + y * y + z * z
    };

    let flat = |ix: usize, iy: usize, iz: usize| -> usize {
        iz * dims[1] * dims[0] + iy * dims[0] + ix
    };

    let mut values = vec![0.0_f64; dims[0] * dims[1] * dims[2]];
    for iz in 0..dims[2] {
        for iy in 0..dims[1] {
            for ix in 0..dims[0] {
                values[flat(ix, iy, iz)] = f(ix, iy, iz);
            }
        }
    }

    // Interior point (2,2,2)
    let (px, py, pz) = (2usize, 2usize, 2usize);
    let lap = (values[flat(px + 1, py, pz)]
        + values[flat(px - 1, py, pz)]
        + values[flat(px, py + 1, pz)]
        + values[flat(px, py - 1, pz)]
        + values[flat(px, py, pz + 1)]
        + values[flat(px, py, pz - 1)]
        - 6.0 * values[flat(px, py, pz)])
        / (h * h);

    assert!(
        (lap - 6.0).abs() < 0.01,
        "Laplacian of x²+y²+z² at (2,2,2) should be 6.0, got {lap:.6}"
    );
}
