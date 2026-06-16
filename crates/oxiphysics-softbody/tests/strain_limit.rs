// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for SVD-based per-triangle strain limiting.

use oxiphysics_softbody::xpbd::strain_limit::{StrainLimitConfig, apply_strain_limiting};

// ─────────────────────────── helper ──────────────────────────────────────────

/// Compute max singular value of the 2×2 deformation gradient for one triangle.
fn max_stretch(positions: &[[f64; 3]], rest_positions: &[[f64; 3]], tri: [usize; 3]) -> f64 {
    let [i0, i1, i2] = tri;
    let q0 = rest_positions[i0];
    let q1 = rest_positions[i1];
    let q2 = rest_positions[i2];
    let e0r = [q1[0] - q0[0], q1[1] - q0[1], q1[2] - q0[2]];
    let e1r = [q2[0] - q0[0], q2[1] - q0[1], q2[2] - q0[2]];

    // Local frame.
    let n0 = norm3(e0r);
    if n0 < 1e-12 {
        return 0.0;
    }
    let t1 = [e0r[0] / n0, e0r[1] / n0, e0r[2] / n0];
    let proj = dot3(e1r, t1);
    let e1_orth = [
        e1r[0] - proj * t1[0],
        e1r[1] - proj * t1[1],
        e1r[2] - proj * t1[2],
    ];
    let n1 = norm3(e1_orth);
    if n1 < 1e-12 {
        return 0.0;
    }
    let t2 = [e1_orth[0] / n1, e1_orth[1] / n1, e1_orth[2] / n1];

    let dm = [
        [dot3(e0r, t1), dot3(e1r, t1)],
        [dot3(e0r, t2), dot3(e1r, t2)],
    ];
    let det = dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0];
    if det.abs() < 1e-12 {
        return 0.0;
    }
    let dm_inv = [
        [dm[1][1] / det, -dm[0][1] / det],
        [-dm[1][0] / det, dm[0][0] / det],
    ];

    let p0 = positions[i0];
    let p1 = positions[i1];
    let p2 = positions[i2];
    let f0 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let f1 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

    let ds = [[dot3(f0, t1), dot3(f1, t1)], [dot3(f0, t2), dot3(f1, t2)]];

    // F = Ds * Dm_inv
    let f_2d = mat2x2_mul(ds, dm_inv);

    // Largest singular value of F via sqrt of largest eigenvalue of F^T F.
    let ft = [[f_2d[0][0], f_2d[1][0]], [f_2d[0][1], f_2d[1][1]]];
    let a = mat2x2_mul(ft, f_2d);
    let aa = a[0][0];
    let b = a[0][1];
    let d = a[1][1];
    let disc = ((aa - d) * 0.5) * ((aa - d) * 0.5) + b * b;
    let lam1 = ((aa + d) * 0.5 + disc.sqrt()).max(0.0);
    lam1.sqrt()
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn mat2x2_mul(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

// ─────────────────────────── tests ───────────────────────────────────────────

/// Equilateral triangle at rest with side length 1 in the XZ plane.
fn equilateral_rest() -> ([[f64; 3]; 3], [[usize; 3]; 1]) {
    let p0 = [0.0_f64, 0.0, 0.0];
    let p1 = [1.0_f64, 0.0, 0.0];
    let p2 = [0.5_f64, 0.0, 3_f64.sqrt() * 0.5];
    ([p0, p1, p2], [[0, 1, 2]])
}

// Test 1: Pre-stretched triangle gets clamped to <= max_stretch.
#[test]
fn test_strain_limiting_clamps_stretch() {
    let (rest, tris) = equilateral_rest();
    let rest_pos: Vec<[f64; 3]> = rest.to_vec();
    // Scale positions by 2x to create 2x stretch.
    let mut positions: Vec<[f64; 3]> = rest_pos
        .iter()
        .map(|p| [p[0] * 2.0, p[1] * 2.0, p[2] * 2.0])
        .collect();
    let inv_masses = vec![1.0_f64, 1.0, 1.0];
    let config = StrainLimitConfig::default();

    // Run several iterations to converge.
    for _ in 0..10 {
        apply_strain_limiting(&mut positions, &rest_pos, &tris, &inv_masses, &config);
    }

    let ms = max_stretch(&positions, &rest_pos, [0, 1, 2]);
    assert!(
        ms <= config.max_stretch + 0.01,
        "max_stretch after limiting = {ms}, expected <= {}",
        config.max_stretch
    );
}

// Test 2: No-op when max_stretch = infinity (effectively disabled).
#[test]
fn test_strain_limiting_noop_when_disabled() {
    let (rest, tris) = equilateral_rest();
    let rest_pos: Vec<[f64; 3]> = rest.to_vec();
    let mut positions: Vec<[f64; 3]> = rest_pos
        .iter()
        .map(|p| [p[0] * 2.0, p[1] * 2.0, p[2] * 2.0])
        .collect();
    let positions_before = positions.clone();
    let inv_masses = vec![1.0_f64, 1.0, 1.0];

    let config = StrainLimitConfig {
        max_stretch: f64::INFINITY,
        min_stretch: 0.0,
        iterations: 1,
    };
    apply_strain_limiting(&mut positions, &rest_pos, &tris, &inv_masses, &config);

    for (i, (a, b)) in positions.iter().zip(positions_before.iter()).enumerate() {
        for k in 0..3 {
            assert!(
                (a[k] - b[k]).abs() < 1e-12,
                "Vertex {i} coord {k} changed when strain limiting should be disabled"
            );
        }
    }
}

// Test 3: Pinned vertex (inv_mass=0) does not move.
#[test]
fn test_strain_limiting_pinned_vertex_unchanged() {
    let (rest, tris) = equilateral_rest();
    let rest_pos: Vec<[f64; 3]> = rest.to_vec();
    // 3x stretch.
    let mut positions: Vec<[f64; 3]> = rest_pos
        .iter()
        .map(|p| [p[0] * 3.0, p[1] * 3.0, p[2] * 3.0])
        .collect();
    let pin_pos = positions[0];
    let inv_masses = vec![0.0_f64, 1.0, 1.0]; // i0 is pinned.
    let config = StrainLimitConfig {
        iterations: 5,
        ..Default::default()
    };

    apply_strain_limiting(&mut positions, &rest_pos, &tris, &inv_masses, &config);

    for k in 0..3 {
        assert!(
            (positions[0][k] - pin_pos[k]).abs() < 1e-12,
            "Pinned vertex 0 coord {k} changed: was {}, now {}",
            pin_pos[k],
            positions[0][k]
        );
    }
}

// Test 4: Triangle at rest -- no correction applied.
#[test]
fn test_strain_limiting_rest_triangle_noop() {
    let (rest, tris) = equilateral_rest();
    let rest_pos: Vec<[f64; 3]> = rest.to_vec();
    let mut positions = rest_pos.clone();
    let inv_masses = vec![1.0_f64, 1.0, 1.0];
    let config = StrainLimitConfig::default();

    apply_strain_limiting(&mut positions, &rest_pos, &tris, &inv_masses, &config);

    for (i, (a, b)) in positions.iter().zip(rest_pos.iter()).enumerate() {
        for k in 0..3 {
            assert!(
                (a[k] - b[k]).abs() < 1e-12,
                "Rest triangle: vertex {i} coord {k} changed (delta={})",
                (a[k] - b[k]).abs()
            );
        }
    }
}

// Test 5: 100x gravity hang test -- 5x5 grid cloth, gravity=980 m/s^2, 200 substeps.
#[test]
fn test_100x_gravity_cloth_hang_strain_limit() {
    const NX: usize = 5;
    const NY: usize = 5;
    const DX: f64 = 0.2;
    const DY: f64 = 0.2;
    const GRAVITY: f64 = 980.0; // 100x g
    const DT: f64 = 1.0 / 600.0;
    const N_SUBSTEPS: usize = 200;

    let n = NX * NY;
    // Flat cloth in XZ plane: x = i*DX, z = j*DY, y = 0.
    let mut rest_pos: Vec<[f64; 3]> = Vec::with_capacity(n);
    for j in 0..NY {
        for i in 0..NX {
            rest_pos.push([i as f64 * DX, 0.0, j as f64 * DY]);
        }
    }
    let mut positions = rest_pos.clone();
    let mut velocities: Vec<[f64; 3]> = vec![[0.0; 3]; n];

    // Triangulate: for each quad (i,j)->(i+1,j)->(i,j+1) and (i+1,j)->(i+1,j+1)->(i,j+1).
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for j in 0..NY - 1 {
        for i in 0..NX - 1 {
            let v00 = j * NX + i;
            let v10 = j * NX + i + 1;
            let v01 = (j + 1) * NX + i;
            let v11 = (j + 1) * NX + i + 1;
            triangles.push([v00, v10, v01]);
            triangles.push([v10, v11, v01]);
        }
    }

    // Pin top row (j=0).
    let mut inv_masses: Vec<f64> = vec![1.0; n];
    for w in inv_masses.iter_mut().take(NX) {
        *w = 0.0;
    }

    let config = StrainLimitConfig {
        max_stretch: 1.05,
        min_stretch: 0.95,
        iterations: 2,
    };

    for _ in 0..N_SUBSTEPS {
        // Symplectic Euler: update velocity with gravity, then position.
        for idx in 0..n {
            if inv_masses[idx] > 0.0 {
                velocities[idx][1] -= GRAVITY * DT;
                positions[idx][0] += velocities[idx][0] * DT;
                positions[idx][1] += velocities[idx][1] * DT;
                positions[idx][2] += velocities[idx][2] * DT;
            }
        }
        // Apply strain limiting.
        apply_strain_limiting(&mut positions, &rest_pos, &triangles, &inv_masses, &config);
        // Update velocities from position corrections (XPBD style).
        // (We don't have prev_positions here, so we just let the limiter constrain stretch.)
    }

    // Assert no NaN.
    for (idx, p) in positions.iter().enumerate() {
        for (k, &coord) in p.iter().enumerate() {
            assert!(
                coord.is_finite(),
                "NaN/Inf in positions[{idx}][{k}] = {}",
                coord
            );
        }
    }

    // Assert max stretch <= 1.05 for all triangles.
    let mut worst = 0.0_f64;
    for &tri in &triangles {
        let s = max_stretch(&positions, &rest_pos, tri);
        if s > worst {
            worst = s;
        }
    }
    assert!(
        worst <= config.max_stretch + 0.01,
        "Worst triangle stretch = {worst}, expected <= {}",
        config.max_stretch
    );
}
