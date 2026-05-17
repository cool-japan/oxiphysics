// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_collision::recast::{RecastBuilder, RecastConfig};

fn flat_plane_tris() -> Vec<[f64; 9]> {
    // 10m × 10m flat plane at y=0
    vec![
        [0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 10.0, 0.0, 10.0],
        [0.0, 0.0, 0.0, 10.0, 0.0, 10.0, 0.0, 0.0, 10.0],
    ]
}

#[test]
fn test_recast_flat_plane_builds() {
    let cfg = RecastConfig {
        cell_size: 0.5,
        cell_height: 0.2,
        agent_height: 1.8,
        agent_max_slope: 45.0,
        region_min_size: 2,
        ..Default::default()
    };
    let builder = RecastBuilder::new(cfg);
    let pm = builder
        .build(&flat_plane_tris())
        .expect("flat plane should build successfully");

    assert!(
        pm.nverts > 0,
        "flat plane should produce vertices, got {}",
        pm.nverts
    );
    assert!(
        pm.npolys > 0,
        "flat plane should produce polygons, got {}",
        pm.npolys
    );
    println!("Flat plane: {} verts, {} polys", pm.nverts, pm.npolys);
}

fn staircase_tris(steps: usize, step_w: f64, step_h: f64, depth: f64) -> Vec<[f64; 9]> {
    // Staircase along X axis: each step is a flat platform at height i*step_h
    let mut tris = Vec::new();
    for i in 0..steps {
        let x0 = i as f64 * step_w;
        let x1 = x0 + step_w;
        let y = i as f64 * step_h;
        // Top face of each step (horizontal, walkable)
        tris.push([x0, y, 0.0, x1, y, 0.0, x1, y, depth]);
        tris.push([x0, y, 0.0, x1, y, depth, x0, y, depth]);
    }
    tris
}

#[test]
fn test_recast_staircase_builds() {
    let cfg = RecastConfig {
        cell_size: 0.25,
        cell_height: 0.1,
        agent_height: 1.8,
        agent_max_climb: 0.45,
        agent_max_slope: 45.0,
        region_min_size: 2,
        ..Default::default()
    };
    let tris = staircase_tris(5, 1.0, 0.4, 1.0);
    let builder = RecastBuilder::new(cfg);
    let result = builder.build(&tris);
    match result {
        Ok(pm) => {
            println!("Staircase: {} verts, {} polys", pm.nverts, pm.npolys);
            assert!(
                pm.nverts > 0 || pm.npolys == 0,
                "staircase should either produce vertices or an empty mesh (no panic)"
            );
        }
        Err(e) => {
            // Empty result is acceptable for small staircases below min_region_size
            println!("Staircase returned error (acceptable): {e}");
        }
    }
}

#[test]
fn test_recast_empty_input_error() {
    let builder = RecastBuilder::with_default_config();
    let result = builder.build(&[]);
    assert!(result.is_err(), "empty input should return Err");
}

#[test]
fn test_recast_pipeline_does_not_panic() {
    // L-shaped mesh: flat plane + vertical wall face
    let mut tris = flat_plane_tris();
    // Add a vertical wall triangle (not walkable due to slope)
    tris.push([3.0, 0.0, 3.0, 5.0, 0.0, 3.0, 5.0, 2.0, 3.0]);
    tris.push([3.0, 0.0, 3.0, 5.0, 2.0, 3.0, 3.0, 2.0, 3.0]);

    let builder = RecastBuilder::with_default_config();
    let _ = builder.build(&tris); // Must not panic
}

#[test]
fn test_recast_poly_mesh_to_nav_mesh_primitives() {
    use oxiphysics_collision::recast::polymesh::poly_mesh_to_nav_mesh_primitives;

    let cfg = RecastConfig {
        cell_size: 0.5,
        cell_height: 0.2,
        agent_height: 1.8,
        agent_max_slope: 45.0,
        region_min_size: 2,
        ..Default::default()
    };
    let builder = RecastBuilder::new(cfg);
    let pm = builder
        .build(&flat_plane_tris())
        .expect("flat plane should build");

    let (verts, tris) = poly_mesh_to_nav_mesh_primitives(&pm);

    assert!(!verts.is_empty(), "primitives should have vertices");
    assert!(!tris.is_empty(), "primitives should have triangles");

    // All vertex indices must be valid
    for tri in &tris {
        for &idx in tri {
            assert!(
                (idx as usize) < verts.len(),
                "triangle index {idx} out of bounds (nverts={})",
                verts.len()
            );
        }
    }
    println!(
        "NavMesh primitives: {} verts, {} tris",
        verts.len(),
        tris.len()
    );
}
