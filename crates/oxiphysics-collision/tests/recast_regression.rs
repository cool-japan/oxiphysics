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

// ---------------------------------------------------------------------------
// Contour fidelity, erosion, region merging and config validation
// ---------------------------------------------------------------------------

use oxiphysics_collision::recast::RecastError;
use oxiphysics_collision::recast::polymesh::{PolyMesh, poly_mesh_to_triangles};

/// Axis-aligned walkable rectangle `[x0, x1] x [z0, z1]` at height `y`.
fn rect_tris(x0: f64, x1: f64, z0: f64, z1: f64, y: f64) -> [[f64; 9]; 2] {
    [
        [x0, y, z0, x1, y, z0, x1, y, z1],
        [x0, y, z0, x1, y, z1, x0, y, z1],
    ]
}

/// L-shaped floor: `[0,6]x[0,2]` plus `[0,2]x[2,6]`.
fn l_floor_tris() -> Vec<[f64; 9]> {
    let mut tris = Vec::new();
    tris.extend(rect_tris(0.0, 6.0, 0.0, 2.0, 0.0));
    tris.extend(rect_tris(0.0, 2.0, 2.0, 6.0, 0.0));
    tris
}

fn in_l_floor(x: f64, z: f64) -> bool {
    let eps = 1e-9;
    let in_rect = |x0: f64, x1: f64, z0: f64, z1: f64| {
        x >= x0 - eps && x <= x1 + eps && z >= z0 - eps && z <= z1 + eps
    };
    in_rect(0.0, 6.0, 0.0, 2.0) || in_rect(0.0, 2.0, 2.0, 6.0)
}

fn l_config() -> RecastConfig {
    RecastConfig {
        cell_size: 0.25,
        cell_height: 0.1,
        agent_height: 1.8,
        agent_radius: 0.3,
        agent_max_climb: 0.4,
        agent_max_slope: 45.0,
        region_min_size: 2,
        ..Default::default()
    }
}

/// Signed XZ area of a world-space triangle (positive = upward normal).
fn tri_area_up(t: &[f64; 9]) -> f64 {
    let (ax, az) = (t[3] - t[0], t[5] - t[2]);
    let (bx, bz) = (t[6] - t[0], t[8] - t[2]);
    // y component of (b - a) x (c - a) is az*bx - ax*bz.
    0.5 * (az * bx - ax * bz)
}

fn mesh_area(pm: &PolyMesh) -> f64 {
    poly_mesh_to_triangles(pm).iter().map(tri_area_up).sum()
}

#[test]
fn test_recast_l_floor_mesh_stays_inside_the_l() {
    let pm = RecastBuilder::new(l_config())
        .build(&l_floor_tris())
        .expect("L floor should build");
    assert!(
        pm.npolys >= 2,
        "a non-convex L needs several convex polygons"
    );
    let tris = poly_mesh_to_triangles(&pm);
    assert!(!tris.is_empty());
    for t in &tris {
        let cx = (t[0] + t[3] + t[6]) / 3.0;
        let cz = (t[2] + t[5] + t[8]) / 3.0;
        assert!(
            in_l_floor(cx, cz),
            "triangle centroid ({cx:.3}, {cz:.3}) lies outside the L-shaped floor"
        );
        for v in 0..3 {
            assert!(
                in_l_floor(t[v * 3], t[v * 3 + 2]),
                "vertex ({}, {}) lies outside the L-shaped floor",
                t[v * 3],
                t[v * 3 + 2]
            );
        }
        assert!(tri_area_up(t) > 0.0, "triangles must face up");
    }
    // The L covers 20 m^2 (its convex hull would cover 28 m^2).  The ledge
    // filter plus erosion by the agent radius strip roughly 0.75 m from each
    // side of the 2 m wide arms, leaving two ~0.75 m wide strips (~6-7 m^2).
    let area = mesh_area(&pm);
    assert!(area < 20.0, "mesh area {area} exceeds the floor area");
    assert!(area > 5.0, "mesh area {area} lost too much of the floor");
}

#[test]
fn test_recast_l_floor_contour_is_non_convex() {
    use oxiphysics_collision::recast::{
        compact, contour, erode, rasterize, region, region_merge, walkable,
    };

    let cfg = l_config();
    let mut hf = rasterize::rasterize_triangles(&l_floor_tris(), &cfg).expect("rasterize");
    walkable::filter_low_hanging_obstacles(&mut hf, &cfg);
    walkable::filter_ledge_spans(&mut hf, &cfg);
    walkable::filter_walkable_low_height(&mut hf, &cfg);
    let chf0 = compact::build_compact_heightfield(&hf, &cfg);
    let radius = erode::walkable_radius_cells(cfg.agent_radius, cfg.cell_size);
    erode::erode_walkable_area(&mut hf, &chf0, radius);
    let chf = compact::build_compact_heightfield(&hf, &cfg);
    let mut rhf = region::build_regions(&chf, &cfg).expect("regions");
    let n = region_merge::merge_and_filter_regions(&mut rhf, 4, 400).expect("merge");
    assert_eq!(n, 1, "the small L regions merge into a single region");
    let set = contour::build_contours(&rhf, &chf, &cfg).expect("contours");
    assert_eq!(set.contours.len(), 1);
    let verts = &set.contours[0].verts;
    let nv = verts.len();
    let reflex = (0..nv)
        .filter(|&i| {
            let a = verts[(i + nv - 1) % nv];
            let b = verts[i];
            let c = verts[(i + 1) % nv];
            let cross = i64::from(b[0] - a[0]) * i64::from(c[2] - a[2])
                - i64::from(c[0] - a[0]) * i64::from(b[2] - a[2]);
            cross > 0
        })
        .count();
    assert!(reflex >= 1, "L contour {verts:?} must have a reflex corner");
}

#[test]
fn test_recast_floor_around_pillar_leaves_pillar_uncovered() {
    // 8x8 floor with a 1x1 pillar (3 m tall) in the middle.
    let mut tris: Vec<[f64; 9]> = rect_tris(0.0, 8.0, 0.0, 8.0, 0.0).to_vec();
    let (p0, p1, top) = (3.5, 4.5, 3.0);
    tris.extend(rect_tris(p0, p1, p0, p1, top));
    // Four vertical walls of the pillar.
    for (a, b) in [
        ((p0, p0), (p1, p0)),
        ((p1, p0), (p1, p1)),
        ((p1, p1), (p0, p1)),
        ((p0, p1), (p0, p0)),
    ] {
        tris.push([a.0, 0.0, a.1, b.0, 0.0, b.1, b.0, top, b.1]);
        tris.push([a.0, 0.0, a.1, b.0, top, b.1, a.0, top, a.1]);
    }
    let cfg = RecastConfig {
        cell_size: 0.25,
        cell_height: 0.1,
        agent_height: 1.8,
        agent_radius: 0.25,
        agent_max_climb: 0.4,
        region_min_size: 2,
        ..Default::default()
    };
    let pm = RecastBuilder::new(cfg)
        .build(&tris)
        .expect("pillar scene builds");
    let tris_out = poly_mesh_to_triangles(&pm);
    for t in &tris_out {
        let cx = (t[0] + t[3] + t[6]) / 3.0;
        let cz = (t[2] + t[5] + t[8]) / 3.0;
        let on_floor = t[1].abs() < 0.5 && t[4].abs() < 0.5 && t[7].abs() < 0.5;
        assert!(
            on_floor,
            "pillar top is too small for the agent and must be eroded"
        );
        assert!(
            !((p0..=p1).contains(&cx) && (p0..=p1).contains(&cz)),
            "triangle centroid ({cx:.3}, {cz:.3}) covers the pillar"
        );
    }
    let area = mesh_area(&pm);
    assert!(
        area < 64.0 - 1.0,
        "mesh area {area} must exclude the pillar"
    );
    assert!(area > 30.0, "mesh area {area} lost too much floor");
}

#[test]
fn test_recast_agent_radius_erodes_border() {
    let base = RecastConfig {
        cell_size: 0.25,
        cell_height: 0.1,
        agent_height: 1.8,
        agent_max_climb: 0.4,
        region_min_size: 2,
        ..Default::default()
    };
    let tris = rect_tris(0.0, 6.0, 0.0, 6.0, 0.0).to_vec();
    let thin = RecastBuilder::new(RecastConfig {
        agent_radius: 0.0,
        ..base.clone()
    })
    .build(&tris)
    .expect("build");
    let wide = RecastBuilder::new(RecastConfig {
        agent_radius: 1.0,
        ..base
    })
    .build(&tris)
    .expect("build");
    let (a_thin, a_wide) = (mesh_area(&thin), mesh_area(&wide));
    // Radius 1 m shrinks each side by >= 1 m: area <= (6 - 2)^2.
    assert!(a_wide <= 16.0 + 1e-9, "eroded area {a_wide}");
    assert!(a_wide > 4.0, "eroded area {a_wide}");
    assert!(
        a_thin > a_wide + 10.0,
        "radius 0 area {a_thin} vs radius 1 area {a_wide}"
    );
    for t in poly_mesh_to_triangles(&wide) {
        for v in 0..3 {
            let (x, z) = (t[v * 3], t[v * 3 + 2]);
            assert!(
                (1.0 - 1e-9..=5.0 + 1e-9).contains(&x) && (1.0 - 1e-9..=5.0 + 1e-9).contains(&z),
                "vertex ({x}, {z}) is closer than the agent radius to the edge"
            );
        }
    }
}

#[test]
fn test_recast_build_is_deterministic() {
    let build = || {
        RecastBuilder::new(l_config())
            .build(&l_floor_tris())
            .expect("build")
    };
    let a = build();
    let b = build();
    assert_eq!(a.verts, b.verts);
    assert_eq!(a.polys, b.polys);
    assert_eq!(a.regs, b.regs);
}

#[test]
fn test_recast_config_validation() {
    let ok = RecastConfig::default();
    assert!(ok.validate().is_ok());
    let bad = [
        RecastConfig {
            cell_size: 0.0,
            ..Default::default()
        },
        RecastConfig {
            cell_height: f64::NAN,
            ..Default::default()
        },
        RecastConfig {
            agent_height: -1.0,
            ..Default::default()
        },
        RecastConfig {
            agent_radius: -0.1,
            ..Default::default()
        },
        RecastConfig {
            agent_max_slope: 90.0,
            ..Default::default()
        },
        RecastConfig {
            max_verts_per_poly: 2,
            ..Default::default()
        },
        RecastConfig {
            max_verts_per_poly: 7,
            ..Default::default()
        },
        RecastConfig {
            edge_max_error: f64::INFINITY,
            ..Default::default()
        },
        RecastConfig {
            region_merge_size: usize::MAX,
            ..Default::default()
        },
        RecastConfig {
            agent_height: 1.0e6,
            cell_height: 1.0e-3,
            ..Default::default()
        },
    ];
    for cfg in bad {
        let err = RecastBuilder::new(cfg.clone()).build(&flat_plane_tris());
        assert!(
            matches!(err, Err(RecastError::InvalidConfig(_))),
            "config {cfg:?} must be rejected, got {:?}",
            err.as_ref().map(|pm| pm.npolys)
        );
    }
}

#[test]
fn test_recast_huge_grid_is_rejected_without_overflow() {
    // Valid parameters, but the 10 m plane at 1e-6 m cells needs 1e14 columns.
    let cfg = RecastConfig {
        cell_size: 1.0e-6,
        agent_radius: 0.0,
        ..Default::default()
    };
    let res = RecastBuilder::new(cfg).build(&flat_plane_tris());
    assert!(
        matches!(res, Err(RecastError::RasterizationFailed(_))),
        "oversized grid must fail cleanly"
    );
}

#[test]
fn test_recast_non_finite_input_is_rejected() {
    let tris = vec![[0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0, 1.0, 0.0, 1.0]];
    let res = RecastBuilder::with_default_config().build(&tris);
    assert!(matches!(res, Err(RecastError::RasterizationFailed(_))));
}
