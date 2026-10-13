// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the watershed region partitioning of the Recast
//! pipeline (`recast::region::build_regions`).

use oxiphysics_collision::recast::compact::{CompactHeightfield, build_compact_heightfield};
use oxiphysics_collision::recast::polymesh::poly_mesh_to_triangles;
use oxiphysics_collision::recast::rasterize::Heightfield;
use oxiphysics_collision::recast::region::{RegionHeightfield, build_regions};
use oxiphysics_collision::recast::{RecastBuilder, RecastConfig};
use std::collections::BTreeSet;
use std::time::Instant;

fn config() -> RecastConfig {
    RecastConfig {
        cell_size: 1.0,
        cell_height: 0.1,
        agent_height: 1.0,
        agent_max_climb: 0.2,
        region_min_size: 2,
        ..Default::default()
    }
}

/// Single-layer floor on a `width x depth` grid; `walkable(x, z)` selects the
/// cells that carry a floor span.
fn floor(
    width: usize,
    depth: usize,
    walkable: impl Fn(usize, usize) -> bool,
) -> CompactHeightfield {
    let mut hf = Heightfield::new(
        width,
        depth,
        [0.0; 3],
        [width as f64, 1.0, depth as f64],
        1.0,
        0.1,
    )
    .expect("heightfield");
    for z in 0..depth {
        for x in 0..width {
            if walkable(x, z) {
                hf.add_span(x, z, 0, 1, true);
            }
        }
    }
    build_compact_heightfield(&hf, &config())
}

/// Distinct region ids of the floor spans in the cells selected by `cells`.
fn regions_in(regions: &RegionHeightfield, cells: impl Fn(usize, usize) -> bool) -> BTreeSet<u16> {
    let chf = &regions.compact;
    let mut found = BTreeSet::new();
    for z in 0..chf.height {
        for x in 0..chf.width {
            if cells(x, z) {
                let cell = &chf.cells[z * chf.width + x];
                found.insert(chf.spans[cell.index as usize].reg);
            }
        }
    }
    found
}

#[test]
fn a_512_by_512_floor_is_partitioned_in_well_under_a_second() {
    // The flood fill used to look up the cell of every span by scanning the
    // whole grid, which made this size take minutes.
    let chf = floor(512, 512, |_, _| true);
    assert_eq!(chf.spans.len(), 512 * 512);
    let started = Instant::now();
    let regions = build_regions(&chf, &RecastConfig::default()).expect("regions");
    let elapsed = started.elapsed();
    assert!(
        elapsed.as_secs_f64() < 1.0,
        "build_regions took {elapsed:?} on a 512 x 512 floor"
    );
    let ids: BTreeSet<u16> = regions.compact.spans.iter().map(|s| s.reg).collect();
    assert_eq!(ids, BTreeSet::from([1]), "an open floor is one basin");
}

#[test]
fn rooms_joined_by_a_doorway_are_separate_regions() {
    let left = |x: usize, z: usize| (1..=12).contains(&x) && (1..=12).contains(&z);
    let right = |x: usize, z: usize| (16..=27).contains(&x) && (1..=12).contains(&z);
    let door = |x: usize, z: usize| (13..=15).contains(&x) && z == 6;
    let chf = floor(29, 14, |x, z| left(x, z) || right(x, z) || door(x, z));
    let regions = build_regions(&chf, &config()).expect("regions");
    let in_left = regions_in(&regions, left);
    let in_right = regions_in(&regions, right);
    assert_eq!(in_left.len(), 1, "left room: {in_left:?}");
    assert_eq!(in_right.len(), 1, "right room: {in_right:?}");
    assert!(in_left.is_disjoint(&in_right));
    assert!(!in_left.contains(&0) && !in_right.contains(&0));
}

#[test]
fn watershed_regions_are_numbered_from_one_and_no_floor_is_left_out() {
    // Pillars cut a 40 x 40 floor into basins.
    let pillar = |x: usize, z: usize| x % 8 == 4 && z % 8 == 4;
    let chf = floor(40, 40, |x, z| !pillar(x, z));
    let regions = build_regions(&chf, &config()).expect("regions");
    assert!(regions.compact.spans.iter().all(|s| s.reg != 0));
    let ids: BTreeSet<u16> = regions.compact.spans.iter().map(|s| s.reg).collect();
    assert_eq!(ids.first().copied(), Some(1));
    assert_eq!(ids.last().copied(), u16::try_from(ids.len()).ok());
}

#[test]
fn region_min_size_drops_floors_too_small_to_become_a_region() {
    let main_floor = |x: usize, z: usize| (1..=8).contains(&x) && (1..=8).contains(&z);
    let speck = |x: usize, z: usize| (12..=13).contains(&x) && (4..=5).contains(&z);
    let chf = floor(16, 10, |x, z| main_floor(x, z) || speck(x, z));
    // 4 x 4 = 16 spans needed per region.
    let cfg = RecastConfig {
        region_min_size: 4,
        ..config()
    };
    let regions = build_regions(&chf, &cfg).expect("regions");
    assert_eq!(regions_in(&regions, speck), BTreeSet::from([0]));
    assert_eq!(regions_in(&regions, main_floor), BTreeSet::from([1]));
    // A smaller threshold keeps the 2 x 2 speck.
    let cfg = RecastConfig {
        region_min_size: 2,
        ..config()
    };
    let regions = build_regions(&chf, &cfg).expect("regions");
    assert!(!regions_in(&regions, speck).contains(&0));
}

#[test]
fn more_than_65535_regions_is_an_error_instead_of_an_overflow() {
    // Walkable cells on a checkerboard never touch, so each one is a region.
    let placed = std::cell::Cell::new(0usize);
    let chf = floor(512, 256, |x, z| {
        if (x + z).is_multiple_of(2) {
            placed.set(placed.get() + 1);
            true
        } else {
            false
        }
    });
    assert_eq!(placed.get(), 65_536);
    let cfg = RecastConfig {
        region_min_size: 0,
        ..config()
    };
    match build_regions(&chf, &cfg) {
        Err(message) => assert!(message.contains("region id overflow"), "{message}"),
        Ok(_) => panic!("65536 regions cannot be numbered with 16-bit ids"),
    }
}

/// Rectangle `[x0, x1] x [z0, z1]` at height 0 as two upward triangles.
fn rect(x0: f64, x1: f64, z0: f64, z1: f64) -> [[f64; 9]; 2] {
    [
        [x0, 0.0, z0, x1, 0.0, z0, x1, 0.0, z1],
        [x0, 0.0, z0, x1, 0.0, z1, x0, 0.0, z1],
    ]
}

/// Whether the point `(x, z)` lies inside (or on the border of) one of the
/// triangles `[x0, y0, z0, x1, y1, z1, x2, y2, z2]`, seen from above.
fn covered(triangles: &[[f64; 9]], x: f64, z: f64) -> bool {
    triangles.iter().any(|t| {
        let side = |ax: f64, az: f64, bx: f64, bz: f64| (bx - ax) * (z - az) - (bz - az) * (x - ax);
        let signs = [
            side(t[0], t[2], t[3], t[5]),
            side(t[3], t[5], t[6], t[8]),
            side(t[6], t[8], t[0], t[2]),
        ];
        signs.iter().all(|&v| v >= -1e-9) || signs.iter().all(|&v| v <= 1e-9)
    })
}

#[test]
fn builder_meshes_two_rooms_and_the_door_between_them() {
    let mut tris = Vec::new();
    tris.extend(rect(0.0, 6.0, 0.0, 6.0));
    tris.extend(rect(8.0, 14.0, 0.0, 6.0));
    tris.extend(rect(6.0, 8.0, 2.5, 3.5));
    let cfg = RecastConfig {
        cell_size: 0.25,
        cell_height: 0.1,
        agent_height: 1.8,
        agent_radius: 0.0,
        agent_max_climb: 0.4,
        region_min_size: 2,
        ..Default::default()
    };
    let mesh = RecastBuilder::new(cfg)
        .build(&tris)
        .expect("two rooms build");
    let triangles = poly_mesh_to_triangles(&mesh);
    assert!(!triangles.is_empty());
    assert!(covered(&triangles, 3.0, 3.0), "left room is meshed");
    assert!(covered(&triangles, 11.0, 3.0), "right room is meshed");
    assert!(
        covered(&triangles, 7.0, 3.0),
        "the door between them is meshed"
    );
    assert!(!covered(&triangles, 7.0, 1.0), "no floor beside the door");
    assert!(!covered(&triangles, 7.0, 5.0), "no floor beside the door");
    // Two 36 m^2 rooms plus a 2 m^2 door, minus the one-cell ledge ring.
    let area: f64 = triangles
        .iter()
        .map(|t| 0.5 * ((t[5] - t[2]) * (t[6] - t[0]) - (t[3] - t[0]) * (t[8] - t[2])))
        .sum();
    assert!(area > 50.0 && area < 74.0, "mesh area {area}");
}
