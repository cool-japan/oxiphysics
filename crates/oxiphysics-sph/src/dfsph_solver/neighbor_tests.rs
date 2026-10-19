// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Tests of the spatial hash and the neighbour lists of `DfsphSolver`: exactness
// at the 10 000-particle scale, correctness after repeated rebuilds with moving
// particles (the cell and list buffers are reused, so stale entries would show
// up here), and the cell-key hasher.

use std::hash::{BuildHasher, Hasher};

use super::*;

fn config(support_radius: f64) -> DfsphConfig {
    DfsphConfig {
        rest_density: 1000.0,
        smoothing_length: 0.5 * support_radius,
        support_radius,
        max_pressure_iterations: 10,
        max_divergence_iterations: 10,
        error_threshold: 1e-3,
        omega: 0.5,
        viscosity: 1e-6,
        gravity: [0.0; 3],
        surface_tension: 0.0,
    }
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Indices of all particles strictly closer than `support` to particle `i`
/// (including `i` itself), by brute force.
fn brute_force_neighbors(positions: &[[f64; 3]], i: usize, support: f64) -> Vec<usize> {
    (0..positions.len())
        .filter(|&j| norm_sq3(sub3(positions[i], positions[j])) < support * support)
        .collect()
}

/// A 22³ = 10 648 particle lattice (spacing 0.05, support 0.13 = 2.6 spacings):
/// every list equals the one obtained by enumerating the lattice offsets with
/// `|offset|² < 2.6² = 6.76`, interior particles have exactly 81 entries
/// (themselves included).
#[test]
fn test_neighbor_lists_of_a_ten_thousand_particle_lattice_are_exact() {
    let n = 22usize;
    let spacing = 0.05;
    let support = 0.13;
    let mut solver = DfsphSolver::new(config(support));
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                solver.add_particle(pos, [0.0; 3], 1.0);
            }
        }
    }
    assert_eq!(solver.num_particles(), 10_648);
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    let index = |ix: usize, iy: usize, iz: usize| (ix * n + iy) * n + iz;
    let reach = 3i64;
    let mut interior = 0;
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let mut expected = Vec::new();
                for dx in -reach..=reach {
                    for dy in -reach..=reach {
                        for dz in -reach..=reach {
                            if dx * dx + dy * dy + dz * dz >= 7 {
                                continue; // 6.76 < 7: offsets with |d|² ≤ 6 only
                            }
                            let (jx, jy, jz) = (ix as i64 + dx, iy as i64 + dy, iz as i64 + dz);
                            if (0..n as i64).contains(&jx)
                                && (0..n as i64).contains(&jy)
                                && (0..n as i64).contains(&jz)
                            {
                                expected.push(index(jx as usize, jy as usize, jz as usize));
                            }
                        }
                    }
                }
                expected.sort_unstable();
                let mut got = solver.neighbor_lists[index(ix, iy, iz)].clone();
                got.sort_unstable();
                assert_eq!(got, expected, "particle ({ix}, {iy}, {iz})");
                if (3..n - 3).contains(&ix) && (3..n - 3).contains(&iy) && (3..n - 3).contains(&iz)
                {
                    assert_eq!(got.len(), 81);
                    interior += 1;
                }
            }
        }
    }
    assert_eq!(interior, 16 * 16 * 16);
}

/// Rebuilding after the particles moved (some across cell boundaries, some far
/// away) must not leave stale entries in the reused cell and list buffers.
#[test]
fn test_rebuilding_after_motion_leaves_no_stale_entries() {
    let mut state = 0xC2B2_AE3D_27D4_EB4F_u64;
    let support = 0.2;
    let mut solver = DfsphSolver::new(config(support));
    for _ in 0..300 {
        let pos = [
            0.9 * lcg_unit(&mut state) - 0.45,
            0.9 * lcg_unit(&mut state) - 0.45,
            0.9 * lcg_unit(&mut state) - 0.45,
        ];
        solver.add_particle(pos, [0.0; 3], 0.01);
    }
    for round in 0..5 {
        solver.build_spatial_hash();
        solver.build_neighbor_lists();
        let positions: Vec<[f64; 3]> = solver.particles.iter().map(|p| p.position).collect();
        assert_eq!(solver.neighbor_lists.len(), positions.len());
        for i in 0..positions.len() {
            let mut got = solver.neighbor_lists[i].clone();
            got.sort_unstable();
            assert_eq!(
                got,
                brute_force_neighbors(&positions, i, support),
                "round {round}, particle {i}"
            );
        }
        // Every hash cell holds at least one particle (empty cells are dropped).
        assert!(solver.spatial_hash.values().all(|cell| !cell.is_empty()));
        let stored: usize = solver.spatial_hash.values().map(Vec::len).sum();
        assert_eq!(stored, positions.len(), "round {round}");
        // Move: small jitters, a few jumps over cell boundaries and one far away.
        for (k, p) in solver.particles.iter_mut().enumerate() {
            for d in 0..3 {
                p.position[d] += 0.15 * (lcg_unit(&mut state) - 0.5);
            }
            if k % 37 == 0 {
                p.position[0] += 1.7;
            }
        }
        solver.particles[0].position = [50.0 + round as f64, -30.0, 12.0];
    }
}

/// Shrinking the particle set also shrinks the neighbour table.
#[test]
fn test_neighbor_table_follows_the_particle_count() {
    let mut solver = DfsphSolver::new(config(0.2));
    for k in 0..20 {
        solver.add_particle([0.04 * k as f64, 0.0, 0.0], [0.0; 3], 0.01);
    }
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    assert_eq!(solver.neighbor_lists.len(), 20);
    solver.particles.truncate(5);
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    assert_eq!(solver.neighbor_lists.len(), 5);
    // Particles 0..5 are 0.04 apart, the farthest pair 0.16: all inside the
    // support radius 0.2 of each other.
    for (i, list) in solver.neighbor_lists.iter().enumerate() {
        let mut sorted = list.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4], "particle {i}");
    }
}

/// The inner vectors are cleared and refilled, not reallocated: with
/// unchanged positions the list allocations stay where they are.
#[test]
fn test_neighbor_lists_are_refilled_in_place() {
    let mut state = 0x1656_67B1_9E37_79F9_u64;
    let mut solver = DfsphSolver::new(config(0.2));
    for _ in 0..200 {
        let pos = [
            0.6 * lcg_unit(&mut state),
            0.6 * lcg_unit(&mut state),
            0.6 * lcg_unit(&mut state),
        ];
        solver.add_particle(pos, [0.0; 3], 0.01);
    }
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    let before: Vec<(usize, usize)> = solver
        .neighbor_lists
        .iter()
        .map(|list| (list.as_ptr() as usize, list.capacity()))
        .collect();
    let reference = solver.neighbor_lists.clone();
    for _ in 0..3 {
        solver.build_spatial_hash();
        solver.build_neighbor_lists();
    }
    let after: Vec<(usize, usize)> = solver
        .neighbor_lists
        .iter()
        .map(|list| (list.as_ptr() as usize, list.capacity()))
        .collect();
    assert_eq!(before, after, "the neighbour vectors were reallocated");
    assert_eq!(solver.neighbor_lists, reference);
}

/// The cell-key hasher is deterministic and separates the 27 cells around a
/// query cell (and a long row of cells) into different hashes.
#[test]
fn test_cell_key_hasher_is_deterministic_and_separates_neighbouring_cells() {
    let build = BuildHasherDefault::<CellKeyHasher>::default();
    let key = (12_i64, -7_i64, 31_i64);
    assert_eq!(build.hash_one(key), build.hash_one(key));
    let mut hashes = Vec::new();
    for dx in -1..=1 {
        for dy in -1..=1 {
            for dz in -1..=1 {
                hashes.push(build.hash_one((key.0 + dx, key.1 + dy, key.2 + dz)));
            }
        }
    }
    let mut unique = hashes.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(unique.len(), 27);
    let row: Vec<u64> = (0..1000_i64)
        .map(|z| build.hash_one((0_i64, 0_i64, z)))
        .collect();
    let mut unique_row = row.clone();
    unique_row.sort_unstable();
    unique_row.dedup();
    assert_eq!(unique_row.len(), 1000);
    // The generic byte path agrees with the integer path for 8-byte chunks.
    let mut by_bytes = CellKeyHasher::default();
    by_bytes.write(&5_i64.to_le_bytes());
    let mut by_integer = CellKeyHasher::default();
    by_integer.write_i64(5);
    assert_eq!(by_bytes.finish(), by_integer.finish());
}
