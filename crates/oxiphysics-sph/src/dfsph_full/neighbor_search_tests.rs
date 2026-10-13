// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Tests of the hashed neighbour search used by the full DFSPH solver
// (`find_neighbors`, `DfsphSimulation`, `DfsphFactorTable`). They live in their
// own file so that `dfsph_full.rs` stays below the 2000-line limit.

use super::*;

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// `n` pseudo-random points in the cube `[-half, half]³` (so negative cell
/// coordinates are exercised as well).
fn random_cloud(n: usize, half: f64, seed: u64) -> Vec<[f64; 3]> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            [
                half * (2.0 * lcg_unit(&mut state) - 1.0),
                half * (2.0 * lcg_unit(&mut state) - 1.0),
                half * (2.0 * lcg_unit(&mut state) - 1.0),
            ]
        })
        .collect()
}

/// Cubic lattice of `n_side³` points with spacing `dx`, every node displaced by
/// up to `±0.3 dx` so that no pair sits exactly on a kernel-support boundary.
fn jittered_lattice(n_side: usize, dx: f64, seed: u64) -> Vec<[f64; 3]> {
    let mut state = seed;
    let mut points = Vec::with_capacity(n_side * n_side * n_side);
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                points.push([
                    (ix as f64 + 0.6 * (lcg_unit(&mut state) - 0.5)) * dx,
                    (iy as f64 + 0.6 * (lcg_unit(&mut state) - 0.5)) * dx,
                    (iz as f64 + 0.6 * (lcg_unit(&mut state) - 0.5)) * dx,
                ]);
            }
        }
    }
    points
}

#[test]
fn test_find_neighbors_equals_brute_force_on_random_cloud() {
    let positions = random_cloud(300, 0.5, 0x9E37_79B9_7F4A_7C15);
    let h = 0.05;
    let hashed = find_neighbors(&positions, h);
    let brute = find_neighbors_brute(&positions, h);
    assert_eq!(hashed, brute);
    let total: usize = hashed.iter().map(Vec::len).sum();
    assert!(
        total > 300,
        "the cloud must be dense enough to have neighbours, got {total} pair entries"
    );
}

#[test]
fn test_find_neighbors_equals_brute_force_on_large_cloud() {
    let positions = random_cloud(1500, 0.4, 0xD1B5_4A32_D192_ED03);
    let h = 0.03;
    assert_eq!(
        find_neighbors(&positions, h),
        find_neighbors_brute(&positions, h)
    );
}

/// A regular lattice has many pairs exactly `2h` apart. The hashed search keeps
/// them (its distance test is inclusive) while the brute-force reference drops
/// them; those are the only differences, and the cubic kernel vanishes there.
#[test]
fn test_find_neighbors_lattice_adds_only_pairs_at_the_support_radius() {
    let spacing = 0.05;
    let h = 0.1;
    let support = 2.0 * h;
    let mut positions = Vec::new();
    for ix in 0..6 {
        for iy in 0..6 {
            for iz in 0..6 {
                positions.push([
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ]);
            }
        }
    }
    let hashed = find_neighbors(&positions, h);
    let brute = find_neighbors_brute(&positions, h);
    let mut boundary_pairs = 0usize;
    for i in 0..positions.len() {
        for j in &brute[i] {
            assert!(
                hashed[i].contains(j),
                "pair ({i}, {j}) inside the support is missing"
            );
        }
        for &j in &hashed[i] {
            if brute[i].contains(&j) {
                continue;
            }
            let r = norm3(sub3(positions[i], positions[j]));
            assert!(
                (r - support).abs() <= 1e-12 * support,
                "extra pair ({i}, {j}) at r = {r} is not on the support radius {support}"
            );
            assert_eq!(cubic_w(r, h), 0.0, "kernel must vanish at the support");
            assert_eq!(
                cubic_grad_w(r, h),
                0.0,
                "gradient must vanish at the support"
            );
            boundary_pairs += 1;
        }
    }
    assert!(
        boundary_pairs > 0,
        "the lattice contains pairs exactly 2h apart"
    );
}

#[test]
fn test_find_neighbors_are_symmetric_sorted_and_exclude_self() {
    let positions = random_cloud(400, 0.3, 0x2545_F491_4F6C_DD1D);
    let lists = find_neighbors(&positions, 0.04);
    assert_eq!(lists.len(), positions.len());
    for (i, list) in lists.iter().enumerate() {
        assert!(
            list.windows(2).all(|w| w[0] < w[1]),
            "list of particle {i} is not strictly ascending: {list:?}"
        );
        assert!(!list.contains(&i), "particle {i} lists itself");
        for &j in list {
            assert!(
                lists[j].contains(&i),
                "neighbour relation not symmetric for ({i}, {j})"
            );
        }
    }
}

#[test]
fn test_find_neighbors_empty_input_and_empty_support() {
    assert!(find_neighbors(&[], 0.1).is_empty());
    let positions = random_cloud(20, 0.1, 3);
    for bad_h in [0.0, -0.1, f64::NAN] {
        let lists = find_neighbors(&positions, bad_h);
        assert_eq!(lists.len(), positions.len());
        assert!(lists.iter().all(Vec::is_empty), "h = {bad_h}");
        assert_eq!(lists, find_neighbors_brute(&positions, bad_h));
    }
}

/// A diverged particle (infinite, NaN or astronomically large coordinates) must
/// neither panic the search nor acquire neighbours, and must not disturb the
/// lists of the well-behaved particles.
#[test]
fn test_find_neighbors_tolerates_non_finite_and_escaped_particles() {
    let mut positions = random_cloud(60, 0.15, 0xA076_1D64_78BD_642F);
    let healthy = positions.len();
    positions.push([f64::INFINITY, 0.0, 0.0]);
    positions.push([0.0, f64::NEG_INFINITY, 0.0]);
    positions.push([f64::NAN, 0.0, 0.0]);
    positions.push([1.0e300, -1.0e300, 1.0e300]);
    positions.push([-1.0e12, 5.0e11, 0.0]);
    let h = 0.05;
    let lists = find_neighbors(&positions, h);
    assert_eq!(lists.len(), positions.len());
    for (i, list) in lists.iter().enumerate().skip(healthy) {
        assert!(list.is_empty(), "escaped particle {i} has neighbours");
    }
    let reference = find_neighbors_brute(&positions[..healthy], h);
    assert_eq!(&lists[..healthy], &reference[..]);
}

/// Far from the origin the cell index saturates; the distance test must still
/// give the exact brute-force neighbour sets.
#[test]
fn test_find_neighbors_far_from_origin_still_matches_brute_force() {
    let offset = 1.0e9;
    let positions: Vec<[f64; 3]> = random_cloud(50, 0.02, 11)
        .into_iter()
        .map(|p| [p[0] + offset, p[1] - offset, p[2] + offset])
        .collect();
    let h = 0.01;
    assert_eq!(
        find_neighbors(&positions, h),
        find_neighbors_brute(&positions, h)
    );
}

/// With both pressure solves disabled the particles move ballistically
/// (`x += v dt`), so pairs cross the kernel-support radius as the cloud
/// expands. After every `advance()` the stored neighbour lists must be exactly
/// the brute-force lists of the new positions.
#[test]
fn test_simulation_neighbors_track_positions_across_steps() {
    let dx = 0.05;
    let h = 0.05; // support 2h = 0.1 = two lattice spacings
    let rho0 = 1000.0;
    let mass = rho0 * dx * dx * dx;
    let dt = 1.0e-3;
    let initial = jittered_lattice(4, dx, 0x5851_F42D_4C95_7F2D);
    let n = initial.len() as f64;
    let centre = {
        let sum = initial.iter().fold([0.0; 3], |acc, p| add3(acc, *p));
        scale3(sum, 1.0 / n)
    };
    let mut state = DfsphState::new(h, rho0, dt);
    for position in &initial {
        state.add_particle(*position, mass);
    }
    // Homogeneous expansion v = 2 (x - c) [1/s].
    for (velocity, position) in state.velocities.iter_mut().zip(initial.iter()) {
        *velocity = scale3(sub3(*position, centre), 2.0);
    }
    let mut sim = DfsphSimulation::new(state, [0.0; 3], 0.0);
    sim.max_div_iter = 0;
    sim.max_density_iter = 0;
    assert_eq!(
        sim.neighbors,
        find_neighbors_brute(&sim.state.positions, h),
        "initial neighbour lists"
    );

    let steps = 60;
    let mut previous = sim.neighbors.clone();
    let mut changed_steps = 0;
    for step in 1..=steps {
        sim.advance();
        assert_eq!(
            sim.neighbors,
            find_neighbors_brute(&sim.state.positions, h),
            "neighbour lists after step {step}"
        );
        if sim.neighbors != previous {
            changed_steps += 1;
        }
        previous = sim.neighbors.clone();
    }
    // The motion really was ballistic ...
    let elapsed = steps as f64 * dt;
    for (position, x0) in sim.state.positions.iter().zip(initial.iter()) {
        let expected = add3(*x0, scale3(sub3(*x0, centre), 2.0 * elapsed));
        assert!(norm3(sub3(*position, expected)) < 1e-12);
    }
    // ... and it moved pairs across the support radius, so the check above
    // exercised changing neighbour sets rather than a frozen one.
    assert!(
        changed_steps >= 3,
        "neighbour sets changed in only {changed_steps} of {steps} steps"
    );
}

/// The factor table is built from the same neighbour sets as the brute-force
/// reference, so its alphas agree exactly.
#[test]
fn test_factor_table_matches_brute_force_alphas() {
    let dx = 0.05;
    let h = 0.1;
    let rho0 = 1000.0;
    let positions = jittered_lattice(4, dx, 0x1405_7B7E_F767_814F);
    let n = positions.len();
    let masses = vec![rho0 * dx * dx * dx; n];
    let densities: Vec<f64> = (0..n).map(|i| rho0 * (0.9 + 0.001 * i as f64)).collect();
    let table = DfsphFactorTable::build(&positions, &masses, &densities, h, rho0);
    assert_eq!(table.len(), n);
    let brute = find_neighbors_brute(&positions, h);
    let mut non_trivial = 0usize;
    for i in 0..n {
        let nb_pos: Vec<[f64; 3]> = brute[i].iter().map(|&j| positions[j]).collect();
        let nb_mass: Vec<f64> = brute[i].iter().map(|&j| masses[j]).collect();
        let expected = compute_alpha_single(positions[i], densities[i], h, &nb_pos, &nb_mass);
        assert_eq!(table.get(i), expected, "alpha of particle {i}");
        if expected > 0.0 {
            non_trivial += 1;
        }
    }
    assert!(non_trivial > n / 2, "most alphas must be non-zero");
}
