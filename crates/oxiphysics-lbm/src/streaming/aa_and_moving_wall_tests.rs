//! Regression tests for the AA-pattern half-steps ([`stream_aa_2d`]) and the
//! half-way moving-wall bounce-back ([`stream_half_way_bb_moving_wall_2d`]).

use super::{
    max_diff_snapshots, snapshot_2d, stream_2d, stream_aa_2d, stream_half_way_bb_2d,
    stream_half_way_bb_moving_wall_2d,
};
use crate::collision::bgk_collide_2d;
use crate::grid::LbmGrid2D;
use crate::lattice::{CS2, D2Q9_WEIGHTS, LatticeType};

/// A D2Q9 grid with deterministic pseudo-random positive populations.
fn random_grid(nx: usize, ny: usize, seed: u64) -> LbmGrid2D {
    let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
    let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
    for plane in grid.f.iter_mut() {
        for value in plane.iter_mut() {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            *value = 0.01 + (state >> 11) as f64 / (1u64 << 53) as f64;
        }
    }
    grid
}

/// A D2Q9 grid without any population.
fn empty_grid(nx: usize, ny: usize) -> LbmGrid2D {
    let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
    for plane in grid.f.iter_mut() {
        plane.fill(0.0);
    }
    grid
}

/// Exchange the planes of every pair of opposite directions.
fn swap_opposite_planes(grid: &mut LbmGrid2D) {
    for i in 0..grid.lattice.q() {
        let opp = grid.lattice.opposite(i);
        if i < opp {
            grid.f.swap(i, opp);
        }
    }
}

fn total_population(grid: &LbmGrid2D) -> f64 {
    grid.f.iter().flatten().sum()
}

/// Total momentum `sum_i f_i c_i` of all populations, read with the standard
/// lattice velocities.
fn total_momentum(grid: &LbmGrid2D) -> [f64; 2] {
    let mut momentum = [0.0_f64; 2];
    for (i, plane) in grid.f.iter().enumerate() {
        let c = grid.lattice.velocity_2d(i);
        let sum: f64 = plane.iter().sum();
        momentum[0] += sum * f64::from(c[0]);
        momentum[1] += sum * f64::from(c[1]);
    }
    momentum
}

// ---------------------------------------------------------------------------
// AA pattern
// ---------------------------------------------------------------------------

#[test]
fn aa_even_step_streams_into_the_reversed_slots() {
    for (nx, ny) in [(1, 1), (2, 3), (7, 5), (8, 6), (13, 11)] {
        let mut aa = random_grid(nx, ny, (nx * 31 + ny) as u64);
        let mut reference = aa.clone();
        stream_aa_2d(&mut aa, true);
        stream_2d(&mut reference);
        swap_opposite_planes(&mut reference);
        let diff = max_diff_snapshots(&snapshot_2d(&aa), &snapshot_2d(&reference));
        assert!(diff < 1e-14, "even step on {nx}x{ny} differs by {diff}");
    }
}

#[test]
fn aa_even_then_odd_advances_the_populations_by_two_streaming_steps() {
    for (nx, ny) in [(1, 1), (2, 3), (7, 5), (8, 6), (13, 11)] {
        let mut aa = random_grid(nx, ny, (nx * 17 + ny * 5) as u64);
        let mut reference = aa.clone();
        stream_aa_2d(&mut aa, true);
        stream_aa_2d(&mut aa, false);
        stream_2d(&mut reference);
        stream_2d(&mut reference);
        let diff = max_diff_snapshots(&snapshot_2d(&aa), &snapshot_2d(&reference));
        assert!(diff < 1e-14, "even+odd on {nx}x{ny} differs by {diff}");
    }
}

#[test]
fn aa_half_steps_move_a_spike_one_cell_each() {
    let (nx, ny) = (9, 7);
    for i in 1..9 {
        let mut grid = empty_grid(nx, ny);
        let c = grid.lattice.velocity_2d(i);
        let opp = grid.lattice.opposite(i);
        let start = grid.idx(3, 2);
        grid.f[i][start] = 1.0;

        stream_aa_2d(&mut grid, true);
        let x1 = (3 + c[0]).rem_euclid(nx as i32) as usize;
        let y1 = (2 + c[1]).rem_euclid(ny as i32) as usize;
        assert_eq!(
            grid.f[opp][grid.idx(x1, y1)],
            1.0,
            "direction {i}: the even step stores the moved population in the slot of {opp}"
        );
        assert_eq!(total_population(&grid), 1.0, "direction {i}: even step");

        stream_aa_2d(&mut grid, false);
        let x2 = (3 + 2 * c[0]).rem_euclid(nx as i32) as usize;
        let y2 = (2 + 2 * c[1]).rem_euclid(ny as i32) as usize;
        assert_eq!(
            grid.f[i][grid.idx(x2, y2)],
            1.0,
            "direction {i}: the odd step restores the standard slot after a second move"
        );
        assert_eq!(total_population(&grid), 1.0, "direction {i}: odd step");
    }
}

#[test]
fn aa_pairs_of_collision_free_half_steps_conserve_mass_and_momentum() {
    let mut grid = random_grid(6, 5, 77);
    let mass = total_population(&grid);
    let momentum = total_momentum(&grid);
    stream_aa_2d(&mut grid, true);
    // In the reversed layout the same populations read as the opposite momentum.
    let reversed = total_momentum(&grid);
    assert!((total_population(&grid) - mass).abs() < 1e-12);
    assert!((reversed[0] + momentum[0]).abs() < 1e-12);
    assert!((reversed[1] + momentum[1]).abs() < 1e-12);
    stream_aa_2d(&mut grid, false);
    assert!((total_population(&grid) - mass).abs() < 1e-12);
    let restored = total_momentum(&grid);
    assert!((restored[0] - momentum[0]).abs() < 1e-12);
    assert!((restored[1] - momentum[1]).abs() < 1e-12);
}

#[test]
fn bgk_collision_commutes_with_the_reversed_layout() {
    // Colliding in the reversed layout (even step -> collide -> odd step) must
    // give the same populations as colliding in the standard layout, because
    // relaxing the slot of c_i towards feq(rho, -u) relaxes the population of
    // c_opp(i) towards feq(rho, u).
    let omega = 1.3;
    let mut reversed = random_grid(6, 5, 5);
    let mut standard = reversed.clone();

    stream_aa_2d(&mut reversed, true);
    bgk_collide_2d(&mut reversed, omega);
    stream_aa_2d(&mut reversed, false);

    stream_2d(&mut standard);
    bgk_collide_2d(&mut standard, omega);
    stream_2d(&mut standard);

    let diff = max_diff_snapshots(&snapshot_2d(&reversed), &snapshot_2d(&standard));
    assert!(
        diff < 1e-13,
        "collision in the reversed layout differs by {diff}"
    );
}

// ---------------------------------------------------------------------------
// Half-way moving-wall bounce-back
// ---------------------------------------------------------------------------

#[test]
fn moving_wall_reflection_subtracts_the_wall_momentum_term() {
    // Fluid node (2, ny - 2) below a wall at y = ny - 1 that moves along +x.
    let (nx, ny) = (5, 6);
    let u_wall = [0.04, 0.0];
    let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
    let src = grid.idx(2, ny - 2);
    let before: Vec<f64> = (0..9).map(|i| 0.05 + 0.01 * i as f64).collect();
    for (plane, &value) in grid.f.iter_mut().zip(&before) {
        plane[src] = value;
    }
    let rho: f64 = before.iter().sum();

    stream_half_way_bb_moving_wall_2d(&mut grid, |_, y| y == 0 || y == ny - 1, u_wall);

    // Links into the wall: north (2), north-east (5), north-west (6).
    let correction = 2.0 * D2Q9_WEIGHTS[5] * rho * u_wall[0] / CS2;
    assert!(correction > 0.0);
    assert!(
        (grid.f[4][src] - before[2]).abs() < 1e-15,
        "a link normal to the wall carries no momentum term"
    );
    assert!(
        (grid.f[7][src] - (before[5] - correction)).abs() < 1e-15,
        "north-east link: {} vs {}",
        grid.f[7][src],
        before[5] - correction
    );
    assert!(
        (grid.f[8][src] - (before[6] + correction)).abs() < 1e-15,
        "north-west link: {} vs {}",
        grid.f[8][src],
        before[6] + correction
    );
    // The +x wall momentum is handed to the fluid: the reflected populations
    // carry more +x momentum than the plain bounce-back by 2 * correction.
    let plain_px = -before[5] + before[6];
    let moving_px = -grid.f[7][src] + grid.f[8][src];
    assert!(
        (moving_px - plain_px - 2.0 * correction).abs() < 1e-15,
        "momentum gained from the wall: {}",
        moving_px - plain_px
    );
}

#[test]
fn moving_wall_with_zero_velocity_is_the_plain_half_way_bounce_back() {
    let (nx, ny) = (7, 6);
    let is_solid = |x: usize, y: usize| y == 0 || y == ny - 1 || (x == 3 && y == 3);
    let mut moving = random_grid(nx, ny, 11);
    let mut plain = moving.clone();
    stream_half_way_bb_moving_wall_2d(&mut moving, is_solid, [0.0, 0.0]);
    stream_half_way_bb_2d(&mut plain, is_solid);
    let diff = max_diff_snapshots(&snapshot_2d(&moving), &snapshot_2d(&plain));
    assert!(
        diff < 1e-15,
        "zero wall velocity must not change anything: {diff}"
    );
}

#[test]
fn walls_moving_along_plus_x_drag_the_fluid_along_plus_x() {
    let (nx, ny) = (4, 10);
    let u_wall = 0.05;
    let omega = 1.0 / 0.8;
    let is_solid = |_x: usize, y: usize| y == 0 || y == ny - 1;
    let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
    for _ in 0..3000 {
        bgk_collide_2d(&mut grid, omega);
        stream_half_way_bb_moving_wall_2d(&mut grid, is_solid, [u_wall, 0.0]);
    }
    grid.compute_macroscopic();

    let mut sum_ux = 0.0;
    let mut cells = 0usize;
    let mut max_deviation = 0.0_f64;
    let mut mass = 0.0;
    for y in 1..ny - 1 {
        for x in 0..nx {
            let k = grid.idx(x, y);
            sum_ux += grid.ux[k];
            cells += 1;
            max_deviation = max_deviation.max((grid.ux[k] - u_wall).abs());
            mass += grid.rho[k];
            assert!(grid.uy[k].abs() < 1e-9, "no wall-normal flow at ({x},{y})");
        }
    }
    let mean_ux = sum_ux / cells as f64;
    assert!(
        mean_ux > 0.0,
        "the fluid must move with the walls, got {mean_ux}"
    );
    assert!(
        (mean_ux - u_wall).abs() < 0.02 * u_wall,
        "mean ux {mean_ux} vs wall velocity {u_wall}"
    );
    // Two walls with the same velocity impose a uniform plug flow.
    assert!(max_deviation < 1e-6, "profile deviates by {max_deviation}");
    assert!(
        (mass / cells as f64 - 1.0).abs() < 1e-12,
        "tangential moving walls must conserve the fluid mass"
    );
}

#[test]
fn walls_moving_along_minus_x_drag_the_fluid_along_minus_x() {
    let (nx, ny) = (4, 8);
    let u_wall = -0.03;
    let omega = 1.0 / 0.7;
    let is_solid = |_x: usize, y: usize| y == 0 || y == ny - 1;
    let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
    for _ in 0..2000 {
        bgk_collide_2d(&mut grid, omega);
        stream_half_way_bb_moving_wall_2d(&mut grid, is_solid, [u_wall, 0.0]);
    }
    grid.compute_macroscopic();
    for y in 1..ny - 1 {
        let k = grid.idx(1, y);
        assert!(
            (grid.ux[k] - u_wall).abs() < 1e-6,
            "row {y}: ux = {} vs wall velocity {u_wall}",
            grid.ux[k]
        );
    }
}
