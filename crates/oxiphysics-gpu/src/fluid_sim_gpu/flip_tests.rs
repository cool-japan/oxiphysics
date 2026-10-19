// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the FLIP/PIC pipeline.
//!
//! The pipeline used to transfer only the x velocity (the y/z components were
//! never read from the grid), computed the FLIP delta against a clone of the
//! grid it had just projected (so it was identically zero), advected particles
//! by `dx` instead of `dt`, and never marked any cell as fluid, so the pressure
//! solve was a no-op.  The tests below pin every one of those stages.

use super::*;

/// Particles on a regular 2 x 2 x 2-per-cell lattice filling the cell range
/// `[lo, hi)` (cell units) of a grid with cell size `dx`.
fn lattice_particles(lo: [usize; 3], hi: [usize; 3], dx: f64) -> Vec<FlipParticle> {
    let mut out = Vec::new();
    for k in lo[2]..hi[2] {
        for j in lo[1]..hi[1] {
            for i in lo[0]..hi[0] {
                for sub in 0..8_usize {
                    let (a, b, c) = (sub & 1, (sub >> 1) & 1, (sub >> 2) & 1);
                    out.push(FlipParticle::new([
                        (i as f64 + 0.25 + 0.5 * a as f64) * dx,
                        (j as f64 + 0.25 + 0.5 * b as f64) * dx,
                        (k as f64 + 0.25 + 0.5 * c as f64) * dx,
                    ]));
                }
            }
        }
    }
    out
}

fn config(grid: [usize; 3], dt: f64, gravity: [f64; 3]) -> FluidSimConfig {
    FluidSimConfig {
        grid_size: grid,
        dx: 0.1,
        dt,
        gravity,
        ..FluidSimConfig::default()
    }
}

// ── staggered interpolation ──────────────────────────────────────────────────

/// Fill the three staggered components of `grid` from the analytic field
/// `field(component, x, y, z)` evaluated at the true face-centre positions.
fn fill_staggered(grid: &mut MacGrid, field: impl Fn(usize, f64, f64, f64) -> f64) {
    let dx = grid.dx;
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                let v = field(
                    0,
                    i as f64 * dx,
                    (j as f64 + 0.5) * dx,
                    (k as f64 + 0.5) * dx,
                );
                grid.set_u(i, j, k, v);
            }
        }
    }
    for k in 0..grid.nz {
        for j in 0..=grid.ny {
            for i in 0..grid.nx {
                let v = field(
                    1,
                    (i as f64 + 0.5) * dx,
                    j as f64 * dx,
                    (k as f64 + 0.5) * dx,
                );
                grid.set_v(i, j, k, v);
            }
        }
    }
    for k in 0..=grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let v = field(
                    2,
                    (i as f64 + 0.5) * dx,
                    (j as f64 + 0.5) * dx,
                    k as f64 * dx,
                );
                grid.set_w(i, j, k, v);
            }
        }
    }
}

/// Trilinear interpolation reproduces affine fields exactly, which only holds
/// when the sample positions assumed by `interp_u/v/w` are the true face
/// centres of the staggered grid.
#[test]
fn staggered_interpolators_reproduce_affine_fields() {
    let mut grid = MacGrid::new(6, 5, 4, 0.25);
    let field = |c: usize, x: f64, y: f64, z: f64| match c {
        0 => 2.0 * x + 3.0 * y - 5.0 * z + 1.0,
        1 => -1.5 * x + 0.5 * y + 4.0 * z - 2.0,
        _ => 7.0 * x - 2.0 * y + 1.0 * z + 0.25,
    };
    fill_staggered(&mut grid, field);
    let points = [
        [0.6, 0.55, 0.45],
        [1.1, 0.9, 0.7],
        [0.3, 0.3, 0.3],
        [1.2, 1.1, 0.6],
    ];
    for [x, y, z] in points {
        let got = [
            grid.interp_u(x, y, z),
            grid.interp_v(x, y, z),
            grid.interp_w(x, y, z),
        ];
        for (c, &value) in got.iter().enumerate() {
            let expected = field(c, x, y, z);
            assert!(
                (value - expected).abs() < 1e-12,
                "component {c} at ({x}, {y}, {z}): {value} vs {expected}"
            );
        }
    }
}

#[test]
fn interpolators_clamp_outside_the_grid() {
    let mut grid = MacGrid::new(3, 3, 3, 0.5);
    fill_staggered(&mut grid, |c, _, _, _| 1.0 + c as f64);
    for c in [-10.0, 10.0] {
        assert!((grid.interp_u(c, c, c) - 1.0).abs() < 1e-12);
        assert!((grid.interp_v(c, c, c) - 2.0).abs() < 1e-12);
        assert!((grid.interp_w(c, c, c) - 3.0).abs() < 1e-12);
    }
}

/// P2G followed by G2P (pure PIC) returns the original velocity of particles
/// that are well inside the grid, on all three components.
#[test]
fn p2g_then_pic_g2p_reproduces_a_linear_velocity_field() {
    let dx = 0.1;
    let mut particles = lattice_particles([0, 0, 0], [10, 10, 10], dx);
    let field = |x: f64, y: f64, z: f64| [2.0 * y + 1.0, 3.0 * x - z + 0.5, x + 4.0 * y + 2.0];
    for p in particles.iter_mut() {
        p.velocity = field(p.position[0], p.position[1], p.position[2]);
    }
    let mut grid = MacGrid::new(10, 10, 10, dx);
    p2g_transfer(&particles, &mut grid);
    let reference = grid.clone();

    let mut moved = particles.clone();
    for p in moved.iter_mut() {
        p.velocity = [0.0; 3];
    }
    // dt = 0: pure PIC resampling without advection.
    g2p_transfer(&mut moved, &grid, &reference, 0.0, 0.0);
    let mut checked = 0;
    for (orig, p) in particles.iter().zip(moved.iter()) {
        let interior = orig
            .position
            .iter()
            .all(|&c| c > 2.0 * dx && c < 10.0 * dx - 2.0 * dx);
        if !interior {
            continue;
        }
        checked += 1;
        for c in 0..3 {
            assert!(
                (p.velocity[c] - orig.velocity[c]).abs() < 1e-10,
                "component {c} at {:?}: {} vs {}",
                orig.position,
                p.velocity[c],
                orig.velocity[c]
            );
        }
    }
    assert!(checked > 100, "only {checked} interior particles checked");
}

// ── G2P ──────────────────────────────────────────────────────────────────────

fn uniform_grid(n: usize, velocity: [f64; 3]) -> MacGrid {
    let mut grid = MacGrid::new(n, n, n, 0.1);
    grid.u.fill(velocity[0]);
    grid.v.fill(velocity[1]);
    grid.w.fill(velocity[2]);
    grid
}

#[test]
fn g2p_blends_pic_and_flip_on_all_three_axes() {
    let grid_new = uniform_grid(8, [1.0, 2.0, 3.0]);
    let grid_old = uniform_grid(8, [0.0, 0.0, 0.0]);
    let v0 = [10.0, 20.0, 30.0];
    for (ratio, expected) in [
        (0.0, [1.0, 2.0, 3.0]),
        (1.0, [11.0, 22.0, 33.0]),
        (0.5, [6.0, 12.0, 18.0]),
    ] {
        let mut particles = vec![FlipParticle::new([0.4, 0.4, 0.4])];
        particles[0].velocity = v0;
        g2p_transfer(&mut particles, &grid_new, &grid_old, ratio, 0.01);
        for (c, &want) in expected.iter().enumerate() {
            assert!(
                (particles[0].velocity[c] - want).abs() < 1e-12,
                "flip_ratio {ratio}, component {c}: {} vs {want}",
                particles[0].velocity[c]
            );
        }
        // The sampled grid velocity is recorded on the particle.
        assert_eq!(particles[0].velocity_grid, [1.0, 2.0, 3.0]);
    }
}

#[test]
fn g2p_advects_by_dt_not_by_dx() {
    let grid = uniform_grid(8, [1.0, -2.0, 0.5]);
    let zero = uniform_grid(8, [0.0; 3]);
    let mut particles = vec![FlipParticle::new([0.4, 0.4, 0.4])];
    let dt = 0.02;
    g2p_transfer(&mut particles, &grid, &zero, 0.0, dt);
    let expected = [0.4 + 1.0 * dt, 0.4 - 2.0 * dt, 0.4 + 0.5 * dt];
    for (c, &want) in expected.iter().enumerate() {
        assert!(
            (particles[0].position[c] - want).abs() < 1e-12,
            "component {c}: {} vs {want}",
            particles[0].position[c]
        );
    }
}

#[test]
fn g2p_keeps_particles_inside_the_domain() {
    let grid = uniform_grid(8, [5.0, -5.0, 5.0]);
    let zero = uniform_grid(8, [0.0; 3]);
    let extent = 0.8;
    let mut particles = vec![FlipParticle::new([0.79, 0.01, 0.79])];
    g2p_transfer(&mut particles, &grid, &zero, 0.0, 0.1);
    let p = &particles[0];
    assert_eq!(p.position[0], extent);
    assert_eq!(p.position[1], 0.0);
    assert_eq!(p.position[2], extent);
    // The outward velocity component is removed, the others are untouched.
    assert_eq!(p.velocity[0], 0.0);
    assert_eq!(p.velocity[1], 0.0);
    assert_eq!(p.velocity[2], 0.0);
}

// ── fluid-cell marking ───────────────────────────────────────────────────────

#[test]
fn mark_fluid_cells_flags_cells_containing_particles() {
    let mut grid = MacGrid::new(4, 4, 4, 1.0);
    grid.p.fill(5.0);
    let idx = |i: usize, j: usize, k: usize| i + 4 * (j + 4 * k);
    // A solid cell keeps its flag even when a particle sits in it.
    let solid_idx = idx(1, 1, 1);
    grid.flags[solid_idx] = CELL_SOLID;
    // A stale fluid flag with no particle becomes air.
    grid.flags[0] = CELL_FLUID;
    let particles = vec![
        FlipParticle::new([3.9, 1.2, 2.5]),
        FlipParticle::new([1.5, 1.5, 1.5]),
        FlipParticle::new([-0.3, 7.0, 100.0]),
    ];
    grid.mark_fluid_cells(&particles);

    assert_eq!(grid.flags[idx(3, 1, 2)], CELL_FLUID);
    assert_eq!(
        grid.flags[idx(0, 3, 3)],
        CELL_FLUID,
        "positions are clamped"
    );
    assert_eq!(grid.flags[solid_idx], CELL_SOLID);
    assert_eq!(grid.flags[0], CELL_AIR);
    let fluid = grid.flags.iter().filter(|&&f| f == CELL_FLUID).count();
    assert_eq!(fluid, 2);
    // Pressure survives only in fluid cells.
    for (cell, (&flag, &p)) in grid.flags.iter().zip(grid.p.iter()).enumerate() {
        if flag == CELL_FLUID {
            assert_eq!(p, 5.0, "cell {cell}");
        } else {
            assert_eq!(p, 0.0, "cell {cell}");
        }
    }
}

#[test]
fn mark_fluid_cells_handles_nan_positions_and_empty_input() {
    let mut grid = MacGrid::new(3, 3, 3, 0.5);
    grid.mark_fluid_cells(&[]);
    assert!(grid.flags.iter().all(|&f| f == CELL_AIR));
    grid.mark_fluid_cells(&[FlipParticle::new([f64::NAN, 0.2, 0.2])]);
    assert_eq!(grid.flags.iter().filter(|&&f| f == CELL_FLUID).count(), 1);
}

// ── one full step ────────────────────────────────────────────────────────────

/// One particle in the middle of the box: nothing but gravity acts on it, so
/// after one step `v = g dt` on every axis, whatever the PIC/FLIP blend.
#[test]
fn single_particle_accelerates_with_gravity_on_every_axis() {
    let dt = 0.01;
    for gravity in [
        [0.0, -9.81, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 0.0, -3.0],
        [1.0, -2.0, 3.0],
    ] {
        for flip_ratio in [0.0, 0.95, 1.0] {
            let cfg = FluidSimConfig {
                flip_ratio,
                ..config([10, 10, 10], dt, gravity)
            };
            let mut sim = FluidSimulation::new(cfg);
            sim.particles.push(FlipParticle::new([0.5, 0.5, 0.5]));
            sim.step();
            let p = &sim.particles[0];
            for c in 0..3 {
                let v = gravity[c] * dt;
                assert!(
                    (p.velocity[c] - v).abs() < 1e-12,
                    "gravity {gravity:?}, ratio {flip_ratio}, component {c}: {} vs {v}",
                    p.velocity[c]
                );
                let x = 0.5 + v * dt;
                assert!(
                    (p.position[c] - x).abs() < 1e-12,
                    "gravity {gravity:?}, component {c}: position {} vs {x}",
                    p.position[c]
                );
            }
        }
    }
}

/// The FLIP delta is measured against the grid right after P2G, so the
/// particle keeps its own velocity and adds only the grid change.  (Vorticity
/// confinement is switched off: the splatted single-particle field has curl.)
#[test]
fn flip_delta_is_taken_against_the_post_p2g_grid() {
    let dt = 0.01;
    let cfg = FluidSimConfig {
        flip_ratio: 1.0,
        vorticity_eps: 0.0,
        ..config([10, 10, 10], dt, [0.0, -9.81, 0.0])
    };
    let mut sim = FluidSimulation::new(cfg);
    // Off the face planes, so every face of the particle's cell receives
    // splat weight and the splatted field is divergence-free in that cell.
    let at = [0.53, 0.47, 0.51];
    sim.particles.push(FlipParticle::new(at));
    sim.particles[0].velocity = [0.3, -0.2, 0.1];
    sim.step();
    let v = sim.particles[0].velocity;
    assert!((v[0] - 0.3).abs() < 1e-12, "{v:?}");
    assert!((v[1] - (-0.2 - 9.81 * dt)).abs() < 1e-12, "{v:?}");
    assert!((v[2] - 0.1).abs() < 1e-12, "{v:?}");
    // grid_prev holds the splatted particle velocity, not the end-of-step one.
    assert!((sim.grid_prev.interp_u(at[0], at[1], at[2]) - 0.3).abs() < 1e-12);
    assert!((sim.grid_prev.interp_v(at[0], at[1], at[2]) + 0.2).abs() < 1e-12);
    assert!((sim.grid.interp_v(at[0], at[1], at[2]) - v[1]).abs() < 1e-12);
}

/// PIC smooths particle velocities through the grid, FLIP keeps them: two
/// particles sharing a cell with opposite velocities.
#[test]
fn flip_preserves_particle_velocity_differences_that_pic_smooths_out() {
    let spread = |flip_ratio: f64| -> f64 {
        let cfg = FluidSimConfig {
            flip_ratio,
            vorticity_eps: 0.0,
            ..config([10, 10, 10], 0.001, [0.0; 3])
        };
        let mut sim = FluidSimulation::new(cfg);
        let mut a = FlipParticle::new([0.52, 0.5, 0.5]);
        let mut b = FlipParticle::new([0.58, 0.5, 0.5]);
        a.velocity = [1.0, 0.0, 0.0];
        b.velocity = [-1.0, 0.0, 0.0];
        sim.particles.push(a);
        sim.particles.push(b);
        sim.step();
        sim.particles[0].velocity[0] - sim.particles[1].velocity[0]
    };
    let flip = spread(1.0);
    let pic = spread(0.0);
    assert!(flip > 1.5, "FLIP lost the velocity difference: {flip}");
    assert!(
        pic < 0.5 * flip,
        "PIC must smooth the difference: pic {pic} vs flip {flip}"
    );
}

/// A water column that fills the box width and rests on the floor: the
/// discrete hydrostatic solution is `p_j = (n - j) rho g dx` per layer, the
/// pressure projection cancels gravity, and the particles stay at rest.
#[test]
fn resting_water_column_builds_hydrostatic_pressure_and_stays_at_rest() {
    let dt = 0.01;
    let n_layers = 3;
    let cfg = FluidSimConfig {
        pressure_iters: 600,
        vorticity_eps: 0.0,
        ..config([4, 6, 4], dt, [0.0, -9.81, 0.0])
    };
    let (rho, dx) = (cfg.density, cfg.dx);
    let mut sim = FluidSimulation::new(cfg);
    sim.particles = lattice_particles([0, 0, 0], [4, n_layers, 4], dx);
    let before: Vec<[f64; 3]> = sim.particles.iter().map(|p| p.position).collect();
    sim.step();

    for k in 0..4 {
        for j in 0..6 {
            for i in 0..4 {
                let idx = i + 4 * (j + 6 * k);
                let p = sim.grid.p[idx];
                if j < n_layers {
                    assert_eq!(sim.grid.flags[idx], CELL_FLUID);
                    let expected = (n_layers - j) as f64 * rho * 9.81 * dx;
                    assert!(
                        p > 0.0 && (p - expected).abs() < 1e-3 * expected,
                        "cell ({i},{j},{k}): p = {p}, hydrostatic {expected}"
                    );
                } else {
                    assert_eq!(p, 0.0, "air cell ({i},{j},{k}) must stay at p = 0");
                }
            }
        }
    }
    // Gravity is balanced: no particle accelerates or drifts.
    for (p, x0) in sim.particles.iter().zip(before.iter()) {
        for c in 0..3 {
            assert!(
                p.velocity[c].abs() < 1e-3 * 9.81 * dt,
                "particle at {x0:?} accelerated: {:?}",
                p.velocity
            );
            assert!((p.position[c] - x0[c]).abs() < 1e-6);
        }
    }
}

/// A block of fluid resting on the floor develops non-zero pressure in every
/// fluid cell after a single default step; a free-falling block does not.
#[test]
fn block_on_floor_has_pressure_but_free_falling_block_has_none() {
    let dx = 0.1;
    let cfg = config([12, 12, 12], 0.01, [0.0, -9.81, 0.0]);

    let mut resting = FluidSimulation::new(cfg.clone());
    resting.particles = lattice_particles([4, 0, 4], [8, 3, 8], dx);
    resting.step();
    let mut fluid_cells = 0;
    for (&flag, &p) in resting.grid.flags.iter().zip(resting.grid.p.iter()) {
        if flag == CELL_FLUID {
            fluid_cells += 1;
            assert!(p > 0.0, "fluid cell without pressure");
        } else {
            assert_eq!(p, 0.0);
        }
    }
    assert_eq!(fluid_cells, 4 * 3 * 4);

    let mut falling = FluidSimulation::new(cfg);
    falling.particles = lattice_particles([4, 5, 4], [8, 8, 8], dx);
    falling.step();
    assert!(falling.grid.p.iter().all(|&p| p.abs() < 1e-9));
}

// ── projection ───────────────────────────────────────────────────────────────

/// Largest `|div|` over the fluid cells.
fn max_fluid_divergence(grid: &mut MacGrid) -> f64 {
    grid.compute_divergence();
    grid.div
        .iter()
        .zip(grid.flags.iter())
        .filter(|&(_, &f)| f == CELL_FLUID)
        .fold(0.0_f64, |m, (&d, _)| m.max(d.abs()))
}

fn divergent_block_grid() -> MacGrid {
    let mut grid = MacGrid::new(8, 8, 8, 0.1);
    for k in 2..6 {
        for j in 2..6 {
            for i in 2..6 {
                let idx = i + 8 * (j + 8 * k);
                grid.flags[idx] = CELL_FLUID;
            }
        }
    }
    // Deterministic, strongly divergent velocity field on every face.
    for (n, u) in grid.u.iter_mut().enumerate() {
        *u = ((n * 37 % 17) as f64 - 8.0) * 0.1;
    }
    for (n, v) in grid.v.iter_mut().enumerate() {
        *v = ((n * 53 % 19) as f64 - 9.0) * 0.1;
    }
    for (n, w) in grid.w.iter_mut().enumerate() {
        *w = ((n * 71 % 23) as f64 - 11.0) * 0.1;
    }
    grid.enforce_solid_boundaries();
    grid
}

#[test]
fn projection_reduces_fluid_divergence_by_more_than_ten() {
    for iterations in [20, 300] {
        let mut grid = divergent_block_grid();
        let before = max_fluid_divergence(&mut grid);
        assert!(before > 1.0, "test setup: divergence {before}");
        grid.jacobi_pressure_solve(1000.0, 0.01, iterations);
        grid.pressure_project(1000.0, 0.01);
        let after = max_fluid_divergence(&mut grid);
        assert!(
            after * 10.0 < before,
            "{iterations} sweeps: divergence {before} -> {after}"
        );
        if iterations >= 300 {
            assert!(after * 1e4 < before, "converged solve: {before} -> {after}");
        }
        // The solve only produces pressure in the fluid.
        for (&flag, &p) in grid.flags.iter().zip(grid.p.iter()) {
            if flag != CELL_FLUID {
                assert_eq!(p, 0.0);
            }
        }
        assert!(grid.p.iter().any(|&p| p != 0.0));
    }
}

/// 1-D chain `air | fluid | fluid | solid` closed by the domain wall: air is a
/// Dirichlet (`p = 0`) neighbour, the solid and the wall are Neumann.  The
/// discrete solution for inflow `a` through the air face and `b` through the
/// middle face is `p1 = s a / dx`, `p2 = s (a + b) / dx` with
/// `s = rho dx^2 / dt`, and the projection removes both fluxes.
#[test]
fn jacobi_uses_dirichlet_air_and_neumann_solid_neighbours() {
    let (rho, dt, dx) = (1000.0, 0.01, 0.1);
    let mut grid = MacGrid::new(4, 1, 1, dx);
    grid.flags = vec![CELL_AIR, CELL_FLUID, CELL_FLUID, CELL_SOLID];
    let (a, b) = (0.7, -0.4);
    grid.set_u(1, 0, 0, a);
    grid.set_u(2, 0, 0, b);
    grid.set_u(3, 0, 0, 0.9); // face into the solid: must be removed
    grid.enforce_solid_boundaries();
    assert_eq!(grid.get_u(3, 0, 0), 0.0);

    grid.jacobi_pressure_solve(rho, dt, 4000);
    let s = rho * dx * dx / dt;
    let (p1, p2) = (s * a / dx, s * (a + b) / dx);
    assert_eq!(grid.p[0], 0.0, "air stays at p = 0");
    assert_eq!(grid.p[3], 0.0, "solid carries no pressure");
    assert!(
        (grid.p[1] - p1).abs() < 1e-6 * p1.abs(),
        "p1 = {}",
        grid.p[1]
    );
    assert!(
        (grid.p[2] - p2).abs() < 1e-6 * p2.abs(),
        "p2 = {}",
        grid.p[2]
    );

    grid.pressure_project(rho, dt);
    assert!(grid.get_u(1, 0, 0).abs() < 1e-8, "{}", grid.get_u(1, 0, 0));
    assert!(grid.get_u(2, 0, 0).abs() < 1e-8, "{}", grid.get_u(2, 0, 0));
    assert_eq!(grid.get_u(3, 0, 0), 0.0, "solid faces are left alone");
}

/// The free function and the grid method are the same solver.
#[test]
fn gpu_jacobi_matches_the_grid_solver() {
    let mut a = divergent_block_grid();
    let mut b = a.clone();
    a.jacobi_pressure_solve(1000.0, 0.01, 40);
    gpu_jacobi_pressure_solve(&mut b, 1000.0, 0.01, 40);
    assert_eq!(a.p, b.p);
}

#[test]
fn solid_boundaries_zero_normal_velocity_on_walls_and_solids() {
    let mut grid = MacGrid::new(3, 3, 3, 0.1);
    grid.u.fill(1.0);
    grid.v.fill(2.0);
    grid.w.fill(3.0);
    // Centre cell (1, 1, 1): index 1 + 3 * 1 + 9 * 1.
    let centre = 1 + 3 + 9;
    grid.flags[centre] = CELL_SOLID;
    grid.enforce_solid_boundaries();
    // Domain walls.
    assert_eq!(grid.get_u(0, 1, 1), 0.0);
    assert_eq!(grid.get_u(3, 1, 1), 0.0);
    assert_eq!(grid.get_v(1, 0, 1), 0.0);
    assert_eq!(grid.get_v(1, 3, 1), 0.0);
    assert_eq!(grid.get_w(1, 1, 0), 0.0);
    assert_eq!(grid.get_w(1, 1, 3), 0.0);
    // Faces of the solid centre cell.
    assert_eq!(grid.get_u(1, 1, 1), 0.0);
    assert_eq!(grid.get_u(2, 1, 1), 0.0);
    assert_eq!(grid.get_v(1, 1, 1), 0.0);
    assert_eq!(grid.get_v(1, 2, 1), 0.0);
    assert_eq!(grid.get_w(1, 1, 1), 0.0);
    assert_eq!(grid.get_w(1, 1, 2), 0.0);
    // Interior faces away from walls and solids keep their value.
    assert_eq!(grid.get_u(1, 0, 0), 1.0);
    assert_eq!(grid.get_v(0, 1, 0), 2.0);
    assert_eq!(grid.get_w(0, 0, 1), 3.0);
}

#[test]
fn body_acceleration_reaches_all_three_components() {
    let mut grid = MacGrid::new(3, 3, 3, 0.1);
    grid.add_body_acceleration([1.0, -2.0, 3.0], 0.5);
    assert!(grid.u.iter().all(|&u| (u - 0.5).abs() < 1e-15));
    assert!(grid.v.iter().all(|&v| (v + 1.0).abs() < 1e-15));
    assert!(grid.w.iter().all(|&w| (w - 1.5).abs() < 1e-15));
}

#[test]
fn copy_velocities_rejects_mismatched_grids() {
    let mut a = MacGrid::new(3, 3, 3, 0.1);
    let mut b = MacGrid::new(4, 3, 3, 0.1);
    b.u.fill(1.0);
    assert!(!a.copy_velocities_from(&b));
    assert!(a.u.iter().all(|&u| u == 0.0));
    let mut c = MacGrid::new(3, 3, 3, 0.1);
    c.u.fill(2.0);
    c.v.fill(3.0);
    c.w.fill(4.0);
    assert!(a.copy_velocities_from(&c));
    assert_eq!(a.u, c.u);
    assert_eq!(a.v, c.v);
    assert_eq!(a.w, c.w);
}

/// `FluidSimConfig::surface_tension` is documented as not read by
/// `FluidSimulation::step`: two simulations that differ only in that
/// coefficient stay bit-identical.
#[test]
fn flip_step_does_not_read_the_surface_tension_coefficient() {
    let run = |surface_tension: f64| -> Vec<([f64; 3], [f64; 3])> {
        let cfg = FluidSimConfig {
            surface_tension,
            ..config([8, 8, 8], 0.01, [0.0, -9.81, 0.0])
        };
        let mut sim = FluidSimulation::new(cfg);
        sim.particles = lattice_particles([2, 0, 2], [5, 3, 5], 0.1);
        for _ in 0..3 {
            sim.step();
        }
        sim.particles
            .iter()
            .map(|p| (p.position, p.velocity))
            .collect()
    };
    let without = run(0.0);
    let with = run(25.0);
    assert!(!without.is_empty());
    assert_eq!(without, with);
}
