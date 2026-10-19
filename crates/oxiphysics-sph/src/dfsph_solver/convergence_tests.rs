// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Iteration-by-iteration behaviour of the divergence-free and constant-density
// solves of `DfsphSolver`.
//
// Each iteration re-derives the per-particle stiffnesses from the *current*
// velocities and applies them once (Jacobi). For two equal particles the
// compressive error is therefore multiplied by exactly `1 − 2ω` per iteration,
// and the error of a compressed block falls with every iteration. The earlier
// code accumulated the stiffnesses and re-applied the running sum in every
// iteration (divergence solve) or read a density that nothing in the loop
// updated (constant-density solve), so neither error decreased.
//
// The recursion lives *inside* one `solve_*` call, so every check below runs a
// fresh solver for `k` iterations of a single call (`max_*_iterations = k`, and
// the tolerance 0 that can never be met) instead of calling a one-iteration
// solve `k` times, which would hide any state carried between the iterations.

use std::time::{Duration, Instant};

use super::*;

/// Rest density of every test (kg/m³).
const REST_DENSITY: f64 = 1000.0;

/// Lattice spacing at which the particle mass `ρ₀ s³` gives the rest density.
const REFERENCE_SPACING: f64 = 0.05;

/// Time step of the solves (s).
const DT: f64 = 1e-3;

/// Configuration of a solver whose `solve_*` calls never exit early: the
/// tolerance 0 can never be met, so a call performs exactly
/// `max_*_iterations` projections (one by default) and returns the residual
/// that is left after the last of them.
fn one_iteration_config(omega: f64, smoothing_length: f64) -> DfsphConfig {
    DfsphConfig {
        rest_density: REST_DENSITY,
        smoothing_length,
        support_radius: 2.0 * smoothing_length,
        max_pressure_iterations: 1,
        max_divergence_iterations: 1,
        error_threshold: 0.0,
        omega,
        viscosity: 0.0,
        gravity: [0.0; 3],
        surface_tension: 0.0,
    }
}

/// Neighbour lists, summation densities and alphas of the current positions.
fn prepare(solver: &mut DfsphSolver) {
    solver.build_spatial_hash();
    solver.build_neighbor_lists();
    solver.update_densities();
    solver.update_alphas();
}

/// Two equal, over-dense particles 0.06 apart (support `0.1`) that approach each
/// other at 1 m/s.
fn pair_solver(omega: f64) -> DfsphSolver {
    let mut solver = DfsphSolver::new(one_iteration_config(omega, 0.05));
    solver.add_particle([0.0; 3], [0.5, 0.0, 0.0], 0.3);
    solver.add_particle([0.06, 0.0, 0.0], [-0.5, 0.0, 0.0], 0.3);
    prepare(&mut solver);
    solver
}

/// Uniformly compressing velocity field `v = −rate (x − c)` of an `n³` lattice
/// with spacing `spacing` centred on `c` (divergence `−3 rate`).
fn compressing_field(n: usize, spacing: f64, rate: f64) -> impl Fn([f64; 3]) -> [f64; 3] {
    let centre = 0.5 * (n as f64 - 1.0) * spacing;
    move |position| {
        [
            -rate * (position[0] - centre),
            -rate * (position[1] - centre),
            -rate * (position[2] - centre),
        ]
    }
}

/// Solver holding an `n³` lattice with spacing `spacing` (particle mass
/// `ρ₀ REFERENCE_SPACING³`, so a spacing below the reference spacing is
/// over-dense) and the velocity field `velocity(position)`.
fn lattice_solver(
    config: DfsphConfig,
    n: usize,
    spacing: f64,
    velocity: impl Fn([f64; 3]) -> [f64; 3],
) -> DfsphSolver {
    let mut solver = DfsphSolver::new(config);
    let mass = REST_DENSITY * REFERENCE_SPACING.powi(3);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let position = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                solver.add_particle(position, velocity(position), mass);
            }
        }
    }
    prepare(&mut solver);
    solver
}

/// Mean compressive divergence (1/s) of the current velocities.
fn divergence_error(solver: &mut DfsphSolver) -> f64 {
    let kernel = WcSphDensity::new(KernelType::WendlandC2, solver.config.smoothing_length);
    solver.divergence_stiffnesses(DT, &kernel).1
}

/// Mean compressive relative density error of the current velocities.
fn density_error(solver: &DfsphSolver) -> f64 {
    let kernel = WcSphDensity::new(KernelType::WendlandC2, solver.config.smoothing_length);
    solver.density_stiffnesses(DT, &kernel).1
}

/// Error that one `solve_divergence_free` call leaves after exactly
/// `iterations` projections of a freshly built solver.
fn divergence_error_after(mut solver: DfsphSolver, iterations: usize) -> f64 {
    solver.config.max_divergence_iterations = iterations;
    let error = solver.solve_divergence_free(DT);
    assert_eq!(solver.stats.divergence_iterations, vec![iterations]);
    assert_eq!(solver.stats.divergence_error_history, vec![error]);
    error
}

/// Error that one `solve_constant_density` call leaves after exactly
/// `iterations` projections of a freshly built solver.
fn density_error_after(mut solver: DfsphSolver, iterations: usize) -> f64 {
    solver.config.max_pressure_iterations = iterations;
    let error = solver.solve_constant_density(DT);
    assert_eq!(solver.stats.pressure_iterations, vec![iterations]);
    assert_eq!(solver.stats.density_error_history, vec![error]);
    error
}

/// With the stiffnesses assigned (not accumulated) the divergence of a pair
/// shrinks by exactly `1 − 2ω` in *every* iteration of one solve. Re-applying
/// the running sum of the stiffnesses (the old behaviour) breaks the recursion
/// from the second iteration on.
#[test]
fn test_divergence_solve_of_a_pair_shrinks_the_error_by_one_minus_two_omega_per_iteration() {
    let omega = 0.25;
    let initial = divergence_error(&mut pair_solver(omega));
    assert!(initial > 1.0, "the pair must be compressing: {initial}");
    for iterations in 1..=8usize {
        let error = divergence_error_after(pair_solver(omega), iterations);
        let expected = initial * (1.0 - 2.0 * omega).powi(iterations as i32);
        assert!(
            (error - expected).abs() < 1e-9 * initial,
            "{iterations} iterations: error {error}, expected {expected}"
        );
    }
}

/// The constant-density solve re-predicts the density from the corrected
/// velocities each iteration, so the compressive error of a pair also shrinks
/// by exactly `1 − 2ω` per iteration (a frozen density gave a constant error and
/// an accumulated pressure an ever larger kick).
#[test]
fn test_constant_density_solve_of_a_pair_shrinks_the_error_by_one_minus_two_omega_per_iteration() {
    let omega = 0.25;
    let initial = density_error(&pair_solver(omega));
    assert!(initial > 0.05, "the pair must be over-dense: {initial}");
    for iterations in 1..=8usize {
        let error = density_error_after(pair_solver(omega), iterations);
        let expected = initial * (1.0 - 2.0 * omega).powi(iterations as i32);
        assert!(
            (error - expected).abs() < 1e-9 * initial,
            "{iterations} iterations: error {error}, expected {expected}"
        );
    }
}

/// A uniformly compressing lattice (`∇·v = −6 s⁻¹`): the mean compressive
/// divergence left by a solve falls strictly with the number of iterations, for
/// a conservative and for the default relaxation factor.
#[test]
fn test_divergence_solve_error_of_a_compressing_lattice_decreases_every_iteration() {
    let n = 6;
    for omega in [0.25, 0.5] {
        let build = || {
            let config = one_iteration_config(omega, 1.3 * REFERENCE_SPACING);
            let field = compressing_field(n, REFERENCE_SPACING, 2.0);
            lattice_solver(config, n, REFERENCE_SPACING, field)
        };
        let initial = divergence_error(&mut build());
        assert!(initial > 1.0, "omega {omega}: initial error {initial}");
        let mut previous = initial;
        for iterations in 1..=15usize {
            let error = divergence_error_after(build(), iterations);
            assert!(
                error > 0.0 && error < previous,
                "omega {omega}, {iterations} iterations: {previous} -> {error}"
            );
            previous = error;
        }
        assert!(
            previous < 0.02 * initial,
            "omega {omega}: only {initial} -> {previous} after 15 iterations"
        );
    }
}

/// A block compressed to about 1.2 ρ₀ (spacing 0.047 instead of 0.05): the
/// compressive density error left by a solve falls strictly with the number of
/// iterations.
#[test]
fn test_constant_density_solve_error_of_a_compressed_block_decreases_every_iteration() {
    let n = 6;
    for omega in [0.25, 0.5] {
        let build = || {
            let config = one_iteration_config(omega, 1.3 * REFERENCE_SPACING);
            lattice_solver(config, n, 0.047, |_| [0.0; 3])
        };
        let initial = density_error(&build());
        assert!(initial > 0.03, "omega {omega}: initial error {initial}");
        let mut previous = initial;
        for iterations in 1..=15usize {
            let error = density_error_after(build(), iterations);
            assert!(
                error > 0.0 && error < previous,
                "omega {omega}, {iterations} iterations: {previous} -> {error}"
            );
            previous = error;
        }
        assert!(
            previous < 0.02 * initial,
            "omega {omega}: only {initial} -> {previous} after 15 iterations"
        );
    }
}

/// With the real tolerance both solves reach it well before the iteration cap
/// (the frozen-density solve always ran into the cap) and record the iteration
/// count they actually used.
#[test]
fn test_compressed_block_converges_before_the_iteration_cap() {
    let n = 6;
    let cap = 50;
    let mut config = one_iteration_config(0.5, 1.3 * REFERENCE_SPACING);
    config.error_threshold = 1e-3;
    config.max_pressure_iterations = cap;
    config.max_divergence_iterations = cap;

    let mut block = lattice_solver(config.clone(), n, 0.047, |_| [0.0; 3]);
    assert!(density_error(&block) > 0.03);
    let error = block.solve_constant_density(DT);
    assert!(error < 1e-3, "density error {error}");
    assert_eq!(block.stats.pressure_iterations.len(), 1);
    let iterations = block.stats.pressure_iterations[0];
    assert!(
        (2..cap).contains(&iterations),
        "constant density: {iterations} iterations"
    );

    let field = compressing_field(n, REFERENCE_SPACING, 2.0);
    let mut lattice = lattice_solver(config, n, REFERENCE_SPACING, field);
    assert!(divergence_error(&mut lattice) > 1.0);
    let error = lattice.solve_divergence_free(DT);
    assert!(error < 1e-3, "divergence error {error}");
    assert_eq!(lattice.stats.divergence_iterations.len(), 1);
    let iterations = lattice.stats.divergence_iterations[0];
    assert!(
        (2..cap).contains(&iterations),
        "divergence-free: {iterations} iterations"
    );
}

/// A complete step of a 22³ = 10 648 particle block (about 70 neighbours per
/// particle, 20 iterations in each solve). Cloning the neighbour table for
/// every particle and iteration made this quadratic in the particle count (tens
/// of gigabytes copied per iteration); the step now takes a few seconds even on
/// a heavily loaded machine, so the generous bound only trips on such a
/// regression.
#[test]
fn test_step_of_a_ten_thousand_particle_block_is_fast_and_reduces_the_divergence() {
    let n = 22;
    let spacing = REFERENCE_SPACING;
    let config = DfsphConfig {
        max_pressure_iterations: 20,
        max_divergence_iterations: 20,
        error_threshold: 1e-6,
        gravity: [0.0, -9.81, 0.0],
        ..one_iteration_config(0.5, 1.3 * spacing)
    };
    // v = −rate (x − c) with 0.05 m/s at the faces of the block.
    let rate = 0.05 / (0.5 * (n as f64 - 1.0) * spacing);
    let mut solver = lattice_solver(config, n, spacing, compressing_field(n, spacing, rate));
    assert_eq!(solver.num_particles(), 10_648);
    let initial = divergence_error(&mut solver);
    assert!(initial > 0.05, "initial compressive divergence {initial}");

    let started = Instant::now();
    let (density_err, divergence_err) = solver.step(DT);
    let elapsed = started.elapsed();
    assert!(
        elapsed < Duration::from_secs(100),
        "one step of 10 648 particles took {elapsed:?}"
    );

    assert!(density_err.is_finite() && density_err >= 0.0);
    assert!(
        divergence_err < 0.1 * initial,
        "divergence error {initial} -> {divergence_err}"
    );
    assert_eq!(solver.stats.pressure_iterations.len(), 1);
    assert_eq!(solver.stats.divergence_iterations.len(), 1);
    assert!((1..=20).contains(&solver.stats.pressure_iterations[0]));
    assert!((1..=20).contains(&solver.stats.divergence_iterations[0]));
    let average_neighbours = solver.neighbor_lists.iter().map(Vec::len).sum::<usize>() as f64
        / solver.num_particles() as f64;
    assert!(
        (60.0..80.0).contains(&average_neighbours),
        "average neighbour count {average_neighbours}"
    );
    for p in &solver.particles {
        assert!(p.position.iter().all(|c| c.is_finite()));
        assert!(p.velocity.iter().all(|c| c.is_finite()));
        assert!(p.density.is_finite() && p.density > 0.0);
    }
}
