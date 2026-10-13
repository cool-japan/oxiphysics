// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the IISPH pressure solver (see the parent module).

use super::*;
use crate::kernel::CubicSplineKernel;

/// Build a larger uniform grid of particles for interior-density tests.
fn uniform_block(spacing: f64, n_per_side: usize, mass: f64) -> (Vec<Vec3>, Vec<Vec3>, Vec<f64>) {
    let mut positions = Vec::new();
    let mut velocities = Vec::new();
    let mut masses = Vec::new();
    for i in 0..n_per_side {
        for j in 0..n_per_side {
            for k in 0..n_per_side {
                positions.push(Vec3::new(
                    i as f64 * spacing,
                    j as f64 * spacing,
                    k as f64 * spacing,
                ));
                velocities.push(Vec3::zeros());
                masses.push(mass);
            }
        }
    }
    (positions, velocities, masses)
}

/// `Σ_{a,b,c ∈ Z} W(spacing |(a, b, c)|, h)` over an infinite cubic lattice.
fn lattice_kernel_sum(spacing: f64, h: f64) -> f64 {
    let kernel = CubicSplineKernel;
    let reach = (2.0 * h / spacing).ceil() as i32;
    let mut sum = 0.0;
    for a in -reach..=reach {
        for b in -reach..=reach {
            for c in -reach..=reach {
                let r = spacing * f64::from(a * a + b * b + c * c).sqrt();
                sum += kernel.w(r, h);
            }
        }
    }
    sum
}

/// Particle mass for which a deep interior particle of a cubic lattice of
/// spacing `rest_spacing` has the SPH density `rho0`.
fn calibrated_mass(rho0: f64, rest_spacing: f64, h: f64) -> f64 {
    rho0 / lattice_kernel_sum(rest_spacing, h)
}

/// Average and maximum density error after applying the accelerations over
/// `dt` (continuity equation evaluated with `v + dt a`), with the same
/// complementarity convention as the solver: particles at zero pressure
/// are only penalised for being too dense.  Independent of the solver's
/// internal displacement vectors.
fn predicted_density_error(
    solution: &IisphSolution,
    positions: &[Vec3],
    velocities: &[Vec3],
    masses: &[f64],
    neighbors: &[Vec<usize>],
    dt: f64,
    rho0: f64,
    h: f64,
) -> (f64, f64) {
    let kernel = CubicSplineKernel;
    let densities = IisphSolver::compute_densities(positions, masses, neighbors, &kernel, h);
    let new_velocity: Vec<Vec3> = velocities
        .iter()
        .zip(solution.accelerations.iter())
        .map(|(v, a)| v + a * dt)
        .collect();
    let mut sum = 0.0;
    let mut max = 0.0_f64;
    for i in 0..positions.len() {
        let mut drho = 0.0;
        for &j in &neighbors[i] {
            let rij = positions[i] - positions[j];
            let r = rij.norm();
            let grad = rij * (kernel.grad_w(r, h) / r);
            drho += masses[j] * (new_velocity[i] - new_velocity[j]).dot(&grad);
        }
        let deviation = (densities[i] + dt * drho - rho0) / rho0;
        let error = if solution.pressures[i] > 0.0 {
            deviation.abs()
        } else {
            deviation.max(0.0)
        };
        sum += error;
        max = max.max(error);
    }
    (sum / positions.len() as f64, max)
}

fn total_momentum_ratio(solution: &IisphSolution, masses: &[f64]) -> f64 {
    let net = solution
        .accelerations
        .iter()
        .zip(masses.iter())
        .fold(Vec3::zeros(), |acc, (a, &m)| acc + a * m);
    let scale: f64 = solution
        .accelerations
        .iter()
        .zip(masses.iter())
        .map(|(a, &m)| a.norm() * m)
        .sum();
    net.norm() / scale.max(1e-300)
}

#[test]
fn test_iisph_solver_creation() {
    let solver = IisphSolver::new(1000.0);
    assert!(
        (solver.density0 - 1000.0).abs() < 1e-10,
        "density0 should be 1000"
    );
    assert!(
        (solver.omega - 0.5).abs() < 1e-10,
        "default omega should be 0.5"
    );
    assert!(solver.max_iter >= 10, "max_iter should be at least 10");
    assert!(solver.tolerance > 0.0, "tolerance should be positive");
    assert!(solver.tolerance < 0.1, "tolerance should be small (<0.1)");
}

#[test]
fn test_iisph_density_computation() {
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let rho0 = 1000.0_f64;
    let mass = rho0 * spacing.powi(3);

    // Use a larger block so interior particles have enough neighbors.
    let (positions, _, masses) = uniform_block(spacing, 5, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let kernel = CubicSplineKernel;

    let densities = IisphSolver::compute_densities(&positions, &masses, &neighbors, &kernel, h);

    // At least some particles should have density close to rho0.
    let interior: Vec<f64> = densities.iter().copied().filter(|&d| d > 700.0).collect();
    assert!(
        !interior.is_empty(),
        "There should be particles with density > 700"
    );
    let avg: f64 = interior.iter().sum::<f64>() / interior.len() as f64;
    assert!(
        (avg - rho0).abs() < 300.0,
        "Average interior density = {avg:.1}, expected near {rho0}"
    );
}

#[test]
fn test_iisph_zero_pressure_stationary() {
    // A uniform, stationary fluid: total momentum change should be zero
    // (Newton's third law — pressure forces are symmetric).
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let rho0 = 1000.0_f64;
    let mass = rho0 * spacing.powi(3);
    let dt = 0.001_f64;

    let (positions, velocities, masses) = uniform_block(spacing, 5, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let kernel = CubicSplineKernel;

    let solver = IisphSolver::new(rho0);
    let (forces, avg_err) =
        solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

    // The avg density error must be finite.
    assert!(avg_err.is_finite(), "avg_density_error should be finite");

    // Total momentum change (sum of all pressure forces * mass) should be ~0.
    let net: Vec3 = forces
        .iter()
        .zip(masses.iter())
        .map(|(f, &m)| *f * m)
        .fold(Vec3::zeros(), |acc, v| acc + v);
    let net_norm = net.norm();
    let total_force: f64 = forces
        .iter()
        .zip(masses.iter())
        .map(|(f, &m)| f.norm() * m)
        .sum();
    let scale = total_force.max(1e-14);
    assert!(
        net_norm / scale < 1e-6,
        "Net pressure momentum should be ~0, got {net_norm:.4e} (scale {scale:.4e})"
    );
}

// --- Ihmsen IISPH: pressure response of compressed fluid ---

#[test]
fn test_compressed_block_reaches_incompressible_state() {
    // 10³ block compressed by 1 % in spacing (about 3 % in density): the
    // solver must remove the compression (average density error < 0.1 %)
    // with positive interior pressure and conserved momentum.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.99 * rest, 10, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);

    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 4000;
    solver.tolerance = 1e-4;
    let untouched = IisphSolver {
        max_iter: 0,
        ..solver.clone()
    }
    .solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert!(
        untouched.metrics.avg_error > 5e-3,
        "compression before the solve: {}",
        untouched.metrics.avg_error
    );
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );

    // the interior is genuinely over-dense before the solve
    let centre = 5 * 100 + 5 * 10 + 5;
    let rho_adv = solution.advected_densities[centre];
    assert!(
        rho_adv > 1.02 * rho0 && rho_adv < 1.05 * rho0,
        "interior advected density {rho_adv} should be ≈ 1.03 ρ₀"
    );

    let metrics = solution.metrics;
    assert!(metrics.converged, "not converged: {metrics:?}");
    assert!(
        metrics.avg_error < 1e-4,
        "average density error {} must be far below 0.1 %",
        metrics.avg_error
    );
    assert!(metrics.iterations > 5, "needs real work: {metrics:?}");

    // end-to-end check from the returned accelerations
    let (avg, max) = predicted_density_error(
        &solution,
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        rho0,
        h,
    );
    assert!(avg < 1e-4, "independent average error {avg}");
    assert!(
        (avg - metrics.avg_error).abs() < 1e-9,
        "{avg} vs {}",
        metrics.avg_error
    );
    assert!((max - metrics.max_error).abs() < 1e-9);

    assert!(
        solution.pressures[centre] > 0.0,
        "interior pressure must be positive"
    );
    assert_eq!(
        solution.pressures[0], 0.0,
        "free corner cannot carry pressure"
    );
    assert!(solution.pressures.iter().all(|&p| p >= 0.0));
    assert!(
        total_momentum_ratio(&solution, &masses) < 1e-9,
        "Σ m a must vanish"
    );

    // surface particles are pushed outward
    let block_centre =
        positions.iter().fold(Vec3::zeros(), |acc, p| acc + p) / positions.len() as f64;
    for corner in [0, 9, 90, 99, 900, 909, 990, 999] {
        let outward = (positions[corner] - block_centre).dot(&solution.accelerations[corner]);
        assert!(outward > 0.0, "corner {corner} must be pushed outward");
    }
}

#[test]
fn test_small_compressed_block_interior_pressure_and_outward_forces() {
    // 5³ block at spacing 0.04 whose masses fill the rest spacing 0.05
    // (h = 0.1): over-dense, so the interior carries pressure and the
    // surface is pushed away from the centre.
    let kernel = CubicSplineKernel;
    let h = 0.1;
    let rho0 = 1000.0;
    let mass = rho0 * 0.05_f64.powi(3);
    let dt = 1e-3;
    let (positions, velocities, masses) = uniform_block(0.04, 5, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 400;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );

    let centre = 2 * 25 + 2 * 5 + 2;
    assert!(
        solution.advected_densities[centre] > rho0,
        "the block must be over-dense"
    );
    assert!(
        solution.pressures[centre] > 0.0,
        "interior pressure must be positive"
    );
    let x_c = positions[centre];
    let radial_sum: f64 = positions
        .iter()
        .zip(solution.accelerations.iter())
        .map(|(x, a)| (x - x_c).dot(a))
        .sum();
    assert!(
        radial_sum > 0.0,
        "net outward push expected, got {radial_sum}"
    );
    for corner in [0, 4, 20, 24, 100, 104, 120, 124] {
        let outward = (positions[corner] - x_c).dot(&solution.accelerations[corner]);
        assert!(outward > 0.0, "corner {corner} must be pushed outward");
    }
    assert!(total_momentum_ratio(&solution, &masses) < 1e-9);
    assert!(
        solution.metrics.avg_error
            < solution
                .advected_densities
                .iter()
                .map(|r| (r - rho0).abs() / rho0)
                .sum::<f64>()
                / positions.len() as f64,
        "the solve must reduce the density error"
    );
}

#[test]
fn test_two_overdense_particles_repel_with_equal_and_opposite_forces() {
    let kernel = CubicSplineKernel;
    let h = 0.1;
    let rho0 = 1000.0;
    let mass = 5.0; // heavy: ρ_i = m (W(0) + W(d)) ≈ 3000 kg/m³
    let dt = 1e-3;
    let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.03, 0.0, 0.0)];
    let velocities = vec![Vec3::zeros(); 2];
    let masses = vec![mass; 2];
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let densities = IisphSolver::compute_densities(&positions, &masses, &neighbors, &kernel, h);
    assert!(densities[0] > 2.0 * rho0 && densities[1] > 2.0 * rho0);

    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 200;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert!(solution.pressures[0] > 0.0 && solution.pressures[1] > 0.0);
    let f0 = solution.accelerations[0];
    let f1 = solution.accelerations[1];
    assert!(f0.x < 0.0, "particle 0 pushed to −x, got {}", f0.x);
    assert!(f1.x > 0.0, "particle 1 pushed to +x, got {}", f1.x);
    assert!((f0.x + f1.x).abs() < 1e-9 * f0.x.abs());
    assert!(f0.y.abs() < 1e-12 * f0.x.abs() && f0.z.abs() < 1e-12 * f0.x.abs());
}

#[test]
fn test_block_at_rest_density_needs_no_pressure() {
    // Calibrated rest lattice: deep interior particles have ρ = ρ₀ exactly,
    // the sparser surface layers cannot be corrected without suction.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(rest, 9, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let solver = IisphSolver::new(rho0);
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        1e-3,
        &kernel,
        h,
        None,
    );
    let centre = 4 * 81 + 4 * 9 + 4;
    assert!((solution.densities[centre] - rho0).abs() < 1e-9 * rho0);
    assert!(solution.metrics.converged);
    assert_eq!(solution.metrics.iterations, 0);
    assert!(solution.pressures.iter().all(|&p| p == 0.0));
    assert!(solution.accelerations.iter().all(|a| a.norm() == 0.0));
}

#[test]
fn test_converging_flow_is_projected_to_divergence_free_compression() {
    // Inward velocity field v = −c (x − x_c) squeezes a rest block by about
    // 3 c dt in density; the pressure acceleration must cancel it.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let c = 10.0;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, _, masses) = uniform_block(rest, 8, mass);
    let centre_x = positions.iter().fold(Vec3::zeros(), |acc, p| acc + p) / positions.len() as f64;
    let velocities: Vec<Vec3> = positions.iter().map(|x| (x - centre_x) * (-c)).collect();
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 3000;
    solver.tolerance = 1e-5;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    let interior = 3 * 64 + 3 * 8 + 3;
    let compression = solution.advected_densities[interior] / rho0 - 1.0;
    assert!(
        compression > 0.02 && compression < 0.04,
        "converging flow compresses by ≈ 3 c dt = 3 %, got {compression}"
    );
    assert!(solution.metrics.converged, "{:?}", solution.metrics);
    let (avg, _) = predicted_density_error(
        &solution,
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        rho0,
        h,
    );
    assert!(avg < 1e-5, "post-projection density error {avg}");
    assert!(solution.pressures[interior] > 0.0);
    assert!(total_momentum_ratio(&solution, &masses) < 1e-9);
}

#[test]
fn test_pressure_scales_inversely_with_dt_squared() {
    // The system matrix is proportional to Δt², so the pressure at every
    // Jacobi iteration scales with 1/Δt² and the iteration counts agree.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.97 * rest, 6, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 1500;
    let coarse = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        2e-3,
        &kernel,
        h,
        None,
    );
    let fine = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        1e-3,
        &kernel,
        h,
        None,
    );
    assert!(coarse.metrics.converged && fine.metrics.converged);
    assert_eq!(coarse.metrics.iterations, fine.metrics.iterations);
    let p_max = coarse.pressures.iter().cloned().fold(0.0, f64::max);
    assert!(p_max > 0.0);
    for (pc, pf) in coarse.pressures.iter().zip(fine.pressures.iter()) {
        assert!(
            (pf - 4.0 * pc).abs() < 1e-9 * 4.0 * p_max,
            "{pf} vs 4 × {pc}"
        );
    }
}

#[test]
fn test_wide_kernel_relaxation_is_lowered_and_converges() {
    // Support of four particle spacings (h = 2 s, as in the other SPH
    // tests): λ_max(D⁻¹A) ≈ 10, so the plain ω = 0.5 Jacobi iteration
    // diverges.  The safeguarded iteration lowers ω and converges.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 2.0 * rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.97 * rest, 8, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 5000;
    solver.tolerance = 1e-4;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert!(
        solution.omega < 0.25,
        "ω must be lowered for a wide kernel, got {}",
        solution.omega
    );
    assert!(solution.metrics.converged, "{:?}", solution.metrics);
    let (avg, _) = predicted_density_error(
        &solution,
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        rho0,
        h,
    );
    assert!(avg < 1e-4, "independent average error {avg}");
    assert!(total_momentum_ratio(&solution, &masses) < 1e-9);
}

#[test]
fn test_narrow_kernel_keeps_configured_relaxation() {
    // Support of two particle spacings: λ_max ≈ 3, ω = 0.5 contracts and is
    // kept (at most mildly reduced); a smaller configured ω is never raised.
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.97 * rest, 8, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 2000;
    solver.tolerance = 1e-4;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert!(
        solution.omega > 0.4 && solution.omega <= 0.5,
        "ω = {}",
        solution.omega
    );
    assert!(solution.metrics.converged, "{:?}", solution.metrics);

    solver.omega = 0.2;
    let slow = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert_eq!(slow.omega, 0.2);
    assert!(slow.metrics.converged, "{:?}", slow.metrics);
}

#[test]
fn test_warm_start_from_converged_pressure_needs_no_iterations() {
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.97 * rest, 6, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 3000;
    solver.tolerance = 1e-6;
    let cold = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert!(cold.metrics.converged);
    assert!(cold.metrics.iterations > 20);

    let warm = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        Some(&cold.pressures),
    );
    assert!(warm.metrics.converged);
    assert_eq!(warm.metrics.iterations, 0);
    for (a, b) in warm.accelerations.iter().zip(cold.accelerations.iter()) {
        assert!((a - b).norm() <= 1e-12 * (1.0 + b.norm()));
    }

    // the convenience wrapper with warm start agrees with the core
    let (forces, avg_err) = solver.solve_pressure_warm_start(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        &cold.pressures,
    );
    assert_eq!(forces.len(), positions.len());
    assert!((avg_err - warm.metrics.avg_error).abs() < 1e-15);

    // a warm start with the wrong length is ignored (zero initial guess)
    let ignored = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        Some(&cold.pressures[..3]),
    );
    assert_eq!(ignored.metrics.iterations, cold.metrics.iterations);
}

#[test]
fn test_solve_rejects_inconsistent_input() {
    let kernel = CubicSplineKernel;
    let solver = IisphSolver::new(1000.0);
    let positions = vec![Vec3::zeros(), Vec3::new(0.05, 0.0, 0.0)];
    let velocities = vec![Vec3::zeros(); 2];
    let masses = vec![1.0; 2];
    let neighbors = vec![vec![1], vec![0]];

    let short_velocities = vec![Vec3::zeros(); 1];
    let rejected = solver.solve(
        &positions,
        &short_velocities,
        &masses,
        &neighbors,
        1e-3,
        &kernel,
        0.1,
        None,
    );
    assert!(!rejected.metrics.converged);
    assert!(rejected.metrics.avg_error.is_infinite());
    assert_eq!(rejected.accelerations.len(), 2);
    assert!(rejected.accelerations.iter().all(|a| a.norm() == 0.0));

    let bad_neighbors = vec![vec![7], vec![0]];
    let rejected = solver.solve(
        &positions,
        &velocities,
        &masses,
        &bad_neighbors,
        1e-3,
        &kernel,
        0.1,
        None,
    );
    assert!(!rejected.metrics.converged);
    assert!(rejected.pressures.iter().all(|&p| p == 0.0));

    let (forces, err) = solver.solve_pressure(
        &positions,
        &velocities,
        &masses[..1],
        &neighbors,
        1e-3,
        &kernel,
        0.1,
    );
    assert_eq!(forces.len(), 2);
    assert!(err.is_infinite());
}

#[test]
fn test_solve_empty_and_isolated_inputs() {
    let kernel = CubicSplineKernel;
    let solver = IisphSolver::new(1000.0);
    let empty = solver.solve(&[], &[], &[], &[], 1e-3, &kernel, 0.1, None);
    assert!(empty.metrics.converged && empty.pressures.is_empty());
    assert_eq!(empty.metrics.avg_error, 0.0);

    // a single heavy particle without neighbours cannot be relieved
    let alone = solver.solve(
        &[Vec3::zeros()],
        &[Vec3::zeros()],
        &[1e3],
        &[Vec::new()],
        1e-3,
        &kernel,
        0.1,
        None,
    );
    assert_eq!(alone.diagonal[0], 0.0);
    assert_eq!(alone.pressures[0], 0.0);
    assert!(alone.metrics.converged);
    assert_eq!(alone.accelerations[0].norm(), 0.0);
}

#[test]
fn test_max_iter_zero_reports_initial_error() {
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.97 * rest, 10, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 0;
    let (forces, metrics) = solver.solve_pressure_with_metrics(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        1e-3,
        &kernel,
        h,
    );
    assert_eq!(metrics.iterations, 0);
    assert!(!metrics.converged);
    assert!(
        metrics.avg_error > 0.01,
        "compressed block starts far from ρ₀"
    );
    assert!(forces.iter().all(|f| f.norm() == 0.0));
}

// --- Diagonal coefficient ---

#[test]
fn test_diagonal_coefficient_is_negative_and_matches_closed_form() {
    // Two particles at distance d: ∇W_12 = W'(d) x̂ (|∇W| = |W'|), so
    // a_11 = −Δt²/ρ_1² (m_2² + m_1 m_2) |W'(d)|².
    let kernel = CubicSplineKernel;
    let h = 0.1;
    let d = 0.07;
    let dt = 2e-3;
    let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(d, 0.0, 0.0)];
    let masses = vec![0.3, 0.2];
    let densities = vec![900.0, 1100.0];
    let a11 =
        compute_diagonal_coefficient(0, &positions, &masses, &densities, &[1], &kernel, h, dt);
    let dw = kernel.grad_w(d, h);
    let expected = -dt * dt / (900.0 * 900.0) * (0.2 * 0.2 + 0.3 * 0.2) * dw * dw;
    assert!(a11 < 0.0);
    assert!(
        (a11 - expected).abs() < 1e-12 * expected.abs(),
        "a_11 = {a11}, closed form {expected}"
    );
    // particle 2 sees particle 1 in the opposite direction with its own mass and density
    let a22 =
        compute_diagonal_coefficient(1, &positions, &masses, &densities, &[0], &kernel, h, dt);
    let expected2 = -dt * dt / (1100.0 * 1100.0) * (0.3 * 0.3 + 0.2 * 0.3) * dw * dw;
    assert!((a22 - expected2).abs() < 1e-12 * expected2.abs());
    // quadratic in dt
    let half = compute_diagonal_coefficient(
        0,
        &positions,
        &masses,
        &densities,
        &[1],
        &kernel,
        h,
        0.5 * dt,
    );
    assert!((half - 0.25 * a11).abs() < 1e-15 * a11.abs());
}

#[test]
fn test_diagonal_coefficient_isolated_and_out_of_range() {
    let kernel = CubicSplineKernel;
    let positions = vec![Vec3::zeros(), Vec3::new(0.05, 0.0, 0.0)];
    let masses = vec![1.0, 1.0];
    let densities = vec![1000.0, 1000.0];
    assert_eq!(
        compute_diagonal_coefficient(0, &positions, &masses, &densities, &[], &kernel, 0.1, 1e-3),
        0.0
    );
    assert_eq!(
        compute_diagonal_coefficient(5, &positions, &masses, &densities, &[1], &kernel, 0.1, 1e-3),
        0.0
    );
    // an out-of-range neighbour index is skipped
    let skipped =
        compute_diagonal_coefficient(0, &positions, &masses, &densities, &[9], &kernel, 0.1, 1e-3);
    assert_eq!(skipped, 0.0);
}

#[test]
fn test_solver_diagonal_matches_public_helper() {
    let kernel = CubicSplineKernel;
    let h = 0.1;
    let (positions, velocities, masses) = uniform_block(0.05, 4, 0.125);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let solver = IisphSolver::new(1000.0);
    let dt = 1e-3;
    let solution = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    for (i, list) in neighbors.iter().enumerate() {
        let expected = compute_diagonal_coefficient(
            i,
            &positions,
            &masses,
            &solution.densities,
            list,
            &kernel,
            h,
            dt,
        );
        assert!(solution.diagonal[i] < 0.0);
        assert!(
            (solution.diagonal[i] - expected).abs() < 1e-12 * expected.abs(),
            "particle {i}: {} vs {expected}",
            solution.diagonal[i]
        );
    }
}

// --- Source term, Jacobi update and relaxation helpers ---

#[test]
fn test_relaxed_jacobi_update_compression_gives_positive_pressure() {
    // Ihmsen's a_ii is negative, the source of a compressed particle is
    // negative: p = ω s / a_ii = 0.5 · (−30) / (−2) = 7.5.
    let p = relaxed_jacobi_update(0.0, -30.0, 0.0, -2.0, 0.5);
    assert!((p - 7.5).abs() < 1e-12, "got {p}");
    // under relaxation keeps a share of the old pressure: 0.5·4 + 0.5·15
    let p = relaxed_jacobi_update(4.0, -30.0, 0.0, -2.0, 0.5);
    assert!((p - 9.5).abs() < 1e-12, "got {p}");
    // the off-diagonal part enters with a minus sign: (−30 − 10)/(−2) = 20
    let p = relaxed_jacobi_update(0.0, -30.0, 10.0, -2.0, 1.0);
    assert!((p - 20.0).abs() < 1e-12, "got {p}");
    let p = relaxed_jacobi_update(0.0, -30.0, -10.0, -2.0, 1.0);
    assert!((p - 10.0).abs() < 1e-12, "got {p}");
}

#[test]
fn test_relaxed_jacobi_update_expanded_particle_clamped_to_zero() {
    // an under-dense particle (source > 0) would need suction
    let p = relaxed_jacobi_update(0.0, 30.0, 0.0, -2.0, 0.5);
    assert_eq!(p, 0.0);
    // negligible diagonal: isolated particle
    assert_eq!(relaxed_jacobi_update(3.0, -30.0, 0.0, 0.0, 0.5), 0.0);
}

#[test]
fn test_estimate_optimal_omega_uses_diagonal_magnitudes() {
    // The IISPH diagonal is negative; the heuristic must not collapse to the default.
    let negative = vec![-0.1, -0.2, -0.3, -0.15, -0.25];
    let positive: Vec<f64> = negative.iter().map(|a| -a).collect();
    let omega = estimate_optimal_omega(&negative, 5);
    assert!((omega - estimate_optimal_omega(&positive, 5)).abs() < 1e-15);
    assert!((0.3..=0.9).contains(&omega), "omega out of range: {omega}");
    // uniform diagonal: best conditioned, largest factor 0.9
    assert!((estimate_optimal_omega(&[-0.2; 4], 4) - 0.9).abs() < 1e-12);
    // strongly non-uniform: ratio 1/100 → 0.5 + 0.4 · 0.01 = 0.504
    assert!((estimate_optimal_omega(&[-1.0, -100.0], 2) - 0.504).abs() < 1e-12);
    // isolated particles are ignored
    assert!((estimate_optimal_omega(&[0.0, -0.2, -0.2], 3) - 0.9).abs() < 1e-12);
    assert!((estimate_optimal_omega(&[0.0, 0.0], 2) - 0.5).abs() < 1e-12);
}

#[test]
fn test_solve_pressure_with_metrics_matches_core() {
    let kernel = CubicSplineKernel;
    let rest = 0.05;
    let h = 1.3 * rest;
    let rho0 = 1000.0;
    let dt = 1e-3;
    let mass = calibrated_mass(rho0, rest, h);
    let (positions, velocities, masses) = uniform_block(0.98 * rest, 5, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let mut solver = IisphSolver::new(rho0);
    solver.max_iter = 1000;
    let (forces, metrics) = solver.solve_pressure_with_metrics(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
    );
    let core = solver.solve(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
        None,
    );
    assert_eq!(forces.len(), positions.len());
    assert!(metrics.avg_error.is_finite());
    assert!(metrics.max_error >= metrics.avg_error);
    assert!(metrics.iterations > 0);
    assert_eq!(metrics.iterations, core.metrics.iterations);
    assert_eq!(metrics.converged, core.metrics.converged);
    for (a, b) in forces.iter().zip(core.accelerations.iter()) {
        assert_eq!(a, b);
    }
    let (forces2, avg2) =
        solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);
    assert_eq!(avg2, metrics.avg_error);
    assert_eq!(forces2.len(), forces.len());
}

// --- Density error metrics tests ---

#[test]
fn test_density_error_metrics_zero() {
    let m = DensityErrorMetrics::zero();
    assert!(m.avg_error.abs() < 1e-15);
    assert!(m.max_error.abs() < 1e-15);
    assert_eq!(m.iterations, 0);
    assert!(m.converged);
}

// --- IISPH timestep adapter tests ---

#[test]
fn test_timestep_adapter_shrink_on_poor_convergence() {
    let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
    let metrics = DensityErrorMetrics {
        avg_error: 0.1,
        max_error: 0.5,
        iterations: 100,
        converged: false,
    };
    let new_dt = adapter.adapt(&metrics);
    assert!(
        new_dt < 0.001,
        "dt should shrink on poor convergence, got {new_dt}"
    );
}

#[test]
fn test_timestep_adapter_grow_on_good_convergence() {
    let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
    adapter.growth_threshold = 2;
    let metrics = DensityErrorMetrics {
        avg_error: 1e-5,
        max_error: 1e-4,
        iterations: 5,
        converged: true,
    };
    // Need multiple good steps
    adapter.adapt(&metrics);
    let new_dt = adapter.adapt(&metrics);
    assert!(
        new_dt > 0.001,
        "dt should grow after consecutive good steps, got {new_dt}"
    );
}

#[test]
fn test_timestep_adapter_reset() {
    let mut adapter = IisphTimestepAdapter::new(0.001, 1e-6, 0.01);
    adapter.good_steps = 10;
    adapter.reset(0.005);
    assert!((adapter.dt - 0.005).abs() < 1e-15);
    assert_eq!(adapter.good_steps, 0);
}

#[test]
fn test_timestep_adapter_respects_bounds() {
    let mut adapter = IisphTimestepAdapter::new(0.001, 1e-4, 0.01);
    let bad = DensityErrorMetrics {
        avg_error: 1.0,
        max_error: 5.0,
        iterations: 100,
        converged: false,
    };
    // Repeatedly shrink
    for _ in 0..100 {
        adapter.adapt(&bad);
    }
    assert!(adapter.dt >= 1e-4, "dt must not go below min_dt");
}

// --- Diagonal coefficient test ---

// --- Source term test ---

#[test]
fn test_compute_source_term() {
    let source = compute_source_term(1050.0, 1000.0);
    assert!(
        (source - (-50.0)).abs() < 1e-12,
        "source should be rho0 - rho_adv"
    );
}

// --- Relaxed Jacobi update test ---

// --- Optimal omega test ---

#[test]
fn test_estimate_optimal_omega_range() {
    let a_ii = vec![0.1, 0.2, 0.3, 0.15, 0.25];
    let omega = estimate_optimal_omega(&a_ii, 5);
    assert!((0.3..=0.9).contains(&omega), "omega out of range: {omega}");
}

#[test]
fn test_estimate_optimal_omega_empty() {
    let omega = estimate_optimal_omega(&[], 0);
    assert!((omega - 0.5).abs() < 1e-12, "Default omega should be 0.5");
}

// --- Solve with metrics test ---

// --- IISPH-FLIP hybrid tests ---

#[test]
fn test_iisph_flip_hybrid_blending() {
    // blend = 0 → pure IISPH (no FLIP), blend = 1 → pure FLIP.
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let rho0 = 1000.0_f64;
    let mass = rho0 * spacing.powi(3);
    let dt = 0.001_f64;

    let (positions, velocities, masses) = uniform_block(spacing, 3, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let kernel = CubicSplineKernel;
    let solver = IisphSolver::new(rho0);

    let (forces_pure, _) =
        solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

    // Apply FLIP blending (blend=0 → identical to pure IISPH).
    let blended = iisph_flip_blend(&forces_pure, &velocities, &velocities, 0.0, masses.len());
    for (i, (f, b)) in forces_pure.iter().zip(blended.iter()).enumerate() {
        let diff = (f - b).norm();
        assert!(
            diff < 1e-12,
            "blend=0 should be identical to IISPH forces at particle {i}: diff={diff}"
        );
    }
}

#[test]
fn test_iisph_flip_hybrid_full_flip() {
    // blend = 1 → velocities = old_velocities + dt * forces (FLIP update).
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let rho0 = 1000.0_f64;
    let mass = rho0 * spacing.powi(3);
    let dt = 0.001_f64;

    let (positions, velocities, masses) = uniform_block(spacing, 3, mass);
    let n = positions.len();
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let kernel = CubicSplineKernel;
    let solver = IisphSolver::new(rho0);

    let (forces, _) =
        solver.solve_pressure(&positions, &velocities, &masses, &neighbors, dt, &kernel, h);

    // With blend=1 we get FLIP: v_new = v_old + dt * force.
    // The blended function returns modified "effective forces" or velocities.
    let blended = iisph_flip_blend(&forces, &velocities, &velocities, 1.0, n);

    // blend=1 returns forces unchanged (pure IISPH pressure forces are kept).
    assert_eq!(
        blended.len(),
        n,
        "blended should have same length as forces"
    );
    for f in &blended {
        assert!(f.norm().is_finite(), "blended forces should be finite");
    }
}

// --- IISPH multiple iterations / error monitoring ---

#[test]
fn test_iisph_error_monitor_tracks_iterations() {
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let rho0 = 1000.0_f64;
    let mass = rho0 * spacing.powi(3);
    let dt = 0.001_f64;

    let (positions, velocities, masses) = uniform_block(spacing, 4, mass);
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
    let kernel = CubicSplineKernel;

    let solver = IisphSolver::new(rho0);
    let mut monitor = IisphErrorMonitor::new();
    let (_, metrics) = solver.solve_pressure_with_metrics(
        &positions,
        &velocities,
        &masses,
        &neighbors,
        dt,
        &kernel,
        h,
    );
    monitor.record(metrics.avg_error, metrics.iterations);

    assert_eq!(monitor.num_steps(), 1);
    assert!(monitor.last_error().is_finite());
}

#[test]
fn test_iisph_error_monitor_multiple_steps() {
    let mut monitor = IisphErrorMonitor::new();
    monitor.record(0.1, 10);
    monitor.record(0.05, 8);
    monitor.record(0.02, 5);
    assert_eq!(monitor.num_steps(), 3);
    assert!((monitor.last_error() - 0.02).abs() < 1e-14);
    assert!(monitor.avg_iterations() > 0.0);
}

// --- IISPH acceleration data structures ---

#[test]
fn test_iisph_neighborhood_cache_lookup() {
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let (positions, _, _) = uniform_block(spacing, 3, 1.0);
    let cache = IisphNeighborhoodCache::build(&positions, 2.0 * h);

    // Every particle should have at least one neighbor.
    for i in 0..positions.len() {
        let neighbors = cache.neighbors(i);
        assert!(!neighbors.is_empty(), "Particle {i} should have neighbors");
    }
}

#[test]
fn test_iisph_neighborhood_cache_no_self_neighbor() {
    let spacing = 0.05_f64;
    let h = 0.1_f64;
    let (positions, _, _) = uniform_block(spacing, 3, 1.0);
    let cache = IisphNeighborhoodCache::build(&positions, 2.0 * h);

    for i in 0..positions.len() {
        for &j in cache.neighbors(i) {
            assert_ne!(i, j, "Particle should not neighbor itself");
        }
    }
}
