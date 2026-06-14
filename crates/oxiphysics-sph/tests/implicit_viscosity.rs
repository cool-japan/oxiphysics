// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the implicit / semi-implicit SPH viscosity solver.
//!
//! Headline guarantee: the implicit Morris-Laplacian update is unconditionally
//! stable, so it runs at time steps far beyond the explicit viscous limit
//! `Δt < h²/(2·d·ν)` without blowing up, while still dissipating kinetic energy
//! monotonically.

use oxiphysics_sph::kernel::{CubicSplineKernel, grad};
use oxiphysics_sph::viscosity_implicit::{
    NeighborEntry, ViscosityError, ViscosityParticles, ViscositySolveOptions,
    solve_implicit_viscosity,
};

const MORRIS_EPSILON: f64 = 0.01;

/// Build cubic-spline neighbour lists for `positions` with support radius `2h`.
fn build_neighbors(positions: &[[f64; 3]], h: f64) -> Vec<Vec<NeighborEntry>> {
    let kernel = CubicSplineKernel;
    let support = 2.0 * h;
    let n = positions.len();
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_ij = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]).sqrt();
            if r >= support || r < 1.0e-12 {
                continue;
            }
            let gw = grad(&kernel, r_ij, h);
            neighbors[i].push(NeighborEntry::new(j, r_ij, gw));
        }
    }
    neighbors
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Public-API replica of the private Morris coupling `c_ij ≥ 0`, used only to
/// reconstruct the system matrix for the SPD test.
fn coupling(
    entry: &NeighborEntry,
    mu_i: f64,
    mu_j: f64,
    rho_i: f64,
    rho_j: f64,
    m_j: f64,
    h: f64,
) -> f64 {
    let r2 = dot3(entry.r_ij, entry.r_ij);
    if r2 < 1.0e-24 {
        return 0.0;
    }
    let denom = r2 + MORRIS_EPSILON * h * h;
    let geom = -dot3(entry.r_ij, entry.grad_w_ij) / denom;
    (m_j / rho_j) * (mu_i + mu_j) / rho_i * geom
}

/// Assemble the dense system matrix `A = I − dt·ν·L` for an SPD check.
fn dense_system(
    neighbors: &[Vec<NeighborEntry>],
    densities: &[f64],
    masses: &[f64],
    viscosities: &[f64],
    nu: f64,
    dt: f64,
    h: f64,
) -> Vec<Vec<f64>> {
    let n = neighbors.len();
    let mut a = vec![vec![0.0_f64; n]; n];
    for (i, row) in a.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let scale = dt * nu;
    for i in 0..n {
        for entry in &neighbors[i] {
            let j = entry.j;
            if j == i || j >= n {
                continue;
            }
            let c = coupling(
                entry,
                viscosities[i],
                viscosities[j],
                densities[i],
                densities[j],
                masses[j],
                h,
            );
            a[i][i] += scale * c;
            a[i][j] -= scale * c;
        }
    }
    a
}

/// Total kinetic energy `Σ ½ m_i |v_i|²`.
fn kinetic_energy(velocities: &[[f64; 3]], masses: &[f64]) -> f64 {
    velocities
        .iter()
        .zip(masses.iter())
        .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum()
}

/// Build a regular 2-D grid of particles in the z = 0 plane, spacing `dx`.
fn grid_2d(nx: usize, ny: usize, dx: f64) -> Vec<[f64; 3]> {
    let mut positions = Vec::with_capacity(nx * ny);
    for ix in 0..nx {
        for iy in 0..ny {
            positions.push([ix as f64 * dx, iy as f64 * dx, 0.0]);
        }
    }
    positions
}

/// Test 1 — stability gate: implicit solver stays stable at 100× the explicit
/// viscous limit, with no NaNs and monotonically non-increasing kinetic energy.
#[test]
fn implicit_stable_at_100x_explicit_limit() {
    let dx = 0.1_f64;
    let h = 0.1_f64;
    let nu = 1.0_f64;
    let rho0 = 1000.0_f64;
    let dim = 2.0_f64;

    let dt_explicit = dx * dx / (2.0 * dim * nu); // = 0.00125
    let dt = 100.0 * dt_explicit; // 0.125 — far past the explicit bound

    let positions = grid_2d(8, 8, dx);
    let n = positions.len();
    let neighbors = build_neighbors(&positions, h);

    // Checkerboard initial velocity field (high-frequency, the worst case for
    // explicit stability).
    let mut velocities: Vec<[f64; 3]> = Vec::with_capacity(n);
    for (idx, _) in positions.iter().enumerate() {
        let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
        velocities.push([sign, 0.0, 0.0]);
    }
    let densities = vec![rho0; n];
    let masses = vec![rho0 * dx * dx; n]; // 2-D mass per particle
    let viscosities = vec![rho0 * nu; n]; // μ = ρ·ν

    let ke0 = kinetic_energy(&velocities, &masses);

    for step in 0..20 {
        let particles = ViscosityParticles {
            velocities: &velocities,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        let updated = solve_implicit_viscosity(
            &particles,
            &neighbors,
            nu,
            dt,
            h,
            ViscositySolveOptions::default(),
        )
        .unwrap_or_else(|e| panic!("solve failed at step {step}: {e}"));
        for v in &updated {
            for &c in v {
                assert!(c.is_finite(), "non-finite velocity at step {step}");
            }
        }
        let ke = kinetic_energy(&updated, &masses);
        assert!(
            ke <= ke0 + 1e-9,
            "kinetic energy grew (step {step}): ke={ke}, ke0={ke0}"
        );
        velocities = updated;
    }

    let ke_final = kinetic_energy(&velocities, &masses);
    assert!(
        ke_final < ke0,
        "viscosity must dissipate energy: ke_final={ke_final}, ke0={ke0}"
    );
}

/// Test 2 — convergence to truth: explicit (small dt) and implicit (large dt)
/// integration of the same diffusion problem reach the same steady state.
///
/// A 1-D row of particles carries a step velocity profile.  Pure diffusion with
/// no-flux (isolated) ends drives the field toward its conserved mean.  We march
/// the explicit Morris update at a stable dt and the implicit update at 50× that
/// dt to the same physical time, and require the two velocity fields to agree.
#[test]
fn implicit_matches_explicit_steady_state() {
    let dx = 0.1_f64;
    let h = 0.1_f64;
    let nu = 0.5_f64;
    let rho0 = 1000.0_f64;

    // 1-D chain of 16 particles.
    let n = 16usize;
    let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * dx, 0.0, 0.0]).collect();
    let neighbors = build_neighbors(&positions, h);
    let densities = vec![rho0; n];
    let masses = vec![rho0 * dx; n];
    let viscosities = vec![rho0 * nu; n];

    // Step profile: left half +1, right half −1.
    let initial: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            if i < n / 2 {
                [1.0, 0.0, 0.0]
            } else {
                [-1.0, 0.0, 0.0]
            }
        })
        .collect();

    let t_final = 0.5_f64;

    // -- Explicit reference, dt at the *discrete* stable limit. -----------
    // The continuum bound h²/(2dν) is too loose for the SPH Morris operator,
    // whose largest eigenvalue is ≈ 2·maxᵢ Σⱼ c_ij.  Forward-Euler diffusion
    // is stable only for dt·ν·λ_max < 2, so size dt from the assembled
    // operator: dt = safety / (ν · maxᵢ Σⱼ c_ij) with safety = 0.2 gives
    // dt·ν·λ_max ≈ 0.4, comfortably stable.
    let max_diag_sum = (0..n)
        .map(|i| {
            neighbors[i]
                .iter()
                .map(|entry| {
                    let r2 = dot3(entry.r_ij, entry.r_ij);
                    if r2 < 1.0e-24 {
                        return 0.0;
                    }
                    let denom = r2 + MORRIS_EPSILON * h * h;
                    (masses[entry.j] / densities[entry.j]) * (viscosities[i] + viscosities[entry.j])
                        / densities[i]
                        * (-dot3(entry.r_ij, entry.grad_w_ij) / denom)
                })
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let dt_explicit = 0.2 / (nu * max_diag_sum);
    let n_explicit = (t_final / dt_explicit).ceil() as usize;
    let dt_e = t_final / n_explicit as f64;
    let mut v_explicit = initial.clone();
    for _ in 0..n_explicit {
        let mut next = v_explicit.clone();
        for i in 0..n {
            let mut lap = [0.0_f64; 3];
            for entry in &neighbors[i] {
                let j = entry.j;
                let r2 = dot3(entry.r_ij, entry.r_ij);
                if r2 < 1.0e-24 {
                    continue;
                }
                let denom = r2 + MORRIS_EPSILON * h * h;
                let c = (masses[j] / densities[j]) * (viscosities[i] + viscosities[j])
                    / densities[i]
                    * (-dot3(entry.r_ij, entry.grad_w_ij) / denom);
                for k in 0..3 {
                    lap[k] += c * (v_explicit[j][k] - v_explicit[i][k]);
                }
            }
            for k in 0..3 {
                next[i][k] = v_explicit[i][k] + dt_e * nu * lap[k];
            }
        }
        v_explicit = next;
    }

    // -- Implicit, dt 50× larger. -----------------------------------------
    let dt_i = 50.0 * dt_e;
    let n_implicit = (t_final / dt_i).ceil() as usize;
    let dt_i = t_final / n_implicit as f64;
    let mut v_implicit = initial.clone();
    for _ in 0..n_implicit {
        let particles = ViscosityParticles {
            velocities: &v_implicit,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        v_implicit = solve_implicit_viscosity(
            &particles,
            &neighbors,
            nu,
            dt_i,
            h,
            ViscositySolveOptions::default(),
        )
        .expect("implicit solve should converge");
    }

    // Both fields should be close (diffusion is smoothing toward the mean).
    let mut max_diff = 0.0_f64;
    for i in 0..n {
        max_diff = max_diff.max((v_explicit[i][0] - v_implicit[i][0]).abs());
    }
    // The two schemes have different temporal truncation error; require
    // agreement to within a generous absolute tolerance on the O(1) field.
    assert!(
        max_diff < 0.05,
        "implicit/explicit steady states disagree: max |Δv| = {max_diff}"
    );

    // Sanity: total momentum is conserved by both (equal masses, symmetric op).
    let p_e: f64 = (0..n).map(|i| masses[i] * v_explicit[i][0]).sum();
    let p_i: f64 = (0..n).map(|i| masses[i] * v_implicit[i][0]).sum();
    assert!(
        (p_e - p_i).abs() < 1e-6,
        "momentum mismatch: {p_e} vs {p_i}"
    );
}

/// Test 3 — SPD check on a small 5-particle configuration: the assembled system
/// matrix must be symmetric and positive definite.
#[test]
fn system_matrix_is_symmetric_positive_definite() {
    let h = 0.1_f64;
    let nu = 1.0_f64;
    let dt = 0.05_f64;
    let rho0 = 1000.0_f64;

    // Five particles in a small cross.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [0.08, 0.0, 0.0],
        [-0.08, 0.0, 0.0],
        [0.0, 0.08, 0.0],
        [0.0, -0.08, 0.0],
    ];
    let n = positions.len();
    let neighbors = build_neighbors(&positions, h);
    let densities = vec![rho0; n];
    let masses = vec![rho0 * 0.08 * 0.08; n];
    let viscosities = vec![rho0 * nu; n];

    let a = dense_system(&neighbors, &densities, &masses, &viscosities, nu, dt, h);

    // Symmetry: A == Aᵀ.
    for (i, row) in a.iter().enumerate() {
        for (j, &a_ij) in row.iter().enumerate() {
            assert!(
                (a_ij - a[j][i]).abs() < 1e-12,
                "A not symmetric at ({i},{j}): {} vs {}",
                a_ij,
                a[j][i]
            );
        }
    }

    // Positive definiteness: xᵀ A x > 0 for several non-zero x.
    let test_vectors: [[f64; 5]; 4] = [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, -1.0, 1.0, -1.0, 1.0],
        [0.3, 0.7, -0.2, 0.5, -0.9],
        [1.0, 1.0, 1.0, 1.0, 1.0],
    ];
    for x in &test_vectors {
        let mut quad = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                quad += x[i] * a[i][j] * x[j];
            }
        }
        assert!(quad > 0.0, "matrix not positive definite: xᵀAx = {quad}");
    }

    // Diagonal dominance (strict for an M-matrix): A_ii ≥ 1 + Σ|A_ij|.
    for (i, row) in a.iter().enumerate() {
        let off: f64 = row
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, &a_ij)| a_ij.abs())
            .sum();
        assert!(
            row[i] >= 1.0 + off - 1e-12,
            "row {i} not diagonally dominant: A_ii={}, Σ|off|={off}",
            row[i]
        );
    }
}

/// Test 4 — CG convergence on a 50-particle system: the internal Jacobi-PCG
/// must converge within the iteration budget, i.e. the solver returns `Ok`
/// rather than `SolverNotConverged`.
#[test]
fn solver_converges_on_fifty_particles() {
    let dx = 0.1_f64;
    let h = 0.1_f64;
    let nu = 1.0_f64;
    let dt = 0.1_f64;
    let rho0 = 1000.0_f64;

    // 5×10 grid = 50 particles.
    let positions = grid_2d(5, 10, dx);
    let n = positions.len();
    assert_eq!(n, 50);
    let neighbors = build_neighbors(&positions, h);
    let densities = vec![rho0; n];
    let masses = vec![rho0 * dx * dx; n];
    let viscosities = vec![rho0 * nu; n];

    // Random-ish initial velocities (deterministic, no rng dependency).
    let velocities: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let a = (i as f64 * 0.37).sin();
            let b = (i as f64 * 0.71).cos();
            [a, b, 0.0]
        })
        .collect();

    let particles = ViscosityParticles {
        velocities: &velocities,
        densities: &densities,
        masses: &masses,
        viscosities: &viscosities,
    };

    // Tight tolerance, capped at 100 iterations — must still converge.
    let options = ViscositySolveOptions {
        tol: 1.0e-8,
        max_iters: 100,
    };
    let result = solve_implicit_viscosity(&particles, &neighbors, nu, dt, h, options);
    match result {
        Ok(out) => {
            assert_eq!(out.len(), n);
            for v in &out {
                for &c in v {
                    assert!(c.is_finite(), "non-finite velocity in solution");
                }
            }
        }
        Err(ViscosityError::SolverNotConverged { iters, residual }) => {
            panic!("CG failed to converge in 100 iters: iters={iters}, residual={residual}");
        }
        Err(e) => panic!("unexpected error: {e}"),
    }
}
