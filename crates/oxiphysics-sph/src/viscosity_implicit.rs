// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit / semi-implicit SPH viscosity solver (Takahashi et al. 2015).
//!
//! The explicit viscous update `v* = vⁿ + Δt·ν·∇²v` is conditionally stable,
//! limited to `Δt < h²/(2·d·ν)` (d = spatial dimension).  For viscous-dominated
//! flows this bound is far smaller than the advective CFL limit, so explicit
//! integration wastes work.  The implicit update instead solves the linear
//! system
//!
//! ```text
//! (I − Δt·ν·L) vⁿ⁺¹ = vⁿ
//! ```
//!
//! where `L` is the SPH viscous Laplacian.  Because `L` is negative
//! semi-definite, `A = I − Δt·ν·L` is symmetric **positive definite** and
//! strictly diagonally dominant for every `Δt > 0`, `ν > 0`, so the update is
//! unconditionally stable and conjugate gradient converges.
//!
//! # Discrete operator (Morris / Brookshaw)
//!
//! The Morris et al. (1997) symmetric second-derivative estimator gives, for
//! particle `i` with neighbour `j` (with `r_ij = x_i − x_j`,
//! `∇W_ij = ∇_i W(r_ij, h)`):
//!
//! ```text
//! (L v)_i = Σ_j c_ij (v_j − v_i),
//! c_ij = (m_j / ρ_j) · (μ_i + μ_j) / ρ_i · (−(r_ij · ∇W_ij)) / (|r_ij|² + ε h²).
//! ```
//!
//! The kernel gradient `∇W_ij = (dW/dr) r_ij / |r_ij|` points opposite to
//! `r_ij` (the kernel decreases outward, `dW/dr < 0`), so `r_ij · ∇W_ij < 0`
//! and the leading `−` makes `c_ij ≥ 0`.  Hence the Laplacian matrix has
//! off-diagonal `L_ij = c_ij ≥ 0` and diagonal `L_ii = −Σ_j c_ij ≤ 0`, and the
//! system matrix is
//!
//! ```text
//! A_ii = 1 + Δt·ν·Σ_j c_ij  > 0   (diagonal, dominant)
//! A_ij = −Δt·ν·c_ij          ≤ 0   (off-diagonal)
//! ```
//!
//! a symmetric M-matrix.  The three velocity components are decoupled and share
//! the same matrix `A`, so we assemble `A` once and run three CG solves with
//! right-hand sides `vⁿ_x`, `vⁿ_y`, `vⁿ_z`.
//!
//! The linear algebra (`CsrMatrix`, [`pcg_solve`]) is reused from the ISPH
//! module — self-contained, pure Rust, no external sparse-solver dependency.

use crate::isph::{pcg_solve, CsrMatrix};

/// Small dimensionless regulariser for the Morris denominator `|r_ij|² + ε h²`.
///
/// Matches the value used across the explicit viscosity discretisations and
/// prevents the singularity as `r_ij → 0`.
const MORRIS_EPSILON: f64 = 0.01;

/// One neighbour interaction for the implicit viscosity stencil.
///
/// Precomputing `r_ij` and `∇W_ij` keeps the assembly independent of any
/// particular kernel or neighbour-search implementation: callers fill these in
/// however they already compute SPH interactions.
#[derive(Debug, Clone, Copy)]
pub struct NeighborEntry {
    /// Index `j` of the neighbouring particle.
    pub j: usize,
    /// Displacement `r_ij = x_i − x_j`.
    pub r_ij: [f64; 3],
    /// Kernel gradient `∇W_ij = ∇_i W(r_ij, h)` (points opposite to `r_ij`).
    pub grad_w_ij: [f64; 3],
}

impl NeighborEntry {
    /// Convenience constructor.
    pub fn new(j: usize, r_ij: [f64; 3], grad_w_ij: [f64; 3]) -> Self {
        Self { j, r_ij, grad_w_ij }
    }
}

/// Per-particle state consumed by [`solve_implicit_viscosity`].
///
/// Structure-of-references rather than owned data: the caller keeps its own
/// SoA / AoS storage and lends read-only slices for the duration of the solve.
#[derive(Debug, Clone, Copy)]
pub struct ViscosityParticles<'a> {
    /// Current velocities `[vx, vy, vz]` per particle (the system right-hand
    /// side).
    pub velocities: &'a [[f64; 3]],
    /// Per-particle density `ρ_i` (kg/m³).
    pub densities: &'a [f64],
    /// Per-particle mass `m_i` (kg).
    pub masses: &'a [f64],
    /// Per-particle **dynamic** viscosity `μ_i = ρ_i · ν_i` (Pa·s).
    pub viscosities: &'a [f64],
}

/// Errors that can arise during an implicit viscosity solve.
#[derive(Debug, thiserror::Error)]
pub enum ViscosityError {
    /// CG failed to reach the requested tolerance within the iteration budget.
    #[error("implicit viscosity CG did not converge: {iters} iters, residual {residual:.3e}")]
    SolverNotConverged {
        /// Iterations performed before giving up.
        iters: usize,
        /// Final relative residual `‖b − Ax‖ / ‖b‖`.
        residual: f64,
    },
    /// The particle list was empty.
    #[error("empty particle list")]
    EmptyParticleList,
    /// Slice lengths (`velocities`, `densities`, `masses`, `viscosities`,
    /// `neighbors`) did not all match.
    #[error("inconsistent particle slice lengths")]
    InconsistentLengths,
    /// `dt` or `nu` was not strictly positive, or a density was non-positive.
    #[error("dt and nu must be > 0 and all densities must be > 0")]
    InvalidConfig,
}

/// Relative residual tolerance and iteration cap for the implicit solve.
#[derive(Debug, Clone, Copy)]
pub struct ViscositySolveOptions {
    /// Relative residual tolerance for each component CG solve.
    pub tol: f64,
    /// Maximum CG iterations per component.
    pub max_iters: usize,
}

impl Default for ViscositySolveOptions {
    fn default() -> Self {
        Self {
            tol: 1.0e-8,
            max_iters: 2000,
        }
    }
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the symmetric Morris coupling `c_ij ≥ 0` for one neighbour entry.
///
/// Returns `0` for a coincident pair (`|r_ij|² ≈ 0`) where the stencil is not
/// defined.  `h` enters only through the `ε h²` regulariser.
fn morris_coupling(
    entry: &NeighborEntry,
    mu_i: f64,
    mu_j: f64,
    rho_i: f64,
    rho_j: f64,
    m_j: f64,
    h: f64,
) -> f64 {
    let r2 = dot3(entry.r_ij, entry.r_ij);
    if r2 < 1.0e-24 || rho_i <= 0.0 || rho_j <= 0.0 {
        return 0.0;
    }
    let denom = r2 + MORRIS_EPSILON * h * h;
    // r_ij · ∇W_ij ≤ 0, so the leading `-` yields c_ij ≥ 0.
    let geom = -dot3(entry.r_ij, entry.grad_w_ij) / denom;
    (m_j / rho_j) * (mu_i + mu_j) / rho_i * geom
}

/// Assemble the implicit-viscosity system matrix `A = I − Δt·ν·L`.
///
/// Off-diagonal `A_ij = −Δt·ν·c_ij`, diagonal `A_ii = 1 + Δt·ν·Σ_j c_ij`.
/// Duplicate neighbour entries for the same `j` are accumulated into a single
/// column so the resulting CSR has one stored value per `(i, j)` pair.  `nu` is
/// the scalar multiplying the whole Laplacian; the per-particle `μ_i`, `μ_j`
/// supply the spatially varying part of the coupling.
fn assemble_system(
    particles: &ViscosityParticles<'_>,
    neighbors: &[Vec<NeighborEntry>],
    nu: f64,
    dt: f64,
    h: f64,
) -> CsrMatrix {
    let n = particles.velocities.len();
    let mut off_diagonals: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut diagonal = vec![1.0_f64; n];
    let scale = dt * nu;
    for i in 0..n {
        let mu_i = particles.viscosities[i];
        let rho_i = particles.densities[i];
        // Accumulate couplings per neighbour column (handles duplicates).
        let mut columns: Vec<(usize, f64)> = Vec::with_capacity(neighbors[i].len());
        for entry in &neighbors[i] {
            let j = entry.j;
            if j == i || j >= n {
                continue;
            }
            let c = morris_coupling(
                entry,
                mu_i,
                particles.viscosities[j],
                rho_i,
                particles.densities[j],
                particles.masses[j],
                h,
            );
            if c == 0.0 {
                continue;
            }
            match columns.iter_mut().find(|(col, _)| *col == j) {
                Some((_, acc)) => *acc += c,
                None => columns.push((j, c)),
            }
        }
        for (j, c) in columns {
            diagonal[i] += scale * c; // A_ii = 1 + dt·ν·Σ c_ij
            off_diagonals[i].push((j, -scale * c)); // A_ij = −dt·ν·c_ij
        }
    }
    CsrMatrix::from_rows(n, &off_diagonals, &diagonal)
}

/// Implicit viscosity update via Jacobi-preconditioned conjugate gradient.
///
/// Solves `(I − Δt·ν·L) vⁿ⁺¹ = vⁿ` for the new velocity field, where `L` is the
/// Morris symmetric SPH viscous Laplacian (see the module documentation).  The
/// three velocity components share the assembled SPD matrix and are solved
/// independently.
///
/// # Arguments
/// - `particles`: current velocities, densities, masses and dynamic
///   viscosities (`μ_i = ρ_i·ν_i`).
/// - `neighbors`: `neighbors[i]` lists every interaction `(j, r_ij, ∇W_ij)` of
///   particle `i`.
/// - `nu`: kinematic viscosity scale `ν` (m²/s) multiplying the Laplacian.
/// - `dt`: time step `Δt` (s).
/// - `h`: smoothing length (only used in the `ε h²` regulariser).
/// - `options`: CG tolerance and iteration cap.
///
/// # Errors
/// - [`ViscosityError::EmptyParticleList`] if there are no particles.
/// - [`ViscosityError::InconsistentLengths`] if the slices disagree in length.
/// - [`ViscosityError::InvalidConfig`] if `dt ≤ 0`, `nu ≤ 0`, or any density is
///   non-positive.
/// - [`ViscosityError::SolverNotConverged`] if any component CG exceeds the
///   iteration budget without reaching `options.tol`.
///
/// # Returns
/// `Ok(velocities)` with the updated `[vx, vy, vz]` for every particle, indexed
/// to match `particles.velocities`.
pub fn solve_implicit_viscosity(
    particles: &ViscosityParticles<'_>,
    neighbors: &[Vec<NeighborEntry>],
    nu: f64,
    dt: f64,
    h: f64,
    options: ViscositySolveOptions,
) -> Result<Vec<[f64; 3]>, ViscosityError> {
    let n = particles.velocities.len();
    if n == 0 {
        return Err(ViscosityError::EmptyParticleList);
    }
    if particles.densities.len() != n
        || particles.masses.len() != n
        || particles.viscosities.len() != n
        || neighbors.len() != n
    {
        return Err(ViscosityError::InconsistentLengths);
    }
    if dt <= 0.0 || nu <= 0.0 || h <= 0.0 {
        return Err(ViscosityError::InvalidConfig);
    }
    if particles.densities.iter().any(|&rho| rho <= 0.0) {
        return Err(ViscosityError::InvalidConfig);
    }

    let a = assemble_system(particles, neighbors, nu, dt, h);

    let mut result = vec![[0.0_f64; 3]; n];
    let mut rhs = vec![0.0_f64; n];
    let mut sol = vec![0.0_f64; n];
    // The x/y/z velocity components are decoupled and share the SPD matrix `a`,
    // so solve them one at a time, scattering each solution back into `result`.
    for comp in 0..3 {
        for (i, (r, s)) in rhs.iter_mut().zip(sol.iter_mut()).enumerate() {
            *r = particles.velocities[i][comp];
            *s = particles.velocities[i][comp]; // warm start
        }
        let (iters, residual) = pcg_solve(&a, &rhs, &mut sol, options.tol, options.max_iters);
        // A zero right-hand side legitimately returns residual 0; only a residual
        // above tolerance (or a non-finite one) is a genuine non-convergence.
        if !residual.is_finite() || residual > options.tol {
            return Err(ViscosityError::SolverNotConverged { iters, residual });
        }
        for (out, &si) in result.iter_mut().zip(sol.iter()) {
            out[comp] = si;
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{grad, CubicSplineKernel};

    /// Build a neighbour list for a set of positions using the cubic-spline
    /// kernel and support radius `2h`.  Shared by several tests.
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

    #[test]
    fn coupling_is_non_negative() {
        // r_ij · ∇W_ij < 0 ⇒ c_ij ≥ 0 for any positive masses/densities.
        let kernel = CubicSplineKernel;
        let h = 0.1;
        let r_ij = [0.08_f64, 0.0, 0.0];
        let gw = grad(&kernel, r_ij, h);
        let entry = NeighborEntry::new(1, r_ij, gw);
        let c = morris_coupling(&entry, 1.0, 1.0, 1000.0, 1000.0, 0.001, h);
        assert!(c >= 0.0, "Morris coupling must be non-negative, got {c}");
        assert!(c > 0.0, "coupling should be strictly positive for r < 2h");
    }

    #[test]
    fn coupling_zero_for_coincident() {
        let entry = NeighborEntry::new(1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let c = morris_coupling(&entry, 1.0, 1.0, 1000.0, 1000.0, 0.001, 0.1);
        assert_eq!(c, 0.0, "coincident particles contribute nothing");
    }

    #[test]
    fn empty_particle_list_errors() {
        let particles = ViscosityParticles {
            velocities: &[],
            densities: &[],
            masses: &[],
            viscosities: &[],
        };
        let neighbors: Vec<Vec<NeighborEntry>> = Vec::new();
        let err = solve_implicit_viscosity(
            &particles,
            &neighbors,
            0.1,
            0.01,
            0.1,
            ViscositySolveOptions::default(),
        );
        assert!(matches!(err, Err(ViscosityError::EmptyParticleList)));
    }

    #[test]
    fn inconsistent_lengths_error() {
        let velocities = vec![[0.0; 3]; 2];
        let densities = vec![1000.0; 2];
        let masses = vec![0.001; 2];
        let viscosities = vec![1.0]; // wrong length
        let particles = ViscosityParticles {
            velocities: &velocities,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        let neighbors = vec![Vec::new(); 2];
        let err = solve_implicit_viscosity(
            &particles,
            &neighbors,
            0.1,
            0.01,
            0.1,
            ViscositySolveOptions::default(),
        );
        assert!(matches!(err, Err(ViscosityError::InconsistentLengths)));
    }

    #[test]
    fn invalid_config_errors() {
        let velocities = vec![[0.0; 3]; 2];
        let densities = vec![1000.0; 2];
        let masses = vec![0.001; 2];
        let viscosities = vec![1.0; 2];
        let particles = ViscosityParticles {
            velocities: &velocities,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        let neighbors = vec![Vec::new(); 2];
        // dt <= 0
        let err = solve_implicit_viscosity(
            &particles,
            &neighbors,
            0.1,
            0.0,
            0.1,
            ViscositySolveOptions::default(),
        );
        assert!(matches!(err, Err(ViscosityError::InvalidConfig)));
        // nu <= 0
        let err = solve_implicit_viscosity(
            &particles,
            &neighbors,
            -1.0,
            0.01,
            0.1,
            ViscositySolveOptions::default(),
        );
        assert!(matches!(err, Err(ViscosityError::InvalidConfig)));
    }

    #[test]
    fn isolated_particle_velocity_unchanged() {
        // No neighbours ⇒ A = I ⇒ vⁿ⁺¹ = vⁿ.
        let velocities = vec![[1.0, -2.0, 3.0]];
        let densities = vec![1000.0];
        let masses = vec![0.001];
        let viscosities = vec![1.0];
        let particles = ViscosityParticles {
            velocities: &velocities,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        let neighbors = vec![Vec::new()];
        let out = solve_implicit_viscosity(
            &particles,
            &neighbors,
            0.1,
            0.01,
            0.1,
            ViscositySolveOptions::default(),
        )
        .expect("solve should succeed");
        for k in 0..3 {
            assert!(
                (out[0][k] - velocities[0][k]).abs() < 1e-12,
                "isolated particle velocity must be unchanged"
            );
        }
    }

    #[test]
    fn two_particles_velocities_relax_together() {
        // Two equal particles with opposite velocities: viscous diffusion must
        // pull them toward their mean (here zero), conserving total momentum.
        let h = 0.1;
        let positions = vec![[0.0, 0.0, 0.0], [0.08, 0.0, 0.0]];
        let velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let densities = vec![1000.0, 1000.0];
        let masses = vec![1.0, 1.0];
        let viscosities = vec![1000.0, 1000.0]; // μ = ρ·ν, ν = 1
        let neighbors = build_neighbors(&positions, h);
        let particles = ViscosityParticles {
            velocities: &velocities,
            densities: &densities,
            masses: &masses,
            viscosities: &viscosities,
        };
        let out = solve_implicit_viscosity(
            &particles,
            &neighbors,
            1.0,
            0.5,
            h,
            ViscositySolveOptions::default(),
        )
        .expect("solve should succeed");
        // Magnitudes shrink toward the mean.
        assert!(out[0][0] < velocities[0][0], "particle 0 vx must decrease");
        assert!(out[1][0] > velocities[1][0], "particle 1 vx must increase");
        // Equal-mass momentum is conserved by the symmetric operator.
        let p_before = masses[0] * velocities[0][0] + masses[1] * velocities[1][0];
        let p_after = masses[0] * out[0][0] + masses[1] * out[1][0];
        assert!(
            (p_before - p_after).abs() < 1e-9,
            "momentum drift {}",
            (p_before - p_after).abs()
        );
    }
}
