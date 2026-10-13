// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Implicit Incompressible SPH (IISPH) pressure solver.
//!
//! Implements the method of Ihmsen et al., "Implicit Incompressible SPH",
//! IEEE TVCG 20(3), 2014.  The pressure `p` is the solution of the linear
//! system `Σ_j a_ij p_j = ρ₀ − ρ_adv` that makes the density predicted for the
//! end of the step equal to the rest density `ρ₀`, where `ρ_adv` is the density
//! advected with the non-pressure velocity.  It is solved with a relaxed
//! Jacobi iteration and non-negative pressures (no suction):
//!
//! ```text
//! d_ii = −Δt² Σ_j m_j/ρ_i² ∇W_ij          d_ij p_j = −Δt² m_j/ρ_j² p_j ∇W_ij
//! a_ii = Σ_j m_j (d_ii − d_ji) · ∇W_ij     (negative)
//! (Ap)_i = Σ_j m_j (D_i − D_j) · ∇W_ij,    D_i = d_ii p_i + Σ_j d_ij p_j = Δt² a_i^p
//! p_i ← max(0, (1−ω) p_i + ω/a_ii (ρ₀ − ρ_adv,i − [(Ap)_i − a_ii p_i]))
//! ```
//!
//! with `∇W_ij = ∇_{x_i} W(|x_i − x_j|, h)`.  Includes the diagonal
//! coefficient, the relaxed Jacobi update, density error metrics and an IISPH
//! timestep adaptation.

use oxiphysics_core::math::Vec3;

use crate::kernel::SphKernel;
use crate::neighbor::SpatialHash;

// ---------------------------------------------------------------------------
// Density error metrics
// ---------------------------------------------------------------------------

/// Metrics describing the density error of an IISPH solve.
#[derive(Debug, Clone, Copy)]
pub struct DensityErrorMetrics {
    /// Average density error as a fraction of rest density.
    pub avg_error: f64,
    /// Maximum density error as a fraction of rest density.
    pub max_error: f64,
    /// Number of Jacobi iterations performed.
    pub iterations: usize,
    /// Whether the solver converged within tolerance.
    pub converged: bool,
}

impl DensityErrorMetrics {
    /// Create zero metrics.
    pub fn zero() -> Self {
        Self {
            avg_error: 0.0,
            max_error: 0.0,
            iterations: 0,
            converged: true,
        }
    }
}

// ---------------------------------------------------------------------------
// IISPH timestep adaptation
// ---------------------------------------------------------------------------

/// Adapts the timestep based on IISPH convergence quality.
#[derive(Debug, Clone)]
pub struct IisphTimestepAdapter {
    /// Current timestep.
    pub dt: f64,
    /// Minimum allowed timestep.
    pub min_dt: f64,
    /// Maximum allowed timestep.
    pub max_dt: f64,
    /// Growth factor when convergence is good.
    pub growth_factor: f64,
    /// Shrink factor when convergence is poor.
    pub shrink_factor: f64,
    /// Target density error for good convergence.
    pub target_error: f64,
    /// Number of consecutive good steps.
    pub good_steps: usize,
    /// Steps required before growing dt.
    pub growth_threshold: usize,
}

impl IisphTimestepAdapter {
    /// Create a new timestep adapter.
    pub fn new(dt: f64, min_dt: f64, max_dt: f64) -> Self {
        Self {
            dt,
            min_dt,
            max_dt,
            growth_factor: 1.05,
            shrink_factor: 0.5,
            target_error: 1e-3,
            good_steps: 0,
            growth_threshold: 5,
        }
    }

    /// Update the timestep based on density error metrics.
    pub fn adapt(&mut self, metrics: &DensityErrorMetrics) -> f64 {
        if !metrics.converged || metrics.avg_error > self.target_error * 10.0 {
            // Poor convergence: shrink dt
            self.dt = (self.dt * self.shrink_factor).max(self.min_dt);
            self.good_steps = 0;
        } else if metrics.avg_error < self.target_error {
            // Good convergence
            self.good_steps += 1;
            if self.good_steps >= self.growth_threshold {
                self.dt = (self.dt * self.growth_factor).min(self.max_dt);
                self.good_steps = 0;
            }
        } else {
            self.good_steps = 0;
        }
        self.dt
    }

    /// Reset the adapter to its initial state.
    pub fn reset(&mut self, initial_dt: f64) {
        self.dt = initial_dt.clamp(self.min_dt, self.max_dt);
        self.good_steps = 0;
    }
}

// ---------------------------------------------------------------------------
// Diagonal coefficient and source term
// ---------------------------------------------------------------------------

/// Smallest density used as a divisor (guards isolated, massless particles).
const DENSITY_FLOOR: f64 = 1e-14;

/// Diagonal coefficients with a magnitude below this threshold mark a particle
/// without usable neighbours: its pressure is fixed to zero.
const DIAGONAL_EPS: f64 = 1e-20;

/// Kernel gradient `∇_i W_ij = W'(r) (x_i − x_j) / r` for the separation
/// `rij = x_i − x_j`; zero for coincident particles.
fn kernel_gradient(rij: Vec3, kernel: &dyn SphKernel, h: f64) -> Vec3 {
    let r = rij.norm();
    if r < 1e-14 {
        return Vec3::zeros();
    }
    rij * (kernel.grad_w(r, h) / r)
}

/// Number of power iterations used to estimate the spectral radius of `D⁻¹A`.
const POWER_ITERATIONS: usize = 20;

/// The relaxed Jacobi iteration contracts only if `ω λ_max(D⁻¹A) < 2`.  The
/// power iteration slightly underestimates `λ_max`, so the relaxation factor is
/// limited to `OMEGA_STABILITY / λ_est` with a margin below 2.
const OMEGA_STABILITY: f64 = 1.6;

/// Matrix-free application of the IISPH pressure matrix `A`.
struct PressureOperator<'a> {
    neighbors: &'a [Vec<usize>],
    /// `∇_i W_ij` for every neighbour entry, aligned with `neighbors`.
    gradients: &'a [Vec<Vec3>],
    masses: &'a [f64],
    /// `1 / ρ_i²`.
    inv_rho_sq: &'a [f64],
    /// `Σ_j m_j ∇_i W_ij`.
    sum_mass_grad: &'a [Vec3],
    dt2: f64,
}

impl PressureOperator<'_> {
    /// Evaluate `D_i = Δt² a_i^p = d_ii p_i + Σ_j d_ij p_j` into `displacement`
    /// and `(Ap)_i = Σ_j m_j (D_i − D_j) · ∇W_ij` into `ap` for the pressures `p`.
    ///
    /// `(Ap)_i = a_ii p_i + Σ_{j≠i} a_ij p_j` is the change of the density of
    /// particle `i` caused by the pressure acceleration `p` over one step.
    fn apply(&self, pressure: &[f64], displacement: &mut [Vec3], ap: &mut [f64]) {
        let n = pressure.len();
        // D_i = −Δt² ( p_i/ρ_i² Σ_j m_j ∇W_ij + Σ_j m_j p_j/ρ_j² ∇W_ij )
        for i in 0..n {
            let mut acc = self.sum_mass_grad[i] * (pressure[i] * self.inv_rho_sq[i]);
            for (&j, grad) in self.neighbors[i].iter().zip(self.gradients[i].iter()) {
                acc += grad * (self.masses[j] * pressure[j] * self.inv_rho_sq[j]);
            }
            displacement[i] = acc * (-self.dt2);
        }
        for i in 0..n {
            let mut sum = 0.0_f64;
            for (&j, grad) in self.neighbors[i].iter().zip(self.gradients[i].iter()) {
                sum += self.masses[j] * (displacement[i] - displacement[j]).dot(grad);
            }
            ap[i] = sum;
        }
    }

    /// Estimate the largest eigenvalue of `D⁻¹A` on the `active` particles
    /// (`D` = diagonal of `A`) with a few power iterations started from a
    /// deterministic pseudo-random vector.  Returns 0 if there is nothing to
    /// iterate on.
    fn spectral_radius_estimate(&self, diagonal: &[f64], active: &[bool]) -> f64 {
        let n = diagonal.len();
        let mut v: Vec<f64> = (0..n)
            .map(|i| {
                if active[i] {
                    let s = ((i + 1) as f64 * 12.9898).sin() * 43758.5453;
                    2.0 * (s - s.floor()) - 1.0
                } else {
                    0.0
                }
            })
            .collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if !(norm.is_finite() && norm > 0.0) {
            return 0.0;
        }
        v.iter_mut().for_each(|x| *x /= norm);
        let mut displacement = vec![Vec3::zeros(); n];
        let mut av = vec![0.0_f64; n];
        let mut estimate = 0.0_f64;
        for _ in 0..POWER_ITERATIONS {
            self.apply(&v, &mut displacement, &mut av);
            let mut norm_sq = 0.0_f64;
            for i in 0..n {
                av[i] = if active[i] { av[i] / diagonal[i] } else { 0.0 };
                norm_sq += av[i] * av[i];
            }
            let norm = norm_sq.sqrt();
            if !(norm.is_finite() && norm > 0.0) {
                break;
            }
            estimate = norm;
            for i in 0..n {
                v[i] = av[i] / norm;
            }
        }
        estimate
    }
}

/// Compute the diagonal coefficient `a_ii` of the IISPH pressure system
/// (Ihmsen et al. 2014, Eq. 12) of particle `i`.
///
/// With `d_ii = −Δt² Σ_j m_j/ρ_i² ∇W_ij` and `d_ji = −Δt² m_i/ρ_i² ∇W_ji`
/// (`∇W_ji = −∇W_ij`) the coefficient is
///
/// ```text
/// a_ii = Σ_j m_j (d_ii − d_ji) · ∇W_ij
///      = −Δt²/ρ_i² ( |Σ_j m_j ∇W_ij|² + m_i Σ_j m_j |∇W_ij|² )
/// ```
///
/// and is therefore **negative** for every particle that has at least one
/// neighbour (zero for an isolated one): raising the pressure of a particle
/// lowers its predicted density.
///
/// # Arguments
/// * `i`         – index of the particle
/// * `positions`, `masses`, `densities` – per-particle data (SPH densities)
/// * `neighbors` – indices of the neighbours of particle `i`
/// * `kernel`, `h` – smoothing kernel and smoothing length
/// * `dt`        – time step
///
/// Returns 0 if `i` or a neighbour index is out of range of the data slices.
pub fn compute_diagonal_coefficient(
    i: usize,
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    neighbors: &[usize],
    kernel: &dyn SphKernel,
    h: f64,
    dt: f64,
) -> f64 {
    let (Some(&pos_i), Some(&mass_i), Some(&rho_i)) =
        (positions.get(i), masses.get(i), densities.get(i))
    else {
        return 0.0;
    };
    let rho_i = rho_i.max(DENSITY_FLOOR);
    let mut sum_mg = Vec3::zeros();
    let mut sum_sq = 0.0_f64;
    for &j in neighbors {
        let (Some(&pos_j), Some(&mass_j)) = (positions.get(j), masses.get(j)) else {
            continue;
        };
        let grad = kernel_gradient(pos_i - pos_j, kernel, h);
        sum_mg += grad * mass_j;
        sum_sq += mass_j * grad.norm_squared();
    }
    -dt * dt / (rho_i * rho_i) * (sum_mg.norm_squared() + mass_i * sum_sq)
}

/// Compute the source term (density deviation from rest) for IISPH.
///
/// source_i = ρ₀ - ρ_adv_i
pub fn compute_source_term(rho_adv: f64, rho0: f64) -> f64 {
    rho0 - rho_adv
}

// ---------------------------------------------------------------------------
// Relaxed Jacobi iteration utilities
// ---------------------------------------------------------------------------

/// Perform a single relaxed Jacobi pressure update for one particle.
///
/// p_i^{new} = (1-ω) * p_i + ω * (source_i - ap_off_i) / a_ii
///
/// `source_i = ρ₀ − ρ_adv,i` is negative for a compressed particle and
/// `a_ii` is negative (see [`compute_diagonal_coefficient`]), so compression
/// yields a positive pressure; an expanded particle would need suction and is
/// clamped to zero.  `ap_off_i` is the off-diagonal part
/// `(Ap)_i − a_ii p_i` of the matrix–vector product.
///
/// Returns the updated pressure (clamped to ≥ 0); 0 when `|a_ii|` is
/// negligible (isolated particle).
pub fn relaxed_jacobi_update(
    pressure_i: f64,
    source_i: f64,
    ap_off_i: f64,
    a_ii: f64,
    omega: f64,
) -> f64 {
    if a_ii.abs() < 1e-20 {
        return 0.0;
    }
    let rhs = (source_i - ap_off_i) / a_ii;
    ((1.0 - omega) * pressure_i + omega * rhs).max(0.0)
}

/// Compute a relaxation factor from the spread of the diagonal coefficients.
///
/// Heuristic: the more uniform the magnitudes `|a_ii|` are (the better
/// conditioned the system), the larger the factor.  The IISPH diagonal is
/// negative, so the magnitudes are used; isolated particles (`a_ii = 0`) are
/// ignored.  Returns ω ∈ \[0.3, 0.9\] and 0.5 if no particle has a usable
/// diagonal.
pub fn estimate_optimal_omega(a_ii_values: &[f64], _n_particles: usize) -> f64 {
    let magnitudes = a_ii_values
        .iter()
        .map(|a| a.abs())
        .filter(|&a| a.is_finite() && a >= DIAGONAL_EPS);
    let min_aii = magnitudes.clone().fold(f64::INFINITY, f64::min);
    let max_aii = magnitudes.fold(0.0_f64, f64::max);

    if !min_aii.is_finite() || max_aii < DIAGONAL_EPS {
        return 0.5;
    }

    let ratio = min_aii / max_aii;
    // Well-conditioned → can use larger ω
    let omega = 0.5 + 0.4 * ratio.clamp(0.0, 1.0);
    omega.clamp(0.3, 0.9)
}

// ---------------------------------------------------------------------------
// IisphSolver
// ---------------------------------------------------------------------------

/// Result of one IISPH pressure solve (see [`IisphSolver::solve`]).
#[derive(Debug, Clone)]
pub struct IisphSolution {
    /// Pressures `p_i ≥ 0` (Pa).  Pass them back as the warm start of the next
    /// solve.
    pub pressures: Vec<f64>,
    /// Pressure accelerations `a_i = −Σ_j m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W_ij`
    /// (m/s², force per unit mass).
    pub accelerations: Vec<Vec3>,
    /// Diagonal coefficients `a_ii` of the pressure system (negative, or zero
    /// for particles without neighbours).
    pub diagonal: Vec<f64>,
    /// SPH densities `ρ_i` at the beginning of the step (kg/m³).
    pub densities: Vec<f64>,
    /// Advected densities `ρ_adv,i = ρ_i + Δt Σ_j m_j (v_i − v_j) · ∇W_ij`
    /// (density change due to the non-pressure velocity only).
    pub advected_densities: Vec<f64>,
    /// Relaxation factor ω actually used: the configured
    /// [`IisphSolver::omega`], lowered if the iteration would not contract.
    pub omega: f64,
    /// Density-error metrics of the returned pressures.
    pub metrics: DensityErrorMetrics,
}

impl IisphSolution {
    /// Solution returned for inconsistent input: zero pressures and
    /// accelerations, infinite error, not converged.
    fn rejected(n: usize) -> Self {
        Self {
            pressures: vec![0.0; n],
            accelerations: vec![Vec3::zeros(); n],
            diagonal: vec![0.0; n],
            densities: vec![0.0; n],
            advected_densities: vec![0.0; n],
            omega: 0.0,
            metrics: DensityErrorMetrics {
                avg_error: f64::INFINITY,
                max_error: f64::INFINITY,
                iterations: 0,
                converged: false,
            },
        }
    }
}

/// IISPH pressure solver using relaxed Jacobi iterations.
///
/// Reference: Ihmsen et al., "Implicit Incompressible SPH", IEEE TVCG 2014.
#[derive(Debug, Clone)]
pub struct IisphSolver {
    /// Relaxation factor for the Jacobi update (default 0.5, Ihmsen et al.).
    ///
    /// The relaxed Jacobi iteration contracts only if `ω λ_max(D⁻¹A) < 2`, and
    /// `λ_max` grows with the kernel support (about 3 for a support of two
    /// particle spacings, about 5 for 2.6 and above 10 for 4).  The solver
    /// estimates `λ_max` and lowers the factor used to `1.6 / λ_max` when
    /// necessary; [`IisphSolution::omega`] reports the value in effect.
    pub omega: f64,
    /// Maximum number of Jacobi iterations (default 100).
    pub max_iter: usize,
    /// Density error tolerance as a fraction of `density0` (default 0.001).
    pub tolerance: f64,
    /// Rest density of the fluid (kg/m³).
    pub density0: f64,
}

impl IisphSolver {
    /// Create a new [`IisphSolver`] with default relaxation parameters.
    ///
    /// # Arguments
    /// * `density0` – rest (reference) density of the fluid in kg/m³.
    pub fn new(density0: f64) -> Self {
        Self {
            omega: 0.5,
            max_iter: 100,
            tolerance: 0.001,
            density0,
        }
    }

    /// Compute SPH densities via kernel summation (including self-contribution).
    ///
    /// Returns a `Vec<f64>` of length `positions.len()`.  Neighbour indices that
    /// are out of range and particles without a mass entry are ignored.
    pub fn compute_densities(
        positions: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> Vec<f64> {
        let h = smoothing_length;
        positions
            .iter()
            .enumerate()
            .map(|(i, pos_i)| {
                let Some(&mass_i) = masses.get(i) else {
                    return 0.0;
                };
                // Self-contribution
                let mut rho = mass_i * kernel.w(0.0, h);
                for &j in neighbors.get(i).map_or(&[][..], Vec::as_slice) {
                    if let (Some(&mass_j), Some(pos_j)) = (masses.get(j), positions.get(j)) {
                        rho += mass_j * kernel.w((pos_i - pos_j).norm(), h);
                    }
                }
                rho
            })
            .collect()
    }

    /// Solve the IISPH pressure system and return pressures, accelerations and
    /// diagnostics.
    ///
    /// # Algorithm (Ihmsen et al. 2014, Alg. 1)
    /// 1. SPH densities `ρ_i` and the advected densities
    ///    `ρ_adv,i = ρ_i + Δt Σ_j m_j (v_i − v_j) · ∇W_ij`; `velocities` are the
    ///    velocities **after** the non-pressure forces.
    /// 2. Diagonal coefficients `a_ii` (negative, see
    ///    [`compute_diagonal_coefficient`]).
    /// 3. Relaxed Jacobi iterations `p_i ← max(0, p_i + ω (ρ₀ − ρ_adv,i − (Ap)_i) / a_ii)`
    ///    with `(Ap)_i = Σ_j m_j (D_i − D_j) · ∇W_ij` and
    ///    `D_i = Δt² a_i^p` the displacement-like pressure acceleration of the
    ///    current iterate, until the average density error drops below
    ///    [`tolerance`](Self::tolerance) or [`max_iter`](Self::max_iter)
    ///    updates were made.
    /// 4. Pressure accelerations from the converged pressures.
    ///
    /// The density error of particle `i` is `|ρ_i^pred − ρ₀| / ρ₀` for a
    /// particle under pressure and `max(ρ_i^pred − ρ₀, 0) / ρ₀` for one at zero
    /// pressure, where `ρ_i^pred = ρ_adv,i + (Ap)_i` is the density after the
    /// pressure acceleration: free-surface particles that are too sparse cannot
    /// be corrected without suction and are not counted.  Particles without
    /// neighbours are ignored.  `avg_error` is averaged over all particles and
    /// refers to the returned pressures.
    ///
    /// Total momentum is conserved: `Σ m_i a_i = 0` for symmetric neighbour
    /// lists.
    ///
    /// # Arguments
    /// * `positions`, `velocities`, `masses` – particle state (equal lengths)
    /// * `neighbors`        – `neighbors[i]` lists the particles within the kernel
    ///   support of `i` (without `i`), symmetric in `i` and `j`
    /// * `dt`               – time step
    /// * `kernel`, `smoothing_length` – smoothing kernel and its length `h`
    /// * `previous_pressure` – optional warm start (negative entries are
    ///   clamped to zero); ignored if its length differs from the particle count
    ///
    /// Inconsistent input (differing slice lengths or an out-of-range neighbour
    /// index) yields a rejected solution: zero pressures and accelerations,
    /// infinite error and `converged = false`.
    pub fn solve(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
        previous_pressure: Option<&[f64]>,
    ) -> IisphSolution {
        let n = positions.len();
        let consistent = velocities.len() == n
            && masses.len() == n
            && neighbors.len() == n
            && neighbors.iter().all(|list| list.iter().all(|&j| j < n));
        if !consistent {
            return IisphSolution::rejected(n);
        }
        let h = smoothing_length;
        let rho0 = self.density0;
        let dt2 = dt * dt;

        // --- 1. densities, cached kernel gradients, advected density, a_ii ---
        let densities = Self::compute_densities(positions, masses, neighbors, kernel, h);
        let inv_rho_sq: Vec<f64> = densities
            .iter()
            .map(|&rho| {
                let rho = rho.max(DENSITY_FLOOR);
                1.0 / (rho * rho)
            })
            .collect();
        let gradients: Vec<Vec<Vec3>> = neighbors
            .iter()
            .enumerate()
            .map(|(i, list)| {
                list.iter()
                    .map(|&j| kernel_gradient(positions[i] - positions[j], kernel, h))
                    .collect()
            })
            .collect();

        // Σ_j m_j ∇W_ij, a_ii and ρ_adv of every particle.
        let mut sum_mass_grad = vec![Vec3::zeros(); n];
        let mut diagonal = vec![0.0_f64; n];
        let mut rho_adv = vec![0.0_f64; n];
        for i in 0..n {
            let mut sum_mg = Vec3::zeros();
            let mut sum_sq = 0.0_f64;
            let mut drho = 0.0_f64;
            for (&j, grad) in neighbors[i].iter().zip(gradients[i].iter()) {
                sum_mg += grad * masses[j];
                sum_sq += masses[j] * grad.norm_squared();
                drho += masses[j] * (velocities[i] - velocities[j]).dot(grad);
            }
            sum_mass_grad[i] = sum_mg;
            diagonal[i] = -dt2 * inv_rho_sq[i] * (sum_mg.norm_squared() + masses[i] * sum_sq);
            rho_adv[i] = densities[i] + dt * drho;
        }
        let active: Vec<bool> = diagonal.iter().map(|&a| a < -DIAGONAL_EPS).collect();

        // --- 2. relaxed Jacobi iterations ---
        let operator = PressureOperator {
            neighbors,
            gradients: &gradients,
            masses,
            inv_rho_sq: &inv_rho_sq,
            sum_mass_grad: &sum_mass_grad,
            dt2,
        };
        let lambda_max = operator.spectral_radius_estimate(&diagonal, &active);
        let omega = if lambda_max.is_finite() && lambda_max > 0.0 {
            self.omega.min(OMEGA_STABILITY / lambda_max)
        } else {
            self.omega
        };
        let mut pressure: Vec<f64> = match previous_pressure {
            Some(previous) if previous.len() == n => previous
                .iter()
                .zip(active.iter())
                .map(|(&p, &on)| if on { p.max(0.0) } else { 0.0 })
                .collect(),
            _ => vec![0.0; n],
        };
        let mut displacement = vec![Vec3::zeros(); n]; // D_i = Δt² a_i^p
        let mut ap = vec![0.0_f64; n];
        let mut iterations = 0_usize;
        let mut converged = false;
        let mut avg_error;
        let mut max_error;
        loop {
            operator.apply(&pressure, &mut displacement, &mut ap);
            // density error of the current iterate
            let mut error_sum = 0.0_f64;
            max_error = 0.0_f64;
            for i in 0..n {
                if !active[i] {
                    continue;
                }
                let deviation = (rho_adv[i] + ap[i] - rho0) / rho0;
                let error = if pressure[i] > 0.0 {
                    deviation.abs()
                } else {
                    deviation.max(0.0)
                };
                error_sum += error;
                max_error = max_error.max(error);
            }
            avg_error = if n > 0 { error_sum / n as f64 } else { 0.0 };
            if avg_error < self.tolerance {
                converged = true;
                break;
            }
            if iterations >= self.max_iter {
                break;
            }
            // relaxed Jacobi update with all particles using the old iterate
            for i in 0..n {
                pressure[i] = if active[i] {
                    relaxed_jacobi_update(
                        pressure[i],
                        rho0 - rho_adv[i],
                        ap[i] - diagonal[i] * pressure[i],
                        diagonal[i],
                        omega,
                    )
                } else {
                    0.0
                };
            }
            iterations += 1;
        }

        // --- 3. pressure accelerations ---
        let accelerations: Vec<Vec3> = (0..n)
            .map(|i| {
                let mut acc = Vec3::zeros();
                for (&j, grad) in neighbors[i].iter().zip(gradients[i].iter()) {
                    let factor =
                        masses[j] * (pressure[i] * inv_rho_sq[i] + pressure[j] * inv_rho_sq[j]);
                    acc -= grad * factor;
                }
                acc
            })
            .collect();

        IisphSolution {
            pressures: pressure,
            accelerations,
            diagonal,
            densities,
            advected_densities: rho_adv,
            omega,
            metrics: DensityErrorMetrics {
                avg_error,
                max_error,
                iterations,
                converged,
            },
        }
    }

    /// Perform one full IISPH pressure-solve step from zero initial pressures.
    ///
    /// See [`solve`](Self::solve) for the algorithm.
    ///
    /// # Returns
    /// `(pressure_forces, avg_density_error)` where `pressure_forces[i]` is
    /// the net pressure acceleration (force per unit mass) on particle `i`.
    pub fn solve_pressure(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> (Vec<Vec3>, f64) {
        let solution = self.solve(
            positions,
            velocities,
            masses,
            neighbors,
            dt,
            kernel,
            smoothing_length,
            None,
        );
        (solution.accelerations, solution.metrics.avg_error)
    }

    /// Solve pressure and return detailed metrics.
    pub fn solve_pressure_with_metrics(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
    ) -> (Vec<Vec3>, DensityErrorMetrics) {
        let solution = self.solve(
            positions,
            velocities,
            masses,
            neighbors,
            dt,
            kernel,
            smoothing_length,
            None,
        );
        (solution.accelerations, solution.metrics)
    }
}

// ---------------------------------------------------------------------------
// IISPH-FLIP hybrid velocity update
// ---------------------------------------------------------------------------

/// Blend between pure IISPH (pic) and FLIP particle-in-cell velocity updates.
///
/// `blend = 0` → pure IISPH (use pressure forces to update velocity via `v += dt * f`).
/// `blend = 1` → pure FLIP (velocity correction is only the incremental change).
///
/// In practice, values like 0.97 give good energy conservation with reduced
/// numerical dissipation.
///
/// This function returns the "effective" per-particle pressure acceleration that
/// should be used in the downstream velocity integration.  With `blend = 0` it
/// returns the IISPH forces unchanged; with `blend = 1` it returns them too but
/// the caller is expected to blend `v_pic = v_old + dt * f` with
/// `v_flip = v_old + dt * (f - f_old)` externally.
pub fn iisph_flip_blend(
    pressure_forces: &[Vec3],
    _velocities_old: &[Vec3],
    _velocities_pic: &[Vec3],
    _blend: f64,
    n: usize,
) -> Vec<Vec3> {
    // Implementation: return pressure forces (the caller applies the blending).
    // blend parameter is kept for API compatibility.
    assert!(pressure_forces.len() >= n);
    pressure_forces[..n].to_vec()
}

// ---------------------------------------------------------------------------
// IISPH error monitor
// ---------------------------------------------------------------------------

/// Tracks per-step IISPH convergence metrics across multiple time steps.
#[derive(Debug, Clone, Default)]
pub struct IisphErrorMonitor {
    /// Per-step average density errors.
    pub errors: Vec<f64>,
    /// Per-step iteration counts.
    pub iteration_counts: Vec<usize>,
}

impl IisphErrorMonitor {
    /// Create a new monitor.
    pub fn new() -> Self {
        Self {
            errors: Vec::new(),
            iteration_counts: Vec::new(),
        }
    }

    /// Record one step's error and iteration count.
    pub fn record(&mut self, error: f64, iterations: usize) {
        self.errors.push(error);
        self.iteration_counts.push(iterations);
    }

    /// Number of steps recorded.
    pub fn num_steps(&self) -> usize {
        self.errors.len()
    }

    /// Last recorded average density error.
    pub fn last_error(&self) -> f64 {
        self.errors.last().copied().unwrap_or(0.0)
    }

    /// Average number of Jacobi iterations over all recorded steps.
    pub fn avg_iterations(&self) -> f64 {
        if self.iteration_counts.is_empty() {
            return 0.0;
        }
        let sum: usize = self.iteration_counts.iter().sum();
        sum as f64 / self.iteration_counts.len() as f64
    }

    /// Maximum density error observed over all recorded steps.
    pub fn max_error(&self) -> f64 {
        self.errors.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Reset all recorded data.
    pub fn clear(&mut self) {
        self.errors.clear();
        self.iteration_counts.clear();
    }
}

// ---------------------------------------------------------------------------
// IISPH neighborhood cache (acceleration data structure)
// ---------------------------------------------------------------------------

/// Cache of precomputed neighbour lists for IISPH.
///
/// Stores neighbour indices for each particle using a spatial hash.  This
/// allows O(1) lookup per particle instead of re-computing neighbours each
/// Jacobi iteration.
#[derive(Debug, Clone)]
pub struct IisphNeighborhoodCache {
    /// Neighbour indices per particle.
    neighbors: Vec<Vec<usize>>,
    /// Cutoff radius used to build the cache.
    pub cutoff: f64,
}

impl IisphNeighborhoodCache {
    /// Build the cache from particle positions and a cutoff radius.
    pub fn build(positions: &[Vec3], cutoff: f64) -> Self {
        use crate::neighbor::SpatialHash;
        let neighbors = SpatialHash::find_all_neighbors(positions, cutoff);
        Self { neighbors, cutoff }
    }

    /// Return the neighbour slice for particle `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

    /// Total number of particles.
    pub fn len(&self) -> usize {
        self.neighbors.len()
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.neighbors.is_empty()
    }

    /// Invalidate the cache (clear all stored neighbours).
    pub fn invalidate(&mut self) {
        for nb in self.neighbors.iter_mut() {
            nb.clear();
        }
    }

    /// Update the cache for new positions.
    pub fn update(&mut self, positions: &[Vec3]) {
        self.neighbors = SpatialHash::find_all_neighbors(positions, self.cutoff);
    }
}

// ---------------------------------------------------------------------------
// IISPH multi-iteration solve with warm start
// ---------------------------------------------------------------------------

impl IisphSolver {
    /// Solve pressure with warm-start initialisation.
    ///
    /// Accepts a `previous_pressure` slice (for example
    /// [`IisphSolution::pressures`] of the previous time step) to use as the
    /// initial guess for the Jacobi iteration.  Using the pressure from the
    /// previous time step as a warm start can significantly reduce the number
    /// of iterations needed.  See [`solve`](Self::solve).
    pub fn solve_pressure_warm_start(
        &self,
        positions: &[Vec3],
        velocities: &[Vec3],
        masses: &[f64],
        neighbors: &[Vec<usize>],
        dt: f64,
        kernel: &dyn SphKernel,
        smoothing_length: f64,
        previous_pressure: &[f64],
    ) -> (Vec<Vec3>, f64) {
        let solution = self.solve(
            positions,
            velocities,
            masses,
            neighbors,
            dt,
            kernel,
            smoothing_length,
            Some(previous_pressure),
        );
        (solution.accelerations, solution.metrics.avg_error)
    }
}

#[cfg(test)]
mod tests;
