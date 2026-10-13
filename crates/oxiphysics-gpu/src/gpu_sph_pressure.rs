// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU SPH pressure solver (CPU mock backend).
//!
//! Implements the Smoothed Particle Hydrodynamics pressure and viscosity
//! force computation pipeline using Poly6, Spiky, and Viscosity kernels.
//! [`pcisph_gpu_correction`] adds the predictive-corrective incompressible
//! (PCISPH) pressure solve on top of it.
//! The GPU dispatch is mocked by plain Rust loops for portability.

use std::f64::consts::PI;

// ── Data structures ──────────────────────────────────────────────────────────

/// GPU SPH pressure solver holding all per-particle state.
#[derive(Debug, Clone)]
pub struct GpuSphPressureSolver {
    /// Total number of simulated particles.
    pub n_particles: usize,
    /// Per-particle pressure values (Pa).
    pub pressure: Vec<f64>,
    /// Per-particle density values (kg/m³).
    pub density: Vec<f64>,
    /// Per-particle positions `[x, y, z]` in metres.
    pub positions: Vec<[f64; 3]>,
    /// Per-particle masses (kg).
    pub masses: Vec<f64>,
    /// Smoothing length h (m).
    pub smoothing_h: f64,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Stiffness constant k for the linearised EOS.
    pub stiffness: f64,
}

impl GpuSphPressureSolver {
    /// Create a new solver with `n` particles, smoothing length `h`,
    /// rest density `rho0`, and stiffness `k`.
    pub fn new(n: usize, h: f64, rho0: f64, k: f64) -> Self {
        Self {
            n_particles: n,
            pressure: vec![0.0; n],
            density: vec![0.0; n],
            positions: vec![[0.0; 3]; n],
            masses: vec![1.0; n],
            smoothing_h: h,
            rest_density: rho0,
            stiffness: k,
        }
    }

    /// Return the number of particles in this solver.
    pub fn particle_count(&self) -> usize {
        self.n_particles
    }

    /// Set the position of particle `i`.
    pub fn set_position(&mut self, i: usize, pos: [f64; 3]) {
        self.positions[i] = pos;
    }

    /// Set the mass of particle `i`.
    pub fn set_mass(&mut self, i: usize, mass: f64) {
        self.masses[i] = mass;
    }

    /// Compute total mass of all particles.
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }

    /// Compute the maximum relative density error: max |ρᵢ − ρ₀| / ρ₀.
    pub fn density_error(&self) -> f64 {
        if self.rest_density <= 0.0 {
            return 0.0;
        }
        self.density
            .iter()
            .map(|&rho| (rho - self.rest_density).abs() / self.rest_density)
            .fold(0.0_f64, f64::max)
    }

    /// Compute density statistics.
    pub fn compute_stats(&self) -> GpuSphStats {
        let max_density = self
            .density
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_density = self.density.iter().cloned().fold(f64::INFINITY, f64::min);
        let mean_pressure = if self.n_particles == 0 {
            0.0
        } else {
            self.pressure.iter().sum::<f64>() / self.n_particles as f64
        };
        let compression_error = self.density_error();
        GpuSphStats {
            max_density,
            min_density,
            mean_pressure,
            compression_error,
        }
    }

    /// Mock GPU density computation: ρᵢ = Σⱼ mⱼ · W_poly6(|rᵢ − rⱼ|, h).
    pub fn gpu_compute_density(&mut self) {
        let n = self.n_particles;
        let h = self.smoothing_h;
        for i in 0..n {
            let mut rho = 0.0f64;
            for j in 0..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                rho += self.masses[j] * kernel_poly6(r, h);
            }
            self.density[i] = rho;
        }
    }

    /// Mock GPU pressure computation: pᵢ = k · (ρᵢ − ρ₀).
    pub fn gpu_compute_pressure(&mut self) {
        let k = self.stiffness;
        let rho0 = self.rest_density;
        for i in 0..self.n_particles {
            self.pressure[i] = k * (self.density[i] - rho0);
        }
    }

    /// Compute pressure force on particle `i`:
    /// F_press_i = −Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇W_spiky.
    pub fn gpu_pressure_force(&self, i: usize) -> [f64; 3] {
        let h = self.smoothing_h;
        let rhoi = self.density[i];
        let pi = self.pressure[i];
        let mut force = [0.0f64; 3];
        if rhoi < 1e-15 {
            return force;
        }
        for j in 0..self.n_particles {
            if i == j {
                continue;
            }
            let rhoj = self.density[j];
            if rhoj < 1e-15 {
                continue;
            }
            let r_vec = [
                self.positions[i][0] - self.positions[j][0],
                self.positions[i][1] - self.positions[j][1],
                self.positions[i][2] - self.positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let grad = kernel_spiky_grad(r_vec, r, h);
            let coeff = -self.masses[j] * (pi / (rhoi * rhoi) + self.pressure[j] / (rhoj * rhoj));
            force[0] += coeff * grad[0];
            force[1] += coeff * grad[1];
            force[2] += coeff * grad[2];
        }
        force
    }

    /// Compute viscosity force on particle `i`:
    /// F_visc_i = μ Σⱼ mⱼ/ρⱼ (vⱼ − vᵢ) ∇²W_visc.
    pub fn gpu_viscosity_force(&self, i: usize, velocities: &[[f64; 3]], mu: f64) -> [f64; 3] {
        let h = self.smoothing_h;
        let mut force = [0.0f64; 3];
        for j in 0..self.n_particles {
            if i == j {
                continue;
            }
            let rhoj = self.density[j];
            if rhoj < 1e-15 {
                continue;
            }
            let r_vec = [
                self.positions[i][0] - self.positions[j][0],
                self.positions[i][1] - self.positions[j][1],
                self.positions[i][2] - self.positions[j][2],
            ];
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let lap = kernel_viscosity_laplacian(r, h);
            let dv = [
                velocities[j][0] - velocities[i][0],
                velocities[j][1] - velocities[i][1],
                velocities[j][2] - velocities[i][2],
            ];
            let coeff = mu * self.masses[j] / rhoj * lap;
            force[0] += coeff * dv[0];
            force[1] += coeff * dv[1];
            force[2] += coeff * dv[2];
        }
        force
    }
}

/// Statistics computed from the pressure solver state.
#[derive(Debug, Clone)]
pub struct GpuSphStats {
    /// Maximum density across all particles.
    pub max_density: f64,
    /// Minimum density across all particles.
    pub min_density: f64,
    /// Mean pressure across all particles.
    pub mean_pressure: f64,
    /// Maximum relative density deviation from rest density.
    pub compression_error: f64,
}

// ── Kernel functions ─────────────────────────────────────────────────────────

/// Poly6 SPH kernel: W_poly6(r, h).
///
/// Returns `315/(64·π·h⁹) · (h²−r²)³` for `r < h`, else `0`.
pub fn kernel_poly6(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r >= h {
        return 0.0;
    }
    let h2 = h * h;
    let diff = h2 - r * r;
    315.0 / (64.0 * PI * h.powi(9)) * diff.powi(3)
}

/// Gradient of the Poly6 kernel with respect to the first particle:
/// ∇W_poly6(r_vec, r, h) with `r_vec = xᵢ − xⱼ` and `r = |r_vec|`.
///
/// Returns `−945/(32·π·h⁹) · (h²−r²)² · r_vec` for `r < h`, else the zero vector.
pub fn kernel_poly6_grad(r_vec: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if h <= 0.0 || r >= h {
        return [0.0; 3];
    }
    let diff = h * h - r * r;
    let coeff = -945.0 / (32.0 * PI * h.powi(9)) * diff * diff;
    [coeff * r_vec[0], coeff * r_vec[1], coeff * r_vec[2]]
}

/// Gradient of the Spiky kernel: ∇W_spiky(r_vec, r, h).
///
/// Returns `−45/(π·h⁶) · (h−r)² · r̂` for `r < h` and `r > 0`, else zero vector.
pub fn kernel_spiky_grad(r_vec: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if h <= 0.0 || r >= h || r < 1e-15 {
        return [0.0; 3];
    }
    let coeff = -45.0 / (PI * h.powi(6)) * (h - r).powi(2) / r;
    [coeff * r_vec[0], coeff * r_vec[1], coeff * r_vec[2]]
}

/// Viscosity kernel Laplacian: ∇²W_visc(r, h).
///
/// Returns `45/(π·h⁶) · (h−r)` for `r < h`, else `0`.
pub fn kernel_viscosity_laplacian(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r >= h {
        return 0.0;
    }
    45.0 / (PI * h.powi(6)) * (h - r)
}

// ── Standalone utilities ─────────────────────────────────────────────────────

/// Summary of a PCISPH pressure-correction run, see [`pcisph_gpu_correction_report`].
#[derive(Debug, Clone)]
pub struct PcisphReport {
    /// Number of pressure-update iterations performed, including updates that were
    /// discarded again (see [`rejected_steps`](Self::rejected_steps)).
    pub iterations: usize,
    /// `true` when the largest relative residual of the returned pressure field is
    /// below the requested tolerance.
    pub converged: bool,
    /// Number of tentative updates that were discarded because they let the
    /// residual 2-norm grow beyond 1.5 times the best value seen; each one restarts
    /// the iteration from the best state with half the update gain.
    pub rejected_steps: usize,
    /// Largest relative residual `max_i r_i / ρ₀` of the zero-pressure start
    /// (entry 0) followed by the state after every accepted update; see
    /// [`pcisph_gpu_correction_report`] for the definition of `r_i`.  The history
    /// is not necessarily monotone: the residual of a healthy PCISPH run may rise
    /// temporarily before it falls.  It is empty when the input was rejected.
    pub error_history: Vec<f64>,
    /// Largest relative residual of the pressure field left in the solver, or
    /// `f64::NAN` when the input was rejected and nothing was computed.
    pub final_error: f64,
    /// Pressure-update gain `δ` obtained from the prototype particle (`0` when the
    /// input was rejected before it was computed).
    pub delta: f64,
}

/// Largest per-axis lattice extent, in particle spacings, that the PCISPH
/// prototype scan is willing to enumerate (bounds its cost to ~10⁶ lattice sites).
const PCISPH_PROTOTYPE_MAX_REACH: f64 = 48.0;

/// A tentative pressure update is discarded when its residual 2-norm exceeds this
/// multiple of the best 2-norm seen so far, which signals an overshooting update.
const PCISPH_DIVERGENCE_FACTOR: f64 = 1.5;

/// Smallest damping factor `ω` the PCISPH update is reduced to after discarded steps.
const PCISPH_MIN_RELAXATION: f64 = 1.0 / 128.0;

/// Consecutive discarded updates at the smallest damping factor after which the
/// iteration gives up (the second one is a deterministic repeat of the first).
const PCISPH_MAX_STALLED_RETRIES: usize = 2;

/// One state of the PCISPH iteration: a pressure field, the densities predicted for
/// it and the residual norms of those densities.
#[derive(Clone)]
struct PcisphState {
    pressure: Vec<f64>,
    rho_star: Vec<f64>,
    max_error: f64,
    l2_error: f64,
}

/// PCISPH pressure-update gain `δ` from a prototype particle with a full neighbourhood.
///
/// The prototype sits at the centre of a cubic lattice whose spacing
/// `s = (m/ρ₀)^(1/3)` reproduces the rest density.  Linearising the predicted
/// density `ρ*ᵢ = Σⱼ mⱼ W(x*ᵢ − x*ⱼ)` for a pressure `p` acting on the particle and
/// its neighbours (`Δx = Δt²·a_p`, `a_p = −Σⱼ mⱼ (pᵢ + pⱼ)/ρ₀² ∇W_spikyᵢⱼ`) gives
///
/// `Δρᵢ = −β p (Σⱼ ∇Wρᵢⱼ · Σⱼ ∇Wpᵢⱼ + Σⱼ ∇Wρᵢⱼ · ∇Wpᵢⱼ)`, `β = 2 (Δt m / ρ₀)²`,
///
/// where `∇Wρ` is the Poly6 (density) gradient and `∇Wp` the Spiky (pressure)
/// gradient.  Requiring `Δρᵢ = −(ρ*ᵢ − ρ₀)` yields `δ = 1 / (β (…))`.
///
/// Returns `None` when the inputs are not positive/finite, when the prototype has no
/// neighbour inside the support (`h` smaller than the lattice spacing, so the
/// response vanishes) or when `h` spans more than [`PCISPH_PROTOTYPE_MAX_REACH`]
/// lattice spacings.
fn pcisph_prototype_delta(h: f64, rest_density: f64, mass: f64, dt: f64) -> Option<f64> {
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if !(positive(h) && positive(rest_density) && positive(mass) && positive(dt)) {
        return None;
    }
    let spacing = (mass / rest_density).cbrt();
    if !positive(spacing) {
        return None;
    }
    let reach_f = (h / spacing).ceil();
    if !reach_f.is_finite() || reach_f > PCISPH_PROTOTYPE_MAX_REACH {
        return None;
    }
    let reach = reach_f as i64;
    let mut sum_density_grad = [0.0f64; 3];
    let mut sum_pressure_grad = [0.0f64; 3];
    let mut pair_dot = 0.0f64;
    for ix in -reach..=reach {
        for iy in -reach..=reach {
            for iz in -reach..=reach {
                if ix == 0 && iy == 0 && iz == 0 {
                    continue;
                }
                let r_vec = [
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                ];
                let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
                if r >= h {
                    continue;
                }
                let grad_density = kernel_poly6_grad(r_vec, r, h);
                let grad_pressure = kernel_spiky_grad(r_vec, r, h);
                for dim in 0..3 {
                    sum_density_grad[dim] += grad_density[dim];
                    sum_pressure_grad[dim] += grad_pressure[dim];
                    pair_dot += grad_density[dim] * grad_pressure[dim];
                }
            }
        }
    }
    let sum_dot = sum_density_grad[0] * sum_pressure_grad[0]
        + sum_density_grad[1] * sum_pressure_grad[1]
        + sum_density_grad[2] * sum_pressure_grad[2];
    let beta = 2.0 * (dt * mass / rest_density).powi(2);
    let denominator = beta * (sum_dot + pair_dot);
    if denominator.is_finite() && denominator > 0.0 {
        Some(1.0 / denominator)
    } else {
        None
    }
}

/// `true` if every component of every 3-vector in `values` is finite.
fn all_components_finite(values: &[[f64; 3]]) -> bool {
    values.iter().all(|v| v.iter().all(|c| c.is_finite()))
}

/// Poly6 densities at the PCISPH-predicted positions
/// `x*ᵢ = xᵢ + Δt (vᵢ + Δt (a_extᵢ + a_pressureᵢ))`.
fn pcisph_predicted_density(
    solver: &GpuSphPressureSolver,
    velocities: &[[f64; 3]],
    ext_acc: &[[f64; 3]],
    pressure_acc: &[[f64; 3]],
    dt: f64,
) -> Vec<f64> {
    let h = solver.smoothing_h;
    let predicted: Vec<[f64; 3]> = (0..solver.n_particles)
        .map(|i| {
            let mut x = solver.positions[i];
            for (dim, coord) in x.iter_mut().enumerate() {
                let v_star = velocities[i][dim] + dt * (ext_acc[i][dim] + pressure_acc[i][dim]);
                *coord += dt * v_star;
            }
            x
        })
        .collect();
    predicted
        .iter()
        .map(|xi| {
            predicted
                .iter()
                .zip(solver.masses.iter())
                .map(|(xj, &mj)| {
                    let dx = xi[0] - xj[0];
                    let dy = xi[1] - xj[1];
                    let dz = xi[2] - xj[2];
                    mj * kernel_poly6((dx * dx + dy * dy + dz * dz).sqrt(), h)
                })
                .sum()
        })
        .collect()
}

/// Complementarity residual of the predicted densities, `(max-norm, 2-norm)`, both
/// relative to `ρ₀`.
///
/// Particles carrying pressure (`pᵢ > 0`) must be at rest density, so they
/// contribute `|ρ*ᵢ − ρ₀|`; pressure-free particles may only be at or below it,
/// so they contribute the compression `max(ρ*ᵢ − ρ₀, 0)`.  Density deficiency at a
/// free surface can never be repaired by a non-negative pressure and is therefore
/// not counted against particles without pressure.
///
/// A NaN density makes both norms NaN instead of being silently ignored, so that
/// every comparison against it fails.
fn pcisph_residual_norms(rho_star: &[f64], pressure: &[f64], rest_density: f64) -> (f64, f64) {
    let mut max_norm = 0.0f64;
    let mut sum_sq = 0.0f64;
    for (&rho, &p) in rho_star.iter().zip(pressure.iter()) {
        let err = rho - rest_density;
        if err.is_nan() {
            return (f64::NAN, f64::NAN);
        }
        let residual = (if p > 0.0 { err.abs() } else { err.max(0.0) }) / rest_density;
        max_norm = max_norm.max(residual);
        sum_sq += residual * residual;
    }
    (max_norm, sum_sq.sqrt())
}

/// Predictive-corrective incompressible SPH (PCISPH, Solenthaler & Pajarola 2009)
/// pressure solve; returns the full [`PcisphReport`].
///
/// Starting from zero pressure, every iteration
///
/// 1. predicts velocities and positions of all particles from the external
///    acceleration plus the pressure acceleration of the current pressure field,
///    `v* = v + Δt (a_ext + a_p)`, `x* = x + Δt v*`;
/// 2. evaluates the Poly6 density `ρ*` at the predicted positions;
/// 3. updates the pressure with `pᵢ ← max(0, pᵢ + ω δ (ρ*ᵢ − ρ₀))`, where `δ` is
///    the gain of a prototype particle with a full neighbourhood (a cubic lattice of
///    spacing `(m̄/ρ₀)^(1/3)`, `m̄` the mean mass, linearised as in the PCISPH paper
///    but with the Poly6 density gradient and the Spiky pressure gradient) and
///    `ω ≤ 1` is a damping factor that starts at `1`;
/// 4. recomputes the pressure accelerations from the new pressure field with
///    [`GpuSphPressureSolver::gpu_pressure_force`] at the *current* positions.
///
/// The prototype gain is exact only for a single particle with a regular, full
/// neighbourhood, so on dense or irregular neighbourhoods the plain fixed-point
/// iteration can overshoot, oscillate or even diverge.  The iteration therefore
/// remembers the best state found so far (smallest residual 2-norm).  An update
/// whose residual 2-norm exceeds 1.5 times that best value is discarded: the
/// iteration restarts from the best state with `ω` halved (down to `1/128`).
/// Updates that merely fail to improve, but stay within that bound, are kept
/// because the residual of a healthy PCISPH run is not monotone.  The run stops as
/// soon as the largest residual drops below `tol`, after `max_iter` iterations, or
/// when updates keep being discarded at the smallest damping factor.  When it stops
/// without meeting `tol` the solver keeps the best state seen, so the returned
/// pressure field is never worse than the zero-pressure start.
///
/// The residual of particle `i` is `|ρ*ᵢ − ρ₀|/ρ₀` if it carries pressure and
/// `max(ρ*ᵢ − ρ₀, 0)/ρ₀` otherwise, i.e. only compression is penalised where no
/// pressure acts, because density deficiency near free surfaces cannot be
/// repaired by a non-negative pressure.
///
/// On return `solver.pressure` holds the corrected pressure (consistent with
/// `solver.gpu_pressure_force`) and `solver.density` the Poly6 density at the
/// unchanged current positions, which `solver.positions` still holds.
///
/// Invalid input — `velocities`/`ext_acc` not holding one entry per particle, any
/// non-finite position, velocity, acceleration or mass, `dt`, `h`, `ρ₀` or the mean
/// mass not finite and positive, or an `h` that spans no neighbour or more than 48
/// particle spacings — leaves the solver untouched and returns an unconverged
/// report with zero iterations.
pub fn pcisph_gpu_correction_report(
    solver: &mut GpuSphPressureSolver,
    velocities: &[[f64; 3]],
    ext_acc: &[[f64; 3]],
    dt: f64,
    max_iter: usize,
    tol: f64,
) -> PcisphReport {
    let n = solver.n_particles;
    let rho0 = solver.rest_density;
    let mut report = PcisphReport {
        iterations: 0,
        converged: false,
        rejected_steps: 0,
        error_history: Vec::new(),
        final_error: f64::NAN,
        delta: 0.0,
    };
    let consistent = velocities.len() == n
        && ext_acc.len() == n
        && solver.positions.len() == n
        && solver.masses.len() == n
        && solver.pressure.len() == n
        && solver.density.len() == n
        && all_components_finite(velocities)
        && all_components_finite(ext_acc)
        && all_components_finite(&solver.positions)
        && solver.masses.iter().all(|m| m.is_finite());
    if !(consistent && dt.is_finite() && dt > 0.0) {
        return report;
    }
    if n == 0 {
        report.error_history.push(0.0);
        report.final_error = 0.0;
        report.converged = 0.0 < tol;
        return report;
    }
    let mean_mass = solver.total_mass() / n as f64;
    let Some(delta) = pcisph_prototype_delta(solver.smoothing_h, rho0, mean_mass, dt) else {
        return report;
    };
    report.delta = delta;

    solver.gpu_compute_density();
    solver.pressure.fill(0.0);
    let zero_acc = vec![[0.0f64; 3]; n];
    let start_rho = pcisph_predicted_density(solver, velocities, ext_acc, &zero_acc, dt);
    let (start_max, start_l2) = pcisph_residual_norms(&start_rho, &solver.pressure, rho0);
    let mut current = PcisphState {
        pressure: vec![0.0f64; n],
        rho_star: start_rho,
        max_error: start_max,
        l2_error: start_l2,
    };
    let mut best = current.clone();
    report.error_history.push(start_max);

    let mut omega = 1.0f64;
    let mut stalled = 0usize;
    while current.max_error >= tol && report.iterations < max_iter {
        report.iterations += 1;
        for ((p_new, &p_old), &rho) in solver
            .pressure
            .iter_mut()
            .zip(current.pressure.iter())
            .zip(current.rho_star.iter())
        {
            *p_new = (p_old + omega * delta * (rho - rho0)).max(0.0);
        }
        let pressure_acc: Vec<[f64; 3]> = (0..n).map(|i| solver.gpu_pressure_force(i)).collect();
        let rho_star = pcisph_predicted_density(solver, velocities, ext_acc, &pressure_acc, dt);
        let (max_error, l2_error) = pcisph_residual_norms(&rho_star, &solver.pressure, rho0);
        if max_error < tol || l2_error <= PCISPH_DIVERGENCE_FACTOR * best.l2_error {
            current = PcisphState {
                pressure: solver.pressure.clone(),
                rho_star,
                max_error,
                l2_error,
            };
            report.error_history.push(max_error);
            stalled = 0;
            if current.l2_error < best.l2_error {
                best = current.clone();
            }
        } else {
            report.rejected_steps += 1;
            current = best.clone();
            if omega > PCISPH_MIN_RELAXATION {
                omega = (0.5 * omega).max(PCISPH_MIN_RELAXATION);
            } else {
                stalled += 1;
                if stalled >= PCISPH_MAX_STALLED_RETRIES {
                    break;
                }
            }
        }
    }

    // `best` never has a larger 2-norm than any accepted state, so it is the right
    // fallback whenever the tolerance was not met.
    let result = if current.max_error < tol {
        current
    } else {
        best
    };
    solver.pressure.copy_from_slice(&result.pressure);
    report.final_error = result.max_error;
    report.converged = result.max_error < tol;
    report
}

/// PCISPH pressure correction: finds the pressure field that removes the
/// predicted compression of the particle set.
///
/// `velocities` and `ext_acc` hold the current velocity and the external
/// (non-pressure, e.g. gravity and viscosity) acceleration of every particle, and
/// `dt` is the time step the correction is made for.  Iterates at most `max_iter`
/// times until the largest relative density residual falls below `tol`; the
/// solver's `pressure` field then holds the corrected pressures.
///
/// Returns `(iterations, converged)`.  See [`pcisph_gpu_correction_report`] for the
/// algorithm, the residual definition and the full diagnostics.
pub fn pcisph_gpu_correction(
    solver: &mut GpuSphPressureSolver,
    velocities: &[[f64; 3]],
    ext_acc: &[[f64; 3]],
    dt: f64,
    max_iter: usize,
    tol: f64,
) -> (usize, bool) {
    let report = pcisph_gpu_correction_report(solver, velocities, ext_acc, dt, max_iter, tol);
    (report.iterations, report.converged)
}

/// WCSPH Tait equation of state.
///
/// Returns `ρ₀·cs²/γ · ((ρ/ρ₀)^γ − 1)`.
pub fn wcsph_tait_eos(rho: f64, rho0: f64, cs: f64, gamma: f64) -> f64 {
    if rho0 <= 0.0 || gamma <= 0.0 {
        return 0.0;
    }
    rho0 * cs * cs / gamma * ((rho / rho0).powf(gamma) - 1.0)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── kernel tests ─────────────────────────────────────────────────────

    #[test]
    fn test_poly6_positive_within_h() {
        let w = kernel_poly6(0.5, 1.0);
        assert!(w > 0.0, "poly6 should be positive for r < h, got {w}");
    }

    #[test]
    fn test_poly6_zero_at_h() {
        let w = kernel_poly6(1.0, 1.0);
        assert_eq!(w, 0.0, "poly6 should be zero at r == h");
    }

    #[test]
    fn test_poly6_zero_beyond_h() {
        let w = kernel_poly6(1.5, 1.0);
        assert_eq!(w, 0.0, "poly6 should be zero for r > h");
    }

    #[test]
    fn test_poly6_at_origin_positive() {
        let w = kernel_poly6(0.0, 1.0);
        assert!(w > 0.0, "poly6 at r=0 should be positive");
    }

    #[test]
    fn test_poly6_decreasing() {
        let w0 = kernel_poly6(0.0, 1.0);
        let w1 = kernel_poly6(0.4, 1.0);
        let w2 = kernel_poly6(0.8, 1.0);
        assert!(w0 > w1 && w1 > w2);
    }

    #[test]
    fn test_poly6_zero_h() {
        assert_eq!(kernel_poly6(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_spiky_grad_zero_outside_h() {
        let g = kernel_spiky_grad([1.0, 0.0, 0.0], 1.0, 0.5);
        assert_eq!(g, [0.0; 3]);
    }

    #[test]
    fn test_spiky_grad_nonzero_within_h() {
        let r_vec = [0.3, 0.0, 0.0];
        let r = 0.3_f64;
        let g = kernel_spiky_grad(r_vec, r, 1.0);
        let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!(mag > 0.0, "spiky grad should be nonzero within h");
    }

    #[test]
    fn test_spiky_grad_points_radially() {
        // For r_vec along x axis, gradient should be along x only
        let r_vec = [0.4, 0.0, 0.0];
        let g = kernel_spiky_grad(r_vec, 0.4, 1.0);
        assert!(g[1].abs() < 1e-15 && g[2].abs() < 1e-15);
    }

    #[test]
    fn test_viscosity_laplacian_positive_within_h() {
        let lap = kernel_viscosity_laplacian(0.5, 1.0);
        assert!(lap > 0.0);
    }

    #[test]
    fn test_viscosity_laplacian_zero_at_h() {
        let lap = kernel_viscosity_laplacian(1.0, 1.0);
        assert_eq!(lap, 0.0);
    }

    #[test]
    fn test_viscosity_laplacian_zero_beyond_h() {
        let lap = kernel_viscosity_laplacian(1.5, 1.0);
        assert_eq!(lap, 0.0);
    }

    // ── wcsph_tait_eos tests ──────────────────────────────────────────────

    #[test]
    fn test_wcsph_zero_at_rest_density() {
        let p = wcsph_tait_eos(1000.0, 1000.0, 100.0, 7.0);
        assert!(
            p.abs() < 1e-6,
            "WCSPH pressure at rest density should be zero, got {p}"
        );
    }

    #[test]
    fn test_wcsph_positive_above_rest() {
        let p = wcsph_tait_eos(1100.0, 1000.0, 100.0, 7.0);
        assert!(p > 0.0, "WCSPH pressure should be positive for rho > rho0");
    }

    #[test]
    fn test_wcsph_negative_below_rest() {
        let p = wcsph_tait_eos(900.0, 1000.0, 100.0, 7.0);
        assert!(p < 0.0, "WCSPH pressure should be negative for rho < rho0");
    }

    #[test]
    fn test_wcsph_zero_rho0() {
        assert_eq!(wcsph_tait_eos(1000.0, 0.0, 100.0, 7.0), 0.0);
    }

    // ── solver constructor tests ──────────────────────────────────────────

    #[test]
    fn test_solver_new_particle_count() {
        let s = GpuSphPressureSolver::new(5, 1.0, 1000.0, 1.0);
        assert_eq!(s.particle_count(), 5);
    }

    #[test]
    fn test_solver_new_initial_pressure_zero() {
        let s = GpuSphPressureSolver::new(3, 1.0, 1000.0, 1.0);
        assert!(s.pressure.iter().all(|&p| p == 0.0));
    }

    #[test]
    fn test_solver_total_mass() {
        let mut s = GpuSphPressureSolver::new(3, 1.0, 1000.0, 1.0);
        s.masses = vec![1.0, 2.0, 3.0];
        assert!((s.total_mass() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_solver_set_position() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_position(0, [1.0, 2.0, 3.0]);
        assert_eq!(s.positions[0], [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_solver_set_mass() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_mass(1, 5.0);
        assert!((s.masses[1] - 5.0).abs() < 1e-12);
    }

    // ── density / pressure compute tests ─────────────────────────────────

    #[test]
    fn test_gpu_compute_density_positive() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        s.gpu_compute_density();
        assert!(s.density[0] > 0.0 && s.density[1] > 0.0);
    }

    #[test]
    fn test_gpu_compute_density_self_contribution() {
        // A single particle should self-contribute via W_poly6(0, h) > 0
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 1.0);
        s.gpu_compute_density();
        assert!(s.density[0] > 0.0);
    }

    #[test]
    fn test_gpu_compute_pressure_nonneg_above_rho0() {
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 200.0);
        s.density[0] = 1100.0; // above rest
        s.gpu_compute_pressure();
        assert!(s.pressure[0] > 0.0);
    }

    #[test]
    fn test_gpu_compute_pressure_zero_at_rest() {
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 200.0);
        s.density[0] = 1000.0; // exactly at rest
        s.gpu_compute_pressure();
        assert!(s.pressure[0].abs() < 1e-10);
    }

    // ── stats tests ───────────────────────────────────────────────────────

    #[test]
    fn test_stats_max_density_ge_min() {
        let mut s = GpuSphPressureSolver::new(3, 1.0, 1000.0, 1.0);
        s.density = vec![900.0, 1000.0, 1100.0];
        let stats = s.compute_stats();
        assert!(stats.max_density >= stats.min_density);
    }

    #[test]
    fn test_stats_mean_pressure() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.pressure = vec![10.0, 20.0];
        let stats = s.compute_stats();
        assert!((stats.mean_pressure - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_density_error_zero_at_rest() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.density = vec![1000.0, 1000.0];
        assert!(s.density_error().abs() < 1e-12);
    }

    // ── pcisph correction tests ───────────────────────────────────────────

    /// Lattice spacing at which unit-mass particles reproduce `rest_density`.
    fn rest_spacing(rest_density: f64) -> f64 {
        (1.0 / rest_density).cbrt()
    }

    /// Cubic block of `nx³` unit-mass particles with the given lattice `spacing`;
    /// particle `(ix, iy, iz)` has index `(ix·nx + iy)·nx + iz`.
    fn lattice_block(nx: usize, spacing: f64, h: f64, rest_density: f64) -> GpuSphPressureSolver {
        let mut s = GpuSphPressureSolver::new(nx * nx * nx, h, rest_density, 1.0);
        for ix in 0..nx {
            for iy in 0..nx {
                for iz in 0..nx {
                    s.set_position(
                        (ix * nx + iy) * nx + iz,
                        [
                            ix as f64 * spacing,
                            iy as f64 * spacing,
                            iz as f64 * spacing,
                        ],
                    );
                }
            }
        }
        s
    }

    fn zero_vectors(n: usize) -> Vec<[f64; 3]> {
        vec![[0.0; 3]; n]
    }

    /// Residual norms `(max, 2-norm)` of the solver's current pressure field,
    /// recomputed from scratch and therefore independent of the bookkeeping inside
    /// the PCISPH loop.
    fn residual_of_current_pressure(
        s: &GpuSphPressureSolver,
        velocities: &[[f64; 3]],
        ext_acc: &[[f64; 3]],
        dt: f64,
    ) -> (f64, f64) {
        let pressure_acc: Vec<[f64; 3]> = (0..s.n_particles)
            .map(|i| s.gpu_pressure_force(i))
            .collect();
        let rho_star = pcisph_predicted_density(s, velocities, ext_acc, &pressure_acc, dt);
        pcisph_residual_norms(&rho_star, &s.pressure, s.rest_density)
    }

    /// Two 3×3×3 rest-lattice blocks that approach each other along x: block A
    /// occupies the lattice columns 0..=2, block B the columns 3..=5, and the
    /// velocities close the gap between them by `closing` spacings within one step.
    /// Particle `(column, iy, iz)` has index `column·9 + iy·3 + iz`.
    fn colliding_blocks(
        h_factor: f64,
        closing: f64,
        dt: f64,
    ) -> (GpuSphPressureSolver, Vec<[f64; 3]>) {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let mut s = GpuSphPressureSolver::new(54, h_factor * spacing, rho0, 1.0);
        let speed = 0.5 * closing * spacing / dt;
        let mut velocities = Vec::with_capacity(54);
        let mut idx = 0;
        for column in 0..6 {
            for iy in 0..3 {
                for iz in 0..3 {
                    s.set_position(
                        idx,
                        [
                            column as f64 * spacing,
                            iy as f64 * spacing,
                            iz as f64 * spacing,
                        ],
                    );
                    velocities.push([if column < 3 { speed } else { -speed }, 0.0, 0.0]);
                    idx += 1;
                }
            }
        }
        (s, velocities)
    }

    #[test]
    fn test_pcisph_returns_iterations_le_max() {
        // Two sparse particles: far below rest density, so no pressure is needed.
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        let zero = zero_vectors(2);
        let (iters, converged) = pcisph_gpu_correction(&mut s, &zero, &zero, 1e-3, 10, 1e-3);
        assert_eq!((iters, converged), (0, true));
        assert!(s.pressure.iter().all(|&p| p == 0.0));
        // rho = m·W(0) + m·W(0.3) = 315/(64π)·(1 + 0.91³) ≈ 2.75.
        let expected = 315.0 / (64.0 * PI) * (1.0 + 0.91_f64.powi(3));
        for &rho in &s.density {
            assert!(
                (rho - expected).abs() < 1e-12,
                "density {rho}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_pcisph_single_particle_density_is_self_contribution() {
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 1.0);
        let zero = zero_vectors(1);
        let (iters, converged) = pcisph_gpu_correction(&mut s, &zero, &zero, 1e-3, 1, 1.0);
        assert_eq!((iters, converged), (0, true));
        // W_poly6(0, h) = 315 / (64 π h³) for the unit mass.
        let expected = 315.0 / (64.0 * PI);
        assert!(
            (s.density[0] - expected).abs() < 1e-12,
            "density {}, expected {expected}",
            s.density[0]
        );
    }

    #[test]
    fn test_pcisph_compressed_block_error_strictly_decreases() {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let dt = 1e-3;
        let tol = 1e-2;
        for &(nx, h_factor, squeeze) in &[(4usize, 2.5, 0.9), (4, 2.0, 0.9), (3, 2.5, 0.9)] {
            let mut s = lattice_block(nx, squeeze * spacing, h_factor * spacing, rho0);
            let n = s.n_particles;
            let zero = zero_vectors(n);
            let report = pcisph_gpu_correction_report(&mut s, &zero, &zero, dt, 100, tol);
            let case = format!("nx={nx} h={h_factor}·s squeeze={squeeze}");
            assert!(report.converged, "{case}: {report:?}");
            assert!(report.delta > 0.0 && report.delta.is_finite(), "{case}");
            assert!(
                report.error_history[0] > 0.15,
                "{case}: a block squeezed by 10 % must start compressed, got {}",
                report.error_history[0]
            );
            assert!(
                report.error_history.len() >= 4,
                "{case}: expected several corrections, got {:?}",
                report.error_history
            );
            assert!(
                report.error_history.windows(2).all(|w| w[1] < w[0]),
                "{case}: the error must strictly decrease, got {:?}",
                report.error_history
            );
            assert!(report.final_error < tol, "{case}: {}", report.final_error);
            assert_eq!(
                report.error_history.last().copied(),
                Some(report.final_error),
                "{case}"
            );
            // The report is not taken on trust: recompute the residual of the
            // pressure field the solver was left with.
            let (max_error, _) = residual_of_current_pressure(&s, &zero, &zero, dt);
            assert!(max_error < tol, "{case}: recomputed residual {max_error}");
            assert!((max_error - report.final_error).abs() < 1e-12, "{case}");
            assert!(
                s.pressure.iter().all(|&p| p.is_finite() && p >= 0.0),
                "{case}: pressures must be finite and non-negative"
            );
            assert!(s.pressure.iter().any(|&p| p > 0.0), "{case}");
        }
    }

    #[test]
    fn test_pcisph_pressure_pushes_a_compressed_block_outwards() {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let nx = 4;
        let mut s = lattice_block(nx, 0.9 * spacing, 2.5 * spacing, rho0);
        let n = s.n_particles;
        let zero = zero_vectors(n);
        let report = pcisph_gpu_correction_report(&mut s, &zero, &zero, 1e-3, 100, 1e-2);
        assert!(report.converged);
        let centre = 0.5 * (nx as f64 - 1.0) * 0.9 * spacing;
        let mut checked = 0;
        for ix in [0, nx - 1] {
            for iy in [0, nx - 1] {
                for iz in [0, nx - 1] {
                    let i = (ix * nx + iy) * nx + iz;
                    let force = s.gpu_pressure_force(i);
                    let outward = [
                        s.positions[i][0] - centre,
                        s.positions[i][1] - centre,
                        s.positions[i][2] - centre,
                    ];
                    let radial: f64 = force.iter().zip(outward.iter()).map(|(f, o)| f * o).sum();
                    assert!(
                        radial > 0.0,
                        "corner {i}: pressure acceleration {force:?} must point away from the centre"
                    );
                    checked += 1;
                }
            }
        }
        assert_eq!(checked, 8);
        // The pressure field of a symmetric block is point-symmetric.
        let peak = s.pressure.iter().copied().fold(0.0_f64, f64::max);
        assert!(peak > 0.0);
        for i in 0..n {
            assert!(
                (s.pressure[i] - s.pressure[n - 1 - i]).abs() <= 1e-9 * peak,
                "pressure of particle {i} differs from its mirror image"
            );
        }
    }

    #[test]
    fn test_pcisph_colliding_blocks_build_pressure_only_at_the_contact() {
        let dt = 1e-3;
        let tol = 1e-2;
        for &(h_factor, closing) in &[(2.5, 0.5), (2.0, 0.5), (2.5, 0.3)] {
            let (mut s, velocities) = colliding_blocks(h_factor, closing, dt);
            let zero = zero_vectors(54);
            let report = pcisph_gpu_correction_report(&mut s, &velocities, &zero, dt, 100, tol);
            let case = format!("h={h_factor}·s closing={closing}·s");
            assert!(report.converged, "{case}: {report:?}");
            assert!(
                report.error_history[0] > 0.05,
                "{case}: the approach must predict compression, got {}",
                report.error_history[0]
            );
            assert!(
                report.error_history.windows(2).all(|w| w[1] < w[0]),
                "{case}: the error must strictly decrease, got {:?}",
                report.error_history
            );
            let (max_error, _) = residual_of_current_pressure(&s, &velocities, &zero, dt);
            assert!(max_error < tol, "{case}: recomputed residual {max_error}");
            // Only the two fully surrounded particles on the contact plane are
            // compressed: every particle with a free-surface deficiency stays at
            // zero pressure, and the two blocks are mirror images of each other.
            let (centre_a, centre_b) = (2 * 9 + 4, 3 * 9 + 4);
            let contact_peak = s.pressure[centre_a];
            assert!(contact_peak > 0.0, "{case}");
            assert!(
                (contact_peak - s.pressure[centre_b]).abs() <= 1e-9 * contact_peak,
                "{case}: the two contact particles must carry equal pressure ({} vs {})",
                contact_peak,
                s.pressure[centre_b]
            );
            for (i, &p) in s.pressure.iter().enumerate() {
                if i != centre_a && i != centre_b {
                    assert!(
                        p <= 1e-9 * contact_peak,
                        "{case}: particle {i} is not in contact but carries pressure {p}"
                    );
                }
            }
            // The contact pressure decelerates the approach: block A is pushed
            // back (-x) and block B forward (+x).
            let a_centre = s.gpu_pressure_force(centre_a);
            let b_centre = s.gpu_pressure_force(centre_b);
            assert!(a_centre[0] < 0.0, "{case}: {a_centre:?}");
            assert!(b_centre[0] > 0.0, "{case}: {b_centre:?}");
            assert!((a_centre[0] + b_centre[0]).abs() <= 1e-9 * a_centre[0].abs());
        }
    }

    #[test]
    fn test_pcisph_external_acceleration_drives_the_prediction() {
        // A lattice at rest density with no velocity is already fine; an external
        // acceleration towards the centre that moves the outermost particles
        // inwards by 0.3 spacings within one step is not.
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let dt = 1e-3;
        let nx = 3;
        let mut s = lattice_block(nx, spacing, 2.5 * spacing, rho0);
        let n = s.n_particles;
        let zero = zero_vectors(n);
        let calm = pcisph_gpu_correction_report(&mut s.clone(), &zero, &zero, dt, 50, 0.1);
        assert!(
            calm.converged && calm.iterations == 0,
            "the rest lattice needs no pressure: {calm:?}"
        );
        let centre = spacing * (nx as f64 - 1.0) * 0.5;
        let accel = 0.3 * spacing / (dt * dt);
        let ext: Vec<[f64; 3]> = s
            .positions
            .iter()
            .map(|p| {
                [
                    -accel * (p[0] - centre) / centre,
                    -accel * (p[1] - centre) / centre,
                    -accel * (p[2] - centre) / centre,
                ]
            })
            .collect();
        let report = pcisph_gpu_correction_report(&mut s, &zero, &ext, dt, 100, 0.05);
        assert!(report.error_history[0] > 0.3, "{report:?}");
        assert!(report.converged, "{report:?}");
        assert!(report.iterations >= 2, "{report:?}");
        let (max_error, _) = residual_of_current_pressure(&s, &zero, &ext, dt);
        assert!(max_error < 0.05, "recomputed residual {max_error}");
    }

    #[test]
    fn test_pcisph_limited_budget_keeps_best_state() {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let dt = 1e-3;
        let mut s = lattice_block(5, 0.7 * spacing, 2.0 * spacing, rho0);
        let n = s.n_particles;
        let zero = zero_vectors(n);
        // Residual of the zero-pressure start, evaluated independently.
        s.gpu_compute_density();
        let start_rho = pcisph_predicted_density(&s, &zero, &zero, &zero, dt);
        let (start_max, start_l2) = pcisph_residual_norms(&start_rho, &s.pressure, rho0);
        assert!(start_max > 1.5, "a 30 % squeeze is heavily compressed");

        let report = pcisph_gpu_correction_report(&mut s, &zero, &zero, dt, 4, 1e-3);
        assert!(
            !report.converged,
            "four iterations cannot remove 190 % error"
        );
        assert_eq!(report.iterations, 4);
        assert!((report.error_history[0] - start_max).abs() < 1e-12);
        let (end_max, end_l2) = residual_of_current_pressure(&s, &zero, &zero, dt);
        assert!(end_l2 < 0.25 * start_l2, "2-norm {start_l2} -> {end_l2}");
        assert!(
            end_max < 0.5 * start_max,
            "max norm {start_max} -> {end_max}"
        );
        assert!((end_max - report.final_error).abs() < 1e-12);
        assert!(s.pressure.iter().all(|&p| p.is_finite() && p >= 0.0));

        // With a budget of zero iterations the solver keeps zero pressure.
        let mut idle = lattice_block(5, 0.7 * spacing, 2.0 * spacing, rho0);
        let report = pcisph_gpu_correction_report(&mut idle, &zero, &zero, dt, 0, 1e-3);
        assert_eq!(report.iterations, 0);
        assert!(!report.converged);
        assert!(idle.pressure.iter().all(|&p| p == 0.0));
        assert!((report.final_error - start_max).abs() < 1e-12);
    }

    #[test]
    fn test_pcisph_damping_tames_an_overshooting_update() {
        // With ~110 neighbours per particle the prototype gain overshoots by more
        // than a factor of two: the undamped update flips between zero pressure
        // and a hugely over-expanded state.  The best-state rollback with a
        // halved gain must still converge.
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let dt = 1e-3;
        let mut s = lattice_block(6, 0.9 * spacing, 3.0 * spacing, rho0);
        let n = s.n_particles;
        let zero = zero_vectors(n);
        let report = pcisph_gpu_correction_report(&mut s, &zero, &zero, dt, 60, 1e-2);
        assert!(report.converged, "{report:?}");
        assert!(report.rejected_steps >= 1, "{report:?}");
        assert!(report.iterations >= report.rejected_steps + 2, "{report:?}");
        let (max_error, _) = residual_of_current_pressure(&s, &zero, &zero, dt);
        assert!(max_error < 1e-2, "recomputed residual {max_error}");
    }

    #[test]
    fn test_pcisph_iteration_is_independent_of_the_time_step() {
        // delta scales with 1/dt² and the position change with dt², so for zero
        // velocities the iteration must not depend on dt at all.
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let build = || lattice_block(4, 0.9 * spacing, 2.5 * spacing, rho0);
        let zero = zero_vectors(64);
        let mut a = build();
        let mut b = build();
        let fine = pcisph_gpu_correction_report(&mut a, &zero, &zero, 1e-3, 100, 1e-2);
        let coarse = pcisph_gpu_correction_report(&mut b, &zero, &zero, 4e-3, 100, 1e-2);
        assert!(
            (fine.delta / coarse.delta - 16.0).abs() < 1e-9 * 16.0,
            "{} vs {}",
            fine.delta,
            coarse.delta
        );
        assert_eq!(fine.iterations, coarse.iterations);
        assert_eq!(fine.error_history.len(), coarse.error_history.len());
        for (x, y) in fine.error_history.iter().zip(coarse.error_history.iter()) {
            assert!((x - y).abs() < 1e-9, "{x} vs {y}");
        }
    }

    #[test]
    fn test_pcisph_wrapper_matches_report() {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let template = lattice_block(3, 0.85 * spacing, 2.5 * spacing, rho0);
        let zero = zero_vectors(27);
        let mut for_report = template.clone();
        let mut for_wrapper = template;
        let report = pcisph_gpu_correction_report(&mut for_report, &zero, &zero, 1e-3, 100, 1e-2);
        let (iterations, converged) =
            pcisph_gpu_correction(&mut for_wrapper, &zero, &zero, 1e-3, 100, 1e-2);
        assert_eq!(iterations, report.iterations);
        assert_eq!(converged, report.converged);
        assert!(iterations > 0 && converged);
        assert_eq!(for_report.pressure, for_wrapper.pressure);
    }

    #[test]
    fn test_pcisph_rejects_invalid_input_without_touching_the_solver() {
        let rho0 = 1000.0;
        let spacing = rest_spacing(rho0);
        let mut s = lattice_block(3, 0.8 * spacing, 2.5 * spacing, rho0);
        s.pressure.fill(7.0);
        s.density.fill(3.0);
        let n = s.n_particles;
        let zero = zero_vectors(n);
        let mut nan_velocity = zero_vectors(n);
        nan_velocity[4][1] = f64::NAN;
        let short = zero_vectors(n - 1);
        let cases = [
            ("short velocities", short.clone(), zero.clone(), 1e-3),
            ("short accelerations", zero.clone(), short, 1e-3),
            ("NaN velocity", nan_velocity, zero.clone(), 1e-3),
            ("zero dt", zero.clone(), zero.clone(), 0.0),
            ("negative dt", zero.clone(), zero.clone(), -1e-3),
            ("NaN dt", zero.clone(), zero.clone(), f64::NAN),
            ("infinite dt", zero.clone(), zero.clone(), f64::INFINITY),
        ];
        for (name, velocities, ext_acc, dt) in cases {
            let report = pcisph_gpu_correction_report(&mut s, &velocities, &ext_acc, dt, 10, 1e-3);
            assert_eq!(report.iterations, 0, "{name}");
            assert!(!report.converged, "{name}");
            assert!(report.final_error.is_nan(), "{name}");
            assert!(report.error_history.is_empty(), "{name}");
            assert!(s.pressure.iter().all(|&p| p == 7.0), "{name}");
            assert!(s.density.iter().all(|&rho| rho == 3.0), "{name}");
        }
        // Non-finite particle data and unusable solver parameters.
        let mut bad_position = s.clone();
        bad_position.positions[2][0] = f64::INFINITY;
        let report = pcisph_gpu_correction_report(&mut bad_position, &zero, &zero, 1e-3, 10, 1e-3);
        assert_eq!(report.iterations, 0);
        assert!(!report.converged);
        let mut bad_mass = s.clone();
        bad_mass.masses[0] = f64::NAN;
        let report = pcisph_gpu_correction_report(&mut bad_mass, &zero, &zero, 1e-3, 10, 1e-3);
        assert!(!report.converged && report.error_history.is_empty());
        let mut tiny_support = s.clone();
        tiny_support.smoothing_h = 0.5 * spacing; // no lattice neighbour inside the support
        let report = pcisph_gpu_correction_report(&mut tiny_support, &zero, &zero, 1e-3, 10, 1e-3);
        assert_eq!(report.iterations, 0);
        assert!(!report.converged && report.delta == 0.0);
        assert!(tiny_support.pressure.iter().all(|&p| p == 7.0));
    }

    #[test]
    fn test_pcisph_empty_solver() {
        let mut s = GpuSphPressureSolver::new(0, 1.0, 1000.0, 1.0);
        let report = pcisph_gpu_correction_report(&mut s, &[], &[], 1e-3, 5, 1e-3);
        assert_eq!(report.iterations, 0);
        assert!(report.converged);
        assert_eq!(report.final_error, 0.0);
        let report = pcisph_gpu_correction_report(&mut s, &[], &[], 1e-3, 5, 0.0);
        assert!(!report.converged, "a zero tolerance can never be met");
    }

    #[test]
    fn test_pcisph_prototype_delta_matches_closed_form_lattice_sum() {
        // For the symmetric prototype lattice Σⱼ ∇Wⱼ vanishes, so
        // δ = 1 / (β Σⱼ ∇Wρⱼ·∇Wpⱼ) with the closed form
        // ∇Wρ·∇Wp = 945/(32π h⁹)(h²−r²)² · 45/(π h⁶)(h−r)² · r.
        let (h, rho0, mass, dt): (f64, f64, f64, f64) = (0.25, 1000.0, 1.0, 2e-3);
        let spacing = (mass / rho0).cbrt();
        let reach = (h / spacing).ceil() as i64;
        let mut pair_sum = 0.0;
        for ix in -reach..=reach {
            for iy in -reach..=reach {
                for iz in -reach..=reach {
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue;
                    }
                    let r = spacing * ((ix * ix + iy * iy + iz * iz) as f64).sqrt();
                    if r >= h {
                        continue;
                    }
                    pair_sum += 945.0 / (32.0 * PI * h.powi(9)) * (h * h - r * r).powi(2) * 45.0
                        / (PI * h.powi(6))
                        * (h - r).powi(2)
                        * r;
                }
            }
        }
        let beta = 2.0 * (dt * mass / rho0).powi(2);
        let expected = 1.0 / (beta * pair_sum);
        let got = pcisph_prototype_delta(h, rho0, mass, dt).expect("valid prototype");
        assert!(
            (got - expected).abs() <= 1e-9 * expected,
            "delta {got}, closed form {expected}"
        );
        // Unusable inputs.
        assert!(pcisph_prototype_delta(h, rho0, mass, 0.0).is_none());
        assert!(pcisph_prototype_delta(h, -1.0, mass, dt).is_none());
        assert!(pcisph_prototype_delta(h, rho0, f64::NAN, dt).is_none());
        assert!(
            pcisph_prototype_delta(0.5 * spacing, rho0, mass, dt).is_none(),
            "support smaller than the lattice spacing has no neighbour"
        );
        assert!(
            pcisph_prototype_delta(100.0 * spacing, rho0, mass, dt).is_none(),
            "support spanning more than 48 spacings is refused"
        );
    }

    #[test]
    fn test_pcisph_predicted_density_uses_the_predicted_positions() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        s.set_mass(0, 1.0);
        s.set_mass(1, 2.0);
        let velocities = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let ext_acc = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let pressure_acc = [[0.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
        let dt = 0.1;
        // x*₁ = 0.3 + 0.1·(−1 + 0.1·(0 − 2)) = 0.18, x*₀ = 0.
        let rho = pcisph_predicted_density(&s, &velocities, &ext_acc, &pressure_acc, dt);
        let w_self = kernel_poly6(0.0, 1.0);
        let w_pair = kernel_poly6(0.18, 1.0);
        assert!((rho[0] - (w_self + 2.0 * w_pair)).abs() < 1e-12);
        assert!((rho[1] - (w_pair + 2.0 * w_self)).abs() < 1e-12);
    }

    #[test]
    fn test_pcisph_residual_norms_complementarity() {
        // (density, pressure): compressed with pressure, over-expanded with
        // pressure, exactly at rest, compressed without pressure, under-dense
        // without pressure (a free surface, not counted).
        let rho = [1500.0, 900.0, 1000.0, 1200.0, 800.0];
        let pressure = [1.0, 1.0, 0.0, 0.0, 0.0];
        let (max_norm, l2_norm) = pcisph_residual_norms(&rho, &pressure, 1000.0);
        assert!((max_norm - 0.5).abs() < 1e-15);
        let expected_l2 = (0.5_f64 * 0.5 + 0.1 * 0.1 + 0.0 + 0.2 * 0.2 + 0.0).sqrt();
        assert!((l2_norm - expected_l2).abs() < 1e-15);
        let (nan_max, nan_l2) = pcisph_residual_norms(&[1000.0, f64::NAN], &[0.0, 0.0], 1000.0);
        assert!(nan_max.is_nan() && nan_l2.is_nan());
    }

    #[test]
    fn test_poly6_grad_matches_finite_difference() {
        let h = 0.8;
        let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let r_vec = [0.3, -0.2, 0.25];
        let analytic = kernel_poly6_grad(r_vec, norm(r_vec), h);
        let eps = 1e-6;
        for dim in 0..3 {
            let mut plus = r_vec;
            let mut minus = r_vec;
            plus[dim] += eps;
            minus[dim] -= eps;
            let numeric =
                (kernel_poly6(norm(plus), h) - kernel_poly6(norm(minus), h)) / (2.0 * eps);
            assert!(
                (analytic[dim] - numeric).abs() <= 1e-6 * numeric.abs(),
                "component {dim}: analytic {} vs finite difference {numeric}",
                analytic[dim]
            );
        }
        // The kernel decreases with distance, so the gradient points towards the
        // neighbour (opposite to r_vec), and it vanishes outside the support.
        let dot: f64 = analytic.iter().zip(r_vec.iter()).map(|(g, r)| g * r).sum();
        assert!(dot < 0.0);
        assert_eq!(kernel_poly6_grad([h, 0.0, 0.0], h, h), [0.0; 3]);
        assert_eq!(kernel_poly6_grad([0.1, 0.0, 0.0], 0.1, 0.0), [0.0; 3]);
        assert_eq!(kernel_poly6_grad([0.0; 3], 0.0, h), [0.0; 3]);
    }

    // ── pressure / viscosity force tests ─────────────────────────────────

    #[test]
    fn test_pressure_force_finite() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 200.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        s.gpu_compute_density();
        s.gpu_compute_pressure();
        let f = s.gpu_pressure_force(0);
        assert!(f.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_viscosity_force_finite() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 200.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        s.gpu_compute_density();
        let vels = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let f = s.gpu_viscosity_force(0, &vels, 0.001);
        assert!(f.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_total_mass_empty() {
        let s = GpuSphPressureSolver::new(0, 1.0, 1000.0, 1.0);
        assert_eq!(s.total_mass(), 0.0);
    }

    #[test]
    fn test_solver_clone() {
        let s = GpuSphPressureSolver::new(3, 1.0, 1000.0, 1.0);
        let s2 = s.clone();
        assert_eq!(s2.particle_count(), 3);
    }

    #[test]
    fn test_stats_compression_error_nonzero() {
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 1.0);
        s.density[0] = 1100.0;
        let stats = s.compute_stats();
        assert!(stats.compression_error > 0.0);
    }
}
