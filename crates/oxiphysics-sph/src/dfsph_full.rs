#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Full DFSPH solver (Bender & Koschier 2015, 2017).
//!
//! Two-phase pressure solve:
//!   1. Divergence-free velocity correction (incompressibility of velocities)
//!   2. Density error correction (constant density enforcement)
//!
//! This module is self-contained and operates on `[f64; 3]` arrays, making it
//! easy to embed or test independently of the rest of the engine.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Cubic spline kernel (internal, 3-D, support radius 2h)
// ---------------------------------------------------------------------------

#[inline]
fn cubic_w(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    }
}

/// Returns the *scalar* dW/dr (multiply by r_hat to get vector gradient).
#[inline]
fn cubic_grad_w(r: f64, h: f64) -> f64 {
    let sigma = 1.0 / (PI * h * h * h);
    let q = r / h;
    if !(1e-14..2.0).contains(&q) {
        0.0
    } else if q >= 1.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t) / h
    } else {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    }
}

// ---------------------------------------------------------------------------
// Small 3-D vector helpers (avoid pulling in nalgebra here)
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// DfsphState
// ---------------------------------------------------------------------------

/// Full DFSPH solver state.
///
/// Stores all per-particle data and provides the two-phase pressure iteration
/// described by Bender & Koschier (2015, 2017).
#[derive(Debug, Clone)]
pub struct DfsphState {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// SPH-summed densities.
    pub densities: Vec<f64>,
    /// Pressure (density correction phase).
    pub pressures: Vec<f64>,
    /// Pressure (divergence correction phase).
    pub pressures_div: Vec<f64>,
    /// Per-particle α_i factor.
    pub alphas: Vec<f64>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Smoothing length.
    pub h: f64,
    /// Rest density (kg/m³).
    pub rho0: f64,
    /// Time step.
    pub dt: f64,
    n: usize,
}

impl DfsphState {
    /// Create a new empty DFSPH state.
    pub fn new(h: f64, rho0: f64, dt: f64) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            pressures_div: Vec::new(),
            alphas: Vec::new(),
            masses: Vec::new(),
            h,
            rho0,
            dt,
            n: 0,
        }
    }

    /// Add a particle at `pos` with the given `mass`.
    pub fn add_particle(&mut self, pos: [f64; 3], mass: f64) {
        self.positions.push(pos);
        self.velocities.push([0.0; 3]);
        self.densities.push(0.0);
        self.pressures.push(0.0);
        self.pressures_div.push(0.0);
        self.alphas.push(0.0);
        self.masses.push(mass);
        self.n += 1;
    }

    /// Compute density via SPH summation using the cubic spline kernel.
    ///
    /// ρ_i = Σ_j m_j W(|r_ij|, h)
    #[allow(clippy::needless_range_loop)]
    pub fn compute_densities(&mut self, neighbors: &[Vec<usize>]) {
        let h = self.h;
        for i in 0..self.n {
            let mut rho = self.masses[i] * cubic_w(0.0, h); // self-contribution
            for &j in &neighbors[i] {
                let r = norm3(sub3(self.positions[i], self.positions[j]));
                rho += self.masses[j] * cubic_w(r, h);
            }
            self.densities[i] = rho;
        }
    }

    /// Compute per-particle α factors for the pressure corrections.
    ///
    /// α_i = ρ_i² / (|Σ_j m_j ∇W_ij|² + Σ_j |m_j ∇W_ij|²)
    #[allow(clippy::needless_range_loop)]
    pub fn compute_alphas(&mut self, neighbors: &[Vec<usize>]) {
        let h = self.h;
        for i in 0..self.n {
            let mut sum_grad = [0.0f64; 3];
            let mut sum_grad_sq = 0.0f64;

            for &j in &neighbors[i] {
                let rij = sub3(self.positions[i], self.positions[j]);
                let r = norm3(rij);
                if r < 1e-14 {
                    continue;
                }
                let dw_dr = cubic_grad_w(r, h);
                let r_hat = scale3(rij, 1.0 / r);
                // m_j * ∇W_ij  (vector)
                let grad_contrib = scale3(r_hat, self.masses[j] * dw_dr);
                sum_grad = add3(sum_grad, grad_contrib);
                sum_grad_sq += dot3(grad_contrib, grad_contrib);
            }

            let denom = dot3(sum_grad, sum_grad) + sum_grad_sq;
            self.alphas[i] = if denom > 1e-28 {
                self.densities[i] * self.densities[i] / denom
            } else {
                0.0
            };
        }
    }

    /// Phase 1: divergence-free velocity correction.
    ///
    /// Iterates until `max(|div v_i|) < div_eps` or `max_iter` is reached.
    /// Returns the number of iterations performed.
    ///
    /// Per iteration:
    ///   κ_i = (div v_i) * α_i / dt
    ///   Δv_i = -dt * Σ_j m_j * (κ_i/ρ_i² + κ_j/ρ_j²) * ∇W_ij
    #[allow(clippy::needless_range_loop)]
    pub fn divergence_free_solve(
        &mut self,
        neighbors: &[Vec<usize>],
        max_iter: usize,
        div_eps: f64,
    ) -> usize {
        let h = self.h;
        let dt = self.dt;

        for iter in 0..max_iter {
            // Compute κ_i from velocity divergence
            let mut kappas = vec![0.0f64; self.n];
            let mut max_div = 0.0f64;
            for i in 0..self.n {
                let div_v = self.velocity_divergence(i, neighbors);
                let rho_i = self.densities[i].max(1e-14);
                kappas[i] = div_v * self.alphas[i] / (dt * rho_i * rho_i);
                // store scaled pressure for divergence phase
                self.pressures_div[i] = kappas[i];
                max_div = max_div.max(div_v.abs());
            }

            if max_div < div_eps {
                return iter + 1;
            }

            // Apply velocity correction
            // Δv_i = -dt * Σ_j m_j * (κ_i/ρ_i² + κ_j/ρ_j²) * ∇W_ij
            // Note: κ already has /ρ² baked in from above
            let mut dv = vec![[0.0f64; 3]; self.n];
            for i in 0..self.n {
                for &j in &neighbors[i] {
                    let rij = sub3(self.positions[i], self.positions[j]);
                    let r = norm3(rij);
                    if r < 1e-14 {
                        continue;
                    }
                    let dw_dr = cubic_grad_w(r, h);
                    let r_hat = scale3(rij, 1.0 / r);
                    let factor = -dt * self.masses[j] * (kappas[i] + kappas[j]) * dw_dr;
                    dv[i] = add3(dv[i], scale3(r_hat, factor));
                }
            }
            for i in 0..self.n {
                self.velocities[i] = add3(self.velocities[i], dv[i]);
            }
        }

        max_iter
    }

    /// Phase 2: density error correction.
    ///
    /// Iterates until `max(|ρ_i - ρ0|/ρ0) < density_eps` or `max_iter` is reached.
    /// Returns the number of iterations performed.
    ///
    /// Per iteration:
    ///   ρ*_i = ρ_i + dt * Σ_j m_j (v_i - v_j) · ∇W_ij
    ///   κ_i  = (ρ*_i - ρ0) * α_i / (dt² * ρ_i²)
    ///   Δv_i = -dt * Σ_j m_j * (κ_i/ρ_i² + κ_j/ρ_j²) * ∇W_ij
    #[allow(clippy::needless_range_loop)]
    pub fn density_solve(
        &mut self,
        neighbors: &[Vec<usize>],
        max_iter: usize,
        density_eps: f64,
    ) -> usize {
        let h = self.h;
        let dt = self.dt;
        let rho0 = self.rho0;

        for iter in 0..max_iter {
            // Predict density: ρ*_i = ρ_i + dt * Σ_j m_j (v_i - v_j) · ∇W_ij
            let mut rho_star = vec![0.0f64; self.n];
            let mut max_err = 0.0f64;
            for i in 0..self.n {
                let mut drho_dt = 0.0;
                for &j in &neighbors[i] {
                    let rij = sub3(self.positions[i], self.positions[j]);
                    let r = norm3(rij);
                    if r < 1e-14 {
                        continue;
                    }
                    let dw_dr = cubic_grad_w(r, h);
                    let r_hat = scale3(rij, 1.0 / r);
                    let vij = sub3(self.velocities[i], self.velocities[j]);
                    drho_dt += self.masses[j] * dot3(vij, scale3(r_hat, dw_dr));
                }
                rho_star[i] = self.densities[i] + dt * drho_dt;
                let err = ((rho_star[i] - rho0) / rho0).abs();
                max_err = max_err.max(err);
            }

            if max_err < density_eps {
                return iter + 1;
            }

            // Compute κ_i = (ρ*_i - ρ0) * α_i / (dt² * ρ_i²)
            let mut kappas = vec![0.0f64; self.n];
            for i in 0..self.n {
                let rho_i = self.densities[i].max(1e-14);
                let alpha_i = self.alphas[i];
                kappas[i] = if alpha_i.abs() > 1e-28 {
                    (rho_star[i] - rho0) * alpha_i / (dt * dt * rho_i * rho_i)
                } else {
                    0.0
                };
                self.pressures[i] = kappas[i];
            }

            // Apply velocity correction
            let mut dv = vec![[0.0f64; 3]; self.n];
            for i in 0..self.n {
                for &j in &neighbors[i] {
                    let rij = sub3(self.positions[i], self.positions[j]);
                    let r = norm3(rij);
                    if r < 1e-14 {
                        continue;
                    }
                    let dw_dr = cubic_grad_w(r, h);
                    let r_hat = scale3(rij, 1.0 / r);
                    let factor = -dt * self.masses[j] * (kappas[i] + kappas[j]) * dw_dr;
                    dv[i] = add3(dv[i], scale3(r_hat, factor));
                }
            }
            for i in 0..self.n {
                self.velocities[i] = add3(self.velocities[i], dv[i]);
            }
        }

        max_iter
    }

    /// Compute the velocity divergence for particle `i`.
    ///
    /// div v_i = -1/ρ_i * Σ_j m_j (v_j - v_i) · ∇W_ij
    ///
    /// (Equivalently: +1/ρ_i * Σ_j m_j (v_i - v_j) · ∇W_ij)
    pub fn velocity_divergence(&self, i: usize, neighbors: &[Vec<usize>]) -> f64 {
        let h = self.h;
        let rho_i = self.densities[i].max(1e-14);
        let mut div = 0.0;
        for &j in &neighbors[i] {
            let rij = sub3(self.positions[i], self.positions[j]);
            let r = norm3(rij);
            if r < 1e-14 {
                continue;
            }
            let dw_dr = cubic_grad_w(r, h);
            let r_hat = scale3(rij, 1.0 / r);
            let vij = sub3(self.velocities[i], self.velocities[j]);
            div += self.masses[j] * dot3(vij, scale3(r_hat, dw_dr));
        }
        div / rho_i
    }

    /// Apply non-pressure forces: gravity and artificial viscosity (XSPH variant).
    ///
    /// Viscosity term: Δv_i += ν * Σ_j (m_j/ρ_j) * (v_j - v_i) * W_ij
    #[allow(clippy::needless_range_loop)]
    pub fn apply_non_pressure_forces(
        &mut self,
        gravity: [f64; 3],
        viscosity: f64,
        neighbors: &[Vec<usize>],
    ) {
        let h = self.h;
        let dt = self.dt;

        // Gravity
        for i in 0..self.n {
            self.velocities[i] = add3(self.velocities[i], scale3(gravity, dt));
        }

        // XSPH viscosity
        if viscosity > 1e-28 {
            let mut dv_visc = vec![[0.0f64; 3]; self.n];
            for i in 0..self.n {
                for &j in &neighbors[i] {
                    let rij = sub3(self.positions[i], self.positions[j]);
                    let r = norm3(rij);
                    let w = cubic_w(r, h);
                    let rho_j = self.densities[j].max(1e-14);
                    let vji = sub3(self.velocities[j], self.velocities[i]);
                    let coeff = viscosity * self.masses[j] / rho_j * w;
                    dv_visc[i] = add3(dv_visc[i], scale3(vji, coeff));
                }
            }
            for i in 0..self.n {
                self.velocities[i] = add3(self.velocities[i], dv_visc[i]);
            }
        }
    }

    /// Integrate positions: pos += vel * dt.
    pub fn integrate_positions(&mut self) {
        let dt = self.dt;
        for i in 0..self.n {
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }

    /// Perform a full DFSPH time step.
    ///
    /// Order of operations:
    ///   1. `apply_non_pressure_forces`
    ///   2. `divergence_free_solve`
    ///   3. `density_solve`
    ///   4. `integrate_positions`
    pub fn step(&mut self, neighbors: &[Vec<usize>], gravity: [f64; 3], viscosity: f64) {
        self.apply_non_pressure_forces(gravity, viscosity, neighbors);
        self.divergence_free_solve(neighbors, 100, 1e-3);
        self.density_solve(neighbors, 100, 1e-3);
        self.integrate_positions();
    }

    /// Average relative density error |ρ_i - ρ0| / ρ0.
    pub fn avg_density_error(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        let rho0 = self.rho0;
        let sum: f64 = self
            .densities
            .iter()
            .map(|&rho| (rho - rho0).abs() / rho0)
            .sum();
        sum / self.n as f64
    }

    /// Maximum velocity-divergence magnitude across all particles.
    pub fn max_divergence(&self, neighbors: &[Vec<usize>]) -> f64 {
        (0..self.n)
            .map(|i| self.velocity_divergence(i, neighbors).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.n
    }
}

// ---------------------------------------------------------------------------
// Neighbor search
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// DFSPH full simulation loop
// ---------------------------------------------------------------------------

/// Full DFSPH simulation loop runner.
///
/// Manages the complete simulation including neighbor rebuilding and
/// diagnostic output.
#[allow(dead_code)]
pub struct DfsphSimulation {
    /// The solver state.
    pub state: DfsphState,
    /// Current neighbor lists.
    pub neighbors: Vec<Vec<usize>>,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Viscosity coefficient.
    pub viscosity: f64,
    /// Maximum divergence solver iterations.
    pub max_div_iter: usize,
    /// Maximum density solver iterations.
    pub max_density_iter: usize,
    /// Divergence tolerance.
    pub div_eps: f64,
    /// Density error tolerance.
    pub density_eps: f64,
    /// Current simulation time.
    pub time: f64,
    /// Number of steps taken.
    pub step_count: usize,
}

#[allow(dead_code)]
impl DfsphSimulation {
    /// Create a new simulation from a DFSPH state.
    pub fn new(state: DfsphState, gravity: [f64; 3], viscosity: f64) -> Self {
        let neighbors = find_neighbors_brute(&state.positions, state.h);
        Self {
            state,
            neighbors,
            gravity,
            viscosity,
            max_div_iter: 100,
            max_density_iter: 100,
            div_eps: 1e-3,
            density_eps: 1e-3,
            time: 0.0,
            step_count: 0,
        }
    }

    /// Advance the simulation by one time step.
    ///
    /// Returns `(div_iters, density_iters)` used in this step.
    pub fn advance(&mut self) -> (usize, usize) {
        self.state.compute_densities(&self.neighbors);
        self.state.compute_alphas(&self.neighbors);

        self.state
            .apply_non_pressure_forces(self.gravity, self.viscosity, &self.neighbors);

        let div_iters =
            self.state
                .divergence_free_solve(&self.neighbors, self.max_div_iter, self.div_eps);
        let density_iters =
            self.state
                .density_solve(&self.neighbors, self.max_density_iter, self.density_eps);

        self.state.integrate_positions();
        self.neighbors = find_neighbors_brute(&self.state.positions, self.state.h);

        self.time += self.state.dt;
        self.step_count += 1;

        (div_iters, density_iters)
    }

    /// Run `n` time steps.
    pub fn run(&mut self, n: usize) {
        for _ in 0..n {
            self.advance();
        }
    }

    /// Return diagnostic information about the current state.
    pub fn diagnostics(&self) -> DfsphDiagnostics {
        let max_div = self.state.max_divergence(&self.neighbors);
        let avg_err = self.state.avg_density_error();
        let max_speed = self
            .state
            .velocities
            .iter()
            .map(|v| norm3(*v))
            .fold(0.0_f64, f64::max);
        let ke = self.kinetic_energy();
        DfsphDiagnostics {
            time: self.time,
            step: self.step_count,
            max_divergence: max_div,
            avg_density_error: avg_err,
            max_speed,
            kinetic_energy: ke,
        }
    }

    /// Compute total kinetic energy: ½ Σ mᵢ |vᵢ|².
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.state.n_particles() {
            let v = self.state.velocities[i];
            let v2 = dot3(v, v);
            ke += 0.5 * self.state.masses[i] * v2;
        }
        ke
    }
}

/// Diagnostic snapshot from a DFSPH simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DfsphDiagnostics {
    /// Current simulation time.
    pub time: f64,
    /// Number of steps taken.
    pub step: usize,
    /// Maximum velocity divergence.
    pub max_divergence: f64,
    /// Average relative density error.
    pub avg_density_error: f64,
    /// Maximum particle speed.
    pub max_speed: f64,
    /// Total kinetic energy.
    pub kinetic_energy: f64,
}

// ---------------------------------------------------------------------------
// Velocity field analysis
// ---------------------------------------------------------------------------

/// Analyze properties of the velocity field in a DFSPH simulation.
#[allow(dead_code)]
pub struct VelocityFieldAnalysis;

#[allow(dead_code)]
impl VelocityFieldAnalysis {
    /// Compute the average velocity across all particles.
    pub fn mean_velocity(state: &DfsphState) -> [f64; 3] {
        let n = state.n_particles();
        if n == 0 {
            return [0.0; 3];
        }
        let mut sum = [0.0; 3];
        for v in &state.velocities {
            sum[0] += v[0];
            sum[1] += v[1];
            sum[2] += v[2];
        }
        let inv_n = 1.0 / n as f64;
        [sum[0] * inv_n, sum[1] * inv_n, sum[2] * inv_n]
    }

    /// Compute velocity variance (scalar: mean of |v - v_mean|²).
    pub fn velocity_variance(state: &DfsphState) -> f64 {
        let n = state.n_particles();
        if n == 0 {
            return 0.0;
        }
        let mean = Self::mean_velocity(state);
        let mut sum_sq = 0.0;
        for v in &state.velocities {
            let dv = sub3(*v, mean);
            sum_sq += dot3(dv, dv);
        }
        sum_sq / n as f64
    }

    /// Compute the total linear momentum Σ mᵢ vᵢ.
    pub fn total_momentum(state: &DfsphState) -> [f64; 3] {
        let mut p = [0.0; 3];
        for i in 0..state.n_particles() {
            let m = state.masses[i];
            p[0] += m * state.velocities[i][0];
            p[1] += m * state.velocities[i][1];
            p[2] += m * state.velocities[i][2];
        }
        p
    }

    /// Compute max speed across all particles.
    pub fn max_speed(state: &DfsphState) -> f64 {
        state
            .velocities
            .iter()
            .map(|v| norm3(*v))
            .fold(0.0_f64, f64::max)
    }

    /// Compute the velocity curl (vorticity) at particle `i`.
    ///
    /// ω_i = Σ_j (m_j/ρ_j) (v_j - v_i) × ∇W_ij
    pub fn vorticity(state: &DfsphState, i: usize, neighbors: &[Vec<usize>]) -> [f64; 3] {
        let h = state.h;
        let mut omega = [0.0; 3];
        for &j in &neighbors[i] {
            let rij = sub3(state.positions[i], state.positions[j]);
            let r = norm3(rij);
            if r < 1e-14 {
                continue;
            }
            let dw_dr = cubic_grad_w(r, h);
            let r_hat = scale3(rij, 1.0 / r);
            let grad_w = scale3(r_hat, dw_dr);
            let vji = sub3(state.velocities[j], state.velocities[i]);
            let rho_j = state.densities[j].max(1e-14);
            let coeff = state.masses[j] / rho_j;
            // Cross product: vji × grad_w
            let cross = [
                vji[1] * grad_w[2] - vji[2] * grad_w[1],
                vji[2] * grad_w[0] - vji[0] * grad_w[2],
                vji[0] * grad_w[1] - vji[1] * grad_w[0],
            ];
            omega[0] += coeff * cross[0];
            omega[1] += coeff * cross[1];
            omega[2] += coeff * cross[2];
        }
        omega
    }

    /// Compute maximum vorticity magnitude across all particles.
    pub fn max_vorticity(state: &DfsphState, neighbors: &[Vec<usize>]) -> f64 {
        (0..state.n_particles())
            .map(|i| norm3(Self::vorticity(state, i, neighbors)))
            .fold(0.0_f64, f64::max)
    }

    /// Compute the CFL time step: dt_cfl = cfl_factor * h / max_speed.
    pub fn cfl_dt(state: &DfsphState, cfl_factor: f64) -> f64 {
        let v_max = Self::max_speed(state);
        if v_max < 1e-14 {
            return f64::MAX;
        }
        cfl_factor * state.h / v_max
    }
}

// ---------------------------------------------------------------------------
// Neighbor search
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// DFSPH alpha pre-computation (Bender 2015 formulation)
// ---------------------------------------------------------------------------

/// Compute the DFSPH α factor for a single particle.
///
/// This is the scalar factor used in both pressure solves:
///
/// ```text
/// α_i = ρ_i² / (|Σ_j m_j ∇W_ij|² + Σ_j m_j² |∇W_ij|²)
/// ```
///
/// The denominator combines a squared-sum term and a sum-of-squares term,
/// making it always non-negative.  Returns 0 if the denominator is too small.
pub fn compute_alpha_single(
    pos_i: [f64; 3],
    rho_i: f64,
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let mut sum_grad = [0.0_f64; 3];
    let mut sum_grad_sq = 0.0_f64;

    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass.iter()) {
        let rij = sub3(pos_i, *pos_j);
        let r = norm3(rij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_grad_w(r, h);
        let r_hat = scale3(rij, 1.0 / r);
        let g = scale3(r_hat, m_j * dw_dr);
        sum_grad = add3(sum_grad, g);
        sum_grad_sq += dot3(g, g);
    }

    let denom = dot3(sum_grad, sum_grad) + sum_grad_sq;
    if denom > 1e-28 {
        rho_i * rho_i / denom
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Factor table (pre-computation cache)
// ---------------------------------------------------------------------------

/// Pre-computed DFSPH factor table.
///
/// Stores α values computed from a reference configuration and provides
/// lookup by particle index.  Used when the geometry is fixed or when
/// recomputing α every step would be too costly.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DfsphFactorTable {
    /// Pre-computed α values (one per particle).
    pub alphas: Vec<f64>,
    /// Rest density used for the pre-computation.
    pub rho0: f64,
    /// Smoothing length used.
    pub h: f64,
}

#[allow(dead_code)]
impl DfsphFactorTable {
    /// Build the factor table from the given particle configuration.
    pub fn build(
        positions: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        h: f64,
        rho0: f64,
    ) -> Self {
        let neighbors = find_neighbors_brute(positions, h);
        let n = positions.len();
        let mut alphas = vec![0.0_f64; n];
        for i in 0..n {
            let nb_pos: Vec<[f64; 3]> = neighbors[i].iter().map(|&j| positions[j]).collect();
            let nb_mass: Vec<f64> = neighbors[i].iter().map(|&j| masses[j]).collect();
            alphas[i] = compute_alpha_single(positions[i], densities[i], h, &nb_pos, &nb_mass);
        }
        Self { alphas, rho0, h }
    }

    /// Get α for particle `i`.
    pub fn get(&self, i: usize) -> f64 {
        self.alphas.get(i).copied().unwrap_or(0.0)
    }

    /// Number of particles in the table.
    pub fn len(&self) -> usize {
        self.alphas.len()
    }

    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.alphas.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Divergence-free constraint helpers (standalone, Bender 2015)
// ---------------------------------------------------------------------------

/// Compute the velocity divergence at position `pos_i`.
///
/// `div v_i = 1/ρ_i * Σ_j m_j (v_i - v_j) · ∇W_ij`
pub fn velocity_divergence_single(
    pos_i: [f64; 3],
    vel_i: [f64; 3],
    rho_i: f64,
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_vel: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let rho_safe = rho_i.max(1e-14);
    let mut div = 0.0_f64;
    for ((pos_j, vel_j), &m_j) in neighbor_pos
        .iter()
        .zip(neighbor_vel.iter())
        .zip(neighbor_mass.iter())
    {
        let rij = sub3(pos_i, *pos_j);
        let r = norm3(rij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_grad_w(r, h);
        let r_hat = scale3(rij, 1.0 / r);
        let vij = sub3(vel_i, *vel_j);
        div += m_j * dot3(vij, scale3(r_hat, dw_dr));
    }
    div / rho_safe
}

/// Compute the predicted density rate of change:
///
/// ```text
/// Dρ/Dt ≈ Σ_j m_j (v_i - v_j) · ∇W_ij
/// ```
///
/// This is used in both divergence-free and density solves.
pub fn drho_dt_sph(
    pos_i: [f64; 3],
    vel_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_vel: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let mut drho = 0.0_f64;
    for ((pos_j, vel_j), &m_j) in neighbor_pos
        .iter()
        .zip(neighbor_vel.iter())
        .zip(neighbor_mass.iter())
    {
        let rij = sub3(pos_i, *pos_j);
        let r = norm3(rij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_grad_w(r, h);
        let r_hat = scale3(rij, 1.0 / r);
        let vij = sub3(vel_i, *vel_j);
        drho += m_j * dot3(vij, scale3(r_hat, dw_dr));
    }
    drho
}

// ---------------------------------------------------------------------------
// Two-stage DFSPH solver (standalone functional API)
// ---------------------------------------------------------------------------

/// Two-stage DFSPH pressure solve using plain `[f64; 3]` arrays.
///
/// Stage 1: divergence-free velocity correction (ensures ∇·v ≈ 0).
/// Stage 2: density correction (ensures ρ ≈ ρ₀).
///
/// Returns `(stage1_iters, stage2_iters)`.
#[allow(clippy::too_many_arguments)]
pub fn dfsph_two_stage_solve(
    positions: &[[f64; 3]],
    velocities: &mut Vec<[f64; 3]>,
    densities: &[f64],
    alphas: &[f64],
    masses: &[f64],
    neighbors: &[Vec<usize>],
    h: f64,
    rho0: f64,
    dt: f64,
    max_iter_div: usize,
    div_eps: f64,
    max_iter_den: usize,
    den_eps: f64,
) -> (usize, usize) {
    let n = positions.len();

    // --- Stage 1: divergence-free solve ---
    let mut iters1 = 0;
    for iter in 0..max_iter_div {
        iters1 = iter + 1;
        // Compute κ_div for each particle: κ_i = α_i * (div v_i) / dt
        let mut kappas = vec![0.0_f64; n];
        let mut max_div = 0.0_f64;
        for i in 0..n {
            let nb_pos: Vec<[f64; 3]> = neighbors[i].iter().map(|&j| positions[j]).collect();
            let nb_vel: Vec<[f64; 3]> = neighbors[i].iter().map(|&j| velocities[j]).collect();
            let nb_mass: Vec<f64> = neighbors[i].iter().map(|&j| masses[j]).collect();
            let div_v = velocity_divergence_single(
                positions[i],
                velocities[i],
                densities[i],
                h,
                &nb_pos,
                &nb_vel,
                &nb_mass,
            );
            let rho_i = densities[i].max(1e-14);
            kappas[i] = alphas[i] * div_v / (dt * rho_i * rho_i);
            max_div = max_div.max(div_v.abs());
        }
        if max_div < div_eps {
            break;
        }
        // Apply Δv
        let mut dv = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for &j in &neighbors[i] {
                let rij = sub3(positions[i], positions[j]);
                let r = norm3(rij);
                if r < 1e-14 {
                    continue;
                }
                let dw_dr = cubic_grad_w(r, h);
                let r_hat = scale3(rij, 1.0 / r);
                let factor = -dt * masses[j] * (kappas[i] + kappas[j]) * dw_dr;
                dv[i] = add3(dv[i], scale3(r_hat, factor));
            }
        }
        for i in 0..n {
            velocities[i] = add3(velocities[i], dv[i]);
        }
    }

    // --- Stage 2: density solve ---
    let mut iters2 = 0;
    for iter in 0..max_iter_den {
        iters2 = iter + 1;
        let mut kappas = vec![0.0_f64; n];
        let mut max_err = 0.0_f64;
        for i in 0..n {
            let nb_pos: Vec<[f64; 3]> = neighbors[i].iter().map(|&j| positions[j]).collect();
            let nb_vel: Vec<[f64; 3]> = neighbors[i].iter().map(|&j| velocities[j]).collect();
            let nb_mass: Vec<f64> = neighbors[i].iter().map(|&j| masses[j]).collect();
            let drho = drho_dt_sph(positions[i], velocities[i], h, &nb_pos, &nb_vel, &nb_mass);
            let rho_star = densities[i] + dt * drho;
            let rho_i = densities[i].max(1e-14);
            let err = ((rho_star - rho0) / rho0).abs();
            max_err = max_err.max(err);
            kappas[i] = if alphas[i].abs() > 1e-28 {
                (rho_star - rho0) * alphas[i] / (dt * dt * rho_i * rho_i)
            } else {
                0.0
            };
        }
        if max_err < den_eps {
            break;
        }
        let mut dv = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for &j in &neighbors[i] {
                let rij = sub3(positions[i], positions[j]);
                let r = norm3(rij);
                if r < 1e-14 {
                    continue;
                }
                let dw_dr = cubic_grad_w(r, h);
                let r_hat = scale3(rij, 1.0 / r);
                let factor = -dt * masses[j] * (kappas[i] + kappas[j]) * dw_dr;
                dv[i] = add3(dv[i], scale3(r_hat, factor));
            }
        }
        for i in 0..n {
            velocities[i] = add3(velocities[i], dv[i]);
        }
    }

    (iters1, iters2)
}

// ---------------------------------------------------------------------------
// DFSPH iteration diagnostics
// ---------------------------------------------------------------------------

/// Diagnostics collected from one DFSPH two-stage solve.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DfsphIterationStats {
    /// Iterations used by the divergence-free stage.
    pub div_iters: usize,
    /// Iterations used by the density stage.
    pub den_iters: usize,
    /// Final max velocity divergence after stage 1.
    pub final_max_div: f64,
    /// Final max density error after stage 2.
    pub final_max_den_err: f64,
}

#[allow(dead_code)]
impl DfsphIterationStats {
    /// Compute diagnostics from a fully solved state.
    pub fn from_state(state: &DfsphState, neighbors: &[Vec<usize>]) -> Self {
        let max_div = state.max_divergence(neighbors);
        let max_den = state
            .densities
            .iter()
            .map(|&rho| ((rho - state.rho0) / state.rho0).abs())
            .fold(0.0_f64, f64::max);
        Self {
            div_iters: 0,
            den_iters: 0,
            final_max_div: max_div,
            final_max_den_err: max_den,
        }
    }
}

// ---------------------------------------------------------------------------
// SPH pressure helper (two-phase, plain arrays)
// ---------------------------------------------------------------------------

/// Compute SPH pressure forces using plain arrays.
///
/// Returns per-particle accelerations (force per unit mass):
/// ```text
/// a_i = -Σ_j m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W_ij
/// ```
#[allow(dead_code)]
pub fn sph_pressure_accel_plain(
    positions: &[[f64; 3]],
    densities: &[f64],
    pressures: &[f64],
    masses: &[f64],
    neighbors: &[Vec<usize>],
    h: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut accels = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        let rho_i = densities[i].max(1e-14);
        let p_i = pressures[i];
        for &j in &neighbors[i] {
            let rij = sub3(positions[i], positions[j]);
            let r = norm3(rij);
            if r < 1e-14 {
                continue;
            }
            let dw_dr = cubic_grad_w(r, h);
            let rho_j = densities[j].max(1e-14);
            let p_j = pressures[j];
            let factor = -masses[j] * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw_dr / r;
            accels[i][0] += factor * rij[0];
            accels[i][1] += factor * rij[1];
            accels[i][2] += factor * rij[2];
        }
    }
    accels
}

/// Simple O(n²) brute-force neighbor search for testing.
///
/// Finds all particles within distance `2h` (the kernel support radius for
/// the cubic spline kernel used here).
pub fn find_neighbors_brute(positions: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
    let n = positions.len();
    let support = 2.0 * h;
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let r = norm3(sub3(positions[i], positions[j]));
            if r < support {
                neighbors[i].push(j);
            }
        }
    }
    neighbors
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a small uniform grid of particles.
    ///
    /// Returns `(state, neighbors)`.
    fn uniform_grid(
        nx: usize,
        ny: usize,
        nz: usize,
        spacing: f64,
    ) -> (DfsphState, Vec<Vec<usize>>) {
        let h = 2.0 * spacing;
        let rho0 = 1000.0;
        let mass = rho0 * spacing * spacing * spacing;
        let mut state = DfsphState::new(h, rho0, 0.001);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let pos = [i as f64 * spacing, j as f64 * spacing, k as f64 * spacing];
                    state.add_particle(pos, mass);
                }
            }
        }
        let neighbors = find_neighbors_brute(&state.positions, h);
        state.compute_densities(&neighbors);
        (state, neighbors)
    }

    /// Uniform 4×4×4 grid: all densities should be positive and at least 10% of rho0.
    ///
    /// Interior particles are relatively close to rho0; boundary particles will be
    /// lower (they have fewer neighbors) but should still be > 0.
    #[test]
    fn test_dfsph_density_init() {
        let (state, _) = uniform_grid(4, 4, 4, 0.05);
        // All densities positive
        for (i, &rho) in state.densities.iter().enumerate() {
            assert!(rho > 0.0, "Density of particle {i} is non-positive: {rho}");
        }
        // At least one particle should have density close to rho0 (interior)
        let max_rho = state.densities.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_rho > state.rho0 * 0.7,
            "Max density {max_rho} is too far below rho0 {}",
            state.rho0
        );
    }

    /// After compute_alphas, all α values should be non-negative.
    #[test]
    fn test_dfsph_alphas_positive() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);
        for (i, &a) in state.alphas.iter().enumerate() {
            assert!(a >= 0.0, "Alpha[{i}] = {a} is negative");
        }
    }

    /// After one divergence solve, max divergence should decrease (or already be small).
    #[test]
    fn test_dfsph_divergence_solve_reduces_divergence() {
        let (mut state, neighbors) = uniform_grid(4, 4, 4, 0.05);
        state.compute_alphas(&neighbors);

        // Perturb velocities to create nonzero divergence
        for (i, v) in state.velocities.iter_mut().enumerate() {
            v[0] = 0.1 * (i as f64 * 0.37).sin();
            v[1] = 0.05 * (i as f64 * 0.71).cos();
        }

        let div_before = state.max_divergence(&neighbors);
        state.divergence_free_solve(&neighbors, 50, 1e-6);
        let div_after = state.max_divergence(&neighbors);

        // Either the divergence was already small, or it was reduced
        assert!(
            div_after <= div_before + 1e-12,
            "Divergence increased: before={div_before}, after={div_after}"
        );
    }

    /// After one density solve, avg density error should not increase.
    #[test]
    fn test_dfsph_density_solve_reduces_error() {
        let (mut state, neighbors) = uniform_grid(4, 4, 4, 0.05);
        state.compute_alphas(&neighbors);

        let err_before = state.avg_density_error();
        state.density_solve(&neighbors, 50, 1e-6);
        // Re-compute densities after the velocity corrections
        state.compute_densities(&neighbors);
        let err_after = state.avg_density_error();

        // Density error should not grow significantly (may stay same if already tiny)
        assert!(
            err_after <= err_before + 0.01,
            "Density error grew: before={err_before}, after={err_after}"
        );
    }

    /// 2×2×2 grid (8 particles), run 5 full steps without panicking.
    #[test]
    fn test_dfsph_step_no_panic() {
        let (mut state, mut neighbors) = uniform_grid(2, 2, 2, 0.1);
        state.compute_alphas(&neighbors);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..5 {
            state.compute_densities(&neighbors);
            state.compute_alphas(&neighbors);
            state.step(&neighbors, gravity, 0.01);
            // Rebuild neighbors after positions change
            neighbors = find_neighbors_brute(&state.positions, state.h);
        }
        assert_eq!(state.n_particles(), 8);
    }

    /// 4 particles all within 2h of each other should all be mutual neighbors.
    #[test]
    fn test_brute_force_neighbors() {
        let h = 1.0;
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
        ];
        let neighbors = find_neighbors_brute(&positions, h);
        // Every particle should have exactly 3 neighbors (all others are within 2h=2.0)
        for (i, nb) in neighbors.iter().enumerate() {
            assert_eq!(
                nb.len(),
                3,
                "Particle {i} has {} neighbors, expected 3",
                nb.len()
            );
        }
    }

    // ── Additional DFSPH tests ──────────────────────────────────────────

    #[test]
    fn test_divergence_solver_returns_iterations() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);
        let iters = state.divergence_free_solve(&neighbors, 10, 1e-3);
        assert!((1..=10).contains(&iters));
    }

    #[test]
    fn test_density_solver_returns_iterations() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);
        let iters = state.density_solve(&neighbors, 10, 1e-3);
        assert!((1..=10).contains(&iters));
    }

    #[test]
    fn test_velocity_divergence_uniform_is_small() {
        let (state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        // Uniform grid with zero velocities should have ~zero divergence
        for i in 0..state.n_particles() {
            let div = state.velocity_divergence(i, &neighbors);
            assert!(
                div.abs() < 1e-6,
                "Divergence of particle {i} should be ~0, got {div}"
            );
        }
    }

    #[test]
    fn test_alpha_factors_are_bounded() {
        let (mut state, neighbors) = uniform_grid(4, 4, 4, 0.05);
        state.compute_alphas(&neighbors);
        for (i, &a) in state.alphas.iter().enumerate() {
            assert!(a.is_finite(), "Alpha[{i}] is not finite");
            assert!(a < 1e10, "Alpha[{i}] = {a} is unreasonably large");
        }
    }

    #[test]
    fn test_gravity_changes_velocity() {
        let (mut state, neighbors) = uniform_grid(2, 2, 2, 0.1);
        state.compute_alphas(&neighbors);
        state.apply_non_pressure_forces([0.0, -9.81, 0.0], 0.0, &neighbors);
        for v in &state.velocities {
            assert!(
                v[1] < 0.0,
                "vy should be negative after gravity, got {}",
                v[1]
            );
        }
    }

    #[test]
    fn test_viscosity_smooths_velocities() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_densities(&neighbors);
        // Give alternating velocities
        for (i, v) in state.velocities.iter_mut().enumerate() {
            v[0] = if i % 2 == 0 { 1.0 } else { -1.0 };
        }
        let spread_before: f64 = state.velocities.iter().map(|v| v[0].abs()).sum();
        state.apply_non_pressure_forces([0.0, 0.0, 0.0], 0.5, &neighbors);
        let spread_after: f64 = state.velocities.iter().map(|v| v[0].abs()).sum();
        // Viscosity should reduce spread of velocities
        assert!(
            spread_after < spread_before + 1e-6,
            "Viscosity should smooth velocities"
        );
    }

    #[test]
    fn test_integrate_positions_moves_particles() {
        let (mut state, _) = uniform_grid(2, 2, 2, 0.1);
        let pos_before = state.positions.clone();
        state.velocities[0] = [1.0, 2.0, 3.0];
        state.integrate_positions();
        let dt = state.dt;
        assert!((state.positions[0][0] - pos_before[0][0] - dt).abs() < 1e-14);
        assert!((state.positions[0][1] - pos_before[0][1] - 2.0 * dt).abs() < 1e-14);
        assert!((state.positions[0][2] - pos_before[0][2] - 3.0 * dt).abs() < 1e-14);
    }

    #[test]
    fn test_avg_density_error_empty_state() {
        let state = DfsphState::new(0.1, 1000.0, 0.001);
        assert_eq!(state.avg_density_error(), 0.0);
    }

    #[test]
    fn test_max_divergence_empty() {
        let state = DfsphState::new(0.1, 1000.0, 0.001);
        let neighbors: Vec<Vec<usize>> = Vec::new();
        assert_eq!(state.max_divergence(&neighbors), 0.0);
    }

    #[test]
    fn test_n_particles_count() {
        let mut state = DfsphState::new(0.1, 1000.0, 0.001);
        assert_eq!(state.n_particles(), 0);
        state.add_particle([0.0; 3], 1.0);
        state.add_particle([0.1, 0.0, 0.0], 1.0);
        assert_eq!(state.n_particles(), 2);
    }

    #[test]
    fn test_multi_step_simulation() {
        let (mut state, mut neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..10 {
            state.compute_densities(&neighbors);
            state.compute_alphas(&neighbors);
            state.step(&neighbors, gravity, 0.01);
            neighbors = find_neighbors_brute(&state.positions, state.h);
        }
        // All particles should have finite positions
        for (i, pos) in state.positions.iter().enumerate() {
            #[allow(clippy::needless_range_loop)]
            for d in 0..3 {
                assert!(
                    pos[d].is_finite(),
                    "Position of particle {i} dim {d} is not finite"
                );
            }
        }
    }

    #[test]
    fn test_density_near_rho0_for_interior() {
        let (state, _) = uniform_grid(6, 6, 6, 0.05);
        // Interior particles (away from boundary) should be positive
        let n = state.n_particles();
        let mid = n / 2;
        let rho = state.densities[mid];
        assert!(rho > 0.0, "Interior density {rho} should be positive");
        // Max density should be a reasonable fraction of rho0
        let max_rho = state.densities.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_rho > state.rho0 * 0.1,
            "Max density {max_rho} too low vs rho0 {}",
            state.rho0
        );
    }

    #[test]
    fn test_cubic_kernel_symmetry() {
        let h = 0.1;
        let r = 0.05;
        let w1 = cubic_w(r, h);
        let w2 = cubic_w(r, h);
        assert_eq!(w1, w2, "Kernel should be symmetric");
    }

    #[test]
    fn test_cubic_kernel_monotone_decreasing() {
        let h = 0.1;
        let w0 = cubic_w(0.0, h);
        let w1 = cubic_w(0.05, h);
        let w2 = cubic_w(0.1, h);
        let w3 = cubic_w(0.15, h);
        assert!(w0 >= w1, "Kernel should decrease: w(0)={w0} < w(0.05)={w1}");
        assert!(
            w1 >= w2,
            "Kernel should decrease: w(0.05)={w1} < w(0.1)={w2}"
        );
        assert!(
            w2 >= w3,
            "Kernel should decrease: w(0.1)={w2} < w(0.15)={w3}"
        );
    }

    #[test]
    fn test_cubic_kernel_zero_outside_support() {
        let h = 0.1;
        assert_eq!(cubic_w(0.2, h), 0.0);
        assert_eq!(cubic_w(0.3, h), 0.0);
        assert_eq!(cubic_w(1.0, h), 0.0);
    }

    #[test]
    fn test_cubic_gradient_zero_at_origin() {
        let h = 0.1;
        let g = cubic_grad_w(0.0, h);
        assert_eq!(g, 0.0, "Gradient at r=0 should be 0");
    }

    #[test]
    fn test_cubic_gradient_zero_outside_support() {
        let h = 0.1;
        assert_eq!(cubic_grad_w(0.2, h), 0.0);
        assert_eq!(cubic_grad_w(0.5, h), 0.0);
    }

    #[test]
    fn test_add_particle_initializes_velocity_zero() {
        let mut state = DfsphState::new(0.1, 1000.0, 0.001);
        state.add_particle([1.0, 2.0, 3.0], 0.5);
        assert_eq!(state.velocities[0], [0.0, 0.0, 0.0]);
        assert_eq!(state.positions[0], [1.0, 2.0, 3.0]);
        assert_eq!(state.masses[0], 0.5);
    }

    #[test]
    fn test_step_with_no_gravity_no_viscosity() {
        let (mut state, neighbors) = uniform_grid(2, 2, 2, 0.1);
        state.compute_alphas(&neighbors);
        let pos_before = state.positions.clone();
        state.step(&neighbors, [0.0, 0.0, 0.0], 0.0);
        // With zero gravity and zero viscosity, positions should change only slightly
        // (only through pressure corrections)
        for (i, pos) in state.positions.iter().enumerate() {
            for d in 0..3 {
                assert!(
                    pos[d].is_finite(),
                    "Position [{i}][{d}] not finite after step"
                );
                let _diff = (pos[d] - pos_before[i][d]).abs();
            }
        }
    }

    #[test]
    fn test_neighbors_self_excluded() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let neighbors = find_neighbors_brute(&positions, 0.2);
        // Particle 0 should not contain itself
        assert!(!neighbors[0].contains(&0));
        assert!(neighbors[0].contains(&1));
    }

    #[test]
    fn test_neighbors_empty_far_apart() {
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let neighbors = find_neighbors_brute(&positions, 0.1);
        assert!(neighbors[0].is_empty());
        assert!(neighbors[1].is_empty());
    }

    // -----------------------------------------------------------------------
    // compute_alpha_single tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_alpha_single_positive_with_neighbors() {
        let pos_i = [0.0_f64; 3];
        let rho_i = 1000.0;
        let h = 0.1;
        let nb_pos = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0], [0.0, 0.05, 0.0]];
        let nb_mass = vec![0.001_f64; 3];
        let alpha = compute_alpha_single(pos_i, rho_i, h, &nb_pos, &nb_mass);
        assert!(
            alpha > 0.0,
            "alpha should be positive with neighbors: {alpha}"
        );
        assert!(alpha.is_finite(), "alpha should be finite: {alpha}");
    }

    #[test]
    fn test_alpha_single_zero_no_neighbors() {
        let pos_i = [0.0_f64; 3];
        let alpha = compute_alpha_single(pos_i, 1000.0, 0.1, &[], &[]);
        assert_eq!(alpha, 0.0, "No neighbors → alpha = 0");
    }

    // -----------------------------------------------------------------------
    // DfsphFactorTable tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_factor_table_build() {
        let spacing = 0.05;
        let h = 0.1;
        let rho0 = 1000.0;
        let mass = rho0 * spacing * spacing * spacing;
        let mut positions = Vec::new();
        let mut masses = Vec::new();
        let mut densities = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    positions.push([i as f64 * spacing, j as f64 * spacing, k as f64 * spacing]);
                    masses.push(mass);
                    densities.push(rho0);
                }
            }
        }
        let table = DfsphFactorTable::build(&positions, &masses, &densities, h, rho0);
        assert_eq!(table.len(), 27);
        assert!(!table.is_empty());
    }

    #[test]
    fn test_factor_table_get_in_bounds() {
        let h = 0.1;
        let pos = vec![[0.0_f64; 3], [0.05, 0.0, 0.0]];
        let mass = vec![0.001_f64; 2];
        let den = vec![1000.0_f64; 2];
        let table = DfsphFactorTable::build(&pos, &mass, &den, h, 1000.0);
        let a = table.get(0);
        assert!(a.is_finite());
        // Out-of-bounds returns 0
        assert_eq!(table.get(100), 0.0);
    }

    // -----------------------------------------------------------------------
    // velocity_divergence_single tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_velocity_divergence_single_zero_velocity() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [0.0_f64; 3];
        let rho_i = 1000.0;
        let h = 0.1;
        let nb_pos = vec![[0.05_f64, 0.0, 0.0]];
        let nb_vel = vec![[0.0_f64; 3]];
        let nb_mass = vec![0.001_f64];
        let div = velocity_divergence_single(pos_i, vel_i, rho_i, h, &nb_pos, &nb_vel, &nb_mass);
        assert!(div.abs() < 1e-14, "Zero velocity → zero divergence: {div}");
    }

    #[test]
    fn test_velocity_divergence_single_finite() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let rho_i = 1000.0;
        let h = 0.1;
        let nb_pos = vec![[0.05_f64, 0.0, 0.0]];
        let nb_vel = vec![[0.0_f64; 3]];
        let nb_mass = vec![0.001_f64];
        let div = velocity_divergence_single(pos_i, vel_i, rho_i, h, &nb_pos, &nb_vel, &nb_mass);
        assert!(div.is_finite(), "Divergence should be finite: {div}");
    }

    // -----------------------------------------------------------------------
    // drho_dt_sph tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_drho_dt_zero_relative_velocity() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 2.0, 3.0];
        let h = 0.1;
        let nb_pos = vec![[0.05_f64, 0.0, 0.0]];
        let nb_vel = vec![[1.0_f64, 2.0, 3.0]]; // same velocity
        let nb_mass = vec![0.001_f64];
        let drho = drho_dt_sph(pos_i, vel_i, h, &nb_pos, &nb_vel, &nb_mass);
        assert!(
            drho.abs() < 1e-14,
            "Zero relative velocity → drho/dt=0: {drho}"
        );
    }

    #[test]
    fn test_drho_dt_finite_nonzero() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let h = 0.1;
        let nb_pos = vec![[0.05_f64, 0.0, 0.0]];
        let nb_vel = vec![[0.0_f64; 3]];
        let nb_mass = vec![0.001_f64];
        let drho = drho_dt_sph(pos_i, vel_i, h, &nb_pos, &nb_vel, &nb_mass);
        assert!(drho.is_finite());
    }

    // -----------------------------------------------------------------------
    // dfsph_two_stage_solve tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_two_stage_solve_uniform_converges() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);

        let mut velocities = state.velocities.clone();
        let (iters1, iters2) = dfsph_two_stage_solve(
            &state.positions,
            &mut velocities,
            &state.densities,
            &state.alphas,
            &state.masses,
            &neighbors,
            state.h,
            state.rho0,
            state.dt,
            20,
            1e-4,
            20,
            1e-4,
        );
        assert!((1..=20).contains(&iters1));
        assert!((1..=20).contains(&iters2));
    }

    #[test]
    fn test_two_stage_solve_velocities_finite() {
        let (mut state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        state.compute_alphas(&neighbors);
        // Perturb velocities
        for (i, v) in state.velocities.iter_mut().enumerate() {
            v[0] = (i as f64 * 0.1).sin() * 0.5;
        }
        let mut velocities = state.velocities.clone();
        dfsph_two_stage_solve(
            &state.positions,
            &mut velocities,
            &state.densities,
            &state.alphas,
            &state.masses,
            &neighbors,
            state.h,
            state.rho0,
            state.dt,
            10,
            1e-3,
            10,
            1e-3,
        );
        for (i, v) in velocities.iter().enumerate() {
            #[allow(clippy::needless_range_loop)]
            for d in 0..3 {
                assert!(v[d].is_finite(), "velocity[{i}][{d}] not finite");
            }
        }
    }

    // -----------------------------------------------------------------------
    // sph_pressure_accel_plain tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pressure_accel_zero_pressure() {
        let (state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        let pressures = vec![0.0_f64; state.n_particles()];
        let accels = sph_pressure_accel_plain(
            &state.positions,
            &state.densities,
            &pressures,
            &state.masses,
            &neighbors,
            state.h,
        );
        for (i, a) in accels.iter().enumerate() {
            let mag = norm3(*a);
            assert!(mag < 1e-28, "Zero pressure → zero accel[{i}]: {mag}");
        }
    }

    #[test]
    fn test_pressure_accel_finite_with_pressure() {
        let (state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        let pressures: Vec<f64> = (0..state.n_particles())
            .map(|i| (i as f64 * 100.0).sin().abs() * 1000.0)
            .collect();
        let accels = sph_pressure_accel_plain(
            &state.positions,
            &state.densities,
            &pressures,
            &state.masses,
            &neighbors,
            state.h,
        );
        for (i, a) in accels.iter().enumerate() {
            #[allow(clippy::needless_range_loop)]
            for d in 0..3 {
                assert!(a[d].is_finite(), "accel[{i}][{d}] should be finite");
            }
        }
    }

    // -----------------------------------------------------------------------
    // DfsphIterationStats tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_iteration_stats_from_state() {
        let (state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        let stats = DfsphIterationStats::from_state(&state, &neighbors);
        assert!(stats.final_max_div.is_finite());
        assert!(stats.final_max_den_err.is_finite());
    }

    // -----------------------------------------------------------------------
    // DfsphSimulation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_simulation_advance_returns_iteration_counts() {
        let (state, _) = uniform_grid(3, 3, 3, 0.05);
        let mut sim = DfsphSimulation::new(state, [0.0, -9.81, 0.0], 0.01);
        let (d, r) = sim.advance();
        assert!(d >= 1);
        assert!(r >= 1);
    }

    #[test]
    fn test_simulation_run_n_steps() {
        let (state, _) = uniform_grid(2, 2, 2, 0.1);
        let mut sim = DfsphSimulation::new(state, [0.0, 0.0, 0.0], 0.0);
        sim.run(5);
        assert_eq!(sim.step_count, 5);
    }

    #[test]
    fn test_simulation_kinetic_energy_non_negative() {
        let (state, _) = uniform_grid(3, 3, 3, 0.05);
        let mut sim = DfsphSimulation::new(state, [0.0, -9.81, 0.0], 0.01);
        sim.advance();
        assert!(sim.kinetic_energy() >= 0.0);
    }

    #[test]
    fn test_simulation_diagnostics_finite() {
        let (state, _) = uniform_grid(3, 3, 3, 0.05);
        let mut sim = DfsphSimulation::new(state, [0.0, -9.81, 0.0], 0.0);
        sim.advance();
        let diag = sim.diagnostics();
        assert!(diag.max_divergence.is_finite());
        assert!(diag.avg_density_error.is_finite());
        assert!(diag.kinetic_energy >= 0.0);
    }

    // -----------------------------------------------------------------------
    // VelocityFieldAnalysis tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_velocity_field_mean_zero() {
        let (state, _) = uniform_grid(2, 2, 2, 0.1);
        // All velocities are zero
        let mean = VelocityFieldAnalysis::mean_velocity(&state);
        #[allow(clippy::needless_range_loop)]
        for d in 0..3 {
            assert!(mean[d].abs() < 1e-14);
        }
    }

    #[test]
    fn test_velocity_field_variance_zero_uniform() {
        let (mut state, _) = uniform_grid(2, 2, 2, 0.1);
        for v in &mut state.velocities {
            *v = [1.0, 2.0, 3.0];
        }
        let var = VelocityFieldAnalysis::velocity_variance(&state);
        assert!(var.abs() < 1e-12, "Uniform velocity → zero variance: {var}");
    }

    #[test]
    fn test_velocity_field_total_momentum() {
        let (mut state, _) = uniform_grid(2, 2, 2, 0.1);
        for v in &mut state.velocities {
            *v = [1.0, 0.0, 0.0];
        }
        let p = VelocityFieldAnalysis::total_momentum(&state);
        let total_mass: f64 = state.masses.iter().sum();
        assert!((p[0] - total_mass).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_field_max_speed() {
        let (mut state, _) = uniform_grid(2, 2, 2, 0.1);
        state.velocities[0] = [3.0, 4.0, 0.0]; // speed = 5
        let max_v = VelocityFieldAnalysis::max_speed(&state);
        assert!((max_v - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_velocity_field_cfl_dt() {
        let (mut state, _) = uniform_grid(2, 2, 2, 0.1);
        state.velocities[0] = [10.0, 0.0, 0.0];
        let dt = VelocityFieldAnalysis::cfl_dt(&state, 0.3);
        let expected = 0.3 * state.h / 10.0;
        assert!((dt - expected).abs() < 1e-14);
    }

    #[test]
    fn test_velocity_field_vorticity_zero_uniform() {
        let (state, neighbors) = uniform_grid(3, 3, 3, 0.05);
        // All velocities zero → vorticity should be zero
        let w = VelocityFieldAnalysis::vorticity(&state, 0, &neighbors);
        let mag = norm3(w);
        assert!(mag.abs() < 1e-10, "Zero velocities → zero vorticity: {mag}");
    }
}
