// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU SPH pressure solver (CPU mock backend).
//!
//! Implements the Smoothed Particle Hydrodynamics pressure and viscosity
//! force computation pipeline using Poly6, Spiky, and Viscosity kernels.
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

/// PCISPH pressure correction loop.
///
/// Iteratively updates pressure until the density error drops below `tol`
/// or `max_iter` iterations are reached. Returns the number of iterations.
pub fn pcisph_gpu_correction(
    solver: &mut GpuSphPressureSolver,
    max_iter: usize,
    tol: f64,
) -> usize {
    for iter in 0..max_iter {
        solver.gpu_compute_density();
        solver.gpu_compute_pressure();
        if solver.density_error() < tol {
            return iter + 1;
        }
    }
    max_iter
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

    #[test]
    fn test_pcisph_returns_iterations_le_max() {
        let mut s = GpuSphPressureSolver::new(2, 1.0, 1000.0, 1.0);
        s.set_position(0, [0.0, 0.0, 0.0]);
        s.set_position(1, [0.3, 0.0, 0.0]);
        let iters = pcisph_gpu_correction(&mut s, 10, 1e-3);
        assert!(iters <= 10);
    }

    #[test]
    fn test_pcisph_updates_density() {
        let mut s = GpuSphPressureSolver::new(1, 1.0, 1000.0, 1.0);
        pcisph_gpu_correction(&mut s, 1, 1.0);
        assert!(s.density[0] > 0.0);
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
