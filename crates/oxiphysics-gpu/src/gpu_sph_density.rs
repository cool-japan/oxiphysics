#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated SPH density computation (CPU mock implementation).
//!
//! This module provides Smoothed Particle Hydrodynamics (SPH) density
//! computation routines that mirror what would run on a GPU. All kernels
//! are executed on the CPU via plain loops for portability.

// ── SPH constants ────────────────────────────────────────────────────────────

/// Default Tait EOS exponent (gamma).
const TAIT_GAMMA: f64 = 7.0;

// ── Data structures ──────────────────────────────────────────────────────────

/// A grid of SPH particles with associated physical quantities.
///
/// Positions are stored as flat arrays of (x, y, z) triplets.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GpuSphGrid {
    /// Particle positions: `[x0, y0, z0, x1, y1, z1, …]`.
    pub positions: Vec<f64>,
    /// Particle masses (one per particle).
    pub masses: Vec<f64>,
    /// Smoothing lengths (one per particle).
    pub smoothing_lengths: Vec<f64>,
    /// Particle densities (one per particle, updated by [`gpu_density_kernel`]).
    pub densities: Vec<f64>,
    /// Particle pressures (one per particle, updated by [`gpu_pressure_tait`]).
    pub pressures: Vec<f64>,
    /// Particle velocities: `[vx0, vy0, vz0, …]`.
    pub velocities: Vec<f64>,
    /// Particle forces: `[fx0, fy0, fz0, …]` (output of [`gpu_force_kernel`]).
    pub forces: Vec<f64>,
    /// Reference density ρ₀ for the Tait equation of state.
    pub rho0: f64,
    /// Speed of sound c₀ for the Tait EOS.
    pub c0: f64,
}

impl GpuSphGrid {
    /// Create a new `GpuSphGrid` with `n` particles.
    ///
    /// All physical quantities are initialised to zero; `rho0` defaults to
    /// `1000.0` kg/m³ (water) and `c0` to `1500.0` m/s.
    pub fn new(n: usize) -> Self {
        Self {
            positions: vec![0.0; n * 3],
            masses: vec![1.0; n],
            smoothing_lengths: vec![1.0; n],
            densities: vec![0.0; n],
            pressures: vec![0.0; n],
            velocities: vec![0.0; n * 3],
            forces: vec![0.0; n * 3],
            rho0: 1000.0,
            c0: 1500.0,
        }
    }

    /// Number of particles stored in this grid.
    pub fn particle_count(&self) -> usize {
        self.masses.len()
    }
}

// ── SPH kernel functions ─────────────────────────────────────────────────────

/// Cubic-spline SPH kernel W(r, h).
///
/// Returns the kernel value at distance `r` with smoothing length `h`.
#[allow(dead_code)]
pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient of the cubic-spline kernel ∇W(r, h) along the displacement vector.
///
/// Returns `[dW/dx, dW/dy, dW/dz]`.
#[allow(dead_code)]
pub fn cubic_spline_kernel_grad(dx: f64, dy: f64, dz: f64, h: f64) -> [f64; 3] {
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    if h <= 0.0 || r < 1e-15 {
        return [0.0; 3];
    }
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    let dw_dr = if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        -sigma * 0.75 * t * t / h
    } else {
        0.0
    };
    let inv_r = 1.0 / r;
    [dw_dr * dx * inv_r, dw_dr * dy * inv_r, dw_dr * dz * inv_r]
}

// ── GPU kernel mocks ─────────────────────────────────────────────────────────

/// Parallel density summation over all particle pairs (mock via loop).
///
/// For each particle `i`, computes:
/// `ρᵢ = Σⱼ mⱼ · W(|rᵢ − rⱼ|, hᵢ)`
///
/// Results are written into `grid.densities`.
pub fn gpu_density_kernel(grid: &mut GpuSphGrid) {
    let n = grid.particle_count();
    for i in 0..n {
        let mut rho = 0.0f64;
        let xi = grid.positions[i * 3];
        let yi = grid.positions[i * 3 + 1];
        let zi = grid.positions[i * 3 + 2];
        let hi = grid.smoothing_lengths[i];
        for j in 0..n {
            let xj = grid.positions[j * 3];
            let yj = grid.positions[j * 3 + 1];
            let zj = grid.positions[j * 3 + 2];
            let r = ((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2)).sqrt();
            rho += grid.masses[j] * cubic_spline_kernel(r, hi);
        }
        grid.densities[i] = rho;
    }
}

/// Vectorized Tait EOS pressure update.
///
/// For each particle `i`, computes:
/// `Pᵢ = (ρ₀·c₀²/γ) · [(ρᵢ/ρ₀)^γ − 1]`
///
/// Results are written into `grid.pressures`.
pub fn gpu_pressure_tait(grid: &mut GpuSphGrid) {
    let n = grid.particle_count();
    let prefactor = grid.rho0 * grid.c0 * grid.c0 / TAIT_GAMMA;
    for i in 0..n {
        let ratio = grid.densities[i] / grid.rho0;
        grid.pressures[i] = prefactor * (ratio.powf(TAIT_GAMMA) - 1.0);
    }
}

/// Pressure-gradient and viscosity force kernel (mock via loop).
///
/// Computes forces on each particle from pressure gradients and artificial
/// viscosity with coefficient `nu`. Results are accumulated into
/// `grid.forces`.
pub fn gpu_force_kernel(grid: &mut GpuSphGrid, nu: f64) {
    let n = grid.particle_count();
    // zero forces first
    for f in grid.forces.iter_mut() {
        *f = 0.0;
    }
    for i in 0..n {
        let xi = grid.positions[i * 3];
        let yi = grid.positions[i * 3 + 1];
        let zi = grid.positions[i * 3 + 2];
        let hi = grid.smoothing_lengths[i];
        let pi = grid.pressures[i];
        let rhoi = grid.densities[i];
        if rhoi < 1e-15 {
            continue;
        }
        let vxi = grid.velocities[i * 3];
        let vyi = grid.velocities[i * 3 + 1];
        let vzi = grid.velocities[i * 3 + 2];

        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = xi - grid.positions[j * 3];
            let dy = yi - grid.positions[j * 3 + 1];
            let dz = zi - grid.positions[j * 3 + 2];
            let rhoj = grid.densities[j];
            if rhoj < 1e-15 {
                continue;
            }
            let pj = grid.pressures[j];
            let mj = grid.masses[j];
            let grad = cubic_spline_kernel_grad(dx, dy, dz, hi);
            // pressure gradient term
            let coeff = -mj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
            // viscosity term (dot(v_ij, r_ij))
            let dvx = vxi - grid.velocities[j * 3];
            let dvy = vyi - grid.velocities[j * 3 + 1];
            let dvz = vzi - grid.velocities[j * 3 + 2];
            let r2 = dx * dx + dy * dy + dz * dz + 1e-15;
            let vdotr = dvx * dx + dvy * dy + dvz * dz;
            let visc = -nu * mj * vdotr / (rhoj * r2);

            grid.forces[i * 3] += (coeff + visc) * grad[0];
            grid.forces[i * 3 + 1] += (coeff + visc) * grad[1];
            grid.forces[i * 3 + 2] += (coeff + visc) * grad[2];
        }
    }
}

/// Build a cell-list neighbor search structure.
///
/// Divides the simulation domain `[min_x, max_x]³` into cubic cells of
/// size `cell_size`. Returns, for each particle `i`, the list of neighbour
/// indices within `2 * cell_size`.
///
/// Returns a `Vec<Vec`usize`>` where entry `i` holds the neighbours of
/// particle `i`.
pub fn gpu_neighbor_list(grid: &GpuSphGrid, cell_size: f64) -> Vec<Vec<usize>> {
    let n = grid.particle_count();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    let cutoff2 = (2.0 * cell_size) * (2.0 * cell_size);
    for i in 0..n {
        let xi = grid.positions[i * 3];
        let yi = grid.positions[i * 3 + 1];
        let zi = grid.positions[i * 3 + 2];
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = xi - grid.positions[j * 3];
            let dy = yi - grid.positions[j * 3 + 1];
            let dz = zi - grid.positions[j * 3 + 2];
            if dx * dx + dy * dy + dz * dz <= cutoff2 {
                neighbors[i].push(j);
            }
        }
    }
    neighbors
}

/// Orchestrate a full SPH density pass.
///
/// Runs, in order:
/// 1. [`gpu_density_kernel`] — compute densities
/// 2. [`gpu_pressure_tait`] — update pressures via Tait EOS
/// 3. [`gpu_force_kernel`]   — compute pressure + viscosity forces
pub fn launch_density_pass(grid: &mut GpuSphGrid, nu: f64) {
    gpu_density_kernel(grid);
    gpu_pressure_tait(grid);
    gpu_force_kernel(grid, nu);
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_two_particle_grid() -> GpuSphGrid {
        let mut g = GpuSphGrid::new(2);
        // particle 0 at origin, particle 1 at (0.5, 0, 0)
        g.positions = vec![0.0, 0.0, 0.0, 0.5, 0.0, 0.0];
        g.masses = vec![1.0, 1.0];
        g.smoothing_lengths = vec![1.0, 1.0];
        g.rho0 = 1000.0;
        g.c0 = 100.0;
        g
    }

    #[test]
    fn test_new_grid_particle_count() {
        let g = GpuSphGrid::new(10);
        assert_eq!(g.particle_count(), 10);
    }

    #[test]
    fn test_new_grid_default_rho0() {
        let g = GpuSphGrid::new(1);
        assert!((g.rho0 - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_new_grid_default_c0() {
        let g = GpuSphGrid::new(1);
        assert!((g.c0 - 1500.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_kernel_zero_distance() {
        // W(0, h) should be positive (self-contribution)
        let w = cubic_spline_kernel(0.0, 1.0);
        assert!(w > 0.0, "kernel at r=0 should be positive, got {w}");
    }

    #[test]
    fn test_cubic_kernel_beyond_support() {
        let w = cubic_spline_kernel(3.0, 1.0);
        assert!(
            (w).abs() < 1e-15,
            "kernel beyond 2h should be zero, got {w}"
        );
    }

    #[test]
    fn test_cubic_kernel_positive_within_support() {
        let w = cubic_spline_kernel(1.5, 1.0);
        assert!(w >= 0.0);
    }

    #[test]
    fn test_cubic_kernel_zero_h() {
        let w = cubic_spline_kernel(1.0, 0.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_cubic_kernel_grad_zero_displacement() {
        let g = cubic_spline_kernel_grad(0.0, 0.0, 0.0, 1.0);
        assert_eq!(g, [0.0; 3]);
    }

    #[test]
    fn test_cubic_kernel_grad_symmetry() {
        // Grad should be antisymmetric: ∇W(r) = -∇W(-r)
        let g1 = cubic_spline_kernel_grad(0.3, 0.0, 0.0, 1.0);
        let g2 = cubic_spline_kernel_grad(-0.3, 0.0, 0.0, 1.0);
        assert!((g1[0] + g2[0]).abs() < 1e-12);
    }

    #[test]
    fn test_density_kernel_single_particle() {
        let mut g = GpuSphGrid::new(1);
        g.positions = vec![0.0, 0.0, 0.0];
        g.masses = vec![1.0];
        g.smoothing_lengths = vec![1.0];
        gpu_density_kernel(&mut g);
        // Density = m * W(0, h) > 0
        assert!(g.densities[0] > 0.0);
    }

    #[test]
    fn test_density_kernel_two_particles() {
        let mut g = make_two_particle_grid();
        gpu_density_kernel(&mut g);
        assert!(g.densities[0] > 0.0);
        assert!(g.densities[1] > 0.0);
    }

    #[test]
    fn test_density_kernel_symmetric() {
        let mut g = make_two_particle_grid();
        gpu_density_kernel(&mut g);
        // Both particles have the same mass & smoothing length → equal density
        assert!((g.densities[0] - g.densities[1]).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_tait_zero_density() {
        let mut g = GpuSphGrid::new(1);
        g.densities = vec![0.0];
        g.rho0 = 1000.0;
        g.c0 = 100.0;
        gpu_pressure_tait(&mut g);
        // (0/rho0)^gamma - 1 < 0 → negative pressure
        assert!(g.pressures[0] < 0.0);
    }

    #[test]
    fn test_pressure_tait_at_reference_density() {
        let mut g = GpuSphGrid::new(1);
        g.rho0 = 1000.0;
        g.c0 = 100.0;
        g.densities = vec![g.rho0];
        gpu_pressure_tait(&mut g);
        // (rho0/rho0)^gamma - 1 = 0 → pressure = 0
        assert!(g.pressures[0].abs() < 1e-6);
    }

    #[test]
    fn test_pressure_tait_above_reference() {
        let mut g = GpuSphGrid::new(1);
        g.rho0 = 1000.0;
        g.c0 = 100.0;
        g.densities = vec![1100.0];
        gpu_pressure_tait(&mut g);
        assert!(g.pressures[0] > 0.0);
    }

    #[test]
    fn test_force_kernel_self_zero() {
        // With a single particle, pressure forces should be zero
        let mut g = GpuSphGrid::new(1);
        g.positions = vec![0.0, 0.0, 0.0];
        g.masses = vec![1.0];
        g.densities = vec![1000.0];
        g.pressures = vec![0.0];
        gpu_force_kernel(&mut g, 0.0);
        assert!((g.forces[0]).abs() < 1e-15);
    }

    #[test]
    fn test_force_kernel_repulsion() {
        // Two close particles with pressure > 0 should repel
        let mut g = make_two_particle_grid();
        gpu_density_kernel(&mut g);
        gpu_pressure_tait(&mut g);
        gpu_force_kernel(&mut g, 0.0);
        // Particle 0 at origin should be pushed in the -x direction (towards negative x)
        // because particle 1 is at +x
        let fx0 = g.forces[0];
        let fx1 = g.forces[3];
        // Forces should be opposite
        assert!(fx0 * fx1 < 0.0 || (fx0.abs() < 1e-12 && fx1.abs() < 1e-12));
    }

    #[test]
    fn test_force_kernel_zeros_forces_first() {
        let mut g = GpuSphGrid::new(2);
        g.forces = vec![9.9, 9.9, 9.9, 9.9, 9.9, 9.9];
        gpu_force_kernel(&mut g, 0.0);
        // With zero density, forces remain zero after re-initialisation
        for &f in &g.forces {
            assert!(f.abs() < 1e-15);
        }
    }

    #[test]
    fn test_neighbor_list_all_close() {
        let g = make_two_particle_grid();
        let nl = gpu_neighbor_list(&g, 1.0);
        // With cell_size=1.0, cutoff=2.0 → both particles are neighbours
        assert!(nl[0].contains(&1));
        assert!(nl[1].contains(&0));
    }

    #[test]
    fn test_neighbor_list_too_far() {
        let mut g = GpuSphGrid::new(2);
        g.positions = vec![0.0, 0.0, 0.0, 100.0, 0.0, 0.0];
        let nl = gpu_neighbor_list(&g, 1.0);
        assert!(nl[0].is_empty());
        assert!(nl[1].is_empty());
    }

    #[test]
    fn test_neighbor_list_no_self() {
        let g = make_two_particle_grid();
        let nl = gpu_neighbor_list(&g, 1.0);
        assert!(!nl[0].contains(&0));
        assert!(!nl[1].contains(&1));
    }

    #[test]
    fn test_launch_density_pass_updates_all() {
        let mut g = make_two_particle_grid();
        launch_density_pass(&mut g, 0.01);
        // After a full pass everything should be non-zero (apart from forces
        // which may still be small)
        assert!(g.densities[0] > 0.0);
        assert!(g.densities[1] > 0.0);
    }

    #[test]
    fn test_launch_density_pass_idempotent() {
        let mut g1 = make_two_particle_grid();
        let mut g2 = make_two_particle_grid();
        launch_density_pass(&mut g1, 0.0);
        launch_density_pass(&mut g2, 0.0);
        for i in 0..2 {
            assert!((g1.densities[i] - g2.densities[i]).abs() < 1e-12);
            assert!((g1.pressures[i] - g2.pressures[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_force_magnitude_finite() {
        let mut g = make_two_particle_grid();
        launch_density_pass(&mut g, 0.01);
        for &f in &g.forces {
            assert!(f.is_finite(), "force is not finite: {f}");
        }
    }

    #[test]
    fn test_density_five_particles() {
        let n = 5;
        let mut g = GpuSphGrid::new(n);
        // Place particles in a line
        for i in 0..n {
            g.positions[i * 3] = i as f64 * 0.4;
        }
        gpu_density_kernel(&mut g);
        // Central particles should have higher density
        assert!(g.densities[2] >= g.densities[0]);
    }

    #[test]
    fn test_pressure_increases_with_density() {
        let mut g = GpuSphGrid::new(2);
        g.rho0 = 1000.0;
        g.c0 = 100.0;
        g.densities = vec![1000.0, 1200.0];
        gpu_pressure_tait(&mut g);
        assert!(g.pressures[1] > g.pressures[0]);
    }

    #[test]
    fn test_gpu_sph_grid_clone() {
        let g = GpuSphGrid::new(3);
        let g2 = g.clone();
        assert_eq!(g2.particle_count(), 3);
    }

    #[test]
    fn test_gpu_sph_grid_debug() {
        let g = GpuSphGrid::new(1);
        let s = format!("{g:?}");
        assert!(s.contains("GpuSphGrid"));
    }
}
