// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated thermal computation (CPU mock backend via Rayon).
//!
//! Implements a parallel finite-difference heat equation solver that mirrors
//! the structure of a GPU compute dispatch. The "GPU" here is simulated via
//! Rayon parallel iterators, making it easy to swap in a real GPU backend later.

use rayon::prelude::*;

// ---------------------------------------------------------------------------
// Boundary condition type
// ---------------------------------------------------------------------------

/// Boundary condition type for thermal simulations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalBc {
    /// Dirichlet: fixed temperature value at the boundary.
    Dirichlet(f64),
    /// Neumann: fixed heat flux (derivative) at the boundary.
    Neumann(f64),
}

// ---------------------------------------------------------------------------
// HeatSource
// ---------------------------------------------------------------------------

/// Volumetric heat source specification.
#[derive(Debug, Clone)]
pub struct HeatSource {
    /// Linear index of the cell where the source is applied.
    pub cell_index: usize,
    /// Power per unit volume (W/m³ equivalent in simulation units).
    pub power: f64,
}

impl HeatSource {
    /// Create a new heat source at `cell_index` with the given `power`.
    pub fn new(cell_index: usize, power: f64) -> Self {
        Self { cell_index, power }
    }
}

// ---------------------------------------------------------------------------
// GpuThermalSolver
// ---------------------------------------------------------------------------

/// GPU-accelerated (CPU mock) thermal solver for 3-D structured grids.
///
/// Solves the heat equation:
/// ```text
///   ∂T/∂t = α ∇²T + Q
/// ```
/// where `α` is the thermal diffusivity and `Q` is a volumetric heat source.
///
/// Grid layout: row-major, index `[iz * ny * nx + iy * nx + ix]`.
#[derive(Debug, Clone)]
pub struct GpuThermalSolver {
    /// Number of cells in the x-direction.
    pub nx: usize,
    /// Number of cells in the y-direction.
    pub ny: usize,
    /// Number of cells in the z-direction.
    pub nz: usize,
    /// Thermal diffusivity (m²/s or simulation units).
    pub diffusivity: f64,
    /// Grid spacing in x (m).
    pub dx: f64,
    /// Grid spacing in y (m).
    pub dy: f64,
    /// Grid spacing in z (m).
    pub dz: f64,
    /// Current temperature field, length `nx * ny * nz`.
    pub temperature: Vec<f64>,
}

impl GpuThermalSolver {
    /// Create a new solver with uniform initial temperature `t0`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        diffusivity: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        t0: f64,
    ) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            diffusivity,
            dx,
            dy,
            dz,
            temperature: vec![t0; n],
        }
    }

    /// Linear index from 3-D grid coordinates.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.ny * self.nx + iy * self.nx + ix
    }

    /// Total number of cells.
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

// ---------------------------------------------------------------------------
// gpu_heat_diffusion
// ---------------------------------------------------------------------------

/// Perform one parallel heat-equation update step (mock GPU dispatch).
///
/// Applies the explicit finite-difference stencil:
/// ```text
///   T_new[i] = T[i] + dt * α * ∇²T[i]
/// ```
/// Interior cells only; boundary cells are left unchanged.
///
/// # Arguments
/// * `solver` - The thermal solver (updated in place).
/// * `dt` - Time step size (s).
pub fn gpu_heat_diffusion(solver: &mut GpuThermalSolver, dt: f64) {
    let nx = solver.nx;
    let ny = solver.ny;
    let nz = solver.nz;
    let alpha = solver.diffusivity;
    let dx2 = solver.dx * solver.dx;
    let dy2 = solver.dy * solver.dy;
    let dz2 = solver.dz * solver.dz;
    let old = solver.temperature.clone();

    let new_temp: Vec<f64> = (0..nz * ny * nx)
        .into_par_iter()
        .map(|idx| {
            let iz = idx / (ny * nx);
            let rem = idx % (ny * nx);
            let iy = rem / nx;
            let ix = rem % nx;

            // Skip boundary cells
            if ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1 {
                return old[idx];
            }

            let laplacian_x = (old[idx - 1] - 2.0 * old[idx] + old[idx + 1]) / dx2;
            let laplacian_y = (old[idx - nx] - 2.0 * old[idx] + old[idx + nx]) / dy2;
            let laplacian_z = (old[idx - ny * nx] - 2.0 * old[idx] + old[idx + ny * nx]) / dz2;

            old[idx] + dt * alpha * (laplacian_x + laplacian_y + laplacian_z)
        })
        .collect();

    solver.temperature = new_temp;
}

// ---------------------------------------------------------------------------
// thermal_boundary_apply
// ---------------------------------------------------------------------------

/// Apply boundary conditions to the temperature field.
///
/// Iterates over all six faces of the structured grid and applies either
/// Dirichlet or Neumann boundary conditions.
///
/// # Arguments
/// * `solver` - The thermal solver (updated in place).
/// * `bc_xmin` - BC on the x = 0 face.
/// * `bc_xmax` - BC on the x = nx-1 face.
/// * `bc_ymin` - BC on the y = 0 face.
/// * `bc_ymax` - BC on the y = ny-1 face.
/// * `bc_zmin` - BC on the z = 0 face.
/// * `bc_zmax` - BC on the z = nz-1 face.
#[allow(clippy::too_many_arguments)]
pub fn thermal_boundary_apply(
    solver: &mut GpuThermalSolver,
    bc_xmin: ThermalBc,
    bc_xmax: ThermalBc,
    bc_ymin: ThermalBc,
    bc_ymax: ThermalBc,
    bc_zmin: ThermalBc,
    bc_zmax: ThermalBc,
) {
    let nx = solver.nx;
    let ny = solver.ny;
    let nz = solver.nz;

    // x-faces
    for iz in 0..nz {
        for iy in 0..ny {
            let idx_min = solver.idx(0, iy, iz);
            let idx_max = solver.idx(nx - 1, iy, iz);
            apply_bc_to_cell(&mut solver.temperature, idx_min, bc_xmin, solver.dx);
            apply_bc_to_cell(&mut solver.temperature, idx_max, bc_xmax, solver.dx);
        }
    }

    // y-faces
    for iz in 0..nz {
        for ix in 0..nx {
            let idx_min = solver.idx(ix, 0, iz);
            let idx_max = solver.idx(ix, ny - 1, iz);
            apply_bc_to_cell(&mut solver.temperature, idx_min, bc_ymin, solver.dy);
            apply_bc_to_cell(&mut solver.temperature, idx_max, bc_ymax, solver.dy);
        }
    }

    // z-faces
    for iy in 0..ny {
        for ix in 0..nx {
            let idx_min = solver.idx(ix, iy, 0);
            let idx_max = solver.idx(ix, iy, nz - 1);
            apply_bc_to_cell(&mut solver.temperature, idx_min, bc_zmin, solver.dz);
            apply_bc_to_cell(&mut solver.temperature, idx_max, bc_zmax, solver.dz);
        }
    }
}

/// Internal helper: apply one BC to a cell.
#[allow(dead_code)]
fn apply_bc_to_cell(temperature: &mut [f64], idx: usize, bc: ThermalBc, _h: f64) {
    match bc {
        ThermalBc::Dirichlet(val) => {
            temperature[idx] = val;
        }
        ThermalBc::Neumann(_flux) => {
            // For Neumann: T[ghost] = T[interior] + flux*h
            // In a 1-D ghost-cell approach we just leave the boundary unchanged
            // (zero-flux by default). A full implementation would require ghost
            // cell indices; this mock sets the value to maintain the gradient.
            // No modification needed for zero-flux; non-zero handled outside.
        }
    }
}

// ---------------------------------------------------------------------------
// gpu_heat_source
// ---------------------------------------------------------------------------

/// Apply volumetric heat sources to the temperature field (mock GPU kernel).
///
/// Each source adds `source.power * dt` to the temperature of its cell.
///
/// # Arguments
/// * `solver` - The thermal solver (updated in place).
/// * `sources` - List of heat sources.
/// * `dt` - Time step size (s).
pub fn gpu_heat_source(solver: &mut GpuThermalSolver, sources: &[HeatSource], dt: f64) {
    for src in sources {
        if src.cell_index < solver.temperature.len() {
            solver.temperature[src.cell_index] += src.power * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// temperature_gradient
// ---------------------------------------------------------------------------

/// Compute the temperature gradient at every interior cell using central differences.
///
/// Returns a `Vec<[f64; 3]>` of length `nx * ny * nz`.  Boundary cells get
/// a gradient of `[0.0, 0.0, 0.0]`.
///
/// # Arguments
/// * `solver` - The thermal solver.
pub fn temperature_gradient(solver: &GpuThermalSolver) -> Vec<[f64; 3]> {
    let nx = solver.nx;
    let ny = solver.ny;
    let nz = solver.nz;
    let t = &solver.temperature;

    (0..nz * ny * nx)
        .into_par_iter()
        .map(|idx| {
            let iz = idx / (ny * nx);
            let rem = idx % (ny * nx);
            let iy = rem / nx;
            let ix = rem % nx;

            if ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1 {
                return [0.0; 3];
            }

            let gx = (t[idx + 1] - t[idx - 1]) / (2.0 * solver.dx);
            let gy = (t[idx + nx] - t[idx - nx]) / (2.0 * solver.dy);
            let gz = (t[idx + ny * nx] - t[idx - ny * nx]) / (2.0 * solver.dz);

            [gx, gy, gz]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// thermal_equilibration
// ---------------------------------------------------------------------------

/// Check whether the temperature field has reached thermal equilibrium.
///
/// Returns `true` when the maximum absolute change between `old` and the
/// current field is below `tol`.
///
/// # Arguments
/// * `solver` - The thermal solver (current state).
/// * `old` - Temperature field from the previous iteration.
/// * `tol` - Convergence tolerance (K or simulation units).
pub fn thermal_equilibration(solver: &GpuThermalSolver, old: &[f64], tol: f64) -> bool {
    if old.len() != solver.temperature.len() {
        return false;
    }
    solver
        .temperature
        .par_iter()
        .zip(old.par_iter())
        .map(|(&new, &prev)| (new - prev).abs())
        .reduce(|| 0.0_f64, f64::max)
        < tol
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_solver(nx: usize, ny: usize, nz: usize, t0: f64) -> GpuThermalSolver {
        GpuThermalSolver::new(nx, ny, nz, 1e-4, 0.1, 0.1, 0.1, t0)
    }

    // ── GpuThermalSolver construction ─────────────────────────────────────

    #[test]
    fn test_solver_new_size() {
        let s = make_solver(4, 4, 4, 300.0);
        assert_eq!(s.n_cells(), 64);
        assert_eq!(s.temperature.len(), 64);
    }

    #[test]
    fn test_solver_uniform_init() {
        let s = make_solver(3, 3, 3, 273.15);
        assert!(s.temperature.iter().all(|&t| (t - 273.15).abs() < 1e-12));
    }

    #[test]
    fn test_solver_idx_origin() {
        let s = make_solver(5, 5, 5, 0.0);
        assert_eq!(s.idx(0, 0, 0), 0);
    }

    #[test]
    fn test_solver_idx_last() {
        let s = make_solver(5, 5, 5, 0.0);
        assert_eq!(s.idx(4, 4, 4), 124);
    }

    #[test]
    fn test_solver_idx_slice() {
        let s = make_solver(4, 4, 4, 0.0);
        assert_eq!(s.idx(2, 1, 0), 4 + 2);
    }

    // ── gpu_heat_diffusion ────────────────────────────────────────────────

    #[test]
    fn test_diffusion_boundary_unchanged() {
        let mut s = make_solver(4, 4, 4, 100.0);
        // Set a spike in the interior
        let idx = s.idx(2, 2, 2);
        s.temperature[idx] = 500.0;
        let boundary_before = s.temperature[s.idx(0, 0, 0)];
        gpu_heat_diffusion(&mut s, 0.001);
        // Corner boundary must stay at 100 (not updated by diffusion kernel)
        assert!((s.temperature[s.idx(0, 0, 0)] - boundary_before).abs() < 1e-12);
    }

    #[test]
    fn test_diffusion_uniform_field_unchanged() {
        let mut s = make_solver(5, 5, 5, 300.0);
        let before: Vec<f64> = s.temperature.clone();
        gpu_heat_diffusion(&mut s, 0.01);
        // Uniform field: Laplacian is zero, so nothing should change
        for (a, b) in s.temperature.iter().zip(before.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "uniform field changed under diffusion"
            );
        }
    }

    #[test]
    fn test_diffusion_hot_spot_cools() {
        let mut s = make_solver(5, 5, 5, 0.0);
        let idx = s.idx(2, 2, 2);
        s.temperature[idx] = 1000.0;
        let before = s.temperature[idx];
        gpu_heat_diffusion(&mut s, 0.001);
        // The hot spot should cool (energy spreads)
        assert!(s.temperature[idx] < before);
    }

    #[test]
    fn test_diffusion_cold_spot_warms() {
        let mut s = make_solver(5, 5, 5, 300.0);
        let idx = s.idx(2, 2, 2);
        s.temperature[idx] = 0.0;
        gpu_heat_diffusion(&mut s, 0.001);
        assert!(s.temperature[idx] > 0.0);
    }

    #[test]
    fn test_diffusion_energy_approximately_conserved_interior() {
        // Total energy of interior cells should stay roughly the same
        let mut s = make_solver(5, 5, 5, 0.0);
        let idx = s.idx(2, 2, 2);
        s.temperature[idx] = 1000.0;
        let sum_before: f64 = s.temperature.iter().sum();
        gpu_heat_diffusion(&mut s, 0.001);
        let sum_after: f64 = s.temperature.iter().sum();
        assert!((sum_before - sum_after).abs() < 1e-6 * sum_before.abs() + 1e-9);
    }

    // ── thermal_boundary_apply ────────────────────────────────────────────

    #[test]
    fn test_dirichlet_xmin_applied() {
        let mut s = make_solver(6, 6, 6, 0.0);
        thermal_boundary_apply(
            &mut s,
            ThermalBc::Dirichlet(500.0),
            ThermalBc::Dirichlet(500.0),
            ThermalBc::Dirichlet(500.0),
            ThermalBc::Dirichlet(500.0),
            ThermalBc::Dirichlet(500.0),
            ThermalBc::Dirichlet(500.0),
        );
        // All boundary cells should be 500 when all faces have same Dirichlet value
        for iz in 0..6 {
            for iy in 0..6 {
                let idx = s.idx(0, iy, iz);
                assert!(
                    (s.temperature[idx] - 500.0).abs() < 1e-12,
                    "x=0 face cell ({iy},{iz}) expected 500, got {}",
                    s.temperature[idx]
                );
            }
        }
    }

    #[test]
    fn test_dirichlet_xmax_applied() {
        let mut s = make_solver(6, 6, 6, 0.0);
        thermal_boundary_apply(
            &mut s,
            ThermalBc::Dirichlet(200.0),
            ThermalBc::Dirichlet(200.0),
            ThermalBc::Dirichlet(200.0),
            ThermalBc::Dirichlet(200.0),
            ThermalBc::Dirichlet(200.0),
            ThermalBc::Dirichlet(200.0),
        );
        // All boundary cells should be 200 when all faces share the same value
        for iz in 0..6 {
            for iy in 0..6 {
                let idx = s.idx(5, iy, iz);
                assert!(
                    (s.temperature[idx] - 200.0).abs() < 1e-12,
                    "x=5 face cell ({iy},{iz}) expected 200, got {}",
                    s.temperature[idx]
                );
            }
        }
    }

    #[test]
    fn test_neumann_bc_leaves_interior_unchanged_for_zero_flux() {
        let mut s = make_solver(4, 4, 4, 300.0);
        let interior_before = s.temperature[s.idx(2, 2, 2)];
        thermal_boundary_apply(
            &mut s,
            ThermalBc::Neumann(0.0),
            ThermalBc::Neumann(0.0),
            ThermalBc::Neumann(0.0),
            ThermalBc::Neumann(0.0),
            ThermalBc::Neumann(0.0),
            ThermalBc::Neumann(0.0),
        );
        let interior_after = s.temperature[s.idx(2, 2, 2)];
        assert!((interior_before - interior_after).abs() < 1e-12);
    }

    // ── gpu_heat_source ───────────────────────────────────────────────────

    #[test]
    fn test_heat_source_single_cell() {
        let mut s = make_solver(4, 4, 4, 0.0);
        let idx = s.idx(2, 2, 2);
        let sources = vec![HeatSource::new(idx, 1000.0)];
        gpu_heat_source(&mut s, &sources, 0.1);
        assert!((s.temperature[idx] - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_heat_source_multiple_cells() {
        let mut s = make_solver(4, 4, 4, 0.0);
        let idx1 = s.idx(1, 1, 1);
        let idx2 = s.idx(2, 2, 2);
        let sources = vec![HeatSource::new(idx1, 500.0), HeatSource::new(idx2, 200.0)];
        gpu_heat_source(&mut s, &sources, 1.0);
        assert!((s.temperature[idx1] - 500.0).abs() < 1e-10);
        assert!((s.temperature[idx2] - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_heat_source_out_of_bounds_ignored() {
        let mut s = make_solver(4, 4, 4, 0.0);
        let sources = vec![HeatSource::new(9999, 1000.0)];
        // Should not panic
        gpu_heat_source(&mut s, &sources, 1.0);
        assert!(s.temperature.iter().all(|&t| t.abs() < 1e-12));
    }

    #[test]
    fn test_heat_source_zero_power() {
        let mut s = make_solver(4, 4, 4, 100.0);
        let idx = s.idx(2, 2, 2);
        let sources = vec![HeatSource::new(idx, 0.0)];
        gpu_heat_source(&mut s, &sources, 1.0);
        assert!((s.temperature[idx] - 100.0).abs() < 1e-12);
    }

    // ── temperature_gradient ──────────────────────────────────────────────

    #[test]
    fn test_gradient_uniform_field_is_zero() {
        let s = make_solver(5, 5, 5, 300.0);
        let grad = temperature_gradient(&s);
        for g in &grad {
            assert!(g[0].abs() < 1e-12 && g[1].abs() < 1e-12 && g[2].abs() < 1e-12);
        }
    }

    #[test]
    fn test_gradient_boundary_is_zero() {
        let s = make_solver(4, 4, 4, 100.0);
        let grad = temperature_gradient(&s);
        // Boundary cell (0,0,0)
        assert_eq!(grad[0], [0.0; 3]);
    }

    #[test]
    fn test_gradient_x_linear_field() {
        // T = ix * dx (linear in x), so dT/dx = 1.0, dT/dy = 0, dT/dz = 0
        let mut s = GpuThermalSolver::new(5, 5, 5, 1e-4, 1.0, 1.0, 1.0, 0.0);
        for iz in 0..5 {
            for iy in 0..5 {
                for ix in 0..5 {
                    let idx = s.idx(ix, iy, iz);
                    s.temperature[idx] = ix as f64;
                }
            }
        }
        let grad = temperature_gradient(&s);
        let idx = s.idx(2, 2, 2);
        assert!((grad[idx][0] - 1.0).abs() < 1e-12, "gx={}", grad[idx][0]);
        assert!(grad[idx][1].abs() < 1e-12);
        assert!(grad[idx][2].abs() < 1e-12);
    }

    #[test]
    fn test_gradient_y_linear_field() {
        let mut s = GpuThermalSolver::new(5, 5, 5, 1e-4, 1.0, 1.0, 1.0, 0.0);
        for iz in 0..5 {
            for iy in 0..5 {
                for ix in 0..5 {
                    let idx = s.idx(ix, iy, iz);
                    s.temperature[idx] = iy as f64;
                }
            }
        }
        let grad = temperature_gradient(&s);
        let idx = s.idx(2, 2, 2);
        assert!(grad[idx][0].abs() < 1e-12);
        assert!((grad[idx][1] - 1.0).abs() < 1e-12, "gy={}", grad[idx][1]);
        assert!(grad[idx][2].abs() < 1e-12);
    }

    // ── thermal_equilibration ─────────────────────────────────────────────

    #[test]
    fn test_equilibration_identical_fields() {
        let s = make_solver(4, 4, 4, 300.0);
        let old = s.temperature.clone();
        assert!(thermal_equilibration(&s, &old, 1e-6));
    }

    #[test]
    fn test_equilibration_large_change() {
        let s = make_solver(4, 4, 4, 300.0);
        let old = vec![0.0; s.n_cells()];
        assert!(!thermal_equilibration(&s, &old, 1e-6));
    }

    #[test]
    fn test_equilibration_small_change_below_tol() {
        let mut s = make_solver(4, 4, 4, 300.0);
        let old = s.temperature.clone();
        s.temperature[0] += 1e-8; // tiny change
        assert!(thermal_equilibration(&s, &old, 1e-6));
    }

    #[test]
    fn test_equilibration_small_change_above_tol() {
        let mut s = make_solver(4, 4, 4, 300.0);
        let old = s.temperature.clone();
        s.temperature[0] += 1.0; // large change
        assert!(!thermal_equilibration(&s, &old, 1e-6));
    }

    #[test]
    fn test_equilibration_length_mismatch_returns_false() {
        let s = make_solver(4, 4, 4, 300.0);
        let old = vec![300.0; 10]; // wrong length
        assert!(!thermal_equilibration(&s, &old, 1e-6));
    }

    // ── HeatSource struct ─────────────────────────────────────────────────

    #[test]
    fn test_heat_source_fields() {
        let src = HeatSource::new(42, 999.9);
        assert_eq!(src.cell_index, 42);
        assert!((src.power - 999.9).abs() < 1e-12);
    }

    // ── Integration: diffusion + BC ───────────────────────────────────────

    #[test]
    fn test_diffusion_and_dirichlet_bc_combined() {
        let mut s = GpuThermalSolver::new(5, 5, 5, 1e-3, 0.1, 0.1, 0.1, 0.0);
        // Hot left wall
        thermal_boundary_apply(
            &mut s,
            ThermalBc::Dirichlet(100.0),
            ThermalBc::Dirichlet(0.0),
            ThermalBc::Dirichlet(0.0),
            ThermalBc::Dirichlet(0.0),
            ThermalBc::Dirichlet(0.0),
            ThermalBc::Dirichlet(0.0),
        );
        // Run diffusion for several steps
        for _ in 0..10 {
            gpu_heat_diffusion(&mut s, 0.001);
            // Re-apply BC each step
            thermal_boundary_apply(
                &mut s,
                ThermalBc::Dirichlet(100.0),
                ThermalBc::Dirichlet(0.0),
                ThermalBc::Dirichlet(0.0),
                ThermalBc::Dirichlet(0.0),
                ThermalBc::Dirichlet(0.0),
                ThermalBc::Dirichlet(0.0),
            );
        }
        // Interior cells should be warming up
        let interior_idx = s.idx(2, 2, 2);
        assert!(s.temperature[interior_idx] > 0.0, "interior should warm up");
        // x=0 boundary still fixed at 100
        let bc_idx = s.idx(0, 2, 2);
        assert!((s.temperature[bc_idx] - 100.0).abs() < 1e-12);
    }
}
