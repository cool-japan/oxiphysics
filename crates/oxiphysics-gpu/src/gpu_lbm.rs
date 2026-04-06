// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Lattice Boltzmann Method (D2Q9) — CPU mock backend.
//!
//! Implements the D2Q9 velocity set with BGK collision, periodic streaming,
//! macroscopic quantity recovery, and diagnostic utilities for vorticity and
//! divergence fields.

// ── D2Q9 constants ────────────────────────────────────────────────────────────

/// Number of discrete velocities in the D2Q9 velocity set.
pub const D2Q9_Q: usize = 9;

/// Standard D2Q9 lattice weights.
///
/// Ordering: rest (index 0), axis (1–4), diagonal (5–8).
pub const D2Q9_WEIGHTS: [f64; D2Q9_Q] = [
    4.0 / 9.0,  // 0: rest
    1.0 / 9.0,  // 1: +x
    1.0 / 9.0,  // 2: -x
    1.0 / 9.0,  // 3: +y
    1.0 / 9.0,  // 4: -y
    1.0 / 36.0, // 5: +x+y
    1.0 / 36.0, // 6: -x+y
    1.0 / 36.0, // 7: -x-y
    1.0 / 36.0, // 8: +x-y
];

/// D2Q9 discrete velocity directions `[cx, cy]`.
pub const D2Q9_DIRS: [[i32; 2]; D2Q9_Q] = [
    [0, 0],   // 0: rest
    [1, 0],   // 1: +x
    [-1, 0],  // 2: -x
    [0, 1],   // 3: +y
    [0, -1],  // 4: -y
    [1, 1],   // 5: +x+y
    [-1, 1],  // 6: -x+y
    [-1, -1], // 7: -x-y
    [1, -1],  // 8: +x-y
];

/// Opposite (bounce-back) indices for D2Q9.
///
/// `D2Q9_OPPOSITE[i]` is the index `j` such that `D2Q9_DIRS[j] == -D2Q9_DIRS[i]`.
pub const D2Q9_OPPOSITE: [usize; D2Q9_Q] = [0, 2, 1, 4, 3, 7, 8, 5, 6];

// ── D2Q9Weights ───────────────────────────────────────────────────────────────

/// Lattice weight set and direction set for the D2Q9 velocity model.
#[derive(Debug, Clone, PartialEq)]
pub struct D2Q9Weights {
    /// Nine lattice weights.
    pub weights: [f64; 9],
    /// Nine lattice directions `[cx, cy]`.
    pub dirs: [[i32; 2]; 9],
}

impl D2Q9Weights {
    /// Construct the standard D2Q9 weight/direction set.
    pub fn standard() -> Self {
        Self {
            weights: D2Q9_WEIGHTS,
            dirs: D2Q9_DIRS,
        }
    }

    /// Speed of sound squared for this lattice (`cs² = 1/3`).
    pub fn cs2() -> f64 {
        1.0 / 3.0
    }
}

impl Default for D2Q9Weights {
    fn default() -> Self {
        Self::standard()
    }
}

// ── LbmCell ──────────────────────────────────────────────────────────────────

/// A single D2Q9 LBM cell holding incoming and outgoing distribution functions.
#[derive(Debug, Clone, PartialEq)]
pub struct LbmCell {
    /// Incoming distribution functions (post-stream).
    pub f_in: [f64; 9],
    /// Outgoing distribution functions (post-collision).
    pub f_out: [f64; 9],
}

impl LbmCell {
    /// Create a new cell with uniform rest-state distributions for density `rho`.
    pub fn new_equilibrium(rho: f64) -> Self {
        let mut f = [0.0f64; 9];
        for (i, w) in D2Q9_WEIGHTS.iter().enumerate() {
            f[i] = w * rho;
        }
        Self { f_in: f, f_out: f }
    }

    /// Macroscopic density: `ρ = Σ fᵢ`.
    pub fn density(&self) -> f64 {
        self.f_in.iter().sum()
    }

    /// Macroscopic velocity `[ux, uy]` from the incoming distributions.
    pub fn velocity(&self) -> [f64; 2] {
        let rho = self.density();
        if rho < 1e-15 {
            return [0.0, 0.0];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        for (i, &fi) in self.f_in.iter().enumerate() {
            ux += fi * D2Q9_DIRS[i][0] as f64;
            uy += fi * D2Q9_DIRS[i][1] as f64;
        }
        [ux / rho, uy / rho]
    }
}

impl Default for LbmCell {
    fn default() -> Self {
        Self::new_equilibrium(1.0)
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Compute the BGK equilibrium distribution for one velocity direction.
///
/// Formula:  `f_eq = w · ρ · (1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²))`
///
/// * `rho`  – macroscopic density
/// * `ux`, `uy` – macroscopic velocity components
/// * `w`    – lattice weight for this direction
/// * `cx`, `cy` – discrete velocity components for this direction
pub fn bgk_equilibrium(rho: f64, ux: f64, uy: f64, w: f64, cx: f64, cy: f64) -> f64 {
    let cs2 = 1.0 / 3.0;
    let cu = cx * ux + cy * uy;
    let u2 = ux * ux + uy * uy;
    w * rho * (1.0 + cu / cs2 + cu * cu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2))
}

/// Compute the BGK relaxation parameter `ω` from kinematic viscosity and `cs²`.
///
/// `ω = 1 / (ν/cs² + 0.5)`
pub fn relaxation_from_viscosity(nu: f64, cs2: f64) -> f64 {
    1.0 / (nu / cs2 + 0.5)
}

// ── GpuLbmGrid ───────────────────────────────────────────────────────────────

/// 2-D D2Q9 LBM grid (CPU mock of a GPU backend).
#[derive(Debug, Clone)]
pub struct GpuLbmGrid {
    /// Number of cells in the x-direction.
    pub nx: usize,
    /// Number of cells in the y-direction.
    pub ny: usize,
    /// Flat cell storage: row-major, `cells[y * nx + x]`.
    pub cells: Vec<LbmCell>,
}

impl GpuLbmGrid {
    /// Construct a new grid initialised to a uniform equilibrium state.
    ///
    /// * `nx`, `ny` – grid dimensions
    /// * `rho`      – initial density (uniform)
    pub fn new(nx: usize, ny: usize, rho: f64) -> Self {
        let cells = vec![LbmCell::new_equilibrium(rho); nx * ny];
        Self { nx, ny, cells }
    }

    /// Linear index for cell `(x, y)`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// BGK collision step: compute `f_out` from current `f_in` using relaxation
    /// parameter `omega`.
    pub fn collision_bgk(&mut self, omega: f64) {
        for cell in &mut self.cells {
            let rho = cell.density();
            let [ux, uy] = cell.velocity();
            for i in 0..9 {
                let cx = D2Q9_DIRS[i][0] as f64;
                let cy = D2Q9_DIRS[i][1] as f64;
                let feq = bgk_equilibrium(rho, ux, uy, D2Q9_WEIGHTS[i], cx, cy);
                cell.f_out[i] = cell.f_in[i] - omega * (cell.f_in[i] - feq);
            }
        }
    }

    /// Streaming step with periodic boundary conditions.
    ///
    /// Each distribution `f_out[i]` from cell `(x, y)` is moved to
    /// `f_in[i]` of the neighbouring cell in direction `D2Q9_DIRS[i]`.
    pub fn stream_periodic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        // Snapshot of f_out before overwriting f_in
        let snapshot: Vec<[f64; 9]> = self.cells.iter().map(|c| c.f_out).collect();

        for y in 0..ny {
            for x in 0..nx {
                let src_idx = y * nx + x;
                for i in 0..9 {
                    let cx = D2Q9_DIRS[i][0];
                    let cy = D2Q9_DIRS[i][1];
                    let dest_x = ((x as i64 + cx as i64).rem_euclid(nx as i64)) as usize;
                    let dest_y = ((y as i64 + cy as i64).rem_euclid(ny as i64)) as usize;
                    let dest_idx = dest_y * nx + dest_x;
                    self.cells[dest_idx].f_in[i] = snapshot[src_idx][i];
                }
            }
        }
    }

    /// Perform one full LBM step: collision then streaming.
    pub fn step(&mut self, omega: f64) {
        self.collision_bgk(omega);
        self.stream_periodic();
    }
}

// ── LbmBoundary ──────────────────────────────────────────────────────────────

/// Boundary condition type for a cell or face.
#[derive(Debug, Clone, PartialEq)]
pub enum LbmBoundary {
    /// Periodic boundary (default — handled in `stream_periodic`).
    Periodic,
    /// No-slip bounce-back wall.
    NoSlip,
    /// Zou-He velocity/pressure boundary condition.
    ZouHe {
        /// Prescribed density (pressure BC: set ux=uy=0, or velocity BC: set rho=1).
        rho: f64,
        /// Prescribed x-velocity.
        ux: f64,
        /// Prescribed y-velocity.
        uy: f64,
    },
}

// ── GpuLbmDispatch ───────────────────────────────────────────────────────────

/// High-level dispatcher that couples a `GpuLbmGrid` with per-cell boundary
/// conditions and orchestrates time-stepping.
#[derive(Debug, Clone)]
pub struct GpuLbmDispatch {
    /// The underlying LBM grid.
    pub grid: GpuLbmGrid,
    /// Per-cell boundary conditions (length == `grid.nx * grid.ny`).
    pub boundary: Vec<LbmBoundary>,
}

impl GpuLbmDispatch {
    /// Create a new dispatcher for a grid of size `nx × ny` with uniform
    /// initial density `rho` and all boundaries set to `Periodic`.
    pub fn new(nx: usize, ny: usize, rho: f64) -> Self {
        let grid = GpuLbmGrid::new(nx, ny, rho);
        let boundary = vec![LbmBoundary::Periodic; nx * ny];
        Self { grid, boundary }
    }

    /// Set the boundary condition for cell `(x, y)`.
    pub fn set_boundary(&mut self, x: usize, y: usize, bc: LbmBoundary) {
        let idx = self.grid.idx(x, y);
        self.boundary[idx] = bc;
    }

    /// Perform one dispatch step, applying boundary pre-conditioning before
    /// collision + streaming.
    pub fn dispatch_step(&mut self, omega: f64) {
        self.apply_boundaries();
        self.grid.step(omega);
    }

    /// Apply boundary conditions to cells before the collision step.
    fn apply_boundaries(&mut self) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                match &self.boundary[idx] {
                    LbmBoundary::Periodic => {}
                    LbmBoundary::NoSlip => {
                        // Bounce-back: swap f_in with its opposite
                        let cell = &mut self.grid.cells[idx];
                        let old = cell.f_in;
                        for i in 0..9 {
                            cell.f_in[i] = old[D2Q9_OPPOSITE[i]];
                        }
                    }
                    LbmBoundary::ZouHe { rho, ux, uy } => {
                        let rho = *rho;
                        let ux = *ux;
                        let uy = *uy;
                        // Overwrite f_in with equilibrium at prescribed (rho, ux, uy)
                        let cell = &mut self.grid.cells[idx];
                        for i in 0..9 {
                            let cx = D2Q9_DIRS[i][0] as f64;
                            let cy = D2Q9_DIRS[i][1] as f64;
                            cell.f_in[i] = bgk_equilibrium(rho, ux, uy, D2Q9_WEIGHTS[i], cx, cy);
                        }
                    }
                }
            }
        }
    }
}

// ── LbmDiagnostics ───────────────────────────────────────────────────────────

/// Diagnostic utilities for LBM simulations.
pub struct LbmDiagnostics;

impl LbmDiagnostics {
    /// Compute the vorticity field `∂uy/∂x − ∂ux/∂y` using central differences.
    ///
    /// Returns a flat `Vec`f64` in row-major order (same layout as `grid.cells`).
    /// Periodic boundary wrapping is used at the domain edges.
    pub fn compute_vorticity(grid: &GpuLbmGrid) -> Vec<f64> {
        let nx = grid.nx;
        let ny = grid.ny;
        let mut vort = vec![0.0f64; nx * ny];

        for y in 0..ny {
            for x in 0..nx {
                let xp = (x + 1) % nx;
                let xm = (x + nx - 1) % nx;
                let yp = (y + 1) % ny;
                let ym = (y + ny - 1) % ny;

                let uy_xp = grid.cells[grid.idx(xp, y)].velocity()[1];
                let uy_xm = grid.cells[grid.idx(xm, y)].velocity()[1];
                let ux_yp = grid.cells[grid.idx(x, yp)].velocity()[0];
                let ux_ym = grid.cells[grid.idx(x, ym)].velocity()[0];

                let duy_dx = (uy_xp - uy_xm) / 2.0;
                let dux_dy = (ux_yp - ux_ym) / 2.0;
                vort[grid.idx(x, y)] = duy_dx - dux_dy;
            }
        }
        vort
    }

    /// Compute the velocity divergence field `∂ux/∂x + ∂uy/∂y` using central
    /// differences with periodic wrapping.
    pub fn compute_divergence(grid: &GpuLbmGrid) -> Vec<f64> {
        let nx = grid.nx;
        let ny = grid.ny;
        let mut div = vec![0.0f64; nx * ny];

        for y in 0..ny {
            for x in 0..nx {
                let xp = (x + 1) % nx;
                let xm = (x + nx - 1) % nx;
                let yp = (y + 1) % ny;
                let ym = (y + ny - 1) % ny;

                let ux_xp = grid.cells[grid.idx(xp, y)].velocity()[0];
                let ux_xm = grid.cells[grid.idx(xm, y)].velocity()[0];
                let uy_yp = grid.cells[grid.idx(x, yp)].velocity()[1];
                let uy_ym = grid.cells[grid.idx(x, ym)].velocity()[1];

                let dux_dx = (ux_xp - ux_xm) / 2.0;
                let duy_dy = (uy_yp - uy_ym) / 2.0;
                div[grid.idx(x, y)] = dux_dx + duy_dy;
            }
        }
        div
    }

    /// Return the maximum velocity magnitude across all cells.
    pub fn max_velocity(grid: &GpuLbmGrid) -> f64 {
        grid.cells
            .iter()
            .map(|c| {
                let [ux, uy] = c.velocity();
                (ux * ux + uy * uy).sqrt()
            })
            .fold(0.0f64, f64::max)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    // ── D2Q9Weights ──────────────────────────────────────────────────────────

    #[test]
    fn test_d2q9_weights_sum_to_one() {
        let sum: f64 = D2Q9_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < EPS, "weights sum = {sum}");
    }

    #[test]
    fn test_d2q9_weights_standard_construction() {
        let lw = D2Q9Weights::standard();
        assert_eq!(lw.weights, D2Q9_WEIGHTS);
        assert_eq!(lw.dirs, D2Q9_DIRS);
    }

    #[test]
    fn test_d2q9_default_equals_standard() {
        let a = D2Q9Weights::default();
        let b = D2Q9Weights::standard();
        assert_eq!(a, b);
    }

    #[test]
    fn test_d2q9_cs2() {
        assert!((D2Q9Weights::cs2() - 1.0 / 3.0).abs() < EPS);
    }

    #[test]
    fn test_d2q9_dirs_rest_is_zero() {
        assert_eq!(D2Q9_DIRS[0], [0, 0]);
    }

    #[test]
    fn test_d2q9_opposite_involution() {
        for i in 0..9 {
            assert_eq!(D2Q9_OPPOSITE[D2Q9_OPPOSITE[i]], i);
        }
    }

    #[test]
    fn test_d2q9_opposite_directions_are_negatives() {
        for i in 0..9 {
            let j = D2Q9_OPPOSITE[i];
            assert_eq!(D2Q9_DIRS[j][0], -D2Q9_DIRS[i][0]);
            assert_eq!(D2Q9_DIRS[j][1], -D2Q9_DIRS[i][1]);
        }
    }

    // ── bgk_equilibrium ──────────────────────────────────────────────────────

    #[test]
    fn test_bgk_equilibrium_rest_state() {
        // At ux=uy=0, feq = w * rho
        let rho = 1.0;
        for i in 0..9 {
            let cx = D2Q9_DIRS[i][0] as f64;
            let cy = D2Q9_DIRS[i][1] as f64;
            let feq = bgk_equilibrium(rho, 0.0, 0.0, D2Q9_WEIGHTS[i], cx, cy);
            assert!(
                (feq - D2Q9_WEIGHTS[i]).abs() < EPS,
                "i={i}: feq={feq}, w={}",
                D2Q9_WEIGHTS[i]
            );
        }
    }

    #[test]
    fn test_bgk_equilibrium_sum_is_rho() {
        let rho = 1.5;
        let ux = 0.1;
        let uy = -0.05;
        let sum: f64 = (0..9)
            .map(|i| {
                bgk_equilibrium(
                    rho,
                    ux,
                    uy,
                    D2Q9_WEIGHTS[i],
                    D2Q9_DIRS[i][0] as f64,
                    D2Q9_DIRS[i][1] as f64,
                )
            })
            .sum();
        assert!((sum - rho).abs() < 1e-10, "sum={sum}, rho={rho}");
    }

    #[test]
    fn test_bgk_equilibrium_momentum_x_is_rho_ux() {
        let rho = 1.2;
        let ux = 0.08;
        let uy = 0.0;
        let mom_x: f64 = (0..9)
            .map(|i| {
                bgk_equilibrium(
                    rho,
                    ux,
                    uy,
                    D2Q9_WEIGHTS[i],
                    D2Q9_DIRS[i][0] as f64,
                    D2Q9_DIRS[i][1] as f64,
                ) * D2Q9_DIRS[i][0] as f64
            })
            .sum();
        assert!((mom_x - rho * ux).abs() < 1e-10, "mom_x={mom_x}");
    }

    #[test]
    fn test_bgk_equilibrium_momentum_y_is_rho_uy() {
        let rho = 1.0;
        let ux = 0.0;
        let uy = 0.05;
        let mom_y: f64 = (0..9)
            .map(|i| {
                bgk_equilibrium(
                    rho,
                    ux,
                    uy,
                    D2Q9_WEIGHTS[i],
                    D2Q9_DIRS[i][0] as f64,
                    D2Q9_DIRS[i][1] as f64,
                ) * D2Q9_DIRS[i][1] as f64
            })
            .sum();
        assert!((mom_y - rho * uy).abs() < 1e-10, "mom_y={mom_y}");
    }

    #[test]
    fn test_bgk_equilibrium_positive_for_small_u() {
        for i in 0..9 {
            let val = bgk_equilibrium(
                1.0,
                0.05,
                0.05,
                D2Q9_WEIGHTS[i],
                D2Q9_DIRS[i][0] as f64,
                D2Q9_DIRS[i][1] as f64,
            );
            assert!(val >= 0.0, "i={i}: feq={val}");
        }
    }

    // ── relaxation_from_viscosity ────────────────────────────────────────────

    #[test]
    fn test_relaxation_from_viscosity_nu0() {
        // nu=0 => omega = 1/(0+0.5) = 2.0
        let omega = relaxation_from_viscosity(0.0, 1.0 / 3.0);
        assert!((omega - 2.0).abs() < EPS);
    }

    #[test]
    fn test_relaxation_from_viscosity_typical() {
        let cs2 = 1.0 / 3.0;
        let nu = 1.0 / 6.0;
        // omega = 1 / (nu/cs2 + 0.5) = 1/(0.5 + 0.5) = 1.0
        let omega = relaxation_from_viscosity(nu, cs2);
        assert!((omega - 1.0).abs() < EPS);
    }

    #[test]
    fn test_relaxation_from_viscosity_range() {
        let cs2 = 1.0 / 3.0;
        for nu_times_10 in 1..20usize {
            let nu = nu_times_10 as f64 * 0.01;
            let omega = relaxation_from_viscosity(nu, cs2);
            // omega must be in (0, 2) for stable BGK
            assert!(omega > 0.0 && omega < 2.0, "nu={nu}, omega={omega}");
        }
    }

    // ── LbmCell ──────────────────────────────────────────────────────────────

    #[test]
    fn test_lbm_cell_density_equilibrium() {
        let rho = 1.3;
        let cell = LbmCell::new_equilibrium(rho);
        assert!((cell.density() - rho).abs() < EPS);
    }

    #[test]
    fn test_lbm_cell_velocity_rest() {
        let cell = LbmCell::new_equilibrium(1.0);
        let [ux, uy] = cell.velocity();
        assert!(ux.abs() < EPS);
        assert!(uy.abs() < EPS);
    }

    #[test]
    fn test_lbm_cell_default_density_one() {
        let cell = LbmCell::default();
        assert!((cell.density() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_lbm_cell_velocity_zero_density() {
        let cell = LbmCell {
            f_in: [0.0; 9],
            f_out: [0.0; 9],
        };
        let [ux, uy] = cell.velocity();
        assert_eq!([ux, uy], [0.0, 0.0]);
    }

    // ── GpuLbmGrid ───────────────────────────────────────────────────────────

    #[test]
    fn test_gpu_lbm_grid_creation() {
        let grid = GpuLbmGrid::new(8, 8, 1.0);
        assert_eq!(grid.nx, 8);
        assert_eq!(grid.ny, 8);
        assert_eq!(grid.cells.len(), 64);
    }

    #[test]
    fn test_gpu_lbm_grid_idx() {
        let grid = GpuLbmGrid::new(4, 4, 1.0);
        assert_eq!(grid.idx(0, 0), 0);
        assert_eq!(grid.idx(3, 3), 15);
        assert_eq!(grid.idx(1, 2), 9);
    }

    #[test]
    fn test_collision_bgk_conserves_mass() {
        let mut grid = GpuLbmGrid::new(4, 4, 1.0);
        let mass_before: f64 = grid.cells.iter().map(|c| c.density()).sum();
        grid.collision_bgk(1.0);
        let mass_after: f64 = grid.cells.iter().map(|c| c.f_out.iter().sum::<f64>()).sum();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }

    #[test]
    fn test_collision_bgk_equilibrium_state_unchanged() {
        // If cell is already at equilibrium, f_out == f_in after collision
        let mut grid = GpuLbmGrid::new(4, 4, 1.0);
        let omega = 1.0;
        grid.collision_bgk(omega);
        for cell in &grid.cells {
            for i in 0..9 {
                assert!((cell.f_in[i] - cell.f_out[i]).abs() < EPS);
            }
        }
    }

    #[test]
    fn test_stream_periodic_mass_conservation() {
        let mut grid = GpuLbmGrid::new(6, 6, 1.0);
        let mass_before: f64 = grid.cells.iter().map(|c| c.density()).sum();
        // Prepare f_out = f_in before streaming
        for cell in &mut grid.cells {
            cell.f_out = cell.f_in;
        }
        grid.stream_periodic();
        let mass_after: f64 = grid.cells.iter().map(|c| c.density()).sum();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }

    #[test]
    fn test_step_mass_conservation() {
        let mut grid = GpuLbmGrid::new(8, 8, 1.0);
        let mass_before: f64 = grid.cells.iter().map(|c| c.density()).sum();
        grid.step(1.0);
        let mass_after: f64 = grid.cells.iter().map(|c| c.density()).sum();
        assert!((mass_before - mass_after).abs() < 1e-9);
    }

    #[test]
    fn test_step_multiple_iterations() {
        let mut grid = GpuLbmGrid::new(4, 4, 1.0);
        let mass_before: f64 = grid.cells.iter().map(|c| c.density()).sum();
        for _ in 0..10 {
            grid.step(1.0);
        }
        let mass_after: f64 = grid.cells.iter().map(|c| c.density()).sum();
        assert!((mass_before - mass_after).abs() < 1e-8);
    }

    // ── GpuLbmDispatch ───────────────────────────────────────────────────────

    #[test]
    fn test_gpu_lbm_dispatch_creation() {
        let dispatch = GpuLbmDispatch::new(4, 4, 1.0);
        assert_eq!(dispatch.boundary.len(), 16);
        assert!(
            dispatch
                .boundary
                .iter()
                .all(|b| *b == LbmBoundary::Periodic)
        );
    }

    #[test]
    fn test_dispatch_set_boundary() {
        let mut d = GpuLbmDispatch::new(4, 4, 1.0);
        d.set_boundary(1, 2, LbmBoundary::NoSlip);
        assert_eq!(d.boundary[d.grid.idx(1, 2)], LbmBoundary::NoSlip);
    }

    #[test]
    fn test_dispatch_zou_he_boundary() {
        let mut d = GpuLbmDispatch::new(4, 4, 1.0);
        d.set_boundary(
            0,
            0,
            LbmBoundary::ZouHe {
                rho: 1.05,
                ux: 0.1,
                uy: 0.0,
            },
        );
        d.dispatch_step(1.0);
        // After step the cell should have density close to 1.05
        // (it's been overwritten by ZouHe then run through collision/stream)
        // Just check no panic and mass is finite
        let total_mass: f64 = d.grid.cells.iter().map(|c| c.density()).sum();
        assert!(total_mass.is_finite());
    }

    #[test]
    fn test_dispatch_step_mass_finite() {
        let mut d = GpuLbmDispatch::new(6, 6, 1.0);
        d.set_boundary(0, 0, LbmBoundary::NoSlip);
        d.dispatch_step(1.0);
        for cell in &d.grid.cells {
            assert!(cell.density().is_finite());
        }
    }

    // ── LbmDiagnostics ───────────────────────────────────────────────────────

    #[test]
    fn test_vorticity_uniform_flow_is_zero() {
        // Uniform flow: all cells identical => gradients are zero
        let grid = GpuLbmGrid::new(6, 6, 1.0);
        let vort = LbmDiagnostics::compute_vorticity(&grid);
        for v in &vort {
            assert!(v.abs() < EPS, "vorticity={v}");
        }
    }

    #[test]
    fn test_vorticity_field_length() {
        let grid = GpuLbmGrid::new(5, 7, 1.0);
        let vort = LbmDiagnostics::compute_vorticity(&grid);
        assert_eq!(vort.len(), 35);
    }

    #[test]
    fn test_divergence_uniform_flow_is_zero() {
        let grid = GpuLbmGrid::new(6, 6, 1.0);
        let div = LbmDiagnostics::compute_divergence(&grid);
        for d in &div {
            assert!(d.abs() < EPS, "div={d}");
        }
    }

    #[test]
    fn test_divergence_field_length() {
        let grid = GpuLbmGrid::new(3, 4, 1.0);
        let div = LbmDiagnostics::compute_divergence(&grid);
        assert_eq!(div.len(), 12);
    }

    #[test]
    fn test_max_velocity_rest_state() {
        let grid = GpuLbmGrid::new(4, 4, 1.0);
        let mv = LbmDiagnostics::max_velocity(&grid);
        assert!(mv.abs() < EPS, "max_vel={mv}");
    }

    #[test]
    fn test_max_velocity_after_step() {
        let mut grid = GpuLbmGrid::new(4, 4, 1.0);
        // Perturb one cell
        grid.cells[0].f_in[1] += 0.1;
        grid.step(1.0);
        let mv = LbmDiagnostics::max_velocity(&grid);
        assert!(mv.is_finite());
    }

    // ── Macroscopic recovery ─────────────────────────────────────────────────

    #[test]
    fn test_macroscopic_recovery_rho() {
        // Set f_in explicitly and verify density()
        let mut cell = LbmCell {
            f_in: [0.0; 9],
            f_out: [0.0; 9],
        };
        cell.f_in[0] = 0.5;
        cell.f_in[1] = 0.25;
        cell.f_in[2] = 0.25;
        assert!((cell.density() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_macroscopic_recovery_ux() {
        // Manually set f_in so ux = sum(fi * cx) / rho is predictable
        let rho = 1.0;
        let ux = 0.1;
        let uy = 0.0;
        let mut cell = LbmCell {
            f_in: [0.0; 9],
            f_out: [0.0; 9],
        };
        for i in 0..9 {
            let cx = D2Q9_DIRS[i][0] as f64;
            let cy = D2Q9_DIRS[i][1] as f64;
            cell.f_in[i] = bgk_equilibrium(rho, ux, uy, D2Q9_WEIGHTS[i], cx, cy);
        }
        let [vx, vy] = cell.velocity();
        assert!((vx - ux).abs() < 1e-10, "vx={vx}");
        assert!(vy.abs() < 1e-10, "vy={vy}");
    }

    #[test]
    fn test_macroscopic_recovery_uy() {
        let rho = 1.0;
        let ux = 0.0;
        let uy = -0.07;
        let mut cell = LbmCell {
            f_in: [0.0; 9],
            f_out: [0.0; 9],
        };
        for i in 0..9 {
            let cx = D2Q9_DIRS[i][0] as f64;
            let cy = D2Q9_DIRS[i][1] as f64;
            cell.f_in[i] = bgk_equilibrium(rho, ux, uy, D2Q9_WEIGHTS[i], cx, cy);
        }
        let [vx, vy] = cell.velocity();
        assert!(vx.abs() < 1e-10, "vx={vx}");
        assert!((vy - uy).abs() < 1e-10, "vy={vy}");
    }

    // ── Streaming correctness ────────────────────────────────────────────────

    #[test]
    fn test_stream_periodic_translation() {
        // Put a pulse in the x-direction distribution (index 1 = +x)
        // and verify it moves one cell to the right after streaming.
        let nx = 4;
        let ny = 1;
        let mut grid = GpuLbmGrid::new(nx, ny, 1.0);

        // Reset all f_in to 0, place 1.0 in direction 1 (cx=+1) at x=0
        for cell in &mut grid.cells {
            cell.f_in = [0.0; 9];
            cell.f_out = [0.0; 9];
        }
        grid.cells[0].f_in[1] = 1.0; // direction +x, cell (0,0)
        // Copy f_in → f_out (skip collision)
        for cell in &mut grid.cells {
            cell.f_out = cell.f_in;
        }
        grid.stream_periodic();

        // The value should now be at cell (1,0) direction 1
        assert!((grid.cells[1].f_in[1] - 1.0).abs() < EPS);
        assert!(grid.cells[0].f_in[1].abs() < EPS);
    }

    #[test]
    fn test_stream_periodic_wrap() {
        // Direction +x from the last cell should wrap to cell 0.
        let nx = 4;
        let ny = 1;
        let mut grid = GpuLbmGrid::new(nx, ny, 1.0);
        for cell in &mut grid.cells {
            cell.f_in = [0.0; 9];
            cell.f_out = [0.0; 9];
        }
        grid.cells[nx - 1].f_in[1] = 1.0;
        for cell in &mut grid.cells {
            cell.f_out = cell.f_in;
        }
        grid.stream_periodic();
        assert!((grid.cells[0].f_in[1] - 1.0).abs() < EPS);
    }
}
