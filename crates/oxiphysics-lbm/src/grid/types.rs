//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[cfg(test)]
#[cfg(test)]
#[allow(unused_imports)]
use super::functions_2::*;
use crate::lattice::CS2;
use crate::lattice::{
    D2Q9_VELOCITIES, D2Q9_WEIGHTS, bgk_d2q9, bounce_back_node_d2q9, equilibrium_d2q9,
    macros_from_d2q9, stream_d2q9_periodic,
};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{C9, C19, OPP9, W9, W19, equilibrium_2d};
use super::grid_extended::{BoundaryNodeType, NodeFlag};

/// Cell-based 3-D LBM grid for D3Q19.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CellGrid3D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
    /// Cells in row-major order: `cells[z * ny * nx + y * nx + x]`.
    pub cells: Vec<LbmCell3D>,
}
#[allow(dead_code)]
impl CellGrid3D {
    /// Create a new 3-D grid initialised to equilibrium with uniform density `rho0`.
    pub fn new(nx: usize, ny: usize, nz: usize, rho0: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            cells: vec![LbmCell3D::at_equilibrium(rho0); nx * ny * nz],
        }
    }
    /// Linear index for cell `(x, y, z)`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Perform one complete BGK step: collision, streaming, macroscopic update.
    pub fn step(&mut self, omega: f64) {
        self.collide(omega);
        self.stream();
        self.update_macroscopic();
    }
    /// BGK collision in-place on fluid cells.
    fn collide(&mut self, omega: f64) {
        for cell in self.cells.iter_mut() {
            if cell.obstacle {
                continue;
            }
            let rho = cell.rho;
            let ux = cell.u[0];
            let uy = cell.u[1];
            let uz = cell.u[2];
            let u_sq = ux * ux + uy * uy + uz * uz;
            for i in 0..19 {
                let cx = C19[i][0] as f64;
                let cy = C19[i][1] as f64;
                let cz = C19[i][2] as f64;
                let eu = cx * ux + cy * uy + cz * uz;
                let feq = W19[i]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                cell.f[i] -= omega * (cell.f[i] - feq);
            }
        }
    }
    /// Periodic pull-streaming for D3Q19.
    fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let f_old: Vec<[f64; 19]> = self.cells.iter().map(|c| c.f).collect();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = z * ny * nx + y * nx + x;
                    for a in 0..19 {
                        let sx = (x as isize - C19[a][0] as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - C19[a][1] as isize).rem_euclid(ny as isize) as usize;
                        let sz = (z as isize - C19[a][2] as isize).rem_euclid(nz as isize) as usize;
                        self.cells[dst].f[a] = f_old[sz * ny * nx + sy * nx + sx][a];
                    }
                }
            }
        }
    }
    /// Recompute macroscopic density and velocity for all fluid cells.
    fn update_macroscopic(&mut self) {
        for cell in self.cells.iter_mut() {
            if cell.obstacle {
                continue;
            }
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for a in 0..19 {
                let fi = cell.f[a];
                rho += fi;
                mx += fi * C19[a][0] as f64;
                my += fi * C19[a][1] as f64;
                mz += fi * C19[a][2] as f64;
            }
            cell.rho = rho;
            if rho.abs() > 1e-15 {
                cell.u = [mx / rho, my / rho, mz / rho];
            } else {
                cell.u = [0.0; 3];
            }
        }
    }
    /// Domain-averaged velocity `[mean_ux, mean_uy, mean_uz]` over fluid cells.
    pub fn mean_velocity(&self) -> [f64; 3] {
        let fluid: Vec<_> = self.cells.iter().filter(|c| !c.obstacle).collect();
        let n = fluid.len() as f64;
        if n < 1.0 {
            return [0.0; 3];
        }
        let sum_x: f64 = fluid.iter().map(|c| c.u[0]).sum();
        let sum_y: f64 = fluid.iter().map(|c| c.u[1]).sum();
        let sum_z: f64 = fluid.iter().map(|c| c.u[2]).sum();
        [sum_x / n, sum_y / n, sum_z / n]
    }
    /// Total mass across all cells.
    pub fn total_mass(&self) -> f64 {
        self.cells.iter().map(|c| c.rho).sum()
    }
}
/// A 2D LBM grid that stores populations per cell as `[f64; 9]` arrays.
///
/// Supports:
/// - BGK collision with spatially uniform ω
/// - Periodic streaming
/// - Bounce-back wall flags
/// - Guo body-force scheme
#[allow(dead_code)]
pub struct CellularGrid2D {
    /// Width in cells.
    pub nx: usize,
    /// Height in cells.
    pub ny: usize,
    /// Population arrays: `pop[y * nx + x]` = `[f_0 .. f_8]`.
    pub pop: Vec<[f64; 9]>,
    /// Wall flags: `wall[y * nx + x]` = true means solid (bounce-back).
    pub wall: Vec<bool>,
    /// BGK relaxation frequency ω = 1/τ.
    pub omega: f64,
}
impl CellularGrid2D {
    /// Create a grid initialised to rest equilibrium (ρ=1, u=0).
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, omega: f64) -> Self {
        let n = nx * ny;
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        Self {
            nx,
            ny,
            pop: vec![feq; n],
            wall: vec![false; n],
            omega,
        }
    }
    /// Linear index for cell `(x, y)`.
    #[inline]
    #[allow(dead_code)]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Mark a cell as solid wall.
    #[allow(dead_code)]
    pub fn set_wall(&mut self, x: usize, y: usize) {
        let i = y * self.nx + x;
        self.wall[i] = true;
    }
    /// Remove wall flag.
    #[allow(dead_code)]
    pub fn clear_wall(&mut self, x: usize, y: usize) {
        let i = y * self.nx + x;
        self.wall[i] = false;
    }
    /// Compute total mass (sum of all populations).
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f64 {
        self.pop.iter().flat_map(|n| n.iter()).sum()
    }
    /// Compute total x-momentum.
    #[allow(dead_code)]
    pub fn total_momentum_x(&self) -> f64 {
        let mut jx = 0.0f64;
        for (idx, node) in self.pop.iter().enumerate() {
            if !self.wall[idx] {
                for q in 0..9 {
                    jx += D2Q9_VELOCITIES[q][0] as f64 * node[q];
                }
            }
        }
        jx
    }
    /// Compute total kinetic energy (½ ρ u²).
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0f64;
        for node in &self.pop {
            let (rho, ux, uy) = macros_from_d2q9(node);
            ke += 0.5 * rho * (ux * ux + uy * uy);
        }
        ke
    }
    /// Get macroscopic variables at cell `(x, y)`.
    #[allow(dead_code)]
    pub fn macros_at(&self, x: usize, y: usize) -> (f64, f64, f64) {
        macros_from_d2q9(&self.pop[self.idx(x, y)])
    }
    /// Apply BGK collision to all fluid cells.
    #[allow(dead_code)]
    pub fn collide(&mut self) {
        for idx in 0..self.pop.len() {
            if !self.wall[idx] {
                let (rho, ux, uy) = macros_from_d2q9(&self.pop[idx]);
                bgk_d2q9(&mut self.pop[idx], rho, ux, uy, self.omega);
            }
        }
    }
    /// Apply BGK collision with Guo body force to all fluid cells.
    ///
    /// Adds the Guo correction `F_i = w_i * (1 - 1/(2τ)) * (e_i - u + (e_i·u) e_i/cs²) · F / cs²`
    #[allow(clippy::too_many_arguments)]
    #[allow(dead_code)]
    pub fn collide_with_force(&mut self, fx: f64, fy: f64) {
        let tau = 1.0 / self.omega;
        let pre = 1.0 - 0.5 / tau;
        for idx in 0..self.pop.len() {
            if !self.wall[idx] {
                let (rho, ux, uy) = macros_from_d2q9(&self.pop[idx]);
                bgk_d2q9(&mut self.pop[idx], rho, ux, uy, self.omega);
                for q in 0..9 {
                    let cx = D2Q9_VELOCITIES[q][0] as f64;
                    let cy = D2Q9_VELOCITIES[q][1] as f64;
                    let eu = cx * ux + cy * uy;
                    let guo = D2Q9_WEIGHTS[q]
                        * pre
                        * ((cx - ux) * fx / CS2
                            + (cy - uy) * fy / CS2
                            + eu * (cx * fx + cy * fy) / (CS2 * CS2));
                    self.pop[idx][q] += guo;
                }
            }
        }
    }
    /// Apply periodic streaming.
    #[allow(dead_code)]
    pub fn stream(&mut self) {
        stream_d2q9_periodic(&mut self.pop, self.nx, self.ny);
    }
    /// Apply half-way bounce-back to all wall cells.
    #[allow(dead_code)]
    pub fn bounce_back(&mut self) {
        for idx in 0..self.pop.len() {
            if self.wall[idx] {
                bounce_back_node_d2q9(&mut self.pop[idx]);
            }
        }
    }
    /// Execute one full LBM step: collide → stream → bounce-back.
    #[allow(dead_code)]
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.bounce_back();
    }
    /// Execute one step with body force.
    #[allow(dead_code)]
    pub fn step_with_force(&mut self, fx: f64, fy: f64) {
        self.collide_with_force(fx, fy);
        self.stream();
        self.bounce_back();
    }
    /// Extract x-velocity profile at column `x` (all y values).
    #[allow(dead_code)]
    pub fn ux_profile(&self, x: usize) -> Vec<f64> {
        (0..self.ny)
            .map(|y| {
                let (_, ux, _) = macros_from_d2q9(&self.pop[self.idx(x, y)]);
                ux
            })
            .collect()
    }
    /// Return max velocity magnitude over all fluid cells.
    #[allow(dead_code)]
    pub fn max_speed(&self) -> f64 {
        self.pop
            .iter()
            .enumerate()
            .filter(|(i, _)| !self.wall[*i])
            .map(|(_, node)| {
                let (_, ux, uy) = macros_from_d2q9(node);
                (ux * ux + uy * uy).sqrt()
            })
            .fold(0.0f64, f64::max)
    }
    /// Initialise all fluid cells to equilibrium at given density and velocity.
    #[allow(dead_code)]
    pub fn set_uniform_equilibrium(&mut self, rho: f64, ux: f64, uy: f64) {
        let feq = equilibrium_d2q9(rho, ux, uy);
        for (idx, node) in self.pop.iter_mut().enumerate() {
            if !self.wall[idx] {
                *node = feq;
            }
        }
    }
}
/// A D2Q9 LBM grid with per-node boundary flags and built-in
/// collide-stream-BC stepping.
///
/// Distribution storage: flat `f[idx * 9 + alpha]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlaggedGrid2D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Relaxation frequency omega = 1/tau.
    pub omega: f64,
    /// Distribution functions (flat, length = nx * ny * 9).
    pub f: Vec<f64>,
    /// Per-node flags.
    pub flags: Vec<NodeFlag>,
}
#[allow(dead_code)]
impl FlaggedGrid2D {
    /// Create a new flagged 2D grid initialised to equilibrium at rest.
    pub fn new(nx: usize, ny: usize, omega: f64, rho0: f64) -> Self {
        let n = nx * ny;
        let mut f = vec![0.0; n * 9];
        for idx in 0..n {
            for alpha in 0..9 {
                f[idx * 9 + alpha] = W9[alpha] * rho0;
            }
        }
        Self {
            nx,
            ny,
            omega,
            f,
            flags: vec![NodeFlag::Fluid; n],
        }
    }
    /// Linear index for cell (x, y).
    #[inline]
    pub(crate) fn cell_idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Mark cell (x, y) as a wall.
    pub fn set_wall(&mut self, x: usize, y: usize) {
        let idx = self.cell_idx(x, y);
        self.flags[idx] = NodeFlag::Wall;
    }
    /// Mark cell (x, y) as an inlet with prescribed macroscopic values.
    pub fn set_inlet(&mut self, x: usize, y: usize, rho: f64, ux: f64, uy: f64) {
        let idx = self.cell_idx(x, y);
        self.flags[idx] = NodeFlag::Inlet {
            rho: rho.to_bits(),
            ux_bits: ux.to_bits(),
            uy_bits: uy.to_bits(),
        };
        for alpha in 0..9 {
            self.f[idx * 9 + alpha] = equilibrium_2d(
                W9[alpha],
                rho,
                ux,
                uy,
                C9[alpha][0] as f64,
                C9[alpha][1] as f64,
            );
        }
    }
    /// Mark cell (x, y) as an outlet.
    pub fn set_outlet(&mut self, x: usize, y: usize) {
        let idx = self.cell_idx(x, y);
        self.flags[idx] = NodeFlag::Outlet;
    }
    /// Perform one full time step: collide, stream, apply boundary conditions.
    pub fn step(&mut self) {
        self.collide_interior();
        self.stream_periodic();
        self.apply_bcs();
    }
    /// BGK collision on fluid nodes only.
    fn collide_interior(&mut self) {
        let n = self.nx * self.ny;
        let omega = self.omega;
        for idx in 0..n {
            if self.flags[idx] != NodeFlag::Fluid {
                continue;
            }
            let base = idx * 9;
            let mut rho = 0.0;
            let mut mx = 0.0;
            let mut my = 0.0;
            for a in 0..9 {
                let fi = self.f[base + a];
                rho += fi;
                mx += fi * C9[a][0] as f64;
                my += fi * C9[a][1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            for a in 0..9 {
                let feq = equilibrium_2d(W9[a], rho, ux, uy, C9[a][0] as f64, C9[a][1] as f64);
                self.f[base + a] -= omega * (self.f[base + a] - feq);
            }
        }
    }
    /// Pull-scheme periodic streaming.
    fn stream_periodic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let dst = (y * nx + x) * 9;
                for a in 0..9 {
                    let sx = (x as isize - C9[a][0] as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - C9[a][1] as isize).rem_euclid(ny as isize) as usize;
                    let src = (sy * nx + sx) * 9 + a;
                    self.f[dst + a] = f_old[src];
                }
            }
        }
    }
    /// Apply boundary conditions for all flagged nodes.
    fn apply_bcs(&mut self) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            match self.flags[idx] {
                NodeFlag::Fluid | NodeFlag::Symmetry => {}
                NodeFlag::Wall => {
                    let base = idx * 9;
                    let mut tmp = [0.0; 9];
                    tmp.copy_from_slice(&self.f[base..base + 9]);
                    for a in 0..9 {
                        self.f[base + a] = tmp[OPP9[a]];
                    }
                }
                NodeFlag::Inlet {
                    rho,
                    ux_bits,
                    uy_bits,
                } => {
                    let r = f64::from_bits(rho);
                    let ux = f64::from_bits(ux_bits);
                    let uy = f64::from_bits(uy_bits);
                    let base = idx * 9;
                    for a in 0..9 {
                        self.f[base + a] =
                            equilibrium_2d(W9[a], r, ux, uy, C9[a][0] as f64, C9[a][1] as f64);
                    }
                }
                NodeFlag::Outlet => {
                    let x = idx % self.nx;
                    let y = idx / self.nx;
                    if x > 0 {
                        let src_idx = y * self.nx + (x - 1);
                        for a in 0..9 {
                            self.f[idx * 9 + a] = self.f[src_idx * 9 + a];
                        }
                    }
                }
            }
        }
    }
    /// Return the density at cell (x, y).
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        let base = self.cell_idx(x, y) * 9;
        let mut rho = 0.0;
        for a in 0..9 {
            rho += self.f[base + a];
        }
        rho
    }
    /// Return the velocity \[ux, uy\] at cell (x, y).
    pub fn velocity_at(&self, x: usize, y: usize) -> [f64; 2] {
        let base = self.cell_idx(x, y) * 9;
        let mut rho = 0.0;
        let mut mx = 0.0;
        let mut my = 0.0;
        for a in 0..9 {
            let fi = self.f[base + a];
            rho += fi;
            mx += fi * C9[a][0] as f64;
            my += fi * C9[a][1] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho]
        } else {
            [0.0, 0.0]
        }
    }
    /// Total mass across all cells.
    pub fn total_mass(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut total = 0.0;
        for idx in 0..n {
            let base = idx * 9;
            for a in 0..9 {
                total += self.f[base + a];
            }
        }
        total
    }
}
/// A single 2-D LBM cell holding D2Q9 distributions and macroscopic fields.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct LbmCell2D {
    /// D2Q9 distribution functions.
    pub f: [f64; 9],
    /// Macroscopic density.
    pub rho: f64,
    /// Macroscopic x-velocity.
    pub ux: f64,
    /// Macroscopic y-velocity.
    pub uy: f64,
    /// If `true` this cell is a solid obstacle (bounce-back).
    pub obstacle: bool,
}
impl LbmCell2D {
    /// Create a cell at equilibrium for density `rho` and zero velocity.
    fn at_equilibrium(rho: f64) -> Self {
        let feq = {
            let u_sq = 0.0_f64;
            let mut f = [0.0_f64; 9];
            for i in 0..9 {
                let eu = 0.0_f64;
                f[i] = W9[i]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
            }
            f
        };
        Self {
            f: feq,
            rho,
            ux: 0.0,
            uy: 0.0,
            obstacle: false,
        }
    }
}
/// Cell-based 2-D LBM grid for D2Q9.
///
/// Stores one `LbmCell2D` per lattice node.  The `step` method performs a
/// complete BGK collision + pull-streaming + macroscopic update.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CellGrid2D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Cells in row-major order: `cells[y * nx + x]`.
    pub cells: Vec<LbmCell2D>,
}
#[allow(dead_code)]
impl CellGrid2D {
    /// Create a new grid initialised to equilibrium with uniform density `rho0`.
    pub fn new(nx: usize, ny: usize, rho0: f64) -> Self {
        Self {
            nx,
            ny,
            cells: vec![LbmCell2D::at_equilibrium(rho0); nx * ny],
        }
    }
    /// Linear index for cell `(x, y)`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Perform one complete BGK step: collision, streaming, macroscopic update.
    pub fn step(&mut self, omega: f64) {
        self.collide(omega);
        self.stream();
        self.update_macroscopic();
    }
    /// BGK collision (in-place, fluid cells only).
    fn collide(&mut self, omega: f64) {
        for cell in self.cells.iter_mut() {
            if cell.obstacle {
                continue;
            }
            let rho = cell.rho;
            let ux = cell.ux;
            let uy = cell.uy;
            let u_sq = ux * ux + uy * uy;
            for i in 0..9 {
                let cx = C9[i][0] as f64;
                let cy = C9[i][1] as f64;
                let eu = cx * ux + cy * uy;
                let feq = W9[i]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                cell.f[i] -= omega * (cell.f[i] - feq);
            }
        }
    }
    /// Periodic pull-streaming.
    fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old: Vec<[f64; 9]> = self.cells.iter().map(|c| c.f).collect();
        for y in 0..ny {
            for x in 0..nx {
                let dst = y * nx + x;
                for a in 0..9 {
                    let sx = (x as isize - C9[a][0] as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - C9[a][1] as isize).rem_euclid(ny as isize) as usize;
                    self.cells[dst].f[a] = f_old[sy * nx + sx][a];
                }
            }
        }
    }
    /// Recompute macroscopic density and velocity for all fluid cells.
    fn update_macroscopic(&mut self) {
        for cell in self.cells.iter_mut() {
            if cell.obstacle {
                continue;
            }
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            for a in 0..9 {
                let fi = cell.f[a];
                rho += fi;
                mx += fi * C9[a][0] as f64;
                my += fi * C9[a][1] as f64;
            }
            cell.rho = rho;
            if rho.abs() > 1e-15 {
                cell.ux = mx / rho;
                cell.uy = my / rho;
            } else {
                cell.ux = 0.0;
                cell.uy = 0.0;
            }
        }
    }
    /// Apply full-way bounce-back to all obstacle cells.
    pub fn apply_bounce_back(&mut self) {
        for cell in self.cells.iter_mut() {
            if !cell.obstacle {
                continue;
            }
            let mut tmp = [0.0_f64; 9];
            tmp.copy_from_slice(&cell.f);
            for a in 0..9 {
                cell.f[a] = tmp[OPP9[a]];
            }
        }
    }
    /// Apply a simplified Zou-He inlet BC on the left wall (`x = 0`).
    ///
    /// Prescribes density `rho_in` and y-velocity `uy` at all cells with `x = 0`.
    /// The x-velocity is computed from the Zou-He relations.
    pub fn apply_zou_he_inlet(&mut self, rho_in: f64, uy: f64) {
        for y in 0..self.ny {
            let k = self.idx(0, y);
            let cell = &mut self.cells[k];
            let f3 = cell.f[3];
            let f6 = cell.f[6];
            let f7 = cell.f[7];
            let f0 = cell.f[0];
            let f2 = cell.f[2];
            let f4 = cell.f[4];
            let ux = 1.0 - (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / rho_in;
            let u_sq = ux * ux + uy * uy;
            for i in 0..9 {
                let cx = C9[i][0] as f64;
                let cy = C9[i][1] as f64;
                let eu = cx * ux + cy * uy;
                let feq = W9[i]
                    * rho_in
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                cell.f[i] = feq;
            }
            cell.rho = rho_in;
            cell.ux = ux;
            cell.uy = uy;
        }
    }
    /// Domain-averaged velocity `(mean_ux, mean_uy)` over all fluid cells.
    pub fn mean_velocity(&self) -> (f64, f64) {
        let fluid: Vec<_> = self.cells.iter().filter(|c| !c.obstacle).collect();
        let n = fluid.len() as f64;
        if n < 1.0 {
            return (0.0, 0.0);
        }
        let sum_ux: f64 = fluid.iter().map(|c| c.ux).sum();
        let sum_uy: f64 = fluid.iter().map(|c| c.uy).sum();
        (sum_ux / n, sum_uy / n)
    }
    /// Total mass for conservation checks.
    pub fn total_mass(&self) -> f64 {
        self.cells.iter().map(|c| c.rho).sum()
    }
}
/// Convenience wrapper for a 2D channel flow driven by a uniform body force.
///
/// The channel has no-slip walls at y=0 and y=ny-1 and is periodic in x.
#[allow(dead_code)]
pub struct ChannelFlow2D {
    /// The underlying 2D cellular grid.
    pub grid: CellularGrid2D,
    /// Body force per unit mass in x-direction.
    pub fx: f64,
}
impl ChannelFlow2D {
    /// Construct a channel with no-slip top and bottom walls.
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, omega: f64, fx: f64) -> Self {
        let mut grid = CellularGrid2D::new(nx, ny, omega);
        for x in 0..nx {
            grid.set_wall(x, 0);
            grid.set_wall(x, ny - 1);
        }
        Self { grid, fx }
    }
    /// Advance the channel by one LBM time step with body force.
    #[allow(dead_code)]
    pub fn step(&mut self) {
        self.grid.step_with_force(self.fx, 0.0);
    }
    /// Compute mean x-velocity over interior cells.
    #[allow(dead_code)]
    pub fn mean_ux(&self) -> f64 {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for y in 1..ny - 1 {
            for x in 0..nx {
                let (_, ux, _) = self.grid.macros_at(x, y);
                sum += ux;
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
    /// Analytical Poiseuille centreline velocity: `u_max = fx * H² / (8 ν)`.
    #[allow(dead_code)]
    pub fn poiseuille_umax(&self) -> f64 {
        let h = (self.grid.ny - 2) as f64;
        let nu = CS2 * (1.0 / self.grid.omega - 0.5);
        self.fx * h * h / (8.0 * nu)
    }
}
/// A single 3-D LBM cell holding D3Q19 distributions and macroscopic fields.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct LbmCell3D {
    /// D3Q19 distribution functions.
    pub f: [f64; 19],
    /// Macroscopic density.
    pub rho: f64,
    /// Macroscopic velocity.
    pub u: [f64; 3],
    /// If `true` this cell is a solid obstacle (bounce-back).
    pub obstacle: bool,
}
impl LbmCell3D {
    /// Create a cell at equilibrium for density `rho` and zero velocity.
    fn at_equilibrium(rho: f64) -> Self {
        let mut f = [0.0_f64; 19];
        for i in 0..19 {
            f[i] = W19[i] * rho;
        }
        Self {
            f,
            rho,
            u: [0.0; 3],
            obstacle: false,
        }
    }
}
/// Full 2-D LBM grid using D2Q9 with typed boundary nodes, BGK collision,
/// pull-scheme streaming, and Zou-He inlet/outlet BCs.
///
/// Distribution layout: flat `f2[cell_idx * 9 + alpha]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LbmGrid2DFull {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Relaxation frequency ω = 1/τ.
    pub omega: f64,
    /// Distribution functions: flat, length `nx * ny * 9`.
    pub f2: Vec<f64>,
    /// Per-node boundary type.
    pub node_type: Vec<BoundaryNodeType>,
}
#[allow(dead_code)]
impl LbmGrid2DFull {
    /// Create a new `LbmGrid2DFull` initialised to equilibrium at rest with `rho0`.
    pub fn new(nx: usize, ny: usize, omega: f64, rho0: f64) -> Self {
        let n = nx * ny;
        let mut f2 = vec![0.0_f64; n * 9];
        for idx in 0..n {
            for alpha in 0..9 {
                f2[idx * 9 + alpha] = W9[alpha] * rho0;
            }
        }
        Self {
            nx,
            ny,
            omega,
            f2,
            node_type: vec![BoundaryNodeType::Fluid; n],
        }
    }
    /// Linear index for cell (x, y).
    #[inline]
    pub fn cell_idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Set cell (x, y) as a solid wall.
    pub fn set_solid(&mut self, x: usize, y: usize) {
        let idx = self.cell_idx(x, y);
        self.node_type[idx] = BoundaryNodeType::Solid;
    }
    /// Set cell (x, y) as an inlet with prescribed (rho, ux, uy).
    pub fn set_inlet(&mut self, x: usize, y: usize, rho: f64, ux: f64, uy: f64) {
        let idx = self.cell_idx(x, y);
        self.node_type[idx] = BoundaryNodeType::Inlet { ux, uy, rho };
        for alpha in 0..9 {
            self.f2[idx * 9 + alpha] = equilibrium_2d(
                W9[alpha],
                rho,
                ux,
                uy,
                C9[alpha][0] as f64,
                C9[alpha][1] as f64,
            );
        }
    }
    /// Set cell (x, y) as an outlet.
    pub fn set_outlet(&mut self, x: usize, y: usize) {
        let idx = self.cell_idx(x, y);
        self.node_type[idx] = BoundaryNodeType::Outlet;
    }
    /// Perform one complete time step: collide → stream → apply_boundary_conditions.
    pub fn step(&mut self) {
        self.collide_bgk();
        self.stream_periodic();
        self.apply_boundary_conditions();
    }
    /// BGK collision on fluid nodes only.
    fn collide_bgk(&mut self) {
        let n = self.nx * self.ny;
        let omega = self.omega;
        for idx in 0..n {
            if self.node_type[idx] != BoundaryNodeType::Fluid {
                continue;
            }
            let base = idx * 9;
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            for a in 0..9 {
                let fi = self.f2[base + a];
                rho += fi;
                mx += fi * C9[a][0] as f64;
                my += fi * C9[a][1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            for a in 0..9 {
                let feq = equilibrium_2d(W9[a], rho, ux, uy, C9[a][0] as f64, C9[a][1] as f64);
                self.f2[base + a] -= omega * (self.f2[base + a] - feq);
            }
        }
    }
    /// Pull-scheme periodic streaming.
    fn stream_periodic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f2.clone();
        for y in 0..ny {
            for x in 0..nx {
                let dst = (y * nx + x) * 9;
                for a in 0..9 {
                    let sx = (x as isize - C9[a][0] as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - C9[a][1] as isize).rem_euclid(ny as isize) as usize;
                    self.f2[dst + a] = f_old[(sy * nx + sx) * 9 + a];
                }
            }
        }
    }
    /// Apply boundary conditions for all flagged nodes.
    ///
    /// - `Solid`: full-way bounce-back.
    /// - `Inlet`: Zou-He left-wall velocity BC.
    /// - `Outlet`: zero-gradient (first-order extrapolation).
    /// - `Fluid`: no action.
    pub fn apply_boundary_conditions(&mut self) {
        let n = self.nx * self.ny;
        let nx = self.nx;
        for idx in 0..n {
            match self.node_type[idx] {
                BoundaryNodeType::Fluid => {}
                BoundaryNodeType::Solid => {
                    let base = idx * 9;
                    let mut tmp = [0.0_f64; 9];
                    tmp.copy_from_slice(&self.f2[base..base + 9]);
                    for a in 0..9 {
                        self.f2[base + a] = tmp[OPP9[a]];
                    }
                }
                BoundaryNodeType::Inlet { ux, uy, rho } => {
                    let base = idx * 9;
                    for a in 0..9 {
                        self.f2[base + a] =
                            equilibrium_2d(W9[a], rho, ux, uy, C9[a][0] as f64, C9[a][1] as f64);
                    }
                }
                BoundaryNodeType::Outlet => {
                    let x = idx % nx;
                    let y = idx / nx;
                    if x > 0 {
                        let src = y * nx + (x - 1);
                        for a in 0..9 {
                            self.f2[idx * 9 + a] = self.f2[src * 9 + a];
                        }
                    }
                }
            }
        }
    }
    /// Extract macroscopic density at cell (x, y).
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        let base = self.cell_idx(x, y) * 9;
        (0..9).map(|a| self.f2[base + a]).sum()
    }
    /// Extract macroscopic velocity \[ux, uy\] at cell (x, y).
    pub fn velocity_at(&self, x: usize, y: usize) -> [f64; 2] {
        let base = self.cell_idx(x, y) * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for a in 0..9 {
            let fi = self.f2[base + a];
            rho += fi;
            mx += fi * C9[a][0] as f64;
            my += fi * C9[a][1] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho]
        } else {
            [0.0; 2]
        }
    }
    /// Total mass (for conservation checks).
    pub fn total_mass(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut total = 0.0_f64;
        for idx in 0..n {
            let base = idx * 9;
            for a in 0..9 {
                total += self.f2[base + a];
            }
        }
        total
    }
    /// Apply a uniform body force by adding a momentum source to all fluid cells.
    ///
    /// Adds `F_x` and `F_y` to the first moments at each fluid node:
    ///
    /// ```text
    /// f_i += w_i * (c_ix * Fx + c_iy * Fy) / cs²
    /// ```
    pub fn apply_body_force(&mut self, fx: f64, fy: f64) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            if self.node_type[idx] != BoundaryNodeType::Fluid {
                continue;
            }
            let base = idx * 9;
            for a in 0..9 {
                let cx = C9[a][0] as f64;
                let cy = C9[a][1] as f64;
                self.f2[base + a] += W9[a] * (cx * fx + cy * fy) / CS2;
            }
        }
    }
}
