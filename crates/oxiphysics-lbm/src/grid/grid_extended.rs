//! Extended grid types: LbmGrid2D, LbmGrid3D, FlatLbmGrid2D/3D, FlaggedGrid3D,
//! FullGrid3D, TrtGrid2D, CellularGrid3D, BoundaryNodeType, NodeFlag.

#![allow(clippy::needless_range_loop)]
use crate::lattice::{CS2, Lattice, LatticeType};
use crate::lattice::{
    bgk_d3q19, bounce_back_node_d2q9, bounce_back_node_d3q19, equilibrium_d2q9, equilibrium_d3q19,
    macros_from_d2q9, macros_from_d3q19, stream_d2q9_periodic, stream_d3q19_periodic,
    trt_collision_d2q9, trt_magic_omega_anti,
};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{C9, C19, OPP9, OPP19, W9, W19, equilibrium_2d, equilibrium_3d};

/// 3D LBM grid for D3Q19 or D3Q27 simulations.
///
/// Distribution functions are stored as `f[q][z * ny * nx + y * nx + x]`.
#[derive(Debug, Clone)]
pub struct LbmGrid3D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
    /// Lattice descriptor.
    pub lattice: Lattice,
    /// Distribution functions: `f[q][z * ny * nx + y * nx + x]`.
    pub f: Vec<Vec<f64>>,
    /// Macroscopic density.
    pub rho: Vec<f64>,
    /// Macroscopic x-velocity.
    pub ux: Vec<f64>,
    /// Macroscopic y-velocity.
    pub uy: Vec<f64>,
    /// Macroscopic z-velocity.
    pub uz: Vec<f64>,
}
impl LbmGrid3D {
    /// Create a new 3D grid initialized to equilibrium at rest with `rho=1`.
    ///
    /// # Panics
    ///
    /// Panics if `lattice_type` is `D2Q9`.
    pub fn new(nx: usize, ny: usize, nz: usize, lattice_type: LatticeType) -> Self {
        assert_ne!(
            lattice_type,
            LatticeType::D2Q9,
            "LbmGrid3D requires 3D lattice"
        );
        let lattice = Lattice::new(lattice_type);
        let n = nx * ny * nz;
        let q = lattice.q();
        let rho = vec![1.0; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let uz = vec![0.0; n];
        let mut f = Vec::with_capacity(q);
        for i in 0..q {
            let w = lattice.weight(i);
            f.push(vec![w; n]);
        }
        Self {
            nx,
            ny,
            nz,
            lattice,
            f,
            rho,
            ux,
            uy,
            uz,
        }
    }
    /// Linear index for cell (x, y, z).
    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Recompute macroscopic density and velocity from distributions.
    pub fn compute_macroscopic(&mut self) {
        let q = self.lattice.q();
        let n = self.nx * self.ny * self.nz;
        for k in 0..n {
            let mut rho_k = 0.0;
            let mut ux_k = 0.0;
            let mut uy_k = 0.0;
            let mut uz_k = 0.0;
            for i in 0..q {
                let fi = self.f[i][k];
                rho_k += fi;
                let c = self.lattice.velocity_3d(i);
                ux_k += fi * c[0] as f64;
                uy_k += fi * c[1] as f64;
                uz_k += fi * c[2] as f64;
            }
            self.rho[k] = rho_k;
            if rho_k.abs() > 1e-15 {
                self.ux[k] = ux_k / rho_k;
                self.uy[k] = uy_k / rho_k;
                self.uz[k] = uz_k / rho_k;
            } else {
                self.ux[k] = 0.0;
                self.uy[k] = 0.0;
                self.uz[k] = 0.0;
            }
        }
    }
    /// Return the density at cell (x, y, z).
    pub fn density_at(&self, x: usize, y: usize, z: usize) -> f64 {
        self.rho[self.idx(x, y, z)]
    }
    /// Return the velocity (ux, uy, uz) at cell (x, y, z).
    pub fn velocity_at(&self, x: usize, y: usize, z: usize) -> (f64, f64, f64) {
        let k = self.idx(x, y, z);
        (self.ux[k], self.uy[k], self.uz[k])
    }
    /// Set macroscopic fields and reinitialize distributions to equilibrium.
    #[allow(clippy::too_many_arguments)]
    pub fn set_equilibrium(
        &mut self,
        x: usize,
        y: usize,
        z: usize,
        rho: f64,
        ux: f64,
        uy: f64,
        uz: f64,
    ) {
        let k = self.idx(x, y, z);
        self.rho[k] = rho;
        self.ux[k] = ux;
        self.uy[k] = uy;
        self.uz[k] = uz;
        let q = self.lattice.q();
        for i in 0..q {
            let w = self.lattice.weight(i);
            let c = self.lattice.velocity_3d(i);
            self.f[i][k] =
                equilibrium_3d(w, rho, ux, uy, uz, c[0] as f64, c[1] as f64, c[2] as f64);
        }
    }
    /// Total density across all cells (for conservation checks).
    pub fn total_density(&self) -> f64 {
        self.rho.iter().sum()
    }
    /// Perform a single collide-stream-macroscopic step (placeholder).
    ///
    /// This is a simplified step for the stub; a real simulation would use
    /// the full collision and streaming modules.
    #[allow(dead_code)]
    pub fn step_placeholder(&mut self, _omega: f64) {}
}
/// 3D LBM grid (D3Q19) with flat distribution storage and tau-based collision.
///
/// Distribution layout: `f[(z * ny * nx + y * nx + x) * 19 + alpha]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlatLbmGrid3D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
    /// Distribution functions: flat, length `nx * ny * nz * 19`.
    pub f: Vec<f64>,
    /// Relaxation time tau.
    pub tau: f64,
}
#[allow(dead_code)]
impl FlatLbmGrid3D {
    /// Create a new 3D grid initialized to equilibrium at rest with `rho = 1`.
    pub fn new(nx: usize, ny: usize, nz: usize, tau: f64) -> Self {
        let n = nx * ny * nz;
        let mut f = vec![0.0_f64; n * 19];
        for idx in 0..n {
            for alpha in 0..19 {
                f[idx * 19 + alpha] = W19[alpha];
            }
        }
        Self { nx, ny, nz, f, tau }
    }
    /// Linear cell index for (x, y, z).
    #[inline]
    pub fn cell_idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Perform one full BGK collision + periodic streaming step.
    pub fn step(&mut self) {
        let omega = 1.0 / self.tau;
        let n = self.nx * self.ny * self.nz;
        for idx in 0..n {
            let base = idx * 19;
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for a in 0..19 {
                let fi = self.f[base + a];
                rho += fi;
                mx += fi * C19[a][0] as f64;
                my += fi * C19[a][1] as f64;
                mz += fi * C19[a][2] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            let uz = if rho.abs() > 1e-15 { mz / rho } else { 0.0 };
            let u_sq = ux * ux + uy * uy + uz * uz;
            for a in 0..19 {
                let cx = C19[a][0] as f64;
                let cy = C19[a][1] as f64;
                let cz = C19[a][2] as f64;
                let eu = cx * ux + cy * uy + cz * uz;
                let feq = W19[a]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                self.f[base + a] -= omega * (self.f[base + a] - feq);
            }
        }
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let f_old = self.f.clone();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = (z * ny * nx + y * nx + x) * 19;
                    for a in 0..19 {
                        let sx = (x as isize - C19[a][0] as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - C19[a][1] as isize).rem_euclid(ny as isize) as usize;
                        let sz = (z as isize - C19[a][2] as isize).rem_euclid(nz as isize) as usize;
                        let src = (sz * ny * nx + sy * nx + sx) * 19 + a;
                        self.f[dst + a] = f_old[src];
                    }
                }
            }
        }
    }
    /// Return macroscopic density at cell (x, y, z).
    pub fn density(&self, x: usize, y: usize, z: usize) -> f64 {
        let base = self.cell_idx(x, y, z) * 19;
        (0..19).map(|a| self.f[base + a]).sum()
    }
    /// Return macroscopic velocity \[ux, uy, uz\] at cell (x, y, z).
    pub fn velocity(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let base = self.cell_idx(x, y, z) * 19;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for a in 0..19 {
            let fi = self.f[base + a];
            rho += fi;
            mx += fi * C19[a][0] as f64;
            my += fi * C19[a][1] as f64;
            mz += fi * C19[a][2] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho, mz / rho]
        } else {
            [0.0; 3]
        }
    }
    /// Set all inlet cells (x = 0) to equilibrium with prescribed `(ux, uy, uz)`.
    pub fn set_inlet_velocity(&mut self, ux: f64, uy: f64, uz: f64) {
        let ny = self.ny;
        let nz = self.nz;
        let u_sq = ux * ux + uy * uy + uz * uz;
        for z in 0..nz {
            for y in 0..ny {
                let idx = self.cell_idx(0, y, z);
                for a in 0..19 {
                    let cx = C19[a][0] as f64;
                    let cy = C19[a][1] as f64;
                    let cz = C19[a][2] as f64;
                    let eu = cx * ux + cy * uy + cz * uz;
                    self.f[idx * 19 + a] = W19[a]
                        * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                }
            }
        }
    }
    /// Apply full-way bounce-back on all cells at y = 0 and y = ny - 1.
    pub fn apply_bounce_back(&mut self) {
        for z in 0..self.nz {
            for x in 0..self.nx {
                for &y in &[0_usize, self.ny - 1] {
                    let base = self.cell_idx(x, y, z) * 19;
                    let mut tmp = [0.0_f64; 19];
                    tmp.copy_from_slice(&self.f[base..base + 19]);
                    for a in 0..19 {
                        self.f[base + a] = tmp[OPP19[a]];
                    }
                }
            }
        }
    }
    /// Total mass across all cells (for conservation checks).
    pub fn total_mass(&self) -> f64 {
        self.f.iter().sum()
    }
}
/// 2D LBM grid that uses the TRT (two-relaxation-time) collision operator.
#[allow(dead_code)]
pub struct TrtGrid2D {
    /// Width in cells.
    pub nx: usize,
    /// Height in cells.
    pub ny: usize,
    /// D2Q9 population arrays: `pop[y * nx + x]`.
    pub pop: Vec<[f64; 9]>,
    /// Wall flags: `true` means solid (bounce-back).
    pub wall: Vec<bool>,
    /// Symmetric relaxation frequency (controls viscosity).
    pub omega_sym: f64,
    /// Anti-symmetric relaxation frequency (controls numerical diffusion).
    pub omega_anti: f64,
}
impl TrtGrid2D {
    /// Create a TRT grid with automatic "magic" anti-symmetric rate.
    #[allow(dead_code)]
    pub fn new_magic(nx: usize, ny: usize, omega_sym: f64) -> Self {
        let omega_anti = trt_magic_omega_anti(omega_sym);
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        Self {
            nx,
            ny,
            pop: vec![feq; nx * ny],
            wall: vec![false; nx * ny],
            omega_sym,
            omega_anti,
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
    /// Compute total mass.
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f64 {
        self.pop.iter().flat_map(|n| n.iter()).sum()
    }
    /// Apply TRT collision to all fluid cells.
    #[allow(dead_code)]
    pub fn collide(&mut self) {
        for idx in 0..self.pop.len() {
            if !self.wall[idx] {
                let (rho, ux, uy) = macros_from_d2q9(&self.pop[idx]);
                trt_collision_d2q9(
                    &mut self.pop[idx],
                    rho,
                    ux,
                    uy,
                    self.omega_sym,
                    self.omega_anti,
                );
            }
        }
    }
    /// Apply periodic streaming.
    #[allow(dead_code)]
    pub fn stream(&mut self) {
        stream_d2q9_periodic(&mut self.pop, self.nx, self.ny);
    }
    /// Apply bounce-back to wall cells.
    #[allow(dead_code)]
    pub fn bounce_back(&mut self) {
        for idx in 0..self.pop.len() {
            if self.wall[idx] {
                bounce_back_node_d2q9(&mut self.pop[idx]);
            }
        }
    }
    /// Execute one TRT LBM step.
    #[allow(dead_code)]
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.bounce_back();
    }
    /// Get macroscopic variables at cell `(x, y)`.
    #[allow(dead_code)]
    pub fn macros_at(&self, x: usize, y: usize) -> (f64, f64, f64) {
        macros_from_d2q9(&self.pop[self.idx(x, y)])
    }
}
/// 2D LBM grid for D2Q9 simulations.
///
/// Distribution functions are stored as `f[q][y][x]`.
/// Macroscopic density and velocity are cached and can be recomputed.
#[derive(Debug, Clone)]
pub struct LbmGrid2D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Lattice descriptor.
    pub lattice: Lattice,
    /// Distribution functions: `f[q][y * nx + x]`.
    pub f: Vec<Vec<f64>>,
    /// Macroscopic density field: `rho[y * nx + x]`.
    pub rho: Vec<f64>,
    /// Macroscopic x-velocity field.
    pub ux: Vec<f64>,
    /// Macroscopic y-velocity field.
    pub uy: Vec<f64>,
}
impl LbmGrid2D {
    /// Create a new 2D grid initialized to equilibrium at rest with `rho=1`.
    ///
    /// # Panics
    ///
    /// Panics if `lattice_type` is not `D2Q9`.
    pub fn new(nx: usize, ny: usize, lattice_type: LatticeType) -> Self {
        assert_eq!(lattice_type, LatticeType::D2Q9, "LbmGrid2D requires D2Q9");
        let lattice = Lattice::new(lattice_type);
        let n = nx * ny;
        let q = lattice.q();
        let rho = vec![1.0; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let mut f = Vec::with_capacity(q);
        for i in 0..q {
            let w = lattice.weight(i);
            f.push(vec![w; n]);
        }
        Self {
            nx,
            ny,
            lattice,
            f,
            rho,
            ux,
            uy,
        }
    }
    /// Linear index for cell (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Recompute macroscopic density and velocity from distributions.
    pub fn compute_macroscopic(&mut self) {
        let q = self.lattice.q();
        let n = self.nx * self.ny;
        for k in 0..n {
            let mut rho_k = 0.0;
            let mut ux_k = 0.0;
            let mut uy_k = 0.0;
            for i in 0..q {
                let fi = self.f[i][k];
                rho_k += fi;
                let c = self.lattice.velocity_2d(i);
                ux_k += fi * c[0] as f64;
                uy_k += fi * c[1] as f64;
            }
            self.rho[k] = rho_k;
            if rho_k.abs() > 1e-15 {
                self.ux[k] = ux_k / rho_k;
                self.uy[k] = uy_k / rho_k;
            } else {
                self.ux[k] = 0.0;
                self.uy[k] = 0.0;
            }
        }
    }
    /// Return the density at cell (x, y).
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        self.rho[self.idx(x, y)]
    }
    /// Return the velocity (ux, uy) at cell (x, y).
    pub fn velocity_at(&self, x: usize, y: usize) -> (f64, f64) {
        let k = self.idx(x, y);
        (self.ux[k], self.uy[k])
    }
    /// Set macroscopic fields and reinitialize distributions to equilibrium.
    pub fn set_equilibrium(&mut self, x: usize, y: usize, rho: f64, ux: f64, uy: f64) {
        let k = self.idx(x, y);
        self.rho[k] = rho;
        self.ux[k] = ux;
        self.uy[k] = uy;
        let q = self.lattice.q();
        for i in 0..q {
            let w = self.lattice.weight(i);
            let c = self.lattice.velocity_2d(i);
            self.f[i][k] = equilibrium_2d(w, rho, ux, uy, c[0] as f64, c[1] as f64);
        }
    }
    /// Total density across all cells (for conservation checks).
    pub fn total_density(&self) -> f64 {
        self.rho.iter().sum()
    }
}
/// A 3D LBM grid (D3Q19) with per-node flags, BGK collision and streaming.
///
/// Supports `Fluid`, `Wall` (bounce-back), `Inlet`, `Outlet`, and `Symmetry` nodes.
/// Distribution layout: flat `f[cell_idx * 19 + alpha]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlaggedGrid3D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
    /// Relaxation frequency omega = 1/tau.
    pub omega: f64,
    /// Distribution functions (flat): `f[cell_idx * 19 + alpha]`.
    pub f: Vec<f64>,
    /// Per-node flags.
    pub flags: Vec<NodeFlag>,
}
#[allow(dead_code)]
impl FlaggedGrid3D {
    /// Create a new 3D grid (D3Q19) initialised to equilibrium at rest.
    pub fn new(nx: usize, ny: usize, nz: usize, omega: f64, rho0: f64) -> Self {
        let n = nx * ny * nz;
        let w19 = crate::lattice::D3Q19_WEIGHTS;
        let mut f = vec![0.0; n * 19];
        for idx in 0..n {
            for a in 0..19 {
                f[idx * 19 + a] = w19[a] * rho0;
            }
        }
        Self {
            nx,
            ny,
            nz,
            omega,
            f,
            flags: vec![NodeFlag::Fluid; n],
        }
    }
    /// Linear index for cell (x, y, z).
    #[inline]
    pub fn cell_idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Mark a cell as a solid wall (bounce-back).
    pub fn set_wall(&mut self, x: usize, y: usize, z: usize) {
        let idx = self.cell_idx(x, y, z);
        self.flags[idx] = NodeFlag::Wall;
    }
    /// Set a cell as an inlet with prescribed macroscopic values.
    pub fn set_inlet(&mut self, x: usize, y: usize, z: usize, rho: f64, ux: f64, uy: f64) {
        let idx = self.cell_idx(x, y, z);
        self.flags[idx] = NodeFlag::Inlet {
            rho: rho.to_bits(),
            ux_bits: ux.to_bits(),
            uy_bits: uy.to_bits(),
        };
        let base = idx * 19;
        for a in 0..19 {
            self.f[base + a] = equilibrium_3d(
                W19[a],
                rho,
                ux,
                uy,
                0.0,
                C19[a][0] as f64,
                C19[a][1] as f64,
                C19[a][2] as f64,
            );
        }
    }
    /// Return the macroscopic density at cell (x, y, z).
    pub fn density_at(&self, x: usize, y: usize, z: usize) -> f64 {
        let base = self.cell_idx(x, y, z) * 19;
        (0..19).map(|a| self.f[base + a]).sum()
    }
    /// Return the macroscopic velocity at cell (x, y, z).
    pub fn velocity_at(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let base = self.cell_idx(x, y, z) * 19;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for a in 0..19 {
            let fi = self.f[base + a];
            rho += fi;
            mx += fi * C19[a][0] as f64;
            my += fi * C19[a][1] as f64;
            mz += fi * C19[a][2] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho, mz / rho]
        } else {
            [0.0; 3]
        }
    }
    /// Total mass across all cells.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().sum()
    }
    /// Initialize all fluid cells to equilibrium at uniform (rho0, ux0, uy0, uz0).
    pub fn initialize_uniform(&mut self, rho0: f64, ux0: f64, uy0: f64, uz0: f64) {
        let n = self.nx * self.ny * self.nz;
        for idx in 0..n {
            if self.flags[idx] == NodeFlag::Fluid || self.flags[idx] == NodeFlag::Symmetry {
                let base = idx * 19;
                for a in 0..19 {
                    self.f[base + a] = equilibrium_3d(
                        W19[a],
                        rho0,
                        ux0,
                        uy0,
                        uz0,
                        C19[a][0] as f64,
                        C19[a][1] as f64,
                        C19[a][2] as f64,
                    );
                }
            }
        }
    }
    /// Perform one full BGK collide → stream → boundary-condition step.
    ///
    /// - `Fluid` / `Symmetry`: standard BGK collision then pull-scheme streaming.
    /// - `Wall`: no-slip full-way bounce-back (populations reflected in place,
    ///   then overwritten by streaming so the net effect is opposite-direction
    ///   reflection at solid nodes).
    /// - `Inlet`: BGK collision during collision pass; overwritten with
    ///   equilibrium at prescribed density/velocity after streaming.
    /// - `Outlet`: BGK collision during collision pass; overwritten with
    ///   zero-gradient (upstream copy) after streaming.
    pub fn step(&mut self) {
        let omega = self.omega;
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let n = nx * ny * nz;

        // --- Collision pass ---
        // For fluid/symmetry/inlet/outlet: BGK.
        // For wall: in-place full-way bounce-back (swap opposite populations).
        let mut f_coll = self.f.clone();
        for idx in 0..n {
            let base = idx * 19;
            match self.flags[idx] {
                NodeFlag::Fluid
                | NodeFlag::Symmetry
                | NodeFlag::Inlet { .. }
                | NodeFlag::Outlet => {
                    let mut rho = 0.0_f64;
                    let mut mx = 0.0_f64;
                    let mut my = 0.0_f64;
                    let mut mz = 0.0_f64;
                    for a in 0..19 {
                        let fi = self.f[base + a];
                        rho += fi;
                        mx += fi * C19[a][0] as f64;
                        my += fi * C19[a][1] as f64;
                        mz += fi * C19[a][2] as f64;
                    }
                    let (ux, uy, uz) = if rho.abs() > 1e-15 {
                        (mx / rho, my / rho, mz / rho)
                    } else {
                        (0.0, 0.0, 0.0)
                    };
                    for a in 0..19 {
                        let feq = equilibrium_3d(
                            W19[a],
                            rho,
                            ux,
                            uy,
                            uz,
                            C19[a][0] as f64,
                            C19[a][1] as f64,
                            C19[a][2] as f64,
                        );
                        f_coll[base + a] = self.f[base + a] - omega * (self.f[base + a] - feq);
                    }
                }
                NodeFlag::Wall => {
                    // Full-way bounce-back: reflect all populations.
                    let mut tmp = [0.0_f64; 19];
                    tmp.copy_from_slice(&self.f[base..base + 19]);
                    for a in 0..19 {
                        f_coll[base + a] = tmp[OPP19[a]];
                    }
                }
            }
        }

        // --- Streaming pass (pull scheme) ---
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = (z * ny * nx + y * nx + x) * 19;
                    for a in 0..19 {
                        let sx = (x as isize - C19[a][0] as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - C19[a][1] as isize).rem_euclid(ny as isize) as usize;
                        let sz = (z as isize - C19[a][2] as isize).rem_euclid(nz as isize) as usize;
                        let src = (sz * ny * nx + sy * nx + sx) * 19 + a;
                        self.f[dst + a] = f_coll[src];
                    }
                }
            }
        }

        // --- Boundary condition enforcement ---
        for idx in 0..n {
            match self.flags[idx] {
                NodeFlag::Fluid | NodeFlag::Symmetry | NodeFlag::Wall => {}
                NodeFlag::Inlet {
                    rho,
                    ux_bits,
                    uy_bits,
                } => {
                    let r = f64::from_bits(rho);
                    let ux = f64::from_bits(ux_bits);
                    let uy = f64::from_bits(uy_bits);
                    let base = idx * 19;
                    for a in 0..19 {
                        self.f[base + a] = equilibrium_3d(
                            W19[a],
                            r,
                            ux,
                            uy,
                            0.0,
                            C19[a][0] as f64,
                            C19[a][1] as f64,
                            C19[a][2] as f64,
                        );
                    }
                }
                NodeFlag::Outlet => {
                    // Zero-gradient: copy from upstream neighbor (idx - nx, same
                    // convention as FullGrid3D::apply_bcs).
                    let upstream_nx = self.nx;
                    if idx >= upstream_nx {
                        let src_idx = idx - upstream_nx;
                        for a in 0..19 {
                            self.f[idx * 19 + a] = self.f[src_idx * 19 + a];
                        }
                    }
                }
            }
        }
    }
}
/// A full 3D LBM grid with BGK collision, streaming, and per-node BCs.
///
/// Distribution layout: flat `f[cell_idx * 19 + alpha]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FullGrid3D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
    /// Relaxation frequency omega = 1/tau.
    pub omega: f64,
    /// Distribution functions: `f[cell_idx * 19 + alpha]`.
    pub f: Vec<f64>,
    /// Per-node flags.
    pub flags: Vec<NodeFlag>,
}
#[allow(dead_code)]
impl FullGrid3D {
    /// Create a new FullGrid3D initialized to equilibrium at rest with `rho0`.
    pub fn new(nx: usize, ny: usize, nz: usize, omega: f64, rho0: f64) -> Self {
        let n = nx * ny * nz;
        let mut f = vec![0.0_f64; n * 19];
        for idx in 0..n {
            for alpha in 0..19 {
                f[idx * 19 + alpha] = W19[alpha] * rho0;
            }
        }
        Self {
            nx,
            ny,
            nz,
            omega,
            f,
            flags: vec![NodeFlag::Fluid; n],
        }
    }
    /// Linear index for cell (x, y, z).
    #[inline]
    pub fn cell_idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Perform one full time step: collide, stream, apply boundary conditions.
    pub fn step(&mut self) {
        self.collide_bgk();
        self.stream_periodic();
        self.apply_bcs();
    }
    /// BGK collision on fluid nodes.
    fn collide_bgk(&mut self) {
        let n = self.nx * self.ny * self.nz;
        let omega = self.omega;
        for idx in 0..n {
            if self.flags[idx] != NodeFlag::Fluid {
                continue;
            }
            let base = idx * 19;
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            let mut mz = 0.0_f64;
            for a in 0..19 {
                let fi = self.f[base + a];
                rho += fi;
                mx += fi * C19[a][0] as f64;
                my += fi * C19[a][1] as f64;
                mz += fi * C19[a][2] as f64;
            }
            let (ux, uy, uz) = if rho.abs() > 1e-15 {
                (mx / rho, my / rho, mz / rho)
            } else {
                (0.0, 0.0, 0.0)
            };
            for a in 0..19 {
                let feq = equilibrium_3d(
                    W19[a],
                    rho,
                    ux,
                    uy,
                    uz,
                    C19[a][0] as f64,
                    C19[a][1] as f64,
                    C19[a][2] as f64,
                );
                self.f[base + a] -= omega * (self.f[base + a] - feq);
            }
        }
    }
    /// Periodic pull-scheme streaming.
    fn stream_periodic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let f_old = self.f.clone();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = (z * ny * nx + y * nx + x) * 19;
                    for a in 0..19 {
                        let sx = (x as isize - C19[a][0] as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - C19[a][1] as isize).rem_euclid(ny as isize) as usize;
                        let sz = (z as isize - C19[a][2] as isize).rem_euclid(nz as isize) as usize;
                        let src = (sz * ny * nx + sy * nx + sx) * 19 + a;
                        self.f[dst + a] = f_old[src];
                    }
                }
            }
        }
    }
    /// Apply boundary conditions for all flagged nodes.
    fn apply_bcs(&mut self) {
        let n = self.nx * self.ny * self.nz;
        for idx in 0..n {
            match self.flags[idx] {
                NodeFlag::Fluid | NodeFlag::Symmetry => {}
                NodeFlag::Wall => {
                    let base = idx * 19;
                    let mut tmp = [0.0_f64; 19];
                    tmp.copy_from_slice(&self.f[base..base + 19]);
                    for a in 0..19 {
                        self.f[base + a] = tmp[OPP19[a]];
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
                    let base = idx * 19;
                    for a in 0..19 {
                        self.f[base + a] = equilibrium_3d(
                            W19[a],
                            r,
                            ux,
                            uy,
                            0.0,
                            C19[a][0] as f64,
                            C19[a][1] as f64,
                            C19[a][2] as f64,
                        );
                    }
                }
                NodeFlag::Outlet => {
                    let nx = self.nx;
                    let ny = self.ny;
                    let total = nx * ny;
                    if idx >= total {
                        let src_idx = idx - nx;
                        for a in 0..19 {
                            self.f[idx * 19 + a] = self.f[src_idx * 19 + a];
                        }
                    }
                }
            }
        }
    }
    /// Mark a cell as a solid wall.
    pub fn set_wall(&mut self, x: usize, y: usize, z: usize) {
        let idx = self.cell_idx(x, y, z);
        self.flags[idx] = NodeFlag::Wall;
    }
    /// Set a cell as an inlet with prescribed macroscopic values.
    pub fn set_inlet(&mut self, x: usize, y: usize, z: usize, rho: f64, ux: f64, uy: f64) {
        let idx = self.cell_idx(x, y, z);
        self.flags[idx] = NodeFlag::Inlet {
            rho: rho.to_bits(),
            ux_bits: ux.to_bits(),
            uy_bits: uy.to_bits(),
        };
        let base = idx * 19;
        for a in 0..19 {
            self.f[base + a] = equilibrium_3d(
                W19[a],
                rho,
                ux,
                uy,
                0.0,
                C19[a][0] as f64,
                C19[a][1] as f64,
                C19[a][2] as f64,
            );
        }
    }
    /// Return the macroscopic density at cell (x, y, z).
    pub fn density_at(&self, x: usize, y: usize, z: usize) -> f64 {
        let base = self.cell_idx(x, y, z) * 19;
        (0..19).map(|a| self.f[base + a]).sum()
    }
    /// Return the macroscopic velocity at cell (x, y, z).
    pub fn velocity_at(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let base = self.cell_idx(x, y, z) * 19;
        let mut rho = 0.0_f64;
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut mz = 0.0;
        for a in 0..19 {
            let fi = self.f[base + a];
            rho += fi;
            mx += fi * C19[a][0] as f64;
            my += fi * C19[a][1] as f64;
            mz += fi * C19[a][2] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho, mz / rho]
        } else {
            [0.0; 3]
        }
    }
    /// Total mass across all cells (mass conservation check).
    pub fn total_mass(&self) -> f64 {
        let n = self.nx * self.ny * self.nz;
        let mut total = 0.0_f64;
        for idx in 0..n {
            let base = idx * 19;
            for a in 0..19 {
                total += self.f[base + a];
            }
        }
        total
    }
    /// Compute maximum velocity magnitude across all cells.
    pub fn max_velocity_magnitude(&self) -> f64 {
        let n = self.nx * self.ny * self.nz;
        let mut max_v = 0.0_f64;
        for idx in 0..n {
            let v = {
                let base = idx * 19;
                let mut rho = 0.0_f64;
                let mut mx = 0.0;
                let mut my = 0.0;
                let mut mz = 0.0;
                for a in 0..19 {
                    let fi = self.f[base + a];
                    rho += fi;
                    mx += fi * C19[a][0] as f64;
                    my += fi * C19[a][1] as f64;
                    mz += fi * C19[a][2] as f64;
                }
                if rho.abs() > 1e-15 {
                    let ux = mx / rho;
                    let uy = my / rho;
                    let uz = mz / rho;
                    (ux * ux + uy * uy + uz * uz).sqrt()
                } else {
                    0.0
                }
            };
            if v > max_v {
                max_v = v;
            }
        }
        max_v
    }
    /// Save checkpoint: return a clone of the distribution function array.
    pub fn checkpoint_save(&self) -> Vec<f64> {
        self.f.clone()
    }
    /// Load checkpoint: restore distribution functions from saved data.
    pub fn checkpoint_load(&mut self, data: &[f64]) {
        assert_eq!(data.len(), self.f.len(), "checkpoint size mismatch");
        self.f.copy_from_slice(data);
    }
    /// Check mass conservation between current state and a saved checkpoint.
    ///
    /// Returns the absolute difference in total mass.
    pub fn mass_conservation_error(&self, checkpoint: &[f64]) -> f64 {
        let saved_mass: f64 = checkpoint.iter().sum();
        (self.total_mass() - saved_mass).abs()
    }
    /// Initialize all cells to equilibrium at uniform (rho0, ux0, uy0, uz0).
    pub fn initialize_uniform(&mut self, rho0: f64, ux0: f64, uy0: f64, uz0: f64) {
        let n = self.nx * self.ny * self.nz;
        for idx in 0..n {
            let base = idx * 19;
            for a in 0..19 {
                self.f[base + a] = equilibrium_3d(
                    W19[a],
                    rho0,
                    ux0,
                    uy0,
                    uz0,
                    C19[a][0] as f64,
                    C19[a][1] as f64,
                    C19[a][2] as f64,
                );
            }
        }
    }
}
/// Typed classification of a 2-D grid node for boundary-condition handling.
///
/// This enum is an alternative to the bit-packed `NodeFlag` above and provides
/// a cleaner API for setting up inlet/outlet/solid boundary conditions.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryNodeType {
    /// Normal fluid node: BGK collision + streaming.
    Fluid,
    /// Solid wall: full-way bounce-back; no fluid flow.
    Solid,
    /// Inlet: Zou-He velocity BC with prescribed (ux, uy) and density rho.
    Inlet {
        /// Prescribed x-velocity.
        ux: f64,
        /// Prescribed y-velocity.
        uy: f64,
        /// Prescribed density.
        rho: f64,
    },
    /// Outlet: Zou-He pressure BC (zero-gradient extrapolation).
    Outlet,
}
/// A 3D LBM grid storing D3Q19 populations per cell.
#[allow(dead_code)]
pub struct CellularGrid3D {
    /// Width in cells (x-direction).
    pub nx: usize,
    /// Height in cells (y-direction).
    pub ny: usize,
    /// Depth in cells (z-direction).
    pub nz: usize,
    /// D3Q19 distribution functions: `pop[z*ny*nx + y*nx + x]`.
    pub pop: Vec<[f64; 19]>,
    /// Wall flags: `true` means solid (bounce-back).
    pub wall: Vec<bool>,
    /// BGK relaxation frequency ω = 1/τ.
    pub omega: f64,
}
impl CellularGrid3D {
    /// Create a 3D grid initialised to rest equilibrium (ρ=1, u=0).
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, nz: usize, omega: f64) -> Self {
        let n = nx * ny * nz;
        let feq = equilibrium_d3q19(1.0, 0.0, 0.0, 0.0);
        Self {
            nx,
            ny,
            nz,
            pop: vec![feq; n],
            wall: vec![false; n],
            omega,
        }
    }
    /// Linear index for cell `(x, y, z)`.
    #[inline]
    #[allow(dead_code)]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.nx * self.ny + y * self.nx + x
    }
    /// Mark a cell as solid wall.
    #[allow(dead_code)]
    pub fn set_wall(&mut self, x: usize, y: usize, z: usize) {
        let i = self.idx(x, y, z);
        self.wall[i] = true;
    }
    /// Compute total mass (sum of all populations).
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f64 {
        self.pop.iter().flat_map(|n| n.iter()).sum()
    }
    /// Get macroscopic variables at cell `(x, y, z)`.
    #[allow(dead_code)]
    pub fn macros_at(&self, x: usize, y: usize, z: usize) -> (f64, f64, f64, f64) {
        macros_from_d3q19(&self.pop[self.idx(x, y, z)])
    }
    /// Apply BGK collision to all fluid cells.
    #[allow(dead_code)]
    pub fn collide(&mut self) {
        for idx in 0..self.pop.len() {
            if !self.wall[idx] {
                let (rho, ux, uy, uz) = macros_from_d3q19(&self.pop[idx]);
                bgk_d3q19(&mut self.pop[idx], rho, ux, uy, uz, self.omega);
            }
        }
    }
    /// Apply periodic streaming.
    #[allow(dead_code)]
    pub fn stream(&mut self) {
        stream_d3q19_periodic(&mut self.pop, self.nx, self.ny, self.nz);
    }
    /// Apply half-way bounce-back to all wall cells.
    #[allow(dead_code)]
    pub fn bounce_back(&mut self) {
        for idx in 0..self.pop.len() {
            if self.wall[idx] {
                bounce_back_node_d3q19(&mut self.pop[idx]);
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
    /// Return maximum velocity magnitude over all fluid cells.
    #[allow(dead_code)]
    pub fn max_speed(&self) -> f64 {
        self.pop
            .iter()
            .enumerate()
            .filter(|(i, _)| !self.wall[*i])
            .map(|(_, node)| {
                let (_, ux, uy, uz) = macros_from_d3q19(node);
                (ux * ux + uy * uy + uz * uz).sqrt()
            })
            .fold(0.0f64, f64::max)
    }
}
/// 2D LBM grid (D2Q9) with flat distribution storage and tau-based collision.
///
/// Distribution layout: `f[y * nx * 9 + x * 9 + alpha]`.
/// The relaxation time `tau` relates to kinematic viscosity via
/// `nu = (tau - 0.5) * cs²`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlatLbmGrid2D {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Distribution functions: flat, length `nx * ny * 9`.
    pub f: Vec<f64>,
    /// Relaxation time tau (= 1/omega).
    pub tau: f64,
}
#[allow(dead_code)]
impl FlatLbmGrid2D {
    /// Create a new grid initialized to equilibrium at rest with `rho=1`.
    pub fn new(nx: usize, ny: usize, tau: f64) -> Self {
        let n = nx * ny;
        let mut f = vec![0.0_f64; n * 9];
        for idx in 0..n {
            for alpha in 0..9 {
                f[idx * 9 + alpha] = W9[alpha];
            }
        }
        Self { nx, ny, f, tau }
    }
    /// Linear cell index for (x, y).
    #[inline]
    pub fn cell_idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Perform one full BGK collision + periodic streaming step.
    pub fn step(&mut self) {
        let omega = 1.0 / self.tau;
        let n = self.nx * self.ny;
        for idx in 0..n {
            let base = idx * 9;
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            for a in 0..9 {
                let fi = self.f[base + a];
                rho += fi;
                mx += fi * C9[a][0] as f64;
                my += fi * C9[a][1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            let u_sq = ux * ux + uy * uy;
            for a in 0..9 {
                let cx = C9[a][0] as f64;
                let cy = C9[a][1] as f64;
                let eu = cx * ux + cy * uy;
                let feq = W9[a]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                self.f[base + a] -= omega * (self.f[base + a] - feq);
            }
        }
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let dst = (y * nx + x) * 9;
                for a in 0..9 {
                    let sx = (x as isize - C9[a][0] as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - C9[a][1] as isize).rem_euclid(ny as isize) as usize;
                    self.f[dst + a] = f_old[(sy * nx + sx) * 9 + a];
                }
            }
        }
    }
    /// Return macroscopic density at cell (x, y).
    pub fn density(&self, x: usize, y: usize) -> f64 {
        let base = self.cell_idx(x, y) * 9;
        (0..9).map(|a| self.f[base + a]).sum()
    }
    /// Return macroscopic velocity \[ux, uy\] at cell (x, y).
    pub fn velocity(&self, x: usize, y: usize) -> [f64; 2] {
        let base = self.cell_idx(x, y) * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for a in 0..9 {
            let fi = self.f[base + a];
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
    /// Set all inlet cells (x = 0) to equilibrium with prescribed `(ux, uy)`.
    ///
    /// Density is fixed to `rho = 1.0` at the inlet.
    pub fn set_inlet_velocity(&mut self, ux: f64, uy: f64) {
        for y in 0..self.ny {
            let idx = self.cell_idx(0, y);
            let u_sq = ux * ux + uy * uy;
            for a in 0..9 {
                let cx = C9[a][0] as f64;
                let cy = C9[a][1] as f64;
                let eu = cx * ux + cy * uy;
                self.f[idx * 9 + a] =
                    W9[a] * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
            }
        }
    }
    /// Apply full-way bounce-back on all cells at y = 0 and y = ny - 1.
    ///
    /// This imposes no-slip walls on the top and bottom boundaries.
    pub fn apply_bounce_back(&mut self) {
        for x in 0..self.nx {
            for &y in &[0_usize, self.ny - 1] {
                let base = self.cell_idx(x, y) * 9;
                let mut tmp = [0.0_f64; 9];
                tmp.copy_from_slice(&self.f[base..base + 9]);
                for a in 0..9 {
                    self.f[base + a] = tmp[OPP9[a]];
                }
            }
        }
    }
    /// Total mass across all cells (for conservation checks).
    pub fn total_mass(&self) -> f64 {
        self.f.iter().sum()
    }
}
/// Classification of a grid node for boundary-condition handling.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeFlag {
    /// Normal fluid node.
    Fluid,
    /// Solid wall (bounce-back).
    Wall,
    /// Inlet with prescribed density and velocity.
    Inlet {
        /// Prescribed density at this inlet node.
        rho: u64,
        /// Prescribed x-velocity (stored as bits of f64).
        ux_bits: u64,
        /// Prescribed y-velocity (stored as bits of f64).
        uy_bits: u64,
    },
    /// Outlet (extrapolation / zero-gradient).
    Outlet,
    /// Symmetry plane.
    Symmetry,
}
