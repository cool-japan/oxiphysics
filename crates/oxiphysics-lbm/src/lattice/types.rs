//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    CS2, D2Q9_OPPOSITES, D2Q9_VELOCITIES, D2Q9_WEIGHTS, D3Q19_OPPOSITES, D3Q19_VELOCITIES,
    D3Q19_WEIGHTS, D3Q27_OPPOSITES, D3Q27_VELOCITIES, D3Q27_WEIGHTS, MRT_M, bgk_collision,
    compute_rho, compute_velocity,
};

/// Guo body-force scheme for LBM.
///
/// Adds a body force `[gx, gy]` (acceleration) to the LBM collision step.
/// The modified distribution: f_i^* = f_i - omega * (f_i - feq) + (1 - omega/2) * F_i * dt
/// where `F_i = w_i * (c_i - u)/cs^2 * g`.
pub struct GuoBodyForce {
    /// Body force acceleration \[gx, gy\] (lattice units/step^2).
    pub g: [f64; 2],
}
impl GuoBodyForce {
    /// Create a new Guo body force scheme.
    pub fn new(gx: f64, gy: f64) -> Self {
        Self { g: [gx, gy] }
    }
    /// Compute the Guo forcing term F_i for a given velocity (ux, uy).
    pub fn forcing_term(&self, ux: f64, uy: f64, alpha: usize) -> f64 {
        let cx = D2Q9_VELOCITIES[alpha][0] as f64;
        let cy = D2Q9_VELOCITIES[alpha][1] as f64;
        let w = D2Q9_WEIGHTS[alpha];
        let cdotu = cx * ux + cy * uy;
        let cidotg = cx * self.g[0] + cy * self.g[1];
        let udotg = ux * self.g[0] + uy * self.g[1];
        w * (cidotg / CS2 + cdotu * cidotg / (CS2 * CS2) - udotg / CS2)
    }
    /// Apply BGK collision with Guo body force to a D2Q9 grid.
    pub fn collide_with_force(&self, f: &mut [f64], n: usize, omega: f64) {
        for idx in 0..n {
            let base = idx * 9;
            let mut rho = 0.0_f64;
            let mut mx = 0.0;
            let mut my = 0.0;
            for a in 0..9 {
                let fi = f[base + a];
                rho += fi;
                mx += fi * D2Q9_VELOCITIES[a][0] as f64;
                my += fi * D2Q9_VELOCITIES[a][1] as f64;
            }
            if rho < 1e-15 {
                continue;
            }
            let ux = mx / rho + 0.5 * self.g[0];
            let uy = my / rho + 0.5 * self.g[1];
            for a in 0..9 {
                let feq = D2Q9Grid::feq(rho, ux, uy, a);
                let fi = self.forcing_term(ux, uy, a);
                f[base + a] = f[base + a] - omega * (f[base + a] - feq) + (1.0 - 0.5 * omega) * fi;
            }
        }
    }
}
/// Typed D2Q9 velocity-set struct.
///
/// Provides the weights and velocity vectors as associated constants,
/// and convenience methods for equilibrium, BGK collision, and streaming.
pub struct D2Q9;
impl D2Q9 {
    /// D2Q9 weights.
    pub const W: [f64; 9] = D2Q9_WEIGHTS;
    /// D2Q9 velocity vectors `[cx, cy]`.
    pub const C: [[i32; 2]; 9] = D2Q9_VELOCITIES;
    /// Compute the D2Q9 equilibrium distribution for density `rho` and velocity (ux, uy).
    ///
    /// Returns `[f_eq_0, ..., f_eq_8]`.
    pub fn equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let mut feq = [0.0_f64; 9];
        let u_sq = ux * ux + uy * uy;
        for (i, fi) in feq.iter_mut().enumerate() {
            let cx = Self::C[i][0] as f64;
            let cy = Self::C[i][1] as f64;
            let eu = cx * ux + cy * uy;
            *fi = Self::W[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        feq
    }
    /// Pull-scheme periodic streaming on a flat `Vec<[f64; 9]>` of `nx * ny` cells.
    pub fn stream(f: &mut [[f64; 9]], nx: usize, ny: usize) {
        let f_old = f.to_vec();
        for y in 0..ny {
            for x in 0..nx {
                let dst = y * nx + x;
                for a in 0..9 {
                    let sx = (x as isize - Self::C[a][0] as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - Self::C[a][1] as isize).rem_euclid(ny as isize) as usize;
                    f[dst][a] = f_old[sy * nx + sx][a];
                }
            }
        }
    }
}
/// Typed D3Q27 velocity-set struct.
///
/// Provides the 27-velocity set for 3D simulations with higher accuracy
/// than D3Q19 for flows with large velocity gradients.
pub struct D3Q27;
impl D3Q27 {
    /// D3Q27 weights.
    pub const W: [f64; 27] = D3Q27_WEIGHTS;
    /// D3Q27 velocity vectors `[cx, cy, cz]`.
    pub const C: [[i32; 3]; 27] = D3Q27_VELOCITIES;
    /// Compute the D3Q27 equilibrium distribution for density `rho` and velocity `u`.
    ///
    /// Returns `[f_eq_0, ..., f_eq_26]`.
    pub fn equilibrium(rho: f64, u: [f64; 3]) -> [f64; 27] {
        let mut feq = [0.0_f64; 27];
        let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        for (i, fi) in feq.iter_mut().enumerate() {
            let cx = Self::C[i][0] as f64;
            let cy = Self::C[i][1] as f64;
            let cz = Self::C[i][2] as f64;
            let eu = cx * u[0] + cy * u[1] + cz * u[2];
            *fi = Self::W[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        feq
    }
    /// Compute macroscopic density from a D3Q27 distribution.
    pub fn compute_rho(f: &[f64; 27]) -> f64 {
        f.iter().sum()
    }
    /// Compute macroscopic velocity from a D3Q27 distribution.
    pub fn compute_velocity(f: &[f64; 27], rho: f64) -> [f64; 3] {
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for (i, &fi) in f.iter().enumerate() {
            mx += fi * Self::C[i][0] as f64;
            my += fi * Self::C[i][1] as f64;
            mz += fi * Self::C[i][2] as f64;
        }
        if rho.abs() > 1e-15 {
            [mx / rho, my / rho, mz / rho]
        } else {
            [0.0, 0.0, 0.0]
        }
    }
    /// Pull-scheme periodic streaming on a flat `Vec<[f64; 27]>`.
    pub fn stream(f: &mut [[f64; 27]], nx: usize, ny: usize, nz: usize) {
        let f_old = f.to_vec();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let dst = z * ny * nx + y * nx + x;
                    for a in 0..27 {
                        let sx =
                            (x as isize - Self::C[a][0] as isize).rem_euclid(nx as isize) as usize;
                        let sy =
                            (y as isize - Self::C[a][1] as isize).rem_euclid(ny as isize) as usize;
                        let sz =
                            (z as isize - Self::C[a][2] as isize).rem_euclid(nz as isize) as usize;
                        f[dst][a] = f_old[sz * ny * nx + sy * nx + sx][a];
                    }
                }
            }
        }
    }
    /// Apply full bounce-back (no-slip wall) to a single D3Q27 node.
    pub fn bounce_back(f: &mut [f64; 27]) {
        let mut tmp = [0.0_f64; 27];
        tmp.copy_from_slice(f);
        for a in 0..27 {
            f[D3Q27_OPPOSITES[a]] = tmp[a];
        }
    }
    /// BGK collision in-place on a single D3Q27 node.
    pub fn bgk(f: &mut [f64; 27], rho: f64, u: [f64; 3], omega: f64) {
        let feq = Self::equilibrium(rho, u);
        for i in 0..27 {
            f[i] -= omega * (f[i] - feq[i]);
        }
    }
}
/// TRT (Two-Relaxation-Time) collision operator for D2Q9.
///
/// Uses two relaxation times: omega_plus (symmetric) and omega_minus (anti-symmetric).
/// TRT is known to have better wall-distance independence than BGK.
#[derive(Debug, Clone)]
pub struct TrtCollision {
    /// Relaxation rate for symmetric part (related to viscosity).
    pub omega_plus: f64,
    /// Relaxation rate for anti-symmetric part.
    pub omega_minus: f64,
}
impl TrtCollision {
    /// Create a TRT collision operator.
    ///
    /// The magic parameter `lambda_trt = (omega_plus * omega_minus) / (omega_plus + omega_minus)`.
    pub fn new(omega_plus: f64, omega_minus: f64) -> Self {
        Self {
            omega_plus,
            omega_minus,
        }
    }
    /// Create TRT with the "magic" relaxation parameter for optimal accuracy.
    ///
    /// Uses `lambda = 3/16` (best for Poiseuille flow).
    pub fn with_magic(omega_plus: f64) -> Self {
        const MAGIC: f64 = 3.0 / 16.0;
        let omega_minus = 1.0 / (MAGIC / (1.0 / omega_plus - 0.5) + 0.5);
        Self {
            omega_plus,
            omega_minus,
        }
    }
    /// Apply TRT collision to all cells of a D2Q9 grid.
    ///
    /// Modifies `f` in-place. `n = nx * ny`.
    pub fn collide(&self, f: &mut [f64], n: usize) {
        for idx in 0..n {
            let base = idx * 9;
            let mut rho = 0.0_f64;
            let mut mx = 0.0;
            let mut my = 0.0;
            for a in 0..9 {
                let fi = f[base + a];
                rho += fi;
                mx += fi * D2Q9_VELOCITIES[a][0] as f64;
                my += fi * D2Q9_VELOCITIES[a][1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            for a in 0..9 {
                let opp = D2Q9_OPPOSITES[a];
                let feq_a = D2Q9Grid::feq(rho, ux, uy, a);
                let feq_opp = D2Q9Grid::feq(rho, ux, uy, opp);
                let fa = f[base + a];
                let fopp = f[base + opp];
                let sym = (fa + fopp) * 0.5;
                let sym_eq = (feq_a + feq_opp) * 0.5;
                let anti = (fa - fopp) * 0.5;
                let anti_eq = (feq_a - feq_opp) * 0.5;
                f[base + a] =
                    fa - self.omega_plus * (sym - sym_eq) - self.omega_minus * (anti - anti_eq);
            }
        }
    }
    /// Kinematic viscosity from omega_plus: nu = cs^2 * (1/omega_plus - 0.5).
    pub fn viscosity(&self) -> f64 {
        CS2 * (1.0 / self.omega_plus - 0.5)
    }
}
/// MRT (Multiple-Relaxation-Time) collision operator for D2Q9.
///
/// Allows independent control of all relaxation rates.
#[derive(Debug, Clone)]
pub struct MrtCollision {
    /// Diagonal relaxation matrix S: 9 relaxation rates for the 9 moments.
    pub s: [f64; 9],
}
impl MrtCollision {
    /// Create an MRT collision operator with given relaxation rates.
    ///
    /// `s[1]` and `s[2]` control bulk viscosity; `s[7]` and `s[8]` control shear viscosity.
    pub fn new(s: [f64; 9]) -> Self {
        Self { s }
    }
    /// Create MRT with all non-conserved modes at `omega`.
    pub fn uniform(omega: f64) -> Self {
        Self {
            s: [0.0, 1.6, 1.2, 0.0, 1.2, 0.0, 1.2, omega, omega],
        }
    }
    /// Transform distribution functions f to moment space: m = M * f.
    pub fn to_moment_space(f: &[f64; 9]) -> [f64; 9] {
        let mut m = [0.0_f64; 9];
        for i in 0..9 {
            for j in 0..9 {
                m[i] += MRT_M[i][j] * f[j];
            }
        }
        m
    }
    /// Transform moments m back to distribution space: f = M^{-1} * m.
    pub fn from_moment_space(m: &[f64; 9]) -> [f64; 9] {
        let norms: [f64; 9] = [9.0, 36.0, 36.0, 6.0, 12.0, 6.0, 12.0, 4.0, 4.0];
        let mut f = [0.0_f64; 9];
        for j in 0..9 {
            for i in 0..9 {
                f[j] += MRT_M[i][j] * m[i] / norms[i];
            }
        }
        f
    }
    /// Apply MRT collision to one cell with given rho, ux, uy.
    pub fn collide_cell(&self, f_in: &[f64; 9], rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let m = Self::to_moment_space(f_in);
        let m_eq = [
            rho,
            -2.0 * rho + 3.0 * rho * (ux * ux + uy * uy),
            rho - 3.0 * rho * (ux * ux + uy * uy),
            rho * ux,
            -rho * ux,
            rho * uy,
            -rho * uy,
            rho * (ux * ux - uy * uy),
            rho * ux * uy,
        ];
        let mut m_out = [0.0_f64; 9];
        for i in 0..9 {
            m_out[i] = m[i] - self.s[i] * (m[i] - m_eq[i]);
        }
        Self::from_moment_space(&m_out)
    }
    /// Apply MRT collision to all cells in a D2Q9 grid.
    pub fn collide(&self, f: &mut [f64], n: usize) {
        for idx in 0..n {
            let base = idx * 9;
            let mut rho = 0.0_f64;
            let mut mx = 0.0;
            let mut my = 0.0;
            let mut f9 = [0.0_f64; 9];
            for a in 0..9 {
                f9[a] = f[base + a];
                rho += f9[a];
                mx += f9[a] * D2Q9_VELOCITIES[a][0] as f64;
                my += f9[a] * D2Q9_VELOCITIES[a][1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            let f_out = self.collide_cell(&f9, rho, ux, uy);
            f[base..base + 9].copy_from_slice(&f_out);
        }
    }
}
/// A lattice descriptor providing weights, velocities, and opposite indices.
#[derive(Debug, Clone)]
pub struct Lattice {
    /// The lattice type.
    pub(super) lattice_type: LatticeType,
    /// Weights (length = q).
    pub(super) weights: Vec<f64>,
    /// 2D velocity vectors (only populated for D2Q9).
    pub(super) velocities_2d: Vec<[i32; 2]>,
    /// 3D velocity vectors (only populated for D3Q19 / D3Q27).
    pub(super) velocities_3d: Vec<[i32; 3]>,
    /// Opposite direction index for each direction.
    pub(super) opposites: Vec<usize>,
}
impl Lattice {
    /// Create a new `Lattice` from the given type.
    pub fn new(lattice_type: LatticeType) -> Self {
        match lattice_type {
            LatticeType::D2Q9 => Self {
                lattice_type,
                weights: D2Q9_WEIGHTS.to_vec(),
                velocities_2d: D2Q9_VELOCITIES.to_vec(),
                velocities_3d: Vec::new(),
                opposites: D2Q9_OPPOSITES.to_vec(),
            },
            LatticeType::D3Q19 => Self {
                lattice_type,
                weights: D3Q19_WEIGHTS.to_vec(),
                velocities_2d: Vec::new(),
                velocities_3d: D3Q19_VELOCITIES.to_vec(),
                opposites: D3Q19_OPPOSITES.to_vec(),
            },
            LatticeType::D3Q27 => Self {
                lattice_type,
                weights: D3Q27_WEIGHTS.to_vec(),
                velocities_2d: Vec::new(),
                velocities_3d: D3Q27_VELOCITIES.to_vec(),
                opposites: D3Q27_OPPOSITES.to_vec(),
            },
        }
    }
    /// Return the lattice type.
    pub fn lattice_type(&self) -> LatticeType {
        self.lattice_type
    }
    /// Number of discrete velocities.
    pub fn q(&self) -> usize {
        self.lattice_type.q()
    }
    /// Number of spatial dimensions.
    pub fn dim(&self) -> usize {
        self.lattice_type.dim()
    }
    /// Weight for direction `i`.
    pub fn weight(&self, i: usize) -> f64 {
        self.weights[i]
    }
    /// Slice of all weights.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    /// 2D velocity vector for direction `i` (only valid for D2Q9).
    ///
    /// # Panics
    ///
    /// Panics if the lattice is not D2Q9.
    pub fn velocity_2d(&self, i: usize) -> [i32; 2] {
        self.velocities_2d[i]
    }
    /// 3D velocity vector for direction `i` (only valid for D3Q19/D3Q27).
    ///
    /// # Panics
    ///
    /// Panics if the lattice is D2Q9.
    pub fn velocity_3d(&self, i: usize) -> [i32; 3] {
        self.velocities_3d[i]
    }
    /// Opposite direction index for direction `i`.
    pub fn opposite(&self, i: usize) -> usize {
        self.opposites[i]
    }
    /// Speed of sound squared in lattice units (1/3).
    pub fn cs2(&self) -> f64 {
        CS2
    }
}
/// A compact D2Q9 lattice Boltzmann grid with a single flat `Vec`f64` for
/// distribution functions.
///
/// Storage layout: `f\[idx * 9 + alpha\]` where `idx = y * nx + x` and
/// `alpha` is the direction index `0..9`.
#[derive(Debug, Clone)]
pub struct D2Q9Grid {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Distribution functions stored as `f\[idx * 9 + alpha\]`.
    pub f: Vec<f64>,
}
impl D2Q9Grid {
    /// Create a new D2Q9 grid initialised to equilibrium at rest with
    /// uniform density `rho0`.
    pub fn new(nx: usize, ny: usize, rho0: f64) -> Self {
        let n = nx * ny;
        let mut f = vec![0.0; n * 9];
        for idx in 0..n {
            for alpha in 0..9 {
                f[idx * 9 + alpha] = D2Q9_WEIGHTS[alpha] * rho0;
            }
        }
        Self { nx, ny, f }
    }
    /// Compute the equilibrium distribution value for a single direction.
    ///
    /// `alpha` is the direction index (0..9).
    ///
    /// ```text
    /// feq_alpha = w_alpha * rho * (1 + (e · u)/cs² + (e · u)²/(2 cs⁴) − u²/(2 cs²))
    /// ```
    pub fn feq(rho: f64, ux: f64, uy: f64, alpha: usize) -> f64 {
        let w = D2Q9_WEIGHTS[alpha];
        let cx = D2Q9_VELOCITIES[alpha][0] as f64;
        let cy = D2Q9_VELOCITIES[alpha][1] as f64;
        let eu = cx * ux + cy * uy;
        let u_sq = ux * ux + uy * uy;
        w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2))
    }
    /// Perform periodic streaming (pull scheme).
    ///
    /// Each cell pulls from its upstream neighbour for each direction,
    /// wrapping at domain boundaries.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let dst = (y * nx + x) * 9;
                for (alpha, vel) in D2Q9_VELOCITIES.iter().enumerate() {
                    let cx = vel[0];
                    let cy = vel[1];
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let src = (sy * nx + sx) * 9 + alpha;
                    self.f[dst + alpha] = f_old[src];
                }
            }
        }
    }
    /// BGK collision operator applied in-place.
    ///
    /// `omega = 1 / tau` is the relaxation frequency.
    pub fn collide_bgk(&mut self, omega: f64) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            let base = idx * 9;
            let mut rho = 0.0;
            let mut mx = 0.0;
            let mut my = 0.0;
            for (alpha, vel) in D2Q9_VELOCITIES.iter().enumerate() {
                let fi = self.f[base + alpha];
                rho += fi;
                mx += fi * vel[0] as f64;
                my += fi * vel[1] as f64;
            }
            let ux = if rho.abs() > 1e-15 { mx / rho } else { 0.0 };
            let uy = if rho.abs() > 1e-15 { my / rho } else { 0.0 };
            for alpha in 0..9 {
                let feq = Self::feq(rho, ux, uy, alpha);
                self.f[base + alpha] -= omega * (self.f[base + alpha] - feq);
            }
        }
    }
    /// Return the density field as a `Vec`f64` of length `nx * ny`.
    pub fn density_field(&self) -> Vec<f64> {
        let n = self.nx * self.ny;
        let mut rho = vec![0.0; n];
        for (idx, r) in rho.iter_mut().enumerate() {
            let base = idx * 9;
            *r = self.f[base..base + 9].iter().sum();
        }
        rho
    }
    /// Return the velocity field as a `Vec<[f64; 2]>` of length `nx * ny`.
    pub fn velocity_field(&self) -> Vec<[f64; 2]> {
        let n = self.nx * self.ny;
        let mut vel = vec![[0.0; 2]; n];
        for (idx, v) in vel.iter_mut().enumerate() {
            let base = idx * 9;
            let mut rho = 0.0;
            let mut mx = 0.0;
            let mut my = 0.0;
            for (alpha, cv) in D2Q9_VELOCITIES.iter().enumerate() {
                let fi = self.f[base + alpha];
                rho += fi;
                mx += fi * cv[0] as f64;
                my += fi * cv[1] as f64;
            }
            if rho.abs() > 1e-15 {
                *v = [mx / rho, my / rho];
            }
        }
        vel
    }
    /// Total mass (sum of all densities).
    pub fn total_mass(&self) -> f64 {
        self.density_field().iter().sum()
    }
}
/// Dimensions of a lattice domain (nx, ny, nz).
///
/// Provides index-conversion helpers for 3D row-major storage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LatticeDimensions {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Number of cells in z-direction.
    pub nz: usize,
}
impl LatticeDimensions {
    /// Create a new set of lattice dimensions.
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self { nx, ny, nz }
    }
    /// Total number of cells in the domain.
    pub fn total_cells(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    /// Convert 3D coordinates to a linear index (row-major: z fastest in ny*nx, then y, then x).
    ///
    /// Layout: `idx = z * ny * nx + y * nx + x`.
    pub fn idx_3d(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.ny * self.nx + y * self.nx + x
    }
    /// Convert a linear index back to 3D coordinates `(x, y, z)`.
    pub fn coords(&self, idx: usize) -> (usize, usize, usize) {
        let z = idx / (self.ny * self.nx);
        let rem = idx % (self.ny * self.nx);
        let y = rem / self.nx;
        let x = rem % self.nx;
        (x, y, z)
    }
}
/// Smagorinsky large-eddy-simulation turbulence model for LBM.
///
/// Adjusts the local relaxation rate omega based on the local strain rate,
/// using an effective viscosity: nu_eff = nu + (C_s * Delta)^2 * |S|.
pub struct SmagorinskyModel {
    /// Smagorinsky constant (typical: 0.1–0.2).
    pub cs: f64,
    /// Grid spacing (lattice units = 1).
    pub delta: f64,
    /// Base (molecular) kinematic viscosity.
    pub nu0: f64,
}
impl SmagorinskyModel {
    /// Create a new Smagorinsky model.
    pub fn new(cs: f64, delta: f64, nu0: f64) -> Self {
        Self { cs, delta, nu0 }
    }
    /// Compute the local effective omega from distribution functions at a cell.
    ///
    /// Returns the locally adjusted omega for BGK collision.
    pub fn effective_omega(&self, f: &[f64], base: usize) -> f64 {
        let mut rho = 0.0_f64;
        let mut mx = 0.0;
        let mut my = 0.0;
        for a in 0..9 {
            let fi = f[base + a];
            rho += fi;
            mx += fi * D2Q9_VELOCITIES[a][0] as f64;
            my += fi * D2Q9_VELOCITIES[a][1] as f64;
        }
        if rho < 1e-15 {
            return 1.0 / (self.nu0 / CS2 + 0.5);
        }
        let ux = mx / rho;
        let uy = my / rho;
        let mut pi_xx = 0.0_f64;
        let mut pi_xy = 0.0_f64;
        let mut pi_yy = 0.0_f64;
        for a in 0..9 {
            let feq_a = D2Q9Grid::feq(rho, ux, uy, a);
            let fneq = f[base + a] - feq_a;
            let cx = D2Q9_VELOCITIES[a][0] as f64;
            let cy = D2Q9_VELOCITIES[a][1] as f64;
            pi_xx += fneq * cx * cx;
            pi_xy += fneq * cx * cy;
            pi_yy += fneq * cy * cy;
        }
        let pi_sq = pi_xx * pi_xx + 2.0 * pi_xy * pi_xy + pi_yy * pi_yy;
        let s_sq = pi_sq.sqrt();
        let tau0 = self.nu0 / CS2 + 0.5;
        let tau_eff = 0.5
            * (tau0
                + (tau0 * tau0 + 18.0 * self.cs * self.cs * self.delta * self.delta * s_sq / rho)
                    .sqrt());
        1.0 / tau_eff
    }
}
/// Supported lattice velocity sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeType {
    /// 2D lattice with 9 velocities.
    D2Q9,
    /// 3D lattice with 19 velocities.
    D3Q19,
    /// 3D lattice with 27 velocities.
    D3Q27,
}
impl LatticeType {
    /// Number of discrete velocities for this lattice type.
    pub fn q(&self) -> usize {
        match self {
            Self::D2Q9 => 9,
            Self::D3Q19 => 19,
            Self::D3Q27 => 27,
        }
    }
    /// Number of spatial dimensions.
    pub fn dim(&self) -> usize {
        match self {
            Self::D2Q9 => 2,
            Self::D3Q19 | Self::D3Q27 => 3,
        }
    }
}
/// Typed D3Q19 velocity-set struct.
pub struct D3Q19;
impl D3Q19 {
    /// D3Q19 weights.
    pub const W: [f64; 19] = D3Q19_WEIGHTS;
    /// D3Q19 velocity vectors `[cx, cy, cz]`.
    pub const C: [[i32; 3]; 19] = D3Q19_VELOCITIES;
    /// Compute the D3Q19 equilibrium distribution for density `rho` and velocity `u`.
    pub fn equilibrium(rho: f64, u: [f64; 3]) -> [f64; 19] {
        let mut feq = [0.0_f64; 19];
        let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        for (i, fi) in feq.iter_mut().enumerate() {
            let cx = Self::C[i][0] as f64;
            let cy = Self::C[i][1] as f64;
            let cz = Self::C[i][2] as f64;
            let eu = cx * u[0] + cy * u[1] + cz * u[2];
            *fi = Self::W[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        feq
    }
}
/// A simple D2Q9 LBM grid storing distributions as `Vec<[f64; 9]>`.
///
/// Each element of `f` corresponds to one cell (row-major: `idx = y * nx + x`).
#[derive(Debug, Clone)]
pub struct LatticeD2Q9Grid {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Distribution functions: `f[y * nx + x]` is a `[f64; 9]` array.
    pub f: Vec<[f64; 9]>,
}
impl LatticeD2Q9Grid {
    /// Create a new grid initialised to equilibrium at rest with uniform density `rho0`.
    pub fn new(nx: usize, ny: usize) -> Self {
        let rho0 = 1.0_f64;
        let feq = D2Q9::equilibrium(rho0, 0.0, 0.0);
        Self {
            nx,
            ny,
            f: vec![feq; nx * ny],
        }
    }
    /// Linear index for cell (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Perform one full BGK step: collision then periodic pull-streaming.
    pub fn step(&mut self, omega: f64) {
        for k in 0..self.nx * self.ny {
            let rho = compute_rho(&self.f[k]);
            let vel = compute_velocity(&self.f[k], rho);
            bgk_collision(&mut self.f[k], rho, vel[0], vel[1], omega);
        }
        D2Q9::stream(&mut self.f, self.nx, self.ny);
    }
    /// Total mass (sum of all densities) — useful for conservation checks.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
    /// Compute the density field as `Vec`f64`.
    pub fn density_field(&self) -> Vec<f64> {
        self.f.iter().map(compute_rho).collect()
    }
}
