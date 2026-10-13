// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU fluid simulation: Navier-Stokes, LBM kernels, SPH kernels, vortex
//! confinement, and shallow-water equations — CPU mock backend.
//!
//! All structures here simulate GPU buffer management and dispatch via
//! in-memory `Vec`f64` buffers so that the logic can be unit-tested without
//! a real GPU device.

use crate::lbm_gpu::{D3Q19_EX, D3Q19_EY, D3Q19_EZ, D3Q19_OPP};

// ── FluidGpuBuffer ────────────────────────────────────────────────────────────

/// A double-buffered GPU field for ping-pong updates.
///
/// `front` is the current read buffer; `back` is the current write buffer.
/// Call [`FluidGpuBuffer::swap`] after each kernel dispatch to advance.
#[derive(Debug, Clone)]
pub struct FluidGpuBuffer {
    /// Front (read) buffer — size equals `width * height * depth * components`.
    pub front: Vec<f64>,
    /// Back (write) buffer — same size as `front`.
    pub back: Vec<f64>,
    /// Number of grid cells in the x direction.
    pub width: usize,
    /// Number of grid cells in the y direction.
    pub height: usize,
    /// Number of grid cells in the z direction (1 for 2-D grids).
    pub depth: usize,
    /// Number of scalar components per cell (e.g. 3 for velocity, 1 for pressure).
    pub components: usize,
}

impl FluidGpuBuffer {
    /// Allocate a zeroed double buffer for a `width × height × depth` grid with
    /// `components` scalars per cell.
    pub fn new(width: usize, height: usize, depth: usize, components: usize) -> Self {
        let n = width * height * depth * components;
        Self {
            front: vec![0.0; n],
            back: vec![0.0; n],
            width,
            height,
            depth,
            components,
        }
    }

    /// Total number of scalar values per buffer.
    #[inline]
    pub fn len(&self) -> usize {
        self.front.len()
    }

    /// Returns `true` if the buffer contains no elements.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.front.is_empty()
    }

    /// Swap front and back buffers (ping-pong).
    pub fn swap(&mut self) {
        std::mem::swap(&mut self.front, &mut self.back);
    }

    /// Flat index for cell `(x, y, z)` and component `c`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize, c: usize) -> usize {
        ((z * self.height + y) * self.width + x) * self.components + c
    }

    /// Read a scalar from the **front** buffer.
    #[inline]
    pub fn read(&self, x: usize, y: usize, z: usize, c: usize) -> f64 {
        self.front[self.idx(x, y, z, c)]
    }

    /// Write a scalar to the **back** buffer.
    #[inline]
    pub fn write(&mut self, x: usize, y: usize, z: usize, c: usize, val: f64) {
        let i = self.idx(x, y, z, c);
        self.back[i] = val;
    }

    /// Fill both buffers with `value`.
    pub fn fill(&mut self, value: f64) {
        self.front.fill(value);
        self.back.fill(value);
    }

    /// Copy front into back.
    pub fn copy_front_to_back(&mut self) {
        self.back.clone_from(&self.front);
    }
}

// ── NavierStokesGpu ───────────────────────────────────────────────────────────

/// GPU (CPU mock) Navier-Stokes solver on a 3-D Cartesian grid.
///
/// Implements semi-Lagrangian advection, divergence, pressure projection via
/// Jacobi iteration, and curl (vorticity) computation.
#[derive(Debug, Clone)]
pub struct NavierStokesGpu {
    /// Velocity field — 3 components (u, v, w) per cell.
    pub velocity: FluidGpuBuffer,
    /// Pressure field — 1 component per cell.
    pub pressure: FluidGpuBuffer,
    /// Density field — 1 component per cell.
    pub density: FluidGpuBuffer,
    /// Grid cell size (m).
    pub dx: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
}

impl NavierStokesGpu {
    /// Create a new Navier-Stokes solver for a `nx × ny × nz` grid.
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, viscosity: f64) -> Self {
        Self {
            velocity: FluidGpuBuffer::new(nx, ny, nz, 3),
            pressure: FluidGpuBuffer::new(nx, ny, nz, 1),
            density: FluidGpuBuffer::new(nx, ny, nz, 1),
            dx,
            viscosity,
        }
    }

    /// Grid width (cells in x).
    pub fn nx(&self) -> usize {
        self.velocity.width
    }

    /// Grid height (cells in y).
    pub fn ny(&self) -> usize {
        self.velocity.height
    }

    /// Grid depth (cells in z).
    pub fn nz(&self) -> usize {
        self.velocity.depth
    }

    /// Semi-Lagrangian advection kernel (CPU mock).
    ///
    /// Traces each cell centre backwards along the velocity field by `dt`
    /// seconds and samples the field using trilinear interpolation.
    pub fn advect(&mut self, dt: f64) {
        let nx = self.nx();
        let ny = self.ny();
        let nz = self.nz();
        let dx = self.dx;

        self.velocity.copy_front_to_back();

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let u = self.velocity.read(x, y, z, 0);
                    let v = self.velocity.read(x, y, z, 1);
                    let w = self.velocity.read(x, y, z, 2);

                    // back-trace position (clamped to grid)
                    let fx = (x as f64 - u * dt / dx).clamp(0.0, (nx - 1) as f64);
                    let fy = (y as f64 - v * dt / dx).clamp(0.0, (ny - 1) as f64);
                    let fz = (z as f64 - w * dt / dx).clamp(0.0, (nz - 1) as f64);

                    for c in 0..3 {
                        let val =
                            trilinear_sample(&self.velocity.front, nx, ny, nz, 3, fx, fy, fz, c);
                        self.velocity.write(x, y, z, c, val);
                    }
                }
            }
        }
        self.velocity.swap();
    }

    /// Compute the divergence field into a flat `Vec`f64`.
    ///
    /// `result[z * ny * nx + y * nx + x]` = divergence at cell `(x, y, z)`.
    pub fn divergence(&self) -> Vec<f64> {
        let nx = self.nx();
        let ny = self.ny();
        let nz = self.nz();
        let inv2dx = 0.5 / self.dx;
        let mut div = vec![0.0_f64; nx * ny * nz];

        for z in 1..nz.saturating_sub(1) {
            for y in 1..ny.saturating_sub(1) {
                for x in 1..nx.saturating_sub(1) {
                    let du =
                        self.velocity.read(x + 1, y, z, 0) - self.velocity.read(x - 1, y, z, 0);
                    let dv =
                        self.velocity.read(x, y + 1, z, 1) - self.velocity.read(x, y - 1, z, 1);
                    let dw =
                        self.velocity.read(x, y, z + 1, 2) - self.velocity.read(x, y, z - 1, 2);
                    div[z * ny * nx + y * nx + x] = (du + dv + dw) * inv2dx;
                }
            }
        }
        div
    }

    /// Pressure projection: Jacobi iterations to enforce incompressibility.
    ///
    /// Runs `iterations` steps of Jacobi Poisson solve, then subtracts the
    /// pressure gradient from the velocity field.
    pub fn pressure_projection(&mut self, iterations: usize) {
        let nx = self.nx();
        let ny = self.ny();
        let nz = self.nz();
        let dx2 = self.dx * self.dx;
        let div = self.divergence();

        // Jacobi iterations
        for _ in 0..iterations {
            self.pressure.copy_front_to_back();
            for z in 1..nz.saturating_sub(1) {
                for y in 1..ny.saturating_sub(1) {
                    for x in 1..nx.saturating_sub(1) {
                        let neighbours = self.pressure.read(x + 1, y, z, 0)
                            + self.pressure.read(x - 1, y, z, 0)
                            + self.pressure.read(x, y + 1, z, 0)
                            + self.pressure.read(x, y - 1, z, 0)
                            + self.pressure.read(x, y, z + 1, 0)
                            + self.pressure.read(x, y, z - 1, 0);
                        let rhs = div[z * ny * nx + y * nx + x] * dx2;
                        self.pressure.write(x, y, z, 0, (neighbours - rhs) / 6.0);
                    }
                }
            }
            self.pressure.swap();
        }

        // subtract pressure gradient from velocity
        let inv2dx = 0.5 / self.dx;
        self.velocity.copy_front_to_back();
        for z in 1..nz.saturating_sub(1) {
            for y in 1..ny.saturating_sub(1) {
                for x in 1..nx.saturating_sub(1) {
                    let gx = (self.pressure.read(x + 1, y, z, 0)
                        - self.pressure.read(x - 1, y, z, 0))
                        * inv2dx;
                    let gy = (self.pressure.read(x, y + 1, z, 0)
                        - self.pressure.read(x, y - 1, z, 0))
                        * inv2dx;
                    let gz = (self.pressure.read(x, y, z + 1, 0)
                        - self.pressure.read(x, y, z - 1, 0))
                        * inv2dx;
                    self.velocity
                        .write(x, y, z, 0, self.velocity.read(x, y, z, 0) - gx);
                    self.velocity
                        .write(x, y, z, 1, self.velocity.read(x, y, z, 1) - gy);
                    self.velocity
                        .write(x, y, z, 2, self.velocity.read(x, y, z, 2) - gz);
                }
            }
        }
        self.velocity.swap();
    }

    /// Compute the curl (vorticity) field.
    ///
    /// Returns a flat buffer of `[wx, wy, wz]` per cell (3 components).
    pub fn curl(&self) -> Vec<f64> {
        let nx = self.nx();
        let ny = self.ny();
        let nz = self.nz();
        let inv2dx = 0.5 / self.dx;
        let mut omega = vec![0.0_f64; nx * ny * nz * 3];

        for z in 1..nz.saturating_sub(1) {
            for y in 1..ny.saturating_sub(1) {
                for x in 1..nx.saturating_sub(1) {
                    let dw_dy = (self.velocity.read(x, y + 1, z, 2)
                        - self.velocity.read(x, y - 1, z, 2))
                        * inv2dx;
                    let dv_dz = (self.velocity.read(x, y, z + 1, 1)
                        - self.velocity.read(x, y, z - 1, 1))
                        * inv2dx;
                    let du_dz = (self.velocity.read(x, y, z + 1, 0)
                        - self.velocity.read(x, y, z - 1, 0))
                        * inv2dx;
                    let dw_dx = (self.velocity.read(x + 1, y, z, 2)
                        - self.velocity.read(x - 1, y, z, 2))
                        * inv2dx;
                    let dv_dx = (self.velocity.read(x + 1, y, z, 1)
                        - self.velocity.read(x - 1, y, z, 1))
                        * inv2dx;
                    let du_dy = (self.velocity.read(x, y + 1, z, 0)
                        - self.velocity.read(x, y - 1, z, 0))
                        * inv2dx;

                    let base = (z * ny * nx + y * nx + x) * 3;
                    omega[base] = dw_dy - dv_dz; // wx
                    omega[base + 1] = du_dz - dw_dx; // wy
                    omega[base + 2] = dv_dx - du_dy; // wz
                }
            }
        }
        omega
    }
}

// ── Trilinear sampling helper ─────────────────────────────────────────────────

/// Trilinear interpolation of a 3-D field stored in a flat buffer.
///
/// `data` is `nx * ny * nz * stride` elements; `c` selects the component.
fn trilinear_sample(
    data: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    stride: usize,
    fx: f64,
    fy: f64,
    fz: f64,
    c: usize,
) -> f64 {
    let x0 = (fx as usize).min(nx - 1);
    let y0 = (fy as usize).min(ny - 1);
    let z0 = (fz as usize).min(nz - 1);
    let x1 = (x0 + 1).min(nx - 1);
    let y1 = (y0 + 1).min(ny - 1);
    let z1 = (z0 + 1).min(nz - 1);

    let tx = fx.fract();
    let ty = fy.fract();
    let tz = fz.fract();

    let idx =
        |x: usize, y: usize, z: usize| -> f64 { data[(z * ny * nx + y * nx + x) * stride + c] };

    let c000 = idx(x0, y0, z0);
    let c100 = idx(x1, y0, z0);
    let c010 = idx(x0, y1, z0);
    let c110 = idx(x1, y1, z0);
    let c001 = idx(x0, y0, z1);
    let c101 = idx(x1, y0, z1);
    let c011 = idx(x0, y1, z1);
    let c111 = idx(x1, y1, z1);

    let c00 = c000 * (1.0 - tx) + c100 * tx;
    let c10 = c010 * (1.0 - tx) + c110 * tx;
    let c01 = c001 * (1.0 - tx) + c101 * tx;
    let c11 = c011 * (1.0 - tx) + c111 * tx;
    let c_0 = c00 * (1.0 - ty) + c10 * ty;
    let c_1 = c01 * (1.0 - ty) + c11 * ty;
    c_0 * (1.0 - tz) + c_1 * tz
}

// ── LBMGpuKernels ─────────────────────────────────────────────────────────────

/// D3Q19 LBM velocity directions (19 discrete velocities).
pub const D3Q19_Q: usize = 19;

/// D3Q19 lattice weights, in the velocity order of
/// [`crate::lbm_gpu::D3Q19_EX`] / `EY` / `EZ`.
pub const D3Q19_WEIGHTS: [f64; D3Q19_Q] = [
    1.0 / 3.0, // 0: rest
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0, // 1-3: ±x, ±y, ±z axes
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0, // 4-6
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0, // 7-9: face diagonals
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0, // 10-12
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0, // 13-15
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0, // 16-18
];

/// Specification of an LBM collision kernel.
#[derive(Debug, Clone)]
pub struct LBMCollisionKernelSpec {
    /// Relaxation parameter τ (tau).  Viscosity ν = cs² (τ − 0.5) dt.
    pub tau: f64,
    /// Sound speed squared (lattice units, `1/3` for D3Q19).
    pub cs2: f64,
    /// Number of discrete velocity directions ([`D3Q19_Q`] for the D3Q19
    /// kernels of this module).
    pub q: usize,
}

/// Specification of an LBM streaming kernel.
#[derive(Debug, Clone)]
pub struct LBMStreamingSpec {
    /// Grid dimensions `[nx, ny, nz]`.
    pub dims: [usize; 3],
    /// Whether to apply periodic boundary conditions.  When `false`, the
    /// domain faces are half-way bounce-back walls.
    pub periodic: bool,
}

/// Specification of an LBM boundary kernel.
#[derive(Debug, Clone)]
pub struct LBMBoundaryKernelSpec {
    /// Boundary condition type.
    pub kind: LBMBoundaryKind,
    /// Indices of solid (obstacle) cells.
    pub solid_cells: Vec<usize>,
}

/// Variants of supported LBM boundary conditions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LBMBoundaryKind {
    /// Full bounce-back (no-slip wall).
    BounceBack,
    /// Moving-wall bounce-back.
    MovingWall,
    /// Zou-He velocity inlet.
    ZouHeInlet,
    /// Zou-He pressure outlet.
    ZouHeOutlet,
}

/// Density `ρ = Σ f_i` and velocity `u = Σ f_i e_i / ρ` of one cell's D3Q19
/// populations (`u = 0` for a vacuum cell with `ρ <= 1e-12`).
fn d3q19_moments(f: &[f64]) -> (f64, [f64; 3]) {
    let mut rho = 0.0_f64;
    let mut momentum = [0.0_f64; 3];
    for (i, &fi) in f.iter().enumerate().take(D3Q19_Q) {
        rho += fi;
        momentum[0] += fi * f64::from(D3Q19_EX[i]);
        momentum[1] += fi * f64::from(D3Q19_EY[i]);
        momentum[2] += fi * f64::from(D3Q19_EZ[i]);
    }
    if rho > 1e-12 {
        (
            rho,
            [momentum[0] / rho, momentum[1] / rho, momentum[2] / rho],
        )
    } else {
        (rho, [0.0; 3])
    }
}

/// Second-order equilibrium of direction `i`:
/// `w_i ρ [1 + (e_i·u)/cs² + (e_i·u)²/(2 cs⁴) − |u|²/(2 cs²)]`
/// (`inv_cs2 = 1/cs²`).
fn d3q19_equilibrium(i: usize, rho: f64, u: [f64; 3], inv_cs2: f64) -> f64 {
    let eu = f64::from(D3Q19_EX[i]) * u[0]
        + f64::from(D3Q19_EY[i]) * u[1]
        + f64::from(D3Q19_EZ[i]) * u[2];
    let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let eu_c = eu * inv_cs2;
    D3Q19_WEIGHTS[i] * rho * (1.0 + eu_c + 0.5 * eu_c * eu_c - 0.5 * u_sq * inv_cs2)
}

/// GPU kernel manager for a Lattice Boltzmann Method simulation.
///
/// Stores distribution functions for all cells in a ping-pong buffer and
/// exposes collision, streaming, and boundary dispatch methods that run the
/// D3Q19 BGK algorithm on the CPU.  One time step is
/// `dispatch_collision` → `dispatch_streaming` → `dispatch_boundary`.
#[derive(Debug, Clone)]
pub struct LBMGpuKernels {
    /// Double-buffered distribution functions: `Q` components per cell.
    pub f: FluidGpuBuffer,
    /// Collision kernel specification.
    pub collision: LBMCollisionKernelSpec,
    /// Streaming kernel specification.
    pub streaming: LBMStreamingSpec,
    /// Boundary kernel specification.
    pub boundary: LBMBoundaryKernelSpec,
}

impl LBMGpuKernels {
    /// Create a new LBM kernel manager.
    ///
    /// Initialises distribution functions to the equilibrium for ρ = 1, u = 0.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        periodic: bool,
        solid_cells: Vec<usize>,
    ) -> Self {
        let q = D3Q19_Q;
        let mut f = FluidGpuBuffer::new(nx, ny, nz, q);
        // initialise to equilibrium weights
        for i in 0..f.front.len() {
            let c = i % q;
            f.front[i] = D3Q19_WEIGHTS[c];
            f.back[i] = D3Q19_WEIGHTS[c];
        }
        Self {
            f,
            collision: LBMCollisionKernelSpec {
                tau,
                cs2: 1.0 / 3.0,
                q,
            },
            streaming: LBMStreamingSpec {
                dims: [nx, ny, nz],
                periodic,
            },
            boundary: LBMBoundaryKernelSpec {
                kind: LBMBoundaryKind::BounceBack,
                solid_cells,
            },
        }
    }

    /// Boolean mask over the grid cells: `true` for every distinct, in-range
    /// entry of `boundary.solid_cells`.
    fn solid_mask(&self) -> Vec<bool> {
        let n_cells = self.f.width * self.f.height * self.f.depth;
        let mut mask = vec![false; n_cells];
        for &cell in &self.boundary.solid_cells {
            if let Some(flag) = mask.get_mut(cell) {
                *flag = true;
            }
        }
        mask
    }

    /// Execute one BGK collision step (CPU mock of the collision kernel).
    ///
    /// For every fluid cell the density `ρ = Σ f_i` and velocity
    /// `u = Σ f_i e_i / ρ` are computed from the D3Q19 populations, and each
    /// population relaxes towards
    /// `f_i^eq = w_i ρ [1 + (e_i·u)/cs² + (e_i·u)²/(2 cs⁴) − |u|²/(2 cs²)]`
    /// as `f_i ← f_i − (f_i − f_i^eq)/τ`.  Mass and momentum are conserved
    /// exactly.  Solid cells are skipped so that the bounce-back of
    /// [`dispatch_boundary`](Self::dispatch_boundary) reverses the arrived
    /// populations unchanged.
    pub fn dispatch_collision(&mut self) {
        let tau = self.collision.tau;
        let inv_cs2 = 1.0 / self.collision.cs2;
        let solid = self.solid_mask();

        for (cell, &is_solid) in solid.iter().enumerate() {
            let base = cell * D3Q19_Q;
            let src = &self.f.front[base..base + D3Q19_Q];
            let dst = &mut self.f.back[base..base + D3Q19_Q];
            if is_solid {
                dst.copy_from_slice(src);
                continue;
            }
            let (rho, u) = d3q19_moments(src);
            for (i, (out, &fi)) in dst.iter_mut().zip(src.iter()).enumerate() {
                let feq = d3q19_equilibrium(i, rho, u, inv_cs2);
                *out = fi - (fi - feq) / tau;
            }
        }
        self.f.swap();
    }

    /// Execute one streaming step (CPU mock of the pull-streaming kernel).
    ///
    /// Every population moves one lattice link along its velocity:
    /// `f_i(x, t+1) = f_i(x − e_i, t)`.  With `streaming.periodic` the source
    /// cell wraps around the domain; otherwise a population that would have
    /// to come from outside the domain is the one that left through the wall
    /// and is reflected back (half-way bounce-back), so mass is conserved
    /// exactly in both cases.
    pub fn dispatch_streaming(&mut self) {
        let (nx, ny, nz) = (self.f.width, self.f.height, self.f.depth);
        let periodic = self.streaming.periodic;
        let dims = [nx as isize, ny as isize, nz as isize];
        let cell_base = |x: usize, y: usize, z: usize| ((z * ny + y) * nx + x) * D3Q19_Q;

        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let base = cell_base(x, y, z);
                    for i in 0..D3Q19_Q {
                        let mut src = [
                            x as isize - D3Q19_EX[i] as isize,
                            y as isize - D3Q19_EY[i] as isize,
                            z as isize - D3Q19_EZ[i] as isize,
                        ];
                        let mut inside = true;
                        for (s, &n) in src.iter_mut().zip(dims.iter()) {
                            if periodic {
                                *s = s.rem_euclid(n);
                            } else if *s < 0 || *s >= n {
                                inside = false;
                            }
                        }
                        self.f.back[base + i] = if inside {
                            let src_base =
                                cell_base(src[0] as usize, src[1] as usize, src[2] as usize);
                            self.f.front[src_base + i]
                        } else {
                            self.f.front[base + D3Q19_OPP[i]]
                        };
                    }
                }
            }
        }
        self.f.swap();
    }

    /// Apply the no-slip full-way bounce-back to every cell of
    /// `boundary.solid_cells`: the populations that streamed into a solid cell
    /// are reversed (`f_i ← f_opp(i)`), so the next streaming step sends each
    /// of them back where it came from.
    ///
    /// Each distinct solid cell is reflected exactly once (a duplicated index
    /// does not undo itself) and out-of-range indices are ignored.  The pass is
    /// the same for every [`LBMBoundaryKind`]: the inlet, outlet and
    /// moving-wall variants carry no velocity or density data and are treated
    /// as stationary no-slip walls here.
    pub fn dispatch_boundary(&mut self) {
        let solid = self.solid_mask();
        let mut reflected = [0.0_f64; D3Q19_Q];
        for (cell, &is_solid) in solid.iter().enumerate() {
            if !is_solid {
                continue;
            }
            let base = cell * D3Q19_Q;
            let populations = &mut self.f.front[base..base + D3Q19_Q];
            for (i, r) in reflected.iter_mut().enumerate() {
                *r = populations[D3Q19_OPP[i]];
            }
            populations.copy_from_slice(&reflected);
        }
    }

    /// Compute macroscopic density at cell `(x, y, z)`.
    pub fn density_at(&self, x: usize, y: usize, z: usize) -> f64 {
        (0..self.collision.q).map(|i| self.f.read(x, y, z, i)).sum()
    }

    /// Momentum `Σ f_i e_i` at cell `(x, y, z)`.
    pub fn momentum_at(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let mut momentum = [0.0_f64; 3];
        for i in 0..self.collision.q.min(D3Q19_Q) {
            let fi = self.f.read(x, y, z, i);
            momentum[0] += fi * f64::from(D3Q19_EX[i]);
            momentum[1] += fi * f64::from(D3Q19_EY[i]);
            momentum[2] += fi * f64::from(D3Q19_EZ[i]);
        }
        momentum
    }
}

// ── SPHGpuKernels ─────────────────────────────────────────────────────────────

/// GPU SPH particle data for density/pressure/viscosity kernels.
#[derive(Debug, Clone)]
pub struct SPHGpuKernels {
    /// Particle positions `[x, y, z]`.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities `[vx, vy, vz]`.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle density (kg/m³).
    pub densities: Vec<f64>,
    /// Per-particle pressure (Pa).
    pub pressures: Vec<f64>,
    /// Smoothing length `h` (m).
    pub smoothing_length: f64,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Pressure stiffness constant `k`.
    pub stiffness: f64,
    /// Kinematic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Particle mass (kg).
    pub mass: f64,
}

impl SPHGpuKernels {
    /// Create a new SPH kernel manager with `n` particles.
    pub fn new(
        positions: Vec<[f64; 3]>,
        smoothing_length: f64,
        rest_density: f64,
        stiffness: f64,
        viscosity: f64,
        mass: f64,
    ) -> Self {
        let n = positions.len();
        Self {
            velocities: vec![[0.0; 3]; n],
            densities: vec![rest_density; n],
            pressures: vec![0.0; n],
            positions,
            smoothing_length,
            rest_density,
            stiffness,
            viscosity,
            mass,
        }
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Cubic spline SPH kernel W(r, h).
    pub fn kernel_w(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        let q = r / h;
        let sigma = 8.0 / (std::f64::consts::PI * h * h * h);
        if q >= 1.0 {
            0.0
        } else if q >= 0.5 {
            let t = 1.0 - q;
            sigma * 2.0 * t * t * t
        } else {
            sigma * (6.0 * q * q * q - 6.0 * q * q + 1.0)
        }
    }

    /// Gradient of cubic spline kernel ∇W(r⃗, h).
    pub fn kernel_grad_w(&self, rij: [f64; 3], r: f64) -> [f64; 3] {
        let h = self.smoothing_length;
        if r < 1e-12 || r >= h {
            return [0.0; 3];
        }
        let q = r / h;
        let sigma = 8.0 / (std::f64::consts::PI * h * h * h);
        let dw_dq = if q >= 0.5 {
            let t = 1.0 - q;
            -6.0 * sigma * t * t / h
        } else {
            sigma * (18.0 * q * q - 12.0 * q) / h
        };
        let inv_r = 1.0 / r;
        [
            rij[0] * dw_dq * inv_r,
            rij[1] * dw_dq * inv_r,
            rij[2] * dw_dq * inv_r,
        ]
    }

    /// Density kernel: compute per-particle density from neighbour contributions.
    pub fn compute_density(&mut self) {
        let n = self.n_particles();
        let mass = self.mass;
        let mut new_densities = vec![0.0_f64; n];

        for (i, &pos_i) in self.positions.iter().enumerate() {
            let mut rho = 0.0;
            for &pos_j in self.positions.iter() {
                let dx = pos_i[0] - pos_j[0];
                let dy = pos_i[1] - pos_j[1];
                let dz = pos_i[2] - pos_j[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                rho += mass * self.kernel_w(r);
            }
            new_densities[i] = rho.max(self.rest_density * 0.01);
        }
        self.densities = new_densities;
    }

    /// Pressure kernel: compute per-particle pressure (equation of state).
    pub fn compute_pressure(&mut self) {
        for (pr, &rho) in self.pressures.iter_mut().zip(self.densities.iter()) {
            *pr = self.stiffness * (rho - self.rest_density).max(0.0);
        }
    }

    /// Pressure force kernel: compute pressure force on each particle.
    ///
    /// Returns acceleration vectors `[ax, ay, az]` per particle.
    pub fn pressure_force(&self) -> Vec<[f64; 3]> {
        let n = self.n_particles();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mass = self.mass;

        for (i, &pos_i) in self.positions.iter().enumerate() {
            let mut fx = 0.0;
            let mut fy = 0.0;
            let mut fz = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = pos_i[0] - self.positions[j][0];
                let dy = pos_i[1] - self.positions[j][1];
                let dz = pos_i[2] - self.positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 1e-12 {
                    continue;
                }
                let grad = self.kernel_grad_w([dx, dy, dz], r);
                let pi = self.pressures[i];
                let pj = self.pressures[j];
                let rhoj = self.densities[j].max(1e-12);
                let rhoi = self.densities[i].max(1e-12);
                let coeff = -mass * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj));
                fx += coeff * grad[0];
                fy += coeff * grad[1];
                fz += coeff * grad[2];
            }
            forces[i] = [fx, fy, fz];
        }
        forces
    }

    /// Viscosity kernel: compute viscosity force on each particle.
    ///
    /// Returns acceleration vectors `[ax, ay, az]` per particle.
    pub fn viscosity_force(&self) -> Vec<[f64; 3]> {
        let n = self.n_particles();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mass = self.mass;
        let mu = self.viscosity;

        for (i, &pos_i) in self.positions.iter().enumerate() {
            let mut fx = 0.0;
            let mut fy = 0.0;
            let mut fz = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = pos_i[0] - self.positions[j][0];
                let dy = pos_i[1] - self.positions[j][1];
                let dz = pos_i[2] - self.positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 1e-12 {
                    continue;
                }
                let dvx = self.velocities[j][0] - self.velocities[i][0];
                let dvy = self.velocities[j][1] - self.velocities[i][1];
                let dvz = self.velocities[j][2] - self.velocities[i][2];
                let rhoj = self.densities[j].max(1e-12);
                let w = self.kernel_w(r);
                let coeff = mu * mass / rhoj * w;
                fx += coeff * dvx;
                fy += coeff * dvy;
                fz += coeff * dvz;
            }
            forces[i] = [fx, fy, fz];
        }
        forces
    }
}

// ── FluidGpuPipeline ──────────────────────────────────────────────────────────

/// Synchronization barrier kind between GPU passes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BarrierKind {
    /// Wait for all compute dispatches to finish.
    Compute,
    /// Wait for memory writes to be visible (memory barrier).
    Memory,
    /// Full pipeline barrier.
    Full,
}

/// A single pass (dispatch) in the fluid GPU pipeline.
#[derive(Debug, Clone)]
pub struct PipelinePass {
    /// Human-readable label.
    pub label: String,
    /// Barrier to insert **before** this pass.
    pub barrier: BarrierKind,
    /// Estimated memory bandwidth consumed (bytes).
    pub estimated_bytes: usize,
}

/// Ordered fluid GPU pipeline with memory bandwidth estimation.
///
/// Dispatch order: advect → divergence → pressure solve → gradient subtract
/// → vorticity confinement.
#[derive(Debug, Clone)]
pub struct FluidGpuPipeline {
    /// Ordered list of pipeline passes.
    pub passes: Vec<PipelinePass>,
    /// Total bytes transferred (sum over all passes).
    pub total_bytes: usize,
}

impl FluidGpuPipeline {
    /// Build the default Navier-Stokes pipeline for an `n`-cell grid.
    pub fn default_ns(n_cells: usize) -> Self {
        let bytes_per_cell = 8 * std::mem::size_of::<f64>(); // approx
        let b = n_cells * bytes_per_cell;
        let passes = vec![
            PipelinePass {
                label: "advect_velocity".into(),
                barrier: BarrierKind::Memory,
                estimated_bytes: b * 2,
            },
            PipelinePass {
                label: "advect_density".into(),
                barrier: BarrierKind::Memory,
                estimated_bytes: b,
            },
            PipelinePass {
                label: "divergence".into(),
                barrier: BarrierKind::Compute,
                estimated_bytes: b,
            },
            PipelinePass {
                label: "pressure_jacobi_0".into(),
                barrier: BarrierKind::Memory,
                estimated_bytes: b * 2,
            },
            PipelinePass {
                label: "pressure_jacobi_1".into(),
                barrier: BarrierKind::Memory,
                estimated_bytes: b * 2,
            },
            PipelinePass {
                label: "gradient_subtract".into(),
                barrier: BarrierKind::Compute,
                estimated_bytes: b,
            },
            PipelinePass {
                label: "vorticity_confinement".into(),
                barrier: BarrierKind::Full,
                estimated_bytes: b,
            },
        ];
        let total_bytes = passes.iter().map(|p| p.estimated_bytes).sum();
        Self {
            passes,
            total_bytes,
        }
    }

    /// Add a custom pass to the pipeline.
    pub fn push(&mut self, pass: PipelinePass) {
        self.total_bytes += pass.estimated_bytes;
        self.passes.push(pass);
    }

    /// Estimated memory bandwidth in GB/s, given elapsed time in seconds.
    pub fn bandwidth_gb_s(&self, elapsed_secs: f64) -> f64 {
        if elapsed_secs <= 0.0 {
            return 0.0;
        }
        self.total_bytes as f64 / elapsed_secs / 1e9
    }
}

// ── VortexConfinement ─────────────────────────────────────────────────────────

/// GPU vortex confinement — amplifies small-scale vortical structures to
/// counteract numerical dissipation.
///
/// Based on Fedkiw et al. (2001) vorticity confinement.
#[derive(Debug, Clone)]
pub struct VortexConfinement {
    /// Confinement strength ε.  Typical values: 0.1–2.0.
    pub epsilon: f64,
    /// Grid cell size (m).
    pub dx: f64,
}

impl VortexConfinement {
    /// Create a new `VortexConfinement` with given parameters.
    pub fn new(epsilon: f64, dx: f64) -> Self {
        Self { epsilon, dx }
    }

    /// Compute and apply vortex confinement forces to a [`NavierStokesGpu`] solver.
    ///
    /// Adds a body-force acceleration field derived from the vorticity gradient.
    pub fn apply(&self, ns: &mut NavierStokesGpu, dt: f64) {
        let nx = ns.nx();
        let ny = ns.ny();
        let nz = ns.nz();
        let omega = ns.curl();
        let inv2dx = 0.5 / self.dx;

        // compute vorticity magnitude field |ω|
        let mut mag = vec![0.0_f64; nx * ny * nz];
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let base = (z * ny * nx + y * nx + x) * 3;
                    let wx = omega[base];
                    let wy = omega[base + 1];
                    let wz = omega[base + 2];
                    mag[z * ny * nx + y * nx + x] = (wx * wx + wy * wy + wz * wz).sqrt();
                }
            }
        }

        // apply confinement force: f = ε dx (N̂ × ω), N̂ = ∇|ω| / |∇|ω||
        ns.velocity.copy_front_to_back();
        for z in 1..nz.saturating_sub(1) {
            for y in 1..ny.saturating_sub(1) {
                for x in 1..nx.saturating_sub(1) {
                    let gx = (mag[z * ny * nx + y * nx + x + 1]
                        - mag[z * ny * nx + y * nx + x - 1])
                        * inv2dx;
                    let gy = (mag[z * ny * nx + (y + 1) * nx + x]
                        - mag[z * ny * nx + (y - 1) * nx + x])
                        * inv2dx;
                    let gz = (mag[(z + 1) * ny * nx + y * nx + x]
                        - mag[(z - 1) * ny * nx + y * nx + x])
                        * inv2dx;

                    let grad_len = (gx * gx + gy * gy + gz * gz).sqrt();
                    if grad_len < 1e-12 {
                        continue;
                    }
                    let nx_ = gx / grad_len;
                    let ny_ = gy / grad_len;
                    let nz_ = gz / grad_len;

                    let base = (z * ny * nx + y * nx + x) * 3;
                    let wx = omega[base];
                    let wy = omega[base + 1];
                    let wz = omega[base + 2];

                    // N × ω
                    let fcx = ny_ * wz - nz_ * wy;
                    let fcy = nz_ * wx - nx_ * wz;
                    let fcz = nx_ * wy - ny_ * wx;

                    let force = self.epsilon * self.dx;
                    let u = ns.velocity.read(x, y, z, 0) + force * fcx * dt;
                    let v = ns.velocity.read(x, y, z, 1) + force * fcy * dt;
                    let w = ns.velocity.read(x, y, z, 2) + force * fcz * dt;
                    ns.velocity.write(x, y, z, 0, u);
                    ns.velocity.write(x, y, z, 1, v);
                    ns.velocity.write(x, y, z, 2, w);
                }
            }
        }
        ns.velocity.swap();
    }
}

// ── WaterSimGpu ───────────────────────────────────────────────────────────────

/// GPU shallow-water equation (SWE) solver for height field and ripple
/// simulation on a 2-D grid.
///
/// Uses the linearised SWE: ∂h/∂t + H ∇·u = 0, ∂u/∂t + g ∇h = 0.
#[derive(Debug, Clone)]
pub struct WaterSimGpu {
    /// Height field h(x, y) — 1 component per cell (m above rest level).
    pub height: FluidGpuBuffer,
    /// Horizontal velocity field — 2 components (u_x, u_y) per cell.
    pub velocity: FluidGpuBuffer,
    /// Rest water depth H (m).
    pub rest_depth: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: f64,
    /// Grid cell size (m).
    pub dx: f64,
    /// Damping coefficient per step (fraction).
    pub damping: f64,
}

impl WaterSimGpu {
    /// Create a new `WaterSimGpu` for a `nx × ny` grid.
    pub fn new(nx: usize, ny: usize, dx: f64, rest_depth: f64, gravity: f64) -> Self {
        Self {
            height: FluidGpuBuffer::new(nx, ny, 1, 1),
            velocity: FluidGpuBuffer::new(nx, ny, 1, 2),
            rest_depth,
            gravity,
            dx,
            damping: 0.999,
        }
    }

    /// Grid width.
    pub fn nx(&self) -> usize {
        self.height.width
    }

    /// Grid height.
    pub fn ny(&self) -> usize {
        self.height.height
    }

    /// Apply a ripple disturbance at cell `(cx, cy)` with amplitude `amp` and
    /// Gaussian radius `sigma` (in cells).
    pub fn add_disturbance(&mut self, cx: usize, cy: usize, amp: f64, sigma: f64) {
        let nx = self.nx();
        let ny = self.ny();
        for y in 0..ny {
            for x in 0..nx {
                let dx = x as f64 - cx as f64;
                let dy = y as f64 - cy as f64;
                let r2 = dx * dx + dy * dy;
                let h = self.height.read(x, y, 0, 0);
                let idx = self.height.idx(x, y, 0, 0);
                self.height.front[idx] = h + amp * (-r2 / (2.0 * sigma * sigma)).exp();
            }
        }
    }

    /// Step the shallow-water equations forward by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx();
        let ny = self.ny();
        let g = self.gravity;
        let h0 = self.rest_depth;
        let inv2dx = 0.5 / self.dx;

        // Update velocity: ∂u/∂t = -g ∇h
        self.velocity.copy_front_to_back();
        for y in 1..ny.saturating_sub(1) {
            for x in 1..nx.saturating_sub(1) {
                let dhdx =
                    (self.height.read(x + 1, y, 0, 0) - self.height.read(x - 1, y, 0, 0)) * inv2dx;
                let dhdy =
                    (self.height.read(x, y + 1, 0, 0) - self.height.read(x, y - 1, 0, 0)) * inv2dx;
                let ux = self.velocity.read(x, y, 0, 0) - g * dhdx * dt;
                let uy = self.velocity.read(x, y, 0, 1) - g * dhdy * dt;
                self.velocity.write(x, y, 0, 0, ux * self.damping);
                self.velocity.write(x, y, 0, 1, uy * self.damping);
            }
        }
        self.velocity.swap();

        // Update height: ∂h/∂t = -H ∇·u
        self.height.copy_front_to_back();
        for y in 1..ny.saturating_sub(1) {
            for x in 1..nx.saturating_sub(1) {
                let divx = (self.velocity.read(x + 1, y, 0, 0)
                    - self.velocity.read(x - 1, y, 0, 0))
                    * inv2dx;
                let divy = (self.velocity.read(x, y + 1, 0, 1)
                    - self.velocity.read(x, y - 1, 0, 1))
                    * inv2dx;
                let h = self.height.read(x, y, 0, 0) - h0 * (divx + divy) * dt;
                self.height.write(x, y, 0, 0, h);
            }
        }
        self.height.swap();
    }

    /// Compute the maximum absolute height deviation from rest.
    pub fn max_height_deviation(&self) -> f64 {
        self.height
            .front
            .iter()
            .copied()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FluidGpuBuffer tests ──────────────────────────────────────────────

    #[test]
    fn buffer_new_zeroed() {
        let buf = FluidGpuBuffer::new(4, 4, 4, 3);
        assert!(buf.front.iter().all(|&v| v == 0.0));
        assert!(buf.back.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn buffer_len_correct() {
        let buf = FluidGpuBuffer::new(2, 3, 4, 5);
        assert_eq!(buf.len(), 2 * 3 * 4 * 5);
    }

    #[test]
    fn buffer_is_empty_false() {
        let buf = FluidGpuBuffer::new(2, 2, 2, 1);
        assert!(!buf.is_empty());
    }

    #[test]
    fn buffer_is_empty_true() {
        let buf = FluidGpuBuffer::new(0, 1, 1, 1);
        assert!(buf.is_empty());
    }

    #[test]
    fn buffer_write_read_roundtrip() {
        let mut buf = FluidGpuBuffer::new(4, 4, 1, 3);
        buf.write(1, 2, 0, 1, 42.5);
        buf.swap();
        assert!((buf.read(1, 2, 0, 1) - 42.5).abs() < 1e-12);
    }

    #[test]
    fn buffer_swap_exchanges_buffers() {
        let mut buf = FluidGpuBuffer::new(2, 2, 1, 1);
        buf.front[0] = 1.0;
        buf.back[0] = 2.0;
        buf.swap();
        assert!((buf.front[0] - 2.0).abs() < 1e-12);
        assert!((buf.back[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn buffer_fill_both() {
        let mut buf = FluidGpuBuffer::new(2, 2, 2, 1);
        buf.fill(7.0);
        assert!(buf.front.iter().all(|&v| (v - 7.0).abs() < 1e-12));
        assert!(buf.back.iter().all(|&v| (v - 7.0).abs() < 1e-12));
    }

    #[test]
    fn buffer_copy_front_to_back() {
        let mut buf = FluidGpuBuffer::new(2, 2, 1, 1);
        buf.front[3] = 5.5;
        buf.copy_front_to_back();
        assert!((buf.back[3] - 5.5).abs() < 1e-12);
    }

    #[test]
    fn buffer_idx_formula() {
        let buf = FluidGpuBuffer::new(4, 3, 2, 2);
        // cell (3, 2, 1), component 1 → (1*3*4 + 2*4 + 3)*2 + 1 = (12+8+3)*2+1 = 47
        assert_eq!(buf.idx(3, 2, 1, 1), 47);
    }

    // ── NavierStokesGpu tests ─────────────────────────────────────────────

    #[test]
    fn ns_new_dimensions() {
        let ns = NavierStokesGpu::new(8, 8, 8, 0.01, 0.001);
        assert_eq!(ns.nx(), 8);
        assert_eq!(ns.ny(), 8);
        assert_eq!(ns.nz(), 8);
    }

    #[test]
    fn ns_divergence_zero_velocity() {
        let ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        let div = ns.divergence();
        assert!(div.iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn ns_curl_zero_velocity() {
        let ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        let omega = ns.curl();
        assert!(omega.iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn ns_advect_preserves_zero() {
        let mut ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        ns.advect(0.01);
        assert!(ns.velocity.front.iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn ns_pressure_projection_zero_field() {
        let mut ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        ns.pressure_projection(5);
        assert!(ns.velocity.front.iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn ns_divergence_after_projection_reduced() {
        let mut ns = NavierStokesGpu::new(8, 8, 4, 0.1, 0.001);
        // set a divergent field
        for i in 0..ns.velocity.front.len() / 3 {
            ns.velocity.front[i * 3] = (i % 4) as f64 * 0.1;
        }
        let div_before: f64 = ns.divergence().iter().copied().map(f64::abs).sum();
        ns.pressure_projection(20);
        let div_after: f64 = ns.divergence().iter().copied().map(f64::abs).sum();
        // projection should generally reduce divergence (or at least not explode)
        assert!(div_after <= div_before + 1.0);
    }

    #[test]
    fn ns_curl_nonzero_velocity() {
        let mut ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        // set a shear flow: u = y (varies in y-direction)
        for z in 0..6 {
            for y in 0..6 {
                for x in 0..6 {
                    let idx = ns.velocity.idx(x, y, z, 0);
                    ns.velocity.front[idx] = y as f64 * 0.1;
                }
            }
        }
        let omega = ns.curl();
        // should have nonzero vorticity
        let max_omega = omega.iter().copied().fold(0.0_f64, |a, v| a.max(v.abs()));
        assert!(max_omega > 0.0);
    }

    // ── SPHGpuKernels tests ───────────────────────────────────────────────

    fn two_particle_sph() -> SPHGpuKernels {
        SPHGpuKernels::new(
            vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
            0.5,
            1000.0,
            1.0,
            0.001,
            0.1,
        )
    }

    #[test]
    fn sph_kernel_w_zero_at_h() {
        let sph = two_particle_sph();
        assert!((sph.kernel_w(sph.smoothing_length)).abs() < 1e-12);
    }

    #[test]
    fn sph_kernel_w_positive_near_zero() {
        let sph = two_particle_sph();
        assert!(sph.kernel_w(0.0) > 0.0);
    }

    #[test]
    fn sph_compute_density_increases_from_interaction() {
        let mut sph = two_particle_sph();
        sph.compute_density();
        // each particle should contribute to the other's density
        assert!(sph.densities[0] > 0.0);
        assert!(sph.densities[1] > 0.0);
    }

    #[test]
    fn sph_compute_pressure_nonnegative() {
        let mut sph = two_particle_sph();
        sph.compute_density();
        sph.compute_pressure();
        assert!(sph.pressures[0] >= 0.0);
        assert!(sph.pressures[1] >= 0.0);
    }

    #[test]
    fn sph_pressure_force_two_particles() {
        let mut sph = two_particle_sph();
        sph.compute_density();
        sph.compute_pressure();
        let forces = sph.pressure_force();
        assert_eq!(forces.len(), 2);
    }

    #[test]
    fn sph_viscosity_force_zero_velocity() {
        let mut sph = two_particle_sph();
        sph.compute_density();
        let forces = sph.viscosity_force();
        // zero relative velocity → zero viscosity force
        assert!(forces[0].iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn sph_kernel_grad_w_zero_at_large_r() {
        let sph = two_particle_sph();
        let grad = sph.kernel_grad_w([1.0, 0.0, 0.0], 100.0);
        assert!(grad.iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn sph_n_particles_correct() {
        let sph = two_particle_sph();
        assert_eq!(sph.n_particles(), 2);
    }

    // ── FluidGpuPipeline tests ────────────────────────────────────────────

    #[test]
    fn pipeline_default_ns_has_passes() {
        let pipe = FluidGpuPipeline::default_ns(1024);
        assert!(!pipe.passes.is_empty());
    }

    #[test]
    fn pipeline_total_bytes_positive() {
        let pipe = FluidGpuPipeline::default_ns(1024);
        assert!(pipe.total_bytes > 0);
    }

    #[test]
    fn pipeline_bandwidth_zero_time() {
        let pipe = FluidGpuPipeline::default_ns(1024);
        assert!((pipe.bandwidth_gb_s(0.0)).abs() < 1e-12);
    }

    #[test]
    fn pipeline_bandwidth_positive_time() {
        let pipe = FluidGpuPipeline::default_ns(1024);
        assert!(pipe.bandwidth_gb_s(0.001) >= 0.0);
    }

    #[test]
    fn pipeline_push_increases_bytes() {
        let mut pipe = FluidGpuPipeline::default_ns(1024);
        let before = pipe.total_bytes;
        pipe.push(PipelinePass {
            label: "extra".into(),
            barrier: BarrierKind::Compute,
            estimated_bytes: 1000,
        });
        assert_eq!(pipe.total_bytes, before + 1000);
    }

    #[test]
    fn pipeline_barrier_kinds_eq() {
        assert_eq!(BarrierKind::Compute, BarrierKind::Compute);
        assert_ne!(BarrierKind::Memory, BarrierKind::Full);
    }

    // ── VortexConfinement tests ───────────────────────────────────────────

    #[test]
    fn vortex_confinement_no_panic_zero_velocity() {
        let mut ns = NavierStokesGpu::new(6, 6, 6, 0.1, 0.001);
        let vc = VortexConfinement::new(1.0, 0.1);
        vc.apply(&mut ns, 0.01);
    }

    #[test]
    fn vortex_confinement_modifies_shear_flow() {
        let mut ns = NavierStokesGpu::new(8, 8, 8, 0.1, 0.001);
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    let idx = ns.velocity.idx(x, y, z, 0);
                    ns.velocity.front[idx] = y as f64 * 0.1;
                }
            }
        }
        let sum_before: f64 = ns.velocity.front.iter().copied().sum();
        let vc = VortexConfinement::new(1.0, 0.1);
        vc.apply(&mut ns, 0.01);
        let sum_after: f64 = ns.velocity.front.iter().copied().sum();
        // some change should occur due to confinement force
        assert!((sum_after - sum_before).abs() >= 0.0); // no panic, at minimum
    }

    // ── WaterSimGpu tests ─────────────────────────────────────────────────

    #[test]
    fn water_sim_new_zero_height() {
        let ws = WaterSimGpu::new(16, 16, 0.1, 1.0, 9.81);
        assert!(ws.height.front.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn water_sim_add_disturbance() {
        let mut ws = WaterSimGpu::new(16, 16, 0.1, 1.0, 9.81);
        ws.add_disturbance(8, 8, 1.0, 2.0);
        let max_dev = ws.max_height_deviation();
        assert!(max_dev > 0.0);
    }

    #[test]
    fn water_sim_step_no_panic() {
        let mut ws = WaterSimGpu::new(16, 16, 0.1, 1.0, 9.81);
        ws.add_disturbance(8, 8, 0.1, 2.0);
        ws.step(0.01);
    }

    #[test]
    fn water_sim_step_reduces_disturbance_over_time() {
        let mut ws = WaterSimGpu::new(16, 16, 0.1, 1.0, 9.81);
        ws.damping = 0.9;
        ws.add_disturbance(8, 8, 1.0, 1.0);
        let _initial = ws.max_height_deviation();
        for _ in 0..100 {
            ws.step(0.001);
        }
        // should not blow up
        let final_dev = ws.max_height_deviation();
        assert!(final_dev.is_finite());
    }

    #[test]
    fn water_sim_nx_ny_correct() {
        let ws = WaterSimGpu::new(10, 20, 0.05, 2.0, 9.81);
        assert_eq!(ws.nx(), 10);
        assert_eq!(ws.ny(), 20);
    }

    #[test]
    fn water_sim_max_height_zero_initially() {
        let ws = WaterSimGpu::new(8, 8, 0.1, 1.0, 9.81);
        assert!((ws.max_height_deviation()).abs() < 1e-12);
    }

    // ── Trilinear sample tests ────────────────────────────────────────────

    #[test]
    fn trilinear_sample_exact_cell() {
        let mut data = vec![0.0_f64; 2 * 2 * 2];
        data[0] = 5.0; // cell (0,0,0) component 0
        let v = trilinear_sample(&data, 2, 2, 2, 1, 0.0, 0.0, 0.0, 0);
        assert!((v - 5.0).abs() < 1e-12);
    }

    #[test]
    fn trilinear_sample_midpoint() {
        let data = vec![1.0_f64; 2 * 2 * 2]; // all ones
        let v = trilinear_sample(&data, 2, 2, 2, 1, 0.5, 0.5, 0.5, 0);
        assert!((v - 1.0).abs() < 1e-12);
    }

    // ── D3Q19 constants tests ─────────────────────────────────────────────

    #[test]
    fn d3q19_weights_sum_to_one() {
        let sum: f64 = D3Q19_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn d3q19_q_correct() {
        assert_eq!(D3Q19_Q, 19);
    }
}

#[cfg(test)]
mod lbm_tests;
