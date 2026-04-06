// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-mock Lattice Boltzmann Method (LBM) compute kernels.
//!
//! Mirrors GPU dispatch layout but executes in pure Rust on the CPU.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── Lattice and collision type enums ─────────────────────────────────────────

/// Supported LBM lattice topologies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeType {
    /// 2-D 9-velocity.
    D2Q9,
    /// 3-D 19-velocity.
    D3Q19,
    /// 3-D 27-velocity.
    D3Q27,
}

/// LBM collision operator type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollisionType {
    /// Single-relaxation-time (Bhatnagar-Gross-Krook).
    Bgk,
    /// Multiple-relaxation-time.
    Mrt,
    /// Two-relaxation-time.
    Trt,
}

// ── LbmKernelConfig ───────────────────────────────────────────────────────────

/// Configuration shared across all LBM kernels.
#[derive(Debug, Clone)]
pub struct LbmKernelConfig {
    /// Lattice topology.
    pub lattice_type: LatticeType,
    /// Collision operator.
    pub collision_type: CollisionType,
    /// Grid size X.
    pub nx: usize,
    /// Grid size Y.
    pub ny: usize,
    /// Grid size Z (1 for 2-D).
    pub nz: usize,
    /// Relaxation time τ.
    pub tau: f64,
    /// Lattice speed of sound squared (1/3 for standard D2Q9).
    pub cs2: f64,
}

impl LbmKernelConfig {
    /// Create a default 2-D 64×64 D2Q9 BGK config.
    pub fn new_d2q9(nx: usize, ny: usize, tau: f64) -> Self {
        Self {
            lattice_type: LatticeType::D2Q9,
            collision_type: CollisionType::Bgk,
            nx,
            ny,
            nz: 1,
            tau,
            cs2: 1.0 / 3.0,
        }
    }

    /// Create a default 3-D D3Q19 config.
    pub fn new_d3q19(nx: usize, ny: usize, nz: usize, tau: f64) -> Self {
        Self {
            lattice_type: LatticeType::D3Q19,
            collision_type: CollisionType::Bgk,
            nx,
            ny,
            nz,
            tau,
            cs2: 1.0 / 3.0,
        }
    }

    /// Kinematic viscosity ν = cs² * (τ - 0.5).
    pub fn kinematic_viscosity(&self) -> f64 {
        self.cs2 * (self.tau - 0.5)
    }

    /// Total number of lattice nodes.
    pub fn total_nodes(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

// ── D2Q9 constants ────────────────────────────────────────────────────────────

/// D2Q9 discrete velocity directions (ex, ey).
pub const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0,-1.0, 0.0, 1.0,-1.0,-1.0, 1.0];
/// D2Q9 ey component.
pub const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0,-1.0, 1.0, 1.0,-1.0,-1.0];
/// D2Q9 equilibrium weights.
pub const D2Q9_W:  [f64; 9] = [
    4.0/9.0,
    1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
];

// ── D3Q19 constants ───────────────────────────────────────────────────────────

/// D3Q19 discrete velocity ex components.
pub const D3Q19_EX: [f64; 19] = [
    0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0,
    1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0,
    0.0, 0.0, 0.0, 0.0,
];
/// D3Q19 ey components.
pub const D3Q19_EY: [f64; 19] = [
    0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0,
    1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0,
    1.0,-1.0, 1.0,-1.0,
];
/// D3Q19 ez components.
pub const D3Q19_EZ: [f64; 19] = [
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,
    1.0, 1.0,-1.0,-1.0,
];
/// D3Q19 weights.
pub const D3Q19_W: [f64; 19] = [
    1.0/3.0,
    1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
];

// ── LbmEquilibriumKernel ──────────────────────────────────────────────────────

/// Kernel for computing Maxwell-Boltzmann equilibrium distributions.
#[derive(Debug, Clone)]
pub struct LbmEquilibriumKernel {
    /// Lattice speed of sound squared.
    pub cs2: f64,
}

impl LbmEquilibriumKernel {
    /// Create a new equilibrium kernel.
    pub fn new(cs2: f64) -> Self {
        Self { cs2 }
    }

    /// Compute equilibrium distribution for D2Q9 at a single node.
    ///
    /// Returns array of 9 f_eq values.
    pub fn compute_feq_d2q9(&self, rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let cs2 = self.cs2;
        let u2 = ux * ux + uy * uy;
        let mut feq = [0.0f64; 9];
        for q in 0..9 {
            let eu = D2Q9_EX[q] * ux + D2Q9_EY[q] * uy;
            feq[q] = D2Q9_W[q] * rho * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2)
                - u2 / (2.0 * cs2));
        }
        feq
    }

    /// Compute equilibrium distribution for D3Q19 at a single node.
    ///
    /// `u` — velocity vector \[ux, uy, uz\].
    pub fn compute_feq_d3q19(&self, rho: f64, u: [f64; 3]) -> [f64; 19] {
        let cs2 = self.cs2;
        let u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        let mut feq = [0.0f64; 19];
        for q in 0..19 {
            let eu = D3Q19_EX[q]*u[0] + D3Q19_EY[q]*u[1] + D3Q19_EZ[q]*u[2];
            feq[q] = D3Q19_W[q] * rho * (1.0 + eu/cs2 + eu*eu/(2.0*cs2*cs2) - u2/(2.0*cs2));
        }
        feq
    }
}

// ── LbmCollisionKernel ────────────────────────────────────────────────────────

/// Kernel for performing LBM collision steps.
#[derive(Debug, Clone)]
pub struct LbmCollisionKernel {
    /// Relaxation time τ.
    pub tau: f64,
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmCollisionKernel {
    /// Create a new collision kernel.
    pub fn new(tau: f64, cs2: f64) -> Self {
        Self { tau, cs2 }
    }

    /// BGK collision: f_post\[q\] = f\[q\] - (f\[q\] - feq\[q\]) / τ.
    pub fn bgk_collision(&self, f: &[f64; 9], feq: &[f64; 9]) -> [f64; 9] {
        let mut f_post = [0.0f64; 9];
        let inv_tau = 1.0 / self.tau;
        for q in 0..9 {
            f_post[q] = f[q] - inv_tau * (f[q] - feq[q]);
        }
        f_post
    }

    /// MRT collision: f_post = f - M^{-1} * S * (m - m_eq).
    ///
    /// For demo, uses a simplified diagonal relaxation.
    /// `m_s` — relaxation rates in moment space (9 values).
    pub fn mrt_collision(
        &self,
        f: &[f64; 9],
        feq: &[f64; 9],
        m_s: &[f64; 9],
    ) -> [f64; 9] {
        // Simplified: treat as independent mode relaxation
        let mut f_post = [0.0f64; 9];
        for q in 0..9 {
            f_post[q] = f[q] - m_s[q] * (f[q] - feq[q]);
        }
        f_post
    }

    /// TRT collision using symmetric (+) and antisymmetric (−) relaxation.
    ///
    /// `lambda_plus`  — relaxation for symmetric part.
    /// `lambda_minus` — relaxation for antisymmetric part.
    pub fn trt_collision(
        &self,
        f: &[f64; 9],
        feq: &[f64; 9],
        lambda_plus: f64,
        lambda_minus: f64,
    ) -> [f64; 9] {
        // Opposite directions: 0↔0, 1↔3, 2↔4, 5↔7, 6↔8
        let opp: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        let mut f_post = [0.0f64; 9];
        for q in 0..9 {
            let q_opp = opp[q];
            let f_sym  = 0.5 * ((f[q] - feq[q]) + (f[q_opp] - feq[q_opp]));
            let f_asym = 0.5 * ((f[q] - feq[q]) - (f[q_opp] - feq[q_opp]));
            f_post[q] = f[q] - lambda_plus * f_sym - lambda_minus * f_asym;
        }
        f_post
    }

    /// Apply BGK collision to an entire D2Q9 field (size nx * ny * 9).
    ///
    /// `field` — flat array \[node_idx * 9 + q\].
    pub fn bgk_collision_field(
        &self,
        field: &mut Vec<f64>,
        feq_field: &[f64],
        n_nodes: usize,
    ) {
        let inv_tau = 1.0 / self.tau;
        for node in 0..n_nodes {
            for q in 0..9 {
                let idx = node * 9 + q;
                field[idx] -= inv_tau * (field[idx] - feq_field[idx]);
            }
        }
    }
}

// ── LbmStreamingKernel ────────────────────────────────────────────────────────

/// Kernel for streaming (propagation) steps.
#[derive(Debug, Clone)]
pub struct LbmStreamingKernel;

impl LbmStreamingKernel {
    /// Create a new streaming kernel.
    pub fn new() -> Self {
        Self
    }

    /// Stream D2Q9 field with periodic boundary conditions.
    ///
    /// `f` — input field \[nx * ny * 9\].
    /// Returns new field after streaming.
    pub fn stream_d2q9(&self, f: &[f64], nx: usize, ny: usize) -> Vec<f64> {
        let n = nx * ny * 9;
        let mut f_new = vec![0.0f64; n];
        let ex = D2Q9_EX.map(|v| v as i64);
        let ey = D2Q9_EY.map(|v| v as i64);
        for iy in 0..ny {
            for ix in 0..nx {
                for q in 0..9 {
                    let src_x = (ix as i64 - ex[q]).rem_euclid(nx as i64) as usize;
                    let src_y = (iy as i64 - ey[q]).rem_euclid(ny as i64) as usize;
                    let dst = (iy * nx + ix) * 9 + q;
                    let src = (src_y * nx + src_x) * 9 + q;
                    f_new[dst] = f[src];
                }
            }
        }
        f_new
    }

    /// Stream D3Q19 field with periodic boundary conditions.
    ///
    /// `f` — input field \[nx * ny * nz * 19\].
    pub fn stream_d3q19(&self, f: &[f64], nx: usize, ny: usize, nz: usize) -> Vec<f64> {
        let n = nx * ny * nz * 19;
        let mut f_new = vec![0.0f64; n];
        let ex = D3Q19_EX.map(|v| v as i64);
        let ey = D3Q19_EY.map(|v| v as i64);
        let ez = D3Q19_EZ.map(|v| v as i64);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    for q in 0..19 {
                        let sx = (ix as i64 - ex[q]).rem_euclid(nx as i64) as usize;
                        let sy = (iy as i64 - ey[q]).rem_euclid(ny as i64) as usize;
                        let sz = (iz as i64 - ez[q]).rem_euclid(nz as i64) as usize;
                        let dst = ((iz * ny + iy) * nx + ix) * 19 + q;
                        let src = ((sz * ny + sy) * nx + sx) * 19 + q;
                        f_new[dst] = f[src];
                    }
                }
            }
        }
        f_new
    }
}

impl Default for LbmStreamingKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── LbmBoundaryKernel ─────────────────────────────────────────────────────────

/// Direction for Zou-He boundary conditions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryDirection {
    /// West (x=0) inlet.
    West,
    /// East (x=nx-1) outlet.
    East,
    /// South (y=0) wall.
    South,
    /// North (y=ny-1) wall.
    North,
}

/// Kernel for LBM boundary conditions.
#[derive(Debug, Clone)]
pub struct LbmBoundaryKernel {
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmBoundaryKernel {
    /// Create a new boundary kernel.
    pub fn new(cs2: f64) -> Self {
        Self { cs2 }
    }

    /// Zou-He velocity boundary: enforce given velocity at a boundary node.
    ///
    /// `f`        — 9-component distribution at the node (modified in-place).
    /// `u_wall`   — prescribed velocity \[ux, uy\].
    /// `direction`— boundary face direction.
    pub fn zou_he_velocity(
        &self,
        f: &mut [f64; 9],
        u_wall: [f64; 2],
        direction: BoundaryDirection,
    ) {
        match direction {
            BoundaryDirection::West => {
                // rho from known distributions
                let rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - u_wall[0]);
                let ux = u_wall[0];
                let uy = u_wall[1];
                f[1] = f[3] + 2.0/3.0 * rho * ux;
                f[5] = f[7] - 0.5*(f[2]-f[4]) + 0.5*rho*uy + 1.0/6.0*rho*ux;
                f[8] = f[6] + 0.5*(f[2]-f[4]) - 0.5*rho*uy + 1.0/6.0*rho*ux;
            }
            BoundaryDirection::East => {
                let rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + u_wall[0]);
                let ux = u_wall[0];
                let uy = u_wall[1];
                f[3] = f[1] - 2.0/3.0 * rho * ux;
                f[7] = f[5] + 0.5*(f[2]-f[4]) - 0.5*rho*uy - 1.0/6.0*rho*ux;
                f[6] = f[8] - 0.5*(f[2]-f[4]) + 0.5*rho*uy - 1.0/6.0*rho*ux;
            }
            BoundaryDirection::South => {
                let rho = (f[0] + f[1] + f[3] + 2.0*(f[4]+f[7]+f[8])) / (1.0 - u_wall[1]);
                let ux = u_wall[0];
                let uy = u_wall[1];
                f[2] = f[4] + 2.0/3.0 * rho * uy;
                f[5] = f[7] - 0.5*(f[1]-f[3]) + 0.5*rho*ux + 1.0/6.0*rho*uy;
                f[6] = f[8] + 0.5*(f[1]-f[3]) - 0.5*rho*ux + 1.0/6.0*rho*uy;
            }
            BoundaryDirection::North => {
                let rho = (f[0]+f[1]+f[3]+2.0*(f[2]+f[5]+f[6])) / (1.0+u_wall[1]);
                let ux = u_wall[0];
                let uy = u_wall[1];
                f[4] = f[2] - 2.0/3.0 * rho * uy;
                f[7] = f[5] + 0.5*(f[1]-f[3]) - 0.5*rho*ux - 1.0/6.0*rho*uy;
                f[8] = f[6] - 0.5*(f[1]-f[3]) + 0.5*rho*ux - 1.0/6.0*rho*uy;
            }
        }
    }

    /// Bounce-back boundary: reverse populations at solid nodes.
    ///
    /// `mask` — boolean slice (true = solid node), length n_nodes.
    /// `f`    — flat field \[n_nodes * 9\].
    pub fn bounce_back(&self, f: &mut Vec<f64>, mask: &[bool], n_nodes: usize) {
        let opp: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        for node in 0..n_nodes {
            if mask[node] {
                let base = node * 9;
                let tmp = [
                    f[base], f[base+1], f[base+2], f[base+3], f[base+4],
                    f[base+5], f[base+6], f[base+7], f[base+8],
                ];
                for q in 0..9 {
                    f[base + q] = tmp[opp[q]];
                }
            }
        }
    }
}

// ── LbmMomentKernel ───────────────────────────────────────────────────────────

/// Kernel for computing macroscopic moments from distributions.
#[derive(Debug, Clone)]
pub struct LbmMomentKernel {
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmMomentKernel {
    /// Create a new moment kernel.
    pub fn new(cs2: f64) -> Self {
        Self { cs2 }
    }

    /// Compute density ρ and velocity u from D2Q9 distribution at a node.
    pub fn compute_rho_u_d2q9(&self, f: &[f64; 9]) -> (f64, [f64; 3]) {
        let rho: f64 = f.iter().sum();
        let ux: f64 = f.iter().zip(D2Q9_EX.iter()).map(|(fi, ei)| fi * ei).sum();
        let uy: f64 = f.iter().zip(D2Q9_EY.iter()).map(|(fi, ei)| fi * ei).sum();
        if rho.abs() < 1e-15 {
            return (rho, [0.0; 3]);
        }
        (rho, [ux / rho, uy / rho, 0.0])
    }

    /// Compute density and velocity from D3Q19 distribution.
    pub fn compute_rho_u_d3q19(&self, f: &[f64; 19]) -> (f64, [f64; 3]) {
        let rho: f64 = f.iter().sum();
        let ux: f64 = f.iter().zip(D3Q19_EX.iter()).map(|(fi, ei)| fi * ei).sum();
        let uy: f64 = f.iter().zip(D3Q19_EY.iter()).map(|(fi, ei)| fi * ei).sum();
        let uz: f64 = f.iter().zip(D3Q19_EZ.iter()).map(|(fi, ei)| fi * ei).sum();
        if rho.abs() < 1e-15 {
            return (rho, [0.0; 3]);
        }
        (rho, [ux / rho, uy / rho, uz / rho])
    }

    /// Compute pressure from density: p = ρ * cs².
    pub fn compute_pressure(&self, rho: f64) -> f64 {
        rho * self.cs2
    }

    /// Compute kinetic energy density: ke = 0.5 * ρ * |u|².
    pub fn kinetic_energy(&self, rho: f64, u: [f64; 3]) -> f64 {
        0.5 * rho * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])
    }
}

// ── LbmTurbulenceKernel ───────────────────────────────────────────────────────

/// Kernel for LBM turbulence sub-grid scale models.
#[derive(Debug, Clone)]
pub struct LbmTurbulenceKernel {
    /// Smagorinsky constant.
    pub c_s: f64,
    /// WALE model constant.
    pub c_w: f64,
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmTurbulenceKernel {
    /// Create a new turbulence kernel.
    pub fn new(c_s: f64, c_w: f64, cs2: f64) -> Self {
        Self { c_s, c_w, cs2 }
    }

    /// Smagorinsky effective relaxation time.
    ///
    /// `tau_0` — molecular τ. `s_bar` — strain rate magnitude |S̄|.
    /// `delta` — filter width (lattice spacing = 1 typically).
    pub fn smagorinsky_tau(&self, tau_0: f64, s_bar: f64, delta: f64) -> f64 {
        let nu_0 = self.cs2 * (tau_0 - 0.5);
        let nu_sgs = (self.c_s * delta).powi(2) * s_bar;
        let nu_eff = nu_0 + nu_sgs;
        0.5 + nu_eff / self.cs2
    }

    /// WALE sub-grid scale viscosity.
    ///
    /// `g_ij` — velocity gradient tensor (3×3, row-major flat).
    pub fn wale_tau(&self, g_ij: &[f64; 9], tau_0: f64, delta: f64) -> f64 {
        // S_ij = (g_ij + g_ji) / 2
        // S^d_ij = 0.5*(g_ij^2 + g_ji^2) - delta_ij/3 * g_kk^2
        let s_d_sq: f64 = {
            let mut sd = [0.0f64; 9];
            for i in 0..3 {
                for j in 0..3 {
                    let g2_ij: f64 = (0..3).map(|k| g_ij[i*3+k] * g_ij[k*3+j]).sum();
                    let g2_ji: f64 = (0..3).map(|k| g_ij[j*3+k] * g_ij[k*3+i]).sum();
                    sd[i*3+j] = 0.5 * (g2_ij + g2_ji);
                }
            }
            let trace = sd[0] + sd[4] + sd[8];
            for i in 0..3 { sd[i*3+i] -= trace / 3.0; }
            sd.iter().map(|x| x*x).sum::<f64>()
        };
        // s_ij s_ij
        let s_sq: f64 = {
            let mut s = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let sij = 0.5 * (g_ij[i*3+j] + g_ij[j*3+i]);
                    s += sij * sij;
                }
            }
            s
        };
        let denom = s_sq.powf(2.5) + s_d_sq.powf(1.25);
        let nu_sgs = if denom > 1e-15 {
            (self.c_w * delta).powi(2) * s_d_sq.powf(1.5) / denom
        } else {
            0.0
        };
        let nu_0 = self.cs2 * (tau_0 - 0.5);
        0.5 + (nu_0 + nu_sgs) / self.cs2
    }

    /// Convert effective τ to relaxation rate ω = 1/τ.
    pub fn omega_from_tau(&self, tau: f64) -> f64 {
        1.0 / tau
    }
}

// ── LbmMultiphaseKernel ───────────────────────────────────────────────────────

/// Kernel for multiphase LBM (Shan-Chen and free-energy approaches).
#[derive(Debug, Clone)]
pub struct LbmMultiphaseKernel {
    /// Shan-Chen coupling parameter G.
    pub g_sc: f64,
    /// Free-energy interface parameter κ.
    pub kappa: f64,
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmMultiphaseKernel {
    /// Create a new multiphase kernel.
    pub fn new(g_sc: f64, kappa: f64, cs2: f64) -> Self {
        Self { g_sc, kappa, cs2 }
    }

    /// Shan-Chen pseudo-potential ψ(ρ) = ρ₀ * exp(-ρ₀/ρ).
    pub fn psi(&self, rho: f64, rho0: f64) -> f64 {
        rho0 * (-rho0 / rho.max(1e-15)).exp()
    }

    /// Shan-Chen body force at a node using nearest-neighbor ψ values.
    ///
    /// `psi_neighbors` — ψ values at 8 neighbors (D2Q9 directions 1..8).
    /// `psi_self`      — ψ at the current node.
    /// Returns force \[fx, fy\].
    pub fn shan_chen_force(
        &self,
        psi_self: f64,
        psi_neighbors: &[f64; 8],
        rho_self: f64,
    ) -> [f64; 2] {
        let mut fx = 0.0;
        let mut fy = 0.0;
        for q in 0..8 {
            let ex = D2Q9_EX[q + 1];
            let ey = D2Q9_EY[q + 1];
            let w = D2Q9_W[q + 1];
            fx -= self.g_sc * psi_self * w * psi_neighbors[q] * ex;
            fy -= self.g_sc * psi_self * w * psi_neighbors[q] * ey;
        }
        let _ = rho_self;
        [fx, fy]
    }

    /// Free-energy pressure for a binary mixture.
    ///
    /// p = ρ * cs² - κ * ρ * ∇²ρ + κ/2 * |∇ρ|²
    /// (simplified form without gradient terms).
    pub fn free_energy_pressure(&self, rho: f64, phi: f64) -> f64 {
        // Landau free energy bulk term
        let a = -0.125; // negative for phase separation
        let b = 0.125;
        rho * self.cs2 + a * phi.powi(2) + b * phi.powi(4)
    }

    /// Surface tension estimate σ = κ * ∫|∂φ/∂x|² dx ≈ κ * amplitude² / width.
    pub fn surface_tension_estimate(&self, amplitude: f64, interface_width: f64) -> f64 {
        self.kappa * amplitude * amplitude / interface_width.max(1e-15)
    }
}

// ── LbmAcousticsKernel ────────────────────────────────────────────────────────

/// Kernel for acoustic extraction from LBM fields.
#[derive(Debug, Clone)]
pub struct LbmAcousticsKernel {
    /// Reference density ρ₀.
    pub rho_ref: f64,
    /// Lattice cs².
    pub cs2: f64,
}

impl LbmAcousticsKernel {
    /// Create a new acoustics kernel.
    pub fn new(rho_ref: f64, cs2: f64) -> Self {
        Self { rho_ref, cs2 }
    }

    /// Extract acoustic pressure fluctuation p' = (ρ - ρ₀) * cs².
    pub fn extract_acoustic_pressure(&self, rho: f64) -> f64 {
        (rho - self.rho_ref) * self.cs2
    }

    /// Compute room impulse response (RIR) at a receiver point from a
    /// time-series of pressure samples.
    ///
    /// Applies a simple reflection model: direct + single reflection.
    /// `p_samples` — time series of pressure at source.
    /// `delay`     — integer delay in samples for reflection.
    /// `gain`      — reflection gain coefficient.
    pub fn compute_rir_receiver(
        &self,
        p_samples: &[f64],
        delay: usize,
        gain: f64,
    ) -> Vec<f64> {
        let n = p_samples.len();
        let mut rir = vec![0.0f64; n];
        for i in 0..n {
            rir[i] += p_samples[i];
            if i >= delay {
                rir[i] += gain * p_samples[i - delay];
            }
        }
        rir
    }

    /// Estimate sound pressure level (SPL) in dB from RMS pressure.
    ///
    /// `p_ref` — reference pressure (typically 20 μPa in air).
    pub fn spl_db(&self, p_rms: f64, p_ref: f64) -> f64 {
        if p_rms < 1e-15 || p_ref < 1e-15 { return f64::NEG_INFINITY; }
        20.0 * (p_rms / p_ref).log10()
    }
}

// ── LbmVorticityKernel ────────────────────────────────────────────────────────

/// Kernel for computing vorticity and flow topology criteria.
#[derive(Debug, Clone)]
pub struct LbmVorticityKernel {
    /// Lattice spacing dx.
    pub dx: f64,
}

impl LbmVorticityKernel {
    /// Create a new vorticity kernel.
    pub fn new(dx: f64) -> Self {
        Self { dx }
    }

    /// 2-D vorticity ω_z = ∂v/∂x - ∂u/∂y using central differences.
    ///
    /// `u`, `v` — velocity component fields of size nx×ny.
    /// Returns vorticity field of size nx×ny.
    pub fn compute_vorticity_2d(
        &self,
        u: &[f64],
        v: &[f64],
        nx: usize,
        ny: usize,
    ) -> Vec<f64> {
        let mut omega = vec![0.0f64; nx * ny];
        let inv_2dx = 1.0 / (2.0 * self.dx);
        for iy in 1..ny-1 {
            for ix in 1..nx-1 {
                let idx = iy * nx + ix;
                let dvdx = (v[iy*nx + (ix+1)] - v[iy*nx + (ix-1)]) * inv_2dx;
                let dudy = (u[(iy+1)*nx + ix] - u[(iy-1)*nx + ix]) * inv_2dx;
                omega[idx] = dvdx - dudy;
            }
        }
        omega
    }

    /// Q-criterion: Q = -0.5 * (S_ij S_ij - Ω_ij Ω_ij).
    ///
    /// `vel_grad` — velocity gradient tensor 3×3, row-major \[du/dx, du/dy, du/dz, ...\].
    pub fn q_criterion_3d(&self, vel_grad: &[f64; 9]) -> f64 {
        let mut s_sq = 0.0;
        let mut omega_sq = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let s_ij = 0.5 * (vel_grad[i*3+j] + vel_grad[j*3+i]);
                let o_ij = 0.5 * (vel_grad[i*3+j] - vel_grad[j*3+i]);
                s_sq += s_ij * s_ij;
                omega_sq += o_ij * o_ij;
            }
        }
        0.5 * (omega_sq - s_sq)
    }

    /// Lambda-2 criterion for vortex identification.
    ///
    /// Returns the second eigenvalue of (S² + Ω²) where S is symmetric and
    /// Ω is antisymmetric part of the velocity gradient.
    /// A simplified trace-based approximation is used here.
    pub fn lambda2_criterion(&self, vel_grad: &[f64; 9]) -> f64 {
        let mut m = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let s_ij = 0.5 * (vel_grad[i*3+j] + vel_grad[j*3+i]);
                let o_ij = 0.5 * (vel_grad[i*3+j] - vel_grad[j*3+i]);
                m[i*3+j] = s_ij * s_ij + o_ij * o_ij;
            }
        }
        // Return middle eigenvalue approximation: use trace/3
        (m[0] + m[4] + m[8]) / 3.0 - (m[0] - m[8]).abs() * 0.5
    }

    /// Helicity density h = u · ω.
    pub fn helicity_density(&self, u: [f64; 3], omega: [f64; 3]) -> f64 {
        u[0]*omega[0] + u[1]*omega[1] + u[2]*omega[2]
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod lbm_tests {
    use super::*;

    #[test]
    fn test_config_kinematic_viscosity() {
        let cfg = LbmKernelConfig::new_d2q9(64, 64, 1.0);
        let nu = cfg.kinematic_viscosity();
        assert!((nu - 1.0/6.0).abs() < 1e-10);
    }

    #[test]
    fn test_config_total_nodes() {
        let cfg = LbmKernelConfig::new_d2q9(8, 8, 1.0);
        assert_eq!(cfg.total_nodes(), 64);
    }

    #[test]
    fn test_d2q9_weights_sum_to_one() {
        let sum: f64 = D2Q9_W.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_d3q19_weights_sum_to_one() {
        let sum: f64 = D3Q19_W.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_feq_d2q9_sums_to_rho() {
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.2, 0.05, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.2).abs() < 1e-10);
    }

    #[test]
    fn test_feq_d3q19_sums_to_rho() {
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d3q19(1.0, [0.01, 0.0, 0.0]);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_feq_d2q9_zero_velocity() {
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.0, 0.0, 0.0);
        // At zero velocity, feq[q] = w[q] * rho
        for q in 0..9 {
            assert!((feq[q] - D2Q9_W[q]).abs() < 1e-12, "q={q}");
        }
    }

    #[test]
    fn test_bgk_collision_at_equilibrium() {
        let coll = LbmCollisionKernel::new(1.0, 1.0/3.0);
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.0, 0.0, 0.0);
        let f_post = coll.bgk_collision(&feq, &feq);
        // At equilibrium, no change
        for q in 0..9 {
            assert!((f_post[q] - feq[q]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_mrt_collision_at_equilibrium() {
        let coll = LbmCollisionKernel::new(1.0, 1.0/3.0);
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.0, 0.0, 0.0);
        let m_s = [1.0/1.0; 9];
        let f_post = coll.mrt_collision(&feq, &feq, &m_s);
        for q in 0..9 {
            assert!((f_post[q] - feq[q]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_trt_collision_at_equilibrium() {
        let coll = LbmCollisionKernel::new(1.0, 1.0/3.0);
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.0, 0.0, 0.0);
        let f_post = coll.trt_collision(&feq, &feq, 1.0, 1.0);
        for q in 0..9 {
            assert!((f_post[q] - feq[q]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_stream_d2q9_conserves_mass() {
        let nx = 4; let ny = 4;
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let mut f = Vec::new();
        for _ in 0..nx*ny {
            let feq = eq.compute_feq_d2q9(1.0, 0.05, 0.0);
            f.extend_from_slice(&feq);
        }
        let total_before: f64 = f.iter().sum();
        let sk = LbmStreamingKernel::new();
        let f_new = sk.stream_d2q9(&f, nx, ny);
        let total_after: f64 = f_new.iter().sum();
        assert!((total_before - total_after).abs() < 1e-10);
    }

    #[test]
    fn test_stream_d3q19_conserves_mass() {
        let nx = 3; let ny = 3; let nz = 3;
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let mut f = Vec::new();
        for _ in 0..nx*ny*nz {
            let feq = eq.compute_feq_d3q19(1.0, [0.01, 0.0, 0.0]);
            f.extend_from_slice(&feq);
        }
        let total_before: f64 = f.iter().sum();
        let sk = LbmStreamingKernel::new();
        let f_new = sk.stream_d3q19(&f, nx, ny, nz);
        let total_after: f64 = f_new.iter().sum();
        assert!((total_before - total_after).abs() < 1e-10);
    }

    #[test]
    fn test_compute_rho_u_rest() {
        let mk = LbmMomentKernel::new(1.0/3.0);
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.2, 0.0, 0.0);
        let (rho, u) = mk.compute_rho_u_d2q9(&feq);
        assert!((rho - 1.2).abs() < 1e-10);
        assert!(u[0].abs() < 1e-12);
        assert!(u[1].abs() < 1e-12);
    }

    #[test]
    fn test_compute_pressure() {
        let mk = LbmMomentKernel::new(1.0/3.0);
        let p = mk.compute_pressure(1.5);
        assert!((p - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_bounce_back_reversal() {
        let bk = LbmBoundaryKernel::new(1.0/3.0);
        let eq = LbmEquilibriumKernel::new(1.0/3.0);
        let feq = eq.compute_feq_d2q9(1.0, 0.05, 0.0);
        let mut f = feq.to_vec();
        let mask = vec![true];
        bk.bounce_back(&mut f, &mask, 1);
        // q=1 and q=3 are swapped, etc.
        assert!((f[1] - feq[3]).abs() < 1e-12);
        assert!((f[3] - feq[1]).abs() < 1e-12);
    }

    #[test]
    fn test_smagorinsky_tau() {
        let tk = LbmTurbulenceKernel::new(0.1, 0.5, 1.0/3.0);
        let tau_eff = tk.smagorinsky_tau(1.0, 0.1, 1.0);
        assert!(tau_eff > 0.5);
    }

    #[test]
    fn test_omega_from_tau() {
        let tk = LbmTurbulenceKernel::new(0.1, 0.5, 1.0/3.0);
        let omega = tk.omega_from_tau(1.0);
        assert!((omega - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shan_chen_psi() {
        let mk = LbmMultiphaseKernel::new(-1.0, 0.01, 1.0/3.0);
        let psi = mk.psi(1.0, 1.0);
        assert!(psi > 0.0);
    }

    #[test]
    fn test_free_energy_pressure() {
        let mk = LbmMultiphaseKernel::new(-1.0, 0.01, 1.0/3.0);
        let p = mk.free_energy_pressure(1.0, 0.5);
        assert!(p.is_finite());
    }

    #[test]
    fn test_extract_acoustic_pressure() {
        let ak = LbmAcousticsKernel::new(1.0, 1.0/3.0);
        let p = ak.extract_acoustic_pressure(1.01);
        assert!((p - 0.01/3.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_rir_receiver() {
        let ak = LbmAcousticsKernel::new(1.0, 1.0/3.0);
        let p = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let rir = ak.compute_rir_receiver(&p, 2, 0.5);
        assert!((rir[0] - 1.0).abs() < 1e-12);
        assert!((rir[2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_spl_db() {
        let ak = LbmAcousticsKernel::new(1.0, 1.0/3.0);
        let spl = ak.spl_db(20e-3, 20e-6);
        assert!((spl - 60.0).abs() < 1e-6);
    }

    #[test]
    fn test_vorticity_2d_uniform_flow() {
        let vk = LbmVorticityKernel::new(1.0);
        let nx = 5; let ny = 5;
        let u = vec![1.0f64; nx * ny];
        let v = vec![0.0f64; nx * ny];
        let omega = vk.compute_vorticity_2d(&u, &v, nx, ny);
        // Uniform flow: vorticity = 0 at interior
        for iy in 1..ny-1 {
            for ix in 1..nx-1 {
                assert!(omega[iy*nx+ix].abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_q_criterion_irrotational() {
        let vk = LbmVorticityKernel::new(1.0);
        // Pure shear: g_ij = [[1,0,0],[0,-1,0],[0,0,0]]
        let g = [1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,0.0];
        let q = vk.q_criterion_3d(&g);
        // S²=2, Ω²=0 => Q = -1
        assert!(q < 0.0);
    }

    #[test]
    fn test_helicity_density() {
        let vk = LbmVorticityKernel::new(1.0);
        let h = vk.helicity_density([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((h - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lattice_type_variants() {
        assert_ne!(LatticeType::D2Q9, LatticeType::D3Q19);
        assert_ne!(LatticeType::D3Q19, LatticeType::D3Q27);
    }

    #[test]
    fn test_collision_type_variants() {
        assert_ne!(CollisionType::Bgk, CollisionType::Mrt);
        assert_ne!(CollisionType::Mrt, CollisionType::Trt);
    }

    #[test]
    fn test_kinetic_energy() {
        let mk = LbmMomentKernel::new(1.0/3.0);
        let ke = mk.kinetic_energy(1.0, [1.0, 0.0, 0.0]);
        assert!((ke - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_wale_tau_zero_gradient() {
        let tk = LbmTurbulenceKernel::new(0.1, 0.5, 1.0/3.0);
        let g = [0.0f64; 9];
        let tau_eff = tk.wale_tau(&g, 1.0, 1.0);
        // Zero gradient -> nu_sgs = 0 -> tau_eff = tau_0
        assert!((tau_eff - 1.0).abs() < 1e-10);
    }
}
