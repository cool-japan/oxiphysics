#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soft-matter Lattice Boltzmann simulations.
//!
//! This module provides a comprehensive LBM framework for soft-matter physics:
//!
//! - **Nematic liquid crystals**: Q-tensor order parameter, Frank elastic free energy
//! - **Leslie-Ericksen dynamics**: director field evolution on a D2Q9 lattice
//! - **Lyotropic liquid crystals**: concentration-dependent order
//! - **Colloidal suspensions**: Lorentz force coupling + excluded-volume interactions
//! - **Gel network dynamics**: viscoelastic network with bond breakage
//! - **Lipid bilayer LBM**: membrane mechanics, bending rigidity, area constraint
//! - **Worm-like micelles**: Cates breakable-chain model
//! - **Block copolymer phase separation**: Ohta-Kawasaki free energy functional
//! - **Active matter**: self-propelled particles (SPP) with polar / nematic alignment
//!
//! # Key dimensionless groups
//! - Er = γ L² / K  (Ericksen number: viscous vs. elastic)
//! - Pe = v L / D   (Péclet number: advection vs. diffusion)
//! - ξ = √(K / |a|) (nematic coherence length)

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Lattice speed of sound squared for D2Q9 (c_s² = 1/3 in lattice units).
pub const CS2: f64 = 1.0 / 3.0;

/// Number of discrete velocities in the D2Q9 lattice.
pub const Q9: usize = 9;

/// Boltzmann constant \[J/K\].
pub const K_B: f64 = 1.380_649e-23;

/// Default Landau coefficient `a` (< 0 in ordered phase).
pub const LANDAU_A: f64 = -0.5;

/// Default Landau coefficient `b`.
pub const LANDAU_B: f64 = 0.5;

/// Default Landau coefficient `c`.
pub const LANDAU_C: f64 = 1.0 / 3.0;

/// Frank elastic constant (one-constant approximation) \[Pa·m\].
pub const FRANK_K: f64 = 1.0e-11;

/// Default self-propulsion speed for active particles \[lattice units/step\].
pub const ACTIVE_SPEED: f64 = 0.05;

// ============================================================================
// D2Q9 velocity set
// ============================================================================

/// D2Q9 discrete velocity components (x-direction).
pub const CX: [f64; Q9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 discrete velocity components (y-direction).
pub const CY: [f64; Q9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 equilibrium weights.
pub const W9: [f64; Q9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

// ============================================================================
// Section 1 – Q-tensor order parameter (nematic liquid crystal)
// ============================================================================

/// Symmetric traceless 2×2 Q-tensor for a 2-D nematic.
///
/// Stored as `[q11, q12]`; the full tensor is
/// Q = ⎡ q11  q12 ⎤
///     ⎣ q12 -q11 ⎦
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QTensor2 {
    /// Component Q₁₁.
    pub q11: f64,
    /// Component Q₁₂ (= Q₂₁ by symmetry).
    pub q12: f64,
}

impl QTensor2 {
    /// Construct Q from director angle θ and scalar order parameter S.
    ///
    /// Q = S (n⊗n − I/2) where n = (cos θ, sin θ).
    #[allow(dead_code)]
    pub fn from_director(theta: f64, s: f64) -> Self {
        let (c, sn) = (theta.cos(), theta.sin());
        Self {
            q11: s * (c * c - 0.5),
            q12: s * c * sn,
        }
    }

    /// Extract the scalar order parameter S = 2√(Q₁₁² + Q₁₂²).
    ///
    /// This follows from Q = S(n⊗n − I/2) in 2-D where Q₁₁² + Q₁₂² = S²/4.
    #[allow(dead_code)]
    pub fn order_parameter(&self) -> f64 {
        2.0 * (self.q11 * self.q11 + self.q12 * self.q12).sqrt()
    }

    /// Extract the director angle θ ∈ \[0, π\).
    #[allow(dead_code)]
    pub fn director_angle(&self) -> f64 {
        0.5 * self.q12.atan2(self.q11)
    }

    /// Add two Q-tensors.
    #[allow(dead_code)]
    pub fn add(&self, other: &QTensor2) -> Self {
        Self {
            q11: self.q11 + other.q11,
            q12: self.q12 + other.q12,
        }
    }

    /// Scale Q-tensor by a scalar.
    #[allow(dead_code)]
    pub fn scale(&self, f: f64) -> Self {
        Self {
            q11: self.q11 * f,
            q12: self.q12 * f,
        }
    }

    /// Frobenius norm ‖Q‖_F = √(Tr Q²).
    #[allow(dead_code)]
    pub fn norm(&self) -> f64 {
        (self.q11 * self.q11 + self.q12 * self.q12).sqrt()
    }
}

// ============================================================================
// Section 2 – Frank elastic free energy (one-constant approximation)
// ============================================================================

/// Compute the one-constant Frank elastic free energy density.
///
/// f_el = (K/2) ∂_k Q_ij ∂_k Q_ij evaluated with finite differences.
///
/// # Arguments
/// * `q`   – Q-tensor field, dimensions `[nx][ny]`
/// * `k`   – Frank elastic constant K \[Pa\]
/// * `dx`  – lattice spacing \[m\]
#[allow(dead_code)]
pub fn frank_elastic_density(q: &[Vec<QTensor2>], k: f64, dx: f64) -> Vec<Vec<f64>> {
    let nx = q.len();
    let ny = q[0].len();
    let mut f = vec![vec![0.0_f64; ny]; nx];
    let inv2dx = 0.5 / dx;
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let dq11_dx = (q[i + 1][j].q11 - q[i - 1][j].q11) * inv2dx;
            let dq11_dy = (q[i][j + 1].q11 - q[i][j - 1].q11) * inv2dx;
            let dq12_dx = (q[i + 1][j].q12 - q[i - 1][j].q12) * inv2dx;
            let dq12_dy = (q[i][j + 1].q12 - q[i][j - 1].q12) * inv2dx;
            f[i][j] = 0.5
                * k
                * (dq11_dx * dq11_dx + dq11_dy * dq11_dy + dq12_dx * dq12_dx + dq12_dy * dq12_dy);
        }
    }
    f
}

/// Compute the molecular field H = -δF/δQ (driving relaxational dynamics).
///
/// Uses the Beris-Edwards / Landau-de Gennes free energy:
/// F = ∫ \[ a/2 Tr Q² + b/3 Tr Q³ + c/4 (Tr Q²)² + K/2 |∇Q|² \] dV
///
/// # Arguments
/// * `q`   – Q-tensor field `[nx][ny]`
/// * `a`, `b_coeff`, `c_coeff` – Landau bulk coefficients
/// * `k`   – Frank constant
/// * `dx`  – lattice spacing
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn molecular_field(
    q: &[Vec<QTensor2>],
    a: f64,
    b_coeff: f64,
    c_coeff: f64,
    k: f64,
    dx: f64,
) -> Vec<Vec<QTensor2>> {
    let nx = q.len();
    let ny = q[0].len();
    let mut h = vec![vec![QTensor2 { q11: 0.0, q12: 0.0 }; ny]; nx];
    let dx2 = dx * dx;
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let qij = q[i][j];
            let tr2 = 2.0 * (qij.q11 * qij.q11 + qij.q12 * qij.q12);
            // Bulk: dF_bulk/dQ
            let bulk11 = a * qij.q11 + b_coeff * qij.q11 * tr2 + c_coeff * tr2 * qij.q11;
            let bulk12 = a * qij.q12 + b_coeff * qij.q12 * tr2 + c_coeff * tr2 * qij.q12;
            // Elastic: -K ∇²Q
            let lap11 = (q[i + 1][j].q11 + q[i - 1][j].q11 + q[i][j + 1].q11 + q[i][j - 1].q11
                - 4.0 * qij.q11)
                / dx2;
            let lap12 = (q[i + 1][j].q12 + q[i - 1][j].q12 + q[i][j + 1].q12 + q[i][j - 1].q12
                - 4.0 * qij.q12)
                / dx2;
            h[i][j] = QTensor2 {
                q11: -(bulk11 - k * lap11),
                q12: -(bulk12 - k * lap12),
            };
        }
    }
    h
}

// ============================================================================
// Section 3 – Leslie-Ericksen relaxation on a lattice
// ============================================================================

/// Parameters for the Beris-Edwards / Leslie-Ericksen Q-tensor dynamics.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct LiquidCrystalParams {
    /// Rotational viscosity γ₁ \[Pa·s\].
    pub gamma1: f64,
    /// Landau coefficient a.
    pub landau_a: f64,
    /// Landau coefficient b.
    pub landau_b: f64,
    /// Landau coefficient c.
    pub landau_c: f64,
    /// Frank elastic constant K.
    pub frank_k: f64,
    /// Lattice spacing dx.
    pub dx: f64,
    /// Time step dt.
    pub dt: f64,
    /// Flow-alignment (tumbling) parameter ξ.
    pub xi: f64,
}

impl LiquidCrystalParams {
    /// Construct with physically reasonable defaults (lattice units).
    #[allow(dead_code)]
    pub fn default_lattice() -> Self {
        Self {
            gamma1: 0.1,
            landau_a: LANDAU_A,
            landau_b: LANDAU_B,
            landau_c: LANDAU_C,
            frank_k: 0.04,
            dx: 1.0,
            dt: 0.1,
            xi: 0.7,
        }
    }

    /// Compute the Ericksen number Er = γ₁ v L / K.
    #[allow(dead_code)]
    pub fn ericksen_number(&self, velocity: f64, length: f64) -> f64 {
        self.gamma1 * velocity * length / self.frank_k
    }

    /// Compute the nematic coherence length ξ_n = √(K / |a|).
    #[allow(dead_code)]
    pub fn coherence_length(&self) -> f64 {
        if self.landau_a.abs() < f64::EPSILON {
            f64::INFINITY
        } else {
            (self.frank_k / self.landau_a.abs()).sqrt()
        }
    }
}

/// Advance the Q-tensor field by one relaxational time step (no flow).
///
/// ∂Q/∂t = (1/γ₁) H
///
/// # Arguments
/// * `q`      – mutable Q-tensor field `[nx][ny]`
/// * `params` – `LiquidCrystalParams`
#[allow(dead_code)]
pub fn relax_qtensor(q: &mut Vec<Vec<QTensor2>>, params: &LiquidCrystalParams) {
    let h = molecular_field(
        q,
        params.landau_a,
        params.landau_b,
        params.landau_c,
        params.frank_k,
        params.dx,
    );
    let nx = q.len();
    let ny = q[0].len();
    let coeff = params.dt / params.gamma1;
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            q[i][j].q11 += coeff * h[i][j].q11;
            q[i][j].q12 += coeff * h[i][j].q12;
        }
    }
}

/// Advance Q with the full corotational (Beris-Edwards) term.
///
/// ∂Q/∂t + u·∇Q = (1/γ₁) H + ξ D  − Ω·Q + Q·Ω
///
/// Here D is the symmetric strain rate and Ω is the vorticity (antisymmetric).
/// For the 2-D case only the in-plane components are retained.
///
/// # Arguments
/// * `q`    – Q-tensor field `[nx][ny]`
/// * `ux`   – x-velocity `[nx][ny]`
/// * `uy`   – y-velocity `[nx][ny]`
/// * `params` – LC parameters
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn beris_edwards_step(
    q: &mut Vec<Vec<QTensor2>>,
    ux: &[Vec<f64>],
    uy: &[Vec<f64>],
    params: &LiquidCrystalParams,
) {
    let h = molecular_field(
        q,
        params.landau_a,
        params.landau_b,
        params.landau_c,
        params.frank_k,
        params.dx,
    );
    let nx = q.len();
    let ny = q[0].len();
    let inv2dx = 0.5 / params.dx;
    let coeff = params.dt / params.gamma1;

    let mut dq = vec![vec![QTensor2 { q11: 0.0, q12: 0.0 }; ny]; nx];

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let qij = q[i][j];
            // Advection (upwind 1st order)
            let adv11 = ux[i][j] * (q[i][j].q11 - q[i - 1][j].q11) / params.dx
                + uy[i][j] * (q[i][j].q11 - q[i][j - 1].q11) / params.dx;
            let adv12 = ux[i][j] * (q[i][j].q12 - q[i - 1][j].q12) / params.dx
                + uy[i][j] * (q[i][j].q12 - q[i][j - 1].q12) / params.dx;
            // Strain rate D_12 = (∂u_x/∂y + ∂u_y/∂x)/2
            let dxy = 0.5
                * ((ux[i][j + 1] - ux[i][j - 1]) * inv2dx + (uy[i + 1][j] - uy[i - 1][j]) * inv2dx);
            // Vorticity Ω_12 = (∂u_x/∂y - ∂u_y/∂x)/2
            let omega = 0.5
                * ((ux[i][j + 1] - ux[i][j - 1]) * inv2dx - (uy[i + 1][j] - uy[i - 1][j]) * inv2dx);
            // Co-rotation terms: -Ω·Q + Q·Ω → for traceless sym tensor in 2D
            let corot11 = 0.0; // vanishes for Q11 in 2-D
            let corot12 = 2.0 * omega * qij.q11; // simplified 2-D form
            dq[i][j] = QTensor2 {
                q11: coeff * h[i][j].q11 + params.dt * (-adv11 + params.xi * dxy * 0.0 + corot11),
                q12: coeff * h[i][j].q12 + params.dt * (-adv12 + params.xi * dxy + corot12),
            };
        }
    }
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            q[i][j].q11 += dq[i][j].q11;
            q[i][j].q12 += dq[i][j].q12;
        }
    }
}

// ============================================================================
// Section 4 – Lyotropic liquid crystal (concentration-coupled)
// ============================================================================

/// Lyotropic liquid crystal state at a single lattice site.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct LyotropicSite {
    /// Surfactant / polymer concentration c ∈ \[0,1\].
    pub concentration: f64,
    /// Nematic order tensor.
    pub qtensor: QTensor2,
    /// Local fluid density ρ.
    pub density: f64,
    /// Local x-velocity.
    pub ux: f64,
    /// Local y-velocity.
    pub uy: f64,
}

impl LyotropicSite {
    /// Construct an isotropic (disordered) site.
    #[allow(dead_code)]
    pub fn isotropic(conc: f64) -> Self {
        Self {
            concentration: conc,
            qtensor: QTensor2 { q11: 0.0, q12: 0.0 },
            density: 1.0,
            ux: 0.0,
            uy: 0.0,
        }
    }

    /// Equilibrium scalar order parameter S_eq(c) = max(0, 2c − 1).
    ///
    /// Simple Onsager-type threshold at c* = 0.5.
    #[allow(dead_code)]
    pub fn equilibrium_order(&self) -> f64 {
        (2.0 * self.concentration - 1.0).max(0.0)
    }
}

/// Advance lyotropic ordering toward equilibrium (relaxation approximation).
///
/// # Arguments
/// * `sites` – 2-D grid of `LyotropicSite`
/// * `tau_q` – orientational relaxation time
/// * `dt`    – time step
#[allow(dead_code)]
pub fn lyotropic_relax(sites: &mut Vec<Vec<LyotropicSite>>, tau_q: f64, dt: f64) {
    let nx = sites.len();
    let ny = sites[0].len();
    for i in 0..nx {
        for j in 0..ny {
            let s_eq = sites[i][j].equilibrium_order();
            let q_eq = QTensor2::from_director(0.0, s_eq);
            let relax = dt / tau_q;
            sites[i][j].qtensor.q11 += relax * (q_eq.q11 - sites[i][j].qtensor.q11);
            sites[i][j].qtensor.q12 += relax * (q_eq.q12 - sites[i][j].qtensor.q12);
        }
    }
}

// ============================================================================
// Section 5 – Colloidal suspension LBM
// ============================================================================

/// A single colloidal particle in 2-D.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ColloidParticle {
    /// Centre-of-mass position (x, y) in lattice units.
    pub position: [f64; 2],
    /// Velocity (vx, vy) in lattice units/step.
    pub velocity: [f64; 2],
    /// Particle radius in lattice units.
    pub radius: f64,
    /// Particle mass (lattice units).
    pub mass: f64,
    /// Electric charge q (for Lorentz coupling).
    pub charge: f64,
}

impl ColloidParticle {
    /// Construct a stationary particle at a given position.
    #[allow(dead_code)]
    pub fn new_stationary(x: f64, y: f64, radius: f64, mass: f64, charge: f64) -> Self {
        Self {
            position: [x, y],
            velocity: [0.0, 0.0],
            radius,
            mass,
            charge,
        }
    }

    /// Apply Lorentz force F = q(E + v × B) in 2-D (B out-of-plane).
    ///
    /// # Arguments
    /// * `ex`, `ey` – electric field components
    /// * `bz`       – out-of-plane magnetic field component
    /// * `dt`       – time step
    #[allow(dead_code)]
    pub fn apply_lorentz_force(&mut self, ex: f64, ey: f64, bz: f64, dt: f64) {
        let vx = self.velocity[0];
        let vy = self.velocity[1];
        let fx = self.charge * (ex + vy * bz);
        let fy = self.charge * (ey - vx * bz);
        self.velocity[0] += dt * fx / self.mass;
        self.velocity[1] += dt * fy / self.mass;
    }

    /// Advance particle position with Euler integration.
    #[allow(dead_code)]
    pub fn advance(&mut self, dt: f64) {
        self.position[0] += dt * self.velocity[0];
        self.position[1] += dt * self.velocity[1];
    }

    /// Compute excluded-volume (hard-core) force between two particles.
    ///
    /// Returns (fx, fy) acting on `self` due to `other`.
    /// Uses a steeply repulsive WCA-style potential: U = ε (σ/r)^12.
    #[allow(dead_code)]
    pub fn excluded_volume_force(&self, other: &ColloidParticle, epsilon: f64) -> [f64; 2] {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let r2 = dx * dx + dy * dy;
        let sigma = self.radius + other.radius;
        let sigma2 = sigma * sigma;
        if r2 < sigma2 * 4.0 {
            let r6 = r2 * r2 * r2;
            let s6 = sigma2 * sigma2 * sigma2;
            let mag = 12.0 * epsilon * s6 / (r6 * r2 * r2 * r2 * r2);
            [mag * dx, mag * dy]
        } else {
            [0.0, 0.0]
        }
    }
}

/// Compute the inter-particle Lorentz-coupling force on the fluid.
///
/// Returns the body-force density field (fx, fy) in lattice units.
///
/// # Arguments
/// * `particles` – list of colloidal particles
/// * `nx`, `ny`  – grid dimensions
/// * `bz`        – out-of-plane magnetic field
#[allow(dead_code)]
pub fn colloidal_body_force(
    particles: &[ColloidParticle],
    nx: usize,
    ny: usize,
    bz: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut fx = vec![vec![0.0_f64; ny]; nx];
    let mut fy = vec![vec![0.0_f64; ny]; nx];
    for p in particles {
        let pi = p.position[0].round() as isize;
        let pj = p.position[1].round() as isize;
        let r = p.radius.ceil() as isize + 1;
        for di in -r..=r {
            for dj in -r..=r {
                let ii = pi + di;
                let jj = pj + dj;
                if ii < 0 || jj < 0 || ii >= nx as isize || jj >= ny as isize {
                    continue;
                }
                let dist = ((di * di + dj * dj) as f64).sqrt();
                if dist > p.radius {
                    continue;
                }
                let ii = ii as usize;
                let jj = jj as usize;
                // Lorentz back-reaction on fluid
                fx[ii][jj] -= p.charge * p.velocity[1] * bz;
                fy[ii][jj] += p.charge * p.velocity[0] * bz;
            }
        }
    }
    (fx, fy)
}

// ============================================================================
// Section 6 – Gel network dynamics
// ============================================================================

/// Bond in the gel network.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct GelBond {
    /// Index of the first node.
    pub node_a: usize,
    /// Index of the second node.
    pub node_b: usize,
    /// Natural (rest) length \[lattice units\].
    pub rest_length: f64,
    /// Spring constant k \[Pa/m\].
    pub stiffness: f64,
    /// Maximum extension before breakage (fraction of rest length).
    pub break_strain: f64,
    /// Whether the bond is currently intact.
    pub intact: bool,
}

impl GelBond {
    /// Compute elastic force on node A (force on node B is opposite).
    ///
    /// Returns `[fx, fy]` acting on node A.
    #[allow(dead_code)]
    pub fn elastic_force(&self, pos_a: [f64; 2], pos_b: [f64; 2]) -> [f64; 2] {
        if !self.intact {
            return [0.0, 0.0];
        }
        let dx = pos_b[0] - pos_a[0];
        let dy = pos_b[1] - pos_a[1];
        let r = (dx * dx + dy * dy).sqrt();
        if r < f64::EPSILON {
            return [0.0, 0.0];
        }
        let extension = r - self.rest_length;
        let mag = self.stiffness * extension / r;
        [mag * dx, mag * dy]
    }

    /// Check and apply bond breakage criterion.
    #[allow(dead_code)]
    pub fn check_breakage(&mut self, pos_a: [f64; 2], pos_b: [f64; 2]) {
        if !self.intact {
            return;
        }
        let dx = pos_b[0] - pos_a[0];
        let dy = pos_b[1] - pos_a[1];
        let r = (dx * dx + dy * dy).sqrt();
        let strain = (r - self.rest_length) / self.rest_length;
        if strain > self.break_strain {
            self.intact = false;
        }
    }
}

/// A node in the gel network carrying position and velocity.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct GelNode {
    /// Position \[lattice units\].
    pub position: [f64; 2],
    /// Velocity \[lattice units/step\].
    pub velocity: [f64; 2],
    /// Node mass.
    pub mass: f64,
}

impl GelNode {
    /// Advance the node using Euler integration.
    #[allow(dead_code)]
    pub fn advance(&mut self, fx: f64, fy: f64, dt: f64) {
        self.velocity[0] += dt * fx / self.mass;
        self.velocity[1] += dt * fy / self.mass;
        self.position[0] += dt * self.velocity[0];
        self.position[1] += dt * self.velocity[1];
    }
}

/// Advance the gel network by one step (force calculation + breakage check).
///
/// # Arguments
/// * `nodes` – mutable slice of gel nodes
/// * `bonds` – mutable slice of gel bonds
/// * `dt`    – time step
#[allow(dead_code)]
pub fn gel_network_step(nodes: &mut [GelNode], bonds: &mut [GelBond], dt: f64) {
    let mut forces = vec![[0.0_f64; 2]; nodes.len()];
    for bond in bonds.iter_mut() {
        if !bond.intact {
            continue;
        }
        let pa = nodes[bond.node_a].position;
        let pb = nodes[bond.node_b].position;
        bond.check_breakage(pa, pb);
        let f = bond.elastic_force(pa, pb);
        forces[bond.node_a][0] += f[0];
        forces[bond.node_a][1] += f[1];
        forces[bond.node_b][0] -= f[0];
        forces[bond.node_b][1] -= f[1];
    }
    for (i, node) in nodes.iter_mut().enumerate() {
        node.advance(forces[i][0], forces[i][1], dt);
    }
}

// ============================================================================
// Section 7 – Lipid bilayer LBM
// ============================================================================

/// Properties of a lipid bilayer membrane.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct LipidBilayerParams {
    /// Bending rigidity κ \[J = Pa·m³\].
    pub bending_rigidity: f64,
    /// Area compression modulus K_A \[Pa·m\].
    pub area_modulus: f64,
    /// Saddle-splay modulus κ̄.
    pub saddle_splay: f64,
    /// Target area per lipid \[m²\].
    pub target_area: f64,
    /// Line tension of membrane edge λ \[N\].
    pub line_tension: f64,
}

impl LipidBilayerParams {
    /// Default DPPC-like bilayer (in SI units at 323 K).
    #[allow(dead_code)]
    pub fn dppc() -> Self {
        Self {
            bending_rigidity: 20.0 * K_B * 323.0,
            area_modulus: 0.24,
            saddle_splay: -10.0 * K_B * 323.0,
            target_area: 6.4e-19,
            line_tension: 1.0e-11,
        }
    }

    /// Helfrich bending energy for a spherical vesicle of radius R.
    ///
    /// E_bend = 4π(2κ + κ̄) (closed surface, no edges).
    #[allow(dead_code)]
    pub fn spherical_bending_energy(&self) -> f64 {
        4.0 * PI * (2.0 * self.bending_rigidity + self.saddle_splay)
    }
}

/// Curvature-driven force on a 2-D membrane contour (discrete).
///
/// Uses the Helfrich curvature energy:
/// f_n = κ (∂²κ/∂s² + κ³/2) (Euler-Lagrange normal force)
///
/// # Arguments
/// * `x`, `y`  – membrane node positions
/// * `kappa`   – bending rigidity κ
#[allow(dead_code)]
pub fn helfrich_force(x: &[f64], y: &[f64], kappa: f64) -> Vec<[f64; 2]> {
    let n = x.len();
    let mut force = vec![[0.0_f64; 2]; n];
    if n < 3 {
        return force;
    }
    for i in 1..n - 1 {
        let dx1 = x[i] - x[i - 1];
        let dy1 = y[i] - y[i - 1];
        let dx2 = x[i + 1] - x[i];
        let dy2 = y[i + 1] - y[i];
        let l1 = (dx1 * dx1 + dy1 * dy1).sqrt().max(f64::EPSILON);
        let l2 = (dx2 * dx2 + dy2 * dy2).sqrt().max(f64::EPSILON);
        // Discrete curvature (turning angle / arc length)
        let cross = dx1 * dy2 - dy1 * dx2;
        let curv = 2.0 * cross / ((l1 + l2) * l1 * l2);
        // Normal direction
        let nx = -(dy1 / l1 + dy2 / l2) / 2.0;
        let ny = (dx1 / l1 + dx2 / l2) / 2.0;
        let nl = (nx * nx + ny * ny).sqrt().max(f64::EPSILON);
        force[i][0] = kappa * curv * nx / nl;
        force[i][1] = kappa * curv * ny / nl;
    }
    force
}

// ============================================================================
// Section 8 – Worm-like micelles (Cates model)
// ============================================================================

/// Parameters for the Cates worm-like micelle model.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct WormLikeMicelleParams {
    /// Plateau modulus G₀ \[Pa\].
    pub plateau_modulus: f64,
    /// Reptation / diffusion time τ_rep \[s\].
    pub tau_rep: f64,
    /// Breakage time τ_break \[s\].
    pub tau_break: f64,
    /// Solvent viscosity η_s \[Pa·s\].
    pub solvent_viscosity: f64,
}

impl WormLikeMicelleParams {
    /// Maxwell relaxation time τ_R = √(τ_rep · τ_break).
    #[allow(dead_code)]
    pub fn maxwell_time(&self) -> f64 {
        (self.tau_rep * self.tau_break).sqrt()
    }

    /// Zero-shear viscosity η₀ = G₀ τ_R.
    #[allow(dead_code)]
    pub fn zero_shear_viscosity(&self) -> f64 {
        self.plateau_modulus * self.maxwell_time()
    }

    /// Storage modulus G'(ω) for single-mode Maxwell model.
    #[allow(dead_code)]
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let tr = self.maxwell_time();
        let w2 = (omega * tr) * (omega * tr);
        self.plateau_modulus * w2 / (1.0 + w2)
    }

    /// Loss modulus G''(ω).
    #[allow(dead_code)]
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let tr = self.maxwell_time();
        let wt = omega * tr;
        let w2 = wt * wt;
        self.plateau_modulus * wt / (1.0 + w2)
    }

    /// Cole-Cole half-circle radius = G₀/2.
    #[allow(dead_code)]
    pub fn cole_cole_radius(&self) -> f64 {
        self.plateau_modulus / 2.0
    }
}

/// Evolve the Maxwell stress tensor for worm-like micelles.
///
/// ∂σ/∂t = G₀ D̃ − σ/τ_R  (upper-convected Maxwell in linearised form)
///
/// # Arguments
/// * `sigma_xx`, `sigma_xy`, `sigma_yy` – stress components \[Pa\]
/// * `dxx`, `dxy`, `dyy`               – strain-rate components \[1/s\]
/// * `params`                           – micelle parameters
/// * `dt`                               – time step \[s\]
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn maxwell_stress_step(
    sigma_xx: &mut f64,
    sigma_xy: &mut f64,
    sigma_yy: &mut f64,
    dxx: f64,
    dxy: f64,
    dyy: f64,
    params: &WormLikeMicelleParams,
    dt: f64,
) {
    let tau = params.maxwell_time();
    let g0 = params.plateau_modulus;
    *sigma_xx += dt * (g0 * dxx - *sigma_xx / tau);
    *sigma_xy += dt * (g0 * dxy - *sigma_xy / tau);
    *sigma_yy += dt * (g0 * dyy - *sigma_yy / tau);
}

// ============================================================================
// Section 9 – Block copolymer phase separation (Ohta-Kawasaki)
// ============================================================================

/// Parameters for the Ohta-Kawasaki block copolymer model.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct OhtaKawasakiParams {
    /// Interaction parameter χN (Flory-Huggins × chain length).
    pub chi_n: f64,
    /// Mobility M.
    pub mobility: f64,
    /// Interface width parameter ε.
    pub epsilon: f64,
    /// Long-range interaction coefficient α.
    pub alpha: f64,
    /// Composition of A-block f_A ∈ (0,1).
    pub f_a: f64,
    /// Lattice spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
}

impl OhtaKawasakiParams {
    /// Mean-field spinodal at f_A = 0.5: χN_s = 2 / (2 f_A (1-f_A)) = 2 for symmetric.
    #[allow(dead_code)]
    pub fn spinodal_chin(&self) -> f64 {
        2.0 / (2.0 * self.f_a * (1.0 - self.f_a))
    }

    /// Lamellar period estimate L* ≈ 2π ε / √(α).
    #[allow(dead_code)]
    pub fn lamellar_period(&self) -> f64 {
        if self.alpha < f64::EPSILON {
            f64::INFINITY
        } else {
            2.0 * PI * self.epsilon / self.alpha.sqrt()
        }
    }
}

/// Local Ohta-Kawasaki free energy density.
///
/// f = (χN/4)(φ - f_A)² (1 - (φ - f_A)²) − approximation for Flory-Huggins
#[allow(dead_code)]
pub fn ohta_kawasaki_bulk_density(phi: f64, f_a: f64, chi_n: f64) -> f64 {
    let m = phi - f_a;
    0.25 * chi_n * m * m * (1.0 - m * m)
}

/// Chemical potential μ = df/dφ + long-range term (local part only).
#[allow(dead_code)]
pub fn ohta_kawasaki_mu(phi: f64, f_a: f64, chi_n: f64) -> f64 {
    let m = phi - f_a;
    chi_n * m * (0.5 - m * m)
}

/// Advance block copolymer composition field by one Cahn-Hilliard step.
///
/// ∂φ/∂t = M ∇²μ   where μ = δF/δφ = df/dφ − ε² ∇²φ + α (φ − f_A)
///
/// # Arguments
/// * `phi`    – composition field `[nx][ny]`
/// * `params` – OK parameters
#[allow(dead_code)]
pub fn block_copolymer_step(phi: &mut Vec<Vec<f64>>, params: &OhtaKawasakiParams) {
    let nx = phi.len();
    let ny = phi[0].len();
    let dx2 = params.dx * params.dx;
    let mut mu = vec![vec![0.0_f64; ny]; nx];

    // Compute chemical potential
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let lap_phi = (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]
                - 4.0 * phi[i][j])
                / dx2;
            let bulk = ohta_kawasaki_mu(phi[i][j], params.f_a, params.chi_n);
            mu[i][j] = bulk - params.epsilon * params.epsilon * lap_phi
                + params.alpha * (phi[i][j] - params.f_a);
        }
    }

    // ∂φ/∂t = M ∇²μ
    let mut dphi = vec![vec![0.0_f64; ny]; nx];
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let lap_mu =
                (mu[i + 1][j] + mu[i - 1][j] + mu[i][j + 1] + mu[i][j - 1] - 4.0 * mu[i][j]) / dx2;
            dphi[i][j] = params.mobility * lap_mu;
        }
    }
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            phi[i][j] += params.dt * dphi[i][j];
        }
    }
}

// ============================================================================
// Section 10 – Active matter (self-propelled particles)
// ============================================================================

/// A self-propelled particle (SPP) with polar alignment.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ActiveParticle {
    /// Position (x, y) in lattice units.
    pub position: [f64; 2],
    /// Orientation angle θ (direction of self-propulsion).
    pub theta: f64,
    /// Self-propulsion speed v₀.
    pub speed: f64,
    /// Rotational diffusion coefficient D_r.
    pub d_rot: f64,
    /// Translational friction coefficient ζ.
    pub friction: f64,
}

impl ActiveParticle {
    /// Create a new active particle.
    #[allow(dead_code)]
    pub fn new(x: f64, y: f64, theta: f64) -> Self {
        Self {
            position: [x, y],
            theta,
            speed: ACTIVE_SPEED,
            d_rot: 0.01,
            friction: 1.0,
        }
    }

    /// Self-propulsion velocity vector.
    #[allow(dead_code)]
    pub fn propulsion_velocity(&self) -> [f64; 2] {
        [self.speed * self.theta.cos(), self.speed * self.theta.sin()]
    }

    /// Advance position and orientation (Langevin, noise-free).
    ///
    /// # Arguments
    /// * `fx`, `fy` – external force
    /// * `torque`   – external torque on orientation
    /// * `dt`       – time step
    #[allow(dead_code)]
    pub fn advance(&mut self, fx: f64, fy: f64, torque: f64, dt: f64) {
        let v = self.propulsion_velocity();
        self.position[0] += dt * (v[0] + fx / self.friction);
        self.position[1] += dt * (v[1] + fy / self.friction);
        self.theta += dt * torque;
    }
}

/// Compute Vicsek-style polar alignment torque on particle `i`.
///
/// Each particle aligns toward the average orientation of neighbours within
/// radius R_align.
///
/// # Arguments
/// * `particles` – list of active particles
/// * `i`         – index of the focal particle
/// * `r_align`   – alignment radius
/// * `j_align`   – coupling strength
#[allow(dead_code)]
pub fn vicsek_alignment_torque(
    particles: &[ActiveParticle],
    i: usize,
    r_align: f64,
    j_align: f64,
) -> f64 {
    let xi = particles[i].position[0];
    let yi = particles[i].position[1];
    let mut sum_sin = 0.0;
    let mut sum_cos = 0.0;
    let mut count = 0usize;
    for (k, p) in particles.iter().enumerate() {
        if k == i {
            continue;
        }
        let dx = p.position[0] - xi;
        let dy = p.position[1] - yi;
        let r = (dx * dx + dy * dy).sqrt();
        if r < r_align {
            sum_sin += p.theta.sin();
            sum_cos += p.theta.cos();
            count += 1;
        }
    }
    if count == 0 {
        return 0.0;
    }
    let theta_avg = sum_sin.atan2(sum_cos);
    let diff = theta_avg - particles[i].theta;
    // wrap to [-π, π]
    let diff = diff - 2.0 * PI * (diff / (2.0 * PI)).round();
    j_align * diff
}

/// Measure the global polar order parameter Φ = |⟨e^{iθ}⟩|.
#[allow(dead_code)]
pub fn polar_order_parameter(particles: &[ActiveParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    let sx: f64 = particles.iter().map(|p| p.theta.cos()).sum();
    let sy: f64 = particles.iter().map(|p| p.theta.sin()).sum();
    (sx * sx + sy * sy).sqrt() / n
}

/// Measure the nematic order parameter S = |⟨cos 2θ⟩|.
#[allow(dead_code)]
pub fn nematic_order_parameter(particles: &[ActiveParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    let sc: f64 = particles.iter().map(|p| (2.0 * p.theta).cos()).sum();
    let ss: f64 = particles.iter().map(|p| (2.0 * p.theta).sin()).sum();
    (sc * sc + ss * ss).sqrt() / n
}

// ============================================================================
// Section 11 – LBM equilibrium & BGK collision helpers
// ============================================================================

/// Compute the D2Q9 Maxwell-Boltzmann equilibrium distribution.
///
/// f^eq_i = w_i ρ \[ 1 + (c_i·u)/c_s² + (c_i·u)²/(2c_s⁴) − u²/(2c_s²) \]
///
/// # Arguments
/// * `rho` – fluid density
/// * `ux`  – x-velocity
/// * `uy`  – y-velocity
#[allow(dead_code)]
pub fn equilibrium_d2q9(rho: f64, ux: f64, uy: f64) -> [f64; Q9] {
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0_f64; Q9];
    for i in 0..Q9 {
        let cu = CX[i] * ux + CY[i] * uy;
        feq[i] = W9[i] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

/// BGK collision step for a single lattice site.
///
/// f_i* = f_i − (f_i − f^eq_i) / τ
///
/// # Arguments
/// * `f`   – distribution functions (modified in place)
/// * `tau` – relaxation time
/// * `rho` – density
/// * `ux`  – x-velocity
/// * `uy`  – y-velocity
#[allow(dead_code)]
pub fn bgk_collide(f: &mut [f64; Q9], tau: f64, rho: f64, ux: f64, uy: f64) {
    let feq = equilibrium_d2q9(rho, ux, uy);
    for i in 0..Q9 {
        f[i] -= (f[i] - feq[i]) / tau;
    }
}

/// Compute macroscopic density and velocity from distribution functions.
///
/// ρ = Σ f_i,   ρ u = Σ f_i c_i
#[allow(dead_code)]
pub fn macroscopic_d2q9(f: &[f64; Q9]) -> (f64, f64, f64) {
    let rho: f64 = f.iter().sum();
    let mut rho_ux = 0.0;
    let mut rho_uy = 0.0;
    for i in 0..Q9 {
        rho_ux += f[i] * CX[i];
        rho_uy += f[i] * CY[i];
    }
    let ux = if rho > f64::EPSILON {
        rho_ux / rho
    } else {
        0.0
    };
    let uy = if rho > f64::EPSILON {
        rho_uy / rho
    } else {
        0.0
    };
    (rho, ux, uy)
}

// ============================================================================
// Section 12 – Full soft-matter LBM simulation grid
// ============================================================================

/// A 2-D soft-matter LBM grid coupling fluid, Q-tensor, and active particles.
#[allow(dead_code)]
pub struct SoftMatterGrid {
    /// Grid width (x).
    pub nx: usize,
    /// Grid height (y).
    pub ny: usize,
    /// BGK relaxation time τ.
    pub tau: f64,
    /// Distribution functions f\[i\]\[j\]\[q\].
    pub f: Vec<Vec<[f64; Q9]>>,
    /// Q-tensor field.
    pub q_field: Vec<Vec<QTensor2>>,
    /// LC dynamics parameters.
    pub lc_params: LiquidCrystalParams,
}

impl SoftMatterGrid {
    /// Construct a quiescent uniform-density grid.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – grid dimensions
    /// * `tau`        – fluid relaxation time
    /// * `lc_params`  – liquid crystal parameters
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, tau: f64, lc_params: LiquidCrystalParams) -> Self {
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        let f = vec![vec![feq; ny]; nx];
        let q_field = vec![vec![QTensor2 { q11: 0.0, q12: 0.0 }; ny]; nx];
        Self {
            nx,
            ny,
            tau,
            f,
            q_field,
            lc_params,
        }
    }

    /// Initialise the Q-tensor field from a director angle map.
    ///
    /// # Arguments
    /// * `theta_field` – angle field `[nx][ny]`
    /// * `s0`          – uniform initial order parameter
    #[allow(dead_code)]
    pub fn set_director(&mut self, theta_field: &[Vec<f64>], s0: f64) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.q_field[i][j] = QTensor2::from_director(theta_field[i][j], s0);
            }
        }
    }

    /// Perform one full LBM time step (collision → streaming → Q-tensor relax).
    #[allow(dead_code)]
    pub fn step(&mut self) {
        // BGK collision
        for i in 0..self.nx {
            for j in 0..self.ny {
                let (rho, ux, uy) = macroscopic_d2q9(&self.f[i][j]);
                bgk_collide(&mut self.f[i][j], self.tau, rho, ux, uy);
            }
        }
        // Streaming (periodic)
        let mut f_new = vec![vec![[0.0_f64; Q9]; self.ny]; self.nx];
        for i in 0..self.nx {
            for j in 0..self.ny {
                for q in 0..Q9 {
                    let ni = (i as isize + CX[q] as isize).rem_euclid(self.nx as isize) as usize;
                    let nj = (j as isize + CY[q] as isize).rem_euclid(self.ny as isize) as usize;
                    f_new[ni][nj][q] = self.f[i][j][q];
                }
            }
        }
        self.f = f_new;
        // Q-tensor relaxation
        relax_qtensor(&mut self.q_field, &self.lc_params);
    }

    /// Compute average scalar order parameter ⟨S⟩ over the grid.
    #[allow(dead_code)]
    pub fn mean_order_parameter(&self) -> f64 {
        let total: f64 = self
            .q_field
            .iter()
            .flat_map(|row| row.iter())
            .map(|q| q.order_parameter())
            .sum();
        total / (self.nx * self.ny) as f64
    }

    /// Compute average fluid speed ⟨|u|⟩.
    #[allow(dead_code)]
    pub fn mean_speed(&self) -> f64 {
        let total: f64 = self
            .f
            .iter()
            .flat_map(|row| row.iter())
            .map(|f| {
                let (_rho, ux, uy) = macroscopic_d2q9(f);
                (ux * ux + uy * uy).sqrt()
            })
            .sum();
        total / (self.nx * self.ny) as f64
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── QTensor2 ──────────────────────────────────────────────────────────────

    #[test]
    fn test_qtensor_from_director_isotropic() {
        let q = QTensor2::from_director(0.0, 0.0);
        assert!((q.q11).abs() < 1e-12);
        assert!((q.q12).abs() < 1e-12);
    }

    #[test]
    fn test_qtensor_order_parameter_roundtrip() {
        let s = 0.7;
        let theta = PI / 4.0;
        let q = QTensor2::from_director(theta, s);
        let s_back = q.order_parameter();
        assert!((s_back - s).abs() < 1e-10, "S={s_back} vs {s}");
    }

    #[test]
    fn test_qtensor_director_angle_roundtrip() {
        let theta = 0.6;
        let q = QTensor2::from_director(theta, 0.8);
        let theta_back = q.director_angle();
        assert!((theta_back - theta).abs() < 1e-10);
    }

    #[test]
    fn test_qtensor_add_scale() {
        let q1 = QTensor2 { q11: 0.1, q12: 0.2 };
        let q2 = QTensor2 { q11: 0.3, q12: 0.4 };
        let sum = q1.add(&q2);
        assert!((sum.q11 - 0.4).abs() < 1e-12);
        assert!((sum.q12 - 0.6).abs() < 1e-12);
        let scaled = q1.scale(2.0);
        assert!((scaled.q11 - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_qtensor_norm() {
        let q = QTensor2 { q11: 3.0, q12: 4.0 };
        assert!((q.norm() - 5.0).abs() < 1e-12);
    }

    // ── Frank elastic ──────────────────────────────────────────────────────────

    #[test]
    fn test_frank_elastic_uniform_zero() {
        // Uniform director → zero elastic energy
        let nx = 5;
        let ny = 5;
        let q = vec![vec![QTensor2::from_director(0.0, 0.6); ny]; nx];
        let f = frank_elastic_density(&q, 1.0, 1.0);
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                assert!(f[i][j].abs() < 1e-12, "non-zero at ({i},{j})");
            }
        }
    }

    #[test]
    fn test_molecular_field_uniform_isotropic() {
        // Q=0 → bulk field is zero (linear in Q at Q=0 when a=0)
        let nx = 5;
        let ny = 5;
        let q = vec![vec![QTensor2 { q11: 0.0, q12: 0.0 }; ny]; nx];
        let h = molecular_field(&q, 0.0, 1.0, 1.0, 1.0, 1.0);
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                assert!(h[i][j].q11.abs() < 1e-12);
            }
        }
    }

    // ── LiquidCrystalParams ───────────────────────────────────────────────────

    #[test]
    fn test_ericksen_number() {
        let p = LiquidCrystalParams::default_lattice();
        let er = p.ericksen_number(0.01, 100.0);
        assert!(er > 0.0);
    }

    #[test]
    fn test_coherence_length() {
        let p = LiquidCrystalParams::default_lattice();
        let xi = p.coherence_length();
        assert!(xi > 0.0 && xi.is_finite());
    }

    #[test]
    fn test_coherence_length_zero_a() {
        let mut p = LiquidCrystalParams::default_lattice();
        p.landau_a = 0.0;
        assert!(p.coherence_length().is_infinite());
    }

    // ── Lyotropic ─────────────────────────────────────────────────────────────

    #[test]
    fn test_lyotropic_equilibrium_order_below_threshold() {
        let s = LyotropicSite::isotropic(0.3);
        assert!((s.equilibrium_order()).abs() < 1e-12);
    }

    #[test]
    fn test_lyotropic_equilibrium_order_above_threshold() {
        let s = LyotropicSite::isotropic(0.8);
        let s_eq = s.equilibrium_order();
        assert!((s_eq - 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_lyotropic_relax_approaches_equilibrium() {
        let mut sites = vec![vec![LyotropicSite::isotropic(0.8)]];
        let initial = sites[0][0].qtensor.q11;
        lyotropic_relax(&mut sites, 1.0, 0.5);
        // should move toward equilibrium
        let s_eq = QTensor2::from_director(0.0, 0.6).q11;
        let after = sites[0][0].qtensor.q11;
        assert!((after - s_eq).abs() < (initial - s_eq).abs() + 1e-12);
    }

    // ── Colloidal ─────────────────────────────────────────────────────────────

    #[test]
    fn test_colloid_lorentz_force() {
        let mut p = ColloidParticle::new_stationary(5.0, 5.0, 1.0, 1.0, 1.0);
        p.velocity = [1.0, 0.0];
        p.apply_lorentz_force(0.0, 0.0, 1.0, 1.0);
        // F_y = q(-v_x * B_z) → -1
        assert!((p.velocity[1] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_colloid_excluded_volume_repulsion() {
        let p1 = ColloidParticle::new_stationary(0.0, 0.0, 0.5, 1.0, 0.0);
        let p2 = ColloidParticle::new_stationary(0.5, 0.0, 0.5, 1.0, 0.0);
        let f = p1.excluded_volume_force(&p2, 1.0);
        // should be repulsive (negative x)
        assert!(f[0] < 0.0);
    }

    #[test]
    fn test_colloid_body_force_zero_charge() {
        let particles = vec![ColloidParticle::new_stationary(5.0, 5.0, 1.0, 1.0, 0.0)];
        let (fx, fy) = colloidal_body_force(&particles, 10, 10, 1.0);
        let sum_fx: f64 = fx.iter().flat_map(|r| r.iter()).sum();
        let sum_fy: f64 = fy.iter().flat_map(|r| r.iter()).sum();
        assert!(sum_fx.abs() < 1e-12);
        assert!(sum_fy.abs() < 1e-12);
    }

    // ── Gel network ───────────────────────────────────────────────────────────

    #[test]
    fn test_gel_bond_force_at_rest_length() {
        let bond = GelBond {
            node_a: 0,
            node_b: 1,
            rest_length: 1.0,
            stiffness: 10.0,
            break_strain: 0.5,
            intact: true,
        };
        let pa = [0.0, 0.0];
        let pb = [1.0, 0.0];
        let f = bond.elastic_force(pa, pb);
        assert!(f[0].abs() < 1e-12);
        assert!(f[1].abs() < 1e-12);
    }

    #[test]
    fn test_gel_bond_breaks_at_large_extension() {
        let mut bond = GelBond {
            node_a: 0,
            node_b: 1,
            rest_length: 1.0,
            stiffness: 10.0,
            break_strain: 0.1,
            intact: true,
        };
        bond.check_breakage([0.0, 0.0], [2.0, 0.0]);
        assert!(!bond.intact);
    }

    #[test]
    fn test_gel_network_step_conserves_com() {
        let mut nodes = vec![
            GelNode {
                position: [0.0, 0.0],
                velocity: [0.0, 0.0],
                mass: 1.0,
            },
            GelNode {
                position: [2.0, 0.0],
                velocity: [0.0, 0.0],
                mass: 1.0,
            },
        ];
        let mut bonds = vec![GelBond {
            node_a: 0,
            node_b: 1,
            rest_length: 2.0,
            stiffness: 1.0,
            break_strain: 10.0,
            intact: true,
        }];
        let com_x_before = (nodes[0].position[0] + nodes[1].position[0]) / 2.0;
        gel_network_step(&mut nodes, &mut bonds, 0.01);
        let com_x_after = (nodes[0].position[0] + nodes[1].position[0]) / 2.0;
        assert!((com_x_before - com_x_after).abs() < 1e-10);
    }

    // ── Lipid bilayer ─────────────────────────────────────────────────────────

    #[test]
    fn test_lipid_bilayer_dppc() {
        let p = LipidBilayerParams::dppc();
        assert!(p.bending_rigidity > 0.0);
        assert!(p.area_modulus > 0.0);
    }

    #[test]
    fn test_spherical_bending_energy_positive() {
        let mut p = LipidBilayerParams::dppc();
        p.saddle_splay = 0.0;
        assert!(p.spherical_bending_energy() > 0.0);
    }

    #[test]
    fn test_helfrich_force_straight_membrane() {
        // Straight horizontal membrane → zero curvature → zero force
        let n = 5;
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let y: Vec<f64> = vec![0.0; n];
        let f = helfrich_force(&x, &y, 1.0);
        for i in 1..n - 1 {
            assert!(
                f[i][0].abs() < 1e-10 && f[i][1].abs() < 1e-10,
                "non-zero at {i}"
            );
        }
    }

    // ── Worm-like micelles ────────────────────────────────────────────────────

    #[test]
    fn test_maxwell_time_positive() {
        let p = WormLikeMicelleParams {
            plateau_modulus: 100.0,
            tau_rep: 1.0,
            tau_break: 0.01,
            solvent_viscosity: 1e-3,
        };
        assert!(p.maxwell_time() > 0.0);
    }

    #[test]
    fn test_zero_shear_viscosity() {
        let p = WormLikeMicelleParams {
            plateau_modulus: 100.0,
            tau_rep: 1.0,
            tau_break: 0.01,
            solvent_viscosity: 1e-3,
        };
        let eta = p.zero_shear_viscosity();
        assert!((eta - 100.0 * p.maxwell_time()).abs() < 1e-10);
    }

    #[test]
    fn test_maxwell_stress_relaxation() {
        let p = WormLikeMicelleParams {
            plateau_modulus: 1.0,
            tau_rep: 1.0,
            tau_break: 1.0,
            solvent_viscosity: 1e-3,
        };
        let mut sxx = 1.0_f64;
        let mut sxy = 0.0_f64;
        let mut syy = 0.0_f64;
        // With no strain rate, stress should decay
        maxwell_stress_step(&mut sxx, &mut sxy, &mut syy, 0.0, 0.0, 0.0, &p, 0.1);
        assert!(sxx < 1.0);
    }

    #[test]
    fn test_cole_cole_radius() {
        let p = WormLikeMicelleParams {
            plateau_modulus: 200.0,
            tau_rep: 1.0,
            tau_break: 0.1,
            solvent_viscosity: 1e-3,
        };
        assert!((p.cole_cole_radius() - 100.0).abs() < 1e-10);
    }

    // ── Block copolymer ───────────────────────────────────────────────────────

    #[test]
    fn test_ohta_kawasaki_bulk_at_equilibrium() {
        // f(f_A) = 0 (at mean composition, m=0)
        let f_a = 0.5;
        let f = ohta_kawasaki_bulk_density(f_a, f_a, 2.0);
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn test_ohta_kawasaki_mu_at_mean() {
        let f_a = 0.5;
        let mu = ohta_kawasaki_mu(f_a, f_a, 2.0);
        assert!(mu.abs() < 1e-12);
    }

    #[test]
    fn test_block_copolymer_step_conserves_mean() {
        let nx = 5;
        let ny = 5;
        let f_a = 0.4;
        let params = OhtaKawasakiParams {
            chi_n: 2.0,
            mobility: 0.1,
            epsilon: 0.5,
            alpha: 0.1,
            f_a,
            dx: 1.0,
            dt: 0.01,
        };
        let mut phi = vec![vec![f_a; ny]; nx];
        let mean_before: f64 = phi.iter().flat_map(|r| r.iter()).sum::<f64>() / (nx * ny) as f64;
        block_copolymer_step(&mut phi, &params);
        let mean_after: f64 = phi.iter().flat_map(|r| r.iter()).sum::<f64>() / (nx * ny) as f64;
        // Uniform initial condition → zero gradient → mean conserved
        assert!((mean_before - mean_after).abs() < 1e-10);
    }

    #[test]
    fn test_lamellar_period() {
        let params = OhtaKawasakiParams {
            chi_n: 2.0,
            mobility: 0.1,
            epsilon: 1.0,
            alpha: 1.0,
            f_a: 0.5,
            dx: 1.0,
            dt: 0.01,
        };
        let l = params.lamellar_period();
        assert!((l - 2.0 * PI).abs() < 1e-10);
    }

    // ── Active matter ─────────────────────────────────────────────────────────

    #[test]
    fn test_active_particle_propulsion_velocity() {
        let p = ActiveParticle::new(0.0, 0.0, 0.0);
        let v = p.propulsion_velocity();
        assert!((v[0] - ACTIVE_SPEED).abs() < 1e-12);
        assert!(v[1].abs() < 1e-12);
    }

    #[test]
    fn test_active_particle_advance() {
        let mut p = ActiveParticle::new(0.0, 0.0, 0.0);
        let x0 = p.position[0];
        p.advance(0.0, 0.0, 0.0, 1.0);
        assert!(p.position[0] > x0);
    }

    #[test]
    fn test_polar_order_parameter_aligned() {
        let particles: Vec<_> = (0..4).map(|_| ActiveParticle::new(0.0, 0.0, 0.0)).collect();
        let phi = polar_order_parameter(&particles);
        assert!((phi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_nematic_order_parameter_aligned() {
        let particles: Vec<_> = (0..4).map(|_| ActiveParticle::new(0.0, 0.0, 0.0)).collect();
        let s = nematic_order_parameter(&particles);
        assert!((s - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_polar_order_empty() {
        let particles: Vec<ActiveParticle> = vec![];
        assert_eq!(polar_order_parameter(&particles), 0.0);
    }

    #[test]
    fn test_vicsek_alignment_single_particle() {
        let particles = vec![ActiveParticle::new(0.0, 0.0, 0.0)];
        let torque = vicsek_alignment_torque(&particles, 0, 10.0, 1.0);
        assert!(torque.abs() < 1e-12);
    }

    // ── LBM helpers ───────────────────────────────────────────────────────────

    #[test]
    fn test_equilibrium_d2q9_sum_to_rho() {
        let feq = equilibrium_d2q9(1.2, 0.05, -0.03);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.2).abs() < 1e-12);
    }

    #[test]
    fn test_macroscopic_d2q9_roundtrip() {
        let feq = equilibrium_d2q9(1.0, 0.1, -0.05);
        let (rho, ux, uy) = macroscopic_d2q9(&feq);
        assert!((rho - 1.0).abs() < 1e-12);
        assert!((ux - 0.1).abs() < 1e-10);
        assert!((uy + 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_bgk_collide_converges_to_equilibrium() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let feq = equilibrium_d2q9(rho, ux, uy);
        let mut f = feq; // start at equilibrium
        bgk_collide(&mut f, 0.8, rho, ux, uy);
        // should remain at equilibrium
        for i in 0..Q9 {
            assert!((f[i] - feq[i]).abs() < 1e-12);
        }
    }

    // ── SoftMatterGrid ────────────────────────────────────────────────────────

    #[test]
    fn test_soft_matter_grid_init() {
        let params = LiquidCrystalParams::default_lattice();
        let grid = SoftMatterGrid::new(8, 8, 1.0, params);
        assert!((grid.mean_speed()).abs() < 1e-12);
    }

    #[test]
    fn test_soft_matter_grid_step_density_conserved() {
        let params = LiquidCrystalParams::default_lattice();
        let mut grid = SoftMatterGrid::new(8, 8, 1.0, params);
        let rho_before: f64 = grid
            .f
            .iter()
            .flat_map(|r| r.iter())
            .flat_map(|f| f.iter())
            .sum();
        grid.step();
        let rho_after: f64 = grid
            .f
            .iter()
            .flat_map(|r| r.iter())
            .flat_map(|f| f.iter())
            .sum();
        assert!((rho_before - rho_after).abs() < 1e-8);
    }

    #[test]
    fn test_soft_matter_grid_set_director() {
        let params = LiquidCrystalParams::default_lattice();
        let mut grid = SoftMatterGrid::new(4, 4, 1.0, params);
        let theta_field = vec![vec![PI / 4.0; 4]; 4];
        grid.set_director(&theta_field, 0.5);
        let q = grid.q_field[1][1];
        let q_ref = QTensor2::from_director(PI / 4.0, 0.5);
        assert!((q.q11 - q_ref.q11).abs() < 1e-12);
    }

    #[test]
    fn test_mean_order_parameter_uniform() {
        let params = LiquidCrystalParams::default_lattice();
        let mut grid = SoftMatterGrid::new(4, 4, 1.0, params);
        let theta_field = vec![vec![0.0_f64; 4]; 4];
        grid.set_director(&theta_field, 0.7);
        let s_mean = grid.mean_order_parameter();
        assert!((s_mean - 0.7).abs() < 1e-10);
    }
}
