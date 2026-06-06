// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluid-structure interaction (FSI) for the FEM crate.
//!
//! Implements:
//! - [`FsiInterface`] — interface tracking between fluid and structural domains
//! - [`AleMapping`] — Arbitrary Lagrangian-Eulerian mesh motion
//! - [`MeshVelocity`] — ALE mesh velocity from structural displacement
//! - [`AddedMass`] — fluid added mass matrix for submerged structures
//! - [`HydroelasticPressure`] — pressure load on structure from fluid
//! - [`FsiConvergenceMonitor`] — FSI coupling convergence monitoring
//! - [`PartitionedFsi`] — partitioned (staggered) coupling strategy
//! - [`MonolithicFsi`] — placeholder for monolithic FSI assembly
//! - [`StructuralDamping`] — Rayleigh proportional damping C = αM + βK
//! - [`VortexInducedVibration`] — VIV lock-in criterion and amplitude model

use std::f64::consts::PI;

// ── Math helpers ─────────────────────────────────────────────────────────────

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn scale_vec(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// § 1  FsiInterface
// ═══════════════════════════════════════════════════════════════════════════

/// A node on the fluid-structure interface.
#[derive(Debug, Clone)]
pub struct InterfaceNode {
    /// Node index in global mesh.
    pub global_id: usize,
    /// Spatial position \[x, y, z\].
    pub position: [f64; 3],
    /// Outward normal at this interface node.
    pub normal: [f64; 3],
    /// Area weight associated with this node \[m²\].
    pub area_weight: f64,
}

/// Tracks the interface between a fluid domain and a structural domain.
///
/// The interface stores the set of nodes/faces that are common to both domains
/// and provides utilities for interpolation and load/displacement transfer.
#[derive(Debug, Clone)]
pub struct FsiInterface {
    /// Interface nodes on the structural side.
    pub structural_nodes: Vec<InterfaceNode>,
    /// Interface nodes on the fluid side (may differ if non-conforming mesh).
    pub fluid_nodes: Vec<InterfaceNode>,
    /// Whether the interface meshes are conforming (matching node-to-node).
    pub conforming: bool,
}

impl FsiInterface {
    /// Create a new FSI interface.
    pub fn new(conforming: bool) -> Self {
        Self {
            structural_nodes: Vec::new(),
            fluid_nodes: Vec::new(),
            conforming,
        }
    }

    /// Add a structural interface node.
    pub fn add_structural_node(&mut self, node: InterfaceNode) {
        self.structural_nodes.push(node);
    }

    /// Add a fluid interface node.
    pub fn add_fluid_node(&mut self, node: InterfaceNode) {
        self.fluid_nodes.push(node);
    }

    /// Total integrated area of the structural interface \[m²\].
    pub fn total_area(&self) -> f64 {
        self.structural_nodes.iter().map(|n| n.area_weight).sum()
    }

    /// Transfer displacements from structural nodes to fluid nodes.
    ///
    /// Uses nearest-neighbour interpolation for conforming meshes.
    /// Returns a flat vector of size `3 * fluid_nodes.len()`.
    pub fn transfer_displacements(&self, structural_disp: &[f64]) -> Vec<f64> {
        let ns = self.structural_nodes.len();
        if ns == 0 {
            return vec![0.0; 3 * self.fluid_nodes.len()];
        }
        let mut result = Vec::with_capacity(3 * self.fluid_nodes.len());
        for fnode in &self.fluid_nodes {
            // Nearest-neighbour
            let mut best = 0usize;
            let mut best_d = f64::MAX;
            for (si, snode) in self.structural_nodes.iter().enumerate() {
                let dx = fnode.position[0] - snode.position[0];
                let dy = fnode.position[1] - snode.position[1];
                let dz = fnode.position[2] - snode.position[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < best_d {
                    best_d = d;
                    best = si;
                }
            }
            let base = best * 3;
            if base + 2 < structural_disp.len() {
                result.push(structural_disp[base]);
                result.push(structural_disp[base + 1]);
                result.push(structural_disp[base + 2]);
            } else {
                result.extend_from_slice(&[0.0, 0.0, 0.0]);
            }
        }
        result
    }

    /// Transfer traction (pressure × normal) from fluid nodes to structural load vector.
    ///
    /// Returns a flat vector of size `3 * structural_nodes.len()`.
    pub fn transfer_traction(&self, fluid_pressure: &[f64]) -> Vec<f64> {
        let nf = self.fluid_nodes.len();
        let mut loads = vec![0.0f64; 3 * self.structural_nodes.len()];
        for (si, snode) in self.structural_nodes.iter().enumerate() {
            // Nearest fluid neighbour
            let mut best = 0usize;
            let mut best_d = f64::MAX;
            for (fi, fnode) in self.fluid_nodes.iter().enumerate() {
                let dx = snode.position[0] - fnode.position[0];
                let dy = snode.position[1] - fnode.position[1];
                let dz = snode.position[2] - fnode.position[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < best_d {
                    best_d = d;
                    best = fi;
                }
            }
            let p = if best < nf {
                fluid_pressure[best.min(fluid_pressure.len().saturating_sub(1))]
            } else {
                0.0
            };
            let base = si * 3;
            loads[base] -= p * snode.normal[0] * snode.area_weight;
            loads[base + 1] -= p * snode.normal[1] * snode.area_weight;
            loads[base + 2] -= p * snode.normal[2] * snode.area_weight;
        }
        loads
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 2  AleMapping
// ═══════════════════════════════════════════════════════════════════════════

/// Blending strategy for ALE mesh motion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AleBlend {
    /// Purely Lagrangian: mesh moves with the material.
    Lagrangian,
    /// Purely Eulerian: mesh is fixed in space.
    Eulerian,
    /// Linear blend: mesh velocity = θ * material velocity.
    Linear(f64),
}

/// Arbitrary Lagrangian-Eulerian (ALE) mesh motion descriptor.
///
/// Stores the mesh node positions and provides blended mesh velocity
/// between Lagrangian (θ=1) and Eulerian (θ=0) limits.
#[derive(Debug, Clone)]
pub struct AleMapping {
    /// Reference positions of mesh nodes (ndof × 3, flattened row-major).
    pub ref_positions: Vec<f64>,
    /// Current positions of mesh nodes.
    pub cur_positions: Vec<f64>,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Blending mode.
    pub blend: AleBlend,
}

impl AleMapping {
    /// Create a new ALE mapping from reference node positions.
    ///
    /// `positions` is a flat array of `[x0, y0, z0, x1, y1, z1, …]` length `3*n`.
    pub fn new(positions: Vec<f64>, blend: AleBlend) -> Self {
        let n_nodes = positions.len() / 3;
        Self {
            cur_positions: positions.clone(),
            ref_positions: positions,
            n_nodes,
            blend,
        }
    }

    /// Update current positions from displacements.
    pub fn update(&mut self, displacements: &[f64]) {
        for ((cur, &ref_pos), &disp) in self
            .cur_positions
            .iter_mut()
            .zip(self.ref_positions.iter())
            .zip(displacements.iter())
        {
            *cur = ref_pos + disp;
        }
    }

    /// Blending factor θ in \[0, 1\].
    pub fn theta(&self) -> f64 {
        match self.blend {
            AleBlend::Lagrangian => 1.0,
            AleBlend::Eulerian => 0.0,
            AleBlend::Linear(t) => t.clamp(0.0, 1.0),
        }
    }

    /// Compute the ALE mesh displacement: θ * structural_displacement.
    pub fn mesh_displacement(&self, structural_disp: &[f64]) -> Vec<f64> {
        let t = self.theta();
        structural_disp.iter().map(|d| t * d).collect()
    }

    /// Convective ALE term: (1 - θ) for scaling convective velocity.
    pub fn convective_factor(&self) -> f64 {
        1.0 - self.theta()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 3  MeshVelocity
// ═══════════════════════════════════════════════════════════════════════════

/// Computes ALE mesh velocity from successive structural displacements.
///
/// w_mesh = (d_{n+1} - d_n) / Δt, scaled by the ALE blending factor.
#[derive(Debug, Clone)]
pub struct MeshVelocity {
    /// Previous displacement vector (flat, length `3*n`).
    pub prev_disp: Vec<f64>,
    /// Current displacement vector.
    pub cur_disp: Vec<f64>,
    /// ALE blending factor θ ∈ \[0, 1\].
    pub theta: f64,
}

impl MeshVelocity {
    /// Create a new `MeshVelocity` tracker for `n_dof` degrees of freedom.
    pub fn new(n_dof: usize, theta: f64) -> Self {
        Self {
            prev_disp: vec![0.0; n_dof],
            cur_disp: vec![0.0; n_dof],
            theta: theta.clamp(0.0, 1.0),
        }
    }

    /// Advance to the next time step with new structural displacement.
    pub fn advance(&mut self, new_disp: Vec<f64>) {
        self.prev_disp = std::mem::replace(&mut self.cur_disp, new_disp);
    }

    /// Compute mesh velocity for a given time step `dt`.
    ///
    /// Returns `θ * (cur - prev) / dt`.
    pub fn velocity(&self, dt: f64) -> Vec<f64> {
        let inv_dt = if dt.abs() > 1e-300 { 1.0 / dt } else { 0.0 };
        self.cur_disp
            .iter()
            .zip(self.prev_disp.iter())
            .map(|(c, p)| self.theta * (c - p) * inv_dt)
            .collect()
    }

    /// Maximum mesh velocity magnitude.
    pub fn max_velocity_norm(&self, dt: f64) -> f64 {
        let v = self.velocity(dt);
        v.chunks(3)
            .map(|c| {
                let s: f64 = c.iter().map(|x| x * x).sum();
                s.sqrt()
            })
            .fold(0.0f64, f64::max)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 4  AddedMass
// ═══════════════════════════════════════════════════════════════════════════

/// Fluid added mass contribution for a submerged structure.
///
/// The added mass matrix is approximated as M_a = C_a * ρ_f * V * I
/// where C_a is the added mass coefficient, ρ_f the fluid density,
/// and V the displaced volume.
#[derive(Debug, Clone)]
pub struct AddedMass {
    /// Added mass coefficient (C_a ≈ 1 for a sphere, 1 for long cylinders).
    pub coefficient: f64,
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Displaced volume \[m³\].
    pub displaced_volume: f64,
}

impl AddedMass {
    /// Create a new `AddedMass` model.
    pub fn new(coefficient: f64, fluid_density: f64, displaced_volume: f64) -> Self {
        Self {
            coefficient,
            fluid_density,
            displaced_volume,
        }
    }

    /// Scalar added mass value \[kg\].
    ///
    /// m_a = C_a * ρ_f * V
    pub fn scalar_mass(&self) -> f64 {
        self.coefficient * self.fluid_density * self.displaced_volume
    }

    /// Assemble a diagonal added mass contribution into a DOF-sized vector.
    ///
    /// Each translational DOF receives an equal share of the added mass.
    /// `n_dof` must be divisible by 3 (3D problem).
    pub fn assemble_diagonal(&self, n_dof: usize) -> Vec<f64> {
        let m_a = self.scalar_mass();
        let n_nodes = if n_dof >= 3 { n_dof / 3 } else { 1 };
        let per_node = if n_nodes > 0 {
            m_a / n_nodes as f64
        } else {
            0.0
        };
        vec![per_node; n_dof]
    }

    /// Effective acceleration with added mass: a_eff = F / (m_struct + m_a).
    pub fn effective_acceleration(&self, force: f64, structural_mass: f64) -> f64 {
        let total = structural_mass + self.scalar_mass();
        if total.abs() < 1e-300 {
            0.0
        } else {
            force / total
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 5  HydroelasticPressure
// ═══════════════════════════════════════════════════════════════════════════

/// Assembles the hydrostatic or hydrodynamic pressure load vector on a structure.
///
/// For each interface face the distributed pressure p is integrated against
/// the shape functions to produce a consistent nodal load vector.
#[derive(Debug, Clone)]
pub struct HydroelasticPressure {
    /// Fluid pressure values at interface nodes \[Pa\].
    pub nodal_pressure: Vec<f64>,
    /// Outward normals at interface nodes (flat: \[nx0,ny0,nz0,…\]).
    pub normals: Vec<f64>,
    /// Area weights at interface nodes \[m²\].
    pub area_weights: Vec<f64>,
}

impl HydroelasticPressure {
    /// Create a new `HydroelasticPressure` descriptor.
    pub fn new(nodal_pressure: Vec<f64>, normals: Vec<f64>, area_weights: Vec<f64>) -> Self {
        Self {
            nodal_pressure,
            normals,
            area_weights,
        }
    }

    /// Assemble the nodal load vector (size `3 * n_nodes`).
    ///
    /// f_i = −p_i · n_i · A_i  (pressure acts inward on the structure)
    pub fn load_vector(&self) -> Vec<f64> {
        let n = self.nodal_pressure.len();
        let mut f = vec![0.0f64; 3 * n];
        for i in 0..n {
            let p = self.nodal_pressure[i];
            let aw = if i < self.area_weights.len() {
                self.area_weights[i]
            } else {
                1.0
            };
            for d in 0..3 {
                let ni = if 3 * i + d < self.normals.len() {
                    self.normals[3 * i + d]
                } else {
                    0.0
                };
                f[3 * i + d] = -p * ni * aw;
            }
        }
        f
    }

    /// Total resultant force vector \[Fx, Fy, Fz\].
    pub fn resultant_force(&self) -> [f64; 3] {
        let f = self.load_vector();
        let mut res = [0.0f64; 3];
        for chunk in f.chunks(3) {
            if chunk.len() == 3 {
                res[0] += chunk[0];
                res[1] += chunk[1];
                res[2] += chunk[2];
            }
        }
        res
    }

    /// RMS pressure over the interface.
    pub fn rms_pressure(&self) -> f64 {
        let n = self.nodal_pressure.len();
        if n == 0 {
            return 0.0;
        }
        let sum_sq: f64 = self.nodal_pressure.iter().map(|p| p * p).sum();
        (sum_sq / n as f64).sqrt()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 6  FsiConvergenceMonitor
// ═══════════════════════════════════════════════════════════════════════════

/// Monitors convergence of a partitioned FSI coupling iteration.
///
/// Tracks displacement residuals and pressure residuals across coupling cycles.
#[derive(Debug, Clone)]
pub struct FsiConvergenceMonitor {
    /// Tolerance for displacement increment norm.
    pub disp_tol: f64,
    /// Tolerance for pressure residual norm.
    pub pressure_tol: f64,
    /// Maximum number of coupling iterations.
    pub max_iter: usize,
    /// History of displacement residual norms.
    pub disp_history: Vec<f64>,
    /// History of pressure residual norms.
    pub pressure_history: Vec<f64>,
}

impl FsiConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(disp_tol: f64, pressure_tol: f64, max_iter: usize) -> Self {
        Self {
            disp_tol,
            pressure_tol,
            max_iter,
            disp_history: Vec::new(),
            pressure_history: Vec::new(),
        }
    }

    /// Record one coupling iteration's residuals.
    pub fn record(&mut self, disp_residual: f64, pressure_residual: f64) {
        self.disp_history.push(disp_residual);
        self.pressure_history.push(pressure_residual);
    }

    /// Check if convergence is achieved.
    pub fn converged(&self) -> bool {
        if let (Some(&dr), Some(&pr)) = (self.disp_history.last(), self.pressure_history.last()) {
            dr < self.disp_tol && pr < self.pressure_tol
        } else {
            false
        }
    }

    /// Check if maximum iterations exceeded.
    pub fn max_iter_reached(&self) -> bool {
        self.disp_history.len() >= self.max_iter
    }

    /// Compute displacement convergence rate (last two iterations).
    pub fn disp_convergence_rate(&self) -> Option<f64> {
        let n = self.disp_history.len();
        if n >= 2 {
            let prev = self.disp_history[n - 2];
            let curr = self.disp_history[n - 1];
            if prev.abs() > 1e-300 {
                Some(curr / prev)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Reset monitor for a new time step.
    pub fn reset(&mut self) {
        self.disp_history.clear();
        self.pressure_history.clear();
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 7  PartitionedFsi
// ═══════════════════════════════════════════════════════════════════════════

/// Coupling scheme for partitioned FSI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CouplingScheme {
    /// Sequential (one-way): fluid → structure, no feedback.
    Sequential,
    /// Staggered (loosely coupled): fluid → structure → fluid, fixed point.
    Staggered,
    /// Strong coupling: iterate until convergence within each time step.
    StrongCoupling,
}

/// Partitioned FSI solver that alternates between fluid and structural sub-solvers.
///
/// Uses Aitken relaxation when `use_aitken = true` to accelerate fixed-point
/// convergence.
#[derive(Debug, Clone)]
pub struct PartitionedFsi {
    /// FSI interface descriptor.
    pub interface: FsiInterface,
    /// Convergence monitor.
    pub monitor: FsiConvergenceMonitor,
    /// Coupling scheme.
    pub scheme: CouplingScheme,
    /// Aitken relaxation factor ω ∈ (0, 1].
    pub omega: f64,
    /// Enable Aitken dynamic relaxation.
    pub use_aitken: bool,
    /// Previous displacement increment (for Aitken update).
    pub prev_delta: Vec<f64>,
}

impl PartitionedFsi {
    /// Create a new partitioned FSI solver.
    pub fn new(
        interface: FsiInterface,
        scheme: CouplingScheme,
        disp_tol: f64,
        pressure_tol: f64,
        max_iter: usize,
    ) -> Self {
        Self {
            interface,
            monitor: FsiConvergenceMonitor::new(disp_tol, pressure_tol, max_iter),
            scheme,
            omega: 0.5,
            use_aitken: false,
            prev_delta: Vec::new(),
        }
    }

    /// Apply Aitken relaxation to the displacement increment.
    ///
    /// Returns the relaxed increment and updates `omega`.
    pub fn aitken_relax(&mut self, delta: &[f64]) -> Vec<f64> {
        if self.prev_delta.is_empty() || self.prev_delta.len() != delta.len() {
            self.prev_delta = delta.to_vec();
            return scale_vec(delta, self.omega);
        }
        // ω_{n+1} = -ω_n * (Δd_prev · (Δd - Δd_prev)) / ‖Δd - Δd_prev‖²
        let diff: Vec<f64> = delta
            .iter()
            .zip(self.prev_delta.iter())
            .map(|(a, b)| a - b)
            .collect();
        let denom = dot(&diff, &diff);
        if denom > 1e-300 {
            let num = dot(&self.prev_delta, &diff);
            self.omega = (-self.omega * num / denom).clamp(0.1, 1.0);
        }
        self.prev_delta = delta.to_vec();
        scale_vec(delta, self.omega)
    }

    /// Check if the current scheme requires inner iterations.
    pub fn requires_inner_loop(&self) -> bool {
        matches!(self.scheme, CouplingScheme::StrongCoupling)
    }

    /// Record one sub-iteration and check convergence.
    pub fn step(&mut self, disp_res: f64, pres_res: f64) -> bool {
        self.monitor.record(disp_res, pres_res);
        self.monitor.converged() || self.monitor.max_iter_reached()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 8  MonolithicFsi
// ═══════════════════════════════════════════════════════════════════════════

/// Placeholder for a monolithic FSI system.
///
/// In a monolithic approach, the coupled fluid-structure system is assembled
/// into a single matrix and solved simultaneously. This struct holds the
/// block structure metadata needed for that assembly.
#[derive(Debug, Clone)]
pub struct MonolithicFsi {
    /// Number of structural DOFs.
    pub n_structural_dof: usize,
    /// Number of fluid pressure DOFs.
    pub n_fluid_pressure_dof: usize,
    /// Number of fluid velocity DOFs.
    pub n_fluid_velocity_dof: usize,
    /// Stabilisation parameter τ for SUPG/GLS.
    pub tau_stab: f64,
}

impl MonolithicFsi {
    /// Create a new `MonolithicFsi` descriptor.
    pub fn new(
        n_structural_dof: usize,
        n_fluid_pressure_dof: usize,
        n_fluid_velocity_dof: usize,
        tau_stab: f64,
    ) -> Self {
        Self {
            n_structural_dof,
            n_fluid_pressure_dof,
            n_fluid_velocity_dof,
            tau_stab,
        }
    }

    /// Total number of DOFs in the monolithic system.
    pub fn total_dof(&self) -> usize {
        self.n_structural_dof + self.n_fluid_pressure_dof + self.n_fluid_velocity_dof
    }

    /// Return the block offsets: (struct_start, p_start, v_start).
    pub fn block_offsets(&self) -> (usize, usize, usize) {
        let s = 0;
        let p = s + self.n_structural_dof;
        let v = p + self.n_fluid_pressure_dof;
        (s, p, v)
    }

    /// Assemble a zero flat matrix (row-major) of size `n × n`.
    ///
    /// This is the placeholder for the full monolithic tangent matrix.
    pub fn zero_matrix(&self) -> Vec<f64> {
        let n = self.total_dof();
        vec![0.0f64; n * n]
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 9  StructuralDamping
// ═══════════════════════════════════════════════════════════════════════════

/// Rayleigh proportional damping: C = α·M + β·K.
///
/// The coefficients α and β can be determined from two modal frequencies
/// and target damping ratios.
#[derive(Debug, Clone)]
pub struct StructuralDamping {
    /// Mass-proportional damping coefficient α \[1/s\].
    pub alpha: f64,
    /// Stiffness-proportional damping coefficient β \[s\].
    pub beta: f64,
}

impl StructuralDamping {
    /// Create a new `StructuralDamping` with given Rayleigh coefficients.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    /// Construct Rayleigh coefficients from two frequencies (rad/s) and
    /// corresponding damping ratios ζ₁, ζ₂.
    ///
    /// α = 2·ω₁·ω₂·(ζ₁·ω₂ − ζ₂·ω₁) / (ω₂² − ω₁²)
    /// β = 2·(ζ₂·ω₂ − ζ₁·ω₁) / (ω₂² − ω₁²)
    pub fn from_frequencies(omega1: f64, omega2: f64, zeta1: f64, zeta2: f64) -> Self {
        let denom = omega2 * omega2 - omega1 * omega1;
        if denom.abs() < 1e-300 {
            return Self::new(0.0, 0.0);
        }
        let alpha = 2.0 * omega1 * omega2 * (zeta1 * omega2 - zeta2 * omega1) / denom;
        let beta = 2.0 * (zeta2 * omega2 - zeta1 * omega1) / denom;
        Self { alpha, beta }
    }

    /// Effective damping ratio at frequency ω: ζ = α/(2ω) + βω/2.
    pub fn damping_ratio(&self, omega: f64) -> f64 {
        if omega.abs() < 1e-300 {
            return 0.0;
        }
        self.alpha / (2.0 * omega) + self.beta * omega / 2.0
    }

    /// Apply damping to produce C·v: C·v = α·M·v + β·K·v.
    ///
    /// `mv` = M·v (pre-computed), `kv` = K·v (pre-computed).
    pub fn damping_force(&self, mv: &[f64], kv: &[f64]) -> Vec<f64> {
        mv.iter()
            .zip(kv.iter())
            .map(|(m, k)| self.alpha * m + self.beta * k)
            .collect()
    }

    /// Assemble damping matrix C = α·M + β·K (flat row-major).
    ///
    /// `mass_mat` and `stiff_mat` are flat row-major matrices of size `n×n`.
    pub fn assemble_matrix(&self, mass_mat: &[f64], stiff_mat: &[f64]) -> Vec<f64> {
        mass_mat
            .iter()
            .zip(stiff_mat.iter())
            .map(|(m, k)| self.alpha * m + self.beta * k)
            .collect()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 10  VortexInducedVibration
// ═══════════════════════════════════════════════════════════════════════════

/// Vortex-induced vibration (VIV) model for a circular cylinder.
///
/// Provides:
/// - Reduced velocity and Strouhal-based vortex shedding frequency
/// - Lock-in criterion
/// - Hartlen-Currie wake oscillator model for amplitude estimation
#[derive(Debug, Clone)]
pub struct VortexInducedVibration {
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Natural frequency of the structure \[Hz\].
    pub fn_struct: f64,
    /// Strouhal number (≈ 0.2 for circular cylinder).
    pub strouhal: f64,
    /// Structural damping ratio ζ.
    pub zeta: f64,
    /// Mass ratio m* = m / (ρ_f * π/4 * D² * L).
    pub mass_ratio: f64,
}

impl VortexInducedVibration {
    /// Create a new VIV model.
    pub fn new(diameter: f64, fn_struct: f64, strouhal: f64, zeta: f64, mass_ratio: f64) -> Self {
        Self {
            diameter,
            fn_struct,
            strouhal,
            zeta,
            mass_ratio,
        }
    }

    /// Reduced velocity: U_r = U / (f_n * D).
    pub fn reduced_velocity(&self, flow_velocity: f64) -> f64 {
        if self.fn_struct.abs() < 1e-300 || self.diameter.abs() < 1e-300 {
            return 0.0;
        }
        flow_velocity / (self.fn_struct * self.diameter)
    }

    /// Vortex shedding frequency \[Hz\]: f_vs = St * U / D.
    pub fn shedding_frequency(&self, flow_velocity: f64) -> f64 {
        if self.diameter.abs() < 1e-300 {
            return 0.0;
        }
        self.strouhal * flow_velocity / self.diameter
    }

    /// Check lock-in: |f_vs - f_n| / f_n < threshold (default 0.1).
    pub fn is_lock_in(&self, flow_velocity: f64) -> bool {
        let f_vs = self.shedding_frequency(flow_velocity);
        if self.fn_struct.abs() < 1e-300 {
            return false;
        }
        (f_vs - self.fn_struct).abs() / self.fn_struct < 0.1
    }

    /// Synchronisation range of reduced velocity (U_r_lower, U_r_upper).
    ///
    /// Typically \[1/(St*(1+Δ)), 1/(St*(1-Δ))\] with Δ = 0.1.
    pub fn lock_in_range(&self) -> (f64, f64) {
        let st = self.strouhal;
        if st.abs() < 1e-300 {
            return (0.0, 0.0);
        }
        (1.0 / (st * 1.1), 1.0 / (st * 0.9))
    }

    /// Peak VIV amplitude estimate (Griffin plot approximation).
    ///
    /// A/D ≈ 1.29 * exp(−0.41 * Sg^0.36), where Sg = 2π * ζ * m*.
    pub fn peak_amplitude(&self) -> f64 {
        let sg = 2.0 * PI * self.zeta * self.mass_ratio;
        1.29 * (-0.41 * sg.powf(0.36)).exp()
    }

    /// Scruton number (mass-damping parameter): Sc = 2π * m* * ζ.
    pub fn scruton_number(&self) -> f64 {
        2.0 * PI * self.mass_ratio * self.zeta
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── FsiInterface ────────────────────────────────────────────────────────

    fn make_node(id: usize, x: f64, y: f64, z: f64) -> InterfaceNode {
        InterfaceNode {
            global_id: id,
            position: [x, y, z],
            normal: [0.0, 0.0, 1.0],
            area_weight: 0.25,
        }
    }

    #[test]
    fn test_interface_total_area() {
        let mut iface = FsiInterface::new(true);
        iface.add_structural_node(make_node(0, 0.0, 0.0, 0.0));
        iface.add_structural_node(make_node(1, 1.0, 0.0, 0.0));
        assert!((iface.total_area() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_interface_empty_area() {
        let iface = FsiInterface::new(false);
        assert_eq!(iface.total_area(), 0.0);
    }

    #[test]
    fn test_transfer_displacements_conforming() {
        let mut iface = FsiInterface::new(true);
        iface.add_structural_node(make_node(0, 0.0, 0.0, 0.0));
        iface.add_structural_node(make_node(1, 1.0, 0.0, 0.0));
        iface.add_fluid_node(make_node(0, 0.05, 0.0, 0.0)); // close to node 0
        let disp = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let result = iface.transfer_displacements(&disp);
        assert_eq!(result.len(), 3);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
        assert!((result[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_transfer_displacements_empty() {
        let iface = FsiInterface::new(true);
        let result = iface.transfer_displacements(&[1.0, 2.0, 3.0]);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_transfer_traction_zero_pressure() {
        let mut iface = FsiInterface::new(true);
        iface.add_structural_node(make_node(0, 0.0, 0.0, 0.0));
        iface.add_fluid_node(make_node(0, 0.0, 0.0, 0.0));
        let loads = iface.transfer_traction(&[0.0]);
        assert_eq!(loads.len(), 3);
        assert!(loads.iter().all(|&x| x.abs() < 1e-15));
    }

    #[test]
    fn test_transfer_traction_unit_pressure() {
        let mut iface = FsiInterface::new(true);
        let mut snode = make_node(0, 0.0, 0.0, 0.0);
        snode.normal = [0.0, 0.0, 1.0];
        snode.area_weight = 1.0;
        iface.add_structural_node(snode);
        iface.add_fluid_node(make_node(0, 0.0, 0.0, 0.0));
        let loads = iface.transfer_traction(&[1000.0]);
        // f_z = -p * n_z * A = -1000 * 1 * 1 = -1000
        assert!((loads[2] - (-1000.0)).abs() < 1e-9);
    }

    // ── AleMapping ──────────────────────────────────────────────────────────

    #[test]
    fn test_ale_theta_lagrangian() {
        let ale = AleMapping::new(vec![0.0; 6], AleBlend::Lagrangian);
        assert_eq!(ale.theta(), 1.0);
    }

    #[test]
    fn test_ale_theta_eulerian() {
        let ale = AleMapping::new(vec![0.0; 6], AleBlend::Eulerian);
        assert_eq!(ale.theta(), 0.0);
    }

    #[test]
    fn test_ale_theta_linear() {
        let ale = AleMapping::new(vec![0.0; 6], AleBlend::Linear(0.7));
        assert!((ale.theta() - 0.7).abs() < 1e-12);
    }

    #[test]
    fn test_ale_mesh_displacement_lagrangian() {
        let ale = AleMapping::new(vec![0.0; 6], AleBlend::Lagrangian);
        let d = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let md = ale.mesh_displacement(&d);
        assert_eq!(md, d);
    }

    #[test]
    fn test_ale_mesh_displacement_eulerian() {
        let ale = AleMapping::new(vec![0.0; 6], AleBlend::Eulerian);
        let d = vec![1.0, 2.0, 3.0];
        let md = ale.mesh_displacement(&d);
        assert!(md.iter().all(|&x| x.abs() < 1e-15));
    }

    #[test]
    fn test_ale_update() {
        let mut ale = AleMapping::new(vec![0.0, 0.0, 0.0], AleBlend::Lagrangian);
        ale.update(&[1.0, 2.0, 3.0]);
        assert!((ale.cur_positions[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_ale_convective_factor() {
        let ale = AleMapping::new(vec![], AleBlend::Linear(0.3));
        assert!((ale.convective_factor() - 0.7).abs() < 1e-12);
    }

    // ── MeshVelocity ────────────────────────────────────────────────────────

    #[test]
    fn test_mesh_velocity_zero_initial() {
        let mv = MeshVelocity::new(6, 1.0);
        let v = mv.velocity(0.01);
        assert!(v.iter().all(|&x| x.abs() < 1e-15));
    }

    #[test]
    fn test_mesh_velocity_after_advance() {
        let mut mv = MeshVelocity::new(3, 1.0);
        mv.advance(vec![0.01, 0.02, 0.03]);
        let v = mv.velocity(0.01);
        assert!((v[0] - 1.0).abs() < 1e-9);
        assert!((v[1] - 2.0).abs() < 1e-9);
        assert!((v[2] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_mesh_velocity_theta_half() {
        let mut mv = MeshVelocity::new(3, 0.5);
        mv.advance(vec![0.02, 0.0, 0.0]);
        let v = mv.velocity(0.01);
        assert!((v[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_mesh_velocity_zero_dt() {
        let mv = MeshVelocity::new(3, 1.0);
        let v = mv.velocity(0.0);
        assert!(v.iter().all(|&x| x.abs() < 1e-15));
    }

    #[test]
    fn test_mesh_velocity_max_norm() {
        let mut mv = MeshVelocity::new(3, 1.0);
        mv.advance(vec![0.03, 0.04, 0.0]);
        let vmax = mv.max_velocity_norm(0.01);
        assert!((vmax - 5.0).abs() < 1e-9);
    }

    // ── AddedMass ───────────────────────────────────────────────────────────

    #[test]
    fn test_added_mass_scalar() {
        let am = AddedMass::new(1.0, 1000.0, 0.5);
        assert!((am.scalar_mass() - 500.0).abs() < 1e-9);
    }

    #[test]
    fn test_added_mass_coefficient() {
        let am = AddedMass::new(0.5, 1000.0, 1.0);
        assert!((am.scalar_mass() - 500.0).abs() < 1e-9);
    }

    #[test]
    fn test_added_mass_diagonal_sum() {
        // n_dof=9 → n_nodes=3; per_node = m_a/3; sum = 9 * (m_a/3) = 3*m_a
        let am = AddedMass::new(1.0, 1000.0, 0.3);
        let d = am.assemble_diagonal(9);
        let total: f64 = d.iter().sum();
        let n_nodes = 3usize;
        let expected = am.scalar_mass() / n_nodes as f64 * 9.0;
        assert!((total - expected).abs() < 1e-6);
    }

    #[test]
    fn test_added_mass_effective_acceleration() {
        let am = AddedMass::new(1.0, 1000.0, 1.0);
        let a = am.effective_acceleration(2000.0, 1000.0);
        assert!((a - 1.0).abs() < 1e-9); // F / (m_struct + m_a) = 2000 / 2000
    }

    #[test]
    fn test_added_mass_zero_density() {
        let am = AddedMass::new(1.0, 0.0, 1.0);
        assert_eq!(am.scalar_mass(), 0.0);
    }

    // ── HydroelasticPressure ────────────────────────────────────────────────

    #[test]
    fn test_hydroelastic_load_vector_length() {
        let hp = HydroelasticPressure::new(
            vec![100.0, 200.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            vec![0.5, 0.5],
        );
        let f = hp.load_vector();
        assert_eq!(f.len(), 6);
    }

    #[test]
    fn test_hydroelastic_load_vector_values() {
        let hp = HydroelasticPressure::new(vec![1000.0], vec![0.0, 0.0, 1.0], vec![1.0]);
        let f = hp.load_vector();
        assert!((f[2] - (-1000.0)).abs() < 1e-9);
    }

    #[test]
    fn test_hydroelastic_resultant_force() {
        let hp = HydroelasticPressure::new(
            vec![500.0, 500.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            vec![1.0, 1.0],
        );
        let res = hp.resultant_force();
        assert!((res[2] - (-1000.0)).abs() < 1e-9);
    }

    #[test]
    fn test_hydroelastic_rms_pressure() {
        let hp = HydroelasticPressure::new(vec![3.0, 4.0], vec![0.0; 6], vec![1.0, 1.0]);
        let rms = hp.rms_pressure();
        // sqrt((9+16)/2) = sqrt(12.5)
        assert!((rms - (12.5f64).sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_hydroelastic_empty() {
        let hp = HydroelasticPressure::new(vec![], vec![], vec![]);
        assert_eq!(hp.rms_pressure(), 0.0);
        assert_eq!(hp.load_vector().len(), 0);
    }

    // ── FsiConvergenceMonitor ───────────────────────────────────────────────

    #[test]
    fn test_monitor_not_converged_initially() {
        let monitor = FsiConvergenceMonitor::new(1e-6, 1e-6, 10);
        assert!(!monitor.converged());
    }

    #[test]
    fn test_monitor_converged() {
        let mut monitor = FsiConvergenceMonitor::new(1e-4, 1e-4, 10);
        monitor.record(1e-5, 1e-5);
        assert!(monitor.converged());
    }

    #[test]
    fn test_monitor_not_converged() {
        let mut monitor = FsiConvergenceMonitor::new(1e-6, 1e-6, 10);
        monitor.record(1e-3, 1e-3);
        assert!(!monitor.converged());
    }

    #[test]
    fn test_monitor_max_iter() {
        let mut monitor = FsiConvergenceMonitor::new(1e-12, 1e-12, 3);
        for _ in 0..3 {
            monitor.record(1.0, 1.0);
        }
        assert!(monitor.max_iter_reached());
    }

    #[test]
    fn test_monitor_convergence_rate() {
        let mut monitor = FsiConvergenceMonitor::new(1e-6, 1e-6, 10);
        monitor.record(1.0, 0.5);
        monitor.record(0.5, 0.25);
        let rate = monitor.disp_convergence_rate().unwrap();
        assert!((rate - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_monitor_reset() {
        let mut monitor = FsiConvergenceMonitor::new(1e-6, 1e-6, 10);
        monitor.record(1.0, 1.0);
        monitor.reset();
        assert!(monitor.disp_history.is_empty());
        assert!(!monitor.converged());
    }

    // ── PartitionedFsi ──────────────────────────────────────────────────────

    #[test]
    fn test_partitioned_fsi_requires_inner_loop() {
        let iface = FsiInterface::new(true);
        let fsi = PartitionedFsi::new(
            iface.clone(),
            CouplingScheme::StrongCoupling,
            1e-6,
            1e-6,
            10,
        );
        assert!(fsi.requires_inner_loop());

        let fsi2 = PartitionedFsi::new(iface, CouplingScheme::Sequential, 1e-6, 1e-6, 10);
        assert!(!fsi2.requires_inner_loop());
    }

    #[test]
    fn test_partitioned_fsi_step_convergence() {
        let iface = FsiInterface::new(true);
        let mut fsi = PartitionedFsi::new(iface, CouplingScheme::Staggered, 1e-4, 1e-4, 10);
        let done = fsi.step(1e-5, 1e-5);
        assert!(done); // converged
    }

    #[test]
    fn test_aitken_relax_first_step() {
        let iface = FsiInterface::new(true);
        let mut fsi = PartitionedFsi::new(iface, CouplingScheme::StrongCoupling, 1e-6, 1e-6, 10);
        fsi.omega = 0.5;
        let delta = vec![2.0, 4.0, 6.0];
        let relaxed = fsi.aitken_relax(&delta);
        // First call: just scale by omega
        assert!((relaxed[0] - 1.0).abs() < 1e-12);
        assert!((relaxed[1] - 2.0).abs() < 1e-12);
    }

    // ── MonolithicFsi ───────────────────────────────────────────────────────

    #[test]
    fn test_monolithic_total_dof() {
        let mfsi = MonolithicFsi::new(30, 10, 20, 0.1);
        assert_eq!(mfsi.total_dof(), 60);
    }

    #[test]
    fn test_monolithic_block_offsets() {
        let mfsi = MonolithicFsi::new(30, 10, 20, 0.1);
        let (s, p, v) = mfsi.block_offsets();
        assert_eq!(s, 0);
        assert_eq!(p, 30);
        assert_eq!(v, 40);
    }

    #[test]
    fn test_monolithic_zero_matrix_size() {
        let mfsi = MonolithicFsi::new(5, 3, 2, 0.1);
        let m = mfsi.zero_matrix();
        assert_eq!(m.len(), 100); // 10*10
        assert!(m.iter().all(|&x| x == 0.0));
    }

    // ── StructuralDamping ───────────────────────────────────────────────────

    #[test]
    fn test_rayleigh_damping_from_frequencies() {
        // ω1=10, ω2=100 rad/s, ζ1=ζ2=0.05
        let sd = StructuralDamping::from_frequencies(10.0, 100.0, 0.05, 0.05);
        // Check damping ratios at both frequencies
        let z1 = sd.damping_ratio(10.0);
        let z2 = sd.damping_ratio(100.0);
        assert!((z1 - 0.05).abs() < 1e-6, "z1={}", z1);
        assert!((z2 - 0.05).abs() < 1e-6, "z2={}", z2);
    }

    #[test]
    fn test_rayleigh_damping_force() {
        let sd = StructuralDamping::new(2.0, 0.01);
        let mv = vec![1.0, 2.0, 3.0];
        let kv = vec![4.0, 5.0, 6.0];
        let f = sd.damping_force(&mv, &kv);
        assert!((f[0] - (2.0 * 1.0 + 0.01 * 4.0)).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_damping_assemble_matrix() {
        let sd = StructuralDamping::new(1.0, 1.0);
        let m = vec![1.0, 0.0, 0.0, 1.0];
        let k = vec![2.0, 0.0, 0.0, 2.0];
        let c = sd.assemble_matrix(&m, &k);
        assert!((c[0] - 3.0).abs() < 1e-12);
        assert!((c[3] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_damping_zero_omega() {
        let sd = StructuralDamping::new(2.0, 0.01);
        assert_eq!(sd.damping_ratio(0.0), 0.0);
    }

    #[test]
    fn test_rayleigh_from_equal_frequencies() {
        let sd = StructuralDamping::from_frequencies(10.0, 10.0, 0.05, 0.05);
        // degenerate case: should return zeros
        assert_eq!(sd.alpha, 0.0);
        assert_eq!(sd.beta, 0.0);
    }

    // ── VortexInducedVibration ──────────────────────────────────────────────

    #[test]
    fn test_viv_reduced_velocity() {
        let viv = VortexInducedVibration::new(0.1, 5.0, 0.2, 0.01, 2.0);
        let ur = viv.reduced_velocity(1.0);
        assert!((ur - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_viv_shedding_frequency() {
        let viv = VortexInducedVibration::new(0.1, 5.0, 0.2, 0.01, 2.0);
        let fvs = viv.shedding_frequency(1.0);
        assert!((fvs - 2.0).abs() < 1e-9); // St * U / D = 0.2 * 1 / 0.1
    }

    #[test]
    fn test_viv_lock_in_yes() {
        // f_vs ≈ f_n when U = f_n * D / St
        let viv = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.01, 2.0);
        // U such that f_vs = f_n: U = f_n * D / St = 2.0 * 0.1 / 0.2 = 1.0
        assert!(viv.is_lock_in(1.0));
    }

    #[test]
    fn test_viv_lock_in_no() {
        let viv = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.01, 2.0);
        // At U=10, f_vs = 0.2*10/0.1 = 20 >> f_n = 2
        assert!(!viv.is_lock_in(10.0));
    }

    #[test]
    fn test_viv_lock_in_range() {
        let viv = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.01, 2.0);
        let (lo, hi) = viv.lock_in_range();
        assert!(lo < hi);
        // At U_r = 1/St = 5.0 lock-in should occur
        let mid = (lo + hi) / 2.0;
        assert!(lo < mid && mid < hi);
    }

    #[test]
    fn test_viv_peak_amplitude_positive() {
        let viv = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.01, 2.0);
        let amp = viv.peak_amplitude();
        assert!(amp > 0.0);
    }

    #[test]
    fn test_viv_scruton_number() {
        let viv = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.05, 10.0);
        let sc = viv.scruton_number();
        let expected = 2.0 * PI * 10.0 * 0.05;
        assert!((sc - expected).abs() < 1e-9);
    }

    #[test]
    fn test_viv_zero_diameter() {
        let viv = VortexInducedVibration::new(0.0, 5.0, 0.2, 0.01, 2.0);
        assert_eq!(viv.reduced_velocity(1.0), 0.0);
        assert_eq!(viv.shedding_frequency(1.0), 0.0);
    }

    #[test]
    fn test_viv_high_damping_low_amplitude() {
        let viv_low = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.001, 1.0);
        let viv_high = VortexInducedVibration::new(0.1, 2.0, 0.2, 0.5, 10.0);
        assert!(viv_low.peak_amplitude() > viv_high.peak_amplitude());
    }
}

// ============================================================================
// § NEW: FluidStructureFEM, Morison equation, ALE velocity, IBm, sloshing
// ============================================================================

/// Fluid-structure interaction FEM domain descriptor.
///
/// Combines fluid and solid sub-domains with an interface.
pub struct FluidStructureFEM {
    /// Number of fluid nodes.
    pub fluid_nodes: usize,
    /// Number of solid nodes.
    pub solid_nodes: usize,
    /// Number of interface nodes shared between fluid and solid.
    pub interface_nodes: usize,
    /// Fluid density (kg/m³).
    pub rho_f: f64,
    /// Solid density (kg/m³).
    pub rho_s: f64,
}

impl FluidStructureFEM {
    /// Create a new FSI FEM domain.
    pub fn new(
        fluid_nodes: usize,
        solid_nodes: usize,
        interface_nodes: usize,
        rho_f: f64,
        rho_s: f64,
    ) -> Self {
        Self {
            fluid_nodes,
            solid_nodes,
            interface_nodes,
            rho_f,
            rho_s,
        }
    }

    /// Total nodes across both domains.
    pub fn total_nodes(&self) -> usize {
        self.fluid_nodes + self.solid_nodes
    }

    /// Density ratio ρ_f / ρ_s.
    pub fn density_ratio(&self) -> f64 {
        if self.rho_s.abs() < 1e-300 {
            return 0.0;
        }
        self.rho_f / self.rho_s
    }
}

/// Added mass coefficient C_m for standard geometries.
///
/// Returns:
/// - `"sphere"` → 0.5
/// - `"cylinder"` → 1.0
/// - `"flat_plate"` → π/4
/// - other → 1.0 (default)
pub fn added_mass_coefficient(geometry: &str) -> f64 {
    match geometry {
        "sphere" => 0.5,
        "cylinder" => 1.0,
        "flat_plate" => PI / 4.0,
        _ => 1.0,
    }
}

/// Morison equation: total force per unit length on a cylindrical structure.
///
/// F = 0.5 * C_d * ρ * D * |u| * u  +  C_m * ρ * π/4 * D² * u_dot
///
/// - `cd`: drag coefficient
/// - `cm`: inertia coefficient
/// - `rho`: fluid density (kg/m³)
/// - `d`: diameter (m)
/// - `u`: relative fluid velocity (m/s)
/// - `u_dot`: fluid acceleration (m/s²)
pub fn morison_force(cd: f64, cm: f64, rho: f64, d: f64, u: f64, u_dot: f64) -> f64 {
    let drag = 0.5 * cd * rho * d * u.abs() * u;
    let inertia = cm * rho * PI / 4.0 * d * d * u_dot;
    drag + inertia
}

/// Vortex-induced vibration frequency using the Strouhal relation.
///
/// f_VIV = St * U / D
///
/// - `u`: free-stream velocity (m/s)
/// - `d`: cylinder diameter (m)
/// - `st`: Strouhal number (≈ 0.2 for circular cylinders)
pub fn vortex_induced_vibration_frequency(u: f64, d: f64, st: f64) -> f64 {
    if d.abs() < 1e-300 {
        return 0.0;
    }
    st * u / d
}

/// Check VIV lock-in condition: |f_nat - f_viv| < tol.
pub fn lock_in_condition(f_nat: f64, f_viv: f64, tol: f64) -> bool {
    (f_nat - f_viv).abs() <= tol
}

/// Fluid-structure mass coupling factor: fluid/solid mass ratio.
///
/// coupling = (ρ_f * V_f) / (ρ_s * V_s)
pub fn fluid_structure_coupling_factor(
    rho_f: f64,
    rho_s: f64,
    volume_f: f64,
    volume_s: f64,
) -> f64 {
    let denom = rho_s * volume_s;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    rho_f * volume_f / denom
}

/// ALE convective velocity: relative velocity of fluid with respect to mesh.
///
/// c_ALE = v_fluid − v_mesh
pub fn arbitrary_lagrangian_eulerian_velocity(v_mesh: [f64; 3], v_fluid: [f64; 3]) -> [f64; 3] {
    [
        v_fluid[0] - v_mesh[0],
        v_fluid[1] - v_mesh[1],
        v_fluid[2] - v_mesh[2],
    ]
}

/// FSI interface conditions.
///
/// Encodes the physical conditions enforced at a fluid-structure interface.
pub struct InterfaceCondition {
    /// No-slip condition: fluid velocity equals structural velocity at interface.
    pub no_slip: bool,
    /// Pressure continuity across the interface.
    pub pressure_continuity: bool,
    /// Traction balance: fluid traction equals structural surface traction.
    pub traction_balance: bool,
}

impl InterfaceCondition {
    /// Create a set of standard FSI interface conditions.
    pub fn new(no_slip: bool, pressure_continuity: bool, traction_balance: bool) -> Self {
        Self {
            no_slip,
            pressure_continuity,
            traction_balance,
        }
    }

    /// Create the standard set of FSI conditions (all enabled).
    pub fn standard() -> Self {
        Self::new(true, true, true)
    }

    /// Returns true if all three conditions are satisfied.
    pub fn all_satisfied(&self) -> bool {
        self.no_slip && self.pressure_continuity && self.traction_balance
    }
}

/// Immersed boundary method (IBM) penalty force density.
///
/// f_IBM = k_ib * (u_target - u_fluid)
///
/// - `target_vel`: desired velocity at the immersed boundary (m/s)
/// - `fluid_vel`: current fluid velocity at the boundary (m/s)
/// - `k_ib`: penalty spring constant (N·s/m⁴)
pub fn ibm_force_density(target_vel: [f64; 3], fluid_vel: [f64; 3], k_ib: f64) -> [f64; 3] {
    [
        k_ib * (target_vel[0] - fluid_vel[0]),
        k_ib * (target_vel[1] - fluid_vel[1]),
        k_ib * (target_vel[2] - fluid_vel[2]),
    ]
}

/// Sloshing natural frequency for the n-th mode in a rectangular tank.
///
/// ω_n = sqrt(n π g / L * tanh(n π h / L))
///
/// - `g`: gravitational acceleration (m/s²)
/// - `h`: liquid depth (m)
/// - `l`: tank length (m)
/// - `n`: mode number (1, 2, …)
pub fn sloshing_frequency(g: f64, h: f64, l: f64, n: usize) -> f64 {
    if l.abs() < 1e-300 || n == 0 {
        return 0.0;
    }
    let k = n as f64 * PI / l;
    (k * g * (k * h).tanh()).sqrt()
}

// ============================================================================
// Tests for new FSI functions
// ============================================================================

#[cfg(test)]
mod fsi_new_tests {
    use super::*;

    #[test]
    fn test_fluid_structure_fem_new() {
        let fsi = FluidStructureFEM::new(100, 50, 20, 1000.0, 7800.0);
        assert_eq!(fsi.fluid_nodes, 100);
        assert_eq!(fsi.solid_nodes, 50);
        assert_eq!(fsi.interface_nodes, 20);
    }

    #[test]
    fn test_fluid_structure_fem_total_nodes() {
        let fsi = FluidStructureFEM::new(100, 50, 20, 1000.0, 7800.0);
        assert_eq!(fsi.total_nodes(), 150);
    }

    #[test]
    fn test_fluid_structure_fem_density_ratio() {
        let fsi = FluidStructureFEM::new(10, 10, 5, 1000.0, 2000.0);
        assert!((fsi.density_ratio() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_fluid_structure_fem_zero_solid_density() {
        let fsi = FluidStructureFEM::new(10, 10, 5, 1000.0, 0.0);
        assert_eq!(fsi.density_ratio(), 0.0);
    }

    #[test]
    fn test_added_mass_sphere() {
        assert!((added_mass_coefficient("sphere") - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_added_mass_cylinder() {
        assert!((added_mass_coefficient("cylinder") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_added_mass_flat_plate() {
        let expected = PI / 4.0;
        assert!((added_mass_coefficient("flat_plate") - expected).abs() < 1e-12);
    }

    #[test]
    fn test_added_mass_default() {
        assert!((added_mass_coefficient("unknown") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_morison_force_pure_drag() {
        // No acceleration → only drag term
        let f = morison_force(1.0, 2.0, 1000.0, 0.5, 1.0, 0.0);
        let expected = 0.5 * 1.0 * 1000.0 * 0.5 * 1.0;
        assert!((f - expected).abs() < 1e-10);
    }

    #[test]
    fn test_morison_force_pure_inertia() {
        // No velocity → only inertia term
        let f = morison_force(1.0, 2.0, 1000.0, 0.5, 0.0, 1.0);
        let expected = 2.0 * 1000.0 * PI / 4.0 * 0.25 * 1.0;
        assert!((f - expected).abs() < 1e-8);
    }

    #[test]
    fn test_morison_force_negative_velocity() {
        // Drag should oppose direction of flow
        let f_pos = morison_force(1.0, 0.0, 1000.0, 0.5, 1.0, 0.0);
        let f_neg = morison_force(1.0, 0.0, 1000.0, 0.5, -1.0, 0.0);
        assert!(f_pos > 0.0);
        assert!(f_neg < 0.0);
        assert!((f_pos + f_neg).abs() < 1e-12);
    }

    #[test]
    fn test_vortex_induced_vibration_frequency_formula() {
        let f = vortex_induced_vibration_frequency(1.0, 0.5, 0.2);
        assert!((f - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_vortex_induced_vibration_frequency_zero_diameter() {
        let f = vortex_induced_vibration_frequency(1.0, 0.0, 0.2);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_vortex_induced_vibration_frequency_proportional_u() {
        let f1 = vortex_induced_vibration_frequency(1.0, 0.5, 0.2);
        let f2 = vortex_induced_vibration_frequency(2.0, 0.5, 0.2);
        assert!((f2 - 2.0 * f1).abs() < 1e-12);
    }

    #[test]
    fn test_lock_in_condition_true() {
        assert!(lock_in_condition(5.0, 5.05, 0.1));
    }

    #[test]
    fn test_lock_in_condition_false() {
        assert!(!lock_in_condition(5.0, 6.0, 0.1));
    }

    #[test]
    fn test_lock_in_condition_exact() {
        assert!(lock_in_condition(5.0, 5.0, 0.0));
    }

    #[test]
    fn test_coupling_factor_formula() {
        let cf = fluid_structure_coupling_factor(1000.0, 7800.0, 1.0, 1.0);
        assert!((cf - 1000.0 / 7800.0).abs() < 1e-10);
    }

    #[test]
    fn test_coupling_factor_zero_solid_volume() {
        let cf = fluid_structure_coupling_factor(1000.0, 7800.0, 1.0, 0.0);
        assert_eq!(cf, 0.0);
    }

    #[test]
    fn test_ale_velocity_zero_mesh() {
        let c = arbitrary_lagrangian_eulerian_velocity([0.0; 3], [1.0, 2.0, 3.0]);
        assert_eq!(c, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_ale_velocity_formula() {
        let c = arbitrary_lagrangian_eulerian_velocity([1.0, 0.5, 0.25], [2.0, 1.5, 1.25]);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ale_velocity_mesh_faster_than_fluid() {
        let c = arbitrary_lagrangian_eulerian_velocity([3.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(c[0] < 0.0);
    }

    #[test]
    fn test_interface_condition_standard() {
        let ic = InterfaceCondition::standard();
        assert!(ic.no_slip);
        assert!(ic.pressure_continuity);
        assert!(ic.traction_balance);
        assert!(ic.all_satisfied());
    }

    #[test]
    fn test_interface_condition_partial() {
        let ic = InterfaceCondition::new(true, false, true);
        assert!(!ic.all_satisfied());
    }

    #[test]
    fn test_ibm_force_density_zero_diff() {
        let f = ibm_force_density([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1000.0);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_ibm_force_density_formula() {
        let f = ibm_force_density([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 500.0);
        assert!((f[0] - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_ibm_force_density_scales_with_k() {
        let f1 = ibm_force_density([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        let f2 = ibm_force_density([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 200.0);
        assert!((f2[0] - 2.0 * f1[0]).abs() < 1e-10);
    }

    #[test]
    fn test_sloshing_frequency_n1() {
        let g = 9.81;
        let h = 1.0;
        let l = 2.0;
        let omega = sloshing_frequency(g, h, l, 1);
        let k = PI / l;
        let expected = (k * g * (k * h).tanh()).sqrt();
        assert!((omega - expected).abs() < 1e-10);
    }

    #[test]
    fn test_sloshing_frequency_zero_length() {
        let omega = sloshing_frequency(9.81, 1.0, 0.0, 1);
        assert_eq!(omega, 0.0);
    }

    #[test]
    fn test_sloshing_frequency_zero_mode() {
        let omega = sloshing_frequency(9.81, 1.0, 2.0, 0);
        assert_eq!(omega, 0.0);
    }

    #[test]
    fn test_sloshing_frequency_higher_mode_larger() {
        let omega1 = sloshing_frequency(9.81, 1.0, 2.0, 1);
        let omega2 = sloshing_frequency(9.81, 1.0, 2.0, 2);
        assert!(omega2 > omega1);
    }

    #[test]
    fn test_sloshing_frequency_positive() {
        let omega = sloshing_frequency(9.81, 0.5, 1.0, 1);
        assert!(omega > 0.0);
    }

    #[test]
    fn test_morison_force_zero_all() {
        let f = morison_force(1.0, 1.0, 1000.0, 0.5, 0.0, 0.0);
        assert_eq!(f, 0.0);
    }
}
