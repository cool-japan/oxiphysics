// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM chassis co-simulation for flexible body vehicle dynamics.
//!
//! Couples a rigid-body raycast vehicle to a reduced-order FEM chassis model
//! so that chassis flexure feeds back into wheel attachment geometry and
//! ride-quality predictions.
//!
//! # Architecture
//!
//! ```text
//!  ┌─────────────────────────────┐
//!  │     RaycastVehicle (rigid)  │  ←— suspension forces, inertial loads
//!  └──────────────┬──────────────┘
//!                 │  ChassisCosimBridge::couple()
//!  ┌──────────────▼──────────────┐
//!  │  ModalChassisModel (FEM)    │  Φᵀ K Φ η = Φᵀ f_ext
//!  │  Craig-Bampton reduction    │
//!  └──────────────┬──────────────┘
//!                 │  deformation at attachment nodes
//!  ┌──────────────▼──────────────┐
//!  │  Wheel attachment offsets   │  δx_wheel = Φ_attach · η
//!  └─────────────────────────────┘
//! ```
//!
//! ## Reduced-order model
//!
//! The chassis deformation is expressed as
//!
//! **u** = Φ **η**
//!
//! where Φ ∈ ℝ^{3N×K} contains the first K normal-mode shapes and **η** ∈ ℝ^K
//! are the modal coordinates.  The modal equations of motion are
//!
//! **M̃** η̈ + **C̃** η̇ + **K̃** η = Φᵀ **f**
//!
//! with diagonal matrices **M̃** = I (mass-normalised), **K̃** = diag(ω_i²),
//! and proportional damping **C̃** = 2 ζ diag(ω_i).
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_vehicle::fem_chassis_cosim::{
//!     ModalChassisModel, ChassisCosimBridge, ChassisAttachmentPoint,
//! };
//!
//! // Build a 3-mode, 12-node chassis
//! let mut model = ModalChassisModel::new(12, 3);
//! // (In practice: populate mode shapes from an FE solver export.)
//! model.set_natural_frequency(0, 45.0);   // 45 Hz first bending mode
//! model.set_natural_frequency(1, 80.0);
//! model.set_natural_frequency(2, 120.0);
//!
//! // Attach front-left wheel at node 2
//! let mut bridge = ChassisCosimBridge::new(model, 0.02);
//! bridge.add_attachment(ChassisAttachmentPoint { node_index: 2, label: "FL".into() });
//! bridge.add_attachment(ChassisAttachmentPoint { node_index: 3, label: "FR".into() });
//! bridge.add_attachment(ChassisAttachmentPoint { node_index: 8, label: "RL".into() });
//! bridge.add_attachment(ChassisAttachmentPoint { node_index: 9, label: "RR".into() });
//!
//! // Step (normally called each vehicle simulation frame)
//! let suspension_forces = vec![[0.0f64; 3]; 4];
//! bridge.step(&suspension_forces, 1.0 / 60.0);
//!
//! // Query deformation at FL attachment
//! let deform = bridge.attachment_deformation("FL");
//! println!("FL chassis flex: {:?}", deform);
//! ```

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── types ─────────────────────────────────────────────────────────────────────

/// Identification label + FEM node index for a wheel or subframe attachment.
#[derive(Debug, Clone)]
pub struct ChassisAttachmentPoint {
    /// Index into the FEM mesh node list.
    pub node_index: usize,
    /// Human-readable label, e.g. "FL", "FR", "RL", "RR", "subframe_front".
    pub label: String,
}

/// Reduced-order modal chassis model.
///
/// Stores the first K normal-mode shapes and integrates the modal equations
/// of motion under externally applied nodal forces.
///
/// Mode shapes are stored column-major in a flat Vec with stride `3 * n_nodes`:
/// `phi[mode * 3 * n_nodes + node * 3 + dof]`.
#[derive(Debug, Clone)]
pub struct ModalChassisModel {
    /// Number of FEM mesh nodes (3 DOFs each: x, y, z).
    pub n_nodes: usize,
    /// Number of retained mode shapes.
    pub n_modes: usize,
    /// Mode shapes Φ ∈ ℝ^{3n_nodes × n_modes}, column-major flat storage.
    pub phi: Vec<f64>,
    /// Natural frequencies ω_i (rad/s) for each mode.
    pub omega: Vec<f64>,
    /// Modal damping ratios ζ_i (dimensionless, default 0.02).
    pub zeta: Vec<f64>,
    /// Current modal coordinates η (size: n_modes).
    pub eta: Vec<f64>,
    /// Current modal velocities η̇ (size: n_modes).
    pub eta_dot: Vec<f64>,
    /// Current modal accelerations η̈ (size: n_modes) — carried between steps for Newmark-β.
    pub eta_ddot: Vec<f64>,
    /// Nodal mass vector (diagonal lumped mass, size: 3 * n_nodes).
    pub lumped_mass: Vec<f64>,
}

impl ModalChassisModel {
    /// Create a new modal chassis model with `n_nodes` nodes and `n_modes` retained modes.
    ///
    /// Mode shapes are initialised to zero; natural frequencies default to 10 Hz.
    /// Lumped mass defaults to 1 kg per DOF; set with [`Self::set_node_mass`] or [`Self::set_uniform_mass`].
    pub fn new(n_nodes: usize, n_modes: usize) -> Self {
        let n_dof = 3 * n_nodes;
        let default_omega = 2.0 * std::f64::consts::PI * 10.0;
        Self {
            n_nodes,
            n_modes,
            phi: vec![0.0; n_modes * n_dof],
            omega: vec![default_omega; n_modes],
            zeta: vec![0.02; n_modes],
            eta: vec![0.0; n_modes],
            eta_dot: vec![0.0; n_modes],
            eta_ddot: vec![0.0; n_modes],
            lumped_mass: vec![1.0; n_dof],
        }
    }

    /// Set the natural frequency for mode `m` in Hz.
    pub fn set_natural_frequency(&mut self, m: usize, freq_hz: f64) {
        assert!(m < self.n_modes, "mode index out of range");
        self.omega[m] = 2.0 * std::f64::consts::PI * freq_hz;
    }

    /// Set the modal damping ratio ζ for mode `m`.
    pub fn set_damping_ratio(&mut self, m: usize, zeta: f64) {
        assert!(m < self.n_modes);
        self.zeta[m] = zeta;
    }

    /// Set a single entry of the mode shape matrix.
    ///
    /// `m` — mode index, `node` — node index, `dof` — 0 (x), 1 (y), 2 (z).
    pub fn set_phi(&mut self, m: usize, node: usize, dof: usize, value: f64) {
        assert!(m < self.n_modes && node < self.n_nodes && dof < 3);
        let idx = m * 3 * self.n_nodes + node * 3 + dof;
        self.phi[idx] = value;
    }

    /// Set uniform nodal mass (total chassis mass evenly distributed).
    pub fn set_uniform_mass(&mut self, total_kg: f64) {
        let per_dof = total_kg / self.lumped_mass.len() as f64;
        for m in &mut self.lumped_mass {
            *m = per_dof;
        }
    }

    /// Set the lumped mass for a specific node (same value for all 3 DOFs).
    pub fn set_node_mass(&mut self, node: usize, mass_kg: f64) {
        assert!(node < self.n_nodes);
        for dof in 0..3 {
            self.lumped_mass[node * 3 + dof] = mass_kg;
        }
    }

    /// Reset modal coordinates, velocities, and accelerations to zero.
    pub fn reset(&mut self) {
        self.eta.fill(0.0);
        self.eta_dot.fill(0.0);
        self.eta_ddot.fill(0.0);
    }

    // ── modal projection ────────────────────────────────────────────────────

    /// Project nodal force vector **f** (size 3*n_nodes) onto modal coordinates.
    ///
    /// Returns Φᵀ M⁻¹ f (mass-normalised modal force) for each mode.
    /// For mass-normalised mode shapes, this reduces to Φᵀ f.
    pub fn modal_force(&self, nodal_force: &[f64]) -> Vec<f64> {
        assert_eq!(nodal_force.len(), 3 * self.n_nodes);
        let n_dof = 3 * self.n_nodes;
        let mut qf = vec![0.0_f64; self.n_modes];
        for (m, qf_m) in qf.iter_mut().enumerate() {
            let phi_m = &self.phi[m * n_dof..(m + 1) * n_dof];
            *qf_m = phi_m.iter().zip(nodal_force).map(|(a, b)| a * b).sum();
        }
        qf
    }

    /// Integrate modal equations of motion for one time step `dt` (seconds)
    /// using Newmark-β average acceleration (β=0.25, γ=0.5 — unconditionally stable).
    ///
    /// The formulation is the standard three-step Newmark algorithm that tracks
    /// acceleration explicitly, ensuring correct convergence to static equilibrium
    /// under constant loading.
    ///
    /// For mass-normalised modes the equations of motion are:
    ///   η̈ + 2ζω η̇ + ω² η = Q   (M̃ = I, K̃ = ω², C̃ = 2ζω)
    ///
    /// `modal_force_vec` — Φᵀ f from [`Self::modal_force`].
    pub fn integrate(&mut self, modal_force_vec: &[f64], dt: f64) {
        assert_eq!(modal_force_vec.len(), self.n_modes);

        // Newmark-β average acceleration: β = 1/4, γ = 1/2
        // Predictor (explicit):
        //   η_pred    = η_n + dt·η̇_n + dt²·(0.5 − β)·η̈_n
        //   η̇_pred   = η̇_n + dt·(1 − γ)·η̈_n
        // Effective stiffness:
        //   k_eff = 1 + γ·c·dt + β·k·dt²   (mass M̃ = 1)
        // Solve for η̈_{n+1}:
        //   k_eff·η̈_{n+1} = Q − c·η̇_pred − k·η_pred
        // Corrector:
        //   η_{n+1}  = η_pred + β·dt²·η̈_{n+1}
        //   η̇_{n+1} = η̇_pred + γ·dt·η̈_{n+1}
        const BETA: f64 = 0.25;
        const GAMMA: f64 = 0.5;

        let iter = self
            .omega
            .iter()
            .zip(self.zeta.iter())
            .zip(modal_force_vec.iter())
            .zip(self.eta.iter_mut())
            .zip(self.eta_dot.iter_mut())
            .zip(self.eta_ddot.iter_mut());
        for (((((omega, zeta), qf), eta), eta_dot), eta_ddot) in iter {
            let c = 2.0 * zeta * omega;
            let k = omega * omega;

            let eta_n = *eta;
            let eta_dot_n = *eta_dot;
            let eta_ddot_n = *eta_ddot;

            // Predictor
            let eta_pred = eta_n + dt * eta_dot_n + dt * dt * (0.5 - BETA) * eta_ddot_n;
            let eta_dot_pred = eta_dot_n + dt * (1.0 - GAMMA) * eta_ddot_n;

            // Effective stiffness (mass = 1 for mass-normalised modes)
            let k_eff = 1.0 + GAMMA * c * dt + BETA * k * dt * dt;

            // Residual → acceleration at n+1
            let rhs = qf - c * eta_dot_pred - k * eta_pred;
            let eta_ddot_n1 = rhs / k_eff;

            // Corrector
            *eta = eta_pred + BETA * dt * dt * eta_ddot_n1;
            *eta_dot = eta_dot_pred + GAMMA * dt * eta_ddot_n1;
            *eta_ddot = eta_ddot_n1;
        }
    }

    /// Reconstruct the full nodal displacement vector **u** = Φ **η**.
    pub fn nodal_displacement(&self) -> Vec<f64> {
        let n_dof = 3 * self.n_nodes;
        let mut u = vec![0.0_f64; n_dof];
        for m in 0..self.n_modes {
            let phi_m = &self.phi[m * n_dof..(m + 1) * n_dof];
            let eta_m = self.eta[m];
            for (ui, &phi_mi) in u.iter_mut().zip(phi_m) {
                *ui += phi_mi * eta_m;
            }
        }
        u
    }

    /// Deformation at a specific node [x, y, z] (metres).
    pub fn node_deformation(&self, node: usize) -> [f64; 3] {
        assert!(node < self.n_nodes);
        let n_dof = 3 * self.n_nodes;
        let mut d = [0.0_f64; 3];
        for (m, eta_m) in self.eta.iter().enumerate() {
            let base = m * n_dof + node * 3;
            for (dof, d_dof) in d.iter_mut().enumerate() {
                *d_dof += self.phi[base + dof] * eta_m;
            }
        }
        d
    }

    /// Maximum absolute nodal displacement across all nodes and DOFs.
    pub fn max_displacement(&self) -> f64 {
        let u = self.nodal_displacement();
        u.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
    }

    /// Root-mean-square of modal coordinates (dimensionless deformation indicator).
    pub fn eta_rms(&self) -> f64 {
        let sum2: f64 = self.eta.iter().map(|v| v * v).sum();
        (sum2 / self.n_modes as f64).sqrt()
    }
}

// ── ChassisCosimBridge ────────────────────────────────────────────────────────

/// Log entry for one co-simulation time step.
#[derive(Debug, Clone)]
pub struct CosimLogEntry {
    /// Simulation time in seconds.
    pub time: f64,
    /// RMS modal deformation at this step.
    pub eta_rms: f64,
    /// Maximum nodal displacement at this step (metres).
    pub max_disp_m: f64,
}

/// Bridge between a rigid-body vehicle model and the [`ModalChassisModel`].
///
/// Each call to [`Self::step`] applies the current suspension attachment forces to
/// the FEM model, integrates the modal equations, and makes deformed attachment
/// positions available for the next vehicle dynamics step.
pub struct ChassisCosimBridge {
    /// The reduced-order FEM chassis model.
    pub model: ModalChassisModel,
    /// Registered attachment points (wheel mounts, subframe brackets).
    pub attachments: Vec<ChassisAttachmentPoint>,
    /// Fixed sub-step size used inside [`Self::step`] for FEM integration (seconds).
    pub sub_dt: f64,
    /// Total accumulated simulation time.
    pub time: f64,
    /// Simulation log.
    pub log: Vec<CosimLogEntry>,
    /// Cached nodal displacement from the last step (size: 3*n_nodes).
    last_u: Vec<f64>,
    /// Scale factor for inertial load contribution (d'Alembert, default 1.0).
    pub inertial_load_scale: f64,
    /// Number of sub-steps per vehicle frame (default 1).
    pub sub_steps: usize,
}

impl ChassisCosimBridge {
    /// Create a new co-simulation bridge.
    ///
    /// `sub_dt` — the FEM integration time step.  Should satisfy the Nyquist
    /// criterion for the highest retained mode: `sub_dt ≤ π / ω_max`.
    pub fn new(model: ModalChassisModel, sub_dt: f64) -> Self {
        let n_dof = 3 * model.n_nodes;
        Self {
            model,
            attachments: Vec::new(),
            sub_dt,
            time: 0.0,
            log: Vec::new(),
            last_u: vec![0.0; n_dof],
            inertial_load_scale: 1.0,
            sub_steps: 1,
        }
    }

    /// Add a wheel or bracket attachment point.
    pub fn add_attachment(&mut self, pt: ChassisAttachmentPoint) {
        self.attachments.push(pt);
    }

    /// Find an attachment by label, returning its index in `self.attachments`.
    pub fn find_attachment(&self, label: &str) -> Option<usize> {
        self.attachments.iter().position(|a| a.label == label)
    }

    /// Deformation [x, y, z] in metres at the attachment labelled `label`.
    ///
    /// Returns `[0.0; 3]` if the label is not registered.
    pub fn attachment_deformation(&self, label: &str) -> [f64; 3] {
        match self.find_attachment(label) {
            Some(idx) => {
                let node = self.attachments[idx].node_index;
                self.model.node_deformation(node)
            }
            None => [0.0; 3],
        }
    }

    /// Deformations for all registered attachments in registration order.
    pub fn all_attachment_deformations(&self) -> Vec<([f64; 3], &str)> {
        self.attachments
            .iter()
            .map(|a| (self.model.node_deformation(a.node_index), a.label.as_str()))
            .collect()
    }

    /// Advance the co-simulation by one vehicle frame.
    ///
    /// # Parameters
    ///
    /// - `suspension_forces`: forces applied at each attachment node
    ///   (same order as [`Self::add_attachment`] calls), each `[fx, fy, fz]` in N.
    /// - `dt`: vehicle frame duration in seconds.
    ///
    /// Internally runs `self.sub_steps` FEM sub-steps of size `self.sub_dt`.
    /// If `dt / sub_steps` differs from `sub_dt`, the last sub-step is shortened
    /// to fill the remainder (accurate Newmark integration).
    pub fn step(&mut self, suspension_forces: &[[f64; 3]], dt: f64) {
        let n_dof = 3 * self.model.n_nodes;

        // Assemble nodal force vector
        let mut f = vec![0.0_f64; n_dof];
        for (i, attach) in self.attachments.iter().enumerate() {
            let node = attach.node_index;
            if node < self.model.n_nodes {
                let force = if i < suspension_forces.len() {
                    suspension_forces[i]
                } else {
                    [0.0; 3]
                };
                for dof in 0..3 {
                    f[node * 3 + dof] += force[dof];
                }
            }
        }

        // Project to modal space
        let qf = self.model.modal_force(&f);

        // Sub-step integration
        let n_sub = self.sub_steps.max(1);
        let actual_sub_dt = dt / n_sub as f64;
        for _ in 0..n_sub {
            self.model.integrate(&qf, actual_sub_dt);
        }

        // Cache nodal displacement
        self.last_u = self.model.nodal_displacement();
        self.time += dt;

        // Log
        self.log.push(CosimLogEntry {
            time: self.time,
            eta_rms: self.model.eta_rms(),
            max_disp_m: self.model.max_displacement(),
        });
    }

    /// Apply a rigid-body inertial load (d'Alembert) distributed across all
    /// chassis nodes.
    ///
    /// `accel_world [ax, ay, az]` — vehicle CoG acceleration in world frame (m/s²).
    /// The inertial force on each node is `f_i = −m_i * a` scaled by
    /// `self.inertial_load_scale`.
    ///
    /// This method assembles the nodal force vector and calls [`Self::step`] with
    /// zero suspension attachment forces.
    pub fn step_inertial(&mut self, accel_world: [f64; 3], dt: f64) {
        let n_dof = 3 * self.model.n_nodes;
        let mut f = vec![0.0_f64; n_dof];
        for node in 0..self.model.n_nodes {
            let m = self.model.lumped_mass[node * 3]; // same for all 3 DOFs
            for dof in 0..3 {
                f[node * 3 + dof] += -self.inertial_load_scale * m * accel_world[dof];
            }
        }
        let qf = self.model.modal_force(&f);
        let n_sub = self.sub_steps.max(1);
        let sub_dt = dt / n_sub as f64;
        for _ in 0..n_sub {
            self.model.integrate(&qf, sub_dt);
        }
        self.last_u = self.model.nodal_displacement();
        self.time += dt;
        self.log.push(CosimLogEntry {
            time: self.time,
            eta_rms: self.model.eta_rms(),
            max_disp_m: self.model.max_displacement(),
        });
    }

    /// Combined step: suspension attachment forces + rigid-body inertial loading.
    ///
    /// The most accurate coupling: applies both wheel-load-induced bending and
    /// rigid-body d'Alembert forces in the same time step.
    pub fn step_coupled(&mut self, suspension_forces: &[[f64; 3]], accel_world: [f64; 3], dt: f64) {
        let n_dof = 3 * self.model.n_nodes;
        let mut f = vec![0.0_f64; n_dof];

        // Suspension loads at attachment nodes
        for (i, attach) in self.attachments.iter().enumerate() {
            let node = attach.node_index;
            if node < self.model.n_nodes {
                let force = if i < suspension_forces.len() {
                    suspension_forces[i]
                } else {
                    [0.0; 3]
                };
                for dof in 0..3 {
                    f[node * 3 + dof] += force[dof];
                }
            }
        }

        // Inertial (d'Alembert) loads at all nodes
        for node in 0..self.model.n_nodes {
            let m = self.model.lumped_mass[node * 3];
            for dof in 0..3 {
                f[node * 3 + dof] -= self.inertial_load_scale * m * accel_world[dof];
            }
        }

        let qf = self.model.modal_force(&f);
        let n_sub = self.sub_steps.max(1);
        let sub_dt = dt / n_sub as f64;
        for _ in 0..n_sub {
            self.model.integrate(&qf, sub_dt);
        }
        self.last_u = self.model.nodal_displacement();
        self.time += dt;
        self.log.push(CosimLogEntry {
            time: self.time,
            eta_rms: self.model.eta_rms(),
            max_disp_m: self.model.max_displacement(),
        });
    }

    /// Reset the chassis model and simulation clock.
    pub fn reset(&mut self) {
        self.model.reset();
        self.time = 0.0;
        self.log.clear();
        self.last_u.fill(0.0);
    }
}

// ── chassis mode-shape generation helpers ─────────────────────────────────────

/// Generate synthetic cantilever-beam mode shapes for a chassis represented as
/// a straight N-node beam along the X axis.
///
/// Mode shapes are mass-normalised sinusoidal bending shapes in the Y (vertical)
/// direction and are stored directly into `model.phi`.  Natural frequencies
/// follow the Euler-Bernoulli beam relation: ω_m ∝ (m+1)² × π² √(EI/ρAL⁴).
///
/// This is useful for unit-testing the co-simulation without importing an FE
/// solver.  For production, supply mode shapes from a real FE analysis.
///
/// # Parameters
/// - `model`        — target modal chassis model (modified in-place)
/// - `length_m`     — total beam length (m)
/// - `ei_over_rho_a` — bending stiffness / (linear mass density) EI/(ρA) in m⁴/s²
pub fn populate_beam_modes(model: &mut ModalChassisModel, length_m: f64, ei_over_rho_a: f64) {
    let n = model.n_nodes;
    let n_dof = 3 * n;
    let l = length_m;

    for m in 0..model.n_modes {
        // Euler-Bernoulli: ω_m = (m+1)² π² / L² √(EI/ρA)
        let beta_l = (m as f64 + 1.0) * std::f64::consts::PI;
        let omega_m = beta_l * beta_l / (l * l) * ei_over_rho_a.sqrt();
        model.omega[m] = omega_m;
        let freq_hz = omega_m / (2.0 * std::f64::consts::PI);
        // clamp unrealistic values
        model.omega[m] = (freq_hz * 2.0 * std::f64::consts::PI).max(1.0);

        // Sinusoidal Y-direction shape: φ_y(x) = sin((m+1)πx/L)
        // Normalise so that Σ φᵢ² / N ≈ 1 (crude mass-normalisation)
        let norm = (n as f64).sqrt();
        for node in 0..n {
            let x = node as f64 * l / (n as f64 - 1.0).max(1.0);
            let shape_y = ((m as f64 + 1.0) * std::f64::consts::PI * x / l).sin() / norm;
            let idx = m * n_dof + node * 3;
            model.phi[idx] = 0.0; // x DOF: no axial
            model.phi[idx + 1] = shape_y; // y DOF: vertical bending
            model.phi[idx + 2] = 0.0; // z DOF: no lateral
        }
    }
}

/// Query the bending stress at node `node` due to modal deformation.
///
/// Uses the Euler-Bernoulli approximation:
/// σ = E · c · κ,  where κ = d²u_y/dx²
///
/// The second derivative is estimated by finite differences across the
/// three adjacent nodes.  Returns 0 at boundary nodes.
///
/// # Parameters
/// - `model`     — modal chassis model
/// - `node`      — node index (1..n_nodes-2 for valid FD stencil)
/// - `e_modulus` — Young's modulus in Pa
/// - `c_dist`    — distance from neutral axis to extreme fibre in m
/// - `dx`        — nodal spacing in m
pub fn bending_stress(
    model: &ModalChassisModel,
    node: usize,
    e_modulus: f64,
    c_dist: f64,
    dx: f64,
) -> f64 {
    if node == 0 || node + 1 >= model.n_nodes {
        return 0.0;
    }
    let u_prev = model.node_deformation(node - 1)[1];
    let u_curr = model.node_deformation(node)[1];
    let u_next = model.node_deformation(node + 1)[1];
    let kappa = (u_next - 2.0 * u_curr + u_prev) / (dx * dx);
    e_modulus * c_dist * kappa.abs()
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modal_model_construction() {
        let mut m = ModalChassisModel::new(10, 3);
        assert_eq!(m.n_nodes, 10);
        assert_eq!(m.n_modes, 3);
        m.set_natural_frequency(0, 45.0);
        let expected = 2.0 * std::f64::consts::PI * 45.0;
        assert!((m.omega[0] - expected).abs() < 1e-6);
    }

    #[test]
    fn test_modal_force_projection() {
        let mut m = ModalChassisModel::new(4, 2);
        // Set identity-like mode shapes: mode 0 → Y of node 0; mode 1 → Y of node 2
        m.set_phi(0, 0, 1, 1.0);
        m.set_phi(1, 2, 1, 1.0);

        let mut f = vec![0.0; 12];
        f[1] = 100.0; // node 0, Y direction
        let qf = m.modal_force(&f);
        assert!((qf[0] - 100.0).abs() < 1e-10);
        assert!(qf[1].abs() < 1e-10);
    }

    #[test]
    fn test_integration_static_load() {
        // Static load: apply constant force and check convergence to static equilibrium.
        // For a mass-normalised modal oscillator:  η̈ + 2ζω η̇ + ω² η = Q
        // Static solution:  η_static = Q / ω²
        //
        // With ζ = 0.02 (lightly damped) the settling time to within 5% of the
        // static value is approximately  t_settle = ln(0.05) / (−ζ·ω) ≈ 2.4 s,
        // i.e. roughly 24 periods at 10 Hz.  We integrate for 25 periods to
        // guarantee convergence within the 5 % tolerance.
        let mut m = ModalChassisModel::new(2, 1);
        m.set_natural_frequency(0, 10.0); // 10 Hz → ω = 2π·10 rad/s
        // zeta defaults to 0.02 — no override needed
        m.set_phi(0, 0, 1, 1.0);

        let omega = m.omega[0];
        let f_modal = vec![1.0_f64]; // unit modal force Q = 1
        let dt = 1e-4;

        // Integrate for 25 natural periods (sufficient for ζ=0.02 to settle to <5 %)
        let n_periods = 25;
        let t_end = n_periods as f64 / 10.0; // 2.5 s
        let n_steps = (t_end / dt) as usize;
        for _ in 0..n_steps {
            m.integrate(&f_modal, dt);
        }

        // Static equilibrium: η = Q / ω²  (K̃ η = Q, K̃ = diag(ω²))
        let eta_static = 1.0 / (omega * omega);

        assert!(
            (m.eta[0] - eta_static).abs() < eta_static * 0.05,
            "eta={:.6} static={:.6}",
            m.eta[0],
            eta_static
        );
    }

    #[test]
    fn test_bridge_step() {
        let mut model = ModalChassisModel::new(10, 2);
        model.set_natural_frequency(0, 30.0);
        model.set_natural_frequency(1, 60.0);
        populate_beam_modes(&mut model, 4.0, 1e6);

        let mut bridge = ChassisCosimBridge::new(model, 1.0 / 1000.0);
        bridge.add_attachment(ChassisAttachmentPoint {
            node_index: 2,
            label: "FL".into(),
        });
        bridge.add_attachment(ChassisAttachmentPoint {
            node_index: 7,
            label: "RR".into(),
        });

        // Apply 5 kN upward at FL, 4 kN upward at RR
        let forces = vec![[0.0, 5000.0, 0.0], [0.0, 4000.0, 0.0]];
        bridge.sub_steps = 10;
        bridge.step(&forces, 1.0 / 60.0);

        let fl_deform = bridge.attachment_deformation("FL");
        println!("FL deformation: {:?}", fl_deform);
        assert!(bridge.log.len() == 1);
        assert!(bridge.time > 0.0);
    }

    #[test]
    fn test_beam_mode_populate() {
        let mut m = ModalChassisModel::new(20, 4);
        populate_beam_modes(&mut m, 4.5, 5e5);
        // All natural frequencies should be positive
        for &omega in &m.omega {
            assert!(omega > 0.0, "omega should be positive, got {}", omega);
        }
        // Mode shapes should be non-trivial
        let phi_sum: f64 = m.phi.iter().map(|v| v.abs()).sum();
        assert!(phi_sum > 0.0);
    }
}
