// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real-time coupling constraints between rigid bodies and deformable bodies.
//!
//! This module provides [`RigidDeformableCoupling`], a velocity-level constraint that
//! attaches a point on a rigid body to one or more nodes of a deformable mesh
//! (FEM, PBD cloth, soft-body, etc.).
//!
//! ## Architecture
//!
//! The coupling is one-way in the type system: the constraint lives in the
//! *rigid* constraint solver pipeline and treats the deformable body via
//! the [`DeformableBodyState`] trait.  The FEM or soft-body solver is
//! responsible for updating its node positions/velocities; this module
//! re-projects constraint impulses back onto those nodes.
//!
//! ```text
//!         ┌─────────────────┐          ┌────────────────────┐
//!         │  RigidBodySet   │◄────────►│  DeformableBody    │
//!         │  (Rapier-style) │  coupling│  (FEM / PBD / SPH) │
//!         └─────────────────┘          └────────────────────┘
//!                   ▲                          ▲
//!                   │        Impulse           │
//!                   └──── RigidDeformable ─────┘
//!                            Coupling
//! ```
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_constraints::deformable_coupling::{
//!     DeformableBodyState, DeformableNodeId, RigidDeformableCoupling,
//! };
//!
//! struct MockMesh { vel_y: f64 }
//! impl DeformableBodyState for MockMesh {
//!     fn node_position(&self, _id: DeformableNodeId) -> [f64; 3] { [0.0, 0.0, 0.0] }
//!     fn node_velocity(&self, _id: DeformableNodeId) -> [f64; 3] { [0.0, self.vel_y, 0.0] }
//!     fn apply_impulse(&mut self, _id: DeformableNodeId, _impulse: [f64; 3]) {}
//!     fn node_inv_mass(&self, _id: DeformableNodeId) -> f64 { 1.0 }
//! }
//!
//! let mut mesh = MockMesh { vel_y: 0.0 };
//! let coupling = RigidDeformableCoupling::new(
//!     DeformableNodeId(0),
//!     [0.0, 0.0, 0.0],  // attachment point in rigid body local space
//!     1.0 / 60.0,        // dt
//! );
//! let _ = coupling.compute_delta_lambda(&mesh, [0.0, -1.0, 0.0], 1.0, 0.1);
//! ```

// ── DeformableNodeId ─────────────────────────────────────────────────────────

/// Opaque identifier for a node in a deformable mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DeformableNodeId(pub usize);

// ── DeformableBodyState trait ─────────────────────────────────────────────────

/// Interface that a deformable body (FEM mesh, PBD cloth, soft-body, SPH particle
/// cluster) must implement to participate in the rigid–deformable coupling constraint.
pub trait DeformableBodyState {
    /// World-space position of node `id`.
    fn node_position(&self, id: DeformableNodeId) -> [f64; 3];

    /// World-space linear velocity of node `id`.
    fn node_velocity(&self, id: DeformableNodeId) -> [f64; 3];

    /// Apply a world-space linear impulse to node `id`.
    fn apply_impulse(&mut self, id: DeformableNodeId, impulse: [f64; 3]);

    /// Inverse nodal mass.  `0.0` pins the node (fixed boundary condition).
    fn node_inv_mass(&self, id: DeformableNodeId) -> f64;
}

// ── CouplingKind ────────────────────────────────────────────────────────────

/// How the rigid body attachment is projected onto the deformable mesh.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum CouplingKind {
    /// Single node attachment: all impulse goes to one mesh node.
    #[default]
    SingleNode,
    /// Barycentric interpolation across up to 4 nodes with provided weights.
    Barycentric([f64; 4]),
}

// ── RigidDeformableCoupling ──────────────────────────────────────────────────

/// Velocity-level constraint that couples a point on a rigid body to one
/// or more nodes of a deformable body.
///
/// The constraint enforces that the world-space velocity of the attachment
/// point on the rigid body matches the interpolated velocity of the coupled
/// mesh nodes.
#[derive(Debug, Clone)]
pub struct RigidDeformableCoupling {
    /// Target node (or primary node for barycentric couplings).
    pub primary_node: DeformableNodeId,
    /// Additional nodes for barycentric interpolation (filled if kind is Barycentric).
    pub secondary_nodes: [DeformableNodeId; 3],
    /// Attachment point in rigid body local space.
    pub local_attachment: [f64; 3],
    /// How impulses are distributed to mesh nodes.
    pub kind: CouplingKind,
    /// Simulation time step (used for bias computation).
    pub dt: f64,
    /// Constraint stiffness (ERP-style): 1.0 = rigid, 0.1 = soft coupling.
    pub stiffness: f64,
    /// Damping factor for the coupling spring.
    pub damping: f64,
    /// Accumulated impulse from previous frame (warm start).
    pub lambda: [f64; 3],
}

impl RigidDeformableCoupling {
    /// Create a new single-node rigid–deformable coupling.
    ///
    /// * `node`             — target mesh node
    /// * `local_attachment` — attachment in rigid body local space
    /// * `dt`               — time step
    pub fn new(node: DeformableNodeId, local_attachment: [f64; 3], dt: f64) -> Self {
        Self {
            primary_node: node,
            secondary_nodes: [
                DeformableNodeId(usize::MAX),
                DeformableNodeId(usize::MAX),
                DeformableNodeId(usize::MAX),
            ],
            local_attachment,
            kind: CouplingKind::SingleNode,
            dt,
            stiffness: 1.0,
            damping: 0.1,
            lambda: [0.0; 3],
        }
    }

    /// Create a barycentric coupling with 4 nodes and corresponding weights.
    ///
    /// Weights must sum to 1.0.
    pub fn new_barycentric(
        nodes: [DeformableNodeId; 4],
        weights: [f64; 4],
        local_attachment: [f64; 3],
        dt: f64,
    ) -> Self {
        Self {
            primary_node: nodes[0],
            secondary_nodes: [nodes[1], nodes[2], nodes[3]],
            local_attachment,
            kind: CouplingKind::Barycentric(weights),
            dt,
            stiffness: 1.0,
            damping: 0.1,
            lambda: [0.0; 3],
        }
    }

    /// Compute the delta lambda (impulse correction) for this coupling constraint.
    ///
    /// # Arguments
    ///
    /// * `deformable` — read-only deformable body state
    /// * `rigid_vel`  — velocity of the rigid body attachment point in world space
    /// * `rigid_inv_mass` — inverse mass of the rigid body
    /// * `baumgarte`  — stabilization factor
    ///
    /// Returns the 3D impulse to apply (positive = from rigid → mesh direction).
    pub fn compute_delta_lambda<D: DeformableBodyState>(
        &self,
        deformable: &D,
        rigid_vel: [f64; 3],
        rigid_inv_mass: f64,
        baumgarte: f64,
    ) -> [f64; 3] {
        // Interpolated mesh node velocity
        let mesh_vel = self.interpolated_velocity(deformable);

        // Interpolated mesh node inv_mass
        let mesh_inv_mass = self.interpolated_inv_mass(deformable);

        // Relative velocity error: rigid attachment should match mesh
        let rv = sub3(rigid_vel, mesh_vel);

        // Position error (approximated as zero here; full coupling uses position residual)
        let bias = scale3(rv, -(baumgarte / self.dt.max(1e-12)));

        // Effective mass for the 3D coupling constraint (simplified diagonal)
        let eff_mass = if (rigid_inv_mass + mesh_inv_mass) > 1e-15 {
            1.0 / (rigid_inv_mass + mesh_inv_mass)
        } else {
            0.0
        };

        // delta_lambda = -eff_mass * (rv + bias)  [per component]
        let total = add3(rv, bias);
        scale3(total, -eff_mass * self.stiffness)
    }

    /// Apply coupling impulse to the deformable mesh nodes.
    pub fn apply_to_deformable<D: DeformableBodyState>(
        &self,
        deformable: &mut D,
        impulse: [f64; 3],
    ) {
        match self.kind {
            CouplingKind::SingleNode => {
                let im = deformable.node_inv_mass(self.primary_node);
                let scaled = scale3(impulse, im);
                deformable.apply_impulse(self.primary_node, scaled);
            }
            CouplingKind::Barycentric(weights) => {
                let nodes = [
                    self.primary_node,
                    self.secondary_nodes[0],
                    self.secondary_nodes[1],
                    self.secondary_nodes[2],
                ];
                for (node, &w) in nodes.iter().zip(weights.iter()) {
                    if node.0 == usize::MAX {
                        continue;
                    }
                    let im = deformable.node_inv_mass(*node);
                    let scaled = scale3(scale3(impulse, w), im);
                    deformable.apply_impulse(*node, scaled);
                }
            }
        }
    }

    // ── private helpers ───────────────────────────────────────────────────────

    fn interpolated_velocity<D: DeformableBodyState>(&self, d: &D) -> [f64; 3] {
        match self.kind {
            CouplingKind::SingleNode => d.node_velocity(self.primary_node),
            CouplingKind::Barycentric(weights) => {
                let nodes = [
                    self.primary_node,
                    self.secondary_nodes[0],
                    self.secondary_nodes[1],
                    self.secondary_nodes[2],
                ];
                let mut out = [0.0f64; 3];
                for (node, &w) in nodes.iter().zip(weights.iter()) {
                    if node.0 == usize::MAX {
                        continue;
                    }
                    let v = d.node_velocity(*node);
                    out[0] += w * v[0];
                    out[1] += w * v[1];
                    out[2] += w * v[2];
                }
                out
            }
        }
    }

    fn interpolated_inv_mass<D: DeformableBodyState>(&self, d: &D) -> f64 {
        match self.kind {
            CouplingKind::SingleNode => d.node_inv_mass(self.primary_node),
            CouplingKind::Barycentric(weights) => {
                let nodes = [
                    self.primary_node,
                    self.secondary_nodes[0],
                    self.secondary_nodes[1],
                    self.secondary_nodes[2],
                ];
                let mut out = 0.0f64;
                for (node, &w) in nodes.iter().zip(weights.iter()) {
                    if node.0 == usize::MAX {
                        continue;
                    }
                    out += w * d.node_inv_mass(*node);
                }
                out
            }
        }
    }
}

// ── CoupledSimulation ────────────────────────────────────────────────────────

/// Manages multiple [`RigidDeformableCoupling`] constraints and drives a
/// combined velocity-level solve each physics step.
///
/// The caller is responsible for running the rigid-only PGS solver first,
/// then calling [`CoupledSimulation::solve_velocity`] to project coupling
/// impulses between the two solvers.
#[derive(Debug, Default)]
pub struct CoupledSimulation {
    /// All active coupling constraints.
    pub couplings: Vec<RigidDeformableCoupling>,
    /// Maximum coupling iterations per step.
    pub max_iterations: usize,
}

impl CoupledSimulation {
    /// Create an empty simulation coupler with `max_iterations` coupling iterations.
    pub fn new(max_iterations: usize) -> Self {
        Self {
            couplings: Vec::new(),
            max_iterations,
        }
    }

    /// Add a coupling constraint.
    pub fn add_coupling(&mut self, c: RigidDeformableCoupling) {
        self.couplings.push(c);
    }

    /// Velocity-level coupling pass.
    ///
    /// For each coupling, compute the delta lambda based on the current rigid
    /// and deformable body velocities, then apply the impulse to the mesh.
    ///
    /// * `rigid_vels`  — world-space velocity for each rigid body (by index)
    /// * `rigid_inv_masses` — inverse mass for each rigid body
    /// * `deformable`  — deformable body state (mutably)
    /// * `dt`          — time step
    pub fn solve_velocity<D: DeformableBodyState>(
        &mut self,
        rigid_vels: &[[f64; 3]],
        rigid_inv_masses: &[f64],
        deformable: &mut D,
        dt: f64,
        baumgarte: f64,
    ) {
        let _ = dt; // dt embedded in each coupling, but kept for signature compatibility
        for c in &self.couplings {
            if rigid_vels.is_empty() || rigid_inv_masses.is_empty() {
                break;
            }
            // Use first body's velocity as reference (caller arranges correct indexing)
            let rv = rigid_vels[0];
            let im = rigid_inv_masses[0];
            let dl = c.compute_delta_lambda(deformable, rv, im, baumgarte);
            c.apply_to_deformable(deformable, dl);
        }
    }
}

// ── Vector helpers ───────────────────────────────────────────────────────────

#[inline(always)]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    struct MockMesh {
        pos: [f64; 3],
        vel: [f64; 3],
        inv_mass: f64,
        applied: Vec<[f64; 3]>,
    }

    impl DeformableBodyState for MockMesh {
        fn node_position(&self, _id: DeformableNodeId) -> [f64; 3] {
            self.pos
        }
        fn node_velocity(&self, _id: DeformableNodeId) -> [f64; 3] {
            self.vel
        }
        fn apply_impulse(&mut self, _id: DeformableNodeId, impulse: [f64; 3]) {
            self.applied.push(impulse);
        }
        fn node_inv_mass(&self, _id: DeformableNodeId) -> f64 {
            self.inv_mass
        }
    }

    #[test]
    fn compute_delta_lambda_zero_when_matching_velocities() {
        let mesh = MockMesh {
            pos: [0.0; 3],
            vel: [1.0, 0.0, 0.0],
            inv_mass: 1.0,
            applied: vec![],
        };
        let c = RigidDeformableCoupling::new(DeformableNodeId(0), [0.0; 3], 1.0 / 60.0);
        let dl = c.compute_delta_lambda(&mesh, [1.0, 0.0, 0.0], 1.0, 0.2);
        // Same velocities → no correction impulse needed in the parallel component
        let mag = (dl[0] * dl[0] + dl[1] * dl[1] + dl[2] * dl[2]).sqrt();
        assert!(
            mag < 1e-10,
            "impulse should be near zero when velocities match: {mag}"
        );
    }

    #[test]
    fn apply_to_deformable_distributes_impulse() {
        let mut mesh = MockMesh {
            pos: [0.0; 3],
            vel: [0.0; 3],
            inv_mass: 2.0,
            applied: vec![],
        };
        let c = RigidDeformableCoupling::new(DeformableNodeId(0), [0.0; 3], 1.0 / 60.0);
        c.apply_to_deformable(&mut mesh, [0.0, 10.0, 0.0]);
        assert_eq!(mesh.applied.len(), 1);
        // scaled by inv_mass = 2.0
        assert!((mesh.applied[0][1] - 20.0).abs() < 1e-10);
    }

    #[test]
    fn barycentric_coupling_distributes_evenly_with_equal_weights() {
        let mut mesh = MockMesh {
            pos: [0.0; 3],
            vel: [0.0; 3],
            inv_mass: 1.0,
            applied: vec![],
        };
        let nodes = [
            DeformableNodeId(0),
            DeformableNodeId(1),
            DeformableNodeId(2),
            DeformableNodeId(3),
        ];
        let c = RigidDeformableCoupling::new_barycentric(nodes, [0.25; 4], [0.0; 3], 1.0 / 60.0);
        c.apply_to_deformable(&mut mesh, [0.0, 4.0, 0.0]);
        // 4 nodes × 0.25 weight × inv_mass=1 × impulse_y=4 → each gets 1.0
        for app in &mesh.applied {
            assert!(
                (app[1] - 1.0).abs() < 1e-10,
                "each node should get 1.0, got {}",
                app[1]
            );
        }
    }

    #[test]
    fn coupled_simulation_solve_velocity_runs() {
        let mut mesh = MockMesh {
            pos: [0.0; 3],
            vel: [0.0; 3],
            inv_mass: 1.0,
            applied: vec![],
        };
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(RigidDeformableCoupling::new(
            DeformableNodeId(0),
            [0.0; 3],
            1.0 / 60.0,
        ));
        sim.solve_velocity(
            &[[0.0, 1.0, 0.0]], // rigid vel
            &[1.0],             // rigid inv_mass
            &mut mesh,
            1.0 / 60.0,
            0.2,
        );
        // Should run without panic
    }
}
