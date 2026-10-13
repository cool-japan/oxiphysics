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
//! The constraint lives in the *rigid* constraint solver pipeline and treats the
//! deformable body via the [`DeformableBodyState`] trait.  The FEM or soft-body
//! solver is responsible for updating its node positions/velocities; this module
//! exchanges equal and opposite impulses between the two sides: the mesh nodes
//! receive `w_i * P` (an impulse, which the mesh converts to a velocity change
//! with its own inverse mass) and the rigid body receives `-P`.
//!
//! The coupled quantity is the point-to-point constraint
//! `C = x_rigid - sum_i(w_i * x_i) = 0` with `Cdot = v_rigid - sum_i(w_i * v_i)`.
//! One relaxation solves `K P = stiffness * (Cdot + beta / dt * C)` with the
//! scalar effective inverse mass `K = inv_mass_rigid + sum_i(w_i^2 * inv_mass_i)`.
//!
//! **Rigid-body model:** the rigid body is a *translating point mass*: only its
//! linear velocity and inverse mass take part and the world attachment point is
//! `body_position + local_attachment` (no rotation, no `omega x r` term).
//! Rotating bodies must supply their own attachment kinematics through
//! [`RigidDeformableCoupling::compute_delta_lambda`] /
//! [`RigidDeformableCoupling::solve_iteration`].
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
//!     CoupledSimulation, DeformableBodyState, DeformableNodeId, RigidDeformableCoupling,
//! };
//!
//! /// A one-node "mesh": an impulse changes its velocity by `inv_mass * impulse`.
//! struct OneNode { pos: [f64; 3], vel: [f64; 3], inv_mass: f64 }
//! impl DeformableBodyState for OneNode {
//!     fn node_position(&self, _id: DeformableNodeId) -> [f64; 3] { self.pos }
//!     fn node_velocity(&self, _id: DeformableNodeId) -> [f64; 3] { self.vel }
//!     fn apply_impulse(&mut self, _id: DeformableNodeId, impulse: [f64; 3]) {
//!         for k in 0..3 { self.vel[k] += self.inv_mass * impulse[k]; }
//!     }
//!     fn node_inv_mass(&self, _id: DeformableNodeId) -> f64 { self.inv_mass }
//! }
//!
//! let mut mesh = OneNode { pos: [0.0; 3], vel: [0.0; 3], inv_mass: 1.0 };
//! let mut sim = CoupledSimulation::new(8);
//! sim.add_coupling(RigidDeformableCoupling::new(
//!     DeformableNodeId(0),
//!     [0.0, 0.0, 0.0], // attachment point relative to the rigid body position
//!     1.0 / 60.0,      // dt
//! ));
//!
//! // Rigid body 0 sits on the mesh node and moves along +y (kinematic: inv_mass 0).
//! let rigid_positions = [[0.0, 0.0, 0.0]];
//! let mut rigid_vels = [[0.0, 1.0, 0.0]];
//! let rigid_inv_masses = [0.0];
//! sim.solve_velocity(
//!     &rigid_positions,
//!     &mut rigid_vels,
//!     &rigid_inv_masses,
//!     &mut mesh,
//!     1.0 / 60.0,
//!     0.0,
//! );
//! // The mesh node now follows the rigid body.
//! assert!((mesh.vel[1] - 1.0).abs() < 1e-12);
//! ```

// ── Constants ────────────────────────────────────────────────────────────────

/// Lower bound used whenever a time step is divided by.
const MIN_DT: f64 = 1e-12;

/// Combined inverse mass below which a coupling cannot move either side.
const MIN_INV_MASS: f64 = 1e-15;

/// Coupling iterations used by [`CoupledSimulation::default`].
const DEFAULT_MAX_ITERATIONS: usize = 10;

/// Largest impulse component that still counts as "converged" in
/// [`CoupledSimulation::solve_velocity`].
const CONVERGED_IMPULSE: f64 = 1e-14;

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

    /// Apply a world-space linear impulse (N·s) to node `id`.
    ///
    /// Callers pass the impulse *unscaled*; the implementation converts it to
    /// a velocity change `inv_mass * impulse` with its own nodal inverse mass.
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
/// mesh nodes, with a Baumgarte term that pulls the interpolated mesh position
/// onto the attachment point.
#[derive(Debug, Clone)]
pub struct RigidDeformableCoupling {
    /// Target node (or primary node for barycentric couplings).
    pub primary_node: DeformableNodeId,
    /// Additional nodes for barycentric interpolation (filled if kind is Barycentric).
    pub secondary_nodes: [DeformableNodeId; 3],
    /// Attachment point relative to the rigid body position (the body is
    /// treated as a non-rotating point mass, see the module documentation).
    pub local_attachment: [f64; 3],
    /// How impulses are distributed to mesh nodes.
    pub kind: CouplingKind,
    /// Index of the rigid body this coupling is attached to: the position into
    /// the body arrays handed to [`CoupledSimulation::solve_velocity`].
    pub rigid_body: usize,
    /// Fallback simulation time step for the Baumgarte term, used by
    /// [`RigidDeformableCoupling::compute_delta_lambda`] and whenever the step
    /// passed to the solver is not a positive finite number.
    pub dt: f64,
    /// Constraint stiffness (ERP-style): 1.0 = rigid, 0.1 = soft coupling.
    pub stiffness: f64,
    /// Damping factor for the coupling spring.
    pub damping: f64,
    /// Total impulse this coupling applied to the mesh during the latest
    /// [`CoupledSimulation::solve_velocity`] call (the rigid body received the
    /// opposite impulse).  It is cleared at the start of every solve.
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
            rigid_body: 0,
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
            rigid_body: 0,
            dt,
            stiffness: 1.0,
            damping: 0.1,
            lambda: [0.0; 3],
        }
    }

    /// Select the rigid body this coupling is attached to: an index into the
    /// arrays passed to [`CoupledSimulation::solve_velocity`].
    pub fn with_rigid_body(mut self, rigid_body: usize) -> Self {
        self.rigid_body = rigid_body;
        self
    }

    /// World-space attachment point for a rigid body located at `body_position`.
    ///
    /// The body is treated as a non-rotating point mass, so the attachment is
    /// `body_position + local_attachment`.
    pub fn attachment_world_position(&self, body_position: [f64; 3]) -> [f64; 3] {
        add3(body_position, self.local_attachment)
    }

    /// Compute the coupling impulse `P` of this constraint.
    ///
    /// The constraint is `C = rigid_attach_pos - x_mesh = 0` with `x_mesh` the
    /// weight-interpolated mesh position.  With `rv = rigid_vel - v_mesh` the
    /// impulse that drives `Cdot + (baumgarte / dt) * C` to zero is
    /// `P = stiffness * (rv + (baumgarte / dt) * C) / K` with the scalar
    /// effective inverse mass `K = rigid_inv_mass + sum_i(w_i^2 * inv_mass_i)`.
    ///
    /// # Arguments
    ///
    /// * `deformable`       — read-only deformable body state
    /// * `rigid_attach_pos` — world-space position of the rigid attachment point
    /// * `rigid_vel`        — world-space velocity of the rigid attachment point
    /// * `rigid_inv_mass`   — inverse mass of the rigid body
    /// * `baumgarte`        — position-error feedback factor (`0` matches
    ///   velocities only)
    ///
    /// Returns the impulse to apply **to the mesh** (distributed over its nodes
    /// as `w_i * P`); the rigid body receives `-P`.  The result is zero when
    /// both sides are immovable (`K = 0`).  The Baumgarte term uses
    /// [`RigidDeformableCoupling::dt`].
    pub fn compute_delta_lambda<D: DeformableBodyState>(
        &self,
        deformable: &D,
        rigid_attach_pos: [f64; 3],
        rigid_vel: [f64; 3],
        rigid_inv_mass: f64,
        baumgarte: f64,
    ) -> [f64; 3] {
        self.delta_lambda_for_step(
            deformable,
            rigid_attach_pos,
            rigid_vel,
            rigid_inv_mass,
            baumgarte,
            self.dt,
        )
    }

    /// Relax this coupling once.
    ///
    /// Computes the impulse `P` (see [`Self::compute_delta_lambda`]) with the
    /// step `dt` (falling back to [`RigidDeformableCoupling::dt`] when `dt` is
    /// not a positive finite number), applies it to the mesh, applies the
    /// reaction to the rigid body (`rigid_vel -= rigid_inv_mass * P`) and adds
    /// `P` to [`RigidDeformableCoupling::lambda`].  Returns `P`.
    pub fn solve_iteration<D: DeformableBodyState>(
        &mut self,
        deformable: &mut D,
        rigid_attach_pos: [f64; 3],
        rigid_vel: &mut [f64; 3],
        rigid_inv_mass: f64,
        baumgarte: f64,
        dt: f64,
    ) -> [f64; 3] {
        let step = if dt.is_finite() && dt > 0.0 {
            dt
        } else {
            self.dt
        };
        let impulse = self.delta_lambda_for_step(
            &*deformable,
            rigid_attach_pos,
            *rigid_vel,
            rigid_inv_mass,
            baumgarte,
            step,
        );
        self.apply_to_deformable(deformable, impulse);
        *rigid_vel = sub3(*rigid_vel, scale3(impulse, rigid_inv_mass.max(0.0)));
        self.lambda = add3(self.lambda, impulse);
        impulse
    }

    /// Distribute a coupling impulse over the mesh nodes.
    ///
    /// Node `i` receives the impulse `w_i * impulse` (`w_i = 1` for a
    /// single-node coupling).  The impulse is passed *unscaled*: converting it
    /// to a velocity change with the nodal inverse mass is the job of
    /// [`DeformableBodyState::apply_impulse`].
    pub fn apply_to_deformable<D: DeformableBodyState>(
        &self,
        deformable: &mut D,
        impulse: [f64; 3],
    ) {
        for (node, weight) in self.weighted_nodes().into_iter().flatten() {
            deformable.apply_impulse(node, scale3(impulse, weight));
        }
    }

    // ── private helpers ───────────────────────────────────────────────────────

    /// Coupling impulse for an explicit time step `dt`.
    fn delta_lambda_for_step<D: DeformableBodyState>(
        &self,
        deformable: &D,
        rigid_attach_pos: [f64; 3],
        rigid_vel: [f64; 3],
        rigid_inv_mass: f64,
        baumgarte: f64,
        dt: f64,
    ) -> [f64; 3] {
        let effective_inv_mass = rigid_inv_mass.max(0.0) + self.attachment_inv_mass(deformable);
        if effective_inv_mass <= MIN_INV_MASS {
            return [0.0; 3];
        }
        // Velocity error Cdot and position error C of the point-to-point constraint.
        let rel_vel = sub3(rigid_vel, self.interpolated_velocity(deformable));
        let pos_err = sub3(rigid_attach_pos, self.interpolated_position(deformable));
        let target = add3(rel_vel, scale3(pos_err, baumgarte / dt.max(MIN_DT)));
        scale3(target, self.stiffness / effective_inv_mass)
    }

    /// Mesh nodes taking part in the coupling together with their weights.
    fn weighted_nodes(&self) -> [Option<(DeformableNodeId, f64)>; 4] {
        match self.kind {
            CouplingKind::SingleNode => [Some((self.primary_node, 1.0)), None, None, None],
            CouplingKind::Barycentric(weights) => {
                let nodes = [
                    self.primary_node,
                    self.secondary_nodes[0],
                    self.secondary_nodes[1],
                    self.secondary_nodes[2],
                ];
                let mut out = [None; 4];
                for ((slot, node), &weight) in out.iter_mut().zip(nodes.iter()).zip(weights.iter())
                {
                    if node.0 != usize::MAX {
                        *slot = Some((*node, weight));
                    }
                }
                out
            }
        }
    }

    /// Weight-interpolated world position of the coupled mesh nodes.
    fn interpolated_position<D: DeformableBodyState>(&self, d: &D) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        for (node, weight) in self.weighted_nodes().into_iter().flatten() {
            out = add3(out, scale3(d.node_position(node), weight));
        }
        out
    }

    /// Weight-interpolated world velocity of the coupled mesh nodes.
    fn interpolated_velocity<D: DeformableBodyState>(&self, d: &D) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        for (node, weight) in self.weighted_nodes().into_iter().flatten() {
            out = add3(out, scale3(d.node_velocity(node), weight));
        }
        out
    }

    /// Mesh-side effective inverse mass `J M^-1 J^T = sum_i(w_i^2 * inv_mass_i)`
    /// of the interpolated attachment point.
    fn attachment_inv_mass<D: DeformableBodyState>(&self, d: &D) -> f64 {
        self.weighted_nodes()
            .into_iter()
            .flatten()
            .map(|(node, weight)| weight * weight * d.node_inv_mass(node))
            .sum()
    }
}

// ── CoupledSimulation ────────────────────────────────────────────────────────

/// Manages multiple [`RigidDeformableCoupling`] constraints and drives a
/// combined velocity-level solve each physics step.
///
/// The caller is responsible for running the rigid-only PGS solver first,
/// then calling [`CoupledSimulation::solve_velocity`] to exchange coupling
/// impulses between the two solvers: the mesh receives every coupling impulse
/// and the matching rigid body receives the reaction, so momentum is conserved.
#[derive(Debug)]
pub struct CoupledSimulation {
    /// All active coupling constraints.
    pub couplings: Vec<RigidDeformableCoupling>,
    /// Maximum coupling iterations per step.
    pub max_iterations: usize,
}

impl Default for CoupledSimulation {
    /// An empty coupler with a budget of 10 coupling iterations per step.
    fn default() -> Self {
        Self::new(DEFAULT_MAX_ITERATIONS)
    }
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

    /// Velocity-level coupling solve (Gauss–Seidel over all couplings).
    ///
    /// Up to [`CoupledSimulation::max_iterations`] sweeps are made over the
    /// couplings; the solve stops early once a sweep applies no impulse.  Each
    /// coupling uses the rigid body selected by its
    /// [`RigidDeformableCoupling::rigid_body`] index: the impulse `P` computed
    /// from the current rigid and mesh state is applied to the mesh and its
    /// reaction is applied to the rigid body
    /// (`rigid_vels[i] -= rigid_inv_masses[i] * P`).  Couplings whose body
    /// index lies outside the supplied arrays are skipped.  Every coupling's
    /// [`RigidDeformableCoupling::lambda`] is reset at the start and holds the
    /// total impulse it applied to the mesh afterwards.
    ///
    /// Rigid bodies are translating point masses (see the module
    /// documentation): the world attachment point of a coupling is
    /// `rigid_positions[i] + local_attachment`.
    ///
    /// * `rigid_positions`  — world-space position of each rigid body (by index)
    /// * `rigid_vels`       — world-space linear velocity of each rigid body
    ///   (by index), updated in place with the reaction impulses
    /// * `rigid_inv_masses` — inverse mass of each rigid body (by index)
    /// * `deformable`       — deformable body state (mutably)
    /// * `dt`               — time step of this solve; the Baumgarte term uses
    ///   `baumgarte / dt` (a non-positive or non-finite value falls back to
    ///   each coupling's own `dt`)
    /// * `baumgarte`        — position-error feedback factor (`0` = match
    ///   velocities only)
    pub fn solve_velocity<D: DeformableBodyState>(
        &mut self,
        rigid_positions: &[[f64; 3]],
        rigid_vels: &mut [[f64; 3]],
        rigid_inv_masses: &[f64],
        deformable: &mut D,
        dt: f64,
        baumgarte: f64,
    ) {
        for c in &mut self.couplings {
            c.lambda = [0.0; 3];
        }
        for _ in 0..self.max_iterations {
            let mut largest_impulse = 0.0_f64;
            for c in &mut self.couplings {
                let body = c.rigid_body;
                let (Some(&position), Some(&inv_mass), Some(vel)) = (
                    rigid_positions.get(body),
                    rigid_inv_masses.get(body),
                    rigid_vels.get_mut(body),
                ) else {
                    continue;
                };
                let attach = c.attachment_world_position(position);
                let impulse = c.solve_iteration(deformable, attach, vel, inv_mass, baumgarte, dt);
                for component in impulse {
                    largest_impulse = largest_impulse.max(component.abs());
                }
            }
            if largest_impulse <= CONVERGED_IMPULSE {
                break;
            }
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

    /// Mesh with per-node state: `apply_impulse` changes the node velocity by
    /// `inv_mass * impulse`, exactly like a real solver would.
    struct MockMesh {
        pos: Vec<[f64; 3]>,
        vel: Vec<[f64; 3]>,
        inv_mass: Vec<f64>,
        applied: Vec<(usize, [f64; 3])>,
    }

    impl MockMesh {
        fn uniform(nodes: usize, inv_mass: f64) -> Self {
            Self {
                pos: vec![[0.0; 3]; nodes],
                vel: vec![[0.0; 3]; nodes],
                inv_mass: vec![inv_mass; nodes],
                applied: Vec::new(),
            }
        }
    }

    impl DeformableBodyState for MockMesh {
        fn node_position(&self, id: DeformableNodeId) -> [f64; 3] {
            self.pos[id.0]
        }
        fn node_velocity(&self, id: DeformableNodeId) -> [f64; 3] {
            self.vel[id.0]
        }
        fn apply_impulse(&mut self, id: DeformableNodeId, impulse: [f64; 3]) {
            self.applied.push((id.0, impulse));
            let im = self.inv_mass[id.0];
            for (v, p) in self.vel[id.0].iter_mut().zip(impulse.iter()) {
                *v += im * p;
            }
        }
        fn node_inv_mass(&self, id: DeformableNodeId) -> f64 {
            self.inv_mass[id.0]
        }
    }

    fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
        for k in 0..3 {
            assert!(
                (actual[k] - expected[k]).abs() <= tol,
                "component {k}: {actual:?} != {expected:?}"
            );
        }
    }

    fn single(node: usize) -> RigidDeformableCoupling {
        RigidDeformableCoupling::new(DeformableNodeId(node), [0.0; 3], 1.0 / 60.0)
    }

    #[test]
    fn compute_delta_lambda_zero_when_state_matches() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        mesh.vel[0] = [1.0, 0.0, 0.0];
        mesh.pos[0] = [0.5, 0.0, 0.0];
        let dl = single(0).compute_delta_lambda(&mesh, [0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.2);
        assert_vec_close(dl, [0.0; 3], 1e-12);
    }

    #[test]
    fn compute_delta_lambda_drags_mesh_to_rigid_velocity() {
        // Kinematic rigid body (inv_mass 0), unit-mass node, no Baumgarte:
        // one relaxation must bring the node exactly to the rigid velocity.
        let mut mesh = MockMesh::uniform(1, 1.0);
        let c = single(0);
        let p = c.compute_delta_lambda(&mesh, [0.0; 3], [1.0, 0.0, 0.0], 0.0, 0.0);
        assert_vec_close(p, [1.0, 0.0, 0.0], 1e-12);
        c.apply_to_deformable(&mut mesh, p);
        assert_vec_close(mesh.vel[0], [1.0, 0.0, 0.0], 1e-12);
    }

    #[test]
    fn baumgarte_does_not_distort_the_velocity_impulse() {
        // Matching positions: the position bias is zero, so the impulse is
        // exactly eff * rv whatever beta and dt are (the old code multiplied
        // it by `1 - beta / dt`).
        let mut mesh = MockMesh::uniform(1, 2.0);
        mesh.vel[0] = [0.0, -1.0, 0.0];
        let c = single(0);
        let p = c.compute_delta_lambda(&mesh, [0.0; 3], [0.0, 2.0, 0.0], 0.0, 0.2);
        // rv = 3, K = 2 -> P = 1.5
        assert_vec_close(p, [0.0, 1.5, 0.0], 1e-12);
    }

    #[test]
    fn baumgarte_pulls_mesh_toward_the_attachment_point() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut c = single(0);
        c.dt = 0.1;
        // Rigid attachment 0.1 m ahead along +x, velocities equal (zero).
        let p = c.compute_delta_lambda(&mesh, [0.1, 0.0, 0.0], [0.0; 3], 0.0, 0.2);
        // P = K^-1 * (beta / dt) * C = 1 * 2 * 0.1
        assert_vec_close(p, [0.2, 0.0, 0.0], 1e-12);
        c.apply_to_deformable(&mut mesh, p);
        assert!(
            mesh.vel[0][0] > 0.0,
            "node must move toward the rigid point"
        );
    }

    #[test]
    fn stiffness_scales_the_impulse() {
        let mesh = MockMesh::uniform(1, 1.0);
        let mut c = single(0);
        c.stiffness = 0.25;
        let p = c.compute_delta_lambda(&mesh, [0.0; 3], [4.0, 0.0, 0.0], 0.0, 0.0);
        assert_vec_close(p, [1.0, 0.0, 0.0], 1e-12);
    }

    #[test]
    fn effective_mass_uses_squared_barycentric_weights() {
        // After one relaxation against a kinematic body the interpolated mesh
        // velocity must equal the rigid velocity exactly; this only holds with
        // K = sum(w^2 * inv_mass) (not sum(w * inv_mass)).
        let mut mesh = MockMesh::uniform(4, 1.0);
        mesh.inv_mass = vec![1.0, 2.0, 0.5, 4.0];
        let weights = [0.5, 0.25, 0.125, 0.125];
        let nodes = [
            DeformableNodeId(0),
            DeformableNodeId(1),
            DeformableNodeId(2),
            DeformableNodeId(3),
        ];
        let c = RigidDeformableCoupling::new_barycentric(nodes, weights, [0.0; 3], 1.0 / 60.0);
        let rigid_vel = [0.0, 0.0, 2.0];
        let p = c.compute_delta_lambda(&mesh, [0.0; 3], rigid_vel, 0.0, 0.0);
        c.apply_to_deformable(&mut mesh, p);
        let mut interpolated = [0.0; 3];
        for (node, w) in weights.iter().enumerate() {
            for (acc, v) in interpolated.iter_mut().zip(mesh.vel[node].iter()) {
                *acc += w * v;
            }
        }
        assert_vec_close(interpolated, rigid_vel, 1e-12);
    }

    #[test]
    fn solve_iteration_applies_reaction_and_conserves_momentum() {
        // Rigid mass 2 at 3 m/s hits a unit-mass node at rest: perfectly
        // inelastic, both end at 2 m/s and the momentum 6 is conserved.
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut c = single(0);
        let mut rigid_vel = [3.0, 0.0, 0.0];
        let p = c.solve_iteration(&mut mesh, [0.0; 3], &mut rigid_vel, 0.5, 0.0, 1.0 / 60.0);
        assert_vec_close(p, [2.0, 0.0, 0.0], 1e-12);
        assert_vec_close(mesh.vel[0], [2.0, 0.0, 0.0], 1e-12);
        assert_vec_close(rigid_vel, [2.0, 0.0, 0.0], 1e-12);
        let momentum = 2.0 * rigid_vel[0] + mesh.vel[0][0];
        assert!((momentum - 6.0).abs() < 1e-12, "momentum = {momentum}");
        assert_vec_close(c.lambda, p, 1e-15);
    }

    #[test]
    fn immovable_pair_produces_no_impulse() {
        let mesh = MockMesh::uniform(1, 0.0);
        let p = single(0).compute_delta_lambda(&mesh, [1.0; 3], [1.0; 3], 0.0, 0.3);
        assert_vec_close(p, [0.0; 3], 0.0);
    }

    #[test]
    fn apply_to_deformable_passes_the_impulse_unscaled() {
        // The trait takes an impulse; the node converts it with its own
        // inverse mass (the mock does `vel += inv_mass * impulse`).
        let mut mesh = MockMesh::uniform(1, 2.0);
        single(0).apply_to_deformable(&mut mesh, [0.0, 10.0, 0.0]);
        assert_eq!(mesh.applied.len(), 1);
        assert_vec_close(mesh.applied[0].1, [0.0, 10.0, 0.0], 1e-15);
        assert_vec_close(mesh.vel[0], [0.0, 20.0, 0.0], 1e-12);
    }

    #[test]
    fn barycentric_apply_splits_impulse_by_weights() {
        let mut mesh = MockMesh::uniform(4, 1.0);
        let nodes = [
            DeformableNodeId(0),
            DeformableNodeId(1),
            DeformableNodeId(2),
            DeformableNodeId(3),
        ];
        let weights = [0.4, 0.3, 0.2, 0.1];
        let c = RigidDeformableCoupling::new_barycentric(nodes, weights, [0.0; 3], 1.0 / 60.0);
        c.apply_to_deformable(&mut mesh, [0.0, 10.0, 0.0]);
        assert_eq!(mesh.applied.len(), 4);
        let mut total = 0.0;
        for (i, (node, impulse)) in mesh.applied.iter().enumerate() {
            assert_eq!(*node, i);
            assert!((impulse[1] - 10.0 * weights[i]).abs() < 1e-12);
            total += impulse[1];
        }
        assert!((total - 10.0).abs() < 1e-12, "impulse must be conserved");
    }

    #[test]
    fn barycentric_coupling_skips_unused_node_slots() {
        // The mock panics on an out-of-range node id, so a skipped slot proves
        // the sentinel is honoured by every code path.
        let mut mesh = MockMesh::uniform(3, 1.0);
        let nodes = [
            DeformableNodeId(0),
            DeformableNodeId(1),
            DeformableNodeId(2),
            DeformableNodeId(usize::MAX),
        ];
        let c =
            RigidDeformableCoupling::new_barycentric(nodes, [0.5, 0.25, 0.25, 0.0], [0.0; 3], 0.1);
        let p = c.compute_delta_lambda(&mesh, [0.0; 3], [1.0, 0.0, 0.0], 0.0, 0.0);
        c.apply_to_deformable(&mut mesh, p);
        assert_eq!(mesh.applied.len(), 3);
    }

    #[test]
    fn attachment_offset_enters_the_position_error() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        mesh.pos[0] = [0.0, 1.0, 0.0];
        let c = RigidDeformableCoupling::new(DeformableNodeId(0), [0.0, 1.0, 0.0], 1.0 / 60.0);
        let attach = c.attachment_world_position([0.0; 3]);
        assert_vec_close(attach, [0.0, 1.0, 0.0], 0.0);
        let p = c.compute_delta_lambda(&mesh, attach, [0.0; 3], 0.0, 0.2);
        assert_vec_close(p, [0.0; 3], 1e-12);
    }

    #[test]
    fn with_rigid_body_selects_the_body() {
        assert_eq!(single(0).rigid_body, 0);
        assert_eq!(single(0).with_rigid_body(3).rigid_body, 3);
    }

    #[test]
    fn solve_velocity_couples_each_coupling_to_its_own_body() {
        // Two nodes, two kinematic bodies with different velocities: every node
        // must follow its own body (the old code used body 0 for everything).
        let mut mesh = MockMesh::uniform(2, 1.0);
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(single(0).with_rigid_body(0));
        sim.add_coupling(single(1).with_rigid_body(1));
        let positions = [[0.0; 3]; 2];
        let mut vels = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let inv_masses = [0.0, 0.0];
        sim.solve_velocity(
            &positions,
            &mut vels,
            &inv_masses,
            &mut mesh,
            1.0 / 60.0,
            0.0,
        );
        assert_vec_close(mesh.vel[0], [1.0, 0.0, 0.0], 1e-12);
        assert_vec_close(mesh.vel[1], [0.0, 2.0, 0.0], 1e-12);
        // Kinematic bodies are unaffected by the reaction.
        assert_vec_close(vels[0], [1.0, 0.0, 0.0], 0.0);
        assert_vec_close(vels[1], [0.0, 2.0, 0.0], 0.0);
    }

    #[test]
    fn solve_velocity_is_two_way_and_conserves_momentum_per_pair() {
        let mut mesh = MockMesh::uniform(2, 1.0);
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(single(0).with_rigid_body(0));
        sim.add_coupling(single(1).with_rigid_body(1));
        let positions = [[0.0; 3]; 2];
        let mut vels = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        sim.solve_velocity(
            &positions,
            &mut vels,
            &inv_masses,
            &mut mesh,
            1.0 / 60.0,
            0.0,
        );
        // Unit masses: each pair ends at its common mean velocity.
        assert_vec_close(vels[0], [0.5, 0.0, 0.0], 1e-12);
        assert_vec_close(mesh.vel[0], [0.5, 0.0, 0.0], 1e-12);
        assert_vec_close(vels[1], [0.0, 1.0, 0.0], 1e-12);
        assert_vec_close(mesh.vel[1], [0.0, 1.0, 0.0], 1e-12);
    }

    #[test]
    fn solve_velocity_iterates_for_a_shared_node() {
        // One node held by two unit-mass bodies (2 m/s and 0 m/s): the exact
        // solution is a common velocity of 2/3 and needs several sweeps.
        let build = |iterations: usize| {
            let mut mesh = MockMesh::uniform(1, 1.0);
            let mut sim = CoupledSimulation::new(iterations);
            sim.add_coupling(single(0).with_rigid_body(0));
            sim.add_coupling(single(0).with_rigid_body(1));
            let positions = [[0.0; 3]; 2];
            let mut vels = [[2.0, 0.0, 0.0], [0.0; 3]];
            sim.solve_velocity(
                &positions,
                &mut vels,
                &[1.0, 1.0],
                &mut mesh,
                1.0 / 60.0,
                0.0,
            );
            (mesh.vel[0][0], vels[0][0], vels[1][0])
        };
        let (node, a, b) = build(200);
        for v in [node, a, b] {
            assert!((v - 2.0 / 3.0).abs() < 1e-9, "v = {v}");
        }
        assert!((node + a + b - 2.0).abs() < 1e-12, "momentum not conserved");
        // A single sweep is visibly short of convergence: max_iterations is honoured.
        let (node1, a1, b1) = build(1);
        let spread = a1.max(b1).max(node1) - a1.min(b1).min(node1);
        assert!(spread > 0.4, "one sweep cannot converge, spread = {spread}");
    }

    #[test]
    fn solve_velocity_zero_iterations_changes_nothing() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut sim = CoupledSimulation::new(0);
        sim.add_coupling(single(0));
        let mut vels = [[1.0, 0.0, 0.0]];
        sim.solve_velocity(&[[0.0; 3]], &mut vels, &[0.0], &mut mesh, 1.0 / 60.0, 0.2);
        assert!(mesh.applied.is_empty());
        assert_vec_close(vels[0], [1.0, 0.0, 0.0], 0.0);
    }

    #[test]
    fn solve_velocity_skips_couplings_without_a_body() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(single(0).with_rigid_body(5));
        let mut vels = [[1.0, 0.0, 0.0]];
        sim.solve_velocity(&[[0.0; 3]], &mut vels, &[0.0], &mut mesh, 1.0 / 60.0, 0.2);
        assert!(mesh.applied.is_empty());
    }

    #[test]
    fn lambda_accumulates_mesh_impulse_and_resets_every_solve() {
        let mut mesh = MockMesh::uniform(1, 0.5);
        let mut sim = CoupledSimulation::new(6);
        sim.add_coupling(single(0));
        let mut vels = [[0.0, 3.0, 0.0]];
        sim.solve_velocity(&[[0.0; 3]], &mut vels, &[0.0], &mut mesh, 1.0 / 60.0, 0.0);
        // Node mass 2 ends at 3 m/s: total impulse 6 along y.
        assert_vec_close(sim.couplings[0].lambda, [0.0, 6.0, 0.0], 1e-12);
        // Already satisfied: the second solve applies nothing and lambda restarts.
        sim.solve_velocity(&[[0.0; 3]], &mut vels, &[0.0], &mut mesh, 1.0 / 60.0, 0.0);
        assert_vec_close(sim.couplings[0].lambda, [0.0; 3], 1e-12);
    }

    #[test]
    fn solve_velocity_uses_the_step_dt_for_the_baumgarte_term() {
        // The coupling was built for 1/60 s but is solved with a 0.1 s step.
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(single(0));
        let mut vels = [[0.0; 3]];
        let positions = [[0.05, 0.0, 0.0]];
        sim.solve_velocity(&positions, &mut vels, &[0.0], &mut mesh, 0.1, 0.2);
        // Target relative speed = (beta / dt) * C = 2 * 0.05.
        assert_vec_close(mesh.vel[0], [0.1, 0.0, 0.0], 1e-12);
    }

    #[test]
    fn solve_velocity_falls_back_to_the_coupling_dt() {
        let mut mesh = MockMesh::uniform(1, 1.0);
        let mut sim = CoupledSimulation::new(4);
        let mut c = single(0);
        c.dt = 0.1;
        sim.add_coupling(c);
        let mut vels = [[0.0; 3]];
        let positions = [[0.05, 0.0, 0.0]];
        sim.solve_velocity(&positions, &mut vels, &[0.0], &mut mesh, 0.0, 0.2);
        assert_vec_close(mesh.vel[0], [0.1, 0.0, 0.0], 1e-12);
    }

    #[test]
    fn coupled_node_converges_onto_a_moving_rigid_body() {
        // Frame loop: the kinematic body translates at 1 m/s, the mesh node
        // starts 1 m behind.  With Baumgarte feedback the position error
        // shrinks by (1 - beta) every frame.
        let dt = 1.0 / 60.0;
        let mut mesh = MockMesh::uniform(1, 1.0);
        mesh.pos[0] = [-1.0, 0.0, 0.0];
        let mut sim = CoupledSimulation::new(4);
        sim.add_coupling(single(0));
        let mut body_pos = [[0.0, 0.0, 0.0]];
        let mut body_vel = [[1.0, 0.0, 0.0]];
        for _ in 0..150 {
            sim.solve_velocity(&body_pos, &mut body_vel, &[0.0], &mut mesh, dt, 0.2);
            for k in 0..3 {
                body_pos[0][k] += dt * body_vel[0][k];
                mesh.pos[0][k] += dt * mesh.vel[0][k];
            }
        }
        let gap = (body_pos[0][0] - mesh.pos[0][0]).abs();
        assert!(gap < 1e-9, "node must track the rigid body, gap = {gap}");
        assert!((mesh.vel[0][0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn default_simulation_has_an_iteration_budget() {
        let sim = CoupledSimulation::default();
        assert!(sim.couplings.is_empty());
        assert_eq!(sim.max_iterations, DEFAULT_MAX_ITERATIONS);
    }
}
