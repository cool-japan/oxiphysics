// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body – soft body coupling for hybrid simulation.
//!
//! Provides:
//! - Minimal rigid body (position, velocity, orientation quaternion, angular velocity)
//! - Soft body mesh (nodes + velocities + masses)
//! - Contact detection and penalty-force response
//! - Embedded rigid body constraints (rigid core inside a soft shell)

#![allow(dead_code)]

// ─── Vector helpers ──────────────────────────────────────────────────────────

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_len(v: [f64; 3]) -> f64 {
    v3_dot(v, v).sqrt()
}

#[inline]
fn v3_norm(v: [f64; 3]) -> [f64; 3] {
    let l = v3_len(v);
    if l < 1e-30 {
        [0.0, 0.0, 0.0]
    } else {
        v3_scale(v, 1.0 / l)
    }
}

// ─── RigidBody ───────────────────────────────────────────────────────────────

/// Minimal rigid body for coupling with a soft body mesh.
///
/// Orientation is stored as a unit quaternion `[qx, qy, qz, qw]`.
/// Inertia tensor is stored as a flat row-major 3×3 matrix.
#[derive(Debug, Clone)]
pub struct RigidBody {
    /// World-space position of the centre of mass \[m\].
    pub position: [f64; 3],
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Orientation quaternion \[qx, qy, qz, qw\] (unit).
    pub orientation: [f64; 4],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Total mass \[kg\].
    pub mass: f64,
    /// Body-space inertia tensor (row-major 3×3) \[kg·m²\].
    pub inertia: [f64; 9],
}

impl RigidBody {
    /// Construct a new `RigidBody` at rest.
    pub fn new(position: [f64; 3], mass: f64, inertia: [f64; 9]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            mass,
            inertia,
        }
    }

    /// Construct a uniform sphere with moment of inertia `2/5 m r²`.
    pub fn sphere(position: [f64; 3], mass: f64, radius: f64) -> Self {
        let i = 2.0 / 5.0 * mass * radius * radius;
        Self::new(position, mass, [i, 0.0, 0.0, 0.0, i, 0.0, 0.0, 0.0, i])
    }

    /// Integrate the rigid body forward by `dt` seconds under a net force and torque.
    pub fn integrate(&mut self, force: [f64; 3], torque: [f64; 3], dt: f64) {
        // Linear motion.
        let a = v3_scale(force, 1.0 / self.mass.max(1e-30));
        self.velocity = v3_add(self.velocity, v3_scale(a, dt));
        self.position = v3_add(self.position, v3_scale(self.velocity, dt));

        // Rotational motion (simple explicit Euler on angular velocity).
        // Use diagonal approximation of inertia inverse.
        let ixx = self.inertia[0].max(1e-30);
        let iyy = self.inertia[4].max(1e-30);
        let izz = self.inertia[8].max(1e-30);
        let alpha = [torque[0] / ixx, torque[1] / iyy, torque[2] / izz];
        self.angular_velocity = v3_add(self.angular_velocity, v3_scale(alpha, dt));

        // Integrate quaternion: q' = q + 0.5 × Ω⊗q × dt
        let [ox, oy, oz] = self.angular_velocity;
        let [qx, qy, qz, qw] = self.orientation;
        let dqx = 0.5 * (qw * ox - qz * oy + qy * oz);
        let dqy = 0.5 * (qz * ox + qw * oy - qx * oz);
        let dqz = 0.5 * (-qy * ox + qx * oy + qw * oz);
        let dqw = 0.5 * (-qx * ox - qy * oy - qz * oz);
        let mut q = [qx + dqx * dt, qy + dqy * dt, qz + dqz * dt, qw + dqw * dt];
        // Normalise quaternion.
        let qlen = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        if qlen > 1e-30 {
            q[0] /= qlen;
            q[1] /= qlen;
            q[2] /= qlen;
            q[3] /= qlen;
        }
        self.orientation = q;
    }
}

// ─── SoftBodyMesh ────────────────────────────────────────────────────────────

/// Soft body represented as a collection of mass-weighted nodes.
#[derive(Debug, Clone)]
pub struct SoftBodyMesh {
    /// Node positions \[m\].
    pub nodes: Vec<[f64; 3]>,
    /// Node velocities \[m/s\].
    pub velocities: Vec<[f64; 3]>,
    /// Node masses \[kg\].
    pub masses: Vec<f64>,
}

impl SoftBodyMesh {
    /// Construct a new `SoftBodyMesh`.
    pub fn new(nodes: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            velocities: vec![[0.0; 3]; n],
            masses,
        }
    }

    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Apply gravity to all nodes for one time step.
    pub fn apply_gravity(&mut self, g: [f64; 3], dt: f64) {
        for (v, _m) in self.velocities.iter_mut().zip(self.masses.iter()) {
            *v = v3_add(*v, v3_scale(g, dt));
        }
        for (p, v) in self.nodes.iter_mut().zip(self.velocities.iter()) {
            *p = v3_add(*p, v3_scale(*v, dt));
        }
    }
}

// ─── ContactPoint ────────────────────────────────────────────────────────────

/// A detected contact between a rigid body surface and a soft body node.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// World-space contact position \[m\].
    pub position: [f64; 3],
    /// Contact normal pointing from rigid surface into soft body \[unit vector\].
    pub normal: [f64; 3],
    /// Penetration depth \[m\] (positive = overlapping).
    pub penetration: f64,
    /// Index of the soft body node involved in this contact.
    pub soft_node: usize,
}

// ─── RigidSoftContact ────────────────────────────────────────────────────────

/// Manages contact between a `RigidBody` and a `SoftBodyMesh`.
///
/// Uses a sphere-vs-nodes broadphase: any soft node inside the rigid sphere
/// is considered a contact.
#[derive(Debug, Clone)]
pub struct RigidSoftContact {
    /// The rigid body acting as the collision shape.
    pub rigid: RigidBody,
    /// The soft body mesh.
    pub soft: SoftBodyMesh,
    /// Contact stiffness coefficient \[N/m\].
    pub contact_stiffness: f64,
    /// Contact damping coefficient \[N·s/m\].
    pub contact_damping: f64,
    /// Radius of the rigid sphere used for broadphase contact detection \[m\].
    pub rigid_radius: f64,
}

impl RigidSoftContact {
    /// Construct a new `RigidSoftContact`.
    pub fn new(
        rigid: RigidBody,
        soft: SoftBodyMesh,
        contact_stiffness: f64,
        contact_damping: f64,
        rigid_radius: f64,
    ) -> Self {
        Self {
            rigid,
            soft,
            contact_stiffness,
            contact_damping,
            rigid_radius,
        }
    }

    /// Detect contact points between the rigid sphere and soft nodes.
    ///
    /// Returns a list of `ContactPoint` for every soft node that penetrates the
    /// rigid sphere.
    pub fn detect_contacts(&self) -> Vec<ContactPoint> {
        let mut contacts = Vec::new();
        for (i, &node_pos) in self.soft.nodes.iter().enumerate() {
            let diff = v3_sub(node_pos, self.rigid.position);
            let dist = v3_len(diff);
            let penetration = self.rigid_radius - dist;
            if penetration > 0.0 {
                let normal = if dist > 1e-12 {
                    v3_scale(diff, 1.0 / dist)
                } else {
                    [0.0, 1.0, 0.0] // fallback normal
                };
                let contact_pos = v3_add(self.rigid.position, v3_scale(normal, self.rigid_radius));
                contacts.push(ContactPoint {
                    position: contact_pos,
                    normal,
                    penetration,
                    soft_node: i,
                });
            }
        }
        contacts
    }

    /// Apply penalty contact forces to soft nodes and impulse to rigid body.
    ///
    /// Each contact contributes a force `F = penalty_force(penetration, stiffness)`.
    pub fn apply_contact_forces(&mut self, contacts: &[ContactPoint], dt: f64) {
        for cp in contacts {
            let idx = cp.soft_node;
            let m_soft = self.soft.masses[idx].max(1e-30);

            // Relative velocity of soft node w.r.t. rigid surface.
            let rel_vel = v3_sub(self.soft.velocities[idx], self.rigid.velocity);
            let rel_vel_normal = v3_dot(rel_vel, cp.normal);

            let f_pen = penalty_force(cp.penetration, self.contact_stiffness);
            let f_damp = -self.contact_damping * rel_vel_normal;
            let f_total = f_pen + f_damp;

            // Apply force to soft node (impulse = F * dt / m).
            let impulse_soft = v3_scale(cp.normal, f_total * dt / m_soft);
            self.soft.velocities[idx] = v3_add(self.soft.velocities[idx], impulse_soft);

            // Reaction on rigid body.
            let m_rigid = self.rigid.mass.max(1e-30);
            let impulse_rigid = v3_scale(cp.normal, -f_total * dt / m_rigid);
            self.rigid.velocity = v3_add(self.rigid.velocity, impulse_rigid);
        }
    }
}

// ─── EmbeddedRigidBody ───────────────────────────────────────────────────────

/// A rigid body embedded inside a soft body mesh.
///
/// Constraint nodes are those soft nodes attached to the rigid body frame.
/// The constraint enforces that each attached node follows the rigid body's
/// motion exactly.
#[derive(Debug, Clone)]
pub struct EmbeddedRigidBody {
    /// The embedded rigid body.
    pub rigid: RigidBody,
    /// Indices of soft mesh nodes that are attached to (rigidly constrained by) this body.
    pub constraint_nodes: Vec<usize>,
    /// Rest offsets of constraint nodes from the rigid body centre, in world space at construction.
    pub rest_offsets: Vec<[f64; 3]>,
}

impl EmbeddedRigidBody {
    /// Construct a new `EmbeddedRigidBody`, recording rest offsets from the mesh.
    pub fn new(rigid: RigidBody, constraint_nodes: Vec<usize>, mesh: &SoftBodyMesh) -> Self {
        let rest_offsets = constraint_nodes
            .iter()
            .map(|&idx| v3_sub(mesh.nodes[idx], rigid.position))
            .collect();
        Self {
            rigid,
            constraint_nodes,
            rest_offsets,
        }
    }

    /// Apply rigid-body constraints: move constrained soft nodes to follow the rigid body.
    ///
    /// The velocity correction blends node velocity toward the rigid body velocity
    /// with a stiffness factor in \[0, 1\].  `dt` is the time step \[s\].
    #[allow(clippy::too_many_arguments)]
    pub fn apply_constraints(&self, soft: &mut SoftBodyMesh, _dt: f64) {
        // Use the stored orientation quaternion to rotate rest offsets.
        let [qx, qy, qz, qw] = self.rigid.orientation;
        for (&node_idx, &offset) in self.constraint_nodes.iter().zip(self.rest_offsets.iter()) {
            // Rotate offset by quaternion.
            let rotated = quat_rotate_vec([qx, qy, qz, qw], offset);
            let target_pos = v3_add(self.rigid.position, rotated);
            soft.nodes[node_idx] = target_pos;
            soft.velocities[node_idx] = self.rigid.velocity;
        }
    }
}

/// Rotate a 3-vector by a unit quaternion `[qx, qy, qz, qw]`.
fn quat_rotate_vec(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    let [vx, vy, vz] = v;
    // t = 2 × (q_xyz × v)
    let tx = 2.0 * (qy * vz - qz * vy);
    let ty = 2.0 * (qz * vx - qx * vz);
    let tz = 2.0 * (qx * vy - qy * vx);
    [
        vx + qw * tx + qy * tz - qz * ty,
        vy + qw * ty + qz * tx - qx * tz,
        vz + qw * tz + qx * ty - qy * tx,
    ]
}

// ─── Free functions ──────────────────────────────────────────────────────────

/// Penalty force \[N\] for a given penetration depth and stiffness.
///
/// `F = stiffness × penetration` (linear spring model, clamped to ≥ 0).
pub fn penalty_force(penetration: f64, stiffness: f64) -> f64 {
    (stiffness * penetration).max(0.0)
}

/// Unit contact normal pointing from `rigid_surface_point` toward `soft_node`.
pub fn contact_normal(soft_node: [f64; 3], rigid_surface_point: [f64; 3]) -> [f64; 3] {
    v3_norm(v3_sub(soft_node, rigid_surface_point))
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn sphere_rigid() -> RigidBody {
        RigidBody::sphere([0.0, 0.0, 0.0], 1.0, 0.5)
    }

    fn simple_mesh(nodes: Vec<[f64; 3]>) -> SoftBodyMesh {
        let n = nodes.len();
        SoftBodyMesh::new(nodes, vec![1.0; n])
    }

    // 1. penalty_force: positive penetration returns positive force.
    #[test]
    fn test_penalty_force_positive() {
        let f = penalty_force(0.05, 1000.0);
        assert!(
            f > 0.0,
            "Penalty force should be positive for penetrating contact"
        );
    }

    // 2. penalty_force: zero penetration returns zero.
    #[test]
    fn test_penalty_force_zero_penetration() {
        let f = penalty_force(0.0, 1000.0);
        assert!(
            f.abs() < EPS,
            "Penalty force should be zero at zero penetration"
        );
    }

    // 3. penalty_force: negative penetration (separation) returns zero.
    #[test]
    fn test_penalty_force_negative_penetration() {
        let f = penalty_force(-0.1, 1000.0);
        assert!(
            f.abs() < EPS,
            "Penalty force should be zero when not penetrating"
        );
    }

    // 4. penalty_force: scales linearly with penetration.
    #[test]
    fn test_penalty_force_linear_scaling() {
        let f1 = penalty_force(0.01, 500.0);
        let f2 = penalty_force(0.02, 500.0);
        assert!(
            (f2 - 2.0 * f1).abs() < EPS,
            "Penalty force must scale linearly"
        );
    }

    // 5. penalty_force: scales linearly with stiffness.
    #[test]
    fn test_penalty_force_stiffness_scaling() {
        let f1 = penalty_force(0.05, 100.0);
        let f2 = penalty_force(0.05, 200.0);
        assert!(
            (f2 - 2.0 * f1).abs() < EPS,
            "Penalty force must scale linearly with stiffness"
        );
    }

    // 6. contact_normal: normal points from rigid surface to soft node.
    #[test]
    fn test_contact_normal_direction() {
        let n = contact_normal([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((n[1] - 1.0).abs() < EPS, "Normal should point in +Y");
    }

    // 7. contact_normal: result is unit length.
    #[test]
    fn test_contact_normal_unit_length() {
        let n = contact_normal([3.0, 4.0, 0.0], [0.0, 0.0, 0.0]);
        let len = v3_len(n);
        assert!(
            (len - 1.0).abs() < 1e-10,
            "Contact normal must be unit length, got {len}"
        );
    }

    // 8. contact_normal: degenerate case (same point) returns zero vector.
    #[test]
    fn test_contact_normal_degenerate() {
        let n = contact_normal([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]);
        assert!(v3_len(n) < EPS, "Degenerate contact normal should be zero");
    }

    // 9. detect_contacts: node outside sphere is not detected.
    #[test]
    fn test_detect_contacts_no_penetration() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[1.0, 0.0, 0.0]]); // distance 1.0 > radius 0.5
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        assert!(
            contacts.is_empty(),
            "Node outside sphere should not be a contact"
        );
    }

    // 10. detect_contacts: node inside sphere is detected.
    #[test]
    fn test_detect_contacts_penetration() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.1, 0.0, 0.0]]); // inside radius 0.5
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        assert_eq!(
            contacts.len(),
            1,
            "Node inside sphere should generate a contact"
        );
    }

    // 11. detect_contacts: penetration depth is correct.
    #[test]
    fn test_detect_contacts_penetration_depth() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        let expected_pen = 0.5 - 0.3;
        assert!(
            (contacts[0].penetration - expected_pen).abs() < 1e-10,
            "Penetration depth mismatch"
        );
    }

    // 12. detect_contacts: contact normal points outward from rigid centre.
    #[test]
    fn test_detect_contacts_normal_outward() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        // Normal should have positive x component (pointing away from origin toward node).
        assert!(
            contacts[0].normal[0] > 0.0,
            "Normal should point toward the penetrating node"
        );
    }

    // 13. detect_contacts: multiple nodes, only inside ones detected.
    #[test]
    fn test_detect_contacts_multiple_nodes() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![
            [0.2, 0.0, 0.0],  // inside
            [0.8, 0.0, 0.0],  // outside
            [-0.1, 0.0, 0.0], // inside
        ]);
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        assert_eq!(contacts.len(), 2, "Exactly two nodes should be in contact");
    }

    // 14. apply_contact_forces: soft node velocity increases in normal direction.
    #[test]
    fn test_apply_contact_forces_velocity_change() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let mut coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        coupling.apply_contact_forces(&contacts, 0.01);
        // Soft node should now have positive x velocity (pushed away from sphere).
        assert!(
            coupling.soft.velocities[0][0] > 0.0,
            "Soft node should be pushed away from rigid body"
        );
    }

    // 15. apply_contact_forces: rigid body receives reaction impulse.
    #[test]
    fn test_apply_contact_forces_rigid_reaction() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let mut coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        coupling.apply_contact_forces(&contacts, 0.01);
        // Rigid body should move in negative x (pushed back).
        assert!(
            coupling.rigid.velocity[0] < 0.0,
            "Rigid body should receive reaction impulse"
        );
    }

    // 16. apply_contact_forces: zero stiffness gives no velocity change.
    #[test]
    fn test_apply_contact_forces_zero_stiffness() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let mut coupling = RigidSoftContact::new(rigid, soft, 0.0, 0.0, 0.5);
        let contacts = coupling.detect_contacts();
        coupling.apply_contact_forces(&contacts, 0.01);
        let vx = coupling.soft.velocities[0][0];
        assert!(
            vx.abs() < EPS,
            "Zero stiffness should give no velocity change"
        );
    }

    // 17. EmbeddedRigidBody: constrained node follows rigid position.
    #[test]
    fn test_embedded_constraint_follows_rigid() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.3);
        let mut soft = simple_mesh(vec![[0.2, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let emb = EmbeddedRigidBody::new(rigid, vec![0], &soft);
        // Move rigid to (1, 0, 0).
        let mut emb2 = emb.clone();
        emb2.rigid.position = [1.0, 0.0, 0.0];
        emb2.apply_constraints(&mut soft, 0.01);
        // Node 0 should be near (1 + offset) = (1 + 0.2, 0, 0) = (1.2, 0, 0).
        assert!(
            (soft.nodes[0][0] - 1.2).abs() < 1e-10,
            "Constrained node should follow rigid body, got {:?}",
            soft.nodes[0]
        );
    }

    // 18. EmbeddedRigidBody: unconstrained node is not moved.
    #[test]
    fn test_embedded_constraint_leaves_free_nodes() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.3);
        let mut soft = simple_mesh(vec![[0.2, 0.0, 0.0], [5.0, 3.0, 1.0]]);
        let emb = EmbeddedRigidBody::new(rigid, vec![0], &soft);
        emb.apply_constraints(&mut soft, 0.01);
        // Node 1 should be unchanged.
        assert!(
            (soft.nodes[1][0] - 5.0).abs() < EPS,
            "Free node should not be moved by constraint"
        );
    }

    // 19. EmbeddedRigidBody: constrained node velocity matches rigid velocity.
    #[test]
    fn test_embedded_constraint_velocity() {
        let mut rigid = RigidBody::sphere([0.0; 3], 1.0, 0.3);
        rigid.velocity = [2.0, 0.0, 0.0];
        let mut soft = simple_mesh(vec![[0.1, 0.0, 0.0]]);
        let emb = EmbeddedRigidBody::new(rigid, vec![0], &soft);
        emb.apply_constraints(&mut soft, 0.01);
        assert!(
            (soft.velocities[0][0] - 2.0).abs() < EPS,
            "Constrained node velocity should match rigid velocity"
        );
    }

    // 20. RigidBody: integrate under gravity moves position downward.
    #[test]
    fn test_rigid_body_integrate_gravity() {
        let mut rb = RigidBody::sphere([0.0, 10.0, 0.0], 1.0, 0.5);
        let gravity_force = [0.0, -9.81, 0.0];
        for _ in 0..60 {
            rb.integrate(gravity_force, [0.0; 3], 1.0 / 60.0);
        }
        assert!(rb.position[1] < 9.0, "Rigid body should fall under gravity");
    }

    // 21. RigidBody: no force → no velocity change.
    #[test]
    fn test_rigid_body_no_force() {
        let mut rb = RigidBody::sphere([0.0; 3], 2.0, 0.5);
        rb.integrate([0.0; 3], [0.0; 3], 0.01);
        assert!(
            v3_len(rb.velocity) < EPS,
            "No force should give no velocity"
        );
    }

    // 22. RigidBody: quaternion stays unit after integration with angular velocity.
    #[test]
    fn test_rigid_body_quaternion_unit() {
        let mut rb = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        rb.angular_velocity = [1.0, 0.5, 0.3];
        for _ in 0..100 {
            rb.integrate([0.0; 3], [0.0; 3], 0.01);
        }
        let [qx, qy, qz, qw] = rb.orientation;
        let qlen = (qx * qx + qy * qy + qz * qz + qw * qw).sqrt();
        assert!(
            (qlen - 1.0).abs() < 1e-6,
            "Quaternion must remain unit length, got {qlen}"
        );
    }

    // 23. SoftBodyMesh: apply_gravity accelerates nodes downward.
    #[test]
    fn test_soft_mesh_apply_gravity() {
        let mut mesh = simple_mesh(vec![[0.0, 5.0, 0.0]]);
        let g = [0.0, -9.81, 0.0];
        for _ in 0..60 {
            mesh.apply_gravity(g, 1.0 / 60.0);
        }
        assert!(mesh.nodes[0][1] < 4.0, "Node should fall under gravity");
    }

    // 24. detect_contacts: node on sphere surface (distance = radius) has zero penetration.
    #[test]
    fn test_detect_contacts_surface_node() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let soft = simple_mesh(vec![[0.5, 0.0, 0.0]]); // exactly on surface
        let coupling = RigidSoftContact::new(rigid, soft, 1000.0, 10.0, 0.5);
        let contacts = coupling.detect_contacts();
        // Penetration = 0.5 - 0.5 = 0, so no contact.
        assert!(
            contacts.is_empty(),
            "Node on surface should not have positive penetration"
        );
    }

    // 25. quat_rotate_vec: identity quaternion leaves vector unchanged.
    #[test]
    fn test_quat_rotate_identity() {
        let v = [1.0, 2.0, 3.0];
        let rv = quat_rotate_vec([0.0, 0.0, 0.0, 1.0], v);
        for i in 0..3 {
            assert!(
                (rv[i] - v[i]).abs() < 1e-12,
                "Identity rotation should not change vector"
            );
        }
    }

    // 26. quat_rotate_vec: 90° rotation about Z swaps x and y.
    #[test]
    fn test_quat_rotate_90_deg_z() {
        // Quaternion for 90° about Z: [0, 0, sin(45°), cos(45°)]
        let s = (std::f64::consts::FRAC_PI_4).sin();
        let c = (std::f64::consts::FRAC_PI_4).cos();
        let q = [0.0, 0.0, s, c];
        let rv = quat_rotate_vec(q, [1.0, 0.0, 0.0]);
        // Expected: [0, 1, 0]
        assert!((rv[0]).abs() < 1e-10, "x should map to 0, got {}", rv[0]);
        assert!((rv[1] - 1.0).abs() < 1e-10, "y should be 1, got {}", rv[1]);
        assert!((rv[2]).abs() < 1e-10, "z should be 0, got {}", rv[2]);
    }

    // 27. contact_normal: negative direction.
    #[test]
    fn test_contact_normal_negative() {
        let n = contact_normal([0.0, -1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((n[1] + 1.0).abs() < EPS, "Normal should point in -Y");
    }

    // 28. EmbeddedRigidBody: multiple constraint nodes.
    #[test]
    fn test_embedded_multiple_constraint_nodes() {
        let rigid = RigidBody::sphere([0.0; 3], 1.0, 0.5);
        let mut soft = simple_mesh(vec![[0.3, 0.0, 0.0], [-0.3, 0.0, 0.0], [0.0, 5.0, 0.0]]);
        let emb = EmbeddedRigidBody::new(rigid, vec![0, 1], &soft);
        emb.apply_constraints(&mut soft, 0.01);
        // Node 2 should be untouched.
        assert!(
            (soft.nodes[2][1] - 5.0).abs() < EPS,
            "Free node should not be constrained"
        );
        // Nodes 0 and 1 should be at rigid.position + offset = offset (rigid still at origin).
        assert!(
            (soft.nodes[0][0] - 0.3).abs() < 1e-10,
            "Constrained node 0 should stay at offset"
        );
        assert!(
            (soft.nodes[1][0] + 0.3).abs() < 1e-10,
            "Constrained node 1 should stay at offset"
        );
    }

    // 29. RigidBody sphere inertia tensor is diagonal with correct value.
    #[test]
    fn test_rigid_body_sphere_inertia() {
        let rb = RigidBody::sphere([0.0; 3], 2.0, 1.0);
        let expected = 2.0 / 5.0 * 2.0 * 1.0 * 1.0;
        assert!((rb.inertia[0] - expected).abs() < EPS);
        assert!((rb.inertia[4] - expected).abs() < EPS);
        assert!((rb.inertia[8] - expected).abs() < EPS);
    }

    // 30. apply_contact_forces: momentum is approximately conserved.
    #[test]
    fn test_contact_momentum_conservation() {
        let rigid = RigidBody::sphere([0.0; 3], 2.0, 0.5);
        let soft = simple_mesh(vec![[0.3, 0.0, 0.0]]);
        let mut coupling = RigidSoftContact::new(rigid, soft, 500.0, 0.0, 0.5);
        let p_soft_before: f64 = coupling
            .soft
            .velocities
            .iter()
            .zip(coupling.soft.masses.iter())
            .map(|(v, m)| v[0] * m)
            .sum();
        let p_rigid_before = coupling.rigid.velocity[0] * coupling.rigid.mass;
        let total_before = p_soft_before + p_rigid_before;

        let contacts = coupling.detect_contacts();
        coupling.apply_contact_forces(&contacts, 0.01);

        let p_soft_after: f64 = coupling
            .soft
            .velocities
            .iter()
            .zip(coupling.soft.masses.iter())
            .map(|(v, m)| v[0] * m)
            .sum();
        let p_rigid_after = coupling.rigid.velocity[0] * coupling.rigid.mass;
        let total_after = p_soft_after + p_rigid_after;

        assert!(
            (total_after - total_before).abs() < 1e-10,
            "Momentum should be conserved: before={total_before}, after={total_after}"
        );
    }

    // 31. SoftBodyMesh: num_nodes returns correct count.
    #[test]
    fn test_soft_mesh_num_nodes() {
        let mesh = simple_mesh(vec![[0.0; 3]; 5]);
        assert_eq!(mesh.num_nodes(), 5);
    }

    // 32. penalty_force: large penetration gives large force.
    #[test]
    fn test_penalty_force_large() {
        let f = penalty_force(1.0, 1e6);
        assert!((f - 1e6).abs() < EPS, "Large penetration * stiffness");
    }
}
