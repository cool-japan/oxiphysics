#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated rigid body batch simulation (CPU mock).
//!
//! Provides structs and functions for simulating many rigid bodies in batch,
//! broadphase collision detection (SAP), and constraint solving via sequential
//! impulse. All computation is done on the CPU as a reference implementation.

// ── Quaternion helpers ────────────────────────────────────────────────────────

/// Rotate vector `v` by unit quaternion `q` (format: \[x, y, z, w\]).
///
/// Uses the sandwich product `q * v * q⁻¹` via the Rodrigues formula.
pub fn quat_rotate(q: [f32; 4], v: [f32; 3]) -> [f32; 3] {
    let (qx, qy, qz, qw) = (q[0], q[1], q[2], q[3]);
    // t = 2 * cross(q.xyz, v)
    let tx = 2.0 * (qy * v[2] - qz * v[1]);
    let ty = 2.0 * (qz * v[0] - qx * v[2]);
    let tz = 2.0 * (qx * v[1] - qy * v[0]);
    // result = v + qw * t + cross(q.xyz, t)
    [
        v[0] + qw * tx + (qy * tz - qz * ty),
        v[1] + qw * ty + (qz * tx - qx * tz),
        v[2] + qw * tz + (qx * ty - qy * tx),
    ]
}

/// Multiply two unit quaternions `a` and `b` (Hamilton product).
///
/// Format for both operands and the result: \[x, y, z, w\].
pub fn quat_mul(a: [f32; 4], b: [f32; 4]) -> [f32; 4] {
    let (ax, ay, az, aw) = (a[0], a[1], a[2], a[3]);
    let (bx, by, bz, bw) = (b[0], b[1], b[2], b[3]);
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}

/// Normalise a quaternion to unit length.
///
/// Returns the identity quaternion `[0,0,0,1]` if the norm is nearly zero.
pub fn quat_normalize(q: [f32; 4]) -> [f32; 4] {
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if norm < 1e-9 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm]
}

/// Integrate a quaternion orientation `q` by angular velocity `omega` (rad/s)
/// over time step `dt` (seconds).
///
/// Uses the first-order approximation `q' = normalize(q + 0.5 * dt * Ω⊗q)`.
pub fn integrate_orientation(q: [f32; 4], omega: [f32; 3], dt: f32) -> [f32; 4] {
    // omega as a pure quaternion [ox, oy, oz, 0]
    let omega_q = [omega[0], omega[1], omega[2], 0.0_f32];
    let dq = quat_mul(omega_q, q);
    let q_new = [
        q[0] + 0.5 * dt * dq[0],
        q[1] + 0.5 * dt * dq[1],
        q[2] + 0.5 * dt * dq[2],
        q[3] + 0.5 * dt * dq[3],
    ];
    quat_normalize(q_new)
}

// ── Core rigid body ───────────────────────────────────────────────────────────

/// A single rigid body stored in GPU-friendly packed arrays of `f32`.
///
/// Orientation is stored as a unit quaternion `[x, y, z, w]`.
/// The inverse inertia tensor is stored in row-major order as a 3×3 matrix.
#[derive(Debug, Clone)]
pub struct GpuRigidBody {
    /// World-space position (metres).
    pub position: [f32; 3],
    /// Linear velocity (m/s).
    pub velocity: [f32; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub orientation: [f32; 4],
    /// Angular velocity in world space (rad/s).
    pub angular_velocity: [f32; 3],
    /// Mass in kilograms.
    pub mass: f32,
    /// Inverse of the body-space inertia tensor, row-major (3×3).
    pub inv_inertia: [f32; 9],
}

impl GpuRigidBody {
    /// Create a new rigid body at rest with the given mass and diagonal inertia.
    ///
    /// `ixx`, `iyy`, `izz` are the principal moments of inertia.
    pub fn new(mass: f32, ixx: f32, iyy: f32, izz: f32) -> Self {
        let safe_inv = |v: f32| if v.abs() > 1e-12 { 1.0 / v } else { 0.0 };
        let inv_inertia = [
            safe_inv(ixx),
            0.0,
            0.0,
            0.0,
            safe_inv(iyy),
            0.0,
            0.0,
            0.0,
            safe_inv(izz),
        ];
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            mass,
            inv_inertia,
        }
    }

    /// Apply an inverse inertia tensor (row-major 3×3) to a vector.
    #[allow(dead_code)]
    fn apply_inv_inertia(&self, v: [f32; 3]) -> [f32; 3] {
        let i = &self.inv_inertia;
        [
            i[0] * v[0] + i[1] * v[1] + i[2] * v[2],
            i[3] * v[0] + i[4] * v[1] + i[5] * v[2],
            i[6] * v[0] + i[7] * v[1] + i[8] * v[2],
        ]
    }
}

// ── Batch integration ─────────────────────────────────────────────────────────

/// A batch of GPU rigid bodies supporting bulk integration.
#[derive(Debug, Clone, Default)]
pub struct GpuRigidBodyBatch {
    /// The rigid bodies in this batch.
    pub bodies: Vec<GpuRigidBody>,
}

impl GpuRigidBodyBatch {
    /// Create an empty batch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a body to the batch and return its index.
    pub fn add(&mut self, body: GpuRigidBody) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(body);
        idx
    }

    /// Integrate all bodies by `dt` seconds under a uniform `gravity` (m/s²).
    ///
    /// Uses explicit Euler integration for linear dynamics and a first-order
    /// quaternion update for rotational dynamics.
    pub fn integrate_all(&mut self, dt: f32, gravity: [f32; 3]) {
        for body in &mut self.bodies {
            // linear: v += g*dt,  p += v*dt
            body.velocity[0] += gravity[0] * dt;
            body.velocity[1] += gravity[1] * dt;
            body.velocity[2] += gravity[2] * dt;
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;
            // rotational
            let new_q = integrate_orientation(body.orientation, body.angular_velocity, dt);
            body.orientation = new_q;
        }
    }

    /// Apply a linear impulse `impulse` (N·s) at world-space `point` to body
    /// at `body_idx`.
    ///
    /// Updates both linear and angular velocity.
    pub fn apply_impulse(&mut self, body_idx: usize, impulse: [f32; 3], point: [f32; 3]) {
        let body = &mut self.bodies[body_idx];
        let inv_mass = if body.mass > 1e-12 {
            1.0 / body.mass
        } else {
            0.0
        };
        // linear impulse
        body.velocity[0] += impulse[0] * inv_mass;
        body.velocity[1] += impulse[1] * inv_mass;
        body.velocity[2] += impulse[2] * inv_mass;
        // torque arm: r = point - position
        let r = [
            point[0] - body.position[0],
            point[1] - body.position[1],
            point[2] - body.position[2],
        ];
        // angular impulse: omega += I⁻¹ * (r × impulse)
        let torque_impulse = cross3f(r, impulse);
        let delta_omega = apply_mat3(body.inv_inertia, torque_impulse);
        body.angular_velocity[0] += delta_omega[0];
        body.angular_velocity[1] += delta_omega[1];
        body.angular_velocity[2] += delta_omega[2];
    }
}

// ── Broadphase SAP ────────────────────────────────────────────────────────────

/// A candidate collision pair from the broadphase, with AABB information.
#[derive(Debug, Clone)]
pub struct BroadphasePairGpu {
    /// Index of the first body.
    pub body_a: usize,
    /// Index of the second body.
    pub body_b: usize,
    /// AABB centre of body A.
    pub aabb_a_center: [f32; 3],
    /// AABB half-extents of body A.
    pub aabb_a_half: [f32; 3],
    /// AABB centre of body B.
    pub aabb_b_center: [f32; 3],
    /// AABB half-extents of body B.
    pub aabb_b_half: [f32; 3],
}

/// Broadphase collision detection for GPU rigid bodies using a sweep-and-prune
/// (SAP) approach.
///
/// Each body is approximated by a sphere; overlapping sphere AABBs are reported
/// as candidate pairs.
#[derive(Debug, Clone, Default)]
pub struct GpuBroadphase {
    /// Bodies to test.
    pub bodies: Vec<GpuRigidBody>,
    /// Bounding sphere radii, one per body.
    pub radii: Vec<f32>,
}

impl GpuBroadphase {
    /// Create an empty broadphase structure.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a body with its bounding sphere radius.
    pub fn add_body(&mut self, body: GpuRigidBody, radius: f32) {
        self.bodies.push(body);
        self.radii.push(radius);
    }

    /// Run a sort-and-sweep on the X axis and return overlapping pairs.
    ///
    /// Only pairs whose AABBs overlap on all three axes are returned.
    pub fn compute_pairs_sap(&self) -> Vec<BroadphasePairGpu> {
        let n = self.bodies.len();
        let mut pairs = Vec::new();

        // Sort by AABB min-x
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            let ax = self.bodies[a].position[0] - self.radii[a];
            let bx = self.bodies[b].position[0] - self.radii[b];
            ax.partial_cmp(&bx).unwrap_or(std::cmp::Ordering::Equal)
        });

        for i in 0..order.len() {
            let ai = order[i];
            let a_max_x = self.bodies[ai].position[0] + self.radii[ai];
            for j in (i + 1)..order.len() {
                let bi = order[j];
                let b_min_x = self.bodies[bi].position[0] - self.radii[bi];
                if b_min_x > a_max_x {
                    break; // sorted: no further pair can overlap on X
                }
                // Check all three axes
                if self.aabb_overlap(ai, bi) {
                    let ra = self.radii[ai];
                    let rb = self.radii[bi];
                    pairs.push(BroadphasePairGpu {
                        body_a: ai,
                        body_b: bi,
                        aabb_a_center: self.bodies[ai].position,
                        aabb_a_half: [ra, ra, ra],
                        aabb_b_center: self.bodies[bi].position,
                        aabb_b_half: [rb, rb, rb],
                    });
                }
            }
        }
        pairs
    }

    /// Test whether two bodies' sphere AABBs overlap on all three axes.
    fn aabb_overlap(&self, a: usize, b: usize) -> bool {
        for k in 0..3 {
            let a_min = self.bodies[a].position[k] - self.radii[a];
            let a_max = self.bodies[a].position[k] + self.radii[a];
            let b_min = self.bodies[b].position[k] - self.radii[b];
            let b_max = self.bodies[b].position[k] + self.radii[b];
            if a_max < b_min || b_max < a_min {
                return false;
            }
        }
        true
    }
}

// ── Contact manifold ──────────────────────────────────────────────────────────

/// Contact manifold between two bodies produced by the narrowphase.
#[derive(Debug, Clone)]
pub struct ContactManifoldGpu {
    /// Index of the first body.
    pub body_a: usize,
    /// Index of the second body.
    pub body_b: usize,
    /// Contact point positions in world space.
    pub contact_points: Vec<[f32; 3]>,
    /// Outward contact normals (pointing from B to A), one per contact point.
    pub normals: Vec<[f32; 3]>,
    /// Penetration depths (positive means overlapping), one per contact point.
    pub penetrations: Vec<f32>,
}

impl ContactManifoldGpu {
    /// Create a new, empty contact manifold between two bodies.
    pub fn new(body_a: usize, body_b: usize) -> Self {
        Self {
            body_a,
            body_b,
            contact_points: Vec::new(),
            normals: Vec::new(),
            penetrations: Vec::new(),
        }
    }

    /// Add a single contact point to the manifold.
    pub fn add_contact(&mut self, point: [f32; 3], normal: [f32; 3], penetration: f32) {
        self.contact_points.push(point);
        self.normals.push(normal);
        self.penetrations.push(penetration);
    }

    /// Number of contact points in the manifold.
    pub fn contact_count(&self) -> usize {
        self.contact_points.len()
    }
}

// ── Sequential impulse solver ─────────────────────────────────────────────────

/// Iterative sequential-impulse constraint solver for rigid body contacts.
#[derive(Debug, Clone, Default)]
pub struct GpuConstraintSolver {
    /// Broadphase candidate pairs.
    pub pairs: Vec<BroadphasePairGpu>,
    /// Contact manifolds, aligned with `pairs`.
    pub manifolds: Vec<ContactManifoldGpu>,
}

impl GpuConstraintSolver {
    /// Create a new solver with no constraints.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a manifold (and its corresponding broadphase pair) to the solver.
    pub fn add_manifold(&mut self, pair: BroadphasePairGpu, manifold: ContactManifoldGpu) {
        self.pairs.push(pair);
        self.manifolds.push(manifold);
    }

    /// Solve all contact constraints via sequential impulses.
    ///
    /// Iterates `iterations` times over every contact point and applies a
    /// non-penetration impulse. A restitution coefficient of 0.3 is used.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_sequential_impulse(
        &self,
        bodies: &mut [GpuRigidBody],
        dt: f32,
        iterations: usize,
    ) {
        let _ = dt; // reserved for future bias/baumgarte
        let restitution = 0.3_f32;
        for _ in 0..iterations {
            for manifold in &self.manifolds {
                let a = manifold.body_a;
                let b = manifold.body_b;
                for c in 0..manifold.contact_count() {
                    let n = manifold.normals[c];
                    // relative velocity at contact
                    let va = bodies[a].velocity;
                    let vb = bodies[b].velocity;
                    let rv = [va[0] - vb[0], va[1] - vb[1], va[2] - vb[2]];
                    let vn = dot3f(rv, n);
                    if vn >= 0.0 {
                        continue; // separating
                    }
                    let inv_ma = if bodies[a].mass > 1e-12 {
                        1.0 / bodies[a].mass
                    } else {
                        0.0
                    };
                    let inv_mb = if bodies[b].mass > 1e-12 {
                        1.0 / bodies[b].mass
                    } else {
                        0.0
                    };
                    let j = -(1.0 + restitution) * vn / (inv_ma + inv_mb);
                    // apply
                    bodies[a].velocity[0] += j * inv_ma * n[0];
                    bodies[a].velocity[1] += j * inv_ma * n[1];
                    bodies[a].velocity[2] += j * inv_ma * n[2];
                    bodies[b].velocity[0] -= j * inv_mb * n[0];
                    bodies[b].velocity[1] -= j * inv_mb * n[1];
                    bodies[b].velocity[2] -= j * inv_mb * n[2];
                }
            }
        }
    }
}

// ── Internal math helpers ─────────────────────────────────────────────────────

/// Cross product of two `[f32; 3]` vectors.
fn cross3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two `[f32; 3]` vectors.
fn dot3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Apply a row-major 3×3 matrix to a vector.
fn apply_mat3(m: [f32; 9], v: [f32; 3]) -> [f32; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f32 = 1e-5;

    fn approx_eq3(a: [f32; 3], b: [f32; 3]) -> bool {
        (a[0] - b[0]).abs() < EPS && (a[1] - b[1]).abs() < EPS && (a[2] - b[2]).abs() < EPS
    }

    fn approx_eq4(a: [f32; 4], b: [f32; 4]) -> bool {
        (a[0] - b[0]).abs() < EPS
            && (a[1] - b[1]).abs() < EPS
            && (a[2] - b[2]).abs() < EPS
            && (a[3] - b[3]).abs() < EPS
    }

    fn quat_norm(q: [f32; 4]) -> f32 {
        (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt()
    }

    // ── quat_rotate ────────────────────────────────────────────────────────

    #[test]
    fn test_quat_rotate_identity() {
        let q = [0.0, 0.0, 0.0, 1.0_f32];
        let v = [1.0, 2.0, 3.0_f32];
        let r = quat_rotate(q, v);
        assert!(approx_eq3(r, v), "identity quat should not rotate: {r:?}");
    }

    #[test]
    fn test_quat_rotate_180_about_z() {
        // 180° about Z: q = [0, 0, 1, 0]
        let q = [0.0, 0.0, 1.0_f32, 0.0];
        let v = [1.0, 0.0, 0.0_f32];
        let r = quat_rotate(q, v);
        assert!(approx_eq3(r, [-1.0, 0.0, 0.0]), "180 Z rotate: {r:?}");
    }

    #[test]
    fn test_quat_rotate_90_about_y() {
        // 90° about Y: q = [0, sin45, 0, cos45]
        let half = std::f32::consts::FRAC_PI_4;
        let q = [0.0, half.sin(), 0.0, half.cos()];
        let v = [1.0, 0.0, 0.0_f32];
        let r = quat_rotate(q, v);
        // Should become ~[0, 0, -1]
        assert!((r[0]).abs() < EPS, "x should be ~0: {}", r[0]);
        assert!((r[1]).abs() < EPS, "y should be ~0: {}", r[1]);
        assert!((r[2] + 1.0).abs() < EPS, "z should be ~-1: {}", r[2]);
    }

    #[test]
    fn test_quat_rotate_preserves_length() {
        let half = std::f32::consts::FRAC_PI_6;
        let q = quat_normalize([half.sin(), 0.0, 0.0, half.cos()]);
        let v = [3.0, 4.0, 0.0_f32];
        let r = quat_rotate(q, v);
        let len_v = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let len_r = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        assert!((len_v - len_r).abs() < EPS, "rotation must preserve length");
    }

    // ── quat_mul ───────────────────────────────────────────────────────────

    #[test]
    fn test_quat_mul_identity() {
        let id = [0.0, 0.0, 0.0, 1.0_f32];
        let q = [0.1, 0.2, 0.3_f32, 0.9];
        let q = quat_normalize(q);
        assert!(approx_eq4(quat_mul(id, q), q));
        assert!(approx_eq4(quat_mul(q, id), q));
    }

    #[test]
    fn test_quat_mul_unit_norm() {
        let a = quat_normalize([1.0, 0.0, 0.0_f32, 1.0]);
        let b = quat_normalize([0.0, 1.0, 0.0_f32, 1.0]);
        let c = quat_mul(a, b);
        assert!(
            (quat_norm(c) - 1.0).abs() < EPS,
            "product must be unit: {}",
            quat_norm(c)
        );
    }

    #[test]
    fn test_quat_mul_double_rotation() {
        // Two 90° rotations about Z => 180° about Z
        let half = std::f32::consts::FRAC_PI_4;
        let q90 = [0.0, 0.0, half.sin(), half.cos()];
        let q180 = quat_mul(q90, q90);
        let v = [1.0, 0.0, 0.0_f32];
        let r = quat_rotate(q180, v);
        assert!((r[0] + 1.0).abs() < EPS, "should be [-1,0,0]: {r:?}");
    }

    // ── integrate_orientation ──────────────────────────────────────────────

    #[test]
    fn test_integrate_orientation_no_rotation() {
        let q = [0.0, 0.0, 0.0, 1.0_f32];
        let q2 = integrate_orientation(q, [0.0; 3], 0.01);
        assert!((quat_norm(q2) - 1.0).abs() < EPS);
        assert!(approx_eq4(q2, q));
    }

    #[test]
    fn test_integrate_orientation_stays_unit() {
        let q = [0.0, 0.0, 0.0, 1.0_f32];
        let omega = [0.0, 0.0, 1.0_f32]; // 1 rad/s about Z
        let mut q_cur = q;
        for _ in 0..100 {
            q_cur = integrate_orientation(q_cur, omega, 0.01);
        }
        assert!((quat_norm(q_cur) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_integrate_orientation_direction() {
        // Small rotation about X: after many steps orientation should shift
        let q = [0.0, 0.0, 0.0, 1.0_f32];
        let omega = [1.0, 0.0, 0.0_f32];
        let q2 = integrate_orientation(q, omega, 0.1);
        // x component should now be non-zero
        assert!(
            q2[0].abs() > 1e-4,
            "qx should be > 0 after rotation: {}",
            q2[0]
        );
    }

    // ── quat_normalize ─────────────────────────────────────────────────────

    #[test]
    fn test_quat_normalize_unit() {
        let q = [1.0_f32, 0.0, 0.0, 0.0];
        assert!(approx_eq4(quat_normalize(q), q));
    }

    #[test]
    fn test_quat_normalize_zero_returns_identity() {
        let q = quat_normalize([0.0; 4]);
        assert!(approx_eq4(q, [0.0, 0.0, 0.0, 1.0]));
    }

    #[test]
    fn test_quat_normalize_scales() {
        let q = [2.0_f32, 0.0, 0.0, 0.0];
        let n = quat_normalize(q);
        assert!((n[0] - 1.0).abs() < EPS);
    }

    // ── GpuRigidBody ──────────────────────────────────────────────────────

    #[test]
    fn test_rigid_body_new_defaults() {
        let b = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        assert_eq!(b.position, [0.0; 3]);
        assert_eq!(b.velocity, [0.0; 3]);
        assert_eq!(b.orientation, [0.0, 0.0, 0.0, 1.0]);
        assert_eq!(b.angular_velocity, [0.0; 3]);
        assert!((b.mass - 1.0).abs() < EPS);
    }

    #[test]
    fn test_rigid_body_inv_inertia_diagonal() {
        let b = GpuRigidBody::new(1.0, 2.0, 4.0, 8.0);
        assert!((b.inv_inertia[0] - 0.5).abs() < EPS, "ixx");
        assert!((b.inv_inertia[4] - 0.25).abs() < EPS, "iyy");
        assert!((b.inv_inertia[8] - 0.125).abs() < EPS, "izz");
    }

    // ── GpuRigidBodyBatch::integrate_all ──────────────────────────────────

    #[test]
    fn test_batch_gravity_integration() {
        let mut batch = GpuRigidBodyBatch::new();
        batch.add(GpuRigidBody::new(1.0, 1.0, 1.0, 1.0));
        batch.integrate_all(1.0, [0.0, -9.81, 0.0]);
        assert!((batch.bodies[0].velocity[1] + 9.81).abs() < EPS);
        assert!((batch.bodies[0].position[1] + 9.81).abs() < EPS);
    }

    #[test]
    fn test_batch_no_gravity_no_motion() {
        let mut batch = GpuRigidBodyBatch::new();
        batch.add(GpuRigidBody::new(1.0, 1.0, 1.0, 1.0));
        batch.integrate_all(1.0, [0.0; 3]);
        assert_eq!(batch.bodies[0].position, [0.0; 3]);
        assert_eq!(batch.bodies[0].velocity, [0.0; 3]);
    }

    #[test]
    fn test_batch_multiple_bodies() {
        let mut batch = GpuRigidBodyBatch::new();
        for _ in 0..5 {
            batch.add(GpuRigidBody::new(1.0, 1.0, 1.0, 1.0));
        }
        batch.integrate_all(0.5, [0.0, -9.81, 0.0]);
        for b in &batch.bodies {
            assert!((b.velocity[1] + 9.81 * 0.5).abs() < EPS);
        }
    }

    #[test]
    fn test_batch_add_returns_index() {
        let mut batch = GpuRigidBodyBatch::new();
        let i0 = batch.add(GpuRigidBody::new(1.0, 1.0, 1.0, 1.0));
        let i1 = batch.add(GpuRigidBody::new(2.0, 1.0, 1.0, 1.0));
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
    }

    // ── GpuRigidBodyBatch::apply_impulse ──────────────────────────────────

    #[test]
    fn test_apply_impulse_linear() {
        let mut batch = GpuRigidBodyBatch::new();
        batch.add(GpuRigidBody::new(2.0, 1.0, 1.0, 1.0));
        // impulse [2,0,0] on mass 2 => dv = [1,0,0]
        batch.apply_impulse(0, [2.0, 0.0, 0.0], [0.0; 3]);
        assert!((batch.bodies[0].velocity[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_apply_impulse_angular() {
        let mut batch = GpuRigidBodyBatch::new();
        let mut b = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b.position = [0.0; 3];
        batch.add(b);
        // Apply Z-axis impulse at offset on X axis → angular velocity change
        batch.apply_impulse(0, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(batch.bodies[0].angular_velocity[2].abs() > 1e-5);
    }

    #[test]
    fn test_apply_impulse_zero_mass() {
        let mut batch = GpuRigidBodyBatch::new();
        batch.add(GpuRigidBody::new(0.0, 1.0, 1.0, 1.0));
        batch.apply_impulse(0, [100.0, 0.0, 0.0], [0.0; 3]);
        assert_eq!(batch.bodies[0].velocity, [0.0; 3]);
    }

    // ── GpuBroadphase ─────────────────────────────────────────────────────

    #[test]
    fn test_broadphase_no_bodies() {
        let bp = GpuBroadphase::new();
        assert!(bp.compute_pairs_sap().is_empty());
    }

    #[test]
    fn test_broadphase_two_overlapping() {
        let mut bp = GpuBroadphase::new();
        let mut b1 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b1.position = [0.0; 3];
        let mut b2 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b2.position = [0.5, 0.0, 0.0];
        bp.add_body(b1, 1.0);
        bp.add_body(b2, 1.0);
        let pairs = bp.compute_pairs_sap();
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn test_broadphase_two_separated() {
        let mut bp = GpuBroadphase::new();
        let mut b1 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b1.position = [0.0; 3];
        let mut b2 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b2.position = [100.0, 0.0, 0.0];
        bp.add_body(b1, 0.5);
        bp.add_body(b2, 0.5);
        let pairs = bp.compute_pairs_sap();
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_broadphase_three_bodies_two_pairs() {
        let mut bp = GpuBroadphase::new();
        for x in [0.0_f32, 1.0, 2.0] {
            let mut b = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
            b.position = [x, 0.0, 0.0];
            bp.add_body(b, 0.8);
        }
        let pairs = bp.compute_pairs_sap();
        // 0-1 and 1-2 overlap, 0-2 may or may not
        assert!(pairs.len() >= 2, "expected >= 2 pairs, got {}", pairs.len());
    }

    #[test]
    fn test_broadphase_pair_indices_valid() {
        let mut bp = GpuBroadphase::new();
        let mut b1 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b1.position = [0.0; 3];
        let mut b2 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b2.position = [0.3, 0.0, 0.0];
        bp.add_body(b1, 0.5);
        bp.add_body(b2, 0.5);
        let pairs = bp.compute_pairs_sap();
        assert_eq!(pairs.len(), 1);
        let p = &pairs[0];
        assert!(p.body_a < 2 && p.body_b < 2 && p.body_a != p.body_b);
    }

    // ── ContactManifoldGpu ────────────────────────────────────────────────

    #[test]
    fn test_manifold_empty() {
        let m = ContactManifoldGpu::new(0, 1);
        assert_eq!(m.contact_count(), 0);
    }

    #[test]
    fn test_manifold_add_contact() {
        let mut m = ContactManifoldGpu::new(0, 1);
        m.add_contact([0.0; 3], [0.0, 1.0, 0.0], 0.01);
        assert_eq!(m.contact_count(), 1);
        assert!((m.penetrations[0] - 0.01).abs() < EPS);
    }

    #[test]
    fn test_manifold_multiple_contacts() {
        let mut m = ContactManifoldGpu::new(0, 1);
        for i in 0..4 {
            m.add_contact([i as f32, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01 * i as f32);
        }
        assert_eq!(m.contact_count(), 4);
    }

    // ── GpuConstraintSolver ───────────────────────────────────────────────

    #[test]
    fn test_solver_no_manifolds() {
        let solver = GpuConstraintSolver::new();
        let mut bodies = vec![GpuRigidBody::new(1.0, 1.0, 1.0, 1.0)];
        solver.solve_sequential_impulse(&mut bodies, 0.01, 10);
        assert_eq!(bodies[0].velocity, [0.0; 3]);
    }

    #[test]
    fn test_solver_separates_colliding_bodies() {
        let mut b1 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b1.velocity = [1.0, 0.0, 0.0];
        let mut b2 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b2.velocity = [-1.0, 0.0, 0.0];
        b2.position = [0.5, 0.0, 0.0];

        let mut solver = GpuConstraintSolver::new();
        let pair = BroadphasePairGpu {
            body_a: 0,
            body_b: 1,
            aabb_a_center: [0.0; 3],
            aabb_a_half: [0.5; 3],
            aabb_b_center: [0.5; 3],
            aabb_b_half: [0.5; 3],
        };
        let mut manifold = ContactManifoldGpu::new(0, 1);
        manifold.add_contact([0.25, 0.0, 0.0], [1.0, 0.0, 0.0], 0.01);
        solver.add_manifold(pair, manifold);

        let mut bodies = vec![b1, b2];
        solver.solve_sequential_impulse(&mut bodies, 0.01, 10);

        // After impulse, relative velocity along normal should be >= 0
        let rv_x = bodies[0].velocity[0] - bodies[1].velocity[0];
        assert!(
            rv_x >= -EPS,
            "bodies should not penetrate further: rv_x={rv_x}"
        );
    }

    #[test]
    fn test_solver_static_body_kinematic() {
        // body_b has zero mass => should not move
        let mut b1 = GpuRigidBody::new(1.0, 1.0, 1.0, 1.0);
        b1.velocity = [0.0, -1.0, 0.0];
        let b2 = GpuRigidBody::new(0.0, 1.0, 1.0, 1.0); // static

        let mut solver = GpuConstraintSolver::new();
        let pair = BroadphasePairGpu {
            body_a: 0,
            body_b: 1,
            aabb_a_center: [0.0; 3],
            aabb_a_half: [0.5; 3],
            aabb_b_center: [0.0, -0.9, 0.0],
            aabb_b_half: [0.5; 3],
        };
        let mut manifold = ContactManifoldGpu::new(0, 1);
        manifold.add_contact([0.0, -0.5, 0.0], [0.0, 1.0, 0.0], 0.1);
        solver.add_manifold(pair, manifold);

        let mut bodies = vec![b1, b2];
        solver.solve_sequential_impulse(&mut bodies, 0.01, 10);

        assert_eq!(bodies[1].velocity, [0.0; 3], "static body must not move");
        assert!(
            bodies[0].velocity[1] >= -EPS,
            "dynamic body should bounce up"
        );
    }

    // ── Internal helpers ──────────────────────────────────────────────────

    #[test]
    fn test_cross3f() {
        let i = [1.0_f32, 0.0, 0.0];
        let j = [0.0_f32, 1.0, 0.0];
        let k = cross3f(i, j);
        assert!(approx_eq3(k, [0.0, 0.0, 1.0]));
    }

    #[test]
    fn test_dot3f() {
        assert!((dot3f([1.0, 2.0, 3.0_f32], [4.0, 5.0, 6.0]) - 32.0).abs() < EPS);
    }

    #[test]
    fn test_apply_mat3_identity() {
        let id = [1.0_f32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let v = [3.0_f32, 4.0, 5.0];
        assert!(approx_eq3(apply_mat3(id, v), v));
    }

    #[test]
    fn test_apply_mat3_scale() {
        let m = [2.0_f32, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0];
        let v = [1.0_f32, 1.0, 1.0];
        assert!(approx_eq3(apply_mat3(m, v), [2.0, 3.0, 4.0]));
    }
}
