// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint types for articulated-body models.
//!
//! Each joint implements the [`Joint`] trait providing:
//! - the joint pose `X_J(q)` of the child joint frame relative to the parent
//!   joint frame (see [`SpatialTransform`] for the convention)
//! - the configuration/velocity bookkeeping (`nq`, `neutral`, `integrate`,
//!   `difference`) for joints whose configuration is not Euclidean
//! - the motion subspace `S` mapping q̇ → spatial velocity
//! - optional Coriolis term `c_J`

use crate::spatial::{SpatialTransform, SpatialVec};

// ─── Joint trait ─────────────────────────────────────────────────────────────

/// A kinematic joint connecting a parent body to a child body.
///
/// Implementors must be `Send + Sync` for multi-threaded use.
pub trait Joint: Send + Sync {
    /// Number of degrees of freedom `nv` (0 for fixed, 1 for revolute/prismatic,
    /// 3 for spherical, 6 for free-floating, …).
    ///
    /// This is the length of the joint's velocity, acceleration and
    /// generalised-force slices.
    fn dof(&self) -> usize;

    /// Number of configuration coordinates `nq` (length of the joint's slice of
    /// `q`).
    ///
    /// Equals [`Joint::dof`] for joints whose configuration space is Euclidean.
    /// Joints that carry a unit quaternion override it (spherical: 4,
    /// free-floating: 7).
    fn nq(&self) -> usize {
        self.dof()
    }

    /// Neutral (reference) configuration of length [`Joint::nq`].
    ///
    /// Zero for Euclidean joints, the identity quaternion for rotational
    /// joints parameterised by a quaternion.
    fn neutral(&self) -> Vec<f64> {
        vec![0.0; self.nq()]
    }

    /// Integrate a constant joint velocity over `dt`: `q_out = q ⊕ v·dt`.
    ///
    /// `q` and `q_out` have length [`Joint::nq`]; `v` has length
    /// [`Joint::dof`]. The default is the Euclidean update `q + v·dt`; joints
    /// with a quaternion use the exponential map.
    fn integrate(&self, q: &[f64], v: &[f64], dt: f64, q_out: &mut [f64]) {
        for ((out, &qk), &vk) in q_out.iter_mut().zip(q.iter()).zip(v.iter()) {
            *out = qk + vk * dt;
        }
    }

    /// Tangent-space difference `v_out = q1 ⊖ q0`, the inverse of
    /// [`Joint::integrate`]: `integrate(q0, v_out, 1.0) == q1`.
    ///
    /// `q0`, `q1` have length [`Joint::nq`]; `v_out` has length [`Joint::dof`].
    fn difference(&self, q0: &[f64], q1: &[f64], v_out: &mut [f64]) {
        for ((out, &a), &b) in v_out.iter_mut().zip(q0.iter()).zip(q1.iter()) {
            *out = b - a;
        }
    }

    /// Joint pose `X_J(q)`: the child joint frame relative to the parent joint
    /// frame (see [`SpatialTransform`] for the convention).
    ///
    /// `q` must have length [`Joint::nq`].
    fn transform(&self, q: &[f64]) -> SpatialTransform;

    /// Motion subspace matrix columns `S_i`.
    ///
    /// The spatial velocity due to the joint is `v_J = Σ_i S_i · q̇_i`.
    /// Returns a `Vec` of length `self.dof()`, each entry a spatial vector.
    fn motion_subspace(&self) -> Vec<SpatialVec>;

    /// Coriolis acceleration `c_J = Ṡ · q̇` (bias term in acceleration propagation).
    ///
    /// For constant-axis joints (revolute, prismatic, fixed) this is always zero.
    /// For variable-geometry joints it may depend on `q` and `q_dot`.
    fn coriolis(&self, _q: &[f64], _q_dot: &[f64]) -> SpatialVec {
        SpatialVec::ZERO
    }

    /// Motion subspace columns evaluated at configuration `q`.
    ///
    /// For simple joints (revolute, prismatic, fixed, free-floating), the motion
    /// subspace is configuration-independent, so this delegates to `motion_subspace()`.
    /// Override for joints where `S` depends on `q` (e.g., universal/cardan joint).
    fn motion_subspace_at(&self, _q: &[f64]) -> Vec<SpatialVec> {
        self.motion_subspace()
    }
}

// ─── Helper: axis normalisation ───────────────────────────────────────────────

fn normalize_axis(axis: [f64; 3]) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if n < 1e-14 {
        [0.0, 0.0, 1.0] // fallback to z-axis
    } else {
        [axis[0] / n, axis[1] / n, axis[2] / n]
    }
}

// ─── RevoluteJoint ───────────────────────────────────────────────────────────

/// 1-DOF revolute joint: rotation about a fixed axis.
///
/// `q[0]` is the rotation angle \[rad\]. The motion subspace is `[axis; 0]`
/// (pure rotation, no translation).
pub struct RevoluteJoint {
    /// Rotation axis in the joint frame (unit vector, normalised on construction).
    pub axis: [f64; 3],
}

impl RevoluteJoint {
    /// Create a revolute joint with the given axis (will be normalised).
    pub fn new(axis: [f64; 3]) -> Self {
        Self {
            axis: normalize_axis(axis),
        }
    }
}

impl Joint for RevoluteJoint {
    fn dof(&self) -> usize {
        1
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let angle = q[0];
        SpatialTransform::from_axis_angle(self.axis, angle)
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        // S = [axis; 0] — pure rotation
        vec![SpatialVec::new(self.axis, [0.0; 3])]
    }
}

// ─── PrismaticJoint ──────────────────────────────────────────────────────────

/// 1-DOF prismatic joint: translation along a fixed axis.
///
/// `q[0]` is the displacement \[m\]. The motion subspace is `[0; axis]`
/// (pure translation, no rotation).
pub struct PrismaticJoint {
    /// Translation axis in the joint frame (unit vector, normalised on construction).
    pub axis: [f64; 3],
}

impl PrismaticJoint {
    /// Create a prismatic joint with the given axis (will be normalised).
    pub fn new(axis: [f64; 3]) -> Self {
        Self {
            axis: normalize_axis(axis),
        }
    }
}

impl Joint for PrismaticJoint {
    fn dof(&self) -> usize {
        1
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let d = q[0];
        let trans = [self.axis[0] * d, self.axis[1] * d, self.axis[2] * d];
        SpatialTransform::from_translation(trans)
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        // S = [0; axis] — pure translation
        vec![SpatialVec::new([0.0; 3], self.axis)]
    }
}

// ─── FixedJoint ──────────────────────────────────────────────────────────────

/// 0-DOF fixed joint: rigid attachment, no relative motion.
pub struct FixedJoint;

impl Joint for FixedJoint {
    fn dof(&self) -> usize {
        0
    }

    fn transform(&self, _q: &[f64]) -> SpatialTransform {
        SpatialTransform::IDENTITY
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        vec![] // no columns
    }
}

// ─── FreeFloatingJoint ───────────────────────────────────────────────────────

/// 6-DOF free-floating joint for floating-base systems.
///
/// **Configuration** (`nq = 7`):
/// - `q[0..3]` = position `[x, y, z]` of the child origin in the parent frame,
/// - `q[3..7]` = unit quaternion `[w, x, y, z]` giving the child orientation
///   relative to the parent (normalised on use).
///
/// **Velocity** (`nv = 6`) is the body-frame twist of the child:
/// - `q_dot[0..3]` = linear velocity of the child origin, in **child** coordinates,
/// - `q_dot[3..6]` = angular velocity, in **child** coordinates.
///
/// The motion subspace `S` is the 6×6 identity (reordered to `[v; ω]`), so the
/// joint velocity is directly the spatial velocity of the child in its own
/// frame. [`Joint::integrate`] applies the SE(3) exponential of the body
/// twist; [`Joint::difference`] is its exact inverse (SE(3) logarithm).
pub struct FreeFloatingJoint;

impl Joint for FreeFloatingJoint {
    fn dof(&self) -> usize {
        6
    }

    fn nq(&self) -> usize {
        7
    }

    fn neutral(&self) -> Vec<f64> {
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    }

    fn integrate(&self, q: &[f64], v: &[f64], dt: f64, q_out: &mut [f64]) {
        let p = [q[0], q[1], q[2]];
        let quat = quat_normalize([q[3], q[4], q[5], q[6]]);
        let rot = quat_to_rotmat_arr(quat);
        let phi = [v[3] * dt, v[4] * dt, v[5] * dt];
        let rho = [v[0] * dt, v[1] * dt, v[2] * dt];
        // Constant body twist: T(dt) = T · exp([ρ; φ]^) with
        // exp = (exp(φ×), V(φ)·ρ).
        let dp_body = so3_left_jacobian_apply(phi, rho);
        let dp = mat3_vec3(rot, dp_body);
        let quat_new = quat_normalize(quat_mul(quat, quat_exp(phi)));
        q_out[0] = p[0] + dp[0];
        q_out[1] = p[1] + dp[1];
        q_out[2] = p[2] + dp[2];
        q_out[3..7].copy_from_slice(&quat_new);
    }

    fn difference(&self, q0: &[f64], q1: &[f64], v_out: &mut [f64]) {
        let quat0 = quat_normalize([q0[3], q0[4], q0[5], q0[6]]);
        let quat1 = quat_normalize([q1[3], q1[4], q1[5], q1[6]]);
        let rot0 = quat_to_rotmat_arr(quat0);
        // Relative pose T0⁻¹·T1 = (R0ᵀR1, R0ᵀ(p1 − p0)).
        let dp_world = [q1[0] - q0[0], q1[1] - q0[1], q1[2] - q0[2]];
        let dp_body = mat3_t_vec3(rot0, dp_world);
        let phi = quat_log(quat_mul(quat_conj(quat0), quat1));
        let rho = so3_left_jacobian_inverse_apply(phi, dp_body);
        v_out[0] = rho[0];
        v_out[1] = rho[1];
        v_out[2] = rho[2];
        v_out[3] = phi[0];
        v_out[4] = phi[1];
        v_out[5] = phi[2];
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let trans = [q[0], q[1], q[2]];
        let rot = quat_to_rotmat_arr(quat_normalize([q[3], q[4], q[5], q[6]]));
        SpatialTransform::from_rotation_translation(rot, trans)
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        // 6 columns of 6×6 identity, ordered: [v_lin(3), ω(3)]
        // Convention: first 3 DOF are linear velocity, last 3 are angular
        vec![
            SpatialVec::new([0.0; 3], [1.0, 0.0, 0.0]), // vx
            SpatialVec::new([0.0; 3], [0.0, 1.0, 0.0]), // vy
            SpatialVec::new([0.0; 3], [0.0, 0.0, 1.0]), // vz
            SpatialVec::new([1.0, 0.0, 0.0], [0.0; 3]), // ωx
            SpatialVec::new([0.0, 1.0, 0.0], [0.0; 3]), // ωy
            SpatialVec::new([0.0, 0.0, 1.0], [0.0; 3]), // ωz
        ]
    }
}

// ─── UniversalJoint ──────────────────────────────────────────────────────────

/// 2-DOF universal (cardan/Hooke) joint: two perpendicular rotation axes.
///
/// Parameterisation:
/// - `q[0]` = rotation angle about `axis1` (in the parent joint frame).
/// - `q[1]` = rotation angle about `axis2` (in the intermediate frame after
///   the first rotation). For a standard cardan joint, `axis1 ⊥ axis2`.
///
/// The joint pose (child relative to parent) is
/// `X_J = R(axis1, q0) · R(axis2, q1)`: rotate about `axis1` first, then about
/// `axis2` expressed in the intermediate (once-rotated) frame.
///
/// **Motion subspace:** `S(q)` depends on `q[1]` — the first column (axis1's
/// contribution in the child frame) is `R(axis2, −q1)·axis1`. This is why the
/// [`Joint::motion_subspace_at`] override is needed. The q-independent fallback
/// `motion_subspace()` returns the columns at `q = [0, 0]`.
///
/// **Bias acceleration:** because `S` varies with `q1`, the joint contributes
/// a velocity-product term `c_J = Ṡ·q̇ = −q̇0·q̇1·(axis2 × S0(q))` (see
/// [`Joint::coriolis`]).
pub struct UniversalJoint {
    /// First rotation axis (in the parent joint frame). Normalised on construction.
    pub axis1: [f64; 3],
    /// Second rotation axis (in the intermediate frame after the first rotation).
    /// Normalised on construction.
    pub axis2: [f64; 3],
}

impl UniversalJoint {
    /// Create a universal joint with the given axes (both are normalised).
    ///
    /// For a standard cardan joint, `axis1` and `axis2` should be perpendicular.
    pub fn new(axis1: [f64; 3], axis2: [f64; 3]) -> Self {
        Self {
            axis1: normalize_axis(axis1),
            axis2: normalize_axis(axis2),
        }
    }
}

impl Joint for UniversalJoint {
    fn dof(&self) -> usize {
        2
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let r1 = SpatialTransform::from_axis_angle(self.axis1, q[0]);
        let r2 = SpatialTransform::from_axis_angle(self.axis2, q[1]);
        // Pose chain parent → intermediate (about axis1) → child (about axis2
        // in the intermediate frame): R_J = R1 · R2.
        r1.compose(&r2)
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        // At q = [0, 0]: columns are just axis1 and axis2 in the (still-aligned) child frame.
        vec![
            SpatialVec::new(self.axis1, [0.0; 3]),
            SpatialVec::new(self.axis2, [0.0; 3]),
        ]
    }

    fn motion_subspace_at(&self, q: &[f64]) -> Vec<SpatialVec> {
        // In the child frame (after both rotations), axis1's column must be
        // rotated back through the second rotation (axis2, angle q[1]):
        //   axis1_child = R(axis2, -q[1]) · axis1
        // axis2 is constant in the child frame (it's the second rotation axis).
        let q1 = if q.len() > 1 { q[1] } else { 0.0 };
        let axis1_child = rotate_vec3_axis_angle(self.axis2, -q1, self.axis1);
        vec![
            SpatialVec::new(axis1_child, [0.0; 3]),
            SpatialVec::new(self.axis2, [0.0; 3]),
        ]
    }

    fn coriolis(&self, q: &[f64], q_dot: &[f64]) -> SpatialVec {
        // S0(q) = R(axis2, −q1)·axis1 ⇒ ∂S0/∂q1 = −axis2 × S0(q); S1 is constant.
        // c_J = Ṡ·q̇ = (∂S0/∂q1·q̇1)·q̇0 = −q̇0·q̇1·(axis2 × S0(q)).
        let q1 = if q.len() > 1 { q[1] } else { 0.0 };
        let qd0 = q_dot.first().copied().unwrap_or(0.0);
        let qd1 = q_dot.get(1).copied().unwrap_or(0.0);
        let axis1_child = rotate_vec3_axis_angle(self.axis2, -q1, self.axis1);
        let a2_x_s0 = [
            self.axis2[1] * axis1_child[2] - self.axis2[2] * axis1_child[1],
            self.axis2[2] * axis1_child[0] - self.axis2[0] * axis1_child[2],
            self.axis2[0] * axis1_child[1] - self.axis2[1] * axis1_child[0],
        ];
        let k = -qd0 * qd1;
        SpatialVec::new([k * a2_x_s0[0], k * a2_x_s0[1], k * a2_x_s0[2]], [0.0; 3])
    }
}

// ─── SphericalJoint ──────────────────────────────────────────────────────────

/// 3-DOF spherical (ball-and-socket) joint for shoulder/hip articulations.
///
/// **Configuration** (`nq = 4`): a unit quaternion `[w, x, y, z]` giving the
/// child orientation relative to the parent (normalised on use). The model's
/// `q` vector carries all four components ([`crate::model::ArticulatedModel::q_slice`]
/// returns them), so every algorithm sees the true orientation.
///
/// **Velocity** (`nv = 3`): the body-frame angular velocity, matching the
/// constant motion subspace `[I₃; 0₃]`. [`Joint::integrate`] applies the
/// quaternion exponential `q ⊗ exp(ω·dt/2)`; [`Joint::difference`] returns the
/// body-frame rotation vector `log(q0⁻¹ ⊗ q1)`.
pub struct SphericalJoint;

impl Joint for SphericalJoint {
    fn dof(&self) -> usize {
        3
    }

    fn nq(&self) -> usize {
        4
    }

    fn neutral(&self) -> Vec<f64> {
        vec![1.0, 0.0, 0.0, 0.0]
    }

    fn integrate(&self, q: &[f64], v: &[f64], dt: f64, q_out: &mut [f64]) {
        let quat = quat_normalize([q[0], q[1], q[2], q[3]]);
        let phi = [v[0] * dt, v[1] * dt, v[2] * dt];
        let quat_new = quat_normalize(quat_mul(quat, quat_exp(phi)));
        q_out[..4].copy_from_slice(&quat_new);
    }

    fn difference(&self, q0: &[f64], q1: &[f64], v_out: &mut [f64]) {
        let quat0 = quat_normalize([q0[0], q0[1], q0[2], q0[3]]);
        let quat1 = quat_normalize([q1[0], q1[1], q1[2], q1[3]]);
        let phi = quat_log(quat_mul(quat_conj(quat0), quat1));
        v_out[..3].copy_from_slice(&phi);
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let rot = quat_to_rotmat_arr(quat_normalize([q[0], q[1], q[2], q[3]]));
        SpatialTransform::from_rotation_translation(rot, [0.0; 3])
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        // [I₃; 0₃]: three angular-only columns in the child frame.
        // Constant because q_dot is body-frame ω, not dq/dt.
        vec![
            SpatialVec::new([1.0, 0.0, 0.0], [0.0; 3]),
            SpatialVec::new([0.0, 1.0, 0.0], [0.0; 3]),
            SpatialVec::new([0.0, 0.0, 1.0], [0.0; 3]),
        ]
    }
}

// ─── HelicalJoint ────────────────────────────────────────────────────────────

/// 1-DOF helical (screw) joint coupling rotation and translation.
///
/// `q[0] = θ` gives both a rotation of `θ` rad about `axis` AND a translation
/// of `pitch_per_rad · θ` metres along `axis`. This models lead-screw actuators.
///
/// The motion subspace `S = [axis; pitch_per_rad · axis]` is constant.
pub struct HelicalJoint {
    /// Rotation/translation axis in the joint frame (unit vector).
    pub axis: [f64; 3],
    /// Axial translation per radian of rotation [m/rad].
    ///
    /// For a lead-screw with p metres per full turn, set `p / (2π)`.
    pub pitch_per_rad: f64,
}

impl HelicalJoint {
    /// Create a helical joint (axis is normalised automatically).
    pub fn new(axis: [f64; 3], pitch_per_rad: f64) -> Self {
        Self {
            axis: normalize_axis(axis),
            pitch_per_rad,
        }
    }

    /// Convenience constructor: specify lead-screw pitch in metres per full turn.
    pub fn from_turn_pitch(axis: [f64; 3], meters_per_turn: f64) -> Self {
        Self::new(axis, meters_per_turn / (2.0 * std::f64::consts::PI))
    }
}

impl Joint for HelicalJoint {
    fn dof(&self) -> usize {
        1
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let theta = q[0];
        let rot = SpatialTransform::from_axis_angle(self.axis, theta);
        let trans = [
            self.axis[0] * theta * self.pitch_per_rad,
            self.axis[1] * theta * self.pitch_per_rad,
            self.axis[2] * theta * self.pitch_per_rad,
        ];
        SpatialTransform::from_rotation_translation(rot.rot, trans)
    }

    fn motion_subspace(&self) -> Vec<SpatialVec> {
        let lin = [
            self.axis[0] * self.pitch_per_rad,
            self.axis[1] * self.pitch_per_rad,
            self.axis[2] * self.pitch_per_rad,
        ];
        vec![SpatialVec::new(self.axis, lin)]
    }
}

// ─── Helper: quaternion / SO(3) / SE(3) utilities ───────────────────────────

/// Unit quaternion `[w, x, y, z]`.
type Quat = [f64; 4];

/// Normalise a quaternion; a (numerically) zero quaternion maps to identity.
fn quat_normalize(q: Quat) -> Quat {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n < 1e-300 {
        [1.0, 0.0, 0.0, 0.0]
    } else {
        [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
    }
}

/// Hamilton product `a ⊗ b`.
fn quat_mul(a: Quat, b: Quat) -> Quat {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// Conjugate (inverse of a unit quaternion).
fn quat_conj(q: Quat) -> Quat {
    [q[0], -q[1], -q[2], -q[3]]
}

/// Quaternion exponential of a rotation vector `φ` (rotation by `|φ|` about
/// `φ/|φ|`): `[cos(θ/2), sin(θ/2)·φ/θ]`.
fn quat_exp(phi: [f64; 3]) -> Quat {
    let theta = (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]).sqrt();
    let half = 0.5 * theta;
    // sin(θ/2)/θ, with its Taylor series near zero.
    let k = if theta < 1e-6 {
        0.5 - theta * theta / 48.0
    } else {
        half.sin() / theta
    };
    [half.cos(), k * phi[0], k * phi[1], k * phi[2]]
}

/// Quaternion logarithm: the rotation vector `φ` with `quat_exp(φ) = ±q`,
/// choosing the shortest rotation (`|φ| ≤ π`).
fn quat_log(q: Quat) -> [f64; 3] {
    let q = quat_normalize(q);
    // q and −q represent the same rotation; pick w ≥ 0 for |θ| ≤ π.
    let (w, v) = if q[0] < 0.0 {
        (-q[0], [-q[1], -q[2], -q[3]])
    } else {
        (q[0], [q[1], q[2], q[3]])
    };
    let s = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let theta = 2.0 * s.atan2(w);
    // θ / sin(θ/2), with its Taylor series near zero (sin(θ/2) = s).
    let k = if s < 1e-8 {
        2.0 / w * (1.0 - s * s / (3.0 * w * w))
    } else {
        theta / s
    };
    [k * v[0], k * v[1], k * v[2]]
}

/// Rotation matrix (child → parent, row-major) of a unit quaternion.
fn quat_to_rotmat_arr(q: Quat) -> [[f64; 3]; 3] {
    quat_to_rotmat(q[0], q[1], q[2], q[3])
}

/// Convert a unit quaternion `(qw, qx, qy, qz)` to a 3×3 rotation matrix.
///
/// The quaternion must be normalised (`qw² + qx² + qy² + qz² = 1`).
/// Returns a row-major rotation matrix.
fn quat_to_rotmat(qw: f64, qx: f64, qy: f64, qz: f64) -> [[f64; 3]; 3] {
    [
        [
            1.0 - 2.0 * (qy * qy + qz * qz),
            2.0 * (qx * qy - qz * qw),
            2.0 * (qx * qz + qy * qw),
        ],
        [
            2.0 * (qx * qy + qz * qw),
            1.0 - 2.0 * (qx * qx + qz * qz),
            2.0 * (qy * qz - qx * qw),
        ],
        [
            2.0 * (qx * qz - qy * qw),
            2.0 * (qy * qz + qx * qw),
            1.0 - 2.0 * (qx * qx + qy * qy),
        ],
    ]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mat3_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn mat3_t_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

/// Apply the SO(3) left Jacobian `V(φ) = I + b·[φ]× + c·[φ]×²` to `ρ`, where
/// `b = (1 − cos θ)/θ²` and `c = (θ − sin θ)/θ³`. This is the translational
/// part of the SE(3) exponential.
fn so3_left_jacobian_apply(phi: [f64; 3], rho: [f64; 3]) -> [f64; 3] {
    let t2 = phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2];
    let theta = t2.sqrt();
    let (b, c) = if theta < 1e-4 {
        (0.5 - t2 / 24.0, 1.0 / 6.0 - t2 / 120.0)
    } else {
        (
            (1.0 - theta.cos()) / t2,
            (theta - theta.sin()) / (t2 * theta),
        )
    };
    let px = cross3(phi, rho);
    let ppx = cross3(phi, px);
    [
        rho[0] + b * px[0] + c * ppx[0],
        rho[1] + b * px[1] + c * ppx[1],
        rho[2] + b * px[2] + c * ppx[2],
    ]
}

/// Apply the inverse SO(3) left Jacobian
/// `V(φ)⁻¹ = I − ½[φ]× + d·[φ]×²`, `d = (1 − θ·sin θ / (2(1 − cos θ)))/θ²`.
fn so3_left_jacobian_inverse_apply(phi: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let t2 = phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2];
    let theta = t2.sqrt();
    let d = if theta < 1e-4 {
        1.0 / 12.0 + t2 / 720.0
    } else {
        (1.0 - theta * theta.sin() / (2.0 * (1.0 - theta.cos()))) / t2
    };
    let px = cross3(phi, p);
    let ppx = cross3(phi, px);
    [
        p[0] - 0.5 * px[0] + d * ppx[0],
        p[1] - 0.5 * px[1] + d * ppx[1],
        p[2] - 0.5 * px[2] + d * ppx[2],
    ]
}

/// Rotate vector `v` about `axis` by `angle` radians (Rodrigues' rotation formula).
fn rotate_vec3_axis_angle(axis: [f64; 3], angle: f64, v: [f64; 3]) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if n < 1e-14 {
        return v;
    }
    let u = [axis[0] / n, axis[1] / n, axis[2] / n];
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    [
        cos_a * v[0] + sin_a * cross[0] + (1.0 - cos_a) * dot * u[0],
        cos_a * v[1] + sin_a * cross[1] + (1.0 - cos_a) * dot * u[1],
        cos_a * v[2] + sin_a * cross[2] + (1.0 - cos_a) * dot * u[2],
    ]
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn test_revolute_at_zero() {
        let j = RevoluteJoint::new([0.0, 0.0, 1.0]);
        let x = j.transform(&[0.0]);
        // At q=0, transform should be identity rotation
        assert!(
            approx_eq(x.rot[0][0], 1.0, 1e-12),
            "R[0][0] should be 1 at q=0"
        );
        assert!(
            approx_eq(x.rot[1][1], 1.0, 1e-12),
            "R[1][1] should be 1 at q=0"
        );
    }

    #[test]
    fn test_revolute_quarter_turn() {
        let j = RevoluteJoint::new([0.0, 0.0, 1.0]);
        let x = j.transform(&[std::f64::consts::FRAC_PI_2]);
        // 90° around z: R·[1,0,0] = [0,1,0]
        let v = [x.rot[0][0], x.rot[1][0], x.rot[2][0]];
        assert!(
            approx_eq(v[0], 0.0, 1e-12),
            "R*x[0] should be 0, got {:.6}",
            v[0]
        );
        assert!(
            approx_eq(v[1], 1.0, 1e-12),
            "R*x[1] should be 1, got {:.6}",
            v[1]
        );
    }

    #[test]
    fn test_prismatic_motion_subspace() {
        let j = PrismaticJoint::new([1.0, 0.0, 0.0]);
        let s = j.motion_subspace();
        assert_eq!(s.len(), 1);
        assert!(
            approx_eq(s[0].linear[0], 1.0, 1e-12),
            "prismatic S linear[0] should be 1"
        );
        assert!(
            approx_eq(s[0].angular[0], 0.0, 1e-12),
            "prismatic S angular should be 0"
        );
    }

    #[test]
    fn test_fixed_joint() {
        let j = FixedJoint;
        assert_eq!(j.dof(), 0);
        let s = j.motion_subspace();
        assert!(s.is_empty(), "fixed joint has no motion subspace columns");
    }
}
