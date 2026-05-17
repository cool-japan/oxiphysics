// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint types for articulated-body models.
//!
//! Each joint implements the [`Joint`] trait providing:
//! - the joint transform `X_J(q)` from parent to child frame
//! - the motion subspace `S` mapping q̇ → spatial velocity
//! - optional Coriolis term `c_J`

use crate::spatial::{SpatialTransform, SpatialVec};

// ─── Joint trait ─────────────────────────────────────────────────────────────

/// A kinematic joint connecting a parent body to a child body.
///
/// Implementors must be `Send + Sync` for multi-threaded use.
pub trait Joint: Send + Sync {
    /// Number of degrees of freedom (0 for fixed, 1 for revolute/prismatic, etc.).
    fn dof(&self) -> usize;

    /// Joint transform `X_J(q)` from parent joint frame to child joint frame.
    ///
    /// `q` must have length `self.dof()`.
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
/// Parameterised as:
/// - `q[0..3]` = translation `[x, y, z]`
/// - `q[3..7]` = unit quaternion `[w, x, y, z]` (7 position variables)
///
/// Velocity is parameterised as:
/// - `q_dot[0..3]` = linear velocity `[vx, vy, vz]` in parent frame
/// - `q_dot[3..6]` = angular velocity `[ωx, ωy, ωz]` in body frame
///
/// The motion subspace `S` is the 6×6 identity (spatial velocity = S · q̇).
///
/// **Note:** Full quaternion-based floating-base dynamics (quaternion kinematics,
/// constraint forces) are complex. This implementation provides correct DOF count
/// and motion subspace for RNEA/ABA, but quaternion-to-rotation decoding uses
/// only `q[3..7]` which must be a unit quaternion supplied by the caller.
pub struct FreeFloatingJoint;

impl Joint for FreeFloatingJoint {
    fn dof(&self) -> usize {
        6
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        // q[0..3] = translation, q[3..7] = quaternion [w, x, y, z]
        let trans = [
            if !q.is_empty() { q[0] } else { 0.0 },
            if q.len() > 1 { q[1] } else { 0.0 },
            if q.len() > 2 { q[2] } else { 0.0 },
        ];
        let rot = if q.len() >= 7 {
            let (qw, qx, qy, qz) = (q[3], q[4], q[5], q[6]);
            quat_to_rotmat(qw, qx, qy, qz)
        } else {
            // Identity rotation if quaternion not supplied
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        };
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
/// The joint transform is `X_J = R(axis2, q1) ∘ R(axis1, q0)`.
///
/// **Motion subspace:** `S(q)` depends on `q[1]` — the first column (axis1's
/// contribution in the child frame) is rotated by `q[1]` about `axis2`. This
/// is why the [`Joint::motion_subspace_at`] override is needed. The q-independent
/// fallback `motion_subspace()` returns the columns at `q = [0, 0]`.
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
        // Apply r1 first (about axis1), then r2 (about axis2 in intermediate frame).
        // compose: r2 ∘ r1 means: first apply r1, then r2.
        r2.compose(&r1)
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
}

// ─── SphericalJoint ──────────────────────────────────────────────────────────

/// 3-DOF spherical (ball-and-socket) joint for shoulder/hip articulations.
///
/// **Configuration parameterisation:** `q ∈ ℝ⁴` stores a unit quaternion
/// `[w, x, y, z]`. This is NOT accessed via `q_slice` (which only returns
/// `dof() = 3` elements); callers must supply quaternion state separately
/// and pass the raw q array when calling [`Joint::transform`].
///
/// **Velocity parameterisation:** `q_dot ∈ ℝ³` is the body-frame angular
/// velocity, matching the 3-column constant motion subspace `[I₃; 0₃]`.
///
/// This mirrors the convention used by [`FreeFloatingJoint`] for its rotational
/// sub-motion.
pub struct SphericalJoint;

impl Joint for SphericalJoint {
    fn dof(&self) -> usize {
        3
    }

    fn transform(&self, q: &[f64]) -> SpatialTransform {
        let rot = if q.len() >= 4 {
            let (qw, qx, qy, qz) = (q[0], q[1], q[2], q[3]);
            quat_to_rotmat(qw, qx, qy, qz)
        } else {
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        };
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

// ─── Helper: quaternion to rotation matrix ────────────────────────────────────

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
