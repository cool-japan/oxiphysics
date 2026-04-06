#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constraint solvers for the OxiPhysics engine.
//!
//! Provides contact constraints (sequential impulse), joint constraints
//! (fixed, revolute, prismatic, ball, spring), PGS and TGS solvers,
//! and an island manager for grouping connected bodies.
#![warn(missing_docs)]
#![allow(ambiguous_glob_reexports)]
#![allow(dead_code)]

mod error;
pub use error::*;

mod traits;
pub use traits::Constraint;

pub mod contact;
pub use contact::ContactConstraint;

pub mod pgs_solver;
pub use pgs_solver::PgsSolver;

pub mod tgs_solver;
pub use tgs_solver::TgsSolver;

pub mod joints;
pub use joints::{
    BallJoint, FixedJoint, GearJoint, MotorJoint, PrismaticJoint, PulleyJoint, RackPinionJoint,
    RevoluteJoint, SpringJoint,
};

pub mod islands;

pub mod island;
pub use island::{Island, IslandManager};

pub mod pbd;

pub mod warm_start;

pub mod pid_motor;

pub mod control_theory;
pub use control_theory::*;

pub mod friction;
pub use friction::*;

pub mod ccd_constraints;
pub use ccd_constraints::{CcdConstraint, CcdConstraintSolver, estimate_linear_toi};

pub mod motor_constraints;
pub use motor_constraints::{
    AngularMotorConstraint, LinearMotorConstraint, MotorConstraint, MotorPid, MotorSolver,
    ServoMotorConstraint,
};

pub mod six_dof_constraint;
pub use six_dof_constraint::{AxisConfig, MotorConfig, SixDofConstraint, compute_jacobian_6dof};

// ── Constraint utility functions ────────────────────────────────────────────

/// Compute the Baumgarte bias velocity for a positional error.
///
/// bias = -(beta / dt) * error
///
/// * `error` - The positional error (scalar, e.g. penetration depth).
/// * `beta` - Baumgarte stabilization factor (typically 0.1 to 0.3).
/// * `dt` - Time step.
#[allow(dead_code)]
pub fn baumgarte_bias(error: f64, beta: f64, dt: f64) -> f64 {
    -(beta / dt) * error
}

/// Compute the effective mass for a 1-DOF constraint between two bodies.
///
/// effective_mass = 1 / (inv_mass_a + inv_mass_b + ...)
///
/// For a simple point constraint along a normal `n`:
/// K = inv_mass_a + inv_mass_b + n^T * (inv_I_a * (r_a x n) x r_a + inv_I_b * (r_b x n) x r_b)
///
/// This simplified version considers only the linear part.
#[allow(dead_code)]
pub fn effective_mass_linear(inv_mass_a: f64, inv_mass_b: f64) -> f64 {
    let sum = inv_mass_a + inv_mass_b;
    if sum.abs() < 1e-15 { 0.0 } else { 1.0 / sum }
}

/// Clamp an impulse to a range `[lo, hi]`.
#[allow(dead_code)]
pub fn clamp_impulse(impulse: f64, lo: f64, hi: f64) -> f64 {
    impulse.max(lo).min(hi)
}

/// Compute the velocity constraint error (Cdot) for a distance constraint.
///
/// Cdot = n . (v_b + omega_b x r_b - v_a - omega_a x r_a)
///
/// All vectors as `[f64; 3]` arrays to avoid nalgebra dependency.
#[allow(dead_code)]
pub fn distance_cdot(
    normal: [f64; 3],
    v_a: [f64; 3],
    omega_a: [f64; 3],
    r_a: [f64; 3],
    v_b: [f64; 3],
    omega_b: [f64; 3],
    r_b: [f64; 3],
) -> f64 {
    let cross_a = cross3(omega_a, r_a);
    let cross_b = cross3(omega_b, r_b);
    let rel = [
        v_b[0] + cross_b[0] - v_a[0] - cross_a[0],
        v_b[1] + cross_b[1] - v_a[1] - cross_a[1],
        v_b[2] + cross_b[2] - v_a[2] - cross_a[2],
    ];
    dot3(normal, rel)
}

/// Simple 3D cross product on arrays.
#[allow(dead_code)]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Simple 3D dot product on arrays.
#[allow(dead_code)]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the angular error between two orientations represented as
/// rotation matrices (3x3 column-major arrays).
///
/// Returns the axis-angle error vector `[ex, ey, ez]` such that rotating
/// body A by this vector (in world space) would align it with body B.
///
/// Uses the skew-symmetric part of R_b * R_a^T.
#[allow(dead_code)]
pub fn angular_error_3x3(rot_a: [[f64; 3]; 3], rot_b: [[f64; 3]; 3]) -> [f64; 3] {
    // R_err = R_b * R_a^T
    let rat = transpose3(rot_a);
    let r_err = mat_mul3(rot_b, rat);
    // Extract axis-angle from skew-symmetric part
    [
        0.5 * (r_err[2][1] - r_err[1][2]),
        0.5 * (r_err[0][2] - r_err[2][0]),
        0.5 * (r_err[1][0] - r_err[0][1]),
    ]
}

#[allow(dead_code)]
fn transpose3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

#[allow(dead_code)]
fn mat_mul3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

/// Compute the constraint force magnitude from a Lagrange multiplier
/// and the time step: F = lambda / dt.
#[allow(dead_code)]
pub fn constraint_force_from_lambda(lambda: f64, dt: f64) -> f64 {
    if dt.abs() < 1e-15 { 0.0 } else { lambda / dt }
}

/// Sequential impulse iteration: apply one impulse step for a 1-DOF
/// constraint and return the impulse applied.
///
/// * `cdot` - Current velocity error.
/// * `eff_mass` - Effective mass.
/// * `accumulated` - Current accumulated impulse (mutated in place).
/// * `lo` / `hi` - Impulse clamp bounds.
///
/// Returns the delta impulse that was applied.
#[allow(dead_code)]
pub fn solve_1dof_impulse(
    cdot: f64,
    eff_mass: f64,
    accumulated: &mut f64,
    lo: f64,
    hi: f64,
) -> f64 {
    let raw = -eff_mass * cdot;
    let old = *accumulated;
    *accumulated = (old + raw).max(lo).min(hi);
    *accumulated - old
}

#[cfg(test)]
mod constraint_util_tests {
    use super::*;

    #[test]
    fn test_baumgarte_bias_positive_error() {
        let bias = baumgarte_bias(0.1, 0.2, 0.01);
        // -(0.2 / 0.01) * 0.1 = -2.0
        assert!((bias + 2.0).abs() < 1e-10, "bias={bias}");
    }

    #[test]
    fn test_baumgarte_bias_zero_error() {
        let bias = baumgarte_bias(0.0, 0.2, 0.01);
        assert!(bias.abs() < 1e-10);
    }

    #[test]
    fn test_effective_mass_linear_equal() {
        let m = effective_mass_linear(1.0, 1.0);
        assert!((m - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_effective_mass_linear_one_static() {
        // If body B is static (inv_mass = 0), effective mass = 1 / inv_mass_a
        let m = effective_mass_linear(0.5, 0.0);
        assert!((m - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_effective_mass_linear_both_static() {
        let m = effective_mass_linear(0.0, 0.0);
        assert!(m.abs() < 1e-10);
    }

    #[test]
    fn test_clamp_impulse_within_range() {
        assert!((clamp_impulse(0.5, 0.0, 1.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_clamp_impulse_below_min() {
        assert!((clamp_impulse(-1.0, 0.0, 1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_clamp_impulse_above_max() {
        assert!((clamp_impulse(2.0, 0.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_cdot_stationary() {
        let n = [1.0, 0.0, 0.0];
        let zero = [0.0, 0.0, 0.0];
        let cdot = distance_cdot(n, zero, zero, zero, zero, zero, zero);
        assert!(cdot.abs() < 1e-10);
    }

    #[test]
    fn test_distance_cdot_relative_velocity() {
        let n = [1.0, 0.0, 0.0];
        let va = [0.0, 0.0, 0.0];
        let vb = [3.0, 0.0, 0.0];
        let zero = [0.0, 0.0, 0.0];
        let cdot = distance_cdot(n, va, zero, zero, vb, zero, zero);
        assert!((cdot - 3.0).abs() < 1e-10, "cdot={cdot}");
    }

    #[test]
    fn test_angular_error_identity() {
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let err = angular_error_3x3(i, i);
        for e in err {
            assert!(e.abs() < 1e-10);
        }
    }

    #[test]
    fn test_angular_error_small_rotation() {
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Small rotation around Z axis by angle theta
        let theta = 0.01_f64;
        let c = theta.cos();
        let s = theta.sin();
        let rz = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let err = angular_error_3x3(i, rz);
        // Error should be approximately [0, 0, theta]
        assert!(err[0].abs() < 1e-6);
        assert!(err[1].abs() < 1e-6);
        assert!((err[2] - theta).abs() < 1e-4, "ez={}", err[2]);
    }

    #[test]
    fn test_constraint_force_from_lambda() {
        let f = constraint_force_from_lambda(10.0, 0.01);
        assert!((f - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_constraint_force_zero_dt() {
        let f = constraint_force_from_lambda(10.0, 0.0);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_solve_1dof_impulse_basic() {
        let mut acc = 0.0;
        let delta = solve_1dof_impulse(1.0, 2.0, &mut acc, -100.0, 100.0);
        // raw = -2.0 * 1.0 = -2.0
        assert!((delta + 2.0).abs() < 1e-10, "delta={delta}");
        assert!((acc + 2.0).abs() < 1e-10, "acc={acc}");
    }

    #[test]
    fn test_solve_1dof_impulse_clamped() {
        let mut acc = 0.0;
        let delta = solve_1dof_impulse(1.0, 2.0, &mut acc, -1.0, 1.0);
        // raw = -2.0, clamped to -1.0
        assert!((delta + 1.0).abs() < 1e-10, "delta={delta}");
        assert!((acc + 1.0).abs() < 1e-10, "acc={acc}");
    }

    #[test]
    fn test_cross3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[0]).abs() < 1e-10);
        assert!((c[1]).abs() < 1e-10);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dot3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot3(a, b) - 32.0).abs() < 1e-10);
    }

    #[test]
    fn test_transpose3() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = transpose3(m);
        assert!((t[0][1] - 4.0).abs() < 1e-10);
        assert!((t[1][0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mat_mul3_identity() {
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let r = mat_mul3(i, m);
        for row in 0..3 {
            for col in 0..3 {
                assert!(
                    (r[row][col] - m[row][col]).abs() < 1e-10,
                    "r[{row}][{col}]={} expected {}",
                    r[row][col],
                    m[row][col]
                );
            }
        }
    }
}

// ── Impulse / Jacobian Utilities ─────────────────────────────────────────────

/// Compute the rotational contribution to effective mass for a single body.
///
/// `r × n` dotted through `I⁻¹`, then dotted with `r × n` again:
///
/// ```text
/// K_rot = (r × n)^T I⁻¹ (r × n)
/// ```
///
/// * `r` - lever arm from centre of mass to contact point (world space).
/// * `n` - constraint axis / contact normal (world space, unit vector).
/// * `inv_inertia_world` - 3×3 world-space inverse inertia tensor (column-major).
#[allow(dead_code)]
pub fn rotational_effective_mass_contribution(
    r: [f64; 3],
    n: [f64; 3],
    inv_inertia_world: [[f64; 3]; 3],
) -> f64 {
    let rxn = cross3(r, n);
    let iinv_rxn = mat_vec_mul3(inv_inertia_world, rxn);
    dot3(rxn, iinv_rxn)
}

/// Multiply a 3×3 matrix by a 3-vector.
#[allow(dead_code)]
fn mat_vec_mul3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Compute the full scalar effective mass K for a contact or joint constraint.
///
/// ```text
/// K = invMa + invMb
///   + (ra × n)^T Ia⁻¹ (ra × n)
///   + (rb × n)^T Ib⁻¹ (rb × n)
/// ```
///
/// Returns `1/K` (the effective mass used to scale impulses), or 0.0 if K ≈ 0.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn effective_mass_full(
    inv_mass_a: f64,
    inv_mass_b: f64,
    r_a: [f64; 3],
    r_b: [f64; 3],
    normal: [f64; 3],
    inv_inertia_a: [[f64; 3]; 3],
    inv_inertia_b: [[f64; 3]; 3],
) -> f64 {
    let k = inv_mass_a
        + inv_mass_b
        + rotational_effective_mass_contribution(r_a, normal, inv_inertia_a)
        + rotational_effective_mass_contribution(r_b, normal, inv_inertia_b);
    if k.abs() < 1e-15 { 0.0 } else { 1.0 / k }
}

/// Compute the relative velocity at a contact point.
///
/// ```text
/// v_rel = (vb + ωb × rb) - (va + ωa × ra)
/// ```
#[allow(dead_code)]
pub fn relative_velocity(
    v_a: [f64; 3],
    omega_a: [f64; 3],
    r_a: [f64; 3],
    v_b: [f64; 3],
    omega_b: [f64; 3],
    r_b: [f64; 3],
) -> [f64; 3] {
    let vel_a = [
        v_a[0] + cross3(omega_a, r_a)[0],
        v_a[1] + cross3(omega_a, r_a)[1],
        v_a[2] + cross3(omega_a, r_a)[2],
    ];
    let vel_b = [
        v_b[0] + cross3(omega_b, r_b)[0],
        v_b[1] + cross3(omega_b, r_b)[1],
        v_b[2] + cross3(omega_b, r_b)[2],
    ];
    [
        vel_b[0] - vel_a[0],
        vel_b[1] - vel_a[1],
        vel_b[2] - vel_a[2],
    ]
}

/// Compute the restitution bias for a contact constraint.
///
/// The bias adds a target velocity based on the incoming relative normal velocity
/// (coefficient-of-restitution model):
///
/// ```text
/// bias = max(0, -e * v_n)   if |v_n| > restitution_threshold
/// ```
#[allow(dead_code)]
pub fn restitution_bias(v_rel_normal: f64, restitution: f64, threshold: f64) -> f64 {
    if v_rel_normal < -threshold {
        -restitution * v_rel_normal
    } else {
        0.0
    }
}

/// Apply an impulse to a body, returning updated linear and angular velocities.
///
/// Linear:  v' = v + J * inv_mass
/// Angular: ω' = ω + I⁻¹ * (r × J)
///
/// `impulse_dir` should be a signed magnitude times the constraint direction.
#[allow(dead_code)]
pub fn apply_impulse(
    vel: [f64; 3],
    omega: [f64; 3],
    inv_mass: f64,
    inv_inertia: [[f64; 3]; 3],
    r: [f64; 3],
    impulse_dir: [f64; 3],
    sign: f64,
) -> ([f64; 3], [f64; 3]) {
    let j = [
        sign * impulse_dir[0],
        sign * impulse_dir[1],
        sign * impulse_dir[2],
    ];
    let new_vel = [
        vel[0] + inv_mass * j[0],
        vel[1] + inv_mass * j[1],
        vel[2] + inv_mass * j[2],
    ];
    let torque = cross3(r, j);
    let delta_omega = mat_vec_mul3(inv_inertia, torque);
    let new_omega = [
        omega[0] + delta_omega[0],
        omega[1] + delta_omega[1],
        omega[2] + delta_omega[2],
    ];
    (new_vel, new_omega)
}

// ── Quaternion Utilities ──────────────────────────────────────────────────────

/// Normalize a quaternion `[x, y, z, w]`.
///
/// Returns the input unchanged if its magnitude is below `1e-15`.
#[allow(dead_code)]
pub fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let len2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if len2 < 1e-30 {
        return q;
    }
    let inv_len = 1.0 / len2.sqrt();
    [
        q[0] * inv_len,
        q[1] * inv_len,
        q[2] * inv_len,
        q[3] * inv_len,
    ]
}

/// Multiply two quaternions `p * q` (Hamilton product).
///
/// Quaternion layout: `[x, y, z, w]`.
#[allow(dead_code)]
pub fn quat_mul(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    let [px, py, pz, pw] = p;
    let [qx, qy, qz, qw] = q;
    [
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
        pw * qw - px * qx - py * qy - pz * qz,
    ]
}

/// Conjugate (inverse for unit quaternions) of `[x, y, z, w]`.
#[allow(dead_code)]
pub fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}

/// Convert a quaternion `[x, y, z, w]` to a 3×3 rotation matrix (row-major).
#[allow(dead_code)]
pub fn quat_to_rot3(q: [f64; 4]) -> [[f64; 3]; 3] {
    let [x, y, z, w] = q;
    let x2 = 2.0 * x * x;
    let y2 = 2.0 * y * y;
    let z2 = 2.0 * z * z;
    let xy = 2.0 * x * y;
    let xz = 2.0 * x * z;
    let yz = 2.0 * y * z;
    let wx = 2.0 * w * x;
    let wy = 2.0 * w * y;
    let wz = 2.0 * w * z;
    [
        [1.0 - y2 - z2, xy - wz, xz + wy],
        [xy + wz, 1.0 - x2 - z2, yz - wx],
        [xz - wy, yz + wx, 1.0 - x2 - y2],
    ]
}

/// Compute the angular velocity from the time derivative of a quaternion.
///
/// `ω = 2 * conj(q) ⊗ dq/dt`  (body-space angular velocity)
#[allow(dead_code)]
pub fn quat_to_angular_velocity(q: [f64; 4], q_dot: [f64; 4]) -> [f64; 3] {
    let q_conj = quat_conjugate(q);
    let omega_q = quat_mul(q_conj, q_dot);
    [2.0 * omega_q[0], 2.0 * omega_q[1], 2.0 * omega_q[2]]
}

// ── Vector Utilities ──────────────────────────────────────────────────────────

/// Normalize a 3-vector. Returns `[0,0,0]` if length < 1e-15.
#[allow(dead_code)]
pub fn vec3_normalize(v: [f64; 3]) -> [f64; 3] {
    let len2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if len2 < 1e-30 {
        return [0.0, 0.0, 0.0];
    }
    let inv = 1.0 / len2.sqrt();
    [v[0] * inv, v[1] * inv, v[2] * inv]
}

/// Length of a 3-vector.
#[allow(dead_code)]
pub fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Scale a 3-vector.
#[allow(dead_code)]
pub fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Add two 3-vectors.
#[allow(dead_code)]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors (a - b).
#[allow(dead_code)]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Compute a tangent vector perpendicular to `n`.
///
/// Picks the most orthogonal world-axis and orthogonalizes against `n`.
#[allow(dead_code)]
pub fn perpendicular_tangent(n: [f64; 3]) -> [f64; 3] {
    let nx_abs = n[0].abs();
    let ny_abs = n[1].abs();
    let nz_abs = n[2].abs();

    let candidate = if nx_abs <= ny_abs && nx_abs <= nz_abs {
        [1.0, 0.0, 0.0]
    } else if ny_abs <= nx_abs && ny_abs <= nz_abs {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };

    let t = cross3(n, candidate);
    vec3_normalize(t)
}

// ── Solver Convergence Helpers ────────────────────────────────────────────────

/// Compute the L2 norm of a residual vector.
#[allow(dead_code)]
pub fn residual_norm(residuals: &[f64]) -> f64 {
    residuals.iter().map(|r| r * r).sum::<f64>().sqrt()
}

/// Compute the L∞ (max) norm of a residual vector.
#[allow(dead_code)]
pub fn residual_max(residuals: &[f64]) -> f64 {
    residuals
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max)
}

/// Determine the over-relaxation factor for SOR (Successive Over-Relaxation).
///
/// Standard PGS uses ω = 1.0 (no relaxation).  Values in (1, 2) provide
/// over-relaxation; values in (0, 1) provide under-relaxation.
#[allow(dead_code)]
pub fn sor_omega(omega: f64) -> f64 {
    omega.clamp(1e-4, 2.0 - 1e-4)
}

// ── Split-Impulse / Pseudo-Velocity Helpers ───────────────────────────────────

/// Compute the split-impulse pseudo-velocity bias.
///
/// Split impulse separates position correction from velocity correction by
/// introducing a pseudo-velocity `v*` that corrects drift without affecting
/// the true velocity.
///
/// ```text
/// v*_bias = -(beta / dt) * max(error - slop, 0)
/// ```
///
/// * `error` - Positional penetration depth (positive = overlapping).
/// * `beta`  - Position stabilization factor (0.1–0.3 typical).
/// * `dt`    - Time step.
/// * `slop`  - Allowed penetration slop before correction kicks in.
#[allow(dead_code)]
pub fn split_impulse_bias(error: f64, beta: f64, dt: f64, slop: f64) -> f64 {
    let correctable = (error - slop).max(0.0);
    -(beta / dt) * correctable
}

// ── Jacobian Utilities ────────────────────────────────────────────────────────

/// Build the 6-component Jacobian row for a point-contact normal constraint.
///
/// The constraint is: `n · (v_b + ω_b × r_b - v_a - ω_a × r_a) = 0`
///
/// Returns `[J_va, J_ωa, J_vb, J_ωb]` as four 3-vectors.
#[allow(dead_code)]
pub fn contact_normal_jacobian(
    n: [f64; 3],
    r_a: [f64; 3],
    r_b: [f64; 3],
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // J_va = -n
    let j_va = [-n[0], -n[1], -n[2]];
    // J_ωa = -(r_a × n)
    let rxn_a = cross3(r_a, n);
    let j_wa = [-rxn_a[0], -rxn_a[1], -rxn_a[2]];
    // J_vb = +n
    let j_vb = n;
    // J_ωb = +(r_b × n)
    let j_wb = cross3(r_b, n);
    (j_va, j_wa, j_vb, j_wb)
}

// ── Interpolation / Blending ──────────────────────────────────────────────────

/// Linear interpolation between two scalars.
#[allow(dead_code)]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Smoothstep (cubic Hermite) interpolation.
///
/// Maps `t` in `[0, 1]` to a smooth curve: `3t² - 2t³`.
#[allow(dead_code)]
pub fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Clamp a value to `[lo, hi]`.
#[allow(dead_code)]
pub fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ── Extended Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod extended_util_tests {
    use super::*;

    #[test]
    fn test_rotational_effective_mass_zero_lever() {
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let contrib = rotational_effective_mass_contribution([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], i);
        assert!(contrib.abs() < 1e-10);
    }

    #[test]
    fn test_rotational_effective_mass_unit() {
        // r = (1,0,0), n = (0,0,1), I⁻¹ = identity
        // r × n = (0,1,0), I⁻¹ * (0,1,0) = (0,1,0), dot = 1
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let contrib = rotational_effective_mass_contribution([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], i);
        assert!((contrib - 1.0).abs() < 1e-10, "contrib={contrib}");
    }

    #[test]
    fn test_effective_mass_full_equal_masses() {
        let i = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let m = effective_mass_full(1.0, 1.0, [0.0; 3], [0.0; 3], [1.0, 0.0, 0.0], i, i);
        // K = 1 + 1 + 0 + 0 = 2, eff_mass = 0.5
        assert!((m - 0.5).abs() < 1e-10, "m={m}");
    }

    #[test]
    fn test_relative_velocity_stationary() {
        let z = [0.0; 3];
        let rv = relative_velocity(z, z, z, z, z, z);
        for c in rv {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_relative_velocity_linear() {
        let z = [0.0; 3];
        let va = [1.0, 0.0, 0.0];
        let vb = [4.0, 0.0, 0.0];
        let rv = relative_velocity(va, z, z, vb, z, z);
        assert!((rv[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_restitution_bias_above_threshold() {
        // v_rel_n = -2.0, e = 0.5, threshold = 0.1
        let bias = restitution_bias(-2.0, 0.5, 0.1);
        assert!((bias - 1.0).abs() < 1e-10, "bias={bias}");
    }

    #[test]
    fn test_restitution_bias_below_threshold() {
        // Small incoming velocity — no restitution applied
        let bias = restitution_bias(-0.05, 0.5, 0.1);
        assert!(bias.abs() < 1e-10, "bias={bias}");
    }

    #[test]
    fn test_restitution_bias_positive_vel() {
        // Bodies separating — no restitution
        let bias = restitution_bias(0.5, 0.5, 0.1);
        assert!(bias.abs() < 1e-10);
    }

    #[test]
    fn test_quat_normalize_unit() {
        let q = [0.0, 0.0, 0.0, 1.0];
        let qn = quat_normalize(q);
        assert!((qn[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quat_normalize_nontrivial() {
        let q = [1.0, 1.0, 1.0, 1.0];
        let qn = quat_normalize(q);
        let len2 = qn[0] * qn[0] + qn[1] * qn[1] + qn[2] * qn[2] + qn[3] * qn[3];
        assert!((len2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quat_mul_identity() {
        let id = [0.0, 0.0, 0.0, 1.0];
        let q = [0.1, 0.2, 0.3, 0.9];
        let qn = quat_normalize(q);
        let r = quat_mul(id, qn);
        for i in 0..4 {
            assert!((r[i] - qn[i]).abs() < 1e-10, "r[{i}]={}", r[i]);
        }
    }

    #[test]
    fn test_quat_conjugate() {
        let q = [0.1, 0.2, 0.3, 0.9274];
        let c = quat_conjugate(q);
        assert!((c[0] + q[0]).abs() < 1e-10);
        assert!((c[3] - q[3]).abs() < 1e-10);
    }

    #[test]
    fn test_quat_to_rot3_identity() {
        let id = [0.0, 0.0, 0.0, 1.0];
        let r = quat_to_rot3(id);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (r[i][j] - expected).abs() < 1e-10,
                    "r[{i}][{j}]={}",
                    r[i][j]
                );
            }
        }
    }

    #[test]
    fn test_vec3_normalize() {
        let v = [3.0, 0.0, 0.0];
        let n = vec3_normalize(v);
        assert!((n[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_normalize_zero() {
        let v = [0.0, 0.0, 0.0];
        let n = vec3_normalize(v);
        for c in n {
            assert!(c.abs() < 1e-10);
        }
    }

    #[test]
    fn test_vec3_len() {
        let v = [3.0, 4.0, 0.0];
        assert!((vec3_len(v) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_perpendicular_tangent() {
        let n = [1.0, 0.0, 0.0];
        let t = perpendicular_tangent(n);
        // t must be perpendicular to n
        assert!(dot3(n, t).abs() < 1e-10, "dot={}", dot3(n, t));
        // t must be unit length
        assert!((vec3_len(t) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_perpendicular_tangent_y() {
        let n = [0.0, 1.0, 0.0];
        let t = perpendicular_tangent(n);
        assert!(dot3(n, t).abs() < 1e-10);
        assert!((vec3_len(t) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_residual_norm_basic() {
        let r = [3.0, 4.0];
        assert!((residual_norm(&r) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_residual_max() {
        let r = [1.0, -5.0, 3.0];
        assert!((residual_max(&r) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_sor_omega_clamp() {
        assert!((sor_omega(1.0) - 1.0).abs() < 1e-10);
        assert!(sor_omega(0.0) > 0.0);
        assert!(sor_omega(3.0) < 2.0);
    }

    #[test]
    fn test_split_impulse_bias_no_slop() {
        let bias = split_impulse_bias(0.1, 0.2, 0.01, 0.0);
        // -(0.2 / 0.01) * 0.1 = -2.0
        assert!((bias + 2.0).abs() < 1e-10, "bias={bias}");
    }

    #[test]
    fn test_split_impulse_bias_within_slop() {
        // error = 0.001 < slop = 0.01 → no correction
        let bias = split_impulse_bias(0.001, 0.2, 0.01, 0.01);
        assert!(bias.abs() < 1e-10, "bias={bias}");
    }

    #[test]
    fn test_contact_normal_jacobian() {
        let n = [1.0, 0.0, 0.0];
        let r_a = [0.0, 1.0, 0.0];
        let r_b = [0.0, -1.0, 0.0];
        let (j_va, _j_wa, j_vb, _j_wb) = contact_normal_jacobian(n, r_a, r_b);
        // J_va = -n
        assert!((j_va[0] + 1.0).abs() < 1e-10);
        // J_vb = +n
        assert!((j_vb[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lerp() {
        assert!((lerp(0.0, 10.0, 0.5) - 5.0).abs() < 1e-10);
        assert!((lerp(0.0, 10.0, 0.0) - 0.0).abs() < 1e-10);
        assert!((lerp(0.0, 10.0, 1.0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_smoothstep() {
        // At t=0 and t=1, smoothstep returns 0 and 1 respectively.
        assert!(smoothstep(0.0).abs() < 1e-10);
        assert!((smoothstep(1.0) - 1.0).abs() < 1e-10);
        // At t=0.5 the value is 0.5.
        assert!((smoothstep(0.5) - 0.5).abs() < 1e-10);
        // Values outside [0,1] are clamped.
        assert!(smoothstep(-1.0).abs() < 1e-10);
        assert!((smoothstep(2.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clamp() {
        assert!((clamp(5.0, 0.0, 10.0) - 5.0).abs() < 1e-10);
        assert!((clamp(-1.0, 0.0, 10.0) - 0.0).abs() < 1e-10);
        assert!((clamp(11.0, 0.0, 10.0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_impulse_linear_only() {
        // Static inertia (zeros) so only linear velocity changes.
        let zero_inertia = [[0.0; 3]; 3];
        let vel = [1.0, 0.0, 0.0];
        let omega = [0.0; 3];
        let r = [0.0; 3];
        let dir = [1.0, 0.0, 0.0];
        let (new_vel, new_omega) = apply_impulse(vel, omega, 0.5, zero_inertia, r, dir, 1.0);
        assert!(
            (new_vel[0] - 1.5).abs() < 1e-10,
            "new_vel[0]={}",
            new_vel[0]
        );
        for c in new_omega {
            assert!(c.abs() < 1e-10);
        }
    }
}
pub mod cable_constraints;
pub mod contact_constraints;
pub mod contact_stability;
pub mod contact_theory;
pub mod energy_management_constraints;
pub mod game_physics_constraints;
pub mod haptic_constraints;
pub mod holonomic_constraints;
pub mod island_solver;
pub mod locomotion;
pub mod locomotion_constraints;
pub mod mechanism;
pub mod motion_planning;
pub mod multi_agent_constraints;
pub mod multibody;
pub mod optimal_control;
pub mod optimization_constraints;
pub mod robot_control;
pub mod servo_constraints;
pub mod swarm_robotics;
pub mod trajectory_optimization;
pub mod variational_constraints;
pub mod wearable_constraints;
