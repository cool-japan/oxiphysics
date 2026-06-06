// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial algebra (Featherstone) for articulated rigid body dynamics.
//!
//! Implements 6-dimensional spatial vectors following Featherstone's notation:
//! velocity `[ω; v]`, force `[τ; f]`, and the 6×6 spatial inertia tensor.
//! All computations use plain `f64` arrays — no external algebra library.

/// A 6-dimensional spatial velocity vector `[ω; v]`.
///
/// The first three components are angular velocity `ω` (rad/s);
/// the last three are linear velocity `v` (m/s).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialVelocity {
    /// Angular velocity components `[ωx, ωy, ωz]` in rad/s.
    pub omega: [f64; 3],
    /// Linear velocity components `[vx, vy, vz]` in m/s.
    pub v: [f64; 3],
}

impl SpatialVelocity {
    /// Create a spatial velocity from angular and linear parts.
    pub fn new(omega: [f64; 3], v: [f64; 3]) -> Self {
        Self { omega, v }
    }

    /// Zero spatial velocity.
    pub fn zero() -> Self {
        Self {
            omega: [0.0; 3],
            v: [0.0; 3],
        }
    }

    /// Return the full 6-vector `[ωx, ωy, ωz, vx, vy, vz]`.
    pub fn as_array(&self) -> [f64; 6] {
        [
            self.omega[0],
            self.omega[1],
            self.omega[2],
            self.v[0],
            self.v[1],
            self.v[2],
        ]
    }
}

/// A 6-dimensional spatial force vector `[τ; f]`.
///
/// The first three components are torque `τ` (N·m);
/// the last three are force `f` (N).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialForce {
    /// Torque components `[τx, τy, τz]` in N·m.
    pub tau: [f64; 3],
    /// Force components `[fx, fy, fz]` in N.
    pub f: [f64; 3],
}

impl SpatialForce {
    /// Create a spatial force from torque and force parts.
    pub fn new(tau: [f64; 3], f: [f64; 3]) -> Self {
        Self { tau, f }
    }

    /// Zero spatial force.
    pub fn zero() -> Self {
        Self {
            tau: [0.0; 3],
            f: [0.0; 3],
        }
    }

    /// Return the full 6-vector `[τx, τy, τz, fx, fy, fz]`.
    pub fn as_array(&self) -> [f64; 6] {
        [
            self.tau[0],
            self.tau[1],
            self.tau[2],
            self.f[0],
            self.f[1],
            self.f[2],
        ]
    }
}

/// A 6×6 spatial inertia tensor (articulated-body inertia).
///
/// Stored in row-major order. Symmetric positive-definite for valid bodies.
#[derive(Debug, Clone, Copy)]
pub struct SpatialInertia {
    /// Row-major 6×6 matrix entries.
    pub data: [[f64; 6]; 6],
}

impl SpatialInertia {
    /// Create from a row-major 6×6 array.
    pub fn new(data: [[f64; 6]; 6]) -> Self {
        Self { data }
    }

    /// Identity spatial inertia (unit mass, unit rotational inertia).
    pub fn identity() -> Self {
        let mut d = [[0.0f64; 6]; 6];
        for (i, row) in d.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        Self { data: d }
    }

    /// Build a rigid-body spatial inertia given scalar mass and a diagonal
    /// rotational inertia tensor `[Ixx, Iyy, Izz]`.
    ///
    /// The body frame is assumed coincident with the centre of mass so the
    /// off-diagonal coupling blocks are zero.
    pub fn from_rigid_body(mass: f64, inertia_diag: [f64; 3]) -> Self {
        // Structure:
        //  | I_rot   0     |
        //  | 0       m·Id  |
        let mut d = [[0.0f64; 6]; 6];
        d[0][0] = inertia_diag[0];
        d[1][1] = inertia_diag[1];
        d[2][2] = inertia_diag[2];
        d[3][3] = mass;
        d[4][4] = mass;
        d[5][5] = mass;
        Self { data: d }
    }

    /// Multiply the spatial inertia by a spatial velocity, yielding a spatial
    /// momentum / force.
    pub fn mul_velocity(&self, vel: &SpatialVelocity) -> SpatialForce {
        let v = vel.as_array();
        let mut result = [0.0f64; 6];
        for (res_i, row) in result.iter_mut().zip(self.data.iter()) {
            for (d_ij, v_j) in row.iter().zip(v.iter()) {
                *res_i += d_ij * v_j;
            }
        }
        SpatialForce {
            tau: [result[0], result[1], result[2]],
            f: [result[3], result[4], result[5]],
        }
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Spatial cross-product of two spatial velocities.
///
/// Given `V1 = [ω1; v1]` and `V2 = [ω2; v2]`, returns:
/// ```text
/// V1 × V2 = [ω1 × ω2 ;  ω1 × v2 + v1 × ω2]
/// ```
/// This is the Lie bracket of two spatial velocities.
pub fn spatial_cross_product(v1: &SpatialVelocity, v2: &SpatialVelocity) -> SpatialVelocity {
    let omega_out = cross3(v1.omega, v2.omega);
    let v_out_a = cross3(v1.omega, v2.v);
    let v_out_b = cross3(v1.v, v2.omega);
    let v_out = [
        v_out_a[0] + v_out_b[0],
        v_out_a[1] + v_out_b[1],
        v_out_a[2] + v_out_b[2],
    ];
    SpatialVelocity::new(omega_out, v_out)
}

/// Spatial dot product (mechanical power): `F · V = τ·ω + f·v`.
///
/// Returns instantaneous power in watts.
pub fn spatial_dot(force: &SpatialForce, vel: &SpatialVelocity) -> f64 {
    dot3(force.tau, vel.omega) + dot3(force.f, vel.v)
}

/// Plücker transform of a spatial velocity from frame A to frame B.
///
/// Given the position of frame B's origin expressed in frame A (`r_ab`,
/// in metres) and the 3×3 rotation matrix `rot` that rotates vectors from
/// frame A coordinates to frame B coordinates (row-major), compute the
/// transformed spatial velocity.
///
/// Transform law:
/// ```text
/// ω_B = R · ω_A
/// v_B = R · (v_A - r_ab × ω_A)
/// ```
pub fn plucker_transform(
    vel: &SpatialVelocity,
    rot: &[[f64; 3]; 3],
    r_ab: [f64; 3],
) -> SpatialVelocity {
    // Rotate omega
    let omega_a = vel.omega;
    let v_a = vel.v;

    let omega_b = mat3_mul_vec(rot, omega_a);

    // v_B = R * (v_A - r_ab × omega_A)
    let r_cross_omega = cross3(r_ab, omega_a);
    let v_shifted = [
        v_a[0] - r_cross_omega[0],
        v_a[1] - r_cross_omega[1],
        v_a[2] - r_cross_omega[2],
    ];
    let v_b = mat3_mul_vec(rot, v_shifted);

    SpatialVelocity::new(omega_b, v_b)
}

fn mat3_mul_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Add two spatial inertias (composite-body algorithm step).
///
/// Composite inertia of a sub-tree is the sum of all constituent inertias.
pub fn composite_inertia(a: &SpatialInertia, b: &SpatialInertia) -> SpatialInertia {
    let mut d = [[0.0f64; 6]; 6];
    for (d_row, (a_row, b_row)) in d.iter_mut().zip(a.data.iter().zip(b.data.iter())) {
        for (d_ij, (a_ij, b_ij)) in d_row.iter_mut().zip(a_row.iter().zip(b_row.iter())) {
            *d_ij = a_ij + b_ij;
        }
    }
    SpatialInertia { data: d }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    fn approx_eq3(a: [f64; 3], b: [f64; 3]) -> bool {
        approx_eq(a[0], b[0]) && approx_eq(a[1], b[1]) && approx_eq(a[2], b[2])
    }

    // ── SpatialVelocity ───────────────────────────────────────────────────

    #[test]
    fn spatial_velocity_zero() {
        let v = SpatialVelocity::zero();
        assert!(approx_eq3(v.omega, [0.0; 3]));
        assert!(approx_eq3(v.v, [0.0; 3]));
    }

    #[test]
    fn spatial_velocity_as_array() {
        let sv = SpatialVelocity::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let arr = sv.as_array();
        assert_eq!(arr, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    // ── SpatialForce ──────────────────────────────────────────────────────

    #[test]
    fn spatial_force_zero() {
        let f = SpatialForce::zero();
        assert!(approx_eq3(f.tau, [0.0; 3]));
        assert!(approx_eq3(f.f, [0.0; 3]));
    }

    #[test]
    fn spatial_force_as_array() {
        let sf = SpatialForce::new([7.0, 8.0, 9.0], [10.0, 11.0, 12.0]);
        let arr = sf.as_array();
        assert_eq!(arr, [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
    }

    // ── spatial_cross_product ─────────────────────────────────────────────

    #[test]
    fn cross_product_zero_velocities() {
        let v1 = SpatialVelocity::zero();
        let v2 = SpatialVelocity::zero();
        let result = spatial_cross_product(&v1, &v2);
        assert!(approx_eq3(result.omega, [0.0; 3]));
        assert!(approx_eq3(result.v, [0.0; 3]));
    }

    #[test]
    fn cross_product_angular_only() {
        // omega1 = [1,0,0], omega2 = [0,1,0], v1=v2=0
        let v1 = SpatialVelocity::new([1.0, 0.0, 0.0], [0.0; 3]);
        let v2 = SpatialVelocity::new([0.0, 1.0, 0.0], [0.0; 3]);
        let res = spatial_cross_product(&v1, &v2);
        // omega_out = [1,0,0] × [0,1,0] = [0,0,1]
        assert!(approx_eq3(res.omega, [0.0, 0.0, 1.0]));
        assert!(approx_eq3(res.v, [0.0; 3]));
    }

    #[test]
    fn cross_product_antisymmetry() {
        let v1 = SpatialVelocity::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let v2 = SpatialVelocity::new([7.0, 8.0, 9.0], [10.0, 11.0, 12.0]);
        let r12 = spatial_cross_product(&v1, &v2);
        let r21 = spatial_cross_product(&v2, &v1);
        // r12 = -r21 for the omega part (cross product antisymmetry)
        for i in 0..3 {
            assert!(approx_eq(r12.omega[i], -r21.omega[i]));
        }
    }

    #[test]
    fn cross_product_linear_velocity_coupling() {
        let v1 = SpatialVelocity::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let v2 = SpatialVelocity::new([0.0; 3], [0.0, 1.0, 0.0]);
        let res = spatial_cross_product(&v1, &v2);
        // omega_out = [0,0,1] × [0,0,0] = [0,0,0]
        // v_out = [0,0,1]×[0,1,0] + [1,0,0]×[0,0,0]
        //       = [-1,0,0] + [0,0,0] = [-1,0,0]
        assert!(approx_eq3(res.omega, [0.0; 3]));
        assert!(approx_eq3(res.v, [-1.0, 0.0, 0.0]));
    }

    // ── spatial_dot ───────────────────────────────────────────────────────

    #[test]
    fn power_zero_for_zero_velocity() {
        let f = SpatialForce::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let v = SpatialVelocity::zero();
        assert!(approx_eq(spatial_dot(&f, &v), 0.0));
    }

    #[test]
    fn power_zero_for_zero_force() {
        let f = SpatialForce::zero();
        let v = SpatialVelocity::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(approx_eq(spatial_dot(&f, &v), 0.0));
    }

    #[test]
    fn power_torque_times_omega() {
        let f = SpatialForce::new([2.0, 0.0, 0.0], [0.0; 3]);
        let v = SpatialVelocity::new([3.0, 0.0, 0.0], [0.0; 3]);
        // P = tau · omega = 2*3 = 6
        assert!(approx_eq(spatial_dot(&f, &v), 6.0));
    }

    #[test]
    fn power_force_times_velocity() {
        let f = SpatialForce::new([0.0; 3], [5.0, 0.0, 0.0]);
        let v = SpatialVelocity::new([0.0; 3], [4.0, 0.0, 0.0]);
        // P = f · v = 5*4 = 20
        assert!(approx_eq(spatial_dot(&f, &v), 20.0));
    }

    #[test]
    fn power_combined() {
        let f = SpatialForce::new([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]);
        let v = SpatialVelocity::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        // P = (1+1+1)*2 + (1+1+1)*3 = 6 + 9 = 15
        assert!(approx_eq(spatial_dot(&f, &v), 15.0));
    }

    // ── plucker_transform ─────────────────────────────────────────────────

    #[test]
    fn plucker_identity_rotation_no_translation() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let vel = SpatialVelocity::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let transformed = plucker_transform(&vel, &identity, [0.0; 3]);
        assert!(approx_eq3(transformed.omega, vel.omega));
        assert!(approx_eq3(transformed.v, vel.v));
    }

    #[test]
    fn plucker_pure_translation_shifts_linear_velocity() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Pure rotation omega = [0,0,1], v = [1,0,0], translate along x by 2
        let vel = SpatialVelocity::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let r_ab = [2.0, 0.0, 0.0];
        let transformed = plucker_transform(&vel, &identity, r_ab);
        // omega unchanged; v_B = R*(v_A - r_ab × omega_A)
        // r_ab × omega = [2,0,0] × [0,0,1] = [0*1-0*0, 0*0-2*1, 2*0-0*0] = [0,-2,0]
        // v_A - r_ab×omega = [1,0,0] - [0,-2,0] = [1,2,0]
        assert!(approx_eq3(transformed.omega, [0.0, 0.0, 1.0]));
        assert!(approx_eq3(transformed.v, [1.0, 2.0, 0.0]));
    }

    #[test]
    fn plucker_90deg_rotation_around_z() {
        // Rotate 90° around Z: R = [[0,-1,0],[1,0,0],[0,0,1]]
        let rot = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let vel = SpatialVelocity::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let transformed = plucker_transform(&vel, &rot, [0.0; 3]);
        // omega_b = R * [1,0,0] = [0,1,0]
        assert!((transformed.omega[0]).abs() < 1e-10);
        assert!((transformed.omega[1] - 1.0).abs() < 1e-10);
        assert!((transformed.omega[2]).abs() < 1e-10);
    }

    // ── SpatialInertia ────────────────────────────────────────────────────

    #[test]
    fn spatial_inertia_identity() {
        let si = SpatialInertia::identity();
        for i in 0..6 {
            for j in 0..6 {
                if i == j {
                    assert!(approx_eq(si.data[i][j], 1.0));
                } else {
                    assert!(approx_eq(si.data[i][j], 0.0));
                }
            }
        }
    }

    #[test]
    fn spatial_inertia_from_rigid_body_diagonal() {
        let si = SpatialInertia::from_rigid_body(2.0, [1.0, 2.0, 3.0]);
        assert!(approx_eq(si.data[0][0], 1.0));
        assert!(approx_eq(si.data[1][1], 2.0));
        assert!(approx_eq(si.data[2][2], 3.0));
        assert!(approx_eq(si.data[3][3], 2.0));
        assert!(approx_eq(si.data[4][4], 2.0));
        assert!(approx_eq(si.data[5][5], 2.0));
        // Off-diagonals zero
        assert!(approx_eq(si.data[0][3], 0.0));
    }

    #[test]
    fn mul_velocity_identity_inertia() {
        let si = SpatialInertia::identity();
        let vel = SpatialVelocity::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let force = si.mul_velocity(&vel);
        assert!(approx_eq3(force.tau, [1.0, 2.0, 3.0]));
        assert!(approx_eq3(force.f, [4.0, 5.0, 6.0]));
    }

    #[test]
    fn mul_velocity_scaled_inertia() {
        let si = SpatialInertia::from_rigid_body(2.0, [3.0, 3.0, 3.0]);
        let vel = SpatialVelocity::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let force = si.mul_velocity(&vel);
        // tau = I * omega = 3 * [1,0,0] = [3,0,0]
        assert!(approx_eq(force.tau[0], 3.0));
        // f = m * v = 2 * [1,0,0] = [2,0,0]
        assert!(approx_eq(force.f[0], 2.0));
    }

    // ── composite_inertia ─────────────────────────────────────────────────

    #[test]
    fn composite_inertia_sums_correctly() {
        let a = SpatialInertia::from_rigid_body(1.0, [1.0, 1.0, 1.0]);
        let b = SpatialInertia::from_rigid_body(2.0, [3.0, 3.0, 3.0]);
        let c = composite_inertia(&a, &b);
        // Mass block: 1+2=3
        assert!(approx_eq(c.data[3][3], 3.0));
        // Rotational: 1+3=4
        assert!(approx_eq(c.data[0][0], 4.0));
    }

    #[test]
    fn composite_inertia_symmetric() {
        let a = SpatialInertia::identity();
        let b = SpatialInertia::identity();
        let c = composite_inertia(&a, &b);
        // All diagonal entries should be 2
        for i in 0..6 {
            assert!(approx_eq(c.data[i][i], 2.0));
        }
    }
}
