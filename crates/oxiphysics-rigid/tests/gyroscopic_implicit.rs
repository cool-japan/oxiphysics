// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the implicit gyroscopic step (Catto, GDC 2015).

use oxiphysics_rigid::gyroscopic::{Gyroscope, gyroscopic_implicit_step, step_gyroscope};

// ---------------------------------------------------------------------------
// Quaternion helper for tests (scalar-last [x, y, z, w] convention)
// ---------------------------------------------------------------------------

fn quat_integrate_test(q: [f64; 4], omega: [f64; 3], dt: f64) -> [f64; 4] {
    // First-order: dq/dt = 0.5 * q ⊗ Ω where Ω = [ωx, ωy, ωz, 0] scalar-last.
    let [qx, qy, qz, qw] = q;
    let [wx, wy, wz] = omega;
    let h = 0.5 * dt;
    let nqx = qx + (qw * wx + qy * wz - qz * wy) * h;
    let nqy = qy + (qw * wy - qx * wz + qz * wx) * h;
    let nqz = qz + (qw * wz + qx * wy - qy * wx) * h;
    let nqw = qw + (-qx * wx - qy * wy - qz * wz) * h;
    let norm = (nqx * nqx + nqy * nqy + nqz * nqz + nqw * nqw).sqrt();
    if norm < 1e-15 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [nqx / norm, nqy / norm, nqz / norm, nqw / norm]
}

/// Rotate a vector by a quaternion [x,y,z,w] (scalar-last): returns R·v (world = R·body).
fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    let r = [
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
    ];
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

/// Rotate a vector by the TRANSPOSE of the rotation matrix of q: v_body = Rᵀ·v_world.
fn quat_rotate_inverse(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    // Conjugate quaternion: [-qx,-qy,-qz,qw] gives Rᵀ
    let r = [
        [
            1.0 - 2.0 * (qy * qy + qz * qz),
            2.0 * (qx * qy + qz * qw),
            2.0 * (qx * qz - qy * qw),
        ],
        [
            2.0 * (qx * qy - qz * qw),
            1.0 - 2.0 * (qx * qx + qz * qz),
            2.0 * (qy * qz + qx * qw),
        ],
        [
            2.0 * (qx * qz + qy * qw),
            2.0 * (qy * qz - qx * qw),
            1.0 - 2.0 * (qx * qx + qy * qy),
        ],
    ];
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// Test 1 — Energy stability gate (Dzhanibekov / T-handle effect)
// ---------------------------------------------------------------------------

/// Simulate torque-free motion of a body spinning about its intermediate
/// principal axis (Dzhanibekov effect) for 60 seconds using the implicit step.
/// Asserts: no NaN/Inf at any step, and kinetic energy does not blow up.
#[test]
fn test_dzhanibekov_energy_stability() {
    // Distinct principal moments: Ix < Iy < Iz
    let inertia = [1.0_f64, 2.0, 3.0];
    let [ix, iy, iz] = inertia;
    let dt = 1.0_f64 / 60.0;
    let n_steps = 3_600_usize; // 60 seconds

    // Initial ω is almost aligned with the intermediate axis (y), which is
    // unstable — causes the Dzhanibekov flip.
    // We work entirely in world frame; initial orientation = identity.
    let mut omega_world = [0.01_f64, 10.0, 0.01];
    let mut q = [0.0_f64, 0.0, 0.0, 1.0]; // identity

    // Compute initial kinetic energy in body frame.
    let ke_initial = {
        let ob = quat_rotate_inverse(q, omega_world);
        0.5 * (ix * ob[0] * ob[0] + iy * ob[1] * ob[1] + iz * ob[2] * ob[2])
    };
    assert!(ke_initial > 0.0, "initial KE must be positive");

    for step in 0..n_steps {
        omega_world = gyroscopic_implicit_step(omega_world, q, inertia, dt);
        q = quat_integrate_test(q, omega_world, dt);

        // Guard: no NaN or Inf.
        for &c in omega_world.iter() {
            assert!(
                c.is_finite(),
                "ω is non-finite at step {step}: ω={omega_world:?}"
            );
        }
        for &c in q.iter() {
            assert!(
                c.is_finite(),
                "quaternion is non-finite at step {step}: q={q:?}"
            );
        }
    }

    // Kinetic energy must stay bounded (implicit step is dissipative).
    let ke_final = {
        let ob = quat_rotate_inverse(q, omega_world);
        0.5 * (ix * ob[0] * ob[0] + iy * ob[1] * ob[1] + iz * ob[2] * ob[2])
    };
    assert!(
        ke_final < ke_initial * 10.0,
        "kinetic energy blew up: initial={ke_initial:.4}, final={ke_final:.4}"
    );
    assert!(
        ke_final > 0.0,
        "kinetic energy went non-positive: {ke_final}"
    );
}

// ---------------------------------------------------------------------------
// Test 2 — Short-term oracle agreement with RK4
// ---------------------------------------------------------------------------

/// Over a short window (1 second, 60 steps), the implicit step and the RK4
/// oracle (step_gyroscope) should produce orientation trajectories that remain
/// close (|q_implicit · q_rk4| > 0.99).
#[test]
fn test_oracle_agreement_with_rk4() {
    let inertia = [1.0_f64, 2.0, 3.0];
    let dt = 1.0_f64 / 60.0;
    let n_steps = 60_usize;

    // Stable initial condition: spin mostly about the major axis (Iz=3 largest).
    // Both representations start at identity orientation.
    let omega_init = [0.1_f64, 0.1, 5.0];

    // --- RK4 path (Gyroscope struct, body frame) ---
    let mut gyro = Gyroscope::new([0.0, 0.0, 1.0], 5.0, inertia);
    gyro.angular_velocity = omega_init; // body frame = world frame at identity

    // --- Implicit path (world frame) ---
    let mut q_impl = [0.0_f64, 0.0, 0.0, 1.0];
    let mut omega_impl = omega_init;

    for _ in 0..n_steps {
        step_gyroscope(&mut gyro, [0.0, 0.0, 0.0], dt);

        omega_impl = gyroscopic_implicit_step(omega_impl, q_impl, inertia, dt);
        q_impl = quat_integrate_test(q_impl, omega_impl, dt);
    }

    // Compare orientations via quaternion dot product.
    let q_rk4 = gyro.orientation; // [x, y, z, w] scalar-last
    let dot =
        (q_rk4[0] * q_impl[0] + q_rk4[1] * q_impl[1] + q_rk4[2] * q_impl[2] + q_rk4[3] * q_impl[3])
            .abs(); // absolute value handles double-cover
    assert!(
        dot > 0.99,
        "orientation diverged between RK4 and implicit: |q·q_rk4|={dot:.6}"
    );
}

// ---------------------------------------------------------------------------
// Test 3 — Off-parity (near-zero ω returns near-zero correction)
// ---------------------------------------------------------------------------

/// When angular velocity is effectively zero, the implicit step is a near-identity
/// transformation (no gyroscopic coupling possible with zero angular momentum).
#[test]
fn test_off_parity_zero_omega() {
    let inertia = [1.0_f64, 2.0, 3.0];
    let dt = 1.0_f64 / 60.0;
    let q = [0.0_f64, 0.0, 0.0, 1.0]; // identity orientation
    let omega_zero = [0.0_f64, 0.0, 0.0];

    let result = gyroscopic_implicit_step(omega_zero, q, inertia, dt);

    // With zero ω, f(0) = 0, so Newton step gives Δω = 0.
    for (i, &c) in result.iter().enumerate() {
        assert!(
            c.abs() < 1e-12,
            "expected near-zero result for zero ω, got result[{i}]={c}"
        );
    }
}

/// For a very small ω (far from any spinning-top regime), the correction
/// should be tiny (second-order in ω).
#[test]
fn test_off_parity_tiny_omega() {
    let inertia = [1.0_f64, 2.0, 3.0];
    let dt = 1.0_f64 / 60.0;
    let q = [0.0_f64, 0.0, 0.0, 1.0]; // identity orientation
    let omega_tiny = [1e-6_f64, 2e-6, 3e-6];

    let result = gyroscopic_implicit_step(omega_tiny, q, inertia, dt);

    // Result should be very close to input (gyroscopic coupling is second-order).
    for i in 0..3 {
        let diff = (result[i] - omega_tiny[i]).abs();
        assert!(
            diff < 1e-10,
            "large correction for tiny ω: diff[{i}]={diff:.2e}"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 4 — Jacobian unit test (analytical vs finite difference)
// ---------------------------------------------------------------------------

/// Verify the 3×3 Jacobian J = I/dt + [ω]×·I − [I·ω]× against finite
/// differences of the residual f(ω₁) at a fixed evaluation point.
#[test]
fn test_jacobian_finite_difference() {
    let [ix, iy, iz] = [1.0_f64, 2.0, 3.0];
    let dt = 0.016_f64;
    let eps = 1e-6_f64;

    // Evaluation point ω₁ (current Newton iterate).
    let omega = [2.0_f64, 3.0, 1.5];
    // Reference point ω₀ (previous step angular velocity — appears only in f, not J).
    let omega0 = [1.5_f64, 2.0, 1.0];

    // Residual function (body frame, diagonal I).
    let residual = |w: [f64; 3]| -> [f64; 3] {
        let [wx, wy, wz] = w;
        let lx = ix * wx;
        let ly = iy * wy;
        let lz = iz * wz;
        [
            ix * (wx - omega0[0]) / dt + wy * lz - wz * ly,
            iy * (wy - omega0[1]) / dt + wz * lx - wx * lz,
            iz * (wz - omega0[2]) / dt + wx * ly - wy * lx,
        ]
    };

    // Analytical Jacobian at `omega`.
    // J[i][k] = ∂f_i/∂ω_k
    let [wx, wy, wz] = omega;
    let j_analytical = [
        [ix / dt, wz * (iz - iy), wy * (iz - iy)],
        [wz * (ix - iz), iy / dt, wx * (ix - iz)],
        [wy * (iy - ix), wx * (iy - ix), iz / dt],
    ];

    let f_base = residual(omega);

    for k in 0..3 {
        let mut omega_perturbed = omega;
        omega_perturbed[k] += eps;
        let f_perturbed = residual(omega_perturbed);

        for i in 0..3 {
            let fd = (f_perturbed[i] - f_base[i]) / eps;
            let an = j_analytical[i][k];
            let diff = (fd - an).abs();
            assert!(
                diff < 1e-4,
                "J[{i}][{k}]: analytical={an:.6}, finite-diff={fd:.6}, diff={diff:.2e}"
            );
        }
    }
}
