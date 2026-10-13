// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

use oxiphysics_articulated::{
    body::RigidBody,
    centroidal::{centroidal_momentum, compute_cmm, compute_com},
    ik::{forward_kinematics, geometric_jacobian},
    joint::{FreeFloatingJoint, RevoluteJoint, SphericalJoint},
    model::ArticulatedModel,
    spatial::{SpatialInertia, SpatialTransform},
};

#[test]
fn test_cmm_dimensions() {
    let mut model = ArticulatedModel::new([0.0; 3]);
    let inertia = SpatialInertia::from_com(
        1.0,
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    for i in 0..4 {
        let parent = if i == 0 { None } else { Some(i - 1) };
        let body = RigidBody::new(
            format!("link{i}"),
            inertia,
            parent,
            SpatialTransform::IDENTITY,
        );
        model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
    }
    let q = vec![0.0; 4];
    let a_g = compute_cmm(&model, &q);
    assert_eq!(a_g.len(), 6, "A_G must have 6 rows");
    assert_eq!(a_g[0].len(), 4, "A_G must have n columns");
}

#[test]
fn test_cmm_nonzero_for_revolute_chain() {
    // A 2-link revolute arm should have non-zero CMM entries
    let mut model = ArticulatedModel::new([0.0; 3]);
    let i0 = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let b0 = RigidBody::new("b0", i0, None, SpatialTransform::IDENTITY);
    model.add_body(b0, Box::new(RevoluteJoint::new([1.0, 0.0, 0.0])));
    let i1 = SpatialInertia::from_com(
        1.0,
        [0.0, 0.0, -0.5],
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    );
    let b1 = RigidBody::new(
        "b1",
        i1,
        Some(0),
        SpatialTransform::from_translation([0.0, 0.0, -1.0]),
    );
    model.add_body(b1, Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])));
    let q = vec![0.2, -0.3];
    let a_g = compute_cmm(&model, &q);
    // At least some entries should be non-zero
    let any_nonzero = a_g.iter().any(|row| row.iter().any(|&v| v.abs() > 1e-10));
    assert!(
        any_nonzero,
        "CMM should have non-zero entries for a physical system"
    );
}

// ─── A_G·q̇ equals the true centroidal momentum ──────────────────────────────

fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Per-body data needed by the momentum oracle: mass, COM (body frame) and
/// rotational inertia about the COM (body frame).
struct BodyData {
    mass: f64,
    com: [f64; 3],
    i_com: [[f64; 3]; 3],
}

/// Σ_i momentum of body i about the system COM, using world velocities from
/// the geometric Jacobian of each body (validated against finite differences
/// of FK, including through quaternion joints, in `ik::tests`).
fn momentum_oracle(
    model: &ArticulatedModel,
    data: &[BodyData],
    q: &[f64],
    qd: &[f64],
) -> ([f64; 3], [f64; 3]) {
    let n = model.total_dof();
    let fk = forward_kinematics(model, q);
    let (com_g, _) = compute_com(model, q);
    let mut k_g = [0.0; 3];
    let mut l = [0.0; 3];
    for (i, bd) in data.iter().enumerate() {
        let jac = geometric_jacobian(model, q, i);
        let mut twist = [0.0; 6];
        for (r, t) in twist.iter_mut().enumerate() {
            *t = (0..n).map(|c| jac[r * n + c] * qd[c]).sum();
        }
        let omega = [twist[0], twist[1], twist[2]];
        let v_origin = [twist[3], twist[4], twist[5]];
        let r_i = fk[i].rot;
        let c_off = mat3_vec(r_i, bd.com);
        let c_world = [
            fk[i].trans[0] + c_off[0],
            fk[i].trans[1] + c_off[1],
            fk[i].trans[2] + c_off[2],
        ];
        let w_c = cross(omega, c_off);
        let v_com = [
            v_origin[0] + w_c[0],
            v_origin[1] + w_c[1],
            v_origin[2] + w_c[2],
        ];
        let p = [bd.mass * v_com[0], bd.mass * v_com[1], bd.mass * v_com[2]];
        // World inertia about the body COM: R·I_com·Rᵀ·ω = R·(I_com·(Rᵀω)).
        let omega_body = [
            r_i[0][0] * omega[0] + r_i[1][0] * omega[1] + r_i[2][0] * omega[2],
            r_i[0][1] * omega[0] + r_i[1][1] * omega[1] + r_i[2][1] * omega[2],
            r_i[0][2] * omega[0] + r_i[1][2] * omega[1] + r_i[2][2] * omega[2],
        ];
        let spin = mat3_vec(r_i, mat3_vec(bd.i_com, omega_body));
        let arm = [
            c_world[0] - com_g[0],
            c_world[1] - com_g[1],
            c_world[2] - com_g[2],
        ];
        let orbital = cross(arm, p);
        for k in 0..3 {
            k_g[k] += spin[k] + orbital[k];
            l[k] += p[k];
        }
    }
    (k_g, l)
}

fn build_offset_chain() -> (ArticulatedModel, Vec<BodyData>) {
    let data = vec![
        BodyData {
            mass: 2.0,
            com: [0.1, 0.0, -0.4],
            i_com: [[0.20, 0.01, 0.0], [0.01, 0.15, 0.02], [0.0, 0.02, 0.05]],
        },
        BodyData {
            mass: 1.2,
            com: [0.0, 0.15, -0.3],
            i_com: [[0.08, 0.0, 0.01], [0.0, 0.07, 0.0], [0.01, 0.0, 0.03]],
        },
        BodyData {
            mass: 0.7,
            com: [0.2, -0.05, 0.1],
            i_com: [[0.02, 0.0, 0.0], [0.0, 0.03, 0.005], [0.0, 0.005, 0.02]],
        },
    ];
    let mut model = ArticulatedModel::new([0.0, 0.0, -9.81]);
    model.add_body(
        RigidBody::new(
            "link0",
            SpatialInertia::from_com(data[0].mass, data[0].com, data[0].i_com),
            None,
            SpatialTransform::from_translation([0.3, -0.2, 1.0]),
        ),
        Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
    );
    let x_t1 = SpatialTransform::from_translation([0.4, 0.0, -0.5])
        .compose(&SpatialTransform::from_axis_angle([1.0, 0.0, 0.0], 0.3));
    model.add_body(
        RigidBody::new(
            "link1",
            SpatialInertia::from_com(data[1].mass, data[1].com, data[1].i_com),
            Some(0),
            x_t1,
        ),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    model.add_body(
        RigidBody::new(
            "link2",
            SpatialInertia::from_com(data[2].mass, data[2].com, data[2].i_com),
            Some(1),
            SpatialTransform::from_translation([0.0, 0.3, -0.6]),
        ),
        Box::new(SphericalJoint),
    );
    (model, data)
}

#[test]
fn test_cmm_times_qdot_equals_centroidal_momentum() {
    let (model, data) = build_offset_chain();
    assert_eq!(model.total_dof(), 5);
    assert_eq!(model.total_nq(), 6);
    let h = 0.35_f64;
    let q = vec![
        0.6,
        -0.8,
        h.cos(),
        h.sin() * 0.36,
        h.sin() * 0.48,
        h.sin() * 0.8,
    ];
    let qd = vec![1.1, -0.7, 0.4, 0.9, -1.3];
    let a_g = compute_cmm(&model, &q);
    let mut h_g = [0.0; 6];
    for (r, row) in a_g.iter().enumerate() {
        h_g[r] = row.iter().zip(qd.iter()).map(|(a, v)| a * v).sum();
    }
    let (k_expected, l_expected) = momentum_oracle(&model, &data, &q, &qd);
    for k in 0..3 {
        assert!(
            (h_g[k] - k_expected[k]).abs() < 1e-10,
            "angular momentum[{k}]: A_G·q̇={:.12} oracle={:.12}",
            h_g[k],
            k_expected[k]
        );
        assert!(
            (h_g[3 + k] - l_expected[k]).abs() < 1e-10,
            "linear momentum[{k}]: A_G·q̇={:.12} oracle={:.12}",
            h_g[3 + k],
            l_expected[k]
        );
    }
    let hv = centroidal_momentum(&model, &q, &qd);
    for k in 0..3 {
        assert!((hv.angular[k] - h_g[k]).abs() < 1e-14);
        assert!((hv.linear[k] - h_g[3 + k]).abs() < 1e-14);
    }
}

#[test]
fn test_cmm_free_floating_base_linear_rows_are_total_mass() {
    // For a floating base, the linear-momentum block w.r.t. the base linear
    // velocity (body frame) is M_total·R_base, independent of the limbs.
    let mut model = ArticulatedModel::new([0.0, 0.0, -9.81]);
    let i_com = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]];
    model.add_body(
        RigidBody::new(
            "base",
            SpatialInertia::from_com(3.0, [0.05, 0.0, 0.0], i_com),
            None,
            SpatialTransform::IDENTITY,
        ),
        Box::new(FreeFloatingJoint),
    );
    model.add_body(
        RigidBody::new(
            "leg",
            SpatialInertia::from_com(1.5, [0.0, 0.0, -0.4], i_com),
            Some(0),
            SpatialTransform::from_translation([0.0, 0.1, -0.2]),
        ),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    let h = 0.5_f64;
    let q = vec![0.2, 0.0, 0.9, h.cos(), 0.0, 0.0, h.sin(), 0.7];
    let a_g = compute_cmm(&model, &q);
    let r_base = forward_kinematics(&model, &q)[0].rot;
    for row in 0..3 {
        for col in 0..3 {
            assert!(
                (a_g[3 + row][col] - 4.5 * r_base[row][col]).abs() < 1e-12,
                "A_G[{}][{col}]={} expected {}",
                3 + row,
                a_g[3 + row][col],
                4.5 * r_base[row][col]
            );
        }
    }
    // Moving only the leg joint leaves total linear momentum consistent with
    // the oracle as well.
    let data = vec![
        BodyData {
            mass: 3.0,
            com: [0.05, 0.0, 0.0],
            i_com,
        },
        BodyData {
            mass: 1.5,
            com: [0.0, 0.0, -0.4],
            i_com,
        },
    ];
    let qd = vec![0.3, -0.2, 0.1, 0.5, 0.4, -0.6, 1.2];
    let hv = centroidal_momentum(&model, &q, &qd);
    let (k_expected, l_expected) = momentum_oracle(&model, &data, &q, &qd);
    for k in 0..3 {
        assert!((hv.angular[k] - k_expected[k]).abs() < 1e-10);
        assert!((hv.linear[k] - l_expected[k]).abs() < 1e-10);
    }
}
