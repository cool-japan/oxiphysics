// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the shape-function gradients behind
//! [`tet_stiffness_matrix`] and [`corot_internal_forces`].
//!
//! Every earlier test used the unit corner tetrahedron, whose `Dm^-1` is the
//! identity: reading the gradients from the columns instead of the rows of
//! `Dm^-1` is invisible there.  The tests below use a skewed tetrahedron and
//! compare against quantities derived independently (reciprocal basis from
//! cross products, analytic linear elasticity of a uniform strain field).

use super::{assemble_global_stiffness, corot_internal_forces, tet_stiffness_matrix};

const LAMBDA: f64 = 1.0e5;
const MU: f64 = 7.0e4;

type Vec3 = [f64; 3];

/// A skewed tetrahedron: `Dm^-1` is neither symmetric nor diagonal.
fn skewed_tet() -> [Vec3; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.1, 0.0],
        [0.2, 0.4, 0.9],
    ]
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mat_vec(m: [[f64; 3]; 3], v: Vec3) -> Vec3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// Volume and the four shape-function gradients, from the reciprocal basis
/// `b_k = (e_{k+1} x e_{k+2}) / (e_1 . e_2 x e_3)` of the edge vectors.
fn volume_and_gradients(p: &[Vec3; 4]) -> (f64, [Vec3; 4]) {
    let e1 = sub(p[1], p[0]);
    let e2 = sub(p[2], p[0]);
    let e3 = sub(p[3], p[0]);
    let triple = dot(e1, cross(e2, e3));
    let scale = |v: Vec3| [v[0] / triple, v[1] / triple, v[2] / triple];
    let g1 = scale(cross(e2, e3));
    let g2 = scale(cross(e3, e1));
    let g3 = scale(cross(e1, e2));
    let g0 = [
        -(g1[0] + g2[0] + g3[0]),
        -(g1[1] + g2[1] + g3[1]),
        -(g1[2] + g2[2] + g3[2]),
    ];
    (triple.abs() / 6.0, [g0, g1, g2, g3])
}

fn stiffness(p: &[Vec3; 4]) -> [f64; 144] {
    tet_stiffness_matrix(p[0], p[1], p[2], p[3], LAMBDA, MU)
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

fn apply(k: &[f64; 144], u: &[f64; 12]) -> [f64; 12] {
    let mut out = [0.0_f64; 12];
    for (i, out_i) in out.iter_mut().enumerate() {
        for (j, u_j) in u.iter().enumerate() {
            *out_i += k[i * 12 + j] * u_j;
        }
    }
    out
}

/// Displacements `u_a = f(p_a)` flattened to the 12-vector layout of `K`.
fn displacement_field(p: &[Vec3; 4], field: impl Fn(Vec3) -> Vec3) -> [f64; 12] {
    let mut u = [0.0_f64; 12];
    for (a, p_a) in p.iter().enumerate() {
        let u_a = field(*p_a);
        u[a * 3..a * 3 + 3].copy_from_slice(&u_a);
    }
    u
}

#[test]
fn stiffness_has_no_strain_energy_for_infinitesimal_rotation() {
    let p = skewed_tet();
    let k = stiffness(&p);
    let omega = [0.3, -0.2, 0.5];
    let u = displacement_field(&p, |x| cross(omega, x));
    let ku = apply(&k, &u);
    let bound = 1.0e-8 * max_abs(&k);
    assert!(
        max_abs(&ku) < bound,
        "an infinitesimal rigid rotation must be in the null space of K: |K u| = {:e}, bound {bound:e}",
        max_abs(&ku)
    );
}

#[test]
fn stiffness_has_no_strain_energy_for_translation() {
    let p = skewed_tet();
    let k = stiffness(&p);
    let u = displacement_field(&p, |_| [0.7, -1.3, 0.4]);
    let ku = apply(&k, &u);
    assert!(max_abs(&ku) < 1.0e-8 * max_abs(&k));
}

#[test]
fn stiffness_energy_matches_linear_elasticity_of_uniform_strain() {
    let p = skewed_tet();
    let k = stiffness(&p);
    let (volume, _) = volume_and_gradients(&p);
    // Symmetric strain tensor with all entries non-zero.
    let a = [
        [0.010, 0.004, -0.002],
        [0.004, -0.006, 0.003],
        [-0.002, 0.003, 0.008],
    ];
    let u = displacement_field(&p, |x| mat_vec(a, x));
    let ku = apply(&k, &u);
    let energy_twice: f64 = u.iter().zip(ku.iter()).map(|(ui, kui)| ui * kui).sum();

    let trace = a[0][0] + a[1][1] + a[2][2];
    let mut a_colon_a = 0.0;
    for row in &a {
        for v in row {
            a_colon_a += v * v;
        }
    }
    // Strain energy of a uniform strain: V (lambda/2 tr(eps)^2 + mu eps:eps).
    let expected_twice = volume * (LAMBDA * trace * trace + 2.0 * MU * a_colon_a);
    assert!(
        (energy_twice - expected_twice).abs() < 1.0e-9 * expected_twice.abs(),
        "u^T K u = {energy_twice}, analytic {expected_twice}"
    );
}

#[test]
fn stiffness_nodal_forces_equal_stress_times_shape_gradient() {
    let p = skewed_tet();
    let k = stiffness(&p);
    let (volume, grads) = volume_and_gradients(&p);
    let a = [
        [0.010, 0.004, -0.002],
        [0.004, -0.006, 0.003],
        [-0.002, 0.003, 0.008],
    ];
    let trace = a[0][0] + a[1][1] + a[2][2];
    let mut sigma = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            sigma[i][j] = 2.0 * MU * a[i][j];
        }
        sigma[i][i] += LAMBDA * trace;
    }
    let u = displacement_field(&p, |x| mat_vec(a, x));
    let ku = apply(&k, &u);
    let scale = max_abs(&ku);
    for (node, grad) in grads.iter().enumerate() {
        let expected = mat_vec(sigma, *grad);
        for d in 0..3 {
            let got = ku[node * 3 + d];
            let want = volume * expected[d];
            assert!(
                (got - want).abs() < 1.0e-9 * scale,
                "node {node} component {d}: K u = {got}, V sigma grad N = {want}"
            );
        }
    }
}

#[test]
fn stiffness_of_unit_corner_tet_is_unchanged() {
    // For Dm = I the row and column readings coincide: the corrected code must
    // reproduce the classic unit-tetrahedron stiffness.
    let p = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let k = stiffness(&p);
    let volume = 1.0 / 6.0;
    // grad N_0 = (-1, -1, -1), so the x column of node 0 in B holds -1 in the
    // exx, exz and exy rows: K[0][0] = V ((lambda + 2 mu) + mu + mu).
    let expected_xx = volume * ((LAMBDA + 2.0 * MU) + MU + MU);
    assert!(
        (k[0] - expected_xx).abs() < 1.0e-9 * expected_xx,
        "K[0][0] = {}, expected {expected_xx}",
        k[0]
    );
}

#[test]
fn assembled_stiffness_has_rigid_rotation_null_space_on_skewed_mesh() {
    // Two skewed tetrahedra sharing a face.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.1, 0.0],
        [0.2, 0.4, 0.9],
        [1.1, 1.0, 1.2],
    ];
    let elements = vec![
        ([0usize, 1, 2, 3], LAMBDA, MU),
        ([1usize, 2, 3, 4], LAMBDA, MU),
    ];
    let k = assemble_global_stiffness(&positions, &elements);
    let dof = 3 * positions.len();
    let omega = [-0.4, 0.25, 0.35];
    let mut u = vec![0.0_f64; dof];
    for (i, x) in positions.iter().enumerate() {
        let u_i = cross(omega, *x);
        u[i * 3..i * 3 + 3].copy_from_slice(&u_i);
    }
    let mut worst = 0.0_f64;
    for row in 0..dof {
        let mut s = 0.0;
        for col in 0..dof {
            s += k[row * dof + col] * u[col];
        }
        worst = worst.max(s.abs());
    }
    assert!(
        worst < 1.0e-8 * max_abs(&k),
        "assembled K u for a rigid rotation = {worst:e}"
    );
}

fn rotation_matrix(axis: Vec3, angle: f64) -> [[f64; 3]; 3] {
    let norm = dot(axis, axis).sqrt();
    let n = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [
            c + n[0] * n[0] * t,
            n[0] * n[1] * t - n[2] * s,
            n[0] * n[2] * t + n[1] * s,
        ],
        [
            n[1] * n[0] * t + n[2] * s,
            c + n[1] * n[1] * t,
            n[1] * n[2] * t - n[0] * s,
        ],
        [
            n[2] * n[0] * t - n[1] * s,
            n[2] * n[1] * t + n[0] * s,
            c + n[2] * n[2] * t,
        ],
    ]
}

#[test]
fn corot_internal_forces_vanish_for_finite_rigid_rotation_of_skewed_tet() {
    let rest = skewed_tet();
    let q = rotation_matrix([1.0, 2.0, 2.0], 0.9);
    let deformed = rest.map(|x| {
        let r = mat_vec(q, x);
        [r[0] + 0.3, r[1] - 0.2, r[2] + 0.1]
    });
    let forces = corot_internal_forces(rest, deformed, LAMBDA, MU);
    let bound = 1.0e-9 * max_abs(&stiffness(&rest));
    for (a, f) in forces.iter().enumerate() {
        for (d, v) in f.iter().enumerate() {
            assert!(
                v.abs() < bound,
                "rigid motion must be force free: node {a} component {d} = {v:e} (bound {bound:e})"
            );
        }
    }
}

#[test]
fn corot_internal_forces_rotate_with_the_element() {
    // x_a = Q (I + A) X_a with A symmetric: polar(F) = Q exactly, so the
    // un-rotated displacement is A X_a and the force is Q (V sigma grad N_a).
    let rest = skewed_tet();
    let q = rotation_matrix([1.0, -1.0, 0.5], 0.6);
    let a = [
        [0.012, -0.003, 0.002],
        [-0.003, 0.007, 0.004],
        [0.002, 0.004, -0.009],
    ];
    let deformed = rest.map(|x| {
        let ax = mat_vec(a, x);
        mat_vec(q, [x[0] + ax[0], x[1] + ax[1], x[2] + ax[2]])
    });
    let forces = corot_internal_forces(rest, deformed, LAMBDA, MU);

    let (volume, grads) = volume_and_gradients(&rest);
    let trace = a[0][0] + a[1][1] + a[2][2];
    let mut sigma = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            sigma[i][j] = 2.0 * MU * a[i][j];
        }
        sigma[i][i] += LAMBDA * trace;
    }
    let expected: Vec<Vec3> = grads
        .iter()
        .map(|grad| {
            let local = mat_vec(sigma, *grad);
            mat_vec(q, [volume * local[0], volume * local[1], volume * local[2]])
        })
        .collect();
    let magnitude = expected.iter().fold(0.0_f64, |acc, f| acc.max(max_abs(f)));
    assert!(magnitude > 1.0, "the strained element must carry force");
    for (node, want) in expected.iter().enumerate() {
        for d in 0..3 {
            assert!(
                (forces[node][d] - want[d]).abs() < 1.0e-8 * magnitude,
                "node {node} component {d}: force {} expected {}",
                forces[node][d],
                want[d]
            );
        }
    }
}
