// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`CorotationalElementRaw`].
//!
//! The previous element polar-decomposed the raw current edge matrix instead of
//! the deformation gradient `F = Ds Dm^-1` and applied a "force proportional to
//! the displacement from the rotated rest shape" with a stiffness of
//! `E * V0`.  On any rest shape other than the unit corner tetrahedron the
//! rotation was not the identity at rest, so a resting element carried internal
//! forces.  Poisson's ratio and the damping coefficient were never read.

use super::{CorotationalElementRaw, FemSoftBodyRaw, SoftMaterial};

type Vec3 = [f64; 3];

const YOUNG: f64 = 2.0e5;
const POISSON: f64 = 0.35;

fn material(poisson: f64, damping: f64) -> SoftMaterial {
    SoftMaterial {
        youngs_modulus: YOUNG,
        poisson_ratio: poisson,
        damping,
    }
}

/// Skewed rest tetrahedron.
fn skewed_rest() -> Vec<Vec3> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.1, 0.0],
        [0.2, 0.4, 0.9],
    ]
}

fn element(poisson: f64, damping: f64) -> CorotationalElementRaw {
    CorotationalElementRaw::new([0, 1, 2, 3], &skewed_rest(), material(poisson, damping))
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

/// Independent shape-function gradients from the reciprocal basis of the edges.
fn volume_and_gradients(p: &[Vec3]) -> (f64, [Vec3; 4]) {
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

fn max_abs(forces: &[Vec3; 4]) -> f64 {
    forces
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

/// Linear-elastic nodal forces `-V sigma(a) grad N_k` of a uniform symmetric
/// strain (or strain-rate) tensor `a`.
fn linear_forces(rest: &[Vec3], a: [[f64; 3]; 3], lambda: f64, mu: f64, scale: f64) -> [Vec3; 4] {
    let (volume, grads) = volume_and_gradients(rest);
    let trace = a[0][0] + a[1][1] + a[2][2];
    let mut sigma = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            sigma[i][j] = 2.0 * mu * a[i][j];
        }
        sigma[i][i] += lambda * trace;
    }
    grads.map(|g| {
        let t = mat_vec(sigma, g);
        [
            -scale * volume * t[0],
            -scale * volume * t[1],
            -scale * volume * t[2],
        ]
    })
}

const STRAIN: [[f64; 3]; 3] = [
    [1.0e-5, 4.0e-6, -2.0e-6],
    [4.0e-6, -6.0e-6, 3.0e-6],
    [-2.0e-6, 3.0e-6, 8.0e-6],
];

#[test]
fn skewed_element_is_force_free_at_rest() {
    let rest = skewed_rest();
    let forces = element(POISSON, 0.0).compute_forces(&rest);
    for (node, f) in forces.iter().enumerate() {
        for (d, v) in f.iter().enumerate() {
            assert!(
                v.abs() < 1.0e-9,
                "node {node} component {d} of the resting element carries force {v:e}"
            );
        }
    }
}

#[test]
fn skewed_element_is_force_free_under_finite_rigid_motion() {
    let q = rotation_matrix([1.0, 2.0, 2.0], 1.1);
    let moved: Vec<Vec3> = skewed_rest()
        .iter()
        .map(|&x| {
            let r = mat_vec(q, x);
            [r[0] + 0.7, r[1] - 1.3, r[2] + 0.2]
        })
        .collect();
    let forces = element(POISSON, 0.0).compute_forces(&moved);
    for (node, f) in forces.iter().enumerate() {
        for (d, v) in f.iter().enumerate() {
            assert!(
                v.abs() < 1.0e-8,
                "node {node} component {d} of the rigidly moved element carries force {v:e}"
            );
        }
    }
}

#[test]
fn rotation_is_the_polar_factor_of_the_deformation_gradient() {
    let rest = skewed_rest();
    let elem = element(POISSON, 0.0);
    let at_rest = elem.compute_rotation(&rest);
    for (i, row) in at_rest.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let want = if i == j { 1.0 } else { 0.0 };
            assert!(
                (v - want).abs() < 1.0e-9,
                "R[{i}][{j}] = {v} at rest (raw-edge decomposition gives a non-identity matrix)"
            );
        }
    }
    let q = rotation_matrix([1.0, 2.0, 2.0], 1.1);
    let moved: Vec<Vec3> = rest.iter().map(|&x| mat_vec(q, x)).collect();
    let r = elem.compute_rotation(&moved);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (r[i][j] - q[i][j]).abs() < 1.0e-9,
                "R[{i}][{j}] = {} expected {}",
                r[i][j],
                q[i][j]
            );
        }
    }
}

#[test]
fn stretch_gives_restoring_force_without_net_force_or_torque() {
    let rest = skewed_rest();
    let elem = element(POISSON, 0.0);
    let mut stretched = rest.clone();
    stretched[1][0] += 0.2;
    let forces = elem.compute_forces(&stretched);
    let pull = sub(stretched[1], rest[1]);
    assert!(
        dot(forces[1], pull) < 0.0,
        "the force on the displaced node must oppose its displacement"
    );
    let scale = max_abs(&forces);
    assert!(scale > 1.0, "a stretched element must carry force");
    let mut net = [0.0_f64; 3];
    let mut torque = [0.0_f64; 3];
    for (x, f) in stretched.iter().zip(forces.iter()) {
        let moment = cross(*x, *f);
        for d in 0..3 {
            net[d] += f[d];
            torque[d] += moment[d];
        }
    }
    for d in 0..3 {
        assert!(net[d].abs() < 1.0e-9 * scale, "net force {net:?}");
        assert!(torque[d].abs() < 1.0e-9 * scale, "net torque {torque:?}");
    }
}

#[test]
fn small_uniform_strain_reproduces_linear_elasticity() {
    let rest = skewed_rest();
    let elem = element(POISSON, 0.0);
    let (lambda, mu) = elem.material.lame_parameters();
    let deformed: Vec<Vec3> = rest
        .iter()
        .map(|&x| {
            let ax = mat_vec(STRAIN, x);
            [x[0] + ax[0], x[1] + ax[1], x[2] + ax[2]]
        })
        .collect();
    let forces = elem.compute_forces(&deformed);
    let expected = linear_forces(&rest, STRAIN, lambda, mu, 1.0);
    let scale = max_abs(&expected);
    assert!(scale > 1.0e-3);
    for k in 0..4 {
        for d in 0..3 {
            assert!(
                (forces[k][d] - expected[k][d]).abs() < 1.0e-3 * scale,
                "node {k} component {d}: {} vs linear {}",
                forces[k][d],
                expected[k][d]
            );
        }
    }
}

#[test]
fn poisson_ratio_and_youngs_modulus_enter_the_response() {
    let rest = skewed_rest();
    let mut stretched = rest.clone();
    stretched[1][0] += 0.1;
    let incompressible_like = element(0.45, 0.0).compute_forces(&stretched);
    let auxetic_free = element(0.0, 0.0).compute_forces(&stretched);
    let difference: f64 = (0..4)
        .flat_map(|k| (0..3).map(move |d| (k, d)))
        .map(|(k, d)| (incompressible_like[k][d] - auxetic_free[k][d]).abs())
        .sum();
    assert!(
        difference > 0.1 * max_abs(&auxetic_free),
        "nu = 0.45 and nu = 0 must give different forces (difference {difference})"
    );

    let mut doubled = element(POISSON, 0.0);
    doubled.material.youngs_modulus *= 2.0;
    let base = element(POISSON, 0.0).compute_forces(&stretched);
    let twice = doubled.compute_forces(&stretched);
    for k in 0..4 {
        for d in 0..3 {
            assert!(
                (twice[k][d] - 2.0 * base[k][d]).abs() < 1.0e-9 * max_abs(&base),
                "forces are linear in E"
            );
        }
    }
}

#[test]
fn lame_parameters_follow_youngs_modulus_and_poisson_ratio() {
    let m = SoftMaterial {
        youngs_modulus: 1000.0,
        poisson_ratio: 0.25,
        damping: 0.0,
    };
    let (lambda, mu) = m.lame_parameters();
    assert!((lambda - 400.0).abs() < 1.0e-9, "lambda = {lambda}");
    assert!((mu - 400.0).abs() < 1.0e-9, "mu = {mu}");
    let incompressible = SoftMaterial {
        poisson_ratio: 0.5,
        ..m
    };
    let (lambda_inc, mu_inc) = incompressible.lame_parameters();
    assert!(lambda_inc.is_finite() && lambda_inc > 0.0);
    assert!((mu_inc - 1000.0 / 3.0).abs() < 1.0);
}

fn deformed_configuration() -> Vec<Vec3> {
    let q = rotation_matrix([1.0, -1.0, 0.5], 0.6);
    skewed_rest()
        .iter()
        .map(|&x| {
            let stretched = [
                x[0] * 1.06 + 0.03 * x[1],
                x[1] * 0.95,
                x[2] * 1.03 - 0.02 * x[0],
            ];
            mat_vec(q, stretched)
        })
        .collect()
}

#[test]
fn damping_forces_vanish_for_rigid_motion_of_a_deformed_element() {
    let elem = element(POISSON, 0.05);
    let positions = deformed_configuration();
    let omega = [0.4, -0.7, 0.3];
    let translation = [0.2, 0.5, -0.9];
    let rigid: Vec<Vec3> = positions
        .iter()
        .map(|&x| {
            let w = cross(omega, x);
            [
                w[0] + translation[0],
                w[1] + translation[1],
                w[2] + translation[2],
            ]
        })
        .collect();
    let deforming = vec![
        [0.3, -0.2, 0.1],
        [-0.1, 0.4, 0.2],
        [0.25, 0.05, -0.3],
        [-0.15, -0.35, 0.45],
    ];
    let reference = max_abs(&elem.damping_forces(&positions, &deforming));
    assert!(reference > 1.0, "a strain rate must produce damping forces");
    let forces = elem.damping_forces(&positions, &rigid);
    assert!(
        max_abs(&forces) < 1.0e-9 * reference,
        "rigid motion must not be damped: {:e} vs {reference:e}",
        max_abs(&forces)
    );
}

#[test]
fn damping_forces_equal_the_viscous_stress_of_the_strain_rate() {
    let rest = skewed_rest();
    let eta = 0.05;
    let elem = element(POISSON, eta);
    let (lambda, mu) = elem.material.lame_parameters();
    let rate = [[0.3, 0.1, -0.2], [0.1, -0.4, 0.15], [-0.2, 0.15, 0.25]];
    let velocities: Vec<Vec3> = rest.iter().map(|&x| mat_vec(rate, x)).collect();
    let forces = elem.damping_forces(&rest, &velocities);
    let expected = linear_forces(&rest, rate, lambda, mu, eta);
    let scale = max_abs(&expected);
    assert!(scale > 1.0);
    for k in 0..4 {
        for d in 0..3 {
            assert!(
                (forces[k][d] - expected[k][d]).abs() < 1.0e-9 * scale,
                "node {k} component {d}: {} vs {}",
                forces[k][d],
                expected[k][d]
            );
        }
    }
}

#[test]
fn damping_forces_dissipate_power_and_need_a_coefficient() {
    let positions = deformed_configuration();
    let velocities = vec![
        [0.3, -0.2, 0.1],
        [-0.1, 0.4, 0.2],
        [0.25, 0.05, -0.3],
        [-0.15, -0.35, 0.45],
    ];
    let damped = element(POISSON, 0.05);
    let forces = damped.damping_forces(&positions, &velocities);
    let power: f64 = forces
        .iter()
        .zip(velocities.iter())
        .map(|(f, v)| dot(*f, *v))
        .sum();
    assert!(power < -1.0e-6, "damping must remove energy: P = {power}");
    let undamped = element(POISSON, 0.0);
    assert_eq!(
        max_abs(&undamped.damping_forces(&positions, &velocities)),
        0.0
    );
}

#[test]
fn stretched_element_accelerates_back_in_one_step() {
    let mut nodes = skewed_rest();
    nodes[1][0] += 0.3;
    let elem = element(POISSON, 0.01);
    let mut body = FemSoftBodyRaw::new(nodes.clone(), vec![elem], 1.0);
    body.step(1.0e-4, [0.0; 3]);
    assert!(
        body.velocities[1][0] < 0.0,
        "the stretched node must be pulled back: v = {:?}",
        body.velocities[1]
    );
}

#[test]
fn free_flight_of_a_skewed_element_is_not_disturbed() {
    let rest = skewed_rest();
    let elem = element(POISSON, 0.05);
    let mut body = FemSoftBodyRaw::new(rest.clone(), vec![elem], 1.0);
    let v0 = [1.0, -2.0, 0.5];
    for v in body.velocities.iter_mut() {
        *v = v0;
    }
    for _ in 0..100 {
        body.step(0.01, [0.0; 3]);
    }
    for (i, x0) in rest.iter().enumerate() {
        for d in 0..3 {
            assert!(
                (body.velocities[i][d] - v0[d]).abs() < 1.0e-8,
                "velocity of node {i} changed in flight"
            );
            assert!(
                (body.nodes[i][d] - (x0[d] + v0[d])).abs() < 1.0e-6,
                "node {i} component {d} left its ballistic path"
            );
        }
    }
}
