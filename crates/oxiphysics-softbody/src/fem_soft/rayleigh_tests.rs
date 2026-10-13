// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the stiffness-proportional part of Rayleigh damping.
//!
//! The damping force is `-beta_K * K v` with `K = -d f_el / d x` built from the
//! **elastic** forces only.  The previous implementation scaled the complete
//! force array (gravity included) by `-beta_K`: it was velocity independent,
//! changed the free-fall acceleration to `(1 - beta_K) g` and dissipated nothing.

use super::{
    CorotFemElement4, CorotFemNode, CorotFemTet, FemSoftBodyVerlet, RayleighDamping,
    VelocityVerletIntegrator, assemble_internal_forces, tet_stiffness_matrix,
};

const YOUNG: f64 = 1000.0;
const POISSON: f64 = 0.3;

fn lame(young: f64, poisson: f64) -> (f64, f64) {
    let lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = young / (2.0 * (1.0 + poisson));
    (lambda, mu)
}

fn unit_tet() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// Skewed tetrahedron: neither `Dm` nor `Dm^-1` is symmetric.
fn skewed_tet() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [0.3, 1.1, 0.0],
        [0.2, 0.4, 0.9],
    ]
}

fn verlet_body(beta_k: f64, gravity: [f64; 3]) -> FemSoftBodyVerlet {
    let pos = unit_tet();
    let elem = CorotFemElement4::new([0, 1, 2, 3], &pos, YOUNG, POISSON);
    FemSoftBodyVerlet::new(pos, vec![1.0; 4], vec![elem], gravity, 1.0e-3, 0.0, beta_k)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn max_abs(values: &[[f64; 3]]) -> f64 {
    values
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

#[test]
fn free_fall_acceleration_is_exactly_g_with_stiffness_damping() {
    let g = [0.0, -9.81, 0.0];
    let mut damped = verlet_body(0.1, g);
    let mut undamped = verlet_body(0.0, g);
    for _ in 0..100 {
        damped.step();
        undamped.step();
    }
    // Rigid translation has no strain rate: K v = 0, so beta_K must not matter.
    for i in 0..4 {
        for d in 0..3 {
            assert!(
                (damped.integrator.prev_accel[i][d] - undamped.integrator.prev_accel[i][d]).abs()
                    < 1.0e-8,
                "acceleration of node {i} component {d}"
            );
            assert!(
                (damped.velocities[i][d] - undamped.velocities[i][d]).abs() < 1.0e-8,
                "velocity of node {i} component {d}"
            );
            assert!(
                (damped.positions[i][d] - undamped.positions[i][d]).abs() < 1.0e-8,
                "position of node {i} component {d}"
            );
        }
        assert!(
            (damped.integrator.prev_accel[i][1] + 9.81).abs() < 1.0e-8,
            "free-fall acceleration {} must be g, not (1 - beta_K) g",
            damped.integrator.prev_accel[i][1]
        );
    }
}

#[test]
fn step_sim_free_fall_ignores_stiffness_damping() {
    let rest_nodes: Vec<CorotFemNode> = unit_tet()
        .iter()
        .map(|&position| CorotFemNode {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass: 1.0,
        })
        .collect();
    let tets = vec![CorotFemTet::new([0, 1, 2, 3], &rest_nodes)];
    let (lambda, mu) = lame(YOUNG, POISSON);
    let rayleigh = RayleighDamping::new(0.0, 0.1);
    let pinned = vec![false; 4];
    let run = |damping: Option<&RayleighDamping>| {
        let mut integrator = VelocityVerletIntegrator::new(4, 1.0e-3);
        let mut positions = unit_tet();
        let mut velocities = vec![[0.0_f64; 3]; 4];
        for _ in 0..50 {
            integrator.step_sim(
                &mut positions,
                &mut velocities,
                &[1.0; 4],
                &tets,
                mu,
                lambda,
                [0.0, -9.81, 0.0],
                damping,
                &pinned,
            );
        }
        (positions, velocities, integrator.prev_accel)
    };
    let (x_damped, v_damped, a_damped) = run(Some(&rayleigh));
    let (x_free, v_free, a_free) = run(None);
    for i in 0..4 {
        for d in 0..3 {
            assert!((x_damped[i][d] - x_free[i][d]).abs() < 1.0e-8);
            assert!((v_damped[i][d] - v_free[i][d]).abs() < 1.0e-8);
            assert!((a_damped[i][d] - a_free[i][d]).abs() < 1.0e-8);
        }
    }
}

#[test]
fn stiffness_damping_dissipates_vibration_energy_monotonically() {
    // Nodes 0..2 are fixed, so the whole motion of node 3 is deformation.
    let make = |beta_k: f64| {
        let mut body = verlet_body(beta_k, [0.0; 3]);
        for node in 0..3 {
            body.pin_node(node);
        }
        body.positions[3] = [0.0, 0.0, 1.2];
        body
    };
    let mut damped = make(0.04);
    let mut undamped = make(0.0);
    let e0 = damped.mechanical_energy();
    assert!(e0 > 0.0, "the displaced node stores strain energy");

    let mut samples = vec![e0];
    for _ in 0..10 {
        for _ in 0..100 {
            damped.step();
            undamped.step();
        }
        samples.push(damped.mechanical_energy());
    }
    for (i, pair) in samples.windows(2).enumerate() {
        assert!(
            pair[1] < pair[0],
            "energy must decrease monotonically: sample {i} = {}, sample {} = {}",
            pair[0],
            i + 1,
            pair[1]
        );
    }
    let undamped_energy = undamped.mechanical_energy();
    assert!(
        (undamped_energy - e0).abs() < 0.05 * e0,
        "without damping the integrator keeps the energy: {e0} -> {undamped_energy}"
    );
    let damped_energy = *samples.last().unwrap_or(&e0);
    assert!(
        damped_energy < 0.05 * undamped_energy,
        "beta_K must dissipate: damped {damped_energy}, undamped {undamped_energy}"
    );
}

#[test]
fn stiffness_velocity_product_matches_linear_stiffness_at_rest() {
    let pos = skewed_tet();
    let element = CorotFemElement4::new([0, 1, 2, 3], &pos, YOUNG, POISSON);
    let velocities = vec![
        [0.3, -0.2, 0.1],
        [-0.1, 0.4, 0.2],
        [0.25, 0.05, -0.3],
        [-0.15, -0.35, 0.45],
    ];
    let elements = vec![element];
    let kv = RayleighDamping::stiffness_velocity_product(&pos, &velocities, |x| {
        assemble_internal_forces(x, &elements)
    });

    let (lambda, mu) = lame(YOUNG, POISSON);
    let k = tet_stiffness_matrix(pos[0], pos[1], pos[2], pos[3], lambda, mu);
    let mut expected = [[0.0_f64; 3]; 4];
    for (row, out) in expected.iter_mut().flatten().enumerate() {
        for (col, v) in velocities.iter().flatten().enumerate() {
            *out += k[row * 12 + col] * v;
        }
    }
    let scale = max_abs(&expected);
    assert!(scale > 1.0, "the velocity field must strain the element");
    for i in 0..4 {
        for d in 0..3 {
            assert!(
                (kv[i][d] - expected[i][d]).abs() < 1.0e-6 * scale,
                "node {i} component {d}: finite-difference K v = {} vs linear K v = {}",
                kv[i][d],
                expected[i][d]
            );
        }
    }
}

#[test]
fn stiffness_velocity_product_ignores_rigid_motion() {
    let pos = skewed_tet();
    let elements = vec![CorotFemElement4::new([0, 1, 2, 3], &pos, YOUNG, POISSON)];
    let force_field = |x: &[[f64; 3]]| assemble_internal_forces(x, &elements);

    let deforming = vec![
        [0.3, -0.2, 0.1],
        [-0.1, 0.4, 0.2],
        [0.25, 0.05, -0.3],
        [-0.15, -0.35, 0.45],
    ];
    let reference = max_abs(&RayleighDamping::stiffness_velocity_product(
        &pos,
        &deforming,
        force_field,
    ));
    assert!(reference > 1.0);

    let translation = vec![[0.7, -1.3, 0.4]; 4];
    let kv_translation =
        RayleighDamping::stiffness_velocity_product(&pos, &translation, force_field);
    assert!(max_abs(&kv_translation) < 1.0e-6 * reference);

    let omega = [0.3, -0.2, 0.5];
    let rotation: Vec<[f64; 3]> = pos.iter().map(|&x| cross(omega, x)).collect();
    let kv_rotation = RayleighDamping::stiffness_velocity_product(&pos, &rotation, force_field);
    assert!(
        max_abs(&kv_rotation) < 1.0e-6 * reference,
        "a corotational element must not damp rigid rotation: {}",
        max_abs(&kv_rotation)
    );
}

#[test]
fn stiffness_velocity_product_of_resting_body_is_zero() {
    let pos = unit_tet();
    let elements = vec![CorotFemElement4::new([0, 1, 2, 3], &pos, YOUNG, POISSON)];
    let kv = RayleighDamping::stiffness_velocity_product(&pos, &[[0.0; 3]; 4], |x| {
        assemble_internal_forces(x, &elements)
    });
    assert_eq!(max_abs(&kv), 0.0);
    let empty = RayleighDamping::stiffness_velocity_product(&[], &[], |_| Vec::new());
    assert!(empty.is_empty());
}

#[test]
fn apply_adds_mass_and_stiffness_terms() {
    let rd = RayleighDamping::new(0.5, 0.25);
    let mut forces = vec![[0.0, 0.0, 5.0]];
    let velocities = vec![[1.0, 2.0, 3.0]];
    let masses = vec![2.0];
    let kv = vec![[10.0, 0.0, -4.0]];
    rd.apply(&mut forces, &velocities, &masses, &kv);
    // f = f0 - alpha m v - beta_K K v
    let expected = [
        -0.5 * 2.0 * 1.0 - 0.25 * 10.0,
        -0.5 * 2.0 * 2.0,
        5.0 - 0.5 * 2.0 * 3.0 + 0.25 * 4.0,
    ];
    for d in 0..3 {
        assert!(
            (forces[0][d] - expected[d]).abs() < 1.0e-12,
            "component {d}: {} vs {}",
            forces[0][d],
            expected[d]
        );
    }
    assert!(rd.uses_stiffness_term());
    assert!(!RayleighDamping::new(1.0, 0.0).uses_stiffness_term());
}
