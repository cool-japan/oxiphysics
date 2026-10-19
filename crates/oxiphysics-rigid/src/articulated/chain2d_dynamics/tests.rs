// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the dynamics of [`ArticulatedChain2D`].
//!
//! The old model measured the gravity lever arm of link `i` without its own
//! length and with the wrong sign, summed the angles over all earlier links
//! ignoring `parent`, treated every joint as an independent rotor
//! (`alpha = tau / I`) and omitted the translation of the links from the
//! kinetic energy.
//!
//! The oracles do not rely on the code under test: the public planar-chain
//! functions of `articulated::functions` (a separate implementation) for serial
//! chains, and finite differences of `forward_kinematics` for arbitrary trees.

use super::*;
use crate::articulated::{
    aba_revolute_chain, coriolis_bias_forces, gravity_compensation_torques,
    joint_space_inertia_matrix_planar, kinetic_energy_planar, rnea_3d_revolute_z,
};
use std::f64::consts::{FRAC_PI_2, PI};

const G: f64 = 9.81;

/// Deterministic generator for the randomised comparisons.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn vector(&mut self, n: usize, lo: f64, hi: f64) -> Vec<f64> {
        (0..n).map(|_| self.range(lo, hi)).collect()
    }
}

/// One link of a test chain.
#[derive(Clone)]
struct Spec {
    mass: f64,
    length: f64,
    /// Moment of inertia about the centre of mass.
    inertia_com: f64,
    parent: Option<usize>,
}

/// Random links: a serial chain, or a random tree whose links attach to a random
/// earlier link (and sometimes start a new root at the base).
fn random_specs(rng: &mut Lcg, n: usize, tree: bool) -> Vec<Spec> {
    (0..n)
        .map(|i| {
            let parent = if i == 0 {
                None
            } else if tree {
                if rng.unit() < 0.15 {
                    None
                } else {
                    Some(((rng.unit() * i as f64) as usize).min(i - 1))
                }
            } else {
                Some(i - 1)
            };
            let mass = rng.range(0.2, 3.0);
            let length = rng.range(0.3, 1.5);
            Spec {
                mass,
                length,
                inertia_com: mass * length * length / 12.0 * rng.range(0.0, 2.0),
                parent,
            }
        })
        .collect()
}

fn build(specs: &[Spec], q: &[f64], qd: &[f64]) -> ArticulatedChain2D {
    let mut chain = ArticulatedChain2D::new();
    for (i, spec) in specs.iter().enumerate() {
        chain.add_link(spec.mass, spec.length, spec.parent);
        let link = &mut chain.links[i];
        // `inertia` is the moment of inertia about the joint (proximal end).
        link.body.inertia = spec.inertia_com + spec.mass * 0.25 * spec.length * spec.length;
        link.joint_angle = q[i];
        link.joint_velocity = qd[i];
    }
    chain
}

fn masses(specs: &[Spec]) -> Vec<f64> {
    specs.iter().map(|s| s.mass).collect()
}

fn lengths(specs: &[Spec]) -> Vec<f64> {
    specs.iter().map(|s| s.length).collect()
}

fn inertias_com(specs: &[Spec]) -> Vec<f64> {
    specs.iter().map(|s| s.inertia_com).collect()
}

fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - e).abs() <= tolerance * (1.0 + e.abs()),
            "{what}[{i}]: {a} vs {e}"
        );
    }
}

// ── independent kinematic oracles ───────────────────────────────────────────

/// Centres of mass and absolute angles from `forward_kinematics`.
fn kinematics(chain: &ArticulatedChain2D) -> (Vec<[f64; 2]>, Vec<f64>) {
    let fk = chain.forward_kinematics([0.0, 0.0], 0.0);
    let mut coms = Vec::new();
    let mut angles = Vec::new();
    for (i, link) in chain.links.iter().enumerate() {
        let origin = match link.parent {
            Some(p) => fk[p].0,
            None => [0.0, 0.0],
        };
        let tip = fk[i].0;
        coms.push([0.5 * (origin[0] + tip[0]), 0.5 * (origin[1] + tip[1])]);
        angles.push(fk[i].1);
    }
    (coms, angles)
}

/// Central-difference derivative of the centres of mass and absolute angles
/// with respect to joint `k`.
fn pose_derivative(chain: &ArticulatedChain2D, k: usize) -> (Vec<[f64; 2]>, Vec<f64>) {
    let h = 1e-6;
    let mut plus = chain.clone();
    plus.links[k].joint_angle += h;
    let mut minus = chain.clone();
    minus.links[k].joint_angle -= h;
    let (com_plus, angle_plus) = kinematics(&plus);
    let (com_minus, angle_minus) = kinematics(&minus);
    let dcom = com_plus
        .iter()
        .zip(&com_minus)
        .map(|(p, m)| [(p[0] - m[0]) / (2.0 * h), (p[1] - m[1]) / (2.0 * h)])
        .collect();
    let dangle = angle_plus
        .iter()
        .zip(&angle_minus)
        .map(|(p, m)| (p - m) / (2.0 * h))
        .collect();
    (dcom, dangle)
}

/// `H_kl = sum_j m_j dc_j/dq_k . dc_j/dq_l + I_c,j dtheta_j/dq_k dtheta_j/dq_l`
/// from the numerical pose Jacobians.
fn numeric_mass_matrix(chain: &ArticulatedChain2D) -> Vec<f64> {
    let n = chain.links.len();
    let derivatives: Vec<_> = (0..n).map(|k| pose_derivative(chain, k)).collect();
    let mut h = vec![0.0; n * n];
    for k in 0..n {
        for l in 0..n {
            for j in 0..n {
                let m = chain.links[j].body.mass;
                let dot = derivatives[k].0[j][0] * derivatives[l].0[j][0]
                    + derivatives[k].0[j][1] * derivatives[l].0[j][1];
                h[k * n + l] +=
                    m * dot + chain.com_inertia(j) * derivatives[k].1[j] * derivatives[l].1[j];
            }
        }
    }
    h
}

/// Gravitational potential energy `-sum m g . c` for any gravity vector.
fn potential(chain: &ArticulatedChain2D, gravity: [f64; 2]) -> f64 {
    let (coms, _) = kinematics(chain);
    chain
        .links
        .iter()
        .zip(&coms)
        .map(|(l, c)| -l.body.mass * (gravity[0] * c[0] + gravity[1] * c[1]))
        .sum()
}

/// Coriolis/centrifugal torques from the Christoffel symbols of the first
/// kind of the (finite-difference) derivative of the mass matrix.
fn numeric_bias(chain: &ArticulatedChain2D) -> Vec<f64> {
    let n = chain.links.len();
    let h = 1e-6;
    let derivative: Vec<Vec<f64>> = (0..n)
        .map(|m| {
            let mut plus = chain.clone();
            plus.links[m].joint_angle += h;
            let mut minus = chain.clone();
            minus.links[m].joint_angle -= h;
            plus.mass_matrix()
                .iter()
                .zip(minus.mass_matrix().iter())
                .map(|(p, q)| (p - q) / (2.0 * h))
                .collect()
        })
        .collect();
    let qd = chain.joint_velocities();
    (0..n)
        .map(|k| {
            let mut c = 0.0;
            for l in 0..n {
                for m in 0..n {
                    c += 0.5
                        * (derivative[m][k * n + l] + derivative[l][k * n + m]
                            - derivative[k][l * n + m])
                        * qd[l]
                        * qd[m];
                }
            }
            c
        })
        .collect()
}

fn mat_vec(h: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|k| (0..n).map(|l| h[k * n + l] * x[l]).sum())
        .collect()
}

// ── gravity ─────────────────────────────────────────────────────────────────

#[test]
fn bent_two_link_chain_needs_the_holding_torques_of_the_audit_example() {
    // q = [0, pi/2], g = (0, -9.81), no acceleration:
    // tau_0 = g (m_0 l_0 / 2 + m_1 l_0), tau_1 = 0.
    let (m0, l0, m1, l1) = (2.0, 1.5, 3.0, 0.8);
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(m0, l0, None);
    chain.add_link(m1, l1, Some(0));
    chain.links[1].joint_angle = FRAC_PI_2;
    let tau = chain.inverse_dynamics(&[0.0, 0.0], [0.0, -G]);
    assert!(
        (tau[0] - G * (m0 * l0 / 2.0 + m1 * l0)).abs() < 1e-10,
        "tau0 = {}",
        tau[0]
    );
    assert!(tau[1].abs() < 1e-10, "tau1 = {}", tau[1]);
    // The same numbers from the dedicated entry point.
    assert_close(&chain.gravity_torques([0.0, -G]), &tau, 1e-12, "G(q)");
}

#[test]
fn horizontal_chain_holding_torques_include_every_distal_lever_arm() {
    let (m, l) = ([1.0, 2.0, 3.0], [1.0, 0.8, 0.6]);
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(m[0], l[0], None);
    chain.add_link(m[1], l[1], Some(0));
    chain.add_link(m[2], l[2], Some(1));
    let tau = chain.gravity_torques([0.0, -G]);
    let expected = [
        G * (m[0] * l[0] / 2.0 + m[1] * (l[0] + l[1] / 2.0) + m[2] * (l[0] + l[1] + l[2] / 2.0)),
        G * (m[1] * l[1] / 2.0 + m[2] * (l[1] + l[2] / 2.0)),
        G * m[2] * l[2] / 2.0,
    ];
    assert_close(&tau, &expected, 1e-12, "static horizontal chain");
}

#[test]
fn single_link_gravity_torque_follows_the_cosine_of_its_angle() {
    let (m, l) = (1.5, 2.0);
    for q in [0.0, 0.5, FRAC_PI_2, 2.0, PI, -1.0, -FRAC_PI_2] {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(m, l, None);
        chain.links[0].joint_angle = q;
        let tau = chain.gravity_torques([0.0, -G]);
        let expected = m * G * l / 2.0 * q.cos();
        assert!(
            (tau[0] - expected).abs() < 1e-12,
            "q = {q}: {} vs {expected}",
            tau[0]
        );
    }
    // A horizontal link needs a positive (counter-clockwise) holding torque:
    // the old sign was inverted.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(m, l, None);
    assert!(chain.gravity_torques([0.0, -G])[0] > 0.0);
    // Gravity along +x pulls a link pointing along +y clockwise, which a
    // counter-clockwise holding torque m g l/2 sin(q) opposes.
    chain.links[0].joint_angle = FRAC_PI_2;
    let sideways = chain.gravity_torques([G, 0.0])[0];
    assert!((sideways - m * G * l / 2.0).abs() < 1e-12, "{sideways}");
}

#[test]
fn gravity_torques_are_the_gradient_of_the_potential_energy() {
    let mut rng = Lcg(0xA11CE);
    for case in 0..120 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, case % 2 == 0);
        let q = rng.vector(n, -3.0, 3.0);
        let chain = build(&specs, &q, &vec![0.0; n]);
        let gravity = [rng.range(-5.0, 5.0), rng.range(-12.0, 12.0)];
        let analytic = chain.gravity_torques(gravity);
        let h = 1e-6;
        let numeric: Vec<f64> = (0..n)
            .map(|k| {
                let mut plus = chain.clone();
                plus.links[k].joint_angle += h;
                let mut minus = chain.clone();
                minus.links[k].joint_angle -= h;
                (potential(&plus, gravity) - potential(&minus, gravity)) / (2.0 * h)
            })
            .collect();
        assert_close(&analytic, &numeric, 1e-6, &format!("case {case} dV/dq"));
    }
}

#[test]
fn gravity_torques_agree_with_the_existing_potential_energy() {
    // `potential_energy` of the chain uses |g| and the centre-of-mass heights.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(2.0, 0.7, Some(0));
    chain.add_link(1.5, 0.9, Some(0));
    chain.links[0].joint_angle = 0.4;
    chain.links[1].joint_angle = -0.8;
    chain.links[2].joint_angle = 1.1;
    let tau = chain.gravity_torques([0.0, -G]);
    let h = 1e-6;
    for (k, tau_k) in tau.iter().enumerate() {
        let mut plus = chain.clone();
        plus.links[k].joint_angle += h;
        let mut minus = chain.clone();
        minus.links[k].joint_angle -= h;
        let slope = (plus.potential_energy([0.0, 0.0], 0.0, G)
            - minus.potential_energy([0.0, 0.0], 0.0, G))
            / (2.0 * h);
        assert!(
            (tau_k - slope).abs() < 1e-6,
            "joint {k}: {tau_k} vs {slope}"
        );
    }
}

// ── mass matrix ─────────────────────────────────────────────────────────────

#[test]
fn mass_matrix_of_one_link_is_its_inertia_about_the_joint() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(3.0, 2.0, None);
    // m l^2 / 3 for the uniform rod created by add_link.
    assert_close(&chain.mass_matrix(), &[4.0], 1e-14, "H");
    // A custom inertia about the joint is used as is.
    chain.links[0].body.inertia = 5.5;
    assert_close(&chain.mass_matrix(), &[5.5], 1e-14, "H");
}

#[test]
fn mass_matrix_of_a_two_link_arm_matches_the_closed_form() {
    let mut rng = Lcg(0xBEEF);
    for _ in 0..50 {
        let (m1, l1) = (rng.range(0.3, 3.0), rng.range(0.3, 1.5));
        let (m2, l2) = (rng.range(0.3, 3.0), rng.range(0.3, 1.5));
        let q2 = rng.range(-PI, PI);
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(m1, l1, None);
        chain.add_link(m2, l2, Some(0));
        chain.links[0].joint_angle = rng.range(-PI, PI);
        chain.links[1].joint_angle = q2;
        let h00 = m1 * l1 * l1 / 3.0 + m2 * (l1 * l1 + l2 * l2 / 3.0 + l1 * l2 * q2.cos());
        let h01 = m2 * l2 * l2 / 3.0 + m2 * l1 * l2 * q2.cos() / 2.0;
        let h11 = m2 * l2 * l2 / 3.0;
        assert_close(&chain.mass_matrix(), &[h00, h01, h01, h11], 1e-12, "H");
    }
}

#[test]
fn mass_matrix_depends_on_the_configuration_and_couples_the_joints() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    let straight = chain.mass_matrix();
    chain.links[1].joint_angle = PI;
    let folded = chain.mass_matrix();
    // The coupling term m l^2/3 + m l^2 cos(q2)/2 changes with q2 ...
    assert!((straight[1] - 5.0 / 6.0).abs() < 1e-12, "{}", straight[1]);
    assert!(
        (folded[1] - (1.0 / 3.0 - 0.5)).abs() < 1e-12,
        "{}",
        folded[1]
    );
    // ... and so does the inertia seen by the first joint.
    assert!(straight[0] > folded[0] + 1.0);
    assert_eq!(straight[1], straight[2]);
    assert_eq!(folded[1], folded[2]);
}

#[test]
fn mass_matrix_matches_the_planar_chain_implementation() {
    let mut rng = Lcg(0xC0FFEE);
    for case in 0..80 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, false);
        let q = rng.vector(n, -3.0, 3.0);
        let chain = build(&specs, &q, &vec![0.0; n]);
        let oracle = joint_space_inertia_matrix_planar(
            n,
            &masses(&specs),
            &lengths(&specs),
            &inertias_com(&specs),
            &q,
        );
        assert_close(
            &chain.mass_matrix(),
            &oracle,
            1e-10,
            &format!("case {case} H"),
        );
    }
}

#[test]
fn mass_matrix_matches_the_numerical_jacobian_of_any_tree() {
    let mut rng = Lcg(0xD00D);
    for case in 0..80 {
        let n = 2 + case % 6;
        let specs = random_specs(&mut rng, n, true);
        let q = rng.vector(n, -3.0, 3.0);
        let chain = build(&specs, &q, &vec![0.0; n]);
        assert_close(
            &chain.mass_matrix(),
            &numeric_mass_matrix(&chain),
            1e-6,
            &format!("case {case} H"),
        );
    }
}

#[test]
fn mass_matrix_is_symmetric_and_positive_definite() {
    let mut rng = Lcg(0xFACE);
    for case in 0..60 {
        let n = 1 + case % 7;
        let specs = random_specs(&mut rng, n, true);
        let q = rng.vector(n, -3.0, 3.0);
        let chain = build(&specs, &q, &vec![0.0; n]);
        let h = chain.mass_matrix();
        for k in 0..n {
            for l in 0..n {
                assert!((h[k * n + l] - h[l * n + k]).abs() < 1e-12, "case {case}");
            }
        }
        assert!(
            DMatrix::from_row_slice(n, n, &h).cholesky().is_some(),
            "case {case}: H is not positive definite"
        );
    }
}

#[test]
fn sibling_branches_do_not_couple() {
    // Link 0 carries two child links; links 1 and 2 hang off its tip.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(2.0, 0.8, Some(0));
    chain.add_link(1.5, 0.6, Some(0));
    chain.links[0].joint_angle = 0.3;
    chain.links[1].joint_angle = 0.9;
    chain.links[2].joint_angle = -0.7;
    let h = chain.mass_matrix();
    // Entries (1, 2) and (2, 1) of the row-major 3x3 matrix.
    assert_eq!(h[5], 0.0);
    assert_eq!(h[7], 0.0);
    assert!(h[1] != 0.0 && h[2] != 0.0, "children couple to the parent");
    // The parent joint sees both subtrees: it is the sum of the contributions.
    let numeric = numeric_mass_matrix(&chain);
    assert_close(&h, &numeric, 1e-6, "tree H");
}

#[test]
fn separate_roots_are_independent_pendulums() {
    // Two links attached to the base: no coupling at all.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(2.0, 0.5, None);
    chain.links[0].joint_angle = 0.4;
    chain.links[1].joint_angle = 1.2;
    let h = chain.mass_matrix();
    assert_eq!(h[1], 0.0);
    assert_eq!(h[2], 0.0);
    assert!((h[0] - 1.0 / 3.0).abs() < 1e-14);
    assert!((h[3] - 2.0 * 0.25 / 3.0).abs() < 1e-14);
    let tau = chain.gravity_torques([0.0, -G]);
    assert!((tau[0] - G * 0.5 * 0.4_f64.cos()).abs() < 1e-12);
    assert!((tau[1] - 2.0 * G * 0.25 * 1.2_f64.cos()).abs() < 1e-12);
}

// ── kinetic energy and Coriolis ─────────────────────────────────────────────

#[test]
fn kinetic_energy_is_half_qd_h_qd_and_matches_the_numerical_velocities() {
    let mut rng = Lcg(0x5EED);
    for case in 0..100 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, case % 2 == 1);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let chain = build(&specs, &q, &qd);
        let energy = chain.total_kinetic_energy();
        let h = chain.mass_matrix();
        let quadratic: f64 = 0.5
            * qd.iter()
                .zip(mat_vec(&h, &qd))
                .map(|(a, b)| a * b)
                .sum::<f64>();
        assert!(
            (energy - quadratic).abs() < 1e-9 * (1.0 + energy),
            "case {case}: {energy} vs 1/2 qd.H.qd = {quadratic}"
        );
        // Independent: velocities of the centres of mass by finite differences.
        let numeric_h = numeric_mass_matrix(&chain);
        let numeric: f64 = 0.5
            * qd.iter()
                .zip(mat_vec(&numeric_h, &qd))
                .map(|(a, b)| a * b)
                .sum::<f64>();
        assert!(
            (energy - numeric).abs() < 1e-5 * (1.0 + energy),
            "case {case}: {energy} vs {numeric}"
        );
    }
}

#[test]
fn kinetic_energy_matches_the_planar_chain_implementation() {
    let mut rng = Lcg(0x7777);
    for case in 0..60 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, false);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let chain = build(&specs, &q, &qd);
        let oracle = kinetic_energy_planar(
            n,
            &masses(&specs),
            &lengths(&specs),
            &inertias_com(&specs),
            &q,
            &qd,
        );
        assert!(
            (chain.total_kinetic_energy() - oracle).abs() < 1e-9 * (1.0 + oracle),
            "case {case}"
        );
    }
}

#[test]
fn kinetic_energy_includes_the_translation_of_the_links() {
    // Joint 0 spins and joint 1 does not: the straight chain rotates rigidly
    // about the base, so the second link translates with its centre of mass at
    // 1.5 m. The old model summed rotor energies of the individual joints and
    // saw only joint 0 (0.67 J).
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    chain.links[0].joint_velocity = 2.0;
    chain.links[1].joint_velocity = 0.0;
    // Rigid rotation of the straight chain about the base: both links rotate at
    // 2 rad/s, T = 1/2 (I_0 + I_1 + m l^2) 4 with I about the joints.
    let rigid = 0.5 * 4.0 * (1.0 / 3.0 + (1.0 / 12.0 + 1.0 * 1.5 * 1.5));
    assert!((chain.total_kinetic_energy() - rigid).abs() < 1e-12);
    let old_rotor_energy = 0.5 * (1.0 / 3.0) * 4.0;
    assert!(chain.total_kinetic_energy() > 2.0 * old_rotor_energy);
}

#[test]
fn coriolis_matches_the_planar_chain_implementation() {
    let mut rng = Lcg(0x1234);
    for case in 0..80 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, false);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let chain = build(&specs, &q, &qd);
        let oracle = coriolis_bias_forces(n, &masses(&specs), &lengths(&specs), &q, &qd);
        assert_close(
            &chain.coriolis_forces(),
            &oracle,
            1e-9,
            &format!("case {case} C"),
        );
    }
}

#[test]
fn coriolis_matches_the_christoffel_symbols_of_any_tree() {
    let mut rng = Lcg(0x4321);
    for case in 0..60 {
        let n = 2 + case % 5;
        let specs = random_specs(&mut rng, n, true);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let chain = build(&specs, &q, &qd);
        assert_close(
            &chain.coriolis_forces(),
            &numeric_bias(&chain),
            1e-5,
            &format!("case {case} C"),
        );
    }
}

#[test]
fn a_single_spinning_link_and_a_rigidly_rotating_straight_chain_have_no_coriolis_force() {
    let mut single = ArticulatedChain2D::new();
    single.add_link(2.0, 1.2, None);
    single.links[0].joint_velocity = 7.0;
    assert!(single.coriolis_forces()[0].abs() < 1e-12);

    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    chain.links[0].joint_velocity = 3.0;
    assert!(chain.coriolis_forces().iter().all(|c| c.abs() < 1e-12));
    // Bending the second joint while the first rotates does produce a torque.
    chain.links[1].joint_angle = 0.8;
    assert!(chain.coriolis_forces()[0].abs() > 1e-3 || chain.coriolis_forces()[1].abs() > 1e-3);
}

// ── inverse dynamics ────────────────────────────────────────────────────────

#[test]
fn inverse_dynamics_of_a_free_single_link_is_the_joint_inertia() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    let tau = chain.inverse_dynamics(&[3.0], [0.0, 0.0]);
    assert!((tau[0] - 1.0).abs() < 1e-12, "tau = {}", tau[0]);
    // With gravity and an angle: tau = I qdd + m g l/2 cos q.
    chain.links[0].joint_angle = 0.6;
    let tau = chain.inverse_dynamics(&[3.0], [0.0, -G]);
    assert!((tau[0] - (1.0 + G * 0.5 * 0.6_f64.cos())).abs() < 1e-12);
}

#[test]
fn inverse_dynamics_matches_the_planar_chain_implementation() {
    let mut rng = Lcg(0x9999);
    for case in 0..100 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, false);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let qdd = rng.vector(n, -8.0, 8.0);
        let gravity = [rng.range(-5.0, 5.0), rng.range(-12.0, 12.0)];
        let chain = build(&specs, &q, &qd);
        let oracle = rnea_3d_revolute_z(
            n,
            &masses(&specs),
            &inertias_com(&specs),
            &lengths(&specs),
            &q,
            &qd,
            &qdd,
            [gravity[0], gravity[1], 0.0],
        );
        assert_close(
            &chain.inverse_dynamics(&qdd, gravity),
            &oracle,
            1e-9,
            &format!("case {case} tau"),
        );
        // The gravity torques equal the dedicated planar-chain function.
        let holding = gravity_compensation_torques(
            n,
            &masses(&specs),
            &lengths(&specs),
            &q,
            [gravity[0], gravity[1], 0.0],
        );
        assert_close(&chain.gravity_torques(gravity), &holding, 1e-9, "G");
    }
}

#[test]
fn inverse_dynamics_is_h_qdd_plus_c_plus_g_for_any_tree() {
    let mut rng = Lcg(0x2468);
    for case in 0..100 {
        let n = 1 + case % 7;
        let specs = random_specs(&mut rng, n, true);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let qdd = rng.vector(n, -8.0, 8.0);
        let gravity = [rng.range(-5.0, 5.0), rng.range(-12.0, 12.0)];
        let chain = build(&specs, &q, &qd);
        let tau = chain.inverse_dynamics(&qdd, gravity);
        let h_qdd = mat_vec(&chain.mass_matrix(), &qdd);
        let c = chain.coriolis_forces();
        let g = chain.gravity_torques(gravity);
        let expected: Vec<f64> = (0..n).map(|k| h_qdd[k] + c[k] + g[k]).collect();
        assert_close(&tau, &expected, 1e-9, &format!("case {case} tau"));
        // And with the numerical mass matrix, Christoffel bias and potential.
        let numeric: Vec<f64> = mat_vec(&numeric_mass_matrix(&chain), &qdd)
            .iter()
            .zip(numeric_bias(&chain))
            .zip(&g)
            .map(|((a, b), c)| a + b + c)
            .collect();
        assert_close(&tau, &numeric, 1e-4, &format!("case {case} numeric tau"));
    }
}

#[test]
fn short_and_long_acceleration_slices_are_padded_or_truncated() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(2.0, 0.5, Some(0));
    chain.links[0].joint_angle = 0.3;
    chain.links[1].joint_angle = -0.4;
    chain.links[0].joint_velocity = 1.0;
    let full = chain.inverse_dynamics(&[2.0, 0.0], [0.0, -G]);
    assert_close(
        &chain.inverse_dynamics(&[2.0], [0.0, -G]),
        &full,
        1e-14,
        "short",
    );
    assert_close(
        &chain.inverse_dynamics(&[2.0, 0.0, 9.0, 9.0], [0.0, -G]),
        &full,
        1e-14,
        "long",
    );
    let none = chain.inverse_dynamics(&[], [0.0, -G]);
    assert_close(
        &none,
        &chain.inverse_dynamics(&[0.0, 0.0], [0.0, -G]),
        1e-14,
        "empty",
    );
}

// ── forward dynamics ────────────────────────────────────────────────────────

#[test]
fn torque_on_the_distal_joint_accelerates_the_proximal_joint_too() {
    // Two uniform unit rods, straight, no gravity: H = [[8/3, 5/6], [5/6, 1/3]]
    // and q'' = H^-1 [0, 1] = [-30/7, 96/7]. The old model applied tau_i / I_i to
    // every joint independently, so the first joint never reacted.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    let qdd = chain
        .forward_dynamics(&[0.0, 1.0], [0.0, 0.0])
        .expect("positive definite inertia matrix");
    assert_close(&qdd, &[-30.0 / 7.0, 96.0 / 7.0], 1e-12, "q''");

    let dt = 1e-3;
    chain.apply_torques(&[0.0, 1.0], dt);
    assert!((chain.links[0].joint_velocity + 30.0 / 7.0 * dt).abs() < 1e-12);
    assert!((chain.links[1].joint_velocity - 96.0 / 7.0 * dt).abs() < 1e-12);
    assert!((chain.links[0].joint_angle - chain.links[0].joint_velocity * dt).abs() < 1e-15);
}

#[test]
fn forward_dynamics_inverts_the_inverse_dynamics_for_any_tree() {
    let mut rng = Lcg(0x1357);
    for case in 0..100 {
        let n = 1 + case % 7;
        let specs = random_specs(&mut rng, n, true);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let qdd = rng.vector(n, -8.0, 8.0);
        let gravity = [rng.range(-5.0, 5.0), rng.range(-12.0, 12.0)];
        let chain = build(&specs, &q, &qd);
        let tau = chain.inverse_dynamics(&qdd, gravity);
        let recovered = chain
            .forward_dynamics(&tau, gravity)
            .expect("H is positive definite");
        assert_close(&recovered, &qdd, 1e-7, &format!("case {case} q''"));
    }
}

#[test]
fn forward_dynamics_matches_the_articulated_body_algorithm() {
    let mut rng = Lcg(0x8642);
    for case in 0..80 {
        let n = 1 + case % 6;
        let specs = random_specs(&mut rng, n, false);
        let q = rng.vector(n, -3.0, 3.0);
        let qd = rng.vector(n, -4.0, 4.0);
        let tau = rng.vector(n, -10.0, 10.0);
        let gravity = [rng.range(-5.0, 5.0), rng.range(-12.0, 12.0)];
        let chain = build(&specs, &q, &qd);
        let oracle = aba_revolute_chain(
            n,
            &inertias_com(&specs),
            &masses(&specs),
            &lengths(&specs),
            &q,
            &qd,
            &tau,
            [gravity[0], gravity[1], 0.0],
        );
        let qdd = chain
            .forward_dynamics(&tau, gravity)
            .expect("H is positive definite");
        assert_close(&qdd, &oracle, 1e-7, &format!("case {case} q''"));
    }
}

#[test]
fn a_released_horizontal_link_falls_with_the_pendulum_acceleration() {
    // Uniform rod pivoting about its end: q'' = -(3 g / 2 l) cos q.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 2.0, None);
    let qdd = chain
        .forward_dynamics(&[0.0], [0.0, -G])
        .expect("single link");
    assert!((qdd[0] + 1.5 * G / 2.0).abs() < 1e-12, "q'' = {}", qdd[0]);
    // And the integration step follows it (semi-implicit Euler).
    let dt = 1e-3;
    chain.apply_torques_with_gravity(&[0.0], [0.0, -G], dt);
    assert!((chain.links[0].joint_velocity + 1.5 * G / 2.0 * dt).abs() < 1e-12);
    assert!(chain.links[0].joint_angle < 0.0);
}

#[test]
fn apply_torques_ignores_gravity() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.apply_torques(&[0.0], 0.1);
    assert_eq!(chain.links[0].joint_velocity, 0.0);
    assert_eq!(chain.links[0].joint_angle, 0.0);
}

#[test]
fn unforced_double_pendulum_conserves_energy_to_the_integration_error() {
    let make = || {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.2, 0.9, None);
        chain.add_link(0.8, 0.7, Some(0));
        chain.links[0].joint_angle = 0.5;
        chain.links[1].joint_angle = -0.3;
        chain.links[0].joint_velocity = 0.1;
        chain.links[1].joint_velocity = 0.05;
        chain
    };
    let energy = |chain: &ArticulatedChain2D| {
        chain.total_kinetic_energy() + chain.potential_energy([0.0, 0.0], 0.0, G)
    };
    let drift = |dt: f64, steps: usize| {
        let mut chain = make();
        let e0 = energy(&chain);
        let mut worst: f64 = 0.0;
        for _ in 0..steps {
            chain.apply_torques_with_gravity(&[0.0, 0.0], [0.0, -G], dt);
            worst = worst.max((energy(&chain) - e0).abs());
        }
        worst
    };
    let coarse = drift(2e-4, 5000);
    let fine = drift(1e-4, 10000);
    assert!(fine < 0.05, "energy drift {fine} J");
    // First-order scheme: halving the step halves the (bounded) energy error.
    assert!(fine < 0.75 * coarse, "coarse {coarse} J vs fine {fine} J");
}

#[test]
fn pendulum_swing_agrees_with_the_closed_form_small_angle_period() {
    // Uniform rod pivoting about its end: omega^2 = 3 g / (2 l) near q = -pi/2.
    let l = 1.0;
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, l, None);
    chain.links[0].joint_angle = -FRAC_PI_2 + 0.01;
    let dt = 1e-4;
    let mut upward_crossings = Vec::new();
    let mut previous = chain.links[0].joint_angle + FRAC_PI_2;
    for step in 0..30_000 {
        chain.apply_torques_with_gravity(&[0.0], [0.0, -G], dt);
        let offset = chain.links[0].joint_angle + FRAC_PI_2;
        if previous < 0.0 && offset >= 0.0 {
            let fraction = -previous / (offset - previous);
            upward_crossings.push((step as f64 + fraction) * dt);
        }
        previous = offset;
    }
    assert!(upward_crossings.len() >= 2, "{upward_crossings:?}");
    let period = upward_crossings[1] - upward_crossings[0];
    let expected = 2.0 * PI / (1.5 * G / l).sqrt();
    assert!(
        (period - expected).abs() < 5e-3 * expected,
        "period {period} vs {expected}"
    );
}

// ── degenerate input ────────────────────────────────────────────────────────

#[test]
fn links_without_inertia_are_locked_by_the_integration_and_an_error_for_the_strict_solver() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(0.0, 1.0, Some(0)); // massless, inertia-less pointer
    let err = chain
        .forward_dynamics(&[1.0, 5.0], [0.0, 0.0])
        .expect_err("singular inertia matrix");
    assert!(
        matches!(err, Error::ForwardDynamicsFailure { .. }),
        "{err:?}"
    );

    chain.apply_torques(&[1.0, 5.0], 0.1);
    // Joint 0 carries I = 1/3: q'' = 3; the massless joint cannot accelerate.
    assert!((chain.links[0].joint_velocity - 0.3).abs() < 1e-9);
    assert!(chain.links[1].joint_velocity.abs() < 1e-12);
    assert!(chain.links[0].joint_angle.is_finite());
}

#[test]
fn inertia_below_the_parallel_axis_term_is_clamped_at_zero() {
    // A rod of mass 1 and length 1 cannot have less than m (l/2)^2 = 0.25 about
    // its end: the inertia about the centre of mass is clamped at zero.
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.links[0].body.inertia = 0.1;
    assert_close(&chain.mass_matrix(), &[0.25], 1e-14, "H");
    chain.links[0].joint_velocity = 2.0;
    assert!((chain.total_kinetic_energy() - 0.5).abs() < 1e-14);
}

#[test]
fn parent_indices_that_do_not_refer_to_an_earlier_link_attach_to_the_base() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, Some(0)); // itself
    chain.add_link(2.0, 0.5, Some(7)); // outside the chain
    chain.add_link(1.5, 0.8, Some(1));
    chain.links[0].joint_angle = 0.4;
    chain.links[1].joint_angle = -0.3;
    chain.links[2].joint_angle = 0.9;
    chain.links[1].joint_velocity = 1.5;

    let mut reference = ArticulatedChain2D::new();
    reference.add_link(1.0, 1.0, None);
    reference.add_link(2.0, 0.5, None);
    reference.add_link(1.5, 0.8, Some(1));
    for i in 0..3 {
        reference.links[i].joint_angle = chain.links[i].joint_angle;
        reference.links[i].joint_velocity = chain.links[i].joint_velocity;
    }

    let gravity = [0.0, -G];
    let qdd = [1.0, -2.0, 0.5];
    assert_close(&chain.mass_matrix(), &reference.mass_matrix(), 1e-14, "H");
    assert_close(
        &chain.inverse_dynamics(&qdd, gravity),
        &reference.inverse_dynamics(&qdd, gravity),
        1e-14,
        "tau",
    );
    assert_close(
        &chain.coriolis_forces(),
        &reference.coriolis_forces(),
        1e-14,
        "C",
    );
    assert!((chain.total_kinetic_energy() - reference.total_kinetic_energy()).abs() < 1e-14);
    assert!(chain.forward_dynamics(&[0.0; 3], gravity).is_ok());
}

#[test]
fn empty_chain_has_empty_dynamics() {
    let mut chain = ArticulatedChain2D::new();
    assert!(chain.mass_matrix().is_empty());
    assert!(chain.gravity_torques([0.0, -G]).is_empty());
    assert!(chain.coriolis_forces().is_empty());
    assert!(chain.inverse_dynamics(&[1.0], [0.0, -G]).is_empty());
    assert_eq!(chain.forward_dynamics(&[1.0], [0.0, -G]), Ok(Vec::new()));
    assert_eq!(chain.total_kinetic_energy(), 0.0);
    chain.apply_torques(&[1.0], 0.1);
    assert_eq!(chain.link_count(), 0);
}

#[test]
fn non_finite_state_gives_nan_accelerations_without_hanging() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    chain.links[1].joint_angle = f64::NAN;
    assert!(chain.forward_dynamics(&[0.0, 0.0], [0.0, -G]).is_err());
    chain.apply_torques_with_gravity(&[0.0, 0.0], [0.0, -G], 0.01);
    assert!(chain.links[0].joint_velocity.is_nan());
}

#[test]
fn missing_torques_default_to_zero() {
    let mut chain = ArticulatedChain2D::new();
    chain.add_link(1.0, 1.0, None);
    chain.add_link(1.0, 1.0, Some(0));
    let short = chain
        .forward_dynamics(&[1.0], [0.0, 0.0])
        .expect("positive definite");
    let full = chain
        .forward_dynamics(&[1.0, 0.0], [0.0, 0.0])
        .expect("positive definite");
    assert_close(&short, &full, 1e-14, "q''");
}
