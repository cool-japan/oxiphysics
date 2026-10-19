// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the multibody dynamics: the mass matrix (CRBA), the
//! gravity torques and the Coriolis vector are checked against closed-form
//! two-link results and against independent finite-difference oracles on a
//! branched, non-planar chain with full inertia tensors.

use super::*;

const GRAVITY: [f64; 3] = [0.0, -9.81, 0.0];

// ── fixtures ────────────────────────────────────────────────────────────

/// Four links, branching (`0 -> 1 -> 3` and `0 -> 2`), non-parallel joint axes,
/// non-zero joint offsets, off-axis centres of mass and inertia tensors with
/// products of inertia.
fn spatial_chain() -> MultibodyChain {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let links = vec![
        Link::with_inertia(
            2.0,
            [[0.30, 0.02, -0.01], [0.02, 0.25, 0.03], [-0.01, 0.03, 0.20]],
            [0.1, 0.25, -0.05],
            [0.0, 0.0, 1.0],
            [0.2, -0.1, 0.4],
        ),
        Link::with_inertia(
            1.5,
            [[0.10, 0.01, 0.0], [0.01, 0.12, -0.02], [0.0, -0.02, 0.08]],
            [0.3, 0.0, 0.1],
            [0.0, 1.0, 0.0],
            [0.0, 0.5, 0.1],
        ),
        Link::with_inertia(
            0.8,
            [[0.05, 0.0, 0.01], [0.0, 0.04, 0.0], [0.01, 0.0, 0.03]],
            [0.0, 0.2, 0.15],
            [1.0, 0.0, 0.0],
            [0.3, 0.2, -0.2],
        ),
        Link::with_inertia(
            1.1,
            [[0.06, -0.01, 0.0], [-0.01, 0.07, 0.01], [0.0, 0.01, 0.05]],
            [0.15, 0.1, -0.05],
            [s, s, 0.0],
            [0.4, 0.0, 0.2],
        ),
    ];
    let mut chain = MultibodyChain::new(links, vec![None, Some(0), Some(0), Some(1)]);
    chain.joint_angles = vec![0.4, -0.7, 1.1, 0.5];
    chain.joint_velocities = vec![0.8, -0.6, 0.5, 1.2];
    chain
}

fn planar(n: usize, angles: &[f64], velocities: &[f64]) -> MultibodyChain {
    let mut chain = build_planar_arm(n, 1.0, 1.0);
    chain.joint_angles = angles.to_vec();
    chain.joint_velocities = velocities.to_vec();
    chain
}

fn com_positions(chain: &MultibodyChain) -> Vec<[f64; 3]> {
    chain
        .forward_kinematics()
        .iter()
        .zip(chain.links.iter())
        .map(|((pos, rot), link)| vec3_add(*pos, mat3_vec_mul(*rot, link.com)))
        .collect()
}

fn potential_energy(chain: &MultibodyChain, gravity: [f64; 3]) -> f64 {
    com_positions(chain)
        .iter()
        .zip(chain.links.iter())
        .map(|(c, link)| -link.mass * vec3_dot(gravity, *c))
        .sum()
}

fn kinetic_energy(chain: &MultibodyChain) -> f64 {
    let m = chain.mass_matrix();
    let qd = &chain.joint_velocities;
    let mut t = 0.0;
    for (row, qi) in m.iter().zip(qd.iter()) {
        for (mij, qj) in row.iter().zip(qd.iter()) {
            t += 0.5 * qi * mij * qj;
        }
    }
    t
}

fn perturbed(chain: &MultibodyChain, joint: usize, delta: f64) -> MultibodyChain {
    let mut copy = chain.clone();
    copy.joint_angles[joint] += delta;
    copy
}

fn max_abs_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Mass matrix from its definition `M = Σ m Jvᵀ Jv + Jωᵀ I_world Jω`, with the
/// centre-of-mass velocity Jacobian obtained by differencing the forward
/// kinematics (independent of the CRBA code under test).
fn mass_matrix_by_definition(chain: &MultibodyChain) -> Vec<Vec<f64>> {
    let n = chain.n_dof();
    let h = 1e-6;
    let jv: Vec<Vec<[f64; 3]>> = (0..n)
        .map(|joint| {
            let plus = com_positions(&perturbed(chain, joint, h));
            let minus = com_positions(&perturbed(chain, joint, -h));
            plus.iter()
                .zip(minus.iter())
                .map(|(p, m)| vec3_scale(vec3_sub(*p, *m), 0.5 / h))
                .collect()
        })
        .collect();
    let poses = chain.forward_kinematics();
    let mut m = vec![vec![0.0; n]; n];
    for (k, link) in chain.links.iter().enumerate() {
        let rot = poses[k].1;
        let inertia_world = mat3_mul(rot, mat3_mul(link.inertia, mat3_transpose(rot)));
        let angular: Vec<[f64; 3]> = chain
            .jacobian(k)
            .iter()
            .map(|col| [col[0], col[1], col[2]])
            .collect();
        for a in 0..n {
            for b in 0..n {
                m[a][b] += link.mass * vec3_dot(jv[a][k], jv[b][k])
                    + vec3_dot(angular[a], mat3_vec_mul(inertia_world, angular[b]));
            }
        }
    }
    m
}

/// Cholesky factorisation succeeds iff the symmetric matrix is positive definite.
fn is_positive_definite(m: &[Vec<f64>]) -> bool {
    let n = m.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let partial: f64 = l[i][..j]
                .iter()
                .zip(l[j][..j].iter())
                .map(|(a, b)| a * b)
                .sum();
            if i == j {
                let d = m[i][i] - partial;
                if d <= 1e-12 {
                    return false;
                }
                l[i][i] = d.sqrt();
            } else {
                l[i][j] = (m[i][j] - partial) / l[j][j];
            }
        }
    }
    true
}

/// Deterministic pseudo-random angles in `[-pi, pi]`.
fn pseudo_random(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let unit = ((*seed >> 11) as f64) / ((1u64 << 53) as f64);
    (unit * 2.0 - 1.0) * PI
}

// ── articulated body inertia ────────────────────────────────────────────

#[test]
fn spatial_inertia_of_a_point_mass_uses_the_parallel_axis_theorem() {
    // Point mass 2 kg at (1, 2, 0) from the reference point.
    let link = Link::new(2.0, [1.0, 2.0, 0.0]);
    let abi = ArticulatedBodyInertia::from_link(&link, mat3_identity(), [0.0; 3]);
    let expected = [[8.0, -4.0, 0.0], [-4.0, 2.0, 0.0], [0.0, 0.0, 10.0]];
    for (row, expected_row) in abi.h.iter().zip(expected.iter()) {
        for (value, want) in row.iter().zip(expected_row.iter()) {
            assert!((value - want).abs() < 1e-12, "{:?}", abi.h);
        }
    }
    assert_eq!(abi.c, [2.0, 4.0, 0.0]);
    // The same mass described through a joint located at (1, 0, 0) relative to
    // the reference point and a centre of mass at (0, 2, 0) from that joint.
    let shifted = ArticulatedBodyInertia::from_link(
        &Link::new(2.0, [0.0, 2.0, 0.0]),
        mat3_identity(),
        [1.0, 0.0, 0.0],
    );
    for (row, expected_row) in shifted.h.iter().zip(abi.h.iter()) {
        for (value, want) in row.iter().zip(expected_row.iter()) {
            assert!((value - want).abs() < 1e-12);
        }
    }
}

#[test]
fn composite_of_two_point_masses_has_no_first_moment() {
    let left = Link::new(1.0, [-1.0, 0.0, 0.0]);
    let right = Link::new(1.0, [1.0, 0.0, 0.0]);
    let total = ArticulatedBodyInertia::from_link(&left, mat3_identity(), [0.0; 3]).add(
        &ArticulatedBodyInertia::from_link(&right, mat3_identity(), [0.0; 3]),
    );
    assert!((total.m - 2.0).abs() < 1e-15);
    assert!(vec3_norm(total.c) < 1e-15);
    assert!((total.h[0][0]).abs() < 1e-15);
    assert!((total.h[1][1] - 2.0).abs() < 1e-15);
    assert!((total.h[2][2] - 2.0).abs() < 1e-15);
}

#[test]
fn spatial_inertia_momentum_matches_a_rigid_rotation_about_the_reference() {
    // 3 kg point mass at (0, 1, 0) rotating about z at 2 rad/s with the
    // reference point at rest: v = w x r = (-2, 0, 0).
    let abi = ArticulatedBodyInertia::from_link(
        &Link::new(3.0, [0.0, 1.0, 0.0]),
        mat3_identity(),
        [0.0; 3],
    );
    let (angular, linear) = abi.momentum([0.0, 0.0, 2.0], [0.0; 3]);
    assert!((angular[2] - 6.0).abs() < 1e-12, "{angular:?}");
    assert!(
        (linear[0] + 6.0).abs() < 1e-12 && linear[1].abs() < 1e-12,
        "{linear:?}"
    );
}

// ── mass matrix ─────────────────────────────────────────────────────────

#[test]
fn single_link_mass_is_inertia_plus_parallel_axis_term() {
    // Uniform rod of length 1 and mass 1 about its end: 1/12 + 1/4 = 1/3.
    let chain = build_planar_arm(1, 1.0, 1.0);
    let m = chain.mass_matrix();
    assert!((m[0][0] - 1.0 / 3.0).abs() < 1e-12, "M = {}", m[0][0]);
}

#[test]
fn two_link_mass_matrix_matches_the_closed_form() {
    // m1 = m2 = 1, l1 = l2 = 1, lc = 0.5, I = 1/12 about the centre of mass.
    for &(q1, q2) in &[
        (0.0, 0.0),
        (0.3, 1.2),
        (-1.0, 2.5),
        (2.0, -0.7),
        (PI, PI / 2.0),
    ] {
        let chain = planar(2, &[q1, q2], &[0.0, 0.0]);
        let m = chain.mass_matrix();
        let c2 = f64::cos(q2);
        let i_c = 1.0 / 12.0;
        let m11 = 0.25 + (1.0 + 0.25 + c2) + 2.0 * i_c;
        let m12 = 0.25 + 0.5 * c2 + i_c;
        let m22 = 0.25 + i_c;
        assert!((m[0][0] - m11).abs() < 1e-12, "M11 {} vs {m11}", m[0][0]);
        assert!((m[0][1] - m12).abs() < 1e-12, "M12 {} vs {m12}", m[0][1]);
        assert!((m[1][0] - m12).abs() < 1e-12, "M21 {} vs {m12}", m[1][0]);
        assert!((m[1][1] - m22).abs() < 1e-12, "M22 {} vs {m22}", m[1][1]);
    }
}

#[test]
fn spatial_mass_matrix_matches_the_definition_from_jacobians() {
    let mut seed = 7;
    for trial in 0..6 {
        let mut chain = spatial_chain();
        if trial > 0 {
            for angle in chain.joint_angles.iter_mut() {
                *angle = pseudo_random(&mut seed);
            }
        }
        let m = chain.mass_matrix();
        let oracle = mass_matrix_by_definition(&chain);
        for (row, oracle_row) in m.iter().zip(oracle.iter()) {
            assert!(
                max_abs_difference(row, oracle_row) < 1e-6,
                "trial {trial}: M = {m:?}\noracle = {oracle:?}"
            );
        }
    }
}

#[test]
fn mass_matrix_is_symmetric_positive_definite_and_respects_the_tree() {
    let mut seed = 99;
    for _ in 0..20 {
        let mut chain = spatial_chain();
        for angle in chain.joint_angles.iter_mut() {
            *angle = pseudo_random(&mut seed);
        }
        let m = chain.mass_matrix();
        for (i, row) in m.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert!((value - m[j][i]).abs() < 1e-14);
            }
        }
        assert!(is_positive_definite(&m), "M = {m:?}");
        // Siblings (1, 2) and the cross-branch pair (2, 3) do not interact.
        assert_eq!(m[1][2], 0.0);
        assert_eq!(m[2][3], 0.0);
        // Ancestor relations do interact.
        assert!(m[0][1].abs() > 1e-6 || m[0][3].abs() > 1e-6);
        assert!(m[1][3].abs() > 1e-9);
    }
    for &(q1, q2, q3) in &[(0.0, 0.0, 0.0), (1.0, -2.0, 0.5), (3.0, 3.0, -3.0)] {
        let chain = planar(3, &[q1, q2, q3], &[0.0; 3]);
        assert!(is_positive_definite(&chain.mass_matrix()));
    }
}

#[test]
fn mass_matrix_does_not_depend_on_where_the_chain_is_in_the_world() {
    let reference = spatial_chain();
    let mut moved = spatial_chain();
    moved.links[0].joint_offset = vec3_add(moved.links[0].joint_offset, [1000.0, -500.0, 250.0]);
    let a = reference.mass_matrix();
    let b = moved.mass_matrix();
    for (ra, rb) in a.iter().zip(b.iter()) {
        assert!(max_abs_difference(ra, rb) < 1e-9, "{a:?} vs {b:?}");
    }
    assert!(max_abs_difference(&reference.coriolis_vector(), &moved.coriolis_vector()) < 1e-8);
    assert!(
        max_abs_difference(
            &reference.gravity_torques(GRAVITY),
            &moved.gravity_torques(GRAVITY)
        ) < 1e-8
    );
}

// ── gravity torques ─────────────────────────────────────────────────────

#[test]
fn two_link_gravity_torques_match_the_holding_torques() {
    // Horizontal arm: q = (pi/2, 0), links along -x.
    let chain = planar(2, &[PI / 2.0, 0.0], &[0.0, 0.0]);
    let g = chain.gravity_torques(GRAVITY);
    assert!((g[0] + 19.62).abs() < 1e-9, "g0 = {}", g[0]);
    assert!((g[1] + 4.905).abs() < 1e-9, "g1 = {}", g[1]);
}

#[test]
fn two_link_gravity_torques_match_the_potential_gradient_for_any_pose() {
    for &(q1, q2) in &[(0.0, 0.0), (0.4, 1.0), (-1.3, 0.8), (2.5, -2.0), (PI, 0.3)] {
        let chain = planar(2, &[q1, q2], &[0.0; 2]);
        let g = chain.gravity_torques(GRAVITY);
        let g1 = -9.81 * ((0.5 + 1.0) * q1.sin() + 0.5 * (q1 + q2).sin());
        let g2 = -9.81 * 0.5 * (q1 + q2).sin();
        assert!((g[0] - g1).abs() < 1e-9, "g0 {} vs {g1}", g[0]);
        assert!((g[1] - g2).abs() < 1e-9, "g1 {} vs {g2}", g[1]);
    }
}

#[test]
fn spatial_gravity_torques_are_the_gradient_of_the_potential_energy() {
    let gravity = [0.3, -9.81, 0.5];
    let mut seed = 11;
    for _ in 0..6 {
        let mut chain = spatial_chain();
        for angle in chain.joint_angles.iter_mut() {
            *angle = pseudo_random(&mut seed);
        }
        let h = 1e-6;
        let gradient: Vec<f64> = (0..chain.n_dof())
            .map(|j| {
                let plus = potential_energy(&perturbed(&chain, j, h), gravity);
                let minus = potential_energy(&perturbed(&chain, j, -h), gravity);
                (plus - minus) / (2.0 * h)
            })
            .collect();
        let g = chain.gravity_torques(gravity);
        assert!(
            max_abs_difference(&g, &gradient) < 1e-6,
            "g = {g:?} vs dV/dq = {gradient:?}"
        );
    }
}

#[test]
fn a_horizontal_link_falls_towards_the_hanging_position() {
    let mut chain = build_planar_arm(1, 1.0, 1.0);
    chain.joint_angles[0] = PI / 2.0;
    // q = pi/2 points along -x; gravity pulls the tip down, i.e. q increases.
    let qdd = chain.forward_dynamics(&[0.0], GRAVITY);
    assert!(
        (qdd[0] - 4.905 * 3.0).abs() < 1e-9,
        "q'' = {} expected 14.715",
        qdd[0]
    );
    chain.step_rk4(&[0.0], GRAVITY, 1e-3);
    assert!(
        chain.joint_angles[0] > PI / 2.0,
        "q = {}",
        chain.joint_angles[0]
    );
    assert!(chain.joint_velocities[0] > 0.0);
}

// ── Coriolis vector ─────────────────────────────────────────────────────

#[test]
fn two_link_coriolis_vector_matches_the_closed_form() {
    for &(q1, q2, v1, v2) in &[
        (0.3, 1.1, 0.7, -0.4),
        (-1.0, 0.6, 2.0, 1.5),
        (2.0, -2.0, -1.2, 0.9),
        (0.0, 0.5, 0.0, 3.0),
    ] {
        let chain = planar(2, &[q1, q2], &[v1, v2]);
        let c = chain.coriolis_vector();
        let h = -0.5 * q2.sin();
        let c1 = h * (2.0 * v1 * v2 + v2 * v2);
        let c2 = -h * v1 * v1;
        assert!((c[0] - c1).abs() < 1e-9, "c1 {} vs {c1}", c[0]);
        assert!((c[1] - c2).abs() < 1e-9, "c2 {} vs {c2}", c[1]);
    }
}

#[test]
fn coriolis_vector_vanishes_at_rest_and_is_quadratic_in_the_velocities() {
    let mut chain = spatial_chain();
    chain.joint_velocities = vec![0.0; 4];
    assert!(chain.coriolis_vector().iter().all(|c| *c == 0.0));

    let mut chain = spatial_chain();
    let c1 = chain.coriolis_vector();
    for v in chain.joint_velocities.iter_mut() {
        *v *= 3.0;
    }
    let c3 = chain.coriolis_vector();
    for (a, b) in c1.iter().zip(c3.iter()) {
        assert!((9.0 * a - b).abs() < 1e-10, "{a} {b}");
    }
}

#[test]
fn spatial_coriolis_vector_matches_the_christoffel_terms_of_the_mass_matrix() {
    let mut seed = 5;
    for _ in 0..4 {
        let mut chain = spatial_chain();
        for angle in chain.joint_angles.iter_mut() {
            *angle = pseudo_random(&mut seed);
        }
        for v in chain.joint_velocities.iter_mut() {
            *v = pseudo_random(&mut seed) * 0.5;
        }
        let n = chain.n_dof();
        let h = 1e-5;
        // dm[k][i][j] = dM_ij / dq_k
        let dm: Vec<Vec<Vec<f64>>> = (0..n)
            .map(|k| {
                let plus = perturbed(&chain, k, h).mass_matrix();
                let minus = perturbed(&chain, k, -h).mass_matrix();
                (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| (plus[i][j] - minus[i][j]) / (2.0 * h))
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let qd = &chain.joint_velocities;
        let mut expected = vec![0.0; n];
        for (i, e) in expected.iter_mut().enumerate() {
            for j in 0..n {
                for k in 0..n {
                    *e += (dm[k][i][j] - 0.5 * dm[i][j][k]) * qd[j] * qd[k];
                }
            }
        }
        let c = chain.coriolis_vector();
        assert!(
            max_abs_difference(&c, &expected) < 1e-5,
            "C = {c:?} vs Christoffel {expected:?}"
        );
    }
}

#[test]
fn coriolis_vector_satisfies_the_power_balance_identity() {
    // qd . C(q, qd) qd = 1/2 qd^T (dM/dt) qd, with dM/dt differenced along qd.
    let chain = spatial_chain();
    let qd = chain.joint_velocities.clone();
    let h = 1e-6;
    let shifted = |sign: f64| {
        let mut copy = chain.clone();
        for (angle, v) in copy.joint_angles.iter_mut().zip(qd.iter()) {
            *angle += sign * h * v;
        }
        copy.mass_matrix()
    };
    let plus = shifted(1.0);
    let minus = shifted(-1.0);
    let mut half_qd_mdot_qd = 0.0;
    for (i, qi) in qd.iter().enumerate() {
        for (j, qj) in qd.iter().enumerate() {
            half_qd_mdot_qd += 0.5 * qi * (plus[i][j] - minus[i][j]) / (2.0 * h) * qj;
        }
    }
    let c = chain.coriolis_vector();
    let power: f64 = qd.iter().zip(c.iter()).map(|(v, ci)| v * ci).sum();
    assert!(
        (power - half_qd_mdot_qd).abs() < 1e-6,
        "qd.C.qd = {power} vs 1/2 qd.Mdot.qd = {half_qd_mdot_qd}"
    );
}

// ── RNEA consistency and time integration ───────────────────────────────

#[test]
fn rnea_agrees_with_mass_matrix_coriolis_and_gravity() {
    let mut seed = 21;
    let negative_gravity = [-0.3, 9.81, -0.5];
    for _ in 0..5 {
        let mut chain = spatial_chain();
        for angle in chain.joint_angles.iter_mut() {
            *angle = pseudo_random(&mut seed);
        }
        for v in chain.joint_velocities.iter_mut() {
            *v = pseudo_random(&mut seed);
        }
        let qdd: Vec<f64> = (0..4).map(|_| pseudo_random(&mut seed)).collect();
        let gravity = [0.3, -9.81, 0.5];
        let from_rnea = chain.rnea(&chain.joint_velocities, &qdd, negative_gravity);
        let from_matrices = chain.inverse_dynamics(&qdd, gravity);
        assert!(
            max_abs_difference(&from_rnea, &from_matrices) < 1e-9,
            "RNEA {from_rnea:?} vs M q'' + C + g {from_matrices:?}"
        );
        let static_rnea = chain.rnea(&[], &[], negative_gravity);
        assert!(max_abs_difference(&static_rnea, &chain.gravity_torques(gravity)) < 1e-10);
    }
}

#[test]
fn forward_and_inverse_dynamics_are_inverse_to_each_other() {
    let chain = spatial_chain();
    let qdd = [0.5, -1.0, 2.0, 0.25];
    let tau = chain.inverse_dynamics(&qdd, GRAVITY);
    let recovered = chain.forward_dynamics(&tau, GRAVITY);
    assert!(max_abs_difference(&recovered, &qdd) < 1e-9, "{recovered:?}");
}

#[test]
fn rk4_conserves_energy_without_torques() {
    for gravity in [GRAVITY, [0.3, -9.81, 0.5]] {
        let mut chain = spatial_chain();
        let energy_0 = kinetic_energy(&chain) + potential_energy(&chain, gravity);
        for _ in 0..1000 {
            chain.step_rk4(&[0.0; 4], gravity, 1e-3);
        }
        let energy = kinetic_energy(&chain) + potential_energy(&chain, gravity);
        assert!(
            (energy - energy_0).abs() < 1e-7,
            "energy drifted from {energy_0} to {energy}"
        );
        assert!(chain.joint_angles.iter().all(|q| q.is_finite()));
        // The motion is not trivial.
        assert!((chain.joint_angles[0] - 0.4).abs() > 0.1);
    }
}

#[test]
fn small_oscillations_of_a_hanging_rod_have_the_physical_period() {
    // Physical pendulum: w^2 = m g lc / (I_c + m lc^2) = 9.81 * 0.5 / (1/3).
    let omega = (9.81_f64 * 0.5 / (1.0 / 3.0)).sqrt();
    let period = 2.0 * PI / omega;
    let amplitude = 0.02;
    let dt = 1e-3;
    let mut chain = build_planar_arm(1, 1.0, 1.0);
    chain.joint_angles[0] = PI + amplitude;

    let half_steps = (0.5 * period / dt).round() as usize;
    for _ in 0..half_steps {
        chain.step_rk4(&[0.0], GRAVITY, dt);
    }
    assert!(
        (chain.joint_angles[0] - (PI - amplitude)).abs() < 2e-4 * amplitude + 1e-5,
        "half period: q - pi = {}",
        chain.joint_angles[0] - PI
    );
    for _ in 0..half_steps {
        chain.step_rk4(&[0.0], GRAVITY, dt);
    }
    assert!(
        (chain.joint_angles[0] - (PI + amplitude)).abs() < 2e-4 * amplitude + 1e-5,
        "full period: q - pi = {}",
        chain.joint_angles[0] - PI
    );
}
