// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the stretch / shear coupling of the [`DiscreteCosseratRod`].
//!
//! Before the coupling the nodes only felt a stretch force along the chord; the
//! bending stiffness acted on the frames alone and never reached the positions
//! (a cantilever did not sag, an intrinsically curved rod stayed straight), and the
//! joint moments were not paired.  The tests check the invariants of the coupled
//! model (momentum balance, exact consistency with the stored energy) and the
//! engineering results it must reproduce (Timoshenko cantilever, constant-curvature
//! arc).

use super::*;

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];
const Z: [f64; 3] = [0.0, 0.0, 1.0];

/// Soft, light rod material: explicit stepping at `dt = 1e-4` is far inside the
/// stability limit (stiffest mode ~ 2e3 rad/s).
fn soft_props() -> CosseratRodProperties {
    CosseratRodProperties::circular(0.01, 1.0e6, 4.0e5, 1000.0)
}

fn rotate_frame(rod: &mut DiscreteCosseratRod, node: usize, axis: [f64; 3], angle: f64) {
    rod.nodes[node].frame = rod.nodes[node].frame.rotate(rodrigues(axis, angle));
}

/// Deterministic pseudo-random numbers in `[-1, 1]` (xorshift).
fn noise(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 52) as f64 * 2.0 - 1.0
}

/// A free (nothing pinned) five-node rod with displaced nodes, rotated frames and
/// a non-trivial intrinsic curvature; anisotropic cross-section (`EI₁ ≠ EI₂`).
fn scrambled_rod(seed: u64, angle: f64, shift: f64) -> DiscreteCosseratRod {
    let props = CosseratRodProperties::rectangular(0.03, 0.02, 2.0e5, 7.0e4, 1100.0);
    let mut rod = DiscreteCosseratRod::straight(4, 0.8, props);
    let mut state = seed;
    for node in rod.nodes.iter_mut() {
        node.is_pinned = false;
        for p in node.position.iter_mut() {
            *p += shift * noise(&mut state);
        }
        let axis = [noise(&mut state), noise(&mut state), noise(&mut state)];
        node.frame = node
            .frame
            .rotate(rodrigues(axis, angle * noise(&mut state)));
    }
    for kappa in rod.intrinsic_kappa.iter_mut() {
        *kappa = [0.2, -0.1, 0.05];
    }
    rod
}

fn sum3(vectors: &[[f64; 3]]) -> [f64; 3] {
    vectors.iter().fold([0.0; 3], |acc, v| add3(acc, *v))
}

fn max_norm(vectors: &[[f64; 3]]) -> f64 {
    vectors.iter().fold(0.0_f64, |acc, v| acc.max(norm3(*v)))
}

// ── Midpoint frame ────────────────────────────────────────────────────────────

fn orthonormality_error(r: [f64; 9]) -> f64 {
    let rtr = mat3_mul(mat3_transpose(r), r);
    let id = mat3_identity();
    let mut err = 0.0_f64;
    for k in 0..9 {
        err = err.max((rtr[k] - id[k]).abs());
    }
    err
}

#[test]
fn midpoint_rotation_halves_the_relative_rotation() {
    let cases: [([f64; 3], f64, [f64; 3], f64); 5] = [
        (X, 0.0, Y, 0.0),
        (Z, 0.3, Z, 1.1),
        ([1.0, 2.0, -0.5], 0.7, [-0.3, 0.4, 1.0], -1.9),
        // Large relative rotations, including one close to pi.
        ([0.2, 1.0, 0.1], 0.1, [0.2, 1.0, 0.1], 3.0),
        ([1.0, 1.0, 1.0], -2.5, [0.0, 0.3, -1.0], 0.6),
    ];
    for (axis_a, angle_a, axis_b, angle_b) in cases {
        let ra = rodrigues(axis_a, angle_a);
        let rb = rodrigues(axis_b, angle_b);
        let rm = midpoint_rotation(ra, rb);
        assert!(orthonormality_error(rm) < 1e-12, "R_m must be a rotation");
        // R_aᵀ R_m == R_mᵀ R_b: the midpoint is equally far from both ends.
        let first_half = mat3_mul(mat3_transpose(ra), rm);
        let second_half = mat3_mul(mat3_transpose(rm), rb);
        for k in 0..9 {
            assert!(
                (first_half[k] - second_half[k]).abs() < 1e-12,
                "case {axis_a:?}/{angle_a}: entry {k}: {} vs {}",
                first_half[k],
                second_half[k]
            );
        }
        // And twice the half rotation is the full relative rotation.
        let relative = mat3_mul(mat3_transpose(ra), rb);
        let doubled = mat3_mul(first_half, second_half);
        for k in 0..9 {
            assert!((doubled[k] - relative[k]).abs() < 1e-12);
        }
    }
}

#[test]
fn midpoint_of_rotations_about_one_axis_is_the_mean_angle() {
    let ra = rodrigues(Y, 0.5);
    let rb = rodrigues(Y, -0.3);
    let rm = midpoint_rotation(ra, rb);
    let expected = rodrigues(Y, 0.1);
    for k in 0..9 {
        assert!((rm[k] - expected[k]).abs() < 1e-12, "entry {k}");
    }
    // The midpoint frame of identical frames is that frame.
    let frame = DirectorTriad::canonical().rotate(rodrigues([1.0, 1.0, 0.0], 0.8));
    let mid = midpoint_frame(&frame, &frame);
    for (m, f) in [(mid.d1, frame.d1), (mid.d2, frame.d2), (mid.d3, frame.d3)] {
        for d in 0..3 {
            assert!((m[d] - f[d]).abs() < 1e-14);
        }
    }
}

// The Darboux vector is the rotation vector over the segment length, exact for a
// frame that turns uniformly: R_b = R_a exp(ds [κ]×) gives back κ.
#[test]
fn darboux_vector_is_exact_for_constant_curvature() {
    let kappa = [1.3, -0.7, 2.1];
    let ds = 0.4;
    let phi = scale3(kappa, ds);
    let frame_a = DirectorTriad::canonical().rotate(rodrigues([0.3, -1.0, 0.5], 0.9));
    let rb = mat3_mul(frame_a.to_rotation_matrix(), rodrigues(phi, norm3(phi)));
    let frame_b = DirectorTriad::from_rotation_matrix(rb);
    let recovered = darboux_vector(&frame_a, &frame_b, ds);
    for d in 0..3 {
        assert!(
            (recovered[d] - kappa[d]).abs() < 1e-12,
            "component {d}: {recovered:?} vs {kappa:?}"
        );
    }
}

// ── SO(3) helpers ─────────────────────────────────────────────────────────────

use super::coupling::{exp_rotation, left_jacobian, left_jacobian_inverse, rotation_vector};

/// Angles that straddle the Taylor-series threshold of the Jacobians and reach pi.
const ANGLES: [f64; 12] = [
    0.0,
    1.0e-12,
    1.0e-6,
    5.0e-3,
    9.99e-3,
    1.01e-2,
    0.3,
    1.0,
    2.5,
    3.1,
    std::f64::consts::PI - 1.0e-3,
    std::f64::consts::PI,
];

#[test]
fn rotation_vector_inverts_the_exponential_map() {
    let axes = [[1.0, 0.0, 0.0], [0.3, -1.0, 0.5], [-0.2, 0.1, 1.0]];
    for axis in axes {
        let unit = scale3(axis, 1.0 / norm3(axis));
        for &angle in &ANGLES[..ANGLES.len() - 1] {
            let phi = scale3(unit, angle);
            let recovered = rotation_vector(exp_rotation(phi));
            for d in 0..3 {
                assert!(
                    (recovered[d] - phi[d]).abs() < 1e-9,
                    "angle {angle}: {recovered:?} vs {phi:?}"
                );
            }
        }
    }
    // pi: the sign of the axis is ambiguous, the rotation is not.
    let half_turn = rodrigues(Z, std::f64::consts::PI);
    let phi = rotation_vector(half_turn);
    assert!((norm3(phi) - std::f64::consts::PI).abs() < 1e-9, "{phi:?}");
}

#[test]
fn left_jacobian_inverse_is_the_inverse_of_the_left_jacobian() {
    for axis in [[1.0, 0.0, 0.0], [0.4, 0.7, -0.5]] {
        let unit = scale3(axis, 1.0 / norm3(axis));
        for &angle in &ANGLES {
            let phi = scale3(unit, angle);
            let product = mat3_mul(left_jacobian(phi), left_jacobian_inverse(phi));
            let id = mat3_identity();
            for k in 0..9 {
                assert!(
                    (product[k] - id[k]).abs() < 1e-10,
                    "angle {angle}: entry {k} = {}",
                    product[k]
                );
            }
        }
    }
}

// The series and the closed forms agree across the threshold (no jump).
#[test]
fn jacobians_are_continuous_across_the_series_threshold() {
    let unit = scale3([0.3, -0.5, 0.8], 1.0 / norm3([0.3, -0.5, 0.8]));
    let below = scale3(unit, 1.0e-2 - 1.0e-9);
    let above = scale3(unit, 1.0e-2 + 1.0e-9);
    for (a, b) in [
        (left_jacobian(below), left_jacobian(above)),
        (left_jacobian_inverse(below), left_jacobian_inverse(above)),
    ] {
        for k in 0..9 {
            assert!(
                (a[k] - b[k]).abs() < 1e-8,
                "entry {k}: {} vs {}",
                a[k],
                b[k]
            );
        }
    }
}

// Defining property of the left Jacobian: exp(φ + δ) = exp(J(φ) δ) exp(φ), and of
// its inverse: log(exp(δ) exp(φ)) = φ + J⁻¹(φ) δ, to first order in δ.
#[test]
fn left_jacobian_relates_rotation_vector_and_left_perturbation() {
    let delta = [2.0e-7, -3.0e-7, 1.0e-7];
    for phi in [[0.2, -0.4, 0.1], [1.1, 0.7, -0.9], [-2.0, 1.0, 0.5]] {
        let direct = exp_rotation(add3(phi, delta));
        let jacobian = mat3_mul_vec(left_jacobian(phi), delta);
        let via_jacobian = mat3_mul(exp_rotation(jacobian), exp_rotation(phi));
        for k in 0..9 {
            assert!(
                (direct[k] - via_jacobian[k]).abs() < 1e-12,
                "phi {phi:?} entry {k}: {} vs {}",
                direct[k],
                via_jacobian[k]
            );
        }
        let perturbed = mat3_mul(exp_rotation(delta), exp_rotation(phi));
        let log = rotation_vector(perturbed);
        let predicted = add3(phi, mat3_mul_vec(left_jacobian_inverse(phi), delta));
        for d in 0..3 {
            assert!(
                (log[d] - predicted[d]).abs() < 1e-12,
                "phi {phi:?} axis {d}: {} vs {}",
                log[d],
                predicted[d]
            );
        }
    }
}

// ── Stretch / shear forces ───────────────────────────────────────────────────

#[test]
fn stretch_shear_strain_is_measured_in_the_midpoint_frame() {
    let props = soft_props();
    let mut rod = DiscreteCosseratRod::straight(1, 0.2, props);
    // The end frames lean to opposite sides of the tangent; their midpoint is the
    // canonical frame, so the straight segment is free of shear.
    rotate_frame(&mut rod, 0, Y, 0.2);
    rotate_frame(&mut rod, 1, Y, -0.2);
    let gamma = rod.stretch_shear(0);
    for (k, g) in gamma.iter().enumerate() {
        assert!(g.abs() < 1e-12, "gamma[{k}] = {g}");
    }
    assert!(norm3(rod.segment_force(0)) < 1e-9);
}

#[test]
fn stretched_segment_draws_its_nodes_together() {
    let props = soft_props();
    let mut rod = DiscreteCosseratRod::straight(1, 0.2, props);
    rod.nodes[1].position = [0.0, 0.0, 0.22];
    let forces = rod.elastic_forces();
    let axial = props.ea * 0.1;
    assert!(
        (forces[0][2] - axial).abs() < 1e-9 * axial,
        "{:?}",
        forces[0]
    );
    assert!(
        (forces[1][2] + axial).abs() < 1e-9 * axial,
        "{:?}",
        forces[1]
    );
    assert!(norm3(sum3(&forces)) < 1e-9 * axial);
    // Stretch alone loads no frame.
    for couple in rod.shear_couple_moments() {
        assert!(norm3(couple) < 1e-9 * axial);
    }
}

#[test]
fn sheared_segment_pushes_nodes_and_turns_frames_back() {
    let props = soft_props();
    let mut rod = DiscreteCosseratRod::straight(1, 0.2, props);
    // The tangent leans towards +d1: gamma_1 = 0.05.
    rod.nodes[1].position = [0.01, 0.0, 0.2];
    let gamma = rod.stretch_shear(0);
    assert!((gamma[0] - 0.05).abs() < 1e-12 && gamma[1].abs() < 1e-12);
    let forces = rod.elastic_forces();
    assert!(
        forces[0][0] > 0.0,
        "node 0 is pushed along +d1: {:?}",
        forces[0]
    );
    assert!(
        forces[1][0] < 0.0,
        "node 1 is pushed along -d1: {:?}",
        forces[1]
    );
    // The directors turn about +d2 (d3 x d1 = d2) towards the tangent, on both nodes.
    let couples = rod.shear_couple_moments();
    assert!(couples[0][1] > 0.0 && couples[1][1] > 0.0, "{couples:?}");
    for d in 0..3 {
        assert!((couples[0][d] - couples[1][d]).abs() < 1e-12 * couples[0][1]);
    }
}

#[test]
fn elastic_loads_split_into_forces_joint_moments_and_couples() {
    let rod = scrambled_rod(0x2545_F491_4F6C_DD1D, 0.15, 0.02);
    let (forces, moments) = rod.elastic_loads();
    let joint = rod.elastic_moments();
    let couples = rod.shear_couple_moments();
    let stretch = rod.elastic_forces();
    for node in 0..rod.nodes.len() {
        for d in 0..3 {
            assert!((forces[node][d] - stretch[node][d]).abs() < 1e-12);
            let sum = joint[node][d] + couples[node][d];
            assert!((moments[node][d] - sum).abs() < 1e-12 * (1.0 + sum.abs()));
        }
    }
    assert!(max_norm(&forces) > 1.0 && max_norm(&moments) > 1e-3);
    // The joint moments alone are pairs of action and reaction.
    assert!(norm3(sum3(&joint)) < 1e-12 * max_norm(&joint));
}

// ── Momentum balance ─────────────────────────────────────────────────────────

fn total_force_and_torque(rod: &DiscreteCosseratRod) -> ([f64; 3], [f64; 3], f64) {
    let (forces, moments) = rod.elastic_loads();
    let scale = max_norm(&forces).max(max_norm(&moments));
    let mut net_torque = sum3(&moments);
    for (node, force) in rod.nodes.iter().zip(forces.iter()) {
        net_torque = add3(net_torque, cross3(node.position, *force));
    }
    (sum3(&forces), net_torque, scale)
}

// A free rod in an arbitrary configuration: the elastic loads add up to no force and
// no torque about the origin.  This is the algebraic form of linear and angular
// momentum conservation.
#[test]
fn elastic_loads_carry_no_net_force_or_torque() {
    for seed in [
        0x9E37_79B9_7F4A_7C15_u64,
        0xD1B5_4A32_D192_ED03,
        0x2545_F491_4F6C_DD1D,
    ] {
        let rod = scrambled_rod(seed, 0.4, 0.05);
        let (net_force, net_torque, scale) = total_force_and_torque(&rod);
        assert!(scale > 1.0, "the configuration must load the rod: {scale}");
        assert!(norm3(net_force) < 1e-9 * scale, "net force {net_force:?}");
        assert!(
            norm3(net_torque) < 1e-9 * scale,
            "net torque {net_torque:?}"
        );
    }
}

// ── Consistency with the stored energy ───────────────────────────────────────

fn energy_with_displaced_node(rod: &DiscreteCosseratRod, node: usize, d: usize, h: f64) -> f64 {
    let mut copy = rod.clone();
    copy.nodes[node].position[d] += h;
    copy.elastic_energy()
}

fn energy_with_rotated_node(rod: &DiscreteCosseratRod, node: usize, axis: [f64; 3], h: f64) -> f64 {
    let mut copy = rod.clone();
    rotate_frame(&mut copy, node, axis, h);
    copy.elastic_energy()
}

/// The nodal loads of `rod` against central finite differences of its elastic
/// energy: forces are `-dW/dx` and moments `-dW/dθ` (global rotation of the frame).
fn check_loads_are_the_energy_gradient(rod: &DiscreteCosseratRod, tolerance: f64) {
    let (forces, moments) = rod.elastic_loads();
    let force_scale = max_norm(&forces);
    let moment_scale = max_norm(&moments);
    assert!(
        force_scale > 0.0 || moment_scale > 0.0,
        "the configuration must load the rod"
    );
    let h = 1e-6;
    for node in 0..rod.nodes.len() {
        for (d, component) in forces[node].iter().enumerate() {
            let slope = (energy_with_displaced_node(rod, node, d, h)
                - energy_with_displaced_node(rod, node, d, -h))
                / (2.0 * h);
            assert!(
                (component + slope).abs() <= tolerance * force_scale.max(1e-12),
                "node {node} axis {d}: force {component} vs -dW/dx {} (scale {force_scale})",
                -slope
            );
        }
        for (d, axis) in [X, Y, Z].into_iter().enumerate() {
            let slope = (energy_with_rotated_node(rod, node, axis, h)
                - energy_with_rotated_node(rod, node, axis, -h))
                / (2.0 * h);
            assert!(
                (moments[node][d] + slope).abs() <= tolerance * moment_scale.max(1e-12),
                "node {node} axis {d}: moment {} vs -dW/dθ {} (scale {moment_scale})",
                moments[node][d],
                -slope
            );
        }
    }
}

// The nodal loads are the exact negative gradient of the elastic energy with respect
// to the node positions and to rotations of the frames, also for large relative
// rotations: the joint moments carry the SO(3) Jacobian of the Darboux vector and
// the couples the exact derivative of the midpoint frame.
#[test]
fn loads_are_the_negative_energy_gradient() {
    for (seed, angle, shift) in [
        (0x9E37_79B9_7F4A_7C15_u64, 0.05, 0.005),
        (0xD1B5_4A32_D192_ED03, 0.4, 0.03),
        (0x2545_F491_4F6C_DD1D, 1.2, 0.05),
    ] {
        check_loads_are_the_energy_gradient(&scrambled_rod(seed, angle, shift), 1e-6);
    }
}

// Each energy on its own: the bending / twist moments of an anisotropic section with
// intrinsic curvature (no stretch or shear stiffness), and the couples of the
// stretch / shear forces (no bending or twist stiffness).  The bending moments are
// two orders of magnitude smaller than the couples, so only the isolated checks
// resolve them.
#[test]
fn bending_and_shear_loads_are_each_an_exact_energy_gradient() {
    for (seed, angle, shift) in [
        (0xD1B5_4A32_D192_ED03_u64, 0.3, 0.02),
        (0x9E37_79B9_7F4A_7C15, 1.5, 0.04),
    ] {
        let mut bending = scrambled_rod(seed, angle, shift);
        bending.props.ea = 0.0;
        bending.props.ga1 = 0.0;
        bending.props.ga2 = 0.0;
        assert!(bending.stretch_shear_energy() == 0.0);
        assert!(bending.bending_twist_energy() > 0.0);
        check_loads_are_the_energy_gradient(&bending, 1e-6);

        let mut shear = scrambled_rod(seed, angle, shift);
        shear.props.ei1 = 0.0;
        shear.props.ei2 = 0.0;
        shear.props.gj = 0.0;
        assert!(shear.bending_twist_energy() == 0.0);
        assert!(shear.stretch_shear_energy() > 0.0);
        check_loads_are_the_energy_gradient(&shear, 1e-6);
    }
}

#[test]
fn elastic_energy_is_the_sum_of_its_parts_and_zero_at_rest() {
    let rod = scrambled_rod(0x2545_F491_4F6C_DD1D, 0.3, 0.02);
    let total = rod.elastic_energy();
    assert!(rod.stretch_shear_energy() > 0.0 && rod.bending_twist_energy() > 0.0);
    assert!(
        (total - rod.stretch_shear_energy() - rod.bending_twist_energy()).abs() < 1e-15 * total
    );
    let straight = DiscreteCosseratRod::straight(5, 1.0, soft_props());
    assert!(straight.elastic_energy().abs() < 1e-18);
    let (forces, moments) = straight.elastic_loads();
    assert!(max_norm(&forces) < 1e-9 && max_norm(&moments) < 1e-12);
}

// ── Bending reaches the positions ─────────────────────────────────────────────

/// Relax a rod under a constant tip force with damped explicit dynamics.
fn relax(rod: &mut DiscreteCosseratRod, tip_force: [f64; 3], dt: f64, steps: usize, damping: f64) {
    let tip = rod.nodes.len() - 1;
    for _ in 0..steps {
        rod.nodes[tip].apply_force(tip_force);
        rod.step(dt, damping);
    }
}

// The tip deflection of a cantilever follows the Timoshenko beam formula
// delta = P L^3 / (3 EI) + P L / (G A) with the bending stiffness of the load axis
// (EI₂ = 4 EI₁ for the 0.04 x 0.02 section), and the bending part scales with 1/EI.
// The shear modulus is low enough for the shear deflection to be 9 - 44 % of the
// total, so a wrong shear coupling cannot hide behind the bending part.
#[test]
fn cantilever_sag_follows_timoshenko_theory_and_scales_with_ei() {
    let length = 0.5;
    let load = 1.0e-2;
    let dt = 2.0e-4;
    // (bending part, shear part, total) of the four runs.
    let mut runs: Vec<(f64, f64)> = Vec::new();
    for scale in [1.0, 2.0] {
        let mut props =
            CosseratRodProperties::rectangular(0.04, 0.02, 1.0e7 * scale, 4.0e4, 1000.0);
        props.rho_a = 1000.0 * 0.04 * 0.02;
        for (axis, direction) in [(0usize, X), (1usize, Y)] {
            let mut rod = DiscreteCosseratRod::straight(10, length, props);
            let ei = if axis == 0 { props.ei2 } else { props.ei1 };
            let bending = load * length.powi(3) / (3.0 * ei);
            let shear = load * length / props.ga1;
            assert!(shear > 0.08 * (bending + shear), "shear must matter");
            let expected = bending + shear;
            relax(&mut rod, scale3(direction, load), dt, 12_000, 8.0);
            let tip = rod.nodes.len() - 1;
            let sag = rod.nodes[tip].position[axis];
            assert!(
                (sag - expected).abs() < 0.01 * expected,
                "scale {scale}, load axis {axis}: tip sag {sag:e}, Timoshenko {expected:e}"
            );
            let other = rod.nodes[tip].position[1 - axis];
            assert!(other.abs() < 1e-3 * sag, "cross deflection {other:e}");
            runs.push((sag - shear, shear));
        }
    }
    // Doubled E (hence EI) halves the bending part; the softer axis bends 4x more.
    assert!(
        (runs[0].0 / runs[2].0 - 2.0).abs() < 0.04,
        "d1 load: {runs:?}"
    );
    assert!(
        (runs[1].0 / runs[3].0 - 2.0).abs() < 0.04,
        "d2 load: {runs:?}"
    );
    assert!(
        (runs[1].0 / runs[0].0 - 4.0).abs() < 0.1,
        "EI2/EI1 = 4: {runs:?}"
    );
}

/// Positions of the polygon that a rod of `segments` segments of length `ds` and
/// uniform relative rotation `θ` about `d₁` per segment relaxes to: each segment
/// follows the midpoint director, `x_{k+1} − x_k = ds · R_m e₃`.
fn arc_polygon(segments: usize, ds: f64, theta: f64) -> Vec<[f64; 3]> {
    let mut points = vec![[0.0; 3]];
    for i in 0..segments {
        let angle = (i as f64 + 0.5) * theta;
        let last = points[i];
        points.push([
            last[0],
            last[1] - ds * angle.sin(),
            last[2] + ds * angle.cos(),
        ]);
    }
    points
}

// A clamped rod with an intrinsic curvature about d1 relaxes to a circular arc: the
// bending stiffness moves the nodes.  (Without the stretch / shear coupling the
// frames turned but the positions stayed on the straight line.)
#[test]
fn intrinsically_curved_rod_relaxes_to_a_circular_arc() {
    let (segments, length, kappa0) = (8usize, 0.2, 3.0);
    let mut rod = DiscreteCosseratRod::straight(segments, length, soft_props());
    for kappa in rod.intrinsic_kappa.iter_mut() {
        *kappa = [kappa0, 0.0, 0.0];
    }
    let ds = length / segments as f64;
    let theta = kappa0 * ds;
    for _ in 0..25_000 {
        rod.step(1.0e-4, 30.0);
    }
    // Every joint has relaxed to the rest curvature and no strain is left ...
    for i in 0..segments {
        let kappa = rod.darboux(i);
        assert!(
            (kappa[0] - kappa0).abs() < 1e-4 * kappa0,
            "segment {i}: {kappa:?}"
        );
        assert!(kappa[1].abs() < 1e-4 && kappa[2].abs() < 1e-4);
        let gamma = rod.stretch_shear(i);
        assert!(norm3(gamma) < 1e-5, "segment {i}: gamma {gamma:?}");
    }
    // ... the nodes sit on the polygon of the midpoint directors, which differs from
    // the exact circle of radius 1/κ₀ only by the chord error (κ₀ ds)² / 24.
    let polygon = arc_polygon(segments, ds, theta);
    for (i, expected) in polygon.iter().enumerate() {
        let position = rod.nodes[i].position;
        assert!(
            norm3(sub3(position, *expected)) < 2e-4 * length,
            "node {i}: {position:?} vs {expected:?}"
        );
    }
    let tip = rod.nodes[segments].position;
    let arc_y = ((kappa0 * length).cos() - 1.0) / kappa0;
    let arc_z = (kappa0 * length).sin() / kappa0;
    assert!(tip[1] < -0.1 * length, "the rod must bend: tip {tip:?}");
    assert!(
        (tip[1] - arc_y).abs() < 2e-3 * length && (tip[2] - arc_z).abs() < 2e-3 * length,
        "tip {tip:?} vs circle ({arc_y}, {arc_z})"
    );
    assert!(tip[0].abs() < 1e-6, "the arc stays in the plane of bending");
}

// Tuning the intrinsic curvature of an `ActiveCosseratRod` bends the rod.
#[test]
fn active_rod_bends_towards_its_target_curvature() {
    let (segments, length, kappa0) = (8usize, 0.2, 3.0);
    let rod = DiscreteCosseratRod::straight(segments, length, soft_props());
    let mut active = ActiveCosseratRod::new(rod, 0.05);
    for i in 0..segments {
        active.set_target_curvature(i, [kappa0, 0.0, 0.0]);
    }
    let straight_tip = active.rod.nodes[segments].position;
    for _ in 0..30_000 {
        active.step(1.0e-4, [0.0; 3], 30.0);
    }
    let tip = active.rod.nodes[segments].position;
    let expected = ((kappa0 * length).cos() - 1.0) / kappa0;
    assert!(
        (tip[1] - expected).abs() < 5e-3 * length,
        "tip y {} vs arc {expected} (straight {})",
        tip[1],
        straight_tip[1]
    );
    assert!(
        tip[2] < straight_tip[2],
        "the chord shortens as the rod bends"
    );
}

// Interior frame rotated against the tangent: bending, shear and stretch all act,
// the whole rod (clamped at node 0) returns to the straight configuration.
#[test]
fn interior_rotation_relaxes_to_a_straight_rod() {
    let mut rod = DiscreteCosseratRod::straight(4, 0.4, soft_props());
    rotate_frame(&mut rod, 2, Y, 0.05);
    let initial_energy = rod.elastic_energy();
    assert!(rod.stretch_shear_energy() > 0.0 && rod.bending_twist_energy() > 0.0);
    for _ in 0..20_000 {
        rod.step(1.0e-4, 10.0);
    }
    assert!(
        rod.elastic_energy() < 1e-6 * initial_energy,
        "elastic energy {initial_energy} -> {}",
        rod.elastic_energy()
    );
    for (i, node) in rod.nodes.iter().enumerate() {
        let expected_z = 0.1 * i as f64;
        assert!(node.position[0].abs() < 1e-5 && node.position[1].abs() < 1e-5);
        assert!((node.position[2] - expected_z).abs() < 1e-5);
    }
}

// ── Free rod dynamics ───────────────────────────────────────────────────────

/// `(linear momentum, angular momentum about the origin, kinetic energy)` of a rod
/// with nodal masses `ρA ds` and rotary inertias `diag(ρI₁, ρI₂, ρJ) ds`.
fn momenta_and_kinetic_energy(rod: &DiscreteCosseratRod) -> ([f64; 3], [f64; 3], f64) {
    let node_length = rod.rest_lengths.iter().sum::<f64>() / rod.num_segments() as f64;
    let mass = rod.props.rho_a * node_length;
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    let mut kinetic = 0.0;
    for node in &rod.nodes {
        linear = add3(linear, scale3(node.velocity, mass));
        angular = add3(angular, scale3(cross3(node.position, node.velocity), mass));
        kinetic += 0.5 * mass * dot3(node.velocity, node.velocity);
        let rotation = node.frame.to_rotation_matrix();
        let omega_body = mat3_mul_vec(mat3_transpose(rotation), node.omega);
        let spin_body = [
            rod.props.rho_i1 * node_length * omega_body[0],
            rod.props.rho_i2 * node_length * omega_body[1],
            rod.props.rho_j * node_length * omega_body[2],
        ];
        angular = add3(angular, mat3_mul_vec(rotation, spin_body));
        kinetic += 0.5 * dot3(omega_body, spin_body);
    }
    (linear, angular, kinetic)
}

/// Largest drift of the linear momentum, of the angular momentum and of the total
/// energy over `time` seconds of undamped free motion of a scrambled rod.
fn free_flight_drift(dt: f64, time: f64) -> (f64, f64, f64, f64) {
    let mut rod = scrambled_rod(0x9E37_79B9_7F4A_7C15, 0.1, 0.01);
    for kappa in rod.intrinsic_kappa.iter_mut() {
        *kappa = [0.0; 3];
    }
    let (_, _, kinetic0) = momenta_and_kinetic_energy(&rod);
    let energy0 = rod.elastic_energy() + kinetic0;
    let mut max_linear = 0.0_f64;
    let mut max_angular = 0.0_f64;
    let mut max_energy = 0.0_f64;
    let mut momentum_scale = 0.0_f64;
    let mass = rod.props.rho_a * rod.rest_lengths.iter().sum::<f64>() / rod.num_segments() as f64;
    for _ in 0..(time / dt).round() as usize {
        rod.step(dt, 0.0);
        let (linear, angular, kinetic) = momenta_and_kinetic_energy(&rod);
        max_linear = max_linear.max(norm3(linear));
        max_angular = max_angular.max(norm3(angular));
        max_energy = max_energy.max(((rod.elastic_energy() + kinetic - energy0) / energy0).abs());
        for node in &rod.nodes {
            momentum_scale = momentum_scale.max(mass * norm3(node.velocity));
        }
    }
    (max_linear, max_angular, max_energy, momentum_scale)
}

// A free rod started at rest conserves its (zero) linear and angular momentum: the
// linear momentum to rounding, because the forces add up to zero exactly; the
// angular momentum up to the first-order integration error of the explicit scheme,
// which halves with the time step.  The total energy oscillates with the same
// first-order accuracy instead of drifting.
#[test]
fn free_rod_conserves_momentum_and_energy_in_time() {
    let (linear_coarse, angular_coarse, energy_coarse, scale) = free_flight_drift(4.0e-5, 0.12);
    let (linear_fine, angular_fine, energy_fine, _) = free_flight_drift(2.0e-5, 0.12);
    assert!(scale > 1e-3, "the rod must move: {scale}");
    assert!(linear_coarse < 1e-12 * scale && linear_fine < 1e-12 * scale);
    assert!(angular_coarse < 1e-4 * scale, "{angular_coarse}");
    assert!(angular_fine < 1e-4 * scale, "{angular_fine}");
    let angular_ratio = angular_coarse / angular_fine;
    assert!(
        (1.5..2.6).contains(&angular_ratio),
        "angular momentum drift must be first order in dt: {angular_coarse} / {angular_fine}"
    );
    assert!(energy_coarse < 3e-2 && energy_fine < 2e-2);
    let energy_ratio = energy_coarse / energy_fine;
    assert!(
        (1.5..2.6).contains(&energy_ratio),
        "energy error must be first order in dt: {energy_coarse} / {energy_fine}"
    );
}
