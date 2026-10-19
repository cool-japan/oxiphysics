// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the neo-Hookean elasticity, the deformed-geometry Maxwell forces and
//! the explicit integrator of the [`ElectroacitiveFemBody`].
//!
//! The earlier body integrated only gravity and damping: the Maxwell force was added
//! to the velocity without the time step, the shear and bulk moduli were never used,
//! and every node of an element received the same force.

use super::*;

const MU: f64 = 1.0e5;
const KAPPA: f64 = 4.0e5;

type Vec3 = [f64; 3];

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

fn mat_vec(m: &[[f64; 3]; 3], v: Vec3) -> Vec3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

/// Rotation matrix about the unit axis `axis` by `angle` (Rodrigues).
fn rotation(axis: Vec3, angle: f64) -> [[f64; 3]; 3] {
    let n = norm(axis);
    let (x, y, z) = (axis[0] / n, axis[1] / n, axis[2] / n);
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [c + x * x * t, x * y * t - z * s, x * z * t + y * s],
        [y * x * t + z * s, c + y * y * t, y * z * t - x * s],
        [z * x * t - y * s, z * y * t + x * s, c + z * z * t],
    ]
}

/// A skewed millimetre-scale tetrahedron (neither `Dm` nor `Dm⁻¹` is symmetric).
const REST: [Vec3; 4] = [
    [0.0, 0.0, 0.0],
    [1.0e-3, 2.0e-4, 1.0e-4],
    [3.0e-4, 1.1e-3, 0.0],
    [2.0e-4, 4.0e-4, 9.0e-4],
];

fn element_and_nodes() -> (ElectroacitiveTetraElement, Vec<ElectroacitiveFemNode>) {
    let nodes: Vec<ElectroacitiveFemNode> = REST
        .iter()
        .map(|&p| ElectroacitiveFemNode::new(p, 1.0e-6))
        .collect();
    let element = ElectroacitiveTetraElement::new([0, 1, 2, 3], &REST, MU, KAPPA, EPSILON_R_VHB);
    (element, nodes)
}

/// Move the nodes to `x = F X + t`.
fn deform(nodes: &mut [ElectroacitiveFemNode], f: &[[f64; 3]; 3], translation: Vec3) {
    for (node, rest) in nodes.iter_mut().zip(REST.iter()) {
        let image = mat_vec(f, *rest);
        for d in 0..3 {
            node.position[d] = image[d] + translation[d];
        }
    }
}

/// A general deformation gradient with stretch, shear and rotation (J ≈ 1.07).
const GENERAL_F: [[f64; 3]; 3] = [[1.08, 0.12, -0.05], [0.04, 0.93, 0.10], [-0.07, 0.03, 1.09]];

fn max_abs(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .flat_map(|f| f.iter())
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

fn sum_forces(forces: &[[f64; 3]]) -> Vec3 {
    let mut sum = [0.0; 3];
    for f in forces {
        for d in 0..3 {
            sum[d] += f[d];
        }
    }
    sum
}

// ── Neo-Hookean elasticity ───────────────────────────────────────────────────

#[test]
fn rest_configuration_is_force_free() {
    let (element, nodes) = element_and_nodes();
    let f = element.deformation_gradient(&nodes);
    for (i, row) in f.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((value - want).abs() < 1e-12, "F[{i}][{j}] = {value}");
        }
    }
    assert!(max_abs(&element.elastic_nodal_forces(&nodes)) < 1e-9);
    assert!(element.elastic_energy(&nodes).abs() < 1e-18);
}

// A finite rigid rotation and translation costs no energy and exerts no force.
#[test]
fn rigid_motion_exerts_no_elastic_force() {
    let (element, mut nodes) = element_and_nodes();
    let r = rotation([0.3, -1.0, 0.6], 1.1);
    deform(&mut nodes, &r, [2.0e-3, -1.0e-3, 5.0e-4]);
    let scale = MU * element.volume_0 / 1.0e-3;
    let forces = element.elastic_nodal_forces(&nodes);
    assert!(
        max_abs(&forces) < 1e-9 * scale,
        "forces {forces:?} (scale {scale:e})"
    );
    assert!(element.elastic_energy(&nodes).abs() < 1e-9 * MU * element.volume_0);
}

// Newton's third law and the balance of angular momentum for a general deformation.
#[test]
fn elastic_forces_have_no_net_force_or_torque() {
    let (element, mut nodes) = element_and_nodes();
    deform(&mut nodes, &GENERAL_F, [1.0e-3, 0.0, -2.0e-3]);
    let forces = element.elastic_nodal_forces(&nodes);
    let scale = max_abs(&forces);
    assert!(
        scale > 1e-6,
        "the deformation must load the element: {scale:e}"
    );
    let net = sum_forces(&forces);
    let mut torque = [0.0; 3];
    for (node, f) in nodes.iter().zip(forces.iter()) {
        let x = node.position;
        torque[0] += x[1] * f[2] - x[2] * f[1];
        torque[1] += x[2] * f[0] - x[0] * f[2];
        torque[2] += x[0] * f[1] - x[1] * f[0];
    }
    for d in 0..3 {
        assert!(net[d].abs() < 1e-9 * scale, "net force {net:?}");
        assert!(
            torque[d].abs() < 1e-9 * scale * 1e-3,
            "net torque {torque:?}"
        );
    }
}

fn displaced_elastic_energy(
    element: &ElectroacitiveTetraElement,
    nodes: &[ElectroacitiveFemNode],
    node: usize,
    d: usize,
    h: f64,
) -> f64 {
    let mut copy = nodes.to_vec();
    copy[node].position[d] += h;
    element.elastic_energy(&copy)
}

// The forces are the exact negative gradient of the strain energy.
#[test]
fn elastic_forces_are_the_negative_energy_gradient() {
    let (element, mut nodes) = element_and_nodes();
    deform(&mut nodes, &GENERAL_F, [0.0; 3]);
    let forces = element.elastic_nodal_forces(&nodes);
    let scale = max_abs(&forces);
    let h = 1.0e-9;
    for (node, force) in forces.iter().enumerate() {
        for (d, component) in force.iter().enumerate() {
            let slope = (displaced_elastic_energy(&element, &nodes, node, d, h)
                - displaced_elastic_energy(&element, &nodes, node, d, -h))
                / (2.0 * h);
            assert!(
                (component + slope).abs() < 1e-6 * scale,
                "node {node} axis {d}: {component} vs -dΨ/dx {}",
                -slope
            );
        }
    }
}

// For a homogeneous deformation x = F X the first moment of the nodal forces
// Σ_k f_k ⊗ X_k equals -V₀ P, with P = μ (F - F⁻ᵀ) + λ ln J F⁻ᵀ.
#[test]
fn homogeneous_deformation_recovers_the_neo_hookean_stress() {
    let (element, mut nodes) = element_and_nodes();
    deform(&mut nodes, &GENERAL_F, [0.0; 3]);
    let forces = element.elastic_nodal_forces(&nodes);
    // F⁻ᵀ from an independent inversion.
    let inverse = invert_3x3(GENERAL_F);
    let det = GENERAL_F[0][0]
        * (GENERAL_F[1][1] * GENERAL_F[2][2] - GENERAL_F[1][2] * GENERAL_F[2][1])
        - GENERAL_F[0][1] * (GENERAL_F[1][0] * GENERAL_F[2][2] - GENERAL_F[1][2] * GENERAL_F[2][0])
        + GENERAL_F[0][2] * (GENERAL_F[1][0] * GENERAL_F[2][1] - GENERAL_F[1][1] * GENERAL_F[2][0]);
    let lambda = KAPPA - 2.0 * MU / 3.0;
    let mut scale = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let inverse_transpose = inverse[j][i];
            let p =
                MU * (GENERAL_F[i][j] - inverse_transpose) + lambda * det.ln() * inverse_transpose;
            let mut moment = 0.0;
            for k in 0..4 {
                moment += forces[k][i] * REST[k][j];
            }
            let expected = -element.volume_0 * p;
            scale = scale.max(expected.abs());
            assert!(
                (moment - expected).abs() < 1e-9 * expected.abs().max(1e-3 * scale.max(1e-12)),
                "moment[{i}][{j}] = {moment:e}, -V0 P = {expected:e}"
            );
        }
    }
}

// The stored moduli are the small-strain shear and bulk moduli: hydrostatic
// stretching gives P = 3 κ e I and simple shear P₀₁ = P₁₀ = μ γ.
#[test]
fn stored_moduli_are_the_small_strain_moduli() {
    let (element, _) = element_and_nodes();
    let e = 1.0e-7;
    let hydrostatic = element.first_piola_stress([
        [1.0 + e, 0.0, 0.0],
        [0.0, 1.0 + e, 0.0],
        [0.0, 0.0, 1.0 + e],
    ]);
    for (i, row) in hydrostatic.iter().enumerate() {
        assert!(
            (row[i] - 3.0 * KAPPA * e).abs() < 1e-5 * 3.0 * KAPPA * e,
            "P[{i}][{i}] = {} vs 3 K e = {}",
            row[i],
            3.0 * KAPPA * e
        );
    }
    let gamma = 1.0e-7;
    let shear = element.first_piola_stress([[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    assert!(
        (shear[0][1] - MU * gamma).abs() < 1e-5 * MU * gamma,
        "{shear:?}"
    );
    assert!(
        (shear[1][0] - MU * gamma).abs() < 1e-5 * MU * gamma,
        "{shear:?}"
    );
    assert!(shear[0][0].abs() < 1e-5 * MU * gamma && shear[2][2].abs() < 1e-5 * MU * gamma);
}

// A node pushed through the opposite face (negative volume) must not produce NaN or
// net force, and the element is pushed back: moving the node along the force
// lowers the clamped energy trend towards positive volume.
#[test]
fn inverted_element_has_finite_forces() {
    let (element, mut nodes) = element_and_nodes();
    // Mirror the apex through the plane of the opposite face.
    let mirrored = [REST[3][0], REST[3][1], -REST[3][2]];
    nodes[3].position = mirrored;
    let f = element.deformation_gradient(&nodes);
    let j = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
        - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
        + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
    assert!(j < 0.0, "the element must be inverted: J = {j}");
    let forces = element.elastic_nodal_forces(&nodes);
    assert!(forces.iter().flat_map(|f| f.iter()).all(|v| v.is_finite()));
    let net = sum_forces(&forces);
    assert!(norm(net) < 1e-9 * max_abs(&forces));
    // The apex is pushed back towards positive z (towards positive volume).
    assert!(forces[3][2] > 0.0, "apex force {:?}", forces[3]);
}

// ── Maxwell forces on the deformed geometry ──────────────────────────────────

fn with_potentials(nodes: &mut [ElectroacitiveFemNode], potentials: [f64; 4]) {
    for (node, phi) in nodes.iter_mut().zip(potentials.iter()) {
        node.potential = *phi;
    }
}

#[test]
fn deformed_maxwell_forces_reduce_to_the_reference_form_at_rest() {
    let (element, mut nodes) = element_and_nodes();
    with_potentials(&mut nodes, [0.0, 1.5e3, -8.0e2, 2.2e3]);
    let reference = element.maxwell_nodal_forces(&nodes);
    let deformed = element.maxwell_nodal_forces_deformed(&nodes);
    let scale = max_abs(&reference);
    assert!(scale > 1e-9);
    for k in 0..4 {
        for d in 0..3 {
            assert!(
                (reference[k][d] - deformed[k][d]).abs() < 1e-9 * scale,
                "node {k} axis {d}: {} vs {}",
                reference[k][d],
                deformed[k][d]
            );
        }
    }
}

fn displaced_electric_energy(
    element: &ElectroacitiveTetraElement,
    nodes: &[ElectroacitiveFemNode],
    node: usize,
    d: usize,
    h: f64,
) -> f64 {
    let mut copy = nodes.to_vec();
    copy[node].position[d] += h;
    element.electric_energy(&copy)
}

// At fixed nodal potentials the electric force is +∂W_e/∂x with W_e = ½ ε v |E|²
// (the sources do work), also for a strongly deformed element.
#[test]
fn deformed_maxwell_forces_are_the_electric_energy_gradient() {
    let (element, mut nodes) = element_and_nodes();
    deform(&mut nodes, &GENERAL_F, [0.0; 3]);
    with_potentials(&mut nodes, [0.0, 1.5e3, -8.0e2, 2.2e3]);
    let forces = element.maxwell_nodal_forces_deformed(&nodes);
    let scale = max_abs(&forces);
    assert!(scale > 1e-9);
    assert!(norm(sum_forces(&forces)) < 1e-9 * scale);
    let h = 1.0e-9;
    for (node, force) in forces.iter().enumerate() {
        for (d, component) in force.iter().enumerate() {
            let slope = (displaced_electric_energy(&element, &nodes, node, d, h)
                - displaced_electric_energy(&element, &nodes, node, d, -h))
                / (2.0 * h);
            assert!(
                (component - slope).abs() < 1e-6 * scale,
                "node {node} axis {d}: {component} vs +dW_e/dx {slope}"
            );
        }
    }
}

// A film that thins at fixed potentials sees a stronger field: compressing the
// element by 20 % along the field raises the force by 1 / 0.8² (E² grows by 1.5625
// while v ∇N_z stays constant), which the reference form cannot show.
#[test]
fn thinner_film_feels_a_stronger_maxwell_force() {
    let h = 1.0e-3;
    let rest = [[0.0, 0.0, 0.0], [h, 0.0, 0.0], [0.0, h, 0.0], [0.0, 0.0, h]];
    let element = ElectroacitiveTetraElement::new([0, 1, 2, 3], &rest, MU, KAPPA, EPSILON_R_VHB);
    let mut nodes: Vec<ElectroacitiveFemNode> = rest
        .iter()
        .map(|&p| ElectroacitiveFemNode::new(p, 1.0e-6))
        .collect();
    nodes[3].potential = 2.0e3;
    let thick = element.maxwell_nodal_forces_deformed(&nodes)[3][2];
    nodes[3].position[2] = 0.8 * h;
    let thin = element.maxwell_nodal_forces_deformed(&nodes)[3][2];
    assert!(thick < 0.0 && thin < 0.0, "the top node is pulled down");
    assert!(
        (thin / thick - 1.0 / 0.64).abs() < 1e-9,
        "force ratio {} vs {}",
        thin / thick,
        1.0 / 0.64
    );
    // The reference form is blind to the deformation.
    let reference = element.maxwell_nodal_forces(&nodes)[3][2];
    assert!((reference / thick - 1.0).abs() < 1e-9);
}

#[test]
fn degenerate_element_carries_no_maxwell_force() {
    let (element, mut nodes) = element_and_nodes();
    with_potentials(&mut nodes, [0.0, 1.0e3, 2.0e3, 3.0e3]);
    // Flatten the apex into the opposite face (its centroid): zero volume.
    let mut centroid = [0.0; 3];
    for node in &nodes[..3] {
        for (c, x) in centroid.iter_mut().zip(node.position.iter()) {
            *c += x / 3.0;
        }
    }
    nodes[3].position = centroid;
    let forces = element.maxwell_nodal_forces_deformed(&nodes);
    assert!(forces.iter().flat_map(|f| f.iter()).all(|v| *v == 0.0));
    assert!(element.electric_energy(&nodes) == 0.0);
}

// ── The body integrator ──────────────────────────────────────────────────────

/// A block of `nx × ny × nz` cubic cells of edge `h` with six tetrahedra per cell.
/// The bottom layer is pinned at potential zero, the potential rises linearly to
/// `voltage` at the free top layer.
fn film(
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
    voltage: f64,
    damping: f64,
) -> ElectroacitiveFemBody {
    film_pinned_at(nx, ny, nz, h, voltage, damping, 0, 1.0e7)
}

/// [`film`] with the layer `pinned_layer` (`0` = bottom, `nz` = top) pinned and the
/// given bulk modulus (the default film is nearly incompressible).
fn film_pinned_at(
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
    voltage: f64,
    damping: f64,
    pinned_layer: usize,
    bulk_modulus: f64,
) -> ElectroacitiveFemBody {
    let node_id = |i: usize, j: usize, k: usize| (k * (ny + 1) + j) * (nx + 1) + i;
    let mut positions = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                positions.push([i as f64 * h, j as f64 * h, k as f64 * h]);
            }
        }
    }
    let mut elements = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let c = |di: usize, dj: usize, dk: usize| node_id(i + di, j + dj, k + dk);
                // Kuhn split along the main diagonal of the cell.
                let corners = [
                    [c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(1, 1, 1)],
                    [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(1, 1, 1)],
                    [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 1, 1)],
                    [c(0, 0, 0), c(0, 1, 0), c(0, 1, 1), c(1, 1, 1)],
                    [c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1)],
                    [c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(1, 1, 1)],
                ];
                for indices in corners {
                    elements.push(ElectroacitiveTetraElement::new(
                        indices,
                        &positions,
                        MU,
                        bulk_modulus,
                        EPSILON_R_VHB,
                    ));
                }
            }
        }
    }
    // Lumped masses from the density 1000 kg/m^3.
    let mut masses = vec![0.0_f64; positions.len()];
    for element in &elements {
        for &node in &element.node_indices {
            masses[node] += 1000.0 * element.volume_0 / 4.0;
        }
    }
    let mut nodes = Vec::new();
    for (index, &position) in positions.iter().enumerate() {
        let layer = index / ((nx + 1) * (ny + 1));
        let mut node = if layer == pinned_layer {
            ElectroacitiveFemNode::new_pinned(position)
        } else {
            ElectroacitiveFemNode::new(position, masses[index])
        };
        node.potential = voltage * layer as f64 / nz as f64;
        nodes.push(node);
    }
    ElectroacitiveFemBody::new(nodes, elements, damping)
}

fn top_layer(nx: usize, ny: usize, nz: usize) -> Vec<usize> {
    let per_layer = (nx + 1) * (ny + 1);
    (nz * per_layer..(nz + 1) * per_layer).collect()
}

// The Maxwell and elastic forces of a body carry no net force.
#[test]
fn body_forces_balance() {
    let mut body = film(2, 2, 1, 1.0e-3, 3.0e3, 0.0);
    // Perturb the free nodes so that both force sets are non-trivial.
    for (i, node) in body.nodes.iter_mut().enumerate() {
        if !node.pinned {
            node.position[0] += 1.0e-5 * ((i % 3) as f64 - 1.0);
            node.position[2] -= 2.0e-5 * (i % 2) as f64;
        }
    }
    for forces in [
        body.elastic_forces(),
        body.maxwell_forces(),
        body.net_forces(),
    ] {
        let scale = max_abs(&forces);
        assert!(scale > 1e-6, "{scale:e}");
        assert!(norm(sum_forces(&forces)) < 1e-9 * scale);
    }
}

// `apply_maxwell_forces` / `apply_elastic_forces` only add to the force buffer: they
// neither accelerate nor move the nodes (the old `apply_maxwell_forces` added F / m
// to the velocities without a time step).
#[test]
fn applying_forces_only_accumulates() {
    let body = film(1, 1, 1, 1.0e-3, 2.0e3, 0.0);
    let offset = [1.0, 2.0, 3.0];
    let mut forces = vec![offset; body.nodes.len()];
    body.apply_maxwell_forces(&mut forces);
    body.apply_elastic_forces(&mut forces);
    let maxwell = body.maxwell_forces();
    let elastic = body.elastic_forces();
    let net = body.net_forces();
    assert!(max_abs(&maxwell) > 1e-6, "the field must load the film");
    for i in 0..body.nodes.len() {
        for d in 0..3 {
            let expected = offset[d] + maxwell[i][d] + elastic[i][d];
            assert!((forces[i][d] - expected).abs() < 1e-12 * expected.abs().max(1.0));
            assert!((net[i][d] - maxwell[i][d] - elastic[i][d]).abs() < 1e-12);
        }
        assert_eq!(body.nodes[i].velocity, [0.0; 3]);
    }
}

// Δv scales with dt: one step of dt gives the same velocity change per time as two
// steps of dt/2.  (The old integrator added F/m to the velocity regardless of dt.)
#[test]
fn velocity_change_is_proportional_to_the_time_step() {
    let reference = film(1, 1, 1, 1.0e-3, 2.0e3, 0.0);
    let dt = 1.0e-7;
    let mut once = reference.clone();
    once.step(dt, [0.0; 3]);
    let mut twice = reference.clone();
    twice.step(0.5 * dt, [0.0; 3]);
    twice.step(0.5 * dt, [0.0; 3]);
    let mut doubled = reference.clone();
    doubled.step(2.0 * dt, [0.0; 3]);
    let forces = reference.net_forces();
    let mut checked = 0;
    for (index, node) in reference.nodes.iter().enumerate() {
        if node.pinned {
            continue;
        }
        let expected = scale3(forces[index], dt / node.mass);
        for (d, expected_d) in expected.iter().enumerate() {
            let v1 = once.nodes[index].velocity[d];
            assert!(
                (v1 - expected_d).abs() < 1e-9 * expected_d.abs().max(1e-12),
                "node {index} axis {d}: v = {v1:e}, F/m dt = {expected_d:e}"
            );
            let v2 = twice.nodes[index].velocity[d];
            assert!(
                (v1 - v2).abs() < 1e-3 * v1.abs().max(1e-12),
                "{v1:e} vs {v2:e}"
            );
            let v_double = doubled.nodes[index].velocity[d];
            assert!(
                (v_double - 2.0 * v1).abs() < 1e-9 * v1.abs().max(1e-12),
                "Δv must double with dt: {v1:e} -> {v_double:e}"
            );
        }
        checked += 1;
    }
    assert!(checked > 0);
}

fn scale3(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Relax a film under `voltage` and return `(thickness change of the top layer,
/// lateral growth of the top layer, largest remaining node force)`.
fn relax_film(voltage: f64) -> (f64, f64, f64) {
    let (nx, ny, nz, h) = (2usize, 2usize, 1usize, 1.0e-3);
    let mut body = film(nx, ny, nz, h, voltage, 2.0e4);
    let top = top_layer(nx, ny, nz);
    for _ in 0..4000 {
        body.step(1.0e-6, [0.0; 3]);
    }
    let thickness: f64 =
        top.iter().map(|&i| body.nodes[i].position[2]).sum::<f64>() / top.len() as f64;
    let width = body.nodes[top[top.len() - 1]].position[0] - body.nodes[top[0]].position[0];
    let forces = body.net_forces();
    let residual = body
        .nodes
        .iter()
        .zip(forces.iter())
        .filter(|(node, _)| !node.pinned)
        .map(|(_, f)| norm(*f))
        .fold(0.0_f64, f64::max);
    (thickness - h, width - nx as f64 * h, residual)
}

// A potential difference across the film squeezes it along the field and spreads
// it sideways; the squeeze grows with the square of the voltage, and the body
// comes to rest at an elastic equilibrium.
#[test]
fn field_squeezes_the_film_and_spreads_it_sideways() {
    let (squeeze_1, spread_1, residual_1) = relax_film(2.0e3);
    let (squeeze_2, spread_2, residual_2) = relax_film(4.0e3);
    assert!(
        squeeze_1 < 0.0 && squeeze_2 < 0.0,
        "{squeeze_1:e} {squeeze_2:e}"
    );
    assert!(
        spread_1 > 0.0 && spread_2 > 0.0,
        "{spread_1:e} {spread_2:e}"
    );
    let ratio = squeeze_2 / squeeze_1;
    assert!(
        (ratio - 4.0).abs() < 0.2,
        "squeeze ratio {ratio} for doubled voltage"
    );
    // Equilibrium: the elastic restoring force balances the Maxwell force.
    let maxwell_scale = EPSILON_0 * EPSILON_R_VHB * (4.0e3_f64 / 1.0e-3).powi(2) * 1.0e-6;
    assert!(
        residual_1 < 1e-3 * maxwell_scale && residual_2 < 1e-3 * maxwell_scale,
        "residual forces {residual_1:e}, {residual_2:e} (Maxwell scale {maxwell_scale:e})"
    );
    // The compression is of the order Maxwell pressure / Young's modulus (3 mu).
    let pressure = EPSILON_0 * EPSILON_R_VHB * (2.0e3_f64 / 1.0e-3).powi(2);
    let strain = -squeeze_1 / 1.0e-3;
    assert!(
        strain > 0.2 * pressure / (3.0 * MU) && strain < 5.0 * pressure / (3.0 * MU),
        "strain {strain:e} vs p / 3 mu = {:e}",
        pressure / (3.0 * MU)
    );
    // Spreading is the Poisson response: far smaller than the squeeze times the
    // aspect ratio and of the same sign convention as an incompressible solid.
    assert!(spread_1 < 2.0 * squeeze_1.abs() * 2.0);
}

// Gravity loads the elastic body: a column hanging from its pinned top stretches
// until the stress balances its weight.
#[test]
fn elastic_forces_hold_a_hanging_column() {
    let (nx, ny, nz, h) = (1usize, 1usize, 4usize, 1.0e-3);
    // A compressible material (nu = 0.35): linear tetrahedra lock volumetrically for
    // nearly incompressible ones and would converge to the beam value only slowly.
    let bulk = 3.0e5;
    let mut body = film_pinned_at(nx, ny, nz, h, 0.0, 2.0e4, nz, bulk);
    let g = [0.0, 0.0, -9.81];
    for _ in 0..60_000 {
        body.step(1.0e-6, g);
    }
    // Self-weight stretch of a column with the Young's modulus E = 9 K mu / (3 K + mu)
    // of the material: the free end sags by rho g L^2 / (2 E).
    let young = 9.0 * bulk * MU / (3.0 * bulk + MU);
    let length = nz as f64 * h;
    let expected = 1000.0 * 9.81 * length * length / (2.0 * young);
    let bottom: f64 = (0..(nx + 1) * (ny + 1))
        .map(|i| body.nodes[i].position[2])
        .sum::<f64>()
        / ((nx + 1) * (ny + 1)) as f64;
    let sag = -bottom;
    assert!(sag > 0.0, "the free end hangs below its rest position");
    assert!(
        (sag - expected).abs() < 0.12 * expected,
        "sag {sag:e} vs rho g L^2 / 2E = {expected:e}"
    );
    // The column is at rest: the elastic forces carry exactly the weight.
    let forces = body.net_forces();
    let weight: f64 = body
        .nodes
        .iter()
        .filter(|node| !node.pinned)
        .map(|node| node.mass * 9.81)
        .sum();
    let pull: f64 = body
        .nodes
        .iter()
        .zip(forces.iter())
        .filter(|(node, _)| !node.pinned)
        .map(|(_, f)| f[2])
        .sum();
    assert!(
        (pull - weight).abs() < 1e-3 * weight,
        "elastic pull {pull:e} vs weight {weight:e}"
    );
}

// Without damping the body conserves the free energy kinetic + elastic - electric
// (the electric energy enters with a minus sign at fixed potentials), up to the
// first-order accuracy of the symplectic Euler scheme: it fluctuates by a small
// fraction of the energy that is actually exchanged between the three forms.
#[test]
fn undamped_body_conserves_the_free_energy() {
    let mut body = film(2, 2, 1, 1.0e-3, 3.0e3, 0.0);
    let free_energy =
        |b: &ElectroacitiveFemBody| b.kinetic_energy() + b.elastic_energy() - b.electric_energy();
    let initial = free_energy(&body);
    let mut max_deviation = 0.0_f64;
    let mut max_kinetic = 0.0_f64;
    let mut max_motion = 0.0_f64;
    for _ in 0..4000 {
        body.step(2.0e-7, [0.0; 3]);
        max_deviation = max_deviation.max((free_energy(&body) - initial).abs());
        max_kinetic = max_kinetic.max(body.kinetic_energy());
        for node in body.nodes.iter().filter(|node| !node.pinned) {
            max_motion = max_motion.max((node.position[2] - 1.0e-3).abs());
        }
    }
    assert!(max_motion > 1.0e-6, "the film must move: {max_motion:e}");
    assert!(max_kinetic > 1.0e-11, "kinetic energy {max_kinetic:e}");
    // The exchanged energy is orders of magnitude below the stored electric energy,
    // so a small relative error of W_e would already swamp it: this is a sharp test.
    assert!(
        max_deviation < 0.02 * max_kinetic,
        "free energy drifts by {max_deviation:e}, exchanged energy {max_kinetic:e}"
    );
}
