// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the tissue dynamics: gravity is an acceleration,
//! damping is mass-proportional, and the tetrahedral elements of
//! [`TissueModel`] (linear elastic) and [`SoftTissueModel`] (neo-Hookean)
//! really exert forces.

use super::elasticity::{
    MAX_STABILITY_SUBSTEPS, STABILITY_SAFETY, TetRest, explicit_stability_limit,
    linear_elastic_stress, stability_substeps,
};
use super::functions::{add3, dot3, norm3, scale3};
use super::{
    CuttingModel, NeoHookeanParams, SimulationStep, SoftTissueModel, SurgicalTool, TissueEdge,
    TissueLayer, TissueModel, TissueNode, ToolGeometry,
};

const G: f64 = 9.81;
const GRAVITY: [f64; 3] = [0.0, -G, 0.0];
const NO_GRAVITY: [f64; 3] = [0.0; 3];

fn skin() -> TissueLayer {
    TissueLayer::new("skin", 200_000.0, 0.45, 0.002, 1100.0, 500.0)
}

fn fat() -> TissueLayer {
    TissueLayer::new("fat", 3_000.0, 0.49, 0.005, 920.0, 50.0)
}

/// Unit tetrahedron with a pinned node 0 and three free nodes of `mass`.
fn unit_tet_tissue(mass: f64) -> TissueModel {
    let mut model = TissueModel::new(vec![skin()], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    model.add_node(TissueNode::new_static([0.0, 0.0, 0.0], 0));
    model.add_node(TissueNode::new([1.0, 0.0, 0.0], mass, 0));
    model.add_node(TissueNode::new([0.0, 1.0, 0.0], mass, 0));
    model.add_node(TissueNode::new([0.0, 0.0, 1.0], mass, 0));
    model.add_element([0, 1, 2, 3]);
    model
}

/// Unit tetrahedron whose base (nodes 0, 1, 2) is pinned; only the apex
/// (node 3) is free, so the structure has no rigid-body mode.
fn base_pinned_tet_tissue(mass: f64) -> TissueModel {
    let mut model = TissueModel::new(vec![skin()], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    model.add_node(TissueNode::new_static([0.0, 0.0, 0.0], 0));
    model.add_node(TissueNode::new_static([1.0, 0.0, 0.0], 0));
    model.add_node(TissueNode::new_static([0.0, 1.0, 0.0], 0));
    model.add_node(TissueNode::new([0.0, 0.0, 1.0], mass, 0));
    model.add_element([0, 1, 2, 3]);
    model
}

/// Two tetrahedra sharing the face (1, 2, 3); the second one lies in the fat
/// layer.  Node 0 belongs to the first tetrahedron only, node 4 to the second.
fn two_layer_tissue() -> TissueModel {
    let mut model = TissueModel::new(vec![skin(), fat()], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    model.add_node(TissueNode::new([0.0, 0.0, 0.0], 0.1, 0));
    model.add_node(TissueNode::new([1.0, 0.0, 0.0], 0.1, 1));
    model.add_node(TissueNode::new([0.0, 1.0, 0.0], 0.1, 0));
    model.add_node(TissueNode::new([0.0, 0.0, 1.0], 0.1, 0));
    model.add_node(TissueNode::new([1.0, 1.0, 1.0], 0.1, 0));
    model.add_element([0, 1, 2, 3]);
    model.add_element([1, 2, 3, 4]);
    model
}

fn total_force(forces: &[[f64; 3]]) -> [f64; 3] {
    forces.iter().fold([0.0; 3], |sum, f| add3(sum, *f))
}

// ---------------------------------------------------------------------------
// TissueModel: gravity, damping and applied forces
// ---------------------------------------------------------------------------

#[test]
fn tissue_free_fall_acceleration_is_independent_of_mass() {
    let mut model = TissueModel::new(Vec::new(), 0.0);
    model.damping_alpha = 0.0;
    for (i, mass) in [1.0, 2.0, 0.01].into_iter().enumerate() {
        model.add_node(TissueNode::new([i as f64, 10.0, 0.0], mass, 0));
    }
    let dt = 0.01;
    model.step(dt, GRAVITY);
    for node in &model.nodes {
        assert!(
            (node.velocity[1] + G * dt).abs() < 1.0e-12,
            "a = g for every mass, got v = {} after one step",
            node.velocity[1]
        );
    }
    for _ in 1..100 {
        model.step(dt, GRAVITY);
    }
    // Symplectic Euler: v_N = -g N dt, y_N = y_0 - g dt^2 N (N + 1) / 2.
    let expected_y = 10.0 - G * dt * dt * 100.0 * 101.0 / 2.0;
    for node in &model.nodes {
        assert!(
            (node.velocity[1] + G).abs() < 1.0e-9,
            "v = {}",
            node.velocity[1]
        );
        assert!(
            (node.position[1] - expected_y).abs() < 1.0e-9,
            "y = {} vs {expected_y}",
            node.position[1]
        );
    }
}

#[test]
fn tissue_gravity_is_not_divided_by_mass() {
    // A light node used to be accelerated with g / m = 981 m/s^2.
    let mut model = TissueModel::new(Vec::new(), 0.0);
    model.damping_alpha = 0.0;
    model.add_node(TissueNode::new([0.0; 3], 0.01, 0));
    model.step(0.1, GRAVITY);
    assert!(
        (model.nodes[0].velocity[1] + 0.981).abs() < 1.0e-12,
        "v = {}",
        model.nodes[0].velocity[1]
    );
}

#[test]
fn tissue_damping_is_mass_proportional() {
    // Force -alpha m v means a = -alpha v, whatever the mass.
    let mut model = TissueModel::new(Vec::new(), 0.0);
    model.damping_alpha = 2.0;
    for mass in [1.0, 5.0] {
        let mut node = TissueNode::new([0.0; 3], mass, 0);
        node.velocity = [1.0, 0.0, 0.0];
        model.add_node(node);
    }
    model.step(0.01, NO_GRAVITY);
    for node in &model.nodes {
        assert!(
            (node.velocity[0] - 0.98).abs() < 1.0e-12,
            "v = {} should be 1 - alpha dt",
            node.velocity[0]
        );
    }
}

#[test]
fn tissue_applied_force_is_divided_by_mass_and_cleared() {
    let mut model = TissueModel::new(Vec::new(), 0.0);
    model.damping_alpha = 0.0;
    let mut node = TissueNode::new([0.0; 3], 2.0, 0);
    node.force = [10.0, 0.0, 0.0];
    model.add_node(node);
    let mut pinned = TissueNode::new_static([1.0, 0.0, 0.0], 0);
    pinned.force = [5.0, 0.0, 0.0];
    model.add_node(pinned);
    model.step(0.1, NO_GRAVITY);
    assert!(
        (model.nodes[0].velocity[0] - 0.5).abs() < 1.0e-12,
        "a = F / m = 5"
    );
    assert_eq!(model.nodes[0].force, [0.0; 3]);
    assert_eq!(
        model.nodes[1].force, [0.0; 3],
        "forces on pinned nodes are cleared too"
    );
    assert_eq!(model.nodes[1].position, [1.0, 0.0, 0.0]);
}

// ---------------------------------------------------------------------------
// TissueModel: elastic element forces
// ---------------------------------------------------------------------------

#[test]
fn tissue_rest_configuration_is_force_free() {
    let model = two_layer_tissue();
    for force in model.internal_forces() {
        assert_eq!(force, [0.0; 3]);
    }
    assert_eq!(model.elastic_energy(), 0.0);
}

#[test]
fn tissue_stretched_node_feels_the_analytic_restoring_force() {
    let mut model = unit_tet_tissue(0.01);
    model.nodes[1].position[0] = 1.1;
    let forces = model.internal_forces();

    // Dm = I, V0 = 1/6, displacement gradient = 0.1 e_x e_x^T:
    // f_1 = -V0 (lambda + 2 mu) 0.1 e_x and f_2 = -V0 lambda 0.1 e_y.
    let skin = skin();
    let (lambda, mu) = (skin.lame_lambda(), skin.shear_modulus());
    let v0 = 1.0 / 6.0;
    let expected_f1 = -v0 * (lambda + 2.0 * mu) * 0.1;
    let expected_f2 = -v0 * lambda * 0.1;
    assert!(
        (forces[1][0] - expected_f1).abs() < 1.0e-9 * expected_f1.abs(),
        "f1.x = {} vs {expected_f1}",
        forces[1][0]
    );
    assert!(forces[1][1].abs() < 1.0e-9 && forces[1][2].abs() < 1.0e-9);
    assert!((forces[2][1] - expected_f2).abs() < 1.0e-9 * expected_f2.abs());
    assert!((forces[3][2] - expected_f2).abs() < 1.0e-9 * expected_f2.abs());
    let net = norm3(total_force(&forces));
    assert!(
        net < 1.0e-9 * expected_f1.abs(),
        "internal forces must sum to zero, got {net}"
    );
}

#[test]
fn tissue_forces_are_the_negative_energy_gradient() {
    let mut model = two_layer_tissue();
    model.pre_stress = 5_000.0;
    // A non-trivial deformation of all five nodes.
    let offsets = [
        [0.02, -0.01, 0.03],
        [0.05, 0.02, -0.04],
        [-0.03, 0.04, 0.01],
        [0.01, -0.05, 0.06],
        [-0.04, 0.03, -0.02],
    ];
    for (node, offset) in model.nodes.iter_mut().zip(offsets.iter()) {
        node.position = add3(node.position, *offset);
    }
    let forces = model.internal_forces();
    let h = 1.0e-6;
    for (node, node_forces) in forces.iter().enumerate() {
        for (axis, &force) in node_forces.iter().enumerate() {
            let original = model.nodes[node].position[axis];
            model.nodes[node].position[axis] = original + h;
            let plus = model.elastic_energy();
            model.nodes[node].position[axis] = original - h;
            let minus = model.elastic_energy();
            model.nodes[node].position[axis] = original;
            let gradient = (plus - minus) / (2.0 * h);
            assert!(
                (gradient + force).abs() <= 1.0e-6 * (1.0 + force.abs()),
                "node {node} axis {axis}: -dE/dx = {} vs force {force}",
                -gradient
            );
        }
    }
}

#[test]
fn tissue_elements_use_the_material_of_their_layer() {
    // Two layers with the same Poisson ratio; the second is four times stiffer.
    let layers = || {
        vec![
            TissueLayer::new("soft", 1.0e4, 0.3, 0.01, 1000.0, 100.0),
            TissueLayer::new("stiff", 4.0e4, 0.3, 0.01, 1000.0, 100.0),
        ]
    };
    let build = |second_layer: usize| {
        let mut model = TissueModel::new(layers(), 0.0);
        model.add_node(TissueNode::new([0.0, 0.0, 0.0], 0.1, 0));
        model.add_node(TissueNode::new([1.0, 0.0, 0.0], 0.1, second_layer));
        model.add_node(TissueNode::new([0.0, 1.0, 0.0], 0.1, 0));
        model.add_node(TissueNode::new([0.0, 0.0, 1.0], 0.1, 0));
        model.add_node(TissueNode::new([1.0, 1.0, 1.0], 0.1, 0));
        model.add_element([0, 1, 2, 3]);
        model.add_element([1, 2, 3, 4]);
        model.nodes[0].position = [-0.05, 0.02, 0.0];
        model.nodes[4].position = [1.0, 1.1, 1.05];
        model
    };
    let uniform = build(0).internal_forces();
    let layered = build(1).internal_forces();
    // Node 0 only belongs to the first element, node 4 only to the second.
    assert_eq!(uniform[0], layered[0], "first element keeps its material");
    for axis in 0..3 {
        assert!(
            (layered[4][axis] - 4.0 * uniform[4][axis]).abs()
                <= 1.0e-9 * uniform[4][axis].abs().max(1.0),
            "second element is four times stiffer on axis {axis}: {} vs {}",
            layered[4][axis],
            uniform[4][axis]
        );
    }
    assert!(norm3(uniform[4]) > 0.0);
}

#[test]
fn tissue_pre_stress_acts_on_the_element_boundary() {
    let mut model = unit_tet_tissue(1.0);
    model.pre_stress = 6_000.0;
    let forces = model.internal_forces();
    // Dm = I: f_k = -V0 sigma_0 e_k for k = 1..3 and f_0 = -sum f_k.
    let expected = -6_000.0 / 6.0;
    assert!((forces[1][0] - expected).abs() < 1.0e-9, "{:?}", forces[1]);
    assert!((forces[2][1] - expected).abs() < 1.0e-9, "{:?}", forces[2]);
    assert!((forces[3][2] - expected).abs() < 1.0e-9, "{:?}", forces[3]);
    for axis in 0..3 {
        assert!(
            (forces[0][axis] + expected).abs() < 1.0e-9,
            "{:?}",
            forces[0]
        );
    }
    assert!(norm3(total_force(&forces)) < 1.0e-9);
}

/// Total (kinetic + elastic) energy of an undamped stretched unit tetrahedron
/// released at rest and integrated for `steps` steps of `dt`: returns the
/// relative peak-to-peak energy fluctuation and the closest approach of the
/// released node to its rest position.
fn stretched_tet_energy_fluctuation(dt: f64, steps: usize) -> (f64, f64) {
    let mut model = unit_tet_tissue(0.01);
    model.nodes[1].position[0] = 1.1;
    let initial_energy = model.elastic_energy();
    assert!(initial_energy > 0.0);
    assert!(
        model.internal_forces()[1][0] < 0.0,
        "the force must pull node 1 back"
    );

    let mut closest_to_rest = f64::INFINITY;
    let mut max_total_energy = initial_energy;
    let mut min_total_energy = initial_energy;
    for _ in 0..steps {
        model.step(dt, NO_GRAVITY);
        closest_to_rest = closest_to_rest.min((model.nodes[1].position[0] - 1.0).abs());
        let kinetic: f64 = model
            .nodes
            .iter()
            .map(|n| 0.5 * n.mass * norm3(n.velocity).powi(2))
            .sum();
        let total = kinetic + model.elastic_energy();
        max_total_energy = max_total_energy.max(total);
        min_total_energy = min_total_energy.min(total);
    }
    (
        (max_total_energy - min_total_energy) / initial_energy,
        closest_to_rest,
    )
}

#[test]
fn tissue_stretched_tet_returns_toward_rest_and_conserves_energy_without_damping() {
    // Symplectic Euler does not conserve the energy exactly: it oscillates
    // around the true value with an amplitude proportional to the step. Halving
    // the step must therefore halve the fluctuation, and for a step far below
    // the stability limit the fluctuation is small.
    let (coarse, closest_coarse) = stretched_tet_energy_fluctuation(5.0e-6, 2000);
    let (fine, closest_fine) = stretched_tet_energy_fluctuation(2.5e-6, 4000);
    assert!(
        closest_coarse < 0.02 && closest_fine < 0.02,
        "the released node must swing back through its rest position, closest approach \
         {closest_coarse} / {closest_fine}"
    );
    assert!(
        coarse < 0.05,
        "undamped energy must stay bounded, fluctuation {coarse}"
    );
    assert!(
        fine < 0.65 * coarse,
        "energy error must shrink with the step: {fine} (dt/2) vs {coarse} (dt)"
    );
}

#[test]
fn tissue_damped_tet_settles_and_energy_decays() {
    let mut model = unit_tet_tissue(0.01);
    model.damping_alpha = 3_000.0;
    model.nodes[1].position[0] = 1.1;
    let initial_energy = model.elastic_energy();
    let initial_offset = 0.1;
    for _ in 0..40 {
        model.step(1.0e-3, NO_GRAVITY);
    }
    let offset = (model.nodes[1].position[0] - 1.0).abs();
    assert!(
        offset < 1.0e-3 * initial_offset,
        "node must settle at its rest position, offset {offset}"
    );
    assert!(model.elastic_energy() < 1.0e-6 * initial_energy);
    let kinetic: f64 = model
        .nodes
        .iter()
        .map(|n| 0.5 * n.mass * norm3(n.velocity).powi(2))
        .sum();
    assert!(
        kinetic < 1.0e-6 * initial_energy,
        "kinetic energy {kinetic}"
    );
}

#[test]
fn tissue_stiffness_proportional_damping_dissipates() {
    let run = |beta: f64| {
        let mut model = unit_tet_tissue(0.01);
        model.damping_beta = beta;
        model.nodes[1].position[0] = 1.1;
        let initial = model.elastic_energy();
        for _ in 0..100 {
            model.step(1.0e-4, NO_GRAVITY);
        }
        let kinetic: f64 = model
            .nodes
            .iter()
            .map(|n| 0.5 * n.mass * norm3(n.velocity).powi(2))
            .sum();
        (kinetic + model.elastic_energy()) / initial
    };
    let undamped = run(0.0);
    let damped = run(2.0e-4);
    assert!(
        damped < 0.5 * undamped,
        "beta damping must remove energy: {damped} vs undamped {undamped}"
    );
}

#[test]
fn tissue_hanging_under_gravity_reaches_the_analytic_sag() {
    // Pure shear of the apex: sigma_yz = mu * delta, f_y = -V0 mu delta, so the
    // weight m g is balanced at delta = -m g / (V0 mu).
    let mass = 0.01;
    let mut model = base_pinned_tet_tissue(mass);
    model.damping_alpha = 3_000.0;
    for _ in 0..100 {
        model.step(1.0e-3, GRAVITY);
    }
    let expected_sag = -mass * G / (skin().shear_modulus() / 6.0);
    let sag = model.nodes[3].position[1];
    assert!(
        (sag - expected_sag).abs() < 1.0e-3 * expected_sag.abs(),
        "apex sag {sag} vs analytic {expected_sag}"
    );
    assert!(norm3(model.nodes[3].velocity) < 1.0e-9);
    let residual = norm3(add3(model.internal_forces()[3], [0.0, -mass * G, 0.0]));
    assert!(residual < 1.0e-3 * mass * G, "force residual {residual}");
}

#[test]
fn tissue_explicit_step_is_stable_far_above_the_naive_limit() {
    // A stiff element integrated at 60 Hz with the default Rayleigh damping:
    // without sub-stepping this blows up.  The base is pinned so that the
    // structure has no zero-energy rigid-body mode that gravity could drive.
    let mut model = TissueModel::new(vec![skin()], 0.0);
    model.add_node(TissueNode::new_static([0.0, 0.0, 0.0], 0));
    model.add_node(TissueNode::new_static([1.0, 0.0, 0.0], 0));
    model.add_node(TissueNode::new_static([0.0, 1.0, 0.0], 0));
    model.add_node(TissueNode::new([0.0, 0.0, 1.0], 0.01, 0));
    model.add_element([0, 1, 2, 3]);
    let dt = 1.0 / 60.0;
    let critical = model.critical_time_step();
    assert!(
        critical.is_finite() && critical < dt / 10.0,
        "critical dt = {critical}"
    );
    for _ in 0..120 {
        model.step(dt, GRAVITY);
    }
    let apex = &model.nodes[3];
    assert!(
        apex.position.iter().all(|c| c.is_finite() && c.abs() < 2.0),
        "unstable: {:?}",
        apex.position
    );
    // After two seconds the apex rests where its weight is carried by the
    // shear of the element: delta = -m g / (V0 mu), V0 = 1/6.
    let expected_sag = -0.01 * G / (skin().shear_modulus() / 6.0);
    assert!(
        (apex.position[1] - expected_sag).abs() < 1.0e-2 * expected_sag.abs(),
        "apex sag {} vs analytic {expected_sag}",
        apex.position[1]
    );
    assert!(
        norm3(apex.velocity) < 1.0e-6,
        "apex must be at rest, v = {:?}",
        apex.velocity
    );
}

#[test]
fn tissue_critical_time_step_scales_with_stiffness() {
    let make = |young: f64| {
        let mut model = TissueModel::new(
            vec![TissueLayer::new("tissue", young, 0.3, 0.01, 1000.0, 100.0)],
            0.0,
        );
        model.damping_alpha = 0.0;
        model.damping_beta = 0.0;
        model.add_node(TissueNode::new_static([0.0; 3], 0));
        model.add_node(TissueNode::new([1.0, 0.0, 0.0], 0.1, 0));
        model.add_node(TissueNode::new([0.0, 1.0, 0.0], 0.1, 0));
        model.add_node(TissueNode::new([0.0, 0.0, 1.0], 0.1, 0));
        model.add_element([0, 1, 2, 3]);
        model
    };
    let soft = make(1.0e4).critical_time_step();
    let stiff = make(4.0e4).critical_time_step();
    assert!(
        (soft / stiff - 2.0).abs() < 1.0e-9,
        "dt_crit ~ 1/sqrt(E): {soft} vs {stiff}"
    );

    // Without elements the only limit is the mass-proportional damping, 2 / alpha.
    let mut free = TissueModel::new(Vec::new(), 0.0);
    free.damping_alpha = 0.0;
    free.damping_beta = 0.0;
    free.add_node(TissueNode::new([0.0; 3], 1.0, 0));
    assert!(free.critical_time_step().is_infinite());
    free.damping_alpha = 4.0;
    assert!((free.critical_time_step() - 0.5).abs() < 1.0e-12);
}

// ---------------------------------------------------------------------------
// TissueModel: severed elements
// ---------------------------------------------------------------------------

#[test]
fn tissue_severed_element_stops_carrying_load() {
    let mut model = two_layer_tissue();
    for node in model.nodes.iter_mut() {
        node.position = add3(node.position, [0.03, -0.02, 0.04]);
    }
    model.nodes[0].position = add3(model.nodes[0].position, [-0.1, 0.0, 0.0]);
    model.nodes[4].position = add3(model.nodes[4].position, [0.0, 0.1, 0.0]);
    let intact = model.internal_forces();
    assert!(norm3(intact[0]) > 0.0 && norm3(intact[4]) > 0.0);

    assert!(!model.is_element_severed(0));
    assert!(model.sever_element(0));
    assert!(model.is_element_severed(0));
    assert!(!model.is_element_severed(1));
    assert!(
        !model.sever_element(7),
        "out-of-range elements cannot be severed"
    );
    assert!(!model.is_element_severed(7));

    let cut = model.internal_forces();
    assert_eq!(
        cut[0], [0.0; 3],
        "node 0 only belongs to the severed element"
    );
    assert!(
        norm3(cut[4]) > 0.0,
        "the other element still pulls on node 4"
    );
    // The stored energy is that of the second element alone: rebuild it as a
    // stand-alone model (rest positions first, then the same deformation).
    let mut second = TissueModel::new(model.layers.clone(), 0.0);
    for node in 1..5 {
        second.add_node(TissueNode::new(
            model.rest_positions[node],
            0.1,
            model.nodes[node].layer_index,
        ));
    }
    second.add_element([0, 1, 2, 3]);
    for node in 1..5 {
        second.nodes[node - 1].position = model.nodes[node].position;
    }
    let expected = second.elastic_energy();
    assert!(expected > 0.0);
    assert!(
        (model.elastic_energy() - expected).abs() < 1.0e-9 * expected,
        "{} vs {expected}",
        model.elastic_energy()
    );
}

#[test]
fn tissue_edges_sever_exactly_the_elements_that_contain_them() {
    let mut model = two_layer_tissue();
    // Edge (0, 1) lies only in element 0; edge (2, 3) is shared by both.
    assert_eq!(model.sever_elements_by_edges([(1, 0)]), 1);
    assert!(model.is_element_severed(0) && !model.is_element_severed(1));
    assert_eq!(
        model.sever_elements_by_edges([(0, 1)]),
        0,
        "already severed"
    );
    assert_eq!(
        model.sever_elements_by_edges([(0, 4)]),
        0,
        "no element contains both nodes"
    );
    assert_eq!(
        model.sever_elements_by_edges([(3, 2)]),
        1,
        "only element 1 is left"
    );
    assert!(model.is_element_severed(1));
    assert_eq!(model.sever_elements_by_edges(std::iter::empty()), 0);
}

#[test]
fn simulation_step_cut_unloads_the_cut_element() {
    let mut tissue = TissueModel::new(vec![skin()], 0.0);
    tissue.damping_alpha = 0.0;
    tissue.damping_beta = 0.0;
    tissue.add_node(TissueNode::new_static([0.0, 0.0, 0.0], 0));
    tissue.add_node(TissueNode::new([0.01, 0.0, 0.0], 0.01, 0));
    tissue.add_node(TissueNode::new([0.0, 0.01, 0.0], 0.01, 0));
    tissue.add_node(TissueNode::new([0.0, 0.0, 0.01], 0.01, 0));
    tissue.add_element([0, 1, 2, 3]);
    // The probe presses on node 1 with enough force to cut the edge (1, 2).
    let tool = SurgicalTool::new("blade", ToolGeometry::Blade, 0.01, [0.0105, 0.0, 0.0]);
    let mut cutting = CuttingModel::new(1.0e-3, 1.0);
    cutting.add_edge(TissueEdge::new(1, 2, 0.01 * 2.0_f64.sqrt()));
    let mut sim = SimulationStep::new(tissue, tool, cutting, 1.0e5, 0.001);

    assert!(!sim.tissue.is_element_severed(0));
    sim.advance(1.0e-3, NO_GRAVITY);
    assert_eq!(
        sim.num_cuts(),
        1,
        "the edge must be severed by the contact force"
    );
    assert!(
        sim.tissue.is_element_severed(0),
        "the element that contains the cut edge no longer carries load"
    );
    // A cut element exerts nothing, even when strongly deformed.
    sim.tissue.nodes[2].position = [0.0, 0.02, 0.0];
    for force in sim.tissue.internal_forces() {
        assert_eq!(force, [0.0; 3]);
    }
}

#[test]
fn simulation_step_applies_gravity_as_acceleration() {
    let mut tissue = TissueModel::new(Vec::new(), 0.0);
    tissue.damping_alpha = 0.0;
    tissue.add_node(TissueNode::new([0.0, 0.0, 0.0], 5.0, 0));
    tissue.add_node(TissueNode::new([1.0, 0.0, 0.0], 0.001, 0));
    let tool = SurgicalTool::new("probe", ToolGeometry::Sphere, 0.001, [10.0; 3]);
    let mut sim = SimulationStep::new(tissue, tool, CuttingModel::new(1.0, 0.5), 1.0e5, 0.001);
    sim.advance(0.01, GRAVITY);
    for node in &sim.tissue.nodes {
        assert!(
            (node.velocity[1] + G * 0.01).abs() < 1.0e-12,
            "v = {}",
            node.velocity[1]
        );
    }
}

// ---------------------------------------------------------------------------
// SoftTissueModel: neo-Hookean forces
// ---------------------------------------------------------------------------

fn neo_tet(mass: f64, alpha: f64, eta: f64) -> SoftTissueModel {
    let mut model = SoftTissueModel::new(alpha);
    let material = NeoHookeanParams::from_young(1.0e4, 0.3, eta);
    model.add_node([0.0, 0.0, 0.0], 0.0, true);
    model.add_node([1.0, 0.0, 0.0], mass, false);
    model.add_node([0.0, 1.0, 0.0], mass, false);
    model.add_node([0.0, 0.0, 1.0], mass, false);
    model.add_element([0, 1, 2, 3], material);
    model
}

fn kinetic_energy(model: &SoftTissueModel) -> f64 {
    model
        .masses
        .iter()
        .zip(model.velocities.iter())
        .zip(model.fixed.iter())
        .filter(|(_, fixed)| !**fixed)
        .map(|((m, v), _)| 0.5 * m * norm3(*v).powi(2))
        .sum()
}

#[test]
fn soft_tissue_rest_configuration_is_force_free() {
    let model = neo_tet(0.01, 0.0, 1.0);
    for force in model.internal_forces() {
        assert!(norm3(force) < 1.0e-9, "{force:?}");
    }
    assert!(model.strain_energy().abs() < 1.0e-12);
    assert!((model.rest_volume(0) - 1.0 / 6.0).abs() < 1.0e-12);
}

#[test]
fn soft_tissue_forces_are_the_negative_energy_gradient() {
    let mut model = neo_tet(0.01, 0.0, 0.0);
    model.nodes[1] = [1.2, 0.1, 0.0];
    model.nodes[2] = [0.05, 0.9, 0.1];
    model.nodes[3] = [0.0, 0.1, 1.1];
    let forces = model.internal_forces();
    assert!(
        norm3(total_force(&forces)) < 1.0e-8,
        "net force must vanish"
    );
    let h = 1.0e-6;
    for (node, node_forces) in forces.iter().enumerate() {
        for (axis, &force) in node_forces.iter().enumerate() {
            let original = model.nodes[node][axis];
            model.nodes[node][axis] = original + h;
            let plus = model.strain_energy();
            model.nodes[node][axis] = original - h;
            let minus = model.strain_energy();
            model.nodes[node][axis] = original;
            let gradient = (plus - minus) / (2.0 * h);
            assert!(
                (gradient + force).abs() <= 1.0e-5 * (1.0 + force.abs()),
                "node {node} axis {axis}: -dE/dx = {} vs force {force}",
                -gradient
            );
        }
    }
}

#[test]
fn soft_tissue_elements_pull_stretched_nodes_back() {
    let mut model = neo_tet(0.01, 0.0, 0.0);
    model.nodes[1][0] = 1.1;
    let forces = model.internal_forces();
    assert!(
        forces[1][0] < 0.0,
        "stretched node must be pulled back, {:?}",
        forces[1]
    );
    assert!(model.strain_energy() > 0.0);
}

#[test]
fn soft_tissue_stretched_tet_returns_toward_rest_and_energy_decays() {
    let run = |alpha: f64| {
        let mut model = neo_tet(0.01, alpha, 0.0);
        model.nodes[1][0] = 1.1;
        let initial = model.strain_energy();
        let mut closest = f64::INFINITY;
        for _ in 0..400 {
            model.step(1.0e-4, NO_GRAVITY);
            closest = closest.min((model.nodes[1][0] - 1.0).abs());
        }
        (
            initial,
            closest,
            model.strain_energy() + kinetic_energy(&model),
        )
    };
    let (initial, closest_undamped, total_undamped) = run(0.0);
    assert!(
        closest_undamped < 0.02,
        "released node must swing back toward rest, closest {closest_undamped}"
    );
    assert!(
        (total_undamped - initial).abs() < 0.1 * initial,
        "undamped energy {total_undamped} vs initial {initial}"
    );
    let (_, _, total_damped) = run(2_000.0);
    assert!(
        total_damped < 0.1 * total_undamped,
        "mass damping must dissipate: {total_damped} vs {total_undamped}"
    );
}

#[test]
fn soft_tissue_viscosity_dissipates_energy() {
    let run = |eta: f64| {
        let mut model = neo_tet(0.01, 0.0, eta);
        model.nodes[1][0] = 1.1;
        for _ in 0..400 {
            model.step(1.0e-4, NO_GRAVITY);
        }
        model.strain_energy() + kinetic_energy(&model)
    };
    let elastic_only = run(0.0);
    let viscous = run(50.0);
    assert!(
        viscous < 0.5 * elastic_only,
        "viscosity must remove energy: {viscous} vs {elastic_only}"
    );
}

#[test]
fn soft_tissue_viscous_force_ignores_rigid_rotation() {
    // v = omega x X is a rigid rotation: velocity gradient is skew, no dissipation.
    let mut model = neo_tet(0.01, 0.0, 200.0);
    let omega = [0.0, 0.0, 3.0];
    let rest = model.rest_positions.clone();
    for (velocity, x) in model.velocities.iter_mut().zip(rest.iter()) {
        *velocity = [
            omega[1] * x[2] - omega[2] * x[1],
            omega[2] * x[0] - omega[0] * x[2],
            omega[0] * x[1] - omega[1] * x[0],
        ];
    }
    for force in model.internal_forces() {
        assert!(
            norm3(force) < 1.0e-9,
            "rigid rotation must not be damped: {force:?}"
        );
    }
    // Pure stretching rate does produce a resisting force.
    let mut stretching = neo_tet(0.01, 0.0, 200.0);
    stretching.velocities[1] = [0.5, 0.0, 0.0];
    assert!(stretching.internal_forces()[1][0] < 0.0);
}

#[test]
fn soft_tissue_inverted_element_exerts_no_force() {
    let mut model = neo_tet(0.01, 0.0, 1.0);
    model.nodes[3] = [0.0, 0.0, -1.0];
    for force in model.internal_forces() {
        assert_eq!(force, [0.0; 3]);
    }
    model.step(0.01, GRAVITY);
    assert!(model.nodes.iter().all(|n| n.iter().all(|c| c.is_finite())));
}

#[test]
fn soft_tissue_free_fall_is_mass_independent() {
    let mut model = SoftTissueModel::new(0.0);
    model.add_node([0.0, 5.0, 0.0], 1.0, false);
    model.add_node([1.0, 5.0, 0.0], 0.01, false);
    model.step(0.01, GRAVITY);
    for velocity in &model.velocities {
        assert!(
            (velocity[1] + 0.0981).abs() < 1.0e-12,
            "v = {}",
            velocity[1]
        );
    }
}

#[test]
fn soft_tissue_hanging_under_gravity_reaches_the_analytic_sag() {
    // Base pinned, apex free: for small shear the neo-Hookean shear modulus is mu,
    // so the apex settles at delta = -m g / (V0 mu).
    //
    // The apex oscillates with omega = sqrt(V0 mu / m) = 253 rad/s; mass
    // damping of 2 omega is critical damping, so it settles within a few
    // hundredths of a second.
    let mass = 0.01;
    let mut model = SoftTissueModel::new(500.0);
    let material = NeoHookeanParams::from_young(1.0e4, 0.3, 0.0);
    let shear_modulus = material.mu;
    model.add_node([0.0, 0.0, 0.0], 0.0, true);
    model.add_node([1.0, 0.0, 0.0], 0.0, true);
    model.add_node([0.0, 1.0, 0.0], 0.0, true);
    model.add_node([0.0, 0.0, 1.0], mass, false);
    model.add_element([0, 1, 2, 3], material);
    for _ in 0..100 {
        model.step(1.0e-3, GRAVITY);
    }
    let expected_sag = -mass * G / (shear_modulus / 6.0);
    let sag = model.nodes[3][1];
    assert!(
        (sag - expected_sag).abs() < 1.0e-4 * expected_sag.abs(),
        "apex sag {sag} vs analytic {expected_sag}"
    );
    assert!(norm3(model.velocities[3]) < 1.0e-6);
    let residual = norm3(add3(model.internal_forces()[3], [0.0, -mass * G, 0.0]));
    assert!(residual < 1.0e-3 * mass * G, "force residual {residual}");
}

#[test]
fn soft_tissue_explicit_step_is_stable_far_above_the_naive_limit() {
    let mut model = neo_tet(0.01, 0.02, 1.0);
    let dt = 0.1;
    let critical = model.critical_time_step();
    assert!(
        critical.is_finite() && critical < dt / 10.0,
        "critical dt = {critical}"
    );
    for _ in 0..50 {
        model.step(dt, GRAVITY);
    }
    for node in &model.nodes {
        assert!(
            node.iter().all(|c| c.is_finite() && c.abs() < 2.0),
            "unstable: {node:?}"
        );
    }
}

#[test]
fn soft_tissue_materials_set_the_stiffness_of_their_element() {
    let build = |young: f64| {
        let mut model = SoftTissueModel::new(0.0);
        let material = NeoHookeanParams::from_young(young, 0.3, 0.0);
        model.add_node([0.0, 0.0, 0.0], 0.0, true);
        model.add_node([1.1, 0.0, 0.0], 0.01, false);
        model.add_node([0.0, 1.0, 0.0], 0.01, false);
        model.add_node([0.0, 0.0, 1.0], 0.01, false);
        model.rest_positions[1] = [1.0, 0.0, 0.0];
        model.add_element([0, 1, 2, 3], material);
        model
    };
    let soft = norm3(build(1.0e3).internal_forces()[1]);
    let stiff = norm3(build(1.0e4).internal_forces()[1]);
    assert!(
        (stiff / soft - 10.0).abs() < 1.0e-6,
        "force scales with the modulus: {soft} vs {stiff}"
    );
}

#[test]
fn soft_tissue_energy_and_density_helpers_agree() {
    let mut model = neo_tet(0.01, 0.0, 0.0);
    model.nodes[1][0] = 1.2;
    let by_density = model.total_energy() * model.rest_volume(0);
    assert!((model.strain_energy() - by_density).abs() < 1.0e-12 * by_density.abs().max(1.0));
    // `element_volume` reports the deformed volume, `rest_volume` the reference one.
    assert!((model.element_volume(0) - 1.2 / 6.0).abs() < 1.0e-12);
    assert!((model.rest_volume(0) - 1.0 / 6.0).abs() < 1.0e-12);
    assert_eq!(model.rest_volume(5), 0.0, "unknown elements have no volume");
}

// ---------------------------------------------------------------------------
// Tetrahedron rest data, nodal forces and the explicit stability estimates
// ---------------------------------------------------------------------------

const UNIT_TET: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Row-major 3x3 matrix times vector.
fn mat3_apply(m: &[f64; 9], v: [f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

/// Tetrahedra of very different shape, size and orientation.
fn sample_tetrahedra() -> Vec<[[f64; 3]; 4]> {
    vec![
        UNIT_TET,
        // Regular tetrahedron with unit edges.
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866_025_403_784_438_6, 0.0],
            [0.5, 0.288_675_134_594_812_9, 0.816_496_580_927_726],
        ],
        // Stretched by a different factor along every axis.
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ],
        // Needle-like and sheared.
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.3, 0.05, 0.0],
            [0.2, 0.1, 0.02],
        ],
        // Centimetre-sized, far from the origin, negatively oriented.
        [
            [10.0, 10.0, 10.0],
            [10.01, 10.0, 10.0],
            [10.0, 10.0, 10.01],
            [10.0, 10.01, 10.0],
        ],
    ]
}

#[test]
fn tet_rest_volume_and_inverse_edge_matrix() {
    let unit = TetRest::from_positions(UNIT_TET).expect("the unit tetrahedron is valid");
    assert!((unit.volume - 1.0 / 6.0).abs() < 1.0e-15);
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    for (got, expected) in unit.dm_inv.iter().zip(identity.iter()) {
        assert!((got - expected).abs() < 1.0e-15, "{:?}", unit.dm_inv);
    }

    let stretched =
        TetRest::from_positions([[0.0; 3], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]])
            .expect("a box-shaped tetrahedron is valid");
    assert!((stretched.volume - 4.0).abs() < 1.0e-14);
    let diagonal = [0.5, 1.0 / 3.0, 0.25];
    for (axis, expected) in diagonal.iter().enumerate() {
        assert!((stretched.dm_inv[4 * axis] - expected).abs() < 1.0e-15);
    }
    assert_eq!(stretched.dm_inv[1], 0.0);

    // The orientation of the nodes does not change the volume.
    let mirrored =
        TetRest::from_positions([[0.0; 3], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .expect("a mirrored tetrahedron is valid");
    assert!((mirrored.volume - 1.0 / 6.0).abs() < 1.0e-15);
}

#[test]
fn tet_rest_rejects_flat_sliver_and_non_finite_tetrahedra() {
    let flat = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    assert!(TetRest::from_positions(flat).is_none());
    let sliver = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.3, 0.3, 1.0e-14],
    ];
    assert!(TetRest::from_positions(sliver).is_none());
    let coincident = [[0.5; 3]; 4];
    assert!(TetRest::from_positions(coincident).is_none());
    for poison in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut broken = UNIT_TET;
        broken[3][2] = poison;
        assert!(TetRest::from_positions(broken).is_none(), "{poison}");
    }
}

#[test]
fn tet_gradient_reproduces_affine_fields() {
    let map = [1.2, -0.3, 0.1, 0.4, 0.9, -0.2, -0.1, 0.25, 1.1];
    let shift = [0.7, -1.3, 2.1];
    for rest_positions in sample_tetrahedra() {
        let rest = TetRest::from_positions(rest_positions).expect("the samples are valid");
        let field = rest_positions.map(|x| add3(mat3_apply(&map, x), shift));
        let gradient = rest.gradient(&field);
        for (got, expected) in gradient.iter().zip(map.iter()) {
            assert!(
                (got - expected).abs() < 1.0e-9 * (1.0 + expected.abs()),
                "{gradient:?} vs {map:?}"
            );
        }
    }
}

#[test]
fn tet_nodal_forces_balance_force_and_moment_for_a_symmetric_stress() {
    let stress = [3.0, 0.5, -0.2, 0.5, -1.0, 0.7, -0.2, 0.7, 2.0];
    for rest_positions in sample_tetrahedra() {
        let rest = TetRest::from_positions(rest_positions).expect("the samples are valid");
        let forces = rest.nodal_forces(&stress);
        let largest = forces.iter().map(|f| norm3(*f)).fold(0.0, f64::max);
        assert!(largest > 0.0);
        assert!(
            norm3(total_force(&forces)) <= 1.0e-12 * largest,
            "net force {:?} for forces {forces:?}",
            total_force(&forces)
        );
        // Torque about the first node: sum of (x_k - x_0) x f_k.
        let mut torque = [0.0; 3];
        for (position, force) in rest_positions.iter().zip(forces.iter()) {
            let arm = [
                position[0] - rest_positions[0][0],
                position[1] - rest_positions[0][1],
                position[2] - rest_positions[0][2],
            ];
            torque = add3(
                torque,
                [
                    arm[1] * force[2] - arm[2] * force[1],
                    arm[2] * force[0] - arm[0] * force[2],
                    arm[0] * force[1] - arm[1] * force[0],
                ],
            );
        }
        let reach = rest_positions[1..]
            .iter()
            .map(|p| {
                norm3([
                    p[0] - rest_positions[0][0],
                    p[1] - rest_positions[0][1],
                    p[2] - rest_positions[0][2],
                ])
            })
            .fold(0.0, f64::max);
        assert!(
            norm3(torque) <= 1.0e-9 * largest * reach,
            "a symmetric stress exerts no net torque, got {torque:?}"
        );
    }
}

#[test]
fn tet_nodal_forces_are_minus_volume_times_stress_times_shape_gradients() {
    // Dm = I for the unit tetrahedron: the shape-function gradients are the unit
    // vectors, so node k = 1..3 receives -V0 * (column k - 1 of the stress).
    let rest = TetRest::from_positions(UNIT_TET).expect("valid");
    let stress = [3.0, 0.5, -0.2, 0.4, -1.0, 0.7, -0.6, 0.8, 2.0];
    let forces = rest.nodal_forces(&stress);
    let v0 = 1.0 / 6.0;
    for axis in 0..3 {
        let column = [stress[axis], stress[3 + axis], stress[6 + axis]];
        for (component, value) in column.iter().enumerate() {
            assert!(
                (forces[axis + 1][component] + v0 * value).abs() < 1.0e-15,
                "node {} component {component}: {} vs {}",
                axis + 1,
                forces[axis + 1][component],
                -v0 * value
            );
        }
    }
    // The stress is linear in the forces.
    let doubled = rest.nodal_forces(&stress.map(|s| 2.0 * s));
    for (single, twice) in forces.iter().zip(doubled.iter()) {
        for component in 0..3 {
            assert!((2.0 * single[component] - twice[component]).abs() < 1.0e-15);
        }
    }
}

#[test]
fn explicit_stability_limit_matches_the_damped_oscillator_formula() {
    assert!(explicit_stability_limit(0.0, 0.0).is_infinite());
    assert!((explicit_stability_limit(50.0, 0.0) - 2.0 / 50.0).abs() < 1.0e-15);
    assert!((explicit_stability_limit(0.0, 8.0) - 0.25).abs() < 1.0e-15);
    let (omega, gamma) = (50.0_f64, 20.0_f64);
    let expected = 2.0 / ((omega * omega + 0.25 * gamma * gamma).sqrt() + 0.5 * gamma);
    assert!((explicit_stability_limit(omega, gamma) - expected).abs() < 1.0e-15);
    // Damping and stiffness both shorten the admissible step.
    assert!(explicit_stability_limit(omega, gamma) < explicit_stability_limit(omega, 0.0));
    assert!(explicit_stability_limit(omega, gamma) < explicit_stability_limit(0.0, gamma));
    // Negative inputs mean "no stiffness" / "no damping".
    assert!(explicit_stability_limit(-1.0, -1.0).is_infinite());
    assert_eq!(
        explicit_stability_limit(-3.0, 8.0),
        explicit_stability_limit(0.0, 8.0)
    );
}

/// Largest `|x|` reached by `v += h (-w^2 x - g v); x += h v` from `x = 1`.
fn oscillator_peak(h: f64, omega: f64, gamma: f64, steps: usize) -> f64 {
    let (mut x, mut v) = (1.0_f64, 0.0_f64);
    let mut peak = 0.0_f64;
    for _ in 0..steps {
        v += h * (-omega * omega * x - gamma * v);
        x += h * v;
        peak = peak.max(x.abs());
        if !peak.is_finite() {
            return f64::INFINITY;
        }
    }
    peak
}

#[test]
fn stability_limit_separates_bounded_from_diverging_oscillations() {
    // Semi-implicit Euler with the damping evaluated at the old velocity, the
    // scheme the tissue models use.
    for (omega, gamma) in [
        (1.0, 0.0),
        (50.0, 0.0),
        (50.0, 20.0),
        (300.0, 1000.0),
        (10.0, 5000.0),
    ] {
        let limit = explicit_stability_limit(omega, gamma);
        let below = oscillator_peak(0.98 * limit, omega, gamma, 20_000);
        let above = oscillator_peak(1.02 * limit, omega, gamma, 20_000);
        assert!(
            below < 10.0,
            "omega {omega}, gamma {gamma}: 98% of the limit must stay bounded, peak {below}"
        );
        assert!(
            above > 1.0e6,
            "omega {omega}, gamma {gamma}: 102% of the limit must diverge, peak {above}"
        );
    }
}

#[test]
fn stability_substeps_respect_the_safety_margin_and_the_cap() {
    // Nothing limits the step.
    assert_eq!(stability_substeps(0.01, f64::INFINITY), 1);
    // The step already fits into STABILITY_SAFETY * critical.
    assert_eq!(stability_substeps(0.01, 0.04), 1);
    // Smallest number of equal sub-steps below the allowed length.
    assert_eq!(stability_substeps(0.01, 0.0015), 14);
    for (dt, critical) in [(0.02, 0.003), (1.0 / 60.0, 1.0e-4), (0.5, 0.7), (3.0, 0.1)] {
        let substeps = stability_substeps(dt, critical);
        assert!((1..MAX_STABILITY_SUBSTEPS).contains(&substeps));
        let allowed = STABILITY_SAFETY * critical;
        assert!(
            dt / substeps as f64 <= allowed * (1.0 + 1.0e-12),
            "{substeps} sub-steps of {dt} exceed {allowed}"
        );
        assert!(
            substeps == 1 || dt / (substeps - 1) as f64 > allowed,
            "{substeps} sub-steps are more than needed for {dt} / {allowed}"
        );
    }
    // Hopeless cases hit the cap instead of looping for ever.
    assert_eq!(stability_substeps(1.0, 1.0e-9), MAX_STABILITY_SUBSTEPS);
    assert_eq!(stability_substeps(1.0, 0.0), MAX_STABILITY_SUBSTEPS);
    assert_eq!(stability_substeps(1.0, f64::NAN), MAX_STABILITY_SUBSTEPS);
    // Nothing to integrate: a single (empty) step.
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(stability_substeps(dt, 1.0e-9), 1, "dt {dt}");
    }
}

/// Largest eigenvalue of the stiffness matrix of a linear-elastic tetrahedron,
/// by power iteration on `u -> -forces(u)` (the six rigid-body modes have
/// eigenvalue zero and are never amplified).
fn largest_stiffness_eigenvalue(rest: &TetRest, lame_lambda: f64, shear_modulus: f64) -> f64 {
    let stiffness_times = |u: &[[f64; 3]; 4]| -> [[f64; 3]; 4] {
        let stress = linear_elastic_stress(&rest.gradient(u), lame_lambda, shear_modulus, 0.0);
        rest.nodal_forces(&stress).map(|f| scale3(f, -1.0))
    };
    let dot = |a: &[[f64; 3]; 4], b: &[[f64; 3]; 4]| -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| x[0] * y[0] + x[1] * y[1] + x[2] * y[2])
            .sum()
    };
    let mut u = [[0.0; 3]; 4];
    for (k, node) in u.iter_mut().enumerate() {
        for (axis, component) in node.iter_mut().enumerate() {
            *component = ((1 + 3 * k + axis) as f64 * 0.7).sin() + 0.1 * k as f64;
        }
    }
    let mut estimate = 0.0;
    for _ in 0..4000 {
        let ku = stiffness_times(&u);
        let norm = dot(&ku, &ku).sqrt();
        if norm == 0.0 {
            break;
        }
        estimate = dot(&u, &ku) / dot(&u, &u);
        u = ku.map(|node| scale3(node, 1.0 / norm));
    }
    estimate
}

#[test]
fn stiffness_row_sums_match_the_blocks_of_the_force_law() {
    // Block (a, b) of the stiffness matrix: column c is minus the force on node
    // a for a unit displacement of node b along axis c.
    let (lame_lambda, shear_modulus) = (7.0e5, 9.0e4);
    for rest_positions in sample_tetrahedra() {
        let rest = TetRest::from_positions(rest_positions).expect("the samples are valid");
        let mut expected = [0.0; 4];
        for b in 0..4 {
            let mut squares = [0.0; 4];
            for axis in 0..3 {
                let mut u = [[0.0; 3]; 4];
                u[b][axis] = 1.0;
                let stress =
                    linear_elastic_stress(&rest.gradient(&u), lame_lambda, shear_modulus, 0.0);
                for (square, force) in squares.iter_mut().zip(rest.nodal_forces(&stress)) {
                    *square += dot3(force, force);
                }
            }
            for (total, square) in expected.iter_mut().zip(squares) {
                *total += f64::sqrt(square);
            }
        }
        let sums = rest.stiffness_row_sums(lame_lambda, shear_modulus);
        for (got, want) in sums.iter().zip(expected.iter()) {
            assert!(
                (got - want).abs() <= 1.0e-9 * want,
                "{rest_positions:?}: {sums:?} vs {expected:?}"
            );
        }
    }
}

#[test]
fn stiffness_row_sums_dominate_the_largest_eigenvalue_of_the_element_matrix() {
    for (young, poisson) in [(2.0e5, 0.45), (1.0e4, 0.3), (5.0e3, 0.0), (3.0e3, 0.49)] {
        let layer = TissueLayer::new("probe", young, poisson, 0.01, 1000.0, 100.0);
        let (lame_lambda, shear_modulus) = (layer.lame_lambda(), layer.shear_modulus());
        for rest_positions in sample_tetrahedra() {
            let rest = TetRest::from_positions(rest_positions).expect("the samples are valid");
            let bound = rest
                .stiffness_row_sums(lame_lambda, shear_modulus)
                .into_iter()
                .fold(0.0, f64::max);
            let eigenvalue = largest_stiffness_eigenvalue(&rest, lame_lambda, shear_modulus);
            assert!(eigenvalue > 0.0, "E {young}, nu {poisson}");
            assert!(
                eigenvalue <= bound * (1.0 + 1.0e-9),
                "E {young}, nu {poisson}, {rest_positions:?}: the bound {bound} is below the \
                 eigenvalue {eigenvalue}"
            );
            assert!(
                bound <= 4.0 * eigenvalue,
                "E {young}, nu {poisson}, {rest_positions:?}: the bound {bound} is far too loose \
                 against the eigenvalue {eigenvalue}"
            );
        }
    }
}

#[test]
fn power_iteration_finds_the_known_eigenvalue_of_the_unit_tetrahedron() {
    // E = 2e5, nu = 0.45 on the unit tetrahedron: lambda_max = 691293.0409...
    // (dense symmetric eigensolver of the 12 x 12 stiffness matrix).
    let layer = TissueLayer::new("skin", 2.0e5, 0.45, 0.002, 1100.0, 500.0);
    let rest = TetRest::from_positions(UNIT_TET).expect("valid");
    let eigenvalue =
        largest_stiffness_eigenvalue(&rest, layer.lame_lambda(), layer.shear_modulus());
    assert!(
        (eigenvalue - 691_293.040_907).abs() < 1.0e-6 * 691_293.040_907,
        "{eigenvalue}"
    );
}

// ---------------------------------------------------------------------------
// Robustness: elements that cannot carry load
// ---------------------------------------------------------------------------

/// Skin tissue with one proper element (nodes 0-3) next to three that cannot
/// carry load: a flat rest shape (nodes 4-7), a tetrahedron of a layer with an
/// infinite Lame parameter (nodes 8-11) and an element that refers to a node
/// that does not exist.  Everything is deformed, so that an element that was
/// wrongly included would exert a force.
fn tissue_with_unusable_elements() -> TissueModel {
    let incompressible = TissueLayer::new("incompressible", 1.0e4, 0.5, 0.01, 1000.0, 100.0);
    let mut model = TissueModel::new(vec![skin(), incompressible], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    let corners = |origin: [f64; 3]| {
        [
            origin,
            [origin[0] + 1.0, origin[1], origin[2]],
            [origin[0], origin[1] + 1.0, origin[2]],
            [origin[0], origin[1], origin[2] + 1.0],
        ]
    };
    for position in corners([0.0, 0.0, 0.0]) {
        model.add_node(TissueNode::new(position, 0.1, 0));
    }
    // Coplanar in the rest configuration.
    for position in [
        [5.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [5.0, 1.0, 0.0],
        [6.0, 1.0, 0.0],
    ] {
        model.add_node(TissueNode::new(position, 0.1, 0));
    }
    for position in corners([9.0, 0.0, 0.0]) {
        model.add_node(TissueNode::new(position, 0.1, 1));
    }
    model.add_element([0, 1, 2, 3]);
    model.add_element([4, 5, 6, 7]);
    model.add_element([8, 9, 10, 11]);
    model.add_element([0, 1, 2, 99]);
    model.nodes[1].position[0] = 1.1;
    model.nodes[7].position[2] = 0.3;
    model.nodes[10].position[1] = 1.2;
    model
}

/// The proper element of [`tissue_with_unusable_elements`] on its own.
fn tissue_with_the_proper_element_only() -> TissueModel {
    let mut model = TissueModel::new(vec![skin()], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    for position in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ] {
        model.add_node(TissueNode::new(position, 0.1, 0));
    }
    model.add_element([0, 1, 2, 3]);
    model.nodes[1].position[0] = 1.1;
    model
}

#[test]
fn tissue_skips_elements_without_nodes_rest_shape_or_finite_material() {
    let mut model = tissue_with_unusable_elements();
    let reference = tissue_with_the_proper_element_only();

    let forces = model.internal_forces();
    let expected = reference.internal_forces();
    assert!(norm3(expected[1]) > 0.0, "the proper element is loaded");
    assert_eq!(forces.len(), 12);
    assert_eq!(forces[..4], expected[..]);
    for (node, force) in forces.iter().enumerate().skip(4) {
        assert_eq!(*force, [0.0; 3], "node {node} belongs to no usable element");
    }
    assert_eq!(model.elastic_energy(), reference.elastic_energy());
    let critical = model.critical_time_step();
    assert!(critical.is_finite());
    assert_eq!(critical, reference.critical_time_step());

    // The nodes of the unusable elements are in free fall.
    let dt = 1.0e-3;
    model.step(dt, GRAVITY);
    for node in 4..12 {
        let velocity = model.nodes[node].velocity;
        assert!(
            (velocity[1] + G * dt).abs() < 1.0e-12 && velocity[0] == 0.0 && velocity[2] == 0.0,
            "node {node}: {velocity:?}"
        );
    }
    assert!(
        model
            .nodes
            .iter()
            .all(|n| n.position.iter().all(|c| c.is_finite()))
    );
}

#[test]
fn tissue_without_layers_has_no_material_and_exerts_no_elastic_force() {
    let mut model = TissueModel::new(Vec::new(), 0.0);
    model.damping_alpha = 0.0;
    for position in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ] {
        model.add_node(TissueNode::new(position, 0.1, 0));
    }
    model.add_element([0, 1, 2, 3]);
    model.nodes[1].position[0] = 1.3;

    assert!(model.internal_forces().iter().all(|f| *f == [0.0; 3]));
    assert_eq!(model.elastic_energy(), 0.0);
    assert_eq!(model.elastic_energy_density(0, &model.rest_positions), 0.0);
    assert!(model.critical_time_step().is_infinite());
    model.step(0.01, GRAVITY);
    for node in &model.nodes {
        assert!((node.velocity[1] + G * 0.01).abs() < 1.0e-12);
    }
}

#[test]
fn tissue_step_of_zero_duration_keeps_the_state_and_clears_the_forces() {
    let mut model = unit_tet_tissue(0.01);
    model.nodes[1].position[0] = 1.1;
    model.nodes[2].force = [3.0, 0.0, 0.0];
    let before: Vec<[f64; 3]> = model.nodes.iter().map(|n| n.position).collect();
    model.step(0.0, GRAVITY);
    for (node, position) in model.nodes.iter().zip(before.iter()) {
        assert_eq!(node.position, *position);
        assert_eq!(node.velocity, [0.0; 3]);
        assert_eq!(node.force, [0.0; 3]);
    }
}

#[test]
fn tissue_strong_stiffness_damping_does_not_destabilise_the_explicit_step() {
    // beta * omega^2 is orders of magnitude above 1 / dt: the damping term
    // limits the step long before the stiffness does.
    let mut model = unit_tet_tissue(0.01);
    model.damping_beta = 0.005;
    model.nodes[1].position[0] = 1.1;
    let initial = model.elastic_energy();
    let critical = model.critical_time_step();
    assert!(critical < 1.0e-5, "critical dt = {critical}");
    for _ in 0..50 {
        model.step(1.0e-3, NO_GRAVITY);
    }
    assert!(model.nodes.iter().all(|n| {
        n.position
            .iter()
            .chain(n.velocity.iter())
            .all(|c| c.is_finite())
    }));
    // The damping makes the node creep back to its rest position with the time
    // constant beta = 5 ms: after ten time constants it has settled.
    let offset = (model.nodes[1].position[0] - 1.0).abs();
    assert!(offset < 1.0e-3 * 0.1, "offset {offset}");
    assert!(model.elastic_energy() < 1.0e-3 * initial);
}

#[test]
fn soft_tissue_skips_elements_without_nodes_rest_shape_or_finite_material() {
    let good = NeoHookeanParams::from_young(1.0e4, 0.3, 0.5);
    let mut broken_material = good.clone();
    broken_material.mu = f64::NAN;

    let mut model = SoftTissueModel::new(0.0);
    model.add_node([0.0, 0.0, 0.0], 0.0, true);
    model.add_node([1.0, 0.0, 0.0], 0.01, false);
    model.add_node([0.0, 1.0, 0.0], 0.01, false);
    model.add_node([0.0, 0.0, 1.0], 0.01, false);
    // Coplanar in the rest configuration.
    model.add_node([5.0, 0.0, 0.0], 0.01, false);
    model.add_node([6.0, 0.0, 0.0], 0.01, false);
    model.add_node([5.0, 1.0, 0.0], 0.01, false);
    model.add_node([6.0, 1.0, 0.0], 0.01, false);
    model.add_element([0, 1, 2, 3], good.clone());
    model.add_element([4, 5, 6, 7], good.clone());
    model.add_element([0, 1, 2, 99], good.clone());
    model.add_element([0, 1, 2, 3], broken_material);
    model.nodes[1][0] = 1.1;
    model.nodes[7][2] = 0.3;
    model.velocities[1] = [0.2, 0.0, 0.0];

    let mut reference = neo_tet(0.01, 0.0, 0.5);
    reference.nodes[1][0] = 1.1;
    reference.velocities[1] = [0.2, 0.0, 0.0];

    let forces = model.internal_forces();
    let expected = reference.internal_forces();
    assert!(norm3(expected[1]) > 0.0, "the proper element is loaded");
    assert_eq!(forces.len(), 8);
    assert_eq!(forces[..4], expected[..]);
    for (node, force) in forces.iter().enumerate().skip(4) {
        assert_eq!(*force, [0.0; 3], "node {node} belongs to no usable element");
    }
    let energy = model.strain_energy();
    assert!(energy.is_finite() && energy > 0.0);
    assert_eq!(energy, reference.strain_energy());
    let critical = model.critical_time_step();
    assert!(critical.is_finite());
    assert_eq!(critical, reference.critical_time_step());

    let dt = 1.0e-3;
    model.step(dt, GRAVITY);
    for node in 4..8 {
        let velocity = model.velocities[node];
        assert!(
            (velocity[1] + G * dt).abs() < 1.0e-12 && velocity[0] == 0.0 && velocity[2] == 0.0,
            "node {node}: {velocity:?}"
        );
    }
    assert!(model.nodes.iter().all(|n| n.iter().all(|c| c.is_finite())));
}

#[test]
fn soft_tissue_ignores_elements_that_outrun_the_per_node_arrays() {
    // A caller that removed a velocity entry by hand leaves node 3 without
    // state; the element that needs it must not be evaluated.
    let mut model = neo_tet(0.01, 0.0, 0.0);
    model.nodes[1][0] = 1.1;
    model.velocities.pop();
    assert!(model.internal_forces().iter().all(|f| *f == [0.0; 3]));
    assert_eq!(model.strain_energy(), 0.0);
    model.step(0.01, GRAVITY);
    assert!(model.nodes.iter().all(|n| n.iter().all(|c| c.is_finite())));
}

#[test]
fn soft_tissue_critical_time_step_shrinks_with_stiffness_and_viscosity() {
    let build = |young: f64, eta: f64, alpha: f64| {
        let mut model = SoftTissueModel::new(alpha);
        let material = NeoHookeanParams::from_young(young, 0.3, eta);
        model.add_node([0.0, 0.0, 0.0], 0.0, true);
        model.add_node([1.0, 0.0, 0.0], 0.1, false);
        model.add_node([0.0, 1.0, 0.0], 0.1, false);
        model.add_node([0.0, 0.0, 1.0], 0.1, false);
        model.add_element([0, 1, 2, 3], material);
        model
    };
    let soft = build(1.0e4, 0.0, 0.0).critical_time_step();
    let stiff = build(4.0e4, 0.0, 0.0).critical_time_step();
    assert!(
        (soft / stiff - 2.0).abs() < 1.0e-9,
        "dt_crit ~ 1/sqrt(E): {soft} vs {stiff}"
    );
    assert!(build(1.0e4, 50.0, 0.0).critical_time_step() < soft);
    assert!(build(1.0e4, 0.0, 100.0).critical_time_step() < soft);
    // Mass damping alone limits the step to 2 / alpha.
    let mut free = SoftTissueModel::new(4.0);
    free.add_node([0.0; 3], 1.0, false);
    assert!((free.critical_time_step() - 0.5).abs() < 1.0e-12);
    free.alpha = 0.0;
    assert!(free.critical_time_step().is_infinite());
}

// ---------------------------------------------------------------------------
// Whole meshes: the critical time step bounds the natural frequencies
// ---------------------------------------------------------------------------

/// Nodes of a regular grid of `n` x `n` x `n` points with spacing `edge`
/// (index `i + n (j + n k)`).
fn block_positions(n: usize, edge: f64) -> Vec<[f64; 3]> {
    let mut positions = Vec::with_capacity(n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                positions.push([i as f64 * edge, j as f64 * edge, k as f64 * edge]);
            }
        }
    }
    positions
}

/// Six tetrahedra per cube around its main diagonal; neighbouring cubes agree
/// on the diagonal of their shared face, so the mesh is conforming.
fn block_tetrahedra(n: usize) -> Vec<[usize; 4]> {
    let index = |i: usize, j: usize, k: usize| i + n * (j + n * k);
    let mut tetrahedra = Vec::new();
    for k in 0..n - 1 {
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let cube = [
                    index(i, j, k),
                    index(i + 1, j, k),
                    index(i + 1, j + 1, k),
                    index(i, j + 1, k),
                    index(i, j, k + 1),
                    index(i + 1, j, k + 1),
                    index(i + 1, j + 1, k + 1),
                    index(i, j + 1, k + 1),
                ];
                for corners in [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    tetrahedra.push(corners.map(|corner| cube[corner]));
                }
            }
        }
    }
    tetrahedra
}

/// Node mass of the block: about one cell of water, varied from node to node so
/// that the bound has to use the individual masses.
fn block_node_mass(index: usize, edge: f64) -> f64 {
    1000.0 * edge.powi(3) * (1.0 + 0.5 * ((7 * index) % 5) as f64 / 4.0)
}

/// A skin block whose face `y = 0` is pinned; no damping.
fn tissue_block(n: usize, edge: f64) -> TissueModel {
    let mut model = TissueModel::new(vec![skin()], 0.0);
    model.damping_alpha = 0.0;
    model.damping_beta = 0.0;
    for (index, position) in block_positions(n, edge).into_iter().enumerate() {
        if position[1] == 0.0 {
            model.add_node(TissueNode::new_static(position, 0));
        } else {
            model.add_node(TissueNode::new(position, block_node_mass(index, edge), 0));
        }
    }
    for tetrahedron in block_tetrahedra(n) {
        model.add_element(tetrahedron);
    }
    model
}

/// A neo-Hookean block whose face `y = 0` is pinned.
fn soft_tissue_block(n: usize, edge: f64, material: &NeoHookeanParams) -> SoftTissueModel {
    let mut model = SoftTissueModel::new(0.0);
    for (index, position) in block_positions(n, edge).into_iter().enumerate() {
        let fixed = position[1] == 0.0;
        model.add_node(position, block_node_mass(index, edge), fixed);
    }
    for tetrahedron in block_tetrahedra(n) {
        model.add_element(tetrahedron, material.clone());
    }
    model
}

/// Largest eigenvalue of `M^-1 A` by power iteration, as the Rayleigh quotient
/// `u'A u / u'M u`.  `apply(u)` returns `A u` (zero for constrained nodes) and
/// `masses` holds the lumped masses, zero for constrained nodes.  The iterate is
/// rescaled to the given largest component `amplitude` every time.
fn largest_generalized_eigenvalue(
    masses: &[f64],
    mut apply: impl FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
    amplitude: f64,
) -> f64 {
    let rescale = |vector: &mut Vec<[f64; 3]>| {
        let peak = vector
            .iter()
            .flat_map(|node| node.iter())
            .fold(0.0_f64, |acc, c| acc.max(c.abs()));
        if peak > 0.0 {
            for node in vector.iter_mut() {
                *node = scale3(*node, amplitude / peak);
            }
        }
    };
    let mut u: Vec<[f64; 3]> = masses
        .iter()
        .enumerate()
        .map(|(i, &mass)| {
            if mass > 0.0 {
                let t = i as f64;
                [(1.7 * t).sin(), (2.3 * t).cos(), (0.9 * t + 0.4).sin()]
            } else {
                [0.0; 3]
            }
        })
        .collect();
    rescale(&mut u);
    let mut estimate = 0.0;
    for _ in 0..400 {
        let au = apply(&u);
        let work: f64 = u.iter().zip(au.iter()).map(|(x, y)| dot3(*x, *y)).sum();
        let inertia: f64 = u
            .iter()
            .zip(masses.iter())
            .map(|(x, &m)| m * dot3(*x, *x))
            .sum();
        estimate = work / inertia;
        let mut next: Vec<[f64; 3]> = au
            .iter()
            .zip(masses.iter())
            .map(|(a, &m)| {
                if m > 0.0 {
                    scale3(*a, 1.0 / m)
                } else {
                    [0.0; 3]
                }
            })
            .collect();
        rescale(&mut next);
        u = next;
    }
    estimate
}

/// Lumped masses of a tissue block, zero for pinned nodes.
fn free_masses(masses: &[f64], fixed: impl Iterator<Item = bool>) -> Vec<f64> {
    masses
        .iter()
        .zip(fixed)
        .map(|(&mass, fixed)| if fixed { 0.0 } else { mass })
        .collect()
}

#[test]
fn tissue_block_critical_time_step_bounds_the_natural_frequencies() {
    let (n, edge) = (4, 0.02);
    let mut model = tissue_block(n, edge);
    let critical = model.critical_time_step();
    assert!(critical.is_finite() && critical > 0.0);
    let omega_sq_bound = (2.0 / critical).powi(2);

    let rest = model.rest_positions.clone();
    let masses = free_masses(
        &model.nodes.iter().map(|node| node.mass).collect::<Vec<_>>(),
        model.nodes.iter().map(|node| node.is_static),
    );
    let omega_sq = largest_generalized_eigenvalue(
        &masses,
        |u| {
            for (node, (rest_position, displacement)) in
                model.nodes.iter_mut().zip(rest.iter().zip(u.iter()))
            {
                node.position = add3(*rest_position, *displacement);
            }
            let forces = model.internal_forces();
            forces
                .iter()
                .zip(masses.iter())
                .map(|(f, &m)| if m > 0.0 { scale3(*f, -1.0) } else { [0.0; 3] })
                .collect()
        },
        1.0e-3 * edge,
    );
    assert!(omega_sq > 0.0);
    assert!(
        omega_sq_bound >= omega_sq * (1.0 - 1.0e-6),
        "the bound {omega_sq_bound} is below the largest squared frequency {omega_sq}"
    );
    assert!(
        omega_sq_bound <= 8.0 * omega_sq,
        "the bound {omega_sq_bound} is far too loose against {omega_sq}"
    );
}

#[test]
fn soft_tissue_block_critical_time_step_bounds_the_natural_frequencies() {
    let (n, edge) = (4, 0.02);
    let material = NeoHookeanParams::from_young(2.0e5, 0.3, 0.0);
    let mut model = soft_tissue_block(n, edge, &material);
    let critical = model.critical_time_step();
    assert!(critical.is_finite() && critical > 0.0);
    let omega_sq_bound = (2.0 / critical).powi(2);

    let rest = model.rest_positions.clone();
    let masses = free_masses(&model.masses, model.fixed.iter().copied());
    // Displacements of one millionth of an edge stay in the linear regime.
    let omega_sq = largest_generalized_eigenvalue(
        &masses,
        |u| {
            for (position, (rest_position, displacement)) in
                model.nodes.iter_mut().zip(rest.iter().zip(u.iter()))
            {
                *position = add3(*rest_position, *displacement);
            }
            let forces = model.internal_forces();
            forces
                .iter()
                .zip(masses.iter())
                .map(|(f, &m)| if m > 0.0 { scale3(*f, -1.0) } else { [0.0; 3] })
                .collect()
        },
        1.0e-6 * edge,
    );
    assert!(omega_sq > 0.0);
    assert!(
        omega_sq_bound >= omega_sq * (1.0 - 1.0e-4),
        "the bound {omega_sq_bound} is below the largest squared frequency {omega_sq}"
    );
    assert!(
        omega_sq_bound <= 8.0 * omega_sq,
        "the bound {omega_sq_bound} is far too loose against {omega_sq}"
    );
}

#[test]
fn soft_tissue_block_viscous_bound_dominates_the_damping_rate() {
    // A purely viscous material (no elasticity): the critical step is 2 / gamma
    // for the largest damping rate gamma of M^-1 C, which is the largest
    // eigenvalue of M^-1 C for velocity fields at the rest configuration.
    let (n, edge) = (4, 0.02);
    let material = NeoHookeanParams {
        mu: 0.0,
        kappa: 0.0,
        eta: 40.0,
    };
    let mut model = soft_tissue_block(n, edge, &material);
    let critical = model.critical_time_step();
    assert!(critical.is_finite() && critical > 0.0);
    let gamma_bound = 2.0 / critical;

    let masses = free_masses(&model.masses, model.fixed.iter().copied());
    let gamma = largest_generalized_eigenvalue(
        &masses,
        |u| {
            model.velocities.clone_from_slice(u);
            let forces = model.internal_forces();
            forces
                .iter()
                .zip(masses.iter())
                .map(|(f, &m)| if m > 0.0 { scale3(*f, -1.0) } else { [0.0; 3] })
                .collect()
        },
        1.0,
    );
    assert!(gamma > 0.0);
    assert!(
        gamma_bound >= gamma * (1.0 - 1.0e-6),
        "the bound {gamma_bound} is below the largest damping rate {gamma}"
    );
    assert!(
        gamma_bound <= 8.0 * gamma,
        "the bound {gamma_bound} is far too loose against {gamma}"
    );
}

#[test]
fn tissue_block_with_default_damping_is_stable_at_sixty_hertz() {
    // A mesh of centimetre-sized elements with the default Rayleigh
    // coefficients: the guaranteed step is reached without hitting the sub-step
    // cap.
    let (n, edge) = (3, 0.04);
    let mut model = tissue_block(n, edge);
    model.damping_alpha = 0.02;
    model.damping_beta = 0.002;
    let dt = 1.0 / 60.0;
    let critical = model.critical_time_step();
    let substeps = dt / (STABILITY_SAFETY * critical);
    assert!(
        substeps < MAX_STABILITY_SUBSTEPS as f64,
        "{substeps} sub-steps per frame would exceed the cap"
    );
    for _ in 0..6 {
        model.step(dt, GRAVITY);
    }
    let extent = edge * (n - 1) as f64;
    for node in &model.nodes {
        assert!(
            node.position
                .iter()
                .all(|c| c.is_finite() && *c > -0.1 * extent && *c < 1.1 * extent),
            "{:?}",
            node.position
        );
        assert!(norm3(node.velocity) < 1.0, "{:?}", node.velocity);
    }
    // The pinned face keeps its place and the free nodes sag downwards.
    let sagging = model
        .nodes
        .iter()
        .zip(model.rest_positions.iter())
        .filter(|(node, _)| !node.is_static)
        .all(|(node, rest)| node.position[1] <= rest[1] + 1.0e-12);
    assert!(sagging);
}
