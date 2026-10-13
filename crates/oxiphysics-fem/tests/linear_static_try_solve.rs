// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for [`LinearStaticAnalysis::try_solve`].
//!
//! `LinearStaticAnalysis::solve` used to be infallible: it returned whatever
//! the conjugate-gradient iteration produced, including the garbage of a
//! singular (unconstrained) system, a solve that never converged, or NaN from a
//! degenerate element.  `try_solve` validates the problem and reports those
//! failures.

use oxiphysics_core::math::Vec3;
use oxiphysics_fem::Error;
use oxiphysics_fem::analysis::{FemError, FemResult, LinearStaticAnalysis};
use oxiphysics_fem::boundary::{DirichletBc, NeumannBc};
use oxiphysics_fem::constitutive::LinearElasticMaterial;
use oxiphysics_fem::mesh::TetrahedralMesh;

const LENGTH: f64 = 1.0;
const WIDTH: f64 = 0.1;
const HEIGHT: f64 = 0.1;

fn steel() -> LinearElasticMaterial {
    LinearElasticMaterial::new(200.0e9, 0.3)
}

fn zero_force() -> Vec3 {
    Vec3::new(0.0, 0.0, 0.0)
}

/// Beam clamped at x = 0 with a downward tip load at x = L.
fn clamped_beam() -> (
    TetrahedralMesh,
    Vec<DirichletBc>,
    Vec<NeumannBc>,
    Vec<usize>,
) {
    let mesh = TetrahedralMesh::generate_beam(LENGTH, WIDTH, HEIGHT, 4, 1, 1);
    let mut clamp = Vec::new();
    for (i, node) in mesh.nodes.iter().enumerate() {
        if node.x.abs() < 1e-10 {
            for dof in 0..3 {
                clamp.push(DirichletBc::new(i, dof, 0.0));
            }
        }
    }
    let tip: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - LENGTH).abs() < 1e-10)
        .map(|(i, _)| i)
        .collect();
    let per_node = -500.0 / tip.len() as f64;
    let loads = tip
        .iter()
        .map(|&i| NeumannBc::new(i, [0.0, per_node, 0.0]))
        .collect();
    (mesh, clamp, loads, tip)
}

fn analysis() -> LinearStaticAnalysis {
    LinearStaticAnalysis {
        max_iter: 20_000,
        tolerance: 1e-10,
    }
}

fn assert_invalid_input<T: std::fmt::Debug>(result: Result<T, FemError>, what: &str) {
    assert!(
        matches!(result, Err(FemError::InvalidInput(_))),
        "{what}: expected InvalidInput, got {result:?}"
    );
}

#[test]
fn clamped_cantilever_solves_with_finite_downward_deflection() {
    let (mesh, clamp, loads, tip) = clamped_beam();
    let result = analysis()
        .try_solve(&mesh, &steel(), &clamp, &loads, &zero_force())
        .expect("a clamped beam is a well-posed problem");
    assert_eq!(result.displacements.len(), mesh.num_nodes());
    assert_eq!(result.stress.len(), mesh.num_elements());
    assert_eq!(result.strain.len(), mesh.num_elements());
    let tip_y = tip
        .iter()
        .map(|&i| result.displacements[i].y)
        .fold(f64::INFINITY, f64::min);
    assert!(tip_y < 0.0 && tip_y.is_finite(), "tip deflection {tip_y}");
    for bc in &clamp {
        assert_eq!(
            result.displacements[bc.node].norm(),
            0.0,
            "clamped node {} moved",
            bc.node
        );
    }
}

#[test]
fn solve_returns_exactly_the_try_solve_result_on_success() {
    let (mesh, clamp, loads, _) = clamped_beam();
    let strict = analysis()
        .try_solve(&mesh, &steel(), &clamp, &loads, &zero_force())
        .expect("well posed");
    let lenient = analysis().solve(&mesh, &steel(), &clamp, &loads, &zero_force());
    assert_eq!(strict.displacements, lenient.displacements);
    assert_eq!(strict.stress, lenient.stress);
    assert_eq!(strict.strain, lenient.strain);
}

#[test]
fn unconstrained_structure_under_gravity_is_an_error() {
    // No Dirichlet BCs: K has six rigid-body modes and the gravity resultant is
    // not balanced, so no displacement field exists.
    let mesh = TetrahedralMesh::generate_beam(LENGTH, WIDTH, HEIGHT, 2, 1, 1);
    let gravity = Vec3::new(0.0, -9.81 * 7800.0, 0.0);
    let result = analysis().try_solve(&mesh, &steel(), &[], &[], &gravity);
    assert_eq!(result.as_ref().err(), Some(&FemError::SolverNotConverged));
}

#[test]
fn solve_flags_a_failed_analysis_with_nan_instead_of_plausible_numbers() {
    let mesh = TetrahedralMesh::generate_beam(LENGTH, WIDTH, HEIGHT, 2, 1, 1);
    let gravity = Vec3::new(0.0, -9.81 * 7800.0, 0.0);
    let result: FemResult = analysis().solve(&mesh, &steel(), &[], &[], &gravity);
    assert_eq!(result.displacements.len(), mesh.num_nodes());
    assert_eq!(result.stress.len(), mesh.num_elements());
    assert_eq!(result.strain.len(), mesh.num_elements());
    assert!(
        result
            .displacements
            .iter()
            .all(|d| d.x.is_nan() && d.y.is_nan() && d.z.is_nan())
    );
    assert!(result.stress.iter().all(|s| s.iter().all(|c| c.is_nan())));
    assert!(result.strain.iter().all(|s| s.iter().all(|c| c.is_nan())));
}

#[test]
fn exhausted_iteration_budget_is_reported_as_not_converged() {
    let (mesh, clamp, loads, _) = clamped_beam();
    let starved = LinearStaticAnalysis {
        max_iter: 1,
        tolerance: 1e-12,
    };
    let result = starved.try_solve(&mesh, &steel(), &clamp, &loads, &zero_force());
    assert_eq!(result.as_ref().err(), Some(&FemError::SolverNotConverged));
}

#[test]
fn empty_mesh_is_invalid_input() {
    let mesh = TetrahedralMesh {
        nodes: Vec::new(),
        elements: Vec::new(),
    };
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &[], &[], &zero_force()),
        "empty mesh",
    );
}

#[test]
fn element_referencing_a_missing_node_is_invalid_input() {
    let mut mesh = TetrahedralMesh::generate_beam(LENGTH, WIDTH, HEIGHT, 1, 1, 1);
    let missing = mesh.num_nodes() + 3;
    mesh.elements[0][2] = missing;
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &[], &[], &zero_force()),
        "dangling element node",
    );
}

#[test]
fn degenerate_element_is_invalid_input() {
    // Four coplanar points: zero volume, the B-matrix would divide by zero.
    let nodes = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
    ];
    let mesh = TetrahedralMesh {
        nodes,
        elements: vec![[0, 1, 2, 3]],
    };
    let fixed = [DirichletBc::new(0, 0, 0.0)];
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &fixed, &[], &zero_force()),
        "flat tetrahedron",
    );

    // A sliver whose volume is far below 1e-12 h^3 is degenerate as well ...
    let sliver = TetrahedralMesh {
        nodes: vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.3, 0.3, 1.0e-15),
        ],
        elements: vec![[0, 1, 2, 3]],
    };
    assert_invalid_input(
        analysis().try_solve(&sliver, &steel(), &fixed, &[], &zero_force()),
        "sliver tetrahedron",
    );
    // ... while a merely small but healthy element is not.
    let small = TetrahedralMesh::generate_beam(1.0e-4, 1.0e-5, 1.0e-5, 2, 1, 1);
    let small_clamp: Vec<DirichletBc> = small
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.x.abs() < 1e-12)
        .flat_map(|(i, _)| (0..3).map(move |dof| DirichletBc::new(i, dof, 0.0)))
        .collect();
    let tiny = analysis().try_solve(
        &small,
        &steel(),
        &small_clamp,
        &[],
        &Vec3::new(0.0, -1.0, 0.0),
    );
    assert!(tiny.is_ok(), "a small element is not degenerate: {tiny:?}");
}

#[test]
fn invalid_material_is_invalid_input() {
    let (mesh, clamp, loads, _) = clamped_beam();
    let cases = [
        (0.0, 0.3, "zero Young's modulus"),
        (-1.0e9, 0.3, "negative Young's modulus"),
        (f64::NAN, 0.3, "NaN Young's modulus"),
        (f64::INFINITY, 0.3, "infinite Young's modulus"),
        (200.0e9, 0.5, "incompressible Poisson ratio"),
        (200.0e9, 0.7, "Poisson ratio above 1/2"),
        (200.0e9, -1.0, "Poisson ratio of -1"),
        (200.0e9, f64::NAN, "NaN Poisson ratio"),
    ];
    for (young, poisson, what) in cases {
        let material = LinearElasticMaterial::new(young, poisson);
        assert_invalid_input(
            analysis().try_solve(&mesh, &material, &clamp, &loads, &zero_force()),
            what,
        );
    }
}

#[test]
fn out_of_range_boundary_conditions_are_invalid_input() {
    let (mesh, clamp, loads, _) = clamped_beam();
    let n = mesh.num_nodes();

    let mut bad_node = clamp.clone();
    bad_node.push(DirichletBc::new(n, 0, 0.0));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &bad_node, &loads, &zero_force()),
        "Dirichlet node out of range",
    );

    // Local DOF 3 would silently alias DOF 0 of the next node.
    let mut bad_dof = clamp.clone();
    bad_dof.push(DirichletBc::new(0, 3, 0.0));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &bad_dof, &loads, &zero_force()),
        "Dirichlet local DOF out of range",
    );

    let mut bad_load = loads.clone();
    bad_load.push(NeumannBc::new(n + 10, [0.0, 1.0, 0.0]));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &clamp, &bad_load, &zero_force()),
        "Neumann node out of range",
    );
}

#[test]
fn non_finite_loads_and_prescribed_values_are_invalid_input() {
    let (mesh, clamp, loads, _) = clamped_beam();

    let mut nan_value = clamp.clone();
    nan_value.push(DirichletBc::new(1, 1, f64::NAN));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &nan_value, &loads, &zero_force()),
        "NaN prescribed displacement",
    );

    let mut inf_load = loads.clone();
    inf_load.push(NeumannBc::new(0, [f64::INFINITY, 0.0, 0.0]));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &clamp, &inf_load, &zero_force()),
        "infinite nodal force",
    );

    assert_invalid_input(
        analysis().try_solve(
            &mesh,
            &steel(),
            &clamp,
            &loads,
            &Vec3::new(0.0, f64::NAN, 0.0),
        ),
        "NaN body force",
    );
}

#[test]
fn conflicting_prescribed_values_are_invalid_input() {
    let (mesh, mut clamp, loads, _) = clamped_beam();
    let first = clamp[0].clone();
    clamp.push(DirichletBc::new(first.node, first.dof, 0.01));
    assert_invalid_input(
        analysis().try_solve(&mesh, &steel(), &clamp, &loads, &zero_force()),
        "two different values for one DOF",
    );
}

#[test]
fn invalid_solver_tolerance_is_invalid_input() {
    let (mesh, clamp, loads, _) = clamped_beam();
    for tolerance in [0.0, -1.0e-10, f64::NAN, f64::INFINITY] {
        let solver = LinearStaticAnalysis {
            max_iter: 1000,
            tolerance,
        };
        assert_invalid_input(
            solver.try_solve(&mesh, &steel(), &clamp, &loads, &zero_force()),
            "solver tolerance",
        );
    }
}

#[test]
fn default_solver_parameters_converge_on_a_small_cantilever() {
    // `LinearStaticAnalysis::new()` must be usable out of the box: a 4 x 1 x 1
    // steel cantilever (1 m x 0.25 m x 0.25 m) under a tip load converges within
    // the default iteration budget, and the result agrees with a tightly
    // converged reference solve.
    let (length, width, height) = (1.0_f64, 0.25_f64, 0.25_f64);
    let mesh = TetrahedralMesh::generate_beam(length, width, height, 4, 1, 1);
    let clamp: Vec<DirichletBc> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.x.abs() < 1e-10)
        .flat_map(|(i, _)| (0..3).map(move |dof| DirichletBc::new(i, dof, 0.0)))
        .collect();
    let tip: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - length).abs() < 1e-10)
        .map(|(i, _)| i)
        .collect();
    let load = -1000.0;
    let loads: Vec<NeumannBc> = tip
        .iter()
        .map(|&i| NeumannBc::new(i, [0.0, load / tip.len() as f64, 0.0]))
        .collect();
    let tip_deflection = |analysis: &LinearStaticAnalysis| {
        let result = analysis
            .try_solve(&mesh, &steel(), &clamp, &loads, &zero_force())
            .expect("a 60-DOF cantilever converges");
        tip.iter()
            .map(|&i| result.displacements[i].y)
            .fold(f64::INFINITY, f64::min)
    };

    let by_default = tip_deflection(&LinearStaticAnalysis::new());
    let reference = tip_deflection(&LinearStaticAnalysis {
        max_iter: 200_000,
        tolerance: 1e-14,
    });
    assert!(by_default < 0.0, "tip deflection {by_default}");
    assert!(
        (by_default - reference).abs() <= 1e-4 * reference.abs(),
        "default solve {by_default:e} vs reference {reference:e}"
    );
    // Displacement-based finite elements are never more flexible than the exact
    // solution, which is close to the Euler-Bernoulli value P L^3 / (3 E I).
    let second_moment = width * height.powi(3) / 12.0;
    let euler_bernoulli = load * length.powi(3) / (3.0 * 200.0e9 * second_moment);
    assert!(
        by_default.abs() < 1.5 * euler_bernoulli.abs(),
        "tip deflection {by_default:e} vs Euler-Bernoulli {euler_bernoulli:e}"
    );
}

#[test]
fn clamped_beam_without_any_load_solves_to_the_zero_field() {
    // A zero right-hand side is a valid problem whose solution is u = 0; it must
    // not be mistaken for a stalled iteration.
    let (mesh, clamp, _, _) = clamped_beam();
    let result = analysis()
        .try_solve(&mesh, &steel(), &clamp, &[], &zero_force())
        .expect("zero load on a clamped beam is well posed");
    assert!(result.displacements.iter().all(|d| d.norm() == 0.0));
    assert!(
        result
            .stress
            .iter()
            .chain(result.strain.iter())
            .all(|voigt| voigt.iter().all(|&c| c == 0.0))
    );
}

#[test]
fn convergence_criterion_is_relative_to_the_load_norm() {
    // With an absolute residual test a 1e-12 N load would "converge" after the
    // first iteration with a wildly inaccurate field; the relative test must
    // return the same shape of solution scaled by the load factor.
    let (mesh, clamp, loads, tip) = clamped_beam();
    let reference = analysis()
        .try_solve(&mesh, &steel(), &clamp, &loads, &zero_force())
        .expect("well posed");
    for scale in [1.0e-12, 1.0e9] {
        let scaled_loads: Vec<NeumannBc> = loads
            .iter()
            .map(|bc| NeumannBc::new(bc.node, bc.force.map(|component| component * scale)))
            .collect();
        let scaled = analysis()
            .try_solve(&mesh, &steel(), &clamp, &scaled_loads, &zero_force())
            .expect("a scaled load is still well posed");
        for &node in &tip {
            let expected = reference.displacements[node].y * scale;
            let actual = scaled.displacements[node].y;
            assert!(
                (actual - expected).abs() <= 1e-4 * expected.abs(),
                "scale {scale:e}, node {node}: {actual:e} vs {expected:e}"
            );
        }
    }
}

#[test]
fn uniaxial_patch_test_recovers_the_exact_linear_field_and_uniform_stress() {
    // E = 1000, nu = 0.25  =>  lambda = mu = 400.  Prescribing u = (eps x, 0, 0)
    // on the whole boundary of the unit cube must reproduce it everywhere with
    // sigma_xx = (lambda + 2 mu) eps = 1200 eps and sigma_yy = sigma_zz = lambda eps.
    let eps = 1.0e-3;
    let mesh = TetrahedralMesh::generate_beam(1.0, 1.0, 1.0, 2, 2, 2);
    let on_boundary = |p: &Vec3| {
        [p.x, p.y, p.z]
            .iter()
            .any(|&c| c.abs() < 1e-12 || (c - 1.0).abs() < 1e-12)
    };
    let mut prescribed = Vec::new();
    for (i, node) in mesh.nodes.iter().enumerate() {
        if on_boundary(node) {
            prescribed.push(DirichletBc::new(i, 0, eps * node.x));
            prescribed.push(DirichletBc::new(i, 1, 0.0));
            prescribed.push(DirichletBc::new(i, 2, 0.0));
        }
    }
    assert!(
        prescribed.len() < 3 * mesh.num_nodes(),
        "the patch must contain interior nodes"
    );
    let material = LinearElasticMaterial::new(1000.0, 0.25);
    let solver = LinearStaticAnalysis {
        max_iter: 10_000,
        tolerance: 1e-14,
    };
    let result = solver
        .try_solve(&mesh, &material, &prescribed, &[], &zero_force())
        .expect("patch test is well posed");

    for (i, node) in mesh.nodes.iter().enumerate() {
        let u = result.displacements[i];
        assert!(
            (u.x - eps * node.x).abs() < 1e-12,
            "u_x at node {i}: {}",
            u.x
        );
        assert!(u.y.abs() < 1e-12 && u.z.abs() < 1e-12, "node {i}: {u:?}");
    }
    for (e, (sigma, strain)) in result.stress.iter().zip(&result.strain).enumerate() {
        assert!(
            (strain[0] - eps).abs() < 1e-12,
            "element {e} strain {strain:?}"
        );
        assert!(
            (sigma[0] - 1200.0 * eps).abs() < 1e-9,
            "element {e} stress {sigma:?}"
        );
        assert!(
            (sigma[1] - 400.0 * eps).abs() < 1e-9,
            "element {e} stress {sigma:?}"
        );
        assert!(
            (sigma[2] - 400.0 * eps).abs() < 1e-9,
            "element {e} stress {sigma:?}"
        );
        for shear in &sigma[3..] {
            assert!(shear.abs() < 1e-9, "element {e} shear stress {sigma:?}");
        }
    }
}

#[test]
fn crate_errors_convert_into_fem_errors() {
    assert_eq!(
        FemError::from(Error::SingularSystem),
        FemError::SingularSystem
    );
    match FemError::from(Error::DofOutOfRange { dof: 77, n: 12 }) {
        FemError::InvalidInput(message) => {
            assert!(
                message.contains("77") && message.contains("12"),
                "{message}"
            );
        }
        other => panic!("unexpected conversion: {other:?}"),
    }
    let converted = FemError::from(Error::InvalidInput("bad".to_string()));
    assert_eq!(converted, FemError::InvalidInput("bad".to_string()));
    assert_eq!(converted.to_string(), "invalid input: bad");
    assert_eq!(
        FemError::from(Error::General("oops".to_string())),
        FemError::InvalidInput("oops".to_string())
    );
}
