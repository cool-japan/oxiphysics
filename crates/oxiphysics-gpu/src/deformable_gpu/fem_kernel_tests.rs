// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`FEMGpuKernel`] (Saint Venant-Kirchhoff tetrahedra).
//!
//! The kernel used to return three edge springs `-mu (|e| - |e0|) / |e0| V0 e/|e|`
//! along the edges `a-b`, `a-c`, `a-d` only: `lambda` and the other three edges
//! did not enter, the forces were not the gradient of the reported energy, a
//! non-orthogonal rest tetrahedron was not even force free at rest, and the
//! unchecked `rest_positions[ia]` panicked for a short rest array.  The tests
//! check the element against its defining energy: rest pose and rigid motions
//! are force free, the forces are the central-difference gradient of the
//! energy, closed forms for uniform and uniaxial stretch, and the guards.

use super::functions::{add3, cross3, dot3, length3, normalize3, scale3, sub3};
use super::types::*;

/// Rest tetrahedron with a 45 degree skew between its edges: `Dm` is not
/// orthogonal, so edge springs would not vanish at rest.  `|det Dm| / 6 = 1/6`.
const SKEWED_REST: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Fixed, asymmetric node displacements.
const DISPLACEMENTS: [[f64; 3]; 4] = [
    [0.02, -0.01, 0.03],
    [0.10, 0.05, -0.07],
    [-0.04, 0.12, 0.02],
    [0.06, -0.08, 0.11],
];

const YOUNGS: f64 = 1.0e5;
const POISSON: f64 = 0.3;

fn kernel_with(rest_volume: f64, youngs: f64, poisson: f64) -> FEMGpuKernel {
    let mut kernel = FEMGpuKernel::new(4);
    kernel.add_element(TetElement::new([0, 1, 2, 3], rest_volume, youngs, poisson));
    kernel
}

fn skewed_kernel() -> FEMGpuKernel {
    kernel_with(1.0 / 6.0, YOUNGS, POISSON)
}

fn displaced(rest: &[[f64; 3]], scale: f64) -> Vec<[f64; 3]> {
    rest.iter()
        .zip(DISPLACEMENTS.iter())
        .map(|(p, d)| add3(*p, scale3(*d, scale)))
        .collect()
}

/// Rodrigues rotation of `p` about the axis through the origin.
fn rotate_about(axis: [f64; 3], angle: f64, p: [f64; 3]) -> [f64; 3] {
    let k = normalize3(axis);
    let (s, c) = angle.sin_cos();
    add3(
        add3(scale3(p, c), scale3(cross3(k, p), s)),
        scale3(k, dot3(k, p) * (1.0 - c)),
    )
}

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: {actual} vs {expected} (tol {tol})"
    );
}

#[test]
fn rest_pose_of_a_skewed_tetrahedron_is_force_and_energy_free() {
    let kernel = skewed_kernel();
    let forces = kernel.compute_forces(&SKEWED_REST, &SKEWED_REST);
    assert_eq!(forces.len(), 4);
    for (i, f) in forces.iter().enumerate() {
        assert!(length3(*f) < 1e-8, "node {i} rest force {f:?}");
    }
    assert!(kernel.compute_strain_energy(&SKEWED_REST, &SKEWED_REST) < 1e-12);
}

#[test]
fn rigid_motion_is_force_and_energy_free() {
    let kernel = skewed_kernel();
    let moved: Vec<[f64; 3]> = SKEWED_REST
        .iter()
        .map(|p| add3(rotate_about([1.0, 2.0, 3.0], 0.9, *p), [4.0, -2.0, 7.5]))
        .collect();
    for (i, f) in kernel
        .compute_forces(&moved, &SKEWED_REST)
        .iter()
        .enumerate()
    {
        assert!(length3(*f) < 1e-6, "node {i} force {f:?} in rigid motion");
    }
    assert!(kernel.compute_strain_energy(&moved, &SKEWED_REST) < 1e-9);
}

/// `f = -dE/dx` by central differences, for every node and axis.
#[test]
fn forces_are_the_negative_energy_gradient() {
    let kernel = skewed_kernel();
    let positions = displaced(&SKEWED_REST, 1.0);
    let forces = kernel.compute_forces(&positions, &SKEWED_REST);
    assert!(
        forces.iter().map(|f| length3(*f)).fold(0.0, f64::max) > 1.0,
        "the deformed element must carry force"
    );
    let h = 1e-6;
    for node in 0..4 {
        for axis in 0..3 {
            let mut plus = positions.clone();
            let mut minus = positions.clone();
            plus[node][axis] += h;
            minus[node][axis] -= h;
            let gradient = (kernel.compute_strain_energy(&plus, &SKEWED_REST)
                - kernel.compute_strain_energy(&minus, &SKEWED_REST))
                / (2.0 * h);
            let tol = 1e-5 * (1.0 + forces[node][axis].abs());
            assert_close(
                forces[node][axis],
                -gradient,
                tol,
                &format!("force[{node}][{axis}]"),
            );
        }
    }
}

/// The gradient property also holds for large strain and for a different
/// material, so it is not tied to the linear regime.
#[test]
fn forces_are_the_negative_energy_gradient_at_large_strain() {
    let kernel = kernel_with(1.0 / 6.0, 3.0e6, 0.42);
    let positions = displaced(&SKEWED_REST, 4.0);
    let forces = kernel.compute_forces(&positions, &SKEWED_REST);
    let h = 1e-6;
    for node in 0..4 {
        for axis in 0..3 {
            let mut plus = positions.clone();
            let mut minus = positions.clone();
            plus[node][axis] += h;
            minus[node][axis] -= h;
            let gradient = (kernel.compute_strain_energy(&plus, &SKEWED_REST)
                - kernel.compute_strain_energy(&minus, &SKEWED_REST))
                / (2.0 * h);
            let tol = 1e-5 * (1.0 + forces[node][axis].abs());
            assert_close(
                forces[node][axis],
                -gradient,
                tol,
                &format!("force[{node}][{axis}]"),
            );
        }
    }
}

#[test]
fn forces_carry_no_net_force_or_torque() {
    let kernel = skewed_kernel();
    let positions = displaced(&SKEWED_REST, 1.0);
    let forces = kernel.compute_forces(&positions, &SKEWED_REST);
    let mut net = [0.0; 3];
    let mut torque = [0.0; 3];
    for (p, f) in positions.iter().zip(forces.iter()) {
        net = add3(net, *f);
        torque = add3(torque, cross3(*p, *f));
    }
    let scale = forces.iter().map(|f| length3(*f)).fold(0.0, f64::max);
    assert!(scale > 1.0);
    assert!(length3(net) < 1e-9 * scale, "net force {net:?}");
    assert!(length3(torque) < 1e-8 * scale, "net torque {torque:?}");
}

/// `F = (1 + eps) I`: `E = (eps + eps^2/2) I`, so the energy is
/// `V0 e^2 (3 mu + 9/2 lambda)` with `e = eps + eps^2/2` -- both Lamé
/// parameters enter, and the forces restore the shape.
#[test]
fn uniform_stretch_energy_matches_the_closed_form() {
    let element = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, 2.0e5, 0.25);
    let (lambda, mu) = element.lame_params();
    let mut kernel = FEMGpuKernel::new(4);
    kernel.add_element(element);
    for eps in [0.01, -0.02, 0.1] {
        let stretched: Vec<[f64; 3]> = SKEWED_REST.iter().map(|p| scale3(*p, 1.0 + eps)).collect();
        let e = eps + 0.5 * eps * eps;
        let expected = (1.0 / 6.0) * e * e * (3.0 * mu + 4.5 * lambda);
        assert_close(
            kernel.compute_strain_energy(&stretched, &SKEWED_REST),
            expected,
            1e-9 * expected,
            &format!("uniform stretch eps = {eps}"),
        );
        // Restoring: the forces point against the displacement from the centroid.
        let forces = kernel.compute_forces(&stretched, &SKEWED_REST);
        let centroid = scale3(
            stretched.iter().fold([0.0; 3], |acc, p| add3(acc, *p)),
            0.25,
        );
        let outward: f64 = stretched
            .iter()
            .zip(forces.iter())
            .map(|(p, f)| dot3(sub3(*p, centroid), *f))
            .sum();
        assert!(
            outward * eps < 0.0,
            "forces must oppose the stretch eps = {eps}: {outward}"
        );
    }
}

/// Uniaxial stretch of the x coordinates, `F = diag(1 + eps, 1, 1)`, is
/// independent of the (non-orthogonal) rest shape: `E = diag(e, 0, 0)` and the
/// energy is `V0 (mu + lambda/2) e^2`.
#[test]
fn uniaxial_stretch_energy_matches_the_closed_form() {
    let element = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, 4.0e5, 0.33);
    let (lambda, mu) = element.lame_params();
    let mut kernel = FEMGpuKernel::new(4);
    kernel.add_element(element);
    for eps in [0.005, 0.05, -0.05] {
        let stretched: Vec<[f64; 3]> = SKEWED_REST
            .iter()
            .map(|p| [p[0] * (1.0 + eps), p[1], p[2]])
            .collect();
        let e = eps + 0.5 * eps * eps;
        let expected = (1.0 / 6.0) * (mu + 0.5 * lambda) * e * e;
        assert_close(
            kernel.compute_strain_energy(&stretched, &SKEWED_REST),
            expected,
            1e-9 * expected,
            &format!("uniaxial stretch eps = {eps}"),
        );
    }
}

/// Stretching the edge `a-b` along x gives restoring forces on its end nodes,
/// and compressing it the opposite sign.
#[test]
fn stretched_and_compressed_edges_are_restored() {
    let kernel = skewed_kernel();
    for (stretch, sign) in [(0.2, -1.0), (-0.2, 1.0)] {
        let mut deformed = SKEWED_REST.to_vec();
        deformed[1][0] += stretch;
        let forces = kernel.compute_forces(&deformed, &SKEWED_REST);
        assert!(sign * forces[1][0] > 0.0, "f_b = {:?}", forces[1]);
        assert!(sign * forces[0][0] < 0.0, "f_a = {:?}", forces[0]);
    }
}

/// At small strain Saint Venant-Kirchhoff and the corotational model are both
/// linear elasticity, so their energies agree to first order in the strain.
#[test]
fn small_strain_energy_agrees_with_the_corotational_element() {
    let kernel = skewed_kernel();
    let corotational = CorotationalElement::from_rest([0, 1, 2, 3], &SKEWED_REST, YOUNGS, POISSON);
    let positions = displaced(&SKEWED_REST, 1e-4);
    let stvk = kernel.compute_strain_energy(&positions, &SKEWED_REST);
    let corot = corotational.strain_energy(&positions);
    assert!(stvk > 0.0 && corot > 0.0);
    assert!(
        (stvk - corot).abs() < 2e-3 * corot,
        "StVK {stvk} vs corotational {corot}"
    );
    // The forces agree as well.
    let f_stvk = kernel.compute_forces(&positions, &SKEWED_REST);
    let f_corot = corotational.compute_nodal_forces(&positions);
    let scale = f_corot.iter().map(|f| length3(*f)).fold(0.0, f64::max);
    for (node, (a, b)) in f_stvk.iter().zip(f_corot.iter()).enumerate() {
        assert!(
            length3(sub3(*a, *b)) < 2e-3 * scale,
            "node {node}: {a:?} vs {b:?}"
        );
    }
}

#[test]
fn elements_sharing_nodes_assemble_additively() {
    // Two tetrahedra glued along the face (1, 2, 3).
    let apex = [1.0, 1.0, 1.0];
    let rest = vec![
        SKEWED_REST[0],
        SKEWED_REST[1],
        SKEWED_REST[2],
        SKEWED_REST[3],
        apex,
    ];
    // Volume of the second tetrahedron (1, 2, 3, 4): |det| / 6.
    let edges = [
        sub3(rest[2], rest[1]),
        sub3(rest[3], rest[1]),
        sub3(rest[4], rest[1]),
    ];
    let volume_b = dot3(edges[0], cross3(edges[1], edges[2])).abs() / 6.0;
    let first = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, 1.0e5, 0.3);
    let second = TetElement::new([1, 2, 3, 4], volume_b, 2.0e5, 0.2);

    let mut both = FEMGpuKernel::new(0);
    both.add_element(first.clone());
    both.add_element(second.clone());
    assert_eq!(both.num_nodes, 5);
    let mut only_first = FEMGpuKernel::new(5);
    only_first.add_element(first);
    let mut only_second = FEMGpuKernel::new(5);
    only_second.add_element(second);

    let mut positions = displaced(&rest[..4], 1.0);
    positions.push(add3(apex, [0.05, -0.03, 0.07]));
    let total = both.compute_forces(&positions, &rest);
    let a = only_first.compute_forces(&positions, &rest);
    let b = only_second.compute_forces(&positions, &rest);
    for node in 0..5 {
        for axis in 0..3 {
            assert_close(
                total[node][axis],
                a[node][axis] + b[node][axis],
                1e-9 * (1.0 + total[node][axis].abs()),
                "assembled force",
            );
        }
    }
    // Node 4 only belongs to the second element, node 0 only to the first.
    assert_eq!(a[4], [0.0; 3]);
    assert_eq!(b[0], [0.0; 3]);
    let energy = both.compute_strain_energy(&positions, &rest);
    let sum = only_first.compute_strain_energy(&positions, &rest)
        + only_second.compute_strain_energy(&positions, &rest);
    assert!(energy > 0.0);
    assert_close(energy, sum, 1e-12 * energy, "total energy");
}

/// A node index beyond `rest_positions` used to panic (`rest_positions[ia]`);
/// every kind of out-of-range element is now skipped for forces and energy.
#[test]
fn out_of_range_elements_are_skipped_instead_of_panicking() {
    let positions = displaced(&SKEWED_REST, 1.0);

    // Node 3 exists in `positions` but not in a three-entry rest array.
    let kernel = skewed_kernel();
    let short_rest = &SKEWED_REST[..3];
    let forces = kernel.compute_forces(&positions, short_rest);
    assert!(forces.iter().all(|f| *f == [0.0; 3]));
    assert_eq!(kernel.compute_strain_energy(&positions, short_rest), 0.0);

    // Node 3 missing from the current positions.
    let forces = kernel.compute_forces(&positions[..3], &SKEWED_REST);
    assert!(forces.iter().all(|f| *f == [0.0; 3]));
    assert_eq!(
        kernel.compute_strain_energy(&positions[..3], &SKEWED_REST),
        0.0
    );

    // An element pushed past `num_nodes` without `add_element`.
    let mut manual = FEMGpuKernel::new(3);
    manual
        .elements
        .push(TetElement::new([0, 1, 2, 3], 1.0 / 6.0, YOUNGS, POISSON));
    let forces = manual.compute_forces(&positions, &SKEWED_REST);
    assert_eq!(forces.len(), 3);
    assert!(forces.iter().all(|f| *f == [0.0; 3]));
    assert_eq!(manual.compute_strain_energy(&positions, &SKEWED_REST), 0.0);

    // A valid element still works next to an invalid one.
    let mut mixed = skewed_kernel();
    mixed
        .elements
        .push(TetElement::new([0, 1, 2, 99], 1.0 / 6.0, YOUNGS, POISSON));
    let valid_only = skewed_kernel();
    assert_eq!(
        mixed.compute_forces(&positions, &SKEWED_REST),
        valid_only.compute_forces(&positions, &SKEWED_REST)
    );
}

#[test]
fn degenerate_or_invalid_elements_carry_no_force() {
    let positions = displaced(&SKEWED_REST, 1.0);
    let flat_rest = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let flat_positions = displaced(&flat_rest, 1.0);
    let flat = skewed_kernel();
    assert_eq!(flat.compute_strain_energy(&flat_positions, &flat_rest), 0.0);
    let forces = flat.compute_forces(&flat_positions, &flat_rest);
    assert!(forces.iter().all(|f| *f == [0.0; 3]));

    // A repeated node is degenerate too.
    let mut repeated = FEMGpuKernel::new(4);
    repeated.add_element(TetElement::new([0, 1, 1, 3], 1.0 / 6.0, YOUNGS, POISSON));
    assert_eq!(
        repeated.compute_strain_energy(&positions, &SKEWED_REST),
        0.0
    );

    // Non-positive and non-finite rest volumes.
    for volume in [0.0, -1.0 / 6.0, f64::NAN, f64::INFINITY] {
        let kernel = kernel_with(volume, YOUNGS, POISSON);
        assert_eq!(
            kernel.compute_strain_energy(&positions, &SKEWED_REST),
            0.0,
            "rest volume {volume}"
        );
        let forces = kernel.compute_forces(&positions, &SKEWED_REST);
        assert!(
            forces.iter().all(|f| *f == [0.0; 3]),
            "rest volume {volume}"
        );
    }
}

#[test]
fn incompressible_poisson_ratio_stays_finite() {
    let element = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, YOUNGS, 0.5);
    let (lambda, mu) = element.lame_params();
    assert!(lambda.is_finite() && mu.is_finite() && lambda > 0.0);
    let mut kernel = FEMGpuKernel::new(4);
    kernel.add_element(element);
    let positions = displaced(&SKEWED_REST, 1.0);
    assert!(
        kernel
            .compute_strain_energy(&positions, &SKEWED_REST)
            .is_finite()
    );
    for f in kernel.compute_forces(&positions, &SKEWED_REST) {
        assert!(f.iter().all(|c| c.is_finite()));
    }
    for f in kernel.compute_forces(&SKEWED_REST, &SKEWED_REST) {
        assert!(f.iter().all(|c| c.is_finite() && c.abs() < 1e-6));
    }
}
