// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the bending stencils of [`XpbdClothMesh`].
//!
//! The old stencils used a wing vertex that lay on the hinge line itself, so
//! every bending constraint saw a zero face normal and did nothing.  The
//! stencils are now derived from the triangles that share an edge.
//!
//! Once the constraints act, the cloth also depends on the hinge projection of
//! [`BendingConstraint`]: the gravity tests below check that bending stiffens
//! the sheet and does so equally for both bending directions.  The hinge-level
//! behaviour (gradients, conservation laws, mirror symmetry) is covered in
//! `constraint::bending_tests`.

use super::types::XpbdClothMesh;
use crate::constraint::{BendingConstraint, SoftConstraint};
use crate::particle::SoftParticle;
use crate::solver::XpbdSolver;
use oxiphysics_core::math::Vec3;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// Undirected edge with sorted endpoints.
fn sorted_edge(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Interior edges (shared by exactly two triangles) of the cloth triangulation.
fn interior_edges(cloth: &XpbdClothMesh) -> BTreeSet<(usize, usize)> {
    let mut counts: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for &[a, b, c] in &cloth.body.triangles {
        for (u, v) in [(a, b), (b, c), (c, a)] {
            *counts.entry(sorted_edge(u, v)).or_insert(0) += 1;
        }
    }
    counts
        .into_iter()
        .filter(|&(_, count)| count == 2)
        .map(|(edge, _)| edge)
        .collect()
}

/// Face normal `(p_wing - p_a) x (p_b - p_a)` of the triangle `(a, b, wing)`.
fn wing_normal(particles: &[SoftParticle], a: usize, b: usize, wing: usize) -> Vec3 {
    let pa = particles[a].position;
    (particles[wing].position - pa).cross(&(particles[b].position - pa))
}

#[test]
fn bending_constraints_cover_every_interior_edge_exactly_once() {
    let (nx, ny) = (5, 4);
    let cloth = XpbdClothMesh::new(nx, ny, 4.0, 3.0, 1.0, 1e-4);
    let expected = interior_edges(&cloth);
    // Horizontal, vertical and diagonal interior edges of the triangulated grid.
    let analytic = (nx - 1) * (ny - 2) + (nx - 2) * (ny - 1) + (nx - 1) * (ny - 1);
    assert_eq!(expected.len(), analytic);
    assert_eq!(cloth.bending_constraints.len(), analytic);
    let hinges: BTreeSet<(usize, usize)> = cloth
        .bending_constraints
        .iter()
        .map(|c| sorted_edge(c.indices[0], c.indices[1]))
        .collect();
    assert_eq!(
        hinges, expected,
        "one hinge per interior edge, no duplicates"
    );
}

#[test]
fn every_stencil_is_a_two_triangle_hinge_with_non_degenerate_wings() {
    let cloth = XpbdClothMesh::new(6, 5, 5.0, 4.0, 1.0, 1e-4);
    assert!(!cloth.bending_constraints.is_empty());
    for constraint in &cloth.bending_constraints {
        let [a, b, c, d] = constraint.indices;
        let mut distinct = [a, b, c, d];
        distinct.sort_unstable();
        assert!(
            distinct.windows(2).all(|w| w[0] != w[1]),
            "stencil {:?} repeats a particle",
            constraint.indices
        );
        // The wings are the apexes of the two triangles that share the hinge.
        let apexes: BTreeSet<usize> = cloth
            .body
            .triangles
            .iter()
            .filter(|tri| tri.contains(&a) && tri.contains(&b))
            .flat_map(|tri| tri.iter().copied().filter(|&v| v != a && v != b))
            .collect();
        assert_eq!(
            apexes,
            BTreeSet::from([c, d]),
            "stencil {:?} is not made of the two triangles sharing its hinge",
            constraint.indices
        );
        // Neither wing may lie on the hinge line (the old, inert stencils did).
        let normal_c = wing_normal(&cloth.body.particles, a, b, c);
        let normal_d = wing_normal(&cloth.body.particles, a, b, d);
        assert!(
            normal_c.norm() > 1e-6,
            "wing {c} of {:?} is on the hinge",
            constraint.indices
        );
        assert!(
            normal_d.norm() > 1e-6,
            "wing {d} of {:?} is on the hinge",
            constraint.indices
        );
    }
}

#[test]
fn flat_cloth_has_rest_angle_pi_and_antiparallel_wing_normals() {
    let cloth = XpbdClothMesh::new(7, 6, 6.0, 5.0, 1.0, 1e-4);
    for constraint in &cloth.bending_constraints {
        assert!(
            (constraint.rest_angle - PI).abs() < 1e-6,
            "flat hinge {:?} must have rest angle pi, got {}",
            constraint.indices,
            constraint.rest_angle
        );
        let [a, b, c, d] = constraint.indices;
        let n_c = wing_normal(&cloth.body.particles, a, b, c);
        let n_d = wing_normal(&cloth.body.particles, a, b, d);
        // The wings sit on opposite sides of the hinge, so the face normals built
        // with the same edge direction are antiparallel.
        let cos = n_c.dot(&n_d) / (n_c.norm() * n_d.norm());
        assert!(
            (cos + 1.0).abs() < 1e-9,
            "hinge {:?}: cos = {cos}",
            constraint.indices
        );
    }
}

#[test]
fn bending_constraint_count_scales_with_grid_size() {
    for (nx, ny) in [(2, 2), (3, 2), (2, 3), (3, 3), (10, 7)] {
        let cloth = XpbdClothMesh::new(nx, ny, 1.0, 1.0, 1.0, 1e-4);
        let analytic =
            (nx - 1) * ny.saturating_sub(2) + nx.saturating_sub(2) * (ny - 1) + (nx - 1) * (ny - 1);
        assert_eq!(
            cloth.bending_constraints.len(),
            analytic,
            "grid {nx} x {ny}"
        );
        assert_eq!(
            cloth.bending_constraints.len(),
            interior_edges(&cloth).len()
        );
    }
}

/// The inert stencils the old `XpbdClothMesh::new` generated: the wing
/// `(i + 1, j)` lies in the middle of the hinge `(i, j) -> (i + 2, j)`, and the
/// `j`-direction twin has the same defect.
fn legacy_collinear_stencils(cloth: &XpbdClothMesh) -> Vec<BendingConstraint> {
    let (nx, ny) = (cloth.nx, cloth.ny);
    let idx = |ix: usize, iy: usize| iy * nx + ix;
    let mut legacy = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if i + 2 < nx && j + 1 < ny {
                legacy.push(BendingConstraint::from_particles(
                    [idx(i, j), idx(i + 2, j), idx(i + 1, j), idx(i + 1, j + 1)],
                    &cloth.body.particles,
                    1e-6,
                ));
            }
            if j + 2 < ny && i + 1 < nx {
                legacy.push(BendingConstraint::from_particles(
                    [idx(i, j), idx(i, j + 2), idx(i, j + 1), idx(i + 1, j + 1)],
                    &cloth.body.particles,
                    1e-6,
                ));
            }
        }
    }
    legacy
}

/// Hang a 9 x 9 cloth pinned along the `iy = 0` edge under a vertical
/// acceleration `gravity_y` with the given bending constraints and return the
/// mean displacement of the free particles *along gravity* once it has settled
/// (positive means the cloth sags in the direction of gravity).
fn settled_mean_sag(
    bending: impl FnOnce(&XpbdClothMesh) -> Vec<BendingConstraint>,
    gravity_y: f64,
) -> f64 {
    let n = 9;
    let mut cloth = XpbdClothMesh::new(n, n, 2.0, 2.0, 1.0, 1e-7);
    for ix in 0..n {
        cloth.pin(ix, 0);
    }
    cloth.bending_constraints = bending(&cloth);
    cloth.body.damping = 0.05;
    let mut constraints: Vec<Box<dyn SoftConstraint>> = cloth
        .distance_constraints
        .iter()
        .cloned()
        .map(|c| Box::new(c) as Box<dyn SoftConstraint>)
        .chain(
            cloth
                .bending_constraints
                .iter()
                .cloned()
                .map(|c| Box::new(c) as Box<dyn SoftConstraint>),
        )
        .collect();
    let gravity = Vec3::new(0.0, gravity_y, 0.0);
    let mut solver = XpbdSolver::new(10);
    for _ in 0..300 {
        cloth.body.apply_force(&gravity);
        solver.solve(&mut cloth.body, &mut constraints, 1.0 / 60.0);
        cloth.body.clear_forces();
    }
    let free: Vec<f64> = cloth
        .body
        .particles
        .iter()
        .filter(|p| !p.is_static())
        .map(|p| p.position.y)
        .collect();
    free.iter().sum::<f64>() / free.len() as f64 * gravity_y.signum()
}

#[test]
fn cloth_pinned_along_one_edge_sags_less_with_bending_enabled() {
    let with_bending = settled_mean_sag(|cloth| cloth.bending_constraints.clone(), -9.81);
    let without_bending = settled_mean_sag(|_| Vec::new(), -9.81);
    // Without any bending resistance the cloth hangs like a curtain.
    assert!(
        without_bending > 0.5,
        "an unsupported cloth must hang, mean sag {without_bending}"
    );
    // The hinge stencils stiffen the cloth.
    assert!(
        with_bending < without_bending,
        "bending must reduce the sag: {with_bending} (bending) vs {without_bending} (none)"
    );
    assert!(
        with_bending > 0.0,
        "gravity still bends the cloth down: {with_bending}"
    );
}

#[test]
fn bending_stiffness_does_not_depend_on_the_bending_direction() {
    // Reversing gravity bends the cloth the other way.  The sheet and its
    // triangulation are mirror symmetric about the plane of the cloth, so a
    // correct hinge constraint gives the mirror image of the same sag.  The old
    // projection pushed along a fixed normal and would stiffen only one sense.
    let down = settled_mean_sag(|cloth| cloth.bending_constraints.clone(), -9.81);
    let up = settled_mean_sag(|cloth| cloth.bending_constraints.clone(), 9.81);
    assert!(
        (down - up).abs() < 1e-9 * (1.0 + down.abs()),
        "sag with gravity down {down} vs gravity up {up}"
    );
    let free_down = settled_mean_sag(|_| Vec::new(), -9.81);
    let free_up = settled_mean_sag(|_| Vec::new(), 9.81);
    assert!((free_down - free_up).abs() < 1e-9 * (1.0 + free_down.abs()));
    assert!(down < free_down && up < free_up);
}

#[test]
fn legacy_collinear_stencils_were_degenerate_at_rest() {
    // Documents the defect: the middle vertex of the old hinge lies exactly on
    // the hinge line, so the first face normal vanishes, `compute_dihedral`
    // falls back to 0 and the constraint has no restoring effect on the flat
    // cloth (rest angle 0, zero violation).
    let cloth = XpbdClothMesh::new(6, 5, 5.0, 4.0, 1.0, 1e-4);
    let legacy = legacy_collinear_stencils(&cloth);
    assert!(!legacy.is_empty());
    for constraint in &legacy {
        let [a, b, c, _] = constraint.indices;
        let normal = wing_normal(&cloth.body.particles, a, b, c);
        assert!(
            normal.norm() < 1e-12,
            "wing {c} of {:?} is not on the hinge",
            constraint.indices
        );
        assert_eq!(constraint.rest_angle, 0.0, "{:?}", constraint.indices);
    }
    // The new stencils are the exact opposite: proper hinges with rest angle pi.
    for constraint in &cloth.bending_constraints {
        assert!((constraint.rest_angle - PI).abs() < 1e-6);
    }
    // A flat cloth stays exactly flat under the new constraints.
    let mut particles = cloth.body.particles.clone();
    let before: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
    for constraint in &cloth.bending_constraints {
        let mut constraint = constraint.clone();
        constraint.project(&mut particles, 1.0 / 60.0);
    }
    for (p, b) in particles.iter().zip(&before) {
        assert!(
            (p.position - b).norm() < 1e-12,
            "a flat cloth must not be bent"
        );
    }
}
