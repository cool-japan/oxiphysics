// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the tetrahedron volume gradients and for the volume
//! constraints that are built on them.
//!
//! The signed volume `V = e1 · (e2 × e3) / 6` is multilinear in the vertex
//! coordinates, so central differences of it are exact up to rounding; every
//! gradient is therefore checked against finite differences for all four
//! vertices and all three axes, and the four gradients of a tetrahedron must
//! sum to zero (the volume is translation invariant).

use super::functions::*;
use super::types::{TetrahedralSoftBody, VolumeParticle};

const FD_STEP: f64 = 1.0e-6;

type Tet = [[f64; 3]; 4];

/// Irregular, well-shaped tetrahedron with positive orientation.
fn skew_tet() -> Tet {
    [
        [0.10, -0.20, 0.30],
        [1.30, 0.10, -0.15],
        [-0.25, 0.95, 0.20],
        [0.20, 0.15, 1.10],
    ]
}

/// The same tetrahedron with two vertices swapped: negative orientation.
fn mirrored(tet: Tet) -> Tet {
    [tet[0], tet[2], tet[1], tet[3]]
}

/// Right-angled unit tetrahedron with positive orientation (volume 1/6).
fn unit_tet() -> Tet {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

fn signed_volume(tet: &Tet) -> f64 {
    tet_volume_signed(tet[0], tet[1], tet[2], tet[3])
}

fn centroid(tet: &Tet) -> [f64; 3] {
    tet_centroid(tet)
}

/// Scale `tet` about its centroid so that its volume becomes `volume_ratio`
/// times the original one.
fn compressed_to(tet: Tet, volume_ratio: f64) -> Tet {
    let c = centroid(&tet);
    let factor = volume_ratio.cbrt();
    tet.map(|p| add3(c, scale3(sub3(p, c), factor)))
}

fn particle(position: [f64; 3], mass: f64) -> VolumeParticle {
    VolumeParticle {
        position,
        prev_position: position,
        mass,
        fixed: false,
    }
}

/// Body whose single tetrahedron `index_order` was built at `rest` (so that its
/// rest volume and orientation are recorded) and is then moved to `current`.
fn body_at(
    rest: Tet,
    current: Tet,
    index_order: [usize; 4],
    masses: [f64; 4],
) -> TetrahedralSoftBody {
    let particles = rest
        .iter()
        .zip(masses)
        .map(|(&position, mass)| particle(position, mass))
        .collect();
    let mut body = TetrahedralSoftBody::new(particles, vec![index_order]);
    for (p, &position) in body.particles.iter_mut().zip(current.iter()) {
        p.position = position;
        p.prev_position = position;
    }
    body
}

fn body_signed_volume(body: &TetrahedralSoftBody) -> f64 {
    let t = body.tetrahedra[0];
    tet_volume_signed(
        body.particles[t[0]].position,
        body.particles[t[1]].position,
        body.particles[t[2]].position,
        body.particles[t[3]].position,
    )
}

fn body_volume_error(body: &TetrahedralSoftBody) -> f64 {
    (body_signed_volume(body).abs() - body.rest_volumes[0]).abs()
}

fn mass_weighted_centre(body: &TetrahedralSoftBody) -> [f64; 3] {
    body.particles
        .iter()
        .fold([0.0; 3], |acc, p| add3(acc, scale3(p.position, p.mass)))
}

// ---------------------------------------------------------------------------
// tet_volume_gradient
// ---------------------------------------------------------------------------

#[test]
fn tet_volume_gradient_unit_tet_has_the_analytic_vertex_gradients() {
    let t = unit_tet();
    let g = tet_volume_gradient(t[0], t[1], t[2], t[3]);
    let sixth = 1.0 / 6.0;
    let expected = [
        [-sixth, -sixth, -sixth],
        [sixth, 0.0, 0.0],
        [0.0, sixth, 0.0],
        [0.0, 0.0, sixth],
    ];
    for (vertex, (got, want)) in g.iter().zip(expected.iter()).enumerate() {
        for (axis, (&a, &b)) in got.iter().zip(want.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1.0e-15,
                "dV/dp{vertex}[{axis}] = {a}, expected {b}"
            );
        }
    }
}

#[test]
fn tet_volume_gradient_matches_central_differences_for_every_vertex_and_axis() {
    for tet in [skew_tet(), mirrored(skew_tet())] {
        let grads = tet_volume_gradient(tet[0], tet[1], tet[2], tet[3]);
        for (vertex, gradient) in grads.iter().enumerate() {
            for (axis, &analytic) in gradient.iter().enumerate() {
                let mut plus = tet;
                plus[vertex][axis] += FD_STEP;
                let mut minus = tet;
                minus[vertex][axis] -= FD_STEP;
                let numeric = (signed_volume(&plus) - signed_volume(&minus)) / (2.0 * FD_STEP);
                assert!(
                    (numeric - analytic).abs() < 1.0e-6,
                    "dV/dp{vertex}[{axis}]: finite difference {numeric}, analytic {analytic}"
                );
            }
        }
    }
}

#[test]
fn tet_volume_gradients_sum_to_zero() {
    for tet in [unit_tet(), skew_tet(), mirrored(skew_tet())] {
        let grads = tet_volume_gradient(tet[0], tet[1], tet[2], tet[3]);
        let sum = grads.iter().fold([0.0; 3], |acc, &g| add3(acc, g));
        for (axis, component) in sum.iter().enumerate() {
            assert!(
                component.abs() < 1.0e-15,
                "gradients of a translation-invariant volume must sum to zero, axis {axis}: {component}"
            );
        }
    }
}

#[test]
fn tet_volume_gradient_of_the_first_vertex_is_the_opposite_face_normal() {
    // dV/dp0 = (p3 - p1) x (p2 - p1) / 6.
    let t = skew_tet();
    let g = tet_volume_gradient(t[0], t[1], t[2], t[3]);
    let direct = scale3(cross3(sub3(t[3], t[1]), sub3(t[2], t[1])), 1.0 / 6.0);
    for (axis, (&a, &b)) in g[0].iter().zip(direct.iter()).enumerate() {
        assert!((a - b).abs() < 1.0e-14, "axis {axis}: {a} vs {b}");
    }
}

// ---------------------------------------------------------------------------
// volume_preservation_gradient
// ---------------------------------------------------------------------------

/// Two tetrahedra of opposite orientation that share the face `(1, 2, 3)`.
fn two_tet_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let t = skew_tet();
    let apex = [1.0, 0.9, 0.3];
    (
        vec![t[0], t[1], t[2], t[3], apex],
        // Tet 0 has positive orientation; tet 1 reuses the shared face with
        // the apex on the other side, which gives it the opposite orientation.
        vec![[0, 1, 2, 3], [4, 1, 2, 3]],
    )
}

#[test]
fn volume_preservation_gradient_matches_central_differences_of_the_energy() {
    let (positions, tetrahedra) = two_tet_mesh();
    let orientation = |tet: &[usize; 4]| {
        tet_volume_signed(
            positions[tet[0]],
            positions[tet[1]],
            positions[tet[2]],
            positions[tet[3]],
        )
    };
    assert!(orientation(&tetrahedra[0]) > 0.1);
    assert!(
        orientation(&tetrahedra[1]) < -0.1,
        "the second tetrahedron has the opposite orientation"
    );
    let rest_volumes = [0.11, 0.07];
    let stiffness = 250.0;
    let gradient = volume_preservation_gradient(&tetrahedra, &rest_volumes, &positions, stiffness);
    assert_eq!(gradient.len(), positions.len());
    for (vertex, vertex_gradient) in gradient.iter().enumerate() {
        for (axis, &analytic) in vertex_gradient.iter().enumerate() {
            let mut plus = positions.clone();
            plus[vertex][axis] += FD_STEP;
            let mut minus = positions.clone();
            minus[vertex][axis] -= FD_STEP;
            let numeric =
                (volume_preservation_energy(&tetrahedra, &rest_volumes, &plus, stiffness)
                    - volume_preservation_energy(&tetrahedra, &rest_volumes, &minus, stiffness))
                    / (2.0 * FD_STEP);
            assert!(
                (numeric - analytic).abs() < 1.0e-5 * (1.0 + analytic.abs()),
                "dE/dx[{vertex}][{axis}]: finite difference {numeric}, analytic {analytic}"
            );
        }
    }
}

#[test]
fn volume_preservation_gradient_matches_finite_differences_for_both_orientations() {
    for tet in [skew_tet(), mirrored(skew_tet())] {
        let positions = tet.to_vec();
        let tetrahedra = [[0usize, 1, 2, 3]];
        let rest_volumes = [0.2];
        let stiffness = 80.0;
        let gradient =
            volume_preservation_gradient(&tetrahedra, &rest_volumes, &positions, stiffness);
        for (vertex, vertex_gradient) in gradient.iter().enumerate() {
            for (axis, &analytic) in vertex_gradient.iter().enumerate() {
                let mut plus = positions.clone();
                plus[vertex][axis] += FD_STEP;
                let mut minus = positions.clone();
                minus[vertex][axis] -= FD_STEP;
                let numeric =
                    (volume_preservation_energy(&tetrahedra, &rest_volumes, &plus, stiffness)
                        - volume_preservation_energy(
                            &tetrahedra,
                            &rest_volumes,
                            &minus,
                            stiffness,
                        ))
                        / (2.0 * FD_STEP);
                assert!(
                    (numeric - analytic).abs() < 1.0e-5 * (1.0 + analytic.abs()),
                    "dE/dx[{vertex}][{axis}]: finite difference {numeric}, analytic {analytic}"
                );
            }
        }
    }
}

#[test]
fn volume_preservation_gradient_of_one_tet_sums_to_zero() {
    for tet in [skew_tet(), mirrored(skew_tet())] {
        let positions = tet.to_vec();
        let gradient = volume_preservation_gradient(&[[0, 1, 2, 3]], &[0.2], &positions, 120.0);
        let sum = gradient.iter().fold([0.0; 3], |acc, &g| add3(acc, g));
        for (axis, component) in sum.iter().enumerate() {
            assert!(
                component.abs() < 1.0e-12,
                "a translation of the whole tetrahedron does not change its energy, axis {axis}: {component}"
            );
        }
        let magnitude: f64 = gradient.iter().map(|&g| len3(g)).sum();
        assert!(
            magnitude > 1.0e-3,
            "the volume error must produce a non-trivial gradient, got {magnitude}"
        );
    }
}

#[test]
fn apply_volume_preservation_descends_the_penalty_energy() {
    let mut body = body_at(unit_tet(), unit_tet(), [0, 1, 2, 3], [1.0; 4]);
    body.particles[3].position = [0.1, 0.05, 2.0];
    let rest_volumes = body.rest_volumes.clone();
    let stiffness = 1000.0;
    let energy = |body: &TetrahedralSoftBody| {
        let positions: Vec<[f64; 3]> = body.particles.iter().map(|p| p.position).collect();
        volume_preservation_energy(&body.tetrahedra, &rest_volumes, &positions, stiffness)
    };
    let mut previous = energy(&body);
    for step in 0..20 {
        body.apply_volume_preservation(&rest_volumes, stiffness);
        let now = energy(&body);
        assert!(
            now < previous,
            "step {step}: the penalty energy must decrease, {previous} -> {now}"
        );
        previous = now;
    }
}

// ---------------------------------------------------------------------------
// TetrahedralSoftBody::solve_volume_constraints
// ---------------------------------------------------------------------------

#[test]
fn new_records_the_rest_orientation_of_every_tetrahedron() {
    let t = skew_tet();
    let particles: Vec<VolumeParticle> = t.iter().map(|&p| particle(p, 1.0)).collect();
    let body = TetrahedralSoftBody::new(particles, vec![[0, 1, 2, 3], [0, 2, 1, 3]]);
    assert_eq!(body.rest_orientations, [1.0, -1.0]);
    assert!((body.rest_volumes[0] - body.rest_volumes[1]).abs() < 1.0e-15);
    assert!(body.rest_volumes.iter().all(|&v| v > 0.0));
}

#[test]
fn volume_projection_of_a_compressed_equal_mass_tet_keeps_the_centroid_and_converges_monotonically()
{
    let rest = skew_tet();
    let mut body = body_at(rest, compressed_to(rest, 0.8), [0, 1, 2, 3], [1.0; 4]);
    let centre_before = centroid(&[
        body.particles[0].position,
        body.particles[1].position,
        body.particles[2].position,
        body.particles[3].position,
    ]);
    let v_rest = body.rest_volumes[0];
    let mut previous = body_volume_error(&body);
    assert!(
        (previous / v_rest - 0.2).abs() < 1.0e-12,
        "the tetrahedron starts 20% compressed"
    );
    for iteration in 0..10 {
        body.solve_volume_constraints();
        let error = body_volume_error(&body);
        assert!(
            error < previous,
            "iteration {iteration}: |V - V0| must decrease, {previous} -> {error}"
        );
        previous = error;
        let centre = centroid(&[
            body.particles[0].position,
            body.particles[1].position,
            body.particles[2].position,
            body.particles[3].position,
        ]);
        for (axis, (&a, &b)) in centre.iter().zip(centre_before.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1.0e-12,
                "iteration {iteration}: the projection injected momentum along axis {axis}: {a} vs {b}"
            );
        }
    }
    assert!(
        previous < 0.002 * v_rest,
        "ten half-stiffness projections remove almost all of a 20% error, left {}",
        previous / v_rest
    );
}

#[test]
fn volume_projection_conserves_linear_momentum_for_unequal_masses() {
    let rest = skew_tet();
    let mut body = body_at(
        rest,
        compressed_to(rest, 0.7),
        [0, 1, 2, 3],
        [1.0, 2.0, 3.0, 4.0],
    );
    let before = mass_weighted_centre(&body);
    for _ in 0..8 {
        body.solve_volume_constraints();
    }
    let after = mass_weighted_centre(&body);
    for (axis, (&a, &b)) in after.iter().zip(before.iter()).enumerate() {
        assert!(
            (a - b).abs() < 1.0e-12,
            "sum of m * x must be invariant, axis {axis}: {a} vs {b}"
        );
    }
}

#[test]
fn volume_projection_restores_a_negatively_oriented_tet() {
    // The rest configuration itself has negative orientation: the tetrahedra
    // were listed with swapped vertices.
    let rest = skew_tet();
    let mut body = body_at(rest, compressed_to(rest, 0.8), [0, 2, 1, 3], [1.0; 4]);
    assert_eq!(body.rest_orientations, [-1.0]);
    assert!(body_signed_volume(&body) < 0.0);
    let v_rest = body.rest_volumes[0];
    let mut previous = body_volume_error(&body);
    for iteration in 0..10 {
        body.solve_volume_constraints();
        let error = body_volume_error(&body);
        assert!(
            error < previous,
            "iteration {iteration}: a negatively oriented tetrahedron must be inflated, {previous} -> {error}"
        );
        previous = error;
    }
    assert!(previous < 0.002 * v_rest, "left {}", previous / v_rest);
    assert!(
        body_signed_volume(&body) < 0.0,
        "the rest orientation must be kept"
    );
}

#[test]
fn volume_projection_recovers_a_fully_inverted_tet() {
    // Reflect the apex through the base plane: |V| equals the rest volume, so
    // a constraint on |V| sees nothing to correct and leaves the tetrahedron
    // inside out. The signed constraint pushes it back.
    let rest = unit_tet();
    let mut inverted = rest;
    inverted[3] = [0.0, 0.0, -1.0];
    let mut body = body_at(rest, inverted, [0, 1, 2, 3], [1.0; 4]);
    body.volume_stiffness = 1.0;
    let v_rest = body.rest_volumes[0];
    assert!(
        (body_signed_volume(&body) + v_rest).abs() < 1.0e-15,
        "the starting configuration is the mirror image of the rest shape"
    );
    let apex_before = body.particles[3].position;
    body.solve_volume_constraints();
    assert_ne!(
        body.particles[3].position, apex_before,
        "an inverted tetrahedron must be moved by the first projection"
    );
    for _ in 0..80 {
        body.solve_volume_constraints();
    }
    let volume = body_signed_volume(&body);
    assert!(
        volume > 0.0,
        "the orientation must be restored, signed volume {volume}"
    );
    assert!(
        (volume - v_rest).abs() < 1.0e-6 * v_rest,
        "restored volume {volume} vs rest {v_rest}"
    );
}

#[test]
fn volume_projection_keeps_fixed_particles_in_place_and_corrects_the_rest() {
    let rest = unit_tet();
    let mut body = body_at(rest, compressed_to(rest, 0.6), [0, 1, 2, 3], [1.0; 4]);
    body.particles[0].fixed = true;
    let anchored = body.particles[0].position;
    let initial_error = body_volume_error(&body);
    for _ in 0..20 {
        body.solve_volume_constraints();
    }
    assert_eq!(body.particles[0].position, anchored);
    assert!(body_volume_error(&body) < 0.05 * initial_error);
}

#[test]
fn tetrahedron_pushed_without_an_orientation_entry_uses_its_current_sign() {
    let rest = skew_tet();
    let mut body = body_at(rest, rest, [0, 1, 2, 3], [1.0; 4]);
    // A second tetrahedron registered by hand, without a rest orientation.
    body.particles.push(particle([0.6, 0.55, -0.9], 1.0));
    body.tetrahedra.push([4, 1, 2, 3]);
    let apex_volume = tet_volume_signed(
        body.particles[4].position,
        body.particles[1].position,
        body.particles[2].position,
        body.particles[3].position,
    );
    body.rest_volumes.push(apex_volume.abs());
    assert_eq!(body.rest_orientations.len(), 1);
    body.particles[4].position = [0.6, 0.55, -1.4];
    let before = body.particles[4].position;
    body.solve_volume_constraints();
    assert_ne!(
        body.particles[4].position, before,
        "the hand-registered tetrahedron must take part in the projection"
    );
    for p in &body.particles {
        assert!(p.position.iter().all(|c| c.is_finite()));
    }
}

// ---------------------------------------------------------------------------
// TetrahedralSoftBody::solve_volume_with_pressure
// ---------------------------------------------------------------------------

#[test]
fn pressurised_projection_inflates_a_negatively_oriented_compressed_tet() {
    let rest = skew_tet();
    let mut body = body_at(rest, compressed_to(rest, 0.8), [0, 2, 1, 3], [1.0; 4]);
    let v_rest = body.rest_volumes[0];
    let initial_error = body_volume_error(&body);
    let mut previous = initial_error;
    for iteration in 0..10 {
        body.solve_volume_with_pressure(0.01);
        let error = body_volume_error(&body);
        assert!(
            error < previous,
            "iteration {iteration}: the tetrahedron must be inflated, {previous} -> {error}"
        );
        previous = error;
    }
    assert!(
        previous < 0.05 * v_rest,
        "left {} of the rest volume",
        previous / v_rest
    );
}

#[test]
fn pressurised_projection_conserves_the_centroid_of_a_free_equal_mass_tet() {
    let rest = skew_tet();
    let mut body = body_at(rest, compressed_to(rest, 0.8), [0, 1, 2, 3], [1.0; 4]);
    let before = mass_weighted_centre(&body);
    for _ in 0..10 {
        body.solve_volume_with_pressure(0.01);
    }
    let after = mass_weighted_centre(&body);
    for (axis, (&a, &b)) in after.iter().zip(before.iter()).enumerate() {
        assert!((a - b).abs() < 1.0e-12, "axis {axis}: {a} vs {b}");
    }
}
