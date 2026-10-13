// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the dihedral [`BendingConstraint`].
//!
//! The hinge `x0 = (0, 0, 0)`, `x1 = (1, 0, 0)` with the wings
//! `x2 = (0.5, 1, z2)` and `x3 = (0.5, -1, z3)` is flat for `z2 = z3 = 0`.
//! Raising both wings (`z > 0`) or lowering both (`z < 0`) folds it by the same
//! amount in opposite senses.  The old projection pushed the wings along a
//! fixed face normal, so it restored one fold direction and deepened the other
//! and injected net momentum; the tests below pin the correct behaviour with
//! finite-difference gradients, mirror symmetry and conservation laws.

use super::functions::SoftConstraint;
use super::types::BendingConstraint;
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};
use std::f64::consts::PI;

fn v(x: Real, y: Real, z: Real) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Positions of the flat or folded hinge described in the module documentation.
fn hinge_points(z2: Real, z3: Real) -> [Vec3; 4] {
    [
        v(0.0, 0.0, 0.0),
        v(1.0, 0.0, 0.0),
        v(0.5, 1.0, z2),
        v(0.5, -1.0, z3),
    ]
}

/// A hinge whose edge is pinned and whose wings are free.
fn pinned_hinge(z2: Real, z3: Real) -> Vec<SoftParticle> {
    let [x0, x1, x2, x3] = hinge_points(z2, z3);
    vec![
        SoftParticle::new_static(x0),
        SoftParticle::new_static(x1),
        SoftParticle::new(x2, 1.0),
        SoftParticle::new(x3, 1.0),
    ]
}

fn free_hinge(points: [Vec3; 4], masses: [Real; 4]) -> Vec<SoftParticle> {
    points
        .iter()
        .zip(masses)
        .map(|(&p, m)| SoftParticle::new(p, m))
        .collect()
}

fn positions(particles: &[SoftParticle]) -> [Vec3; 4] {
    [
        particles[0].position,
        particles[1].position,
        particles[2].position,
        particles[3].position,
    ]
}

/// Signed fold angle of the first four particles.
fn fold_of(particles: &[SoftParticle]) -> Real {
    let p = positions(particles);
    BendingConstraint::hinge_fold(&p[0], &p[1], &p[2], &p[3])
        .expect("the test hinge is not degenerate")
        .angle
}

/// The acos formula that `compute_dihedral` used before it was rewritten.
fn legacy_dihedral(p: &[Vec3; 4]) -> Real {
    let e = p[1] - p[0];
    let n1 = (p[2] - p[0]).cross(&e);
    let n2 = (p[3] - p[0]).cross(&e);
    let cos_a = (n1 / n1.norm()).dot(&(n2 / n2.norm()));
    cos_a.clamp(-1.0, 1.0).acos()
}

/// Hinges in general position, flat, slightly and strongly folded.
fn generic_hinges() -> Vec<[Vec3; 4]> {
    vec![
        hinge_points(0.0, 0.0),
        hinge_points(0.3, 0.3),
        hinge_points(-0.4, -0.4),
        hinge_points(0.5, -0.1),
        [
            v(0.1, -0.2, 0.3),
            v(1.3, 0.1, -0.1),
            v(0.4, 1.2, 0.5),
            v(0.9, -0.9, 0.7),
        ],
        // Both wings on nearly the same side: a strongly folded hinge.
        [
            v(0.0, 0.0, 0.0),
            v(1.0, 0.2, 0.0),
            v(0.5, 1.0, 0.1),
            v(0.4, 0.6, 0.6),
        ],
    ]
}

fn rest_constraint(points: [Vec3; 4], compliance: Real) -> BendingConstraint {
    let particles = free_hinge(points, [1.0; 4]);
    BendingConstraint::from_particles([0, 1, 2, 3], &particles, compliance)
}

#[test]
fn fold_gradient_matches_central_differences() {
    let h = 1e-6;
    for (case, points) in generic_hinges().into_iter().enumerate() {
        let fold = BendingConstraint::hinge_fold(&points[0], &points[1], &points[2], &points[3])
            .expect("non-degenerate hinge");
        assert!(
            fold.angle.abs() < 3.0,
            "case {case}: stay off the branch cut"
        );
        for k in 0..4 {
            for axis in 0..3 {
                let mut plus = points;
                let mut minus = points;
                plus[k][axis] += h;
                minus[k][axis] -= h;
                let angle = |p: &[Vec3; 4]| {
                    BendingConstraint::hinge_fold(&p[0], &p[1], &p[2], &p[3])
                        .expect("non-degenerate hinge")
                        .angle
                };
                let numeric = (angle(&plus) - angle(&minus)) / (2.0 * h);
                let analytic = fold.gradient[k][axis];
                assert!(
                    (numeric - analytic).abs() < 1e-6 * (1.0 + analytic.abs()),
                    "case {case}, particle {k}, axis {axis}: analytic {analytic} vs numeric {numeric}"
                );
            }
        }
    }
}

#[test]
fn fold_gradient_is_translation_and_rotation_invariant() {
    for (case, points) in generic_hinges().into_iter().enumerate() {
        let fold = BendingConstraint::hinge_fold(&points[0], &points[1], &points[2], &points[3])
            .expect("non-degenerate hinge");
        let force = fold.gradient.iter().fold(Vec3::zeros(), |acc, g| acc + g);
        assert!(force.norm() < 1e-12, "case {case}: net force {force:?}");
        let torque = fold
            .gradient
            .iter()
            .zip(&points)
            .fold(Vec3::zeros(), |acc, (g, x)| acc + x.cross(g));
        assert!(torque.norm() < 1e-12, "case {case}: net torque {torque:?}");
    }
}

#[test]
fn fold_angle_is_signed_zero_when_flat_and_smooth_through_it() {
    assert!(fold_of(&free_hinge(hinge_points(0.0, 0.0), [1.0; 4])).abs() < 1e-15);
    // Raising both wings lowers phi, lowering them raises it, the tilted flat
    // plane z = 0.3 y keeps phi at zero.
    let up = fold_of(&free_hinge(hinge_points(0.3, 0.3), [1.0; 4]));
    let down = fold_of(&free_hinge(hinge_points(-0.3, -0.3), [1.0; 4]));
    let tilted = fold_of(&free_hinge(hinge_points(0.3, -0.3), [1.0; 4]));
    assert!(up < -0.5 && down > 0.5, "up {up}, down {down}");
    assert!((up + down).abs() < 1e-12, "the two folds mirror each other");
    assert!(tilted.abs() < 1e-12, "a tilted flat sheet is not folded");
    // Small folds scale linearly: phi ~ -(z2 + z3) for unit wing heights.
    let tiny = fold_of(&free_hinge(hinge_points(1e-7, 2e-7), [1.0; 4]));
    assert!((tiny + 3e-7).abs() < 1e-12, "tiny fold {tiny}");
    // Monotone across the flat state.
    let samples: Vec<Real> = (-10..=10)
        .map(|i| {
            fold_of(&free_hinge(
                hinge_points(0.05 * i as Real, 0.05 * i as Real),
                [1.0; 4],
            ))
        })
        .collect();
    assert!(samples.windows(2).all(|w| w[0] > w[1]), "{samples:?}");
}

#[test]
fn unsigned_dihedral_keeps_the_legacy_convention() {
    // The rest angle of a flat hinge is pi and folding (in either direction)
    // lowers it; the rewritten function agrees with the old acos formula.
    for (case, points) in generic_hinges().into_iter().enumerate() {
        let new =
            BendingConstraint::compute_dihedral(&points[0], &points[1], &points[2], &points[3]);
        let old = legacy_dihedral(&points);
        assert!(
            (new - old).abs() < 1e-6,
            "case {case}: new {new} vs legacy {old}"
        );
    }
    let flat = hinge_points(0.0, 0.0);
    let angle = BendingConstraint::compute_dihedral(&flat[0], &flat[1], &flat[2], &flat[3]);
    assert!((angle - PI).abs() < 1e-15, "flat hinge: {angle}");
    // Close to flat, acos loses half its digits; atan2 does not.
    let nearly = hinge_points(1e-9, 1e-9);
    let angle = BendingConstraint::compute_dihedral(&nearly[0], &nearly[1], &nearly[2], &nearly[3]);
    assert!((angle - (PI - 2e-9)).abs() < 1e-15, "nearly flat: {angle}");
    assert_eq!(rest_constraint(flat, 0.0).rest_angle, PI);
}

#[test]
fn degenerate_hinges_have_no_fold() {
    let ok = hinge_points(0.2, 0.1);
    let fold = |p: [Vec3; 4]| BendingConstraint::hinge_fold(&p[0], &p[1], &p[2], &p[3]);
    assert!(fold(ok).is_some());
    // Wing on the hinge line.
    let mut on_line = ok;
    on_line[2] = v(0.5, 0.0, 0.0);
    assert!(fold(on_line).is_none());
    on_line[2] = ok[2];
    on_line[3] = v(2.0, 0.0, 0.0);
    assert!(fold(on_line).is_none());
    // Zero-length hinge.
    let mut collapsed = ok;
    collapsed[1] = collapsed[0];
    assert!(fold(collapsed).is_none());
    // Non-finite coordinates.
    let mut broken = ok;
    broken[2] = v(f64::NAN, 1.0, 0.0);
    assert!(fold(broken).is_none());
    broken[2] = v(0.5, f64::INFINITY, 0.0);
    assert!(fold(broken).is_none());
    assert_eq!(
        BendingConstraint::compute_dihedral(&on_line[0], &on_line[1], &on_line[2], &on_line[3]),
        0.0
    );
}

#[test]
fn a_fold_is_restored_whichever_way_it_is_bent() {
    for (z2, z3) in [
        (0.3, 0.3),
        (-0.3, -0.3),
        (0.5, 0.1),
        (-0.5, -0.1),
        (0.2, 0.0),
        (0.0, -0.2),
        (0.4, -0.1),
        (-0.4, 0.1),
    ] {
        let mut particles = pinned_hinge(z2, z3);
        let start = positions(&particles);
        let mut constraint = BendingConstraint::new([0, 1, 2, 3], PI, 0.0);
        let before = fold_of(&particles).abs();
        assert!(before > 0.05, "({z2}, {z3}) must start folded");
        constraint.project(&mut particles, 1.0 / 60.0);
        let after = fold_of(&particles).abs();
        assert!(
            after < 0.25 * before,
            "({z2}, {z3}): one projection must remove most of the fold, {before} -> {after}"
        );
        for _ in 0..10 {
            constraint.reset_lambda();
            constraint.project(&mut particles, 1.0 / 60.0);
        }
        assert!(
            fold_of(&particles).abs() < 1e-9,
            "({z2}, {z3}) converges to flat, fold {}",
            fold_of(&particles)
        );
        let end = positions(&particles);
        assert_eq!(start[0], end[0], "pinned hinge end x0 moved");
        assert_eq!(start[1], end[1], "pinned hinge end x1 moved");
    }
}

#[test]
fn the_projection_is_mirror_symmetric() {
    // Reflecting the hinge through the z = 0 plane must reflect the result.
    for (z2, z3) in [(0.3, 0.3), (0.5, 0.1), (0.4, -0.2)] {
        let mut up = pinned_hinge(z2, z3);
        let mut down = pinned_hinge(-z2, -z3);
        let mut c_up = BendingConstraint::new([0, 1, 2, 3], PI, 1e-6);
        let mut c_down = c_up.clone();
        for _ in 0..3 {
            c_up.project(&mut up, 1.0 / 60.0);
            c_down.project(&mut down, 1.0 / 60.0);
        }
        for (a, b) in up.iter().zip(&down) {
            assert!(
                (a.position.x - b.position.x).abs() < 1e-12
                    && (a.position.y - b.position.y).abs() < 1e-12
                    && (a.position.z + b.position.z).abs() < 1e-12,
                "({z2}, {z3}): {:?} vs mirrored {:?}",
                a.position,
                b.position
            );
        }
    }
}

#[test]
fn a_projection_conserves_linear_and_angular_momentum() {
    let masses = [1.0, 2.0, 3.0, 4.0];
    for (case, points) in generic_hinges().into_iter().enumerate() {
        let mut particles = free_hinge(points, masses);
        // Ask for a different angle than the current one so that it acts.
        let current =
            BendingConstraint::compute_dihedral(&points[0], &points[1], &points[2], &points[3]);
        let mut constraint = BendingConstraint::new([0, 1, 2, 3], current - 0.2, 1e-5);
        constraint.project(&mut particles, 1.0 / 120.0);
        let mut linear = Vec3::zeros();
        let mut angular = Vec3::zeros();
        let mut moved: Real = 0.0;
        for (k, p) in particles.iter().enumerate() {
            let displacement = p.position - points[k];
            moved += displacement.norm();
            linear += displacement * masses[k];
            angular += points[k].cross(&(displacement * masses[k]));
        }
        assert!(moved > 1e-6, "case {case}: the constraint must act");
        assert!(
            linear.norm() < 1e-12,
            "case {case}: linear momentum {linear:?}"
        );
        assert!(
            angular.norm() < 1e-12,
            "case {case}: angular momentum {angular:?}"
        );
    }
}

#[test]
fn a_hinge_converges_to_its_rest_state_with_all_four_particles_free() {
    for (case, rest_points) in generic_hinges().into_iter().enumerate() {
        let mut constraint = rest_constraint(rest_points, 0.0);
        let rest_angle = constraint.rest_angle;
        // Disturb the rest shape by a smooth deformation and relax it.
        let mut particles = free_hinge(rest_points, [1.0, 2.0, 1.0, 3.0]);
        for (k, p) in particles.iter_mut().enumerate() {
            let phase = k as Real + case as Real;
            p.position += v(
                0.03 * phase.sin(),
                0.02 * (2.0 * phase).cos(),
                0.05 * (3.0 * phase).sin(),
            );
        }
        for _ in 0..200 {
            constraint.reset_lambda();
            constraint.project(&mut particles, 1.0 / 60.0);
        }
        let p = positions(&particles);
        let angle = BendingConstraint::compute_dihedral(&p[0], &p[1], &p[2], &p[3]);
        assert!(
            (angle - rest_angle).abs() < 1e-9,
            "case {case}: angle {angle} vs rest {rest_angle}"
        );
    }
}

#[test]
fn the_rest_state_is_a_fixed_point() {
    for (case, points) in generic_hinges().into_iter().enumerate() {
        let mut constraint = rest_constraint(points, 1e-4);
        let mut particles = free_hinge(points, [1.0, 2.0, 3.0, 4.0]);
        constraint.project(&mut particles, 1.0 / 60.0);
        for (k, p) in particles.iter().enumerate() {
            assert!(
                (p.position - points[k]).norm() < 1e-12,
                "case {case}: particle {k} moved at rest"
            );
        }
        assert_eq!(constraint.lambda, 0.0, "case {case}");
    }
}

#[test]
fn a_creased_rest_angle_is_reached_from_either_side() {
    let crease = 0.8;
    let rest_angle = PI - crease;
    // From the flat state the hinge folds towards the lowered-wings side (the
    // tie at exactly zero is broken towards positive phi) and stops at |phi| = crease.
    let mut flat = pinned_hinge(0.0, 0.0);
    let mut constraint = BendingConstraint::new([0, 1, 2, 3], rest_angle, 0.0);
    for _ in 0..30 {
        constraint.reset_lambda();
        constraint.project(&mut flat, 1.0 / 60.0);
    }
    assert!(
        (fold_of(&flat).abs() - crease).abs() < 1e-9,
        "{}",
        fold_of(&flat)
    );
    // A slightly folded hinge keeps its sense and deepens to the crease angle.
    for sign in [1.0, -1.0] {
        let mut particles = pinned_hinge(-0.1 * sign, -0.1 * sign);
        let initial = fold_of(&particles);
        let mut constraint = BendingConstraint::new([0, 1, 2, 3], rest_angle, 0.0);
        for _ in 0..30 {
            constraint.reset_lambda();
            constraint.project(&mut particles, 1.0 / 60.0);
        }
        let end = fold_of(&particles);
        assert_eq!(
            end.signum(),
            initial.signum(),
            "the fold must keep its sense"
        );
        assert!((end.abs() - crease).abs() < 1e-9, "sign {sign}: {end}");
    }
    // Over-folded: the hinge opens back to the crease angle.
    let mut particles = pinned_hinge(-1.5, -1.5);
    let mut constraint = BendingConstraint::new([0, 1, 2, 3], rest_angle, 0.0);
    assert!(fold_of(&particles).abs() > crease);
    for _ in 0..30 {
        constraint.reset_lambda();
        constraint.project(&mut particles, 1.0 / 60.0);
    }
    assert!((fold_of(&particles).abs() - crease).abs() < 1e-9);
}

#[test]
fn compliance_softens_the_correction_and_multipliers_accumulate() {
    let correction = |compliance: Real| -> (Real, BendingConstraint) {
        let mut particles = pinned_hinge(0.3, 0.3);
        let mut constraint = BendingConstraint::new([0, 1, 2, 3], PI, compliance);
        constraint.project(&mut particles, 1.0 / 60.0);
        (fold_of(&particles).abs(), constraint)
    };
    let initial = fold_of(&pinned_hinge(0.3, 0.3)).abs();
    let (stiff, stiff_constraint) = correction(0.0);
    let (medium, _) = correction(1e-4);
    let (soft, mut soft_constraint) = correction(1e-1);
    assert!(
        stiff < medium && medium < soft && soft < initial,
        "{stiff} < {medium} < {soft} < {initial}"
    );
    assert!(stiff_constraint.lambda != 0.0 && soft_constraint.lambda != 0.0);
    // The multiplier of the first projection is -C / (W + alpha_tilde).
    let particles = pinned_hinge(0.3, 0.3);
    let p = positions(&particles);
    let fold = BendingConstraint::hinge_fold(&p[0], &p[1], &p[2], &p[3]).expect("hinge");
    let w: Real = fold.gradient[2].norm_squared() + fold.gradient[3].norm_squared();
    let c = (PI - fold.angle.abs()) - PI;
    let alpha_tilde = 1e-1 / (1.0 / 60.0_f64).powi(2);
    assert!(
        (soft_constraint.lambda - (-c / (w + alpha_tilde))).abs() < 1e-12,
        "lambda {} vs {}",
        soft_constraint.lambda,
        -c / (w + alpha_tilde)
    );
    soft_constraint.reset_lambda();
    assert_eq!(soft_constraint.lambda, 0.0);
    let boxed: &mut dyn SoftConstraint = &mut soft_constraint;
    boxed.reset_lambda();
    assert_eq!(soft_constraint.lambda, 0.0);
}

#[test]
fn invalid_input_leaves_the_particles_untouched() {
    let folded = pinned_hinge(0.3, 0.3);
    let snapshot = positions(&folded);
    let unchanged = |particles: &[SoftParticle]| {
        positions(particles)
            .iter()
            .zip(&snapshot)
            .all(|(a, b)| a == b)
    };
    let mut constraint = BendingConstraint::new([0, 1, 2, 3], PI, 0.0);
    for bad_dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut particles = folded.clone();
        constraint.project(&mut particles, bad_dt);
        assert!(unchanged(&particles), "dt = {bad_dt}");
    }
    assert_eq!(constraint.lambda, 0.0);

    let mut out_of_range = BendingConstraint::new([0, 1, 2, 9], PI, 0.0);
    let mut particles = folded.clone();
    out_of_range.project(&mut particles, 1.0 / 60.0);
    assert!(unchanged(&particles));

    // A wing on the hinge line, the situation of the old collinear stencils.
    let mut particles = folded.clone();
    particles[2].position = v(0.5, 0.0, 0.0);
    let on_line = positions(&particles);
    constraint.project(&mut particles, 1.0 / 60.0);
    assert_eq!(positions(&particles), on_line);

    // Everything pinned: nothing can move and no multiplier accumulates.
    let mut frozen: Vec<SoftParticle> = folded
        .iter()
        .map(|p| SoftParticle::new_static(p.position))
        .collect();
    constraint.project(&mut frozen, 1.0 / 60.0);
    assert!(unchanged(&frozen));
    assert_eq!(constraint.lambda, 0.0);
}
