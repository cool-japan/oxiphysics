// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the XPBD Lagrange-multiplier reset (dyn-constraint
//! solvers) and for the Gauss-Seidel over-relaxation.
//!
//! The reference scenario is a particle of mass `m` hanging from a pinned
//! anchor through a compliant distance constraint.  XPBD has a closed-form
//! static solution for it: the constraint force balances gravity, so the
//! constraint is stretched by exactly `m * g * compliance`.  A solver that
//! carries the Lagrange multiplier over from sub-step to sub-step keeps
//! softening the constraint and the stretch grows without bound.

use super::*;
use crate::constraint::DistanceConstraint;
use crate::pbd_constraints::{AnchorConstraint, ConstraintSolver};
use oxiphysics_core::math::Vec3;

const MASS: Real = 1.0;
const GRAVITY: Real = 9.81;
const COMPLIANCE: Real = 1.0e-3;
const REST_LENGTH: Real = 1.0;
const FRAME_DT: Real = 1.0 / 60.0;
/// Static equilibrium stretch `m * g * compliance` of the hanging particle.
const EXPECTED_STRETCH: Real = MASS * GRAVITY * COMPLIANCE;

/// Pinned anchor at the origin plus a particle hanging `REST_LENGTH` below it.
fn hanging_particles() -> Vec<SoftParticle> {
    let mut particles = vec![
        SoftParticle::new_static(Vec3::zeros()),
        SoftParticle::new(Vec3::new(0.0, -REST_LENGTH, 0.0), MASS),
    ];
    particles[1].external_force = Vec3::new(0.0, -MASS * GRAVITY, 0.0);
    particles
}

fn hanging_constraints() -> Vec<Box<dyn SoftConstraint>> {
    vec![Box::new(DistanceConstraint::new(
        0,
        1,
        REST_LENGTH,
        COMPLIANCE,
    ))]
}

/// Elongation of the hanging distance constraint.
fn stretch(particles: &[SoftParticle]) -> Real {
    (particles[1].position - particles[0].position).norm() - REST_LENGTH
}

fn assert_static_stretch(label: &str, value: Real) {
    assert!(
        (value - EXPECTED_STRETCH).abs() <= 0.05 * EXPECTED_STRETCH,
        "{label}: stretch {value:.6e} should be within 5% of m*g*compliance = {EXPECTED_STRETCH:.6e}"
    );
}

// ---------------------------------------------------------------------------
// Lagrange multiplier reset
// ---------------------------------------------------------------------------

#[test]
fn xpbd_solver_hanging_particle_settles_at_static_stretch_and_stays() {
    let mut body = SoftBody::from_particles(hanging_particles());
    let mut constraints = hanging_constraints();
    let mut solver = XpbdSolver::with_iterations(4, 1);
    // Never fall asleep: sleeping would freeze (and thereby hide) any drift.
    solver.sleep_threshold = 0.0;

    let mut stretch_after_1_frame = 0.0;
    let mut stretch_after_120_frames = 0.0;
    for frame in 1..=600 {
        solver.solve(&mut body, &mut constraints, FRAME_DT);
        if frame == 1 {
            stretch_after_1_frame = stretch(&body.particles);
        }
        if frame == 120 {
            stretch_after_120_frames = stretch(&body.particles);
        }
    }
    let stretch_after_600_frames = stretch(&body.particles);

    // One frame is far from equilibrium, but the stretch must already be the
    // compliant one rather than a stale-multiplier artefact.
    assert!(
        stretch_after_1_frame > 0.0 && stretch_after_1_frame < 2.0 * EXPECTED_STRETCH,
        "stretch after one frame {stretch_after_1_frame:.6e} must lie between 0 and twice the static stretch"
    );
    assert_static_stretch("after 120 frames", stretch_after_120_frames);
    assert_static_stretch("after 600 frames", stretch_after_600_frames);
    assert!(
        (stretch_after_600_frames - stretch_after_120_frames).abs() < 1.0e-4 * REST_LENGTH,
        "stretch must not drift: {stretch_after_120_frames:.6e} -> {stretch_after_600_frames:.6e}"
    );
}

#[test]
fn xpbd_solver_stretch_scales_with_compliance_and_mass() {
    // Static stretch = m * g * compliance, so doubling either doubles it.
    for (mass, compliance) in [(1.0, 5.0e-4), (2.0, 1.0e-3), (0.5, 4.0e-3)] {
        let mut particles = vec![
            SoftParticle::new_static(Vec3::zeros()),
            SoftParticle::new(Vec3::new(0.0, -REST_LENGTH, 0.0), mass),
        ];
        particles[1].external_force = Vec3::new(0.0, -mass * GRAVITY, 0.0);
        let mut body = SoftBody::from_particles(particles);
        let mut constraints: Vec<Box<dyn SoftConstraint>> = vec![Box::new(
            DistanceConstraint::new(0, 1, REST_LENGTH, compliance),
        )];
        let mut solver = XpbdSolver::with_iterations(4, 1);
        solver.sleep_threshold = 0.0;
        for _ in 0..600 {
            solver.solve(&mut body, &mut constraints, FRAME_DT);
        }
        let expected = mass * GRAVITY * compliance;
        let got = stretch(&body.particles);
        assert!(
            (got - expected).abs() <= 0.05 * expected,
            "m={mass}, compliance={compliance}: stretch {got:.6e} vs expected {expected:.6e}"
        );
    }
}

#[test]
fn solve_with_ordering_hanging_particle_settles_for_every_ordering() {
    for ordering in [
        ConstraintOrdering::Sequential,
        ConstraintOrdering::Alternating,
        ConstraintOrdering::Shuffled,
    ] {
        let mut body = SoftBody::from_particles(hanging_particles());
        let mut constraints = hanging_constraints();
        let mut solver = XpbdSolver::with_iterations(4, 2);
        solver.sleep_threshold = 0.0;
        for _ in 0..600 {
            solve_with_ordering(&mut solver, &mut body, &mut constraints, FRAME_DT, ordering);
        }
        assert_static_stretch(&format!("{ordering:?}"), stretch(&body.particles));
    }
}

#[test]
fn constraint_solver_hanging_particle_settles_at_static_stretch() {
    let mut particles = hanging_particles();
    // `ConstraintSolver` takes gravity as an acceleration, not as a force.
    particles[1].external_force = Vec3::zeros();
    let mut constraints = hanging_constraints();
    let solver = ConstraintSolver::new(1, 4);
    let gravity = Vec3::new(0.0, -GRAVITY, 0.0);
    for _ in 0..600 {
        solver.solve(&mut particles, &mut constraints, gravity, FRAME_DT);
    }
    assert_static_stretch("ConstraintSolver", stretch(&particles));
}

#[test]
fn constraint_solver_anchor_constraint_settles_at_static_offset() {
    // A compliant world-space pin behaves like a spring to the anchor point.
    let mut particles = vec![SoftParticle::new(Vec3::zeros(), MASS)];
    let mut constraints: Vec<Box<dyn SoftConstraint>> = vec![Box::new(AnchorConstraint::new(
        0,
        Vec3::zeros(),
        COMPLIANCE,
    ))];
    let solver = ConstraintSolver::new(1, 4);
    let gravity = Vec3::new(0.0, -GRAVITY, 0.0);
    for _ in 0..600 {
        solver.solve(&mut particles, &mut constraints, gravity, FRAME_DT);
    }
    assert_static_stretch("AnchorConstraint", particles[0].position.norm());
}

#[test]
fn constraint_solver_project_only_is_repeatable() {
    // `project_only` is one sub-step's projection phase: running it twice from
    // the same particle state must give the same result, i.e. the first call
    // may not leave a multiplier behind that the second call picks up.
    let make_particles = || {
        vec![
            SoftParticle::new_static(Vec3::zeros()),
            SoftParticle::new(Vec3::new(0.0, -2.0, 0.0), MASS),
        ]
    };
    let mut constraints = hanging_constraints();
    let solver = ConstraintSolver::new(3, 1);

    let mut first = make_particles();
    solver.project_only(&mut first, &mut constraints, FRAME_DT);
    let mut second = make_particles();
    solver.project_only(&mut second, &mut constraints, FRAME_DT);

    assert!(
        (first[1].position - second[1].position).norm() < 1.0e-14,
        "stale multiplier leaked between project_only calls: {:?} vs {:?}",
        first[1].position,
        second[1].position
    );
}

#[test]
fn gauss_seidel_solve_resets_multipliers_on_every_call() {
    let make_particles = || {
        let mut particles = vec![
            SoftParticle::new_static(Vec3::zeros()),
            SoftParticle::new(Vec3::new(0.0, -1.2, 0.0), MASS),
        ];
        for p in &mut particles {
            p.prev_position = p.position;
        }
        particles
    };
    let mut constraints = hanging_constraints();
    let gs = GaussSeidelSolver::new(3);

    let mut first = make_particles();
    gs.solve(&mut first, &mut constraints, FRAME_DT);
    let mut second = make_particles();
    gs.solve(&mut second, &mut constraints, FRAME_DT);

    assert!(
        (first[1].position - second[1].position).norm() < 1.0e-14,
        "stale multiplier leaked between solve calls: {:?} vs {:?}",
        first[1].position,
        second[1].position
    );
    // The compliant constraint is only partially satisfied after one sub-step.
    let remaining = stretch(&first);
    assert!(
        remaining > 0.0 && remaining < 0.2,
        "compliance should leave part of the 0.2 violation: remaining stretch {remaining}"
    );
}

#[test]
fn dyn_reset_lambda_matches_inherent_reset() {
    // Two identical constraints, one driven through `dyn SoftConstraint`.
    let mut particles_a = hanging_particles();
    let mut particles_b = hanging_particles();
    particles_a[1].position.y = -1.5;
    particles_b[1].position.y = -1.5;

    let mut inherent = DistanceConstraint::new(0, 1, REST_LENGTH, COMPLIANCE);
    let mut boxed: Box<dyn SoftConstraint> =
        Box::new(DistanceConstraint::new(0, 1, REST_LENGTH, COMPLIANCE));

    for _ in 0..4 {
        inherent.reset_lambda();
        inherent.project(&mut particles_a, FRAME_DT);
        boxed.reset_lambda();
        boxed.project(&mut particles_b, FRAME_DT);
    }
    assert!(
        (particles_a[1].position - particles_b[1].position).norm() < 1.0e-15,
        "dyn reset must behave like the inherent reset"
    );
}

// ---------------------------------------------------------------------------
// Gauss-Seidel over-relaxation
// ---------------------------------------------------------------------------

#[test]
fn sor_never_rescales_inertial_displacement() {
    // One unconstrained particle that has already moved away from its
    // sub-step start position: no constraint touches it, so no sweep may move
    // it, whatever the relaxation factor.
    for omega in [0.5, 1.0, 1.5, 1.9] {
        let mut particles = vec![SoftParticle::new(Vec3::new(1.0, 2.0, 3.0), MASS)];
        particles[0].prev_position = Vec3::new(0.9, 2.1, 3.0);
        let mut constraints: Vec<Box<dyn SoftConstraint>> = Vec::new();
        GaussSeidelSolver::with_sor(10, omega).solve(&mut particles, &mut constraints, FRAME_DT);
        assert_eq!(
            particles[0].position,
            Vec3::new(1.0, 2.0, 3.0),
            "omega={omega}: inertial displacement must be left untouched"
        );
    }
}

#[test]
fn sor_identical_results_for_unconstrained_particle() {
    let run = |omega: Real| {
        let mut particles = vec![SoftParticle::new(Vec3::new(0.3, -0.2, 0.7), MASS)];
        particles[0].prev_position = Vec3::new(0.0, 0.0, 0.0);
        let mut constraints: Vec<Box<dyn SoftConstraint>> = Vec::new();
        GaussSeidelSolver::with_sor(5, omega).solve(&mut particles, &mut constraints, FRAME_DT);
        particles[0].position
    };
    assert_eq!(run(1.0), run(1.5));
}

/// Initial elongation of [`violated_link`] (2.4 m apart, 1 m rest length).
const LINK_INITIAL_ERROR: Real = 1.4;

/// A pinned particle and a dynamic one 2.4 m away, joined by a rigid 1 m link.
///
/// The error never reaches the degenerate coincident state for relaxation
/// factors in `(0, 2)`: after the first sweep `|error| <= 1.4 * |1 - omega| < 1`.
fn violated_link() -> (Vec<SoftParticle>, Vec<Box<dyn SoftConstraint>>) {
    let mut particles = vec![
        SoftParticle::new_static(Vec3::zeros()),
        SoftParticle::new(Vec3::new(1.0 + LINK_INITIAL_ERROR, 0.0, 0.0), MASS),
    ];
    for p in &mut particles {
        p.prev_position = p.position;
    }
    let constraints: Vec<Box<dyn SoftConstraint>> =
        vec![Box::new(DistanceConstraint::new(0, 1, 1.0, 0.0))];
    (particles, constraints)
}

#[test]
fn sor_sweeps_converge_geometrically_for_valid_omega() {
    // A rigid link is satisfied exactly by one projection, so a sweep with
    // relaxation omega multiplies the remaining error by (1 - omega).
    for (omega, sweeps, tolerance) in [(1.0, 1, 1.0e-12), (1.5, 40, 1.0e-9), (0.5, 40, 1.0e-9)] {
        let (mut particles, mut constraints) = violated_link();
        GaussSeidelSolver::try_with_sor(sweeps, omega)
            .expect("omega inside (0, 2) is valid")
            .solve(&mut particles, &mut constraints, FRAME_DT);
        let length = (particles[1].position - particles[0].position).norm();
        assert!(
            (length - 1.0).abs() < tolerance,
            "omega={omega}: link length {length} should converge to 1"
        );
    }
}

#[test]
fn sor_error_follows_one_minus_omega_per_sweep() {
    let (mut particles, mut constraints) = violated_link();
    let omega = 1.5;
    let gs = GaussSeidelSolver::with_sor(1, omega);
    let mut previous_error = LINK_INITIAL_ERROR;
    for sweep in 1..=6 {
        gs.solve(&mut particles, &mut constraints, FRAME_DT);
        let error = (particles[1].position - particles[0].position).norm() - 1.0;
        let expected = previous_error * (1.0 - omega).abs();
        assert!(
            (error.abs() - expected).abs() < 1.0e-12,
            "sweep {sweep}: |error| {} vs expected {expected}",
            error.abs()
        );
        previous_error = error.abs();
    }
}

#[test]
fn sor_keeps_pinned_particles_fixed() {
    let (mut particles, mut constraints) = violated_link();
    let pinned = particles[0].position;
    GaussSeidelSolver::with_sor(8, 1.7).solve(&mut particles, &mut constraints, FRAME_DT);
    assert_eq!(particles[0].position, pinned);
}

#[test]
fn sor_matches_plain_gauss_seidel_when_omega_is_one() {
    let (mut plain, mut plain_constraints) = violated_link();
    let (mut relaxed, mut relaxed_constraints) = violated_link();
    GaussSeidelSolver::new(4).solve(&mut plain, &mut plain_constraints, FRAME_DT);
    GaussSeidelSolver::with_sor(4, 1.0).solve(&mut relaxed, &mut relaxed_constraints, FRAME_DT);
    assert_eq!(plain[1].position, relaxed[1].position);
}

#[test]
fn try_with_sor_accepts_only_open_interval_zero_to_two() {
    for omega in [0.1, 0.5, 1.0, 1.5, 1.99] {
        let gs = GaussSeidelSolver::try_with_sor(7, omega).expect("valid omega");
        assert_eq!(gs.iterations, 7);
        assert_eq!(gs.omega, omega);
    }
    for omega in [
        0.0,
        -0.5,
        2.0,
        2.5,
        Real::NAN,
        Real::INFINITY,
        Real::NEG_INFINITY,
    ] {
        match GaussSeidelSolver::try_with_sor(7, omega) {
            Err(Error::InvalidParameter(message)) => {
                assert!(
                    message.contains("omega"),
                    "message should name the parameter: {message}"
                );
            }
            other => panic!("omega={omega} must be rejected, got {other:?}"),
        }
    }
}
