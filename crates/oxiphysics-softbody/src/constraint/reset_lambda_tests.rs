// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Per-constraint tests for [`SoftConstraint::reset_lambda`].
//!
//! The XPBD solvers that drive constraints as `Box<dyn SoftConstraint>` can
//! only reach the Lagrange multiplier through the trait.  Every implementor
//! that carries a multiplier must therefore override the (no-op) default.
//!
//! The check is behavioural: a constraint whose multiplier history was wiped
//! through the trait object has to project exactly like a freshly built one.

use super::SoftConstraint;
use crate::constraint as soft;
use crate::particle::SoftParticle;
use crate::pbd_constraints as pbd;
use oxiphysics_core::math::{Real, Vec3};

const DT: Real = 1.0 / 60.0;

fn particles_at(points: &[[Real; 3]]) -> Vec<SoftParticle> {
    points
        .iter()
        .map(|&[x, y, z]| SoftParticle::new(Vec3::new(x, y, z), 1.0))
        .collect()
}

/// Largest per-particle position difference between two particle sets.
fn max_deviation(a: &[SoftParticle], b: &[SoftParticle]) -> Real {
    a.iter()
        .zip(b.iter())
        .map(|(pa, pb)| (pa.position - pb.position).norm())
        .fold(0.0, Real::max)
}

/// What [`probe`] measured about one constraint.
struct Probe {
    /// How far a single fresh projection moves the particles.
    movement: Real,
    /// How far the projection of a constraint with a stale multiplier
    /// deviates from the fresh projection (the multiplier's influence).
    stale_deviation: Real,
}

/// Verifies that `reset_lambda` called through `dyn SoftConstraint` restores
/// fresh-constraint behaviour, and reports how much the scenario exercised the
/// multiplier.
fn probe<C, F>(make: F, start: &[SoftParticle]) -> Probe
where
    C: SoftConstraint + Clone,
    F: Fn() -> C,
{
    // Reference: a fresh constraint projects once.
    let mut fresh = make();
    let mut expected = start.to_vec();
    fresh.project(&mut expected, DT);

    // A second constraint accumulates multiplier history over several
    // projections of its own particle copy.
    let mut used = make();
    let mut scratch = start.to_vec();
    for _ in 0..4 {
        used.project(&mut scratch, DT);
    }

    // Without a reset the history leaks into the next projection.
    let mut stale = used.clone();
    let mut stale_particles = start.to_vec();
    stale.project(&mut stale_particles, DT);

    // Reset through the trait object, as the solvers do.
    {
        let dynamic: &mut dyn SoftConstraint = &mut used;
        dynamic.reset_lambda();
    }
    let mut actual = start.to_vec();
    used.project(&mut actual, DT);

    let residual = max_deviation(&actual, &expected);
    assert!(
        residual < 1.0e-12,
        "projection after reset_lambda deviates from a fresh constraint by {residual:e}"
    );

    Probe {
        movement: max_deviation(&expected, start),
        stale_deviation: max_deviation(&stale_particles, &expected),
    }
}

/// `probe` plus the requirement that the scenario really exercises the
/// multiplier (so a missing override would be caught).
fn probe_sensitive<C, F>(make: F, start: &[SoftParticle])
where
    C: SoftConstraint + Clone,
    F: Fn() -> C,
{
    let result = probe(make, start);
    assert!(
        result.movement > 1.0e-6,
        "scenario must violate the constraint (movement {:e})",
        result.movement
    );
    assert!(
        result.stale_deviation > 1.0e-9,
        "a stale multiplier must change the projection (deviation {:e}), otherwise the test is vacuous",
        result.stale_deviation
    );
}

fn flat_hinge() -> Vec<SoftParticle> {
    particles_at(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, -1.0, 0.0],
    ])
}

fn folded_hinge() -> Vec<SoftParticle> {
    particles_at(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.5],
        [0.5, -1.0, 0.0],
    ])
}

fn rest_tet() -> Vec<SoftParticle> {
    particles_at(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
}

fn stretched_tet() -> Vec<SoftParticle> {
    particles_at(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.5],
    ])
}

fn rest_triangle() -> Vec<SoftParticle> {
    particles_at(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
}

fn stretched_triangle() -> Vec<SoftParticle> {
    particles_at(&[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.4, 0.0]])
}

/// Closed tetrahedron surface (outward winding) used by the volume constraints.
fn tet_faces() -> Vec<[usize; 3]> {
    vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
}

// ---------------------------------------------------------------------------
// constraint::* implementors
// ---------------------------------------------------------------------------

#[test]
fn reset_distance_constraint() {
    let start = particles_at(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    probe_sensitive(|| soft::DistanceConstraint::new(0, 1, 1.0, 1.0e-2), &start);
}

#[test]
fn reset_bending_constraint() {
    let rest = flat_hinge();
    probe_sensitive(
        || soft::BendingConstraint::from_particles([0, 1, 2, 3], &rest, 1.0e-2),
        &folded_hinge(),
    );
}

#[test]
fn reset_volume_constraint() {
    let rest_volume = 1.0 / 6.0;
    probe_sensitive(
        || soft::VolumeConstraint::new([0, 1, 2, 3], rest_volume, 1.0e-2),
        &stretched_tet(),
    );
}

#[test]
fn reset_area_constraint() {
    probe_sensitive(
        || soft::AreaConstraint::new([0, 1, 2], 0.5, 1.0e-2),
        &stretched_triangle(),
    );
}

#[test]
fn reset_lra_constraint() {
    let start = particles_at(&[[2.0, 0.0, 0.0]]);
    probe_sensitive(
        || soft::LraConstraint::new(0, Vec3::zeros(), 1.0, 1.0e-2),
        &start,
    );
}

#[test]
fn reset_rigid_binding_constraint() {
    let start = particles_at(&[[1.0, 0.0, 0.0]]);
    probe_sensitive(
        || soft::RigidBindingConstraint::new(0, Vec3::zeros(), 1.0e-2),
        &start,
    );
}

#[test]
fn reset_motor_constraint() {
    let start = particles_at(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    probe_sensitive(
        || soft::MotorConstraint::new(0, 1, 1.0, 0.5, 1.0e-2),
        &start,
    );
}

#[test]
fn reset_surface_tension_constraint() {
    let rest = rest_triangle();
    probe_sensitive(
        || soft::SurfaceTensionConstraint::from_particles([0, 1, 2], &rest, 1.0, 1.0e-2),
        &stretched_triangle(),
    );
}

#[test]
fn reset_isometric_bending_constraint() {
    let rest = flat_hinge();
    probe(
        || soft::IsometricBendingConstraint::from_particles([0, 1, 2, 3], &rest, 1.0, 1.0e-2),
        &folded_hinge(),
    );
}

#[test]
fn reset_isometric_bending_constraint_v2() {
    let rest = flat_hinge();
    probe(
        || soft::IsometricBendingConstraintV2::from_particles([0, 1, 2, 3], &rest, 1.0, 1.0e-2),
        &folded_hinge(),
    );
}

#[test]
fn reset_balloon_pressure_constraint() {
    let start = particles_at(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    probe(
        || soft::BalloonPressureConstraint::new([0, 1, 2], 1.0, 1.0, 1.0e-2),
        &start,
    );
}

#[test]
fn reset_neo_hookean_constraint() {
    let rest = rest_tet();
    probe(
        || soft::NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, 1.0e3, 1.0e3, 1.0e-4),
        &stretched_tet(),
    );
}

#[test]
fn reset_strain_limiting_constraint() {
    let rest = rest_triangle();
    probe(
        || soft::StrainLimitingConstraint::from_particles([0, 1, 2], &rest, 1.1, 0.9, 1.0e-3),
        &stretched_triangle(),
    );
}

#[test]
fn constraints_without_multiplier_keep_default_reset() {
    // Collision and shape matching carry no Lagrange multiplier: the default
    // no-op must leave their behaviour untouched.
    let below_plane = particles_at(&[[0.0, -1.0, 0.0]]);
    let collision = probe(
        || soft::CollisionConstraint::new(0, Vec3::new(0.0, 1.0, 0.0), Vec3::zeros()),
        &below_plane,
    );
    assert!(
        collision.movement > 0.5,
        "particle must be pushed out of the plane"
    );
    assert!(
        collision.stale_deviation < 1.0e-15,
        "collision projection has no history"
    );

    let rest = rest_tet();
    let shape = probe(
        || soft::ShapeMatchingConstraint::from_particles(vec![0, 1, 2, 3], &rest, 0.5),
        &stretched_tet(),
    );
    assert!(
        shape.movement > 1.0e-6,
        "shape matching must move the particles"
    );
    assert!(
        shape.stale_deviation < 1.0e-15,
        "shape matching projection has no history"
    );
}

// ---------------------------------------------------------------------------
// pbd_constraints::* implementors
// ---------------------------------------------------------------------------

#[test]
fn reset_pbd_isometric_bending_constraint() {
    let rest = flat_hinge();
    probe(
        || pbd::IsometricBendingConstraint::new([0, 1, 2, 3], &rest, 1.0e-2),
        &folded_hinge(),
    );
}

#[test]
fn reset_pbd_cst_strain_constraint() {
    let rest = rest_triangle();
    probe(
        || pbd::CSTStrainConstraint::new([0, 1, 2], &rest, 1.0e-2),
        &stretched_triangle(),
    );
}

#[test]
fn reset_pbd_neo_hookean_constraint() {
    let rest = rest_tet();
    probe(
        || pbd::NeoHookeanConstraint::new([0, 1, 2, 3], &rest, 1.0e4, 0.3, 1.0e-4),
        &stretched_tet(),
    );
}

#[test]
fn reset_pbd_volume_conservation_constraint() {
    let rest_volume = pbd::VolumeConservationConstraint::compute_volume(&tet_faces(), &rest_tet());
    probe_sensitive(
        || pbd::VolumeConservationConstraint::new(tet_faces(), rest_volume, 1.0e-2),
        &stretched_tet(),
    );
}

#[test]
fn reset_pbd_cable_constraint() {
    let start = particles_at(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    probe_sensitive(|| pbd::CableConstraint::new(0, 1, 1.0, 1.0e-2), &start);
}

#[test]
fn reset_pbd_rigid_body_constraint() {
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let start = particles_at(&[[3.0, 0.0, 0.0]]);
    probe_sensitive(
        || {
            pbd::RigidBodyConstraint::new(
                0,
                Vec3::zeros(),
                identity,
                Vec3::new(1.0, 0.0, 0.0),
                1.0e-2,
            )
        },
        &start,
    );
}

#[test]
fn reset_pbd_sewing_constraint() {
    let start = particles_at(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.5, 0.0],
    ]);
    probe_sensitive(
        || pbd::SewingConstraint::new(vec![[0, 1], [2, 3]], 1.0e-2),
        &start,
    );
}

#[test]
fn reset_pbd_pressure_constraint() {
    probe_sensitive(
        || pbd::PressureConstraint::new(tet_faces(), 1.0 / 6.0, 1.5, 1.0e-2),
        &rest_tet(),
    );
}

#[test]
fn reset_pbd_anchor_constraint() {
    let start = particles_at(&[[1.0, 0.0, 0.0]]);
    probe_sensitive(
        || pbd::AnchorConstraint::new(0, Vec3::zeros(), 1.0e-2),
        &start,
    );
}
