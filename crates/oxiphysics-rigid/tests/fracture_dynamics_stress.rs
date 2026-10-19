// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for `AggregateBody` / `FractureSimulation::step`.
//!
//! The nominal stress is `|F| / A` over the body's cross-section (it used to be
//! `F / (mass * 0.001)`), the applied force accelerates the body, and the
//! fragments of a broken body take the body's own density and extents.

use oxiphysics_rigid::Error;
use oxiphysics_rigid::fracture_dynamics::{
    AggregateBody, BreakableJoint, CrackPlane, FractureSimulation, FractureThreshold, FragmentSet,
};

const DENSITY: f64 = 2500.0;

fn body_with_area(mass: f64, area: f64) -> AggregateBody {
    AggregateBody::new(
        0,
        [0.0; 3],
        vec![0.5 * mass, 0.5 * mass],
        DENSITY,
        [0.1; 3],
        area,
    )
    .expect("valid aggregate body")
}

/// Number of fracture events after a single step under `force` newtons.
fn events_after_one_step(mass: f64, area: f64, force: f64, stress_threshold: f64) -> usize {
    let mut sim = FractureSimulation::new();
    let mut body = body_with_area(mass, area);
    // The joint gives way as soon as it is loaded at all.
    body.add_joint(BreakableJoint::new(0, 0, 1, 1.0, 1.0e9, 1.0e9));
    sim.add_body(body);
    let threshold = FractureThreshold::new(stress_threshold, [1.0, 0.0, 0.0], 4);
    sim.step(0.01, [force, 0.0, 0.0], &threshold)
}

// ── stress ──────────────────────────────────────────────────────────────────

#[test]
fn stress_is_force_magnitude_over_cross_section() {
    let body = body_with_area(2.0, 0.5);
    // |(3, 4, 0)| = 5 N over 0.5 m² is 10 Pa.
    assert!((body.stress([3.0, 4.0, 0.0]) - 10.0).abs() < 1.0e-12);
    assert_eq!(body.stress([0.0; 3]), 0.0);
}

#[test]
fn doubling_the_area_halves_the_stress() {
    let thin = body_with_area(3.0, 0.02);
    let wide = body_with_area(3.0, 0.04);
    let force = [120.0, -50.0, 10.0];
    let ratio = thin.stress(force) / wide.stress(force);
    assert!((ratio - 2.0).abs() < 1.0e-12, "ratio = {ratio}");
}

#[test]
fn stress_does_not_depend_on_the_mass() {
    let light = body_with_area(1.0, 0.04);
    let heavy = body_with_area(100.0, 0.04);
    let force = [10.0, 0.0, 0.0];
    assert_eq!(light.stress(force), heavy.stress(force));
    assert!((light.stress(force) - 250.0).abs() < 1.0e-12);
}

#[test]
fn a_section_without_area_is_infinitely_stressed_by_any_force() {
    let mut body = body_with_area(1.0, 0.04);
    body.cross_section_area = 0.0;
    assert_eq!(body.stress([1.0e-9, 0.0, 0.0]), f64::INFINITY);
    assert_eq!(body.stress([0.0; 3]), 0.0);
}

// ── fracture decision ───────────────────────────────────────────────────────

#[test]
fn fracture_decision_follows_the_stress_not_the_mass() {
    // 100 N over 0.05 m² is 2000 Pa (above the 1000 Pa threshold); over 0.2 m²
    // it is 500 Pa (below).
    for mass in [1.0, 100.0] {
        assert_eq!(
            events_after_one_step(mass, 0.05, 100.0, 1000.0),
            1,
            "mass {mass}: 2000 Pa must fracture"
        );
        assert_eq!(
            events_after_one_step(mass, 0.2, 100.0, 1000.0),
            0,
            "mass {mass}: 500 Pa must not fracture"
        );
    }
}

#[test]
fn equal_force_per_area_gives_the_same_decision_for_any_mass() {
    // Both bodies carry 5000 Pa, the first from 50 N on 0.01 m² and the second
    // from 500 N on 0.1 m².
    let light = events_after_one_step(1.0, 0.01, 50.0, 4000.0);
    let heavy = events_after_one_step(100.0, 0.1, 500.0, 4000.0);
    assert_eq!(light, 1);
    assert_eq!(heavy, 1);
    let light_safe = events_after_one_step(1.0, 0.01, 50.0, 6000.0);
    let heavy_safe = events_after_one_step(100.0, 0.1, 500.0, 6000.0);
    assert_eq!(light_safe, 0);
    assert_eq!(heavy_safe, 0);
}

// ── applied force ───────────────────────────────────────────────────────────

fn unbreakable_threshold() -> FractureThreshold {
    FractureThreshold::new(1.0e18, [1.0, 0.0, 0.0], 4)
}

#[test]
fn applied_force_changes_the_velocity() {
    let mut sim = FractureSimulation::new();
    sim.gravity = [0.0; 3];
    sim.add_body(body_with_area(4.0, 0.04));
    // a = F / m = 8 / 4 = 2 m/s^2 for 0.5 s.
    sim.step(0.5, [8.0, 0.0, 0.0], &unbreakable_threshold());
    let body = &sim.bodies[0];
    assert!(
        (body.velocity[0] - 1.0).abs() < 1.0e-12,
        "{:?}",
        body.velocity
    );
    assert_eq!(body.velocity[1], 0.0);
    assert_eq!(body.velocity[2], 0.0);
    // Semi-implicit Euler: the position advances with the updated velocity.
    assert!(
        (body.position[0] - 0.5).abs() < 1.0e-12,
        "{:?}",
        body.position
    );
}

#[test]
fn gravity_is_mass_independent_while_the_applied_acceleration_is_not() {
    let mut sim = FractureSimulation::new();
    sim.add_body(body_with_area(1.0, 0.04));
    sim.add_body(body_with_area(100.0, 0.04));
    sim.step(0.5, [8.0, 0.0, 0.0], &unbreakable_threshold());
    let (light, heavy) = (&sim.bodies[0], &sim.bodies[1]);
    let free_fall = -9.81 * 0.5;
    assert!((light.velocity[1] - free_fall).abs() < 1.0e-12);
    assert!((heavy.velocity[1] - free_fall).abs() < 1.0e-12);
    assert!(
        (light.velocity[0] - 4.0).abs() < 1.0e-12,
        "{:?}",
        light.velocity
    );
    assert!(
        (heavy.velocity[0] - 0.04).abs() < 1.0e-12,
        "{:?}",
        heavy.velocity
    );
}

// ── fragments ───────────────────────────────────────────────────────────────

#[test]
fn fragments_take_the_density_and_extents_of_the_broken_body() {
    let half_extents = [0.5, 0.2, 0.05];
    let density = 1800.0;
    let mut body = AggregateBody::new(
        7,
        [1.0, 2.0, 3.0],
        vec![3.0, 3.0],
        density,
        half_extents,
        0.02,
    )
    .expect("valid aggregate body");
    body.add_joint(BreakableJoint::new(0, 0, 1, 1.0, 1.0e9, 1.0e9));
    let mut sim = FractureSimulation::new();
    sim.add_body(body);
    let threshold = FractureThreshold::new(10.0, [0.0, 1.0, 0.0], 8);

    // 100 N over 0.02 m² is 5000 Pa, well above the 10 Pa threshold.
    let events = sim.step(0.01, [100.0, 0.0, 0.0], &threshold);
    assert_eq!(events, 1);
    assert_eq!(sim.fragment_sets.len(), 1);

    let set = &sim.fragment_sets[0];
    assert_eq!(set.len(), 8);
    assert!((set.total_mass() - 6.0).abs() < 1.0e-12);
    for fragment in &set.fragments {
        assert!(
            (fragment.volume - fragment.mass / density).abs() < 1.0e-15,
            "volume {} must be mass / {density}",
            fragment.volume
        );
        for (coordinate, extent) in fragment.seed.iter().zip(half_extents.iter()) {
            assert!(
                coordinate.abs() <= extent + 1.0e-12,
                "seed {:?} outside the body's extents",
                fragment.seed
            );
        }
    }
    let spread = |axis: usize| {
        set.fragments
            .iter()
            .map(|f| f.seed[axis].abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(
        spread(0) > 0.4,
        "x spread {} should follow the 0.5 m extent",
        spread(0)
    );
    assert!(
        spread(2) <= 0.05 + 1.0e-12,
        "z spread {} must stay thin",
        spread(2)
    );
}

#[test]
fn voronoi_split_fragment_volume_follows_the_supplied_density() {
    let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0; 3]);
    let set = FragmentSet::voronoi_split(5.0, 1250.0, [0.1; 3], 5, plane, 0.0, [0.0; 3]);
    assert_eq!(set.len(), 5);
    for fragment in &set.fragments {
        // Each fragment is 1 kg of 1250 kg/m^3 material.
        assert!(
            (fragment.volume - 8.0e-4).abs() < 1.0e-15,
            "{}",
            fragment.volume
        );
    }
}

// ── construction ────────────────────────────────────────────────────────────

fn build(masses: Vec<f64>, density: f64, extents: [f64; 3], area: f64) -> Result<(), Error> {
    AggregateBody::new(0, [0.0; 3], masses, density, extents, area).map(|_| ())
}

#[test]
fn invalid_masses_are_rejected() {
    let ok = [0.1; 3];
    assert!(matches!(
        build(vec![], DENSITY, ok, 0.04),
        Err(Error::NonPositiveMass { .. })
    ));
    assert!(matches!(
        build(vec![1.0, -2.0], DENSITY, ok, 0.04),
        Err(Error::NonPositiveMass { mass }) if mass == -2.0
    ));
    assert!(matches!(
        build(vec![f64::NAN], DENSITY, ok, 0.04),
        Err(Error::NonPositiveMass { .. })
    ));
    assert!(build(vec![1.0, 2.0], DENSITY, ok, 0.04).is_ok());
}

#[test]
fn invalid_geometry_is_rejected() {
    let ok = [0.1; 3];
    for area in [0.0, -0.04, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(build(vec![1.0], DENSITY, ok, area), Err(Error::General(_))),
            "area {area} must be rejected"
        );
    }
    for density in [0.0, -1.0, f64::NAN] {
        assert!(
            matches!(build(vec![1.0], density, ok, 0.04), Err(Error::General(_))),
            "density {density} must be rejected"
        );
    }
    assert!(matches!(
        build(vec![1.0], DENSITY, [0.1, 0.0, 0.1], 0.04),
        Err(Error::General(_))
    ));
    assert!(matches!(
        build(vec![1.0], DENSITY, [0.1, 0.1, f64::NAN], 0.04),
        Err(Error::General(_))
    ));
}
