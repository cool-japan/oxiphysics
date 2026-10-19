// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the Hertz normal damping sign in `DemContact`.
//!
//! The closing velocity `v_n = (v_a - v_b) . n`, with `n` pointing from `a` to
//! `b`, is positive while two particles approach.  The dashpot has to add to the
//! elastic repulsion while the particles close and take away from it while they
//! separate.  With the opposite sign a collision injects kinetic energy.

use oxiphysics_rigid::granular::{DemContact, DemSimulation, GranularParticle};

const RADIUS: f64 = 0.01;
const DENSITY: f64 = 2500.0;
const KN: f64 = 1.0e6;
const DT: f64 = 1.0e-6;
const STEPS: usize = 8000;

/// A contact law with only the normal Hertz spring and dashpot switched on.
fn normal_only(cn: f64) -> DemContact {
    DemContact {
        kn: KN,
        cn,
        kt: 0.0,
        ct: 0.0,
        mu: 0.0,
    }
}

/// Two equal spheres on the x axis, `a` on the left, approaching head-on.
fn head_on_pair(contact: DemContact, speed_a: f64, speed_b: f64) -> DemSimulation {
    let mut sim = DemSimulation::new(contact, 0.0);
    let mut a = GranularParticle::new([-(RADIUS + 5.0e-4), 0.0, 0.0], RADIUS, DENSITY);
    let mut b = GranularParticle::new([RADIUS + 5.0e-4, 0.0, 0.0], RADIUS, DENSITY);
    a.velocity = [speed_a, 0.0, 0.0];
    b.velocity = [speed_b, 0.0, 0.0];
    sim.add_particle(a);
    sim.add_particle(b);
    sim
}

/// Run the collision to completion and return `(ke_before, ke_after, sim)`.
fn collide(contact: DemContact, speed_a: f64, speed_b: f64) -> (f64, f64, DemSimulation) {
    let mut sim = head_on_pair(contact, speed_a, speed_b);
    let ke_before = sim.total_kinetic_energy();
    for _ in 0..STEPS {
        sim.step(DT);
    }
    let ke_after = sim.total_kinetic_energy();
    (ke_before, ke_after, sim)
}

fn separation(sim: &DemSimulation) -> f64 {
    sim.particles[1].position[0] - sim.particles[0].position[0]
}

fn total_momentum_x(sim: &DemSimulation) -> f64 {
    sim.particles.iter().map(|p| p.mass * p.velocity[0]).sum()
}

#[test]
fn glass_bead_head_on_collision_loses_kinetic_energy() {
    let (before, after, sim) = collide(DemContact::glass_beads(), 1.0, -1.0);
    assert!(before > 0.0);
    assert!(
        after < before,
        "a damped collision must dissipate energy: {before} -> {after}"
    );
    assert!(
        after < 0.5 * before,
        "glass-bead damping should remove most of the energy: {before} -> {after}"
    );
    assert!(
        separation(&sim) > 2.0 * RADIUS,
        "the pair must have separated again"
    );
}

#[test]
fn undamped_hertz_collision_is_elastic() {
    let (before, after, _) = collide(normal_only(0.0), 1.0, -1.0);
    assert!(
        (after - before).abs() < 1.0e-5 * before,
        "cn = 0 must conserve energy: {before} -> {after}"
    );
}

#[test]
fn energy_loss_grows_monotonically_with_the_damping_coefficient() {
    let mut previous = 1.0;
    for cn in [0.5, 2.0, 5.0, 20.0, 50.0] {
        let (before, after, _) = collide(normal_only(cn), 1.0, -1.0);
        let ratio = after / before;
        assert!(
            ratio < 1.0,
            "cn = {cn}: energy ratio {ratio} must be below 1"
        );
        assert!(
            ratio < previous,
            "cn = {cn}: energy ratio {ratio} must fall below the previous {previous}"
        );
        previous = ratio;
    }
}

#[test]
fn damped_collision_conserves_linear_momentum_and_pushes_the_pair_apart() {
    let (_, _, sim) = collide(normal_only(5.0), 1.0, 0.0);
    let mass = sim.particles[0].mass;
    let momentum = total_momentum_x(&sim);
    assert!(
        (momentum - mass).abs() < 1.0e-9 * mass,
        "momentum {momentum} must stay equal to the initial {mass}"
    );
    let (va, vb) = (sim.particles[0].velocity[0], sim.particles[1].velocity[0]);
    assert!(
        vb > va,
        "after the collision the struck sphere leads: va = {va}, vb = {vb}"
    );
    assert!(separation(&sim) > 2.0 * RADIUS);
}

#[test]
fn closing_contact_adds_the_dashpot_to_the_elastic_repulsion() {
    let contact = normal_only(50.0);
    let mut a = GranularParticle::new([0.0, 0.0, 0.0], RADIUS, DENSITY);
    let mut b = GranularParticle::new([0.0199, 0.0, 0.0], RADIUS, DENSITY);
    a.velocity = [0.25, 0.0, 0.0];
    b.velocity = [-0.25, 0.0, 0.0];
    let overlap = a.radius + b.radius - 0.0199;
    let elastic = KN * overlap.powf(1.5);
    let closing_speed = 0.5;

    contact.apply_contact(&mut a, &mut b);

    let expected = elastic + 50.0 * closing_speed;
    assert!(
        (a.force[0] + expected).abs() < 1.0e-9 * expected,
        "force on a is {} but the repulsion should be {expected} (elastic {elastic})",
        a.force[0]
    );
    assert!(
        (a.force[0] + b.force[0]).abs() < 1.0e-9 * expected,
        "Newton's third law"
    );
    assert!(
        expected > elastic,
        "closing must be stiffer than the quasi-static spring"
    );
}

#[test]
fn separating_contact_subtracts_the_dashpot_and_never_pulls() {
    let contact = normal_only(50.0);
    let overlap_gap = 0.0199;

    // Slowly separating: the repulsion is reduced but still pushes apart.
    let mut a = GranularParticle::new([0.0, 0.0, 0.0], RADIUS, DENSITY);
    let mut b = GranularParticle::new([overlap_gap, 0.0, 0.0], RADIUS, DENSITY);
    a.velocity = [-0.005, 0.0, 0.0];
    b.velocity = [0.005, 0.0, 0.0];
    let overlap = a.radius + b.radius - overlap_gap;
    let elastic = KN * overlap.powf(1.5);
    contact.apply_contact(&mut a, &mut b);
    let expected = elastic - 50.0 * 0.01;
    assert!(expected > 0.0, "test set-up: repulsion must survive");
    assert!(
        (a.force[0] + expected).abs() < 1.0e-9 * elastic,
        "force on a is {} but the reduced repulsion should be {expected}",
        a.force[0]
    );

    // Quickly separating: the dashpot would exceed the spring, so the contact
    // force is clamped at zero instead of pulling the particles together.
    let mut a = GranularParticle::new([0.0, 0.0, 0.0], RADIUS, DENSITY);
    let mut b = GranularParticle::new([overlap_gap, 0.0, 0.0], RADIUS, DENSITY);
    a.velocity = [-0.25, 0.0, 0.0];
    b.velocity = [0.25, 0.0, 0.0];
    contact.apply_contact(&mut a, &mut b);
    assert_eq!(a.force, [0.0; 3], "no tensile contact force");
    assert_eq!(b.force, [0.0; 3], "no tensile contact force");
}
