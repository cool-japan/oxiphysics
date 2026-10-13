// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the CCD event response.
//!
//! `CcdEvent::body_a` / `body_b` carry user-defined IDs while the response must
//! address the bodies by slice index (`index_a` / `index_b`).  Resolving by ID
//! used to silently do nothing for IDs outside the slice and hit the wrong
//! bodies for permuted IDs.

use oxiphysics_rigid::ccd::{CcdEvent, CcdPipeline, SubstepIntegrator, SweptBody};

fn sphere(x: f64, vx: f64, inv_mass: f64) -> SweptBody {
    SweptBody::new([x, 0.0, 0.0], [vx, 0.0, 0.0], 0.5, inv_mass)
}

fn head_on_pair() -> Vec<SweptBody> {
    vec![sphere(-4.5, 1.0, 1.0), sphere(4.5, -1.0, 1.0)]
}

#[test]
fn events_carry_user_ids_and_slice_indices() {
    let bodies = head_on_pair();
    let events = CcdPipeline::new(10.0).detect_tunneling(&bodies, &[10, 20]);
    assert_eq!(events.len(), 1);
    let ev = &events[0];
    assert_eq!((ev.body_a, ev.body_b), (10, 20));
    assert_eq!((ev.index_a, ev.index_b), (0, 1));
    assert!((ev.toi - 4.0).abs() < 1e-12, "toi = {}", ev.toi);
}

#[test]
fn ids_shorter_than_the_slice_fall_back_to_indices_instead_of_panicking() {
    let bodies = vec![
        sphere(-4.5, 1.0, 1.0),
        sphere(4.5, -1.0, 1.0),
        sphere(40.0, 0.0, 1.0),
    ];
    let events = CcdPipeline::new(10.0).detect_tunneling(&bodies, &[7]);
    assert_eq!(events.len(), 1);
    assert_eq!((events[0].body_a, events[0].body_b), (7, 1));
    assert_eq!((events[0].index_a, events[0].index_b), (0, 1));
    let events = CcdPipeline::new(10.0).detect_tunneling(&bodies, &[]);
    assert_eq!((events[0].body_a, events[0].body_b), (0, 1));
}

#[test]
fn run_with_stats_resolves_head_on_spheres_whose_ids_exceed_the_slice() {
    let mut bodies = head_on_pair();
    let stats = CcdPipeline::new(10.0).run_with_stats(&mut bodies, &[10, 20], 1.0);
    assert_eq!(stats.events_found, 1);
    assert_eq!(stats.impulses_applied, 1);
    assert_eq!(stats.pairs_tested, 1);
    assert!((stats.earliest_toi - 4.0).abs() < 1e-9);
    assert!((stats.time_covered - 10.0).abs() < 1e-6);
    // Elastic head-on bounce: velocities are exchanged ...
    assert!(
        (bodies[0].vel[0] + 1.0).abs() < 1e-12,
        "{:?}",
        bodies[0].vel
    );
    assert!(
        (bodies[1].vel[0] - 1.0).abs() < 1e-12,
        "{:?}",
        bodies[1].vel
    );
    // ... and the bodies keep moving for the rest of the step: they touch at
    // x = -/+0.5 after 4 s and then fly apart for another 6 s.
    assert!((bodies[0].pos[0] + 6.5).abs() < 1e-6, "{:?}", bodies[0].pos);
    assert!((bodies[1].pos[0] - 6.5).abs() < 1e-6, "{:?}", bodies[1].pos);
}

#[test]
fn run_with_stats_hits_the_right_bodies_for_every_id_assignment() {
    // Unequal masses make a mix-up observable: the heavy body (inverse mass
    // 0.5) at x = -4.5 moves at +1, the light body at x = 4.5 moves at -1.
    // Elastic outcome: heavy -> -1/3, light -> +5/3.
    for ids in [vec![0, 1], vec![1, 0], vec![10, 20], vec![5], vec![]] {
        let mut bodies = vec![sphere(-4.5, 1.0, 0.5), sphere(4.5, -1.0, 1.0)];
        let stats = CcdPipeline::new(10.0).run_with_stats(&mut bodies, &ids, 1.0);
        assert_eq!(stats.impulses_applied, 1, "ids {ids:?}");
        assert!(
            (bodies[0].vel[0] + 1.0 / 3.0).abs() < 1e-12,
            "ids {ids:?}: heavy body velocity {:?}",
            bodies[0].vel
        );
        assert!(
            (bodies[1].vel[0] - 5.0 / 3.0).abs() < 1e-12,
            "ids {ids:?}: light body velocity {:?}",
            bodies[1].vel
        );
        let momentum = 2.0 * bodies[0].vel[0] + bodies[1].vel[0];
        assert!((momentum - 1.0).abs() < 1e-12, "momentum {momentum}");
    }
}

#[test]
fn only_the_colliding_pair_is_touched_when_other_bodies_are_present() {
    // The colliding pair is at slice indices 1 and 2; slice index 0 is a
    // bystander that must keep its velocity.
    let mut bodies = vec![
        sphere(100.0, 0.25, 1.0),
        sphere(-4.5, 1.0, 1.0),
        sphere(4.5, -1.0, 1.0),
    ];
    let stats = CcdPipeline::new(10.0).run_with_stats(&mut bodies, &[100, 200, 300], 1.0);
    assert_eq!(stats.pairs_tested, 3);
    assert_eq!(stats.impulses_applied, 1);
    assert_eq!(bodies[0].vel, [0.25, 0.0, 0.0]);
    assert!((bodies[1].vel[0] + 1.0).abs() < 1e-12);
    assert!((bodies[2].vel[0] - 1.0).abs() < 1e-12);
}

#[test]
fn separating_touching_pair_is_not_an_event() {
    // The spheres touch but already fly apart: nothing to resolve and the
    // whole step must still be integrated.
    let mut bodies = vec![sphere(-0.5, -1.0, 1.0), sphere(0.5, 1.0, 1.0)];
    let stats = CcdPipeline::new(2.0).run_with_stats(&mut bodies, &[0, 1], 1.0);
    assert_eq!(stats.events_found, 0);
    assert_eq!(stats.impulses_applied, 0);
    assert!((stats.time_covered - 2.0).abs() < 1e-12);
    assert!((bodies[0].pos[0] + 2.5).abs() < 1e-12);
    assert!((bodies[1].pos[0] - 2.5).abs() < 1e-12);
    assert_eq!(bodies[0].vel, [-1.0, 0.0, 0.0]);
}

#[test]
fn resolve_event_addresses_bodies_by_slice_index() {
    let integrator = SubstepIntegrator::new(1);
    let event = CcdEvent {
        toi: 0.0,
        body_a: 900,
        body_b: 901,
        index_a: 1,
        index_b: 0,
        contact_normal: [1.0, 0.0, 0.0],
        contact_point: [0.0; 3],
    };
    // Slice index 1 is "A" and sits on the +x side moving towards slice index 0.
    let mut bodies = vec![sphere(-0.5, 1.0, 1.0), sphere(0.5, -1.0, 1.0)];
    assert!(integrator.resolve_event(&mut bodies, &event, 1.0));
    assert!((bodies[0].vel[0] + 1.0).abs() < 1e-12);
    assert!((bodies[1].vel[0] - 1.0).abs() < 1e-12);
}

#[test]
fn resolve_event_ignores_invalid_or_separating_events() {
    let integrator = SubstepIntegrator::new(1);
    let base = CcdEvent {
        toi: 0.0,
        body_a: 0,
        body_b: 1,
        index_a: 0,
        index_b: 1,
        contact_normal: [-1.0, 0.0, 0.0],
        contact_point: [0.0; 3],
    };
    // Index outside the slice.
    let mut bodies = head_on_pair();
    let out_of_range = CcdEvent {
        index_b: 5,
        ..base.clone()
    };
    assert!(!integrator.resolve_event(&mut bodies, &out_of_range, 1.0));
    // Same body twice.
    let same = CcdEvent {
        index_b: 0,
        ..base.clone()
    };
    assert!(!integrator.resolve_event(&mut bodies, &same, 1.0));
    // Separating pair (normal points the other way): no impulse.
    let flipped = CcdEvent {
        contact_normal: [1.0, 0.0, 0.0],
        ..base.clone()
    };
    assert!(!integrator.resolve_event(&mut bodies, &flipped, 1.0));
    // Two immovable bodies.
    let mut statics = vec![sphere(-0.5, 1.0, 0.0), sphere(0.5, -1.0, 0.0)];
    assert!(!integrator.resolve_event(&mut statics, &base, 1.0));
    assert_eq!(bodies, head_on_pair());
}
