// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for samplers that draw from user-supplied floating-point
//! bounds.
//!
//! `rand` 0.10 panics when asked for a value from an empty range, so a joint
//! pinned at `lower == upper` or a flat search box used to abort the process.
//! Such bounds are legitimate; the sampler has to return the bound itself.

use oxiphysics_rigid::robotic_systems::{DhParam, JointType, KinematicChain, WorkspaceAnalysis};
use oxiphysics_rigid::swarm_dynamics::ParticleSwarmOptimizer;
use std::f64::consts::PI;

/// Planar two-link arm (link lengths 1 and 0.5) with the given joint limits.
fn planar_arm(limits: [(f64, f64); 2]) -> KinematicChain {
    let params = vec![
        DhParam::new(1.0, 0.0, 0.0, 0.0),
        DhParam::new(0.5, 0.0, 0.0, 0.0),
    ];
    let mut chain = KinematicChain::new(
        params,
        vec![JointType::Revolute, JointType::Revolute],
        vec![0.0, 0.0],
    );
    chain.limits = limits.to_vec();
    chain
}

// ── workspace sampling ──────────────────────────────────────────────────────

#[test]
fn joint_pinned_at_its_limit_keeps_the_workspace_on_one_circle() {
    let pinned = 0.3;
    let chain = planar_arm([(pinned, pinned), (-PI, PI)]);
    let analysis = WorkspaceAnalysis::new(200);
    let samples = analysis.sample_workspace(&chain);
    assert_eq!(samples.len(), 200);
    // The first joint cannot move, so the end effector stays on the circle of
    // radius 0.5 around the tip of the first link.
    let tip = [pinned.cos(), pinned.sin()];
    for s in &samples {
        let radius = ((s[0] - tip[0]).powi(2) + (s[1] - tip[1]).powi(2)).sqrt();
        assert!((radius - 0.5).abs() < 1.0e-9, "sample {s:?} off the circle");
        assert!(s[2].abs() < 1.0e-12);
    }
    // The second joint does sweep: the samples are not all the same point.
    assert!(
        samples
            .iter()
            .any(|s| (s[0] - samples[0][0]).abs() > 1.0e-3)
    );
}

#[test]
fn fully_pinned_arm_has_a_single_workspace_point() {
    let chain = planar_arm([(0.2, 0.2), (-0.4, -0.4)]);
    let analysis = WorkspaceAnalysis::new(50);
    let samples = analysis.sample_workspace(&chain);
    assert_eq!(samples.len(), 50);
    for s in &samples {
        assert_eq!(*s, samples[0]);
    }
    // The reach and the bounding box follow from that single point.
    let reach = analysis.max_reach(&chain);
    let norm = (samples[0][0].powi(2) + samples[0][1].powi(2) + samples[0][2].powi(2)).sqrt();
    assert!((reach - norm).abs() < 1.0e-12);
    let (lo, hi) = analysis.bounding_box(&chain);
    assert_eq!(lo, hi);
}

#[test]
fn pinned_and_inverted_limits_do_not_panic_the_other_analyses() {
    // An inverted interval has no interior either; the lower bound is used.
    let chain = planar_arm([(0.1, 0.1), (1.0, -1.0)]);
    let analysis = WorkspaceAnalysis::new(64);
    let singular = analysis.detect_singularities(&chain, 1.0e3);
    assert!(!singular.is_empty());
    for q in &singular {
        assert_eq!(q.len(), 2);
    }
    let fraction = analysis.dexterous_workspace_fraction(&chain, -1.0);
    assert_eq!(fraction, 1.0);
}

// ── particle swarm ──────────────────────────────────────────────────────────

#[test]
fn flat_search_box_places_every_particle_on_the_bound() {
    let swarm = ParticleSwarmOptimizer::new(5, 3, 2.0, 2.0);
    assert_eq!(swarm.particles.len(), 5);
    for particle in &swarm.particles {
        assert_eq!(particle.position, [2.0; 3]);
        assert!(particle.velocity.iter().all(|v| *v == 0.0));
    }
}

#[test]
fn proper_search_box_still_samples_inside_the_bounds() {
    let swarm = ParticleSwarmOptimizer::new(40, 4, -1.0, 3.0);
    let mut spread = 0.0_f64;
    for particle in &swarm.particles {
        for &x in &particle.position {
            assert!((-1.0..3.0).contains(&x), "position {x} outside the box");
            spread = spread.max(x.abs());
        }
        for &v in &particle.velocity {
            assert!(
                v.abs() <= 0.4 + 1.0e-12,
                "velocity {v} beyond 10 % of the box"
            );
        }
    }
    assert!(spread > 0.5, "positions should be spread over the box");
}

#[test]
fn swarm_with_a_flat_box_still_runs() {
    let mut swarm = ParticleSwarmOptimizer::new(6, 2, 1.0, 1.0);
    swarm.run(5, &|x: &[f64]| x.iter().map(|v| v * v).sum::<f64>());
    assert_eq!(swarm.iteration, 5);
    assert!(swarm.global_best_value.is_finite());
}
