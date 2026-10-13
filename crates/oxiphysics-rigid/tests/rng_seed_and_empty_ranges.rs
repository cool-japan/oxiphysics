// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the seeded fracture generators and for random sampling
//! over degenerate (empty) floating-point ranges.
//!
//! `rand` 0.10 panics when asked for a value from an empty range, so zero
//! scatter, a flat bounding box or a pinned joint used to abort the process.

use oxiphysics_rigid::crowd_simulation::spawn_pedestrians_in_rect;
use oxiphysics_rigid::destruction_rigid::{
    DebrisConfig, FractureMaterial, VoronoiCell, execute_voronoi_destruction, generate_planar_cuts,
    generate_voronoi_sites, spawn_debris,
};
use oxiphysics_rigid::neural_control::{ModelPredictiveControl, MotionPlanning};
use oxiphysics_rigid::swarm_physics::{BoidForces, SwarmSimulator};

fn unit_cell(index: usize) -> VoronoiCell {
    let x = index as f64;
    VoronoiCell {
        site_index: index,
        vertices: vec![
            [x, 0.0, 0.0],
            [x + 1.0, 0.0, 0.0],
            [x, 1.0, 0.0],
            [x, 0.0, 1.0],
        ],
        volume: 1.0,
        centroid: [x + 0.5, 0.5, 0.5],
    }
}

fn cells(count: usize) -> Vec<VoronoiCell> {
    (0..count).map(unit_cell).collect()
}

// ── seeded generators ───────────────────────────────────────────────────────

#[test]
fn voronoi_sites_are_reproducible_for_a_seed() {
    let a = generate_voronoi_sites([0.0; 3], [2.0, 3.0, 4.0], 24, 1234);
    let b = generate_voronoi_sites([0.0; 3], [2.0, 3.0, 4.0], 24, 1234);
    assert_eq!(a.len(), 24);
    for (sa, sb) in a.iter().zip(b.iter()) {
        assert_eq!(sa.position, sb.position, "same seed must give same sites");
    }
}

#[test]
fn voronoi_sites_differ_between_seeds_and_stay_inside_the_box() {
    let min = [-1.0, 0.0, 2.0];
    let max = [1.0, 5.0, 3.0];
    let a = generate_voronoi_sites(min, max, 32, 1);
    let b = generate_voronoi_sites(min, max, 32, 2);
    let differing = a
        .iter()
        .zip(b.iter())
        .filter(|(sa, sb)| sa.position != sb.position)
        .count();
    assert!(
        differing >= 30,
        "seeds 1 and 2 agree on {} sites",
        32 - differing
    );
    for site in a.iter().chain(b.iter()) {
        for k in 0..3 {
            assert!(
                site.position[k] >= min[k] && site.position[k] < max[k],
                "site {:?} outside the box",
                site.position
            );
        }
    }
}

#[test]
fn planar_cuts_are_reproducible_for_a_seed() {
    let centre = [0.5, -1.0, 2.0];
    let a = generate_planar_cuts(centre, 12, 77);
    let b = generate_planar_cuts(centre, 12, 77);
    let c = generate_planar_cuts(centre, 12, 78);
    assert_eq!(a.len(), 12);
    for (pa, pb) in a.iter().zip(b.iter()) {
        assert_eq!(pa.normal, pb.normal);
        assert_eq!(pa.point, centre);
    }
    let differing = a
        .iter()
        .zip(c.iter())
        .filter(|(pa, pc)| pa.normal != pc.normal)
        .count();
    assert!(
        differing >= 10,
        "different seeds should give different planes"
    );
    for plane in &a {
        let n = plane.normal;
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-12,
            "plane normal must be a unit vector"
        );
    }
}

#[test]
fn debris_velocities_are_reproducible_for_a_seed() {
    let material = FractureMaterial::default();
    let config = DebrisConfig {
        seed: 4242,
        ..DebrisConfig::default()
    };
    let a = spawn_debris(&cells(8), &material, [1.0, 0.0, 0.0], &config);
    let b = spawn_debris(&cells(8), &material, [1.0, 0.0, 0.0], &config);
    assert_eq!(a.len(), 8);
    for (pa, pb) in a.iter().zip(b.iter()) {
        assert_eq!(pa.fragment.velocity, pb.fragment.velocity);
        assert_eq!(pa.fragment.angular_velocity, pb.fragment.angular_velocity);
    }
}

#[test]
fn debris_velocities_differ_between_seeds() {
    let material = FractureMaterial::default();
    let low = DebrisConfig {
        seed: 1,
        ..DebrisConfig::default()
    };
    let high = DebrisConfig {
        seed: 2,
        ..DebrisConfig::default()
    };
    let a = spawn_debris(&cells(8), &material, [0.0; 3], &low);
    let b = spawn_debris(&cells(8), &material, [0.0; 3], &high);
    let differing = a
        .iter()
        .zip(b.iter())
        .filter(|(pa, pb)| pa.fragment.velocity != pb.fragment.velocity)
        .count();
    assert_eq!(
        differing, 8,
        "every fragment should get a different scatter"
    );
}

#[test]
fn voronoi_destruction_is_reproducible_per_body() {
    let material = FractureMaterial::concrete();
    let run = |body_id: usize| {
        execute_voronoi_destruction(
            body_id, [0.0; 3], [1.0; 3], [0.5; 3], 500.0, [0.0; 3], &material, 6, 10,
        )
    };
    let first = run(3);
    let second = run(3);
    let other = run(4);
    let volumes = |event: &oxiphysics_rigid::destruction_rigid::DestructionEvent| {
        event
            .fragments
            .iter()
            .map(|f| f.volume)
            .collect::<Vec<f64>>()
    };
    assert_eq!(volumes(&first), volumes(&second));
    assert_ne!(volumes(&first), volumes(&other));
}

// ── empty-range robustness ──────────────────────────────────────────────────

#[test]
fn zero_velocity_scatter_does_not_panic_and_keeps_the_parent_velocity() {
    let material = FractureMaterial::default();
    let config = DebrisConfig {
        velocity_scatter: 0.0,
        angular_scatter: 0.0,
        ..DebrisConfig::default()
    };
    let parent = [0.5, -0.25, 2.0];
    let debris = spawn_debris(&cells(5), &material, parent, &config);
    assert_eq!(debris.len(), 5);
    for particle in &debris {
        assert_eq!(particle.fragment.velocity, parent);
        assert_eq!(particle.fragment.angular_velocity, [0.0; 3]);
    }
}

#[test]
fn flat_aabb_places_sites_on_the_flat_plane() {
    let sites = generate_voronoi_sites([0.0, 1.5, 0.0], [2.0, 1.5, 2.0], 16, 9);
    assert_eq!(sites.len(), 16);
    for site in &sites {
        assert_eq!(site.position[1], 1.5, "flat axis must keep its coordinate");
        assert!((0.0..2.0).contains(&site.position[0]));
        assert!((0.0..2.0).contains(&site.position[2]));
    }
    let point = generate_voronoi_sites([3.0; 3], [3.0; 3], 4, 0);
    assert!(point.iter().all(|s| s.position == [3.0; 3]));
}

#[test]
fn zero_sized_swarm_box_does_not_panic() {
    let swarm = SwarmSimulator::new(4, 0.0, BoidForces::default_weights());
    assert_eq!(swarm.agent_count(), 4);
    for agent in &swarm.agents {
        assert_eq!(agent.pos, [0.0; 3]);
    }
    let spread = SwarmSimulator::new(16, 5.0, BoidForces::default_weights());
    for agent in &spread.agents {
        assert!(agent.pos[0].abs() < 5.0 && agent.pos[2].abs() < 5.0);
    }
}

#[test]
fn flat_pedestrian_spawn_region_does_not_panic() {
    let crowd = spawn_pedestrians_in_rect(6, 1.0, 1.0, 0.0, 4.0, [10.0, 0.0], 70.0);
    assert_eq!(crowd.len(), 6);
    for pedestrian in &crowd {
        assert_eq!(pedestrian.position[0], 1.0);
        assert!((0.0..4.0).contains(&pedestrian.position[1]));
    }
    let point = spawn_pedestrians_in_rect(3, 2.0, 2.0, -1.0, -1.0, [5.0, 5.0], 70.0);
    assert!(point.iter().all(|p| p.position == [2.0, -1.0]));
}

#[test]
fn pinned_joint_bounds_do_not_panic_in_random_config() {
    let planner = MotionPlanning::new(
        3,
        vec![0.0, -1.0, 0.5],
        vec![0.0, 1.0, 0.5],
        0.1,
        10,
        0.05,
        false,
    );
    for _ in 0..50 {
        let q = planner.random_config();
        assert_eq!(q.len(), 3);
        assert_eq!(q[0], 0.0, "joint pinned at 0");
        assert!((-1.0..1.0).contains(&q[1]));
        assert_eq!(q[2], 0.5, "joint pinned at 0.5");
    }
}

#[test]
fn zero_action_bound_does_not_panic_in_mpc_planning() {
    let mpc = ModelPredictiveControl::new(2, 2, &[4], 3, 6, 0.95, 0.0);
    let cost = |state: &[f64]| state.iter().map(|x| x * x).sum::<f64>();
    let first = mpc.plan(&[0.5, -0.5], cost);
    assert_eq!(first, vec![0.0, 0.0]);
    let refined = mpc.cem_plan(&[0.5, -0.5], cost, 2);
    assert_eq!(refined, vec![0.0, 0.0]);
}
