// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the physics pipeline.
//!
//! This file holds the shared scene helpers and the tests of the pipeline
//! plumbing (statistics, builder, checkpoints, callbacks, profiling). The
//! physics of collision detection, the contact solver, continuous collision
//! detection and the SPH-rigid coupling are tested in the child modules.

mod ccd_tests;
mod collision_tests;
mod coupling_tests;
mod solver_tests;

use super::*;
use oxiphysics_core::{BodyHandle, Transform};
use oxiphysics_geometry::{BoxShape, Shape, Sphere};
use oxiphysics_rigid::{Collider, RigidBody};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

/// 60 Hz step used by most tests.
pub(super) const DT: Real = 1.0 / 60.0;

/// A configuration whose bodies never fall asleep.
pub(super) fn awake_config() -> PhysicsConfig {
    PhysicsConfig {
        time_before_sleep: 1.0e9,
        linear_sleep_threshold: 0.0,
        angular_sleep_threshold: 0.0,
        ..PhysicsConfig::default()
    }
}

/// Bodies and colliders of a test scene.
pub(super) struct Scene {
    pub(super) bodies: RigidBodySet,
    pub(super) colliders: ColliderSet,
}

impl Scene {
    pub(super) fn new() -> Self {
        Self {
            bodies: RigidBodySet::new(),
            colliders: ColliderSet::new(),
        }
    }

    /// Insert a body with one collider.
    pub(super) fn add(
        &mut self,
        body: RigidBody,
        shape: Arc<dyn Shape>,
        friction: Real,
        restitution: Real,
    ) -> (BodyHandle, ColliderHandle) {
        let body_handle = self.bodies.insert(body);
        let collider_handle = self.colliders.insert(
            Collider::new(shape)
                .with_body(body_handle)
                .with_friction(friction)
                .with_restitution(restitution),
        );
        (body_handle, collider_handle)
    }

    /// Insert a static box (a static body with a box collider).
    pub(super) fn static_box(
        &mut self,
        centre: Vec3,
        half_extents: Vec3,
        friction: Real,
        restitution: Real,
    ) -> (BodyHandle, ColliderHandle) {
        let mut body = RigidBody::new_static();
        body.transform = Transform::from_position(centre);
        self.add(
            body,
            Arc::new(BoxShape::new(half_extents)),
            friction,
            restitution,
        )
    }

    /// Insert a dynamic body of the given mass at `position`.
    pub(super) fn dynamic(
        &mut self,
        mass: Real,
        position: Vec3,
        shape: Arc<dyn Shape>,
        friction: Real,
        restitution: Real,
    ) -> (BodyHandle, ColliderHandle) {
        let mut body = RigidBody::new(mass);
        body.transform = Transform::from_position(position);
        self.add(body, shape, friction, restitution)
    }

    /// Advance the scene by one 60 Hz step.
    pub(super) fn step(&mut self, pipeline: &mut PhysicsPipeline) {
        pipeline.step(DT, &mut self.bodies, &self.colliders);
    }

    /// Position of a body.
    pub(super) fn position(&self, body: BodyHandle) -> Vec3 {
        self.bodies
            .get(body)
            .map(|b| b.transform.position)
            .unwrap_or_default()
    }

    /// Linear velocity of a body.
    pub(super) fn velocity(&self, body: BodyHandle) -> Vec3 {
        self.bodies
            .get(body)
            .map(|b| b.velocity)
            .unwrap_or_default()
    }
}

/// A unit sphere shape of the given radius.
pub(super) fn sphere(radius: Real) -> Arc<dyn Shape> {
    Arc::new(Sphere::new(radius))
}

/// A box shape with the given half-extents.
pub(super) fn cuboid(half_x: Real, half_y: Real, half_z: Real) -> Arc<dyn Shape> {
    Arc::new(BoxShape::new(Vec3::new(half_x, half_y, half_z)))
}

#[test]
fn test_pipeline_gravity() {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();

    let mut body = RigidBody::new(1.0);
    body.transform = Transform::from_position(Vec3::new(0.0, 10.0, 0.0));
    body.linear_damping = 0.0;
    let h = bodies.insert(body);

    // Step for 1 second
    for _ in 0..60 {
        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    }

    let body = bodies.get(h).unwrap();
    assert!(
        body.transform.position.y < 6.0,
        "y={}",
        body.transform.position.y
    );
    assert!(body.velocity.y < -5.0, "vy={}", body.velocity.y);
}

#[test]
fn test_pipeline_sphere_on_floor() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    let (sphere_body, _) = scene.dynamic(1.0, Vec3::new(0.0, 5.0, 0.0), sphere(0.5), 0.5, 0.0);
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );

    for _ in 0..300 {
        scene.step(&mut pipeline);
    }

    // The sphere comes to rest on the floor: centre one radius above y = 0.
    let p = scene.position(sphere_body);
    assert!(
        (p.y - 0.5).abs() < 0.02,
        "sphere should rest on the floor at y = 0.5, got y = {}",
        p.y
    );
    assert!(scene.velocity(sphere_body).norm() < 1e-3);
}

#[test]
fn test_pipeline_no_colliders() {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();

    bodies.insert(RigidBody::new(1.0));
    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
}

#[test]
fn test_pipeline_stats_step_count() {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    assert_eq!(pipeline.stats.step_count, 0);
    for i in 1..=5u64 {
        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        assert_eq!(pipeline.stats.step_count, i);
    }
}

#[test]
fn test_pipeline_stats_timing_accumulates() {
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    scene.dynamic(1.0, Vec3::new(0.0, 0.49, 0.0), sphere(0.5), 0.5, 0.0);
    let mut pipeline = PhysicsPipeline::with_config(awake_config());

    for _ in 0..50 {
        scene.step(&mut pipeline);
    }

    // The collision stages of a CCD-enabled step did real work every step.
    let timings = &pipeline.stats.timings;
    assert!(timings.broadphase > Duration::ZERO);
    assert!(timings.narrowphase > Duration::ZERO);
    assert!(timings.solver > Duration::ZERO);
    assert_eq!(timings.coupling, Duration::ZERO, "no coupling installed");
    assert_eq!(
        timings.total(),
        timings.forces
            + timings.coupling
            + timings.broadphase
            + timings.narrowphase
            + timings.solver
            + timings.integration
            + timings.sleep_check
    );
    for stage in ["broadphase", "narrowphase", "solver"] {
        assert!(
            pipeline.stats.avg_stage_time(stage) > Duration::ZERO,
            "stage {stage} has no timing"
        );
    }
    assert_eq!(pipeline.stats.avg_stage_time("unknown"), Duration::ZERO);

    pipeline.stats.timings.reset();
    assert_eq!(pipeline.stats.timings.total(), Duration::ZERO);
}

#[test]
fn test_pipeline_stats_count_pairs_and_manifolds_without_ccd() {
    let config = PhysicsConfig {
        ccd_enabled: false,
        ..awake_config()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    scene.dynamic(1.0, Vec3::new(0.0, 0.49, 0.0), sphere(0.5), 0.5, 0.0);
    scene.dynamic(1.0, Vec3::new(5.0, 3.0, 0.0), sphere(0.5), 0.5, 0.0);

    scene.step(&mut pipeline);
    assert_eq!(pipeline.stats.last_broadphase_pairs, 1);
    assert_eq!(pipeline.stats.last_contact_manifolds, 1);
    assert_eq!(pipeline.stats.last_ccd_events, 0);
}

#[test]
fn test_pipeline_stats_reset() {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    assert_eq!(pipeline.stats.step_count, 1);

    pipeline.reset_stats();
    assert_eq!(pipeline.stats.step_count, 0);
}

#[test]
fn test_pipeline_builder() {
    let pipeline = PhysicsPipeline::builder()
        .gravity(Vec3::new(0.0, -9.81, 0.0))
        .ccd_enabled(false)
        .solver_iterations(8_u32)
        .max_ccd_substeps(3)
        .profiling(true)
        .finish();

    assert!((pipeline.config.gravity.y + 9.81).abs() < 1e-12);
    assert!(!pipeline.config.ccd_enabled);
    assert_eq!(pipeline.config.solver_iterations, 8_u32);
    assert_eq!(pipeline.max_ccd_substeps, 3);
    assert!(pipeline.profiling_enabled);
    assert!(pipeline.sph_rigid_coupling().is_none());
}

#[test]
fn test_pipeline_default_ccd_budget() {
    let pipeline = PhysicsPipeline::new();
    assert_eq!(pipeline.max_ccd_substeps, 8);
    assert!(pipeline.config.ccd_enabled);
}

#[test]
fn test_pipeline_types_stay_send_and_sync() {
    // The pipeline is handed to worker threads and stored in shared state by
    // callers; its (private) solver and collision state must not take that
    // away.
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<PhysicsPipeline>();
    assert_send_sync::<PipelineBuilder>();
    assert_send_sync::<PipelineCheckpoint>();
    assert_send_sync::<PipelineStats>();
    assert_send_sync::<WarmStartCache>();
    assert_send_sync::<SphRigidCouplingParams>();
}

#[test]
fn test_pipeline_config_builder() {
    let cfg = PipelineConfigBuilder::new()
        .gravity(Vec3::new(0.0, -9.81, 0.0))
        .solver_iterations(12_u32)
        .ccd_enabled(true)
        .build();

    assert!((cfg.gravity.y + 9.81).abs() < 1e-12);
    assert_eq!(cfg.solver_iterations, 12_u32);
    assert!(cfg.ccd_enabled);
}

#[test]
fn test_pipeline_checkpoint_restore() {
    let mut pipeline = PhysicsPipeline::with_config(awake_config());
    let mut scene = Scene::new();
    scene.static_box(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(10.0, 0.5, 10.0),
        0.5,
        0.0,
    );
    scene.dynamic(
        1.0,
        Vec3::new(0.0, 0.49, 0.0),
        cuboid(0.5, 0.5, 0.5),
        0.5,
        0.0,
    );

    // Run a few steps and take a checkpoint.
    for _ in 0..5 {
        scene.step(&mut pipeline);
    }
    let checkpoint = pipeline.checkpoint();
    let time_at_checkpoint = checkpoint.time;
    let cached_at_checkpoint = checkpoint.contact_cache.contact_count();
    assert!(
        cached_at_checkpoint > 0,
        "the resting box has warm-start data"
    );

    // Run more steps, then lose the cache.
    for _ in 0..10 {
        scene.step(&mut pipeline);
    }
    assert!(pipeline.time > time_at_checkpoint);
    pipeline.contact_cache.clear();
    assert!(pipeline.warm_start_cache().is_empty());

    // Restore: time and the warm-start cache come back.
    pipeline.restore_checkpoint(checkpoint);
    assert!((pipeline.time - time_at_checkpoint).abs() < 1e-14);
    assert_eq!(
        pipeline.warm_start_cache().contact_count(),
        cached_at_checkpoint
    );
}

#[test]
fn test_stage_callback_fires() {
    let counter = Arc::new(AtomicUsize::new(0));
    let counter_clone = Arc::clone(&counter);

    // Without CCD a step reports forces, broadphase, narrowphase, solver,
    // integration and sleep_check.
    let config = PhysicsConfig {
        ccd_enabled: false,
        ..PhysicsConfig::default()
    };
    let mut pipeline = PhysicsPipeline::with_config(config);
    pipeline.set_stage_callback(Box::new(move |_stage, _dur| {
        counter_clone.fetch_add(1, Ordering::Relaxed);
    }));

    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    assert_eq!(counter.load(Ordering::Relaxed), 6);

    pipeline.clear_stage_callback();
    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    assert_eq!(counter.load(Ordering::Relaxed), 6, "callback was removed");
}

#[test]
fn test_stage_callback_reports_the_ccd_search_as_well() {
    let seen = Arc::new(std::sync::Mutex::new(Vec::<String>::new()));
    let seen_clone = Arc::clone(&seen);

    let mut pipeline = PhysicsPipeline::new();
    pipeline.set_stage_callback(Box::new(move |stage, _dur| {
        if let Ok(mut stages) = seen_clone.lock() {
            stages.push(stage.to_owned());
        }
    }));
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));
    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);

    // The contact pass first, then the swept CCD search, then the integration.
    let stages = seen.lock().map(|s| s.clone()).unwrap_or_default();
    assert_eq!(
        stages,
        vec![
            "forces",
            "broadphase",
            "narrowphase",
            "solver",
            "broadphase",
            "narrowphase",
            "integration",
            "sleep_check",
        ]
    );
}

#[test]
fn test_profiling_records_stages() {
    let mut pipeline = PhysicsPipeline::new();
    pipeline.enable_profiling(Some(128));

    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);

    let profiles = pipeline.drain_profiles();
    assert!(
        !profiles.is_empty(),
        "profiling should record at least one stage"
    );
    // All profiles should be for step 0 (step_count was 0 when recorded, incremented after).
    assert!(profiles.iter().all(|p| p.step == 0));
    assert!(
        pipeline.drain_profiles().is_empty(),
        "draining empties the history"
    );
}

#[test]
fn test_profiling_history_is_bounded() {
    let mut pipeline = PhysicsPipeline::new();
    pipeline.enable_profiling(Some(4));
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    for _ in 0..10 {
        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    }
    let profiles = pipeline.drain_profiles();
    assert_eq!(profiles.len(), 4);
    // The ring buffer keeps the most recent records.
    assert!(profiles.iter().all(|p| p.step == 9));
}

#[test]
fn test_profiling_with_zero_limit_keeps_nothing_and_does_not_panic() {
    let mut pipeline = PhysicsPipeline::new();
    pipeline.enable_profiling(Some(0));
    let mut bodies = RigidBodySet::new();
    let colliders = ColliderSet::new();
    bodies.insert(RigidBody::new(1.0));

    for _ in 0..3 {
        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
    }
    assert!(pipeline.drain_profiles().is_empty());
    assert!(
        pipeline.stats.timings.forces > Duration::ZERO,
        "stage timings are accumulated independently of the history"
    );
}

#[test]
fn test_total_kinetic_energy_zero_static() {
    let mut bodies = RigidBodySet::new();
    let mut floor = RigidBody::new_static();
    floor.transform = Transform::from_position(Vec3::zeros());
    bodies.insert(floor);

    let ek = total_kinetic_energy(&bodies);
    assert!(
        ek.abs() < 1e-14,
        "static body should contribute zero kinetic energy"
    );
}

#[test]
fn test_total_kinetic_energy_moving_body() {
    let mut bodies = RigidBodySet::new();
    let mut body = RigidBody::new(2.0);
    body.velocity = Vec3::new(3.0, 4.0, 0.0); // speed = 5, KE = 0.5*2*25 = 25
    bodies.insert(body);

    let ek = total_kinetic_energy(&bodies);
    assert!((ek - 25.0).abs() < 1e-9, "KE should be 25 J, got {}", ek);
}

#[test]
fn test_count_awake_bodies() {
    let mut bodies = RigidBodySet::new();
    bodies.insert(RigidBody::new(1.0)); // dynamic, awake by default
    let mut floor = RigidBody::new_static();
    floor.transform = Transform::from_position(Vec3::zeros());
    bodies.insert(floor);
    let mut asleep = RigidBody::new(1.0);
    asleep.state = BodyState::Sleeping;
    bodies.insert(asleep);

    let awake = count_awake_bodies(&bodies);
    assert_eq!(awake, 1, "only the awake dynamic body is counted");
}
