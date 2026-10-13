// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Full physics simulation pipeline.
//!
//! Orchestrates the complete rigid body simulation. One
//! [`PhysicsPipeline::step`](crate::pipeline::PhysicsPipeline::step) runs these
//! stages:
//!
//! 1. **forces** – integrate gravity / external forces and damping into the
//!    velocities
//! 2. **coupling** – optional SPH-rigid buoyancy
//!    ([`SphRigidCouplingParams`](crate::pipeline::SphRigidCouplingParams))
//! 3. **broadphase** – rotation-aware world AABBs and candidate pairs
//! 4. **narrowphase** – contact manifolds from the narrow-phase dispatcher
//! 5. **solver** – sequential-impulse contact solve with accumulated impulses,
//!    warm starting and split-impulse position correction
//! 6. **integration** – integrate velocities into positions and orientations;
//!    with CCD enabled the fast pairs are handled first by time-of-impact
//!    sub-steps (see below), which repeat the broadphase and narrowphase
//!    stages in swept form
//! 7. **sleep_check** – put slow bodies to sleep
//!
//! Sensor colliders never receive contact forces; the pairs they overlap with
//! are reported through
//! [`PhysicsPipeline::last_sensor_overlaps`](crate::pipeline::PhysicsPipeline::last_sensor_overlaps).
//!
//! With continuous collision detection enabled the contacts are solved first,
//! exactly as in a step without CCD, so resting contacts are handled in every
//! frame whether or not a fast body produces an event, and the sub-steps move
//! the bodies with their *solved* velocities. The integration then resolves up
//! to [`PhysicsPipeline::max_ccd_substeps`](crate::pipeline::PhysicsPipeline::max_ccd_substeps)
//! times of impact of fast pairs, each with its own sub-step, and integrates
//! what is left of the step.
//!
//! Extended with:
//! - Pipeline statistics (step count, timing per stage)
//! - Pipeline configuration builder
//! - SPH-rigid buoyancy coupling applied by the pipeline
//! - Pipeline checkpointing and replay
//! - Debug callback hooks
//! - Per-stage performance profiling

use oxiphysics_collision::NarrowPhaseDispatcher;
use oxiphysics_collision::broadphase::SweepAndPrune;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_core::{ColliderHandle, PhysicsConfig};
use oxiphysics_rigid::{BodyState, ColliderSet, RigidBodySet};
use std::collections::VecDeque;
use std::time::{Duration, Instant};

mod ccd;
mod collision;
mod coupling;
mod solver;
#[cfg(test)]
mod tests;

use ccd::follows_target;
use collision::{ColliderEntry, collect_entries};
pub use coupling::{SphRigidCouplingParams, apply_buoyancy_impulse};
pub use solver::WarmStartCache;
use solver::{ContactPair, SolverSettings, solve_contacts};

// ---------------------------------------------------------------------------
// Pipeline statistics
// ---------------------------------------------------------------------------

/// Per-stage timing accumulated over multiple steps.
#[derive(Debug, Clone, Default)]
pub struct StageTimings {
    /// Accumulated wall-clock time spent in the force-integration stage.
    pub forces: Duration,
    /// Accumulated wall-clock time spent applying SPH-rigid coupling forces.
    pub coupling: Duration,
    /// Accumulated wall-clock time spent in the broadphase stage.
    pub broadphase: Duration,
    /// Accumulated wall-clock time spent in the narrowphase / CCD stage.
    pub narrowphase: Duration,
    /// Accumulated wall-clock time spent solving constraints / contacts.
    pub solver: Duration,
    /// Accumulated wall-clock time spent integrating positions.
    pub integration: Duration,
    /// Accumulated wall-clock time spent on the sleep-check stage.
    pub sleep_check: Duration,
}

impl StageTimings {
    /// Reset all accumulated timings to zero.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Total accumulated time across all stages.
    pub fn total(&self) -> Duration {
        self.forces
            + self.coupling
            + self.broadphase
            + self.narrowphase
            + self.solver
            + self.integration
            + self.sleep_check
    }
}

/// Rolling statistics collected by [`PhysicsPipeline`].
#[derive(Debug, Clone, Default)]
pub struct PipelineStats {
    /// Number of simulation steps executed so far.
    pub step_count: u64,
    /// Number of broadphase collision pairs found in the last step.
    pub last_broadphase_pairs: usize,
    /// Number of contact manifolds generated in the last step.
    pub last_contact_manifolds: usize,
    /// Number of time-of-impact events resolved by continuous collision
    /// detection in the last step.
    pub last_ccd_events: usize,
    /// Accumulated per-stage timings (wall-clock, not simulated time).
    pub timings: StageTimings,
}

impl PipelineStats {
    /// Average time (wall-clock) per stage per step, or zero if no steps taken.
    pub fn avg_stage_time(&self, stage: &str) -> Duration {
        if self.step_count == 0 {
            return Duration::ZERO;
        }
        let n = self.step_count as u32;
        match stage {
            "forces" => self.timings.forces / n,
            "coupling" => self.timings.coupling / n,
            "broadphase" => self.timings.broadphase / n,
            "narrowphase" => self.timings.narrowphase / n,
            "solver" => self.timings.solver / n,
            "integration" => self.timings.integration / n,
            "sleep_check" => self.timings.sleep_check / n,
            _ => Duration::ZERO,
        }
    }
}

// ---------------------------------------------------------------------------
// Debug callback hook type
// ---------------------------------------------------------------------------

/// Callback invoked after each simulation stage.
///
/// Receives a stage name and the elapsed wall-clock duration for that stage.
pub type StageCallback = Box<dyn Fn(&str, Duration) + Send + Sync>;

// ---------------------------------------------------------------------------
// Pipeline configuration builder
// ---------------------------------------------------------------------------

/// Builder for [`PhysicsConfig`] with a fluent API.
///
/// # Example
/// ```
/// use oxiphysics::core::math::Vec3;
/// use oxiphysics::pipeline::PipelineConfigBuilder;
///
/// let config = PipelineConfigBuilder::new()
///     .gravity(Vec3::new(0.0, -9.81, 0.0))
///     .solver_iterations(10)
///     .ccd_enabled(true)
///     .build();
/// assert_eq!(config.solver_iterations, 10);
/// assert!(config.ccd_enabled);
/// ```
#[derive(Debug, Clone, Default)]
pub struct PipelineConfigBuilder {
    config: PhysicsConfig,
}

impl PipelineConfigBuilder {
    /// Start with default configuration.
    pub fn new() -> Self {
        Self {
            config: PhysicsConfig::default(),
        }
    }

    /// Set the gravity vector.
    pub fn gravity(mut self, g: Vec3) -> Self {
        self.config.gravity = g;
        self
    }

    /// Set the number of sequential-impulse solver iterations per step.
    pub fn solver_iterations(mut self, n: u32) -> Self {
        self.config.solver_iterations = n;
        self
    }

    /// Enable or disable continuous collision detection (CCD).
    pub fn ccd_enabled(mut self, enabled: bool) -> Self {
        self.config.ccd_enabled = enabled;
        self
    }

    /// Set the linear velocity threshold below which a body may sleep.
    pub fn linear_sleep_threshold(mut self, v: Real) -> Self {
        self.config.linear_sleep_threshold = v;
        self
    }

    /// Set the angular velocity threshold below which a body may sleep.
    pub fn angular_sleep_threshold(mut self, v: Real) -> Self {
        self.config.angular_sleep_threshold = v;
        self
    }

    /// Set the time (seconds of simulated time) a body must be slow before sleeping.
    pub fn time_before_sleep(mut self, t: Real) -> Self {
        self.config.time_before_sleep = t;
        self
    }

    /// Build the [`PhysicsConfig`].
    pub fn build(self) -> PhysicsConfig {
        self.config
    }
}

// ---------------------------------------------------------------------------
// Checkpoint / replay
// ---------------------------------------------------------------------------

/// A snapshot of a pipeline's mutable state for checkpointing and replay.
///
/// This stores everything needed to restore the pipeline to a known-good
/// state: the simulation time and the contact warmstarting cache.
/// Body and collider state must be checkpointed separately by the caller.
#[derive(Debug, Clone)]
pub struct PipelineCheckpoint {
    /// Simulation time at the point this checkpoint was taken.
    pub time: Real,
    /// Snapshot of the warmstarting contact cache.
    pub contact_cache: WarmStartCache,
}

// ---------------------------------------------------------------------------
// Per-stage profiling record
// ---------------------------------------------------------------------------

/// A single per-step per-stage profiling record.
#[derive(Debug, Clone)]
pub struct StageProfile {
    /// Step index when this record was taken.
    pub step: u64,
    /// Simulated time at the start of the step.
    pub sim_time: Real,
    /// Wall-clock duration of the stage.
    pub duration: Duration,
    /// Name of the stage.
    pub stage: &'static str,
}

// ---------------------------------------------------------------------------
// PhysicsPipeline
// ---------------------------------------------------------------------------

/// A full physics simulation pipeline.
pub struct PhysicsPipeline {
    /// Configuration.
    pub config: PhysicsConfig,
    /// Simulation time.
    pub time: Real,
    /// Maximum number of time-of-impact events a CCD step resolves; what is
    /// left of the step after that is integrated without further impact
    /// handling.
    pub max_ccd_substeps: usize,
    /// Broadphase algorithm.
    broadphase: SweepAndPrune,
    /// Narrow-phase algorithm registry (B-to-A contact normals).
    narrowphase: NarrowPhaseDispatcher,
    /// Warmstarting cache: accumulated impulses per contact point from the previous step.
    contact_cache: WarmStartCache,
    /// Sensor overlaps found in the last step.
    sensor_overlaps: Vec<(ColliderHandle, ColliderHandle)>,
    /// Rolling statistics collected during simulation.
    pub stats: PipelineStats,
    /// Optional debug callback invoked after each stage.
    stage_callback: Option<StageCallback>,
    /// Per-stage profile history (bounded ring-buffer of the last N records).
    profile_history: VecDeque<StageProfile>,
    /// Maximum number of stage profile records to retain.
    profile_history_limit: usize,
    /// Whether per-stage profiling is enabled.
    pub profiling_enabled: bool,
    /// Active SPH-rigid coupling parameters, if any.
    sph_rigid_coupling: Option<SphRigidCouplingParams>,
}

impl PhysicsPipeline {
    /// Create a new pipeline with default config.
    pub fn new() -> Self {
        Self::with_config(PhysicsConfig::default())
    }

    /// Create with custom config.
    pub fn with_config(config: PhysicsConfig) -> Self {
        Self {
            config,
            time: 0.0,
            max_ccd_substeps: ccd::DEFAULT_MAX_CCD_SUBSTEPS,
            broadphase: SweepAndPrune::default(),
            narrowphase: NarrowPhaseDispatcher::default(),
            contact_cache: WarmStartCache::new(),
            sensor_overlaps: Vec::new(),
            stats: PipelineStats::default(),
            stage_callback: None,
            profile_history: VecDeque::new(),
            profile_history_limit: 1024,
            profiling_enabled: false,
            sph_rigid_coupling: None,
        }
    }

    /// Use a [`PipelineConfigBuilder`] to construct and configure a pipeline.
    ///
    /// # Example
    /// ```
    /// use oxiphysics::core::math::Vec3;
    /// use oxiphysics::pipeline::PhysicsPipeline;
    ///
    /// let pipeline = PhysicsPipeline::builder()
    ///     .gravity(Vec3::new(0.0, -9.81, 0.0))
    ///     .ccd_enabled(true)
    ///     .finish();
    /// assert!(pipeline.config.ccd_enabled);
    /// ```
    pub fn builder() -> PipelineBuilder {
        PipelineBuilder::new()
    }

    /// Register a debug callback that is invoked after each simulation stage.
    ///
    /// The callback receives the stage name and the elapsed wall-clock duration.
    pub fn set_stage_callback(&mut self, cb: StageCallback) {
        self.stage_callback = Some(cb);
    }

    /// Remove any previously registered stage callback.
    pub fn clear_stage_callback(&mut self) {
        self.stage_callback = None;
    }

    /// Enable per-stage profiling with an optional ring-buffer size limit.
    ///
    /// When `limit` is `None` the default of 1024 records is used. A limit of
    /// zero keeps no records (the per-stage timings in [`PipelineStats`] are
    /// still accumulated).
    pub fn enable_profiling(&mut self, limit: Option<usize>) {
        self.profiling_enabled = true;
        if let Some(l) = limit {
            self.profile_history_limit = l;
        }
        while self.profile_history.len() > self.profile_history_limit {
            self.profile_history.pop_front();
        }
    }

    /// Disable per-stage profiling.
    pub fn disable_profiling(&mut self) {
        self.profiling_enabled = false;
    }

    /// Drain and return all accumulated profile records, leaving an empty history.
    pub fn drain_profiles(&mut self) -> Vec<StageProfile> {
        self.profile_history.drain(..).collect()
    }

    /// Configure SPH-rigid coupling.
    ///
    /// From the next [`step`](Self::step) on, the buoyancy of the fluid region
    /// described by `params` is applied to the rigid bodies that dip into it
    /// (see [`SphRigidCouplingParams`]).
    ///
    /// # Example
    /// ```
    /// use oxiphysics::core::math::Vec3;
    /// use oxiphysics::pipeline::{PhysicsPipeline, SphRigidCouplingParams};
    ///
    /// let mut pipeline = PhysicsPipeline::new();
    /// pipeline.set_sph_rigid_coupling(SphRigidCouplingParams {
    ///     fluid_density: 1000.0,
    ///     coupling_alpha: 1.0,
    ///     region_half_extents: Vec3::new(2.0, 1.0, 2.0),
    ///     region_centre: Vec3::zeros(),
    /// });
    /// assert!(pipeline.sph_rigid_coupling().is_some());
    ///
    /// pipeline.clear_sph_rigid_coupling();
    /// assert!(pipeline.sph_rigid_coupling().is_none());
    /// ```
    pub fn set_sph_rigid_coupling(&mut self, params: SphRigidCouplingParams) {
        self.sph_rigid_coupling = Some(params);
    }

    /// Disable SPH-rigid coupling.
    pub fn clear_sph_rigid_coupling(&mut self) {
        self.sph_rigid_coupling = None;
    }

    /// The SPH-rigid coupling currently installed, if any.
    pub fn sph_rigid_coupling(&self) -> Option<&SphRigidCouplingParams> {
        self.sph_rigid_coupling.as_ref()
    }

    /// Collider pairs in which a sensor overlapped another collider during the
    /// last [`step`](Self::step).
    ///
    /// Sensors take part in collision detection but never in the contact
    /// solve: a sensor does not push bodies and bodies do not push it. Each
    /// entry is `(collider, collider)` in collider-set order; at least one of
    /// the two is a sensor. Pairs between bodies that can neither move (for
    /// example a static sensor and static geometry) are not reported. The list
    /// is rebuilt on every step.
    pub fn last_sensor_overlaps(&self) -> &[(ColliderHandle, ColliderHandle)] {
        &self.sensor_overlaps
    }

    /// The accumulated contact impulses that warm-start the next step.
    pub fn warm_start_cache(&self) -> &WarmStartCache {
        &self.contact_cache
    }

    /// Take a checkpoint of the pipeline's mutable state.
    ///
    /// The caller is responsible for saving the corresponding [`RigidBodySet`]
    /// and [`ColliderSet`] state if full scene restoration is required.
    pub fn checkpoint(&self) -> PipelineCheckpoint {
        PipelineCheckpoint {
            time: self.time,
            contact_cache: self.contact_cache.clone(),
        }
    }

    /// Restore the pipeline to a previously captured [`PipelineCheckpoint`].
    ///
    /// Statistics accumulated after the checkpoint point are not rolled back.
    pub fn restore_checkpoint(&mut self, ckpt: PipelineCheckpoint) {
        self.time = ckpt.time;
        self.contact_cache = ckpt.contact_cache;
    }

    /// Reset the accumulated [`PipelineStats`] to zero.
    pub fn reset_stats(&mut self) {
        self.stats = PipelineStats::default();
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Record a single stage duration into stats + history.
    fn record_stage(&mut self, stage: &'static str, duration: Duration) {
        // Update rolling timing sums.
        match stage {
            "forces" => self.stats.timings.forces += duration,
            "coupling" => self.stats.timings.coupling += duration,
            "broadphase" => self.stats.timings.broadphase += duration,
            "narrowphase" => self.stats.timings.narrowphase += duration,
            "solver" => self.stats.timings.solver += duration,
            "integration" => self.stats.timings.integration += duration,
            "sleep_check" => self.stats.timings.sleep_check += duration,
            _ => {}
        }

        if self.profiling_enabled && self.profile_history_limit > 0 {
            let rec = StageProfile {
                step: self.stats.step_count,
                sim_time: self.time,
                duration,
                stage,
            };
            while self.profile_history.len() >= self.profile_history_limit {
                self.profile_history.pop_front();
            }
            self.profile_history.push_back(rec);
        }
    }

    /// Invoke the stage callback if one is registered.
    fn fire_callback(&self, stage: &str, duration: Duration) {
        if let Some(cb) = &self.stage_callback {
            cb(stage, duration);
        }
    }

    /// Record a stage duration and invoke the stage callback.
    fn record_and_fire(&mut self, stage: &'static str, duration: Duration) {
        self.record_stage(stage, duration);
        self.fire_callback(stage, duration);
    }

    // -----------------------------------------------------------------------
    // Public step
    // -----------------------------------------------------------------------

    /// Run one full simulation step.
    ///
    /// Pipeline:
    /// 1. Apply gravity + external forces -> velocity
    /// 2. SPH-rigid buoyancy, if a coupling region is installed
    /// 3. Broadphase with rotation-aware AABBs and narrowphase contacts
    /// 4. Sequential-impulse contact solve (the full contact set, also in
    ///    steps that contain CCD events)
    /// 5. Integrate positions; with CCD enabled, the times of impact of fast
    ///    pairs are resolved first, each with its own sub-step (see the module
    ///    documentation)
    /// 6. Sleep check
    pub fn step(&mut self, dt: Real, bodies: &mut RigidBodySet, colliders: &ColliderSet) {
        let gravity = self.config.gravity;
        self.sensor_overlaps.clear();
        self.stats.last_ccd_events = 0;
        let entries = collect_entries(colliders);

        // Stage 1: forces
        let t0 = Instant::now();
        for (_, body) in bodies.iter_mut() {
            body.integrate_forces(dt, &gravity);
        }
        self.record_and_fire("forces", t0.elapsed());

        // Stage 2: SPH-rigid coupling
        if let Some(params) = self.sph_rigid_coupling.clone() {
            let t_coupling = Instant::now();
            self.apply_sph_rigid_coupling(&params, dt, bodies, &entries);
            self.record_and_fire("coupling", t_coupling.elapsed());
        }

        // Stages 3-5: broadphase, narrowphase and the contact solve at the
        // positions of the start of the step. This is the same for every step,
        // with or without CCD.
        self.contact_pass(dt, bodies, &entries);

        // Stage 6: integration. With CCD the fast pairs are resolved by
        // time-of-impact sub-steps first; `interval` is what is left of the
        // step afterwards.
        let interval = if self.config.ccd_enabled {
            self.resolve_ccd_events(dt, bodies, &entries)
        } else {
            dt
        };
        self.integrate(dt, interval, bodies);

        // Stage: sleep check
        let t_sleep = Instant::now();
        for (_, body) in bodies.iter_mut() {
            body.check_sleep(
                dt,
                self.config.linear_sleep_threshold,
                self.config.angular_sleep_threshold,
                self.config.time_before_sleep,
            );
        }
        self.record_and_fire("sleep_check", t_sleep.elapsed());

        self.time += dt;
        self.stats.step_count += 1;
    }

    /// Contact pass: broadphase, narrowphase and sequential-impulse solve at
    /// the current positions.
    ///
    /// `step_dt` is the length of the whole step (it scales the warm-start
    /// impulses). Positions are left alone except for the split-impulse
    /// separation of overlapping bodies; the solved velocities are integrated
    /// by [`integrate`](Self::integrate).
    fn contact_pass(
        &mut self,
        step_dt: Real,
        bodies: &mut RigidBodySet,
        entries: &[ColliderEntry<'_>],
    ) {
        // Broadphase + narrowphase
        let mut detection = self.detect_collisions(entries, bodies);
        self.stats.last_broadphase_pairs = detection.broadphase_pairs;
        self.stats.last_contact_manifolds = detection.contacts.len();
        self.record_and_fire("broadphase", detection.broadphase_time);
        self.record_and_fire("narrowphase", detection.narrowphase_time);
        self.sensor_overlaps = std::mem::take(&mut detection.sensor_overlaps);

        // Solver
        let t_solver = Instant::now();
        let pairs: Vec<ContactPair<'_>> = detection
            .contacts
            .iter()
            .filter_map(|record| {
                let entry_a = entries.get(record.a)?;
                let entry_b = entries.get(record.b)?;
                Some(ContactPair {
                    handle_a: entry_a.handle,
                    handle_b: entry_b.handle,
                    body_a: entry_a.collider.body_handle,
                    body_b: entry_b.collider.body_handle,
                    friction: (entry_a.collider.friction + entry_b.collider.friction) * 0.5,
                    restitution: (entry_a.collider.restitution + entry_b.collider.restitution)
                        * 0.5,
                    manifold: &record.manifold,
                })
            })
            .collect();
        solve_contacts(
            &mut self.contact_cache,
            &pairs,
            bodies,
            SolverSettings {
                iterations: self.config.solver_iterations,
                step_dt,
            },
        );
        self.record_and_fire("solver", t_solver.elapsed());
    }

    /// Integrate the velocities of all bodies over `interval`.
    ///
    /// `interval` is the part of the step of length `step_dt` that is still to
    /// be integrated: the whole step without CCD, the remainder after the CCD
    /// sub-steps otherwise.
    ///
    /// A body that follows a scripted kinematic target stood still while the
    /// step was split into CCD sub-steps and covers its whole way to the target
    /// now, so the velocity it ends up with is the displacement over the whole
    /// step. Integrating it over the remainder only would inflate that velocity
    /// by `step_dt / interval`.
    fn integrate(&mut self, step_dt: Real, interval: Real, bodies: &mut RigidBodySet) {
        if step_dt > 0.0 {
            let t_integration = Instant::now();
            for (_, body) in bodies.iter_mut() {
                let span = if follows_target(body) {
                    step_dt
                } else {
                    interval
                };
                if span > 0.0 {
                    body.integrate_velocity(span);
                }
            }
            self.record_and_fire("integration", t_integration.elapsed());
        }
    }
}

impl Default for PhysicsPipeline {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PipelineBuilder – fluent pipeline factory
// ---------------------------------------------------------------------------

/// Builder that constructs and configures a [`PhysicsPipeline`].
///
/// # Example
/// ```
/// use oxiphysics::core::math::Vec3;
/// use oxiphysics::pipeline::PhysicsPipeline;
///
/// let pipeline = PhysicsPipeline::builder()
///     .gravity(Vec3::new(0.0, -9.81, 0.0))
///     .ccd_enabled(true)
///     .profiling(true)
///     .finish();
/// assert!(pipeline.profiling_enabled);
/// ```
pub struct PipelineBuilder {
    config_builder: PipelineConfigBuilder,
    profiling: bool,
    profile_limit: Option<usize>,
    sph_coupling: Option<SphRigidCouplingParams>,
    max_ccd_substeps: Option<usize>,
}

impl PipelineBuilder {
    /// Start a new builder with defaults.
    pub fn new() -> Self {
        Self {
            config_builder: PipelineConfigBuilder::new(),
            profiling: false,
            profile_limit: None,
            sph_coupling: None,
            max_ccd_substeps: None,
        }
    }

    /// Set gravity.
    pub fn gravity(mut self, g: Vec3) -> Self {
        self.config_builder = self.config_builder.gravity(g);
        self
    }

    /// Enable or disable CCD.
    pub fn ccd_enabled(mut self, enabled: bool) -> Self {
        self.config_builder = self.config_builder.ccd_enabled(enabled);
        self
    }

    /// Set solver iterations.
    pub fn solver_iterations(mut self, n: u32) -> Self {
        self.config_builder = self.config_builder.solver_iterations(n);
        self
    }

    /// Set the maximum number of time-of-impact events a CCD step resolves
    /// (see [`PhysicsPipeline::max_ccd_substeps`]).
    pub fn max_ccd_substeps(mut self, n: usize) -> Self {
        self.max_ccd_substeps = Some(n);
        self
    }

    /// Enable per-stage profiling with an optional ring-buffer size limit.
    pub fn profiling(mut self, enabled: bool) -> Self {
        self.profiling = enabled;
        self
    }

    /// Set the maximum number of profile records to retain.
    pub fn profile_limit(mut self, limit: usize) -> Self {
        self.profile_limit = Some(limit);
        self
    }

    /// Configure SPH-rigid coupling.
    pub fn sph_rigid_coupling(mut self, params: SphRigidCouplingParams) -> Self {
        self.sph_coupling = Some(params);
        self
    }

    /// Consume the builder and produce a [`PhysicsPipeline`].
    pub fn finish(self) -> PhysicsPipeline {
        let config = self.config_builder.build();
        let mut pipeline = PhysicsPipeline::with_config(config);
        if let Some(n) = self.max_ccd_substeps {
            pipeline.max_ccd_substeps = n;
        }
        if self.profiling {
            pipeline.enable_profiling(self.profile_limit);
        }
        if let Some(params) = self.sph_coupling {
            pipeline.set_sph_rigid_coupling(params);
        }
        pipeline
    }
}

impl Default for PipelineBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Diagnostics helpers
// ---------------------------------------------------------------------------

/// Estimate the total kinetic energy of all dynamic bodies in a set.
pub fn total_kinetic_energy(bodies: &RigidBodySet) -> Real {
    let mut ek = 0.0;
    for (_, body) in bodies.iter() {
        if body.is_dynamic() {
            let m = if body.inverse_mass > 1e-15 {
                1.0 / body.inverse_mass
            } else {
                0.0
            };
            ek += 0.5 * m * body.velocity.norm_squared();
        }
    }
    ek
}

/// Count how many dynamic bodies are currently awake (not sleeping).
pub fn count_awake_bodies(bodies: &RigidBodySet) -> usize {
    bodies
        .iter()
        .filter(|(_, b)| b.is_dynamic() && b.state != BodyState::Sleeping)
        .count()
}
