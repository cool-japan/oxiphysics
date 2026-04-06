// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Full physics simulation pipeline.
//!
//! Orchestrates the complete rigid body simulation:
//! 1. Apply gravity / external forces
//! 2. Broadphase collision detection
//! 3. Narrowphase contact generation
//! 4. Solve constraints
//! 5. Integrate velocities + positions
//! 6. Update broadphase
//! 7. Sleep check
//!
//! Extended with:
//! - Pipeline statistics (step count, timing per stage)
//! - Pipeline configuration builder
//! - Multi-physics coupling stubs (rigid + fluid, SPH + rigid)
//! - Pipeline checkpointing and replay
//! - Debug callback hooks
//! - Per-stage performance profiling

use oxiphysics_collision::broadphase::SweepAndPrune;
use oxiphysics_collision::ccd::{ToiResult, time_of_impact};
use oxiphysics_collision::{BroadPhase, CollisionPair, ContactManifold, NarrowPhaseDispatcher};
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_core::{Aabb, PhysicsConfig, Transform};
use oxiphysics_rigid::{BodyState, ColliderSet, RigidBodySet};
use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Key for the contact cache: sorted pair of collider indices.
type ContactKey = (usize, usize);

/// Minimum relative speed (m/step) for CCD to be attempted on a pair.
const CCD_SPEED_THRESHOLD: Real = 0.5;

// ---------------------------------------------------------------------------
// Pipeline statistics
// ---------------------------------------------------------------------------

/// Per-stage timing accumulated over multiple steps.
#[derive(Debug, Clone, Default)]
pub struct StageTimings {
    /// Accumulated wall-clock time spent in the force-integration stage.
    pub forces: Duration,
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
    /// Accumulated per-stage timings (wall-clock, not simulated time).
    pub timings: StageTimings,
}

impl PipelineStats {
    /// Average time (wall-clock) per stage per step, or zero if no steps taken.
    #[allow(dead_code)]
    pub fn avg_stage_time(&self, stage: &str) -> Duration {
        if self.step_count == 0 {
            return Duration::ZERO;
        }
        let n = self.step_count as u32;
        match stage {
            "forces" => self.timings.forces / n,
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
/// ```ignore
/// let config = PipelineConfigBuilder::new()
///     .gravity(Vec3::new(0.0, -9.81, 0.0))
///     .solver_iterations(10)
///     .ccd_enabled(true)
///     .build();
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
    pub contact_cache: HashMap<ContactKey, f64>,
}

// ---------------------------------------------------------------------------
// Multi-physics coupling descriptors
// ---------------------------------------------------------------------------

/// Descriptor for a region where an SPH fluid couples with rigid bodies.
///
/// In a full multi-physics simulation the pipeline can apply buoyancy forces
/// to rigid bodies based on the local SPH fluid density and vice-versa push
/// SPH particles away from rigid body surfaces.  This struct holds the
/// parameters that govern that exchange.
#[derive(Debug, Clone)]
pub struct SphRigidCouplingParams {
    /// Average fluid density in the coupling region \[kg/m³\].
    pub fluid_density: Real,
    /// Coupling strength coefficient (0 = no coupling, 1 = full).
    pub coupling_alpha: Real,
    /// Half-extents of the axis-aligned coupling region.
    pub region_half_extents: Vec3,
    /// Centre of the coupling region in world space.
    pub region_centre: Vec3,
}

impl Default for SphRigidCouplingParams {
    fn default() -> Self {
        Self {
            fluid_density: 1000.0,
            coupling_alpha: 1.0,
            region_half_extents: Vec3::new(5.0, 5.0, 5.0),
            region_centre: Vec3::zeros(),
        }
    }
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
    /// Broadphase algorithm.
    broadphase: SweepAndPrune,
    /// Warmstarting cache: accumulated normal impulse per contact pair from the previous frame.
    #[allow(dead_code)]
    contact_cache: HashMap<ContactKey, f64>,
    /// Rolling statistics collected during simulation.
    pub stats: PipelineStats,
    /// Optional debug callback invoked after each stage.
    #[allow(dead_code)]
    stage_callback: Option<StageCallback>,
    /// Per-stage profile history (bounded ring-buffer of the last N records).
    #[allow(dead_code)]
    profile_history: Vec<StageProfile>,
    /// Maximum number of stage profile records to retain.
    #[allow(dead_code)]
    profile_history_limit: usize,
    /// Whether per-stage profiling is enabled.
    pub profiling_enabled: bool,
    /// Active SPH-rigid coupling parameters, if any.
    #[allow(dead_code)]
    sph_rigid_coupling: Option<SphRigidCouplingParams>,
}

impl PhysicsPipeline {
    /// Create a new pipeline with default config.
    pub fn new() -> Self {
        Self {
            config: PhysicsConfig::default(),
            time: 0.0,
            broadphase: SweepAndPrune::default(),
            contact_cache: HashMap::new(),
            stats: PipelineStats::default(),
            stage_callback: None,
            profile_history: Vec::new(),
            profile_history_limit: 1024,
            profiling_enabled: false,
            sph_rigid_coupling: None,
        }
    }

    /// Create with custom config.
    pub fn with_config(config: PhysicsConfig) -> Self {
        Self {
            config,
            time: 0.0,
            broadphase: SweepAndPrune::default(),
            contact_cache: HashMap::new(),
            stats: PipelineStats::default(),
            stage_callback: None,
            profile_history: Vec::new(),
            profile_history_limit: 1024,
            profiling_enabled: false,
            sph_rigid_coupling: None,
        }
    }

    /// Use a [`PipelineConfigBuilder`] to construct and configure a pipeline.
    ///
    /// # Example
    /// ```ignore
    /// let pipeline = PhysicsPipeline::builder()
    ///     .gravity(Vec3::new(0.0, -9.81, 0.0))
    ///     .ccd_enabled(true)
    ///     .finish();
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
    /// When `limit` is `None` the default of 1024 records is used.
    pub fn enable_profiling(&mut self, limit: Option<usize>) {
        self.profiling_enabled = true;
        if let Some(l) = limit {
            self.profile_history_limit = l;
        }
    }

    /// Disable per-stage profiling.
    pub fn disable_profiling(&mut self) {
        self.profiling_enabled = false;
    }

    /// Drain and return all accumulated profile records, leaving an empty history.
    pub fn drain_profiles(&mut self) -> Vec<StageProfile> {
        std::mem::take(&mut self.profile_history)
    }

    /// Configure SPH-rigid coupling.
    pub fn set_sph_rigid_coupling(&mut self, params: SphRigidCouplingParams) {
        self.sph_rigid_coupling = Some(params);
    }

    /// Disable SPH-rigid coupling.
    pub fn clear_sph_rigid_coupling(&mut self) {
        self.sph_rigid_coupling = None;
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
            "broadphase" => self.stats.timings.broadphase += duration,
            "narrowphase" => self.stats.timings.narrowphase += duration,
            "solver" => self.stats.timings.solver += duration,
            "integration" => self.stats.timings.integration += duration,
            "sleep_check" => self.stats.timings.sleep_check += duration,
            _ => {}
        }

        // Fire the debug callback (borrow checked workaround: re-check below).
        if self.profiling_enabled {
            let rec = StageProfile {
                step: self.stats.step_count,
                sim_time: self.time,
                duration,
                stage,
            };
            if self.profile_history.len() >= self.profile_history_limit {
                self.profile_history.remove(0);
            }
            self.profile_history.push(rec);
        }
    }

    /// Invoke the stage callback if one is registered.
    fn fire_callback(&self, stage: &str, duration: Duration) {
        if let Some(cb) = &self.stage_callback {
            cb(stage, duration);
        }
    }

    // -----------------------------------------------------------------------
    // Public step
    // -----------------------------------------------------------------------

    /// Run one full simulation step.
    ///
    /// Pipeline (CCD enabled):
    /// 1. Apply gravity + external forces -> velocity
    /// 2. Broadphase with swept AABBs
    /// 3. CCD: find earliest TOI among candidate pairs
    /// 4. If TOI found at t < 1.0: sub-step to t*dt, apply impulse, continue (1-t)*dt
    /// 5. Else: normal narrowphase + impulse solve
    /// 6. Integrate positions
    /// 7. Sleep check
    pub fn step(&mut self, dt: Real, bodies: &mut RigidBodySet, colliders: &ColliderSet) {
        let gravity = self.config.gravity;

        // Stage 1: forces
        let t0 = Instant::now();
        for (_, body) in bodies.iter_mut() {
            body.integrate_forces(dt, &gravity);
        }
        let d_forces = t0.elapsed();
        self.record_stage("forces", d_forces);
        self.fire_callback("forces", d_forces);

        if self.config.ccd_enabled {
            self.step_with_ccd(dt, bodies, colliders);
        } else {
            self.step_normal(dt, bodies, colliders);
        }

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
        let d_sleep = t_sleep.elapsed();
        self.record_stage("sleep_check", d_sleep);
        self.fire_callback("sleep_check", d_sleep);

        self.time += dt;
        self.stats.step_count += 1;
    }

    /// Normal (non-CCD) step: narrowphase + solve + integrate.
    fn step_normal(&mut self, dt: Real, bodies: &mut RigidBodySet, colliders: &ColliderSet) {
        // Broadphase + narrowphase
        let t_np = Instant::now();
        let manifolds = self.detect_collisions(bodies, colliders);
        let d_np = t_np.elapsed();
        self.stats.last_contact_manifolds = manifolds.len();
        self.record_stage("narrowphase", d_np);
        self.fire_callback("narrowphase", d_np);

        // Solver
        let t_sol = Instant::now();
        self.solve_contacts(&manifolds, bodies, colliders, dt);
        let d_sol = t_sol.elapsed();
        self.record_stage("solver", d_sol);
        self.fire_callback("solver", d_sol);

        // Integration
        let t_int = Instant::now();
        for (_, body) in bodies.iter_mut() {
            body.integrate_velocity(dt);
        }
        let d_int = t_int.elapsed();
        self.record_stage("integration", d_int);
        self.fire_callback("integration", d_int);
    }

    /// CCD-aware step.
    fn step_with_ccd(&mut self, dt: Real, bodies: &mut RigidBodySet, colliders: &ColliderSet) {
        // Broadphase + CCD
        let t_bp = Instant::now();
        let ccd_hits = self.detect_collisions_ccd(dt, bodies, colliders);
        let d_bp = t_bp.elapsed();
        self.stats.last_broadphase_pairs = ccd_hits.len();
        self.record_stage("broadphase", d_bp);
        self.fire_callback("broadphase", d_bp);

        // Narrowphase: pick earliest TOI
        let t_np = Instant::now();
        let earliest = ccd_hits
            .iter()
            .filter(|(toi_result, _, _)| toi_result.toi < 1.0)
            .min_by(|(a, _, _), (b, _, _)| {
                a.toi
                    .partial_cmp(&b.toi)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        let d_np = t_np.elapsed();
        self.record_stage("narrowphase", d_np);
        self.fire_callback("narrowphase", d_np);

        if let Some((toi_result, col_a, col_b)) = earliest {
            let t = toi_result.toi.clamp(0.0, 1.0);
            let dt_first = t * dt;
            let dt_rest = (1.0 - t) * dt;

            // Sub-step 1: integrate positions to TOI
            let t_int = Instant::now();
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt_first);
            }
            let d_int = t_int.elapsed();
            self.record_stage("integration", d_int);
            self.fire_callback("integration", d_int);

            // Apply collision impulse at TOI
            let t_sol = Instant::now();
            self.apply_ccd_impulse(toi_result, *col_a, *col_b, bodies, colliders);
            let d_sol = t_sol.elapsed();
            self.record_stage("solver", d_sol);
            self.fire_callback("solver", d_sol);

            // Sub-step 2: integrate remaining time
            if dt_rest > 1e-10 {
                let t_int2 = Instant::now();
                for (_, body) in bodies.iter_mut() {
                    body.integrate_velocity(dt_rest);
                }
                let d_int2 = t_int2.elapsed();
                self.record_stage("integration", d_int2);
                self.fire_callback("integration", d_int2);
            }
        } else {
            // No CCD event — fall back to normal narrowphase.
            let t_np2 = Instant::now();
            let manifolds = self.detect_collisions(bodies, colliders);
            let d_np2 = t_np2.elapsed();
            self.stats.last_contact_manifolds = manifolds.len();
            self.record_stage("narrowphase", d_np2);
            self.fire_callback("narrowphase", d_np2);

            let t_sol = Instant::now();
            self.solve_contacts(&manifolds, bodies, colliders, dt);
            let d_sol = t_sol.elapsed();
            self.record_stage("solver", d_sol);
            self.fire_callback("solver", d_sol);

            let t_int = Instant::now();
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt);
            }
            let d_int = t_int.elapsed();
            self.record_stage("integration", d_int);
            self.fire_callback("integration", d_int);
        }
    }

    /// For each broadphase pair, compute swept AABBs and run CCD.
    /// Returns `(ToiResult, collider_index_a, collider_index_b)` for hits.
    #[allow(dead_code)]
    fn detect_collisions_ccd(
        &self,
        dt: Real,
        bodies: &RigidBodySet,
        colliders: &ColliderSet,
    ) -> Vec<(ToiResult, usize, usize)> {
        let collider_data: Vec<_> = colliders.iter().collect();
        if collider_data.is_empty() {
            return Vec::new();
        }

        // Build swept (start + end) transforms and swept AABBs for broadphase.
        let swept_aabbs: Vec<Aabb> = collider_data
            .iter()
            .map(|(_, col)| {
                let body = col.body_handle.and_then(|h| bodies.get(h));

                let body_t_start = body.map(|b| b.transform.clone()).unwrap_or_default();

                // Predict end position based on current velocity * dt.
                let body_t_end = body
                    .map(|b| {
                        let end_pos = b.transform.position + b.velocity * dt;
                        Transform::from_position(end_pos)
                    })
                    .unwrap_or_default();

                let local_aabb = col.shape.bounding_box();

                // AABB at start position.
                let world_t_start = col.world_transform(&body_t_start);
                let aabb_start = Aabb::new(
                    local_aabb.min + world_t_start.position,
                    local_aabb.max + world_t_start.position,
                );

                // AABB at end position.
                let world_t_end = col.world_transform(&body_t_end);
                let aabb_end = Aabb::new(
                    local_aabb.min + world_t_end.position,
                    local_aabb.max + world_t_end.position,
                );

                // Swept AABB = union of start and end.
                aabb_start.merge(&aabb_end)
            })
            .collect();

        let pairs = self.broadphase.find_pairs(&swept_aabbs);

        let mut results = Vec::new();

        for pair in &pairs {
            let (_, col_a) = &collider_data[pair.a];
            let (_, col_b) = &collider_data[pair.b];

            if col_a.is_sensor && col_b.is_sensor {
                continue;
            }

            let body_a = col_a.body_handle.and_then(|h| bodies.get(h));
            let body_b = col_b.body_handle.and_then(|h| bodies.get(h));

            // At least one must be dynamic and moving.
            let vel_a = body_a.map(|b| b.velocity).unwrap_or(Vec3::zeros());
            let vel_b = body_b.map(|b| b.velocity).unwrap_or(Vec3::zeros());
            let rel_speed = (vel_a - vel_b).norm() * dt;
            if rel_speed < CCD_SPEED_THRESHOLD {
                continue;
            }

            // Build start/end transforms for each collider.
            let body_t_a_start = body_a.map(|b| b.transform.clone()).unwrap_or_default();
            let body_t_a_end = body_a
                .map(|b| {
                    let end_pos = b.transform.position + b.velocity * dt;
                    Transform::from_position(end_pos)
                })
                .unwrap_or_default();

            let body_t_b_start = body_b.map(|b| b.transform.clone()).unwrap_or_default();
            let body_t_b_end = body_b
                .map(|b| {
                    let end_pos = b.transform.position + b.velocity * dt;
                    Transform::from_position(end_pos)
                })
                .unwrap_or_default();

            let t_a_start = col_a.world_transform(&body_t_a_start);
            let t_a_end = col_a.world_transform(&body_t_a_end);
            let t_b_start = col_b.world_transform(&body_t_b_start);
            let t_b_end = col_b.world_transform(&body_t_b_end);

            if let Some(toi) = time_of_impact(
                col_a.shape.as_ref(),
                &t_a_start,
                &t_a_end,
                col_b.shape.as_ref(),
                &t_b_start,
                &t_b_end,
            ) {
                results.push((toi, pair.a, pair.b));
            }
        }

        results
    }

    /// Apply a collision impulse derived from a CCD TOI result.
    fn apply_ccd_impulse(
        &self,
        toi_result: &ToiResult,
        col_idx_a: usize,
        col_idx_b: usize,
        bodies: &mut RigidBodySet,
        colliders: &ColliderSet,
    ) {
        let collider_data: Vec<_> = colliders.iter().collect();
        let (_, col_a) = &collider_data[col_idx_a];
        let (_, col_b) = &collider_data[col_idx_b];

        let ha = col_a.body_handle;
        let hb = col_b.body_handle;

        let restitution = (col_a.restitution + col_b.restitution) * 0.5;

        // Use the CCD normal and the midpoint of witness points as contact point.
        let normal = -toi_result.normal;
        let cp = (toi_result.witness_a + toi_result.witness_b) * 0.5;

        let (vel_a, omega_a, inv_mass_a, inv_inertia_a, pos_a, dyn_a) = ha
            .and_then(|h| bodies.get(h))
            .map(|b| {
                (
                    b.velocity,
                    b.angular_velocity,
                    b.inverse_mass,
                    b.world_inverse_inertia,
                    b.transform.position,
                    b.is_dynamic(),
                )
            })
            .unwrap_or((
                Vec3::zeros(),
                Vec3::zeros(),
                0.0,
                oxiphysics_core::math::Mat3::zeros(),
                Vec3::zeros(),
                false,
            ));

        let (vel_b, omega_b, inv_mass_b, inv_inertia_b, pos_b, dyn_b) = hb
            .and_then(|h| bodies.get(h))
            .map(|b| {
                (
                    b.velocity,
                    b.angular_velocity,
                    b.inverse_mass,
                    b.world_inverse_inertia,
                    b.transform.position,
                    b.is_dynamic(),
                )
            })
            .unwrap_or((
                Vec3::zeros(),
                Vec3::zeros(),
                0.0,
                oxiphysics_core::math::Mat3::zeros(),
                Vec3::zeros(),
                false,
            ));

        let r_a = cp - pos_a;
        let r_b = cp - pos_b;

        let vcp_a = vel_a + omega_a.cross(&r_a);
        let vcp_b = vel_b + omega_b.cross(&r_b);
        let rel_vel = vcp_a - vcp_b;
        let normal_vel = rel_vel.dot(&normal);

        // Only apply impulse if bodies are approaching.
        if normal_vel >= 0.0 {
            return;
        }

        let rn_a = r_a.cross(&normal);
        let rn_b = r_b.cross(&normal);
        let ang_a = rn_a.dot(&(inv_inertia_a * rn_a));
        let ang_b = rn_b.dot(&(inv_inertia_b * rn_b));
        let eff_mass_denom = inv_mass_a + inv_mass_b + ang_a + ang_b;

        if eff_mass_denom < 1e-10 {
            return;
        }

        let j_normal = -(1.0 + restitution) * normal_vel / eff_mass_denom;
        let normal_impulse = normal * j_normal;

        if let Some(h) = ha
            && let Some(body) = bodies.get_mut(h)
            && dyn_a
        {
            body.velocity += normal_impulse * body.inverse_mass;
            body.angular_velocity += body.world_inverse_inertia * r_a.cross(&normal_impulse);
            body.wake_up();
        }
        if let Some(h) = hb
            && let Some(body) = bodies.get_mut(h)
            && dyn_b
        {
            body.velocity -= normal_impulse * body.inverse_mass;
            body.angular_velocity -= body.world_inverse_inertia * r_b.cross(&normal_impulse);
            body.wake_up();
        }
    }

    /// Detect collisions and return contact manifolds.
    fn detect_collisions(
        &self,
        bodies: &RigidBodySet,
        colliders: &ColliderSet,
    ) -> Vec<(ContactManifold, usize, usize)> {
        // Collect collider data with world transforms and AABBs
        let collider_data: Vec<_> = colliders.iter().collect();
        if collider_data.is_empty() {
            return Vec::new();
        }

        let aabbs: Vec<Aabb> = collider_data
            .iter()
            .map(|(_, col)| {
                let body_t = col
                    .body_handle
                    .and_then(|h| bodies.get(h))
                    .map(|b| b.transform.clone())
                    .unwrap_or_default();
                let world_t = col.world_transform(&body_t);
                let local_aabb = col.shape.bounding_box();
                Aabb::new(
                    local_aabb.min + world_t.position,
                    local_aabb.max + world_t.position,
                )
            })
            .collect();

        // Broadphase
        let pairs = self.broadphase.find_pairs(&aabbs);

        // Narrowphase
        let mut manifolds = Vec::new();
        for pair in &pairs {
            let (_, col_a) = &collider_data[pair.a];
            let (_, col_b) = &collider_data[pair.b];

            // Skip sensor-sensor pairs
            if col_a.is_sensor && col_b.is_sensor {
                continue;
            }

            let body_t_a = col_a
                .body_handle
                .and_then(|h| bodies.get(h))
                .map(|b| b.transform.clone())
                .unwrap_or_default();
            let body_t_b = col_b
                .body_handle
                .and_then(|h| bodies.get(h))
                .map(|b| b.transform.clone())
                .unwrap_or_default();

            let t_a = col_a.world_transform(&body_t_a);
            let t_b = col_b.world_transform(&body_t_b);

            if let Some(manifold) = NarrowPhaseDispatcher::generate_contacts(
                col_a.shape.as_ref(),
                &t_a,
                col_b.shape.as_ref(),
                &t_b,
                CollisionPair::new(pair.a, pair.b),
            ) {
                manifolds.push((manifold, pair.a, pair.b));
            }
        }
        manifolds
    }

    /// Sequential impulse solver for contacts with full angular response and warmstarting.
    fn solve_contacts(
        &mut self,
        manifolds: &[(ContactManifold, usize, usize)],
        bodies: &mut RigidBodySet,
        colliders: &ColliderSet,
        _dt: Real,
    ) {
        let collider_data: Vec<_> = colliders.iter().collect();

        // --- Warmstarting ---
        let mut active_keys: HashMap<ContactKey, f64> = HashMap::new();

        for (manifold, col_idx_a, col_idx_b) in manifolds {
            let key: ContactKey = if col_idx_a <= col_idx_b {
                (*col_idx_a, *col_idx_b)
            } else {
                (*col_idx_b, *col_idx_a)
            };

            let (_, col_a) = &collider_data[*col_idx_a];
            let (_, col_b) = &collider_data[*col_idx_b];
            let ha = col_a.body_handle;
            let hb = col_b.body_handle;

            if let Some(&cached_impulse) = self.contact_cache.get(&key) {
                let n_contacts = manifold.contacts.len();
                if n_contacts == 0 {
                    continue;
                }
                let impulse_per_contact = cached_impulse / n_contacts as f64;

                for contact in &manifold.contacts {
                    let cp = contact.point();
                    let normal = contact.normal;

                    let (inv_mass_a, _inv_inertia_a, pos_a, dyn_a) = ha
                        .and_then(|h| bodies.get(h))
                        .map(|b| {
                            (
                                b.inverse_mass,
                                b.world_inverse_inertia,
                                b.transform.position,
                                b.is_dynamic(),
                            )
                        })
                        .unwrap_or((
                            0.0,
                            oxiphysics_core::math::Mat3::zeros(),
                            Vec3::zeros(),
                            false,
                        ));

                    let (inv_mass_b, _inv_inertia_b, pos_b, dyn_b) = hb
                        .and_then(|h| bodies.get(h))
                        .map(|b| {
                            (
                                b.inverse_mass,
                                b.world_inverse_inertia,
                                b.transform.position,
                                b.is_dynamic(),
                            )
                        })
                        .unwrap_or((
                            0.0,
                            oxiphysics_core::math::Mat3::zeros(),
                            Vec3::zeros(),
                            false,
                        ));

                    let r_a = cp - pos_a;
                    let r_b = cp - pos_b;
                    let warm_impulse = normal * impulse_per_contact;

                    if let Some(h) = ha
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_a
                    {
                        body.velocity += warm_impulse * inv_mass_a;
                        body.angular_velocity +=
                            body.world_inverse_inertia * r_a.cross(&warm_impulse);
                        body.wake_up();
                    }
                    if let Some(h) = hb
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_b
                    {
                        body.velocity -= warm_impulse * inv_mass_b;
                        body.angular_velocity -=
                            body.world_inverse_inertia * r_b.cross(&warm_impulse);
                        body.wake_up();
                    }
                }
            }

            active_keys.insert(key, 0.0);
        }

        let solver_iterations = self.config.solver_iterations;

        // --- PGS iterations ---
        for iteration in 0..solver_iterations {
            let is_last = iteration + 1 == solver_iterations;

            for (manifold, col_idx_a, col_idx_b) in manifolds {
                let key: ContactKey = if col_idx_a <= col_idx_b {
                    (*col_idx_a, *col_idx_b)
                } else {
                    (*col_idx_b, *col_idx_a)
                };

                let (_, col_a) = &collider_data[*col_idx_a];
                let (_, col_b) = &collider_data[*col_idx_b];

                let ha = col_a.body_handle;
                let hb = col_b.body_handle;

                let restitution = (col_a.restitution + col_b.restitution) * 0.5;
                let friction = (col_a.friction + col_b.friction) * 0.5;

                for contact in &manifold.contacts {
                    let cp = contact.point();

                    let (vel_a, omega_a, inv_mass_a, inv_inertia_a, pos_a, dyn_a) = ha
                        .and_then(|h| bodies.get(h))
                        .map(|b| {
                            (
                                b.velocity,
                                b.angular_velocity,
                                b.inverse_mass,
                                b.world_inverse_inertia,
                                b.transform.position,
                                b.is_dynamic(),
                            )
                        })
                        .unwrap_or((
                            Vec3::zeros(),
                            Vec3::zeros(),
                            0.0,
                            oxiphysics_core::math::Mat3::zeros(),
                            Vec3::zeros(),
                            false,
                        ));
                    let (vel_b, omega_b, inv_mass_b, inv_inertia_b, pos_b, dyn_b) = hb
                        .and_then(|h| bodies.get(h))
                        .map(|b| {
                            (
                                b.velocity,
                                b.angular_velocity,
                                b.inverse_mass,
                                b.world_inverse_inertia,
                                b.transform.position,
                                b.is_dynamic(),
                            )
                        })
                        .unwrap_or((
                            Vec3::zeros(),
                            Vec3::zeros(),
                            0.0,
                            oxiphysics_core::math::Mat3::zeros(),
                            Vec3::zeros(),
                            false,
                        ));

                    let r_a = cp - pos_a;
                    let r_b = cp - pos_b;

                    let vcp_a = vel_a + omega_a.cross(&r_a);
                    let vcp_b = vel_b + omega_b.cross(&r_b);
                    let rel_vel = vcp_a - vcp_b;
                    let normal_vel = rel_vel.dot(&contact.normal);

                    if normal_vel > 0.0 {
                        continue;
                    }

                    let rn_a = r_a.cross(&contact.normal);
                    let rn_b = r_b.cross(&contact.normal);
                    let ang_a = rn_a.dot(&(inv_inertia_a * rn_a));
                    let ang_b = rn_b.dot(&(inv_inertia_b * rn_b));
                    let eff_mass_denom = inv_mass_a + inv_mass_b + ang_a + ang_b;

                    if eff_mass_denom < 1e-10 {
                        continue;
                    }

                    let j_normal = -(1.0 + restitution) * normal_vel / eff_mass_denom;
                    let normal_impulse = contact.normal * j_normal;

                    if is_last {
                        *active_keys.entry(key).or_insert(0.0) += j_normal;
                    }

                    if let Some(h) = ha
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_a
                    {
                        body.velocity += normal_impulse * body.inverse_mass;
                        body.angular_velocity +=
                            body.world_inverse_inertia * r_a.cross(&normal_impulse);
                        body.wake_up();
                    }
                    if let Some(h) = hb
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_b
                    {
                        body.velocity -= normal_impulse * body.inverse_mass;
                        body.angular_velocity -=
                            body.world_inverse_inertia * r_b.cross(&normal_impulse);
                        body.wake_up();
                    }

                    // Tangent friction impulse
                    let tangent_vel = rel_vel - contact.normal * normal_vel;
                    let tangent_speed = tangent_vel.norm();
                    if tangent_speed > 1e-10 {
                        let tangent = tangent_vel / tangent_speed;
                        let rt_a = r_a.cross(&tangent);
                        let rt_b = r_b.cross(&tangent);
                        let ang_t_a = rt_a.dot(&(inv_inertia_a * rt_a));
                        let ang_t_b = rt_b.dot(&(inv_inertia_b * rt_b));
                        let eff_mass_t = inv_mass_a + inv_mass_b + ang_t_a + ang_t_b;
                        let j_tangent = if eff_mass_t > 1e-10 {
                            (-tangent_speed / eff_mass_t).max(-friction * j_normal)
                        } else {
                            0.0
                        };
                        let friction_impulse = tangent * j_tangent;

                        if let Some(h) = ha
                            && let Some(body) = bodies.get_mut(h)
                            && dyn_a
                        {
                            body.velocity += friction_impulse * body.inverse_mass;
                            body.angular_velocity +=
                                body.world_inverse_inertia * r_a.cross(&friction_impulse);
                        }
                        if let Some(h) = hb
                            && let Some(body) = bodies.get_mut(h)
                            && dyn_b
                        {
                            body.velocity -= friction_impulse * body.inverse_mass;
                            body.angular_velocity -=
                                body.world_inverse_inertia * r_b.cross(&friction_impulse);
                        }
                    }

                    // Baumgarte position correction
                    let slop = 0.005;
                    let baumgarte = 0.2;
                    let corr_mag = baumgarte * (contact.depth - slop).max(0.0) / eff_mass_denom;
                    let correction = contact.normal * corr_mag;

                    if let Some(h) = ha
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_a
                    {
                        body.transform.position += correction * body.inverse_mass;
                    }
                    if let Some(h) = hb
                        && let Some(body) = bodies.get_mut(h)
                        && dyn_b
                    {
                        body.transform.position -= correction * body.inverse_mass;
                    }
                }
            }
        }

        self.contact_cache = active_keys;
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
/// ```ignore
/// let pipeline = PhysicsPipeline::builder()
///     .gravity(Vec3::new(0.0, -9.81, 0.0))
///     .ccd_enabled(true)
///     .profiling(true)
///     .finish();
/// ```
pub struct PipelineBuilder {
    config_builder: PipelineConfigBuilder,
    profiling: bool,
    profile_limit: Option<usize>,
    sph_coupling: Option<SphRigidCouplingParams>,
}

impl PipelineBuilder {
    /// Start a new builder with defaults.
    pub fn new() -> Self {
        Self {
            config_builder: PipelineConfigBuilder::new(),
            profiling: false,
            profile_limit: None,
            sph_coupling: None,
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
// Multi-physics coupling utilities
// ---------------------------------------------------------------------------

/// Apply buoyancy from a uniform fluid region to all dynamic bodies.
///
/// This is a simplified coupling: each dynamic body whose centroid falls
/// inside the coupling region receives an upward buoyancy force proportional
/// to its bounding volume (approximated as a sphere) and the fluid density.
///
/// For production use, couple an actual [`oxiphysics_sph`] simulation by
/// querying SPH densities at each body's centroid.
///
/// `g_magnitude` – magnitude of gravitational acceleration \[m/s²\].
#[allow(dead_code)]
pub fn apply_buoyancy_impulse(
    bodies: &mut RigidBodySet,
    coupling: &SphRigidCouplingParams,
    effective_radius: Real,
    g_magnitude: Real,
    dt: Real,
) {
    let half = coupling.region_half_extents;
    let centre = coupling.region_centre;
    let rho_f = coupling.fluid_density;
    let alpha = coupling.coupling_alpha;

    // Volume of a sphere: (4/3) π r³
    let submerged_vol = (4.0 / 3.0) * std::f64::consts::PI * effective_radius.powi(3);
    let buoyancy_force_mag = rho_f * g_magnitude * submerged_vol * alpha;

    for (_, body) in bodies.iter_mut() {
        if !body.is_dynamic() {
            continue;
        }
        let p = body.transform.position;
        // Check if inside coupling region.
        if (p.x - centre.x).abs() <= half.x
            && (p.y - centre.y).abs() <= half.y
            && (p.z - centre.z).abs() <= half.z
        {
            // Buoyancy acts upward (+Y).
            let impulse = Vec3::new(0.0, buoyancy_force_mag * dt, 0.0);
            body.velocity += impulse * body.inverse_mass;
        }
    }
}

/// Estimate the total kinetic energy of all dynamic bodies in a set.
#[allow(dead_code)]
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
#[allow(dead_code)]
pub fn count_awake_bodies(bodies: &RigidBodySet) -> usize {
    bodies
        .iter()
        .filter(|(_, b)| b.is_dynamic() && b.state != BodyState::Sleeping)
        .count()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::Transform;
    use oxiphysics_geometry::{BoxShape, Shape, Sphere};
    use oxiphysics_rigid::{Collider, RigidBody};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

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
        let mut pipeline = PhysicsPipeline::new();
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        let mut sphere_body = RigidBody::new(1.0);
        sphere_body.transform = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        sphere_body.linear_damping = 0.0;
        let sh = bodies.insert(sphere_body);

        let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(
            Collider::new(sphere_shape)
                .with_body(sh)
                .with_restitution(0.5),
        );

        let mut floor = RigidBody::new_static();
        floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let fh = bodies.insert(floor);

        let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
        colliders.insert(Collider::new(floor_shape).with_body(fh));

        for _ in 0..300 {
            pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        }

        let sphere = bodies.get(sh).unwrap();
        assert!(
            sphere.transform.position.y > -1.0 && sphere.transform.position.y < 5.0,
            "Sphere y={} should be near floor",
            sphere.transform.position.y
        );
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
    fn test_pipeline_ccd_no_tunnel() {
        let config = PhysicsConfig {
            ccd_enabled: true,
            gravity: Vec3::zeros(),
            ..PhysicsConfig::default()
        };

        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        let mut body_a = RigidBody::new(1.0);
        body_a.transform = Transform::from_position(Vec3::new(-5.0, 0.0, 0.0));
        body_a.velocity = Vec3::new(100.0, 0.0, 0.0);
        body_a.linear_damping = 0.0;
        let ha = bodies.insert(body_a);

        let shape_a: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(Collider::new(shape_a).with_body(ha).with_restitution(1.0));

        let mut body_b = RigidBody::new(1.0);
        body_b.transform = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        body_b.velocity = Vec3::new(-100.0, 0.0, 0.0);
        body_b.linear_damping = 0.0;
        let hb = bodies.insert(body_b);

        let shape_b: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(Collider::new(shape_b).with_body(hb).with_restitution(1.0));

        let dt = 0.1_f64;
        let x_a_initial = bodies.get(ha).unwrap().transform.position.x;
        let x_b_initial = bodies.get(hb).unwrap().transform.position.x;

        pipeline.step(dt, &mut bodies, &colliders);

        let vel_a_after = bodies.get(ha).unwrap().velocity.x;
        let vel_b_after = bodies.get(hb).unwrap().velocity.x;
        let x_a_after = bodies.get(ha).unwrap().transform.position.x;
        let x_b_after = bodies.get(hb).unwrap().transform.position.x;

        assert!(
            vel_a_after < 0.0,
            "Sphere A velocity should be negative after collision, got {}",
            vel_a_after
        );
        assert!(
            vel_b_after > 0.0,
            "Sphere B velocity should be positive after collision, got {}",
            vel_b_after
        );
        assert!(
            x_a_after <= x_b_after + 1.0 + 1e-3,
            "Sphere A (x={}) tunneled past sphere B (x={})",
            x_a_after,
            x_b_after
        );
        assert!(
            (x_a_after - x_a_initial).abs() > 0.01,
            "Sphere A did not move"
        );
        assert!(
            (x_b_after - x_b_initial).abs() > 0.01,
            "Sphere B did not move"
        );
    }

    #[test]
    fn test_warmstarting_converges_faster() {
        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 1,
            ccd_enabled: false,
            ..PhysicsConfig::default()
        };

        let make_scene = |cfg: PhysicsConfig| {
            let mut bodies = RigidBodySet::new();
            let mut colliders = ColliderSet::new();

            let mut floor = RigidBody::new_static();
            floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
            let fh = bodies.insert(floor);
            let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
            colliders.insert(Collider::new(floor_shape).with_body(fh));

            let mut boxb = RigidBody::new(1.0);
            boxb.transform = Transform::from_position(Vec3::new(0.0, 0.48, 0.0));
            boxb.linear_damping = 0.0;
            let bh = bodies.insert(boxb);
            let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
            colliders.insert(Collider::new(box_shape).with_body(bh));

            (PhysicsPipeline::with_config(cfg), bodies, colliders, bh)
        };

        let (mut pipe_cold, mut bodies_cold, colliders_cold, bh_cold) = make_scene(config.clone());
        pipe_cold.step(1.0 / 60.0, &mut bodies_cold, &colliders_cold);
        pipe_cold.contact_cache.clear();
        let vy_before_cold = bodies_cold.get(bh_cold).unwrap().velocity.y;
        pipe_cold.step(1.0 / 60.0, &mut bodies_cold, &colliders_cold);
        let vy_after_cold = bodies_cold.get(bh_cold).unwrap().velocity.y;
        let residual_cold = vy_after_cold.abs();

        let (mut pipe_warm, mut bodies_warm, colliders_warm, bh_warm) = make_scene(config.clone());
        pipe_warm.step(1.0 / 60.0, &mut bodies_warm, &colliders_warm);
        let _vy_before_warm = bodies_warm.get(bh_warm).unwrap().velocity.y;
        pipe_warm.step(1.0 / 60.0, &mut bodies_warm, &colliders_warm);
        let vy_after_warm = bodies_warm.get(bh_warm).unwrap().velocity.y;
        let residual_warm = vy_after_warm.abs();

        assert!(
            !pipe_warm.contact_cache.is_empty(),
            "contact_cache should be populated after a contact step"
        );
        assert!(
            residual_warm <= residual_cold + 1e-9,
            "warmstarted residual ({}) should be <= cold residual ({}); vy_before_cold={}",
            residual_warm,
            residual_cold,
            vy_before_cold,
        );
    }

    #[test]
    fn test_warmstart_cache_cleared_on_separation() {
        let config = PhysicsConfig {
            gravity: Vec3::zeros(),
            solver_iterations: 4,
            ccd_enabled: false,
            ..PhysicsConfig::default()
        };

        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        let mut floor = RigidBody::new_static();
        floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let fh = bodies.insert(floor);
        let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
        colliders.insert(Collider::new(floor_shape).with_body(fh));

        let mut boxb = RigidBody::new(1.0);
        boxb.transform = Transform::from_position(Vec3::new(0.0, 0.48, 0.0));
        boxb.linear_damping = 0.0;
        let bh = bodies.insert(boxb);
        let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
        colliders.insert(Collider::new(box_shape).with_body(bh));

        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        assert!(
            !pipeline.contact_cache.is_empty(),
            "cache should be non-empty while bodies are in contact"
        );

        bodies.get_mut(bh).unwrap().transform.position.y = 100.0;
        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);

        assert!(
            pipeline.contact_cache.is_empty(),
            "cache should be empty after bodies separate; found {:?}",
            pipeline.contact_cache
        );
    }

    // -----------------------------------------------------------------------
    // New tests: statistics, builder, checkpoint, callbacks
    // -----------------------------------------------------------------------

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
        let mut pipeline = PhysicsPipeline::new();
        let mut bodies = RigidBodySet::new();
        let colliders = ColliderSet::new();
        bodies.insert(RigidBody::new(1.0));

        for _ in 0..10 {
            pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        }

        // At least some non-zero timing should have been recorded.
        let total = pipeline.stats.timings.total();
        assert!(total >= Duration::ZERO);
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
            .profiling(true)
            .finish();

        assert!((pipeline.config.gravity.y + 9.81).abs() < 1e-12);
        assert!(!pipeline.config.ccd_enabled);
        assert_eq!(pipeline.config.solver_iterations, 8_u32);
        assert!(pipeline.profiling_enabled);
    }

    #[test]
    fn test_pipeline_checkpoint_restore() {
        let mut pipeline = PhysicsPipeline::new();
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        let mut floor = RigidBody::new_static();
        floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let fh = bodies.insert(floor);
        let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
        colliders.insert(Collider::new(floor_shape).with_body(fh));

        let mut boxb = RigidBody::new(1.0);
        boxb.transform = Transform::from_position(Vec3::new(0.0, 0.48, 0.0));
        let bh = bodies.insert(boxb);
        let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
        colliders.insert(Collider::new(box_shape).with_body(bh));

        // Run a few steps and take a checkpoint.
        for _ in 0..5 {
            pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        }
        let checkpoint = pipeline.checkpoint();
        let time_at_checkpoint = checkpoint.time;

        // Run more steps.
        for _ in 0..10 {
            pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
        }
        assert!(pipeline.time > time_at_checkpoint);

        // Restore.
        pipeline.restore_checkpoint(checkpoint);
        assert!((pipeline.time - time_at_checkpoint).abs() < 1e-14);
    }

    #[test]
    fn test_stage_callback_fires() {
        let counter = Arc::new(AtomicUsize::new(0));
        let counter_clone = Arc::clone(&counter);

        let mut pipeline = PhysicsPipeline::new();
        pipeline.set_stage_callback(Box::new(move |_stage, _dur| {
            counter_clone.fetch_add(1, Ordering::Relaxed);
        }));

        let mut bodies = RigidBodySet::new();
        let colliders = ColliderSet::new();
        bodies.insert(RigidBody::new(1.0));

        pipeline.step(1.0 / 60.0, &mut bodies, &colliders);

        // We expect at least one callback (forces + sleep_check at minimum).
        assert!(counter.load(Ordering::Relaxed) >= 2);
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
    }

    #[test]
    fn test_sph_rigid_coupling_params_default() {
        let params = SphRigidCouplingParams::default();
        assert!((params.fluid_density - 1000.0).abs() < 1e-12);
        assert!((params.coupling_alpha - 1.0).abs() < 1e-12);
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

        let awake = count_awake_bodies(&bodies);
        assert_eq!(awake, 1, "only the dynamic body should be awake");
    }

    #[test]
    fn test_apply_buoyancy_impulse() {
        let coupling = SphRigidCouplingParams {
            fluid_density: 1000.0,
            coupling_alpha: 1.0,
            region_half_extents: Vec3::new(10.0, 10.0, 10.0),
            region_centre: Vec3::zeros(),
        };

        let mut bodies = RigidBodySet::new();
        let mut body = RigidBody::new(1.0);
        body.velocity = Vec3::zeros();
        body.transform = Transform::from_position(Vec3::zeros());
        let bh = bodies.insert(body);

        apply_buoyancy_impulse(&mut bodies, &coupling, 0.5, 9.81, 1.0 / 60.0);

        let vy = bodies.get(bh).unwrap().velocity.y;
        assert!(vy > 0.0, "buoyancy should push body upward, got vy={}", vy);
    }
}
