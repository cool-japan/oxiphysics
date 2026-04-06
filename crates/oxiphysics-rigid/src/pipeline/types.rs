//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
use std::time::Instant;

/// Broadphase AABB overlap detector: O(n²) all-pairs.
///
/// Returns a list of `(i, j)` index pairs where the bounding AABBs of bodies
/// `i` and `j` overlap.
pub struct BroadphaseDetector;
impl BroadphaseDetector {
    /// Detect all overlapping pairs in `bodies` using AABB tests.
    pub fn detect(bodies: &[BodySnapshot]) -> Vec<(usize, usize)> {
        let n = bodies.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if !bodies[i].active && !bodies[j].active {
                    continue;
                }
                let (mn_a, mx_a) = body_aabb(bodies[i].position, bodies[i].inv_mass);
                let (mn_b, mx_b) = body_aabb(bodies[j].position, bodies[j].inv_mass);
                if aabb_overlap(mn_a, mx_a, mn_b, mx_b) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
}
/// Interpolates body state between two physics sub-steps for smooth rendering.
///
/// The physics pipeline may run at a fixed timestep (e.g. 240 Hz) while the
/// renderer runs at a different rate (e.g. 60 Hz).  The interpolator blends
/// between the previous and current physics state using the sub-step
/// remainder `alpha ∈ [0, 1]`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SubStepInterpolator {
    /// Previous-frame body states.
    pub prev_states: Vec<BodySnapshot>,
    /// Current-frame body states.
    pub curr_states: Vec<BodySnapshot>,
}
impl SubStepInterpolator {
    /// Create a new interpolator with empty state buffers.
    pub fn new() -> Self {
        Self {
            prev_states: Vec::new(),
            curr_states: Vec::new(),
        }
    }
    /// Record the current states as `prev`, and update `curr`.
    pub fn advance(&mut self, new_states: Vec<BodySnapshot>) {
        self.prev_states = std::mem::replace(&mut self.curr_states, new_states);
    }
    /// Linearly interpolate position for body at index `i` with blend `alpha`.
    ///
    /// `alpha = 0` → previous state, `alpha = 1` → current state.
    ///
    /// Returns `None` if `i` is out of range in either buffer.
    pub fn interpolate_position(&self, i: usize, alpha: f64) -> Option<[f64; 3]> {
        let p = self.prev_states.get(i)?;
        let c = self.curr_states.get(i)?;
        let lerp = |a: f64, b: f64| a + (b - a) * alpha;
        Some([
            lerp(p.position[0], c.position[0]),
            lerp(p.position[1], c.position[1]),
            lerp(p.position[2], c.position[2]),
        ])
    }
    /// Linearly interpolate velocity for body at index `i`.
    pub fn interpolate_velocity(&self, i: usize, alpha: f64) -> Option<[f64; 3]> {
        let p = self.prev_states.get(i)?;
        let c = self.curr_states.get(i)?;
        let lerp = |a: f64, b: f64| a + (b - a) * alpha;
        Some([
            lerp(p.velocity[0], c.velocity[0]),
            lerp(p.velocity[1], c.velocity[1]),
            lerp(p.velocity[2], c.velocity[2]),
        ])
    }
    /// Number of bodies tracked.
    pub fn body_count(&self) -> usize {
        self.curr_states.len()
    }
}
/// Physics pipeline with fine-grained per-phase profiling.
///
/// Each call to [`ProfiledPipeline::step`] records [`PhaseTimings`] for
/// the most recent step and accumulates running totals.
pub struct ProfiledPipeline {
    /// Underlying pipeline configuration.
    pub config: PhysicsPipelineConfig,
    /// Sleep timer state (indexed like `bodies`).
    pub(super) sleep_timers: Vec<f64>,
    /// Timings from the most recent step.
    pub last_timings: PhaseTimings,
    /// Cumulative timings since creation or last [`ProfiledPipeline::reset_stats`].
    pub cumulative_timings: PhaseTimings,
    /// Number of profiled steps taken.
    pub step_count: u64,
}
impl ProfiledPipeline {
    /// Create a new `ProfiledPipeline` from config.
    pub fn new(config: PhysicsPipelineConfig) -> Self {
        Self {
            config,
            sleep_timers: Vec::new(),
            last_timings: PhaseTimings::default(),
            cumulative_timings: PhaseTimings::default(),
            step_count: 0,
        }
    }
    /// Reset accumulated statistics (does not affect simulation state).
    pub fn reset_stats(&mut self) {
        self.cumulative_timings = PhaseTimings::default();
        self.step_count = 0;
    }
    fn sync_timers(&mut self, n: usize) {
        if self.sleep_timers.len() < n {
            self.sleep_timers.resize(n, 0.0);
        }
    }
    /// Step the simulation, measuring per-phase time.
    pub fn step(&mut self, bodies: &mut Vec<BodySnapshot>) -> StepReport {
        let step_start = Instant::now();
        let dt = self.config.dt / self.config.sub_steps as f64;
        self.sync_timers(bodies.len());
        let mut report = StepReport::default();
        let mut timings = PhaseTimings::default();
        for _ in 0..self.config.sub_steps {
            let step = PipelineStep::new(PhysicsPipelineConfig {
                dt,
                ..self.config.clone()
            });
            let t0 = Instant::now();
            step.apply_gravity(bodies);
            timings.gravity_us += t0.elapsed().as_micros() as u64;
            let t1 = Instant::now();
            let pairs = BroadphaseDetector::detect(bodies);
            timings.broadphase_us += t1.elapsed().as_micros() as u64;
            report.num_collisions += pairs.len();
            let t2 = Instant::now();
            let contacts = NarrowphaseDetector::generate_contacts(bodies, &pairs);
            timings.narrowphase_us += t2.elapsed().as_micros() as u64;
            report.num_contacts += contacts.len();
            let t3 = Instant::now();
            let solver = SequentialImpulseSolver::new(self.config.solver_iterations, 0.3);
            solver.solve(bodies, &contacts);
            timings.solver_us += t3.elapsed().as_micros() as u64;
            let t4 = Instant::now();
            step.integrate_positions(bodies);
            timings.integrate_us += t4.elapsed().as_micros() as u64;
            let t5 = Instant::now();
            let sleep_mgr = SleepManager::new(0.01, 0.01, 0.5);
            sleep_mgr.update(bodies, &mut self.sleep_timers, dt);
            timings.sleep_us += t5.elapsed().as_micros() as u64;
        }
        timings.total_us = step_start.elapsed().as_micros() as u64;
        report.num_sleeping = SleepManager::sleeping_count(bodies);
        report.step_time_us = timings.total_us;
        self.cumulative_timings.gravity_us += timings.gravity_us;
        self.cumulative_timings.broadphase_us += timings.broadphase_us;
        self.cumulative_timings.narrowphase_us += timings.narrowphase_us;
        self.cumulative_timings.solver_us += timings.solver_us;
        self.cumulative_timings.integrate_us += timings.integrate_us;
        self.cumulative_timings.sleep_us += timings.sleep_us;
        self.cumulative_timings.total_us += timings.total_us;
        self.step_count += 1;
        self.last_timings = timings;
        report
    }
    /// Average step time in microseconds over all steps since last reset.
    pub fn avg_step_us(&self) -> f64 {
        if self.step_count == 0 {
            0.0
        } else {
            self.cumulative_timings.total_us as f64 / self.step_count as f64
        }
    }
}
/// Event queue that accumulates [`PhysicsEvent`]s during a simulation step.
///
/// After calling [`EventAwarePipeline::step`] the caller drains `events` to
/// process collision begin/end and sleep/wake notifications.
#[derive(Debug, Clone, Default)]
pub struct EventQueue {
    /// Events accumulated since the last drain.
    pub events: Vec<PhysicsEvent>,
}
impl EventQueue {
    /// Create a new empty event queue.
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }
    /// Push a new event.
    pub fn push(&mut self, e: PhysicsEvent) {
        self.events.push(e);
    }
    /// Drain (consume and return) all events.
    pub fn drain(&mut self) -> Vec<PhysicsEvent> {
        std::mem::take(&mut self.events)
    }
    /// Number of pending events.
    pub fn len(&self) -> usize {
        self.events.len()
    }
    /// Returns `true` if there are no pending events.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
    /// Count events of a given kind.
    pub fn count_kind(&self, kind: PhysicsEventKind) -> usize {
        self.events.iter().filter(|e| e.kind == kind).count()
    }
}
/// A constant external force that is applied to all dynamic bodies every
/// sub-step (e.g. a wind force).
pub struct ConstantForceHook {
    /// Force vector in world space (N).
    pub force: [f64; 3],
}
/// A velocity-level contact constraint between two bodies.
///
/// Used in the velocity-impulse step: each contact constraint enforces
/// non-penetration via normal impulse and Coulomb friction via tangent
/// impulses.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VelocityContactConstraint {
    /// Index of body A in the snapshot array.
    pub body_a: usize,
    /// Index of body B in the snapshot array (`usize::MAX` for static).
    pub body_b: usize,
    /// Contact normal pointing from B to A.
    pub normal: [f64; 3],
    /// Relative velocity along the normal before impulse application.
    pub v_rel_n: f64,
    /// Penetration depth (negative = overlap).
    pub penetration: f64,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Coefficient of friction.
    pub friction: f64,
    /// Effective mass along the normal direction.
    pub effective_mass: f64,
    /// Accumulated normal impulse (for clamping).
    pub lambda_n: f64,
}
impl VelocityContactConstraint {
    /// Create a new velocity contact constraint.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        body_a: usize,
        body_b: usize,
        normal: [f64; 3],
        v_rel_n: f64,
        penetration: f64,
        restitution: f64,
        friction: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
    ) -> Self {
        let effective_mass = 1.0 / (inv_mass_a + inv_mass_b).max(1e-30);
        Self {
            body_a,
            body_b,
            normal,
            v_rel_n,
            penetration,
            restitution,
            friction,
            effective_mass,
            lambda_n: 0.0,
        }
    }
    /// Compute the Baumgarte velocity bias for position stabilization.
    ///
    /// ```text
    /// b = β/Δt · max(0, −d − slop)
    /// ```
    ///
    /// where `β` is the Baumgarte coefficient (typically 0.2) and `slop` is
    /// the allowed penetration (typically 0.005 m).
    pub fn baumgarte_bias(&self, dt: f64, beta: f64, slop: f64) -> f64 {
        let d = self.penetration;
        (beta / dt) * ((-d - slop).max(0.0))
    }
    /// Compute the normal impulse magnitude for one solver iteration.
    ///
    /// Returns the change in accumulated impulse `Δλ` (clamped).
    pub fn solve_normal(&mut self, v_rel_n_current: f64, bias: f64) -> f64 {
        let target_velocity = -self.restitution * self.v_rel_n.min(0.0);
        let delta_v = target_velocity - v_rel_n_current + bias;
        let delta_lambda = self.effective_mass * delta_v;
        let old_lambda = self.lambda_n;
        self.lambda_n = (self.lambda_n + delta_lambda).max(0.0);
        self.lambda_n - old_lambda
    }
    /// Maximum friction impulse magnitude (Coulomb cone).
    pub fn max_friction_impulse(&self) -> f64 {
        self.friction * self.lambda_n.abs()
    }
    /// Returns `true` if the contact is in the Coulomb friction cone.
    pub fn in_friction_cone(&self, tangent_impulse: f64) -> bool {
        tangent_impulse.abs() <= self.max_friction_impulse()
    }
}
/// Orchestrates a complete simulation step:
///
/// 1. Apply gravity → update velocities.
/// 2. Broadphase AABB overlap detection.
/// 3. Narrowphase contact generation (sphere–sphere).
/// 4. Sequential-impulse constraint solving.
/// 5. CCD TOI handling (advance to earliest TOI, resolve, then continue).
/// 6. Integrate positions.
/// 7. Island-based sleeping.
pub struct PhysicsPipeline {
    /// Pipeline configuration.
    pub config: PhysicsPipelineConfig,
    /// Per-body sleep timer state (indexed like `bodies` slice).
    pub(super) sleep_timers: Vec<f64>,
}
impl PhysicsPipeline {
    /// Creates a new `PhysicsPipeline` from the given config.
    pub fn new(config: PhysicsPipelineConfig) -> Self {
        Self {
            config,
            sleep_timers: Vec::new(),
        }
    }
    /// Ensures the internal sleep timer buffer matches the length of `bodies`.
    fn sync_timers(&mut self, n: usize) {
        if self.sleep_timers.len() < n {
            self.sleep_timers.resize(n, 0.0);
        }
    }
    /// Executes one full simulation step over `bodies` and returns a
    /// [`StepReport`].
    ///
    /// # Arguments
    /// * `bodies` — slice of body snapshots (mutated in-place).
    pub fn step(&mut self, bodies: &mut Vec<BodySnapshot>) -> StepReport {
        let start = Instant::now();
        let dt = self.config.dt / self.config.sub_steps as f64;
        self.sync_timers(bodies.len());
        let mut report = StepReport::default();
        for _ in 0..self.config.sub_steps {
            let step = PipelineStep::new(PhysicsPipelineConfig {
                dt,
                ..self.config.clone()
            });
            step.apply_gravity(bodies);
            let pairs = BroadphaseDetector::detect(bodies);
            report.num_collisions += pairs.len();
            let contacts = NarrowphaseDetector::generate_contacts(bodies, &pairs);
            report.num_contacts += contacts.len();
            if !contacts.is_empty() {
                self.handle_ccd(bodies, &contacts, dt);
            }
            let solver = SequentialImpulseSolver::new(self.config.solver_iterations, 0.3);
            solver.solve(bodies, &contacts);
            step.integrate_positions(bodies);
            let sleep_mgr = SleepManager::new(0.01, 0.01, 0.5);
            sleep_mgr.update(bodies, &mut self.sleep_timers, dt);
        }
        report.num_sleeping = SleepManager::sleeping_count(bodies);
        report.step_time_us = start.elapsed().as_micros() as u64;
        report
    }
    /// Minimal CCD: find the earliest contact (by penetration depth) and apply
    /// a separating impulse before the solver runs.
    fn handle_ccd(&self, bodies: &mut Vec<BodySnapshot>, contacts: &[Contact], _dt: f64) {
        let deepest = contacts.iter().max_by(|a, b| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(c) = deepest {
            let n = c.normal;
            let ia = c.idx_a;
            let ib = c.idx_b;
            let inv_ma = bodies[ia].inv_mass;
            let inv_mb = bodies[ib].inv_mass;
            let va = bodies[ia].velocity;
            let vb = bodies[ib].velocity;
            let rel_v = sub3(va, vb);
            let vn = dot3(rel_v, n);
            if vn < 0.0 {
                let denom = inv_ma + inv_mb;
                if denom > 1e-30 {
                    let j = -vn / denom;
                    if inv_ma > 0.0 {
                        let dv = scale3(n, j * inv_ma);
                        bodies[ia].velocity = add3(bodies[ia].velocity, dv);
                    }
                    if inv_mb > 0.0 {
                        let dv = scale3(n, j * inv_mb);
                        bodies[ib].velocity = sub3(bodies[ib].velocity, dv);
                    }
                }
            }
        }
    }
}
/// Executes a single physics pipeline step.
pub struct PipelineStep {
    /// Pipeline configuration.
    pub config: PhysicsPipelineConfig,
}
impl PipelineStep {
    /// Creates a new `PipelineStep` with the given config.
    pub fn new(config: PhysicsPipelineConfig) -> Self {
        Self { config }
    }
    /// Applies gravity to all active dynamic bodies (inv_mass > 0).
    pub fn apply_gravity(&self, bodies: &mut [BodySnapshot]) {
        let g = self.config.gravity;
        let dt = self.config.dt;
        for body in bodies.iter_mut() {
            if body.active && body.inv_mass > 0.0 {
                body.velocity[0] += g[0] * dt;
                body.velocity[1] += g[1] * dt;
                body.velocity[2] += g[2] * dt;
            }
        }
    }
    /// Integrates positions using semi-implicit Euler: position += velocity * dt.
    pub fn integrate_positions(&self, bodies: &mut [BodySnapshot]) {
        let dt = self.config.dt;
        for body in bodies.iter_mut() {
            if body.active {
                body.position[0] += body.velocity[0] * dt;
                body.position[1] += body.velocity[1] * dt;
                body.position[2] += body.velocity[2] * dt;
            }
        }
    }
    /// Runs one full pipeline step: apply gravity, then integrate positions.
    pub fn step(&self, bodies: &mut [BodySnapshot]) {
        self.apply_gravity(bodies);
        self.integrate_positions(bodies);
    }
}
/// A contact manifold manages a set of [`VelocityContactConstraint`]s.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ContactManifold {
    /// All contact constraints in this manifold.
    pub contacts: Vec<VelocityContactConstraint>,
}
impl ContactManifold {
    /// Create a new empty manifold.
    pub fn new() -> Self {
        Self {
            contacts: Vec::new(),
        }
    }
    /// Add a contact.
    pub fn add(&mut self, c: VelocityContactConstraint) {
        self.contacts.push(c);
    }
    /// Clear all contacts.
    pub fn clear(&mut self) {
        self.contacts.clear();
    }
    /// Number of contacts.
    pub fn len(&self) -> usize {
        self.contacts.len()
    }
    /// Returns `true` if there are no contacts.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }
    /// Total accumulated normal impulse across all contacts.
    pub fn total_normal_impulse(&self) -> f64 {
        self.contacts.iter().map(|c| c.lambda_n).sum()
    }
    /// Maximum penetration depth across all contacts.
    pub fn max_penetration(&self) -> f64 {
        self.contacts
            .iter()
            .map(|c| -c.penetration)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
/// Data for a single simulation island (a set of connected bodies).
#[derive(Debug, Clone)]
pub struct IslandData {
    /// Body IDs belonging to this island.
    pub body_ids: Vec<u64>,
    /// Contact pairs within this island.
    pub contact_pairs: Vec<(u64, u64)>,
    /// Whether this island is active (non-sleeping).
    pub active: bool,
}
/// Cache of warm-start entries indexed by (body_a, body_b) pair.
///
/// The cache is rebuilt each frame by matching current contacts against the
/// previous frame's contact set.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct WarmStartCache {
    /// All cached entries from the previous frame.
    pub entries: Vec<WarmStartEntry>,
}
impl WarmStartCache {
    /// Create a new empty warm-start cache.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }
    /// Look up the warm-start entry for a given body pair (order-independent).
    pub fn find(&self, a: u64, b: u64) -> Option<&WarmStartEntry> {
        self.entries
            .iter()
            .find(|e| (e.body_a == a && e.body_b == b) || (e.body_a == b && e.body_b == a))
    }
    /// Insert or update the entry for the given body pair.
    pub fn upsert(&mut self, entry: WarmStartEntry) {
        let a = entry.body_a;
        let b = entry.body_b;
        if let Some(e) = self
            .entries
            .iter_mut()
            .find(|e| (e.body_a == a && e.body_b == b) || (e.body_a == b && e.body_b == a))
        {
            *e = entry;
        } else {
            self.entries.push(entry);
        }
    }
    /// Remove entries for bodies that are no longer in contact.
    ///
    /// `active_pairs` is the set of (a, b) pairs present in the current
    /// contact manifold.
    pub fn prune(&mut self, active_pairs: &[(u64, u64)]) {
        self.entries.retain(|e| {
            active_pairs
                .iter()
                .any(|&(a, b)| (e.body_a == a && e.body_b == b) || (e.body_a == b && e.body_b == a))
        });
    }
    /// Clear all cached entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Scale all cached impulses by `factor`.
    pub fn scale_all(&mut self, factor: f64) {
        for e in &mut self.entries {
            e.scale(factor);
        }
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// Drives a physics pipeline at a fixed simulation timestep regardless of
/// the variable render frame rate.
///
/// A time accumulator absorbs the variable render `dt` and fires one physics
/// step for every complete `fixed_dt` that accumulates.  Any leftover time
/// carries forward to the next frame.
///
/// # Usage
/// ```ignore
/// let mut stepper = DeterministicStepper::new(1.0 / 240.0, 8);
/// loop {
///     let render_dt = measure_frame_time();
///     stepper.advance(render_dt, &mut bodies, &mut pipeline);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct DeterministicStepper {
    /// Fixed simulation timestep (seconds).
    pub fixed_dt: f64,
    /// Maximum number of physics steps allowed per rendered frame (prevents
    /// spiral-of-death when the simulation falls behind).
    pub max_steps_per_frame: usize,
    /// Accumulated time not yet consumed by simulation steps.
    pub accumulator: f64,
    /// Total simulation time elapsed so far (seconds).
    pub sim_time: f64,
    /// Total number of physics steps executed.
    pub total_steps: u64,
}
impl DeterministicStepper {
    /// Create a new stepper with the given fixed timestep and per-frame step cap.
    pub fn new(fixed_dt: f64, max_steps_per_frame: usize) -> Self {
        Self {
            fixed_dt,
            max_steps_per_frame,
            accumulator: 0.0,
            sim_time: 0.0,
            total_steps: 0,
        }
    }
    /// Feed `render_dt` seconds into the accumulator and fire physics steps
    /// on `pipeline` / `bodies` until the accumulator is drained.
    ///
    /// Returns the number of sub-steps taken this frame and the interpolation
    /// remainder `alpha ∈ [0, 1)` for rendering.
    pub fn advance(
        &mut self,
        render_dt: f64,
        bodies: &mut Vec<BodySnapshot>,
        pipeline: &mut PhysicsPipeline,
    ) -> (usize, f64) {
        self.accumulator += render_dt;
        let mut steps_taken = 0;
        while self.accumulator >= self.fixed_dt && steps_taken < self.max_steps_per_frame {
            let old_dt = pipeline.config.dt;
            pipeline.config.dt = self.fixed_dt;
            pipeline.step(bodies);
            pipeline.config.dt = old_dt;
            self.accumulator -= self.fixed_dt;
            self.sim_time += self.fixed_dt;
            self.total_steps += 1;
            steps_taken += 1;
        }
        let alpha = self.accumulator / self.fixed_dt;
        (steps_taken, alpha.min(1.0))
    }
    /// Reset the accumulator and simulation time.
    pub fn reset(&mut self) {
        self.accumulator = 0.0;
        self.sim_time = 0.0;
        self.total_steps = 0;
    }
    /// Returns `true` if the stepper has accumulated enough time to take at
    /// least one step.
    pub fn is_ready(&self) -> bool {
        self.accumulator >= self.fixed_dt
    }
}
/// A linear velocity damping hook (drag), applied as `v *= 1 − coeff * dt`.
pub struct DampingHook {
    /// Damping coefficient (0 = no damping, 1 = full stop in 1 s).
    pub coeff: f64,
}
/// A hook that records the total kinetic energy after each sub-step.
pub struct EnergyMonitorHook {
    /// Kinetic energy samples, one per sub-step (appended).
    pub samples: std::sync::Mutex<Vec<f64>>,
}
impl EnergyMonitorHook {
    /// Create a new energy monitor.
    pub fn new() -> Self {
        Self {
            samples: std::sync::Mutex::new(Vec::new()),
        }
    }
    /// All recorded kinetic energy samples.
    pub fn snapshot(&self) -> Vec<f64> {
        self.samples
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .clone()
    }
}
/// Cached impulse for a contact pair, used to warm-start the sequential
/// impulse solver on the next frame.
///
/// Warm-starting reuses the previous frame's accumulated impulse as the
/// initial guess, which reduces the number of solver iterations needed for
/// convergence in stacked/persistent contacts.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WarmStartEntry {
    /// First body ID.
    pub body_a: u64,
    /// Second body ID.
    pub body_b: u64,
    /// Contact normal (world frame).
    pub normal: [f64; 3],
    /// Accumulated normal impulse from the previous frame.
    pub lambda_n: f64,
    /// Accumulated tangent impulse (friction), component 1.
    pub lambda_t1: f64,
    /// Accumulated tangent impulse (friction), component 2.
    pub lambda_t2: f64,
}
impl WarmStartEntry {
    /// Create a zero-initialized warm-start entry.
    pub fn new(body_a: u64, body_b: u64, normal: [f64; 3]) -> Self {
        Self {
            body_a,
            body_b,
            normal,
            lambda_n: 0.0,
            lambda_t1: 0.0,
            lambda_t2: 0.0,
        }
    }
    /// Scale all impulses by `factor` (e.g. for timestep ratio adjustment).
    pub fn scale(&mut self, factor: f64) {
        self.lambda_n *= factor;
        self.lambda_t1 *= factor;
        self.lambda_t2 *= factor;
    }
    /// Clamp normal impulse to be non-negative (contacts only push, not pull).
    pub fn clamp_normal(&mut self) {
        if self.lambda_n < 0.0 {
            self.lambda_n = 0.0;
        }
    }
}
/// Narrowphase contact generator: sphere–sphere, sphere–box, box–box.
///
/// Every body is represented as a sphere whose radius is derived from its
/// mass.  This provides a unified contact generation pass without requiring
/// explicit shape metadata in `BodySnapshot`.
pub struct NarrowphaseDetector;
impl NarrowphaseDetector {
    /// Generate contacts for each candidate pair `(i, j)` in `pairs`.
    pub fn generate_contacts(bodies: &[BodySnapshot], pairs: &[(usize, usize)]) -> Vec<Contact> {
        let mut contacts = Vec::new();
        for &(i, j) in pairs {
            if let Some(c) = Self::sphere_sphere(bodies, i, j) {
                contacts.push(c);
            }
        }
        contacts
    }
    /// Sphere–sphere contact test between bodies at indices `i` and `j`.
    fn sphere_sphere(bodies: &[BodySnapshot], i: usize, j: usize) -> Option<Contact> {
        let ba = &bodies[i];
        let bb = &bodies[j];
        let ra = Self::bounding_radius(ba.inv_mass);
        let rb = Self::bounding_radius(bb.inv_mass);
        let diff = sub3(ba.position, bb.position);
        let dist_sq = dot3(diff, diff);
        let sum_r = ra + rb;
        if dist_sq >= sum_r * sum_r {
            return None;
        }
        let dist = dist_sq.sqrt();
        let depth = sum_r - dist;
        let normal = normalize3(diff);
        let contact_point = add3(bb.position, scale3(normal, rb));
        Some(Contact {
            idx_a: i,
            idx_b: j,
            normal,
            depth,
            contact_point,
        })
    }
    /// Bounding sphere radius heuristic from inverse mass.
    fn bounding_radius(inv_mass: f64) -> f64 {
        if inv_mass > 0.0 {
            (1.0 / inv_mass).cbrt() * 0.3
        } else {
            0.5_f64
        }
        .max(0.1)
    }
}
/// Sequential-impulse constraint solver: resolves velocity-level contact
/// constraints for a list of [`Contact`]s over a slice of [`BodySnapshot`]s.
pub struct SequentialImpulseSolver {
    /// Number of Gauss-Seidel iterations.
    pub iterations: usize,
    /// Coefficient of restitution (bounciness, 0..1).
    pub restitution: f64,
}
impl SequentialImpulseSolver {
    /// Creates a solver with the given iteration count and restitution.
    pub fn new(iterations: usize, restitution: f64) -> Self {
        Self {
            iterations,
            restitution,
        }
    }
    /// Resolves all contacts in-place.
    pub fn solve(&self, bodies: &mut [BodySnapshot], contacts: &[Contact]) {
        for _ in 0..self.iterations {
            for c in contacts {
                self.resolve_one(bodies, c);
            }
        }
    }
    fn resolve_one(&self, bodies: &mut [BodySnapshot], c: &Contact) {
        let n = c.normal;
        let ia = c.idx_a;
        let ib = c.idx_b;
        let va = bodies[ia].velocity;
        let vb = bodies[ib].velocity;
        let inv_ma = bodies[ia].inv_mass;
        let inv_mb = bodies[ib].inv_mass;
        let rel_v = sub3(va, vb);
        let vn = dot3(rel_v, n);
        if vn >= 0.0 {
            return;
        }
        let denom = inv_ma + inv_mb;
        if denom < 1e-30 {
            return;
        }
        let j = -(1.0 + self.restitution) * vn / denom;
        if inv_ma > 0.0 {
            let dv = scale3(n, j * inv_ma);
            bodies[ia].velocity = add3(bodies[ia].velocity, dv);
        }
        if inv_mb > 0.0 {
            let dv = scale3(n, j * inv_mb);
            bodies[ib].velocity = sub3(bodies[ib].velocity, dv);
        }
    }
}
/// Manages body sleep state based on velocity thresholds.
pub struct SleepManager {
    /// Linear velocity threshold (m/s).
    pub linear_threshold: f64,
    /// Angular velocity threshold (rad/s).
    pub angular_threshold: f64,
    /// Time below threshold required to sleep (seconds).
    pub time_before_sleep: f64,
}
impl SleepManager {
    /// Creates a new `SleepManager` with the given thresholds.
    pub fn new(linear_threshold: f64, angular_threshold: f64, time_before_sleep: f64) -> Self {
        Self {
            linear_threshold,
            angular_threshold,
            time_before_sleep,
        }
    }
    /// Advances sleep timers and sets `active = false` for bodies that have
    /// been below threshold long enough.
    ///
    /// `timers` must have the same length as `bodies`.
    pub fn update(&self, bodies: &mut [BodySnapshot], timers: &mut [f64], dt: f64) {
        for (body, timer) in bodies.iter_mut().zip(timers.iter_mut()) {
            if body.inv_mass == 0.0 {
                continue;
            }
            let v_sq = dot3(body.velocity, body.velocity);
            let w_sq = dot3(body.ang_velocity, body.ang_velocity);
            let lt = self.linear_threshold;
            let at = self.angular_threshold;
            if v_sq < lt * lt && w_sq < at * at {
                *timer += dt;
                if *timer >= self.time_before_sleep {
                    body.active = false;
                }
            } else {
                *timer = 0.0;
                body.active = true;
            }
        }
    }
    /// Returns the number of sleeping bodies (`active == false` and dynamic).
    pub fn sleeping_count(bodies: &[BodySnapshot]) -> usize {
        bodies
            .iter()
            .filter(|b| !b.active && b.inv_mass > 0.0)
            .count()
    }
}
/// Configuration for the physics pipeline.
#[derive(Debug, Clone)]
pub struct PhysicsPipelineConfig {
    /// Gravitational acceleration vector \[x, y, z\].
    pub gravity: [f64; 3],
    /// Time step in seconds.
    pub dt: f64,
    /// Number of sub-steps per frame.
    pub sub_steps: usize,
    /// Number of solver iterations per sub-step.
    pub solver_iterations: usize,
}
impl PhysicsPipelineConfig {
    /// Creates a `PhysicsPipelineConfig` with physics defaults:
    /// gravity = \[0, -9.81, 0\], dt = 1/60, sub_steps = 1, solver_iterations = 10.
    pub fn new() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            sub_steps: 1,
            solver_iterations: 10,
        }
    }
}
/// A snapshot of a rigid body's state for pipeline processing.
#[derive(Debug, Clone)]
pub struct BodySnapshot {
    /// Unique body identifier.
    pub id: u64,
    /// World-space position.
    pub position: [f64; 3],
    /// Linear velocity.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub ang_velocity: [f64; 3],
    /// Inverse mass (0 = static/kinematic).
    pub inv_mass: f64,
    /// Whether the body is actively simulated.
    pub active: bool,
}
/// Detailed timing breakdown for one physics pipeline step.
///
/// Each field measures the wall-clock time in *microseconds* spent in the
/// corresponding pipeline phase.
#[derive(Debug, Clone, Default)]
pub struct PhaseTimings {
    /// Time spent applying gravity to all bodies.
    pub gravity_us: u64,
    /// Time spent in broadphase AABB detection.
    pub broadphase_us: u64,
    /// Time spent in narrowphase contact generation.
    pub narrowphase_us: u64,
    /// Time spent running the sequential-impulse solver.
    pub solver_us: u64,
    /// Time spent integrating positions.
    pub integrate_us: u64,
    /// Time spent in sleep management.
    pub sleep_us: u64,
    /// Total wall-clock time for the step.
    pub total_us: u64,
}
impl PhaseTimings {
    /// Returns the fraction of time spent in the solver vs. total step time.
    pub fn solver_fraction(&self) -> f64 {
        if self.total_us == 0 {
            0.0
        } else {
            self.solver_us as f64 / self.total_us as f64
        }
    }
    /// Returns the fraction of time spent in broadphase + narrowphase.
    pub fn collision_fraction(&self) -> f64 {
        if self.total_us == 0 {
            0.0
        } else {
            (self.broadphase_us + self.narrowphase_us) as f64 / self.total_us as f64
        }
    }
}
/// The type of a physics event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicsEventKind {
    /// Two bodies began to overlap (first contact).
    CollisionBegin,
    /// Two bodies ceased to overlap (contact lost).
    CollisionEnd,
    /// A body entered sleep state.
    BodySlept,
    /// A body was woken from sleep.
    BodyWoke,
}
/// A physics pipeline that emits [`PhysicsEvent`]s for collision begin/end
/// and body sleep/wake transitions.
pub struct EventAwarePipeline {
    /// Pipeline configuration.
    pub config: PhysicsPipelineConfig,
    /// Per-body sleep timers.
    pub(super) sleep_timers: Vec<f64>,
    /// Previous-frame set of overlapping body-ID pairs for delta detection.
    pub(super) prev_contact_pairs: std::collections::HashSet<(u64, u64)>,
    /// Previous-frame set of sleeping body IDs for sleep/wake detection.
    pub(super) prev_sleeping: std::collections::HashSet<u64>,
    /// Accumulated simulation time (seconds).
    pub sim_time: f64,
}
impl EventAwarePipeline {
    /// Create a new `EventAwarePipeline`.
    pub fn new(config: PhysicsPipelineConfig) -> Self {
        Self {
            config,
            sleep_timers: Vec::new(),
            prev_contact_pairs: std::collections::HashSet::new(),
            prev_sleeping: std::collections::HashSet::new(),
            sim_time: 0.0,
        }
    }
    fn sync_timers(&mut self, n: usize) {
        if self.sleep_timers.len() < n {
            self.sleep_timers.resize(n, 0.0);
        }
    }
    /// Run one simulation step and collect events into `queue`.
    pub fn step(&mut self, bodies: &mut Vec<BodySnapshot>, queue: &mut EventQueue) -> StepReport {
        let start = Instant::now();
        let dt = self.config.dt / self.config.sub_steps as f64;
        self.sync_timers(bodies.len());
        let mut report = StepReport::default();
        let mut current_contacts: std::collections::HashSet<(u64, u64)> =
            std::collections::HashSet::new();
        for _ in 0..self.config.sub_steps {
            let step_cfg = PipelineStep::new(PhysicsPipelineConfig {
                dt,
                ..self.config.clone()
            });
            step_cfg.apply_gravity(bodies);
            let pairs = BroadphaseDetector::detect(bodies);
            report.num_collisions += pairs.len();
            let contacts = NarrowphaseDetector::generate_contacts(bodies, &pairs);
            report.num_contacts += contacts.len();
            for c in &contacts {
                let ia = bodies[c.idx_a].id;
                let ib = bodies[c.idx_b].id;
                let key = if ia <= ib { (ia, ib) } else { (ib, ia) };
                current_contacts.insert(key);
            }
            let solver = SequentialImpulseSolver::new(self.config.solver_iterations, 0.3);
            solver.solve(bodies, &contacts);
            step_cfg.integrate_positions(bodies);
            let sleep_mgr = SleepManager::new(0.01, 0.01, 0.5);
            sleep_mgr.update(bodies, &mut self.sleep_timers, dt);
        }
        self.sim_time += self.config.dt;
        let t = self.sim_time;
        for &(a, b) in &current_contacts {
            if !self.prev_contact_pairs.contains(&(a, b)) {
                queue.push(PhysicsEvent::collision(
                    PhysicsEventKind::CollisionBegin,
                    a,
                    b,
                    t,
                ));
            }
        }
        for &(a, b) in &self.prev_contact_pairs {
            if !current_contacts.contains(&(a, b)) {
                queue.push(PhysicsEvent::collision(
                    PhysicsEventKind::CollisionEnd,
                    a,
                    b,
                    t,
                ));
            }
        }
        let current_sleeping: std::collections::HashSet<u64> = bodies
            .iter()
            .filter(|b| !b.active && b.inv_mass > 0.0)
            .map(|b| b.id)
            .collect();
        for &id in &current_sleeping {
            if !self.prev_sleeping.contains(&id) {
                queue.push(PhysicsEvent::body(PhysicsEventKind::BodySlept, id, t));
            }
        }
        for &id in &self.prev_sleeping {
            if !current_sleeping.contains(&id) {
                queue.push(PhysicsEvent::body(PhysicsEventKind::BodyWoke, id, t));
            }
        }
        self.prev_contact_pairs = current_contacts;
        self.prev_sleeping = current_sleeping;
        report.num_sleeping = SleepManager::sleeping_count(bodies);
        report.step_time_us = start.elapsed().as_micros() as u64;
        report
    }
}
/// A single physics event carrying the event kind and the affected body IDs.
#[derive(Debug, Clone)]
pub struct PhysicsEvent {
    /// What happened.
    pub kind: PhysicsEventKind,
    /// Primary body ID.
    pub body_a: u64,
    /// Secondary body ID (for pair events; `u64::MAX` for single-body events).
    pub body_b: u64,
    /// Simulation time at which the event occurred.
    pub sim_time: f64,
}
impl PhysicsEvent {
    /// Construct a collision event.
    pub fn collision(kind: PhysicsEventKind, a: u64, b: u64, t: f64) -> Self {
        Self {
            kind,
            body_a: a,
            body_b: b,
            sim_time: t,
        }
    }
    /// Construct a single-body event.
    pub fn body(kind: PhysicsEventKind, id: u64, t: f64) -> Self {
        Self {
            kind,
            body_a: id,
            body_b: u64::MAX,
            sim_time: t,
        }
    }
}
/// Manages simulation islands built via union-find.
pub struct IslandManager {
    /// All islands produced by the last build.
    pub islands: Vec<IslandData>,
}
impl IslandManager {
    /// Creates an empty `IslandManager`.
    pub fn new() -> Self {
        Self {
            islands: Vec::new(),
        }
    }
    /// Groups `body_ids` into connected islands using `contact_pairs` as edges.
    pub fn build_islands(body_ids: &[u64], contact_pairs: &[(u64, u64)]) -> Self {
        let n = body_ids.len();
        if n == 0 {
            return Self {
                islands: Vec::new(),
            };
        }
        let mut id_to_idx = std::collections::HashMap::new();
        for (i, &id) in body_ids.iter().enumerate() {
            id_to_idx.insert(id, i);
        }
        let edges: Vec<(usize, usize)> = contact_pairs
            .iter()
            .filter_map(|(a, b)| {
                let ia = id_to_idx.get(a)?;
                let ib = id_to_idx.get(b)?;
                Some((*ia, *ib))
            })
            .collect();
        let mut parent = union_find_build(n, &edges);
        let roots: Vec<usize> = (0..n).map(|i| union_find_root(&mut parent, i)).collect();
        let mut root_to_island: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        let mut islands: Vec<IslandData> = Vec::new();
        for (i, &root) in roots.iter().enumerate() {
            let island_idx = root_to_island.entry(root).or_insert_with(|| {
                let idx = islands.len();
                islands.push(IslandData {
                    body_ids: Vec::new(),
                    contact_pairs: Vec::new(),
                    active: true,
                });
                idx
            });
            islands[*island_idx].body_ids.push(body_ids[i]);
        }
        for &(a, b) in contact_pairs {
            if let (Some(&ia), Some(&ib)) = (id_to_idx.get(&a), id_to_idx.get(&b)) {
                let root = union_find_root(&mut parent, ia);
                if let Some(&island_idx) = root_to_island.get(&root) {
                    let root_b = union_find_root(&mut parent, ib);
                    if root == root_b {
                        islands[island_idx].contact_pairs.push((a, b));
                    }
                }
            }
        }
        Self { islands }
    }
    /// Merges island `b` into island `a`.
    pub fn merge_islands(&mut self, a: usize, b: usize) {
        if a == b || b >= self.islands.len() || a >= self.islands.len() {
            return;
        }
        let island_b = self.islands.remove(b);
        let target = if a > b { a - 1 } else { a };
        self.islands[target].body_ids.extend(island_b.body_ids);
        self.islands[target]
            .contact_pairs
            .extend(island_b.contact_pairs);
    }
    /// Returns the number of islands.
    pub fn island_count(&self) -> usize {
        self.islands.len()
    }
    /// Returns the body IDs in island `i`.
    pub fn bodies_in_island(&self, i: usize) -> &[u64] {
        &self.islands[i].body_ids
    }
}
/// Statistics produced by one [`PhysicsPipeline::step`] call.
#[derive(Debug, Clone, Default)]
pub struct StepReport {
    /// Number of contact pairs generated by the narrowphase.
    pub num_contacts: usize,
    /// Number of broadphase collision candidate pairs.
    pub num_collisions: usize,
    /// Number of sleeping bodies at the end of the step.
    pub num_sleeping: usize,
    /// Wall-clock time of the step in microseconds.
    pub step_time_us: u64,
}
/// A contact generated by the narrowphase.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Index of the first body in the slice.
    pub idx_a: usize,
    /// Index of the second body in the slice.
    pub idx_b: usize,
    /// Contact normal (from B toward A, unit length).
    pub normal: [f64; 3],
    /// Penetration depth (metres, positive = overlap).
    pub depth: f64,
    /// World-space contact point.
    pub contact_point: [f64; 3],
}
