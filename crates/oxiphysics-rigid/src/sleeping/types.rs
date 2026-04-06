//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg, clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
/// A single rigid body tracked by the sleep manager.
#[derive(Debug, Clone)]
pub struct SleepingBody {
    /// Application-defined body identifier.
    pub body_id: u32,
    /// Current sleep lifecycle state.
    pub state: SleepState,
    /// Linear velocity in world space (m/s).
    pub linear_velocity: [f64; 3],
    /// Angular velocity in world space (rad/s).
    pub angular_velocity: [f64; 3],
    /// World-space position.
    pub position: [f64; 3],
    /// World-space orientation as a unit quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// Accumulated sleep timer in seconds.
    pub sleep_timer: f64,
    /// Mass (kg), used for energy-based sleeping.
    pub mass: f64,
    /// Scalar moment of inertia (kg*m^2), used for energy-based sleeping.
    pub inertia: f64,
}
impl SleepingBody {
    /// Magnitude of the linear velocity vector.
    pub fn linear_speed(&self) -> f64 {
        let [vx, vy, vz] = self.linear_velocity;
        (vx * vx + vy * vy + vz * vz).sqrt()
    }
    /// Magnitude of the angular velocity vector.
    pub fn angular_speed(&self) -> f64 {
        let [wx, wy, wz] = self.angular_velocity;
        (wx * wx + wy * wy + wz * wz).sqrt()
    }
    /// Approximate kinetic energy: `0.5 * mass * v^2 + 0.5 * inertia * w^2`.
    ///
    /// `mass` is in kg; `inertia` is a scalar moment of inertia (kg*m^2).
    pub fn kinetic_energy(&self, mass: f64, inertia: f64) -> f64 {
        let v = self.linear_speed();
        let w = self.angular_speed();
        0.5 * mass * v * v + 0.5 * inertia * w * w
    }
    /// Kinetic energy using the body's own mass and inertia fields.
    pub fn kinetic_energy_self(&self) -> f64 {
        self.kinetic_energy(self.mass, self.inertia)
    }
    /// Check if this body is below the velocity thresholds.
    pub fn is_below_threshold(&self, linear_threshold: f64, angular_threshold: f64) -> bool {
        self.linear_speed() < linear_threshold && self.angular_speed() < angular_threshold
    }
    /// Check if this body is below the energy threshold.
    pub fn is_below_energy_threshold(&self, energy_threshold: f64) -> bool {
        self.kinetic_energy_self() < energy_threshold
    }
    /// Apply gradual wake-up damping to velocities.
    pub fn apply_wake_damping(&mut self, damping: f64) {
        for v in &mut self.linear_velocity {
            *v *= damping;
        }
        for w in &mut self.angular_velocity {
            *w *= damping;
        }
    }
}
/// A body that is completely frozen in place.
///
/// A frozen body participates in collision detection but is never integrated.
/// It is equivalent to an infinite-mass static object, but retains its
/// original mass for density-based queries.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FrozenBody {
    /// Body ID.
    pub id: u32,
    /// World-space position.
    pub position: [f64; 3],
    /// Orientation as a unit quaternion `[w, x, y, z]`.
    pub orientation: [f64; 4],
    /// Mass (kg) — retained for informational use only.
    pub mass: f64,
    /// Whether the body is currently frozen.
    pub frozen: bool,
}
impl FrozenBody {
    /// Create a new frozen body.
    pub fn new(id: u32, position: [f64; 3], mass: f64) -> Self {
        Self {
            id,
            position,
            orientation: [1.0, 0.0, 0.0, 0.0],
            mass,
            frozen: true,
        }
    }
    /// Freeze the body (disable dynamics).
    pub fn freeze(&mut self) {
        self.frozen = true;
    }
    /// Unfreeze the body (re-enable dynamics).
    pub fn unfreeze(&mut self) {
        self.frozen = false;
    }
    /// Effective inverse mass: infinite mass when frozen (returns 0).
    pub fn inv_mass(&self) -> f64 {
        if self.frozen {
            0.0
        } else {
            1.0 / self.mass.max(f64::EPSILON)
        }
    }
    /// Teleport the body to a new position (allowed even when frozen).
    pub fn set_position(&mut self, pos: [f64; 3]) {
        self.position = pos;
    }
    /// Set orientation as a unit quaternion.
    pub fn set_orientation(&mut self, q: [f64; 4]) {
        self.orientation = q;
    }
}
/// Snapshot statistics about the sleep system.
#[derive(Debug, Clone)]
pub struct SleepStats {
    /// Total number of bodies tracked.
    pub total: usize,
    /// Number of bodies currently sleeping.
    pub sleeping: usize,
    /// Number of bodies currently active.
    pub active: usize,
    /// Number of bodies in the drowsy transition state.
    pub drowsy: usize,
    /// Total kinetic energy of all bodies.
    pub total_energy: f64,
    /// Average kinetic energy per body (0 if no bodies).
    pub avg_energy: f64,
}
impl SleepStats {
    /// Compute statistics from the current state of `manager`.
    pub fn compute(manager: &SleepManager) -> Self {
        let mut stats = Self {
            total: manager.bodies.len(),
            sleeping: 0,
            active: 0,
            drowsy: 0,
            total_energy: 0.0,
            avg_energy: 0.0,
        };
        for body in &manager.bodies {
            match body.state {
                SleepState::Active => stats.active += 1,
                SleepState::Drowsy { .. } => stats.drowsy += 1,
                SleepState::Sleeping => stats.sleeping += 1,
            }
            stats.total_energy += body.kinetic_energy_self();
        }
        if stats.total > 0 {
            stats.avg_energy = stats.total_energy / stats.total as f64;
        }
        stats
    }
    /// Fraction of bodies that are sleeping (0.0 - 1.0).  Returns 0.0 for an
    /// empty manager.
    pub fn sleep_ratio(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.sleeping as f64 / self.total as f64
    }
    /// Fraction of bodies that are active (0.0 - 1.0).
    pub fn active_ratio(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.active as f64 / self.total as f64
    }
}
/// A body in pseudo-sleep: still simulated but at a reduced time step multiple.
///
/// Pseudo-sleeping allows slowly moving bodies to be simulated at a fraction of
/// the full physics rate, saving CPU while avoiding the sharp discontinuity of
/// hard sleeping.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PseudoSleepBody {
    /// Body ID.
    pub id: u32,
    /// Full simulation time step (s).
    pub dt_full: f64,
    /// Fraction of full rate at which pseudo-sleeping bodies are updated.
    /// E.g. 0.25 means one update per 4 full ticks.
    pub rate_fraction: f64,
    /// Internal tick counter.
    pub(super) tick: u32,
    /// Whether pseudo-sleep is active.
    pub active: bool,
}
impl PseudoSleepBody {
    /// Create a new pseudo-sleep body.
    pub fn new(id: u32, dt_full: f64, rate_fraction: f64) -> Self {
        Self {
            id,
            dt_full,
            rate_fraction: rate_fraction.clamp(0.01, 1.0),
            tick: 0,
            active: true,
        }
    }
    /// Advance one tick. Returns `Some(dt_pseudo)` if an update should be
    /// performed this tick, or `None` otherwise.
    pub fn tick(&mut self) -> Option<f64> {
        if !self.active {
            return Some(self.dt_full);
        }
        self.tick += 1;
        let period = (1.0 / self.rate_fraction).round() as u32;
        if self.tick >= period {
            self.tick = 0;
            Some(self.dt_full / self.rate_fraction)
        } else {
            None
        }
    }
    /// Effective time step for one issued update.
    pub fn effective_dt(&self) -> f64 {
        self.dt_full / self.rate_fraction
    }
    /// Disable pseudo-sleep (full simulation rate).
    pub fn disable(&mut self) {
        self.active = false;
        self.tick = 0;
    }
    /// Enable pseudo-sleep.
    pub fn enable(&mut self) {
        self.active = true;
        self.tick = 0;
    }
}
/// Higher-level sleeping system managing energy trackers per body.
pub struct SleepingSystem {
    /// Per-body energy trackers keyed by body ID.
    pub trackers: std::collections::HashMap<u32, EnergyTracker>,
    /// Number of frames to observe before sleeping.
    pub sleep_delay: usize,
    /// Impulse magnitude above which a body should be woken.
    pub wake_impulse_threshold: f64,
}
impl SleepingSystem {
    /// Create a new sleeping system.
    pub fn new(sleep_delay: usize, wake_impulse_threshold: f64) -> Self {
        Self {
            trackers: std::collections::HashMap::new(),
            sleep_delay,
            wake_impulse_threshold,
        }
    }
    /// Update the energy tracker for a body.
    pub fn update_body(&mut self, id: u32, ke: f64, pe: f64) {
        let tracker = self
            .trackers
            .entry(id)
            .or_insert_with(|| EnergyTracker::new(self.sleep_delay.max(1), 0.01));
        tracker.push(ke + pe);
    }
    /// Returns true if the body qualifies for sleep.
    pub fn should_sleep(&self, id: u32) -> bool {
        match self.trackers.get(&id) {
            None => false,
            Some(tracker) => tracker.count >= tracker.max_history && tracker.is_below_threshold(),
        }
    }
    /// Reset the tracker for a body (e.g. on collision or impulse).
    pub fn wake_body(&mut self, id: u32) {
        if let Some(tracker) = self.trackers.get_mut(&id) {
            tracker.count = 0;
            tracker.write_pos = 0;
            for v in &mut tracker.history {
                *v = 0.0;
            }
        }
    }
}
/// Records the sleep/wake event history for a body.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SleepEventType {
    /// Simulation time of the event (s).
    pub time: f64,
    /// Whether this was a sleep event (true) or wake event (false).
    pub is_sleep: bool,
    /// Body ID this event belongs to.
    pub body_id: u32,
}
/// Tunable parameters for the sleep system.
#[derive(Debug, Clone)]
pub struct SleepConfig {
    /// Maximum linear speed (m/s) below which a body is considered still.
    pub linear_threshold: f64,
    /// Maximum angular speed (rad/s) below which a body is considered still.
    pub angular_threshold: f64,
    /// Number of consecutive sub-threshold frames required before a body
    /// transitions from `Drowsy` to `Sleeping`.
    pub drowsy_frames: u32,
    /// Energy threshold for the energy-based sleep criterion (J).
    /// If set to a positive value, energy-based sleeping is used.
    pub energy_threshold: f64,
    /// Gradual wake-up damping factor (0..1). When waking, velocities are
    /// scaled by this factor to avoid sudden jumps.
    pub wake_damping: f64,
    /// Timer-based sleep: seconds of stillness required before sleeping.
    /// If > 0, used instead of frame counting.
    pub sleep_time_seconds: f64,
}
impl SleepConfig {
    /// Sensible defaults: 0.1 m/s, 0.1 rad/s, 60 frames.
    pub fn default() -> Self {
        Self {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            drowsy_frames: 60,
            energy_threshold: 0.0,
            wake_damping: 1.0,
            sleep_time_seconds: 0.0,
        }
    }
    /// Create a config tuned for fast sleeping (useful for stacking scenarios).
    pub fn fast() -> Self {
        Self {
            linear_threshold: 0.05,
            angular_threshold: 0.05,
            drowsy_frames: 10,
            energy_threshold: 0.0,
            wake_damping: 1.0,
            sleep_time_seconds: 0.0,
        }
    }
    /// Create a config that uses energy-based sleep criterion.
    pub fn energy_based(energy_threshold: f64) -> Self {
        Self {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            drowsy_frames: 60,
            energy_threshold,
            wake_damping: 1.0,
            sleep_time_seconds: 0.0,
        }
    }
    /// Create a config that uses timer-based sleeping.
    pub fn timer_based(seconds: f64) -> Self {
        Self {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            drowsy_frames: 60,
            energy_threshold: 0.0,
            wake_damping: 1.0,
            sleep_time_seconds: seconds,
        }
    }
}
/// A recurring sleep/wake pattern for bodies that should be periodically
/// activated (e.g., a motor that runs for a period then rests).
///
/// Bodies matching the pattern are forced awake during "active" windows
/// and allowed to sleep during "rest" windows.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SleepWakePattern {
    /// Duration of the active window (s).
    pub active_duration: f64,
    /// Duration of the rest window (s).
    pub rest_duration: f64,
    /// Current phase timer (s).
    pub phase_timer: f64,
    /// Whether the pattern is currently in the active phase.
    pub is_active: bool,
}
#[allow(dead_code)]
impl SleepWakePattern {
    /// Create a new sleep/wake pattern.
    pub fn new(active_duration: f64, rest_duration: f64) -> Self {
        Self {
            active_duration,
            rest_duration,
            phase_timer: 0.0,
            is_active: true,
        }
    }
    /// Advance the pattern timer by `dt` and update the active phase.
    pub fn step(&mut self, dt: f64) {
        self.phase_timer += dt;
        if self.is_active && self.phase_timer >= self.active_duration {
            self.is_active = false;
            self.phase_timer = 0.0;
        } else if !self.is_active && self.phase_timer >= self.rest_duration {
            self.is_active = true;
            self.phase_timer = 0.0;
        }
    }
    /// Current duty cycle (active fraction).
    pub fn duty_cycle(&self) -> f64 {
        let period = self.active_duration + self.rest_duration;
        if period < 1e-30 {
            return 0.0;
        }
        self.active_duration / period
    }
    /// Apply to a manager: force bodies in `indices` awake if active, else let them sleep.
    pub fn apply_to_manager(&self, manager: &mut SleepManager, indices: &[usize]) {
        if self.is_active {
            for &idx in indices {
                if idx < manager.bodies.len() && manager.bodies[idx].state == SleepState::Sleeping {
                    manager.bodies[idx].state = SleepState::Active;
                }
            }
        }
    }
}
/// Gradually damps body velocity as it approaches the sleep threshold to avoid
/// "snap to sleep" artefacts.
///
/// When a body is in the `Drowsy` state its velocity is scaled by a factor
/// that decreases toward zero over the drowsy window.
#[derive(Debug, Clone)]
pub struct VelocityDamper {
    /// Damping factor applied per frame when body is drowsy ∈ (0, 1].
    /// Values close to 1.0 give slow, smooth damping.
    pub frame_damping: f64,
    /// Whether the damper is active.
    pub enabled: bool,
}
impl VelocityDamper {
    /// Create a velocity damper.
    pub fn new(frame_damping: f64) -> Self {
        Self {
            frame_damping: frame_damping.clamp(0.0, 1.0),
            enabled: true,
        }
    }
    /// Identity damper that does not change velocity.
    pub fn identity() -> Self {
        Self {
            frame_damping: 1.0,
            enabled: false,
        }
    }
    /// Apply damping to linear and angular velocities.
    ///
    /// Returns `(new_linear_vel, new_angular_vel)`.
    pub fn apply(&self, linear: [f64; 3], angular: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        if !self.enabled {
            return (linear, angular);
        }
        let f = self.frame_damping;
        (
            [linear[0] * f, linear[1] * f, linear[2] * f],
            [angular[0] * f, angular[1] * f, angular[2] * f],
        )
    }
    /// Apply in-place to a `SleepingBody` that is in a `Drowsy` state.
    ///
    /// Only acts when the body is `Drowsy`; leaves `Active` and `Sleeping`
    /// bodies unchanged.
    pub fn apply_to_drowsy(&self, body: &mut SleepingBody) {
        if !self.enabled {
            return;
        }
        if matches!(body.state, SleepState::Drowsy { .. }) {
            let f = self.frame_damping;
            body.linear_velocity[0] *= f;
            body.linear_velocity[1] *= f;
            body.linear_velocity[2] *= f;
            body.angular_velocity[0] *= f;
            body.angular_velocity[1] *= f;
            body.angular_velocity[2] *= f;
        }
    }
    /// Apply to all drowsy bodies in a manager.
    pub fn apply_to_manager(&self, manager: &mut SleepManager) {
        if !self.enabled {
            return;
        }
        for body in &mut manager.bodies {
            self.apply_to_drowsy(body);
        }
    }
}
/// Sleep threshold with hysteresis: a lower threshold for entering sleep and
/// a higher threshold for waking.
///
/// This prevents bodies from rapidly toggling between sleeping and active when
/// their speed is near the threshold.
#[derive(Debug, Clone)]
pub struct SleepHysteresis {
    /// Threshold below which a body moves toward sleep.
    pub sleep_threshold: f64,
    /// Threshold above which a sleeping body wakes.
    pub wake_threshold: f64,
}
impl SleepHysteresis {
    /// Create a hysteresis band.  `wake_threshold` must be ≥ `sleep_threshold`.
    pub fn new(sleep_threshold: f64, wake_threshold: f64) -> Self {
        assert!(
            wake_threshold >= sleep_threshold,
            "wake_threshold must be >= sleep_threshold"
        );
        Self {
            sleep_threshold,
            wake_threshold,
        }
    }
    /// Default hysteresis with a 2× wake/sleep ratio.
    pub fn default_band(base_threshold: f64) -> Self {
        Self::new(base_threshold, base_threshold * 2.0)
    }
    /// Determine whether a body should go to sleep given `speed`.
    ///
    /// Only the *linear* speed is checked here; angular velocity should be
    /// tested separately.
    pub fn should_sleep(&self, speed: f64) -> bool {
        speed < self.sleep_threshold
    }
    /// Determine whether a sleeping body should wake given `speed`.
    pub fn should_wake(&self, speed: f64) -> bool {
        speed > self.wake_threshold
    }
    /// Update the sleep state of a body based on its current speed.
    ///
    /// Returns the new state:
    /// - If currently sleeping and speed > wake_threshold → `Active`
    /// - If currently active  and speed < sleep_threshold → `Drowsy { frames: 0 }`
    /// - Otherwise unchanged.
    pub fn update_state(&self, current: &SleepState, speed: f64) -> SleepState {
        match current {
            SleepState::Sleeping => {
                if self.should_wake(speed) {
                    SleepState::Active
                } else {
                    SleepState::Sleeping
                }
            }
            SleepState::Active => {
                if self.should_sleep(speed) {
                    SleepState::Drowsy { frames: 0 }
                } else {
                    SleepState::Active
                }
            }
            SleepState::Drowsy { frames } => {
                if !self.should_sleep(speed) {
                    SleepState::Active
                } else {
                    SleepState::Drowsy { frames: *frames }
                }
            }
        }
    }
    /// Width of the hysteresis band.
    pub fn band_width(&self) -> f64 {
        self.wake_threshold - self.sleep_threshold
    }
}
/// Adapts sleep thresholds based on observed energy statistics.
///
/// Tracks a rolling mean energy of all bodies and adjusts thresholds to
/// reduce false positives (waking sleeping bodies due to numeric noise) while
/// keeping simulation quality high.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SleepThresholdAdapter {
    /// Current linear speed threshold (m/s).
    pub linear_threshold: f64,
    /// Current angular speed threshold (rad/s).
    pub angular_threshold: f64,
    /// Adaptation rate (fraction per step).
    pub adapt_rate: f64,
    /// Minimum allowed threshold (floor).
    pub min_threshold: f64,
    /// Maximum allowed threshold (ceiling).
    pub max_threshold: f64,
    /// Recent mean kinetic energy.
    pub(super) mean_ke: f64,
    /// Exponential moving average coefficient.
    pub(super) ema_alpha: f64,
}
impl SleepThresholdAdapter {
    /// Create a new adapter with given initial thresholds.
    pub fn new(linear: f64, angular: f64, min_t: f64, max_t: f64) -> Self {
        Self {
            linear_threshold: linear,
            angular_threshold: angular,
            adapt_rate: 0.01,
            min_threshold: min_t,
            max_threshold: max_t,
            mean_ke: 0.0,
            ema_alpha: 0.1,
        }
    }
    /// Update the adapter with the current scene mean kinetic energy.
    ///
    /// Thresholds are nudged upward if mean KE is high (active scene) and
    /// downward if the scene is quiet.
    pub fn update(&mut self, scene_mean_ke: f64) {
        self.mean_ke = self.ema_alpha * scene_mean_ke + (1.0 - self.ema_alpha) * self.mean_ke;
        let target = (self.mean_ke * 0.001)
            .sqrt()
            .clamp(self.min_threshold, self.max_threshold);
        let delta = self.adapt_rate * (target - self.linear_threshold);
        self.linear_threshold =
            (self.linear_threshold + delta).clamp(self.min_threshold, self.max_threshold);
        self.angular_threshold =
            (self.angular_threshold + delta).clamp(self.min_threshold, self.max_threshold);
    }
    /// Check if the given speed qualifies a body as "still" under current thresholds.
    pub fn is_still(&self, lin_speed: f64, ang_speed: f64) -> bool {
        lin_speed < self.linear_threshold && ang_speed < self.angular_threshold
    }
}
/// Propagates wake events through a contact graph: when a body collides with
/// a sleeping body, the sleeping body is woken.
///
/// Performs a BFS/DFS through the contact graph starting from all active bodies
/// and wakes any connected sleeping body.
#[derive(Debug, Clone)]
pub struct WakePropagator {
    /// Maximum propagation depth (to avoid cascade waking the whole scene).
    pub max_depth: u32,
}
impl WakePropagator {
    /// Create a wake propagator.
    pub fn new(max_depth: u32) -> Self {
        Self { max_depth }
    }
    /// Propagate wake events through `contacts` (pairs of body indices).
    ///
    /// Active bodies cause adjacent sleeping bodies to wake.  The process
    /// repeats up to `max_depth` times (BFS layers).
    ///
    /// Returns the number of bodies that were woken.
    pub fn propagate(&self, manager: &mut SleepManager, contacts: &[(usize, usize)]) -> usize {
        let mut woken = 0usize;
        for _depth in 0..self.max_depth {
            let mut any_new_wake = false;
            for &(a, b) in contacts {
                if a >= manager.bodies.len() || b >= manager.bodies.len() {
                    continue;
                }
                let a_active = manager.bodies[a].state == SleepState::Active;
                let b_active = manager.bodies[b].state == SleepState::Active;
                let a_sleeping = manager.bodies[a].state == SleepState::Sleeping;
                let b_sleeping = manager.bodies[b].state == SleepState::Sleeping;
                if a_active && b_sleeping {
                    manager.bodies[b].state = SleepState::Active;
                    manager.bodies[b].sleep_timer = 0.0;
                    woken += 1;
                    any_new_wake = true;
                }
                if b_active && a_sleeping {
                    manager.bodies[a].state = SleepState::Active;
                    manager.bodies[a].sleep_timer = 0.0;
                    woken += 1;
                    any_new_wake = true;
                }
            }
            if !any_new_wake {
                break;
            }
        }
        woken
    }
}
/// Adaptive sleep threshold that adjusts based on scene-wide activity.
///
/// The threshold is raised when many bodies are active (fast-moving scene)
/// and lowered when few bodies are active (stable scene), allowing the
/// simulation to put more bodies to sleep in quiet phases.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AdaptiveSleepThreshold {
    /// Minimum threshold (m/s or rad/s).
    pub min_threshold: f64,
    /// Maximum threshold (m/s or rad/s).
    pub max_threshold: f64,
    /// Current threshold.
    pub current_threshold: f64,
    /// Adaptation rate (how fast to adjust per tick).
    pub adaptation_rate: f64,
    /// Target fraction of sleeping bodies.
    pub target_sleep_fraction: f64,
}
#[allow(dead_code)]
impl AdaptiveSleepThreshold {
    /// Create a new adaptive threshold.
    pub fn new(min_t: f64, max_t: f64, adaptation_rate: f64, target_fraction: f64) -> Self {
        Self {
            min_threshold: min_t,
            max_threshold: max_t,
            current_threshold: 0.5 * (min_t + max_t),
            adaptation_rate,
            target_sleep_fraction: target_fraction,
        }
    }
    /// Update the threshold based on the actual sleep fraction.
    ///
    /// If the actual fraction is below the target, raise the threshold
    /// (put more bodies to sleep). If above, lower it.
    pub fn update(&mut self, actual_sleep_fraction: f64) {
        let error = actual_sleep_fraction - self.target_sleep_fraction;
        self.current_threshold -= self.adaptation_rate * error;
        self.current_threshold = self
            .current_threshold
            .clamp(self.min_threshold, self.max_threshold);
    }
    /// Apply the current threshold to a SleepConfig, modifying both linear and angular thresholds.
    pub fn apply_to_config(&self, config: &mut SleepConfig) {
        config.linear_threshold = self.current_threshold;
        config.angular_threshold = self.current_threshold;
    }
}
/// Evaluate the sleep decision for an island of bodies.
///
/// An island sleeps only if *every* body in it is below the sleep thresholds.
/// If any body is active, all sleeping bodies in the island are woken.
#[allow(dead_code)]
pub struct IslandSleepEvaluator {
    /// Energy threshold below which a body is considered still (J).
    pub energy_threshold: f64,
}
impl IslandSleepEvaluator {
    /// Create a new evaluator.
    pub fn new(energy_threshold: f64) -> Self {
        Self { energy_threshold }
    }
    /// Evaluate island sleep decision from per-body kinetic energies.
    pub fn evaluate(&self, energies: &[f64]) -> IslandSleepDecision {
        if energies.is_empty() {
            return IslandSleepDecision::SleepAll;
        }
        let all_below = energies.iter().all(|&e| e < self.energy_threshold);
        let any_below = energies.iter().any(|&e| e < self.energy_threshold);
        if all_below {
            IslandSleepDecision::SleepAll
        } else if !any_below {
            IslandSleepDecision::WakeAll
        } else {
            IslandSleepDecision::Mixed
        }
    }
    /// Number of bodies below threshold.
    pub fn count_sleeping_candidates(&self, energies: &[f64]) -> usize {
        energies
            .iter()
            .filter(|&&e| e < self.energy_threshold)
            .count()
    }
}
/// Tracks how long each body has been sleeping.
#[derive(Debug, Clone)]
pub struct SleepDurationTracker {
    /// Map from body ID → total accumulated sleep time (s).
    pub sleep_duration: std::collections::HashMap<u64, f64>,
    /// Map from body ID → time at which current sleep period started.
    pub(super) sleep_start: std::collections::HashMap<u64, f64>,
}
impl SleepDurationTracker {
    /// Create a new tracker.
    pub fn new() -> Self {
        Self {
            sleep_duration: std::collections::HashMap::new(),
            sleep_start: std::collections::HashMap::new(),
        }
    }
    /// Notify the tracker that body `id` fell asleep at `time`.
    pub fn record_sleep_start(&mut self, id: u64, time: f64) {
        self.sleep_start.insert(id, time);
    }
    /// Notify the tracker that body `id` woke up at `time`.
    ///
    /// Accumulates the elapsed sleep time.
    pub fn record_wake(&mut self, id: u64, time: f64) {
        if let Some(&start) = self.sleep_start.get(&id) {
            let elapsed = (time - start).max(0.0);
            *self.sleep_duration.entry(id).or_insert(0.0) += elapsed;
            self.sleep_start.remove(&id);
        }
    }
    /// Total accumulated sleep time for body `id`.
    pub fn total_sleep_time(&self, id: u64) -> f64 {
        *self.sleep_duration.get(&id).unwrap_or(&0.0)
    }
    /// Whether body `id` is currently sleeping (started but not yet woken).
    pub fn is_currently_sleeping(&self, id: u64) -> bool {
        self.sleep_start.contains_key(&id)
    }
    /// Current sleep duration for a body that is still sleeping.
    pub fn current_sleep_duration(&self, id: u64, now: f64) -> f64 {
        if let Some(&start) = self.sleep_start.get(&id) {
            (now - start).max(0.0)
        } else {
            0.0
        }
    }
    /// Reset all tracking data.
    pub fn reset(&mut self) {
        self.sleep_duration.clear();
        self.sleep_start.clear();
    }
}
/// Lifecycle state of a rigid body with respect to the sleep system.
#[derive(Debug, Clone, PartialEq)]
pub enum SleepState {
    /// Body is fully simulated every frame.
    Active,
    /// Body has been below the velocity threshold for `frames` consecutive
    /// frames and is a candidate for sleeping.
    Drowsy {
        /// Number of consecutive sub-threshold frames so far.
        frames: u32,
    },
    /// Body is fully at rest; the integrator may skip it.
    Sleeping,
}
/// Sleep management at the island (connected-component) level.
///
/// Bodies in the same contact island must all be still before any of them
/// can sleep; otherwise a sleeping body can be "stuck" while a touching body
/// pushes it.
pub struct IslandSleepManager {
    /// Sleep manager that owns the bodies.
    pub manager: SleepManager,
}
impl IslandSleepManager {
    /// Create a new island sleep manager.
    pub fn new(config: SleepConfig) -> Self {
        Self {
            manager: SleepManager::new(config),
        }
    }
    /// Partition `bodies` into connected components using Union-Find on
    /// `contact_pairs`.  Returns a list of islands, each island being a list
    /// of indices into `bodies`.
    pub fn group_into_islands(
        bodies: &[SleepingBody],
        contact_pairs: &[(usize, usize)],
    ) -> Vec<Vec<usize>> {
        let n = bodies.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for &(a, b) in contact_pairs {
            let ra = find(&mut parent, a);
            let rb = find(&mut parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
        let mut map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for i in 0..n {
            let root = find(&mut parent, i);
            map.entry(root).or_default().push(i);
        }
        map.into_values().collect()
    }
    /// Return `true` if every body in `island` is below both velocity
    /// thresholds in `config`.
    pub fn island_can_sleep(
        island: &[usize],
        bodies: &[SleepingBody],
        config: &SleepConfig,
    ) -> bool {
        island.iter().all(|&i| {
            if config.energy_threshold > 0.0 {
                bodies[i].is_below_energy_threshold(config.energy_threshold)
            } else {
                bodies[i].linear_speed() < config.linear_threshold
                    && bodies[i].angular_speed() < config.angular_threshold
            }
        })
    }
    /// Immediately set every body in `island` to `Sleeping`.
    pub fn sleep_island(&mut self, island: &[usize]) {
        for &i in island {
            if i < self.manager.bodies.len() {
                self.manager.bodies[i].state = SleepState::Sleeping;
                self.manager.bodies[i].sleep_timer = 0.0;
            }
        }
    }
    /// Wake all bodies in `island`.
    pub fn wake_island(&mut self, island: &[usize]) {
        for &i in island {
            if i < self.manager.bodies.len() {
                self.manager.bodies[i].state = SleepState::Active;
                self.manager.bodies[i].sleep_timer = 0.0;
            }
        }
    }
    /// Process all islands: sleep islands that can sleep, wake those that cannot.
    pub fn process_islands(&mut self, contact_pairs: &[(usize, usize)]) {
        let islands = Self::group_into_islands(&self.manager.bodies, contact_pairs);
        let config = self.manager.config.clone();
        for island in &islands {
            if Self::island_can_sleep(island, &self.manager.bodies, &config) {
                for &i in island {
                    if i < self.manager.bodies.len() {
                        self.manager.bodies[i].state = SleepState::Sleeping;
                    }
                }
            } else {
                for &i in island {
                    if i < self.manager.bodies.len()
                        && self.manager.bodies[i].state == SleepState::Sleeping
                    {
                        self.manager.bodies[i].state = SleepState::Active;
                        self.manager.bodies[i].sleep_timer = 0.0;
                    }
                }
            }
        }
    }
}
impl IslandSleepManager {
    /// Propagate wake signals through the island graph.
    ///
    /// Any body that was woken externally will propagate its wake signal
    /// to all bodies in the same island transitively. This ensures that
    /// if one body in a resting stack is hit, all connected bodies wake.
    pub fn propagate_wake(&mut self, contact_pairs: &[(usize, usize)]) {
        let woke: Vec<usize> = self
            .manager
            .bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| b.state != SleepState::Sleeping)
            .map(|(i, _)| i)
            .collect();
        if woke.is_empty() {
            return;
        }
        let islands = Self::group_into_islands(&self.manager.bodies, contact_pairs);
        for island in &islands {
            let has_woken = island.iter().any(|&i| woke.contains(&i));
            if has_woken {
                for &i in island {
                    if i < self.manager.bodies.len()
                        && self.manager.bodies[i].state == SleepState::Sleeping
                    {
                        self.manager.bodies[i].state = SleepState::Active;
                        self.manager.bodies[i].sleep_timer = 0.0;
                    }
                }
            }
        }
    }
    /// Island-aware tick: run tick on manager, then enforce island consistency.
    ///
    /// After the standard tick, if any body in an island is Active,
    /// wake all Sleeping members of that island.
    pub fn tick_islands(&mut self, contact_pairs: &[(usize, usize)]) {
        self.manager.tick();
        self.propagate_wake(contact_pairs);
    }
    /// Count the number of distinct islands.
    pub fn island_count(&self, contact_pairs: &[(usize, usize)]) -> usize {
        Self::group_into_islands(&self.manager.bodies, contact_pairs).len()
    }
    /// Return indices of all sleeping bodies.
    pub fn sleeping_indices(&self) -> Vec<usize> {
        self.manager
            .bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| b.state == SleepState::Sleeping)
            .map(|(i, _)| i)
            .collect()
    }
    /// Return indices of all active bodies.
    pub fn active_indices(&self) -> Vec<usize> {
        self.manager
            .bodies
            .iter()
            .enumerate()
            .filter(|(_, b)| b.state != SleepState::Sleeping)
            .map(|(i, _)| i)
            .collect()
    }
    /// Force all bodies in the manager to sleep immediately.
    pub fn force_all_to_sleep(&mut self) {
        for b in &mut self.manager.bodies {
            b.state = SleepState::Sleeping;
        }
    }
}
/// Result of an island-level sleep evaluation.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum IslandSleepDecision {
    /// All bodies in the island should sleep.
    SleepAll,
    /// All bodies in the island should be woken.
    WakeAll,
    /// Mixed — individual bodies have different states.
    Mixed,
}
/// Tracks kinetic energy over time to estimate when a body will naturally
/// decay below the sleep threshold.
#[allow(dead_code)]
pub struct EnergyDecayTracker {
    /// Body ID.
    pub id: u32,
    /// Most recent energy (J).
    pub current_energy: f64,
    /// Previous energy (J).
    pub prev_energy: f64,
    /// Estimated time constant τ (s) for exponential decay.
    pub tau: f64,
    /// Time since last update (s).
    pub elapsed: f64,
}
impl EnergyDecayTracker {
    /// Create a new decay tracker.
    pub fn new(id: u32) -> Self {
        Self {
            id,
            current_energy: 0.0,
            prev_energy: 0.0,
            tau: f64::INFINITY,
            elapsed: 0.0,
        }
    }
    /// Update with a new energy measurement.
    pub fn update(&mut self, energy: f64, dt: f64) {
        self.prev_energy = self.current_energy;
        self.current_energy = energy;
        self.elapsed += dt;
        if self.prev_energy > 1e-30 && energy > 1e-30 && dt > 1e-15 {
            let ratio = energy / self.prev_energy;
            if ratio > 0.0 && ratio < 1.0 {
                self.tau = -dt / ratio.ln();
            }
        }
    }
    /// Estimated time until energy falls below `threshold` (s).
    /// Returns `None` if energy is already below threshold or tau is infinite.
    pub fn time_to_threshold(&self, threshold: f64) -> Option<f64> {
        if self.current_energy <= threshold {
            return None;
        }
        if self.tau.is_infinite() || self.tau <= 0.0 {
            return None;
        }
        let t = self.tau * (self.current_energy / threshold).ln();
        Some(t.max(0.0))
    }
    /// Returns true if energy is decaying (current < previous).
    pub fn is_decaying(&self) -> bool {
        self.current_energy < self.prev_energy
    }
}
/// Ring-buffer energy tracker for determining if a body should sleep.
pub struct EnergyTracker {
    /// Energy history (ring buffer).
    pub history: Vec<f64>,
    /// Maximum number of history entries.
    pub max_history: usize,
    /// Energy threshold; all values must be below this to qualify for sleep.
    pub threshold: f64,
    /// Current write position in the ring buffer.
    pub(super) write_pos: usize,
    /// Number of valid entries.
    pub(super) count: usize,
}
impl EnergyTracker {
    /// Create a new energy tracker.
    pub fn new(max_history: usize, threshold: f64) -> Self {
        Self {
            history: vec![0.0; max_history.max(1)],
            max_history: max_history.max(1),
            threshold,
            write_pos: 0,
            count: 0,
        }
    }
    /// Push a new energy value into the ring buffer.
    pub fn push(&mut self, energy: f64) {
        self.history[self.write_pos] = energy;
        self.write_pos = (self.write_pos + 1) % self.max_history;
        if self.count < self.max_history {
            self.count += 1;
        }
    }
    /// Returns true if all values in the buffer are below the threshold.
    pub fn is_below_threshold(&self) -> bool {
        if self.count == 0 {
            return true;
        }
        self.history[..self.count]
            .iter()
            .all(|&e| e < self.threshold)
    }
    /// Mean energy of the values in the buffer. Returns 0.0 if empty.
    pub fn mean_energy(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        let sum: f64 = self.history[..self.count].iter().sum();
        sum / self.count as f64
    }
}
/// Predicts when a body will transition to sleep, based on exponential
/// velocity decay under damping.
///
/// Models: `v(t) = v0 * exp(-damping * t)`
/// Sleep occurs when `v(t) < threshold`, i.e. `t > -ln(threshold/v0) / damping`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SleepPredictor {
    /// Linear damping coefficient (1/s).
    pub linear_damping: f64,
    /// Angular damping coefficient (1/s).
    pub angular_damping: f64,
    /// Velocity threshold (m/s) below which sleep is triggered.
    pub linear_threshold: f64,
    /// Angular velocity threshold (rad/s).
    pub angular_threshold: f64,
}
#[allow(dead_code)]
impl SleepPredictor {
    /// Create a new sleep predictor.
    pub fn new(
        linear_damping: f64,
        angular_damping: f64,
        linear_threshold: f64,
        angular_threshold: f64,
    ) -> Self {
        Self {
            linear_damping,
            angular_damping,
            linear_threshold,
            angular_threshold,
        }
    }
    /// Predict time until linear speed falls below threshold.
    ///
    /// Returns `f64::INFINITY` if the body won't reach the threshold (no damping),
    /// and `0.0` if already below threshold.
    pub fn time_to_linear_sleep(&self, current_speed: f64) -> f64 {
        if current_speed <= self.linear_threshold {
            return 0.0;
        }
        if self.linear_damping < 1e-30 {
            return f64::INFINITY;
        }
        (current_speed / self.linear_threshold).ln() / self.linear_damping
    }
    /// Predict time until angular speed falls below threshold.
    pub fn time_to_angular_sleep(&self, current_angular_speed: f64) -> f64 {
        if current_angular_speed <= self.angular_threshold {
            return 0.0;
        }
        if self.angular_damping < 1e-30 {
            return f64::INFINITY;
        }
        (current_angular_speed / self.angular_threshold).ln() / self.angular_damping
    }
    /// Predict the total time until a body sleeps (max of linear and angular).
    pub fn time_to_sleep(&self, body: &SleepingBody) -> f64 {
        let t_lin = self.time_to_linear_sleep(body.linear_speed());
        let t_ang = self.time_to_angular_sleep(body.angular_speed());
        t_lin.max(t_ang)
    }
    /// Predict the velocity of a body at time `t` from now.
    pub fn predicted_linear_speed(&self, current_speed: f64, t: f64) -> f64 {
        current_speed * (-self.linear_damping * t).exp()
    }
    /// Predict the angular speed of a body at time `t` from now.
    pub fn predicted_angular_speed(&self, current_angular_speed: f64, t: f64) -> f64 {
        current_angular_speed * (-self.angular_damping * t).exp()
    }
}
/// Classification of a body's current motion level.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ActivityLevel {
    /// Body is fully at rest (below sleep threshold).
    Resting,
    /// Body is barely moving (drowsy zone).
    Minimal,
    /// Body has moderate velocity.
    Moderate,
    /// Body is moving fast.
    High,
}
impl ActivityLevel {
    /// Classify from linear speed (m/s) using standard thresholds.
    pub fn from_speed(speed: f64) -> Self {
        if speed < 0.05 {
            ActivityLevel::Resting
        } else if speed < 0.5 {
            ActivityLevel::Minimal
        } else if speed < 3.0 {
            ActivityLevel::Moderate
        } else {
            ActivityLevel::High
        }
    }
    /// Classify from kinetic energy using configurable thresholds.
    pub fn from_kinetic_energy(
        ke: f64,
        rest_thresh: f64,
        mod_thresh: f64,
        high_thresh: f64,
    ) -> Self {
        if ke < rest_thresh {
            ActivityLevel::Resting
        } else if ke < mod_thresh {
            ActivityLevel::Minimal
        } else if ke < high_thresh {
            ActivityLevel::Moderate
        } else {
            ActivityLevel::High
        }
    }
    /// Whether this level requires full simulation.
    pub fn needs_simulation(&self) -> bool {
        !matches!(self, ActivityLevel::Resting)
    }
}
/// Records sleep/wake transitions for diagnostics.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct SleepHistory {
    /// Recorded events in order.
    pub events: Vec<SleepEventType>,
}
#[allow(dead_code)]
impl SleepHistory {
    /// Create an empty history.
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }
    /// Record a sleep event for a body.
    pub fn record_sleep(&mut self, time: f64, body_id: u32) {
        self.events.push(SleepEventType {
            time,
            is_sleep: true,
            body_id,
        });
    }
    /// Record a wake event for a body.
    pub fn record_wake(&mut self, time: f64, body_id: u32) {
        self.events.push(SleepEventType {
            time,
            is_sleep: false,
            body_id,
        });
    }
    /// Total number of sleep events recorded.
    pub fn sleep_count(&self) -> usize {
        self.events.iter().filter(|e| e.is_sleep).count()
    }
    /// Total number of wake events recorded.
    pub fn wake_count(&self) -> usize {
        self.events.iter().filter(|e| !e.is_sleep).count()
    }
    /// Events for a specific body.
    pub fn events_for(&self, body_id: u32) -> Vec<&SleepEventType> {
        self.events
            .iter()
            .filter(|e| e.body_id == body_id)
            .collect()
    }
    /// Clear all recorded events.
    pub fn clear(&mut self) {
        self.events.clear();
    }
}
/// Manages sleep/wake transitions for a collection of rigid bodies.
pub struct SleepManager {
    /// All tracked bodies.
    pub bodies: Vec<SleepingBody>,
    /// Configuration shared by all bodies.
    pub config: SleepConfig,
}
impl SleepManager {
    /// Create a new manager with the given configuration.
    pub fn new(config: SleepConfig) -> Self {
        Self {
            bodies: Vec::new(),
            config,
        }
    }
    /// Register a new body at `pos` and return its index.
    pub fn add_body(&mut self, id: u32, pos: [f64; 3]) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(SleepingBody {
            body_id: id,
            state: SleepState::Active,
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            position: pos,
            orientation: [0.0, 0.0, 0.0, 1.0],
            sleep_timer: 0.0,
            mass: 1.0,
            inertia: 1.0,
        });
        idx
    }
    /// Register a new body with mass and inertia for energy-based sleep.
    pub fn add_body_with_mass(&mut self, id: u32, pos: [f64; 3], mass: f64, inertia: f64) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(SleepingBody {
            body_id: id,
            state: SleepState::Active,
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            position: pos,
            orientation: [0.0, 0.0, 0.0, 1.0],
            sleep_timer: 0.0,
            mass,
            inertia,
        });
        idx
    }
    /// Overwrite the velocity of the body at `idx`.
    pub fn update_body_velocity(&mut self, idx: usize, lin_vel: [f64; 3], ang_vel: [f64; 3]) {
        let b = &mut self.bodies[idx];
        b.linear_velocity = lin_vel;
        b.angular_velocity = ang_vel;
    }
    /// Overwrite the position and orientation of the body at `idx`.
    pub fn update_body_position(&mut self, idx: usize, pos: [f64; 3], orient: [f64; 4]) {
        let b = &mut self.bodies[idx];
        b.position = pos;
        b.orientation = orient;
    }
    /// Check if a body should sleep based on the configured criterion.
    #[allow(dead_code)]
    fn should_sleep(&self, body: &SleepingBody) -> bool {
        if self.config.energy_threshold > 0.0 {
            body.is_below_energy_threshold(self.config.energy_threshold)
        } else {
            body.is_below_threshold(self.config.linear_threshold, self.config.angular_threshold)
        }
    }
    /// Advance the sleep state machine by one simulation frame.
    ///
    /// Transition rules:
    /// * `Active`  + below threshold -> `Drowsy { frames: 1 }`
    /// * `Drowsy`  + below threshold + frames+1 < drowsy_frames -> `Drowsy { frames+1 }`
    /// * `Drowsy`  + below threshold + frames+1 >= drowsy_frames -> `Sleeping`
    /// * `Drowsy`  + above threshold -> `Active`
    /// * `Sleeping` -> stays `Sleeping` (use [`wake`](Self::wake) to revive)
    pub fn tick(&mut self) {
        let needed = self.config.drowsy_frames;
        for body in &mut self.bodies {
            let below = if self.config.energy_threshold > 0.0 {
                body.is_below_energy_threshold(self.config.energy_threshold)
            } else {
                body.linear_speed() < self.config.linear_threshold
                    && body.angular_speed() < self.config.angular_threshold
            };
            body.state = match body.state {
                SleepState::Active => {
                    if below {
                        SleepState::Drowsy { frames: 1 }
                    } else {
                        SleepState::Active
                    }
                }
                SleepState::Drowsy { frames } => {
                    if below {
                        let next = frames + 1;
                        if next >= needed {
                            SleepState::Sleeping
                        } else {
                            SleepState::Drowsy { frames: next }
                        }
                    } else {
                        SleepState::Active
                    }
                }
                SleepState::Sleeping => SleepState::Sleeping,
            };
        }
    }
    /// Advance the sleep state machine using a timer-based criterion.
    ///
    /// `dt` is the simulation time step in seconds.
    pub fn tick_timer(&mut self, dt: f64) {
        let sleep_time = self.config.sleep_time_seconds;
        if sleep_time <= 0.0 {
            self.tick();
            return;
        }
        for body in &mut self.bodies {
            let below = if self.config.energy_threshold > 0.0 {
                body.is_below_energy_threshold(self.config.energy_threshold)
            } else {
                body.linear_speed() < self.config.linear_threshold
                    && body.angular_speed() < self.config.angular_threshold
            };
            match body.state {
                SleepState::Sleeping => {}
                _ => {
                    if below {
                        body.sleep_timer += dt;
                        if body.sleep_timer >= sleep_time {
                            body.state = SleepState::Sleeping;
                        } else {
                            body.state = SleepState::Drowsy { frames: 0 };
                        }
                    } else {
                        body.sleep_timer = 0.0;
                        body.state = SleepState::Active;
                    }
                }
            }
        }
    }
    /// Force the body at `idx` to the `Active` state regardless of velocity.
    pub fn wake(&mut self, idx: usize) {
        self.bodies[idx].state = SleepState::Active;
        self.bodies[idx].sleep_timer = 0.0;
    }
    /// Wake with gradual damping: sets body active and scales velocities.
    pub fn wake_gradual(&mut self, idx: usize) {
        let damping = self.config.wake_damping;
        self.bodies[idx].state = SleepState::Active;
        self.bodies[idx].sleep_timer = 0.0;
        self.bodies[idx].apply_wake_damping(damping);
    }
    /// Wake every sleeping body whose position is within `radius` of `pos`.
    pub fn wake_nearby(&mut self, pos: [f64; 3], radius: f64) {
        let r2 = radius * radius;
        for body in &mut self.bodies {
            if body.state == SleepState::Sleeping {
                let [bx, by, bz] = body.position;
                let [px, py, pz] = pos;
                let dx = bx - px;
                let dy = by - py;
                let dz = bz - pz;
                if dx * dx + dy * dy + dz * dz <= r2 {
                    body.state = SleepState::Active;
                    body.sleep_timer = 0.0;
                }
            }
        }
    }
    /// Wake nearby bodies with gradual damping.
    pub fn wake_nearby_gradual(&mut self, pos: [f64; 3], radius: f64) {
        let r2 = radius * radius;
        let damping = self.config.wake_damping;
        for body in &mut self.bodies {
            if body.state == SleepState::Sleeping {
                let [bx, by, bz] = body.position;
                let [px, py, pz] = pos;
                let dx = bx - px;
                let dy = by - py;
                let dz = bz - pz;
                if dx * dx + dy * dy + dz * dz <= r2 {
                    body.state = SleepState::Active;
                    body.sleep_timer = 0.0;
                    body.apply_wake_damping(damping);
                }
            }
        }
    }
    /// Return the number of sleeping bodies.
    pub fn sleeping_count(&self) -> usize {
        self.bodies
            .iter()
            .filter(|b| b.state == SleepState::Sleeping)
            .count()
    }
    /// Return the number of active bodies.
    pub fn active_count(&self) -> usize {
        self.bodies
            .iter()
            .filter(|b| b.state == SleepState::Active)
            .count()
    }
    /// Return the total kinetic energy of all bodies.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(|b| b.kinetic_energy_self()).sum()
    }
    /// Update linear and angular thresholds at runtime.
    pub fn set_thresholds(&mut self, linear: f64, angular: f64) {
        self.config.linear_threshold = linear;
        self.config.angular_threshold = angular;
    }
}
