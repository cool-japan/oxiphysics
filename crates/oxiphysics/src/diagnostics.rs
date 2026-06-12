// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Diagnostics and performance monitoring for OxiPhysics simulations.
//!
//! This module provides tools for measuring timing, tracking energy
//! conservation, monitoring constraint solver convergence, and
//! collecting memory usage estimates across all physics subsystems.
//!
//! # Quick-start
//!
//! ```no_run
//! use oxiphysics::diagnostics::{PhysicsTimer, SimulationDiagnostics};
//!
//! let mut timer = PhysicsTimer::new();
//! timer.start();
//! // ... run simulation step ...
//! timer.stop();
//!
//! let mut diag = SimulationDiagnostics::default();
//! diag.step_count += 1;
//! diag.n_rigid_bodies = 10;
//! ```

use std::collections::VecDeque;
use std::time::{Duration, Instant};

// ---------------------------------------------------------------------------
// PhysicsTimer
// ---------------------------------------------------------------------------

/// High-resolution wall-clock timer for profiling physics stages.
///
/// Supports lap timing, split recording, and elapsed-millisecond queries.
pub struct PhysicsTimer {
    /// Instant at which the last `start()` was called.
    start_instant: Option<Instant>,
    /// Accumulated duration from completed `start`/`stop` pairs.
    accumulated: Duration,
    /// Lap durations recorded via [`PhysicsTimer::lap`].
    lap_times: Vec<Duration>,
    /// Split durations recorded via [`PhysicsTimer::split_times`].
    splits: Vec<Duration>,
    /// Whether the timer is currently running.
    running: bool,
}

impl Default for PhysicsTimer {
    fn default() -> Self {
        Self::new()
    }
}

impl PhysicsTimer {
    /// Create a new stopped timer with no accumulated time.
    pub fn new() -> Self {
        Self {
            start_instant: None,
            accumulated: Duration::ZERO,
            lap_times: Vec::new(),
            splits: Vec::new(),
            running: false,
        }
    }

    /// Start the timer.  Has no effect if already running.
    pub fn start(&mut self) {
        if !self.running {
            self.start_instant = Some(Instant::now());
            self.running = true;
        }
    }

    /// Stop the timer and accumulate the elapsed time.  Has no effect if not running.
    pub fn stop(&mut self) {
        if self.running {
            if let Some(t0) = self.start_instant.take() {
                self.accumulated += t0.elapsed();
            }
            self.running = false;
        }
    }

    /// Return the total accumulated time (including any currently-running interval).
    pub fn elapsed(&self) -> Duration {
        if self.running && let Some(t0) = self.start_instant {
            return self.accumulated + t0.elapsed();
        }
        self.accumulated
    }

    /// Return the total accumulated time in milliseconds.
    pub fn elapsed_ms(&self) -> f64 {
        self.elapsed().as_secs_f64() * 1000.0
    }

    /// Record a lap: stores the current elapsed time as a lap split and
    /// continues running.
    pub fn lap(&mut self) {
        let t = self.elapsed();
        self.lap_times.push(t);
    }

    /// Return a snapshot of all lap times recorded so far.
    pub fn split_times(&self) -> &[Duration] {
        &self.lap_times
    }

    /// Record a split time (alias for [`PhysicsTimer::lap`], stored separately).
    pub fn split(&mut self) {
        let t = self.elapsed();
        self.splits.push(t);
    }

    /// Return a snapshot of all split times recorded so far.
    pub fn splits(&self) -> &[Duration] {
        &self.splits
    }

    /// Reset the timer to its initial state.
    pub fn reset(&mut self) {
        *self = Self::new();
    }

    /// Return whether the timer is currently running.
    pub fn is_running(&self) -> bool {
        self.running
    }
}

// ---------------------------------------------------------------------------
// SimulationDiagnostics
// ---------------------------------------------------------------------------

/// Snapshot of simulation-wide statistics at a single simulation step.
#[derive(Debug, Clone, Default)]
pub struct SimulationDiagnostics {
    /// Cumulative number of simulation steps executed.
    pub step_count: u64,
    /// Elapsed simulated time (s).
    pub sim_time: f64,
    /// Elapsed wall-clock time since the simulation was created (s).
    pub wall_time: f64,
    /// Number of active rigid bodies.
    pub n_rigid_bodies: usize,
    /// Number of active constraint edges (joints + contacts combined).
    pub n_constraints: usize,
    /// Number of active contact points.
    pub n_contacts: usize,
    /// Number of sleeping rigid bodies.
    pub n_sleeping: usize,
}

impl SimulationDiagnostics {
    /// Create a zeroed diagnostics snapshot.
    pub fn new() -> Self {
        Self::default()
    }

    /// Fraction of bodies currently sleeping (0.0 if no bodies).
    pub fn sleeping_fraction(&self) -> f64 {
        if self.n_rigid_bodies == 0 {
            0.0
        } else {
            self.n_sleeping as f64 / self.n_rigid_bodies as f64
        }
    }

    /// Simulated-time-to-wall-time ratio (real-time factor).
    ///
    /// A value > 1.0 means faster-than-real-time; < 1.0 means slower.
    pub fn real_time_factor(&self) -> f64 {
        if self.wall_time < 1e-15 {
            0.0
        } else {
            self.sim_time / self.wall_time
        }
    }
}

// ---------------------------------------------------------------------------
// EnergyMonitor
// ---------------------------------------------------------------------------

/// Tracks kinetic, potential, and total mechanical energy per step.
///
/// Also records the momentum conservation error to detect instabilities.
#[derive(Debug, Clone, Default)]
pub struct EnergyMonitor {
    /// Kinetic energy (J) at each recorded step.
    pub kinetic: Vec<f64>,
    /// Potential energy (J) at each recorded step.
    pub potential: Vec<f64>,
    /// Total mechanical energy (J) at each recorded step.
    pub total: Vec<f64>,
    /// Linear momentum magnitude conservation error (kg m s⁻¹) per step.
    pub momentum_error: Vec<f64>,
}

impl EnergyMonitor {
    /// Create a new empty energy monitor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record energy and momentum conservation data for one step.
    pub fn record(&mut self, kinetic: f64, potential: f64, momentum_err: f64) {
        let total = kinetic + potential;
        self.kinetic.push(kinetic);
        self.potential.push(potential);
        self.total.push(total);
        self.momentum_error.push(momentum_err);
    }

    /// Return the maximum energy drift relative to the first recorded value.
    ///
    /// Returns `0.0` if fewer than two steps have been recorded.
    pub fn max_energy_drift(&self) -> f64 {
        if self.total.len() < 2 {
            return 0.0;
        }
        let e0 = self.total[0].abs().max(1e-30);
        self.total
            .iter()
            .map(|e| (e - self.total[0]).abs() / e0)
            .fold(0.0_f64, f64::max)
    }

    /// Return the average kinetic energy across all recorded steps.
    pub fn avg_kinetic(&self) -> f64 {
        if self.kinetic.is_empty() {
            return 0.0;
        }
        self.kinetic.iter().sum::<f64>() / self.kinetic.len() as f64
    }

    /// Return the average potential energy across all recorded steps.
    pub fn avg_potential(&self) -> f64 {
        if self.potential.is_empty() {
            return 0.0;
        }
        self.potential.iter().sum::<f64>() / self.potential.len() as f64
    }

    /// Clear all recorded data.
    pub fn reset(&mut self) {
        *self = Self::new();
    }

    /// Return the number of recorded steps.
    pub fn len(&self) -> usize {
        self.total.len()
    }

    /// Return `true` if no steps have been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.total.is_empty()
    }
}

// ---------------------------------------------------------------------------
// ConstraintDiagnostics
// ---------------------------------------------------------------------------

/// Constraint solver convergence diagnostics.
#[derive(Debug, Clone, Default)]
pub struct ConstraintDiagnostics {
    /// Number of solver iterations actually taken for each step.
    pub solver_iterations: Vec<usize>,
    /// Final solver residual norm at each step.
    pub residuals: Vec<f64>,
    /// Warm-start quality metric in \[0, 1\] (higher = better initial guess).
    pub warm_start_quality: f64,
}

impl ConstraintDiagnostics {
    /// Create a new empty constraint diagnostics object.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record the solver result for one step.
    pub fn record(&mut self, iters: usize, residual: f64) {
        self.solver_iterations.push(iters);
        self.residuals.push(residual);
    }

    /// Return the average number of solver iterations.
    pub fn avg_iterations(&self) -> f64 {
        if self.solver_iterations.is_empty() {
            return 0.0;
        }
        self.solver_iterations.iter().sum::<usize>() as f64
            / self.solver_iterations.len() as f64
    }

    /// Return the maximum residual across all recorded steps.
    pub fn max_residual(&self) -> f64 {
        self.residuals
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }

    /// Return the number of recorded steps.
    pub fn len(&self) -> usize {
        self.solver_iterations.len()
    }

    /// Return `true` if no steps have been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.solver_iterations.is_empty()
    }
}

// ---------------------------------------------------------------------------
// CollisionDiagnostics
// ---------------------------------------------------------------------------

/// Collision detection performance diagnostics.
#[derive(Debug, Clone, Default)]
pub struct CollisionDiagnostics {
    /// Number of candidate pairs produced by the broadphase in the last step.
    pub broadphase_pairs: usize,
    /// Number of narrowphase overlap tests performed in the last step.
    pub narrowphase_checks: usize,
    /// Number of contact manifolds generated in the last step.
    pub contacts_generated: usize,
    /// Number of CCD sub-steps taken in the last step.
    pub ccd_steps: usize,
}

impl CollisionDiagnostics {
    /// Create a new zeroed collision diagnostics object.
    pub fn new() -> Self {
        Self::default()
    }

    /// Ratio of narrowphase checks to broadphase pairs (false-positive rate + 1).
    ///
    /// Returns `0.0` if no broadphase pairs were reported.
    pub fn narrowphase_ratio(&self) -> f64 {
        if self.broadphase_pairs == 0 {
            0.0
        } else {
            self.narrowphase_checks as f64 / self.broadphase_pairs as f64
        }
    }

    /// Reset all counters to zero.
    pub fn reset(&mut self) {
        *self = Self::new();
    }
}

// ---------------------------------------------------------------------------
// MemoryStats
// ---------------------------------------------------------------------------

/// Rough memory usage estimates (bytes) per physics subsystem.
///
/// These are programmer estimates rather than OS-reported figures.
#[derive(Debug, Clone, Default)]
pub struct MemoryStats {
    /// Estimated bytes used by rigid body data.
    pub rigid_bytes: usize,
    /// Estimated bytes used by collision data structures.
    pub collision_bytes: usize,
    /// Estimated bytes used by constraint data.
    pub constraint_bytes: usize,
    /// Estimated bytes used by LBM grid data.
    pub lbm_bytes: usize,
    /// Estimated bytes used by SPH particle data.
    pub sph_bytes: usize,
    /// Estimated bytes used by FEM mesh and stiffness data.
    pub fem_bytes: usize,
}

impl MemoryStats {
    /// Create a zeroed memory stats object.
    pub fn new() -> Self {
        Self::default()
    }

    /// Total estimated memory across all subsystems (bytes).
    pub fn total_bytes(&self) -> usize {
        self.rigid_bytes
            + self.collision_bytes
            + self.constraint_bytes
            + self.lbm_bytes
            + self.sph_bytes
            + self.fem_bytes
    }

    /// Total estimated memory in megabytes.
    pub fn total_mb(&self) -> f64 {
        self.total_bytes() as f64 / (1024.0 * 1024.0)
    }
}

// ---------------------------------------------------------------------------
// PerformanceSummary
// ---------------------------------------------------------------------------

/// Aggregates all diagnostic subsystems into a single summary.
#[derive(Debug, Clone, Default)]
pub struct PerformanceSummary {
    /// High-level simulation statistics.
    pub simulation: SimulationDiagnostics,
    /// Energy tracking data.
    pub energy: EnergyMonitor,
    /// Constraint solver diagnostics.
    pub constraints: ConstraintDiagnostics,
    /// Collision detection diagnostics.
    pub collision: CollisionDiagnostics,
    /// Memory usage estimates.
    pub memory: MemoryStats,
}

impl PerformanceSummary {
    /// Create a new zeroed performance summary.
    pub fn new() -> Self {
        Self::default()
    }

    /// Print a human-readable summary to stdout.
    pub fn print_summary(&self) {
        println!("=== OxiPhysics Performance Summary ===");
        println!("Steps          : {}", self.simulation.step_count);
        println!("Sim time       : {:.4} s", self.simulation.sim_time);
        println!("Wall time      : {:.4} s", self.simulation.wall_time);
        println!("Real-time factor: {:.2}x", self.simulation.real_time_factor());
        println!("Rigid bodies   : {}", self.simulation.n_rigid_bodies);
        println!("Sleeping frac  : {:.1}%", self.simulation.sleeping_fraction() * 100.0);
        println!("Contacts       : {}", self.simulation.n_contacts);
        if !self.constraints.is_empty() {
            println!("Avg solver iters: {:.1}", self.constraints.avg_iterations());
            println!("Max residual    : {:.2e}", self.constraints.max_residual());
        }
        if !self.energy.is_empty() {
            println!("Avg KE         : {:.4e} J", self.energy.avg_kinetic());
            println!("Energy drift   : {:.2e}", self.energy.max_energy_drift());
        }
        println!("Memory (total) : {:.2} MB", self.memory.total_mb());
        println!("======================================");
    }

    /// Serialise the summary to a JSON-like string.
    ///
    /// This is a hand-rolled serialiser that avoids pulling in `serde` as a
    /// hard dependency of the facade crate.
    pub fn to_json(&self) -> String {
        let sim = &self.simulation;
        let mem = &self.memory;
        format!(
            r#"{{
  "step_count": {},
  "sim_time": {},
  "wall_time": {},
  "real_time_factor": {},
  "n_rigid_bodies": {},
  "n_constraints": {},
  "n_contacts": {},
  "n_sleeping": {},
  "sleeping_fraction": {},
  "avg_kinetic_energy": {},
  "max_energy_drift": {},
  "avg_solver_iterations": {},
  "max_residual": {},
  "broadphase_pairs": {},
  "narrowphase_checks": {},
  "contacts_generated": {},
  "ccd_steps": {},
  "memory_total_mb": {}
}}"#,
            sim.step_count,
            sim.sim_time,
            sim.wall_time,
            sim.real_time_factor(),
            sim.n_rigid_bodies,
            sim.n_constraints,
            sim.n_contacts,
            sim.n_sleeping,
            sim.sleeping_fraction(),
            self.energy.avg_kinetic(),
            self.energy.max_energy_drift(),
            self.constraints.avg_iterations(),
            self.constraints.max_residual(),
            self.collision.broadphase_pairs,
            self.collision.narrowphase_checks,
            self.collision.contacts_generated,
            self.collision.ccd_steps,
            mem.total_mb(),
        )
    }
}

// ---------------------------------------------------------------------------
// StepLogger
// ---------------------------------------------------------------------------

/// Per-step diagnostics logger with a rolling buffer.
pub struct StepLogger {
    /// The rolling buffer of per-step diagnostics.
    buffer: VecDeque<SimulationDiagnostics>,
    /// Maximum number of entries to keep.
    capacity: usize,
    /// Accumulated total wall-clock time across all logged steps.
    total_wall_time: f64,
    /// Number of steps ever logged (including ones that fell off the buffer).
    total_logged: u64,
}

impl StepLogger {
    /// Create a new step logger with the given rolling-buffer capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: VecDeque::with_capacity(capacity),
            capacity,
            total_wall_time: 0.0,
            total_logged: 0,
        }
    }

    /// Log diagnostics for a single simulation step.
    pub fn log_step(&mut self, diag: SimulationDiagnostics) {
        self.total_wall_time += diag.wall_time;
        self.total_logged += 1;
        if self.buffer.len() >= self.capacity {
            self.buffer.pop_front();
        }
        self.buffer.push_back(diag);
    }

    /// Return up to the last `n` diagnostics snapshots (oldest first).
    pub fn get_last_n(&self, n: usize) -> Vec<&SimulationDiagnostics> {
        let skip = if self.buffer.len() > n {
            self.buffer.len() - n
        } else {
            0
        };
        self.buffer.iter().skip(skip).collect()
    }

    /// Return the average step wall-clock time in milliseconds.
    ///
    /// Uses all steps ever logged, not just the ones in the rolling buffer.
    pub fn average_step_time_ms(&self) -> f64 {
        if self.total_logged == 0 {
            0.0
        } else {
            (self.total_wall_time / self.total_logged as f64) * 1000.0
        }
    }

    /// Return the number of entries currently in the rolling buffer.
    pub fn len(&self) -> usize {
        self.buffer.len()
    }

    /// Return `true` if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }

    /// Return the total number of steps ever logged.
    pub fn total_logged(&self) -> u64 {
        self.total_logged
    }

    /// Clear the rolling buffer (does not reset `total_logged`).
    pub fn clear(&mut self) {
        self.buffer.clear();
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PhysicsTimer ---

    #[test]
    fn test_timer_new_not_running() {
        let t = PhysicsTimer::new();
        assert!(!t.is_running());
    }

    #[test]
    fn test_timer_start_sets_running() {
        let mut t = PhysicsTimer::new();
        t.start();
        assert!(t.is_running());
    }

    #[test]
    fn test_timer_stop_clears_running() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.stop();
        assert!(!t.is_running());
    }

    #[test]
    fn test_timer_elapsed_zero_when_not_started() {
        let t = PhysicsTimer::new();
        assert_eq!(t.elapsed(), Duration::ZERO);
    }

    #[test]
    fn test_timer_elapsed_ms_zero_when_not_started() {
        let t = PhysicsTimer::new();
        assert!((t.elapsed_ms() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_timer_elapsed_increases_after_start() {
        let mut t = PhysicsTimer::new();
        t.start();
        // Just check that elapsed is non-negative; no actual sleep needed.
        let e = t.elapsed_ms();
        assert!(e >= 0.0);
    }

    #[test]
    fn test_timer_lap_records_time() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.lap();
        t.stop();
        assert_eq!(t.split_times().len(), 1);
    }

    #[test]
    fn test_timer_multiple_laps() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.lap();
        t.lap();
        t.lap();
        t.stop();
        assert_eq!(t.split_times().len(), 3);
    }

    #[test]
    fn test_timer_reset() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.lap();
        t.stop();
        t.reset();
        assert!(!t.is_running());
        assert_eq!(t.split_times().len(), 0);
        assert_eq!(t.elapsed(), Duration::ZERO);
    }

    #[test]
    fn test_timer_split_separate_from_lap() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.lap();
        t.split();
        t.stop();
        assert_eq!(t.split_times().len(), 1); // lap
        assert_eq!(t.splits().len(), 1); // split
    }

    #[test]
    fn test_timer_double_start_no_op() {
        let mut t = PhysicsTimer::new();
        t.start();
        let _e1 = t.elapsed();
        t.start(); // second start should be ignored
        assert!(t.is_running());
    }

    #[test]
    fn test_timer_double_stop_no_op() {
        let mut t = PhysicsTimer::new();
        t.start();
        t.stop();
        t.stop(); // second stop should be ignored
        assert!(!t.is_running());
    }

    // --- SimulationDiagnostics ---

    #[test]
    fn test_simulation_diag_defaults() {
        let d = SimulationDiagnostics::default();
        assert_eq!(d.step_count, 0);
        assert!((d.sim_time - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_sleeping_fraction_no_bodies() {
        let d = SimulationDiagnostics::default();
        assert!((d.sleeping_fraction() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_sleeping_fraction_half() {
        let d = SimulationDiagnostics {
            n_rigid_bodies: 10,
            n_sleeping: 5,
            ..Default::default()
        };
        assert!((d.sleeping_fraction() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_real_time_factor_zero_wall_time() {
        let d = SimulationDiagnostics::default();
        assert!((d.real_time_factor() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_real_time_factor_one() {
        let d = SimulationDiagnostics {
            sim_time: 1.0,
            wall_time: 1.0,
            ..Default::default()
        };
        assert!((d.real_time_factor() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_real_time_factor_faster() {
        let d = SimulationDiagnostics {
            sim_time: 2.0,
            wall_time: 1.0,
            ..Default::default()
        };
        assert!((d.real_time_factor() - 2.0).abs() < 1e-12);
    }

    // --- EnergyMonitor ---

    #[test]
    fn test_energy_monitor_empty() {
        let m = EnergyMonitor::new();
        assert!(m.is_empty());
        assert_eq!(m.len(), 0);
    }

    #[test]
    fn test_energy_monitor_record() {
        let mut m = EnergyMonitor::new();
        m.record(10.0, 5.0, 0.0);
        assert_eq!(m.len(), 1);
        assert!((m.total[0] - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_avg_kinetic() {
        let mut m = EnergyMonitor::new();
        m.record(10.0, 0.0, 0.0);
        m.record(20.0, 0.0, 0.0);
        assert!((m.avg_kinetic() - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_avg_potential() {
        let mut m = EnergyMonitor::new();
        m.record(0.0, 4.0, 0.0);
        m.record(0.0, 8.0, 0.0);
        assert!((m.avg_potential() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_max_drift_less_than_two() {
        let mut m = EnergyMonitor::new();
        m.record(10.0, 5.0, 0.0);
        assert!((m.max_energy_drift() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_max_drift_stable() {
        let mut m = EnergyMonitor::new();
        m.record(10.0, 5.0, 0.0);
        m.record(10.0, 5.0, 0.0);
        assert!((m.max_energy_drift() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_monitor_reset() {
        let mut m = EnergyMonitor::new();
        m.record(10.0, 5.0, 0.0);
        m.reset();
        assert!(m.is_empty());
    }

    // --- ConstraintDiagnostics ---

    #[test]
    fn test_constraint_diag_empty() {
        let d = ConstraintDiagnostics::new();
        assert!(d.is_empty());
        assert!((d.avg_iterations() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_constraint_diag_record() {
        let mut d = ConstraintDiagnostics::new();
        d.record(10, 0.001);
        assert_eq!(d.len(), 1);
        assert!((d.avg_iterations() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_constraint_diag_avg_iterations() {
        let mut d = ConstraintDiagnostics::new();
        d.record(8, 0.01);
        d.record(12, 0.02);
        assert!((d.avg_iterations() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_constraint_diag_max_residual() {
        let mut d = ConstraintDiagnostics::new();
        d.record(5, 0.01);
        d.record(5, 0.05);
        d.record(5, 0.02);
        assert!((d.max_residual() - 0.05).abs() < 1e-12);
    }

    // --- CollisionDiagnostics ---

    #[test]
    fn test_collision_diag_defaults() {
        let d = CollisionDiagnostics::new();
        assert_eq!(d.broadphase_pairs, 0);
        assert!((d.narrowphase_ratio() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_collision_diag_narrowphase_ratio() {
        let d = CollisionDiagnostics {
            broadphase_pairs: 100,
            narrowphase_checks: 50,
            ..Default::default()
        };
        assert!((d.narrowphase_ratio() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_collision_diag_reset() {
        let mut d = CollisionDiagnostics {
            broadphase_pairs: 100,
            ..Default::default()
        };
        d.reset();
        assert_eq!(d.broadphase_pairs, 0);
    }

    // --- MemoryStats ---

    #[test]
    fn test_memory_stats_total_zero() {
        let m = MemoryStats::new();
        assert_eq!(m.total_bytes(), 0);
        assert!((m.total_mb() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_memory_stats_total() {
        let m = MemoryStats {
            rigid_bytes: 1024,
            collision_bytes: 512,
            ..Default::default()
        };
        assert_eq!(m.total_bytes(), 1536);
    }

    #[test]
    fn test_memory_stats_mb() {
        let m = MemoryStats {
            rigid_bytes: 1024 * 1024,
            ..Default::default()
        };
        assert!((m.total_mb() - 1.0).abs() < 1e-6);
    }

    // --- PerformanceSummary ---

    #[test]
    fn test_performance_summary_defaults() {
        let s = PerformanceSummary::new();
        assert_eq!(s.simulation.step_count, 0);
    }

    #[test]
    fn test_performance_summary_to_json_contains_step_count() {
        let mut s = PerformanceSummary::new();
        s.simulation.step_count = 42;
        let json = s.to_json();
        assert!(json.contains("42"));
    }

    #[test]
    fn test_performance_summary_to_json_is_string() {
        let s = PerformanceSummary::new();
        let json = s.to_json();
        assert!(!json.is_empty());
        assert!(json.starts_with('{'));
        assert!(json.ends_with('}'));
    }

    #[test]
    fn test_performance_summary_print_does_not_panic() {
        let mut s = PerformanceSummary::new();
        s.simulation.step_count = 100;
        s.simulation.sim_time = 1.0;
        s.simulation.wall_time = 1.1;
        s.simulation.n_rigid_bodies = 5;
        s.energy.record(10.0, 5.0, 0.0);
        s.constraints.record(8, 0.001);
        // Just ensure it doesn't panic
        s.print_summary();
    }

    // --- StepLogger ---

    #[test]
    fn test_step_logger_empty() {
        let logger = StepLogger::new(10);
        assert!(logger.is_empty());
        assert_eq!(logger.len(), 0);
        assert_eq!(logger.total_logged(), 0);
    }

    #[test]
    fn test_step_logger_log_one() {
        let mut logger = StepLogger::new(10);
        let d = SimulationDiagnostics {
            step_count: 1,
            wall_time: 0.016,
            ..Default::default()
        };
        logger.log_step(d);
        assert_eq!(logger.len(), 1);
        assert_eq!(logger.total_logged(), 1);
    }

    #[test]
    fn test_step_logger_rolling_overflow() {
        let mut logger = StepLogger::new(3);
        for i in 0..5 {
            let d = SimulationDiagnostics {
                step_count: i,
                ..Default::default()
            };
            logger.log_step(d);
        }
        assert_eq!(logger.len(), 3); // Only last 3 kept
        assert_eq!(logger.total_logged(), 5); // But total is 5
    }

    #[test]
    fn test_step_logger_get_last_n() {
        let mut logger = StepLogger::new(10);
        for i in 0..5 {
            let d = SimulationDiagnostics {
                step_count: i,
                ..Default::default()
            };
            logger.log_step(d);
        }
        let last2 = logger.get_last_n(2);
        assert_eq!(last2.len(), 2);
        assert_eq!(last2[1].step_count, 4); // Last one
    }

    #[test]
    fn test_step_logger_average_step_time_ms() {
        let mut logger = StepLogger::new(10);
        for _ in 0..4 {
            let d = SimulationDiagnostics {
                wall_time: 0.01, // 10 ms in seconds
                ..Default::default()
            };
            logger.log_step(d);
        }
        // average_step_time_ms = (total_wall_time / total_logged) * 1000
        // = (0.04 / 4) * 1000 = 10.0 ms
        let avg = logger.average_step_time_ms();
        assert!((avg - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_step_logger_average_empty() {
        let logger = StepLogger::new(10);
        assert!((logger.average_step_time_ms() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_step_logger_clear() {
        let mut logger = StepLogger::new(10);
        logger.log_step(SimulationDiagnostics::default());
        logger.log_step(SimulationDiagnostics::default());
        logger.clear();
        assert!(logger.is_empty());
        assert_eq!(logger.total_logged(), 2); // Still tracks total
    }

    #[test]
    fn test_step_logger_get_last_n_more_than_available() {
        let mut logger = StepLogger::new(10);
        logger.log_step(SimulationDiagnostics::default());
        let last = logger.get_last_n(100);
        assert_eq!(last.len(), 1);
    }
}
