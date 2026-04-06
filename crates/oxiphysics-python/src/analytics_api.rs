#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics analytics Python bindings.
//!
//! Provides Python-friendly types for measuring and analysing simulation
//! quality: energy tracking, momentum conservation, collision statistics,
//! performance profiling, and full simulation reports.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ─────────────────────────────────────────────────────────────────────────────
// Snapshot
// ─────────────────────────────────────────────────────────────────────────────

/// A snapshot of global physics state at one time step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySnapshot {
    /// Simulation time (seconds).
    pub time: f64,
    /// Total kinetic energy (J).
    pub kinetic_energy: f64,
    /// Total potential energy (J).
    pub potential_energy: f64,
    /// Total energy (KE + PE).
    pub total_energy: f64,
    /// Linear momentum magnitude `|p|`.
    pub linear_momentum: f64,
    /// Angular momentum magnitude `|L|`.
    pub angular_momentum: f64,
    /// Number of active contacts.
    pub contact_count: usize,
    /// Number of awake bodies.
    pub awake_count: usize,
    /// Number of sleeping bodies.
    pub sleeping_count: usize,
    /// Step wall-clock time in milliseconds.
    pub step_ms: f64,
}

impl Default for PySnapshot {
    fn default() -> Self {
        Self {
            time: 0.0,
            kinetic_energy: 0.0,
            potential_energy: 0.0,
            total_energy: 0.0,
            linear_momentum: 0.0,
            angular_momentum: 0.0,
            contact_count: 0,
            awake_count: 0,
            sleeping_count: 0,
            step_ms: 0.0,
        }
    }
}

impl PySnapshot {
    /// Compute the total energy.
    pub fn compute_total(&mut self) {
        self.total_energy = self.kinetic_energy + self.potential_energy;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyPhysicsAnalytics
// ─────────────────────────────────────────────────────────────────────────────

/// Master analytics object for a simulation.
///
/// Holds snapshot history and delegates to sub-trackers.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyPhysicsAnalytics {
    /// All recorded snapshots.
    pub snapshots: Vec<PySnapshot>,
    /// Maximum number of snapshots to retain.
    pub max_snapshots: usize,
    /// Current step index.
    pub step: u64,
}

impl PyPhysicsAnalytics {
    /// Create analytics with a given history capacity.
    pub fn new(max_snapshots: usize) -> Self {
        Self {
            snapshots: Vec::new(),
            max_snapshots,
            step: 0,
        }
    }

    /// Record a new snapshot.
    pub fn snapshot(&mut self, snap: PySnapshot) {
        if self.max_snapshots > 0 && self.snapshots.len() >= self.max_snapshots {
            self.snapshots.remove(0);
        }
        self.snapshots.push(snap);
        self.step += 1;
    }

    /// Advance step counter.
    pub fn step(&mut self) {
        self.step += 1;
    }

    /// Compute aggregate metrics over all snapshots.
    pub fn compute_metrics(&self) -> PyAggregateMetrics {
        if self.snapshots.is_empty() {
            return PyAggregateMetrics::default();
        }
        let n = self.snapshots.len();
        let mean_ke = self.snapshots.iter().map(|s| s.kinetic_energy).sum::<f64>() / n as f64;
        let mean_pe = self
            .snapshots
            .iter()
            .map(|s| s.potential_energy)
            .sum::<f64>()
            / n as f64;
        let max_ke = self
            .snapshots
            .iter()
            .map(|s| s.kinetic_energy)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_ke = self
            .snapshots
            .iter()
            .map(|s| s.kinetic_energy)
            .fold(f64::INFINITY, f64::min);
        let mean_step_ms = self.snapshots.iter().map(|s| s.step_ms).sum::<f64>() / n as f64;
        let max_step_ms = self
            .snapshots
            .iter()
            .map(|s| s.step_ms)
            .fold(f64::NEG_INFINITY, f64::max);
        PyAggregateMetrics {
            mean_ke,
            mean_pe,
            max_ke,
            min_ke,
            mean_step_ms,
            max_step_ms,
            sample_count: n,
        }
    }

    /// Most recent snapshot.
    pub fn latest(&self) -> Option<&PySnapshot> {
        self.snapshots.last()
    }

    /// True if no snapshots recorded.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }
}

/// Aggregate metrics derived from snapshot history.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyAggregateMetrics {
    /// Mean kinetic energy.
    pub mean_ke: f64,
    /// Mean potential energy.
    pub mean_pe: f64,
    /// Maximum kinetic energy observed.
    pub max_ke: f64,
    /// Minimum kinetic energy observed.
    pub min_ke: f64,
    /// Mean step wall-clock time (ms).
    pub mean_step_ms: f64,
    /// Maximum step wall-clock time (ms).
    pub max_step_ms: f64,
    /// Number of snapshots in this aggregate.
    pub sample_count: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// PyEnergyTracker
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks kinetic, potential, and total energy over time.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyEnergyTracker {
    /// Time series of (time, KE, PE, total).
    pub history: Vec<(f64, f64, f64, f64)>,
    /// Maximum history size.
    pub capacity: usize,
    /// Jump threshold for anomaly detection (energy change ratio).
    pub anomaly_threshold: f64,
    /// Detected energy anomalies `(step_index, delta_energy)`.
    pub anomalies: Vec<(usize, f64)>,
}

impl PyEnergyTracker {
    /// Create a new tracker.
    pub fn new(capacity: usize, anomaly_threshold: f64) -> Self {
        Self {
            history: Vec::new(),
            capacity,
            anomaly_threshold,
            anomalies: Vec::new(),
        }
    }

    /// Record a step's energy.
    pub fn record(&mut self, time: f64, ke: f64, pe: f64) {
        let total = ke + pe;
        if let Some(&(_t0, _ke0, _pe0, prev_total)) = self.history.last()
            && prev_total.abs() > 1e-10
        {
            let delta = (total - prev_total).abs() / prev_total.abs();
            if delta > self.anomaly_threshold {
                self.anomalies
                    .push((self.history.len(), total - prev_total));
            }
        }
        if self.capacity > 0 && self.history.len() >= self.capacity {
            self.history.remove(0);
        }
        self.history.push((time, ke, pe, total));
    }

    /// Mean kinetic energy.
    pub fn mean_ke(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        self.history.iter().map(|&(_, ke, _, _)| ke).sum::<f64>() / self.history.len() as f64
    }

    /// Mean total energy.
    pub fn mean_total(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        self.history
            .iter()
            .map(|&(_, _, _, total)| total)
            .sum::<f64>()
            / self.history.len() as f64
    }

    /// Latest total energy (or 0).
    pub fn latest_total(&self) -> f64 {
        self.history.last().map_or(0.0, |&(_, _, _, t)| t)
    }

    /// True if any anomalies were detected.
    pub fn has_anomalies(&self) -> bool {
        !self.anomalies.is_empty()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyMomentumTracker
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks linear and angular momentum conservation.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyMomentumTracker {
    /// History of `(time, lin_mom_mag, ang_mom_mag)`.
    pub history: Vec<(f64, f64, f64)>,
    /// Initial linear momentum for conservation check.
    pub initial_linear: f64,
    /// Initial angular momentum for conservation check.
    pub initial_angular: f64,
    /// Tolerance for conservation check.
    pub tolerance: f64,
}

impl PyMomentumTracker {
    /// Create with default tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }

    /// Record a step.
    pub fn record(&mut self, time: f64, linear_mag: f64, angular_mag: f64) {
        if self.history.is_empty() {
            self.initial_linear = linear_mag;
            self.initial_angular = angular_mag;
        }
        self.history.push((time, linear_mag, angular_mag));
    }

    /// True if linear momentum is conserved within tolerance.
    pub fn linear_conserved(&self) -> bool {
        if self.history.is_empty() {
            return true;
        }
        self.history.iter().all(|&(_, lm, _)| {
            (lm - self.initial_linear).abs() <= self.tolerance * (self.initial_linear.abs() + 1.0)
        })
    }

    /// True if angular momentum is conserved within tolerance.
    pub fn angular_conserved(&self) -> bool {
        if self.history.is_empty() {
            return true;
        }
        self.history.iter().all(|&(_, _, am)| {
            (am - self.initial_angular).abs() <= self.tolerance * (self.initial_angular.abs() + 1.0)
        })
    }

    /// Maximum linear momentum deviation.
    pub fn max_linear_deviation(&self) -> f64 {
        self.history
            .iter()
            .map(|&(_, lm, _)| (lm - self.initial_linear).abs())
            .fold(0.0_f64, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyCollisionStats
// ─────────────────────────────────────────────────────────────────────────────

/// Statistics about collisions over the simulation.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyCollisionStats {
    /// Total collision count.
    pub total_collisions: u64,
    /// Per-body-pair collision count (body_a_id, body_b_id, count).
    pub pair_counts: Vec<(u32, u32, u64)>,
    /// Impulse histogram (n_bins buckets in \[0, max_impulse\]).
    pub impulse_histogram: Vec<u64>,
    /// Maximum impulse bin value.
    pub max_impulse: f64,
    /// Penetration depth statistics.
    pub depth_min: f64,
    /// Maximum penetration depth recorded.
    pub depth_max: f64,
    /// Mean penetration depth.
    pub depth_mean: f64,
    depth_sum: f64,
    depth_count: u64,
}

impl PyCollisionStats {
    /// Create with given histogram resolution.
    pub fn new(n_bins: usize, max_impulse: f64) -> Self {
        Self {
            max_impulse,
            impulse_histogram: vec![0; n_bins],
            depth_min: f64::INFINITY,
            depth_max: f64::NEG_INFINITY,
            ..Default::default()
        }
    }

    /// Record a collision event.
    pub fn record_collision(&mut self, body_a: u32, body_b: u32, impulse: f64, depth: f64) {
        self.total_collisions += 1;
        // Pair counts
        let pair = self
            .pair_counts
            .iter_mut()
            .find(|(a, b, _)| (*a == body_a && *b == body_b) || (*a == body_b && *b == body_a));
        match pair {
            Some((_, _, cnt)) => *cnt += 1,
            None => self.pair_counts.push((body_a, body_b, 1)),
        }
        // Histogram
        let n = self.impulse_histogram.len();
        if n > 0 && self.max_impulse > 0.0 {
            let bin = ((impulse / self.max_impulse) * n as f64) as usize;
            let bin = bin.min(n - 1);
            self.impulse_histogram[bin] += 1;
        }
        // Depth stats
        if depth < self.depth_min {
            self.depth_min = depth;
        }
        if depth > self.depth_max {
            self.depth_max = depth;
        }
        self.depth_sum += depth;
        self.depth_count += 1;
        self.depth_mean = self.depth_sum / self.depth_count as f64;
    }

    /// Reset all stats.
    pub fn reset(&mut self) {
        self.total_collisions = 0;
        self.pair_counts.clear();
        for b in &mut self.impulse_histogram {
            *b = 0;
        }
        self.depth_min = f64::INFINITY;
        self.depth_max = f64::NEG_INFINITY;
        self.depth_mean = 0.0;
        self.depth_sum = 0.0;
        self.depth_count = 0;
    }

    /// Most collided pair (if any).
    pub fn most_collided_pair(&self) -> Option<(u32, u32, u64)> {
        self.pair_counts.iter().max_by_key(|&&(_, _, c)| c).copied()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySleepStats
// ─────────────────────────────────────────────────────────────────────────────

/// Statistics about body sleep/wake behaviour.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PySleepStats {
    /// Total wake events observed.
    pub total_wake_events: u64,
    /// Total sleep events observed.
    pub total_sleep_events: u64,
    /// History of (time, sleeping_fraction) tuples.
    pub sleeping_fraction_history: Vec<(f64, f64)>,
    /// Capacity limit.
    pub capacity: usize,
}

impl PySleepStats {
    /// Create with history capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            ..Default::default()
        }
    }

    /// Record a step with body counts.
    pub fn record(&mut self, time: f64, sleeping: usize, total: usize, wakes: u64, sleeps: u64) {
        self.total_wake_events += wakes;
        self.total_sleep_events += sleeps;
        let frac = if total == 0 {
            0.0
        } else {
            sleeping as f64 / total as f64
        };
        if self.capacity > 0 && self.sleeping_fraction_history.len() >= self.capacity {
            self.sleeping_fraction_history.remove(0);
        }
        self.sleeping_fraction_history.push((time, frac));
    }

    /// Mean sleeping fraction.
    pub fn mean_sleeping_fraction(&self) -> f64 {
        if self.sleeping_fraction_history.is_empty() {
            return 0.0;
        }
        self.sleeping_fraction_history
            .iter()
            .map(|&(_, f)| f)
            .sum::<f64>()
            / self.sleeping_fraction_history.len() as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyConstraintStats
// ─────────────────────────────────────────────────────────────────────────────

/// Solver and constraint quality statistics.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyConstraintStats {
    /// Mean solver residual per step.
    pub mean_residual: f64,
    /// Maximum residual observed.
    pub max_residual: f64,
    /// Mean solver iterations used.
    pub mean_iterations: f64,
    /// Warm-start hit rate (fraction of constraints using warm start).
    pub warm_start_rate: f64,
    /// Number of samples.
    sample_count: u64,
    residual_sum: f64,
    iterations_sum: f64,
}

impl PyConstraintStats {
    /// Create empty stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a solver step.
    pub fn record(&mut self, residual: f64, iterations: u32, warm_start_fraction: f64) {
        self.sample_count += 1;
        self.residual_sum += residual;
        self.iterations_sum += iterations as f64;
        self.mean_residual = self.residual_sum / self.sample_count as f64;
        self.mean_iterations = self.iterations_sum / self.sample_count as f64;
        self.warm_start_rate = warm_start_fraction; // latest value
        if residual > self.max_residual {
            self.max_residual = residual;
        }
    }

    /// True if the solver is converging well.
    pub fn is_converging(&self, threshold: f64) -> bool {
        self.mean_residual < threshold
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyPerformanceTracker
// ─────────────────────────────────────────────────────────────────────────────

/// Breakdown of time spent per pipeline phase (milliseconds).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyStepTimings {
    /// Broadphase collision detection (ms).
    pub broad_ms: f64,
    /// Narrowphase collision detection (ms).
    pub narrow_ms: f64,
    /// Constraint solver (ms).
    pub solver_ms: f64,
    /// Integration (ms).
    pub integration_ms: f64,
    /// Total step time (ms).
    pub total_ms: f64,
}

impl PyStepTimings {
    /// Compute total as sum of phases.
    pub fn compute_total(&mut self) {
        self.total_ms = self.broad_ms + self.narrow_ms + self.solver_ms + self.integration_ms;
    }
}

/// Tracks per-step performance timings.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyPerformanceTracker {
    /// History of per-step timings.
    pub history: Vec<PyStepTimings>,
    /// Capacity limit.
    pub capacity: usize,
}

impl PyPerformanceTracker {
    /// Create with capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            history: Vec::new(),
        }
    }

    /// Record a step's timings.
    pub fn record(&mut self, mut timings: PyStepTimings) {
        timings.compute_total();
        if self.capacity > 0 && self.history.len() >= self.capacity {
            self.history.remove(0);
        }
        self.history.push(timings);
    }

    /// Rolling average total step time.
    pub fn mean_total_ms(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        self.history.iter().map(|t| t.total_ms).sum::<f64>() / self.history.len() as f64
    }

    /// Peak total step time.
    pub fn peak_total_ms(&self) -> f64 {
        self.history
            .iter()
            .map(|t| t.total_ms)
            .fold(0.0_f64, f64::max)
    }

    /// Export timings as CSV string.
    pub fn export_csv(&self) -> String {
        let mut out = String::from("broad_ms,narrow_ms,solver_ms,integration_ms,total_ms\n");
        for t in &self.history {
            out.push_str(&format!(
                "{:.4},{:.4},{:.4},{:.4},{:.4}\n",
                t.broad_ms, t.narrow_ms, t.solver_ms, t.integration_ms, t.total_ms
            ));
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PySimulationReport
// ─────────────────────────────────────────────────────────────────────────────

/// A full simulation quality report.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PySimulationReport {
    /// Total steps simulated.
    pub total_steps: u64,
    /// Total simulation time (s).
    pub total_time: f64,
    /// Energy summary.
    pub energy: PyEnergyTracker,
    /// Momentum summary.
    pub momentum: PyMomentumTracker,
    /// Collision stats.
    pub collision: PyCollisionStats,
    /// Sleep stats.
    pub sleep: PySleepStats,
    /// Constraint stats.
    pub constraint: PyConstraintStats,
    /// Performance stats.
    pub performance: PyPerformanceTracker,
    /// List of anomalies as human-readable strings.
    pub anomaly_list: Vec<String>,
}

impl PySimulationReport {
    /// Create a new report.
    pub fn new() -> Self {
        Self {
            energy: PyEnergyTracker::new(1000, 0.5),
            momentum: PyMomentumTracker::new(0.01),
            collision: PyCollisionStats::new(32, 1000.0),
            sleep: PySleepStats::new(1000),
            constraint: PyConstraintStats::new(),
            performance: PyPerformanceTracker::new(1000),
            ..Default::default()
        }
    }

    /// Serialise report to JSON string.
    pub fn as_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Append an anomaly message.
    pub fn flag_anomaly(&mut self, msg: impl Into<String>) {
        self.anomaly_list.push(msg.into());
    }

    /// Run a check and flag if energy anomalies exist.
    pub fn check_energy(&mut self) {
        if self.energy.has_anomalies() {
            self.flag_anomaly(format!(
                "Energy anomalies detected: {} jumps",
                self.energy.anomalies.len()
            ));
        }
    }

    /// Run a check and flag if momentum is not conserved.
    pub fn check_momentum(&mut self) {
        if !self.momentum.linear_conserved() {
            self.flag_anomaly(format!(
                "Linear momentum not conserved; max deviation = {:.4}",
                self.momentum.max_linear_deviation()
            ));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyBenchmark
// ─────────────────────────────────────────────────────────────────────────────

/// Benchmarks a physics simulation over N steps.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyBenchmark {
    /// Number of steps run.
    pub steps_run: u64,
    /// Total wall-clock time for all steps (ms, simulated).
    pub total_ms: f64,
    /// Throughput in steps per second.
    pub throughput: f64,
    /// Body count used in the benchmark.
    pub body_count: usize,
    /// Constraint count used.
    pub constraint_count: usize,
}

impl PyBenchmark {
    /// Create a benchmark for the given body and constraint counts.
    pub fn new(body_count: usize, constraint_count: usize) -> Self {
        Self {
            body_count,
            constraint_count,
            ..Default::default()
        }
    }

    /// Run `n_steps` of a mock simulation with artificial cost.
    pub fn run(&mut self, n_steps: u64, step_cost_ms: f64) {
        self.steps_run = n_steps;
        // Mock: cost scales linearly with body+constraint count
        let scale = 1.0 + (self.body_count + self.constraint_count) as f64 / 1000.0;
        self.total_ms = n_steps as f64 * step_cost_ms * scale;
        if self.total_ms > 0.0 {
            self.throughput = self.steps_run as f64 / (self.total_ms / 1000.0);
        }
    }

    /// Throughput in steps/s.
    pub fn steps_per_second(&self) -> f64 {
        self.throughput
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PyDataExporter
// ─────────────────────────────────────────────────────────────────────────────

/// Exports simulation state to CSV or JSON.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PyDataExporter {
    /// Recorded rows: (time, body_id, px, py, pz, vx, vy, vz, ke).
    rows: Vec<(f64, u32, f64, f64, f64, f64, f64, f64, f64)>,
}

impl PyDataExporter {
    /// Create an empty exporter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a body state.
    #[allow(clippy::too_many_arguments)]
    pub fn record(
        &mut self,
        time: f64,
        body_id: u32,
        position: [f64; 3],
        velocity: [f64; 3],
        ke: f64,
    ) {
        self.rows.push((
            time,
            body_id,
            position[0],
            position[1],
            position[2],
            velocity[0],
            velocity[1],
            velocity[2],
            ke,
        ));
    }

    /// Export all rows to CSV.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("time,body_id,px,py,pz,vx,vy,vz,ke\n");
        for &(t, id, px, py, pz, vx, vy, vz, ke) in &self.rows {
            out.push_str(&format!(
                "{t:.6},{id},{px:.6},{py:.6},{pz:.6},{vx:.6},{vy:.6},{vz:.6},{ke:.6}\n"
            ));
        }
        out
    }

    /// Export all rows to JSON.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.rows)
    }

    /// Number of recorded rows.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Clear all rows.
    pub fn clear(&mut self) {
        self.rows.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- PySnapshot ---

    #[test]
    fn test_snapshot_default() {
        let s = PySnapshot::default();
        assert_eq!(s.time, 0.0);
        assert_eq!(s.kinetic_energy, 0.0);
    }

    #[test]
    fn test_snapshot_compute_total() {
        let mut s = PySnapshot {
            kinetic_energy: 3.0,
            potential_energy: 2.0,
            ..Default::default()
        };
        s.compute_total();
        assert!((s.total_energy - 5.0).abs() < 1e-12);
    }

    // --- PyPhysicsAnalytics ---

    #[test]
    fn test_analytics_snapshot() {
        let mut a = PyPhysicsAnalytics::new(10);
        for i in 0..5 {
            let s = PySnapshot {
                time: i as f64,
                kinetic_energy: i as f64,
                ..Default::default()
            };
            a.snapshot(s);
        }
        assert_eq!(a.snapshots.len(), 5);
    }

    #[test]
    fn test_analytics_evicts_old() {
        let mut a = PyPhysicsAnalytics::new(3);
        for i in 0..5 {
            let s = PySnapshot {
                kinetic_energy: i as f64,
                ..Default::default()
            };
            a.snapshot(s);
        }
        assert_eq!(a.snapshots.len(), 3);
    }

    #[test]
    fn test_analytics_metrics() {
        let mut a = PyPhysicsAnalytics::new(100);
        for i in 1..=4 {
            let s = PySnapshot {
                kinetic_energy: i as f64,
                ..Default::default()
            };
            a.snapshot(s);
        }
        let m = a.compute_metrics();
        assert!((m.mean_ke - 2.5).abs() < 1e-10);
        assert!((m.max_ke - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_analytics_is_empty() {
        let a = PyPhysicsAnalytics::new(10);
        assert!(a.is_empty());
    }

    // --- PyEnergyTracker ---

    #[test]
    fn test_energy_tracker_record() {
        let mut t = PyEnergyTracker::new(100, 0.5);
        t.record(0.0, 10.0, 5.0);
        t.record(0.1, 11.0, 5.0);
        assert_eq!(t.history.len(), 2);
    }

    #[test]
    fn test_energy_tracker_anomaly_detect() {
        let mut t = PyEnergyTracker::new(100, 0.1);
        t.record(0.0, 10.0, 0.0);
        t.record(0.1, 1000.0, 0.0); // huge jump
        assert!(t.has_anomalies());
    }

    #[test]
    fn test_energy_tracker_no_anomaly() {
        let mut t = PyEnergyTracker::new(100, 0.5);
        t.record(0.0, 10.0, 0.0);
        t.record(0.1, 10.1, 0.0);
        assert!(!t.has_anomalies());
    }

    #[test]
    fn test_energy_tracker_mean_ke() {
        let mut t = PyEnergyTracker::new(100, 1.0);
        t.record(0.0, 4.0, 0.0);
        t.record(0.1, 6.0, 0.0);
        assert!((t.mean_ke() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_energy_tracker_eviction() {
        let mut t = PyEnergyTracker::new(2, 10.0);
        t.record(0.0, 1.0, 0.0);
        t.record(0.1, 2.0, 0.0);
        t.record(0.2, 3.0, 0.0);
        assert_eq!(t.history.len(), 2);
    }

    // --- PyMomentumTracker ---

    #[test]
    fn test_momentum_tracker_conserved() {
        let mut t = PyMomentumTracker::new(0.01);
        t.record(0.0, 5.0, 2.0);
        t.record(0.1, 5.001, 2.001);
        assert!(t.linear_conserved());
    }

    #[test]
    fn test_momentum_tracker_not_conserved() {
        let mut t = PyMomentumTracker::new(0.001);
        t.record(0.0, 5.0, 2.0);
        t.record(0.1, 6.0, 2.0); // large jump
        assert!(!t.linear_conserved());
    }

    #[test]
    fn test_momentum_tracker_deviation() {
        let mut t = PyMomentumTracker::new(0.5);
        t.record(0.0, 10.0, 0.0);
        t.record(0.1, 10.5, 0.0);
        assert!((t.max_linear_deviation() - 0.5).abs() < 1e-10);
    }

    // --- PyCollisionStats ---

    #[test]
    fn test_collision_stats_record() {
        let mut s = PyCollisionStats::new(10, 100.0);
        s.record_collision(0, 1, 50.0, 0.01);
        assert_eq!(s.total_collisions, 1);
    }

    #[test]
    fn test_collision_stats_pair_count() {
        let mut s = PyCollisionStats::new(10, 100.0);
        s.record_collision(0, 1, 10.0, 0.01);
        s.record_collision(0, 1, 20.0, 0.02);
        let pair = s.most_collided_pair().unwrap();
        assert_eq!(pair.2, 2);
    }

    #[test]
    fn test_collision_stats_depth() {
        let mut s = PyCollisionStats::new(10, 100.0);
        s.record_collision(0, 1, 10.0, 0.05);
        s.record_collision(2, 3, 10.0, 0.10);
        assert!((s.depth_mean - 0.075).abs() < 1e-10);
    }

    #[test]
    fn test_collision_stats_reset() {
        let mut s = PyCollisionStats::new(10, 100.0);
        s.record_collision(0, 1, 10.0, 0.01);
        s.reset();
        assert_eq!(s.total_collisions, 0);
    }

    // --- PySleepStats ---

    #[test]
    fn test_sleep_stats_record() {
        let mut s = PySleepStats::new(100);
        s.record(0.0, 3, 10, 0, 3);
        s.record(0.1, 5, 10, 2, 2);
        assert_eq!(s.total_wake_events, 2);
    }

    #[test]
    fn test_sleep_stats_mean_fraction() {
        let mut s = PySleepStats::new(100);
        s.record(0.0, 5, 10, 0, 0);
        s.record(0.1, 3, 10, 0, 0);
        // (0.5 + 0.3) / 2 = 0.4
        assert!((s.mean_sleeping_fraction() - 0.4).abs() < 1e-10);
    }

    // --- PyConstraintStats ---

    #[test]
    fn test_constraint_stats_record() {
        let mut s = PyConstraintStats::new();
        s.record(0.001, 5, 0.8);
        s.record(0.002, 7, 0.9);
        assert!((s.mean_iterations - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_constraint_stats_converging() {
        let mut s = PyConstraintStats::new();
        s.record(0.0001, 3, 1.0);
        assert!(s.is_converging(0.01));
    }

    // --- PyPerformanceTracker ---

    #[test]
    fn test_perf_tracker_mean() {
        let mut t = PyPerformanceTracker::new(100);
        t.record(PyStepTimings {
            broad_ms: 1.0,
            narrow_ms: 2.0,
            solver_ms: 3.0,
            integration_ms: 0.5,
            total_ms: 0.0,
        });
        t.record(PyStepTimings {
            broad_ms: 2.0,
            narrow_ms: 1.0,
            solver_ms: 2.0,
            integration_ms: 0.5,
            total_ms: 0.0,
        });
        // totals: 6.5 and 5.5 -> mean 6.0
        assert!((t.mean_total_ms() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_perf_tracker_csv() {
        let mut t = PyPerformanceTracker::new(10);
        t.record(PyStepTimings {
            broad_ms: 1.0,
            narrow_ms: 1.0,
            solver_ms: 1.0,
            integration_ms: 1.0,
            total_ms: 0.0,
        });
        let csv = t.export_csv();
        assert!(csv.contains("broad_ms"));
    }

    // --- PySimulationReport ---

    #[test]
    fn test_report_as_json() {
        let mut r = PySimulationReport::new();
        r.total_steps = 100;
        let json = r.as_json().unwrap();
        assert!(json.contains("total_steps"));
    }

    #[test]
    fn test_report_anomaly() {
        let mut r = PySimulationReport::new();
        r.energy.record(0.0, 10.0, 0.0);
        r.energy.record(0.1, 1000.0, 0.0); // triggers anomaly
        r.check_energy();
        assert!(!r.anomaly_list.is_empty());
    }

    // --- PyBenchmark ---

    #[test]
    fn test_benchmark_run() {
        let mut b = PyBenchmark::new(100, 50);
        b.run(1000, 1.0);
        assert_eq!(b.steps_run, 1000);
        assert!(b.throughput > 0.0);
    }

    #[test]
    fn test_benchmark_throughput() {
        let mut b = PyBenchmark::new(0, 0);
        b.run(1000, 1.0); // 1000 steps × 1ms/step = 1000ms total = 1s
        // throughput ~ 1000 steps/s
        assert!((b.steps_per_second() - 1000.0).abs() < 1.0);
    }

    // --- PyDataExporter ---

    #[test]
    fn test_exporter_csv() {
        let mut e = PyDataExporter::new();
        e.record(0.0, 1, [0.0, 1.0, 2.0], [0.1, 0.2, 0.3], 5.0);
        let csv = e.to_csv();
        assert!(csv.contains("time,body_id"));
        assert!(csv.contains("0.000000"));
    }

    #[test]
    fn test_exporter_json() {
        let mut e = PyDataExporter::new();
        e.record(0.0, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0);
        let json = e.to_json().unwrap();
        assert!(json.contains("["));
    }

    #[test]
    fn test_exporter_row_count() {
        let mut e = PyDataExporter::new();
        for i in 0..5 {
            e.record(i as f64, i as u32, [0.0; 3], [0.0; 3], 0.0);
        }
        assert_eq!(e.row_count(), 5);
    }

    #[test]
    fn test_exporter_clear() {
        let mut e = PyDataExporter::new();
        e.record(0.0, 0, [0.0; 3], [0.0; 3], 0.0);
        e.clear();
        assert_eq!(e.row_count(), 0);
    }
}
