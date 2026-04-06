// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly analytics bridge.
//!
//! Provides types for tracking physics simulation analytics including energy,
//! performance, collision statistics, and constraint residuals. Designed for
//! serialisation across the WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// AnalyticsConfig
// ---------------------------------------------------------------------------

/// Configuration for the analytics subsystem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalyticsConfig {
    /// Maximum number of history entries to retain.
    pub max_history: usize,
    /// Number of simulation steps between samples.
    pub sample_interval: usize,
    /// Export format string (e.g. "csv", "json").
    pub export_format: String,
}

impl Default for AnalyticsConfig {
    fn default() -> Self {
        AnalyticsConfig {
            max_history: 1000,
            sample_interval: 1,
            export_format: "json".to_string(),
        }
    }
}

impl AnalyticsConfig {
    /// Create a new `AnalyticsConfig` with the given settings.
    pub fn new(max_history: usize, sample_interval: usize, export_format: &str) -> Self {
        AnalyticsConfig {
            max_history,
            sample_interval,
            export_format: export_format.to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// TimeSeriesData
// ---------------------------------------------------------------------------

/// A labelled time-series of (timestamp, value) pairs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeSeriesData {
    /// Timestamps (seconds).
    pub timestamps: Vec<f64>,
    /// Scalar values corresponding to each timestamp.
    pub values: Vec<f64>,
    /// Human-readable label for this series.
    pub label: String,
}

impl TimeSeriesData {
    /// Create a new empty `TimeSeriesData` with the given label.
    pub fn new(label: &str) -> Self {
        TimeSeriesData {
            timestamps: Vec::new(),
            values: Vec::new(),
            label: label.to_string(),
        }
    }

    /// Append a (timestamp, value) pair.
    pub fn push(&mut self, t: f64, v: f64) {
        self.timestamps.push(t);
        self.values.push(v);
    }

    /// Compute the arithmetic mean of the stored values.
    ///
    /// Returns `0.0` if empty.
    pub fn mean(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }

    /// Compute the population variance of the stored values.
    ///
    /// Returns `0.0` if empty.
    pub fn variance(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        let m = self.mean();
        self.values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / self.values.len() as f64
    }

    /// Return the minimum stored value, or `f64::NAN` if empty.
    pub fn min(&self) -> f64 {
        self.values.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Return the maximum stored value, or `f64::NEG_INFINITY` if empty.
    pub fn max(&self) -> f64 {
        self.values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Return a flat `Vec`f64` suitable for passing across the WASM boundary.
    ///
    /// Layout: `\[t0, v0, t1, v1, ...\]`
    pub fn to_js_array(&self) -> Vec<f64> {
        self.timestamps
            .iter()
            .zip(self.values.iter())
            .flat_map(|(t, v)| [*t, *v])
            .collect()
    }

    /// Number of data points stored.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns `true` if no data points have been recorded.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

// ---------------------------------------------------------------------------
// PhysicsMetrics
// ---------------------------------------------------------------------------

/// A snapshot of key physics quantities at a single simulation step.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PhysicsMetrics {
    /// Kinetic energy (J).
    pub kinetic_energy: f64,
    /// Potential energy (J).
    pub potential_energy: f64,
    /// Total mechanical energy (J).
    pub total_energy: f64,
    /// Linear momentum vector [px, py, pz] (kg·m/s).
    pub momentum: [f64; 3],
    /// Angular momentum vector [Lx, Ly, Lz] (kg·m²/s).
    pub angular_momentum: [f64; 3],
}

impl PhysicsMetrics {
    /// Create a new `PhysicsMetrics` snapshot.
    pub fn new(
        kinetic_energy: f64,
        potential_energy: f64,
        momentum: [f64; 3],
        angular_momentum: [f64; 3],
    ) -> Self {
        PhysicsMetrics {
            kinetic_energy,
            potential_energy,
            total_energy: kinetic_energy + potential_energy,
            momentum,
            angular_momentum,
        }
    }
}

// ---------------------------------------------------------------------------
// PerformanceTracker
// ---------------------------------------------------------------------------

/// Tracks per-frame timing to compute FPS and frame-time statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceTracker {
    /// Per-frame wall-clock durations in milliseconds.
    pub frame_times: Vec<f64>,
}

impl PerformanceTracker {
    /// Create a new empty `PerformanceTracker`.
    pub fn new() -> Self {
        PerformanceTracker {
            frame_times: Vec::new(),
        }
    }

    /// Record a new frame with duration `dt` milliseconds.
    pub fn push_frame(&mut self, dt: f64) {
        self.frame_times.push(dt);
    }

    /// Compute the mean frames-per-second over all recorded frames.
    ///
    /// Returns `0.0` if no frames have been recorded.
    pub fn fps(&self) -> f64 {
        if self.frame_times.is_empty() {
            return 0.0;
        }
        let avg_ms = self.avg_frame_ms();
        if avg_ms <= 0.0 {
            return 0.0;
        }
        1000.0 / avg_ms
    }

    /// Compute the mean frame duration in milliseconds.
    ///
    /// Returns `0.0` if no frames have been recorded.
    pub fn avg_frame_ms(&self) -> f64 {
        if self.frame_times.is_empty() {
            return 0.0;
        }
        self.frame_times.iter().sum::<f64>() / self.frame_times.len() as f64
    }

    /// Return the longest recorded frame duration in milliseconds.
    ///
    /// Returns `0.0` if no frames have been recorded.
    pub fn worst_frame_ms(&self) -> f64 {
        self.frame_times.iter().cloned().fold(0.0_f64, f64::max)
    }
}

impl Default for PerformanceTracker {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EnergyTracker
// ---------------------------------------------------------------------------

/// Tracks a history of `PhysicsMetrics` to monitor energy conservation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyTracker {
    /// Ordered list of recorded metrics snapshots.
    pub history: Vec<PhysicsMetrics>,
}

impl EnergyTracker {
    /// Create a new empty `EnergyTracker`.
    pub fn new() -> Self {
        EnergyTracker {
            history: Vec::new(),
        }
    }

    /// Append a metrics snapshot.
    pub fn push(&mut self, m: PhysicsMetrics) {
        self.history.push(m);
    }

    /// Compute energy drift: `(E_last - E_first) / E_first`.
    ///
    /// Returns `0.0` if fewer than two entries are recorded or if `E_first` is
    /// zero.
    pub fn energy_drift(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        let e0 = self.history[0].total_energy;
        if e0 == 0.0 {
            return 0.0;
        }
        let e_last = self
            .history
            .last()
            .expect("collection should not be empty")
            .total_energy;
        (e_last - e0) / e0.abs()
    }

    /// Return the maximum total energy seen across all recorded snapshots.
    ///
    /// Returns `0.0` if history is empty.
    pub fn max_energy(&self) -> f64 {
        self.history
            .iter()
            .map(|m| m.total_energy)
            .fold(0.0_f64, f64::max)
    }

    /// Return a `Vec<\[f64; 2\]>` of `\[index, total_energy\]` pairs for plotting.
    pub fn plot_data(&self) -> Vec<[f64; 2]> {
        self.history
            .iter()
            .enumerate()
            .map(|(i, m)| [i as f64, m.total_energy])
            .collect()
    }
}

impl Default for EnergyTracker {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CollisionStats
// ---------------------------------------------------------------------------

/// Aggregated collision-detection statistics for a simulation step.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct CollisionStats {
    /// Number of broad-phase collision pairs.
    pub n_broad: usize,
    /// Number of narrow-phase collision tests.
    pub n_narrow: usize,
    /// Number of confirmed contacts.
    pub n_contacts: usize,
}

impl CollisionStats {
    /// Create a new zeroed `CollisionStats`.
    pub fn new() -> Self {
        CollisionStats::default()
    }

    /// Update stats with counts from the latest simulation step.
    pub fn update(&mut self, broad: usize, narrow: usize, contacts: usize) {
        self.n_broad = broad;
        self.n_narrow = narrow;
        self.n_contacts = contacts;
    }

    /// Compute narrow-phase / broad-phase efficiency ratio.
    ///
    /// Returns `0.0` when no broad-phase pairs exist.
    pub fn efficiency_ratio(&self) -> f64 {
        if self.n_broad == 0 {
            return 0.0;
        }
        self.n_narrow as f64 / self.n_broad as f64
    }
}

// ---------------------------------------------------------------------------
// ConstraintResiduals
// ---------------------------------------------------------------------------

/// Tracks per-iteration solver residuals to monitor convergence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConstraintResiduals {
    /// Residual value recorded at each solver iteration.
    pub residuals: Vec<f64>,
}

impl ConstraintResiduals {
    /// Create a new empty `ConstraintResiduals` tracker.
    pub fn new() -> Self {
        ConstraintResiduals {
            residuals: Vec::new(),
        }
    }

    /// Record the residual at solver iteration `_iter`.
    pub fn push_residual(&mut self, _iter: usize, residual: f64) {
        self.residuals.push(residual);
    }

    /// Returns `true` if the last recorded residual is below `tol`.
    ///
    /// Returns `false` if no residuals have been recorded.
    pub fn converged(&self, tol: f64) -> bool {
        match self.residuals.last() {
            Some(&r) => r < tol,
            None => false,
        }
    }

    /// Return the full residual history as a `Vec`f64`.
    pub fn residual_history(&self) -> Vec<f64> {
        self.residuals.clone()
    }

    /// Clear all recorded residuals.
    pub fn clear(&mut self) {
        self.residuals.clear();
    }
}

impl Default for ConstraintResiduals {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// AnalyticsDashboard
// ---------------------------------------------------------------------------

/// Top-level analytics dashboard aggregating all tracker sub-systems.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalyticsDashboard {
    /// Configuration governing sampling and export behaviour.
    pub config: AnalyticsConfig,
    /// Energy conservation tracker.
    pub energy: EnergyTracker,
    /// Frame-time performance tracker.
    pub performance: PerformanceTracker,
    /// Collision detection statistics.
    pub collisions: CollisionStats,
    /// Constraint solver residuals.
    pub residuals: ConstraintResiduals,
    /// Kinetic energy time series.
    pub ke_series: TimeSeriesData,
    /// Potential energy time series.
    pub pe_series: TimeSeriesData,
    /// Current simulation time (seconds).
    pub sim_time: f64,
}

impl AnalyticsDashboard {
    /// Create a new `AnalyticsDashboard` with default configuration.
    pub fn new() -> Self {
        AnalyticsDashboard {
            config: AnalyticsConfig::default(),
            energy: EnergyTracker::new(),
            performance: PerformanceTracker::new(),
            collisions: CollisionStats::new(),
            residuals: ConstraintResiduals::new(),
            ke_series: TimeSeriesData::new("kinetic_energy"),
            pe_series: TimeSeriesData::new("potential_energy"),
            sim_time: 0.0,
        }
    }

    /// Record kinetic energy `ke` and potential energy `pe` at the current
    /// simulation time.
    pub fn update_energy(&mut self, ke: f64, pe: f64) {
        let metrics = PhysicsMetrics::new(ke, pe, [0.0; 3], [0.0; 3]);
        self.energy.push(metrics);
        self.ke_series.push(self.sim_time, ke);
        self.pe_series.push(self.sim_time, pe);
    }

    /// Record a frame with wall-clock duration `dt` milliseconds and advance
    /// simulation time accordingly.
    pub fn update_performance(&mut self, dt: f64) {
        self.performance.push_frame(dt);
        self.sim_time += dt / 1000.0;
    }

    /// Update collision statistics for the current step.
    pub fn update_collisions(&mut self, b: usize, n: usize, c: usize) {
        self.collisions.update(b, n, c);
    }

    /// Produce a JSON summary string of the current dashboard state.
    pub fn summary_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
}

impl Default for AnalyticsDashboard {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- AnalyticsConfig ----

    #[test]
    fn test_analytics_config_default() {
        let cfg = AnalyticsConfig::default();
        assert_eq!(cfg.max_history, 1000);
        assert_eq!(cfg.sample_interval, 1);
        assert_eq!(cfg.export_format, "json");
    }

    #[test]
    fn test_analytics_config_new() {
        let cfg = AnalyticsConfig::new(500, 5, "csv");
        assert_eq!(cfg.max_history, 500);
        assert_eq!(cfg.sample_interval, 5);
        assert_eq!(cfg.export_format, "csv");
    }

    // ---- TimeSeriesData ----

    #[test]
    fn test_time_series_empty() {
        let ts = TimeSeriesData::new("test");
        assert!(ts.is_empty());
        assert_eq!(ts.len(), 0);
        assert_eq!(ts.mean(), 0.0);
        assert_eq!(ts.variance(), 0.0);
    }

    #[test]
    fn test_time_series_push_len() {
        let mut ts = TimeSeriesData::new("test");
        ts.push(0.0, 1.0);
        ts.push(1.0, 2.0);
        assert_eq!(ts.len(), 2);
        assert!(!ts.is_empty());
    }

    #[test]
    fn test_time_series_mean() {
        let mut ts = TimeSeriesData::new("x");
        ts.push(0.0, 2.0);
        ts.push(1.0, 4.0);
        ts.push(2.0, 6.0);
        assert!((ts.mean() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_variance() {
        let mut ts = TimeSeriesData::new("x");
        ts.push(0.0, 2.0);
        ts.push(1.0, 4.0);
        ts.push(2.0, 6.0);
        // mean=4, var = ((2-4)^2 + (4-4)^2 + (6-4)^2) / 3 = 8/3
        let expected = 8.0 / 3.0;
        assert!((ts.variance() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_min_max() {
        let mut ts = TimeSeriesData::new("x");
        ts.push(0.0, 3.0);
        ts.push(1.0, 1.0);
        ts.push(2.0, 7.0);
        assert!((ts.min() - 1.0).abs() < 1e-10);
        assert!((ts.max() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_to_js_array() {
        let mut ts = TimeSeriesData::new("x");
        ts.push(0.5, 1.5);
        ts.push(1.5, 2.5);
        let arr = ts.to_js_array();
        assert_eq!(arr.len(), 4);
        assert!((arr[0] - 0.5).abs() < 1e-10);
        assert!((arr[1] - 1.5).abs() < 1e-10);
        assert!((arr[2] - 1.5).abs() < 1e-10);
        assert!((arr[3] - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_time_series_single_value_variance() {
        let mut ts = TimeSeriesData::new("x");
        ts.push(0.0, 5.0);
        assert!((ts.variance() - 0.0).abs() < 1e-10);
    }

    // ---- PhysicsMetrics ----

    #[test]
    fn test_physics_metrics_total_energy() {
        let m = PhysicsMetrics::new(10.0, 5.0, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert!((m.total_energy - 15.0).abs() < 1e-10);
        assert!((m.kinetic_energy - 10.0).abs() < 1e-10);
        assert!((m.potential_energy - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_physics_metrics_momentum() {
        let m = PhysicsMetrics::new(1.0, 1.0, [3.0, 4.0, 5.0], [0.0; 3]);
        assert!((m.momentum[0] - 3.0).abs() < 1e-10);
        assert!((m.momentum[1] - 4.0).abs() < 1e-10);
        assert!((m.momentum[2] - 5.0).abs() < 1e-10);
    }

    // ---- PerformanceTracker ----

    #[test]
    fn test_performance_tracker_empty() {
        let pt = PerformanceTracker::new();
        assert_eq!(pt.fps(), 0.0);
        assert_eq!(pt.avg_frame_ms(), 0.0);
        assert_eq!(pt.worst_frame_ms(), 0.0);
    }

    #[test]
    fn test_performance_tracker_fps() {
        let mut pt = PerformanceTracker::new();
        pt.push_frame(16.666);
        pt.push_frame(16.666);
        let fps = pt.fps();
        assert!(fps > 59.0 && fps < 61.0);
    }

    #[test]
    fn test_performance_tracker_avg_frame_ms() {
        let mut pt = PerformanceTracker::new();
        pt.push_frame(10.0);
        pt.push_frame(20.0);
        pt.push_frame(30.0);
        assert!((pt.avg_frame_ms() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_performance_tracker_worst_frame() {
        let mut pt = PerformanceTracker::new();
        pt.push_frame(5.0);
        pt.push_frame(50.0);
        pt.push_frame(10.0);
        assert!((pt.worst_frame_ms() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_performance_tracker_default() {
        let pt = PerformanceTracker::default();
        assert!(pt.frame_times.is_empty());
    }

    // ---- EnergyTracker ----

    #[test]
    fn test_energy_tracker_empty() {
        let et = EnergyTracker::new();
        assert_eq!(et.energy_drift(), 0.0);
        assert_eq!(et.max_energy(), 0.0);
        assert!(et.plot_data().is_empty());
    }

    #[test]
    fn test_energy_tracker_push() {
        let mut et = EnergyTracker::new();
        et.push(PhysicsMetrics::new(10.0, 5.0, [0.0; 3], [0.0; 3]));
        assert_eq!(et.history.len(), 1);
    }

    #[test]
    fn test_energy_tracker_drift() {
        let mut et = EnergyTracker::new();
        et.push(PhysicsMetrics::new(10.0, 0.0, [0.0; 3], [0.0; 3]));
        et.push(PhysicsMetrics::new(11.0, 0.0, [0.0; 3], [0.0; 3]));
        // drift = (11-10)/10 = 0.1
        assert!((et.energy_drift() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_energy_tracker_drift_zero_initial() {
        let mut et = EnergyTracker::new();
        et.push(PhysicsMetrics::new(0.0, 0.0, [0.0; 3], [0.0; 3]));
        et.push(PhysicsMetrics::new(1.0, 0.0, [0.0; 3], [0.0; 3]));
        assert_eq!(et.energy_drift(), 0.0);
    }

    #[test]
    fn test_energy_tracker_max_energy() {
        let mut et = EnergyTracker::new();
        et.push(PhysicsMetrics::new(5.0, 0.0, [0.0; 3], [0.0; 3]));
        et.push(PhysicsMetrics::new(20.0, 0.0, [0.0; 3], [0.0; 3]));
        et.push(PhysicsMetrics::new(10.0, 0.0, [0.0; 3], [0.0; 3]));
        assert!((et.max_energy() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_energy_tracker_plot_data() {
        let mut et = EnergyTracker::new();
        et.push(PhysicsMetrics::new(3.0, 2.0, [0.0; 3], [0.0; 3]));
        et.push(PhysicsMetrics::new(4.0, 1.0, [0.0; 3], [0.0; 3]));
        let pd = et.plot_data();
        assert_eq!(pd.len(), 2);
        assert!((pd[0][0] - 0.0).abs() < 1e-10);
        assert!((pd[0][1] - 5.0).abs() < 1e-10);
        assert!((pd[1][0] - 1.0).abs() < 1e-10);
        assert!((pd[1][1] - 5.0).abs() < 1e-10);
    }

    // ---- CollisionStats ----

    #[test]
    fn test_collision_stats_default() {
        let cs = CollisionStats::new();
        assert_eq!(cs.n_broad, 0);
        assert_eq!(cs.n_narrow, 0);
        assert_eq!(cs.n_contacts, 0);
    }

    #[test]
    fn test_collision_stats_update() {
        let mut cs = CollisionStats::new();
        cs.update(100, 40, 10);
        assert_eq!(cs.n_broad, 100);
        assert_eq!(cs.n_narrow, 40);
        assert_eq!(cs.n_contacts, 10);
    }

    #[test]
    fn test_collision_stats_efficiency_ratio() {
        let mut cs = CollisionStats::new();
        cs.update(200, 50, 10);
        assert!((cs.efficiency_ratio() - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_collision_stats_efficiency_ratio_zero_broad() {
        let cs = CollisionStats::new();
        assert_eq!(cs.efficiency_ratio(), 0.0);
    }

    // ---- ConstraintResiduals ----

    #[test]
    fn test_constraint_residuals_empty() {
        let cr = ConstraintResiduals::new();
        assert!(!cr.converged(1e-6));
        assert!(cr.residual_history().is_empty());
    }

    #[test]
    fn test_constraint_residuals_push() {
        let mut cr = ConstraintResiduals::new();
        cr.push_residual(0, 1.0);
        cr.push_residual(1, 0.5);
        cr.push_residual(2, 0.01);
        assert_eq!(cr.residuals.len(), 3);
    }

    #[test]
    fn test_constraint_residuals_converged_true() {
        let mut cr = ConstraintResiduals::new();
        cr.push_residual(0, 1.0);
        cr.push_residual(1, 1e-8);
        assert!(cr.converged(1e-6));
    }

    #[test]
    fn test_constraint_residuals_converged_false() {
        let mut cr = ConstraintResiduals::new();
        cr.push_residual(0, 0.5);
        assert!(!cr.converged(1e-6));
    }

    #[test]
    fn test_constraint_residuals_clear() {
        let mut cr = ConstraintResiduals::new();
        cr.push_residual(0, 1.0);
        cr.clear();
        assert!(cr.residuals.is_empty());
    }

    #[test]
    fn test_constraint_residuals_history_clone() {
        let mut cr = ConstraintResiduals::new();
        cr.push_residual(0, 2.0);
        cr.push_residual(1, 1.0);
        let hist = cr.residual_history();
        assert_eq!(hist.len(), 2);
        assert!((hist[0] - 2.0).abs() < 1e-10);
        assert!((hist[1] - 1.0).abs() < 1e-10);
    }

    // ---- AnalyticsDashboard ----

    #[test]
    fn test_dashboard_new() {
        let db = AnalyticsDashboard::new();
        assert!((db.sim_time - 0.0).abs() < 1e-10);
        assert!(db.energy.history.is_empty());
    }

    #[test]
    fn test_dashboard_update_energy() {
        let mut db = AnalyticsDashboard::new();
        db.update_energy(10.0, 5.0);
        assert_eq!(db.energy.history.len(), 1);
        assert_eq!(db.ke_series.len(), 1);
        assert_eq!(db.pe_series.len(), 1);
    }

    #[test]
    fn test_dashboard_update_performance() {
        let mut db = AnalyticsDashboard::new();
        db.update_performance(16.666);
        assert_eq!(db.performance.frame_times.len(), 1);
        assert!(db.sim_time > 0.0);
    }

    #[test]
    fn test_dashboard_update_collisions() {
        let mut db = AnalyticsDashboard::new();
        db.update_collisions(50, 20, 5);
        assert_eq!(db.collisions.n_broad, 50);
    }

    #[test]
    fn test_dashboard_summary_json_valid() {
        let mut db = AnalyticsDashboard::new();
        db.update_energy(5.0, 2.0);
        let json = db.summary_json();
        assert!(!json.is_empty());
        assert!(json.contains("kinetic_energy") || json.contains("ke_series"));
    }

    #[test]
    fn test_dashboard_default() {
        let db = AnalyticsDashboard::default();
        assert_eq!(db.config.max_history, 1000);
    }

    #[test]
    fn test_dashboard_multiple_updates() {
        let mut db = AnalyticsDashboard::new();
        for i in 0..10 {
            db.update_energy(i as f64, (10 - i) as f64);
            db.update_performance(16.0);
        }
        assert_eq!(db.energy.history.len(), 10);
        assert_eq!(db.performance.frame_times.len(), 10);
    }

    #[test]
    fn test_time_series_min_empty() {
        let ts = TimeSeriesData::new("empty");
        // min on empty returns INFINITY
        assert!(ts.min().is_infinite());
    }

    #[test]
    fn test_time_series_max_empty() {
        let ts = TimeSeriesData::new("empty");
        // max on empty returns NEG_INFINITY
        assert!(ts.max().is_infinite());
    }

    #[test]
    fn test_performance_tracker_single_frame() {
        let mut pt = PerformanceTracker::new();
        pt.push_frame(20.0);
        assert!((pt.avg_frame_ms() - 20.0).abs() < 1e-10);
        assert!((pt.fps() - 50.0).abs() < 1e-10);
    }
}
