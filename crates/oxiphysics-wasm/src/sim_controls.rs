// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly simulation control APIs.
//!
//! This module provides pure-Rust structs for controlling physics simulations,
//! recording replay data, profiling performance, and exporting simulation state.
//! All types are designed for ergonomic use from JavaScript via wasm-bindgen.

#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// BroadphaseType
// ---------------------------------------------------------------------------

/// Broadphase algorithm selection for collision detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum BroadphaseType {
    /// Sweep and prune along a single axis.
    SweepAndPrune,
    /// Bounding volume hierarchy tree.
    #[default]
    Bvh,
    /// Simple O(n²) brute-force (useful for small scenes).
    BruteForce,
    /// Spatial hashing grid.
    SpatialHash,
}

// ---------------------------------------------------------------------------
// SimulationConfig
// ---------------------------------------------------------------------------

/// Configuration parameters for the simulation.
///
/// Controls the integration timestep, gravity vector, solver iterations,
/// sleeping threshold, and broadphase algorithm.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    /// Fixed timestep in seconds (e.g. 1/60).
    pub timestep: f64,
    /// Gravity vector \[x, y, z\] in m/s².
    pub gravity: [f64; 3],
    /// Number of constraint solver iterations per step.
    pub iterations: u32,
    /// Linear velocity threshold below which a body may sleep (m/s).
    pub sleeping_threshold_linear: f64,
    /// Angular velocity threshold below which a body may sleep (rad/s).
    pub sleeping_threshold_angular: f64,
    /// Broadphase algorithm to use.
    pub broadphase: BroadphaseType,
    /// Maximum number of sub-steps per call to `step`.
    pub max_substeps: u32,
    /// Whether to allow bodies to sleep.
    pub sleeping_enabled: bool,
    /// CCD (continuous collision detection) enabled.
    pub ccd_enabled: bool,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        SimulationConfig {
            timestep: 1.0 / 60.0,
            gravity: [0.0, -9.81, 0.0],
            iterations: 10,
            sleeping_threshold_linear: 0.01,
            sleeping_threshold_angular: 0.01,
            broadphase: BroadphaseType::Bvh,
            max_substeps: 4,
            sleeping_enabled: true,
            ccd_enabled: false,
        }
    }
}

impl SimulationConfig {
    /// Create a new `SimulationConfig` with the given timestep and gravity.
    pub fn new(timestep: f64, gravity: [f64; 3]) -> Self {
        SimulationConfig {
            timestep,
            gravity,
            ..Default::default()
        }
    }

    /// Set the number of solver iterations.
    pub fn with_iterations(mut self, iterations: u32) -> Self {
        self.iterations = iterations;
        self
    }

    /// Set sleeping thresholds (linear and angular).
    pub fn with_sleeping(mut self, linear: f64, angular: f64) -> Self {
        self.sleeping_threshold_linear = linear;
        self.sleeping_threshold_angular = angular;
        self
    }

    /// Set the broadphase algorithm.
    pub fn with_broadphase(mut self, broadphase: BroadphaseType) -> Self {
        self.broadphase = broadphase;
        self
    }

    /// Enable or disable CCD.
    pub fn with_ccd(mut self, enabled: bool) -> Self {
        self.ccd_enabled = enabled;
        self
    }

    /// Validate configuration values.
    pub fn validate(&self) -> Result<(), String> {
        if self.timestep <= 0.0 {
            return Err("timestep must be positive".to_string());
        }
        if self.iterations == 0 {
            return Err("iterations must be at least 1".to_string());
        }
        if self.max_substeps == 0 {
            return Err("max_substeps must be at least 1".to_string());
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// SimulationState
// ---------------------------------------------------------------------------

/// The current execution state of the simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SimulationState {
    /// Simulation is actively advancing every frame.
    Running,
    /// Simulation is paused; no integration is performed.
    #[default]
    Paused,
    /// Simulation advances exactly one timestep and then pauses.
    Stepping,
    /// Simulation is playing back recorded frames in reverse.
    Rewinding,
    /// Simulation is running and every frame is being recorded.
    Recording,
}

impl SimulationState {
    /// Returns `true` if integration should occur this frame.
    pub fn is_active(&self) -> bool {
        matches!(
            self,
            SimulationState::Running | SimulationState::Stepping | SimulationState::Recording
        )
    }

    /// Returns `true` if the simulation is recording frames.
    pub fn is_recording(&self) -> bool {
        *self == SimulationState::Recording
    }

    /// Human-readable label for the state.
    pub fn label(&self) -> &'static str {
        match self {
            SimulationState::Running => "running",
            SimulationState::Paused => "paused",
            SimulationState::Stepping => "stepping",
            SimulationState::Rewinding => "rewinding",
            SimulationState::Recording => "recording",
        }
    }
}

// ---------------------------------------------------------------------------
// SimulationStats
// ---------------------------------------------------------------------------

/// Aggregate statistics collected over the lifetime of the simulation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SimulationStats {
    /// Total number of integration steps completed.
    pub step_count: u64,
    /// Total simulated time in seconds.
    pub elapsed_time: f64,
    /// Current number of rigid bodies in the scene.
    pub body_count: u32,
    /// Current number of active constraints.
    pub constraint_count: u32,
    /// Number of simulation islands (disconnected groups) this step.
    pub island_count: u32,
    /// Number of active collision pairs detected this step.
    pub collision_pairs: u32,
    /// Peak memory used by the physics system (bytes).
    pub peak_memory_bytes: u64,
    /// Number of bodies currently sleeping.
    pub sleeping_bodies: u32,
    /// Number of CCD sub-steps taken this step.
    pub ccd_substeps: u32,
}

impl SimulationStats {
    /// Create a zeroed `SimulationStats`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Advance step count and elapsed time.
    pub fn advance(&mut self, dt: f64) {
        self.step_count += 1;
        self.elapsed_time += dt;
    }

    /// Reset all counters while keeping lifetime step_count/elapsed_time.
    pub fn reset_frame_counters(&mut self) {
        self.collision_pairs = 0;
        self.island_count = 0;
        self.ccd_substeps = 0;
    }
}

// ---------------------------------------------------------------------------
// StepResult
// ---------------------------------------------------------------------------

/// Per-step diagnostic result returned by `SimulationController::step`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StepResult {
    /// Number of new contact manifolds generated this step.
    pub new_contacts: u32,
    /// Number of contact manifolds resolved (separated) this step.
    pub resolved_contacts: u32,
    /// Number of bodies that transitioned to sleep.
    pub sleeping_bodies: u32,
    /// Number of active simulation islands this step.
    pub active_islands: u32,
    /// Wall-clock time spent in physics this step (ms).
    pub perf_ms: f64,
    /// Whether any constraint was violated beyond the tolerance.
    pub constraint_violation: bool,
}

impl StepResult {
    /// Create a zeroed `StepResult`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Return `true` if any performance-relevant event occurred.
    pub fn has_events(&self) -> bool {
        self.new_contacts > 0 || self.sleeping_bodies > 0 || self.constraint_violation
    }
}

// ---------------------------------------------------------------------------
// SimulationController
// ---------------------------------------------------------------------------

/// High-level controller that owns a `SimulationConfig` and drives the loop.
///
/// `SimulationController` does not depend on wasm-bindgen and can be used
/// in native unit tests.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationController {
    /// Active configuration.
    pub config: SimulationConfig,
    /// Current execution state.
    pub state: SimulationState,
    /// Accumulated statistics.
    pub stats: SimulationStats,
    /// Accumulated time not yet consumed by fixed steps.
    pub accumulator: f64,
}

impl SimulationController {
    /// Create a new controller with the given config. Initial state is `Paused`.
    pub fn new(config: SimulationConfig) -> Self {
        SimulationController {
            config,
            state: SimulationState::Paused,
            stats: SimulationStats::new(),
            accumulator: 0.0,
        }
    }

    /// Advance the simulation by `frame_time` seconds (wall-clock delta).
    ///
    /// Returns a `StepResult` summarising this frame's activity.
    pub fn step(&mut self, frame_time: f64) -> StepResult {
        let mut result = StepResult::new();

        if !self.state.is_active() {
            return result;
        }

        self.accumulator += frame_time;
        let dt = self.config.timestep;
        let mut substeps = 0u32;

        while self.accumulator >= dt && substeps < self.config.max_substeps {
            self.stats.advance(dt);
            self.accumulator -= dt;
            substeps += 1;
            result.active_islands = self.stats.island_count;
        }

        // If we were in Stepping mode, pause after one logical frame.
        if self.state == SimulationState::Stepping {
            self.state = SimulationState::Paused;
        }

        result.perf_ms = frame_time * 1_000.0 * 0.01; // mock
        result
    }

    /// Pause the simulation.
    pub fn pause(&mut self) {
        self.state = SimulationState::Paused;
    }

    /// Resume from paused.
    pub fn resume(&mut self) {
        self.state = SimulationState::Running;
    }

    /// Advance exactly one fixed-timestep and pause.
    pub fn step_once(&mut self) {
        self.state = SimulationState::Stepping;
        self.accumulator += self.config.timestep;
    }

    /// Begin recording frames.
    pub fn start_recording(&mut self) {
        self.state = SimulationState::Recording;
    }

    /// Begin rewinding recorded frames.
    pub fn start_rewind(&mut self) {
        self.state = SimulationState::Rewinding;
    }

    /// Reset simulation to t=0.
    pub fn reset(&mut self) {
        self.stats = SimulationStats::new();
        self.accumulator = 0.0;
        self.state = SimulationState::Paused;
    }

    /// Get a reference to the current statistics.
    pub fn get_stats(&self) -> &SimulationStats {
        &self.stats
    }

    /// Update body and constraint counts in stats.
    pub fn update_counts(&mut self, bodies: u32, constraints: u32) {
        self.stats.body_count = bodies;
        self.stats.constraint_count = constraints;
    }
}

// ---------------------------------------------------------------------------
// ReplayBuffer
// ---------------------------------------------------------------------------

/// A single recorded physics frame.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayFrame {
    /// Simulation time at this frame (s).
    pub time: f64,
    /// Flat array of body positions: \[x0,y0,z0, x1,y1,z1, ...\].
    pub positions: Vec<f64>,
    /// Flat array of body linear velocities.
    pub velocities: Vec<f64>,
    /// Flat array of body quaternions \[qx,qy,qz,qw, ...\].
    pub orientations: Vec<f64>,
    /// Frame index.
    pub frame_index: u64,
}

impl ReplayFrame {
    /// Create a new empty frame at the given time.
    pub fn new(time: f64, frame_index: u64) -> Self {
        ReplayFrame {
            time,
            positions: Vec::new(),
            velocities: Vec::new(),
            orientations: Vec::new(),
            frame_index,
        }
    }

    /// Number of bodies encoded in this frame.
    pub fn body_count(&self) -> usize {
        self.positions.len() / 3
    }
}

/// Ring-buffer for recording and playing back simulation frames.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayBuffer {
    /// All recorded frames in chronological order.
    pub frames: Vec<ReplayFrame>,
    /// Maximum number of frames to retain.
    pub capacity: usize,
    /// Current playback cursor index.
    pub cursor: usize,
    /// Whether the buffer is currently in record mode.
    pub recording: bool,
}

impl ReplayBuffer {
    /// Create a new `ReplayBuffer` with the given capacity.
    pub fn new(capacity: usize) -> Self {
        ReplayBuffer {
            frames: Vec::with_capacity(capacity),
            capacity,
            cursor: 0,
            recording: false,
        }
    }

    /// Record a new frame. If at capacity the oldest frame is discarded.
    pub fn record_state(&mut self, frame: ReplayFrame) {
        if self.frames.len() >= self.capacity {
            self.frames.remove(0);
        }
        self.frames.push(frame);
    }

    /// Return the frame at the current playback cursor, or `None`.
    pub fn playback_frame(&self) -> Option<&ReplayFrame> {
        self.frames.get(self.cursor)
    }

    /// Seek playback cursor to the frame nearest to `time`.
    pub fn seek_to_time(&mut self, time: f64) -> usize {
        if self.frames.is_empty() {
            self.cursor = 0;
            return 0;
        }
        let mut best = 0;
        let mut best_diff = f64::MAX;
        for (i, f) in self.frames.iter().enumerate() {
            let diff = (f.time - time).abs();
            if diff < best_diff {
                best_diff = diff;
                best = i;
            }
        }
        self.cursor = best;
        best
    }

    /// Advance playback cursor by one frame; returns `false` if at end.
    pub fn advance_cursor(&mut self) -> bool {
        if self.cursor + 1 < self.frames.len() {
            self.cursor += 1;
            true
        } else {
            false
        }
    }

    /// Rewind playback cursor by one frame; returns `false` if at start.
    pub fn rewind_cursor(&mut self) -> bool {
        if self.cursor > 0 {
            self.cursor -= 1;
            true
        } else {
            false
        }
    }

    /// Export all frames to a JSON string.
    pub fn export_frames(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }

    /// Import frames from a JSON string, replacing current buffer.
    pub fn import_frames(&mut self, json: &str) -> Result<(), String> {
        match serde_json::from_str::<ReplayBuffer>(json) {
            Ok(buf) => {
                *self = buf;
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    /// Clear all frames.
    pub fn clear(&mut self) {
        self.frames.clear();
        self.cursor = 0;
    }

    /// Total number of recorded frames.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Returns `true` if no frames have been recorded.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

// ---------------------------------------------------------------------------
// PhysicsProfiler
// ---------------------------------------------------------------------------

/// Per-step timing breakdown for the physics pipeline.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PhysicsProfilerSample {
    /// Time spent in broadphase collision detection (ms).
    pub broadphase_ms: f64,
    /// Time spent in narrowphase collision detection (ms).
    pub narrowphase_ms: f64,
    /// Time spent in constraint solver (ms).
    pub solver_ms: f64,
    /// Time spent in integration (ms).
    pub integration_ms: f64,
    /// Total physics time this step (ms).
    pub total_ms: f64,
    /// Step index at which this sample was taken.
    pub step_index: u64,
}

impl PhysicsProfilerSample {
    /// Compute total from sub-phases.
    pub fn compute_total(&mut self) {
        self.total_ms =
            self.broadphase_ms + self.narrowphase_ms + self.solver_ms + self.integration_ms;
    }
}

/// Profiler that accumulates timing history over multiple steps.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicsProfiler {
    /// Ring buffer of per-step samples.
    pub history: Vec<PhysicsProfilerSample>,
    /// Capacity of the history ring buffer.
    pub history_capacity: usize,
    /// Running average total_ms.
    pub avg_total_ms: f64,
    /// Peak total_ms seen.
    pub peak_total_ms: f64,
}

impl PhysicsProfiler {
    /// Create a new profiler with a history ring buffer of the given size.
    pub fn new(history_capacity: usize) -> Self {
        PhysicsProfiler {
            history: Vec::with_capacity(history_capacity),
            history_capacity,
            avg_total_ms: 0.0,
            peak_total_ms: 0.0,
        }
    }

    /// Record a new sample.
    pub fn record(&mut self, mut sample: PhysicsProfilerSample) {
        sample.compute_total();
        if sample.total_ms > self.peak_total_ms {
            self.peak_total_ms = sample.total_ms;
        }
        // Update running average.
        let n = (self.history.len() + 1) as f64;
        self.avg_total_ms = (self.avg_total_ms * (n - 1.0) + sample.total_ms) / n;
        if self.history.len() >= self.history_capacity {
            self.history.remove(0);
        }
        self.history.push(sample);
    }

    /// Return the most recent sample, if any.
    pub fn latest(&self) -> Option<&PhysicsProfilerSample> {
        self.history.last()
    }

    /// Reset all history.
    pub fn reset(&mut self) {
        self.history.clear();
        self.avg_total_ms = 0.0;
        self.peak_total_ms = 0.0;
    }

    /// Return the number of recorded samples.
    pub fn sample_count(&self) -> usize {
        self.history.len()
    }
}

// ---------------------------------------------------------------------------
// SimulationExporter
// ---------------------------------------------------------------------------

/// Export format selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExportFormat {
    /// JSON text format.
    Json,
    /// Binary (CBOR-like) format.
    Binary,
    /// VTK legacy ASCII format for ParaView.
    Vtk,
    /// CSV comma-separated values.
    Csv,
}

/// Serialised snapshot of one simulation frame for export.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    /// Simulation time (s).
    pub time: f64,
    /// Flat body positions \[x,y,z, ...\].
    pub positions: Vec<f64>,
    /// Flat body velocities \[vx,vy,vz, ...\].
    pub velocities: Vec<f64>,
    /// Flat body quaternions \[qx,qy,qz,qw, ...\].
    pub orientations: Vec<f64>,
    /// Flat angular velocities \[wx,wy,wz, ...\].
    pub angular_velocities: Vec<f64>,
    /// Body masses.
    pub masses: Vec<f64>,
}

impl SimulationSnapshot {
    /// Create an empty snapshot at the given time.
    pub fn new(time: f64) -> Self {
        SimulationSnapshot {
            time,
            positions: Vec::new(),
            velocities: Vec::new(),
            orientations: Vec::new(),
            angular_velocities: Vec::new(),
            masses: Vec::new(),
        }
    }

    /// Number of bodies in the snapshot.
    pub fn body_count(&self) -> usize {
        self.masses.len()
    }
}

/// Utility for exporting and importing simulation state.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SimulationExporter {
    /// Cached snapshots (e.g., for VTK series).
    pub snapshots: Vec<SimulationSnapshot>,
}

impl SimulationExporter {
    /// Create a new empty exporter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Export a single snapshot to JSON.
    pub fn export_json(&self, snapshot: &SimulationSnapshot) -> String {
        serde_json::to_string_pretty(snapshot).unwrap_or_default()
    }

    /// Export a single snapshot to a simple binary format (little-endian f64 stream).
    pub fn export_binary(&self, snapshot: &SimulationSnapshot) -> Vec<u8> {
        let mut buf = Vec::new();
        // Header: time as 8 bytes
        buf.extend_from_slice(&snapshot.time.to_le_bytes());
        let n = snapshot.positions.len() as u64;
        buf.extend_from_slice(&n.to_le_bytes());
        for &v in &snapshot.positions {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &snapshot.velocities {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// Export a snapshot in VTK legacy ASCII format.
    pub fn export_vtk(&self, snapshot: &SimulationSnapshot, frame_idx: usize) -> String {
        let n = snapshot.positions.len() / 3;
        let mut out = format!(
            "# vtk DataFile Version 3.0\nFrame {frame_idx} t={:.6}\nASCII\nDATASET UNSTRUCTURED_GRID\n",
            snapshot.time
        );
        out.push_str(&format!("POINTS {n} double\n"));
        for i in 0..n {
            let x = snapshot.positions.get(i * 3).copied().unwrap_or(0.0);
            let y = snapshot.positions.get(i * 3 + 1).copied().unwrap_or(0.0);
            let z = snapshot.positions.get(i * 3 + 2).copied().unwrap_or(0.0);
            out.push_str(&format!("{x:.6} {y:.6} {z:.6}\n"));
        }
        out
    }

    /// Export all cached snapshots as a VTK series (one string per frame).
    pub fn export_vtk_series(&self) -> Vec<String> {
        self.snapshots
            .iter()
            .enumerate()
            .map(|(i, s)| self.export_vtk(s, i))
            .collect()
    }

    /// Import a snapshot from JSON.
    pub fn import_state(&self, json: &str) -> Result<SimulationSnapshot, String> {
        serde_json::from_str(json).map_err(|e| e.to_string())
    }

    /// Add a snapshot to the internal series.
    pub fn push_snapshot(&mut self, snapshot: SimulationSnapshot) {
        self.snapshots.push(snapshot);
    }

    /// Clear cached snapshots.
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }
}

// ---------------------------------------------------------------------------
// TimeIntegrator
// ---------------------------------------------------------------------------

/// Integration scheme for advancing body state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum IntegratorType {
    /// Forward (explicit) Euler.
    Euler,
    /// Semi-implicit Euler (symplectic Euler).
    #[default]
    SemiImplicit,
    /// Velocity Verlet.
    Verlet,
    /// Fourth-order Runge-Kutta.
    Rk4,
}

impl IntegratorType {
    /// Human-readable name.
    pub fn name(&self) -> &'static str {
        match self {
            IntegratorType::Euler => "euler",
            IntegratorType::SemiImplicit => "semi_implicit",
            IntegratorType::Verlet => "verlet",
            IntegratorType::Rk4 => "rk4",
        }
    }

    /// Whether the integrator is symplectic (energy-conserving in limit).
    pub fn is_symplectic(&self) -> bool {
        matches!(self, IntegratorType::SemiImplicit | IntegratorType::Verlet)
    }
}

/// Configuration for the time integrator.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeIntegrator {
    /// Selected integration scheme.
    pub integrator_type: IntegratorType,
    /// Number of sub-steps per timestep.
    pub sub_steps: u32,
    /// Damping coefficient applied during integration.
    pub damping: f64,
    /// Angular damping coefficient.
    pub angular_damping: f64,
}

impl Default for TimeIntegrator {
    fn default() -> Self {
        TimeIntegrator {
            integrator_type: IntegratorType::SemiImplicit,
            sub_steps: 1,
            damping: 0.0,
            angular_damping: 0.0,
        }
    }
}

impl TimeIntegrator {
    /// Create with a specific scheme.
    pub fn new(integrator_type: IntegratorType) -> Self {
        TimeIntegrator {
            integrator_type,
            ..Default::default()
        }
    }

    /// Set sub-step count.
    pub fn with_substeps(mut self, n: u32) -> Self {
        self.sub_steps = n.max(1);
        self
    }

    /// Set linear damping.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping.clamp(0.0, 1.0);
        self
    }

    /// Integrate a 1D position/velocity pair using the configured scheme.
    ///
    /// Returns `(new_position, new_velocity)`.
    pub fn integrate_1d(&self, pos: f64, vel: f64, acc: f64, dt: f64) -> (f64, f64) {
        let sub_dt = dt / self.sub_steps as f64;
        let mut p = pos;
        let mut v = vel;
        for _ in 0..self.sub_steps {
            (p, v) = match self.integrator_type {
                IntegratorType::Euler => (p + v * sub_dt, v + acc * sub_dt),
                IntegratorType::SemiImplicit => {
                    let v2 = v + acc * sub_dt;
                    (p + v2 * sub_dt, v2)
                }
                IntegratorType::Verlet => {
                    let p2 = p + v * sub_dt + 0.5 * acc * sub_dt * sub_dt;
                    let v2 = v + acc * sub_dt;
                    (p2, v2)
                }
                IntegratorType::Rk4 => {
                    // Simplified RK4 for constant acceleration.
                    let k1v = acc;
                    let k2v = acc;
                    let k3v = acc;
                    let k4v = acc;
                    let k1p = v;
                    let k2p = v + 0.5 * sub_dt * k1v;
                    let k3p = v + 0.5 * sub_dt * k2v;
                    let k4p = v + sub_dt * k3v;
                    let v2 = v + (sub_dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
                    let p2 = p + (sub_dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
                    (p2, v2)
                }
            };
            v *= 1.0 - self.damping * sub_dt;
        }
        (p, v)
    }
}

// ---------------------------------------------------------------------------
// DebugDrawConfig
// ---------------------------------------------------------------------------

/// A packed RGBA colour as \[r, g, b, a\] each in \[0, 1\].
pub type Color4 = [f32; 4];

/// Configuration for debug visualisation overlays.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DebugDrawConfig {
    /// Draw axis-aligned bounding boxes.
    pub show_aabbs: bool,
    /// Draw contact points and normals.
    pub show_contacts: bool,
    /// Draw joint/constraint anchors and axes.
    pub show_joints: bool,
    /// Draw velocity arrows on bodies.
    pub show_velocities: bool,
    /// Draw force arrows on bodies.
    pub show_forces: bool,
    /// Colour used for contact visualisation.
    pub contact_color: Color4,
    /// Colour used for AABB visualisation.
    pub aabb_color: Color4,
    /// Colour used for joint visualisation.
    pub joint_color: Color4,
    /// Colour used for velocity arrows.
    pub velocity_color: Color4,
    /// Scale factor for arrow lengths.
    pub arrow_scale: f32,
    /// Whether to show island IDs as text.
    pub show_island_ids: bool,
    /// Whether to show body sleep state.
    pub show_sleep_state: bool,
}

impl Default for DebugDrawConfig {
    fn default() -> Self {
        DebugDrawConfig {
            show_aabbs: false,
            show_contacts: true,
            show_joints: true,
            show_velocities: false,
            show_forces: false,
            contact_color: [1.0, 0.0, 0.0, 1.0],
            aabb_color: [0.0, 1.0, 0.0, 0.5],
            joint_color: [0.0, 0.5, 1.0, 1.0],
            velocity_color: [1.0, 1.0, 0.0, 1.0],
            arrow_scale: 0.1,
            show_island_ids: false,
            show_sleep_state: false,
        }
    }
}

impl DebugDrawConfig {
    /// Create a config with all overlays disabled.
    pub fn disabled() -> Self {
        DebugDrawConfig {
            show_aabbs: false,
            show_contacts: false,
            show_joints: false,
            show_velocities: false,
            show_forces: false,
            ..Default::default()
        }
    }

    /// Enable all overlays.
    pub fn all_enabled() -> Self {
        DebugDrawConfig {
            show_aabbs: true,
            show_contacts: true,
            show_joints: true,
            show_velocities: true,
            show_forces: true,
            ..Default::default()
        }
    }

    /// Returns `true` if any overlay is enabled.
    pub fn any_enabled(&self) -> bool {
        self.show_aabbs
            || self.show_contacts
            || self.show_joints
            || self.show_velocities
            || self.show_forces
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
#[cfg(test)]
mod tests {
    use super::*;

    // --- SimulationConfig ---

    #[test]
    fn test_default_config_is_valid() {
        let cfg = SimulationConfig::default();
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn test_config_zero_timestep_is_invalid() {
        let cfg = SimulationConfig {
            timestep: 0.0,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_config_negative_timestep_is_invalid() {
        let cfg = SimulationConfig {
            timestep: -0.01,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_config_zero_iterations_is_invalid() {
        let cfg = SimulationConfig {
            iterations: 0,
            ..Default::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn test_config_builder_chain() {
        let cfg = SimulationConfig::new(1.0 / 120.0, [0.0, -9.81, 0.0])
            .with_iterations(20)
            .with_sleeping(0.05, 0.05)
            .with_broadphase(BroadphaseType::SweepAndPrune)
            .with_ccd(true);
        assert_eq!(cfg.iterations, 20);
        assert!(cfg.ccd_enabled);
        assert_eq!(cfg.broadphase, BroadphaseType::SweepAndPrune);
    }

    #[test]
    fn test_config_serialization_roundtrip() {
        let cfg = SimulationConfig::default();
        let json = serde_json::to_string(&cfg).unwrap();
        let cfg2: SimulationConfig = serde_json::from_str(&json).unwrap();
        assert!((cfg.timestep - cfg2.timestep).abs() < 1e-15);
        assert_eq!(cfg.iterations, cfg2.iterations);
    }

    // --- SimulationState ---

    #[test]
    fn test_state_running_is_active() {
        assert!(SimulationState::Running.is_active());
    }

    #[test]
    fn test_state_paused_is_not_active() {
        assert!(!SimulationState::Paused.is_active());
    }

    #[test]
    fn test_state_recording_is_recording() {
        assert!(SimulationState::Recording.is_recording());
    }

    #[test]
    fn test_state_labels() {
        assert_eq!(SimulationState::Running.label(), "running");
        assert_eq!(SimulationState::Paused.label(), "paused");
        assert_eq!(SimulationState::Stepping.label(), "stepping");
        assert_eq!(SimulationState::Rewinding.label(), "rewinding");
        assert_eq!(SimulationState::Recording.label(), "recording");
    }

    // --- SimulationController ---

    #[test]
    fn test_controller_starts_paused() {
        let ctrl = SimulationController::new(SimulationConfig::default());
        assert_eq!(ctrl.state, SimulationState::Paused);
    }

    #[test]
    fn test_controller_no_step_while_paused() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.step(1.0);
        assert_eq!(ctrl.stats.step_count, 0);
    }

    #[test]
    fn test_controller_steps_when_running() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.resume();
        ctrl.step(1.0 / 60.0);
        assert!(ctrl.stats.step_count >= 1);
    }

    #[test]
    fn test_controller_pause_resume() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.resume();
        assert_eq!(ctrl.state, SimulationState::Running);
        ctrl.pause();
        assert_eq!(ctrl.state, SimulationState::Paused);
    }

    #[test]
    fn test_controller_reset_clears_stats() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.resume();
        ctrl.step(1.0);
        assert!(ctrl.stats.step_count > 0);
        ctrl.reset();
        assert_eq!(ctrl.stats.step_count, 0);
        assert_eq!(ctrl.state, SimulationState::Paused);
    }

    #[test]
    fn test_controller_step_once_pauses_after() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.step_once();
        ctrl.step(ctrl.config.timestep * 2.0);
        assert_eq!(ctrl.state, SimulationState::Paused);
    }

    #[test]
    fn test_controller_update_counts() {
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.update_counts(100, 50);
        assert_eq!(ctrl.stats.body_count, 100);
        assert_eq!(ctrl.stats.constraint_count, 50);
    }

    // --- StepResult ---

    #[test]
    fn test_step_result_has_events() {
        let mut r = StepResult::new();
        assert!(!r.has_events());
        r.new_contacts = 1;
        assert!(r.has_events());
    }

    // --- ReplayBuffer ---

    #[test]
    fn test_replay_buffer_record_and_playback() {
        let mut buf = ReplayBuffer::new(10);
        let frame = ReplayFrame::new(0.0, 0);
        buf.record_state(frame);
        assert_eq!(buf.len(), 1);
        assert!(buf.playback_frame().is_some());
    }

    #[test]
    fn test_replay_buffer_capacity_evicts_oldest() {
        let mut buf = ReplayBuffer::new(3);
        for i in 0..5u64 {
            buf.record_state(ReplayFrame::new(i as f64, i));
        }
        assert_eq!(buf.len(), 3);
        assert_eq!(buf.frames[0].frame_index, 2);
    }

    #[test]
    fn test_replay_buffer_seek() {
        let mut buf = ReplayBuffer::new(10);
        for i in 0..5u64 {
            buf.record_state(ReplayFrame::new(i as f64 * 0.1, i));
        }
        let idx = buf.seek_to_time(0.25);
        assert_eq!(idx, 2); // 0.2 is nearest
    }

    #[test]
    fn test_replay_buffer_export_import() {
        let mut buf = ReplayBuffer::new(5);
        buf.record_state(ReplayFrame::new(1.0, 0));
        let json = buf.export_frames();
        let mut buf2 = ReplayBuffer::new(5);
        buf2.import_frames(&json).unwrap();
        assert_eq!(buf2.len(), 1);
    }

    #[test]
    fn test_replay_buffer_cursor_advance_rewind() {
        let mut buf = ReplayBuffer::new(5);
        buf.record_state(ReplayFrame::new(0.0, 0));
        buf.record_state(ReplayFrame::new(0.1, 1));
        assert!(buf.advance_cursor());
        assert_eq!(buf.cursor, 1);
        assert!(buf.rewind_cursor());
        assert_eq!(buf.cursor, 0);
        assert!(!buf.rewind_cursor()); // already at start
    }

    // --- PhysicsProfiler ---

    #[test]
    fn test_profiler_records_samples() {
        let mut prof = PhysicsProfiler::new(10);
        let s = PhysicsProfilerSample {
            broadphase_ms: 1.0,
            narrowphase_ms: 2.0,
            solver_ms: 3.0,
            integration_ms: 0.5,
            ..Default::default()
        };
        prof.record(s);
        assert_eq!(prof.sample_count(), 1);
        assert!((prof.latest().unwrap().total_ms - 6.5).abs() < 1e-10);
    }

    #[test]
    fn test_profiler_peak_tracking() {
        let mut prof = PhysicsProfiler::new(10);
        for total in [5.0f64, 10.0, 3.0] {
            prof.record(PhysicsProfilerSample {
                solver_ms: total,
                ..Default::default()
            });
        }
        assert!((prof.peak_total_ms - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_profiler_reset() {
        let mut prof = PhysicsProfiler::new(10);
        prof.record(PhysicsProfilerSample::default());
        prof.reset();
        assert_eq!(prof.sample_count(), 0);
        assert_eq!(prof.peak_total_ms, 0.0);
    }

    // --- SimulationExporter ---

    #[test]
    fn test_exporter_json_roundtrip() {
        let exp = SimulationExporter::new();
        let mut snap = SimulationSnapshot::new(1.5);
        snap.positions = vec![1.0, 2.0, 3.0];
        snap.masses = vec![1.0];
        let json = exp.export_json(&snap);
        let snap2 = exp.import_state(&json).unwrap();
        assert!((snap2.time - 1.5).abs() < 1e-10);
        assert_eq!(snap2.positions.len(), 3);
    }

    #[test]
    fn test_exporter_binary_roundtrip_header() {
        let exp = SimulationExporter::new();
        let snap = SimulationSnapshot::new(2.72);
        let bytes = exp.export_binary(&snap);
        // First 8 bytes = time as little-endian f64
        let time = f64::from_le_bytes(bytes[0..8].try_into().unwrap());
        assert!((time - 2.72).abs() < 1e-10);
    }

    #[test]
    fn test_exporter_vtk_contains_header() {
        let exp = SimulationExporter::new();
        let mut snap = SimulationSnapshot::new(0.0);
        snap.positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let vtk = exp.export_vtk(&snap, 0);
        assert!(vtk.contains("vtk DataFile"));
        assert!(vtk.contains("POINTS 2"));
    }

    #[test]
    fn test_exporter_vtk_series() {
        let mut exp = SimulationExporter::new();
        exp.push_snapshot(SimulationSnapshot::new(0.0));
        exp.push_snapshot(SimulationSnapshot::new(1.0));
        let series = exp.export_vtk_series();
        assert_eq!(series.len(), 2);
    }

    // --- TimeIntegrator ---

    #[test]
    fn test_integrator_euler_free_fall() {
        let integr = TimeIntegrator::new(IntegratorType::Euler);
        let (pos, vel) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        assert!((vel - (-9.81)).abs() < 1e-10);
        assert!((pos - 0.0).abs() < 1e-10); // Euler: pos updates first
    }

    #[test]
    fn test_integrator_semi_implicit_free_fall() {
        let integr = TimeIntegrator::new(IntegratorType::SemiImplicit);
        let (pos, vel) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        assert!((vel - (-9.81)).abs() < 1e-10);
        // position should be non-zero with semi-implicit
        assert!(pos < 0.0);
    }

    #[test]
    fn test_integrator_verlet() {
        let integr = TimeIntegrator::new(IntegratorType::Verlet);
        let (pos, _vel) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        // 0 + 0*1 + 0.5*(-9.81)*1^2 = -4.905
        assert!((pos - (-4.905)).abs() < 1e-6);
    }

    #[test]
    fn test_integrator_rk4() {
        let integr = TimeIntegrator::new(IntegratorType::Rk4);
        let (_pos, vel) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        assert!((vel - (-9.81)).abs() < 1e-10);
    }

    #[test]
    fn test_integrator_substeps() {
        let integr = TimeIntegrator::new(IntegratorType::SemiImplicit).with_substeps(10);
        let (pos1, _) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        // Both methods must produce negative position under gravity
        assert!(pos1 < 0.0);
        // Final velocity should be approximately g*t = -9.81 m/s
        let (_, vel1) = integr.integrate_1d(0.0, 0.0, -9.81, 1.0);
        assert!((vel1 - (-9.81)).abs() < 1e-10);
    }

    #[test]
    fn test_integrator_symplectic_flag() {
        assert!(IntegratorType::SemiImplicit.is_symplectic());
        assert!(IntegratorType::Verlet.is_symplectic());
        assert!(!IntegratorType::Euler.is_symplectic());
        assert!(!IntegratorType::Rk4.is_symplectic());
    }

    // --- DebugDrawConfig ---

    #[test]
    fn test_debug_draw_disabled() {
        let cfg = DebugDrawConfig::disabled();
        assert!(!cfg.any_enabled());
    }

    #[test]
    fn test_debug_draw_all_enabled() {
        let cfg = DebugDrawConfig::all_enabled();
        assert!(cfg.any_enabled());
        assert!(cfg.show_aabbs);
        assert!(cfg.show_contacts);
        assert!(cfg.show_joints);
    }

    #[test]
    fn test_debug_draw_serialization() {
        let cfg = DebugDrawConfig::default();
        let json = serde_json::to_string(&cfg).unwrap();
        let cfg2: DebugDrawConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(cfg.show_contacts, cfg2.show_contacts);
    }

    // --- BroadphaseType ---

    #[test]
    fn test_broadphase_default() {
        assert_eq!(BroadphaseType::default(), BroadphaseType::Bvh);
    }

    // --- SimulationStats ---

    #[test]
    fn test_stats_advance() {
        let mut stats = SimulationStats::new();
        stats.advance(1.0 / 60.0);
        assert_eq!(stats.step_count, 1);
        assert!((stats.elapsed_time - 1.0 / 60.0).abs() < 1e-15);
    }
}
