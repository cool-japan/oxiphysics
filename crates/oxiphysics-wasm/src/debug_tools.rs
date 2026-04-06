// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly debug and inspection tools.
//!
//! Provides debug draw lists, a physics inspector, a performance HUD,
//! runtime assertions, a replay controller, a scene explorer, a logger,
//! and a test harness—all designed for easy JavaScript interop.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ─────────────────────────────────────────────────────────────────────────────
// Colour helpers
// ─────────────────────────────────────────────────────────────────────────────

/// RGBA colour represented as four u8 components.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rgba {
    /// Red channel (0–255).
    pub r: u8,
    /// Green channel (0–255).
    pub g: u8,
    /// Blue channel (0–255).
    pub b: u8,
    /// Alpha channel (0–255).
    pub a: u8,
}

impl Rgba {
    /// Construct from components.
    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }

    /// Opaque red.
    pub fn red() -> Self {
        Self::new(255, 0, 0, 255)
    }
    /// Opaque green.
    pub fn green() -> Self {
        Self::new(0, 255, 0, 255)
    }
    /// Opaque blue.
    pub fn blue() -> Self {
        Self::new(0, 0, 255, 255)
    }
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(255, 255, 255, 255)
    }
    /// Opaque yellow.
    pub fn yellow() -> Self {
        Self::new(255, 255, 0, 255)
    }
    /// Transparent black.
    pub fn transparent() -> Self {
        Self::new(0, 0, 0, 0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmDebugConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Runtime debug visualisation configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmDebugConfig {
    /// Draw axis-aligned bounding boxes.
    pub show_aabbs: bool,
    /// Draw velocity vectors.
    pub show_velocities: bool,
    /// Draw contact points.
    pub show_contacts: bool,
    /// Draw joint anchors and axes.
    pub show_joints: bool,
    /// Colour for AABB outlines.
    pub aabb_color: Rgba,
    /// Colour for velocity arrows.
    pub velocity_color: Rgba,
    /// Colour for contact normals.
    pub contact_color: Rgba,
    /// Scale factor for velocity arrows.
    pub velocity_scale: f64,
}

impl Default for WasmDebugConfig {
    fn default() -> Self {
        Self {
            show_aabbs: false,
            show_velocities: false,
            show_contacts: false,
            show_joints: false,
            aabb_color: Rgba::green(),
            velocity_color: Rgba::yellow(),
            contact_color: Rgba::red(),
            velocity_scale: 0.1,
        }
    }
}

impl WasmDebugConfig {
    /// Enable all visual overlays.
    pub fn all_enabled() -> Self {
        Self {
            show_aabbs: true,
            show_velocities: true,
            show_contacts: true,
            show_joints: true,
            ..Default::default()
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmDebugDraw primitives
// ─────────────────────────────────────────────────────────────────────────────

/// A single debug draw call.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmDrawCall {
    /// Draw a line segment.
    Line {
        /// Start point.
        start: [f64; 3],
        /// End point.
        end: [f64; 3],
        /// Line colour.
        color: Rgba,
    },
    /// Draw a wire-frame sphere.
    Sphere {
        /// Centre of the sphere.
        center: [f64; 3],
        /// Sphere radius.
        radius: f64,
        /// Colour.
        color: Rgba,
    },
    /// Draw a wire-frame axis-aligned box.
    Box {
        /// Minimum corner.
        min: [f64; 3],
        /// Maximum corner.
        max: [f64; 3],
        /// Colour.
        color: Rgba,
    },
    /// Draw an arrow (line + arrowhead).
    Arrow {
        /// Tail of the arrow.
        from: [f64; 3],
        /// Tip of the arrow.
        to: [f64; 3],
        /// Colour.
        color: Rgba,
    },
    /// Draw a text label.
    Text {
        /// World-space position.
        position: [f64; 3],
        /// Text to display.
        label: String,
        /// Colour.
        color: Rgba,
    },
}

/// A list of debug draw calls for a single frame.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmDebugDrawList {
    /// All draw calls.
    pub calls: Vec<WasmDrawCall>,
}

impl WasmDebugDrawList {
    /// Create an empty draw list.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a draw call.
    pub fn push(&mut self, call: WasmDrawCall) {
        self.calls.push(call);
    }

    /// Number of draw calls.
    pub fn len(&self) -> usize {
        self.calls.len()
    }

    /// True if no draw calls.
    pub fn is_empty(&self) -> bool {
        self.calls.is_empty()
    }

    /// Clear all draw calls.
    pub fn clear(&mut self) {
        self.calls.clear();
    }

    /// Serialise to JSON for JavaScript consumption.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(self)
    }
}

/// Debug draw helper with fluent API.
#[derive(Debug, Clone, Default)]
pub struct WasmDebugDraw {
    /// Accumulated draw list.
    pub list: WasmDebugDrawList,
}

impl WasmDebugDraw {
    /// Create a new draw helper.
    pub fn new() -> Self {
        Self::default()
    }

    /// Draw a line segment.
    pub fn draw_line(&mut self, start: [f64; 3], end: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Line { start, end, color });
    }

    /// Draw a sphere.
    pub fn draw_sphere(&mut self, center: [f64; 3], radius: f64, color: Rgba) {
        self.list.push(WasmDrawCall::Sphere {
            center,
            radius,
            color,
        });
    }

    /// Draw an AABB box.
    pub fn draw_box(&mut self, min: [f64; 3], max: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Box { min, max, color });
    }

    /// Draw an arrow.
    pub fn draw_arrow(&mut self, from: [f64; 3], to: [f64; 3], color: Rgba) {
        self.list.push(WasmDrawCall::Arrow { from, to, color });
    }

    /// Draw a text label.
    pub fn draw_text(&mut self, position: [f64; 3], label: impl Into<String>, color: Rgba) {
        self.list.push(WasmDrawCall::Text {
            position,
            label: label.into(),
            color,
        });
    }

    /// Flush (clone) the draw list and clear the internal one.
    pub fn flush(&mut self) -> WasmDebugDrawList {
        let out = self.list.clone();
        self.list.clear();
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmPhysicsInspector
// ─────────────────────────────────────────────────────────────────────────────

/// A body state snapshot returned by the inspector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyInspectResult {
    /// Body handle.
    pub id: u32,
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Angular velocity `[wx, wy, wz]`.
    pub ang_vel: [f64; 3],
    /// Kinetic energy (J).
    pub kinetic_energy: f64,
    /// Whether sleeping.
    pub sleeping: bool,
    /// Whether static.
    pub is_static: bool,
}

/// A joint state snapshot.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JointInspectResult {
    /// Joint handle.
    pub id: u32,
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
    /// Joint type name.
    pub joint_type: String,
    /// Reaction force `[fx, fy, fz]`.
    pub reaction_force: [f64; 3],
    /// Motor enabled.
    pub motor_enabled: bool,
}

/// Inspects simulation objects and returns JSON-serialisable snapshots.
#[derive(Debug, Clone, Default)]
pub struct WasmPhysicsInspector {
    /// Mock body registry: `(id, position, velocity, ang_vel, sleeping, static)`.
    #[allow(clippy::type_complexity)]
    bodies: Vec<(u32, [f64; 3], [f64; 3], [f64; 3], bool, bool)>,
    /// Mock joint registry: `(id, body_a, body_b, type_name)`.
    joints: Vec<(u32, u32, u32, String)>,
}

impl WasmPhysicsInspector {
    /// Create a new inspector.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a body for inspection.
    pub fn register_body(
        &mut self,
        id: u32,
        position: [f64; 3],
        velocity: [f64; 3],
        ang_vel: [f64; 3],
        sleeping: bool,
        is_static: bool,
    ) {
        self.bodies.retain(|(bid, _, _, _, _, _)| *bid != id);
        self.bodies
            .push((id, position, velocity, ang_vel, sleeping, is_static));
    }

    /// Register a joint for inspection.
    pub fn register_joint(
        &mut self,
        id: u32,
        body_a: u32,
        body_b: u32,
        type_name: impl Into<String>,
    ) {
        self.joints.retain(|(jid, _, _, _)| *jid != id);
        self.joints.push((id, body_a, body_b, type_name.into()));
    }

    /// Inspect a body; returns JSON string or `None` if not found.
    pub fn inspect_body(&self, id: u32) -> Option<String> {
        let &(_id, pos, vel, ang, sleeping, is_static) =
            self.bodies.iter().find(|(bid, _, _, _, _, _)| *bid == id)?;
        let ke = {
            let v2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
            0.5 * 1.0 * v2 // assume mass = 1 for mock
        };
        let result = BodyInspectResult {
            id,
            position: pos,
            velocity: vel,
            ang_vel: ang,
            kinetic_energy: ke,
            sleeping,
            is_static,
        };
        serde_json::to_string_pretty(&result).ok()
    }

    /// Inspect a joint; returns JSON or `None`.
    pub fn inspect_joint(&self, id: u32) -> Option<String> {
        let &(jid, body_a, body_b, ref jtype) =
            self.joints.iter().find(|(jid, _, _, _)| *jid == id)?;
        let result = JointInspectResult {
            id: jid,
            body_a,
            body_b,
            joint_type: jtype.clone(),
            reaction_force: [0.0; 3],
            motor_enabled: false,
        };
        serde_json::to_string_pretty(&result).ok()
    }

    /// Inspect an island (mock: returns body count JSON).
    pub fn inspect_island(&self, island_id: u32) -> String {
        format!(
            "{{\"island_id\":{},\"body_count\":{}}}",
            island_id,
            self.bodies.len()
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmPerformanceHud
// ─────────────────────────────────────────────────────────────────────────────

/// Per-frame timing record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameTiming {
    /// Frame index.
    pub frame: u64,
    /// Frame duration (ms).
    pub ms: f64,
}

/// Rolling-average performance HUD.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmPerformanceHud {
    /// History of frame timings.
    pub history: Vec<FrameTiming>,
    /// Capacity limit.
    pub capacity: usize,
    /// Current frame counter.
    pub frame: u64,
}

impl WasmPerformanceHud {
    /// Create with capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            ..Default::default()
        }
    }

    /// Record a frame timing.
    pub fn record_frame(&mut self, ms: f64) {
        if self.capacity > 0 && self.history.len() >= self.capacity {
            self.history.remove(0);
        }
        self.history.push(FrameTiming {
            frame: self.frame,
            ms,
        });
        self.frame += 1;
    }

    /// Rolling average frame time.
    pub fn rolling_average(&self) -> f64 {
        if self.history.is_empty() {
            return 0.0;
        }
        self.history.iter().map(|f| f.ms).sum::<f64>() / self.history.len() as f64
    }

    /// Peak frame time.
    pub fn peak_ms(&self) -> f64 {
        self.history.iter().map(|f| f.ms).fold(0.0_f64, f64::max)
    }

    /// Export CSV string.
    pub fn export_csv(&self) -> String {
        let mut out = String::from("frame,ms\n");
        for f in &self.history {
            out.push_str(&format!("{},{:.4}\n", f.frame, f.ms));
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmPhysicsAssert
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a physics assertion check.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssertResult {
    /// Whether the assertion passed.
    pub passed: bool,
    /// Human-readable message.
    pub message: String,
}

impl AssertResult {
    fn pass(msg: impl Into<String>) -> Self {
        Self {
            passed: true,
            message: msg.into(),
        }
    }
    fn fail(msg: impl Into<String>) -> Self {
        Self {
            passed: false,
            message: msg.into(),
        }
    }
}

/// Runtime physics assertion checker.
#[derive(Debug, Clone, Default)]
pub struct WasmPhysicsAssert {
    /// History of all assertion results.
    pub results: Vec<AssertResult>,
}

impl WasmPhysicsAssert {
    /// Create a new asserter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Assert that energy is below `max_energy`.
    pub fn assert_energy_lt(&mut self, energy: f64, max_energy: f64) -> bool {
        if energy < max_energy {
            self.results.push(AssertResult::pass(format!(
                "energy {energy:.4} < {max_energy:.4} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "energy {energy:.4} >= {max_energy:.4} FAIL"
            )));
            false
        }
    }

    /// Assert that penetration depth is below epsilon.
    pub fn assert_penetration_lt(&mut self, depth: f64, eps: f64) -> bool {
        if depth < eps {
            self.results.push(AssertResult::pass(format!(
                "depth {depth:.6} < {eps:.6} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "depth {depth:.6} >= {eps:.6} FAIL"
            )));
            false
        }
    }

    /// Assert that velocity magnitude is below `max_vel`.
    pub fn assert_velocity_lt(&mut self, vel_mag: f64, max_vel: f64) -> bool {
        if vel_mag < max_vel {
            self.results.push(AssertResult::pass(format!(
                "vel {vel_mag:.4} < {max_vel:.4} OK"
            )));
            true
        } else {
            self.results.push(AssertResult::fail(format!(
                "vel {vel_mag:.4} >= {max_vel:.4} FAIL"
            )));
            false
        }
    }

    /// Number of failing assertions.
    pub fn failure_count(&self) -> usize {
        self.results.iter().filter(|r| !r.passed).count()
    }

    /// True if all assertions passed.
    pub fn all_passed(&self) -> bool {
        self.failure_count() == 0
    }

    /// Clear results.
    pub fn clear(&mut self) {
        self.results.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmReplayController
// ─────────────────────────────────────────────────────────────────────────────

/// A recorded frame in the replay buffer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReplayFrame {
    /// Frame index.
    pub index: usize,
    /// Simulation time (s).
    pub time: f64,
    /// Serialised state (body positions as JSON).
    pub state_json: String,
}

/// Controls replay of a recorded simulation.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmReplayController {
    /// Recorded frames.
    pub frames: Vec<ReplayFrame>,
    /// Current playback position.
    pub current_frame: usize,
    /// Whether playback is active.
    pub playing: bool,
}

impl WasmReplayController {
    /// Create a new controller.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a frame.
    pub fn record_frame(&mut self, time: f64, state_json: impl Into<String>) {
        let index = self.frames.len();
        self.frames.push(ReplayFrame {
            index,
            time,
            state_json: state_json.into(),
        });
    }

    /// Step forward one frame.
    pub fn step_forward(&mut self) -> Option<&ReplayFrame> {
        if self.current_frame + 1 < self.frames.len() {
            self.current_frame += 1;
        }
        self.frames.get(self.current_frame)
    }

    /// Step backward one frame.
    pub fn step_backward(&mut self) -> Option<&ReplayFrame> {
        if self.current_frame > 0 {
            self.current_frame -= 1;
        }
        self.frames.get(self.current_frame)
    }

    /// Start playback.
    pub fn play(&mut self) {
        self.playing = true;
    }

    /// Pause playback.
    pub fn pause(&mut self) {
        self.playing = false;
    }

    /// Seek to a specific frame index.
    pub fn seek_to_frame(&mut self, frame: usize) -> Option<&ReplayFrame> {
        self.current_frame = frame.min(self.frames.len().saturating_sub(1));
        self.frames.get(self.current_frame)
    }

    /// Current frame count.
    pub fn total_frames(&self) -> usize {
        self.frames.len()
    }

    /// Current playback frame.
    pub fn current(&self) -> Option<&ReplayFrame> {
        self.frames.get(self.current_frame)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmSceneExplorer
// ─────────────────────────────────────────────────────────────────────────────

/// A node in the scene hierarchy.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneNode {
    /// Node (body) id.
    pub id: u32,
    /// Optional name.
    pub name: Option<String>,
    /// Parent id (None = root).
    pub parent: Option<u32>,
    /// Child node ids.
    pub children: Vec<u32>,
}

impl SceneNode {
    /// Create a new scene node.
    pub fn new(id: u32, name: Option<String>, parent: Option<u32>) -> Self {
        Self {
            id,
            name,
            parent,
            children: Vec::new(),
        }
    }
}

/// Enumerates and queries all scene objects.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WasmSceneExplorer {
    /// Scene nodes.
    pub nodes: Vec<SceneNode>,
    /// Joint records `(joint_id, body_a, body_b)`.
    pub joints: Vec<(u32, u32, u32)>,
}

impl WasmSceneExplorer {
    /// Create an empty explorer.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a body node.
    pub fn add_node(&mut self, id: u32, name: Option<String>, parent: Option<u32>) {
        // Link to parent
        if let Some(pid) = parent
            && let Some(p) = self.nodes.iter_mut().find(|n| n.id == pid)
        {
            p.children.push(id);
        }
        self.nodes.push(SceneNode::new(id, name, parent));
    }

    /// Register a joint.
    pub fn add_joint(&mut self, joint_id: u32, body_a: u32, body_b: u32) {
        self.joints.push((joint_id, body_a, body_b));
    }

    /// Enumerate all body ids.
    pub fn enumerate_bodies(&self) -> Vec<u32> {
        self.nodes.iter().map(|n| n.id).collect()
    }

    /// Enumerate all joint ids.
    pub fn enumerate_joints(&self) -> Vec<u32> {
        self.joints.iter().map(|&(jid, _, _)| jid).collect()
    }

    /// Get hierarchy as JSON.
    pub fn get_hierarchy(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.nodes)
    }

    /// Find a body by name; returns id if found.
    pub fn find_by_name(&self, name: &str) -> Option<u32> {
        self.nodes.iter().find_map(|n| {
            if n.name.as_deref() == Some(name) {
                Some(n.id)
            } else {
                None
            }
        })
    }

    /// Root nodes (nodes with no parent).
    pub fn roots(&self) -> Vec<u32> {
        self.nodes
            .iter()
            .filter(|n| n.parent.is_none())
            .map(|n| n.id)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmPhysicsLogger
// ─────────────────────────────────────────────────────────────────────────────

/// Log level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LogLevel {
    /// Error — always shown.
    Error = 0,
    /// Warning.
    Warn = 1,
    /// Informational.
    Info = 2,
    /// Verbose debug output.
    Debug = 3,
}

/// A single log entry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LogEntry {
    /// Severity.
    pub level: LogLevel,
    /// Category tag.
    pub category: String,
    /// Log message.
    pub message: String,
    /// Simulation time when logged.
    pub time: f64,
}

impl LogEntry {
    fn new(
        level: LogLevel,
        category: impl Into<String>,
        message: impl Into<String>,
        time: f64,
    ) -> Self {
        Self {
            level,
            category: category.into(),
            message: message.into(),
            time,
        }
    }
}

/// Ring-buffer logger with category filtering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmPhysicsLogger {
    /// Ring buffer of entries.
    pub entries: Vec<LogEntry>,
    /// Ring buffer capacity.
    pub capacity: usize,
    /// Minimum level to log.
    pub min_level: LogLevel,
    /// Optional category filter (`None` = all categories).
    pub category_filter: Option<String>,
}

impl WasmPhysicsLogger {
    /// Create a logger.
    pub fn new(capacity: usize, min_level: LogLevel) -> Self {
        Self {
            entries: Vec::new(),
            capacity,
            min_level,
            category_filter: None,
        }
    }

    /// Filter to a specific category.
    pub fn set_category_filter(&mut self, category: impl Into<String>) {
        self.category_filter = Some(category.into());
    }

    /// Remove category filter.
    pub fn clear_filter(&mut self) {
        self.category_filter = None;
    }

    /// Log a message.
    pub fn log(
        &mut self,
        level: LogLevel,
        category: impl Into<String>,
        message: impl Into<String>,
        time: f64,
    ) {
        if level > self.min_level {
            return;
        }
        let cat: String = category.into();
        if let Some(ref f) = self.category_filter
            && &cat != f
        {
            return;
        }
        if self.capacity > 0 && self.entries.len() >= self.capacity {
            self.entries.remove(0);
        }
        self.entries.push(LogEntry::new(level, cat, message, time));
    }

    /// Shorthand error log.
    pub fn error(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Error, cat, msg, time);
    }

    /// Shorthand warn log.
    pub fn warn(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Warn, cat, msg, time);
    }

    /// Shorthand info log.
    pub fn info(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Info, cat, msg, time);
    }

    /// Shorthand debug log.
    pub fn debug(&mut self, cat: impl Into<String>, msg: impl Into<String>, time: f64) {
        self.log(LogLevel::Debug, cat, msg, time);
    }

    /// Entries at or above a given level.
    pub fn entries_at_level(&self, level: LogLevel) -> Vec<&LogEntry> {
        self.entries.iter().filter(|e| e.level <= level).collect()
    }

    /// Flush all entries.
    pub fn flush(&mut self) -> Vec<LogEntry> {
        std::mem::take(&mut self.entries)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WasmTestHarness
// ─────────────────────────────────────────────────────────────────────────────

/// Expected state for a single body.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExpectedBodyState {
    /// Body id.
    pub id: u32,
    /// Expected position tolerance per axis.
    pub position: [f64; 3],
    /// Position tolerance.
    pub position_tol: f64,
    /// Expected velocity.
    pub velocity: [f64; 3],
    /// Velocity tolerance.
    pub velocity_tol: f64,
}

/// Result of a test scenario.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct TestResult {
    /// Whether the test passed.
    pub passed: bool,
    /// Messages from assertions.
    pub messages: Vec<String>,
}

impl TestResult {
    fn pass() -> Self {
        Self {
            passed: true,
            messages: Vec::new(),
        }
    }

    fn fail(messages: Vec<String>) -> Self {
        Self {
            passed: false,
            messages,
        }
    }
}

/// Runs a physics scenario and asserts expected state.
#[derive(Debug, Clone, Default)]
pub struct WasmTestHarness {
    /// Baseline states (name -> JSON snapshot).
    baselines: Vec<(String, String)>,
    /// Last test result.
    pub last_result: Option<TestResult>,
}

impl WasmTestHarness {
    /// Create a new harness.
    pub fn new() -> Self {
        Self::default()
    }

    /// Save a named baseline state.
    pub fn save_baseline(&mut self, name: impl Into<String>, snapshot_json: impl Into<String>) {
        let n = name.into();
        self.baselines.retain(|(bn, _)| bn != &n);
        self.baselines.push((n, snapshot_json.into()));
    }

    /// Load a baseline by name.
    pub fn load_baseline(&self, name: &str) -> Option<&str> {
        self.baselines
            .iter()
            .find_map(|(n, s)| if n == name { Some(s.as_str()) } else { None })
    }

    /// Assert that a set of actual body states matches expectations.
    /// Each actual state is `(id, position, velocity)`.
    #[allow(clippy::needless_range_loop)]
    pub fn assert_expected(
        &mut self,
        actual: &[(u32, [f64; 3], [f64; 3])],
        expected: &[ExpectedBodyState],
    ) -> bool {
        let mut messages = Vec::new();
        for exp in expected {
            if let Some(&(_id, pos, vel)) = actual.iter().find(|&&(id, _, _)| id == exp.id) {
                // Check position
                for i in 0..3 {
                    let diff = (pos[i] - exp.position[i]).abs();
                    if diff > exp.position_tol {
                        messages.push(format!(
                            "Body {} pos[{i}] diff={diff:.6} > tol={:.6}",
                            exp.id, exp.position_tol
                        ));
                    }
                }
                // Check velocity
                for i in 0..3 {
                    let diff = (vel[i] - exp.velocity[i]).abs();
                    if diff > exp.velocity_tol {
                        messages.push(format!(
                            "Body {} vel[{i}] diff={diff:.6} > tol={:.6}",
                            exp.id, exp.velocity_tol
                        ));
                    }
                }
            } else {
                messages.push(format!("Body {} not found in actual state", exp.id));
            }
        }
        let passed = messages.is_empty();
        self.last_result = Some(if passed {
            TestResult::pass()
        } else {
            TestResult::fail(messages)
        });
        passed
    }

    /// Diff current state against a saved baseline.
    pub fn diff_against_baseline(&self, name: &str, current_json: &str) -> Option<String> {
        let baseline = self.load_baseline(name)?;
        if baseline == current_json {
            Some("MATCH".to_string())
        } else {
            Some(format!(
                "DIFFER\nBaseline: {baseline}\nCurrent: {current_json}"
            ))
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- Rgba ---

    #[test]
    fn test_rgba_constructors() {
        assert_eq!(Rgba::red(), Rgba::new(255, 0, 0, 255));
        assert_eq!(Rgba::transparent(), Rgba::new(0, 0, 0, 0));
    }

    // --- WasmDebugConfig ---

    #[test]
    fn test_debug_config_default() {
        let cfg = WasmDebugConfig::default();
        assert!(!cfg.show_aabbs);
        assert!(!cfg.show_contacts);
    }

    #[test]
    fn test_debug_config_all_enabled() {
        let cfg = WasmDebugConfig::all_enabled();
        assert!(cfg.show_aabbs);
        assert!(cfg.show_velocities);
        assert!(cfg.show_contacts);
        assert!(cfg.show_joints);
    }

    // --- WasmDebugDraw ---

    #[test]
    fn test_draw_line() {
        let mut d = WasmDebugDraw::new();
        d.draw_line([0.0; 3], [1.0, 0.0, 0.0], Rgba::red());
        assert_eq!(d.list.len(), 1);
    }

    #[test]
    fn test_draw_sphere() {
        let mut d = WasmDebugDraw::new();
        d.draw_sphere([0.0; 3], 1.0, Rgba::blue());
        assert_eq!(d.list.len(), 1);
    }

    #[test]
    fn test_draw_box() {
        let mut d = WasmDebugDraw::new();
        d.draw_box([-1.0; 3], [1.0; 3], Rgba::green());
        assert_eq!(d.list.len(), 1);
    }

    #[test]
    fn test_draw_arrow() {
        let mut d = WasmDebugDraw::new();
        d.draw_arrow([0.0; 3], [0.0, 1.0, 0.0], Rgba::yellow());
        assert_eq!(d.list.len(), 1);
    }

    #[test]
    fn test_draw_text() {
        let mut d = WasmDebugDraw::new();
        d.draw_text([0.0; 3], "hello", Rgba::white());
        assert_eq!(d.list.len(), 1);
    }

    #[test]
    fn test_draw_flush_clears() {
        let mut d = WasmDebugDraw::new();
        d.draw_line([0.0; 3], [1.0, 0.0, 0.0], Rgba::red());
        let flushed = d.flush();
        assert_eq!(flushed.len(), 1);
        assert_eq!(d.list.len(), 0);
    }

    #[test]
    fn test_draw_list_to_json() {
        let mut list = WasmDebugDrawList::new();
        list.push(WasmDrawCall::Line {
            start: [0.0; 3],
            end: [1.0, 0.0, 0.0],
            color: Rgba::red(),
        });
        let json = list.to_json().unwrap();
        assert!(json.contains("Line"));
    }

    // --- WasmPhysicsInspector ---

    #[test]
    fn test_inspector_body() {
        let mut insp = WasmPhysicsInspector::new();
        insp.register_body(1, [1.0, 2.0, 3.0], [0.1, 0.0, 0.0], [0.0; 3], false, false);
        let json = insp.inspect_body(1).unwrap();
        assert!(json.contains("position"));
    }

    #[test]
    fn test_inspector_body_not_found() {
        let insp = WasmPhysicsInspector::new();
        assert!(insp.inspect_body(999).is_none());
    }

    #[test]
    fn test_inspector_joint() {
        let mut insp = WasmPhysicsInspector::new();
        insp.register_joint(0, 1, 2, "Hinge");
        let json = insp.inspect_joint(0).unwrap();
        assert!(json.contains("Hinge"));
    }

    #[test]
    fn test_inspector_island() {
        let insp = WasmPhysicsInspector::new();
        let s = insp.inspect_island(0);
        assert!(s.contains("island_id"));
    }

    // --- WasmPerformanceHud ---

    #[test]
    fn test_hud_record_and_avg() {
        let mut hud = WasmPerformanceHud::new(10);
        hud.record_frame(16.0);
        hud.record_frame(20.0);
        assert!((hud.rolling_average() - 18.0).abs() < 1e-10);
    }

    #[test]
    fn test_hud_eviction() {
        let mut hud = WasmPerformanceHud::new(2);
        hud.record_frame(1.0);
        hud.record_frame(2.0);
        hud.record_frame(3.0);
        assert_eq!(hud.history.len(), 2);
    }

    #[test]
    fn test_hud_csv() {
        let mut hud = WasmPerformanceHud::new(5);
        hud.record_frame(10.0);
        let csv = hud.export_csv();
        assert!(csv.contains("frame,ms"));
    }

    // --- WasmPhysicsAssert ---

    #[test]
    fn test_assert_energy_pass() {
        let mut a = WasmPhysicsAssert::new();
        assert!(a.assert_energy_lt(5.0, 10.0));
        assert!(a.all_passed());
    }

    #[test]
    fn test_assert_energy_fail() {
        let mut a = WasmPhysicsAssert::new();
        assert!(!a.assert_energy_lt(15.0, 10.0));
        assert_eq!(a.failure_count(), 1);
    }

    #[test]
    fn test_assert_penetration() {
        let mut a = WasmPhysicsAssert::new();
        assert!(a.assert_penetration_lt(0.001, 0.01));
        assert!(!a.assert_penetration_lt(0.1, 0.01));
    }

    #[test]
    fn test_assert_velocity() {
        let mut a = WasmPhysicsAssert::new();
        assert!(a.assert_velocity_lt(5.0, 100.0));
        assert!(!a.assert_velocity_lt(200.0, 100.0));
    }

    // --- WasmReplayController ---

    #[test]
    fn test_replay_record_seek() {
        let mut r = WasmReplayController::new();
        r.record_frame(0.0, "{}");
        r.record_frame(0.1, "{\"x\":1}");
        assert_eq!(r.total_frames(), 2);
        let f = r.seek_to_frame(1).unwrap();
        assert_eq!(f.index, 1);
    }

    #[test]
    fn test_replay_step_forward_backward() {
        let mut r = WasmReplayController::new();
        r.record_frame(0.0, "a");
        r.record_frame(0.1, "b");
        r.step_forward();
        assert_eq!(r.current_frame, 1);
        r.step_backward();
        assert_eq!(r.current_frame, 0);
    }

    #[test]
    fn test_replay_play_pause() {
        let mut r = WasmReplayController::new();
        r.play();
        assert!(r.playing);
        r.pause();
        assert!(!r.playing);
    }

    // --- WasmSceneExplorer ---

    #[test]
    fn test_scene_explorer_add_find() {
        let mut e = WasmSceneExplorer::new();
        e.add_node(1, Some("player".to_string()), None);
        e.add_node(2, Some("bullet".to_string()), Some(1));
        assert_eq!(e.find_by_name("player"), Some(1));
        assert_eq!(e.find_by_name("missing"), None);
    }

    #[test]
    fn test_scene_explorer_roots() {
        let mut e = WasmSceneExplorer::new();
        e.add_node(1, None, None);
        e.add_node(2, None, Some(1));
        let roots = e.roots();
        assert_eq!(roots, vec![1]);
    }

    #[test]
    fn test_scene_explorer_hierarchy_json() {
        let mut e = WasmSceneExplorer::new();
        e.add_node(1, Some("root".to_string()), None);
        let json = e.get_hierarchy().unwrap();
        assert!(json.contains("root"));
    }

    // --- WasmPhysicsLogger ---

    #[test]
    fn test_logger_levels() {
        let mut log = WasmPhysicsLogger::new(100, LogLevel::Info);
        log.info("physics", "step done", 0.0);
        log.debug("physics", "verbose", 0.0); // below min_level=Info => not logged
        assert_eq!(log.entries.len(), 1);
    }

    #[test]
    fn test_logger_category_filter() {
        let mut log = WasmPhysicsLogger::new(100, LogLevel::Debug);
        log.set_category_filter("collision");
        log.info("physics", "not me", 0.0);
        log.info("collision", "this one", 0.0);
        assert_eq!(log.entries.len(), 1);
    }

    #[test]
    fn test_logger_flush() {
        let mut log = WasmPhysicsLogger::new(10, LogLevel::Debug);
        log.error("sys", "oh no", 0.0);
        let drained = log.flush();
        assert_eq!(drained.len(), 1);
        assert_eq!(log.entries.len(), 0);
    }

    // --- WasmTestHarness ---

    #[test]
    fn test_harness_assert_pass() {
        let mut h = WasmTestHarness::new();
        let actual = vec![(1, [0.0, 5.0, 0.0], [0.0; 3])];
        let expected = vec![ExpectedBodyState {
            id: 1,
            position: [0.0, 5.0, 0.0],
            position_tol: 0.01,
            velocity: [0.0; 3],
            velocity_tol: 0.01,
        }];
        assert!(h.assert_expected(&actual, &expected));
    }

    #[test]
    fn test_harness_assert_fail() {
        let mut h = WasmTestHarness::new();
        let actual = vec![(1, [0.0, 100.0, 0.0], [0.0; 3])];
        let expected = vec![ExpectedBodyState {
            id: 1,
            position: [0.0, 5.0, 0.0],
            position_tol: 0.01,
            velocity: [0.0; 3],
            velocity_tol: 0.01,
        }];
        assert!(!h.assert_expected(&actual, &expected));
        let result = h.last_result.unwrap();
        assert!(!result.passed);
    }

    #[test]
    fn test_harness_baseline() {
        let mut h = WasmTestHarness::new();
        h.save_baseline("t1", "{\"x\":0}");
        let diff = h.diff_against_baseline("t1", "{\"x\":0}").unwrap();
        assert_eq!(diff, "MATCH");
    }
}
