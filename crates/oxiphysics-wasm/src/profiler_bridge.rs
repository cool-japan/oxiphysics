// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the hierarchical scoped profiler.
//!
//! Exposes a JSON-oriented surface around frame-level profiling with named
//! scopes, call counts, and nanosecond timing, suitable for use across the
//! WASM boundary.
//!
//! # Clock
//!
//! All timings come from [`crate::wasm_helpers::now_ms`] (`performance.now()`
//! on `wasm32`, a monotonic `Instant` baseline natively).  `std::time::Instant`
//! is deliberately not used directly: on `wasm32-unknown-unknown` it panics,
//! which aborts the whole WebAssembly instance.  When a target exposes no
//! clock at all, the affected timings are reported as `0` and
//! [`WasmFrameReport::clock_available`] is `false` instead of a made-up number.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, now_ms, to_js_value};

/// Nanoseconds between two [`now_ms`] samples.
///
/// Returns `None` when either endpoint is missing (the target has no clock) or
/// the interval is not finite.  A negative interval (a clock that stepped
/// backwards) is clamped to zero.
fn elapsed_ns(start: Option<f64>, end: Option<f64>) -> Option<u64> {
    let elapsed_ms = end? - start?;
    if !elapsed_ms.is_finite() {
        return None;
    }
    Some((elapsed_ms.max(0.0) * 1.0e6) as u64)
}

/// Serde default for [`WasmFrameReport::clock_available`]: reports produced
/// before the field existed were always measured with a real clock.
fn default_clock_available() -> bool {
    true
}

// ---------------------------------------------------------------------------
// WasmScopeNode
// ---------------------------------------------------------------------------

/// A node in the profiling tree, populated by [`WasmProfilerSession::end_frame`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmScopeNode {
    pub name: String,
    pub call_count: u64,
    pub total_ns: u64,
    pub max_ns: u64,
    pub children: Vec<WasmScopeNode>,
}

// ---------------------------------------------------------------------------
// WasmFrameReport
// ---------------------------------------------------------------------------

/// Per-frame profiling summary produced by [`WasmProfilerSession::end_frame`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmFrameReport {
    /// Total wall-clock nanoseconds for the frame (`0` when the frame was not
    /// begun with [`WasmProfilerSession::begin_frame`] or no clock exists).
    pub frame_ns: u64,
    /// `false` when the platform exposed no clock for at least one timing
    /// sample of this frame; every unmeasurable duration is then reported as
    /// `0` rather than fabricated.  Always `true` on browsers, workers and
    /// native targets.
    #[serde(default = "default_clock_available")]
    pub clock_available: bool,
    /// Root of the scope tree (name = `"frame"`).
    pub root: WasmScopeNode,
}

impl WasmFrameReport {
    /// Serialise this report as a JSON string.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }

    /// Render all scope nodes depth-first as CSV.
    ///
    /// Header: `name,call_count,total_ns,max_ns`
    pub fn to_csv(&self) -> String {
        let mut out = String::from("name,call_count,total_ns,max_ns\n");
        Self::node_to_csv(&self.root, &mut out);
        out
    }

    fn node_to_csv(node: &WasmScopeNode, out: &mut String) {
        out.push_str(&format!(
            "{},{},{},{}\n",
            node.name, node.call_count, node.total_ns, node.max_ns
        ));
        for child in &node.children {
            Self::node_to_csv(child, out);
        }
    }

    /// Render the scope tree as folded-stacks format (flamegraph-compatible).
    pub fn to_folded_stacks(&self) -> String {
        let mut out = String::new();
        Self::node_to_folded(&self.root, "frame", &mut out);
        out
    }

    fn node_to_folded(node: &WasmScopeNode, path: &str, out: &mut String) {
        out.push_str(&format!("{} {}\n", path, node.total_ns));
        for child in &node.children {
            let child_path = format!("{};{}", path, child.name);
            Self::node_to_folded(child, &child_path, out);
        }
    }
}

// ---------------------------------------------------------------------------
// Internal arena types
// ---------------------------------------------------------------------------

struct ArenaNode {
    name: String,
    parent: Option<usize>,
    call_count: u64,
    total_ns: u64,
    max_ns: u64,
}

// ---------------------------------------------------------------------------
// WasmScopeGuard — RAII scope guard
// ---------------------------------------------------------------------------

/// RAII guard returned by [`WasmProfilerSession::scope`].
///
/// On drop, records elapsed nanoseconds and pops the scope off the stack.
pub struct WasmScopeGuard<'a> {
    session: &'a mut WasmProfilerSession,
    /// Clock sample taken when the scope opened (`None` without a clock).
    start: Option<f64>,
    node_idx: usize,
}

impl Drop for WasmScopeGuard<'_> {
    fn drop(&mut self) {
        let elapsed_ns = self.session.elapsed_since(self.start);
        self.session.record_scope(self.node_idx, elapsed_ns);
        self.session.stack.pop();
    }
}

// ---------------------------------------------------------------------------
// WasmProfilerSession
// ---------------------------------------------------------------------------

/// WASM profiler session — manages one frame at a time.
///
/// Call [`Self::begin_frame`], open scopes with [`Self::scope`], then finalise with
/// [`Self::end_frame`] to get a [`WasmFrameReport`].
#[wasm_bindgen]
pub struct WasmProfilerSession {
    arena: Vec<ArenaNode>,
    stack: Vec<usize>,
    /// Parallel stack of scope start clock samples (milliseconds from
    /// [`now_ms`], `None` when the target has no clock) used by the JS-facing
    /// `push_scope_js` / `pop_scope_js` API.  The Rust-side RAII `scope()`
    /// API does not push to this stack.
    scope_starts: Vec<Option<f64>>,
    /// Clock sample taken by [`Self::begin_frame`], in milliseconds.
    frame_start: Option<f64>,
    /// Set when any clock sample of the current frame was unavailable; it
    /// becomes `!clock_available` in the next [`WasmFrameReport`].
    clock_missed: bool,
}

impl WasmProfilerSession {
    /// Create a new profiler session.
    pub fn new() -> Self {
        WasmProfilerSession {
            arena: Vec::new(),
            stack: Vec::new(),
            scope_starts: Vec::new(),
            frame_start: None,
            clock_missed: false,
        }
    }

    /// Take one clock sample, remembering when the target has no clock.
    fn sample_clock(&mut self) -> Option<f64> {
        let now = now_ms();
        if now.is_none() {
            self.clock_missed = true;
        }
        now
    }

    /// Nanoseconds elapsed since the clock sample `start`, or `0` when the
    /// interval could not be measured (which flags the frame as clockless).
    fn elapsed_since(&mut self, start: Option<f64>) -> u64 {
        let end = self.sample_clock();
        match elapsed_ns(start, end) {
            Some(ns) => ns,
            None => {
                self.clock_missed = true;
                0
            }
        }
    }

    /// Fold one finished call of arena node `node_idx` into its statistics.
    fn record_scope(&mut self, node_idx: usize, elapsed_ns: u64) {
        if let Some(node) = self.arena.get_mut(node_idx) {
            node.total_ns = node.total_ns.saturating_add(elapsed_ns);
            node.max_ns = node.max_ns.max(elapsed_ns);
            node.call_count += 1;
        }
    }

    /// Begin a new profiling frame (clears previous frame data).
    pub fn begin_frame(&mut self) {
        self.arena.clear();
        self.stack.clear();
        self.scope_starts.clear();
        self.clock_missed = false;
        self.frame_start = self.sample_clock();
    }

    /// Open a named scope.  Returns a guard that records time on drop.
    pub fn scope(&mut self, name: &str) -> WasmScopeGuard<'_> {
        let parent = self.stack.last().copied();
        // Find or create a node with this name under the current parent
        let node_idx = self
            .arena
            .iter()
            .position(|n| n.name == name && n.parent == parent)
            .unwrap_or_else(|| {
                let idx = self.arena.len();
                self.arena.push(ArenaNode {
                    name: name.to_string(),
                    parent,
                    call_count: 0,
                    total_ns: 0,
                    max_ns: 0,
                });
                idx
            });
        self.stack.push(node_idx);
        let start = self.sample_clock();
        WasmScopeGuard {
            session: self,
            start,
            node_idx,
        }
    }

    /// Finalise the current frame and return a [`WasmFrameReport`].
    ///
    /// Any scope left open (a `push_scope` without its `pop_scope`) is
    /// discarded together with the frame, so the next frame starts from an
    /// empty scope stack.
    pub fn end_frame(&mut self) -> WasmFrameReport {
        let frame_ns = match self.frame_start.take() {
            Some(start) => self.elapsed_since(Some(start)),
            None => 0,
        };
        let clock_available = !self.clock_missed;
        self.clock_missed = false;

        // Build scope tree from arena
        let root_children: Vec<WasmScopeNode> = self
            .arena
            .iter()
            .enumerate()
            .filter(|(_, n)| n.parent.is_none())
            .map(|(i, _)| self.build_node(i))
            .collect();

        let root = WasmScopeNode {
            name: "frame".to_string(),
            call_count: 1,
            total_ns: frame_ns,
            max_ns: frame_ns,
            children: root_children,
        };

        self.arena.clear();
        self.stack.clear();
        self.scope_starts.clear();

        WasmFrameReport {
            frame_ns,
            clock_available,
            root,
        }
    }

    fn build_node(&self, idx: usize) -> WasmScopeNode {
        let node = &self.arena[idx];
        let children: Vec<WasmScopeNode> = self
            .arena
            .iter()
            .enumerate()
            .filter(|(_, n)| n.parent == Some(idx))
            .map(|(i, _)| self.build_node(i))
            .collect();
        WasmScopeNode {
            name: node.name.clone(),
            call_count: node.call_count,
            total_ns: node.total_ns,
            max_ns: node.max_ns,
            children,
        }
    }

    /// Number of arena nodes (useful for testing that `end_frame` clears state).
    pub fn node_count(&self) -> usize {
        self.arena.len()
    }
}

impl Default for WasmProfilerSession {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// wasm-bindgen JavaScript API
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmProfilerSession {
    /// Create a new profiler session (JS constructor).
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> WasmProfilerSession {
        WasmProfilerSession::new()
    }

    /// Begin a new profiling frame (clears previous frame data).
    #[wasm_bindgen(js_name = "begin_frame")]
    pub fn begin_frame_js(&mut self) {
        self.begin_frame();
    }

    /// Open a named scope.
    ///
    /// The Rust RAII [`Self::scope`] API is not exposed because lifetime-bearing
    /// guards are incompatible with the WASM ABI.  Instead, JavaScript must
    /// pair every `push_scope` with a matching `pop_scope`.
    #[wasm_bindgen(js_name = "push_scope")]
    pub fn push_scope_js(&mut self, name: &str) {
        let parent = self.stack.last().copied();
        let node_idx = self
            .arena
            .iter()
            .position(|n| n.name == name && n.parent == parent)
            .unwrap_or_else(|| {
                let idx = self.arena.len();
                self.arena.push(ArenaNode {
                    name: name.to_string(),
                    parent,
                    call_count: 0,
                    total_ns: 0,
                    max_ns: 0,
                });
                idx
            });
        self.stack.push(node_idx);
        let start = self.sample_clock();
        self.scope_starts.push(start);
    }

    /// Close the most recently opened scope.
    ///
    /// No-op when the scope stack is empty.
    #[wasm_bindgen(js_name = "pop_scope")]
    pub fn pop_scope_js(&mut self) {
        let (Some(node_idx), Some(start)) = (self.stack.pop(), self.scope_starts.pop()) else {
            return;
        };
        let elapsed_ns = self.elapsed_since(start);
        self.record_scope(node_idx, elapsed_ns);
    }

    /// Finalise the current frame and return the report as a `JsValue`.
    ///
    /// # Errors
    ///
    /// Returns a `JsValue` error string if serialisation fails.
    #[wasm_bindgen(js_name = "end_frame")]
    pub fn end_frame_js(&mut self) -> std::result::Result<JsValue, JsValue> {
        let report = self.end_frame();
        to_js_value(&report)
    }

    /// Finalise the current frame and return the report as a JSON string.
    #[wasm_bindgen(js_name = "end_frame_json")]
    pub fn end_frame_json_js(&mut self) -> String {
        self.end_frame().to_json()
    }

    /// Finalise the current frame and return the report as CSV.
    #[wasm_bindgen(js_name = "end_frame_csv")]
    pub fn end_frame_csv_js(&mut self) -> String {
        self.end_frame().to_csv()
    }

    /// Finalise the current frame and return the report as folded-stacks
    /// (flamegraph-compatible) text.
    #[wasm_bindgen(js_name = "end_frame_folded_stacks")]
    pub fn end_frame_folded_js(&mut self) -> String {
        self.end_frame().to_folded_stacks()
    }

    /// Number of arena nodes currently allocated.
    #[wasm_bindgen(js_name = "node_count")]
    pub fn node_count_js(&self) -> u32 {
        self.node_count() as u32
    }

    /// Current scope-stack depth (useful for debugging mismatched push/pop).
    #[wasm_bindgen(js_name = "scope_depth")]
    pub fn scope_depth_js(&self) -> u32 {
        self.stack.len() as u32
    }
}

// ---------------------------------------------------------------------------
// JSON helpers used by the wasm-bindgen API
// ---------------------------------------------------------------------------

/// Parse a [`WasmFrameReport`] from its JSON form, returning a `JsValue`
/// error on failure.
///
/// This module-level helper is exposed for JS callers that wish to round-trip
/// a frame report through JSON (for example after IPC).
#[wasm_bindgen(js_name = "parse_wasm_frame_report")]
pub fn parse_wasm_frame_report(json: &str) -> std::result::Result<JsValue, JsValue> {
    let report: WasmFrameReport = serde_json::from_str(json).map_err(err_to_jsvalue)?;
    to_js_value(&report)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn one_frame<F: FnOnce(&mut WasmProfilerSession)>(f: F) -> WasmFrameReport {
        let mut s = WasmProfilerSession::new();
        s.begin_frame();
        f(&mut s);
        s.end_frame()
    }

    #[test]
    fn test_profiler_bridge_instantiation() {
        let report = one_frame(|_| {});
        let json = report.to_json();
        assert!(!json.is_empty(), "to_json should return non-empty string");
    }

    #[test]
    fn test_profiler_bridge_nested_scopes() {
        // Sequential sibling scopes (can't hold two &mut self borrows simultaneously).
        let report = one_frame(|s| {
            let _a = s.scope("a");
            // _a is dropped here (end of statement) before _b is created.
        });
        let report2 = one_frame(|s| {
            let _b = s.scope("b");
        });
        assert_eq!(report.root.name, "frame");
        assert_eq!(report.root.children.len(), 1);
        assert_eq!(report.root.children[0].name, "a");
        assert_eq!(report2.root.children.len(), 1);
        assert_eq!(report2.root.children[0].name, "b");
    }

    #[test]
    fn test_profiler_bridge_sibling_scopes() {
        let report = one_frame(|s| {
            {
                let _x = s.scope("x");
            }
            {
                let _y = s.scope("y");
            }
        });
        assert_eq!(report.root.children.len(), 2);
        let names: Vec<&str> = report
            .root
            .children
            .iter()
            .map(|c| c.name.as_str())
            .collect();
        assert!(names.contains(&"x"));
        assert!(names.contains(&"y"));
    }

    #[test]
    fn test_profiler_bridge_repeated_scope_accumulates() {
        let report = one_frame(|s| {
            for _ in 0..3 {
                let _z = s.scope("z");
            }
        });
        assert_eq!(report.root.children.len(), 1);
        let z = &report.root.children[0];
        assert_eq!(z.call_count, 3);
    }

    #[test]
    fn test_profiler_bridge_csv_format() {
        let report = one_frame(|s| {
            let _a = s.scope("work");
        });
        let csv = report.to_csv();
        let mut lines = csv.lines();
        let header = lines.next().expect("CSV should have header");
        assert_eq!(header, "name,call_count,total_ns,max_ns");
    }

    #[test]
    fn test_profiler_bridge_folded_stacks() {
        let report = one_frame(|s| {
            let _a = s.scope("physics");
        });
        let folded = report.to_folded_stacks();
        assert!(!folded.is_empty());
        assert!(folded.contains("frame"));
    }

    #[test]
    fn test_profiler_bridge_begin_clears_previous() {
        let mut s = WasmProfilerSession::new();
        s.begin_frame();
        {
            let _g = s.scope("scope1");
        }
        let _ = s.end_frame();
        s.begin_frame();
        assert_eq!(s.node_count(), 0, "begin_frame should clear arena");
    }

    // ---- wasm-safe clock regression (F5717 / F5719) -----------------------

    /// Busy-wait until `ms` milliseconds have elapsed on the profiler's clock.
    fn spin_for_ms(ms: f64) {
        let start = now_ms().expect("native clock is always available");
        while now_ms().expect("native clock is always available") - start < ms {
            std::hint::spin_loop();
        }
    }

    #[test]
    fn test_elapsed_ns_converts_milliseconds() {
        assert_eq!(elapsed_ns(Some(1.0), Some(3.5)), Some(2_500_000));
        assert_eq!(elapsed_ns(Some(10.0), Some(10.0)), Some(0));
    }

    #[test]
    fn test_elapsed_ns_without_clock_is_unmeasurable() {
        assert_eq!(elapsed_ns(None, Some(1.0)), None);
        assert_eq!(elapsed_ns(Some(1.0), None), None);
        assert_eq!(elapsed_ns(None, None), None);
        assert_eq!(elapsed_ns(Some(f64::NAN), Some(1.0)), None);
        assert_eq!(elapsed_ns(Some(0.0), Some(f64::INFINITY)), None);
    }

    #[test]
    fn test_elapsed_ns_backwards_clock_clamps_to_zero() {
        assert_eq!(elapsed_ns(Some(5.0), Some(4.0)), Some(0));
    }

    #[test]
    fn test_profiler_raii_scope_measures_real_time() {
        let report = one_frame(|s| {
            let _guard = s.scope("work");
            spin_for_ms(2.0);
        });
        assert!(report.clock_available);
        let work = &report.root.children[0];
        assert_eq!(work.call_count, 1);
        assert!(
            work.total_ns >= 2_000_000,
            "a 2 ms scope was recorded as {} ns",
            work.total_ns
        );
        assert_eq!(work.max_ns, work.total_ns);
        assert!(
            report.frame_ns >= work.total_ns,
            "frame {} ns shorter than its scope {} ns",
            report.frame_ns,
            work.total_ns
        );
        assert_eq!(report.root.total_ns, report.frame_ns);
    }

    #[test]
    fn test_profiler_js_push_pop_measures_real_time() {
        // The JS-facing begin_frame / push_scope / pop_scope / end_frame path
        // is the one that trapped on wasm32 when it used `Instant::now`.
        let mut s = WasmProfilerSession::new_js();
        s.begin_frame_js();
        s.push_scope_js("physics");
        assert_eq!(s.scope_depth_js(), 1);
        spin_for_ms(1.0);
        s.push_scope_js("narrowphase");
        spin_for_ms(1.0);
        s.pop_scope_js();
        s.pop_scope_js();
        assert_eq!(s.scope_depth_js(), 0);
        let report = s.end_frame();
        assert!(report.clock_available);
        let physics = &report.root.children[0];
        assert_eq!(physics.name, "physics");
        assert_eq!(physics.call_count, 1);
        assert!(physics.total_ns >= 2_000_000, "got {}", physics.total_ns);
        let narrow = &physics.children[0];
        assert_eq!(narrow.name, "narrowphase");
        assert!(narrow.total_ns >= 1_000_000, "got {}", narrow.total_ns);
        assert!(narrow.total_ns <= physics.total_ns);
    }

    #[test]
    fn test_profiler_end_frame_discards_unbalanced_scopes() {
        let mut s = WasmProfilerSession::new();
        s.begin_frame();
        s.push_scope_js("leaked");
        assert_eq!(s.scope_starts.len(), 1);
        let _ = s.end_frame();
        assert_eq!(s.scope_depth_js(), 0, "scope stack survived end_frame");
        assert!(
            s.scope_starts.is_empty(),
            "stale scope start times survived end_frame"
        );
        // A stray pop after the frame ended must be a harmless no-op.
        s.pop_scope_js();
        assert_eq!(s.scope_depth_js(), 0);
        // The next frame starts clean and times its own scope.
        s.begin_frame();
        s.push_scope_js("fresh");
        s.pop_scope_js();
        let report = s.end_frame();
        assert_eq!(report.root.children.len(), 1);
        assert_eq!(report.root.children[0].name, "fresh");
        assert_eq!(report.root.children[0].call_count, 1);
    }

    #[test]
    fn test_profiler_missing_clock_is_reported_not_fabricated() {
        let mut s = WasmProfilerSession::new();
        s.begin_frame();
        s.push_scope_js("x");
        // Emulate a clockless target: `now_ms()` returned `None` at push time.
        s.scope_starts[0] = None;
        spin_for_ms(1.0);
        s.pop_scope_js();
        let report = s.end_frame();
        assert!(
            !report.clock_available,
            "an unmeasured scope must clear clock_available"
        );
        let x = &report.root.children[0];
        assert_eq!(x.call_count, 1, "the call is still counted");
        assert_eq!(x.total_ns, 0, "unmeasured time must be 0, not invented");
        assert_eq!(x.max_ns, 0);

        // The flag is per frame: the next, fully measured frame is available.
        let report = one_frame(|s| {
            let _g = s.scope("y");
        });
        assert!(report.clock_available);
    }

    #[test]
    fn test_profiler_end_frame_without_begin_reports_zero_frame() {
        let mut s = WasmProfilerSession::new();
        let report = s.end_frame();
        assert_eq!(report.frame_ns, 0);
        assert_eq!(report.root.total_ns, 0);
        assert!(report.root.children.is_empty());
        assert!(report.clock_available);
    }

    #[test]
    fn test_frame_report_json_clock_available_roundtrip_and_default() {
        let report = one_frame(|s| {
            let _g = s.scope("a");
        });
        let json = report.to_json();
        assert!(json.contains("\"clock_available\":true"), "json = {json}");
        let back: WasmFrameReport = serde_json::from_str(&json).expect("round trip");
        assert!(back.clock_available);
        assert_eq!(back.root.children.len(), 1);

        // Reports serialised before the field existed still parse, and are
        // treated as measured with a real clock.
        let legacy = r#"{"frame_ns":42,"root":{"name":"frame","call_count":1,
            "total_ns":42,"max_ns":42,"children":[]}}"#;
        let old: WasmFrameReport = serde_json::from_str(legacy).expect("legacy json");
        assert_eq!(old.frame_ns, 42);
        assert!(old.clock_available);

        let unmeasured = r#"{"frame_ns":0,"clock_available":false,"root":{"name":"frame",
            "call_count":1,"total_ns":0,"max_ns":0,"children":[]}}"#;
        let flagged: WasmFrameReport = serde_json::from_str(unmeasured).expect("flagged json");
        assert!(!flagged.clock_available);
    }
}
