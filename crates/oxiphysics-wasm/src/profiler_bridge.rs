// Copyright (C) 2026 COOLJAPAN OU (Team Kitasan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly bridge for the hierarchical scoped profiler.
//!
//! Exposes a JSON-oriented surface around frame-level profiling with named
//! scopes, call counts, and nanosecond timing, suitable for use across the
//! WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use serde::{Deserialize, Serialize};
use std::time::Instant;
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

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
    /// Total wall-clock nanoseconds for the frame.
    pub frame_ns: u64,
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
    start: Instant,
    node_idx: usize,
}

impl Drop for WasmScopeGuard<'_> {
    fn drop(&mut self) {
        let elapsed_ns = self.start.elapsed().as_nanos() as u64;
        if let Some(node) = self.session.arena.get_mut(self.node_idx) {
            node.total_ns += elapsed_ns;
            if elapsed_ns > node.max_ns {
                node.max_ns = elapsed_ns;
            }
            node.call_count += 1;
        }
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
    /// Parallel stack of scope start `Instant`s used by the JS-facing
    /// `push_scope_js` / `pop_scope_js` API.  The Rust-side RAII `scope()`
    /// API does not push to this stack.
    scope_starts: Vec<Instant>,
    frame_start: Option<Instant>,
}

impl WasmProfilerSession {
    /// Create a new profiler session.
    pub fn new() -> Self {
        WasmProfilerSession {
            arena: Vec::new(),
            stack: Vec::new(),
            scope_starts: Vec::new(),
            frame_start: None,
        }
    }

    /// Begin a new profiling frame (clears previous frame data).
    pub fn begin_frame(&mut self) {
        self.arena.clear();
        self.stack.clear();
        self.scope_starts.clear();
        self.frame_start = Some(Instant::now());
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
        WasmScopeGuard {
            session: self,
            start: Instant::now(),
            node_idx,
        }
    }

    /// Finalise the current frame and return a [`WasmFrameReport`].
    pub fn end_frame(&mut self) -> WasmFrameReport {
        let frame_ns = self
            .frame_start
            .take()
            .map_or(0, |s| s.elapsed().as_nanos() as u64);

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

        WasmFrameReport { frame_ns, root }
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
        self.scope_starts.push(Instant::now());
    }

    /// Close the most recently opened scope.
    ///
    /// No-op when the scope stack is empty.
    #[wasm_bindgen(js_name = "pop_scope")]
    pub fn pop_scope_js(&mut self) {
        let (Some(node_idx), Some(start)) = (self.stack.pop(), self.scope_starts.pop()) else {
            return;
        };
        let elapsed_ns = start.elapsed().as_nanos() as u64;
        if let Some(node) = self.arena.get_mut(node_idx) {
            node.total_ns += elapsed_ns;
            if elapsed_ns > node.max_ns {
                node.max_ns = elapsed_ns;
            }
            node.call_count += 1;
        }
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
}
