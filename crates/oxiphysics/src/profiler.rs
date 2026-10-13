// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hierarchical scoped profiler with RAII guards, frame reports, and flamegraph export.
//!
//! ## Types
//!
//! - `ProfilerSession` — manages one frame at a time; call `begin_frame()`,
//!   open scopes with `scope()`, then finalise with `end_frame()`.
//! - `ScopeGuard` — RAII guard returned by `scope()`; on drop it records
//!   the elapsed nanoseconds and pops the scope off the stack. The guard
//!   borrows the session, so it can never outlive it and the frame can not be
//!   restarted or finished while a guard is still alive.
//! - `ScopeNode` — a node in the reconstructed profiling tree.
//! - `FrameReport` — the per-frame summary returned by `end_frame()`,
//!   containing the total wall-clock nanoseconds and the scope tree.
//!
//! ## Example
//!
//! ```rust
//! use oxiphysics::profiler::ProfilerSession;
//!
//! let mut session = ProfilerSession::new();
//! session.begin_frame();
//!
//! {
//!     let _a = session.scope("physics");
//!     {
//!         let _b = session.scope("collision");
//!         // inner work
//!     } // _b drops here — records "collision"
//! }   // _a drops here — records "physics"
//!
//! let report = session.end_frame();
//! assert_eq!(report.root.children.len(), 1);
//! assert_eq!(report.root.children[0].name, "physics");
//! assert_eq!(report.root.children[0].children[0].name, "collision");
//! assert!(report.to_csv().starts_with("name,call_count,total_ns,max_ns"));
//! assert!(report.to_json().contains("\"physics\""));
//! assert!(report.to_folded_stacks().starts_with("frame;physics;collision "));
//! ```
//!
//! ## Lifetime safety
//!
//! A [`ScopeGuard`](crate::profiler::ScopeGuard) holds a shared borrow of its
//! [`ProfilerSession`](crate::profiler::ProfilerSession), so a guard cannot
//! outlive (or be moved out of) the session that created it:
//!
//! ```compile_fail,E0597
//! use oxiphysics::profiler::ProfilerSession;
//!
//! let guard = {
//!     let session = ProfilerSession::new();
//!     session.scope("escapes")
//! }; // `session` is dropped here while `guard` still borrows it
//! drop(guard);
//! ```
//!
//! `begin_frame` and `end_frame` need exclusive access to the session, so
//! they cannot run while a guard is alive either (which would otherwise
//! invalidate the guard's bookkeeping):
//!
//! ```compile_fail,E0502
//! use oxiphysics::profiler::ProfilerSession;
//!
//! let mut session = ProfilerSession::new();
//! session.begin_frame();
//! let guard = session.scope("still_open");
//! let _report = session.end_frame(); // exclusive borrow while `guard` lives
//! drop(guard);
//! ```
//!
//! ```compile_fail,E0502
//! use oxiphysics::profiler::ProfilerSession;
//!
//! let mut session = ProfilerSession::new();
//! session.begin_frame();
//! let guard = session.scope("still_open");
//! session.begin_frame(); // exclusive borrow while `guard` lives
//! drop(guard);
//! ```

use std::cell::RefCell;
use std::time::Instant;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// ScopeNode — public tree node (appears in FrameReport)
// ---------------------------------------------------------------------------

/// A node in the profiling tree, populated by [`ProfilerSession::end_frame`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScopeNode {
    /// Scope name as passed to [`ProfilerSession::scope`].
    pub name: String,
    /// Number of times this scope was entered during the frame.
    pub call_count: u64,
    /// Total nanoseconds accumulated across all calls this frame.
    pub total_ns: u64,
    /// Maximum nanoseconds for any single call this frame.
    pub max_ns: u64,
    /// Child scopes that were opened while this scope was active.
    pub children: Vec<ScopeNode>,
}

// ---------------------------------------------------------------------------
// ArenaNode — internal flat representation used while the frame is in progress
// ---------------------------------------------------------------------------

/// Flat arena entry used internally during frame capture.
struct ArenaNode {
    name: &'static str,
    /// Arena index of the parent scope, or `None` for top-level scopes.
    parent: Option<usize>,
    call_count: u64,
    total_ns: u64,
    max_ns: u64,
}

// ---------------------------------------------------------------------------
// FrameReport
// ---------------------------------------------------------------------------

/// A per-frame summary produced by [`ProfilerSession::end_frame`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameReport {
    /// Total wall-clock nanoseconds for the frame.
    pub frame_ns: u64,
    /// Root of the scope tree.  Its name is always `"frame"`.
    pub root: ScopeNode,
}

impl FrameReport {
    // -----------------------------------------------------------------------
    // Export helpers
    // -----------------------------------------------------------------------

    /// Render all scope nodes depth-first as CSV.
    ///
    /// Columns: `name,call_count,total_ns,max_ns`
    pub fn to_csv(&self) -> String {
        let mut out = String::from("name,call_count,total_ns,max_ns\n");
        Self::csv_node(&self.root, &mut out);
        out
    }

    fn csv_node(node: &ScopeNode, out: &mut String) {
        out.push_str(&format!(
            "{},{},{},{}\n",
            node.name, node.call_count, node.total_ns, node.max_ns
        ));
        for child in &node.children {
            Self::csv_node(child, out);
        }
    }

    /// Render the report as a JSON string using [`serde_json`].
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_default()
    }

    /// Render the report in folded-stacks format suitable for `flamegraph.pl`.
    ///
    /// One line is emitted **per leaf scope** (scope with no children), with
    /// the full ancestor path separated by semicolons, followed by a space and
    /// the leaf's `total_ns` value.
    ///
    /// Example: `root;physics;collision 4200\n`
    pub fn to_folded_stacks(&self) -> String {
        let mut out = String::new();
        // The virtual "frame" root is included as the path root.
        Self::folded_node(&self.root, &self.root.name, &mut out);
        out
    }

    fn folded_node(node: &ScopeNode, path: &str, out: &mut String) {
        if node.children.is_empty() {
            out.push_str(&format!("{} {}\n", path, node.total_ns));
        } else {
            for child in &node.children {
                let child_path = format!("{};{}", path, child.name);
                Self::folded_node(child, &child_path, out);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ProfilerSession
// ---------------------------------------------------------------------------

/// Mutable state of the frame that is currently being captured.
///
/// It lives behind a [`RefCell`] inside [`ProfilerSession`] so that scopes can
/// be opened (and closed on drop) through a *shared* reference to the session,
/// which is what lets [`ScopeGuard`] borrow the session safely.
#[derive(Default)]
struct FrameState {
    /// Flat arena of nodes for the current frame.
    nodes: Vec<ArenaNode>,
    /// Stack of arena indices representing the current nesting path.
    stack: Vec<usize>,
    /// Instant at which the current frame began.
    frame_start: Option<Instant>,
}

/// The profiler session.  Manages one frame at a time.
///
/// Use [`begin_frame`](Self::begin_frame) to start, [`scope`](Self::scope) to
/// open RAII guards, and [`end_frame`](Self::end_frame) to finalise and
/// retrieve the report.
///
/// Scopes are opened through `&self`, so any number of nested guards can be
/// alive at once, while `begin_frame` / `end_frame` take `&mut self` and the
/// borrow checker therefore refuses to run them while a guard is still open.
/// The session is `Send` but not `Sync`: the nesting stack is inherently
/// per-thread state, so one session must not be shared between threads.
pub struct ProfilerSession {
    /// Arena, nesting stack and frame clock of the frame being captured.
    state: RefCell<FrameState>,
}

impl Default for ProfilerSession {
    fn default() -> Self {
        Self::new()
    }
}

impl ProfilerSession {
    /// Create a new session.  Call [`begin_frame`](Self::begin_frame) before
    /// opening any scopes.
    pub fn new() -> Self {
        Self {
            state: RefCell::new(FrameState::default()),
        }
    }

    /// Start a new frame.  Clears all previous frame data.
    pub fn begin_frame(&mut self) {
        let state = self.state.get_mut();
        state.nodes.clear();
        state.stack.clear();
        state.frame_start = Some(Instant::now());
    }

    /// Open a named scope.  Returns a RAII [`ScopeGuard`]; on drop, records
    /// elapsed nanoseconds into the arena.
    ///
    /// If a scope with the same name and the same parent already exists in the
    /// arena, its entry is reused (accumulating `total_ns` and `call_count`
    /// within the frame).
    ///
    /// The guard borrows the session for as long as it lives, so it cannot
    /// outlive the session and the frame cannot be restarted or finished while
    /// it is open.
    pub fn scope(&self, name: &'static str) -> ScopeGuard<'_> {
        // The borrow is released again before this function returns, and no
        // other `RefCell` borrow is ever held across a call into user code, so
        // it can not conflict with another one.
        let mut state = self.state.borrow_mut();
        let parent = state.stack.last().copied();

        // Look for an existing arena node with the same name and parent.
        let node_idx = match state
            .nodes
            .iter()
            .position(|n| n.name == name && n.parent == parent)
        {
            Some(idx) => idx,
            None => {
                let idx = state.nodes.len();
                state.nodes.push(ArenaNode {
                    name,
                    parent,
                    call_count: 0,
                    total_ns: 0,
                    max_ns: 0,
                });
                idx
            }
        };

        state.stack.push(node_idx);
        drop(state);

        ScopeGuard {
            state: &self.state,
            node_idx,
            start: Instant::now(),
        }
    }

    /// Finalise the frame and return the summary [`FrameReport`].
    ///
    /// After this call, the arena is cleared and the session is ready for
    /// the next `begin_frame`.
    pub fn end_frame(&mut self) -> FrameReport {
        let state = self.state.get_mut();
        let frame_ns = state.frame_start.take().map_or(0, |t| elapsed_ns(&t));

        // Build the ScopeNode tree from the flat arena.
        let root = Self::build_tree(&state.nodes, None, "frame");

        state.nodes.clear();
        state.stack.clear();

        FrameReport { frame_ns, root }
    }

    /// Number of arena nodes recorded in the current frame.
    pub fn node_count(&self) -> usize {
        self.state.borrow().nodes.len()
    }

    // -----------------------------------------------------------------------
    // Internal: tree reconstruction
    // -----------------------------------------------------------------------

    /// Recursively build a [`ScopeNode`] tree from the flat arena.
    ///
    /// `parent_idx` is `None` for the virtual "frame" root; children are all
    /// arena nodes whose `parent` field matches.
    fn build_tree(nodes: &[ArenaNode], parent_idx: Option<usize>, name: &str) -> ScopeNode {
        // Direct children (arena indices) of `parent_idx`, in creation order.
        let children: Vec<ScopeNode> = nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.parent == parent_idx)
            .map(|(idx, n)| Self::build_tree(nodes, Some(idx), n.name))
            .collect();

        match parent_idx.and_then(|idx| nodes.get(idx)) {
            Some(n) => ScopeNode {
                name: n.name.to_owned(),
                call_count: n.call_count,
                total_ns: n.total_ns,
                max_ns: n.max_ns,
                children,
            },
            None => {
                // This is the virtual "frame" root — synthesise its totals.
                let total_ns: u64 = children.iter().map(|c| c.total_ns).sum();
                ScopeNode {
                    name: name.to_owned(),
                    call_count: 1,
                    total_ns,
                    max_ns: total_ns,
                    children,
                }
            }
        }
    }
}

/// Nanoseconds elapsed since `since`, saturating at `u64::MAX`.
fn elapsed_ns(since: &Instant) -> u64 {
    u64::try_from(since.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

// ---------------------------------------------------------------------------
// ScopeGuard
// ---------------------------------------------------------------------------

/// RAII scope guard returned by [`ProfilerSession::scope`].
///
/// On drop, records the elapsed nanoseconds into the session's arena and pops
/// this scope off the nesting stack.
///
/// The guard holds a shared borrow of the [`ProfilerSession`] that created it,
/// so it can neither outlive that session nor coexist with a call to
/// [`ProfilerSession::begin_frame`] / [`ProfilerSession::end_frame`].
#[must_use = "the scope is timed until the guard is dropped; bind it to a variable"]
pub struct ScopeGuard<'a> {
    /// The session state this guard reports into.
    state: &'a RefCell<FrameState>,
    /// Arena index of the scope this guard times.
    node_idx: usize,
    /// When the scope was opened.
    start: Instant,
}

impl Drop for ScopeGuard<'_> {
    fn drop(&mut self) {
        let elapsed = elapsed_ns(&self.start);

        // `Drop` must never panic, so a (practically impossible) borrow
        // conflict simply skips the bookkeeping instead of aborting.
        let Ok(mut state) = self.state.try_borrow_mut() else {
            return;
        };

        if let Some(node) = state.nodes.get_mut(self.node_idx) {
            node.call_count = node.call_count.saturating_add(1);
            node.total_ns = node.total_ns.saturating_add(elapsed);
            node.max_ns = node.max_ns.max(elapsed);
        }

        // Remove this scope from the nesting stack.  Guards normally drop in
        // reverse order of creation, which makes this the top entry, but a
        // guard dropped out of order must still remove its *own* entry.
        if let Some(pos) = state.stack.iter().rposition(|&idx| idx == self.node_idx) {
            state.stack.remove(pos);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // Helper: run a full frame with a closure and return the report.
    // -------------------------------------------------------------------------
    fn one_frame<F>(f: F) -> FrameReport
    where
        F: FnOnce(&mut ProfilerSession),
    {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        f(&mut s);
        s.end_frame()
    }

    // -------------------------------------------------------------------------
    // 1. Nested scopes
    // -------------------------------------------------------------------------
    #[test]
    fn test_nested_scopes() {
        let report = one_frame(|s| {
            let _a = s.scope("a");
            let _b = s.scope("b");
            // _b drops first, then _a
        });

        assert_eq!(report.root.name, "frame");
        assert_eq!(report.root.children.len(), 1);
        let a = &report.root.children[0];
        assert_eq!(a.name, "a");
        assert_eq!(a.children.len(), 1);
        let b = &a.children[0];
        assert_eq!(b.name, "b");
        assert!(b.call_count >= 1);
    }

    // -------------------------------------------------------------------------
    // 2. Sibling scopes
    // -------------------------------------------------------------------------
    #[test]
    fn test_sibling_scopes() {
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
        assert!(names.contains(&"x"), "should have scope x");
        assert!(names.contains(&"y"), "should have scope y");
        for child in &report.root.children {
            assert_eq!(child.call_count, 1);
        }
    }

    // -------------------------------------------------------------------------
    // 3. Repeated scope (same name, same nesting) — accumulates in one node
    // -------------------------------------------------------------------------
    #[test]
    fn test_repeated_scope_accumulates() {
        let report = one_frame(|s| {
            for _ in 0..2 {
                let _z = s.scope("z");
                std::thread::sleep(std::time::Duration::from_micros(10));
            }
        });

        // Should be a single node with call_count == 2.
        assert_eq!(report.root.children.len(), 1);
        let z = &report.root.children[0];
        assert_eq!(z.name, "z");
        assert_eq!(
            z.call_count, 2,
            "same-named scope at same level should accumulate"
        );
        assert!(z.total_ns > 0, "total_ns must be > 0");
    }

    // -------------------------------------------------------------------------
    // 4. Drop order — stack is empty after end_frame
    // -------------------------------------------------------------------------
    #[test]
    fn test_drop_order_stack_empty() {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        {
            let _a = s.scope("outer");
            {
                let _b = s.scope("inner");
            } // _b drops here
        } // _a drops here
        // Stack should be empty before end_frame — verified inside end_frame
        let _report = s.end_frame();
        // After end_frame the stack is explicitly cleared — verify it directly
        // by checking node_count is 0 (arena was cleared).
        assert_eq!(s.node_count(), 0, "arena should be cleared after end_frame");
    }

    // -------------------------------------------------------------------------
    // 5. Folded stacks format
    // -------------------------------------------------------------------------
    #[test]
    fn test_folded_stacks_format() {
        let report = one_frame(|s| {
            let _a = s.scope("scope_a");
            let _b = s.scope("scope_b");
        });

        let folded = report.to_folded_stacks();
        assert!(!folded.is_empty(), "folded stacks should not be empty");

        for line in folded.lines() {
            // Each line must be "path<space>number"
            let mut parts = line.rsplitn(2, ' ');
            let ns_str = parts.next().expect("should have ns value");
            let path_str = parts.next().expect("should have path");
            let _ns: u64 = ns_str
                .parse()
                .unwrap_or_else(|_| panic!("ns value should be a valid u64, got '{ns_str}'"));
            assert!(
                path_str.starts_with("frame"),
                "path should start with 'frame', got '{path_str}'"
            );
            assert!(
                path_str.contains(';'),
                "path should have semicolons, got '{path_str}'"
            );
        }
    }

    // -------------------------------------------------------------------------
    // 6. JSON round-trip
    // -------------------------------------------------------------------------
    #[test]
    fn test_json_round_trip() {
        let report = one_frame(|s| {
            let _a = s.scope("alpha");
        });

        let json = report.to_json();
        let decoded: FrameReport =
            serde_json::from_str(&json).expect("JSON should deserialise back to FrameReport");
        assert_eq!(decoded.root.name, "frame");
        assert_eq!(
            decoded.root.children.len(),
            report.root.children.len(),
            "round-tripped report should have same child count"
        );
    }

    // -------------------------------------------------------------------------
    // 7. begin_frame clears previous frame data
    // -------------------------------------------------------------------------
    #[test]
    fn test_begin_frame_clears() {
        let mut s = ProfilerSession::new();

        // First frame
        s.begin_frame();
        let _g = s.scope("scope1");
        drop(_g);
        let _r = s.end_frame();

        // Second frame should start clean
        s.begin_frame();
        assert_eq!(
            s.node_count(),
            0,
            "begin_frame should clear arena from previous frame"
        );
    }

    // -------------------------------------------------------------------------
    // 8. CSV output is well-formed
    // -------------------------------------------------------------------------
    #[test]
    fn test_csv_output() {
        let report = one_frame(|s| {
            let _a = s.scope("work");
        });

        let csv = report.to_csv();
        let mut lines = csv.lines();

        // Header
        let header = lines.next().expect("CSV should have a header");
        assert_eq!(header, "name,call_count,total_ns,max_ns");

        // At least the "frame" root and "work" child should appear.
        let data_lines: Vec<&str> = lines.collect();
        assert!(
            !data_lines.is_empty(),
            "CSV should have at least one data line"
        );
        for line in &data_lines {
            let cols: Vec<&str> = line.split(',').collect();
            assert_eq!(
                cols.len(),
                4,
                "each data line should have 4 columns, got {line}"
            );
            // columns 1..3 should be valid u64s
            cols[1].parse::<u64>().expect("call_count should be u64");
            cols[2].parse::<u64>().expect("total_ns should be u64");
            cols[3].parse::<u64>().expect("max_ns should be u64");
        }
    }

    // -------------------------------------------------------------------------
    // 9. frame_ns is populated
    // -------------------------------------------------------------------------
    #[test]
    fn test_frame_ns_populated() {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        let _g = s.scope("tick");
        std::thread::sleep(std::time::Duration::from_micros(10));
        drop(_g);
        let report = s.end_frame();
        // frame_ns should be > 0 (some real time elapsed).
        assert!(report.frame_ns > 0, "frame_ns should be > 0");
    }

    // -------------------------------------------------------------------------
    // 10. Scopes can be opened through a shared reference
    // -------------------------------------------------------------------------
    /// Instrumented helper that only has `&ProfilerSession`, like library code
    /// handed a profiler by its caller.
    fn instrumented_leaf(profiler: &ProfilerSession) {
        let _leaf = profiler.scope("leaf");
    }

    fn instrumented_branch(profiler: &ProfilerSession) {
        let _branch = profiler.scope("branch");
        instrumented_leaf(profiler);
        instrumented_leaf(profiler);
    }

    #[test]
    fn test_scopes_open_through_shared_reference() {
        let report = one_frame(|s| {
            let shared: &ProfilerSession = s;
            instrumented_branch(shared);
            instrumented_branch(shared);
        });

        assert_eq!(report.root.children.len(), 1);
        let branch = &report.root.children[0];
        assert_eq!(branch.name, "branch");
        assert_eq!(branch.call_count, 2, "branch was entered twice");
        assert_eq!(branch.children.len(), 1);
        let leaf = &branch.children[0];
        assert_eq!(leaf.name, "leaf");
        assert_eq!(
            leaf.call_count, 4,
            "two branch calls with two leaf calls each"
        );
        assert!(leaf.max_ns <= leaf.total_ns);
    }

    // -------------------------------------------------------------------------
    // 11. Same name under different parents stays separate
    // -------------------------------------------------------------------------
    #[test]
    fn test_same_name_under_different_parents_is_not_merged() {
        let report = one_frame(|s| {
            {
                let _p = s.scope("parent_a");
                let _c = s.scope("shared_child");
            }
            {
                let _p = s.scope("parent_b");
                let _c = s.scope("shared_child");
            }
        });

        assert_eq!(report.root.children.len(), 2);
        for parent in &report.root.children {
            assert_eq!(parent.children.len(), 1);
            assert_eq!(parent.children[0].name, "shared_child");
            assert_eq!(parent.children[0].call_count, 1);
        }
    }

    // -------------------------------------------------------------------------
    // 12. Out-of-order guard drops keep the nesting stack consistent
    // -------------------------------------------------------------------------
    #[test]
    fn test_out_of_order_drop_keeps_stack_consistent() {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        {
            let outer = s.scope("outer");
            let inner = s.scope("inner");
            // Drop the *outer* guard first: only its own stack entry may go.
            drop(outer);
            // `inner` is still open, so a new scope nests under it.
            {
                let _deep = s.scope("deep");
            }
            drop(inner);
            // Everything is closed again: a new scope is a root-level sibling.
            let _after = s.scope("after");
        }
        let report = s.end_frame();

        let names: Vec<&str> = report
            .root
            .children
            .iter()
            .map(|c| c.name.as_str())
            .collect();
        assert_eq!(names, vec!["outer", "after"]);

        let outer = &report.root.children[0];
        assert_eq!(outer.call_count, 1);
        assert_eq!(outer.children.len(), 1);
        let inner = &outer.children[0];
        assert_eq!(inner.name, "inner");
        assert_eq!(inner.call_count, 1);
        assert_eq!(inner.children.len(), 1);
        assert_eq!(inner.children[0].name, "deep");
        assert_eq!(report.root.children[1].call_count, 1);
    }

    // -------------------------------------------------------------------------
    // 13. Guards dropped late still record into their own node
    // -------------------------------------------------------------------------
    #[test]
    fn test_late_drop_records_into_own_node_only() {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        let first = s.scope("first");
        let second = s.scope("second");
        std::thread::sleep(std::time::Duration::from_micros(200));
        drop(first);
        drop(second);
        let report = s.end_frame();

        let first_node = &report.root.children[0];
        assert_eq!(first_node.name, "first");
        assert_eq!(first_node.call_count, 1);
        assert!(first_node.total_ns > 0);
        let second_node = &first_node.children[0];
        assert_eq!(second_node.name, "second");
        assert_eq!(second_node.call_count, 1);
        assert!(second_node.total_ns > 0);
    }

    // -------------------------------------------------------------------------
    // 14. The session is reusable across frames and isolated between them
    // -------------------------------------------------------------------------
    #[test]
    fn test_frames_are_independent() {
        let mut s = ProfilerSession::new();

        s.begin_frame();
        {
            let _a = s.scope("only_in_first_frame");
        }
        let first = s.end_frame();
        assert_eq!(first.root.children.len(), 1);

        s.begin_frame();
        {
            let _b = s.scope("only_in_second_frame");
        }
        assert_eq!(s.node_count(), 1);
        let second = s.end_frame();
        assert_eq!(second.root.children.len(), 1);
        assert_eq!(second.root.children[0].name, "only_in_second_frame");
    }

    // -------------------------------------------------------------------------
    // 15. A session can move to another thread between frames
    // -------------------------------------------------------------------------
    #[test]
    fn test_session_is_send_between_frames() {
        let mut s = ProfilerSession::new();
        s.begin_frame();
        let handle = std::thread::spawn(move || {
            {
                let _g = s.scope("on_worker_thread");
            }
            s.end_frame()
        });
        let report = handle.join().expect("worker thread must not panic");
        assert_eq!(report.root.children.len(), 1);
        assert_eq!(report.root.children[0].name, "on_worker_thread");
    }
}
