// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly debug and inspection tools.
//!
//! Provides debug draw lists, a physics inspector, a performance HUD,
//! runtime assertions, a replay controller, a scene explorer, a logger,
//! and a test harness — all designed for easy JavaScript interop. Each
//! sub-module exposes a focused set of types annotated for `wasm-bindgen`.

use wasm_bindgen::prelude::*;

pub mod assertions;
pub mod colors;
pub mod config;
pub mod draw;
pub mod harness;
pub mod inspector;
pub mod logger;
pub mod performance;
pub mod replay;
pub mod scene;

pub use assertions::{AssertResult, WasmPhysicsAssert};
pub use colors::Rgba;
pub use config::WasmDebugConfig;
pub use draw::{WasmDebugDraw, WasmDebugDrawList, WasmDrawCall};
pub use harness::{ExpectedBodyState, TestResult, WasmTestHarness};
pub use inspector::{BodyInspectResult, JointInspectResult, WasmPhysicsInspector};
pub use logger::{LogEntry, LogLevel, WasmPhysicsLogger};
pub use performance::{FrameTiming, WasmPerformanceHud};
pub use replay::{ReplayFrame, WasmReplayController};
pub use scene::{SceneNode, WasmSceneExplorer};

/// Convert a JavaScript `Float64Array` slice into a `[f64; 3]`, returning a
/// descriptive `JsValue` error if the slice is not exactly three components.
pub(crate) fn vec3_from_slice(slice: &[f64], name: &str) -> Result<[f64; 3], JsValue> {
    if slice.len() != 3 {
        return Err(JsValue::from_str(&format!(
            "{name} must have exactly 3 components, got {}",
            slice.len()
        )));
    }
    Ok([slice[0], slice[1], slice[2]])
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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
        let json = list.to_json().expect("to_json must succeed");
        assert!(json.contains("Line"));
    }

    // --- WasmPhysicsInspector ---

    #[test]
    fn test_inspector_body() {
        let mut insp = WasmPhysicsInspector::new();
        insp.register_body(1, [1.0, 2.0, 3.0], [0.1, 0.0, 0.0], [0.0; 3], false, false);
        let json = insp.inspect_body(1).expect("body 1 should exist");
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
        let json = insp.inspect_joint(0).expect("joint 0 should exist");
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
        let f = r.seek_to_frame(1).expect("frame 1 should exist");
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
        let json = e.get_hierarchy().expect("hierarchy JSON must serialise");
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
        let result = h.last_result.clone().expect("last_result should be set");
        assert!(!result.passed);
    }

    #[test]
    fn test_harness_baseline() {
        let mut h = WasmTestHarness::new();
        h.save_baseline("t1", "{\"x\":0}");
        let diff = h
            .diff_against_baseline("t1", "{\"x\":0}")
            .expect("baseline t1 must exist");
        assert_eq!(diff, "MATCH");
    }
}
