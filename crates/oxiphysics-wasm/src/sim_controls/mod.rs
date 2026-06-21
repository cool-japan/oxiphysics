// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly simulation control APIs.
//!
//! This module provides pure-Rust structs for controlling physics simulations,
//! recording replay data, profiling performance, and exporting simulation state.
//! All types are designed for ergonomic use from JavaScript via wasm-bindgen.

pub mod config;
pub mod controller;
pub mod debug_draw;
pub mod exporter;
pub mod integrator;
pub mod profiler;
pub mod replay;
pub mod state;

// Re-export BroadphaseType / SolverType from physics_config (F1 canonical).
pub use crate::physics_config::{BroadphaseType, SolverType};

// Flatten public API — all types stay accessible as `crate::sim_controls::Foo`.
pub use config::SimulationConfig;
pub use controller::SimulationController;
pub use debug_draw::{Color4, DebugDrawConfig};
pub use exporter::{ExportFormat, SimulationExporter, SimulationSnapshot};
pub use integrator::{IntegratorType, TimeIntegrator, integrator_is_symplectic, integrator_name};
pub use profiler::{PhysicsProfiler, PhysicsProfilerSample};
pub use replay::{ReplayBuffer, ReplayFrame};
pub use state::{
    SimulationState, SimulationStats, StepResult, sim_state_is_active, sim_state_is_recording,
    sim_state_label,
};

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
        let json = serde_json::to_string(&cfg).expect("serialize SimulationConfig");
        let cfg2: SimulationConfig =
            serde_json::from_str(&json).expect("deserialize SimulationConfig");
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

    #[test]
    fn test_controller_perf_ms_is_real_measurement_not_input_scaling() {
        // The old fabrication set perf_ms = frame_time * 1000 * 0.01, i.e. a
        // deterministic 1% of the input. Verify the reported value is a real
        // wall-clock measurement instead: on the host it must be finite and
        // non-negative, and it must NOT equal that fabricated formula (the
        // chance a real elapsed time coincides exactly is negligible).
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        ctrl.resume();
        let frame_time = 1.0 / 60.0;
        let r = ctrl.step(frame_time);
        let fabricated = frame_time * 1_000.0 * 0.01;
        assert!(
            r.perf_ms.is_finite() && r.perf_ms >= 0.0,
            "perf_ms should be a real non-negative measurement on host, got {}",
            r.perf_ms
        );
        assert!(
            (r.perf_ms - fabricated).abs() > f64::EPSILON,
            "perf_ms still matches the discarded fabricated formula ({fabricated})"
        );
    }

    #[test]
    fn test_controller_perf_ms_zero_while_paused() {
        // A paused controller does no physics work; perf_ms must be exactly
        // zero, never a fabricated fraction of the frame time.
        let mut ctrl = SimulationController::new(SimulationConfig::default());
        let r = ctrl.step(1.0);
        assert_eq!(r.perf_ms, 0.0);
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
        buf2.import_frames(&json).expect("import_frames");
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
        let latest = prof.latest().expect("latest sample");
        assert!((latest.total_ms - 6.5).abs() < 1e-10);
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
        let snap2 = exp.import_state(&json).expect("import_state");
        assert!((snap2.time - 1.5).abs() < 1e-10);
        assert_eq!(snap2.positions.len(), 3);
    }

    #[test]
    fn test_exporter_binary_roundtrip_header() {
        let exp = SimulationExporter::new();
        let snap = SimulationSnapshot::new(2.72);
        let bytes = exp.export_binary(&snap);
        // First 8 bytes = time as little-endian f64
        let header: [u8; 8] = bytes[0..8].try_into().expect("8-byte header");
        let time = f64::from_le_bytes(header);
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
        let json = serde_json::to_string(&cfg).expect("serialize DebugDrawConfig");
        let cfg2: DebugDrawConfig =
            serde_json::from_str(&json).expect("deserialize DebugDrawConfig");
        assert_eq!(cfg.show_contacts, cfg2.show_contacts);
    }

    // --- BroadphaseType ---

    #[test]
    fn test_broadphase_default() {
        // The canonical default (from physics_config::BroadphaseType) is SweepAndPrune.
        assert_eq!(BroadphaseType::default(), BroadphaseType::SweepAndPrune);
    }

    // --- SimulationStats ---

    #[test]
    fn test_stats_advance() {
        let mut stats = SimulationStats::new();
        stats.advance(1.0 / 60.0);
        assert_eq!(stats.step_count, 1);
        assert!((stats.elapsed_time - 1.0 / 60.0).abs() < 1e-15);
    }

    // --- Free helper functions ---

    #[test]
    fn test_sim_state_helpers() {
        assert!(sim_state_is_active(SimulationState::Running));
        assert!(!sim_state_is_active(SimulationState::Paused));
        assert!(sim_state_is_recording(SimulationState::Recording));
        assert_eq!(sim_state_label(SimulationState::Stepping), "stepping");
    }

    #[test]
    fn test_integrator_helpers() {
        assert_eq!(integrator_name(IntegratorType::Euler), "euler");
        assert!(integrator_is_symplectic(IntegratorType::SemiImplicit));
        assert!(!integrator_is_symplectic(IntegratorType::Rk4));
    }
}
