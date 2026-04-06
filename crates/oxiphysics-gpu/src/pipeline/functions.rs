//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WorldState;

/// Mock broad-phase: returns number of potential collision pairs.
#[allow(dead_code)]
pub(super) fn run_broadphase(world: &WorldState) -> usize {
    let n = world.body_count();
    n.saturating_sub(1) * n / 2
}
/// Mock constraint solve: returns number of solved constraints (= collision pairs).
#[allow(dead_code)]
pub(super) fn run_constraint_solve(world: &WorldState) -> usize {
    let n = world.body_count();
    n.saturating_sub(1) * n / 2
}
/// Mock integration: semi-implicit Euler with gravity `[0, -9.81, 0]`.
#[allow(dead_code)]
pub(super) fn run_integration(world: &mut WorldState, dt: f64) {
    let n = world.body_count();
    for i in 0..n {
        let inv_m = world.inverse_masses[i];
        if inv_m == 0.0 {
            continue;
        }
        if world.velocities.len() > i * 3 + 1 {
            world.velocities[i * 3 + 1] -= 9.81 * dt;
        }
        if world.positions.len() >= (i + 1) * 3 && world.velocities.len() >= (i + 1) * 3 {
            world.positions[i * 3] += world.velocities[i * 3] * dt;
            world.positions[i * 3 + 1] += world.velocities[i * 3 + 1] * dt;
            world.positions[i * 3 + 2] += world.velocities[i * 3 + 2] * dt;
        }
    }
}
/// Mock post-process: no-op for an empty world; just validates invariants.
#[allow(dead_code)]
pub(super) fn run_postprocess(_world: &mut WorldState) {}
#[cfg(test)]
mod tests {

    use crate::pipeline::AsyncComputeQueue;

    use crate::pipeline::BarrierSet;
    use crate::pipeline::BufferUsage;

    use crate::pipeline::ComputePipeline;
    use crate::pipeline::CpuBuffer;
    use crate::pipeline::DispatchBatch;

    use crate::pipeline::GpuMemoryPool;

    use crate::pipeline::PhysicsPipeline;
    use crate::pipeline::PipelineBuilder;
    use crate::pipeline::PipelineConfig;
    use crate::pipeline::PipelineProfiler;
    use crate::pipeline::PipelineStage;

    use crate::pipeline::PipelineStats;

    use crate::pipeline::ResourceBarrier;
    use crate::pipeline::ResourceHandle;

    use crate::pipeline::WorldState;
    #[test]
    fn test_pipeline_workgroups() {
        let p = ComputePipeline::new("test", "", "main");
        assert_eq!(p.workgroups_needed(200), [4, 1, 1]);
        assert_eq!(p.workgroups_needed(128), [2, 1, 1]);
        assert_eq!(p.workgroups_needed(1), [1, 1, 1]);
    }
    #[test]
    fn test_cpu_buffer_zeros() {
        let buf = CpuBuffer::new_zeros("zeros", 10, BufferUsage::Storage);
        assert_eq!(buf.len(), 10);
        assert!(buf.data.iter().all(|&v| v == 0.0));
        assert!(!buf.is_empty());
    }
    #[test]
    fn test_dispatch_batch_bind() {
        let pipeline = ComputePipeline::new("batch", "", "main");
        let mut batch = DispatchBatch::new(pipeline, 64);
        batch.bind(CpuBuffer::new_zeros("a", 64, BufferUsage::Storage));
        batch.bind(CpuBuffer::new_zeros("b", 64, BufferUsage::StorageReadOnly));
        batch.bind(CpuBuffer::new_f32("c", vec![1.0; 4], BufferUsage::Uniform));
        assert_eq!(batch.bindings.len(), 3);
    }
    #[test]
    fn default_config_has_all_stages() {
        let cfg = PipelineConfig::new();
        for stage in PipelineStage::all_in_order() {
            assert!(cfg.is_enabled(stage), "stage {stage:?} should be enabled");
        }
        assert_eq!(cfg.substeps, 1);
        assert!(!cfg.use_gpu);
    }
    #[test]
    fn config_disable_stage() {
        let mut cfg = PipelineConfig::new();
        cfg.enabled_stages
            .retain(|&s| s != PipelineStage::PostProcess);
        assert!(!cfg.is_enabled(PipelineStage::PostProcess));
        assert!(cfg.is_enabled(PipelineStage::BroadPhase));
    }
    #[test]
    fn builder_pattern_substeps_and_gpu_flag() {
        let p = PipelineBuilder::new().substeps(4).use_gpu(true).build();
        assert_eq!(p.config.substeps, 4);
        assert!(p.config.use_gpu);
    }
    #[test]
    fn builder_disable_stage() {
        let p = PipelineBuilder::new()
            .disable_stage(PipelineStage::NarrowPhase)
            .build();
        assert!(!p.config.is_enabled(PipelineStage::NarrowPhase));
        assert!(p.config.is_enabled(PipelineStage::BroadPhase));
    }
    #[test]
    fn builder_enable_stage_idempotent() {
        let p = PipelineBuilder::new()
            .enable_stage(PipelineStage::BroadPhase)
            .enable_stage(PipelineStage::BroadPhase)
            .build();
        let count = p
            .config
            .enabled_stages
            .iter()
            .filter(|&&s| s == PipelineStage::BroadPhase)
            .count();
        assert_eq!(count, 1);
    }
    #[test]
    fn stage_order_is_canonical() {
        let order = PipelineStage::all_in_order();
        for w in order.windows(2) {
            assert!(
                w[0] < w[1],
                "stage order violated: {:?} should be < {:?}",
                w[0],
                w[1]
            );
        }
    }
    #[test]
    fn empty_world_step_no_panic() {
        let mut pipeline = PhysicsPipeline::new(PipelineConfig::new());
        let mut world = WorldState::default();
        let stats = pipeline.step(&mut world, 1.0 / 60.0);
        assert_eq!(stats.collision_pairs, 0);
        assert_eq!(stats.solved_constraints, 0);
    }
    #[test]
    fn stats_accumulate_sums_fields() {
        let mut a = PipelineStats {
            broadphase_ms: 1.0,
            narrowphase_ms: 2.0,
            constraint_ms: 3.0,
            integration_ms: 4.0,
            postprocess_ms: 5.0,
            total_time_ms: 6.0,
            collision_pairs: 10,
            solved_constraints: 5,
        };
        let b = PipelineStats {
            broadphase_ms: 0.5,
            narrowphase_ms: 0.5,
            constraint_ms: 0.5,
            integration_ms: 0.5,
            postprocess_ms: 0.5,
            total_time_ms: 1.0,
            collision_pairs: 3,
            solved_constraints: 2,
        };
        a.accumulate(&b);
        assert!((a.broadphase_ms - 1.5).abs() < 1e-12);
        assert_eq!(a.collision_pairs, 13);
        assert_eq!(a.solved_constraints, 7);
    }
    #[test]
    fn pipeline_cumulative_stats_grow() {
        let mut pipeline = PhysicsPipeline::new(PipelineConfig::new());
        let mut world = WorldState {
            positions: vec![0.0; 6],
            velocities: vec![0.0; 6],
            inverse_masses: vec![1.0, 1.0],
        };
        pipeline.step(&mut world, 0.016);
        pipeline.step(&mut world, 0.016);
        assert!(pipeline.stats.total_time_ms >= 0.0);
    }
    #[test]
    fn integration_stage_moves_body_downward() {
        let cfg = PipelineBuilder::new()
            .disable_stage(PipelineStage::BroadPhase)
            .disable_stage(PipelineStage::NarrowPhase)
            .disable_stage(PipelineStage::ConstraintSolve)
            .disable_stage(PipelineStage::PostProcess)
            .build()
            .config;
        let mut pipeline = PhysicsPipeline::new(cfg);
        let mut world = WorldState {
            positions: vec![0.0, 10.0, 0.0],
            velocities: vec![0.0, 0.0, 0.0],
            inverse_masses: vec![1.0],
        };
        pipeline.step(&mut world, 1.0);
        assert!(
            world.positions[1] < 10.0,
            "y should decrease due to gravity"
        );
    }
    #[test]
    fn barrier_valid_order() {
        let b = ResourceBarrier::new(
            PipelineStage::BroadPhase,
            PipelineStage::NarrowPhase,
            "pair_buffer",
        );
        assert!(b.is_valid_order());
    }
    #[test]
    fn barrier_invalid_reverse_order() {
        let b = ResourceBarrier::new(PipelineStage::NarrowPhase, PipelineStage::BroadPhase, "bad");
        assert!(!b.is_valid_order());
    }
    #[test]
    fn barrier_set_validate_catches_bad() {
        let mut bs = BarrierSet::new();
        bs.add(ResourceBarrier::new(
            PipelineStage::BroadPhase,
            PipelineStage::NarrowPhase,
            "ok",
        ));
        bs.add(ResourceBarrier::new(
            PipelineStage::PostProcess,
            PipelineStage::Integration,
            "backwards",
        ));
        assert_eq!(bs.validate().len(), 1);
    }
    #[test]
    fn barrier_set_filter_by_stage() {
        let mut bs = BarrierSet::new();
        bs.add(ResourceBarrier::new(
            PipelineStage::BroadPhase,
            PipelineStage::NarrowPhase,
            "pair_buf",
        ));
        bs.add(ResourceBarrier::new(
            PipelineStage::NarrowPhase,
            PipelineStage::ConstraintSolve,
            "contact_buf",
        ));
        assert_eq!(bs.barriers_from(PipelineStage::BroadPhase).len(), 1);
        assert_eq!(bs.barriers_to(PipelineStage::ConstraintSolve).len(), 1);
        assert_eq!(bs.len(), 2);
    }
    #[test]
    fn async_queue_submit_and_flush() {
        let mut q = AsyncComputeQueue::new();
        assert!(q.is_idle());
        let p = ComputePipeline::new("test", "", "main");
        q.submit(DispatchBatch::new(p.clone(), 64));
        q.submit(DispatchBatch::new(p, 128));
        assert_eq!(q.pending(), 2);
        let executed = q.flush();
        assert_eq!(executed, 2);
        assert!(q.is_idle());
        assert_eq!(q.total_enqueued, 2);
        assert_eq!(q.total_executed, 2);
    }
    #[test]
    fn async_queue_flush_empty_returns_zero() {
        let mut q = AsyncComputeQueue::new();
        assert_eq!(q.flush(), 0);
    }
    #[test]
    fn async_queue_pending_decrements_on_flush() {
        let mut q = AsyncComputeQueue::new();
        let p = ComputePipeline::new("t", "", "main");
        for _ in 0..5 {
            q.submit(DispatchBatch::new(p.clone(), 32));
        }
        assert_eq!(q.pending(), 5);
        q.flush();
        assert_eq!(q.pending(), 0);
    }
    #[test]
    fn profiler_record_and_summary() {
        let mut prof = PipelineProfiler::new();
        prof.record("broadphase", 1.0);
        prof.record("broadphase", 3.0);
        let (mean, _std, n) = prof.summary("broadphase").unwrap();
        assert!((mean - 2.0).abs() < 1e-10);
        assert_eq!(n, 2);
    }
    #[test]
    fn profiler_summary_unknown_stage_is_none() {
        let prof = PipelineProfiler::new();
        assert!(prof.summary("nonexistent").is_none());
    }
    #[test]
    fn profiler_total_samples() {
        let mut prof = PipelineProfiler::new();
        prof.record("a", 1.0);
        prof.record("a", 2.0);
        prof.record("b", 3.0);
        assert_eq!(prof.total_samples(), 3);
    }
    #[test]
    fn profiler_stage_names_sorted() {
        let mut prof = PipelineProfiler::new();
        prof.record("zzz", 1.0);
        prof.record("aaa", 1.0);
        prof.record("mmm", 1.0);
        let names = prof.stage_names();
        assert_eq!(names, vec!["aaa", "mmm", "zzz"]);
    }
    #[test]
    fn profiler_reset_clears_all() {
        let mut prof = PipelineProfiler::new();
        prof.record("broadphase", 1.5);
        prof.reset();
        assert_eq!(prof.total_samples(), 0);
        assert!(prof.stage_names().is_empty());
    }
    #[test]
    fn profiler_stddev_uniform_values() {
        let mut prof = PipelineProfiler::new();
        for _ in 0..4 {
            prof.record("stage", 5.0);
        }
        let (mean, std, _) = prof.summary("stage").unwrap();
        assert!((mean - 5.0).abs() < 1e-10);
        assert!(std < 1e-10, "stddev should be ~0 for identical samples");
    }
    #[test]
    fn pool_alloc_basic() {
        let mut pool = GpuMemoryPool::new(1024);
        let h = pool.alloc(256).expect("alloc should succeed");
        assert_eq!(h.1, 256);
        assert_eq!(pool.allocated, 256);
        assert_eq!(pool.free_space(), 768);
    }
    #[test]
    fn pool_alloc_and_free_round_trip() {
        let mut pool = GpuMemoryPool::new(1024);
        let h = pool.alloc(512).unwrap();
        pool.free(h.0, h.1).expect("free should succeed");
        assert!(pool.is_fully_free());
        assert_eq!(pool.fragmentation_count(), 1);
    }
    #[test]
    fn pool_alloc_exhaustion_returns_none() {
        let mut pool = GpuMemoryPool::new(100);
        let _ = pool.alloc(100).unwrap();
        assert!(
            pool.alloc(1).is_none(),
            "pool exhausted → alloc should fail"
        );
    }
    #[test]
    fn pool_multiple_allocs_and_frees() {
        let mut pool = GpuMemoryPool::new(1024);
        let h1 = pool.alloc(256).unwrap();
        let h2 = pool.alloc(256).unwrap();
        pool.free(h1.0, h1.1).unwrap();
        pool.free(h2.0, h2.1).unwrap();
        assert!(pool.is_fully_free());
    }
    #[test]
    fn pool_alloc_buffer_returns_correct_size() {
        let mut pool = GpuMemoryPool::new(2048);
        let (buf, handle) = pool
            .alloc_buffer("positions", 512, BufferUsage::Storage)
            .expect("alloc_buffer should succeed");
        assert_eq!(buf.len(), 512);
        assert_eq!(handle.1, 512);
    }
    #[test]
    fn resource_handle_from_alloc() {
        let handle = ResourceHandle::from_alloc((64, 128));
        assert_eq!(handle.offset, 64);
        assert_eq!(handle.size, 128);
    }
    #[test]
    fn pool_coalesces_free_blocks() {
        let mut pool = GpuMemoryPool::new(300);
        let h1 = pool.alloc(100).unwrap();
        let h2 = pool.alloc(100).unwrap();
        let h3 = pool.alloc(100).unwrap();
        pool.free(h1.0, h1.1).unwrap();
        pool.free(h2.0, h2.1).unwrap();
        pool.free(h3.0, h3.1).unwrap();
        assert_eq!(pool.fragmentation_count(), 1, "all blocks should coalesce");
    }
}
#[cfg(test)]
mod extended_pipeline_tests {

    use crate::pipeline::BarrierOptimizer;

    use crate::pipeline::ComputeOverlapScheduler;
    use crate::pipeline::ComputePipeline;

    use crate::pipeline::DispatchBatch;
    use crate::pipeline::FrameGraph;
    use crate::pipeline::FrameGraphPass;
    use crate::pipeline::GpuMemoryPool;
    use crate::pipeline::MultiQueueRecorder;

    use crate::pipeline::PipelineStage;
    use crate::pipeline::PipelineStatistics;
    use crate::pipeline::PipelineStats;
    use crate::pipeline::QueueType;
    use crate::pipeline::ResourceAliasingTracker;
    use crate::pipeline::ResourceBarrier;

    use crate::pipeline::StageTimer;
    use crate::pipeline::TimestampQuery;
    use crate::pipeline::TimestampQuerySet;
    use crate::pipeline::WorldState;
    #[test]
    fn timestamp_query_elapsed() {
        let q = TimestampQuery::new("broadphase", 1.0, 3.5);
        assert!((q.elapsed_ms() - 2.5).abs() < 1e-10);
    }
    #[test]
    fn timestamp_query_zero_duration() {
        let q = TimestampQuery::new("instant", 5.0, 5.0);
        assert!((q.elapsed_ms()).abs() < 1e-10);
    }
    #[test]
    fn query_set_slowest_pass() {
        let mut qs = TimestampQuerySet::new();
        qs.record(TimestampQuery::new("a", 0.0, 1.0));
        qs.record(TimestampQuery::new("b", 0.0, 5.0));
        qs.record(TimestampQuery::new("c", 0.0, 2.0));
        let slowest = qs.slowest_pass().unwrap();
        assert_eq!(slowest.label, "b");
    }
    #[test]
    fn query_set_total_elapsed() {
        let mut qs = TimestampQuerySet::new();
        qs.record(TimestampQuery::new("x", 0.0, 1.0));
        qs.record(TimestampQuery::new("y", 0.0, 2.0));
        assert!((qs.total_elapsed_ms() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn query_set_empty_slowest_none() {
        let qs = TimestampQuerySet::new();
        assert!(qs.slowest_pass().is_none());
    }
    #[test]
    fn query_set_clear() {
        let mut qs = TimestampQuerySet::new();
        qs.record(TimestampQuery::new("a", 0.0, 1.0));
        qs.clear();
        assert!(qs.queries().is_empty());
    }
    #[test]
    fn pipeline_stats_arithmetic_intensity() {
        let s = PipelineStatistics {
            cs_invocations: 1024,
            workgroups_dispatched: 16,
            flops: 2048,
            bytes_read: 512,
            bytes_written: 512,
        };
        let ai = s.arithmetic_intensity();
        assert!((ai - 2.0).abs() < 1e-10, "ai = {ai}");
    }
    #[test]
    fn pipeline_stats_zero_bytes() {
        let s = PipelineStatistics::default();
        assert!((s.arithmetic_intensity()).abs() < 1e-10);
    }
    #[test]
    fn pipeline_stats_bandwidth_utilisation() {
        let s = PipelineStatistics {
            bytes_read: 500_000_000,
            bytes_written: 500_000_000,
            ..Default::default()
        };
        let util = s.bandwidth_utilization(1_000_000_000.0, 1.0);
        assert!((util - 1.0).abs() < 1e-6, "util = {util}");
    }
    #[test]
    fn multi_queue_recorder_submit_and_flush() {
        let mut rec = MultiQueueRecorder::new();
        let p = ComputePipeline::new("t", "", "main");
        rec.submit(DispatchBatch::new(p.clone(), 64), QueueType::Main);
        rec.submit(DispatchBatch::new(p.clone(), 64), QueueType::AsyncCompute);
        rec.submit(DispatchBatch::new(p, 64), QueueType::Transfer);
        assert_eq!(rec.pending_total(), 3);
        let flushed = rec.flush_all();
        assert_eq!(flushed, 3);
        assert_eq!(rec.pending_total(), 0);
    }
    #[test]
    fn multi_queue_recorder_total_recorded() {
        let mut rec = MultiQueueRecorder::new();
        let p = ComputePipeline::new("t", "", "main");
        for _ in 0..5 {
            rec.submit(DispatchBatch::new(p.clone(), 32), QueueType::Main);
        }
        assert_eq!(rec.total_recorded, 5);
    }
    #[test]
    fn barrier_optimizer_removes_duplicates() {
        let barriers = vec![
            ResourceBarrier::new(PipelineStage::BroadPhase, PipelineStage::NarrowPhase, "buf"),
            ResourceBarrier::new(PipelineStage::BroadPhase, PipelineStage::NarrowPhase, "buf"),
        ];
        let opt = BarrierOptimizer::optimize(&barriers);
        assert_eq!(opt.len(), 1, "duplicate should be removed");
        assert_eq!(BarrierOptimizer::savings(&barriers), 1);
    }
    #[test]
    fn barrier_optimizer_preserves_different_resources() {
        let barriers = vec![
            ResourceBarrier::new(
                PipelineStage::BroadPhase,
                PipelineStage::NarrowPhase,
                "buf_a",
            ),
            ResourceBarrier::new(
                PipelineStage::BroadPhase,
                PipelineStage::NarrowPhase,
                "buf_b",
            ),
        ];
        let opt = BarrierOptimizer::optimize(&barriers);
        assert_eq!(opt.len(), 2);
    }
    #[test]
    fn barrier_optimizer_empty_input() {
        let opt = BarrierOptimizer::optimize(&[]);
        assert!(opt.is_empty());
        assert_eq!(BarrierOptimizer::savings(&[]), 0);
    }
    #[test]
    fn aliasing_tracker_detects_shared_allocation() {
        let mut t = ResourceAliasingTracker::new();
        t.track("shadow_map", 0, 1024);
        t.track("gbuffer_depth", 0, 1024);
        assert!(t.are_aliased("shadow_map", "gbuffer_depth"));
    }
    #[test]
    fn aliasing_tracker_no_alias() {
        let mut t = ResourceAliasingTracker::new();
        t.track("a", 0, 512);
        t.track("b", 512, 512);
        assert!(!t.are_aliased("a", "b"));
    }
    #[test]
    fn aliasing_tracker_aliases_for() {
        let mut t = ResourceAliasingTracker::new();
        t.track("x", 100, 200);
        t.track("y", 100, 200);
        let aliases = t.aliases_for(100, 200);
        assert_eq!(aliases.len(), 2);
    }
    #[test]
    fn aliasing_tracker_counts() {
        let mut t = ResourceAliasingTracker::new();
        t.track("a", 0, 256);
        t.track("b", 256, 256);
        t.track("c", 0, 256);
        assert_eq!(t.allocation_count(), 2);
        assert_eq!(t.total_resource_registrations(), 3);
    }
    #[test]
    fn overlap_scheduler_critical_and_background() {
        let mut sched = ComputeOverlapScheduler::new();
        let p = ComputePipeline::new("t", "", "main");
        sched.submit_critical(DispatchBatch::new(p.clone(), 64));
        sched.submit_background(DispatchBatch::new(p, 64));
        assert!(sched.has_pending());
        let n = sched.end_frame();
        assert_eq!(n, 2);
        assert!(!sched.has_pending());
        assert_eq!(sched.critical_count, 0);
    }
    #[test]
    fn overlap_scheduler_empty_frame() {
        let mut sched = ComputeOverlapScheduler::new();
        assert!(!sched.has_pending());
        let n = sched.end_frame();
        assert_eq!(n, 0);
    }
    #[test]
    fn frame_graph_writers_and_readers() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FrameGraphPass::new("shadow_pass", QueueType::Main).writes("shadow_map"));
        fg.add_pass(FrameGraphPass::new("lighting_pass", QueueType::Main).reads("shadow_map"));
        let writers = fg.writers_of("shadow_map");
        assert_eq!(writers.len(), 1);
        assert_eq!(writers[0].name, "shadow_pass");
        let readers = fg.readers_of("shadow_map");
        assert_eq!(readers.len(), 1);
        assert_eq!(readers[0].name, "lighting_pass");
    }
    #[test]
    fn frame_graph_validate_valid_deps() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FrameGraphPass::new("pass_a", QueueType::Main));
        fg.add_pass(FrameGraphPass::new("pass_b", QueueType::Main).depends_on("pass_a"));
        assert!(
            fg.validate_dependencies().is_empty(),
            "valid deps should produce no errors"
        );
    }
    #[test]
    fn frame_graph_validate_invalid_dep() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FrameGraphPass::new("pass_b", QueueType::Main).depends_on("nonexistent"));
        let errors = fg.validate_dependencies();
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("nonexistent"));
    }
    #[test]
    fn frame_graph_async_pass_count() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FrameGraphPass::new("main1", QueueType::Main));
        fg.add_pass(FrameGraphPass::new("async1", QueueType::AsyncCompute));
        fg.add_pass(FrameGraphPass::new("async2", QueueType::AsyncCompute));
        assert_eq!(fg.async_pass_count(), 2);
    }
    #[test]
    fn frame_graph_pass_names_order() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FrameGraphPass::new("z", QueueType::Main));
        fg.add_pass(FrameGraphPass::new("a", QueueType::Main));
        let names = fg.pass_names();
        assert_eq!(names, vec!["z", "a"]);
    }
    #[test]
    fn pool_free_out_of_bounds_returns_err() {
        let mut pool = GpuMemoryPool::new(100);
        assert!(
            pool.free(200, 10).is_err(),
            "out-of-bounds free should fail"
        );
    }
    #[test]
    fn pool_fragmentation_after_partial_free() {
        let mut pool = GpuMemoryPool::new(300);
        let h1 = pool.alloc(100).unwrap();
        let _h2 = pool.alloc(100).unwrap();
        let _h3 = pool.alloc(100).unwrap();
        pool.free(h1.0, h1.1).unwrap();
        assert!(pool.fragmentation_count() >= 1);
    }
    #[test]
    fn stage_timer_records_elapsed() {
        let mut timer = StageTimer::start();
        timer.stop();
        assert!(timer.elapsed_ms >= 0.0, "elapsed should be non-negative");
    }
    #[test]
    fn pipeline_stats_stage_total() {
        let s = PipelineStats {
            broadphase_ms: 1.0,
            narrowphase_ms: 2.0,
            constraint_ms: 3.0,
            integration_ms: 4.0,
            postprocess_ms: 5.0,
            ..Default::default()
        };
        assert!((s.stage_total_ms() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn world_state_body_count() {
        let w = WorldState {
            positions: vec![0.0; 9],
            velocities: vec![0.0; 9],
            inverse_masses: vec![1.0, 1.0, 1.0],
        };
        assert_eq!(w.body_count(), 3);
    }
}
