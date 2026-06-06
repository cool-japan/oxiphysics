//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
use std::time::Instant;

/// Schedules compute dispatches to maximise GPU utilisation through overlapping.
///
/// Classifies dispatches as either "latency-critical" (must finish before the
/// next frame) or "background" (can run async on a secondary queue).
#[derive(Debug)]
pub struct ComputeOverlapScheduler {
    /// Number of latency-critical dispatches submitted this frame.
    pub critical_count: usize,
    /// Number of background dispatches submitted this frame.
    pub background_count: usize,
    pub(super) recorder: MultiQueueRecorder,
}
impl ComputeOverlapScheduler {
    /// Create a new scheduler.
    pub fn new() -> Self {
        Self {
            critical_count: 0,
            background_count: 0,
            recorder: MultiQueueRecorder::new(),
        }
    }
    /// Submit a latency-critical dispatch (→ main queue).
    pub fn submit_critical(&mut self, batch: DispatchBatch) {
        self.critical_count += 1;
        self.recorder.submit(batch, QueueType::Main);
    }
    /// Submit a background dispatch (→ async compute queue).
    pub fn submit_background(&mut self, batch: DispatchBatch) {
        self.background_count += 1;
        self.recorder.submit(batch, QueueType::AsyncCompute);
    }
    /// Execute all pending work and reset frame counters.
    pub fn end_frame(&mut self) -> usize {
        let n = self.recorder.flush_all();
        self.critical_count = 0;
        self.background_count = 0;
        n
    }
    /// Returns `true` if there is any pending work.
    pub fn has_pending(&self) -> bool {
        self.recorder.pending_total() > 0
    }
}
/// Configuration for the high-level physics pipeline.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Which stages are active.  Stages not in this list are skipped.
    pub enabled_stages: Vec<PipelineStage>,
    /// Number of sub-steps per call to [`PhysicsPipeline::step`].
    pub substeps: u32,
    /// Whether a GPU backend should be preferred over the CPU fallback.
    pub use_gpu: bool,
}
impl PipelineConfig {
    /// Default config: all stages enabled, 1 substep, CPU backend.
    pub fn new() -> Self {
        Self {
            enabled_stages: PipelineStage::all_in_order().to_vec(),
            substeps: 1,
            use_gpu: false,
        }
    }
    /// Returns `true` if `stage` is in the enabled list.
    pub fn is_enabled(&self, stage: PipelineStage) -> bool {
        self.enabled_stages.contains(&stage)
    }
}
/// Orchestrates a complete simulation step through all enabled pipeline stages.
///
/// The stage execution order is guaranteed:
/// `BroadPhase → NarrowPhase → ConstraintSolve → Integration → PostProcess`
pub struct PhysicsPipeline {
    /// Pipeline configuration.
    pub config: PipelineConfig,
    /// Cumulative stats across all steps since creation.
    pub stats: PipelineStats,
}
impl PhysicsPipeline {
    /// Create a new pipeline with the given configuration.
    pub fn new(config: PipelineConfig) -> Self {
        Self {
            config,
            stats: PipelineStats::default(),
        }
    }
    /// Run one full simulation step (all enabled stages, for `substeps`
    /// iterations) over `world_state` with time step `dt`.
    ///
    /// Returns the [`PipelineStats`] for this step.
    pub fn step(&mut self, world_state: &mut WorldState, dt: f64) -> PipelineStats {
        let step_start = Instant::now();
        let mut step_stats = PipelineStats::default();
        let sub_dt = if self.config.substeps > 0 {
            dt / self.config.substeps as f64
        } else {
            dt
        };
        for _ in 0..self.config.substeps.max(1) {
            let sub_stats = self.run_stages(world_state, sub_dt);
            step_stats.accumulate(&sub_stats);
        }
        step_stats.total_time_ms = step_start.elapsed().as_secs_f64() * 1000.0;
        self.stats.accumulate(&step_stats);
        step_stats
    }
    /// Run all enabled stages once for a single sub-step.
    fn run_stages(&self, world_state: &mut WorldState, dt: f64) -> PipelineStats {
        let mut stats = PipelineStats::default();
        for stage in PipelineStage::all_in_order() {
            if !self.config.is_enabled(stage) {
                continue;
            }
            let mut timer = StageTimer::start();
            match stage {
                PipelineStage::BroadPhase => {
                    let pairs = run_broadphase(world_state);
                    stats.collision_pairs += pairs as u32;
                }
                PipelineStage::NarrowPhase => {}
                PipelineStage::ConstraintSolve => {
                    let solved = run_constraint_solve(world_state);
                    stats.solved_constraints += solved as u32;
                }
                PipelineStage::Integration => {
                    run_integration(world_state, dt);
                }
                PipelineStage::PostProcess => {
                    run_postprocess(world_state);
                }
            }
            timer.stop();
            match stage {
                PipelineStage::BroadPhase => stats.broadphase_ms += timer.elapsed_ms,
                PipelineStage::NarrowPhase => stats.narrowphase_ms += timer.elapsed_ms,
                PipelineStage::ConstraintSolve => stats.constraint_ms += timer.elapsed_ms,
                PipelineStage::Integration => stats.integration_ms += timer.elapsed_ms,
                PipelineStage::PostProcess => stats.postprocess_ms += timer.elapsed_ms,
            }
        }
        stats
    }
}
/// GPU compute pipeline descriptor (wgpu-agnostic, usable in unit tests).
#[derive(Debug, Clone)]
pub struct ComputePipeline {
    /// Human-readable label for debugging.
    pub label: String,
    /// Raw WGSL shader source.
    pub shader_source: String,
    /// Name of the entry-point function inside the shader.
    pub entry_point: String,
    /// Workgroup size declared in the shader `[X, Y, Z]`.
    pub workgroup_size: [u32; 3],
}
impl ComputePipeline {
    /// Create a new pipeline descriptor.
    ///
    /// `workgroup_size` defaults to `[64, 1, 1]` (matching the shaders in
    /// [`crate::shaders`]).
    pub fn new(label: &str, shader: &str, entry_point: &str) -> Self {
        Self {
            label: label.to_owned(),
            shader_source: shader.to_owned(),
            entry_point: entry_point.to_owned(),
            workgroup_size: [64, 1, 1],
        }
    }
    /// Compute the dispatch grid dimensions needed to cover `n_items` work
    /// items.
    ///
    /// Returns `[ceil(n_items / workgroup_size[0\]), 1, 1]`.
    pub fn workgroups_needed(&self, n_items: u32) -> [u32; 3] {
        let x = n_items.div_ceil(self.workgroup_size[0]);
        [x, 1, 1]
    }
}
/// Describes a memory/execution barrier between two pipeline stages.
///
/// On a real GPU this would translate to a `vkCmdPipelineBarrier` or
/// `wgpu::CommandEncoder::insert_debug_marker`.  Here it is a pure data
/// type used to validate that the user has not forgotten to insert barriers
/// between dependent passes.
#[derive(Debug, Clone, PartialEq)]
pub struct ResourceBarrier {
    /// The stage that *writes* the resource.
    pub src_stage: PipelineStage,
    /// The stage that *reads* the resource.
    pub dst_stage: PipelineStage,
    /// Human-readable name of the protected resource.
    pub resource_name: String,
}
impl ResourceBarrier {
    /// Create a new barrier descriptor.
    pub fn new(src: PipelineStage, dst: PipelineStage, name: &str) -> Self {
        Self {
            src_stage: src,
            dst_stage: dst,
            resource_name: name.to_owned(),
        }
    }
    /// Returns `true` if the barrier is between stages in the correct
    /// dependency order (src must come before dst).
    pub fn is_valid_order(&self) -> bool {
        self.src_stage < self.dst_stage
    }
}
/// An opaque handle to a pool-allocated GPU resource.
///
/// Stores the `(offset, size)` pair returned by [`GpuMemoryPool::alloc`] so
/// that the caller can later free the resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResourceHandle {
    /// Offset into the pool's backing store.
    pub offset: usize,
    /// Size of the allocation in `f32` elements.
    pub size: usize,
}
impl ResourceHandle {
    /// Create a handle from a raw `(offset, size)` pair.
    pub fn from_alloc(alloc: (usize, usize)) -> Self {
        Self {
            offset: alloc.0,
            size: alloc.1,
        }
    }
}
/// A simple frame-graph pass descriptor.
///
/// Frame graphs describe the rendering / compute workload as a DAG of passes,
/// enabling automatic barrier insertion and resource aliasing.
#[derive(Debug, Clone)]
pub struct FrameGraphPass {
    /// Unique name for this pass.
    pub name: String,
    /// Names of resources this pass reads.
    pub reads: Vec<String>,
    /// Names of resources this pass writes.
    pub writes: Vec<String>,
    /// Names of passes that must complete before this one.
    pub dependencies: Vec<String>,
    /// Which queue this pass runs on.
    pub queue: QueueType,
}
impl FrameGraphPass {
    /// Create a new pass with no dependencies.
    pub fn new(name: impl Into<String>, queue: QueueType) -> Self {
        Self {
            name: name.into(),
            reads: Vec::new(),
            writes: Vec::new(),
            dependencies: Vec::new(),
            queue,
        }
    }
    /// Add a resource read.
    pub fn reads(mut self, resource: impl Into<String>) -> Self {
        self.reads.push(resource.into());
        self
    }
    /// Add a resource write.
    pub fn writes(mut self, resource: impl Into<String>) -> Self {
        self.writes.push(resource.into());
        self
    }
    /// Add a pass dependency.
    pub fn depends_on(mut self, pass: impl Into<String>) -> Self {
        self.dependencies.push(pass.into());
        self
    }
}
/// A minimal frame graph that can validate pass dependencies.
#[derive(Debug, Default)]
pub struct FrameGraph {
    pub(super) passes: Vec<FrameGraphPass>,
}
impl FrameGraph {
    /// Create an empty frame graph.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a pass.
    pub fn add_pass(&mut self, pass: FrameGraphPass) {
        self.passes.push(pass);
    }
    /// Return all pass names in insertion order.
    pub fn pass_names(&self) -> Vec<&str> {
        self.passes.iter().map(|p| p.name.as_str()).collect()
    }
    /// Return all passes that write to `resource`.
    pub fn writers_of(&self, resource: &str) -> Vec<&FrameGraphPass> {
        self.passes
            .iter()
            .filter(|p| p.writes.iter().any(|w| w == resource))
            .collect()
    }
    /// Return all passes that read from `resource`.
    pub fn readers_of(&self, resource: &str) -> Vec<&FrameGraphPass> {
        self.passes
            .iter()
            .filter(|p| p.reads.iter().any(|r| r == resource))
            .collect()
    }
    /// Validate that all declared dependencies reference existing passes.
    /// Returns a list of invalid dependency references.
    pub fn validate_dependencies(&self) -> Vec<String> {
        let names: std::collections::HashSet<&str> =
            self.passes.iter().map(|p| p.name.as_str()).collect();
        let mut errors = Vec::new();
        for pass in &self.passes {
            for dep in &pass.dependencies {
                if !names.contains(dep.as_str()) {
                    errors.push(format!("{}: unknown dependency '{}'", pass.name, dep));
                }
            }
        }
        errors
    }
    /// Count how many passes run on the async compute queue.
    pub fn async_pass_count(&self) -> usize {
        self.passes
            .iter()
            .filter(|p| p.queue == QueueType::AsyncCompute)
            .count()
    }
}
/// A distinct stage in the physics simulation pipeline.
///
/// Stages run in dependency order:
/// `BroadPhase → NarrowPhase → ConstraintSolve → Integration → PostProcess`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PipelineStage {
    /// Broad-phase collision detection (AABB / BVH).
    BroadPhase,
    /// Narrow-phase contact generation.
    NarrowPhase,
    /// Constraint / impulse solver.
    ConstraintSolve,
    /// Velocity and position integration.
    Integration,
    /// Post-processing (sleeping, user callbacks, …).
    PostProcess,
}
impl PipelineStage {
    /// Ordered list of all stages, from first to last.
    pub fn all_in_order() -> [PipelineStage; 5] {
        [
            PipelineStage::BroadPhase,
            PipelineStage::NarrowPhase,
            PipelineStage::ConstraintSolve,
            PipelineStage::Integration,
            PipelineStage::PostProcess,
        ]
    }
}
/// A single compute dispatch pass: one pipeline plus its bound buffers.
#[derive(Debug)]
pub struct DispatchBatch {
    /// The pipeline to dispatch.
    pub pipeline: ComputePipeline,
    /// Buffers bound to this pass, in binding-slot order.
    pub bindings: Vec<CpuBuffer>,
    /// The dispatch grid dimensions `[X, Y, Z]`.
    pub dispatch_dims: [u32; 3],
}
impl DispatchBatch {
    /// Create a new dispatch batch, computing the required workgroup grid from
    /// `workitems`.
    pub fn new(pipeline: ComputePipeline, workitems: u32) -> Self {
        let dispatch_dims = pipeline.workgroups_needed(workitems);
        Self {
            pipeline,
            bindings: Vec::new(),
            dispatch_dims,
        }
    }
    /// Append `buffer` to the binding list (next available slot).
    pub fn bind(&mut self, buffer: CpuBuffer) {
        self.bindings.push(buffer);
    }
}
/// Models an async compute queue: a FIFO of [`DispatchBatch`]es that can be
/// submitted independently of the graphics queue.
///
/// On a GPU this enables overlapping async compute with graphics work.
/// Here it is a simple in-memory queue used to test scheduling logic.
pub struct AsyncComputeQueue {
    pub(super) queue: std::collections::VecDeque<DispatchBatch>,
    /// Total number of batches ever enqueued (for profiling).
    pub total_enqueued: usize,
    /// Total number of batches ever drained (executed).
    pub total_executed: usize,
}
impl AsyncComputeQueue {
    /// Create an empty async compute queue.
    pub fn new() -> Self {
        Self {
            queue: std::collections::VecDeque::new(),
            total_enqueued: 0,
            total_executed: 0,
        }
    }
    /// Enqueue a dispatch batch.
    pub fn submit(&mut self, batch: DispatchBatch) {
        self.total_enqueued += 1;
        self.queue.push_back(batch);
    }
    /// Execute (drain) all pending batches and return the count executed.
    ///
    /// In a real backend this would flush the command buffer to the GPU.
    /// Here it simply clears the queue and updates counters.
    pub fn flush(&mut self) -> usize {
        let n = self.queue.len();
        self.total_executed += n;
        self.queue.clear();
        n
    }
    /// Number of batches currently waiting in the queue.
    pub fn pending(&self) -> usize {
        self.queue.len()
    }
    /// Returns `true` when there are no pending batches.
    pub fn is_idle(&self) -> bool {
        self.queue.is_empty()
    }
}
/// A dispatched batch annotated with its target queue.
#[derive(Debug)]
pub struct MultiQueueBatch {
    /// The dispatch batch.
    pub batch: DispatchBatch,
    /// Which queue to submit to.
    pub queue: QueueType,
    /// Semaphore/fence dependency: must wait for this frame's signal.
    pub wait_frame: u64,
}
/// Minimal world state passed to [`PhysicsPipeline::step`].
///
/// Uses pure f64 arrays (no nalgebra) compatible with GPU buffer uploads.
#[derive(Debug, Clone, Default)]
pub struct WorldState {
    /// Flat position array `[x0, y0, z0, x1, y1, z1, ...]`.
    pub positions: Vec<f64>,
    /// Flat velocity array `[vx0, vy0, vz0, ...]`.
    pub velocities: Vec<f64>,
    /// Inverse masses `[inv_m0, inv_m1, ...]`.
    pub inverse_masses: Vec<f64>,
}
impl WorldState {
    /// Number of bodies in the world state.
    pub fn body_count(&self) -> usize {
        self.inverse_masses.len()
    }
}
/// Records per-stage timing samples and provides basic statistics.
///
/// Each call to [`PipelineProfiler::record`] appends one sample for the
/// named stage.  [`PipelineProfiler::summary`] returns mean ± stddev.
#[derive(Debug, Default)]
pub struct PipelineProfiler {
    pub(super) samples: std::collections::HashMap<String, Vec<f64>>,
}
impl PipelineProfiler {
    /// Create a new, empty profiler.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a timing sample (milliseconds) for `stage_name`.
    pub fn record(&mut self, stage_name: &str, ms: f64) {
        self.samples
            .entry(stage_name.to_owned())
            .or_default()
            .push(ms);
    }
    /// Returns `(mean_ms, stddev_ms, sample_count)` for `stage_name`, or
    /// `None` if no samples have been recorded.
    pub fn summary(&self, stage_name: &str) -> Option<(f64, f64, usize)> {
        let v = self.samples.get(stage_name)?;
        if v.is_empty() {
            return None;
        }
        let n = v.len() as f64;
        let mean = v.iter().sum::<f64>() / n;
        let variance = v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
        Some((mean, variance.sqrt(), v.len()))
    }
    /// Return all stage names that have recorded samples.
    pub fn stage_names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.samples.keys().map(String::as_str).collect();
        names.sort_unstable();
        names
    }
    /// Total number of samples across all stages.
    pub fn total_samples(&self) -> usize {
        self.samples.values().map(Vec::len).sum()
    }
    /// Clear all recorded samples.
    pub fn reset(&mut self) {
        self.samples.clear();
    }
}
/// Models multiple compute queues (graphics queue + async compute queue).
///
/// On Vulkan/Metal/DX12 certain hardware can overlap work on the graphics
/// queue and an async compute queue.  This mock records which queue each
/// dispatch was submitted to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueueType {
    /// Main graphics/compute queue.
    Main,
    /// Async compute queue (can overlap with graphics).
    AsyncCompute,
    /// Transfer/copy queue.
    Transfer,
}
/// Wraps an [`Instant`] and records the elapsed duration of a stage.
pub struct StageTimer {
    pub(super) start: Instant,
    /// Elapsed time recorded by [`StageTimer::stop`].
    pub elapsed_ms: f64,
}
impl StageTimer {
    /// Create a new timer and start the clock.
    pub fn start() -> Self {
        Self {
            start: Instant::now(),
            elapsed_ms: 0.0,
        }
    }
    /// Stop the timer and record elapsed time in milliseconds.
    pub fn stop(&mut self) {
        self.elapsed_ms = self.start.elapsed().as_secs_f64() * 1000.0;
    }
}
/// Low-level GPU pipeline statistics for one dispatch (mock).
///
/// Mirrors the statistics provided by `VkQueryPool` with
/// `VK_QUERY_TYPE_PIPELINE_STATISTICS` or Metal/D3D equivalents.
#[derive(Debug, Clone, Default)]
pub struct PipelineStatistics {
    /// Number of invocations of the compute shader.
    pub cs_invocations: u64,
    /// Number of work-groups dispatched.
    pub workgroups_dispatched: u64,
    /// Estimated FLOP count (mock value).
    pub flops: u64,
    /// Estimated bytes read from global memory.
    pub bytes_read: u64,
    /// Estimated bytes written to global memory.
    pub bytes_written: u64,
}
impl PipelineStatistics {
    /// Compute the arithmetic intensity (FLOPs per byte).
    pub fn arithmetic_intensity(&self) -> f64 {
        let bytes = self.bytes_read + self.bytes_written;
        if bytes == 0 {
            return 0.0;
        }
        self.flops as f64 / bytes as f64
    }
    /// Estimated memory bandwidth utilization as a fraction of `peak_bw_bytes_s`.
    pub fn bandwidth_utilization(&self, peak_bw_bytes_s: f64, elapsed_s: f64) -> f64 {
        if peak_bw_bytes_s <= 0.0 || elapsed_s <= 0.0 {
            return 0.0;
        }
        let used = (self.bytes_read + self.bytes_written) as f64 / elapsed_s;
        (used / peak_bw_bytes_s).min(1.0)
    }
}
/// Multi-queue command recorder: tracks which batches go to which queue.
#[derive(Debug, Default)]
pub struct MultiQueueRecorder {
    /// Pending batches on the main queue.
    pub main_queue: Vec<DispatchBatch>,
    /// Pending batches on the async compute queue.
    pub async_queue: Vec<DispatchBatch>,
    /// Pending batches on the transfer queue.
    pub transfer_queue: Vec<DispatchBatch>,
    /// Total batches ever recorded.
    pub total_recorded: usize,
}
impl MultiQueueRecorder {
    /// Create an empty recorder.
    pub fn new() -> Self {
        Self::default()
    }
    /// Submit a batch to the specified queue.
    pub fn submit(&mut self, batch: DispatchBatch, queue: QueueType) {
        self.total_recorded += 1;
        match queue {
            QueueType::Main => self.main_queue.push(batch),
            QueueType::AsyncCompute => self.async_queue.push(batch),
            QueueType::Transfer => self.transfer_queue.push(batch),
        }
    }
    /// Flush all queues and return total batches executed.
    pub fn flush_all(&mut self) -> usize {
        let n = self.main_queue.len() + self.async_queue.len() + self.transfer_queue.len();
        self.main_queue.clear();
        self.async_queue.clear();
        self.transfer_queue.clear();
        n
    }
    /// Number of pending batches across all queues.
    pub fn pending_total(&self) -> usize {
        self.main_queue.len() + self.async_queue.len() + self.transfer_queue.len()
    }
}
/// A simple slab-style GPU memory pool that sub-allocates [`CpuBuffer`]s from
/// a fixed-capacity backing store.
///
/// The pool pre-allocates a large buffer and hands out non-overlapping slices
/// (as [`CpuBuffer`] views backed by offsets).  On a real GPU this avoids per-
/// allocation overhead from `vkAllocateMemory`.
///
/// This CPU-side mock tracks allocations as `(offset, size)` pairs and
/// simulates fragmentation/free-list behaviour.
#[derive(Debug)]
pub struct GpuMemoryPool {
    /// Total capacity of the pool in `f32` elements.
    pub capacity: usize,
    /// Currently allocated `f32` elements.
    pub allocated: usize,
    /// Free-list entries `(offset_in_elements, size_in_elements)`.
    pub(super) free_list: Vec<(usize, usize)>,
}
impl GpuMemoryPool {
    /// Create a pool with the given capacity in `f32` elements.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            allocated: 0,
            free_list: vec![(0, capacity)],
        }
    }
    /// Attempt to allocate `size` `f32` elements from the pool.
    ///
    /// Returns `Some((offset, size))` on success, `None` if the pool is full.
    /// Uses a first-fit strategy.
    pub fn alloc(&mut self, size: usize) -> Option<(usize, usize)> {
        for i in 0..self.free_list.len() {
            let (off, avail) = self.free_list[i];
            if avail >= size {
                let alloc_off = off;
                if avail == size {
                    self.free_list.remove(i);
                } else {
                    self.free_list[i] = (off + size, avail - size);
                }
                self.allocated += size;
                return Some((alloc_off, size));
            }
        }
        None
    }
    /// Free a previously allocated block `(offset, size)`.
    ///
    /// Returns `Err` if the block was never allocated (invalid free).
    pub fn free(&mut self, offset: usize, size: usize) -> Result<(), &'static str> {
        if offset + size > self.capacity {
            return Err("block out of bounds");
        }
        if self.allocated < size {
            return Err("double-free: allocated count underflow");
        }
        self.allocated -= size;
        self.free_list.push((offset, size));
        self.free_list.sort_by_key(|&(off, _)| off);
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for &(off, sz) in &self.free_list {
            if let Some(last) = merged.last_mut()
                && last.0 + last.1 == off
            {
                last.1 += sz;
                continue;
            }
            merged.push((off, sz));
        }
        self.free_list = merged;
        Ok(())
    }
    /// Remaining free `f32` elements.
    pub fn free_space(&self) -> usize {
        self.capacity - self.allocated
    }
    /// Returns `true` if the pool has no outstanding allocations.
    pub fn is_fully_free(&self) -> bool {
        self.allocated == 0
    }
    /// Number of fragmented free-list entries (1 = perfectly contiguous).
    pub fn fragmentation_count(&self) -> usize {
        self.free_list.len()
    }
    /// Allocate a named [`CpuBuffer`] from the pool, returning the buffer
    /// and the pool allocation handle `(offset, size)`.
    pub fn alloc_buffer(
        &mut self,
        label: &str,
        n: usize,
        usage: BufferUsage,
    ) -> Option<(CpuBuffer, (usize, usize))> {
        let handle = self.alloc(n)?;
        let buf = CpuBuffer::new_zeros(label, n, usage);
        Some((buf, handle))
    }
}
/// A set of [`TimestampQuery`]s collected during one frame / step.
#[derive(Debug, Clone, Default)]
pub struct TimestampQuerySet {
    pub(super) queries: Vec<TimestampQuery>,
}
impl TimestampQuerySet {
    /// Create an empty query set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a timestamp query.
    pub fn record(&mut self, query: TimestampQuery) {
        self.queries.push(query);
    }
    /// Return all queries.
    pub fn queries(&self) -> &[TimestampQuery] {
        &self.queries
    }
    /// Find the query with the longest elapsed time.
    pub fn slowest_pass(&self) -> Option<&TimestampQuery> {
        self.queries.iter().max_by(|a, b| {
            a.elapsed_ms()
                .partial_cmp(&b.elapsed_ms())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
    /// Sum of all elapsed times.
    pub fn total_elapsed_ms(&self) -> f64 {
        self.queries.iter().map(|q| q.elapsed_ms()).sum()
    }
    /// Clear all recorded queries.
    pub fn clear(&mut self) {
        self.queries.clear();
    }
}
/// Per-step statistics produced by [`PhysicsPipeline::step`].
#[derive(Debug, Clone, Default)]
pub struct PipelineStats {
    /// Time spent in the broad-phase stage (ms).
    pub broadphase_ms: f64,
    /// Time spent in the narrow-phase stage (ms).
    pub narrowphase_ms: f64,
    /// Time spent in the constraint-solve stage (ms).
    pub constraint_ms: f64,
    /// Time spent in the integration stage (ms).
    pub integration_ms: f64,
    /// Time spent in the post-process stage (ms).
    pub postprocess_ms: f64,
    /// Total wall-clock time for the step (ms).
    pub total_time_ms: f64,
    /// Number of collision pairs found in broad-phase.
    pub collision_pairs: u32,
    /// Number of constraints solved.
    pub solved_constraints: u32,
}
impl PipelineStats {
    /// Accumulate another stats snapshot into `self` (used for sub-step sums).
    pub fn accumulate(&mut self, other: &PipelineStats) {
        self.broadphase_ms += other.broadphase_ms;
        self.narrowphase_ms += other.narrowphase_ms;
        self.constraint_ms += other.constraint_ms;
        self.integration_ms += other.integration_ms;
        self.postprocess_ms += other.postprocess_ms;
        self.total_time_ms += other.total_time_ms;
        self.collision_pairs += other.collision_pairs;
        self.solved_constraints += other.solved_constraints;
    }
    /// Sum of all per-stage times.
    pub fn stage_total_ms(&self) -> f64 {
        self.broadphase_ms
            + self.narrowphase_ms
            + self.constraint_ms
            + self.integration_ms
            + self.postprocess_ms
    }
}
/// Optimises a sequence of resource barriers by removing redundant ones.
///
/// Two barriers are redundant when a later barrier subsumes an earlier one
/// (same `src_stage → dst_stage` pair for the same resource).
#[derive(Debug, Default)]
pub struct BarrierOptimizer;
impl BarrierOptimizer {
    /// Deduplicate `barriers`: for each `(src, dst, resource)` triple,
    /// keep only the last occurrence.  Returns the optimised set.
    pub fn optimize(barriers: &[ResourceBarrier]) -> BarrierSet {
        let mut seen: std::collections::HashMap<
            (PipelineStage, PipelineStage, &str),
            &ResourceBarrier,
        > = std::collections::HashMap::new();
        for b in barriers {
            seen.insert((b.src_stage, b.dst_stage, b.resource_name.as_str()), b);
        }
        let mut out = BarrierSet::new();
        for b in seen.values() {
            out.add(ResourceBarrier::new(
                b.src_stage,
                b.dst_stage,
                &b.resource_name,
            ));
        }
        out
    }
    /// Count how many barriers would be removed from `before` to produce
    /// the optimised set.
    pub fn savings(barriers: &[ResourceBarrier]) -> usize {
        let optimized = Self::optimize(barriers);
        barriers.len().saturating_sub(optimized.len())
    }
}
/// Fluent builder for [`PhysicsPipeline`].
///
/// ```no_run
/// use oxiphysics_gpu::pipeline::{PipelineBuilder, PipelineStage};
///
/// let pipeline = PipelineBuilder::new()
///     .substeps(2)
///     .use_gpu(false)
///     .disable_stage(PipelineStage::PostProcess)
///     .build();
///
/// assert_eq!(pipeline.config.substeps, 2);
/// assert!(!pipeline.config.is_enabled(PipelineStage::PostProcess));
/// ```
pub struct PipelineBuilder {
    pub(super) config: PipelineConfig,
}
impl PipelineBuilder {
    /// Create a builder with a default config (all stages, 1 substep, CPU).
    pub fn new() -> Self {
        Self {
            config: PipelineConfig::new(),
        }
    }
    /// Set the number of sub-steps.
    pub fn substeps(mut self, n: u32) -> Self {
        self.config.substeps = n;
        self
    }
    /// Set whether to prefer a GPU backend.
    pub fn use_gpu(mut self, gpu: bool) -> Self {
        self.config.use_gpu = gpu;
        self
    }
    /// Add a stage to the enabled list (idempotent).
    pub fn enable_stage(mut self, stage: PipelineStage) -> Self {
        if !self.config.enabled_stages.contains(&stage) {
            self.config.enabled_stages.push(stage);
            self.config.enabled_stages.sort();
        }
        self
    }
    /// Remove a stage from the enabled list.
    pub fn disable_stage(mut self, stage: PipelineStage) -> Self {
        self.config.enabled_stages.retain(|&s| s != stage);
        self
    }
    /// Consume the builder and produce a [`PhysicsPipeline`].
    pub fn build(self) -> PhysicsPipeline {
        PhysicsPipeline::new(self.config)
    }
}
/// A collection of [`ResourceBarrier`]s that form a complete barrier schedule
/// for one pipeline configuration.
#[derive(Debug, Clone, Default)]
pub struct BarrierSet {
    pub(super) barriers: Vec<ResourceBarrier>,
}
impl BarrierSet {
    /// Create an empty barrier set.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a barrier.
    pub fn add(&mut self, barrier: ResourceBarrier) {
        self.barriers.push(barrier);
    }
    /// Return all barriers whose `src_stage` matches `stage`.
    pub fn barriers_from(&self, stage: PipelineStage) -> Vec<&ResourceBarrier> {
        self.barriers
            .iter()
            .filter(|b| b.src_stage == stage)
            .collect()
    }
    /// Return all barriers whose `dst_stage` matches `stage`.
    pub fn barriers_to(&self, stage: PipelineStage) -> Vec<&ResourceBarrier> {
        self.barriers
            .iter()
            .filter(|b| b.dst_stage == stage)
            .collect()
    }
    /// Number of barriers registered.
    pub fn len(&self) -> usize {
        self.barriers.len()
    }
    /// Returns `true` if no barriers are registered.
    pub fn is_empty(&self) -> bool {
        self.barriers.is_empty()
    }
    /// Validate that every barrier has `src < dst` (i.e., no backwards
    /// dependencies).  Returns the list of invalid barriers.
    pub fn validate(&self) -> Vec<&ResourceBarrier> {
        self.barriers
            .iter()
            .filter(|b| !b.is_valid_order())
            .collect()
    }
}
/// A timestamp query pair `(begin, end)` for a named pass.
///
/// On a real GPU this maps to `VkQueryPool` or `wgpu::QuerySet` writes at the
/// start and end of a render/compute pass.  Here we store wall-clock
/// `f64` timestamps in milliseconds.
#[derive(Debug, Clone)]
pub struct TimestampQuery {
    /// Human-readable label for this pass.
    pub label: String,
    /// Begin timestamp in milliseconds (relative to some epoch).
    pub begin_ms: f64,
    /// End timestamp in milliseconds.
    pub end_ms: f64,
}
impl TimestampQuery {
    /// Create a new timestamp query pair.
    pub fn new(label: impl Into<String>, begin_ms: f64, end_ms: f64) -> Self {
        Self {
            label: label.into(),
            begin_ms,
            end_ms,
        }
    }
    /// Elapsed duration in milliseconds.
    pub fn elapsed_ms(&self) -> f64 {
        self.end_ms - self.begin_ms
    }
}
/// Tracks resource aliasing: two logical resources that share the same physical
/// memory (pool allocation) at non-overlapping lifetimes.
///
/// This is important for reducing peak GPU memory usage.
#[derive(Debug, Default)]
pub struct ResourceAliasingTracker {
    /// Maps physical allocation handle `(offset, size)` to a list of logical
    /// resource names currently using it.
    pub(super) aliases: std::collections::HashMap<(usize, usize), Vec<String>>,
}
impl ResourceAliasingTracker {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record that `resource_name` uses physical allocation `(offset, size)`.
    pub fn track(&mut self, resource_name: impl Into<String>, offset: usize, size: usize) {
        self.aliases
            .entry((offset, size))
            .or_default()
            .push(resource_name.into());
    }
    /// Return all resource names sharing the physical allocation `(offset, size)`.
    pub fn aliases_for(&self, offset: usize, size: usize) -> &[String] {
        self.aliases
            .get(&(offset, size))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
    /// Return `true` if two named resources share any physical allocation.
    pub fn are_aliased(&self, a: &str, b: &str) -> bool {
        for names in self.aliases.values() {
            if names.contains(&a.to_string()) && names.contains(&b.to_string()) {
                return true;
            }
        }
        false
    }
    /// Total number of physical allocations tracked.
    pub fn allocation_count(&self) -> usize {
        self.aliases.len()
    }
    /// Total number of logical resource registrations across all allocations.
    pub fn total_resource_registrations(&self) -> usize {
        self.aliases.values().map(Vec::len).sum()
    }
}
/// Simulated GPU buffer backed by CPU memory for unit testing without a GPU.
#[derive(Debug, Clone)]
pub struct CpuBuffer {
    /// Human-readable label for debugging.
    pub label: String,
    /// Buffer contents as 32-bit floats.
    pub data: Vec<f32>,
    /// Intended usage of this buffer.
    pub usage: BufferUsage,
}
impl CpuBuffer {
    /// Create a buffer pre-filled with the given data.
    pub fn new_f32(label: &str, data: Vec<f32>, usage: BufferUsage) -> Self {
        Self {
            label: label.to_owned(),
            data,
            usage,
        }
    }
    /// Create a buffer of `n` zeros.
    pub fn new_zeros(label: &str, n: usize, usage: BufferUsage) -> Self {
        Self {
            label: label.to_owned(),
            data: vec![0.0_f32; n],
            usage,
        }
    }
    /// Number of `f32` elements in the buffer.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Returns `true` if the buffer contains no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
/// Describes how a [`CpuBuffer`] will be used by the pipeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BufferUsage {
    /// Read-write storage buffer.
    Storage,
    /// Uniform (read-only, small) buffer.
    Uniform,
    /// Read-only storage buffer.
    StorageReadOnly,
}
