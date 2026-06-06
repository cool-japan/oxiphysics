// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU workload scheduler.
//!
//! Provides task graphs, topological scheduling, resource barriers, async
//! compute simulation, frame graphs, and timestamp queries — all CPU-side.

use std::collections::{HashMap, HashSet, VecDeque};

// ---------------------------------------------------------------------------
// TaskPriority
// ---------------------------------------------------------------------------

/// Priority level for a compute task.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum TaskPriority {
    /// Must complete within the current frame; highest urgency.
    RealTime = 4,
    /// High-importance work that should not be deferred.
    High = 3,
    /// Standard work.
    #[default]
    Normal = 2,
    /// Can be deferred to future frames if time is tight.
    Low = 1,
    /// Background processing with no frame deadline.
    Background = 0,
}

// ---------------------------------------------------------------------------
// ComputeTask
// ---------------------------------------------------------------------------

/// A single GPU compute dispatch.
#[derive(Debug, Clone)]
pub struct ComputeTask {
    /// Unique task name.
    pub name: String,
    /// Workgroup size in (x, y, z).
    pub workgroup_size: [u32; 3],
    /// Dispatch count in (x, y, z).
    pub dispatch_count: [u32; 3],
    /// Names of tasks that must complete before this one.
    pub dependencies: Vec<String>,
    /// Priority for the scheduler.
    pub priority: TaskPriority,
    /// Estimated execution time in milliseconds.
    pub estimated_ms: f64,
}

impl ComputeTask {
    /// Create a simple 1-D compute task.
    pub fn new_1d(name: impl Into<String>, dispatch_x: u32) -> Self {
        Self {
            name: name.into(),
            workgroup_size: [64, 1, 1],
            dispatch_count: [dispatch_x, 1, 1],
            dependencies: vec![],
            priority: TaskPriority::Normal,
            estimated_ms: 1.0,
        }
    }

    /// Create a 2-D compute task.
    pub fn new_2d(name: impl Into<String>, dispatch_x: u32, dispatch_y: u32) -> Self {
        Self {
            name: name.into(),
            workgroup_size: [8, 8, 1],
            dispatch_count: [dispatch_x, dispatch_y, 1],
            dependencies: vec![],
            priority: TaskPriority::Normal,
            estimated_ms: 1.0,
        }
    }

    /// Total number of workgroup invocations.
    pub fn total_workgroups(&self) -> u64 {
        self.dispatch_count[0] as u64
            * self.dispatch_count[1] as u64
            * self.dispatch_count[2] as u64
    }

    /// Total number of shader invocations.
    pub fn total_invocations(&self) -> u64 {
        self.total_workgroups()
            * self.workgroup_size[0] as u64
            * self.workgroup_size[1] as u64
            * self.workgroup_size[2] as u64
    }

    /// Add a dependency by name.
    pub fn depends_on(mut self, dep: impl Into<String>) -> Self {
        self.dependencies.push(dep.into());
        self
    }

    /// Set priority.
    pub fn with_priority(mut self, priority: TaskPriority) -> Self {
        self.priority = priority;
        self
    }

    /// Set estimated execution time.
    pub fn with_estimated_ms(mut self, ms: f64) -> Self {
        self.estimated_ms = ms;
        self
    }
}

// ---------------------------------------------------------------------------
// TaskGraph
// ---------------------------------------------------------------------------

/// A directed acyclic graph of compute tasks.
#[derive(Debug, Clone, Default)]
pub struct TaskGraph {
    /// All tasks keyed by name.
    tasks: HashMap<String, ComputeTask>,
}

impl TaskGraph {
    /// Create an empty task graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a task to the graph.  Replaces any existing task with the same name.
    pub fn add_task(&mut self, task: ComputeTask) {
        self.tasks.insert(task.name.clone(), task);
    }

    /// Remove a task by name.
    pub fn remove_task(&mut self, name: &str) {
        self.tasks.remove(name);
    }

    /// Number of tasks.
    pub fn len(&self) -> usize {
        self.tasks.len()
    }

    /// True when there are no tasks.
    pub fn is_empty(&self) -> bool {
        self.tasks.is_empty()
    }

    /// Topological sort using Kahn's algorithm.
    ///
    /// Returns `Ok(order)` where `order` is a valid execution order, or
    /// `Err(cycle)` naming one task involved in a cycle.
    pub fn topological_sort(&self) -> Result<Vec<String>, String> {
        // Build adjacency and in-degree maps
        let mut in_degree: HashMap<&str, usize> = HashMap::new();
        let mut rev: HashMap<&str, Vec<&str>> = HashMap::new(); // task -> tasks that depend on it

        for (name, task) in &self.tasks {
            in_degree.entry(name.as_str()).or_insert(0);
            for dep in &task.dependencies {
                if !self.tasks.contains_key(dep.as_str()) {
                    // Unknown dependency — skip
                    continue;
                }
                // dep -> name (name depends on dep)
                rev.entry(dep.as_str()).or_default().push(name.as_str());
                *in_degree.entry(name.as_str()).or_insert(0) += 1;
            }
        }

        let mut queue: VecDeque<&str> = in_degree
            .iter()
            .filter(|(_, d)| **d == 0)
            .map(|(&n, _)| n)
            .collect();

        // Sort for determinism
        let mut queue_vec: Vec<&str> = queue.drain(..).collect();
        queue_vec.sort();
        queue.extend(queue_vec);

        let mut order = Vec::new();
        while let Some(name) = queue.pop_front() {
            order.push(name.to_owned());
            if let Some(dependents) = rev.get(name) {
                let mut next: Vec<&str> = dependents
                    .iter()
                    .filter_map(|&d| {
                        let deg = in_degree.get_mut(d)?;
                        *deg -= 1;
                        if *deg == 0 { Some(d) } else { None }
                    })
                    .collect();
                next.sort();
                queue.extend(next);
            }
        }

        if order.len() != self.tasks.len() {
            // Find a node still in a cycle
            let cycle_node = self
                .tasks
                .keys()
                .find(|n| !order.contains(*n))
                .cloned()
                .unwrap_or_else(|| "unknown".to_owned());
            Err(cycle_node)
        } else {
            Ok(order)
        }
    }

    /// Compute the critical path (longest chain by estimated_ms).
    ///
    /// Returns the list of task names on the critical path.
    pub fn critical_path(&self) -> Vec<String> {
        let order = match self.topological_sort() {
            Ok(o) => o,
            Err(_) => return vec![],
        };

        // Compute earliest finish times
        let mut eft: HashMap<&str, f64> = HashMap::new();
        let mut pred: HashMap<&str, &str> = HashMap::new();

        for name in &order {
            let task = &self.tasks[name.as_str()];
            let dep_max = task
                .dependencies
                .iter()
                .filter_map(|d| eft.get(d.as_str()).copied())
                .fold(0.0f64, f64::max);
            let ef = dep_max + task.estimated_ms;
            eft.insert(name.as_str(), ef);
            // Track predecessor that provides dep_max
            if let Some(best_pred) = task
                .dependencies
                .iter()
                .filter_map(|d| {
                    let t = eft.get(d.as_str()).copied()?;
                    Some((d.as_str(), t))
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(d, _)| d)
            {
                pred.insert(name.as_str(), best_pred);
            }
        }

        // Walk back from the task with maximum eft
        let end = order.iter().max_by(|a, b| {
            eft.get(a.as_str())
                .unwrap_or(&0.0)
                .partial_cmp(eft.get(b.as_str()).unwrap_or(&0.0))
                .expect("operation should succeed")
        });

        let mut path = Vec::new();
        let mut cur = match end {
            Some(s) => s.as_str(),
            None => return vec![],
        };
        loop {
            path.push(cur.to_owned());
            match pred.get(cur) {
                Some(&p) => cur = p,
                None => break,
            }
        }
        path.reverse();
        path
    }

    /// Check whether the graph contains a cycle.
    pub fn has_cycle(&self) -> bool {
        self.topological_sort().is_err()
    }
}

// ---------------------------------------------------------------------------
// ResourceBarrier
// ---------------------------------------------------------------------------

/// Type of resource access ordering barrier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarrierType {
    /// Read-after-write: a reader must wait for a preceding writer.
    ReadAfterWrite,
    /// Write-after-read: a writer must wait for a preceding reader.
    WriteAfterRead,
    /// Write-after-write: serialise two writers.
    WriteAfterWrite,
}

/// A resource barrier between two tasks.
#[derive(Debug, Clone)]
pub struct ResourceBarrier {
    /// Name of the task that writes / produces.
    pub producer: String,
    /// Name of the task that reads / consumes.
    pub consumer: String,
    /// The kind of hazard this barrier prevents.
    pub barrier_type: BarrierType,
    /// Resource being protected (e.g. buffer name).
    pub resource: String,
}

impl ResourceBarrier {
    /// Create a read-after-write barrier.
    pub fn raw(
        producer: impl Into<String>,
        consumer: impl Into<String>,
        resource: impl Into<String>,
    ) -> Self {
        Self {
            producer: producer.into(),
            consumer: consumer.into(),
            barrier_type: BarrierType::ReadAfterWrite,
            resource: resource.into(),
        }
    }

    /// Create a write-after-read barrier.
    pub fn war(
        producer: impl Into<String>,
        consumer: impl Into<String>,
        resource: impl Into<String>,
    ) -> Self {
        Self {
            producer: producer.into(),
            consumer: consumer.into(),
            barrier_type: BarrierType::WriteAfterRead,
            resource: resource.into(),
        }
    }
}

// ---------------------------------------------------------------------------
// TaskScheduler
// ---------------------------------------------------------------------------

/// Schedules a task graph into an execution order.
#[derive(Debug, Default)]
pub struct TaskScheduler {
    /// Barriers to inject between tasks.
    pub barriers: Vec<ResourceBarrier>,
}

impl TaskScheduler {
    /// Create a new scheduler.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a resource barrier.
    pub fn add_barrier(&mut self, barrier: ResourceBarrier) {
        self.barriers.push(barrier);
    }

    /// Schedule a task graph.
    ///
    /// Returns the topological execution order, or an error if a cycle exists.
    pub fn schedule(&self, graph: &TaskGraph) -> Result<Vec<String>, String> {
        graph.topological_sort()
    }

    /// Schedule and group independent tasks into parallel batches.
    ///
    /// Each inner `Vec` contains tasks that can run concurrently.
    pub fn batch_schedule(&self, graph: &TaskGraph) -> Result<Vec<Vec<String>>, String> {
        let order = self.schedule(graph)?;
        let tasks = &graph.tasks;

        // Compute the depth of each task (longest dependency chain)
        let mut depth: HashMap<&str, usize> = HashMap::new();
        for name in &order {
            let task = &tasks[name.as_str()];
            let d = task
                .dependencies
                .iter()
                .filter_map(|dep| depth.get(dep.as_str()).copied())
                .max()
                .map(|m| m + 1)
                .unwrap_or(0);
            depth.insert(name.as_str(), d);
        }

        let max_depth = depth.values().copied().max().unwrap_or(0);
        let mut batches: Vec<Vec<String>> = vec![vec![]; max_depth + 1];
        for name in &order {
            let d = *depth.get(name.as_str()).unwrap_or(&0);
            batches[d].push(name.clone());
        }
        Ok(batches)
    }
}

// ---------------------------------------------------------------------------
// WorkloadBalancer
// ---------------------------------------------------------------------------

/// Splits large dispatches across frames to stay within a time budget.
#[derive(Debug, Clone)]
pub struct WorkloadBalancer {
    /// GPU time budget per frame in milliseconds.
    pub budget_ms: f64,
    /// Accumulated pending tasks with their estimated costs.
    pending: Vec<(ComputeTask, f64)>,
}

impl WorkloadBalancer {
    /// Create a new balancer with the given budget.
    pub fn new(budget_ms: f64) -> Self {
        Self {
            budget_ms,
            pending: vec![],
        }
    }

    /// Submit a task for scheduling.
    pub fn submit(&mut self, task: ComputeTask) {
        let cost = task.estimated_ms;
        self.pending.push((task, cost));
    }

    /// Extract tasks that fit within this frame's budget.
    ///
    /// Higher-priority tasks are selected first.  Returns the tasks to
    /// execute this frame.
    pub fn extract_frame_work(&mut self) -> Vec<ComputeTask> {
        // Sort by descending priority then descending estimated_ms
        self.pending.sort_by(|a, b| {
            b.0.priority
                .cmp(&a.0.priority)
                .then(b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
        });

        let mut remaining = self.budget_ms;
        let mut this_frame = Vec::new();
        let mut leftover = Vec::new();

        for (task, cost) in self.pending.drain(..) {
            if cost <= remaining || this_frame.is_empty() {
                remaining -= cost;
                this_frame.push(task);
            } else {
                leftover.push((task, cost));
            }
        }
        self.pending = leftover;
        this_frame
    }

    /// Number of pending tasks.
    pub fn pending_count(&self) -> usize {
        self.pending.len()
    }
}

// ---------------------------------------------------------------------------
// AsyncCompute
// ---------------------------------------------------------------------------

/// State of an async compute task.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AsyncState {
    /// Waiting to be dispatched.
    Pending,
    /// Currently running.
    Running,
    /// Completed successfully.
    Done,
    /// Failed with an error message.
    Failed(String),
}

/// A promise-like result from an async compute submission.
#[derive(Debug, Clone)]
pub struct AsyncResult {
    /// Task name.
    pub name: String,
    /// Current state.
    pub state: AsyncState,
    /// Simulated output data (bytes).
    pub output: Vec<u8>,
}

impl AsyncResult {
    /// True when the task has finished (successfully or not).
    pub fn is_complete(&self) -> bool {
        matches!(self.state, AsyncState::Done | AsyncState::Failed(_))
    }
}

/// Simulated async compute queue.
#[derive(Debug, Default)]
pub struct AsyncCompute {
    /// All submitted tasks.
    results: Vec<AsyncResult>,
}

impl AsyncCompute {
    /// Create a new async compute queue.
    pub fn new() -> Self {
        Self::default()
    }

    /// Submit a task for async execution.  Returns an index into the result list.
    pub fn submit(&mut self, task: &ComputeTask) -> usize {
        let idx = self.results.len();
        self.results.push(AsyncResult {
            name: task.name.clone(),
            state: AsyncState::Pending,
            output: vec![],
        });
        idx
    }

    /// Advance all pending tasks by one simulated tick.
    ///
    /// - Pending → Running
    /// - Running → Done (with placeholder output)
    pub fn tick(&mut self) {
        for r in &mut self.results {
            match r.state {
                AsyncState::Pending => r.state = AsyncState::Running,
                AsyncState::Running => {
                    r.state = AsyncState::Done;
                    r.output = vec![0u8; 4]; // placeholder
                }
                _ => {}
            }
        }
    }

    /// Get the result for submission index `idx`.
    pub fn poll(&self, idx: usize) -> Option<&AsyncResult> {
        self.results.get(idx)
    }

    /// Drain completed results.
    pub fn drain_completed(&mut self) -> Vec<AsyncResult> {
        let mut done = Vec::new();
        let mut remaining = Vec::new();
        for r in self.results.drain(..) {
            if r.is_complete() {
                done.push(r);
            } else {
                remaining.push(r);
            }
        }
        self.results = remaining;
        done
    }
}

// ---------------------------------------------------------------------------
// PipelineBarrier
// ---------------------------------------------------------------------------

/// Describes a memory barrier stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineStage {
    /// Beginning of the pipeline.
    Top,
    /// Vertex shading.
    Vertex,
    /// Fragment / pixel shading.
    Fragment,
    /// Compute dispatch.
    Compute,
    /// Transfer / copy operations.
    Transfer,
    /// Color attachment output.
    ColorAttachment,
    /// Shader read.
    ShaderRead,
    /// End of pipeline.
    Bottom,
}

/// A pipeline memory barrier between passes.
#[derive(Debug, Clone)]
pub struct PipelineBarrier {
    /// Stage that must complete before the barrier.
    pub src_stage: PipelineStage,
    /// Stage that must wait after the barrier.
    pub dst_stage: PipelineStage,
    /// Human-readable label.
    pub label: String,
    /// Whether this is a color attachment → shader read transition.
    pub color_to_shader_read: bool,
}

impl PipelineBarrier {
    /// Create a color attachment output → shader read barrier.
    pub fn color_attachment_to_shader_read(label: impl Into<String>) -> Self {
        Self {
            src_stage: PipelineStage::ColorAttachment,
            dst_stage: PipelineStage::ShaderRead,
            label: label.into(),
            color_to_shader_read: true,
        }
    }

    /// Create a compute → compute barrier (for storage buffer hazards).
    pub fn compute_to_compute(label: impl Into<String>) -> Self {
        Self {
            src_stage: PipelineStage::Compute,
            dst_stage: PipelineStage::Compute,
            label: label.into(),
            color_to_shader_read: false,
        }
    }

    /// True when the barrier crosses a compute → read hazard.
    pub fn is_compute_read_hazard(&self) -> bool {
        self.src_stage == PipelineStage::Compute
            && matches!(
                self.dst_stage,
                PipelineStage::ShaderRead | PipelineStage::Fragment
            )
    }
}

// ---------------------------------------------------------------------------
// GpuTimestampQuery
// ---------------------------------------------------------------------------

/// A single GPU timestamp query pair.
#[derive(Debug, Clone)]
pub struct GpuTimestampQuery {
    /// Label for this query.
    pub label: String,
    /// Simulated start time (nanoseconds).
    pub start_ns: u64,
    /// Simulated end time (nanoseconds).
    pub end_ns: u64,
    /// Whether `begin` has been called.
    active: bool,
}

impl GpuTimestampQuery {
    /// Create a new timestamp query.
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            start_ns: 0,
            end_ns: 0,
            active: false,
        }
    }

    /// Record the start timestamp.
    pub fn begin(&mut self, now_ns: u64) {
        self.start_ns = now_ns;
        self.active = true;
    }

    /// Record the end timestamp.
    pub fn end(&mut self, now_ns: u64) {
        self.end_ns = now_ns;
        self.active = false;
    }

    /// Elapsed time in microseconds.
    pub fn elapsed_us(&self) -> f64 {
        (self.end_ns.saturating_sub(self.start_ns)) as f64 / 1_000.0
    }

    /// Elapsed time in milliseconds.
    pub fn elapsed_ms(&self) -> f64 {
        self.elapsed_us() / 1_000.0
    }

    /// True when a `begin` is outstanding.
    pub fn is_active(&self) -> bool {
        self.active
    }
}

/// A pool of timestamp query pairs for profiling a frame.
#[derive(Debug, Default)]
pub struct TimestampPool {
    /// All queries.
    queries: Vec<GpuTimestampQuery>,
}

impl TimestampPool {
    /// Create an empty pool.
    pub fn new() -> Self {
        Self::default()
    }

    /// Allocate and begin a new timestamp query.  Returns its index.
    pub fn begin(&mut self, label: impl Into<String>, now_ns: u64) -> usize {
        let mut q = GpuTimestampQuery::new(label);
        q.begin(now_ns);
        let idx = self.queries.len();
        self.queries.push(q);
        idx
    }

    /// End the query at `idx`.
    pub fn end(&mut self, idx: usize, now_ns: u64) {
        if let Some(q) = self.queries.get_mut(idx) {
            q.end(now_ns);
        }
    }

    /// Get elapsed ms for query `idx`.
    pub fn elapsed_ms(&self, idx: usize) -> f64 {
        self.queries.get(idx).map(|q| q.elapsed_ms()).unwrap_or(0.0)
    }

    /// Total elapsed ms across all finished queries.
    pub fn total_ms(&self) -> f64 {
        self.queries
            .iter()
            .filter(|q| !q.is_active())
            .map(|q| q.elapsed_ms())
            .sum()
    }

    /// Reset all queries.
    pub fn reset(&mut self) {
        self.queries.clear();
    }
}

// ---------------------------------------------------------------------------
// FrameGraph
// ---------------------------------------------------------------------------

/// A transient resource in the frame graph.
#[derive(Debug, Clone)]
pub struct FrameResource {
    /// Resource name.
    pub name: String,
    /// Size in bytes.
    pub size: usize,
    /// First pass index that uses this resource.
    pub first_use: usize,
    /// Last pass index that uses this resource.
    pub last_use: usize,
    /// Allocated byte offset (set during aliasing).
    pub offset: usize,
}

/// A render pass in the frame graph.
#[derive(Debug, Clone)]
pub struct FramePass {
    /// Pass name.
    pub name: String,
    /// Resources read by this pass.
    pub reads: Vec<String>,
    /// Resources written by this pass.
    pub writes: Vec<String>,
    /// Pipeline barriers to inject before this pass.
    pub barriers: Vec<PipelineBarrier>,
}

impl FramePass {
    /// Create a new frame pass.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            reads: vec![],
            writes: vec![],
            barriers: vec![],
        }
    }

    /// Declare a resource read.
    pub fn reads(mut self, res: impl Into<String>) -> Self {
        self.reads.push(res.into());
        self
    }

    /// Declare a resource write.
    pub fn writes(mut self, res: impl Into<String>) -> Self {
        self.writes.push(res.into());
        self
    }

    /// Add a pipeline barrier.
    pub fn barrier(mut self, b: PipelineBarrier) -> Self {
        self.barriers.push(b);
        self
    }
}

/// A full-frame resource graph with transient resource aliasing.
#[derive(Debug, Default)]
pub struct FrameGraph {
    /// All passes in submission order.
    passes: Vec<FramePass>,
    /// Declared transient resources.
    resources: HashMap<String, FrameResource>,
}

impl FrameGraph {
    /// Create an empty frame graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a render pass.
    pub fn add_pass(&mut self, pass: FramePass) {
        let idx = self.passes.len();
        // Track resource lifetimes
        for res in pass.reads.iter().chain(pass.writes.iter()) {
            let e = self.resources.entry(res.clone()).or_insert(FrameResource {
                name: res.clone(),
                size: 0,
                first_use: idx,
                last_use: idx,
                offset: 0,
            });
            // Update first_use if this is earlier (handles declare_resource + add_pass order)
            if idx < e.first_use {
                e.first_use = idx;
            }
            if idx > e.last_use {
                e.last_use = idx;
            }
        }
        self.passes.push(pass);
    }

    /// Declare a transient resource with its size.
    pub fn declare_resource(&mut self, name: impl Into<String>, size: usize) {
        let name = name.into();
        let e = self.resources.entry(name.clone()).or_insert(FrameResource {
            name: name.clone(),
            size: 0,
            first_use: usize::MAX,
            last_use: 0,
            offset: 0,
        });
        e.size = size;
    }

    /// Run a simple aliasing pass: resources that do not overlap in lifetime
    /// share the same memory offset.
    pub fn alias_resources(&mut self) {
        // Greedy aliasing: sort by first_use, then assign offsets
        let names: Vec<String> = {
            let mut v: Vec<String> = self.resources.keys().cloned().collect();
            v.sort();
            v
        };

        // Track which offsets are "free" at each pass
        let mut allocations: Vec<(usize, usize, usize)> = Vec::new(); // (offset, end_pass, size)
        let pass_count = self.passes.len();

        for name in &names {
            if let Some(res) = self.resources.get_mut(name) {
                if res.first_use > pass_count {
                    continue;
                }
                // Find a free slot
                let mut found = None;
                for (off, end, sz) in &mut allocations {
                    if *end < res.first_use && *sz >= res.size {
                        found = Some(*off);
                        *end = res.last_use;
                        break;
                    }
                }
                if let Some(off) = found {
                    res.offset = off;
                } else {
                    let off: usize = allocations.iter().map(|(o, _, s)| o + s).max().unwrap_or(0);
                    res.offset = off;
                    let (last_use, size) = (res.last_use, res.size);
                    allocations.push((off, last_use, size));
                }
            }
        }
    }

    /// Compute the peak memory required (maximum end of any allocation).
    pub fn peak_memory(&self) -> usize {
        self.resources
            .values()
            .map(|r| r.offset + r.size)
            .max()
            .unwrap_or(0)
    }

    /// Number of passes.
    pub fn pass_count(&self) -> usize {
        self.passes.len()
    }

    /// Get all barriers for a given pass by index.
    pub fn barriers_for_pass(&self, idx: usize) -> &[PipelineBarrier] {
        self.passes
            .get(idx)
            .map(|p| p.barriers.as_slice())
            .unwrap_or(&[])
    }

    /// Collect all pipeline barriers across the frame in order.
    pub fn all_barriers(&self) -> Vec<&PipelineBarrier> {
        self.passes.iter().flat_map(|p| p.barriers.iter()).collect()
    }

    /// Find all resources used by pass at index `idx`.
    pub fn resources_for_pass(&self, idx: usize) -> Vec<&str> {
        if let Some(pass) = self.passes.get(idx) {
            pass.reads
                .iter()
                .chain(pass.writes.iter())
                .map(|s| s.as_str())
                .collect::<HashSet<_>>()
                .into_iter()
                .collect()
        } else {
            vec![]
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- TaskPriority tests ---

    #[test]
    fn test_priority_ordering() {
        assert!(TaskPriority::RealTime > TaskPriority::High);
        assert!(TaskPriority::High > TaskPriority::Normal);
        assert!(TaskPriority::Normal > TaskPriority::Low);
        assert!(TaskPriority::Low > TaskPriority::Background);
    }

    // --- ComputeTask tests ---

    #[test]
    fn test_compute_task_invocations_1d() {
        let t = ComputeTask::new_1d("particles", 100);
        // dispatch 100 × 1 × 1, workgroup 64 × 1 × 1 → 6400 invocations
        assert_eq!(t.total_invocations(), 6400);
    }

    #[test]
    fn test_compute_task_invocations_2d() {
        let t = ComputeTask::new_2d("shadows", 8, 8);
        // dispatch 8×8×1, workgroup 8×8×1 → 8*8*8*8 = 4096
        assert_eq!(t.total_invocations(), 4096);
    }

    #[test]
    fn test_compute_task_depends_on() {
        let t = ComputeTask::new_1d("B", 1).depends_on("A");
        assert!(t.dependencies.contains(&"A".to_owned()));
    }

    #[test]
    fn test_compute_task_priority() {
        let t = ComputeTask::new_1d("t", 1).with_priority(TaskPriority::High);
        assert_eq!(t.priority, TaskPriority::High);
    }

    // --- TaskGraph tests ---

    #[test]
    fn test_task_graph_topo_sort_simple() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("A", 1));
        g.add_task(ComputeTask::new_1d("B", 1).depends_on("A"));
        g.add_task(ComputeTask::new_1d("C", 1).depends_on("B"));
        let order = g.topological_sort().unwrap();
        let pos: HashMap<&str, usize> = order
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        assert!(pos["A"] < pos["B"]);
        assert!(pos["B"] < pos["C"]);
    }

    #[test]
    fn test_task_graph_topo_sort_diamond() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("A", 1));
        g.add_task(ComputeTask::new_1d("B", 1).depends_on("A"));
        g.add_task(ComputeTask::new_1d("C", 1).depends_on("A"));
        g.add_task(ComputeTask::new_1d("D", 1).depends_on("B").depends_on("C"));
        let order = g.topological_sort().unwrap();
        assert_eq!(order.len(), 4);
    }

    #[test]
    fn test_task_graph_cycle_detection() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("A", 1).depends_on("B"));
        g.add_task(ComputeTask::new_1d("B", 1).depends_on("A"));
        assert!(g.has_cycle());
    }

    #[test]
    fn test_task_graph_critical_path() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("A", 1).with_estimated_ms(1.0));
        g.add_task(
            ComputeTask::new_1d("B", 1)
                .depends_on("A")
                .with_estimated_ms(2.0),
        );
        g.add_task(ComputeTask::new_1d("C", 1).with_estimated_ms(10.0));
        let cp = g.critical_path();
        // C alone is the critical path (10ms vs A+B = 3ms)
        assert!(cp.contains(&"C".to_owned()));
    }

    #[test]
    fn test_task_graph_empty_topo() {
        let g = TaskGraph::new();
        let order = g.topological_sort().unwrap();
        assert!(order.is_empty());
    }

    // --- TaskScheduler tests ---

    #[test]
    fn test_scheduler_schedule() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("X", 1));
        g.add_task(ComputeTask::new_1d("Y", 1).depends_on("X"));
        let sched = TaskScheduler::new();
        let order = sched.schedule(&g).unwrap();
        assert_eq!(order.len(), 2);
    }

    #[test]
    fn test_scheduler_batch_schedule() {
        let mut g = TaskGraph::new();
        g.add_task(ComputeTask::new_1d("A", 1));
        g.add_task(ComputeTask::new_1d("B", 1));
        g.add_task(ComputeTask::new_1d("C", 1).depends_on("A").depends_on("B"));
        let sched = TaskScheduler::new();
        let batches = sched.batch_schedule(&g).unwrap();
        // A and B should be in the same batch (depth 0)
        assert!(batches[0].len() >= 2);
        // C in a later batch
        assert!(batches.len() >= 2);
    }

    // --- ResourceBarrier tests ---

    #[test]
    fn test_resource_barrier_raw() {
        let b = ResourceBarrier::raw("write_task", "read_task", "position_buffer");
        assert_eq!(b.barrier_type, BarrierType::ReadAfterWrite);
        assert_eq!(b.resource, "position_buffer");
    }

    #[test]
    fn test_resource_barrier_war() {
        let b = ResourceBarrier::war("reader", "writer", "depth");
        assert_eq!(b.barrier_type, BarrierType::WriteAfterRead);
    }

    // --- WorkloadBalancer tests ---

    #[test]
    fn test_workload_balancer_respects_budget() {
        let mut wb = WorkloadBalancer::new(10.0);
        wb.submit(ComputeTask::new_1d("A", 1).with_estimated_ms(3.0));
        wb.submit(ComputeTask::new_1d("B", 1).with_estimated_ms(4.0));
        wb.submit(ComputeTask::new_1d("C", 1).with_estimated_ms(6.0));
        let frame = wb.extract_frame_work();
        let total: f64 = frame.iter().map(|t| t.estimated_ms).sum();
        // At most budget + one overflow (for non-empty guarantee)
        assert!(total <= 10.0 + 6.0);
    }

    #[test]
    fn test_workload_balancer_priority_order() {
        let mut wb = WorkloadBalancer::new(5.0);
        wb.submit(
            ComputeTask::new_1d("low", 1)
                .with_priority(TaskPriority::Low)
                .with_estimated_ms(2.0),
        );
        wb.submit(
            ComputeTask::new_1d("rt", 1)
                .with_priority(TaskPriority::RealTime)
                .with_estimated_ms(2.0),
        );
        let frame = wb.extract_frame_work();
        // RealTime should be first
        assert_eq!(frame[0].name, "rt");
    }

    #[test]
    fn test_workload_balancer_pending_count() {
        let mut wb = WorkloadBalancer::new(1.0);
        for i in 0..5 {
            wb.submit(ComputeTask::new_1d(format!("t{i}"), 1).with_estimated_ms(1.0));
        }
        wb.extract_frame_work();
        assert!(wb.pending_count() < 5);
    }

    // --- AsyncCompute tests ---

    #[test]
    fn test_async_compute_submit_poll() {
        let mut ac = AsyncCompute::new();
        let task = ComputeTask::new_1d("sim", 64);
        let idx = ac.submit(&task);
        let r = ac.poll(idx).unwrap();
        assert_eq!(r.state, AsyncState::Pending);
    }

    #[test]
    fn test_async_compute_tick_to_done() {
        let mut ac = AsyncCompute::new();
        let task = ComputeTask::new_1d("sim", 1);
        let idx = ac.submit(&task);
        ac.tick(); // Pending → Running
        ac.tick(); // Running → Done
        assert_eq!(ac.poll(idx).unwrap().state, AsyncState::Done);
    }

    #[test]
    fn test_async_compute_drain_completed() {
        let mut ac = AsyncCompute::new();
        let t = ComputeTask::new_1d("t", 1);
        ac.submit(&t);
        ac.tick();
        ac.tick();
        let done = ac.drain_completed();
        assert_eq!(done.len(), 1);
        assert!(ac.poll(0).is_none()); // drained
    }

    // --- PipelineBarrier tests ---

    #[test]
    fn test_pipeline_barrier_color_to_shader_read() {
        let b = PipelineBarrier::color_attachment_to_shader_read("gbuffer");
        assert!(b.color_to_shader_read);
        assert_eq!(b.src_stage, PipelineStage::ColorAttachment);
        assert_eq!(b.dst_stage, PipelineStage::ShaderRead);
    }

    #[test]
    fn test_pipeline_barrier_compute_to_compute() {
        let b = PipelineBarrier::compute_to_compute("particles");
        assert_eq!(b.src_stage, PipelineStage::Compute);
        assert!(!b.is_compute_read_hazard()); // dst is also Compute
    }

    #[test]
    fn test_pipeline_barrier_compute_read_hazard() {
        let b = PipelineBarrier {
            src_stage: PipelineStage::Compute,
            dst_stage: PipelineStage::ShaderRead,
            label: "test".to_owned(),
            color_to_shader_read: false,
        };
        assert!(b.is_compute_read_hazard());
    }

    // --- GpuTimestampQuery tests ---

    #[test]
    fn test_timestamp_query_elapsed() {
        let mut q = GpuTimestampQuery::new("render");
        q.begin(1_000_000); // 1 ms in ns
        q.end(2_000_000); // 2 ms in ns
        assert!((q.elapsed_ms() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_timestamp_query_is_active() {
        let mut q = GpuTimestampQuery::new("x");
        assert!(!q.is_active());
        q.begin(0);
        assert!(q.is_active());
        q.end(100);
        assert!(!q.is_active());
    }

    #[test]
    fn test_timestamp_pool_total() {
        let mut pool = TimestampPool::new();
        let i0 = pool.begin("a", 0);
        pool.end(i0, 1_000_000);
        let i1 = pool.begin("b", 0);
        pool.end(i1, 2_000_000);
        let total = pool.total_ms();
        assert!((total - 3.0).abs() < 1e-6, "total={total}");
    }

    #[test]
    fn test_timestamp_pool_reset() {
        let mut pool = TimestampPool::new();
        pool.begin("x", 0);
        pool.reset();
        assert!((pool.total_ms()).abs() < 1e-10);
    }

    // --- FrameGraph tests ---

    #[test]
    fn test_frame_graph_add_pass() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FramePass::new("gbuffer").writes("color").writes("depth"));
        fg.add_pass(
            FramePass::new("lighting")
                .reads("color")
                .reads("depth")
                .writes("hdr"),
        );
        assert_eq!(fg.pass_count(), 2);
    }

    #[test]
    fn test_frame_graph_resource_lifetime() {
        let mut fg = FrameGraph::new();
        fg.declare_resource("color", 1024 * 1024 * 4);
        fg.add_pass(FramePass::new("p0").writes("color"));
        fg.add_pass(FramePass::new("p1").reads("color"));
        let res = &fg.resources["color"];
        assert_eq!(res.first_use, 0);
        assert_eq!(res.last_use, 1);
    }

    #[test]
    fn test_frame_graph_aliasing() {
        let mut fg = FrameGraph::new();
        fg.declare_resource("A", 1024);
        fg.declare_resource("B", 1024);
        fg.add_pass(FramePass::new("p0").writes("A"));
        fg.add_pass(FramePass::new("p1").reads("A"));
        fg.add_pass(FramePass::new("p2").writes("B"));
        fg.alias_resources();
        // B's lifetime starts after A ends, so they may share memory
        let peak = fg.peak_memory();
        assert!(peak > 0);
    }

    #[test]
    fn test_frame_graph_barriers() {
        let mut fg = FrameGraph::new();
        fg.add_pass(
            FramePass::new("render")
                .barrier(PipelineBarrier::color_attachment_to_shader_read("test")),
        );
        let barriers = fg.barriers_for_pass(0);
        assert_eq!(barriers.len(), 1);
    }

    #[test]
    fn test_frame_graph_all_barriers() {
        let mut fg = FrameGraph::new();
        fg.add_pass(FramePass::new("p0").barrier(PipelineBarrier::compute_to_compute("c0")));
        fg.add_pass(FramePass::new("p1").barrier(PipelineBarrier::compute_to_compute("c1")));
        assert_eq!(fg.all_barriers().len(), 2);
    }
}
